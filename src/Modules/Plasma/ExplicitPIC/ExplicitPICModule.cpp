#include "Modules/Plasma/IO/PicSpeciesSeries.h"
#include "Modules/Plasma/IO/PicCaseInputLoader.h"
#include "Modules/Plasma/IO/PicStepSeries.h"
#include "Modules/Plasma/IO/PicOutputManifest.h"
#include "Modules/Plasma/Diagnostics/PicStepDiagnostics.h"
#include "Modules/Plasma/Diagnostics/PicParticleTrajectories.h"
#include "Modules/Plasma/ExplicitPIC/ExplicitPICModule.h"
#include "Modules/Plasma/Physics/PicFieldState.h"
#include "Modules/Plasma/Physics/PicCurrentWeights.h"
#include "Modules/Plasma/IO/PicCheckpointStore.h"
#include "Modules/Plasma/IO/PicFieldOutput.h"
#include "Modules/Plasma/IO/PicCheckpointReader.h"
#include "Modules/Plasma/Particles/PicSpeciesState.h"
#include "Modules/Plasma/Boundary/PicBoundarySetup.h"
#include "Modules/Plasma/Boundary/PicFieldBoundaryState.h"
#include "Modules/Plasma/Boundary/PicEmbeddedGeometry.h"
#include "Modules/Plasma/Boundary/PicAmrGeometry.h"
#include "Modules/Plasma/Physics/PicDomain.h"
#include "Modules/Plasma/IO/PicCheckpointWriter.h"

#include <algorithm>
#include <cctype>
#include <array>
#include <cmath>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "Core/Mesh/FaceSegmentGraph.h"
#include "Modules/Plasma/Physics/PicMeshRefinement.h"
#include "Core/Mesh/PeriodicFaceUtils.h"
#include "Core/Mesh/QuadtreeMeshRefiner.h"
#include "Core/Field/CellField.h"
#include "Core/Field/VertField.h"
#include "Core/IO/DictionaryGeometry.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/IO/ConsoleControl.h"
#include "Core/Base/Structures.h"
#include "Modules/Plasma/Particles/ParticleType.h"
#include "Modules/Plasma/Particles/ParticleManager.h"
#include "Modules/Plasma/ExplicitPIC/Deposit/CouplingWeightsCache.h"
#include "Modules/Plasma/Physics/PicChargeDeposition.h"
#include "Modules/Plasma/ExplicitPIC/Pusher/ParticlePusherLeapfrog.h"
#include "Modules/Plasma/IO/PicPhysicsDiagnostics.h"
#include "Modules/Plasma/Diagnostics/TimeSeriesIO.h"

#include "Modules/Plasma/Physics/PhysicalConstants.h"
#include "Modules/Plasma/ExplicitPIC/Diagnostics/PicDiagnostics.h"
#include "Modules/Plasma/IO/PicProfiles.h"
#include "Solvers/Poisson/CutCellMetrics.h"
#include "Modules/Plasma/Physics/PicPoissonFieldSolve.h"
#include "Modules/Plasma/Physics/PicSecondaryEmission.h"
#include "Modules/Plasma/Physics/PicParticleSources.h"
#include "Modules/Plasma/Physics/PicSpawnManagement.h"
#include "Solvers/SolversConfig.h"
#include "Core/Base/Base.h"
#include "Core/Base/Interpolation.h"
#include "Core/Geometry/ImmersedGeometry.h"
#include "Modules/Plasma/Boundary/PicMovingGeometry.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Units/UnitsContext.h"
#include "Core/Units/NormalizationContext.h"
#include "Modules/Plasma/Physics/PicCollisions.h"
#include "Modules/Plasma/IO/PicCollisionDiagnostics.h"

using namespace std;
namespace
{



using Solvers::Poisson::CutCellMetrics;


} // namespace

ExplicitPICModule::ExplicitPICModule(const PlasmaRunInputs& inputs)
    : cfg_(inputs)
{
}

int ExplicitPICModule::withSession(PicCaseInput input, const PicSessionDriver& driver)
{
    if (!driver) throw std::invalid_argument("PIC session driver is required");
    const PicConfig& cfg = input.pic;
    const auto& dictPath = input.source;
    const auto& pmCfg = cfg.particleManagement;
    const auto& ckCfg = cfg.checkpoint;
    double dt = cfg.dt;
    int steps = cfg.steps;

    auto& dictGeometry = input.geometry;
    const auto& movingEBConfig = input.movingBoundary;
    const auto& motionSchedule = input.gridMotion;
    const bool haveGeometry = input.haveGeometry;
    const auto boundarySetup = PicBoundarySetup::prepareExplicit(cfg, dictGeometry, haveGeometry);
    const auto& domainBCMapped = boundarySetup.fields.domain;
    const auto& embeddedBCMapped = boundarySetup.fields.embedded;
    const auto& embeddedBCByIdMapped = boundarySetup.fields.embeddedById;
    const auto& lineToBoundaryId = boundarySetup.fields.lineToBoundaryId;
    const int idLeft = boundarySetup.fields.domainIds[0], idRight = boundarySetup.fields.domainIds[1];
    const int idBottom = boundarySetup.fields.domainIds[2], idTop = boundarySetup.fields.domainIds[3];
    const bool periodicX = boundarySetup.periodicX, periodicY = boundarySetup.periodicY;
    const auto& seeCfg = boundarySetup.particles.see;
    const auto& seeBoundaryIds = boundarySetup.particles.embeddedIds;

#ifdef _OPENMP
    // Allow PIC dict to set OpenMP threads (single source of truth).
    if (cfg.ompThreads > 0)
    {
        omp_set_dynamic(0);
        omp_set_num_threads(cfg.ompThreads);
    }
#endif


    using Clock = std::chrono::steady_clock;
    auto secondsBetween = [](const Clock::time_point& a, const Clock::time_point& b) -> double {
        return std::chrono::duration_cast<std::chrono::duration<double>>(b - a).count();
    };


    // O4-lite: configure console verbosity from pic.dict
    ConsoleControl::setConfig(cfg.output.console.level, cfg.output.console.every);
    const bool consoleDebug = ConsoleControl::enabled(ConsoleControl::Level::kDebug);

#ifdef _OPENMP
    {
        const int nThreads = omp_get_max_threads();
        std::cout << "[ExplicitPIC] OpenMP threads = " << nThreads << std::endl;
    }
#endif

    if (consoleDebug)
    {
        std::cout << "[ExplicitPIC] Macro-step pipeline (deposit->Poisson->E->gather->push)" << std::endl;
        std::cout << "  dt    = " << std::setprecision(16) << dt << std::endl;
        std::cout << "  steps = " << steps << std::endl;
        std::cout << "  species = " << cfg.species.size() << std::endl;
    }

    auto& dictMesh = input.mesh;
    const bool meshBlocksPresent = dictMesh.geometryRefinement.present || dictMesh.physicsRefinement.present;
    const bool meshGeoEnabled = meshBlocksPresent ? dictMesh.geometryRefinement.enabled : false;
    const bool geometryHasEmbeddedBoundaries = haveGeometry && !dictGeometry.immersedBoundaryIds().empty();
    const bool geoRefinementRequested = meshGeoEnabled &&
        (!dictMesh.geometryRefinement.boundaryRules.empty() ||
         !dictMesh.geometryRefinement.phaseRules.empty() ||
         !dictMesh.geometryRefinement.regionRules.empty());

    const auto checkpointStore = PicCheckpointStore::forRun(PicIntegrator::Explicit, input, cfg_.collisions.get());
    std::optional<PicCheckpointReader::PreparedExplicitRestart> preparedRestart;
    int outputResumeStep = 0;
    if (ckCfg.restart)
    {
        preparedRestart = PicCheckpointReader::prepareExplicitRestart(
            checkpointStore, ckCfg.restartStep, cfg_.collisions.get(), ckCfg.saveRng);
        if (!preparedRestart) return -1;
        outputResumeStep = preparedRestart->outputResumeStep();
    }

    // Output directory
    const PicOutputDirectory outputDirectory(PicOutputDirectory::Mode::Explicit, cfg.output.diagnostics.baseDir);
    const auto& outputRoot = outputDirectory.root;
    const auto& fieldsDir = outputDirectory.fields;
    const auto& seriesDir = outputDirectory.series;
    const auto& metaDir = outputDirectory.meta;
    outputDirectory.create();
    outputDirectory.prepareFields(outputResumeStep);
    const std::string seeReportPath = (metaDir / "see_report.txt").string();
    const std::string overviewPrimaryPath = (seriesDir / "overview_step.csv").string();
    const std::string overviewDetailPath = (seriesDir / "overview_detail_step.csv").string();
    const std::string speciesPrimaryPath = (seriesDir / "species_step.csv").string();
    const std::string speciesBudgetPath = (seriesDir / "species_budget_step.csv").string();
    const std::string boundaryPrimaryPath = (seriesDir / "boundary_step.csv").string();
    const std::string boundaryEnergyPath = (seriesDir / "boundary_energy_step.csv").string();
    const std::string processStepPath = (seriesDir / "process_step.csv").string();
    const std::string geometryCouplingPath = (seriesDir / "geometry_coupling_step.csv").string();
    const std::string impactReportPath = (seriesDir / "impact_step.csv").string();
    const std::string particleMgmtReportPath = (seriesDir / "particle_management_step.csv").string();
    const std::string fluxReportPath = (seriesDir / "boundary_flux_step.csv").string();
    const std::string powerReportPath = (seriesDir / "power_step.csv").string();
    const std::string rfPhaseReportPath = (seriesDir / "rf_phase_step.csv").string();
    const std::string weightStatsPath = (seriesDir / "particle_weight_stats_step.csv").string();
    const auto& diagnosticsRoot = outputDirectory.diagnostics;
    const PicOutputCadence outputCadence(PicOutputCadence::Mode::Explicit, cfg, steps);
    const int fieldsEvery = outputCadence.fieldsEvery();
    const int seriesEvery = outputCadence.seriesEvery();
    const int diagnosticsBaseEvery = outputCadence.diagnosticsEvery();
    PicParticleTrajectories trajectories;
    trajectories.configure(cfg.output.particleTrajectory, seriesEvery, seriesDir, "ExplicitPIC");
    PicSpeciesSeries speciesSeries(PicSpeciesSeries::Mode::Explicit, [&]() {
        std::vector<PicSpeciesSeries::Species> columns;
        columns.reserve(cfg.species.size());
        for (const auto& sp : cfg.species) columns.push_back({sp.name, sp.q});
        return columns;
    }());


    outputDirectory.removeObsoleteReports();

    PicOutputManifest::write(
        {cfg, outputRoot, dictPath, cfg_.units.get(),
         cfg_.normalization ? &cfg_.normalization->config() : nullptr},
        outputCadence);

    // Periodic closure diagnostics (optional)
    const std::filesystem::path periodicDiagPath = seriesDir / "periodic_step.csv";
    PicStepSeries stepSeries(PicStepSeries::Mode::Explicit);
    if (cfg.output.periodicDiagnostics.enabled)
        stepSeries.openPeriodic(periodicDiagPath, outputResumeStep);

    PicProfiles profiles(PicFieldMode::Explicit, fieldsDir, outputResumeStep, fieldsEvery,
        {cfg.output.profile1D.enabled, cfg.output.profile1D.nbinsX,
         cfg.output.profile1DAvg.enabled, cfg.output.profile1DAvg.avgWindowSteps});

    // ---------------------------------------------------------------------
    // Current-waveform (PB8.1): capture area-averaged normal field on a
    // selected domain boundary side (for displacement-current diagnostics).
    //
    // Note: we only compute <En> here; the current waveform writer and
    // J_disp/J_cond synthesis are added in later PB8.1 steps.
    // ---------------------------------------------------------------------
    const bool currentWaveEnabled = cfg.output.currentWaveform.enabled;
    const std::filesystem::path currentWavePath = seriesDir / "waveform_step.csv";
    if (currentWaveEnabled)
        stepSeries.openWaveform(currentWavePath, outputResumeStep, seriesEvery);

    // Source/sink current-balance diagnostics (optional)
    const bool currentBalanceEnabled = cfg.output.currentBalance.enabled;
    const std::filesystem::path currentBalancePath = seriesDir / "current_balance_step.csv";

    // EB near-wall RMS(|E|) diagnostics (optional)
    const bool ebBandRMSEnabled = cfg.output.ebBandRMS.enabled;
    const int ebBandLayers = std::max(0, cfg.output.ebBandRMS.bandLayers);
    const std::filesystem::path ebBandRMSPath = seriesDir / "eb_band_rms_step.csv";

    auto lowerCopy = [](std::string s) -> std::string {
        for (char& c : s)
            c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        return s;
    };


    auto boundaryIdFromDomainSide = [&](ExplicitPIC::Pusher::DomainSide side) -> int {
        switch (side)
        {
        case ExplicitPIC::Pusher::DomainSide::kLeft:
            return idLeft;
        case ExplicitPIC::Pusher::DomainSide::kRight:
            return idRight;
        case ExplicitPIC::Pusher::DomainSide::kBottom:
            return idBottom;
        case ExplicitPIC::Pusher::DomainSide::kTop:
            return idTop;
        default:
            return 0;
        }
    };



    PicElectrodeCurrent waveform(cfg.output.currentWaveform.side, domainBCMapped);
    const Boundary waveformBoundary = waveform.boundary();
    const int waveformSideIndex = waveform.sideIndex();
    if (currentWaveEnabled && waveformSideIndex < 0)
        std::cerr << "[ExplicitPIC] WARNING: CurrentWaveform enabled but cannot resolve side='"
                  << cfg.output.currentWaveform.side << "' to a domain boundary (or no Dirichlet sides for auto).\n";

    // Per-step normal field on domain boundary segments (same order as faceGraph.segments())
    std::vector<double> EnBySegmentForWaveform;











    // Collisions stats cadence (collisions.dict: stats.every). Default is 1 (every step).
    const int collisionsStatsEvery = (cfg_.collisions != nullptr) ? cfg_.collisions->stats_every : 1;
    stepSeries.openPerformance(seriesDir / "perf_step.csv", outputResumeStep);

    // (Gather+push fine-grained timing CSV removed)

    // SEE / wall-impact report.
    // Only create this file when SEE is configured to avoid redundant writes in non-SEE runs.
    if (seeCfg.enabled) stepSeries.openSecondaryEmission(seeReportPath);

    ExplicitPicCollisionReports collisionReports(metaDir, collisionsStatsEvery);

    // Domain impact diagnostic report (counts by side + unknown).
    stepSeries.openImpact(impactReportPath, outputResumeStep);

    // Particle management report (merge/split/thinning diagnostics).
    if (pmCfg.reportEnabled)
        stepSeries.openParticleManagement(particleMgmtReportPath, outputResumeStep);

    const bool writeSI = cfg_.writeSIOutputs && (cfg_.units != nullptr);

    // Boundary electrical work (Dirichlet electrodes): dE_drive = -\int_{t^n}^{t^{n+1}} \sum_{\Gamma_D} V\,I\,dt.
    // Here I is the total current flowing out of the domain through the boundary (conduction+displacement),
    // estimated from wall charge outflow (conduction) and d<En>/dt (displacement) per side.
    const bool driveWorkEnabled = writeSI && (cfg_.normalization != nullptr);
    PicBoundaryWork boundaryWork(false);

    if (currentBalanceEnabled)
        stepSeries.openCurrentBalance(currentBalancePath, outputResumeStep);

    if (ebBandRMSEnabled)
        stepSeries.openEmbeddedBand(ebBandRMSPath, outputResumeStep);

    std::unique_ptr<PicCollisions> mccState;
    const MCCCollisionOperator* mccOp = nullptr;

    PhysicsDiagnostics::PicDiagnosticsOutput physicsOutput(
        PhysicsDiagnostics::PicDiagnosticsOutput::Mode::Explicit, cfg.output.diagnostics,
        diagnosticsRoot, diagnosticsBaseEvery, outputResumeStep, cfg.species);
    const auto& consoleElectrodes = physicsOutput.consoleElectrodes();
    const auto& electrodeGroups = consoleElectrodes.groups();

    // Electrode flux/current report (per species, per side).
    speciesSeries.openFlux(fluxReportPath, outputResumeStep, writeSI);

    speciesSeries.openBoundary(boundaryPrimaryPath, outputResumeStep);

    speciesSeries.openPrimary(speciesPrimaryPath, outputResumeStep);

    stepSeries.openOverview(overviewPrimaryPath, outputResumeStep);

    speciesSeries.openBudget(speciesBudgetPath, outputResumeStep);

    speciesSeries.openBoundaryEnergy(boundaryEnergyPath, outputResumeStep);

    stepSeries.openProcess(processStepPath, outputResumeStep);

    stepSeries.openGeometry(geometryCouplingPath, outputResumeStep);

    // Power/energy report (global).
    stepSeries.openPower(powerReportPath, outputResumeStep, writeSI);

    // RF phase report (time-varying Dirichlet boundaries).
    stepSeries.openRfPhase(rfPhaseReportPath, outputResumeStep, writeSI);

    // Particle weight / per-cell occupancy stats (diagnostic).
    speciesSeries.openWeightStatistics(weightStatsPath, outputResumeStep);

    // ---------------------------------------------------------------------
    // Checkpoint / restart
    // ---------------------------------------------------------------------
    std::optional<PicCheckpointReader> restart;
    if (preparedRestart)
    {
        restart = preparedRestart->read(dt);
        if (!restart) return -1;
    }
    const bool restartReady = restart.has_value();

    // ---------------------------------------------------------------------
    // Mesh
    // ---------------------------------------------------------------------


    ExplicitPicAmrSetup amrSetup(input);
    const auto& amrCfg = amrSetup.config();
    const auto& meshRefinerDict = amrSetup.refinerMesh();

    PicDomain domain(dictMesh, meshRefinerDict, PicDomain::Mode::Explicit);
    auto& mesh = domain.mesh();
    if (!domain.prepareMesh(restartReady ? &*restart : nullptr)) return -1;
    auto& meshRefiner = domain.refiner();

    // ---------------------------------------------------------------------
    // Embedded geometry (optional)
    // ---------------------------------------------------------------------
    auto& embeddedGeometry = domain.embedded();
    PicAmrGeometry amrGeometry(PicAmrGeometry::Mode::Explicit, {});
    constexpr bool kEmbeddedAbsorbCrossing = false;
    const auto& embeddedLevelSetIndices = embeddedGeometry.levelSetIndices();
    const auto& embeddedBoundaryIds = embeddedGeometry.boundaryIds();
    const auto& cutMetrics = embeddedGeometry.metrics();
    const auto& cutAlphaByCellIndex = embeddedGeometry.alphaByCell();
    PicMovingGeometry movingGeometry(PicMovingGeometry::Mode::Explicit);
    const auto initialDomain = domain.initializeForRun(haveGeometry ? &dictGeometry : nullptr,
        geometryHasEmbeddedBoundaries, geoRefinementRequested,
        haveGeometry && (movingEBConfig.enabled || motionSchedule.enabled), lineToBoundaryId);
    if (!initialDomain) return -1;
    const bool useEmbeddedGeometry = initialDomain->embeddedEnabled;
    const bool geometryRefinementActive = initialDomain->refinementActive;
    auto* const iGeometry = domain.geometry();
    auto& faceGraph = domain.faceGraph();

    auto computeEmbeddedSurfaceAreaM2 = [&]() {
        return embeddedGeometry.surfaceAreaM2(iGeometry, dictGeometry, lineToBoundaryId,
            cfg_.units.get(), cfg_.normalization.get());
    };

    auto boundaryIdFromPosition = [&](const Vector2D& pos, int cellLabelHint) {
        return embeddedGeometry.boundaryIdAt(mesh, iGeometry, lineToBoundaryId, pos, cellLabelHint);
    };

    if (!movingGeometry.initialize(mesh, iGeometry, haveGeometry ? &dictGeometry : nullptr,
            movingEBConfig, cfg_.units.get(), cfg_.normalization.get(),
            restartReady ? std::optional<std::filesystem::path>(restart->directory()) : std::nullopt))
        return -1;

    if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
    {
        std::cout << "  mesh: verts=" << mesh.vertIndexNum()
                  << " cells=" << mesh.cellIndexNum()
                  << " leafSegments=" << faceGraph.segments().size() << std::endl;
        if (embeddedGeometry.ready())
        {
            std::ostringstream lsList;
            lsList << "{";
            bool first = true;
            for (int lsI : embeddedLevelSetIndices)
            {
                if (!first) lsList << ",";
                first = false;
                lsList << lsI;
            }
            lsList << "}";
            std::cout << "  embedded: levelSets=" << (embeddedLevelSetIndices.empty() ? std::string("{") + std::to_string(embeddedGeometry.levelSetIndex()) + "}" : lsList.str())
                      << " cutMetrics=" << (embeddedGeometry.metricsReady() ? "on" : "off")
                      << (embeddedGeometry.usesComposite() ? " (composite)" : "")
                      << std::endl;
        }
    }

    if (embeddedBCMapped.enabled && !embeddedGeometry.ready())
    {
        std::cerr << "[ExplicitPIC][EB] WARNING: BC/Embedded enabled but EmbeddedGeometry is not active; "
                  << "embedded BC will be ignored.\n";
    }

    amrSetup.activateForMesh(mesh.maxLevel(), geometryRefinementActive, std::cerr);

    if (amrSetup.enabled() && amrCfg.report.enabled)
        stepSeries.openAmrIndicators(amrCfg.report.path, outputResumeStep);

    // ---------------------------------------------------------------------
    // Fields (cell-centered storage; output interpolates to vertices)
    // ---------------------------------------------------------------------
    PicFieldState fields(mesh);
    const auto& phi = fields.potential();
    const auto& rhs = fields.poissonRhs();
    const auto& rho = fields.chargeDensity(); // Charge density; Q is deposited charge amount.
    const auto& Q = fields.chargeAmount();
    const auto& Ecell = fields.electricField();

    ExplicitPicRefinement refinement(mesh);
    ExplicitPicChargeDeposition chargeDeposition(fields);
    const auto& amrNeInstant = chargeDeposition.electronAmount();
    if (consoleDebug)
        std::cerr << "[ExplicitPIC] init: fields ready.\n";

    // ---------------------------------------------------------------------
    // Particles
    // ---------------------------------------------------------------------
    PicSpeciesState particleState(cfg.species.size(), cfg.randomSeed);
    auto& species = particleState.species();
    const auto& macroWeight = particleState.weights();

    // Deterministic RNG unless user provides a seed.
    auto& rng = particleState.rng();

    const double x0 = mesh.minCoord.x;
    const double x1 = mesh.maxCoord.x;
    const double y0 = mesh.minCoord.y;
    const double y1 = mesh.maxCoord.y;

    // ---------------------------------------------------------------------
    // Current-waveform (PB8.1): constants for electrode area and dt in SI
    // used to compute conduction current density from wall outflow.
    // ---------------------------------------------------------------------
    const auto& scalesCW = cfg_.normalization->scales();
    const double T0_s_cw = (scalesCW.T0_s > 0.0) ? scalesCW.T0_s : 1.0;
    const double L0_m_cw = (scalesCW.L0_m > 0.0) ? scalesCW.L0_m : 1.0;
    const double dt_si_forWaveform = T0_s_cw * dt;
    double electrodeArea_m2_forWaveform = 0.0;
    if (currentWaveEnabled && waveformSideIndex >= 0)
    {
        const double width_m = (cfg_.units != nullptr) ? cfg_.units->LengthCodeToSI(x1 - x0)
                                                     : (x1 - x0) * L0_m_cw;
        const double height_m = (cfg_.units != nullptr) ? cfg_.units->LengthCodeToSI(y1 - y0)
                                                      : (y1 - y0) * L0_m_cw;
        const double slab_m = (scalesCW.slabThickness_m > 0.0) ? scalesCW.slabThickness_m : 1.0;
        // 2D slab electrode area (m^2): LR -> height*slab, BT -> width*slab
        electrodeArea_m2_forWaveform = (waveformSideIndex < 2) ? (height_m * slab_m) : (width_m * slab_m);
        if (electrodeArea_m2_forWaveform <= 0.0)
        {
            std::cerr << "[ExplicitPIC] WARNING: CurrentWaveform electrode area is non-positive (area="
                      << electrodeArea_m2_forWaveform << " m^2). Check normalization slabThickness.\n";
        }
    }


    // Small margin to avoid sampling exactly on the boundary.
    const double mx = 1e-12 + 1e-6 * (x1 - x0);
    const double my = 1e-12 + 1e-6 * (y1 - y0);


    auto isInsideEmbeddedPoint = [&](const Vector2D& pos, int cellLabelHint) {
        return embeddedGeometry.containsPoint(mesh, iGeometry, pos, cellLabelHint);
    };

    std::vector<long long> stepParticleMgmtMacroLoss(species.size(), 0ll);
    std::vector<double> stepParticleMgmtPhysLoss(species.size(), 0.0);

    auto cullParticlesInsideEmbedded = [&](const char* tag,
                                           const std::vector<unsigned char>* speciesMask = nullptr,
                                           std::vector<long long>* macroLossBySpecies = nullptr,
                                           std::vector<double>* physLossBySpecies = nullptr)
    {
        particleState.cullEmbedded(domain, tag, speciesMask, macroLossBySpecies, physLossBySpecies);
    };

    const auto particleStart = particleState.initializeForRun(cfg.species, cfg_.normalization->scales(),
        domain, restartReady ? &*restart : nullptr, ckCfg.saveRng);
    if (!particleStart) return -1;
    const int startStep = *particleStart;
    const bool didRestart = restartReady;

    if (startStep >= steps)
    {
        std::cout << "[ExplicitPIC][Checkpoint] NOTE: startStep >= steps ("
                  << startStep << " >= " << steps << "), nothing to run.\n";
        return 0;
    }

    trajectories.assignMissingIds(species);
    trajectories.resetOutput();

    const PicCheckpointWriter checkpointWriter(checkpointStore, rng, ckCfg.saveRng, ckCfg.keep, "ExplicitPIC");


    // Map species name -> index for fast lookup (used by MCC spawn merging).
    particleState.rebuildIndex();
    const auto& speciesIndexByName = particleState.index();
    if (consoleDebug)
        std::cerr << "[ExplicitPIC] init: species index map ready.\n";

    PicParticleSources particleSources(PicParticleSources::Mode::Explicit,
        cfg.species, speciesIndexByName, haveGeometry ? &dictGeometry : nullptr,
        mesh, embeddedGeometry, iGeometry);

    if (consoleDebug && !particleSources.empty())
        std::cerr << "[ExplicitPIC] init: source emitters ready.\n";

    amrSetup.selectElectron(particleState, std::cerr);
    const int amrElectronIndex = amrSetup.electronIndex();
    if (!amrSetup.enabled()) stepSeries.closeAmrIndicators();
    if (consoleDebug)
        std::cerr << "[ExplicitPIC] init: AMR electron selection done.\n";

    const double amrTe_eV = amrSetup.temperatureEV(cfg.species);

    // ---------------------------------------------------------------------
    // SEE wiring: ions impacting configured boundary sides emit
    // secondary electrons with a 3D hemisphere angular distribution.
    // ---------------------------------------------------------------------
    if (consoleDebug)
        std::cerr << "[ExplicitPIC] init: SEE config begin.\n";
    const auto& scales = cfg_.normalization->scales();

    auto sideIndex = [](ExplicitPIC::Pusher::DomainSide side) -> int {
        switch (side)
        {
        case ExplicitPIC::Pusher::DomainSide::kLeft:   return 0;
        case ExplicitPIC::Pusher::DomainSide::kRight:  return 1;
        case ExplicitPIC::Pusher::DomainSide::kBottom: return 2;
        case ExplicitPIC::Pusher::DomainSide::kTop:    return 3;
        default: return -1;
        }
    };

    PicSecondaryEmission emission(seeCfg, PicSecondaryEmission::Mode::Explicit);
    emission.initialize(seeBoundaryIds, species, macroWeight, speciesIndexByName, scales.L0_m, scales.v0_mps);
    const bool seeEnabled = emission.enabled();
    const std::size_t seeEmittedIndex = emission.emittedIndex();
    if (consoleDebug)
        std::cerr << "[ExplicitPIC] init: SEE config done.\n";

    if (consoleDebug)
        std::cerr << "[ExplicitPIC] init: collisions config begin.\n";
    if (cfg_.collisions && cfg_.collisions->enabled)
    {
        std::string error;
        mccState = PicCollisions::create(
            {cfg_.gas, cfg_.crossSections, cfg_.reactions, cfg_.collisions}, cfg.species,
            cfg_.normalization ? cfg_.normalization->scales().v0_mps : 1.0, error,
            restart ? &*restart : nullptr, ckCfg.saveRng);
        if (!mccState)
        {
            std::cerr << "[ExplicitPIC] ERROR: " << error << std::endl;
            return -1;
        }
        mccOp = &mccState->operation();

        physicsOutput.prepareExplicitReactions(cfg_.reactions.get(), mccOp);

        if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
            std::cout << mccOp->SummaryString() << std::endl;

        writeExplicitPicCollisionInitialization(*mccOp, *cfg_.gas, cfg_.units.get(),
            cfg_.normalization ? cfg_.normalization->scales() : NormalizationScales{},
            mesh, cfg.species, dt, collisionReports.reportPath());
    }
    if (consoleDebug)
        std::cerr << "[ExplicitPIC] init: collisions config done.\n";

    // ---------------------------------------------------------------------
    // Coupling width hp: choose a conservative constant based on the
    // smallest leaf-cell size (paperA uses hp = local leaf-cell size).
    // ---------------------------------------------------------------------
    PicCurrentWeights::Policy coupling(mesh, cfg.deposition, PicCurrentWeights::Policy::Mode::Explicit, periodicX, periodicY);
    if (consoleDebug)
        std::cerr << "[ExplicitPIC] init: hp computed.\n";

    if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
        std::cout << "  hpConst(min leaf) = " << coupling.minimumLeafWidth() << " (local hp by host cell)\n";

    // ---------------------------------------------------------------------
    // Solvers / operators
    // ---------------------------------------------------------------------
    PicPoissonFieldSolve poisson(fields, input.poisson);
    PoissonFVSolverCellCentered::SolveReport poissonReport;
    if (consoleDebug)
        std::cerr << "[ExplicitPIC] init: poisson options ready.\n";


    particleState.configurePush(domain, useEmbeddedGeometry, periodicX, periodicY, kEmbeddedAbsorbCrossing);
    if (consoleDebug)
        std::cerr << "[ExplicitPIC] init: pusher options ready.\n";

    // ---------------------------------------------------------------------
    // Time loop skeleton
    // ---------------------------------------------------------------------
    const int consoleEvery = std::max(1, cfg.output.console.every);
    PicExplicitConsole consoleReport(cfg, startStep, steps);

    auto dumpStep = [&](int step)
    {
        PicFieldOutput::writeExplicit(step, dt,
            {mesh, faceGraph, embeddedGeometry, phi, rho, Ecell, species, macroWeight,
             cfg_.units.get(), cfg_.normalization.get(), coupling, cfg.output.fields},
            profiles, cfg_.writeSIOutputs);
    };
    if (consoleDebug)
        std::cerr << "[ExplicitPIC] init: dumpStep ready.\n";

    ExplicitPIC::Diagnostics::PicDiagnostics::Result lastDiag;
    bool lastDiagValid = false;
    // Power and series windows keep their independent output cadences.
    PicPowerWindow powerWindow;
    PicDiagnosticWindow seriesWindow(species.size());

    auto printDiagnostics = [&](int stepForPrint,
                                double totalGridCharge,
                                double totalMacroCharge,
                                const PoissonFVSolverCellCentered::SolveReport& rep,
                                const PicDomainBC& domainBC,
                                const PicEmbeddedBC& embeddedBC,
                                const PicEmbeddedBCMap& embeddedBCById)
    {
        const bool wantSeriesStep =
            outputCadence.series(stepForPrint);
        const bool wantDumpStep = (stepForPrint == 0 || outputCadence.fields(stepForPrint));
        const bool wantConsoleDebug =
            ConsoleControl::enabled(ConsoleControl::Level::kDebug) &&
            ConsoleControl::shouldPrintInfoSummary(stepForPrint);
        if (!wantSeriesStep && !wantDumpStep && !wantConsoleDebug)
            return;

        const auto sample = ExplicitPIC::Diagnostics::PicDiagnostics::sampleExplicit(
            mesh, faceGraph, domainBC, rhs, phi, Q, totalMacroCharge,
            species, macroWeight, *cfg_.normalization,
            embeddedGeometry.metricsReady() ? &cutMetrics : nullptr, embeddedBC, embeddedBCById);
        const auto& diag = sample.values;

        lastDiag = diag;
        lastDiagValid = true;

        if (wantConsoleDebug)
            PicStepSeries::writeExplicitDiagnostic(std::cout, stepForPrint,
                totalGridCharge, totalMacroCharge, sample, rep.message);

        // Write diagnostics output at the same cadence as other outputs.
        // - CSV summary
        if (outputCadence.fields(stepForPrint))
        {
            // CSV scalar diagnostics
            const double tNow = static_cast<double>(stepForPrint) * dt;
            stepSeries.writeDetail(
                overviewDetailPath, stepForPrint, tNow, diag,
                cfg_.units.get(), cfg_.writeSIOutputs, startStep);
        }
    };

    bool amrAdaptedThisStep = false;
    double amrAdapt_s_step = 0.0;
    auto updateAmrIndicators = [&](double t_si, int stepForPrint, bool adaptedThisStep)
    {
        if (!amrSetup.enabled())
            return;

        const bool useDebye = amrCfg.debye.enabled;
        const bool useEField = amrCfg.indicator.useEField;
        const bool useGeom = geometryRefinementActive;
        if (!useDebye && !useEField && !useGeom)
            return;

        if ((useDebye || useEField) && (cfg_.units == nullptr || cfg_.normalization == nullptr))
            return;
        if (useDebye && (amrElectronIndex < 0 || amrElectronIndex >= static_cast<int>(species.size())))
            return;

        if (!refinement.update(amrCfg, {amrNeInstant, Ecell,
                embeddedGeometry.metricsReady() ? &cutAlphaByCellIndex : nullptr,
                *cfg_.units, *cfg_.normalization}, amrTe_eV, stepForPrint, steps, meshRefiner)) return;
        const auto& amrStats = refinement.report();
        if (stepSeries.amrIndicatorsReady() && amrCfg.report.every > 0 &&
            ((stepForPrint % amrCfg.report.every) == 0 || stepForPrint == steps))
            stepSeries.writeAmrIndicators(stepForPrint, t_si, amrStats, refinement.pending(), adaptedThisStep);

        if (useDebye && ConsoleControl::enabled(ConsoleControl::Level::kDebug) &&
            ConsoleControl::shouldPrintPeriodic(stepForPrint, ConsoleControl::g_config.every))
            PicStepSeries::writeAmrDebyeWarning(std::cout, amrStats);
    };

    // Helper: build caches -> deposit -> rhs -> solve -> E
    auto buildSolveAtCurrentPositions = [&](PicCurrentWeights& weights,
                                           PoissonFVSolverCellCentered::SolveReport& rep,
                                           bool printDiag,
                                           int stepForPrint,
                                           double t_si,
                                           PicExplicitStepTiming* stepTiming)
    {
        PicFieldBoundaryState fieldBoundaries(domainBCMapped, embeddedBCMapped, embeddedBCByIdMapped);
        const bool needEnBySegment = (currentWaveEnabled && waveformSideIndex >= 0) || driveWorkEnabled;
        poisson.advanceExplicit(faceGraph,
            {species, macroWeight, *cfg_.normalization, coupling, embeddedGeometry, iGeometry,
             t_si, stepForPrint, amrSetup.enabled(), amrElectronIndex,
             needEnBySegment ? &EnBySegmentForWaveform : nullptr},
            weights, chargeDeposition, fieldBoundaries, rep, stepTiming);
        const auto& domainBC = fieldBoundaries.domain();
        const auto& embeddedBC = fieldBoundaries.embedded();
        const auto& embeddedBCById = fieldBoundaries.embeddedById();

        const auto tFieldSample0 = Clock::now();
        // PB8.1-4: area-weighted average normal field on selected electrode boundary.
        if (currentWaveEnabled && waveformSideIndex >= 0)
        {
            waveform.sampleNormalFieldCode(faceGraph, EnBySegmentForWaveform);
        }

        // Boundary En averages per side for electrical-work accounting.
        if (driveWorkEnabled && !EnBySegmentForWaveform.empty())
        {
            boundaryWork.sampleNormalFieldsSI(faceGraph, EnBySegmentForWaveform, *cfg_.units);
        }
        if (stepTiming)
            stepTiming->evalE_s += secondsBetween(tFieldSample0, Clock::now());

        const auto tAmrInd0 = Clock::now();
        updateAmrIndicators(t_si, stepForPrint, amrAdaptedThisStep);
        if (stepTiming)
            stepTiming->amr_indicator_s = secondsBetween(tAmrInd0, Clock::now());

        if (printDiag)
            printDiagnostics(stepForPrint, chargeDeposition.totals().grid, chargeDeposition.totals().particles, rep, domainBC, embeddedBC, embeddedBCById);
    };



    auto doAmrAdapt = [&](int stepForPrint) -> bool
    {
        if (!refinement.adaptIfNeeded(amrCfg, stepForPrint,
                {amrSetup.enabled(), useEmbeddedGeometry, geometryRefinementActive,
                 ConsoleControl::enabled(ConsoleControl::Level::kInfo),
                 domain, amrGeometry, movingGeometry, lineToBoundaryId,
                 fields, chargeDeposition, particleState, poisson, coupling}))
            return false;
        const auto& adaptation = refinement.adaptationReport();
        amrAdapt_s_step = adaptation.seconds;
        if (ConsoleControl::enabled(ConsoleControl::Level::kInfo))
            PicStepSeries::writeAmrConsole(std::cout, stepForPrint, mesh.cellIndexNum(),
                coupling.minimumLeafWidth(), adaptation);
        return true;
    };

    if (!didRestart)
    {
        // -----------------------------------------------------------------
        // Initialize leapfrog with a backward half-kick using E(x^0).
        // We first build fields at t=0, then apply: v^{-1/2} = v^{0} - kappa_s E dt/2.
        // -----------------------------------------------------------------
        PicExplicitStepTiming initFieldTiming;
        double initHalfKick_s = 0.0;
        double initDump_s = 0.0;
        const auto tInit0 = Clock::now();

        try
        {
            PicCurrentWeights initialWeights(mesh);
            const auto& cachesInit = initialWeights.caches();
            buildSolveAtCurrentPositions(initialWeights, poissonReport, /*printDiag=*/true, /*stepForPrint=*/0,
                                         /*t_si=*/0.0, &initFieldTiming);
            const auto tInitHalfKick0 = Clock::now();
            particleState.initializeBackwardHalfStep(mesh, cachesInit, Ecell, dt, iGeometry);
            initHalfKick_s = secondsBetween(tInitHalfKick0, Clock::now());
        }
        catch (const std::exception& e)
        {
            std::cerr << "[ExplicitPIC] ERROR during half-kick initialization: " << e.what() << std::endl;
            return -3;
        }

        // Dump at t=0 after fields are initialized (and velocities shifted to v^{-1/2}).
        {
            const auto tDump0 = Clock::now();
            if (outputCadence.fields(0))
                dumpStep(0);
            initDump_s = secondsBetween(tDump0, Clock::now());
        }


        // Write timing rows for the initialization (step 0).
        const double initTotal_s = secondsBetween(tInit0, Clock::now());
        if (stepSeries.performanceReady())
        {
            const std::uint64_t nP = countPicMacroParticles(species);
            stepSeries.writePerformance(0, 0.0, {initFieldTiming,
                {poissonReport.timeAssemble_s, poissonReport.timeFactor_s,
                 poissonReport.timeSolve_s, poissonReport.timeTotal_s},
                0.0, 0.0, initHalfKick_s, initDump_s, initTotal_s,
                nP, poissonReport.N, poissonReport.nnz, poissonReport.usedSuperLU});
        }
    }

    PicCollisions::Stage mccStage(mccState.get(), speciesIndexByName, macroWeight);
    const auto& mccSpawnPending = mccStage.pending();
    const auto& mccStepStats = mccStage.batchStats();
    const auto& mccCumulativeStats = mccStage.accumulatedStats();
    ParticleSpawnBuffer sourceSpawnPending;
    const PicSpawnManagement spawnManagement(PicSpawnManagement::Mode::Explicit, pmCfg, speciesIndexByName, macroWeight);

    PicCurrentWeights currentWeights(mesh);
    const auto& caches = currentWeights.caches();
    const auto tLoopStart = Clock::now();
    if (consoleDebug)
        std::cerr << "[ExplicitPIC] init: entering time loop.\n";
    const auto advanceStep = [&](const PicStepTime& stepTime) -> int
    {
        const int step = stepTime.index;
        if (consoleDebug && step == startStep)
            std::cerr << "[ExplicitPIC] loop: step begin.\n";
        const auto tStep0 = Clock::now();
        PicExplicitStepTiming fieldTiming;
        stepParticleMgmtMacroLoss.assign(species.size(), 0ll);
        stepParticleMgmtPhysLoss.assign(species.size(), 0.0);
        seriesWindow.beginFromParticles(stepTime.begin, species, macroWeight);
        // One explicit PIC macro-step:
        //   build weights -> deposit -> Poisson -> evaluate E -> matched gather -> push -> output

        amrAdapt_s_step = 0.0;
        amrAdaptedThisStep = false;
        if (amrSetup.enabled() && refinement.pending())
            amrAdaptedThisStep = doAmrAdapt(step + 1);

        try
        {
            const double t_code_step = stepTime.end;
            const double t_si_step = scales.T0_s * t_code_step;
            movingGeometry.beginStep();
            if (consoleDebug && step == startStep)
                std::cerr << "[ExplicitPIC] loop: buildSolve begin.\n";
            buildSolveAtCurrentPositions(currentWeights, poissonReport, /*printDiag=*/true, /*stepForPrint=*/(step + 1),
                                         /*t_si=*/t_si_step, &fieldTiming);
            if (consoleDebug && step == startStep)
                std::cerr << "[ExplicitPIC] loop: buildSolve done.\n";

            // PB8.2: accumulate 1D steady averages using the same particle-to-grid
            // caches as the field solve (i.e., consistent with deposition/gather).
            profiles.accumulate(step + 1,
                {mesh, cfg_.units.get(), cfg_.normalization.get(), phi, Ecell, species, macroWeight,
                 (embeddedGeometry.metricsReady() && cutAlphaByCellIndex.size() == static_cast<size_t>(mesh.cellIndexNum()))
                     ? &cutAlphaByCellIndex : nullptr}, caches);

            // PB8.1-6: driven electrode voltage for current-waveform diagnostics.
            // Re-evaluate timeDirichlet at the current SI time (t_si_step).
            if (currentWaveEnabled && waveformSideIndex >= 0)
            {
                waveform.sampleFieldSI(domainBCMapped, t_si_step, cfg_.units.get(), scales.E0_Vpm);
            }

            // Dirichlet side voltages for boundary electrical-work accounting.
            if (driveWorkEnabled)
            {
                boundaryWork.sampleVoltagesSI(domainBCMapped, *cfg_.units, t_si_step);
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "[ExplicitPIC] ERROR during field solve at step " << (step + 1)
                      << ": " << e.what() << std::endl;
            return -2;
        }

        const auto tGatherPush0 = Clock::now();

        mccStage.beginBatch();
        sourceSpawnPending.Clear();

        physicsOutput.beginStep(step + 1, cfg_.reactions.get(), mccOp);
        const bool doElectrodeDiag = physicsOutput.sampleElectrodes();
        const bool doEnergyDiag = physicsOutput.sampleEnergy();

        PicBoundaryLosses boundaryLosses(species.size(), scales.v0_mps, periodicX, periodicY, PicBoundaryLosses::Mode::Explicit);
        const auto& wallOutChargeE = boundaryLosses.totals().domainChargeE;
        const auto& wallOutChargeI = boundaryLosses.totals().domainChargeI;
        const auto& wallOutMacroBySpecies = boundaryLosses.totals().domainMacro;
        const auto& wallOutPhysBySpecies = boundaryLosses.totals().domainPhysical;
        const auto& wallAbsEnergyJBySpecies = boundaryLosses.totals().domainEnergy;
        const auto& ebOutPhysBySpecies = boundaryLosses.totals().embeddedPhysical;
        const auto& ebOutMacroBySpecies = boundaryLosses.totals().embeddedMacro;
        const auto& ebAbsEnergyJBySpecies = boundaryLosses.totals().embeddedEnergy;
        const auto& ebAbsorbStatsStep = boundaryLosses.totals().embeddedCounts;
        std::uint64_t ebResidual = 0;
        const bool ebStatsActive =
            embeddedGeometry.ready() && embeddedGeometry.absorbing() && ConsoleControl::shouldPrintInfoSummary(step + 1);
        const bool ebImpactEnergyActive = embeddedGeometry.ready() && embeddedGeometry.absorbing() && stepSeries.powerReady();
        const bool ebImpactChargeActive =
            embeddedGeometry.ready() && embeddedGeometry.absorbing() &&
            (ebStatsActive || currentBalanceEnabled || ebImpactEnergyActive || static_cast<bool>(movingGeometry.controller()));
        // (6) Matched gather + push per species + wall diagnostics/SEE
        const auto& seeStep = emission.stats();
        PicExplicitCurrentBalance currentBalance;

        std::vector<Particle> seeSpawn;
        seeSpawn.reserve(256);
        std::uniform_real_distribution<double> unif01(0.0, 1.0);

        emission.beginStep(step + 1);

        if (consoleDebug && step == startStep)
            std::cerr << "[ExplicitPIC] loop: gather/push begin.\n";

        const double v0_mps = scales.v0_mps;

        double dE_inj_J_step = 0.0;
        const auto& dE_abs_domain_J_step = boundaryLosses.totals().domainEnergyJ;
        const auto& dE_abs_embedded_J_step = boundaryLosses.totals().embeddedEnergyJ;
        const auto& qOutEmbedded_C_byBoundaryId_step = boundaryLosses.totals().embeddedChargeById;

        const double t_code = stepTime.end;
        const double t_si = scales.T0_s * t_code;
        particleState.beginPush();
        for (std::size_t si = 0; si < species.size(); ++si)
        {
            auto& sp = species[si];
            if (consoleDebug && step == startStep)
                std::cerr << "[ExplicitPIC] loop: species " << sp.name() << " begin.\n";
            const auto& cache = caches[si];
            const double w_default = (si < macroWeight.size()) ? macroWeight[si] : 1.0;
            if (sp.particle.empty())
                continue;

            const auto& pushed = particleState.pushSpecies(si, mesh, cache, Ecell, dt,
                iGeometry, ebImpactChargeActive);
            const auto& impacts = pushed.domain;
            const auto& ebImpacts = pushed.embedded;
            const auto* ebImpactsPtr = ebImpactChargeActive ? &ebImpacts : nullptr;

            if (!impacts.empty())
            {
                for (const auto& imp : impacts.impacts())
                {
                    const int idx = sideIndex(imp.side);
                    const auto impact = boundaryLosses.recordDomain(si, sp, imp.particle, w_default, idx);
                    if (!impact.accepted) continue;
                    const double w = impact.weight;

                    if (physicsOutput.includesDomainElectrodes())
                    {
                        const double qC = sp.charge() * w;
                        const double neutralW = (sp.charge() == 0.0) ? w : 0.0;
                        const int boundaryId = boundaryIdFromDomainSide(imp.side);
                        physicsOutput.recordElectrode(boundaryId, qC, neutralW);
                    }
                    if (doEnergyDiag)
                        physicsOutput.observeEnergy(sp.name(), sp.charge(), "domain",
                                                 ParticleManager::kineticEnergyEV(imp.particle, sp.mass(), v0_mps), w);

                }
            }
            emission.processExplicit({mesh, iGeometry, dictGeometry, embeddedGeometry},
                species, si, w_default, impacts, ebImpactsPtr ? &ebImpacts : nullptr,
                {{x0, y0}, {x1, y1}, {mx, my}}, lineToBoundaryId,
                seeSpawn, rng, unif01);

            // (7) MCC collisions applied after push (optional)
            // Add MCC fine-grained timing (prescale/apply/postscale/spawnproc)
            // and per-step+species counts (mcc_collisions, mcc_spawned_kept).
            if (!mccStage.applyExplicit(mesh, sp, macroWeight[si],
                    cfg_.normalization->scales().T0_s * dt, rng, step + 1)) return -4;
            if (consoleDebug && step == startStep)
                std::cerr << "[ExplicitPIC] loop: species " << sp.name() << " done.\n";

            if (ebImpactsPtr && !ebImpacts.empty())
            {
                for (const auto& imp : ebImpacts.impacts())
                {
                    if (movingGeometry.controller() && iGeometry)
                    {
                        movingGeometry.ingestImpact(mesh, *iGeometry,
                            embeddedGeometry.metricsReady() ? &cutMetrics : nullptr, lineToBoundaryId,
                            imp, static_cast<int>(si), sp.mass(), sp.charge(), w_default, v0_mps, t_si);
                    }
                    const auto impact = boundaryLosses.recordEmbedded(si, sp, imp.particle, w_default);
                    const double w = impact.weight;
                    const int bid = (driveWorkEnabled || doElectrodeDiag || ebStatsActive || !electrodeGroups.empty())
                                        ? boundaryIdFromPosition(imp.pos, imp.cellLabel)
                                        : 0;

                    boundaryLosses.recordEmbeddedBoundary(bid, impact, imp.crossing, ebStatsActive);
                    if (bid > 0)
                    {
                        const double qC = sp.charge() * w;
                        const double neutralW = (sp.charge() == 0.0) ? w : 0.0;
                        physicsOutput.recordElectrode(bid, qC, neutralW);
                    }
                    if (doEnergyDiag)
                        physicsOutput.observeEnergy(sp.name(), sp.charge(), "embedded",
                                                 ParticleManager::kineticEnergyEV(imp.particle, sp.mass(), v0_mps), w);

                }
            }
        }
        if (consoleDebug && step == startStep)
            std::cerr << "[ExplicitPIC] loop: gather/push done.\n";

        // Source injection (boundary emission).
        const auto& sc_src = cfg_.normalization->scales();
        particleSources.emit({step, sc_src.T0_s * dt, sc_src.v0_mps, sc_src.slabThickness_m},
            mesh, embeddedGeometry, iGeometry, rng, unif01, sourceSpawnPending);

        // Particle management (spawn thinning + merge/split) cadence.
        if (consoleDebug && step == startStep)
            std::cerr << "[ExplicitPIC] loop: particle mgmt begin.\n";
        const auto tPm0 = Clock::now();
        double pm_s = 0.0;
        const bool pmStep = pmCfg.enabled && (pmCfg.every > 0) && (((step + 1) % pmCfg.every) == 0);

        const auto accepted = mccStage.acceptExplicit(spawnManagement,
            {mesh, species, rng, unif01, v0_mps, {idLeft, idRight, idBottom, idTop}, pmStep},
            emission, seeSpawn, particleSources, sourceSpawnPending);
        dE_inj_J_step = accepted.injectionEnergyJ;
        const auto& spawnTouchedSpecies = accepted.touchedSpecies;
        if (consoleDebug && step == startStep)
            std::cerr << "[ExplicitPIC] loop: spawn commit done.\n";

        // Ensure newly spawned particles are not inside embedded solids.
        bool hasSpawnTouchedSpecies = false;
        for (const unsigned char touched : spawnTouchedSpecies)
        {
            if (touched != 0u)
            {
                hasSpawnTouchedSpecies = true;
                break;
            }
        }
        if (hasSpawnTouchedSpecies)
            cullParticlesInsideEmbedded("spawn", &spawnTouchedSpecies,
                                        &stepParticleMgmtMacroLoss, &stepParticleMgmtPhysLoss);
        if (consoleDebug && step == startStep)
            std::cerr << "[ExplicitPIC] loop: spawn cull done.\n";

        const double v0_mps_energyCutoff = cfg_.normalization->scales().v0_mps;
        const std::size_t energyCutoffRemovedTotal = ParticleManager::pruneByEnergyCutoff(
            species, cfg.species, macroWeight, v0_mps_energyCutoff,
            &stepParticleMgmtMacroLoss, &stepParticleMgmtPhysLoss);
        if (energyCutoffRemovedTotal > 0 && ConsoleControl::enabled(ConsoleControl::Level::kDebug))
        {
            std::cout << "[ExplicitPIC][EnergyCutoff] step=" << (step + 1)
                      << " removed=" << energyCutoffRemovedTotal << std::endl;
        }

        // Particle management (merge/split) after spawn commit.
        const bool pmReport = pmCfg.reportEnabled && stepSeries.particleManagementReady()
                              && (pmCfg.reportEvery > 0)
                              && (((step + 1) % pmCfg.reportEvery) == 0 || (step + 1) == steps);
        const bool pmConsole = pmCfg.enabled
                               && (pmCfg.reportEvery > 0)
                               && (((step + 1) % pmCfg.reportEvery) == 0 || (step + 1) == steps);
        ParticleManagementStats pmTotals;
        if (pmReport)
        {
            if (pmCfg.thinningEnabled)
            {
                stepSeries.writeSpawnManagement(step + 1, t_si, "mcc", accepted.collision, accepted.keepProbability);
                stepSeries.writeSpawnManagement(step + 1, t_si, "see", accepted.secondary, accepted.keepProbability);
                stepSeries.writeSpawnManagement(step + 1, t_si, "source", accepted.source, accepted.keepProbability);
            }
        }
        if (consoleDebug && step == startStep)
            std::cerr << "[ExplicitPIC] loop: pm report done.\n";

        const bool pmApply = pmStep && (pmCfg.mergeEnabled || pmCfg.splitEnabled);
        if (pmApply)
        {
            ParticleManager::Step management(pmCfg, mesh, faceGraph, embeddedGeometry, ParticleManager::Step::Mode::Explicit);
            ParticleManager::Step::Observer report;
            if (pmReport)
                report = [&](const ParticleType& sp, const ParticleManagementReport& values) {
                    stepSeries.writeParticleManagement(step + 1, t_si, sp.name(), values);
                };
            pmTotals = management.applyAll(species, macroWeight, scales.v0_mps, rng, pmReport, report);
        }
        if (consoleDebug && step == startStep)
            std::cerr << "[ExplicitPIC] loop: pm apply done.\n";

        // Ensure particle management does not leave particles inside solids.
        if (pmApply)
            cullParticlesInsideEmbedded("pm", nullptr,
                                        &stepParticleMgmtMacroLoss, &stepParticleMgmtPhysLoss);
        if (pmReport)
            stepSeries.flushParticleManagement();

        if (pmConsole)
        {
            std::cout << "[PM] step=" << (step + 1)
                      << " total_s=" << std::setprecision(6) << std::fixed << pm_s
                      << " keepProb=" << std::setprecision(6) << accepted.keepProbability
                      << " spawn(mcc/see in=" << accepted.collision.in << "/" << accepted.secondary.in
                      << " kept=" << accepted.collision.kept << "/" << accepted.secondary.kept << ")"
                      << " mergePairs=" << pmTotals.mergedPairs
                      << " split=" << pmTotals.splitCount
                      << " splitEmergency=" << pmTotals.splitEmergencyCount
                      << " splitCountBased=" << pmTotals.splitCountBasedCount;
            if (!pmStep)
                std::cout << " (pmStep off)";
            std::cout << std::defaultfloat << std::setprecision(6) << std::endl;
        }
        pm_s = secondsBetween(tPm0, Clock::now());
        if (consoleDebug && step == startStep)
            std::cerr << "[ExplicitPIC] loop: particle mgmt done.\n";

        // Reports and console summary (SEE/ion/neutral/collisions timing + info prints).
        if (consoleDebug && step == startStep)
            std::cerr << "[ExplicitPIC] loop: reports begin.\n";
        const auto emissionEnergy = summarizePicEmissionEnergy(
            seeStep.emittedPhys, seeStep.emittedEnergyJ, seeStep.emittedEnergyJ2);

        if (stepSeries.secondaryEmissionReady() && seeEnabled)
            stepSeries.writeSecondaryEmission(step + 1, t_si, seeStep, emissionEnergy,
                ((step + 1) % consoleEvery) == 0 || (step + 1) == steps);

        if (stepSeries.impactReady() && (((step + 1) % consoleEvery) == 0 || (step + 1) == steps))
            stepSeries.writeImpact(step + 1, t_si, boundaryLosses.totals());

        collisionReports.writeStep(step + 1, t_si, particleState, scales.v0_mps, mccOp, mccStepStats);

        // Accumulate into run total after finishing all projectile Apply() calls.
        mccStage.accumulateBatch();

        if (ebStatsActive && embeddedGeometry.ready())
        {
            std::uint64_t count = 0;
            for (const auto& sp : species)
            {
                for (const auto& p : sp.particle)
                {
                    if (isInsideEmbeddedPoint(p.position, p.inCellL))
                        ++count;
                }
            }
            ebResidual = count;
        }

        // PB8.1-5: conduction current density on selected electrode boundary from wall outflow.
        // J_cond = sum_s (q_s * N_out_phys_s) / (A_electrode * dt)
        waveform.sampleConductionSI(species, wallOutPhysBySpecies,
                electrodeArea_m2_forWaveform, dt_si_forWaveform);

        // Boundary electrical work per macro step from Dirichlet electrodes.
        // dE_drive = - sum_sides V_side * I_side * dt, where I_side = J_total * A_side.
        boundaryWork.resetStep(dt_si_forWaveform);
        if (driveWorkEnabled && dt_si_forWaveform > 0.0)
        {
            const double slab_m = (scalesCW.slabThickness_m > 0.0) ? scalesCW.slabThickness_m : 1.0;
            const double width_m = cfg_.units->LengthCodeToSI(x1 - x0);
            const double height_m = cfg_.units->LengthCodeToSI(y1 - y0);
            const double area_lr_m2 = std::max(0.0, height_m * slab_m);
            const double area_bt_m2 = std::max(0.0, width_m * slab_m);

            const PicBoundaryWork::Injection injection{accepted.sourceChargeC, seeStep.injectedCharge};
            boundaryWork.addDomain({area_lr_m2, area_lr_m2, area_bt_m2, area_bt_m2},
                species, wallOutPhysBySpecies, &injection);

            boundaryWork.addEmbedded(qOutEmbedded_C_byBoundaryId_step, embeddedBCMapped, embeddedBCByIdMapped, *cfg_.units, t_si);

            boundaryWork.finishStep();
        }

        // Accumulate interval energy terms for the next power-report row.
        // These are per-macro-step (this `step`) contributions.
        powerWindow.add({boundaryWork.current().energyJ, dE_inj_J_step,
            dE_abs_domain_J_step, dE_abs_embedded_J_step, mccStepStats.energy_loss_J_total,
            seeStep.emittedEnergyJ, seeStep.ionImpactEnergyJCathode});
        seriesWindow.addEnergy({boundaryWork.current().energyJ, dE_inj_J_step,
            dE_abs_domain_J_step, dE_abs_embedded_J_step, mccStepStats.energy_loss_J_total,
            seeStep.emittedEnergyJ, seeStep.ionImpactEnergyJCathode});
        seriesWindow.addBoundary({wallOutMacroBySpecies, wallOutPhysBySpecies, wallAbsEnergyJBySpecies,
            ebOutPhysBySpecies, ebOutMacroBySpecies, ebAbsEnergyJBySpecies,
            stepParticleMgmtMacroLoss, stepParticleMgmtPhysLoss});
        seriesWindow.addSpawns(PicDiagnosticWindow::Spawn::Source, sourceSpawnPending, speciesIndexByName, macroWeight);
        seriesWindow.addSpawns(PicDiagnosticWindow::Spawn::Collision, mccSpawnPending, speciesIndexByName, macroWeight);
        if (seeEnabled && seeEmittedIndex != static_cast<std::size_t>(-1))
            seriesWindow.addSecondary(seeEmittedIndex, seeSpawn, macroWeight);
        seriesWindow.addReactionEnergy(mccStepStats.per_channel_energy_loss_J);

        // Publish the selected electrode current every step, independently of CSV cadence.
        waveform.finishStep(dt_si_forWaveform);
        if (waveformSideIndex >= 0 && stepSeries.waveformDue(step + 1))
        {
            const double t_si_out = (cfg_.units != nullptr) ? cfg_.units->TimeCodeToSI(t_code) : t_si;
            const auto& sample = waveform.current();
            stepSeries.writeWaveform(step + 1, t_si_out,
                sample.voltageV, sample.conduction, sample.displacement);
        }

        // Periodic closure diagnostics: monitor max phi jump across periodic mates.
        if (stepSeries.periodicReady() && cfg.output.periodicDiagnostics.enabled
            && outputCadence.periodicSeries(step + 1))
        {
            const double t_si = (cfg_.units != nullptr) ? cfg_.units->TimeCodeToSI(t_code) : t_code;
            const auto jumps_code = measurePicPeriodicJumps(mesh, faceGraph, phi, periodicX, periodicY,
                embeddedGeometry.metricsReady() ? &cutMetrics : nullptr, PicPeriodicMatePolicy::AnyMate);
            const double jx_V = (cfg_.units != nullptr) ? cfg_.units->PotentialCodeToSI(jumps_code.maxJumpX)
                                                       : jumps_code.maxJumpX;
            const double jy_V = (cfg_.units != nullptr) ? cfg_.units->PotentialCodeToSI(jumps_code.maxJumpY)
                                                       : jumps_code.maxJumpY;
            stepSeries.writePeriodic(step + 1, t_si, jx_V, jy_V, periodicX, periodicY);

            if (cfg.output.periodicDiagnostics.printConsole
                && ConsoleControl::enabled(ConsoleControl::Level::kInfo)
                && ConsoleControl::shouldPrintInfoSummary(step + 1))
            {
                std::cout << "[ExplicitPIC][periodic] step=" << (step + 1)
                          << " phiJumpX=" << jx_V << " V"
                          << " phiJumpY=" << jy_V << " V" << std::endl;
            }
        }

        physicsOutput.writeReactions(step + 1, t_si, scales.T0_s * dt, mccStepStats);

        if (doEnergyDiag)
        {
            observeExplicitPicBulkEnergy(species, speciesIndexByName, macroWeight, v0_mps,
                [&](const std::string& name, double charge, double energyEV, double weight) {
                    physicsOutput.observeEnergy(name, charge, "bulk", energyEV, weight);
                });
            physicsOutput.writeEnergy(step + 1, t_si);
        }

        physicsOutput.writeElectrodes(step + 1, t_si, scales.T0_s * dt);

        // End of gather+push segment for this step.
        const auto tAfterGatherPush = Clock::now();
        const double gatherPush_s = secondsBetween(tGatherPush0, tAfterGatherPush);

        if (currentBalanceEnabled && stepSeries.currentBalanceReady() && outputCadence.periodicSeries(step + 1))
        {
            currentBalance = summarizeExplicitPicCurrentBalance(sourceSpawnPending, speciesIndexByName,
                species, macroWeight, wallOutPhysBySpecies, ebOutPhysBySpecies, scales.T0_s * dt);
            stepSeries.writeCurrentBalance(step + 1, t_si, currentBalance);
        }

        if (ebBandRMSEnabled && stepSeries.embeddedBandReady() && outputCadence.periodicSeries(step + 1))
        {
            const double e0 = (cfg_.units != nullptr) ? cfg_.units->scales().E0_Vpm : scales.E0_Vpm;
            stepSeries.writeEmbeddedBand(step + 1, t_si,
                summarizeExplicitPicEmbeddedBandField(mesh, faceGraph,
                    embeddedGeometry.metricsReady() ? &cutMetrics : nullptr, ebBandLayers, Ecell, e0));
            stepSeries.flushEmbeddedBand();
        }

        const bool seriesOutStep = outputCadence.series(step + 1);
        if (seriesOutStep)
        {
            const double t_series_si = scales.T0_s * t_code;
            const double dt_series_si = std::max(1e-300, t_series_si - scales.T0_s * seriesWindow.startTimeCode());
            const double width_si = cfg_.units->LengthCodeToSI(x1 - x0);
            const double height_si = cfg_.units->LengthCodeToSI(y1 - y0);
            const double slab_m = scales.slabThickness_m;
            const double area_lr_m2 = std::max(0.0, height_si * slab_m);
            const double area_bt_m2 = std::max(0.0, width_si * slab_m);

            double fieldEnergy = 0.0;
            if (lastDiagValid)
            {
                if (writeSI)
                    fieldEnergy = lastDiag.energyW * cfg_.units->scales().Phi0_V;
                else
                    fieldEnergy = lastDiag.energyW;
            }

            const auto kinetic = summarizeExplicitPicKineticEnergy(species, macroWeight, scales.v0_mps);
            const double kin_e = kinetic.electronJ, kin_i = kinetic.ionJ, kin_n = kinetic.neutralJ;
            const double kin_total = kinetic.total();
            const double totalEnergy_J = fieldEnergy + kin_total;
            const double jCondOut = waveform.current().ready ? waveform.current().conduction : std::numeric_limits<double>::quiet_NaN();
            const double jDispOut = waveform.current().ready ? waveform.current().displacement : std::numeric_limits<double>::quiet_NaN();
            const double jTotOut = waveform.current().ready ? (waveform.current().conduction + waveform.current().displacement)
                                                   : std::numeric_limits<double>::quiet_NaN();

            if (stepSeries.overviewReady())
            {
                stepSeries.writeOverview((step + 1), t_series_si, PicStepSeries::ExplicitOverview{
                    (lastDiagValid ? lastDiag.rG_L2 : 0.0),
                    (lastDiagValid ? lastDiag.rG_Linf : 0.0),
                    (lastDiagValid ? lastDiag.rG_L2_per_active : 0.0),
                    (lastDiagValid ? lastDiag.rG_L2_rel_rhs : 0.0),
                    (lastDiagValid ? lastDiag.rG_Linf_rel_rhs : 0.0),
                    (lastDiagValid ? lastDiag.S_eb_rel_rhs : 0.0),
                    {
                        fieldEnergy,
                        kin_e,
                        kin_i,
                        kin_n,
                        kin_total,
                        totalEnergy_J,
                        jCondOut,
                        jDispOut,
                        jTotOut,
                        (lastDiagValid ? lastDiag.macro_total : 0.0),
                        (lastDiagValid ? lastDiag.phys_total : 0.0),
                        (lastDiagValid ? lastDiag.dQrel : 0.0)
                    }
                });
            }

            if (speciesSeries.primaryReady() && lastDiagValid)
            {
                const std::size_t nSpeciesDiag = std::min(lastDiag.speciesNames.size(),
                    std::min(lastDiag.speciesMacroCount.size(),
                             std::min(lastDiag.speciesPhysCount.size(), lastDiag.speciesNmean_1pm3.size())));
                speciesSeries.writePrimary((step + 1), t_series_si, {nSpeciesDiag,
                    lastDiag.speciesMacroCount, lastDiag.speciesPhysCount, lastDiag.speciesNmean_1pm3,
                    lastDiag.electron_mean_energy_eV, lastDiag.Te_e_mean_eV});
            }

            if (speciesSeries.boundaryReady())
            {
                speciesSeries.beginBoundary((step + 1), t_series_si);
                for (std::size_t si = 0; si < species.size(); ++si)
                {
                    const auto& sp = species[si];
                    for (int s = 0; s < 4; ++s)
                    {
                        const double area = (s < 2) ? area_lr_m2 : area_bt_m2;
                        const double n_macro = static_cast<double>(seriesWindow.totals().wallOutMacroBySpecies[si][static_cast<std::size_t>(s)]);
                        const double n_phys = seriesWindow.totals().wallOutPhysBySpecies[si][static_cast<std::size_t>(s)];
                        speciesSeries.addBoundary(PicBoundaryFlux<double>::sample(
                            n_macro, n_phys, sp.charge(), area, dt_series_si));
                    }
                }
                speciesSeries.endBoundary();
            }

            if (speciesSeries.boundaryEnergyReady())
            {
                speciesSeries.beginBoundaryEnergy((step + 1), t_series_si);
                for (std::size_t si = 0; si < species.size(); ++si)
                {
                    for (int s = 0; s < 4; ++s)
                    {
                        const long long macroOut = static_cast<long long>(seriesWindow.totals().wallOutMacroBySpecies[si][static_cast<std::size_t>(s)]);
                        const double physOut = seriesWindow.totals().wallOutPhysBySpecies[si][static_cast<std::size_t>(s)];
                        const double energyJ = seriesWindow.totals().wallAbsEnergyJBySpecies[si][static_cast<std::size_t>(s)];
                        const double area_m2 = (s < 2) ? area_lr_m2 : area_bt_m2;
                        speciesSeries.addBoundaryEnergy(PicBoundaryEnergy::sample(
                            {macroOut, physOut}, energyJ, area_m2, dt_series_si));
                    }

                    const long long macroOutEmbedded = seriesWindow.totals().ebOutMacroBySpecies[si];
                    const double physOutEmbedded = seriesWindow.totals().ebOutPhysBySpecies[si];
                    const double energyEmbeddedJ = seriesWindow.totals().ebAbsEnergyJBySpecies[si];
                    const double areaEmbedded_m2 = computeEmbeddedSurfaceAreaM2();
                    speciesSeries.addBoundaryEnergy(PicBoundaryEnergy::sample(
                        {macroOutEmbedded, physOutEmbedded}, energyEmbeddedJ, areaEmbedded_m2, dt_series_si));
                }
                speciesSeries.endBoundaryEnergy();
            }

            if (speciesSeries.budgetReady())
            {
                speciesSeries.beginBudget((step + 1), t_series_si);
                for (std::size_t si = 0; si < species.size(); ++si)
                {
                    const double wDefault = (si < macroWeight.size()) ? macroWeight[si] : 1.0;
                    const auto end = PicSpeciesCount::sample(species[si], wDefault);
                    speciesSeries.addBudget(seriesWindow.budget(si, end));
                }
                speciesSeries.endBudget();
            }

            if (stepSeries.processReady())
            {
                const auto processSummary = summarizePicProcesses({
                    jCondOut,
                    jDispOut,
                    seriesWindow.totals().energy.driveJ,
                    seriesWindow.totals().energy.injectionJ,
                    seriesWindow.totals().energy.domainAbsorptionJ,
                    seriesWindow.totals().energy.embeddedAbsorptionJ,
                    seriesWindow.totals().energy.reactionLossJ,
                    dt_series_si,
                    species.size(),
                    seriesWindow.totals().wallAbsEnergyJBySpecies,
                    seriesWindow.totals().ebAbsEnergyJBySpecies,
                    seriesWindow.totals().reactionEnergyJByChannel
                });

                stepSeries.writeProcess((step + 1), t_series_si, {
                    seriesWindow.totals().energy.driveJ,
                    seriesWindow.totals().energy.injectionJ,
                    seriesWindow.totals().energy.reactionLossJ,
                    seriesWindow.totals().energy.domainAbsorptionJ,
                    seriesWindow.totals().energy.embeddedAbsorptionJ,
                    seriesWindow.totals().energy.seeEmissionJ,
                    seriesWindow.totals().energy.ionImpactJ,
                    processSummary.displacementFraction,
                    processSummary.dominantSource,
                    processSummary.dominantSourceFraction,
                    processSummary.dominantBoundary,
                    processSummary.dominantBoundaryPowerW,
                    processSummary.dominantBoundaryFraction,
                    processSummary.dominantReactionChannel,
                    processSummary.dominantReactionPowerW,
                    processSummary.dominantReactionFraction
                });
            }

            if (stepSeries.geometryReady())
            {
                const auto geometrySummary = summarizePicGeometry(embeddedGeometry.metricsReady(), cutAlphaByCellIndex,
                    mesh.cellIndexNum(), movingGeometry.controller(), step + 1);

                stepSeries.writeGeometry((step + 1), t_series_si, {
                    (movingGeometry.updated() ? 1 : 0),
                    geometrySummary.activeCutCellCount,
                    computeEmbeddedSurfaceAreaM2(),
                    geometrySummary.movingRegionCount,
                    geometrySummary.movingLevelsetCount,
                    geometrySummary.movingRegionTotalAreaM2,
                    geometrySummary.movingRegionTotalFlux,
                    geometrySummary.movingRegionMeanEnergyEV,
                    geometrySummary.movingRegionMeanCos,
                    geometrySummary.movingRegionMeanAbsSpeedMps,
                    geometrySummary.movingRegionMaxAbsSpeedMps,
                    geometrySummary.geometryGsaMm,
                    geometrySummary.geometryGadMm,
                    geometrySummary.geometryDaMm,
                    fieldEnergy,
                    seriesWindow.totals().energy.domainAbsorptionJ,
                    seriesWindow.totals().energy.embeddedAbsorptionJ,
                    seriesWindow.totals().energy.reactionLossJ
                });
            }

            seriesWindow.finishReport();
        }

        trajectories.sample(species, rng, step + 1, t_code * scales.T0_s, scales.L0_m, scales.v0_mps);

        const bool outStep = outputCadence.fields(step + 1);
        if (outStep)
        {
            const double t_out_code = t_code;
            const double t_out_si = scales.T0_s * t_code;
            const double dt_out = writeSI ? (scales.T0_s * dt) : dt;
            const double width_code = x1 - x0;
            const double height_code = y1 - y0;
            const double slab_m = scales.slabThickness_m;
            const double width_si = writeSI ? cfg_.units->LengthCodeToSI(width_code) : width_code;
            const double height_si = writeSI ? cfg_.units->LengthCodeToSI(height_code) : height_code;
            const double area_lr = (writeSI ? height_si * slab_m : height_si);
            const double area_bt = (writeSI ? width_si * slab_m : width_si);

            if (speciesSeries.fluxReady())
            {
                static const char* kSideName[4] = {"L", "R", "B", "T"};
                for (std::size_t si = 0; si < species.size(); ++si)
                {
                    const auto& sp = species[si];
                    for (int s = 0; s < 4; ++s)
                    {
                        const double area = (s < 2) ? area_lr : area_bt;
                        const double n_macro = static_cast<double>(wallOutMacroBySpecies[si][static_cast<std::size_t>(s)]);
                        const double n_phys = wallOutPhysBySpecies[si][static_cast<std::size_t>(s)];
                        const auto flux = PicBoundaryFlux<double>::sample(n_macro, n_phys, sp.charge(), area, dt_out);
                        if (writeSI)
                        {
                                speciesSeries.writeFlux((step + 1), t_out_si, kSideName[s], sp.name(),
                                    flux);
                        }
                        else
                        {
                            speciesSeries.writeFlux((step + 1), t_out_code, kSideName[s], sp.name(),
                                    flux);
                        }
                    }
                }
                if (outputCadence.explicitDetailFlush(step + 1))
                {
                    speciesSeries.flushFlux();
                }
            }

            if (stepSeries.powerReady())
            {
                double rG_L2 = 0.0;
                double rG_Linf = 0.0;
                double S_cf = 0.0;
                double S_eb = 0.0;
                double S_cf_p50 = 0.0;
                double S_cf_p90 = 0.0;
                double S_cf_p99 = 0.0;
                double S_eb_p50 = 0.0;
                double S_eb_p90 = 0.0;
                double S_eb_p99 = 0.0;
                if (lastDiagValid)
                {
                    rG_L2 = lastDiag.rG_L2;
                    rG_Linf = lastDiag.rG_Linf;
                    S_cf = lastDiag.S_cf;
                    S_eb = lastDiag.S_eb;
                    S_cf_p50 = lastDiag.S_cf_p50;
                    S_cf_p90 = lastDiag.S_cf_p90;
                    S_cf_p99 = lastDiag.S_cf_p99;
                    S_eb_p50 = lastDiag.S_eb_p50;
                    S_eb_p90 = lastDiag.S_eb_p90;
                    S_eb_p99 = lastDiag.S_eb_p99;
                }

                double fieldEnergy = 0.0;
                if (lastDiagValid)
                {
                    if (writeSI)
                        fieldEnergy = lastDiag.energyW * cfg_.units->scales().Phi0_V;
                    else
                        fieldEnergy = lastDiag.energyW;
                }

                const double t_ref = writeSI ? t_out_si : t_out_code;
                const auto kinetic = summarizeExplicitPicKineticEnergy(species, macroWeight, scales.v0_mps);
                const double kin_e = kinetic.electronJ, kin_i = kinetic.ionJ, kin_n = kinetic.neutralJ;
                const double kin_total = kinetic.total();

                const auto power = powerWindow.sample(t_ref, fieldEnergy, kin_total);
                const double jCondOut = waveform.current().ready ? waveform.current().conduction : std::numeric_limits<double>::quiet_NaN();
                const double jDispOut = waveform.current().ready ? waveform.current().displacement : std::numeric_limits<double>::quiet_NaN();
                const double jTotOut = waveform.current().ready ? (waveform.current().conduction + waveform.current().displacement)
                                                       : std::numeric_limits<double>::quiet_NaN();

                stepSeries.writePower(step + 1, t_ref, {
                    rG_L2, rG_Linf, lastDiagValid ? lastDiag.rG_L2_per_active : 0.0,
                    lastDiagValid ? lastDiag.rG_L2_rel_rhs : 0.0, lastDiagValid ? lastDiag.rG_Linf_rel_rhs : 0.0,
                    S_cf, S_eb, lastDiagValid ? lastDiag.S_eb_rel_rhs : 0.0,
                    {S_cf_p50, S_cf_p90, S_cf_p99}, {S_eb_p50, S_eb_p90, S_eb_p99},
                    fieldEnergy, kin_e, kin_i, kin_n, kin_total, power.balance,
                    powerWindow.transfers(), power.potentialCurrentW, power.potentialEnergyJ});
                powerWindow.finishReport();
                if (outputCadence.explicitDetailFlush(step + 1))
                    stepSeries.flushPower();
            }

            if (speciesSeries.weightStatisticsReady())
            {
                for (std::size_t si = 0; si < species.size(); ++si)
                    speciesSeries.writeWeightStatistics(step + 1, t_out_si, species[si].name(),
                        summarizePicParticleWeights(mesh, species[si], si < macroWeight.size() ? macroWeight[si] : 1.0));
                if (outputCadence.explicitDetailFlush(step + 1))
                    speciesSeries.flushWeightStatistics();
            }

            if (stepSeries.rfPhaseReady())
            {
                stepSeries.writeRfPhase(step + 1, writeSI ? t_out_si : t_out_code,
                    summarizePicRfPhases(domainBCMapped, embeddedBCMapped, t_out_si,
                        writeSI ? cfg_.units.get() : nullptr));
                if (outputCadence.explicitDetailFlush(step + 1))
                    stepSeries.flushRfPhase();
            }
        }
        if (consoleDebug && step == startStep)
            std::cerr << "[ExplicitPIC] loop: reports done.\n";
        if (movingGeometry.advance(step + 1, t_si, scales.T0_s * dt, mesh, iGeometry))
        {
            movingGeometry.refreshAfterUpdate(step + 1,
                {mesh, *iGeometry, embeddedGeometry, faceGraph, lineToBoundaryId, useEmbeddedGeometry,
                 nullptr, poisson, nullptr});
            refinement.requestAdaptation();
        }
        const auto tOut0 = Clock::now();
        if (outputCadence.fields(step + 1))
        {
            if (consoleDebug && step == startStep)
                std::cerr << "[ExplicitPIC] loop: output begin.\n";
            if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
                std::cout << "[ExplicitPIC] step " << (step + 1) << "/" << steps
                          << " (t=" << (stepTime.end) << ")" << std::endl;
            dumpStep(step + 1);
            if (mccOp) collisionReports.append(step + 1, cfg_.units->TimeCodeToSI(stepTime.end),
                particleState, scales.v0_mps, mccStepStats, mccCumulativeStats);
            if (consoleDebug && step == startStep)
                std::cerr << "[ExplicitPIC] loop: output done.\n";
        }
        const double output_s = secondsBetween(tOut0, Clock::now());
        const double total_s = secondsBetween(tStep0, Clock::now());
        const std::uint64_t nP = countPicMacroParticles(species);
        const int stepsCompleted = std::max(1, step + 1 - startStep);
        const int stepsRemain = std::max(0, steps - (step + 1));
        const double elapsedLoop_s = secondsBetween(tLoopStart, Clock::now());
        const double eta_s = (stepsCompleted > 0)
                                 ? (elapsedLoop_s / static_cast<double>(stepsCompleted)) * static_cast<double>(stepsRemain)
                                 : std::numeric_limits<double>::quiet_NaN();
        const PicStepSeries::ExplicitPerformance performance{fieldTiming,
            {poissonReport.timeAssemble_s, poissonReport.timeFactor_s,
             poissonReport.timeSolve_s, poissonReport.timeTotal_s},
            amrAdapt_s_step, pm_s, gatherPush_s, output_s, total_s,
            nP, poissonReport.N, poissonReport.nnz, poissonReport.usedSuperLU};
        consoleReport.writeStep(step + 1, t_code, dt, eta_s, performance,
            {{mesh, phi, cfg_.units.get(), particleState, mccStepStats, boundaryLosses.totals(),
              chargeDeposition.totals().grid, chargeDeposition.totals().particles,
              x1 - x0, y1 - y0, scales.slabThickness_m, ebStatsActive},
             consoleElectrodes, lastDiagValid ? &lastDiag : nullptr,
             {seeStep.ionImpactMacroCathode, seeStep.emittedMacro, seeStep.ionImpactPhysCathode,
              seeStep.emittedPhys, emissionEnergy.meanEV, emissionEnergy.stdEV},
             {static_cast<long long>(sourceSpawnPending.size()), currentBalance.sourceC, currentBalance.domainC, currentBalance.embeddedC,
              seeEnabled, useEmbeddedGeometry && embeddedGeometry.ready(), mccOp != nullptr,
              amrAdaptedThisStep, outputCadence.fields(step + 1), ebResidual}}, std::cout);
        if (stepSeries.performanceReady() && (((step + 1) % consoleEvery) == 0 || (step + 1) == steps))
        {
            stepSeries.writePerformance(step + 1, t_code, performance);
        }
        if (consoleDebug && step == startStep)
            std::cerr << "[ExplicitPIC] loop: timing done.\n";



        return 0;
    };
    const auto completeStep = [&](const PicStepTime& stepTime, bool saveCheckpoint) -> int
    {
        const int step = stepTime.index;
        if (consoleDebug && step == startStep)
            std::cerr << "[ExplicitPIC] loop: checkpoint begin.\n";
        if (saveCheckpoint) checkpointWriter.write(step + 1, stepTime.end, stepTime.dt,
            {domain, movingGeometry, species, mccState.get()});
        if (consoleDebug && step == startStep)
            std::cerr << "[ExplicitPIC] loop: checkpoint done.\n";
        return 0;
    };
    const int sessionResult = drivePicStepSession(startStep, steps, dt, advanceStep, driver, completeStep);
    if (sessionResult != 0) return sessionResult;

    // Ensure final dump even if steps is not a multiple of output.every.
    if (outputCadence.explicitFinalFields())
    {
        dumpStep(steps);
        if (mccOp) collisionReports.append(steps, cfg_.units->TimeCodeToSI(steps * dt),
            particleState, scales.v0_mps, mccStepStats, mccCumulativeStats);
    }

    trajectories.finish();

    if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
        std::cout << "[ExplicitPIC] Macro-step pipeline finished." << std::endl;
    return 0;
}
