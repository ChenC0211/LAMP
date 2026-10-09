#include "Modules/Plasma/IO/PicSpeciesSeries.h"
#include "Modules/Plasma/IO/PicCaseInputLoader.h"
#include "Modules/Plasma/IO/PicStepSeries.h"
#include "Modules/Plasma/IO/PicOutputManifest.h"
#include "Modules/Plasma/Diagnostics/PicStepDiagnostics.h"
#include "Modules/Plasma/Diagnostics/PicParticleTrajectories.h"
#include "Modules/Plasma/ImplicitPIC/ImplicitPICModule.h"
#include "Modules/Plasma/IO/PicCheckpointStore.h"
#include "Modules/Plasma/IO/PicFieldOutputQueue.h"
#include "Modules/Plasma/IO/PicFieldOutput.h"
#include "Modules/Plasma/IO/PicCheckpointReader.h"
#include "Modules/Plasma/Particles/PicSpeciesState.h"
#include "Modules/Plasma/Boundary/PicBoundarySetup.h"
#include "Modules/Plasma/Boundary/PicFieldBoundaryState.h"
#include "Modules/Plasma/Boundary/PicEmbeddedGeometry.h"
#include "Modules/Plasma/Boundary/PicAmrGeometry.h"
#include "Modules/Plasma/Physics/PicDomain.h"
#include "Modules/Plasma/Boundary/PicScheduledGeometry.h"
#include "Modules/Plasma/IO/PicCheckpointWriter.h"
#include "Modules/Plasma/Physics/PicFieldState.h"
#include "Modules/Plasma/Physics/PicCurrentWeights.h"
#include "Modules/Plasma/ImplicitPIC/ImplicitPicPolarization.h"

#include <array>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <cstdint>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif


#include "Core/Field/CellField.h"
#include "Core/Field/HEdgeField.h"
#include "Core/Field/VEdgeField.h"
#include "Core/Field/VertField.h"
#include "Core/Geometry/ImmersedGeometry.h"
#include "Modules/Plasma/Boundary/PicMovingGeometry.h"
#include "Modules/Plasma/Boundary/LegacyGridMotionEvolution.h"
#include "Core/IO/ConsoleControl.h"
#include "Core/IO/DictionaryGeometry.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Modules/Plasma/Physics/PicMeshRefinement.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/QuadtreeMeshRefiner.h"
#include "Modules/Plasma/Particles/ParticleManager.h"
#include "Core/Base/Interpolation.h"
#include "Core/Units/UnitsContext.h"
#include "Core/Units/NormalizationContext.h"

#include "Modules/Plasma/ExplicitPIC/Deposit/CouplingWeightsCache.h"

#include "Modules/Plasma/Diagnostics/PicParticleFields.h"
#include "Modules/Plasma/IO/PicProfiles.h"
#include "Modules/Plasma/ExplicitPIC/Geometry/LevelSetQuery.h"
#include "Modules/Plasma/ExplicitPIC/Diagnostics/PicDiagnostics.h"
#include "Modules/Plasma/ExplicitPIC/Pusher/ParticlePusherLeapfrog.h"
#include "Modules/Plasma/IO/PicPhysicsDiagnostics.h"
#include "Modules/Plasma/Diagnostics/TimeSeriesIO.h"
#include "Modules/Plasma/Physics/PicCollisions.h"
#include "Modules/Plasma/Physics/PhysicalConstants.h"

#include "Modules/Plasma/Physics/PicPoissonFieldSolve.h"
#include "Modules/Plasma/Physics/PicSecondaryEmission.h"
#include "Modules/Plasma/Physics/PicParticleSources.h"
#include "Modules/Plasma/Physics/PicSpawnManagement.h"
#include "Modules/Plasma/ImplicitPIC/ImplicitPicFieldSolve.h"
#include "Solvers/Poisson/CutCellMetrics.h"
#include "Solvers/SolversConfig.h"

#include "Modules/Plasma/ImplicitPIC/Orbit/OrbitAveragedWeightsBuilder.h"
#include "Modules/Plasma/ImplicitPIC/Orbit/ImplicitPicOrbit.h"

using Solvers::Poisson::CutCellMetrics;

static std::string toLowerCopy(std::string s)
{
    for (auto& c : s)
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}



ImplicitPICModule::ImplicitPICModule(const PlasmaRunInputs& inputs)
    : cfg_(inputs)
{
}

int ImplicitPICModule::withSession(PicCaseInput input, const PicSessionDriver& driver)
{
    if (!driver) throw std::invalid_argument("PIC session driver is required");
    const PicConfig& cfg = input.pic;
    const auto& overrides = input.implicitOverrides;
    const auto& dictPath = input.source;
    const auto& pmCfg = cfg.particleManagement;
    const auto& ckCfg = cfg.checkpoint;
    double dt = cfg.dt;
    int steps = cfg.steps;
    ConsoleControl::setConfig(cfg.output.console.level, cfg.output.console.every);
    const bool consoleInfo = ConsoleControl::enabled(ConsoleControl::Level::kInfo);
    const bool consoleDebug = ConsoleControl::enabled(ConsoleControl::Level::kDebug);

#ifdef _OPENMP
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

    if (!(dt > 0.0) || steps <= 0)
    {
        std::cerr << "[ImplicitPIC] ERROR: invalid dt/steps in " << dictPath << " (dt=" << dt
                  << ", steps=" << steps << ")\n";
        return -2;
    }

    if (!cfg_.units || !cfg_.normalization)
    {
        std::cerr << "[ImplicitPIC] ERROR: Units/Normalization context is missing.\n";
        return -2;
    }

    auto& dictGeometry = input.geometry;
    const auto& movingEBConfig = input.movingBoundary;
    const auto& motionSchedule = input.gridMotion;
    const bool haveGeometry = input.haveGeometry;
    auto& dictMesh = input.mesh;
    const bool meshBlocksPresent = dictMesh.geometryRefinement.present || dictMesh.physicsRefinement.present;
    const bool meshGeoEnabled = meshBlocksPresent ? dictMesh.geometryRefinement.enabled : false;
    const bool geometryHasEmbeddedBoundaries = haveGeometry && !dictGeometry.immersedBoundaryIds().empty();
    const bool geoRefinementRequested = meshGeoEnabled &&
        (!dictMesh.geometryRefinement.boundaryRules.empty() ||
         !dictMesh.geometryRefinement.phaseRules.empty() ||
         !dictMesh.geometryRefinement.regionRules.empty());

    const auto checkpointStore = PicCheckpointStore::forRun(PicIntegrator::Implicit, input, cfg_.collisions.get());

    std::optional<PicCheckpointReader> restart;
    if (ckCfg.restart)
    {
        restart = PicCheckpointReader::select(checkpointStore, ckCfg.restartStep, dt, "ImplicitPIC");
        if (!restart) return -1;
    }
    const bool restartReady = restart.has_value();

    // ---------------------------------------------------------------------
    // Geometry.dict (optional) + mesh.dict
    // ---------------------------------------------------------------------



    PicDomain domain(dictMesh, dictMesh, PicDomain::Mode::Implicit);
    auto& mesh = domain.mesh();
    if (!domain.prepareMesh(restartReady ? &*restart : nullptr,
        overrides.initGeomSnapshot.value_or(false)
            ? std::optional<std::filesystem::path>(overrides.initGeomSnapshotPath.value_or("./cache/init_geom_snapshot"))
            : std::nullopt)) return -1;

    const auto fieldBoundarySetup = PicBoundarySetup::mapFields(cfg, dictGeometry, haveGeometry,
        PicBoundarySetup::Mode::Implicit);
    const auto& domainBCMapped = fieldBoundarySetup.domain;
    const auto& embeddedBCMapped = fieldBoundarySetup.embedded;
    const auto& embeddedBCByIdMapped = fieldBoundarySetup.embeddedById;
    const auto& lineToBoundaryId = fieldBoundarySetup.lineToBoundaryId;

    // ---------------------------------------------------------------------
    // Embedded geometry + cut metrics + AMR refiner wiring
    // ---------------------------------------------------------------------
    auto& embeddedGeometry = domain.embedded();
    constexpr bool kEmbeddedAbsorbCrossing = false;
    const auto& embeddedLevelSetIndices = embeddedGeometry.levelSetIndices();
    const auto& embeddedBoundaryIds = embeddedGeometry.boundaryIds();
    const auto& cutMetrics = embeddedGeometry.metrics();
    const auto& cutAlphaByCellIndex = embeddedGeometry.alphaByCell();
    PicMovingGeometry movingGeometry(PicMovingGeometry::Mode::Implicit);
    auto& meshRefiner = domain.refiner();
    const auto initialDomain = domain.initializeForRun(haveGeometry ? &dictGeometry : nullptr,
        geometryHasEmbeddedBoundaries, geoRefinementRequested,
        haveGeometry && (movingEBConfig.enabled || motionSchedule.enabled), lineToBoundaryId,
        movingEBConfig.enabled ? nullptr : &motionSchedule);
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

    // Runtime BC state (time-varying Dirichlet values are updated in-place each step).
    PicFieldBoundaryState fieldBoundaries(domainBCMapped, embeddedBCMapped, embeddedBCByIdMapped);
    const auto& domainBCUsed = fieldBoundaries.domain();
    const auto& embeddedBCUsed = fieldBoundaries.embedded();
    const auto& embeddedBCByIdUsed = fieldBoundaries.embeddedById();

    // Periodicity switches for the PIC loop.
    const bool periodicX = (domainBCUsed.left.type == PicBCType::Periodic) &&
                           (domainBCUsed.right.type == PicBCType::Periodic);
    const bool periodicY = (domainBCUsed.bottom.type == PicBCType::Periodic) &&
                           (domainBCUsed.top.type == PicBCType::Periodic);

    if (consoleInfo)
    {
        std::cout << "[ImplicitPIC] Loaded config from " << dictPath << "\n"
                  << "[ImplicitPIC] dt=" << std::setprecision(16) << dt
                  << " steps=" << steps
                  << " periodicX=" << (periodicX ? "true" : "false")
                  << " periodicY=" << (periodicY ? "true" : "false")
                  << " EB=" << ((useEmbeddedGeometry && embeddedGeometry.ready()) ? "on" : "off")
                  << std::defaultfloat << std::setprecision(6)
                  << std::endl;
    }

#ifdef _OPENMP
    if (consoleInfo)
        std::cout << "[ImplicitPIC] OpenMP threads = " << omp_get_max_threads() << std::endl;
#endif

    // ---------------------------------------------------------------------
    // Fields
    // ---------------------------------------------------------------------
    PicFieldState fields(mesh);
    const auto& phi = fields.potential();
    const auto& rhs = fields.poissonRhs();
    const auto& rho = fields.chargeDensity(); // Charge density; Q is deposited charge amount.
    const auto& Q = fields.chargeAmount();
    const auto& Ecell = fields.electricField();

    // ---------------------------------------------------------------------
    // Particles
    // ---------------------------------------------------------------------
    PicSpeciesState particleState(cfg.species.size(), cfg.randomSeed);
    auto& species = particleState.species();
    const auto& macroWeight = particleState.weights();

    auto& rng = particleState.rng();

    const double x0 = mesh.minCoord.x;
    const double x1 = mesh.maxCoord.x;
    const double y0 = mesh.minCoord.y;
    const double y1 = mesh.maxCoord.y;
    const double mx = 1e-12 + 1e-6 * (x1 - x0);
    const double my = 1e-12 + 1e-6 * (y1 - y0);
    std::uniform_real_distribution<double> unif01(0.0, 1.0);

    const NormalizationContext& norm = *cfg_.normalization;
    const auto& sc = norm.scales();
    const double v0_mps = std::max(1e-300, sc.v0_mps);
    const double L0_m = std::max(1e-300, sc.L0_m);
    const double T0_s = std::max(1e-300, sc.T0_s);
    const double Phi0_V = sc.Phi0_V;
    const double dt_si = dt * T0_s;

    std::vector<long long> stepParticleMgmtMacroLoss;
    std::vector<double> stepParticleMgmtPhysLoss;

    auto cullParticlesInsideEmbedded = [&](const char* tag,
                                          std::vector<long long>* macroLossBySpecies = nullptr,
                                          std::vector<double>* physLossBySpecies = nullptr) {
        return particleState.cullEmbedded(domain, tag, nullptr, macroLossBySpecies, physLossBySpecies);
    };

    const auto particleStart = particleState.initializeForRun(cfg.species, sc,
        domain, restartReady ? &*restart : nullptr, ckCfg.saveRng);
    if (!particleStart) return -1;
    const int startStep = *particleStart;
    const bool didRestart = restartReady;

    // Map species name -> index for fast lookup (used by MCC spawn merging).
    particleState.rebuildIndex();
    const auto& speciesIndexByName = particleState.index();

    // ---------------------------------------------------------------------
    // Minimal SEE wiring for embedded boundaryIds.
    // ---------------------------------------------------------------------
    const auto particleBoundarySetup = PicBoundarySetup::prepareImplicitParticles(cfg, dictGeometry);
    const auto& seeCfg = particleBoundarySetup.see;
    const auto& seeBoundaryIds = particleBoundarySetup.embeddedIds;

    const bool seeDebug = overrides.seeDebug.value_or(false);
    PicSecondaryEmission emission(seeCfg, PicSecondaryEmission::Mode::Implicit);
    emission.initialize(seeBoundaryIds, species, macroWeight, speciesIndexByName, L0_m, v0_mps, seeDebug);
    const bool seeEnabled = emission.enabled();

    // ---------------------------------------------------------------------
    // Optional boundary source injection (geometry-driven).
    // ---------------------------------------------------------------------
    PicParticleSources particleSources(PicParticleSources::Mode::Implicit,
        cfg.species, speciesIndexByName, haveGeometry ? &dictGeometry : nullptr,
        mesh, embeddedGeometry, iGeometry, useEmbeddedGeometry);

    ParticleSpawnBuffer sourceSpawnPending;
    ParticleSpawnBuffer seeSpawnPending;

    using SpawnPmKindStats = PicSpawnManagement::Stats;
    const PicSpawnManagement spawnManagement(PicSpawnManagement::Mode::Implicit, pmCfg, speciesIndexByName, macroWeight);

    using Spawn = PicStepParticleBalance::Spawn;

    const double v0_mps_energyCutoff = cfg_.normalization->scales().v0_mps;


    // ---------------------------------------------------------------------
    // Optional MCC wiring (Strang splitting: half-step before + after push).
    // ---------------------------------------------------------------------
    std::unique_ptr<PicCollisions> mccState;
    const MCCCollisionOperator* mccOp = nullptr;

    if (cfg_.collisions && cfg_.collisions->enabled)
    {
        std::string error;
        mccState = PicCollisions::create(
            {cfg_.gas, cfg_.crossSections, cfg_.reactions, cfg_.collisions}, cfg.species,
            cfg_.normalization ? cfg_.normalization->scales().v0_mps : 1.0, error,
            restart ? &*restart : nullptr, ckCfg.saveRng);
        if (!mccState)
        {
            std::cerr << "[ImplicitPIC] ERROR: " << error << std::endl;
            return -1;
        }
        mccOp = &mccState->operation();

        if (cfg_.collisions->enabled)
            std::cout << "[ImplicitPIC][MCC] enabled" << std::endl;
    }


    PicCollisions::Stage mccStage(mccState.get(), speciesIndexByName, macroWeight);
    const auto& mccSpawnPending = mccStage.pending();
    const auto& mccStepStats = mccStage.accumulatedStats();


    // ---------------------------------------------------------------------
    // Coupling width hp
    // ---------------------------------------------------------------------
    PicCurrentWeights::Policy coupling(mesh, cfg.deposition, PicCurrentWeights::Policy::Mode::Implicit, periodicX, periodicY);
    const auto& hp = coupling.halfWidths();

    // ---------------------------------------------------------------------
    // Solvers / operators
    // ---------------------------------------------------------------------
    auto poissonSetup = input.poisson;
    // Debug hardening: the implicit macro-step currently calls Poisson twice per step
    // with potentially changing coefficients later. Disable caching for now to
    // avoid any cache-signature edge cases while bringing the pipeline up.
    poissonSetup.cacheEnabled = false;
    // Small cut-cell treatment: cells with alpha < threshold get a Laplacian
    // averaging stencil instead of the degenerate flux equation.  This avoids
    // phi≈0 artifacts at EB aperture corners where beta_f → 0.
    poissonSetup.alphaMinThreshold = 0.05;
    PicPoissonFieldSolve poisson(fields, poissonSetup);
    const auto& poissonOpt = poisson.options();
    ImplicitPicFieldSolve polPoisson(fields, poissonOpt, cfg.polarizationPoisson,
        domain.hasScheduledMotion(), consoleInfo, overrides.polarization);

    const bool useCutMetrics = useEmbeddedGeometry && embeddedGeometry.ready() && embeddedGeometry.metricsReady();
    const CutCellMetrics* cutMetricsPtr = useCutMetrics ? &cutMetrics : nullptr;
    const PicEmbeddedBC* embeddedBCPtr = (useCutMetrics && embeddedBCUsed.enabled) ? &embeddedBCUsed : nullptr;
    const PicEmbeddedBCMap* embeddedBCByIdPtr = (useCutMetrics && !embeddedBCByIdUsed.empty()) ? &embeddedBCByIdUsed : nullptr;
    ImmersedGeometry<QuadtreeMesh>* iGeometryPtr = (useEmbeddedGeometry && embeddedGeometry.ready() && iGeometry) ? iGeometry : nullptr;
    if (!movingGeometry.initialize(mesh, iGeometry, haveGeometry ? &dictGeometry : nullptr,
            movingEBConfig, cfg_.units.get(), cfg_.normalization.get(),
            restartReady ? std::optional<std::filesystem::path>(restart->directory()) : std::nullopt))
        return -1;

    particleState.configurePush(domain, useEmbeddedGeometry, periodicX, periodicY, kEmbeddedAbsorbCrossing);

    const PicPoissonFieldSolve::ImplicitRebuildInput fieldRebuild{
        species, macroWeight, norm, coupling, embeddedGeometry, useCutMetrics};

    // ---------------------------------------------------------------------
    // Init: build fields at t=0 then backward half-kick -> v^{-1/2}
    // ---------------------------------------------------------------------
    if (!didRestart)
    {
        PicCurrentWeights rebuiltWeights(mesh);
        poisson.rebuildImplicit({faceGraph, domainBCUsed, cutMetricsPtr, embeddedBCPtr, embeddedBCByIdPtr},
            fieldRebuild, rebuiltWeights);
        const auto& caches0 = rebuiltWeights.caches();

        particleState.initializeBackwardHalfStep(mesh, caches0, Ecell, dt, iGeometryPtr);

        if (consoleInfo)
            std::cout << "[ImplicitPIC] init: Poisson + backward half-kick done" << std::endl;
    }
    else if (consoleInfo)
    {
        std::cout << "[ImplicitPIC][Checkpoint] restart state loaded; skip backward half-kick init" << std::endl;
    }

    PicParticleTrajectories trajectories;
    trajectories.assignMissingIds(species);

    // ---------------------------------------------------------------------
    // Diagnostics outputs
    // ---------------------------------------------------------------------
    const PicOutputCadence outputCadence(PicOutputCadence::Mode::Implicit, cfg, steps);
    const int fieldsEvery = outputCadence.fieldsEvery();
    const int seriesEvery = outputCadence.seriesEvery();
    const int amountEvery = outputCadence.amountEvery();
    PicSpeciesSeries speciesSeries(PicSpeciesSeries::Mode::Implicit, [&]() {
        std::vector<PicSpeciesSeries::Species> columns;
        columns.reserve(species.size());
        for (const auto& sp : species) columns.push_back({sp.name(), sp.charge()});
        return columns;
    }());


    const PicOutputDirectory outputDirectory(PicOutputDirectory::Mode::Implicit, cfg.output.diagnostics.baseDir);
    const auto& outputRoot = outputDirectory.root;
    const auto& fieldsDir = outputDirectory.fields;
    const auto& seriesDir = outputDirectory.series;
    const auto& metaDir = outputDirectory.meta;
    const auto& diagnosticsRoot = outputDirectory.diagnostics;
    const int diagnosticsBaseEvery = outputCadence.diagnosticsEvery();
    trajectories.configure(cfg.output.particleTrajectory, seriesEvery, seriesDir, "ImplicitPIC");
    // restart->metadata().nextStep is the next loop index to execute, while step-tagged outputs
    // are written with the completed physical step number (step + 1). Preserve the
    // checkpoint-aligned outputs and only prune data strictly after the restored state.
    const int outputResumeStep = restartReady ? std::max(0, restart->metadata().nextStep + 1) : 0;
    outputDirectory.create();
    trajectories.resetOutput();
    if (restartReady)
    {
        std::cout << "[ImplicitPIC][Output] Restart resume checkpoint step=" << restart->metadata().nextStep
                  << "; preserve outputs <= " << restart->metadata().nextStep
                  << ", prune outputs >= " << outputResumeStep << std::endl;
    }
    outputDirectory.prepareFields(outputResumeStep);
    PicStepSeries stepSeries(PicStepSeries::Mode::Implicit);
    stepSeries.openPerformance(seriesDir / "perf_step.csv", outputResumeStep);

    stepSeries.openAmrAdaptation(seriesDir / "amr_adapt_step.csv", outputResumeStep);
    stepSeries.openAmrLevelSets(seriesDir / "amr_levelset_diag.csv", outputResumeStep);

    stepSeries.openPhysics(seriesDir / "physics_step.csv", outputResumeStep);

    const std::string overviewPrimaryPath = (seriesDir / "overview_step.csv").string();
    stepSeries.openOverview(overviewPrimaryPath, outputResumeStep);

    const std::string processStepPath = (seriesDir / "process_step.csv").string();
    stepSeries.openProcess(processStepPath, outputResumeStep);

    const std::string geometryCouplingPath = (seriesDir / "geometry_coupling_step.csv").string();
    stepSeries.openGeometry(geometryCouplingPath, outputResumeStep);

    const std::string speciesDetailPath = (seriesDir / "species_detail_step.csv").string();
    speciesSeries.openDetail(speciesDetailPath, outputResumeStep);

    const std::string speciesPrimaryPath = (seriesDir / "species_step.csv").string();
    speciesSeries.openPrimary(speciesPrimaryPath, outputResumeStep);

    const std::string speciesBudgetPath = (seriesDir / "species_budget_step.csv").string();
    speciesSeries.openBudget(speciesBudgetPath, outputResumeStep);

    if (pmCfg.reportEnabled)
        stepSeries.openParticleManagement(seriesDir / "particle_management_step.csv", outputResumeStep);

    if (amountEvery > 0)
    {
        stepSeries.openAmount(seriesDir / "amount_step.csv", outputResumeStep);
        stepSeries.openAmountSpecies(seriesDir / "amount_species_step.csv", outputResumeStep);
    }

    stepSeries.openCouplingFallback(seriesDir / "coupling_fallback_step.csv", outputResumeStep);

    const std::string boundaryFluxPath = (seriesDir / "boundary_flux_step.csv").string();
    speciesSeries.openFlux(boundaryFluxPath, outputResumeStep, true);

    const std::string boundaryPrimaryPath = (seriesDir / "boundary_step.csv").string();
    speciesSeries.openBoundary(boundaryPrimaryPath, outputResumeStep);

    const std::string boundaryEnergyPath = (seriesDir / "boundary_energy_step.csv").string();
    speciesSeries.openBoundaryEnergy(boundaryEnergyPath, outputResumeStep);

    const bool currentBalanceEnabled = cfg.output.currentBalance.enabled;
    std::filesystem::path currentBalancePath = seriesDir / "current_balance_step.csv";

    if (currentBalanceEnabled)
        stepSeries.openCurrentBalance(currentBalancePath, outputResumeStep);

    PhysicsDiagnostics::PicDiagnosticsOutput physicsOutput(
        PhysicsDiagnostics::PicDiagnosticsOutput::Mode::Implicit, cfg.output.diagnostics,
        diagnosticsRoot, diagnosticsBaseEvery, outputResumeStep, cfg.species, cfg_.reactions.get(), mccOp);
    const auto& consoleElectrodes = physicsOutput.consoleElectrodes();
    const auto& electrodeGroups = consoleElectrodes.groups();

    outputDirectory.removeObsoleteReports();

    const bool geometryMotionEnabled =
        haveGeometry &&
        motionSchedule.enabled &&
        !movingEBConfig.enabled &&
        iGeometry &&
        !iGeometry->levelSets.empty();
    const int geometryMetricsEvery =
        std::max(1, motionSchedule.outputEvery > 0 ? motionSchedule.outputEvery : seriesEvery);
    if (geometryMotionEnabled)
        stepSeries.openGeometryMetrics(seriesDir, motionSchedule.metricsPath, outputResumeStep, geometryMetricsEvery);

    if (cfg.output.ebBandRMS.enabled)
        stepSeries.openEmbeddedBand(seriesDir / "eb_band_rms_step.csv", outputResumeStep);

    // Periodic closure diagnostics (optional).
    const std::filesystem::path periodicDiagPath = seriesDir / "periodic_step.csv";

    if (cfg.output.periodicDiagnostics.enabled)
        stepSeries.openPeriodic(periodicDiagPath, outputResumeStep);

    PicOutputManifest::write(
        {cfg, outputRoot, dictPath, cfg_.units.get(),
         cfg_.normalization ? &cfg_.normalization->config() : nullptr},
        outputCadence);

    const bool periodicDiagDebug = overrides.periodicDiagDebug.value_or(false);

    PicProfiles profiles(PicFieldMode::Implicit, fieldsDir, outputResumeStep, fieldsEvery,
        {cfg.output.profile1D.enabled, cfg.output.profile1D.nbinsX,
         cfg.output.profile1DAvg.enabled, cfg.output.profile1DAvg.avgWindowSteps});

    // Current waveform diagnostics (optional).
    const bool currentWaveEnabled = cfg.output.currentWaveform.enabled;
    const std::filesystem::path currentWavePath = seriesDir / "waveform_step.csv";
    if (currentWaveEnabled)
        stepSeries.openWaveform(currentWavePath, outputResumeStep, seriesEvery);

    PicElectrodeCurrent waveform(cfg.output.currentWaveform.side, domainBCUsed);
    const Boundary waveformBoundary = waveform.boundary();
    const int waveformSideIndex = waveform.sideIndex();
    if (currentWaveEnabled && waveformSideIndex < 0)
        std::cerr << "[ImplicitPIC] WARNING: CurrentWaveform enabled but cannot resolve side='"
                  << cfg.output.currentWaveform.side << "' to a domain boundary (or no Dirichlet sides for auto).\n";

    // Per-step normal field on domain boundary segments (same order as faceGraph.segments()).
    std::vector<double> EnBySegmentForWaveform;








    // Boundary electrical work (Dirichlet electrodes): dE_drive = -\int \sum_{\Gamma_D} V I dt.
    const bool driveWorkEnabled = (cfg_.units != nullptr) && (cfg_.normalization != nullptr);
    PicBoundaryWork boundaryWork(true);

    double electrodeArea_m2_forWaveform = 0.0;
    const double dt_si_forWaveform = (cfg_.units != nullptr) ? cfg_.units->TimeCodeToSI(dt) : dt_si;
    if (currentWaveEnabled && (cfg_.units != nullptr) && (cfg_.normalization != nullptr) && waveformSideIndex >= 0)
    {
        const double slabCfg_m = cfg_.normalization->config().slabThickness_m;
        const double slab_m = (slabCfg_m > 0.0) ? slabCfg_m : 1.0;
        const double Lx_m = cfg_.units->LengthCodeToSI(mesh.maxCoord.x - mesh.minCoord.x);
        const double Ly_m = cfg_.units->LengthCodeToSI(mesh.maxCoord.y - mesh.minCoord.y);
        if (waveformBoundary == Boundary::west || waveformBoundary == Boundary::east)
            electrodeArea_m2_forWaveform = std::max(0.0, Ly_m * slab_m);
        else if (waveformBoundary == Boundary::south || waveformBoundary == Boundary::north)
            electrodeArea_m2_forWaveform = std::max(0.0, Lx_m * slab_m);
    }

    const std::string implicitPowerReportPath = (seriesDir / "power_step.csv").string();
    stepSeries.openPower(implicitPowerReportPath, outputResumeStep);
    PicStepEnergyHistory energyHistory;

    PicImplicitConsole consoleReport(cfg, overrides, startStep, steps, fieldsEvery, consoleInfo);
    const int depositFillOmpMinParticles =
        std::max(1, overrides.depositFillOmpMinParticles.value_or(4096));

    // ---------------------------------------------------------------------
    // Main macro-step loop
    // ---------------------------------------------------------------------
    PicCurrentWeights currentWeights(mesh);
    const auto& cachesN = currentWeights.caches();
    ImplicitPicOrbit orbit(mesh, species.size(), cfg, overrides, cfg_.units.get(), periodicX, periodicY);
    const auto& cachesBar = orbit.averagedCaches();
    ImplicitPicPolarization polarization(mesh, species.size(), faceGraph.segments().size());
    const auto& Qstar = polarization.chargeAmount();
    const auto& rhoStar = polarization.chargeDensity();
    const auto& Nstar = polarization.numberAmounts();
    const auto& muCell_Fpm = polarization.susceptibilitySI();
    const auto& kappaByCell = polarization.cellCoefficients();
    const auto& kappaBySeg = polarization.faceCoefficients();

    ImplicitPicFieldHistory fieldHistory(fields);
    const auto& EcellNew = fieldHistory.nextElectricField();
    const auto& Ehalf = fieldHistory.halfElectricField();

    if (restartReady)
    {
        PicCurrentWeights restartWeights(mesh);
        if (poisson.restoreImplicit(restart->directory(), {faceGraph, domainBCUsed, cutMetricsPtr, embeddedBCPtr, embeddedBCByIdPtr},
                fieldRebuild, restartWeights, fieldHistory))
        {
            if (consoleInfo)
                std::cout << "[ImplicitPIC][Checkpoint] field state restored from "
                          << PicCheckpointStore::fieldName << std::endl;
        }
        else
        {
            std::cerr << "[ImplicitPIC][Checkpoint] WARNING: field state missing; rebuilt phi/E from restart particles."
                      << std::endl;
        }
    }

    const int cutcellStateDumpStep = overrides.cutcellStateDumpStep.value_or(-1);
    const std::string cutcellStateDumpPathRaw =
        overrides.cutcellStateDumpPath.value_or(std::string());
    const bool cutcellStateDumpOnlyActive =
        overrides.cutcellStateDumpOnlyActive.value_or(true);
    PicScheduledGeometry scheduledGeometry(domain.takeScheduledMotion(),
        {mesh, iGeometry, embeddedGeometry, faceGraph, lineToBoundaryId, movingGeometry.controller(),
         meshRefiner, poisson, polPoisson},
        {geometryMotionEnabled, useEmbeddedGeometry, geometryRefinementActive,
         [&](int step, const std::string& message) { stepSeries.writeGeometryTrace(step, message); }});
    stepSeries.openGeometryTrace(overrides.movingTracePath.value_or(std::string()),
        overrides.movingTraceStep.value_or(-1), overrides.movingTraceAll.value_or(false));

    PicFieldOutputQueue packedFieldOutput(mesh, cfg_.units.get(), cfg.output.fields,
        PicFieldOutputQueue::Options::fromImplicit(overrides));

    // ---------------------------------------------------------------------
    // AMR (mesh.dict PhysicsRefinement/GeometryRefinement)
    // ---------------------------------------------------------------------
    const auto amrSetup = resolveImplicitPicAmr(input, particleState.electronIndex(), geometryRefinementActive);
    const auto& refinementOptions = amrSetup.indicators;
    const bool amrEnabled = refinementOptions.enabled;
    const int amrAdaptEvery = amrSetup.adaptEvery;
    const int amrIndicatorEvery = refinementOptions.indicatorEvery;
    packedFieldOutput.configureBatching(amrEnabled);

    ImplicitPicRefinement refinement(mesh);



    const PicCheckpointWriter checkpointWriter(checkpointStore, rng, ckCfg.saveRng, ckCfg.keep, "ImplicitPIC");

    amrSetup.writeReport(std::cout, std::cerr, consoleInfo);

    auto updateAmrTargetsIfNeeded = [&](int stepForPrint) {
        refinement.update(refinementOptions, {Nstar, Ecell, useCutMetrics ? &cutAlphaByCellIndex : nullptr,
                          cfg_.units.get(), norm}, stepForPrint, meshRefiner);
    };


    PicAmrGeometry amrGeometry(PicAmrGeometry::Mode::Implicit, amrSetup.geometry);

    auto doAmrAdaptIfNeeded = [&](int stepForPrint, PicImplicitStepTiming* stepTiming = nullptr) -> bool {
        if (!refinement.adaptIfNeeded(stepForPrint,
                {amrEnabled, amrSetup.fastRemap, useEmbeddedGeometry, geometryRefinementActive,
                 amrSetup.resetTransients, amrAdaptEvery, domain, amrGeometry, movingGeometry,
                 lineToBoundaryId, fields, fieldHistory, polarization, particleState,
                 poisson, coupling, currentWeights, orbit,
                 {faceGraph, domainBCUsed, cutMetricsPtr, embeddedBCPtr, embeddedBCByIdPtr}, fieldRebuild,
                 [&](const QuadtreeMesh& adaptedMesh) {
                     packedFieldOutput.refreshCoordinates(adaptedMesh);
                 }})) return false;
        const auto& adaptation = refinement.adaptationReport();
        const auto& geometryRefresh = amrGeometry.report();
        stepSeries.writeAmrAdaptation(stepForPrint, adaptation, geometryRefresh);
        stepSeries.writeAmrLevelSets(stepForPrint, iGeometry, geometryRefresh);
        if (ConsoleControl::enabled(ConsoleControl::Level::kInfo))
            PicStepSeries::writeAmrConsole(std::cout, stepForPrint, adaptation, geometryRefresh);
        
        // Record subphase timings if stepTiming provided
        if (stepTiming)
        {
            stepTiming->amr_faceGraph_s = adaptation.faceGraphS;
            stepTiming->amr_cutMetrics_s = adaptation.embeddedMetricsS;
            stepTiming->amr_reinitLevelset_s = geometryRefresh.geometryRefreshS;
            stepTiming->amr_geom_initializeLevelSets_s = geometryRefresh.geomInitializeLevelSetsS;
            stepTiming->amr_geom_pdeReinitLevelset_s = geometryRefresh.geomPdeReinitS;
            stepTiming->amr_geom_refreshDerived_s = geometryRefresh.geomRefreshDerivedS;
            stepTiming->amr_geom_checkCutInfo_s = geometryRefresh.geomCheckCutInfoS;
            stepTiming->amr_geom_rebuildComposite_s = geometryRefresh.geomRebuildCompositeS;
            stepTiming->amr_geom_markPhases_s = geometryRefresh.geomMarkPhasesS;
            stepTiming->amr_fieldRebuild_s = adaptation.fieldRebuildS;
        }
        return true;
    };

    auto applyScheduledGeometry = [&](int step, double time_s, bool force) {
        if (scheduledGeometry.advance(step, time_s, force) && amrEnabled)
        {
            refinement.requestAdaptation();
        }
    };

    if (startStep >= steps)
    {
        std::cout << "[ImplicitPIC][Checkpoint] NOTE: startStep >= steps ("
                  << startStep << " >= " << steps << "), nothing to run.\n";
        return 0;
    }

    if (!(restartReady && scheduledGeometry.hasSchedule()))
    {
        applyScheduledGeometry(
            startStep,
            (cfg_.units != nullptr) ? cfg_.units->TimeCodeToSI(static_cast<double>(startStep) * dt)
                                    : (static_cast<double>(startStep) * dt_si),
            /*force=*/true);
    }
    stepSeries.writeGeometryMetrics(
        startStep,
        (cfg_.units != nullptr) ? cfg_.units->TimeCodeToSI(static_cast<double>(startStep) * dt)
                                : (static_cast<double>(startStep) * dt_si),
        scheduledGeometry.current(),
        /*force=*/true);

    const auto advanceStep = [&](const PicStepTime& stepTime) -> int
    {
        const int step = stepTime.index;
        try
        {
        stepParticleMgmtMacroLoss.assign(species.size(), 0ll);
        stepParticleMgmtPhysLoss.assign(species.size(), 0.0);
        const auto tStep0 = Clock::now();
        PicImplicitStepTiming stepTiming;
        PicStepParticleBalance particleBalance(species, macroWeight);
        const auto& sourceSpawns = particleBalance.spawns(Spawn::Source);
        const auto& mccSpawns = particleBalance.spawns(Spawn::Collision);
        const auto& secondarySpawns = particleBalance.spawns(Spawn::Secondary);
        const auto& stepCharge = particleBalance.current();
        PicBoundaryLosses boundaryLosses(species.size(), v0_mps, periodicX, periodicY, PicBoundaryLosses::Mode::Implicit);
        const auto& dE_abs_domain_J = boundaryLosses.totals().domainEnergyJ;
        const auto& dE_abs_embedded_J = boundaryLosses.totals().embeddedEnergyJ;
        emission.beginStep(step + 1);
        const auto& seeStep = emission.stats();
        PicBoundaryLosses::ImplicitWall wall;
        const auto& wallOutMacroBySpeciesStep = boundaryLosses.totals().domainMacro;
        const auto& wallOutPhysBySpeciesStep = boundaryLosses.totals().domainPhysical;
        const auto& wallAbsEnergyJBySpeciesStep = boundaryLosses.totals().domainEnergy;
        const auto& ebOutMacroBySpeciesStep = boundaryLosses.totals().embeddedMacro;
        const auto& ebOutPhysBySpeciesStep = boundaryLosses.totals().embeddedPhysical;
        const auto& ebAbsEnergyJBySpeciesStep = boundaryLosses.totals().embeddedEnergy;
        const auto& embeddedChargeOutByBoundaryIdStep_C = boundaryLosses.totals().embeddedChargeById;
        SpawnPmKindStats pmMccSpawnStats;
        SpawnPmKindStats pmSeeSpawnStats;
        SpawnPmKindStats pmSourceSpawnStats;
        seeSpawnPending.Clear();
        if (consoleDebug && ((step % std::max(1, cfg.output.console.every)) == 0))
            std::cout << "[ImplicitPIC] stepBegin=" << (step + 1) << "/" << steps << std::endl;

        const double t_code_step = stepTime.end;
        const double t_si_step = (cfg_.units != nullptr) ? cfg_.units->TimeCodeToSI(t_code_step)
                                                         : (step + 1) * dt_si;
        movingGeometry.beginStep();
        const bool reportStep = outputCadence.series(step + 1);
        const bool pmStep = pmCfg.enabled && (pmCfg.every > 0) && (((step + 1) % pmCfg.every) == 0);
        const bool pmReport = pmCfg.reportEnabled && stepSeries.particleManagementReady()
                              && (pmCfg.reportEvery > 0)
                              && ((((step + 1) % pmCfg.reportEvery) == 0) || ((step + 1) == steps));
        const bool pmConsole = pmCfg.enabled
                               && (pmCfg.reportEvery > 0)
                               && ((((step + 1) % pmCfg.reportEvery) == 0) || ((step + 1) == steps));
        const bool summaryStepDefault = ConsoleControl::shouldPrintInfoSummary(step + 1);

        applyScheduledGeometry(step, t_si_step, /*force=*/false);

        // Apply time-varying Dirichlet BC values before the implicit field solve.
        fieldBoundaries.update(t_si_step);

        if (mccOp)
            mccStage.clearAccumulated();

        physicsOutput.beginStep(step + 1);
        const bool doElectrodeDiag = physicsOutput.sampleElectrodes();
        const bool doEnergyDiag = physicsOutput.sampleEnergy();

        const PicSpawnManagement::ImplicitStage spawnStage{
            mesh, species, rng, unif01, particleBalance, sc.v0_mps, pmStep};

        // (0) MCC half-step before field solve/predictor (Strang splitting).
        const auto tMccPre0 = Clock::now();
        mccStage.advanceImplicitHalf(spawnManagement, spawnStage, 0.5 * dt_si, &pmMccSpawnStats);
        if (mccOp)
            mccStage.accumulateBatch();
        stepTiming.mcc_pre_s += secondsBetween(tMccPre0, Clock::now());

        // (1) Build weights at x^n.
        if (consoleDebug && step == 0)
            std::cout << "[ImplicitPIC] stage: buildCachesAtCurrentPositions" << std::endl;
        const auto tBuildWeightsN0 = Clock::now();
        stepTiming.buildWeightsN_cwFallbackCount = 0u;
        stepTiming.buildWeightsN_particleCount = 0u;
        const auto currentWeightCounts = currentWeights.rebuild(species, coupling, embeddedGeometry, useCutMetrics);
        stepTiming.buildWeightsN_cwFallbackCount = currentWeightCounts.fallback;
        stepTiming.buildWeightsN_particleCount = currentWeightCounts.particles;
        stepTiming.buildWeightsN_s += secondsBetween(tBuildWeightsN0, Clock::now());
        if (consoleDebug && step == 0)
            std::cout << "[ImplicitPIC] stage: gatherEn" << std::endl;
        if (step == 0)
            ImplicitPicOrbit::validateCurrentCaches(species, cachesN);

        // (2) Gather E^n.
        const auto tGatherEn0 = Clock::now();
        orbit.gather(species, cachesN, Ecell, step);
        stepTiming.gatherEn_s += secondsBetween(tGatherEn0, Clock::now());

        // (3) Predictor segments for each species.
        if (consoleDebug && step == 0)
            std::cout << "[ImplicitPIC] stage: predictorSegments" << std::endl;

        const auto orbitSpeciesPerf = orbit.build(species, hp, dt, embeddedGeometry, useCutMetrics,
            particleState.pushOptions().absorbAtEmbedded || particleState.pushOptions().absorbEmbeddedCrossing, {step, consoleInfo, consoleDebug});

        for (const auto& perf : orbitSpeciesPerf)
        {
            stepTiming.predictorPrep_s += perf.predictorPrep_s;
            stepTiming.predictorSegment_s += perf.predictorSegment_s;
            stepTiming.predictor_segments_total += perf.predictor_segments_total;
            stepTiming.predictor_fast_exit_hits += perf.predictor_fast_exit_hits;
            stepTiming.predictor_bisection_fallback_hits += perf.predictor_bisection_fallback_hits;
            stepTiming.predictor_coarsen_hits += perf.predictor_coarsen_hits;
            stepTiming.predictor_coarsen_parent_climb += perf.predictor_coarsen_parent_climb;

            stepTiming.orbitAvgWeights_s += perf.orbitAvgWeights_s;
            stepTiming.orbitAvg_fused_mode = std::max(stepTiming.orbitAvg_fused_mode, perf.orbitAvg_fused_mode);
            stepTiming.orbitAvg_fused_eval_s += perf.orbitAvg_fused_eval_s;
            stepTiming.orbitAvg_fused_accum_s += perf.orbitAvg_fused_accum_s;
            stepTiming.orbitAvg_sampleFill_s += perf.orbitAvg_sampleFill_s;
            stepTiming.orbitAvg_sampleCacheBuild_s += perf.orbitAvg_sampleCacheBuild_s;
            stepTiming.orbitAvg_localReduce_s += perf.orbitAvg_localReduce_s;
            stepTiming.orbitAvg_concat_s += perf.orbitAvg_concat_s;
            stepTiming.orbitAvg_sampleCwFallbackCount += perf.orbitAvg_sampleCwFallbackCount;
            stepTiming.orbitAvg_sampleCount += perf.orbitAvg_sampleCount;
            stepTiming.orbitAvg_sampleCache_nnz += perf.orbitAvg_sampleCache_nnz;
            stepTiming.orbitAvg_out_nnz += perf.orbitAvg_out_nnz;
            stepTiming.orbitAvg_adaptAcceptCount += perf.orbitAvg_adaptAcceptCount;
            stepTiming.orbitAvg_adaptRefineCount += perf.orbitAvg_adaptRefineCount;
            stepTiming.orbitAvg_adaptHitMaxDepthCount += perf.orbitAvg_adaptHitMaxDepthCount;
            stepTiming.orbitAvg_adaptAcceptedDepthSum += perf.orbitAvg_adaptAcceptedDepthSum;
        }
        const bool needEnBySegment = (currentWaveEnabled && waveformSideIndex >= 0) || driveWorkEnabled;
        const auto polResult = polPoisson.advance(
            {faceGraph, domainBCUsed, cutMetricsPtr, embeddedBCPtr, embeddedBCByIdPtr},
            {species, macroWeight, cachesBar, norm, dt_si, step, depositFillOmpMinParticles,
             periodicX, periodicY, consoleDebug, needEnBySegment ? &EnBySegmentForWaveform : nullptr},
            polarization, poisson, fieldHistory, particleBalance, stepTiming);
        const auto& polRep = polResult.report;
        if (cutcellStateDumpStep > 0 && step + 1 == cutcellStateDumpStep)
            PicFieldOutput::writeCutcellState(outputRoot, cutcellStateDumpPathRaw, step + 1, t_si_step,
                cutcellStateDumpOnlyActive,
                {mesh, rhoStar, kappaByCell, phi, EcellNew, cutMetricsPtr, embeddedBCPtr, embeddedBCByIdPtr,
                 iGeometryPtr ? &iGeometryPtr->levelSets : nullptr, cfg_.units.get(), cfg_.normalization.get(),
                 poissonOpt.embeddedBCMissingPolicy == PoissonFVSolverCellCentered::Options::EmbeddedBCMissingPolicy::Skip});

        // Current-waveform helper values from the field solve: <En> and driven voltage.
        const auto tWavePeriodicSample0 = Clock::now();
        if (currentWaveEnabled && waveformSideIndex >= 0)
        {
            // Area-weighted average normal field on selected electrode boundary.
            waveform.sampleNormalFieldCode(faceGraph, EnBySegmentForWaveform);

            // Driven electrode voltage.
            waveform.sampleFieldSI(domainBCUsed, t_si_step, cfg_.units.get());
        }

        // Boundary En averages per side for electrical-work accounting.
        if (driveWorkEnabled && !EnBySegmentForWaveform.empty())
        {
            boundaryWork.sampleNormalFieldsSI(faceGraph, EnBySegmentForWaveform, *cfg_.units);

            boundaryWork.sampleVoltagesSI(domainBCUsed, *cfg_.units);
        }

        // Cache periodic phi-jump diagnostics from the same field-solve snapshot
        // used for waveform <En> so we can write them later in the step in the
        // same order as ExplicitPIC (after current-waveform bookkeeping).
        bool periodicDiagPendingWrite = false;
        double periodicDiagJx_V = 0.0;
        double periodicDiagJy_V = 0.0;
        if (stepSeries.periodicReady() && cfg.output.periodicDiagnostics.enabled
            && outputCadence.periodicSeries(step + 1))
        {
            const auto jumps_code = measurePicPeriodicJumps(mesh, faceGraph, phi, periodicX, periodicY,
                cutMetricsPtr, PicPeriodicMatePolicy::CleanMate, periodicDiagDebug);
            periodicDiagJx_V = (cfg_.units != nullptr) ? cfg_.units->PotentialCodeToSI(jumps_code.maxJumpX)
                                                       : jumps_code.maxJumpX;
            periodicDiagJy_V = (cfg_.units != nullptr) ? cfg_.units->PotentialCodeToSI(jumps_code.maxJumpY)
                                                       : jumps_code.maxJumpY;
            if (periodicDiagDebug)
            {
                if (!jumps_code.debugX.empty())
                    std::cout << "[ImplicitPIC][periodic][debug][X] step=" << (step + 1)
                              << " " << jumps_code.debugX << std::endl;
                if (!jumps_code.debugY.empty())
                    std::cout << "[ImplicitPIC][periodic][debug][Y] step=" << (step + 1)
                              << " " << jumps_code.debugY << std::endl;
            }
            periodicDiagPendingWrite = true;
        }
        stepTiming.waveformPeriodicSample_s += secondsBetween(tWavePeriodicSample0, Clock::now());

        // Profile1D writers (Turner-style x-binned diagnostics).
        const auto tProfiles0 = Clock::now();
        const bool hasAlphaOut = useCutMetrics && (cutAlphaByCellIndex.size() == static_cast<std::size_t>(mesh.cellIndexNum()));
        const std::vector<double>* alphaPtrOut = hasAlphaOut ? &cutAlphaByCellIndex : nullptr;

        const PicProfiles::Sample profileSample{mesh, cfg_.units.get(), cfg_.normalization.get(),
            phi, EcellNew, species, macroWeight, alphaPtrOut};
        profiles.accumulate(step + 1, profileSample, cachesN);

        if (profiles.instantaneousDue(step + 1))
        {
            auto buildCacheForSpeciesCopy = [&](std::size_t si, CouplingWeightsCache& cache) {
                if (si < cachesN.size())
                    cache = cachesN[si];
            };
            profiles.writeInstantaneous(step + 1, t_code_step, t_si_step, profileSample,
                rhoStar, buildCacheForSpeciesCopy, cfg_.writeSIOutputs);
        }
        trajectories.sample(species, rng, step + 1, t_si_step, cfg_.units->scales().L0_m, cfg_.units->scales().v0_mps);
        stepTiming.profiles_s += secondsBetween(tProfiles0, Clock::now());

        // (9) Gather E^{n+1/2} using orbit-averaged overlaps and permanent push.
        const auto tGatherPush0 = Clock::now();
        particleState.pushAll(mesh, cachesBar, Ehalf, dt, iGeometryPtr);
        if (movingGeometry.controller() && iGeometry)
        {
            for (std::size_t si = 0; si < species.size(); ++si)
            {
                const double wDefault = (si < macroWeight.size()) ? macroWeight[si] : 1.0;
                for (const auto& imp : particleState.pushResult(si).embedded.impacts())
                    movingGeometry.ingestImpact(mesh, *iGeometry, cutMetricsPtr, lineToBoundaryId,
                        imp, static_cast<int>(si), species[si].mass(), species[si].charge(),
                        wDefault, v0_mps, t_si_step);
            }
        }
        stepTiming.gatherPush_s += secondsBetween(tGatherPush0, Clock::now());

        // Domain outflow / flux diagnostics.
        const auto tFluxWave0 = Clock::now();
        {
            auto sideIndexFromDomainSide = [](ExplicitPIC::Pusher::DomainSide s) -> int {
                using DS = ExplicitPIC::Pusher::DomainSide;
                switch (s)
                {
                case DS::kLeft:   return 0;
                case DS::kRight:  return 1;
                case DS::kBottom: return 2;
                case DS::kTop:    return 3;
                default:          return -1;
                }
            };

            auto sideNameFromIndex = [](int idx) -> const char* {
                switch (idx)
                {
                case 0: return "left";
                case 1: return "right";
                case 2: return "bottom";
                case 3: return "top";
                default: return "unknown";
                }
            };

            auto& wallOutMacroBySpecies = wallOutMacroBySpeciesStep;
            auto& wallOutPhysBySpecies = wallOutPhysBySpeciesStep;


            for (std::size_t si = 0; si < species.size(); ++si)
            {
                const double w_default = (si < macroWeight.size()) ? macroWeight[si] : 1.0;
                for (const auto& imp : particleState.pushResult(si).domain.impacts())
                {
                    const int sidx = sideIndexFromDomainSide(imp.side);
                    const auto impact = boundaryLosses.recordDomain(si, species[si], imp.particle, w_default, sidx);
                    if (!impact.accepted) continue;
                    const double w = impact.weight;
                    particleBalance.addDomain(impact.chargeC);

                    if (physicsOutput.includesDomainElectrodes())
                    {
                        const double qC = species[si].charge() * w;
                        const double neutralW = (species[si].charge() == 0.0) ? w : 0.0;
                        physicsOutput.recordElectrode(sidx + 1, qC, neutralW);
                    }
                    if (doEnergyDiag)
                        physicsOutput.observeEnergy(species[si].name(), species[si].charge(), "domain",
                                                 ParticleManager::kineticEnergyEV(imp.particle, species[si].mass(), v0_mps), w);

                }

                for (const auto& imp : particleState.pushResult(si).embedded.impacts())
                {
                    const auto impact = boundaryLosses.recordEmbedded(si, species[si], imp.particle, w_default);
                    const double w = impact.weight;
                    particleBalance.addEmbedded(impact.chargeC);
                    const int boundaryId = (driveWorkEnabled || seeEnabled || doElectrodeDiag || doEnergyDiag || !electrodeGroups.empty())
                                               ? boundaryIdFromPosition(imp.pos, imp.cellLabel)
                                               : 0;

                    boundaryLosses.recordEmbeddedBoundary(boundaryId, impact, imp.crossing, false);
                    if (boundaryId > 0)
                    {
                        const double qC = species[si].charge() * w;
                        const double neutralW = (species[si].charge() == 0.0) ? w : 0.0;
                        physicsOutput.recordElectrode(boundaryId, qC, neutralW);
                    }
                    if (doEnergyDiag)
                        physicsOutput.observeEnergy(species[si].name(), species[si].charge(), "embedded",
                                                 ParticleManager::kineticEnergyEV(imp.particle, species[si].mass(), v0_mps), w);

                    const double impactEnergyOneJ = impact.energyPerParticleJ;

                    emission.processImplicitEmbedded({mesh, iGeometry, dictGeometry, embeddedGeometry},
                        imp, boundaryId, species, si, w, impactEnergyOneJ, seeSpawnPending, rng, unif01);
                }
            }


            wall = boundaryLosses.summarizeImplicitWall(species, mesh, cfg_.units.get(),
                cfg_.normalization.get(), dt_si);

            // Write boundary flux detail (long form) plus wide primary boundary summaries.
            if (reportStep)
            {
                for (std::size_t si = 0; si < species.size(); ++si)
                {
                    for (int sidx = 0; sidx < 4; ++sidx)
                    {
                        const std::uint64_t macroOut = wallOutMacroBySpecies[si][static_cast<std::size_t>(sidx)];
                        const double physOut = wallOutPhysBySpecies[si][static_cast<std::size_t>(sidx)];
                        if (macroOut == 0 && !(physOut > 0.0))
                            continue;

                        const double A = wall.areaM2[sidx];
                        const double dt_out = std::max(1e-300, dt_si);

                        if (speciesSeries.fluxReady())
                        {
                            speciesSeries.writeFlux((step + 1), t_si_step, sideNameFromIndex(sidx), species[si].name(),
                                    PicBoundaryFlux<std::uint64_t>::sample(macroOut, physOut, species[si].charge(), A, dt_out));
                        }
                    }
                }

                if (speciesSeries.boundaryReady())
                {
                    speciesSeries.beginBoundary((step + 1), t_si_step);
                    for (std::size_t si = 0; si < species.size(); ++si)
                    {
                        for (int sidx = 0; sidx < 4; ++sidx)
                        {
                            const double A = wall.areaM2[sidx];
                            const double dt_out = std::max(1e-300, dt_si);
                            const std::uint64_t macroOut = wallOutMacroBySpecies[si][static_cast<std::size_t>(sidx)];
                            const double physOut = wallOutPhysBySpecies[si][static_cast<std::size_t>(sidx)];
                            speciesSeries.addBoundary(PicBoundaryFlux<std::uint64_t>::sample(macroOut, physOut, species[si].charge(), A, dt_out));
                        }
                    }
                    speciesSeries.endBoundary();
                }

                if (speciesSeries.fluxReady())
                    speciesSeries.flushFlux();
            }

            if (speciesSeries.boundaryEnergyReady() && reportStep)
            {
                const double dt_out = std::max(1e-300, dt_si);
                speciesSeries.beginBoundaryEnergy((step + 1), t_si_step);
                for (std::size_t si = 0; si < species.size(); ++si)
                {
                    for (int sidx = 0; sidx < 4; ++sidx)
                    {
                        const long long macroOut = static_cast<long long>(wallOutMacroBySpecies[si][static_cast<std::size_t>(sidx)]);
                        const double physOut = wallOutPhysBySpecies[si][static_cast<std::size_t>(sidx)];
                        const double energyJ = wallAbsEnergyJBySpeciesStep[si][static_cast<std::size_t>(sidx)];

                        const double area_m2 = wall.areaM2[sidx];
                        speciesSeries.addBoundaryEnergy(PicBoundaryEnergy::sample(
                            {macroOut, physOut}, energyJ, area_m2, dt_out));
                    }

                    const long long macroOutEmbedded = ebOutMacroBySpeciesStep[si];
                    const double physOutEmbedded = ebOutPhysBySpeciesStep[si];
                    const double energyEmbeddedJ = ebAbsEnergyJBySpeciesStep[si];

                    const double areaEmbedded_m2 = computeEmbeddedSurfaceAreaM2();
                    speciesSeries.addBoundaryEnergy(PicBoundaryEnergy::sample(
                        {macroOutEmbedded, physOutEmbedded}, energyEmbeddedJ, areaEmbedded_m2, dt_out));
                }
                speciesSeries.endBoundaryEnergy();
            }

            // Conduction current density on selected electrode boundary from wall outflow.
            waveform.sampleConductionSI(species, wallOutPhysBySpecies,
                electrodeArea_m2_forWaveform, dt_si_forWaveform);

            // Boundary electrical work per macro step from Dirichlet electrodes.
            boundaryWork.resetStep(dt_si_forWaveform);
            if (driveWorkEnabled && dt_si_forWaveform > 0.0)
            {
                boundaryWork.addDomain(wall.areaM2, species, wallOutPhysBySpecies);

                boundaryWork.addEmbedded(embeddedChargeOutByBoundaryIdStep_C, fieldBoundaries, *cfg_.units);

                boundaryWork.finishStep();
            }
        }
        stepTiming.fluxWaveform_s += secondsBetween(tFluxWave0, Clock::now());

        spawnManagement.acceptImplicit(spawnStage, seeSpawnPending,
            PicSpawnManagement::ImplicitKind::Secondary, &pmSeeSpawnStats);

        // (10) MCC half-step after the permanent push (Strang splitting).
        const auto tMccPost0 = Clock::now();
        mccStage.advanceImplicitHalf(spawnManagement, spawnStage, 0.5 * dt_si, &pmMccSpawnStats);
        stepTiming.mcc_post_s += secondsBetween(tMccPost0, Clock::now());

        // Publish the selected electrode current every step, independently of CSV cadence.
        waveform.finishStep(dt_si_forWaveform);
        if (waveformSideIndex >= 0 && stepSeries.waveformDue(step + 1))
        {
            const auto& sample = waveform.current();
            stepSeries.writeWaveform(step + 1, t_si_step,
                sample.voltageV, sample.conduction, sample.displacement);
        }

        // Periodic closure diagnostics: write after waveform bookkeeping so the
        // output ordering matches ExplicitPIC while preserving the field-solve
        // snapshot sampled above.
        if (periodicDiagPendingWrite)
        {
            stepSeries.writePeriodic(step + 1, t_si_step,
                periodicDiagJx_V, periodicDiagJy_V, periodicX, periodicY, reportStep);
        }

        if (stepSeries.embeddedBandReady() && outputCadence.periodicSeries(step + 1))
        {
            const double E0 = (cfg_.units != nullptr) ? cfg_.units->scales().E0_Vpm : 1.0;
            stepSeries.writeEmbeddedBand(step + 1, t_si_step,
                summarizePicEmbeddedBandField(mesh, faceGraph, useCutMetrics, cutAlphaByCellIndex,
                    cfg.output.ebBandRMS.bandLayers, EcellNew, E0));
        }

        // Boundary source injection (emission).
        const auto tSource0 = Clock::now();
        if (!particleSources.empty())
        {
            sourceSpawnPending.Clear();

            particleSources.emit({step, dt_si, v0_mps, sc.slabThickness_m},
                mesh, embeddedGeometry, iGeometry, rng, unif01, sourceSpawnPending);

            spawnManagement.acceptImplicit(spawnStage, sourceSpawnPending,
                PicSpawnManagement::ImplicitKind::Source, &pmSourceSpawnStats);
        }
        stepTiming.source_s += secondsBetween(tSource0, Clock::now());

        if (mccOp)
        {
            mccStage.accumulateBatch();
            const int every = (cfg_.collisions ? cfg_.collisions->stats_every : 0);
            if (consoleDebug && every > 0 && ((step + 1) % every == 0))
                std::cout << "[ImplicitPIC][MCC] " << mccStepStats.SummaryString();
        }

        const std::size_t energyCutoffRemovedTotal = ParticleManager::pruneByEnergyCutoff(
            species, cfg.species, macroWeight, v0_mps_energyCutoff,
            &stepParticleMgmtMacroLoss, &stepParticleMgmtPhysLoss);
        if (energyCutoffRemovedTotal > 0 && ConsoleControl::enabled(ConsoleControl::Level::kDebug))
        {
            std::cout << "[ImplicitPIC][EnergyCutoff] step=" << (step + 1)
                      << " removed=" << energyCutoffRemovedTotal << std::endl;
        }

        ParticleManagementStats pmTotals;
        const bool pmApply = pmStep && (pmCfg.mergeEnabled || pmCfg.splitEnabled);
        const auto tPm0 = Clock::now();
        if (pmReport)
        {
            stepSeries.writeSpawnManagement(step + 1, t_si_step, "mcc", pmMccSpawnStats, pmMccSpawnStats.reportKeepProb());
            stepSeries.writeSpawnManagement(step + 1, t_si_step, "see", pmSeeSpawnStats, pmSeeSpawnStats.reportKeepProb());
            stepSeries.writeSpawnManagement(step + 1, t_si_step, "source", pmSourceSpawnStats, pmSourceSpawnStats.reportKeepProb());
        }
        if (pmApply)
        {
            ParticleManager::Step management(pmCfg, mesh, faceGraph, embeddedGeometry, ParticleManager::Step::Mode::Implicit, useCutMetrics);
            ParticleManager::Step::Observer report;
            if (pmReport)
                report = [&](const ParticleType& sp, const ParticleManagementReport& values) {
                    stepSeries.writeParticleManagement(step + 1, t_si_step, sp.name(), values);
                };
            pmTotals = management.applyAll(species, macroWeight, sc.v0_mps, rng, pmReport, report,
                &stepParticleMgmtMacroLoss, &stepParticleMgmtPhysLoss);

            if (pmTotals.splitCount > 0)
                cullParticlesInsideEmbedded("pm", &stepParticleMgmtMacroLoss, &stepParticleMgmtPhysLoss);
        }

        if (pmReport)
            stepSeries.flushParticleManagement();

        stepTiming.pm_s += secondsBetween(tPm0, Clock::now());

        if (pmConsole)
        {
            std::cout << "[ImplicitPIC][PM] step=" << (step + 1)
                      << " total_s=" << std::setprecision(6) << std::fixed << stepTiming.pm_s
                      << " mergePairs=" << pmTotals.mergedPairs
                      << " split=" << pmTotals.splitCount
                      << " splitEmergency=" << pmTotals.splitEmergencyCount
                      << " splitCountBased=" << pmTotals.splitCountBasedCount;
            if (!pmStep)
                std::cout << " (pmStep off)";
            std::cout << std::defaultfloat << std::setprecision(6) << std::endl;
        }

        // Rotate fields.
        fieldHistory.acceptNext();

        // (11) Diagnostics: assembled Gauss residual for var-coeff elliptic + energy summaries.
        const auto tDiag0 = Clock::now();
        {
            const auto particles = summarizeImplicitPicParticles(species, macroWeight, v0_mps,
                doEnergyDiag ? PicBulkEnergyObserver([&](const std::string& name, double charge, double energyEV, double weight) {
                    physicsOutput.observeEnergy(name, charge, "bulk", energyEV, weight);
                }) : PicBulkEnergyObserver{});
            const double Qmacro = particles.macroCharge;
            const double kinetic_e_J = particles.kinetic.electronJ, kinetic_ion_J = particles.kinetic.ionJ,
                         kinetic_neu_J = particles.kinetic.neutralJ;

            particleBalance.finish(Qmacro);

            const auto diag = ExplicitPIC::Diagnostics::PicDiagnostics::sampleImplicit(
                mesh, faceGraph, domainBCUsed, rhoStar, phi, Qstar, Qmacro,
                species, macroWeight, norm, cutMetricsPtr, embeddedBCPtr, embeddedBCByIdPtr,
                kappaBySeg, kappaByCell);
            consoleReport.observeDiagnostics(diag);

            const int stepOut = step + 1;
            const double t_code = stepTime.end;
            const double t_si = t_code * T0_s;

            const double fieldEnergy_J = diag.energyW * Phi0_V;
            const double kineticTotal_J = kinetic_e_J + kinetic_ion_J + kinetic_neu_J;
            const auto energy = energyHistory.sample(dt_si, fieldEnergy_J, kineticTotal_J,
                boundaryWork.current().energyJ, particleBalance.injectionEnergyJ(), dE_abs_domain_J, dE_abs_embedded_J);
            const double totalEnergy_J = energy.totalJ;

            physicsOutput.writeReactions(stepOut, t_si, dt_si, mccStepStats);

            physicsOutput.writeEnergy(stepOut, t_si);

            physicsOutput.writeElectrodes(stepOut, t_si, dt_si);

            if (currentBalanceEnabled && stepSeries.currentBalanceReady() &&
                (outputCadence.series(stepOut)))
            {
                stepSeries.writeCurrentBalance(stepOut, t_si, stepCharge);
            }

            stepSeries.writeGeometryMetrics(stepOut, t_si, scheduledGeometry.current(), /*force=*/reportStep || stepOut == steps);

            if (stepSeries.physicsReady() && outputCadence.series(stepOut))
            {
                const double jCondOut = waveform.current().ready ? waveform.current().conduction : std::numeric_limits<double>::quiet_NaN();
                const double jDispOut = waveform.current().ready ? waveform.current().displacement : std::numeric_limits<double>::quiet_NaN();
                const double jTotOut = waveform.current().ready ? (waveform.current().conduction + waveform.current().displacement)
                                                       : std::numeric_limits<double>::quiet_NaN();

                if (stepSeries.overviewReady())
                {
                    stepSeries.writeOverview(stepOut, t_si, PicStepSeries::ImplicitOverview{
                        diag.rG_L2,
                        diag.rG_Linf,
                        {
                            fieldEnergy_J,
                            kinetic_e_J,
                            kinetic_ion_J,
                            kinetic_neu_J,
                            kineticTotal_J,
                            totalEnergy_J,
                            jCondOut,
                            jDispOut,
                            jTotOut,
                            diag.macro_total,
                            diag.phys_total,
                            stepCharge.relative
                        }
                    });
                }

                if (stepSeries.processReady())
                {
                    const auto processSummary = summarizePicProcesses({
                        jCondOut,
                        jDispOut,
                        boundaryWork.current().energyJ,
                        particleBalance.injectionEnergyJ(),
                        dE_abs_domain_J,
                        dE_abs_embedded_J,
                        mccStepStats.energy_loss_J_total,
                        dt_si,
                        species.size(),
                        wallAbsEnergyJBySpeciesStep,
                        ebAbsEnergyJBySpeciesStep,
                        mccStepStats.per_channel_energy_loss_J
                    });

                    stepSeries.writeProcess(stepOut, t_si, {
                        boundaryWork.current().energyJ,
                        particleBalance.injectionEnergyJ(),
                        mccStepStats.energy_loss_J_total,
                        dE_abs_domain_J,
                        dE_abs_embedded_J,
                        secondarySpawns.energyJ,
                        seeStep.ionImpactEnergyJCathode,
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
                    const auto geometrySummary = summarizePicGeometry(useCutMetrics, cutAlphaByCellIndex,
                        mesh.cellIndexNum(), movingGeometry.controller(), stepOut);

                    stepSeries.writeGeometry(stepOut, t_si, {
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
                        fieldEnergy_J,
                        dE_abs_domain_J,
                        dE_abs_embedded_J,
                        mccStepStats.energy_loss_J_total
                    });
                }

                stepSeries.writePhysics(stepOut, t_si, {diag.rG_L2, diag.rG_Linf,
                    {fieldEnergy_J, kinetic_e_J, kinetic_ion_J, kinetic_neu_J, kineticTotal_J,
                     totalEnergy_J, jCondOut, jDispOut, jTotOut, diag.macro_total, diag.phys_total,
                     stepCharge.relative}}, stepCharge);
            }

            if (speciesSeries.detailReady() && outputCadence.series(stepOut))
            {
                const std::size_t nSpeciesDiag =
                    std::min(diag.speciesMacroCount.size(),
                             std::min(diag.speciesPhysCount.size(), diag.speciesNmean_1pm3.size()));
                speciesSeries.writeDetail(stepOut, t_si, {nSpeciesDiag,
                    diag.speciesMacroCount, diag.speciesPhysCount, diag.speciesNmean_1pm3,
                    diag.electron_mean_energy_eV, diag.Te_e_mean_eV});

                if (speciesSeries.primaryReady())
                {
                    speciesSeries.writePrimary(stepOut, t_si, {nSpeciesDiag,
                    diag.speciesMacroCount, diag.speciesPhysCount, diag.speciesNmean_1pm3,
                    diag.electron_mean_energy_eV, diag.Te_e_mean_eV});
                }

                if (speciesSeries.budgetReady())
                {
                    speciesSeries.beginBudget(stepOut, t_si);
                    for (std::size_t si = 0; si < species.size(); ++si)
                    {
                        const double wDefault = (si < macroWeight.size()) ? macroWeight[si] : 1.0;
                        const auto end = PicSpeciesCount::sample(species[si], wDefault);

                        speciesSeries.addBudget(summarizePicSpeciesBudget({
                            particleBalance.begin(si),
                            sourceSpawns.species[si],
                            mccSpawns.species[si],
                            secondarySpawns.species[si],
                            wallOutMacroBySpeciesStep[si], wallOutPhysBySpeciesStep[si],
                            {ebOutMacroBySpeciesStep[si], ebOutPhysBySpeciesStep[si]},
                            end,
                            {(si < stepParticleMgmtMacroLoss.size()) ? stepParticleMgmtMacroLoss[si] : 0ll,
                             (si < stepParticleMgmtPhysLoss.size()) ? stepParticleMgmtPhysLoss[si] : 0.0}}));
                    }
                    speciesSeries.endBudget();
                }
            }

            if (outputCadence.amounts(stepOut))
            {
                if (stepSeries.amountReady())
                    stepSeries.writeAmount(stepOut, t_si, stepCharge);

                if (stepSeries.amountSpeciesReady())
                {
                    for (std::size_t si = 0; si < species.size(); ++si)
                    {
                        const double nPhys = (si < diag.speciesPhysCount.size()) ? diag.speciesPhysCount[si] : 0.0;
                        stepSeries.writeAmountSpecies(stepOut, t_si, species[si].name(),
                            summarizePicSpeciesAmount(si < Nstar.size() ? &Nstar[si] : nullptr,
                                mesh.cellIndexNum(), nPhys));
                    }
                }
            }

            if (stepSeries.powerReady() && reportStep)
            {
                stepSeries.writePower(stepOut, t_si, {
                    diag.rG_L2, diag.rG_Linf, 0.0, 0.0, 0.0, diag.S_cf, diag.S_eb, 0.0,
                    {diag.S_cf_p50, diag.S_cf_p90, diag.S_cf_p99}, {diag.S_eb_p50, diag.S_eb_p90, diag.S_eb_p99},
                    fieldEnergy_J, kinetic_e_J, kinetic_ion_J, kinetic_neu_J, kineticTotal_J, energy,
                    {boundaryWork.current().energyJ, particleBalance.injectionEnergyJ(), dE_abs_domain_J, dE_abs_embedded_J},
                    boundaryWork.current().potentialCurrentW, boundaryWork.current().potentialCurrentW * dt_si_forWaveform});
            }
        }
        stepTiming.diagnostics_s += secondsBetween(tDiag0, Clock::now());

        // Optional packed field outputs (muCell/muFace).
        const auto tPackedOut0 = Clock::now();
        PicFieldOutputQueue::Timing packedTiming;
        if (outputCadence.fields(step + 1)) packedTiming = packedFieldOutput.writeImplicit(step + 1,
            {faceGraph, embeddedGeometry, iGeometry, dictGeometry.levelSetName, cutMetricsPtr, useCutMetrics,
             phi, rhoStar, Ecell, Nstar, species, macroWeight, cfg_.normalization.get(), coupling,
             muCell_Fpm, kappaByCell, kappaBySeg, periodicX, periodicY});
        stepTiming.packedFieldCompute_s += packedTiming.compute_s;
        stepTiming.packedFieldIo_s += packedTiming.io_s;
        stepTiming.packedFieldOutput_s += secondsBetween(tPackedOut0, Clock::now());

        if (packedFieldOutput.batchingEnabled() && amrEnabled && (amrAdaptEvery > 0) &&
            (((step + 1) % std::max(1, amrAdaptEvery)) == 0))
        {
            const double flushS = packedFieldOutput.flush();
            stepTiming.packedFieldIo_s += flushS;
            stepTiming.packedFieldOutput_s += flushS;
        }

        if (movingGeometry.advance(step + 1, t_si_step, dt_si, mesh, iGeometry))
        {
            cullParticlesInsideEmbedded("movingEB",
                                        &stepParticleMgmtMacroLoss, &stepParticleMgmtPhysLoss);
            movingGeometry.refreshAfterUpdate(step + 1,
                {mesh, *iGeometry, embeddedGeometry, faceGraph, lineToBoundaryId, useEmbeddedGeometry,
                 geometryRefinementActive ? &meshRefiner : nullptr, poisson, &polPoisson});
            if (amrEnabled)
            {
                refinement.requestAdaptation();
            }
        }

        // AMR: update targets and optionally adapt after finishing this macro-step.
        const auto tAmr0 = Clock::now();
        const bool isAmrIndicatorStep =
            amrEnabled && (amrIndicatorEvery > 0) &&
            (((step + 1) % std::max(1, amrIndicatorEvery)) == 0);
        const auto tAmrIndicator0 = Clock::now();
        updateAmrTargetsIfNeeded(step + 1);
        stepTiming.amr_indicator_s += secondsBetween(tAmrIndicator0, Clock::now());
        const bool amrSkippedByGateStep = refinement.skippedByGate();
        const auto tAmrAdapt0 = Clock::now();
        const bool didAmrAdaptStep = doAmrAdaptIfNeeded(step + 1, &stepTiming);
        stepTiming.amr_adapt_s += secondsBetween(tAmrAdapt0, Clock::now());
        if (didAmrAdaptStep)
        {
            poisson.clearCache();
            polPoisson.resetAfterMeshChange(step + 1);
        }
        stepTiming.amr_s += secondsBetween(tAmr0, Clock::now());

        stepTiming.total_s = secondsBetween(tStep0, Clock::now());
        const auto timingSummary = summarizePicImplicitTiming(stepTiming);
        const std::uint64_t nP = consoleReport.writeStep(step, t_si_step, dt_si, summaryStepDefault,
            stepTiming, timingSummary,
            {mesh, phi, cfg_.units.get(), species, macroWeight, stepCharge, boundaryLosses.totals(),
             {wall.electronCurrent, wall.ionCurrent, wall.electronCharge, wall.ionCharge},
             consoleElectrodes, waveform.current(), mccStepStats,
             {isAmrIndicatorStep, amrSkippedByGateStep, didAmrAdaptStep, useEmbeddedGeometry && embeddedGeometry.ready(), mccOp != nullptr,
              static_cast<long long>(sourceSpawnPending.size()), static_cast<long long>(mccSpawnPending.size())}}, std::cout);

        if (stepSeries.couplingFallbackReady())
            stepSeries.writeCouplingFallback(step + 1, t_si_step, stepTiming, timingSummary);
        if (stepSeries.performanceReady() && outputCadence.series(step + 1))
        {
            const bool isOutputStep =
                outputCadence.fields(step + 1);
            const bool isAmrAdaptStep =
                amrEnabled && (amrAdaptEvery > 0) &&
                (((step + 1) % std::max(1, amrAdaptEvery)) == 0) &&
                didAmrAdaptStep;
            const char* polSolverMode = polRep.cgFallbackToDirect ? "cg_fallback"
                                       : (polRep.usedCG ? "cg" : "direct");
            stepSeries.writePerformance(step + 1, t_si_step, stepTiming, timingSummary,
                {nP, isOutputStep, isAmrAdaptStep, amrSkippedByGateStep,
                 polRep.N, polRep.nnz, polSolverMode});
        }

        if (outputCadence.series(step + 1))
        {
            stepSeries.flushPerformance();
            stepSeries.flushPhysics();
            speciesSeries.flushDetail();
            stepSeries.flushAmounts();
            if (speciesSeries.fluxReady())
                speciesSeries.flushFlux();
            if (stepSeries.currentBalanceReady())
                stepSeries.flushCurrentBalance();
            stepSeries.flushEmbeddedBand();
            if (stepSeries.periodicReady())
                stepSeries.flushPeriodic();
            if (stepSeries.powerReady())
                stepSeries.flushPower();
        }


        }
        catch (const std::exception& e)
        {
            std::cerr << "[ImplicitPIC] FATAL: exception at step=" << (step + 1)
                      << " : " << e.what() << std::endl;
            return -3;
        }
        catch (...)
        {
            std::cerr << "[ImplicitPIC] FATAL: unknown exception at step=" << (step + 1) << std::endl;
            return -4;
        }
        return 0;
    };
    const auto completeStep = [&](const PicStepTime& stepTime, bool saveCheckpoint) -> int
    {
        const int step = stepTime.index;
        try
        {
        if (saveCheckpoint)
        {
            const int nextStep = step + 1;
            checkpointWriter.write(nextStep, stepTime.end, stepTime.dt,
                {domain, movingGeometry, species, mccState.get(), &fieldHistory});
        }

        if (consoleDebug && (step % std::max(1, cfg.output.console.every)) == 0)
            std::cout << "[ImplicitPIC] step=" << (step + 1) << "/" << steps << std::endl;
        }
        catch (const std::exception& e)
        {
            std::cerr << "[ImplicitPIC] FATAL: exception at step=" << (step + 1)
                      << " : " << e.what() << std::endl;
            return -3;
        }
        catch (...)
        {
            std::cerr << "[ImplicitPIC] FATAL: unknown exception at step=" << (step + 1) << std::endl;
            return -4;
        }
        return 0;
    };
    const int sessionResult = drivePicStepSession(startStep, steps, dt, advanceStep, driver, completeStep);
    if (sessionResult != 0) return sessionResult;

    (void)packedFieldOutput.flush();

    trajectories.finish();

    return 0;
}
