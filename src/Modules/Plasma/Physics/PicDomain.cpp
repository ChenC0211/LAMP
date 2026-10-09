#include "Modules/Plasma/Physics/PicDomain.h"
#include "Modules/Plasma/Boundary/MovingEBController.h"
#include "Modules/Plasma/IO/PicCheckpointStore.h"
#include "Modules/Plasma/IO/PicCheckpointReader.h"
#include "Core/Mesh/QuadtreeMeshRefiner.h"
#include "Core/Field/HEdgeField.h"
#include "Core/Field/VEdgeField.h"
#include "Core/IO/ConsoleControl.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace
{
using Clock = std::chrono::steady_clock;
double secondsBetween(Clock::time_point begin, Clock::time_point end)
{
    return std::chrono::duration_cast<std::chrono::duration<double>>(end - begin).count();
}
}
PicDomain::PicDomain(const DictionaryMesh& meshInput,
    const DictionaryMesh& refinementInput, Mode mode)
    : mesh_(meshInput), refinementInput_(refinementInput), mode_(mode),
      embedded_(mode), refiner_(refinementInput_) {}

bool PicDomain::prepareMesh(const PicCheckpointReader* restart,
    std::optional<std::filesystem::path> initialSnapshot)
{
    if (meshPreparationAttempted_ || initializationAttempted_)
        throw std::logic_error("PIC mesh preparation may only be attempted once before geometry initialization");
    if (mode_ == Mode::Explicit && initialSnapshot)
        throw std::invalid_argument("Initial geometry snapshots are only supported by implicit PIC");
    meshPreparationAttempted_ = true;
    initialSnapshot_ = std::move(initialSnapshot);
    if (restart)
    {
        startupRestore_.directory = restart->directory();
        startupRestore_.nextStep = restart->metadata().nextStep;
        if (!restart->readMesh(mesh_, mode_ == Mode::Explicit
            ? PicCheckpointReader::MeshCache::Optional : PicCheckpointReader::MeshCache::Required)) return false;
    }
    else if (initialSnapshot_)
    {
        std::error_code ec;
        const std::filesystem::path meshCachePath = *initialSnapshot_ / "mesh_cache.bin";
        if (std::filesystem::exists(meshCachePath, ec))
        {
            const std::string meshPrefix = PicCheckpointStore::pathPrefix(*initialSnapshot_);
            try
            {
                mesh_.readCache(meshPrefix);
                initialSnapshotLoaded_ = true;
                if (ConsoleControl::enabled(ConsoleControl::Level::kInfo))
                    std::cout << "[ImplicitPIC][init] geometry snapshot loaded: "
                              << initialSnapshot_->string() << std::endl;
            }
            catch (const std::exception& e)
            {
                std::cerr << "[ImplicitPIC][init] WARNING: failed to load geometry snapshot mesh: "
                          << e.what() << std::endl;
                initialSnapshotLoaded_ = false;
            }
        }
        if (initialSnapshotLoaded_) startupRestore_.directory = *initialSnapshot_;
    }
    meshPrepared_ = true;
    return true;
}

std::optional<PicDomain::InitialState> PicDomain::initializeForRun(DictionaryGeometry* input,
    bool useEmbeddedGeometry, bool useGeometryRefinement, bool geometryChangesDuringRun,
    const std::vector<int>& lineToBoundaryId, const LegacyGridMotionConfig* motionSchedule)
{
    if (!meshPrepared_ || initializationAttempted_)
        throw std::logic_error("PIC run geometry requires a prepared mesh and may only be initialized once");
    const bool implicit = mode_ == Mode::Implicit;
    const bool consoleInfo = ConsoleControl::enabled(ConsoleControl::Level::kInfo);
    const bool consoleDebug = !implicit && ConsoleControl::enabled(ConsoleControl::Level::kDebug);
    const bool needGeometry = useEmbeddedGeometry || useGeometryRefinement;
    const bool restartReady = startupRestore_.nextStep.has_value();
    const char* prefix = implicit ? "[ImplicitPIC]" : "[ExplicitPIC]";
    const auto logInitStage = [&](const char* stage, double seconds) {
        if (!implicit || !consoleInfo) return;
        std::cout << "[ImplicitPIC][init] " << stage
                  << " done in " << std::fixed << std::setprecision(3) << seconds << " s"
                  << std::defaultfloat << std::setprecision(6) << std::endl;
    };
    if (implicit && consoleInfo && input && needGeometry)
        std::cout << "[ImplicitPIC][init] preparing embedded geometry..." << std::endl;
    {
        const auto begin = Clock::now();
        if (!initializeGeometry(input, useEmbeddedGeometry, useGeometryRefinement, startupRestore_, motionSchedule))
            return std::nullopt;
        if (input && needGeometry) logInitStage("geometry state setup", secondsBetween(begin, Clock::now()));
    }
    bool geometryRefinementActive = false;
    const auto attachRefinement = [&]() {
        if (useGeometryRefinement)
        {
            geometryRefinementActive = refiner_.attachGeometry(input, geometry_.get());
            if (geometryRefinementActive && refiner_.phaseActive() && geometry_)
                geometry_->markPhases(mesh_);
            if (!geometryRefinementActive) useGeometryRefinement = false;
        }
    };
    // Explicit attaches phase refinement before the composite; implicit does
    // so afterwards. Preserve both orders, including startup diagnostics.
    if (!implicit) attachRefinement();
    if (embedded_.ready())
    {
        const auto begin = Clock::now();
        embedded_.rebuildComposite(mesh_, geometry_.get());
        logInitStage("embedded composite", secondsBetween(begin, Clock::now()));
        if (useEmbeddedGeometry && consoleInfo)
        {
            std::ostringstream boundaryList, levelSetList;
            boundaryList << "{";
            for (std::size_t i = 0; i < embedded_.boundaryIds().size(); ++i)
            {
                if (i > 0) boundaryList << ",";
                boundaryList << embedded_.boundaryIds()[i];
            }
            boundaryList << "}";
            levelSetList << "{";
            for (std::size_t i = 0; i < embedded_.levelSetIndices().size(); ++i)
            {
                if (i > 0) levelSetList << ",";
                levelSetList << (embedded_.levelSetIndices()[i] + 1);
            }
            levelSetList << "}";
            std::cout << prefix << "[EB] resolved boundaryIds=" << boundaryList.str()
                      << " levelSetIndices(1-based)=" << levelSetList.str()
                      << " insideIsNegative=" << (embedded_.insideIsNegative() ? "true" : "false")
                      << " absorbModel=BC/Particle(default=absorb)"
                      << " absorbCrossing=false alphaSampleN=" << PicEmbeddedGeometry::alphaSampleCount
                      << std::endl;
        }
    }
    if (implicit) attachRefinement();
    if ((!implicit || (!restartReady && !initialSnapshotLoaded_)) &&
        geometryRefinementActive && geometry_ && mesh_.maxInitialIterations() > 0)
        refineInitial(input, nullptr, useEmbeddedGeometry, refiner_.phaseActive());

    if (implicit && !restartReady && initialSnapshot_ && !initialSnapshotLoaded_ && geometry_)
    {
        std::error_code ec;
        std::filesystem::create_directories(*initialSnapshot_, ec);
        if (!ec)
        {
            const std::string snapshotPrefix = PicCheckpointStore::pathPrefix(*initialSnapshot_);
            try
            {
                mesh_.writeCache(snapshotPrefix);
                writeGeometryCache(snapshotPrefix);
                if (consoleInfo)
                    std::cout << "[ImplicitPIC][init] geometry snapshot saved: "
                              << initialSnapshot_->string() << std::endl;
            }
            catch (const std::exception& e)
            {
                std::cerr << "[ImplicitPIC][init] WARNING: failed to save geometry snapshot: "
                          << e.what() << std::endl;
            }
        }
    }
    if (geometryRefinementActive && !geometryChangesDuringRun)
    {
        if (consoleInfo)
            std::cout << prefix << "[AMR] static geometry detected; runtime geometry-driven refinement is disabled after startup mesh setup"
                      << std::endl;
        refiner_.clearGeometry();
        geometryRefinementActive = false;
    }
    {
        const auto begin = Clock::now();
        if (consoleDebug) std::cerr << "[FaceSegmentGraph] build begin.\n";
        rebuildFaceGraph();
        if (consoleDebug) std::cerr << "[FaceSegmentGraph] build done.\n";
        logInitStage("face graph", secondsBetween(begin, Clock::now()));
    }
    if ((!implicit || useEmbeddedGeometry) && embedded_.ready())
    {
        const auto begin = Clock::now();
        if (consoleDebug) std::cerr << "[CutCellMetrics] build begin.\n";
        embedded_.rebuildMetrics(mesh_, geometry_.get(), faceGraph_, lineToBoundaryId, nullptr);
        if (consoleDebug) std::cerr << "[CutCellMetrics] build done.\n";
        logInitStage("cut metrics", secondsBetween(begin, Clock::now()));
    }
    return InitialState{useEmbeddedGeometry, geometryRefinementActive};
}

void PicDomain::rebuildFaceGraph()
{
    if (embedded_.ready() && !embedded_.levelSetIndices().empty())
        faceGraph_.build(mesh_, [this](const QuadtreeMesh& mesh, const Cell& cell) {
            return embedded_.activeCell(mesh, geometry_.get(), cell);
        });
    else
        faceGraph_.build(mesh_);
}

std::unique_ptr<LegacyGridMotionEvolution> PicDomain::takeScheduledMotion()
{
    return std::move(scheduledMotion_);
}

void PicDomain::writeGeometryCache(const std::string& prefix) const
{
    if (geometry_) geometry_->wirteCache(prefix);
}

bool PicDomain::initializeGeometry(DictionaryGeometry* input,
    bool& useEmbeddedGeometry, bool& useGeometryRefinement, const Restore& restore,
    const LegacyGridMotionConfig* motionSchedule)
{
    if (initializationAttempted_) throw std::logic_error("PIC geometry initialization may only be attempted once");
    initializationAttempted_ = true;
    auto& mesh = mesh_;
    const auto mode = mode_;
    auto& embeddedGeometry = embedded_;
    auto& iGeometry = geometry_;
    auto& scheduledMotion = scheduledMotion_;
    const bool needGeometry = useEmbeddedGeometry || useGeometryRefinement;
    const bool haveGeometry = input != nullptr;
    const char* prefix = mode == Mode::Explicit ? "[ExplicitPIC]" : "[ImplicitPIC]";
    const bool consoleDebug = mode == Mode::Explicit && ConsoleControl::enabled(ConsoleControl::Level::kDebug);
    if (needGeometry)
    {
        if (!haveGeometry)
        {
            if (useEmbeddedGeometry)
                std::cerr << prefix << "[EB] WARNING: EmbeddedGeometry enabled but geometry.dict could not be read. Disabling EB.\n";
            if (useGeometryRefinement)
                std::cerr << prefix << "[AMR] WARNING: GeometryRefinement enabled but geometry.dict could not be read. Disabling geometry refinement.\n";
            useEmbeddedGeometry = false;
            useGeometryRefinement = false;
        }
    }
    if (needGeometry && haveGeometry)
    {
        auto& dictGeometry = *input;
        iGeometry = std::make_unique<ImmersedGeometry<QuadtreeMesh>>(dictGeometry, mesh);
        if (consoleDebug)
            std::cerr << "[ImmersedGeometry] initializeLevelSets begin.\n";
        iGeometry->initializeLevelSets(mesh);
        if (consoleDebug)
            std::cerr << "[ImmersedGeometry] initializeLevelSets done.\n";
        if (consoleDebug)
            std::cerr << "[ImmersedGeometry] reinitializeLevelSets begin.\n";
        if (motionSchedule && motionSchedule->enabled && !iGeometry->levelSets.empty())
        {
            LegacyGridMotionEvolution::Reference reference;
            reference.lines = iGeometry->lines;
            for (const auto& levelSet : iGeometry->levelSets)
                reference.innerPoints.push_back(levelSet.innerPoint());
            auto boundaryLines = [&](int id) -> std::vector<int> {
                const auto found = dictGeometry.boundaryById.find(id);
                return found == dictGeometry.boundaryById.end() ? std::vector<int>{} : found->second.lineIds;
            };
            LegacyGridMotion::Bindings lines{
                boundaryLines(10), boundaryLines(11), boundaryLines(12),
                boundaryLines(13), boundaryLines(14), boundaryLines(15)};
            LegacyGridMotionEvolution::InnerPointBindings innerPoints;
            if (reference.innerPoints.size() >= 6u)
                innerPoints = {{2, 3}, {4, 5}};
            scheduledMotion = std::make_unique<LegacyGridMotionEvolution>(
                *motionSchedule, std::move(lines), std::move(reference), std::move(innerPoints));
        }
        if (restore.directory.has_value())
        {
            const std::string prefix = PicCheckpointStore::pathPrefix(*restore.directory);
            const std::filesystem::path geomCachePath = (*restore.directory) / "iGeometry_cache.bin";
            std::error_code geomEc;
            if (std::filesystem::exists(geomCachePath, geomEc))
                iGeometry->readCache(mesh, prefix);
            else
                iGeometry->refreshLevelSetDerivedData(mesh, 0.0, /*reinitializeMovingLevelSets=*/false);
        }
        else
        {
            iGeometry->refreshLevelSetDerivedData(mesh, 0.0, /*reinitializeMovingLevelSets=*/false);
        }
        if (consoleDebug)
            std::cerr << "[ImmersedGeometry] reinitializeLevelSets done.\n";
        if (restore.nextStep.has_value() && scheduledMotion)
        {
            try
            {
                scheduledMotion->restore(std::max(0, *restore.nextStep), iGeometry->lines);
                // Inner points are input metadata, not part of the geometry cache.
                for (std::size_t lsI = 0; lsI < std::min(iGeometry->levelSets.size(), scheduledMotion->innerPointCount()); ++lsI)
                    iGeometry->levelSets[lsI].setInnerPoint(scheduledMotion->innerPoint(lsI, scheduledMotion->current()));
            }
            catch (const std::exception& error)
            {
                std::cerr << "[ImplicitPIC][MotionSchedule] ERROR: " << error.what() << '\n';
                return false;
            }
        }
        if (useEmbeddedGeometry)
        {
            if (consoleDebug)
                std::cerr << "[ImmersedGeometry] checkCutInfo begin.\n";
            iGeometry->checkCutInfo(mesh, /*isInit=*/true);
            if (consoleDebug)
                std::cerr << "[ImmersedGeometry] checkCutInfo done.\n";
        }

        const int nLS = static_cast<int>(iGeometry->levelSets.size());
        if (nLS <= 0 && mode == Mode::Explicit)
        {
            if (useEmbeddedGeometry)
                std::cerr << prefix << "[EB] WARNING: geometry.dict defines no level sets. Disabling EB.\n";
            if (useGeometryRefinement)
                std::cerr << prefix << "[AMR] WARNING: geometry.dict defines no level sets. Disabling geometry refinement.\n";
            useEmbeddedGeometry = false;
            useGeometryRefinement = false;
            iGeometry.reset();
        }
        else if (useEmbeddedGeometry && nLS > 0)
        {
            embeddedGeometry.resolveSelection(dictGeometry, nLS);
            if (mode == Mode::Explicit && !embeddedGeometry.ready()) useEmbeddedGeometry = false;
        }
    }

    return true;
}

void PicDomain::refineInitial(const DictionaryGeometry* input, MovingEBController* movingEb,
    bool useEmbeddedGeometry, bool phaseActive)
{
    if (!geometry_) throw std::logic_error("PIC startup refinement requires initialized geometry");
    auto& mesh = mesh_;
    const auto mode = mode_;
    auto& embeddedGeometry = embedded_;
    Geometry* iGeometry = geometry_.get();
    const bool verbose = ConsoleControl::enabled(mode == Mode::Explicit
        ? ConsoleControl::Level::kDebug : ConsoleControl::Level::kInfo);
    const char* prefix = mode == Mode::Explicit ? "[ExplicitPIC]" : "[ImplicitPIC]";
    const auto tInitRefine0 = Clock::now();
    std::ofstream initRefineCsv;
    if (mode == Mode::Implicit)
    {
        const std::filesystem::path initRefineCsvPath =
            std::filesystem::path("./output") / "meta" / "init_geometry_refine.csv";
        std::filesystem::create_directories(initRefineCsvPath.parent_path());
        initRefineCsv.open(initRefineCsvPath.string(), std::ios::out | std::ios::trunc);
        if (mode == Mode::Implicit && initRefineCsv)
        {
            initRefineCsv << "iter,cells_before,cells_after,build_target_s,adapt_s,movingeb_onMeshChanged_s,"
                          << "post_check_s,post_converged,reinit_levelset_s,cutinfo_composite_s,mark_phases_s,iter_total_s\n";
            initRefineCsv << std::setprecision(16) << std::fixed;
            initRefineCsv.flush();
        }
    }
    auto initTargetHasPendingChange = [&](const CellField<double>& target) -> bool {
        const int nCell = mesh.cellIndexNum();
        if (target.size() != nCell)
            return true;
        for (int ci = 0; ci < nCell; ++ci)
        {
            const int cLbl = mesh.cellIndexToLabel[ci];
            const Cell& c = mesh.cells(cLbl);
            const double v = target.value[ci];
            if (!std::isfinite(v))
                return true;
            const int aimLevel = static_cast<int>(std::lround(v));
            if (aimLevel != c.level())
                return true;
        }
        return false;
    };
    DictionaryMesh meshRefinerGeoDict = refinementInput_;
    meshRefinerGeoDict.physicsRefinement.enabled = false;
    meshRefinerGeoDict.physicsRefinement.indicators.clear();
    QuadtreeMeshRefiner meshRefinerGeo(meshRefinerGeoDict);
    if (input)
        meshRefinerGeo.attachGeometry(input, iGeometry);

    const int initIters = std::max(0, mesh.maxInitialIterations());
    CellField<double> preTarget(mesh, 0.0);
    CellField<double> preIndicator(mesh, 0.0);
    CellField<double> postTarget(mesh, 0.0);
    CellField<double> postIndicator(mesh, 0.0);
    for (int it = 0; it < initIters; ++it)
    {
        const auto tPreIter0 = Clock::now();
        if (mode == Mode::Explicit && verbose)
            std::cerr << "[ExplicitPIC][AMR] pre-iter " << (it + 1) << "/" << initIters << std::endl;
        const int cellsBefore = mesh.cellIndexNum();
        const auto tBuildTarget0 = Clock::now();
        meshRefinerGeo.buildTargetLevels(mesh, preTarget, &preIndicator, nullptr);
        const double buildTargetS = secondsBetween(tBuildTarget0, Clock::now());
        if (!initTargetHasPendingChange(preTarget))
        {
            if (verbose)
            {
                std::cout << prefix << "[init] initial geometry refinement converged at iter "
                          << (it + 1) << "/" << initIters
                          << " in " << std::fixed << std::setprecision(3)
                          << secondsBetween(tPreIter0, Clock::now()) << " s"
                          << std::defaultfloat << std::setprecision(6) << std::endl;
            }
            break;
        }

        std::vector<VertField<double>> levelSets;
        std::vector<VertField<double>*> vertFields;
        std::vector<HEdgeField<double>*> hEdgeFields;
        std::vector<VEdgeField<double>*> vEdgeFields;
        std::vector<CellField<double>*> cellFields;
        cellFields.push_back(&preTarget);
        for (auto& ls : iGeometry->levelSets)
            vertFields.push_back(&ls);

        const auto tAdapt0 = Clock::now();
        mesh.adapt(levelSets, preIndicator, vertFields, hEdgeFields, vEdgeFields, cellFields, &preTarget);
        const double adaptS = secondsBetween(tAdapt0, Clock::now());
        double movingEbMeshChangedS = 0.0;
        if (movingEb)
        {
            const auto tMovingEb0 = Clock::now();
            movingEb->onMeshChanged(mesh);
            movingEbMeshChangedS = secondsBetween(tMovingEb0, Clock::now());
        }

        const auto tReinit0 = Clock::now();
        iGeometry->initializeLevelSets(mesh);
        iGeometry->refreshLevelSetDerivedData(mesh, 0.0, /*reinitializeMovingLevelSets=*/false);
        const double reinitLevelsetS = secondsBetween(tReinit0, Clock::now());

        double cutInfoCompositeS = 0.0;
        if ((mode == Mode::Explicit || useEmbeddedGeometry) && embeddedGeometry.ready())
        {
            const auto tCutInfoComposite0 = Clock::now();
            iGeometry->checkCutInfo(mesh, /*isInit=*/true);
            embeddedGeometry.rebuildComposite(mesh, iGeometry);
            cutInfoCompositeS = secondsBetween(tCutInfoComposite0, Clock::now());
        }

        double markPhasesS = 0.0;
        if (phaseActive)
        {
            const auto tMarkPhases0 = Clock::now();
            iGeometry->markPhases(mesh);
            markPhasesS = secondsBetween(tMarkPhases0, Clock::now());
        }

        const auto tPostCheck0 = Clock::now();
        meshRefinerGeo.buildTargetLevels(mesh, postTarget, &postIndicator, nullptr);
        const bool postConverged = !initTargetHasPendingChange(postTarget);
        const double postCheckS = secondsBetween(tPostCheck0, Clock::now());

        const double iterTotalS = secondsBetween(tPreIter0, Clock::now());
        if (mode == Mode::Implicit && initRefineCsv)
        {
            initRefineCsv << (it + 1) << ','
                          << cellsBefore << ','
                          << mesh.cellIndexNum() << ','
                          << buildTargetS << ','
                          << adaptS << ','
                          << movingEbMeshChangedS << ','
                          << postCheckS << ','
                          << (postConverged ? 1 : 0) << ','
                          << reinitLevelsetS << ','
                          << cutInfoCompositeS << ','
                          << markPhasesS << ','
                          << iterTotalS
                          << "\n";
            initRefineCsv.flush();
        }

        if (verbose)
        {
            std::cout << prefix << "[init] geometry refine iter "
                      << (it + 1) << "/" << initIters
                      << " done in " << std::fixed << std::setprecision(3)
                      << iterTotalS << " s"
                      << " cells=" << mesh.cellIndexNum()
                      << std::defaultfloat << std::setprecision(6) << std::endl;
        }

        if (postConverged)
        {
            if (verbose)
            {
                std::cout << prefix << "[init] initial geometry refinement converged after adapt at iter "
                          << (it + 1) << "/" << initIters
                          << std::endl;
            }
            break;
        }
    }
    if (mode == Mode::Implicit && verbose)
    {
        std::cout << "[ImplicitPIC][init] initial geometry refinement"
                  << " done in " << std::fixed << std::setprecision(3)
                  << secondsBetween(tInitRefine0, Clock::now()) << " s"
                  << std::defaultfloat << std::setprecision(6) << std::endl;
    }
}
