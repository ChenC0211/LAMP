#include "Modules/Plasma/Boundary/PicScheduledGeometry.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/QuadtreeMeshRefiner.h"
#include "Modules/Plasma/Physics/PicPoissonFieldSolve.h"
#include "Modules/Plasma/ImplicitPIC/ImplicitPicFieldSolve.h"
#include <algorithm>
#include <utility>

PicScheduledGeometry::PicScheduledGeometry(std::unique_ptr<LegacyGridMotionEvolution> motion,
    Bindings bindings, Options options)
    : motion_(std::move(motion)), bindings_(bindings), options_(std::move(options))
{
    if (motion_) current_ = motion_->current();
}

void PicScheduledGeometry::trace(int step, const std::string& message)
{
    if (options_.trace) options_.trace(step, message);
}

bool PicScheduledGeometry::advance(int step, double time_s, bool force)
{
    auto& mesh = bindings_.mesh;
    auto* iGeometry = bindings_.geometry;
    auto* scheduledMotion = motion_.get();
    auto& embeddedGeometry = bindings_.embedded;
    auto& faceGraph = bindings_.graph;
    auto& meshRefiner = bindings_.refiner;
    auto& poisson = bindings_.poisson;
    auto& polPoisson = bindings_.polarization;
    const auto isActivePlasmaCell = [&](const QuadtreeMesh& m, const Cell& c) {
        return embeddedGeometry.activeCell(m, iGeometry, c);
    };
    if (!options_.enabled || !iGeometry || !scheduledMotion)
    {
        current_ = {};
        return false;
    }
    LegacyGridMotionState state = scheduledMotion->prepare(step, force);
    trace(step,
                   "applyScheduledGeometry begin"
                       " force=" + std::to_string(force ? 1 : 0) +
                       " stage=" + state.stageName +
                       " changed=" + std::to_string(state.changed ? 1 : 0));

    if (!force && !state.changed)
    {
        current_ = state;
        return false;
    }

    scheduledMotion->deformReference(iGeometry->lines, state);
    trace(step, "translated boundaries");
    trace(step, "applied accel recession");
    for (std::size_t lsI = 0; lsI < std::min(iGeometry->levelSets.size(), scheduledMotion->innerPointCount()); ++lsI)
        iGeometry->levelSets[lsI].setInnerPoint(scheduledMotion->innerPoint(lsI, state));
    trace(step, "updated level-set inner points");

    // For prescribed rigid shifts, moved polylines are the geometry truth.
    // Rebuild signed distance from the moved lines directly, then refresh
    // derivatives/cut data without PDE redistancing of moving level sets.
    trace(step, "initializeLevelSets begin");
    iGeometry->initializeLevelSets(mesh);
    trace(step, "initializeLevelSets end");
    trace(step, "refreshLevelSetDerivedData begin");
    iGeometry->refreshLevelSetDerivedData(mesh, time_s, /*reinitializeMovingLevelSets=*/false);
    trace(step, "refreshLevelSetDerivedData end");
    if (options_.embeddedEnabled && embeddedGeometry.ready())
    {
        trace(step, "checkCutInfo begin");
        iGeometry->checkCutInfo(mesh, /*isInit=*/true);
        trace(step, "checkCutInfo end");
        trace(step, "rebuildEmbeddedComposite begin");
        embeddedGeometry.rebuildComposite(mesh, iGeometry);
        trace(step, "rebuildEmbeddedComposite end");
    }
    if (options_.refinementActive && meshRefiner.phaseActive())
    {
        trace(step, "markPhases begin");
        iGeometry->markPhases(mesh);
        trace(step, "markPhases end");
    }

    if (embeddedGeometry.ready() && !embeddedGeometry.levelSetIndices().empty())
    {
        trace(step, "faceGraph.build(active) begin");
        faceGraph.build(mesh, isActivePlasmaCell);
        trace(step, "faceGraph.build(active) end");
    }
    else
    {
        trace(step, "faceGraph.build(all) begin");
        faceGraph.build(mesh);
        trace(step, "faceGraph.build(all) end");
    }

    if (options_.embeddedEnabled && embeddedGeometry.ready())
    {
        trace(step, "rebuildEmbeddedMetrics begin");
        embeddedGeometry.rebuildMetrics(mesh, iGeometry, faceGraph, bindings_.lineToBoundaryId, bindings_.moving);
        trace(step, "rebuildEmbeddedMetrics end");
    }

    // A moved cut boundary changes coefficients even when sparsity is unchanged.
    // Invalidate the same solver caches as the MovingEB and AMR paths.
    poisson.clearCache();
    polPoisson.clearDirectCache();
    scheduledMotion->accept(state, iGeometry->lines);
    state = scheduledMotion->current();
    current_ = state;

    // Initialization retains the original refiner setup. A resumed run keeps
    // the cached geometry; only actual time-loop motion refreshes polylines.
    if (options_.refinementActive && !force)
        meshRefiner.refreshBoundaryPolylines();

    trace(step, "applyScheduledGeometry end");
    return true;
}
