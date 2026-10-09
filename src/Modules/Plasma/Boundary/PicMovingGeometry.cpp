#include "Modules/Plasma/Boundary/PicMovingGeometry.h"
#include "Modules/Plasma/ExplicitPIC/Pusher/ParticlePusherLeapfrog.h"
#include "Core/Mesh/QuadtreeMeshRefiner.h"
#include "Modules/Plasma/Physics/PicPoissonFieldSolve.h"
#include "Modules/Plasma/ImplicitPIC/ImplicitPicFieldSolve.h"
#include <sstream>
#include <stdexcept>

bool PicMovingGeometry::initialize(QuadtreeMesh& mesh, Geometry* geometry,
    const DictionaryGeometry* dictionary, const MovingEBConfig& config,
    const UnitsContext* units, const NormalizationContext* normalization,
    const std::optional<std::filesystem::path>& restartDirectory)
{
    if (!dictionary || !geometry || !config.enabled) return true;
    controller_ = std::make_unique<MovingEBController>(*dictionary, config, units, normalization);
    restartEnabled_ = config.restart;
    bool restarted = false;
    if (restartDirectory && restartEnabled_)
    {
        restarted = controller_->readCache(mesh, *restartDirectory);
        if (!restarted)
        {
            const std::filesystem::path cachePath = controller_->resolvedCachePath(*restartDirectory);
            std::cerr << (mode_ == Mode::Explicit ? "[ExplicitPIC]" : "[ImplicitPIC]")
                      << "[Checkpoint] ERROR: failed to restore MovingEB cache: "
                      << cachePath.string() << std::endl;
            return false;
        }
    }
    controller_->initialize(mesh, *geometry, restarted);
    return true;
}

void PicMovingGeometry::beginStep()
{
    if (controller_) controller_->beginStep();
    impactContext_.reset();
}

void PicMovingGeometry::onMeshChanged(const QuadtreeMesh& mesh)
{
    if (controller_) controller_->onMeshChanged(mesh);
    impactContext_.reset();
}

void PicMovingGeometry::ingestImpact(const QuadtreeMesh& mesh, const Geometry& geometry,
    const PicEmbeddedGeometry::Metrics* metrics, const std::vector<int>& lineToBoundaryId,
    const ExplicitPIC::Pusher::EmbeddedImpact& impact, int speciesId,
    double mass, double charge, double defaultWeight, double velocityScaleMps, double timeSI)
{
    if (!controller_) return;
    if (!impactContext_)
        impactContext_ = std::make_unique<MovingEBController::ImpactResolveContext>(
            MovingEBController::buildImpactResolveContext(geometry, lineToBoundaryId));
    WallImpactSample sample;
    MovingEBController::resolveWallImpact(mesh, geometry, metrics, lineToBoundaryId,
        impactContext_.get(), impact.pos, impact.cellLabel, impact.particle,
        speciesId, mass, charge, defaultWeight, velocityScaleMps, timeSI, impact.crossing, sample);
    controller_->ingest(sample);
}

bool PicMovingGeometry::advance(int step, double timeSI, double dtSI,
    QuadtreeMesh& mesh, Geometry* geometry)
{
    impactContext_.reset(); // no lookup may survive geometry mutation
    updated_ = controller_ && geometry
        ? controller_->maybeUpdate(step, timeSI, dtSI, mesh, *geometry) : false;
    return updated_;
}

void PicMovingGeometry::verifyCutMetrics(int step, const Refresh& bindings) const
{
    if (!controller_->lastUpdateHadActiveNormalErosion()) return;
    int invalidGammaBoundary = 0;
    int invalidGammaLine = 0;
    int zeroNormal = 0;
    for (const auto& cm : bindings.embedded.metrics().cellMetricsByIndex())
    {
        if (!cm.isCut) continue;
        if (cm.gammaBoundaryId <= 0) ++invalidGammaBoundary;
        if (cm.gammaLineIndex < 0 || cm.gammaLineIndex >= static_cast<int>(bindings.geometry.lines.size()))
            ++invalidGammaLine;
        if (!(cm.n_gamma.mod2() > 1e-20)) ++zeroNormal;
    }
    if (invalidGammaBoundary > 0 || invalidGammaLine > 0 || zeroNormal > 0)
    {
        std::ostringstream oss;
        oss << "[ImplicitPIC][MovingEB] invalid cut metrics after erosion update at step=" << step
            << " invalidGammaBoundary=" << invalidGammaBoundary
            << " invalidGammaLine=" << invalidGammaLine
            << " zeroNormal=" << zeroNormal;
        throw std::runtime_error(oss.str());
    }
}

void PicMovingGeometry::refreshAfterUpdate(int step, Refresh b)
{
    if (!updated_) return;
    const bool rebuildEmbedded = b.embedded.ready() && (mode_ == Mode::Explicit || b.embeddedEnabled);
    if (rebuildEmbedded) b.embedded.rebuildComposite(b.mesh, &b.geometry);
    if (b.embedded.ready() && !b.embedded.levelSetIndices().empty())
        b.graph.build(b.mesh, [&](const QuadtreeMesh& m, const Cell& c) {
            return b.embedded.activeCell(m, &b.geometry, c);
        });
    else
        b.graph.build(b.mesh);
    if (rebuildEmbedded)
    {
        b.embedded.rebuildMetrics(b.mesh, &b.geometry, b.graph, b.lineToBoundaryId, controller_.get());
        if (mode_ == Mode::Implicit) verifyCutMetrics(step, b);
    }
    if (mode_ == Mode::Implicit && b.refiner)
    {
        // Exact rigid shifts move the reference lines; advected fields leave
        // them stale, so geometry-band AMR must fall back to level-set distance.
        if (controller_->linesInCurrentGeometryFrame()) b.refiner->refreshBoundaryPolylines();
        else b.refiner->clearBoundaryPolylines();
    }
    b.poisson.clearCache();
    if (mode_ == Mode::Implicit && b.polarization) b.polarization->clearDirectCache();
}

void PicMovingGeometry::writeCache(const std::filesystem::path& directory) const
{
    if (controller_ && restartEnabled_) controller_->writeCache(directory);
}
