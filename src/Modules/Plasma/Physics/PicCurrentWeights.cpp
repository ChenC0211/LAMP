#include "Modules/Plasma/Physics/PicCurrentWeights.h"
#include "Modules/Plasma/Particles/ParticleType.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Modules/Plasma/Boundary/PicEmbeddedGeometry.h"
#include "Modules/Plasma/ExplicitPIC/PicConfig.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

PicCurrentWeights::Policy::Policy(const QuadtreeMesh& mesh, const DepositionConfig& config,
    Mode mode, bool periodicX, bool periodicY)
    : mesh_(mesh), mode_(mode), halfWidths_(1)
{
    base_.periodicX = periodicX;
    base_.periodicY = periodicY;
    base_.enableCompositeBand = config.compositeBand;
    base_.enableMomentExactCut = config.momentExactCut;
    if (mode_ == Mode::Explicit)
    {
        base_.debugParticleId = config.debugParticleId >= 0
            ? static_cast<std::size_t>(config.debugParticleId) : std::numeric_limits<std::size_t>::max();
        base_.dumpWeightsOnce = config.dumpWeightsOnce;
        base_.dumpWeightsPath = config.dumpWeightsPath;
    }
    refreshAfterMeshChange();
}

void PicCurrentWeights::Policy::refreshAfterMeshChange()
{
    double minH = std::numeric_limits<double>::max();
    for (int ci = 0; ci < mesh_.cellIndexNum(); ++ci)
    {
        const int cLbl = mesh_.cellIndexToLabel[ci];
        const Cell& c = mesh_.cells(cLbl);
        if (c.son[0] != -1) continue;
        minH = std::min(minH, std::min(c.size().x, c.size().y));
    }
    if (!std::isfinite(minH) || minH <= 0.0) minH = 1.0;
    halfWidths_[0] = minH;
}

CouplingWeightsCache::BuildOptions PicCurrentWeights::Policy::physicalOptions(
    const PicEmbeddedGeometry& geometry, bool useCutMetrics) const
{
    auto options = base_;
    if (mode_ == Mode::Explicit)
    {
        options.levelSetI = geometry.usesComposite() ? -1 : geometry.levelSetIndex();
        options.levelSetVert = geometry.ready() ? geometry.levelSet() : nullptr;
        options.insideIsNegative = geometry.insideIsNegative();
        options.alphaByCellIndex = geometry.metricsReady() &&
            geometry.alphaByCell().size() == static_cast<std::size_t>(mesh_.cellIndexNum())
            ? &geometry.alphaByCell() : nullptr;
        options.cutByAlpha = geometry.usesComposite();
    }
    else
    {
        if (useCutMetrics)
        {
            options.levelSetI = geometry.levelSetIndex();
            options.levelSetVert = geometry.levelSet();
            options.insideIsNegative = geometry.insideIsNegative();
            options.alphaByCellIndex = &geometry.alphaByCell();
            options.cutByAlpha = true;
        }
        options.recordHpStats = false;
        options.recordPerParticleHp = false;
        options.recordPerParticleCompositeBand = false;
    }
    return options;
}

CouplingWeightsCache::BuildOptions PicCurrentWeights::Policy::diagnosticOptions(
    const PicEmbeddedGeometry& geometry, bool useCutMetrics) const
{
    auto options = base_;
    options.levelSetI = geometry.usesComposite() ? -1 : geometry.levelSetIndex();
    options.levelSetVert = geometry.ready() ? geometry.levelSet() : nullptr;
    options.insideIsNegative = geometry.insideIsNegative();
    const bool haveMetrics = mode_ == Mode::Explicit ? geometry.metricsReady() : useCutMetrics;
    options.alphaByCellIndex = haveMetrics &&
        geometry.alphaByCell().size() == static_cast<std::size_t>(mesh_.cellIndexNum())
        ? &geometry.alphaByCell() : nullptr;
    options.cutByAlpha = geometry.usesComposite();
    if (mode_ == Mode::Explicit)
    {
        options.debugParticleId = std::numeric_limits<std::size_t>::max();
        options.dumpWeightsOnce = false;
        options.dumpWeightsPath.clear();
    }
    return options;
}

PicCurrentWeights::Counts PicCurrentWeights::rebuild(const std::vector<ParticleType>& species,
    const Policy& policy, const PicEmbeddedGeometry& geometry, bool useCutMetrics)
{
    if (&mesh_ != &policy.mesh_) throw std::invalid_argument("PIC coupling policy belongs to a different mesh");
    return rebuild(species, policy.halfWidths(), policy.physicalOptions(geometry, useCutMetrics),
        policy.mode_ == Policy::Mode::Implicit);
}

void PicCurrentWeights::buildDiagnostic(const Policy& policy, const PicEmbeddedGeometry& geometry,
    const ParticleType& species, CouplingWeightsCache& cache, bool useCutMetrics)
{
    buildDiagnostic(policy.mesh_, species, policy.halfWidths(), policy.diagnosticOptions(geometry, useCutMetrics), cache);
}

void PicCurrentWeights::buildSpecies(const ParticleType& species, const std::vector<double>& halfWidths,
    CouplingWeightsCache::BuildOptions options, CouplingWeightsCache& cache, bool parallelPositionCopy)
{
    const int n = static_cast<int>(species.particle.size());
    positions_.resize(static_cast<std::size_t>(n));
    hostLabels_.resize(static_cast<std::size_t>(n));
#ifdef _OPENMP
    const int threshold = std::max(1, CouplingWeightsCache::BuildOptions{}.ompMinParticles);
    #pragma omp parallel for if(parallelPositionCopy && n >= threshold) schedule(static)
#endif
    for (int p = 0; p < n; ++p)
    {
        const auto& particle = species.particle[static_cast<std::size_t>(p)];
        positions_[static_cast<std::size_t>(p)] = particle.position;
        hostLabels_[static_cast<std::size_t>(p)] = particle.inCellL;
    }
    options.hostCellLabels = &hostLabels_;
    cache.build(mesh_, positions_, halfWidths, options);
}

PicCurrentWeights::Counts PicCurrentWeights::rebuild(const std::vector<ParticleType>& species,
    const std::vector<double>& halfWidths, CouplingWeightsCache::BuildOptions options, bool parallelPositionCopy)
{
    caches_.clear();
    caches_.resize(species.size());
    Counts counts;
    for (std::size_t si = 0; si < species.size(); ++si)
    {
        auto& cache = caches_[si];
        buildSpecies(species[si], halfWidths, options, cache, parallelPositionCopy);
        counts.fallback += static_cast<std::uint64_t>(cache.lastBuildFallbackCount());
        counts.particles += static_cast<std::uint64_t>(cache.lastBuildParticleCount());
    }
    return counts;
}

void PicCurrentWeights::buildDiagnostic(const QuadtreeMesh& mesh, const ParticleType& species,
    const std::vector<double>& halfWidths, CouplingWeightsCache::BuildOptions options,
    CouplingWeightsCache& cache)
{
    PicCurrentWeights workspace(mesh);
    workspace.buildSpecies(species, halfWidths, options, cache, false);
}
