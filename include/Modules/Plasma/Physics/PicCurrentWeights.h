#pragma once

#include "Modules/Plasma/ExplicitPIC/Deposit/CouplingWeightsCache.h"
#include <cstdint>

class QuadtreeMesh;
class ParticleType;
class PicEmbeddedGeometry;
struct DepositionConfig;

// Current-position coupling, bound to the session mesh. Rebuild is explicit:
// particle motion/reordering and mesh changes invalidate the borrowed caches.
// Policy owns configured coupling support; geometry is sampled at each build.
// Rebuild timing remains explicit at the original integration/output stages.
class PicCurrentWeights
{
public:
    class Policy
    {
    public:
        enum class Mode { Explicit, Implicit };
        Policy(const QuadtreeMesh& mesh, const DepositionConfig& config,
            Mode mode, bool periodicX, bool periodicY);
        Policy(const Policy&) = delete;
        Policy& operator=(const Policy&) = delete;
        // A single minimum-leaf fallback width, also borrowed by the implicit
        // orbit builder. Refresh updates this stable vector in place.
        const std::vector<double>& halfWidths() const { return halfWidths_; }
        double minimumLeafWidth() const { return halfWidths_[0]; }
        void refreshAfterMeshChange();
        // Returned options borrow the current geometry only until its next
        // rebuild. Explicit mode retains its own metrics/size gate; implicit
        // mode uses the caller's original useCutMetrics decision.
        CouplingWeightsCache::BuildOptions physicalOptions(const PicEmbeddedGeometry& geometry,
            bool useCutMetrics = false) const;
        CouplingWeightsCache::BuildOptions diagnosticOptions(const PicEmbeddedGeometry& geometry,
            bool useCutMetrics = false) const;
    private:
        friend class ImplicitPicRefinement;
        friend class ExplicitPicRefinement;
        friend class PicCurrentWeights;
        friend class PicPoissonFieldSolve;
        const QuadtreeMesh& mesh_;
        Mode mode_;
        CouplingWeightsCache::BuildOptions base_;
        std::vector<double> halfWidths_;
    };
    struct Counts { std::uint64_t fallback = 0, particles = 0; };
    explicit PicCurrentWeights(const QuadtreeMesh& mesh) : mesh_(mesh) {}
    PicCurrentWeights(const PicCurrentWeights&) = delete;
    PicCurrentWeights& operator=(const PicCurrentWeights&) = delete;

    Counts rebuild(const std::vector<ParticleType>& species, const std::vector<double>& halfWidths,
        CouplingWeightsCache::BuildOptions options, bool parallelPositionCopy = false);
    Counts rebuild(const std::vector<ParticleType>& species, const Policy& policy,
        const PicEmbeddedGeometry& geometry, bool useCutMetrics = false);
    const std::vector<CouplingWeightsCache>& caches() const { return caches_; }
    void clear() { caches_.clear(); }

    // Independent diagnostic cache; never replaces physical caches. Options
    // are copied and host labels are borrowed only during the synchronous build.
    static void buildDiagnostic(const QuadtreeMesh& mesh, const ParticleType& species,
        const std::vector<double>& halfWidths, CouplingWeightsCache::BuildOptions options,
        CouplingWeightsCache& cache);
    static void buildDiagnostic(const Policy& policy, const PicEmbeddedGeometry& geometry,
        const ParticleType& species, CouplingWeightsCache& cache, bool useCutMetrics = false);

private:
    friend class ImplicitPicRefinement;
    friend class PicPoissonFieldSolve;
    void buildSpecies(const ParticleType& species, const std::vector<double>& halfWidths,
        CouplingWeightsCache::BuildOptions options, CouplingWeightsCache& cache, bool parallelPositionCopy);
    const QuadtreeMesh& mesh_;
    std::vector<CouplingWeightsCache> caches_;
    // Species are built sequentially; one input workspace serves all of them.
    std::vector<Vector2D> positions_;
    std::vector<int> hostLabels_;
};
