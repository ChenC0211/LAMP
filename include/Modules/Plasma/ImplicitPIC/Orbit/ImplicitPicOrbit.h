#pragma once

#include "Modules/Plasma/ImplicitPIC/Orbit/OrbitAveragedWeightsBuilder.h"
#include "Modules/Plasma/ImplicitPIC/Orbit/PredictorOrbitSegmenter.h"
#include <cstdint>

class QuadtreeMesh;
class ParticleType;
class PicEmbeddedGeometry;
class UnitsContext;
struct PicConfig;
struct ImplicitPicOverrides;
template<class T> class CellField;

// Owns E^n gathers, predictor endpoints/segments and orbit-averaged caches.
// Build borrows immutable particles/geometry for the call; it never pushes
// particles. Outputs are invalidated by the next build/reset or mesh change.
// On failure, the caller stops the step; partial buffers are not rolled back.
class ImplicitPicOrbit
{
public:
    struct Options
    {
        ImplicitPIC::Orbit::PredictorOrbitSegmenter::Options segmenter;
        CouplingWeightsCache::BuildOptions weights;
        ImplicitPIC::Orbit::OrbitAveragedWeightsBuilderOptions averaging;
        bool parallelSpecies = false;
    };
    struct Trace { int step = 0; bool info = false, debug = false; };
    struct SpeciesPerf
    {
        double predictorPrep_s = 0.0;
        double predictorSegment_s = 0.0;
        std::uint64_t predictor_segments_total = 0u;
        std::uint64_t predictor_fast_exit_hits = 0u;
        std::uint64_t predictor_bisection_fallback_hits = 0u;
        std::uint64_t predictor_coarsen_hits = 0u;
        std::uint64_t predictor_coarsen_parent_climb = 0u;

        double orbitAvgWeights_s = 0.0;
        int orbitAvg_fused_mode = 0;
        double orbitAvg_fused_eval_s = 0.0;
        double orbitAvg_fused_accum_s = 0.0;
        double orbitAvg_sampleFill_s = 0.0;
        double orbitAvg_sampleCacheBuild_s = 0.0;
        double orbitAvg_localReduce_s = 0.0;
        double orbitAvg_concat_s = 0.0;
        std::uint64_t orbitAvg_sampleCwFallbackCount = 0u;
        std::uint64_t orbitAvg_sampleCount = 0u;
        std::uint64_t orbitAvg_sampleCache_nnz = 0u;
        std::uint64_t orbitAvg_out_nnz = 0u;
        std::uint64_t orbitAvg_adaptAcceptCount = 0u;
        std::uint64_t orbitAvg_adaptRefineCount = 0u;
        std::uint64_t orbitAvg_adaptHitMaxDepthCount = 0u;
        std::uint64_t orbitAvg_adaptAcceptedDepthSum = 0u;
    };


    ImplicitPicOrbit(const QuadtreeMesh& mesh, std::size_t speciesCount);
    // Snapshot configured integration/segmentation controls. Geometry is bound
    // for each build, so no pointers into rebuilt cut metrics are retained.
    ImplicitPicOrbit(const QuadtreeMesh& mesh, std::size_t speciesCount,
        const PicConfig& config, const ImplicitPicOverrides& overrides,
        const UnitsContext* units, bool periodicX, bool periodicY);
    ImplicitPicOrbit(const ImplicitPicOrbit&) = delete;
    ImplicitPicOrbit& operator=(const ImplicitPicOrbit&) = delete;

    static void validateCurrentCaches(const std::vector<ParticleType>& species,
        const std::vector<CouplingWeightsCache>& caches);
    void gather(const std::vector<ParticleType>& species,
        const std::vector<CouplingWeightsCache>& cachesN, const CellField<Vector2D>& Ecell, int step);
    std::vector<SpeciesPerf> build(const std::vector<ParticleType>& species,
        const std::vector<double>& hp, double dt, const Options& options, Trace trace);
    std::vector<SpeciesPerf> build(const std::vector<ParticleType>& species,
        const std::vector<double>& hp, double dt, const PicEmbeddedGeometry& geometry,
        bool useCutMetrics, bool enableEmbeddedIntercept, Trace trace);
    // Returned options borrow geometry until its next rebuild. The caller's
    // original cut/interception gates are explicit, independent of readiness.
    Options optionsFor(const PicEmbeddedGeometry& geometry,
        bool useCutMetrics, bool enableEmbeddedIntercept) const;
    void resetScratch();
    const std::vector<CouplingWeightsCache>& averagedCaches() const { return cachesBar; }

private:
    friend class ImplicitPicRefinement;
    const QuadtreeMesh& mesh;
    Options configuredOptions_;
    std::vector<CouplingWeightsCache> cachesBar;
    std::vector<ImplicitPIC::Orbit::OrbitAveragedWeightsBuilderWorkspace> orbitAvgWorkspaces;
    std::vector<std::vector<Vector2D>> gatherScratch, predictorX0Scratch, predictorX1Scratch;
    std::vector<std::vector<int>> predictorHostScratch, segOffsetsScratch;
    std::vector<std::vector<ImplicitPIC::Orbit::OrbitSegment>> segsScratch;
};
