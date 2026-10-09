#pragma once

#include <cstddef>
#include <string>
#include <vector>

class BackgroundMesh;
class CouplingWeightsCache;

#include "Modules/Plasma/ExplicitPIC/Deposit/CouplingWeightsCache.h"

namespace ImplicitPIC
{
namespace Orbit
{

struct OrbitSegment;

struct OrbitAveragedWeightsBuilderOptions
{
    bool deterministic = true;

    // If true, renormalize per particle to enforce sum Wbar = 1.
    bool renormalize = true;

    // If true, evaluate/accumulate per-sample weights directly into per-particle
    // accumulators (avoid global sample-cache CSR materialization).
    bool useFusedPath = false;

    // -----------------------------------------------------------------
    // Adaptive time quadrature (whole-step integration)
    // -----------------------------------------------------------------
    bool adaptiveEnabled = false;
    double adaptiveTol = 1e-3;
    int adaptiveMaxRefineDepth = 8;
    // Currently supported: "gauss2_4".
    std::string adaptiveRule = "gauss2_4";

    // If true, use a cheap refinement criterion based on host-cell label
    // changes ("补刀" mode). When enabled, accepted leaf intervals are
    // integrated with Gauss2 (no Gauss4 error estimator).
    bool adaptiveCheapRefine = false;

    // Optional extra split at max depth for cheapRefine crossing intervals.
    bool adaptiveCrossSplitAtMaxDepth = false;
    int adaptiveCrossSplitMaxIters = 12;
    // Sheath gate in code length units.
    double adaptiveCrossSplitSheathPad = 0.0;
    // Per-particle limit for max-depth extra splits.
    int adaptiveCrossSplitMaxSplits = 1;

    // -----------------------------------------------------------------
    // Hybrid mode: segment-based integration with adaptive "patch" only
    // on segments that appear to cross leaf cells.
    // -----------------------------------------------------------------
    bool segmentPatchEnabled = false;
    double segmentPatchTol = 1e-3;
    int segmentPatchMaxRefineDepth = 8;
};

struct OrbitAveragedWeightsBuilderTiming
{
    int fusedMode = 0;
    double fusedEval_s = 0.0;
    double fusedAccum_s = 0.0;
    double sampleFill_s = 0.0;
    double sampleCacheBuild_s = 0.0;
    double localReduce_s = 0.0;
    double concat_s = 0.0;
    double total_s = 0.0;
    std::size_t sampleCount = 0u;
    std::size_t sampleCacheFallbackCount = 0u;
    std::size_t sampleCacheNnz = 0u;
    std::size_t outNnz = 0u;

    // Adaptive diagnostics (counts are totals across all particles).
    std::size_t adaptAcceptCount = 0u;
    std::size_t adaptRefineCount = 0u;
    std::size_t adaptHitMaxDepthCount = 0u;

    // Sum of depths for intervals that were accepted (including hitMaxDepth accept).
    // This enables reporting an average accepted depth.
    std::size_t adaptAcceptedDepthSum = 0u;
};

struct OrbitAveragedWeightsBuilderWorkspace
{
    // Sample-stage scratch (reused across macro-steps).
    std::vector<std::size_t> sampleOffsetsByParticle;
    std::vector<Vector2D> samplePos;
    std::vector<int> sampleHostLabel;
    std::vector<double> sampleCoeff;
    CouplingWeightsCache sampleCache;

    // Reduce-stage scratch.
    std::vector<int> particleNnz;
    std::vector<int> outOffsets;
    std::vector<std::vector<int>> reducedCellByParticle;
    std::vector<std::vector<double>> reducedWByParticle;

    // Parallel reduce/concat scratch (chunked by thread).
    std::vector<std::vector<int>> reducedCellByThread;
    std::vector<std::vector<double>> reducedWByThread;
    std::vector<std::vector<int>> reducedOffsetsByThread;
    std::vector<std::size_t> reducedThreadBegin;
    std::vector<std::size_t> reducedThreadEnd;
    std::vector<std::size_t> reducedChunkNnzHint;
};

// Build orbit-averaged overlaps Wbar in CouplingWeightsCache format by integrating
// instantaneous CouplingWeightsCache weights along each particle's segmented orbit.
//
// Quadrature: 2-point Gauss per segment, and the average is normalized by dt.
// Output CSR uses cell *indices* (compatible with CellField storage).
bool buildOrbitAveragedCacheGauss2(const BackgroundMesh& mesh,
                                  const std::vector<double>& hp,
                                  const CouplingWeightsCache::BuildOptions& buildOpt,
                                  std::size_t numParticles,
                                  const std::vector<int>& segmentOffsets,
                                  const std::vector<OrbitSegment>& segments,
                                  const std::vector<int>& hostCellLabelByParticle,
                                  double dt,
                                  CouplingWeightsCache& out,
                                  OrbitAveragedWeightsBuilderWorkspace* workspace,
                                  const OrbitAveragedWeightsBuilderOptions& opt = {},
                                  OrbitAveragedWeightsBuilderTiming* timing = nullptr);

bool buildOrbitAveragedCacheGauss2(const BackgroundMesh& mesh,
                                  const std::vector<double>& hp,
                                  const CouplingWeightsCache::BuildOptions& buildOpt,
                                  std::size_t numParticles,
                                  const std::vector<int>& segmentOffsets,
                                  const std::vector<OrbitSegment>& segments,
                                  const std::vector<int>& hostCellLabelByParticle,
                                  double dt,
                                  CouplingWeightsCache& out,
                                  const OrbitAveragedWeightsBuilderOptions& opt = {},
                                  OrbitAveragedWeightsBuilderTiming* timing = nullptr);

// Build orbit-averaged overlaps Wbar by adaptive time quadrature over the whole
// predictor interval [t_n, t_{n+1}]. The orbit is assumed linear from x0 to x1.
//
// Quadrature: Gauss2 vs Gauss4 error estimate; recursively bisect when the
// L1 difference |I4 - I2| exceeds adaptiveTol.
//
// Note: This path is intended to reduce samples after AMR without depending on
// event-based segmentation.
bool buildOrbitAveragedCacheAdaptiveGauss24(const BackgroundMesh& mesh,
                                           const std::vector<double>& hp,
                                           const CouplingWeightsCache::BuildOptions& buildOpt,
                                           std::size_t numParticles,
                                           const std::vector<Vector2D>& x0,
                                           const std::vector<Vector2D>& x1,
                                           const std::vector<int>& hostCellLabelByParticle,
                                           double dt,
                                           CouplingWeightsCache& out,
                                           OrbitAveragedWeightsBuilderWorkspace* workspace,
                                           const OrbitAveragedWeightsBuilderOptions& opt = {},
                                           OrbitAveragedWeightsBuilderTiming* timing = nullptr);

bool buildOrbitAveragedCacheAdaptiveGauss24(const BackgroundMesh& mesh,
                                           const std::vector<double>& hp,
                                           const CouplingWeightsCache::BuildOptions& buildOpt,
                                           std::size_t numParticles,
                                           const std::vector<Vector2D>& x0,
                                           const std::vector<Vector2D>& x1,
                                           const std::vector<int>& hostCellLabelByParticle,
                                           double dt,
                                           CouplingWeightsCache& out,
                                           const OrbitAveragedWeightsBuilderOptions& opt = {},
                                           OrbitAveragedWeightsBuilderTiming* timing = nullptr);

} // namespace Orbit
} // namespace ImplicitPIC
