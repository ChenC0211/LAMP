#pragma once

#include <cstdint>
#include <vector>

#include "Core/Base/Vector2D.h"

class BackgroundMesh;

template<typename T>
class VertField;

namespace ImplicitPIC
{
namespace Orbit
{

struct OrbitSegment
{
    // Cell label that contains the segment interior.
    // Default behavior uses leaf cells; when Options::segmentCellMaxLevel >= 0,
    // this may be a non-leaf (coarsened) cell label.
    int cellLabel = -1;
    double t0 = 0.0;        // segment start time (code)
    double t1 = 0.0;        // segment end time (code)
    Vector2D x0;            // segment start position (code)
    Vector2D x1;            // segment end position (code)
};

class PredictorOrbitSegmenter
{
public:
    struct Stats
    {
        std::uint64_t fastExitHits = 0u;
        std::uint64_t bisectionFallbackHits = 0u;

        // segmentCellMaxLevel diagnostics:
        // How often the host cell label was clamped to a coarser (ancestor) cell.
        std::uint64_t coarsenHits = 0u;
        // Total number of parent links climbed across all clamps.
        std::uint64_t coarsenParentClimb = 0u;
    };

    struct Options
    {
        // Explicit debug controls; the segmenter never reads process state.
        bool trace = false;
        // Use position lookup for both the initial host and subsequent crossings.
        bool ignoreHostHints = false;
        bool periodicX = false;
        bool periodicY = false;

        // Tolerances for bisection when detecting leaf-cell crossings.
        double timeTolRel = 1e-12;
        double timeTolAbs = 0.0;
        int maxBisectionIter = 24;
        int maxSegmentsPerParticle = 512;

        // Optional performance/robustness control:
        // If segmentCellMaxLevel >= 0, segment predictor orbits using cells
        // coarsened up to this maximum level instead of always segmenting on leaf
        // cell boundaries.
        //
        // Rationale: after AMR, leaf sizes can drop sharply, inflating
        // segments/samples per particle and dominating step cost.
        //
        // Notes:
        // - This is an approximation (fewer split points than leaf-accurate).
        // - Embedded-boundary interception splitting is leaf-based; if enabled,
        //   callers should leave this at -1.
        int segmentCellMaxLevel = -1;

        // Position wrap tolerance (half-open interval [min,max)).
        double wrapEps = 1e-12;

        // Mesh in-cell epsilon (passed to isInsideCell).
        double inCellEps = 0.0;

        // ------------------------------------------------------------
        // Embedded boundary (EB) interception splitting (cut leaves)
        // ------------------------------------------------------------
        // If enabled, and if a cut-cell indicator is available, split (and terminate)
        // predictor segments at the first intersection with the embedded boundary
        // inside the current leaf cell.
        bool enableEmbeddedIntercept = false;

        // Vertex level-set field used to build a planar EB approximation per cut leaf.
        // When null, EB interception is disabled.
        const VertField<double>* embeddedLevelSetVert = nullptr;

        // Level-set sign convention: insideIsNegative=true means solid is psi < 0,
        // plasma/fluid is psi >= 0.
        bool embeddedInsideIsNegative = true;

        // Optional cut-fraction alpha by *cell index*. Used to detect cut leaves.
        // If null, EB interception is disabled.
        const std::vector<double>* alphaByCellIndex = nullptr;

        // Cut-cell detection threshold for alpha-based cut flag.
        double cutAlphaEps = 1e-12;

        // Hardening: when alphaByCellIndex is available, terminate segmentation
        // immediately if the active host cell is fully solid (alpha <= eps).
        // This avoids producing orbit segments entirely inside solid regions,
        // which otherwise create many zero-weight coupling stencils.
        bool terminateOnSolidCell = true;
        double solidAlphaEps = 1e-12;
    };

    // Segment straight-line predictor orbits x(t) between x0[p] and x1[p] into
    // leaf-cell-constant segments.
    //
    // Outputs are in CSR-like form:
    //   segmentOffsets[p]..segmentOffsets[p+1]-1 are segments for particle p.
    static void segmentLinearOrbits(const BackgroundMesh& mesh,
                                    const std::vector<Vector2D>& x0,
                                    const std::vector<Vector2D>& x1,
                                    const std::vector<int>* hostCellHints,
                                    double dt,
                                    std::vector<int>& segmentOffsets,
                                    std::vector<OrbitSegment>& segments,
                                    const Options& opt = Options{},
                                    Stats* stats = nullptr);
};

} // namespace Orbit
} // namespace ImplicitPIC
