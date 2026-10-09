#pragma once

#include <string>
#include <vector>
#include <functional>

#include "Core/Base/Vector2D.h"
#include "Core/Base/Structures.h"

class QuadtreeMesh;
class Cell;
class FaceSegmentGraph;
template<typename T>
class VertField;

template<typename Mesh>
class ImmersedGeometry;

namespace Solvers
{
namespace Poisson
{

// ============================================================================
// CutCellMetrics
// ----------------------------------------------------------------------------
// A lightweight, non-intrusive cache for cut-cell related metrics needed by
// FV/EB Poisson and deposition.
//
// The current implementation provides a conservative baseline build:
//   alpha_K   = 1
//   beta_f    = 1
//   A_gamma,K = 0
//   d_gamma,K = 0
// and two CSV debug dump helpers.
//
// Higher-fidelity alpha/beta/gamma metrics can be built from:
//   - Level-set (psi, gradpsi)
//   - Edge::cutLength() and cut flags
//   - ImmersedGeometry cut-cell utilities (if available)
// ----------------------------------------------------------------------------
// Notes on indexing:
//   - Per-cell metrics are stored by *cell index* (0..mesh.cellIndexNum()-1),
//     consistent with FaceSegmentGraph::cellOffsets().
//   - Per-segment beta is stored aligned with FaceSegmentGraph::segments().
// ============================================================================

class CutCellMetrics
{
public:
    using CellActiveFn = std::function<bool(const QuadtreeMesh&, const Cell&)>;

    struct CellMetric
    {
        double alpha = 1.0;        // volume/area fraction of plasma in cell
        double A_gamma = 0.0;      // EB interface length inside cell (2D)
        double d_gamma = 0.0;      // center-to-EB distance along normal
        Vector2D n_gamma = Vector2D(0.0, 0.0); // EB normal (optional)
        bool isCut = false;        // cut flag for the chosen level-set

        // Optional: geometry boundary mapping for embedded interfaces.
        // gammaBoundaryId is the 1-based Boundary ID from geometry.dict that the
        // cut-cell's local EB segment is closest to.
        // 0 means "unknown / not mapped".
        int gammaBoundaryId = 0;

        // Matched global line index in geometry.dict (0-based), if any.
        // -1 means "unknown".
        int gammaLineIndex = -1;
    };

    struct BuildOptions
    {
        int levelSetI = 0;
        CellActiveFn isActive = nullptr; // if null, infer from FaceSegmentGraph

        // alpha_K approximation control
        // If a cut-cell volume/area fraction utility is not available, we
        // approximate alpha_K by sub-cell sampling of the level-set.
        //   inside = (psi <= 0) when insideIsNegative=true.
        int alphaSampleN = 6;        // sampling resolution per dimension (N x N)
        bool insideIsNegative = true;

        // Optional: map geometry line index (0-based) -> 1-based geometry boundaryId.
        // If provided, CutCellMetrics will attempt to tag each cut-cell's
        // local EB segment with the closest geometry line and propagate its
        // boundaryId into CellMetric::gammaBoundaryId.
        //
        // This enables multi-boundaryId EB BCs by allowing downstream solvers
        // to pick the proper BC based on gammaBoundaryId.
        const std::vector<int>* lineToBoundaryId = nullptr;

        // Optional: override level-set field and gradients (for composite EB).
        // If psiOverride is provided, psiDxOverride/psiDyOverride must be provided too.
        const VertField<double>* psiOverride = nullptr;
        const VertField<double>* psiDxOverride = nullptr;
        const VertField<double>* psiDyOverride = nullptr;

        // Optional: override level-set interface lines (SignInt list).
        // Useful when psiOverride represents a composite of multiple level sets.
        const std::vector<SignInt>* linesOverride = nullptr;

        // Optional: map a discretized geometry-line point from the static
        // geometry frame into the current geometry frame before gamma-tagging.
        // This is used by MovingEB so gammaBoundaryId / gammaLineIndex can be
        // matched against the current advected interface instead of the
        // original static line coordinates.
        std::function<Vector2D(int, const Vector2D&)> transformLinePoint = nullptr;
    };

    // Workaround for GCC: avoid BuildOptions{} as a default argument inside the enclosing class.
    static BuildOptions defaultBuildOptions() { return BuildOptions{}; }

    CutCellMetrics() = default;
    ~CutCellMetrics() = default;

    void clear();

    // Build default metrics for selected cells/segments.
    // This function does NOT modify mesh/geometry.
    void build(
        const QuadtreeMesh& mesh,
        ImmersedGeometry<QuadtreeMesh>& iGeometry,
        const FaceSegmentGraph& faceGraph,
        const BuildOptions& opt = defaultBuildOptions());

    // Accessors
    const std::vector<CellMetric>& cellMetricsByIndex() const { return cell_; }
    const std::vector<double>& betaBySegment() const { return beta_; }
    const std::vector<unsigned char>& activeMaskByCellIndex() const { return active_; }
    unsigned long long revision() const { return revision_; }

    bool isActiveCellIndex(int cellIndex) const;

    // Debug dumps (CSV)
    void dumpCellsCSV(const QuadtreeMesh& mesh, const std::string& filePath,
                      bool onlyActive = true) const;

    void dumpSegmentsCSV(const QuadtreeMesh& mesh, const FaceSegmentGraph& faceGraph,
                         const std::string& filePath,
                         bool onlyBoundary = false) const;

private:
    std::vector<CellMetric> cell_;           // size = mesh.cellIndexNum()
    std::vector<double> beta_;              // size = faceGraph.segments().size()
    std::vector<unsigned char> active_;     // size = mesh.cellIndexNum(), 1 if active
    int levelSetI_ = 0;
    unsigned long long revision_ = 0;
};

} // namespace Poisson
} // namespace Solvers
