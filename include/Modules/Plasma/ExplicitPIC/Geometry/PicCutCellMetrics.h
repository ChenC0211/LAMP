#pragma once

#include <string>
#include <vector>
#include <functional>

#include "Core/Base/Vector2D.h"

class QuadtreeMesh;
class Cell;
class FaceSegmentGraph;

template<typename Mesh>
class ImmersedGeometry;

namespace ExplicitPIC
{
namespace Geometry
{

// ============================================================================
// PicCutCellMetrics
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

class PicCutCellMetrics
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
    };

    // Workaround for GCC: avoid BuildOptions{} as a default argument inside the enclosing class.
    static BuildOptions defaultBuildOptions() { return BuildOptions{}; }

    PicCutCellMetrics() = default;
    ~PicCutCellMetrics() = default;

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
};

} // namespace Geometry
} // namespace ExplicitPIC
