#include "Core/Mesh/FaceSegmentGraph.h"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <stdexcept>

#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/Cell.h"
#include "Core/Mesh/Edge.h"

namespace
{
// Recursively collect leaf-edge labels under an edge label.
// On a 2:1 balanced quadtree, depth is at most 1 for face splitting, but we
// keep recursion for robustness.
static void collectLeafEdges(const QuadtreeMesh& mesh, const int edgeLabel,
                             const bool isHorizontal, std::vector<int>& out)
{
    const Edge& e = isHorizontal ? mesh.hEdges(edgeLabel) : mesh.vEdges(edgeLabel);
    if (e.son[0] == -1 && e.son[1] == -1)
    {
        out.push_back(edgeLabel);
        return;
    }
    for (int k = 0; k < 2; ++k)
    {
        if (e.son[k] != -1)
        {
            collectLeafEdges(mesh, e.son[k], isHorizontal, out);
        }
    }
}

static void sortByTangentialCoord(const QuadtreeMesh& mesh, const bool isHorizontal,
                                 std::vector<int>& edgeLabels)
{
    std::sort(edgeLabels.begin(), edgeLabels.end(), [&](int a, int b) {
        const Edge& ea = isHorizontal ? mesh.hEdges(a) : mesh.vEdges(a);
        const Edge& eb = isHorizontal ? mesh.hEdges(b) : mesh.vEdges(b);
        const Vector2D ca = ea.midpointCoord();
        const Vector2D cb = eb.midpointCoord();
        return isHorizontal ? (ca.x < cb.x) : (ca.y < cb.y);
    });
}
}

void FaceSegmentGraph::clear()
{
    segments_.clear();
    cellOffsets_.clear();
}

bool FaceSegmentGraph::isLeafCell(const Cell& cell)
{
    // A leaf cell has no children.
    return (cell.son[0] == -1 && cell.son[1] == -1 && cell.son[2] == -1 && cell.son[3] == -1);
}

std::pair<size_t, size_t> FaceSegmentGraph::segmentRangeForCellIndex(int cellIndex) const
{
    if (cellIndex < 0 || static_cast<size_t>(cellIndex + 1) >= cellOffsets_.size())
    {
        return {0, 0};
    }
    return {cellOffsets_[cellIndex], cellOffsets_[cellIndex + 1]};
}

void FaceSegmentGraph::build(const QuadtreeMesh& mesh, CellActiveFn isActive)
{
    clear();

    const int nCell = mesh.cellIndexNum();
    cellOffsets_.assign(static_cast<size_t>(nCell) + 1u, 0u);

    size_t cursor = 0u;
    cellOffsets_[0] = 0u;
    // Scratch belongs to this build call, not to a cell side. Reuse capacity
    // while preserving recursive split-edge expansion and tangential ordering.
    std::vector<int> segEdgeLabels;
    segEdgeLabels.reserve(2);

    for (int ci = 0; ci < nCell; ++ci)
    {
        cellOffsets_[ci] = cursor;

        const int Klabel = mesh.cellIndexToLabel[ci];
        const Cell& K = mesh.cells(Klabel);

        if (!isLeafCell(K))
        {
            cellOffsets_[ci + 1] = cursor;
            continue;
        }
        if (isActive && !isActive(mesh, K))
        {
            cellOffsets_[ci + 1] = cursor;
            continue;
        }

        // Enumerate 4 sides.
        for (int s = 0; s < 4; ++s)
        {
            const FaceSide side = static_cast<FaceSide>(s);
            const bool isHorizontal = (side == FaceSide::South || side == FaceSide::North);

            // Pull edge labels from cellToHEdge/cellToVEdge.
            int baseEdgeLabels[2] = {-1, -1};
            if (isHorizontal)
            {
                const int row = (side == FaceSide::South) ? 0 : 1;
                baseEdgeLabels[0] = K.cellToHEdge[row][0];
                baseEdgeLabels[1] = K.cellToHEdge[row][1];
            }
            else
            {
                const int col = (side == FaceSide::West) ? 0 : 1;
                baseEdgeLabels[0] = K.cellToVEdge[col][0];
                baseEdgeLabels[1] = K.cellToVEdge[col][1];
            }

            segEdgeLabels.clear();

            // Expand each base edge into leaf segments using Edge::son[].
            for (int j = 0; j < 2; ++j)
            {
                if (baseEdgeLabels[j] == -1) continue;

                // In most cases on a 2:1 balanced mesh, baseEdgeLabels[j] is already a leaf.
                // We still expand for safety.
                collectLeafEdges(mesh, baseEdgeLabels[j], isHorizontal, segEdgeLabels);
            }

            if (segEdgeLabels.empty())
            {
                continue;
            }

            // Ensure a deterministic order along tangential direction.
            sortByTangentialCoord(mesh, isHorizontal, segEdgeLabels);

            // For each segment, create an oriented FaceSegment owned by K.
            for (int eLabel : segEdgeLabels)
            {
                const Edge& e = isHorizontal ? mesh.hEdges(eLabel) : mesh.vEdges(eLabel);

                FaceSegment seg;
                seg.K = Klabel;
                seg.side = side;
                seg.isHorizontal = isHorizontal;
                seg.edgeLabel = eLabel;

                // Outward normal is axis-aligned on this quadtree mesh.
                seg.n = faceSideToNormal(side);

                // Neighbor across the edge depends on the side and the master/slave convention.
                // Mesh conventions (BackgroundMesh::findNeighbourCells):
                //   - Horizontal edge: masterCell is north, slaveCell is south.
                //   - Vertical edge:   masterCell is east,  slaveCell is west.
                if (isHorizontal)
                {
                    seg.L = (side == FaceSide::South) ? e.slaveCell : e.masterCell;
                }
                else
                {
                    seg.L = (side == FaceSide::West) ? e.slaveCell : e.masterCell;
                }

                seg.boundary = (seg.L == -1) ? faceSideToBoundary(side) : Boundary::non;

                // Geometry
                seg.Af = e.length();
                seg.faceCenter = e.midpointCoord();

                // Distance metric df
                if (seg.L != -1)
                {
                    const Cell& Lc = mesh.cells(seg.L);
                    Vector2D d = Lc.centerCoord() - K.centerCoord();
                    seg.df = std::fabs(d.dotProduct(seg.n));
                }
                else
                {
                    Vector2D d = seg.faceCenter - K.centerCoord();
                    seg.df = std::fabs(d.dotProduct(seg.n));
                }

                // EB diagnostics if cut information exists on the edge.
                seg.edgeCut = e.cut();
                seg.edgeLength = e.length();
                seg.edgeCutLength = e.cutLength();

                segments_.push_back(seg);
                ++cursor;
            }
        }

        cellOffsets_[ci + 1] = cursor;
    }

    // Ensure last entry is set.
    cellOffsets_[static_cast<size_t>(nCell)] = cursor;
}

void FaceSegmentGraph::dumpDebugCSV(const QuadtreeMesh& mesh, const std::string& filePath,
                                   bool onlyBoundary) const
{
    std::ofstream ofs(filePath);
    if (!ofs)
    {
        throw std::runtime_error("FaceSegmentGraph::dumpDebugCSV: cannot open file: " + filePath);
    }

    ofs << "K_label,K_index,K_level,side,L_label,boundary,Af,nx,ny,df,edgeLabel,isHorizontal,faceCx,faceCy,edgeCut,edgeLen,edgeCutLen\n";
    ofs << std::setprecision(17);

    for (const auto& seg : segments_)
    {
        if (onlyBoundary && seg.boundary == Boundary::non)
        {
            continue;
        }

        const Cell& K = mesh.cells(seg.K);
        ofs << seg.K << ','
            << K.index() << ','
            << K.level() << ','
            << static_cast<int>(seg.side) << ','
            << seg.L << ','
            << static_cast<int>(seg.boundary) << ','
            << seg.Af << ','
            << seg.n.x << ','
            << seg.n.y << ','
            << seg.df << ','
            << seg.edgeLabel << ','
            << (seg.isHorizontal ? 1 : 0) << ','
            << seg.faceCenter.x << ','
            << seg.faceCenter.y << ','
            << (seg.edgeCut ? 1 : 0) << ','
            << seg.edgeLength << ','
            << seg.edgeCutLength
            << '\n';
    }
}
