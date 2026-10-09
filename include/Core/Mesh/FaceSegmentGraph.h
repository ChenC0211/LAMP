#pragma once

#include <vector>
#include <string>
#include <functional>
#include <utility>

#include "Core/Base/Vector2D.h"
#include "Core/Base/Enumeration.h"

// Forward declarations
class QuadtreeMesh;
class Cell;

// Face order follows Cell::neiCell indexing convention used in the mesh:
//   0: west, 1: south, 2: east, 3: north.
enum class FaceSide : int
{
    West  = 0,
    South = 1,
    East  = 2,
    North = 3
};

inline Boundary faceSideToBoundary(const FaceSide side)
{
    switch (side)
    {
    case FaceSide::West:  return Boundary::west;
    case FaceSide::South: return Boundary::south;
    case FaceSide::East:  return Boundary::east;
    case FaceSide::North: return Boundary::north;
    default:              return Boundary::non;
    }
}

inline Vector2D faceSideToNormal(const FaceSide side)
{
    switch (side)
    {
    case FaceSide::West:  return Vector2D(-1.0, 0.0);
    case FaceSide::South: return Vector2D(0.0, -1.0);
    case FaceSide::East:  return Vector2D(1.0, 0.0);
    case FaceSide::North: return Vector2D(0.0, 1.0);
    default:              return Vector2D(0.0, 0.0);
    }
}

// A single oriented face segment associated with an owner cell K.
//
// - K is a leaf cell label (plasma-active leaf chosen by the caller).
// - L is the neighbor cell label across this segment, or -1 if domain boundary.
// - boundary is Boundary::non for interior segments; otherwise set to west/south/east/north.
// - Af is the segment length (2D); currently the geometric edge length.
// - n is the outward unit normal of cell K (axis-aligned on quadtree).
// - df is the center-to-center spacing along n when L exists; otherwise center-to-face distance.
struct FaceSegment
{
    int K = -1;
    int L = -1;
    Boundary boundary = Boundary::non;

    FaceSide side = FaceSide::West;
    bool isHorizontal = false;
    int edgeLabel = -1;

    double Af = 0.0;
    Vector2D n = Vector2D(0.0, 0.0);
    double df = 0.0;

    Vector2D faceCenter = Vector2D(0.0, 0.0);

    // Embedded-boundary diagnostics (optional; populated if available)
    bool edgeCut = false;
    double edgeLength = 0.0;
    double edgeCutLength = 0.0;
};

// ----------------------------------------------------------------------------
// FaceSegmentGraph
//
// Build oriented face segments for each (selected) leaf cell.
// The graph stores segments in a flat array, and a per-cell offset table for
// quick access.
// ----------------------------------------------------------------------------
class FaceSegmentGraph
{
public:
    // The predicate decides whether a leaf cell participates.
    // If predicate is empty (nullptr), all leaf cells are included.
    using CellActiveFn = std::function<bool(const QuadtreeMesh&, const Cell&)>;

    FaceSegmentGraph() = default;
    ~FaceSegmentGraph() = default;

    void clear();

    // Build from mesh with an optional active-cell predicate.
    void build(const QuadtreeMesh& mesh, CellActiveFn isActive = nullptr);

    // Flat segment array.
    const std::vector<FaceSegment>& segments() const { return segments_; }

    // Offsets by *cell index* in QuadtreeMesh (not label). Range is
    // [cellOffsets[i], cellOffsets[i+1]).
    const std::vector<size_t>& cellOffsets() const { return cellOffsets_; }

    // Convenience: returns [begin,end) indices in segments() for a cell index.
    std::pair<size_t, size_t> segmentRangeForCellIndex(int cellIndex) const;

    // Debug dump: one line per segment.
    // File format: CSV with header.
    void dumpDebugCSV(const QuadtreeMesh& mesh, const std::string& filePath,
                      bool onlyBoundary = false) const;

private:
    std::vector<FaceSegment> segments_;
    std::vector<size_t> cellOffsets_; // size = mesh.cellIndexNum()+1

    static bool isLeafCell(const Cell& cell);
};
