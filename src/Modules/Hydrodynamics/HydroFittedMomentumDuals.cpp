#include "Modules/Hydrodynamics/Detail/HydroFittedMomentumDuals.h"

#include "Modules/Hydrodynamics/HydroFittedDomain.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/CartesianMeasure.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace HydroFittedDetail
{
namespace
{
using Face = NativeLeafFace;
using Volume = MomentumDuals::Volume;
using Segment = MomentumDuals::Segment;

void require(bool condition, const char* message)
{
    if (!condition)
        throw std::invalid_argument(std::string("[HydroFittedMomentumDuals] ") + message);
}
double coordinate(Vector2D point, int axis) { return axis == 0 ? point.x : point.y; }
Vector2D pointOnPlane(int normal, double plane, double tangent)
{
    return normal == 0 ? Vector2D(plane, tangent) : Vector2D(tangent, plane);
}
double volumeMeasure(Vector2D low, Vector2D high, bool rz)
{
    return (high.x - low.x) * (high.y - low.y)
        * (rz ? 0.5 * CartesianMeasure::twoPi * (high.x + low.x) : 1.0);
}
double segmentMeasure(const Segment& segment, bool rz)
{
    const int normal = segment.normalAxis;
    const double length = coordinate(segment.end, 1 - normal) - coordinate(segment.begin, 1 - normal);
    return length * (rz ? 0.5 * CartesianMeasure::twoPi * (segment.begin.x + segment.end.x) : 1.0);
}
bool nearMeasure(double a, double b)
{
    return std::isfinite(a) && std::isfinite(b)
        && std::abs(a - b) <= 512.0 * std::numeric_limits<double>::epsilon() * std::max(std::abs(a), std::abs(b));
}
bool active(const Face& face) { return face.negativeCell >= 0 || face.positiveCell >= 0; }
const Cell& cellAt(const QuadtreeMesh& mesh, int index)
{
    return mesh.cells(mesh.cellIndexToLabel[index]);
}
const Edge& edgeAt(const QuadtreeMesh& mesh, const Face& face)
{
    return face.horizontal ? mesh.hEdges(face.edgeLabel) : mesh.vEdges(face.edgeLabel);
}

// At most two native component faces on each side of a Cartesian leaf.
// Each owns its actual half-cell rectangle, split tangentially if necessary.
struct HalfPiece
{
    int native = -1;
    Vector2D lower{0.0, 0.0}, upper{0.0, 0.0};
};
struct HalfSide
{
    std::array<HalfPiece, 2> pieces;
    int count = 0;
};
using CellPieces = std::array<std::array<HalfSide, 2>, 2>; // component, normal side

class Builder
{
public:
    Builder(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
        std::vector<Volume>& volumes, std::vector<Segment>& segments)
        : mesh_(mesh), domain_(domain), volumes_(volumes), segments_(segments),
          pieces_(domain.cellIndexByRow().size())
    {
    }
    void run()
    {
        domain_.validateMeshIdentity(mesh_);
        volumes_.resize(domain_.nativeFaceCount());
        segments_.reserve(5 * domain_.nativeFaceCount());
        makePieces();
        insideCells();
        acrossPrimalFaces();
        validateCoverage();
    }

private:
    const QuadtreeMesh& mesh_;
    const NativeLeafDomain& domain_;
    std::vector<Volume>& volumes_;
    std::vector<Segment>& segments_;
    std::vector<CellPieces> pieces_; // Compact active-cell indexing only

    CellPieces& pieces(int cell) { return pieces_[domain_.cellRowByIndex()[cell]]; }
    const CellPieces& pieces(int cell) const { return pieces_[domain_.cellRowByIndex()[cell]]; }

    void append(int negative, int positive, int normal, double plane, double low, double high,
        int primal, int owner)
    {
        require(low < high && std::isfinite(low) && std::isfinite(high) && std::isfinite(plane),
            "nonpositive or nonfinite common segment");
        require(negative != positive && (negative >= 0 || positive >= 0)
            && (primal >= 0) != (owner >= 0), "invalid segment ownership");
        Segment segment;
        segment.negativeNative = negative;
        segment.positiveNative = positive;
        segment.normalAxis = normal;
        segment.begin = pointOnPlane(normal, plane, low);
        segment.end = pointOnPlane(normal, plane, high);
        segment.primalNative = primal;
        segment.ownerCell = owner;
        segment.measure = segmentMeasure(segment, domain_.axisymmetric());
        require(std::isfinite(segment.measure) && segment.measure >= 0.0
            && (segment.measure > 0.0 || (domain_.axisymmetric() && normal == 0 && plane == 0.0)),
            "invalid common-segment area");
        if (negative >= 0 && positive >= 0)
            require(domain_.face(mesh_,negative).horizontal == domain_.face(mesh_,positive).horizontal,
                "a momentum segment connects different velocity components");
        const int id = static_cast<int>(segments_.size());
        segments_.push_back(segment);
        if (negative >= 0) volumes_[negative].sides[2 * normal + 1].push_back(id);
        if (positive >= 0) volumes_[positive].sides[2 * normal].push_back(id);
    }

    void makePieces()
    {
        for (int index : domain_.cellIndexByRow())
        {
            const auto& cell = cellAt(mesh_, index);
            const auto low = mesh_.verts(cell.cellToVert[0]).coord();
            const auto high = mesh_.verts(cell.cellToVert[2]).coord();
            const auto center = cell.centerCoord();
            auto& local = pieces(index);
            for (int component = 0; component < 2; ++component)
            {
                double total = 0.0;
                for (int side = 0; side < 2; ++side)
                {
                    auto& half = local[component][side];
                    for (int part = 0; part < 2; ++part)
                    {
                        const int label = component == 0 ? cell.cellToVEdge[side][part] : cell.cellToHEdge[side][part];
                        if (label < 0) continue;
                        const auto& edge = component == 0 ? mesh_.vEdges(label) : mesh_.hEdges(label);
                        const int native = edge.index() + (component == 0 ? 0 : mesh_.vEdgeIndexNum());
                        const auto face = domain_.face(mesh_,native);
                        require((side == 0 ? face.positiveCell : face.negativeCell) == index,
                            "half-cell owner does not match its native face");
                        const auto begin = mesh_.verts(edge.edgeToVert[0]).coord();
                        const auto end = mesh_.verts(edge.edgeToVert[1]).coord();
                        auto& piece = half.pieces[half.count++];
                        piece.native = native;
                        piece.lower = component == 0 ? Vector2D(side == 0 ? low.x : center.x, begin.y)
                            : Vector2D(begin.x, side == 0 ? low.y : center.y);
                        piece.upper = component == 0 ? Vector2D(side == 0 ? center.x : high.x, end.y)
                            : Vector2D(end.x, side == 0 ? center.y : high.y);
                        auto& volume = volumes_[native];
                        if (!volume.active)
                        {
                            volume.active = true;
                            volume.lower = piece.lower;
                            volume.upper = piece.upper;
                        }
                        else
                        {
                            require(coordinate(volume.lower, 1 - component) == coordinate(piece.lower, 1 - component)
                                && coordinate(volume.upper, 1 - component) == coordinate(piece.upper, 1 - component),
                                "native face owners disagree on their tangential interval");
                            volume.lower = Vector2D(std::min(volume.lower.x, piece.lower.x), std::min(volume.lower.y, piece.lower.y));
                            volume.upper = Vector2D(std::max(volume.upper.x, piece.upper.x), std::max(volume.upper.y, piece.upper.y));
                        }
                        const double measure = volumeMeasure(piece.lower, piece.upper, domain_.axisymmetric());
                        require(std::isfinite(measure) && measure > 0.0, "invalid half-cell volume");
                        volume.measure += measure;
                        total += measure;
                    }
                    require(half.count > 0, "missing native component face");
                    // Exact endpoint coverage, independent of roundoff in area sums.
                    std::sort(half.pieces.begin(), half.pieces.begin() + half.count,
                        [component](const HalfPiece& a, const HalfPiece& b) {
                            return coordinate(a.lower, 1 - component) < coordinate(b.lower, 1 - component);
                        });
                    double cursor = coordinate(low, 1 - component);
                    for (int part = 0; part < half.count; ++part)
                    {
                        const auto& piece = half.pieces[part];
                        require(coordinate(piece.lower, 1 - component) == cursor, "gap or overlap in half-cell pieces");
                        cursor = coordinate(piece.upper, 1 - component);
                    }
                    require(cursor == coordinate(high, 1 - component), "half-cell pieces do not cover the full side");
                }
                require(nearMeasure(total, volumeMeasure(low, high, domain_.axisymmetric())),
                    "momentum pieces do not partition the active cell volume");
            }
        }
    }

    void insideCells()
    {
        for (int index : domain_.cellIndexByRow())
        {
            const auto& center = cellAt(mesh_, index).centerCoord();
            const auto& local = pieces(index);
            for (int component = 0; component < 2; ++component)
            {
                const auto& negative = local[component][0];
                const auto& positive = local[component][1];
                for (int i = 0; i < negative.count; ++i)
                    for (int j = 0; j < positive.count; ++j)
                    {
                        const auto& a = negative.pieces[i];
                        const auto& b = positive.pieces[j];
                        const double low = std::max(coordinate(a.lower, 1 - component), coordinate(b.lower, 1 - component));
                        const double high = std::min(coordinate(a.upper, 1 - component), coordinate(b.upper, 1 - component));
                        if (low < high)
                            append(a.native, b.native, component, coordinate(center, component), low, high, -1, index);
                    }
                for (int side = 0; side < 2; ++side)
                {
                    const auto& half = local[component][side];
                    if (half.count != 2) continue;
                    const auto& a = half.pieces[0];
                    const auto& b = half.pieces[1];
                    append(a.native, b.native, 1 - component, coordinate(a.upper, 1 - component),
                        coordinate(a.lower, component), coordinate(a.upper, component), -1, index);
                }
            }
        }
    }

    struct Interval { int native; double low, high; };
    struct Intervals { std::array<Interval, 4> values; int count = 0; };
    Intervals touching(const Face& cross, int owner, int side) const
    {
        Intervals result;
        if (owner < 0) return result;
        const int normal = cross.horizontal ? 1 : 0;
        const int component = 1 - normal;
        const auto& edge = edgeAt(mesh_, cross);
        const double start = coordinate(mesh_.verts(edge.edgeToVert[0]).coord(), component);
        const double finish = coordinate(mesh_.verts(edge.edgeToVert[1]).coord(), component);
        const double plane = coordinate(cross.center, normal);
        for (const auto& half : pieces(owner)[component])
            for (int i = 0; i < half.count; ++i)
            {
                const auto& piece = half.pieces[i];
                if (coordinate(side == 0 ? piece.lower : piece.upper, normal) != plane) continue;
                const double low = std::max(start, coordinate(piece.lower, component));
                const double high = std::min(finish, coordinate(piece.upper, component));
                if (low < high) result.values[result.count++] = {piece.native, low, high};
            }
        require(result.count > 0, "no momentum piece touches an active primal cross-face");
        return result;
    }

    void acrossPrimalFaces()
    {
        for (const auto face : domain_.faces(mesh_))
        {
            if (!active(face)) continue;
            const int normal = face.horizontal ? 1 : 0;
            const double plane = coordinate(face.center, normal);
            // Tangential momentum volumes on either side of this primal face.
            const auto negative = touching(face, face.negativeCell, 1);
            const auto positive = touching(face, face.positiveCell, 0);
            if (negative.count > 0 && positive.count > 0)
            {
                for (int i = 0; i < negative.count; ++i)
                    for (int j = 0; j < positive.count; ++j)
                    {
                        const auto& a = negative.values[i];
                        const auto& b = positive.values[j];
                        const double low = std::max(a.low, b.low), high = std::min(a.high, b.high);
                        if (low < high) append(a.native, b.native, normal, plane, low, high, face.nativeRow, -1);
                    }
            }
            else
            {
                const bool outwardPositive = face.negativeCell >= 0;
                const auto& intervals = outwardPositive ? negative : positive;
                for (int i = 0; i < intervals.count; ++i)
                {
                    const auto& part = intervals.values[i];
                    append(outwardPositive ? part.native : -1, outwardPositive ? -1 : part.native,
                        normal, plane, part.low, part.high, face.nativeRow, -1);
                }
                // The normal momentum half-volume has its own exterior side.
                const auto& edge = edgeAt(mesh_, face);
                append(outwardPositive ? face.nativeRow : -1, outwardPositive ? -1 : face.nativeRow,
                    normal, plane, coordinate(mesh_.verts(edge.edgeToVert[0]).coord(), 1 - normal),
                    coordinate(mesh_.verts(edge.edgeToVert[1]).coord(), 1 - normal), face.nativeRow, -1);
            }
        }
    }

    void validateCoverage()
    {
        std::vector<int> uses(segments_.size(), 0);
        for (const auto face : domain_.faces(mesh_))
        {
            auto& volume = volumes_[face.nativeRow];
            require(volume.active == active(face), "momentum volume active mask differs from the fluid domain");
            if (!volume.active) continue;
            require(nearMeasure(volume.measure, volumeMeasure(volume.lower, volume.upper, domain_.axisymmetric())),
                "native momentum volume has a gap or overlap");
            for (int side = 0; side < 4; ++side)
            {
                const int normal = side / 2, tangent = 1 - normal;
                const double plane = coordinate(side % 2 == 0 ? volume.lower : volume.upper, normal);
                auto& ids = volume.sides[side];
                std::sort(ids.begin(), ids.end(), [&](int a, int b) {
                    return coordinate(segments_[a].begin, tangent) < coordinate(segments_[b].begin, tangent);
                });
                double cursor = coordinate(volume.lower, tangent);
                for (int id : ids)
                {
                    const auto& segment = segments_[id];
                    require(segment.normalAxis == normal && coordinate(segment.begin, normal) == plane
                        && coordinate(segment.end, normal) == plane
                        && coordinate(segment.begin, tangent) == cursor,
                        "gap, overlap or misplaced common momentum segment");
                    require((side % 2 == 0 ? segment.positiveNative : segment.negativeNative) == face.nativeRow,
                        "common momentum segment has the wrong outward orientation");
                    cursor = coordinate(segment.end, tangent);
                    ++uses[id];
                }
                require(cursor == coordinate(volume.upper, tangent), "momentum side is not fully covered");
            }
        }
        for (std::size_t id = 0; id < segments_.size(); ++id)
        {
            const auto& segment = segments_[id];
            const int expected = (segment.negativeNative >= 0 ? 1 : 0) + (segment.positiveNative >= 0 ? 1 : 0);
            require(uses[id] == expected, "common momentum segment is missing a reciprocal incidence");
        }
    }
};
}

MomentumDuals MomentumDuals::build(const QuadtreeMesh& mesh, const HydroFittedDomain& domain)
{
    domain.validateMeshIdentity(mesh);
    return build(mesh,domain.geometry());
}
MomentumDuals MomentumDuals::build(const QuadtreeMesh& mesh, const NativeLeafDomain& domain)
{
    MomentumDuals result;
    Builder(mesh, domain, result.volumes_, result.segments_).run();
    return result;
}
}
