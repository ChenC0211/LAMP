#include "Core/Mesh/NativeLeafDomain.h"

#include "Core/Mesh/CartesianLeafMesh.h"
#include "Core/Mesh/CartesianMeasure.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/NativeMeshIdentity.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/UniformLeafGrid.h"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace
{
void require(bool condition, const char* message)
{
    if (!condition) throw std::invalid_argument(std::string("[NativeLeafDomain] ") + message);
}

bool sameSegment(const FaceSegment& a, const FaceSegment& b)
{
    return a.K == b.K && a.L == b.L && a.boundary == b.boundary && a.side == b.side &&
        a.isHorizontal == b.isHorizontal && a.edgeLabel == b.edgeLabel &&
        a.Af == b.Af && a.df == b.df && a.n.x == b.n.x && a.n.y == b.n.y &&
        a.faceCenter.x == b.faceCenter.x && a.faceCenter.y == b.faceCenter.y;
}
}

NativeLeafDomain NativeLeafDomain::build(const QuadtreeMesh& mesh, const FaceSegmentGraph& graph,
    bool axisymmetric)
{
    NativeLeafDomain result;
    result.identity_ = std::make_shared<NativeMeshIdentity>(mesh);
    for (int i = 1; i < mesh.cellIndexNum(); ++i)
        result.uniformLeaves_ &= mesh.cells(mesh.cellIndexToLabel[i]).level()
            == mesh.cells(mesh.cellIndexToLabel[0]).level();
    if (result.uniformLeaves_)
    {
        const UniformLeafGrid lattice(mesh, CartesianGridScope::UniformLeaves);
        result.spacing_ = lattice.spacing();
    }
    else result.spacing_ = validateCartesianLeafMesh(mesh).minimumSpacing;
    if (axisymmetric) require(mesh.minCoord.x >= 0.0, "RZ requires a nonnegative radial domain");
    for (int i = 0; i < mesh.cellIndexNum(); ++i)
    {
        const auto& cell = mesh.cells(mesh.cellIndexToLabel[i]);
        const double volume = CartesianMeasure::cell(cell.area(), cell.centerCoord().x, axisymmetric);
        require(std::isfinite(volume) && volume > 0.0, "unrepresentable cell volume");
    }
    if (axisymmetric && mesh.minCoord.x == 0.0)
        for (int j = 0; j < mesh.vEdgeIndexNum(); ++j)
        {
            const auto& edge = mesh.vEdges(mesh.vEdgeIndexToLabel[j]);
            if (edge.slaveCell >= 0) continue;
            require(edge.midpointCoord().x == 0.0 && mesh.verts(edge.edgeToVert[0]).coord().x == 0.0 &&
                    mesh.verts(edge.edgeToVert[1]).coord().x == 0.0,
                    "the radial axis must retain its exact zero coordinate");
        }
    result.axisymmetric_ = axisymmetric;
    const int cells = mesh.cellIndexNum(), nv = mesh.vEdgeIndexNum(), nh = mesh.hEdgeIndexNum();
    require(nv <= std::numeric_limits<int>::max() - nh, "native face count overflows indices");
    result.verticalFaces_ = nv; result.horizontalFaces_ = nh;
    const auto& segments = graph.segments();
    const auto& offsets = graph.cellOffsets();
    require(offsets.size() == static_cast<std::size_t>(cells) + 1 && offsets.front() == 0 &&
            offsets.back() == segments.size(), "graph offsets/count mismatch");
    result.cellRowByIndex_.assign(cells, -1);
    for (int ci = 0; ci < cells; ++ci)
    {
        require(offsets[ci] <= offsets[ci+1] && offsets[ci+1] <= segments.size(), "invalid graph cell range");
        if (offsets[ci] == offsets[ci+1]) continue;
        require(offsets[ci+1] - offsets[ci] >= 4 && offsets[ci+1] - offsets[ci] <= 8,
            "active 2:1 leaf must supply all real side segments");
        result.cellRowByIndex_[ci] = static_cast<int>(result.cellIndexByRow_.size());
        result.cellIndexByRow_.push_back(ci);
        const auto& cell = mesh.cells(mesh.cellIndexToLabel[ci]);
        result.cellVolumes_.push_back(CartesianMeasure::cell(cell.area(), cell.centerCoord().x, axisymmetric));
    }
    require(!result.cellIndexByRow_.empty(), "empty selected native leaf domain");
    FaceSegmentGraph actual;
    actual.build(mesh, [&](const QuadtreeMesh&, const Cell& cell) { return result.cellRowByIndex_[cell.index()] >= 0; });
    require(actual.cellOffsets() == offsets && actual.segments().size() == segments.size(),
        "graph is not the complete graph of the selected actual leaf cells");

    // Preserve the original exact oriented-segment and incidence checks. Only
    // physical boundary classification has been removed from this geometry pass.
    std::vector<int> incidence(result.nativeFaceCount(), 0);
    for (std::size_t sid = 0; sid < segments.size(); ++sid)
    {
        const auto& segment = segments[sid];
        require(sameSegment(segment, actual.segments()[sid]), "stale/altered oriented graph identity or metric");
        const auto& edge = segment.isHorizontal ? mesh.hEdges(segment.edgeLabel) : mesh.vEdges(segment.edgeLabel);
        const int row = edge.index() + (segment.isHorizontal ? nv : 0);
        const bool interior = result.activeCellIndex(mesh, edge.slaveCell) >= 0
            && result.activeCellIndex(mesh, edge.masterCell) >= 0;
        require(++incidence[row] <= (interior ? 2 : 1), "duplicate oriented physical face");
    }
    for (int row = 0; row < result.nativeFaceCount(); ++row)
    {
        const auto& edge = row < nv ? mesh.vEdges(mesh.vEdgeIndexToLabel[row]) : mesh.hEdges(mesh.hEdgeIndexToLabel[row-nv]);
        const bool negative = result.activeCellIndex(mesh, edge.slaveCell) >= 0;
        const bool positive = result.activeCellIndex(mesh, edge.masterCell) >= 0;
        require(incidence[row] == static_cast<int>(negative) + static_cast<int>(positive),
            "missing oriented active physical face");
    }
    return result;
}

int NativeLeafDomain::activeCellIndex(const QuadtreeMesh& mesh, int label) const
{
    if (label < 0) return -1;
    const int index = mesh.cells(label).index();
    return cellRowByIndex_[index] >= 0 ? index : -1;
}

NativeLeafFace NativeLeafDomain::face(const QuadtreeMesh& mesh, int nativeRow) const
{
    require(identity_ && mesh.cellIndexNum() == static_cast<int>(cellRowByIndex_.size())
        && mesh.vEdgeIndexNum() == verticalFaces_ && mesh.hEdgeIndexNum() == horizontalFaces_,
        "face query requires a built domain and matching native dimensions");
    require(nativeRow >= 0 && nativeRow < nativeFaceCount(), "native face index out of range");
    NativeLeafFace face;
    face.nativeRow = nativeRow;
    face.horizontal = nativeRow >= verticalFaces_;
    face.edgeLabel = face.horizontal ? mesh.hEdgeIndexToLabel[nativeRow-verticalFaces_] : mesh.vEdgeIndexToLabel[nativeRow];
    const auto& edge = face.horizontal ? mesh.hEdges(face.edgeLabel) : mesh.vEdges(face.edgeLabel);
    face.center = edge.midpointCoord();
    const Cell* fine = nullptr;
    for (int label : {edge.slaveCell, edge.masterCell}) if (label >= 0)
        if (!fine || mesh.cells(label).level() > fine->level()) fine = &mesh.cells(label);
    require(fine != nullptr, "native face has no owner");
    face.normalSpacing = uniformLeaves_ ? (face.horizontal ? spacing_.y : spacing_.x)
        : (face.horizontal ? fine->size().y : fine->size().x);
    face.tangentSpacing = uniformLeaves_ ? (face.horizontal ? spacing_.x : spacing_.y) : edge.length();
    face.area = CartesianMeasure::face(edge.length(), face.center.x, axisymmetric_);
    face.negativeCell = activeCellIndex(mesh, edge.slaveCell);
    face.positiveCell = activeCellIndex(mesh, edge.masterCell);
    if ((face.negativeCell >= 0) != (face.positiveCell >= 0))
        face.outwardSign = face.negativeCell >= 0 ? 1.0 : -1.0;
    return face;
}

void NativeLeafDomain::validateMeshIdentity(const QuadtreeMesh& mesh) const
{
    require(static_cast<bool>(identity_), "domain has not been built");
    identity_->validate(mesh);
}

Vector2D NativeLeafDomain::spacing() const
{
    require(uniformLeaves_, "mixed-level domain has no single global stencil spacing");
    return spacing_;
}
