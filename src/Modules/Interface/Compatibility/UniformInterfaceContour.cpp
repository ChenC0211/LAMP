#include "Modules/Interface/Compatibility/UniformInterfaceContour.h"
#include "Core/Geometry/ImmersedGeometry.h"
#include "Core/Geometry/LevelSetCurvature.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Base/Interpolation.h"
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace UniformContourDetail
{
[[noreturn]] inline void fail(const std::string& stage, const std::string& reason)
{
    throw std::runtime_error("[SharpInterface:Hydro " + stage + "] " + reason);
}

template<typename Mesh>
bool validVertex(const Mesh& mesh, int label)
{
    if (label < 0 || label >= mesh.vertLabelCapacity()) return false;
    const int index = mesh.verts(label).index();
    return index >= 0 && index < mesh.vertIndexNum() && mesh.vertIndexToLabel[index] == label;
}

template<typename Mesh>
bool validCell(const Mesh& mesh, int label)
{
    if (label < 0 || label >= mesh.cellLabelCapacity()) return false;
    const int index = mesh.cells(label).index();
    return index >= 0 && index < mesh.cellIndexNum() && mesh.cellIndexToLabel[index] == label;
}

}

UniformInterfaceContour::UniformInterfaceContour(const QuadtreeMesh& mesh,
    int levelSetIndex, bool axisymmetric, std::vector<InterfaceNode> nodes, Report report)
    : mesh_(&mesh), identity_(mesh), levelSetIndex_(levelSetIndex), axisymmetric_(axisymmetric),
      nodes_(std::move(nodes)), report_(report)
{}

void UniformInterfaceContour::validateMesh(const QuadtreeMesh& mesh) const
{
    if (mesh_ != &mesh || nodes_.empty())
        throw std::invalid_argument("uniform contour belongs to a different mesh or was moved from");
    identity_.validate(mesh);
}

UniformInterfaceContour UniformInterfaceContour::build(QuadtreeMesh& mesh,
    const ImmersedGeometry<QuadtreeMesh>& geometry, int levelSetIndex, bool axisymmetric)
{
    using namespace UniformContourDetail;
    if (levelSetIndex < 0 || levelSetIndex >= static_cast<int>(geometry.levelSets.size()))
        fail("geometry", "selected level-set index is out of range");
    if (mesh.Nx < 4 || mesh.Ny < 4 || mesh.cellIndexNum() != mesh.Nx * mesh.Ny)
        fail("geometry", "only a fixed uniform mesh with at least 4 by 4 cells is supported");
    const double hx = mesh.cellSize.x;
    const double hy = mesh.cellSize.y;
    if (!(hx > 0.0) || !(hy > 0.0) || !std::isfinite(hx) || !std::isfinite(hy))
        fail("geometry", "invalid uniform mesh spacing");
    if (axisymmetric && (!(mesh.minCoord.x >= 0.0) || !std::isfinite(mesh.maxCoord.x)))
        fail("geometry", "RZ requires a nonnegative radial domain");
    const double tolerance = 64.0 * std::numeric_limits<double>::epsilon();
    for (int i = 0; i < mesh.cellIndexNum(); ++i)
    {
        const int label = mesh.cellIndexToLabel[i];
        if (!validCell(mesh, label)) fail("geometry", "invalid cell mapping");
        const auto& cell = mesh.cells(label);
        if (std::abs(cell.size().x - hx) > tolerance * hx
            || std::abs(cell.size().y - hy) > tolerance * hy)
            fail("geometry", "coarse-fine/AMR cells are unsupported by the checked interface path");
        for (int corner = 0; corner < 4; ++corner)
            if (!validVertex(mesh, cell.cellToVert[corner])) fail("geometry", "invalid cell vertex mapping");
        for (int side = 0; side < 2; ++side)
            if (cell.cellToVEdge[side][1] != -1 || cell.cellToHEdge[side][1] != -1)
                fail("geometry", "split cell faces are unsupported by the checked interface path");
    }

    const auto& phi = geometry.levelSets[levelSetIndex];
    const auto& gx = geometry.LS_dx(levelSetIndex);
    const auto& gy = geometry.LS_dy(levelSetIndex);
    const auto& gxx = geometry.LS_dxx(levelSetIndex);
    const auto& gyy = geometry.LS_dyy(levelSetIndex);
    const auto& gxy = geometry.LS_dxy(levelSetIndex);
    auto checkField = [&](const auto& field)
    {
        if (field.size() != mesh.vertIndexNum() || !field.value)
            fail("geometry", "level set and derivative fields must be refreshed on the current mesh");
        for (int i = 0; i < field.size(); ++i)
            if (!std::isfinite(field.value[i])) fail("geometry", "level set or derivative field is non-finite");
    };
    checkField(phi); checkField(gx); checkField(gy);
    checkField(gxx); checkField(gyy); checkField(gxy);

    const double coordinateTolerance = tolerance * std::max({1.0, std::abs(mesh.minCoord.x),
        std::abs(mesh.maxCoord.x), std::abs(mesh.minCoord.y), std::abs(mesh.maxCoord.y)});
    auto inside = [&](const Vector2D& point)
    {
        return std::isfinite(point.x) && std::isfinite(point.y)
            && point.x > mesh.minCoord.x + coordinateTolerance && point.x < mesh.maxCoord.x - coordinateTolerance
            && point.y > mesh.minCoord.y + coordinateTolerance && point.y < mesh.maxCoord.y - coordinateTolerance;
    };
    std::vector<InterfaceNode> nodes;
    nodes.reserve(mesh.vertIndexNum());
    std::unordered_map<int, int> nodeForNegativeVertex;
    auto nodeForVertex = [&](int label)
    {
        const auto existing = nodeForNegativeVertex.find(label);
        if (existing != nodeForNegativeVertex.end()) return existing->second;
        const auto& vertex = mesh.verts(label);
        const int index = vertex.index();
        const double gradient = std::hypot(gx.value[index], gy.value[index]);
        if (!std::isfinite(gradient) || gradient <= 1e-12)
            fail("geometry", "degenerate level-set gradient at a cut-edge vertex");
        // Same negative-face-vertex node construction as the legacy method,
        // with one gradient-scaled projection for a not-yet-perfect signed distance.
        const Vector2D normal(gx.value[index] / gradient, gy.value[index] / gradient);
        const Vector2D point = vertex.coord() - normal * (phi.value[index] / gradient);
        if (!inside(point))
            fail("geometry", "interface touches/leaves the domain or RZ axis; only closed interior interfaces are supported");
        const int nodeIndex = static_cast<int>(nodes.size());
        nodes.emplace_back();
        auto& node = nodes.back();
        node.setIndex(nodeIndex); node.setLabel(nodeIndex);
        node.setRelatedNegVertLabel(label);
        node.setCoord(point); node.setPreCoord(point); node.setPrePreCoord(point);
        node.setNormal(normal); node.setDs(0.0);
        node.setCurvature(0.0); node.setCurvatureInPlane(0.0);
        nodeForNegativeVertex.emplace(label, nodeIndex);
        return nodeIndex;
    };
    auto cutNode = [&](const Edge& edge)
    {
        if (!validVertex(mesh, edge.edgeToVert[0]) || !validVertex(mesh, edge.edgeToVert[1]))
            fail("geometry", "invalid cut-edge vertex mapping");
        const double a = phi.value[mesh.verts(edge.edgeToVert[0]).index()];
        const double b = phi.value[mesh.verts(edge.edgeToVert[1]).index()];
        if (a == 0.0 && b == 0.0)
            fail("geometry", "level set coincides with a mesh edge; this contour topology is unsupported");
        if ((a <= 0.0) == (b <= 0.0)) return -1;
        if (edge.masterCell == -1 || edge.slaveCell == -1)
            fail("geometry", "open/contact-line interface on a domain boundary is unsupported");
        return nodeForVertex(edge.edgeToVert[a <= 0.0 ? 0 : 1]);
    };
    std::vector<int> vNode(mesh.vEdgeIndexNum(), -1), hNode(mesh.hEdgeIndexNum(), -1);
    for (int i = 0; i < mesh.vEdgeIndexNum(); ++i)
    {
        const int label = mesh.vEdgeIndexToLabel[i];
        if (label < 0 || label >= mesh.vEdgeLabelCapacity() || mesh.vEdges(label).index() != i)
            fail("geometry", "invalid vertical-edge mapping");
        vNode[i] = cutNode(mesh.vEdges(label));
    }
    for (int i = 0; i < mesh.hEdgeIndexNum(); ++i)
    {
        const int label = mesh.hEdgeIndexToLabel[i];
        if (label < 0 || label >= mesh.hEdgeLabelCapacity() || mesh.hEdges(label).index() != i)
            fail("geometry", "invalid horizontal-edge mapping");
        hNode[i] = cutNode(mesh.hEdges(label));
    }
    if (nodes.size() < 3) fail("geometry", "no resolved closed interface (at least three nodes are required)");

    std::vector<std::vector<int>> adjacent(nodes.size());
    for (int i = 0; i < mesh.cellIndexNum(); ++i)
    {
        const auto& cell = mesh.cells(mesh.cellIndexToLabel[i]);
        std::vector<int> endpoints;
        auto append = [&](int node)
        {
            if (node >= 0 && std::find(endpoints.begin(), endpoints.end(), node) == endpoints.end())
                endpoints.push_back(node);
        };
        for (int side = 0; side < 2; ++side)
        {
            const int vLabel = cell.cellToVEdge[side][0], hLabel = cell.cellToHEdge[side][0];
            if (vLabel < 0 || vLabel >= mesh.vEdgeLabelCapacity()
                || hLabel < 0 || hLabel >= mesh.hEdgeLabelCapacity())
                fail("geometry", "invalid cell-to-edge mapping");
            const int vi = mesh.vEdges(vLabel).index(), hi = mesh.hEdges(hLabel).index();
            if (vi < 0 || vi >= static_cast<int>(vNode.size()) || hi < 0 || hi >= static_cast<int>(hNode.size()))
                fail("geometry", "invalid cell-to-edge index");
            append(vNode[vi]); append(hNode[hi]);
        }
        if (endpoints.size() > 2)
            fail("geometry", "ambiguous multiple interface segments in one cell; refine the initial uniform mesh");
        if (endpoints.size() == 2)
            for (int side = 0; side < 2; ++side)
            {
                auto& neighbors = adjacent[endpoints[side]];
                const int other = endpoints[1 - side];
                if (std::find(neighbors.begin(), neighbors.end(), other) == neighbors.end()) neighbors.push_back(other);
            }
    }
    for (const auto& neighbors : adjacent)
        if (neighbors.size() != 2)
            fail("geometry", "disconnected/open/branched interface node topology; only one resolved closed contour is supported");
    std::vector<bool> visited(nodes.size(), false);
    int previous = -1, current = 0;
    for (std::size_t visitedCount = 0; visitedCount < nodes.size(); ++visitedCount)
    {
        if (visited[current]) fail("geometry", "multiple disconnected interface contours are unsupported");
        visited[current] = true;
        const int next = adjacent[current][0] == previous ? adjacent[current][1] : adjacent[current][0];
        nodes[current].setUpNodeLabel(next);
        nodes[next].setDownNodeLabel(current);
        previous = current;
        current = next;
    }
    if (current != 0) fail("geometry", "interface contour does not close");

    // Reconstruct all planar derivatives from one local scalar polynomial.
    // Keep the existing negative-vertex projection, contour, ds and force
    // port unchanged. Rebuild the snapshot whenever phi changes. RZ retains
    // its separately validated derivative/azimuthal-curvature path for now.
    std::optional<LevelSetCurvature::Reconstruction<QuadtreeMesh>> planarCurvature;
    if (!axisymmetric)
        planarCurvature.emplace(mesh, phi, 4, LevelSetCurvature::Basis::TensorProduct);

    auto sample = [&](const auto& field, const Cell& cell, const Vector2D& point)
    {
        const auto lower = mesh.verts(cell.cellToVert[0]).coord();
        const double x = std::clamp((point.x - lower.x) / hx, 0.0, 1.0);
        const double y = std::clamp((point.y - lower.y) / hy, 0.0, 1.0);
        const double values[4] = {field.value[mesh.verts(cell.cellToVert[0]).index()],
            field.value[mesh.verts(cell.cellToVert[1]).index()], field.value[mesh.verts(cell.cellToVert[2]).index()],
            field.value[mesh.verts(cell.cellToVert[3]).index()]};
        return (1.0 - y) * ((1.0 - x) * values[0] + x * values[1])
            + y * ((1.0 - x) * values[3] + x * values[2]);
    };
    Report report;
    report.nodeCount = static_cast<int>(nodes.size());
    for (auto& node : nodes)
    {
        const auto point = node.coord();
        const double upLength = (nodes[node.upNodeLabel()].coord() - point).mod();
        const double downLength = (nodes[node.downNodeLabel()].coord() - point).mod();
        if (!std::isfinite(upLength) || !std::isfinite(downLength) || upLength <= tolerance * hx || downLength <= tolerance * hy)
            fail("geometry", "coincident/non-finite neighboring interface nodes");
        node.setDs(0.5 * (upLength + downLength));
        node.relatedNodes[0] = node.label();
        node.relatedNodes[1] = node.upNodeLabel();
        node.relatedNodes[2] = node.downNodeLabel();
        const int cellLabel = mesh.findCellLabel(point);
        if (!validCell(mesh, cellLabel)) fail("geometry", "interface node has no valid interpolation cell");
        const auto& cell = mesh.cells(cellLabel);
        double nx, ny, curvature;
        if (planarCurvature)
        {
            const auto reconstructed = planarCurvature->evaluate(point);
            nx = reconstructed.nx; ny = reconstructed.ny;
            curvature = reconstructed.curvature;
        }
        else
        {
            const double dx = sample(gx, cell, point), dy = sample(gy, cell, point);
            const double gradient = std::hypot(dx, dy);
            if (!std::isfinite(gradient) || gradient <= 1e-12)
                fail("geometry", "degenerate interpolated level-set gradient");
            nx = dx / gradient; ny = dy / gradient;
            curvature = -(ny * ny * sample(gxx, cell, point) - 2.0 * nx * ny * sample(gxy, cell, point)
                + nx * nx * sample(gyy, cell, point)) / gradient;
        }
        const double totalCurvature = curvature - (axisymmetric ? nx / point.x : 0.0);
        if (!std::isfinite(curvature) || !std::isfinite(totalCurvature))
            fail("geometry", "non-finite interface curvature");
        node.setNormal(Vector2D(nx, ny));
        node.setCurvatureInPlane(curvature); node.setCurvature(totalCurvature);
        const double measure = node.ds() * (axisymmetric ? 2.0 * Pi * point.x : 1.0);
        if (!std::isfinite(measure) || measure <= 0.0) fail("geometry", "invalid interface line/surface measure");
        report.surfaceMeasure += measure;
        report.maxAbsCurvature = std::max(report.maxAbsCurvature, std::abs(totalCurvature));

        int entry = 0;
        double weightSum = 0.0;
        for (int row = -2; row <= 2; ++row)
            for (int column = -2; column <= 2; ++column)
            {
                const Vector2D center(cell.centerCoord().x + column * hx, cell.centerCoord().y + row * hy);
                if (!inside(center)) continue;
                const int label = mesh.findCellLabel(center);
                if (!validCell(mesh, label)) fail("geometry", "invalid spreading support cell");
                const auto& support = mesh.cells(label);
                const double weight = Interpolation::diracInterpolation(point, support.centerCoord(), Vector2D(hx, hy));
                if (!std::isfinite(weight) || weight < 0.0) fail("geometry", "invalid Peskin spreading weight");
                if (weight == 0.0) continue;
                node.cellStencil[entry] = support.index();
                node.cellStencilCoeff[entry++] = weight;
                weightSum += weight;
            }
        if (!std::isfinite(weightSum) || weightSum <= 0.0)
            fail("geometry", "interface node has no nonzero spreading support");
        for (int j = 0; j < entry; ++j) node.cellStencilCoeff[j] /= weightSum;
    }
    if (!std::isfinite(report.surfaceMeasure)) fail("geometry", "integrated interface measure overflow");
    return UniformInterfaceContour(mesh, levelSetIndex, axisymmetric, std::move(nodes), report);
}
