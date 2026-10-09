#include "Core/Mesh/CartesianLeafLocator.h"
#include "Core/Mesh/CartesianLeafMesh.h"
#include "Core/Mesh/QuadtreeMesh.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

CartesianLeafLocator::CartesianLeafLocator(const QuadtreeMesh &mesh)
{
    validateCartesianLeafMesh(mesh);
    if (mesh.cellIndexNum() <= 0 || mesh.cellIndexNum() > std::numeric_limits<int>::max() / 2)
        throw std::invalid_argument("CartesianLeafLocator leaf count exceeds node indexing capacity");
    std::vector<Node> leaves(mesh.cellIndexNum());
    for (int i = 0; i < mesh.cellIndexNum(); ++i)
    {
        const auto &cell = mesh.cells(mesh.cellIndexToLabel[i]);
        const auto sw = mesh.verts(cell.cellToVert[0]).coord(), se = mesh.verts(cell.cellToVert[1]).coord();
        const auto ne = mesh.verts(cell.cellToVert[2]).coord(), nw = mesh.verts(cell.cellToVert[3]).coord();
        if (sw.x != nw.x || se.x != ne.x || sw.y != se.y || nw.y != ne.y || !(sw.x < ne.x) || !(sw.y < ne.y))
            throw std::invalid_argument("CartesianLeafLocator requires actual rectangular leaf corners");
        leaves[i] = {sw, ne, -1, -1, i};
    }
    nodes_.reserve(2 * leaves.size() - 1);
    build(leaves, 0, static_cast<int>(leaves.size()));
}
int CartesianLeafLocator::build(std::vector<Node> &leaves, int begin, int end)
{
    const int id = static_cast<int>(nodes_.size());
    Node node = leaves[begin];
    for (int i = begin + 1; i < end; ++i)
    {
        node.lower.x = std::min(node.lower.x, leaves[i].lower.x);
        node.lower.y = std::min(node.lower.y, leaves[i].lower.y);
        node.upper.x = std::max(node.upper.x, leaves[i].upper.x);
        node.upper.y = std::max(node.upper.y, leaves[i].upper.y);
    }
    nodes_.push_back(node);
    if (end - begin == 1)
        return id;
    const bool x = node.upper.x - node.lower.x >= node.upper.y - node.lower.y;
    const int middle = begin + (end - begin) / 2;
    std::nth_element(leaves.begin() + begin, leaves.begin() + middle, leaves.begin() + end,
                     [&](const Node &a, const Node &b) {
                         const double ac =
                             x ? 0.5 * a.lower.x + 0.5 * a.upper.x : 0.5 * a.lower.y + 0.5 * a.upper.y;
                         const double bc =
                             x ? 0.5 * b.lower.x + 0.5 * b.upper.x : 0.5 * b.lower.y + 0.5 * b.upper.y;
                         return ac == bc ? a.index < b.index : ac < bc;
                     });
    const int left = build(leaves, begin, middle), right = build(leaves, middle, end);
    nodes_[id].index = -1;
    nodes_[id].left = left;
    nodes_[id].right = right;
    return id;
}
void CartesianLeafLocator::visit(int id, Vector2D p, IncidentCells &result) const
{
    const auto &node = nodes_[id];
    if (p.x < node.lower.x || p.x > node.upper.x || p.y < node.lower.y || p.y > node.upper.y)
        return;
    if (node.index >= 0)
    {
        if (result.count == static_cast<int>(result.indices.size()))
            throw std::runtime_error("CartesianLeafLocator found overlapping incident leaves");
        result.indices[result.count++] = node.index;
    }
    else
    {
        visit(node.left, p, result);
        visit(node.right, p, result);
    }
}
CartesianLeafLocator::IncidentCells CartesianLeafLocator::incident(Vector2D point) const
{
    if (!std::isfinite(point.x) || !std::isfinite(point.y))
        throw std::invalid_argument("CartesianLeafLocator requires a finite query");
    IncidentCells result;
    if (!nodes_.empty())
        visit(0, point, result);
    std::sort(result.indices.begin(), result.indices.begin() + result.count);
    return result;
}

void CartesianLeafLocator::visitBox(int id, Vector2D lower, Vector2D upper, std::vector<int> &result) const
{
    const auto &node = nodes_[id];
    if (upper.x < node.lower.x || upper.y < node.lower.y || lower.x > node.upper.x || lower.y > node.upper.y)
        return;
    if (node.index >= 0) result.push_back(node.index);
    else
    {
        visitBox(node.left, lower, upper, result);
        visitBox(node.right, lower, upper, result);
    }
}
std::vector<int> CartesianLeafLocator::intersectingBox(Vector2D lower, Vector2D upper) const
{
    if (!std::isfinite(lower.x) || !std::isfinite(lower.y) || !std::isfinite(upper.x) ||
        !std::isfinite(upper.y) || lower.x > upper.x || lower.y > upper.y)
        throw std::invalid_argument("CartesianLeafLocator requires an ordered finite closed box");
    std::vector<int> result;
    if (!nodes_.empty()) visitBox(0, lower, upper, result);
    std::sort(result.begin(), result.end());
    return result;
}
