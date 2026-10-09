#include "Core/Mesh/NativeMeshIdentity.h"
#include "Core/Mesh/QuadtreeMesh.h"

#include <cmath>
#include <cstring>
#include <stdexcept>
#include <string>

namespace
{
void require(bool condition, const char* message)
{
    if (!condition) throw std::invalid_argument(std::string("[NativeMeshIdentity] ") + message);
}

// Exact ordered words, not a fingerprint. Exclude phase/cut flags and the
// mutable supCellIndex solver cache; retain connectivity used by native stencils.
std::vector<std::uint64_t> identityWords(const QuadtreeMesh& mesh)
{
    require(mesh.cellIndexNum() > 0 && mesh.vertIndexNum() > 0 &&
            mesh.vEdgeIndexNum() > 0 && mesh.hEdgeIndexNum() > 0 &&
            mesh.cellIndexToLabel && mesh.vertIndexToLabel &&
            mesh.vEdgeIndexToLabel && mesh.hEdgeIndexToLabel &&
            mesh.cellLabelCapacity() >= mesh.cellIndexNum() &&
            mesh.vertLabelCapacity() >= mesh.vertIndexNum() &&
            mesh.vEdgeLabelCapacity() >= mesh.vEdgeIndexNum() &&
            mesh.hEdgeLabelCapacity() >= mesh.hEdgeIndexNum(), "missing or undersized indexed mesh storage");
    std::vector<std::uint64_t> words;
    const auto integer = [&](int value) { words.push_back(static_cast<std::uint64_t>(static_cast<std::int64_t>(value))); };
    const auto scalar = [&](double value) {
        require(std::isfinite(value), "non-finite mesh coordinate/metric");
        std::uint64_t bits = 0;
        static_assert(sizeof(bits) == sizeof(value), "binary64 mesh identity requires 64-bit double");
        std::memcpy(&bits, &value, sizeof(bits));
        words.push_back(bits);
    };
    const auto point = [&](Vector2D value) { scalar(value.x); scalar(value.y); };
    integer(mesh.cellIndexNum()); integer(mesh.vertIndexNum());
    integer(mesh.vEdgeIndexNum()); integer(mesh.hEdgeIndexNum());
    integer(mesh.Nx); integer(mesh.Ny);
    point(mesh.minCoord); point(mesh.maxCoord);
    bool mixedLevels=false;
    for (int i = 0; i < mesh.cellIndexNum(); ++i)
    {
        const int label = mesh.cellIndexToLabel[i];
        require(label >= 0 && label < mesh.cellLabelCapacity(), "invalid indexed cell label");
        const auto& cell = mesh.cells(label);
        integer(label); integer(cell.label()); integer(cell.index());
        integer(cell.level()); integer(cell.parent); integer(static_cast<int>(cell.boundary()));
        mixedLevels |= cell.level()!=mesh.cells(mesh.cellIndexToLabel[0]).level();
        point(cell.centerCoord()); point(cell.size());
        for (int son : cell.son) integer(son);
        for (int vertex : cell.cellToVert) integer(vertex);
        for (const auto& side : cell.cellToVEdge) for (int edge : side) integer(edge);
        for (const auto& side : cell.cellToHEdge) for (int edge : side) integer(edge);
        for (const auto& side : cell.neiCell) for (int neighbor : side) integer(neighbor);
    }
    // Adaptive RZ rows follow parent.son[] to native sibling faces. Leaf words
    // alone cannot identify that operator after an ancestor-only mutation.
    if (mixedLevels)
    {
        std::vector<unsigned char> seen(mesh.cellLabelCapacity(),0);
        for (int i=0;i<mesh.cellIndexNum();++i)
        {
            int label=mesh.cells(mesh.cellIndexToLabel[i]).parent;
            while (label>=0)
            {
                require(label<mesh.cellLabelCapacity(), "invalid native ancestor label");
                if (seen[label]) break;
                seen[label]=1;
                const auto& cell=mesh.cells(label);
                integer(label); integer(cell.label()); integer(cell.level()); integer(cell.parent);
                point(cell.centerCoord()); point(cell.size());
                for (int son:cell.son) integer(son);
                for (int vertex:cell.cellToVert) integer(vertex);
                for (const auto& side:cell.cellToVEdge) for (int edge:side) integer(edge);
                for (const auto& side:cell.cellToHEdge) for (int edge:side) integer(edge);
                label=cell.parent;
            }
        }
    }
    for (int i = 0; i < mesh.vertIndexNum(); ++i)
    {
        const int label = mesh.vertIndexToLabel[i];
        require(label >= 0 && label < mesh.vertLabelCapacity(), "invalid indexed vertex label");
        const auto& vertex = mesh.verts(label);
        integer(label); integer(vertex.label()); integer(vertex.index()); point(vertex.coord());
        integer(vertex.parent); integer(vertex.son[0]); integer(static_cast<int>(vertex.boundary()));
        integer(static_cast<int>(vertex.type()));
        integer(vertex.supCell);
        for (int value : vertex.supVert) integer(value);
        for (int value : vertex.vertToCell) integer(value);
        for (int value : vertex.vertToVEdge) integer(value);
        for (int value : vertex.vertToHEdge) integer(value);
        for (int value : vertex.neiVertices) integer(value);
    }
    for (bool horizontal : {false, true})
    {
        const int count = horizontal ? mesh.hEdgeIndexNum() : mesh.vEdgeIndexNum();
        const int capacity = horizontal ? mesh.hEdgeLabelCapacity() : mesh.vEdgeLabelCapacity();
        const auto labels = horizontal ? mesh.hEdgeIndexToLabel : mesh.vEdgeIndexToLabel;
        for (int i = 0; i < count; ++i)
        {
            const int label = labels[i];
            require(label >= 0 && label < capacity, "invalid indexed edge label");
            const auto& edge = horizontal ? mesh.hEdges(label) : mesh.vEdges(label);
            integer(label); integer(edge.label()); integer(edge.index()); integer(edge.parent);
            integer(edge.son[0]); integer(edge.son[1]); integer(edge.masterCell); integer(edge.slaveCell);
            integer(edge.edgeToVert[0]); integer(edge.edgeToVert[1]);
            integer(edge.isHorizontal() ? 1 : 0); integer(static_cast<int>(edge.boundary()));
            point(edge.midpointCoord()); scalar(edge.length());
        }
    }
    return words;
}

} // namespace

NativeMeshIdentity::NativeMeshIdentity(const QuadtreeMesh& mesh) : words_(identityWords(mesh)) {}

void NativeMeshIdentity::validate(const QuadtreeMesh& mesh) const
{
    require(identityWords(mesh) == words_, "mesh geometry/topology changed since identity construction");
}

