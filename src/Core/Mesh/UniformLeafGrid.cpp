#include "Core/Mesh/UniformLeafGrid.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include "Core/Mesh/QuadtreeMesh.h"

namespace
{
void require(bool condition, const char* message)
{
    if (!condition) throw std::invalid_argument(std::string("[UniformLeafGrid] ") + message);
}

Boundary boundaryAt(bool west, bool south, bool east, bool north)
{
    if (south) return west ? Boundary::southWest : east ? Boundary::southEast : Boundary::south;
    if (north) return west ? Boundary::northWest : east ? Boundary::northEast : Boundary::north;
    return west ? Boundary::west : east ? Boundary::east : Boundary::non;
}

Boundary cellBoundaryAt(bool west, bool south, bool east, bool north)
{
    // Match native metadata also when one cell spans both opposite walls.
    Boundary result = west ? Boundary::west : Boundary::non;
    if (east) result = Boundary::east;
    if (south) result = result == Boundary::west ? Boundary::southWest
        : result == Boundary::east ? Boundary::southEast : Boundary::south;
    if (north) result = result == Boundary::west ? Boundary::northWest
        : result == Boundary::east ? Boundary::northEast : Boundary::north;
    return result;
}
} // namespace

UniformLeafGrid::UniformLeafGrid(const QuadtreeMesh& mesh, CartesianGridScope scope)
{
    require(scope == CartesianGridScope::RootOnly || scope == CartesianGridScope::UniformLeaves,
            "invalid Cartesian grid scope");
    require(mesh.cellIndexNum() > 0 && mesh.vertIndexNum() > 0 &&
            mesh.vEdgeIndexNum() > 0 && mesh.hEdgeIndexNum() > 0 &&
            mesh.cellIndexToLabel && mesh.vertIndexToLabel &&
            mesh.vEdgeIndexToLabel && mesh.hEdgeIndexToLabel &&
            mesh.cellLabelCapacity() >= mesh.cellIndexNum() &&
            mesh.vertLabelCapacity() >= mesh.vertIndexNum() &&
            mesh.vEdgeLabelCapacity() >= mesh.vEdgeIndexNum() &&
            mesh.hEdgeLabelCapacity() >= mesh.hEdgeIndexNum(), "missing or undersized indexed mesh storage");
    const bool rootOnly = scope == CartesianGridScope::RootOnly;
    // Preserve the Q4 root-grid tolerance; the extracted fitted-leaf contract
    // retains its existing, explicitly different metadata tolerance.
    const double toleranceFactor = rootOnly ? 64.0 : 256.0;
    const auto near = [toleranceFactor](double a, double b, double scale) {
        return std::isfinite(a) && std::isfinite(b) && std::isfinite(scale) &&
            std::abs(a - b) <= toleranceFactor * std::numeric_limits<double>::epsilon() * scale;
    };
    const int firstLabel = mesh.cellIndexToLabel[0];
    require(firstLabel >= 0 && firstLabel < mesh.cellLabelCapacity(), "invalid first indexed cell label");
    const auto lower = mesh.minCoord, upper = mesh.maxCoord;
    require(std::isfinite(lower.x) && std::isfinite(lower.y) &&
            std::isfinite(upper.x) && std::isfinite(upper.y) &&
            upper.x > lower.x && upper.y > lower.y, "invalid or empty domain");
    lower_ = lower; upper_ = upper;
    const auto& first = mesh.cells(firstLabel);
    require(first.level() >= 0, "invalid leaf level");
    level_ = first.level();
    spacing_ = rootOnly ? mesh.cellSize : first.size();
    const auto h = spacing_;
    const double sx = std::max({std::abs(lower.x), std::abs(upper.x), h.x});
    const double sy = std::max({std::abs(lower.y), std::abs(upper.y), h.y});
    const double resolutionFactor = rootOnly ? 256.0 : 1024.0;
    require(std::isfinite(h.x) && std::isfinite(h.y) && h.x > resolutionFactor * std::numeric_limits<double>::epsilon() * sx &&
            h.y > resolutionFactor * std::numeric_limits<double>::epsilon() * sy,
            "leaf spacing must be positive and numerically resolved");
    const auto dimension = [&](double length, double spacing, double scale) {
        const double count = std::round(length / spacing);
        require(std::isfinite(count) && count >= 1.0 && count <= mesh.cellIndexNum() &&
                near(count * spacing, length, scale), "domain is not a complete uniform leaf lattice");
        return static_cast<int>(count);
    };
    nx_ = dimension(upper.x - lower.x, h.x, sx);
    ny_ = dimension(upper.y - lower.y, h.y, sy);
    if (rootOnly)
        require(nx_ == mesh.Nx && ny_ == mesh.Ny && level_ == 0,
                "RootOnly requires the original level-zero root dimensions and spacing");
    const long long nx = nx_, ny = ny_;
    require(nx * ny == mesh.cellIndexNum() && (nx + 1) * (ny + 1) == mesh.vertIndexNum() &&
            (nx + 1) * ny == mesh.vEdgeIndexNum() && nx * (ny + 1) == mesh.hEdgeIndexNum(),
            "requires complete same-level leaf cells/vertices/native faces, not adaptive or omitted geometry");
    cells_.assign(mesh.cellIndexNum(), -1);
    vertices_.assign(mesh.vertIndexNum(), -1);
    vertical_.assign(mesh.vEdgeIndexNum(), -1);
    horizontal_.assign(mesh.hEdgeIndexNum(), -1);
    const auto coordinateIndex = [&near](double x, double lo, double spacing, double offset, int limit, double scale) {
        const double index = std::round((x - lo) / spacing - offset);
        require(std::isfinite(index) && index >= 0.0 && index <= limit &&
                near(x, lo + (index + offset) * spacing, scale), "geometry is off the actual leaf lattice");
        return static_cast<int>(index);
    };
    const auto put = [](std::vector<int>& labels, std::size_t slot, int label) {
        require(slot < labels.size() && labels[slot] == -1, "duplicate leaf lattice location");
        labels[slot] = label;
    };
    for (int i = 0; i < mesh.cellIndexNum(); ++i)
    {
        const int label = mesh.cellIndexToLabel[i];
        require(label >= 0 && label < mesh.cellLabelCapacity(), "invalid indexed cell label");
        const auto& cell = mesh.cells(label);
        require(cell.label() == label && cell.index() == i && cell.level() == first.level() &&
                near(cell.size().x, h.x, h.x) && near(cell.size().y, h.y, h.y),
                "invalid cell identity or nonuniform leaf size/level");
        for (int son : cell.son) require(son == -1, "indexed cells must be actual leaves");
        require(cell.parent >= -1 && cell.parent < mesh.cellLabelCapacity(),
                "invalid leaf parent label");
        if (rootOnly) require(cell.parent == -1, "RootOnly forbids cell refinement ancestry");
        const auto center = cell.centerCoord();
        const int x = coordinateIndex(center.x, lower.x, h.x, 0.5, nx_ - 1, sx);
        const int y = coordinateIndex(center.y, lower.y, h.y, 0.5, ny_ - 1, sy);
        put(cells_, static_cast<std::size_t>(y) * nx_ + x, label);
        require(cell.boundary() == cellBoundaryAt(x == 0, y == 0, x == nx_ - 1, y == ny_ - 1),
                "cell boundary metadata disagrees with the leaf lattice");
    }
    for (int i = 0; i < mesh.vertIndexNum(); ++i)
    {
        const int label = mesh.vertIndexToLabel[i];
        require(label >= 0 && label < mesh.vertLabelCapacity(), "invalid indexed vertex label");
        const auto& vertex = mesh.verts(label);
        require(vertex.index() == i && vertex.label() == label && vertex.supCell == -1 &&
                vertex.supVert[0] == -1 && vertex.supVert[1] == -1 && vertex.type() == TVertType::non,
                "invalid vertex identity or hanging-node topology");
        if (rootOnly) require(vertex.parent == -1 && vertex.son[0] == -1,
                "RootOnly forbids vertex refinement ancestry");
        const auto p = vertex.coord();
        const int x = coordinateIndex(p.x, lower.x, h.x, 0.0, nx_, sx);
        const int y = coordinateIndex(p.y, lower.y, h.y, 0.0, ny_, sy);
        put(vertices_, static_cast<std::size_t>(y) * (nx_ + 1) + x, label);
        require(vertex.boundary() == boundaryAt(x == 0, y == 0, x == nx_, y == ny_),
                "vertex boundary metadata disagrees with the leaf lattice");
    }
    const auto cellAt = [&](int x, int y) {
        return x < 0 || y < 0 || x >= nx_ || y >= ny_ ? -1 : cells_[static_cast<std::size_t>(y) * nx_ + x];
    };
    const auto vertexAt = [&](int x, int y) { return vertices_[static_cast<std::size_t>(y) * (nx_ + 1) + x]; };
    for (bool horizontal : {false, true})
    {
        const int count = horizontal ? mesh.hEdgeIndexNum() : mesh.vEdgeIndexNum();
        const auto labels = horizontal ? mesh.hEdgeIndexToLabel : mesh.vEdgeIndexToLabel;
        for (int i = 0; i < count; ++i)
        {
            const int label = labels[i];
            require(label >= 0 && label < (horizontal ? mesh.hEdgeLabelCapacity() : mesh.vEdgeLabelCapacity()),
                    "invalid indexed edge label");
            const auto& edge = horizontal ? mesh.hEdges(label) : mesh.vEdges(label);
            require(edge.index() == i && edge.label() == label && edge.isHorizontal() == horizontal &&
                    edge.son[0] == -1 && edge.son[1] == -1 &&
                    near(edge.length(), horizontal ? h.x : h.y, horizontal ? h.x : h.y),
                    "invalid native edge identity/length or split face");
            require(edge.parent >= -1 &&
                    edge.parent < (horizontal ? mesh.hEdgeLabelCapacity() : mesh.vEdgeLabelCapacity()),
                    "invalid edge parent label");
            if (rootOnly) require(edge.parent == -1, "RootOnly forbids edge refinement ancestry");
            const auto p = edge.midpointCoord();
            const int x = coordinateIndex(p.x, lower.x, h.x, horizontal ? 0.5 : 0.0, nx_ - (horizontal ? 1 : 0), sx);
            const int y = coordinateIndex(p.y, lower.y, h.y, horizontal ? 0.0 : 0.5, ny_ - (horizontal ? 0 : 1), sy);
            const int negative = horizontal ? cellAt(x, y - 1) : cellAt(x - 1, y);
            const int positive = cellAt(x, y);
            require(edge.slaveCell == negative && edge.masterCell == positive &&
                    edge.edgeToVert[0] == vertexAt(x, y) &&
                    edge.edgeToVert[1] == vertexAt(x + (horizontal ? 1 : 0), y + (horizontal ? 0 : 1)),
                    "native edge owner/endpoint incidence disagrees with actual leaf lattice");
            const Boundary boundary = horizontal
                ? (y == 0 ? Boundary::south : y == ny_ ? Boundary::north : Boundary::non)
                : (x == 0 ? Boundary::west : x == nx_ ? Boundary::east : Boundary::non);
            require(edge.boundary() == boundary, "edge boundary metadata disagrees with the leaf lattice");
            auto& target = horizontal ? horizontal_ : vertical_;
            put(target, static_cast<std::size_t>(y) * (horizontal ? nx_ : nx_ + 1) + x, label);
        }
    }
    for (int y = 0; y < ny_; ++y)
        for (int x = 0; x < nx_; ++x)
        {
            const auto& cell = mesh.cells(cellAt(x, y));
            const std::array<int, 4> vertices{{vertexAt(x,y), vertexAt(x+1,y), vertexAt(x+1,y+1), vertexAt(x,y+1)}};
            const std::array<int, 4> neighbors{{cellAt(x-1,y),cellAt(x,y-1),cellAt(x+1,y),cellAt(x,y+1)}};
            const int west = vertical_[static_cast<std::size_t>(y)*(nx_+1)+x];
            const int east = vertical_[static_cast<std::size_t>(y)*(nx_+1)+x+1];
            const int south = horizontal_[static_cast<std::size_t>(y)*nx_+x];
            const int north = horizontal_[static_cast<std::size_t>(y+1)*nx_+x];
            require(cell.cellToVEdge[0][0] == west && cell.cellToVEdge[1][0] == east &&
                    cell.cellToHEdge[0][0] == south && cell.cellToHEdge[1][0] == north,
                    "cell/native edge incidence disagrees with leaf lattice");
            for (int side = 0; side < 4; ++side)
                require(cell.cellToVert[side] == vertices[side] && cell.neiCell[side][0] == neighbors[side] &&
                        cell.neiCell[side][1] == -1, "cell corner/neighbor incidence disagrees with leaf lattice");
            for (int side = 0; side < 2; ++side)
                require(cell.cellToVEdge[side][1] == -1 && cell.cellToHEdge[side][1] == -1,
                        "split cell faces are unsupported in this same-level port");
        }
    for (int y = 0; y <= ny_; ++y)
        for (int x = 0; x <= nx_; ++x)
        {
            const auto& vertex = mesh.verts(vertexAt(x,y));
            const std::array<int, 4> adjacent{{cellAt(x-1,y-1),cellAt(x,y-1),cellAt(x,y),cellAt(x-1,y)}};
            for (int quadrant = 0; quadrant < 4; ++quadrant)
                require(vertex.vertToCell[quadrant] == adjacent[quadrant], "vertex/cell incidence disagrees with leaf lattice");
            const std::array<int, 4> neighbors{{x == 0 ? -1 : vertexAt(x-1,y),
                y == 0 ? -1 : vertexAt(x,y-1), x == nx_ ? -1 : vertexAt(x+1,y),
                y == ny_ ? -1 : vertexAt(x,y+1)}};
            for (int side = 0; side < 4; ++side)
                require(vertex.neiVertices[side] == neighbors[side],
                        "vertex neighbor incidence disagrees with leaf lattice");
            require(vertex.vertToVEdge[0] == (y == 0 ? -1 : vertical_[static_cast<std::size_t>(y-1)*(nx_+1)+x]) &&
                    vertex.vertToVEdge[1] == (y == ny_ ? -1 : vertical_[static_cast<std::size_t>(y)*(nx_+1)+x]) &&
                    vertex.vertToHEdge[0] == (x == 0 ? -1 : horizontal_[static_cast<std::size_t>(y)*nx_+x-1]) &&
                    vertex.vertToHEdge[1] == (x == nx_ ? -1 : horizontal_[static_cast<std::size_t>(y)*nx_+x]),
                    "vertex/native edge incidence disagrees with leaf lattice");
        }
}
