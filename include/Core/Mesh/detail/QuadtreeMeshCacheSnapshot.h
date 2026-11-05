#pragma once

// Private staging storage for QMC2. This does not serialize a second format or
// recompute cut/derivative data: accepted records retain their stored values.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "Core/Mesh/MeshCacheIO.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/DyadicGridCoordinates.h"

namespace QuadtreeCacheDetail
{
inline void require(bool condition, const std::string& reason)
{
    if (!condition) throw std::runtime_error(reason);
}

inline std::uint64_t addBytes(std::uint64_t a, std::uint64_t b)
{
    require(b <= (std::numeric_limits<std::uint64_t>::max)() - a, "serialized byte count overflow");
    return a + b;
}
inline std::uint64_t multiplyBytes(std::uint64_t a, std::uint64_t b)
{
    require(a == 0 || b <= (std::numeric_limits<std::uint64_t>::max)() / a, "serialized byte count overflow");
    return a * b;
}

struct Snapshot
{
    // Entity order is vertex, horizontal edge, vertical edge, cell.
    int levelSets = 0, maxLevel = 0;
    std::array<int, 4> count{}, capacity{}, ptr{};
    std::array<std::unique_ptr<int[]>, 4> indexToLabel, labelSpace;
    std::unique_ptr<Vert[]> vertices;
    std::unique_ptr<Edge[]> hEdges, vEdges;
    std::unique_ptr<Cell[]> cells;
    std::vector<int> minNegative, minPositive, maxLocal, cellsPerLevel;

    void read(std::istream& input, MeshCacheIO::ReadBudget& budget)
    {
        char magic[4]{};
        const char expected[4] = {'Q', 'M', 'C', '2'};
        MeshCacheIO::readRawArray(input, magic, 4);
        require(MeshCacheIO::magicEquals(magic, expected), "incompatible mesh cache magic (expected QMC2)");
        std::uint32_t version = 0;
        MeshCacheIO::readRaw(input, version);
        require(version == 2, "unsupported mesh cache version " + std::to_string(version));
        MeshCacheIO::readRaw(input, levelSets);
        require(levelSets >= 0, "level-set count must be nonnegative");
        for (int& value : count)
        {
            MeshCacheIO::readRaw(input, value);
            require(value > 0, "active entity count must be positive");
        }
        std::uint64_t mappingBytes = 0;
        for (int value : count) mappingBytes = addBytes(mappingBytes, multiplyBytes(value, sizeof(int)));
        MeshCacheIO::ReadBudget::requireBytes(input, mappingBytes);
        for (int kind = 0; kind < 4; ++kind)
        {
            indexToLabel[kind] = std::make_unique<int[]>(count[kind]);
            MeshCacheIO::readRawArray(input, indexToLabel[kind].get(), count[kind]);
        }
        MeshCacheIO::readRaw(input, maxLevel);
        require(maxLevel >= 0, "maximum refinement level must be nonnegative");
        for (int& value : capacity) MeshCacheIO::readRaw(input, value);
        for (int& value : ptr) MeshCacheIO::readRaw(input, value);
        for (int kind = 0; kind < 4; ++kind)
            require(capacity[kind] >= count[kind] && ptr[kind] >= 0 && ptr[kind] <= capacity[kind],
                    "invalid entity capacity/high-water mark");

        // Exact minimum serialized bytes (empty history vectors) per record.
        // Capacity, not active count, is serialized by the existing QMC2 writer.
        const auto n = static_cast<std::uint64_t>(levelSets);
        const std::uint64_t vertexBytes = addBytes(28 * sizeof(int) + sizeof(Vector2D) + sizeof(bool) +
            sizeof(Boundary) + sizeof(TVertType) + 8 * sizeof(DiscCoeff),
            multiplyBytes(n, 2 * sizeof(unsigned char) + sizeof(Vector2D) + sizeof(int)));
        const std::uint64_t edgeBytes = addBytes(12 * sizeof(int) + sizeof(Vector2D) + 2 * sizeof(double) +
            3 * sizeof(bool) + sizeof(Boundary), multiplyBytes(n, sizeof(unsigned char)));
        const std::uint64_t cellBytes = addBytes(34 * sizeof(int) + 2 * sizeof(Vector2D) +
            3 * sizeof(bool) + sizeof(Boundary), multiplyBytes(n, sizeof(unsigned char)));
        const std::array<std::uint64_t, 4> minimumBytes{{vertexBytes, edgeBytes, edgeBytes, cellBytes}};
        std::uint64_t requiredBytes = 4 * sizeof(int); // trailing vector lengths
        for (int kind = 0; kind < 4; ++kind)
            requiredBytes = addBytes(requiredBytes,
                multiplyBytes(capacity[kind], addBytes(minimumBytes[kind], sizeof(int))));
        require(requiredBytes <= budget.remaining(), "entity capacities/level-set count exceed remaining cache file bytes");
        for (int kind = 0; kind < 4; ++kind)
        {
            labelSpace[kind] = std::make_unique<int[]>(capacity[kind]);
            MeshCacheIO::readRawArray(input, labelSpace[kind].get(), capacity[kind]);
        }
        vertices = std::make_unique<Vert[]>(capacity[0]);
        hEdges = std::make_unique<Edge[]>(capacity[1]);
        vEdges = std::make_unique<Edge[]>(capacity[2]);
        cells = std::make_unique<Cell[]>(capacity[3]);
        for (int i = 0; i < capacity[0]; ++i) MeshCacheIO::readVert(input, vertices[i], levelSets);
        for (int i = 0; i < capacity[1]; ++i) MeshCacheIO::readEdge(input, hEdges[i], levelSets);
        for (int i = 0; i < capacity[2]; ++i) MeshCacheIO::readEdge(input, vEdges[i], levelSets);
        for (int i = 0; i < capacity[3]; ++i) MeshCacheIO::readCell(input, cells[i], levelSets);
        MeshCacheIO::readVector(input, minNegative);
        MeshCacheIO::readVector(input, minPositive);
        MeshCacheIO::readVector(input, maxLocal);
        MeshCacheIO::readVector(input, cellsPerLevel);
        require(budget.remaining() == 0, "unexpected trailing bytes in QMC2 cache");
        // Preserve the old reader's padding contract, but do it before commit.
        if (minNegative.size() < n) minNegative.resize(static_cast<std::size_t>(n), 0);
        if (minPositive.size() < n) minPositive.resize(static_cast<std::size_t>(n), 0);
        if (maxLocal.size() < n) maxLocal.resize(static_cast<std::size_t>(n), 0);
    }

    template<class Entity>
    void validateBookkeeping(const Entity* objects, int kind, const char* name) const
    {
        int leafIndex = 0, freeIndex = 0, lastLive = 0;
        for (int label = 0; label < capacity[kind]; ++label)
        {
            const auto& object = objects[label];
            require(object.label() == -1 || object.label() == label,
                    std::string(name) + " label does not match its storage slot");
            if (object.label() == -1)
            {
                // Only this prefix of the free-label storage is initialized.
                require(labelSpace[kind][freeIndex++] == label, std::string(name) + " free-label prefix is inconsistent");
                continue;
            }
            lastLive = label + 1;
            if (object.son[0] == -1)
            {
                require(leafIndex < count[kind] && indexToLabel[kind][leafIndex] == label && object.index() == leafIndex,
                        std::string(name) + " active index mapping is inconsistent");
                ++leafIndex;
            }
        }
        require(leafIndex == count[kind] && ptr[kind] >= lastLive,
                std::string(name) + " active count/high-water mark is inconsistent");
    }

    // Called only after the complete topology validation below. Root labels are
    // storage identities, not row-major lattice addresses. Recover each root's
    // already-validated location, then derive every descendant from its parent.
    // No coordinate, cut data or history is repaired when reading an old cache.
    void validateDyadicCoordinates(const QuadtreeMesh& mesh) const
    {
        using LAMP::Mesh::DyadicGridCoordinates;
        using Index = std::uint64_t;
        struct Address
        {
            Index x = 0, y = 0;
            unsigned level = 0;
            bool assigned = false;
        };
        const DyadicGridCoordinates grid(mesh.minCoord, mesh.maxCoord, mesh.Nx, mesh.Ny);
        const std::string incompatible =
            ": cached geometry does not match the canonical dyadic lattice; rebuild the mesh cache from the case";
        auto samePoint = [](const Vector2D& a, const Vector2D& b)
        { return a.x == b.x && a.y == b.y; };
        auto doubled = [](Index value, Index offset)
        {
            require(value <= ((std::numeric_limits<Index>::max)() - offset) / 2,
                    "cached dyadic address overflow");
            return 2 * value + offset;
        };
        auto normalized = [](Address address)
        {
            while (address.level > 0 && (address.x & 1u) == 0 && (address.y & 1u) == 0)
            {
                address.x /= 2;
                address.y /= 2;
                --address.level;
            }
            return address;
        };
        auto atLevel = [](Index value, unsigned oldLevel, unsigned newLevel)
        {
            require(newLevel >= oldLevel && newLevel - oldLevel < 64,
                    "cached dyadic address has an unsupported level span");
            const unsigned shift = newLevel - oldLevel;
            require(value <= (std::numeric_limits<Index>::max)() >> shift,
                    "cached dyadic address overflow");
            return value << shift;
        };

        std::vector<Address> cellAddress(static_cast<std::size_t>(capacity[3]));
        std::vector<Address> vertexAddress(static_cast<std::size_t>(capacity[0]));
        std::vector<int> pending;
        pending.reserve(static_cast<std::size_t>(capacity[3]));
        for (int label = 0; label < capacity[3]; ++label)
        {
            const auto& cell = cells[label];
            if (cell.label() == -1 || cell.parent != -1) continue;
            const double x = (cell.centerCoord().x - mesh.minCoord.x) / mesh.cellSize.x - 0.5;
            const double y = (cell.centerCoord().y - mesh.minCoord.y) / mesh.cellSize.y - 0.5;
            // The preceding root-location checks made these rounded indices
            // nonnegative, in range, unique, and complete over the base grid.
            cellAddress[label] = {static_cast<Index>(std::round(x)), static_cast<Index>(std::round(y)), 0, true};
            pending.push_back(label);
        }
        for (std::size_t next = 0; next < pending.size(); ++next)
        {
            const int label = pending[next];
            const auto& cell = cells[label];
            const auto address = cellAddress[label];
            const std::string context = "cell " + std::to_string(label);
            require(address.level == static_cast<unsigned>(cell.level()), context + ": dyadic level/topology mismatch");
            const auto center = grid.point(doubled(address.x, 1), doubled(address.y, 1), address.level + 1);
            require(samePoint(cell.centerCoord(), center) && samePoint(cell.size(), grid.cellSize(address.level)),
                    context + incompatible);
            for (int corner = 0; corner < 4; ++corner)
            {
                const Index ox = corner == 1 || corner == 2 ? 1 : 0;
                const Index oy = corner >= 2 ? 1 : 0;
                require(address.x < (std::numeric_limits<Index>::max)() &&
                        address.y < (std::numeric_limits<Index>::max)(), "cached dyadic corner address overflow");
                Address point{address.x + ox, address.y + oy, address.level, true};
                const int vertexLabel = cell.cellToVert[corner];
                require(samePoint(vertices[vertexLabel].coord(), grid.point(point.x, point.y, point.level)),
                        context + ": corner vertex " + std::to_string(vertexLabel) + incompatible);
                point = normalized(point);
                auto& saved = vertexAddress[vertexLabel];
                require(!saved.assigned || (saved.x == point.x && saved.y == point.y && saved.level == point.level),
                        context + ": shared corner has inconsistent dyadic addresses");
                saved = point;

                if (cell.son[0] == -1) continue;
                const int child = cell.son[corner];
                require(!cellAddress[child].assigned, context + ": child is reachable from multiple dyadic addresses");
                cellAddress[child] = {doubled(address.x, ox), doubled(address.y, oy), address.level + 1, true};
                pending.push_back(child);
            }
        }
        for (int label = 0; label < capacity[3]; ++label)
            require(cells[label].label() == -1 || cellAddress[label].assigned,
                    "live cell is unreachable from the cached root lattice");
        for (int label = 0; label < capacity[0]; ++label)
            require(vertices[label].label() == -1 || vertexAddress[label].assigned,
                    "live vertex has no cached dyadic cell-corner incidence");

        // Endpoints already have topology-derived addresses. This includes
        // inactive parent edges and both pieces at a hanging-node face. Their
        // midpoint must be the SAME lattice point used by a cell center/vertex,
        // not a fresh floating average of rounded endpoint coordinates.
        for (int kind : {1, 2})
        {
            const bool horizontal = kind == 1;
            const auto* edges = horizontal ? hEdges.get() : vEdges.get();
            for (int label = 0; label < capacity[kind]; ++label)
            {
                const auto& edge = edges[label];
                if (edge.label() == -1) continue;
                const std::string context = std::string(horizontal ? "horizontal edge " : "vertical edge ") +
                    std::to_string(label);
                const auto a = vertexAddress[edge.edgeToVert[0]];
                const auto b = vertexAddress[edge.edgeToVert[1]];
                const unsigned level = std::max(a.level, b.level);
                const Index ax = atLevel(a.x, a.level, level), ay = atLevel(a.y, a.level, level);
                const Index bx = atLevel(b.x, b.level, level), by = atLevel(b.y, b.level, level);
                require(horizontal ? (ay == by && bx > ax) : (ax == bx && by > ay),
                        context + ": endpoints do not form a positive dyadic edge");
                const Index span = horizontal ? bx - ax : by - ay;
                require((span & (span - 1)) == 0, context + ": edge spans a non-dyadic interval");
                unsigned spanLevel = 0;
                for (Index value = span; value > 1; value >>= 1) ++spanLevel;
                require(spanLevel <= level, context + ": edge is larger than a background-grid edge");
                require(ax <= (std::numeric_limits<Index>::max)() - bx &&
                        ay <= (std::numeric_limits<Index>::max)() - by, "cached dyadic midpoint address overflow");
                const auto midpoint = grid.point(ax + bx, ay + by, level + 1);
                const auto size = grid.cellSize(level - spanLevel);
                require(samePoint(edge.midpointCoord(), midpoint) &&
                        edge.length() == (horizontal ? size.x : size.y), context + incompatible);
            }
        }
    }

    void validate(const QuadtreeMesh& mesh) const
    {
        validateBookkeeping(vertices.get(), 0, "vertex");
        validateBookkeeping(hEdges.get(), 1, "horizontal edge");
        validateBookkeeping(vEdges.get(), 2, "vertical edge");
        validateBookkeeping(cells.get(), 3, "cell");
        require(mesh.Nx > 0 && mesh.Ny > 0 && std::isfinite(mesh.cellSize.x) && mesh.cellSize.x > 0.0 &&
                std::isfinite(mesh.cellSize.y) && mesh.cellSize.y > 0.0 &&
                std::isfinite(mesh.minCoord.x) && std::isfinite(mesh.minCoord.y) &&
                std::isfinite(mesh.maxCoord.x) && std::isfinite(mesh.maxCoord.y), "invalid receiving mesh domain");
        const double sx = std::max({std::abs(mesh.minCoord.x), std::abs(mesh.maxCoord.x), mesh.cellSize.x});
        const double sy = std::max({std::abs(mesh.minCoord.y), std::abs(mesh.maxCoord.y), mesh.cellSize.y});
        const double eps = 256.0 * std::numeric_limits<double>::epsilon();
        auto near = [&](double a, double b, double scale)
        { return std::isfinite(a) && std::isfinite(b) && std::abs(a - b) <= eps * scale; };
        auto pointNear = [&](const Vector2D& a, const Vector2D& b)
        { return near(a.x, b.x, sx) && near(a.y, b.y, sy); };
        auto inDomain = [&](const Vector2D& point)
        {
            return std::isfinite(point.x) && std::isfinite(point.y) &&
                point.x >= mesh.minCoord.x - eps * sx && point.x <= mesh.maxCoord.x + eps * sx &&
                point.y >= mesh.minCoord.y - eps * sy && point.y <= mesh.maxCoord.y + eps * sy;
        };
        auto optionalLabel = [&](int label, int kind, const std::string& context)
        { require(label == -1 || (label >= 0 && label < capacity[kind]), context + ": topology label out of range"); };
        auto vertex = [&](int label, const std::string& context) -> const Vert&
        {
            require(label >= 0 && label < capacity[0] && vertices[label].label() == label,
                    context + ": missing/live vertex label out of range");
            return vertices[label];
        };
        for (int label = 0; label < capacity[0]; ++label)
        {
            const auto& value = vertices[label];
            if (value.label() == -1) continue;
            const std::string context = "vertex " + std::to_string(label);
            require(inDomain(value.coord()), context + ": coordinate is non-finite or outside receiving domain");
            for (int item : value.vertToHEdge) optionalLabel(item, 1, context);
            for (int item : value.vertToVEdge) optionalLabel(item, 2, context);
            for (int item : value.vertToCell) optionalLabel(item, 3, context);
            for (int item : value.neiVertices) optionalLabel(item, 0, context);
            for (int item : value.supVert) optionalLabel(item, 0, context);
            optionalLabel(value.supCell, 3, context);
            optionalLabel(value.parent, 0, context);
            optionalLabel(value.son[0], 0, context);
            require(value.parent != label && value.son[0] != label, context + ": self-referencing parent/child");
            require(value.relatedVertsNum >= 0 && value.relatedVertsNum <= 6,
                    context + ": invalid related vertex count");
            for (int i = 0; i < value.relatedVertsNum; ++i) optionalLabel(value.relatedVerts[i], 0, context);
        }

        for (int kind : {1, 2})
        {
            const bool horizontal = kind == 1;
            const auto* edges = horizontal ? hEdges.get() : vEdges.get();
            for (int label = 0; label < capacity[kind]; ++label)
            {
                const auto& edge = edges[label];
                if (edge.label() == -1) continue;
                const std::string context = std::string(horizontal ? "horizontal edge " : "vertical edge ") + std::to_string(label);
                const auto a = vertex(edge.edgeToVert[0], context).coord();
                const auto b = vertex(edge.edgeToVert[1], context).coord();
                const double length = horizontal ? b.x - a.x : b.y - a.y;
                require(edge.isHorizontal() == horizontal && inDomain(edge.midpointCoord()) &&
                        std::isfinite(edge.length()) && edge.length() > 0.0 && length > 0.0 &&
                        (horizontal ? near(a.y, b.y, sy) : near(a.x, b.x, sx)),
                        context + ": invalid axis/endpoints/positive length");
                // Include inactive parent edges: coarsening can reactivate them.
                // Never repair this metric alone while retaining old cut data.
                require(near(edge.length(), length, horizontal ? sx : sy) &&
                        pointNear(edge.midpointCoord(), a + (b - a) * 0.5),
                        context + ": stored length/midpoint disagrees with endpoints (possibly an old rectangular-grid cache)");
                optionalLabel(edge.masterCell, 3, context);
                optionalLabel(edge.slaveCell, 3, context);
                optionalLabel(edge.parent, kind, context);
                for (int child : edge.son) optionalLabel(child, kind, context);
                require((edge.son[0] == -1) == (edge.son[1] == -1), context + ": partial child pair");
                if (edge.parent != -1)
                {
                    const auto& parent = edges[edge.parent];
                    require(edge.parent != label && parent.label() == edge.parent &&
                            (parent.son[0] == label || parent.son[1] == label), context + ": invalid parent backlink");
                    require(near(2.0 * edge.length(), parent.length(), horizontal ? sx : sy),
                            context + ": child/parent length mismatch");
                }
                if (edge.son[0] != -1)
                {
                    require(edge.son[0] != label && edge.son[1] != label && edge.son[0] != edge.son[1],
                            context + ": cyclic/duplicate children");
                    const auto& first = edges[edge.son[0]];
                    const auto& second = edges[edge.son[1]];
                    require(first.label() == edge.son[0] && second.label() == edge.son[1] &&
                            first.parent == label && second.parent == label, context + ": invalid child backlink");
                    require(first.edgeToVert[0] == edge.edgeToVert[0] && second.edgeToVert[1] == edge.edgeToVert[1] &&
                            first.edgeToVert[1] == second.edgeToVert[0] &&
                            near(first.length(), 0.5 * edge.length(), horizontal ? sx : sy) &&
                            near(second.length(), 0.5 * edge.length(), horizontal ? sx : sy),
                            context + ": children do not geometrically bisect their parent");
                }
            }
        }

        require(static_cast<std::uint64_t>(maxLevel) + 1u == cellsPerLevel.size(),
                "per-level cell counts do not match maximum refinement level");
        std::vector<int> observedCounts(cellsPerLevel.size(), 0);
        const auto baseCount = static_cast<std::uint64_t>(mesh.Nx) * mesh.Ny;
        require(baseCount <= static_cast<std::uint64_t>(capacity[3]), "receiving domain grid is larger than cached storage");
        std::vector<unsigned char> roots(static_cast<std::size_t>(baseCount), 0);
        for (int label = 0; label < capacity[3]; ++label)
        {
            const auto& cell = cells[label];
            if (cell.label() == -1) continue;
            const std::string context = "cell " + std::to_string(label);
            require(cell.level() >= 0 && cell.level() <= maxLevel, context + ": refinement level out of range");
            const Vector2D expectedSize(std::ldexp(mesh.cellSize.x, -cell.level()), std::ldexp(mesh.cellSize.y, -cell.level()));
            require(expectedSize.x > 0.0 && expectedSize.y > 0.0 &&
                    near(cell.size().x, expectedSize.x, expectedSize.x) &&
                    near(cell.size().y, expectedSize.y, expectedSize.y) &&
                    std::isfinite(cell.area()) && cell.area() > 0.0 && inDomain(cell.centerCoord()),
                    context + ": size/level/domain mismatch");
            for (int corner = 0; corner < 4; ++corner)
            {
                const Vector2D expected = cell.centerCoord() + Vector2D(
                    (corner == 1 || corner == 2 ? 0.5 : -0.5) * cell.size().x,
                    (corner >= 2 ? 0.5 : -0.5) * cell.size().y);
                require(inDomain(expected) && pointNear(vertex(cell.cellToVert[corner], context).coord(), expected),
                        context + ": corners disagree with center/size/domain");
                for (int part = 0; part < 2; ++part) optionalLabel(cell.neiCell[corner][part], 3, context);
                optionalLabel(cell.son[corner], 3, context);
            }
            for (int side = 0; side < 2; ++side)
                for (int part = 0; part < 2; ++part)
                {
                    optionalLabel(cell.cellToHEdge[side][part], 1, context);
                    optionalLabel(cell.cellToVEdge[side][part], 2, context);
                }
            optionalLabel(cell.parent, 3, context);
            const bool leaf = cell.son[0] == -1;
            for (int child : cell.son) require((child == -1) == leaf, context + ": incomplete child quartet");
            if (leaf) ++observedCounts[static_cast<std::size_t>(cell.level())];
            if (cell.parent == -1)
            {
                require(cell.level() == 0, context + ": refined cell has no parent");
                const double x = (cell.centerCoord().x - mesh.minCoord.x) / mesh.cellSize.x - 0.5;
                const double y = (cell.centerCoord().y - mesh.minCoord.y) / mesh.cellSize.y - 0.5;
                require(x > -0.5 && x < mesh.Nx - 0.5 && y > -0.5 && y < mesh.Ny - 0.5,
                        context + ": root center outside receiving background lattice");
                const int ix = static_cast<int>(std::round(x)), iy = static_cast<int>(std::round(y));
                require(pointNear(cell.centerCoord(), mesh.minCoord + Vector2D((ix + 0.5) * mesh.cellSize.x,
                        (iy + 0.5) * mesh.cellSize.y)), context + ": root center off receiving background lattice");
                auto& occupied = roots[static_cast<std::size_t>(iy) * mesh.Nx + ix];
                require(!occupied, context + ": duplicate root location");
                occupied = 1;
            }
            else
            {
                const auto& parent = cells[cell.parent];
                require(cell.parent != label && parent.label() == cell.parent &&
                        parent.level() == cell.level() - 1, context + ": invalid parent/level backlink");
                bool found = false;
                for (int child : parent.son) found = found || child == label;
                require(found, context + ": parent does not contain child");
            }
            if (!leaf)
            {
                for (int corner = 0; corner < 4; ++corner)
                {
                    const int childLabel = cell.son[corner];
                    const auto& child = cells[childLabel];
                    require(childLabel != label && child.label() == childLabel && child.parent == label &&
                            child.level() == cell.level() + 1, context + ": invalid child quartet/backlink");
                    const Vector2D expected = cell.centerCoord() + Vector2D(
                        (corner == 1 || corner == 2 ? 0.25 : -0.25) * cell.size().x,
                        (corner >= 2 ? 0.25 : -0.25) * cell.size().y);
                    require(pointNear(child.centerCoord(), expected), context + ": child center does not quarter parent");
                }
            }
        }
        // Active incidence is stronger than a label-range check. Parent edges
        // retain historical owner information and are intentionally excluded.
        for (int label = 0; label < capacity[3]; ++label)
        {
            const auto& cell = cells[label];
            if (cell.label() == -1 || cell.son[0] != -1) continue;
            const std::string context = "active cell " + std::to_string(label);
            for (bool horizontal : {false, true})
            {
                const auto* edges = horizontal ? hEdges.get() : vEdges.get();
                const int kind = horizontal ? 1 : 2;
                for (int side = 0; side < 2; ++side)
                {
                    const int* face = horizontal ? cell.cellToHEdge[side] : cell.cellToVEdge[side];
                    require(face[0] >= 0, context + ": required first face edge is missing");
                    const int pieces = face[1] == -1 ? 1 : 2;
                    require(pieces == 1 || face[0] != face[1], context + ": duplicate split face edge");
                    Vector2D firstPoint, lastPoint;
                    for (int part = 0; part < pieces; ++part)
                    {
                        const int edgeLabel = face[part];
                        require(edgeLabel >= 0 && edgeLabel < capacity[kind] &&
                                edges[edgeLabel].label() == edgeLabel && edges[edgeLabel].son[0] == -1,
                                context + ": required face does not reference a live leaf edge");
                        const auto& edge = edges[edgeLabel];
                        require((side == 0 ? edge.masterCell : edge.slaveCell) == label,
                                context + ": face/edge owner incidence is inconsistent");
                        const auto start = vertices[edge.edgeToVert[0]].coord();
                        const auto finish = vertices[edge.edgeToVert[1]].coord();
                        if (part == 0) firstPoint = start;
                        else require(edge.edgeToVert[0] == edges[face[0]].edgeToVert[1],
                                     context + ": split face pieces do not meet in order");
                        lastPoint = finish;
                        const int direction = horizontal ? (side == 0 ? 1 : 3) : (side == 0 ? 0 : 2);
                        const int neighbor = side == 0 ? edge.slaveCell : edge.masterCell;
                        require(cell.neiCell[direction][part] == neighbor,
                                context + ": face/neighbor incidence is inconsistent");
                    }
                    const auto center = cell.centerCoord();
                    const auto size = cell.size();
                    const Vector2D expectedStart = horizontal ?
                        center + Vector2D(-0.5 * size.x, (side == 0 ? -0.5 : 0.5) * size.y) :
                        center + Vector2D((side == 0 ? -0.5 : 0.5) * size.x, -0.5 * size.y);
                    const Vector2D expectedEnd = horizontal ? expectedStart + Vector2D(size.x, 0.0) :
                                                                            expectedStart + Vector2D(0.0, size.y);
                    require(pointNear(firstPoint, expectedStart) && pointNear(lastPoint, expectedEnd),
                            context + ": face edges do not cover the correct geometric side");
                }
            }
        }
        for (int kind : {1, 2})
        {
            const bool horizontal = kind == 1;
            const auto* edges = horizontal ? hEdges.get() : vEdges.get();
            for (int label = 0; label < capacity[kind]; ++label)
            {
                const auto& edge = edges[label];
                if (edge.label() == -1 || edge.son[0] != -1) continue;
                const std::string context = std::string(horizontal ? "active horizontal edge " : "active vertical edge ") +
                    std::to_string(label);
                require(edge.masterCell != -1 || edge.slaveCell != -1, context + ": no live owner");
                require(edge.masterCell == -1 || edge.masterCell != edge.slaveCell, context + ": duplicate owners");
                const auto& start = vertices[edge.edgeToVert[0]];
                const auto& finish = vertices[edge.edgeToVert[1]];
                require(start.son[0] == -1 && finish.son[0] == -1 &&
                        (horizontal ? start.vertToHEdge[1] : start.vertToVEdge[1]) == label &&
                        (horizontal ? finish.vertToHEdge[0] : finish.vertToVEdge[0]) == label &&
                        start.neiVertices[horizontal ? 2 : 3] == edge.edgeToVert[1] &&
                        finish.neiVertices[horizontal ? 0 : 1] == edge.edgeToVert[0],
                        context + ": endpoint edge/necessary-neighbor incidence is inconsistent");
                for (int side = 0; side < 2; ++side)
                {
                    const int owner = side == 0 ? edge.masterCell : edge.slaveCell;
                    if (owner == -1)
                    {
                        const double position = horizontal ? edge.midpointCoord().y : edge.midpointCoord().x;
                        const double limit = horizontal ? (side == 0 ? mesh.maxCoord.y : mesh.minCoord.y) :
                                                         (side == 0 ? mesh.maxCoord.x : mesh.minCoord.x);
                        require(near(position, limit, horizontal ? sy : sx), context + ": missing interior owner");
                        continue;
                    }
                    require(owner >= 0 && owner < capacity[3] && cells[owner].label() == owner && cells[owner].son[0] == -1,
                            context + ": owner is not a live leaf cell");
                    const auto& cell = cells[owner];
                    const int* face = horizontal ? cell.cellToHEdge[side] : cell.cellToVEdge[side];
                    require(face[0] == label || face[1] == label, context + ": owner has no matching face backlink");
                }
            }
        }
        // Infer necessary neighbors from actual edge incidence, not cached
        // Boundary/TVertType: static refineBoxes historically precedes type
        // classification. A true T-junction may therefore have type=non here.
        // Classification refresh remains a separate caller responsibility.
        for (int label = 0; label < capacity[0]; ++label)
        {
            const auto& value = vertices[label];
            if (value.label() == -1 || value.son[0] != -1) continue;
            const std::string context = "active vertex " + std::to_string(label);
            int degree = 0;
            for (int direction = 0; direction < 4; ++direction)
            {
                const bool horizontal = direction == 0 || direction == 2;
                const int side = direction >= 2 ? 1 : 0;
                const int edgeLabel = horizontal ? value.vertToHEdge[side] : value.vertToVEdge[side];
                const int neighbor = value.neiVertices[direction];
                require((edgeLabel == -1) == (neighbor == -1), context + ": edge/neighbor missingness is inconsistent");
                if (edgeLabel == -1) continue;
                ++degree;
                const auto* edges = horizontal ? hEdges.get() : vEdges.get();
                const int kind = horizontal ? 1 : 2;
                require(edgeLabel >= 0 && edgeLabel < capacity[kind] && edges[edgeLabel].label() == edgeLabel &&
                        edges[edgeLabel].son[0] == -1 && neighbor >= 0 && neighbor < capacity[0] &&
                        vertices[neighbor].label() == neighbor && vertices[neighbor].son[0] == -1,
                        context + ": necessary neighbor/edge is not live");
                const auto& edge = edges[edgeLabel];
                require(edge.edgeToVert[side == 1 ? 0 : 1] == label &&
                        edge.edgeToVert[side == 1 ? 1 : 0] == neighbor,
                        context + ": directed edge/neighbor endpoints are inconsistent");
            }
            const auto point = value.coord();
            const bool verticalBoundary = near(point.x, mesh.minCoord.x, sx) || near(point.x, mesh.maxCoord.x, sx);
            const bool horizontalBoundary = near(point.y, mesh.minCoord.y, sy) || near(point.y, mesh.maxCoord.y, sy);
            const int boundaryCount = static_cast<int>(verticalBoundary) + static_cast<int>(horizontalBoundary);
            require(boundaryCount == 2 ? degree == 2 : (boundaryCount == 1 ? degree == 3 : degree == 3 || degree == 4),
                    context + ": incomplete geometric neighbor incidence (not a domain corner/boundary or valid T-junction)");
        }
        require(std::all_of(roots.begin(), roots.end(), [](unsigned char value) { return value != 0; }),
                "cache root cells do not cover the complete receiving domain");
        require(observedCounts == cellsPerLevel, "cached per-level leaf counts are inconsistent");
        for (const auto* levels : {&minNegative, &minPositive, &maxLocal})
            for (int level : *levels)
                require(level >= 0 && level <= maxLevel, "cached refinement configuration exceeds supported levels");
        validateDyadicCoordinates(mesh);
    }
};
} // namespace QuadtreeCacheDetail
