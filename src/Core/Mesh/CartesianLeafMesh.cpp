#include "Core/Mesh/CartesianLeafMesh.h"
#include "Core/Mesh/DyadicGridCoordinates.h"
#include "Core/Mesh/QuadtreeMesh.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
void require(bool value, const char *message)
{
    if (!value)
        throw std::invalid_argument(std::string("[CartesianLeafMesh] ") + message);
}
bool near(double a, double b, double scale)
{
    return std::isfinite(a) && std::isfinite(b) &&
           std::abs(a - b) <= 256.0 * std::numeric_limits<double>::epsilon() * scale;
}
Boundary boundaryAt(bool w, bool s, bool e, bool n)
{
    if (s)
        return w ? Boundary::southWest : e ? Boundary::southEast : Boundary::south;
    if (n)
        return w ? Boundary::northWest : e ? Boundary::northEast : Boundary::north;
    return w ? Boundary::west : e ? Boundary::east : Boundary::non;
}
Boundary cellBoundaryAt(bool w, bool s, bool e, bool n)
{
    Boundary result = w ? Boundary::west : Boundary::non;
    if (e)
        result = Boundary::east;
    if (s)
        result = result == Boundary::west   ? Boundary::southWest
                 : result == Boundary::east ? Boundary::southEast
                                            : Boundary::south;
    if (n)
        result = result == Boundary::west   ? Boundary::northWest
                 : result == Boundary::east ? Boundary::northEast
                                            : Boundary::north;
    return result;
}
} // namespace

CartesianLeafMetrics validateCartesianLeafMesh(const QuadtreeMesh &mesh)
{
    const int nc = mesh.cellIndexNum(), nv = mesh.vertIndexNum();
    const int ne = mesh.vEdgeIndexNum(), nh = mesh.hEdgeIndexNum();
    require(nc > 0 && nv > 0 && ne > 0 && nh > 0 && mesh.cellIndexToLabel && mesh.vertIndexToLabel &&
                mesh.vEdgeIndexToLabel && mesh.hEdgeIndexToLabel,
            "missing native storage");
    require(mesh.cellLabelCapacity() >= nc && mesh.vertLabelCapacity() >= nv &&
                mesh.vEdgeLabelCapacity() >= ne && mesh.hEdgeLabelCapacity() >= nh,
            "undersized native storage");
    const long long roots = static_cast<long long>(mesh.Nx) * mesh.Ny;
    require(mesh.Nx > 0 && mesh.Ny > 0 && roots <= nc, "invalid root lattice dimensions");
    const LAMP::Mesh::DyadicGridCoordinates coordinates(mesh.minCoord, mesh.maxCoord, mesh.Nx, mesh.Ny);
    const auto cellAt = [&](int label) -> const Cell & {
        require(label >= 0 && label < mesh.cellLabelCapacity(), "invalid cell reference");
        const auto &c = mesh.cells(label);
        require(c.label() == label, "cell label identity mismatch");
        return c;
    };
    const auto vertexAt = [&](int label) -> const Vert & {
        require(label >= 0 && label < mesh.vertLabelCapacity(), "invalid vertex reference");
        const auto &v = mesh.verts(label);
        require(v.label() == label && v.index() >= 0 && v.index() < nv &&
                    mesh.vertIndexToLabel[v.index()] == label,
                "vertex index identity mismatch");
        return v;
    };
    std::vector<unsigned char> leaf(mesh.cellLabelCapacity(), 0), visited(leaf.size(), 0);
    CartesianLeafMetrics result;
    for (int i = 0; i < nc; ++i)
    {
        const auto &c = cellAt(mesh.cellIndexToLabel[i]);
        require(c.index() == i && !leaf[c.label()], "duplicate/misindexed leaf cell");
        for (int son : c.son)
            require(son == -1, "indexed cell is not a leaf");
        leaf[c.label()] = 1;
        const auto address = LAMP::Mesh::cellAddress(mesh, c);
        const auto h = coordinates.cellSize(address.level);
        const auto center = coordinates.point(2 * address.i + 1, 2 * address.j + 1, address.level + 1);
        const double sx = std::max({std::abs(mesh.minCoord.x), std::abs(mesh.maxCoord.x), h.x});
        const double sy = std::max({std::abs(mesh.minCoord.y), std::abs(mesh.maxCoord.y), h.y});
        require(h.x > 1024 * std::numeric_limits<double>::epsilon() * sx &&
                    h.y > 1024 * std::numeric_limits<double>::epsilon() * sy,
                "unresolved leaf spacing");
        require(near(c.size().x, h.x, h.x) && near(c.size().y, h.y, h.y) &&
                    near(c.centerCoord().x, center.x, sx) && near(c.centerCoord().y, center.y, sy),
                "leaf geometry disagrees with native dyadic ancestry");
        const std::array<Vector2D, 4> corners{{coordinates.point(address.i, address.j, address.level),
                                               coordinates.point(address.i + 1, address.j, address.level),
                                               coordinates.point(address.i + 1, address.j + 1, address.level),
                                               coordinates.point(address.i, address.j + 1, address.level)}};
        for (int k = 0; k < 4; ++k)
        {
            const auto p = vertexAt(c.cellToVert[k]).coord();
            require(near(p.x, corners[k].x, sx) && near(p.y, corners[k].y, sy), "wrong native cell corner");
        }
        require(c.boundary() ==
                    cellBoundaryAt(corners[0].x == mesh.minCoord.x, corners[0].y == mesh.minCoord.y,
                                   corners[2].x == mesh.maxCoord.x, corners[2].y == mesh.maxCoord.y),
                "stale cell boundary metadata");
        if (i == 0)
            result.minimumSpacing = h;
        else
        {
            result.minimumSpacing.x = std::min(result.minimumSpacing.x, h.x);
            result.minimumSpacing.y = std::min(result.minimumSpacing.y, h.y);
            result.mixedLevels |= c.level() != mesh.cells(mesh.cellIndexToLabel[0]).level();
        }
    }
    // Traverse the original forest, not an area-sum surrogate: no omitted,
    // overlapping, duplicated or disconnected indexed leaf is accepted.
    std::vector<int> pending;
    for (int root = 0; root < roots; ++root)
        pending.push_back(root);
    int reached = 0;
    while (!pending.empty())
    {
        const int label = pending.back();
        pending.pop_back();
        const auto &c = cellAt(label);
        require(!visited[label], "cyclic/duplicate quadtree child");
        visited[label] = 1;
        if (c.son[0] == -1)
        {
            require(leaf[label] != 0, "omitted leaf in root coverage");
            ++reached;
            for (int son : c.son)
                require(son == -1, "partial child array");
        }
        else
            for (int son : c.son)
            {
                const auto &child = cellAt(son);
                require(child.parent == label && child.level() == c.level() + 1, "invalid child ancestry");
                pending.push_back(son);
            }
    }
    require(reached == nc, "indexed leaves do not cover the original roots exactly");
    const auto edgeAt = [&](bool horizontal, int label) -> const Edge & {
        require(label >= 0 && label < (horizontal ? mesh.hEdgeLabelCapacity() : mesh.vEdgeLabelCapacity()),
                "invalid edge reference");
        const auto &e = horizontal ? mesh.hEdges(label) : mesh.vEdges(label);
        const int count = horizontal ? nh : ne;
        const auto mapping = horizontal ? mesh.hEdgeIndexToLabel : mesh.vEdgeIndexToLabel;
        require(e.label() == label && e.index() >= 0 && e.index() < count && mapping[e.index()] == label &&
                    e.isHorizontal() == horizontal && e.son[0] == -1 && e.son[1] == -1,
                "invalid real leaf edge");
        return e;
    };
    std::vector<std::array<int, 4>> incident(nv, {{-1, -1, -1, -1}}), quadrants(nv, {{-1, -1, -1, -1}});
    for (bool horizontal : {false, true})
        for (int i = 0; i < (horizontal ? nh : ne); ++i)
        {
            const auto &e =
                edgeAt(horizontal, (horizontal ? mesh.hEdgeIndexToLabel : mesh.vEdgeIndexToLabel)[i]);
            require(e.index() == i, "misindexed native edge");
            const auto &a = vertexAt(e.edgeToVert[0]);
            const auto &b = vertexAt(e.edgeToVert[1]);
            const auto p = a.coord(), q = b.coord(), m = e.midpointCoord();
            const double length = horizontal ? q.x - p.x : q.y - p.y;
            require(e.masterCell >= 0 || e.slaveCell >= 0, "edge has no leaf owner");
            int edgeLevel = -1;
            for (int label : {e.masterCell, e.slaveCell})
                if (label != -1)
                {
                    const auto& owner = cellAt(label);
                    require(leaf[label] != 0, "edge references a non-leaf owner");
                    edgeLevel = std::max(edgeLevel, owner.level());
                }
            const auto nominal = coordinates.cellSize(static_cast<unsigned>(edgeLevel));
            // Endpoints have already been checked against the same dyadic
            // ancestry as the cells. Their rounded difference loses digits
            // for small leaves far from zero; it is not the nominal metric.
            // Check the stored length against the finest owner's dyadic size,
            // keeping the original relative metric tolerance and all incidence
            // and coordinate checks (no coordinate-scaled length tolerance).
            const double nominalLength = horizontal ? nominal.x : nominal.y;
            const double scale =
                std::max({std::abs(p.x), std::abs(p.y), std::abs(q.x), std::abs(q.y), length});
            require(length > 0 && near(nominalLength, e.length(), e.length()) &&
                        (horizontal ? p.y == q.y : p.x == q.x) && near(m.x, .5 * (p.x + q.x), scale) &&
                        near(m.y, .5 * (p.y + q.y), scale),
                    "invalid edge geometry");
            if (e.masterCell >= 0 && e.slaveCell >= 0)
                require(std::abs(cellAt(e.masterCell).level() - cellAt(e.slaveCell).level()) <= 1,
                        "non-2:1 face");
            const Boundary expected = e.slaveCell < 0    ? (horizontal ? Boundary::south : Boundary::west)
                                      : e.masterCell < 0 ? (horizontal ? Boundary::north : Boundary::east)
                                                         : Boundary::non;
            require(e.boundary() == expected, "stale edge boundary metadata");
            if (e.slaveCell < 0)
                require(horizontal ? m.y == mesh.minCoord.y : m.x == mesh.minCoord.x,
                        "interior missing negative owner");
            if (e.masterCell < 0)
                require(horizontal ? m.y == mesh.maxCoord.y : m.x == mesh.maxCoord.x,
                        "interior missing positive owner");
            const int forward = horizontal ? 2 : 3, backward = horizontal ? 0 : 1;
            require(incident[a.index()][forward] == -1 && incident[b.index()][backward] == -1,
                    "duplicate vertex edge direction");
            incident[a.index()][forward] = e.label();
            incident[b.index()][backward] = e.label();
        }
    for (int i = 0; i < nc; ++i)
    {
        const auto &c = cellAt(mesh.cellIndexToLabel[i]);
        for (int side = 0; side < 4; ++side)
        {
            const bool horizontal = side == 1 || side == 3;
            const int slot = side < 2 ? 0 : 1;
            const auto &labels = horizontal ? c.cellToHEdge[slot] : c.cellToVEdge[slot];
            require(labels[0] >= 0, "missing complete cell side");
            double cursor =
                horizontal ? vertexAt(c.cellToVert[0]).coord().x : vertexAt(c.cellToVert[0]).coord().y;
            const double end =
                horizontal ? vertexAt(c.cellToVert[2]).coord().x : vertexAt(c.cellToVert[2]).coord().y;
            const double normal = horizontal ? vertexAt(c.cellToVert[slot ? 2 : 0]).coord().y
                                             : vertexAt(c.cellToVert[slot ? 2 : 0]).coord().x;
            for (int part = 0; part < 2; ++part)
            {
                if (labels[part] == -1)
                {
                    require(c.neiCell[side][part] == -1, "neighbor without side segment");
                    continue;
                }
                const auto &e = edgeAt(horizontal, labels[part]);
                require((slot ? e.slaveCell : e.masterCell) == c.label() &&
                            c.neiCell[side][part] == (slot ? e.masterCell : e.slaveCell),
                        "wrong cell/edge/neighbor incidence");
                const auto p = vertexAt(e.edgeToVert[0]).coord(), q = vertexAt(e.edgeToVert[1]).coord();
                require((horizontal ? p.y : p.x) == normal && (horizontal ? p.x : p.y) == cursor,
                        "gapped/overlapping cell side");
                cursor = horizontal ? q.x : q.y;
                for (int label : e.edgeToVert)
                {
                    const auto &v = vertexAt(label);
                    const auto point = v.coord();
                    const double dx = c.centerCoord().x - point.x, dy = c.centerCoord().y - point.y;
                    const bool midX = near(dx, 0.0, c.size().x), midY = near(dy, 0.0, c.size().y);
                    for (int quadrant = 0; quadrant < 4; ++quadrant)
                    {
                        const bool right = quadrant == 1 || quadrant == 2, up = quadrant >= 2;
                        if ((!midX && right != (dx > 0)) || (!midY && up != (dy > 0)))
                            continue;
                        auto &owner = quadrants[v.index()][quadrant];
                        require(owner == -1 || owner == c.label(), "overlapping cells at vertex");
                        owner = c.label();
                    }
                }
            }
            require(cursor == end, "incomplete cell side coverage");
        }
    }
    for (int i = 0; i < nv; ++i)
    {
        const auto &v = vertexAt(mesh.vertIndexToLabel[i]);
        const auto p = v.coord();
        require(v.index() == i, "misindexed native vertex");
        const Boundary boundary = boundaryAt(p.x == mesh.minCoord.x, p.y == mesh.minCoord.y,
                                             p.x == mesh.maxCoord.x, p.y == mesh.maxCoord.y);
        require(v.boundary() == boundary && v.vertToHEdge[0] == incident[i][0] &&
                    v.vertToVEdge[0] == incident[i][1] && v.vertToHEdge[1] == incident[i][2] &&
                    v.vertToVEdge[1] == incident[i][3],
                "stale vertex boundary/edge metadata");
        int missing = -1, missingCount = 0;
        for (int k = 0; k < 4; ++k)
        {
            require(v.vertToCell[k] == quadrants[i][k], "wrong vertex cell quadrant");
            if (incident[i][k] < 0)
            {
                missing = k;
                ++missingCount;
            }
        }
        if (boundary == Boundary::non && missingCount)
        {
            require(missingCount == 1, "invalid interior hanging vertex degree");
            const int owner = quadrants[i][missing]; // W/S/E/N -> SW/SE/NE/NW
            require(v.supCell == owner && owner >= 0, "missing hanging coarse-cell identity");
            const auto &coarse = cellAt(owner);
            const std::array<TVertType, 4> types{
                {TVertType::west, TVertType::south, TVertType::east, TVertType::north}};
            const int first = missing == 0 ? 0 : missing == 1 ? 0 : missing == 2 ? 1 : 3;
            const int last = missing == 0 ? 3 : missing == 1 ? 1 : 2;
            require(v.type() == types[missing] && v.supVert[0] == coarse.cellToVert[first] &&
                        v.supVert[1] == coarse.cellToVert[last],
                    "stale hanging interpolation metadata");
        }
        else
            require(v.type() == TVertType::non && v.supCell == -1 && v.supVert[0] == -1 && v.supVert[1] == -1,
                    "unexpected hanging metadata on a regular/boundary vertex");
    }
    return result;
}
