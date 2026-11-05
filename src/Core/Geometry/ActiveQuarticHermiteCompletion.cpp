#include "Core/Geometry/ActiveQuarticHermiteCompletion.h"

#include "Core/Geometry/Detail/ExactExpansion.h"
#include "Core/Mesh/NativeMeshIdentity.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/UniformLeafGrid.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

namespace LAMP::Geometry {
namespace {
using Completion = ActiveQuarticHermiteCompletion;
using Weight = Completion::SourceWeight;
using Kind = Weight::Kind;
using EdgeKey = std::array<int, 2>;

void require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::invalid_argument("[ActiveQuarticHermiteCompletion] " + message);
}

// Attribute delegated topology guards to this owner while preserving the
// original validation reason. Allocation and unrelated exceptions propagate.
template<class F>
auto delegatedValidation(const char* stage, F&& action) -> decltype(action())
{
    try {
        return action();
    } catch (const std::invalid_argument& error) {
        throw std::invalid_argument(std::string("[ActiveQuarticHermiteCompletion] ")
            + stage + ": " + error.what());
    }
}

struct SourceKey {
    Kind kind = Kind::Vertex;
    int index = -1, component = 0;
    bool operator<(const SourceKey& other) const
    {
        return std::make_tuple(static_cast<int>(kind), index, component)
            < std::make_tuple(static_cast<int>(other.kind), other.index, other.component);
    }
};

struct IntegerTerm {
    SourceKey source;
    std::int64_t numerator = 0;
};

// Selectors retain their source value exactly. Derived rows retain the integer
// operator before division, with a true source reference for constant removal.
struct OperatorRow {
    std::vector<IntegerTerm> terms;
    SourceKey reference;
    long double denominator = 1.0L;
    bool derived = false;
};

double sourceValue(const ActiveQuarticHermitePatch& patch, const SourceKey& key)
{
    switch (key.kind) {
    case Kind::Vertex: {
        const auto& node = patch.nodes.at(static_cast<std::size_t>(key.index));
        switch (key.component) {
        case 0: return node.value;
        case 1: return node.dx;
        case 2: return node.dy;
        case 3: return node.dxy;
        default: break;
        }
        break;
    }
    case Kind::EdgeMidpoint: {
        const auto& edge = patch.edgeMidpoints.at(static_cast<std::size_t>(key.index));
        if (key.component == 0) return edge.value;
        if (key.component == 1) return edge.normalDerivative;
        break;
    }
    case Kind::CellCenter:
        if (key.component == 0)
            return patch.cellCenters.at(static_cast<std::size_t>(key.index));
        break;
    }
    throw std::logic_error("[ActiveQuarticHermiteCompletion] invalid retained source identity");
}

double evaluate(const OperatorRow& row, const ActiveQuarticHermitePatch& patch)
{
    if (row.terms.empty()) return 0.0;
    if (!row.derived) {
        require(row.terms.size() == 1 && row.terms.front().numerator == 1
            && row.denominator == 1.0L, "invalid retained selector");
        return sourceValue(patch, row.terms.front().source);
    }
    const long double reference = sourceValue(patch, row.reference);
    long double sum = 0.0L, correction = 0.0L;
    for (const auto& term : row.terms) {
        const long double difference =
            static_cast<long double>(sourceValue(patch, term.source)) - reference;
        const long double contribution =
            static_cast<long double>(term.numerator) * difference;
        const long double next = sum + contribution;
        // Fixed source ordering and compensated contraction apply identically
        // to the completed values and to admissible datum rates.
        correction += std::abs(sum) >= std::abs(contribution)
            ? (sum - next) + contribution : (contribution - next) + sum;
        sum = next;
    }
    const long double physical = (sum + correction) / row.denominator;
    const double result = static_cast<double>(physical);
    require(std::isfinite(physical) && std::isfinite(result)
        && (physical == 0.0L || result != 0.0),
        "completed derivative is not finitely representable");
    return result;
}

void validateFinite(const ActiveQuarticHermitePatch& patch)
{
    for (const auto& node : patch.nodes)
        require(std::isfinite(node.value) && std::isfinite(node.dx)
            && std::isfinite(node.dy) && std::isfinite(node.dxy),
            "nonfinite formal vertex data");
    for (const auto& edge : patch.edgeMidpoints)
        require(std::isfinite(edge.value) && std::isfinite(edge.normalDerivative),
            "nonfinite formal edge data");
    for (double value : patch.cellCenters)
        require(std::isfinite(value), "nonfinite actual cell-center data");
}

void validateLayout(const ActiveQuarticHermitePatch& patch,
    const ActiveQuarticHermitePatch& source)
{
    require(patch.activeCells == source.activeCells && patch.nodes.size() == source.nodes.size()
        && patch.edgeMidpoints.size() == source.edgeMidpoints.size()
        && patch.cellCenters.size() == source.cellCenters.size()
        && patch.cornerSectors.empty(), "formal patch layout differs from retained source");
    for (std::size_t edge = 0; edge < source.edgeMidpoints.size(); ++edge)
        require(patch.edgeMidpoints[edge].vertices == source.edgeMidpoints[edge].vertices,
            "formal edge order or native identity differs from retained source");
    validateFinite(patch);
}

// Integer coefficients for quartic interpolation differentiated twice, at
// target 0..4 within a consecutive five-point half-lattice window.
constexpr std::array<std::array<std::int64_t, 5>, 5> secondNumerators{{
    {{35, -104, 114, -56, 11}},
    {{11, -20, 6, 4, -1}},
    {{-1, 16, -30, 16, -1}},
    {{-1, 4, 6, -20, 11}},
    {{11, -56, 114, -104, 35}}
}};

struct Stencil {
    int begin = 0;
    std::array<std::int64_t, 5> numerator{};
};

Stencil secondStencil(int target, int count)
{
    require(count >= 5 && target >= 0 && target < count,
        "missing five-point half-lattice support");
    const int begin = std::clamp(target - 2, 0, count - 5);
    return {begin, secondNumerators.at(static_cast<std::size_t>(target - begin))};
}

std::vector<Weight> mergedPullback(
    const std::vector<Completion::Row>& rows,
    const std::vector<std::pair<std::size_t, long double>>& coefficients)
{
    std::map<SourceKey, long double> accumulated;
    for (const auto& coefficient : coefficients) {
        require(coefficient.first < rows.size() && std::isfinite(coefficient.second),
            "invalid completed-row pullback coefficient");
        if (coefficient.second == 0.0L) continue;
        for (const auto& weight : rows[coefficient.first].sourceWeights) {
            const SourceKey key{weight.kind, weight.index, weight.component};
            const long double contribution = coefficient.second * weight.coefficient;
            require(std::isfinite(contribution), "nonfinite formal-source pullback");
            const long double sum = accumulated[key] + contribution;
            require(std::isfinite(sum), "nonfinite merged formal-source pullback");
            accumulated[key] = sum;
        }
    }
    std::vector<Weight> result;
    result.reserve(accumulated.size());
    for (const auto& entry : accumulated)
        if (entry.second != 0.0L)
            result.push_back({entry.first.kind, entry.first.index,
                entry.first.component, entry.second});
    return result;
}
} // namespace

struct ActiveQuarticHermiteCompletion::Data {
    NativeMeshIdentity mesh;
    ActiveQuarticHermitePatch source;
    std::vector<SymmetryPlane> symmetry;
    std::vector<Row> rows;
    std::vector<OperatorRow> operators;
    std::vector<std::array<std::size_t, 49>> cellRows;
    std::vector<std::array<bool, 2>> reflectedVertices;
    std::vector<bool> reflectedEdges;

    Data(const QuadtreeMesh& native, const ActiveQuarticHermitePatch& patch,
        const std::vector<SymmetryPlane>& planes)
        : mesh(delegatedValidation("mesh identity construction",
              [&] { return NativeMeshIdentity(native); })),
          source(patch), symmetry(planes)
    {
        Detail::ExactExpansion::validateEnvironment();
        const auto grid = delegatedValidation("uniform leaf layout",
            [&] { return UniformLeafGrid(native, CartesianGridScope::UniformLeaves); });
        const int nx = grid.nx(), ny = grid.ny();
        require(nx >= 2 && ny >= 2, "requires at least two native cells in each axis");
        require(nx <= (std::numeric_limits<int>::max() - 1) / 2
            && ny <= (std::numeric_limits<int>::max() - 1) / 2,
            "half-lattice dimensions exceed index range");
        require(source.activeCells.size() == static_cast<std::size_t>(native.cellIndexNum())
            && source.nodes.size() == static_cast<std::size_t>(native.vertIndexNum())
            && source.cellCenters.size() == static_cast<std::size_t>(native.cellIndexNum())
            && source.edgeMidpoints.size() ==
                static_cast<std::size_t>(native.hEdgeIndexNum())
                    + static_cast<std::size_t>(native.vEdgeIndexNum()),
            "requires complete actual native formal data");
        require(source.edgeMidpoints.size() <=
            static_cast<std::size_t>(std::numeric_limits<int>::max()),
            "formal edge source ordinals exceed index range");
        require(source.cornerSectors.empty(), "corner-sector completion is unsupported");
        for (auto active : source.activeCells)
            require(active == 1, "requires all actual native cells active");
        validateFinite(source);

        std::vector<int> vertices(grid.vertices().size()), cells(grid.cells().size());
        for (std::size_t i = 0; i < vertices.size(); ++i)
            vertices[i] = native.verts(grid.vertices()[i]).index();
        for (std::size_t i = 0; i < cells.size(); ++i)
            cells[i] = native.cells(grid.cells()[i]).index();
        const auto vertex = [&](int x, int y) {
            return vertices.at(static_cast<std::size_t>(y) * (nx + 1) + x);
        };
        const auto cell = [&](int x, int y) {
            return cells.at(static_cast<std::size_t>(y) * nx + x);
        };

        struct Face {
            bool horizontal = false;
            std::size_t latticeSlot = 0;
            int sourceOrdinal = -1;
        };
        std::map<EdgeKey, Face> actualFaces;
        std::vector<int> horizontal(grid.horizontal().size(), -1);
        std::vector<int> vertical(grid.vertical().size(), -1);
        const auto insertFace = [&](const Edge& edge, bool horizontalFace, std::size_t slot) {
            const EdgeKey key{{
                native.verts(edge.edgeToVert[0]).index(),
                native.verts(edge.edgeToVert[1]).index()
            }};
            require(actualFaces.emplace(key, Face{horizontalFace, slot, -1}).second,
                "duplicate actual canonical native face");
        };
        for (std::size_t i = 0; i < grid.horizontal().size(); ++i)
            insertFace(native.hEdges(grid.horizontal()[i]), true, i);
        for (std::size_t i = 0; i < grid.vertical().size(); ++i)
            insertFace(native.vEdges(grid.vertical()[i]), false, i);
        for (std::size_t i = 0; i < source.edgeMidpoints.size(); ++i) {
            const auto& key = source.edgeMidpoints[i].vertices;
            require(key[0] >= 0 && key[1] >= 0
                && key[0] < native.vertIndexNum() && key[1] < native.vertIndexNum(),
                "formal edge contains a foreign native vertex");
            const auto found = actualFaces.find(key);
            require(found != actualFaces.end(),
                "formal edge is not an actual canonical complete native face");
            auto& face = found->second;
            require(face.sourceOrdinal == -1, "duplicate formal native face identity");
            face.sourceOrdinal = static_cast<int>(i);
            (face.horizontal ? horizontal : vertical).at(face.latticeSlot) = face.sourceOrdinal;
        }
        for (const auto& face : actualFaces)
            require(face.second.sourceOrdinal >= 0, "missing actual native face data");

        reflectedVertices.assign(source.nodes.size(), {{false, false}});
        reflectedEdges.assign(source.edgeMidpoints.size(), false);
        const auto coordinate = [](Vector2D point, int axis) {
            return axis == 0 ? point.x : point.y;
        };
        for (std::size_t i = 0; i < symmetry.size(); ++i) {
            const auto& plane = symmetry[i];
            require((plane.normalAxis == 0 || plane.normalAxis == 1)
                && std::isfinite(plane.coordinate), "invalid scalar symmetry plane");
            require(plane.coordinate == coordinate(native.minCoord, plane.normalAxis)
                || plane.coordinate == coordinate(native.maxCoord, plane.normalAxis),
                "scalar symmetry must be an exact outer native section");
            for (std::size_t j = 0; j < i; ++j)
                require(!(plane == symmetry[j]), "duplicate scalar symmetry plane");
            const bool lower =
                plane.coordinate == coordinate(native.minCoord, plane.normalAxis);
            for (int y = 0; y <= ny; ++y)
                for (int x = 0; x <= nx; ++x) {
                    const bool onPlane = plane.normalAxis == 0
                        ? x == (lower ? 0 : nx) : y == (lower ? 0 : ny);
                    if (!onPlane) continue;
                    const int v = vertex(x, y);
                    require(coordinate(native.verts(grid.vertices()[
                        static_cast<std::size_t>(y) * (nx + 1) + x]).coord(),
                        plane.normalAxis) == plane.coordinate,
                        "declared scalar symmetry lacks an exact native boundary vertex");
                    reflectedVertices[v][plane.normalAxis] = true;
                }
            bool found = false;
            for (std::size_t e = 0; e < source.edgeMidpoints.size(); ++e) {
                const auto& key = source.edgeMidpoints[e].vertices;
                if (!reflectedVertices[key[0]][plane.normalAxis]
                    || !reflectedVertices[key[1]][plane.normalAxis]) continue;
                // Masks can include the opposite wall from an earlier
                // declaration; the actual coordinate selects this section.
                const auto first = native.verts(native.vertIndexToLabel[key[0]]).coord();
                const auto last = native.verts(native.vertIndexToLabel[key[1]]).coord();
                if (coordinate(first, plane.normalAxis) != plane.coordinate
                    || coordinate(last, plane.normalAxis) != plane.coordinate) continue;
                reflectedEdges[e] = true;
                found = true;
            }
            require(found, "scalar symmetry plane has no actual boundary trace");
        }
        validateSymmetry(source);

        const long double dx = static_cast<long double>(grid.spacing().x) * 0.5L;
        const long double dy = static_cast<long double>(grid.spacing().y) * 0.5L;
        const long double denominatorX = 12.0L * dx * dx;
        const long double denominatorY = 12.0L * dy * dy;
        const long double denominatorXY = 144.0L * dx * dx * dy * dy;
        require(std::isfinite(denominatorX) && denominatorX > 0.0L
            && std::isfinite(denominatorY) && denominatorY > 0.0L
            && std::isfinite(denominatorXY) && denominatorXY > 0.0L,
            "half-lattice derivative scales are unresolved");

        const auto valueKey = [&](int x, int y) -> SourceKey {
            if ((x & 1) == 0 && (y & 1) == 0)
                return {Kind::Vertex, vertex(x / 2, y / 2), 0};
            if ((y & 1) == 0)
                return {Kind::EdgeMidpoint, horizontal.at(
                    static_cast<std::size_t>(y / 2) * nx + x / 2), 0};
            if ((x & 1) == 0)
                return {Kind::EdgeMidpoint, vertical.at(
                    static_cast<std::size_t>(y / 2) * (nx + 1) + x / 2), 0};
            return {Kind::CellCenter, cell(x / 2, y / 2), 0};
        };
        const auto gxKey = [&](int x, int y) -> SourceKey {
            require((x & 1) == 0, "gx source requires an integer-x line");
            if ((y & 1) == 0) return {Kind::Vertex, vertex(x / 2, y / 2), 1};
            return {Kind::EdgeMidpoint, vertical.at(
                static_cast<std::size_t>(y / 2) * (nx + 1) + x / 2), 1};
        };
        const auto gyKey = [&](int x, int y) -> SourceKey {
            require((y & 1) == 0, "gy source requires an integer-y line");
            if ((x & 1) == 0) return {Kind::Vertex, vertex(x / 2, y / 2), 2};
            return {Kind::EdgeMidpoint, horizontal.at(
                static_cast<std::size_t>(y / 2) * nx + x / 2), 1};
        };
        const auto selector = [](SourceKey key) {
            OperatorRow result;
            result.terms.push_back({key, 1});
            return result;
        };
        const auto second = [&](int x, int y, int axis, const auto& sourceKey) {
            OperatorRow result;
            result.derived = true;
            result.reference = sourceKey(x, y);
            result.denominator = axis == 0 ? denominatorX : denominatorY;
            const Stencil stencil = secondStencil(axis == 0 ? x : y,
                axis == 0 ? 2 * nx + 1 : 2 * ny + 1);
            for (int j = 0; j < 5; ++j)
                result.terms.push_back({sourceKey(axis == 0 ? stencil.begin + j : x,
                    axis == 1 ? stencil.begin + j : y), stencil.numerator[j]});
            return result;
        };
        const auto tensorSecond = [&](int x, int y) {
            OperatorRow result;
            result.derived = true;
            result.reference = valueKey(x, y);
            result.denominator = denominatorXY;
            const auto sx = secondStencil(x, 2 * nx + 1);
            const auto sy = secondStencil(y, 2 * ny + 1);
            // One unrounded tensor contraction to actual formal V sources.
            for (int j = 0; j < 5; ++j)
                for (int i = 0; i < 5; ++i)
                    result.terms.push_back({valueKey(sx.begin + i, sy.begin + j),
                        sx.numerator[i] * sy.numerator[j]});
            return result;
        };

        rows.reserve(source.nodes.size() * 9 + source.edgeMidpoints.size() * 3
            + source.cellCenters.size());
        operators.reserve(rows.capacity());
        std::vector<std::array<int, 2>> vertexLattice(source.nodes.size());
        std::vector<std::array<int, 2>> edgeLattice(source.edgeMidpoints.size());
        for (int y = 0; y <= ny; ++y)
            for (int x = 0; x <= nx; ++x)
                vertexLattice[vertex(x, y)] = {{2 * x, 2 * y}};
        for (int y = 0; y <= ny; ++y)
            for (int x = 0; x < nx; ++x)
                edgeLattice[horizontal[static_cast<std::size_t>(y) * nx + x]]
                    = {{2 * x + 1, 2 * y}};
        for (int y = 0; y < ny; ++y)
            for (int x = 0; x <= nx; ++x)
                edgeLattice[vertical[static_cast<std::size_t>(y) * (nx + 1) + x]]
                    = {{2 * x, 2 * y + 1}};

        for (int v = 0; v < native.vertIndexNum(); ++v) {
            const int x = vertexLattice[v][0], y = vertexLattice[v][1];
            for (int component = 0; component < 9; ++component) {
                OperatorRow op;
                const bool constrained =
                    (component == 1 || component == 7) && reflectedVertices[v][0]
                    || (component == 2 || component == 6) && reflectedVertices[v][1]
                    || component == 3 && (reflectedVertices[v][0] || reflectedVertices[v][1]);
                if (component < 4) op = selector({Kind::Vertex, v, component});
                else if (constrained) op.derived = true; // exact added odd row
                else if (component == 4) op = second(x, y, 0, valueKey);
                else if (component == 5) op = second(x, y, 1, valueKey);
                else if (component == 6) op = second(x, y, 0, gyKey);
                else if (component == 7) op = second(x, y, 1, gxKey);
                else op = tensorSecond(x, y);
                appendRow(Kind::Vertex, v, component, constrained, std::move(op));
            }
        }
        for (std::size_t edge = 0; edge < source.edgeMidpoints.size(); ++edge) {
            const int x = edgeLattice[edge][0], y = edgeLattice[edge][1];
            for (int component = 0; component < 3; ++component) {
                OperatorRow op = component < 2
                    ? selector({Kind::EdgeMidpoint, static_cast<int>(edge), component})
                    : second(x, y, (x & 1) == 0 ? 0 : 1, valueKey);
                appendRow(Kind::EdgeMidpoint, static_cast<int>(edge), component,
                    component == 1 && reflectedEdges[edge], std::move(op));
            }
        }
        for (int c = 0; c < native.cellIndexNum(); ++c)
            appendRow(Kind::CellCenter, c, 0, false,
                selector({Kind::CellCenter, c, 0}));

        const auto vertexRow = [](int v, int component) {
            return static_cast<std::size_t>(v) * 9 + component;
        };
        const std::size_t edgeBase = source.nodes.size() * 9;
        const std::size_t centerBase = edgeBase + source.edgeMidpoints.size() * 3;
        const auto edgeRow = [&](int e, int component) {
            return edgeBase + static_cast<std::size_t>(e) * 3 + component;
        };
        // Derivative orders for [V0,Vmid,V1,D10,D11,D20,D21].
        constexpr std::array<int, 7> location{{0, 1, 2, 0, 2, 0, 2}};
        constexpr std::array<int, 7> order{{0, 0, 0, 1, 1, 2, 2}};
        constexpr int vertexComponent[3][3] = {
            {0, 1, 4}, {2, 3, 6}, {5, 7, 8}
        }; // [y derivative][x derivative]
        cellRows.resize(source.cellCenters.size());
        for (int y = 0; y < ny; ++y)
            for (int x = 0; x < nx; ++x) {
                auto& local = cellRows[cell(x, y)];
                for (int j = 0; j < 7; ++j)
                    for (int i = 0; i < 7; ++i) {
                        const int px = 2 * x + location[i], py = 2 * y + location[j];
                        std::size_t row = 0;
                        if ((px & 1) == 0 && (py & 1) == 0)
                            row = vertexRow(vertex(px / 2, py / 2),
                                vertexComponent[order[j]][order[i]]);
                        else if ((px & 1) == 0)
                            row = edgeRow(vertical[static_cast<std::size_t>(py / 2)
                                * (nx + 1) + px / 2], order[i]);
                        else if ((py & 1) == 0)
                            row = edgeRow(horizontal[static_cast<std::size_t>(py / 2)
                                * nx + px / 2], order[j]);
                        else row = centerBase + static_cast<std::size_t>(cell(x, y));
                        local[static_cast<std::size_t>(7 * j + i)] = row;
                    }
            }
        delegatedValidation("final mesh binding", [&] { mesh.validate(native); });
        Detail::ExactExpansion::validateEnvironment();
    }

    void validateSymmetry(const ActiveQuarticHermitePatch& patch) const
    {
        for (std::size_t v = 0; v < reflectedVertices.size(); ++v) {
            const auto reflection = reflectedVertices[v];
            const auto& node = patch.nodes[v];
            require(!reflection[0] || node.dx == 0.0,
                "scalar x-symmetry requires an exactly zero vertex normal datum");
            require(!reflection[1] || node.dy == 0.0,
                "scalar y-symmetry requires an exactly zero vertex normal datum");
            require(!(reflection[0] || reflection[1]) || node.dxy == 0.0,
                "scalar symmetry requires an exactly zero vertex mixed datum");
        }
        for (std::size_t e = 0; e < reflectedEdges.size(); ++e)
            require(!reflectedEdges[e] || patch.edgeMidpoints[e].normalDerivative == 0.0,
                "scalar symmetry requires an exactly zero edge normal datum");
    }

    void appendRow(Kind kind, int index, int component,
        bool constrained, OperatorRow op)
    {
        Row row;
        row.kind = kind;
        row.index = index;
        row.component = component;
        row.symmetryConstrained = constrained;
        require(std::isfinite(op.denominator) && op.denominator > 0.0L,
            "invalid retained row scale");
        std::map<SourceKey, std::int64_t> numerators;
        for (const auto& term : op.terms)
            numerators[term.source] += term.numerator;
        std::int64_t sum = 0;
        for (const auto& entry : numerators) {
            sum += entry.second;
            if (entry.second == 0) continue;
            const long double coefficient =
                static_cast<long double>(entry.second) / op.denominator;
            require(std::isfinite(coefficient) && coefficient != 0.0L,
                "physical source coefficient is not finitely representable");
            row.sourceWeights.push_back({entry.first.kind, entry.first.index,
                entry.first.component, coefficient});
        }
        require(!op.derived || sum == 0,
            "derived integer row does not annihilate a constant source");
        row.value = evaluate(op, source);
        rows.push_back(std::move(row));
        operators.push_back(std::move(op));
    }
};

ActiveQuarticHermiteCompletion::ActiveQuarticHermiteCompletion(
    std::shared_ptr<const Data> data) : data_(std::move(data)) {}

const ActiveQuarticHermiteCompletion::Data& ActiveQuarticHermiteCompletion::data() const
{
    if (!data_)
        throw std::runtime_error("[ActiveQuarticHermiteCompletion] unbound moved-from owner");
    return *data_;
}

ActiveQuarticHermiteCompletion ActiveQuarticHermiteCompletion::fromUniformPatch(
    const QuadtreeMesh& mesh, const ActiveQuarticHermitePatch& patch,
    const std::vector<SymmetryPlane>& symmetryPlanes)
{
    return ActiveQuarticHermiteCompletion(
        std::make_shared<const Data>(mesh, patch, symmetryPlanes));
}

const ActiveQuarticHermitePatch& ActiveQuarticHermiteCompletion::sourcePatch() const
{
    return data().source;
}

const std::vector<ActiveQuarticHermiteCompletion::Row>&
ActiveQuarticHermiteCompletion::rows() const
{
    return data().rows;
}

const std::array<std::size_t, 49>&
ActiveQuarticHermiteCompletion::cellRows(int nativeCellIndex) const
{
    const auto& retained = data();
    require(nativeCellIndex >= 0
        && static_cast<std::size_t>(nativeCellIndex) < retained.cellRows.size(),
        "foreign native cell row query");
    return retained.cellRows[static_cast<std::size_t>(nativeCellIndex)];
}

std::array<double, 49>
ActiveQuarticHermiteCompletion::physicalCellDofs(int nativeCellIndex) const
{
    const auto& retained = data();
    const auto& selected = cellRows(nativeCellIndex);
    std::array<double, 49> result{};
    for (std::size_t i = 0; i < result.size(); ++i)
        result[i] = retained.rows[selected[i]].value;
    return result;
}

const std::vector<ActiveQuarticHermiteCompletion::SymmetryPlane>&
ActiveQuarticHermiteCompletion::symmetryPlanes() const
{
    return data().symmetry;
}

std::vector<double> ActiveQuarticHermiteCompletion::applyRates(
    const ActiveQuarticHermitePatch& formalRates) const
{
    Detail::ExactExpansion::validateEnvironment();
    const auto& retained = data();
    validateLayout(formalRates, retained.source);
    retained.validateSymmetry(formalRates);
    std::vector<double> result;
    result.reserve(retained.operators.size());
    for (const auto& row : retained.operators)
        result.push_back(evaluate(row, formalRates));
    Detail::ExactExpansion::validateEnvironment();
    return result;
}

std::vector<ActiveQuarticHermiteCompletion::SourceWeight>
ActiveQuarticHermiteCompletion::pullbackCell(int nativeCellIndex,
    const std::array<long double, 49>& completedWeights) const
{
    Detail::ExactExpansion::validateEnvironment();
    const auto& retained = data();
    const auto& selected = cellRows(nativeCellIndex);
    std::vector<std::pair<std::size_t, long double>> coefficients;
    coefficients.reserve(selected.size());
    for (std::size_t i = 0; i < selected.size(); ++i)
        coefficients.emplace_back(selected[i], completedWeights[i]);
    return mergedPullback(retained.rows, coefficients);
}

std::vector<ActiveQuarticHermiteCompletion::SourceWeight>
ActiveQuarticHermiteCompletion::pullbackRows(
    const std::vector<long double>& completedWeights) const
{
    Detail::ExactExpansion::validateEnvironment();
    const auto& retained = data();
    require(completedWeights.size() == retained.rows.size(),
        "completed-row pullback has a foreign layout");
    std::vector<std::pair<std::size_t, long double>> coefficients;
    coefficients.reserve(completedWeights.size());
    for (std::size_t i = 0; i < completedWeights.size(); ++i)
        coefficients.emplace_back(i, completedWeights[i]);
    return mergedPullback(retained.rows, coefficients);
}

void ActiveQuarticHermiteCompletion::validateFor(
    const QuadtreeMesh& mesh, const ActiveQuarticHermitePatch& patch,
    const std::vector<SymmetryPlane>& symmetryPlanes) const
{
    const auto& retained = data();
    delegatedValidation("mesh binding", [&] { retained.mesh.validate(mesh); });
    validateLayout(patch, retained.source);
    require(symmetryPlanes == retained.symmetry,
        "scalar symmetry declarations differ from retained owner");
    for (std::size_t v = 0; v < patch.nodes.size(); ++v) {
        const auto& actual = patch.nodes[v];
        const auto& expected = retained.source.nodes[v];
        require(actual.value == expected.value && actual.dx == expected.dx
            && actual.dy == expected.dy && actual.dxy == expected.dxy,
            "formal vertex values differ from retained source");
    }
    for (std::size_t e = 0; e < patch.edgeMidpoints.size(); ++e) {
        const auto& actual = patch.edgeMidpoints[e];
        const auto& expected = retained.source.edgeMidpoints[e];
        require(actual.value == expected.value
            && actual.normalDerivative == expected.normalDerivative,
            "formal edge values differ from retained source");
    }
    require(patch.cellCenters == retained.source.cellCenters,
        "actual cell-center values differ from retained source");
    retained.validateSymmetry(patch);
}
} // namespace LAMP::Geometry
