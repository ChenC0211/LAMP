#include "Core/Geometry/LevelSetInterfacePanels.h"
#include "Core/Geometry/Detail/LevelSetInterfaceChainOrder.h"
#include "Core/Geometry/Detail/LevelSetLocalProjectionSearch.h"
#include "Core/Geometry/Detail/LevelSetInterfaceQuadrature.h"
#include "Core/Geometry/Detail/ExactExpansion.h"
#include "Core/Mesh/CartesianMeasure.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/UniformLeafGrid.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace
{
using Snapshot = LevelSetHermiteGeometry::Snapshot;
using Frame = LevelSetHermiteGeometry::Frame;
using Panel = LevelSetInterfacePanels::Panel;
using PanelGraph = LevelSetInterfacePanels::PanelGraph;
using EndpointKind = LevelSetInterfacePanels::EndpointKind;

void require(bool ok, const std::string& reason)
{
    if (!ok) throw std::runtime_error("[LevelSetInterfacePanels] " + reason);
}
bool finite(Vector2D p) { return std::isfinite(p.x) && std::isfinite(p.y); }
bool same(Vector2D a, Vector2D b) { return a.x == b.x && a.y == b.y; }
using LAMP::Geometry::Detail::InterfaceQuadrature::checked;
using LAMP::Geometry::Detail::InterfaceQuadrature::offAxisWeight;
Vector2D midpoint(Vector2D a, Vector2D b)
{
    return {checked((static_cast<long double>(a.x) + b.x) / 2.0L, "midpoint overflow"),
        checked((static_cast<long double>(a.y) + b.y) / 2.0L, "midpoint overflow")};
}
double dot(Vector2D a, Vector2D b) { return a.x * b.x + a.y * b.y; }
Vector2D tangent(const Frame& f) { return {-f.normal.y, f.normal.x}; }

struct PhysicalEdge
{
    int label = -1, index = -1, i = -1, j = -1;
    bool horizontal = false, cut = false;
    Vector2D begin, end, crossing;
    double crossingParameter = std::numeric_limits<double>::quiet_NaN();
    std::vector<int> cells; // Incident cut-cell Cartesian indices.
};

struct Grid
{
    int nx, ny;
    Vector2D lower, upper, h;
    std::vector<int> cellLabels, vertexLabels, xEdges, yEdges;
    std::vector<PhysicalEdge> edges;
    std::vector<std::uint8_t> active;

    bool cellActive(int cell) const { return cell >= 0 && active[static_cast<std::size_t>(cell)] != 0; }
    std::array<int, 2> adjacent(const PhysicalEdge& edge) const
    {
        return {{edge.horizontal ? (edge.j > 0 ? edge.i + nx * (edge.j - 1) : -1)
                : (edge.i > 0 ? edge.i - 1 + nx * edge.j : -1),
            edge.horizontal ? (edge.j < ny ? edge.i + nx * edge.j : -1)
                : (edge.i < nx ? edge.i + nx * edge.j : -1)}};
    }
    bool activeBoundary(Vector2D point) const
    {
        for (const auto& edge : edges)
        {
            const auto cells = adjacent(edge);
            if (cellActive(cells[0]) == cellActive(cells[1])) continue;
            if (edge.horizontal ? point.y == edge.begin.y && point.x >= edge.begin.x && point.x <= edge.end.x
                : point.x == edge.begin.x && point.y >= edge.begin.y && point.y <= edge.end.y) return true;
        }
        return false;
    }

    double coordinate(int index, bool x) const
    {
        const int count = x ? nx : ny;
        if (index == count) return x ? upper.x : upper.y;
        return (x ? lower.x : lower.y) + index * (x ? h.x : h.y);
    }
    Vector2D vertex(int i, int j) const { return {coordinate(i, true), coordinate(j, false)}; }
    int cell(int i, int j) const { return cellLabels[i + nx * j]; }

    Grid(const QuadtreeMesh& mesh, const Snapshot& geometry, CartesianGridScope scope)
        : nx(geometry.nx()), ny(geometry.ny()), lower(geometry.minCoord()),
          upper(geometry.maxCoord()), h(geometry.cellSize())
    {
        const UniformLeafGrid layout(mesh, scope);
        require(nx >= 4 && ny >= 4 && nx == layout.nx() && ny == layout.ny()
            && same(lower, layout.minCoord()) && same(upper, layout.maxCoord()) && same(h, layout.spacing()),
            "mesh and retained snapshot must have the same complete Cartesian layout");
        const long long nxf = (static_cast<long long>(nx) + 1) * ny;
        const long long nyf = static_cast<long long>(nx) * (ny + 1);
        require(nxf + nyf <= std::numeric_limits<int>::max(),
            "native panel edge count exceeds supported index range");
        // Keep the panel geometry's original precision/metric gates in addition
        // to the shared topology validator; accepting leaves is not a relaxed
        // curve-resolution or integration contract.
        const double tolerance = 256.0 * std::numeric_limits<double>::epsilon()
            * std::max({h.x, h.y, std::abs(lower.x), std::abs(lower.y), std::abs(upper.x), std::abs(upper.y)});
        require(std::isfinite(tolerance) && tolerance < std::min(h.x, h.y) / 16.0,
            "Cartesian coordinates are unresolved");
        const auto metric = [](double a, double b)
        { return std::isfinite(a) && a > 0.0 && std::abs(static_cast<long double>(a) - b)
            <= 64.0L * std::numeric_limits<double>::epsilon() * b; };
        cellLabels = layout.cells();
        vertexLabels = layout.vertices();
        active.resize(cellLabels.size());
        for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i)
            active[static_cast<std::size_t>(j) * nx + i] = geometry.cellActive(i, j);
        for (int label : cellLabels)
        {
            const auto& cell = mesh.cells(label);
            require(metric(cell.size().x, h.x) && metric(cell.size().y, h.y), "invalid uniform leaf cell");
        }
        xEdges.assign(static_cast<std::size_t>(nxf), -1);
        yEdges.assign(static_cast<std::size_t>(nyf), -1);
        // Invert only the already-validated maps. Preserve original native-edge
        // traversal/seed order instead of changing it to Cartesian label order.
        std::vector<int> xPositions(static_cast<std::size_t>(nxf));
        std::vector<int> yPositions(static_cast<std::size_t>(nyf));
        for (std::size_t p = 0; p < layout.vertical().size(); ++p)
            xPositions[mesh.vEdges(layout.vertical()[p]).index()] = static_cast<int>(p);
        for (std::size_t p = 0; p < layout.horizontal().size(); ++p)
            yPositions[mesh.hEdges(layout.horizontal()[p]).index()] = static_cast<int>(p);
        const auto inspectEdge = [&](const Edge& edge, int label, int index, bool horizontal)
        {
            require(metric(edge.length(), horizontal ? h.x : h.y), "invalid uniform edge");
            const int position = (horizontal ? yPositions : xPositions)[index];
            const int width = horizontal ? nx : nx + 1;
            const int i = position % width, j = position / width;
            PhysicalEdge record;
            record.label = label; record.index = index; record.i = i; record.j = j; record.horizontal = horizontal;
            record.begin = vertex(i, j); record.end = vertex(i + (horizontal ? 1 : 0), j + (horizontal ? 0 : 1));
            // Snapshot-defined grid coordinates identify shared physical edges;
            // layout validation is not a query-point snapping rule.
            const int id = static_cast<int>(edges.size());
            (horizontal ? yEdges : xEdges)[position] = id;
            edges.push_back(std::move(record));
        };
        edges.reserve(static_cast<std::size_t>(nxf + nyf));
        for (int k = 0; k < mesh.vEdgeIndexNum(); ++k)
        {
            const int label = mesh.vEdgeIndexToLabel[k];
            inspectEdge(mesh.vEdges(label), label, k, false);
        }
        for (int k = 0; k < mesh.hEdgeIndexNum(); ++k)
        {
            const int label = mesh.hEdgeIndexToLabel[k];
            inspectEdge(mesh.hEdges(label), label, k, true);
        }
    }
};

// Only actual sign-changing interior roots are accepted. Exact-zero endpoints
// are unsupported even if the other part of the edge has a definite phase.
// The public partitioner itself rejects unresolved tangencies/small intervals.
std::vector<Vector2D> crossings(const Snapshot& geometry, Vector2D a, Vector2D b)
{
    require(!same(a, b), "zero-length root query");
    require(geometry.value(a) != 0.0 && geometry.value(b) != 0.0,
        "vertex/edge-endpoint zero or along-edge interface is unsupported");
    const auto pieces = geometry.partitionSegment(a, b);
    require(!pieces.empty() && same(pieces.front().begin, a) && same(pieces.back().end, b),
        "incomplete public segment partition");
    std::vector<Vector2D> result;
    for (std::size_t k = 1; k < pieces.size(); ++k)
    {
        require(same(pieces[k - 1].end, pieces[k].begin)
            && pieces[k - 1].negative != pieces[k].negative, "unresolved phase partition boundary");
        const auto p = pieces[k].begin;
        require(finite(p) && !same(p, a) && !same(p, b), "unresolved endpoint crossing");
        result.push_back(p);
    }
    return result;
}

using LAMP::Geometry::Detail::InterfaceQuadrature::gaussX;
using LAMP::Geometry::Detail::InterfaceQuadrature::gaussW;

struct CurveGraph
{
    const Snapshot& geometry;
    Vector2D lower, upper, begin, end;
    bool x;
    double direction;
    bool cellLocal = false;
    int cellI = -1, cellJ = -1;
    Vector2D spacing{1.0, 1.0};

    long double physicalSpan() const
    {
        const long double span = std::abs(static_cast<long double>(x ? end.x : end.y) - (x ? begin.x : begin.y));
        return cellLocal ? span * (x ? spacing.x : spacing.y) : span;
    }

    Frame atCoordinate(double coordinate, bool checkOrientation = true) const
    {
        const double first = x ? begin.x : begin.y, last = x ? end.x : end.y;
        require(std::isfinite(coordinate) && coordinate > std::min(first, last)
            && coordinate < std::max(first, last),
            "graph quadrature coordinate is not representably interior");
        const Vector2D a = x ? Vector2D{coordinate, lower.y} : Vector2D{lower.x, coordinate};
        const Vector2D b = x ? Vector2D{coordinate, upper.y} : Vector2D{upper.x, coordinate};
        const auto roots = cellLocal ? geometry.intersectCellSegment(cellI, cellJ, a, b) : crossings(geometry, a, b);
        require(roots.size() == 1, "cell curve is not a single resolved graph at a quadrature section");
        Frame result = cellLocal ? geometry.jetAtCellLocal(cellI, cellJ, roots.front()) : geometry.jet(roots.front());
        if (cellLocal)
            require(result.hasCellLocalPoint && result.cellI == cellI && result.cellJ == cellJ
                && same(result.cellLocalPoint, roots.front()) && roots.front().x > 0.0 && roots.front().x < 1.0
                && roots.front().y > 0.0 && roots.front().y < 1.0,
                "quadrature root is not strictly inside its cell-local domain");
        else require(result.point.x > lower.x && result.point.x < upper.x
            && result.point.y > lower.y && result.point.y < upper.y, "quadrature root is not strictly inside its cell");
        if (checkOrientation) metric(result);
        return result;
    }

    Frame at(double fraction, bool checkOrientation = true) const
    {
        const char* stage = "graph coordinate";
        try
        {
            require(fraction > 0.0 && fraction < 1.0, "graph quadrature parameter must be interior");
            const double first = x ? begin.x : begin.y, last = x ? end.x : end.y;
            const double coordinate = checked(static_cast<long double>(first)
                + fraction * (static_cast<long double>(last) - first), "graph coordinate overflow");
            stage = "transverse root, local jet and oriented metric";
            return atCoordinate(coordinate, checkOrientation);
        }
        catch (const std::exception& error)
        {
            std::ostringstream context;
            context << std::setprecision(17) << "[LevelSetInterfacePanels::CurveGraph] cell=("
                << cellI << ',' << cellJ << "), stage=" << stage
                << ", coordinates=" << (cellLocal ? "cell-local" : "physical")
                << ", begin=(" << begin.x << ',' << begin.y << "), end=(" << end.x << ',' << end.y
                << "), chosen=" << (x ? 'x' : 'y') << ", fraction=" << fraction << ": " << error.what();
            throw std::runtime_error(context.str());
        }
    }
    double metric(const Frame& f) const
    {
        const double component = x ? -f.normal.y : f.normal.x;
        // No near-zero normal division and no fallback to a chord or another
        // scalar. This is a sampled regular-graph check, not a new interval
        // certificate; the Snapshot already enforces its resolved-cell scope.
        require(std::isfinite(component) && component * direction > 256.0 * std::numeric_limits<double>::epsilon(),
            "chosen cell graph has an unresolved or reversed tangent");
        return checked(1.0L / std::abs(component), "graph arc metric overflow");
    }
    double integrate(double finish, bool physical = false) const
    {
        require(finish > 0.0 && finish <= 1.0, "invalid local arclength interval");
        const long double span = physicalSpan();
        long double sum = 0.0L;
        for (std::size_t k = 0; k < gaussX.size(); ++k)
        {
            const double fraction = checked(static_cast<long double>(finish) * (1.0L + gaussX[k]) / 2.0L,
                "quadrature parameter underflow");
            if (physical)
            {
                const auto frame = at(fraction);
                sum += static_cast<long double>(gaussW[k]) * metric(frame) * offAxisWeight(frame);
            }
            else
                sum += static_cast<long double>(gaussW[k]) * metric(at(fraction));
        }
        const double result = checked(span * finish * sum / 2.0L, "unrepresentable local arclength");
        require(result > 0.0, "nonpositive local arclength");
        return result;
    }
};

// Every explicitly sector-enabled declared-chain panel supplies local endpoints.
// The edge-root parameters remain authoritative even when a strictly interior
// quadrature point rounds to the same global double as a cell boundary.
CurveGraph curveGraph(const Grid& grid, const Snapshot& geometry, int cell,
    Vector2D begin, Vector2D end, bool x, double direction,
    const std::array<Vector2D, 2>* localEndpoints)
{
    const int i = cell % grid.nx, j = cell / grid.nx;
    if (!localEndpoints)
    {
        CurveGraph result{geometry, grid.vertex(i,j), grid.vertex(i+1,j+1), begin, end, x, direction};
        result.cellI = i; result.cellJ = j;
        return result;
    }
    for (const auto point : *localEndpoints)
        require(finite(point) && point.x >= 0.0 && point.x <= 1.0 && point.y >= 0.0 && point.y <= 1.0,
            "invalid retained cell-local panel endpoint");
    CurveGraph result{geometry, {0.0,0.0}, {1.0,1.0}, (*localEndpoints)[0], (*localEndpoints)[1], x, direction};
    result.cellLocal = true; result.cellI = i; result.cellJ = j; result.spacing = grid.h;
    return result;
}

bool samePanelSide(const Grid& grid, int cartesianCell, Vector2D begin, Vector2D end)
{
    const int i = cartesianCell % grid.nx, j = cartesianCell / grid.nx;
    const auto lower = grid.vertex(i, j), upper = grid.vertex(i + 1, j + 1);
    return (begin.x == end.x && (begin.x == lower.x || begin.x == upper.x))
        || (begin.y == end.y && (begin.y == lower.y || begin.y == upper.y));
}

Frame graphMidpointFrame(const Grid& grid, const Snapshot& geometry,
    int cartesianCell, Vector2D begin, Vector2D end,
    const std::array<Vector2D, 2>* localEndpoints = nullptr)
{
    const auto difference = localEndpoints ? (*localEndpoints)[1]-(*localEndpoints)[0] : end-begin;
    const bool x = localEndpoints ? std::abs(difference.x) >= std::abs(difference.y)
        : std::abs(difference.x/grid.h.x) >= std::abs(difference.y/grid.h.y);
    const double delta = x ? difference.x : difference.y;
    require(finite(difference) && std::isfinite(delta) && delta != 0.0, "unresolved graph midpoint span");
    const CurveGraph graph = curveGraph(grid, geometry, cartesianCell, begin, end, x,
        delta > 0.0 ? 1.0 : -1.0, localEndpoints);
    // A same-side chord need not have a stationary foot on this curve branch.
    // A short cross-edge arc can also lose its transverse interior coordinate
    // when a physical midpoint is rounded. Resolve the original scalar at the
    // graph-coordinate midpoint, retaining local coordinates when supplied.
    // Direction comes from this actual interior point; makePanel subsequently
    // checks every original oriented metric.
    auto frame = graph.at(0.5, false);
    if (localEndpoints)
    {
        const auto query = midpoint((*localEndpoints)[0], (*localEndpoints)[1]);
        frame.stationarityResidual = checked(std::abs(
            (static_cast<long double>(frame.cellLocalPoint.x)-query.x)*grid.h.x*(-frame.normal.y)
            + (static_cast<long double>(frame.cellLocalPoint.y)-query.y)*grid.h.y*frame.normal.x),
            "unrepresentable local graph midpoint stationarity residual");
    }
    else
    {
        const auto query = midpoint(begin, end);
        frame.stationarityResidual = checked(std::abs(
            (static_cast<long double>(frame.point.x)-query.x)*(-frame.normal.y)
            + (static_cast<long double>(frame.point.y)-query.y)*frame.normal.x),
            "unrepresentable graph midpoint stationarity residual");
    }
    return frame;
}

Panel makePanel(const Grid& grid, const Snapshot& geometry, int cartesianCell,
    int solverCell, Vector2D begin, Vector2D end, bool axisymmetric, bool axisEndpoints,
    bool declared = false, const std::array<Vector2D, 2>* localEndpoints = nullptr)
{
    const int i = cartesianCell % grid.nx, j = cartesianCell / grid.nx;
    Panel panel;
    panel.cellIndex = solverCell; panel.begin = begin; panel.end = end;
    if (localEndpoints)
    {
        require(declared && geometry.hasCornerSectors(),
            "cell-local panel path requires an explicit sector declared-chain panel");
        panel.beginFrame = geometry.jetAtCellLocal(i, j, (*localEndpoints)[0]);
        panel.endFrame = geometry.jetAtCellLocal(i, j, (*localEndpoints)[1]);
    }
    else { panel.beginFrame = geometry.jet(begin); panel.endFrame = geometry.jet(end); }
    panel.frame = localEndpoints || (declared && samePanelSide(grid, cartesianCell, begin, end))
        ? graphMidpointFrame(grid, geometry, cartesianCell, begin, end, localEndpoints)
        : geometry.projectedFrame(midpoint(begin, end));
    if (axisymmetric)
    {
        const auto endpoint = [&](const Frame& frame)
        {
            if (axisEndpoints && frame.point.x == 0.0)
                require(frame.normal.x == 0.0, "axis endpoint is not radially symmetric");
            else offAxisWeight(frame);
        };
        endpoint(panel.beginFrame);
        endpoint(panel.endFrame);
        offAxisWeight(panel.frame);
    }
    const auto lower = grid.vertex(i, j), upper = grid.vertex(i + 1, j + 1);
    if (localEndpoints)
        require(panel.frame.hasCellLocalPoint && panel.frame.cellI == i && panel.frame.cellJ == j
            && panel.frame.cellLocalPoint.x > 0.0 && panel.frame.cellLocalPoint.x < 1.0
            && panel.frame.cellLocalPoint.y > 0.0 && panel.frame.cellLocalPoint.y < 1.0,
            "cell-local midpoint did not stay in the panel's cell");
    else require(panel.frame.point.x > lower.x && panel.frame.point.x < upper.x
        && panel.frame.point.y > lower.y && panel.frame.point.y < upper.y,
        "chord-midpoint projection did not stay in the panel's cell");
    const Vector2D difference = end - begin;
    const Vector2D localDifference = localEndpoints ? (*localEndpoints)[1]-(*localEndpoints)[0] : Vector2D{};
    const double orientation = localEndpoints ? checked(
        static_cast<long double>(localDifference.x)*grid.h.x*(-panel.frame.normal.y)
        + static_cast<long double>(localDifference.y)*grid.h.y*panel.frame.normal.x,
        "unrepresentable local panel orientation") : dot(difference, tangent(panel.frame));
    require(finite(difference) && std::isfinite(orientation) && orientation > 0.0,
        "ordered panel disagrees with the negative-to-positive tangent orientation");
    // Prefer the larger normalized endpoint span; do not select an axis from
    // an analytic reference, curve fit, or shape-specific rule.
    panel.coordinateIsX = localEndpoints ? std::abs(localDifference.x) >= std::abs(localDifference.y)
        : std::abs(difference.x / grid.h.x) >= std::abs(difference.y / grid.h.y);
    const double delta = localEndpoints ? (panel.coordinateIsX ? localDifference.x : localDifference.y)
        : (panel.coordinateIsX ? difference.x : difference.y);
    require(delta != 0.0 && std::isfinite(delta), "zero graph-coordinate span");
    const CurveGraph graph = curveGraph(grid, geometry, cartesianCell, begin, end, panel.coordinateIsX,
        delta > 0.0 ? 1.0 : -1.0, localEndpoints);
    graph.metric(panel.beginFrame); graph.metric(panel.endFrame); graph.metric(panel.frame);
    const long double span = localEndpoints ? graph.physicalSpan() : std::abs(static_cast<long double>(delta));
    long double measure = 0.0L;
    long double physicalMeasure = 0.0L;
    for (std::size_t k = 0; k < gaussX.size(); ++k)
    {
        const double fraction = (1.0 + gaussX[k]) / 2.0;
        LevelSetInterfacePanels::QuadratureNode node;
        node.frame = graph.at(fraction);
        node.weight = checked(span * gaussW[k] * graph.metric(node.frame) / 2.0L, "unrepresentable arc quadrature weight");
        require(node.weight > 0.0, "nonpositive arc quadrature weight");
        node.localArclength = graph.integrate(fraction);
        require(panel.quadratureNodes.empty() || node.localArclength > panel.quadratureNodes.back().localArclength,
            "quadrature arclength order is unresolved");
        node.physicalWeight = axisymmetric
            ? checked(static_cast<long double>(offAxisWeight(node.frame)) * node.weight,
                "unrepresentable physical quadrature weight") : node.weight;
        node.localPhysicalMeasure = axisymmetric ? graph.integrate(fraction, true) : node.localArclength;
        require(node.physicalWeight > 0.0 && (panel.quadratureNodes.empty()
            || node.localPhysicalMeasure > panel.quadratureNodes.back().localPhysicalMeasure),
            "physical quadrature weight/order is unresolved");
        measure += node.weight;
        if (axisymmetric) physicalMeasure += node.physicalWeight;
        panel.quadratureNodes.push_back(std::move(node));
    }
    panel.measure = checked(measure, "panel measure overflow");
    require(panel.measure > panel.quadratureNodes.back().localArclength, "panel endpoint arclength is not ordered");
    panel.physicalMeasure = axisymmetric ? checked(physicalMeasure, "panel physical measure overflow") : panel.measure;
    require(panel.physicalMeasure > panel.quadratureNodes.back().localPhysicalMeasure,
        "panel endpoint physical measure is not ordered");
    return panel;
}

PanelGraph retainedGraph(const Grid& grid, int cell, const Panel& panel,
    Vector2D localBegin, Vector2D localEnd)
{
    for (const auto point : {localBegin, localEnd})
        require(finite(point) && point.x >= 0.0 && point.x <= 1.0 && point.y >= 0.0 && point.y <= 1.0,
            "invalid retained cell-local graph endpoint");
    const double begin = panel.coordinateIsX ? localBegin.x : localBegin.y;
    const double end = panel.coordinateIsX ? localEnd.x : localEnd.y;
    require(begin != end, "retained graph has no resolved local coordinate span");
    return {cell % grid.nx, cell / grid.nx, localBegin, localEnd, grid.h, panel.coordinateIsX};
}

struct DeclaredChain
{
    std::vector<Panel> panels;
    std::vector<PanelGraph> graphs;
    std::vector<int> mapping;
    EndpointKind beginKind = EndpointKind::None, endKind = EndpointKind::None;
};

DeclaredChain declaredChain(const QuadtreeMesh& mesh, const Grid& grid, const Snapshot& snapshot,
    const std::vector<int>& signs, bool axisymmetric, bool axisContact,
    const std::vector<Vector2D>& contacts, bool retainGraphs = false)
{
    require(!contacts.empty() && contacts.size() <= 2, "an open chain accepts one or two declared boundary endpoints");
    for (std::size_t k = 0; k < contacts.size(); ++k)
    {
        const auto p = contacts[k];
        require(finite(p) && p.x >= grid.lower.x && p.x <= grid.upper.x
            && p.y >= grid.lower.y && p.y <= grid.upper.y, "declared endpoint lies outside the physical domain");
        const int walls = static_cast<int>(p.x == grid.lower.x) + static_cast<int>(p.x == grid.upper.x)
            + static_cast<int>(p.y == grid.lower.y) + static_cast<int>(p.y == grid.upper.y);
        if (snapshot.restrictedDomain())
            require(grid.activeBoundary(p), "declared endpoint is not incident to an actual active-domain boundary edge");
        else require(walls == 1, "declared endpoints must lie on one physical wall, not a domain corner or in the interior");
        require(!axisymmetric || p.x > 0.0, "r=0 endpoints use radial symmetry, not non-axis boundary declarations");
        for (std::size_t j = 0; j < k; ++j) require(!same(p, contacts[j]), "duplicate declared boundary endpoint");
    }
    struct Node
    {
        Vector2D point;
        int vertex = -1;
        int edge = -1;
        double edgeParameter = 0.0;
        EndpointKind kind = EndpointKind::None;
        std::vector<int> cells;
    };
    std::vector<Node> nodes;
    std::vector<int> vertexNodes(grid.vertexLabels.size(), -1);
    std::vector<std::vector<int>> edgeNodes(grid.edges.size());
    std::vector<bool> matched(contacts.size(), false);
    const bool sectors = snapshot.hasCornerSectors();
    const auto endpointKind = [&](Vector2D p) {
        const bool wall = snapshot.restrictedDomain() ? grid.activeBoundary(p)
            : p.x == grid.lower.x || p.x == grid.upper.x || p.y == grid.lower.y || p.y == grid.upper.y;
        if (!wall) return EndpointKind::None;
        if (axisContact && p.x == 0.0 && p.y > grid.lower.y && p.y < grid.upper.y) return EndpointKind::Axis;
        for (std::size_t k = 0; k < contacts.size(); ++k)
            if (same(p, contacts[k])) { matched[k] = true; return EndpointKind::DeclaredBoundary; }
        throw std::runtime_error("[LevelSetInterfacePanels] interface touches an undeclared physical boundary point");
    };
    for (std::size_t e = 0; e < grid.edges.size(); ++e)
    {
        const auto& edge = grid.edges[e];
        const int low = edge.horizontal ? (edge.j > 0 ? edge.i + grid.nx * (edge.j - 1) : -1)
            : (edge.i > 0 ? edge.i - 1 + grid.nx * edge.j : -1);
        const int high = edge.horizontal ? (edge.j < grid.ny ? edge.i + grid.nx * edge.j : -1)
            : (edge.i < grid.nx ? edge.i + grid.nx * edge.j : -1);
        if (!grid.cellActive(low) && !grid.cellActive(high)) continue;
        if ((grid.cellActive(low) && signs[low] != 0) || (grid.cellActive(high) && signs[high] != 0)) continue;
        // The public root query retains endpoint/tangent/coincident information
        // which phase partitions intentionally omit. Every physical edge is
        // queried once, with exact grid-vertex identity shared by all edges.
        const auto roots = snapshot.intersectGridEdge(edge.begin, edge.end);
        require(!roots.coincident, "an interface coincident with a physical grid edge is unsupported");
        if (sectors) require(roots.parameters.size() == roots.points.size(),
            "sector edge roots require their original local parameters");
        else if (retainGraphs) require(roots.parameters.size() == roots.points.size(),
            "retained edge roots require their original local parameters");
        for (std::size_t root = 0; root < roots.points.size(); ++root)
        {
            const auto point = roots.points[root];
            require(finite(point), "non-finite grid-edge intersection");
            int vertex = -1;
            if (same(point, edge.begin)) vertex = edge.i + (grid.nx + 1) * edge.j;
            else if (same(point, edge.end)) vertex = edge.i + (edge.horizontal ? 1 : 0)
                + (grid.nx + 1) * (edge.j + (edge.horizontal ? 0 : 1));
            else
            {
                const auto frame = snapshot.jet(point);
                const double derivative = edge.horizontal ? frame.normal.x : frame.normal.y;
                require(std::abs(derivative) > 256.0 * std::numeric_limits<double>::epsilon(),
                    "tangent or unresolved interior grid-edge root is unsupported");
            }
            if (vertex >= 0 && vertexNodes[vertex] >= 0)
            {
                edgeNodes[e].push_back(vertexNodes[vertex]);
                require(same(nodes[vertexNodes[vertex]].point, point), "inconsistent shared vertex intersection");
                continue;
            }
            Node node;
            node.point = point; node.vertex = vertex; node.kind = endpointKind(point);
            if ((sectors || retainGraphs) && vertex < 0)
            {
                node.edge = static_cast<int>(e); node.edgeParameter = roots.parameters[root];
                require(std::isfinite(node.edgeParameter) && node.edgeParameter > 0.0 && node.edgeParameter < 1.0,
                    "interior sector root lost its strict local edge parameter");
            }
            const int id = static_cast<int>(nodes.size());
            edgeNodes[e].push_back(id);
            if (vertex >= 0) vertexNodes[vertex] = id;
            nodes.push_back(std::move(node));
        }
    }
    require(std::all_of(matched.begin(), matched.end(), [](bool value) { return value; }),
        "a declared boundary endpoint does not exactly match a retained grid-edge root");
    std::vector<std::array<int, 2>> connections(grid.cellLabels.size(), {{-1, -1}});
    std::size_t cutCells = 0;
    for (int j = 0; j < grid.ny; ++j) for (int i = 0; i < grid.nx; ++i)
    {
        const int cell = i + grid.nx * j;
        if (!grid.cellActive(cell)) continue;
        const int sides[4] = {grid.xEdges[i + (grid.nx + 1) * j], grid.yEdges[i + grid.nx * j],
            grid.xEdges[i + 1 + (grid.nx + 1) * j], grid.yEdges[i + grid.nx * (j + 1)]};
        const int corners[4] = {i + (grid.nx + 1) * j, i + 1 + (grid.nx + 1) * j,
            i + 1 + (grid.nx + 1) * (j + 1), i + (grid.nx + 1) * (j + 1)};
        std::vector<int> local;
        for (const int edge : sides)
            for (const int node : edgeNodes[edge])
                if (std::find(local.begin(), local.end(), node) == local.end()) local.push_back(node);
        unsigned excludedCorners = 0;
        local.erase(std::remove_if(local.begin(), local.end(), [&](int node) {
            if (nodes[node].vertex < 0) return false;
            int corner = -1;
            for (int k = 0; k < 4; ++k) if (nodes[node].vertex == corners[k]) corner = k;
            require(corner >= 0, "cell root names an unrelated Cartesian vertex");
            const auto derivative = sectors ? snapshot.zeroVertexGradientSignsInCell(i, j, nodes[node].point)
                : snapshot.zeroVertexGradientSigns(nodes[node].point);
            const int xIn = (corner == 0 || corner == 3 ? 1 : -1) * derivative[0];
            const int yIn = (corner == 0 || corner == 1 ? 1 : -1) * derivative[1];
            // Both STRICT inward partials have the same sign. Continuity of
            // this cell's gradient excludes a zero germ in its inward cone;
            // an explicit corner sector need not share another cell's jet.
            // Remove only this cell incidence, never the global root or any
            // other edge root. Zero/tangent partials get no such exemption.
            if (xIn == 0 || xIn != yIn) return false;
            excludedCorners |= 1u << corner;
            return true;
        }), local.end());
        if (local.empty())
        {
            require(signs[cell] != 0 || (excludedCorners != 0
                && snapshot.closedCellSignExceptVertices(mesh.cells(grid.cellLabels[cell]).centerCoord(), excludedCorners) != 0),
                "uncertified no-root cell in declared open-chain geometry");
            continue;
        }
        if (local.size() == 1)
        {
            // A shared zero corner can touch this cell without an interior
            // segment. Skip it ONLY with the authoritative closed-cell sign
            // certificate excluding that one corner, never from root count,
            // a point sample, or the approximate regular-cell partial test.
            require(nodes[local.front()].vertex >= 0, "an isolated nonvertex cell contact is unsupported");
            unsigned mask = 0;
            for (int corner = 0; corner < 4; ++corner)
                if (nodes[local.front()].vertex == corners[corner]) mask = 1u << corner;
            require(mask != 0 && snapshot.closedCellSignExceptVertices(
                mesh.cells(grid.cellLabels[cell]).centerCoord(), mask | excludedCorners) != 0,
                "a zero-corner-only cell lacks a strict sign certificate away from that corner");
            continue;
        }
        require(local.size() == 2 && signs[cell] == 0,
            "branch/multiple segments or inconsistent crossings in a declared-chain cell");
        connections[cell] = {{local[0], local[1]}};
        nodes[local[0]].cells.push_back(cell); nodes[local[1]].cells.push_back(cell);
        ++cutCells;
    }
    require(cutCells >= 3, "one resolved interface with at least three panels is required");
    std::vector<int> ends;
    for (std::size_t k = 0; k < nodes.size(); ++k)
    {
        const auto& node = nodes[k];
        if (node.kind != EndpointKind::None)
        {
            require(node.cells.size() == 1, "declared/axis boundary root is not a unique open-chain endpoint");
            ends.push_back(static_cast<int>(k));
        }
        else require(node.cells.size() == 2, "interior intersection has a branch, isolated contact or incomplete incidence");
    }
    require(ends.size() == 2, "a declared open interface must have exactly two boundary endpoints");
    const auto otherNode = [&](int cell, int node) {
        const auto pair = connections[cell];
        require(pair[0] == node || pair[1] == node, "broken declared-chain cell incidence");
        return pair[0] == node ? pair[1] : pair[0];
    };
    const auto localPoint = [&](int node, int cell) {
        const int i = cell % grid.nx, j = cell / grid.nx;
        const auto& root = nodes[node];
        if (root.vertex >= 0)
            return Vector2D{static_cast<double>(root.vertex % (grid.nx+1) - i),
                static_cast<double>(root.vertex / (grid.nx+1) - j)};
        require(root.edge >= 0 && root.edge < static_cast<int>(grid.edges.size()),
            "missing retained local edge-root identity");
        const auto& edge = grid.edges[root.edge];
        return edge.horizontal ? Vector2D{edge.i-i+root.edgeParameter, static_cast<double>(edge.j-j)}
            : Vector2D{static_cast<double>(edge.i-i), edge.j-j+root.edgeParameter};
    };
    const auto endpointIdentity = [&](int node, int cell) {
        const int i=cell%grid.nx,j=cell/grid.nx;
        const int sides[4]={grid.xEdges[i+(grid.nx+1)*j],grid.yEdges[i+grid.nx*j],
            grid.xEdges[i+1+(grid.nx+1)*j],grid.yEdges[i+grid.nx*(j+1)]};
        for(int side=0;side<4;++side)
        {
            const auto& inventory=edgeNodes[sides[side]];
            const auto found=std::find(inventory.begin(),inventory.end(),node);
            if(found!=inventory.end())
                return PanelGraph::Endpoint{side,static_cast<int>(found-inventory.begin())};
        }
        throw std::runtime_error("[LevelSetInterfacePanels] graph endpoint lost original edge-root identity");
    };
    const auto panelContext = [&](int cell, const char* stage, int beginNode, int endNode,
        const std::array<Vector2D, 2>* localEndpoints) {
        std::ostringstream context;
        context << std::setprecision(17) << "[LevelSetInterfacePanels] cellIndex="
            << mesh.cells(grid.cellLabels[cell]).index() << ", cell=(" << cell % grid.nx << ','
            << cell / grid.nx << "), stage=" << stage << ", physicalBegin=("
            << nodes[beginNode].point.x << ',' << nodes[beginNode].point.y << "), physicalEnd=("
            << nodes[endNode].point.x << ',' << nodes[endNode].point.y << ')';
        if (localEndpoints)
            context << ", localBegin=(" << (*localEndpoints)[0].x << ',' << (*localEndpoints)[0].y
                << "), localEnd=(" << (*localEndpoints)[1].x << ',' << (*localEndpoints)[1].y << ')';
        return context.str();
    };
    const auto outgoingOrientation = [&](int node) {
        const int cell = nodes[node].cells.front(), end = otherNode(cell, node);
        const bool sameSide = samePanelSide(grid, cell, nodes[node].point, nodes[end].point);
        const bool local = sectors;
        std::array<Vector2D, 2> endpoints{};
        bool retainedEndpoints = false;
        const char* stage = "orientation endpoint coordinates";
        try
        {
            if (local)
            {
                endpoints = {{localPoint(node,cell),localPoint(end,cell)}};
                retainedEndpoints = true;
            }
            stage = "orientation frame";
            const auto frame = local || sameSide
                ? graphMidpointFrame(grid, snapshot, cell, nodes[node].point, nodes[end].point, local ? &endpoints : nullptr)
                : snapshot.projectedFrame(midpoint(nodes[node].point, nodes[end].point));
            stage = "endpoint orientation";
            const double orientation = local ? checked(
                (static_cast<long double>(endpoints[1].x)-endpoints[0].x)*grid.h.x*(-frame.normal.y)
                + (static_cast<long double>(endpoints[1].y)-endpoints[0].y)*grid.h.y*frame.normal.x,
                "unrepresentable local endpoint orientation") : dot(nodes[end].point - nodes[node].point, tangent(frame));
            require(std::isfinite(orientation) && orientation != 0.0, "unresolved declared-chain endpoint orientation");
            return orientation;
        }
        catch (const std::exception& error)
        {
            throw std::runtime_error(panelContext(cell, stage, node, end, retainedEndpoints ? &endpoints : nullptr)
                + ": " + error.what());
        }
    };
    const int start = outgoingOrientation(ends.front()) > 0.0 ? ends.front() : ends.back();
    require(outgoingOrientation(start) > 0.0, "boundary endpoints do not orient one consistent open chain");
    DeclaredChain result;
    result.panels.reserve(cutCells);
    if (retainGraphs) result.graphs.reserve(cutCells);
    result.mapping.assign(static_cast<std::size_t>(mesh.cellIndexNum()), -1);
    std::vector<std::vector<int>> nodeCells;
    for(const auto& node:nodes) nodeCells.push_back(node.cells);
    const auto visits=LAMP::Geometry::Detail::orderInterfaceChain(
        connections,nodeCells,start,nodes[start].cells.front(),cutCells,false);
    int currentNode = start, currentCell = nodes[start].cells.front();
    for (std::size_t p = 0; p < cutCells; ++p)
    {
        currentCell=visits[p].cell; currentNode=visits[p].begin;
        const int solverCell = mesh.cells(grid.cellLabels[currentCell]).index();
        require(result.mapping[solverCell] == -1, "declared chain visits a cell twice or closes early");
        const int end = otherNode(currentCell, currentNode);
        const bool local = sectors;
        std::array<Vector2D, 2> endpoints{};
        bool retainedEndpoints = false;
        const char* stage = "panel endpoint coordinates";
        try
        {
            if (local)
            {
                endpoints = {{localPoint(currentNode,currentCell),localPoint(end,currentCell)}};
                retainedEndpoints = true;
            }
            stage = "panel frames and GL8 quadrature";
            result.panels.push_back(makePanel(grid, snapshot, currentCell, solverCell,
                nodes[currentNode].point, nodes[end].point, axisymmetric, true, true, local ? &endpoints : nullptr));
            if (retainGraphs)
            {
                stage = "retained panel graph";
                result.graphs.push_back(retainedGraph(grid, currentCell, result.panels.back(),
                    localPoint(currentNode, currentCell), localPoint(end, currentCell)));
                result.graphs.back().beginRoot=endpointIdentity(currentNode,currentCell);
                result.graphs.back().endRoot=endpointIdentity(end,currentCell);
            }
        }
        catch (const std::exception& error)
        {
            throw std::runtime_error(panelContext(currentCell, stage, currentNode, end,
                retainedEndpoints ? &endpoints : nullptr) + ": " + error.what());
        }
        result.mapping[solverCell] = static_cast<int>(p);
        currentNode = end;
    }
    require(result.panels.size() == cutCells && currentNode != start
        && nodes[currentNode].kind != EndpointKind::None
        && same(result.panels.front().begin, nodes[start].point)
        && same(result.panels.back().end, nodes[currentNode].point),
        "interface graph did not form one complete declared open chain");
    result.beginKind = nodes[start].kind; result.endKind = nodes[currentNode].kind;
    return result;
}
}

LevelSetInterfacePanels::LevelSetInterfacePanels(const QuadtreeMesh& mesh, const Snapshot& snapshot, bool axisymmetric,
    std::vector<Vector2D> boundaryContacts, CartesianGridScope scope)
    : axisymmetric_(axisymmetric), boundaryContacts_(std::move(boundaryContacts))
{
    initialize(mesh, snapshot, scope);
}

LevelSetInterfacePanels::LevelSetInterfacePanels(const QuadtreeMesh& mesh,
    std::shared_ptr<const Snapshot> snapshot, bool axisymmetric,
    std::vector<Vector2D> boundaryContacts, CartesianGridScope scope)
    : axisymmetric_(axisymmetric), boundaryContacts_(std::move(boundaryContacts)), source_(std::move(snapshot))
{
    require(static_cast<bool>(source_), "retained panel construction requires a nonnull immutable snapshot");
    initialize(mesh, *source_, scope);
}

void LevelSetInterfacePanels::initialize(const QuadtreeMesh& mesh, const Snapshot& snapshot, CartesianGridScope scope)
{
    const bool axisymmetric = axisymmetric_;
    Grid grid(mesh, snapshot, scope);
    require(!axisymmetric || grid.lower.x >= 0.0,
        "axisymmetric panels require a nonnegative radial domain");
    const bool axisContact = axisymmetric && grid.lower.x == 0.0 && snapshot.radialAxisSymmetry();
    // Querying the public scalar also checks the supported FP environment
    // before root work; the caller's floating-point controls are never changed.
    std::vector<int> signs(grid.cellLabels.size());
    for (std::size_t k = 0; k < signs.size(); ++k)
        if (grid.active[k]) signs[k] = snapshot.closedCellSign(mesh.cells(grid.cellLabels[k]).centerCoord());
    if (!boundaryContacts_.empty())
    {
        auto result = declaredChain(mesh, grid, snapshot, signs, axisymmetric, axisContact, boundaryContacts_,
            static_cast<bool>(source_));
        panels_ = std::move(result.panels); cellToPanel_ = std::move(result.mapping);
        graphs_ = std::move(result.graphs);
        closed_ = false; beginEndpointKind_ = result.beginKind; endEndpointKind_ = result.endKind;
        return;
    }
    require(!snapshot.restrictedDomain(), "active-domain panels require explicit open-chain endpoint declarations");
    for (auto& edge : grid.edges)
    {
        // A strict CLOSED-cell sign excludes all roots on this physical edge.
        // Uncertified edges are partitioned only here, once, never per owner.
        const int low = edge.horizontal ? (edge.j > 0 ? edge.i + grid.nx * (edge.j - 1) : -1)
            : (edge.i > 0 ? edge.i - 1 + grid.nx * edge.j : -1);
        const int high = edge.horizontal ? (edge.j < grid.ny ? edge.i + grid.nx * edge.j : -1)
            : (edge.i < grid.nx ? edge.i + grid.nx * edge.j : -1);
        if ((low >= 0 && signs[low] != 0) || (high >= 0 && signs[high] != 0)) continue;
        const auto roots = crossings(snapshot, edge.begin, edge.end);
        require(roots.size() <= 1, "multiple crossings on a physical edge are unsupported");
        if (roots.empty()) continue;
        require((low >= 0 && high >= 0) || (axisContact && !edge.horizontal && edge.i == 0
            && low == -1 && high >= 0 && roots.front().x == 0.0),
            "interface contacts an unsupported physical domain boundary");
        edge.cut = true; edge.crossing = roots.front();
        if (source_)
        {
            // The result-only path above keeps its original partition/root
            // arithmetic. Retention obtains the original local root parameter
            // from the same public scalar and verifies exact representative
            // identity; it never inverts the rounded global crossing point.
            const auto retained = snapshot.intersectGridEdge(edge.begin, edge.end);
            require(!retained.coincident && retained.points.size() == 1 && retained.parameters.size() == 1
                && same(retained.points.front(), edge.crossing),
                "retained local root disagrees with the original panel edge inventory");
            edge.crossingParameter = retained.parameters.front();
            require(std::isfinite(edge.crossingParameter) && edge.crossingParameter > 0.0
                && edge.crossingParameter < 1.0, "retained graph edge parameter is not a strict interior root");
        }
    }
    std::vector<std::array<int, 2>> connections(grid.cellLabels.size(), {{-1, -1}});
    std::size_t cutCells = 0;
    for (int j = 0; j < grid.ny; ++j) for (int i = 0; i < grid.nx; ++i)
    {
        const int cell = i + grid.nx * j;
        const int sides[4] = {grid.xEdges[i + (grid.nx + 1) * j], grid.yEdges[i + grid.nx * j],
            grid.xEdges[i + 1 + (grid.nx + 1) * j], grid.yEdges[i + grid.nx * (j + 1)]};
        int count = 0;
        for (int edge : sides) if (grid.edges[edge].cut)
        {
            require(count < 2, "branch/multiple interface segments in one cell are unsupported");
            connections[cell][count++] = edge;
            grid.edges[edge].cells.push_back(cell);
        }
        if (count == 0)
        {
            require(signs[cell] != 0, "uncertified no-crossing cell: hidden/near-critical zero topology is unsupported");
            continue;
        }
        require(count == 2 && signs[cell] == 0, "cut cell does not have exactly two consistent crossings");
        ++cutCells;
    }
    require(cutCells >= 3, "one resolved interface with at least three panels is required");
    int startEdge = -1;
    std::vector<int> axisEnds;
    for (std::size_t e = 0; e < grid.edges.size(); ++e) if (grid.edges[e].cut)
    {
        const auto& edge = grid.edges[e];
        if (edge.cells.size() == 1)
        {
            require(axisContact && !edge.horizontal && edge.i == 0 && edge.crossing.x == 0.0,
                "open interface endpoint is not on the symmetry axis");
            axisEnds.push_back(static_cast<int>(e));
        }
        else require(edge.cells.size() == 2, "interface intersection graph is not degree two");
        if (startEdge < 0) startEdge = static_cast<int>(e);
    }
    require(axisEnds.empty() || axisEnds.size() == 2, "axisymmetric open interface requires exactly two axis endpoints");
    closed_ = axisEnds.empty();
    if (!closed_) startEdge = axisEnds.front();
    require(startEdge >= 0, "interface graph is empty");
    int firstCell = grid.edges[startEdge].cells.front();
    const auto otherEdge = [&](int cell, int edge)
    {
        const auto pair = connections[cell];
        require(pair[0] == edge || pair[1] == edge, "broken ordered graph incidence");
        return pair[0] == edge ? pair[1] : pair[0];
    };
    const int firstEnd = otherEdge(firstCell, startEdge);
    const auto firstMid = snapshot.projectedFrame(midpoint(grid.edges[startEdge].crossing, grid.edges[firstEnd].crossing));
    const double orientation = dot(grid.edges[firstEnd].crossing - grid.edges[startEdge].crossing, tangent(firstMid));
    require(std::isfinite(orientation) && orientation != 0.0, "unresolved initial interface orientation");
    if (orientation < 0.0)
    {
        if (closed_) firstCell = grid.edges[startEdge].cells.back();
        else
        {
            startEdge = axisEnds.back();
            firstCell = grid.edges[startEdge].cells.front();
            const int end = otherEdge(firstCell, startEdge);
            const auto mid = snapshot.projectedFrame(midpoint(grid.edges[startEdge].crossing, grid.edges[end].crossing));
            require(dot(grid.edges[end].crossing-grid.edges[startEdge].crossing, tangent(mid)) > 0.0,
                "axis endpoints do not define a consistently oriented open chain");
        }
    }
    std::vector<Panel> panels;
    panels.reserve(cutCells);
    if (source_) graphs_.reserve(cutCells);
    std::vector<int> mapping(static_cast<std::size_t>(mesh.cellIndexNum()), -1);
    int currentCell = firstCell, currentEdge = startEdge;
    for (std::size_t p = 0; p < cutCells; ++p)
    {
        const int solverCell = mesh.cells(grid.cellLabels[currentCell]).index();
        require(mapping[solverCell] == -1, "interface graph closes early or visits a cell twice");
        const int endEdge = otherEdge(currentCell, currentEdge);
        require(endEdge >= 0, "missing panel end intersection");
        try
        {
            panels.push_back(makePanel(grid, snapshot, currentCell, solverCell,
                grid.edges[currentEdge].crossing, grid.edges[endEdge].crossing, axisymmetric, !closed_));
            if (source_)
            {
                const auto localRoot = [&](int edgeId) {
                    const auto& edge = grid.edges[edgeId];
                    const int i = currentCell % grid.nx, j = currentCell / grid.nx;
                    return edge.horizontal
                        ? Vector2D{edge.i-i+edge.crossingParameter, static_cast<double>(edge.j-j)}
                        : Vector2D{static_cast<double>(edge.i-i), edge.j-j+edge.crossingParameter};
                };
                graphs_.push_back(retainedGraph(grid, currentCell, panels.back(),
                    localRoot(currentEdge), localRoot(endEdge)));
                const auto endpointIdentity = [&](int edgeId) {
                    const auto& edge=grid.edges[edgeId];
                    const int i=currentCell%grid.nx,j=currentCell/grid.nx;
                    // This no-declaration path requires exactly one regular
                    // interior crossing per cut edge; there are no corner roots.
                    return PanelGraph::Endpoint{edge.horizontal ? (edge.j==j?1:3) : (edge.i==i?0:2),0};
                };
                graphs_.back().beginRoot=endpointIdentity(currentEdge);
                graphs_.back().endRoot=endpointIdentity(endEdge);
            }
        }
        catch (const std::exception& e)
        {
            throw std::runtime_error("[LevelSetInterfacePanels] cellIndex=" + std::to_string(solverCell) + ": " + e.what());
        }
        mapping[solverCell] = static_cast<int>(p);
        const auto& incident = grid.edges[endEdge].cells;
        if (!closed_ && incident.size() == 1)
        {
            require(p + 1 == cutCells && endEdge != startEdge &&
                (endEdge == axisEnds.front() || endEdge == axisEnds.back()),
                "axis chain ends early or leaves disconnected interface components");
            currentEdge = endEdge;
            break;
        }
        require(incident.size() == 2, "broken interior interface incidence");
        currentCell = incident[0] == currentCell ? incident[1] : incident[0];
        currentEdge = endEdge;
        if (p + 1 < cutCells)
            require(currentEdge != startEdge, "multiple disconnected interface components are unsupported");
    }
    if (closed_)
        require(currentEdge == startEdge && currentCell == firstCell
            && same(panels.back().end, panels.front().begin), "interface graph did not form one closed ordered loop");
    else
        require(panels.size() == cutCells && panels.front().begin.x == 0.0 && panels.back().end.x == 0.0
            && currentEdge != startEdge, "interface graph did not form one axis-to-axis ordered chain");
    panels_ = std::move(panels);
    cellToPanel_ = std::move(mapping);
    if (!closed_) beginEndpointKind_ = endEndpointKind_ = EndpointKind::Axis;
}

const LevelSetInterfacePanels::PanelGraph& LevelSetInterfacePanels::panelGraph(int panelIndex) const
{
    require(static_cast<bool>(source_), "subarc queries require explicitly retained source construction");
    require(panelIndex >= 0 && static_cast<std::size_t>(panelIndex) < panels_.size()
        && graphs_.size() == panels_.size(), "invalid retained panel graph index");
    const auto& chart = graphs_[panelIndex];
    require(source_->cellActive(chart.cellI, chart.cellJ), "retained panel source no longer supplies its active cell");
    return chart;
}

std::vector<LevelSetInterfacePanels::ArcSample> LevelSetInterfacePanels::subarcQuadrature(int panelIndex,
    double localCoordinateBegin, double localCoordinateEnd) const
{
    const auto& chart = panelGraph(panelIndex);
    const double first = chart.coordinateIsX ? chart.localBegin.x : chart.localBegin.y;
    const double last = chart.coordinateIsX ? chart.localEnd.x : chart.localEnd.y;
    for (double coordinate : {localCoordinateBegin, localCoordinateEnd})
        require(std::isfinite(coordinate) && coordinate >= std::min(first, last)
            && coordinate <= std::max(first, last), "subarc endpoint lies outside its retained local graph");
    std::vector<ArcSample> result;
    if (localCoordinateBegin == localCoordinateEnd) return result;
    CurveGraph graph{*source_, {0.0,0.0}, {1.0,1.0}, chart.localBegin, chart.localEnd,
        chart.coordinateIsX, last > first ? 1.0 : -1.0};
    graph.cellLocal = true; graph.cellI = chart.cellI; graph.cellJ = chart.cellJ; graph.spacing = chart.cellSize;
    if (chart.coordinateIsX) { graph.begin.x = localCoordinateBegin; graph.end.x = localCoordinateEnd; }
    else { graph.begin.y = localCoordinateBegin; graph.end.y = localCoordinateEnd; }
    const long double span = graph.physicalSpan();
    return LAMP::Geometry::Detail::InterfaceQuadrature::sampleGraph(span,axisymmetric_,
        [&](double fraction) { return graph.at(fraction); },
        [&](ArcSample& node) { node.localGraphCoordinate=chart.coordinateIsX?
            node.frame.cellLocalPoint.x:node.frame.cellLocalPoint.y; },
        [&](const Frame& frame) { return graph.metric(frame); },
        [&](const Frame& frame) { return source_->scalarJetAtCellLocal(chart.cellI,chart.cellJ,frame.cellLocalPoint); });
}

namespace
{
namespace ProjectionExpansion = LAMP::Geometry::Detail::ExactExpansion;
using Expansion = ProjectionExpansion::Expansion;
using ProjectionQuery = LevelSetInterfacePanels::LocalProjectionQuery;
using ProjectionResult = LevelSetInterfacePanels::LocalProjection;

struct ProjectionSample
{
    ProjectionResult result;
    double derivative = 0.0; // Signed derivative along the oriented interface.
};

// Root identity is topological. In particular, two different edge roots may
// have the same rounded physical display coordinate and must remain distinct.
struct ProjectionRoot
{
    int kind = -1, i = -1, j = -1, ordinal = -1; // vertex / vertical / horizontal
    bool operator==(const ProjectionRoot& other) const
    { return kind == other.kind && i == other.i && j == other.j && ordinal == other.ordinal; }
};

ProjectionRoot projectionRoot(const PanelGraph& graph, bool begin)
{
    const auto point = begin ? graph.localBegin : graph.localEnd;
    const auto root = begin ? graph.beginRoot : graph.endRoot;
    require(root.side >= 0 && root.side < 4 && root.rootOrdinal >= 0,
        "local projection requires the original endpoint root identity");
    require((root.side == 0 && point.x == 0.0) || (root.side == 1 && point.y == 0.0)
        || (root.side == 2 && point.x == 1.0) || (root.side == 3 && point.y == 1.0),
        "retained endpoint root side disagrees with its local coordinate");
    if ((point.x == 0.0 || point.x == 1.0) && (point.y == 0.0 || point.y == 1.0))
        return {0, graph.cellI + (point.x == 1.0), graph.cellJ + (point.y == 1.0), -1};
    if (root.side == 0 || root.side == 2)
        return {1, graph.cellI + (root.side == 2), graph.cellJ, root.rootOrdinal};
    return {2, graph.cellI, graph.cellJ + (root.side == 3), root.rootOrdinal};
}

class LocalProjectionContext
{
public:
    // Physical-query metrics only. Load projection has its own original
    // boundary KKT policy; this constructor does NOT enable accept(), a
    // material radius, motion allowance, or the material branch-box guard.
    LocalProjectionContext(const Snapshot& source, Vector2D physicalQuery)
        : source_(source), spacing_(source.cellSize())
    {
        ProjectionExpansion::validateEnvironment();
        require(finite(physicalQuery), "physical projection query must be finite");
        targetX_ = ProjectionExpansion::fromDouble(physicalQuery.x);
        targetY_ = ProjectionExpansion::fromDouble(physicalQuery.y);
    }

    LocalProjectionContext(const Snapshot& source, const ProjectionQuery& query)
        : source_(source), spacing_(source.cellSize())
    {
        ProjectionExpansion::validateEnvironment();
        require(query.cellI >= 0 && query.cellI < source.nx()
            && query.cellJ >= 0 && query.cellJ < source.ny()
            && finite(query.cellLocalPoint) && query.cellLocalPoint.x >= 0.0 && query.cellLocalPoint.x <= 1.0
            && query.cellLocalPoint.y >= 0.0 && query.cellLocalPoint.y <= 1.0 && finite(query.displacement),
            "local projection requires a finite closed-cell anchor and displacement");
        const double h = std::min(spacing_.x, spacing_.y);
        radius = h / 4.0;
        tolerance = 1e-9 * h;
        const auto lower = source.minCoord(), upper = source.maxCoord();
        // Same branch-box margin as the existing material projection. It is
        // never used to merge roots, clip a coordinate, or accept stationarity.
        boxTolerance = 512.0 * std::numeric_limits<double>::epsilon()
            * std::max({spacing_.x, spacing_.y, std::abs(lower.x), std::abs(lower.y),
                std::abs(upper.x), std::abs(upper.y)});
        require(std::isfinite(radius) && radius > 0.0 && std::isfinite(tolerance) && tolerance > 0.0
            && std::isfinite(boxTolerance) && std::hypot(query.displacement.x, query.displacement.y) <= h / 2.0,
            "local material displacement exceeds half a source cell");
        targetX_ = ProjectionExpansion::add(component(query.cellI, query.cellLocalPoint.x, true),
            ProjectionExpansion::fromDouble(query.displacement.x));
        targetY_ = ProjectionExpansion::add(component(query.cellJ, query.cellLocalPoint.y, false),
            ProjectionExpansion::fromDouble(query.displacement.y));
    }

    ProjectionSample evaluate(int panelIndex, const PanelGraph& graph, Frame frame,
        bool endpoint = false, bool begin = false) const
    {
        require((frame.hasCellLocalPoint || frame.hasCellChartPoint)
            && frame.cellI == graph.cellI && frame.cellJ == graph.cellJ,
            "local projection lost its actual source-jet provenance");
        const auto dx = frameDifference(frame, true);
        const auto dy = frameDifference(frame, false);
        const auto tangential = ProjectionExpansion::add(
            ProjectionExpansion::multiply(dx, ProjectionExpansion::fromDouble(-frame.normal.y)),
            ProjectionExpansion::multiply(dy, ProjectionExpansion::fromDouble(frame.normal.x)));
        ProjectionSample sample;
        sample.derivative = ProjectionExpansion::estimate(tangential);
        auto& result = sample.result;
        result.frame = std::move(frame);
        result.frame.stationarityResidual = std::abs(sample.derivative);
        result.panelIndex = panelIndex;
        result.endpoint = endpoint; result.atPanelBegin = begin;
        result.distance = std::hypot(ProjectionExpansion::estimate(dx), ProjectionExpansion::estimate(dy));
        result.inwardDerivative = endpoint ? (begin ? sample.derivative : -sample.derivative) : 0.0;
        result.kktResidual = endpoint ? std::max(0.0, -result.inwardDerivative) : std::abs(sample.derivative);
        if (endpoint) result.root = begin ? graph.beginRoot : graph.endRoot;
        require(std::isfinite(result.distance) && std::isfinite(result.frame.levelSetResidual)
            && std::isfinite(result.frame.stationarityResidual) && std::isfinite(result.kktResidual),
            "non-finite local projection residual");
        return sample;
    }

    ProjectionSample endpoint(int index, const PanelGraph& graph, bool begin) const
    {
        projectionRoot(graph, begin); // Validate the original identity, never infer it from distance.
        return evaluate(index, graph, source_.jetAtCellLocal(graph.cellI, graph.cellJ,
            begin ? graph.localBegin : graph.localEnd), true, begin);
    }

    double boxDistance(const PanelGraph& graph, Vector2D low, Vector2D high) const
    {
        const auto distance = [&](int index, double first, double last, bool x)
        {
            const auto a = difference(index, first, x), b = difference(index, last, x);
            if (ProjectionExpansion::sign(a) > 0) return ProjectionExpansion::estimate(a);
            if (ProjectionExpansion::sign(b) < 0) return -ProjectionExpansion::estimate(b);
            return 0.0;
        };
        return std::hypot(distance(graph.cellI, low.x, high.x, true),
            distance(graph.cellJ, low.y, high.y, false));
    }

    double targetCoordinate(const PanelGraph& graph) const
    {
        const bool x = graph.coordinateIsX;
        const int index = x ? graph.cellI : graph.cellJ;
        // Convert only a candidate section seed, never an endpoint identity.
        // An unrepresentable interior seed is omitted; no clipping is used.
        return ProjectionExpansion::estimate(ProjectionExpansion::subtract(x ? targetX_ : targetY_,
            ProjectionExpansion::fromDouble(coordinate(index, x)))) / (x ? spacing_.x : spacing_.y);
    }

    int compareDistance(const Frame& first, const Frame& second) const
    {
        return ProjectionExpansion::sign(ProjectionExpansion::subtract(
            squaredDistance(first), squaredDistance(second)));
    }

    bool samePosition(const Frame& first, const Frame& second) const
    {
        return ProjectionExpansion::sign(ProjectionExpansion::subtract(
                frameDifference(first, true), frameDifference(second, true))) == 0
            && ProjectionExpansion::sign(ProjectionExpansion::subtract(
                frameDifference(first, false), frameDifference(second, false))) == 0;
    }

    void accept(const ProjectionResult& result) const
    {
        require(result.distance <= radius, "local material projection exceeds a quarter source cell");
        require(result.frame.levelSetResidual <= tolerance && result.kktResidual <= tolerance,
            "local material projection exceeds the original level-set/stationarity or endpoint KKT tolerance");
    }

    double radius = 0.0, tolerance = 0.0, boxTolerance = 0.0;

private:
    const Snapshot& source_;
    Vector2D spacing_;
    Expansion targetX_, targetY_;

    double coordinate(int index, bool x) const
    {
        const int count = x ? source_.nx() : source_.ny();
        if (index == count) return x ? source_.maxCoord().x : source_.maxCoord().y;
        return (x ? source_.minCoord().x : source_.minCoord().y) + index * (x ? spacing_.x : spacing_.y);
    }
    Expansion component(int index, double local, bool x) const
    {
        // Match Snapshot::jetAtCellLocal's boundary convention. Adjacent cell
        // origins need not differ by exactly nominal h on a nonbinary grid.
        if (local == 1.0) return ProjectionExpansion::fromDouble(coordinate(index + 1, x));
        return ProjectionExpansion::add(ProjectionExpansion::fromDouble(coordinate(index, x)),
            ProjectionExpansion::multiply(ProjectionExpansion::fromDouble(x ? spacing_.x : spacing_.y),
                ProjectionExpansion::fromDouble(local)));
    }
    Expansion difference(int index, double local, bool x) const
    { return ProjectionExpansion::subtract(component(index, local, x), x ? targetX_ : targetY_); }
    Expansion frameDifference(const Frame& frame, bool x) const
    {
        if (!frame.hasCellChartPoint)
            return difference(x ? frame.cellI : frame.cellJ,
                x ? frame.cellLocalPoint.x : frame.cellLocalPoint.y, x);
        const auto chart = LevelSetHermiteGeometry::cellChartPoint(frame);
        const auto local = x ? chart.x : chart.y;
        const int index = x ? frame.cellI : frame.cellJ;
        if (!local.fromUpper) return difference(index, local.offset, x);
        // Exact cell boundaries keep the legacy actual-grid convention.
        // An interior upper-relative offset is NEVER first rounded to 1-tail.
        if (local.offset == 0.0) return difference(index, 1.0, x);
        if (local.offset == 1.0) return difference(index, 0.0, x);
        const auto spacing = ProjectionExpansion::fromDouble(x ? spacing_.x : spacing_.y);
        const auto position = ProjectionExpansion::add(
            ProjectionExpansion::fromDouble(coordinate(index, x)),
            ProjectionExpansion::subtract(spacing,
                ProjectionExpansion::multiply(spacing,ProjectionExpansion::fromDouble(local.offset))));
        return ProjectionExpansion::subtract(position,x ? targetX_ : targetY_);
    }
    Expansion squaredDistance(const Frame& frame) const
    {
        const auto dx = frameDifference(frame, true);
        const auto dy = frameDifference(frame, false);
        return ProjectionExpansion::add(ProjectionExpansion::multiply(dx, dx),
            ProjectionExpansion::multiply(dy, dy));
    }
};

CurveGraph projectionGraph(const Snapshot& source, const PanelGraph& chart)
{
    const double first = chart.coordinateIsX ? chart.localBegin.x : chart.localBegin.y;
    const double last = chart.coordinateIsX ? chart.localEnd.x : chart.localEnd.y;
    CurveGraph result{source, {0.0, 0.0}, {1.0, 1.0}, chart.localBegin, chart.localEnd,
        chart.coordinateIsX, last > first ? 1.0 : -1.0};
    result.cellLocal = true; result.cellI = chart.cellI; result.cellJ = chart.cellJ; result.spacing = chart.cellSize;
    return result;
}

double projectionCoordinate(const Frame& frame, const PanelGraph& graph)
{
    if (!frame.hasCellChartPoint)
        return graph.coordinateIsX ? frame.cellLocalPoint.x : frame.cellLocalPoint.y;
    const auto point = LevelSetHermiteGeometry::cellChartPoint(frame);
    return LAMP::Geometry::cellChartCoordinateValue(graph.coordinateIsX ? point.x : point.y);
}

namespace SharedProjection = LAMP::Geometry::Detail::LocalInterfaceProjection;
using SharedProjection::projectionNeighbor;

// A sign bracket is in increasing ORIENTED arclength, not necessarily increasing
// chosen coordinate. Resolve it in the original parameter interval. Merely
// being within the physical residual tolerance does not turn a tiny panel's
// arbitrary quadrature node into the material projection.
template<class Evaluate>
ProjectionSample stationaryBracket(const LocalProjectionContext& context, int index,
    const PanelGraph& chart, ProjectionSample first, ProjectionSample last,
    double stationaryTolerance, const Evaluate& evaluate)
{
    return SharedProjection::stationaryMinimum(std::move(first),std::move(last),stationaryTolerance,
        [&](const ProjectionSample& sample) { return projectionCoordinate(sample.result.frame,chart); },
        [&](double coordinate) { return context.evaluate(index,chart,evaluate(coordinate)); });
}


class UniformProjectionLayout
{
public:
    using Sample = ProjectionSample;
    UniformProjectionLayout(const Snapshot& source,const std::vector<PanelGraph>& graphs,bool closed)
        : source_(source),graphs_(graphs),closed_(closed) {}
    int count() const { return static_cast<int>(graphs_.size()); }
    bool closed() const { return closed_; }
    ProjectionRoot identity(int index,bool begin) const { return projectionRoot(graphs_[index],begin); }
    Sample endpoint(const LocalProjectionContext& context,int index,bool begin) const
    { return context.endpoint(index,graphs_[index],begin); }
    double boxDistance(const LocalProjectionContext& context,int index) const
    { return context.boxDistance(graphs_[index],{0.0,0.0},{1.0,1.0}); }
    std::vector<Sample> samples(const LocalProjectionContext& context,int index) const
    {
        const auto& chart=graphs_[index];
        const auto curve=projectionGraph(source_,chart);
        std::vector<Sample> samples;
        samples.push_back(context.endpoint(index,chart,true));
        for (double abscissa : gaussX)
            samples.push_back(context.evaluate(index,chart,curve.at((1.0+abscissa)/2.0)));
        samples.push_back(context.endpoint(index,chart,false));
        const double a=projectionCoordinate(samples.front().result.frame,chart);
        const double b=projectionCoordinate(samples.back().result.frame,chart);
        const double seed=context.targetCoordinate(chart);
        if (std::isfinite(seed) && seed>std::min(a,b) && seed<std::max(a,b))
            samples.push_back(context.evaluate(index,chart,curve.atCoordinate(seed)));
        std::sort(samples.begin(),samples.end(),[&](const auto& left,const auto& right) {
            return curve.direction*projectionCoordinate(left.result.frame,chart)
                < curve.direction*projectionCoordinate(right.result.frame,chart);
        });
        return samples;
    }
    Sample minimum(const LocalProjectionContext& context,int index,Sample first,Sample last) const
    {
        const auto& chart=graphs_[index];
        const auto curve=projectionGraph(source_,chart);
        return stationaryBracket(context,index,chart,std::move(first),std::move(last),context.tolerance,
            [&](double coordinate) { return curve.atCoordinate(coordinate); });
    }
private:
    const Snapshot& source_;
    const std::vector<PanelGraph>& graphs_;
    bool closed_;
};

// A load foot may be much nearer a cell side than 1-local can represent.
// Keep its section root in the endpoint-relative chart all the way to the
// jet. The retained graph's chosen coordinate still has its original double
// parameterization; an unresolved chosen parameter rejects rather than snaps.
Frame loadSectionFrame(const Snapshot& source, const PanelGraph& graph, double coordinate)
{
    const double first = graph.coordinateIsX ? graph.localBegin.x : graph.localBegin.y;
    const double last = graph.coordinateIsX ? graph.localEnd.x : graph.localEnd.y;
    require(std::isfinite(coordinate) && coordinate > std::min(first,last)
        && coordinate < std::max(first,last), "load graph coordinate is not strictly interior");
    const auto roots = source.intersectCellSection(graph.cellI,graph.cellJ,
        graph.coordinateIsX,{coordinate,false});
    require(roots.size() == 1, "load section is not one resolved retained curve graph");
    require(LAMP::Geometry::cellChartInterior(roots.front()),
        "load section root is not inside its retained cell chart");
    auto frame = source.jetAtCellChart(graph.cellI,graph.cellJ,roots.front());
    require(frame.hasCellChartPoint && !frame.hasCellLocalPoint
        && projectionCoordinate(frame,graph) == coordinate, "load section lost its original chosen coordinate");
    projectionGraph(source,graph).metric(frame); // Original oriented regular-graph check.
    return frame;
}

class UniformLoadProjectionLayout
{
public:
    using Sample = ProjectionSample;
    UniformLoadProjectionLayout(const Snapshot& source, const std::vector<PanelGraph>& graphs, bool closed)
        : source_(source), graphs_(graphs), closed_(closed) {}
    int count() const { return static_cast<int>(graphs_.size()); }
    bool closed() const { return closed_; }
    double endpointDistance(const LocalProjectionContext& context,int index,bool begin) const
    { return endpoint(context,index,begin).result.distance; }
    ProjectionRoot identity(int index, bool begin) const { return projectionRoot(graphs_[index],begin); }
    Sample endpoint(const LocalProjectionContext& context, int index, bool begin) const
    { return context.endpoint(index,graphs_[index],begin); }
    double boxDistance(const LocalProjectionContext& context, int index) const
    { return context.boxDistance(graphs_[index],{0.0,0.0},{1.0,1.0}); }
    std::vector<Sample> samples(const LocalProjectionContext& context, int index) const
    {
        const auto& graph = graphs_[index];
        const auto curve = projectionGraph(source_,graph);
        const auto evaluate = [&](double coordinate) { return loadSectionFrame(source_,graph,coordinate); };
        std::vector<ProjectionSample> samples;
        samples.push_back(context.endpoint(index,graph,true));
        for (double abscissa : gaussX)
        {
            const double first = graph.coordinateIsX ? graph.localBegin.x : graph.localBegin.y;
            const double last = graph.coordinateIsX ? graph.localEnd.x : graph.localEnd.y;
            const double coordinate = checked(static_cast<long double>(first)
                + ((1.0+abscissa)/2.0)*(static_cast<long double>(last)-first), "load graph sample overflow");
            samples.push_back(context.evaluate(index,graph,evaluate(coordinate)));
        }
        samples.push_back(context.endpoint(index,graph,false));
        const double a = projectionCoordinate(samples.front().result.frame,graph);
        const double b = projectionCoordinate(samples.back().result.frame,graph);
        const double seed = context.targetCoordinate(graph);
        if (std::isfinite(seed) && seed > std::min(a,b) && seed < std::max(a,b))
            samples.push_back(context.evaluate(index,graph,evaluate(seed)));
        std::sort(samples.begin(),samples.end(),[&](const auto& left,const auto& right) {
            return curve.direction*projectionCoordinate(left.result.frame,graph)
                < curve.direction*projectionCoordinate(right.result.frame,graph);
        });
        return samples;
    }
    Sample loadMinimum(const LocalProjectionContext& context, int index, Sample first, Sample last,
        double tolerance) const
    {
        const auto& graph=graphs_[index];
        return stationaryBracket(context,index,graph,std::move(first),std::move(last),tolerance,
            [&](double coordinate) { return loadSectionFrame(source_,graph,coordinate); });
    }
private:
    const Snapshot& source_;
    const std::vector<PanelGraph>& graphs_;
    bool closed_;
};

std::vector<ProjectionResult> retainedLoadCandidates(const Snapshot& source,
    const std::vector<PanelGraph>& graphs, bool closed, const LocalProjectionContext& context,
    double boundaryTolerance)
{
    return SharedProjection::loadCandidates(context,UniformLoadProjectionLayout(source,graphs,closed),
        boundaryTolerance);
}
} // namespace

LevelSetInterfacePanels::LocalProjection LevelSetInterfacePanels::projectLocalEndpoint(
    const LocalProjectionQuery& query, int panelIndex, bool atPanelBegin) const
{
    const auto& graph = panelGraph(panelIndex);
    const LocalProjectionContext context(*source_, query);
    auto result = context.endpoint(panelIndex, graph, atPanelBegin).result;
    context.accept(result);
    SharedProjection::guardBranches(context,UniformProjectionLayout(*source_,graphs_,closed_),result);
    return result;
}

LevelSetInterfacePanels::LocalProjection LevelSetInterfacePanels::projectLocal(const LocalProjectionQuery& query) const
{
    require(source_ && !graphs_.empty() && graphs_.size()==panels_.size(),
        "local projection requires explicitly retained source construction");
    const LocalProjectionContext context(*source_,query);
    return SharedProjection::search(context,UniformProjectionLayout(*source_,graphs_,closed_));
}

LevelSetInterfacePanels::LocalProjection LevelSetInterfacePanels::projectedFrame(Vector2D query) const
{
    if (nativeChain_) return nativeChain_->projectedFrame(query);
    require(source_ && !graphs_.empty() && graphs_.size() == panels_.size(),
        "load projection requires explicitly retained source construction");
    // This private inventory runs only the Snapshot's original eight closest
    // seeds. No arithmetic/geometry exceptions are swallowed to obtain an
    // endpoint fallback, and the source validates actual query activity.
    const auto stationary = source_->stationaryCellLocalFrames(query);
    const LocalProjectionContext context(*source_, query);
    const auto h = source_->cellSize(), lower = source_->minCoord(), upper = source_->maxCoord();
    const long double length = std::max(h.x, h.y);
    const long double coordinateScale = std::max({length,
        std::abs(static_cast<long double>(lower.x)), std::abs(static_cast<long double>(lower.y)),
        std::abs(static_cast<long double>(upper.x)), std::abs(static_cast<long double>(upper.y))});
    const long double tolerance = 2048.0L * std::numeric_limits<double>::epsilon() * coordinateScale;
    require(std::isfinite(tolerance) && tolerance > 0.0L
        && tolerance <= std::sqrt(static_cast<long double>(std::numeric_limits<double>::epsilon()))
            * std::min(h.x, h.y), "coordinate precision cannot resolve the requested boundary projection");

    LocalProjection best;
    bool found = false;
    const auto consider = [&](LocalProjection candidate)
    {
        const int order = found ? context.compareDistance(candidate.frame, best.frame) : -1;
        // Preserve the boundary selector's endpoint preference only for the
        // SAME actual local position. Equal rounded display coordinates are
        // not an endpoint/root identity, nor an excuse to merge a short arc.
        if (!found || order < 0 || (order == 0 && candidate.endpoint
            && context.samePosition(candidate.frame, best.frame)))
        {
            best = std::move(candidate);
            found = true;
        }
    };
    for (const auto& frame : stationary.frames)
    {
        require(frame.hasCellLocalPoint, "stationary load candidate lost local provenance");
        const auto panel = std::find_if(graphs_.begin(), graphs_.end(), [&](const auto& graph)
            { return graph.cellI == frame.cellI && graph.cellJ == frame.cellJ; });
        // In particular, a pure wet-sector corner may have a tiny Newton
        // residual but has no free-curve panel incidence. It is not a load
        // point. Other original converged seeds are still inspected.
        if (panel == graphs_.end()) continue;
        const int index = static_cast<int>(panel - graphs_.begin());
        const auto local = frame.cellLocalPoint;
        const double q = panel->coordinateIsX ? local.x : local.y;
        const double first = panel->coordinateIsX ? panel->localBegin.x : panel->localBegin.y;
        const double last = panel->coordinateIsX ? panel->localEnd.x : panel->localEnd.y;
        const bool originalRoot = same(local, panel->localBegin) || same(local, panel->localEnd);
        const bool interior = finite(local) && local.x > 0.0 && local.x < 1.0
            && local.y > 0.0 && local.y < 1.0 && q > std::min(first, last) && q < std::max(first, last);
        if (!originalRoot && !interior) continue; // No global bbox, tolerance or endpoint snapping.
        const auto actual = source_->jetAtCellLocal(frame.cellI, frame.cellJ, local);
        require(same(frame.point, actual.point) && same(frame.normal, actual.normal)
            && frame.curvature == actual.curvature, "stationary load candidate does not match its actual source jet");
        auto candidate = context.evaluate(index, *panel, frame).result;
        if (candidate.frame.levelSetResidual > tolerance || candidate.frame.stationarityResidual > tolerance) continue;
        // A stationary shared-panel root is not a constrained open boundary.
        // The actual open endpoints below carry their original root metadata.
        consider(std::move(candidate));
    }
    if (!closed_)
    {
        const int indices[2] = {0, static_cast<int>(graphs_.size()-1)};
        const EndpointKind kinds[2] = {beginEndpointKind_, endEndpointKind_};
        for (int end = 0; end < 2; ++end)
        {
            require(kinds[end] == EndpointKind::DeclaredBoundary || kinds[end] == EndpointKind::Axis,
                "load projection encountered an undeclared open-chain endpoint");
            const int index = indices[end];
            auto candidate = context.endpoint(index, panelGraph(index), end == 0).result;
            require(candidate.frame.levelSetResidual <= tolerance,
                "original local endpoint has an unresolved scalar residual");
            if (kinds[end] == EndpointKind::Axis)
                require(axisymmetric_ && candidate.frame.point.x == 0.0 && candidate.frame.normal.x == 0.0,
                    "original load-projection pole lost radial symmetry");
            if (candidate.inwardDerivative < -tolerance) continue;
            candidate.frame.attemptedSeeds = stationary.attemptedSeeds;
            candidate.frame.convergedSeeds = stationary.convergedSeeds;
            consider(std::move(candidate));
        }
    }
    if (!found)
        for (auto candidate : retainedLoadCandidates(*source_,graphs_,closed_,context,
            checked(tolerance,"unrepresentable load boundary tolerance")))
        {
            candidate.frame.attemptedSeeds = stationary.attemptedSeeds;
            candidate.frame.convergedSeeds = stationary.convergedSeeds;
            consider(std::move(candidate));
        }
    if (!found)
    {
        std::ostringstream message;
        message << std::setprecision(17)
            << "load projection has no retained-curve stationary or admissible original endpoint candidate"
            << "; query=(" << query.x << ',' << query.y << "), tolerance=" << tolerance
            << ", attemptedSeeds=" << stationary.attemptedSeeds
            << ", convergedSeeds=" << stationary.convergedSeeds;
        for (const auto& frame : stationary.frames)
            message << "; seedFrame cell=(" << frame.cellI << ',' << frame.cellJ
                << "), local=(" << frame.cellLocalPoint.x << ',' << frame.cellLocalPoint.y
                << "), phiResidual=" << frame.levelSetResidual
                << ", stationaryResidual=" << frame.stationarityResidual;
        if (!closed_)
            for (int end = 0; end < 2; ++end)
            {
                const int index = end == 0 ? 0 : static_cast<int>(graphs_.size()-1);
                const auto endpoint = context.endpoint(index,panelGraph(index),end == 0).result;
                message << "; endpoint=" << end << ", cell=(" << endpoint.frame.cellI
                    << ',' << endpoint.frame.cellJ << "), local=(" << endpoint.frame.cellLocalPoint.x
                    << ',' << endpoint.frame.cellLocalPoint.y << "), inwardDerivative=" << endpoint.inwardDerivative;
            }
        require(false,message.str());
    }
    return best;
}
