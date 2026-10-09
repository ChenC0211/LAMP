#include "Modules/Hydrodynamics/HydroInterfaceLoadGeometry.h"
#include "Core/Geometry/LevelSetPhaseFraction.h"
#include "Core/Geometry/LevelSetInterfacePanels.h"
#include "Core/Field/VertField.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/UniformLeafGrid.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace
{
void require(bool ok, const char* reason)
{
    if (!ok) throw std::runtime_error(std::string("[HydroInterfaceLoadGeometry] ") + reason);
}
bool finite(Vector2D point) { return std::isfinite(point.x) && std::isfinite(point.y); }
UniformLeafGrid validatedGrid(const QuadtreeMesh& mesh, CartesianGridScope scope)
{
    try { return UniformLeafGrid(mesh, scope); }
    catch (const std::invalid_argument& error)
    {
        // Keep the load-geometry rejection port while sharing topology checks.
        throw std::runtime_error(std::string("[HydroInterfaceLoadGeometry] ") + error.what());
    }
}
struct Source
{
    const std::vector<double>& fraction; // Same immutable Core diagnostic fraction.
    std::vector<double> indicator; // Empty for PL: reuse fraction without a second definition.
    std::vector<HydroInterfaceLoadGeometry::Face> faces;
};
Source inspectSource(const QuadtreeMesh& mesh, const VertField<double>& phi,
    HydroInterfaceLoadGeometry::Options options, const InterfaceGeometry& geometry)
{
    const auto* snapshot = &geometry.snapshot();
    using Layer = HydroInterfaceLoadGeometry::Layer;
    require(options.layer == Layer::GeometricPL || options.layer == Layer::Regularized, "unknown load layer");
    const bool regularized = options.layer == Layer::Regularized;
    const bool boundaryContact = !options.normalVelocityBoundaries.empty() || !options.boundaryContacts.empty();
    require(!boundaryContact || (regularized && !options.includeTangentialLoads
        && !options.normalVelocityBoundaries.empty() && !options.boundaryContacts.empty()),
        "boundary contact requires regularized normal loads, declared contacts and constrained normal-velocity intervals");
    require(!options.axisymmetric || regularized,
        "axisymmetric load geometry requires Regularized; RZ geometric volume fractions are unsupported");
    require(std::isfinite(options.halfWidth) && (regularized ? options.halfWidth > 0.0 : options.halfWidth == 0.0),
        "regularized load layer requires positive finite half-width; geometricPL requires zero");
    require(!regularized || snapshot, "regularized load layer requires the shared scalar snapshot");
    require(!options.includeTangentialLoads || snapshot, "tangential load geometry requires the shared scalar snapshot");
    const UniformLeafGrid grid = validatedGrid(mesh, options.gridScope);
    const int nx = grid.nx(), ny = grid.ny();
    require(nx >= 4 && ny >= 4, "at least 4x4 uniform Cartesian cells are required");
    const long long cells = static_cast<long long>(nx) * ny;
    const long long vertices = (static_cast<long long>(nx) + 1) * (ny + 1);
    const long long xFaces = (static_cast<long long>(nx) + 1) * ny;
    const long long yFaces = static_cast<long long>(nx) * (ny + 1);
    require(phi.value && phi.size() == vertices, "level-set storage must match the mesh vertices");
    const double hx = grid.spacing().x, hy = grid.spacing().y, area = hx * hy;
    require(finite(mesh.minCoord) && finite(mesh.maxCoord)
        && finite(grid.spacing()) && hx > 0. && hy > 0.
        && mesh.maxCoord.x > mesh.minCoord.x && mesh.maxCoord.y > mesh.minCoord.y
        && std::isfinite(area) && area > 0., "invalid Cartesian domain, spacing or dual area");
    require(!options.axisymmetric || mesh.minCoord.x >= 0.0,
        "axisymmetric load geometry requires nonnegative domain radius");
    const double coordinateScale = std::max({hx, hy, std::abs(mesh.minCoord.x), std::abs(mesh.minCoord.y),
        std::abs(mesh.maxCoord.x), std::abs(mesh.maxCoord.y)});
    const double tolerance = 256. * std::numeric_limits<double>::epsilon() * coordinateScale;
    require(std::isfinite(tolerance) && tolerance < std::min(hx, hy) / 16.,
        "Cartesian spacing is unresolved at this coordinate scale");
    const auto near = [&](double actual, long double expected)
    {
        return std::isfinite(actual) && std::isfinite(expected)
            && std::abs(static_cast<long double>(actual) - expected) <= tolerance;
    };
    const auto metricNear = [](double actual, double expected)
    {
        return std::isfinite(actual) && actual > 0.
            && std::abs(static_cast<long double>(actual) - expected)
                <= 64.L * std::numeric_limits<double>::epsilon() * expected;
    };
    require(near(mesh.maxCoord.x, static_cast<long double>(mesh.minCoord.x) + nx * hx)
        && near(mesh.maxCoord.y, static_cast<long double>(mesh.minCoord.y) + ny * hy),
        "domain extents disagree with Cartesian spacing/counts");
    const auto& cellGrid = grid.cells();
    // Reverse logical positions preserve the original native-index face
    // iteration order without reconstructing a second coordinate lattice.
    std::vector<int> xPositionByIndex(static_cast<std::size_t>(xFaces));
    std::vector<int> yPositionByIndex(static_cast<std::size_t>(yFaces));
    for (std::size_t p = 0; p < grid.vertical().size(); ++p)
        xPositionByIndex[mesh.vEdges(grid.vertical()[p]).index()] = static_cast<int>(p);
    for (std::size_t p = 0; p < grid.horizontal().size(); ++p)
        yPositionByIndex[mesh.hEdges(grid.horizontal()[p]).index()] = static_cast<int>(p);
    for (int k = 0; k < mesh.vertIndexNum(); ++k)
    {
        const int label = mesh.vertIndexToLabel[k];
        require(label >= 0 && label < mesh.vertLabelCapacity(), "vertex label outside allocated storage");
        const auto& vertex = mesh.verts(label);
        require(vertex.label() == label && vertex.index() == k && phi.label(k) == label, "vertex index-label mapping mismatch");
        require(std::isfinite(phi.value[k]), "non-finite source level set");
    }
    // All actual labels and native incidences were checked by UniformLeafGrid.
    const auto cellAt = [&](int i, int j) { return cellGrid[i + nx * j]; };
    for (std::size_t k = 0; k < options.normalVelocityBoundaries.size(); ++k)
    {
        const auto& patch = options.normalVelocityBoundaries[k];
        require(patch.side == Boundary::west || patch.side == Boundary::east
            || patch.side == Boundary::south || patch.side == Boundary::north, "invalid constrained boundary side");
        require(!(options.axisymmetric && mesh.minCoord.x == 0.0 && patch.side == Boundary::west),
            "the symmetry axis is not a prescribed physical-wall interval");
        const bool horizontal = patch.side == Boundary::south || patch.side == Boundary::north;
        const double lo = horizontal ? mesh.minCoord.x : mesh.minCoord.y;
        const double hi = horizontal ? mesh.maxCoord.x : mesh.maxCoord.y;
        const double spacing = horizontal ? hx : hy;
        require(std::isfinite(patch.start) && std::isfinite(patch.end)
            && patch.start >= lo && patch.end <= hi && patch.start < patch.end,
            "invalid constrained normal-velocity interval");
        for (double coordinate : {patch.start, patch.end})
            require(near(coordinate, static_cast<long double>(lo)
                + std::round((static_cast<long double>(coordinate)-lo)/spacing)*spacing),
                "constrained boundary endpoints must align with native faces");
        for (std::size_t j = 0; j < k; ++j)
            require(patch.side != options.normalVelocityBoundaries[j].side
                || patch.start >= options.normalVelocityBoundaries[j].end
                || patch.end <= options.normalVelocityBoundaries[j].start,
                "constrained boundary intervals overlap");
    }
    const auto covered = [&](Boundary side, double begin, double end)
    {
        if (options.axisymmetric && snapshot && snapshot->radialAxisSymmetry() && side == Boundary::west) return true;
        for (const auto& patch : options.normalVelocityBoundaries)
            if (patch.side == side && (patch.start <= begin || near(patch.start, begin))
                && (patch.end >= end || near(patch.end, end))) return true;
        return false;
    };
    const auto constrainedCell = [&](const Cell& cell)
    {
        const auto lower = mesh.verts(cell.cellToVert[0]).coord();
        const auto upper = mesh.verts(cell.cellToVert[2]).coord();
        const auto boundary = cell.boundary();
        const bool west = boundary == Boundary::west || boundary == Boundary::southWest || boundary == Boundary::northWest;
        const bool east = boundary == Boundary::east || boundary == Boundary::southEast || boundary == Boundary::northEast;
        const bool south = boundary == Boundary::south || boundary == Boundary::southWest || boundary == Boundary::southEast;
        const bool north = boundary == Boundary::north || boundary == Boundary::northWest || boundary == Boundary::northEast;
        // At a domain corner, BOTH physical sides must be constrained (or the
        // west side must be the RZ axis); one declared side is not sufficient.
        return (!west || covered(Boundary::west, lower.y, upper.y))
            && (!east || covered(Boundary::east, lower.y, upper.y))
            && (!south || covered(Boundary::south, lower.x, upper.x))
            && (!north || covered(Boundary::north, lower.x, upper.x));
    };
    for (const auto point : options.boundaryContacts)
    {
        require(finite(point), "nonfinite boundary contact");
        const bool allowed = (point.x == mesh.minCoord.x && covered(Boundary::west, point.y, point.y))
            || (point.x == mesh.maxCoord.x && covered(Boundary::east, point.y, point.y))
            || (point.y == mesh.minCoord.y && covered(Boundary::south, point.x, point.x))
            || (point.y == mesh.maxCoord.y && covered(Boundary::north, point.x, point.x));
        require(allowed, "interface contact is not on a declared constrained normal-velocity boundary");
    }
    for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i)
    {
        const auto& cell = mesh.cells(cellAt(i, j));
        require(metricNear(cell.size().x, hx) && metricNear(cell.size().y, hy), "nonuniform cell dimensions");
    }
    Source result{geometry.positiveFraction()};
    if (regularized)
    {
        result.indicator.resize(static_cast<std::size_t>(cells));
        for (int k = 0; k < mesh.cellIndexNum(); ++k)
        {
            const auto& cell = mesh.cells(mesh.cellIndexToLabel[k]);
            const auto point = cell.centerCoord();
            // A sign-only certificate is insufficient for a finite load band.
            // Certify the entire boundary cell, not just a saturated center.
            const bool axisOnly = options.axisymmetric && snapshot->radialAxisSymmetry()
                && cell.boundary() == Boundary::west;
            if (cell.boundary() != Boundary::non && !axisOnly)
                require(snapshot->closedCellOutsideBand(point, options.halfWidth)
                    || (boundaryContact && constrainedCell(cell)),
                    "regularized load band is not certified separated from boundary cells");
            result.indicator[k] = LevelSetPhaseFraction::regularizedPositiveIndicator(
                snapshot->value(point), options.halfWidth);
        }
    }
    const auto& indicator = regularized ? result.indicator : result.fraction;
    std::vector<double> gradientX, gradientY;
    if (options.includeTangentialLoads)
    {
        // A missing cross-component face at a physical wall may be extended
        // by zero only on a certified pure H plateau. Regularized mode has
        // already required the stronger CLOSED-cell band-separation proof.
        // These queries and all extra arithmetic are absent in default mode.
        for (int k = 0; k < mesh.cellIndexNum(); ++k)
        {
            const auto& cell = mesh.cells(mesh.cellIndexToLabel[k]);
            if (cell.boundary() == Boundary::non) continue;
            // At the symmetry axis, G_r(H)=0 follows from the even scalar;
            // it needs no pure-phase extension. Axial-wall corners still do.
            if (options.axisymmetric && snapshot->radialAxisSymmetry()
                && cell.boundary() == Boundary::west) continue;
            const int sign = snapshot->closedCellSign(cell.centerCoord());
            require(sign != 0 && indicator[k] == (sign > 0 ? 1.0 : 0.0),
                "tangential load zero extension requires strictly separated boundary cells on an exact H plateau");
        }
        gradientX.assign(static_cast<std::size_t>(xFaces), 0.0);
        gradientY.assign(static_cast<std::size_t>(yFaces), 0.0);
        const auto pressureGradient = [&](int highLabel, int lowLabel, bool horizontal)
        {
            const auto& high = mesh.cells(highLabel);
            const auto& low = mesh.cells(lowLabel);
            const double difference = indicator[high.index()] - indicator[low.index()];
            const double gradient = difference / (horizontal ? high.size().y : high.size().x);
            require(std::isfinite(gradient) && (difference == 0.0 || gradient != 0.0),
                "transverse source pressure gradient is unrepresentable");
            return gradient;
        };
        for (int j = 0; j < ny; ++j)
            for (int i = 1; i < nx; ++i)
                gradientX[i + (nx + 1) * j] = pressureGradient(cellAt(i, j), cellAt(i - 1, j), false);
        for (int j = 1; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
                gradientY[i + nx * j] = pressureGradient(cellAt(i, j), cellAt(i, j - 1), true);
    }
    auto& active = result.faces;
    active.reserve(static_cast<std::size_t>(xFaces + yFaces));
    const auto inspectFace = [&](const Edge& edge, int label, int k, bool horizontal)
    {
        require(edge.label() == label && edge.index() == k && edge.isHorizontal() == horizontal,
            "edge/output index-label mapping or orientation mismatch");
        require(metricNear(edge.length(), horizontal ? hx : hy), "edge length disagrees with Cartesian metric");
        const auto point = edge.midpointCoord();
        const int position = horizontal ? yPositionByIndex[k] : xPositionByIndex[k];
        const int stride = horizontal ? nx : nx + 1;
        const int i = position % stride, j = position / stride;
        const int master = edge.masterCell, slave = edge.slaveCell;
        if (master == -1 || slave == -1) return;
        const auto& plus = mesh.cells(master);
        const auto& minus = mesh.cells(slave);
        const auto plusPoint = plus.centerCoord(), minusPoint = minus.centerCoord();
        const double spacing = horizontal ? plus.size().y : plus.size().x;
        require(near(horizontal ? plusPoint.y - minusPoint.y : plusPoint.x - minusPoint.x, spacing),
            "interior center spacing disagrees with the pressure gradient stencil");
        // Production pressureFaceData uses this same master-minus-slave
        // difference divided by the uniform cell size, not cell interpolation.
        const double difference = indicator[plus.index()] - indicator[minus.index()];
        double transverse = 0.0;
        if (options.includeTangentialLoads)
        {
            // Four surrounding faces of the other Cartesian component.
            // At x(i,j+1/2): y(i-1/2,j), y(i+1/2,j), and both at j+1.
            // At y(i+1/2,j): x(i,j-1/2), x(i+1,j-1/2), and both at j+1/2.
            const std::array<double, 4> surrounding = horizontal
                ? std::array<double, 4>{{gradientX[i + (nx + 1) * (j - 1)],
                    gradientX[i + 1 + (nx + 1) * (j - 1)],
                    gradientX[i + (nx + 1) * j], gradientX[i + 1 + (nx + 1) * j]}}
                : std::array<double, 4>{{gradientY[i - 1 + nx * j], gradientY[i + nx * j],
                    gradientY[i - 1 + nx * (j + 1)], gradientY[i + nx * (j + 1)]}};
            long double sum = 0.0L;
            for (double value : surrounding) sum += value;
            const long double average = sum / 4.0L;
            transverse = static_cast<double>(average);
            require(std::isfinite(sum) && std::isfinite(transverse)
                && (sum == 0.0L || average != 0.0L) && (average == 0.0L || transverse != 0.0),
                "transverse pressure-gradient average overflow or underflow");
        }
        if (difference == 0. && transverse == 0.) return;
        const double gradient = difference / spacing;
        require(std::isfinite(gradient) && (difference == 0. || gradient != 0.), "fraction gradient is unrepresentable");
        const double dualArea = plus.area();
        require(std::isfinite(dualArea) && dualArea > 0., "invalid active pressure-face dual area");
        HydroInterfaceLoadGeometry::Face face;
        face.point = point; face.edgeIndex = k; face.horizontal = horizontal;
        face.lowCellIndex = minus.index(); face.highCellIndex = plus.index();
        face.gradient = gradient; face.dualArea = dualArea;
        face.transverseGradient = transverse;
        face.dualMeasure = dualArea;
        if (options.axisymmetric)
        {
            require(point.x > 0.0, "active axisymmetric load face must have strictly positive radius");
            face.dualMeasure = dualArea * (2.0 * std::acos(-1.0) * point.x);
            require(std::isfinite(face.dualMeasure) && face.dualMeasure > 0.0,
                "axisymmetric load dual measure overflow or underflow");
        }
        active.push_back(face);
    };
    for (int k = 0; k < mesh.vEdgeIndexNum(); ++k)
    {
        const int label = mesh.vEdgeIndexToLabel[k];
        require(label >= 0 && label < mesh.vEdgeLabelCapacity(), "vertical edge label outside allocated storage");
        inspectFace(mesh.vEdges(label), label, k, false);
    }
    for (int k = 0; k < mesh.hEdgeIndexNum(); ++k)
    {
        const int label = mesh.hEdgeIndexToLabel[k];
        require(label >= 0 && label < mesh.hEdgeLabelCapacity(), "horizontal edge label outside allocated storage");
        inspectFace(mesh.hEdges(label), label, k, true);
    }
    return result;
}
} // namespace

struct HydroInterfaceLoadGeometry::Impl
{
    Options options;
    InterfaceGeometry geometry;
    const LevelSetHermiteGeometry::Snapshot& snapshot;
    Source source;
    Impl(const QuadtreeMesh& inputMesh, const VertField<double>& phi, Options inputOptions)
        : options(inputOptions),
          geometry(inputMesh, phi, inputOptions.axisymmetric, inputOptions.gridScope, inputOptions.reconstruction),
          snapshot(geometry.snapshot()), source(inspectSource(inputMesh, phi, options, geometry))
    {
        if (!options.boundaryContacts.empty())
        {
            // Validate the actual open zero-curve topology before accepting
            // endpoint candidates. Declarations do not invent a curve or pin it.
            const LevelSetInterfacePanels chain(inputMesh, snapshot, options.axisymmetric, options.boundaryContacts, options.gridScope);
            require(!chain.closed(), "declared boundary contacts require an open interface chain");
        }
        for (auto& face : source.faces)
        {
            if (options.boundaryContacts.empty()) face.frame = snapshot.projectedFrame(face.point);
            else
            {
                const auto projection = snapshot.boundaryProjectedFrame(face.point, options.boundaryContacts);
                face.frame = projection.frame;
                face.boundaryEndpointProjection = projection.endpoint;
                face.endpointKktResidual = projection.kktResidual;
            }
            face.curvatureTangentialDerivative = snapshot.reconstruction() == LevelSetHermiteGeometry::Reconstruction::EnrichedQuartic
                ? snapshot.curvatureTangentialDerivative(face.frame.point)
                : snapshot.reconstructedCurvatureTangentialDerivative(face.frame.point);
            require(finite(face.frame.point) && finite(face.frame.normal)
                && std::isfinite(face.frame.curvature) && std::isfinite(face.curvatureTangentialDerivative),
                "non-finite active interface frame");
            face.surfaceCurvature = face.frame.curvature;
            if (options.axisymmetric)
            {
                require(face.frame.point.x >= 0.0, "negative reconstructed interface radius");
                const bool atAxis = face.frame.point.x == 0.0;
                require(!atAxis || (snapshot.radialAxisSymmetry() && face.frame.normal.x == 0.0),
                    "axis contact requires a regular even scalar with zero radial normal");
                // For a regular surface of revolution, lim(n_r/r)=kappa_m
                // at a pole. This is the limit, not a radius floor.
                const double azimuthalCurvature = atAxis ? face.frame.curvature
                    : face.frame.normal.x / face.frame.point.x;
                require(std::isfinite(azimuthalCurvature)
                    && (face.frame.normal.x == 0.0 || azimuthalCurvature != 0.0),
                    "axisymmetric azimuthal curvature overflow or underflow");
                face.surfaceCurvature = face.frame.curvature + azimuthalCurvature;
                require(std::isfinite(face.surfaceCurvature), "non-finite axisymmetric surface curvature");
            }
        }
    }
};

HydroInterfaceLoadGeometry::HydroInterfaceLoadGeometry(const QuadtreeMesh& mesh, const VertField<double>& phi)
    : HydroInterfaceLoadGeometry(mesh, phi, Options{}) {}
HydroInterfaceLoadGeometry::HydroInterfaceLoadGeometry(const QuadtreeMesh& mesh, const VertField<double>& phi, Options options)
    : impl_(std::make_unique<const Impl>(mesh, phi, options))
{
    std::vector<InterfaceTraceGeometry::Point> queries;
    queries.reserve(impl_->source.faces.size());
    for (const auto& face : impl_->source.faces)
        queries.push_back({face.frame, face.curvatureTangentialDerivative, face.endpointKktResidual});
    traceGeometry_ = std::make_shared<const InterfaceTraceGeometry>(impl_->geometry.source(), std::move(queries));
}
HydroInterfaceLoadGeometry::~HydroInterfaceLoadGeometry() = default;
HydroInterfaceLoadGeometry::HydroInterfaceLoadGeometry(HydroInterfaceLoadGeometry&&) noexcept = default;
HydroInterfaceLoadGeometry& HydroInterfaceLoadGeometry::operator=(HydroInterfaceLoadGeometry&&) noexcept = default;
const HydroInterfaceLoadGeometry::Impl& HydroInterfaceLoadGeometry::checked() const
{
    require(static_cast<bool>(impl_), "moved-from interface load geometry"); return *impl_;
}
const InterfaceGeometry& HydroInterfaceLoadGeometry::geometry() const { return checked().geometry; }
const LevelSetHermiteGeometry::Snapshot& HydroInterfaceLoadGeometry::snapshot() const { return checked().snapshot; }
const std::vector<double>& HydroInterfaceLoadGeometry::positiveFraction() const { return checked().source.fraction; }
const std::vector<double>& HydroInterfaceLoadGeometry::loadIndicator() const
{
    const auto& state = checked();
    return state.options.layer == Layer::GeometricPL ? state.source.fraction : state.source.indicator;
}
const HydroInterfaceLoadGeometry::Options& HydroInterfaceLoadGeometry::options() const { return checked().options; }
const std::vector<HydroInterfaceLoadGeometry::Face>& HydroInterfaceLoadGeometry::faces() const { return checked().source.faces; }
const std::shared_ptr<const InterfaceTraceGeometry>& HydroInterfaceLoadGeometry::traceGeometry() const
{ checked(); return traceGeometry_; }
void HydroInterfaceLoadGeometry::validateSource(const QuadtreeMesh& mesh, const VertField<double>& phi) const
{
    checked().geometry.validateSource(mesh, phi);
}
