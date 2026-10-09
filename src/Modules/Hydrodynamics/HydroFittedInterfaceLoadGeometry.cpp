#include "Modules/Hydrodynamics/HydroFittedInterfaceLoadGeometry.h"
#include "Modules/Hydrodynamics/HydroFittedDomain.h"
#include "Modules/Hydrodynamics/HydroNativePressureStencil.h"
#include "Core/Geometry/LevelSetInterfacePanels.h"
#include "Core/Geometry/LevelSetPhaseFraction.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/UniformLeafGrid.h"

#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>

namespace Quantum::Hydrodynamics
{
namespace
{
using Snapshot = LevelSetHermiteGeometry::Snapshot;
using Face = HydroFittedDomain::NativeFace;
using BoundaryType = HydroFittedDomain::BoundaryType;

void require(bool condition, const char *message)
{
    if (!condition)
        throw std::invalid_argument(std::string("[HydroFittedInterfaceLoadGeometry] ") + message);
}

template <class Query>
auto checkedQuery(const char *kind, int index, Vector2D point, Query query) -> decltype(query())
{
    try
    {
        return query();
    }
    catch (const std::exception &error)
    {
        std::ostringstream message;
        message << std::setprecision(17) << "[HydroFittedInterfaceLoadGeometry] " << kind << " " << index
                << " at (" << point.x << ',' << point.y << "): " << error.what();
        throw std::runtime_error(message.str());
    }
}

double checked(long double value, const char *message)
{
    const double result = static_cast<double>(value);
    require(std::isfinite(value) && std::isfinite(result) && (value == 0.0L || result != 0.0), message);
    return result;
}

bool pressure(BoundaryType type)
{
    return type == BoundaryType::PressureReservoir || type == BoundaryType::PressureOutlet;
}

// The boundary gradient is evaluated at the ACTUAL native midpoint, from its
// unique active incident cell. Match the retained scalar's cell-origin map,
// including its exact local==1 convention; do not use (p-min)/h-index, clamp
// a coordinate, or recover a near-wall point from an unrelated global jet.
Vector2D wallLocalPoint(const Snapshot &snapshot, const std::array<int, 2> &cell, const Face &face)
{
    const auto lower = snapshot.minCoord(), upper = snapshot.maxCoord(), h = snapshot.cellSize();
    const auto component = [&](int index, int count, double minimum, double maximum, double spacing,
                               double point) {
        const double origin = minimum + index * spacing;
        const double end = index + 1 == count ? maximum : minimum + (index + 1) * spacing;
        if (point == origin)
            return 0.0;
        if (point == end)
            return 1.0;
        const double local = static_cast<double>((static_cast<long double>(point) - origin) / spacing);
        require(std::isfinite(local) && local > 0.0 && local < 1.0,
                "actual wall midpoint is outside its active scalar cell");
        require(origin + spacing * local == point,
                "actual wall midpoint failed its same-cell physical round trip");
        return local;
    };
    const Vector2D local{component(cell[0], snapshot.nx(), lower.x, upper.x, h.x, face.center.x),
                         component(cell[1], snapshot.ny(), lower.y, upper.y, h.y, face.center.y)};
    require((face.horizontal ? local.y : local.x) == (face.negativeCell >= 0 ? 1.0 : 0.0),
            "actual wall does not coincide with its active scalar cell boundary");
    return local;
}

double wallIndicatorGradient(const Snapshot::ScalarJet2 &jet, bool horizontal, double halfWidth)
{
    // H is exactly constant outside the compact layer. Inside, evaluate the
    // derivative of the EXISTING quintic H without subtracting 1-s*s close
    // to an endpoint. This derivative adds no second interface thickness.
    if (jet.value <= -halfWidth || jet.value >= halfWidth)
        return 0.0;
    const long double tail = (static_cast<long double>(halfWidth) - std::abs(jet.value)) / halfWidth;
    const long double squareFactor = tail * (2.0L - tail);
    const long double derivative = (15.0L / 16.0L) * squareFactor * squareFactor / halfWidth;
    require(std::isfinite(derivative) && derivative > 0.0L,
            "compact wall indicator derivative overflow or underflow");
    const double component = horizontal ? jet.dy : jet.dx;
    const long double gradient = derivative * component;
    require(component == 0.0 || gradient != 0.0L, "nonzero wall indicator gradient underflow");
    return checked(gradient, "nonfinite or underflowed wall indicator gradient");
}

supCoeff nativeTransverseStencil(const QuadtreeMesh &mesh, const HydroFittedDomain &domain,
                                 const Face &target, double &maxAffineError)
{
    // Use the existing linear fit, on the fixed REAL face incidence of the
    // two active owners. No virtual uniform grid or solid scalar is sampled.
    supCoeff result;
    Vector2D scale{0, 0};
    for (int owner : {target.negativeCell, target.positiveCell})
    {
        require(owner >= 0 && domain.cellRowByIndex()[owner] >= 0,
                "native transverse interpolation attempted to consume a solid cell");
        const auto &cell = mesh.cells(mesh.cellIndexToLabel[owner]);
        scale.x = std::max(scale.x, cell.size().x);
        scale.y = std::max(scale.y, cell.size().y);
        for (int side = 0; side < 2; ++side)
            for (int part = 0; part < 2; ++part)
            {
                const int label =
                    target.horizontal ? cell.cellToVEdge[side][part] : cell.cellToHEdge[side][part];
                if (label < 0)
                    continue;
                const int row = target.horizontal ? mesh.vEdges(label).index()
                                                  : mesh.vEdgeIndexNum() + mesh.hEdges(label).index();
                require(row >= 0 && row < static_cast<int>(domain.faces().size()),
                        "native transverse interpolation has an invalid face index");
                const auto &face = domain.faces()[row];
                require(face.negativeCell == owner || face.positiveCell == owner,
                        "native transverse source is not incident on its active owner");
                result.ind.push_back(row);
            }
    }
    std::sort(result.ind.begin(), result.ind.end());
    result.ind.erase(std::unique(result.ind.begin(), result.ind.end()), result.ind.end());
    require(result.ind.size() >= 3 && scale.x > 0.0 && scale.y > 0.0,
            "native transverse fit lacks real support");
    std::vector<Vector2D> points;
    for (int row : result.ind)
    {
        const auto point = domain.faces()[row].center;
        points.push_back({(point.x - target.center.x) / scale.x, (point.y - target.center.y) / scale.y});
    }
    // Check rank before the legacy least-squares inverse, which still has an
    // exit(0) branch for insufficient support. A missing closure must fail.
    long double maxArea = 0.0L;
    for (std::size_t i = 1; i < points.size(); ++i)
        for (std::size_t j = i + 1; j < points.size(); ++j)
            maxArea = std::max(maxArea, std::abs((static_cast<long double>(points[i].x) - points[0].x) *
                                                     (points[j].y - points[0].y) -
                                                 (static_cast<long double>(points[j].x) - points[0].x) *
                                                     (points[i].y - points[0].y)));
    require(maxArea > 0.0L, "singular native transverse stencil");
    result.coe = Interpolation::leastSquareInterpolationCoeff(points, {0, 0});
    require(result.coe.size() == points.size(), "incomplete native transverse weights");
    long double one = 0, x = 0, y = 0;
    for (std::size_t i = 0; i < points.size(); ++i)
    {
        require(std::isfinite(result.coe[i]), "nonfinite native transverse coefficient");
        one += result.coe[i];
        x += static_cast<long double>(result.coe[i]) * points[i].x;
        y += static_cast<long double>(result.coe[i]) * points[i].y;
    }
    const double error = static_cast<double>(std::max({std::abs(one - 1), std::abs(x), std::abs(y)}));
    require(error <= 1e-10, "native transverse interpolation lost affine consistency");
    maxAffineError = std::max(maxAffineError, error);
    return result;
}
} // namespace

HydroFittedInterfaceLoadGeometry HydroFittedInterfaceLoadGeometry::build(
    const QuadtreeMesh &mesh, const HydroFittedDomain &domain, std::shared_ptr<const Snapshot> snapshot,
    const LevelSetInterfacePanels &retainedPanels, double halfWidth,
    FittedTransverseGradientClosure transverseClosure, bool nativeLeaves)
{
    require(snapshot && std::isfinite(halfWidth) && halfWidth > 0.0,
            "a retained scalar and positive finite load half-width are required");
    require(transverseClosure == FittedTransverseGradientClosure::ActiveHermiteWallGradient,
            "unknown fitted transverse-gradient closure");
    require(domain.axisymmetric() && snapshot->radialAxisSymmetry() && snapshot->restrictedDomain(),
            "fitted interface loads currently require an explicit active RZ scalar");
    require(retainedPanels.usesSource(*snapshot) && retainedPanels.axisymmetric() &&
                !retainedPanels.closed() && !retainedPanels.panels().empty(),
            "retained open RZ panels must use the identical shared scalar");
    require(snapshot->usesNativeLeaves() == nativeLeaves && retainedPanels.usesNativeLeaves() == nativeLeaves,
            "scalar and panels must match the explicit native/same-level load entry");
    domain.validateMeshIdentity(mesh);
    const int cellCount = mesh.cellIndexNum(), verticalCount = mesh.vEdgeIndexNum();
    const int nativeCount = verticalCount + mesh.hEdgeIndexNum();
    require(domain.faces().size() == static_cast<std::size_t>(nativeCount) &&
                domain.cellRowByIndex().size() == static_cast<std::size_t>(cellCount) &&
                retainedPanels.cellToPanel().size() == static_cast<std::size_t>(cellCount),
            "incomplete actual native domain or panel layout");
    std::vector<std::array<int, 2>> logicalCell(cellCount);
    if (nativeLeaves)
    {
        require(snapshot->nativeCellCount() == static_cast<std::size_t>(cellCount),
                "native scalar and mesh cell counts disagree");
        for (int k = 0; k < cellCount; ++k)
        {
            const auto &cell = mesh.cells(mesh.cellIndexToLabel[k]);
            const auto bounds = snapshot->nativeCellBounds(k);
            const auto low = mesh.verts(cell.cellToVert[0]).coord(),
                       high = mesh.verts(cell.cellToVert[2]).coord();
            require(bounds[0].x == low.x && bounds[0].y == low.y && bounds[1].x == high.x &&
                        bounds[1].y == high.y,
                    "native scalar and actual leaf bounds disagree");
            require(snapshot->nativeCellActive(k) == (domain.cellRowByIndex()[k] >= 0),
                    "scalar and fitted pressure domains must have identical activity");
        }
    }
    else
    {
        const UniformLeafGrid grid(mesh, CartesianGridScope::UniformLeaves);
        require(
            snapshot->nx() == grid.nx() && snapshot->ny() == grid.ny() &&
                snapshot->minCoord().x == grid.minCoord().x && snapshot->minCoord().y == grid.minCoord().y &&
                snapshot->maxCoord().x == grid.maxCoord().x && snapshot->maxCoord().y == grid.maxCoord().y &&
                snapshot->cellSize().x == grid.spacing().x && snapshot->cellSize().y == grid.spacing().y,
            "snapshot and fitted mesh layouts disagree");
        for (int j = 0; j < grid.ny(); ++j)
            for (int i = 0; i < grid.nx(); ++i)
            {
                const int k = mesh.cells(grid.cells()[static_cast<std::size_t>(j) * grid.nx() + i]).index();
                logicalCell[k] = {{i, j}};
                require(snapshot->cellActive(i, j) == (domain.cellRowByIndex()[k] >= 0),
                        "scalar and fitted pressure domains must have identical activity");
            }
    }
    for (std::size_t p = 0; p < retainedPanels.panels().size(); ++p)
    {
        const int cell = retainedPanels.panels()[p].cellIndex;
        require(cell >= 0 && cell < cellCount && domain.cellRowByIndex()[cell] >= 0 &&
                    retainedPanels.cellToPanel()[cell] == static_cast<int>(p),
                "retained interface graph does not match the actual active native cell");
        if (nativeLeaves)
            require(retainedPanels.nativeChain().graphs()[p].nativeCellIndex == cell,
                    "retained native graph and panel cell identities disagree");
        else
        {
            const auto &graph = retainedPanels.panelGraph(static_cast<int>(p));
            require(graph.cellI == logicalCell[cell][0] && graph.cellJ == logicalCell[cell][1],
                    "retained interface graph does not match its logical cell");
        }
    }

    HydroFittedInterfaceLoadGeometry result;
    result.halfWidth = halfWidth;
    result.indicator.assign(cellCount, std::numeric_limits<double>::quiet_NaN());
    Eigen::VectorXd h(domain.cellIndexByRow().size());
    for (std::size_t row = 0; row < domain.cellIndexByRow().size(); ++row)
    {
        const int cell = domain.cellIndexByRow()[row];
        const auto point = mesh.cells(mesh.cellIndexToLabel[cell]).centerCoord();
        // Preserve the original fitted capillary indicator evaluation and
        // pressure-load arithmetic, rather than switching the primary H to
        // a new local evaluation while adding transverse support.
        h[static_cast<Eigen::Index>(row)] = result.indicator[cell] =
            LevelSetPhaseFraction::regularizedPositiveIndicator(
                checkedQuery("cell indicator", cell, point, [&] { return snapshot->value(point); }),
                halfWidth);
    }
    const Eigen::VectorXd pressureLoad = domain.pressureGradientLoad() * h;
    require(pressureLoad.allFinite(), "nonfinite fitted C(H)");
    std::vector<double> sourceGradient(nativeCount, std::numeric_limits<double>::quiet_NaN());
    std::vector<unsigned char> pressurePlateau(nativeCount, 0);
    double identityScale = 0.0;
    std::map<std::pair<int, int>, supCoeff> centers;
    const auto childCenter = [&](const Cell &cell, int child) -> const supCoeff & {
        const auto key = std::make_pair(cell.index(), child);
        const auto found = centers.find(key);
        if (found != centers.end())
            return found->second;
        auto stencil = HydroNativePressureStencil::childCenter(mesh, cell, child);
        for (int column : stencil.ind)
            require(column >= 0 && column < cellCount && domain.cellRowByIndex()[column] >= 0,
                    "native load pressure interpolation consumes an inactive column");
        return centers.emplace(key, std::move(stencil)).first->second;
    };
    for (const auto &face : domain.faces())
    {
        if (face.velocityRow < 0)
            continue;
        const int low = face.negativeCell, high = face.positiveCell;
        if (low < 0 || high < 0)
        {
            const int owner = low >= 0 ? low : high;
            require(owner >= 0 && pressure(face.boundary.type), "unsupported free boundary force row");
            const auto center = mesh.cells(mesh.cellIndexToLabel[owner]).centerCoord();
            require(checkedQuery("pressure port owner", owner, center,
                                 [&] { return snapshot->closedCellOutsideBand(center, halfWidth); }) &&
                        (result.indicator[owner] == 0.0 || result.indicator[owner] == 1.0),
                    "interface load band reaches a pressure port; boundary force/phase closure is required");
            pressurePlateau[face.nativeRow] = 1;
            continue; // A certified constant H, not an arbitrary port force.
        }
        const double integrated = pressureLoad[face.velocityRow];
        double independent = face.area * (result.indicator[high] - result.indicator[low]);
        if (nativeLeaves &&
            mesh.cells(mesh.cellIndexToLabel[low]).level() != mesh.cells(mesh.cellIndexToLabel[high]).level())
        {
            long double sum = 0;
            const auto &edge = face.horizontal ? mesh.hEdges(face.edgeLabel) : mesh.vEdges(face.edgeLabel);
            HydroNativePressureStencil::visitGradient(
                mesh, edge, face.horizontal, childCenter, [&](int column, double coefficient) {
                    sum += static_cast<long double>(face.dualVolume * coefficient) * result.indicator[column];
                });
            independent = checked(sum, "unrepresentable independent native pressure load");
        }
        result.maxGradientIdentityError =
            std::max(result.maxGradientIdentityError, std::abs(integrated - independent));
        identityScale = std::max({identityScale, std::abs(integrated), std::abs(independent)});
        const double gradient = integrated / face.dualVolume;
        require(std::isfinite(gradient) && (integrated == 0.0 || gradient != 0.0),
                "unrepresentable fitted pressure-gradient density");
        sourceGradient[face.nativeRow] = gradient;
    }
    require(result.maxGradientIdentityError <= 512.0 * std::numeric_limits<double>::epsilon() * identityScale,
            "actual C(H) disagrees with independently assembled native pressure-face load");

    const auto source = [&](int native) {
        require(native >= 0 && native < nativeCount, "transverse source has no actual native face");
        double &value = sourceGradient[native];
        if (std::isfinite(value))
            return value;
        const auto &face = domain.faces()[native];
        const int owner = face.negativeCell >= 0 ? face.negativeCell : face.positiveCell;
        require(owner >= 0 && ((face.negativeCell >= 0) != (face.positiveCell >= 0)),
                "missing transverse source is not a unique active-side boundary face");
        if (pressure(face.boundary.type))
        {
            require(pressurePlateau[native] != 0,
                    "pressure source lacks its closed-cell load-band certificate");
            ++result.pressurePlateauSources;
            value = 0.0;
            return value;
        }
        require(face.velocityRow < 0 && (face.boundary.type == BoundaryType::NoSlipWall ||
                                         face.boundary.type == BoundaryType::Symmetry),
                "unsupported transverse scalar boundary closure");
        return checkedQuery("transverse boundary source", native, face.center, [&] {
            const auto ij = logicalCell[owner];
            Vector2D local;
            Snapshot::ScalarJet2 jet;
            if (nativeLeaves)
            {
                const auto bounds = snapshot->nativeCellBounds(owner);
                const auto side = bounds[face.negativeCell >= 0 ? 1 : 0];
                require((face.horizontal ? face.center.y == side.y : face.center.x == side.x),
                        "actual wall does not coincide with its native scalar leaf boundary");
                jet = snapshot->scalarJetInNativeCell(owner, face.center);
            }
            else
            {
                local = wallLocalPoint(*snapshot, ij, face);
                jet = snapshot->scalarJetAtCellLocal(ij[0], ij[1], local);
            }
            if (face.boundary.type == BoundaryType::Symmetry)
            {
                require(!face.horizontal && face.center.x == 0.0 && snapshot->minCoord().x == 0.0 &&
                            (nativeLeaves ? snapshot->nativeCellBounds(owner)[0].x == 0.0
                                          : ij[0] == 0 && local.x == 0.0) &&
                            jet.dx == 0.0,
                        "only the true even-scalar radial axis has a zero symmetry source");
                ++result.axisGradientSources;
                value = 0.0;
            }
            else
            {
                value = wallIndicatorGradient(jet, face.horizontal, halfWidth);
                ++result.wallGradientSources;
                if (value != 0.0)
                    ++result.nonzeroWallGradientSources;
                result.maxAbsWallGradient = std::max(result.maxAbsWallGradient, std::abs(value));
            }
            return value;
        });
    };
    const auto transverseFace = [&](int owner, bool horizontal, int side) {
        require(owner >= 0 && owner < cellCount && domain.cellRowByIndex()[owner] >= 0,
                "transverse interpolation attempted to consume a solid cell");
        const auto &cell = mesh.cells(mesh.cellIndexToLabel[owner]);
        const int label = horizontal ? cell.cellToHEdge[side][0] : cell.cellToVEdge[side][0];
        require(label >= 0, "transverse interpolation lacks an actual native cell face");
        return horizontal ? verticalCount + mesh.hEdges(label).index() : mesh.vEdges(label).index();
    };
    for (const auto &native : domain.faces())
    {
        if (native.velocityRow < 0 || pressurePlateau[native.nativeRow])
            continue;
        const int low = native.negativeCell, high = native.positiveCell;
        require(low >= 0 && high >= 0, "unclassified free interface-load row");
        // Preserve the original full-domain four-face SUM ORDER as well as
        // its 1/4 average: y row W/E below, W/E above; x row S left/right,
        // N left/right. Actual cell incidence, not velocity constraints,
        // selects every transverse scalar-gradient source.
        double transverse = 0;
        if (nativeLeaves)
        {
            const auto stencil =
                nativeTransverseStencil(mesh, domain, native, result.maxTransverseAffineError);
            long double sum = 0;
            for (std::size_t i = 0; i < stencil.ind.size(); ++i)
                sum += static_cast<long double>(stencil.coe[i]) * source(stencil.ind[i]);
            transverse = checked(sum, "nonfinite or underflowed native transverse gradient");
            ++result.nativeTransverseFits;
        }
        else
        {
            const std::array<int, 4> neighbors =
                native.horizontal
                    ? std::array<int, 4>{{transverseFace(low, false, 0), transverseFace(low, false, 1),
                                          transverseFace(high, false, 0), transverseFace(high, false, 1)}}
                    : std::array<int, 4>{{transverseFace(low, true, 0), transverseFace(high, true, 0),
                                          transverseFace(low, true, 1), transverseFace(high, true, 1)}};
            long double sum = 0.0L;
            for (int neighbor : neighbors)
                sum += source(neighbor);
            const long double average = sum / 4.0L;
            require(std::isfinite(sum) && (sum == 0.0L || average != 0.0L),
                    "transverse pressure-gradient average overflow or underflow");
            transverse = checked(average, "transverse pressure-gradient average overflow or underflow");
        }
        const double gradient = sourceGradient[native.nativeRow];
        require(std::isfinite(gradient), "free interface-load row lacks its primary fitted gradient");
        if (gradient == 0.0 && transverse == 0.0)
            continue;
        if (gradient == 0.0)
            ++result.transverseOnlyFaces;
        result.maxAbsTransverseGradient = std::max(result.maxAbsTransverseGradient, std::abs(transverse));
        HydroInterfaceLoadGeometry::Face face;
        face.horizontal = native.horizontal;
        face.edgeIndex = native.horizontal ? native.nativeRow - verticalCount : native.nativeRow;
        face.lowCellIndex = low;
        face.highCellIndex = high;
        face.point = native.center;
        face.gradient = gradient;
        face.transverseGradient = transverse;
        face.dualMeasure = native.dualVolume;
        face.dualArea = native.normalSpacing * native.tangentSpacing;
        result.faces_.push_back(face);
    }

    // The pressure-compatible support and its ordering are now complete.
    // Projection/curvature at distinct faces share only immutable geometry;
    // each task writes one separate result slot, not a field reduction.
    const auto populateFrame = [&](std::size_t index) {
        auto& face = result.faces_[index];
        const int nativeRow = face.horizontal ? verticalCount + face.edgeIndex : face.edgeIndex;
        const auto projection = checkedQuery("native load-face projection", nativeRow, face.point,
                                             [&] { return retainedPanels.projectedFrame(face.point); });
        face.frame = projection.frame;
        face.boundaryEndpointProjection = projection.endpoint;
        face.endpointKktResidual = projection.kktResidual;
        require(face.frame.hasCellLocalPoint != face.frame.hasCellChartPoint && face.frame.point.x >= 0.0 &&
                    std::isfinite(face.frame.curvature),
                "load projection must retain exactly one valid local or boundary-chart RZ frame");
        const bool axis = face.frame.point.x == 0.0;
        require(!axis || face.frame.normal.x == 0.0, "axis curvature requires the actual even-scalar pole");
        face.surfaceCurvature =
            face.frame.curvature + (axis ? face.frame.curvature : face.frame.normal.x / face.frame.point.x);
        require(std::isfinite(face.surfaceCurvature), "nonfinite RZ surface curvature");
        face.curvatureTangentialDerivative =
            checkedQuery("load-frame curvature derivative", nativeRow, face.frame.point, [&] {
                if (nativeLeaves)
                    return snapshot->curvatureTangentialDerivativeAtNativeCellChart(
                        face.frame.nativeCellIndex, face.frame.cellChartPoint);
                return face.frame.hasCellChartPoint
                           ? snapshot->curvatureTangentialDerivativeAtCellChart(
                                 face.frame.cellI, face.frame.cellJ, face.frame.cellChartPoint)
                           : snapshot->curvatureTangentialDerivativeAtCellLocal(
                                 face.frame.cellI, face.frame.cellJ, face.frame.cellLocalPoint);
            });
        require(std::isfinite(face.curvatureTangentialDerivative),
                "nonfinite retained local curvature derivative");
    };
    // Match the other native geometry kernels: validate a real query on the
    // caller before starting workers, capture failures, and rethrow in face
    // order after joining. No incomplete geometry escapes this function.
    if (!result.faces_.empty()) populateFrame(0);
    std::vector<std::exception_ptr> failures(result.faces_.size());
#pragma omp parallel for schedule(dynamic, 8) if(nativeLeaves && result.faces_.size() > 32)
    for (std::ptrdiff_t index = 1; index < static_cast<std::ptrdiff_t>(result.faces_.size()); ++index)
    {
        try { populateFrame(static_cast<std::size_t>(index)); }
        catch (...) { failures[index] = std::current_exception(); }
    }
    for (std::size_t index = 0; index < result.faces_.size(); ++index)
    {
        if (failures[index]) std::rethrow_exception(failures[index]);
        result.endpointProjections += result.faces_[index].boundaryEndpointProjection ? 1 : 0;
    }
    std::vector<InterfaceTraceGeometry::Point> points;
    points.reserve(result.faces_.size());
    for (const auto& face : result.faces_)
        points.push_back({face.frame, face.curvatureTangentialDerivative, face.endpointKktResidual});
    result.traceGeometry_ = std::make_shared<const InterfaceTraceGeometry>(snapshot, std::move(points));
    result.meshIdentity_ = domain.geometry().meshIdentity();
    return result;
}

void HydroFittedInterfaceLoadGeometry::validate(const QuadtreeMesh& mesh, const HydroFittedDomain& domain,
    const std::shared_ptr<const InterfaceTraceGeometry>& traceGeometry) const
{
    require(traceGeometry_ && traceGeometry == traceGeometry_ && meshIdentity_
        && meshIdentity_ == domain.geometry().meshIdentity(), "load geometry, trace sampling and receiving domain must share their actual source");
    domain.validateMeshIdentity(mesh);
    require(traceGeometry_->points().size() == faces_.size(), "incomplete retained load sampling");
}

HydroFittedInterfaceLoadGeometry buildFittedInterfaceLoadGeometry(
    const QuadtreeMesh &mesh, const HydroFittedDomain &domain, std::shared_ptr<const Snapshot> snapshot,
    const LevelSetInterfacePanels &panels, double halfWidth, FittedTransverseGradientClosure closure)
{
    return HydroFittedInterfaceLoadGeometry::build(mesh, domain, std::move(snapshot), panels, halfWidth, closure, false);
}

HydroFittedInterfaceLoadGeometry buildNativeFittedInterfaceLoadGeometry(
    const QuadtreeMesh &mesh, const HydroFittedDomain &domain, std::shared_ptr<const Snapshot> snapshot,
    const LevelSetInterfacePanels &panels, double halfWidth, FittedTransverseGradientClosure closure)
{
    return HydroFittedInterfaceLoadGeometry::build(mesh, domain, std::move(snapshot), panels, halfWidth, closure, true);
}
} // namespace Quantum::Hydrodynamics
