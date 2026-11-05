#include "Modules/Electric/SurfaceChargeField.h"
#include "Core/Field/StaggeredVelocitySampler.h"
#include "Core/Geometry/CellChartPoint.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/CartesianMeasure.h"
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <exception>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace LAMP::Electric
{
namespace
{
void require(bool condition, const char* reason)
{ if (!condition) throw std::runtime_error(reason); }
bool same(Vector2D a, Vector2D b) { return a.x == b.x && a.y == b.y; }
std::vector<double> measures(const LevelSetInterfacePanels& geometry)
{
    std::vector<double> result;
    for (const auto& panel : geometry.panels())
        result.push_back(geometry.axisymmetric() ? panel.physicalMeasure : panel.measure);
    return result;
}
}

std::array<double,3> SurfaceChargeDensityStencil::densityWeightsAtArclength(double s) const
{
    require(std::isfinite(panelMeasure) && panelMeasure>0 && std::isfinite(scale) && scale>0
        && std::isfinite(origin) && origin>=0 && origin<=panelMeasure,
        "invalid surface-charge density stencil coordinates");
    require(std::isfinite(s) && s>=0 && s<=panelMeasure,
        "surface-charge density stencil arclength lies outside its panel");
    const double u=(s-origin)/scale;
    require(std::isfinite(u),"surface-charge density stencil coordinate overflow");
    std::array<double,3> result;
    for(std::size_t col=0;col<3;++col) {
        const double a=polynomialWeights[col],b=polynomialWeights[3+col],c=polynomialWeights[6+col];
        require(std::isfinite(a) && std::isfinite(b) && std::isfinite(c),
            "nonfinite surface-charge density stencil coefficient");
        if(evenAtAxis) {
            const double square=u*u;result[col]=a+square*(b+square*c);
        } else result[col]=a+u*(b+u*c);
        require(std::isfinite(result[col]),"surface-charge density stencil value overflow");
    }
    return result;
}

SurfaceChargeField::SurfaceChargeField(const QuadtreeMesh& mesh,
    const LevelSetHermiteGeometry::Snapshot& snapshot, double initialDensity, bool axisymmetric,
    std::vector<Vector2D> boundaryContacts)
    : panels_(mesh, snapshot, axisymmetric, std::move(boundaryContacts)), charge_(measures(panels_),
          std::vector<double>(panels_.panels().size(), initialDensity), panels_.closed())
{ initialize(mesh); }

SurfaceChargeField::SurfaceChargeField(const QuadtreeMesh& mesh,
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> snapshot, double initialDensity,
    SurfaceChargeGeometryOptions options)
    : source_(std::move(snapshot)), panels_(mesh, source_, options.axisymmetric,
          std::move(options.boundaryContacts), options.scope), charge_(measures(panels_),
          std::vector<double>(panels_.panels().size(), initialDensity), panels_.closed()), scope_(options.scope)
{ initialize(mesh); }

SurfaceChargeField::SurfaceChargeField(const QuadtreeMesh &mesh,
                                       LevelSetInterfacePanels::NativeGraphChain chain, double initialDensity)
    : SurfaceChargeField(mesh,LevelSetInterfacePanels(std::move(chain)),initialDensity)
{}

SurfaceChargeField::SurfaceChargeField(const QuadtreeMesh &mesh,
                                       LevelSetInterfacePanels panels, double initialDensity)
    : source_(panels.nativeChain().source()), panels_(std::move(panels)),
      charge_(measures(panels_), std::vector<double>(panels_.panels().size(), initialDensity),
              panels_.closed())
{
    initialize(mesh);
    validateNativeRemapGrid(mesh, *source_);
}

void SurfaceChargeField::validateNativeRemapGrid(const QuadtreeMesh &mesh,
                                                 const LevelSetHermiteGeometry::Snapshot &next) const
{
    require(source_ && source_->usesNativeLeaves() && next.usesNativeLeaves(),
            "native material remap requires two native scalar sources");
    require(source_->nativeCellCount() == next.nativeCellCount() &&
                source_->nativeCellCount() == static_cast<std::size_t>(mesh.cellIndexNum()) &&
                retainedCellLabels_.size() == static_cast<std::size_t>(mesh.cellIndexNum()) &&
                retainedVertexLabels_.size() == static_cast<std::size_t>(mesh.vertIndexNum()) &&
                same(source_->minCoord(), next.minCoord()) && same(source_->maxCoord(), next.maxCoord()) &&
                same(source_->minCoord(), mesh.minCoord) && same(source_->maxCoord(), mesh.maxCoord) &&
                source_->radialAxisSymmetry() == next.radialAxisSymmetry(),
            "native material remap changed its fixed grid/domain");
    for (int c = 0; c < mesh.cellIndexNum(); ++c)
    {
        require(retainedCellLabels_[c] == mesh.cellIndexToLabel[c],
                "native material remap changed cell labels");
        const auto a = source_->nativeCellBounds(c), b = next.nativeCellBounds(c);
        const auto &cell = mesh.cells(mesh.cellIndexToLabel[c]);
        require(same(a[0], b[0]) && same(a[1], b[1]) && same(a[0], mesh.verts(cell.cellToVert[0]).coord()) &&
                    same(a[1], mesh.verts(cell.cellToVert[2]).coord()) &&
                    source_->nativeCellActive(c) == next.nativeCellActive(c),
                "native material remap changed leaf bounds or the active/solid mask");
    }
    for (int v = 0; v < mesh.vertIndexNum(); ++v)
        require(retainedVertexLabels_[v] == mesh.vertIndexToLabel[v],
                "native material remap changed vertex labels");
}

double SurfaceChargeField::nativeFrameFraction(std::size_t panel,
                                               const LevelSetHermiteGeometry::Frame &frame) const
{
    return panels_.nativeChain().graphFraction(static_cast<int>(panel), frame.cellChartPoint);
}

void SurfaceChargeField::initialize(const QuadtreeMesh &mesh)
{
    const auto &panels = panels_.panels();
    const std::size_t count = panels.size();
    coordinateTolerance_ =
        512.0 * std::numeric_limits<double>::epsilon() *
        std::max({mesh.cellSize.x, mesh.cellSize.y, std::abs(mesh.minCoord.x), std::abs(mesh.minCoord.y),
                  std::abs(mesh.maxCoord.x), std::abs(mesh.maxCoord.y)});
    if (source_)
    {
        if (panels_.usesNativeLeaves())
        {
            for (int c = 0; c < mesh.cellIndexNum(); ++c)
                retainedCellLabels_.push_back(mesh.cellIndexToLabel[c]);
            for (int v = 0; v < mesh.vertIndexNum(); ++v)
                retainedVertexLabels_.push_back(mesh.vertIndexToLabel[v]);
        }
        else
        {
            const UniformLeafGrid grid(mesh, scope_);
            retainedCellLabels_ = grid.cells();
            retainedVertexLabels_ = grid.vertices();
        }
    }
    mappings_.resize(count);
    for (std::size_t i = 0; i < count; ++i)
    {
        const auto &panel = panels[i];
        auto &mapping = mappings_[i];
        double reciprocalCondition = std::numeric_limits<double>::quiet_NaN();
        double inverseResidual = std::numeric_limits<double>::quiet_NaN();
        mapping.neighbor = {{(i + count - 1) % count, i, (i + 1) % count}};
        if (!panels_.closed())
        {
            const std::size_t first = i == 0 ? 0 : i + 1 == count ? count - 3 : i - 1;
            mapping.neighbor = {{first, first + 1, first + 2}};
        }
        try
        {
            const bool axisBegin =
                i == 0 && panels_.beginEndpointKind() == LevelSetInterfacePanels::EndpointKind::Axis;
            const bool axisEnd =
                i + 1 == count && panels_.endEndpointKind() == LevelSetInterfacePanels::EndpointKind::Axis;
            const auto &cell = mesh.cells(mesh.cellIndexToLabel[panel.cellIndex]);
            mapping.low = cell.centerCoord() - cell.size() * .5;
            mapping.high = cell.centerCoord() + cell.size() * .5;
            mapping.coordinateBegin = panel.coordinateIsX ? panel.begin.x : panel.begin.y;
            mapping.coordinateSpan =
                (panel.coordinateIsX ? panel.end.x : panel.end.y) - mapping.coordinateBegin;
            double coordinateScale = 1.0;
            if (panels_.usesNativeLeaves())
            {
                const auto &chain = panels_.nativeChain();
                const auto bounds = source_->nativeCellBounds(panel.cellIndex);
                mapping.low = bounds[0];
                mapping.high = bounds[1];
                mapping.coordinateBegin = 0;
                mapping.coordinateSpan = 1;
                coordinateScale = chain.physicalCoordinateSpan(static_cast<int>(i));
            }
            else if (source_)
            {
                const auto &graph = panels_.panelGraph(static_cast<int>(i));
                require(panelByCell_.emplace(std::make_pair(graph.cellI, graph.cellJ), i).second,
                        "multiple charge panels share one retained cell-local graph");
                mapping.coordinateBegin = graph.coordinateIsX ? graph.localBegin.x : graph.localBegin.y;
                mapping.coordinateSpan =
                    (graph.coordinateIsX ? graph.localEnd.x : graph.localEnd.y) - mapping.coordinateBegin;
                coordinateScale = graph.coordinateIsX ? graph.cellSize.x : graph.cellSize.y;
                require(std::isfinite(coordinateScale) && coordinateScale > 0.0,
                        "invalid retained charge graph physical coordinate scale");
            }
            require(std::isfinite(mapping.coordinateSpan) && mapping.coordinateSpan != 0.0,
                    "surface-charge panel lacks a resolved monotone coordinate");
            auto appendArc = [&](const LevelSetHermiteGeometry::Frame &savedFrame, double arclength,
                                 double physicalMeasure, int endpoint, double nativeFraction) {
                require(source_ || !savedFrame.hasCellChartPoint,
                        "legacy global charge quadrature cannot consume a boundary-chart Frame");
                const auto frame = source_ ? retainedFrame(i, savedFrame, endpoint) : savedFrame;
                const double tangentCoordinate = panel.coordinateIsX ? -frame.normal.y : frame.normal.x;
                const double derivative =
                    source_ ? static_cast<double>(std::abs(static_cast<long double>(mapping.coordinateSpan) *
                                                           coordinateScale / tangentCoordinate))
                            : std::abs(mapping.coordinateSpan / tangentCoordinate);
                double coordinate;
                if (panels_.usesNativeLeaves())
                    coordinate = nativeFraction;
                else if (source_ && frame.hasCellChartPoint)
                {
                    const auto &graph = panels_.panelGraph(static_cast<int>(i));
                    const double last = graph.coordinateIsX ? graph.localEnd.x : graph.localEnd.y;
                    coordinate = Geometry::cellChartIntervalFraction(
                        graph.coordinateIsX ? frame.cellChartPoint.x : frame.cellChartPoint.y,
                        mapping.coordinateBegin, last);
                }
                else
                {
                    const auto point = source_ ? frame.cellLocalPoint : frame.point;
                    const double u = panel.coordinateIsX ? point.x : point.y;
                    coordinate = (u - mapping.coordinateBegin) / mapping.coordinateSpan;
                }
                require(std::isfinite(derivative) && derivative > 0.0 && std::isfinite(coordinate),
                        "surface-charge arc mapping is singular");
                mapping.arc.push_back({coordinate, arclength, derivative});
                if (panels_.axisymmetric())
                {
                    const double physicalDerivative =
                        CartesianMeasure::radialWeight(frame.point.x) * derivative;
                    const bool axisEndpoint =
                        frame.point.x == 0.0 &&
                        ((axisBegin && arclength == 0.0) || (axisEnd && arclength == panel.measure));
                    require(std::isfinite(physicalDerivative) &&
                                (physicalDerivative > 0.0 || (axisEndpoint && physicalDerivative == 0.0)) &&
                                std::isfinite(physicalMeasure) && physicalMeasure >= 0.0,
                            "surface-charge physical measure mapping is singular");
                    mapping.physical.push_back({coordinate, physicalMeasure, physicalDerivative});
                }
            };
            appendArc(panel.beginFrame, 0.0, 0.0, -1, 0);
            for (const auto &node : panel.quadratureNodes)
                appendArc(node.frame, node.localArclength, node.localPhysicalMeasure, 0,
                          node.nativeGraphFraction);
            appendArc(panel.endFrame, panel.measure, panel.physicalMeasure, 1, 1);
            measureInterpolation_.rationalArcIntervals+=prepareMapping(mapping.arc, "surface-charge arc mapping is not strictly ordered",
                            "surface-charge arc interpolant is not monotone");
            if (panels_.axisymmetric())
                measureInterpolation_.rationalAreaIntervals+=prepareMapping(
                    mapping.physical, "surface-charge physical measure mapping is not strictly ordered",
                    "surface-charge physical measure interpolant is not monotone", axisBegin, axisEnd);

            mapping.evenAtAxis = axisBegin || axisEnd;
            mapping.polynomialOrigin = axisBegin ? 0.0 : axisEnd ? panel.measure : .5 * panel.measure;
            mapping.scale =
                std::max({panels[mapping.neighbor[0]].measure, panels[mapping.neighbor[1]].measure,
                          panels[mapping.neighbor[2]].measure});
            Eigen::Matrix3d moments;
            for (int row = 0; row < 3; ++row)
            {
                const auto &other = panels[mapping.neighbor[row]];
                double offset = row == 0 ? -other.measure : row == 1 ? 0.0 : panel.measure;
                if (!panels_.closed())
                {
                    offset = 0.0;
                    for (std::size_t j = mapping.neighbor[row]; j < i; ++j)
                        offset -= panels[j].measure;
                    for (std::size_t j = i; j < mapping.neighbor[row]; ++j)
                        offset += panels[j].measure;
                }
                long double first = 0.0L, second = 0.0L;
                for (const auto &node : other.quadratureNodes)
                {
                    const long double x =
                        (offset + node.localArclength - mapping.polynomialOrigin) / mapping.scale;
                    const double weight = panels_.axisymmetric() ? node.physicalWeight : node.weight;
                    if (mapping.evenAtAxis)
                    {
                        const long double square = x * x;
                        first += weight * square;
                        second += weight * square * square;
                    }
                    else
                    {
                        first += weight * x;
                        second += weight * x * x;
                    }
                }
                const double measure = panels_.axisymmetric() ? other.physicalMeasure : other.measure;
                moments(row, 0) = 1.0;
                moments(row, 1) = static_cast<double>(first / measure);
                moments(row, 2) = static_cast<double>(second / measure);
            }
            Eigen::FullPivLU<Eigen::Matrix3d> factor(moments);
            const bool invertible = moments.allFinite() && factor.isInvertible();
            if (invertible)
                reciprocalCondition = factor.rcond();
            require(invertible && reciprocalCondition > 1e-12,
                    "unresolved neighboring charge panels for conservative quadratic reconstruction");
            const Eigen::Matrix3d inverse = factor.inverse();
            if (inverse.allFinite())
                inverseResidual = (moments * inverse - Eigen::Matrix3d::Identity()).norm();
            require(inverse.allFinite() && inverseResidual < 1e-9,
                    "charge reconstruction moment inverse failed");
            for (int row = 0; row < 3; ++row)
                for (int col = 0; col < 3; ++col)
                    mapping.inverseMoments[3 * row + col] = inverse(row, col);
        }
        catch (const std::exception &error)
        {
            std::ostringstream context;
            context << std::setprecision(17) << "[SurfaceChargeField] panel=" << i
                << ", nativeCell=" << panel.cellIndex << ", neighbors=[";
            for (int k = 0; k < 3; ++k)
            {
                const auto index = mapping.neighbor[k];
                if (k) context << ';';
                context << index << ":cell=" << panels[index].cellIndex
                    << ",length=" << panels[index].measure << ",area=" << panels[index].physicalMeasure;
            }
            context << "], scale=" << mapping.scale << ", rcond=" << reciprocalCondition
                << ", inverseResidual=" << inverseResidual << ": " << error.what();
            throw std::runtime_error(context.str());
        }
    }
}

void SurfaceChargeField::validateLocalFrame(std::size_t panelIndex,
    const LevelSetHermiteGeometry::Frame& frame) const
{
    require(static_cast<bool>(source_), "cell-local charge query requires a retained scalar");
    if (panels_.usesNativeLeaves())
    {
        const int cell=panels_.panels().at(panelIndex).cellIndex;
        require(frame.hasCellChartPoint && !frame.hasCellLocalPoint && frame.cellI==-1 && frame.cellJ==-1
            && frame.nativeCellIndex==cell,"native charge frame lost its actual leaf/chart identity");
        Geometry::validateCellChart(frame.cellChartPoint);
        const auto actual=source_->jetAtNativeCellChart(cell,frame.cellChartPoint);
        require(same(frame.point,actual.point) && same(frame.normal,actual.normal)
            && frame.curvature==actual.curvature,"native charge frame differs from its retained scalar jet");
        const auto bounds=source_->nativeCellBounds(cell);
        const double h=std::min(bounds[1].x-bounds[0].x,bounds[1].y-bounds[0].y);
        require(std::isfinite(frame.levelSetResidual) && frame.levelSetResidual>=0 && frame.levelSetResidual<=1e-9*h
            && std::isfinite(actual.levelSetResidual) && actual.levelSetResidual<=1e-9*h,
            "native charge frame is off the retained interface");
        return;
    }
    const auto& graph = panels_.panelGraph(static_cast<int>(panelIndex));
    require(frame.hasCellLocalPoint != frame.hasCellChartPoint,
        "charge query requires exactly one retained local or boundary-chart provenance");
    require(frame.cellI == graph.cellI && frame.cellJ == graph.cellJ,
        "charge frame lost its retained cell identity");
    if (frame.hasCellChartPoint) Geometry::validateCellChart(frame.cellChartPoint);
    else
    {
        const auto local = frame.cellLocalPoint;
        require(std::isfinite(local.x) && std::isfinite(local.y) && local.x >= 0.0 && local.x <= 1.0
            && local.y >= 0.0 && local.y <= 1.0, "charge frame is outside its retained local cell");
    }
    // A cell id alone does not authorize a different one-sided corner jet.
    // Reevaluate the immutable source at the unchanged local coordinate.
    // Projection iteration/stationarity data are intentionally not compared:
    // they describe how the caller reached the point, not the scalar jet.
    const auto actual = frame.hasCellChartPoint
        ? source_->jetAtCellChart(graph.cellI,graph.cellJ,frame.cellChartPoint)
        : source_->jetAtCellLocal(graph.cellI,graph.cellJ,frame.cellLocalPoint);
    require(actual.hasCellLocalPoint == frame.hasCellLocalPoint && actual.hasCellChartPoint == frame.hasCellChartPoint
        && same(frame.point, actual.point) && same(frame.normal, actual.normal)
        && frame.curvature == actual.curvature,
        "charge frame does not match the retained scalar's local/chart jet/sector");
    const double h = std::min(graph.cellSize.x, graph.cellSize.y);
    require(std::isfinite(frame.levelSetResidual) && frame.levelSetResidual >= 0.0 && frame.levelSetResidual <= 1e-9*h
        && std::isfinite(actual.levelSetResidual) && actual.levelSetResidual >= 0.0 && actual.levelSetResidual <= 1e-9*h,
        "charge frame is off the retained interface");
}

LevelSetHermiteGeometry::Frame SurfaceChargeField::retainedFrame(std::size_t panelIndex,
    const LevelSetHermiteGeometry::Frame& frame, int endpoint) const
{
    if (panels_.usesNativeLeaves())
    {
        validateLocalFrame(panelIndex,frame);
        if (endpoint!=0)
        {
            const auto actual=panels_.nativeChain().frameAt(static_cast<int>(panelIndex),endpoint<0?0:1);
            require(Geometry::cellChartEqual(frame.cellChartPoint,actual.cellChartPoint),
                "native charge endpoint differs from its original root identity");
        }
        return frame;
    }
    const auto& graph = panels_.panelGraph(static_cast<int>(panelIndex));
    if (frame.hasCellLocalPoint || frame.hasCellChartPoint)
    {
        validateLocalFrame(panelIndex, frame);
        if (endpoint != 0)
            require(frame.hasCellChartPoint
                    ? Geometry::cellChartEqual(frame.cellChartPoint,endpoint < 0 ? graph.localBegin : graph.localEnd)
                    : same(frame.cellLocalPoint,endpoint < 0 ? graph.localBegin : graph.localEnd),
                "charge endpoint differs from its retained edge-root parameter");
        return frame;
    }
    require(!source_->hasCornerSectors(), "sector charge quadrature lost its original local provenance");
    const double nan = std::numeric_limits<double>::quiet_NaN();
    Vector2D local{nan,nan}, actual{nan,nan};
    const char* stage = "legacy local coordinate";
    try
    {
        if (endpoint != 0)
            local = endpoint < 0 ? graph.localBegin : graph.localEnd;
        else
        {
            // Legacy non-sector panels used physical quadrature points. Preserve
            // THOSE points and their saved GL8 integrals, not a second local GL8
            // sequence paired with old lengths. Invert the SAME cell-origin
            // map as jetAtCellLocal. Dividing relative to the domain origin
            // and then subtracting a large cell index can lose enough bits to
            // move the saved physical point by an ULP. Both implementations
            // disable FP contraction; the origin below has the same two-step
            // binary64 construction as Snapshot::coordinate. Still require a
            // strict interior and exact physical round trip, with no clipping.
            // Sector sub-ULP points never enter this compatibility branch.
            const auto lower = source_->minCoord();
            const Vector2D origin{lower.x+graph.cellI*graph.cellSize.x, lower.y+graph.cellJ*graph.cellSize.y};
            local = {static_cast<double>((static_cast<long double>(frame.point.x)-origin.x)/graph.cellSize.x),
                static_cast<double>((static_cast<long double>(frame.point.y)-origin.y)/graph.cellSize.y)};
            require(std::isfinite(local.x) && std::isfinite(local.y) && local.x > 0.0 && local.x < 1.0
                && local.y > 0.0 && local.y < 1.0,
                "legacy quadrature point has no resolved interior in its retained local cell");
        }
        stage = "legacy local jet";
        const auto result = source_->jetAtCellLocal(graph.cellI, graph.cellJ, local);
        actual = result.point;
        stage = "legacy quadrature physical round trip";
        if (endpoint == 0)
            require(same(result.point, frame.point), "legacy charge quadrature point failed its same-cell physical round trip");
        // Endpoints instead retain the ORIGINAL graph edge-root parameter.
        // GridEdge's interpolation between represented physical vertices and
        // the cell's nominal-h map can display that same root one ULP apart.
        // Never invert the old display point to change the endpoint parameter.
        // Its zero/full saved measure is attached to the same ordered root
        // identity; this does not replace any interior GL8 point or integral.
        // The new local endpoint Frame therefore need not have the old global
        // display point, but must still pass the same-source check below.
        // An interior quadrature point now has the local jet at its unchanged
        // physical position. Its final bits need not match the old domain-
        // normalized jet's differently rounded local coordinate. Do not weaken
        // validateLocalFrame: all consumer/sector frames still have to match
        // this immutable scalar's local jet exactly. Saved Panel nodes,
        // quadrature weights and cumulative GL8 measures are not modified.
        stage = "retained local frame validation";
        validateLocalFrame(panelIndex, result);
        return result;
    }
    catch (const std::exception& error)
    {
        const auto lower = source_->minCoord();
        std::ostringstream reason;
        reason << std::setprecision(17) << "[SurfaceChargeField::retainedFrame] stage=" << stage
            << ", endpoint=" << endpoint << ", cell=(" << graph.cellI << ',' << graph.cellJ
            << "), h=(" << graph.cellSize.x << ',' << graph.cellSize.y << "), cellOrigin=("
            << lower.x+graph.cellI*graph.cellSize.x << ',' << lower.y+graph.cellJ*graph.cellSize.y
            << "), saved=(" << frame.point.x << ',' << frame.point.y << "), local=("
            << local.x << ',' << local.y << "), actual=(" << actual.x << ',' << actual.y << "): " << error.what();
        throw std::runtime_error(reason.str());
    }
}

SurfaceChargeSource SurfaceChargeField::source() const
{
    return source(charge_);
}

SurfaceChargeSource SurfaceChargeField::source(const std::vector<double>& integratedCharge) const
{
    return source(SurfaceChargeTransport::fromIntegratedCharge(
        charge_.panelMeasure(),integratedCharge,panels_.closed()));
}

void SurfaceChargeField::assignIntegratedCharge(const std::vector<double>& integratedCharge)
{
    auto candidate = SurfaceChargeTransport::fromIntegratedCharge(
        charge_.panelMeasure(),integratedCharge,panels_.closed());
    charge_ = std::move(candidate);
}

SurfaceChargeSource SurfaceChargeField::source(const SurfaceChargeTransport& charge) const
{
    SurfaceChargeSource result;
    result.axisymmetric = panels_.axisymmetric();
    result.cellIntegratedCharge.assign(panels_.usesNativeLeaves()
        ?panels_.nativeChain().cellToGraphs().size():panels_.cellToPanel().size(), 0.0);
    const auto& panels = panels_.panels();
    const auto density = charge.density();
    std::vector<Polynomial> polynomial(panels.size());
    for (std::size_t i = 0; i < panels.size(); ++i)
    {
        result.cellIntegratedCharge[panels[i].cellIndex] += charge.panelCharge()[i];
        const auto& mapping = mappings_[i];
        polynomial[i]=reconstructPolynomial(i,density[i],{{density[mapping.neighbor[0]],
            density[mapping.neighbor[1]],density[mapping.neighbor[2]]}});
    }
    result.evaluate = [this, polynomial = std::move(polynomial)](const LevelSetHermiteGeometry::Frame& frame)
    {
        require(source_ || !frame.hasCellChartPoint,
            "legacy global charge source cannot consume a boundary-chart Frame");
        const auto location = source_ ? (panels_.usesNativeLeaves() ? locateNativeCached(frame) : locate(frame))
                                      : locate(frame.point);
        return evaluate(location, polynomial);
    };
    return result;
}

SurfaceChargeField::Polynomial SurfaceChargeField::reconstructPolynomial(std::size_t panel,
    double panelDensity,const std::array<double,3>& neighborDensity) const
{
    const auto& mapping=mappings_.at(panel);
    // Keep the original difference arithmetic and accumulation order. This is
    // shared by the legacy current-Q evaluation and the unit-Q stencil columns.
    const double difference[3]={neighborDensity[0]-panelDensity,
        neighborDensity[1]-panelDensity,neighborDensity[2]-panelDensity};
    double coefficients[3]={};
    for(int row=0;row<3;++row) {
        long double value=0.0L;
        for(int col=0;col<3;++col)value+=mapping.inverseMoments[3*row+col]*difference[col];
        coefficients[row]=static_cast<double>(value);
        require(std::isfinite(coefficients[row]),"charge reconstruction coefficient overflow");
    }
    const Polynomial result={panelDensity+coefficients[0],coefficients[1],coefficients[2],
        mapping.scale,mapping.polynomialOrigin,mapping.evenAtAxis};
    require(std::isfinite(result.a),"charge reconstruction constant overflow");
    return result;
}

SurfaceChargeDensityStencil SurfaceChargeField::densityStencil(std::size_t panel) const
{
    const auto& geometry=panels_.panels().at(panel);
    const auto& mapping=mappings_.at(panel);
    const auto& measure=charge_.panelMeasure();
    SurfaceChargeDensityStencil result;
    result.neighbor=mapping.neighbor;result.origin=mapping.polynomialOrigin;
    result.scale=mapping.scale;result.panelMeasure=geometry.measure;result.evenAtAxis=mapping.evenAtAxis;
    require(std::isfinite(result.panelMeasure) && result.panelMeasure>0
        && std::isfinite(result.scale) && result.scale>0
        && std::isfinite(result.origin) && result.origin>=0 && result.origin<=result.panelMeasure,
        "invalid retained surface-charge stencil coordinates");
    const double ownMeasure=measure.at(panel);
    require(std::isfinite(ownMeasure) && ownMeasure>0,"invalid physical charge panel measure");
    for(std::size_t col=0;col<3;++col) {
        const std::size_t unitPanel=result.neighbor[col];
        for(std::size_t previous=0;previous<col;++previous)
            require(result.neighbor[previous]!=unitPanel,"duplicate surface-charge density stencil neighbor");
        const double unitMeasure=measure.at(unitPanel);
        require(std::isfinite(unitMeasure) && unitMeasure>0,
            "invalid physical neighboring charge panel measure");
        const double unitDensity=1.0/unitMeasure;
        require(std::isfinite(unitDensity) && unitDensity>0,"unit-Q density stencil is unrepresentable");
        std::array<double,3> neighborDensity{};
        for(std::size_t neighbor=0;neighbor<3;++neighbor)
            if(result.neighbor[neighbor]==unitPanel)neighborDensity[neighbor]=unitDensity;
        const auto polynomial=reconstructPolynomial(panel,panel==unitPanel?unitDensity:0.,neighborDensity);
        result.polynomialWeights[col]=polynomial.a;
        result.polynomialWeights[3+col]=polynomial.b;
        result.polynomialWeights[6+col]=polynomial.c;
    }
    return result;
}

std::size_t SurfaceChargeField::prepareMapping(std::vector<MeasureNode>& nodes,
    const char* orderReason, const char* monotonicityReason, bool zeroBeginDerivative, bool zeroEndDerivative)
{
    std::size_t rational=0;
    for (std::size_t j = 1; j < nodes.size(); ++j)
    {
        require(nodes[j].coordinate > nodes[j - 1].coordinate &&
            nodes[j].measure > nodes[j - 1].measure, orderReason);
        auto& a = nodes[j - 1]; const auto& b = nodes[j];
        const bool zeroA = zeroBeginDerivative && j == 1 && a.derivative == 0.0;
        const bool zeroB = zeroEndDerivative && j + 1 == nodes.size() && b.derivative == 0.0;
        require(std::isfinite(a.derivative) && std::isfinite(b.derivative) &&
            (a.derivative > 0.0 || zeroA) && (b.derivative > 0.0 || zeroB), monotonicityReason);
        const double slope = (b.measure-a.measure)/(b.coordinate-a.coordinate);
        require(std::isfinite(slope) && slope>0,monotonicityReason);
        const double qa = 3*(a.derivative+b.derivative-2*slope);
        const double qb = 2*(3*slope-2*a.derivative-b.derivative);
        double minimum = std::min(a.derivative, b.derivative);
        bool interiorNonpositive=false;
        if (qa > 0.0)
        {
            const double t = -qb/(2*qa);
            if (t > 0.0 && t < 1.0)
            {
                const double interior = (qa*t+qb)*t+a.derivative;
                require(std::isfinite(interior), monotonicityReason);
                interiorNonpositive=interior<=0;
                minimum = std::min(minimum, interior);
            }
        }
        require(std::isfinite(qa) && std::isfinite(qb) && std::isfinite(minimum),monotonicityReason);
        a.rationalToNext=interiorNonpositive || !(minimum>0 || ((zeroA || zeroB) && minimum==0));
        if(a.rationalToNext)++rational;
    }
    return rational;
}

double SurfaceChargeField::interpolateMapping(const std::vector<MeasureNode>& nodes, double coordinate)
{
    const auto upper = std::upper_bound(nodes.begin(), nodes.end(), coordinate,
        [](double u, const MeasureNode& node) { return u < node.coordinate; });
    const std::size_t right = std::min<std::size_t>(nodes.size() - 1,
        std::max<std::size_t>(1, static_cast<std::size_t>(upper - nodes.begin())));
    const auto& a = nodes[right - 1]; const auto& b = nodes[right];
    const double interval = b.coordinate - a.coordinate, t = (coordinate - a.coordinate) / interval;
    if(a.rationalToNext) {
        // Gregory--Delbourgo monotone rational quadratic: retain BOTH saved
        // cumulative integrals and actual endpoint metrics, without a slope
        // limiter. Positive secant delta and nonnegative endpoint derivatives
        // give f' = delta^2*(db*t^2+2*delta*t*(1-t)+da*(1-t)^2)/den^2 > 0
        // in the interval. A zero endpoint derivative is allowed only at a
        // declared RZ pole, as checked above. Adjacent cubic pieces share
        // these same endpoint values/derivatives, so the map remains C1.
        // Formula: https://doi.org/10.1093/imanum/2.2.123
        require(std::isfinite(t) && t>=0 && t<=1,"rational measure query left its interval");
        if(coordinate==a.coordinate)return a.measure;
        if(coordinate==b.coordinate)return b.measure;
        const double delta=(b.measure-a.measure)/interval;
        const double scale=std::max({delta,a.derivative,b.derivative});
        const double s=delta/scale,da=a.derivative/scale,db=b.derivative/scale,u=1-t;
        require(s>0 && (a.derivative==0 || da>0) && (b.derivative==0 || db>0),
            "rational measure metric ratio underflow");
        // Complementary positive numerators avoid subtracting the large
        // endpoint metrics in the denominator or 1-f near the right knot.
        const double left=s*t*t+da*t*u,right=s*u*u+db*t*u,den=left+right;
        require(std::isfinite(den) && den>0 && left>=0 && right>=0,
            "rational measure interval is unrepresentable");
        return t<=.5?a.measure+(b.measure-a.measure)*(left/den)
            :b.measure-(b.measure-a.measure)*(right/den);
    }
    return (2*t*t*t - 3*t*t + 1)*a.measure + (t*t*t - 2*t*t + t)*interval*a.derivative +
        (-2*t*t*t + 3*t*t)*b.measure + (t*t*t - t*t)*interval*b.derivative;
}

SurfaceChargeField::Location SurfaceChargeField::locate(Vector2D point) const
{
    require(std::isfinite(point.x) && std::isfinite(point.y), "nonfinite charge evaluation point");
    const auto& panels = panels_.panels();
    for (std::size_t i = 0; i < panels.size(); ++i)
    {
        const auto& mapping = mappings_[i];
        if (point.x < mapping.low.x - coordinateTolerance_ ||
            point.x > mapping.high.x + coordinateTolerance_ ||
            point.y < mapping.low.y - coordinateTolerance_ || point.y > mapping.high.y + coordinateTolerance_)
            continue;
        const double rawCoordinate =
            ((panels[i].coordinateIsX ? point.x : point.y) - mapping.coordinateBegin) /
            mapping.coordinateSpan;
        const double tolerance = coordinateTolerance_ / std::abs(mapping.coordinateSpan);
        if (rawCoordinate < -tolerance || rawCoordinate > 1.0 + tolerance)
            continue;
        // Only roundoff at a shared endpoint may be clipped. The mapping is a
        // cubic Hermite interpolation of saved GL8 arc integrals, not a new
        // projection or a claim of exact arclength on the implicit curve.
        const double coordinate = std::clamp(rawCoordinate, 0.0, 1.0);
        const double s = interpolateMapping(mapping.arc, coordinate);
        require(std::isfinite(s) && s >= -coordinateTolerance_ &&
                    s <= panels[i].measure + coordinateTolerance_,
                "charge arclength interpolation left its panel");
        const double arclength = std::clamp(s, 0.0, panels[i].measure);
        double physicalMeasure = arclength;
        if (panels_.axisymmetric())
        {
            physicalMeasure = interpolateMapping(mapping.physical, coordinate);
            const double tolerance = coordinateTolerance_ * CartesianMeasure::radialWeight(mapping.high.x);
            require(std::isfinite(tolerance) && std::isfinite(physicalMeasure) &&
                        physicalMeasure >= -tolerance &&
                        physicalMeasure <= panels[i].physicalMeasure + tolerance,
                    "charge physical measure interpolation left its panel");
            physicalMeasure = std::clamp(physicalMeasure, 0.0, panels[i].physicalMeasure);
        }
        return {i, arclength, physicalMeasure};
    }
    throw std::runtime_error("interface point does not belong to the retained charge panels");
}

SurfaceChargeField::Location SurfaceChargeField::locate(const LevelSetHermiteGeometry::Frame &frame) const
{
    require(frame.hasCellLocalPoint != frame.hasCellChartPoint,
            "retained charge source requires exactly one local/chart Frame; global fallback is unsupported");
    std::size_t i;
    if(panels_.usesNativeLeaves()) {
        require(frame.hasCellChartPoint,"native charge source requires a retained chart Frame");
        i=panels_.nativeChain().graphAt(frame.nativeCellIndex,frame.cellChartPoint);
    } else {
        const auto found=panelByCell_.find({frame.cellI,frame.cellJ});
        require(found!=panelByCell_.end(),"no retained charge panel for Frame cell");
        i=found->second;
    }
    const auto &panel = panels_.panels()[i];
    try
    {
        validateLocalFrame(i, frame);
        double coordinate;
        if (panels_.usesNativeLeaves())
            coordinate = nativeFrameFraction(i, frame);
        else
        {
            const auto &graph = panels_.panelGraph(static_cast<int>(i));
            const double first = graph.coordinateIsX ? graph.localBegin.x : graph.localBegin.y;
            const double last = graph.coordinateIsX ? graph.localEnd.x : graph.localEnd.y;
            if (frame.hasCellChartPoint)
            {
                const auto point = LevelSetHermiteGeometry::cellChartPoint(frame);
                const auto q = graph.coordinateIsX ? point.x : point.y;
                const bool atFirst = Geometry::compareCellChartCoordinate(q, {first, false}) == 0;
                const bool atLast = Geometry::compareCellChartCoordinate(q, {last, false}) == 0;
                if (atFirst || atLast)
                {
                    require(Geometry::cellChartEqual(point, atFirst ? graph.localBegin : graph.localEnd),
                            "charge chart is not the original retained root at this graph coordinate");
                    coordinate = atFirst ? 0.0 : 1.0;
                }
                else
                {
                    require(Geometry::cellChartInterior(point),
                            "interior charge chart is outside its retained cell");
                    // Exact chart subtraction/division is public Core geometry.
                    // Terminal 0/1 rounding still rejects; it is NOT an endpoint.
                    coordinate = Geometry::cellChartIntervalFraction(q, first, last);
                    require(std::isfinite(coordinate) && coordinate > 0.0 && coordinate < 1.0,
                            "interior charge chart has no represented graph parameter");
                }
            }
            else
            {
                const double q = graph.coordinateIsX ? frame.cellLocalPoint.x : frame.cellLocalPoint.y;
                require(q >= std::min(first, last) && q <= std::max(first, last),
                        "charge frame lies outside its retained graph-coordinate interval");
                if (q == first || q == last)
                {
                    require(same(frame.cellLocalPoint, q == first ? graph.localBegin : graph.localEnd),
                            "charge frame is not the retained endpoint at this graph coordinate");
                    coordinate = q == first ? 0.0 : 1.0;
                }
                else
                {
                    require(frame.cellLocalPoint.x > 0.0 && frame.cellLocalPoint.x < 1.0 &&
                                frame.cellLocalPoint.y > 0.0 && frame.cellLocalPoint.y < 1.0,
                            "interior charge frame is not inside its retained local cell");
                    coordinate = static_cast<double>((static_cast<long double>(q) - first) /
                                                     (static_cast<long double>(last) - first));
                    require(std::isfinite(coordinate) && coordinate > 0.0 && coordinate < 1.0,
                            "interior charge frame has no represented graph parameter");
                }
            }
        }
        return locateFraction(i,coordinate);
    }
    catch (const std::exception &error)
    {
        std::ostringstream reason;
        reason << std::setprecision(17) << "[SurfaceChargeField] panel=" << i
               << ", nativeCell=" << panel.cellIndex << ", Frame cell=(" << frame.cellI << ',' << frame.cellJ
               << ')';
        if (frame.hasCellChartPoint)
            reason << ", boundaryChart=((offset=" << frame.cellChartPoint.x.offset
                   << ",upper=" << frame.cellChartPoint.x.fromUpper
                   << "),(offset=" << frame.cellChartPoint.y.offset
                   << ",upper=" << frame.cellChartPoint.y.fromUpper << "))";
        else
            reason << ", local=(" << frame.cellLocalPoint.x << ',' << frame.cellLocalPoint.y << ')';
        reason << ": " << error.what();
        throw std::runtime_error(reason.str());
    }
}

SurfaceChargeField::Location SurfaceChargeField::locateNativeCached(
    const LevelSetHermiteGeometry::Frame& frame) const
{
    require(frame.hasCellChartPoint && !frame.hasCellLocalPoint && frame.cellI == -1 && frame.cellJ == -1,
        "native charge frame lost its actual leaf/chart identity");
    // This public Core check also rejects unsupported rounding/underflow
    // modes on cache hits. Do not bypass the exact-geometry environment gate.
    Geometry::validateCellChart(frame.cellChartPoint);
    const auto bits = [](double value) {
        std::uint64_t result;
        static_assert(sizeof(result) == sizeof(value));
        std::memcpy(&result, &value, sizeof(result));
        return result;
    };
    // All inputs consumed by native validateLocalFrame/graphFraction belong
    // to the key. No padding bytes or NaN floating-point ordering are used.
    // Stationarity/iteration diagnostics and the unused local representation
    // are deliberately excluded, exactly as in the original validation.
    const NativeLocationKey key{{static_cast<std::uint64_t>(frame.nativeCellIndex),
        bits(frame.cellChartPoint.x.offset), frame.cellChartPoint.x.fromUpper,
        bits(frame.cellChartPoint.y.offset), frame.cellChartPoint.y.fromUpper,
        bits(frame.point.x), bits(frame.point.y), bits(frame.normal.x), bits(frame.normal.y),
        bits(frame.curvature), bits(frame.levelSetResidual)}};
    {
        const std::lock_guard<std::mutex> lock(nativeLocationMutex_);
        const auto found = nativeLocations_.find(key);
        if (found != nativeLocations_.end()) return found->second;
    }
    // A miss executes the unchanged full source-jet, endpoint identity and
    // saved GL8 arclength/area checks. Failed queries are never inserted.
    const auto location = locate(frame);
    {
        const std::lock_guard<std::mutex> lock(nativeLocationMutex_);
        // Bound memory for arbitrary caller query sequences. Reaching this
        // capacity only stops insertion; it never drops or approximates a
        // query. Field destruction/remapping releases the entire cache.
        constexpr std::size_t maximumCachedLocations = 4096;
        if (nativeLocations_.size() < maximumCachedLocations)
            nativeLocations_.emplace(key, location);
    }
    return location;
}

SurfaceChargeField::Location SurfaceChargeField::locateFraction(std::size_t i, double coordinate) const
{
    require(std::isfinite(coordinate) && coordinate >= 0 && coordinate <= 1,
            "material measure query has an invalid graph fraction");
    const auto &panel = panels_.panels().at(i);
    const auto &mapping = mappings_.at(i);
    // Saved GL8 integrals with the admitted cubic or monotone rational map.
    // No clamping or normalization of an unresolved subarc or surface area.
    const double arclength = interpolateMapping(mapping.arc, coordinate);
    require(std::isfinite(arclength) && arclength >= 0 && arclength <= panel.measure,
            "cell-local charge arclength interpolation left its panel");
    const double physicalMeasure =
        panels_.axisymmetric() ? interpolateMapping(mapping.physical, coordinate) : arclength;
    require(std::isfinite(physicalMeasure) && physicalMeasure >= 0 &&
                physicalMeasure <= panel.physicalMeasure,
            "cell-local charge physical measure interpolation left its panel");
    return {i, arclength, physicalMeasure};
}

SurfaceChargeJet SurfaceChargeField::evaluate(const Location &location,
                                              const std::vector<Polynomial> &polynomial) const
{
    const auto &p = polynomial[location.panel];
    const double x = (location.arclength - p.origin) / p.scale;
    SurfaceChargeJet result;
    if (p.evenAtAxis)
    {
        const double square = x * x;
        result = {p.a + square * (p.b + square * p.c), x * (2 * p.b + 4 * square * p.c) / p.scale,
                  ((2 * p.b + 12 * square * p.c) / p.scale) / p.scale};
    }
    else
        result = {p.a + x * (p.b + x * p.c), (p.b + 2 * x * p.c) / p.scale, (2 * p.c / p.scale) / p.scale};
    require(std::isfinite(result.density) && std::isfinite(result.firstDerivative) &&
                std::isfinite(result.secondDerivative),
            "charge jet overflow");
    return result;
}

SurfaceChargeFieldStep SurfaceChargeField::advanceNormalCurrent(const std::vector<double> &cellCurrent,
                                                                double currentScale,
                                                                double relativeBulkTolerance, double dt)
{
    return advance(cellCurrent,
                   std::vector<double>(panels_.panels().size() + (panels_.closed() ? 0 : 1), 0.0),
                   currentScale, relativeBulkTolerance, dt);
}

SurfaceChargeFieldStep SurfaceChargeField::advance(const std::vector<double> &cellCurrent,
                                                   const std::vector<double> &relativeAdvectiveMeasureFlux,
                                                   double currentScale, double relativeBulkTolerance,
                                                   double dt)
{
    SurfaceChargeFieldStep result;
    const auto current = panelCurrent(cellCurrent, currentScale, relativeBulkTolerance, result);
    result.transport = charge_.advance(charge_.panelMeasure(), relativeAdvectiveMeasureFlux, current, dt);
    return result;
}

SurfaceChargeFieldStep SurfaceChargeField::advance(const std::vector<double>& cellCurrent,
    const std::vector<double>& relativeAdvectiveMeasureFlux,
    double currentScale, double relativeBulkTolerance, double dt,
    const SurfaceChargeBoundaryState& boundary)
{
    SurfaceChargeFieldStep result;
    const auto current = panelCurrent(cellCurrent, currentScale, relativeBulkTolerance, result);
    if (!panels_.closed())
    {
        require(relativeAdvectiveMeasureFlux.size() == panels_.panels().size() + 1,
            "open surface-charge field requires one measure flux per shared boundary");
        require((panels_.beginEndpointKind() != LevelSetInterfacePanels::EndpointKind::Axis ||
                relativeAdvectiveMeasureFlux.front() == 0.0) &&
            (panels_.endEndpointKind() != LevelSetInterfacePanels::EndpointKind::Axis ||
                relativeAdvectiveMeasureFlux.back() == 0.0),
            "surface-charge symmetry-axis endpoint requires exactly zero measure flux");
    }
    result.transport = charge_.advance(charge_.panelMeasure(), relativeAdvectiveMeasureFlux, current, dt, boundary);
    return result;
}

std::vector<double> SurfaceChargeField::panelCurrent(const std::vector<double>& cellCurrent,
    double currentScale, double relativeBulkTolerance, SurfaceChargeFieldStep& result) const
{
    const bool native=panels_.usesNativeLeaves();
    require(cellCurrent.size() == (native?panels_.nativeChain().cellToGraphs().size():panels_.cellToPanel().size()),
        "charge/current cell count mismatch");
    require(std::isfinite(currentScale) && currentScale >= 0.0 && std::isfinite(relativeBulkTolerance) &&
        relativeBulkTolerance > 0.0 && relativeBulkTolerance < 1.0, "invalid bulk current tolerance/scale");
    std::vector<double> panelCurrent(panels_.panels().size(), 0.0);
    long double bulk = 0.0L, bulkAbs = 0.0L;
    for (std::size_t i = 0; i < cellCurrent.size(); ++i)
    {
        require(std::isfinite(cellCurrent[i]), "nonfinite conduction cell current");
        if(native) {
            const auto& pieces=panels_.nativeChain().cellToGraphs()[i];
            if(pieces.empty()) {bulk+=cellCurrent[i];bulkAbs+=std::abs(cellCurrent[i]);continue;}
            if(pieces.size()==1) {panelCurrent[pieces.front()]=cellCurrent[i];continue;}
            // The conduction solver supplies one FV current integral per
            // cell. Its explicit subpanel closure is uniform normal-current
            // density over that cell's retained physical interface measure.
            long double area=0,assigned=0;int largest=pieces.front();
            for(int p:pieces) {
                area+=panels_.panels()[p].physicalMeasure;
                if(panels_.panels()[p].physicalMeasure>panels_.panels()[largest].physicalMeasure)largest=p;
            }
            require(std::isfinite(area) && area>0,"native current has no positive interface measure");
            for(int p:pieces) {
                panelCurrent[p]=static_cast<double>(static_cast<long double>(cellCurrent[i])
                    *panels_.panels()[p].physicalMeasure/area);
                assigned+=panelCurrent[p];
            }
            // Close summation roundoff on the largest piece, not a vanishing
            // pole panel. No non-interface/bulk current is redistributed.
            panelCurrent[largest]+=static_cast<double>(static_cast<long double>(cellCurrent[i])-assigned);
        } else {
            const int panel=panels_.cellToPanel()[i];
            if(panel>=0)panelCurrent[panel]=cellCurrent[i];
            else {bulk+=cellCurrent[i];bulkAbs+=std::abs(cellCurrent[i]);}
        }
    }
    result.unassignedBulkCurrent = static_cast<double>(bulk);
    result.absoluteBulkCurrent = static_cast<double>(bulkAbs);
    require(std::isfinite(result.unassignedBulkCurrent) && std::isfinite(result.absoluteBulkCurrent) &&
        bulkAbs <= static_cast<long double>(relativeBulkTolerance)*currentScale,
        "non-interface conduction current exceeds roundoff allowance; bulk charge is not supported");
    // Bulk roundoff is reported, never moved into an arbitrary interface cell.
    return panelCurrent;
}

std::vector<double> SurfaceChargeField::mapNormalCurrent(const std::vector<double>& cellCurrent,
    double currentScale, double relativeBulkTolerance, SurfaceChargeFieldStep& report) const
{
    SurfaceChargeFieldStep staged;
    auto current=panelCurrent(cellCurrent,currentScale,relativeBulkTolerance,staged);
    report.unassignedBulkCurrent=staged.unassignedBulkCurrent;
    report.absoluteBulkCurrent=staged.absoluteBulkCurrent;
    return current;
}

std::unique_ptr<SurfaceChargeField> SurfaceChargeField::advanceMoving(const QuadtreeMesh& mesh,
    const LevelSetHermiteGeometry::Snapshot& nextSnapshot,
    const std::function<Vector2D(Vector2D)>& forwardMaterialPoint,
    const std::vector<double>& cellCurrent, double currentScale,
    double relativeBulkTolerance, double dt, SurfaceChargeFieldStep& report) const
{
    require(!source_, "moving cell-local surface charge requires a retained-chart remap; legacy remapping cannot discard its geometry mode");
    require(panels_.boundaryContacts().empty(),
        "moving declared-boundary surface charge requires a boundary-aware remap; fixed contacts are not material endpoints");
    require(static_cast<bool>(forwardMaterialPoint), "moving charge requires a material trajectory");
    const bool closed = panels_.closed();
    auto next = std::make_unique<SurfaceChargeField>(mesh, nextSnapshot, 0.0, panels_.axisymmetric());
    require(next->panels_.closed() == closed, "moving surface charge cannot change interface topology");
    require(closed || (panels_.axisymmetric() && nextSnapshot.radialAxisSymmetry()),
        "moving open charge requires radially symmetric axis geometry");
    SurfaceChargeFieldStep result;
    const auto current = panelCurrent(cellCurrent, currentScale, relativeBulkTolerance, result);
    auto staged = SurfaceChargeTransport::fromIntegratedCharge(charge_.panelMeasure(), charge_.panelCharge(), closed);
    result.transport = staged.advance(staged.panelMeasure(), std::vector<double>(current.size() + (closed ? 0 : 1), 0.0), current, dt);
    if (!result.transport.success) throw std::runtime_error("moving charge " +
        result.transport.stage + ": " + result.transport.reason);
    const auto nextMeasures = measures(next->panels_);
    std::vector<double> boundaries(nextMeasures.size()+1, 0.0);
    for (std::size_t j = 0; j < nextMeasures.size(); ++j)
    {
        boundaries[j+1] = boundaries[j] + nextMeasures[j];
        require(std::isfinite(boundaries[j+1]) && boundaries[j+1] > boundaries[j],
            panels_.axisymmetric() ? "unresolved new-interface physical measure" : "unresolved new-interface arclength");
    }
    const double period = boundaries.back();
    const double h = std::min(mesh.cellSize.x, mesh.cellSize.y);
    const auto& oldPanels = panels_.panels();
    const auto& newPanels = next->panels_.panels();
    std::vector<double> mapped(oldPanels.size()+1);
    std::vector<double> endpointCoordinates(oldPanels.size() + (closed ? 0 : 1));
    for (std::size_t i = 0; i < endpointCoordinates.size(); ++i)
    {
        const bool beginPole = !closed && i == 0;
        const bool endPole = !closed && i == oldPanels.size();
        const auto oldPoint = endPole ? oldPanels.back().end : oldPanels[i].begin;
        const auto& oldFrame = endPole ? oldPanels.back().endFrame : oldPanels[i].beginFrame;
        const auto moved = forwardMaterialPoint(oldPoint);
        const double displacement = std::hypot(moved.x-oldPoint.x, moved.y-oldPoint.y);
        require(std::isfinite(displacement) && displacement <= .5*h,
            "material interface displacement exceeds half a cell; reduce dt");
        LevelSetHermiteGeometry::Frame frame;
        if (beginPole || endPole)
        {
            require(oldPoint.x == 0.0 && moved.x == 0.0 &&
                moved.y > mesh.minCoord.y && moved.y < mesh.maxCoord.y,
                "material pole trajectory must stay exactly on the interior r=0 axis");
            // The chain's oriented endpoint is the actual shared-snapshot
            // axis crossing, not a fitted pole or a clipped radial trajectory.
            frame = beginPole ? newPanels.front().beginFrame : newPanels.back().endFrame;
            require(frame.point.x == 0.0 && frame.normal.x == 0.0,
                "new material pole is not a radially symmetric axis crossing");
            frame.stationarityResidual = std::abs(-(frame.point.x-moved.x)*frame.normal.y +
                (frame.point.y-moved.y)*frame.normal.x);
        }
        else
        {
            require(closed || moved.x > 0.0, "interior material endpoint reached or crossed the symmetry axis");
            frame = nextSnapshot.projectedFrame(moved);
        }
        const double projection = std::hypot(frame.point.x-moved.x, frame.point.y-moved.y);
        require(std::isfinite(projection) && projection <= .25*h &&
            std::isfinite(frame.levelSetResidual) && frame.levelSetResidual <= 1e-9*h &&
            std::isfinite(frame.stationarityResidual) && frame.stationarityResidual <= 1e-9*h,
            "material endpoint cannot be locally projected to the new interface");
        const auto location = next->locate(frame.point);
        require((!beginPole || (location.panel == 0 && location.localPhysicalMeasure == 0.0)) &&
            (!endPole || (location.panel+1 == newPanels.size() &&
                location.localPhysicalMeasure == newPanels.back().physicalMeasure)),
            "material pole does not match the same oriented new-chain endpoint");
        // Exclude a nearby nonadjacent panel branch from the projection's
        // correction ball. This is a local resolved-motion guard, not a
        // globally closest-point certificate for arbitrary topology.
        for (std::size_t j = 0; j < next->mappings_.size(); ++j)
        {
            if (j == location.panel || (closed
                ? (j == (location.panel+1)%nextMeasures.size() || (j+1)%nextMeasures.size() == location.panel)
                : (j == location.panel+1 || j+1 == location.panel))) continue;
            const auto& box = next->mappings_[j];
            const double dx = std::max({box.low.x-moved.x, 0.0, moved.x-box.high.x});
            const double dy = std::max({box.low.y-moved.y, 0.0, moved.y-box.high.y});
            require(std::hypot(dx,dy) > projection + next->coordinateTolerance_,
                "material projection is ambiguous near a nonadjacent interface branch");
        }
        require(oldFrame.normal.x*frame.normal.x + oldFrame.normal.y*frame.normal.y > 0.0,
            "material interface normal reversed during a step");
        endpointCoordinates[i] = boundaries[location.panel] + location.localPhysicalMeasure;
        if (closed)
        {
            if (endpointCoordinates[i] == period) endpointCoordinates[i] = 0.0;
            require(endpointCoordinates[i] >= 0.0 && endpointCoordinates[i] < period,
                panels_.axisymmetric() ? "mapped endpoint is outside periodic physical measure" :
                    "mapped endpoint is outside periodic arclength");
        }
        else
            require((beginPole && endpointCoordinates[i] == 0.0) || (endPole && endpointCoordinates[i] == period) ||
                (!beginPole && !endPole && endpointCoordinates[i] > 0.0 && endpointCoordinates[i] < period),
                "mapped endpoint is outside its open physical-area interval");
        result.maxMaterialDisplacement = std::max(result.maxMaterialDisplacement, displacement);
        result.maxProjectionDistance = std::max(result.maxProjectionDistance, projection);
    }
    if (closed)
    {
        mapped[0] = endpointCoordinates[0];
        long double windingLength = 0.0L;
        for (std::size_t i = 0; i < oldPanels.size(); ++i)
        {
            double delta = endpointCoordinates[(i+1)%oldPanels.size()] - endpointCoordinates[i];
            if (delta < 0.0) delta += period;
            require(std::isfinite(delta) && delta > 0.0, "moving interface boundaries merged or folded");
            windingLength += delta;
            mapped[i+1] = mapped[i] + delta;
            require(mapped[i+1] > mapped[i], panels_.axisymmetric() ?
                "material panel physical measure is unresolved" : "material panel arclength is unresolved");
        }
        require(std::abs(windingLength-period) <= 512.0*std::numeric_limits<double>::epsilon()*period,
            "material interface mapping does not have one positive winding");
        mapped.back() = mapped.front()+period; // Same shared periodic endpoint, not a charge correction.
    }
    else
    {
        mapped = endpointCoordinates; // All N+1 material endpoints are distinct; never wrap.
        require(mapped.front() == 0.0 && mapped.back() == period, "open material map does not cover both axis endpoints");
        for (std::size_t i = 1; i < mapped.size(); ++i)
            require(mapped[i] > mapped[i-1], "open material interface boundaries merged or folded");
    }
    const auto remap = staged.remap(mapped, boundaries);
    if (!remap.success) throw std::runtime_error("surface charge remap " + remap.stage + ": " + remap.reason);
    // The overlap coordinates use cumulative physical measures; restore the original
    // GL8 measures without round-tripping or changing ANY integrated charge.
    next->charge_ = SurfaceChargeTransport::fromIntegratedCharge(nextMeasures, staged.panelCharge(), closed);
    result.remapped = true;
    result.previousPanelCount = oldPanels.size();
    result.remapConservationResidual = remap.conservationResidual;
    auto& q = result.transport;
    q.totalChargeAfter = remap.totalChargeAfter;
    q.conservationResidual = q.totalChargeAfter-q.totalChargeBefore-dt*q.integratedNormalCurrent;
    q.relativeConservationResidual = q.conservationScale > 0.0 ? std::abs(q.conservationResidual)/q.conservationScale : 0.0;
    q.maxAbsDensity = 0.0;
    for (double value : next->charge_.density()) q.maxAbsDensity = std::max(q.maxAbsDensity, std::abs(value));
    require(std::isfinite(q.conservationResidual) && std::isfinite(q.relativeConservationResidual) &&
        q.relativeConservationResidual <= 4096.0*std::numeric_limits<double>::epsilon(),
        "moving surface-charge update does not conserve integrated charge");
    report = std::move(result);
    return next;
}

SurfaceChargeField::MaterialRemapPlan SurfaceChargeField::prepareMaterialRemap(
    const QuadtreeMesh &mesh, std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> nextSnapshot,
    const std::vector<Vector2D> &fixedMaterialContacts,
    const std::function<Vector2D(Vector2D)> &forwardMaterialPoint,
    const LevelSetInterfacePanels* preparedPanels) const
{
    require(source_ && nextSnapshot, "cell-local material remap requires old/new retained snapshots");
    require(!panels_.closed(), "cell-local material remap currently requires an open chain");
    require(static_cast<bool>(forwardMaterialPoint), "moving charge requires a material trajectory");
    const bool native = panels_.usesNativeLeaves();
    if (preparedPanels)
    {
        require(native && preparedPanels->usesNativeLeaves() && preparedPanels->usesSource(*nextSnapshot)
            && preparedPanels->axisymmetric()==panels_.axisymmetric(),
            "prepared material panels must use the same native scalar and symmetry");
        const auto& contacts=preparedPanels->boundaryContacts();
        require(contacts.size()==panels_.boundaryContacts().size(),
            "prepared material panels changed declared contacts");
        for (std::size_t i=0;i<contacts.size();++i)
            require(same(contacts[i],panels_.boundaryContacts()[i]),
                "prepared material panels changed declared contact identity/order");
    }
    if (native)
        validateNativeRemapGrid(mesh, *nextSnapshot);
    else
    {
        require(!nextSnapshot->usesNativeLeaves(), "Cartesian material remap cannot change source layout");
        require(source_->nx() == nextSnapshot->nx() && source_->ny() == nextSnapshot->ny() &&
                    same(source_->minCoord(), nextSnapshot->minCoord()) &&
                    same(source_->maxCoord(), nextSnapshot->maxCoord()) &&
                    same(source_->cellSize(), nextSnapshot->cellSize()) &&
                    source_->radialAxisSymmetry() == nextSnapshot->radialAxisSymmetry() &&
                    source_->restrictedDomain() == nextSnapshot->restrictedDomain(),
                "cell-local material remap requires the same old/new scalar grid and domain kind");
        const UniformLeafGrid grid(mesh, scope_);
        require(grid.cells() == retainedCellLabels_ && grid.vertices() == retainedVertexLabels_,
                "cell-local material remap cannot change retained native cell/vertex labels");
        for (int j = 0; j < source_->ny(); ++j)
            for (int i = 0; i < source_->nx(); ++i)
                require(source_->cellActive(i, j) == nextSnapshot->cellActive(i, j),
                        "cell-local material remap cannot change the active/solid cell mask");
        const auto &oldPanels = panels_.panels();
        for (std::size_t i = 0; i < oldPanels.size(); ++i)
        {
            const auto &graph = panels_.panelGraph(static_cast<int>(i));
            const int label = grid.cells()[static_cast<std::size_t>(graph.cellJ) * grid.nx() + graph.cellI];
            require(mesh.cells(label).index() == oldPanels[i].cellIndex,
                    "cell-local material remap changed a retained panel's native solver index");
        }
    }
    const auto &oldPanels = panels_.panels();
    const auto &declared = panels_.boundaryContacts();
    require(fixedMaterialContacts.size() == declared.size(),
            "every declared contact requires an explicit fixed-material identity");
    for (std::size_t i = 0; i < fixedMaterialContacts.size(); ++i)
    {
        const auto point = fixedMaterialContacts[i];
        require(std::isfinite(point.x) && std::isfinite(point.y) &&
                    std::find_if(declared.begin(), declared.end(),
                                 [&](Vector2D p) { return same(p, point); }) != declared.end(),
                "fixed material contact differs from the declared boundary geometry");
        for (std::size_t j = 0; j < i; ++j)
            require(!same(point, fixedMaterialContacts[j]), "duplicate fixed material contact identity");
    }

    SurfaceChargeGeometryOptions options;
    options.axisymmetric = panels_.axisymmetric();
    options.boundaryContacts = declared;
    options.scope = scope_;
    auto next = preparedPanels ? std::make_unique<SurfaceChargeField>(mesh,*preparedPanels,0.0)
        : native ? std::make_unique<SurfaceChargeField>(
                     mesh,
                     LevelSetInterfacePanels::discoverNativeGraphs(mesh, std::move(nextSnapshot),
                                                                   panels_.axisymmetric(), declared),
                     0.0)
               : std::make_unique<SurfaceChargeField>(mesh, std::move(nextSnapshot), 0.0, std::move(options));
    require(!next->panels_.closed() && next->panels_.beginEndpointKind() == panels_.beginEndpointKind() &&
                next->panels_.endEndpointKind() == panels_.endEndpointKind(),
            "moving surface charge cannot change open-chain endpoint topology");
    SurfaceChargeFieldStep result;
    const auto &newPanels = next->panels_.panels();
    const auto spacing = native ? panels_.nativeChain().minimumActiveSpacing() : source_->cellSize();
    const double h = std::min(spacing.x, spacing.y);
    std::vector<SurfaceChargeTransport::PanelCoordinate> mapped(oldPanels.size() + 1);
    struct EndpointMotion
    {
        LevelSetHermiteGeometry::Frame oldFrame;
        Vector2D displacement;
        double motion;
        bool begin, end, fixed, axis;
    };
    std::vector<EndpointMotion> endpoints;
    endpoints.reserve(mapped.size());
    std::size_t fixedContactsVisited = 0;
    for (std::size_t i = 0; i < mapped.size(); ++i)
    {
        const bool begin = i == 0, end = i == oldPanels.size();
        const std::size_t oldPanel = end ? oldPanels.size() - 1 : i;
        const auto oldFrame =
            retainedFrame(oldPanel, end ? oldPanels.back().endFrame : oldPanels[i].beginFrame, end ? 1 : -1);
        require(native ? (oldFrame.hasCellChartPoint && !oldFrame.hasCellLocalPoint)
                       : (oldFrame.hasCellLocalPoint && !oldFrame.hasCellChartPoint),
                "material endpoint projection does not accept boundary-chart anchors; no implicit "
                "local-coordinate conversion");
        const auto kind = begin ? panels_.beginEndpointKind()
                          : end ? panels_.endEndpointKind()
                                : LevelSetInterfacePanels::EndpointKind::None;
        const bool fixed = kind == LevelSetInterfacePanels::EndpointKind::DeclaredBoundary;
        const bool axis = kind == LevelSetInterfacePanels::EndpointKind::Axis;
        // A fixed declaration is the caller's physical trajectory identity.
        // Every other trajectory is sampled at the retained root's display
        // point; only the resulting displacement is added to its LOCAL anchor.
        const auto trajectoryPoint =
            fixed ? (begin ? oldPanels.front().begin : oldPanels.back().end) : oldFrame.point;
        const auto moved =
            forwardMaterialPoint(trajectoryPoint); // Exactly one trajectory evaluation per old boundary.
        require(std::isfinite(moved.x) && std::isfinite(moved.y), "nonfinite material endpoint trajectory");
        const Vector2D displacement = moved - trajectoryPoint;
        const double motion = std::hypot(displacement.x, displacement.y);
        require(std::isfinite(motion) && motion <= .5 * h,
                "material interface displacement exceeds half a cell; reduce dt");
        if (fixed)
        {
            require(std::find_if(fixedMaterialContacts.begin(), fixedMaterialContacts.end(),
                                 [&](Vector2D p) { return same(p, trajectoryPoint); }) !=
                            fixedMaterialContacts.end() &&
                        same(moved, trajectoryPoint),
                    "fixed material contact trajectory must remain exactly at its declaration");
            require(same(begin ? newPanels.front().begin : newPanels.back().end, trajectoryPoint),
                    "new interface moved a fixed declared material endpoint");
            ++fixedContactsVisited;
        }
        else if (axis)
            require(panels_.axisymmetric() && oldFrame.point.x == 0.0 && moved.x == 0.0 &&
                        moved.y > source_->minCoord().y && moved.y < source_->maxCoord().y,
                    "material pole trajectory must stay exactly on the interior r=0 axis");
        else if (panels_.axisymmetric())
            require(moved.x > 0.0, "interior material endpoint reached or crossed the symmetry axis");

        endpoints.push_back({oldFrame, displacement, motion, begin, end, fixed, axis});
    }
    // Keep arbitrary material-trajectory callbacks on the caller, exactly
    // once per endpoint and in their original order. Only the subsequent
    // independent geometry queries consume the immutable new interface in
    // parallel. Each task owns its mapping and diagnostic slots.
    std::vector<double> projectionDistances(endpoints.size());
    const auto projectEndpoint = [&](std::size_t i) {
        const auto& endpoint = endpoints[i];
        const auto& oldFrame = endpoint.oldFrame;
        const auto displacement = endpoint.displacement;
        const bool begin = endpoint.begin, end = endpoint.end;
        const bool fixed = endpoint.fixed, axis = endpoint.axis;
        LevelSetInterfacePanels::LocalProjection projection;
        const int endpointPanel = begin ? 0 : static_cast<int>(newPanels.size() - 1);
        if (native)
        {
            const LevelSetInterfacePanels::NativeProjectionQuery query{oldFrame.nativeCellIndex,
                                                                       oldFrame.cellChartPoint, displacement};
            const auto &chain = next->panels_.nativeChain();
            projection = (fixed || axis) ? chain.projectLocalEndpoint(query, endpointPanel, begin)
                                         : chain.projectLocal(query);
        }
        else
        {
            const LevelSetInterfacePanels::LocalProjectionQuery query{oldFrame.cellI, oldFrame.cellJ,
                                                                      oldFrame.cellLocalPoint, displacement};
            projection = (fixed || axis) ? next->panels_.projectLocalEndpoint(query, endpointPanel, begin)
                                         : next->panels_.projectLocal(query);
        }
        const auto &frame = projection.frame;
        require(std::isfinite(projection.distance) && projection.distance >= 0.0 &&
                    projection.distance <= .25 * h && std::isfinite(frame.levelSetResidual) &&
                    frame.levelSetResidual >= 0.0 && frame.levelSetResidual <= 1e-9 * h &&
                    std::isfinite(frame.stationarityResidual) && frame.stationarityResidual >= 0.0 &&
                    frame.stationarityResidual <= 1e-9 * h && std::isfinite(projection.kktResidual) &&
                    projection.kktResidual >= 0.0 && projection.kktResidual <= 1e-9 * h,
                "material endpoint cannot be locally projected to the new retained interface");
        Location location;
        if (native)
        {
            require(projection.panelIndex >= 0, "native material projection lost its panel index");
            next->validateLocalFrame(projection.panelIndex, frame);
            const auto retained =
                next->panels_.nativeChain().frameAt(projection.panelIndex, projection.nativeGraphFraction);
            require(Geometry::cellChartEqual(retained.cellChartPoint, frame.cellChartPoint),
                    "native projection fraction differs from its retained chart");
            // Consume the parameter that actually generated the projection,
            // not a quotient recovered from its rounded terminal chart.
            location = next->locateFraction(projection.panelIndex, projection.nativeGraphFraction);
        }
        else
            location = next->locate(frame);
        require(projection.panelIndex >= 0 &&
                    static_cast<std::size_t>(projection.panelIndex) == location.panel,
                "material projection lost its retained panel identity");
        const auto &locatedPanel = newPanels[location.panel];
        if (projection.endpoint)
        {
            LevelSetInterfacePanels::PanelGraph::Endpoint root;
            if (native)
            {
                const auto &graph = next->panels_.nativeChain().graphs().at(projection.panelIndex);
                const auto endpoint = projection.atPanelBegin ? graph.begin : graph.end;
                root = {endpoint.side, endpoint.rootOrdinal};
            }
            else
            {
                const auto &graph = next->panels_.panelGraph(projection.panelIndex);
                root = projection.atPanelBegin ? graph.beginRoot : graph.endRoot;
            }
            require(projection.root.side == root.side && projection.root.rootOrdinal == root.rootOrdinal,
                    "material projection lost its original graph edge-root identity");
            require(projection.atPanelBegin
                        ? location.arclength == 0.0 && location.localPhysicalMeasure == 0.0
                        : location.arclength == locatedPanel.measure &&
                              location.localPhysicalMeasure == locatedPanel.physicalMeasure,
                    "material projection root identity disagrees with its endpoint measures");
        }
        else
            // A rounded cumulative measure is not an endpoint identity. In
            // particular, a positive RZ tail can be below the panel total's
            // resolution near the axis; never normalize it to a shared root.
            require(location.arclength > 0.0 && location.arclength < locatedPanel.measure &&
                        location.localPhysicalMeasure > 0.0 &&
                        location.localPhysicalMeasure < locatedPanel.physicalMeasure,
                    "interior material projection has an unrepresentable positive subarc length or physical "
                    "subarea");
        if (begin || end)
        {
            const bool correctEndpoint =
                begin ? location.panel == 0 && location.localPhysicalMeasure == 0.0
                      : location.panel + 1 == newPanels.size() &&
                            location.localPhysicalMeasure == newPanels.back().physicalMeasure;
            require(projection.endpoint && projection.atPanelBegin == begin && correctEndpoint,
                    "material endpoint does not match the same oriented new-chain endpoint");
            if (axis)
                require(frame.point.x == 0.0 && frame.normal.x == 0.0,
                        "new material pole is not a radially symmetric axis crossing");
        }
        else
            require(!(location.panel == 0 && location.localPhysicalMeasure == 0.0) &&
                        !(location.panel + 1 == newPanels.size() &&
                          location.localPhysicalMeasure == newPanels.back().physicalMeasure),
                    "interior material boundary reached an open-chain endpoint");

        // Both public local projection entry points own the common whole-cell-
        // box/nonadjacent-branch guard and its retained-coordinate distances.
        // Do not add a second global-display bbox projection policy here.
        require(oldFrame.normal.x * frame.normal.x + oldFrame.normal.y * frame.normal.y > 0.0,
                "material interface normal reversed during a step");
        mapped[i] = {location.panel, location.localPhysicalMeasure};
        projectionDistances[i] = projection.distance;
    };
    if (!endpoints.empty()) projectEndpoint(0);
    std::vector<std::exception_ptr> failures(endpoints.size());
#pragma omp parallel for schedule(dynamic, 1) if(native && endpoints.size() > 8)
    for (std::ptrdiff_t i = 1; i < static_cast<std::ptrdiff_t>(endpoints.size()); ++i)
    {
        try { projectEndpoint(static_cast<std::size_t>(i)); }
        catch (...) { failures[i] = std::current_exception(); }
    }
    for (std::size_t i = 0; i < endpoints.size(); ++i)
    {
        if (failures[i]) std::rethrow_exception(failures[i]);
        result.maxMaterialDisplacement = std::max(result.maxMaterialDisplacement, endpoints[i].motion);
        result.maxProjectionDistance = std::max(result.maxProjectionDistance, projectionDistances[i]);
    }
    require(fixedContactsVisited == fixedMaterialContacts.size(),
            "a declared fixed material contact was not traced");
    return {std::move(next), std::move(mapped), std::move(result)};
}

std::unique_ptr<SurfaceChargeField> SurfaceChargeField::advanceMoving(
    const QuadtreeMesh &mesh, std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> nextSnapshot,
    const std::vector<Vector2D> &fixedMaterialContacts,
    const std::function<Vector2D(Vector2D)> &forwardMaterialPoint, const std::vector<double> &cellCurrent,
    double currentScale, double relativeBulkTolerance, double dt, SurfaceChargeFieldStep &report,
    const LevelSetInterfacePanels* preparedPanels) const
{
    auto plan =
        prepareMaterialRemap(mesh, std::move(nextSnapshot), fixedMaterialContacts, forwardMaterialPoint,preparedPanels);
    auto &next = plan.field;
    auto &result = plan.report;
    const auto current = panelCurrent(cellCurrent, currentScale, relativeBulkTolerance, result);
    auto staged =
        SurfaceChargeTransport::fromIntegratedCharge(charge_.panelMeasure(), charge_.panelCharge(), false);
    result.transport =
        staged.advance(staged.panelMeasure(), std::vector<double>(current.size() + 1, 0.0), current, dt);
    if (!result.transport.success)
        throw std::runtime_error("moving charge " + result.transport.stage + ": " + result.transport.reason);
    const auto nextMeasures = measures(next->panels_);
    const auto remap = staged.remap(plan.boundaries, nextMeasures);
    if (!remap.success)
        throw std::runtime_error("surface charge local remap " + remap.stage + ": " + remap.reason);
    next->charge_ = SurfaceChargeTransport::fromIntegratedCharge(nextMeasures, staged.panelCharge(), false);
    result.remapped = true;
    result.previousPanelCount = panels_.panels().size();
    result.remapConservationResidual = remap.conservationResidual;
    auto &chargeReport = result.transport;
    chargeReport.totalChargeAfter = remap.totalChargeAfter;
    chargeReport.conservationResidual = chargeReport.totalChargeAfter - chargeReport.totalChargeBefore -
                                        dt * chargeReport.integratedNormalCurrent;
    chargeReport.relativeConservationResidual =
        chargeReport.conservationScale > 0.0
            ? std::abs(chargeReport.conservationResidual) / chargeReport.conservationScale
            : 0.0;
    chargeReport.maxAbsDensity = 0.0;
    for (double value : next->charge_.density())
        chargeReport.maxAbsDensity = std::max(chargeReport.maxAbsDensity, std::abs(value));
    require(std::isfinite(chargeReport.conservationResidual) &&
                std::isfinite(chargeReport.relativeConservationResidual) &&
                chargeReport.relativeConservationResidual <= 4096.0 * std::numeric_limits<double>::epsilon(),
            "moving surface-charge update does not conserve integrated charge");
    report = std::move(result);
    return std::move(next);
}

struct SurfaceChargeMaterialTransfer::Data
{
    Data(LevelSetInterfacePanels from,LevelSetInterfacePanels to,
        std::vector<SurfaceChargeTransport::PanelCoordinate> points,
        SurfaceChargeFieldStep diagnostics,int index,double begin,double end)
        :source(std::move(from)),target(std::move(to)),boundaries(std::move(points)),
         report(std::move(diagnostics)),step(index),sourceTime(begin),targetTime(end){}
    LevelSetInterfacePanels source,target;
    std::vector<SurfaceChargeTransport::PanelCoordinate> boundaries;
    SurfaceChargeFieldStep report;
    int step;
    double sourceTime,targetTime;
};
SurfaceChargeMaterialTransfer::SurfaceChargeMaterialTransfer(std::shared_ptr<const Data> data)
    :data_(std::move(data)){}
const SurfaceChargeMaterialTransfer::Data& SurfaceChargeMaterialTransfer::checked() const
{
    if(!data_)throw std::logic_error("[SurfaceChargeMaterialTransfer] moved-from transfer");
    return *data_;
}
const LevelSetInterfacePanels& SurfaceChargeMaterialTransfer::sourcePanels() const{return checked().source;}
const LevelSetInterfacePanels& SurfaceChargeMaterialTransfer::targetPanels() const{return checked().target;}
const std::vector<SurfaceChargeTransport::PanelCoordinate>&
SurfaceChargeMaterialTransfer::materialBoundaries() const{return checked().boundaries;}
const SurfaceChargeFieldStep& SurfaceChargeMaterialTransfer::geometryReport() const{return checked().report;}
int SurfaceChargeMaterialTransfer::step() const{return checked().step;}
double SurfaceChargeMaterialTransfer::sourceTime() const{return checked().sourceTime;}
double SurfaceChargeMaterialTransfer::targetTime() const{return checked().targetTime;}
double SurfaceChargeMaterialTransfer::dt() const{return checked().targetTime-checked().sourceTime;}

SurfaceChargeMaterialTransfer SurfaceChargeField::prepareMaterialTransfer(const QuadtreeMesh& mesh,
    const LevelSetInterfacePanels& target,const std::function<Vector2D(Vector2D)>& forwardMaterialPoint,
    int step,double sourceTime,double targetTime) const
{
    require(step>0 && std::isfinite(sourceTime) && std::isfinite(targetTime)
        && sourceTime>=0 && targetTime>sourceTime && std::isfinite(targetTime-sourceTime),
        "native material transfer has no positive represented physical interval");
    require(source_ && panels_.usesNativeLeaves() && target.usesNativeLeaves()
        && target.nativeChain().source(),"native material transfer requires actual retained native panels");
    require(panels_.boundaryContacts().empty() && target.boundaryContacts().empty(),
        "native material transfer requires a separately implemented material law for declared contacts");
    require(panels_.nativeChain().symmetryPlanes()==target.nativeChain().symmetryPlanes(),
        "native material transfer changed reflection sections");
    auto plan=prepareMaterialRemap(mesh,target.nativeChain().source(),{},forwardMaterialPoint,&target);
    // The original map's coverage/order/overlap gates are checked already at
    // preparation using zero integrals; no caller Q or history is published.
    const auto targetMeasures=measures(plan.field->panels_);
    auto probe=SurfaceChargeTransport::fromIntegratedCharge(charge_.panelMeasure(),
        std::vector<double>(charge_.panelCharge().size(),0.0),false);
    const auto coverage=probe.remap(plan.boundaries,targetMeasures);
    require(coverage.success,coverage.reason.c_str());
    auto data=std::make_shared<SurfaceChargeMaterialTransfer::Data>(panels_,plan.field->panels_,
        std::move(plan.boundaries),std::move(plan.report),step,sourceTime,targetTime);
    return SurfaceChargeMaterialTransfer(std::move(data));
}

SurfaceChargeMaterialRemap SurfaceChargeField::remapPrepared(const QuadtreeMesh& mesh,
    const SurfaceChargeMaterialTransfer& transfer,const std::vector<std::vector<double>>& history) const
{
    const auto& plan=transfer.checked();
    require(source_ && panels_.usesNativeLeaves() && plan.source.usesSource(*source_),
        "native material transfer source does not own this Field geometry");
    require(!panels_.closed() && !plan.source.closed() && !plan.target.closed()
        && plan.source.axisymmetric()==panels_.axisymmetric()
        && plan.target.axisymmetric()==panels_.axisymmetric()
        && plan.source.panels().size()==panels_.panels().size()
        && plan.source.beginEndpointKind()==panels_.beginEndpointKind()
        && plan.source.endEndpointKind()==panels_.endEndpointKind(),
        "native material transfer changed source partition, topology or metric");
    for(std::size_t i=0;i<panels_.panels().size();++i) {
        const auto& a=panels_.panels()[i];const auto& b=plan.source.panels()[i];
        require(a.cellIndex==b.cellIndex && a.coordinateIsX==b.coordinateIsX
            && a.measure==b.measure && a.physicalMeasure==b.physicalMeasure
            && a.beginFrame.hasCellChartPoint && b.beginFrame.hasCellChartPoint
            && a.endFrame.hasCellChartPoint && b.endFrame.hasCellChartPoint
            && Geometry::cellChartEqual(a.beginFrame.cellChartPoint,b.beginFrame.cellChartPoint)
            && Geometry::cellChartEqual(a.endFrame.cellChartPoint,b.endFrame.cellChartPoint),
            "native material transfer source panel partition differs from charge storage");
    }
    validateNativeRemapGrid(mesh,*plan.target.nativeChain().source());
    SurfaceChargeMaterialRemap result;
    result.field=std::make_unique<SurfaceChargeField>(mesh,plan.target,0.0);
    const auto& targetMeasures=result.field->charge_.panelMeasure();
    auto current=SurfaceChargeTransport::fromIntegratedCharge(charge_.panelMeasure(),charge_.panelCharge(),false);
    result.report=plan.report;
    result.report.transport=current.remap(plan.boundaries,targetMeasures);
    require(result.report.transport.success,result.report.transport.reason.c_str());
    result.field->charge_=SurfaceChargeTransport::fromIntegratedCharge(targetMeasures,current.panelCharge(),false);
    for(const auto& values:history) {
        auto previous=SurfaceChargeTransport::fromIntegratedCharge(charge_.panelMeasure(),values,false);
        auto report=previous.remap(plan.boundaries,targetMeasures);
        require(report.success,report.reason.c_str());
        result.history.push_back(previous.panelCharge());result.historyReports.push_back(std::move(report));
    }
    result.report.remapped=true;result.report.previousPanelCount=panels_.panels().size();
    result.report.remapConservationResidual=result.report.transport.conservationResidual;
    return result;
}

SurfaceChargeMaterialRemap SurfaceChargeField::remapMoving(
    const QuadtreeMesh &mesh, std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> nextSnapshot,
    const std::vector<Vector2D> &fixedMaterialContacts,
    const std::function<Vector2D(Vector2D)> &forwardMaterialPoint,
    const std::vector<std::vector<double>> &history,
    const LevelSetInterfacePanels* preparedPanels) const
{
    // Validate all time-level arrays before tracing endpoints. Each copy is
    // staged; neither an input history nor this Field can be partially updated.
    std::vector<SurfaceChargeTransport> stagedHistory;
    stagedHistory.reserve(history.size());
    for (const auto &values : history)
        stagedHistory.push_back(
            SurfaceChargeTransport::fromIntegratedCharge(charge_.panelMeasure(), values, panels_.closed()));
    auto plan = prepareMaterialRemap(mesh,std::move(nextSnapshot),fixedMaterialContacts,forwardMaterialPoint,preparedPanels);
    const auto nextMeasures = measures(plan.field->panels_);
    auto staged = SurfaceChargeTransport::fromIntegratedCharge(charge_.panelMeasure(),charge_.panelCharge(),false);
    const auto remap = staged.remap(plan.boundaries,nextMeasures);
    if (!remap.success) throw std::runtime_error("surface charge material remap " + remap.stage + ": " + remap.reason);
    plan.field->charge_ = SurfaceChargeTransport::fromIntegratedCharge(nextMeasures,staged.panelCharge(),false);
    SurfaceChargeMaterialRemap result;
    result.report = std::move(plan.report);
    result.report.transport = remap;
    result.report.remapped = true;
    result.report.previousPanelCount = panels_.panels().size();
    result.report.remapConservationResidual = remap.conservationResidual;
    for (auto& previous : stagedHistory)
    {
        auto previousReport = previous.remap(plan.boundaries,nextMeasures);
        if (!previousReport.success)
            throw std::runtime_error("surface charge history remap " + previousReport.stage + ": " + previousReport.reason);
        result.history.push_back(previous.panelCharge());
        result.historyReports.push_back(std::move(previousReport));
    }
    result.field = std::move(plan.field);
    return result;
}

std::vector<double> SurfaceChargeField::advectiveMeasureFlux(const StaggeredVelocitySampler& velocity) const
{
    require(!panels_.usesNativeLeaves(),
        "native charge requires native wall-aware velocity fluxes, not the Cartesian staggered sampler");
    require(panels_.boundaryContacts().empty(),
        "declared-boundary surface charge requires wall-aware velocities and explicit measure fluxes; sampler extrapolation is unsupported");
    const auto& panels = panels_.panels();
    std::vector<double> result(panels.size() + (panels_.closed() ? 0 : 1), 0.0);
    const auto finitePoint = [](Vector2D p) { return std::isfinite(p.x) && std::isfinite(p.y); };
    const auto samePoint = [](Vector2D a, Vector2D b) { return a.x == b.x && a.y == b.y; };
    const std::size_t sharedCount = panels_.closed() ? panels.size() : panels.size() - 1;
    // The axis has zero physical boundary measure. Do not query either
    // staggered component at the poles, where its interpolation grid may end.
    require(panels_.closed() || (panels_.axisymmetric() && panels.front().begin.x == 0.0 && panels.back().end.x == 0.0),
        "open surface-charge endpoints must lie on the symmetry axis");
    for (std::size_t i = 0; i < sharedCount; ++i)
    {
        const auto& panel = panels[i];
        const auto& next = panels[(i + 1) % panels.size()];
        const auto& frame = panel.endFrame;
        const std::size_t face = panels_.closed() ? i : i + 1;
        // The retained graph uses one shared crossing value for both ends.
        require(finitePoint(frame.point) && finitePoint(next.beginFrame.point) &&
            samePoint(frame.point, panel.end) && samePoint(next.beginFrame.point, next.begin) &&
            samePoint(panel.end, next.begin), "surface-charge advective face is not a shared ordered endpoint");
        require(finitePoint(frame.normal) &&
            std::abs(std::hypot(frame.normal.x, frame.normal.y) - 1.0) <=
                512.0 * std::numeric_limits<double>::epsilon(),
            "surface-charge advective face normal is invalid");
        const auto u = velocity.sample(frame.point);
        require(finitePoint(u), "nonfinite sampled surface-charge velocity");
        result[face] = -u.x * frame.normal.y + u.y * frame.normal.x;
        if (panels_.axisymmetric())
        {
            const double tangentialVelocity = result[face];
            const double radialWeight = CartesianMeasure::radialWeight(frame.point.x);
            result[face] = radialWeight * tangentialVelocity;
            require(radialWeight > 0.0 && (tangentialVelocity == 0.0 || result[face] != 0.0),
                "surface-charge advective physical measure flux is unresolved");
        }
        require(std::isfinite(result[face]), "surface-charge advective measure flux overflow");
    }
    return result;
}
}
