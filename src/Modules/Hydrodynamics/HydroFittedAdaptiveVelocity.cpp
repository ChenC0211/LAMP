#include "Modules/Hydrodynamics/Detail/HydroFittedAdaptiveVelocity.h"
#include "Core/Field/HEdgeField.h"
#include "Core/Field/VEdgeField.h"
#include "Core/Mesh/CartesianLeafLocator.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/UniformLeafGrid.h"
#include "Modules/Hydrodynamics/Detail/HydroVelocityHermiteBasis.h"
#include "Modules/Hydrodynamics/HydroFittedDomain.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>

namespace HydroFittedDetail
{
namespace
{
constexpr long double twoPi = 6.283185307179586476925286766559L;
using Trace = HydroFittedVelocitySampler::VertexTrace;
using Type = HydroFittedDomain::BoundaryType;
void require(bool valid, const char *message)
{
    if (!valid)
        throw std::invalid_argument(std::string("[HydroFittedAdaptiveVelocity] ") + message);
}
double finite(long double value)
{
    const double result = static_cast<double>(value);
    require(std::isfinite(value) && std::isfinite(result) && (value == 0.0L || result != 0.0),
            "unrepresentable velocity/metric");
    return result;
}
double unit(long double value)
{
    constexpr double tolerance = 128.0 * std::numeric_limits<double>::epsilon();
    require(std::isfinite(value) && value >= -tolerance && value <= 1.0L + tolerance,
            "inverse coordinate left its retained containing interval");
    // Only inverse-coordinate roundoff after exact physical containment.
    return static_cast<double>(std::max(0.0L, std::min(1.0L, value)));
}
struct CurveJet
{
    long double value = 0.0L, derivative = 0.0L, second = 0.0L;
};
} // namespace

struct AdaptiveVelocity::Impl
{
    struct Face
    {
        HydroFittedFaceFlux flux;
        int vertices[2];
        bool noSlip;
        Face(HydroFittedFaceFlux f, int a, int b, bool wall)
            : flux(std::move(f)), vertices{a, b}, noSlip(wall)
        {
        }
    };
    struct Cell
    {
        bool active = false, axis = false;
        unsigned walls = 0;
        Vector2D lower, upper;
        double width = 0.0, height = 0.0;
        std::array<int, 4> vertices{};         // SW,SE,NE,NW.
        std::array<std::vector<int>, 4> sides; // W,E,S,N.
        long double q[4]{};
    };
    CartesianLeafLocator locator;
    bool axisymmetric;
    std::vector<Vector2D> vertices;
    std::vector<Trace> traces;
    std::vector<std::optional<Face>> faces;
    std::vector<Cell> cells;

    Impl(const QuadtreeMesh &mesh, const NativeLeafDomain &domain, const VEdgeField<double> &ux,
         const HEdgeField<double> &uy, const std::vector<Trace> &input,
         const HydroFittedDomain *fitted)
        : locator(mesh), axisymmetric(domain.axisymmetric()), traces(input)
    {
        domain.validateMeshIdentity(mesh);
        if (!fitted)
        {
            require(domain.cellIndexByRow().size() == static_cast<std::size_t>(mesh.cellIndexNum()),
                    "complete velocity requires every native leaf active");
        }
        require(!fitted || axisymmetric, "explicit fitted trace reconstruction requires RZ");
        require(input.size() == static_cast<std::size_t>(mesh.vertIndexNum()),
                "vertex trace layout mismatch");
        const int nv = mesh.vEdgeIndexNum(), nh = mesh.hEdgeIndexNum();
        require(ux.size() == nv && uy.size() == nh && ux.value && uy.value,
                "native face field layout mismatch");
        vertices.resize(mesh.vertIndexNum());
        for (int i = 0; i < mesh.vertIndexNum(); ++i)
            vertices[i] = mesh.verts(mesh.vertIndexToLabel[i]).coord();
        faces.resize(domain.nativeFaceCount());
        for (const auto native : domain.faces(mesh))
        {
            if (native.negativeCell < 0 && native.positiveCell < 0)
                continue;
            const auto &edge =
                native.horizontal ? mesh.hEdges(native.edgeLabel) : mesh.vEdges(native.edgeLabel);
            const int a = mesh.verts(edge.edgeToVert[0]).index(), b = mesh.verts(edge.edgeToVert[1]).index();
            const HydroFittedFaceFlux::Chart chart{vertices[a], native.center, vertices[b]};
            require(native.horizontal
                        ? (chart.begin.y == chart.midpoint.y && chart.midpoint.y == chart.end.y)
                        : (chart.begin.x == chart.midpoint.x && chart.midpoint.x == chart.end.x),
                    "trace-driven reconstruction requires actual straight rectangular faces");
            const int index = native.horizontal ? native.nativeRow - nv : native.nativeRow;
            require(index >= 0 && index < (native.horizontal ? nh : nv) &&
                        (native.horizontal ? uy.label(index) : ux.label(index)) == native.edgeLabel,
                    "stale native face labels");
            const double speed = native.horizontal ? uy.value[index] : ux.value[index];
            require(std::isfinite(speed), "nonfinite active native face speed");
            const auto *physical = fitted ? &fitted->faces().at(native.nativeRow) : nullptr;
            const bool noSlip = physical && physical->boundary.type == Type::NoSlipWall;
            if ((physical && physical->velocityRow < 0) || (axisymmetric && !native.horizontal && chart.begin.x == 0.0))
                require(speed == 0.0, "wall/axis normal face speed is not zero");
            for (int k : {a, b})
            {
                const auto &t = traces[k];
                require(std::isfinite(t.velocity.x) && std::isfinite(t.velocity.y) &&
                            std::isfinite(t.axialDerivative),
                        "nonfinite active vertex trace");
                if (noSlip)
                    require(t.velocity.x == 0.0 && t.velocity.y == 0.0, "no-slip vertex trace is not zero");
                if (axisymmetric && vertices[k].x == 0.0)
                    require(t.velocity.x == 0.0, "axis vertex radial speed is not zero");
            }
            faces[native.nativeRow].emplace(
                HydroFittedFaceFlux::fromEndpointVelocities(chart, native.horizontal,
                                                            static_cast<long double>(native.area) * speed,
                                                            traces[a].velocity, traces[b].velocity, axisymmetric),
                a, b, noSlip);
        }
        cells.resize(mesh.cellIndexNum());
        for (int index : domain.cellIndexByRow())
        {
            const auto &native = mesh.cells(mesh.cellIndexToLabel[index]);
            auto &c = cells[index];
            c.active = true;
            for (int k = 0; k < 4; ++k)
                c.vertices[k] = mesh.verts(native.cellToVert[k]).index();
            c.lower = vertices[c.vertices[0]];
            c.upper = vertices[c.vertices[2]];
            c.width = c.upper.x - c.lower.x;
            c.height = c.upper.y - c.lower.y;
            c.axis = axisymmetric && c.lower.x == 0.0;
            require((!axisymmetric || c.lower.x >= 0.0) && c.width > 0.0 && c.height > 0.0, "invalid active XY/RZ cell metrics");
            const long double volume =
                axisymmetric ? 0.5L * twoPi * (static_cast<long double>(c.upper.x) + c.lower.x) * c.width * c.height
                    : static_cast<long double>(c.width)*c.height;
            const double reference = domain.cellVolumes()[domain.cellRowByIndex()[index]];
            require(std::isfinite(volume) &&
                        std::abs(volume - reference) <=
                            1e-12L * std::max(volume, static_cast<long double>(reference)),
                    "actual reconstruction volume disagrees with native XY/RZ volume");
            for (int side = 0; side < 4; ++side)
            {
                const bool horizontal = side >= 2;
                const auto &labels = horizontal ? native.cellToHEdge[side % 2] : native.cellToVEdge[side % 2];
                for (int part = 0; part < 2; ++part)
                {
                    const int label = labels[part];
                    if (label >= 0)
                    {
                        const auto &edge = horizontal ? mesh.hEdges(label) : mesh.vEdges(label);
                        const int row = edge.index() + (horizontal ? nv : 0);
                        require(faces[row].has_value(), "active cell lacks a retained face");
                        c.sides[side].push_back(row);
                        c.q[side] += faces[row]->flux.nativeFlux();
                    }
                }
                require(!c.sides[side].empty(), "active cell lacks a complete side");
                const bool wall = faces[c.sides[side].front()]->noSlip;
                for (int row : c.sides[side])
                    require(faces[row]->noSlip == wall,
                            "one cell side mixes no-slip and non-wall traces");
                if (wall)
                    c.walls |= 1u << side;
            }
            require((c.walls & 3u) == 0 || (c.walls & 12u) == 0,
                    "perpendicular no-slip walls require a nonconstant divergence lifting");
        }
    }
    double radialCoordinate(const Cell &c, double r) const
    {
        const long double a = (static_cast<long double>(r) - c.lower.x) / c.width;
        const double result = unit(c.axis ? a * a : a);
        require(r == c.lower.x || result != 0.0, "unresolved positive radial coordinate");
        return result;
    }
    double coordinate(const Cell &c, int side, Vector2D p) const
    {
        return side >= 2 ? radialCoordinate(c, p.x)
                         : unit((static_cast<long double>(p.y) - c.lower.y) / c.height);
    }
    long double radialMeasure(const Cell &c, double r) const
    {
        if(!axisymmetric) return c.width;
        return c.axis ? 0.5L * twoPi * c.width * c.width : twoPi * r * c.width;
    }
    long double radialWeight(double r) const { return axisymmetric ? twoPi*r : 1.0L; }
    long double mixed(const Cell &c, Vector2D p, const Trace &trace) const
    {
        if (((c.walls & 1u) && p.x == c.lower.x) || ((c.walls & 2u) && p.x == c.upper.x))
            return c.q[3] - c.q[2];
        if (((c.walls & 4u) && p.y == c.lower.y) || ((c.walls & 8u) && p.y == c.upper.y))
            return -(c.q[1] - c.q[0]);
        return c.q[3] - c.q[2] - radialMeasure(c, p.x) * c.height * trace.axialDerivative;
    }
    long double normalDerivativeValue(const Cell &c, int side, int vertex) const
    {
        const auto p = vertices[vertex];
        const auto &trace = traces[vertex];
        const long double t = coordinate(c, side, p);
        if (side >= 2)
            return radialWeight(p.x) * c.height * trace.velocity.x - ((1.0L - t) * c.q[0] + t * c.q[1]);
        return (1.0L - t) * c.q[2] + t * c.q[3] - radialMeasure(c, p.x) * trace.velocity.y;
    }
    struct SidePosition
    {
        const Face *face = nullptr;
        long double preceding = 0.0L;
        double begin = 0, end = 1;
    };
    SidePosition sidePosition(const Cell &c, int side, double t) const
    {
        SidePosition result;
        for (int row : c.sides[side])
        {
            const auto &face = *faces[row];
            const double begin = coordinate(c, side, face.flux.chart().begin),
                         end = coordinate(c, side, face.flux.chart().end);
            if (t >= begin && t <= end)
            {
                result.face = &face;
                result.begin = begin;
                result.end = end;
                return result;
            }
            result.preceding += face.flux.nativeFlux();
        }
        throw std::invalid_argument("adaptive velocity side query has no actual native face");
    }
    template<bool Differentiate>
    CurveJet boundary(const Cell &c, int side, double t, const SidePosition &part) const
    {
        const bool horizontal = side >= 2;
        const auto &chart = part.face->flux.chart();
        const double mid = coordinate(c, side, chart.midpoint);
        const int half = t > mid ? 1 : 0;
        const auto a = half ? chart.midpoint : chart.begin, b = half ? chart.end : chart.midpoint;
        const double lo = coordinate(c, side, a), hi = coordinate(c, side, b);
        double local = 0.0;
        if (t == hi)
            local = 1.0;
        else if (t != lo)
        {
            if (horizontal && c.axis)
            {
                const long double radius = c.width * std::sqrt(static_cast<long double>(t));
                local = unit(static_cast<long double>(c.width) * c.width * (t - lo) /
                             ((radius + a.x) * (b.x - a.x)));
            }
            else
                local = unit((static_cast<long double>(t) - lo) / (hi - lo));
        }
        long double primitive = part.preceding;
        if (half)
            primitive += part.face->flux.integrate(0, 0, 1);
        primitive += part.face->flux.integrate(half, 0, local);
        const double speed = part.face->flux.normalVelocity(half, local);
        const double radius =
            horizontal ? (c.axis ? c.width * std::sqrt(t) : c.lower.x + c.width * t) : chart.begin.x;
        const long double density =
            horizontal ? radialMeasure(c, radius) * speed : radialWeight(radius) * c.height * speed;
        // The endpoint values are identically zero by the full-face flux
        // constraints, not rounded subtraction of two equal large primitives.
        const long double value = (t == 0.0 || t == 1.0)
                                      ? 0.0L
                                      : (horizontal ? c.q[side] * t - primitive : primitive - c.q[side] * t);
        long double second=0;
        if constexpr(Differentiate) {
            if(horizontal && c.axis && a.x==0.0) {
                // On the first axis half, P(l) is even and l^2=ratio*t.
                // Cancel the leading l^2 term analytically, including r=0.
                const long double ratio=static_cast<long double>(c.width)*c.width/
                    (static_cast<long double>(b.x)*b.x);
                for(const auto& term:part.face->flux.densityPolynomial(half))
                    for(int k=3;k<=5;k+=2)
                        second+=ratio*ratio*.25L*term.factor*term.weight*(k-1)*
                            term.coefficients[k]*std::pow(static_cast<long double>(local),k-3);
            } else {
                long double firstScale=1.0L/(hi-lo),secondScale=0;
                if(horizontal && c.axis) {
                    const long double r=static_cast<long double>(c.width)*std::sqrt(static_cast<long double>(t));
                    firstScale=static_cast<long double>(c.width)*c.width/(2*r*(b.x-a.x));
                    secondScale=-firstScale*static_cast<long double>(c.width)*c.width/(2*r*r);
                }
                const bool reverse=local>.5;const long double offset=reverse?1.0L-local:local;
                second=part.face->flux.densityDerivative(half,reverse,offset)*firstScale*firstScale+
                    part.face->flux.density(half,reverse,offset)*secondScale;
            }
            if(horizontal)second=-second;
        }
        return {value, horizontal ? c.q[side] - density : density - c.q[side], second};
    }
    template<bool Differentiate>
    CurveJet normalDerivative(const Cell &c, int side, double t, const SidePosition &part) const
    {
        if (c.walls & (1u << side))
            return side >= 2 ? CurveJet{-((1.0L - t) * c.q[0] + t * c.q[1]), -(c.q[1] - c.q[0])}
                             : CurveJet{(1.0L - t) * c.q[2] + t * c.q[3], c.q[3] - c.q[2]};
        const auto &chart = part.face->flux.chart();
        const bool horizontal = side >= 2, faceAxis = axisymmetric && horizontal && chart.begin.x == 0.0;
        const long double length = horizontal ? static_cast<long double>(chart.end.x) - chart.begin.x
                                              : static_cast<long double>(chart.end.y) - chart.begin.y;
        const long double radius = horizontal ? (c.axis ? c.width * std::sqrt(static_cast<long double>(t))
                                                        : c.lower.x + static_cast<long double>(c.width) * t)
                                              : chart.begin.x;
        double parameter = 0.0;
        if (t == part.end)
            parameter = 1.0;
        else if (t != part.begin)
            parameter = unit(horizontal && c.axis && !faceAxis
                                 ? static_cast<long double>(c.width) * c.width * (t - part.begin) /
                                       ((radius + chart.begin.x) * length)
                                 : (static_cast<long double>(t) - part.begin) / (part.end - part.begin));
        const Basis basis(parameter);
        CurveJet physical;
        for (int k = 0; k < 2; ++k)
        {
            const int vertex = part.face->vertices[k];
            const auto point = vertices[vertex];
            const auto &trace = traces[vertex];
            const long double mixedJet = mixed(c, point, trace);
            long double value = trace.velocity.y, slope = length * (c.q[3] - c.q[2] - mixedJet) /
                                                          (radialMeasure(c, point.x) * c.height);
            if (horizontal)
            {
                value = radialWeight(point.x) * trace.velocity.x;
                // Shared NATIVE face coordinate, not the owner's whole-cell
                // coordinate. An off-axis fine face can border an axis coarse
                // cell: cubic-in-r and cubic-in-r^2 must not give two traces.
                const long double ratio =
                    faceAxis
                        ? length * length / (static_cast<long double>(c.width) * c.width)
                        : (c.axis ? 2.0L * point.x * length / (static_cast<long double>(c.width) * c.width)
                                  : length / c.width);
                slope = (mixedJet + c.q[1] - c.q[0]) * ratio / c.height;
            }
            physical.value += basis.h[k] * value + basis.g[k] * slope;
            physical.derivative += basis.dh[k] * value + basis.dg[k] * slope;
            if constexpr(Differentiate)physical.second+=basis.ddh[k]*value+basis.ddg[k]*slope;
        }
        if (!horizontal)
            return {(1.0L - t) * c.q[2] + t * c.q[3] -
                        radialMeasure(c, static_cast<double>(radius)) * physical.value,
                    c.q[3] - c.q[2] -
                        radialMeasure(c, static_cast<double>(radius)) * physical.derivative * c.height /
                            length,
                    -radialMeasure(c, static_cast<double>(radius))*physical.second*c.height*c.height/(length*length)};
        const long double derivativeScale =
            faceAxis ? static_cast<long double>(c.width) * c.width / (length * length)
                     : (c.axis ? static_cast<long double>(c.width) * c.width / (2.0L * radius * length)
                               : c.width / length);
        const long double secondScale=horizontal && c.axis && !faceAxis
            ?-derivativeScale*static_cast<long double>(c.width)*c.width/(2*radius*radius):0;
        return {c.height * physical.value - ((1.0L - t) * c.q[0] + t * c.q[1]),
                c.height * physical.derivative * derivativeScale - (c.q[1] - c.q[0]),
                c.height*(physical.second*derivativeScale*derivativeScale+physical.derivative*secondScale)};
    }
    template<bool Differentiate>
    std::array<long double, 5> correction(const Cell &c, double u, double v) const
    {
        // Hermite Boolean sum: extend both pairs of piecewise sides, then
        // subtract their shared corner Hermite patch counted twice. This is
        // the existing streamfunction correction with nonzero side functions
        // at hanging vertices; every correction remains a curl.
        const Basis a(u), b(v);
        long double dx = 0.0L, dy = 0.0L, dxx=0, dxy=0, dyy=0;
        for (int side = 0; side < 4; ++side)
        {
            const bool horizontal = side >= 2;
            const auto position = sidePosition(c, side, horizontal ? u : v);
            const auto f = boundary<Differentiate>(c, side, horizontal ? u : v, position);
            const auto g = normalDerivative<Differentiate>(c, side, horizontal ? u : v, position);
            const int k = side % 2;
            if (horizontal)
            {
                dx += b.h[k] * f.derivative + b.g[k] * g.derivative;
                dy += b.dh[k] * f.value + b.dg[k] * g.value;
                if constexpr(Differentiate) {
                    dxx+=b.h[k]*f.second+b.g[k]*g.second;
                    dxy+=b.dh[k]*f.derivative+b.dg[k]*g.derivative;
                    dyy+=b.ddh[k]*f.value+b.ddg[k]*g.value;
                }
            }
            else
            {
                dx += a.dh[k] * f.value + a.dg[k] * g.value;
                dy += a.h[k] * f.derivative + a.g[k] * g.derivative;
                if constexpr(Differentiate) {
                    dxx+=a.ddh[k]*f.value+a.ddg[k]*g.value;
                    dxy+=a.dh[k]*f.derivative+a.dg[k]*g.derivative;
                    dyy+=a.h[k]*f.second+a.g[k]*g.second;
                }
            }
        }
        for (int y = 0; y < 2; ++y)
            for (int x = 0; x < 2; ++x)
            {
                const int vertex = c.vertices[y ? (x ? 2 : 3) : x];
                const long double cx = normalDerivativeValue(c, x, vertex),
                                  cy = normalDerivativeValue(c, 2 + y, vertex);
                const long double cxy = mixed(c, vertices[vertex], traces[vertex]);
                dx -= cx * a.dg[x] * b.h[y] + cy * a.dh[x] * b.g[y] + cxy * a.dg[x] * b.g[y];
                dy -= cx * a.g[x] * b.dh[y] + cy * a.h[x] * b.dg[y] + cxy * a.g[x] * b.dg[y];
                if constexpr(Differentiate) {
                    dxx-=cx*a.ddg[x]*b.h[y]+cy*a.ddh[x]*b.g[y]+cxy*a.ddg[x]*b.g[y];
                    dxy-=cx*a.dg[x]*b.dh[y]+cy*a.dh[x]*b.dg[y]+cxy*a.dg[x]*b.dg[y];
                    dyy-=cx*a.g[x]*b.ddh[y]+cy*a.h[x]*b.ddg[y]+cxy*a.g[x]*b.ddg[y];
                }
            }
        return {{dx,dy,dxx,dxy,dyy}};
    }
    const Cell &selected(Vector2D p) const
    {
        require(!axisymmetric || p.x >= 0.0, "negative RZ query radius");
        const auto incident = locator.incident(p);
        for (int k = 0; k < incident.count; ++k)
            if (cells[incident.indices[k]].active)
                return cells[incident.indices[k]];
        throw std::invalid_argument("adaptive velocity query has no actual active incident cell");
    }
    Vector2D sample(Vector2D p) const
    {
        const auto &c = selected(p);
        const double u = radialCoordinate(c, p.x),
                     v = unit((static_cast<long double>(p.y) - c.lower.y) / c.height);
        if (((c.walls & 1u) && p.x == c.lower.x) || ((c.walls & 2u) && p.x == c.upper.x) ||
            ((c.walls & 4u) && p.y == c.lower.y) || ((c.walls & 8u) && p.y == c.upper.y))
            return {0.0, 0.0};
        const auto psi = correction<false>(c, u, v);
        const long double radial = (1.0L - u) * c.q[0] + u * c.q[1] + psi[1];
        const long double axial = (1.0L - v) * c.q[2] + v * c.q[3] - psi[0];
        return {axisymmetric && p.x == 0.0 ? 0.0 : finite(radial / (radialWeight(p.x) * c.height)),
                finite(axial / radialMeasure(c, p.x))};
    }
    HydroFittedVelocitySampler::FirstJet firstJet(Vector2D p,bool volumeCoordinates=false) const
    {
        const auto& c=selected(p);const double u=radialCoordinate(c,p.x),v=unit((static_cast<long double>(p.y)-c.lower.y)/c.height);
        const auto psi=correction<true>(c,u,v);const long double dqX=c.q[1]-c.q[0],dqY=c.q[3]-c.q[2];
        const long double radial=(1.0L-u)*c.q[0]+u*c.q[1]+psi[1],axial=(1.0L-v)*c.q[2]+v*c.q[3]-psi[0];
        HydroFittedVelocitySampler::FirstJet out;out.nativeCell=static_cast<int>(&c-cells.data());
        out.value={axisymmetric && p.x==0?0:finite(radial/(radialWeight(p.x)*c.height)),finite(axial/radialMeasure(c,p.x))};
        if(axisymmetric && volumeCoordinates) {
            out.value.x=finite(radial/(twoPi*c.height));
            if(c.axis) {
                out.dx={finite((dqX+psi[3])/(radialMeasure(c,p.x)*c.height)),
                    finite(-2*psi[2]/(radialMeasure(c,p.x)*c.width*c.width))};
            } else {
                out.dx={finite((dqX+psi[3])/(radialMeasure(c,p.x)*c.height)),
                    finite(-psi[2]/(radialMeasure(c,p.x)*c.width*p.x)-out.value.y/(p.x*p.x))};
            }
            out.dy={finite(psi[4]/(twoPi*c.height*c.height)),finite((dqY-psi[3])/(radialMeasure(c,p.x)*c.height))};
        } else if(!axisymmetric) {
            out.dx={finite((dqX+psi[3])/(c.width*c.height)),finite(-psi[2]/(c.width*c.width))};
            out.dy={finite(psi[4]/(c.height*c.height)),finite((dqY-psi[3])/(c.width*c.height))};
        } else if(c.axis && p.x==0) {
            out.dx={finite((dqX+psi[3])/(twoPi*c.width*c.width*c.height)),0};
            out.dy={0,finite((dqY-psi[3])/(radialMeasure(c,p.x)*c.height))};
        } else {
            const long double coordinateDerivative=c.axis?2.0L*p.x/(c.width*c.width):1.0L/c.width;
            out.dx={finite((dqX+psi[3])*coordinateDerivative/(radialWeight(p.x)*c.height)-out.value.x/p.x),
                finite(-psi[2]*coordinateDerivative/radialMeasure(c,p.x)-(c.axis?0:out.value.y/p.x))};
            out.dy={finite(psi[4]/(radialWeight(p.x)*c.height*c.height)),finite((dqY-psi[3])/(radialMeasure(c,p.x)*c.height))};
        }
        out.oneSided=u==0 || u==1 || v==0 || v==1;
        for(int side=0;side<4;++side){const double t=side>=2?u:v;const auto position=sidePosition(c,side,t);
            out.oneSided=out.oneSided || t==position.begin || t==position.end || t==coordinate(c,side,position.face->flux.chart().midpoint);}
        if(((c.walls&1u) && p.x==c.lower.x) || ((c.walls&2u) && p.x==c.upper.x) ||
            ((c.walls&4u) && p.y==c.lower.y) || ((c.walls&8u) && p.y==c.upper.y))out.value={0,0};
        return out;
    }
    double divergence(Vector2D p) const
    {
        const auto &c = selected(p);
        return finite(((c.q[1] - c.q[0]) + (c.q[3] - c.q[2])) / (radialMeasure(c, p.x) * c.height));
    }
};

AdaptiveVelocity::AdaptiveVelocity(const QuadtreeMesh &mesh, const HydroFittedDomain &domain,
                                   const VEdgeField<double> &ux, const HEdgeField<double> &uy,
                                   const std::vector<Trace> &traces)
    : impl_(std::make_unique<const Impl>(mesh, domain.geometry(), ux, uy, traces, &domain))
{
}
AdaptiveVelocity::AdaptiveVelocity(const QuadtreeMesh &mesh, const NativeLeafDomain &domain,
                                   const VEdgeField<double> &ux, const HEdgeField<double> &uy,
                                   const std::vector<Trace> &traces)
    : impl_(std::make_unique<const Impl>(mesh, domain, ux, uy, traces, nullptr))
{
}
AdaptiveVelocity::~AdaptiveVelocity() = default;
Vector2D AdaptiveVelocity::sample(Vector2D p) const
{
    return impl_->sample(p);
}
double AdaptiveVelocity::divergence(Vector2D p) const
{
    return impl_->divergence(p);
}
HydroFittedVelocitySampler::FirstJet AdaptiveVelocity::firstJet(Vector2D p) const
{
    return impl_->firstJet(p);
}
HydroFittedVelocitySampler::FirstJet AdaptiveVelocity::volumeFirstJet(Vector2D p) const
{
    return impl_->firstJet(p,true);
}
const HydroFittedFaceFlux &AdaptiveVelocity::face(int native) const
{
    require(native >= 0 && static_cast<std::size_t>(native) < impl_->faces.size() && impl_->faces[native],
            "native face has no retained active flux trace");
    return impl_->faces[native]->flux;
}
} // namespace HydroFittedDetail
