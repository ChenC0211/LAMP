#include "Modules/Hydrodynamics/HydroFittedFaceFlux.h"
#include "Core/Geometry/LevelSetMappedSegment.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <initializer_list>
#include <limits>
#include <stdexcept>
#include <string>
#include <sstream>
#include <iomanip>

namespace
{
void require(bool condition, const char *reason)
{
    if (!condition)
        throw std::invalid_argument(std::string("[HydroFittedFaceFlux] ") + reason);
}
bool finite(Vector2D point)
{
    return std::isfinite(point.x) && std::isfinite(point.y);
}
// P(w) is the actual half-face positive-normal flux primitive measured from
// the selected endpoint. In a reversed chart P(w)=F(end)-F(end-w), so P'(w)
// keeps the original +r/+z flux sign. Keep these coefficient expressions shared
// with the legacy anchored integral, including the axis xi^2 basis.
std::array<long double, 7> fluxPrimitiveCoefficients(bool radialPolynomial, int half, bool reverseHalf,
                                                     long double q, long double d0, long double d1,
                                                     long double sign)
{
    std::array<long double, 7> p{};
    const int degree = radialPolynomial ? 6 : 3;
    if (!radialPolynomial)
    {
        if (half == 0 && !reverseHalf)
        {
            p[1] = 0.5L * (q + sign * d0);
            p[2] = -0.25L * sign * (2.0L * d0 + d1);
            p[3] = 0.125L * sign * (d0 + d1);
        }
        else if (half == 1 && reverseHalf)
        {
            p[1] = 0.5L * (q + sign * d1);
            p[2] = -0.25L * sign * (2.0L * d1 + d0);
            p[3] = 0.125L * sign * (d0 + d1);
        }
        else
        {
            p[1] = 0.5L * q - 0.125L * sign * (d0 + d1);
            p[2] = 0.125L * sign * (d1 - d0);
            p[3] = 0.125L * sign * (d0 + d1);
        }
    }
    else
    {
        // F(xi)=A*xi^2+B*xi^4+C*xi^6, because the axis
        // streamfunction is cubic in xi^2, not cubic in xi.
        const long double a = q - d0, b = 2.0L * d0 + d1, c = -d0 - d1;
        if (half == 0 && !reverseHalf)
        {
            p[2] = a / 4.0L;
            p[4] = b / 16.0L;
            p[6] = c / 64.0L;
        }
        else if (half == 1 && reverseHalf)
        {
            // First reverse in s=xi^2, then use 1-s=2*v-v^2,
            // v=(1-halfParameter)/2. Factoring at xi=1 retains
            // an exactly zero endpoint flux without large-value subtraction.
            const long double ae = q - d1, be = 2.0L * d1 + d0;
            p[1] = ae;
            p[2] = -ae / 4.0L + be;
            p[3] = -be / 2.0L + c;
            p[4] = (be - 12.0L * c) / 16.0L;
            p[5] = 3.0L * c / 16.0L;
            p[6] = -c / 64.0L;
        }
        else
        {
            p[1] = a / 2.0L + b / 4.0L + 3.0L * c / 32.0L;
            p[2] = a / 4.0L + 3.0L * b / 8.0L + 15.0L * c / 64.0L;
            p[3] = b / 4.0L + 5.0L * c / 16.0L;
            p[4] = b / 16.0L + 15.0L * c / 64.0L;
            p[5] = 3.0L * c / 32.0L;
            p[6] = c / 64.0L;
        }
    }
    if (half == 0 && reverseHalf)
        for (int k = 2; k <= degree; k += 2)
            p[k] = -p[k];
    for (long double coefficient : p)
        require(std::isfinite(coefficient), "nonfinite flux primitive coefficient");
    return p;
}

} // namespace

HydroFittedFaceFlux::HydroFittedFaceFlux(Chart chart, bool isHorizontal, bool radial, long double flux,
                                         long double derivativeBegin, long double derivativeEnd, bool axisymmetric)
    : chart_(chart), radialPolynomial(radial), axisymmetric_(axisymmetric), q(flux), d0(derivativeBegin),
      d1(derivativeEnd), sign(isHorizontal ? -1.0L : 1.0L)
{
    require(finite(chart.begin) && finite(chart.midpoint) && finite(chart.end)
                && (!axisymmetric || (chart.begin.x >= 0.0 && chart.midpoint.x >= 0.0 && chart.end.x >= 0.0)),
            "invalid retained XY/RZ face coordinates");
    require(isHorizontal ? (chart.begin.x < chart.midpoint.x && chart.midpoint.x < chart.end.x)
                         : (chart.begin.y < chart.midpoint.y && chart.midpoint.y < chart.end.y),
            "face chart is not ordered in the positive direction");
    require(!radial || (axisymmetric && isHorizontal && chart.begin.x == 0.0),
            "radial polynomial requires a horizontal face incident on the axis");
    require(std::isfinite(q) && std::isfinite(d0) && std::isfinite(d1), "nonfinite retained face flux/jet");
    if (axisymmetric && !isHorizontal && chart.begin.x == 0.0 && chart.midpoint.x == 0.0 && chart.end.x == 0.0)
        require(q == 0.0L && d0 == 0.0L && d1 == 0.0L, "radial axis cannot carry normal flux");
}

HydroFittedFaceFlux HydroFittedFaceFlux::fromEndpointVelocities(Chart chart, bool horizontal, long double q,
                                                                Vector2D beginVelocity, Vector2D endVelocity, bool axisymmetric)
{
    require(finite(beginVelocity) && finite(endVelocity), "nonfinite endpoint velocity");
    if (axisymmetric && !horizontal && chart.begin.x == 0.0 && chart.midpoint.x == 0.0 && chart.end.x == 0.0)
        require(beginVelocity.x == 0.0 && endVelocity.x == 0.0, "radial axis requires zero normal velocity");
    constexpr long double twoPi = 6.283185307179586476925286766559L;
    const Vector2D lo = (chart.midpoint - chart.begin) * 2.0, hi = (chart.end - chart.midpoint) * 2.0;
    const bool radial = axisymmetric && horizontal && chart.begin.x == 0.0;
    const auto density = [&](Vector2D p, Vector2D tangent, Vector2D velocity) {
        const long double normal = horizontal ? static_cast<long double>(tangent.x) * velocity.y -
                                                    static_cast<long double>(tangent.y) * velocity.x
                                              : static_cast<long double>(tangent.y) * velocity.x -
                                                    static_cast<long double>(tangent.x) * velocity.y;
        return axisymmetric ? twoPi * p.x * normal : normal;
    };
    long double first = density(chart.begin, lo, beginVelocity), last = density(chart.end, hi, endVelocity);
    if (radial)
    {
        require(beginVelocity.x == 0.0 && lo.y == 0.0,
                "axis endpoint requires radial symmetry and orthogonal chart");
        // Analytic xi^2 limit; never divide by r at the axis.
        first = 0.5L * twoPi * lo.x * lo.x * beginVelocity.y;
        last *= 0.5L;
    }
    return HydroFittedFaceFlux(chart, horizontal, radial, q, horizontal ? q - first : first - q,
                               horizontal ? q - last : last - q, axisymmetric);
}

Quantum::Geometry::UnitPolynomialRoots HydroFittedFaceFlux::densityRoots(int half) const
{
    try {return Quantum::Geometry::isolateUnitPolynomial(densityPolynomial(half));}
    catch(const std::exception& error) {
        std::ostringstream message;message<<std::setprecision(21)<<error.what()<<"; physical face flux half="<<half
            <<" q="<<q<<" d0="<<d0<<" d1="<<d1<<" radial="<<radialPolynomial<<" sign="<<sign
            <<" begin=("<<chart_.begin.x<<','<<chart_.begin.y<<") end=("<<chart_.end.x<<','<<chart_.end.y<<')';
        throw std::runtime_error(message.str());
    }
}
std::vector<Quantum::Geometry::UnitPolynomialTerm> HydroFittedFaceFlux::densityPolynomial(int half) const
{
    require(half==0 || half==1,"invalid flux root half-face");
    std::vector<Quantum::Geometry::UnitPolynomialTerm> terms;
    for (int component=0;component<3;++component) {
        const auto p=fluxPrimitiveCoefficients(radialPolynomial,half,false,component==0?1:0,
            component==1?1:0,component==2?1:0,sign);
        Quantum::Geometry::UnitPolynomialTerm term;term.factor=component==0?q:component==1?d0:d1;
        for (int k=1;k<=6;++k) term.coefficients[k-1]=static_cast<double>(k*p[k]);
        terms.push_back(term);
    }
    return terms;
}
long double HydroFittedFaceFlux::density(int half, bool fromEnd, long double offset) const
{
    require((half == 0 || half == 1) && std::isfinite(offset) && offset >= 0 && offset <= 1,
            "invalid flux-density half-face coordinate");
    const auto p = fluxPrimitiveCoefficients(radialPolynomial, half, fromEnd, q, d0, d1, sign);
    const int degree = radialPolynomial ? 6 : 3;
    long double value = degree * p[degree];
    for (int k = degree - 1; k >= 1; --k)
        value = value * offset + k * p[k];
    require(std::isfinite(value), "nonfinite normal-flux density");
    return value;
}

long double HydroFittedFaceFlux::densityDerivative(int half, bool fromEnd, long double offset) const
{
    require((half==0 || half==1) && std::isfinite(offset) && offset>=0 && offset<=1,
            "invalid flux-density derivative coordinate");
    const auto p=fluxPrimitiveCoefficients(radialPolynomial,half,fromEnd,q,d0,d1,sign);
    const int degree=radialPolynomial?6:3;long double value=degree*(degree-1)*p[degree];
    for(int k=degree-1;k>=2;--k)value=value*offset+k*(k-1)*p[k];
    require(std::isfinite(value),"nonfinite flux-density derivative");
    return fromEnd?-value:value;
}

double HydroFittedFaceFlux::normalVelocity(int half, double parameter) const
{
    require((half==0 || half==1) && std::isfinite(parameter) && parameter>=0.0 && parameter<=1.0,
        "invalid normal-velocity face coordinate");
    const auto a=half==0?chart_.begin:chart_.midpoint;
    const auto b=half==0?chart_.midpoint:chart_.end;
    const long double dx=static_cast<long double>(b.x)-a.x,dy=static_cast<long double>(b.y)-a.y;
    const long double length=std::hypot(dx,dy),t=parameter;
    require(std::isfinite(length) && length>0.0L,"unresolved physical half-face length");
    long double velocity=0.0L;
    if(radialPolynomial && half==0)
    {
        // P'(t)/r: cancel t analytically before evaluation, including t=0.
        const auto p=fluxPrimitiveCoefficients(true,0,false,q,d0,d1,sign);
        const long double t2=t*t;
        velocity=((6.0L*p[6]*t2+4.0L*p[4])*t2+2.0L*p[2])
            /(6.283185307179586476925286766559L*dx*length);
    }
    else
    {
        const bool reverse=t>0.5L;
        const auto p=fluxPrimitiveCoefficients(radialPolynomial,half,reverse,q,d0,d1,sign);
        const int degree=radialPolynomial?6:3;
        const long double w=reverse?1.0L-t:t;
        long double density=degree*p[degree];
        for(int k=degree-1;k>=1;--k) density=density*w+k*p[k];
        const long double radius=(1.0L-t)*a.x+t*b.x;
        if(!axisymmetric_) velocity=density/length;
        else if(radius==0.0L)
            require(q==0.0L && d0==0.0L && d1==0.0L,"singular nonzero axis face density");
        else velocity=density/(6.283185307179586476925286766559L*radius*length);
    }
    const double result=static_cast<double>(velocity);
    require(std::isfinite(velocity) && std::isfinite(result) && (velocity==0.0L || result!=0.0),
        "unrepresentable normal face velocity");
    return result;
}

long double HydroFittedFaceFlux::integrate(int half, double anchor0, double anchor1, double sBegin,
                                           double sEnd) const
{
    require(half == 0 || half == 1, "flux half-edge must be zero or one");
    for (double t : {anchor0, anchor1, sBegin, sEnd})
        require(std::isfinite(t) && t >= 0.0 && t <= 1.0,
                "flux parameters must lie in the closed unit chart");
    if (anchor0 == anchor1 || sBegin == sEnd)
        return 0.0L;
    // Keep the small independent parameter near zero where possible.
    // This is a chart reparameterization, not movement of a physical point.
    long double a0 = anchor0, a1 = anchor1, s0 = sBegin, s1 = sEnd;
    if (s0 + s1 > 1.0L)
    {
        std::swap(a0, a1);
        s0 = 1.0L - s0;
        s1 = 1.0L - s1;
    }
    const long double width = a1 - a0;
    const bool reverseHalf = a0 + width * (0.5L * (s0 + s1)) > 0.5L;
    // P is the half-edge flux primitive, expanded at its nearer endpoint.
    // For a reversed half chart P(w)=F(end)-F(end-w), so P'(w) is still
    // positive-normal flux. Constant primitive values never get subtracted.
    const int degree = radialPolynomial ? 6 : 3;
    const auto p = fluxPrimitiveCoefficients(radialPolynomial, half, reverseHalf, q, d0, d1, sign);

    const long double origin = reverseHalf ? 1.0L - a0 : a0;
    const long double slope = reverseHalf ? -width : width;
    // Analytically compose P'(origin+slope*s), rather than rounding the
    // two composed t endpoints. The integration scale keeps both interval
    // factors independent, including sub-ULP anchored phase intervals.
    std::array<long double, 6> density{};
    density[0] = degree * p[degree];
    int currentDegree = 0;
    for (int k = degree - 2; k >= 0; --k)
    {
        for (int j = currentDegree + 1; j > 0; --j)
            density[j] = origin * density[j] + slope * density[j - 1];
        density[0] = origin * density[0] + (k + 1) * p[k + 1];
        ++currentDegree;
    }
    for (long double coefficient : density)
        require(std::isfinite(coefficient), "nonfinite anchored flux polynomial");
    // Factored divided difference of the composed antiderivative. This
    // never computes P(tEnd)-P(tBegin), even when each P is large or its
    // global parameter endpoints would coincide after rounding.
    long double atBegin = density[degree - 1] / degree, divided = 0.0L;
    for (int k = degree - 1; k >= 0; --k)
    {
        divided = divided * s1 + atBegin;
        atBegin = atBegin * s0 + (k ? density[k - 1] / k : 0.0L);
    }
    const long double parameterMeasure = width * (s1 - s0);
    require(std::isfinite(parameterMeasure) && parameterMeasure != 0.0L,
            "anchored flux interval is unrepresentable");
    const long double result = parameterMeasure * divided;
    // Windows long double has binary64 precision. Near a flux-density zero,
    // independent half-chart compositions can otherwise lose partition
    // additivity through cancellation. Reuse Core's exact polynomial basis
    // only for this poorly conditioned evaluation; the physical polynomial,
    // roots and integration endpoints are unchanged.
    if (std::abs(divided)<std::max({std::abs(q),std::abs(d0),std::abs(d1)})/32)
        return Quantum::Geometry::integrateUnitPolynomial(densityPolynomial(half),anchor0,anchor1,sBegin,sEnd);
    require(std::isfinite(divided) && std::isfinite(result) && (divided == 0.0L || result != 0.0L),
            "integrated face flux is unrepresentable");
    return result;
}

long double HydroFittedFaceFlux::integrate(int half, const LevelSetMappedSegment &segment) const
{
    require(half == 0 || half == 1, "flux half-edge must be zero or one");
    const auto &chart = chart_;
    const auto a = segment.begin(), b = segment.end();
    const auto expectedBegin = half == 0 ? chart.begin : chart.midpoint;
    const auto expectedEnd = half == 0 ? chart.midpoint : chart.end;
    require(a.x == expectedBegin.x && a.y == expectedBegin.y && b.x == expectedEnd.x && b.y == expectedEnd.y,
            "mapped segment differs from the actual positive-direction flux half-face");

    const int degree = radialPolynomial ? 6 : 3;
    long double sum = 0.0L, correction = 0.0L;
    for (const auto &range : segment.ranges())
    {
        // The source retains exact domain coverage. Check its terminal
        // floating values without reconstructing offset+length endpoints.
        require(std::isfinite(range.offset) && std::isfinite(range.length) && range.offset >= 0.0L &&
                    range.offset <= 1.0L && range.length > 0.0L && range.length <= 1.0L,
                "invalid mapped flux range");
        if (!range.negative)
            continue;
        const auto p = fluxPrimitiveCoefficients(radialPolynomial, half, range.fromEnd, q, d0, d1, sign);
        // Translate P'(offset+s) into powers of s. In a fromEnd chart its
        // derivative already has the original positive-normal flux sign.
        std::array<long double, 6> density{};
        density[0] = degree * p[degree];
        int currentDegree = 0;
        for (int k = degree - 2; k >= 0; --k)
        {
            for (int j = currentDegree + 1; j > 0; --j)
                density[j] = range.offset * density[j] + density[j - 1];
            density[0] = range.offset * density[0] + (k + 1) * p[k + 1];
            ++currentDegree;
        }
        for (long double coefficient : density)
            require(std::isfinite(coefficient), "nonfinite translated mapped flux polynomial");
        long double average = density[degree - 1] / degree;
        for (int k = degree - 2; k >= 0; --k)
            average = average * range.length + density[k] / (k + 1);
        const long double value = range.length * average;
        require(std::isfinite(average) && std::isfinite(value) && (average == 0.0L || value != 0.0L),
                "mapped face flux is unrepresentable");
        const long double corrected = value - correction, next = sum + corrected;
        require(std::isfinite(corrected) && std::isfinite(next), "nonfinite mapped face flux sum");
        correction = (next - sum) - corrected;
        sum = next;
    }
    return sum;
}
