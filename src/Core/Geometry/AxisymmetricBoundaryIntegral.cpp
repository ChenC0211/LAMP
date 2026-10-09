#include "Core/Geometry/AxisymmetricBoundaryIntegral.h"
#include "Core/Geometry/CartesianBoundaryIntegral.h"
#include "Core/Geometry/LevelSetMappedSegment.h"
#include "Core/Mesh/CartesianMeasure.h"

#include <cmath>
#include <stdexcept>

namespace CartesianBoundaryIntegral
{
long double yVolume(const LevelSetHermiteGeometry::ParametricPiece& segment,
    bool axisymmetric, double axialReference)
{
    const auto a = segment.anchorBegin, b = segment.anchorEnd;
    const auto unit = [](double t) { return std::isfinite(t) && t >= 0.0 && t <= 1.0; };
    if (!std::isfinite(a.x) || !std::isfinite(a.y) || !std::isfinite(b.x)
        || !std::isfinite(b.y) || (axisymmetric && (a.x < 0.0 || b.x < 0.0))
        || !std::isfinite(axialReference) || !unit(segment.tBegin) || !unit(segment.tEnd))
        throw std::invalid_argument("[AxisymmetricBoundaryIntegral] invalid segment, radius or parameter");
    const long double dr = static_cast<long double>(b.x)-a.x;
    const long double dz = static_cast<long double>(b.y)-a.y;
    const long double dt = static_cast<long double>(segment.tEnd)-segment.tBegin;
    const long double z0 = static_cast<long double>(a.y)-axialReference;
    if (!std::isfinite(dr) || !std::isfinite(dz) || !std::isfinite(z0))
        throw std::overflow_error("[AxisymmetricBoundaryIntegral] unrepresentable anchor differences");
    if (dr == 0.0L || dt == 0.0L) return 0.0L;
    // Translate to the interval's beginning, retaining its separate width.
    // fma also retains a sub-global-ULP offset from the axial reference.
    const long double r = std::fma(dr, static_cast<long double>(segment.tBegin), static_cast<long double>(a.x));
    const long double z = std::fma(dz, static_cast<long double>(segment.tBegin), z0);
    const long double average = axisymmetric ? r*z + 0.5L*dt*(dr*z+dz*r) + (dt*dt/3.0L)*dr*dz : z+0.5L*dt*dz;
    const long double result = -(axisymmetric ? static_cast<long double>(CartesianMeasure::twoPi) : 1.0L)*dr*dt*average;
    if (!std::isfinite(r) || !std::isfinite(z) || !std::isfinite(average) || !std::isfinite(result)
        || (average != 0.0L && result == 0.0L))
        throw std::overflow_error("[AxisymmetricBoundaryIntegral] unrepresentable polynomial moment");
    return result;
}

long double yVolume(const LevelSetMappedSegment& segment, bool axisymmetric, double axialReference)
{
    const auto a = segment.begin(), b = segment.end();
    if (!std::isfinite(a.x) || !std::isfinite(a.y) || !std::isfinite(b.x)
        || !std::isfinite(b.y) || (axisymmetric && (a.x < 0.0 || b.x < 0.0)) || !std::isfinite(axialReference))
        throw std::invalid_argument("[AxisymmetricBoundaryIntegral] invalid mapped segment or radius");
    const long double dr = static_cast<long double>(b.x)-a.x;
    const long double dz = static_cast<long double>(b.y)-a.y;
    const long double za = static_cast<long double>(a.y)-axialReference;
    const long double zb = static_cast<long double>(b.y)-axialReference;
    if (!std::isfinite(dr) || !std::isfinite(dz) || !std::isfinite(za) || !std::isfinite(zb))
        throw std::overflow_error("[AxisymmetricBoundaryIntegral] unrepresentable mapped anchor differences");
    long double sum = 0.0L, correction = 0.0L;
    for (const auto& range : segment.ranges())
    {
        // Domain coverage belongs to the retained mapped segment. Its offset
        // and width are separately rounded terminal values, so do not recreate
        // an endpoint (even for an independent coverage comparison) here.
        if (!std::isfinite(range.offset) || !std::isfinite(range.length)
            || range.offset < 0.0L || range.offset > 1.0L
            || range.length <= 0.0L || range.length > 1.0L)
            throw std::invalid_argument("[AxisymmetricBoundaryIntegral] invalid mapped range");
        if (!range.negative || dr == 0.0L) continue;
        const long double chartDr = range.fromEnd ? -dr : dr;
        const long double chartDz = range.fromEnd ? -dz : dz;
        const long double r = std::fma(chartDr, range.offset,
            static_cast<long double>(range.fromEnd ? b.x : a.x));
        const long double z = std::fma(chartDz, range.offset, range.fromEnd ? zb : za);
        const long double width = range.length;
        const long double average = axisymmetric ? r*z + 0.5L*width*(chartDr*z+chartDz*r)
            + (width*width/3.0L)*chartDr*chartDz : z+0.5L*width*chartDz;
        // dr remains in the original input direction, also in a fromEnd chart.
        const long double moment = -(axisymmetric ? static_cast<long double>(CartesianMeasure::twoPi) : 1.0L)*dr*width*average;
        if (!std::isfinite(r) || !std::isfinite(z) || !std::isfinite(average) || !std::isfinite(moment)
            || (average != 0.0L && moment == 0.0L))
            throw std::overflow_error("[AxisymmetricBoundaryIntegral] unrepresentable mapped polynomial moment");
        const long double corrected = moment-correction, next = sum+corrected;
        if (!std::isfinite(corrected) || !std::isfinite(next))
            throw std::overflow_error("[AxisymmetricBoundaryIntegral] nonfinite mapped moment sum");
        correction = (next-sum)-corrected;
        sum = next;
    }
    return sum;
}
}

namespace AxisymmetricBoundaryIntegral
{
long double axialVolume(const LevelSetHermiteGeometry::ParametricPiece& segment, double reference)
{ return CartesianBoundaryIntegral::yVolume(segment,true,reference); }
long double axialVolume(const LevelSetMappedSegment& segment, double reference)
{ return CartesianBoundaryIntegral::yVolume(segment,true,reference); }
}
