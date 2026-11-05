#pragma once

#include "Core/Geometry/Detail/LevelSetHermiteCell.h"
#include <algorithm>
#include <cmath>
#include <limits>

namespace LAMP::Geometry::Detail
{
namespace
{
namespace HermiteDifferential
{
// Original differential geometry of one retained scalar cell. No mesh
// selection, distance assumption, reference shape, projection or RZ pole law.
// The caller supplies the original checked terminal conversion so diagnostic
// prefixes and success-path arithmetic remain unchanged in the legacy engine.
inline bool usableGradient(const LevelSetHermiteScalar::Sample &s, long double gradient,
                           long double localVariationScale, Vector2D h)
{
    const long double derivativeScale = std::max(
        {gradient, localVariationScale,
         std::max({std::abs(static_cast<long double>(s.dxx)), std::abs(static_cast<long double>(s.dxy)),
                   std::abs(static_cast<long double>(s.dyy))}) *
             std::max(h.x, h.y)});
    return std::isfinite(gradient) && std::isfinite(derivativeScale) &&
           gradient > 512.0L * std::numeric_limits<double>::epsilon() * derivativeScale;
}
template <class Frame, class Checked>
Frame frame(Vector2D point, const LevelSetHermiteScalar::Sample &s, long double gradient, Checked checked)
{
    Frame result;
    result.point = point;
    result.normal = Vector2D(checked(s.dx / gradient, "non-finite normal x"),
                             checked(s.dy / gradient, "non-finite normal y"));
    const long double nx = result.normal.x, ny = result.normal.y;
    result.curvature = checked((ny * ny * s.dxx - 2.0L * nx * ny * s.dxy + nx * nx * s.dyy) / gradient,
                               "non-finite Hermite curvature");
    result.levelSetResidual =
        checked(std::abs(static_cast<long double>(s.value)) / gradient, "non-finite level-set residual");
    return result;
}
template <class Checked>
double curvatureDerivative(const LevelSetHermiteScalar::Sample &s,
                           const LevelSetHermiteScalar::ThirdDerivatives &third, long double gradient,
                           Checked checked)
{
    const long double nx = s.dx / gradient, ny = s.dy / gradient;
    const long double tx = -ny, ty = nx;
    const long double htt = tx * tx * s.dxx + 2.0L * tx * ty * s.dxy + ty * ty * s.dyy;
    const long double hnt = nx * (s.dxx * tx + s.dxy * ty) + ny * (s.dxy * tx + s.dyy * ty);
    const long double ttt = third.dxxx * tx * tx * tx + 3.0L * third.dxxy * tx * tx * ty +
                            3.0L * third.dxyy * tx * ty * ty + third.dyyy * ty * ty * ty;
    // g=|grad(phi)|, kappa=t.H.t/g. No signed-distance assumption:
    // D_t kappa = T[t,t,t]/g - 3*kappa*(n.H.t)/g.
    return checked(ttt / gradient - 3.0L * (htt / gradient) * (hnt / gradient),
                   "non-finite Hermite tangential curvature derivative");
}
} // namespace HermiteDifferential
} // namespace
} // namespace LAMP::Geometry::Detail
