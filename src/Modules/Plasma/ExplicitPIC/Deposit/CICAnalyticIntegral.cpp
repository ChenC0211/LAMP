#include "Modules/Plasma/ExplicitPIC/Deposit/CICAnalyticIntegral.h"

#include <algorithm>
#include <cmath>

namespace CICAnalytic {

namespace {
// Antiderivatives for the two linear branches on [-1,0] and [0,1].
// For t in [-1,0]: s(t) = 1 + t.
// For t in [0,1]:  s(t) = 1 - t.
inline double Fneg(const double t) { return t + 0.5 * t * t; }
inline double Fpos(const double t) { return t - 0.5 * t * t; }
} // namespace

double integrateHatNormalized(double u0, double u1)
{
    if (!std::isfinite(u0) || !std::isfinite(u1))
        return 0.0;

    // Ensure u0 <= u1
    double u = std::min(u0, u1);
    double v = std::max(u0, u1);
    if (u == v)
        return 0.0;

    // No overlap with support [-1,1]
    if (v <= -1.0 || u >= 1.0)
        return 0.0;

    // Clamp to support
    u = std::max(u, -1.0);
    v = std::min(v, 1.0);

    // Integrate piecewise.
    // On [-1,0]: s(t)=1+t. On [0,1]: s(t)=1-t.
    if (v <= 0.0)
    {
        return Fneg(v) - Fneg(u);
    }
    if (u >= 0.0)
    {
        return Fpos(v) - Fpos(u);
    }

    // Split across zero.
    return (Fneg(0.0) - Fneg(u)) + (Fpos(v) - Fpos(0.0));
}

double integrate1D(double a, double b, double xp, double hp)
{
    if (!(hp > 0.0) || !std::isfinite(hp) || !std::isfinite(xp))
        return 0.0;

    if (!std::isfinite(a) || !std::isfinite(b))
        return 0.0;

    const double u0 = (a - xp) / hp;
    const double u1 = (b - xp) / hp;

    // Physical-length integral: dx = hp dxi
    return hp * integrateHatNormalized(u0, u1);
}

double integrate1DNormalized(double a, double b, double xp, double hp)
{
    if (!(hp > 0.0) || !std::isfinite(hp) || !std::isfinite(xp))
        return 0.0;

    if (!std::isfinite(a) || !std::isfinite(b))
        return 0.0;

    const double u0 = (a - xp) / hp;
    const double u1 = (b - xp) / hp;
    return integrateHatNormalized(u0, u1);
}

double rectWeight2D(double xL, double xR,
                    double yL, double yR,
                    double xp, double yp,
                    double hp)
{
    if (!(hp > 0.0) || !std::isfinite(hp))
        return 0.0;

    const double Ix = integrate1D(xL, xR, xp, hp);
    const double Iy = integrate1D(yL, yR, yp, hp);

    const double denom = hp * hp;
    const double W = (Ix * Iy) / denom;

    // Numerical guard: enforce non-negativity.
    return (W < 0.0) ? 0.0 : W;
}

bool partitionOfUnitySelfCheck(double xp, double yp, double hp,
                              double dx, double dy,
                              double tol,
                              double* outSum)
{
    if (!(hp > 0.0) || !(dx > 0.0) || !(dy > 0.0))
    {
        if (outSum) *outSum = 0.0;
        return false;
    }

    // Build a tiled rectangle that fully covers the compact support.
    // Support box: [xp-hp, xp+hp] x [yp-hp, yp+hp]. We pad by one cell on each side.
    const double xMin = (xp - hp) - dx;
    const double xMax = (xp + hp) + dx;
    const double yMin = (yp - hp) - dy;
    const double yMax = (yp + hp) + dy;

    const std::size_t nx = static_cast<std::size_t>(std::ceil((xMax - xMin) / dx)) + 1;
    const std::size_t ny = static_cast<std::size_t>(std::ceil((yMax - yMin) / dy)) + 1;

    double sumW = 0.0;
    for (std::size_t j = 0; j < ny; ++j)
    {
        const double yL = yMin + static_cast<double>(j) * dy;
        const double yR = yL + dy;
        for (std::size_t i = 0; i < nx; ++i)
        {
            const double xL = xMin + static_cast<double>(i) * dx;
            const double xR = xL + dx;
            sumW += rectWeight2D(xL, xR, yL, yR, xp, yp, hp);
        }
    }

    if (outSum)
        *outSum = sumW;

    return std::abs(sumW - 1.0) <= tol;
}

} // namespace CICAnalytic
