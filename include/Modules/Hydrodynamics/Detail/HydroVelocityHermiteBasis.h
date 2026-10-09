#pragma once

// Original fitted-velocity cubic Hermite basis, shared by the same-level and
// piecewise-side reconstructions. Expressions and double arithmetic retained.
namespace HydroFittedDetail
{
struct Basis
{
    double h[2], g[2], dh[2], dg[2], ddh[2], ddg[2];
    explicit Basis(double t)
    {
        const double t2 = t * t, t3 = t2 * t;
        h[0] = 2 * t3 - 3 * t2 + 1;
        h[1] = -2 * t3 + 3 * t2;
        g[0] = t3 - 2 * t2 + t;
        g[1] = t3 - t2;
        dh[0] = 6 * t2 - 6 * t;
        dh[1] = -dh[0];
        dg[0] = 3 * t2 - 4 * t + 1;
        dg[1] = 3 * t2 - 2 * t;
        ddh[0] = 12 * t - 6;
        ddh[1] = -ddh[0];
        ddg[0] = 6 * t - 4;
        ddg[1] = 6 * t - 2;
    }
};
} // namespace HydroFittedDetail
