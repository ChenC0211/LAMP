#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>

// ----------------------------------------------------------------------------
// CIC analytic integration utilities (explicit electrostatic PIC)
//
// 1D hat / CIC shape function:
//   s(\xi) = max(0, 1 - |\xi|)
//
// 2D kernel used by this module:
//   S(x,y) = s((x-xp)/hp) * s((y-yp)/hp) / hp^2
//
// For a rectangular cell K = [xL,xR] x [yL,yR], the integrated weight is
//
//   W_{pK} = \int_K S(x,y) dA
//         = ( \int_{xL}^{xR} s((x-xp)/hp) dx )
//           ( \int_{yL}^{yR} s((y-yp)/hp) dy ) / hp^2
//         = Ix * Iy / hp^2.
//
// Notes:
// - Ix, Iy here are *physical-length* integrals of the unnormalized hat.
// - W_{pK} is dimensionless and lies in [0,1].
// - Outside the support (|x-xp|>=hp) the 1D integrand is zero.
// ----------------------------------------------------------------------------

namespace CICAnalytic {

// 1D hat / CIC shape function.
inline double hat1D(const double xi)
{
    const double a = std::abs(xi);
    return (a >= 1.0) ? 0.0 : (1.0 - a);
}

// Analytic integral of hat1D over a *normalized* interval [u0,u1] in xi-space:
//   I = \int_{u0}^{u1} s(xi) dxi, where s(xi)=max(0,1-|xi|).
//
// The interval endpoints may be in any order; the function returns a
// non-negative integral.
//
// This normalized integral is dimensionless and bounded by 1.
// In particular, \int_{-1}^{1} s(xi) dxi = 1.

double integrateHatNormalized(double u0, double u1);

// 1D analytic interval integral in physical space (x-units):
//   Ix = \int_a^b s((x-xp)/hp) dx.
//
// Returns 0 if hp<=0 or the interval does not overlap the support.
// The integral is always non-negative.

double integrate1D(double a, double b, double xp, double hp);

// Normalized 1D integral:
//   Ix_norm = (1/hp) * \int_a^b s((x-xp)/hp) dx
//           = \int_{(a-xp)/hp}^{(b-xp)/hp} s(xi) dxi.
//
// This is dimensionless and is equal to integrate1D(a,b,xp,hp)/hp.

double integrate1DNormalized(double a, double b, double xp, double hp);

// 2D rectangular weight for this kernel:
//   W = Ix * Iy / hp^2.

double rectWeight2D(double xL, double xR,
                    double yL, double yR,
                    double xp, double yp,
                    double hp);

// Small self-check helper: tile a rectangle covering the full CIC support
// with uniform rectangles of size (dx,dy), and verify partition-of-unity:
//   sum_{cells} W_{pK} \approx 1 (uncut case).
//
// This function is framework-free and can be used by tiny unit tests.
// Returns true if |sum-1| <= tol. If outSum is non-null, writes the sum.

bool partitionOfUnitySelfCheck(double xp, double yp, double hp,
                              double dx, double dy,
                              double tol,
                              double* outSum);

} // namespace CICAnalytic
