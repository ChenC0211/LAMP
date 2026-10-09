#pragma once

#include <cmath>
#include <stdexcept>

namespace DielectricPolarizationForce
{
// Charge-free (q=0) polarization coefficient. dn is the common
// oriented face displacement and tangent is the reconstructed tangential E.
// Multiplication by the receiving gradient gives the Maxwell Tnn jump per
// distance for constant Dn and Et in a planar layer.
// General fields use a discrete polarization approximation, not a sharp
// two-sided interface trace or a discrete energy-variation certificate.
// Do not add a second Maxwell interface traction to the same electric load.
// All inputs must use mutually consistent caller-scaled units; the caller
// applies its scalar gradient, mechanical normalization and overflow checks.
inline double gradientCoefficient(double dn, double tangent, double epsilonLow, double epsilonHigh)
{
    if (!std::isfinite(epsilonLow) || !(epsilonLow > 0.0)
        || !std::isfinite(epsilonHigh) || !(epsilonHigh > 0.0)
        || !std::isfinite(dn) || !std::isfinite(tangent))
        throw std::invalid_argument("dielectric polarization force requires positive finite permittivities and finite inputs");
    const double normal2 = (dn / epsilonLow) * (dn / epsilonHigh);
    return -0.5 * (tangent * tangent + normal2);
}
}
