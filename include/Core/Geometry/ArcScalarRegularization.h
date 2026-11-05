#pragma once
#include <vector>

namespace LAMP::Geometry::ArcScalarRegularization
{
// Closed ordered arc, one scalar at each knot. edgeLengths[i] is the positive
// arc length from i to (i+1)%n, not a Euclidean proximity connection. The caller
// owns branch/endpoint/support validity and supplies a common physical scale h.
// Minimize ||sqrt(M)(u-f)||^2 + h^6 ||sqrt(E) B3 u||^2,
// with dual arc mass M, central-edge lengths E and B3 = six times the
// four-knot third divided difference. B3 annihilates quadratics even on
// irregular knots. At least four closed knots are required. Sparse augmented QR
// avoids forming normal equations. No input scalar or geometry is changed.
// On uniform spacing ds the Fourier gain is
// 1/(1+64*(h/ds)^6*sin(theta/2)^6), positive at every wave number.
// This is a scalar reconstruction operator, not a coupled stability guarantee.
struct Result {
    std::vector<double> values;
    double weightedMeanChange=0;
    double relativeNormalResidual=0;
    double dataEnergy=0,derivativeEnergy=0;
};
Result closedThirdDerivative(const std::vector<double>& edgeLengths,
    const std::vector<double>& values,double h);
enum class Endpoint { FixedValue, FixedEven };
// Finite interval with the original scalar values fixed at BOTH endpoints.
// A FixedValue end uses only physical knots: no reflection, imposed slope,
// or wrap to the other end. The endpoint third-difference stencil integrates
// to the physical boundary. This is the discrete variational boundary of
// the roughness functional, not a contact-angle model. FixedEven additionally
// reflects/ties knots across that end; useful for a regular axis limit.
// At least four physical knots, with n-1 positive lengths, are required.
Result openThirdDerivative(const std::vector<double>& edgeLengths,
    const std::vector<double>& values,double h,Endpoint begin,Endpoint end);
// Open interval with EVEN reflection at BOTH endpoints. Endpoint values are
// hard constraints from the caller's existing endpoint estimator. Reflection
// pairs share one unknown, so symmetry is imposed algebraically, not by
// averaging an unconstrained solve. n values require n-1 positive edge lengths.
// This is not a boundary model for a general wall, nozzle or contact line.
Result evenOpenThirdDerivative(const std::vector<double>& edgeLengths,
    const std::vector<double>& values,double h);
}
