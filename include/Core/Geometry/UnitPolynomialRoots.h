#pragma once
#include <array>
#include <vector>

namespace LAMP::Geometry
{
struct UnitPolynomialRoots
{
    struct Root { double point=0,lower=0,upper=0; };
    std::vector<Root> roots;
    bool identicallyZero=false;
};
struct UnitPolynomialTerm
{
    std::array<double,12> coefficients{}; // Ascending powers, including native C2 face moments through degree eleven.
    long double factor=0;
    double weight=1; // Kept separate for exact weighted combinations of traces.
};
// Sum of weight*factor*polynomial terms, degree <=11, on [0,1]. Reuses the existing
// bounded, outward-enclosed polynomial engine. No small-coefficient deletion
// or unresolved multiple-root merging. Endpoint roots retain their identity.
// Form the coefficient sums/products in expansion arithmetic so algebraic
// endpoint cancellation is not lost before isolation. Factors must have a
// finite exact binary64 expansion; wider long doubles are decomposed.
UnitPolynomialRoots isolateUnitPolynomial(const std::vector<UnitPolynomialTerm>& terms);
// Evaluate the same retained combination in expansion arithmetic. In
// particular, classification near an isolated root must not lose cancellation
// by separately rounding each weighted term before summation.
long double evaluateUnitPolynomial(const std::vector<UnitPolynomialTerm>& terms,double parameter);
// Outward Bernstein range of the SAME exact weighted polynomial on [a,b].
// Combine cancelling terms before bounding; no point sampling or coefficient
// deletion. Returned bounds are in polynomial-density units.
std::array<double,2> encloseUnitPolynomial(const std::vector<UnitPolynomialTerm>& terms,double a=0,double b=1);
// Integral on t=a0+(a1-a0)*s, s in [s0,s1]. Keep the anchored endpoints,
// coefficient combinations and antiderivative cancellation in expansions;
// only the final division by a common integer denominator is rounded.
long double integrateUnitPolynomial(const std::vector<UnitPolynomialTerm>& terms,
    double a0,double a1,double s0,double s1);
// Degree-exact positive Gauss integration of q(t)*B_k^degree((half+t)/2).
// The degree <=11 normal trace is evaluated in expansion arithmetic; its
// factored nonnegative work weight has degree <=20. No high-degree root
// isolation, ill-conditioned power expansion or coefficient fitting.
struct UnitPolynomialIntegral;
UnitPolynomialIntegral integrateUnitPolynomialBernstein(const std::vector<UnitPolynomialTerm>& terms,
    int degree,int k,int half,double lo,double hi);
// All k=0,...,degree on the same interval and the same positive rule.
// Share the exact trace evaluation; retain each weight's value/error sum.
std::vector<UnitPolynomialIntegral> integrateUnitPolynomialBernsteinBasis(const std::vector<UnitPolynomialTerm>& terms,
    int degree,int half,double lo,double hi);
// Same primitive on a terminal mapped-segment range. Keep offset and length
// separate even when their sum rounds to offset; fromEnd means s=1-offset-q.
long double integrateUnitPolynomialRange(const std::vector<UnitPolynomialTerm>& terms,
    double a0,double a1,bool fromEnd,long double offset,long double length);
struct UnitPolynomialIntegral {long double value=0,absoluteError=0;};
// Integral of max(p,0), or its exact-in-time linear interpolation average.
// Bernstein sign bounds classify whole intervals. Unresolved sign changes
// are bounded in INTEGRAL units; no representable global root is required.
// On strictly opposite endpoint signs, smooth reversal quadrature uses the
// same monitored Gauss rule as the face-time primitive. The returned error
// includes unresolved-sign bounds and the accumulated quadrature estimate.
UnitPolynomialIntegral integratePositiveUnitPolynomial(const std::vector<UnitPolynomialTerm>&,
    double lo,double hi,long double absoluteTolerance);
// SAME positive primitive on t=a0+(a1-a0)*s, with s=offset+q or
// 1-offset-q, 0<=q<=length. Exact polynomial translation/scaling preserves
// ranges whose two floating endpoints would coincide. Returns integral in
// original t units; the tolerance has those same integral units.
UnitPolynomialIntegral integratePositiveUnitPolynomialRange(const std::vector<UnitPolynomialTerm>&,
    double a0,double a1,bool fromEnd,long double offset,long double length,long double absoluteTolerance);
UnitPolynomialIntegral integratePositiveTimeLinearUnitPolynomial(const std::vector<UnitPolynomialTerm>& begin,
    const std::vector<UnitPolynomialTerm>& end,double lo,double hi,long double absoluteTolerance);
}
