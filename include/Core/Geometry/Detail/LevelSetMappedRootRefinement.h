#pragma once

#include "Core/Geometry/Detail/LevelSetMappedCoordinates.h"
#include "Core/Geometry/Detail/PolynomialSegmentRoots.h"
#include <algorithm>
#include <array>
#include <cstddef>
#include <stdexcept>
#include <string>

namespace Quantum::Geometry::Detail
{
namespace
{
namespace MappedHermiteRootRefinement
{
namespace M = MappedHermiteSegment;
namespace Exact = ExactExpansion;
using Ratio = M::Ratio;
using Polynomial = PolynomialSegmentRoots::ExpandedPolynomial;

inline void require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(std::string("[MappedHermiteRootRefinement] ")+message);
}

struct RefinedRoot
{
    // New retained child chart, possibly directed from upper to lower.
    // The three doubles belong ONLY to that chart, never the original
    // segment or a rounded global coordinate.
    M::Chart chart;
    double parameter=0.0, lower=0.0, upper=0.0;
    bool reversed=false;

    // Same original forward parameter space as the input bracket. The
    // enclosure is non-widening and strictly narrower on at least one side.
    Ratio originalParameter, originalLower, originalUpper;

    // Nonzero sign of the scalar derivative with respect to the original
    // increasing parameter. Together with the single-root inventory and
    // nested bracket this retains the parent's simple-root identity.
    int derivativeSign=0;
};

inline Exact::Expansion absolute(Exact::Expansion value)
{
    return Exact::sign(value)<0?Exact::subtract(Exact::Expansion{},value):value;
}

template<std::size_t Count>
inline int endpointDerivativeSign(const std::array<Exact::Expansion,Count>& polynomial, bool upper)
{
    if (!upper) return Exact::sign(polynomial[1]);
    Exact::Expansion derivative;
    for (std::size_t k=1;k<Count;++k)
        derivative=Exact::add(derivative,Exact::multiply(polynomial[k],M::number(static_cast<double>(k))));
    return Exact::sign(derivative);
}

// Relocate and refine ONE previously isolated simple root. The caller owns
// the parent root identity and its complete, ordered, positive-width bracket
// [lower,upper]. restriction(a,b) must rebuild the SAME retained cell scalar
// on a+s*(b-a), with a strictly positive clearing factor. It must not refit
// values, snap endpoints, drop coefficients, or map a/b through double.
//
// One unchanged endpoint-relative isolation is performed in the child chart.
// Its complete inventory must contain exactly one root. An exact endpoint
// root is allowed only when its exact derivative is nonzero. All unresolved
// multiplicity, zero restrictions, widening, missing progress and arithmetic
// scale failures reject explicitly. This helper supplies no iteration loop,
// alternate threshold, fallback, global-coordinate approximation or policy
// for deciding whether a consumer needs another refinement.
template<class Restriction>
inline RefinedRoot refine(const Ratio& lower, const Ratio& upper, Restriction&& restriction)
{
    Exact::validateEnvironment();
    require(Exact::sign(lower.denominator)>0 && Exact::sign(upper.denominator)>0,
        "input bracket denominators must be positive");
    require(M::compare(lower,upper)<0, "refinement requires an ordered positive-width bracket");

    RefinedRoot result;
    result.chart={lower,upper};
    auto polynomial=restriction(result.chart.begin,result.chart.end);
    if (Exact::sign(Exact::subtract(absolute(Exact::evaluate(polynomial,1.0)),absolute(polynomial[0])))<0)
    {
        result.reversed=true;
        std::swap(result.chart.begin,result.chart.end);
        polynomial=restriction(result.chart.begin,result.chart.end);
    }

    const auto isolated=PolynomialSegmentRoots::isolateEndpointRelative(polynomial);
    require(!isolated.identicallyZero, "parent bracket restricted to an identically zero polynomial");
    require(isolated.roots.size()==1 && isolated.enclosures.size()==1,
        "refined bracket does not contain exactly one isolated root");
    result.parameter=isolated.roots.front();
    result.lower=isolated.enclosures.front().lo;
    result.upper=isolated.enclosures.front().hi;
    require(result.lower>=0.0 && result.upper<=1.0 && result.lower<=result.parameter
        && result.parameter<=result.upper, "invalid child root representative or enclosure");

    int childDerivativeSign=0;
    if (result.parameter==0.0 || result.parameter==1.0)
    {
        const bool atUpper=result.parameter==1.0;
        const int endpointSign=Exact::sign(atUpper?Exact::evaluate(polynomial,1.0):polynomial[0]);
        require(endpointSign==0, "child root collapsed onto a nonzero endpoint");
        childDerivativeSign=endpointDerivativeSign(polynomial,atUpper);
        require(childDerivativeSign!=0, "refined endpoint root is not simple");
    }
    else
    {
        // Complete isolation certifies this interior root simple. There are
        // no endpoint roots in a one-root inventory, so the two nonzero
        // endpoint phases must be opposite. Their crossing fixes P'(root)
        // without a rounded-coordinate derivative evaluation.
        const int firstSign=Exact::sign(polynomial[0]);
        const int lastSign=Exact::sign(Exact::evaluate(polynomial,1.0));
        require(firstSign!=0 && lastSign==-firstSign,
            "single interior-root inventory disagrees with endpoint phases");
        childDerivativeSign=-firstSign;
    }
    result.derivativeSign=result.reversed?-childDerivativeSign:childDerivativeSign;

    result.originalParameter=M::between(result.chart.begin,result.chart.end,result.parameter);
    result.originalLower=M::between(result.chart.begin,result.chart.end,result.lower);
    result.originalUpper=M::between(result.chart.begin,result.chart.end,result.upper);
    if (M::compare(result.originalLower,result.originalUpper)>0)
        std::swap(result.originalLower,result.originalUpper);
    require(M::compare(result.originalLower,lower)>=0 && M::compare(result.originalUpper,upper)<=0,
        "refined enclosure widened beyond its parent bracket");
    require(M::compare(result.originalLower,result.originalParameter)<=0
        && M::compare(result.originalParameter,result.originalUpper)<=0,
        "refined representative left its original-parameter enclosure");
    require(M::compare(result.originalLower,lower)>0 || M::compare(result.originalUpper,upper)<0,
        "child isolation did not strictly narrow the parent bracket");
    return result;
}
}
}
}
