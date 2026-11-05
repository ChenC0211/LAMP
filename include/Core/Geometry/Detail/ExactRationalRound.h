#pragma once
#include "Core/Geometry/Detail/ExactFaceTimeRational.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

// Implementation-only, TU-local conversion certificate. This preserves the
// existing exact algebra and rounds its retained quotient once, to the actual
// binary long-double format in round-to-nearest, ties-to-even mode. Unsupported
// exponent windows, a correctly rounded nonzero becoming zero, and overflow
// are explicit failures; no phase floor or clipping is part of conversion.
namespace LAMP::Geometry::Detail
{
namespace
{
namespace ExactRationalRound
{
namespace Exact = ExactExpansion;
using E = Exact::Expansion;
using R = FaceTime::FaceTimeRational;

[[noreturn]] inline void fail(const char* reason)
{ throw std::runtime_error(std::string("[ExactRationalRound] ") + reason); }
inline void require(bool condition, const char* reason)
{ if (!condition) fail(reason); }

// Retain every bit of a native long double as binary64 components. Native
// exponents outside the binary64 component window are deliberately unsupported.
inline E scalar(long double input)
{
    require(std::isfinite(input), "non-finite native scalar");
    E result;
    long double remainder = input;
    constexpr int pieces = (std::numeric_limits<long double>::digits + 52) / 53 + 1;
    for (int i = 0; i < pieces && remainder != 0; ++i)
    {
        const double piece = static_cast<double>(remainder);
        require(std::isfinite(piece) && piece != 0,
            "native scalar bits exceed binary64 component exponent window");
        result = Exact::add(result, Exact::fromDouble(piece));
        const long double next = remainder - static_cast<long double>(piece);
        require(next != remainder, "native scalar decomposition did not advance");
        remainder = next;
    }
    require(remainder == 0, "native scalar decomposition did not retain every bit");
    return result;
}

// Compare two positive retained ratios. Both exact cross products receive ONE
// common power-of-two scale. Its allocation to the two complete homogeneous
// ratios also respects every input bit, so a tiny candidate is never multiplied
// directly into an unrelated wide denominator. The two leading guard bits
// allow the exact products and their subtraction to remain finite.
inline int compare(const R& a, const R& b)
{
    Exact::ExponentRange an, ad, bn, bd, wholeA, wholeB;
    an.include(a.n); ad.include(a.d); bn.include(b.n); bd.include(b.d);
    wholeA.include(a.n); wholeA.include(a.d);
    wholeB.include(b.n); wholeB.include(b.d);
    const int low = std::min(an.lowest + bd.lowest, bn.lowest + ad.lowest);
    const int high = std::max(an.highest + bd.highest, bn.highest + ad.highest);
    const int aLow = -1074 - wholeA.lowest, aHigh = 1023 - wholeA.highest;
    const int bLow = -1074 - wholeB.lowest, bHigh = 1023 - wholeB.highest;
    const int commonLow = std::max(-1074 - low, aLow + bLow);
    const int commonHigh = std::min(1021 - high, aHigh + bHigh);
    require(commonLow <= commonHigh,
        "exact midpoint certificate exceeds binary64 product exponent window");
    const int common = std::clamp(-(low + high) / 2, commonLow, commonHigh);
    const int leftLow = std::max(aLow, common - bHigh);
    const int leftHigh = std::min(aHigh, common - bLow);
    require(leftLow <= leftHigh, "exact midpoint certificate has no homogeneous scale");
    const int left = std::clamp(common / 2, leftLow, leftHigh);
    const int right = common - left;
    const E lhs = Exact::multiply(Exact::scalePow2(a.n, left),
        Exact::scalePow2(b.d, right));
    const E rhs = Exact::multiply(Exact::scalePow2(b.n, right),
        Exact::scalePow2(a.d, left));
    return Exact::sign(Exact::subtract(lhs, rhs));
}

inline R midpoint(long double left, long double right)
{
    E a = scalar(left), b = scalar(right), d = Exact::fromDouble(2);
    Exact::ExponentRange range; range.include(a); range.include(b); range.include(d);
    const int shift = range.centeredShift();
    // Scale BEFORE adding: max+previous must not overflow, and min-subnormal/2
    // must not disappear before the exact homogeneous ratio is constructed.
    a = Exact::scalePow2(a, shift); b = Exact::scalePow2(b, shift);
    d = Exact::scalePow2(d, shift);
    return R(Exact::add(a, b), std::move(d));
}

inline bool evenSignificand(long double value)
{
    if (value == 0) return true;
    int exponent = 0;
    const long double fraction = std::frexp(value, &exponent);
    int precision = std::numeric_limits<long double>::digits;
    if (exponent < std::numeric_limits<long double>::min_exponent)
        precision -= std::numeric_limits<long double>::min_exponent - exponent;
    require(precision > 0, "native subnormal has no retained significand");
    const long double integer = std::scalbn(fraction, precision);
    require(std::isfinite(integer) && integer == std::trunc(integer),
        "native significand parity could not be certified");
    return std::fmod(integer, static_cast<long double>(2)) == 0;
}

inline R overflowMidpoint()
{
    constexpr int exponent = std::numeric_limits<long double>::max_exponent
        - std::numeric_limits<long double>::digits - 2;
    const long double quarterUlp = std::scalbn(static_cast<long double>(1), exponent);
    require(std::isfinite(quarterUlp) && quarterUlp > 0,
        "native overflow midpoint exceeds component exponent window");
    // (max/2 + ulp(max)/4)/(1/2) = max + ulp(max)/2. The virtual
    // upper neighbor is 2^max_exponent; it is never materialized as a scalar.
    return R(Exact::add(scalar(std::numeric_limits<long double>::max() / 2),
        scalar(quarterUlp)), Exact::fromDouble(0.5));
}

inline long double approximation(const E& value)
{
    long double result = 0;
    for (double component : value.components()) result += component;
    require(std::isfinite(result) && result > 0,
        "positive retained operand has no finite native approximation");
    return result;
}

inline long double quotient(const E& numerator, const E& denominator)
{
    static_assert(std::numeric_limits<long double>::radix == 2,
        "exact quotient rounding requires a binary native format");
    static_assert(std::numeric_limits<long double>::has_denorm == std::denorm_present,
        "exact quotient rounding requires native gradual underflow");
    Exact::validateActualEnvironment();
    const int denominatorSign = Exact::sign(denominator);
    require(denominatorSign != 0, "zero retained denominator");
    const int numeratorSign = Exact::sign(numerator);
    if (numeratorSign == 0) return 0;
    const int outputSign = numeratorSign * denominatorSign;
    const R exact(numeratorSign < 0 ? Exact::subtract(E{}, numerator) : numerator,
        denominatorSign < 0 ? Exact::subtract(E{}, denominator) : denominator);
    const long double seed = approximation(exact.n) / approximation(exact.d);
    require(!std::isnan(seed) && seed >= 0, "invalid native quotient seed");
    long double candidate = seed;
    if (candidate == 0) candidate = std::numeric_limits<long double>::denorm_min();
    if (!std::isfinite(candidate)) candidate = std::numeric_limits<long double>::max();
    constexpr unsigned maximumNeighborSteps = 16;
    for (unsigned step = 0; step < maximumNeighborSteps; ++step)
    {
        const bool even = evenSignificand(candidate);
        const long double lower = std::nextafter(candidate, static_cast<long double>(0));
        const int lowerComparison = compare(exact, midpoint(lower, candidate));
        if (lowerComparison < 0 || (lowerComparison == 0 && !even))
        {
            require(lower != 0, "correctly rounded nonzero quotient would become zero");
            candidate = lower;
            continue;
        }
        const long double upper = std::nextafter(candidate,
            std::numeric_limits<long double>::infinity());
        const int upperComparison = compare(exact,
            std::isfinite(upper) ? midpoint(candidate, upper) : overflowMidpoint());
        if (upperComparison > 0 || (upperComparison == 0 && !even))
        {
            require(std::isfinite(upper), "correctly rounded quotient would overflow");
            candidate = upper;
            continue;
        }
        Exact::validateActualEnvironment();
        return outputSign < 0 ? -candidate : candidate;
    }
    fail("native quotient seed exceeded finite neighbor certificate budget");
}
}
}
}