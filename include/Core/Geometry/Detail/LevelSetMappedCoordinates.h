#pragma once

#include "Core/Geometry/Detail/ExactExpansion.h"

namespace Quantum::Geometry::Detail
{
namespace
{
namespace MappedHermiteSegment
{
namespace Exact = ExactExpansion;
using Expansion = Exact::Expansion;

inline void require(bool valid, const char* reason)
{
    if (!valid) throw std::runtime_error(std::string("[LevelSetMappedSegment] ") + reason);
}
inline Expansion number(double value) { return Exact::fromDouble(value); }
inline Expansion difference(double a, double b) { return Exact::subtract(number(a), number(b)); }

// Retained exact-rational coordinates, shared by face partitions and curve
// quadrature. Only terminal integration converts a ratio to floating point.
// The expansion engine's finite exponent range remains checked explicitly.
struct Ratio
{
    Expansion numerator, denominator = number(1.0);
};
inline Ratio ratio(Expansion numerator, Expansion denominator)
{
    const int sign = Exact::sign(denominator);
    require(sign != 0, "zero chart denominator");
    if (sign < 0)
    {
        numerator = Exact::subtract(Expansion{}, numerator);
        denominator = Exact::subtract(Expansion{}, denominator);
    }
    if (Exact::sign(numerator) == 0) return {};
    if (Exact::sign(Exact::subtract(numerator, denominator)) == 0) return {number(1.0), number(1.0)};
    return {std::move(numerator), std::move(denominator)};
}
inline Ratio constant(double value) { return {number(value), number(1.0)}; }
inline bool sameDenominator(const Ratio& a, const Ratio& b)
{
    return Exact::sign(Exact::subtract(a.denominator,b.denominator)) == 0;
}
inline int compare(const Ratio& a, const Ratio& b)
{
    if (sameDenominator(a,b)) return Exact::sign(Exact::subtract(a.numerator,b.numerator));
    return Exact::sign(Exact::subtract(Exact::multiply(a.numerator, b.denominator),
        Exact::multiply(b.numerator, a.denominator)));
}
inline Ratio add(const Ratio& a, const Ratio& b)
{
    if (sameDenominator(a,b)) return ratio(Exact::add(a.numerator,b.numerator),a.denominator);
    return ratio(Exact::add(Exact::multiply(a.numerator, b.denominator),
        Exact::multiply(b.numerator, a.denominator)), Exact::multiply(a.denominator, b.denominator));
}
inline Ratio subtract(const Ratio& a, const Ratio& b)
{
    if (sameDenominator(a,b)) return ratio(Exact::subtract(a.numerator,b.numerator),a.denominator);
    return ratio(Exact::subtract(Exact::multiply(a.numerator, b.denominator),
        Exact::multiply(b.numerator, a.denominator)), Exact::multiply(a.denominator, b.denominator));
}
inline Ratio between(const Ratio& a, const Ratio& b, double parameter)
{
    if (sameDenominator(a,b))
        return ratio(Exact::add(a.numerator,Exact::multiply(Exact::subtract(b.numerator,a.numerator),number(parameter))),
            a.denominator);
    const auto origin = Exact::multiply(a.numerator, b.denominator);
    const auto delta = Exact::subtract(Exact::multiply(b.numerator, a.denominator), origin);
    return ratio(Exact::add(origin, Exact::multiply(delta, number(parameter))),
        Exact::multiply(a.denominator, b.denominator));
}
inline long double approximate(const Ratio& value)
{
    const long double result = static_cast<long double>(Exact::estimate(value.numerator))
        / Exact::estimate(value.denominator);
    require(std::isfinite(result) && (Exact::sign(value.numerator) == 0 || result != 0.0L),
        "unrepresentable integration ratio");
    return result;
}
struct Chart { Ratio begin, end; };
}
}
}
