#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "Core/Geometry/CellChartPoint.h"
#include "Core/Base/Vector2D.h"
#include "Core/Geometry/Detail/ExactExpansion.h"

// Shared single-cell Hermite arithmetic. No mesh, node reconstruction, root
// selection, boundary law or IO. Views borrow retained scaled nodal doubles and
// their original normalized cache; the caller owns that one scalar authority.
// Keep TU-local linkage and the existing strict floating-point boundary.
namespace LAMP::Geometry::Detail
{
namespace
{
namespace LevelSetHermiteScalar
{
struct Sample
{
    double value = 0.0;
    double dx = 0.0, dy = 0.0;
    double dxx = 0.0, dxy = 0.0, dyy = 0.0;
};

// Physical third derivatives of ONE selected bicubic cell polynomial. These
// need not agree across cell boundaries: the shared Hermite scalar is C1,
// not C2 or C3, and a declared sector seam is only C0. This separate result/API
// leaves the legacy Sample path intact.
struct ThirdDerivatives
{
    double dxxx = 0.0, dxxy = 0.0, dxyy = 0.0, dyyy = 0.0;
};

namespace Detail
{
inline void require(bool condition, const char *reason)
{
    if (!condition)
        throw std::runtime_error(std::string("[LevelSetHermiteScalarCandidate] ") + reason);
}

inline double checked(long double value, const char *reason)
{
    const double result = static_cast<double>(value);
    require(std::isfinite(value) && std::isfinite(result), reason);
    return result;
}

// Ordering: endpoint values at 0/1, then endpoint derivatives at 0/1.
// Derivatives here are with respect to the dimensionless cell coordinate.
struct Basis
{
    std::array<long double, 4> value, first, second;
    explicit Basis(double coordinate)
    {
        const long double x = coordinate, square = x * x;
        value = {{(2.0L * x - 3.0L) * square + 1.0L, (-2.0L * x + 3.0L) * square, ((x - 2.0L) * x + 1.0L) * x,
                  (x - 1.0L) * square}};
        first = {{6.0L * square - 6.0L * x, -6.0L * square + 6.0L * x, 3.0L * square - 4.0L * x + 1.0L,
                  3.0L * square - 2.0L * x}};
        second = {{12.0L * x - 6.0L, 6.0L - 12.0L * x, 6.0L * x - 4.0L, 6.0L * x - 2.0L}};
    }
    explicit Basis(CellChartCoordinate coordinate) : Basis(coordinate.offset)
    {
        if (!coordinate.fromUpper)
            return;
        // h00(1-a)=h01(a), h10(1-a)=-h11(a), and conversely.
        // Odd derivatives also change sign under a=1-u. No 1-a is formed.
        const auto lowerValue = value, lowerFirst = first, lowerSecond = second;
        constexpr std::array<int, 4> reflected{{1, 0, 3, 2}};
        constexpr std::array<int, 4> parity{{1, 1, -1, -1}};
        for (int k = 0; k < 4; ++k)
        {
            value[k] = parity[k] * lowerValue[reflected[k]];
            first[k] = -parity[k] * lowerFirst[reflected[k]];
            second[k] = parity[k] * lowerSecond[reflected[k]];
        }
    }
};

// Allocation-free outward interval arithmetic for filters only. Public filter
// entry points check the floating-point environment ONCE. Each rounded
// operation is enclosed by adjacent binary64 values; overflow/indeterminate
// intermediate ranges become the whole real line, never NaN or a false sign.
// Inputs to point() are already validated nodal/local values or small integers.
struct Outward
{
    using Interval = ExactExpansion::Interval;
    static Interval whole()
    {
        return {-std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
    }
    static Interval point(double value)
    {
        return {value, value};
    }
    static bool zero(Interval value)
    {
        return value.lo == 0.0 && value.hi == 0.0;
    }
    static bool finite(Interval value)
    {
        return std::isfinite(value.lo) && std::isfinite(value.hi) && value.lo <= value.hi;
    }
    static Interval negate(Interval value)
    {
        return {-value.hi, -value.lo};
    }
    static double down(double value)
    {
        return std::nextafter(value, -std::numeric_limits<double>::infinity());
    }
    static double up(double value)
    {
        return std::nextafter(value, std::numeric_limits<double>::infinity());
    }
    static Interval add(Interval a, Interval b)
    {
        if (!finite(a) || !finite(b))
            return whole();
        if (zero(a))
            return b;
        if (zero(b))
            return a;
        const double lo = a.lo + b.lo, hi = a.hi + b.hi;
        if (!std::isfinite(lo) || !std::isfinite(hi))
            return whole();
        return {down(lo), up(hi)};
    }
    static Interval multiply(Interval a, Interval b)
    {
        if (!finite(a) || !finite(b))
            return whole();
        if (zero(a) || zero(b))
            return {};
        if (a.lo == 1.0 && a.hi == 1.0)
            return b;
        if (b.lo == 1.0 && b.hi == 1.0)
            return a;
        if (a.lo == -1.0 && a.hi == -1.0)
            return negate(b);
        if (b.lo == -1.0 && b.hi == -1.0)
            return negate(a);
        // Monotonicity selects the SAME extrema of the four endpoint
        // products. In particular, nonnegative Bernstein weights need only
        // two products. Keep the original adjacent-binary64 enclosure and
        // overflow rule; only two intervals straddling zero need all four.
        double lo, hi;
        if (a.lo >= 0.0)
        {
            if (b.lo >= 0.0) { lo = a.lo * b.lo; hi = a.hi * b.hi; }
            else if (b.hi <= 0.0) { lo = a.hi * b.lo; hi = a.lo * b.hi; }
            else { lo = a.hi * b.lo; hi = a.hi * b.hi; }
        }
        else if (a.hi <= 0.0)
        {
            if (b.lo >= 0.0) { lo = a.lo * b.hi; hi = a.hi * b.lo; }
            else if (b.hi <= 0.0) { lo = a.hi * b.hi; hi = a.lo * b.lo; }
            else { lo = a.lo * b.hi; hi = a.lo * b.lo; }
        }
        else
        {
            if (b.lo >= 0.0) { lo = a.lo * b.hi; hi = a.hi * b.hi; }
            else if (b.hi <= 0.0) { lo = a.hi * b.lo; hi = a.lo * b.lo; }
            else
            {
                lo = std::min(a.lo * b.hi, a.hi * b.lo);
                hi = std::max(a.lo * b.lo, a.hi * b.hi);
            }
        }
        if (!std::isfinite(lo) || !std::isfinite(hi))
            return whole();
        return {down(lo), up(hi)};
    }
    // A positive represented divisor; Bernstein controls use integers 3 and 9.
    static Interval dividePositive(Interval value, double divisor)
    {
        if (!finite(value) || !std::isfinite(divisor) || !(divisor > 0.0))
            return whole();
        if (zero(value) || divisor == 1.0)
            return value;
        const double lo = value.lo / divisor, hi = value.hi / divisor;
        if (!std::isfinite(lo) || !std::isfinite(hi))
            return whole();
        return {down(lo), up(hi)};
    }
    static std::array<Interval, 4> basis(double coordinate)
    {
        // Endpoint identities are exact properties of the same integer
        // Hermite basis, not a near-endpoint or shape-dependent shortcut.
        if (coordinate == 0.0)
            return {{point(1.0), Interval{}, Interval{}, Interval{}}};
        if (coordinate == 1.0)
            return {{Interval{}, point(1.0), Interval{}, Interval{}}};
        const std::array<std::array<int, 4>, 4> coefficients{
            {{{1, 0, -3, 2}}, {{0, 0, 3, -2}}, {{0, 1, -2, 1}}, {{0, 0, -1, 1}}}};
        const Interval argument = point(coordinate);
        std::array<Interval, 4> result;
        for (int row = 0; row < 4; ++row)
        {
            result[row] = point(coefficients[row][3]);
            for (int power = 2; power >= 0; --power)
                result[row] = add(multiply(result[row], argument), point(coefficients[row][power]));
        }
        return result;
    }
};
} // namespace Detail

struct Node
{
    double value = 0.0, x = 0.0, y = 0.0, xy = 0.0;
};
struct CellData
{
    double offset = 0.0, scale = 0.0;
    std::array<double, 16> data{};
};
using ExpandedPolynomial = std::array<ExactExpansion::Expansion, 7>;
struct RationalSegmentCoordinate
{
    ExactExpansion::Expansion numeratorBegin, numeratorEnd, denominator;
};

// Original geometry candidate test, shared by logical and native layouts.
// These are floating resolution checks, not exact zero-set certificates.
template<class Sampler>
inline std::array<long double,16> sampledBernstein(Vector2D h, Sampler sample)
{
    std::array<long double,16> b{};
    for (int y=0;y<2;++y) for (int x=0;x<2;++x)
    {
        const auto s=sample(x,y);
        const int bx=3*x,by=3*y,sx=x?-1:1,sy=y?-1:1;
        const long double dx=static_cast<long double>(s.dx)*h.x/3.0L;
        const long double dy=static_cast<long double>(s.dy)*h.y/3.0L;
        const long double dxy=static_cast<long double>(s.dxy)*h.x*h.y/9.0L;
        b[4*by+bx]=s.value;
        b[4*by+bx+sx]=s.value+sx*dx;
        b[4*(by+sy)+bx]=s.value+sy*dy;
        b[4*(by+sy)+bx+sx]=s.value+sx*dx+sy*dy+sx*sy*dxy;
    }
    for (auto v:b) Detail::require(std::isfinite(v),"non-finite cell Bernstein coefficient");
    return b;
}
// Same outward tensor blossom and convex hull for cubic/quartic controls.
// The degree changes only array extents; endpoint identities and arithmetic
// order are retained. This does not certify a sign if the bound contains zero.
namespace Detail
{
template<std::size_t Count>
inline std::array<ExactExpansion::Interval, Count> tensorBernsteinRestriction(
    std::array<ExactExpansion::Interval, Count> controls, Vector2D localLower, Vector2D localUpper)
{
    static_assert(Count == 16 || Count == 25, "supported cubic/quartic tensor controls");
    constexpr int Side = Count == 16 ? 4 : 5, Degree = Side - 1;
    using Interval = ExactExpansion::Interval;
    using Arithmetic = Outward;
    const auto blend = [](Interval a, Interval b, double t) {
        // These are exact Bernstein endpoint identities. In particular
        // they neither widen an unchanged endpoint nor consume a control
        // multiplied by an exact zero.
        if (t == 0.0)
            return a;
        if (t == 1.0)
            return b;
        const auto argument = Arithmetic::point(t);
        const auto complement = Arithmetic::add(Arithmetic::point(1.0), Arithmetic::negate(argument));
        return Arithmetic::add(Arithmetic::multiply(complement, a), Arithmetic::multiply(argument, b));
    };
    const auto restrictLine = [&](const auto &line, double lower, double upper) {
        auto restricted = line;
        for (int k = 0; k < Side; ++k)
        {
            auto work = line;
            // The kth control of p(lower+(upper-lower)*s) is the
            // symmetric degree-D blossom at (lower^(D-k), upper^k).
            // This avoids division by a possibly tiny subinterval length
            // and remains valid when lower==upper.
            for (int stage = 0; stage < Degree; ++stage)
            {
                const double argument = stage < Degree - k ? lower : upper;
                for (int entry = 0; entry < Degree - stage; ++entry)
                    work[entry] = blend(work[entry], work[entry + 1], argument);
            }
            restricted[k] = work[0];
        }
        return restricted;
    };
    if (localLower.x != 0.0 || localUpper.x != 1.0)
        for (int y = 0; y < Side; ++y)
        {
            std::array<Interval, Side> line;
            for (int x = 0; x < Side; ++x)
                line[x] = controls[Side * y + x];
            const auto restricted = restrictLine(line, localLower.x, localUpper.x);
            for (int x = 0; x < Side; ++x)
                controls[Side * y + x] = restricted[x];
        }
    if (localLower.y != 0.0 || localUpper.y != 1.0)
        for (int x = 0; x < Side; ++x)
        {
            std::array<Interval, Side> line;
            for (int y = 0; y < Side; ++y)
                line[y] = controls[Side * y + x];
            const auto restricted = restrictLine(line, localLower.y, localUpper.y);
            for (int y = 0; y < Side; ++y)
                controls[Side * y + x] = restricted[y];
        }
    return controls;
}
template<std::size_t Count>
inline ExactExpansion::Interval tensorBernsteinBounds(
    std::array<ExactExpansion::Interval, Count> controls, Vector2D localLower, Vector2D localUpper)
{
    using Interval = ExactExpansion::Interval;
    using Arithmetic = Outward;
    controls = tensorBernsteinRestriction(controls, localLower, localUpper);
    Interval result{std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};
    for (auto control : controls)
    {
        if (!Arithmetic::finite(control))
            return Arithmetic::whole();
        result.lo = std::min(result.lo, control.lo);
        result.hi = std::max(result.hi, control.hi);
    }
    return result;
}
} // namespace Detail

template<std::size_t Count>
inline bool resolvedCartesianPartial(const std::array<long double,Count>& b, long double guard,
    int selectedAxis = -1)
{
    static_assert(Count == 16 || Count == 25, "supported cubic/quartic tensor controls");
    constexpr int Side = Count == 16 ? 4 : 5, Degree = Side - 1;
    bool regular=false;
    for (int axis=0;axis<2;++axis)
    {
        if (selectedAxis >= 0 && axis != selectedAxis) continue;
        long double low=std::numeric_limits<long double>::infinity(),high=-low;
        for (int y=0;y<(axis?Degree:Side);++y) for (int x=0;x<(axis?Side:Degree);++x)
        {
            const auto d=static_cast<long double>(Degree)*(b[Side*(y+(axis?1:0))+x+(axis?0:1)]-b[Side*y+x]);
            low=std::min(low,d); high=std::max(high,d);
        }
        regular=regular || low>(2.0L*Degree)*guard || high<-(2.0L*Degree)*guard;
    }
    return regular;
}

class CellView
{
  public:
    using Interval = ExactExpansion::Interval;
    CellView(std::array<const Node *, 4> corners, const CellData &cache, double hx, double hy)
        : corners_(corners), cache_(cache), hx_(hx), hy_(hy)
    {
    }

    Sample evaluateCellBasis(const Detail::Basis &bx, const Detail::Basis &by) const
    {
        const auto &cell = cache_;
        if (cell.scale == 0.0)
        {
            Sample result;
            result.value = cell.offset;
            return result;
        }
        long double value = 0.0L, dx = 0.0L, dy = 0.0L;
        long double dxx = 0.0L, dxy = 0.0L, dyy = 0.0L;
        for (int y = 0; y < 4; ++y)
            for (int x = 0; x < 4; ++x)
            {
                const long double coefficient = cell.data[static_cast<std::size_t>(4 * y + x)];
                value += coefficient * bx.value[x] * by.value[y];
                dx += coefficient * bx.first[x] * by.value[y];
                dy += coefficient * bx.value[x] * by.first[y];
                dxx += coefficient * bx.second[x] * by.value[y];
                dxy += coefficient * bx.first[x] * by.first[y];
                dyy += coefficient * bx.value[x] * by.second[y];
            }
        Sample result;
        result.value = Detail::checked(static_cast<long double>(cell.offset) + cell.scale * value,
                                       "non-finite reconstructed scalar");
        result.dx = Detail::checked(cell.scale * dx / hx_, "non-finite reconstructed dx");
        result.dy = Detail::checked(cell.scale * dy / hy_, "non-finite reconstructed dy");
        result.dxx = Detail::checked(cell.scale * dxx / hx_ / hx_, "non-finite reconstructed dxx");
        result.dxy = Detail::checked(cell.scale * dxy / hx_ / hy_, "non-finite reconstructed dxy");
        result.dyy = Detail::checked(cell.scale * dyy / hy_ / hy_, "non-finite reconstructed dyy");
        return result;
    }

    ThirdDerivatives thirdDerivativesCellBasis(const Detail::Basis &bx, const Detail::Basis &by) const
    {
        const auto &cell = cache_;
        if (cell.scale == 0.0)
            return {};
        // Third derivatives of {h00,h01,h10,h11}, in the existing ordering.
        // These constants are unchanged by the reflected basis: both its
        // row permutation/sign and the third-order chain-rule sign apply.
        const std::array<long double, 4> third{{12.0L, -12.0L, 6.0L, 6.0L}};
        long double dxxx = 0.0L, dxxy = 0.0L, dxyy = 0.0L, dyyy = 0.0L;
        for (int y = 0; y < 4; ++y)
            for (int x = 0; x < 4; ++x)
            {
                const long double coefficient = cell.data[static_cast<std::size_t>(4 * y + x)];
                dxxx += coefficient * third[x] * by.value[y];
                dxxy += coefficient * bx.second[x] * by.first[y];
                dxyy += coefficient * bx.first[x] * by.second[y];
                dyyy += coefficient * bx.value[x] * third[y];
            }
        auto physical = [&](long double derivative, double first, double second, double thirdScale) {
            long double value = cell.scale * derivative;
            Detail::require(std::isfinite(value) && (derivative == 0.0L || value != 0.0L),
                            "third derivative scaling overflow or underflow");
            for (double length : {first, second, thirdScale})
            {
                const long double divided = value / length;
                Detail::require(std::isfinite(divided) && (value == 0.0L || divided != 0.0L),
                                "third derivative physical scaling overflow or underflow");
                value = divided;
            }
            const double result = Detail::checked(value, "non-finite reconstructed third derivative");
            Detail::require(value == 0.0L || result != 0.0, "third derivative conversion underflow");
            return result;
        };
        ThirdDerivatives result;
        result.dxxx = physical(dxxx, hx_, hx_, hx_);
        result.dxxy = physical(dxxy, hx_, hx_, hy_);
        result.dxyy = physical(dxyy, hx_, hy_, hy_);
        result.dyyy = physical(dyyy, hy_, hy_, hy_);
        return result;
    }

    std::array<double, 7> segmentPolynomial(const Vector2D &localBegin, const Vector2D &localEnd) const
    {
        auto valid = [](const Vector2D &p) {
            return std::isfinite(p.x) && std::isfinite(p.y) && p.x >= 0.0 && p.x <= 1.0 && p.y >= 0.0 &&
                   p.y <= 1.0;
        };
        Detail::require(valid(localBegin) && valid(localEnd),
                        "segment endpoints must lie in the closed local unit cell");
        const auto &cell = cache_;
        std::array<double, 7> result{};
        if (cell.scale == 0.0)
        {
            result[0] = cell.offset;
            return result;
        }
        auto restrictedBasis = [](double start, double end) {
            // Rows match Detail::Basis: endpoint values 0/1, derivatives 0/1.
            const std::array<std::array<int, 4>, 4> coefficients{
                {{{1, 0, -3, 2}}, {{0, 0, 3, -2}}, {{0, 1, -2, 1}}, {{0, 0, -1, 1}}}};
            std::array<std::array<long double, 4>, 4> powers{}, basis{};
            powers[0][0] = 1.0L;
            const long double delta = static_cast<long double>(end) - start;
            for (int power = 1; power <= 3; ++power)
                for (int k = 0; k <= power; ++k)
                {
                    if (k < power)
                        powers[power][k] += start * powers[power - 1][k];
                    if (k > 0)
                        powers[power][k] += delta * powers[power - 1][k - 1];
                }
            for (int row = 0; row < 4; ++row)
                for (int power = 0; power <= 3; ++power)
                    for (int k = 0; k <= power; ++k)
                        basis[row][k] += coefficients[row][power] * powers[power][k];
            return basis;
        };
        const auto bx = restrictedBasis(localBegin.x, localEnd.x);
        const auto by = restrictedBasis(localBegin.y, localEnd.y);
        std::array<long double, 7> coefficients{};
        for (int y = 0; y < 4; ++y)
            for (int x = 0; x < 4; ++x)
                for (int v = 0; v < 4; ++v)
                    for (int u = 0; u < 4; ++u)
                        coefficients[u + v] +=
                            static_cast<long double>(cell.data[4 * y + x]) * bx[x][u] * by[y][v];
        for (int k = 0; k < 7; ++k)
        {
            const long double value = cell.scale * coefficients[k] + (k == 0 ? cell.offset : 0.0L);
            result[k] = Detail::checked(value, "non-finite restricted scalar polynomial coefficient");
            Detail::require(value == 0.0L || result[k] != 0.0, "restricted scalar coefficient underflow");
        }
        return result;
    }

    ExpandedPolynomial segmentPolynomialExact(const Vector2D &localBegin, const Vector2D &localEnd) const
    {
        namespace Exact = ExactExpansion;
        Exact::validateEnvironment();
        auto valid = [](const Vector2D &p) {
            return std::isfinite(p.x) && std::isfinite(p.y) && p.x >= 0.0 && p.x <= 1.0 && p.y >= 0.0 &&
                   p.y <= 1.0;
        };
        Detail::require(valid(localBegin) && valid(localEnd),
                        "segment endpoints must lie in the closed local unit cell");
        const auto bx = restrictedBasisExact(localBegin.x, localEnd.x);
        const auto by = restrictedBasisExact(localBegin.y, localEnd.y);
        ExpandedPolynomial result{};
        for (int y = 0; y < 4; ++y)
            for (int x = 0; x < 4; ++x)
            {
                const auto &node = cellNode(x % 2, y % 2);
                const double dof = y < 2 ? (x < 2 ? node.value : node.x) : (x < 2 ? node.y : node.xy);
                if (dof == 0.0)
                    continue; // EXACT zero, not a magnitude threshold.
                const auto value = Exact::fromDouble(dof);
                for (int v = 0; v < 4; ++v)
                {
                    if (Exact::sign(by[y][v]) == 0)
                        continue;
                    for (int u = 0; u < 4; ++u)
                    {
                        if (Exact::sign(bx[x][u]) == 0)
                            continue;
                        const auto term = Exact::multiply(value, Exact::multiply(bx[x][u], by[y][v]));
                        result[u + v] = Exact::add(result[u + v], term);
                    }
                }
            }
        return result;
    }

    ExpandedPolynomial segmentPolynomialCleared(const RationalSegmentCoordinate &x,
                                                const RationalSegmentCoordinate &y) const
    {
        namespace Exact = ExactExpansion;
        Exact::validateEnvironment();
        const auto validate = [](const RationalSegmentCoordinate &coordinate) {
            Detail::require(Exact::sign(coordinate.denominator) > 0,
                            "cleared restriction denominator must be strictly positive");
            for (const auto *numerator : {&coordinate.numeratorBegin, &coordinate.numeratorEnd})
                Detail::require(Exact::sign(*numerator) >= 0 &&
                                    Exact::sign(Exact::subtract(coordinate.denominator, *numerator)) >= 0,
                                "cleared restriction endpoints must lie in the closed local unit cell");
        };
        validate(x);
        validate(y);
        const auto bx = restrictedBasisCleared(x), by = restrictedBasisCleared(y);
        ExpandedPolynomial result{};
        for (int row = 0; row < 4; ++row)
            for (int column = 0; column < 4; ++column)
            {
                const auto &node = cellNode(column % 2, row % 2);
                const double dof =
                    row < 2 ? (column < 2 ? node.value : node.x) : (column < 2 ? node.y : node.xy);
                if (dof == 0.0)
                    continue;
                const auto value = Exact::fromDouble(dof);
                for (int v = 0; v < 4; ++v)
                {
                    if (Exact::sign(by[row][v]) == 0)
                        continue;
                    for (int u = 0; u < 4; ++u)
                    {
                        if (Exact::sign(bx[column][u]) == 0)
                            continue;
                        const auto term = Exact::multiply(value, Exact::multiply(bx[column][u], by[row][v]));
                        result[u + v] = Exact::add(result[u + v], term);
                    }
                }
            }
        return result;
    }

    Interval cellBounds(Vector2D localLower = {0.0, 0.0}, Vector2D localUpper = {1.0, 1.0}) const
    {
        ExactExpansion::validateEnvironment();
        Detail::require(std::isfinite(localLower.x) && std::isfinite(localLower.y) &&
                            std::isfinite(localUpper.x) && std::isfinite(localUpper.y) &&
                            localLower.x >= 0.0 && localLower.y >= 0.0 && localUpper.x <= 1.0 &&
                            localUpper.y <= 1.0 && localLower.x <= localUpper.x &&
                            localLower.y <= localUpper.y,
                        "bounds box must be an ordered closed subset of the local unit cell");
        using Arithmetic = Detail::Outward;
        std::array<Interval, 16> controls{};
        for (int y = 0; y < 2; ++y)
            for (int x = 0; x < 2; ++x)
            {
                const auto &node = cellNode(x, y);
                const int bx = 3 * x, by = 3 * y, sx = x ? -1 : 1, sy = y ? -1 : 1;
                const auto f = Arithmetic::point(node.value);
                auto dx = Arithmetic::dividePositive(Arithmetic::point(node.x), 3.0);
                auto dy = Arithmetic::dividePositive(Arithmetic::point(node.y), 3.0);
                auto dxy = Arithmetic::dividePositive(Arithmetic::point(node.xy), 9.0);
                if (sx < 0)
                    dx = Arithmetic::negate(dx);
                if (sy < 0)
                    dy = Arithmetic::negate(dy);
                if (sx * sy < 0)
                    dxy = Arithmetic::negate(dxy);
                controls[4 * by + bx] = f;
                controls[4 * by + bx + sx] = Arithmetic::add(f, dx);
                controls[4 * (by + sy) + bx] = Arithmetic::add(f, dy);
                controls[4 * (by + sy) + bx + sx] =
                    Arithmetic::add(Arithmetic::add(Arithmetic::add(f, dx), dy), dxy);
            }
        return Detail::tensorBernsteinBounds(controls, localLower, localUpper);
    }

    int cellSignExceptVertices(unsigned allowedCornerMask) const
    {
        namespace Exact = ExactExpansion;
        Exact::validateEnvironment();
        Detail::require((allowedCornerMask & ~15u) == 0u, "invalid allowed-corner mask");
        std::array<int, 16> signs{};
        const auto scaled = [](double value, double factor) {
            return Exact::multiply(Exact::fromDouble(value), Exact::fromDouble(factor));
        };
        for (int y = 0; y < 2; ++y)
            for (int x = 0; x < 2; ++x)
            {
                const auto &node = cellNode(x, y);
                const int bx = 3 * x, by = 3 * y, sx = x ? -1 : 1, sy = y ? -1 : 1;
                const auto f = scaled(node.value, 9.0);
                const auto dx = scaled(node.x, 3.0 * sx), dy = scaled(node.y, 3.0 * sy);
                const auto dxy = scaled(node.xy, static_cast<double>(sx * sy));
                signs[4 * by + bx] = Exact::sign(f);
                signs[4 * by + bx + sx] = Exact::sign(Exact::add(f, dx));
                signs[4 * (by + sy) + bx] = Exact::sign(Exact::add(f, dy));
                signs[4 * (by + sy) + bx + sx] =
                    Exact::sign(Exact::add(Exact::add(Exact::add(f, dx), dy), dxy));
            }
        return certifyBernsteinSigns(signs, allowedCornerMask);
    }

    // Shared Bernstein certificate; common positive denominators do not
    // affect signs. Native leaves retain expanded nodal numerators here.
    static int expandedCellSignExceptVertices(const std::array<ExactExpansion::Expansion, 16> &dofs,
                                              unsigned allowedCornerMask)
    {
        Detail::require((allowedCornerMask & ~15u) == 0u, "invalid allowed-corner mask");
        return certifyBernsteinSigns(expandedBernsteinSigns(dofs), allowedCornerMask);
    }

    static std::array<ExactExpansion::Expansion,16> bernsteinControlsCleared(
        const std::array<ExactExpansion::Expansion,16>& dofs)
    {
        namespace Exact = ExactExpansion;
        std::array<Exact::Expansion,16> controls;
        for (int y = 0; y < 2; ++y)
            for (int x = 0; x < 2; ++x)
            {
                const int bx = 3 * x, by = 3 * y, sx = x ? -1 : 1, sy = y ? -1 : 1;
                const auto f = Exact::multiply(dofs[4 * y + x], Exact::fromDouble(9.0));
                const auto dx = Exact::multiply(dofs[4 * y + x + 2], Exact::fromDouble(3.0 * sx));
                const auto dy = Exact::multiply(dofs[4 * (y + 2) + x], Exact::fromDouble(3.0 * sy));
                const auto dxy = Exact::multiply(dofs[4 * (y + 2) + x + 2], Exact::fromDouble(sx * sy));
                controls[4 * by + bx] = f;
                controls[4 * by + bx + sx] = Exact::add(f, dx);
                controls[4 * (by + sy) + bx] = Exact::add(f, dy);
                controls[4 * (by + sy) + bx + sx] = Exact::add(Exact::add(Exact::add(f, dx), dy), dxy);
            }
        return controls; // Nine times the homogeneous Bernstein numerator.
    }

    static std::array<int,16> expandedBernsteinSigns(const std::array<ExactExpansion::Expansion,16>& dofs)
    {
        const auto controls=bernsteinControlsCleared(dofs);std::array<int,16> signs{};
        for (std::size_t i=0;i<signs.size();++i) signs[i]=ExactExpansion::sign(controls[i]);
        return signs;
    }

    template<std::size_t Count>
    static int certifyBernsteinSigns(const std::array<int, Count> &signs, unsigned allowedCornerMask)
    {
        static_assert(Count == 16 || Count == 25, "supported cubic/quartic tensor controls");
        constexpr int Side = Count == 16 ? 4 : 5;
        int sign = 0;
        for (int value : signs)
            if (value != 0)
            {
                if (sign != 0 && sign != value)
                    return 0;
                sign = value;
            }
        if (sign == 0)
            return 0;
        const int corners[4] = {0, Side-1, Side*Side-1, Side*(Side-1)}; // SW, SE, NE, NW.
        for (unsigned k = 0; k < 4; ++k)
            if (signs[corners[k]] == 0 && (allowedCornerMask & (1u << k)) == 0u)
                return 0;
        for (int edge = 0; edge < 4; ++edge)
        {
            bool nonzero = false;
            for (int k = 0; k < Side; ++k)
            {
                const int control = edge == 0 ? k : edge == 1 ? Side*k+Side-1
                    : edge == 2 ? Side*(Side-1)+k : Side*k;
                nonzero = nonzero || signs[control] != 0;
            }
            if (!nonzero)
                return 0;
        }
        // All Bernstein bases are positive in the open cell/edge interiors.
        // Weak common sign and one strict control exclude zeros there; only
        // the explicitly allowed corners can remain zero.
        return sign;
    }

    Interval evaluateEnclosureCell(double localX, double localY) const
    {
        ExactExpansion::validateEnvironment();
        Detail::require(std::isfinite(localX) && std::isfinite(localY) && localX >= 0.0 && localX <= 1.0 &&
                            localY >= 0.0 && localY <= 1.0,
                        "local cell coordinates must be finite and in [0,1]");
        using Arithmetic = Detail::Outward;
        const auto bx = Arithmetic::basis(localX), by = Arithmetic::basis(localY);
        Interval result{};
        for (int y = 0; y < 4; ++y)
            for (int x = 0; x < 4; ++x)
            {
                const auto &node = cellNode(x % 2, y % 2);
                const double dof = y < 2 ? (x < 2 ? node.value : node.x) : (x < 2 ? node.y : node.xy);
                if (dof == 0.0 || Arithmetic::zero(bx[x]) || Arithmetic::zero(by[y]))
                    continue;
                result = Arithmetic::add(
                    result, Arithmetic::multiply(Arithmetic::point(dof), Arithmetic::multiply(bx[x], by[y])));
            }
        return result;
    }

    static CellData makeCell(const std::array<const Node *, 4> &corners)
    {
        CellData cell;
        cell.offset = corners[0]->value;
        for (int y = 0; y < 2; ++y)
            for (int x = 0; x < 2; ++x)
            {
                const auto &node = *corners[static_cast<std::size_t>(2 * y + x)];
                cell.data[4 * y + x] = Detail::checked(static_cast<long double>(node.value) - cell.offset,
                                                       "unrepresentable local scalar difference");
                cell.data[4 * y + x + 2] = node.x;
                cell.data[4 * (y + 2) + x] = node.y;
                cell.data[4 * (y + 2) + x + 2] = node.xy;
            }
        for (double value : cell.data)
            cell.scale = std::max(cell.scale, std::abs(value));
        Detail::require(std::isfinite(cell.scale), "non-finite cell normalization scale");
        if (cell.scale != 0.0)
            for (double &value : cell.data)
            {
                const double normalized = value / cell.scale;
                Detail::require(std::isfinite(normalized) && (value == 0.0 || normalized != 0.0),
                                "Hermite data normalization underflow; unsupported numeric scale");
                value = normalized;
            }
        return cell;
    }

  private:
    std::array<const Node *, 4> corners_;
    const CellData &cache_;
    double hx_, hy_;
    const Node &cellNode(int x, int y) const
    {
        return *corners_[2 * y + x];
    }

  public:
    // The leaf owner reuses these exact restrictions when constraining a
    // hanging node. They are the existing Hermite basis, not another fit.
    using ExpandedBasis = std::array<std::array<ExactExpansion::Expansion, 4>, 4>;
    static ExpandedBasis restrictedBasisCleared(const RationalSegmentCoordinate &coordinate,
                                                unsigned derivativeOrder = 0)
    {
        namespace Exact = ExactExpansion;
        // Homogenize each existing integer basis coefficient c_k u^k as
        // c_k N(t)^k D^(3-k). Accumulate powers, not sampled/interpolated
        // values; exact product tails and unsupported-range failures remain.
        Detail::require(derivativeOrder <= 3, "Hermite basis derivative order exceeds degree three");
        std::array<std::array<int, 4>, 4> coefficients{
            {{{1, 0, -3, 2}}, {{0, 0, 3, -2}}, {{0, 1, -2, 1}}, {{0, 0, -1, 1}}}};
        for (unsigned order = 0; order < derivativeOrder; ++order)
            for (auto &row : coefficients)
            {
                for (int k = 0; k < 3; ++k)
                    row[k] = (k + 1) * row[k + 1];
                row[3] = 0;
            }
        ExpandedBasis powers{}, basis{};
        std::array<Exact::Expansion, 4> denominatorPowers;
        denominatorPowers[0] = Exact::fromDouble(1.0);
        for (int power = 1; power <= 3; ++power)
            denominatorPowers[power] = Exact::multiply(denominatorPowers[power - 1], coordinate.denominator);
        const auto delta = Exact::subtract(coordinate.numeratorEnd, coordinate.numeratorBegin);
        // A point or axis-aligned restriction has exactly zero parameter
        // degree in this direction. Keep the constant coefficient's original
        // arithmetic; the omitted higher coefficients are exact empty zeros.
        const int parameterDegree = Exact::sign(delta) == 0 ? 0 : 3;
        powers[0][0] = Exact::fromDouble(1.0);
        for (int power = 1; power <= 3; ++power)
            for (int k = 0; k <= std::min(power, parameterDegree); ++k)
            {
                if (k < power)
                    powers[power][k] = Exact::add(
                        powers[power][k], Exact::multiply(coordinate.numeratorBegin, powers[power - 1][k]));
                if (k > 0)
                    powers[power][k] =
                        Exact::add(powers[power][k], Exact::multiply(delta, powers[power - 1][k - 1]));
            }
        for (int row = 0; row < 4; ++row)
            for (int power = 0; power <= 3; ++power)
            {
                if (coefficients[row][power] == 0)
                    continue;
                const auto factor =
                    Exact::multiply(Exact::fromDouble(static_cast<double>(coefficients[row][power])),
                                    denominatorPowers[3 - power]);
                for (int k = 0; k <= std::min(power, parameterDegree); ++k)
                    basis[row][k] = Exact::add(basis[row][k], Exact::multiply(factor, powers[power][k]));
            }
        return basis;
    }

    static ExpandedBasis restrictedBasisExact(double start, double end)
    {
        namespace Exact = ExactExpansion;
        // The SAME four Hermite integer polynomials as Detail::Basis and the
        // compatibility segmentPolynomial entry point; no interpolated fit.
        const std::array<std::array<int, 4>, 4> coefficients{
            {{{1, 0, -3, 2}}, {{0, 0, 3, -2}}, {{0, 1, -2, 1}}, {{0, 0, -1, 1}}}};
        ExpandedBasis powers{}, basis{};
        const auto origin = Exact::fromDouble(start);
        const auto delta = Exact::subtract(Exact::fromDouble(end), origin);
        powers[0][0] = Exact::fromDouble(1.0);
        for (int power = 1; power <= 3; ++power)
            for (int k = 0; k <= power; ++k)
            {
                if (k < power)
                    powers[power][k] =
                        Exact::add(powers[power][k], Exact::multiply(origin, powers[power - 1][k]));
                if (k > 0)
                    powers[power][k] =
                        Exact::add(powers[power][k], Exact::multiply(delta, powers[power - 1][k - 1]));
            }
        for (int row = 0; row < 4; ++row)
            for (int power = 0; power <= 3; ++power)
            {
                if (coefficients[row][power] == 0)
                    continue;
                const auto coefficient = Exact::fromDouble(static_cast<double>(coefficients[row][power]));
                for (int k = 0; k <= power; ++k)
                    basis[row][k] = Exact::add(basis[row][k], Exact::multiply(coefficient, powers[power][k]));
            }
        return basis;
    }
};
} // namespace LevelSetHermiteScalar
} // namespace
} // namespace LAMP::Geometry::Detail
