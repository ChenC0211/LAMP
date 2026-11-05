#pragma once

// Shared SW--NE piecewise-linear phase-area snapshot. Corner order is
// SW, SE, NE, NW; phi<=0 is the negative phase. This is NOT transported VOF.
// The historical arithmetic, guards and exception strings are retained.
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace LevelSetPhaseFraction
{
inline void require(bool condition, const char* reason)
{
    if (!condition) throw std::runtime_error(std::string("[LevelSetPhaseFraction] ") + reason);
}

namespace Detail
{
struct Vertex
{
    double x, y, phi;
};

inline double negativeArea(const std::array<Vertex, 3>& triangle)
{
    require(!(triangle[0].phi == 0.0 && triangle[1].phi == 0.0 && triangle[2].phi == 0.0),
        "ambiguous zero-valued level-set triangle");
    std::vector<Vertex> clipped;
    clipped.reserve(4);
    for (int i = 0; i < 3; ++i)
    {
        const auto& start = triangle[i];
        const auto& end = triangle[(i + 1) % 3];
        const bool insideStart = start.phi <= 0.0, insideEnd = end.phi <= 0.0;
        if (insideStart != insideEnd)
        {
            const double difference = start.phi - end.phi;
            require(std::isfinite(difference) && difference != 0.0, "unrepresentable level-set edge difference");
            const double fraction = start.phi / difference;
            clipped.push_back({start.x + fraction * (end.x - start.x),
                start.y + fraction * (end.y - start.y), 0.0});
        }
        if (insideEnd) clipped.push_back(end);
    }
    double cross = 0.0;
    for (std::size_t i = 0; i < clipped.size(); ++i)
    {
        const auto& a = clipped[i];
        const auto& b = clipped[(i + 1) % clipped.size()];
        cross += a.x * b.y - b.x * a.y;
    }
    const double area = 0.5 * std::abs(cross);
    require(std::isfinite(area) && area >= 0.0 && area <= 0.5 + 32.0 * std::numeric_limits<double>::epsilon(),
        "invalid clipped triangle area");
    return area;
}

}

inline double positiveFraction(const std::array<double, 4>& phi)
{
    for (double value : phi) require(std::isfinite(value), "non-finite nodal level set");
    const double negative = Detail::negativeArea({{{0.0, 0.0, phi[0]}, {1.0, 0.0, phi[1]}, {1.0, 1.0, phi[2]}}})
        + Detail::negativeArea({{{0.0, 0.0, phi[0]}, {1.0, 1.0, phi[2]}, {0.0, 1.0, phi[3]}}});
    require(std::isfinite(negative) && negative >= 0.0 && negative <= 1.0 + 64.0 * std::numeric_limits<double>::epsilon(),
        "invalid geometric phase fraction");
    return std::clamp(1.0 - negative, 0.0, 1.0);
}

// Compact C2 load indicator, separate from the geometric PL area fraction.
// Its derivative is 15*(1-s*s)^2/(16*halfWidth) inside |s|<1,
// nonnegative with unit integral. Width has the units of phi; it is a physical
// distance half-width only when phi is a resolved signed-distance field.
inline double regularizedPositiveIndicator(double phi, double halfWidth)
{
    if (!std::isfinite(phi) || !std::isfinite(halfWidth) || halfWidth <= 0.0)
        throw std::runtime_error("[LevelSetPhaseFraction] finite phi and positive finite load half-width required");
    if (phi <= -halfWidth) return 0.0;
    if (phi >= halfWidth) return 1.0;
    const long double s = static_cast<long double>(phi) / halfWidth;
    // Factored lower tail avoids subtracting nearly equal halves near -1.
    // This is the same quintic .5+s*(15-10*s*s+3*s*s*s*s)/16.
    const long double z = (1.0L - std::abs(s)) / 2.0L;
    const long double lower = z * z * z * (10.0L + z * (-15.0L + 6.0L * z));
    const double value = static_cast<double>(s < 0.0L ? lower : 1.0L - lower);
    if (!std::isfinite(value) || value < 0.0 || value > 1.0)
        throw std::runtime_error("[LevelSetPhaseFraction] regularized load indicator outside [0,1]");
    return value;
}

} // namespace LevelSetPhaseFraction
