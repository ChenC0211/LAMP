#pragma once

#include <cmath>
#include <stdexcept>

// Measures on a Cartesian (x,y) grid, optionally revolved with x=r and y=z.
// Input lengths/areas remain the actual 2D geometry. These helpers neither
// change coordinates nor clip a negative radius onto the symmetry axis.
namespace CartesianMeasure
{
inline constexpr double twoPi = 6.283185307179586476925286766559005768;

inline double radialWeight(double radius)
{
    if (!std::isfinite(radius) || radius < 0.0)
        throw std::runtime_error("[CartesianMeasure] radius must be finite and nonnegative");
    if (radius == 0.0) return 0.0;
    const double weight = twoPi * radius;
    if (!std::isfinite(weight) || weight <= 0.0)
        throw std::runtime_error("[CartesianMeasure] radial weight is unrepresentable");
    return weight;
}

inline double face(double planarLength, double radius, bool axisymmetric)
{
    if (!axisymmetric) return planarLength;
    if (!std::isfinite(planarLength) || planarLength <= 0.0)
        throw std::runtime_error("[CartesianMeasure] face length must be finite and positive");
    const double weight = radialWeight(radius);
    if (weight == 0.0) return 0.0; // Exactly zero area only on the actual axis.
    const double result = weight * planarLength;
    if (!std::isfinite(result) || result <= 0.0)
        throw std::runtime_error("[CartesianMeasure] revolved face area is unrepresentable");
    return result;
}

inline double cell(double planarArea, double centerRadius, bool axisymmetric)
{
    if (!axisymmetric) return planarArea;
    if (!std::isfinite(centerRadius) || centerRadius <= 0.0)
        throw std::runtime_error("[CartesianMeasure] revolved cell center must have positive radius");
    // For a rectangular (r,z) cell the midpoint rule integrates 2*pi*r
    // exactly in real arithmetic, including the first cell adjacent to r=0.
    return face(planarArea, centerRadius, true);
}
} // namespace CartesianMeasure
