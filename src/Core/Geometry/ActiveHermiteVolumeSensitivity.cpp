#include "Core/Geometry/ActiveHermiteVolumeSensitivity.h"
#include "Core/Geometry/Detail/ExactExpansion.h"

#include <cmath>
#include <stdexcept>

namespace ActiveHermiteVolumeSensitivity
{
namespace
{
void require(bool condition, const char* reason)
{
    if (!condition) throw std::invalid_argument(reason);
}

long double multiply(long double a, long double b)
{
    const long double result = a * b;
    if (!std::isfinite(result) || (a != 0.0L && b != 0.0L && result == 0.0L))
        throw std::overflow_error("[ActiveHermiteVolumeSensitivity] unrepresentable Hermite basis product");
    return result;
}

// The same tensor cubic Hermite basis as the shared scalar, written with
// endpoint tails factored so a small positive endpoint basis is not cancelled.
struct Basis
{
    std::array<long double, 2> value, derivative;
    explicit Basis(double coordinate)
    {
        const long double t = coordinate, u = 1.0L - t;
        const long double t2 = multiply(t, t), u2 = multiply(u, u);
        value = {{multiply(u2, 1.0L + 2.0L * t), multiply(t2, 1.0L + 2.0L * u)}};
        derivative = {{multiply(t, u2), -multiply(t2, u)}};
    }
};

long double weightedBasis(const std::array<long double, 5>& factors, long double gradient)
{
    // Accumulate exponents separately: tiny arc weights must not disappear
    // in an intermediate product before division by the actual gradient.
    long double mantissa = 1.0L;
    int exponent = 0;
    for (long double factor : factors)
    {
        if (factor == 0.0L) return 0.0L; // An exact endpoint basis zero.
        int part = 0;
        mantissa *= std::frexp(factor, &part);
        exponent += part;
    }
    int denominatorExponent = 0;
    mantissa /= std::frexp(gradient, &denominatorExponent);
    const long double result = -std::scalbn(mantissa, exponent - denominatorExponent);
    if (!std::isfinite(result) || result == 0.0L)
        throw std::overflow_error("[ActiveHermiteVolumeSensitivity] unrepresentable weighted sensitivity");
    return result;
}
}

Coefficients assemble(int cellI, int cellJ, Vector2D cellSize,
    const std::vector<LevelSetInterfacePanels::ArcSample>& samples)
{
    LAMP::Geometry::Detail::ExactExpansion::validateEnvironment();
    require(cellI >= 0 && cellJ >= 0 && std::isfinite(cellSize.x) && cellSize.x > 0.0
        && std::isfinite(cellSize.y) && cellSize.y > 0.0,
        "[ActiveHermiteVolumeSensitivity] invalid cell or physical spacing");
    Coefficients result{}, correction{};
    constexpr std::array<std::array<int, 2>, 4> corners{{{{0, 0}}, {{1, 0}}, {{1, 1}}, {{0, 1}}}};
    const auto unit = [](double t) { return std::isfinite(t) && t >= 0.0 && t <= 1.0; };
    for (const auto& sample : samples)
    {
        const auto& frame = sample.frame;
        require(frame.hasCellLocalPoint && frame.cellI == cellI && frame.cellJ == cellJ
            && unit(frame.cellLocalPoint.x) && unit(frame.cellLocalPoint.y)
            && unit(sample.localGraphCoordinate)
            && std::isfinite(frame.point.x) && std::isfinite(frame.point.y),
            "[ActiveHermiteVolumeSensitivity] sample has no matching cell-local provenance");
        require(std::isfinite(sample.weight) && sample.weight > 0.0
            && std::isfinite(sample.physicalWeight) && sample.physicalWeight > 0.0
            && std::isfinite(sample.gradientMagnitude) && sample.gradientMagnitude > 0.0,
            "[ActiveHermiteVolumeSensitivity] sample measure/gradient must be positive and finite");
        const Basis x(frame.cellLocalPoint.x), y(frame.cellLocalPoint.y);
        for (std::size_t corner = 0; corner < corners.size(); ++corner)
            for (int dof = 0; dof < 4; ++dof)
            {
                const int ix = corners[corner][0], iy = corners[corner][1];
                const long double term = weightedBasis({{sample.physicalWeight,
                    dof & 1 ? x.derivative[ix] : x.value[ix],
                    dof & 2 ? y.derivative[iy] : y.value[iy],
                    dof & 1 ? cellSize.x : 1.0L, dof & 2 ? cellSize.y : 1.0L}}, sample.gradientMagnitude);
                const std::size_t column = 4 * corner + dof;
                const long double corrected = term - correction[column];
                const long double next = result[column] + corrected;
                if (!std::isfinite(corrected) || !std::isfinite(next))
                    throw std::overflow_error("[ActiveHermiteVolumeSensitivity] nonfinite sensitivity sum");
                correction[column] = (next - result[column]) - corrected;
                result[column] = next;
            }
    }
    return result;
}
}
