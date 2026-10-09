#pragma once

#include "Core/Geometry/LevelSetInterfacePanels.h"
#include "Core/Mesh/CartesianMeasure.h"
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace Quantum::Geometry::Detail
{
namespace
{
namespace InterfaceQuadrature
{
inline void require(bool valid, const char *message)
{
    if (!valid)
        throw std::runtime_error(std::string("[LevelSetInterfacePanels] ") + message);
}
inline double checked(long double value, const char *reason)
{
    const double result = static_cast<double>(value);
    require(std::isfinite(value) && std::isfinite(result) && (value == 0.0L || result != 0.0), reason);
    return result;
}
inline double offAxisWeight(const ::LevelSetHermiteGeometry::Frame &frame)
{
    require(std::isfinite(frame.point.x) && frame.point.x > 0.0,
            "axisymmetric panel interiors and non-axis endpoints require r>0");
    return CartesianMeasure::radialWeight(frame.point.x);
}
constexpr std::array<double, 8> gaussX{{-.96028985649753623168, -.79666647741362673959,
                                        -.52553240991632898582, -.18343464249564980494, .18343464249564980494,
                                        .52553240991632898582, .79666647741362673959, .96028985649753623168}};
constexpr std::array<double, 8> gaussW{{.10122853629037625915, .22238103445337447054, .31370664587788728734,
                                        .36268378337836198297, .36268378337836198297, .31370664587788728734,
                                        .22238103445337447054, .10122853629037625915}};

// Original subarc GL8 operation order, shared by the two source layouts.
// Frame and coordinate callbacks preserve their own explicit provenance;
// this routine neither reconstructs a scalar nor projects a query point.
template <class At, class Coordinate, class Metric, class Gradient>
std::vector<LevelSetInterfacePanels::ArcSample> sampleGraph(long double span, bool axisymmetric, At at,
                                                            Coordinate coordinate, Metric metric,
                                                            Gradient gradient)
{
    require(std::isfinite(span) && span > 0.0L, "subarc has no representable positive physical span");
    std::vector<LevelSetInterfacePanels::ArcSample> result;
    result.reserve(gaussX.size());
    for (std::size_t k = 0; k < gaussX.size(); ++k)
    {
        LevelSetInterfacePanels::ArcSample node;
        const double fraction = (1.0 + gaussX[k]) / 2.0;
        node.frame = at(fraction);
        coordinate(node);
        node.weight =
            checked(span * gaussW[k] * metric(node.frame) / 2.0L, "unrepresentable subarc quadrature weight");
        node.physicalWeight = axisymmetric
                                  ? checked(static_cast<long double>(offAxisWeight(node.frame)) * node.weight,
                                            "unrepresentable subarc physical quadrature weight")
                                  : node.weight;
        const auto raw = gradient(node.frame);
        node.gradientMagnitude =
            checked(std::hypot(static_cast<long double>(raw.dx), static_cast<long double>(raw.dy)),
                    "unrepresentable subarc scalar gradient");
        require(node.weight > 0.0 && node.physicalWeight > 0.0 && node.gradientMagnitude > 0.0,
                "subarc quadrature requires positive finite weights and scalar gradient");
        result.push_back(std::move(node));
    }
    return result;
}
} // namespace InterfaceQuadrature
} // namespace
} // namespace Quantum::Geometry::Detail
