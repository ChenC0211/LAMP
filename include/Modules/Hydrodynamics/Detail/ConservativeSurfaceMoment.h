#pragma once

#include "Core/Geometry/CellChartPoint.h"
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace Quantum::Hydrodynamics::Detail
{
// Unit physical-face-flux basis on one of its two owner cells, divided by
// center spacing. In RZ the radial basis has constant cylindrical divergence.
// Endpoint-relative distances preserve small positive integration weights;
// these terminal arithmetic values never determine interface topology.
inline long double surfaceFaceBasis(Geometry::CellChartPoint chart,
    const std::array<Vector2D,2>& bounds,long double volume,double spacing,
    bool horizontal,bool negativeOwner,bool axisymmetric)
{
    Geometry::validateCellChart(chart);
    if(!(volume>0 && spacing>0))throw std::invalid_argument("invalid surface moment volume/spacing");
    const auto coordinate=horizontal?chart.y:chart.x;
    const long double width=horizontal?bounds[1].y-bounds[0].y:bounds[1].x-bounds[0].x;
    const long double distance=width*(coordinate.fromUpper==negativeOwner
        ?1.L-coordinate.offset:coordinate.offset);
    long double value=distance/volume/spacing;
    if(axisymmetric && !horizontal) {
        const long double radius=chart.x.fromUpper?bounds[1].x-width*chart.x.offset
            :bounds[0].x+width*chart.x.offset;
        if(!(radius>0))throw std::invalid_argument("radial surface moment needs an interior quadrature radius");
        value*=(radius+(negativeOwner?bounds[0].x:bounds[1].x))/(2*radius);
    }
    if(!(value>=0) || !std::isfinite(value))throw std::invalid_argument("invalid surface moment basis");
    return value;
}
struct SurfaceMomentNode {long double weight=0;double normal=0,tangent=0;};
struct SurfaceMomentBalance
{
    std::vector<long double> normal,tangent;
    long double normalDefect=0,tangentDefect=0;
};
// Minimum weighted L2 correction of the two oriented moments. A constant
// traction then uses exactly the caller's pressure-compatible C gradients.
// This algebra does not certify variable-curvature accuracy or global force
// conservation; both remain properties of the geometry and physical model.
inline SurfaceMomentBalance balanceSurfaceMoments(const std::vector<SurfaceMomentNode>& nodes,
    long double normalTarget,long double tangentTarget)
{
    long double area=0,normal=0,tangent=0;
    for(const auto& q:nodes) {
        if(!(q.weight>=0) || !std::isfinite(q.weight) || !std::isfinite(q.normal) || !std::isfinite(q.tangent))
            throw std::invalid_argument("invalid oriented surface moment");
        area+=q.weight;normal+=q.weight*q.normal;tangent+=q.weight*q.tangent;
    }
    if(!(area>0) || !std::isfinite(area) || !std::isfinite(normalTarget) || !std::isfinite(tangentTarget))
        throw std::invalid_argument("empty/nonfinite surface moment support");
    SurfaceMomentBalance result;result.normalDefect=normalTarget-normal;result.tangentDefect=tangentTarget-tangent;
    for(const auto& q:nodes) {
        result.normal.push_back(q.weight*(q.normal+result.normalDefect/area));
        result.tangent.push_back(q.weight*(q.tangent+result.tangentDefect/area));
    }
    return result;
}
}
