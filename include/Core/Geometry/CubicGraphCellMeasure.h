#pragma once
#include "Core/Geometry/PlanarCellCut.h"
#include <array>

namespace LAMP::Geometry::CubicGraphCellMeasure
{
// Rotated scaled graph n=sum(c[k]*t^k), liquid n<=H(t).
// t=(-ny,nx).(x-origin)/scale, n=normal.(x-origin)/scale.
// Roots reuse the bounded polynomial isolator: unresolved contacts reject.
// Polynomial volume moments are integrated in local slab coordinates.
// Root enclosures refer to the assembled floating polynomial; this is not
// a certificate for the combined frame-conversion and volume-rounding error.
struct Result {long double volume=0;std::array<long double,4> derivative{};};
Result evaluate(const PlanarCellCut::Rectangle&,PlanarCellCut::Point origin,
    Vector2D normal,long double scale,const std::array<long double,4>& coefficients,bool axisymmetric);
}
