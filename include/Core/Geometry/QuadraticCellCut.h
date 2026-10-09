#pragma once
#include "Core/Geometry/PlanarCellCut.h"
#include <optional>

namespace Quantum::Geometry::QuadraticCellCut
{
using Rectangle=PlanarCellCut::Rectangle;
using Point=PlanarCellCut::Point;
// Cell-centred graph n = offset + slope*t + second*t*t/2.
// axis is the normal coordinate (0=x/r, 1=y/z); sign*(n-H)<=0 is liquid.
// Optional frameNormal replaces the coordinate axes by its unit normal n and
// t=(-n.y,n.x). The same sign*(n-H)<=0 convention applies. Coefficients remain
// physical derivatives about the cell centre; RZ still measures 2*pi*r dA.
struct Shape {int axis=1,sign=1;long double slope=0,second=0;
    std::optional<Vector2D> frameNormal;};
struct Basis {Point tangent,normal;};
Basis basis(const Shape&);
struct Reconstruction {Shape shape;long double offset=0,volume=0,residual=0;};
// Same centred phase polynomial used by interval and temporal-event queries.
long double implicitValue(const Rectangle&,const Reconstruction&,Point);
// Rotated, scaled graph n=c[0]+c[1]*t+c[2]*t^2, liquid below the graph.
// (t,n) are coordinates in the tangent=(-ny,nx)/normal frame, divided
// by scale. The supplied normal is normalized internally. The derivatives
// are with respect to these three dimensionless coefficients, including
// the actual XY/RZ cell measure and moving edge intersections.
struct GraphMeasure {long double volume=0;std::array<long double,3> derivative{};};
GraphMeasure graphMeasure(const Rectangle&,Point origin,Vector2D normal,long double scale,
    const std::array<long double,3>& coefficients,bool axisymmetric);
// Closed polynomial moments after floating quadratic root splitting; this is
// not an interval certificate. The inverse enforces the original PLIC absolute
// volume tolerance and never changes the requested fraction.
Reconstruction reconstruct(const Rectangle&,Shape,long double fraction,bool axisymmetric);
// Project an entirely interior single graph onto a parabola preserving its
// two endpoint heights and the supplied XY/RZ liquid volume. Heights are
// absolute physical normal coordinates at the tangent bounds. Returns no
// graph when these data admit no contained regular parabola; no clipping.
std::optional<Reconstruction> projectEndpointVolume(const Rectangle&,int axis,int sign,
    long double lowerHeight,long double upperHeight,long double liquidVolume,bool axisymmetric);
PlanarCellCut::Overlap overlapInFrame(const Rectangle&,const std::vector<Point>& polygon,Point origin,
    const Reconstruction&,bool axisymmetric);
// Liquid portions of a retained chord interval, in its original [0,1]
// parameter. Caller has already clipped that interval to the donor rectangle.
std::vector<std::array<long double,2>> intervals(const Rectangle&,const Reconstruction&,
    Point begin,Point end,long double lower=0,long double upper=1);
}
