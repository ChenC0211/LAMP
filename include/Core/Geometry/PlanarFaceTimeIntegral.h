#pragma once
#include "Core/Geometry/PlanarCellCut.h"
#include <cstddef>
#include <vector>

namespace LAMP::Geometry::PlanarFaceTime
{
using Point=PlanarCellCut::Point;
using AnchoredPoint=PlanarCellCut::AnchoredPoint;
struct Reflection {int signX=1,signY=1;Point translation{};};
struct Bounds {long double x0=0,x1=0,y0=0,y1=0;};
struct Gross
{
    long double halfPositive=0,halfNegative=0,fullPositive=0,fullNegative=0;
    std::size_t temporalEvents=0,temporalIntervals=0;
    unsigned maximumRootRefinements=0;
    long double maximumRelativeUncertainty=0;
};
// Original actual points and complete physical clock are retained before any
// differences/products. Gross signs use canonical +x/+y native normals and
// therefore do not depend on face endpoint order. Only axis-aligned faces.
Gross integratePolygon(const std::vector<Point>& polygon,Point faceBegin,Point faceEnd,
    Point velocityBegin,Point velocityEnd,long double timeBegin,long double timeEnd,
    Point frameOrigin,bool axisymmetric,unsigned eventLimit=512,Reflection reflection={});
// Each polygon coordinate is the exact sum of its public anchor and offset.
// Face, clock, frame and reflection retain the original Point contract; the
// shared engine retains its existing final scalar conversion semantics.
Gross integrateAnchoredPolygon(const std::vector<AnchoredPoint>& polygon,Point faceBegin,Point faceEnd,
    Point velocityBegin,Point velocityEnd,long double timeBegin,long double timeEnd,
    Point frameOrigin,bool axisymmetric,unsigned eventLimit=512,Reflection reflection={});
// Native total measure uses the same exact linear reversal events.
Gross integrateNormal(long double velocityBegin,long double velocityEnd,
    long double timeBegin,long double timeEnd,long double nativeArea);
// Outward enclosure of the actual face's quadratic departure trajectory.
Bounds departureBounds(Point faceBegin,Point faceEnd,Point velocityBegin,Point velocityEnd,
    long double timeBegin,long double timeEnd);
// Exact support guard, independent of outward query rounding. The permitted
// one-reflection domain is [2*lower-upper,2*upper-lower] in each coordinate.
void validateDepartureSupport(Point faceBegin,Point faceEnd,Point velocityBegin,Point velocityEnd,
    long double timeBegin,long double timeEnd,Bounds originalDomainBounds);
}
