#include "Core/Geometry/PlanarCellCut.h"
#include "Core/Geometry/Detail/ExactExpansion.h"
#include "Core/Geometry/Detail/ExactRationalRound.h"
#include "Core/Geometry/Detail/FiniteVolumeClosure.h"
#include "Core/Mesh/CartesianMeasure.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace LAMP::Geometry::PlanarCellCut
{
namespace
{
namespace Exact=LAMP::Geometry::Detail::ExactExpansion;
using Expanded=Exact::Expansion;
// Expand the ACTUAL published coordinate, including long-double tail bits on
// platforms where Point has more precision than binary64. No rational edge
// intersection or separately represented geometry is measured here.
Expanded exactNumber(long double value)
{
    Expanded result;
    while (value!=0) {
        const double part=static_cast<double>(value);
        if (!std::isfinite(part) || part==0)
            throw std::runtime_error("[PlanarCellCut] coordinate outside exact arithmetic scale");
        Exact::accumulate(result,Exact::fromDouble(part));value-=part;
    }
    return result;
}
Expanded exactCross(Point a,Point b,Point c)
{
    const auto ax=exactNumber(a.x),ay=exactNumber(a.y);
    const auto bx=Exact::subtract(exactNumber(b.x),ax),by=Exact::subtract(exactNumber(b.y),ay);
    const auto cx=Exact::subtract(exactNumber(c.x),ax),cy=Exact::subtract(exactNumber(c.y),ay);
    return Exact::subtract(Exact::multiply(bx,cy),Exact::multiply(by,cx));
}
Expanded positiveTwiceArea(Point a,Point b,Point c,int orientation)
{
    auto twice=exactCross(a,b,c);
    if (orientation<0) twice=Exact::subtract(Expanded{},twice);
    if (Exact::sign(twice)<0)
        throw std::runtime_error("[PlanarCellCut] polygon fan has a negative actual triangle");
    return twice;
}
Expanded coordinateSum(long double a,long double b,long double c)
{ return Exact::add(Exact::add(exactNumber(a),exactNumber(b)),exactNumber(c)); }
long double dividedValue(const Expanded& numerator,long double denominator)
{
    long double value=0;
    for (double part:numerator.components()) value+=part;
    value/=denominator;
    if (!std::isfinite(value) || (value==0 && Exact::sign(numerator)!=0))
        throw std::runtime_error("[PlanarCellCut] exact polygon moment is unrepresentable");
    return value;
}
struct Shape { long double hx, hy, radius, nx, ny; };
Shape shape(const Rectangle& cell, Vector2D normal, bool rz)
{
    if (!std::isfinite(cell.lower.x) || !std::isfinite(cell.lower.y)
        || !std::isfinite(cell.upper.x) || !std::isfinite(cell.upper.y)
        || cell.upper.x <= cell.lower.x || cell.upper.y <= cell.lower.y
        || (rz && cell.lower.x < 0) || !std::isfinite(normal.x) || !std::isfinite(normal.y))
        throw std::invalid_argument("[PlanarCellCut] invalid rectangle/normal");
    const long double length = std::hypot(static_cast<long double>(normal.x),
        static_cast<long double>(normal.y));
    if (!(length > 0) || !std::isfinite(length))
        throw std::invalid_argument("[PlanarCellCut] zero/invalid normal");
    return {(static_cast<long double>(cell.upper.x)-cell.lower.x)/2,
        (static_cast<long double>(cell.upper.y)-cell.lower.y)/2,
        static_cast<long double>(cell.lower.x)+(static_cast<long double>(cell.upper.x)-cell.lower.x)/2,
        normal.x/length, normal.y/length};
}
long double full(const Shape& s, bool rz)
{
    const long double value = 4*s.hx*s.hy*(rz ? CartesianMeasure::twoPi*s.radius : 1.L);
    if (!std::isfinite(value) || value <= 0)
        throw std::invalid_argument("[PlanarCellCut] unrepresentable capacity");
    return value;
}
long double cut(const Shape& s, long double offset, bool rz)
{
    if (!std::isfinite(offset)) throw std::invalid_argument("[PlanarCellCut] invalid offset");
    const long double extent = std::abs(s.nx)*s.hx+std::abs(s.ny)*s.hy;
    if (offset <= -extent) return 0;
    if (offset >= extent) return full(s,rz);
    const std::array<Point,4> corners{{{-s.hx,-s.hy},{s.hx,-s.hy},{s.hx,s.hy},{-s.hx,s.hy}}};
    std::array<Point,8> polygon{};
    int count=0;
    for (int k=0;k<4;++k) {
        const auto a=corners[k], b=corners[(k+1)%4];
        const auto da=s.nx*a.x+s.ny*a.y-offset, db=s.nx*b.x+s.ny*b.y-offset;
        if (da<=0) polygon[count++]=a;
        if ((da<0 && db>0) || (da>0 && db<0)) {
            const auto t=da/(da-db);
            polygon[count++]={a.x+t*(b.x-a.x),a.y+t*(b.y-a.y)};
        }
    }
    long double result=0;
    // Fan triangles avoid subtracting large global-coordinate moments.
    for (int k=1;k+1<count;++k) {
        const auto a=polygon[0],b=polygon[k],c=polygon[k+1];
        const auto area=((b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x))/2;
        result+=area*(rz ? CartesianMeasure::twoPi*(s.radius+(a.x+b.x+c.x)/3) : 1.L);
    }
    if (!std::isfinite(result) || result<0)
        throw std::runtime_error("[PlanarCellCut] invalid clipped measure");
    return result;
}
}
long double capacity(const Rectangle& cell, bool rz) { return full(shape(cell,{1,0},rz),rz); }
long double volume(const Rectangle& cell, Vector2D normal, long double offset, bool rz)
{ return cut(shape(cell,normal,rz),offset,rz); }
Reconstruction reconstruct(const Rectangle& cell, Vector2D normal, long double fraction, bool rz)
{
    if (!std::isfinite(fraction) || fraction<0 || fraction>1)
        throw std::invalid_argument("[PlanarCellCut] fraction outside [0,1]");
    const auto s=shape(cell,normal,rz);
    const auto total=full(s,rz), target=fraction*total;
    const auto extent=std::abs(s.nx)*s.hx+std::abs(s.ny)*s.hy;
    long double lo=-extent,hi=extent;
    long double offset=fraction==0?lo:(fraction==1?hi:0),actual=cut(s,offset,rz);
    if (fraction>0 && fraction<1) for (int iteration=0;iteration<96;++iteration) {
        offset=lo+(hi-lo)/2;
        actual=cut(s,offset,rz);
        if (actual==target || offset==lo || offset==hi) break;
        if (actual<target) lo=offset; else hi=offset;
    }
    if (std::abs(actual-target)>relativeVolumeResolution*total)
        throw std::runtime_error("[PlanarCellCut] unresolved PLIC volume constraint");
    return {{static_cast<double>(s.nx),static_cast<double>(s.ny)},offset,actual,actual-target};
}
CornerReconstruction reconstructThroughCorner(const Rectangle& cell,RectangleCorner corner,
    long double fraction,bool rz,bool negativeAlongHorizontalEdge)
{
    if (!std::isfinite(fraction) || fraction<0 || fraction>1)
        throw std::invalid_argument("[PlanarCellCut] corner fraction outside [0,1]");
    int sx=0,sy=0;
    switch (corner) {
    case RectangleCorner::LowerLeft: sx=-1;sy=-1;break;
    case RectangleCorner::LowerRight: sx=1;sy=-1;break;
    case RectangleCorner::UpperRight: sx=1;sy=1;break;
    case RectangleCorner::UpperLeft: sx=-1;sy=1;break;
    default: throw std::invalid_argument("[PlanarCellCut] invalid rectangle corner");
    }
    const auto rectangleShape=shape(cell,{1,0},rz);
    const auto total=full(rectangleShape,rz),target=fraction*total;
    const Point localCorner{sx*rectangleShape.hx,sy*rectangleShape.hy};
    const bool mixed=fraction>0 && fraction<1;
    if (mixed && (!(target>0) || !(target<total) || !((1-fraction)*total>0)))
        throw std::runtime_error("[PlanarCellCut] unrepresentable positive corner target");
    const auto candidate=[&](double ratio,bool horizontalRatio) {
        const double sign=negativeAlongHorizontalEdge?1.:-1.;
        // These two ratio charts are the same n(t)=(sx*(1-t),-sy*t)
        // family. They retain a small component without subtracting it from 1.
        const Vector2D raw{sign*sx*(horizontalRatio?ratio:1.),
            -sign*sy*(horizontalRatio?1.:ratio)};
        const auto length=std::hypot(static_cast<long double>(raw.x),static_cast<long double>(raw.y));
        const Vector2D normal{static_cast<double>(raw.x/length),static_cast<double>(raw.y/length)};
        const auto s=shape(cell,normal,rz); // Normalize the FINAL stored double.
        const auto offset=s.nx*localCorner.x+s.ny*localCorner.y;
        const auto actual=cut(s,offset,rz);
        return CornerReconstruction{{normal,offset,actual,actual-target},localCorner,corner};
    };
    const auto hasBothPhases=[&](const CornerReconstruction& value) {
        const auto& p=value.plane;
        return p.volume>0 && cut(shape(cell,{-p.normal.x,-p.normal.y},rz),-p.offset,rz)>0;
    };
    if (!mixed) {
        // t=0 is full for the true branch; t=1 is empty. Pure geometry
        // can coincide with an edge and does not establish a physical contact.
        const bool horizontalRatio=negativeAlongHorizontalEdge==(fraction==0);
        return candidate(0,horizontalRatio);
    }
    auto best=candidate(1,true);
    bool bestMixed=hasBothPhases(best);
    const auto consider=[&](const CornerReconstruction& value) {
        if ((!bestMixed || std::abs(value.plane.residual)<std::abs(best.plane.residual)) && hasBothPhases(value)) {
            best=value;bestMixed=true;
        }
    };
    const bool horizontalRatio=negativeAlongHorizontalEdge
        ?target<=best.plane.volume:target>best.plane.volume;
    const bool increasing=horizontalRatio==negativeAlongHorizontalEdge;
    const auto below=[&](const CornerReconstruction& value) {
        return increasing?value.plane.volume<target:value.plane.volume>target;
    };
    // First bracket the binary exponent. This covers the actual double normal
    // range without thousands of halvings or a hidden minimum component.
    int lowerExponent=std::numeric_limits<double>::min_exponent
        -std::numeric_limits<double>::digits-1,upperExponent=0;
    double lowerRatio=0,upperRatio=1;
    while (upperExponent-lowerExponent>1 && best.plane.residual!=0) {
        const int exponent=lowerExponent+(upperExponent-lowerExponent)/2;
        const double ratio=std::scalbn(1.,exponent);
        const auto value=candidate(ratio,horizontalRatio);consider(value);
        if (below(value)) {lowerExponent=exponent;lowerRatio=ratio;}
        else {upperExponent=exponent;upperRatio=ratio;}
    }
    // Within one exponent bin, bisection terminates at adjacent actual doubles.
    for (int iteration=0;iteration<64 && best.plane.residual!=0;++iteration) {
        const double ratio=lowerRatio+(upperRatio-lowerRatio)/2;
        if (ratio==lowerRatio || ratio==upperRatio) break;
        const auto value=candidate(ratio,horizontalRatio);consider(value);
        if (below(value)) lowerRatio=ratio;else upperRatio=ratio;
    }
    if (!bestMixed || std::abs(best.plane.residual)>relativeVolumeResolution*total)
        throw std::runtime_error("[PlanarCellCut] unresolved corner PLIC volume constraint");
    // The same final normal/offset are deliberately returned without another
    // double conversion or offset inverse, preserving the exact corner zero.
    return best;
}
static std::vector<Point> clipHalfPlane(const std::vector<Point>& source,long double nx,long double ny,long double d)
{
    std::vector<Point> result;result.reserve(source.size()+1);
    for (std::size_t i=0;i<source.size();++i) {
        const auto a=source[i],b=source[(i+1)%source.size()];
        const auto da=nx*a.x+ny*a.y-d,db=nx*b.x+ny*b.y-d;
        if (da<=0) result.push_back(a);
        if ((da<0 && db>0) || (da>0 && db<0)) {
            const auto t=da/(da-db);result.push_back({a.x+t*(b.x-a.x),a.y+t*(b.y-a.y)});
        }
    }
    return result;
}
static int polygonOrientation(const Point* input,std::size_t count)
{
    for (std::size_t i=0;i<count;++i) if (!std::isfinite(input[i].x) || !std::isfinite(input[i].y))
        throw std::invalid_argument("[PlanarCellCut] nonfinite polygon");
    if (count<3) return 0;
    int orientation=0;
    // The original all-vertices/edge proof excludes concavity and crossings.
    for (std::size_t i=0;i<count;++i) {
        const auto a=input[i],b=input[(i+1)%count];
        for (std::size_t j=0;j<count;++j) {
            const auto c=input[j];
            const int sign=Exact::sign(exactCross(a,b,c));
            if (sign && orientation && sign!=orientation)
                throw std::invalid_argument("[PlanarCellCut] polygon is not convex and cyclic");
            if (sign) orientation=sign;
        }
    }
    return orientation;
}
static std::vector<Point> clipRectangle(const std::vector<Point>& input,Point lower,Point upper)
{
    if(input.size()==4){
        const auto a=input[0],b=input[1],c=input[2],d=input[3];
        // Canonical LL/LR/UR/UL is already a convex CCW rectangle. If wholly
        // inside, the general four clips would preserve these exact vertices.
        if(a.x==d.x&&b.x==c.x&&a.y==b.y&&c.y==d.y
            &&std::isfinite(a.x)&&std::isfinite(a.y)&&std::isfinite(b.x)&&std::isfinite(c.y)
            &&a.x<b.x&&a.y<c.y&&a.x>=lower.x&&a.y>=lower.y&&c.x<=upper.x&&c.y<=upper.y)
            return input;
    }
    const int orientation=polygonOrientation(input.data(),input.size());
    if (!orientation) return {};
    auto polygon=input;
    if (orientation<0) std::reverse(polygon.begin(),polygon.end());
    polygon=clipHalfPlane(polygon,1,0,upper.x);polygon=clipHalfPlane(polygon,-1,0,-lower.x);
    polygon=clipHalfPlane(polygon,0,1,upper.y);polygon=clipHalfPlane(polygon,0,-1,-lower.y);
    return polygon;
}
std::vector<Point> clipPolygonInFrame(const Rectangle& cell,const std::vector<Point>& polygon,Point origin)
{
    shape(cell,{1,0},false);
    if (!std::isfinite(origin.x) || !std::isfinite(origin.y))
        throw std::invalid_argument("[PlanarCellCut] invalid polygon frame origin");
    return clipRectangle(polygon,{static_cast<long double>(cell.lower.x)-origin.x,static_cast<long double>(cell.lower.y)-origin.y},
        {static_cast<long double>(cell.upper.x)-origin.x,static_cast<long double>(cell.upper.y)-origin.y});
}
static Overlap overlapAt(const Shape& s,const std::vector<Point>& input,Point center,
    Point lower,Point upper,long double radius,long double offset,bool rz)
{
    if (!std::isfinite(offset)) throw std::invalid_argument("[PlanarCellCut] invalid plane offset");
    const auto polygon=clipRectangle(input,lower,upper);
    const auto measure=[&](const std::vector<Point>& p) {
        long double value=0;
        for (std::size_t i=1;i+1<p.size();++i) {
            const auto a=p[0],b=p[i],c=p[i+1];
            const auto area=((b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x))/2;
            value+=area*(rz?CartesianMeasure::twoPi*(radius+(a.x+b.x+c.x)/3):1.L);
        }
        if (!std::isfinite(value) || value<0)
            throw std::runtime_error("[PlanarCellCut] invalid swept moment");
        return value;
    };
    const auto total=measure(polygon),liquid=measure(clipHalfPlane(polygon,s.nx,s.ny,offset+s.nx*center.x+s.ny*center.y));
    return {total,liquid};
}
Overlap overlap(const Rectangle& cell,const std::vector<Point>& input,Vector2D normal,long double offset,bool rz)
{
    const auto s=shape(cell,normal,rz);
    return overlapAt(s,input,{0,0},{-s.hx,-s.hy},{s.hx,s.hy},s.radius,offset,rz);
}
Overlap overlapInFrame(const Rectangle& cell,const std::vector<Point>& input,Point origin,
    Vector2D normal,long double offset,bool rz)
{
    if (!std::isfinite(origin.x) || !std::isfinite(origin.y))
        throw std::invalid_argument("[PlanarCellCut] invalid polygon frame origin");
    const auto s=shape(cell,normal,rz);
    const Point center{s.radius-origin.x,
        static_cast<long double>(cell.lower.y)+(static_cast<long double>(cell.upper.y)-cell.lower.y)/2-origin.y};
    const Point lower{static_cast<long double>(cell.lower.x)-origin.x,static_cast<long double>(cell.lower.y)-origin.y};
    const Point upper{static_cast<long double>(cell.upper.x)-origin.x,static_cast<long double>(cell.upper.y)-origin.y};
    return overlapAt(s,input,center,lower,upper,origin.x,offset,rz);
}
std::optional<std::array<Point,2>> interfaceSegment(const Rectangle& cell,Vector2D normal,long double offset)
{
    const auto s=shape(cell,normal,false);
    if (!std::isfinite(offset)) throw std::invalid_argument("[PlanarCellCut] invalid segment offset");
    const std::array<Point,4> corners{{{-s.hx,-s.hy},{s.hx,-s.hy},{s.hx,s.hy},{-s.hx,s.hy}}};
    std::vector<Point> points;
    const auto append=[&](Point p) {
        for (auto q:points) if (q.x==p.x && q.y==p.y) return;
        points.push_back(p);
    };
    for (int i=0;i<4;++i) {
        const auto a=corners[i],b=corners[(i+1)%4];
        const auto da=s.nx*a.x+s.ny*a.y-offset,db=s.nx*b.x+s.ny*b.y-offset;
        if (da==0) append(a);
        if ((da<0 && db>0) || (da>0 && db<0)) {
            const auto t=da/(da-db);append({a.x+t*(b.x-a.x),a.y+t*(b.y-a.y)});
        }
    }
    if (points.size()>2) throw std::runtime_error("[PlanarCellCut] inconsistent plane/rectangle intersections");
    if (points.size()!=2 || (points[0].x==points[1].x && points[0].y==points[1].y)) return std::nullopt;
    return std::array<Point,2>{{points[0],points[1]}};
}

namespace
{
bool equalPoint(Point a,Point b) { return a.x==b.x && a.y==b.y; }
Point crossing(Point a,Point b,long double va,long double vb)
{
    // Use the closer endpoint and a ratio <=1: no fa-fb overflow, no loss of
    // a resolvable tiny phase through interpolation from the distant endpoint.
    const auto aa=std::abs(va),bb=std::abs(vb);
    if (aa>bb || (aa==bb && (a.x>b.x || (a.x==b.x && a.y>b.y)))) {
        std::swap(a,b);std::swap(va,vb);
    }
    const auto ratio=std::abs(va)/std::abs(vb);
    const auto t=ratio/(1+ratio);
    const Point p{std::fma(t,b.x-a.x,a.x),std::fma(t,b.y-a.y,a.y)};
    if (!std::isfinite(p.x) || !std::isfinite(p.y) || t<=0 || t>=1
        || equalPoint(p,a) || equalPoint(p,b))
        throw std::runtime_error("[PlanarCellCut] line/edge intersection is unrepresentable");
    return p;
}
void validateOrigin(Point origin)
{
    if (!std::isfinite(origin.x) || !std::isfinite(origin.y))
        throw std::invalid_argument("[PlanarCellCut] invalid polygon frame origin");
}
long double measurePoints(const Point* points,std::size_t count,Point origin,bool rz)
{
    validateOrigin(origin);
    const int orientation=polygonOrientation(points,count);
    const auto radialOrigin=rz?exactNumber(origin.x):Expanded{};
    if (rz) for (std::size_t i=0;i<count;++i)
        if (Exact::sign(Exact::add(radialOrigin,exactNumber(points[i].x)))<0)
            throw std::invalid_argument("[PlanarCellCut] negative/nonfinite polygon radius");
    if (!orientation) return 0;
    Expanded numerator;
    for (std::size_t i=1;i+1<count;++i) {
        const auto a=points[0],b=points[i],c=points[i+1];
        auto term=positiveTwiceArea(a,b,c,orientation);
        if (rz) {
            // Integrate the physical radius of these actual Point vertices.
            // Keeping the factor six outside avoids rounded area/centroid
            // intermediates and retains a thin polygon's represented volume.
            const auto radiusSum=Exact::add(Exact::multiply(radialOrigin,Exact::fromDouble(3)),
                coordinateSum(a.x,b.x,c.x));
            if (Exact::sign(radiusSum)<0 || (Exact::sign(term)>0 && Exact::sign(radiusSum)==0))
                throw std::runtime_error("[PlanarCellCut] unrepresentable polygon triangle radius");
            term=Exact::multiply(term,radiusSum);
        }
        Exact::accumulate(numerator,term);
    }
    if (rz) numerator=Exact::multiply(numerator,exactNumber(CartesianMeasure::twoPi));
    const auto value=dividedValue(numerator,rz?6.L:2.L);
    if (value<0)
        throw std::runtime_error("[PlanarCellCut] invalid physical polygon measure");
    return value;
}
}
Point linearZeroCrossing(Point a,Point b,long double va,long double vb)
{
    if (!std::isfinite(a.x) || !std::isfinite(a.y) || !std::isfinite(b.x) || !std::isfinite(b.y)
        || !std::isfinite(va) || !std::isfinite(vb) || equalPoint(a,b)
        || !((va<0 && vb>0) || (va>0 && vb<0)))
        throw std::invalid_argument("[PlanarCellCut] crossing requires distinct finite endpoints and opposite signs");
    return crossing(a,b,va,vb);
}
PolygonMoments polygonMoments(const std::vector<Point>& points)
{
    const int orientation=polygonOrientation(points.data(),points.size());
    PolygonMoments result;
    if (!orientation) return result;
    Expanded twiceArea,sixFirstX,sixFirstY;
    for (std::size_t i=1;i+1<points.size();++i) {
        const auto a=points[0],b=points[i],c=points[i+1];
        const auto twice=positiveTwiceArea(a,b,c,orientation);
        Exact::accumulate(twiceArea,twice);
        Exact::accumulate(sixFirstX,Exact::multiply(twice,coordinateSum(a.x,b.x,c.x)));
        Exact::accumulate(sixFirstY,Exact::multiply(twice,coordinateSum(a.y,b.y,c.y)));
    }
    result.area=dividedValue(twiceArea,2);
    result.firstX=dividedValue(sixFirstX,6);
    result.firstY=dividedValue(sixFirstY,6);
    return result;
}
long double polygonMeasureInFrame(const std::vector<Point>& polygon,Point origin,bool rz)
{ return measurePoints(polygon.data(),polygon.size(),origin,rz); }
std::vector<Point> clipConvexPolygon(const std::vector<Point>& input,Point normal,long double offset)
{
    if (!std::isfinite(normal.x) || !std::isfinite(normal.y) || (normal.x==0 && normal.y==0)
        || !std::isfinite(offset))
        throw std::invalid_argument("[PlanarCellCut] invalid polygon cutting plane");
    const int orientation=polygonOrientation(input.data(),input.size());
    if (!orientation) return {};
    std::vector<Point> result;result.reserve(input.size()+1);
    for (std::size_t i=0;i<input.size();++i) {
        const auto a=input[i],b=input[(i+1)%input.size()];
        const auto da=normal.x*a.x+normal.y*a.y-offset,db=normal.x*b.x+normal.y*b.y-offset;
        if (!std::isfinite(da) || !std::isfinite(db))
            throw std::invalid_argument("[PlanarCellCut] unrepresentable polygon plane value");
        if (da<=0) result.push_back(a);
        if ((da<0 && db>0) || (da>0 && db<0)) result.push_back(crossing(a,b,da,db));
    }
    if (orientation<0) std::reverse(result.begin(),result.end());
    return result;
}
static LinearTriangleCut triangleCutData(const std::array<Point,3>& triangle,const std::array<long double,3>& scalar,bool withInterface)
{
    const int orientation=polygonOrientation(triangle.data(),triangle.size());
    if (!orientation) throw std::invalid_argument("[PlanarCellCut] degenerate triangle");
    LinearTriangleCut result;
    unsigned negatives=0,positives=0;
    for (unsigned i=0;i<3;++i) {
        if (!std::isfinite(scalar[i])) throw std::invalid_argument("[PlanarCellCut] nonfinite triangle scalar");
        if (scalar[i]<0) ++negatives;
        if (scalar[i]>0) ++positives;
        if (scalar[i]==0) result.zeroVertexMask|=1u<<i;
    }
    for (unsigned i=0;i<3;++i) if (scalar[i]==0 && scalar[(i+1)%3]==0) result.zeroEdgeMask|=1u<<i;
    if (!negatives && !positives) return result; // Positive area, undefined phase.
    result.kind=negatives?(positives?TriangleRegion::Cut:TriangleRegion::Negative):TriangleRegion::Positive;
    std::array<Point,2> zeros{};std::size_t zeroCount=0;
    const auto appendZero=[&](Point p) {
        for (std::size_t j=0;j<zeroCount;++j) if (equalPoint(zeros[j],p)) return;
        if (zeroCount==zeros.size()) throw std::runtime_error("[PlanarCellCut] triangle zero-line topology is inconsistent");
        zeros[zeroCount++]=p;
    };
    for (unsigned i=0;i<3;++i) {
        const auto a=triangle[i],b=triangle[(i+1)%3];
        const auto va=scalar[i],vb=scalar[(i+1)%3];
        if (va<=0) result.negativeVertices[result.negativeCount++]=a;
        if (withInterface && va==0) appendZero(a);
        if ((va<0 && vb>0) || (va>0 && vb<0)) {
            const auto p=crossing(a,b,va,vb);
            result.negativeVertices[result.negativeCount++]=p;if (withInterface) appendZero(p);
        }
    }
    if (zeroCount==2) {
        const auto scale=std::max({std::abs(scalar[0]),std::abs(scalar[1]),std::abs(scalar[2])});
        int exponent=0;std::frexp(scale,&exponent);
        const auto scaled=[&](long double value) {return Exact::scalePow2(exactNumber(value),-exponent);};
        const auto s0=scaled(scalar[0]);
        const auto d1=Exact::subtract(scaled(scalar[1]),s0),d2=Exact::subtract(scaled(scalar[2]),s0);
        const auto x0=exactNumber(triangle[0].x),y0=exactNumber(triangle[0].y);
        const auto dx1=Exact::subtract(exactNumber(triangle[1].x),x0),dy1=Exact::subtract(exactNumber(triangle[1].y),y0);
        const auto dx2=Exact::subtract(exactNumber(triangle[2].x),x0),dy2=Exact::subtract(exactNumber(triangle[2].y),y0);
        auto gx=Exact::subtract(Exact::multiply(d1,dy2),Exact::multiply(d2,dy1));
        auto gy=Exact::subtract(Exact::multiply(dx1,d2),Exact::multiply(dx2,d1));
        if (orientation<0) {gx=Exact::subtract(Expanded{},gx);gy=Exact::subtract(Expanded{},gy);}
        const auto nx=dividedValue(gx,1),ny=dividedValue(gy,1);
        const auto length=std::hypot(nx,ny);
        if (!std::isfinite(length) || length<=0)
            throw std::runtime_error("[PlanarCellCut] triangle interface normal is unrepresentable");
        result.unitNormal={nx/length,ny/length};
        const auto zx=Exact::subtract(exactNumber(zeros[1].x),exactNumber(zeros[0].x));
        const auto zy=Exact::subtract(exactNumber(zeros[1].y),exactNumber(zeros[0].y));
        const int tangent=Exact::sign(Exact::subtract(Exact::multiply(gx,zy),Exact::multiply(gy,zx)));
        if (!tangent)
            throw std::runtime_error("[PlanarCellCut] triangle interface direction is unrepresentable");
        if (tangent<0) std::swap(zeros[0],zeros[1]);
        result.interface=zeros;
    }
    return result;
}
LinearTriangleCut cutTriangle(const std::array<Point,3>& triangle,const std::array<long double,3>& scalar)
{ return triangleCutData(triangle,scalar,true); }
long double negativeMeasure(const LinearTriangleCut& result,Point origin,bool rz)
{
    if (result.kind==TriangleRegion::AllZeroScalar)
        throw std::domain_error("[PlanarCellCut] positive-area all-zero scalar has no phase measure");
    if (result.negativeCount>result.negativeVertices.size())
        throw std::invalid_argument("[PlanarCellCut] invalid fixed triangle polygon count");
    const auto value=measurePoints(result.negativeVertices.data(),result.negativeCount,origin,rz);
    if ((result.kind==TriangleRegion::Cut || result.kind==TriangleRegion::Negative) && value<=0)
        throw std::runtime_error("[PlanarCellCut] positive-area negative phase measure is unrepresentable");
    return value;
}

long double negativeMeasure(const std::array<Point,3>& triangle,const std::array<long double,3>& scalar,
    Point origin,bool rz)
{
    validateOrigin(origin);
    if (rz) for (const auto p:triangle) {
        const auto r=origin.x+p.x;
        if (!std::isfinite(r) || r<0)
            throw std::invalid_argument("[PlanarCellCut] negative/nonfinite source triangle radius");
    }
    return negativeMeasure(triangleCutData(triangle,scalar,false),origin,rz);
}

// Explicit anchored-site operators. The original Point operators above are
// retained unchanged; all predicates below read the same public actual sites.
namespace
{
namespace AnchoredRound=LAMP::Geometry::Detail::ExactRationalRound;
struct AnchoredExactPoint {Expanded x,y;};
void anchoredRequire(bool ok,const char* message)
{if(!ok)throw std::invalid_argument(std::string("[PlanarCellCut/Anchored] ")+message);}
Expanded anchoredNegate(const Expanded& value){return Exact::subtract(Expanded{},value);}
int anchoredCompare(const Expanded& a,const Expanded& b){return Exact::sign(Exact::subtract(a,b));}
AnchoredExactPoint anchoredExact(AnchoredPoint p)
{
    return {Exact::add(AnchoredRound::scalar(p.anchor.x),AnchoredRound::scalar(p.offset.x)),
        Exact::add(AnchoredRound::scalar(p.anchor.y),AnchoredRound::scalar(p.offset.y))};
}
AnchoredExactPoint anchoredSubtract(const AnchoredExactPoint& a,const AnchoredExactPoint& b)
{return {Exact::subtract(a.x,b.x),Exact::subtract(a.y,b.y)};}
Expanded anchoredCrossExact(AnchoredPoint a,AnchoredPoint b,AnchoredPoint c)
{
    const auto origin=anchoredExact(a),u=anchoredSubtract(anchoredExact(b),origin),v=anchoredSubtract(anchoredExact(c),origin);
    return Exact::subtract(Exact::multiply(u.x,v.y),Exact::multiply(u.y,v.x));
}
bool anchoredSame(AnchoredPoint a,AnchoredPoint b)
{
    const auto x=anchoredExact(a),y=anchoredExact(b);
    return anchoredCompare(x.x,y.x)==0 && anchoredCompare(x.y,y.y)==0;
}
bool anchoredLexGreater(AnchoredPoint a,AnchoredPoint b)
{
    const auto x=anchoredExact(a),y=anchoredExact(b);const int first=anchoredCompare(x.x,y.x);
    return first>0 || (first==0 && anchoredCompare(x.y,y.y)>0);
}
// Each first factor and every unmultiplied value receive ONE power of two.
// Include all retained least bits BEFORE products; reserve three carry bits
// for the four signed numerator terms and the positive denominator sum.
int anchoredProductShift(const std::vector<std::pair<Expanded,Expanded>>& products,
    const std::vector<Expanded>& unmultiplied)
{
    int lowest=std::numeric_limits<int>::max(),highest=std::numeric_limits<int>::min();
    for(const auto& value:unmultiplied){Exact::ExponentRange range;range.include(value);
        if(range.lowest!=std::numeric_limits<int>::max()){lowest=std::min(lowest,range.lowest);highest=std::max(highest,range.highest);}}
    for(const auto& product:products){Exact::ExponentRange a,b;a.include(product.first);b.include(product.second);
        if(a.lowest==std::numeric_limits<int>::max() || b.lowest==std::numeric_limits<int>::max())continue;
        lowest=std::min(lowest,a.lowest+b.lowest);highest=std::max(highest,a.highest+b.highest+1);}
    if(lowest==std::numeric_limits<int>::max())return 0;
    const int lower=-1074-lowest,upper=1020-highest;
    if(lower>upper)throw std::runtime_error("[PlanarCellCut/Anchored] exact weighted products exceed arithmetic window");
    return std::clamp(-(lowest+highest)/2,lower,upper);
}
long double anchoredOffsetCoordinate(long double aOffset,long double aAnchor,
    long double bOffset,long double bAnchor,const Expanded& weightA,const Expanded& weightB)
{
    const auto ao=AnchoredRound::scalar(aOffset),aa=AnchoredRound::scalar(aAnchor);
    const auto bo=AnchoredRound::scalar(bOffset),ba=AnchoredRound::scalar(bAnchor);
    const int shift=anchoredProductShift({{weightB,ao},{weightA,bo},{weightA,ba},{weightA,aa}},{weightA,weightB});
    const auto wa=Exact::scalePow2(weightA,shift),wb=Exact::scalePow2(weightB,shift);
    // Distribute the anchor difference BEFORE products, avoiding an otherwise
    // unnecessary overflow of b.anchor-a.anchor. The expression is unchanged.
    auto numerator=Exact::add(Exact::multiply(wb,ao),Exact::multiply(wa,bo));
    numerator=Exact::add(numerator,Exact::multiply(wa,ba));
    numerator=Exact::subtract(numerator,Exact::multiply(wa,aa));
    return AnchoredRound::quotient(numerator,Exact::add(wa,wb));
}
AnchoredPoint anchoredCrossing(AnchoredPoint a,AnchoredPoint b,Expanded valueA,Expanded valueB)
{
    const int signA=Exact::sign(valueA),signB=Exact::sign(valueB);
    anchoredRequire(signA*signB<0 && !anchoredSame(a,b),"crossing requires distinct sites and strictly opposite scalars");
    auto wa=signA<0?anchoredNegate(valueA):valueA,wb=signB<0?anchoredNegate(valueB):valueB;
    const int nearer=anchoredCompare(wa,wb);
    if(nearer>0 || (nearer==0 && anchoredLexGreater(a,b))){std::swap(a,b);std::swap(wa,wb);}
    AnchoredPoint result;result.anchor=a.anchor;
    result.offset={anchoredOffsetCoordinate(a.offset.x,a.anchor.x,b.offset.x,b.anchor.x,wa,wb),
        anchoredOffsetCoordinate(a.offset.y,a.anchor.y,b.offset.y,b.anchor.y,wa,wb)};
    if(anchoredSame(result,a) || anchoredSame(result,b))
        throw std::runtime_error("[PlanarCellCut/Anchored] actual line/edge intersection collapsed onto an endpoint");
    return result;
}
int anchoredOrientation(const AnchoredPoint* input,std::size_t count)
{
    for(std::size_t i=0;i<count;++i)anchoredExact(input[i]);
    if(count<3)return 0;
    int orientation=0;
    for(std::size_t i=0;i<count;++i)for(std::size_t j=0;j<count;++j){
        const int sign=Exact::sign(anchoredCrossExact(input[i],input[(i+1)%count],input[j]));
        anchoredRequire(!sign || !orientation || sign==orientation,"polygon is not convex and cyclic");
        if(sign)orientation=sign;
    }
    return orientation;
}
Expanded anchoredPositiveArea(AnchoredPoint a,AnchoredPoint b,AnchoredPoint c,int orientation)
{
    auto twice=anchoredCrossExact(a,b,c);if(orientation<0)twice=anchoredNegate(twice);
    anchoredRequire(Exact::sign(twice)>=0,"polygon fan has a negative actual triangle");return twice;
}
long double anchoredDivide(const Expanded& numerator,long double denominator)
{return AnchoredRound::quotient(numerator,AnchoredRound::scalar(denominator));}
long double anchoredMeasure(const AnchoredPoint* points,std::size_t count,Point origin,bool rz)
{
    validateOrigin(origin);const int orientation=anchoredOrientation(points,count);
    const auto radialOrigin=rz?AnchoredRound::scalar(origin.x):Expanded{};
    if(rz)for(std::size_t i=0;i<count;++i)
        anchoredRequire(Exact::sign(Exact::add(radialOrigin,anchoredExact(points[i]).x))>=0,"negative physical polygon radius");
    if(!orientation)return 0;
    Expanded numerator;
    for(std::size_t i=1;i+1<count;++i){
        auto term=anchoredPositiveArea(points[0],points[i],points[i+1],orientation);
        if(rz){
            // Form physical radii from the actual site before multiplying;
            // no translated native Point or rounded centroid is constructed.
            const auto ra=Exact::add(radialOrigin,anchoredExact(points[0]).x);
            const auto rb=Exact::add(radialOrigin,anchoredExact(points[i]).x);
            const auto rc=Exact::add(radialOrigin,anchoredExact(points[i+1]).x);
            const auto radii=Exact::add(Exact::add(ra,rb),rc);
            if(Exact::sign(term)>0 && Exact::sign(radii)==0)
                throw std::runtime_error("[PlanarCellCut/Anchored] positive-area physical radius is unrepresentable");
            term=Exact::multiply(term,radii);
        }
        Exact::accumulate(numerator,term);
    }
    if(rz)numerator=Exact::multiply(numerator,AnchoredRound::scalar(CartesianMeasure::twoPi));
    const auto result=anchoredDivide(numerator,rz?6.L:2.L);
    anchoredRequire(result>=0,"negative actual physical polygon measure");return result;
}
void anchoredValidateBox(const Rectangle& box)
{
    anchoredRequire(std::isfinite(box.lower.x) && std::isfinite(box.lower.y) && std::isfinite(box.upper.x) &&
        std::isfinite(box.upper.y) && box.upper.x>box.lower.x && box.upper.y>box.lower.y,"invalid rectangle");
}
std::vector<AnchoredPoint> anchoredClipValues(const std::vector<AnchoredPoint>& polygon,
    const std::vector<Expanded>& values)
{
    std::vector<AnchoredPoint> result;result.reserve(polygon.size()+1);
    for(std::size_t i=0;i<polygon.size();++i){const auto j=(i+1)%polygon.size();
        const int a=Exact::sign(values[i]),b=Exact::sign(values[j]);
        if(a<=0)result.push_back(polygon[i]);
        if(a*b<0)result.push_back(anchoredCrossing(polygon[i],polygon[j],values[i],values[j]));
    }
    return result;
}
AnchoredTriangleCut anchoredTriangleData(const std::array<AnchoredPoint,3>& triangle,
    const std::array<long double,3>& scalar,bool withInterface)
{
    const int orientation=anchoredOrientation(triangle.data(),triangle.size());
    anchoredRequire(orientation!=0,"degenerate triangle");AnchoredTriangleCut result;
    unsigned negatives=0,positives=0;
    for(unsigned i=0;i<3;++i){anchoredRequire(std::isfinite(scalar[i]),"nonfinite triangle scalar");
        negatives+=scalar[i]<0;positives+=scalar[i]>0;if(scalar[i]==0)result.zeroVertexMask|=1u<<i;}
    for(unsigned i=0;i<3;++i)if(scalar[i]==0 && scalar[(i+1)%3]==0)result.zeroEdgeMask|=1u<<i;
    if(!negatives && !positives)return result;
    result.kind=negatives?(positives?TriangleRegion::Cut:TriangleRegion::Negative):TriangleRegion::Positive;
    std::array<AnchoredPoint,2> zeros{};std::size_t zeroCount=0;
    const auto appendZero=[&](AnchoredPoint p){
        for(std::size_t j=0;j<zeroCount;++j)if(anchoredSame(zeros[j],p))return;
        if(zeroCount==zeros.size())throw std::runtime_error("[PlanarCellCut/Anchored] inconsistent triangle zero topology");
        zeros[zeroCount++]=p;
    };
    for(unsigned i=0;i<3;++i){const unsigned j=(i+1)%3;
        if(scalar[i]<=0)result.negativeVertices[result.negativeCount++]=triangle[i];
        if(withInterface && scalar[i]==0)appendZero(triangle[i]);
        if((scalar[i]<0 && scalar[j]>0) || (scalar[i]>0 && scalar[j]<0)){
            const auto p=anchoredCrossing(triangle[i],triangle[j],AnchoredRound::scalar(scalar[i]),AnchoredRound::scalar(scalar[j]));
            result.negativeVertices[result.negativeCount++]=p;if(withInterface)appendZero(p);}
    }
    if(zeroCount==2){
        std::array<Expanded,3> values{{AnchoredRound::scalar(scalar[0]),AnchoredRound::scalar(scalar[1]),AnchoredRound::scalar(scalar[2])}};
        Exact::ExponentRange valueRange;for(const auto& v:values)valueRange.include(v);
        const int valueShift=valueRange.centeredShift();for(auto& v:values)v=Exact::scalePow2(v,valueShift);
        const auto ds1=Exact::subtract(values[1],values[0]),ds2=Exact::subtract(values[2],values[0]);
        const auto p0=anchoredExact(triangle[0]),d1=anchoredSubtract(anchoredExact(triangle[1]),p0),d2=anchoredSubtract(anchoredExact(triangle[2]),p0);
        auto gx=Exact::subtract(Exact::multiply(ds1,d2.y),Exact::multiply(ds2,d1.y));
        auto gy=Exact::subtract(Exact::multiply(ds2,d1.x),Exact::multiply(ds1,d2.x));
        if(orientation<0){gx=anchoredNegate(gx);gy=anchoredNegate(gy);}
        Exact::ExponentRange gradientRange;gradientRange.include(gx);gradientRange.include(gy);
        const int gradientShift=gradientRange.centeredShift();gx=Exact::scalePow2(gx,gradientShift);gy=Exact::scalePow2(gy,gradientShift);
        auto tangent=anchoredSubtract(anchoredExact(zeros[1]),anchoredExact(zeros[0]));
        const int direction=Exact::sign(Exact::subtract(Exact::multiply(gx,tangent.y),Exact::multiply(gy,tangent.x)));
        if(!direction)throw std::runtime_error("[PlanarCellCut/Anchored] unrepresentable interface direction");
        if(direction<0){std::swap(zeros[0],zeros[1]);tangent.x=anchoredNegate(tangent.x);tangent.y=anchoredNegate(tangent.y);}
        // The cofactor chooses the negative-left orientation. The metric
        // normal is perpendicular to this SAME public actual-site chord.
        Exact::ExponentRange chordRange;chordRange.include(tangent.x);chordRange.include(tangent.y);
        const int chordShift=chordRange.centeredShift();
        const auto dx=anchoredDivide(Exact::scalePow2(tangent.x,chordShift),1);
        const auto dy=anchoredDivide(Exact::scalePow2(tangent.y,chordShift),1);
        const auto length=std::hypot(dx,dy);
        if(!std::isfinite(length) || length<=0)throw std::runtime_error("[PlanarCellCut/Anchored] unrepresentable interface normal");
        result.unitNormal={dy/length,-dx/length};result.interface=zeros;
    }
    return result;
}
}
Point approximateAnchoredPoint(AnchoredPoint point)
{
    Exact::validateEnvironment();const auto value=anchoredExact(point);
    return {anchoredDivide(value.x,1),anchoredDivide(value.y,1)};
}
Point anchoredDifference(AnchoredPoint a,AnchoredPoint b)
{
    Exact::validateEnvironment();const auto difference=anchoredSubtract(anchoredExact(a),anchoredExact(b));
    return {anchoredDivide(difference.x,1),anchoredDivide(difference.y,1)};
}
Vector2D anchoredDifferenceDouble(AnchoredPoint query,AnchoredPoint frameOrigin,Point localCenter)
{
    Exact::validateActualEnvironment();
    anchoredRequire(std::isfinite(query.anchor.x)&&std::isfinite(query.anchor.y)
        &&std::isfinite(query.offset.x)&&std::isfinite(query.offset.y)
        &&std::isfinite(frameOrigin.anchor.x)&&std::isfinite(frameOrigin.anchor.y)
        &&std::isfinite(frameOrigin.offset.x)&&std::isfinite(frameOrigin.offset.y)
        &&std::isfinite(localCenter.x)&&std::isfinite(localCenter.y),"nonfinite binary64 frame input");
    const auto coordinate=[](long double qa,long double qo,long double fa,long double fo,long double center){
        if constexpr(std::numeric_limits<long double>::radix==2
            &&std::numeric_limits<long double>::digits==std::numeric_limits<double>::digits
            &&std::numeric_limits<long double>::min_exponent==std::numeric_limits<double>::min_exponent
            &&std::numeric_limits<long double>::max_exponent==std::numeric_limits<double>::max_exponent){
            if(qa==fa&&fo==0){
                // Equal anchors cancel exactly. Only one binary64 subtraction
                // remains; no three-term partially rounded frame is formed.
                const double value=static_cast<double>(qo)-static_cast<double>(center);
                anchoredRequire(std::isfinite(value),"final binary64 frame overflow");
                return value==0?0.:value;
            }
        }
        Detail::FiniteVolume::Dyadic difference(qa);difference.add(qo);
        difference.add(fa,-1);difference.add(fo,-1);difference.add(center,-1);
        return difference.roundedDouble();
    };
    return {coordinate(query.anchor.x,query.offset.x,frameOrigin.anchor.x,frameOrigin.offset.x,localCenter.x),
        coordinate(query.anchor.y,query.offset.y,frameOrigin.anchor.y,frameOrigin.offset.y,localCenter.y)};
}
bool sameAnchoredPoint(AnchoredPoint a,AnchoredPoint b)
{Exact::validateEnvironment();return anchoredSame(a,b);}
AnchoredPoint anchoredZeroCrossing(AnchoredPoint a,AnchoredPoint b,long double valueA,long double valueB)
{
    Exact::validateEnvironment();anchoredRequire(std::isfinite(valueA) && std::isfinite(valueB),"nonfinite crossing scalar");
    return anchoredCrossing(a,b,AnchoredRound::scalar(valueA),AnchoredRound::scalar(valueB));
}
AnchoredTriangleCut cutAnchoredTriangle(const std::array<AnchoredPoint,3>& triangle,const std::array<long double,3>& scalar)
{Exact::validateEnvironment();return anchoredTriangleData(triangle,scalar,true);}
long double negativeAnchoredMeasure(const AnchoredTriangleCut& cut,Point origin,bool rz)
{
    Exact::validateEnvironment();if(cut.kind==TriangleRegion::AllZeroScalar)
        throw std::domain_error("[PlanarCellCut/Anchored] positive-area all-zero scalar has no phase measure");
    anchoredRequire(cut.negativeCount<=cut.negativeVertices.size(),"invalid fixed triangle polygon count");
    const auto value=anchoredMeasure(cut.negativeVertices.data(),cut.negativeCount,origin,rz);
    if((cut.kind==TriangleRegion::Cut || cut.kind==TriangleRegion::Negative) && value<=0)
        throw std::runtime_error("[PlanarCellCut/Anchored] positive-area negative phase is unrepresentable");
    return value;
}
long double negativeAnchoredMeasure(const std::array<AnchoredPoint,3>& triangle,const std::array<long double,3>& scalar,Point origin,bool rz)
{
    Exact::validateEnvironment();validateOrigin(origin);
    if(rz){const auto radialOrigin=AnchoredRound::scalar(origin.x);
        for(const auto& point:triangle)anchoredRequire(Exact::sign(Exact::add(radialOrigin,anchoredExact(point).x))>=0,
            "negative source triangle physical radius");}
    return negativeAnchoredMeasure(anchoredTriangleData(triangle,scalar,false),origin,rz);
}
PolygonMoments anchoredPolygonMoments(const std::vector<AnchoredPoint>& polygon)
{
    Exact::validateEnvironment();const int orientation=anchoredOrientation(polygon.data(),polygon.size());
    PolygonMoments result;if(!orientation)return result;Expanded twiceArea,sixX,sixY;
    for(std::size_t i=1;i+1<polygon.size();++i){
        const auto term=anchoredPositiveArea(polygon[0],polygon[i],polygon[i+1],orientation);
        const auto a=anchoredExact(polygon[0]),b=anchoredExact(polygon[i]),c=anchoredExact(polygon[i+1]);
        Exact::accumulate(twiceArea,term);
        Exact::accumulate(sixX,Exact::multiply(term,Exact::add(Exact::add(a.x,b.x),c.x)));
        Exact::accumulate(sixY,Exact::multiply(term,Exact::add(Exact::add(a.y,b.y),c.y)));
    }
    result.area=anchoredDivide(twiceArea,2);result.firstX=anchoredDivide(sixX,6);result.firstY=anchoredDivide(sixY,6);return result;
}
long double anchoredPolygonMeasureInFrame(const std::vector<AnchoredPoint>& polygon,Point origin,bool rz)
{Exact::validateEnvironment();return anchoredMeasure(polygon.data(),polygon.size(),origin,rz);}
std::vector<AnchoredPoint> clipAnchoredPolygonInFrame(const Rectangle& box,const std::vector<AnchoredPoint>& input,Point origin)
{
    Exact::validateEnvironment();anchoredValidateBox(box);validateOrigin(origin);
    const int orientation=anchoredOrientation(input.data(),input.size());if(!orientation)return {};
    auto polygon=input;if(orientation<0)std::reverse(polygon.begin(),polygon.end());
    const auto ox=AnchoredRound::scalar(origin.x),oy=AnchoredRound::scalar(origin.y);
    const std::array<Expanded,4> bounds{{AnchoredRound::scalar(box.upper.x),AnchoredRound::scalar(box.lower.x),
        AnchoredRound::scalar(box.upper.y),AnchoredRound::scalar(box.lower.y)}};
    for(unsigned side=0;side<4 && !polygon.empty();++side){std::vector<Expanded> values;values.reserve(polygon.size());
        for(const auto& p:polygon){const auto actual=anchoredExact(p);
            const auto world=Exact::add(side<2?actual.x:actual.y,side<2?ox:oy);
            values.push_back(side%2==0?Exact::subtract(world,bounds[side]):Exact::subtract(bounds[side],world));}
        polygon=anchoredClipValues(polygon,values);
    }
    return polygon;
}
std::vector<AnchoredPoint> clipAnchoredConvexPolygon(const std::vector<AnchoredPoint>& input,Point normal,long double offset)
{
    Exact::validateEnvironment();anchoredRequire(std::isfinite(normal.x) && std::isfinite(normal.y) &&
        (normal.x!=0 || normal.y!=0) && std::isfinite(offset),"invalid polygon cutting plane");
    const int orientation=anchoredOrientation(input.data(),input.size());if(!orientation)return {};
    const auto nx=AnchoredRound::scalar(normal.x),ny=AnchoredRound::scalar(normal.y),d=AnchoredRound::scalar(offset);
    std::vector<AnchoredExactPoint> points;std::vector<std::pair<Expanded,Expanded>> products;
    for(const auto& p:input){points.push_back(anchoredExact(p));products.push_back({nx,points.back().x});products.push_back({ny,points.back().y});}
    const int shift=anchoredProductShift(products,{nx,ny,d});
    const auto sx=Exact::scalePow2(nx,shift),sy=Exact::scalePow2(ny,shift),sd=Exact::scalePow2(d,shift);
    std::vector<Expanded> values;values.reserve(points.size());
    for(const auto& p:points)values.push_back(Exact::subtract(Exact::add(Exact::multiply(sx,p.x),Exact::multiply(sy,p.y)),sd));
    auto result=anchoredClipValues(input,values);if(orientation<0)std::reverse(result.begin(),result.end());return result;
}
bool anchoredPolygonInsideBox(const std::vector<AnchoredPoint>& polygon,Point origin,const Rectangle& box)
{
    Exact::validateEnvironment();anchoredValidateBox(box);validateOrigin(origin);anchoredOrientation(polygon.data(),polygon.size());
    const auto ox=AnchoredRound::scalar(origin.x),oy=AnchoredRound::scalar(origin.y);
    const auto x0=AnchoredRound::scalar(box.lower.x),x1=AnchoredRound::scalar(box.upper.x),y0=AnchoredRound::scalar(box.lower.y),y1=AnchoredRound::scalar(box.upper.y);
    for(const auto& p:polygon){const auto actual=anchoredExact(p);const auto x=Exact::add(actual.x,ox),y=Exact::add(actual.y,oy);
        if(anchoredCompare(x,x0)<0 || anchoredCompare(x,x1)>0 || anchoredCompare(y,y0)<0 || anchoredCompare(y,y1)>0)return false;}
    return true;
}

}
