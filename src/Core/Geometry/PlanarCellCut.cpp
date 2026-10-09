#include "Core/Geometry/PlanarCellCut.h"
#include "Core/Mesh/CartesianMeasure.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace Quantum::Geometry::PlanarCellCut
{
namespace
{
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
static std::vector<Point> clipRectangle(const std::vector<Point>& input,Point lower,Point upper)
{
    for (const auto p:input) if (!std::isfinite(p.x) || !std::isfinite(p.y))
        throw std::invalid_argument("[PlanarCellCut] nonfinite swept polygon");
    if (input.size()<3) return {};
    int orientation=0;
    // Every vertex must lie on the same side of every oriented polygon edge.
    // This also excludes crossed polygons; repeated/collinear vertices are OK.
    for (std::size_t i=0;i<input.size();++i) {
        const auto a=input[i],b=input[(i+1)%input.size()];
        for (const auto c:input) {
            const auto cross=(b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x);
            if (!std::isfinite(cross)) throw std::invalid_argument("[PlanarCellCut] polygon overflow");
            const int sign=cross>0?1:cross<0?-1:0;
            if (sign && orientation && sign!=orientation)
                throw std::invalid_argument("[PlanarCellCut] polygon is not convex and cyclic");
            if (sign) orientation=sign;
        }
    }
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
}
