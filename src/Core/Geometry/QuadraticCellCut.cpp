#include "Core/Geometry/QuadraticCellCut.h"
#include "Core/Mesh/CartesianMeasure.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace Quantum::Geometry::QuadraticCellCut
{
namespace
{
using Poly=std::array<long double,5>;
void require(bool ok,const char* message) {if (!ok) throw std::invalid_argument(message);}
long double coordinate(Point p,int axis) {return axis?p.y:p.x;}
Point center(const Rectangle& b)
{return {static_cast<long double>(b.lower.x)+(static_cast<long double>(b.upper.x)-b.lower.x)/2,
    static_cast<long double>(b.lower.y)+(static_cast<long double>(b.upper.y)-b.lower.y)/2};}
void validate(Shape s)
{require((s.axis==0 || s.axis==1) && (s.sign==-1 || s.sign==1)
    && std::isfinite(s.slope) && std::isfinite(s.second),"[QuadraticCellCut] invalid graph shape");
    if(s.frameNormal) {const auto n=*s.frameNormal;const auto norm=std::hypot(n.x,n.y);
        require(std::isfinite(norm) && norm>0,"[QuadraticCellCut] invalid normal frame");}}
long double value(const Poly& p,long double x)
{long double result=p[4];for (int k=3;k>=0;--k) result=result*x+p[k];return result;}
void roots(long double a,long double b,long double c,long double lo,long double hi,std::vector<long double>& cuts)
{
    require(std::isfinite(a) && std::isfinite(b) && std::isfinite(c),"[QuadraticCellCut] nonfinite root polynomial");
    const auto add=[&](long double x){if (x>lo && x<hi) cuts.push_back(x);};
    if (a==0) {if (b!=0) add(-c/b);return;}
    const auto scale=std::max({std::abs(a),std::abs(b),std::abs(c)});a/=scale;b/=scale;c/=scale;
    const auto product=4*a*c;
    const auto discriminant=std::fma(b,b,-product)-std::fma(4*a,c,-product);
    if (discriminant<0) return;
    const auto q=-.5L*(b+std::copysign(std::sqrt(discriminant),b));
    if (q==0) add(-b/(2*a));else {add(q/a);add(c/q);}
}
// Integrate in a local interval coordinate, avoiding subtraction of large
// primitive values. The physical XY/RZ integrands have degree at most four.
long double integrate(const Poly& p,long double lo,long double hi)
{
    const auto middle=lo+(hi-lo)/2,half=(hi-lo)/2;
    const auto second=p[2]+3*p[3]*middle+6*p[4]*middle*middle;
    return 2*half*(value(p,middle)+second*half*half/3+p[4]*half*half*half*half/5);
}
Poly measure(const Poly& upper,const Poly& lower,bool rz,int axis,long double radius)
{
    Poly difference{};for (int k=0;k<5;++k) difference[k]=upper[k]-lower[k];
    if (!rz) return difference;
    Poly result{};
    if (axis==1) {
        for (int k=0;k<5;++k) {result[k]+=CartesianMeasure::twoPi*radius*difference[k];
            if (k+1<5) result[k+1]+=CartesianMeasure::twoPi*difference[k];}
    } else {
        Poly sum{};for (int k=0;k<5;++k) sum[k]=upper[k]+lower[k];sum[0]+=2*radius;
        for (int i=0;i<5;++i) for (int j=0;i+j<5;++j)
            result[i+j]+=(CartesianMeasure::twoPi/2)*difference[i]*sum[j];
    }
    return result;
}
GraphMeasure measureGraphPolygon(const Rectangle&,const std::vector<Point>&,Point polygonOrigin,
    Point origin,Vector2D normal,long double scale,const std::array<long double,3>&,bool rz,int sign,bool* full);
long double polygonVolume(const Rectangle& box,const std::vector<Point>& polygon,Point origin,
    const Reconstruction& reconstruction,bool rz,bool* entirelyLiquid=nullptr)
{
    if (entirelyLiquid) *entirelyLiquid=true;
    const auto s=reconstruction.shape;const auto c=center(box);
    if(s.frameNormal) {
        const long double scale=std::max(static_cast<long double>(box.upper.x)-box.lower.x,
            static_cast<long double>(box.upper.y)-box.lower.y);
        return measureGraphPolygon(box,polygon,origin,c,*s.frameNormal,scale,
            {{reconstruction.offset/scale,s.slope,s.second*scale/2}},rz,s.sign,entirelyLiquid).volume;
    }
    std::vector<Point> local;std::vector<long double> slabs;
    for (auto p:polygon) {p.x-=c.x-origin.x;p.y-=c.y-origin.y;
        const Point q{coordinate(p,1-s.axis),coordinate(p,s.axis)};local.push_back(q);slabs.push_back(q.x);}
    std::sort(slabs.begin(),slabs.end());slabs.erase(std::unique(slabs.begin(),slabs.end()),slabs.end());
    const Poly graph{{reconstruction.offset,s.slope,s.second/2,0,0}};long double volume=0;
    for (std::size_t k=1;k<slabs.size();++k) {
        const auto lo=slabs[k-1],hi=slabs[k],middle=lo+(hi-lo)/2;
        std::vector<Poly> edges;
        for (std::size_t e=0;e<local.size();++e) {
            const auto a=local[e],b=local[(e+1)%local.size()];
            if (a.x==b.x || middle<=std::min(a.x,b.x) || middle>=std::max(a.x,b.x)) continue;
            const auto slope=(b.y-a.y)/(b.x-a.x);edges.push_back({{a.y-slope*a.x,slope,0,0,0}});
        }
        if (edges.empty()) continue;
        require(edges.size()==2,"[QuadraticCellCut] invalid convex slice");
        if (value(edges[0],middle)>value(edges[1],middle)) std::swap(edges[0],edges[1]);
        std::vector<long double> cuts{lo,hi};
        for (const auto& edge:edges) roots(graph[2],graph[1]-edge[1],graph[0]-edge[0],lo,hi,cuts);
        std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
        for (std::size_t i=1;i<cuts.size();++i) {
            const auto a=cuts[i-1],b=cuts[i],mid=a+(b-a)/2;
            auto lower=edges[0],upper=edges[1];
            if (s.sign>0) {
                if (value(graph,mid)<=value(lower,mid)) {if (entirelyLiquid) *entirelyLiquid=false;continue;}
                if (value(graph,mid)<value(upper,mid)) {upper=graph;if (entirelyLiquid) *entirelyLiquid=false;}
            } else {
                if (value(graph,mid)>=value(upper,mid)) {if (entirelyLiquid) *entirelyLiquid=false;continue;}
                if (value(graph,mid)>value(lower,mid)) {lower=graph;if (entirelyLiquid) *entirelyLiquid=false;}
            }
            volume+=integrate(measure(upper,lower,rz,s.axis,c.x),a,b);
        }
    }
    require(std::isfinite(volume) && volume>=0,"[QuadraticCellCut] invalid curved volume");return volume;
}
}
Basis basis(const Shape& s)
{
    validate(s);
    if(!s.frameNormal)return s.axis?Basis{{1,0},{0,1}}:Basis{{0,1},{1,0}};
    const auto normal=*s.frameNormal;const auto length=std::hypot(static_cast<long double>(normal.x),static_cast<long double>(normal.y));
    const Point n{normal.x/length,normal.y/length};return {{-n.y,n.x},n};
}
long double implicitValue(const Rectangle& box,const Reconstruction& graph,Point point)
{
    require(std::isfinite(box.lower.x) && std::isfinite(box.lower.y) && std::isfinite(box.upper.x)
        && std::isfinite(box.upper.y) && box.lower.x<box.upper.x && box.lower.y<box.upper.y
        && std::isfinite(graph.offset) && std::isfinite(point.x) && std::isfinite(point.y),
        "[QuadraticCellCut] invalid phase query");
    const auto b=basis(graph.shape);const auto c=center(box);point.x-=c.x;point.y-=c.y;
    const auto t=point.x*b.tangent.x+point.y*b.tangent.y,n=point.x*b.normal.x+point.y*b.normal.y;
    return graph.shape.sign*(n-graph.offset-graph.shape.slope*t-graph.shape.second*t*t/2);
}
namespace {
GraphMeasure measureGraphPolygon(const Rectangle& box,const std::vector<Point>& polygon,Point polygonOrigin,
    Point origin,Vector2D normal,long double scale,const std::array<long double,3>& coefficients,bool rz,int sign,bool* fullOutput)
{
    const auto capacity=PlanarCellCut::capacity(box,rz);
    const auto length=std::hypot(static_cast<long double>(normal.x),static_cast<long double>(normal.y));
    require(std::isfinite(scale) && scale>0 && std::isfinite(length) && length>0
        && std::isfinite(origin.x) && std::isfinite(origin.y),"[QuadraticCellCut] invalid rotated graph frame");
    for (auto c:coefficients) require(std::isfinite(c),"[QuadraticCellCut] nonfinite rotated graph coefficient");
    const Point n{normal.x/length,normal.y/length},t{-n.y,n.x};
    std::vector<Point> local;local.reserve(polygon.size());
    std::vector<long double> slabs;
    for (auto point:polygon) {
        const auto x=(point.x-(origin.x-polygonOrigin.x))/scale,y=(point.y-(origin.y-polygonOrigin.y))/scale;
        local.push_back({x*t.x+y*t.y,x*n.x+y*n.y});slabs.push_back(local.back().x);
    }
    std::sort(slabs.begin(),slabs.end());slabs.erase(std::unique(slabs.begin(),slabs.end()),slabs.end());
    GraphMeasure result;bool full=true,empty=true;
    for (std::size_t k=1;k<slabs.size();++k) {
        const auto lo=slabs[k-1],hi=slabs[k],width=hi-lo;
        // Integrate in the slab's own [0,1] parameter. Nearly aligned frames
        // can create sub-ULP-width slivers: no rounded interior abscissa or
        // enormous global edge intercept is needed to represent their area.
        const Poly graph{{(coefficients[2]*lo+coefficients[1])*lo+coefficients[0],
            (coefficients[1]+2*coefficients[2]*lo)*width,coefficients[2]*width*width,0,0}};
        std::vector<Poly> edges;
        for (std::size_t i=0;i<local.size();++i) {
            const auto a=local[i],b=local[(i+1)%local.size()];
            if (a.x==b.x || lo<std::min(a.x,b.x) || hi>std::max(a.x,b.x)) continue;
            edges.push_back({{a.y+(lo-a.x)/(b.x-a.x)*(b.y-a.y),width/(b.x-a.x)*(b.y-a.y),0,0,0}});
        }
        if (edges.empty()) continue;
        require(edges.size()==2,"[QuadraticCellCut] invalid rotated cell slice");
        if (value(edges[0],.5L)>value(edges[1],.5L)) std::swap(edges[0],edges[1]);
        std::vector<long double> cuts{0,1};
        for (const auto& edge:edges) roots(graph[2],graph[1]-edge[1],graph[0]-edge[0],0,1,cuts);
        std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
        for (std::size_t j=1;j<cuts.size();++j) {
            const auto a=cuts[j-1],b=cuts[j],mid=a+(b-a)/2;
            auto lower=edges[0],upper=edges[1];bool cut=false;
            if(sign>0) {
                if(value(graph,mid)<=value(lower,mid)){full=false;continue;}
                cut=value(graph,mid)<value(upper,mid);if(cut)upper=graph;
            } else {
                if(value(graph,mid)>=value(upper,mid)){full=false;continue;}
                cut=value(graph,mid)>value(lower,mid);if(cut)lower=graph;
            }
            empty=false;if(cut)full=false;
            Poly difference{},density{};
            for(int i=0;i<5;++i) difference[i]=upper[i]-lower[i];
            if (rz) {
                Poly meanRadius{};
                for(int i=0;i<5;++i) meanRadius[i]=scale*n.x*(upper[i]+lower[i])/2;
                meanRadius[0]+=origin.x+scale*t.x*lo;meanRadius[1]+=scale*t.x*width;
                for(int i=0;i<5;++i) for(int m=0;i+m<5;++m)
                    density[i+m]+=CartesianMeasure::twoPi*difference[i]*meanRadius[m];
            } else density=difference;
            result.volume+=scale*scale*width*integrate(density,a,b);
            if (cut) {
                // Endpoint variations cancel because the two slice formulas
                // agree where the graph meets the rectangle boundary.
                Poly radius{};radius[0]=1;
                if (rz) {
                    radius={};for(int i=0;i<3;++i) radius[i]=scale*n.x*graph[i];
                    radius[0]+=origin.x+scale*t.x*lo;radius[1]+=scale*t.x*width;
                    for(auto& q:radius) q*=CartesianMeasure::twoPi;
                }
                const std::array<Poly,3> basis{{Poly{{1,0,0,0,0}},Poly{{lo,width,0,0,0}},
                    Poly{{lo*lo,2*lo*width,width*width,0,0}}}};
                for(int m=0;m<3;++m) {
                    Poly sensitivity{};
                    for(int i=0;i<5;++i) for(int p=0;i+p<5;++p) sensitivity[i+p]+=radius[i]*basis[m][p];
                    result.derivative[m]+=sign*scale*scale*width*integrate(sensitivity,a,b);
                }
            }
        }
    }
    if(fullOutput)*fullOutput=full;if(empty)result.volume=0;
    const auto allowance=PlanarCellCut::relativeVolumeResolution*capacity;
    require(std::isfinite(result.volume) && result.volume>=-allowance && result.volume<=capacity+allowance,
        "[QuadraticCellCut] invalid rotated graph volume");
    for(auto d:result.derivative) require(std::isfinite(d),"[QuadraticCellCut] nonfinite graph sensitivity");
    return result;
}
}
GraphMeasure graphMeasure(const Rectangle& box,Point origin,Vector2D normal,long double scale,
    const std::array<long double,3>& coefficients,bool rz)
{
    const std::vector<Point> corners{{box.lower.x,box.lower.y},{box.upper.x,box.lower.y},
        {box.upper.x,box.upper.y},{box.lower.x,box.upper.y}};
    bool full=false;auto result=measureGraphPolygon(box,corners,{0,0},origin,normal,scale,coefficients,rz,1,&full);
    if(full)result.volume=PlanarCellCut::capacity(box,rz);return result;
}
PlanarCellCut::Overlap overlapInFrame(const Rectangle& box,const std::vector<Point>& polygon,Point origin,
    const Reconstruction& reconstruction,bool rz)
{
    validate(reconstruction.shape);require(std::isfinite(reconstruction.offset),"[QuadraticCellCut] invalid offset");
    const auto total=PlanarCellCut::overlapInFrame(box,polygon,origin,{1,0},0,rz).total;
    const auto clipped=PlanarCellCut::clipPolygonInFrame(box,polygon,origin);
    bool entirelyLiquid=false;
    const auto integrated=polygonVolume(box,clipped,origin,reconstruction,rz,&entirelyLiquid);
    // A region classified wholly liquid shares the exact total metric. This
    // avoids comparing different integration roundoff for the same domain;
    // genuinely cut regions retain their independently integrated volume.
    const auto liquid=entirelyLiquid?total:integrated;
    require(liquid<=total+PlanarCellCut::relativeVolumeResolution*PlanarCellCut::capacity(box,rz),
        "[QuadraticCellCut] curved subvolume exceeds donor overlap");
    // No phase clipping: the caller retains the same strict shared-face bounds.
    return {total,liquid};
}
Reconstruction reconstruct(const Rectangle& box,Shape shape,long double fraction,bool rz)
{
    validate(shape);require(std::isfinite(fraction) && fraction>=0 && fraction<=1,"[QuadraticCellCut] invalid fraction");
    const auto capacity=PlanarCellCut::capacity(box,rz);const auto c=center(box);
    const long double hx=(static_cast<long double>(box.upper.x)-box.lower.x)/2,hy=(static_cast<long double>(box.upper.y)-box.lower.y)/2;
    const auto frame=basis(shape);
    const auto ht=std::abs(frame.tangent.x)*hx+std::abs(frame.tangent.y)*hy,
        hn=std::abs(frame.normal.x)*hx+std::abs(frame.normal.y)*hy;
    const auto extent=hn+std::abs(shape.slope)*ht+std::abs(shape.second)*ht*ht/2;
    const std::vector<Point> polygon{{-hx,-hy},{hx,-hy},{hx,hy},{-hx,hy}};
    Reconstruction result{shape};long double lo=-extent,hi=extent;const auto target=fraction*capacity;
    for (int i=0;i<96;++i) {
        result.offset=lo+(hi-lo)/2;bool full=false;result.volume=polygonVolume(box,polygon,c,result,rz,&full);
        if(full)result.volume=capacity;
        if (result.volume==target || result.offset==lo || result.offset==hi) break;
        if ((result.volume<target)==(shape.sign>0)) lo=result.offset;else hi=result.offset;
    }
    result.residual=result.volume-target;
    require(std::abs(result.residual)<=PlanarCellCut::relativeVolumeResolution*capacity,
        "[QuadraticCellCut] unresolved curved volume constraint");return result;
}
std::optional<Reconstruction> projectEndpointVolume(const Rectangle& box,int axis,int sign,
    long double lowerHeight,long double upperHeight,long double liquidVolume,bool rz)
{
    validate({axis,sign,0,0});const auto capacity=PlanarCellCut::capacity(box,rz);const auto c=center(box);
    require(std::isfinite(lowerHeight) && std::isfinite(upperHeight) && std::isfinite(liquidVolume)
        && liquidVolume>=0 && liquidVolume<=capacity,"[QuadraticCellCut] invalid projection data");
    const long double width=axis?static_cast<long double>(box.upper.x)-box.lower.x:
        static_cast<long double>(box.upper.y)-box.lower.y;
    const long double nlo=axis?box.lower.y:box.lower.x,nhi=axis?box.upper.y:box.upper.x;
    if(!(liquidVolume>0 && liquidVolume<capacity) || lowerHeight<nlo || lowerHeight>nhi
        || upperHeight<nlo || upperHeight>nhi) return std::nullopt;
    const auto average=(lowerHeight+upperHeight)/2,slope=(upperHeight-lowerHeight)/width;
    // H(t)=average+slope*t+k*(t*t/width^2-1/4). Endpoints stay fixed.
    long double k=0;
    if(rz && axis==0) {
        const auto meanSquare=sign>0?nlo*nlo+liquidVolume/((CartesianMeasure::twoPi/2)*width):
            nhi*nhi-liquidVolume/((CartesianMeasure::twoPi/2)*width);
        const auto difference=(meanSquare-average*average)-slope*slope*width*width/12;
        const auto discriminant=100*average*average+120*difference;
        if(!(discriminant>0) || !(average>0)) return std::nullopt;
        // Physical root continuous with k=0, without subtracting close roots.
        k=-60*difference/(10*average+std::sqrt(discriminant));
    } else {
        const auto measure=rz?CartesianMeasure::twoPi*c.x*width:width;
        const auto mean=sign>0?nlo+liquidVolume/measure:nhi-liquidVolume/measure;
        const auto linearMean=average+(rz?slope*width*width/(12*c.x):0);
        k=6*(linearMean-mean);
    }
    const auto middle=average-k/4,second=2*k/(width*width);
    if(!std::isfinite(middle) || !std::isfinite(second)) return std::nullopt;
    if(second!=0) {
        const auto stationary=-slope/second;
        if(std::abs(stationary)<width/2) {
            const auto value=middle+slope*stationary+second*stationary*stationary/2;
            if(value<nlo || value>nhi) return std::nullopt;
        }
    }
    Reconstruction result{{axis,sign,slope,second},middle-coordinate(c,axis)};
    const long double hx=(static_cast<long double>(box.upper.x)-box.lower.x)/2,
        hy=(static_cast<long double>(box.upper.y)-box.lower.y)/2;
    result.volume=polygonVolume(box,{{-hx,-hy},{hx,-hy},{hx,hy},{-hx,hy}},c,result,rz);
    result.residual=result.volume-liquidVolume;
    if(std::abs(result.residual)>PlanarCellCut::relativeVolumeResolution*capacity) return std::nullopt;
    return result;
}
std::vector<std::array<long double,2>> intervals(const Rectangle& box,const Reconstruction& graph,Point begin,Point end,
    long double lower,long double upper)
{
    validate(graph.shape);require(lower>=0 && upper<=1 && lower<=upper,"[QuadraticCellCut] invalid chord interval");
    const auto c=center(box);begin.x-=c.x;begin.y-=c.y;end.x-=c.x;end.y-=c.y;
    const auto s=graph.shape;const auto frame=basis(s);
    const auto t=begin.x*frame.tangent.x+begin.y*frame.tangent.y,
        dt=(end.x-begin.x)*frame.tangent.x+(end.y-begin.y)*frame.tangent.y;
    const auto n=begin.x*frame.normal.x+begin.y*frame.normal.y,
        dn=(end.x-begin.x)*frame.normal.x+(end.y-begin.y)*frame.normal.y;
    const auto a=-s.second*dt*dt/2,b=dn-s.slope*dt-s.second*t*dt,d=n-graph.offset-s.slope*t-s.second*t*t/2;
    std::vector<long double> cuts{lower,upper};roots(a,b,d,lower,upper,cuts);
    std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
    std::vector<std::array<long double,2>> result;
    for (std::size_t i=1;i<cuts.size();++i) {
        const auto mid=cuts[i-1]+(cuts[i]-cuts[i-1])/2;
        if (s.sign*((a*mid+b)*mid+d)<=0) result.push_back({{cuts[i-1],cuts[i]}});
    }
    return result;
}
}
