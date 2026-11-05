#include "Core/Geometry/QuadraticSurfaceMeasure.h"
#include "Core/Geometry/QuadraticSurfacePreparation.h"
#include "Core/Geometry/QuadraticFixedCornerSensitivity.h"
#include "Core/Geometry/Detail/QuadraticGraphExact.h"
#include "Core/Geometry/Detail/ExactLinearNativeFilter.h"
#include "Core/Mesh/CartesianMeasure.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <chrono>
#include <limits>
#include <optional>
#include <stdexcept>

namespace LAMP::Geometry::QuadraticSurfaceMeasure
{
namespace
{
void require(bool ok, const char* message) { if(!ok) throw std::invalid_argument(message); }
using SlopeRoundingFilter=Detail::ExactLinearNative::Filter;
struct Frame
{
    Point centre, tangent, normal;
    long double offset, slope, second;
    int sign;
    std::shared_ptr<const Detail::QuadraticGraphExact> exact;
    long double metric=1;
    SlopeRoundingFilter slopeFilter;
    long double slopeAt(long double t) const
    {
        if(!exact)return slope+second*t;
        long double result=0;if(slopeFilter.evaluate(t,result))return result;
        return Detail::QuadraticExact::rounded(Detail::graphSlope(*exact,Detail::QuadraticExact::scalar(t)));
    }
    Point point(long double t) const
    {
        if(exact){for(std::size_t i=0;i<exact->portParameters.size();++i)if(Detail::QuadraticExact::rounded(exact->portParameters[i])==t)return exact->ports[i];return Detail::QuadraticExact::native(Detail::graphPosition(*exact,Detail::QuadraticExact::scalar(t)));}
        const auto h=std::fma(second*t/2,t,std::fma(slope,t,offset));
        return {centre.x+tangent.x*t+normal.x*h,centre.y+tangent.y*t+normal.y*h};
    }
};
Frame frame(const Rectangle& source,const Graph& graph,bool rz)
{
    (void)PlanarCellCut::capacity(source,rz);
    require(std::isfinite(graph.offset),"[QuadraticSurfaceMeasure] nonfinite offset");
    Detail::validateGraph(graph);
    if(graph.exact){const auto&e=*graph.exact;Frame result{Detail::QuadraticExact::native(e.origin),Detail::QuadraticExact::native(e.tangent),Detail::QuadraticExact::native(e.normal),graph.offset,graph.shape.slope,graph.shape.second,graph.shape.sign,graph.exact,std::sqrt(Detail::QuadraticExact::rounded(e.determinant))};
        result.slopeFilter=SlopeRoundingFilter::prepare(e.q[1],2*e.q[2]);return result;}
    const auto b=QuadraticCellCut::basis(graph.shape);
    return {{source.lower.x+(static_cast<long double>(source.upper.x)-source.lower.x)/2,
        source.lower.y+(static_cast<long double>(source.upper.y)-source.lower.y)/2},
        b.tangent,b.normal,graph.offset,graph.shape.slope,graph.shape.second,graph.shape.sign};
}
struct Primitive { long double zero,first,second; };
Primitive primitive(long double z)
{
    const auto h=std::hypot(1.0L,z),a=std::asinh(z);
    return {(z*h+a)/2,h*h*h/3,(z*h*(2*z*z+1)-a)/8};
}
struct Pair { long double length=0,area=0; };
constexpr std::array<long double,4> gaussX{{.18343464249564980494L,.52553240991632898582L,
    .79666647741362673959L,.96028985649753623168L}};
constexpr std::array<long double,4> gaussW{{.36268378337836198297L,.31370664587788728734L,
    .22238103445337447054L,.10122853629037625915L}};
// Fixed positive rule for the SAME arc integral in the well-conditioned
// interval coordinate. Only conditioning selects this evaluation; no curve
// coefficient, phase or curvature is deleted in the nearly straight limit.
Pair gauss(const Frame& f,long double lo,long double hi,bool rz)
{
    const auto mid=lo+(hi-lo)/2,half=(hi-lo)/2;
    Pair result;
    for(std::size_t i=0;i<gaussX.size();++i) for(int sign:{-1,1})
    {
        const auto t=mid+sign*half*gaussX[i];
        const auto length=half*gaussW[i]*f.metric*std::hypot(1.0L,f.slopeAt(t));
        const auto radius=rz?f.point(t).x:0.L;
        require(!rz || radius>=0,"[QuadraticSurfaceMeasure] arc crosses negative radius");
        result.length+=length;result.area+=rz?CartesianMeasure::twoPi*radius*length:length;
    }
    return result;
}
struct Derivatives { std::array<long double,3> value{}, absolute{}; };
Derivatives gaussDerivatives(const Frame& f,long double lo,long double hi,bool rz)
{
    const auto mid=lo+(hi-lo)/2,half=(hi-lo)/2;
    Derivatives result;
    for(std::size_t i=0;i<gaussX.size();++i) for(int sign:{-1,1})
    {
        const auto t=mid+sign*half*gaussX[i],p=f.slopeAt(t),g=std::hypot(1.0L,p);
        const auto radius=rz?f.point(t).x:1;
        const auto multiplier=rz?static_cast<long double>(CartesianMeasure::twoPi):1;
        const std::array<long double,3> height{{1,t,t*t/2}},slope{{0,1,t}};
        for(int j=0;j<3;++j)
        {
            const auto derivative=f.metric*multiplier*((rz?f.normal.x*height[j]*g:0)+radius*p*slope[j]/g);
            const auto weighted=half*gaussW[i]*derivative;
            result.value[j]+=weighted;result.absolute[j]+=std::abs(weighted);
        }
    }
    return result;
}
std::array<long double,3> derivativeIntegral(const Frame& f,long double lo,long double hi,bool rz,unsigned depth=0)
{
    const auto mid=lo+(hi-lo)/2;
    const auto whole=gaussDerivatives(f,lo,hi,rz),a=gaussDerivatives(f,lo,mid,rz),b=gaussDerivatives(f,mid,hi,rz);
    std::array<long double,3> result{};bool accepted=true;
    for(int j=0;j<3;++j)
    {
        result[j]=a.value[j]+b.value[j];
        accepted &= std::abs(result[j]-whole.value[j])<=512*std::numeric_limits<long double>::epsilon()*(a.absolute[j]+b.absolute[j]);
    }
    if(accepted)return result;
    require(depth<16 && mid!=lo && mid!=hi,"[QuadraticSurfaceMeasure] unresolved area sensitivity");
    const auto left=derivativeIntegral(f,lo,mid,rz,depth+1),right=derivativeIntegral(f,mid,hi,rz,depth+1);
    for(int j=0;j<3;++j) result[j]=left[j]+right[j];return result;
}
Integral stableIntegral(const Frame& f,long double lo,long double hi,bool rz,unsigned depth=0)
{
    const auto whole=gauss(f,lo,hi,rz);
    const auto mid=lo+(hi-lo)/2;
    const auto left=gauss(f,lo,mid,rz),right=gauss(f,mid,hi,rz);
    const Pair split{left.length+right.length,left.area+right.area};
    const auto dl=std::abs(split.length-whole.length),da=std::abs(split.area-whole.area);
    constexpr auto tolerance=512*std::numeric_limits<long double>::epsilon();
    if(dl<=tolerance*split.length && da<=tolerance*split.area)
        return {split.length,split.area,da};
    require(depth<16 && mid!=lo && mid!=hi,"[QuadraticSurfaceMeasure] unresolved arc integral");
    const auto a=stableIntegral(f,lo,mid,rz,depth+1),b=stableIntegral(f,mid,hi,rz,depth+1);
    return {a.length+b.length,a.physicalMeasure+b.physicalMeasure,
        a.quadratureErrorEstimate+b.quadratureErrorEstimate};
}
Integral integral(const Frame& f,long double lo,long double hi,bool rz,
    std::optional<long double> radiusLo=std::nullopt,std::optional<long double> radiusHi=std::nullopt)
{
    require(std::isfinite(lo) && std::isfinite(hi) && hi>=lo,"[QuadraticSurfaceMeasure] invalid tangent interval");
    if(lo==hi) return {};
    if(rz)
    {
        // Only evaluate() supplies these values, from an actual clipping-side
        // root identity. Untagged integrate() retains its original strict test.
        const auto firstRadius=radiusLo.value_or(f.point(lo).x),lastRadius=radiusHi.value_or(f.point(hi).x);
        auto least=std::min(firstRadius,lastRadius);
        if(f.normal.x*f.second>0)
        {
            const auto stationary=-(f.tangent.x+f.normal.x*f.slope)/(f.normal.x*f.second);
            if(stationary>lo && stationary<hi) least=std::min(least,f.point(stationary).x);
        }
        require(least>=0,"[QuadraticSurfaceMeasure] arc crosses negative radius");
    }
    const auto mid=lo+(hi-lo)/2,half=(hi-lo)/2;
    const auto slope=f.slopeAt(mid);
    if(f.second==0)
    {
        const auto length=(hi-lo)*f.metric*std::hypot(1.0L,slope),radius=f.point(mid).x;
        require(!rz || (radiusLo.value_or(f.point(lo).x)>=0 && radiusHi.value_or(f.point(hi).x)>=0),"[QuadraticSurfaceMeasure] negative-radius line");
        return {length,rz?CartesianMeasure::twoPi*radius*length:length,0};
    }
    const auto variation=std::abs(f.second)*half;
    // Avoid catastrophic cancellation of primitive differences and 1/k^3
    // weighted moments. The positive local rule keeps nonzero k unchanged.
    if(variation<.5L*std::hypot(1.0L,slope)) return stableIntegral(f,lo,hi,rz);
    const auto a=primitive(slope-f.second*half),b=primitive(slope+f.second*half);
    const auto d0=b.zero-a.zero,d1=b.first-a.first,d2=b.second-a.second;
    const auto i0=d0/f.second;
    Integral result{f.metric*i0,f.metric*i0,0};
    if(rz)
    {
        const auto i1=(d1-slope*d0)/(f.second*f.second);
        const auto i2=(d2-2*slope*d1+slope*slope*d0)/(f.second*f.second*f.second);
        const auto r0=f.point(mid).x,r1=f.tangent.x+f.normal.x*slope,r2=f.normal.x*f.second/2;
        result.physicalMeasure=f.metric*CartesianMeasure::twoPi*(r0*i0+r1*i1+r2*i2);
    }
    require(std::isfinite(result.length) && result.length>0 && std::isfinite(result.physicalMeasure)
        && result.physicalMeasure>=0,"[QuadraticSurfaceMeasure] unrepresentable surface measure");
    return result;
}
struct Event { long double t=0; unsigned mask=0; };
void roots(long double a,long double b,long double c,long double lo,long double hi,
    unsigned mask,std::vector<Event>& events)
{
    require(std::isfinite(a) && std::isfinite(b) && std::isfinite(c),"[QuadraticSurfaceMeasure] nonfinite side equation");
    const auto add=[&](long double t){if(t>=lo && t<=hi) events.push_back({t,mask});};
    if(a==0) { if(b!=0) add(-c/b);return; }
    const auto scale=std::max({std::abs(a),std::abs(b),std::abs(c)});
    a/=scale;b/=scale;c/=scale;
    const auto product=4*a*c;
    const auto discriminant=std::fma(b,b,-product)-std::fma(4*a,c,-product);
    if(discriminant<0) return;
    const auto q=-.5L*(b+std::copysign(std::sqrt(discriminant),b));
    if(q==0) add(-b/(2*a));else {add(q/a);add(c/q);}
}
}
Sample sample(const Rectangle& source,const Graph& graph,long double t,bool rz)
{
    const auto f=frame(source,graph,rz);
    require(std::isfinite(t),"[QuadraticSurfaceMeasure] nonfinite tangent");
    const auto slope=f.slopeAt(t),g=std::hypot(1.0L,slope);
    Sample result;
    result.position=f.point(t);
    result.unitNormal={f.sign*(f.normal.x-slope*f.tangent.x)/(f.metric*g),
        f.sign*(f.normal.y-slope*f.tangent.y)/(f.metric*g)};
    result.meridionalCurvature=-f.sign*f.second/(f.metric*g*g*g);
    result.curvature=result.meridionalCurvature;
    if(rz)
    {
        require(result.position.x>=0,"[QuadraticSurfaceMeasure] negative sample radius");
        if(result.position.x>0) result.curvature+=result.unitNormal.x/result.position.x;
        else
        {
            require(std::abs(result.unitNormal.x)<=64*std::numeric_limits<long double>::epsilon() &&
                f.tangent.x+f.normal.x*slope!=0,"[QuadraticSurfaceMeasure] nonregular axis intersection");
            result.curvature+=result.meridionalCurvature;
        }
    }
    require(std::isfinite(result.position.x) && std::isfinite(result.position.y) && std::isfinite(result.curvature),
        "[QuadraticSurfaceMeasure] nonfinite graph sample");
    return result;
}
Integral integrate(const Rectangle& source,const Graph& graph,long double lo,long double hi,bool rz)
{ return integral(frame(source,graph,rz),lo,hi,rz); }
namespace
{
Result evaluateImpl(const Rectangle& source,const Graph& graph,const Rectangle& region,bool rz,unsigned ownership,MaterialCurves* materialCurves,SurfacePreparationTiming* timing=nullptr)
{
    const auto f=frame(source,graph,rz);(void)PlanarCellCut::capacity(region,rz);
    require(region.lower.x>=source.lower.x && region.lower.y>=source.lower.y &&
        region.upper.x<=source.upper.x && region.upper.y<=source.upper.y && (ownership&~All)==0,
        "[QuadraticSurfaceMeasure] region outside original cell or invalid ownership");
    if(graph.exact){
        Result result;std::vector<Detail::GraphEvent> events;
        const auto eventStart=timing?std::chrono::steady_clock::now():std::chrono::steady_clock::time_point{};
        const auto arcs=Detail::graphArcs(*graph.exact,region,ownership,materialCurves?&events:nullptr);
        if(timing)timing->eventSeconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-eventStart).count();
        for(const auto&arc:arcs){
            Patch p;p.lower=Detail::QuadraticExact::rounded(arc.first.parameter);p.upper=Detail::QuadraticExact::rounded(arc.last.parameter);
            require(p.upper>p.lower,"[QuadraticSurfaceMeasure] native interval cannot represent distinct exact events");
            p.begin=arc.first.position;p.end=arc.last.position;p.beginBoundary=arc.first.boundary;p.endBoundary=arc.last.boundary;
            const auto integralStart=timing?std::chrono::steady_clock::now():std::chrono::steady_clock::time_point{};
            p.integral=integral(f,p.lower,p.upper,rz,p.begin.x,p.end.x);
            if(timing)timing->areaIntegralSeconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-integralStart).count();
            result.integral.length+=p.integral.length;result.integral.physicalMeasure+=p.integral.physicalMeasure;result.integral.quadratureErrorEstimate+=p.integral.quadratureErrorEstimate;result.patches.push_back(p);
            if(materialCurves){
                const auto curveStart=timing?std::chrono::steady_clock::now():std::chrono::steady_clock::time_point{};
                const Detail::GraphEvent *first=nullptr,*last=nullptr;
                for(const auto&event:events){if(!event.boundary)continue;
                    if(Detail::QuadraticExact::rounded(event.parameter)==p.lower){require(!first,"[QuadraticCellCut] ambiguous native material lower event");first=&event;}
                    if(Detail::QuadraticExact::rounded(event.parameter)==p.upper){require(!last,"[QuadraticCellCut] ambiguous native material upper event");last=&event;}}
                require(first&&last&&first->parameter<last->parameter,"[QuadraticCellCut] material interval has no original side-root identities");
                materialCurves->push_back(Detail::graphArcMaterialData(*graph.exact,{*first,*last}));
                if(timing)timing->materialCurveSeconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-curveStart).count();
            }
        }return result;
    }
    long double lo=std::numeric_limits<long double>::infinity(),hi=-lo;
    for(auto x:{region.lower.x,region.upper.x}) for(auto y:{region.lower.y,region.upper.y})
    {
        const auto t=(x-f.centre.x)*f.tangent.x+(y-f.centre.y)*f.tangent.y;
        lo=std::min(lo,t);hi=std::max(hi,t);
    }
    std::vector<Event> events{{lo,0},{hi,0}};
    unsigned coincident=0;
    const auto side=[&](long double n,long double t,long double centre,long double border,unsigned mask)
    {
        const auto a=n*f.second/2,b=t+n*f.slope,c=centre+n*f.offset-border;
        if(a==0 && b==0 && c==0) coincident|=mask;
        else roots(a,b,c,lo,hi,mask,events);
    };
    side(f.normal.x,f.tangent.x,f.centre.x,region.lower.x,Left);
    side(f.normal.x,f.tangent.x,f.centre.x,region.upper.x,Right);
    side(f.normal.y,f.tangent.y,f.centre.y,region.lower.y,Bottom);
    side(f.normal.y,f.tangent.y,f.centre.y,region.upper.y,Top);
    Result result;
    if(coincident&~ownership) return result;
    std::sort(events.begin(),events.end(),[](auto a,auto b){return a.t<b.t;});
    std::vector<Event> cuts;
    for(const auto e:events)
        if(!cuts.empty() && cuts.back().t==e.t) cuts.back().mask|=e.mask;
        else cuts.push_back(e);
    for(std::size_t i=1;i<cuts.size();++i)
    {
        const auto a=cuts[i-1],b=cuts[i];
        const auto p=f.point(a.t+(b.t-a.t)/2);
        if(p.x<region.lower.x || p.x>region.upper.x || p.y<region.lower.y || p.y>region.upper.y) continue;
        // Adjacent inside intervals are one connected arc, including a
        // retained tangential touch. No close roots are merged or snapped.
        if(!result.patches.empty() && result.patches.back().upper==a.t)
        {
            auto& patch=result.patches.back();patch.upper=b.t;patch.end=f.point(b.t);patch.endBoundary=b.mask|coincident;
        }
        else result.patches.push_back({a.t,b.t,f.point(a.t),f.point(b.t),a.mask|coincident,b.mask|coincident,{}});
    }
    for(auto& p:result.patches)
    {
        const auto provenRadius=[&](unsigned mask)->std::optional<long double>{
            if(mask&Left)return region.lower.x;if(mask&Right)return region.upper.x;return std::nullopt;};
        p.integral=integral(f,p.lower,p.upper,rz,provenRadius(p.beginBoundary),provenRadius(p.endBoundary));
        result.integral.length+=p.integral.length;
        result.integral.physicalMeasure+=p.integral.physicalMeasure;
        result.integral.quadratureErrorEstimate+=p.integral.quadratureErrorEstimate;
    }
    return result;
}
Sensitivity differentiateImpl(const Rectangle& source,const Graph& graph,const Rectangle& region,bool rz,unsigned ownership,MaterialCurves* materialCurves,SurfacePreparationTiming* timing=nullptr)
{
    const auto f=frame(source,graph,rz);
    Sensitivity result;result.surface=evaluateImpl(source,graph,region,rz,ownership,materialCurves,timing);
    const auto derivativeStart=timing?std::chrono::steady_clock::now():std::chrono::steady_clock::time_point{};
    for(const auto& patch:result.surface.patches)
    {
        const auto inside=derivativeIntegral(f,patch.lower,patch.upper,rz);
        for(int j=0;j<3;++j) result.derivative[j]+=inside[j];
        const auto boundary=[&](long double t,unsigned mask,int endpointSign)
        {
            require(mask!=0 && (mask&(mask-1))==0,"[QuadraticSurfaceMeasure] area inverse has a corner/coincident endpoint");
            const auto n=(mask&(Left|Right))?f.normal.x:f.normal.y;
            const auto tangent=(mask&(Left|Right))?f.tangent.x:f.tangent.y;
            const auto p=f.slopeAt(t),denominator=tangent+n*p;
            require(denominator!=0,"[QuadraticSurfaceMeasure] tangential endpoint derivative singular");
            // The active radial side is the exact endpoint constraint. Its
            // r=0 boundary density vanishes even if the displayed raw root
            // coordinate has rounded just outside the physical rectangle.
            const auto radius=(mask&Left)?static_cast<long double>(region.lower.x):(mask&Right)?static_cast<long double>(region.upper.x):f.point(t).x;
            const auto density=f.metric*std::hypot(1.0L,p)*(rz?CartesianMeasure::twoPi*radius:1);
            const std::array<long double,3> basis{{1,t,t*t/2}};
            for(int j=0;j<3;++j) result.derivative[j]-=endpointSign*density*n*basis[j]/denominator;
        };
        boundary(patch.lower,patch.beginBoundary,-1);boundary(patch.upper,patch.endBoundary,1);
    }
    if(timing)timing->derivativeSeconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-derivativeStart).count();
    for(auto d:result.derivative)require(std::isfinite(d),"[QuadraticSurfaceMeasure] nonfinite area sensitivity");
    return result;
}
Sensitivity differentiateWithFixedCornerEndpointsImpl(const Rectangle& source,const Graph& graph,const Rectangle& region,bool rz,unsigned ownership,MaterialCurves* materialCurves,SurfacePreparationTiming* timing=nullptr)
{
    require(graph.exact!=nullptr,"[QuadraticSurfaceMeasure/FixedCorner] exact graph required");
    const auto f=frame(source,graph,rz);
    Sensitivity result;result.surface=evaluateImpl(source,graph,region,rz,ownership,materialCurves,timing);
    const auto derivativeStart=timing?std::chrono::steady_clock::now():std::chrono::steady_clock::time_point{};
    require(result.surface.patches.size()==1,"[QuadraticSurfaceMeasure/FixedCorner] one ordinary arc required");
    const auto& patch=result.surface.patches.front();
    const auto corner=[](unsigned mask){return (mask&(Left|Right))!=0&&(mask&(Bottom|Top))!=0;};
    require(corner(patch.beginBoundary)||corner(patch.endBoundary),"[QuadraticSurfaceMeasure/FixedCorner] actual corner endpoint required");
    const auto inside=derivativeIntegral(f,patch.lower,patch.upper,rz);
    for(int j=0;j<3;++j)result.derivative[j]+=inside[j];
    const auto boundary=[&](long double t,Point point,unsigned mask,int endpointSign)
    {
        if(corner(mask)){
            const bool left=(mask&Left)!=0,right=(mask&Right)!=0,bottom=(mask&Bottom)!=0,top=(mask&Top)!=0;
            require(left!=right&&bottom!=top&&point.x==(left?region.lower.x:region.upper.x)&&
                point.y==(bottom?region.lower.y:region.upper.y),"[QuadraticSurfaceMeasure/FixedCorner] invalid actual corner event");
            // Only the declared constrained pullback is meaningful here:
            // delta F(corner)=0 keeps this same physical parameter fixed.
            return;
        }
        require(mask!=0&&(mask&(mask-1))==0,"[QuadraticSurfaceMeasure] area inverse has a corner/coincident endpoint");
        const auto n=(mask&(Left|Right))?f.normal.x:f.normal.y;
        const auto tangent=(mask&(Left|Right))?f.tangent.x:f.tangent.y;
        const auto p=f.slopeAt(t),denominator=tangent+n*p;
        require(denominator!=0,"[QuadraticSurfaceMeasure] tangential endpoint derivative singular");
        const auto radius=(mask&Left)?static_cast<long double>(region.lower.x):(mask&Right)?static_cast<long double>(region.upper.x):f.point(t).x;
        const auto density=f.metric*std::hypot(1.0L,p)*(rz?CartesianMeasure::twoPi*radius:1);
        const std::array<long double,3> basis{{1,t,t*t/2}};
        for(int j=0;j<3;++j)result.derivative[j]-=endpointSign*density*n*basis[j]/denominator;
    };
    boundary(patch.lower,patch.begin,patch.beginBoundary,-1);boundary(patch.upper,patch.end,patch.endBoundary,1);
    if(timing)timing->derivativeSeconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-derivativeStart).count();
    for(auto d:result.derivative)require(std::isfinite(d),"[QuadraticSurfaceMeasure] nonfinite area sensitivity");
    return result;
}
}
Result evaluate(const Rectangle& source,const Graph& graph,const Rectangle& region,bool rz,unsigned ownership)
{return evaluateImpl(source,graph,region,rz,ownership,nullptr);}
struct SurfaceArcPreparation::Data
{
    Rectangle source,region;Graph graph;bool rz=false;Result surface;
};
SurfaceArcPreparation::SurfaceArcPreparation(std::shared_ptr<const Data> data):data_(std::move(data)){}
SurfaceArcPreparation SurfaceArcPreparation::prepare(const Rectangle& source,const Graph& graph,const Rectangle& region,bool rz)
{
    auto surface=evaluateImpl(source,graph,region,rz,All,nullptr);
    return SurfaceArcPreparation(std::make_shared<Data>(Data{source,region,graph,rz,std::move(surface)}));
}
const Result& SurfaceArcPreparation::surface()const noexcept{return data_->surface;}
Result SurfaceArcPreparation::evaluate(const Rectangle& source,const Graph& graph,const Rectangle& region,bool rz,unsigned ownership)const
{
    // Preserve the original public physical/region guards before reusing events.
    const auto f=frame(source,graph,rz);(void)PlanarCellCut::capacity(region,rz);
    require(region.lower.x>=source.lower.x&&region.lower.y>=source.lower.y&&
        region.upper.x<=source.upper.x&&region.upper.y<=source.upper.y&&(ownership&~All)==0,
        "[QuadraticSurfaceMeasure] region outside original cell or invalid ownership");
    const auto sameBox=[](const Rectangle& a,const Rectangle& b){return a.lower.x==b.lower.x&&a.lower.y==b.lower.y&&a.upper.x==b.upper.x&&a.upper.y==b.upper.y;};
    const auto& a=graph.shape;const auto& b=data_->graph.shape;
    require(sameBox(source,data_->source)&&sameBox(region,data_->region)&&graph.exact==data_->graph.exact&&
        graph.offset==data_->graph.offset&&a.axis==b.axis&&a.sign==b.sign&&a.slope==b.slope&&a.second==b.second&&
        a.frameNormal.has_value()==b.frameNormal.has_value()&&(!a.frameNormal||
            (a.frameNormal->x==b.frameNormal->x&&a.frameNormal->y==b.frameNormal->y)),
        "[QuadraticSurfaceMeasure/ArcPreparation] different original graph owner, view or boxes");
    // Unencoded graphs have no exact polynomial owner; retain their old strict path.
    if(!graph.exact)return evaluateImpl(source,graph,region,rz,ownership,nullptr);
    // Ownership only suppresses an arc identically coincident with a side.
    // Test the exact coordinate polynomial, never endpoint distance or masks.
    for(int side=0;side<4;++side)if(!(ownership&(1u<<side))){const int axis=side/2;
        const auto border=axis?(side%2?region.upper.y:region.lower.y):(side%2?region.upper.x:region.lower.x);
        const auto p=Detail::coordinatePolynomial(*graph.exact,axis,border);
        if(p[0]==0&&p[1]==0&&p[2]==0)return {};
    }
    if(rz==data_->rz)return data_->surface;
    Result result;
    for(const auto& retained:data_->surface.patches){auto patch=retained;
        patch.integral=integral(f,patch.lower,patch.upper,rz,patch.begin.x,patch.end.x);
        result.integral.length+=patch.integral.length;result.integral.physicalMeasure+=patch.integral.physicalMeasure;
        result.integral.quadratureErrorEstimate+=patch.integral.quadratureErrorEstimate;result.patches.push_back(patch);
    }
    return result;
}
Sensitivity differentiate(const Rectangle& source,const Graph& graph,const Rectangle& region,bool rz,unsigned ownership)
{return differentiateImpl(source,graph,region,rz,ownership,nullptr);}
Sensitivity differentiateWithFixedCornerEndpoints(const Rectangle& source,const Graph& graph,const Rectangle& region,bool rz,unsigned ownership)
{return differentiateWithFixedCornerEndpointsImpl(source,graph,region,rz,ownership,nullptr);}
PreparedSurface prepareWithMaterialCurves(const Rectangle& source,const Graph& graph,const Rectangle& region,
    bool rz,bool withDerivative,bool fixedCornerEndpoints,unsigned ownership)
{return prepareWithTiming(source,graph,region,rz,withDerivative,fixedCornerEndpoints,true,ownership);}
PreparedSurface prepareWithTiming(const Rectangle& source,const Graph& graph,const Rectangle& region,
    bool rz,bool withDerivative,bool fixedCornerEndpoints,bool retainMaterialCurves,unsigned ownership)
{
    const auto start=std::chrono::steady_clock::now();PreparedSurface result;
    auto* curves=retainMaterialCurves?&result.materialCurves:nullptr;
    if(withDerivative)result.sensitivity=fixedCornerEndpoints?
        differentiateWithFixedCornerEndpointsImpl(source,graph,region,rz,ownership,curves,&result.timing):
        differentiateImpl(source,graph,region,rz,ownership,curves,&result.timing);
    else result.sensitivity.surface=evaluateImpl(source,graph,region,rz,ownership,curves,&result.timing);
    // Historical unencoded callers retain the old lifted-graph bridge.
    if(retainMaterialCurves&&!graph.exact){const auto curveStart=std::chrono::steady_clock::now();
        for(const auto&patch:result.sensitivity.surface.patches)
            result.materialCurves.push_back(QuadraticCellCut::materialCurveData(source,graph,patch.lower,patch.upper));
        result.timing.materialCurveSeconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-curveStart).count();}
    result.timing.totalSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    return result;
}
struct MomentPlan::Data {std::vector<MomentPatch> patches;bool rz=false;};
MomentPlan::MomentPlan(std::shared_ptr<const Data>data):data_(std::move(data)){}
const std::vector<MomentPatch>& MomentPlan::patches()const noexcept{return data_->patches;}
bool MomentPlan::axisymmetric()const noexcept{return data_->rz;}
MomentPlan MomentPlan::prepare(const Rectangle&source,const Graph&graph,const Rectangle&region,bool rz,unsigned ownership)
{
    require(graph.exact!=nullptr,"[QuadraticSurfaceMeasure/MomentPlan] exact graph required; legacy chart unsupported");
    Detail::validateGraph(graph);(void)QuadraticCellCut::basis(graph.shape);
    (void)PlanarCellCut::capacity(source,rz);(void)PlanarCellCut::capacity(region,rz);
    require(region.lower.x>=source.lower.x&&region.lower.y>=source.lower.y&&region.upper.x<=source.upper.x&&region.upper.y<=source.upper.y&&(ownership&~All)==0,
        "[QuadraticSurfaceMeasure/MomentPlan] invalid source/clip/ownership");
    namespace E=Detail::QuadraticExact;using R=E::Scalar;const auto&g=*graph.exact;
    auto data=std::make_shared<Data>();data->rz=rz;const auto metric=std::sqrt(E::rounded(g.determinant));
    const auto endpointDerivative=[&](const Detail::GraphEvent&event){const auto mask=event.boundary;
        require(mask!=0&&(mask&(mask-1))==0,"[QuadraticSurfaceMeasure/MomentPlan] corner/coincident endpoint derivative");
        const bool radial=(mask&(Left|Right))!=0;const R normal=radial?g.normal.x:g.normal.y,tangent=radial?g.tangent.x:g.tangent.y;
        const R t=event.parameter,denominator=tangent+normal*Detail::graphSlope(g,t);
        require(denominator!=0,"[QuadraticSurfaceMeasure/MomentPlan] tangential endpoint derivative singular");
        const std::array<R,3> basis{{1,t,t*t/2}};std::array<R,3> result{};for(int k=0;k<3;++k)result[k]=-normal*basis[k]/denominator;return result;};
    for(const auto&arc:Detail::graphArcs(g,region,ownership)){
        const R a=arc.first.parameter,b=arc.last.parameter,w=b-a,ma=Detail::graphSlope(g,a),mb=Detail::graphSlope(g,b);
        const auto da=endpointDerivative(arc.first),db=endpointDerivative(arc.last);const R qa=(g.q[2]*a+g.q[1])*a+g.q[0];
        const std::array<R,3> radius{{g.origin.x+g.tangent.x*a+g.normal.x*qa,w*(g.tangent.x+g.normal.x*ma),g.normal.x*g.q[2]*w*w}};
        MomentPatch p;p.lower=E::rounded(a);p.upper=E::rounded(b);p.begin=arc.first.position;p.end=arc.last.position;
        p.beginBoundary=arc.first.boundary;p.endBoundary=arc.last.boundary;p.m0=E::rounded(ma);p.m1=E::rounded(mb);p.scale=metric*E::rounded(w);
        require(p.upper>p.lower&&p.scale>0,"[QuadraticSurfaceMeasure/MomentPlan] unrepresentable exact patch");
        for(int j=0;j<3;++j)p.radius[j]=E::rounded(radius[j]);
        const std::array<R,3> height{{1,a,a*a/2}},slopeAtA{{0,1,a}},slopeAtB{{0,1,b}};
        for(int k=0;k<3;++k){const R dw=db[k]-da[k],dma=slopeAtA[k]+2*g.q[2]*da[k],dmb=slopeAtB[k]+2*g.q[2]*db[k];
            const R dqa=height[k]+ma*da[k];p.m0Derivative[k]=E::rounded(dma);p.m1Derivative[k]=E::rounded(dmb);p.scaleDerivative[k]=metric*E::rounded(dw);
            p.radiusDerivative[0][k]=E::rounded(g.tangent.x*da[k]+g.normal.x*dqa);
            p.radiusDerivative[1][k]=E::rounded(dw*(g.tangent.x+g.normal.x*ma)+w*g.normal.x*dma);
            p.radiusDerivative[2][k]=E::rounded(g.normal.x*((k==2?R(1)/2:R(0))*w*w+2*g.q[2]*w*dw));}
        data->patches.push_back(p);
    }return MomentPlan(std::move(data));
}
Sensitivity MomentPlan::apply(const std::vector<MomentValues>&values)const
{
    require(values.size()==data_->patches.size(),"[QuadraticSurfaceMeasure/MomentPlan] moment patch count mismatch");Sensitivity result;
    for(std::size_t i=0;i<values.size();++i){const auto&v=values[i];const auto&p=data_->patches[i];
        for(int j=0;j<3;++j)require(std::isfinite(v.moments[j])&&v.moments[j]>0&&std::isfinite(v.firstDerivative[j])&&std::isfinite(v.lastDerivative[j]),
            "[QuadraticSurfaceMeasure/MomentPlan] invalid same-model moment/partial");
        long double radial=0;for(int j=0;j<3;++j)radial+=p.radius[j]*v.moments[j];
        const auto length=p.scale*v.moments[0],area=data_->rz?CartesianMeasure::twoPi*p.scale*radial:length;
        require(std::isfinite(length)&&length>0&&std::isfinite(area)&&area>=0,"[QuadraticSurfaceMeasure/MomentPlan] invalid moment-model area");
        Patch patch{p.lower,p.upper,p.begin,p.end,p.beginBoundary,p.endBoundary,{length,area,0}};result.surface.patches.push_back(patch);
        result.surface.integral.length+=length;result.surface.integral.physicalMeasure+=area;
        for(int k=0;k<3;++k){long double derivative=0;
            if(!data_->rz)derivative=p.scaleDerivative[k]*v.moments[0]+p.scale*(v.firstDerivative[0]*p.m0Derivative[k]+v.lastDerivative[0]*p.m1Derivative[k]);
            else {long double dr=0;for(int j=0;j<3;++j)dr+=p.radiusDerivative[j][k]*v.moments[j]+p.radius[j]*(v.firstDerivative[j]*p.m0Derivative[k]+v.lastDerivative[j]*p.m1Derivative[k]);
                derivative=CartesianMeasure::twoPi*(p.scaleDerivative[k]*radial+p.scale*dr);}
            require(std::isfinite(derivative),"[QuadraticSurfaceMeasure/MomentPlan] nonfinite moment-model derivative");result.derivative[k]+=derivative;}
    }
    require(std::isfinite(result.surface.integral.length)&&std::isfinite(result.surface.integral.physicalMeasure),"[QuadraticSurfaceMeasure/MomentPlan] total moment-model area overflow");
    for(auto derivative:result.derivative)require(std::isfinite(derivative),"[QuadraticSurfaceMeasure/MomentPlan] total moment-model derivative overflow");
    return result;
}

}
