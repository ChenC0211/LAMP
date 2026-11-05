#include "Core/Geometry/PolynomialVelocitySlab.h"
#include "Core/Geometry/Detail/PolynomialSegmentRoots.h"
#include "Core/Geometry/Detail/ExactFaceTimeRational.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace LAMP::Geometry
{
namespace
{
namespace Arithmetic=Detail::PolynomialSegmentRoots::Detail;
using Interval=Arithmetic::Interval;
void require(bool ok,const char* reason)
{if(!ok)throw std::invalid_argument(std::string("[PolynomialVelocitySlab] ")+reason);}
void certify(bool ok,const char* reason)
{if(!ok)throw PolynomialRk2EnclosureUnresolved(reason);}
bool finite(Vector2D p){return std::isfinite(p.x) && std::isfinite(p.y);}
template<class N> struct Ops;
template<> struct Ops<double>
{
    static double point(double x){return x;}
    static double add(double a,double b){return a+b;}
    static double multiply(double a,double b){return a*b;}
    static double divide(double a,double b){return a/b;}
};
template<> struct Ops<Interval>
{
    static Interval point(double x){return {x,x};}
    static Interval add(Interval a,Interval b){return Arithmetic::add(a,b);}
    static Interval multiply(Interval a,Interval b){return Arithmetic::multiply(a,b);}
    static Interval divide(Interval a,double b){return Arithmetic::divide(a,b);}
};
template<class N> struct Jet {N value{},dx{},dy{};};
template<class N> Jet<N> component(const std::array<Vector2D,16>& coefficients,
    N x,N y,Vector2D scale,bool vertical)
{
    using O=Ops<N>;Jet<N> result;
    for(int i=3;i>=0;--i){N value{},dy{};
        for(int j=3;j>=0;--j){const auto c=coefficients[4*i+j];
            dy=O::add(O::multiply(dy,y),value);
            value=O::add(O::multiply(value,y),O::point(vertical?c.y:c.x));}
        result.dx=O::add(O::multiply(result.dx,x),result.value);
        result.value=O::add(O::multiply(result.value,x),value);
        result.dy=O::add(O::multiply(result.dy,x),dy);}
    result.dx=O::divide(result.dx,scale.x);result.dy=O::divide(result.dy,scale.y);return result;
}
template<class N> N normalized(N x,double origin,double scale)
{using O=Ops<N>;return O::divide(O::add(x,O::point(-origin)),scale);}
double mix(double a,double b,double theta)
{return theta==.5?(a+b)*.5:a*(1-theta)+b*theta;}
Interval hull(Interval a,Interval b)
{return {std::min(a.lo,b.lo),std::max(a.hi,b.hi)};}
Interval mix(Interval a,Interval b,Interval theta,bool previousEndpoint,bool arrivalEndpoint)
{
    const Interval complement=Arithmetic::add({1,1},{-theta.hi,-theta.lo});
    auto result=Arithmetic::add(Arithmetic::multiply(a,complement),Arithmetic::multiply(b,theta));
    // Include the actual evaluation's distinct arithmetic branches as well as
    // the real polynomial. Outward bounds must not rely on a sampled branch.
    if(theta.lo<=.5 && theta.hi>=.5)result=hull(result,Arithmetic::multiply(Arithmetic::add(a,b),{.5,.5}));
    if(previousEndpoint)result=hull(result,a);if(arrivalEndpoint)result=hull(result,b);
    return result;
}
PolynomialVelocityRange range(Interval a)
{require(std::isfinite(a.lo) && std::isfinite(a.hi) && a.lo<=a.hi,"nonfinite interval arithmetic");return {a.lo,a.hi};}
void validateTimes(double begin,double end,double oldTime,double newTime)
{require(std::isfinite(begin) && std::isfinite(end) && oldTime<=begin && begin<=end && end<=newTime,"time query outside slab");}
}

PolynomialRk2EnclosureUnresolved::PolynomialRk2EnclosureUnresolved(const std::string& reason)
    :std::invalid_argument("[PolynomialVelocitySlab] "+reason){}
PolynomialVelocitySlab::PolynomialVelocitySlab(double begin,double end,const PolynomialVelocityData& data)
    :begin_(begin),end_(end),data_(data)
{
    Detail::ExactExpansion::validateEnvironment();
    require(std::isfinite(begin) && std::isfinite(end) && end>begin && std::isfinite(end-begin),"invalid positive time slab");
    require(finite(data_.origin) && finite(data_.scale) && data_.scale.x>0 && data_.scale.y>0,"invalid coordinate basis");
    for(const auto* coefficients:{&data_.previous,&data_.arrival})
        for(const auto c:*coefficients)require(finite(c),"nonfinite velocity coefficient");
}
LevelSetBoundaryTrace::VelocityJet PolynomialVelocitySlab::evaluate(Vector2D point,double time) const
{
    Detail::ExactExpansion::validateEnvironment();validateTimes(time,time,begin_,end_);
    require(finite(point),"nonfinite spatial query");
    const double x=normalized(point.x,data_.origin.x,data_.scale.x),y=normalized(point.y,data_.origin.y,data_.scale.y);
    require(std::isfinite(x) && std::isfinite(y),"nonfinite normalized coordinate");
    const auto a=component(data_.previous,x,y,data_.scale,false),b=component(data_.previous,x,y,data_.scale,true);
    const auto c=component(data_.arrival,x,y,data_.scale,false),d=component(data_.arrival,x,y,data_.scale,true);
    const double theta=(time-begin_)/(end_-begin_);
    const auto blend=[&](double previous,double arrival){return time==begin_?previous:time==end_?arrival:mix(previous,arrival,theta);};
    const LevelSetBoundaryTrace::VelocityJet result{{blend(a.value,c.value),blend(b.value,d.value)},
        {blend(a.dx,c.dx),blend(a.dy,c.dy),blend(b.dx,d.dx),blend(b.dy,d.dy)}};
    require(finite(result.value) && LevelSetBoundaryTrace::Detail::finite(result.derivative),"nonfinite velocity evaluation");
    return result;
}
PolynomialVelocityBounds PolynomialVelocitySlab::enclose(Vector2D lower,Vector2D upper,double begin,double end) const
{
    Detail::ExactExpansion::validateEnvironment();validateTimes(begin,end,begin_,end_);
    require(finite(lower) && finite(upper) && lower.x<=upper.x && lower.y<=upper.y,"invalid closed spatial box");
    const auto x=normalized(Interval{lower.x,upper.x},data_.origin.x,data_.scale.x),
        y=normalized(Interval{lower.y,upper.y},data_.origin.y,data_.scale.y);
    const auto a=component(data_.previous,x,y,data_.scale,false),b=component(data_.previous,x,y,data_.scale,true);
    const auto c=component(data_.arrival,x,y,data_.scale,false),d=component(data_.arrival,x,y,data_.scale,true);
    const auto numerator=Arithmetic::add({begin,end},{-begin_,-begin_});
    const auto denominator=Arithmetic::add({end_,end_},{-begin_,-begin_});
    require(denominator.lo>0 && std::isfinite(denominator.hi),"unresolved time normalization");
    const Interval reciprocal{Arithmetic::down(1/denominator.hi),Arithmetic::up(1/denominator.lo)};
    const auto theta=Arithmetic::multiply(numerator,reciprocal);
    const auto blend=[&](Interval previous,Interval arrival){return range(mix(previous,arrival,theta,begin==begin_,end==end_));};
    return {{{blend(a.value,c.value),blend(b.value,d.value)}},
        {{blend(a.dx,c.dx),blend(a.dy,c.dy),blend(b.dx,d.dx),blend(b.dy,d.dy)}}};
}

namespace
{
// Uses the existing outward arithmetic, then intersects only with exact sign
// and monotonicity facts also obeyed by finite round-to-nearest operations.
struct Number
{
    Interval value{};
    Number()=default;
    Number(double x):value{x,x}{}
    explicit Number(Interval x):value(x){}
};
Number operator+(Number a,Number b)
{
    auto r=Arithmetic::add(a.value,b.value);
    if(a.value.lo>=0 && b.value.lo>=0)r.lo=std::max(0.,r.lo);
    if(a.value.hi<=0 && b.value.hi<=0)r.hi=std::min(0.,r.hi);
    return Number(r);
}
Number operator-(Number a){return Number(Interval{-a.value.hi,-a.value.lo});}
Number operator-(Number a,Number b)
{
    auto r=(a+-b).value;
    if(b.value.lo>=0)r.hi=std::min(r.hi,a.value.hi);
    if(b.value.hi<=0)r.lo=std::max(r.lo,a.value.lo);
    return Number(r);
}
Number operator*(Number a,Number b)
{
    auto r=Arithmetic::multiply(a.value,b.value);
    if((a.value.lo>=0 && b.value.lo>=0) || (a.value.hi<=0 && b.value.hi<=0))r.lo=std::max(0.,r.lo);
    if((a.value.lo>=0 && b.value.hi<=0) || (a.value.hi<=0 && b.value.lo>=0))r.hi=std::min(0.,r.hi);
    return Number(r);
}
template<> struct Ops<Number>
{
    static Number point(double x){return Number(x);}
    static Number add(Number a,Number b){return a+b;}
    static Number multiply(Number a,Number b){return a*b;}
    static Number divide(Number a,double b){
        auto r=Arithmetic::divide(a.value,b);
        if(a.value.lo>=0)r.lo=std::max(0.,r.lo);
        if(a.value.hi<=0)r.hi=std::min(0.,r.hi);
        return Number(r);
    }
};
template<> Number normalized(Number x,double origin,double scale)
{
    auto difference=(x-Number(origin)).value;
    // Rounded subtraction is monotone; equality is an exact zero, including
    // a nonzero chart origin and signed zeros. No near-zero threshold.
    if(x.value.lo>=origin)difference.lo=std::max(0.,difference.lo);
    if(x.value.hi<=origin)difference.hi=std::min(0.,difference.hi);
    return Ops<Number>::divide(Number(difference),scale);
}
Number mixSigned(Number a,Number b,Number theta,bool previousEndpoint,bool arrivalEndpoint)
{
    auto complement=(Number(1)-theta).value;
    complement.lo=std::max(0.,complement.lo);complement.hi=std::min(1.,complement.hi);
    auto result=(a*Number(complement)+b*theta).value;
    if(theta.value.lo<=.5 && theta.value.hi>=.5)result=hull(result,((a+b)*Number(.5)).value);
    if(previousEndpoint)result=hull(result,a.value);if(arrivalEndpoint)result=hull(result,b.value);
    return Number(result);
}
struct BoxVector {Number x,y;};
BoxVector operator-(BoxVector a,BoxVector b){return {a.x-b.x,a.y-b.y};}
BoxVector operator+(BoxVector a,BoxVector b){return {a.x+b.x,a.y+b.y};}
BoxVector operator*(BoxVector a,Number b){return {a.x*b,a.y*b};}
struct BoxMatrix {Number xx,xy,yx,yy;};
BoxMatrix operator+(BoxMatrix a,BoxMatrix b){return {a.xx+b.xx,a.xy+b.xy,a.yx+b.yx,a.yy+b.yy};}
BoxMatrix operator*(BoxMatrix a,Number b){return {a.xx*b,a.xy*b,a.yx*b,a.yy*b};}
Number boxed(PolynomialVelocityRange a){return Number(Interval{a.lower,a.upper});}
BoxVector values(const PolynomialVelocityBounds& a){return {boxed(a.value[0]),boxed(a.value[1])};}
BoxMatrix derivatives(const PolynomialVelocityBounds& a)
{return {boxed(a.derivative[0]),boxed(a.derivative[1]),boxed(a.derivative[2]),boxed(a.derivative[3])};}
Vector2D lower(BoxVector a){return {range(a.x.value).lower,range(a.y.value).lower};}
Vector2D upper(BoxVector a){return {range(a.x.value).upper,range(a.y.value).upper};}
bool inside(BoxVector a,const LevelSetBoundaryTrace::Domain& domain)
{const auto lo=lower(a),hi=upper(a);return lo.x>=domain.minCoord.x && lo.y>=domain.minCoord.y && hi.x<=domain.maxCoord.x && hi.y<=domain.maxCoord.y;}
bool exactZero(Number a){return a.value.lo==0 && a.value.hi==0;}
}
PolynomialVelocityBounds PolynomialVelocitySlab::encloseSignPreserving(Vector2D lo,Vector2D hi,double first,double last) const
{
    auto result=enclose(lo,hi,first,last); // original validation and outer certificate
    const auto x=normalized(Number(Interval{lo.x,hi.x}),data_.origin.x,data_.scale.x),
        y=normalized(Number(Interval{lo.y,hi.y}),data_.origin.y,data_.scale.y);
    const auto a=component(data_.previous,x,y,data_.scale,false),b=component(data_.previous,x,y,data_.scale,true),
        c=component(data_.arrival,x,y,data_.scale,false),d=component(data_.arrival,x,y,data_.scale,true);
    const auto numerator=Arithmetic::add({first,last},{-begin_,-begin_}),
        denominator=Arithmetic::add({end_,end_},{-begin_,-begin_});
    const Interval reciprocal{Arithmetic::down(1/denominator.hi),Arithmetic::up(1/denominator.lo)};
    auto theta=Arithmetic::multiply(numerator,reciprocal);
    // Both the real ratio and the actual rounded evaluation are in [0,1]
    // on this slab, by monotonicity of the identical endpoint expressions.
    theta.lo=std::max(0.,theta.lo);theta.hi=std::min(1.,theta.hi);
    const auto blend=[&](PolynomialVelocityRange outer,Number previous,Number arrival){
        const auto inner=mixSigned(previous,arrival,Number(theta),first==begin_,last==end_).value;
        return range({std::max(outer.lower,inner.lo),std::min(outer.upper,inner.hi)});
    };
    result.value={blend(result.value[0],a.value,c.value),blend(result.value[1],b.value,d.value)};
    result.derivative={blend(result.derivative[0],a.dx,c.dx),blend(result.derivative[1],a.dy,c.dy),
        blend(result.derivative[2],b.dx,d.dx),blend(result.derivative[3],b.dy,d.dy)};
    return result;
}
PolynomialReachableSourceBox PolynomialRk2Reachability::cover(Vector2D lo,Vector2D hi) const
{
    Detail::ExactExpansion::validateEnvironment();
    require(finite(lo) && finite(hi) && lo.x<=hi.x && lo.y<=hi.y
        && lo.x>=lower_.x && lo.y>=lower_.y && hi.x<=upper_.x && hi.y<=upper_.y,
        "reachability input outside closed physical domain");
    const auto extent=[](double a,double b,double radius,double wall0,double wall1){
        if(radius==0)return Interval{a,b}; // exact unchanged coordinate
        return Interval{std::max(wall0,(Number(a)-Number(radius)).value.lo),
            std::min(wall1,(Number(b)+Number(radius)).value.hi)};
    };
    const auto x=extent(lo.x,hi.x,radius_.x,lower_.x,upper_.x),y=extent(lo.y,hi.y,radius_.y,lower_.y,upper_.y);
    return {{x.lo,y.lo},{x.hi,y.hi},{{radius_.x>0 && x.lo==lower_.x,radius_.y>0 && y.lo==lower_.y,
        radius_.x>0 && x.hi==upper_.x,radius_.y>0 && y.hi==upper_.y}}};
}
PolynomialRk2Reachability PolynomialVelocitySlab::encloseRk2Reachability(
    const LevelSetBoundaryTrace::Domain& domain,double last) const
{
    Detail::ExactExpansion::validateEnvironment();validateTimes(begin_,last,begin_,end_);
    require(last>begin_,"reachability requires positive duration");
    const double duration=last-begin_;
    const LevelSetBoundaryTrace::DifferentiableVelocity current=[this,last](Vector2D p){return evaluate(p,last);},
        previous=[this](Vector2D p){return evaluate(p,begin_);};
    const LevelSetBoundaryTrace::Detail::Tracer<true> validation(domain,duration,current,previous);(void)validation;
    require(!domain.activeSegment,"reachability cannot certify opaque active coverage");
    certify(begin_==0 || (begin_>0 && begin_<=std::numeric_limits<double>::max()*.5 && last<=2*begin_),
        "reachability cannot certify physical boundary history time");
    auto speed=encloseTightTime(domain.minCoord,domain.maxCoord,begin_,last);
    const Number epsilon(std::numeric_limits<double>::epsilon()),eta(std::numeric_limits<double>::denorm_min());
    const auto nativeMix=[&](PolynomialVelocityRange r){
        if(r.lower==0 && r.upper==0)return r;
        const Number magnitude(std::max(std::abs(r.lower),std::abs(r.upper)));
        const auto error=(Number(8)*(epsilon*magnitude+eta)).value.hi;
        auto low=(Number(r.lower)-Number(error)).value.lo,high=(Number(r.upper)+Number(error)).value.hi;
        if(r.lower>=0)low=std::max(0.,low);if(r.upper<=0)high=std::min(0.,high);
        // Also ensures the actual native (current+previous)*.5 branch and
        // subsequent event coordinate differences cannot overflow.
        certify(std::max(std::abs(low),std::abs(high))<=std::numeric_limits<double>::max()/16,
            "reachability unresolved intermediate magnitude");
        return range({low,high});
    };
    for(auto& r:speed.value)r=nativeMix(r);for(auto& r:speed.derivative)r=nativeMix(r);
    const double ux=std::max(std::abs(speed.value[0].lower),std::abs(speed.value[0].upper)),
        uy=std::max(std::abs(speed.value[1].lower),std::abs(speed.value[1].upper));
    const double limit=.5*std::min(domain.cellSize.x,domain.cellSize.y);
    const auto travel=((Number(ux)+Number(uy))*Number(duration)).value.hi;
    certify(std::isfinite(travel) && travel<limit && std::hypot(ux,uy)*duration<=limit,
        "reachability cannot certify original half-cell speed guard");
    certify(!domain.radialAxisSymmetry || (ux==0 && speed.derivative[1].lower==0 && speed.derivative[1].upper==0),
        "reachability requires exact zero radial velocity/axis derivative");
    const auto radius=[&](double u,double a,double b){
        if(u==0)return 0.;
        const Number U(u),D(duration),A(std::max(std::abs(a),std::abs(b)));
        const Number magnitude=A+U*D;
        certify(magnitude.value.hi<=std::numeric_limits<double>::max()/16,"reachability unresolved coordinate magnitude");
        // N=128 native steps: time sum error <=2(N+2)(eps*D+eta).
        // Each target/midpoint <=4 eps*(A+UD)+4 eta*(1+D), plus
        // 16 such units for final bracket interpolation. Bisection trials
        // restart at the same point and do not add trajectory increments.
        const Number time=D+Number(260)*(epsilon*D+eta);
        const Number error=Number(532)*(epsilon*magnitude+eta*(Number(1)+D));
        return range((U*time+error).value).upper;
    };
    return {domain.minCoord,domain.maxCoord,{radius(ux,domain.minCoord.x,domain.maxCoord.x),
        radius(uy,domain.minCoord.y,domain.maxCoord.y)}};
}
PolynomialRk2Bounds PolynomialVelocitySlab::encloseInteriorRk2(const LevelSetBoundaryTrace::Domain& domain,
    Vector2D lo,Vector2D hi,double arrivalBegin,double arrivalEnd) const
{return encloseInteriorRk2Impl(domain,lo,hi,arrivalBegin,arrivalEnd,false);}
PolynomialRk2Bounds PolynomialVelocitySlab::encloseInteriorRk2TightTime(const LevelSetBoundaryTrace::Domain& domain,
    Vector2D lo,Vector2D hi,double arrivalBegin,double arrivalEnd) const
{return encloseInteriorRk2Impl(domain,lo,hi,arrivalBegin,arrivalEnd,true);}
PolynomialVelocityBounds PolynomialVelocitySlab::encloseTightTime(Vector2D a,Vector2D b,double first,double last) const
{
    auto result=encloseSignPreserving(a,b,first,last);
    const auto old=encloseSignPreserving(a,b,begin_,begin_),end=encloseSignPreserving(a,b,end_,end_);
    const auto tighten=[](PolynomialVelocityRange original,PolynomialVelocityRange previous,PolynomialVelocityRange arrival){
        const double low=std::min(previous.lower,arrival.lower),high=std::max(previous.upper,arrival.upper);
        if(low==0 && high==0)return PolynomialVelocityRange{0,0};
        // Actual computed theta is in [0,1] by monotonicity of the
        // identical positive endpoint subtraction/division. The real
        // convex combination is in [low,high]. Rounded (1-theta), two
        // products and one sum cost < 8*epsilon*maxabs + 8*denorm_min.
        // This also covers (a+b)*.5 and direct endpoint branches. The
        // original interval already rejects nonfinite intermediates.
        const double magnitude=std::max(std::abs(low),std::abs(high));
        const double relative=8*std::numeric_limits<double>::epsilon(),absolute=8*std::numeric_limits<double>::denorm_min();
        const auto error=Arithmetic::add(Arithmetic::multiply({relative,relative},{magnitude,magnitude}),{absolute,absolute});
        const auto lower=Arithmetic::add({low,low},{-error.hi,-error.hi}),upper=Arithmetic::add({high,high},{error.hi,error.hi});
        return range({std::max(original.lower,lower.lo),std::min(original.upper,upper.hi)});
    };
    for(int k=0;k<2;++k)result.value[k]=tighten(result.value[k],old.value[k],end.value[k]);
    for(int k=0;k<4;++k)result.derivative[k]=tighten(result.derivative[k],old.derivative[k],end.derivative[k]);
    return result;
}
Boundary PolynomialVelocitySlab::encloseInitialBoundaryHit(const LevelSetBoundaryTrace::Domain& domain,
    Vector2D lo,Vector2D hi,double first,double last) const
{
    Detail::ExactExpansion::validateEnvironment();validateTimes(first,last,begin_,end_);
    require(last>begin_ && finite(lo) && finite(hi) && lo.x<=hi.x && lo.y<=hi.y,
        "initial boundary enclosure requires a closed box and positive duration");
    const double duration=last-begin_;
    const LevelSetBoundaryTrace::DifferentiableVelocity current=[this,last](Vector2D p){return evaluate(p,last);},
        previous=[this](Vector2D p){return evaluate(p,begin_);};
    const LevelSetBoundaryTrace::Detail::Tracer<true> validation(domain,duration,current,previous);(void)validation;
    require(!domain.activeSegment,"initial boundary enclosure cannot certify opaque active coverage");
    const BoxVector input{Number(Interval{lo.x,hi.x}),Number(Interval{lo.y,hi.y})};
    require(inside(input,domain),"initial boundary enclosure input outside physical domain");
    const auto speed=encloseTightTime(lo,hi,first,last);
    const double x=std::max(std::abs(speed.value[0].lower),std::abs(speed.value[0].upper)),
        y=std::max(std::abs(speed.value[1].lower),std::abs(speed.value[1].upper));
    const double travel=((Number(x)+Number(y))*Number(duration)).value.hi,
        limit=.5*std::min(domain.cellSize.x,domain.cellSize.y);
    certify(std::isfinite(travel) && travel<limit && std::hypot(x,y)*duration<=limit,
        "initial boundary enclosure cannot certify original half-cell speed guard");
    if(domain.radialAxisSymmetry && lo.x==0){
        const auto axis=encloseTightTime({0,lo.y},{0,hi.y},first,last);
        certify(axis.value[0].lower==0 && axis.value[0].upper==0
            && axis.derivative[1].lower==0 && axis.derivative[1].upper==0,
            "initial boundary enclosure cannot certify exact radial axis velocity/derivative");
    }
    // Mirrors outwardBacktraceSide: one side must be certain for EVERY
    // point, and every other side must be impossible, including corner points.
    const std::array<Boundary,4> sides{Boundary::west,Boundary::south,Boundary::east,Boundary::north};
    const std::array<bool,4> certain{
        lo.x==domain.minCoord.x && hi.x==lo.x && speed.value[0].lower>0,
        lo.y==domain.minCoord.y && hi.y==lo.y && speed.value[1].lower>0,
        hi.x==domain.maxCoord.x && lo.x==hi.x && speed.value[0].upper<0,
        hi.y==domain.maxCoord.y && lo.y==hi.y && speed.value[1].upper<0};
    const std::array<bool,4> possible{
        lo.x==domain.minCoord.x && speed.value[0].upper>0,
        lo.y==domain.minCoord.y && speed.value[1].upper>0,
        hi.x==domain.maxCoord.x && speed.value[0].lower<0,
        hi.y==domain.maxCoord.y && speed.value[1].lower<0};
    for(int k=0;k<4;++k)if(certain[k]){
        bool unique=true;for(int j=0;j<4;++j)if(j!=k && possible[j])unique=false;
        if(unique)return sides[k];
    }
    throw PolynomialRk2EnclosureUnresolved("initial boundary enclosure cannot certify a unique immediate hit");
}
PolynomialRk2Bounds PolynomialVelocitySlab::encloseInteriorRk2Impl(const LevelSetBoundaryTrace::Domain& domain,
    Vector2D lo,Vector2D hi,double arrivalBegin,double arrivalEnd,bool tightTime) const
{
    Detail::ExactExpansion::validateEnvironment();validateTimes(arrivalBegin,arrivalEnd,begin_,end_);
    require(finite(lo) && finite(hi) && lo.x<=hi.x && lo.y<=hi.y,"RK2 enclosure invalid closed spatial box");
    const double fullDuration=arrivalEnd-begin_;
    const LevelSetBoundaryTrace::DifferentiableVelocity current=[this,arrivalEnd](Vector2D p){return evaluate(p,arrivalEnd);},
        previous=[this](Vector2D p){return evaluate(p,begin_);};
    // Reuse the actual tracer's domain/time/coordinate-scale validation. No
    // trajectory is run or substituted for the subsequent range propagation.
    const LevelSetBoundaryTrace::Detail::Tracer<true> validation(domain,fullDuration,current,previous);(void)validation;
    require(!domain.activeSegment,"RK2 enclosure cannot certify opaque active coverage");
    const auto query=[&](Vector2D a,Vector2D b,double first,double last){
        return tightTime?encloseTightTime(a,b,first,last):enclose(a,b,first,last);
    };
    const BoxVector input{Number(Interval{lo.x,hi.x}),Number(Interval{lo.y,hi.y})};
    require(inside(input,domain),"RK2 enclosure input outside physical domain");
    if(arrivalEnd==begin_)return {lo,hi,lo,hi,lo,hi,{{{1,1},{0,0},{0,0},{1,1}}},0};
    auto d=(Number(Interval{arrivalBegin,arrivalEnd})-Number(begin_)).value;
    d.lo=std::max(0.,d.lo);const Number duration(d); // ordered times imply d>=0, including rounded subtraction
    double maximumTravel=0;const double limit=.5*std::min(domain.cellSize.x,domain.cellSize.y);
    const auto speedCheck=[&](const PolynomialVelocityBounds& velocity){
        const double x=std::max(std::abs(velocity.value[0].lower),std::abs(velocity.value[0].upper)),
            y=std::max(std::abs(velocity.value[1].lower),std::abs(velocity.value[1].upper));
        const double travel=((Number(x)+Number(y))*duration).value.hi;
        const double nativeTravel=std::hypot(x,y)*fullDuration;
        certify(std::isfinite(travel) && travel<limit && std::isfinite(nativeTravel) && nativeTravel<=limit,
            "RK2 enclosure cannot certify original half-cell speed guard");
        maximumTravel=std::max(maximumTravel,travel);
    };
    const auto axisCheck=[&](BoxVector points,double a,double b){
        if(!domain.radialAxisSymmetry || points.x.value.lo!=0)return;
        const auto axial=query({0,points.y.value.lo},{0,points.y.value.hi},a,b);
        certify(exactZero(boxed(axial.value[0])) && exactZero(boxed(axial.derivative[1])),
            "RK2 enclosure cannot certify exact radial axis velocity/derivative");
    };
    const auto start=query(lo,hi,arrivalBegin,arrivalEnd);speedCheck(start);axisCheck(input,arrivalBegin,arrivalEnd);
    certify((lo.x!=domain.minCoord.x || start.value[0].upper<=0)
        && (hi.x!=domain.maxCoord.x || start.value[0].lower>=0)
        && (lo.y!=domain.minCoord.y || start.value[1].upper<=0)
        && (hi.y!=domain.maxCoord.y || start.value[1].lower>=0),"RK2 enclosure cannot exclude initial boundary hit");
    const auto half=LevelSetBoundaryTrace::Detail::rk2Midpoint(input,values(start),duration);
    certify(inside(half,domain),"RK2 enclosure midpoint domain unresolved");
    const auto halfLower=lower(half),halfUpper=upper(half);
    const auto newMiddle=query(halfLower,halfUpper,arrivalBegin,arrivalEnd),oldMiddle=query(halfLower,halfUpper,begin_,begin_);
    speedCheck(newMiddle);speedCheck(oldMiddle);axisCheck(half,arrivalBegin,arrivalEnd);axisCheck(half,begin_,begin_);
    auto middle=(values(newMiddle)+values(oldMiddle))*Number(.5);
    auto middleDerivative=(derivatives(newMiddle)+derivatives(oldMiddle))*Number(.5);
    // If half of a positive subnormal duration rounds to zero, the actual
    // tracer selects its elapsed==0 endpoint branch, before the mean branch.
    // Include that branch (also harmless for the zero-duration identity).
    if(d.lo<=std::numeric_limits<double>::denorm_min()){
        const auto unite=[](Number a,Number b){return Number(hull(a.value,b.value));};
        const auto v=values(newMiddle);const auto j=derivatives(newMiddle);
        middle={unite(middle.x,v.x),unite(middle.y,v.y)};
        middleDerivative={unite(middleDerivative.xx,j.xx),unite(middleDerivative.xy,j.xy),
            unite(middleDerivative.yx,j.yx),unite(middleDerivative.yy,j.yy)};
    }
    PolynomialVelocityBounds middleVelocity{};middleVelocity.value={range(middle.x.value),range(middle.y.value)};speedCheck(middleVelocity);
    const auto target=LevelSetBoundaryTrace::Detail::rk2Target(input,middle,duration);
    certify(inside(target,domain),"RK2 enclosure target domain unresolved");
    const auto noNewWall=[](Number from,Number to,Number speed,double minimum,double maximum){
        if(exactZero(speed))return true; // target equals input, including finite rounded arithmetic
        // Original target = from - speed*duration, duration>=0. Rounded
        // multiplication/subtraction preserve this weak coordinate order,
        // including signed zero and underflow. A point strictly inside cannot
        // newly touch the wall in the opposite direction; a point already on
        // that wall is not a touchesNewWall event. Keep the strict spatial
        // proof for directions whose sign has not been certified.
        return (speed.value.hi<=0 || from.value.hi<=minimum || to.value.lo>minimum)
            && (speed.value.lo>=0 || from.value.lo>=maximum || to.value.hi<maximum);
    };
    certify(noNewWall(input.x,target.x,middle.x,domain.minCoord.x,domain.maxCoord.x)
        && noNewWall(input.y,target.y,middle.y,domain.minCoord.y,domain.maxCoord.y),"RK2 enclosure new wall contact unresolved");
    const auto jacobian=LevelSetBoundaryTrace::Detail::rk2Jacobian(duration,derivatives(start),middleDerivative);
    return {lo,hi,halfLower,halfUpper,lower(target),upper(target),
        {{range(jacobian.xx.value),range(jacobian.xy.value),range(jacobian.yx.value),range(jacobian.yy.value)}},maximumTravel};
}

namespace
{
Interval interval(PolynomialVelocityRange a){return {a.lower,a.upper};}
Interval nonnegative(Interval a)
{a.lo=std::max(0.,a.lo);range(a);return a;}
void validIntegral(const PolynomialFaceIntegralBounds& b)
{
    for(const auto v:b.volumes)require(std::isfinite(v.lower) && std::isfinite(v.upper) && v.lower>=0 && v.lower<=v.upper,"invalid gross integral interval");
    for(int k=0;k<2;++k)require(b.volumes[k+2].lower<=b.volumes[k].lower && b.volumes[k+2].upper<=b.volumes[k].upper,"gross liquid/total interval ordering");
}
std::array<Interval,5> powers(Interval a)
{std::array<Interval,5> p{};p[0]={1,1};for(int i=1;i<5;++i)p[i]=Arithmetic::multiply(p[i-1],a);return p;}
std::array<Interval,5> averagePowers(Interval a,Interval b)
{
    const auto pa=powers(a),pb=powers(b);std::array<Interval,5> moments{};moments[0]={1,1};
    for(int j=1;j<5;++j){for(int k=0;k<=j;++k)moments[j]=Arithmetic::add(moments[j],Arithmetic::multiply(pa[j-k],pb[k]));
        moments[j]=Arithmetic::divide(moments[j],double(j+1));}return moments;
}
Interval differenceInterval(double a,double b)
{return Arithmetic::add({a,a},{-b,-b});}
Interval contentInterval(long double lo,long double hi)
{
    require(std::isfinite(lo) && std::isfinite(hi) && lo>=0 && lo<=hi && hi<=1,"invalid bounded content interval");
    // Preserve exact endpoints/pure weights; other long-double weights may
    // need an outward conversion to this binary64 arithmetic engine.
    const double a=static_cast<double>(lo),b=static_cast<double>(hi);
    return {lo==0?0:lo==1?1:std::max(0.,Arithmetic::down(a)),hi==1?1:hi==0?0:std::min(1.,Arithmetic::up(b))};
}
}
PolynomialFaceIntegralBounds PolynomialFaceIntegralBounds::added(const PolynomialFaceIntegralBounds& other) const
{
    Detail::ExactExpansion::validateEnvironment();validIntegral(*this);validIntegral(other);PolynomialFaceIntegralBounds result;
    for(int i=0;i<4;++i)result.volumes[i]=range(nonnegative(Arithmetic::add(interval(volumes[i]),interval(other.volumes[i]))));
    validIntegral(result);return result;
}
PolynomialFaceIntegralBounds PolynomialFaceIntegralBounds::intersected(const PolynomialFaceIntegralBounds& other) const
{
    Detail::ExactExpansion::validateEnvironment();validIntegral(*this);validIntegral(other);PolynomialFaceIntegralBounds result;
    for(int i=0;i<4;++i)result.volumes[i]={std::max(volumes[i].lower,other.volumes[i].lower),std::min(volumes[i].upper,other.volumes[i].upper)};
    validIntegral(result);return result;
}
double PolynomialFaceIntegralBounds::errorBound() const
{
    Detail::ExactExpansion::validateEnvironment();validIntegral(*this);Interval total{};
    for(const auto v:volumes)if(v.upper!=v.lower){const auto width=Arithmetic::up(v.upper-v.lower);total=Arithmetic::add(total,{width,width});}
    return range(total).upper;
}
PolynomialFaceIntegralBounds PolynomialFaceIntegralBounds::withContentRange(long double lo,long double hi) const
{
    Detail::ExactExpansion::validateEnvironment();validIntegral(*this);const auto rho=contentInterval(lo,hi);auto result=*this;
    for(int k=0;k<2;++k){const auto total=interval(volumes[k]);auto liquid=nonnegative(Arithmetic::multiply(total,rho));
        liquid.hi=std::min(liquid.hi,total.hi);result.volumes[k+2]=range(liquid);}
    validIntegral(result);return result;
}
std::array<double,4> PolynomialFaceIntegralBounds::representative() const
{
    Detail::ExactExpansion::validateEnvironment();validIntegral(*this);std::array<double,4> values{};
    for(int i=0;i<4;++i){const auto b=volumes[i];
        // Coordinatewise monotone in both endpoints, hence preserving the
        // liquid/total ordering. Endpoint selection only handles underflow.
        values[i]=std::max(b.lower,std::min(b.upper,b.lower*.5+b.upper*.5));}
    require(values[2]<=values[0] && values[3]<=values[1],"gross representative outside total subset");return values;
}
std::optional<double> polynomialFaceMidpoint(double lo,double hi)
{
    Detail::ExactExpansion::validateEnvironment();require(std::isfinite(lo) && std::isfinite(hi) && lo<hi && std::isfinite(hi-lo),"invalid face subdivision extent");
    const double middle=lo+(hi-lo)*.5;if(middle<=lo || middle>=hi)return std::nullopt;return middle;
}
PolynomialFaceIntegralBounds PolynomialVelocitySlab::encloseNormalFaceIntegral(Vector2D a,Vector2D b,
    double begin,double end,bool radial,long double contentLower,long double contentUpper) const
{
    Detail::ExactExpansion::validateEnvironment();validateTimes(begin,end,begin_,end_);
    require(begin<end && finite(a) && finite(b),"invalid positive face-time extent");
    const bool horizontal=a.y==b.y && a.x<b.x,vertical=a.x==b.x && a.y<b.y;
    require(horizontal || vertical,"face integral requires an ascending axis-aligned segment");
    require(!radial || (a.x>=0 && b.x>=0),"negative radius in face integral");
    contentInterval(contentLower,contentUpper);
    bool zeroNormal=true;for(const auto* coefficients:{&data_.previous,&data_.arrival})
        for(const auto c:*coefficients)zeroNormal=zeroNormal && (horizontal?c.y:c.x)==0;
    if(zeroNormal)return {};
    const auto xa=normalized(Interval{a.x,a.x},data_.origin.x,data_.scale.x),xb=normalized(Interval{b.x,b.x},data_.origin.x,data_.scale.x),
        ya=normalized(Interval{a.y,a.y},data_.origin.y,data_.scale.y),yb=normalized(Interval{b.y,b.y},data_.origin.y,data_.scale.y);
    const auto px=powers(xa),py=powers(ya),mx=averagePowers(xa,xb),my=averagePowers(ya,yb);
    const auto spatial=[&](const std::array<Vector2D,16>& coefficients){Interval sum{};
        for(int i=0;i<4;++i)for(int j=0;j<4;++j){const auto c=coefficients[4*i+j];const double value=horizontal?c.y:c.x;if(value==0)continue;
            auto x=horizontal?mx[i]:px[i];
            if(horizontal && radial)x=Arithmetic::add(Arithmetic::multiply({data_.origin.x,data_.origin.x},mx[i]),
                Arithmetic::multiply({data_.scale.x,data_.scale.x},mx[i+1]));
            const auto y=horizontal?py[j]:my[j];sum=Arithmetic::add(sum,Arithmetic::multiply({value,value},Arithmetic::multiply(x,y)));}return sum;};
    const auto denominator=differenceInterval(end_,begin_);require(denominator.lo>0,"unresolved face time normalization");
    const Interval reciprocal{Arithmetic::down(1/denominator.hi),Arithmetic::up(1/denominator.lo)};
    auto theta=Arithmetic::multiply(Arithmetic::add(differenceInterval(begin,begin_),differenceInterval(end,begin_)),{.5,.5});
    theta=Arithmetic::multiply(theta,reciprocal);theta.lo=std::max(0.,theta.lo);theta.hi=std::min(1.,theta.hi);
    auto complement=Arithmetic::add({1,1},{-theta.hi,-theta.lo});complement.lo=std::max(0.,complement.lo);complement.hi=std::min(1.,complement.hi);
    const auto mean=Arithmetic::add(Arithmetic::multiply(spatial(data_.previous),complement),
        Arithmetic::multiply(spatial(data_.arrival),theta));
    const auto length=differenceInterval(horizontal?b.x:b.y,horizontal?a.x:a.y),duration=differenceInterval(end,begin);
    auto factor=Arithmetic::multiply(length,duration),measure=factor;
    if(radial){constexpr double pi=3.141592653589793238462643383279502884;
        const auto twoPi=Arithmetic::multiply({2,2},{Arithmetic::down(pi),Arithmetic::up(pi)});
        factor=Arithmetic::multiply(factor,twoPi);if(vertical)factor=Arithmetic::multiply(factor,{a.x,a.x});
        const auto meanRadius=Arithmetic::multiply(Arithmetic::add({a.x,a.x},{b.x,b.x}),{.5,.5});
        measure=Arithmetic::multiply(Arithmetic::multiply(measure,twoPi),nonnegative(meanRadius));}
    const auto net=Arithmetic::multiply(mean,factor);range(net);measure=nonnegative(measure);
    const auto velocity=enclose(a,b,begin,end).value[horizontal?1:0];Interval positive{},negative{};
    if(velocity.lower>=0)positive=nonnegative(net);
    else if(velocity.upper<=0)negative=nonnegative({-net.hi,-net.lo});
    else{
        positive={std::max(0.,net.lo),nonnegative(Arithmetic::multiply(measure,{0,velocity.upper})).hi};
        negative={std::max(0.,-net.hi),nonnegative(Arithmetic::multiply(measure,{0,-velocity.lower})).hi};
        positive.hi=std::min(positive.hi,Arithmetic::add(negative,net).hi);
        negative.hi=std::min(negative.hi,Arithmetic::add(positive,{-net.hi,-net.lo}).hi);
    }
    PolynomialFaceIntegralBounds result;result.volumes[0]=range(positive);result.volumes[1]=range(negative);
    for(int k=0;k<2;++k)result.volumes[k+2]={0,result.volumes[k].upper};
    return result.withContentRange(contentLower,contentUpper);
}

namespace
{
struct Affine
{
    Interval center{},tangent{},time{},remainder{};
    Affine()=default;
    Affine(double point):center{point,point}{}
};
bool exact(Interval a,double point){return a.lo==point && a.hi==point;}
bool exact(const Affine& a,double point)
{return exact(a.center,point) && exact(a.tangent,0) && exact(a.time,0) && exact(a.remainder,0);}
Interval negative(Interval a){return {-a.hi,-a.lo};}
double affineMagnitude(Interval a){range(a);return std::max(std::abs(a.lo),std::abs(a.hi));}
Interval variableRange(const Affine& a)
{return Arithmetic::add(Arithmetic::add(Arithmetic::multiply(a.tangent,{-1,1}),Arithmetic::multiply(a.time,{-1,1})),a.remainder);}
Interval affineRange(const Affine& a)
{const auto result=Arithmetic::add(a.center,variableRange(a));range(result);return result;}
void addSymmetric(Affine& a,Interval radius)
{range(radius);require(radius.hi>=0,"negative affine rounding radius");a.remainder=Arithmetic::add(a.remainder,{-radius.hi,radius.hi});}
void roundOperation(Affine& a,Interval value)
{
    const double m=affineMagnitude(value);if(m==0)return;
    constexpr double eps=std::numeric_limits<double>::epsilon(),tiny=std::numeric_limits<double>::denorm_min();
    addSymmetric(a,Arithmetic::add(Arithmetic::multiply({eps,eps},{m,m}),{tiny,tiny}));
}
Affine operator+(const Affine& a,const Affine& b)
{
    if(exact(a,0))return b;if(exact(b,0))return a;
    Affine out;out.center=Arithmetic::add(a.center,b.center);out.tangent=Arithmetic::add(a.tangent,b.tangent);
    out.time=Arithmetic::add(a.time,b.time);out.remainder=Arithmetic::add(a.remainder,b.remainder);
    roundOperation(out,Arithmetic::add(affineRange(a),affineRange(b)));return out;
}
Affine operator-(const Affine& a)
{Affine out;out.center=negative(a.center);out.tangent=negative(a.tangent);out.time=negative(a.time);out.remainder=negative(a.remainder);return out;}
Affine operator-(const Affine& a,const Affine& b){return a+-b;}
Affine operator*(const Affine& a,const Affine& b)
{
    if(exact(a,0) || exact(b,0))return {};if(exact(a,1))return b;if(exact(b,1))return a;
    Affine out;out.center=Arithmetic::multiply(a.center,b.center);
    out.tangent=Arithmetic::add(Arithmetic::multiply(a.center,b.tangent),Arithmetic::multiply(b.center,a.tangent));
    out.time=Arithmetic::add(Arithmetic::multiply(a.center,b.time),Arithmetic::multiply(b.center,a.time));
    out.remainder=Arithmetic::add(Arithmetic::add(Arithmetic::multiply(a.center,b.remainder),Arithmetic::multiply(b.center,a.remainder)),
        Arithmetic::multiply(variableRange(a),variableRange(b)));
    roundOperation(out,Arithmetic::multiply(affineRange(a),affineRange(b)));return out;
}
Affine divideAffine(const Affine& a,double denominator)
{
    require(std::isfinite(denominator) && denominator>0,"invalid affine divisor");if(exact(a,0) || denominator==1)return a;
    Affine out;out.center=Arithmetic::divide(a.center,denominator);out.tangent=Arithmetic::divide(a.tangent,denominator);
    out.time=Arithmetic::divide(a.time,denominator);out.remainder=Arithmetic::divide(a.remainder,denominator);
    roundOperation(out,Arithmetic::divide(affineRange(a),denominator));return out;
}
template<> struct Ops<Affine>
{
    static Affine point(double x){return Affine(x);}
    static Affine add(const Affine& a,const Affine& b){return a+b;}
    static Affine multiply(const Affine& a,const Affine& b){return a*b;}
    static Affine divide(const Affine& a,double b){return divideAffine(a,b);}
};
struct AffineVector {Affine x,y;};
AffineVector operator+(const AffineVector& a,const AffineVector& b){return {a.x+b.x,a.y+b.y};}
AffineVector operator-(const AffineVector& a,const AffineVector& b){return {a.x-b.x,a.y-b.y};}
AffineVector operator*(const AffineVector& a,const Affine& b){return {a.x*b,a.y*b};}
Affine affineParameter(double lo,double hi,bool temporal)
{
    Affine a;a.center=Arithmetic::multiply(Arithmetic::add({lo,lo},{hi,hi}),{.5,.5});
    const auto radius=Arithmetic::multiply(differenceInterval(hi,lo),{.5,.5});
    if(temporal)a.time=radius;else a.tangent=radius;affineRange(a);return a;
}
Affine mixAffine(const Affine& a,const Affine& b,const Affine& theta)
{
    if(exact(a,0) && exact(b,0))return {};
    if(exact(theta,0))return a;if(exact(theta,1))return b;
    auto out=a*(Affine(1)-theta)+b*theta;const auto t=affineRange(theta);
    if(t.lo<=.5 && t.hi>=.5){
        // Contains the real convex combination, so only the alternative
        // original (a+b)*.5 branch's RN error must additionally be covered.
        const double m=std::max(affineMagnitude(affineRange(a)),affineMagnitude(affineRange(b)));
        constexpr double eps=4*std::numeric_limits<double>::epsilon(),tiny=4*std::numeric_limits<double>::denorm_min();
        addSymmetric(out,Arithmetic::add(Arithmetic::multiply({eps,eps},{m,m}),{tiny,tiny}));
    }
    return out;
}
double affineRepresentative(Interval a)
{range(a);return std::max(a.lo,std::min(a.hi,a.lo*.5+a.hi*.5));}
PolynomialAffineCoordinate exportAffine(const Affine& a)
{
    PolynomialAffineCoordinate out;out.center=affineRepresentative(a.center);out.tangent=affineRepresentative(a.tangent);out.time=affineRepresentative(a.time);
    auto residual=Arithmetic::add(a.remainder,Arithmetic::add(a.center,{-out.center,-out.center}));
    residual=Arithmetic::add(residual,Arithmetic::multiply(Arithmetic::add(a.tangent,{-out.tangent,-out.tangent}),{-1,1}));
    residual=Arithmetic::add(residual,Arithmetic::multiply(Arithmetic::add(a.time,{-out.time,-out.time}),{-1,1}));
    out.remainder=range(residual);return out;
}
Affine importAffine(const PolynomialAffineCoordinate& coordinate)
{
    require(std::isfinite(coordinate.center) && std::isfinite(coordinate.tangent) && std::isfinite(coordinate.time),"nonfinite affine coefficient");
    Affine out;out.center={coordinate.center,coordinate.center};out.tangent={coordinate.tangent,coordinate.tangent};
    out.time={coordinate.time,coordinate.time};out.remainder=interval(coordinate.remainder);range(out.remainder);return out;
}
Affine divideAffineInterval(const Affine& a,Interval denominator)
{
    range(denominator);require(denominator.lo>0,"unresolved positive affine chart width");
    const auto divide=[&](Interval n){return hull(Arithmetic::divide(n,denominator.lo),Arithmetic::divide(n,denominator.hi));};
    Affine out;out.center=divide(a.center);out.tangent=divide(a.tangent);out.time=divide(a.time);out.remainder=divide(a.remainder);
    affineRange(out);return out;
}
Interval parameterAt(double point,double lo,double hi)
{
    require(std::isfinite(point) && std::isfinite(lo) && std::isfinite(hi) && lo<hi && point>=lo && point<=hi,"affine point outside parameter extent");
    auto denominator=differenceInterval(hi,lo);
    // A positive difference of finite binary64 endpoints is at least one
    // denorm_min as a real number. Direct endpoint division avoids an
    // overflowing reciprocal for a subnormal but valid parameter extent.
    denominator.lo=std::max(denominator.lo,std::numeric_limits<double>::denorm_min());range(denominator);
    const auto numerator=Arithmetic::add(Arithmetic::multiply(differenceInterval(point,lo),{2,2}),negative(denominator));
    auto q=hull(Arithmetic::divide(numerator,denominator.lo),Arithmetic::divide(numerator,denominator.hi));
    q.lo=std::max(-1.,q.lo);q.hi=std::min(1.,q.hi);range(q);return q;
}
}
PolynomialAffineRk2FaceBounds PolynomialVelocitySlab::encloseAffineFaceRk2(const LevelSetBoundaryTrace::Domain& domain,
    Vector2D a,Vector2D b,double first,double last) const
{
    Detail::ExactExpansion::validateEnvironment();validateTimes(first,last,begin_,end_);
    const bool horizontal=a.y==b.y && a.x<b.x,vertical=a.x==b.x && a.y<b.y;
    require(first<last && finite(a) && finite(b) && (horizontal || vertical),"affine RK2 requires positive ascending face/time extent");
    PolynomialAffineRk2FaceBounds result;result.faceBegin=a;result.faceEnd=b;result.timeBegin=first;result.timeEnd=last;
    result.branch=encloseInteriorRk2TightTime(domain,a,b,first,last);
    const AffineVector input=horizontal?AffineVector{affineParameter(a.x,b.x,false),Affine(a.y)}:
        AffineVector{Affine(a.x),affineParameter(a.y,b.y,false)};
    const auto time=affineParameter(first,last,true),duration=time-Affine(begin_),theta=divideAffine(duration,end_-begin_);
    const auto spatial=[&](const std::array<Vector2D,16>& coefficients,const AffineVector& p){
        const auto x=normalized(p.x,data_.origin.x,data_.scale.x),y=normalized(p.y,data_.origin.y,data_.scale.y);
        return AffineVector{component(coefficients,x,y,data_.scale,false).value,component(coefficients,x,y,data_.scale,true).value};
    };
    const auto arrival=[&](const AffineVector& p){const auto old=spatial(data_.previous,p),end=spatial(data_.arrival,p);
        return AffineVector{mixAffine(old.x,end.x,theta),mixAffine(old.y,end.y,theta)};};
    const auto start=arrival(input),half=LevelSetBoundaryTrace::Detail::rk2Midpoint(input,start,duration);
    const auto newMiddle=arrival(half),oldMiddle=spatial(data_.previous,half),middle=(newMiddle+oldMiddle)*Affine(.5);
    auto target=LevelSetBoundaryTrace::Detail::rk2Target(input,middle,duration);
    if(affineRange(duration).lo<=std::numeric_limits<double>::denorm_min()){
        const auto tinyBranch=[&](Affine& out,const Affine& position,const Affine& newSpeed,const Affine& meanSpeed){
            if(exact(newSpeed,0) && exact(meanSpeed,0))return;
            // dt/2==0 can select arrival-only instead of the mean. Then
            // dt<=denorm_min. Bound the difference of both product/subtract
            // sequences, including both original rounding errors.
            const double x=affineMagnitude(affineRange(position)),u=affineMagnitude(affineRange(newSpeed)),v=affineMagnitude(affineRange(meanSpeed));
            constexpr double eps=8*std::numeric_limits<double>::epsilon(),tiny=8*std::numeric_limits<double>::denorm_min();
            const auto sum=Arithmetic::add({1,1},Arithmetic::add({u,u},{v,v}));
            addSymmetric(out,Arithmetic::add(Arithmetic::multiply({eps,eps},{x,x}),Arithmetic::multiply({tiny,tiny},sum)));
        };
        tinyBranch(target.x,input.x,newMiddle.x,middle.x);tinyBranch(target.y,input.y,newMiddle.y,middle.y);
    }
    result.departure={exportAffine(target.x),exportAffine(target.y)};return result;
}
std::array<PolynomialVelocityRange,2> PolynomialAffineRk2FaceBounds::encloseDeparture(Vector2D point,double time) const
{
    Detail::ExactExpansion::validateEnvironment();const bool horizontal=faceBegin.y==faceEnd.y && faceBegin.x<faceEnd.x,
        vertical=faceBegin.x==faceEnd.x && faceBegin.y<faceEnd.y;
    require((horizontal || vertical) && finite(point) && (horizontal?point.y==faceBegin.y:point.x==faceBegin.x),"affine query not on retained face");
    const auto qs=parameterAt(horizontal?point.x:point.y,horizontal?faceBegin.x:faceBegin.y,horizontal?faceEnd.x:faceEnd.y),
        qt=parameterAt(time,timeBegin,timeEnd);std::array<PolynomialVelocityRange,2> result{};
    for(int k=0;k<2;++k){const auto& a=departure[k];require(std::isfinite(a.center) && std::isfinite(a.tangent) && std::isfinite(a.time),"nonfinite affine coefficient");
        auto value=Arithmetic::add({a.center,a.center},Arithmetic::add(Arithmetic::multiply({a.tangent,a.tangent},qs),Arithmetic::multiply({a.time,a.time},qt)));
        value=Arithmetic::add(value,interval(a.remainder));value.lo=std::max(value.lo,k?branch.departureLower.y:branch.departureLower.x);
        value.hi=std::min(value.hi,k?branch.departureUpper.y:branch.departureUpper.x);result[k]=range(value);
    }
    return result;
}
PolynomialVelocityRange PolynomialAffineCoordinate::enclose(PolynomialVelocityRange tangentParameter,
    PolynomialVelocityRange timeParameter) const
{
    Detail::ExactExpansion::validateEnvironment();const auto a=importAffine(*this);
    const auto qs=interval(tangentParameter),qt=interval(timeParameter);
    range(qs);range(qt);require(qs.lo>=-1 && qs.hi<=1 && qt.lo>=-1 && qt.hi<=1,"affine parameter outside unit square");
    return range(Arithmetic::add(a.center,Arithmetic::add(Arithmetic::add(Arithmetic::multiply(a.tangent,qs),Arithmetic::multiply(a.time,qt)),a.remainder)));
}
PolynomialAffineCoordinate encloseTensorQuarticAffine(const std::array<std::array<double,2>,25>& coefficients,
    Vector2D lo,Vector2D hi,const std::array<PolynomialAffineCoordinate,2>& physical)
{
    Detail::ExactExpansion::validateEnvironment();require(finite(lo) && finite(hi) && lo.x<hi.x && lo.y<hi.y,"invalid affine Q4 physical cell bounds");
    for(const auto c:coefficients)range({c[0],c[1]});
    const auto x=divideAffineInterval(importAffine(physical[0])-Affine(lo.x),differenceInterval(hi.x,lo.x)),
        y=divideAffineInterval(importAffine(physical[1])-Affine(lo.y),differenceInterval(hi.y,lo.y));
    Affine out;
    for(int j=4;j>=0;--j){Affine row;
        for(int i=4;i>=0;--i){Affine coefficient;coefficient.center={coefficients[5*j+i][0],coefficients[5*j+i][1]};row=row*x+coefficient;}
        out=out*y+row;
    }
    return exportAffine(out);
}
PolynomialFaceIntegralBounds PolynomialVelocitySlab::encloseNormalFaceIntegralWithAffinePhase(Vector2D a,Vector2D b,
    double begin,double end,bool radial,const std::vector<PolynomialAffineCoordinate>& forms,
    long double liquidLower,long double liquidUpper,long double gasLower,long double gasUpper) const
{
    Detail::ExactExpansion::validateEnvironment();require(!forms.empty(),"empty affine phase coverage");
    contentInterval(liquidLower,liquidUpper);contentInterval(gasLower,gasUpper);
    require(gasLower<=liquidLower && gasUpper<=liquidUpper,"unordered affine phase weights");
    auto result=encloseNormalFaceIntegral(a,b,begin,end,radial,gasLower,liquidUpper);
    const auto& base=forms.front();importAffine(base);Interval residual{};bool first=true;
    const auto coefficientDifference=[](double x,double y){return x==y?Interval{}:differenceInterval(x,y);};
    for(const auto& f:forms){importAffine(f);
        auto delta=Arithmetic::add(coefficientDifference(f.center,base.center),Arithmetic::add(
            Arithmetic::multiply(coefficientDifference(f.tangent,base.tangent),{-1,1}),
            Arithmetic::multiply(coefficientDifference(f.time,base.time),{-1,1})));
        delta=Arithmetic::add(delta,interval(f.remainder));range(delta);residual=first?delta:hull(residual,delta);first=false;
    }
    namespace E=Detail::ExactExpansion;namespace F=Detail::FaceTime;using R=F::FaceTimeRational;
    const auto exact=[](double x){return E::fromDouble(x);};
    const auto difference=[&](double x,double y){return E::subtract(exact(x),exact(y));};
    const bool horizontal=a.y==b.y;
    const R length(horizontal?difference(b.x,a.x):difference(b.y,a.y),exact(1)),duration(difference(end,begin),exact(1));
    const auto subset=[&](bool gas){
        // Upper phi<0 certifies liquid; -lower phi<0 certifies gas.
        const auto sign=exact(gas?-1:1),s=E::multiply(sign,exact(base.tangent)),t=E::multiply(sign,exact(base.time));
        const auto c=E::multiply(sign,E::add(exact(base.center),exact(gas?residual.lo:residual.hi)));
        const auto moments=F::negativeHalfPlane(E::subtract(E::subtract(c,s),t),E::multiply(exact(2),s),E::multiply(exact(2),t));
        require(E::sign(moments.area.n)>=0 && E::sign((R(1)-moments.area).n)>=0
            && E::sign(moments.first.n)>=0 && E::sign((moments.area-moments.first).n)>=0,"invalid exact affine subset moments");
        const auto weight=radial?R(a.x)*moments.area+R(difference(b.x,a.x),exact(1))*moments.first:moments.area;
        const auto measure=length*duration*weight;
        require(E::sign(measure.n)>=0,"negative exact affine subset measure");
        if(E::sign(measure.n)==0)return Interval{};
        const auto n=E::enclosure(measure.n),d=E::enclosure(measure.d);
        require(d.lo>0 && std::isfinite(d.hi),"unresolved affine subset denominator");
        auto value=nonnegative(hull(Arithmetic::divide({n.lo,n.hi},d.lo),Arithmetic::divide({n.lo,n.hi},d.hi)));
        if(radial){constexpr double pi=3.141592653589793238462643383279502884;
            value=nonnegative(Arithmetic::multiply(value,Arithmetic::multiply({2,2},{Arithmetic::down(pi),Arithmetic::up(pi)})));}
        range(value);return value;
    };
    const auto liquidMeasure=subset(false),gasMeasure=subset(true);
    const auto velocityFirst=enclose(a,b,begin,begin).value[horizontal?1:0],velocityLast=enclose(a,b,end,end).value[horizontal?1:0];
    auto velocity=hull(interval(velocityFirst),interval(velocityLast));
    // The same real linear time interpolation lies in the endpoint hull.
    // Also retain the existing 8-epsilon/8-tiny envelope for original RN
    // mixing branches, including the distinct midpoint mean.
    const double speedMagnitude=std::max(std::abs(velocity.lo),std::abs(velocity.hi));
    const auto speedError=Arithmetic::add(Arithmetic::multiply({8*std::numeric_limits<double>::epsilon(),8*std::numeric_limits<double>::epsilon()},
        {speedMagnitude,speedMagnitude}),{8*std::numeric_limits<double>::denorm_min(),8*std::numeric_limits<double>::denorm_min()});
    velocity=Arithmetic::add(velocity,{-speedError.hi,speedError.hi});range(velocity);
    const auto wl0=contentInterval(liquidLower,liquidLower),wl1=contentInterval(liquidUpper,liquidUpper),
        wg0=contentInterval(gasLower,gasLower),wg1=contentInterval(gasUpper,gasUpper);
    const auto liquidGain=nonnegative(Arithmetic::add(wl0,{-wg0.hi,-wg0.lo})),gasGain=nonnegative(Arithmetic::add(wl1,{-wg1.hi,-wg1.lo}));
    for(int k=0;k<2;++k){const auto total=interval(result.volumes[k]);
        const double minimumSpeed=std::max(0.,k?-velocity.hi:velocity.lo);
        // Lowering a certified subset bound to T.lower preserves validity
        // and the existing componentwise liquid/total ordering convention.
        const double knownLiquid=std::min(total.lo,nonnegative(Arithmetic::multiply({minimumSpeed,minimumSpeed},liquidMeasure)).lo),
            knownGas=std::min(total.lo,nonnegative(Arithmetic::multiply({minimumSpeed,minimumSpeed},gasMeasure)).lo);
        const auto low=Arithmetic::add(Arithmetic::multiply(wg0,total),Arithmetic::multiply(liquidGain,{knownLiquid,knownLiquid}));
        const auto removed=Arithmetic::multiply(gasGain,{knownGas,knownGas}),high=Arithmetic::add(Arithmetic::multiply(wl1,total),{-removed.hi,-removed.lo});
        result.volumes[k+2]={std::max(result.volumes[k+2].lower,std::min(total.lo,std::max(0.,low.lo))),
            std::min(result.volumes[k+2].upper,std::max(0.,high.hi))};
    }
    validIntegral(result);return result;
}
}
