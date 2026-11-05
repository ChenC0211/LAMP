#include "Core/Geometry/PlanarFaceTimeIntegral.h"
#include "Core/Geometry/Detail/ExactFaceTimeRational.h"
#include "Core/Mesh/CartesianMeasure.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace LAMP::Geometry::PlanarFaceTime
{
namespace
{
namespace E=Detail::ExactExpansion;
namespace F=Detail::FaceTime;
using X=E::Expansion;
using R=F::FaceTimeRational;
using Quadratic=std::array<X,3>;
using Polynomial=std::array<X,6>;
constexpr unsigned refinementLimit=256;
constexpr long double relativeUncertainty=32*std::numeric_limits<long double>::epsilon();
void require(bool ok,const char* why)
{if(!ok)throw std::invalid_argument(std::string("[PlanarFaceTime] ")+why);}
X scalar(long double value)
{
    require(std::isfinite(value),"nonfinite original input");X result;
    constexpr unsigned count=(std::numeric_limits<long double>::digits+52)/53+1;
    for(unsigned k=0;value!=0&&k<count;++k){
        const double part=static_cast<double>(value);
        require(std::isfinite(part)&&part!=0,"original input outside exact arithmetic scale");
        E::accumulate(result,E::fromDouble(part));value-=static_cast<long double>(part);}
    require(value==0,"original input has an unsupported exact tail");return result;
}
X neg(const X& a){return E::subtract(X{},a);}
X times(const X& a,double b){return E::multiply(a,E::fromDouble(b));}
int compare(const X& a,const X& b){return E::sign(E::subtract(a,b));}
int compare(const R& a,const R& b)
{return E::sign(E::subtract(E::multiply(a.n,b.d),E::multiply(b.n,a.d)));}
void normalize(Quadratic& p)
{
    E::ExponentRange range;for(const auto& c:p)range.include(c);
    const int shift=range.centeredShift();for(auto& c:p)c=E::scalePow2(c,shift);
}
X at(const Quadratic& p,const X& t)
{return E::add(p[0],E::multiply(t,E::add(p[1],E::multiply(t,p[2]))));}
X at(const Quadratic& p,const R& t)
{
    return E::add(E::multiply(p[0],E::multiply(t.d,t.d)),
        E::add(E::multiply(p[1],E::multiply(t.n,t.d)),E::multiply(p[2],E::multiply(t.n,t.n))));
}
int slope(const Quadratic& p,const R& t)
{return E::sign(E::add(E::multiply(p[1],t.d),times(E::multiply(p[2],t.n),2)));}
struct Root
{
    Quadratic polynomial{};
    R rational{};
    int branch=0; // 0 exact rational; -1/+1 roots of a positive-leading quadratic.
    X lo=E::fromDouble(0),hi=E::fromDouble(1);
    unsigned refinements=0;
};
Root constant(double value)
{Root r;r.rational=R(value);r.lo=r.hi=E::fromDouble(value);return r;}
int compare(const Root& r,const X& t)
{
    if(!r.branch)return E::sign(E::subtract(r.rational.n,E::multiply(t,r.rational.d)));
    const int derivative=E::sign(E::add(r.polynomial[1],times(E::multiply(r.polynomial[2],t),2)));
    if(r.branch<0){if(derivative>=0)return -1;return E::sign(at(r.polynomial,t));}
    if(derivative<=0)return 1;return -E::sign(at(r.polynomial,t));
}
bool equal(const Root& a,const Root& b)
{
    if(!a.branch&&!b.branch)return compare(a.rational,b.rational)==0;
    const auto rationalQuadratic=[](const Root& rational,const Root& quadratic){
        return E::sign(at(quadratic.polynomial,rational.rational))==0&&
            slope(quadratic.polynomial,rational.rational)==quadratic.branch;};
    if(!a.branch)return rationalQuadratic(a,b);
    if(!b.branch)return rationalQuadratic(b,a);
    Quadratic eliminated{};
    for(unsigned k=0;k<3;++k)eliminated[k]=E::subtract(
        E::multiply(a.polynomial[k],b.polynomial[2]),E::multiply(b.polynomial[k],a.polynomial[2]));
    if(E::sign(eliminated[0])==0&&E::sign(eliminated[1])==0)return a.branch==b.branch;
    if(E::sign(eliminated[1])==0)return false;
    Root candidate;candidate.rational=R(neg(eliminated[0]),eliminated[1]);
    return rationalQuadratic(candidate,a)&&rationalQuadratic(candidate,b);
}
void refine(Root& root)
{
    if(compare(root.lo,root.hi)==0)return;
    require(root.refinements<refinementLimit,"exact temporal root refinement budget exhausted");
    const auto middle=E::scalePow2(E::add(root.lo,root.hi),-1);
    require(compare(middle,root.lo)>0&&compare(middle,root.hi)<0,"unrepresentable exact temporal midpoint");
    const int side=compare(root,middle);++root.refinements;
    if(side==0)root.lo=root.hi=middle;else if(side>0)root.lo=middle;else root.hi=middle;
}
double approximation(const X& value)
{
    double result=0;for(double part:value.components())result+=part;return result;
}
void seed(Root& root)
{
    double guess=0;
    if(!root.branch)guess=approximation(root.rational.n)/approximation(root.rational.d);
    else {
        const auto& p=root.polynomial;
        const double a=approximation(p[0]),b=approximation(p[1]),c=approximation(p[2]);
        const auto discriminant=E::subtract(E::multiply(p[1],p[1]),times(E::multiply(p[0],p[2]),4));
        const double disc=approximation(discriminant);
        if(!(disc>0&&std::isfinite(disc)&&std::isfinite(a)&&std::isfinite(b)&&c>0&&std::isfinite(c)))return;
        const double q=-.5*(b+std::copysign(std::sqrt(disc),b));
        if(q==0)return;const double first=q/c,second=a/q;
        guess=root.branch<0?std::min(first,second):std::max(first,second);
    }
    if(!(guess>=0&&guess<=1&&std::isfinite(guess)))return;
    double lo=guess,hi=guess;
    for(unsigned i=0;i<8;++i){lo=std::max(0.,std::nextafter(lo,-std::numeric_limits<double>::infinity()));
        hi=std::min(1.,std::nextafter(hi,std::numeric_limits<double>::infinity()));}
    const auto low=E::fromDouble(lo),high=E::fromDouble(hi);
    const int left=compare(root,low),right=compare(root,high);
    // Floating roots only propose a bracket; both endpoints are certified
    // against the original exact polynomial and the selected root identity.
    if(left>=0&&right<=0){root.lo=low;root.hi=high;
        if(left==0)root.hi=root.lo;else if(right==0)root.lo=root.hi;}
}
void addRoots(Quadratic polynomial,std::vector<Root>& roots,unsigned limit)
{
    normalize(polynomial);
    const auto append=[&](Root root){
        const int first=compare(root,E::fromDouble(0)),last=compare(root,E::fromDouble(1));
        if(first>0&&last<0){seed(root);roots.push_back(std::move(root));
            require(roots.size()<=limit,"exact temporal event budget exhausted");}};
    const int c=E::sign(polynomial[2]),b=E::sign(polynomial[1]);
    if(c==0){if(b!=0){Root root;root.rational=R(neg(polynomial[0]),polynomial[1]);append(std::move(root));}return;}
    if(c<0)for(auto& coefficient:polynomial)coefficient=neg(coefficient);
    const auto discriminant=E::subtract(E::multiply(polynomial[1],polynomial[1]),
        times(E::multiply(polynomial[0],polynomial[2]),4));
    const int sign=E::sign(discriminant);if(sign<0)return;
    if(sign==0){Root root;root.rational=R(neg(polynomial[1]),times(polynomial[2],2));append(std::move(root));return;}
    for(int branch:{-1,1}){Root root;root.polynomial=polynomial;root.branch=branch;append(std::move(root));}
}
int order(Root& a,Root& b)
{
    if(equal(a,b))return 0;
    for(;;){
        if(compare(a.hi,b.lo)<=0)return -1;
        if(compare(b.hi,a.lo)<=0)return 1;
        const auto wa=E::subtract(a.hi,a.lo),wb=E::subtract(b.hi,b.lo);
        if(compare(wa,wb)>=0)refine(a);else refine(b);
    }
}
void sortEvents(std::vector<Root>& roots)
{
    std::vector<Root> sorted;
    for(auto& root:roots){
        std::size_t i=0;bool duplicate=false;
        for(;i<sorted.size();++i){const int side=order(root,sorted[i]);
            if(side==0){duplicate=true;break;}if(side<0)break;}
        if(!duplicate)sorted.insert(sorted.begin()+i,std::move(root));
    }
    roots=std::move(sorted);
}
struct ExactPoint {X x,y;};
ExactPoint point(Point p){return {scalar(p.x),scalar(p.y)};}
ExactPoint point(AnchoredPoint p)
{
    return {E::add(scalar(p.anchor.x),scalar(p.offset.x)),
        E::add(scalar(p.anchor.y),scalar(p.offset.y))};
}
ExactPoint subtract(const ExactPoint& a,const ExactPoint& b){return {E::subtract(a.x,b.x),E::subtract(a.y,b.y)};}
X dot(const ExactPoint& a,const ExactPoint& b){return E::add(E::multiply(a.x,b.x),E::multiply(a.y,b.y));}
X cross(const ExactPoint& a,const ExactPoint& b){return E::subtract(E::multiply(a.x,b.y),E::multiply(a.y,b.x));}
int orientation(const std::vector<ExactPoint>& polygon)
{
    int sign=0;
    for(std::size_t i=0;i<polygon.size();++i){const auto edge=subtract(polygon[(i+1)%polygon.size()],polygon[i]);
        for(std::size_t j=0;j<polygon.size();++j){const int next=E::sign(cross(edge,subtract(polygon[j],polygon[i])));
            if(next!=0){require(sign==0||sign==next,"concave or crossed original polygon");sign=next;}}}
    return sign;
}
struct Bound {Quadratic numerator{};X denominator=E::fromDouble(1);};
void normalize(Bound& b){E::centerHomogeneous(b.numerator,b.denominator);}
Quadratic difference(const Bound& a,const Bound& b)
{
    Quadratic p{};for(unsigned k=0;k<3;++k)p[k]=E::subtract(
        E::multiply(a.numerator[k],b.denominator),E::multiply(b.numerator[k],a.denominator));
    return p;
}
int compare(const Bound& a,const Bound& b,const X& t)
{return E::sign(at(difference(a,b),t));}
struct Density {Polynomial numerator{};X denominator=E::fromDouble(1);int degree=0;};
Density density(const Bound& lower,const Bound& upper,const Quadratic& normal,int sign,
    const X& radius,const X& radialDelta,bool rz)
{
    Density result;Polynomial wet{};
    auto differenceNumerator=difference(upper,lower);
    const auto denominator=E::multiply(upper.denominator,lower.denominator);
    if(!rz){for(unsigned k=0;k<3;++k)wet[k]=differenceNumerator[k];result.denominator=denominator;}
    else {
        const auto upperSquare=E::multiply(upper.denominator,upper.denominator);
        const auto lowerSquare=E::multiply(lower.denominator,lower.denominator);
        for(unsigned k=0;k<3;++k)wet[k]=E::multiply(E::multiply(radius,differenceNumerator[k]),denominator);
        if(E::sign(radialDelta)!=0)for(unsigned i=0;i<3;++i)for(unsigned j=0;j<3;++j){
            const auto term=E::subtract(E::multiply(E::multiply(upper.numerator[i],upper.numerator[j]),lowerSquare),
                E::multiply(E::multiply(lower.numerator[i],lower.numerator[j]),upperSquare));
            E::accumulate(wet[i+j],E::scalePow2(E::multiply(radialDelta,term),-1));}
        result.denominator=E::multiply(denominator,denominator);
    }
    for(unsigned i=0;i<5;++i)for(unsigned j=0;j<2;++j){
        const auto term=E::multiply(wet[i],normal[j]);
        E::accumulate(result.numerator[i+j],sign>0?term:neg(term));}
    E::centerHomogeneous(result.numerator,result.denominator);
    result.degree=5;while(result.degree>0&&E::sign(result.numerator[result.degree])==0)--result.degree;
    return result;
}
int choose(int n,int k)
{int result=1;for(int j=1;j<=k;++j)result=result*(n-j+1)/j;return result;}
Polynomial local(const Density& density,const X& begin,const X& width)
{
    Polynomial a{},w{},result{};a[0]=w[0]=E::fromDouble(1);
    for(int k=1;k<=density.degree;++k){a[k]=E::multiply(a[k-1],begin);w[k]=E::multiply(w[k-1],width);}
    for(int j=0;j<=density.degree;++j){
        X coefficient;for(int k=j;k<=density.degree;++k)E::accumulate(coefficient,
            times(E::multiply(density.numerator[k],a[k-j]),static_cast<double>(choose(k,j))));
        result[j]=E::multiply(coefficient,w[j]);}
    return result;
}
X integralNumerator(const Density& density,const X& begin,const X& end)
{
    const auto width=E::subtract(end,begin);require(E::sign(width)>=0,"inverted exact integration interval");
    if(E::sign(width)==0)return X{};
    const auto coefficients=local(density,begin,width);X numerator;
    for(int k=0;k<=density.degree;++k)E::accumulate(numerator,times(coefficients[k],60./(k+1)));
    return E::multiply(width,numerator);
}
X scaledMaximum(const Density& density,const X& begin,const X& end)
{
    const auto coefficients=local(density,begin,E::subtract(end,begin));X maximum;
    for(int i=0;i<=density.degree;++i){X b;
        for(int j=0;j<=i;++j)E::accumulate(b,times(coefficients[j],
            static_cast<double>(choose(i,j)*(120/choose(density.degree,j)))));
        if(compare(b,maximum)>0)maximum=std::move(b);}
    return maximum; // 120*max Bernstein coefficient, nonnegative upper bound.
}
X stripNumerator(const Density& density,const Root& left,const Root& right)
{
    X numerator;
    if(compare(left.lo,left.hi)!=0)E::accumulate(numerator,
        E::multiply(E::subtract(left.hi,left.lo),scaledMaximum(density,left.lo,left.hi)));
    if(compare(right.lo,right.hi)!=0)E::accumulate(numerator,
        E::multiply(E::subtract(right.hi,right.lo),scaledMaximum(density,right.lo,right.hi)));
    return numerator;
}
struct PositiveSum
{
    long double value=0,error=0;
    void add(long double x){require(x>=0&&std::isfinite(x),"invalid positive gross contribution");
        const auto y=x-error,t=value+y;require(std::isfinite(t),"overflowing positive gross sum");
        error=(t-value)-y;value=t;}
};
long double certifiedIntegral(const Density& density,Root& left,Root& right,const R& units,Gross& report)
{
    bool zero=true;for(const auto& c:density.numerator)zero&=E::sign(c)==0;if(zero)return 0;
    for(;;){
        if(compare(left.hi,right.lo)>=0){refine(left);refine(right);continue;}
        // Both quantities retain the SAME density denominator D:
        // I=Ni/(60D), E=Ns/(120D). Cross-multiplying separate rational
        // objects would duplicate D and can overflow a representable object.
        const auto inside=integralNumerator(density,left.hi,right.lo);
        require(E::sign(inside)>0,"positive branch has nonpositive exact integral");
        const auto error=stripNumerator(density,left,right);
        const auto allowed=E::multiply(scalar(2*relativeUncertainty),inside);
        if(compare(error,allowed)<=0){
            const R midpoint(E::add(times(inside,4),error),times(density.denominator,240));
            const auto physical=midpoint*units;
            const long double value=F::faceTimeValue(physical);
            require(value>0,"positive geometric integral is unrepresentable");
            const auto relative=E::sign(error)==0?0:F::faceTimeValue(R(error,times(inside,2)));
            report.maximumRelativeUncertainty=std::max(report.maximumRelativeUncertainty,relative);
            return value;
        }
        refine(left);refine(right);
    }
}
void addGross(Gross& result,std::array<PositiveSum,4>& sums,long double value,int normalSign,bool half)
{
    if(normalSign>0){sums[2].add(value);if(half)sums[0].add(value);}
    else {sums[3].add(value);if(half)sums[1].add(value);}
    result.halfPositive=sums[0].value;result.halfNegative=sums[1].value;
    result.fullPositive=sums[2].value;result.fullNegative=sums[3].value;
}
X duration(long double begin,long double end)
{const auto dt=E::subtract(scalar(end),scalar(begin));require(E::sign(dt)>0,"nonpositive complete time interval");return dt;}
struct Axis
{
    R low{},high{};
};
Axis displacement(const X& first,const X& last,const X& dt)
{
    // One homogeneous velocity scale preserves the exact trajectory while
    // avoiding an unrepresentable physical velocity square at its extremum.
    E::ExponentRange range;range.include(first);range.include(last);
    const int shift=range.centeredShift();
    const auto normalizedFirst=E::scalePow2(first,shift),normalizedLast=E::scalePow2(last,shift);
    const auto delta=E::subtract(normalizedLast,normalizedFirst);
    const R units(E::scalePow2(E::fromDouble(1),-shift),E::fromDouble(1));
    const R finish=R(E::scalePow2(E::multiply(dt,E::add(normalizedFirst,normalizedLast)),-1),E::fromDouble(1))*units;
    Axis result;if(compare(finish,result.low)<0)result.low=finish;if(compare(finish,result.high)>0)result.high=finish;
    if(E::sign(delta)!=0){const R time(neg(normalizedFirst),delta);
        if(compare(time,R(0))>0&&compare(time,R(1))<0){
            const R critical=R(neg(E::multiply(dt,E::multiply(normalizedFirst,normalizedFirst))),times(delta,2))*units;
            if(compare(critical,result.low)<0)result.low=critical;if(compare(critical,result.high)>0)result.high=critical;}}
    return result;
}
std::array<R,4> departure(Point faceBegin,Point faceEnd,Point velocityBegin,Point velocityEnd,
    long double timeBegin,long double timeEnd)
{
    const auto a=point(faceBegin),b=point(faceEnd),v0=point(velocityBegin),v1=point(velocityEnd);
    const auto dt=duration(timeBegin,timeEnd);
    const auto x=displacement(v0.x,v1.x,dt),y=displacement(v0.y,v1.y,dt);
    const auto minX=compare(a.x,b.x)<0?a.x:b.x,maxX=compare(a.x,b.x)>0?a.x:b.x;
    const auto minY=compare(a.y,b.y)<0?a.y:b.y,maxY=compare(a.y,b.y)>0?a.y:b.y;
    return {{R(minX,E::fromDouble(1))-x.high,R(maxX,E::fromDouble(1))-x.low,
        R(minY,E::fromDouble(1))-y.high,R(maxY,E::fromDouble(1))-y.low}};
}
long double outward(const R& value,bool upper)
{
    const auto n=E::enclosure(value.n),d=E::enclosure(value.d);
    require(d.lo>0&&std::isfinite(d.hi),"unresolved outward rational denominator");
    const long double candidates[4]={static_cast<long double>(n.lo)/d.lo,static_cast<long double>(n.lo)/d.hi,
        static_cast<long double>(n.hi)/d.lo,static_cast<long double>(n.hi)/d.hi};
    const auto result=upper?*std::max_element(candidates,candidates+4):*std::min_element(candidates,candidates+4);
    require(std::isfinite(result),"unrepresentable departure bound");
    return std::nextafter(result,upper?std::numeric_limits<long double>::infinity():-std::numeric_limits<long double>::infinity());
}
template<class Site>
Gross integratePolygonImpl(const std::vector<Site>& polygon,Point faceBegin,Point faceEnd,
    Point velocityBegin,Point velocityEnd,long double timeBegin,long double timeEnd,
    Point frameOrigin,bool rz,unsigned eventLimit,Reflection reflection)
{
    E::validateEnvironment();require(eventLimit>=3&&eventLimit<=4096,"invalid temporal event budget");
    require((reflection.signX==1||reflection.signX==-1)&&(reflection.signY==1||reflection.signY==-1),"invalid exact reflection");
    const auto a=point(faceBegin),b=point(faceEnd),d=subtract(b,a),v0=point(velocityBegin),v1=point(velocityEnd);
    const auto dt=duration(timeBegin,timeEnd);
    const auto origin=point(frameOrigin),translation=point(reflection.translation);
    const bool horizontal=E::sign(d.y)==0;
    require((horizontal&&E::sign(d.x)!=0)||(!horizontal&&E::sign(d.x)==0&&E::sign(d.y)!=0),"face must be nondegenerate and axis aligned");
    if(rz)require(E::sign(E::add(origin.x,a.x))>=0&&E::sign(E::add(origin.x,b.x))>=0,"negative Eulerian face radius");
    std::vector<ExactPoint> points;points.reserve(polygon.size());
    for(auto original:polygon){auto p=point(original);
        p.x=E::add(reflection.signX>0?p.x:neg(p.x),translation.x);
        p.y=E::add(reflection.signY>0?p.y:neg(p.y),translation.y);points.push_back(std::move(p));}
    const int cyclic=orientation(points);if(cyclic==0||points.size()<3)return {};
    if(cyclic<0)std::reverse(points.begin(),points.end());
    Bound zero,one;one.numerator[0]=E::fromDouble(1);
    std::vector<Bound> lower{zero},upper{one};std::vector<Quadratic> predicates;
    const auto velocityDelta=subtract(v1,v0);
    for(std::size_t i=0;i<points.size();++i){
        const auto edge=subtract(points[(i+1)%points.size()],points[i]);const ExactPoint n{edge.y,neg(edge.x)};
        if(E::sign(n.x)==0&&E::sign(n.y)==0)continue;
        const auto lambda=dot(n,d);
        Quadratic time{{dot(n,subtract(points[i],a)),E::multiply(dt,dot(n,v0)),
            E::scalePow2(E::multiply(dt,dot(n,velocityDelta)),-1)}};
        if(E::sign(lambda)==0){normalize(time);predicates.push_back(std::move(time));}
        else {Bound bound;bound.numerator=std::move(time);bound.denominator=lambda;
            if(E::sign(lambda)<0){bound.denominator=neg(lambda);for(auto& c:bound.numerator)c=neg(c);}
            normalize(bound);(E::sign(lambda)>0?upper:lower).push_back(std::move(bound));}
    }
    Quadratic normal{{horizontal?v0.y:v0.x,horizontal?velocityDelta.y:velocityDelta.x,X{}}};
    std::vector<Root> events{constant(0),constant(.5),constant(1)};
    std::vector<Bound> bounds=lower;bounds.insert(bounds.end(),upper.begin(),upper.end());
    for(std::size_t i=0;i<bounds.size();++i)for(std::size_t j=0;j<i;++j)addRoots(difference(bounds[i],bounds[j]),events,eventLimit);
    for(const auto& p:predicates)addRoots(p,events,eventLimit);addRoots(normal,events,eventLimit);sortEvents(events);
    const auto length=horizontal?d.x:d.y;
    const auto positiveLength=E::sign(length)<0?neg(length):length;
    R units(dt,E::fromDouble(1));units=units*R(positiveLength,E::fromDouble(1));
    if(rz)units=units*R(scalar(CartesianMeasure::twoPi),E::fromDouble(1));
    Gross result;result.temporalEvents=events.size();std::array<PositiveSum,4> sums{};
    const auto radius=E::add(origin.x,a.x);
    for(std::size_t i=1;i<events.size();++i){auto& left=events[i-1];auto& right=events[i];
        while(compare(left.hi,right.lo)>=0){refine(left);refine(right);}
        const auto sample=E::scalePow2(E::add(left.hi,right.lo),-1);
        if(std::any_of(predicates.begin(),predicates.end(),[&](const Quadratic& p){return E::sign(at(p,sample))<0;}))continue;
        std::size_t l=0,u=0;for(std::size_t j=1;j<lower.size();++j)if(compare(lower[j],lower[l],sample)>0)l=j;
        for(std::size_t j=1;j<upper.size();++j)if(compare(upper[j],upper[u],sample)<0)u=j;
        if(compare(upper[u],lower[l],sample)<=0)continue;
        const int sign=E::sign(at(normal,sample));if(sign==0)continue;
        const auto branch=density(lower[l],upper[u],normal,sign,radius,d.x,rz);
        const auto value=certifiedIntegral(branch,left,right,units,result);++result.temporalIntervals;
        const bool half=compare(right,E::fromDouble(.5))<=0;
        addGross(result,sums,value,sign,half);
    }
    for(const auto& event:events)result.maximumRootRefinements=std::max(result.maximumRootRefinements,event.refinements);
    E::validateEnvironment();return result;
}
} // query-private arithmetic

Gross integratePolygon(const std::vector<Point>& polygon,Point faceBegin,Point faceEnd,
    Point velocityBegin,Point velocityEnd,long double timeBegin,long double timeEnd,
    Point frameOrigin,bool rz,unsigned eventLimit,Reflection reflection)
{
    return integratePolygonImpl(polygon,faceBegin,faceEnd,velocityBegin,velocityEnd,
        timeBegin,timeEnd,frameOrigin,rz,eventLimit,reflection);
}
Gross integrateAnchoredPolygon(const std::vector<AnchoredPoint>& polygon,Point faceBegin,Point faceEnd,
    Point velocityBegin,Point velocityEnd,long double timeBegin,long double timeEnd,
    Point frameOrigin,bool rz,unsigned eventLimit,Reflection reflection)
{
    return integratePolygonImpl(polygon,faceBegin,faceEnd,velocityBegin,velocityEnd,
        timeBegin,timeEnd,frameOrigin,rz,eventLimit,reflection);
}
Gross integrateNormal(long double first,long double last,long double timeBegin,long double timeEnd,long double area)
{
    E::validateEnvironment();const auto dt=duration(timeBegin,timeEnd),a=scalar(first),b=scalar(last),measure=scalar(area);
    require(E::sign(measure)>=0,"negative native face area");if(E::sign(measure)==0)return {};
    Quadratic normal{{a,E::subtract(b,a),X{}}};std::vector<Root> events{constant(0),constant(.5),constant(1)};
    addRoots(normal,events,512);sortEvents(events);Gross result;result.temporalEvents=events.size();
    const R units=R(dt,E::fromDouble(1))*R(measure,E::fromDouble(1));std::array<PositiveSum,4> sums{};
    for(std::size_t i=1;i<events.size();++i){auto& left=events[i-1];auto& right=events[i];
        while(compare(left.hi,right.lo)>=0){refine(left);refine(right);}
        const auto sample=E::scalePow2(E::add(left.hi,right.lo),-1);const int sign=E::sign(at(normal,sample));if(sign==0)continue;
        Density branch;for(unsigned k=0;k<2;++k)branch.numerator[k]=sign>0?normal[k]:neg(normal[k]);
        branch.degree=E::sign(normal[1])==0?0:1;E::centerHomogeneous(branch.numerator,branch.denominator);
        const auto value=certifiedIntegral(branch,left,right,units,result);++result.temporalIntervals;
        addGross(result,sums,value,sign,compare(right,E::fromDouble(.5))<=0);
    }
    for(const auto& event:events)result.maximumRootRefinements=std::max(result.maximumRootRefinements,event.refinements);
    E::validateEnvironment();return result;
}
Bounds departureBounds(Point a,Point b,Point v0,Point v1,long double begin,long double end)
{
    E::validateEnvironment();const auto values=departure(a,b,v0,v1,begin,end);
    Bounds result{outward(values[0],false),outward(values[1],true),outward(values[2],false),outward(values[3],true)};
    E::validateEnvironment();return result;
}
void validateDepartureSupport(Point a,Point b,Point v0,Point v1,long double begin,long double end,Bounds domain)
{
    E::validateEnvironment();const auto values=departure(a,b,v0,v1,begin,end);
    const auto x0=scalar(domain.x0),x1=scalar(domain.x1),y0=scalar(domain.y0),y1=scalar(domain.y1);
    require(compare(x0,x1)<0&&compare(y0,y1)<0,"invalid original departure domain");
    require(compare(values[0],R(E::subtract(times(x0,2),x1),E::fromDouble(1)))>=0&&
        compare(values[1],R(E::subtract(times(x1,2),x0),E::fromDouble(1)))<=0&&
        compare(values[2],R(E::subtract(times(y0,2),y1),E::fromDouble(1)))>=0&&
        compare(values[3],R(E::subtract(times(y1,2),y0),E::fromDouble(1)))<=0,
        "actual departure requires more than one wall reflection");
    E::validateEnvironment();
}
}



