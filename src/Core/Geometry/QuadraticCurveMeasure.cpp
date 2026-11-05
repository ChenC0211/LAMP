#include "Core/Geometry/QuadraticCurveMeasure.h"
#include "Core/Geometry/QuadraticCurveDifference.h"
#include "Core/Geometry/UnitPolynomialRoots.h"
#include "Core/Geometry/Detail/QuadraticGraphExact.h"
#include "Core/Geometry/Detail/ExactLinearNativeFilter.h"
#include "Core/Mesh/CartesianMeasure.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <mutex>
#include <stdexcept>

namespace LAMP::Geometry::Detail
{
struct QuadraticDerivativeNativeCache
{
    mutable std::once_flag firstOnce,secondOnce;
    mutable std::array<ExactLinearNative::Filter,2> first;
    // Keep the ORIGINAL long-double rounding. Casting it to Point's double
    // belongs to the query environment, not to the cache creation time.
    mutable std::array<long double,2> second{};
    mutable bool secondAvailable=false;
    void prepareFirst(const QuadraticParametricData& data) const
    {
        std::call_once(firstOnce,[&]{
            first[0]=ExactLinearNative::Filter::prepare(data.coefficients[1].x,2*data.coefficients[2].x);
            first[1]=ExactLinearNative::Filter::prepare(data.coefficients[1].y,2*data.coefficients[2].y);
        });
    }
    void prepareSecond(const QuadraticParametricData& data) const
    {
        std::call_once(secondOnce,[&]{
            try{
                second={{QuadraticExact::rounded(2*data.coefficients[2].x),
                    QuadraticExact::rounded(2*data.coefficients[2].y)}};
                secondAvailable=true;
            }catch(const std::invalid_argument&){
                // Do not prepare D1, or replace an original D2 query error.
                // The ordinary exact path below retains its exception.
            }
        });
    }
    static const QuadraticDerivativeNativeCache* get(const QuadraticCurveMeasure::Curve& curve) noexcept
    {return curve.derivativeNative_.get();}
};
}

namespace LAMP::Geometry::QuadraticCurveMeasure
{
namespace
{
namespace E=Detail::QuadraticExact;
void require(bool ok,const char* message){if(!ok)throw std::invalid_argument(message);}
bool finite(Point p){return std::isfinite(p.x)&&std::isfinite(p.y);}
double coordinate(Point p,int axis){return axis?p.y:p.x;}
Point nativePoint(const E::Vector& p){const auto q=E::native(p);return {static_cast<double>(q.x),static_cast<double>(q.y)};}
E::Vector exact(Point p){return {E::scalar(p.x),E::scalar(p.y)};}
E::Scalar coordinate(const E::Vector& p,int axis){return axis?p.y:p.x;}
E::Polynomial polynomial(const Detail::QuadraticParametricData& data,int axis,double level=0)
{return {{coordinate(data.coefficients[0],axis)-E::scalar(level),coordinate(data.coefficients[1],axis),coordinate(data.coefficients[2],axis),0,0}};}

Parametric::PolynomialArc<3> arc(const Curve& curve)
{
    for(const auto p:curve.coefficients)require(finite(p),"[QuadraticCurveMeasure] nonfinite coefficient");
    Parametric::PolynomialArc<3> result;result.a=0;result.b=1;result.p=curve.coefficients;return result;
}
void validate(const Curve& curve)
{
    if(curve.usesExactBasis())require(curve.exactViewUnchanged(),"[QuadraticCurveMeasure] stale coefficient view on exact curve");
    else if(const auto* e=curve.endpointBasis())require(finite(e->first)&&finite(e->last)&&finite(e->bubble),
        "[QuadraticCurveMeasure] nonfinite endpoint basis");
    else (void)arc(curve);
}
void parameter(double s)
{require(std::isfinite(s)&&s>=0&&s<=1,"[QuadraticCurveMeasure] parameter outside [0,1]");}
Point value(const Curve& curve,double s)
{
    if(curve.usesExactBasis()){
        const auto& data=*curve.exactData();
        // Views of the SAME algebraic side events, never spatial snapping.
        if(s==0||s==1){const auto p=data.endpoints[s==0?0:1];return {static_cast<double>(p.x),static_cast<double>(p.y)};}
        return nativePoint(Detail::parametricPosition(data,E::scalar(s)));
    }
    if(const auto* e=curve.endpointBasis()){
        const double q=1-s;
        // The endpoint factors remain explicit throughout the interval.
        // In particular r_last=0 has no sum of power coefficients at s=1.
        if(s<=.5)return {q*e->first.x+s*(e->last.x+q*e->bubble.x),
            q*e->first.y+s*(e->last.y+q*e->bubble.y)};
        return {q*(e->first.x+s*e->bubble.x)+s*e->last.x,
            q*(e->first.y+s*e->bubble.y)+s*e->last.y};
    }
    return Parametric::localValue(curve.coefficients,s);
}
Point differentiated(const Curve& curve,double s,int order)
{
    if(order>=3)return {};
    if(curve.usesExactBasis()){
        // Public queries validate the diagnostic view first; integrate/bounds
        // validate once before their internal immutable derivative samples.
        if(const auto* cached=Detail::QuadraticDerivativeNativeCache::get(curve);
            cached&&Detail::ExactLinearNative::slopeFilterEnvironment()){
            if(order==2){cached->prepareSecond(*curve.exactData());
                if(cached->secondAvailable&&Detail::ExactLinearNative::slopeFilterEnvironment())
                    return {static_cast<double>(cached->second[0]),static_cast<double>(cached->second[1])};}
            if(order==1){cached->prepareFirst(*curve.exactData());long double x=0,y=0;
                if(cached->first[0].evaluate(s,x)&&cached->first[1].evaluate(s,y))
                    return {static_cast<double>(x),static_cast<double>(y)};}
        }
        return nativePoint(Detail::parametricDerivative(*curve.exactData(),E::scalar(s),order));
    }
    if(const auto* e=curve.endpointBasis()){
        if(order==2)return {-2*e->bubble.x,-2*e->bubble.y};
        return {(e->last.x-e->first.x)+(1-2*s)*e->bubble.x,
            (e->last.y-e->first.y)+(1-2*s)*e->bubble.y};
    }
    return Parametric::localValue(curve.coefficients,s,order);
}
std::vector<UnitPolynomialTerm> endpointTerms(const Curve::EndpointBasis& e,int axis,double level,bool firstDerivative=false)
{
    // Exact combinations of the SAME retained binary64 endpoint/bubble data.
    // Combining into three rounded power coefficients would lose an axis root.
    std::vector<UnitPolynomialTerm> terms(4);
    terms[0].factor=coordinate(e.first,axis);terms[1].factor=coordinate(e.last,axis);
    terms[2].factor=coordinate(e.bubble,axis);terms[3].factor=level;
    if(firstDerivative){
        terms[0].coefficients[0]=-1;terms[1].coefficients[0]=1;
        terms[2].coefficients[0]=1;terms[2].coefficients[1]=-2;
    }else{
        terms[0].coefficients[0]=1;terms[0].coefficients[1]=-1;
        terms[1].coefficients[1]=1;terms[2].coefficients[1]=1;terms[2].coefficients[2]=-1;
    }
    terms[3].coefficients[0]=-1;return terms;
}
std::vector<double> stationary(const Curve& curve,int axis)
{
    const auto roots=isolateUnitPolynomial(endpointTerms(*curve.endpointBasis(),axis,0,true));
    std::vector<double> result;
    for(const auto& root:roots.roots)if(root.point>0&&root.point<1)result.push_back(root.point);
    return result;
}
bool finitePositionEnvelope(const Curve& curve)
{
    const double limit=std::numeric_limits<double>::max()/16;
    const auto bounded=[limit](Point p){return finite(p)&&std::abs(p.x)<=limit&&std::abs(p.y)<=limit;};
    if(!curve.usesExactBasis()){
        // Use the retained endpoint authority, never its lossy power view.
        // For s in [0,1], these formulas and legacy Horner stay below 4*limit.
        if(const auto* e=curve.endpointBasis())return bounded(e->first)&&bounded(e->last)&&bounded(e->bubble);
        for(const auto p:curve.coefficients)if(!bounded(p))return false;
        return true;
    }
    // integrate() has already validated the fromGraph certified, unchanged
    // native views. The real rational coefficient is bounded by
    // 2*abs(view)+DBL_MIN, including underflow in the final native cast.
    // Thus all exact polynomial values stay far below DBL_MAX; exact
    // arithmetic itself has no fixed-width intermediate overflow.
    for(const auto p:curve.coefficients)if(!bounded(p))return false;
    const auto& data=*curve.exactData();
    for(const auto p:data.endpoints)if(!std::isfinite(p.x)||!std::isfinite(p.y)||
        std::abs(p.x)>limit||std::abs(p.y)>limit)return false;
    const double minimum=std::numeric_limits<double>::min();
    for(int axis=0;axis<2;++axis){
        const double a=std::abs(coordinate(curve.coefficients[0],axis)),
            b=std::abs(coordinate(curve.coefficients[1],axis)),
            c=std::abs(coordinate(curve.coefficients[2],axis));
        if(a==0&&b==0&&c==0){
            // Zero native views alone do not prove an exact zero polynomial.
            for(const auto& p:data.coefficients)if((axis?p.y:p.x)!=0)return false;
        }else{
            // Also preserve rounded()'s nonzero-underflow rejection: this
            // conservative dominance bound gives abs(p(s)) > 8*DBL_MIN.
            // Near-zero/cancelling coordinates retain the original value().
            if(!(a>=32*minimum&&a/8>b+c+minimum))return false;
        }
    }
    return true;
}
Integral gauss(const Curve& curve,double lo,double hi,bool rz)
{
    const auto& rule=Parametric::gaussEight();const auto half=(hi-lo)/2,mid=lo+half;
    Integral result;const bool positionCertified=!rz&&finitePositionEnvelope(curve);
    for(std::size_t i=0;i<rule.x.size();++i){
        const auto t=mid+half*rule.x[i];const auto p=positionCertified?Point{}:value(curve,t),d=differentiated(curve,t,1);
        require(finite(p)&&finite(d),"[QuadraticCurveMeasure] nonfinite quadrature sample");
        const long double length=half*rule.w[i]*std::hypot(d.x,d.y);
        require(!rz||p.x>=0,"[QuadraticCurveMeasure] negative radius");
        result.length+=length;result.physicalMeasure+=rz?CartesianMeasure::twoPi*p.x*length:length;
    }
    return result;
}
Integral adaptive(const Curve& curve,double lo,double hi,bool rz,long double tolerance,unsigned remaining)
{
    const auto mid=lo+(hi-lo)/2;
    const auto whole=gauss(curve,lo,hi,rz),a=gauss(curve,lo,mid,rz),b=gauss(curve,mid,hi,rz);
    Integral out{a.length+b.length,a.physicalMeasure+b.physicalMeasure,0};
    out.quadratureErrorEstimate=std::abs(out.physicalMeasure-whole.physicalMeasure);
    if(std::abs(out.length-whole.length)<=tolerance*out.length&&out.quadratureErrorEstimate<=tolerance*out.physicalMeasure)return out;
    require(remaining>0&&mid!=lo&&mid!=hi,"[QuadraticCurveMeasure] quadrature budget exhausted");
    const auto left=adaptive(curve,lo,mid,rz,tolerance,remaining-1),right=adaptive(curve,mid,hi,rz,tolerance,remaining-1);
    return {left.length+right.length,left.physicalMeasure+right.physicalMeasure,left.quadratureErrorEstimate+right.quadratureErrorEstimate};
}
// Authoritative coefficients of the SAME polynomial. Endpoint-basis powers
// are assembled exactly here; mutable diagnostic coefficients are unused.
std::array<E::Vector,3> measureCoefficients(const Curve& curve)
{
    if(curve.usesExactBasis())return curve.exactData()->coefficients;
    if(const auto* e=curve.endpointBasis()){
        const auto first=exact(e->first),last=exact(e->last),bubble=exact(e->bubble);
        return {{first,E::add(E::subtract(last,first),bubble),E::multiply(bubble,-1)}};
    }
    return {{exact(curve.coefficients[0]),exact(curve.coefficients[1]),exact(curve.coefficients[2])}};
}
E::Vector measureDerivative(const std::array<E::Vector,3>& c,const E::Scalar& s)
{return E::add(c[1],E::multiply(c[2],2*s));}
void validateMeasurePolynomial(const std::array<E::Vector,3>& c,bool rz)
{
    if(rz){
        require(c[0].x>=0&&c[0].x+c[1].x+c[2].x>=0,"[QuadraticCurveMeasure] radial curve crosses the axis");
        if(c[2].x!=0){const E::Scalar s=-c[1].x/(2*c[2].x);
            if(s>0&&s<1)require((c[2].x*s+c[1].x)*s+c[0].x>=0,"[QuadraticCurveMeasure] radial curve crosses the axis");}
    }
    if(c[2].x!=0||c[2].y!=0){
        const E::Scalar s=c[2].x!=0?-c[1].x/(2*c[2].x):-c[1].y/(2*c[2].y);
        if(s>=0&&s<=1){const auto d=measureDerivative(c,s);
            require(d.x!=0||d.y!=0,"[QuadraticCurveMeasure] nonregular material parameter");}
    }
}
struct NativeMeasureVector {long double x=0,y=0;};
NativeMeasureVector nativeMeasureVector(const E::Vector& p){return {E::rounded(p.x),E::rounded(p.y)};}
struct PairedMeasurePolynomial
{
    std::array<NativeMeasureVector,2> before,after,difference,sum;
    std::array<long double,3> radiusDifference{},radiusSum{};
};
PairedMeasurePolynomial pairedMeasurePolynomial(const Curve& before,const Curve& after,bool rz)
{
    const auto a=measureCoefficients(before),b=measureCoefficients(after);
    validateMeasurePolynomial(a,rz);validateMeasurePolynomial(b,rz);
    PairedMeasurePolynomial p;
    for(unsigned i=0;i<2;++i){
        const auto da=E::multiply(a[i+1],i+1),db=E::multiply(b[i+1],i+1);
        p.before[i]=nativeMeasureVector(da);p.after[i]=nativeMeasureVector(db);
        // Subtract/add the SAME exact coefficients once before rounding.
        p.difference[i]=nativeMeasureVector(E::subtract(db,da));p.sum[i]=nativeMeasureVector(E::add(db,da));
    }
    for(unsigned i=0;i<3;++i){p.radiusDifference[i]=E::rounded(b[i].x-a[i].x);p.radiusSum[i]=E::rounded(b[i].x+a[i].x);}
    return p;
}
NativeMeasureVector pairedDerivative(const std::array<NativeMeasureVector,2>& c,long double s)
{return {std::fma(c[1].x,s,c[0].x),std::fma(c[1].y,s,c[0].y)};}
long double pairedRadius(const std::array<long double,3>& c,long double s)
{return std::fma(std::fma(c[2],s,c[1]),s,c[0]);}
long double pairedMeasureDensity(const PairedMeasurePolynomial& p,double parameter,bool rz)
{
    const auto a=pairedDerivative(p.before,parameter),b=pairedDerivative(p.after,parameter);
    const auto speed0=std::hypot(a.x,a.y),speed1=std::hypot(b.x,b.y),speedSum=speed0+speed1;
    require(std::isfinite(speedSum),"[QuadraticCurveMeasure] nonfinite paired material speed");
    if(speedSum==0)return 0; // Both constant curves have zero physical measure.
    const auto difference=pairedDerivative(p.difference,parameter),sum=pairedDerivative(p.sum,parameter);
    const auto deltaSpeed=std::fma(difference.x,sum.x,difference.y*sum.y)/speedSum;
    long double delta=deltaSpeed;
    if(rz){
        const auto radialDelta=pairedRadius(p.radiusDifference,parameter),radialSum=pairedRadius(p.radiusSum,parameter);
        delta=(CartesianMeasure::twoPi/2)*std::fma(radialDelta,speedSum,radialSum*deltaSpeed);
    }
    require(std::isfinite(delta),"[QuadraticCurveMeasure] nonfinite paired measure difference");return delta;
}
struct PairedMeasureIntegral {long double delta=0,absoluteIntegrand=0,error=0;};
PairedMeasureIntegral pairedGauss(const PairedMeasurePolynomial& p,double lo,double hi,bool rz)
{
    const auto& rule=Parametric::gaussEight();const auto half=(hi-lo)/2,mid=lo+half;
    PairedMeasureIntegral result;long double correction=0;
    for(std::size_t i=0;i<rule.x.size();++i){
        const auto density=pairedMeasureDensity(p,mid+half*rule.x[i],rz);
        const long double term=half*rule.w[i]*density,adjusted=term-correction,next=result.delta+adjusted;
        correction=(next-result.delta)-adjusted;result.delta=next;result.absoluteIntegrand+=std::abs(term);
    }
    require(std::isfinite(result.delta)&&std::isfinite(result.absoluteIntegrand),"[QuadraticCurveMeasure] unrepresentable paired quadrature");return result;
}
PairedMeasureIntegral pairedAdaptive(const PairedMeasurePolynomial& p,double lo,double hi,bool rz,long double tolerance,unsigned remaining)
{
    const auto mid=lo+(hi-lo)/2;
    const auto whole=pairedGauss(p,lo,hi,rz),a=pairedGauss(p,lo,mid,rz),b=pairedGauss(p,mid,hi,rz);
    PairedMeasureIntegral out{a.delta+b.delta,a.absoluteIntegrand+b.absoluteIntegrand,0};
    out.error=std::abs(out.delta-whole.delta);
    if(out.error<=tolerance*std::max(std::abs(out.delta),out.absoluteIntegrand))return out;
    require(remaining>0&&mid!=lo&&mid!=hi,"[QuadraticCurveMeasure] paired quadrature budget exhausted");
    const auto left=pairedAdaptive(p,lo,mid,rz,tolerance,remaining-1),right=pairedAdaptive(p,mid,hi,rz,tolerance,remaining-1);
    return {left.delta+right.delta,left.absoluteIntegrand+right.absoluteIntegrand,left.error+right.error};
}
}
Curve Curve::fromGraph(std::shared_ptr<const Detail::QuadraticParametricData> data)
{
    require(static_cast<bool>(data),"[QuadraticCurveMeasure] missing exact polynomial data");
    Curve result;result.exact_=std::move(data);
    for(unsigned i=0;i<3;++i){result.coefficients[i]=nativePoint(result.exact_->coefficients[i]);
        require(finite(result.coefficients[i]),"[QuadraticCurveMeasure] nonfinite exact coefficient view");}
    result.exactView_=result.coefficients;
    result.derivativeNative_=std::make_shared<Detail::QuadraticDerivativeNativeCache>();return result;
}
bool Curve::exactViewUnchanged() const noexcept
{
    for(unsigned i=0;i<3;++i)if(coefficients[i].x!=exactView_[i].x||coefficients[i].y!=exactView_[i].y)return false;
    return true;
}
bool samePolynomial(const Curve& a,const Curve& b)
{
    validate(a);validate(b);
    if(a.usesExactBasis()!=b.usesExactBasis()||a.usesEndpointBasis()!=b.usesEndpointBasis())return false;
    if(a.usesExactBasis()){
        if(a.exactData()==b.exactData())return true;const auto& x=*a.exactData();const auto& y=*b.exactData();
        for(unsigned i=0;i<3;++i)if(x.coefficients[i].x!=y.coefficients[i].x||x.coefficients[i].y!=y.coefficients[i].y)return false;
        for(unsigned i=0;i<2;++i)if(x.endpoints[i].x!=y.endpoints[i].x||x.endpoints[i].y!=y.endpoints[i].y||x.endpointBoundary[i]!=y.endpointBoundary[i])return false;
        return true;
    }
    if(const auto* x=a.endpointBasis()){const auto& y=*b.endpointBasis();return x->first.x==y.first.x&&x->first.y==y.first.y&&
        x->last.x==y.last.x&&x->last.y==y.last.y&&x->bubble.x==y.bubble.x&&x->bubble.y==y.bubble.y;}
    for(unsigned i=0;i<3;++i)if(a.coefficients[i].x!=b.coefficients[i].x||a.coefficients[i].y!=b.coefficients[i].y)return false;
    return true;
}
Curve affineMapped(const Curve& source,Point sourceAnchor,Point imageAnchor,const std::array<double,4>& matrix)
{
    validate(source);require(finite(sourceAnchor)&&finite(imageAnchor),"[QuadraticCurveMeasure] nonfinite affine anchor");
    for(auto v:matrix)require(std::isfinite(v),"[QuadraticCurveMeasure] nonfinite affine derivative");
    if(matrix==std::array<double,4>{{1,0,0,1}}&&sourceAnchor.x==imageAnchor.x&&sourceAnchor.y==imageAnchor.y)return source;
    const auto mapped=[&](Point p){const auto d=p-sourceAnchor;return Point{imageAnchor.x+matrix[0]*d.x+matrix[1]*d.y,imageAnchor.y+matrix[2]*d.x+matrix[3]*d.y};};
    if(!source.usesExactBasis())return Curve::fromSamples(mapped(value(source,0)),mapped(value(source,.5)),mapped(value(source,1)));
    auto data=std::make_shared<Detail::QuadraticParametricData>(*source.exactData());
    const auto a=exact(sourceAnchor),b=exact(imageAnchor);
    const std::array<E::Scalar,4> f{{E::scalar(matrix[0]),E::scalar(matrix[1]),E::scalar(matrix[2]),E::scalar(matrix[3])}};
    const auto linear=[&](const E::Vector& p){return E::Vector{f[0]*p.x+f[1]*p.y,f[2]*p.x+f[3]*p.y};};
    for(unsigned i=0;i<3;++i)data->coefficients[i]=i?linear(data->coefficients[i]):E::add(b,linear(E::subtract(data->coefficients[i],a)));
    for(unsigned i=0;i<2;++i)data->endpoints[i]=E::native(E::add(b,linear(E::subtract(E::vector(source.exactData()->endpoints[i]),a))));
    data->endpointBoundary={{0,0}};return Curve::fromGraph(std::move(data));
}
Curve displacedSamples(const Curve& source,Point first,Point middle,Point last)
{
    validate(source);require(finite(first)&&finite(middle)&&finite(last),"[QuadraticCurveMeasure] nonfinite displacement sample");
    if(first.x==0&&first.y==0&&middle.x==0&&middle.y==0&&last.x==0&&last.y==0)return source;
    if(!source.usesExactBasis())return Curve::fromSamples(value(source,0)+first,value(source,.5)+middle,value(source,1)+last);
    auto data=std::make_shared<Detail::QuadraticParametricData>(*source.exactData());
    const auto a=exact(first),m=exact(middle),b=exact(last);
    const std::array<E::Vector,3> delta{{a,E::subtract(E::subtract(E::multiply(m,4),E::multiply(a,3)),b),
        E::subtract(E::multiply(E::add(a,b),2),E::multiply(m,4))}};
    for(unsigned i=0;i<3;++i)data->coefficients[i]=E::add(data->coefficients[i],delta[i]);
    data->endpoints[0]=E::native(E::add(E::vector(source.exactData()->endpoints[0]),a));
    data->endpoints[1]=E::native(E::add(E::vector(source.exactData()->endpoints[1]),b));
    data->endpointBoundary={{0,0}};return Curve::fromGraph(std::move(data));
}
Curve Curve::fromSamples(Point first,Point middle,Point last)
{
    require(finite(first)&&finite(middle)&&finite(last),"[QuadraticCurveMeasure] nonfinite original sample");
    Curve result;const Point bubble{2*((middle.x-first.x)+(middle.x-last.x)),
        2*((middle.y-first.y)+(middle.y-last.y))};
    require(finite(bubble),"[QuadraticCurveMeasure] unrepresentable endpoint bubble");
    result.endpoints_=EndpointBasis{first,last,bubble};
    result.coefficients[0]=first;
    result.coefficients[1]={last.x-first.x+bubble.x,last.y-first.y+bubble.y};
    result.coefficients[2]={-bubble.x,-bubble.y};
    return result;
}
Point position(const Curve& curve,double s)
{parameter(s);validate(curve);const auto p=value(curve,s);require(finite(p),"[QuadraticCurveMeasure] nonfinite position");return p;}
Point derivative(const Curve& curve,double s,int order)
{
    parameter(s);validate(curve);require(order>0,"[QuadraticCurveMeasure] derivative order must be positive");
    const auto d=differentiated(curve,s,order);require(finite(d),"[QuadraticCurveMeasure] nonfinite derivative");return d;
}
Bounds bounds(const Curve& curve)
{
    validate(curve);
    if(curve.usesExactBasis()){
        const auto a=value(curve,0),b=value(curve,1);Bounds result{{std::min(a.x,b.x),std::min(a.y,b.y)},
            {std::max(a.x,b.x),std::max(a.y,b.y)}};
        const auto& data=*curve.exactData();
        for(int axis=0;axis<2;++axis){const E::Polynomial d{{coordinate(data.coefficients[1],axis),2*coordinate(data.coefficients[2],axis),0,0,0}};
            for(const auto& t:E::roots(d,0,1))if(t>0&&t<1){const auto p=nativePoint(Detail::parametricPosition(data,t));
                if(axis){result.lower.y=std::min(result.lower.y,p.y);result.upper.y=std::max(result.upper.y,p.y);}
                else{result.lower.x=std::min(result.lower.x,p.x);result.upper.x=std::max(result.upper.x,p.x);}}}
        return result;
    }
    if(curve.usesEndpointBasis()){
        const auto a=value(curve,0),b=value(curve,1);Bounds result{{std::min(a.x,b.x),std::min(a.y,b.y)},
            {std::max(a.x,b.x),std::max(a.y,b.y)}};
        for(int axis=0;axis<2;++axis)for(double s:stationary(curve,axis)){
            const auto p=value(curve,s);require(finite(p),"[QuadraticCurveMeasure] nonfinite stationary point");
            if(axis){result.lower.y=std::min(result.lower.y,p.y);result.upper.y=std::max(result.upper.y,p.y);}
            else{result.lower.x=std::min(result.lower.x,p.x);result.upper.x=std::max(result.upper.x,p.x);}
        }
        return result;
    }
    const auto a=arc(curve);const auto x=Parametric::range(a.coordinate(0)),y=Parametric::range(a.coordinate(1));
    require(std::isfinite(x[0])&&std::isfinite(x[1])&&std::isfinite(y[0])&&std::isfinite(y[1]),"[QuadraticCurveMeasure] nonfinite bounds");
    return {{x[0],y[0]},{x[1],y[1]}};
}
std::vector<double> coordinateRoots(const Curve& curve,int axis,double level,double lo,double hi)
{
    validate(curve);parameter(lo);parameter(hi);
    require((axis==0||axis==1)&&std::isfinite(level)&&lo<=hi,"[QuadraticCurveMeasure] invalid coordinate-root query");
    std::vector<double> result;
    if(curve.usesExactBasis()){
        const auto& data=*curve.exactData();const auto p=polynomial(data,axis,level);
        for(const auto& t:E::roots(p,E::scalar(lo),E::scalar(hi),{0,1}))result.push_back(static_cast<double>(E::rounded(t)));
        const auto mask=axis?(QuadraticSurfaceMeasure::Bottom|QuadraticSurfaceMeasure::Top):(QuadraticSurfaceMeasure::Left|QuadraticSurfaceMeasure::Right);
        for(unsigned i=0;i<2;++i)if((data.endpointBoundary[i]&mask)&&
            (axis?data.endpoints[i].y:data.endpoints[i].x)==level&&i>=lo&&i<=hi)result.push_back(static_cast<double>(i));
    }else if(const auto* e=curve.endpointBasis()){
        const auto found=isolateUnitPolynomial(endpointTerms(*e,axis,level));
        for(const auto& root:found.roots)if(root.point>=lo&&root.point<=hi)result.push_back(root.point);
    }else{
        auto p=arc(curve).coordinate(axis);p[0]-=level;result=Parametric::roots(p,lo,hi);
        if(Parametric::polynomial(p,lo)==0)result.push_back(lo);
        if(Parametric::polynomial(p,hi)==0)result.push_back(hi);
    }
    std::sort(result.begin(),result.end());result.erase(std::unique(result.begin(),result.end()),result.end());return result;
}
std::vector<std::array<double,2>> intervals(const Curve& curve,const Rectangle& box,unsigned ownership)
{
    (void)PlanarCellCut::capacity(box,false);validate(curve);
    require((ownership&~QuadraticSurfaceMeasure::All)==0,"[QuadraticCurveMeasure] invalid side ownership");
    if(curve.usesExactBasis()){
        const auto& data=*curve.exactData();
        for(int axis=0;axis<2;++axis)if(coordinate(data.coefficients[1],axis)==0&&coordinate(data.coefficients[2],axis)==0){
            const auto p=coordinate(data.coefficients[0],axis);const auto lower=axis?box.lower.y:box.lower.x,upper=axis?box.upper.y:box.upper.x;
            const auto lowMask=axis?QuadraticSurfaceMeasure::Bottom:QuadraticSurfaceMeasure::Left;
            const auto highMask=axis?QuadraticSurfaceMeasure::Top:QuadraticSurfaceMeasure::Right;
            if((p==E::scalar(lower)&&!(ownership&lowMask))||(p==E::scalar(upper)&&!(ownership&highMask)))return {};
        }
    }else if(!curve.usesEndpointBasis()){
        const auto a=arc(curve);const auto& p=a.p;
        if(p[1].x==0&&p[2].x==0&&((p[0].x==box.lower.x&&!(ownership&QuadraticSurfaceMeasure::Left))||
            (p[0].x==box.upper.x&&!(ownership&QuadraticSurfaceMeasure::Right))))return {};
        if(p[1].y==0&&p[2].y==0&&((p[0].y==box.lower.y&&!(ownership&QuadraticSurfaceMeasure::Bottom))||
            (p[0].y==box.upper.y&&!(ownership&QuadraticSurfaceMeasure::Top))))return {};
    }else{
        const auto& e=*curve.endpointBasis();
        for(int axis=0;axis<2;++axis)if(coordinate(e.first,axis)==coordinate(e.last,axis)&&coordinate(e.bubble,axis)==0){
            const auto p=coordinate(e.first,axis),lower=coordinate(Point{box.lower.x,box.lower.y},axis),upper=coordinate(Point{box.upper.x,box.upper.y},axis);
            const auto lowMask=axis?QuadraticSurfaceMeasure::Bottom:QuadraticSurfaceMeasure::Left;
            const auto highMask=axis?QuadraticSurfaceMeasure::Top:QuadraticSurfaceMeasure::Right;
            if((p==lower&&!(ownership&lowMask))||(p==upper&&!(ownership&highMask)))return {};
        }
    }
    std::vector<std::array<double,2>> split;
    if(!curve.usesEndpointBasis()&&!curve.usesExactBasis())split=Parametric::rectangleIntervals(arc(curve),{box.lower.x,box.upper.x,box.lower.y,box.upper.y});
    else{
        std::vector<double> cuts{0,1};
        for(int axis=0;axis<2;++axis)for(double edge:axis?std::array<double,2>{{box.lower.y,box.upper.y}}:std::array<double,2>{{box.lower.x,box.upper.x}}){
            const auto found=coordinateRoots(curve,axis,edge);cuts.insert(cuts.end(),found.begin(),found.end());}
        std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
        for(std::size_t i=1;i<cuts.size();++i){
            const auto t=cuts[i-1]+(cuts[i]-cuts[i-1])/2;
            if(curve.usesExactBasis()){
                const auto p=Detail::parametricPosition(*curve.exactData(),E::scalar(t));
                if(p.x>=E::scalar(box.lower.x)&&p.x<=E::scalar(box.upper.x)&&p.y>=E::scalar(box.lower.y)&&p.y<=E::scalar(box.upper.y))split.push_back({cuts[i-1],cuts[i]});
            }else{const auto p=value(curve,t);
                if(p.x>=box.lower.x&&p.x<=box.upper.x&&p.y>=box.lower.y&&p.y<=box.upper.y)split.push_back({cuts[i-1],cuts[i]});}}
    }
    std::vector<std::array<double,2>> result;
    for(const auto interval:split){
        if(!result.empty()&&result.back()[1]==interval[0])result.back()[1]=interval[1];
        else result.push_back(interval);}
    return result;
}
Integral integrate(const Curve& curve,double lo,double hi,bool rz,long double tolerance,unsigned depth)
{
    validate(curve);
    require(std::isfinite(lo)&&std::isfinite(hi)&&lo>=0&&hi<=1&&lo<=hi&&
        std::isfinite(tolerance)&&tolerance>0&&depth<=32,"[QuadraticCurveMeasure] invalid quadrature input");
    if(lo==hi)return {};
    if(rz)require(bounds(curve).lower.x>=0,"[QuadraticCurveMeasure] radial curve crosses the axis");
    if(curve.usesExactBasis()){
        const auto& c=curve.exactData()->coefficients;
        if(c[1].x==0&&c[1].y==0&&c[2].x==0&&c[2].y==0)return {};
        const auto second=E::multiply(c[2],2);
        if(second.x!=0||second.y!=0){const E::Scalar t=second.x!=0?-c[1].x/second.x:-c[1].y/second.y;
            if(t>=E::scalar(lo)&&t<=E::scalar(hi)){const auto d=Detail::parametricDerivative(*curve.exactData(),t,1);
                require(d.x!=0||d.y!=0,"[QuadraticCurveMeasure] nonregular exact material parameter");}}
    }else if(!curve.usesEndpointBasis()){
        const auto p=curve.coefficients;
        if(p[1].x==0&&p[1].y==0&&p[2].x==0&&p[2].y==0)return {};
        const auto square=Parametric::dot(p[2],p[2]);
        if(square>0){const auto stationary=-Parametric::dot(p[1],p[2])/(2*square);
            if(stationary>=lo&&stationary<=hi){const auto d=differentiated(curve,stationary,1);require(d.x!=0||d.y!=0,"[QuadraticCurveMeasure] nonregular material parameter");}}
    }else{
        const auto& e=*curve.endpointBasis();
        if(e.first.x==e.last.x&&e.first.y==e.last.y&&e.bubble.x==0&&e.bubble.y==0)return {};
        const auto first=differentiated(curve,0,1),second=differentiated(curve,0,2);
        require(finite(first)&&finite(second),"[QuadraticCurveMeasure] nonfinite material derivative");
        const long double square=static_cast<long double>(second.x)*second.x+static_cast<long double>(second.y)*second.y;
        if(square>0){const auto s=-(static_cast<long double>(first.x)*second.x+static_cast<long double>(first.y)*second.y)/square;
            if(s>=lo&&s<=hi){const auto d=differentiated(curve,static_cast<double>(s),1);require(d.x!=0||d.y!=0,"[QuadraticCurveMeasure] nonregular material parameter");}}
    }
    return adaptive(curve,lo,hi,rz,tolerance,depth);
}
PhysicalMeasureDifference physicalMeasureDifference(const Curve& before,const Curve& after,bool rz,long double tolerance,unsigned depth)
{
    validate(before);validate(after);
    require(std::isfinite(tolerance)&&tolerance>0&&depth<=32,"[QuadraticCurveMeasure] invalid paired quadrature input");
    const auto polynomial=pairedMeasurePolynomial(before,after,rz);
    const auto result=pairedAdaptive(polynomial,0,1,rz,tolerance,depth);
    return {result.delta,result.error};
}
}