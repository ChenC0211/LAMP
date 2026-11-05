#pragma once
// Private same-polynomial numerics. No mesh, state, phase inventory or solver owner.
#include "Core/Geometry/QuadraticCellCut.h"
#ifndef BOOST_MP_STANDALONE
#define BOOST_MP_STANDALONE
#endif
#include <boost/multiprecision/cpp_int.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace LAMP::Geometry::Detail {
namespace QuadraticExact {
using Integer=boost::multiprecision::cpp_int;
using Scalar=boost::multiprecision::number<boost::multiprecision::cpp_rational_backend,boost::multiprecision::et_off>;
struct Vector { Scalar x=0,y=0; };
using Polynomial=std::array<Scalar,5>;
inline void require(bool ok,const char* message) {if(!ok)throw std::invalid_argument(message);}
inline Scalar scalar(long double v) {
    require(std::isfinite(v),"[QuadraticExact] nonfinite native scalar");if(v==0)return 0;
    int e=0;auto f=std::frexp(std::abs(v),&e);Integer bits=0;
    for(int i=0;i<std::numeric_limits<long double>::digits;++i){f*=2;bits<<=1;if(f>=1){++bits;f-=1;}}
    require(f==0,"[QuadraticExact] unretained native bits");Scalar r(bits);e-=std::numeric_limits<long double>::digits;
    if(e>=0)r*=Scalar(Integer(1)<<e);else r/=Scalar(Integer(1)<<-e);return v<0?-r:r;
}
inline bool even(long double v) {int e=0;auto f=std::frexp(v,&e);for(int i=1;i<std::numeric_limits<long double>::digits;++i)f=std::fmod(f*2,1.L);return std::fmod(f*2,2.L)==0;}
inline long double rounded(const Scalar& x) {
    if(x==0)return 0;const bool negative=x<0;const Scalar v=negative?-x:x;
    // For normal values away from the overflow binade, round the original
    // retained ratio by one integer quotient/remainder. This proves the same
    // nearest-even native view without constructing rational midpoint cells.
    // The original path below retains all subnormal/overflow boundary behavior.
    const Integer numeratorValue=numerator(v),denominatorValue=denominator(v);
    auto exponent=static_cast<long long>(boost::multiprecision::msb(numeratorValue))-
        static_cast<long long>(boost::multiprecision::msb(denominatorValue));
    if(exponent>=std::numeric_limits<long double>::min_exponent-1LL&&
       exponent<std::numeric_limits<long double>::max_exponent-1LL){
        if(exponent>=0){if(numeratorValue<(denominatorValue<<exponent))--exponent;}
        else if((numeratorValue<<(-exponent))<denominatorValue)--exponent;
        if(exponent>=std::numeric_limits<long double>::min_exponent-1LL){
            Integer n=numeratorValue,d=denominatorValue,significand,remainder;
            const auto shift=std::numeric_limits<long double>::digits-1LL-exponent;
            if(shift>=0)n<<=shift;else d<<=-shift;
            boost::multiprecision::divide_qr(n,d,significand,remainder);
            const Integer twice=remainder<<1;
            if(twice>d||(twice==d&&boost::multiprecision::bit_test(significand,0)))++significand;
            const auto result=std::scalbn(significand.convert_to<long double>(),static_cast<int>(-shift));
            require(std::isfinite(result)&&result>=std::numeric_limits<long double>::min(),
                "[QuadraticExact] normal integer rounding left its certified binade");
            return negative?-result:result;
        }
    }
    auto c=v.convert_to<long double>();if(c==0)c=std::numeric_limits<long double>::denorm_min();
    if(!std::isfinite(c))c=std::numeric_limits<long double>::max();
    for(int k=0;k<16;++k){const auto l=std::nextafter(c,0.L),u=std::nextafter(c,std::numeric_limits<long double>::infinity());
        const auto lm=(scalar(l)+scalar(c))/2;const bool e=even(c);
        if(v<lm||(v==lm&&!e)){require(l>0,"[QuadraticExact] native nonzero underflow");c=l;continue;}
        const Scalar um=std::isfinite(u)?(scalar(c)+scalar(u))/2:scalar(c)+scalar(std::scalbn(1.L,std::numeric_limits<long double>::max_exponent-std::numeric_limits<long double>::digits-1));
        if(v>um||(v==um&&!e)){require(std::isfinite(u),"[QuadraticExact] native overflow");c=u;continue;}
        return negative?-c:c;
    }throw std::invalid_argument("[QuadraticExact] native rounding seed outside certified neighborhood");
}
inline Vector vector(PlanarCellCut::Point p) {return {scalar(p.x),scalar(p.y)};}
inline PlanarCellCut::Point native(const Vector& p) {return {rounded(p.x),rounded(p.y)};}
inline Vector add(const Vector&a,const Vector&b){return {a.x+b.x,a.y+b.y};}
inline Vector subtract(const Vector&a,const Vector&b){return {a.x-b.x,a.y-b.y};}
inline Vector multiply(const Vector&a,const Scalar&s){return {a.x*s,a.y*s};}
inline Scalar dot(const Vector&a,const Vector&b){return a.x*b.x+a.y*b.y;}
inline Scalar value(const Polynomial&p,const Scalar&t){Scalar r=p[4];for(int i=3;i>=0;--i)r=r*t+p[i];return r;}
inline Polynomial product(const Polynomial&a,const Polynomial&b){Polynomial r{};for(int i=0;i<5;++i)for(int j=0;i+j<5;++j)r[i+j]+=a[i]*b[j];return r;}
inline Scalar integral(const Polynomial&p,const Scalar&a,const Scalar&b){Scalar x=1,y=1,r=0;for(int i=0;i<5;++i){x*=a;y*=b;r+=p[i]*(y-x)/(i+1);}return r;}
// An explicitly bounded dyadic approximation of the SAME real radical. Its
// precision exceeds the native root model; integer products never overflow.
inline Scalar radical(const Scalar& x) {
    require(x>=0,"[QuadraticExact] negative radical");if(x==0)return 0;
    Integer n=numerator(x),d=denominator(x);const int exponent=static_cast<int>(boost::multiprecision::msb(n))-static_cast<int>(boost::multiprecision::msb(d));
    const int shift=std::numeric_limits<long double>::digits+48-exponent/2;
    if(shift>=0)n<<=2*shift;else d<<=-2*shift;
    const Integer q=n/d;const Integer r=boost::multiprecision::sqrt(q);Scalar result(r);
    if(shift>=0)result/=Scalar(Integer(1)<<shift);else result*=Scalar(Integer(1)<<-shift);return result;
}
inline std::vector<Scalar> roots(const Polynomial&p,const Scalar&lo,const Scalar&hi,const std::vector<Scalar>& known={}) {
    std::vector<Scalar> r;const auto put=[&](const Scalar&t){if(t>=lo&&t<=hi)r.push_back(t);};
    if(p[2]==0){if(p[1]!=0)put(-p[0]/p[1]);return r;}
    // Exact factory events factor the side polynomial before any radical.
    for(const auto&t:known)if(value(p,t)==0){put(t);put(-p[1]/p[2]-t);std::sort(r.begin(),r.end());r.erase(std::unique(r.begin(),r.end()),r.end());return r;}
    const Scalar disc=p[1]*p[1]-4*p[2]*p[0];if(disc<0)return r;
    if(disc==0)put(-p[1]/(2*p[2]));else {const Scalar d=radical(disc);const Scalar q=-(p[1]+(p[1]>=0?d:-d))/2;if(q==0)put(-p[1]/(2*p[2]));else {put(q/p[2]);put(p[0]/q);}}
    std::sort(r.begin(),r.end());r.erase(std::unique(r.begin(),r.end()),r.end());return r;
}
}
// Opaque proof produced only after the original fixed-port geometry gates.
struct QuadraticFixedPortCertificate;
struct QuadraticGraphExact {
    using R=QuadraticExact::Scalar;using V=QuadraticExact::Vector;
    V origin,tangent,normal;std::array<R,3> q{};R determinant=1;
    QuadraticCellCut::Shape nativeShape;long double nativeOffset=0;
    std::vector<R> portParameters;std::vector<PlanarCellCut::Point> ports;
    std::shared_ptr<const QuadraticFixedPortCertificate> fixedPortCertificate;
};
struct QuadraticParametricData {
    std::array<QuadraticExact::Vector,3> coefficients{};
    // These are original ports or actual side-root event views. An active
    // coordinate is the boundary constraint, never a near-coordinate snap.
    std::array<PlanarCellCut::Point,2> endpoints{};
    std::array<unsigned,2> endpointBoundary{};
};
inline QuadraticExact::Vector parametricPosition(const QuadraticParametricData&d,const QuadraticExact::Scalar&s) {
    using namespace QuadraticExact;return add(d.coefficients[0],multiply(add(d.coefficients[1],multiply(d.coefficients[2],s)),s));
}
inline QuadraticExact::Vector parametricDerivative(const QuadraticParametricData&d,const QuadraticExact::Scalar&s,int order=1) {
    using namespace QuadraticExact;if(order==0)return parametricPosition(d,s);if(order==1)return add(d.coefficients[1],multiply(d.coefficients[2],2*s));if(order==2)return multiply(d.coefficients[2],2);return {};
}
inline QuadraticExact::Vector graphPosition(const QuadraticGraphExact&g,const QuadraticExact::Scalar&t) {
    using namespace QuadraticExact;return add(g.origin,add(multiply(g.tangent,t),multiply(g.normal,(g.q[2]*t+g.q[1])*t+g.q[0])));
}
inline QuadraticExact::Scalar graphSlope(const QuadraticGraphExact&g,const QuadraticExact::Scalar&t) {return g.q[1]+2*g.q[2]*t;}
inline QuadraticExact::Polynomial coordinatePolynomial(const QuadraticGraphExact&g,int axis,long double border) {
    using namespace QuadraticExact;const auto n=axis?g.normal.y:g.normal.x,t=axis?g.tangent.y:g.tangent.x,o=axis?g.origin.y:g.origin.x;
    return {{o+n*g.q[0]-scalar(border),t+n*g.q[1],n*g.q[2],0,0}};
}
inline void validateGraph(const QuadraticCellCut::Reconstruction&g) {
    if(!g.exact)return;const auto&s=g.exact->nativeShape;const auto&t=g.shape;
    QuadraticExact::require(s.axis==t.axis&&s.slope==t.slope&&s.second==t.second&&s.frameNormal.has_value()==t.frameNormal.has_value()&&
        (!s.frameNormal||(s.frameNormal->x==t.frameNormal->x&&s.frameNormal->y==t.frameNormal->y))&&g.offset==g.exact->nativeOffset,
        "[QuadraticExact] stale coefficient/frame view on encoded graph");
}
struct GraphEvent {QuadraticExact::Scalar parameter;unsigned boundary=0;PlanarCellCut::Point position;};
struct GraphArc {GraphEvent first,last;};
inline std::shared_ptr<const QuadraticParametricData> graphArcMaterialData(const QuadraticGraphExact&g,const GraphArc&arc) {
    using namespace QuadraticExact;const auto&first=arc.first;const auto&last=arc.last;
    require(first.parameter<last.parameter,"[QuadraticCellCut] material interval has no original side-root identities");
    const Scalar lo=first.parameter,width=last.parameter-lo;
    auto d=std::make_shared<QuadraticParametricData>();d->coefficients[0]=graphPosition(g,lo);
    d->coefficients[1]=multiply(add(g.tangent,multiply(g.normal,g.q[1]+2*g.q[2]*lo)),width);
    d->coefficients[2]=multiply(g.normal,g.q[2]*width*width);d->endpoints={first.position,last.position};d->endpointBoundary={first.boundary,last.boundary};return d;
}
inline std::vector<GraphEvent> graphBoundaryEvents(const QuadraticGraphExact&g,const PlanarCellCut::Rectangle&b,unsigned ownership=15) {
    using namespace QuadraticExact;Scalar lo,hi;bool init=false;
    for(auto x:{b.lower.x,b.upper.x})for(auto y:{b.lower.y,b.upper.y}){const Vector d=subtract(vector({x,y}),g.origin);const Scalar t=dot(d,g.tangent)/g.determinant;if(!init){lo=hi=t;init=true;}else{lo=std::min(lo,t);hi=std::max(hi,t);}}
    std::vector<GraphEvent> e{{lo,0,{}},{hi,0,{}}};unsigned coincident=0;
    for(int i=0;i<4;++i){const int axis=i/2;const auto border=axis?(i%2?b.upper.y:b.lower.y):(i%2?b.upper.x:b.lower.x);const auto p=coordinatePolynomial(g,axis,border);const unsigned mask=1u<<i;
        if(p[0]==0&&p[1]==0&&p[2]==0){coincident|=mask;continue;}
        for(const auto&t:roots(p,lo,hi,g.portParameters))e.push_back({t,mask,{}});}
    if(coincident&~ownership)return {};
    std::sort(e.begin(),e.end(),[](const auto&a,const auto&b){return a.parameter<b.parameter;});std::vector<GraphEvent> cuts;
    for(const auto&a:e){if(!cuts.empty()&&cuts.back().parameter==a.parameter)cuts.back().boundary|=a.boundary;else cuts.push_back(a);}
    for(auto&a:cuts){a.boundary|=coincident;a.position=native(graphPosition(g,a.parameter));
        // Native representation of an actual algebraic side identity.
        if(a.boundary&1)a.position.x=b.lower.x;if(a.boundary&2)a.position.x=b.upper.x;if(a.boundary&4)a.position.y=b.lower.y;if(a.boundary&8)a.position.y=b.upper.y;
        for(std::size_t k=0;k<g.portParameters.size();++k)if(a.parameter==g.portParameters[k])a.position=g.ports[k];}
    return cuts;
}
inline std::vector<GraphArc> graphArcs(const QuadraticGraphExact&g,const PlanarCellCut::Rectangle&b,unsigned ownership=15,std::vector<GraphEvent>* retainedEvents=nullptr) {
    using namespace QuadraticExact;auto cuts=graphBoundaryEvents(g,b,ownership);
    std::vector<GraphArc> arcs;for(std::size_t i=1;i<cuts.size();++i){const auto a=cuts[i-1],c=cuts[i];if(a.parameter==c.parameter)continue;
        const auto p=graphPosition(g,(a.parameter+c.parameter)/2);
        if(p.x<scalar(b.lower.x)||p.x>scalar(b.upper.x)||p.y<scalar(b.lower.y)||p.y>scalar(b.upper.y))continue;
        if(!arcs.empty()&&arcs.back().last.parameter==a.parameter)arcs.back().last=c;else arcs.push_back({a,c});}
    if(retainedEvents)*retainedEvents=std::move(cuts);return arcs;
}
}
