#include "Core/Geometry/QuadraticChordEvents.h"
#include "Core/Geometry/Detail/PolynomialSegmentRoots.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <iomanip>

namespace LAMP::Geometry::QuadraticChordEvents
{
namespace {
namespace E=Detail::ExactExpansion;
using Polynomial=Detail::PolynomialSegmentRoots::ExpandedPolynomial;
void require(bool ok){if(!ok) throw std::invalid_argument("[QuadraticChordEvents] invalid finite sweep or graph");}
void validate(const Sweep& s){for(auto p:{s.begin0,s.end0,s.begin1,s.end1}) require(std::isfinite(p.x)&&std::isfinite(p.y));}
E::Expansion scalar(long double value){
    require(std::isfinite(value));E::Expansion result;
    for(int k=0;value!=0 && k<4;++k) {const double v=static_cast<double>(value);
        require(std::isfinite(v)&&v!=0);result=E::add(result,E::fromDouble(v));value-=static_cast<long double>(v);}
    require(value==0);return result;
}
Polynomial line(long double first,long double last){Polynomial r{};r[0]=scalar(first);r[1]=E::subtract(scalar(last),r[0]);return r;}
Polynomial add(const Polynomial& a,const Polynomial& b){Polynomial r{};
    for(int k=0;k<7;++k) r[k]=E::add(a[k],b[k]);return r;}
Polynomial subtract(const Polynomial& a,const Polynomial& b){Polynomial r{};
    for(int k=0;k<7;++k) r[k]=E::subtract(a[k],b[k]);return r;}
Polynomial multiply(const Polynomial& a,const Polynomial& b){Polynomial r{};
    for(int i=0;i<7;++i) for(int j=0;i+j<7;++j) r[i+j]=E::add(r[i+j],E::multiply(a[i],b[j]));return r;}
Polynomial scaled(const Polynomial& p,long double factor){Polynomial r{};const auto v=scalar(factor);
    for(int k=0;k<7;++k) r[k]=E::multiply(p[k],v);return r;}
UnitPolynomialRoots isolateEvents(const Polynomial& p)
{
    namespace R=Detail::PolynomialSegmentRoots;
    E::validateEnvironment();
    int degree=6;while(degree>0 && E::sign(p[degree])==0) --degree;
    const auto add=E::add,sub=E::subtract,mul=E::multiply;
    const auto scaled=[&](const E::Expansion& x,double v){return mul(x,E::fromDouble(v));};
    std::vector<R::ExpandedPolynomial> factors;
    const auto linear=[&](E::Expansion constant,E::Expansion slope){R::ExpandedPolynomial q{};
        q[0]=std::move(constant);q[1]=std::move(slope);factors.push_back(std::move(q));};
    // Discriminants multiply as many as four coefficients. A rotated frame
    // may give a perfectly ordinary polynomial long exact tails, especially
    // close to an axis. Scale ALL coefficients by one exact positive power
    // of two before forming these products; roots and exact multiplicities
    // are unchanged. Never truncate a tail or declare a small term zero.
    auto predicate=p;bool possibleMultiple=degree==2 || degree==3;
    if(degree==2 || degree==3) {
        int least=1024,highest=-1074;
        for(int k=0;k<=degree;++k) for(double value:p[k].components()) {
            least=std::min(least,E::Detail::leastBitExponent(value));
            highest=std::max(highest,std::ilogb(std::abs(value))+2);
        }
        // A nonzero outward enclosure certifies a simple polynomial without
        // multiplying every exact tail four times. This is especially useful
        // for almost straight rotated donors. Keep the ORIGINAL polynomial
        // for root isolation; the filter never rounds away coefficients or
        // calls an uncertain discriminant zero.
        if(least-highest>=-1074) {
            namespace I=R::Detail;
            std::array<I::Interval,4> c{};
            for(int k=0;k<=degree;++k) {
                const auto bound=E::enclosure(E::scalePow2(p[k],-highest));
                c[k]={bound.lo,bound.hi};
            }
            const auto plus=I::add,product=I::multiply;
            const auto times=[&](I::Interval x,double n){return product(x,{n,n});};
            const auto minus=[&](I::Interval x,I::Interval y){return plus(x,{-y.hi,-y.lo});};
            auto bound=minus(product(c[1],c[1]),times(product(c[2],c[0]),4));
            if(degree==3) {
                const auto a=c[3],b=c[2],v=c[1],d=c[0];
                bound=minus(plus(times(product(product(a,b),product(v,d)),18),product(product(b,b),product(v,v))),
                    plus(plus(times(product(product(product(b,b),b),d),4),times(product(a,product(product(v,v),v)),4)),
                        times(product(product(a,a),product(d,d)),27)));
            }
            possibleMultiple=I::sign(bound)==0;
        }
        const int lower=(degree==2?-500:-250)-least,upper=(degree==2?480:240)-highest;
        // Deflation is optional. If its global discriminant is unrepresentable,
        // the original strict isolator can still certify a regular polynomial
        // on [0,1] (e.g. a multiple root outside the sweep). An unresolved
        // interior multiple root still throws there; no event is discarded.
        if(lower>upper) possibleMultiple=false;
        if(possibleMultiple) {
            const int shift=std::clamp(0,lower,upper);
            if(shift!=0) for(int k=0;k<=degree;++k) predicate[k]=scaled(p[k],std::ldexp(1.,shift));
        }
    }
    // A temporal multiple event is admissible here: it is a subdivision seed,
    // not a regular spatial interface chart. Deflate ONLY algebraically exact
    // zero discriminants; retain the strict existing engine for near-roots.
    if(possibleMultiple && degree==2 && E::sign(sub(mul(predicate[1],predicate[1]),scaled(mul(predicate[2],predicate[0]),4)))==0)
        linear(predicate[1],scaled(predicate[2],2));
    if(possibleMultiple && degree==3) {
        const auto a=predicate[3],b=predicate[2],c=predicate[1],d=predicate[0];
        const auto discriminant=sub(add(scaled(mul(mul(a,b),mul(c,d)),18),mul(mul(b,b),mul(c,c))),
            add(add(scaled(mul(mul(mul(b,b),b),d),4),scaled(mul(a,mul(mul(c,c),c)),4)),scaled(mul(mul(a,a),mul(d,d)),27)));
        if(E::sign(discriminant)==0) {
            const auto delta=sub(mul(b,b),scaled(mul(a,c),3));
            if(E::sign(delta)==0) linear(b,scaled(a,3));
            else {
                linear(sub(mul(b,c),scaled(mul(a,d),9)),scaled(delta,2));
                linear(sub(add(mul(mul(b,b),b),scaled(mul(mul(a,a),d),9)),scaled(mul(mul(a,b),c),4)),mul(a,delta));
            }
        }
    }
    if(factors.empty()) factors.push_back(p);
    UnitPolynomialRoots result;
    for(const auto& q:factors) {
        R::Result roots;
        try {roots=R::isolateEndpointRelative(q);}
        catch(const std::runtime_error& error) {
            std::ostringstream message;message<<"[QuadraticChordEvents] polynomial:"<<std::setprecision(17);
            for(const auto& coefficient:q) {message<<" [";for(double value:coefficient.components()) message<<value<<',';message<<']';}
            message<<"; "<<error.what();throw std::runtime_error(message.str());
        }
        result.identicallyZero|=roots.identicallyZero;
        for(std::size_t k=0;k<roots.roots.size();++k)
            result.roots.push_back({roots.roots[k],roots.enclosures[k].lo,roots.enclosures[k].hi});
    }
    std::sort(result.roots.begin(),result.roots.end(),[](const auto& a,const auto& b){return a.point<b.point;});
    for(std::size_t k=1;k<result.roots.size();++k) if(!(result.roots[k-1].upper<result.roots[k].lower)) {
        std::ostringstream message;message<<std::setprecision(17)<<"[QuadraticChordEvents] overlapping deflated event roots: "
            <<result.roots[k-1].lower<<','<<result.roots[k-1].upper<<" and "<<result.roots[k].lower<<','<<result.roots[k].upper;
        throw std::invalid_argument(message.str());
    }
    return result;
}
Point interpolate(Point a,Point b,double t){return {std::fma(static_cast<long double>(t),b.x-a.x,a.x),
    std::fma(static_cast<long double>(t),b.y-a.y,a.y)};}
struct Frame {long double n0,v0,dn,dv;};
Frame frame(const PlanarCellCut::Rectangle& box,const QuadraticCellCut::Reconstruction& curve,Point a,Point b){
    const auto s=curve.shape;
    require(std::isfinite(box.lower.x)&&std::isfinite(box.lower.y)&&std::isfinite(box.upper.x)&&std::isfinite(box.upper.y)
        && box.lower.x<box.upper.x && box.lower.y<box.upper.y);
    require((s.axis==0 || s.axis==1) && (s.sign==1 || s.sign==-1) && std::isfinite(s.slope)
        && std::isfinite(s.second) && std::isfinite(curve.offset));
    const long double cx=(static_cast<long double>(box.lower.x)+box.upper.x)/2;
    const long double cy=(static_cast<long double>(box.lower.y)+box.upper.y)/2;
    if(s.frameNormal) {
        const auto basis=QuadraticCellCut::basis(s);
        const auto dot=[](Point p,Point q){return p.x*q.x+p.y*q.y;};
        const Point first{a.x-cx,a.y-cy},delta{b.x-a.x,b.y-a.y};
        return {dot(first,basis.normal),dot(first,basis.tangent),dot(delta,basis.normal),dot(delta,basis.tangent)};
    }
    return s.axis?Frame{a.y-cy,a.x-cx,b.y-a.y,b.x-a.x}:Frame{a.x-cx,a.y-cy,b.x-a.x,b.y-a.y};
}
}
UnitPolynomialRoots intersections(const PlanarCellCut::Rectangle& box,
    const QuadraticCellCut::Reconstruction& curve,Point a,Point b)
{
    validate({a,b,a,b});frame(box,curve,a,b);E::validateEnvironment();
    const auto cx=scaled(add(line(box.lower.x,box.lower.x),line(box.upper.x,box.upper.x)),.5);
    const auto cy=scaled(add(line(box.lower.y,box.lower.y),line(box.upper.y,box.upper.y)),.5);
    const auto x=subtract(line(a.x,b.x),cx),y=subtract(line(a.y,b.y),cy);
    auto n=curve.shape.axis?y:x,v=curve.shape.axis?x:y;
    if(curve.shape.frameNormal) {const auto b=QuadraticCellCut::basis(curve.shape);
        n=add(scaled(x,b.normal.x),scaled(y,b.normal.y));v=add(scaled(x,b.tangent.x),scaled(y,b.tangent.y));}
    return isolateEvents(subtract(subtract(subtract(n,line(curve.offset,curve.offset)),scaled(v,curve.shape.slope)),
        scaled(multiply(v,v),curve.shape.second/2)));
}
long double vertexValue(Point a,Point b,Point p){return (b.x-a.x)*(p.y-a.y)-(b.y-a.y)*(p.x-a.x);}
UnitPolynomialRoots vertex(const Sweep& s,Point p)
{
    validate(s);require(std::isfinite(p.x)&&std::isfinite(p.y));E::validateEnvironment();
    const auto ax=line(s.begin0.x,s.begin1.x),ay=line(s.begin0.y,s.begin1.y);
    const auto dx=subtract(line(s.end0.x,s.end1.x),ax),dy=subtract(line(s.end0.y,s.end1.y),ay);
    const auto px=subtract(line(p.x,p.x),ax),py=subtract(line(p.y,p.y),ay);
    auto result=isolateEvents(subtract(multiply(dx,py),multiply(dy,px)));
    result.roots.erase(std::remove_if(result.roots.begin(),result.roots.end(),[&](const auto& r){
        const auto a=interpolate(s.begin0,s.begin1,r.point),b=interpolate(s.end0,s.end1,r.point);
        const auto dx=b.x-a.x,dy=b.y-a.y,length=dx*dx+dy*dy;
        if(length==0) return true;const auto t=((p.x-a.x)*dx+(p.y-a.y)*dy)/length;return t<0 || t>1;
    }),result.roots.end());return result;
}
long double tangencyValue(const PlanarCellCut::Rectangle& box,const QuadraticCellCut::Reconstruction& curve,Point a,Point b)
{
    const auto f=frame(box,curve,a,b);const auto d=f.dn-curve.shape.slope*f.dv;
    // Discriminant after cancelling its equal fourth-degree terms. This
    // factored form also retains the small offset of a near-tangent chord.
    return d*d+2*curve.shape.second*f.dv*(f.dv*(f.n0-curve.offset)-f.v0*f.dn);
}
UnitPolynomialRoots tangencies(const PlanarCellCut::Rectangle& box,
    const QuadraticCellCut::Reconstruction& curve,const Sweep& s)
{
    validate(s);E::validateEnvironment();const auto f=frame(box,curve,s.begin0,s.end0),g=frame(box,curve,s.begin1,s.end1);
    if(curve.shape.second==0 || (f.dv==0 && g.dv==0)) return {};
    const auto ax=line(s.begin0.x,s.begin1.x),ay=line(s.begin0.y,s.begin1.y);
    const auto dx=subtract(line(s.end0.x,s.end1.x),ax),dy=subtract(line(s.end0.y,s.end1.y),ay);
    const auto cx=scaled(add(line(box.lower.x,box.lower.x),line(box.upper.x,box.upper.x)),.5);
    const auto cy=scaled(add(line(box.lower.y,box.lower.y),line(box.upper.y,box.upper.y)),.5);
    auto dn=curve.shape.axis?dy:dx,dv=curve.shape.axis?dx:dy;
    auto n0=subtract(subtract(curve.shape.axis?ay:ax,curve.shape.axis?cy:cx),line(curve.offset,curve.offset));
    auto v0=subtract(curve.shape.axis?ax:ay,curve.shape.axis?cx:cy);
    if(curve.shape.frameNormal) {const auto b=QuadraticCellCut::basis(curve.shape);const auto x=subtract(ax,cx),y=subtract(ay,cy);
        dn=add(scaled(dx,b.normal.x),scaled(dy,b.normal.y));dv=add(scaled(dx,b.tangent.x),scaled(dy,b.tangent.y));
        n0=subtract(add(scaled(x,b.normal.x),scaled(y,b.normal.y)),line(curve.offset,curve.offset));
        v0=add(scaled(x,b.tangent.x),scaled(y,b.tangent.y));}
    const auto d=subtract(dn,scaled(dv,curve.shape.slope));
    auto result=isolateEvents(add(multiply(d,d),scaled(multiply(dv,subtract(multiply(dv,n0),multiply(v0,dn))),2*curve.shape.second)));
    result.roots.erase(std::remove_if(result.roots.begin(),result.roots.end(),[&](const auto& r){
        const auto a=interpolate(s.begin0,s.begin1,r.point),b=interpolate(s.end0,s.end1,r.point);
        const auto q=frame(box,curve,a,b);const auto A=-curve.shape.second*q.dv*q.dv/2;
        const auto B=q.dn-curve.shape.slope*q.dv-curve.shape.second*q.v0*q.dv;
        if(A==0) return true;const auto t=-B/(2*A);return !(t>0 && t<1);
    }),result.roots.end());return result;
}
}
