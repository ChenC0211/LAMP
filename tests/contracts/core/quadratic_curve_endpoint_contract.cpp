#include "Core/Geometry/QuadraticCurveMeasure.h" // Public header stands alone.
#include "Core/Geometry/QuadraticCurveDifference.h"
#include "Core/Mesh/CartesianMeasure.h"
#include <algorithm>
#include <cmath>
#include <cfenv>
#include <string>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace
{
namespace C=LAMP::Geometry::QuadraticCurveMeasure;
using Point=C::Point;
unsigned assertions=0,guards=0;
void check(bool ok,const char* why){++assertions;if(!ok)throw std::runtime_error(why);}
void near(long double a,long double b,long double scale,const char* why)
{check(std::isfinite(a)&&std::isfinite(b)&&scale>0,"invalid comparison");check(std::abs(a-b)<=1e-12L*scale,why);}
template<class F>void rejects(F&& f,const char* why)
{bool caught=false;try{f();}catch(const std::exception&){caught=true;}check(caught,why);++guards;}
struct ReferencePoint {long double x=0,y=0;};
ReferencePoint independent(Point first,Point middle,Point last,long double s,int order)
{
    // Direct Lagrange evaluation of the three original samples, independently
    // of the production endpoint/bubble representation and its roots.
    long double a,b,c;
    if(order==0){a=(1-s)*(1-2*s);b=4*s*(1-s);c=s*(2*s-1);}
    else if(order==1){a=4*s-3;b=4-8*s;c=4*s-1;}
    else{a=4;b=-8;c=4;}
    return {a*first.x+b*middle.x+c*last.x,a*first.y+b*middle.y+c*last.y};
}
long double referenceArea(Point first,Point middle,Point last,bool rz)
{
    constexpr long double x[5]={-.90617984593866399280L,-.53846931010568309104L,0,.53846931010568309104L,.90617984593866399280L};
    constexpr long double w[5]={.23692688505618908751L,.47862867049936646804L,.56888888888888888889L,.47862867049936646804L,.23692688505618908751L};
    long double result=0;
    for(int panel=0;panel<32;++panel){const long double mid=(panel+.5L)/32,half=1.L/64;
        for(int q=0;q<5;++q){const auto s=mid+half*x[q];const auto p=independent(first,middle,last,s,0),d=independent(first,middle,last,s,1);
            check(!rz||p.x>=0,"independent positive interior RZ samples");
            result+=half*w[q]*std::hypot(d.x,d.y)*(rz?CartesianMeasure::twoPi*p.x:1);}}
    return result;
}
void endpointCases()
{
    const Point first{.1,.2},middle{.15,.4},last{0,.6};
    const auto curve=C::Curve::fromSamples(first,middle,last);
    check(curve.usesEndpointBasis(),"explicit three-sample endpoint authority");
    check(C::position(curve,0).x==first.x&&C::position(curve,0).y==first.y&&
        C::position(curve,1).x==last.x&&C::position(curve,1).y==last.y,"both original endpoints exactly retained");
    // Historical imageCurve arithmetic fails on these genuinely nonnegative
    // samples. This is the motivating last-axis regression, not an axis clamp.
    C::Curve old;old.coefficients[0]=first;
    old.coefficients[2]={2*((last.x-first.x)-2*(middle.x-first.x)),2*((last.y-first.y)-2*(middle.y-first.y))};
    old.coefficients[1]={last.x-first.x-old.coefficients[2].x,last.y-first.y-old.coefficients[2].y};
    check(C::position(old,1).x<0,"original double power assembly exposes the negative last-axis residual");
    rejects([&]{C::integrate(old,0,1,true);},"legacy genuinely negative power curve remains strict");
    const auto reverse=C::Curve::fromSamples(last,middle,first);
    check(C::position(reverse,0).x==0&&C::position(curve,1).x==0,"axis-first and axis-last are both exact");
    for(int k=0;k<=32;++k){const double s=static_cast<double>(k)/32;const auto p=C::position(curve,s),r=C::position(reverse,1-s);
        const auto expected=independent(first,middle,last,s,0);near(p.x,expected.x,1,"same physical x polynomial");near(p.y,expected.y,1,"same physical y polynomial");
        near(p.x,r.x,1,"parameter reversal x");near(p.y,r.y,1,"parameter reversal y");
        for(int order=1;order<=2;++order){const auto d=C::derivative(curve,s,order);const auto e=independent(first,middle,last,s,order);
            near(d.x,e.x,1,"same-basis x derivative");near(d.y,e.y,1,"same-basis y derivative");}
        check(p.x>=0&&r.x>=0,"both stable orientations retain nonnegative physical radius");}
    const auto b=C::bounds(curve),br=C::bounds(reverse);check(b.lower.x==0&&br.lower.x==0,"bounds use authoritative axis endpoints");
    near(b.upper.x,.15625L,1,"interior radial stationary maximum retained");
    for(bool rz:{false,true}){const auto full=C::integrate(curve,0,1,rz),reversed=C::integrate(reverse,0,1,rz);
        near(full.physicalMeasure,referenceArea(first,middle,last,rz),1,"positive independent curve integral");
        near(full.physicalMeasure,reversed.physicalMeasure,1,"same integral in reverse parameter direction");}
    const auto endRoots=C::coordinateRoots(curve,0,0);check(endRoots.size()==1&&endRoots[0]==1,"exact last-axis clipping root retained");
    const auto firstRoots=C::coordinateRoots(reverse,0,0);check(firstRoots.size()==1&&firstRoots[0]==0,"exact first-axis clipping root retained");
    const auto crossings=C::coordinateRoots(curve,0,.1);check(crossings.size()==2&&crossings.front()==0,"endpoint and interior coordinate roots both retained");
    near(crossings.back(),.75L,1,"independent interior coordinate root");
    const auto whole=C::intervals(curve,{{0,.2},{.2,.6}});check(whole.size()==1&&whole.front()[0]==0&&whole.front()[1]==1,"closed radial clipping loses no terminal interval");
    const auto lower=C::intervals(curve,{{0,.2},{.2,.4}}),upper=C::intervals(curve,{{0,.4},{.2,.6}});
    check(lower.size()==1&&upper.size()==1&&lower[0][1]==upper[0][0],"adjacent restrictions share the same stable root");
    near(lower[0][1],.5L,1,"independent half-interval clipping");
    const auto l=C::integrate(curve,lower[0][0],lower[0][1],true),u=C::integrate(curve,upper[0][0],upper[0][1],true);
    near(l.physicalMeasure+u.physicalMeasure,C::integrate(curve,0,1,true).physicalMeasure,1,"adjacent true RZ restrictions are additive");
}
void pairedDifferenceCases()
{
    const auto straight=C::Curve::fromSamples({0,0},{.5,0},{1,0});
    const double a=std::ldexp(1.,-30);
    const auto tiny=C::Curve::fromSamples({0,0},{.5,a/4},{1,0});
    const long double a2=static_cast<long double>(a)*a;
    // Integral_0^1 [sqrt(1+a^2*(1-2s)^2)-1] ds has the alternating
    // binomial series a^2/6-a^4/40+a^6/112-5*a^8/1152+... .
    // For 0<a<=1, truncation after a^4 has 0<R<=a^6/112.
    // This reference never subtracts two nearly equal square roots/lengths.
    const long double reference=a2/6-a2*a2/40,remainderBound=a2*a2*a2/112;
    for(bool rz:{false,true}){
        const long double weight=rz?CartesianMeasure::twoPi/2:1;
        // r=s and the speed increment is invariant under s -> 1-s,
        // hence Integral 2*pi*s*deltaSpeed ds = pi*deltaLength.
        const auto forward=C::physicalMeasureDifference(straight,tiny,rz);
        const auto backward=C::physicalMeasureDifference(tiny,straight,rz);
        const auto identity=C::physicalMeasureDifference(tiny,tiny,rz);
        const auto expected=weight*reference;
        check(std::isfinite(forward.delta)&&forward.delta>0,"paired tiny physical difference stays positive");
        check(std::isfinite(forward.quadratureErrorEstimate)&&forward.quadratureErrorEstimate>=0,
            "paired tiny quadrature error is a finite nonnegative estimate");
        check(std::abs(forward.delta-expected)<=1e-12L*expected+weight*remainderBound,
            "paired tiny XY/RZ difference agrees with independent series and next-term bound");
        check(identity.delta==0&&identity.quadratureErrorEstimate==0,"paired identity has exact zero difference and error");
        near(forward.delta,-backward.delta,expected,"paired physical difference is antisymmetric at its own scale");
    }
    const auto normal=C::Curve::fromSamples({0,0},{.5,.25},{1,0});
    const long double normalReference=(std::sqrt(2.L)+std::asinh(1.L))/2-1;
    for(bool rz:{false,true}){
        const auto expected=(rz?CartesianMeasure::twoPi/2:1)*normalReference;
        near(C::physicalMeasureDifference(straight,normal,rz).delta,expected,expected,
            "normal-scale paired physical difference retains original relative accuracy");
    }
    rejects([&]{C::physicalMeasureDifference(straight,normal,false,0);},"paired nonpositive tolerance is rejected");
    rejects([&]{C::physicalMeasureDifference(straight,normal,false,std::numeric_limits<long double>::quiet_NaN());},
        "paired nonfinite tolerance is rejected");
    rejects([&]{C::physicalMeasureDifference(straight,normal,false,1e-13L,33);},"paired depth beyond original upper bound is rejected");
    rejects([&]{C::physicalMeasureDifference(straight,normal,false,1e-30L,0);},"paired unresolved normal-scale quadrature obeys original depth budget");
}
void exactGraphCases()
{
    namespace Q=LAMP::Geometry::QuadraticCellCut;
    namespace S=LAMP::Geometry::QuadraticSurfaceMeasure;
    // Actual repair8 cell29 shared ports and eta. Its native centered offset
    // cannot retain the small correction at the original own-phase/root gates.
    const Q::Rectangle box{{.625,.375},{.75,.5}};
    const Q::Point first{.656249999702048115,.5},last{.71875,.5};
    constexpr long double second=2.80626114169708585e-5L,h=.125L,eta=second*h/2;
    const auto graph=Q::fromFixedPorts(box,{1,1,0,0,{}},first,last,eta);
    check(Q::implicitValue(box,graph,first)==0&&Q::implicitValue(box,graph,last)==0,"original fixed ports are exact roots of the same graph");
    const auto surface=S::evaluate(box,graph,box,false);
    check(surface.patches.size()==1,"actual thin cap is one connected arc");
    const auto& patch=surface.patches[0];
    const auto curve=C::Curve::fromGraph(Q::materialCurveData(box,graph,patch.lower,patch.upper));
    check(curve.usesExactBasis(),"Core exact material bridge is explicitly retained");
    const auto a=C::position(curve,0),b=C::position(curve,1);
    check(a.x==static_cast<double>(first.x)&&a.y==static_cast<double>(first.y)&&
        b.x==static_cast<double>(last.x)&&b.y==static_cast<double>(last.y),"both original binary fixed ports survive bit exactly, no re-fit or snap");
    const long double span=last.x-first.x,k=eta/h;
    // Reuse this actual exact-payload graph at a non-cancellation scale.
    // The original graph and all its existing gates below stay unchanged.
    const auto bent=C::displacedSamples(curve,{0,0},{0,static_cast<double>(span/8)},{0,0});
    check(bent.usesExactBasis(),"paired reference retains the same exact-payload curve authority");
    for(bool rz:{false,true}){
        const auto before=C::integrate(curve,0,1,rz,1e-15L),after=C::integrate(bent,0,1,rz,1e-15L);
        const auto paired=C::physicalMeasureDifference(curve,bent,rz);
        const auto reference=after.physicalMeasure-before.physicalMeasure;
        check(reference>0,"independent exact-payload reference has a non-cancellation physical difference");
        check(std::abs(paired.delta-reference)<=1e-12L*reference+before.quadratureErrorEstimate+after.quadratureErrorEstimate,
            "paired exact-payload physical difference agrees with separate true physical integrals");
    }
    for(double s:{0.,.25,.5,.75,1.}){
        const auto point=C::position(curve,s),d=C::derivative(curve,s),dd=C::derivative(curve,s,2);
        const auto original=S::sample(box,graph,patch.lower+(patch.upper-patch.lower)*s,false);
        near(point.x,original.position.x,span,"same original graph radial coordinate");
        near(point.y,original.position.y,1,"same original graph axial coordinate");
        near(d.x,span,span,"exact bridge original radial derivative");
        near(d.y,k*span*span*(2*s-1),k*span*span,"exact bridge derivative of the same thin bubble");
        near(dd.x,0,span,"exact bridge zero radial second derivative");
        near(dd.y,2*k*span*span,k*span*span,"exact bridge retains the thin quadratic coefficient");
    }
    const auto roof=C::coordinateRoots(curve,1,.5);
    check(roof.size()==2&&roof.front()==0&&roof.back()==1,"both coordinate-boundary roots keep exact original material parameters");
    const auto bounds=C::bounds(curve);
    check(bounds.lower.x==a.x&&bounds.upper.x==b.x&&bounds.upper.y==.5,"exact cap bounds retain coordinate boundary and port identity");
    check(bounds.lower.y<.5,"a genuinely nonzero thin bubble, not a straight chord");
    const auto middle=C::position(curve,.5).x;const auto roots=C::coordinateRoots(curve,0,middle);
    check(roots.size()==1,"single exact radial restriction root");
    near(roots[0],(static_cast<long double>(middle)-first.x)/span,1,"independent linear radial root");
    const auto left=C::intervals(curve,{{box.lower.x,box.lower.y},{middle,box.upper.y}});
    const auto right=C::intervals(curve,{{middle,box.lower.y},{box.upper.x,box.upper.y}});
    check(left.size()==1&&right.size()==1&&left[0][0]==0&&right[0][1]==1&&left[0][1]==right[0][0],
        "adjacent restrictions reuse exactly the same material root");
    for(const bool rz:{false,true}){
        auto gas=graph;gas.shape.sign=-gas.shape.sign;
        const auto minor=Q::graphMeasure(box,gas,rz);
        const auto expected=k*span*span*span/6*(rz?CartesianMeasure::twoPi*(first.x/2+last.x/2):1.L);
        check(minor.volume>0,"same encoded quadratic retains its positive tiny gas phase");
        near(minor.volume,expected,expected,"independent positive XY/RZ cap measure at the original own-phase gate");
        const auto measured=C::integrate(curve,0,1,rz);const auto original=S::evaluate(box,graph,box,rz);
        near(measured.physicalMeasure,original.integral.physicalMeasure,measured.physicalMeasure,"exact material curve and same graph physical area");
        const auto lower=C::integrate(curve,left[0][0],left[0][1],rz),upper=C::integrate(curve,right[0][0],right[0][1],rz);
        near(lower.physicalMeasure+upper.physicalMeasure,measured.physicalMeasure,measured.physicalMeasure,
            "same exact parameter restrictions add back in actual XY/RZ area");
    }
}

void exactDerivativeCacheCases()
{
    namespace Q=LAMP::Geometry::QuadraticCellCut;namespace S=LAMP::Geometry::QuadraticSurfaceMeasure;
    const Q::Rectangle box{{0,0},{1,1}};Q::Shape shape;shape.axis=1;
    const auto make=[&](Q::Point first,Q::Point last,long double eta){
        const auto graph=Q::fromFixedPorts(box,shape,first,last,eta);
        const auto surface=S::evaluate(box,graph,box,false);check(surface.patches.size()==1,"cache fixture has one actual source arc");
        const auto& patch=surface.patches[0];return C::Curve::fromGraph(Q::materialCurveData(box,graph,patch.lower,patch.upper));
    };
    // x=s, y=1/4+s/2+s^2/4 gives exact dyadic reference derivatives.
    const auto curve=make({0,.25},{1,1},.25L);
    const auto identical=[](double a,double b){return a==b&&std::signbit(a)==std::signbit(b);};
    for(const auto s:{0.,.125,.25,.5,.75,1.}){
        const auto dd=C::derivative(curve,s,2),d=C::derivative(curve,s);
        check(identical(d.x,1)&&identical(d.y,.5+.5*s),"cached first derivative retains independent dyadic polynomial bits");
        check(identical(dd.x,0)&&identical(dd.y,.5),"cached second derivative retains its original exact constant");
    }
    const auto tie=std::nextafter(.5,1.);const auto beyond=std::nextafter(tie,1.);
    check(C::derivative(curve,tie).y==.75,"first derivative exact native tie keeps nearest-even result");
    check(C::derivative(curve,beyond).y==std::nextafter(.75,1.),"first derivative beyond tie keeps next native result");
    const auto cancelling=make({0,.25},{1,.25},-.5L); // y=1/4+s/2-s^2/2.
    for(const auto s:{std::nextafter(.5,0.),.5,std::nextafter(.5,1.)})
        check(identical(C::derivative(cancelling,s).y,.5-s),"strong derivative cancellation keeps original rational result");
    auto stale=curve;stale.coefficients[0].x=std::nextafter(stale.coefficients[0].x,1.);
    for(int order:{1,2}){bool guarded=false;try{(void)C::derivative(stale,.25,order);}
        catch(const std::exception& e){guarded=std::string(e.what()).find("stale coefficient view")!=std::string::npos;}
        check(guarded,"stale diagnostic view must reject before either derivative cache");}
    const auto currentRound=std::fegetround();
    if(currentRound==FE_TONEAREST){
        // The exact tie is representable in extended long double; binary64
        // long double instead rounds it to .75 BEFORE Point's native cast.
        volatile long double referenceTie=std::numeric_limits<long double>::digits>
            std::numeric_limits<double>::digits?.75L+std::ldexp(1.L,-std::numeric_limits<double>::digits-1):.75L;
        for(int mode:{FE_DOWNWARD,FE_UPWARD})if(std::fesetround(mode)==0){
            struct Restore{int value;~Restore(){std::fesetround(value);}}restore{currentRound};
            volatile double expected=static_cast<double>(referenceTie);
            check(identical(C::derivative(curve,tie).y,expected),"long-lived curve must honor original query-environment conversion");
            check(C::derivative(curve,.25,2).y==.5,"long-lived second derivative preserves nonnearest fallback");
        }
    }
    const auto dome=make({0,0},{1,0},-1.L); // y=s-s^2.
    const auto scale=.75*std::numeric_limits<double>::max();
    const auto large=C::affineMapped(dome,{0,0},{0,0},{{1,0,0,scale}});
    check(C::position(large,0).y==0&&std::isfinite(C::derivative(large,.25).y),"order-two cache failure must not reject curve creation or valid first derivative");
    rejects([&]{(void)C::derivative(large,.25,2);},"unrepresentable second derivative retains its original query-time failure");
}

void strictAndLegacy()
{
    const auto cross=C::Curve::fromSamples({.1,.2},{-.05,.4},{.1,.6});check(C::bounds(cross).lower.x<0,"actual interior negative stationary radius retained");
    rejects([&]{C::integrate(cross,0,1,true);},"genuine interior cross-axis motion cannot use endpoint exception");
    const auto negative=C::Curve::fromSamples({-.01,.2},{.1,.4},{.2,.6});
    rejects([&]{C::integrate(negative,0,1,true);},"genuine negative endpoint cannot use endpoint exception");
    const auto axis=C::Curve::fromSamples({0,.2},{0,.4},{0,.6});
    const auto all=C::intervals(axis,{{0,0},{1,1}});
    check(all.size()==1&&all[0][0]==0&&all[0][1]==1,"true axis-coincident interval retained with ownership");
    check(C::intervals(axis,{{0,0},{1,1}},LAMP::Geometry::QuadraticSurfaceMeasure::All&
        ~LAMP::Geometry::QuadraticSurfaceMeasure::Left).empty(),"excluded original side ownership remains exact");
    check(C::integrate(axis,0,1,true).physicalMeasure==0,"true radial-axis curve has exact zero physical area");
    C::Curve legacy;legacy.coefficients={{{.25,.25},{-.25,.5},{0,0}}};
    check(!legacy.usesEndpointBasis()&&C::position(legacy,1).x==0,"direct coefficients retain legacy construction");
    near(C::integrate(legacy,0,1,false).length,std::hypot(.25L,.5L),1,"legacy line measure retained");
    const auto zero=C::derivative(axis,.3,3);check(zero.x==0&&zero.y==0,"quadratic higher derivative is exact zero");
    rejects([&]{C::derivative(axis,.5,0);},"nonpositive derivative order");
    rejects([&]{C::coordinateRoots(axis,2,0);},"invalid coordinate axis");
    rejects([&]{C::Curve::fromSamples({0,0},{std::numeric_limits<double>::infinity(),0},{1,1});},"nonfinite original middle sample");
}
}
int main()
{
    try{endpointCases();strictAndLegacy();exactGraphCases();pairedDifferenceCases();exactDerivativeCacheCases();std::cout<<"quadratic curve endpoint contract PASS assertions="<<assertions<<" guards="<<guards<<'\n';return 0;}
    catch(const std::exception& error){std::cerr<<"quadratic curve endpoint contract FAIL "<<error.what()<<'\n';return 1;}
}