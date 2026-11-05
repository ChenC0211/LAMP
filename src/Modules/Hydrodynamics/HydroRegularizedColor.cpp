#include "Modules/Hydrodynamics/HydroRegularizedColor.h"
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Geometry/Detail/ExactExpansion.h"
#include "Core/Geometry/PolynomialVelocitySlab.h"
#include "Core/Mesh/CartesianMeasure.h"
#include <algorithm>
#include <cmath>
#include <exception>
#include <queue>
#include <stdexcept>
#include <string>
#include <utility>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace LAMP::Hydrodynamics
{
namespace
{
constexpr double pi=3.1415926535897932384626433832795;
void require(bool ok,const std::string& why)
{if(!ok)throw std::invalid_argument("[HydroRegularizedColor] "+why);}
void width(double e){require(std::isfinite(e)&&e>0,"epsilon must be a finite positive physical length");}
double lowerH(double t)
{
    // H(-1+a)=a/2*(1-sinc(pi*a)). This positive endpoint expansion avoids
    // cancellation and does not clip a negative approximation to zero.
    const double a=1+t,z=pi*a;
    if(a<.25){double term=z*z/6,sum=term;
        for(int k=2;k<=12;++k){term*=-z*z/((2.*k)*(2*k+1));sum+=term;}
        return .5*a*sum;}
    return .5*(a-std::sin(z)/pi);
}
long double volume(Vector2D a,Vector2D b,bool rz)
{return (static_cast<long double>(b.x)-a.x)*(static_cast<long double>(b.y)-a.y)*(rz?CartesianMeasure::twoPi*(static_cast<long double>(a.x)+b.x)/2:1.L);}
struct BoxBounds
{
    double phiLower=0,phiUpper=0,gradientUpper=0,gradientLower=0;
    double cLower=0,cUpper=1;int pure=0;bool allBand=false;
};
using Range=LAMP::Geometry::PolynomialVelocityRange;
double up(double x){return std::nextafter(x,std::numeric_limits<double>::infinity());}
double down(double x){return std::nextafter(x,-std::numeric_limits<double>::infinity());}
Range add(Range a,Range b)
{return {std::nextafter(a.lower+b.lower,-std::numeric_limits<double>::infinity()),
    std::nextafter(a.upper+b.upper,std::numeric_limits<double>::infinity())};}
Range multiply(Range a,double v)
{const double x=a.lower*v,y=a.upper*v;return {std::nextafter(std::min(x,y),-std::numeric_limits<double>::infinity()),
    std::nextafter(std::max(x,y),std::numeric_limits<double>::infinity())};}
Range positiveProduct(Range a,Range b)
{return {std::max(0.,down(a.lower*b.lower)),up(a.upper*b.upper)};}
Range positiveDivide(Range a,double d)
{return {std::max(0.,down(a.lower/d)),up(a.upper/d)};}
double hBound(double s,double e,bool upper)
{
    if(s<=-e)return 0;if(s>=e)return 1;
    if(s>0)return upper?std::min(1.,up(1-hBound(-s,e,false))):std::max(0.,down(1-hBound(-s,e,true)));
    const double t=s/e;const Range a{std::max(0.,down(1+down(t))),std::min(1.,up(1+up(t)))};
    const auto z=positiveProduct(a,{down(pi),up(pi)}),z2=positiveProduct(z,z);
    auto term=positiveDivide(z2,6);Range sum=term,even{},odd{};
    // For 0<=a<=1 the alternating terms decrease. Directed sums through
    // 16/17 terms enclose the complete analytic series without libm sine.
    for(int k=2;k<=17;++k){term=positiveDivide(positiveProduct(term,z2),(2.*k)*(2*k+1));
        sum=add(sum,multiply(term,k%2?1.:-1.));if(k==16)even=sum;if(k==17)odd=sum;}
    const Range series{std::max(0.,even.lower),std::max(0.,odd.upper)};
    const auto result=positiveProduct(positiveDivide(a,2),series);
    // This intersection applies only to a mathematical enclosure; computed
    // force-color values are never clipped.
    return upper?std::min(1.,result.upper):std::max(0.,result.lower);
}
double normBound(double x,double y,bool upper)
{
    const double scale=std::max(x,y);if(scale==0)return 0;
    const double a=upper?up(x/scale):std::max(0.,down(x/scale));
    const double b=upper?up(y/scale):std::max(0.,down(y/scale));
    const double sum=upper?up(up(a*a)+up(b*b)):std::max(0.,down(down(a*a)+down(b*b)));
    const double root=upper?up(std::sqrt(sum)):std::max(0.,down(std::sqrt(sum)));
    return upper?up(scale*root):std::max(0.,down(scale*root));
}
Range derivativeRange(const LevelSetHermiteGeometry::Snapshot& scalar,int cell,
    const std::array<LAMP::Geometry::PolynomialAffineCoordinate,2>& physical,int axis)
{
    const auto b=scalar.nativeCellBounds(cell);const double step=(axis==0?b[1].x-b[0].x:b[1].y-b[0].y)/4;
    require(step>0&&std::isfinite(12*step),"unrepresentable color polynomial derivative scale");
    Range center{},tangent{},time{},remainder{};
    constexpr double shifts[4]={1,-1,2,-2},weights[4]={8,-8,-1,1};
    for(int k=0;k<4;++k){auto coords=physical;const double old=coords[axis].center,delta=step*shifts[k],sum=old+delta;
        // TwoSum retains the exact offset, even when it is below the global
        // coordinate ULP. Synthetic offsets refer only to this polynomial's
        // mathematical extension, never neighboring physical ghost data.
        const double z=sum-old,error=(old-(sum-z))+(delta-z);coords[axis].center=sum;
        coords[axis].remainder=add(coords[axis].remainder,{error,error});
        const auto f=scalar.nativeQuarticAffineEnclosure(cell,coords);
        center=add(center,multiply({f.center,f.center},weights[k]));
        tangent=add(tangent,multiply({f.tangent,f.tangent},weights[k]));
        time=add(time,multiply({f.time,f.time},weights[k]));
        remainder=add(remainder,multiply(f.remainder,weights[k]));}
    // For EVERY degree<=4 polynomial, [8(P(x+a)-P(x-a))-
    // (P(x+2a)-P(x-2a))]/(12a) is EXACTLY P'(x). Common affine
    // parameters are retained; all roundoff/nonlinear remainder stays outward.
    // This is not a sampled finite-difference approximation or a Q4 refit.
    const double variation=std::nextafter(std::max(std::abs(tangent.lower),std::abs(tangent.upper))+
        std::max(std::abs(time.lower),std::abs(time.upper)),std::numeric_limits<double>::infinity());
    const auto numerator=add(add(center,{-variation,variation}),remainder);
    const double dl=std::nextafter(12*step,0.),du=std::nextafter(12*step,std::numeric_limits<double>::infinity());
    require(dl>0,"unresolved color derivative denominator");
    const std::array<double,4> values{{numerator.lower/dl,numerator.lower/du,numerator.upper/dl,numerator.upper/du}};
    return {std::nextafter(*std::min_element(values.begin(),values.end()),-std::numeric_limits<double>::infinity()),
        std::nextafter(*std::max_element(values.begin(),values.end()),std::numeric_limits<double>::infinity())};
}
// Enclose the actual native subbox without sampling or extending a
// different macro piece. The source endpoints are inside this leaf;
// outward quotient tails may be intersected with its exact [0,1] chart.
std::array<Vector2D,2> macroLocalBox(const LevelSetHermiteGeometry::Snapshot& scalar,int cell,Vector2D a,Vector2D b)
{
    namespace E=LAMP::Geometry::Detail::ExactExpansion;
    const auto bounds=scalar.nativeCellBounds(cell);
    std::array<Vector2D,2> result;
    for(int axis=0;axis<2;++axis) {
        const double lo=axis?bounds[0].y:bounds[0].x,hi=axis?bounds[1].y:bounds[1].x;
        const double first=axis?a.y:a.x,last=axis?b.y:b.x;
        require(first>=lo&&last<=hi&&first<=last,"macro color box leaves native support");
        const auto d=E::enclosure(E::subtract(E::fromDouble(hi),E::fromDouble(lo)));
        const auto n0=E::enclosure(E::subtract(E::fromDouble(first),E::fromDouble(lo)));
        const auto n1=E::enclosure(E::subtract(E::fromDouble(last),E::fromDouble(lo)));
        require(d.lo>0,"macro color chart width unresolved");
        const double a0=std::max(0.,std::nextafter(n0.lo/d.hi,-std::numeric_limits<double>::infinity()));
        const double a1=std::min(1.,std::nextafter(n1.hi/d.lo,std::numeric_limits<double>::infinity()));
        (axis?result[0].y:result[0].x)=a0;(axis?result[1].y:result[1].x)=a1;
    }
    return result;
}
BoxBounds enclose(const LevelSetHermiteGeometry::Snapshot& scalar,int cell,Vector2D a,Vector2D b,double e,double inputScale)
{
    using Form=LAMP::Geometry::PolynomialAffineCoordinate;
    const Vector2D center=a+(b-a)*.5;
    const double rx=std::nextafter(std::max(center.x-a.x,b.x-center.x),std::numeric_limits<double>::infinity()),
                 ry=std::nextafter(std::max(center.y-a.y,b.y-center.y),std::numeric_limits<double>::infinity());
    require(rx>0&&ry>0,"unrepresentable color integration box");
    const std::array<Form,2> coords{{{center.x,rx,0,{0,0}},{center.y,0,ry,{0,0}}}};
    using PRange=LAMP::Geometry::PolynomialVelocityRange;
    PRange p;Range gx,gy;
    if(scalar.usesC2Macro()) {
        const auto local=macroLocalBox(scalar,cell,a,b);
        const auto value=scalar.nativeC2MacroScalarBounds(cell,local[0],local[1]);
        const auto x=scalar.nativeC2MacroScalarBounds(cell,local[0],local[1],1,0);
        const auto y=scalar.nativeC2MacroScalarBounds(cell,local[0],local[1],0,1);
        p={value[0],value[1]};gx={x[0],x[1]};gy={y[0],y[1]};
    } else {
        const auto f=scalar.nativeQuarticAffineEnclosure(cell,coords);
        p=f.enclose({-1,1},{-1,1});
        gx=derivativeRange(scalar,cell,coords,0);gy=derivativeRange(scalar,cell,coords,1);
    }
    const double round=2048*std::numeric_limits<double>::epsilon()*inputScale;
    const auto upper=[](Range r){return std::max(std::abs(r.lower),std::abs(r.upper));};
    const auto lower=[](Range r){return r.lower>0?r.lower:r.upper<0?-r.upper:0.;};
    const double ux=up(upper(gx)+round),uy=up(upper(gy)+round);
    const double lx=std::max(0.,down(lower(gx)-round)),ly=std::max(0.,down(lower(gy)-round));
    BoxBounds out;out.phiLower=p.lower;out.phiUpper=p.upper;
    out.gradientUpper=normBound(ux,uy,true);
    out.gradientLower=normBound(lx,ly,false);
    require(std::isfinite(out.gradientUpper),"unrepresentable color gradient bound");
    const double outerBand=up(e*out.gradientUpper);
    if(p.lower>=outerBand&&p.lower>0){out.pure=1;out.cLower=out.cUpper=0;return out;}
    if(p.upper<=-outerBand&&p.upper<0){out.pure=-1;out.cLower=out.cUpper=1;return out;}
    if(p.lower>=0){out.cUpper=hBound(-std::max(0.,down(p.lower/out.gradientUpper)),e,true);
        out.cLower=out.gradientLower>0?hBound(-up(p.upper/out.gradientLower),e,false):0;}
    else if(p.upper<=0){out.cLower=hBound(std::max(0.,down(-p.upper/out.gradientUpper)),e,false);
        out.cUpper=out.gradientLower>0?hBound(up(-p.lower/out.gradientLower),e,true):1;}
    out.allBand=out.gradientLower>0&&std::max(std::abs(p.lower),std::abs(p.upper))<down(e*out.gradientLower);
    return out;
}
struct Rule {long double liquid=0,gas=0,band=0;bool allZero=true,allOne=true;};
constexpr double x8[8]={-.96028985649753623,-.79666647741362674,-.52553240991632899,-.18343464249564980,
 .18343464249564980,.52553240991632899,.79666647741362674,.96028985649753623};
constexpr double w8[8]={.10122853629037626,.22238103445337447,.31370664587788729,.36268378337836198,
 .36268378337836198,.31370664587788729,.22238103445337447,.10122853629037626};
constexpr double x16[16]={-.98940093499164993,-.94457502307323258,-.86563120238783174,-.75540440835500303,
 -.61787624440264375,-.45801677765722739,-.28160355077925891,-.09501250983763744,
 .09501250983763744,.28160355077925891,.45801677765722739,.61787624440264375,
 .75540440835500303,.86563120238783174,.94457502307323258,.98940093499164993};
constexpr double w16[16]={.02715245941175409,.06225352393864789,.09515851168249278,.12462897125553387,
 .14959598881657673,.16915651939500254,.18260341504492359,.18945061045506850,
 .18945061045506850,.18260341504492359,.16915651939500254,.14959598881657673,
 .12462897125553387,.09515851168249278,.06225352393864789,.02715245941175409};
double query(const LevelSetHermiteGeometry::Snapshot& scalar,int cell,Vector2D p,double e,double scale,double resolution,
    HydroRegularizedColorCellReport* report,bool* inBand=nullptr)
{
    if(inBand)*inBand=false;
    const auto j=scalar.scalarJetInNativeCell(cell,p);const double g=std::hypot(j.dx,j.dy);
    require(std::isfinite(j.value)&&std::isfinite(g),"nonfinite color scalar/gradient in cell "+std::to_string(cell));
    if(j.value>0&&j.value>=e*g)return 0;
    if(j.value<0&&-j.value>=e*g)return 1;
    require(g>resolution*scale,"degenerate/unresolved gradient inside color force band, cell "+std::to_string(cell));
    if(inBand)*inBand=true;
    if(report)report->minimumResolvedGradientRatio=std::min(report->minimumResolvedGradientRatio,g/scale);
    return HydroRegularizedColor::heaviside(-j.value/g,e);
}
class Integrator
{
public:
    Integrator(const LevelSetHermiteGeometry::Snapshot& s,int c,bool axis,HydroRegularizedColorConfig f,double scale,
        double cellWidth,HydroRegularizedColorCellReport& r):scalar(s),cell(c),rz(axis),cfg(f),inputScale(scale),h(cellWidth),report(r){}
    Rule integrate(Vector2D a,Vector2D b,const BoxBounds& initialBounds)
    {
        struct Panel{Vector2D a,b;int depth=0;Rule rule;long double error=0;BoxBounds bound;bool active=true,unresolved=false;};
        std::vector<Panel> panels;
        // Certification is a separate mandatory condition, never a fabricated
        // quadrature error. Unresolved boxes have priority over numerical error.
        std::priority_queue<std::pair<std::pair<bool,long double>,std::size_t>> pending;
        const auto evaluate=[&](Vector2D lo,Vector2D hi,int depth,const BoxBounds* retainedBounds=nullptr){
            Panel p;p.a=lo;p.b=hi;p.depth=depth;
            if(retainedBounds)p.bound=*retainedBounds;
            else{
                require(report.queries+report.enclosureQueries+(scalar.usesC2Macro()?3:9)<=cfg.maxQueriesPerCell,
                    "color geometry certification query limit: unresolved gradient/band, cell "+std::to_string(cell));
                p.bound=enclose(scalar,cell,lo,hi,cfg.epsilon,inputScale);report.enclosureQueries+=scalar.usesC2Macro()?3:9;
            }
            const auto cap=volume(lo,hi,rz);report.deepestPanel=std::max(report.deepestPanel,depth);
            if(p.bound.pure)p.rule=p.bound.pure<0?Rule{cap,0,0,false,true}:Rule{0,cap,0,true,false};
            else{
                require(!p.bound.allBand||cfg.epsilon>=h/2,"color force band underresolved by native cell "+std::to_string(cell));
                p.unresolved=p.bound.gradientLower<=up(cfg.relativeGradientResolution*inputScale)||cfg.epsilon<h/2;
                if(p.unresolved)return p;
                const auto coarse=gauss(lo,hi,x8,w8,8);p.rule=gauss(lo,hi,x16,w16,16);
                p.error=std::abs(p.rule.liquid-coarse.liquid);
                if(p.rule.allZero||p.rule.allOne)p.error=std::max(p.error,cap*(p.bound.cUpper-p.bound.cLower));}
            return p;};
        // build() already certified this exact cell/epsilon/input-scale box.
        panels.push_back(evaluate(a,b,0,&initialBounds));pending.push({{panels[0].unresolved,panels[0].error},0});
        std::size_t unresolved=panels[0].unresolved;
        const long double budget=cfg.relativeTolerance*volume(a,b,rz);long double error=panels[0].error;
        while(error>budget||unresolved){const auto id=pending.top().second;pending.pop();auto p=panels[id];
            require(p.depth<cfg.maxDepth,p.unresolved?"color unresolved gradient/band certification depth exhausted, cell "+std::to_string(cell):
                "color cubature depth exhausted, cell "+std::to_string(cell));
            const int axis=(p.b.x-p.a.x)>=(p.b.y-p.a.y)?0:1;const auto mid=p.a+(p.b-p.a)*.5;
            Vector2D left=p.b,right=p.a;
            if(axis==0){require(mid.x>p.a.x&&mid.x<p.b.x,"unrepresentable color x subdivision");left.x=right.x=mid.x;}
            else{require(mid.y>p.a.y&&mid.y<p.b.y,"unrepresentable color y subdivision");left.y=right.y=mid.y;}
            panels[id].active=false;error-=p.error;unresolved-=p.unresolved;
            for(const auto box:std::array<std::array<Vector2D,2>,2>{{{p.a,left},{right,p.b}}}){
                auto q=evaluate(box[0],box[1],p.depth+1);error+=q.error;unresolved+=q.unresolved;
                pending.push({{q.unresolved,q.error},panels.size()});panels.push_back(std::move(q));}
        }
        Rule sum;
        for(const auto& p:panels)if(p.active){++report.acceptedPanels;sum.liquid+=p.rule.liquid;sum.gas+=p.rule.gas;
            sum.band+=p.rule.band;sum.allZero=sum.allZero&&p.rule.allZero;sum.allOne=sum.allOne&&p.rule.allOne;
            report.absoluteErrorEstimate+=p.error;
            if(!p.bound.pure){require(!p.unresolved,"unresolved color panel accepted");
                report.minimumCertifiedGradientRatio=std::min(report.minimumCertifiedGradientRatio,p.bound.gradientLower/inputScale);
                const auto cap=volume(p.a,p.b,rz);report.transitionVolumeEstimate+=p.rule.band;
                if(p.bound.allBand)report.transitionVolumeLower+=cap;report.transitionVolumeUpper+=cap;}}
        require(report.absoluteErrorEstimate<=budget*(1+64*std::numeric_limits<double>::epsilon()),"color aggregate error ledger exceeds budget");
        return sum;
    }
private:
    Rule gauss(Vector2D a,Vector2D b,const double* nodes,const double* weights,int n)
    {
        require(report.queries+report.enclosureQueries+static_cast<std::size_t>(n*n)<=cfg.maxQueriesPerCell,"color cubature query limit, cell "+std::to_string(cell));
        Rule out;long double total=0;
        for(int i=0;i<n;++i)for(int j=0;j<n;++j){const Vector2D p{a.x+(b.x-a.x)*(nodes[i]+1)/2,a.y+(b.y-a.y)*(nodes[j]+1)/2};
            bool inBand=false;
            const auto c=query(scalar,cell,p,cfg.epsilon,inputScale,cfg.relativeGradientResolution,&report,&inBand);++report.queries;
            const long double w=weights[i]*weights[j]*(rz?CartesianMeasure::radialWeight(p.x):1.);
            out.liquid+=w*c;out.gas+=w*(1-c);total+=w;
            if(inBand){require(cfg.epsilon>=h/2,"color force band underresolved by native cell "+std::to_string(cell));out.band+=w;}
            out.allZero=out.allZero&&c==0;out.allOne=out.allOne&&c==1;}
        require(total>0,"zero positive color quadrature measure");
        const auto cap=volume(a,b,rz);
        // Standard positive average rule normalized to the TRUE box volume.
        // Companion sums make constants/boundedness exact without clipping;
        // these are force-color weights, never material inventory corrections.
        const auto sum=out.liquid+out.gas;out.liquid=cap*(out.liquid/sum);out.gas=cap-out.liquid;out.band=cap*out.band/total;
        require(out.liquid>=0&&out.gas>=0,"color positive cubature lost boundedness");return out;
    }
    const LevelSetHermiteGeometry::Snapshot& scalar;int cell;bool rz;HydroRegularizedColorConfig cfg;
    double inputScale,h;HydroRegularizedColorCellReport& report;
};
}
struct HydroRegularizedColor::Data
{
    Data(HydroPhaseGeometry p,HydroRegularizedColorConfig c):phase(std::move(p)),cfg(c){}
    HydroPhaseGeometry phase;HydroRegularizedColorConfig cfg;
    std::vector<double> color;std::vector<HydroRegularizedColorCellReport> cells;HydroRegularizedColorReport report;
};
HydroRegularizedColor::HydroRegularizedColor(std::shared_ptr<const Data> d):data_(std::move(d)){}
double HydroRegularizedColor::heaviside(double s,double e)
{
    width(e);require(std::isfinite(s),"nonfinite signed length");if(s<=-e)return 0;if(s>=e)return 1;
    const double t=s/e;const double value=t<=0?lowerH(t):1-lowerH(-t);
    require(value>=0&&value<=1,"regularized kernel lost boundedness");return value;
}
double HydroRegularizedColor::delta(double s,double e)
{
    width(e);require(std::isfinite(s),"nonfinite delta length");if(std::abs(s)>=e)return 0;
    const double value=std::cos((pi/2)*(s/e)),result=value*value/e;
    require(std::isfinite(result),"unrepresentable regularized delta");return result;
}
HydroRegularizedColor HydroRegularizedColor::build(const QuadtreeMesh& mesh,HydroPhaseGeometry phase,HydroRegularizedColorConfig cfg)
{
    width(cfg.epsilon);require(std::isfinite(cfg.relativeTolerance)&&cfg.relativeTolerance>0&&cfg.relativeTolerance<1,
        "invalid color relative integration tolerance");
    require(std::isfinite(cfg.relativeGradientResolution)&&cfg.relativeGradientResolution>0&&cfg.relativeGradientResolution<1,
        "invalid color relative gradient resolution");
    require(cfg.maxDepth>=0&&cfg.maxDepth<=40&&cfg.maxQueriesPerCell>=320,"invalid color work limits");
    phase.validateFor(mesh,phase.geometryIdentity(),phase.timeLayer());
    require(phase.publication().state().primary().transport==Interface::HermiteGalsTransport::FirstJet,"color requires Pure FirstJet publication");
    auto d=std::make_shared<Data>(std::move(phase),cfg);const auto& layout=d->phase.layout();
    const auto& scalar=d->phase.publication().state().geometry().snapshot();const bool rz=layout.domain().axisymmetric();
    d->color.resize(layout.pressureCount());d->cells.resize(layout.pressureCount());
    std::vector<std::exception_ptr> failures(layout.pressureCount());
    const auto prepareCell=[&](int p){
        try{
            LAMP::Geometry::Detail::ExactExpansion::validateActualEnvironment();
            const int cell=layout.nativeCell(p);const auto b=scalar.nativeCellBounds(cell);
            const double hx=b[1].x-b[0].x,hy=b[1].y-b[0].y;auto& r=d->cells[p];
            const auto& phaseMeasure=d->phase.cellVolumes()[p];
            if(phaseMeasure.liquid>0&&phaseMeasure.gas>0)
                require(cfg.epsilon>=std::max(hx,hy)/2,"color force band underresolved by cut native cell "+std::to_string(cell));
            using Form=LAMP::Geometry::PolynomialAffineCoordinate;
            const auto center=b[0]+(b[1]-b[0])*.5;
            const auto range=scalar.usesC2Macro()?scalar.nativeC2MacroScalarBounds(cell,{0,0},{1,1}):[&]{
                const auto form=scalar.nativeQuarticAffineEnclosure(cell,std::array<Form,2>{{{center.x,hx/2,0,{0,0}},{center.y,0,hy/2,{0,0}}}});
                const auto bound=form.enclose({-1,1},{-1,1});return std::array<double,2>{{bound.lower,bound.upper}};}();
            r.gradientScale=std::max(std::abs(range[0]),std::abs(range[1]))/std::min(hx,hy);
            require(std::isfinite(r.gradientScale)&&r.gradientScale>0,"unresolved/zero color scalar scale in cell "+std::to_string(cell));
            const auto bounds=enclose(scalar,cell,b[0],b[1],cfg.epsilon,r.gradientScale);r.certifiedPure=bounds.pure!=0;
            r.enclosureQueries+=scalar.usesC2Macro()?4:10; // Actual input-scale and scalar/gradient enclosure calls.
            Integrator integral(scalar,cell,rz,cfg,r.gradientScale,std::max(hx,hy),r);const auto result=integral.integrate(b[0],b[1],bounds);
            const auto total=result.liquid+result.gas;require(total>0,"unresolved physical color cell volume");
            const double color=static_cast<double>(result.liquid/total);require(color>=0&&color<=1,"integrated color outside [0,1]");d->color[p]=color;
        }catch(...){failures[p]=std::current_exception();}
    };
#ifdef _OPENMP
    const bool runParallel=layout.pressureCount()>1&&!omp_in_parallel()&&omp_get_max_threads()>1;
#pragma omp parallel for schedule(dynamic,1) if(runParallel)
#endif
    for(int p=0;p<layout.pressureCount();++p)prepareCell(p);
    // Distinct cells own all result/report/exception slots. Retain native
    // pressure order for both the first failure and every report operation.
    for(int p=0;p<layout.pressureCount();++p){
        if(failures[p])std::rethrow_exception(failures[p]);
        const auto& r=d->cells[p];const double color=d->color[p];
        const auto cap=d->phase.metrics().cellVolumes()[p];const double difference=std::abs(color-d->phase.liquidFractions()[p]);
        d->report.physicalVolume+=cap;d->report.sourceQueries+=r.queries;d->report.enclosureQueries+=r.enclosureQueries;d->report.pureCells+=r.certifiedPure;
        d->report.absoluteErrorEstimate+=r.absoluteErrorEstimate;d->report.transitionVolumeEstimate+=r.transitionVolumeEstimate;
        d->report.transitionVolumeLower+=r.transitionVolumeLower;d->report.transitionVolumeUpper+=r.transitionVolumeUpper;
        d->report.colorMinusSharpL1+=cap*difference;d->report.colorMinusSharpMax=std::max(d->report.colorMinusSharpMax,difference);
    }
    return HydroRegularizedColor(std::move(d));
}
const HydroPhaseGeometry& HydroRegularizedColor::phaseGeometry() const noexcept{return data_->phase;}
const HydroRegularizedColorConfig& HydroRegularizedColor::config() const noexcept{return data_->cfg;}
const std::vector<double>& HydroRegularizedColor::colors() const noexcept{return data_->color;}
const std::vector<HydroRegularizedColorCellReport>& HydroRegularizedColor::cellReports() const noexcept{return data_->cells;}
const HydroRegularizedColorReport& HydroRegularizedColor::report() const noexcept{return data_->report;}
Interface::HermiteGalsIdentity HydroRegularizedColor::geometryIdentity() const{return phaseGeometry().geometryIdentity();}
Interface::HermiteGalsTimeLayer HydroRegularizedColor::timeLayer() const{return phaseGeometry().timeLayer();}
void HydroRegularizedColor::validateFor(const QuadtreeMesh& mesh,Interface::HermiteGalsIdentity id,Interface::HermiteGalsTimeLayer clock) const
{phaseGeometry().validateFor(mesh,id,clock);}
void HydroRegularizedColor::validateFor(const QuadtreeMesh& mesh,const HydroPhaseGeometry& p) const
{
    phaseGeometry().validateFor(mesh,p.geometryIdentity(),p.timeLayer());p.validateFor(mesh,p.geometryIdentity(),p.timeLayer());
    require(phaseGeometry().layout().sameIdentity(p.layout()) && &phaseGeometry().metrics().divergence().matrix()==&p.metrics().divergence().matrix(),
        "different color layout/divergence source");
}
double HydroRegularizedColor::pointColor(int p,Vector2D point) const
{
    require(p>=0&&p<phaseGeometry().layout().pressureCount(),"color pressure index out of range");
    return query(phaseGeometry().publication().state().geometry().snapshot(),phaseGeometry().layout().nativeCell(p),point,
        config().epsilon,cellReports()[p].gradientScale,config().relativeGradientResolution,nullptr);
}
}
