#include "Core/Geometry/QuadraticCellCut.h"
#include "Core/Geometry/QuadraticNativeSingleArc.h"
#include "Core/Geometry/QuadraticSurfaceMeasure.h"
#include "Core/Mesh/CartesianMeasure.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
namespace Q=LAMP::Geometry::QuadraticCellCut;
namespace S=LAMP::Geometry::QuadraticSurfaceMeasure;
void check(bool b,const char*m){if(!b)throw std::runtime_error(m);}
void near(long double a,long double b,const char*m){check(b>0&&std::abs(a-b)<=1e-12L*b,m);}
void cachedSelectedPhase(const Q::Rectangle& box,const Q::FixedPortVolumeResult& solved,bool rz)
{
    const auto independent=Q::graphMeasure(box,solved.graph,rz);
    check(solved.graph.volume==independent.volume,"cached selected-phase volume differs from independent measure");
    for(int k=0;k<3;++k)check(solved.selectedPhaseDerivative[k]==independent.derivative[k],
        "cached normalized selected-phase derivative differs from independent measure");
}
template<class F>void rejects(F f){bool bad=false;try{f();}catch(const std::invalid_argument&){bad=true;}check(bad,"stale encoded polynomial accepted");}
void affineNormalSweepCases()
{
    // An actual dyadic coarse-cell scale and its two 2:1 native subfaces.
    const Q::Rectangle donor{{.25,.25},{.5,.5}};Q::Shape shape;shape.axis=1;
    const auto graph=Q::fromFixedPorts(donor,shape,{.25,31.L/64},{.5,31.L/64},0);
    const Q::Point begin{.5,.25},middle{.5,.375},end{.5,.5};
    const long double duration=1.L/1024,mean=1.25L*duration,variation=.25L*duration,H=.25L,f=15.L/16;
    const auto affine=Q::linearSignedNormalSweepMeasure(donor,graph,0,begin,end,mean,variation,true,false);
    const auto frozen=Q::axisStripMeasure(donor,graph,0,begin,end,mean,true,false);
    const auto primitive=[&](long double s){return (mean-variation)*s+variation*s*s;};
    near(affine.total,H*mean,"affine unique mean gross");
    near(affine.liquid,H*primitive(f),"affine independently integrated liquid");
    near(affine.gas,H*(primitive(1)-primitive(f)),"affine independently integrated gas");
    check(std::abs((affine.liquid-frozen.liquid)+15.L/4096*duration)<=1e-12L*(15.L/4096*duration),
        "affine phase covariance is not -15/4096");
    const auto lower=Q::linearSignedNormalSweepMeasure(donor,graph,0,begin,end,mean,variation,true,false,
        std::array<Q::Point,2>{{begin,middle}});
    const auto upper=Q::linearSignedNormalSweepMeasure(donor,graph,0,begin,end,mean,variation,true,false,
        std::array<Q::Point,2>{{end,middle}});
    near(lower.total+upper.total,affine.total,"affine actual subface gross additivity");
    near(lower.liquid+upper.liquid,affine.liquid,"affine actual subface liquid additivity");
    near(lower.gas+upper.gas,affine.gas,"affine actual subface gas additivity");
    const auto reversed=Q::linearSignedNormalSweepMeasure(donor,graph,0,end,begin,mean,-variation,true,false);
    check(reversed.total==affine.total&&reversed.liquid==affine.liquid&&reversed.gas==affine.gas,"reversed original face changed affine law");
    const auto unsignedSweep=Q::linearNormalSweepMeasure(donor,graph,0,begin,end,mean-variation,mean+variation,true,false);
    check(unsignedSweep.total==affine.total&&unsignedSweep.liquid==affine.liquid&&unsignedSweep.gas==affine.gas,"endpoint and factored affine laws disagree");
    // Same exact polygon with true RZ radius: vertical donor-back-edge radius
    // varies with depth, so its physical measure includes -d(s)^2/2.
    const auto radial=Q::linearSignedNormalSweepMeasure(donor,graph,0,begin,end,mean,variation,true,true);
    const auto radialPrimitive=[&](long double s){const auto a=mean-variation,b=2*variation;
        return .5L*(a*s+b*s*s/2)-(a*a*s+a*b*s*s+b*b*s*s*s/3)/2;};
    near(radial.total,CartesianMeasure::twoPi*H*radialPrimitive(1),"RZ affine wedge total radial moment");
    near(radial.liquid,CartesianMeasure::twoPi*H*radialPrimitive(f),"RZ affine wedge independent liquid");
    near(radial.gas,CartesianMeasure::twoPi*H*(radialPrimitive(1)-radialPrimitive(f)),"RZ affine wedge independent gas");
    // Signed zero is at EXACT parameter 3/7, not a binary world coordinate.
    const long double signedMean=1.L/512,signedVariation=7.L/512;
    const auto positive=Q::linearSignedNormalSweepMeasure(donor,graph,0,begin,end,signedMean,signedVariation,true,false);
    const auto negative=Q::linearSignedNormalSweepMeasure(donor,graph,0,{.25,.25},{.25,.5},signedMean,signedVariation,false,false);
    near(positive.total,H*(signedMean+signedVariation)*(signedMean+signedVariation)/(4*signedVariation),"exact signed positive triangle");
    near(negative.total,H*(signedVariation-signedMean)*(signedVariation-signedMean)/(4*signedVariation),"exact signed negative triangle");
    near(positive.total-negative.total,H*signedMean,"signed gross keeps native mean");
    const auto positiveLo=Q::linearSignedNormalSweepMeasure(donor,graph,0,begin,end,signedMean,signedVariation,true,false,
        std::array<Q::Point,2>{{begin,middle}});
    const auto positiveHi=Q::linearSignedNormalSweepMeasure(donor,graph,0,begin,end,signedMean,signedVariation,true,false,
        std::array<Q::Point,2>{{middle,end}});
    near(positiveLo.total+positiveHi.total,positive.total,"signed-zero subface query additivity");
    check(positive.liquid>0&&positive.gas>0&&negative.liquid>0&&negative.gas==0,"signed donor phases were replaced by a complement");
    // Mean survives even when native mean +/- variation and world zero would
    // both round away this very small triangle on the lower half-query.
    const auto tinyMean=std::ldexp(1.L,-80),largeVariation=1.L/64;
    const auto tiny=Q::linearSignedNormalSweepMeasure(donor,graph,0,begin,end,tinyMean,largeVariation,true,false,
        std::array<Q::Point,2>{{begin,middle}});
    near(tiny.total,H*tinyMean*tinyMean/(4*largeVariation),"factored mean lost below endpoint/world ulp");
    check(tiny.liquid==tiny.total&&tiny.gas==0,"thin factored triangle lost original phase");
    const auto absent=Q::linearSignedNormalSweepMeasure(donor,graph,0,begin,end,-mean,0,true,false);
    check(absent.total==0&&absent.liquid==0&&absent.gas==0,"absent signed part gained stock");
    const auto empty=Q::linearNormalSweepMeasure(donor,graph,0,begin,end,0,0,true,false);
    check(empty.total==0&&empty.liquid==0&&empty.gas==0,"zero endpoint sweep gained stock");
    for(bool rz:{false,true})for(bool direction:{false,true}){
        const auto x=direction?.5L:.25L;const Q::Point a{x,.25},b{x,.5};const auto depth=1.L/4096;
        const auto old=Q::axisStripMeasure(donor,graph,0,a,b,depth,direction,rz);
        const auto same=Q::linearSignedNormalSweepMeasure(donor,graph,0,a,b,direction?depth:-depth,0,direction,rz);
        check(old.total==same.total&&old.liquid==same.liquid&&old.gas==same.gas,"constant-depth original strip compatibility");
        check(same.coveredPhases==0,"coverage unexpectedly paid for whole-donor integration");
        const auto covered=Q::linearNormalSweepMeasure(donor,graph,0,a,b,.25,.25,direction,rz,std::nullopt,true);
        check(covered.coveredPhases==3,"explicit independent donor phase coverage");}
    auto flipped=graph;flipped.shape.sign=-1;
    const auto opposite=Q::linearSignedNormalSweepMeasure(donor,flipped,0,begin,end,mean,variation,true,false);
    check(opposite.liquid==affine.gas&&opposite.gas==affine.liquid,"affine phase orientation changed geometry");
    const Q::Rectangle axisBox{{0,0},{.125,.25}};const auto axis=Q::fromAxisPorts(axisBox,1,{0,.1},{.125,.15});
    const Q::Point axisA{0,.25},axisB{.125,.25};const auto axisDepth=1.L/64;
    const auto axisOld=Q::axisStripMeasure(axisBox,axis,1,axisA,axisB,axisDepth,true,true);
    const auto axisSame=Q::linearNormalSweepMeasure(axisBox,axis,1,axisA,axisB,axisDepth,axisDepth,true,true);
    check(axisOld.total==axisSame.total&&axisOld.liquid==axisSame.liquid&&axisOld.gas==axisSame.gas,"axis donor constant sweep mismatch");
    const auto axisAffine=Q::linearSignedNormalSweepMeasure(axisBox,axis,1,axisA,axisB,1.L/64,1.L/128,true,true);
    near(axisAffine.total,CartesianMeasure::twoPi*.125L*.125L*((1.L/64)/2+(1.L/128)/6),"axis RZ affine radius-weighted measure");
    rejects([&]{(void)Q::linearNormalSweepMeasure(donor,graph,0,begin,end,-mean,mean,true,false);});
    rejects([&]{(void)Q::linearSignedNormalSweepMeasure(donor,graph,0,begin,end,1,0,true,false);});
    rejects([&]{(void)Q::linearSignedNormalSweepMeasure(donor,graph,0,begin,end,mean,std::numeric_limits<long double>::infinity(),true,false);});
    rejects([&]{(void)Q::linearSignedNormalSweepMeasure(donor,graph,0,begin,end,mean,variation,true,false,
        std::array<Q::Point,2>{{begin,{.5,.75}}});});
    rejects([&]{(void)Q::linearSignedNormalSweepMeasure(donor,graph,0,begin,end,mean,variation,true,false,
        std::array<Q::Point,2>{{begin,{.375,.375}}});});
    auto legacy=graph;legacy.exact.reset();rejects([&]{(void)Q::linearSignedNormalSweepMeasure(donor,legacy,0,begin,end,mean,variation,true,false);});
    auto stale=graph;stale.offset=std::nextafter(stale.offset,1.L);rejects([&]{(void)Q::linearSignedNormalSweepMeasure(donor,stale,0,begin,end,mean,variation,true,false);});
}
int main(){try{
    affineNormalSweepCases();
    // Actual repair8 cap scale; the raw shared ports remain authority even
    // when a centered offset cannot represent their small bubble correction.
    const Q::Rectangle box{{.25,0},{.375,.125}};const long double half=.03125L,h=.125L,c2=2.80626114169708585e-5L/2;
    const Q::Point first{.3125L-half,.125L},last{.3125L+half,.125L};Q::Shape shape;shape.axis=1;
    auto graph=Q::fromFixedPorts(box,shape,first,last,c2*h);check(graph.exact!=nullptr,"missing exact polynomial");
    check(Q::implicitValue(box,graph,first)==0&&Q::implicitValue(box,graph,last)==0,"factory ports changed exact polynomial");
    auto gas=graph;gas.shape.sign=-1;
    for(bool rz:{false,true}){
        const long double expected=4*c2*half*half*half/3*(rz?CartesianMeasure::twoPi*.3125L:1);
        near(Q::graphMeasure(box,gas,rz).volume,expected,"tiny cap own-phase measure");
        if(!rz){namespace N=LAMP::Geometry::QuadraticNativeSingleArc;N::Request request;
            request.box=box;request.shape=graph.shape;request.ports={{first,last}};request.targetPhaseVolume=expected;
            request.volumeLiquid=false;request.branchEta=c2*h;
            const auto prepared=N::prepare(request);
            check(!prepared.value&&(prepared.reason==N::Reason::NearEvent||prepared.reason==N::Reason::UncertainRoundoff)&&!prepared.detail.empty(),
                "same-face cap must preselect its named strict chart, not an uncertified native trial");}
        const auto surface=S::evaluate(box,graph,box,rz);check(surface.patches.size()==1,"cap arc topology");const auto&p=surface.patches.front();
        check(p.begin.x==first.x&&p.begin.y==first.y&&p.end.x==last.x&&p.end.y==last.y,"factory side event lost port bits");
        check(Q::materialCurveData(box,graph,p.lower,p.upper)!=nullptr,"same event material restriction");
        long double sum=0;for(int x=0;x<2;++x)for(int y=0;y<2;++y){const Q::Rectangle child{{.25+x*.0625,y*.0625},{.3125+x*.0625,.0625+y*.0625}};const auto shifted=Q::recenter(box,child,gas);sum+=Q::graphMeasure(child,shifted,rz).volume;}
        near(sum,expected,"exact recenter cap additivity");
        const auto strip=Q::axisStripMeasure(box,graph,1,{.25,.125},{.375,.125},std::ldexp(1.L,-40),true,rz);
        check(strip.liquid>=0&&strip.gas>0&&std::abs(strip.liquid+strip.gas-strip.total)<=4*std::numeric_limits<long double>::epsilon()*strip.total,"encoded thin strip phases");
        // The original API means positive OUTWARD transport: upper-side
        // donor material extends backward; lower-side negative flow extends
        // forward. Keep the actual near plane and depth separate in both.
        for(bool positive:{false,true}){
            const long double nearX=positive?.375L:.25L,depth=std::ldexp(1.L,-40);
            const auto side=Q::axisStripMeasure(box,graph,0,{nearX,0},{nearX,.125L},depth,positive,rz);
            const long double radius=nearX+(positive?-depth:depth)/2;
            const long double gross=depth*.125L*(rz?CartesianMeasure::twoPi*radius:1);
            near(side.total,gross,"directional donor strip gross");
            near(side.liquid,gross,"directional pure donor phase");check(side.gas==0,"directional donor gained gas");
        }
    }
    Q::Shape rotated;rotated.frameNormal=Vector2D{.27488060480878,.9614783685034};
    const Q::Point a{.25,.0625},b{.375,.0625};auto framed=Q::fromFixedPorts(box,rotated,a,b,1e-5L);
    check(Q::implicitValue(box,framed,a)==0&&Q::implicitValue(box,framed,b)==0,"nonunit native frame inverse inconsistent");
    const auto arcs=S::evaluate(box,framed,box,false);check(arcs.patches.size()==1,"framed arc topology");const auto&arc=arcs.patches.front();
    check((arc.begin.x==a.x&&arc.end.x==b.x)||(arc.begin.x==b.x&&arc.end.x==a.x),"framed same source endpoints");
    check(Q::materialCurveData(box,framed,arc.lower,arc.upper)!=nullptr,"framed material identity");
    {namespace N=LAMP::Geometry::QuadraticNativeSingleArc;N::Request request;
        request.box=box;request.shape=framed.shape;request.ports={{a,b}};
        request.targetPhaseVolume=Q::graphMeasure(box,framed,false).volume;request.branchEta=1e-5L;
        const auto prepared=N::prepare(request);
        check(!prepared.value&&(prepared.reason==N::Reason::NearEvent||prepared.reason==N::Reason::UncertainRoundoff||prepared.reason==N::Reason::Unsupported)&&!prepared.detail.empty(),
            "equal-height framed ports must preselect a named strict chart before native use");}
    const Q::Rectangle axisBox{{0,0},{.125,.25}};const auto axis=Q::fromAxisPorts(axisBox,1,{0,.1},{.125,.15});const auto axisSurface=S::evaluate(axisBox,axis,axisBox,true);
    check(axisSurface.patches.size()==1&&axisSurface.patches.front().begin.x==0,"axis root identity");
    const auto axisSample=S::sample(axisBox,axis,axisSurface.patches.front().lower,true);check(axisSample.position.x==0&&axisSample.unitNormal.x==0,"same-curve regular axis sample");
    // Explicit analytic volume API: the selected eta remains exact payload
    // data. XY and N.r==0 RZ preselect a linear equation, not an iteration.
    for(bool rz:{false,true}){
        const auto target=Q::graphMeasure(box,gas,rz).volume;
        const auto recovered=Q::fixedPortVolumeGraph(box,gas.shape,first,last,target,rz,c2*h);
        cachedSelectedPhase(box,recovered,rz);
        const auto reversedRecovered=Q::fixedPortVolumeGraph(box,gas.shape,last,first,target,rz,c2*h);
        cachedSelectedPhase(box,reversedRecovered,rz);
        check(reversedRecovered.eta==recovered.eta&&reversedRecovered.graph.volume==recovered.graph.volume&&
            reversedRecovered.selectedPhaseDerivative==recovered.selectedPhaseDerivative,"reversed ordered cap ports changed the same exact phase geometry");
        check(!recovered.quadratic&&recovered.realRoots==1&&recovered.branchCandidates==1,"analytic linear branch classification");
        near(recovered.graph.volume,target,"analytic cap original own-scale certificate");check(recovered.ownPhaseResidual<=1e-12L,"analytic cap own-phase residual");
        const auto arc=S::evaluate(box,recovered.graph,box,rz);check(arc.patches.size()==1&&arc.patches[0].begin.x==first.x&&arc.patches[0].end.x==last.x,"analytic cap fixed port events");
    }
    const auto framedTarget=Q::graphMeasure(box,framed,true).volume;
    const auto recoveredFrame=Q::fixedPortVolumeGraph(box,framed.shape,a,b,framedTarget,true,1e-5L);
    cachedSelectedPhase(box,recoveredFrame,true);
    const auto reversedFrame=Q::fixedPortVolumeGraph(box,framed.shape,b,a,framedTarget,true,1e-5L);
    cachedSelectedPhase(box,reversedFrame,true);
    check(reversedFrame.eta==recoveredFrame.eta&&reversedFrame.graph.volume==recoveredFrame.graph.volume&&
        reversedFrame.selectedPhaseDerivative==recoveredFrame.selectedPhaseDerivative,"reversed ordered nonunit RZ ports changed the same exact phase geometry");
    check(recoveredFrame.quadratic&&recoveredFrame.realRoots==2&&recoveredFrame.branchCandidates==1,"RZ determinant-aware quadratic branch");
    near(recoveredFrame.graph.volume,framedTarget,"RZ analytic frame volume");
    Q::Shape radial;radial.axis=0;const Q::Point lowPort{.3125L,0},highPort{.3125L,.125L};
    const auto radialGraph=Q::fromFixedPorts(box,radial,lowPort,highPort,.01L);const auto radialTarget=Q::graphMeasure(box,radialGraph,true).volume;
    const auto radialRecovered=Q::fixedPortVolumeGraph(box,radial,lowPort,highPort,radialTarget,true,.01L);
    cachedSelectedPhase(box,radialRecovered,true);
    check(radialRecovered.quadratic,"coordinate axis=0 confused with physical symmetry axis");
    rejects([&]{(void)Q::fixedPortVolumeGraph(box,radial,lowPort,highPort,radialTarget,true,12.5L);});
    rejects([&]{(void)Q::fixedPortVolumeGraph(axisBox,shape,{0,.1},{.125,.15},.001L,true,0);});
    rejects([&]{(void)Q::fixedPortVolumeGraph(box,shape,first,last,0,false,0);});
    auto stale=graph;stale.offset=std::nextafter(stale.offset,std::numeric_limits<long double>::infinity());rejects([&]{(void)Q::graphMeasure(box,stale,false);});
    // Legacy lift preserves Q's raw-dot implicit phase, including frame D.
    Q::Reconstruction legacy;legacy.shape=rotated;legacy.offset=.0125;legacy.shape.slope=.1;legacy.shape.second=.3;
    const auto lifted=Q::recenter(box,box,legacy);for(Q::Point p: {Q::Point{.27,.02},Q::Point{.36,.11}}){const auto old=Q::implicitValue(box,legacy,p),now=Q::implicitValue(box,lifted,p);check((old<0)==(now<0),"legacy lift phase changed");}
    std::cout<<"quadratic exact graph contract PASS\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
