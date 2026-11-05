#include "Modules/Interface/VofSmoothGeometry.h"
#include "Modules/Interface/Detail/VofSmoothPhaseMeasure.h"
#include "Modules/Interface/GeometricVofInitialization.h"
#include "Modules/Interface/GeometricVofTransport.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Geometry/Detail/FiniteVolumeClosure.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
using namespace LAMP::Interface;
using Point=LAMP::Geometry::PlanarCellCut::Point;
unsigned assertions=0;
void check(bool value,const char* message)
{ ++assertions;if(!value)throw std::runtime_error(message); }
long double distance(Point a,Point b) { return std::hypot(a.x-b.x,a.y-b.y); }
std::unique_ptr<QuadtreeMesh> meshFor()
{
    DictionaryMesh d;d.Nx=d.Ny=8;d.maxCellNum=512;d.maxInitialIterations=0;
    return std::make_unique<QuadtreeMesh>(d);
}
template<class F> void rejects(F&& f,VofSmoothGeometryFailure reason,const char* message)
{
    bool caught=false;
    try { f(); }catch(const VofSmoothGeometryError& e){caught=e.reason()==reason;}
    check(caught,message);
}
void stripClose(long double actual,long double expected,const char* message)
{check(std::isfinite(actual)&&std::isfinite(expected)&&expected>=0,message);
    check(expected==0?actual==0:std::abs(actual-expected)<=1e-12L*expected,message);}
void anchoredPhaseMeasure()
{
    // Original full03 step-2 cell490's actual binary64 controls. Its world
    // sample-minus-face arithmetic gave gas=1.1741248985913821e-9 and failed
    // the unchanged 1e-12 own-stock gate, although the represented curve fits.
    VofSmoothArc arc;arc.nativeCell=490;arc.graphAxis=1;
    arc.controlPoints={Point{.3125,.46885206539810631},Point{.31250665282528861,.46882575014295869},
        Point{.31248317496955885,.46880103008542312},Point{.31253035462965201,.4687763113347026},
        Point{.31253702294307995,.46875}};
    GeometricVofPlicCell cell;cell.rectangle={{.3125,.46875},{.34375,.5}};
    const long double stock=1.1741248985892801e-9;
    // Independent Fraction integration of these binary64 Bernstein controls:
    // A_gas=1666985748332338611185717/1419768672255616929676307590021120.
    constexpr long double exactArea=1.1741248985892628637687137344026419987112891218989e-9L;
    // V_gas/pi=-integral((x^2-xlo^2)*y' dt), using the physical RZ measure.
    constexpr long double exactRzOverPi=7.3385693125639993814462835486005196572166875737168e-10L;
    const auto originalControls=arc.controlPoints;const auto originalBox=cell.rectangle;
    const auto measured=LAMP::Interface::Detail::measureSmoothArcPhases(arc,cell,false);
    stripClose(measured[1],exactArea,"actual-control XY microphase differs from independent exact integral");
    stripClose(measured[1],stock,"actual full03 curve still fails the original own-stock gate");
    check(measured[0]>0&&measured[1]>0,"actual-control phase measure lost a positive phase");
    std::array<std::vector<VofSmoothWeightedPoint>,2> quadrature;
    const auto exported=LAMP::Interface::Detail::measureSmoothArcPhases(arc,cell,false,&quadrature);
    check(exported==measured,"root-fit and final-quadrature kernels disagree");
    for(unsigned phase=0;phase<2;++phase){long double sum=0;
        check(!quadrature[phase].empty(),"actual-control positive phase has no physical quadrature");
        for(const auto& q:quadrature[phase]){check(q.weight>0&&std::isfinite(q.weight),"nonpositive actual-control phase weight");sum+=q.weight;}
        check(sum==measured[phase],"actual-control physical quadrature sum differs from measured phase");}
    const auto rz=LAMP::Interface::Detail::measureSmoothArcPhases(arc,cell,true);
    stripClose(rz[1]/std::acos(-1.L),exactRzOverPi,"actual-control RZ microphase lost the true radial Jacobian");
    // Subtract the original exact cell-face anchors. Each difference is
    // representable; translation changes neither the actual XY shape nor stock.
    for(auto& p:arc.controlPoints){p.x-=originalBox.lower.x;p.y-=originalBox.lower.y;}
    cell.rectangle={{0,0},{originalBox.upper.x-originalBox.lower.x,originalBox.upper.y-originalBox.lower.y}};
    const auto translated=LAMP::Interface::Detail::measureSmoothArcPhases(arc,cell,false);
    check(translated==measured,"exact face-origin translation changes XY phase measure");
    for(unsigned k=0;k<5;++k)check(arc.controlPoints[k].x+originalBox.lower.x==originalControls[k].x
        &&arc.controlPoints[k].y+originalBox.lower.y==originalControls[k].y,"actual-control fixture translation changed represented curve");
    std::cout<<"ANCHORED_PHASE_MEASURE,PASS actual_full03_curve=1 original_own_phase_gate=1 exact_XY_RZ_reference=1 translated_origin=1\n";
}
void directSuccess(const VofSmoothArc& original,const GeometricVofPlicCell& cell,
    const LAMP::Interface::Detail::SmoothDirectXYResult& result,const LAMP::Interface::Detail::SmoothDirectXYWork& work)
{
    check(work.candidateConstructions==1&&work.geometryCertifications==1&&work.phaseQuadratureBuilds==1,
        "direct XY success did not construct/certify/integrate exactly once");
    check(work.fitPhaseMeasurements==0&&work.volumeNewtonSteps==0&&work.domainLeaves==256,
        "direct XY used volume fitting or changed the fixed full-depth domain");
    check(result.arc.nodes==original.nodes&&result.arc.graphAxis==original.graphAxis,"direct XY changed source nodes/graph chart");
    for(unsigned i:{0U,1U,3U,4U})check(result.arc.controlPoints[i].x==original.controlPoints[i].x
        &&result.arc.controlPoints[i].y==original.controlPoints[i].y,"direct XY changed a complete endpoint jet");
    const auto analytic=LAMP::Interface::Detail::analyticSmoothArcPhasesXY(result.arc,cell);
    for(unsigned phase=0;phase<2;++phase){stripClose(analytic[phase],cell.inventory[phase],"direct XY analytic area misses own stock");
        stripClose(result.measured[phase],cell.inventory[phase],"direct XY independent quadrature misses own stock");
        long double sum=0;for(const auto& point:result.phaseQuadrature[phase]){
            check(point.weight>0&&std::isfinite(point.weight),"direct XY published a nonpositive physical weight");sum+=point.weight;}
        check(sum==result.measured[phase],"direct XY did not retain its single physical quadrature");}
}
void directAlgebra()
{
    VofSmoothArc arc;arc.nativeCell=7;arc.graphAxis=0;arc.nodes={{13,19}};
    arc.controlPoints={Point{0,.25},Point{.25,.375},Point{.5,.625},Point{.75,.625},Point{1,.75}};
    GeometricVofPlicCell cell;cell.rectangle={{0,0},{1,1}};cell.inventory={{.5,.5}};
    const auto original=arc;const auto stocks=cell.inventory;
    const auto before=LAMP::Interface::Detail::analyticSmoothArcPhasesXY(arc,cell);
    // Independent exact-degree Bernstein reference for this x(t)=t fixture.
    check(before[1]==(.25L+.375L+.625L+.625L+.75L)/5,"analytic direct phase area differs from polynomial mean");
    auto moved=arc;moved.controlPoints[2].x+=1.L/32;moved.controlPoints[2].y-=1.L/64;
    const auto after=LAMP::Interface::Detail::analyticSmoothArcPhasesXY(moved,cell);
    const auto kx=5.L+4*(.75L-.25L),ky=5*(.75L-.25L)+4*(.625L-.375L);
    const auto delta=(ky/32+kx/64)/35;
    check(std::abs((after[0]-before[0])-delta)<=8*std::numeric_limits<long double>::epsilon(),
        "direct XY affine P2 gradient disagrees with independent Green identity");
    LAMP::Interface::Detail::SmoothDirectXYWork work;
    const auto result=LAMP::Interface::Detail::constructDirectSmoothArcXY(arc,cell,1e-12L,&work);
    directSuccess(arc,cell,result,work);
    check(cell.inventory==stocks&&arc.controlPoints[2].x==original.controlPoints[2].x
        &&arc.controlPoints[2].y==original.controlPoints[2].y,"direct XY modified immutable source inputs");
    // Independent stocks need not sum EXACTLY to capacity. Here the true
    // minority target is outside the major-phase band, while their joint
    // feasible interval is nonempty. This must positively construct a curve,
    // rather than treating a representability rejection as successful repair.
    auto compatible=cell;compatible.inventory={{.6000000000007L,.4L}};
    auto horizontal=arc;horizontal.controlPoints={Point{0,.4},Point{.25,.4},Point{.5,.4},Point{.75,.4},Point{1,.4}};
    const auto compatibleStocks=compatible.inventory;
    LAMP::Geometry::Detail::FiniteVolume::Dyadic capacityDefect(compatible.inventory[0]);
    capacityDefect.add(compatible.inventory[1]);capacityDefect.add(1,-1);
    check(capacityDefect.compare(LAMP::Geometry::Detail::FiniteVolume::Dyadic{})!=0,
        "dual-stock counterexample unexpectedly has exact capacity closure");
    const auto initial=LAMP::Interface::Detail::analyticSmoothArcPhasesXY(horizontal,compatible);
    stripClose(initial[1],compatible.inventory[1],"dual-stock counterexample lost its true-minor target");
    check(std::abs(initial[0]-compatible.inventory[0])>1e-12L*compatible.inventory[0],
        "true-minor target no longer violates the independent major band");
    stripClose(.3999999999998L,compatible.inventory[1],"independent compatible minority witness is outside its band");
    stripClose(.6000000000002L,compatible.inventory[0],"independent compatible major witness is outside its band");
    const auto joint=LAMP::Interface::Detail::constructDirectSmoothArcXY(horizontal,compatible,1e-12L,&work);
    directSuccess(horizontal,compatible,joint,work);
    check(joint.measured[1]<compatible.inventory[1],"joint-band constructor retained the infeasible exact-minority target");
    const auto reserve=8*std::numeric_limits<long double>::epsilon()*std::min(compatible.inventory[0],compatible.inventory[1]);
    for(unsigned phase=0;phase<2;++phase)check(std::abs(joint.measured[phase]-compatible.inventory[phase])
        <1e-12L*compatible.inventory[phase]-reserve,"joint-band candidate has no positive native acceptance margin");
    check(compatible.inventory==compatibleStocks,"joint-band selection modified the independent stocks");
    auto impossible=cell;impossible.inventory={{.75,.75}};
    rejects([&]{LAMP::Interface::Detail::constructDirectSmoothArcXY(arc,impossible,1e-12L,&work);},
        VofSmoothGeometryFailure::InfeasibleVolume,"incompatible independent stock bands were accepted");
    check(work.candidateConstructions==0&&work.fitPhaseMeasurements==0&&work.volumeNewtonSteps==0,
        "incompatible stock bands generated an actual curve or fitting work");
    // A horizontal family has no x area gradient. Its y gradient is 57/280;
    // the dyadic target below requires a non-dyadic y, so exact world-native
    // representability fails at this explicitly strict test tolerance.
    auto quantized=arc;quantized.controlPoints={Point{0,.5},Point{15.L/64,.5},Point{.5,.5},Point{49.L/64,.5},Point{1,.5}};
    auto narrow=cell;narrow.inventory={{.5L-std::ldexp(1.L,-40),.5L+std::ldexp(1.L,-40)}};
    rejects([&]{LAMP::Interface::Detail::constructDirectSmoothArcXY(quantized,narrow,1e-25L,&work);},
        VofSmoothGeometryFailure::InfeasibleVolume,"unrepresentable world-native own-area target was silently accepted");
    check(work.candidateConstructions==1&&work.geometryCertifications==1&&work.phaseQuadratureBuilds==0,
        "unrepresentable world target used repeated candidates or quadrature fitting");
    std::cout<<"DIRECT_XY_ALGEBRA,PASS affine_area=1 fixed_domain=1 one_candidate=1 fit_measures=0 newton=0 original_gates=1 unrepresentable_rejected=1\n";
}
void directHistoricalAttempt(const char* name,const VofSmoothArc& arc,const GeometricVofPlicCell& cell)
{
    const auto original=arc;const auto stocks=cell.inventory;LAMP::Interface::Detail::SmoothDirectXYWork work;
    bool accepted=false;
    try{const auto result=LAMP::Interface::Detail::constructDirectSmoothArcXY(arc,cell,1e-12L,&work);
        directSuccess(arc,cell,result,work);accepted=true;}
    catch(const VofSmoothGeometryError& error){check(error.reason()==VofSmoothGeometryFailure::InfeasibleVolume,
        "historical direct XY attempt produced an undeclared failure");
        check(work.candidateConstructions<=1&&work.geometryCertifications<=1&&work.phaseQuadratureBuilds<=1,
            "historical direct XY rejection retried a candidate, certificate or quadrature");}
    check(work.fitPhaseMeasurements==0&&work.volumeNewtonSteps==0,"historical direct XY attempt retained old volume fitting");
    check(cell.inventory==stocks,"historical direct XY attempt changed independent stocks");
    for(unsigned i=0;i<5;++i)check(arc.controlPoints[i].x==original.controlPoints[i].x&&arc.controlPoints[i].y==original.controlPoints[i].y,
        "historical direct XY attempt changed the immutable failure fixture");
    std::cout<<name<<" direct_world_status="<<(accepted?"qualified":"typed_rejected")
        <<" candidates="<<work.candidateConstructions<<" fit_measures=0 newton=0\n";
}
void boundedCenter()
{
    // The original full04 exterior fixture remains a negative geometric
    // witness. Its old reflection recovery is retired, not silently deleted.
    VofSmoothArc arc;arc.nativeCell=688;arc.graphAxis=0;arc.nodes={{17,23}};
    arc.controlPoints={Point{.53125,.68746348256573353},Point{.53122421508755702,.68747001738934288},
        Point{.53120003370030622,.68751986392393072},Point{.53117585105893783,.6874934800397563},Point{.5311500623841422,.6875}};
    GeometricVofPlicCell cell;cell.rectangle={{.5,.65625},{.53125,.6875}};cell.inventory={{.00097656143059681919,1.0694031807979259e-9}};
    unsigned depth=0;check(!LAMP::Interface::Detail::certifySmoothArc(arc,cell.rectangle,depth),"full04 true exterior curve no longer rejects");
    check(arc.sample(.705L).y>cell.rectangle.upper.y,"full04 negative fixture no longer has a true exterior point");
    const auto prior=LAMP::Interface::Detail::analyticSmoothArcPhasesXY(arc,cell);
    for(unsigned p=0;p<2;++p)stripClose(prior[p],cell.inventory[p],"full04 failure was misclassified as an area defect");
    directHistoricalAttempt("DIRECT_FULL04",arc,cell);
    directAlgebra();
}
void retainedNativeVolume()
{
    // Preserve source22's actual two world-native values as failure evidence;
    // the new constructor consumes neither a history nor adjacent candidates.
    VofSmoothArc arc;arc.nativeCell=688;arc.graphAxis=0;arc.nodes={{17,23}};
    arc.controlPoints={Point{.53125,.68746323751933613},Point{.53122397520595555,.68746983324052824},
        Point{.53119955155227361,.68751886518362093},Point{.53117512662059529,.68749341942511688},Point{.53114909799256171,.6875}};
    GeometricVofPlicCell cell;cell.rectangle={{.5,.65625},{.53125,.6875}};cell.inventory={{.00097656139047658643,1.1095234135522208e-9}};
    const auto bad=LAMP::Interface::Detail::measureSmoothArcPhases(arc,cell,false);
    auto good=arc;good.controlPoints[2].y=.68751886518362082;
    const auto prior=LAMP::Interface::Detail::measureSmoothArcPhases(good,cell,false);
    check(std::abs(bad[1]-cell.inventory[1])>1e-12L*cell.inventory[1],"source22 worse native point lost the original failure");
    stripClose(prior[1],cell.inventory[1],"source22 qualified prior native point no longer qualifies");
    directHistoricalAttempt("DIRECT_SOURCE22",arc,cell);
}
void centerNativeVolume()
{
    VofSmoothArc source37;source37.nativeCell=688;source37.graphAxis=0;source37.nodes={{17,23}};
    source37.controlPoints={Point{.53125,.68746300057850207},Point{.53122380286799598,.68746963979889952},
        Point{.5311992164054451,.68752518990972533},Point{.53117462864764076,.68749337613097428},Point{.53114842762987713,.6875}};
    GeometricVofPlicCell cell;cell.rectangle={{.5,.65625},{.53125,.6875}};cell.inventory={{.00097656150078457031,9.992154295198383e-10}};
    unsigned depth=0;check(!LAMP::Interface::Detail::certifySmoothArc(source37,cell.rectangle,depth),"source37 true exterior curve passed");
    auto oldSelected=source37;oldSelected.controlPoints[2]={.53122482737886667,.68751501389748315};
    check(LAMP::Interface::Detail::certifySmoothArc(oldSelected,cell.rectangle,depth),"source37 historical selected point lost geometry");
    const auto failed=LAMP::Interface::Detail::measureSmoothArcPhases(oldSelected,cell,false);
    check(std::abs(failed[1]-cell.inventory[1])>1e-12L*cell.inventory[1],"source37 original own-phase failure was erased");
    directHistoricalAttempt("DIRECT_SOURCE37",source37,cell);
    auto source40=source37;source40.controlPoints={Point{.53125,.68746379478805242},Point{.53122456816916253,.68747023983728972},
        Point{.53120075243988107,.68752566082546929},Point{.53117693549129563,.6874935694019354},Point{.53115150000254563,.6875}};
    cell.inventory={{.00097656157127389465,9.2872610510478794e-10}};
    check(!LAMP::Interface::Detail::certifySmoothArc(source40,cell.rectangle,depth),"source40 original exterior curve passed");
    directHistoricalAttempt("DIRECT_SOURCE40",source40,cell);
}

std::array<long double,2> independentGraphStrip(const VofSmoothArc& arc,
    const LAMP::Geometry::PlanarCellCut::Rectangle& box,bool upper,long double volume)
{
    // Independent reference: ordinary Bezier bisection plus degree-exact
    // positive GL4 of the original curve. Only moderate strip widths use
    // world coordinates here; the production tiny-strip test below does not.
    check(arc.graphAxis==0,"reference expects an x graph");
    const auto h=static_cast<long double>(box.upper.y)-box.lower.y;
    const auto x=upper?static_cast<long double>(box.upper.x)-volume/h:
        static_cast<long double>(box.lower.x)+volume/h;
    const bool increasing=arc.sample(1).x>arc.sample(0).x;long double low=0,high=1;
    for(unsigned k=0;k<90;++k){const auto middle=(low+high)/2;if(middle==low||middle==high)break;
        if((arc.sample(middle).x<x)==increasing)low=middle;else high=middle;}
    const auto cut=(low+high)/2;const bool first=(upper?!increasing:increasing);
    const auto begin=first?0.L:cut,end=first?cut:1.L;
    constexpr std::array<long double,4> node{{.069431844202973712388L,.33000947820757186760L,
        .66999052179242813240L,.93056815579702628761L}};
    constexpr std::array<long double,4> weight{{.17392742256872692869L,.32607257743127307131L,
        .32607257743127307131L,.17392742256872692869L}};
    const bool liquidLower=!increasing;std::array<long double,2> result{};
    for(unsigned k=0;k<4;++k){const auto parameter=(1-node[k])*begin+node[k]*end;
        const auto p=arc.sample(parameter),d=arc.derivative(parameter);
        const auto jacobian=(end-begin)*weight[k]*std::abs(d.x);
        result[liquidLower?0:1]+=jacobian*(p.y-box.lower.y);
        result[liquidLower?1:0]+=jacobian*(box.upper.y-p.y);}
    return result;
}
void nativeStrips()
{
    auto mesh=meshFor();VofSmoothGeometryConfig config;config.allowBoundaryEndpoints=true;
    std::size_t collapsedWorldEndpoints=0;
    for(const bool curved:{false,true}){
        const auto source=GeometricVofInitialization::prepare(*mesh,false,[=](Vector2D p){
            return p.y-(curved?.421+.031*std::sin(2*std::acos(-1.)*p.x):.423);});
        const auto geometry=VofSmoothGeometry::prepare(*mesh,source,config);
        const auto original=source.plicCells();
        for(const auto& cell:geometry.cells()){
            const auto capacity=LAMP::Geometry::PlanarCellCut::capacity(cell.rectangle,false);
            const auto q=.37L*capacity;
            for(unsigned axis=0;axis<2;++axis){
                const auto lower=geometry.measureNativeStrip(*mesh,source,cell.nativeCell,axis,false,q);
                const auto upper=geometry.measureNativeStrip(*mesh,source,cell.nativeCell,axis,true,capacity-q);
                for(unsigned p=0;p<2;++p){
                    stripClose(lower.phase[p]+upper.phase[p],cell.inventory[p],"opposite smooth strips lost an independent phase");
                    check(lower.rawPhase[p]>=0&&upper.rawPhase[p]>=0,"negative independently measured strip phase");
                    check(lower.rawPhase[p]==0?lower.absoluteError[p]==0:lower.absoluteError[p]<=1e-12L*lower.rawPhase[p],
                        "strip has no own-phase integral certificate");}
                const auto empty=geometry.measureNativeStrip(*mesh,source,cell.nativeCell,axis,false,0);
                check(empty.phase==std::array<long double,2>{{0,0}}&&empty.rawPhase==empty.phase,"zero strip is not exact");
                const auto full=geometry.measureNativeStrip(*mesh,source,cell.nativeCell,axis,true,capacity);
                for(unsigned p=0;p<2;++p)stripClose(full.phase[p],cell.inventory[p],"full smooth strip differs from its actual own-phase stock");
            }
            if(cell.arc<0)continue;
            const auto& arc=geometry.arcs()[static_cast<std::size_t>(cell.arc)];
            check(arc.graphAxis==0,"mild graph fixture unexpectedly changed graph axis");
            for(const bool upper:{false,true}){
                const auto actual=geometry.measureNativeStrip(*mesh,source,cell.nativeCell,0,upper,q);
                const auto expected=independentGraphStrip(arc,cell.rectangle,upper,q);
                for(unsigned p=0;p<2;++p)stripClose(actual.phase[p],expected[p],"curved strip differs from independent Bezier/GL4 reference");
                if(!curved){const auto h=static_cast<long double>(cell.rectangle.upper.y)-cell.rectangle.lower.y;
                    stripClose(actual.phase[0],q*((.423L-cell.rectangle.lower.y)/h),"flat strip differs from the independent line area");
                    stripClose(actual.phase[1],q*((cell.rectangle.upper.y-.423L)/h),"flat strip lost independent gas area");}
                const auto tiny=capacity*1e-20L;
                const auto micro=geometry.measureNativeStrip(*mesh,source,cell.nativeCell,0,upper,tiny);
                const auto face=upper?cell.rectangle.upper.x:cell.rectangle.lower.x;
                const auto depth=tiny/(static_cast<long double>(cell.rectangle.upper.y)-cell.rectangle.lower.y);
                if(face!=0){check(static_cast<double>(face+(upper?-depth:depth))==face,"microstrip does not exercise a collapsed world endpoint");
                    ++collapsedWorldEndpoints;}
                const auto endpoint=std::abs(arc.sample(0).x-face)<std::abs(arc.sample(1).x-face)?arc.sample(0):arc.sample(1);
                const bool liquidLower=arc.sample(1).x<arc.sample(0).x;
                const auto h=static_cast<long double>(cell.rectangle.upper.y)-cell.rectangle.lower.y;
                const auto lowerHeight=endpoint.y-cell.rectangle.lower.y;
                const std::array<long double,2> expectedMicro=liquidLower?
                    std::array<long double,2>{{tiny*(lowerHeight/h),tiny*((h-lowerHeight)/h)}}:
                    std::array<long double,2>{{tiny*((h-lowerHeight)/h),tiny*(lowerHeight/h)}};
                for(unsigned p=0;p<2;++p){check(micro.phase[p]>0,"positive microscopic smooth strip phase was erased");
                    stripClose(micro.phase[p],expectedMicro[p],"face-anchored microscopic strip lost its actual phase fraction");}
                stripClose(micro.phase[0]+micro.phase[1],tiny,"microscopic shared FT was not retained");
            }
        }
        const auto foreign=GeometricVofInitialization::prepare(*mesh,false,[](Vector2D p){return p.y-.423;});
        rejects([&]{geometry.measureNativeStrip(*mesh,foreign,0,0,false,0);},VofSmoothGeometryFailure::ForeignSource,
            "strip accepted a different actual source owner");
        for(const auto bad:{-1.L,std::numeric_limits<long double>::infinity(),std::numeric_limits<long double>::quiet_NaN(),1.L})
            rejects([&]{geometry.measureNativeStrip(*mesh,source,0,0,false,bad);},VofSmoothGeometryFailure::InvalidCurve,
                "strip accepted an invalid volume range");
        rejects([&]{geometry.measureNativeStrip(*mesh,source,0,2,false,0);},VofSmoothGeometryFailure::InvalidCurve,"strip accepted a foreign normal axis");
        rejects([&]{geometry.measureNativeStrip(*mesh,source,-1,0,false,0);},VofSmoothGeometryFailure::ForeignSource,"strip accepted a foreign cell");
        for(std::size_t ci=0;ci<original.size();++ci)check(source.plicCells()[ci].inventory==original[ci].inventory,
            "strip query changed actual independent stocks");
    }
    check(collapsedWorldEndpoints>0,"no genuine collapsed world endpoint was tested");
    const auto radial=GeometricVofInitialization::prepare(*mesh,true,[](Vector2D p){return p.y-.423;});
    const auto radialGeometry=VofSmoothGeometry::prepare(*mesh,radial,config);
    rejects([&]{radialGeometry.measureNativeStrip(*mesh,radial,0,0,false,0);},VofSmoothGeometryFailure::UnsupportedMesh,
        "XY strip API silently accepted RZ");
    std::cout<<"smooth_native_strip_contract PASS assertions="<<assertions<<'\n';
}
void requestedStripVolumeResidual()
{
    namespace FV=LAMP::Geometry::Detail::FiniteVolume;
    const double h=.1,d=.03,flux=.003;
    auto residual=FV::Dyadic(flux);residual.add(FV::Dyadic::product({h,d}),-1);
    check(residual.compare(FV::Dyadic{})>0&&flux/h==d,"exact FT/h counterexample disappeared");
    stripClose(residual.rounded(),6.938893903907235e-21L,"incorrect independent exact strip remainder");
    DictionaryMesh input;input.Nx=input.Ny=4;input.minCoord={0,0};input.maxCoord={.4,.4};
    input.maxCellNum=256;input.maxInitialIterations=0;QuadtreeMesh mesh(input);
    const auto source=GeometricVofInitialization::prepare(mesh,false,[](Vector2D p){return p.y-.07;});
    VofSmoothGeometryConfig config;config.allowBoundaryEndpoints=true;
    const auto geometry=VofSmoothGeometry::prepare(mesh,source,config);const auto& cell=geometry.cells().front();
    check(cell.arc>=0&&cell.rectangle.upper.x-cell.rectangle.lower.x==h,"non-dyadic strip fixture has wrong cell");
    const auto& arc=geometry.arcs().at(static_cast<std::size_t>(cell.arc));
    bool positiveWitness=false;
    for(const auto parameter:{0.L,1.L}){const auto point=arc.sample(parameter);
        FV::Dyadic distance(cell.rectangle.upper.y);distance.add(point.y,-1);
        if(distance.compare(FV::Dyadic(d))==0)positiveWitness=true;}
    check(positiveWitness,"non-dyadic strip fixture lost its exact positive endpoint witness");
    // Continuity at that actual endpoint proves a positive liquid set. The
    // previous RN depth implementation returned exact zero. A backend unable
    // to certify its tiny integral must issue a typed failure instead.
    try{const auto measure=geometry.measureNativeStrip(mesh,source,cell.nativeCell,1,true,flux);
        check(measure.phase[0]>0&&measure.rawPhase[0]>0,"RN depth erased an exactly positive far phase");
        check(measure.absoluteError[0]<=1e-12L*measure.rawPhase[0],"positive residual phase lacks its own error certificate");
        stripClose(measure.phase[0]+measure.phase[1],flux,"exact requested strip FT was changed");
        bool flat=true;for(const auto control:arc.controlPoints)flat=flat&&control.y==.07;
        if(flat)stripClose(measure.rawPhase[0],residual.rounded(),"flat residual phase differs from independent exact FT-h*d");
    }catch(const VofSmoothGeometryError& error){check(error.reason()==VofSmoothGeometryFailure::VolumeResidual,
        "non-dyadic positive strip failed for an unrelated reason");
        std::cout<<"smooth_strip_exact_FT typed_refusal="<<error.what()<<'\n';}
}
void examine(const QuadtreeMesh& mesh,const GeometricVofPublication& source,bool curved)
{
    const auto stocks=source.plicCells();
    VofSmoothGeometryConfig config;config.allowBoundaryEndpoints=true;
    const auto geometry=VofSmoothGeometry::prepare(mesh,source,config);
    geometry.validateFor(mesh,source);
    check(geometry.sourceIdentity()==source.identity(),"source identity changed");
    check(geometry.sourceTime().step==source.timeLayer().step
        &&geometry.sourceTime().begin==source.timeLayer().begin
        &&geometry.sourceTime().end==source.timeLayer().end
        &&geometry.sourceTime().position==source.timeLayer().position,"full source clock changed");
    const auto copy=geometry;
    check(copy.sameInstance(geometry),"copy lost actual smooth owner");
    const auto other=VofSmoothGeometry::prepare(mesh,source,config);
    check(!other.sameInstance(geometry),"separate derived geometry forged same owner");
    check(!geometry.arcs().empty(),"missing actual B1 curve");
    check(geometry.report().maximumOwnPhaseRelativeResidual<=1e-12L,"own-phase volume certificate failed");
    const auto& direct=geometry.report();
    if(!source.domain().axisymmetric())check(direct.xyDirectCandidates==direct.mixedCells
        &&direct.xyGeometryCertifications==direct.mixedCells&&direct.xyPhaseQuadratureBuilds==direct.mixedCells,
        "actual XY prepare did not consume the one-candidate direct constructor");
    else check(direct.xyDirectCandidates==0&&direct.xyGeometryCertifications==0&&direct.xyPhaseQuadratureBuilds==0,
        "RZ entered the XY direct construction policy");
    check(direct.xyFitPhaseMeasurements==0&&direct.xyVolumeNewtonSteps==0,"actual XY prepare retained volume fitting work");
    std::vector<std::vector<std::pair<std::size_t,unsigned>>> incidence(geometry.nodes().size());
    bool nonStraight=false,nonzeroCorrection=false;
    long double total=0,firstX=0;
    for(const auto& cell:geometry.cells()){
        for(unsigned p=0;p<2;++p){
            long double sum=0;
            for(const auto& q:cell.phaseQuadrature[p]){
                check(q.weight>0&&std::isfinite(q.weight),"nonpositive phase point");
                check(q.point.x>=cell.rectangle.lower.x&&q.point.x<=cell.rectangle.upper.x
                    &&q.point.y>=cell.rectangle.lower.y&&q.point.y<=cell.rectangle.upper.y,"phase point outside cell");
                sum+=q.weight;total+=q.weight;firstX+=q.weight*q.point.x;
                check(geometry.phaseAt(cell.nativeCell,q.point)==p,"actual smooth phase query disagrees with its positive phase quadrature");
            }
            check(cell.inventory[p]==stocks[static_cast<std::size_t>(cell.nativeCell)].inventory[p],"authoritative inventory changed");
            check(cell.inventory[p]==0?sum==0:std::abs(sum-cell.inventory[p])<=1e-12L*cell.inventory[p],"quadrature does not match own phase");
        }
        if(cell.arc<0)continue;
        const auto& arc=geometry.arcs()[static_cast<std::size_t>(cell.arc)];
        nonzeroCorrection=nonzeroCorrection||std::abs(arc.volumeCorrection)>1e-9L;
        for(unsigned end=0;end<2;++end){
            incidence[arc.nodes[end]].push_back({static_cast<std::size_t>(cell.arc),end});
            check(distance(arc.sample(end),geometry.nodes()[arc.nodes[end]].point)==0,"curve endpoint differs from actual shared node");
            const auto d=arc.derivative(end),t=geometry.nodes()[arc.nodes[end]].tangent;
            check(std::abs(d.x*t.y-d.y*t.x)<=1e-12L*std::hypot(d.x,d.y),"shared tangent direction differs");
            check(d.x*t.x+d.y*t.y>0,"shared tangent is reversed");
        }
        const auto chord=Point{arc.sample(1).x-arc.sample(0).x,arc.sample(1).y-arc.sample(0).y};
        const auto d=arc.derivative(.5L);
        nonStraight=nonStraight||std::abs(chord.x*d.y-chord.y*d.x)>1e-8L;
        const auto quarter=arc.sample(.25L),start=arc.sample(0);
        nonStraight=nonStraight||std::abs(chord.x*(quarter.y-start.y)-chord.y*(quarter.x-start.x))>1e-8L;
        const auto normal=arc.normal(.4L),tangent=arc.derivative(.4L);
        check(std::abs(normal.x*tangent.x+normal.y*tangent.y)<1e-12L,"normal/tangent inconsistency");
        check(std::abs(std::hypot(normal.x,normal.y)-1)<1e-12L,"nonunit normal");
        rejects([&]{geometry.phaseAt(cell.nativeCell,arc.sample(.5L));},VofSmoothGeometryFailure::InvalidCurve,
            "point on actual curve silently acquired one phase");
        if(!curved)for(const auto t:{.17L,.5L,.83L}){
            check(std::abs(arc.sample(t).y-.423L)<1e-12L,"coordinate-flat inventory did not remain a flat curve");
            const auto n=arc.normal(t);
            check(std::abs(n.x)<1e-12L&&std::abs(n.y-1)<1e-12L,"coordinate-flat inventory acquired a transverse normal");
        }
    }
    for(const auto& list:incidence)check(list.size()==1||list.size()==2,"branch or unowned node");
    const bool rz=source.domain().axisymmetric();const auto pi=std::acos(-1.L);
    check(std::abs(total-(rz?pi:1))<1e-12L,"actual XY/RZ domain measure mismatch");
    check(std::abs(firstX-(rz?2*pi/3:.5L))<1e-12L,"phase quadrature radial first moment mismatch");
    if(curved){check(nonStraight,"curved inventory silently became straight PLIC");check(nonzeroCorrection,"curved inventory has no actual quartic correction");}
    rejects([&]{VofSmoothGeometry::prepare(mesh,source);},VofSmoothGeometryFailure::BoundaryEndpoint,"unqualified boundary silently accepted");
    check(source.plicCells()[0].inventory==stocks[0].inventory,"rejection changed accepted VOF stock");
    std::cout<<"geometry rz="<<rz<<" curved="<<curved<<" arcs="<<geometry.arcs().size()
        <<" certificate_depth="<<geometry.report().maximumCertificateSubdivisionDepth
        <<" own_phase_residual="<<static_cast<double>(geometry.report().maximumOwnPhaseRelativeResidual)<<'\n';
}
void circle()
{
    DictionaryMesh input;input.Nx=input.Ny=16;input.maxCellNum=2048;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);
    const auto source=GeometricVofInitialization::prepare(mesh,false,[](Vector2D p){return std::hypot(p.x-.487,p.y-.513)-.18;});
    const auto geometry=VofSmoothGeometry::prepare(mesh,source);
    std::vector<unsigned> incidence(geometry.nodes().size());long double radialError=0;
    for(const auto& arc:geometry.arcs()){
        ++incidence[arc.nodes[0]];++incidence[arc.nodes[1]];
        for(const auto t:{0.L,.25L,.5L,.75L,1.L}){const auto p=arc.sample(t);
            radialError=std::max(radialError,std::abs(std::hypot(p.x-.487L,p.y-.513L)-.18L));}
    }
    for(const auto degree:incidence)check(degree==2,"closed circle acquired a missing or branching actual node");
    check(geometry.report().boundaryEndpoints==0,"closed circle acquired exterior endpoints");
    check(geometry.report().maximumOwnPhaseRelativeResidual<=1e-12L,"closed circle own-phase volume certificate failed");
    for(const auto& cell:geometry.cells())for(unsigned phase=0;phase<2;++phase)
        for(const auto& q:cell.phaseQuadrature[phase])check(geometry.phaseAt(cell.nativeCell,q.point)==phase,
            "closed actual curve phase classification differs from its own integration domain");
    std::cout<<"circle arcs="<<geometry.arcs().size()<<" own_phase_residual="
        <<static_cast<double>(geometry.report().maximumOwnPhaseRelativeResidual)
        <<" independent_sample_radial_error="<<static_cast<double>(radialError)<<'\n';
}
void cellReceipt(const QuadtreeMesh& mesh,const GeometricVofPublication& publication,int native,
    const GeometricVofTransport* transport)
{
    const auto& cell=publication.plicCells().at(static_cast<std::size_t>(native));
    const auto capacity=LAMP::Geometry::PlanarCellCut::capacity(cell.rectangle,false);
    std::cout<<"  cell="<<native<<" box="<<cell.rectangle.lower.x<<','<<cell.rectangle.lower.y<<':'
        <<cell.rectangle.upper.x<<','<<cell.rectangle.upper.y<<" stock="<<cell.inventory[0]<<','<<cell.inventory[1]
        <<" measured="<<cell.measured[0]<<','<<cell.measured[1]<<" capacity="<<capacity
        <<" minor_relative_capacity="<<std::min(cell.inventory[0],cell.inventory[1])/capacity
        <<" resolved="<<cell.resolved<<" normal="<<cell.normal.x<<','<<cell.normal.y<<" offset="<<cell.offset<<'\n';
    for(const auto f:publication.domain().faces(mesh))if(f.negativeCell==native||f.positiveCell==native){
        const int other=f.negativeCell==native?f.positiveCell:f.negativeCell;
        std::cout<<"    face="<<f.nativeRow<<" cells="<<f.negativeCell<<','<<f.positiveCell
            <<" horizontal="<<f.horizontal<<" center="<<f.center.x<<','<<f.center.y;
        if(other>=0){const auto& adjacent=publication.plicCells().at(static_cast<std::size_t>(other));
            std::cout<<" adjacent_stock="<<adjacent.inventory[0]<<','<<adjacent.inventory[1]
                <<" adjacent_measured="<<adjacent.measured[0]<<','<<adjacent.measured[1]
                <<" adjacent_resolved="<<adjacent.resolved;}
        if(transport){const auto& faces=transport->faces();
            const auto found=std::find_if(faces.begin(),faces.end(),[&](const auto& flux){return flux.nativeRow==f.nativeRow;});
            check(found!=faces.end(),"movement receipt lacks its actual native flux");
            for(const auto half:{true,false}){const auto& flux=half?found->half:found->full;
                std::cout<<(half?" half=(":" full=(")<<"total="<<flux.total<<",rawLiquid="<<flux.rawLiquid
                    <<",liquid="<<flux.liquid<<",gas="<<flux.gas<<",limiter="<<flux.limiter<<')';}}
        std::cout<<'\n';
    }
}
void geometryReceipt(const QuadtreeMesh& mesh,const GeometricVofPublication& publication,
    const char* initializer,const char* stage,unsigned n,const GeometricVofTransport* transport=nullptr)
{
    const auto clock=publication.timeLayer();const std::vector<int> watch=n==16?std::vector<int>{133,165}:std::vector<int>{688};
    std::cout<<"B1_G1_RECEIPT initializer="<<initializer<<" n="<<n<<" stage="<<stage
        <<" source_identity="<<publication.identity().value()<<" clock=("<<clock.step<<','<<clock.begin<<','
        <<clock.end<<','<<static_cast<int>(clock.position)<<")\n";
    for(const auto native:watch)cellReceipt(mesh,publication,native,transport);
    try{const auto geometry=VofSmoothGeometry::prepare(mesh,publication);
        std::cout<<"  geometry=accepted arcs="<<geometry.arcs().size()<<" own_phase_residual="
            <<geometry.report().maximumOwnPhaseRelativeResidual<<'\n';}
    catch(const VofSmoothGeometryError& error){
        std::cout<<"  geometry=typed_rejected reason="<<static_cast<int>(error.reason())<<" native="<<error.nativeCell()
            <<" detail="<<error.what()<<'\n';
        if(error.nativeCell()>=0&&std::find(watch.begin(),watch.end(),error.nativeCell())==watch.end())
            cellReceipt(mesh,publication,error.nativeCell(),transport);
    }
}
long double receiptTaper(long double x)
{
    if(x<=0||x>=1)return 0;if(x>=.2L&&x<=.8L)return 1;
    const auto z=x<.2L?x/.2L:(1-x)/.2L;return z*z*(3-2*z);
}
void movementReceipts()
{
    const auto flags=std::cout.flags();const auto precision=std::cout.precision();
    std::cout<<std::scientific<<std::setprecision(std::numeric_limits<long double>::max_digits10);
    for(const unsigned n:{16U,32U}){
        DictionaryMesh input;input.Nx=input.Ny=n;input.maxCellNum=4*n*n;input.maxInitialIterations=0;
        QuadtreeMesh mesh(input);
        const auto squared=GeometricVofInitialization::prepare(mesh,false,[](Vector2D p){
            const auto x=p.x-.487,y=p.y-.513;return x*x+y*y-.18*.18;});
        geometryReceipt(mesh,squared,"squared_circle_initial_only","Initial",n);
        const auto source=GeometricVofInitialization::prepare(mesh,false,[](Vector2D p){return std::hypot(p.x-.487,p.y-.513)-.18;});
        const auto acceptedCells=source.plicCells();const auto identity=source.identity();const auto clock=source.timeLayer();
        geometryReceipt(mesh,source,"signed_distance_circle","Initial",n);
        const auto psi=[](long double x,long double y){return .04L*y*receiptTaper(x)*receiptTaper(y);};
        std::vector<double> velocity(static_cast<std::size_t>(source.domain().nativeFaceCount()));
        for(const auto face:source.domain().faces(mesh)){const long double x=face.center.x,y=face.center.y,h=face.area/2;
            velocity[static_cast<std::size_t>(face.nativeRow)]=static_cast<double>(face.horizontal?
                -(psi(x+h,y)-psi(x-h,y))/face.area:(psi(x,y+h)-psi(x,y-h))/face.area);}
        const GeometricVofTimeLayer target{1,0,.025,GeometricVofTimeLayer::Position::End};
        const auto slab=GeometricVofAdvectorSlab::prepare(mesh,source,target,velocity,velocity);
        const auto transport=GeometricVofTransport::prepareLocalPlicSweep(mesh,source,slab,1,1);
        const auto& report=transport.report();std::cout<<"  transport_receipt n="<<n<<" dt="<<target.end-target.begin
            <<" local_unresolved_donors="<<report.localUnresolvedDonors<<" limited_faces="<<report.limitedFaces
            <<" maximum_flux_correction="<<report.maximumFluxCorrection<<'\n';
        geometryReceipt(mesh,transport.midpoint(),"signed_distance_circle","Midpoint",n,&transport);
        geometryReceipt(mesh,transport.endpoint(),"signed_distance_circle","End",n,&transport);
        const auto now=source.timeLayer();check(source.identity()==identity&&now.step==clock.step&&now.begin==clock.begin
            &&now.end==clock.end&&now.position==clock.position,"movement receipt changed accepted source identity/full clock");
        for(std::size_t ci=0;ci<acceptedCells.size();++ci)
            check(source.plicCells()[ci].inventory==acceptedCells[ci].inventory,"movement receipt changed accepted independent stocks");
    }
    std::cout.flags(flags);std::cout.precision(precision);
    std::cout<<"B1 movement receipts emitted; acceptance_status=not_a_motion_contract\n";
}
}
int main(int argc,char** argv)
{
    try {
        if(argc==2&&std::string(argv[1])=="--movement-receipts"){movementReceipts();return 0;}
        if(argc==2&&std::string(argv[1])=="--strip-only"){nativeStrips();requestedStripVolumeResidual();return 0;}
        if(argc==2&&std::string(argv[1])=="--phase-measure-only"){anchoredPhaseMeasure();return 0;}
        if(argc==2&&std::string(argv[1])=="--direct-xy-only"){boundedCenter();retainedNativeVolume();centerNativeVolume();return 0;}
        if(argc==2&&std::string(argv[1])=="--bounded-center-only"){boundedCenter();return 0;}
        if(argc==2&&std::string(argv[1])=="--native-volume-only"){retainedNativeVolume();return 0;}
        if(argc==2&&std::string(argv[1])=="--center-native-volume-only"){centerNativeVolume();return 0;}
        check(argc==1,"unknown geometry contract argument");
        auto mesh=meshFor();
        for(const bool rz:{false,true}){
            const auto initial=GeometricVofInitialization::prepare(*mesh,rz,[](Vector2D p){return p.y-.423;});
            std::vector<std::array<long double,2>> frozenStock;
            for(const auto& cell:initial.plicCells())frozenStock.push_back(cell.inventory);
            const GeometricVofTimeLayer clock{3,.2,.25,GeometricVofTimeLayer::Position::Midpoint};
            const auto plane=GeometricVofPublication::prepareLocalPlicInventory(*mesh,initial.domain(),clock,frozenStock,{});
            examine(*mesh,plane,false);
            // Only initialization queries the independently specified graph.
            // The production reconstruction receives the resulting VOF owner.
            const auto curved=GeometricVofInitialization::prepare(*mesh,rz,[](Vector2D p){
                return p.y-(.421+.031*std::sin(2*std::acos(-1.)*p.x));});
            examine(*mesh,curved,true);
            VofSmoothGeometryConfig config;config.allowBoundaryEndpoints=true;
            const auto geometry=VofSmoothGeometry::prepare(*mesh,plane,config);
            rejects([&]{geometry.validateFor(*mesh,curved);},VofSmoothGeometryFailure::ForeignSource,"foreign source publication accepted");
        }
        // A real independently retained microscopic phase has no resolved
        // source plane; it must reject, never be erased or converted to PLIC.
        FaceSegmentGraph graph;graph.build(*mesh);auto domain=NativeLeafDomain::build(*mesh,graph,false);
        std::vector<std::array<long double,2>> stock(static_cast<std::size_t>(mesh->cellIndexNum()));
        for(std::size_t i=0;i<stock.size();++i){const auto row=domain.cellRowByIndex()[i];stock[i]={0,domain.cellVolumes()[static_cast<std::size_t>(row)]};}
        const auto capacity=stock[0][1];stock[0]={capacity*1e-18L,capacity};
        const auto unresolved=GeometricVofPublication::prepareLocalPlicInventory(*mesh,domain,{},stock,{});
        rejects([&]{VofSmoothGeometry::prepare(*mesh,unresolved);},VofSmoothGeometryFailure::UnresolvedCell,"unresolved minority phase was silently removed");
        check(unresolved.plicCells()[0].inventory[0]==stock[0][0],"unresolved failure modified inventory");
        auto mixedMesh=meshFor();DictionaryMesh::RefineBox refine;refine.level=1;refine.min={.5,0};refine.max={1,1};
        mixedMesh->applyRefineBoxes({refine});
        const auto mixedSource=GeometricVofInitialization::prepare(*mixedMesh,false,[](Vector2D){return 1.;});
        rejects([&]{VofSmoothGeometry::prepare(*mixedMesh,mixedSource);},VofSmoothGeometryFailure::UnsupportedMesh,
            "unqualified coarse-fine geometry silently accepted");
        DictionaryMesh cornerConfig;cornerConfig.Nx=cornerConfig.Ny=2;cornerConfig.maxCellNum=64;cornerConfig.maxInitialIterations=0;
        QuadtreeMesh cornerMesh(cornerConfig);
        const auto cornerSource=GeometricVofInitialization::prepare(cornerMesh,false,[](Vector2D p){return p.x+p.y-1.;});
        rejects([&]{VofSmoothGeometry::prepare(cornerMesh,cornerSource);},VofSmoothGeometryFailure::CornerEndpoint,
            "native-corner crossing silently acquired an arbitrary branch");
        circle();nativeStrips();requestedStripVolumeResidual();anchoredPhaseMeasure();boundedCenter();retainedNativeVolume();centerNativeVolume();
        std::cout<<"vof_smooth_geometry_contract PASS assertions="<<assertions<<'\n';return 0;
    } catch(const std::exception& e){std::cerr<<"vof_smooth_geometry_contract FAIL "<<e.what()<<'\n';return 1;}
}
