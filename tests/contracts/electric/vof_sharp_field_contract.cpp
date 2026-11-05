#include "Modules/Electric/ElectricVofSharpField.h"
#include "Modules/Electric/ElectricSurfaceChargeTransport.h"
#include "Modules/Interface/GeometricVofInitialization.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Geometry/PlanarCellCut.h"
#include "Core/Base/ScalarExpression.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>

namespace
{
using namespace LAMP::Electric;
using namespace LAMP::Interface;
unsigned assertions=0;
void check(bool ok,const char* reason){++assertions;if(!ok)throw std::runtime_error(reason);}
void close(double a,double b,double tolerance,const char* reason)
{if(std::abs(a-b)>tolerance*(1+std::abs(b)))std::cerr<<reason<<" actual="<<a<<" expected="<<b<<" tolerance="<<tolerance<<'\n';
 check(std::abs(a-b)<=tolerance*(1+std::abs(b)),reason);}
template<class F>void rejects(F f,const char* reason)
{bool rejected=false;try{f();}catch(const std::exception&){rejected=true;}check(rejected,reason);}
ElectricVofSharpProblem problem(bool rz,double beta,double gl,double gg,double a)
{
    ElectricVofSharpProblem p;
    for(auto& b:p.boundaries)b.value=[=](Vector2D x){return beta*x.x+(x.y<=a?gl:gg)*(x.y-a);};
    if(rz)p.boundaries[0].type=ElectricVofSharpBoundary::Type::HomogeneousNeumann;
    return p;
}
void layered(QuadtreeMesh& mesh,bool rz)
{
    constexpr double a=.423,gl=.8,qadv=.2,dt=.07;
    const double beta=rz?0:.37;
    const auto vof=GeometricVofInitialization::prepare(mesh,rz,[=](Vector2D x){return x.y-a;});
    VofSmoothGeometryConfig gc;gc.allowBoundaryEndpoints=true;
    const auto g=VofSmoothGeometry::prepare(mesh,vof,gc);
    auto p=problem(rz,beta,gl,1.4,a);p.epsilon={{2.3,.9}};p.conductivity={{.7,2.1}};
    const double density=p.epsilon[0]*gl-p.epsilon[1]*1.4;
    const auto q=ElectricVofSurfaceCharge::fromDensity(mesh,g,[=](std::size_t,double){return density;});
    const auto field=ElectricVofSharpField::prepareElectrostatic(mesh,g,q,p);
    check(field.chargeCandidate().sameInstance(q),"electrostatic snapshot replaced independent charge owner");
    check(field.report().usedSuperLU&&field.report().interfaceFaces==g.arcs().size(),"not a real global B1 SuperLU solve");
    if(!rz){
        namespace Cut=LAMP::Geometry::PlanarCellCut;
        std::size_t framedQueries=0;const auto saved=q.moments();
        for(const auto& cell:g.cells())if(cell.rectangle.lower.x!=0||cell.rectangle.lower.y!=0)
            for(unsigned phase=0;phase<2;++phase){const auto& points=cell.phaseQuadrature[phase];
                if(points.empty())continue;
                for(const auto index:{std::size_t{0},points.size()/2,points.size()-1}){
                    const auto point=points[index].point;
                    const auto e=field.sampleElectricField(cell.nativeCell,phase,
                        {static_cast<double>(point.x),static_cast<double>(point.y)});
                    const Cut::AnchoredPoint world{{static_cast<double>(point.x),static_cast<double>(point.y)},{0,0}};
                    const Cut::AnchoredPoint fromZero{{0,0},world.anchor};
                    const Cut::AnchoredPoint fromCell{{cell.rectangle.lower.x,cell.rectangle.lower.y},
                        {world.anchor.x-cell.rectangle.lower.x,world.anchor.y-cell.rectangle.lower.y}};
                    check(Cut::sameAnchoredPoint(world,fromZero)&&Cut::sameAnchoredPoint(world,fromCell),
                        "anchored field fixture changed its actual physical site");
                    const auto atZero=field.sampleElectricFieldAt(cell.nativeCell,phase,fromZero);
                    const auto atCell=field.sampleElectricFieldAt(cell.nativeCell,phase,fromCell);
                    check(e.x==atZero.x&&e.y==atZero.y&&e.x==atCell.x&&e.y==atCell.y,
                        "equivalent anchored field queries changed the stored P3 sample");
                    close(e.x,-beta,2e-8,"framed XY volume query tangential field");
                    close(e.y,phase==0?-gl:-1.4,2e-8,"framed XY volume query normal field");++framedQueries;
                }
            }
        check(framedQueries>0,"XY field did not query a nonzero actual element origin");
        const auto actual=field.candidateTime(),expected=g.sourceTime();
        check(actual.step==expected.step&&actual.begin==expected.begin&&actual.end==expected.end
            &&actual.position==expected.position,"framed field changed the genuine complete clock");
        check(field.geometry().sameInstance(g)&&field.sourceCharge().sameInstance(q)
            &&field.chargeCandidate().sameInstance(q)&&q.moments()==saved,"framed query changed geometry/charge owners or state");
        check(field.report().maximumRelativeGaussResidual<=1e-10
            &&field.report().maximumRelativeNativeBalance<=1e-10
            &&field.report().relativeGlobalGaussResidual<=1e-10,"framed XY field changed original conservation gates");
        for(const auto& cell:g.cells())if(cell.inventory[0]>0||cell.inventory[1]>0){
            const unsigned phase=cell.inventory[0]>0?0:1;
            const Vector2D outside{std::nextafter(cell.rectangle.upper.x,std::numeric_limits<double>::infinity()),
                cell.rectangle.lower.y};
            rejects([&]{field.sampleElectricField(cell.nativeCell,phase,outside);},"framed query accepted a world point outside actual native cell");
            break;
        }
        bool checkedTiny=false;
        for(const auto& cell:g.cells())if(cell.rectangle.lower.x>0&&cell.rectangle.lower.y>0){
            const auto& r=cell.rectangle;const unsigned phase=cell.inventory[0]>0?0:1;
            const long double x=(static_cast<long double>(r.lower.x)+r.upper.x)/2,
                y=(static_cast<long double>(r.lower.y)+r.upper.y)/2,tiny=std::ldexp(1.L,-80);
            const std::array<Cut::AnchoredPoint,4> outside{{
                {{r.lower.x,y},{-tiny,0}},{{r.upper.x,y},{tiny,0}},
                {{x,r.lower.y},{0,-tiny}},{{x,r.upper.y},{0,tiny}}}};
            for(const auto& site:outside){const auto rounded=Cut::approximateAnchoredPoint(site);
                check(rounded.x==site.anchor.x&&rounded.y==site.anchor.y,
                    "tiny exterior field query no longer rounds back to its native face");
                rejects([&]{field.sampleElectricFieldAt(cell.nativeCell,phase,site);},
                    "anchored field query lost a tiny independent exterior face distance");}
            checkedTiny=true;break;
        }
        check(checkedTiny,"anchored field fixture lacks four nonzero-coordinate native faces");
        std::cout<<"framed XY actual_volume_queries="<<framedQueries<<" owners/clock/Gauss=PASS\n";
    }
    double maximum=0;
    for(std::size_t arc=0;arc<g.arcs().size();++arc){
        for(double t:{.17,.5,.83}){const auto value=field.interfaceField(arc,t);
            close(value.normalElectricField[0],-gl,2e-8,"liquid sharp normal field");
            close(value.normalElectricField[1],-1.4,2e-8,"gas sharp normal field");
            const double exactEt=-beta*value.tangent.x;
            close(value.tangentialElectricField,exactEt,2e-8,"shared nonzero tangential trace");
            close(p.epsilon[1]*value.normalElectricField[1]-p.epsilon[0]*value.normalElectricField[0],
                value.surfaceChargeDensity,1e-9,"independent Q Gauss jump");
            const auto point=g.arcs()[arc].sample(t);
            for(unsigned phase=0;phase<2;++phase){const auto e=field.sampleElectricField(g.arcs()[arc].nativeCell,phase,
                    {static_cast<double>(point.x),static_cast<double>(point.y)});
                close(e.x,-beta,2e-8,"recovered native P3 tangential field");
                close(e.y,phase==0?-gl:-1.4,2e-8,"recovered native P3 normal field");}
            maximum=std::max(maximum,std::abs(value.tangentialElectricField-exactEt));}}
    // Nonproportional epsilon/K and an independently specified initial Q.
    // Target boundary data and analytic layer slopes solve the same BE system.
    const double al=p.epsilon[0]+dt*p.conductivity[0],ag=p.epsilon[1]+dt*p.conductivity[1];
    const double gg=(al*gl-qadv)/ag;
    auto transient=problem(rz,beta,gl,gg,a);transient.epsilon=p.epsilon;transient.conductivity=p.conductivity;
    const auto initial=ElectricVofSurfaceCharge::fromDensity(mesh,g,[=](std::size_t,double){return qadv;});
    LAMP::TimeLayer end{1,0,dt,LAMP::TimeLayer::Position::End};
    const auto ohmic=ElectricVofSharpField::prepareOhmic(mesh,g,initial,end,transient);
    check(ohmic.report().relativeGlobalChargeUpdateResidual<=1e-10,"independent boundary current charge ledger");
    check(ohmic.report().relativeGlobalGaussResidual<=1e-10,"independent global Gauss ledger");
    const auto& updated=ohmic.chargeCandidate();check(!updated.sameInstance(initial),"Ohmic solve did not issue candidate owner");
    check(updated.timeLayer().position==LAMP::TimeLayer::Position::End&&updated.timeLayer().end==dt,"candidate full clock missing");
    const double exactQ=p.epsilon[0]*gl-p.epsilon[1]*gg;
    for(std::size_t arc=0;arc<g.arcs().size();++arc){const auto value=ohmic.interfaceField(arc,.37);
        close(value.normalElectricField[0],-gl,2e-8,"implicit liquid field");
        close(value.normalElectricField[1],-gg,2e-8,"implicit gas field");
        close(value.surfaceChargeDensity,exactQ,2e-8,"conductive surface charge candidate");
        const auto start=g.arcs()[arc].sample(0),finish=g.arcs()[arc].sample(1);
        const double length=static_cast<double>(std::hypot(finish.x-start.x,finish.y-start.y));
        const double measure=rz?std::acos(-1.)*static_cast<double>(finish.x+start.x)*length:length;
        close(updated.integratedCharge()[arc],exactQ*measure,2e-8,"actual XY/RZ integrated surface charge");}
    check(initial.integratedCharge()==ElectricVofSurfaceCharge::fromDensity(mesh,g,[=](std::size_t,double){return qadv;}).integratedCharge(),
        "field changed accepted charge state");
    const auto foreign=VofSmoothGeometry::prepare(mesh,vof,gc);
    rejects([&]{q.validateFor(mesh,foreign,q.timeLayer());},"same-size foreign smooth owner accepted");
    rejects([&]{ElectricVofSurfaceCharge::fromMoments(mesh,g,end,q.moments());},"arbitrary initial clock relabel accepted");
    auto invalid=transient;invalid.epsilon[0]=0;
    rejects([&]{ElectricVofSharpField::prepareOhmic(mesh,g,initial,end,invalid);},"invalid material accepted");
    auto skipped=end;skipped.step=2;
    rejects([&]{ElectricVofSharpField::prepareOhmic(mesh,g,initial,skipped,transient);},"noncontiguous candidate step accepted");
    auto wrongBegin=end;wrongBegin.begin=.01;
    rejects([&]{ElectricVofSharpField::prepareOhmic(mesh,g,initial,wrongBegin,transient);},"different actual Begin accepted");
    std::vector<std::array<long double,2>> stocks;for(const auto& cell:vof.plicCells())stocks.push_back(cell.inventory);
    const auto largest=std::numeric_limits<std::int64_t>::max();
    const LAMP::TimeLayer finalClock{largest,0,1,LAMP::TimeLayer::Position::End};
    const auto finalVof=GeometricVofPublication::prepareLocalPlicInventory(mesh,vof.domain(),finalClock,stocks,{});
    const auto finalGeometry=VofSmoothGeometry::prepare(mesh,finalVof,gc);
    const auto finalCharge=ElectricVofSurfaceCharge::fromDensity(mesh,finalGeometry,[=](std::size_t,double){return qadv;});
    const LAMP::TimeLayer overflow{largest,1,2,LAMP::TimeLayer::Position::End};
    rejects([&]{ElectricVofSharpField::prepareOhmic(mesh,finalGeometry,finalCharge,overflow,transient);},"source step overflow accepted");
    auto shifted=p;for(auto& b:shifted.boundaries){const auto provider=b.value;
        b.value=[provider](Vector2D x){return provider(x)+1e6;};}
    const auto gauge=ElectricVofSharpField::prepareElectrostatic(mesh,g,q,shifted);
    for(std::size_t arc=0;arc<g.arcs().size();++arc){const auto value=gauge.interfaceField(arc,.5);
        close(value.normalElectricField[0],-gl,2e-8,"large constant voltage changed physical field");
        close(value.surfaceChargeDensity,density,2e-8,"large voltage offset changed independent surface charge");}
    std::cout<<"layered rz="<<rz<<" elements="<<field.report().elements<<" dofs="<<field.report().skeletonDofs
        <<" Et_error="<<maximum<<" gauss="<<ohmic.report().maximumRelativeGaussResidual<<'\n';
}
void anchoredQueryKernel()
{
    namespace Cut=LAMP::Geometry::PlanarCellCut;
    DictionaryMesh input;input.Nx=input.Ny=8;input.maxCellNum=512;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);
    const auto vof=GeometricVofInitialization::prepare(mesh,false,[](Vector2D point){return point.y-.423;});
    VofSmoothGeometryConfig config;config.allowBoundaryEndpoints=true;
    const auto geometry=VofSmoothGeometry::prepare(mesh,vof,config);
    const auto charge=ElectricVofSurfaceCharge::fromDensity(mesh,geometry,[](std::size_t,double){return 0.;});
    const auto before=charge.moments();constexpr double center=.5625;
    ElectricVofSharpProblem p;p.epsilon={{1,1}};
    for(auto& boundary:p.boundaries)boundary.value=[=](Vector2D point){return (point.x-center)*(point.y-center);};
    const auto field=ElectricVofSharpField::prepareElectrostatic(mesh,geometry,charge,p);
    int native=-1;unsigned phase=0;
    for(const auto& cell:geometry.cells())if(cell.rectangle.lower.x==.5&&cell.rectangle.lower.y==.5){
        native=cell.nativeCell;phase=cell.inventory[0]>0?0:1;
        check(cell.inventory[1-phase]==0,"anchored harmonic fixture lost its pure symmetric cell");break;
    }
    check(native>=0,"anchored harmonic fixture lacks its selected actual cell");
    const long double tiny=std::ldexp(1.L,-80);
    const Cut::AnchoredPoint origin{{center,center},{0,0}},shifted{{center,center},{tiny,0}};
    check(Cut::approximateAnchoredPoint(shifted).x==center&&!Cut::sameAnchoredPoint(origin,shifted),
        "anchored harmonic fixture does not retain a sub-world-ulp displacement");
    const auto initial=field.sampleElectricFieldAt(native,phase,origin);
    const auto displaced=field.sampleElectricFieldAt(native,phase,shifted);
    const Cut::AnchoredPoint shiftedSwapped{shifted.offset,shifted.anchor};
    check(Cut::sameAnchoredPoint(shifted,shiftedSwapped),"harmonic swapped-anchor fixture changed its actual site");
    const auto displacedSwapped=field.sampleElectricFieldAt(native,phase,shiftedSwapped);
    check(displaced.x==displacedSwapped.x&&displaced.y==displacedSwapped.y,
        "actual P3 field changed under equivalent complete-frame decomposition");
    const auto world=field.sampleElectricField(native,phase,{center,center});
    check(world.x==initial.x&&world.y==initial.y,"world field query did not consume the anchored kernel");
    // For phi=(x-c)(y-c), the physical increment is delta Ey=-delta x.
    // Include rounding of the actual recovered baseline, rather than using
    // the much larger global field error as a tolerance for this increment.
    const long double budget=1e-8L*tiny+16*std::numeric_limits<double>::epsilon()
        *(std::abs(static_cast<long double>(initial.y))+std::abs(static_cast<long double>(displaced.y)));
    check(budget<tiny/4,"actual harmonic field lacks resolution for its retained query displacement");
    const long double increment=static_cast<long double>(displaced.y)-initial.y;
    check(std::abs(increment+tiny)<=budget,"anchored P3 basis discarded the retained query offset");
    check(field.geometry().sameInstance(geometry)&&field.sourceCharge().sameInstance(charge)
        &&field.chargeCandidate().sameInstance(charge)&&charge.moments()==before,
        "anchored harmonic query changed accepted sources or field owners");
    check(field.report().relativeLinearResidual<=1e-11&&field.report().maximumRelativeGaussResidual<=1e-10
        &&field.report().maximumRelativeNativeBalance<=1e-10&&field.report().relativeGlobalGaussResidual<=1e-10,
        "anchored harmonic field changed an original physical/linear gate");
    rejects([&]{field.sampleElectricFieldAt(-1,phase,origin);},"anchored query accepted an invalid native cell");
    rejects([&]{field.sampleElectricFieldAt(native,2,origin);},"anchored query accepted an invalid phase");
    auto nonfinite=origin;nonfinite.offset.x=std::numeric_limits<long double>::infinity();
    rejects([&]{field.sampleElectricFieldAt(native,phase,nonfinite);},"anchored query accepted a nonfinite site component");
    std::cout<<std::setprecision(17)<<"anchored P3 actual_delta_Ey="<<increment<<" exact="<<-tiny
        <<" native_budget="<<budget<<" zero_world_displacement=1\n";
}
void curved(QuadtreeMesh& mesh,bool rz)
{
    const auto vof=GeometricVofInitialization::prepare(mesh,rz,[](Vector2D x){
        return x.y-(.421+.031*std::sin(2*std::acos(-1.)*x.x));});
    VofSmoothGeometryConfig gc;gc.allowBoundaryEndpoints=true;
    const auto g=VofSmoothGeometry::prepare(mesh,vof,gc);
    auto p=problem(rz,0,1,1,0);
    const auto q=ElectricVofSurfaceCharge::fromDensity(mesh,g,[](std::size_t,double){return 0.;});
    const auto field=ElectricVofSharpField::prepareElectrostatic(mesh,g,q,p);
    bool realCurve=false;double error=0;
    for(std::size_t a=0;a<g.arcs().size();++a){realCurve=realCurve||std::abs(g.arcs()[a].volumeCorrection)>1e-9;
        const auto value=field.interfaceField(a,.5);
        check(std::isfinite(value.tangentialElectricField),"nonfinite actual curved shared trace");
        for(unsigned phase=0;phase<2;++phase){error=std::max(error,std::hypot(value.electricField[phase].x,value.electricField[phase].y+1));}
        close(value.normalElectricField[1]-value.normalElectricField[0],value.surfaceChargeDensity,1e-9,
            "curved equilibrated Gauss jump");}
    check(realCurve,"curved B1 test became PLIC");
    check(error<.08,"curved harmonic field error exceeds frozen direct test");
    check(field.report().maximumRelativeNativeBalance<=1e-10,"curved actual native flux imbalance");
    std::cout<<"curved rz="<<rz<<" actual_E_error="<<error<<" gauss="<<field.report().maximumRelativeGaussResidual<<'\n';
}
void obliqueSmallCut(QuadtreeMesh& mesh)
{
    namespace Cut=LAMP::Geometry::PlanarCellCut;
    const double a=.423,slope=.37,beta=.37,gl=.8,gg=1.4,length=std::hypot(1.,slope);
    const Vector2D normal{-slope/length,1/length},tangent{1/length,slope/length};
    FaceSegmentGraph graph;graph.build(mesh);auto domain=NativeLeafDomain::build(mesh,graph,false);
    std::vector<std::array<long double,2>> stocks(mesh.cellIndexNum());
    for(int cell=0;cell<mesh.cellIndexNum();++cell){const auto& c=mesh.cells(mesh.cellIndexToLabel[cell]);
        const auto lo=mesh.verts(c.cellToVert[0]).coord(),hi=mesh.verts(c.cellToVert[2]).coord();
        const long double cx=(static_cast<long double>(lo.x)+hi.x)/2,cy=(static_cast<long double>(lo.y)+hi.y)/2;
        const long double offset=-(cy-a-slope*cx)/std::hypot(1.L,static_cast<long double>(slope));
        stocks[cell][0]=Cut::volume({lo,hi},{-slope,1},offset,false);
        stocks[cell][1]=Cut::volume({lo,hi},{slope,-1},-offset,false);}
    const auto publication=GeometricVofPublication::prepareLocalPlicInventory(mesh,std::move(domain),
        {0,0,0,LAMP::TimeLayer::Position::Initial},stocks,{});
    long double minimumFraction=1;
    for(const auto& cell:publication.plicCells())for(const auto stock:cell.inventory)if(stock>0)
        minimumFraction=std::min(minimumFraction,stock/Cut::capacity(cell.rectangle,false));
    check(minimumFraction<.005L,"declared oblique fixture did not retain a true small-cut phase");
    VofSmoothGeometryConfig gc;gc.allowBoundaryEndpoints=true;
    const auto geometry=VofSmoothGeometry::prepare(mesh,publication,gc);
    check(!geometry.arcs().empty(),"small-cut fixture lost its actual G1 interface");
    auto p=problem(false,0,gl,gg,a);p.epsilon={{2.3,.9}};p.conductivity={{.7,2.1}};
    for(auto& b:p.boundaries)b.value=[=](Vector2D point){const double d=(point.y-a-slope*point.x)/length;
        return beta*(tangent.x*point.x+tangent.y*point.y)+(d<=0?gl:gg)*d;};
    const double density=p.epsilon[0]*gl-p.epsilon[1]*gg;
    const auto charge=ElectricVofSurfaceCharge::fromDensity(mesh,geometry,[=](std::size_t,double){return density;});
    const auto before=charge.moments();const auto field=ElectricVofSharpField::prepareElectrostatic(mesh,geometry,charge,p);
    auto shifted=p;for(auto& b:shifted.boundaries){const auto provider=b.value;b.value=[provider](Vector2D point){return provider(point)+128;};}
    const auto gauge=ElectricVofSharpField::prepareElectrostatic(mesh,geometry,charge,shifted);
    check(field.geometry().sameInstance(geometry)&&gauge.geometry().sameInstance(geometry)
        &&field.sourceCharge().sameInstance(charge)&&gauge.sourceCharge().sameInstance(charge),"small-cut field/gauge lost actual owners");
    check(field.report().usedSuperLU&&gauge.report().usedSuperLU,"small-cut field/gauge did not execute the actual global solve");
    double maximumFieldError=0;
    for(std::size_t arc=0;arc<geometry.arcs().size();++arc)for(double t:{.17,.5,.83}){
        const auto point=geometry.arcs()[arc].sample(t);const Vector2D x{static_cast<double>(point.x),static_cast<double>(point.y)};
        close(x.y-a-slope*x.x,0,2e-8,"small-cut actual G1 plane");
        const auto value=field.interfaceField(arc,t),offsetValue=gauge.interfaceField(arc,t);
        close(value.normalElectricField[0],-gl,2e-8,"small-cut liquid normal field");
        close(value.normalElectricField[1],-gg,2e-8,"small-cut gas normal field");
        close(value.tangentialElectricField,-beta*(tangent.x*value.tangent.x+tangent.y*value.tangent.y),2e-8,"small-cut common tangent trace");
        close(p.epsilon[1]*value.normalElectricField[1]-p.epsilon[0]*value.normalElectricField[0],density,1e-9,"small-cut independent Gauss jump");
        close(value.surfaceChargeDensity,density,2e-8,"small-cut independent surface charge density");
        close(offsetValue.tangentialElectricField,value.tangentialElectricField,2e-8,"small-cut constant gauge tangent field");
        close(offsetValue.surfaceChargeDensity,value.surfaceChargeDensity,2e-8,"small-cut constant gauge charge density");
        for(unsigned phase=0;phase<2;++phase){const double grad=phase?gg:gl;
            const auto e=field.sampleElectricField(geometry.arcs()[arc].nativeCell,phase,x);
            const auto shiftedE=gauge.sampleElectricField(geometry.arcs()[arc].nativeCell,phase,x);
            const double ex=-beta*tangent.x-grad*normal.x,ey=-beta*tangent.y-grad*normal.y;
            close(e.x,ex,2e-8,"small-cut actual P3 field x");close(e.y,ey,2e-8,"small-cut actual P3 field y");
            close(shiftedE.x,e.x,2e-8,"small-cut constant gauge P3 field x");close(shiftedE.y,e.y,2e-8,"small-cut constant gauge P3 field y");
            close(offsetValue.normalElectricField[phase],value.normalElectricField[phase],2e-8,"small-cut constant gauge normal field");
            maximumFieldError=std::max(maximumFieldError,std::hypot(e.x-ex,e.y-ey));}}
    check(charge.moments()==before&&field.chargeCandidate().sameInstance(charge)&&gauge.chargeCandidate().sameInstance(charge),
        "small-cut electrostatic/gauge changed independent accepted charge");
    const double qadv=.2,dt=.07;
    const double implicitGas=((p.epsilon[0]+dt*p.conductivity[0])*gl-qadv)/(p.epsilon[1]+dt*p.conductivity[1]);
    auto transient=p;for(auto& b:transient.boundaries)b.value=[=](Vector2D point){const double d=(point.y-a-slope*point.x)/length;
        return beta*(tangent.x*point.x+tangent.y*point.y)+(d<=0?gl:implicitGas)*d;};
    const auto initial=ElectricVofSurfaceCharge::fromDensity(mesh,geometry,[=](std::size_t,double){return qadv;});
    const auto sourceMoments=initial.moments();const LAMP::TimeLayer end{1,0,dt,LAMP::TimeLayer::Position::End};
    const auto ohmic=ElectricVofSharpField::prepareOhmic(mesh,geometry,initial,end,transient);
    auto transientGauge=transient;for(auto& b:transientGauge.boundaries){const auto provider=b.value;
        b.value=[provider](Vector2D point){return provider(point)+128;};}
    const auto ohmicGauge=ElectricVofSharpField::prepareOhmic(mesh,geometry,initial,end,transientGauge);
    check(initial.moments()==sourceMoments&&ohmic.sourceCharge().sameInstance(initial)&&ohmicGauge.sourceCharge().sameInstance(initial),
        "small-cut Ohmic/gauge replaced or changed accepted Qadv");
    const double nextDensity=p.epsilon[0]*gl-p.epsilon[1]*implicitGas;
    for(std::size_t arc=0;arc<geometry.arcs().size();++arc){const auto value=ohmic.interfaceField(arc,.37),offsetValue=ohmicGauge.interfaceField(arc,.37);
        close(value.normalElectricField[0],-gl,2e-8,"small-cut implicit liquid field");
        close(value.normalElectricField[1],-implicitGas,2e-8,"small-cut implicit gas field");
        close(value.surfaceChargeDensity,nextDensity,2e-8,"small-cut independent Ohmic candidate charge");
        close(p.epsilon[1]*value.normalElectricField[1]-p.epsilon[0]*value.normalElectricField[0],nextDensity,1e-9,"small-cut implicit physical Gauss jump");
        close(value.tangentialElectricField,-beta*(tangent.x*value.tangent.x+tangent.y*value.tangent.y),2e-8,"small-cut implicit common tangent field");
        close(offsetValue.surfaceChargeDensity,value.surfaceChargeDensity,2e-8,"small-cut Ohmic constant gauge candidate charge");
        for(unsigned phase=0;phase<2;++phase)close(offsetValue.normalElectricField[phase],value.normalElectricField[phase],2e-8,
            "small-cut Ohmic constant gauge normal field");
        close(offsetValue.tangentialElectricField,value.tangentialElectricField,2e-8,"small-cut Ohmic constant gauge tangential field");
        const auto begin=geometry.arcs()[arc].sample(0),finish=geometry.arcs()[arc].sample(1);
        close(ohmic.chargeCandidate().integratedCharge()[arc],nextDensity*static_cast<double>(std::hypot(finish.x-begin.x,finish.y-begin.y)),
            2e-8,"small-cut Ohmic actual integrated interface charge");}
    check(ohmic.report().maximumRelativeGaussResidual<=1e-10&&ohmic.report().relativeGlobalGaussResidual<=1e-10
        &&ohmic.report().relativeGlobalChargeUpdateResidual<=1e-10&&ohmicGauge.report().maximumRelativeGaussResidual<=1e-10
        &&ohmicGauge.report().relativeGlobalChargeUpdateResidual<=1e-10,"small-cut Ohmic/gauge original physical conservation gates");
    std::cout<<"small-cut minimum_phase_fraction="<<minimumFraction<<" E_error="<<maximumFieldError
        <<" gauss="<<field.report().maximumRelativeGaussResidual<<" gauge_gauss="<<gauge.report().maximumRelativeGaussResidual
        <<" ohmic_gauss="<<ohmic.report().maximumRelativeGaussResidual<<" ohmic_gauge_gauss="<<ohmicGauge.report().maximumRelativeGaussResidual<<'\n';
}
void sameOperatorDefectCorrection()
{
    DictionaryMesh config;config.Nx=config.Ny=32;config.maxCellNum=4096;config.maxInitialIterations=0;QuadtreeMesh mesh(config);
    const LAMP::ScalarExpression::Constants constants{{"centerX",.487},{"centerY",.513},{"radius",.18},{"chargeAmplitude",.03}};
    const LAMP::ScalarExpression phi("(x-centerX)^2+(y-centerY)^2-radius^2",constants,LAMP::ScalarExpression::Coordinates::Planar);
    const LAMP::ScalarExpression density("chargeAmplitude*(x-centerX)",constants,LAMP::ScalarExpression::Coordinates::Planar);
    const auto vof=GeometricVofInitialization::prepare(mesh,false,[&](Vector2D point){return phi.sample(point.x,point.y).value;});
    std::vector<std::array<long double,2>> sourceStocks;
    for(const auto& cell:vof.plicCells())sourceStocks.push_back(cell.inventory);
    const auto sourceGeometry=VofSmoothGeometry::prepare(mesh,vof,{});
    const auto source=ElectricVofSurfaceCharge::fromDensity(mesh,sourceGeometry,[&](std::size_t arc,double t){
        const auto point=sourceGeometry.arcs()[arc].sample(t);return density.sample(static_cast<double>(point.x),static_cast<double>(point.y)).value;});
    const auto before=source.moments();const LAMP::TimeLayer end{1,0,.0001,LAMP::TimeLayer::Position::End};
    std::vector<double> zero(static_cast<std::size_t>(vof.domain().nativeFaceCount()));
    const auto slab=GeometricVofAdvectorSlab::prepare(mesh,vof,end,zero,zero);
    const auto transport=GeometricVofTransport::prepareLocalPlicSweep(mesh,vof,slab,2,1);
    const auto target=VofSmoothGeometry::prepare(mesh,transport.midpoint(),{});
    const auto moved=ElectricSurfaceChargeTransport::prepareStage(mesh,source,transport,target);
    ElectricVofSharpProblem p;p.epsilon={{2.3,.9}};p.conductivity={{.7,2.1}};
    p.boundaries[0].value=[](Vector2D){return -.02;};p.boundaries[2].value=[](Vector2D){return .02;};
    p.boundaries[1].type=p.boundaries[3].type=ElectricVofSharpBoundary::Type::HomogeneousNeumann;
    const auto field=ElectricVofSharpField::prepareOhmic(mesh,moved,p);const auto& r=field.report();
    long double minimumFraction=1;for(const auto& cell:target.cells())for(const auto stock:cell.inventory)if(stock>0)
        minimumFraction=std::min(minimumFraction,stock/LAMP::Geometry::PlanarCellCut::capacity(cell.rectangle,false));
    // Report before the assertions, so either legitimate production branch
    // leaves its actual initial/final residuals in the validation receipt.
    std::cout<<std::setprecision(17)<<"same-C small-cut="<<minimumFraction<<" corrections="<<r.defectCorrections
        <<" residual_measured="<<r.defectResidualMeasured<<" correction_seconds="<<r.defectCorrectionSeconds
        <<" uncondensed_before="<<r.uncondensedResidualBefore<<" uncondensed_after="<<r.uncondensedResidualAfter
        <<" initial_Gamma="<<r.initialMaximumRelativeGaussResidual<<" final_Gamma="<<r.maximumRelativeGaussResidual
        <<" final_local="<<r.maximumRelativeLocalResidual<<" final_native="<<r.maximumRelativeNativeBalance
        <<" final_update="<<r.maximumRelativeChargeUpdateResidual<<std::endl;
    check(minimumFraction<1e-5L,"same-operator regression lost the actual small phase");
    check(std::isfinite(r.initialMaximumRelativeGaussResidual)&&r.initialMaximumRelativeGaussResidual>=0,
        "same-operator initial physical Gamma residual is unavailable");
    if(r.defectCorrections==0){
        check(r.initialMaximumRelativeGaussResidual<=1e-10
            &&r.initialMaximumRelativeGaussResidual==r.maximumRelativeGaussResidual,
            "uncorrected same-operator solve did not already satisfy the original Gamma gate");
        check(!r.defectResidualMeasured&&r.uncondensedResidualBefore==0&&r.uncondensedResidualAfter==0
            &&r.defectCorrectionSeconds==0,"uncorrected solve forged defect-correction diagnostics");
    }else{
        check(r.defectCorrections<=2&&r.defectResidualMeasured,"same-operator correction exceeded its original bounded policy");
        check(std::isfinite(r.uncondensedResidualBefore)&&std::isfinite(r.uncondensedResidualAfter)
            &&r.uncondensedResidualBefore>0&&r.uncondensedResidualAfter>=0
            &&r.uncondensedResidualAfter<r.uncondensedResidualBefore,
            "same-C defect correction did not reduce the original measured uncondensed residual");
        check(std::isfinite(r.defectCorrectionSeconds)&&r.defectCorrectionSeconds>=0,"invalid actual correction timing");
    }
    check(r.relativeLinearResidual<=1e-11&&r.maximumRelativeGaussResidual<=1e-10
        &&r.maximumRelativeLocalResidual<=1e-10&&r.maximumRelativeNativeBalance<=1e-10
        &&r.maximumRelativeChargeUpdateResidual<=1e-10
        &&r.relativeGlobalGaussResidual<=1e-10&&r.relativeGlobalChargeUpdateResidual<=1e-10,
        "same-operator correction failed an unchanged original gate");
    check(source.moments()==before&&field.sourceCharge().sameInstance(source)&&field.implicitInputCharge().sameInstance(moved.chargeCandidate())
        &&field.advection()&&field.advection()->interfaceTransport()->sameInstance(transport),"same-operator correction changed accepted Q or actual transport provenance");
    check(field.geometry().sameInstance(target)&&field.chargeCandidate().geometry().sameInstance(target),
        "same-operator field or candidate replaced the actual target geometry");
    for(std::size_t cell=0;cell<sourceStocks.size();++cell)
        check(vof.plicCells()[cell].inventory==sourceStocks[cell],"same-operator solve changed accepted independent phase stocks");
}
}
int main(int argc,char** argv)
{
    try{DictionaryMesh input;input.Nx=input.Ny=8;input.maxCellNum=512;input.maxInitialIterations=0;
        QuadtreeMesh mesh(input);
        if(argc==2&&std::string(argv[1])=="--anchored-query-only"){
            anchoredQueryKernel();std::cout<<"vof_sharp_field_contract anchored query PASS assertions="<<assertions<<'\n';return 0;
        }
        if(argc==2&&std::string(argv[1])=="--framed-xy-only"){
            layered(mesh,false);std::cout<<"vof_sharp_field_contract framed XY PASS assertions="<<assertions<<'\n';return 0;
        }
        check(argc==1,"usage: vof_sharp_field_contract [--framed-xy-only | --anchored-query-only]");
        for(bool rz:{false,true}){layered(mesh,rz);curved(mesh,rz);}obliqueSmallCut(mesh);sameOperatorDefectCorrection();anchoredQueryKernel();
        std::cout<<"vof_sharp_field_contract PASS assertions="<<assertions<<'\n';return 0;
    }catch(const std::exception& e){std::cerr<<"vof_sharp_field_contract FAIL "<<e.what()<<'\n';return 1;}
}
