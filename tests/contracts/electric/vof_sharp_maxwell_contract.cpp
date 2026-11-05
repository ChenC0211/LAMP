#include "Modules/Electric/ElectricVofSharpMaxwellStress.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/CartesianMeasure.h"
#include "Core/Geometry/PlanarCellCut.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

namespace
{
using namespace LAMP::Interface;
using namespace LAMP::Electric;
unsigned checks=0,scenarios=0;
void check(bool ok,const char* message){++checks;if(!ok)throw std::runtime_error(message);}
void near(double a,double b,double scale,const char* message)
{check(std::isfinite(a)&&std::abs(a-b)<=1e-8*scale,message);}
void scenario(bool rz,bool curved)
{
    DictionaryMesh config;config.Nx=config.Ny=4;config.maxCellNum=256;config.maxInitialIterations=0;
    if(rz){config.minCoord={1,0};config.maxCoord={2,1};}
    QuadtreeMesh mesh(config);FaceSegmentGraph graph;graph.build(mesh);
    auto domain=NativeLeafDomain::build(mesh,graph,rz);
    const double center=rz?1.5:.5,height=.37;
    const auto shape=[&](double x){return height+(curved?.08*(x-center)*(x-center):0);};
    // Independent fixed Gauss integration of the prescribed quadratic graph.
    constexpr std::array<double,4> gx{{-.8611363115940525752,-.3399810435848562648,.3399810435848562648,.8611363115940525752}};
    constexpr std::array<double,4> gw{{.3478548451374538574,.6521451548625461426,.6521451548625461426,.3478548451374538574}};
    std::vector<std::array<long double,2>> stocks(mesh.cellIndexNum());
    for(int cell=0;cell<mesh.cellIndexNum();++cell){const auto& c=mesh.cells(mesh.cellIndexToLabel[cell]);
        const auto lo=mesh.verts(c.cellToVert[0]).coord(),hi=mesh.verts(c.cellToVert[2]).coord();
        const double minimum=height,maximum=std::max(shape(lo.x),shape(hi.x));
        const long double capacity=LAMP::Geometry::PlanarCellCut::capacity({lo,hi},rz);
        if(hi.y<=minimum){stocks[cell]={capacity,0};continue;}
        if(lo.y>=maximum){stocks[cell]={0,capacity};continue;}
        for(unsigned q=0;q<gx.size();++q){const double x=(lo.x+hi.x)/2+(hi.x-lo.x)*gx[q]/2;
            const long double factor=(hi.x-lo.x)*gw[q]/2*(rz?CartesianMeasure::twoPi*x:1.L);
            stocks[cell][0]+=factor*std::max(0.,std::min(shape(x),hi.y)-lo.y);
            stocks[cell][1]+=factor*std::max(0.,hi.y-std::max(shape(x),lo.y));}}
    auto publication=GeometricVofPublication::prepareLocalPlicInventory(mesh,std::move(domain),{0,0,0,LAMP::TimeLayer::Position::Initial},stocks,{});
    VofSmoothGeometryConfig gc;gc.allowBoundaryEndpoints=true;
    auto geometry=VofSmoothGeometry::prepare(mesh,publication,gc);
    const double enLiquid=.4,enGas=.7,et=rz?0:.3;
    ElectricVofSharpProblem problem;problem.epsilon={2,5};
    for(auto& b:problem.boundaries)b.value=[=](Vector2D p){return -(p.y-height)*(p.y<height?enLiquid:enGas)-et*p.x;};
    auto charge=ElectricVofSurfaceCharge::fromDensity(mesh,geometry,[&](std::size_t arc,double t){
        if(!curved)return problem.epsilon[1]*enGas-problem.epsilon[0]*enLiquid;
        const auto n=geometry.arcs()[arc].normal(t);
        return static_cast<double>((problem.epsilon[1]-problem.epsilon[0])*(et*n.x+enLiquid*n.y));});
    const auto before=charge.moments();
    auto field=ElectricVofSharpField::prepareElectrostatic(mesh,geometry,charge,problem);
    auto stress=ElectricVofSharpMaxwellStress::prepare(mesh,field);stress.validateFor(mesh,field);
    check(stress.axisymmetric()==rz,"wrong physical stress metric");
    check(field.sourceCharge().sameInstance(charge),"stress field lost accepted charge owner");
    check(field.geometry().sameInstance(geometry),"stress field lost actual smooth geometry owner");
    check(charge.moments()==before,"read-only stress changed surface charge");
    namespace Cut=LAMP::Geometry::PlanarCellCut;
    std::size_t anchoredSamples=0;
    for(const auto& cell:geometry.cells())for(unsigned phase=0;phase<2;++phase)
        if(!cell.phaseQuadrature[phase].empty()){
            const auto node=cell.phaseQuadrature[phase].front().point;
            const Vector2D world{static_cast<double>(node.x),static_cast<double>(node.y)};
            const Cut::AnchoredPoint site{{world.x,world.y},{0,0}},fromZero{{0,0},{world.x,world.y}};
            const Cut::AnchoredPoint fromCell{{cell.rectangle.lower.x,cell.rectangle.lower.y},
                {static_cast<long double>(world.x)-cell.rectangle.lower.x,static_cast<long double>(world.y)-cell.rectangle.lower.y}};
            check(Cut::sameAnchoredPoint(site,fromZero)&&Cut::sameAnchoredPoint(site,fromCell),
                "anchored stress fixture changed its physical site");
            const auto legacy=stress.sample(cell.nativeCell,phase,world);
            const auto a=stress.sampleAt(cell.nativeCell,phase,fromZero),b=stress.sampleAt(cell.nativeCell,phase,fromCell);
            for(const auto* sample:{&a,&b}){
                check(sample->electricField.x==legacy.electricField.x&&sample->electricField.y==legacy.electricField.y,
                    "anchored stress changed the actual stored electric field");
                check(sample->tensor.xx==legacy.tensor.xx&&sample->tensor.xy==legacy.tensor.xy
                    &&sample->tensor.yy==legacy.tensor.yy&&sample->tensor.hoop==legacy.tensor.hoop,
                    "anchored stress changed its material tensor kernel");
            }
            ++anchoredSamples;
        }
    check(anchoredSamples>0,"stress had no actual anchored production queries");
    bool checkedTiny=false;
    for(const auto& cell:geometry.cells())if(cell.rectangle.lower.x>0&&cell.rectangle.lower.y>0){
        const auto& r=cell.rectangle;const unsigned phase=cell.inventory[0]>0?0:1;
        const long double x=(static_cast<long double>(r.lower.x)+r.upper.x)/2,
            y=(static_cast<long double>(r.lower.y)+r.upper.y)/2,tiny=std::ldexp(1.L,-80);
        const std::array<Cut::AnchoredPoint,4> sites{{
            {{r.lower.x,y},{-tiny,0}},{{r.upper.x,y},{tiny,0}},
            {{x,r.lower.y},{0,-tiny}},{{x,r.upper.y},{0,tiny}}}};
        for(const auto& site:sites){const auto rounded=Cut::approximateAnchoredPoint(site);
            check(rounded.x==site.anchor.x&&rounded.y==site.anchor.y,"tiny stress fixture no longer folds to its cell face");
            bool rejectedSite=false;try{stress.sampleAt(cell.nativeCell,phase,site);}
            catch(const std::exception&){rejectedSite=true;}
            check(rejectedSite,"anchored stress accepted a tiny exterior site");}
        checkedTiny=true;break;
    }
    check(checkedTiny,"stress fixture lacks four positive-coordinate native faces");
    for(std::size_t arc=0;arc<geometry.arcs().size();++arc)for(double t:{.13,.5,.87}){
        const auto s=stress.interfaceSample(arc,t);
        check(s.relativeNormalIdentityResidual<=1e-12,"normal tensor identity");
        check(s.relativeTangentialJumpResidual<=1e-10,"charge/tangential stress jump");
        for(unsigned phase=0;phase<2;++phase){const auto e=s.field.electricField[phase];const double eps=problem.epsilon[phase];
            near(s.tensor[phase].hoop,-.5*eps*(e.x*e.x+e.y*e.y),eps*(e.x*e.x+e.y*e.y),"missing real hoop component");}
        if(!curved){near(s.field.normalElectricField[0],enLiquid,1,"liquid normal field analytic error");
            near(s.field.normalElectricField[1],enGas,1,"gas normal field analytic error");
            // Arc orientation has liquid on its left and may oppose +x.
            near(std::abs(s.field.tangentialElectricField),std::abs(et),1,"common tangential field analytic error");}}
    auto fresh=ElectricVofSharpField::prepareElectrostatic(mesh,geometry,charge,problem);
    bool rejected=false;try{stress.validateFor(mesh,fresh);}catch(const std::invalid_argument&){rejected=true;}
    check(rejected,"same-value fresh field owner was accepted");++scenarios;
}
void obliqueNearZeroTangential()
{
    namespace Cut=LAMP::Geometry::PlanarCellCut;
    DictionaryMesh config;config.Nx=config.Ny=4;config.maxCellNum=256;config.maxInitialIterations=0;
    QuadtreeMesh mesh(config);FaceSegmentGraph graph;graph.build(mesh);
    auto domain=NativeLeafDomain::build(mesh,graph,false);
    const Vector2D normal{.6,.8},tangent{-.8,.6};const double height=.84375,enLiquid=.4,enGas=.7;
    const long double norm=std::hypot(static_cast<long double>(normal.x),static_cast<long double>(normal.y));
    std::vector<std::array<long double,2>> stocks(mesh.cellIndexNum());
    for(int cell=0;cell<mesh.cellIndexNum();++cell){const auto& c=mesh.cells(mesh.cellIndexToLabel[cell]);
        const auto lo=mesh.verts(c.cellToVert[0]).coord(),hi=mesh.verts(c.cellToVert[2]).coord();
        const long double cx=(static_cast<long double>(lo.x)+hi.x)/2,cy=(static_cast<long double>(lo.y)+hi.y)/2;
        const long double offset=-(normal.x*cx+normal.y*(cy-height))/norm;
        // Independent complementary plane integrals; neither stock is capacity minus the other.
        stocks[cell][0]=Cut::volume({lo,hi},normal,offset,false);
        stocks[cell][1]=Cut::volume({lo,hi},{-normal.x,-normal.y},-offset,false);}
    auto publication=GeometricVofPublication::prepareLocalPlicInventory(mesh,std::move(domain),
        {0,0,0,LAMP::TimeLayer::Position::Initial},stocks,{});
    VofSmoothGeometryConfig gc;gc.allowBoundaryEndpoints=true;
    const auto geometry=VofSmoothGeometry::prepare(mesh,publication,gc);
    check(!geometry.arcs().empty(),"oblique fixture has no actual interface");
    ElectricVofSharpProblem problem;problem.epsilon={2,5};
    const double qs=problem.epsilon[1]*enGas-problem.epsilon[0]*enLiquid;
    const auto charge=ElectricVofSurfaceCharge::fromDensity(mesh,geometry,[=](std::size_t,double){return qs;});
    const auto before=charge.moments();
    for(double et:{0.,1e-12}){
        for(auto& b:problem.boundaries)b.value=[=](Vector2D p){
            const double d=normal.x*p.x+normal.y*(p.y-height),s=tangent.x*p.x+tangent.y*(p.y-height);
            return -(d<=0?enLiquid:enGas)*d-et*s;};
        const auto field=ElectricVofSharpField::prepareElectrostatic(mesh,geometry,charge,problem);
        const auto stress=ElectricVofSharpMaxwellStress::prepare(mesh,field);stress.validateFor(mesh,field);
        check(field.geometry().sameInstance(geometry)&&field.sourceCharge().sameInstance(charge),"oblique stress lost actual field/charge owner");
        check(field.report().usedSuperLU&&field.report().maximumRelativeGaussResidual<=1e-10,"oblique field lost actual solve/Gauss gate");
        double maximumEtError=0,maximumTangentialResidual=0;
        for(std::size_t arc=0;arc<geometry.arcs().size();++arc)for(double t:{.13,.5,.87}){
            const auto point=geometry.arcs()[arc].sample(t);
            near(static_cast<double>(normal.x*point.x+normal.y*(point.y-height)),0,1,"oblique G1 is not the prescribed plane");
            const auto value=stress.interfaceSample(arc,t);const auto& e=value.field;
            near(e.normal.x,normal.x,1,"oblique actual normal x");near(e.normal.y,normal.y,1,"oblique actual normal y");
            near(e.tangent.x,tangent.x,1,"oblique actual tangent x");near(e.tangent.y,tangent.y,1,"oblique actual tangent y");
            near(e.normalElectricField[0],enLiquid,1,"oblique liquid normal field");
            near(e.normalElectricField[1],enGas,1,"oblique gas normal field");
            maximumEtError=std::max(maximumEtError,std::abs(e.tangentialElectricField-et));
            check(std::abs(e.tangentialElectricField-et)<=1e-10,"oblique actual trace is not near the prescribed zero tangential field");
            near(problem.epsilon[1]*e.normalElectricField[1]-problem.epsilon[0]*e.normalElectricField[0],qs,1,"oblique independent charge jump");
            long double contracted=0,budget=std::abs(static_cast<long double>(e.surfaceChargeDensity)*e.tangentialElectricField);
            for(unsigned phase=0;phase<2;++phase){const double en=phase?enGas:enLiquid,eps=problem.epsilon[phase];
                const double ex=en*normal.x+et*tangent.x,ey=en*normal.y+et*tangent.y;
                const auto& a=value.tensor[phase];const double scale=eps*(ex*ex+ey*ey);
                near(a.xx,.5*eps*(ex*ex-ey*ey),scale,"oblique analytic Txx");
                near(a.xy,eps*ex*ey,scale,"oblique analytic Txy");
                near(a.yy,.5*eps*(ey*ey-ex*ex),scale,"oblique analytic Tyy");
                // Independent long-double contraction of the actual returned tensor.
                const long double pxx=static_cast<long double>(e.tangent.x)*a.xx*e.normal.x;
                const long double pxy=static_cast<long double>(e.tangent.x)*a.xy*e.normal.y;
                const long double pyx=static_cast<long double>(e.tangent.y)*a.xy*e.normal.x;
                const long double pyy=static_cast<long double>(e.tangent.y)*a.yy*e.normal.y;
                contracted+=(phase?1:-1)*((pxx+pxy)+(pyx+pyy));
                budget+=std::abs(pxx)+std::abs(pxy)+std::abs(pyx)+std::abs(pyy)
                    +std::abs(static_cast<long double>(eps)*e.normalElectricField[phase]*e.tangentialElectricField);}
            const long double expected=static_cast<long double>(e.surfaceChargeDensity)*e.tangentialElectricField;
            check(budget>0&&std::abs(contracted-expected)<=1e-10L*budget,"oblique full physical traction budget");
            near(value.tangentialJump,qs*et,1,"oblique analytic tangential traction");
            near(value.normalJump,.5*(problem.epsilon[1]*enGas*enGas-problem.epsilon[0]*enLiquid*enLiquid
                -(problem.epsilon[1]-problem.epsilon[0])*et*et),1,"oblique analytic normal traction");
            check(value.relativeNormalIdentityResidual<=1e-12&&value.relativeTangentialJumpResidual<=1e-10,"oblique unchanged stress residual gates");
            maximumTangentialResidual=std::max(maximumTangentialResidual,value.relativeTangentialJumpResidual);}
        check(charge.moments()==before,"oblique read-only field/stress changed accepted charge");
        std::cout<<"oblique Et="<<et<<" maximum_Et_error="<<maximumEtError
            <<" tangential_residual="<<maximumTangentialResidual<<'\n';++scenarios;
    }
}
}
int main()
{
    try{scenario(false,false);scenario(false,true);scenario(true,false);obliqueNearZeroTangential();
        std::cout<<"VOF_SHARP_MAXWELL_PASS scenarios="<<scenarios<<" checks="<<checks<<'\n';return 0;}
    catch(const std::exception& e){std::cerr<<"VOF sharp Maxwell FAIL: "<<e.what()<<'\n';return 1;}
}
