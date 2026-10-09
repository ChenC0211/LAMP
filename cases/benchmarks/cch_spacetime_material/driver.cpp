// H08a research consumer. Analytic references never enter production geometry
// or material integrals. Long windows are deliberately outside default CTest.
#include "Modules/Hydrodynamics/CchDualMaterial.h"
#include "interface_first_jet_fixture.h"
#include <filesystem>
#include <fstream>

namespace {
using namespace Quantum::Hydrodynamics;
CchPhysicalDivergence metrics(const QuadtreeMesh& mesh,bool rz)
{FaceSegmentGraph graph;graph.build(mesh);return CchPhysicalDivergence::build(mesh,CchIntegralDivergence::build(mesh,CchDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,rz))));}
CchNativeValues velocity(const QuadtreeMesh& mesh,const CchDofLayout& l,bool rz)
{std::vector<double> values;for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(mesh,l.nativeFace(u));values.push_back(f.horizontal?(rz?.2:0):(rz?0:.2));}return {l,values,std::vector<double>(l.pressureCount())};}
long double exactAtom(const CchAdvectorAtom& a,const Case& s,double begin,double end)
{
    if(a.normalAxis!=(s.radial?1:0))return 0;
    // Integrating a translating cross-section in time is exactly a swept
    // Cartesian/RZ box volume. No time quadrature or production phase roots.
    Vector2D lo=a.begin,hi=a.end;
    if(s.radial){lo.y=a.begin.y-s.speed*end;hi.y=a.begin.y-s.speed*begin;}
    else {lo.x=a.begin.x-s.speed*end;hi.x=a.begin.x-s.speed*begin;}
    return circleCell(lo,hi,{s.radial?0:s.cx,s.cy},s.radius,s.radial);
}
struct Observation {long double flux=0,midpointFlux=0,defect=0,geometry=0,dualX=0,dualY=0;};
Observation window(std::ofstream& rows,std::ofstream& atoms,std::ofstream& cells,bool rz,int amr,int n,int steps,bool polynomial)
{
    auto mesh=makeMesh(n,amr);Case shape;shape.pureFirstJet=true;shape.radial=rz;shape.motion=1;shape.nonpolynomial=!polynomial;
    const double endTime=.04,rho=1e6;auto original=initial(*mesh,shape);
    const auto metric=metrics(*mesh,rz);const auto g=CchMomentumGeometry::build(*mesh,metric);
    const auto r=CchAdvectorReconstruction::prepare(*mesh,g);const auto search=CchCharacteristicGeometry::prepare(*mesh,r);
    const auto values=velocity(*mesh,g.layout(),rz);Observation out;
    CchSpaceTimeConfig config;config.characteristic.maxStep=.01;
    const HermiteGalsBoundaryFirstJetQuery boundary=[shape](Vector2D p,HermiteGalsTimeLayer t,Boundary){const auto j=shape.profile(t.time())(p);return HermiteGalsBoundaryFirstJet{j.value,{j.dx,j.dy}};};
    std::size_t samples=0,panels=0;double estimated=0;
    for(int step=0;step<steps;++step){
        const double begin=endTime*step/steps,end=endTime*(step+1)/steps,dt=end-begin;
        const auto a=r.evaluate(*mesh,values,original.state().identity(),original.state().timeLayer());
        const auto b=r.evaluate(*mesh,values,original.state().identity(),{step+1,begin,end,Position::End});
        const auto source=search.linear(*mesh,a,b);
        const auto material=CchMaterialTransport::prepareSpaceTime(*mesh,original,source,rho,1,config,boundary);
        const auto mass=[&](const CchPhaseGeometry& p){return CchMomentumMass::build(*mesh,g,p,rho,1);};
        const auto dual=CchDualMaterial::prepare(*mesh,material,mass(material.beginGeometry()),mass(material.midpointRate().geometry()),mass(material.endGeometry()));
        const auto& report=material.spaceTimeReport();samples+=report.geometrySamples;panels+=report.acceptedPanels;estimated=std::max(estimated,static_cast<double>(report.maxRelativePhaseErrorEstimate));
        long double fluxError=0,midpointError=0,work=0,geometryError=0;
        for(std::size_t i=0;i<r.atoms().size();++i){const auto& atom=r.atoms()[i];const auto exact=exactAtom(atom,shape,begin,end);
            const auto actual=material.phaseAtomIntegrals()[i].liquid,mid=dt*material.midpointRate().atomRates()[i].liquid;
            if(atom.primalVelocity>=0){fluxError+=std::abs(actual-exact);midpointError+=std::abs(mid-exact);work+=atom.area*.2*dt;}
            atoms<<rz<<','<<amr<<','<<n<<','<<steps<<','<<polynomial<<','<<step<<','<<i<<','<<atom.primalVelocity<<','<<actual<<','<<exact<<','<<mid<<'\n';
        }
        for(int p=0;p<g.layout().pressureCount();++p){const auto& cell=mesh->cells(mesh->cellIndexToLabel[g.layout().nativeCell(p)]);
            const auto lo=mesh->verts(cell.cellToVert[0]).coord(),hi=mesh->verts(cell.cellToVert[2]).coord();
            const Vector2D center{rz?0:shape.cx+shape.speed*end,shape.cy+(rz?shape.speed*end:0)};
            const auto exact=circleCell(lo,hi,center,shape.radius,rz),actual=material.endGeometry().cellVolumes()[p].liquid;
            geometryError+=std::abs(actual-exact);const auto& c=material.cells()[p];
            cells<<rz<<','<<amr<<','<<n<<','<<steps<<','<<polynomial<<','<<step<<','<<p<<','<<actual<<','<<exact<<','<<c.massChange<<','<<c.outwardIntegral<<','<<c.integratedDefect<<','<<c.defectRate<<'\n';
        }
        out.flux+=fluxError;out.midpointFlux+=midpointError;out.defect+=material.report().l1IntegratedDefect;
        out.geometry=geometryError;out.dualX+=dual.report().components[0].l1IntegratedDefect;out.dualY+=dual.report().components[1].l1IntegratedDefect;
        rows<<rz<<','<<amr<<','<<n<<','<<steps<<','<<polynomial<<','<<step<<','<<dt<<','<<report.geometrySamples<<','<<report.acceptedPanels<<','<<report.maxRelativePhaseErrorEstimate
            <<','<<fluxError<<','<<midpointError<<','<<work<<','<<material.report().l1IntegratedDefect<<','<<material.report().l1IntegratedDefect/dt
            <<','<<dual.report().components[0].l1IntegratedDefect/dt<<','<<dual.report().components[1].l1IntegratedDefect/dt<<','<<geometryError<<std::endl;
        original=material.endGeometry().publication();
    }
    // Physical transport scale is rho*speed*projected interface area, fixed
    // under h/dt refinement. It does not grow with the number of mesh faces.
    const long double physicalRate=rho*.2*(rz?pi*shape.radius*shape.radius:2*shape.radius);
    std::cout<<"SPACETIME_WINDOW,"<<rz<<','<<amr<<','<<n<<','<<steps<<','<<polynomial<<','<<samples<<','<<panels<<','<<estimated
        <<','<<out.flux<<','<<out.midpointFlux<<','<<out.defect/endTime<<','<<out.defect/endTime/physicalRate
        <<','<<out.dualX/endTime/physicalRate<<','<<out.dualY/endTime/physicalRate<<','<<out.geometry<<std::endl;
    return out;
}
}
int main(int argc,char** argv)
{
    try{
        const std::string mode=argc==1?"full":argv[1];check(argc<=2&&(mode=="full"||mode=="probe"||mode=="space"||mode=="time"),"usage: cch spacetime benchmark [full|probe|space|time]");
        std::filesystem::create_directories("output/series");std::ofstream rows("output/series/intervals.csv"),atoms("output/series/atoms.csv"),cells("output/series/cells.csv");
        check(rows&&atoms&&cells,"cannot create research output");std::cout<<std::setprecision(17);rows<<std::setprecision(17);atoms<<std::setprecision(17);cells<<std::setprecision(17);
        rows<<"rz,amr,n,steps,polynomial,step,dt,samples,panels,estimated_relative_error,flux_l1,midpoint_flux_l1,face_work,R_integrated_l1,R_rate_l1,dual_x_R_rate,dual_y_R_rate,geometry_l1\n";
        atoms<<"rz,amr,n,steps,polynomial,step,atom,primal_face,liquid_integral,analytic_integral,midpoint_integral\n";
        cells<<"rz,amr,n,steps,polynomial,step,cell,liquid_volume,analytic_volume,mass_change,outward_integral,R_integral,R_rate\n";
        if(mode=="probe"){for(bool rz:{false,true})window(rows,atoms,cells,rz,0,4,1,false);}
        else {
            if(mode!="time")for(bool rz:{false,true})for(int amr:{0,1})for(int n:{4,8,16})window(rows,atoms,cells,rz,amr,n,1,false);
            if(mode!="space")for(bool rz:{false,true})for(int amr:{0,1})for(int steps:{1,2,4,8})window(rows,atoms,cells,rz,amr,8,steps,false);
        }
        // Successful execution is data completion only, not H08a acceptance.
        std::cout<<"SPACETIME_DATA_COMPLETE,"<<mode<<std::endl;return 0;
    }catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}
}
