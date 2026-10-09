// H12 research driver: actual constant load, physical mass/viscosity and saddle.
// Fixed geometry only; this is not a coupled time integrator.
#include "Modules/Hydrodynamics/CchInterfaceNormalLoad.h"
#include "Modules/Hydrodynamics/CchMomentumGeometry.h"
#include "Modules/Hydrodynamics/CchPhaseViscosity.h"
#include "Modules/Hydrodynamics/CchSinglePhaseSaddle.h"
#include "interface_first_jet_fixture.h"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <map>

namespace {
using namespace Quantum::Hydrodynamics;
using Matrix=CchSaddleVelocityBlock::Matrix;
using Clock=std::chrono::steady_clock;
double elapsed(Clock::time_point a){return std::chrono::duration<double>(Clock::now()-a).count();}
double relative(double e,double scale){check(std::isfinite(e)&&std::isfinite(scale),"nonfinite audit value");if(scale==0){check(e==0,"nonzero unscaled error");return 0;}return std::abs(e)/scale;}
std::ofstream open(const char* name,const char* header){std::ofstream f(std::string("output/series/")+name);check(bool(f),"cannot open series");f<<std::setprecision(17)<<header<<'\n';return f;}
std::unique_ptr<QuadtreeMesh> studyMesh(int n,bool amr,bool rz){
    DictionaryMesh d;d.Nx=d.Ny=n;d.maxCellNum=32768;d.maxInitialIterations=0;auto m=std::make_unique<QuadtreeMesh>(d);
    if(amr){DictionaryMesh::RefineBox b;b.level=1;b.min={rz?0:.25,.25};b.max={rz?.5:.75,.75};m->applyRefineBoxes({b});}return m;
}
}
int main(int argc,char** argv){try{
    const bool probe=argc==2&&std::string(argv[1])=="probe";check(argc<=2&&(argc==1||probe||std::string(argv[1])=="full"),"usage: constant stress benchmark [full|probe]");
    std::filesystem::create_directories("output/series");
    auto geometries=open("geometries.csv","geometry,rz,amr,n,center,radius,cells,faces,cut_cells,coarse_faces,geometry_seconds");
    auto cells=open("cells.csv","geometry,p,native,level,x0,y0,x1,y1,volume,liquid,gas,raw_liquid,alpha");
    auto faces=open("faces.csv","geometry,u,native,horizontal,x0,y0,x1,y1,minus,plus,area,free");
    auto masses=open("masses.csv","geometry,rho_l,u,dual_liquid,dual_gas,mass");
    auto viscous=open("viscosity.csv","geometry,mu_l,row,col,value");
    auto fields=open("pressure.csv","solve,p,value");
    auto velocity=open("velocity.csv","solve,u,value,force,pressure_load,inertia,reaction,diagonal");
    auto summaries=open("solves.csv","solve,geometry,rho_l,mu_l,jump,dt,load_error,velocity_error,pressure_error,balance_error,algebraic,momentum,continuity,mean,refinements,prepare_seconds,solve_seconds,status");
    auto failures=open("failures.csv","solve,geometry,rho_l,mu_l,jump,reason");
    int geometryId=0,solveId=0,failed=0,blocks=0;const auto start=Clock::now();
    const std::vector<int> resolutions=probe?std::vector<int>{8}:std::vector<int>{8,16,32};
    const std::vector<double> positions=probe?std::vector<double>{.4}:std::vector<double>{.4,.44,.48};
    for(bool rz:{false,true})for(bool amr:{false,true})for(int n:resolutions)for(double center:positions){
        const auto gs=Clock::now();++geometryId;auto mesh=studyMesh(n,amr,rz);Case shape;shape.pureFirstJet=true;shape.radial=rz;shape.radius=.19;shape.cx=rz?0:center;shape.cy=rz?center:.5;
        FaceSegmentGraph graph;graph.build(*mesh);const auto l=CchDofLayout::build(*mesh,NativeLeafDomain::build(*mesh,graph,rz));
        const auto metrics=CchPhysicalDivergence::build(*mesh,CchIntegralDivergence::build(*mesh,l));
        const auto phase=CchPhaseGeometry::build(*mesh,metrics,initial(*mesh,shape));check(phase.publication().state().primary().phaseVolumes.empty(),"unexpected Interface inventory");
        const auto pressure=CchPressureLoad::build(metrics.divergence());const auto dual=CchMomentumGeometry::build(*mesh,metrics);
        const auto quadrature=HermiteGalsPhaseQuadrature::build(*mesh,phase.publication());
        std::map<int,long double> raw;for(const auto& c:phase.publication().volumes().cells())raw[c.cell]=c.rawLiquid;
        double sumLiquid=0,sumVolume=0;for(int p=0;p<l.pressureCount();++p){const auto& c=mesh->cells(mesh->cellIndexToLabel[l.nativeCell(p)]);const auto a=mesh->verts(c.cellToVert[0]).coord(),b=mesh->verts(c.cellToVert[2]).coord();
            cells<<geometryId<<','<<p<<','<<l.nativeCell(p)<<','<<c.level()<<','<<a.x<<','<<a.y<<','<<b.x<<','<<b.y<<','<<metrics.cellVolumes()[p]<<','<<phase.cellVolumes()[p].liquid<<','<<phase.cellVolumes()[p].gas<<','<<raw.at(l.nativeCell(p))<<','<<phase.liquidFractions()[p]<<'\n';
            sumLiquid+=metrics.cellVolumes()[p]*phase.liquidFractions()[p];sumVolume+=metrics.cellVolumes()[p];}
        std::vector<int> minus,plus;std::vector<double> area;int coarse=0;
        for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(*mesh,l.nativeFace(u));const auto& e=f.horizontal?mesh->hEdges(f.edgeLabel):mesh->vEdges(f.edgeLabel);const auto a=mesh->verts(e.edgeToVert[0]).coord(),b=mesh->verts(e.edgeToVert[1]).coord();
            minus.push_back(f.negativeCell<0?-1:l.pressureDof(f.negativeCell));plus.push_back(f.positiveCell<0?-1:l.pressureDof(f.positiveCell));
            area.push_back((f.horizontal?std::abs(b.x-a.x):std::abs(b.y-a.y))*(rz?double(2*pi)*(f.horizontal?(a.x+b.x)/2:a.x):1));
            const bool free=minus.back()>=0&&plus.back()>=0&&area.back()>0;
            if(free&&mesh->cells(mesh->cellIndexToLabel[f.negativeCell]).level()!=mesh->cells(mesh->cellIndexToLabel[f.positiveCell]).level())++coarse;
            faces<<geometryId<<','<<u<<','<<l.nativeFace(u)<<','<<f.horizontal<<','<<a.x<<','<<a.y<<','<<b.x<<','<<b.y<<','<<minus.back()<<','<<plus.back()<<','<<metrics.faceMeasures()[u]<<','<<free<<'\n';}
        check(!amr||coarse>0,"missing actual coarse/fine faces");
        geometries<<geometryId<<','<<rz<<','<<amr<<','<<n<<','<<center<<",0.19,"<<l.pressureCount()<<','<<l.velocityCount()<<','<<phase.report().cutCells<<','<<coarse<<','<<elapsed(gs)<<'\n';geometries.flush();
        std::vector<CchPhaseViscosity> viscosities;for(double mu:{1.,1e-3,1e-6,0.}){viscosities.push_back(CchPhaseViscosity::build(*mesh,phase,quadrature,mu,.01*mu));const auto& k=viscosities.back().viscosity().matrix();
            check(mu!=0||k.nonZeros()==0,"nonzero inviscid operator");for(int row=0;row<k.rows();++row)for(Matrix::InnerIterator it(k,row);it;++it)viscous<<geometryId<<','<<mu<<','<<row<<','<<it.col()<<','<<it.value()<<'\n';}
        for(double density:{1.,1e3,1e6}){const auto mass=CchMomentumMass::build(*mesh,dual,phase,density,1);
            for(int u=0;u<l.velocityCount();++u)masses<<geometryId<<','<<density<<','<<u<<','<<mass.phaseVolumes()[u].liquid<<','<<mass.phaseVolumes()[u].gas<<','<<mass.masses()[u]<<'\n';
            for(const auto& viscosity:viscosities){const double mu=viscosity.liquidViscosity();const auto ps=Clock::now();Matrix a=viscosity.viscosity().matrix();for(int u=0;u<l.velocityCount();++u)a.coeffRef(u,u)+=mass.masses()[u]/.125;
                const auto saddle=CchSinglePhaseSaddle::buildClosed(pressure,CchSaddleVelocityBlock::build(l,a));check(saddle.pressureGaugeCount()==1,"full domain disconnected");++blocks;const double prepare=elapsed(ps);
                for(double jump:{-1.75,2.5}){++solveId;const auto ss=Clock::now();std::cout<<"H12_BEGIN,"<<solveId<<','<<geometryId<<','<<density<<','<<mu<<','<<jump<<std::endl;
                    try{const auto normal=CchInterfaceNormalLoad::constant(*mesh,pressure,phase,jump);auto rhs=normal.balanceLoads();for(double& f:rhs)f=-f;
                        const auto result=saddle.solve(*mesh,CchSaddleRhs(l,rhs,std::vector<double>(l.velocityCount())));const double duration=elapsed(ss);
                        const auto& u=result.values().normalVelocity();const auto& p=result.values().pressure();const auto& force=normal.balanceLoads();
                        const auto cp=pressure.apply(*mesh,result.values());const Eigen::Map<const Eigen::VectorXd> ev(u.data(),u.size());const Eigen::VectorXd inertia=a*ev;
                        double eLoad=0,sLoad=0,eVelocity=0,sVelocity=0,ePressure=0,eBalance=0,sBalance=0;
                        for(int k=0;k<l.pressureCount();++k){fields<<solveId<<','<<k<<','<<p[k]<<'\n';ePressure=std::max(ePressure,std::abs(p[k]-jump*(sumLiquid/sumVolume-phase.liquidFractions()[k])));}
                        for(int k=0;k<l.velocityCount();++k){const double exact=jump*area[k]*((plus[k]<0?0:phase.liquidFractions()[plus[k]])-(minus[k]<0?0:phase.liquidFractions()[minus[k]]));
                            eLoad=std::hypot(eLoad,force[k]-exact);sLoad=std::hypot(sLoad,exact);eVelocity=std::max(eVelocity,std::abs(u[k]));sVelocity=std::max(sVelocity,std::abs(force[k])/a.coeff(k,k));
                            const bool free=minus[k]>=0&&plus[k]>=0&&area[k]>0;const double reaction=inertia[k]+cp[k]+force[k];
                            if(free){eBalance=std::hypot(eBalance,cp[k]+force[k]);sBalance=std::hypot(sBalance,std::abs(cp[k])+std::abs(force[k]));}else check(u[k]==0,"nonzero fixed boundary");
                            velocity<<solveId<<','<<k<<','<<u[k]<<','<<force[k]<<','<<cp[k]<<','<<inertia[k]<<','<<reaction<<','<<a.coeff(k,k)<<'\n';}
                        const double le=relative(eLoad,sLoad),ve=relative(eVelocity,sVelocity),pe=relative(ePressure,std::abs(jump)),be=relative(eBalance,sBalance);
                        const bool pass=le<=1e-12&&ve<=1e-11&&pe<=1e-10&&be<=1e-11;const auto& r=result.report();
                        summaries<<solveId<<','<<geometryId<<','<<density<<','<<mu<<','<<jump<<",0.125,"<<le<<','<<ve<<','<<pe<<','<<be<<','<<r.algebraicRelativeResidual<<','<<r.freeMomentumDefect<<','<<r.maxCellContinuityDefect<<','<<r.maxPressureMeanDefect<<','<<r.linearRefinements<<','<<prepare<<','<<duration<<','<<(pass?"PASS":"FAIL")<<'\n';summaries.flush();
                        if(!pass)++failed;std::cout<<"H12_RESULT,"<<solveId<<','<<ve<<','<<pe<<','<<r.maxCellContinuityDefect<<','<<r.linearRefinements<<','<<(pass?"PASS":"FAIL")<<std::endl;
                    }catch(const std::exception& e){++failed;failures<<solveId<<','<<geometryId<<','<<density<<','<<mu<<','<<jump<<",\""<<e.what()<<"\"\n";failures.flush();std::cout<<"H12_FAILED,"<<solveId<<','<<e.what()<<std::endl;}
                }
            }
        }
        for(auto* out:{&cells,&faces,&masses,&viscous,&fields,&velocity}){out->flush();check(bool(*out),"series write failure");}
    }
    check(geometryId==(probe?4:36)&&blocks==(probe?48:432)&&solveId==(probe?96:864),"incomplete parameter product");
    std::cout<<"H12_COMPLETE,"<<geometryId<<','<<blocks<<','<<solveId<<','<<failed<<','<<elapsed(start)<<std::endl;return failed?1:0;
}catch(const std::exception& e){std::cerr<<"H12 BENCHMARK FAIL: "<<e.what()<<'\n';return 1;}}
