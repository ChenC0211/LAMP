#include "Modules/Hydrodynamics/HydroBalancedCsfLoad.h"
#include "Modules/Hydrodynamics/HydroVofStep.h"
#include "Modules/Interface/GeometricVofInitialization.h"
#include "Modules/Interface/GeometricVofMeshTransfer.h"
#include "Modules/Interface/GeometricVofTransport.h"
#include "Core/Field/NativeTargetLevelField.h"
#include "Core/Mesh/NativeMeshIdentity.h"
#include "Core/Mesh/DyadicGridCoordinates.h"
#include "Core/Mesh/CartesianMeasure.h"
#include "Core/Geometry/ActiveQuarticHermiteCompletion.h"
#include "Core/Field/CellField.h"
#include "Core/Field/VertField.h"
#include "Core/Field/HEdgeField.h"
#include "Core/Field/VEdgeField.h"
#include "Modules/Hydrodynamics/HydroMomentumGeometry.h"
#include "Modules/Hydrodynamics/HydroPhaseViscosity.h"
#include "Modules/Hydrodynamics/HydroSinglePhaseSaddle.h"
#include "Modules/Hydrodynamics/HydroSpatialAdvector.h"
#include "Modules/Hydrodynamics/Detail/HydroRegularAdvectorGeometry.h"
#include "interface_first_jet_fixture.h"
#include <cstring> // F25 retained raw independent stock bytes.
#include <cstdint>
#include <map>
#include <numeric>
#include <chrono>
#include <Eigen/Cholesky>
#include <Eigen/LU>
#include <Eigen/Eigenvalues>

namespace
{
using namespace LAMP::Hydrodynamics;
using Matrix=HydroSaddleVelocityBlock::Matrix;
int scenes=0,solves=0,guards=0,negativeControls=0,coarseFaces=0,axisFaces=0,wallFaces=0;
double peakLoad=0,peakVelocity=0,peakPressure=0,peakBalance=0,peakWall=0,peakContinuity=0;
double relative(double error,double scale)
{
    check(std::isfinite(error)&&std::isfinite(scale)&&scale>=0,"invalid independent residual scale");
    if(scale==0){check(error==0,"nonzero error with exactly zero scale");return 0;}
    return std::abs(error)/scale;
}
template<class F>void rejection(F f,const char* expected)
{rejects(f,expected);++guards;}
HydroPhysicalDivergence metricsFor(const QuadtreeMesh& mesh,bool rz,int selection=0)
{
    FaceSegmentGraph graph;graph.build(mesh,[selection](const QuadtreeMesh&,const Cell& cell){const auto p=cell.centerCoord();
        return selection==0 || (selection==1?!(p.x>.375&&p.x<.625&&p.y>.375&&p.y<.625):
            selection==2?p.x<.375||p.x>.625:p.x>.25);});
    return HydroPhysicalDivergence::build(mesh,HydroIntegralDivergence::build(mesh,
        HydroDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,rz))));
}
HermiteGalsPublication published(const QuadtreeMesh& mesh,const Case& shape,double scale=1)
{
    const auto original=shape.profile(0);const Profile f=[=](Vector2D p){auto q=original(p);
        q.value*=scale;q.dx*=scale;q.dy*=scale;q.dxy*=scale;return q;};
    const auto state=stateFor(mesh,f,shape.radial,{},[](int,Vector2D,Vector2D){return 0.;},shape.symmetry(mesh),HermiteGalsTransport::FirstJet);
    HermiteGalsCandidate candidate(mesh,{},state.primary());const auto mixed=closeHermiteGalsMixedAmr(mesh,candidate,shape.symmetry(mesh));
    HermiteGalsPublicationConfig cfg;cfg.volume.symmetryPlanes=shape.symmetry(mesh);const auto id=candidate.identity();
    return publishHermiteGalsFirstJetGeometry(mesh,std::move(candidate),mixed,{id,centers(mesh,f)},cfg);
}
HermiteGalsPublication purePublication(const QuadtreeMesh& mesh,bool rz,double sign)
{
    const Profile f=[sign](Vector2D){return Jet{sign,0,0,0};};
    const auto state=stateFor(mesh,f,rz,{},[](int,Vector2D,Vector2D){return 0.;},{},HermiteGalsTransport::FirstJet);
    HermiteGalsCandidate candidate(mesh,{},state.primary());const auto mixed=closeHermiteGalsMixedAmr(mesh,candidate);
    const auto id=candidate.identity();return publishHermiteGalsFirstJetGeometry(mesh,std::move(candidate),mixed,{id,centers(mesh,f)});
}
HydroRegularizedColor forceColor(const QuadtreeMesh& mesh,const HydroPhaseGeometry& phase,double epsilon=.375)
{HydroRegularizedColorConfig cfg;cfg.epsilon=epsilon;return HydroRegularizedColor::build(mesh,phase,cfg);}

// Independent reference uses native incidence/physical edge coordinates,
// including every actual AMR subface. It never reads C/D coefficients.
struct Reference
{
    std::vector<int> minus,plus,component;std::vector<double> area,volume;int components=0;
    Reference(const QuadtreeMesh& mesh,const HydroPhaseGeometry& phase)
    {
        const auto& l=phase.layout();const int np=l.pressureCount();component.resize(np);std::iota(component.begin(),component.end(),0);
        const auto root=[&](int p){while(component[p]!=p)p=component[p];return p;};
        for(int p=0;p<np;++p){const auto& cell=mesh.cells(mesh.cellIndexToLabel[l.nativeCell(p)]);
            const auto lo=mesh.verts(cell.cellToVert[0]).coord(),hi=mesh.verts(cell.cellToVert[2]).coord();
            const double v=(hi.x-lo.x)*(hi.y-lo.y)*(l.domain().axisymmetric()?double(pi)*(lo.x+hi.x):1);
            volume.push_back(v);close(v,phase.metrics().cellVolumes()[p],1e-14*v,"true native cell volume");}
        for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(mesh,l.nativeFace(u));
            const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
            const auto lo=mesh.verts(edge.edgeToVert[0]).coord(),hi=mesh.verts(edge.edgeToVert[1]).coord();
            const double length=f.horizontal?std::abs(hi.x-lo.x):std::abs(hi.y-lo.y);
            area.push_back(length*(l.domain().axisymmetric()?double(2*pi)*(f.horizontal?(lo.x+hi.x)/2:lo.x):1));
            close(area.back(),phase.metrics().faceMeasures()[u],1e-14*area.back(),"true native subface area");
            minus.push_back(f.negativeCell<0?-1:l.pressureDof(f.negativeCell));plus.push_back(f.positiveCell<0?-1:l.pressureDof(f.positiveCell));
            if(minus.back()>=0&&plus.back()>=0){component[root(plus.back())]=root(minus.back());
                if(mesh.cells(mesh.cellIndexToLabel[f.negativeCell]).level()!=mesh.cells(mesh.cellIndexToLabel[f.positiveCell]).level())++coarseFaces;}
            if(area.back()==0)++axisFaces;
        }
        std::map<int,int> numbered;for(int p=0;p<np;++p){const int r=root(p);if(!numbered.count(r))numbered[r]=components++;}
        std::vector<int> final(np);for(int p=0;p<np;++p)final[p]=numbered[root(p)];component=std::move(final);
    }
    std::vector<double> gradient(const std::vector<double>& p)const
    {
        std::vector<double> result(area.size());for(std::size_t u=0;u<area.size();++u)
            result[u]=area[u]*((plus[u]<0?0:p[plus[u]])-(minus[u]<0?0:p[minus[u]]));return result;
    }
    std::vector<double> meanByCell(const std::vector<double>& c)const
    {
        std::vector<long double> sum(components),capacity(components);
        for(std::size_t p=0;p<c.size();++p){sum[component[p]]+=volume[p]*static_cast<long double>(c[p]);capacity[component[p]]+=volume[p];}
        std::vector<double> mean(c.size());for(std::size_t p=0;p<c.size();++p)mean[p]=static_cast<double>(sum[component[p]]/capacity[component[p]]);return mean;
    }
};
void algebra(const QuadtreeMesh& mesh,const HydroPressureLoad& pressure,const HydroRegularizedColor& color,
    const Reference& ref,double sigma,double kappa)
{
    const auto load=HydroBalancedCsfLoad::constantCurvature(mesh,pressure,color,sigma,kappa);
    load.validateFor(mesh,pressure,color);const double beta=sigma*kappa;
    check(load.surfaceTension()==sigma&&load.constantCurvature()==kappa,"lost explicit constant stress parameters");
    check(load.geometryIdentity()==color.geometryIdentity(),"CSF geometry identity changed");
    const auto gradient=ref.gradient(color.colors());double error=0,scale=0;
    for(std::size_t u=0;u<gradient.size();++u){error=std::hypot(error,load.colorGradient()[u]-gradient[u]);scale=std::hypot(scale,gradient[u]);}
    peakLoad=std::max(peakLoad,relative(error,scale));check(relative(error,scale)<=1e-12,"independent all-row color gradient");
    error=scale=0;for(std::size_t u=0;u<gradient.size();++u){const double expected=beta*gradient[u];
        error=std::hypot(error,load.momentumRhsLoad()[u]-expected);scale=std::hypot(scale,expected);
        if(sigma==0||kappa==0)check(load.momentumRhsLoad()[u]==0&&!std::signbit(load.momentumRhsLoad()[u]),"zero sigma/kappa not exact positive zero");
        if(ref.area[u]==0)check(load.colorGradient()[u]==0&&load.momentumRhsLoad()[u]==0,"axis load nonzero");}
    peakLoad=std::max(peakLoad,relative(error,scale));check(relative(error,scale)<=1e-12,"independent CSF RHS sign/area");
    auto potential=color.colors();for(double& c:potential)c*=beta;
    const auto cp=pressure.apply(mesh,HydroNativeValues(color.phaseGeometry().layout(),std::vector<double>(gradient.size()),potential));
    error=0;for(std::size_t u=0;u<gradient.size();++u)error=std::hypot(error,cp[u]-load.momentumRhsLoad()[u]);
    check(relative(error,scale)<=1e-12,"all-row nongauge constant-curvature identity");
    // Virtual work uses independent physical cell/face incidence.
    std::vector<long double> div(color.colors().size());long double work=0,workScale=0;
    for(std::size_t u=0;u<gradient.size();++u){const double velocity=std::sin(.31*(u+1));const auto term=velocity*load.momentumRhsLoad()[u];
        work+=term;workScale+=std::abs(term);if(ref.minus[u]>=0)div[ref.minus[u]]+=ref.area[u]*velocity;if(ref.plus[u]>=0)div[ref.plus[u]]-=ref.area[u]*velocity;}
    for(std::size_t p=0;p<potential.size();++p){const auto term=potential[p]*div[p];work+=term;workScale+=std::abs(term);}
    check(relative(static_cast<double>(work),static_cast<double>(workScale))<=1e-12,"C=-D_I transpose virtual work");
}
void solveStatic(const QuadtreeMesh& mesh,const HydroPhaseGeometry& phase,const HydroRegularizedColor& color,
    const HydroPressureLoad& pressure,const HydroMomentumGeometry& dual,const Matrix& viscosity,
    const Reference& ref,const char* tag,int grid,int selection,double rho,double mu,bool wrongSign=false,double positiveCurvature=5)
{
    ++scenes;const auto& l=phase.layout();const auto mass=HydroMomentumMass::build(mesh,dual,phase,rho,1);
    Matrix a=viscosity;if(mu==0)check(a.nonZeros()==0,"zero viscosity has artificial floor");
    for(int u=0;u<l.velocityCount();++u)a.coeffRef(u,u)+=mass.masses()[u]/.125;
    const auto saddle=HydroSinglePhaseSaddle::buildClosed(pressure,HydroSaddleVelocityBlock::build(l,a));
    check(saddle.pressureGaugeCount()==ref.components,"independent connected-component gauge");
    std::vector<bool> free(l.velocityCount(),false);for(int u:saddle.freeVelocityDofs())free[u]=true;
    for(int u=0;u<l.velocityCount();++u)check(free[u]==(ref.minus[u]>=0&&ref.plus[u]>=0&&ref.area[u]>0),"native free-row classification");
    const auto means=ref.meanByCell(color.colors()),meanGradient=ref.gradient(means);
    for(const auto parameters:std::array<std::array<double,2>,3>{{{.7,-2.5},{2.5,0},{.5,positiveCurvature}}}){
        const double sigma=parameters[0],kappa=parameters[1],beta=sigma*kappa;
        algebra(mesh,pressure,color,ref,sigma,kappa);
        const auto load=HydroBalancedCsfLoad::constantCurvature(mesh,pressure,color,sigma,kappa);
        const auto result=saddle.solve(mesh,HydroSaddleRhs(l,load.momentumRhsLoad(),std::vector<double>(l.velocityCount())));++solves;
        std::vector<double> expected(l.pressureCount()),au(l.velocityCount());
        for(int p=0;p<l.pressureCount();++p)expected[p]=beta*(color.colors()[p]-means[p]);
        for(int column=0;column<a.outerSize();++column)for(Matrix::InnerIterator it(a,column);it;++it)
            au[it.row()]+=it.value()*result.values().normalVelocity()[it.col()];
        const auto cp=ref.gradient(result.values().pressure());
        double pressureError=0,velocityMax=0,velocityScale=0,balanceError=0,balanceScale=0,wallError=0,wallScale=0;
        for(int p=0;p<l.pressureCount();++p){pressureError=std::max(pressureError,std::abs(result.values().pressure()[p]-expected[p]));
            std::cout<<"BF02_CELL,"<<scenes<<','<<tag<<','<<grid<<','<<selection<<','<<rho<<','<<mu<<','<<sigma<<','<<kappa<<','<<p<<','
                <<l.nativeCell(p)<<','<<ref.component[p]<<','<<ref.volume[p]<<','<<color.config().epsilon<<','<<color.colors()[p]<<','
                <<phase.liquidFractions()[p]<<','<<result.values().pressure()[p]<<','<<expected[p]<<','<<means[p]<<'\n';}
        for(int u=0;u<l.velocityCount();++u){const double reaction=au[u]+cp[u]-load.momentumRhsLoad()[u];
            const double expectedReaction=-beta*meanGradient[u];
            velocityMax=std::max(velocityMax,std::abs(result.values().normalVelocity()[u]));
            velocityScale=std::max(velocityScale,std::abs(load.momentumRhsLoad()[u])/a.coeff(u,u));
            if(free[u]){balanceError=std::hypot(balanceError,cp[u]-load.momentumRhsLoad()[u]);
                balanceScale=std::hypot(balanceScale,std::abs(cp[u])+std::abs(load.momentumRhsLoad()[u]));}
            else{check(result.values().normalVelocity()[u]==0,"fixed-row velocity not exact zero");
                wallError=std::hypot(wallError,reaction-expectedReaction);wallScale=std::hypot(wallScale,std::abs(reaction)+std::abs(expectedReaction));++wallFaces;}
            if(ref.area[u]==0){check(load.momentumRhsLoad()[u]==0&&cp[u]==0,"zero-area axis CSF/pressure load not exact zero");
                // Prescribed axis u is zero. Finite viscosity couples its row
                // to roundoff-level free velocities; Au is an actual reaction,
                // assessed with every wall row at the unchanged wall gate.
                if(mu==0)check(au[u]==0,"inviscid axis inertia not exact zero");}
            std::cout<<"BF02_FACE,"<<scenes<<','<<tag<<','<<grid<<','<<selection<<','<<rho<<','<<mu<<','<<sigma<<','<<kappa<<','<<u<<','
                <<ref.minus[u]<<','<<ref.plus[u]<<','<<ref.area[u]<<','<<free[u]<<','<<(ref.area[u]==0)<<','<<load.colorGradient()[u]<<','
                <<load.momentumRhsLoad()[u]<<','<<cp[u]<<','<<(cp[u]-load.momentumRhsLoad()[u])<<','<<au[u]<<','<<reaction<<','<<expectedReaction<<'\n';}
        const double velocity=relative(velocityMax,velocityScale),pError=relative(pressureError,std::abs(beta));
        const double balance=relative(balanceError,balanceScale),wall=relative(wallError,wallScale);
        check(velocity<=1e-11,"static velocity H12 gate");check(pError<=1e-10,"physical pressure and true-volume gauge H12 gate");
        check(balance<=1e-11,"free-row static balance H12 gate");check(wall<=1e-11,"physical fixed-row reaction");
        peakVelocity=std::max(peakVelocity,velocity);peakPressure=std::max(peakPressure,pError);peakBalance=std::max(peakBalance,balance);
        peakWall=std::max(peakWall,wall);peakContinuity=std::max(peakContinuity,result.report().maxCellContinuityDefect);
        if(beta>0){const auto high=std::max_element(color.colors().begin(),color.colors().end()),low=std::min_element(color.colors().begin(),color.colors().end());
            const int i=static_cast<int>(high-color.colors().begin()),j=static_cast<int>(low-color.colors().begin());
            if(ref.component[i]==ref.component[j]&&*high-*low>1e-6){check(result.values().pressure()[i]>result.values().pressure()[j],"convex liquid pressure sign reversed");
                check(std::abs(-result.values().pressure()[i]-expected[i])+std::abs(-result.values().pressure()[j]-expected[j])>1e-10*beta,
                    "pressure sign negative control invisible");++negativeControls;}
            if(wrongSign){auto badRhs=load.momentumRhsLoad();for(double& value:badRhs)value=-value;
                const auto wrong=saddle.solve(mesh,HydroSaddleRhs(l,badRhs,std::vector<double>(l.velocityCount())));++solves;
                double defect=0;for(int p=0;p<l.pressureCount();++p)defect=std::max(defect,std::abs(wrong.values().pressure()[p]-expected[p]));
                check(defect/std::abs(beta)>1e-10,"actual wrong-RHS sign negative control invisible");
                if(ref.component[i]==ref.component[j]&&*high-*low>1e-6)
                    check(wrong.values().pressure()[i]<wrong.values().pressure()[j],"actual wrong RHS failed to reverse liquid pressure");++negativeControls;}
        }
        std::cout<<"BF02_STATIC,"<<scenes<<','<<tag<<','<<rho<<','<<mu<<','<<sigma<<','<<kappa<<','<<velocity<<','<<pError<<','<<balance<<','<<wall<<','
            <<result.report().linearRefinements<<",PASS\n";
    }
}
void inspect(const QuadtreeMesh& mesh,const HydroPhaseGeometry& phase,const HydroRegularizedColor& color,
    const char* tag,int grid=0,int selection=0,bool densityMatrix=true,double positiveCurvature=5)
{
    check(phase.publication().state().primary().phaseVolumes.empty(),"CSF introduced conservative Interface inventory");
    const auto sharpBefore=phase.liquidFractions();const auto pressure=HydroPressureLoad::build(phase.metrics().divergence());
    const Reference ref(mesh,phase);const auto dual=HydroMomentumGeometry::build(mesh,phase.metrics());
    const auto quadrature=HermiteGalsPhaseQuadrature::build(mesh,phase.publication());
    algebra(mesh,pressure,color,ref,0,std::numeric_limits<double>::max());
    algebra(mesh,pressure,color,ref,std::numeric_limits<double>::max(),0);
    for(double mu:{0.,1e-3}){const auto viscosity=HydroPhaseViscosity::build(mesh,phase,quadrature,mu,.01*mu);
        for(double rho:{1.,1e3,1e6}){if(!densityMatrix&&rho!=1e3)continue;
            solveStatic(mesh,phase,color,pressure,dual,viscosity.viscosity().matrix(),ref,tag,grid,selection,rho,mu,
                densityMatrix&&grid==0&&rho==1&&mu==0,positiveCurvature);}}
    check(phase.liquidFractions()==sharpBefore,"CSF changed sharp-material fractions");
}
void core()
{
    for(bool rz:{false,true})for(int amr:{0,1,2}){auto mesh=makeMesh(4,amr);Case shape;shape.pureFirstJet=true;shape.radial=rz;shape.radius=rz?.34:.19;
        const auto phase=HydroPhaseGeometry::build(*mesh,metricsFor(*mesh,rz),published(*mesh,shape));const auto color=forceColor(*mesh,phase);
        // The constant-curvature drop control uses its independent analytic
        // physical curvature: 1/R in XY and 2/R for the RZ sphere.
        inspect(*mesh,phase,color,rz?"sphere":"circle",amr,0,true,(rz?2.:1.)/shape.radius);}
    std::uint32_t random=170120u;
    for(bool rz:{false,true})for(int amr:{0,1,2}){random=1664525u*random+1013904223u;
        const double offset=.31+.38*static_cast<double>(random)/std::numeric_limits<std::uint32_t>::max();
        auto mesh=makeMesh(4,amr);Case shape;shape.pureFirstJet=true;shape.radial=rz;shape.planar=true;shape.planeOffset=offset;
        const auto met=metricsFor(*mesh,rz);const auto phase=HydroPhaseGeometry::build(*mesh,met,published(*mesh,shape));const auto color=forceColor(*mesh,phase);
        const auto scaledPhase=HydroPhaseGeometry::build(*mesh,met,published(*mesh,shape,7));const auto scaled=forceColor(*mesh,scaledPhase);
        for(std::size_t p=0;p<color.colors().size();++p)close(color.colors()[p],scaled.colors()[p],5e-9,"positive phi scaling changed CSF color");
        const auto pressure=HydroPressureLoad::build(met.divergence());const auto original=HydroBalancedCsfLoad::constantCurvature(*mesh,pressure,color,.5,5);
        const auto altered=HydroBalancedCsfLoad::constantCurvature(*mesh,pressure,scaled,.5,5);
        const Reference ref(*mesh,phase);for(std::size_t u=0;u<ref.area.size();++u)
            close(original.momentumRhsLoad()[u],altered.momentumRhsLoad()[u],2.5*1e-8*ref.area[u],"scaled true-area CSF gradient");
        inspect(*mesh,phase,color,rz?"axial-plane":"plane",amr,0,false);
        std::cout<<"BF02_PLANE_OFFSET,"<<rz<<','<<amr<<','<<offset<<",PASS\n";}
}
void domains()
{
    for(bool rz:{false,true})for(int selection:{1,2}){auto mesh=makeMesh(8);Case shape;shape.pureFirstJet=true;shape.radial=rz;shape.planar=true;shape.planeOffset=.43;
        if(rz&&selection==2){shape.planar=false;shape.radius=.34;}
        const auto phase=HydroPhaseGeometry::build(*mesh,metricsFor(*mesh,rz,selection),published(*mesh,shape));
        const auto color=forceColor(*mesh,phase,.1875);inspect(*mesh,phase,color,selection==1?"hole":rz?"disconnected-sphere":"disconnected",0,selection,false,
            rz&&selection==2?2/shape.radius:5);
        const Reference ref(*mesh,phase);check(ref.components==(selection==1?1:2),"selected connected-component witness");
        if(selection==2){const auto mean=ref.meanByCell(color.colors());
            check(*std::max_element(mean.begin(),mean.end())-*std::min_element(mean.begin(),mean.end())>1e-4,"disconnected gauge lacked distinct component means");}}
    for(bool rz:{false,true})for(double sign:{-1.,1.}){auto mesh=makeMesh(4,1);
        const auto phase=HydroPhaseGeometry::build(*mesh,metricsFor(*mesh,rz),purePublication(*mesh,rz,sign));const auto color=forceColor(*mesh,phase);
        inspect(*mesh,phase,color,sign<0?"pure-liquid":"pure-gas",1,0,false);
        const auto load=HydroBalancedCsfLoad::constantCurvature(*mesh,HydroPressureLoad::build(phase.metrics().divergence()),color,.5,5);
        const Reference ref(*mesh,phase);for(std::size_t u=0;u<ref.area.size();++u)if(ref.minus[u]>=0&&ref.plus[u]>=0)
            check(load.colorGradient()[u]==0&&load.momentumRhsLoad()[u]==0,"pure internal/free CSF not exact zero");}
    auto mesh=makeMesh(4);Case shape;shape.pureFirstJet=true;shape.radial=true;shape.planar=true;shape.planeOffset=.43;
    const auto phase=HydroPhaseGeometry::build(*mesh,metricsFor(*mesh,true,3),published(*mesh,shape));
    inspect(*mesh,phase,forceColor(*mesh,phase),"annular",0,3,false);
}
void provenance()
{
    auto mesh=makeMesh(4);Case shape;shape.pureFirstJet=true;shape.planar=true;shape.planeOffset=.43;
    const auto phase=HydroPhaseGeometry::build(*mesh,metricsFor(*mesh,false),published(*mesh,shape));
    const auto pressure=HydroPressureLoad::build(phase.metrics().divergence());const auto color=forceColor(*mesh,phase);
    const auto load=HydroBalancedCsfLoad::constantCurvature(*mesh,pressure,color,.5,5),copy=load;
    check(&load.momentumRhsLoad()==&copy.momentumRhsLoad() && &load.colorGradient()==&copy.colorGradient(),"CSF copy duplicated immutable arrays");
    load.validateFor(*mesh,pressure,color);copy.validateFor(*mesh,pressure,color);
    const auto colorCopy=color;load.validateFor(*mesh,pressure,colorCopy);
    check(&load.color().colors()==&color.colors(),"CSF copied color data into a second owner");
    check(&load.color().phaseGeometry().publication().state().geometry().snapshot()==&phase.publication().state().geometry().snapshot(),"CSF rebuilt scalar source");
    const auto saved=load.momentumRhsLoad();
    rejection([&]{load.validateFor(*mesh,pressure,HermiteGalsIdentity{},color.timeLayer());},"identity");
    rejection([&]{load.validateFor(*mesh,pressure,published(*mesh,shape).state().identity(),color.timeLayer());},"identity");
    for(const auto clock:{HermiteGalsTimeLayer{1,0,0,Position::Initial},HermiteGalsTimeLayer{0,0,0,Position::End},
        HermiteGalsTimeLayer{0,0,1,Position::Initial},HermiteGalsTimeLayer{1,-1,1,Position::Midpoint}})
        rejection([&]{load.validateFor(*mesh,pressure,color.geometryIdentity(),clock);},"time");
    auto nanClock=color.timeLayer();nanClock.end=std::numeric_limits<double>::quiet_NaN();
    rejection([&]{load.validateFor(*mesh,pressure,color.geometryIdentity(),nanClock);},"time");
    const auto newColor=forceColor(*mesh,phase),otherWidth=forceColor(*mesh,phase,.5);
    rejection([&]{load.validateFor(*mesh,pressure,newColor);},"color owner");
    rejection([&]{load.validateFor(*mesh,pressure,otherWidth);},"color owner");
    const auto otherPressure=HydroPressureLoad::build(HydroIntegralDivergence::build(*mesh,phase.layout()));
    rejection([&]{HydroBalancedCsfLoad::constantCurvature(*mesh,otherPressure,color,.5,5);},"divergence owner");
    rejection([&]{load.validateFor(*mesh,otherPressure,color);},"pressure owner");
    const auto rebuiltPressure=HydroPressureLoad::build(phase.metrics().divergence());
    const auto rebuilt=HydroBalancedCsfLoad::constantCurvature(*mesh,rebuiltPressure,color,.5,5);
    rebuilt.validateFor(*mesh,rebuiltPressure,color);check(rebuilt.momentumRhsLoad()==load.momentumRhsLoad(),"same authoritative D changed CSF construction");
    rejection([&]{load.validateFor(*mesh,rebuiltPressure,color);},"pressure owner");
    const auto foreign=HydroPressureLoad::build(metricsFor(*mesh,false).divergence());
    rejection([&]{HydroBalancedCsfLoad::constantCurvature(*mesh,foreign,color,.5,5);},"layout");
    for(double sigma:{-1.,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()})
        rejection([&]{HydroBalancedCsfLoad::constantCurvature(*mesh,pressure,color,sigma,5);},"surface tension");
    for(double kappa:{std::numeric_limits<double>::infinity(),-std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()})
        rejection([&]{HydroBalancedCsfLoad::constantCurvature(*mesh,pressure,color,.5,kappa);},"curvature");
    rejection([&]{HydroBalancedCsfLoad::constantCurvature(*mesh,pressure,color,std::numeric_limits<double>::max(),2);},"coefficient");
    rejection([&]{HydroBalancedCsfLoad::constantCurvature(*mesh,pressure,color,std::numeric_limits<double>::denorm_min(),.5);},"coefficient");
    rejection([&]{HydroBalancedCsfLoad::constantCurvature(*mesh,pressure,color,1e-200,1e-200);},"coefficient");
    auto target=makeMesh(4,1);
    rejection([&]{load.validateFor(*target,pressure,color);},"mesh");
    for(bool rz:{false,true}){shape.radial=rz;const auto original=published(*mesh,shape);HermiteGalsPublicationConfig cfg;cfg.volume.symmetryPlanes=shape.symmetry(*mesh);
        const auto transfer=transferHermiteGalsFirstJetMesh(*mesh,original,*target,HermiteGalsMeshTransferKind::ExactRefinement,cfg);
        check(transfer.report().refinedCells>0&&transfer.publication().state().primary().phaseVolumes.empty(),"actual FirstJet regrid transfer missing");
        const auto nextPhase=HydroPhaseGeometry::build(*target,metricsFor(*target,rz),transfer.publication());const auto nextColor=forceColor(*target,nextPhase);
        const auto nextPressure=HydroPressureLoad::build(nextPhase.metrics().divergence());
        const auto next=HydroBalancedCsfLoad::constantCurvature(*target,nextPressure,nextColor,.5,5);next.validateFor(*target,nextPressure,nextColor);
        algebra(*target,nextPressure,nextColor,Reference(*target,nextPhase),.5,5);
        rejection([&]{next.validateFor(*target,nextPressure,original.state().identity(),nextColor.timeLayer());},"identity");}
    DictionaryMesh::RefineBox box;box.min={.25,.25};box.max={.75,.75};box.level=1;mesh->applyRefineBoxes({box});
    rejection([&]{load.validateFor(*mesh,pressure,color);},"mesh");check(saved==load.momentumRhsLoad(),"failed validation mutated retained load");
}

// Fixed-geometry t=0 acceleration probe. This isolates the spatial CSF,
// sharp dual inertia and native pressure projection from AC and advection.
// The sharp-cell color load is a diagnostic comparator, never a state update.

// Short geometric-CSF operator contract. References remain test-only; every
// candidate force and static response uses the actual production factory.
void geometricOnly()
{
    int localCases=0,loadedFaces=0;double peakCurvature=0;
    for(int scene=0;scene<4;++scene){
        const bool rz=scene==1,ellipse=scene>=2;const int n=scene==3?4:8,amr=scene==3?1:0;
        auto mesh=makeMesh(n,amr);Case shape;shape.pureFirstJet=true;shape.radial=rz;shape.radius=rz?.34:.19;
        const auto publication=[&]{
            if(!ellipse)return published(*mesh,shape);
            const Profile f=[](Vector2D p){const double x=p.x-.5,y=p.y-.5,a=.168,b=.15238095238095239;
                return Jet{x*x/(a*a)+y*y/(b*b)-1,2*x/(a*a),2*y/(b*b),0};};
            const auto state=stateFor(*mesh,f,false,{},[](int,Vector2D,Vector2D){return 0.;},{},HermiteGalsTransport::FirstJet);
            HermiteGalsCandidate candidate(*mesh,{},state.primary());const auto mixed=closeHermiteGalsMixedAmr(*mesh,candidate);
            const auto id=candidate.identity();return publishHermiteGalsFirstJetGeometry(*mesh,std::move(candidate),mixed,{id,centers(*mesh,f)});
        }();
        const auto phase=HydroPhaseGeometry::build(*mesh,metricsFor(*mesh,rz),publication);
        const auto& l=phase.layout();const auto pressure=HydroPressureLoad::build(phase.metrics().divergence());
        const auto originalAlpha=phase.liquidFractions();
        const auto load=HydroBalancedCsfLoad::geometricCapillary(*mesh,pressure,phase,.7);
        load.validateFor(*mesh,pressure,phase.geometryIdentity(),phase.timeLayer());
        check(load.discretization()==HydroCsfDiscretization::GeometricPhaseFraction&&load.isCapillary()
            &&&load.phaseGeometry().liquidFractions()==&phase.liquidFractions()
            &&&load.phaseGeometry().publication().state().geometry().snapshot()==&phase.publication().state().geometry().snapshot(),
            "geometric CSF replaced same phase/publication owner");
        rejects([&]{load.color();},"no regularized color");
        const Reference ref(*mesh,phase);const auto gradient=ref.gradient(originalAlpha);
        std::vector<long double> area(l.pressureCount()),moment(l.pressureCount());
        // Independent composite rule on the original cell graphs. It does
        // not use the factory's accepted panels or its cell means.
        for(std::size_t c=0;c<publication.topology().curves.size();++c)
            for(std::size_t g=0;g<publication.topology().curves[c].graphs().size();++g)
                for(int part=0;part<32;++part)
                    for(const auto& q:publication.quadrature(c,g,double(part)/32,double(part+1)/32)){
                        const int p=l.pressureDof(q.sample.frame.nativeCellIndex);if(p<0)continue;
                        check(q.physicalWeight>0&&q.sample.geometryIdentity==phase.geometryIdentity(),"independent geometric surface source");
                        area[p]+=q.physicalWeight;moment[p]+=static_cast<long double>(q.physicalWeight)*q.sample.curvature;
                    }
        for(int u=0;u<l.velocityCount();++u){
            const double scale=std::abs(gradient[u])+std::abs(load.colorGradient()[u]);
            check(relative(load.colorGradient()[u]-gradient[u],scale)<=1e-12,"independent all-row geometric C alpha action");
            const auto& support=load.supportReports()[u];
            if(support.curvatureDefined){
                long double mean=0;int count=0;
                for(int p:{ref.minus[u],ref.plus[u]})if(p>=0&&area[p]>0){mean+=moment[p]/area[p];++count;}
                check(count>0&&support.radius==0&&support.totalWeight==count,"geometric face local equal-cell support");
                mean/=count;const double curvature=load.curvatureAverages()[u];
                const double error=relative(curvature-double(mean),std::max(std::abs(curvature),std::abs(double(mean))));
                check(error<=1e-11,"independent adjacent-cell arithmetic curvature average");peakCurvature=std::max(peakCurvature,error);
                check(curvature>=support.curvatureMin&&curvature<=support.curvatureMax,"geometric positive curvature range");
                for(const auto& panel:load.quadraturePanels(u)){
                    const int cell=publication.topology().curves[panel.component].graphs()[panel.graph].nativeCellIndex;
                    check((ref.minus[u]>=0&&cell==l.nativeCell(ref.minus[u]))
                        ||(ref.plus[u]>=0&&cell==l.nativeCell(ref.plus[u])),"geometric force borrowed nonadjacent curvature");
                }
                const double expected=.7*curvature*gradient[u];
                check(relative(load.momentumRhsLoad()[u]-expected,std::abs(expected))<=1e-12,"geometric physical RHS sign/source");
                ++loadedFaces;
            }else check(load.momentumRhsLoad()[u]==0,"geometric force used undefined curvature");
            if(ref.area[u]==0)check(load.colorGradient()[u]==0&&load.momentumRhsLoad()[u]==0,"geometric axis physical zero measure");
        }
        if(!ellipse){
            const double kappa=(rz?2.:1.)/shape.radius;
            const auto dual=HydroMomentumGeometry::build(*mesh,phase.metrics());
            const auto mass=HydroMomentumMass::build(*mesh,dual,phase,100,1);
            Matrix block(l.velocityCount(),l.velocityCount());
            for(int u=0;u<l.velocityCount();++u)block.coeffRef(u,u)=mass.masses()[u]/.125;
            block.makeCompressed();
            const auto saddle=HydroSinglePhaseSaddle::buildClosed(pressure,HydroSaddleVelocityBlock::build(l,block));
            const auto result=saddle.solve(*mesh,HydroSaddleRhs(l,load.momentumRhsLoad(),std::vector<double>(l.velocityCount())));
            double umax=0,uscale=0,pError=0;const auto mean=ref.meanByCell(originalAlpha);
            for(int u=0;u<l.velocityCount();++u){umax=std::max(umax,std::abs(result.values().normalVelocity()[u]));
                uscale=std::max(uscale,std::abs(load.momentumRhsLoad()[u])/block.coeff(u,u));}
            for(int p=0;p<l.pressureCount();++p)pError=std::max(pError,
                std::abs(result.values().pressure()[p]-.7*kappa*(originalAlpha[p]-mean[p])));
            check(relative(umax,uscale)<=1e-11&&relative(pError,.7*kappa)<=1e-10,"geometric real-curvature static balance");
            std::cout<<"GEOMETRIC_CSF_STATIC,"<<rz<<','<<n<<",rho,100:1,velocity_relative,"<<relative(umax,uscale)
                <<",pressure_relative,"<<relative(pError,.7*kappa)<<",PASS\n";
        }
        const auto zero=HydroBalancedCsfLoad::geometricCapillary(*mesh,pressure,phase,0);
        check(zero.curvatureReport().projectionSkippedSigmaZero&&zero.curvatureReport().sourceNodes==0,"geometric sigma-zero sampled curvature");
        for(double rhs:zero.momentumRhsLoad())check(rhs==0&&!std::signbit(rhs),"geometric sigma-zero not exact positive zero");
        const auto fresh=HydroPressureLoad::build(phase.metrics().divergence());
        rejects([&]{load.validateFor(*mesh,fresh,phase.geometryIdentity(),phase.timeLayer());},"pressure owner");
        auto limited=HydroCsfCurvatureConfig{};limited.maxSurfaceQueriesPerFace=8;
        rejects([&]{HydroBalancedCsfLoad::geometricCapillary(*mesh,pressure,phase,.7,limited);},"surface query limit");
        check(phase.liquidFractions()==originalAlpha,"geometric CSF mutated sharp material");
        ++localCases;
    }
    check(localCases==4&&loadedFaces>0,"geometric CSF locality/source coverage missing");
    std::cout<<"GEOMETRIC_CSF_ONLY,"<<localCases<<",loaded_faces,"<<loadedFaces<<",curvature_relative,"<<peakCurvature<<",PASS\n";
}


// Independent fixed composite surface quadrature and actual evaluated field.
// It does not consume the factory panels, moments, or symbolic weights when
// checking native virtual work. The constant anchor flux term is retained for
// arbitrary (including non-solenoidal) native test vectors.
void adjointOnly()
{
    int cases=0,staticCases=0;double peakPower=0,peakReference=0,peakWeights=0;
    bool observedVariableWork=false;
    for(int scene=0;scene<4;++scene){
        const bool rz=scene==1||scene==3,ellipse=scene>=2;constexpr int n=8;constexpr double sigma=.7;
        auto mesh=makeMesh(n);Case shape;shape.pureFirstJet=true;shape.radial=rz;shape.radius=rz?.34:.19;
        const auto publication=[&]{
            if(!ellipse)return published(*mesh,shape);
            const double a=rz?.34:.168,b=rz?.30:.15238095238095239,cx=rz?0:.5;
            const Profile f=[=](Vector2D p){const double x=p.x-cx,y=p.y-.5;
                return Jet{x*x/(a*a)+y*y/(b*b)-1,2*x/(a*a),2*y/(b*b),0};};
            // Radial FirstJet owns the r=0 axis identity automatically,
            // just as published(Case radial=true) above. Do not redeclare it.
            const std::vector<LevelSetInterfacePanels::SymmetryPlane> planes;
            const auto state=stateFor(*mesh,f,rz,{},[](int,Vector2D,Vector2D){return 0.;},planes,HermiteGalsTransport::FirstJet);
            HermiteGalsCandidate candidate(*mesh,{},state.primary());const auto mixed=closeHermiteGalsMixedAmr(*mesh,candidate,planes);
            const auto id=candidate.identity();HermiteGalsPublicationConfig cfg;cfg.volume.symmetryPlanes=planes;
            return publishHermiteGalsFirstJetGeometry(*mesh,std::move(candidate),mixed,{id,centers(*mesh,f)},cfg);
        }();
        const auto phase=HydroPhaseGeometry::build(*mesh,metricsFor(*mesh,rz),publication);const auto alpha=phase.liquidFractions();
        const auto& l=phase.layout();const auto pressure=HydroPressureLoad::build(phase.metrics().divergence());
        const auto geometry=HydroMomentumGeometry::build(*mesh,phase.metrics());
        std::vector<int> walls;
        for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(*mesh,l.nativeFace(u));
            if(f.negativeCell<0||f.positiveCell<0)walls.push_back(u);}
        const auto law=HydroAdvectorBoundary::impermeable(*mesh,geometry,walls);
        const auto rec=HydroAdvectorReconstruction::prepare(*mesh,geometry,law);
        const auto load=HydroBalancedCsfLoad::adjointCapillary(*mesh,pressure,phase,sigma,rec,&law);
        load.validateFor(*mesh,pressure,phase.geometryIdentity(),phase.timeLayer(),rec,&law);
        check(load.discretization()==HydroCsfDiscretization::GeometricAdjoint&&load.reconstruction().sameIdentity(rec)
            &&load.boundary()&&load.boundary()->sameIdentity(law),"adjoint lost true reconstruction/boundary owner");
        check(&load.phaseGeometry().cellVolumes()==&phase.cellVolumes()
            &&&load.phaseGeometry().publication().state().geometry().snapshot()==&publication.state().geometry().snapshot(),
            "adjoint replaced actual material/publication owner");
        rejects([&]{load.curvatureAverages();},"face curvature multiplier");
        rejects([&]{load.color();},"no regularized color");
        std::vector<double> velocity(l.velocityCount());
        for(int u=0;u<l.velocityCount();++u)velocity[u]=std::sin(.731*(u+1))+.3*std::cos(.17*(u+1)*(u+3));
        for(int u:walls)velocity[u]=0;
        const auto field=rec.evaluate(*mesh,HydroNativeValues(l,velocity,std::vector<double>(l.pressureCount())),
            phase.geometryIdentity(),phase.timeLayer(),law);
        long double area=0,kappaMoment=0,flux=0,sharpPower=0,workScale=0;
        std::vector<std::optional<LAMP::Hydrodynamics::Detail::HydroRegularCellStencil>> stencils(l.pressureCount());
        for(std::size_t c=0;c<publication.topology().curves.size();++c)
            for(std::size_t g=0;g<publication.topology().curves[c].graphs().size();++g)
                for(int part=0;part<64;++part)
                    for(const auto& q:publication.quadrature(c,g,double(part)/64,double(part+1)/64)){
                        const auto& f=q.sample.frame;const int p=l.pressureDof(f.nativeCellIndex);
                        check(p>=0&&q.physicalWeight>0&&q.sample.geometryIdentity==phase.geometryIdentity(),"adjoint independent surface owner");
                        const auto actual=field.cellVelocity(p,f.point);
                        const long double nv=static_cast<long double>(f.normal.x)*actual.x+static_cast<long double>(f.normal.y)*actual.y;
                        area+=q.physicalWeight;kappaMoment+=static_cast<long double>(q.physicalWeight)*q.sample.curvature;
                        flux+=q.physicalWeight*nv;sharpPower-=sigma*q.physicalWeight*q.sample.curvature*nv;
                        workScale+=std::abs(sigma*q.physicalWeight*q.sample.curvature*nv);
                        // Separate exact-source symbolic application check; the
                        // independent power integral above used actual evaluate.
                        if(!stencils[p])stencils[p]=rec.regularGeometry().cellStencil(*mesh,p,&law);
                        long double x=0,y=0,scale=0;
                        for(const auto& row:rec.regularGeometry().cellJetWeights(*stencils[p],f.point)){
                            x+=static_cast<long double>(row.jet[0])*velocity[row.velocityDof];
                            y+=static_cast<long double>(row.jet[1])*velocity[row.velocityDof];
                            scale+=std::abs(row.jet[0]*velocity[row.velocityDof])+std::abs(row.jet[1]*velocity[row.velocityDof]);
                        }
                        const auto error=relative(std::hypot(double(x)-actual.x,double(y)-actual.y),
                            std::max(1.,double(scale)));
                        check(error<=1e-12,"adjoint symbolic weights differ from actual source");peakWeights=std::max(peakWeights,error);
                    }
        check(area>0,"adjoint independent surface area missing");
        const double k0=double(kappaMoment/area);std::vector<double> potential(l.pressureCount());
        for(int p=0;p<l.pressureCount();++p){
            check(std::isfinite(load.referenceCurvatures()[p]),"adjoint connected interface reference undefined");
            const auto error=relative(load.referenceCurvatures()[p]-k0,std::max(std::abs(k0),std::abs(load.referenceCurvatures()[p])));
            check(error<=1e-11,"adjoint reference not actual component positive surface mean");peakReference=std::max(peakReference,error);
            potential[p]=sigma*load.referenceCurvatures()[p]*alpha[p];
        }
        const Reference independent(*mesh,phase);const auto base=independent.gradient(potential);
        long double nativePower=0,basePower=0,residualPower=0,nativeScale=0;
        for(int u=0;u<l.velocityCount();++u){
            close(load.referencePressureLoad()[u],base[u],1e-12*std::max(1.,std::abs(base[u])),"adjoint independent native reference pressure");
            check(load.momentumRhsLoad()[u]==load.referencePressureLoad()[u]+load.adjointResidualLoad()[u],
                "adjoint actual force is not exact base plus residual");
            nativePower+=static_cast<long double>(velocity[u])*load.momentumRhsLoad()[u];
            basePower+=static_cast<long double>(velocity[u])*load.referencePressureLoad()[u];
            residualPower+=static_cast<long double>(velocity[u])*load.adjointResidualLoad()[u];
            nativeScale+=std::abs(static_cast<long double>(velocity[u])*load.momentumRhsLoad()[u]);
            for(const auto& panel:load.quadraturePanels(u)){
                const int p=l.pressureDof(publication.topology().curves[panel.component].graphs()[panel.graph].nativeCellIndex);
                const auto q=publication.sample(panel.component,panel.graph,(panel.begin+panel.end)/2);
                if(!stencils[p])stencils[p]=rec.regularGeometry().cellStencil(*mesh,p,&law);
                bool exists=false;for(const auto& row:rec.regularGeometry().cellJetWeights(*stencils[p],q.frame.point))
                    if(row.velocityDof==u)exists=true;
                check(exists,"adjoint borrowed curvature outside same source symbolic support");
            }
        }
        // For arbitrary U: native power = sharp surface work plus the actual
        // reference pressure power and reference-curvature surface flux.
        const long double expected=sharpPower+basePower+sigma*load.referenceCurvatures()[0]*flux;
        const double powerError=relative(double(nativePower-expected),double(std::max(1.L,workScale+nativeScale+std::abs(basePower))));
        check(powerError<=1e-11,"adjoint independent native/surface virtual work");peakPower=std::max(peakPower,powerError);
        if(ellipse&&std::abs(residualPower)>1e-7*std::max(1.L,workScale))observedVariableWork=true;
        if(!ellipse){
            const double exact=(rz?2.:1.)/shape.radius;
            const auto mass=HydroMomentumMass::build(*mesh,geometry,phase,100,1);
            Matrix block(l.velocityCount(),l.velocityCount());
            for(int u=0;u<l.velocityCount();++u)block.coeffRef(u,u)=mass.masses()[u]/.125;
            block.makeCompressed();const auto saddle=HydroSinglePhaseSaddle::buildClosed(pressure,HydroSaddleVelocityBlock::build(l,block));
            const auto solved=saddle.solve(*mesh,HydroSaddleRhs(l,load.momentumRhsLoad(),std::vector<double>(l.velocityCount())));
            double umax=0,uscale=0,pError=0;const auto mean=independent.meanByCell(alpha);
            for(int u=0;u<l.velocityCount();++u){umax=std::max(umax,std::abs(solved.values().normalVelocity()[u]));
                uscale=std::max(uscale,std::abs(load.momentumRhsLoad()[u])/block.coeff(u,u));}
            for(int p=0;p<l.pressureCount();++p)pError=std::max(pError,std::abs(solved.values().pressure()[p]-sigma*exact*(alpha[p]-mean[p])));
            check(relative(umax,uscale)<=1e-11&&relative(pError,sigma*exact)<=1e-10,"adjoint actual circle/sphere static balance");
            ++staticCases;std::cout<<"ADJOINT_CSF_STATIC,"<<rz<<','<<n<<",rho,100:1,velocity_relative,"<<relative(umax,uscale)
                <<",pressure_relative,"<<relative(pError,sigma*exact)<<",PASS\n";
        }
        const auto zero=HydroBalancedCsfLoad::adjointCapillary(*mesh,pressure,phase,0,rec,&law);
        check(zero.curvatureReport().projectionSkippedSigmaZero&&zero.curvatureReport().sourceNodes==0,"adjoint sigma zero sampled surface");
        for(double value:zero.momentumRhsLoad())check(value==0&&!std::signbit(value),"adjoint sigma zero not literal positive zero");
        auto changed=phase.timeLayer();changed.step=1;changed.begin=0;changed.end=.125;changed.position=Position::End;
        rejects([&]{load.validateFor(*mesh,pressure,phase.geometryIdentity(),changed,rec,&law);},"full time layer");
        const auto fresh=HydroAdvectorReconstruction::prepare(*mesh,geometry,law);
        rejects([&]{load.validateFor(*mesh,pressure,phase.geometryIdentity(),phase.timeLayer(),fresh,&law);},"reconstruction owner");
        const auto otherLaw=HydroAdvectorBoundary::impermeable(*mesh,geometry,walls);
        rejects([&]{load.validateFor(*mesh,pressure,phase.geometryIdentity(),phase.timeLayer(),rec,&otherLaw);},"boundary owner");
        auto limited=HydroCsfCurvatureConfig{};limited.maxSurfaceQueriesPerFace=8;
        rejects([&]{HydroBalancedCsfLoad::adjointCapillary(*mesh,pressure,phase,sigma,rec,&law,limited);},"surface query limit");
        if(scene==2){const auto modal=rec.withTangentialModes(*mesh,phase);
            check(!modal.tangentialModes().empty(),"adjoint tangential guard fixture missing");
            rejects([&]{HydroBalancedCsfLoad::adjointCapillary(*mesh,pressure,phase,sigma,modal,&law);},"independent tangential modes");}
        check(phase.liquidFractions()==alpha&&phase.publication().state().primary().phaseVolumes.empty(),
            "adjoint modified raw geometry or introduced a phase inventory");
        ++cases;std::cout<<"ADJOINT_CSF_POWER,"<<scene<<','<<rz<<','<<ellipse<<",native,"<<double(nativePower)
            <<",sharp,"<<double(sharpPower)<<",base,"<<double(basePower)<<",reference_flux,"<<double(sigma*load.referenceCurvatures()[0]*flux)
            <<",residual_work,"<<double(residualPower)<<",relative_error,"<<powerError<<",PASS\n";
    }
    check(cases==4&&staticCases==2&&observedVariableWork,"adjoint constant/variable curvature physical coverage missing");
    std::cout<<"ADJOINT_CSF_ONLY,"<<cases<<",static_cases,"<<staticCases<<",power_relative,"<<peakPower
        <<",reference_relative,"<<peakReference<<",same_source_weights_relative,"<<peakWeights<<",PASS\n";
}

// Independent primary-rate application uses the actual evaluated field/J,
// then a scalar rate. It never uses native symbolic B_u. C2 macro rates
// are completed/interpolated FORWARD, independently of the load's full E^T
// source pullback. The original C1 fixture retains its exact analytic
// Hessian jump/donor reference; neither branch injects analytic curvature.
void projectedAdjointOnly(bool macro=false)
{
    constexpr int n=8;constexpr double sigma=.7;
    int cases=0;double peakPower=0;bool observedResidual=false,observedNonlocal=false;
    for(int scene=0;scene<4;++scene){
        auto mesh=makeMesh(n);const bool jumped=scene==1||scene==2,symmetry=scene==3;const double direction=scene==2?-1.:1.;
        const Profile f=[=](Vector2D p){const double x=p.x-(symmetry?0:.5),y=p.y-.5;
            const double factor=jumped&&x>0?1.35:1.;
            return Jet{factor*x*x+(symmetry?1.2:1.)*y*y-.19*.19,2*factor*x,2*(symmetry?1.2:1.)*y,0};};
        const std::vector<LevelSetInterfacePanels::SymmetryPlane> planes=symmetry?
            std::vector<LevelSetInterfacePanels::SymmetryPlane>{{0,0}}:std::vector<LevelSetInterfacePanels::SymmetryPlane>{};
        const auto initial=stateFor(*mesh,f,false,{},[](int,Vector2D,Vector2D){return 0.;},planes,HermiteGalsTransport::FirstJet);
        HermiteGalsCandidate candidate(*mesh,{},initial.primary());const auto closure=closeHermiteGalsMixedAmr(*mesh,candidate,planes);
        const auto id=candidate.identity();
        HermiteGalsPublicationConfig cfg;cfg.volume.symmetryPlanes=planes;
        if(macro)cfg.scalarRepresentation=HermiteGalsScalarRepresentation::C2MacroQuartic;
        const auto pub=publishHermiteGalsFirstJetGeometry(*mesh,std::move(candidate),closure,{id,centers(*mesh,f)},cfg);
        const auto phase=HydroPhaseGeometry::build(*mesh,metricsFor(*mesh,false),pub);const auto& l=phase.layout();
        const auto pressure=HydroPressureLoad::build(phase.metrics().divergence());
        const auto geometry=HydroMomentumGeometry::build(*mesh,phase.metrics());
        std::vector<int> walls;std::vector<double> velocity(l.velocityCount());
        for(int u=0;u<l.velocityCount();++u){const auto face=l.domain().face(*mesh,l.nativeFace(u));
            if(face.negativeCell<0||face.positiveCell<0)walls.push_back(u);
            // X coordinate alone: exact quadratic face means, hence the same
            // P2 field u_x=direction*x*(1-x), u_y=0, with literal closed walls.
            if(l.normalAxis(u)==HydroDofLayout::NormalAxis::First)
                velocity[u]=direction*face.center.x*(1-face.center.x);
            else if(symmetry)velocity[u]=face.center.y*(1-face.center.y)*(1+.3*face.center.x);
        }
        for(int u:walls)velocity[u]=0;
        const auto law=HydroAdvectorBoundary::impermeable(*mesh,geometry,walls);
        const auto rec=HydroAdvectorReconstruction::prepare(*mesh,geometry,law);
        const auto field=rec.evaluate(*mesh,HydroNativeValues(l,velocity,std::vector<double>(l.pressureCount())),
            phase.geometryIdentity(),phase.timeLayer(),law);
        const auto mixed=HermiteGalsAmrClosureSensitivity::buildUniformFirstJet(*mesh,pub.state(),planes);
        const auto load=HydroBalancedCsfLoad::projectedMotionCapillary(*mesh,pressure,phase,sigma,field,&mixed);
        load.validateFor(*mesh,pressure,phase.geometryIdentity(),phase.timeLayer(),field,&mixed);
        check(load.discretization()==HydroCsfDiscretization::GeometricProjectedAdjoint
            &&load.projectedMotionAdvector().identity()==field.identity()
            &&load.mixedSensitivity().primaryIdentity()==pub.state().identity(),"projected owner provenance");
        rejects([&]{load.validateFor(*mesh,pressure,phase.geometryIdentity(),phase.timeLayer(),rec,&law);},
            "actual motion direction owner");
        rejects([&]{load.curvatureAverages();},"face curvature multiplier");
        const auto& scalar=pub.state().geometry().snapshot();const auto& primary=pub.state().primary();
        check(scalar.usesC2Macro()==macro,"projected fixture has a different scalar representation");
        const auto vertexPoint=[&](int v){return mesh->verts(mesh->vertIndexToLabel[v]).coord();};
        // Closed-form upwind scalar Hessian, independent of production donor
        // selection. Actual field/J is read from the retained physical source.
        const auto rateAt=[&](Vector2D q){
            int cell=-1;for(int c:scalar.nativeCellCandidates(q,q))if(scalar.nativeCellActive(c)&&(cell<0||c<cell))cell=c;
            check(cell>=0,"independent projected point outside source");
            const auto u=field.cellVelocityJet(l.pressureDof(cell),q).jet;
            const auto analytic=f(q);const auto actual=scalar.scalarJetInNativeCell(cell,q);
            const double gx=macro?actual.dx:analytic.dx,gy=macro?actual.dy:analytic.dy;
            // C1 retains an independent exact backward-donor Hessian. Macro
            // is C2: every neighboring native/half chart gives the same H,
            // so the actual scalar jet is the unique physical rate source.
            const double hxx=macro?actual.dxx:2*(jumped&&(q.x>.5||(q.x==.5&&u.value.x<0))?1.35:1.);
            const double hxy=macro?actual.dxy:0,hyy=macro?actual.dyy:2*(symmetry?1.2:1.);
            auto rates=std::array<long double,3>{-static_cast<long double>(u.value.x)*gx-static_cast<long double>(u.value.y)*gy,
                -static_cast<long double>(u.derivative.xx)*gx-static_cast<long double>(u.derivative.yx)*gy-hxx*u.value.x-hxy*u.value.y,
                -static_cast<long double>(u.derivative.xy)*gx-static_cast<long double>(u.derivative.yy)*gy-hxy*u.value.x-hyy*u.value.y};
            if(symmetry&&q.x==0)rates[1]=0; // Actual SymmetryFirstJet primary constraint.
            return rates;
        };
        std::map<int,std::array<long double,4>> vertices;
        std::vector<std::array<long double,2>> edges(primary.edges.size());
        std::map<int,long double> centerRates;
        for(const auto& v:primary.vertices){const auto a=rateAt(vertexPoint(v.vertex));
            vertices[v.vertex]={a[0],a[1],a[2],0};}
        for(std::size_t k=0;k<primary.edges.size();++k){const auto& e=primary.edges[k];
            const auto a=vertexPoint(e.vertices[0]),b=vertexPoint(e.vertices[1]);const auto r=rateAt(a+(b-a)*.5);
            edges[k]={r[0],r[a.y==b.y?2:1]};}
        for(std::size_t k=0;k<primary.vertices.size();++k){
            long double rate=0;for(const auto& term:mixed.rows()[k]){
                const auto value=term.component==HermiteGalsPrimaryDerivative::EdgeNormal?edges.at(term.primaryOrdinal)[1]:
                    vertices.at(primary.vertices.at(term.primaryOrdinal).vertex)[term.component==HermiteGalsPrimaryDerivative::VertexDx?1:2];
                rate+=static_cast<long double>(term.weight)*value;
            }vertices[primary.vertices[k].vertex][3]=rate;
        }
        for(int cell:primary.geometryCells){const auto b=scalar.nativeCellBounds(cell);
            centerRates[cell]=rateAt(b[0]+(b[1]-b[0])*.5)[0];}
        std::optional<LevelSetHermiteGeometry::Snapshot> completedRate;
        if(macro){
            LAMP::Geometry::ActiveQuarticHermitePatch ratePatch;
            ratePatch.activeCells.assign(mesh->cellIndexNum(),1);
            ratePatch.nodes.resize(mesh->vertIndexNum());ratePatch.cellCenters.resize(mesh->cellIndexNum());
            for(const auto& entry:vertices){const auto& r=entry.second;
                ratePatch.nodes.at(entry.first)={double(r[0]),double(r[1]),double(r[2]),double(r[3])};}
            for(std::size_t e=0;e<primary.edges.size();++e)
                ratePatch.edgeMidpoints.push_back({primary.edges[e].vertices,double(edges[e][0]),double(edges[e][1])});
            for(const auto& entry:centerRates)ratePatch.cellCenters.at(entry.first)=double(entry.second);
            const auto completion=LAMP::Geometry::ActiveQuarticHermiteCompletion::fromUniformPatch(*mesh,ratePatch,planes);
            completedRate.emplace(LevelSetHermiteGeometry::Snapshot::fromUniformC2Completion(*mesh,completion));
        }
        const auto scalarRate=[&](const LevelSetHermiteGeometry::Frame& frame){
            const auto chart=LevelSetHermiteGeometry::cellChartPoint(frame);
            long double dot=0;
            if(macro){
                dot=completedRate->scalarJetAtNativeCellChart(frame.nativeCellIndex,chart).value;
                using Kind=LAMP::Geometry::ActiveQuarticHermiteSensitivity::Kind;
                for(const auto& term:scalar.nativeC2MacroDifferentialSensitivity(frame.nativeCellIndex,chart,{1,0,0,0,0,0}))
                    observedNonlocal=observedNonlocal||(term.kind==Kind::CellCenter&&term.index!=frame.nativeCellIndex&&term.coefficient!=0);
            }else for(const auto& term:scalar.nativeQuarticDifferentialSensitivity(frame.nativeCellIndex,
                chart,{1,0,0,0,0,0})){
                using Kind=LAMP::Geometry::ActiveQuarticHermiteSensitivity::Kind;
                const long double value=term.kind==Kind::Vertex?vertices.at(term.index)[term.component]:
                    term.kind==Kind::EdgeMidpoint?edges.at(term.index)[term.component]:centerRates.at(term.index);
                dot+=term.coefficient*value;
            }
            const auto j=scalar.scalarJetAtNativeCellChart(frame.nativeCellIndex,LevelSetHermiteGeometry::cellChartPoint(frame));
            return -dot/std::hypot(static_cast<long double>(j.dx),static_cast<long double>(j.dy));
        };
        long double surfacePower=0,projectedFlux=0,scale=0;
        for(std::size_t c=0;c<pub.topology().curves.size();++c)
            for(std::size_t g=0;g<pub.topology().curves[c].graphs().size();++g)
                for(int part=0;part<64;++part)
                    for(const auto& q:pub.quadrature(c,g,double(part)/64,double(part+1)/64)){
                        const auto rate=scalarRate(q.sample.frame);
                        projectedFlux+=q.physicalWeight*rate;const auto work=-sigma*q.physicalWeight*q.sample.curvature*rate;
                        surfacePower+=work;scale+=std::abs(work);
                    }
        long double power=0,basePower=0,residualPower=0;
        for(int u=0;u<l.velocityCount();++u){
            power+=static_cast<long double>(velocity[u])*load.momentumRhsLoad()[u];
            basePower+=static_cast<long double>(velocity[u])*load.referencePressureLoad()[u];
            residualPower+=static_cast<long double>(velocity[u])*load.adjointResidualLoad()[u];
            scale+=std::abs(static_cast<long double>(velocity[u])*load.momentumRhsLoad()[u]);
            check(load.momentumRhsLoad()[u]==load.referencePressureLoad()[u]+load.adjointResidualLoad()[u],
                "projected physical load base/residual identity");
        }
        const auto expected=surfacePower+basePower+sigma*load.referenceCurvatures()[0]*projectedFlux;
        const double error=relative(double(power-expected),double(std::max(1.L,scale+std::abs(basePower))));
        check(error<=1e-11,"independent projected primary/scalar surface virtual work");peakPower=std::max(peakPower,error);
        if(jumped&&std::abs(residualPower)>1e-7)observedResidual=true;
        if(scene==0){
            const auto mass=HydroMomentumMass::buildCellAverageProjection(*mesh,geometry,phase,100,1);
            Matrix block(l.velocityCount(),l.velocityCount());
            for(int u=0;u<l.velocityCount();++u)block.coeffRef(u,u)=mass.masses()[u]/.125;
            block.makeCompressed();const auto saddle=HydroSinglePhaseSaddle::buildClosed(pressure,HydroSaddleVelocityBlock::build(l,block));
            const auto solved=saddle.solve(*mesh,HydroSaddleRhs(l,load.momentumRhsLoad(),std::vector<double>(l.velocityCount())));
            double umax=0,uscale=0;
            for(int u=0;u<l.velocityCount();++u){umax=std::max(umax,std::abs(solved.values().normalVelocity()[u]));
                uscale=std::max(uscale,std::abs(load.momentumRhsLoad()[u])/block.coeff(u,u));}
            check(relative(umax,uscale)<=1e-11,"projected constant circle static balance");
        }
        if(symmetry){
            bool nonzeroRaw=false;
            for(const auto& v:primary.vertices)if(vertexPoint(v.vertex).x==0){
                const auto q=vertexPoint(v.vertex);const auto u=field.cellVelocityJet(l.pressureDof(scalar.nativeCellCandidates(q,q).front()),q).jet;
                const auto g=f(q);const double raw=-u.derivative.yx*g.dy;
                nonzeroRaw=nonzeroRaw||std::abs(raw)>1e-5;
                check(vertices.at(v.vertex)[1]==0,"symmetry normal-firstjet rate not literal zero");
            }
            check(nonzeroRaw,"wall-rate fixture did not exercise actual scalar symmetry override");
        }
        const auto zero=HydroBalancedCsfLoad::projectedMotionCapillary(*mesh,pressure,phase,0,field,nullptr);
        zero.validateFor(*mesh,pressure,phase.geometryIdentity(),phase.timeLayer(),field,nullptr);
        check(zero.curvatureReport().projectionSkippedSigmaZero&&zero.curvatureReport().sourceNodes==0,"projected sigma zero work");
        rejects([&]{zero.mixedSensitivity();},"no projected mixed owner");
        for(double value:zero.momentumRhsLoad())check(value==0&&!std::signbit(value),"projected sigma zero literal positive zero");
        rejects([&]{HydroBalancedCsfLoad::projectedMotionCapillary(*mesh,pressure,phase,sigma,field,nullptr);},"lacks mixed closure");
        const auto otherField=rec.evaluate(*mesh,HydroNativeValues(l,velocity,std::vector<double>(l.pressureCount())),
            phase.geometryIdentity(),phase.timeLayer(),law);
        rejects([&]{load.validateFor(*mesh,pressure,phase.geometryIdentity(),phase.timeLayer(),otherField,&mixed);},"direction owner");
        auto changed=phase.timeLayer();changed.step=1;changed.begin=0;changed.end=.125;changed.position=Position::End;
        rejects([&]{load.validateFor(*mesh,pressure,phase.geometryIdentity(),changed,field,&mixed);},"full time layer");
        const auto open=rec.evaluate(*mesh,HydroNativeValues(l,velocity,std::vector<double>(l.pressureCount())),
            phase.geometryIdentity(),phase.timeLayer());
        rejects([&]{HydroBalancedCsfLoad::projectedMotionCapillary(*mesh,pressure,phase,sigma,open,&mixed);},"closed boundary");
        // Exact tangent tie: H differs, but H*u is invariant when the native
        // field points along the join. Zero velocity trivially shares that limit.
        if(jumped&&!macro){
            const Vector2D q{.5,.375};long double baselineX=0,baselineY=0;bool first=true,differentH=false;
            for(int cell:scalar.nativeCellCandidates(q,q)){const auto j=scalar.scalarJetInNativeCell(cell,q);
                const long double x=j.dxy*.375*(1-.375),y=j.dyy*.375*(1-.375);
                if(first){baselineX=x;baselineY=y;first=false;}
                else{close(double(x-baselineX),0,1e-12,"C1 tangential Hessian-action tie");
                    close(double(y-baselineY),0,1e-12,"C1 tangential Hessian-action tie");}
                if(j.dxx>2.5)differentH=true;
            }check(!first&&differentH,"directional donor fixture lacks real Hessian jump");
        }
        ++cases;std::cout<<(macro?"PROJECTED_MACRO_POWER,":"PROJECTED_ADJOINT_POWER,")<<scene<<",native,"<<double(power)<<",projected_surface,"<<double(surfacePower)
            <<",base,"<<double(basePower)<<",reference_flux,"<<double(sigma*load.referenceCurvatures()[0]*projectedFlux)
            <<",residual,"<<double(residualPower)<<",relative,"<<error<<",PASS\n";
    }
    check(cases==4&&observedResidual,"projected directional variable-curvature coverage");
    check(!macro||observedNonlocal,"macro projected rate lacks neighboring center source coverage");
    std::cout<<(macro?"PROJECTED_MACRO_ONLY,":"PROJECTED_ADJOINT_ONLY,")<<cases<<",power_relative,"<<peakPower<<",PASS\n";
}


// HF is a new force discretization. Its circle values are approximate and are
// never tested against the old exact same-scalar curvature/static-drop gates.
std::unique_ptr<QuadtreeMesh> heightMesh(int n,Vector2D origin={0,0})
{
    DictionaryMesh cfg;cfg.Nx=cfg.Ny=n;cfg.minCoord=origin;cfg.maxCoord=origin+Vector2D{1,1};
    cfg.maxCellNum=16384;cfg.maxInitialIterations=0;return std::make_unique<QuadtreeMesh>(cfg);
}
HermiteGalsPublication heightPublication(const QuadtreeMesh& mesh,const Profile& profile,
    HermiteGalsTimeLayer clock={},bool macro=false)
{
    const auto seed=stateFor(mesh,profile,false,clock,[](int,Vector2D,Vector2D){return 0.;},{},HermiteGalsTransport::FirstJet);
    HermiteGalsCandidate candidate(mesh,clock,seed.primary());const auto mixed=closeHermiteGalsMixedAmr(mesh,candidate);
    HermiteGalsPublicationConfig cfg;if(macro)cfg.scalarRepresentation=HermiteGalsScalarRepresentation::C2MacroQuartic;
    const auto id=candidate.identity();return publishHermiteGalsFirstJetGeometry(mesh,std::move(candidate),mixed,{id,centers(mesh,profile)},cfg);
}
HydroCsfCurvatureConfig heightConfig()
{HydroCsfCurvatureConfig cfg;cfg.forceMethod=HydroCsfForceMethod::HeightFunction;return cfg;}
void heightForceAlgebra(const QuadtreeMesh& mesh,const HydroPressureLoad& pressure,
    const HydroPhaseGeometry& phase,const HydroBalancedCsfLoad& load,double sigma)
{
    const auto& layout=phase.layout();load.validateFor(mesh,pressure,phase.geometryIdentity(),phase.timeLayer());
    check(load.discretization()==HydroCsfDiscretization::GeometricHeightFunction&&load.isCapillary(),"HF actual discretization kind");
    check(&load.phaseGeometry().cellVolumes()==&phase.cellVolumes()&&&load.phaseGeometry().liquidFractions()==&phase.liquidFractions()
        &&&load.phaseGeometry().publication().state().geometry().snapshot()==&phase.publication().state().geometry().snapshot(),
        "HF replaced retained volume/scalar owner");
    check(load.geometryIdentity()==phase.geometryIdentity()&&load.timeLayer().step==phase.timeLayer().step
        &&load.timeLayer().begin==phase.timeLayer().begin&&load.timeLayer().end==phase.timeLayer().end
        &&load.timeLayer().position==phase.timeLayer().position,"HF source full clock changed");
    check(load.heightFunctionCellCurvatures().size()==static_cast<std::size_t>(layout.pressureCount()),"HF pressure-cell curvature layout");
    const Reference incidence(mesh,phase);const auto gradient=incidence.gradient(phase.liquidFractions());
    const auto actualC=pressure.apply(mesh,HydroNativeValues(layout,std::vector<double>(layout.velocityCount()),phase.liquidFractions()));
    double error=0,scale=0;int loaded=0;
    for(int u=0;u<layout.velocityCount();++u){
        error=std::hypot(error,load.colorGradient()[u]-gradient[u]);scale=std::hypot(scale,gradient[u]);
        close(load.colorGradient()[u],actualC[u],1e-13*std::max(std::abs(actualC[u]),std::abs(gradient[u])),"HF uses actual authoritative C alpha");
        check(load.quadraturePanels(u).empty()&&!load.supportReports()[u].curvatureErrorEstimateAvailable,
            "HF manufactured quadrature panels/curvature error certificate");
        const auto& support=load.supportReports()[u];
        if(support.curvatureDefined){
            long double sum=0;int count=0;double lo=std::numeric_limits<double>::infinity(),hi=-lo;
            for(int p:{incidence.minus[u],incidence.plus[u]})if(p>=0&&std::isfinite(load.heightFunctionCellCurvatures()[p])){
                const double k=load.heightFunctionCellCurvatures()[p];sum+=k;++count;lo=std::min(lo,k);hi=std::max(hi,k);}
            check(count>0&&support.radius==0,"HF borrowed nonadjacent height curvature");
            const double k=load.curvatureAverages()[u];close(k,sum/count,1e-13*std::max(1.,std::abs(k)),"HF adjacent-cell equal mean");
            check(k>=lo-1e-14&&k<=hi+1e-14,"HF face curvature outside actual adjacent range");
            const double expected=sigma*k*gradient[u];close(load.momentumRhsLoad()[u],expected,
                1e-12*std::abs(expected),"HF physical RHS sigma kappa C alpha sign/area");++loaded;
        }else check(load.momentumRhsLoad()[u]==0,"HF used undefined force curvature");
    }
    check(relative(error,scale)<=1e-12,"HF independent native incidence gradient");
    if(phase.report().cutCells&&sigma>0)check(loaded>0,"HF has no force support on resolved interface");
}
void heightFourierSymbol()
{
    // Independent analytic column means, including the cell-average sinc.
    // This establishes the proposed compact derivative's restoring symbol;
    // the actual nonlinear factory is tested separately below.
    constexpr long double h=.03125L;int modes=0;
    for(int mode:{1,4,8,12,16}){
        const long double theta=2*pi*mode/32,mean=std::sin(theta/2)/(theta/2);
        const long double hm=mean*std::cos(theta),h0=mean,hp=hm;
        const long double second=(hp-2*h0+hm)/(h*h);
        const long double symbol=-4*std::pow(std::sin(theta/2),2)*mean/(h*h);
        close(second,symbol,1e-12*std::abs(symbol),"HF independent three-column Fourier symbol");
        check(second<0&&(-second)*h0>0,"HF high-frequency curvature reverses restoring sign");++modes;
    }
    std::cout<<"HEIGHT_FUNCTION_SYMBOL,modes,"<<modes<<",through_Nyquist,PASS\n";
}
void heightSmallContract()
{
    heightFourierSymbol();auto mesh=makeMesh(8);const auto phase=HydroPhaseGeometry::build(*mesh,metricsFor(*mesh,false),purePublication(*mesh,false,1));
    const auto pressure=HydroPressureLoad::build(phase.metrics().divergence());
    const auto load=HydroBalancedCsfLoad::heightFunctionCapillary(*mesh,pressure,phase,0,heightConfig());
    heightForceAlgebra(*mesh,pressure,phase,load,0);
    check(load.curvatureReport().projectionSkippedSigmaZero&&load.curvatureReport().heightFunctionCells==0
        &&load.curvatureReport().heightFunctionColumns==0&&load.curvatureReport().heightFunctionCellQueries==0,"HF sigma-zero evaluated height data");
    for(double value:load.momentumRhsLoad())check(value==0&&!std::signbit(value),"HF sigma-zero not exact positive zero");
    rejects([&]{load.referenceCurvatures();},"no adjoint reference");
    rejects([&]{load.projectedMotionAdvector();},"no projected motion");
    std::cout<<"HEIGHT_FUNCTION_SMALL,pure_cells,"<<phase.layout().pressureCount()<<",PASS\n";
}
Profile heightCapsule(bool normalX,bool perturb=false)
{
    return [=](Vector2D p){
        const double t=(normalX?p.y:p.x)-.5,n=(normalX?p.x:p.y)-.5;
        const double a=perturb?.25:.0625,b=perturb?.1868:.32;
        double eta=0,etaPrime=0;
        if(perturb){
            constexpr double amplitude=.01/64,wave=32*double(pi);const double r=std::abs(t);
            double taper=0,taperPrime=0;
            if(r<=.125)taper=1;
            else if(r<.1875){const double s=(r-.125)/.0625;taper=1-3*s*s+2*s*s*s;
                taperPrime=(-6*s+6*s*s)/.0625*(t>0?1:-1);}
            eta=amplitude*taper*std::cos(wave*t);
            etaPrime=amplitude*(taperPrime*std::cos(wave*t)-taper*wave*std::sin(wave*t));
        }
        const double cap=std::max(std::abs(t)-a,0.),dn=2*(n-eta);
        const double dt=2*cap*(t>0?1:-1)-dn*etaPrime;
        return normalX?Jet{(n-eta)*(n-eta)+cap*cap-b*b,dn,dt,-2*etaPrime}
            :Jet{(n-eta)*(n-eta)+cap*cap-b*b,dt,dn,-2*etaPrime};
    };
}
void heightPlanePatches()
{
    for(bool normalX:{false,true}){
        auto mesh=makeMesh(32);const auto phase=HydroPhaseGeometry::build(*mesh,metricsFor(*mesh,false),heightPublication(*mesh,heightCapsule(normalX)));
        const auto pressure=HydroPressureLoad::build(phase.metrics().divergence());
        const auto start=std::chrono::steady_clock::now();
        const auto load=HydroBalancedCsfLoad::heightFunctionCapillary(*mesh,pressure,phase,.7,heightConfig());
        const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        heightForceAlgebra(*mesh,pressure,phase,load,.7);int straightCells=0;double peak=0;
        for(int p=0;p<phase.layout().pressureCount();++p){const auto center=mesh->cells(mesh->cellIndexToLabel[phase.layout().nativeCell(p)]).centerCoord();
            const double tangent=(normalX?center.y:center.x)-.5,normal=(normalX?center.x:center.y)-.5;
            const double k=load.heightFunctionCellCurvatures()[p];
            if(std::abs(tangent)<=1./64&&std::abs(normal)>.2&&std::isfinite(k)){
                peak=std::max(peak,std::abs(k));check(std::abs(k)<=1e-10*32,"HF closed straight patch acquired curvature");++straightCells;}
        }
        check(straightCells>=4,"HF straight patch normal direction missing");
        std::cout<<"HEIGHT_FUNCTION_PLANE_PATCH,normal_axis,"<<(normalX?'x':'y')<<",cells,"<<straightCells
            <<",curvature_max,"<<peak<<",factory_seconds,"<<seconds<<",PASS\n";
    }
}
struct HeightCircleError {double rms=0,maximum=0;};
HeightCircleError heightCircle(int n,Vector2D origin,int flip,bool provenance)
{
    auto mesh=heightMesh(n,origin);const Vector2D center=origin+Vector2D{.487,.521};constexpr double radius=.32,sigma=.7;
    const Profile f=[=](Vector2D p){const auto d=p-center;return Jet{flip*(d.x*d.x+d.y*d.y-radius*radius),flip*2*d.x,flip*2*d.y,0};};
    const HermiteGalsTimeLayer time{7,.002,.0022,Position::Midpoint};
    const auto phase=HydroPhaseGeometry::build(*mesh,metricsFor(*mesh,false),heightPublication(*mesh,f,time));
    const auto pressure=HydroPressureLoad::build(phase.metrics().divergence());const auto beforeAlpha=phase.liquidFractions();
    const auto start=std::chrono::steady_clock::now();
    const auto load=HydroBalancedCsfLoad::heightFunctionCapillary(*mesh,pressure,phase,sigma,heightConfig());
    const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    heightForceAlgebra(*mesh,pressure,phase,load,sigma);HeightCircleError error;long double squared=0;int count=0;
    for(double k:load.heightFunctionCellCurvatures())if(std::isfinite(k)){
        const double e=std::abs(k-flip/radius)*radius;error.maximum=std::max(error.maximum,e);squared+=e*e;++count;
        check(flip*k>0,"HF circle phase flip lost curvature sign");}
    check(count>0,"HF circle has no measured cell curvature");error.rms=std::sqrt(static_cast<double>(squared/count));
    check(error.maximum<(n==32?.08:.025),"HF independent analytic circle curvature accuracy");
    const auto& report=load.curvatureReport();check(report.heightFunctionCells==static_cast<std::size_t>(count)
        &&report.heightFunctionColumns>=3*report.heightFunctionCells&&report.heightFunctionCellQueries>0,"HF compact cell/column work count");
    check(report.heightFunctionXCells>0&&report.heightFunctionYCells>0
        &&report.heightFunctionXCells+report.heightFunctionYCells==report.heightFunctionCells,"HF circle normal X/Y direction coverage");
    const Reference incidence(*mesh,phase);long double radialWork=0,radialScale=0;int restoringFaces=0;
    for(int u=0;u<phase.layout().velocityCount();++u)if(incidence.minus[u]>=0&&incidence.plus[u]>=0&&load.supportReports()[u].curvatureDefined){
        const auto face=phase.layout().domain().face(*mesh,phase.layout().nativeFace(u));
        const double outward=phase.layout().normalAxis(u)==HydroDofLayout::NormalAxis::First?face.center.x-center.x:face.center.y-center.y;
        const long double term=load.momentumRhsLoad()[u]*outward;radialWork+=term;radialScale+=std::abs(term);++restoringFaces;
        check(term<=1e-13*std::max(1.L,radialScale),"HF circle force points away from interface restoring direction");
    }
    check(restoringFaces>0&&radialWork<0,"HF circle radial restoring force missing");
    check(phase.liquidFractions()==beforeAlpha,"HF mutated source alpha");
    // Analytic volumes are comparison only and never supply the factory.
    long double volumeError=0;
    for(int p=0;p<phase.layout().pressureCount();++p){const auto bounds=phase.publication().state().geometry().snapshot().nativeCellBounds(phase.layout().nativeCell(p));
        const auto v=circleCell(bounds[0],bounds[1],center,radius,false);
        const auto exact=flip>0?v:phase.metrics().cellVolumes()[p]-v;volumeError+=std::abs(phase.cellVolumes()[p].liquid-exact);}
    check(volumeError<=1e-10,"HF analytic reference/source circle geometry mismatch");
    if(provenance){
        const auto saved=load.momentumRhsLoad();
        const auto fresh=HydroPressureLoad::build(phase.metrics().divergence());
        rejects([&]{load.validateFor(*mesh,fresh,phase.geometryIdentity(),time);},"pressure owner");
        rejects([&]{load.validateFor(*mesh,pressure,HermiteGalsIdentity{},time);},"identity");
        rejects([&]{load.validateFor(*mesh,pressure,heightPublication(*mesh,f,time).state().identity(),time);},"identity");
        auto changed=time;changed.begin-=.0001;changed.end+=.0001;
        rejects([&]{load.validateFor(*mesh,pressure,phase.geometryIdentity(),changed);},"time");
        changed=time;changed.step++;rejects([&]{load.validateFor(*mesh,pressure,phase.geometryIdentity(),changed);},"time");
        changed=time;changed.position=Position::End;rejects([&]{load.validateFor(*mesh,pressure,phase.geometryIdentity(),changed);},"time");
        changed=time;changed.end=std::numeric_limits<double>::quiet_NaN();rejects([&]{load.validateFor(*mesh,pressure,phase.geometryIdentity(),changed);},"time");
        const auto wrong=HydroPressureLoad::build(metricsFor(*mesh,false).divergence());
        rejects([&]{HydroBalancedCsfLoad::heightFunctionCapillary(*mesh,wrong,phase,sigma,heightConfig());},"layout");
        DictionaryMesh::RefineBox box;box.min={.25,.25};box.max={.75,.75};box.level=1;mesh->applyRefineBoxes({box});
        rejects([&]{load.validateFor(*mesh,pressure,phase.geometryIdentity(),time);},"mesh");
        check(load.momentumRhsLoad()==saved,"HF failed provenance check mutated retained force");
    }
    std::cout<<"HEIGHT_FUNCTION_CIRCLE,n,"<<n<<",origin_x,"<<origin.x<<",origin_y,"<<origin.y<<",phase_sign,"<<flip
        <<",cells,"<<count<<",relative_rms,"<<error.rms<<",relative_max,"<<error.maximum<<",radial_force_work,"<<double(radialWork)
        <<",columns,"<<report.heightFunctionColumns<<",cell_queries,"<<report.heightFunctionCellQueries
        <<",alternate_directions,"<<report.heightFunctionAlternateDirections<<",normal_x_cells,"<<report.heightFunctionXCells<<",normal_y_cells,"<<report.heightFunctionYCells<<",factory_seconds,"<<seconds<<",PASS\n";
    return error;
}
void heightActualHighFrequency()
{
    auto mesh=makeMesh(64);const auto phase=HydroPhaseGeometry::build(*mesh,metricsFor(*mesh,false),heightPublication(*mesh,heightCapsule(false,true)));
    const auto pressure=HydroPressureLoad::build(phase.metrics().divergence());
    const auto load=HydroBalancedCsfLoad::heightFunctionCapillary(*mesh,pressure,phase,.7,heightConfig());
    heightForceAlgebra(*mesh,pressure,phase,load,.7);long double work=0,scale=0;int rows=0;
    const Reference incidence(*mesh,phase);
    for(int u=0;u<phase.layout().velocityCount();++u)if(phase.layout().normalAxis(u)==HydroDofLayout::NormalAxis::Second
        &&incidence.minus[u]>=0&&incidence.plus[u]>=0&&load.supportReports()[u].curvatureDefined){
        const auto face=phase.layout().domain().face(*mesh,phase.layout().nativeFace(u));const double t=face.center.x-.5;
        if(std::abs(t)>.08)continue;
        const double eta=.01/64*std::cos(32*double(pi)*t);const long double term=load.momentumRhsLoad()[u]*eta;
        work+=term;scale+=std::abs(term);++rows;
    }
    check(rows>=4&&scale>0&&work<-.05L*scale,"HF actual high-frequency force loses restoring correlation");
    std::cout<<"HEIGHT_FUNCTION_ACTUAL_MODE,theta_pi_over_2,rows,"<<rows<<",restoring_work,"<<double(work)
        <<",absolute_work,"<<double(scale)<<",PASS\n";
}
void heightUnsupportedScopes()
{
    for(int scope=0;scope<4;++scope){
        auto mesh=makeMesh(8,scope==1?1:0);const bool rz=scope==0;
        const auto pub=scope==3?heightPublication(*mesh,[](Vector2D){return Jet{1,0,0,0};},{},true):purePublication(*mesh,rz,1);
        const auto phase=HydroPhaseGeometry::build(*mesh,metricsFor(*mesh,rz,scope==2?1:0),pub);
        const auto pressure=HydroPressureLoad::build(phase.metrics().divergence());
        rejects([&]{HydroBalancedCsfLoad::heightFunctionCapillary(*mesh,pressure,phase,0,heightConfig());},
            scope==0?"XY":scope==1?"uniform":scope==2?"complete":"SingleQuartic");
    }
    auto mesh=makeMesh(8);Case plane;plane.pureFirstJet=true;plane.planar=true;plane.planeOffset=.43;
    const auto phase=HydroPhaseGeometry::build(*mesh,metricsFor(*mesh,false),published(*mesh,plane));
    const auto pressure=HydroPressureLoad::build(phase.metrics().divergence());
    rejects([&]{HydroBalancedCsfLoad::heightFunctionCapillary(*mesh,pressure,phase,.7,heightConfig());},"height");
    std::cout<<"HEIGHT_FUNCTION_UNSUPPORTED,rz_amr_partial_macro_wall,PASS\n";
}
// The increment force pairs the compatible first-jet generator with a fixed
// accepted GL8 rule. This is a residual work identity, not a whole-energy or
// adaptive-traction certificate. The reference below evaluates actual motion
// first and applies a forward D4 primary/scalar rate without the sparse load.
void heightFunctionIncrementOnly()
{
    constexpr int n=32;constexpr double sigma=.7;
    auto mesh=makeMesh(n);
    const Profile profile=[](Vector2D p){const double x=p.x-.487,y=p.y-.521;
        return Jet{x*x/(.31*.31)+y*y/(.26*.26)-1,2*x/(.31*.31),2*y/(.26*.26),0};};
    const auto seed=stateFor(*mesh,profile,false,{},[](int,Vector2D,Vector2D){return 0.;},{},HermiteGalsTransport::FirstJet);
    auto primarySeed=seed.primary();primarySeed.gradientUpdate=HermiteGalsGradientUpdate::CompatibleIncrement;
    HermiteGalsCandidate candidate(*mesh,{},primarySeed);const auto closure=closeHermiteGalsMixedAmr(*mesh,candidate);
    HermiteGalsPublicationConfig publicationConfig;publicationConfig.gradientUpdate=HermiteGalsGradientUpdate::CompatibleIncrement;
    const auto id=candidate.identity();
    const auto publication=publishHermiteGalsFirstJetGeometry(*mesh,std::move(candidate),closure,{id,centers(*mesh,profile)},publicationConfig);
    const auto phase=HydroPhaseGeometry::build(*mesh,metricsFor(*mesh,false),publication);const auto& layout=phase.layout();
    const auto pressure=HydroPressureLoad::build(phase.metrics().divergence());
    const auto geometry=HydroMomentumGeometry::build(*mesh,phase.metrics());
    std::vector<int> walls;std::vector<double> velocity(layout.velocityCount());
    for(int u=0;u<layout.velocityCount();++u){const auto face=layout.domain().face(*mesh,layout.nativeFace(u));
        if(face.negativeCell<0||face.positiveCell<0)walls.push_back(u);
        const double x=face.center.x,y=face.center.y;
        velocity[u]=layout.normalAxis(u)==HydroDofLayout::NormalAxis::First?
            std::sin(2*double(pi)*x)*std::pow(std::sin(double(pi)*y),2):
            -.8*std::pow(std::sin(double(pi)*x),2)*std::sin(2*double(pi)*y);
    }
    for(int u:walls)velocity[u]=0;
    const auto boundary=HydroAdvectorBoundary::impermeable(*mesh,geometry,walls);
    const auto reconstruction=HydroAdvectorReconstruction::prepare(*mesh,geometry,boundary);
    const auto motion=reconstruction.evaluate(*mesh,HydroNativeValues(layout,velocity,std::vector<double>(layout.pressureCount())),
        phase.geometryIdentity(),phase.timeLayer(),boundary);
    const auto mixed=HermiteGalsAmrClosureSensitivity::buildUniformFirstJet(*mesh,publication.state());
    auto config=heightConfig();config.forceMethod=HydroCsfForceMethod::HeightFunctionIncrement;
    const auto load=HydroBalancedCsfLoad::heightFunctionIncrementCapillary(*mesh,pressure,phase,sigma,motion,&mixed,config);
    load.validateFor(*mesh,pressure,phase.geometryIdentity(),phase.timeLayer(),motion,&mixed);
    check(load.discretization()==HydroCsfDiscretization::GeometricHeightFunctionIncrement,"increment force kind");
    const auto& scalar=publication.state().geometry().snapshot();const auto& primary=publication.state().primary();
    const auto vertexPoint=[&](int v){return mesh->verts(mesh->vertIndexToLabel[v]).coord();};
    std::map<std::pair<double,double>,long double> pointRates;
    const auto phiRate=[&](Vector2D q){int native=-1;
        for(int c:scalar.nativeCellCandidates(q,q))if(scalar.nativeCellActive(c)&&layout.pressureDof(c)>=0&&(native<0||c<native))native=c;
        check(native>=0,"forward increment rate outside actual source");
        const auto g=scalar.scalarJetNative(q);const auto u=motion.cellVelocityJet(layout.pressureDof(native),q).jet.value;
        return -static_cast<long double>(u.x)*g.dx-static_cast<long double>(u.y)*g.dy;};
    const auto retainPoint=[&](Vector2D q){pointRates.emplace(std::make_pair(q.x,q.y),phiRate(q));};
    for(const auto& v:primary.vertices)retainPoint(vertexPoint(v.vertex));
    for(const auto& e:primary.edges){const auto a=vertexPoint(e.vertices[0]),b=vertexPoint(e.vertices[1]);retainPoint(a+(b-a)*.5);}
    for(int c:primary.geometryCells)retainPoint(mesh->cells(mesh->cellIndexToLabel[c]).centerCoord());
    std::vector<double> xs,ys;
    for(const auto& entry:pointRates){xs.push_back(entry.first.first);ys.push_back(entry.first.second);}
    for(auto* axis:{&xs,&ys}){std::sort(axis->begin(),axis->end());axis->erase(std::unique(axis->begin(),axis->end()),axis->end());}
    check(xs.size()==2*n+1&&ys.size()==2*n+1&&pointRates.size()==xs.size()*ys.size(),"actual complete half-grid rate support");
    // Physical five-node Lagrange derivative, independently formed with
    // unscaled abscissae and target-value subtraction (including one-sided ends).
    const auto derivative=[&](Vector2D q,int axis){const auto& coordinates=axis==0?xs:ys;
        const double x=axis==0?q.x:q.y;const int target=static_cast<int>(std::lower_bound(coordinates.begin(),coordinates.end(),x)-coordinates.begin());
        check(target<static_cast<int>(coordinates.size())&&coordinates[target]==x,"D4 target is not an actual physical node");
        const int begin=std::clamp(target-2,0,static_cast<int>(coordinates.size())-5);
        long double result=0;const auto reference=pointRates.at({q.x,q.y});
        for(int j=0;j<5;++j)if(begin+j!=target){long double numerator=1,denominator=1;
            for(int k=0;k<5;++k)if(k!=j){denominator*=static_cast<long double>(coordinates[begin+j])-coordinates[begin+k];
                if(begin+k!=target)numerator*=static_cast<long double>(x)-coordinates[begin+k];}
            const auto key=axis==0?std::make_pair(coordinates[begin+j],q.y):std::make_pair(q.x,coordinates[begin+j]);
            result+=numerator/denominator*(pointRates.at(key)-reference);
        }return result;};
    std::map<int,std::array<long double,4>> vertexRates;std::vector<std::array<long double,2>> edgeRates(primary.edges.size());
    std::map<int,long double> centerRates;
    for(const auto& v:primary.vertices){const auto q=vertexPoint(v.vertex);
        vertexRates[v.vertex]={pointRates.at({q.x,q.y}),derivative(q,0),derivative(q,1),0};}
    for(std::size_t k=0;k<primary.edges.size();++k){const auto& e=primary.edges[k];
        const auto a=vertexPoint(e.vertices[0]),b=vertexPoint(e.vertices[1]),q=a+(b-a)*.5;
        edgeRates[k]={pointRates.at({q.x,q.y}),derivative(q,a.y==b.y?1:0)};}
    for(std::size_t k=0;k<primary.vertices.size();++k){long double rate=0;
        for(const auto& term:mixed.rows()[k])rate+=static_cast<long double>(term.weight)*
            (term.component==HermiteGalsPrimaryDerivative::EdgeNormal?edgeRates.at(term.primaryOrdinal)[1]:
                vertexRates.at(primary.vertices.at(term.primaryOrdinal).vertex)[term.component==HermiteGalsPrimaryDerivative::VertexDx?1:2]);
        vertexRates[primary.vertices[k].vertex][3]=rate;}
    for(int c:primary.geometryCells){const auto q=mesh->cells(mesh->cellIndexToLabel[c]).centerCoord();centerRates[c]=pointRates.at({q.x,q.y});}
    const auto normalRate=[&](const HermiteGalsSurfaceSample& sample){const auto chart=LevelSetHermiteGeometry::cellChartPoint(sample.frame);
        long double dot=0;using Kind=LAMP::Geometry::ActiveQuarticHermiteSensitivity::Kind;
        for(const auto& term:scalar.nativeQuarticDifferentialSensitivity(sample.frame.nativeCellIndex,chart,{1,0,0,0,0,0})){
            const long double value=term.kind==Kind::Vertex?vertexRates.at(term.index)[term.component]:
                term.kind==Kind::EdgeMidpoint?edgeRates.at(term.index)[term.component]:centerRates.at(term.index);
            dot+=term.coefficient*value;}
        const auto g=scalar.scalarJetAtNativeCellChart(sample.frame.nativeCellIndex,chart);
        return -dot/std::hypot(static_cast<long double>(g.dx),static_cast<long double>(g.dy));};
    long double area=0,deltaMoment=0;std::size_t nodes=0;
    for(std::size_t c=0;c<publication.topology().curves.size();++c)
        for(std::size_t g=0;g<publication.topology().curves[c].graphs().size();++g)
            for(const auto& q:publication.quadrature(c,g)){const int p=layout.pressureDof(q.sample.frame.nativeCellIndex);
                const double hf=load.heightFunctionCellCurvatures().at(p);check(std::isfinite(hf),"cached increment owner lacks HF base curvature");
                area+=q.physicalWeight;deltaMoment+=q.physicalWeight*(static_cast<long double>(q.sample.curvature)-hf);++nodes;}
    check(area>0&&nodes==load.curvatureReport().incrementNodes,"increment accepted-rule node coverage");
    const long double mean=deltaMoment/area;long double expectedResidual=0,powerScale=0;double peakRate=0;
    for(std::size_t c=0;c<publication.topology().curves.size();++c)
        for(std::size_t g=0;g<publication.topology().curves[c].graphs().size();++g)
            for(const auto& q:publication.quadrature(c,g)){
                const int p=layout.pressureDof(q.sample.frame.nativeCellIndex);const long double expected=normalRate(q.sample);
                long double sparse=0;for(const auto& term:load.incrementNormalRateAtSample(q.sample)){
                    check(term.first>=0&&term.first<layout.velocityCount(),"increment rate source DOF outside native layout");
                    sparse+=term.second*velocity[term.first];}
                const double defect=relative(double(sparse-expected),std::max(1.,std::abs(double(expected))));
                check(defect<=1e-10,"sparse increment rate differs from independent actual-motion D4 forward rate");peakRate=std::max(peakRate,defect);
                const long double delta=static_cast<long double>(q.sample.curvature)-load.heightFunctionCellCurvatures().at(p);
                const long double work=-sigma*q.physicalWeight*(delta-mean)*expected;
                expectedResidual+=work;powerScale+=std::abs(work);
            }
    const Reference incidence(*mesh,phase);auto potential=phase.liquidFractions();
    for(int p=0;p<layout.pressureCount();++p){close(load.heightFunctionIncrementMeanCurvatureDifference()[p],double(mean),
            1e-13*std::max(1.,std::abs(double(mean))),"independent same-component GL8 delta mean");potential[p]*=sigma*double(mean);}
    const auto reference=incidence.gradient(potential);
    const auto base=HydroBalancedCsfLoad::heightFunctionCapillary(*mesh,pressure,phase,sigma,heightConfig());
    long double residualPower=0,basePower=0,referencePower=0,totalPower=0;
    for(int u=0;u<layout.velocityCount();++u){const double r=load.heightFunctionIncrementResidualLoad()[u],a=load.heightFunctionIncrementReferenceLoad()[u];
        close(a,reference[u],1e-12*std::max(1.,std::abs(reference[u])),"increment mean pressure load uses independent native C incidence");
        close(load.momentumRhsLoad()[u],base.momentumRhsLoad()[u]+(a+r),1e-13*std::max(1.,std::abs(load.momentumRhsLoad()[u])),
            "actual HF base plus reference/residual native decomposition");
        check(!load.supportReports()[u].curvatureErrorEstimateAvailable,"increment force claims old adaptive traction certificate");
        residualPower+=static_cast<long double>(velocity[u])*r;basePower+=static_cast<long double>(velocity[u])*base.momentumRhsLoad()[u];
        referencePower+=static_cast<long double>(velocity[u])*a;totalPower+=static_cast<long double>(velocity[u])*load.momentumRhsLoad()[u];
        powerScale+=std::abs(static_cast<long double>(velocity[u])*r);
    }
    const double powerDefect=relative(double(residualPower-expectedResidual),double(std::max(1.L,powerScale)));
    check(powerDefect<=1e-10&&std::abs(residualPower)>1e-7,"increment residual actual-motion GL8 adjoint work");
    close(double(totalPower),double(basePower+referencePower+expectedResidual),1e-10*std::max(1.,std::abs(double(totalPower))),"total increment native work decomposition");
    const auto zero=HydroBalancedCsfLoad::heightFunctionIncrementCapillary(*mesh,pressure,phase,0,motion,nullptr,config);
    zero.validateFor(*mesh,pressure,phase.geometryIdentity(),phase.timeLayer(),motion,nullptr);
    check(zero.curvatureReport().projectionSkippedSigmaZero&&zero.curvatureReport().incrementNodes==0
        &&zero.curvatureReport().incrementRatePoints==0&&zero.curvatureReport().heightFunctionCells==0,"sigma-zero increment performed geometric/rate work");
    for(const auto* vector:{&zero.momentumRhsLoad(),&zero.heightFunctionIncrementReferenceLoad(),&zero.heightFunctionIncrementResidualLoad()})
        for(double value:*vector)check(value==0&&!std::signbit(value),"sigma-zero increment is not literal positive zero");
    rejects([&]{load.curvatureAverages();},"face curvature multiplier");
    rejects([&]{HydroBalancedCsfLoad::heightFunctionIncrementCapillary(*mesh,pressure,phase,sigma,motion,nullptr,config);},"mixed-closure owner");
    const auto otherMotion=reconstruction.evaluate(*mesh,HydroNativeValues(layout,velocity,std::vector<double>(layout.pressureCount())),
        phase.geometryIdentity(),phase.timeLayer(),boundary);
    rejects([&]{load.validateFor(*mesh,pressure,phase.geometryIdentity(),phase.timeLayer(),otherMotion,&mixed);},"direction owner");
    auto changed=phase.timeLayer();changed.step=1;changed.begin=0;changed.end=.125;changed.position=Position::End;
    rejects([&]{load.validateFor(*mesh,pressure,phase.geometryIdentity(),changed,motion,&mixed);},"full time layer");
    const auto oldPublication=heightPublication(*mesh,profile);
    const auto oldPhase=HydroPhaseGeometry::build(*mesh,phase.metrics(),oldPublication);
    rejects([&]{HydroBalancedCsfLoad::heightFunctionIncrementCapillary(*mesh,pressure,oldPhase,sigma,motion,&mixed,config);},"compatible-increment publication");
    std::cout<<"HEIGHT_FUNCTION_INCREMENT,nodes,"<<nodes<<",rate_relative,"<<peakRate<<",work_relative,"<<powerDefect
        <<",mean_delta,"<<double(mean)<<",residual_power,"<<double(residualPower)<<",PASS\n";
}

// Independent affine-section integral of the published PLIC plane. The
// endpoint reference never calls the production clipping/box/HF helpers.
struct LocalPlicEndpointReference
{
    struct Result {bool valid=false,opposite=false;double curvature=0;std::array<long double,3> heights{};};
    const HydroPhaseGeometry& phase;
    // Independent Euclidean construction on the three height positions. This
    // uses triangle area and circumradius, not the producer's secant/hypot
    // expression. Integral heights need not lie on the analytic source circle.
    static long double heightCircumcircle(const std::array<long double,3>& heights,long double spacing)
    {
        check(std::isfinite(spacing)&&spacing>0,"height-circle reference has invalid physical spacing");
        const std::array<std::array<long double,2>,3> p{{
            {-spacing,heights[0]-heights[1]},{0,0},{spacing,heights[2]-heights[1]}}};
        const auto edge=[](const auto& a,const auto& b){return std::hypot(b[0]-a[0],b[1]-a[1]);};
        const auto signedTwiceArea=(p[1][0]-p[0][0])*(p[2][1]-p[0][1])-(p[1][1]-p[0][1])*(p[2][0]-p[0][0]);
        if(signedTwiceArea==0)return 0;
        const auto area=std::abs(signedTwiceArea)/2;
        const auto radius=edge(p[0],p[1])*edge(p[1],p[2])*edge(p[2],p[0])/(4*area);
        check(std::isfinite(radius)&&radius>0,"height-circle reference has no finite geometric circumradius");
        return (signedTwiceArea<0?1.L:-1.L)/radius;
    }
    long double declaredCompleteCurvature(const Result& chart,long double spacing) const
    {
        const auto& owner=phase.vofPublication();
        return owner.localPlicResolution()==LAMP::Interface::GeometricVofLocalResolution::NativeInterfaceSupport&&!owner.domain().axisymmetric()
            ?heightCircumcircle(chart.heights,spacing):chart.curvature;
    }
    static long double planeArea(const LAMP::Interface::GeometricVofPlicCell& cell,Vector2D lo,Vector2D hi)
    {
        const long double capacity=static_cast<long double>(hi.x-lo.x)*(hi.y-lo.y);
        if(cell.measured[0]==0)return 0;if(cell.measured[1]==0)return capacity;
        const auto& b=cell.rectangle;const long double ox=static_cast<long double>(b.lower.x)+(static_cast<long double>(b.upper.x)-b.lower.x)/2;
        const long double oy=static_cast<long double>(b.lower.y)+(static_cast<long double>(b.upper.y)-b.lower.y)/2;
        const auto norm=std::hypot(static_cast<long double>(cell.normal.x),static_cast<long double>(cell.normal.y));
        // The publication stores its physical offset along the UNIT normal.
        const auto nx=cell.normal.x/norm,ny=cell.normal.y/norm;const bool vertical=std::abs(ny)>=std::abs(nx);
        const long double t0=vertical?static_cast<long double>(lo.x)-ox:static_cast<long double>(lo.y)-oy;
        const long double t1=vertical?static_cast<long double>(hi.x)-ox:static_cast<long double>(hi.y)-oy;
        const long double v0=vertical?static_cast<long double>(lo.y)-oy:static_cast<long double>(lo.x)-ox;
        const long double v1=vertical?static_cast<long double>(hi.y)-oy:static_cast<long double>(hi.x)-ox;
        const auto major=vertical?ny:nx,minor=vertical?nx:ny;std::vector<long double> cuts{t0,t1};
        if(minor!=0)for(auto edge:{v0,v1}){const auto crossing=(cell.offset-major*edge)/minor;
            if(crossing>t0&&crossing<t1)cuts.push_back(crossing);}
        std::sort(cuts.begin(),cuts.end());
        const auto thickness=[&](long double t){const auto intercept=(cell.offset-minor*t)/major;
            return std::max(0.L,std::min(v1-v0,major>0?intercept-v0:v1-intercept));};
        long double area=0,correction=0;
        for(std::size_t i=1;i<cuts.size();++i){const auto term=(cuts[i]-cuts[i-1])*(thickness(cuts[i-1])+thickness(cuts[i]))/2-correction;
            const auto next=area+term;correction=(next-area)-term;area=next;}return area;
    }
    static Vector2D segmentMidpoint(const LAMP::Interface::GeometricVofPlicCell& cell)
    {
        check(cell.resolved&&cell.measured[0]>0&&cell.measured[1]>0,"reference segment requires actual resolved mixed geometry");
        const auto& b=cell.rectangle;const auto center=(b.lower+b.upper)*.5;
        const auto norm=std::hypot(static_cast<long double>(cell.normal.x),static_cast<long double>(cell.normal.y));
        const long double nx=cell.normal.x/norm,ny=cell.normal.y/norm;
        const std::array<Vector2D,4> corners{{b.lower,{b.upper.x,b.lower.y},b.upper,{b.lower.x,b.upper.y}}};
        std::vector<Vector2D> endpoints;
        const auto add=[&](Vector2D p){for(const auto q:endpoints)if(p.x==q.x&&p.y==q.y)return;endpoints.push_back(p);};
        const auto phi=[&](Vector2D p){return nx*(static_cast<long double>(p.x)-center.x)+ny*(static_cast<long double>(p.y)-center.y)-cell.offset;};
        for(int edge=0;edge<4;++edge){const auto a=corners[edge],z=corners[(edge+1)%4];const auto fa=phi(a),fz=phi(z);
            if(fa==0)add(a);
            if((fa<0&&fz>0)||(fa>0&&fz<0)){const auto t=fa/(fa-fz);
                add({static_cast<double>(a.x+t*(static_cast<long double>(z.x)-a.x)),
                    static_cast<double>(a.y+t*(static_cast<long double>(z.y)-a.y))});}}
        check(endpoints.size()==2,"reference mixed plane must have exactly two actual rectangle-edge intersections");
        return (endpoints[0]+endpoints[1])*.5;
    }
    Result estimate(int pressure,int axis,int lower=-4,int upper=4,const Vector2D* anchor=nullptr,double fineSpacing=0) const
    {
        Result result;const auto& cells=phase.vofPublication().plicCells();const auto& source=cells[phase.layout().nativeCell(pressure)];const auto& b=source.rectangle;
        const double component=axis?source.normal.y:source.normal.x;if(component==0)return result;
        const int outward=component>0?1:-1,tangent=1-axis;
        const double dh=fineSpacing>0?fineSpacing:(axis?b.upper.y-b.lower.y:b.upper.x-b.lower.x),
            spacing=fineSpacing>0?fineSpacing:(tangent?b.upper.y-b.lower.y:b.upper.x-b.lower.x);
        bool ends=true,monotone=true,coveredAll=true;const int count=upper-lower+1;
        for(int column=-1;column<=1;++column){std::vector<long double> fractions(count);long double compensation=0;
            for(int line=lower;line<=upper;++line){Vector2D center=anchor?*anchor:(b.lower+b.upper)*.5;
                (axis?center.y:center.x)+=line*dh;(tangent?center.y:center.x)+=column*spacing;
                const Vector2D half=axis?Vector2D{spacing/2,dh/2}:Vector2D{dh/2,spacing/2};const auto lo=center-half,hi=center+half;
                const long double capacity=static_cast<long double>(hi.x-lo.x)*(hi.y-lo.y);long double liquid=0,covered=0,correction=0;
                for(const auto& cell:cells){const auto& cb=cell.rectangle;
                    const Vector2D a{std::max(lo.x,cb.lower.x),std::max(lo.y,cb.lower.y)},z{std::min(hi.x,cb.upper.x),std::min(hi.y,cb.upper.y)};
                    if(a.x>=z.x||a.y>=z.y)continue;
                    if(cell.measured[0]>0&&cell.measured[1]>0&&(!cell.resolved||outward*(axis?cell.normal.y:cell.normal.x)<=0))result.opposite=true;
                    covered+=static_cast<long double>(z.x-a.x)*(z.y-a.y);const auto term=planeArea(cell,a,z)-correction,next=liquid+term;
                    correction=(next-liquid)-term;liquid=next;}
                coveredAll=coveredAll&&std::abs(covered-capacity)<=1e-12L*capacity;
                fractions[line-lower]=liquid/capacity;
            }
            ends=ends&&fractions[outward>0?0:count-1]>=1-1e-6L&&fractions[outward>0?count-1:0]<=1e-6L;
            long double previous=1;
            for(int step=0;step<count;++step){const auto fraction=fractions[outward>0?step:count-1-step];
                monotone=monotone&&fraction<=previous+1e-6L;previous=fraction;
                const auto term=fraction*dh-compensation,next=result.heights[column+1]+term;
                compensation=(next-result.heights[column+1])-term;result.heights[column+1]=next;}
        }
        // All columns use the SAME lower/upper physical box planes. Independent
        // stopping at a different pure endpoint per column would change these differences.
        const auto first=(result.heights[2]-result.heights[0])/(2*spacing);
        const auto second=(result.heights[2]-2*result.heights[1]+result.heights[0])/(spacing*spacing),metric=std::hypot(1.L,first);
        result.curvature=static_cast<double>(-second/(metric*metric*metric));
        result.valid=coveredAll&&ends&&monotone&&!result.opposite&&std::isfinite(result.curvature);return result;
    }
    long double arithmeticBound(int pressure) const
    {const auto& b=phase.vofPublication().plicCells()[phase.layout().nativeCell(pressure)].rectangle;
        return 4096.L*std::numeric_limits<double>::epsilon()/std::min(b.upper.x-b.lower.x,b.upper.y-b.lower.y);}
};

void localPlicHeightFunction()
{
    const auto localBegin=std::chrono::steady_clock::now();
    auto mesh=makeMesh(32,1);const Vector2D center{.487,.521};constexpr double radius=.32,sigma=.7;
    const auto metrics=metricsFor(*mesh,false);const auto& layout=metrics.layout();
    check(!layout.domain().uniformLeaves(),"local PLIC HF witness retains actual one-level 2:1 AMR");
    std::vector<double> alpha(mesh->cellIndexNum());
    for(int ci=0;ci<mesh->cellIndexNum();++ci){const auto& c=mesh->cells(mesh->cellIndexToLabel[ci]);
        const auto lo=mesh->verts(c.cellToVert[0]).coord(),hi=mesh->verts(c.cellToVert[2]).coord();
        alpha[ci]=static_cast<double>(circleCell(lo,hi,center,radius,false)/((hi.x-lo.x)*(hi.y-lo.y)));}
    // One named positive stock below this PLIC plane's representable geometry.
    // The selected existing coarse leaf and its neighbors are interior pure gas.
    std::vector<std::array<long double,2>> stocks(alpha.size());int microCell=-1;
    for(int ci=0;ci<mesh->cellIndexNum();++ci){const auto& c=mesh->cells(mesh->cellIndexToLabel[ci]);
        const auto lo=mesh->verts(c.cellToVert[0]).coord();const auto capacity=metrics.cellVolumes()[layout.pressureDof(ci)];
        stocks[ci]={static_cast<long double>(alpha[ci])*capacity,(1.L-static_cast<long double>(alpha[ci]))*capacity};
        if(lo.x==.0625&&lo.y==.0625){check(alpha[ci]==0,"named stock-only micro leaf was not originally pure gas");
            microCell=ci;stocks[ci]={capacity*1e-50L,capacity};alpha[ci]=static_cast<double>(stocks[ci][0]/capacity);}}
    check(microCell>=0,"existing local PLIC AMR scene lacks named interior pure-gas leaf");
    const auto retainedStocks=stocks;
    const HermiteGalsTimeLayer clock{7,.002,.0022,Position::Midpoint};
    const auto publication=LAMP::Interface::GeometricVofPublication::prepareLocalPlicInventory(*mesh,layout.domain(),clock,stocks,{});
    const auto phase=HydroPhaseGeometry::build(*mesh,metrics,publication);const auto pressure=HydroPressureLoad::build(metrics.divergence());
    const auto load=HydroBalancedCsfLoad::localPlicHeightFunctionCapillary(*mesh,pressure,phase,sigma,heightConfig());
    load.validateVofFor(*mesh,pressure,publication.identity(),clock);
    check(load.vofIdentity()==publication.identity()&&load.phaseGeometry().model()==HydroPhaseGeometry::Model::GeometricVof&&
        &load.phaseGeometry().cellVolumes()==&phase.cellVolumes(),"local PLIC HF retains its actual typed VOF midpoint owner");
    const Reference incidence(*mesh,phase);std::vector<double> actualColor(layout.pressureCount());std::size_t actualMixed=0;
    for(int p=0;p<layout.pressureCount();++p){const auto volumes=phase.cellVolumes()[p];
        actualColor[p]=volumes.liquid==0?0:volumes.gas==0?1:static_cast<double>(volumes.liquid/metrics.cellVolumes()[p]);
        actualMixed+=volumes.liquid>0&&volumes.gas>0;}
    const auto independentGradient=incidence.gradient(actualColor);
    const auto authoritativeGradient=pressure.apply(*mesh,HydroNativeValues(layout,std::vector<double>(layout.velocityCount()),actualColor));
    const auto inventoryGradient=incidence.gradient(phase.liquidFractions());
    const auto microPressure=layout.pressureDof(microCell);const auto& micro=publication.plicCells()[microCell];
    check(micro.inventory[0]>0&&micro.inventory[1]>0&&micro.measured[0]==0&&micro.measured[1]==metrics.cellVolumes()[microPressure]&&!micro.resolved,
        "named positive micro stock must retain its explicitly unresolved pure actual PLIC geometry");
    check(actualColor[microPressure]==0&&phase.liquidFractions()[microPressure]>0
        &&!std::isfinite(load.heightFunctionCellCurvatures()[microPressure]),
        "actual pure micro leaf must retain positive stock alpha without fabricated interface curvature");
    std::size_t cells=0,loaded=0,microRows=0;double maximum=0;long double sum=0;
    for(double k:load.heightFunctionCellCurvatures())if(std::isfinite(k)){++cells;maximum=std::max(maximum,std::abs(k));sum+=k;}
    check(cells==actualMixed&&cells+1==publication.report().mixedCells&&cells>0&&maximum>0&&sum>0,
        "local PLIC circle has genuinely measured finite nonzero HF curvature on every actual mixed AMR leaf");
    for(int u=0;u<layout.velocityCount();++u){
        close(load.colorGradient()[u],independentGradient[u],1e-12*std::max(1.,std::abs(independentGradient[u])),
            "local PLIC HF uses independent actual AMR incidence C measured color");
        close(load.colorGradient()[u],authoritativeGradient[u],1e-13*std::max(1.,std::abs(authoritativeGradient[u])),
            "local PLIC HF measured color uses authoritative retained pressure C");
        const auto& support=load.supportReports()[u];
        if(incidence.minus[u]==microPressure||incidence.plus[u]==microPressure){
            check(incidence.minus[u]>=0&&incidence.plus[u]>=0&&inventoryGradient[u]!=0,
                "named stock-only leaf requires a real interior nonzero inventory-color gradient witness");
            check(independentGradient[u]==0&&authoritativeGradient[u]==0&&load.colorGradient()[u]==0
                &&load.momentumRhsLoad()[u]==0&&!support.curvatureDefined,
                "actual pure-color row must remain exact zero without requiring a stock-only curvature");++microRows;}
        check(!support.curvatureErrorEstimateAvailable&&load.quadraturePanels(u).empty(),"local PLIC HF does not fabricate quadrature curvature certificates");
        if(support.curvatureDefined){long double k=0;int count=0;
            for(int p:{incidence.minus[u],incidence.plus[u]})if(p>=0&&std::isfinite(load.heightFunctionCellCurvatures()[p])){k+=load.heightFunctionCellCurvatures()[p];++count;}
            check(count>0,"local PLIC force only consumes actual adjacent mixed-cell heights");
            close(load.curvatureAverages()[u],k/count,1e-12*std::max(1.L,std::abs(k/count)),"local PLIC actual adjacent HF mean");
            close(load.momentumRhsLoad()[u],sigma*(k/count)*independentGradient[u],
                1e-12*std::max(1.L,std::abs(sigma*(k/count)*independentGradient[u])),"local PLIC actual force is sigma kappa C alpha");++loaded;
        }else check(load.momentumRhsLoad()[u]==0,"local PLIC has no unsupported force coefficient");
    }
    check(microRows==4,"named interior stock-only micro leaf must exercise its four actual native rows");
    check(loaded>0&&load.curvatureReport().amrFaces>0&&load.curvatureReport().heightFunctionCells==cells,
        "local PLIC HF consumes actual hanging subfaces and local height columns");
    std::vector<double> potential(layout.pressureCount());for(int p=0;p<layout.pressureCount();++p)potential[p]=sigma/radius*actualColor[p];
    const auto constantMode=pressure.apply(*mesh,HydroNativeValues(layout,std::vector<double>(layout.velocityCount()),potential));
    for(int u=0;u<layout.velocityCount();++u)close(constantMode[u],sigma/radius*independentGradient[u],
        1e-12*std::max(1.,std::abs(constantMode[u])),"constant-curvature pressure mode uses the same actual AMR C measured color");
    const auto zero=HydroBalancedCsfLoad::localPlicHeightFunctionCapillary(*mesh,pressure,phase,0,heightConfig());
    check(zero.curvatureReport().projectionSkippedSigmaZero&&zero.curvatureReport().heightFunctionCells==0,
        "zero sigma skips local PLIC box/height curvature queries");
    for(double f:zero.momentumRhsLoad())check(f==0,"zero sigma local PLIC force is exact zero");
    const auto foreign=LAMP::Interface::GeometricVofPublication::prepareLocalPlic(*mesh,layout.domain(),clock,alpha);
    rejects([&]{load.validateVofFor(*mesh,pressure,foreign.identity(),clock);},"owner");
    const auto otherPressure=HydroPressureLoad::build(metrics.divergence());
    rejects([&]{load.validateVofFor(*mesh,otherPressure,publication.identity(),clock);},"pressure owner");
    auto otherClock=clock;++otherClock.step;rejects([&]{load.validateVofFor(*mesh,pressure,publication.identity(),otherClock);},"clock");
    check(publication.alpha()==alpha,"local PLIC HF never mutates authoritative P0 stocks");
    for(std::size_t ci=0;ci<stocks.size();++ci)check(publication.plicCells()[ci].inventory==retainedStocks[ci],
        "local PLIC HF removed or changed independently retained phase stock");
    check(load.curvatureReport().heightFunctionExtendedCells==0&&load.curvatureReport().heightFunctionExtendedSlices==0,
        "existing AMR circle must retain its successful original nine-slice paths without extensions");
    LocalPlicEndpointReference originalReference{phase};bool originalPriority=false;
    for(int p=0;p<layout.pressureCount()&&!originalPriority;++p)if(std::isfinite(load.heightFunctionCellCurvatures()[p])){
        const auto n=publication.plicCells()[layout.nativeCell(p)].normal;
        if(std::abs(n.x)<.3||std::abs(n.y)<.3)continue;
        const auto x=originalReference.estimate(p,0),y=originalReference.estimate(p,1);const auto budget=originalReference.arithmeticBound(p);
        if(x.valid&&y.valid&&std::abs(x.curvature-y.curvature)>16*budget){
            const auto& preferred=std::abs(n.x)>=std::abs(n.y)?x:y;
            close(load.heightFunctionCellCurvatures()[p],preferred.curvature,budget,
                "original preferred nine-slice estimate wins before any endpoint extension");originalPriority=true;}}
    check(originalPriority,"existing shifted AMR circle lacks an independent original-axis priority witness");
    check(load.curvatureReport().heightFunctionRefinedCells==0&&load.curvatureReport().heightFunctionRefinedSlices==0,
        "original successful native nine-slice paths must not prepare refined frames");
    check(load.curvatureReport().heightFunctionFitCells==0&&load.curvatureReport().heightFunctionFitSegments==0,
        "original successful native nine-slice paths must not use the quadratic fit");

    // Reuse the same mesh, metrics, circle-area input and public publication.
    // R=7.15h leaves a small diagonal source at (-5.5h,-5.5h): its
    // outer transverse strip reaches 7h. The line4 interval [-2h,-h]
    // is not pure inside, while line5 [-h,0] is pure inside. Thus BOTH
    // original axes lack an end, and exactly one shared ring repairs it.
    double fineH=1;for(const auto& cell:publication.plicCells())fineH=std::min(fineH,cell.rectangle.upper.x-cell.rectangle.lower.x);
    const Vector2D smallCenter{.5,.5};
    const auto smallCircle=[&](double r){std::vector<double> input(alpha.size());
        for(int ci=0;ci<mesh->cellIndexNum();++ci){const auto& c=mesh->cells(mesh->cellIndexToLabel[ci]);
            const auto lo=mesh->verts(c.cellToVert[0]).coord(),hi=mesh->verts(c.cellToVert[2]).coord();
            input[ci]=static_cast<double>(circleCell(lo,hi,smallCenter,r,false)/((hi.x-lo.x)*(hi.y-lo.y)));}
        return LAMP::Interface::GeometricVofPublication::prepareLocalPlic(*mesh,layout.domain(),clock,input);};
    const auto namedDiagonal=[&](const auto& source,double distance){int found=-1;
        for(int p=0;p<layout.pressureCount();++p){const auto& b=source.plicCells()[layout.nativeCell(p)].rectangle;const auto c=(b.lower+b.upper)*.5;
            if(c.x==smallCenter.x-distance*fineH&&c.y==smallCenter.y-distance*fineH)found=p;}
        check(found>=0,"same AMR mesh lacks the named fine diagonal source");return found;};
    const auto extendedPublication=smallCircle(7.15*fineH);const auto extendedPhase=HydroPhaseGeometry::build(*mesh,metrics,extendedPublication);
    const auto extended=HydroBalancedCsfLoad::localPlicHeightFunctionCapillary(*mesh,pressure,extendedPhase,sigma,heightConfig());
    const int diagonal=namedDiagonal(extendedPublication,5.5);const auto& diagonalCell=extendedPublication.plicCells()[layout.nativeCell(diagonal)];
    check(diagonalCell.measured[0]>0&&diagonalCell.measured[1]>0,"named one-ring source must have actual mixed PLIC geometry");
    LocalPlicEndpointReference extendedReference{extendedPhase};
    const auto oldX=extendedReference.estimate(diagonal,0),oldY=extendedReference.estimate(diagonal,1);
    const auto ringX=extendedReference.estimate(diagonal,0,-4,5),ringY=extendedReference.estimate(diagonal,1,-4,5);
    check(!oldX.valid&&!oldY.valid&&ringX.valid&&ringY.valid,
        "named actual published circle must require exactly one endpoint ring after BOTH original axes fail");
    const int preferred=std::abs(diagonalCell.normal.x)>=std::abs(diagonalCell.normal.y)?0:1;
    const auto& expected=preferred==0?ringX:ringY;
    close(extended.heightFunctionCellCurvatures()[diagonal],expected.curvature,extendedReference.arithmeticBound(diagonal),
        "one-ring curvature uses independent actual-plane columns with a common physical origin");
    check(extended.curvatureReport().heightFunctionExtendedCells>0&&extended.curvatureReport().heightFunctionExtendedSlices>=3&&
        extended.curvatureReport().heightFunctionExtendedCells<=extended.curvatureReport().heightFunctionCells,
        "actual endpoint extension work and adopted cell counts must be reported");
    check(extended.curvatureReport().heightFunctionRefinedCells==0&&extended.curvatureReport().heightFunctionRefinedSlices==0,
        "original successful CT endpoint rings must retain priority over refined frames");
    check(extended.curvatureReport().heightFunctionFitCells==0&&extended.curvatureReport().heightFunctionFitSegments==0,
        "original successful CT endpoint rings must not use the quadratic fit");
    // A smaller circle places the opposite interface inside these columns.
    // HF remains unsupported: the new declared fit is a separate last resort,
    // and does not change the independent original HF topology qualification.
    const auto narrowPublication=smallCircle(4.8*fineH);const auto narrowPhase=HydroPhaseGeometry::build(*mesh,metrics,narrowPublication);
    const int narrowDiagonal=namedDiagonal(narrowPublication,3.5);LocalPlicEndpointReference narrowReference{narrowPhase};
    const auto narrowX=narrowReference.estimate(narrowDiagonal,0),narrowY=narrowReference.estimate(narrowDiagonal,1);
    check(narrowX.opposite&&narrowY.opposite&&!narrowX.valid&&!narrowY.valid,
        "named narrow-circle columns must expose the actual opposite branch");
    const auto fitted=HydroBalancedCsfLoad::localPlicHeightFunctionCapillary(*mesh,pressure,narrowPhase,sigma,heightConfig());
    fitted.validateVofFor(*mesh,pressure,narrowPublication.identity(),clock);
    const auto& fittedReport=fitted.curvatureReport();const auto& fittedSource=narrowPublication.plicCells()[layout.nativeCell(narrowDiagonal)];
    check(fittedSource.rectangle.upper.x-fittedSource.rectangle.lower.x==fineH&&
        fittedSource.rectangle.upper.y-fittedSource.rectangle.lower.y==fineH,
        "independently selected fit witness must be a native fine leaf, ineligible for CU refinement");
    // Both original charts contain an actual opposite branch, so CT cannot
    // repair them by adding slices; the finest source has no CU path. This
    // independently locates a cell that the successful factory MUST fit.
    // Its oracle is the analytic circle +1/R, not a second copy of the LS fit.
    const double analyticNarrow=1/(4.8*fineH),fittedKappa=fitted.heightFunctionCellCurvatures()[narrowDiagonal];
    const double fittedRelative=std::abs(fittedKappa-analyticNarrow)/analyticNarrow;
    check(std::isfinite(fittedKappa)&&fittedKappa>0&&fittedRelative<.35,
        "declared low-order PLIC fit must have the circle's positive curvature and less than 35 percent relative error");
    std::size_t finiteNarrow=0;for(double k:fitted.heightFunctionCellCurvatures())finiteNarrow+=std::isfinite(k);
    check(fittedReport.heightFunctionFitCells>0&&fittedReport.heightFunctionFitSegments>=3*fittedReport.heightFunctionFitCells&&
        fittedReport.heightFunctionCells+fittedReport.heightFunctionFitCells==finiteNarrow&&
        fittedReport.heightFunctionXCells+fittedReport.heightFunctionYCells==fittedReport.heightFunctionCells,
        "HF and adopted-fit cell counts must be separate and cover the actual finite curvature field");
    check(fitted.vofIdentity()==narrowPublication.identity()&&&fitted.phaseGeometry().cellVolumes()==&narrowPhase.cellVolumes(),
        "quadratic fit must retain the actual same-clock VOF geometry owner");
    std::cout<<"LOCAL_PLIC_HEIGHT_FIT,witness_native,"<<layout.nativeCell(narrowDiagonal)
        <<",hf_opposite_charts,2,native_fine,1,kappa,"<<fittedKappa<<",analytic_circle_kappa,"<<analyticNarrow
        <<",relative_error,"<<fittedRelative<<",fit_cells,"<<fittedReport.heightFunctionFitCells
        <<",fit_segments,"<<fittedReport.heightFunctionFitSegments<<",PASS\n";

    // One isolated actual mixed segment, plus a pure-liquid left neighbor,
    // gives a genuine nonzero alpha gradient and resolved PLIC orientation.
    // Every other leaf is pure gas: there is exactly ONE actual segment,
    // neither a left liquid HF end nor three native fit segments.
    std::vector<double> isolatedAlpha(alpha.size());int isolatedCell=-1,leftPure=-1;
    for(int ci=0;ci<mesh->cellIndexNum();++ci){const auto& b=publication.plicCells()[ci].rectangle;
        if(b.lower.x==.5&&b.lower.y==.5&&b.upper.x-b.lower.x==fineH)isolatedCell=ci;
        if(b.upper.x==.5&&b.lower.y==.5&&b.upper.x-b.lower.x==fineH)leftPure=ci;}
    check(isolatedCell>=0&&leftPure>=0,"same fine AMR region lacks the named isolated segment and pure neighbor");
    isolatedAlpha[isolatedCell]=.5;isolatedAlpha[leftPure]=1;
    const auto isolatedPublication=LAMP::Interface::GeometricVofPublication::prepareLocalPlic(*mesh,layout.domain(),clock,isolatedAlpha);
    std::size_t isolatedSegments=0;for(const auto& c:isolatedPublication.plicCells())isolatedSegments+=c.measured[0]>0&&c.measured[1]>0;
    const auto& isolatedSource=isolatedPublication.plicCells()[isolatedCell];
    check(isolatedSegments==1&&isolatedSource.resolved&&isolatedSource.measured[0]>0&&isolatedSource.measured[1]>0,
        "insufficient-fit witness must contain one genuinely resolved actual segment, not an unresolved or stock-only interface");
    const auto isolatedPhase=HydroPhaseGeometry::build(*mesh,metrics,isolatedPublication);
    rejects([&]{HydroBalancedCsfLoad::localPlicHeightFunctionCapillary(*mesh,pressure,isolatedPhase,sigma,heightConfig());},"segments=1");
    std::cout<<"LOCAL_PLIC_HEIGHT_ENDPOINTS,old9_priority,"<<originalPriority<<",one_ring_diagonal,"<<layout.nativeCell(diagonal)
        <<",extended_cells,"<<extended.curvatureReport().heightFunctionExtendedCells
        <<",extended_slices,"<<extended.curvatureReport().heightFunctionExtendedSlices<<",hf_opposite_branch_rejected,1,PASS\n";
    // One translated version of the existing R=7.15h analytic circle places
    // a small actual cut at a COARSE leaf's upper-right corner on x=.25.
    // Its native nine-slice width is 18 fine h and already contains the
    // opposite branch, so adding CT end rings cannot repair either axis.
    const Vector2D coarseCenter{.25+5*fineH,.375+5*fineH};
    std::vector<double> coarseInput(alpha.size());
    for(int ci=0;ci<mesh->cellIndexNum();++ci){const auto& c=mesh->cells(mesh->cellIndexToLabel[ci]);
        const auto lo=mesh->verts(c.cellToVert[0]).coord(),hi=mesh->verts(c.cellToVert[2]).coord();
        coarseInput[ci]=static_cast<double>(circleCell(lo,hi,coarseCenter,7.15*fineH,false)/((hi.x-lo.x)*(hi.y-lo.y)));}
    const auto coarsePublication=LAMP::Interface::GeometricVofPublication::prepareLocalPlic(*mesh,layout.domain(),clock,coarseInput);
    const auto coarsePhase=HydroPhaseGeometry::build(*mesh,metrics,coarsePublication);
    int coarseCorner=-1;for(int p=0;p<layout.pressureCount();++p){const auto& b=coarsePublication.plicCells()[layout.nativeCell(p)].rectangle;
        if(b.upper.x==.25&&b.upper.y==.375&&b.upper.x-b.lower.x==2*fineH&&b.upper.y-b.lower.y==2*fineH)coarseCorner=p;}
    check(coarseCorner>=0,"same true 2:1 AMR mesh lacks the named coarse corner leaf");
    const auto& coarseCell=coarsePublication.plicCells()[layout.nativeCell(coarseCorner)];
    check(coarseCell.resolved&&coarseCell.measured[0]>0&&coarseCell.measured[1]>0,"named coarse corner must retain an actual resolved cut");
    bool actualHangingNeighbor=false;
    for(int u=0;u<layout.velocityCount();++u){const auto face=layout.domain().face(*mesh,layout.nativeFace(u));
        const int native=layout.nativeCell(coarseCorner);
        const int neighbor=face.negativeCell==native?face.positiveCell:face.positiveCell==native?face.negativeCell:-1;
        if(neighbor<0)continue;const auto& b=coarsePublication.plicCells()[neighbor].rectangle;
        actualHangingNeighbor=actualHangingNeighbor||(b.upper.x-b.lower.x==fineH&&b.upper.y-b.lower.y==fineH);}
    check(actualHangingNeighbor,"named coarse source must share an actual native subface with a fine leaf");
    LocalPlicEndpointReference coarseReference{coarsePhase};
    const auto coarseX=coarseReference.estimate(coarseCorner,0),coarseY=coarseReference.estimate(coarseCorner,1);
    check(!coarseX.valid&&!coarseY.valid&&coarseX.opposite&&coarseY.opposite,
        "both native nine-slice frames already contain opposite branches which bounded end rings cannot remove");
    const auto segment=LocalPlicEndpointReference::segmentMidpoint(coarseCell);
    const auto refinedX=coarseReference.estimate(coarseCorner,0,-4,4,&segment,fineH);
    const auto refinedY=coarseReference.estimate(coarseCorner,1,-4,4,&segment,fineH);
    const int refinedFirst=std::abs(coarseCell.normal.x)>=std::abs(coarseCell.normal.y)?0:1;
    const auto& refinedPreferred=refinedFirst==0?refinedX:refinedY;
    const auto& refinedExpected=refinedPreferred.valid?refinedPreferred:(refinedFirst==0?refinedY:refinedX);
    check(refinedExpected.valid,"named coarse actual-segment anchor must have a valid fixed fine nine-slice chart");
    const auto refined=HydroBalancedCsfLoad::localPlicHeightFunctionCapillary(*mesh,pressure,coarsePhase,sigma,heightConfig());
    refined.validateVofFor(*mesh,pressure,coarsePublication.identity(),clock);
    close(refined.heightFunctionCellCurvatures()[coarseCorner],refinedExpected.curvature,
        4096.L*std::numeric_limits<double>::epsilon()/fineH,
        "actual coarse cut uses independent same-owner fine columns at the segment midpoint with their common physical origin");
    check(refined.curvatureReport().heightFunctionRefinedCells>0&&refined.curvatureReport().heightFunctionRefinedSlices>=27&&
        refined.curvatureReport().heightFunctionRefinedCells<=refined.curvatureReport().heightFunctionCells,
        "actual adopted fine-frame cells and measured fine virtual slices must be reported");
    check(refined.curvatureReport().heightFunctionFitCells==0&&refined.curvatureReport().heightFunctionFitSegments==0,
        "original successful CU coarse-cut frames must not use the quadratic fit");
    std::cout<<"LOCAL_PLIC_HEIGHT_REFINED,coarse_corner,"<<layout.nativeCell(coarseCorner)
        <<",segment_x,"<<segment.x<<",segment_y,"<<segment.y
        <<",native_opposite_branches,2,refined_cells,"<<refined.curvatureReport().heightFunctionRefinedCells
        <<",refined_slices,"<<refined.curvatureReport().heightFunctionRefinedSlices<<",PASS\n";
    // One CLOSED 45-degree interface, fully inside the existing fine region.
    // Long sides |dx+dy|=2h stay straight throughout the fixed 5h source box;
    // their closing corners |dx-dy|=12h are far outside it. Exact half-plane
    // clipping initializes alpha, without touching a physical boundary.
    const Vector2D diamondCenter{.5+.5*fineH,.5+.5*fineH};
    std::vector<double> diamondAlpha(alpha.size());
    const std::array<std::array<long double,3>,4> diamondPlanes{{
        {1,1,2*fineH},{-1,-1,2*fineH},{1,-1,12*fineH},{-1,1,12*fineH}}};
    for(int ci=0;ci<mesh->cellIndexNum();++ci){const auto& b=publication.plicCells()[ci].rectangle;
        using P=std::array<long double,2>;
        const long double x0=static_cast<long double>(b.lower.x)-diamondCenter.x,x1=static_cast<long double>(b.upper.x)-diamondCenter.x;
        const long double y0=static_cast<long double>(b.lower.y)-diamondCenter.y,y1=static_cast<long double>(b.upper.y)-diamondCenter.y;
        std::vector<P> polygon{{x0,y0},{x1,y0},{x1,y1},{x0,y1}};
        for(const auto plane:diamondPlanes){if(polygon.empty())break;std::vector<P> clipped;clipped.reserve(8);
            const auto value=[&](const P& p){return plane[0]*p[0]+plane[1]*p[1]-plane[2];};
            auto a=polygon.back();auto fa=value(a);
            for(const auto z:polygon){const auto fz=value(z);
                if((fa<=0)!=(fz<=0)){const auto t=fa/(fa-fz);clipped.push_back({a[0]+t*(z[0]-a[0]),a[1]+t*(z[1]-a[1])});}
                if(fz<=0)clipped.push_back(z);a=z;fa=fz;}polygon=std::move(clipped);}
        long double twiceArea=0;
        for(std::size_t k=0;k<polygon.size();++k){const auto a=polygon[k],z=polygon[(k+1)%polygon.size()];twiceArea+=a[0]*z[1]-a[1]*z[0];}
        diamondAlpha[ci]=static_cast<double>(std::abs(twiceArea)/(2*(x1-x0)*(y1-y0)));
        check(diamondAlpha[ci]>=0&&diamondAlpha[ci]<=1,"independent bounded diamond fraction outside [0,1]");}
    const auto diamondPublication=LAMP::Interface::GeometricVofPublication::prepareLocalPlic(*mesh,layout.domain(),clock,diamondAlpha);
    const auto diamondPhase=HydroPhaseGeometry::build(*mesh,metrics,diamondPublication);
    int diamondSource=-1;for(int p=0;p<layout.pressureCount();++p){const auto& b=diamondPublication.plicCells()[layout.nativeCell(p)].rectangle;
        const auto c=(b.lower+b.upper)*.5;if(c.x==diamondCenter.x+fineH&&c.y==diamondCenter.y+fineH)diamondSource=p;}
    check(diamondSource>=0,"same AMR scene lacks named straight diamond source");
    const auto diamondNative=layout.nativeCell(diamondSource);const auto& diamondCell=diamondPublication.plicCells()[diamondNative];
    check(diamondCell.resolved&&diamondCell.measured[0]>0&&diamondCell.measured[1]>0&&
        diamondCell.rectangle.upper.x-diamondCell.rectangle.lower.x==fineH,
        "diamond witness must be a real resolved fine source, ineligible for CU");
    unsigned faceMixedNeighbors=0;
    for(int u=0;u<layout.velocityCount();++u){const auto face=layout.domain().face(*mesh,layout.nativeFace(u));
        const int neighbor=face.negativeCell==diamondNative?face.positiveCell:face.positiveCell==diamondNative?face.negativeCell:-1;
        if(neighbor>=0){const auto& c=diamondPublication.plicCells()[neighbor];faceMixedNeighbors+=c.measured[0]>0&&c.measured[1]>0;}}
    check(faceMixedNeighbors==0,"45-degree diamond source must have only corner-contact mixed neighbors, so its old face cluster has one segment");
    LocalPlicEndpointReference diamondReference{diamondPhase};
    const auto diamondX=diamondReference.estimate(diamondSource,0),diamondY=diamondReference.estimate(diamondSource,1);
    check(!diamondX.valid&&!diamondY.valid&&diamondX.opposite&&diamondY.opposite,
        "closed narrow diamond must independently reject both original HF charts, which CT cannot repair");
    const auto diamondMid=LocalPlicEndpointReference::segmentMidpoint(diamondCell);
    const auto nNorm=std::hypot(static_cast<long double>(diamondCell.normal.x),static_cast<long double>(diamondCell.normal.y));
    const long double nx=diamondCell.normal.x/nNorm,ny=diamondCell.normal.y/nNorm;
    unsigned localStraightSegments=0;
    const Vector2D localLo{diamondMid.x-2.5*fineH,diamondMid.y-2.5*fineH},localHi{diamondMid.x+2.5*fineH,diamondMid.y+2.5*fineH};
    for(const auto& c:diamondPublication.plicCells()){if(!(c.resolved&&c.measured[0]>0&&c.measured[1]>0)||c.normal.x*nx+c.normal.y*ny<=0)continue;
        const Vector2D lo{std::max(localLo.x,c.rectangle.lower.x),std::max(localLo.y,c.rectangle.lower.y)},
            hi{std::min(localHi.x,c.rectangle.upper.x),std::min(localHi.y,c.rectangle.upper.y)};
        if(lo.x>=hi.x||lo.y>=hi.y)continue;
        const auto capacity=static_cast<long double>(hi.x-lo.x)*(hi.y-lo.y),liquid=LocalPlicEndpointReference::planeArea(c,lo,hi);
        if(!(liquid>0&&liquid<capacity))continue;
        const auto midpoint=LocalPlicEndpointReference::segmentMidpoint(c);
        const auto norm=std::hypot(static_cast<long double>(c.normal.x),static_cast<long double>(c.normal.y));
        close(nx*(static_cast<long double>(midpoint.x)-diamondMid.x)+ny*(static_cast<long double>(midpoint.y)-diamondMid.y),0,
            512.L*std::numeric_limits<double>::epsilon()*fineH,"fixed-box same-direction actual diamond segments share one straight plane");
        close(nx*c.normal.y/norm-ny*c.normal.x/norm,0,512.L*std::numeric_limits<double>::epsilon(),
            "straight diamond segment directions are independently parallel");++localStraightSegments;}
    check(localStraightSegments>=3,"same fixed 5h box must contain enough actual straight corner-contact segments");
    const auto diamondLoad=HydroBalancedCsfLoad::localPlicHeightFunctionCapillary(*mesh,pressure,diamondPhase,sigma,heightConfig());
    diamondLoad.validateVofFor(*mesh,pressure,diamondPublication.identity(),clock);
    close(diamondLoad.heightFunctionCellCurvatures()[diamondSource],0,
        32768.L*std::numeric_limits<double>::epsilon()/fineH,
        "bounded corner-contact quadratic fallback must recover analytic zero curvature on the straight diamond face");
    check(diamondLoad.curvatureReport().heightFunctionFitCells>0&&diamondLoad.curvatureReport().heightFunctionFitSegments>=localStraightSegments,
        "actual corner-contact fallback use must retain existing fit work counters");
    std::cout<<"LOCAL_PLIC_HEIGHT_CORNER_CONTACT,native,"<<diamondNative<<",face_mixed_neighbors,"<<faceMixedNeighbors
        <<",fixed_box_straight_segments,"<<localStraightSegments<<",kappa,"<<diamondLoad.heightFunctionCellCurvatures()[diamondSource]
        <<",analytic_kappa,0,fit_cells,"<<diamondLoad.curvatureReport().heightFunctionFitCells
        <<",fit_segments,"<<diamondLoad.curvatureReport().heightFunctionFitSegments<<",PASS\n";
    const auto wholeLocalSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-localBegin).count();
    std::cout<<"LOCAL_PLIC_HEIGHT_AMR,cells,"<<cells<<",loaded,"<<loaded<<",max_abs_curvature,"<<maximum
        <<",stock_only_micro_cell,"<<microCell<<",stock_only_micro_rows,"<<microRows
        <<",retained_micro_liquid,"<<micro.inventory[0]<<",actual_micro_liquid,"<<micro.measured[0]
        <<",whole_local_seconds,"<<wholeLocalSeconds<<",PASS (no dynamic accuracy claim)\n";
}

std::unique_ptr<QuadtreeMesh> localPlicGradedMesh()
{
    DictionaryMesh cfg;cfg.Nx=cfg.Ny=16;cfg.maxCellNum=65536;cfg.maxInitialIterations=0;cfg.maxLevel_local={3};
    auto mesh=std::make_unique<QuadtreeMesh>(cfg);
    // Use the production target-only Core adapt path. Static refineBoxes
    // alone does not create a balanced graded transition from level0 to3.
    const auto inside=[](Vector2D p,Vector2D lo,Vector2D hi){
        return p.x>=lo.x&&p.x<=hi.x&&p.y>=lo.y&&p.y<=hi.y;};
    for(int sweep=0;sweep<=5;++sweep){
        CellField<double> targets(*mesh,0.),gradient(*mesh,0.);
        for(int i=0;i<mesh->cellIndexNum();++i){const auto p=mesh->cells(mesh->cellIndexToLabel[i]).centerCoord();
            targets.value[i]=inside(p,{.25,.25},{.75,.8125})?3:
                inside(p,{.1875,.1875},{.8125,.875})?2:inside(p,{.125,.125},{.875,.9375})?1:0;}
        if(!mesh->mayAdaptToTargets(targets))return mesh;
        check(sweep<5,"graded HF production Core target adaptation did not settle");
        std::vector<VertField<double>> levelSets;std::vector<VertField<double>*> vertices;
        std::vector<HEdgeField<double>*> horizontal;std::vector<VEdgeField<double>*> vertical;
        std::vector<CellField<double>*> cells{&targets};
        mesh->adapt(levelSets,gradient,vertices,horizontal,vertical,cells,&targets);
    }
    throw std::logic_error("unreachable graded HF fixture state");
}
long double localPlicCapsuleArea(Vector2D lo,Vector2D hi)
{
    // Independent exact strip + two clipped analytic half disks. The two
    // straight sides are real closed-interface plane patches, not a wall.
    constexpr double a=.0625,b=.18;const Vector2D center{.5,.5};
    const double x0=std::max(lo.x,center.x-a),x1=std::min(hi.x,center.x+a);
    const double y0=std::max(lo.y,center.y-b),y1=std::min(hi.y,center.y+b);
    long double area=x1>x0&&y1>y0?static_cast<long double>(x1-x0)*(y1-y0):0;
    if(lo.x<center.x-a)area+=circleCell(lo,{std::min(hi.x,center.x-a),hi.y},{center.x-a,center.y},b,false);
    if(hi.x>center.x+a)area+=circleCell({std::max(lo.x,center.x+a),lo.y},hi,{center.x+a,center.y},b,false);
    return area;
}
void localPlicGradedHeightFunction()
{
    constexpr double sigma=.7,radius=.2;const Vector2D center{.487,.521};
    const LAMP::Interface::GeometricVofTimeLayer clock{7,.002,.0022,LAMP::TimeLayer::Position::Midpoint};
    for(bool circle:{false,true}){
        auto mesh=localPlicGradedMesh();const auto metrics=metricsFor(*mesh,false);const auto& layout=metrics.layout();
        std::set<int> levels;double smallest=1,largest=0;
        for(int ci=0;ci<mesh->cellIndexNum();++ci){const auto& c=mesh->cells(mesh->cellIndexToLabel[ci]);levels.insert(c.level());
            smallest=std::min(smallest,c.size().x);largest=std::max(largest,c.size().x);}
        check(levels==std::set<int>({0,1,2,3})&&largest==8*smallest,"graded HF requires actual levels0..3 and far coarse8h");
        std::size_t hanging=0;
        for(int u=0;u<layout.velocityCount();++u){const auto f=layout.domain().face(*mesh,layout.nativeFace(u));
            if(f.negativeCell<0||f.positiveCell<0)continue;
            const int a=mesh->cells(mesh->cellIndexToLabel[f.negativeCell]).level(),b=mesh->cells(mesh->cellIndexToLabel[f.positiveCell]).level();
            check(std::abs(a-b)<=1,"graded HF fixture violates actual shared-subface local2:1");hanging+=a!=b;}
        check(hanging>0,"graded HF fixture has no real hanging subfaces");
        std::vector<double> alpha(mesh->cellIndexNum());
        for(int ci=0;ci<mesh->cellIndexNum();++ci){const auto& c=mesh->cells(mesh->cellIndexToLabel[ci]);
            const auto lo=mesh->verts(c.cellToVert[0]).coord(),hi=mesh->verts(c.cellToVert[2]).coord();
            const auto liquid=circle?circleCell(lo,hi,center,radius,false):localPlicCapsuleArea(lo,hi);
            alpha[ci]=static_cast<double>(liquid/(static_cast<long double>(hi.x-lo.x)*(hi.y-lo.y)));
            check(alpha[ci]>=0&&alpha[ci]<=1,"independent graded HF input outside true capacity");}
        const auto publication=LAMP::Interface::GeometricVofPublication::prepareLocalPlic(*mesh,layout.domain(),clock,alpha);
        const auto phase=HydroPhaseGeometry::build(*mesh,metrics,publication);const auto pressure=HydroPressureLoad::build(metrics.divergence());
        std::vector<std::array<long double,2>> stocks;for(const auto& c:publication.plicCells())stocks.push_back(c.inventory);
        const auto id=publication.identity();const int leaves=mesh->cellIndexNum();
        const auto load=HydroBalancedCsfLoad::localPlicHeightFunctionCapillary(*mesh,pressure,phase,sigma,heightConfig());
        load.validateVofFor(*mesh,pressure,id,clock);double maximum=0,straightMaximum=0;std::size_t mixed=0,straight=0;
        for(int p=0;p<layout.pressureCount();++p){const int native=layout.nativeCell(p);const auto& c=publication.plicCells()[native];
            const bool actualMixed=c.measured[0]>0&&c.measured[1]>0;const double k=load.heightFunctionCellCurvatures()[p];
            if(!actualMixed){check(!std::isfinite(k),"graded pure cell fabricated HF curvature");continue;}
            ++mixed;check(mesh->cells(mesh->cellIndexToLabel[native]).level()==3&&std::isfinite(k),"graded actual interface must be finest with finite HF");
            if(circle){check(k>0,"graded circle lost analytic curvature sign");maximum=std::max(maximum,std::abs(k-1/radius)*radius);}
            else{const auto midpoint=(c.rectangle.lower+c.rectangle.upper)*.5;
                if(std::abs(midpoint.x-.5)<=1./128&&std::abs(midpoint.y-.5)>.15){
                    straightMaximum=std::max(straightMaximum,std::abs(k));
                    close(k,0,1e-10/smallest,"graded closed straight plane patch acquired curvature");++straight;}}
        }
        check(mixed>0&&load.curvatureReport().heightFunctionCells==mixed&&load.curvatureReport().amrFaces>0,"graded true HF work/source witness missing");
        if(circle)check(maximum<.025,"graded HF independent analytic circle curvature accuracy (old fine-scale gate)");
        else check(straight>=4,"graded closed straight plane witnesses missing");
        const Reference incidence(*mesh,phase);std::vector<double> color(layout.pressureCount());
        for(int p=0;p<layout.pressureCount();++p){const auto& c=publication.plicCells()[layout.nativeCell(p)];
            color[p]=c.measured[0]==0?0:c.measured[1]==0?1:static_cast<double>(c.measured[0]/metrics.cellVolumes()[p]);}
        const auto gradient=incidence.gradient(color);
        for(int u=0;u<layout.velocityCount();++u){close(load.colorGradient()[u],gradient[u],1e-12*std::max(1.,std::abs(gradient[u])),"graded actual C measured-color source");
            if(load.supportReports()[u].curvatureDefined)close(load.momentumRhsLoad()[u],sigma*load.curvatureAverages()[u]*gradient[u],
                1e-12*std::max(1.,std::abs(load.momentumRhsLoad()[u])),"graded unchanged sigma*kappa*Ccolor force");}
        auto wrong=clock;++wrong.step;rejects([&]{load.validateVofFor(*mesh,pressure,id,wrong);},"clock");
        const auto foreign=LAMP::Interface::GeometricVofPublication::prepareLocalPlic(*mesh,layout.domain(),clock,alpha);
        rejects([&]{load.validateVofFor(*mesh,pressure,foreign.identity(),clock);},"owner");
        if(circle){auto wallAlpha=alpha;
            for(int ci=0;ci<mesh->cellIndexNum();++ci){const auto& c=mesh->cells(mesh->cellIndexToLabel[ci]);
                const auto lo=mesh->verts(c.cellToVert[0]).coord(),hi=mesh->verts(c.cellToVert[2]).coord();
                wallAlpha[ci]=static_cast<double>(circleCell(lo,hi,{.05,.5},radius,false)/(static_cast<long double>(hi.x-lo.x)*(hi.y-lo.y)));}
            const auto wall=LAMP::Interface::GeometricVofPublication::prepareLocalPlic(*mesh,layout.domain(),clock,wallAlpha);
            const auto wallPhase=HydroPhaseGeometry::build(*mesh,metrics,wall);
            rejects([&]{HydroBalancedCsfLoad::localPlicHeightFunctionCapillary(*mesh,pressure,wallPhase,sigma,heightConfig());},"physical-contact/boundary-cell support unsupported");}
        check(publication.identity()==id&&publication.alpha()==alpha&&mesh->cellIndexNum()==leaves,"graded HF/failed checks changed source owner or mesh");
        for(std::size_t i=0;i<stocks.size();++i)check(publication.plicCells()[i].inventory==stocks[i],"graded HF/failed checks changed independent source stocks");
        std::cout<<"LOCAL_PLIC_HEIGHT_GRADED,shape,"<<(circle?"circle":"closed_plane_patch")<<",cells,"<<leaves<<",levels,4,far_h_over_finest,8"
            <<",mixed,"<<mixed<<",straight_cells,"<<straight<<",straight_max,"<<straightMaximum<<",circle_relative_max,"<<maximum
            <<",queries,"<<load.curvatureReport().heightFunctionCellQueries<<",PASS (no dynamic accuracy claim)\n";
    }
}



// RZ-C1 independent physical references. No production cut, reconstruction,
// meridional-box port or curvature helper is used by these phase integrals.
constexpr long double rzPhysicalPi=static_cast<long double>(CartesianMeasure::twoPi)/2;
using RzStocks=std::array<long double,2>;
RzStocks rzSphereColumn(Vector2D lo,Vector2D hi,long double center,long double radius)
{
    const long double r0=lo.x,r1=hi.x,a=lo.y,b=hi.y,rr=radius*radius;
    check(r0>=0&&r1>r0&&b>a&&radius>0,"invalid independent sphere column");
    std::vector<long double> cuts{a,b};
    for(auto r:{r0,r1})if(rr>=r*r){const auto z=std::sqrt(rr-r*r);
        for(auto point:{center-z,center+z})if(point>a&&point<b)cuts.push_back(point);}
    std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
    RzStocks result{};
    for(std::size_t i=1;i<cuts.size();++i){const auto lower=cuts[i-1],upper=cuts[i],length=upper-lower;
        const auto t=(lower+upper)/2-center,s=rr-t*t;
        long double liquid=0,gas=0;
        if(s<=r0*r0)gas=(r1-r0)*(r1+r0);
        else if(s>=r1*r1)liquid=(r1-r0)*(r1+r0);
        else{
            // Integrate the quadratic cross section around this interval's
            // midpoint: avg((z-center)^2)=t^2+length^2/12.
            const auto meanSquare=t*t+length*length/12;
            liquid=rr-r0*r0-meanSquare;gas=r1*r1-rr+meanSquare;
        }
        check(liquid>=0&&gas>=0,"independent sphere phase integral became negative");
        result[0]+=rzPhysicalPi*length*liquid;result[1]+=rzPhysicalPi*length*gas;
    }
    const auto capacity=rzPhysicalPi*(r1-r0)*(r1+r0)*(b-a);
    check(std::abs(result[0]+result[1]-capacity)<=1e-12L*capacity,"independent sphere two-phase physical closure");
    return result;
}
RzStocks rzCapsuleColumn(Vector2D lo,Vector2D hi,long double lowerCenter,long double upperCenter,long double radius)
{
    // Exact geometric classification against the axial centre segment:
    // a pure whole rectangle needs one whole-box integral, not a sum of
    // separately rounded sphere/straight pieces. No stock is clipped.
    check(lo.x>=0&&hi.x>lo.x&&hi.y>lo.y&&lowerCenter<=upperCenter&&radius>0,
        "invalid independent capsule rectangle");
    const auto distance=[&](long double z){return z<lowerCenter?lowerCenter-z:z>upperCenter?z-upperCenter:0.L;};
    const long double minimumZ=hi.y<lowerCenter?lowerCenter-hi.y:lo.y>upperCenter?lo.y-upperCenter:0.L;
    const long double maximumZ=std::max(distance(lo.y),distance(hi.y)),rr=radius*radius;
    const long double r0=lo.x,r1=hi.x;
    const auto wholeCapacity=rzPhysicalPi*(r1-r0)*(r1+r0)*(static_cast<long double>(hi.y)-lo.y);
    if(r0*r0+minimumZ*minimumZ>=rr)return {0,wholeCapacity};
    if(r1*r1+maximumZ*maximumZ<=rr)return {wholeCapacity,0};
    RzStocks result{};
    const auto add=[&](RzStocks x){result[0]+=x[0];result[1]+=x[1];};
    if(lo.y<lowerCenter){auto upper=hi;upper.y=static_cast<double>(std::min(static_cast<long double>(hi.y),lowerCenter));
        add(rzSphereColumn(lo,upper,lowerCenter,radius));}
    const auto a=std::max(static_cast<long double>(lo.y),lowerCenter),b=std::min(static_cast<long double>(hi.y),upperCenter);
    if(b>a){const long double r0=lo.x,r1=hi.x;long double liquid=0,gas=0;
        if(r1<=radius)liquid=(r1-r0)*(r1+r0);
        else if(r0>=radius)gas=(r1-r0)*(r1+r0);
        else{liquid=(radius-r0)*(radius+r0);gas=(r1-radius)*(r1+radius);}
        add({rzPhysicalPi*(b-a)*liquid,rzPhysicalPi*(b-a)*gas});}
    if(hi.y>upperCenter){auto lower=lo;lower.y=static_cast<double>(std::max(static_cast<long double>(lo.y),upperCenter));
        add(rzSphereColumn(lower,hi,upperCenter,radius));}
    const auto capacity=rzPhysicalPi*(static_cast<long double>(hi.x)-lo.x)*(static_cast<long double>(hi.x)+lo.x)*(static_cast<long double>(hi.y)-lo.y);
    check(std::abs(result[0]+result[1]-capacity)<=1e-12L*capacity,"independent capsule two-phase physical closure");return result;
}
template<class F>std::vector<RzStocks> rzAnalyticStocks(const QuadtreeMesh& mesh,F integral)
{
    std::vector<RzStocks> stocks;stocks.reserve(mesh.cellIndexNum());
    for(int ci=0;ci<mesh.cellIndexNum();++ci){const auto& cell=mesh.cells(mesh.cellIndexToLabel[ci]);
        stocks.push_back(integral(mesh.verts(cell.cellToVert[0]).coord(),mesh.verts(cell.cellToVert[2]).coord()));}
    return stocks;
}
void rzAnalyticHeightJets()
{
    constexpr long double radius=.23L,center=.513L;
    const auto verify=[&](HydroCsfHeightJet jet,long double expected){
        const auto actual=HydroBalancedCsfLoad::axisymmetricHeightCurvatures(jet);
        check(std::abs(actual[0]-expected)<=1e-12L/radius&&std::abs(actual[1]-expected)<=1e-12L/radius,
            "independent sphere analytic jet meridional/hoop sign or magnitude");};
    for(int flip:{-1,1}){
        // r=g(z), with H measured from the common liquid-side origin toward s.
        const long double z=.3L*radius,g=std::sqrt(radius*radius-z*z);
        HydroCsfHeightJet radial;radial.heightAxis=0;radial.outward=flip;radial.origin=flip>0?0:2*radius;
        radial.height=flip*(g-radial.origin);radial.first=flip*(-z/g);
        radial.second=flip*(-radius*radius/(g*g*g));radial.tangentCoordinate=center+z;verify(radial,flip/radius);
        for(int cap:{-1,1})for(long double r:{0.L,.3L*radius}){
            const auto h=std::sqrt(radius*radius-r*r),surface=center+cap*h;
            HydroCsfHeightJet axial;axial.heightAxis=1;axial.outward=flip*cap;
            axial.origin=axial.outward>0?center-2*radius:center+2*radius;axial.height=axial.outward*(surface-axial.origin);
            axial.first=axial.outward*(-cap*r/h);axial.second=axial.outward*(-cap*radius*radius/(h*h*h));
            axial.tangentCoordinate=r;axial.axisRegular=r==0;verify(axial,flip/radius);
        }
        HydroCsfHeightJet cylinder;cylinder.heightAxis=0;cylinder.outward=flip;cylinder.origin=flip>0?0:2*radius;
        cylinder.height=radius;const auto k=HydroBalancedCsfLoad::axisymmetricHeightCurvatures(cylinder);
        check(k[0]==0&&std::abs(k[1]-flip/radius)<=1e-12L/radius,"independent cylinder jet has wrong 0+1/R or phase sign");
    }
    HydroCsfHeightJet bad;bad.heightAxis=1;bad.height=1;bad.tangentCoordinate=0;bad.second=-1/radius;
    rejects([&]{HydroBalancedCsfLoad::axisymmetricHeightCurvatures(bad);},"even regular axial height jet");
    bad.axisRegular=true;bad.first=.1L;
    rejects([&]{HydroBalancedCsfLoad::axisymmetricHeightCurvatures(bad);},"even regular axial height jet");
    bad={};bad.height=.1L;bad.origin=-.2L;
    rejects([&]{HydroBalancedCsfLoad::axisymmetricHeightCurvatures(bad);},"invalid interface radius");
    std::cout<<"RZ_C1_ANALYTIC_JETS,SPHERE_2_OVER_R,CYLINDER_1_OVER_R,PHASE_FLIPS,AXIS_LIMIT,PASS\n";
}
struct RzHeightError{double maximum=0,rms=0,hmin=std::numeric_limits<double>::infinity(),hmax=0;
    long double effectiveHmin=std::numeric_limits<long double>::infinity(),effectiveHmax=0;
    std::size_t mixed=0,axisMixed=0,reflected=0,jets=0,fit=0,interfaceHanging=0,coarseAligned=0;};
void rzActualForceAlgebra(const QuadtreeMesh& mesh,const HydroPhaseGeometry& phase,const HydroPressureLoad& pressure,
    const HydroBalancedCsfLoad& load,double sigma,double constantKappa,bool checkCylinder=false,double lower=.32,double upper=.68)
{
    const auto& layout=phase.layout();Reference ref(mesh,phase);std::vector<double> color(layout.pressureCount());
    for(int p=0;p<layout.pressureCount();++p){const auto& cell=phase.vofPublication().plicCells()[layout.nativeCell(p)];
        color[p]=cell.measured[0]==0?0:cell.measured[1]==0?1:static_cast<double>(cell.measured[0]/ref.volume[p]);}
    const auto gradient=ref.gradient(color);auto potential=color;for(auto& p:potential)p*=sigma*constantKappa;
    const auto cp=pressure.apply(mesh,HydroNativeValues(layout,std::vector<double>(layout.velocityCount()),potential));
    double gradientError=0,gradientScale=0,forceError=0,forceScale=0,constantError=0,constantScale=0,straightError=0,straightScale=0;
    std::size_t straightActive=0;
    for(int u=0;u<layout.velocityCount();++u){
        gradientError=std::hypot(gradientError,load.colorGradient()[u]-gradient[u]);gradientScale=std::hypot(gradientScale,gradient[u]);
        const auto expected=load.supportReports()[u].curvatureDefined?sigma*load.curvatureAverages()[u]*gradient[u]:0;
        forceError=std::hypot(forceError,load.momentumRhsLoad()[u]-expected);forceScale=std::hypot(forceScale,expected);
        constantError=std::hypot(constantError,cp[u]-sigma*constantKappa*gradient[u]);constantScale=std::hypot(constantScale,cp[u]);
        if(ref.area[u]==0)check(cp[u]==0&&load.colorGradient()[u]==0&&load.momentumRhsLoad()[u]==0,"actual RZ axis has nonzero C/capillary load");
        if(checkCylinder){const auto f=layout.domain().face(mesh,layout.nativeFace(u));
            if(f.center.y>lower+.1&&f.center.y<upper-.1){
                straightError=std::hypot(straightError,load.momentumRhsLoad()[u]-cp[u]);straightScale=std::hypot(straightScale,cp[u]);
                straightActive+=cp[u]!=0;}}
    }
    check(relative(gradientError,gradientScale)<=1e-12&&relative(forceError,forceScale)<=1e-12,"actual RZ same-plane color/force C incidence");
    check(relative(constantError,constantScale)<=1e-12,"constant-curvature physical RZ C balance");
    if(checkCylinder)check(straightActive>=4&&relative(straightError,straightScale)<=1e-12,"actual cylinder middle HF load fails constant 1/R C balance");
}


RzHeightError rzSphereHeightScene(int n,bool amr)
{
    constexpr double radius=.23,center=.513,sigma=.7;
    auto mesh=heightMesh(n);if(amr){DictionaryMesh::RefineBox box;box.level=1;box.min={0,.5};box.max={.4,.85};mesh->applyRefineBoxes({box});}
    const auto metrics=metricsFor(*mesh,true);const auto& layout=metrics.layout();
    const auto stocks=rzAnalyticStocks(*mesh,[&](Vector2D a,Vector2D b){return rzSphereColumn(a,b,center,radius);});
    long double liquid=0;for(const auto& x:stocks)liquid+=x[0];
    const long double volume=4*rzPhysicalPi*radius*radius*radius/3;
    check(std::abs(liquid-volume)<=1e-12L*volume,"independent sphere physical volume is not 4*rzPhysicalPi*R^3/3");
    const LAMP::TimeLayer clock{11,.01,.0102,LAMP::TimeLayer::Position::Midpoint};
    const auto publication=LAMP::Interface::GeometricVofPublication::prepareLocalPlicInventory(*mesh,layout.domain(),clock,stocks,{});
    const auto alpha=publication.alpha();const auto id=publication.identity();const int leaves=mesh->cellIndexNum();
    const auto phase=HydroPhaseGeometry::build(*mesh,metrics,publication);const auto pressure=HydroPressureLoad::build(metrics.divergence());
    const auto load=HydroBalancedCsfLoad::localPlicHeightFunctionCapillary(*mesh,pressure,phase,sigma,heightConfig());
    load.validateVofFor(*mesh,pressure,id,clock);
    check(load.vofIdentity()==id&&&load.phaseGeometry().cellVolumes()==&phase.cellVolumes(),"RZ sphere load replaced actual phase owner");
    const auto& jets=load.heightFunctionCellJets();check(jets.size()==static_cast<std::size_t>(layout.pressureCount()),"RZ sphere jet layout mismatch");
    RzHeightError error;long double squared=0;std::set<int> interfaceLevels;
    for(int p=0;p<layout.pressureCount();++p){const int native=layout.nativeCell(p);const auto& cell=publication.plicCells()[native];
        const auto k=load.heightFunctionCellCurvatures()[p];const bool mixed=cell.resolved&&cell.measured[0]>0&&cell.measured[1]>0;
        if(!mixed){check(!std::isfinite(k)&&!jets[p],"pure/unresolved RZ source fabricated a height jet/curvature");continue;}
        ++error.mixed;interfaceLevels.insert(mesh->cells(mesh->cellIndexToLabel[native]).level());
        const double h=std::max(cell.rectangle.upper.x-cell.rectangle.lower.x,cell.rectangle.upper.y-cell.rectangle.lower.y);
        error.hmin=std::min(error.hmin,h);error.hmax=std::max(error.hmax,h);
        check(std::isfinite(k)&&k>0,"actual sphere HF lacks finite positive liquid-to-gas curvature");
        const auto relativeError=std::abs(k-2/radius)/(2/radius);
        if(relativeError>error.maximum){std::cout<<"RZ_C1_WORST_DIAGNOSTIC,n="<<n<<",amr="<<amr<<",native="<<native<<",relative="<<relativeError<<",k="<<k<<",r0="<<cell.rectangle.lower.x<<",r1="<<cell.rectangle.upper.x<<",z0="<<cell.rectangle.lower.y<<",z1="<<cell.rectangle.upper.y<<",nx="<<cell.normal.x<<",ny="<<cell.normal.y;
            if(jets[p]){const auto& j=*jets[p];const auto parts=HydroBalancedCsfLoad::axisymmetricHeightCurvatures(j);std::cout<<",axis="<<j.heightAxis<<",outward="<<j.outward<<",height="<<j.height<<",first="<<j.first<<",second="<<j.second<<",origin="<<j.origin<<",tangent="<<j.tangentCoordinate<<",meridional="<<parts[0]<<",hoop="<<parts[1]<<",reflected="<<j.axisReflected<<",normal_spacing="<<j.normalSpacing<<",tangent_spacing="<<j.tangentSpacing<<",coarse_aligned="<<j.coarseAligned;}else std::cout<<",fit=1";std::cout<<std::endl;}
        error.maximum=std::max(error.maximum,relativeError);squared+=relativeError*relativeError;
        if(jets[p]){++error.jets;
            const auto& jet=*jets[p];
            check(std::isfinite(jet.normalSpacing)&&jet.normalSpacing>0&&std::isfinite(jet.tangentSpacing)&&jet.tangentSpacing>0,
                "actual RZ qualified jet lacks its effective physical chart spacings");
            const auto effective=std::max(jet.normalSpacing,jet.tangentSpacing);
            error.effectiveHmin=std::min(error.effectiveHmin,effective);error.effectiveHmax=std::max(error.effectiveHmax,effective);
            error.coarseAligned+=jet.coarseAligned;
            const auto parts=HydroBalancedCsfLoad::axisymmetricHeightCurvatures(*jets[p]);
            check(std::abs(parts[0]+parts[1]-k)<=1e-12L*(2/radius),"actual RZ curvature does not consume its retained same-chart height jet");
            error.reflected+=jets[p]->axisReflected;
        }else ++error.fit;
        if(cell.rectangle.lower.x==0){++error.axisMixed;
            check(jets[p]&&jets[p]->heightAxis==1&&jets[p]->axisReflected,"true sphere axis cap was not qualified by reflected axial HF");}
    }
    for(int u=0;u<layout.velocityCount();++u){const auto f=layout.domain().face(*mesh,layout.nativeFace(u));
        if(f.negativeCell<0||f.positiveCell<0)continue;
        const auto& a=mesh->cells(mesh->cellIndexToLabel[f.negativeCell]);const auto& b=mesh->cells(mesh->cellIndexToLabel[f.positiveCell]);
        check(std::abs(a.level()-b.level())<=1,"RZ sphere fixture violates actual shared-face2:1");
        const auto& ca=publication.plicCells()[f.negativeCell];const auto& cb=publication.plicCells()[f.positiveCell];
        error.interfaceHanging+=a.level()!=b.level()&&ca.measured[0]>0&&ca.measured[1]>0&&cb.measured[0]>0&&cb.measured[1]>0;
    }
    check(error.mixed>0&&error.axisMixed>=2&&error.reflected>=2&&error.jets>0,"actual sphere/axis-reflection witnesses missing");
    if(amr)check(interfaceLevels==std::set<int>({0,1})&&error.interfaceHanging>0,"RZ sphere did not cross a real interface hanging face");
    check(load.curvatureReport().heightFunctionCells+load.curvatureReport().heightFunctionFitCells==error.mixed,"RZ sphere height/fit report lost actual mixed source count");
    error.rms=std::sqrt(static_cast<double>(squared/error.mixed));
    check(load.curvatureReport().heightFunctionCoarseAlignedCells==error.coarseAligned,
        "actual final coarse-aligned jet count differs from retained report");
    check(load.curvatureReport().heightFunctionCoarseAlignedCharts>=error.coarseAligned,
        "actual final coarse-aligned cells exceed prequalified directional charts");
    // Diagnostics precede accuracy rejection; no PASS token is printed here.
    // Native-source h and effective chart h are distinct measured quantities.
    std::cout<<"RZ_C1_SPHERE_DIAGNOSTIC,"<<n<<",amr,"<<amr<<",leaves,"<<leaves<<",mixed,"<<error.mixed<<",axis_mixed,"<<error.axisMixed
        <<",reflected,"<<error.reflected<<",jet,"<<error.jets<<",fit,"<<error.fit<<",interface_hanging,"<<error.interfaceHanging
        <<",relative_max,"<<error.maximum<<",relative_rms,"<<error.rms<<",source_h_over_R_min,"<<error.hmin/radius
        <<",source_h_over_R_max,"<<error.hmax/radius<<",source_R_over_h_min,"<<radius/error.hmax
        <<",effective_h_over_R_min,"<<error.effectiveHmin/radius<<",effective_h_over_R_max,"<<error.effectiveHmax/radius
        <<",coarse_aligned_charts,"<<load.curvatureReport().heightFunctionCoarseAlignedCharts<<",coarse_aligned_cells,"<<error.coarseAligned<<std::endl;
    if(amr)check(load.curvatureReport().heightFunctionCoarseAlignedCharts>0&&error.coarseAligned>0,
        "actual AMR sphere never consumed its prequalified coarse-aligned height chart");
    // New RZ finite-grid accuracy gates, frozen before the first direct run.
    // Old XY exact/algebra/analytic thresholds elsewhere are unchanged.
    check(error.maximum<(amr?.06:.15),"RZ sphere independent analytic total-curvature finite-grid gate");
    rzActualForceAlgebra(*mesh,phase,pressure,load,sigma,2/radius);
    auto wrongClock=clock;++wrongClock.step;rejects([&]{load.validateVofFor(*mesh,pressure,id,wrongClock);},"clock");
    const auto zero=HydroBalancedCsfLoad::localPlicHeightFunctionCapillary(*mesh,pressure,phase,0,heightConfig());
    check(zero.curvatureReport().projectionSkippedSigmaZero,"RZ sigma0 did not skip HF queries");
    for(double x:zero.momentumRhsLoad())check(x==0,"RZ sigma0 force is nonzero");
    for(const auto& jet:zero.heightFunctionCellJets())check(!jet,"RZ sigma0 fabricated retained height jets");
    check(publication.identity()==id&&publication.alpha()==alpha&&mesh->cellIndexNum()==leaves,"RZ sphere read/failure changed source owner/mesh");
    for(std::size_t ci=0;ci<stocks.size();++ci)check(publication.plicCells()[ci].inventory==stocks[ci],"RZ sphere HF changed independent input stocks");
    std::cout<<"RZ_C1_SPHERE,"<<n<<",amr,"<<amr<<",leaves,"<<leaves<<",mixed,"<<error.mixed<<",axis_mixed,"<<error.axisMixed
        <<",reflected,"<<error.reflected<<",jet,"<<error.jets<<",fit,"<<error.fit<<",interface_hanging,"<<error.interfaceHanging
        <<",relative_max,"<<error.maximum<<",relative_rms,"<<error.rms<<",h_over_R_min,"<<error.hmin/radius
        <<",h_over_R_max,"<<error.hmax/radius<<",R_over_h_min,"<<radius/error.hmax
        <<",effective_h_over_R_min,"<<error.effectiveHmin/radius<<",effective_h_over_R_max,"<<error.effectiveHmax/radius
        <<",coarse_aligned_charts,"<<load.curvatureReport().heightFunctionCoarseAlignedCharts<<",coarse_aligned_cells,"<<error.coarseAligned<<",PASS\n";return error;
}
void rzCylinderHeightScene()
{
    constexpr int n=64;constexpr double radius=.14,lower=.32,upper=.68,sigma=.7;
    auto mesh=heightMesh(n);const auto metrics=metricsFor(*mesh,true);const auto& layout=metrics.layout();
    const auto stocks=rzAnalyticStocks(*mesh,[&](Vector2D a,Vector2D b){return rzCapsuleColumn(a,b,lower,upper,radius);});
    long double liquid=0;for(const auto& x:stocks)liquid+=x[0];
    const long double volume=rzPhysicalPi*radius*radius*(upper-lower)+4*rzPhysicalPi*radius*radius*radius/3;
    check(std::abs(liquid-volume)<=1e-12L*volume,"independent closed capsule volume oracle");
    const LAMP::TimeLayer clock{12,.01,.0102,LAMP::TimeLayer::Position::Midpoint};
    const auto publication=LAMP::Interface::GeometricVofPublication::prepareLocalPlicInventory(*mesh,layout.domain(),clock,stocks,{});
    const auto phase=HydroPhaseGeometry::build(*mesh,metrics,publication);const auto pressure=HydroPressureLoad::build(metrics.divergence());
    const auto load=HydroBalancedCsfLoad::localPlicHeightFunctionCapillary(*mesh,pressure,phase,sigma,heightConfig());
    load.validateVofFor(*mesh,pressure,publication.identity(),clock);std::size_t straight=0;
    for(int p=0;p<layout.pressureCount();++p){const auto& cell=publication.plicCells()[layout.nativeCell(p)];
        if(!(cell.resolved&&cell.measured[0]>0&&cell.measured[1]>0))continue;
        if(cell.rectangle.lower.y>lower+.1&&cell.rectangle.upper.y<upper-.1){const auto k=load.heightFunctionCellCurvatures()[p];
            check(std::isfinite(k)&&std::abs(k-1/radius)<=1e-12/radius,"closed capsule straight patch is not independent cylinder1/R");
            const auto& jet=load.heightFunctionCellJets()[p];check(jet&&jet->heightAxis==0,"straight cylinder patch lacks radial height chart");
            const auto parts=HydroBalancedCsfLoad::axisymmetricHeightCurvatures(*jet);
            check(std::abs(parts[0])<=1e-12L/radius&&std::abs(parts[1]-1/radius)<=1e-12L/radius,"straight cylinder has spurious meridional curvature or wrong hoop");++straight;}}
    check(straight>=4,"closed capsule has no resolved cylinder middle patch");
    rzActualForceAlgebra(*mesh,phase,pressure,load,sigma,1/radius,true,lower,upper);
    for(std::size_t ci=0;ci<stocks.size();++ci)check(publication.plicCells()[ci].inventory==stocks[ci],"cylinder HF changed independent input stocks");
    std::cout<<"RZ_C1_CYLINDER,"<<n<<",straight_mixed,"<<straight<<",kappa,"<<1/radius<<",C_balance,1e-12,PASS\n";
}
void rzPhysicalWallNegatives()
{
    const LAMP::TimeLayer clock{13,.01,.0102,LAMP::TimeLayer::Position::Midpoint};
    for(bool annulus:{false,true}){
        auto mesh=heightMesh(32,annulus?Vector2D{.125,0}:Vector2D{0,0});
        const auto metrics=metricsFor(*mesh,true);const auto& layout=metrics.layout();
        const double center=annulus?.513:.04,radius=.23;
        const auto stocks=rzAnalyticStocks(*mesh,[&](Vector2D a,Vector2D b){return rzSphereColumn(a,b,center,radius);});
        const auto publication=LAMP::Interface::GeometricVofPublication::prepareLocalPlicInventory(*mesh,layout.domain(),clock,stocks,{});
        const auto phase=HydroPhaseGeometry::build(*mesh,metrics,publication);const auto pressure=HydroPressureLoad::build(metrics.divergence());
        std::size_t contact=0;for(const auto& c:publication.plicCells())if(c.measured[0]>0&&c.measured[1]>0)
            contact+=annulus?c.rectangle.lower.x==mesh->minCoord.x:c.rectangle.lower.y==mesh->minCoord.y;
        check(contact>0,"RZ physical-wall negative lacks actual mixed wall cells");
        rejects([&]{HydroBalancedCsfLoad::localPlicHeightFunctionCapillary(*mesh,pressure,phase,.7,heightConfig());},"physical-contact/boundary-cell support unsupported");
        for(std::size_t ci=0;ci<stocks.size();++ci)check(publication.plicCells()[ci].inventory==stocks[ci],"failed physical-wall request changed RZ independent stocks");
        std::cout<<"RZ_C1_PHYSICAL_WALL_NEGATIVE,"<<(annulus?"nonzero_inner_radius":"axial_wall")<<",mixed_contact,"<<contact<<",PASS\n";
    }
}
void localPlicRzHeightFunction()
{
    rzAnalyticHeightJets();const auto coarse=rzSphereHeightScene(32,false),fine=rzSphereHeightScene(64,true);
    check(fine.rms<.8*coarse.rms,"independent RZ sphere RMS curvature error does not improve on refinement");
    rzCylinderHeightScene();rzPhysicalWallNegatives();
    std::cout<<"RZ_C1_LOCAL_PLIC_HEIGHT_RESULT,PASS (module capability only; no Step/App/full-window claim)\n";
}

void heightFunctionOnly()
{
    heightSmallContract();heightPlanePatches();localPlicHeightFunction();
    const auto coarse=heightCircle(32,{0,0},1,true),fine=heightCircle(64,{0,0},1,false);
    const auto shiftedCoarse=heightCircle(32,{-.25,1.5},-1,false),shiftedFine=heightCircle(64,{-.25,1.5},-1,false);
    for(const auto pair:{std::array<HeightCircleError,2>{coarse,fine},std::array<HeightCircleError,2>{shiftedCoarse,shiftedFine}}){
        check(pair[1].rms<.65*pair[0].rms&&pair[1].maximum<.65*pair[0].maximum,"HF independent analytic circle errors do not decrease on refinement");
        std::cout<<"HEIGHT_FUNCTION_REFINEMENT,rms_order,"<<std::log2(pair[0].rms/pair[1].rms)
            <<",maximum_order,"<<std::log2(pair[0].maximum/pair[1].maximum)<<",PASS\n";
    }
    heightActualHighFrequency();heightUnsupportedScopes();
    std::cout<<"HEIGHT_FUNCTION_ONLY,PASS (force estimator; no dynamic stability claim)\n";
}

void initialDropAcceleration()
{
    constexpr int n=64;constexpr double sigma=341.642,rhoL=100,rhoG=1;
    for(bool circle:{true,false}){
        auto mesh=makeMesh(n);const double a=circle?.16:.168,b=circle?.16:.15238095238095239;
        const Profile f=[=](Vector2D p){const double x=p.x-.5,y=p.y-.5;
            return Jet{x*x/(a*a)+y*y/(b*b)-1,2*x/(a*a),2*y/(b*b),0};};
        const auto state=stateFor(*mesh,f,false,{},[](int,Vector2D,Vector2D){return 0.;},{},HermiteGalsTransport::FirstJet);
        HermiteGalsCandidate candidate(*mesh,{},state.primary());
        const auto mixed=closeHermiteGalsMixedAmr(*mesh,candidate);const auto id=candidate.identity();
        const auto publication=publishHermiteGalsFirstJetGeometry(*mesh,std::move(candidate),mixed,{id,centers(*mesh,f)});
        const auto phase=HydroPhaseGeometry::build(*mesh,metricsFor(*mesh,false),publication);
        const auto& l=phase.layout();const auto pressure=HydroPressureLoad::build(phase.metrics().divergence());
        const auto dual=HydroMomentumGeometry::build(*mesh,phase.metrics());
        const auto mass=HydroMomentumMass::build(*mesh,dual,phase,rhoL,rhoG);
        HydroRegularizedColorConfig colorConfig;colorConfig.epsilon=1.5/n;colorConfig.relativeTolerance=1e-7;
        const auto color=HydroRegularizedColor::build(*mesh,phase,colorConfig);
        HydroCsfCurvatureConfig curvatureConfig;curvatureConfig.relativeQuadratureTolerance=1e-7;
        const auto load=HydroBalancedCsfLoad::capillary(*mesh,pressure,color,sigma,curvatureConfig);
        load.validateFor(*mesh,pressure,color);
        const std::vector<double> zero(l.velocityCount());
        const auto sharpGradient=pressure.apply(*mesh,HydroNativeValues(l,zero,phase.liquidFractions()));
        std::vector<double> sharpRhs(l.velocityCount());
        for(int u=0;u<l.velocityCount();++u){
            if(sharpGradient[u]!=0){
                check(load.supportReports()[u].curvatureDefined,"initial probe comparator outside curvature support");
                sharpRhs[u]=sigma*load.curvatureAverages()[u]*sharpGradient[u];
            }
        }
        Matrix block(l.velocityCount(),l.velocityCount());
        for(int u=0;u<l.velocityCount();++u)block.coeffRef(u,u)=mass.masses()[u];
        block.makeCompressed();
        const auto saddle=HydroSinglePhaseSaddle::buildClosed(pressure,HydroSaddleVelocityBlock::build(l,block));
        const auto current=saddle.solve(*mesh,HydroSaddleRhs(l,load.momentumRhsLoad(),zero));
        const auto sharp=saddle.solve(*mesh,HydroSaddleRhs(l,sharpRhs,zero));
        // Popinet (2009) p.5841: same-level c_f=(c_L+c_R)/2, then
        // rho_f=rho(c_f). Frozen diagnostic only; the production mass owner
        // remains the true sharp dual mass. No equal-density replacement.
        std::vector<double> primalMass(l.velocityCount());
        Matrix primalBlock(l.velocityCount(),l.velocityCount());
        for(int u=0;u<l.velocityCount();++u){
            const auto face=l.domain().face(*mesh,l.nativeFace(u));
            double cf=0;int count=0;
            for(int cell:{face.negativeCell,face.positiveCell})if(cell>=0){
                cf+=phase.liquidFractions()[l.pressureDof(cell)];++count;
            }
            check(count>0,"initial acceleration probe has no adjacent pressure cell");
            cf/=count;primalMass[u]=(rhoG+(rhoL-rhoG)*cf)*double(dual.volumes()[u].volume);
            primalBlock.coeffRef(u,u)=primalMass[u];
        }
        primalBlock.makeCompressed();
        const auto primalSaddle=HydroSinglePhaseSaddle::buildClosed(pressure,HydroSaddleVelocityBlock::build(l,primalBlock));
        const auto currentPrimal=primalSaddle.solve(*mesh,HydroSaddleRhs(l,load.momentumRhsLoad(),zero));
        const auto sharpPrimal=primalSaddle.solve(*mesh,HydroSaddleRhs(l,sharpRhs,zero));
        std::vector<bool> free(l.velocityCount(),false);for(int u:saddle.freeVelocityDofs())free[u]=true;
        check(current.report().algebraicRelativeResidual<=1e-11&&sharp.report().algebraicRelativeResidual<=1e-11,
            "initial acceleration saddle did not retain original algebraic gate");
        double currentPeak=0,sharpPeak=0,currentGasPeak=0,sharpGasPeak=0;
        double currentPrimalPeak=0,sharpPrimalPeak=0,currentPrimalGasPeak=0,sharpPrimalGasPeak=0;int gasCurrent=0,gasSharp=0;
        std::cout<<"CSF_INITIAL_FACE_HEADER,shape,dof,x,y,axis,free,volume,mass,dual_liquid_volume,"
            "gradient_epsilon,gradient_sharp,kappa_epsilon,rhs_epsilon,rhs_sharp,acc_epsilon,acc_sharp,mass_primal_average,acc_epsilon_primal_mass,acc_sharp_primal_mass\n";
        for(int u=0;u<l.velocityCount();++u){
            const auto& box=dual.volumes()[u];const auto face=l.domain().face(*mesh,l.nativeFace(u));
            const auto& edge=face.horizontal?mesh->hEdges(face.edgeLabel):mesh->vEdges(face.edgeLabel);
            const auto point=(mesh->verts(edge.edgeToVert[0]).coord()+mesh->verts(edge.edgeToVert[1]).coord())*.5;
            const int axis=l.normalAxis(u)==HydroDofLayout::NormalAxis::First?0:1;
            const double ca=current.values().normalVelocity()[u],sa=sharp.values().normalVelocity()[u];
            const double cpa=currentPrimal.values().normalVelocity()[u],spa=sharpPrimal.values().normalVelocity()[u];
            if(free[u]){currentPeak=std::max(currentPeak,std::abs(ca));sharpPeak=std::max(sharpPeak,std::abs(sa));
                currentPrimalPeak=std::max(currentPrimalPeak,std::abs(cpa));sharpPrimalPeak=std::max(sharpPrimalPeak,std::abs(spa));
                if(mass.phaseVolumes()[u].liquid==0){
                    currentGasPeak=std::max(currentGasPeak,std::abs(ca));sharpGasPeak=std::max(sharpGasPeak,std::abs(sa));
                    currentPrimalGasPeak=std::max(currentPrimalGasPeak,std::abs(cpa));sharpPrimalGasPeak=std::max(sharpPrimalGasPeak,std::abs(spa));
                    if(load.momentumRhsLoad()[u]!=0)++gasCurrent;if(sharpRhs[u]!=0)++gasSharp;
                }
            }
            std::cout<<"CSF_INITIAL_FACE,"<<(circle?"circle":"ellipse")<<','<<u<<','<<point.x<<','<<point.y<<','<<axis<<','<<free[u]<<','
                <<double(box.volume)<<','<<mass.masses()[u]<<','<<double(mass.phaseVolumes()[u].liquid)<<','
                <<load.colorGradient()[u]<<','<<sharpGradient[u]<<',';
            if(load.supportReports()[u].curvatureDefined)std::cout<<load.curvatureAverages()[u];else std::cout<<"NA";
            std::cout<<','<<load.momentumRhsLoad()[u]<<','<<sharpRhs[u]<<','<<ca<<','<<sa<<','<<primalMass[u]<<','<<cpa<<','<<spa<<'\n';
        }
        std::cout<<"CSF_INITIAL_SUMMARY,"<<(circle?"circle":"ellipse")<<','<<n<<','<<phase.geometryIdentity().value()<<','
            <<phase.timeLayer().begin<<','<<phase.timeLayer().end<<','<<color.config().epsilon<<','
            <<currentPeak<<','<<sharpPeak<<','<<currentGasPeak<<','<<sharpGasPeak<<','<<gasCurrent<<','<<gasSharp<<','
            <<current.report().algebraicRelativeResidual<<','<<sharp.report().algebraicRelativeResidual<<','
            <<currentPrimalPeak<<','<<sharpPrimalPeak<<','<<currentPrimalGasPeak<<','<<sharpPrimalGasPeak<<','
            <<currentPrimal.report().algebraicRelativeResidual<<','<<sharpPrimal.report().algebraicRelativeResidual<<'\n';
    }
    std::cout<<"CSF_INITIAL_ACCELERATION_DIAGNOSIS,COMPLETED (no dynamic accuracy claim)\n";
}



// One fixed, zero-flow geometry diagnosis. Fresh-final init is an independent
// control, never a production transfer/resampling operation or accuracy gate.
using RzProbePublication=LAMP::Interface::GeometricVofPublication;
using RzProbeTransfer=LAMP::Interface::GeometricVofMeshTransfer;
std::vector<RzStocks> rzProbeRetained(const RzProbePublication& p)
{std::vector<RzStocks> result;for(const auto& c:p.plicCells())result.push_back(c.inventory);return result;}
RzStocks rzProbeTotals(const RzProbePublication& p)
{
    RzStocks result{},correction{};for(const auto& c:p.plicCells())for(int phase=0;phase<2;++phase){
        const auto x=c.inventory[phase]-correction[phase],next=result[phase]+x;correction[phase]=(next-result[phase])-x;result[phase]=next;}
    return result;
}
RzProbePublication rzProbeZeroEnd(const QuadtreeMesh& mesh,const RzProbePublication& initial)
{
    using namespace LAMP::Interface;
    const LAMP::TimeLayer end{1,0,1e-5,LAMP::TimeLayer::Position::End};
    const std::vector<double> zero(initial.domain().nativeFaceCount(),0);
    const auto slab=GeometricVofAdvectorSlab::prepare(mesh,initial,end,zero,zero);
    const auto transport=GeometricVofTransport::prepareLocalPlicSweep(mesh,initial,slab,1000,1.226);
    check(rzProbeRetained(initial)==rzProbeRetained(transport.endpoint()),"RZ geometry probe zero slab changed raw stocks");
    for(const auto& face:transport.faces())for(const auto* f:{&face.half,&face.full})
        check(f->total==0&&f->liquid==0&&f->gas==0&&f->mass==0,"RZ geometry probe zero slab produced a flux");
    return transport.endpoint();
}
void rzProbeAdapt(QuadtreeMesh& mesh,bool coarsen)
{
    int sweeps=0;while(true){LAMP::Mesh::NativeTargetLevelField target(mesh,0.),gradient(mesh,0.);
        for(int i=0;i<mesh.cellIndexNum();++i){const auto& c=mesh.cells(mesh.cellIndexToLabel[i]);
            const auto lo=mesh.verts(c.cellToVert[0]).coord(),hi=mesh.verts(c.cellToVert[2]).coord();
            const bool inBand=std::max(lo.x,0.)<std::min(hi.x,9./32)&&std::max(lo.y,7./32)<std::min(hi.y,25./32);
            target.value[i]=coarsen?0:inBand?2:0;}
        if(!mesh.mayAdaptToTargets(target))break;check(sweeps<mesh.maxLevel()+2,"RZ geometry probe target failed original allocated-level bound");
        const int before=mesh.cellIndexNum();std::vector<VertField<double>> levelSets;std::vector<VertField<double>*> vertices;
        std::vector<HEdgeField<double>*> horizontal;std::vector<VEdgeField<double>*> vertical;
        std::vector<LAMP::Mesh::NativeTargetLevelField*> cells{&target};
        mesh.adapt(levelSets,gradient,vertices,horizontal,vertical,cells,&target);++sweeps;
        check(mesh.cellIndexNum()!=before,"RZ geometry probe adaptation made no actual progress");
    }
    check(sweeps>0,"RZ geometry probe did not change real topology");
    std::array<std::size_t,3> levels{};for(int i=0;i<mesh.cellIndexNum();++i){const int level=mesh.cells(mesh.cellIndexToLabel[i]).level();
        check(level>=0&&level<=2,"RZ geometry probe escaped0..2");++levels[level];}
    FaceSegmentGraph graph;graph.build(mesh);std::size_t hanging=0;
    for(const auto& f:graph.segments())if(f.L>=0){const int jump=std::abs(mesh.cells(f.K).level()-mesh.cells(f.L).level());
        check(jump<=1,"RZ geometry probe violates real shared-face2:1");hanging+=jump==1;}
    if(coarsen)check(mesh.cellIndexNum()==32*32&&levels[0]==32*32,"RZ geometry probe did not return to actual base grid");
    else check(levels[0]>0&&levels[1]>0&&levels[2]>0&&hanging>0,"RZ geometry probe lacks true graded0/1/2 leaves");
    std::cout<<"RZ_REGRID_HF_MESH,"<<(coarsen?"coarsen":"refine")<<",sweeps,"<<sweeps<<",leaves,"<<mesh.cellIndexNum()
        <<",levels0/1/2,"<<levels[0]<<'/'<<levels[1]<<'/'<<levels[2]<<",actual_hanging_faces,"<<hanging<<'\n';
}
void rzProbeTransferLedger(const RzProbePublication& source,const RzProbeTransfer& transfer)
{
    check(transfer.sourceIdentity()==source.identity(),"RZ geometry transfer has foreign source");
    const auto before=rzProbeTotals(source),after=rzProbeTotals(transfer.candidate());
    for(int phase=0;phase<2;++phase)check(std::abs(after[phase]-before[phase])<=1e-12L*before[phase],"RZ geometry transfer own-phase total changed beyond original gate");
    const auto& r=transfer.report();check(r.localPlic&&r.maximumRelativeDonorPhaseResidual[0]<=1e-12L&&r.maximumRelativeDonorPhaseResidual[1]<=1e-12L,
        "RZ geometry transfer lost original independent donor gate");
    std::cout<<"RZ_REGRID_HF_TRANSFER,refined_donors,"<<r.refinedDonors<<",refined_targets,"<<r.refinedTargets<<",coarsened_targets,"<<r.coarsenedTargets
        <<",source_L,"<<before[0]<<",source_G,"<<before[1]<<",target_L,"<<after[0]<<",target_G,"<<after[1]
        <<",relative_L,"<<std::abs(after[0]-before[0])/before[0]<<",relative_G,"<<std::abs(after[1]-before[1])/before[1]<<'\n';
}
void rzProbeHeight(const char* label,const QuadtreeMesh& mesh,const RzProbePublication& publication,
    const RzProbeTransfer* transfer=nullptr,const RzProbePublication* donor=nullptr)
{
    const auto id=publication.identity();const auto raw=rzProbeRetained(publication);const NativeMeshIdentity identity(mesh);
    const auto metrics=HydroPhysicalDivergence::build(mesh,HydroIntegralDivergence::build(mesh,HydroDofLayout::build(mesh,publication.domain())));
    const auto phase=HydroPhaseGeometry::build(mesh,metrics,publication);const auto pressure=HydroPressureLoad::build(metrics.divergence());
    const auto load=HydroBalancedCsfLoad::localPlicHeightFunctionCapillary(mesh,pressure,phase,.7,heightConfig());
    load.validateVofFor(mesh,pressure,id,publication.timeLayer());rzActualForceAlgebra(mesh,phase,pressure,load,.7,10.);
    check(!transfer||(donor&&transfer->sourceIdentity()==donor->identity()&&transfer->candidate().identity()==id),"RZ geometry probe parent association mismatch");
    long double maximum=0,squared=0,absolute=0,nativeMin=std::numeric_limits<long double>::infinity(),nativeMax=0;
    long double effectiveMin=std::numeric_limits<long double>::infinity(),effectiveMax=0,parentOverChart=0,normalDifference=0;
    double signedMin=std::numeric_limits<double>::infinity(),signedMax=-signedMin;std::size_t mixed=0,negative=0,nativeCharts=0,ancestorCharts=0,otherMeasuredCharts=0,fits=0;
    int worst=-1;const auto& layout=metrics.layout();const auto& jets=load.heightFunctionCellJets();
    for(int p=0;p<layout.pressureCount();++p){const int ci=layout.nativeCell(p);const auto& c=publication.plicCells()[ci];
        if(!(c.resolved&&c.measured[0]>0&&c.measured[1]>0))continue;const auto k=load.heightFunctionCellCurvatures()[p];
        check(std::isfinite(k),"RZ geometry probe actual mixed cell has no finite curvature");++mixed;negative+=k<0;
        signedMin=std::min(signedMin,k);signedMax=std::max(signedMax,k);absolute=std::max(absolute,std::abs(static_cast<long double>(k)));
        const long double error=std::abs(static_cast<long double>(k)-10)/10;squared+=error*error;if(worst<0||error>maximum){maximum=error;worst=p;}
        const long double dx=c.rectangle.upper.x-c.rectangle.lower.x,dy=c.rectangle.upper.y-c.rectangle.lower.y,h=std::max(dx,dy);
        nativeMin=std::min(nativeMin,h);nativeMax=std::max(nativeMax,h);long double chart=h;
        if(jets[p]){const auto& jet=*jets[p];check(jet.normalSpacing>0&&jet.tangentSpacing>0,"RZ geometry probe missing public actual chart spacing");
            chart=std::max(jet.normalSpacing,jet.tangentSpacing);effectiveMin=std::min(effectiveMin,std::min(jet.normalSpacing,jet.tangentSpacing));effectiveMax=std::max(effectiveMax,chart);
            if(jet.coarseAligned)++ancestorCharts;
            else if(jet.normalSpacing==(jet.heightAxis?dy:dx)&&jet.tangentSpacing==(jet.heightAxis?dx:dy))++nativeCharts;
            else ++otherMeasuredCharts;
        }else ++fits; // Existing public successful-fit class, not a claimed HF jet.
        if(transfer)for(auto oi:transfer->nativeOverlaps().byTarget()[ci]){const auto& o=transfer->overlaps()[oi];const auto& old=donor->plicCells()[o.sourceCell];
            const long double oldH=std::max(old.rectangle.upper.x-old.rectangle.lower.x,old.rectangle.upper.y-old.rectangle.lower.y);
            parentOverChart=std::max(parentOverChart,oldH/chart);
            if(old.resolved&&old.measured[0]>0&&old.measured[1]>0)normalDifference=std::max(normalDifference,std::hypot(static_cast<long double>(c.normal.x)-old.normal.x,static_cast<long double>(c.normal.y)-old.normal.y));}
    }
    check(mixed>0&&worst>=0,"RZ geometry probe has no actual mixed sphere cells");const auto total=rzProbeTotals(publication);
    std::cout<<"RZ_REGRID_HF_GEOMETRY,"<<label<<",owner,"<<id.value()<<",leaves,"<<mesh.cellIndexNum()<<",mixed,"<<mixed
        <<",native_h_min/max,"<<nativeMin<<'/'<<nativeMax<<",actual_jet_spacing_min/max,";
    if(nativeCharts+ancestorCharts+otherMeasuredCharts)std::cout<<effectiveMin<<'/'<<effectiveMax;else std::cout<<"NA/NA";
    std::cout<<",kappa_signed_min/max,"<<signedMin<<'/'<<signedMax<<",kappa_abs_max,"<<absolute<<",sphere_exact_kappa,10,relative_error_max,"<<maximum
        <<",relative_error_RMS,"<<std::sqrt(squared/mixed)<<",negative_cells,"<<negative<<",native_chart_estimated,"<<nativeCharts
        <<",explicit_ancestor_chart,"<<ancestorCharts<<",other_public_measured_chart,"<<otherMeasuredCharts<<",fixed_box_fit,"<<fits
        <<",parent_h_over_chart_max,";if(transfer)std::cout<<parentOverChart;else std::cout<<"NA";
    std::cout<<",target_parent_normal_difference_max,";if(transfer)std::cout<<normalDifference;else std::cout<<"NA";
    std::cout<<",stock_L,"<<total[0]<<",stock_G,"<<total[1]<<",HF_extended,"<<load.curvatureReport().heightFunctionExtendedCells
        <<",HF_refined_frame,"<<load.curvatureReport().heightFunctionRefinedCells<<",coarse_aligned_chart_attempts,"<<load.curvatureReport().heightFunctionCoarseAlignedCharts<<'\n';
    const int ci=layout.nativeCell(worst);const auto& c=publication.plicCells()[ci];
    std::cout<<"RZ_REGRID_HF_WORST,"<<label<<",native,"<<ci<<",rectangle,"<<c.rectangle.lower.x<<'/'<<c.rectangle.lower.y<<'/'<<c.rectangle.upper.x<<'/'<<c.rectangle.upper.y
        <<",normal,"<<c.normal.x<<'/'<<c.normal.y<<",offset,"<<c.offset<<",Q,"<<c.inventory[0]<<'/'<<c.inventory[1]<<",kappa,"<<load.heightFunctionCellCurvatures()[worst];
    if(jets[worst]){const auto& jet=*jets[worst];const auto km=HydroBalancedCsfLoad::axisymmetricHeightCurvatures(jet);
        std::cout<<",jet_axis,"<<jet.heightAxis<<",origin,"<<jet.origin<<",tangent_coordinate,"<<jet.tangentCoordinate<<",height/first/second,"<<jet.height<<'/'<<jet.first<<'/'<<jet.second
            <<",normal/tangent_spacing,"<<jet.normalSpacing<<'/'<<jet.tangentSpacing<<",coarse_aligned,"<<jet.coarseAligned<<",meridional/hoop,"<<km[0]<<'/'<<km[1];
    }else std::cout<<",jet,NA_existing_rotated_fit";std::cout<<'\n';
    if(transfer)for(auto oi:transfer->nativeOverlaps().byTarget()[ci]){const auto& o=transfer->overlaps()[oi];const auto& old=donor->plicCells()[o.sourceCell];
        std::cout<<"RZ_REGRID_HF_WORST_PARENT,"<<label<<",native,"<<o.sourceCell<<",kind,"<<static_cast<int>(o.kind)<<",rectangle,"<<old.rectangle.lower.x<<'/'<<old.rectangle.lower.y<<'/'<<old.rectangle.upper.x<<'/'<<old.rectangle.upper.y
            <<",normal,"<<old.normal.x<<'/'<<old.normal.y<<",offset,"<<old.offset<<",actual_phase_correction,"<<o.phaseCorrection[0]<<'/'<<o.phaseCorrection[1]<<",roundoff_tail,"<<o.phaseRoundoffClosure[0]<<'/'<<o.phaseRoundoffClosure[1]<<'\n';}
    identity.validate(mesh);check(publication.identity()==id&&rzProbeRetained(publication)==raw,"RZ geometry probe HF changed true source mesh/stocks/owner");
}
void rzRegridHeightGeometryProbe()
{
    using namespace LAMP::Interface;
    DictionaryMesh input;input.Nx=input.Ny=32;input.minCoord={0,0};input.maxCoord={1,1};input.maxCellNum=65536;input.maxInitialIterations=0;input.maxLevel_local={2};
    QuadtreeMesh coarse(input);GeometricVofInitializationConfig initial;initial.subcellsPerAxis=8;
    const auto phi=[](Vector2D p){const auto z=p.y-.5;return p.x*p.x+z*z-.2*.2;};
    const auto source=rzProbeZeroEnd(coarse,GeometricVofInitialization::prepare(coarse,true,phi,initial));
    const auto sourceStocks=rzProbeRetained(source);const auto sourceId=source.identity();const NativeMeshIdentity sourceIdentity(coarse);
    rzProbeHeight("coarse_init_zero_end",coarse,source);
    auto fine=coarse.clone();rzProbeAdapt(*fine,false);FaceSegmentGraph fineGraph;fineGraph.build(*fine);
    const auto refined=GeometricVofMeshTransfer::prepare(coarse,source,*fine,NativeLeafDomain::build(*fine,fineGraph,true),source.timeLayer());
    rzProbeTransferLedger(source,refined);rzProbeHeight("refine_from_coarse_plane_stocks",*fine,refined.candidate(),&refined,&source);
    const auto fineStocks=rzProbeRetained(refined.candidate());const auto fineId=refined.candidate().identity();const NativeMeshIdentity fineIdentity(*fine);
    const auto control=rzProbeZeroEnd(*fine,GeometricVofInitialization::prepare(*fine,true,phi,initial));
    rzProbeHeight("fresh_final_mesh_diagnostic_control",*fine,control);
    const auto transferredTotal=rzProbeTotals(refined.candidate()),freshTotal=rzProbeTotals(control);
    std::cout<<"RZ_REGRID_HF_CONTROL_VOLUME_DIFFERENCE,fresh_minus_transferred_L,"<<freshTotal[0]-transferredTotal[0]<<",fresh_minus_transferred_G,"<<freshTotal[1]-transferredTotal[1]<<'\n';
    auto back=fine->clone();rzProbeAdapt(*back,true);FaceSegmentGraph backGraph;backGraph.build(*back);
    const auto coarsened=GeometricVofMeshTransfer::prepare(*fine,refined.candidate(),*back,NativeLeafDomain::build(*back,backGraph,true),source.timeLayer());
    rzProbeTransferLedger(refined.candidate(),coarsened);rzProbeHeight("coarsen_after_refine",*back,coarsened.candidate(),&coarsened,&refined.candidate());
    sourceIdentity.validate(coarse);fineIdentity.validate(*fine);
    check(source.identity()==sourceId&&rzProbeRetained(source)==sourceStocks&&refined.candidate().identity()==fineId&&rzProbeRetained(refined.candidate())==fineStocks,
        "RZ geometry comparison changed its source owner/stock/topology");
    std::cout<<"RZ_REGRID_HF_GEOMETRY_DIAGNOSTIC_COMPLETED,no_NS_solve,no_accuracy_PASS,no_dynamic_mesh_replay_claim\n";
}


// D22: a force closure is authorized by the SAME immutable PLIC publication.
// These are independently declared input stocks; below-resolution neighbours
// retain their nonzero minor stocks while their actual declared geometry is pure.
struct IsolatedFragmentStocks
{
    std::vector<std::array<long double,2>> values;
    int source=-1;
};
IsolatedFragmentStocks isolatedFragmentStocks(const QuadtreeMesh& mesh,
    const HydroPhysicalDivergence& metrics,unsigned phase,int variant=0)
{
    check(phase<2,"isolated-fragment fixture phase");
    const auto& layout=metrics.layout();constexpr double h=1./16;
    const Vector2D lower{variant==3?0.:.5,.5},center=lower+Vector2D{h/2,h/2};
    IsolatedFragmentStocks out;out.values.resize(mesh.cellIndexNum());
    for(int ci=0;ci<mesh.cellIndexNum();++ci){
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[ci]);
        const auto lo=mesh.verts(cell.cellToVert[0]).coord(),hi=mesh.verts(cell.cellToVert[2]).coord();
        check(hi.x-lo.x==h&&hi.y-lo.y==h,"isolated-fragment fixture must retain actual uniform native cells");
        const auto site=cell.centerCoord();const int x=static_cast<int>(std::llround((site.x-center.x)/h));
        const int y=static_cast<int>(std::llround((site.y-center.y)/h));
        long double fraction=0;
        if(x==0&&y==0){fraction=1.L/32;out.source=ci;}
        else if(std::abs(x)<=1&&std::abs(y)<=1)fraction=std::ldexp(static_cast<long double>(4+x+2*y),-50);
        if((variant==1&&x==1&&y==0)||(variant==2&&x==1&&y==1))fraction=1;
        const auto capacity=metrics.cellVolumes()[layout.pressureDof(ci)];
        // Each phase is supplied explicitly; no returned stock/plane is changed.
        const long double minor=capacity*fraction,major=capacity*(1.L-fraction);
        out.values[ci]=phase==0?std::array<long double,2>{minor,major}:std::array<long double,2>{major,minor};
    }
    check(out.source>=0,"isolated-fragment fixture lacks its actual source rectangle");return out;
}
std::vector<std::array<long double,8>> isolatedFragmentSnapshot(
    const LAMP::Interface::GeometricVofPublication& publication)
{
    std::vector<std::array<long double,8>> result;
    for(const auto& p:publication.plicCells())result.push_back({p.inventory[0],p.inventory[1],
        p.measured[0],p.measured[1],p.normal.x,p.normal.y,p.offset,p.resolved?1.L:0.L});
    return result;
}
template<class F>void isolatedFragmentReject(F f,const char* description)
{
    bool rejected=false;try{f();}catch(const std::exception&){rejected=true;}
    check(rejected,description);
}
void localPlicIsolatedFragmentClosure()
{
    using namespace LAMP::Interface;
    constexpr double sigma=.7;const GeometricVofTimeLayer clock{22,.01,.0102,Position::Midpoint};
    std::size_t positive=0,negative=0,sourceGuards=0,activeClosureFaces=0;
    for(unsigned minorPhase:{0u,1u}){
        auto mesh=heightMesh(16);const auto metrics=metricsFor(*mesh,false);const auto& layout=metrics.layout();
        const auto input=isolatedFragmentStocks(*mesh,metrics,minorPhase);
        const auto publication=GeometricVofPublication::prepareLocalPlicInventory(*mesh,layout.domain(),clock,input.values,{});
        const auto id=publication.identity();const auto before=isolatedFragmentSnapshot(publication);
        const auto& cell=publication.plicCells()[input.source];
        std::cout<<"D22_ISOLATED_FIXTURE,phase,"<<minorPhase<<",native,"<<input.source<<",normal,"<<cell.normal.x<<'/'<<cell.normal.y
            <<",offset,"<<cell.offset<<",resolved,"<<cell.resolved<<",stock,"<<cell.inventory[0]<<'/'<<cell.inventory[1]<<'\n';
        check(cell.resolved&&cell.measured[0]>0&&cell.measured[1]>0,
            "isolated-fragment witness must be actual resolved geometry, not only a raw micro stock");
        check(std::abs(cell.normal.x)>.05&&std::abs(cell.normal.y)>.05,
            "isolated-fragment fixture did not produce its intended actual diagonal plane");
        const auto certificate=publication.isolatedLocalPlicFragment(*mesh,input.source);
        check(certificate.has_value(),"closed subcell component did not receive owner-issued qualification");
        certificate->validateFor(*mesh,publication);
        const auto time=certificate->sourceTime();
        check(certificate->sourceIdentity()==id&&certificate->nativeCell()==input.source&&certificate->phase()==minorPhase,
            "isolated-fragment certificate changed its owner/cell/phase");
        check(time.step==clock.step&&time.begin==clock.begin&&time.end==clock.end&&time.position==clock.position,
            "isolated-fragment certificate lost its full source clock");
        const long double delta=std::min(cell.rectangle.upper.x-cell.rectangle.lower.x,cell.rectangle.upper.y-cell.rectangle.lower.y);
        check(certificate->delta()==delta&&certificate->diameter()>0&&certificate->diameter()<delta,
            "isolated-fragment certificate lacks its strict actual diameter bound");
        check(certificate->coveredNativeFaces()>0&&certificate->checkedCornerCells()>0,
            "isolated-fragment qualification omitted real native-face/corner checks");
        // The new certificate is neither transferable to a new owner nor to a
        // new full clock, even when every supplied stock is unchanged.
        const auto other=GeometricVofPublication::prepareLocalPlicInventory(*mesh,layout.domain(),clock,input.values,{});
        isolatedFragmentReject([&]{certificate->validateFor(*mesh,other);},"isolated-fragment foreign owner accepted");++sourceGuards;
        auto nextClock=clock;nextClock.end+=.0001;
        const auto next=GeometricVofPublication::prepareLocalPlicInventory(*mesh,layout.domain(),nextClock,input.values,{});
        isolatedFragmentReject([&]{certificate->validateFor(*mesh,next);},"isolated-fragment foreign full-clock owner accepted");++sourceGuards;
        isolatedFragmentReject([&]{publication.validateFor(*mesh,id,nextClock);},"publication full-clock mismatch accepted");++sourceGuards;
        auto foreignMesh=heightMesh(32); // A different actual native geometry/connectivity identity.
        isolatedFragmentReject([&]{certificate->validateFor(*foreignMesh,publication);},"isolated-fragment foreign mesh accepted");++sourceGuards;
        const auto phase=HydroPhaseGeometry::build(*mesh,metrics,publication);
        const auto pressure=HydroPressureLoad::build(metrics.divergence());
        const auto load=HydroBalancedCsfLoad::localPlicHeightFunctionCapillary(*mesh,pressure,phase,sigma,heightConfig());
        load.validateVofFor(*mesh,pressure,id,clock);const int p=layout.pressureDof(input.source);
        check(load.heightFunctionCellCurvatureMethods()[p]==HydroCsfCellCurvatureMethod::UnderresolvedIsolatedFragmentZero,
            "isolated-fragment zero force lost its explicit cell-method tag");
        check(std::isnan(load.heightFunctionCellCurvatures()[p])&&load.heightFunctionCellEffectiveCurvatures()[p]==0,
            "isolated-fragment force closure masqueraded as a physical zero curvature");
        check(!load.heightFunctionCellJets()[p],"isolated-fragment closure fabricated a qualified height jet");
        check(load.curvatureReport().unresolvedIsolatedFragmentCells==1&&load.curvatureReport().heightFunctionCells==0&&
            load.curvatureReport().heightFunctionFitCells==0,"isolated-fragment closure mixed its counter with successful HF/fit estimates");
        std::vector<double> measuredColor(layout.pressureCount());
        for(int q=0;q<layout.pressureCount();++q){const auto volumes=phase.cellVolumes()[q];
            measuredColor[q]=volumes.liquid==0?0:volumes.gas==0?1:static_cast<double>(volumes.liquid/metrics.cellVolumes()[q]);}
        const Reference independent(*mesh,phase);const auto gradient=independent.gradient(measuredColor);
        double error=0,scale=0;std::size_t active=0;
        for(int u=0;u<layout.velocityCount();++u){error=std::hypot(error,gradient[u]-load.colorGradient()[u]);scale=std::hypot(scale,gradient[u]);
            check(std::isfinite(load.momentumRhsLoad()[u])&&load.momentumRhsLoad()[u]==0,
                "isolated-fragment effective zero produced a nonzero/nonfinite actual load");
            const auto face=layout.domain().face(*mesh,layout.nativeFace(u));
            if(load.colorGradient()[u]!=0&&(face.negativeCell==input.source||face.positiveCell==input.source)){
                ++active;check(load.supportReports()[u].underresolvedFragmentClosure,
                    "active isolated-fragment face lost the explicit closure provenance");
                check(load.heightFunctionFaceCurvatures()[u]==0,"active isolated-fragment face multiplier is not the declared effective zero");}
        }
        check(active>0&&relative(error,scale)<=1e-12,"isolated-fragment closure changed same-plane measured color/native incidence");
        activeClosureFaces+=active;
        check(publication.identity()==id&&isolatedFragmentSnapshot(publication)==before,
            "isolated-fragment query/force changed original micro stocks, planes or geometry");
        for(std::size_t ci=0;ci<input.values.size();++ci)check(publication.plicCells()[ci].inventory==input.values[ci],
            "isolated-fragment closure erased or rewrote an independent input stock");
        ++positive;std::cout<<"D22_ISOLATED_CLOSURE,phase,"<<minorPhase<<",diameter,"<<certificate->diameter()<<",delta,"<<delta
            <<",covered_native_faces,"<<certificate->coveredNativeFaces()<<",corner_cells_checked,"<<certificate->checkedCornerCells()
            <<",active_faces,"<<active<<",physical_kappa,NA,effective_kappa,0,method,UnderresolvedIsolatedFragmentZero\n";
    }
    // Each negative uses the genuine public inventory/normal reconstruction.
    // Connected-face and corner-only contact are not absence certificates.
    for(int variant:{1,2,3,4}){
        auto mesh=heightMesh(16);const bool rz=variant==4;const auto metrics=metricsFor(*mesh,rz);const auto& layout=metrics.layout();
        const auto input=isolatedFragmentStocks(*mesh,metrics,0,variant==4?0:variant);
        const auto publication=GeometricVofPublication::prepareLocalPlicInventory(*mesh,layout.domain(),clock,input.values,{});
        const auto before=isolatedFragmentSnapshot(publication);const auto& cell=publication.plicCells()[input.source];
        check(cell.resolved&&cell.measured[0]>0&&cell.measured[1]>0,"isolated-fragment negative fixture lacks real measured mixed geometry");
        check(!publication.isolatedLocalPlicFragment(*mesh,input.source),
            "connected/wall/corner/RZ geometry received an isolated XY absence certificate");
        check(isolatedFragmentSnapshot(publication)==before,"negative isolated-fragment query modified its source");
        ++negative;std::cout<<"D22_ISOLATED_REJECT,variant,"<<variant<<",native,"<<input.source<<",normal,"<<cell.normal.x<<'/'<<cell.normal.y
            <<",measured,"<<cell.measured[0]<<'/'<<cell.measured[1]<<"\n";
    }
    check(positive==2&&negative==4&&sourceGuards==8&&activeClosureFaces>0,"isolated-fragment evidence incomplete");
    std::cout<<"D22_ISOLATED_FRAGMENT_PASS,positive,"<<positive<<",negative,"<<negative<<",owner_clock_mesh_guards,"<<sourceGuards
        <<",active_closure_faces,"<<activeClosureFaces<<"\n";
}
// D23: a two-bin position deficit, NOT a singleton diameter certificate.
// The public inventory constructor owns the actual PLIC normal/offset.
struct TwoCellFragmentStocks
{
    std::vector<std::array<long double,2>> values;
    std::array<int,2> cells{{-1,-1}};
    int witness=-1;
};
TwoCellFragmentStocks twoCellFragmentStocks(const QuadtreeMesh& mesh,
    const HydroPhysicalDivergence& metrics,unsigned phase,int variant=0)
{
    check(phase<2,"two-cell fragment fixture phase");
    constexpr double h=1./16;const double x0=variant==3?0.:.5,y0=.5;
    TwoCellFragmentStocks out;out.values.resize(mesh.cellIndexNum());
    const auto& layout=metrics.layout();
    for(int ci=0;ci<mesh.cellIndexNum();++ci){
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[ci]);
        const auto lo=mesh.verts(cell.cellToVert[0]).coord(),hi=mesh.verts(cell.cellToVert[2]).coord();
        check(hi.x-lo.x==h&&hi.y-lo.y==h,"two-cell fragment fixture must retain actual uniform native cells");
        const int x=static_cast<int>(std::llround((lo.x-x0)/h)),y=static_cast<int>(std::llround((lo.y-y0)/h));
        long double fraction=0;
        if(x==0&&(y==0||y==1)){
            fraction=variant==1?3.L/8:std::ldexp(1.L,y==0?-10:-11);
            out.cells[y]=ci;
        }
        // Genuine exterior face connection, diagonal corner connection, or
        // a branch that no longer has matching shared-face traces.
        if((variant==2&&x==1&&y==0)||(variant==4&&x==1&&y==2)||(variant==5&&x==0&&y==-1)){
            fraction=1;out.witness=ci;
        }
        const long double capacity=metrics.cellVolumes()[layout.pressureDof(ci)];
        const long double minor=capacity*fraction,major=capacity*(1.L-fraction);
        out.values[ci]=phase==0?std::array<long double,2>{{minor,major}}:std::array<long double,2>{{major,minor}};
    }
    check(out.cells[0]>=0&&out.cells[1]>=0,"two-cell fixture lacks actual neighboring source leaves");
    return out;
}
using TwoCellPoint=std::array<long double,2>;
// Independent one-plane polygon clip. No publication topology/packing helper,
// production box integrator or certificate internals enter this reference.
std::vector<TwoCellPoint> twoCellActualPolygon(const LAMP::Interface::GeometricVofPublication& publication,
    int native,unsigned phase)
{
    const auto& cell=publication.plicCells()[native];const auto& b=cell.rectangle;
    if(cell.measured[phase]==0)return {};
    std::vector<TwoCellPoint> polygon{{b.lower.x,b.lower.y},{b.upper.x,b.lower.y},{b.upper.x,b.upper.y},{b.lower.x,b.upper.y}};
    if(cell.measured[1-phase]==0)return polygon;
    const long double cx=(static_cast<long double>(b.lower.x)+b.upper.x)/2,cy=(static_cast<long double>(b.lower.y)+b.upper.y)/2;
    const auto level=[&](TwoCellPoint p){return (phase==0?1.L:-1.L)*
        (static_cast<long double>(cell.normal.x)*(p[0]-cx)+static_cast<long double>(cell.normal.y)*(p[1]-cy)-cell.offset);};
    std::vector<TwoCellPoint> clipped;auto a=polygon.back();auto fa=level(a);
    for(auto z:polygon){const auto fz=level(z);
        if((fa<=0)!=(fz<=0)){const auto t=fa/(fa-fz);clipped.push_back({a[0]+t*(z[0]-a[0]),a[1]+t*(z[1]-a[1])});}
        if(fz<=0)clipped.push_back(z);a=z;fa=fz;}
    return clipped;
}
std::array<long double,4> twoCellActualBounds(const LAMP::Interface::GeometricVofPublication& publication,
    const std::array<int,2>& cells,unsigned phase)
{
    std::array<long double,4> bounds{{std::numeric_limits<long double>::infinity(),std::numeric_limits<long double>::infinity(),
        -std::numeric_limits<long double>::infinity(),-std::numeric_limits<long double>::infinity()}};
    for(int ci:cells){const auto polygon=twoCellActualPolygon(publication,ci,phase);
        check(polygon.size()>=3,"two-cell actual phase polygon is empty/degenerate");
        for(const auto& p:polygon){bounds[0]=std::min(bounds[0],p[0]);bounds[1]=std::min(bounds[1],p[1]);
            bounds[2]=std::max(bounds[2],p[0]);bounds[3]=std::max(bounds[3],p[1]);}}
    return bounds;
}
bool twoCellContainsCorner(const LAMP::Interface::GeometricVofPublication& publication,int native,
    unsigned phase,TwoCellPoint point)
{
    const auto& c=publication.plicCells()[native];
    if(c.measured[phase]==0)return false;if(c.measured[1-phase]==0)return true;
    const long double cx=(static_cast<long double>(c.rectangle.lower.x)+c.rectangle.upper.x)/2;
    const long double cy=(static_cast<long double>(c.rectangle.lower.y)+c.rectangle.upper.y)/2;
    return (phase==0?1.L:-1.L)*(static_cast<long double>(c.normal.x)*(point[0]-cx)+
        static_cast<long double>(c.normal.y)*(point[1]-cy)-c.offset)<=0;
}
void localPlicTwoCellFragmentClosure()
{
    using namespace LAMP::Interface;
    constexpr double sigma=.7,h=1./16;
    const GeometricVofTimeLayer clock{23,.02,.0202,Position::Midpoint};
    unsigned positive=0,queries=0,negative=0,ownerGuards=0,activeClosureFaces=0;
    long double maximumPackingRatio=0;
    for(unsigned phaseIndex:{0u,1u}){
        auto mesh=heightMesh(16);const auto metrics=metricsFor(*mesh,false);const auto& layout=metrics.layout();
        const auto input=twoCellFragmentStocks(*mesh,metrics,phaseIndex);
        const auto publication=GeometricVofPublication::prepareLocalPlicInventory(*mesh,layout.domain(),clock,input.values,{});
        const auto id=publication.identity();const auto before=isolatedFragmentSnapshot(publication);
        const auto alpha=publication.alpha();const auto* storage=publication.plicCells().data();
        const auto actual=twoCellActualBounds(publication,input.cells,phaseIndex);
        const auto width=actual[2]-actual[0],height=actual[3]-actual[1];
        check(width>0&&height>0&&width/2+height<h&&std::hypot(width,height)>h,
            "two-cell fixture must satisfy two-bin packing while FAILING the old singleton diameter condition");
        const auto& lower=publication.plicCells()[input.cells[0]];const auto& upper=publication.plicCells()[input.cells[1]];
        check(lower.resolved&&upper.resolved&&lower.measured[0]>0&&lower.measured[1]>0&&upper.measured[0]>0&&upper.measured[1]>0,
            "two-cell positive lacks two true resolved mixed source planes");
        check(lower.rectangle.upper.y==upper.rectangle.lower.y&&lower.rectangle.lower.x==upper.rectangle.lower.x
            &&lower.rectangle.upper.x==upper.rectangle.upper.x&&std::abs(lower.normal.x)<1e-12&&std::abs(upper.normal.x)<1e-12
            &&lower.normal.y*upper.normal.y<-.99,"two-cell positive is not the actual opposite horizontal cap pair");
        const auto y=static_cast<long double>(lower.rectangle.upper.y);
        for(const auto x:{lower.rectangle.lower.x,lower.rectangle.upper.x})
            check(twoCellContainsCorner(publication,input.cells[0],phaseIndex,{x,y})
                &&twoCellContainsCorner(publication,input.cells[1],phaseIndex,{x,y}),
                "two-cell positive lacks actual full-side common phase trace");
        std::array<int,2> canonical=input.cells;if(canonical[0]>canonical[1])std::swap(canonical[0],canonical[1]);
        const auto& anchor=publication.plicCells()[canonical[0]].rectangle;
        const std::array<long double,2> origin{{(static_cast<long double>(anchor.lower.x)+anchor.upper.x)/2,
            (static_cast<long double>(anchor.lower.y)+anchor.upper.y)/2}};
        std::array<long double,2> firstLower{},firstUpper{};unsigned firstAxis=0;
        for(int native:input.cells){
            check(!publication.isolatedLocalPlicFragment(*mesh,native),"two-cell component falsely passed original singleton qualification");
            const auto proof=publication.twoCellLocalPlicPositionDeficit(*mesh,native);
            check(proof.has_value(),"two-cell closed subgrid stripe did not receive actual owner-issued packing proof");
            proof->validateFor(*mesh,publication);
            check(proof->sourceIdentity()==id&&proof->nativeCell()==native&&proof->nativeCells()==canonical&&proof->phase()==phaseIndex,
                "two-cell proof lost requested source/canonical component/phase/owner");
            const auto time=proof->sourceTime();
            check(time.step==clock.step&&time.begin==clock.begin&&time.end==clock.end&&time.position==clock.position,
                "two-cell proof changed the full source interval/position");
            const auto lo=proof->bboxLower(),hi=proof->bboxUpper();const auto dx=hi[0]-lo[0],dy=hi[1]-lo[1];
            check(dx>0&&dy>0&&proof->delta()==h&&proof->maximumIndependentPositions()==2&&proof->splitAxis()<2,
                "two-cell proof changed actual delta or independent-position interpretation");
            const auto bound=proof->splitAxis()==0?dx/2+dy:dx+dy/2;
            check(proof->twoBinL1Bound()>0&&proof->twoBinL1Bound()<h&&proof->twoBinL1Bound()>=bound,
                "two-cell proof relaxed strict L1 packing or returned a downward bound");
            const auto allowance=32*std::numeric_limits<long double>::epsilon()*h;
            check(lo[0]+origin[0]<=actual[0]+allowance&&lo[1]+origin[1]<=actual[1]+allowance
                &&hi[0]+origin[0]>=actual[2]-allowance&&hi[1]+origin[1]>=actual[3]-allowance
                &&std::abs(dx-width)<=allowance&&std::abs(dy-height)<=allowance,
                "owner AABB does not enclose the independently clipped SAME actual PLIC component");
            check(proof->coveredNativeFaces()>=7&&proof->checkedCornerCells()>0,
                "two-cell proof omitted actual perimeter/shared-face/corner incidence");
            if(queries%2==0){firstLower=lo;firstUpper=hi;firstAxis=proof->splitAxis();}
            else check(firstLower==lo&&firstUpper==hi&&firstAxis==proof->splitAxis(),
                "two source queries disagree on canonical component coordinates or packing");
            maximumPackingRatio=std::max(maximumPackingRatio,proof->twoBinL1Bound()/h);
            const auto other=GeometricVofPublication::prepareLocalPlicInventory(*mesh,layout.domain(),clock,input.values,{});
            isolatedFragmentReject([&]{proof->validateFor(*mesh,other);},"two-cell proof accepted foreign same-numeric owner");++ownerGuards;
            auto nextClock=clock;++nextClock.step;
            const auto next=GeometricVofPublication::prepareLocalPlicInventory(*mesh,layout.domain(),nextClock,input.values,{});
            isolatedFragmentReject([&]{proof->validateFor(*mesh,next);},"two-cell proof accepted foreign fullclock owner");++ownerGuards;
            isolatedFragmentReject([&]{publication.validateFor(*mesh,id,nextClock);},"two-cell publication ignored equal-time step mismatch");++ownerGuards;
            auto foreignMesh=heightMesh(32);
            isolatedFragmentReject([&]{proof->validateFor(*foreignMesh,publication);},"two-cell proof accepted different native mesh");++ownerGuards;
            ++queries;
        }
        const auto phase=HydroPhaseGeometry::build(*mesh,metrics,publication);
        const auto pressure=HydroPressureLoad::build(metrics.divergence());
        const auto load=HydroBalancedCsfLoad::localPlicHeightFunctionCapillary(*mesh,pressure,phase,sigma,heightConfig());
        load.validateVofFor(*mesh,pressure,id,clock);
        check(load.curvatureReport().unresolvedIsolatedFragmentCells==2&&load.curvatureReport().heightFunctionCells==0
            &&load.curvatureReport().heightFunctionFitCells==0,"two-cell closure masqueraded as successful physical HF/fit");
        for(int native:input.cells){const int p=layout.pressureDof(native);
            check(load.heightFunctionCellCurvatureMethods()[p]==HydroCsfCellCurvatureMethod::UnderresolvedIsolatedFragmentZero
                &&std::isnan(load.heightFunctionCellCurvatures()[p])&&load.heightFunctionCellEffectiveCurvatures()[p]==0
                &&!load.heightFunctionCellJets()[p],"two-cell closure lost method3 effective zero / physical NA distinction");}
        std::vector<double> color(layout.pressureCount());
        for(int p=0;p<layout.pressureCount();++p){const auto c=phase.cellVolumes()[p];color[p]=c.liquid==0?0:c.gas==0?1:static_cast<double>(c.liquid/metrics.cellVolumes()[p]);}
        const Reference independent(*mesh,phase);const auto gradient=independent.gradient(color);
        double error=0,scale=0;unsigned active=0;
        for(int u=0;u<layout.velocityCount();++u){
            error=std::hypot(error,gradient[u]-load.colorGradient()[u]);scale=std::hypot(scale,gradient[u]);
            check(std::isfinite(load.momentumRhsLoad()[u])&&load.momentumRhsLoad()[u]==0,
                "two-cell effective zero created nonzero/nonfinite actual capillary load");
            const auto face=layout.domain().face(*mesh,layout.nativeFace(u));
            const bool adjacent=std::find(input.cells.begin(),input.cells.end(),face.negativeCell)!=input.cells.end()
                ||std::find(input.cells.begin(),input.cells.end(),face.positiveCell)!=input.cells.end();
            if(adjacent&&load.colorGradient()[u]!=0){++active;check(load.supportReports()[u].underresolvedFragmentClosure
                &&load.heightFunctionFaceCurvatures()[u]==0,"two-cell active face lost effective closure provenance");}
        }
        check(active>0&&relative(error,scale)<=1e-12,"two-cell closure changed original same-plane measured color incidence");
        activeClosureFaces+=active;
        check(publication.identity()==id&&publication.plicCells().data()==storage&&publication.alpha()==alpha
            &&isolatedFragmentSnapshot(publication)==before,"two-cell proof/CSF mutated source owner/stocks/normal/offset");
        for(std::size_t ci=0;ci<input.values.size();++ci)check(publication.plicCells()[ci].inventory==input.values[ci],
            "two-cell effective closure erased or repaired an independent stock");
        ++positive;std::cout<<"D23_TWO_CELL_CLOSURE,phase,"<<phaseIndex<<",sources,"<<input.cells[0]<<'/'<<input.cells[1]
            <<",width,"<<width<<",height,"<<height<<",independent_split_L1,"<<width/2+height<<",delta,"<<h
            <<",physical_kappa,NA,effective_kappa,0,active_faces,"<<active<<'\n';
    }
    for(int variant:{1,2,3,4,5,6}){
        auto mesh=heightMesh(16);const bool rz=variant==6;const auto metrics=metricsFor(*mesh,rz);
        const auto input=twoCellFragmentStocks(*mesh,metrics,0,variant==6?0:variant);
        const auto publication=GeometricVofPublication::prepareLocalPlicInventory(*mesh,metrics.layout().domain(),clock,input.values,{});
        const auto before=isolatedFragmentSnapshot(publication);
        if(variant==1){const auto b=twoCellActualBounds(publication,input.cells,0);const auto dx=b[2]-b[0],dy=b[3]-b[1];
            check(dx/2+dy>=h&&dx+dy/2>=h,"thick-stripe negative does not actually violate BOTH strict packing alternatives");}
        if(variant==2||variant==4||variant==5)check(input.witness>=0&&publication.plicCells()[input.witness].measured[0]>0,
            "two-cell topology negative lacks its actual external same-phase source");
        if(variant==2){
            const auto polygon=twoCellActualPolygon(publication,input.cells[0],0);
            const auto x=static_cast<long double>(publication.plicCells()[input.cells[0]].rectangle.upper.x);
            long double lo=std::numeric_limits<long double>::infinity(),hi=-lo;
            for(const auto& p:polygon)if(p[0]==x){lo=std::min(lo,p[1]);hi=std::max(hi,p[1]);}
            check(hi>lo,"third-cell negative lacks a positive actual shared-face phase trace");
        }
        if(variant==5){const auto& c=publication.plicCells()[input.cells[0]].rectangle;
            const TwoCellPoint mid{{(static_cast<long double>(c.lower.x)+c.upper.x)/2,c.upper.y}};
            check(twoCellContainsCorner(publication,input.cells[0],0,mid)!=twoCellContainsCorner(publication,input.cells[1],0,mid),
                "shared-trace negative does not actually disagree on its common face");
        }
        if(variant==4){const auto& c=publication.plicCells()[input.cells[1]].rectangle;
            check(twoCellContainsCorner(publication,input.cells[1],0,{c.upper.x,c.upper.y}),
                "corner negative does not actually contact the declared outside diagonal at the source corner");}
        for(int native:input.cells)check(!publication.twoCellLocalPlicPositionDeficit(*mesh,native),
            "strict-packing/third-cell/wall/corner/mismatching-trace/RZ negative received a two-cell absence certificate");
        check(isolatedFragmentSnapshot(publication)==before,"negative two-cell proof query changed native source");
        ++negative;std::cout<<"D23_TWO_CELL_REJECT,variant,"<<variant<<",sources,"<<input.cells[0]<<'/'<<input.cells[1]<<'\n';
    }
    check(positive==2&&queries==4&&negative==6&&ownerGuards==16&&activeClosureFaces>0,"D23 two-cell evidence incomplete");
    std::cout<<"D23_TWO_CELL_FRAGMENT_PASS,positive,"<<positive<<",source_queries,"<<queries<<",negative,"<<negative
        <<",owner_clock_mesh_guards,"<<ownerGuards<<",maximum_two_bin_ratio,"<<maximumPackingRatio
        <<",active_closure_faces,"<<activeClosureFaces<<",NO_DYNAMIC_ACCURACY_CLAIM\n";
}
namespace F7SouthEven90
{
using Boundary=HydroVofSouthEven90Boundary;
using Vof=LAMP::Interface::GeometricVofPublication;
using Clock=LAMP::TimeLayer;
bool f7ClockEqual(Clock a,Clock b){return a.step==b.step&&a.begin==b.begin&&a.end==b.end&&a.position==b.position;}
unsigned directScenes=0,boundaryGuards=0,straightFaces=0,splitSouth=0,providerStages=0;
double maximumXY=0,maximumRZ=0,constantFree=0,constantReaction=0;
std::vector<int> southRows(const QuadtreeMesh& mesh,const Vof& pub)
{
    std::vector<int> result;
    for(int row=0;row<pub.domain().nativeFaceCount();++row){const auto f=pub.domain().face(mesh,row);
        if(f.horizontal&&f.negativeCell<0&&f.area>0&&f.center.y==mesh.minCoord.y)result.push_back(row);}
    check(!result.empty(),"F7 actual south positive-area rows missing");return result;
}
template<class F>void reject(F action,const char* why)
{
    bool caught=false;try{action();}catch(const std::exception&){caught=true;}
    check(caught,why);++boundaryGuards;
}
long double xyLiquid(Vector2D lo,Vector2D hi,bool column)
{
    constexpr long double radius=.25L,center=.5L;const long double capacity=static_cast<long double>(hi.x-lo.x)*(hi.y-lo.y);
    const long double ycenter=column?.5L:0.L;
    long double liquid=0;
    if(column&&lo.y<.5){
        const long double upper=std::min(static_cast<long double>(hi.y),.5L),left=std::max(static_cast<long double>(lo.x),center-radius),right=std::min(static_cast<long double>(hi.x),center+radius);
        if(upper>lo.y&&right>left)liquid+=(upper-lo.y)*(right-left);}
    if(hi.y>ycenter){
        const Vector2D lower{lo.x,static_cast<double>(std::max(static_cast<long double>(lo.y),ycenter))};
        const long double far=std::max(std::abs(static_cast<long double>(lo.x)-center),std::abs(static_cast<long double>(hi.x)-center));
        const long double zmax=hi.y-ycenter;
        if(far*far+zmax*zmax<=radius*radius)liquid+=static_cast<long double>(hi.x-lo.x)*(hi.y-lower.y);
        else liquid+=circleCell(lower,hi,{static_cast<double>(center),static_cast<double>(ycenter)},static_cast<double>(radius),false);}
    check(liquid>=0&&liquid<=capacity,"F7 independent half-circle/column XY integral");return liquid;
}
std::vector<RzStocks> inputStocks(const QuadtreeMesh& mesh,bool rz,bool column)
{
    std::vector<RzStocks> result;
    for(int ci=0;ci<mesh.cellIndexNum();++ci){const auto& c=mesh.cells(mesh.cellIndexToLabel[ci]);
        const auto lo=mesh.verts(c.cellToVert[0]).coord(),hi=mesh.verts(c.cellToVert[2]).coord();
        if(rz)result.push_back(column?rzCapsuleColumn(lo,hi,0,.5,.25):rzSphereColumn(lo,hi,0,.25));
        else{const auto liquid=xyLiquid(lo,hi,column);const long double capacity=static_cast<long double>(hi.x-lo.x)*(hi.y-lo.y);result.push_back({liquid,capacity-liquid});}}
    return result;
}
void algebra(const QuadtreeMesh& mesh,const HydroPhaseGeometry& phase,const HydroPressureLoad& pressure,
    const HydroBalancedCsfLoad& load,double expected)
{
    const auto& l=phase.layout();const Reference ref(mesh,phase);std::vector<double> color(l.pressureCount());
    for(int p=0;p<l.pressureCount();++p){const auto& cell=phase.vofPublication().plicCells()[l.nativeCell(p)];
        color[p]=cell.measured[0]==0?0:cell.measured[1]==0?1:static_cast<double>(cell.measured[0]/ref.volume[p]);}
    const auto gradient=ref.gradient(color);double eg=0,sg=0,ef=0,sf=0;
    for(int u=0;u<l.velocityCount();++u){
        eg=std::hypot(eg,load.colorGradient()[u]-gradient[u]);sg=std::hypot(sg,gradient[u]);
        const auto& support=load.supportReports()[u];const double want=support.curvatureDefined?.7*load.curvatureAverages()[u]*gradient[u]:0;
        ef=std::hypot(ef,load.momentumRhsLoad()[u]-want);sf=std::hypot(sf,want);
        if(ref.area[u]==0)check(load.colorGradient()[u]==0&&load.momentumRhsLoad()[u]==0,"F7 true axis C/CSF not zero");
    }
    check(relative(eg,sg)<=1e-12&&relative(ef,sf)<=1e-12,"F7 original independent C incidence/actual HF force algebra");
    // Constant-kappa reference is an ORIGINAL C mode, not a claim that the
    // approximate half-circle estimator has exact pointwise constant curvature.
    auto potential=color;for(auto& p:potential)p*=.7*expected;
    const auto rhs=pressure.apply(mesh,HydroNativeValues(l,std::vector<double>(l.velocityCount()),potential));
    const auto means=ref.meanByCell(color),meanGradient=ref.gradient(means);
    Matrix A(l.velocityCount(),l.velocityCount());A.setIdentity();
    const auto saddle=HydroSinglePhaseSaddle::buildClosed(pressure,HydroSaddleVelocityBlock::build(l,A));
    const auto result=saddle.solve(mesh,HydroSaddleRhs(l,rhs,std::vector<double>(l.velocityCount())));
    const auto cp=ref.gradient(result.values().pressure());std::vector<bool> free(l.velocityCount());
    for(int u:saddle.freeVelocityDofs())free[u]=true;
    double fe=0,fs=0,re=0,rs=0;bool reactionWitness=false;
    for(int u=0;u<l.velocityCount();++u){
        if(free[u]){fe=std::hypot(fe,cp[u]-rhs[u]);fs=std::hypot(fs,std::abs(cp[u])+std::abs(rhs[u]));}
        else{const double actual=result.values().normalVelocity()[u]+cp[u]-rhs[u],want=-.7*expected*meanGradient[u];
            re=std::hypot(re,actual-want);rs=std::hypot(rs,std::abs(actual)+std::abs(want));reactionWitness=reactionWitness||want!=0;}}
    constantFree=std::max(constantFree,relative(fe,fs));constantReaction=std::max(constantReaction,relative(re,rs));
    check(relative(fe,fs)<=1e-12&&relative(re,rs)<=1e-12,"F7 original constant C free equilibrium/fixed reaction");
    if(expected!=0)check(reactionWitness,"F7 fixed physical boundary reaction witness missing");
    std::vector<long double> div(l.pressureCount());long double work=0,scale=0;
    for(int u=0;u<l.velocityCount();++u){const double velocity=std::sin(.13*(u+1));const long double term=velocity*rhs[u];work+=term;scale+=std::abs(term);
        if(ref.minus[u]>=0)div[ref.minus[u]]+=ref.area[u]*velocity;if(ref.plus[u]>=0)div[ref.plus[u]]-=ref.area[u]*velocity;}
    for(int p=0;p<l.pressureCount();++p){const long double term=potential[p]*div[p];work+=term;scale+=std::abs(term);}
    check(relative(static_cast<double>(work),static_cast<double>(scale))<=1e-12,"F7 original C=-D_I transpose power");
}
void scene(bool rz,bool column)
{
    const int n=rz&&!column?64:32;
    auto mesh=heightMesh(n);
    if(rz&&!column){DictionaryMesh::RefineBox box;box.level=1;box.min={0,0};box.max={.375,.125};mesh->applyRefineBoxes({box});}
    const auto metrics=metricsFor(*mesh,rz);const auto& layout=metrics.layout();const auto stocks=inputStocks(*mesh,rz,column);
    const Clock clock{41,0,.0001,Clock::Position::Midpoint};
    const auto pub=Vof::prepareLocalPlicInventory(*mesh,layout.domain(),clock,stocks,{});
    const auto saved=isolatedFragmentSnapshot(pub);const auto phase=HydroPhaseGeometry::build(*mesh,metrics,pub);
    const auto pressure=HydroPressureLoad::build(metrics.divergence());const auto rows=southRows(*mesh,pub);
    const auto boundary=Boundary::prepare(*mesh,pub,rows);boundary.validateFor(*mesh,pub);
    check(boundary.source().identity()==pub.identity()&&f7ClockEqual(boundary.source().timeLayer(),clock)
        &&boundary.coordinate()==0&&boundary.nativeRows().size()==rows.size(),"F7 actual immutable source/clock/plane/coverage");
    std::vector<double> lengths;std::vector<std::array<double,2>> spans,merged;
    for(int row:boundary.nativeRows()){const auto face=pub.domain().face(*mesh,row);
        const auto& edge=mesh->hEdges(face.edgeLabel);const double a=mesh->verts(edge.edgeToVert[0]).coord().x,b=mesh->verts(edge.edgeToVert[1]).coord().x;
        spans.push_back({std::min(a,b),std::max(a,b)});lengths.push_back(std::abs(b-a));}
    std::sort(spans.begin(),spans.end());
    for(const auto& span:spans){
        if(merged.empty()||span[0]>merged.back()[1])merged.push_back(span);
        else{check(span[0]==merged.back()[1],"F7 actual native South spans overlap");merged.back()[1]=span[1];}}
    check(boundary.intervals()==merged,"F7 declaration intervals differ from independent exact union of actual split South spans");
    if(rz&&!column){check(!pub.domain().uniformLeaves()&&*std::max_element(lengths.begin(),lengths.end())==2*(*std::min_element(lengths.begin(),lengths.end())),"F7 real 2:1/split south boundary missing");++splitSouth;}
    const auto load=HydroBalancedCsfLoad::localPlicHeightFunctionCapillary(*mesh,pressure,phase,.7,boundary,heightConfig());
    check(load.southEven90Boundary()&&load.southEven90Boundary()->sameIdentity(boundary),"F7 CSF lost exact boundary certificate owner");
    load.validateVofFor(*mesh,pressure,pub.identity(),clock);
    const double expected=column?(rz?4:0):(rz?8:4);unsigned mixed=0,reflected=0,straight=0;double maximum=0;
    for(int p=0;p<layout.pressureCount();++p){const auto& cell=pub.plicCells()[layout.nativeCell(p)];
        if(cell.measured[0]>0&&cell.measured[1]>0){
            const double k=load.heightFunctionCellCurvatures()[p];check(std::isfinite(k),"F7 declared boundary returned unavailable physical mixed-cell curvature");++mixed;
            if(!column)maximum=std::max(maximum,std::abs(k-expected)/expected);
            const auto& jet=load.heightFunctionCellJets()[p];if(jet)reflected+=jet->southEven90Reflected;
        }else check(!std::isfinite(load.heightFunctionCellCurvatures()[p]),"F7 pure-cell physical kappa must stay unavailable");
    }
    if(column){
        for(int u=0;u<layout.velocityCount();++u){const auto face=layout.domain().face(*mesh,layout.nativeFace(u));const auto& report=load.supportReports()[u];
            if(face.center.y<.2&&load.colorGradient()[u]!=0&&face.negativeCell>=0&&face.positiveCell>=0){
                check(report.coincidentPhaseHeight&&report.curvatureDefined,"F7 actual purejump face lacks its same-plane height query");
                check(std::abs(load.curvatureAverages()[u]-expected)<=1e-12/.25
                    &&std::abs(load.momentumRhsLoad()[u]-.7*expected*load.colorGradient()[u])<=1e-12*std::abs(.7*expected*load.colorGradient()[u]),
                    "F7 vertical plane/cylinder original constant face C force");++straight;}}
        check(straight>0&&load.curvatureReport().pureJumpHeightFaces>0,"F7 aligned purejump witness missing");straightFaces+=straight;
    }else{
        check(mixed>0&&(!rz||reflected>0)&&load.curvatureReport().southEven90BoundaryCells>0
            &&load.curvatureReport().southEven90ReflectedSlices>0,"F7 actual source reflection missing (RZ additionally requires its published jets)");
        check(maximum<=(rz?.06:.15),"F7 predeclared half-geometry curvature budget exceeded");
        if(rz)maximumRZ=std::max(maximumRZ,maximum);else maximumXY=std::max(maximumXY,maximum);
    }
    algebra(*mesh,phase,pressure,load,expected);
    const auto foreign=Vof::prepareLocalPlicInventory(*mesh,layout.domain(),clock,stocks,{});
    reject([&]{boundary.validateFor(*mesh,foreign);},"F7 equal-valued foreign publication accepted boundary owner");
    reject([&]{HydroBalancedCsfLoad::localPlicHeightFunctionCapillary(*mesh,pressure,HydroPhaseGeometry::build(*mesh,metrics,foreign),.7,boundary,heightConfig());},"F7 load consumed foreign boundary owner");
    auto wrongClock=clock;++wrongClock.step;
    reject([&]{load.validateVofFor(*mesh,pressure,pub.identity(),wrongClock);},"F7 full-clock mismatch accepted load");
    auto wrongSide=rows;
    for(int row=0;row<pub.domain().nativeFaceCount();++row){const auto f=pub.domain().face(*mesh,row);
        if(f.horizontal&&f.positiveCell<0&&f.area>0){wrongSide.front()=row;break;}}
    reject([&]{Boundary::prepare(*mesh,pub,wrongSide);},"F7 non-south actual native row accepted");
    if(!column){
        reject([&]{HydroBalancedCsfLoad::localPlicHeightFunctionCapillary(*mesh,pressure,phase,.7,heightConfig());},"F7 undeclared physical south wall received implicit contact law");
        std::vector<int> incomplete;
        for(int row:rows){const auto f=pub.domain().face(*mesh,row);const double contact=rz?.25:.75;
            if(std::abs(f.center.x-contact)>.1)incomplete.push_back(row);}
        reject([&]{const auto partial=Boundary::prepare(*mesh,pub,incomplete);
            HydroBalancedCsfLoad::localPlicHeightFunctionCapillary(*mesh,pressure,phase,.7,partial,heightConfig());},"F7 incomplete actual contact-chart coverage accepted");
    }
    for(std::size_t i=0;i<stocks.size();++i)check(pub.plicCells()[i].inventory==stocks[i],"F7 load/proof/negative changed either independent stock");
    check(isolatedFragmentSnapshot(pub)==saved,"F7 source alpha/normal/offset/clock changed");
    ++directScenes;std::cout<<"F7_SOUTH_EVEN90_GEOMETRY,"<<(rz?"RZ":"XY")<<','<<(column?"capped_column":"half_sphere")
        <<",leaves,"<<mesh->cellIndexNum()<<",mixed,"<<mixed<<",reflected_jets,"<<reflected<<",purejump_straight,"<<straight<<",maximum_relative,"<<maximum<<'\n';
}
void providerStep()
{
    auto mesh=heightMesh(32);const auto metrics=metricsFor(*mesh,true);const auto& layout=metrics.layout();
    const auto stocks=inputStocks(*mesh,true,false);
    const Clock initial{0,0,0,Clock::Position::Initial},target{1,0,.00001,Clock::Position::End};
    const auto pub=Vof::prepareLocalPlicInventory(*mesh,layout.domain(),initial,stocks,{});
    const auto phase=HydroPhaseGeometry::build(*mesh,metrics,pub);
    auto accepted=HydroTwoPhaseState::createVof(*mesh,phase,HydroNativeValues(layout,
        std::vector<double>(layout.velocityCount()),std::vector<double>(layout.pressureCount())));
    const auto identity=accepted;const auto source=isolatedFragmentSnapshot(pub);const auto oldU=accepted.values().normalVelocity(),oldP=accepted.values().pressure();
    std::vector<Boundary> views;HydroVofAcConfig config;config.liquidDensity=2;config.gasDensity=1;config.surfaceTension=.7;
    config.momentumDensityRepresentation=HydroMomentumMass::DensityRepresentation::CellAverageProjection;config.curvature=heightConfig();
    config.southEven90Boundary=[&](const QuadtreeMesh& grid,const Vof& actual){
        check(actual.timeLayer().position==Clock::Position::Midpoint&&actual.identity()!=pub.identity(),"F7 provider reused Initial publication");
        views.push_back(Boundary::prepare(grid,actual,southRows(grid,actual)));return views.back();};
    const auto step=HydroVofStep::prepare(*mesh,accepted,target,config);
    check(!views.empty()&&views.size()==step.report().stages.size()&&step.hasCapillaryLoad(),"F7 provider not consumed by every actual positive-sigma stage");
    const auto& load=step.capillaryLoad();const auto* view=load.southEven90Boundary();
    check(view&&view->sameIdentity(views.back())&&view->source().identity()==step.transport().midpoint().identity()
        &&f7ClockEqual(view->source().timeLayer(),step.transport().midpoint().timeLayer()),"F7 actual final Midpoint load/provider/transport owner mismatch");
    check(accepted.sameIdentity(identity)&&accepted.values().normalVelocity()==oldU&&accepted.values().pressure()==oldP
        &&isolatedFragmentSnapshot(pub)==source,"F7 real Step mutated accepted source before commit");
    auto bad=config;const auto stale=Boundary::prepare(*mesh,pub,southRows(*mesh,pub));
    bad.southEven90Boundary=[stale](const QuadtreeMesh&,const Vof&){return stale;};
    reject([&]{HydroVofStep::prepare(*mesh,accepted,target,bad);},"F7 stage accepted stale Initial provider certificate");
    check(accepted.sameIdentity(identity)&&isolatedFragmentSnapshot(pub)==source,"F7 failed stale provider changed source");
    step.commit(*mesh,accepted);check(accepted.timeLayer().step==1&&f7ClockEqual(accepted.timeLayer(),target),"F7 real provider Step did not commit actual End");
    for(const auto& report:step.report().stages)check(report.saddle.algebraicRelativeResidual<=1e-11
        &&report.saddle.freeMomentumDefect<=1e-11&&report.saddle.maxCellContinuityDefect<=1e-11,"F7 original actual Step Saddle gates");
    providerStages+=static_cast<unsigned>(views.size());
    std::cout<<"F7_SOUTH_EVEN90_ACTUAL_STEP,accepted,1,AC_stages,"<<views.size()<<",last_mid_geometry,"<<view->source().identity().value()<<'\n';
}
void run()
{
    for(bool rz:{false,true}){scene(rz,true);scene(rz,false);}providerStep();
    check(directScenes==4&&straightFaces>0&&splitSouth==1&&providerStages>0&&boundaryGuards>=20,"F7 required direct/provider/negative witnesses incomplete");
    std::cout<<"F7_SOUTH_EVEN90_PASS,direct,"<<directScenes<<",guards,"<<boundaryGuards<<",purejump_faces,"<<straightFaces
        <<",provider_stages,"<<providerStages<<",XY_half_max,"<<maximumXY<<",RZ_half_max,"<<maximumRZ
        <<",constant_C_free,"<<constantFree<<",fixed_reaction,"<<constantReaction<<",NO_GENERAL_CONTACT_ANGLE_OR_PINNED_LINE_CLAIM\n";
}
}

}

// BEGIN F25 ACTUAL VERTEX-STAR FRAGMENT CONTRACT
namespace F25VertexStarFragmentContract {
using Publication=LAMP::Interface::GeometricVofPublication;
using Clock=LAMP::TimeLayer;
using Point=TwoCellPoint;
unsigned positive=0,queries=0,negative=0,ownerGuards=0,stepStages=0,circleCells=0;
struct Input {
    std::vector<std::array<long double,2>> stocks;
    std::vector<int> cells;
    int emptyStock=-1,external=-1;
    Vector2D vertex;
};
bool equalClock(Clock a,Clock b){return a.step==b.step&&a.begin==b.begin&&a.end==b.end&&a.position==b.position;}
template<class F>void guard(F action,const char* why){isolatedFragmentReject(action,why);++ownerGuards;}
Input inputFor(const QuadtreeMesh& mesh,const HydroPhysicalDivergence& metrics,unsigned phase,int count,
    bool reflect=false,long double fraction=0,Vector2D vertex={.5,.5},int variant=0)
{
    check(phase<2&&(count==3||count==4),"F25 invalid declarative fixture");
    Input out;out.vertex=vertex;out.stocks.resize(mesh.cellIndexNum());const auto& layout=metrics.layout();
    for(int ci=0;ci<mesh.cellIndexNum();++ci){
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[ci]);
        const auto lo=mesh.verts(cell.cellToVert[0]).coord(),hi=mesh.verts(cell.cellToVert[2]).coord();
        const bool touchX=lo.x==vertex.x||hi.x==vertex.x,touchY=lo.y==vertex.y||hi.y==vertex.y;
        const int qx=lo.x>=vertex.x?1:0,qy=lo.y>=vertex.y?1:0;
        const int logical=2*(qx^(reflect?1:0))+qy;
        long double f=0;
        if(touchX&&touchY){
            if(variant==2)f=(logical==1||logical==2)?1.L/1024:0; // Diagonal-only point contact.
            else if(count==4)f=fraction>0?fraction:1.L/2048;
            else f=logical==1?1.L/1024:logical==3?1.L/2048:logical==2?1.L/4096:1e-50L;
            if(f>=1.L/4096)out.cells.push_back(ci);else if(f>0)out.emptyStock=ci;
        }
        if(variant==1&&hi.x==vertex.x-1./16&&lo.y==vertex.y){f=1;out.external=ci;}
        if(variant==3)f=1.L/2048; // Genuine unresolved measured-mixed plateau.
        const long double capacity=metrics.cellVolumes()[layout.pressureDof(ci)];
        const long double minor=capacity*f,major=capacity*(1-f);
        out.stocks[ci]=phase==0?std::array<long double,2>{{minor,major}}:std::array<long double,2>{{major,minor}};
    }
    std::sort(out.cells.begin(),out.cells.end());return out;
}
struct Shape {
    std::array<long double,4> bbox{{std::numeric_limits<long double>::infinity(),std::numeric_limits<long double>::infinity(),
        -std::numeric_limits<long double>::infinity(),-std::numeric_limits<long double>::infinity()}};
    long double maximumPair=0,midpointPair=0;
    std::vector<Point> vertices,midpoints;
};
Shape actualShape(const Publication& pub,const std::vector<int>& cells,unsigned phase)
{
    Shape out;
    for(int ci:cells){
        const auto& c=pub.plicCells()[ci];
        check(c.resolved&&c.measured[0]>0&&c.measured[1]>0,"F25 declared carrier must be genuinely resolved measured-mixed");
        const auto polygon=twoCellActualPolygon(pub,ci,phase);check(polygon.size()>=3,"F25 empty actual clipped polygon");
        for(auto p:polygon){out.vertices.push_back(p);out.bbox[0]=std::min(out.bbox[0],p[0]);out.bbox[1]=std::min(out.bbox[1],p[1]);
            out.bbox[2]=std::max(out.bbox[2],p[0]);out.bbox[3]=std::max(out.bbox[3],p[1]);}
        const auto segment=LAMP::Geometry::PlanarCellCut::interfaceSegment(c.rectangle,c.normal,c.offset);
        check(bool(segment),"F25 actual carrier lacks its real PLIC segment");
        const auto& b=c.rectangle;const long double cx=(static_cast<long double>(b.lower.x)+b.upper.x)/2,
            cy=(static_cast<long double>(b.lower.y)+b.upper.y)/2;
        out.midpoints.push_back({cx+((*segment)[0].x+(*segment)[1].x)/2,cy+((*segment)[0].y+(*segment)[1].y)/2});
    }
    for(auto a:out.vertices)for(auto b:out.vertices)out.maximumPair=std::max(out.maximumPair,std::hypot(a[0]-b[0],a[1]-b[1]));
    for(auto a:out.midpoints)for(auto b:out.midpoints)out.midpointPair=std::max(out.midpointPair,std::hypot(a[0]-b[0],a[1]-b[1]));
    return out;
}
long double l1(const Shape& s){return s.bbox[2]-s.bbox[0]+s.bbox[3]-s.bbox[1];}
long double commonTrace(const Publication& pub,int first,int second,unsigned phase)
{
    const auto& a=pub.plicCells()[first].rectangle;const auto& b=pub.plicCells()[second].rectangle;
    const bool vertical=(a.upper.x==b.lower.x||b.upper.x==a.lower.x)&&a.lower.y==b.lower.y&&a.upper.y==b.upper.y;
    const bool horizontal=(a.upper.y==b.lower.y||b.upper.y==a.lower.y)&&a.lower.x==b.lower.x&&a.upper.x==b.upper.x;
    if(!vertical&&!horizontal)return 0;
    const long double fixed=vertical?(a.upper.x==b.lower.x?a.upper.x:a.lower.x):(a.upper.y==b.lower.y?a.upper.y:a.lower.y);
    const auto interval=[&](int native){std::array<long double,2> span{{std::numeric_limits<long double>::infinity(),-std::numeric_limits<long double>::infinity()}};
        for(auto p:twoCellActualPolygon(pub,native,phase))if((vertical?p[0]:p[1])==fixed){const auto x=vertical?p[1]:p[0];span[0]=std::min(span[0],x);span[1]=std::max(span[1],x);}return span;};
    const auto x=interval(first),y=interval(second);
    return std::max(0.L,std::min(x[1],y[1])-std::max(x[0],y[0]));
}
void positiveTraceConnected(const Publication& pub,const std::vector<int>& cells,unsigned phase)
{
    std::vector<bool> reached(cells.size());reached[0]=true;unsigned edges=0;
    for(std::size_t i=0;i<cells.size();++i)for(std::size_t j=i+1;j<cells.size();++j)edges+=commonTrace(pub,cells[i],cells[j],phase)>0;
    for(std::size_t pass=0;pass<cells.size();++pass)for(std::size_t i=0;i<cells.size();++i)if(reached[i])
        for(std::size_t j=0;j<cells.size();++j)if(commonTrace(pub,cells[i],cells[j],phase)>0)reached[j]=true;
    check(edges>=cells.size()-1&&std::all_of(reached.begin(),reached.end(),[](bool value){return value;}),
        "F25 independent actual polygon traces do not connect all declared carriers by positive length");
}
void unchanged(const Publication& pub,const Input& input,const std::vector<std::array<long double,8>>& before)
{
    check(isolatedFragmentSnapshot(pub)==before,"F25 certificate/load/negative changed a retained plane, geometry or stock");
    for(std::size_t ci=0;ci<input.stocks.size();++ci)
        check(std::memcmp(pub.plicCells()[ci].inventory.data(),input.stocks[ci].data(),sizeof(long double)*2)==0,
            "F25 changed raw independent inventory bytes");
}
void positiveScene(int count,unsigned minor,bool reflect)
{
    constexpr double h=1./16,sigma=.7;const Vector2D origin=reflect?Vector2D{-.25,1.5}:Vector2D{0,0};
    auto mesh=heightMesh(16,origin);const auto metrics=metricsFor(*mesh,false);const auto& layout=metrics.layout();
    const auto input=inputFor(*mesh,metrics,minor,count,reflect,0,origin+Vector2D{.5,.5});
    check(input.cells.size()==static_cast<std::size_t>(count),"F25 missing declared three/four carriers");
    const Clock clock{25,.03,.0302,Clock::Position::Midpoint};
    const auto pub=Publication::prepareLocalPlicInventory(*mesh,layout.domain(),clock,input.stocks,{});
    const auto before=isolatedFragmentSnapshot(pub);const auto alpha=pub.alpha();const auto* storage=pub.plicCells().data();
    const auto shape=actualShape(pub,input.cells,minor);positiveTraceConnected(pub,input.cells,minor);
    check(shape.maximumPair>0&&l1(shape)<h&&shape.midpointPair<h,
        "F25 true polygon/segment positions must be below one actual Delta");
    if(count==3){check(input.emptyStock>=0,"F25 three-carrier source lacks retained tiny raw stock in fourth quadrant");
        const auto& c=pub.plicCells()[input.emptyStock];
        check(c.inventory[minor]>0&&c.inventory[1-minor]>0&&c.measured[minor]==0&&!c.resolved,
            "F25 fourth quadrant must be measured Empty without deleting its positive independent stock");}
    std::array<long double,2> canonicalLo{},canonicalHi{};
    for(int native:input.cells){
        check(!pub.isolatedLocalPlicFragment(*mesh,native)&&!pub.twoCellLocalPlicPositionDeficit(*mesh,native),
            "F25 three/four actual carriers must fail the old singleton/two-cell component routes");
        const auto proof=pub.vertexStarLocalPlicFragment(*mesh,native);
        check(bool(proof),"F25 closed tiny actual vertex star received no owner-issued qualification");
        proof->validateFor(*mesh,pub);
        check(proof->sourceIdentity()==pub.identity()&&proof->nativeCell()==native&&proof->nativeCells()==input.cells
            &&proof->phase()==minor&&equalClock(proof->sourceTime(),clock),"F25 certificate lost canonical actual component/source/fullclock");
        const auto vertex=proof->starVertex();
        check(vertex.x==input.vertex.x&&vertex.y==input.vertex.y&&proof->delta()==h,
            "F25 certificate changed actual common vertex or native Delta");
        const auto lo=proof->bboxLower(),hi=proof->bboxUpper();const auto& anchor=pub.plicCells()[input.cells.front()].rectangle;
        const long double cx=(static_cast<long double>(anchor.lower.x)+anchor.upper.x)/2,
            cy=(static_cast<long double>(anchor.lower.y)+anchor.upper.y)/2;
        const auto allowance=64*std::numeric_limits<long double>::epsilon()*h;
        check(lo[0]+cx<=shape.bbox[0]+allowance&&lo[1]+cy<=shape.bbox[1]+allowance
            &&hi[0]+cx>=shape.bbox[2]-allowance&&hi[1]+cy>=shape.bbox[3]-allowance,
            "F25 owner bounds do not enclose independent actual clipped polygons");
        check(proof->diameter()>0&&proof->diameter()<h&&proof->diameter()+allowance>=l1(shape)
            &&proof->diameter()+allowance>=shape.maximumPair,"F25 bound is downward or lacks strict native-resolution margin");
        check(proof->coveredNativeFaces()>=static_cast<std::size_t>(4*count-(count==4?4:2))&&proof->checkedCornerCells()>0,
            "F25 proof omitted actual native-face/perimeter/corner incidence");
        if(native==input.cells.front()){canonicalLo=lo;canonicalHi=hi;
            const auto other=Publication::prepareLocalPlicInventory(*mesh,layout.domain(),clock,input.stocks,{});
            guard([&]{proof->validateFor(*mesh,other);},"F25 certificate accepted another actual same-numeric owner");
            auto changed=clock;++changed.step;const auto stepOwner=Publication::prepareLocalPlicInventory(*mesh,layout.domain(),changed,input.stocks,{});
            guard([&]{proof->validateFor(*mesh,stepOwner);},"F25 certificate accepted same interval/time but another step");
            changed=clock;changed.end+=.0001;const auto intervalOwner=Publication::prepareLocalPlicInventory(*mesh,layout.domain(),changed,input.stocks,{});
            guard([&]{proof->validateFor(*mesh,intervalOwner);},"F25 certificate accepted another full interval");
            guard([&]{pub.validateFor(*mesh,pub.identity(),changed);},"F25 source accepted wrong fullclock");
            auto foreign=heightMesh(32,origin);check(foreign->cellIndexNum()!=mesh->cellIndexNum(),"F25 foreign-mesh guard needs genuinely different native topology");
            guard([&]{proof->validateFor(*foreign,pub);},"F25 certificate accepted another actual mesh");
        }else check(lo==canonicalLo&&hi==canonicalHi,"F25 carrier queries disagree on canonical whole-component bounds");
        ++queries;
    }
    const auto phase=HydroPhaseGeometry::build(*mesh,metrics,pub);const auto pressure=HydroPressureLoad::build(metrics.divergence());
    const auto load=HydroBalancedCsfLoad::localPlicHeightFunctionCapillary(*mesh,pressure,phase,sigma,heightConfig());
    load.validateVofFor(*mesh,pressure,pub.identity(),clock);
    check(load.curvatureReport().unresolvedVertexStarFragmentCells==static_cast<std::size_t>(count)
        &&load.curvatureReport().unresolvedIsolatedFragmentCells==static_cast<std::size_t>(count)
        &&load.curvatureReport().heightFunctionCells==0&&load.curvatureReport().heightFunctionFitCells==0,
        "F25 new subset counter must not be added to existing total or successful HF/fit counts");
    for(int native:input.cells){const int p=layout.pressureDof(native);
        check(load.heightFunctionCellCurvatureMethods()[p]==HydroCsfCellCurvatureMethod::UnderresolvedIsolatedFragmentZero
            &&std::isnan(load.heightFunctionCellCurvatures()[p])&&load.heightFunctionCellEffectiveCurvatures()[p]==0
            &&!load.heightFunctionCellJets()[p],"F25 effective closure falsely supplies physical kappa or a qualified jet");}
    std::vector<double> color(layout.pressureCount());
    for(int p=0;p<layout.pressureCount();++p){const auto v=phase.cellVolumes()[p];color[p]=v.liquid==0?0:v.gas==0?1:static_cast<double>(v.liquid/metrics.cellVolumes()[p]);}
    const Reference independent(*mesh,phase);const auto gradient=independent.gradient(color);unsigned active=0,pureBoundary=0;
    for(int u=0;u<layout.velocityCount();++u){close(load.colorGradient()[u],gradient[u],1e-12*std::max(1.,std::abs(gradient[u])),"F25 retained measured native color/incidence");
        check(std::isfinite(load.momentumRhsLoad()[u])&&load.momentumRhsLoad()[u]==0,"F25 effective zero created actual force");
        if(gradient[u]!=0){
            const auto face=layout.domain().face(*mesh,layout.nativeFace(u));const auto& support=load.supportReports()[u];
            const bool adjacent=std::find(input.cells.begin(),input.cells.end(),face.negativeCell)!=input.cells.end()
                ||std::find(input.cells.begin(),input.cells.end(),face.positiveCell)!=input.cells.end();
            if(adjacent){check(support.underresolvedFragmentClosure&&load.heightFunctionFaceCurvatures()[u]==0,
                "F25 active carrier-adjacent force row lost declared closure provenance");++active;}
            else{check((face.negativeCell<0||face.positiveCell<0)&&support.status==HydroCsfSupportStatus::PureBoundaryReaction
                &&!support.curvatureDefined&&!support.underresolvedFragmentClosure,
                "F25 noninterface active row must retain its pure physical boundary reaction");++pureBoundary;}
        }}
    check(active>0,"F25 consumer needs a nontrivial actual carrier-adjacent native gradient");
    check(minor==0?pureBoundary==0:pureBoundary>0,"F25 complement must retain its actual pure-liquid boundary reaction witness");
    check(pub.plicCells().data()==storage&&pub.alpha()==alpha,"F25 certificate/load changed actual storage or alpha");
    unchanged(pub,input,before);++positive;
    std::cout<<"F25_VERTEX_STAR_POSITIVE,carriers,"<<count<<",minor_phase,"<<minor<<",reflected_shifted,"<<reflect
        <<",actual_polygon_max_pair_over_delta,"<<shape.maximumPair/h<<",actual_midpoint_max_pair_over_delta,"<<shape.midpointPair/h
        <<",physical_kappa,NA,effective_kappa,0,total_closure_cells,"<<count<<",subset_vertex_star_cells,"<<count<<",pure_boundary_reaction_rows,"<<pureBoundary<<'\n';
}
void negatives()
{
    constexpr double h=1./16;const Clock clock{25,.03,.0302,Clock::Position::Midpoint};
    // The near-boundary family straddles the exact strict L1 criterion with a
    // margin larger than root/enclosure arithmetic; no epsilon softens <Delta.
    for(bool above:{false,true}){
        auto mesh=heightMesh(16);const auto metrics=metricsFor(*mesh,false);
        const auto input=inputFor(*mesh,metrics,0,4,false,(1.L/32)*(1+(above?1.L:-1.L)/1048576));
        const auto pub=Publication::prepareLocalPlicInventory(*mesh,metrics.layout().domain(),clock,input.stocks,{});
        const auto before=isolatedFragmentSnapshot(pub);const auto shape=actualShape(pub,input.cells,0);
        check(above?l1(shape)>h:l1(shape)<h,"F25 independent public geometry did not straddle the strict boundary");
        for(int native:input.cells){const auto proof=pub.vertexStarLocalPlicFragment(*mesh,native);
            if(above)check(!proof,"F25 >=Delta whole bound received a strict certificate");
            else check(proof&&proof->diameter()<h,"F25 just-below-Delta genuine closed star was rejected");}
        unchanged(pub,input,before);++negative;
        std::cout<<"F25_VERTEX_STAR_BOUNDARY,above,"<<above<<",actual_bbox_L1_over_delta,"<<l1(shape)/h<<'\n';
    }
    for(int variant:{1,2,3}){
        auto mesh=heightMesh(16);const auto metrics=metricsFor(*mesh,false);
        const auto input=inputFor(*mesh,metrics,0,3,false,0,{.5,.5},variant);
        const auto pub=Publication::prepareLocalPlicInventory(*mesh,metrics.layout().domain(),clock,input.stocks,{});
        const auto before=isolatedFragmentSnapshot(pub);
        if(variant==1){check(input.external>=0&&pub.plicCells()[input.external].measured[1]==0,"F25 external continuation lacks a genuine full same-phase leaf");}
        if(variant==2){check(input.cells.size()==2,"F25 point-contact negative needs only diagonal carriers");
            for(int native:input.cells)check(twoCellContainsCorner(pub,native,0,{.5,.5}),"F25 diagonal actual polygons do not touch their common point");}
        if(variant==3){const auto& c=pub.plicCells()[input.cells.front()];
            check(!c.resolved&&c.measured[0]>0&&c.measured[1]>0,"F25 unknown negative must be actual measured-mixed, not declared geometry Empty");}
        if(variant==1){
            int source=-1;for(int native:input.cells){const auto& b=pub.plicCells()[native].rectangle;if(b.upper.x==.5&&b.lower.y==.5)source=native;}
            check(source>=0&&commonTrace(pub,source,input.external,0)>0,"F25 external negative lacks true positive continuation outside the proposed tiny star");
            check(!pub.vertexStarLocalPlicFragment(*mesh,source),"F25 external continuation was certified");
        }else for(int native:input.cells)check(!pub.vertexStarLocalPlicFragment(*mesh,native),"F25 point-only/unknown geometry was certified");
        unchanged(pub,input,before);++negative;
    }
    for(int variant:{0,1}){auto mesh=heightMesh(16);const bool rz=variant==1;const auto metrics=metricsFor(*mesh,rz);
        const Vector2D vertex=variant==0?Vector2D{0,.5}:Vector2D{.5,.5};
        const auto input=inputFor(*mesh,metrics,0,4,false,0,vertex);
        const auto pub=Publication::prepareLocalPlicInventory(*mesh,metrics.layout().domain(),clock,input.stocks,{});
        check(!input.cells.empty(),"F25 wall/RZ negative has no actual carrier");
        for(int native:input.cells)check(!pub.vertexStarLocalPlicFragment(*mesh,native),"F25 physical wall/RZ received an interior XY certificate");++negative;}
    {
        auto mesh=heightMesh(16);DictionaryMesh::RefineBox b;b.level=1;b.min={.5,.5};b.max={.5625,.5625};mesh->applyRefineBoxes({b});
        const auto metrics=metricsFor(*mesh,false);const auto input=inputFor(*mesh,metrics,0,4);
        unsigned split=0;for(int row=0;row<metrics.layout().domain().nativeFaceCount();++row){const auto f=metrics.layout().domain().face(*mesh,row);
            if(f.negativeCell>=0&&f.positiveCell>=0&&mesh->cells(mesh->cellIndexToLabel[f.negativeCell]).level()!=mesh->cells(mesh->cellIndexToLabel[f.positiveCell]).level())++split;}
        check(split>0,"F25 mixed-level negative needs actual native split faces");
        const auto pub=Publication::prepareLocalPlicInventory(*mesh,metrics.layout().domain(),clock,input.stocks,{});
        for(int native:input.cells)check(!pub.vertexStarLocalPlicFragment(*mesh,native),"F25 mixed-level vertex star was certified");++negative;
    }
    {
        auto mesh=heightMesh(16);const auto metrics=metricsFor(*mesh,false);const auto input=inputFor(*mesh,metrics,0,3);
        FaceSegmentGraph graph;graph.build(*mesh,[](const QuadtreeMesh&,const Cell& c){return c.centerCoord().x>.5;});
        const auto open=NativeLeafDomain::build(*mesh,graph,false);
        check(open.cellIndexByRow().size()<input.stocks.size(),"F25 open selected-domain fixture did not remove native coverage");
        guard([&]{Publication::prepareLocalPlicInventory(*mesh,open,clock,input.stocks,{});},"F25 incomplete/open source bypassed complete-publication guard");++negative;
        const auto fan=Publication::prepare(*mesh,metrics.layout().domain(),clock,std::vector<double>(mesh->cellIndexNum(),0));
        check(!fan.vertexStarLocalPlicFragment(*mesh,input.cells.front()),"F25 non-LocalPlic owner received a vertex-star certificate");++negative;
    }
}
void resolvedCircle()
{
    auto mesh=heightMesh(32);const auto metrics=metricsFor(*mesh,false);const auto& layout=metrics.layout();
    const Vector2D center{.487,.521};constexpr double radius=.25;std::vector<double> alpha(mesh->cellIndexNum());
    for(int ci=0;ci<mesh->cellIndexNum();++ci){const auto& c=mesh->cells(mesh->cellIndexToLabel[ci]);
        const auto lo=mesh->verts(c.cellToVert[0]).coord(),hi=mesh->verts(c.cellToVert[2]).coord();
        alpha[ci]=static_cast<double>(circleCell(lo,hi,center,radius,false)/((hi.x-lo.x)*(hi.y-lo.y)));}
    const Clock clock{25,.03,.0302,Clock::Position::Midpoint};
    const auto pub=Publication::prepareLocalPlic(*mesh,layout.domain(),clock,alpha);const auto saved=isolatedFragmentSnapshot(pub);
    const auto phase=HydroPhaseGeometry::build(*mesh,metrics,pub);const auto pressure=HydroPressureLoad::build(metrics.divergence());
    const auto load=HydroBalancedCsfLoad::localPlicHeightFunctionCapillary(*mesh,pressure,phase,.7,heightConfig());
    for(std::size_t p=0;p<load.heightFunctionCellCurvatures().size();++p)if(std::isfinite(load.heightFunctionCellCurvatures()[p])){
        check(load.heightFunctionCellCurvatures()[p]!=0&&load.heightFunctionCellCurvatureMethods()[p]!=HydroCsfCellCurvatureMethod::UnderresolvedIsolatedFragmentZero,
            "F25 resolved main circle was globally zeroed");++circleCells;}
    check(circleCells>0&&load.curvatureReport().unresolvedVertexStarFragmentCells==0
        &&load.curvatureReport().unresolvedIsolatedFragmentCells==0,"F25 resolved circle acquired an underresolved component closure");
    check(isolatedFragmentSnapshot(pub)==saved,"F25 circle regression changed actual stocks/planes");
}
void actualStep()
{
    auto mesh=heightMesh(16);const auto metrics=metricsFor(*mesh,false);const auto& layout=metrics.layout();
    const auto input=inputFor(*mesh,metrics,0,3);const Clock initial{0,0,0,Clock::Position::Initial},end{1,0,1e-6,Clock::Position::End};
    const auto pub=Publication::prepareLocalPlicInventory(*mesh,layout.domain(),initial,input.stocks,{});
    const auto phase=HydroPhaseGeometry::build(*mesh,metrics,pub);
    auto accepted=HydroTwoPhaseState::createVof(*mesh,phase,HydroNativeValues(layout,
        std::vector<double>(layout.velocityCount()),std::vector<double>(layout.pressureCount())));
    const auto source=accepted;const auto saved=isolatedFragmentSnapshot(pub);HydroVofAcConfig config;
    config.liquidDensity=2;config.gasDensity=1;config.surfaceTension=.7;config.curvature=heightConfig();
    config.momentumDensityRepresentation=HydroMomentumMass::DensityRepresentation::CellAverageProjection;
    const auto step=HydroVofStep::prepare(*mesh,accepted,end,config);check(step.hasCapillaryLoad(),"F25 actual positive-sigma Step omitted capillary preparation");
    const auto& load=step.capillaryLoad();const auto& mid=step.transport().midpoint();
    check(load.vofIdentity()==mid.identity()&&equalClock(load.timeLayer(),mid.timeLayer())
        &&load.curvatureReport().unresolvedVertexStarFragmentCells==3,"F25 actual Step did not consume genuine transported Midpoint star");
    for(const auto& report:step.report().stages){check(report.saddle.algebraicRelativeResidual<=1e-11
        &&report.saddle.freeMomentumDefect<=1e-11&&report.saddle.maxCellContinuityDefect<=1e-11
        &&report.momentumPowerRelativeDefect<=1e-11&&report.energy.closureRelativeResidual<=1e-11L,"F25 actual Step changed original Saddle/power/energy gates");++stepStages;}
    check(accepted.sameIdentity(source)&&isolatedFragmentSnapshot(pub)==saved,"F25 prepare changed the accepted source");
    step.commit(*mesh,accepted);check(equalClock(accepted.timeLayer(),end),"F25 actual Step did not commit the declared End");
}
void run()
{
    for(int count:{3,4})for(unsigned minor:{0u,1u})for(bool reflect:{false,true})positiveScene(count,minor,reflect);
    negatives();resolvedCircle();actualStep();
    check(positive==8&&queries==28&&negative==10&&ownerGuards==41&&circleCells>0&&stepStages>0,"F25 required finite qualification/consumer/Step evidence incomplete");
    std::cout<<"F25_VERTEX_STAR_FRAGMENT_PASS,positive,"<<positive<<",carrier_queries,"<<queries<<",boundary_and_negative_scenes,"<<negative
        <<",owner_clock_mesh_guards,"<<ownerGuards<<",resolved_circle_nonzero_cells,"<<circleCells<<",actual_AC_stages,"<<stepStages
        <<",physical_kappa,NA,effective_kappa,0,subset_counter_not_additive,stocks_unchanged,NO_1811_REPLAY_OR_PHYSICAL_BREAKUP_CLAIM\n";
}
} // namespace F25VertexStarFragmentContract
// END F25 ACTUAL VERTEX-STAR FRAGMENT CONTRACT

namespace PinnedCornerHfContract
{
namespace I=LAMP::Interface;
using Vof=I::GeometricVofPublication;using Clock=LAMP::TimeLayer;
enum class Shape { Hemicap,Cylinder,Plane };
unsigned scenes=0,negative=0,actualStages=0,remaps=0;double maxForce=0;
std::unique_ptr<QuadtreeMesh> nozzle(int scale,bool xy=false,bool refine=false)
{
    DictionaryMesh cfg;cfg.Nx=(xy?64:40)*scale;cfg.Ny=128*scale;cfg.minCoord={0,0};cfg.maxCoord=xy?Vector2D{2,4}:Vector2D{5,16};
    cfg.maxCellNum=500000;cfg.maxInitialIterations=0;auto mesh=std::make_unique<QuadtreeMesh>(cfg);
    if(refine){DictionaryMesh::RefineBox b;b.level=1;b.min=cfg.minCoord;b.max=cfg.maxCoord;mesh->applyRefineBoxes({b});}return mesh;
}
NativeLeafDomain selected(const QuadtreeMesh& mesh,bool rz,bool hole=false)
{
    FaceSegmentGraph graph;graph.build(mesh,[hole](const QuadtreeMesh&,const Cell& c){const auto p=c.centerCoord();
        return !(p.x>1&&p.y<3)&&!(hole&&p.x>.5&&p.x<.75&&p.y>3.125&&p.y<3.375);});
    return NativeLeafDomain::build(mesh,graph,rz);
}
HydroPhysicalDivergence metrics(const QuadtreeMesh& mesh,const NativeLeafDomain& d)
{return HydroPhysicalDivergence::build(mesh,HydroIntegralDivergence::build(mesh,HydroDofLayout::build(mesh,d)));}
std::pair<int,int> legs(const QuadtreeMesh& mesh,const NativeLeafDomain& d)
{
    int wet=-1,dry=-1;
    for(const auto f:d.faces(mesh)){
        if((f.negativeCell>=0)==(f.positiveCell>=0))continue;
        const auto& e=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
        const auto a=mesh.verts(e.edgeToVert[0]).coord(),b=mesh.verts(e.edgeToVert[1]).coord();
        if(!((a.x==1&&a.y==3)||(b.x==1&&b.y==3)))continue;
        if(f.horizontal&&f.center.x>1&&f.center.y==3)dry=f.nativeRow;
        if(!f.horizontal&&f.center.x==1&&f.center.y<3)wet=f.nativeRow;
    }
    check(wet>=0&&dry>=0,"HFP missing true selected nozzle wall legs");return {wet,dry};
}
std::vector<RzStocks> stocks(const QuadtreeMesh& mesh,const NativeLeafDomain& d,Shape shape)
{
    const bool rz=d.axisymmetric();std::vector<RzStocks> values(mesh.cellIndexNum());
    for(int ci:d.cellIndexByRow()){
        const auto& c=mesh.cells(mesh.cellIndexToLabel[ci]);const auto lo=mesh.verts(c.cellToVert[0]).coord(),hi=mesh.verts(c.cellToVert[2]).coord();
        const long double capacity=LAMP::Geometry::PlanarCellCut::capacity({lo,hi},rz);
        if(rz){
            const double base=shape==Shape::Cylinder?4:3;
            if(hi.y<=base){values[ci]=hi.x<=1?RzStocks{capacity,0}:RzStocks{0,capacity};continue;}
            if(shape==Shape::Plane){values[ci]={0,capacity};continue;}
            values[ci]=rzSphereColumn(lo,hi,base,1);
        }else{
            const double base=shape==Shape::Cylinder?3.5:shape==Shape::Plane?2.75:3;
            long double liquid=0;
            if(lo.y<base){const double top=std::min(base,hi.y),left=std::max(.5,lo.x),right=std::min(1.,hi.x);
                if(top>lo.y&&right>left)liquid=static_cast<long double>(top-lo.y)*(right-left);}
            if(hi.y>base){const Vector2D lower{lo.x,std::max(lo.y,base)};
                if(shape==Shape::Plane){
                    if(lo.x<.75)liquid+=circleCell(lower,{std::min(hi.x,.75),hi.y},{.75,base},.25,false);
                    const double left=std::max(.75,lo.x),right=std::min(1.,hi.x),top=std::min(3.,hi.y);
                    if(right>left&&top>lower.y)liquid+=static_cast<long double>(right-left)*(top-lower.y);
                }else liquid+=circleCell(lower,hi,{.75,base},.25,false);
            }
            check(liquid>=0&&liquid<=capacity,"HFP independent XY source integration outside capacity");values[ci]={liquid,capacity-liquid};
        }
    }
    return values;
}
Vof pinned(const QuadtreeMesh& mesh,const NativeLeafDomain& d,const std::vector<RzStocks>& input,Clock clock)
{
    const auto raw=Vof::prepareLocalPlicInventory(mesh,d,clock,input,{});const auto rows=legs(mesh,d);
    return Vof::withPinnedCorner(mesh,raw,{1,3},rows.first,rows.second);
}
HydroBalancedCsfLoad loadFor(const QuadtreeMesh& mesh,const Vof& pub)
{
    const auto met=metrics(mesh,pub.domain());const auto phase=HydroPhaseGeometry::build(mesh,met,pub);
    const auto pressure=HydroPressureLoad::build(met.divergence());
    if(pub.domain().axisymmetric())return HydroBalancedCsfLoad::localPlicHeightFunctionCapillary(mesh,pressure,phase,.01,heightConfig());
    std::vector<int> southRows;
    for(const auto f:pub.domain().faces(mesh))if(f.horizontal&&f.negativeCell<0&&f.positiveCell>=0&&
        f.outwardSign==-1&&f.area>0&&f.center.y==mesh.minCoord.y)southRows.push_back(f.nativeRow);
    const auto south=HydroVofSouthEven90Boundary::prepare(mesh,pub,southRows);
    return HydroBalancedCsfLoad::localPlicHeightFunctionCapillary(mesh,pressure,phase,.01,south,heightConfig());
}
void unchanged(const Vof& pub,const std::vector<RzStocks>& input)
{for(std::size_t ci=0;ci<input.size();++ci)check(pub.plicCells()[ci].inventory==input[ci],"HFP read-only force changed independent stocks");}
void algebra(const QuadtreeMesh& mesh,const HydroBalancedCsfLoad& load)
{
    const auto& phase=load.phaseGeometry();const auto& l=phase.layout();const auto& local=phase.vofPublication().plicCells();
    std::vector<double> color(l.pressureCount());
    for(int p=0;p<l.pressureCount();++p){const auto& c=local[l.nativeCell(p)];color[p]=c.measured[0]==0?0:c.measured[1]==0?1:static_cast<double>(c.measured[0]/phase.metrics().cellVolumes()[p]);}
    const auto expected=load.pressureLoad().apply(mesh,HydroNativeValues(l,std::vector<double>(l.velocityCount()),color));
    double error=0,scale=0;
    for(int u=0;u<l.velocityCount();++u){check(expected[u]==load.colorGradient()[u],"HFP color gradient lost its original pressure owner");
        const auto& report=load.supportReports()[u];const double force=report.curvatureDefined?.01*load.curvatureAverages()[u]*expected[u]:0;
        error=std::hypot(error,load.momentumRhsLoad()[u]-force);scale=std::hypot(scale,force);}
    maxForce=std::max(maxForce,relative(error,scale));check(relative(error,scale)<=1e-12,"HFP changed original same-face C force algebra");
}
template<class F>void reject(F action,const char* label)
{bool caught=false;try{action();}catch(const std::exception&){caught=true;}check(caught,label);++negative;}
void arithmetic()
{
    std::array<long double,3> error{};int n=0;
    for(long double h:{.125L,.0625L,.03125L}){
        std::array<long double,3> d{};
        for(int i=0;i<3;++i){const long double a=i*h,b=(i+1)*h;
            const auto primitive=[](long double x){return .2L*x*x/2-.4L*x*x*x/3+.3L*x*x*x*x/4+.1L*x*x*x*x*x/5;};
            d[i]=(primitive(b)-primitive(a))/h;}
        const long double second=(-19*d[0]+8*d[1]-d[2])/(3*h*h),x=h/2,actual=-.8L+1.8L*x+1.2L*x*x;
        error[n++]=std::abs(second-actual);
    }
    check(error[0]>0&&error[1]/error[0]>.24L&&error[1]/error[0]<.26L&&error[2]/error[1]>.24L&&error[2]/error[1]<.26L,
        "HFP exact analytic cellmean cubic second derivative is not second order");
    std::cout<<"HFP_ANALYTIC_CELLMEAN_P3,second_errors,"<<error[0]<<','<<error[1]<<','<<error[2]<<",NO_ACTUAL_PLIC_CLAIM\n";
}
double curve(int scale,bool rz)
{
    auto mesh=nozzle(scale,!rz);const auto d=selected(*mesh,rz);const auto input=stocks(*mesh,d,Shape::Hemicap);
    const Clock clock{3,0,.0001,Clock::Position::End};const auto pub=pinned(*mesh,d,input,clock);const auto saved=isolatedFragmentSnapshot(pub);
    const auto load=loadFor(*mesh,pub);algebra(*mesh,load);unchanged(pub,input);
    check(load.curvatureReport().finiteInventoryHeightColumns==0&&load.curvatureReport().finiteInventoryHeightQueries==0&&
        load.curvatureReport().finiteInventoryHeightCertificates.empty(),"HFP already-qualified source changed its original height path");
    check(isolatedFragmentSnapshot(pub)==saved,"HFP source geometry/clock mutated by force");
    check(load.curvatureReport().pinnedCornerHeightCells>0&&load.curvatureReport().pinnedCornerHeightColumns>=3,
        "HFP curved fixture never consumed the sealed pin P3 chart");
    double error=0;unsigned constrained=0;
    for(const auto& jet:load.heightFunctionCellJets())if(jet&&jet->pinnedCornerConstrained){
        ++constrained;const auto parts=rz?HydroBalancedCsfLoad::axisymmetricHeightCurvatures(*jet):std::array<long double,2>{
            -jet->second/std::pow(std::hypot(1.L,jet->first),3),0};
        error=std::max(error,static_cast<double>(std::abs(parts[0]+parts[1]-(rz?2:4))));
        check(jet->pinnedCornerBin>=0&&jet->pinnedCornerBin<=1,"HFP pin curvature escaped its first/second bin");
    }
    check(constrained>0&&std::isfinite(error),"HFP actual PLIC pin curvature unavailable");
    const auto foreign=pinned(*mesh,d,input,clock);
    reject([&]{load.validateVofFor(*mesh,load.pressureLoad(),foreign.identity(),clock);},"HFP stale source identity accepted force");
    auto wrong=clock;++wrong.step;reject([&]{load.validateVofFor(*mesh,load.pressureLoad(),pub.identity(),wrong);},"HFP stale full clock accepted force");
    reject([&]{HydroPhaseGeometry::build(*mesh,metricsFor(*mesh,rz),pub);},"HFP complete pressure map accepted selected publication");
    ++scenes;std::cout<<"HFP_ACTUAL_PLIC_CURVE,"<<(rz?"RZ_original_cap":"XY_strip_cap")<<",scale,"<<scale<<",pin_jets,"<<constrained<<",kappa_error,"<<error<<'\n';return error;
}
void pure(bool rz,Shape shape)
{
    auto mesh=nozzle(1,!rz);const auto d=selected(*mesh,rz);const auto input=stocks(*mesh,d,shape);
    const auto pub=pinned(*mesh,d,input,{4,0,.0001,Clock::Position::End});
    check(pub.pinnedCorner()->kind()==Vof::PinCertificate::Kind::NativeCoincident,"HFP pure pin manufactured a mixed facet");
    const auto load=loadFor(*mesh,pub);algebra(*mesh,load);unchanged(pub,input);unsigned witness=0;
    const auto& l=load.phaseGeometry().layout();
    for(int p=0;p<l.pressureCount();++p){const auto& cell=pub.plicCells()[l.nativeCell(p)];
        if(cell.inventory[0]==0||cell.inventory[1]==0)check(!std::isfinite(load.heightFunctionCellCurvatures()[p]),"HFP pure pressure cell acquired mixed curvature");}
    for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(*mesh,l.nativeFace(u));const auto& report=load.supportReports()[u];
        const bool nearby=shape==Shape::Cylinder?!f.horizontal&&f.center.x==1&&f.center.y>3&&f.center.y<3+(rz?.125:.03125):
            f.horizontal&&f.center.y==3&&f.center.x<1&&f.center.x>1-(rz?.125:.03125);
        if(nearby&&f.negativeCell>=0&&f.positiveCell>=0&&load.colorGradient()[u]!=0){
            const double expected=rz&&shape==Shape::Cylinder?1:0;
            check(report.coincidentPhaseHeight&&report.curvatureDefined&&std::abs(load.curvatureAverages()[u]-expected)<=1e-11,
                "HFP actual coincident cylinder/plane curvature failed");++witness;}
    }
    check(witness>0,"HFP pure certified face witness missing");++scenes;
    std::cout<<"HFP_NATIVE_COINCIDENT,"<<(rz?"RZ":"XY")<<','<<(shape==Shape::Cylinder?"vertical_cylinder":"horizontal_plane")<<",faces,"<<witness<<'\n';
}
HydroVofPrescribedFlow flow(const QuadtreeMesh& mesh,const HydroDofLayout& l)
{
    using Mode=I::GeometricVofAdvectorSlab::Boundary::Mode;HydroVofPrescribedFlow result;
    result.pressureDrivenReservoirs=true;result.pressureOutletBackflowModel=HydroVofPressureOutletMomentumModel::InteriorContinuationZeroShear;
    result.normalVelocity.assign(l.velocityCount(),0);result.rulesByNativeRow.assign(l.domain().nativeFaceCount(),Mode::Closed);
    for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(mesh,l.nativeFace(u));
        if(f.outwardSign==0)result.rulesByNativeRow[f.nativeRow]=Mode::Interior;
        else if(f.horizontal&&(f.center.y==0||f.center.y==16)){
            const bool supply=f.center.y==0;result.rulesByNativeRow[f.nativeRow]=supply?Mode::PureLiquidReservoir:Mode::PureGasReservoir;
            result.pressureOutlets.push_back({u,supply?.125:0});
        }else if((f.horizontal&&f.center.y==3&&f.center.x>1)||(!f.horizontal&&f.center.x==1&&f.center.y<3))
            result.stationaryNoSlipVelocityDofs.push_back(u);
    }
    check(!result.stationaryNoSlipVelocityDofs.empty(),"HFP actual no-slip wall mask missing");return result;
}
void stepAndAmr()
{
    auto mesh=nozzle(1);const auto d=selected(*mesh,true);const auto input=stocks(*mesh,d,Shape::Hemicap);
    const auto pub=pinned(*mesh,d,input,{0,0,0,Clock::Position::Initial});const auto met=metrics(*mesh,d);const auto& l=met.layout();
    const auto phase=HydroPhaseGeometry::build(*mesh,met,pub);
    auto accepted=HydroTwoPhaseState::createVofPrescribedFlow(*mesh,phase,
        HydroNativeValues(l,std::vector<double>(l.velocityCount()),std::vector<double>(l.pressureCount())),flow(*mesh,l));
    const auto before=accepted;const auto saved=isolatedFragmentSnapshot(pub);HydroVofAcConfig cfg;
    cfg.liquidDensity=2;cfg.gasDensity=1;cfg.liquidViscosity=.01;cfg.gasViscosity=.005;cfg.surfaceTension=.01;cfg.curvature=heightConfig();
    const Clock end{1,0,1./65536,Clock::Position::End};
    const auto sourcePin=*pub.pinnedCorner();
    for(const auto wall:std::array<std::pair<const char*,int>,2>{{{"wet",sourcePin.wetNativeRow()},{"dry",sourcePin.dryNativeRow()}}}){
        const auto f=d.face(*mesh,wall.second);const int ci=f.negativeCell>=0?f.negativeCell:f.positiveCell;
        std::cerr<<"HFP_OPERATION_SOURCE,"<<wall.first<<",cell,"<<ci<<",L,"<<pub.plicCells()[ci].inventory[0]
            <<",G,"<<pub.plicCells()[ci].inventory[1]<<",pin_source,"<<pub.identity().value()<<'\n';}
    std::cerr<<"HFP_OPERATION_BEGIN,preparePrescribedFlow,target_step,"<<end.step<<",begin,"<<end.begin<<",end,"<<end.end<<'\n'<<std::flush;
    const auto step=[&]{try{return HydroVofStep::preparePrescribedFlow(*mesh,accepted,end,cfg);}
        catch(const std::exception& error){std::cerr<<"HFP_OPERATION_FAIL,preparePrescribedFlow,"<<error.what()<<'\n'<<std::flush;throw;}}();
    check(step.hasCapillaryLoad()&&!step.report().stages.empty(),"HFP actual positive-sigma Hydro stage missing");
    const auto& load=step.capillaryLoad();const auto& mid=step.transport().midpoint();
    check(load.phaseGeometry().vofPublication().identity()==mid.identity()&&mid.pinnedCorner()&&
        load.curvatureReport().pinnedCornerHeightCells>0&&load.phaseGeometry().timeLayer().position==Clock::Position::Midpoint,
        "HFP actual stage did not consume its actual pinned Midpoint owner");
    check(step.frozenMidpoint().viscosity().viscosity().strain().stationaryNoSlipVelocityDofs()==accepted.vofPrescribedFlow()->stationaryNoSlipVelocityDofs,
        "HFP actual positive-sigma stage lost no-slip strain metadata");
    for(const auto& r:step.report().stages){check(r.saddle.algebraicRelativeResidual<=1e-11&&r.saddle.freeMomentumDefect<=1e-11&&
        r.saddle.maxCellContinuityDefect<=1e-11&&r.momentumPowerRelativeDefect<=1e-11&&r.energy.closureRelativeResidual<=1e-11L,
        "HFP changed original physical stage gates");++actualStages;}
    check(accepted.sameIdentity(before)&&isolatedFragmentSnapshot(pub)==saved,"HFP candidate stage mutated accepted source");
    std::cerr<<"HFP_OPERATION_BEGIN,commitActualEnd\n"<<std::flush;
    step.commit(*mesh,accepted);check(accepted.timeLayer().step==1&&accepted.phaseGeometry().vofPublication().pinnedCorner(),"HFP physical step failed to commit pin");
    auto fine=nozzle(1,false,true);const auto td=selected(*fine,true);const auto& endPub=accepted.phaseGeometry().vofPublication();
    std::cerr<<"HFP_OPERATION_BEGIN,sameEndRefine,source,"<<endPub.identity().value()<<'\n'<<std::flush;
    const auto refined=I::GeometricVofMeshTransfer::prepare(*mesh,endPub,*fine,td,end);
    check(refined.report().refinedTargets>0&&refined.candidate().identity()!=endPub.identity(),"HFP missing real same-End refinement");
    const auto fineLoad=loadFor(*fine,refined.candidate());check(fineLoad.curvatureReport().pinnedCornerHeightCells>0,"HFP refined fresh HF lost pin");++remaps;
    std::cerr<<"HFP_OPERATION_BEGIN,sameEndCoarsen,source,"<<refined.candidate().identity().value()<<'\n'<<std::flush;
    const auto coarse=I::GeometricVofMeshTransfer::prepare(*fine,refined.candidate(),*mesh,d,end);
    check(coarse.report().coarsenedTargets>0&&coarse.candidate().identity()!=refined.candidate().identity(),"HFP missing real same-End coarsening");
    const auto coarseLoad=loadFor(*mesh,coarse.candidate());check(coarseLoad.curvatureReport().pinnedCornerHeightCells>0,"HFP coarsened fresh HF lost pin");++remaps;
    std::cout<<"HFP_ACTUAL_PRESSURE_NOSLIP_STAGE,AC_stages,"<<actualStages<<",same_End_refine_coarsen,"<<remaps<<'\n';
}
int atLower(const QuadtreeMesh& mesh,Vector2D lower)
{
    for(int native=0;native<mesh.cellIndexNum();++native){const auto& cell=mesh.cells(mesh.cellIndexToLabel[native]);
        const auto p=mesh.verts(cell.cellToVert[0]).coord();if(p.x==lower.x&&p.y==lower.y)return native;}
    check(false,"finite-height fixture lost its actual native rectangle");return -1;
}
template<class F>void finiteReject(F action,const char* expected,const char* label)
{
    bool caught=false;std::string actual;
    try{action();}catch(const std::invalid_argument& error){caught=true;actual=error.what();}
    catch(const std::exception& error){throw std::runtime_error(std::string(label)+
        "; wrong exception type; actual_what={"+error.what()+"}");}
    catch(...){throw std::runtime_error(std::string(label)+"; wrong non-standard exception type");}
    if(!caught)throw std::runtime_error(std::string(label)+"; no invalid_argument was thrown");
    if(actual.find(expected)==std::string::npos)throw std::runtime_error(std::string(label)+
        "; expected_guard={"+expected+"}; actual_what={"+actual+"}");
    ++negative;
}
void finiteInventory(bool rz)
{
    auto mesh=nozzle(1,!rz);const auto d=selected(*mesh,rz);auto input=stocks(*mesh,d,Shape::Cylinder);
    const double h=rz?.125:.03125;
    // These are two genuine stock-mixed, measured-pure publication cells:
    // one pure-G end gains an unknown positive L area; a separate pure-L
    // beginning loses an unknown positive G area. No normal is hand-written.
    const int gas=atLower(*mesh,{1+3*h,3}),liquid=atLower(*mesh,{1-5*h,3+h});
    constexpr long double minor=5.9065793622204327e-112L;
    check(input[gas][0]==0&&input[liquid][1]==0,"finite-height fixture does not start at independent pure ends");
    input[gas][0]=minor;input[liquid][1]=minor;
    const Clock clock{7,0,.0001,Clock::Position::End};const auto pub=pinned(*mesh,d,input,clock);
    const auto saved=isolatedFragmentSnapshot(pub);
    for(int native:{gas,liquid}){const auto& c=pub.plicCells()[native];
        check(!c.resolved&&c.inventory[0]>0&&c.inventory[1]>0&&
            (c.measured[0]==0||c.measured[1]==0),"finite-height fixture did not use actual published unresolved pure geometry");}
    const auto load=loadFor(*mesh,pub);algebra(*mesh,load);unchanged(pub,input);
    check(isolatedFragmentSnapshot(pub)==saved,"finite-height certificate changed source stock/geometry/clock");
    const auto& report=load.curvatureReport();
    check(report.finiteInventoryHeightColumns>0&&report.finiteInventoryHeightQueries>0&&
        !report.finiteInventoryHeightCertificates.empty(),"finite-height actual load did not publish its successful height certificate");
    const long double ulp=std::nextafter(1.L,2.L)-1.L;std::size_t columns=0,queries=0;
    for(const auto& receipt:report.finiteInventoryHeightCertificates){
        check(receipt.heightAxis==0&&receipt.outward==1&&receipt.nativeInverseSpacing==1.L/receipt.tangentSpacing&&
            receipt.queryCount>0,"finite-height receipt lost its actual radial chart coefficients");
        // Independent planar-column integral of the exact cylinder r=1.
        // This uniform binary mesh makes all endpoint differences/products
        // exact; both unknown-area bounds are far inside the same RN bin.
        for(int column=0;column<3;++column)if(receipt.certified[column]){
            long double thickness=0;
            for(int line=-4;line<=4;++line){
                const long double lo=receipt.center.x+(line-.5L)*receipt.normalSpacing,hi=lo+receipt.normalSpacing;
                if(lo<1)thickness+=std::max(0.L,std::min(1.L,hi)-lo);
            }
            const long double independent=receipt.origin+thickness;
            check(independent==1&&receipt.physicalHeights[column]==independent&&!std::signbit(receipt.physicalHeights[column]),
                "finite-height canonical RN differs from the independent actual cylinder integral");
            ++columns;
        }
        queries+=receipt.queryCount;
    }
    for(int native:{gas,liquid}){const auto& r=pub.plicCells()[native].rectangle;
        const long double upper=rz?minor/(static_cast<long double>(CartesianMeasure::twoPi)*r.lower.x*h):minor/h;
        check(upper>0&&upper<ulp/8,"finite-height independent unknown-area bound is not inside its declared RN bin");}
    check(columns==report.finiteInventoryHeightColumns&&queries==report.finiteInventoryHeightQueries,
        "finite-height typed receipt does not account for its actual successful columns/queries");
    unsigned forceFaces=0;const auto& l=load.phaseGeometry().layout();
    for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(*mesh,l.nativeFace(u));
        if(!f.horizontal&&f.center.x==1&&f.center.y>3&&f.center.y<3+2*h&&load.colorGradient()[u]!=0){
            check(load.supportReports()[u].coincidentPhaseHeight&&std::abs(load.curvatureAverages()[u]-(rz?1:0))<=1e-11,
                "finite-height cylinder changed the original physical curvature gate");++forceFaces;}}
    check(forceFaces>0,"finite-height unique RN was not consumed by the actual positive-sigma force");
    ++scenes;std::cout<<"HFP_FINITE_INVENTORY_HEIGHT,"<<(rz?"RZ":"XY")<<",columns,"<<columns<<",queries,"<<queries
        <<",canonical_height,1,minor_geometry_resolved,0,source_unchanged,1\n";
}
void finiteInventoryNegatives()
{
    auto mesh=nozzle(1);const auto d=selected(*mesh,true);const auto pure=stocks(*mesh,d,Shape::Cylinder);
    const Clock clock{8,0,.0001,Clock::Position::End};const int gas=atLower(*mesh,{1.375,3});
    const long double ulp=std::nextafter(1.L,2.L)-1.L,h=.125L,rlo=1.375L,rhi=1.5L;
    auto straddling=pure;
    straddling[gas][0]=static_cast<long double>(CartesianMeasure::twoPi)*((rlo+rhi)/2)*h*(ulp/2);
    const auto uncertain=pinned(*mesh,d,straddling,clock);const auto saved=isolatedFragmentSnapshot(uncertain);
    check(!uncertain.plicCells()[gas].resolved&&uncertain.plicCells()[gas].measured[0]==0,
        "finite-height RN-straddle control was not a genuine unresolved pure placeholder");
    const auto q=straddling[gas][0];
    check(q/(static_cast<long double>(CartesianMeasure::twoPi)*rhi*h)<ulp/2&&
        q/(static_cast<long double>(CartesianMeasure::twoPi)*rlo*h)>ulp/2,
        "finite-height independent RZ interval does not straddle its native halfway boundary");
    finiteReject([&]{loadFor(*mesh,uncertain);},"unique native nearest-even height",
        "finite-height certificate accepted a genuine RN-straddling interval");
    check(isolatedFragmentSnapshot(uncertain)==saved,"RN-straddle rejection changed independent source geometry");
    unchanged(uncertain,straddling);
    // A genuine symmetric isolated mixed stock has zero reconstructed
    // node-alpha gradient, but nonzero measured phases. It is not a pure
    // boundary-placeholder and must retain the original unresolved guard.
    auto nonpure=pure;const int mixed=atLower(*mesh,{1.375,3.125});
    const auto capacity=nonpure[mixed][1];nonpure[mixed]={capacity/2,capacity/2};
    const auto unresolved=pinned(*mesh,d,nonpure,clock);
    check(!unresolved.plicCells()[mixed].resolved&&unresolved.plicCells()[mixed].measured[0]>0&&unresolved.plicCells()[mixed].measured[1]>0,
        "finite-height non-placeholder control has no genuine unresolved mixed geometry");
    finiteReject([&]{loadFor(*mesh,unresolved);},"local PLIC height-function requires a resolved mixed PLIC source:",
        "finite-height certificate accepted unresolved non-placeholder geometry");
    unchanged(unresolved,nonpure);
    // This real coarser radial chart includes r=0 in its first L end.
    // A positive unknown gas stock there has no finite positive-r bound.
    DictionaryMesh cfg;cfg.Nx=20;cfg.Ny=128;cfg.minCoord={0,0};cfg.maxCoord={5,16};
    cfg.maxCellNum=100000;cfg.maxInitialIterations=0;QuadtreeMesh axisMesh(cfg);const auto axisDomain=selected(axisMesh,true);
    auto axisInput=stocks(axisMesh,axisDomain,Shape::Cylinder);const int axis=atLower(axisMesh,{0,3});
    check(axisInput[axis][1]==0,"finite-height axis control does not begin as pure liquid");
    axisInput[axis][1]=5.9065793622204327e-112L;
    const auto axisPub=pinned(axisMesh,axisDomain,axisInput,clock);
    check(!axisPub.plicCells()[axis].resolved&&axisPub.plicCells()[axis].measured[1]==0,
        "finite-height axis control has no actual pure-L placeholder");
    finiteReject([&]{loadFor(axisMesh,axisPub);},"positive actual radial lower bound",
        "finite-height certificate invented a bound at the actual RZ axis");
    unchanged(axisPub,axisInput);
    std::cout<<"HFP_FINITE_INVENTORY_REJECTS,RN_straddle,1,non_placeholder,1,axis,1\n";
}
void negatives()
{
    auto mesh=nozzle(1);const auto d=selected(*mesh,true);const auto input=stocks(*mesh,d,Shape::Hemicap);
    const Clock end{2,0,.0001,Clock::Position::End};const auto pub=pinned(*mesh,d,input,end);const auto saved=isolatedFragmentSnapshot(pub);
    const auto hd=selected(*mesh,true,true);const auto bad=pinned(*mesh,hd,stocks(*mesh,hd,Shape::Hemicap),end);
    reject([&]{loadFor(*mesh,bad);},"HFP invented phase support in a selected solid coverage hole");
    auto reversed=input;
    for(int ci:d.cellIndexByRow()){const auto p=mesh->cells(mesh->cellIndexToLabel[ci]).centerCoord();
        if(p.x>.5&&p.x<.625&&p.y>3.125&&p.y<3.25){const auto c=pub.cells()[ci].capacity;reversed[ci]={0,c};}}
    const auto branch=pinned(*mesh,d,reversed,end);reject([&]{loadFor(*mesh,branch);},"HFP accepted ambiguous multiple radial phase transitions");
    const auto raw=Vof::prepareLocalPlicInventory(*mesh,d,end,input,{});
    reject([&]{loadFor(*mesh,raw);},"HFP unsealed corner gained an implicit pin boundary law");
    auto remote=input;
    for(int ci:d.cellIndexByRow()){const auto p=mesh->cells(mesh->cellIndexToLabel[ci]).centerCoord();const auto c=pub.cells()[ci].capacity;
        if(p.y>3&&p.y<3.125&&p.x>1.625&&p.x<1.75)remote[ci]={c,0};
        if(p.y>3&&p.y<3.125&&p.x>1.75&&p.x<1.875)remote[ci]={.37L*c,.63L*c};}
    const auto far=pinned(*mesh,d,remote,end);bool targeted=false;
    try{loadFor(*mesh,far);}catch(const std::invalid_argument& error){
        check(std::string(error.what()).find("wall-blocked")!=std::string::npos,"HFP remote fragment did not reach its intended pin-window guard");targeted=true;}
    check(targeted,"HFP same-bin remote interface borrowed the pin's fixed normal window");++negative;
    check(isolatedFragmentSnapshot(pub)==saved,"HFP failed support probes mutated original source");unchanged(pub,input);
}
void run()
{
    arithmetic();std::array<double,3> errors{};int i=0;for(int scale:{1,2,4})errors[i++]=curve(scale,true);
    check(errors[2]<errors[0]&&errors[2]<.1,"HFP actual PLIC/HF pin curvature did not converge within its declared engineering budget");
    curve(1,false);curve(2,false);for(bool rz:{false,true}){pure(rz,Shape::Cylinder);pure(rz,Shape::Plane);}
    finiteInventory(false);finiteInventory(true);finiteInventoryNegatives();negatives();stepAndAmr();
    check(actualStages>0&&remaps==2&&negative>=18,"HFP actual consumer or negative witnesses missing");
    std::cout<<"HFP_PINNED_CORNER_PASS,scenes,"<<scenes<<",guards,"<<negative<<",AC_stages,"<<actualStages<<",remaps,"<<remaps
        <<",actual_RZ_kappa_errors,"<<errors[0]<<','<<errors[1]<<','<<errors[2]<<",C_force,"<<maxForce<<",NO_LONG_WINDOW_CLAIM\n";
}
}
namespace DyadicWidthContract
{
namespace I=LAMP::Interface;
std::unique_ptr<QuadtreeMesh> actualMesh(bool perturb)
{
    DictionaryMesh cfg;cfg.Nx=40;cfg.Ny=80;cfg.minCoord={0,0};cfg.maxCoord={5,16};
    cfg.maxCellNum=100000;cfg.maxInitialIterations=0;cfg.maxLevel_local={4};
    auto mesh=std::make_unique<QuadtreeMesh>(cfg);
    DictionaryMesh::RefineBox coarse;coarse.level=3;coarse.min={0,4.4};coarse.max={.25,4.8};
    DictionaryMesh::RefineBox fine;fine.level=4;fine.min={0,4.6};fine.max={.125,4.7};
    DictionaryMesh::RefineBox mixedPatch;mixedPatch.level=4;mixedPatch.min={0,2.75};mixedPatch.max={.375,3.35};
    // Use the same Core target/prerequisite/adapt route as a real graded
    // tree. applyRefineBoxes alone does not complete the box's 2:1 halo.
    unsigned sweeps=0;
    for(;;){
        LAMP::Mesh::NativeTargetLevelField target(*mesh,0.),gradient(*mesh,0.);
        const auto inBox=[](Vector2D q,const DictionaryMesh::RefineBox& box){
            return q.x>=box.min.x&&q.x<=box.max.x&&q.y>=box.min.y&&q.y<=box.max.y;};
        for(int ci=0;ci<mesh->cellIndexNum();++ci){const auto& cell=mesh->cells(mesh->cellIndexToLabel[ci]);const auto q=cell.centerCoord();
            target.value[ci]=std::max(cell.level(),inBox(q,mixedPatch)?mixedPatch.level:inBox(q,fine)?fine.level:inBox(q,coarse)?coarse.level:0);}
        mesh->completeRefinementPrerequisites(target);
        if(!mesh->mayAdaptToTargets(target))break;
        check(sweeps<static_cast<unsigned>(mesh->maxLevel()+2),"dyadic-width balanced target closure did not terminate");
        const int before=mesh->cellIndexNum();
        std::vector<VertField<double>> fields;std::vector<VertField<double>*> vertices;
        std::vector<HEdgeField<double>*> horizontal;std::vector<VEdgeField<double>*> vertical;
        std::vector<LAMP::Mesh::NativeTargetLevelField*> cells{&target};
        mesh->adapt(fields,gradient,vertices,horizontal,vertical,cells,&target);++sweeps;
        check(mesh->cellIndexNum()>before,"dyadic-width balanced refinement made no actual progress");
    }
    check(sweeps>0,"dyadic-width fixture did not refine the real native mesh");
    if(perturb){
        unsigned changed=0;
        // Shift this one actual horizontal vertex layer by one binary64 ulp,
        // retaining rectangular leaf faces but breaking the canonical chart.
        for(int v=0;v<mesh->vertIndexNum();++v){auto& vertex=mesh->verts(mesh->vertIndexToLabel[v]);auto q=vertex.coord();
            if(q.y==4.6125){q.y=std::nextafter(q.y,std::numeric_limits<double>::infinity());vertex.setCoord(q);++changed;}}
        check(changed>0,"dyadic-width negative control did not perturb any actual corner");
    }
    return mesh;
}
void run()
{
    auto mesh=actualMesh(false);const auto metrics=metricsFor(*mesh,true);const auto& layout=metrics.layout();
    const LAMP::TimeLayer clock{1,0,1e-6,LAMP::TimeLayer::Position::Midpoint};
    const auto pub=I::GeometricVofPublication::prepareLocalPlic(*mesh,layout.domain(),clock,std::vector<double>(mesh->cellIndexNum(),0));
    unsigned widthFailures=0,witnesses=0;
    const double roundoff=64*std::numeric_limits<double>::epsilon();
    for(const auto f:layout.domain().faces(*mesh)){
        if(f.negativeCell<0||f.positiveCell<0)continue;
        const auto& a=mesh->cells(mesh->cellIndexToLabel[f.negativeCell]);const auto& b=mesh->cells(mesh->cellIndexToLabel[f.positiveCell]);
        check(std::abs(a.level()-b.level())<=1,"dyadic-width fixture is not a real locally 2:1 tree");
        const auto& ar=pub.plicCells()[f.negativeCell].rectangle;const auto& br=pub.plicCells()[f.positiveCell].rectangle;
        const auto ax=ar.upper.x-ar.lower.x,ay=ar.upper.y-ar.lower.y,bx=br.upper.x-br.lower.x,by=br.upper.y-br.lower.y;
        const bool oldGate=std::max(ax,bx)<=2*std::min(ax,bx)*(1+roundoff)&&std::max(ay,by)<=2*std::min(ay,by)*(1+roundoff);
        if(!oldGate)++widthFailures;
        const auto& small=a.level()>b.level()?ar:br;const auto& large=a.level()>b.level()?br:ar;
        if(!oldGate&&std::min(a.level(),b.level())==3&&std::max(a.level(),b.level())==4&&
            small.lower.y==4.6000000000000005&&small.upper.y==4.6125&&large.lower.y==small.lower.y&&large.upper.y==4.625){
            check(small.upper.y-small.lower.y==.01249999999999929&&large.upper.y-large.lower.y==.024999999999999467,
                "real 40x80/L4 z4.6 canonical endpoint witness changed");
            ++witnesses;
        }
    }
    check(widthFailures>0&&witnesses>0,"dyadic-width fixture never exercises the actual failed guard");
    const auto phase=HydroPhaseGeometry::build(*mesh,metrics,pub);const auto pressure=HydroPressureLoad::build(metrics.divergence());
    const auto load=HydroBalancedCsfLoad::localPlicHeightFunctionCapillary(*mesh,pressure,phase,1,heightConfig());
    load.validateVofFor(*mesh,pressure,pub.identity(),clock);
    for(double value:load.momentumRhsLoad())check(value==0,"pure actual HF qualification invented a capillary load");
    for(const auto& cell:pub.plicCells())check(cell.inventory[0]==0,"HF dyadic qualification changed independent source stock");
    // Real closed RZ material with a straight cylinder through the exact
    // failing z3.025 layer. Both physical phase stocks use independent integrals.
    constexpr long double radius=.12109375L;
    const auto stocks=rzAnalyticStocks(*mesh,[radius](Vector2D a,Vector2D b){return rzCapsuleColumn(a,b,2.98L,3.08L,radius);});
    const auto mixedPub=I::GeometricVofPublication::prepareLocalPlicInventory(*mesh,layout.domain(),clock,stocks,{});
    int self=-1;
    for(int native=0;native<mesh->cellIndexNum();++native){const auto& r=mixedPub.plicCells()[native].rectangle;
        if(r.lower.x==.1171875&&r.upper.x==.125&&r.lower.y==3.0250000000000004&&r.upper.y==3.0375000000000001)self=native;}
    check(self>=0,"canonical-query test lost its actual mixed z3.025 native leaf");
    const auto& cell=mixedPub.plicCells()[self];const auto& r=cell.rectangle;const double actualDy=r.upper.y-r.lower.y;
    check(cell.resolved&&cell.measured[0]>0&&cell.measured[1]>0,"canonical-query source is not an actual resolved mixed plane");
    volatile double storedCenter=(r.lower.y+r.upper.y)*.5,storedHalf=actualDy*.5;
    const double oldLower=storedCenter-storedHalf,oldUpper=storedCenter+storedHalf,oldWidth=oldUpper-oldLower;
    check(oldUpper==3.0374999999999996&&oldUpper!=r.upper.y&&actualDy>oldWidth*(1+roundoff),
        "canonical-query witness no longer exercises original self partial/whole rejection");
    const LAMP::Mesh::DyadicGridCoordinates chart(mesh->minCoord,mesh->maxCoord,mesh->Nx,mesh->Ny);
    const auto address=LAMP::Mesh::cellAddress(*mesh,mesh->cells(mesh->cellIndexToLabel[self]));
    check(chart.point(address.i,address.j,address.level).y==r.lower.y&&
        chart.point(address.i+1,address.j+1,address.level).y==r.upper.y,"actual mixed source lost exact Core chart endpoints");
    const auto mixedPhase=HydroPhaseGeometry::build(*mesh,metrics,mixedPub);
    const auto mixedLoad=HydroBalancedCsfLoad::localPlicHeightFunctionCapillary(*mesh,pressure,mixedPhase,1,heightConfig());
    mixedLoad.validateVofFor(*mesh,pressure,mixedPub.identity(),clock);
    const int selfPressure=layout.pressureDof(self);const auto& jet=mixedLoad.heightFunctionCellJets()[selfPressure];
    check(jet&&jet->heightAxis==0&&!jet->coarseAligned&&!jet->pinnedCornerConstrained&&jet->tangentSpacing==actualDy&&
        mixedLoad.heightFunctionCellCurvatureMethods()[selfPressure]==HydroCsfCellCurvatureMethod::HeightFunction,
        "actual mixed self did not consume its original-spacing complete radial height chart");
    check(std::abs(mixedLoad.heightFunctionCellCurvatures()[selfPressure]-1/static_cast<double>(radius))<1e-6,
        "actual mixed cylinder lost its independent physical hoop curvature");
    unsigned mixedForce=0;for(double value:mixedLoad.momentumRhsLoad()){check(std::isfinite(value),"canonical-query force is nonfinite");mixedForce+=value!=0;}
    check(mixedForce>0,"canonical-query positive-sigma actual mixed load was not consumed");
    for(std::size_t native=0;native<stocks.size();++native)check(mixedPub.plicCells()[native].inventory==stocks[native],
        "canonical-query read changed an independent physical phase stock");
    // Closed real RZ source whose upper cap crosses the original 40x80
    // near-axis z4..4.2 cell. A reflected radial tangent must retain exactly
    // the same normal-Y endpoint as the two nonnegative tangent columns.
    DictionaryMesh axisCfg;axisCfg.Nx=40;axisCfg.Ny=80;axisCfg.minCoord={0,0};axisCfg.maxCoord={5,16};
    axisCfg.maxCellNum=100000;axisCfg.maxInitialIterations=0;axisCfg.maxLevel_local={4};
    QuadtreeMesh axisMesh(axisCfg);const auto axisMetrics=metricsFor(axisMesh,true);const auto& axisLayout=axisMetrics.layout();
    const auto axisStocks=rzAnalyticStocks(axisMesh,[](Vector2D a,Vector2D b){return rzSphereColumn(a,b,2.93L,1.11L);});
    // Supply the independent analytic fraction to the Interface owner, whose
    // strict actual capacity initializes both stocks; a separate rounded RZ
    // integral is not that owner's exact whole-cell capacity.
    std::vector<double> axisAlpha;axisAlpha.reserve(axisStocks.size());
    for(const auto& stock:axisStocks)axisAlpha.push_back(static_cast<double>(stock[0]/(stock[0]+stock[1])));
    const auto axisPub=I::GeometricVofPublication::prepareLocalPlic(axisMesh,axisLayout.domain(),clock,axisAlpha);
    std::vector<RzStocks> axisOwnedStocks;axisOwnedStocks.reserve(axisStocks.size());
    for(const auto& cell:axisPub.plicCells())axisOwnedStocks.push_back(cell.inventory);
    int axisSource=-1;
    for(int native=0;native<axisMesh.cellIndexNum();++native){const auto& box=axisPub.plicCells()[native].rectangle;
        if(box.lower.x==0&&box.upper.x==.125&&box.lower.y==4&&box.upper.y==4.2)axisSource=native;}
    check(axisSource>=0,"canonical-axis test lost its actual near-axis z4..4.2 source");
    const auto& axisCell=axisPub.plicCells()[axisSource];const auto& axisBox=axisCell.rectangle;
    check(axisCell.resolved&&axisCell.measured[0]>0&&axisCell.measured[1]>0&&std::abs(axisCell.normal.y)>std::abs(axisCell.normal.x),
        "canonical-axis source is not a resolved mixed normal-Y cap");
    volatile double oldAxisCenter=(axisBox.lower.y+axisBox.upper.y)*.5,oldAxisDy=axisBox.upper.y-axisBox.lower.y;
    volatile double oldAxisShift=-4*oldAxisDy,oldAxisQueryCenter=oldAxisCenter+oldAxisShift;
    const double oldAxisOrigin=oldAxisQueryCenter-oldAxisDy/2;
    const LAMP::Mesh::DyadicGridCoordinates axisChart(axisMesh.minCoord,axisMesh.maxCoord,axisMesh.Nx,axisMesh.Ny);
    const double actualAxisOrigin=axisChart.point(0,16,0).y;double threeUlps=oldAxisOrigin;
    for(int step=0;step<3;++step)threeUlps=std::nextafter(threeUlps,std::numeric_limits<double>::infinity());
    check(oldAxisOrigin==3.1999999999999988&&threeUlps==actualAxisOrigin&&actualAxisOrigin!=oldAxisOrigin,
        "canonical-axis test no longer exposes the three-ULP reflected-column origin mismatch");
    const auto axisPhase=HydroPhaseGeometry::build(axisMesh,axisMetrics,axisPub);
    const auto axisPressure=HydroPressureLoad::build(axisMetrics.divergence());
    const auto axisLoad=HydroBalancedCsfLoad::localPlicHeightFunctionCapillary(axisMesh,axisPressure,axisPhase,1,heightConfig());
    axisLoad.validateVofFor(axisMesh,axisPressure,axisPub.identity(),clock);
    const int axisP=axisLayout.pressureDof(axisSource);const auto& axisJet=axisLoad.heightFunctionCellJets()[axisP];
    check(axisJet&&axisJet->heightAxis==1&&axisJet->axisReflected&&!axisJet->coarseAligned&&!axisJet->pinnedCornerConstrained&&
        axisJet->origin==actualAxisOrigin&&axisJet->normalSpacing==axisBox.upper.y-axisBox.lower.y&&axisJet->tangentSpacing==.125,
        "near-axis actual normal-Y height lost its original nine-bin canonical reflected origin");
    check(std::isfinite(axisLoad.heightFunctionCellCurvatures()[axisP])&&axisLoad.heightFunctionCellCurvatures()[axisP]>0,
        "near-axis actual closed sphere has no finite positive physical curvature");
    unsigned axisForce=0;for(double value:axisLoad.momentumRhsLoad()){check(std::isfinite(value),"canonical-axis force is nonfinite");axisForce+=value!=0;}
    check(axisForce>0,"canonical-axis positive-sigma mixed force was not consumed");
    for(std::size_t native=0;native<axisOwnedStocks.size();++native)check(axisPub.plicCells()[native].inventory==axisOwnedStocks[native],
        "canonical-axis read changed either Interface-owned physical phase stock");
    auto bad=actualMesh(true);const auto badMetrics=metricsFor(*bad,true);
    const auto badPub=I::GeometricVofPublication::prepareLocalPlic(*bad,badMetrics.layout().domain(),clock,std::vector<double>(bad->cellIndexNum(),0));
    const auto badPhase=HydroPhaseGeometry::build(*bad,badMetrics,badPub);const auto badPressure=HydroPressureLoad::build(badMetrics.divergence());
    bool rejected=false;
    try{(void)HydroBalancedCsfLoad::localPlicHeightFunctionCapillary(*bad,badPressure,badPhase,1,heightConfig());}
    catch(const std::invalid_argument& error){rejected=std::string(error.what()).find("canonical dyadic corners")!=std::string::npos;}
    check(rejected,"HF accepted a perturbed actual corner instead of authenticating its exact dyadic chart");
    std::cout<<"HF_DYADIC_WIDTH_PASS,old_width_failures,"<<widthFailures<<",actual_z4.6_witnesses,"<<witnesses
        <<",pure_positive_sigma_consumer,1,mixed_self_partial_witness,1,mixed_original_spacing_height,1,mixed_force_faces,"<<mixedForce<<",axis_three_ulp_witness,1,axis_reflected_normal_y_height,1,axis_force_faces,"<<axisForce<<",actual_corner_rejection,1,NO_PHYSICAL_PARAMETER_CHANGE\n";
}
}
namespace HfWorkspaceContract
{
void run()
{
    // An empty active-height set must retain the original dense public outputs.
    // Both phases and both measures exercise pure rows before any mixed slot.
    for(bool rz:{false,true})for(double alpha:{0.,1.}){
        auto mesh=makeMesh(16,1);const auto metrics=metricsFor(*mesh,rz);const auto& layout=metrics.layout();
        const LAMP::TimeLayer clock{3,.001,.0011,LAMP::TimeLayer::Position::Midpoint};
        const auto pub=LAMP::Interface::GeometricVofPublication::prepareLocalPlic(*mesh,layout.domain(),clock,
            std::vector<double>(mesh->cellIndexNum(),alpha));
        const auto source=isolatedFragmentSnapshot(pub);
        const auto phase=HydroPhaseGeometry::build(*mesh,metrics,pub);const auto pressure=HydroPressureLoad::build(metrics.divergence());
        const auto load=HydroBalancedCsfLoad::localPlicHeightFunctionCapillary(*mesh,pressure,phase,.7,heightConfig());
        load.validateVofFor(*mesh,pressure,pub.identity(),clock);
        check(load.heightFunctionCellCurvatures().size()==static_cast<std::size_t>(layout.pressureCount()),
            "pure HF workspace lost its dense public pressure rows");
        for(int p=0;p<layout.pressureCount();++p){
            check(!std::isfinite(load.heightFunctionCellCurvatures()[p])&&!std::isfinite(load.heightFunctionCellEffectiveCurvatures()[p])
                &&!load.heightFunctionCellJets()[p]
                &&load.heightFunctionCellCurvatureMethods()[p]==HydroCsfCellCurvatureMethod::Unavailable,
                "pure HF workspace fabricated a mixed-cell height or jet");}
        for(double value:load.momentumRhsLoad())check(value==0,"pure HF workspace invented capillary force");
        const auto& report=load.curvatureReport();
        check(report.heightFunctionCells==0&&report.heightFunctionColumns==0&&report.heightFunctionCellQueries==0,
            "pure HF workspace consumed an actual height query");
        check(isolatedFragmentSnapshot(pub)==source,"pure HF workspace changed the source inventories or geometry");
    }
    // Reuse independent analytic/load checks, sparse mixed/native row mappings,
    // graded 2:1 charts, stock-only microphase, true RZ axis and old wall guards.
    localPlicHeightFunction();localPlicGradedHeightFunction();localPlicRzHeightFunction();
    std::cout<<"LOCAL_PLIC_WORKSPACE_PASS,pure_scenes,4,analytic_graded_RZ_original_groups,1,NO_PARALLEL_OR_TIMING_CLAIM\n";
}
}
namespace NativeSupportHeightChartContract
{
void mixedHeightChartGuard()
{
    namespace I=LAMP::Interface;constexpr long double h=1.L/64,bottom=31.75L*h,sourceX=32.5L*h,A=1.L/6,B=-5.L/16;
    constexpr double sigma=.7;const auto graph=[=](long double s){return h*(33+A+B*s+s*s);};
    const auto leftTop=graph(-1.5L),rightTop=graph(1.5L),leftRadius=(leftTop-bottom)/2,rightRadius=(rightTop-bottom)/2;
    const auto polynomialArea=[&](long double x0,long double x1,long double y0,long double y1){
        if(!(x1>x0&&y1>y0))return 0.L;std::vector<long double> sites{x0,x1};
        for(auto y:{y0,y1}){const auto discriminant=B*B-4*(33+A-y/h);if(discriminant<0)continue;
            for(auto s:{(-B-std::sqrt(discriminant))/2,(-B+std::sqrt(discriminant))/2}){const auto x=sourceX+h*s;if(x>x0&&x<x1)sites.push_back(x);}}
        std::sort(sites.begin(),sites.end());long double area=0;
        for(std::size_t k=1;k<sites.size();++k){const auto lo=sites[k-1],hi=sites[k],mid=graph(((lo+hi)/2-sourceX)/h);
            if(mid<=y0)continue;if(mid>=y1){area+=(hi-lo)*(y1-y0);continue;}
            const auto s0=(lo-sourceX)/h,s1=(hi-sourceX)/h;
            area+=h*h*((33+A-y0/h)*(s1-s0)+B*(s1*s1-s0*s0)/2+(s1*s1*s1-s0*s0*s0)/3);}
        return area;};
    const auto independentArea=[&](Vector2D lo,Vector2D hi){long double area=0;
        if(lo.x<.25)area+=circleCell(lo,{std::min(hi.x,.25),hi.y},{.25,static_cast<double>(bottom+leftRadius)},static_cast<double>(leftRadius),false);
        if(hi.x>.75)area+=circleCell({std::max(lo.x,.75),lo.y},hi,{.75,static_cast<double>(bottom+rightRadius)},static_cast<double>(rightRadius),false);
        const auto x0=std::max(static_cast<long double>(lo.x),.25L),x1=std::min(static_cast<long double>(hi.x),.75L),y0=std::max(static_cast<long double>(lo.y),bottom),y1=static_cast<long double>(hi.y);
        if(!(x1>x0&&y1>y0))return area;const auto join0=sourceX-1.5L*h,join1=sourceX+1.5L*h;
        const auto strip=[&](long double l,long double r,long double top){return r>l?std::max(0.L,std::min(y1,top)-y0)*(r-l):0.L;};
        area+=strip(x0,std::min(x1,join0),leftTop);area+=polynomialArea(std::max(x0,join0),std::min(x1,join1),y0,y1);
        area+=strip(std::max(x0,join1),x1,rightTop);return area;};
    auto mesh=heightMesh(64);FaceSegmentGraph graphDomain;graphDomain.build(*mesh);const auto domain=NativeLeafDomain::build(*mesh,graphDomain,false);
    std::map<std::pair<int,int>,int> grid;std::vector<std::array<long double,2>> stocks(mesh->cellIndexNum());
    for(int ci=0;ci<mesh->cellIndexNum();++ci){const auto& cell=mesh->cells(mesh->cellIndexToLabel[ci]);
        const auto lo=mesh->verts(cell.cellToVert[0]).coord(),hi=mesh->verts(cell.cellToVert[2]).coord();
        const auto capacity=static_cast<long double>(hi.x-lo.x)*(hi.y-lo.y),liquid=independentArea(lo,hi);
        check(liquid>=0&&liquid<=capacity,"mixed-chart invalid setup: independent closed polynomial-band stock escaped capacity");stocks[ci]={liquid,capacity-liquid};
        grid[{static_cast<int>(std::llround(lo.x/h)),static_cast<int>(std::llround(lo.y/h))}]=ci;}
    const int source=grid.at({32,33});const I::GeometricVofTimeLayer clock{51,.002,.0022,LAMP::TimeLayer::Position::Midpoint};I::GeometricVofConfig config;
    config.localResolution=I::GeometricVofLocalResolution::NativeInterfaceSupport;
    const auto publication=I::GeometricVofPublication::prepareLocalPlicInventory(*mesh,domain,clock,stocks,config);
    const auto saved=isolatedFragmentSnapshot(publication);const auto id=publication.identity();const NativeMeshIdentity meshIdentity(*mesh);
    const auto mask=publication.localPlicSubgridMask(*mesh);const auto& sourceCell=publication.plicCells()[source];
    check(!mask[source]&&sourceCell.resolved&&sourceCell.measured[0]>0&&sourceCell.measured[1]>0&&sourceCell.normal.y>0,
        "mixed-chart invalid setup: steep-V source is not an actual unmasked upward two-phase PLIC cut");
    const auto metrics=HydroPhysicalDivergence::build(*mesh,HydroIntegralDivergence::build(*mesh,HydroDofLayout::build(*mesh,publication.domain())));
    const auto phase=HydroPhaseGeometry::build(*mesh,metrics,publication);const auto pressure=HydroPressureLoad::build(metrics.divergence());const auto& layout=phase.layout();
    const int row=layout.pressureDof(source);LocalPlicEndpointReference independent{phase};const auto midpoint=independent.segmentMidpoint(sourceCell);
    check(!independent.estimate(row,0).valid&&!independent.estimate(row,1).valid&&independent.estimate(row,1).opposite,
        "mixed-chart invalid setup: full original source charts did not fail on the actual remote bottom branch");
    std::array<long double,3> worldHeight{};const std::array<long double,3> exact{{34.5625L,33.25L,33.9375L}};
    for(int column=-1;column<=1;++column){const auto& first=publication.plicCells()[grid.at({32+column,32})];
        const auto& last=publication.plicCells()[grid.at({32+column,36})];
        check(first.measured[1]==0&&last.measured[0]==0,"mixed-chart invalid setup: local top columns lack complete liquid/gas ends");
        long double liquid=0,previous=1;
        for(int y=32;y<=36;++y){const int ci=grid.at({32+column,y});const auto& q=publication.plicCells()[ci];
            const auto fraction=q.measured[0]/(h*h);check(fraction<=previous+1e-6L,"mixed-chart invalid setup: local top column is not ordered");previous=fraction;
            if(q.measured[0]>0&&q.measured[1]>0)check(q.resolved&&(mask[ci]||q.normal.y>0),"mixed-chart invalid setup: a local top column has an unrelated resolved branch");
            liquid+=q.measured[0];}
        worldHeight[column+1]=32*h+liquid/h;
        const auto originalBudget=5e-12L+256*std::numeric_limits<double>::epsilon()*40;
        close(worldHeight[column+1]/h,exact[column+1],originalBudget,"mixed-chart actual independent top height matches the analytic stock integral");}
    const auto left=(worldHeight[1]-worldHeight[0])/h,right=(worldHeight[2]-worldHeight[1])/h;
    check(left<0&&right>0&&std::abs(left)>1&&std::abs((left+right)/2)<1,
        "mixed-chart invalid setup: real local height points lack the frozen opposed-steep cancellation");
    // In the alternative x chart, row32 is liquid throughout native radius4;
    // row33 also has positive liquid in every native cell. Neither can supply
    // a pure gas end. Thus the rejected y triple cannot quietly be replaced by
    // three independent same-source x-height points.
    for(int x=28;x<=36;++x){check(publication.plicCells()[grid.at({x,32})].measured[1]==0,
            "mixed-chart invalid setup: lower alternative x strip unexpectedly reaches a gas end");
        check(publication.plicCells()[grid.at({x,33})].measured[0]>1e-6L*h*h,
            "mixed-chart invalid setup: central alternative x strip unexpectedly reaches a gas end");}
    const auto load=HydroBalancedCsfLoad::localPlicHeightFunctionCapillary(*mesh,pressure,phase,sigma,heightConfig());
    load.validateVofFor(*mesh,pressure,id,clock);const auto& methods=load.heightFunctionCellCurvatureMethods();const auto& physical=load.heightFunctionCellCurvatures();
    unsigned noCompleteDonors=0;const double half=2.5*static_cast<double>(h);
    for(int p=0;p<layout.pressureCount();++p){const auto& q=publication.plicCells()[layout.nativeCell(p)];
        if(q.measured[0]<=0||q.measured[1]<=0||static_cast<long double>(q.normal.x)*sourceCell.normal.x+static_cast<long double>(q.normal.y)*sourceCell.normal.y<=0)continue;
        if(q.rectangle.upper.x<=midpoint.x-half||q.rectangle.lower.x>=midpoint.x+half||q.rectangle.upper.y<=midpoint.y-half||q.rectangle.lower.y>=midpoint.y+half)continue;
        check(!independent.estimate(p,0).valid&&!independent.estimate(p,1).valid&&!independent.estimate(p,0,-8,8).valid&&!independent.estimate(p,1,-8,8).valid,
            "mixed-chart invalid setup: fixed5h neighborhood acquired an independently complete-HF donor");
        check(methods[p]!=HydroCsfCellCurvatureMethod::HeightFunction,"mixed-chart invalid setup: immutable qualified-HF donor would hide the intended mixed-point failure");++noCompleteDonors;}
    check(noCompleteDonors>=3,"mixed-chart invalid setup: fixed5h donor exclusion lacks actual segment support");
    const auto k=physical[row];
    if(methods[row]!=HydroCsfCellCurvatureMethod::SegmentFit||!std::isfinite(k))std::cerr<<"MIXED_CHART_ACTUAL_METHOD_DIAGNOSTIC,native,"<<source<<",method,"<<static_cast<int>(methods[row])
        <<",kappa,"<<k<<",normal,"<<sourceCell.normal.x<<'/'<<sourceCell.normal.y<<",world_H_over_h,"<<worldHeight[0]/h<<'/'<<worldHeight[1]/h<<'/'<<worldHeight[2]/h<<'\n';
    check(methods[row]==HydroCsfCellCurvatureMethod::SegmentFit&&std::isfinite(k)&&load.heightFunctionCellEffectiveCurvatures()[row]==k,
        "mixed-chart opposed-steep local points bypassed qualification instead of reaching the existing physical method2 fallback");
    const Reference incidence(*mesh,phase);std::vector<double> color(layout.pressureCount());
    for(int p=0;p<layout.pressureCount();++p){const auto& q=phase.cellVolumes()[p];color[p]=q.liquid==0?0:q.gas==0?1:static_cast<double>(q.liquid/metrics.cellVolumes()[p]);}
    const auto gradient=incidence.gradient(color);const auto authoritative=pressure.apply(*mesh,HydroNativeValues(layout,std::vector<double>(layout.velocityCount()),color));unsigned forceRows=0;
    for(int u=0;u<layout.velocityCount();++u){close(load.colorGradient()[u],gradient[u],1e-13*std::max(1.,std::abs(gradient[u])),"mixed-chart original native physical incidence gradient");
        close(load.colorGradient()[u],authoritative[u],1e-13*std::max(1.,std::abs(authoritative[u])),"mixed-chart original authoritative pressure C");
        const int a=incidence.minus[u],b=incidence.plus[u];if((a!=row&&b!=row)||gradient[u]==0)continue;long double sum=0;unsigned count=0;
        for(int p:{a,b})if(p>=0&&std::isfinite(physical[p])){sum+=physical[p];++count;}
        check(count>0&&load.supportReports()[u].curvatureDefined,"mixed-chart source lost its actual physical force support");
        const auto expected=sigma*(sum/count)*authoritative[u];close(load.momentumRhsLoad()[u],expected,1e-12L*std::max(1e-30L,std::abs(expected)),"mixed-chart same-C adjacent physical curvature force");++forceRows;}
    check(forceRows>0,"mixed-chart lacks an actual source force row");const auto actualClock=load.timeLayer();
    check(actualClock.step==clock.step&&actualClock.begin==clock.begin&&actualClock.end==clock.end&&actualClock.position==clock.position,"mixed-chart changed its full source clock");meshIdentity.validate(*mesh);
    for(std::size_t ci=0;ci<stocks.size();++ci)check(publication.plicCells()[ci].inventory==stocks[ci],"mixed-chart changed independent closed-band phase stocks");
    check(publication.identity()==id&&isolatedFragmentSnapshot(publication)==saved,"mixed-chart changed actual source planes/measurement/inventory/resolution/owner");
    std::cout<<"NATIVE_SUPPORT_MIXED_HEIGHT_CHART_GUARD_PASS,native,"<<source<<",actual_world_H_over_h,"<<worldHeight[0]/h<<'/'<<worldHeight[1]/h<<'/'<<worldHeight[2]/h
        <<",independent_opposed_secants,"<<left<<'/'<<right<<",centered_slope,"<<(left+right)/2<<",normal,"<<sourceCell.normal.x<<'/'<<sourceCell.normal.y
        <<",actual_method,2,physical_kappa,"<<k<<",no_complete_HF_candidates,"<<noCompleteDonors<<",source_force_rows,"<<forceRows
        <<",raw_stocks_planes_clock_unchanged,1,no_fitted_value_or_curvature_magnitude_gate,1\n";
}
void mixedHeightPoints()
{
    namespace I=LAMP::Interface;constexpr double h=1./64,a=.30,b=1.3*h,sigma=.7;const Vector2D center{.5,.5};
    unsigned scenarios=0;
    for(unsigned complement:{0u,1u}){const double orientation=complement?-1.:1.;
    auto mesh=heightMesh(64);FaceSegmentGraph graph;graph.build(*mesh);const auto domain=NativeLeafDomain::build(*mesh,graph,false);
    std::vector<std::array<long double,2>> stocks(mesh->cellIndexNum());std::map<std::pair<int,int>,int> grid;long double liquidTotal=0;
    // Independent analytic circle primitive under an affine change of
    // coordinates. Exact pure boxes use their actual native capacity, rather
    // than acquire coordinate-transform roundoff as artificial minor stock.
    for(int ci=0;ci<mesh->cellIndexNum();++ci){const auto& cell=mesh->cells(mesh->cellIndexToLabel[ci]);
        const auto lo=mesh->verts(cell.cellToVert[0]).coord(),hi=mesh->verts(cell.cellToVert[2]).coord();
        const long double capacity=static_cast<long double>(hi.x-lo.x)*(hi.y-lo.y);
        const Vector2D qlo{(lo.x-center.x)/a,(lo.y-center.y)/b},qhi{(hi.x-center.x)/a,(hi.y-center.y)/b};
        const long double farX=std::max(std::abs(qlo.x),std::abs(qhi.x)),farY=std::max(std::abs(qlo.y),std::abs(qhi.y));
        const long double nearX=qlo.x>0?qlo.x:qhi.x<0?-qhi.x:0,nearY=qlo.y>0?qlo.y:qhi.y<0?-qhi.y:0;
        const long double liquid=farX*farX+farY*farY<=1?capacity:nearX*nearX+nearY*nearY>=1?0:
            static_cast<long double>(a)*b*circleCell(qlo,qhi,{0,0},1,false);
        check(liquid>=0&&liquid<=capacity,"mixed-height invalid setup: independent ellipse stock escaped actual native capacity");
        stocks[ci][complement]=liquid;stocks[ci][1-complement]=capacity-liquid;liquidTotal+=liquid;
        grid[{static_cast<int>(std::llround(lo.x/h)),static_cast<int>(std::llround(lo.y/h))}]=ci;}
    check(std::abs(liquidTotal-pi*a*b)<=1e-12L,"mixed-height invalid setup: independent whole ellipse area");
    const int source=grid.at({32,33});I::GeometricVofConfig config;config.localResolution=I::GeometricVofLocalResolution::NativeInterfaceSupport;
    const I::GeometricVofTimeLayer clock{49+complement,.002,.0022,LAMP::TimeLayer::Position::Midpoint};
    const auto baseline=I::GeometricVofPublication::prepareLocalPlicInventory(*mesh,domain,clock,stocks,{});
    const auto publication=I::GeometricVofPublication::prepareLocalPlicInventory(*mesh,domain,clock,stocks,config);
    const auto saved=isolatedFragmentSnapshot(publication);const auto id=publication.identity();const NativeMeshIdentity meshIdentity(*mesh);
    check(saved==isolatedFragmentSnapshot(baseline),"mixed-height law selection changed raw stock or actual PLIC geometry");
    const auto mask=publication.localPlicSubgridMask(*mesh);const auto& sourceCell=publication.plicCells()[source];
    check(!mask[source]&&sourceCell.resolved&&sourceCell.measured[0]>0&&sourceCell.measured[1]>0&&orientation*sourceCell.normal.y>0,
        "mixed-height invalid setup: named top ellipse source is not an owner-supported actual two-phase PLIC cut");
    const auto metrics=HydroPhysicalDivergence::build(*mesh,HydroIntegralDivergence::build(*mesh,HydroDofLayout::build(*mesh,publication.domain())));
    const auto phase=HydroPhaseGeometry::build(*mesh,metrics,publication);const auto pressure=HydroPressureLoad::build(metrics.divergence());const auto& layout=phase.layout();
    const int row=layout.pressureDof(source);LocalPlicEndpointReference independent{phase};const auto midpoint=independent.segmentMidpoint(sourceCell);
    const auto nativeX=independent.estimate(row,0),nativeY=independent.estimate(row,1);
    check(!nativeX.valid&&!nativeY.valid&&nativeY.opposite,"mixed-height invalid setup: the original complete charts did not reject the remote bottom branch");
    // These three independently certified short columns include one full
    // inside-phase box, the actual top cut and one full outside-phase box. Their
    // positions are separated by at least h and bracket the source tangent.
    std::array<Vector2D,3> points{};
    for(int column=-1;column<=1;++column){const auto& lower=publication.plicCells()[grid.at({32+column,32})];
        const auto& cut=publication.plicCells()[grid.at({32+column,33})];const auto& upper=publication.plicCells()[grid.at({32+column,34})];
        check(lower.measured[1-complement]==0&&upper.measured[complement]==0&&cut.resolved&&cut.measured[0]>0&&cut.measured[1]>0&&orientation*cut.normal.y>0,
            "mixed-height invalid setup: a near-top short column lacks complete ordered same-branch phase ends");
        const auto liquidHeight=(lower.measured[0]+cut.measured[0]+upper.measured[0])/h;
        const auto origin=complement?center.y+3*h:center.y;
        const auto physicalHeight=origin+orientation*liquidHeight;
        const auto insideHeight=(lower.measured[complement]+cut.measured[complement]+upper.measured[complement])/h;
        close(physicalHeight,center.y+insideHeight,3e-12L*h,"mixed-height phase-swapped short columns retain their correct physical origin");
        points[column+1]={(cut.rectangle.lower.x+cut.rectangle.upper.x)/2,static_cast<double>(physicalHeight)};
    }
    const Vector2D tangent{-sourceCell.normal.y,sourceCell.normal.x};
    const auto signedTangent=[&](Vector2D point){return (static_cast<long double>(point.x)-midpoint.x)*tangent.x+(static_cast<long double>(point.y)-midpoint.y)*tangent.y;};
    check(signedTangent(points[0])*signedTangent(points[2])<0,"mixed-height invalid setup: independent short-column points do not span the actual source tangent");
    for(unsigned first=0;first<3;++first)for(unsigned second=first+1;second<3;++second)
        check(std::hypot(points[first].x-points[second].x,points[first].y-points[second].y)>=.5*h,
            "mixed-height invalid setup: independent integral-height positions are not half-cell separated");
    const auto load=HydroBalancedCsfLoad::localPlicHeightFunctionCapillary(*mesh,pressure,phase,sigma,heightConfig());
    load.validateVofFor(*mesh,pressure,id,clock);const auto& methods=load.heightFunctionCellCurvatureMethods();const auto& physical=load.heightFunctionCellCurvatures();
    const double k=physical[row];check(methods[row]==HydroCsfCellCurvatureMethod::MixedHeightFit&&std::isfinite(k)&&orientation*k>0&&load.heightFunctionCellEffectiveCurvatures()[row]==k,
        "mixed-height did not use explicit physical method4 at the independently eligible thin-ellipse source");
    // Analytic graph curvature at the ACTUAL segment's tangential position,
    // independent of a fitted normal, polynomial coefficients or residual.
    const long double x=(static_cast<long double>(midpoint.x)-center.x)/a,ratio=static_cast<long double>(b)/a;
    const auto metricSquared=1-x*x+ratio*ratio*x*x,exact=(static_cast<long double>(b)/(a*a))/(metricSquared*std::sqrt(metricSquared));
    const auto error=std::abs(static_cast<long double>(k)-orientation*exact)/exact;
    check(error<.025L,"mixed-height actual thin-ellipse curvature exceeds the frozen 2.5 percent independent analytic gate");
    // A conservative superset of the original fixed5h donor box: every
    // aligned actual mixed native rectangle touching it must still fail both
    // complete charts, including all four possible end-extension rings.
    unsigned noFullHeightDonors=0;const auto half=2.5*h;
    for(int p=0;p<layout.pressureCount();++p){const auto& q=publication.plicCells()[layout.nativeCell(p)];
        if(q.measured[0]<=0||q.measured[1]<=0||static_cast<long double>(q.normal.x)*sourceCell.normal.x+static_cast<long double>(q.normal.y)*sourceCell.normal.y<=0)continue;
        if(q.rectangle.upper.x<=midpoint.x-half||q.rectangle.lower.x>=midpoint.x+half||q.rectangle.upper.y<=midpoint.y-half||q.rectangle.lower.y>=midpoint.y+half)continue;
        check(!independent.estimate(p,0).valid&&!independent.estimate(p,1).valid&&!independent.estimate(p,0,-8,8).valid&&!independent.estimate(p,1,-8,8).valid,
            "mixed-height invalid setup: the actual same-direction donor neighborhood contains an independently complete height chart");
        check(methods[p]!=HydroCsfCellCurvatureMethod::HeightFunction,"mixed-height invalid setup: an actual qualified full-HF donor would select the earlier constant fallback");++noFullHeightDonors;
    }
    check(noFullHeightDonors>=3,"mixed-height invalid setup: missing independently checked donor-neighborhood coverage");
    const Reference incidence(*mesh,phase);std::vector<double> color(layout.pressureCount());
    for(int p=0;p<layout.pressureCount();++p){const auto& q=phase.cellVolumes()[p];color[p]=q.liquid==0?0:q.gas==0?1:static_cast<double>(q.liquid/metrics.cellVolumes()[p]);}
    const auto gradient=incidence.gradient(color);const auto authoritative=pressure.apply(*mesh,HydroNativeValues(layout,std::vector<double>(layout.velocityCount()),color));unsigned loaded=0;
    for(int u=0;u<layout.velocityCount();++u){close(load.colorGradient()[u],gradient[u],1e-13*std::max(1.,std::abs(gradient[u])),"mixed-height retained independent actual face incidence C");
        close(load.colorGradient()[u],authoritative[u],1e-13*std::max(1.,std::abs(authoritative[u])),"mixed-height retained authoritative pressure C");
        const int first=incidence.minus[u],second=incidence.plus[u];if(first!=row&&second!=row)continue;
        if(gradient[u]==0){check(load.momentumRhsLoad()[u]==0,"mixed-height source has force without its actual color gradient");continue;}
        long double sum=0;unsigned count=0;for(int p:{first,second})if(p>=0&&std::isfinite(physical[p])){sum+=physical[p];++count;}
        check(count>0&&load.supportReports()[u].curvatureDefined,"mixed-height source lost its actual adjacent physical force support");
        const auto expected=sigma*(sum/count)*authoritative[u];close(load.momentumRhsLoad()[u],expected,1e-12L*std::max(1e-30L,std::abs(expected)),"mixed-height actual source force preserves sigma physical mean-kappa C");++loaded;
    }
    check(loaded>0,"mixed-height has no nonzero-gradient actual source force witness");const auto actualClock=load.timeLayer();
    check(actualClock.step==clock.step&&actualClock.begin==clock.begin&&actualClock.end==clock.end&&actualClock.position==clock.position,"mixed-height changed the complete force-owner clock");meshIdentity.validate(*mesh);
    for(std::size_t ci=0;ci<stocks.size();++ci)check(publication.plicCells()[ci].inventory==stocks[ci],"mixed-height changed independent raw phase stocks");
    check(publication.identity()==id&&isolatedFragmentSnapshot(publication)==saved&&isolatedFragmentSnapshot(baseline)==saved,"mixed-height changed actual PLIC planes, measurements, stocks or source identity");
    std::cout<<"NATIVE_SUPPORT_MIXED_HEIGHT_POINTS_SCENE_PASS,complement,"<<complement<<",native,"<<source<<",a,"<<a<<",b_over_h,"<<b/h<<",method,4,short_independent_points,3,fullHF_failed_candidates,"<<noFullHeightDonors
        <<",actual_segment_midpoint,"<<midpoint.x<<'/'<<midpoint.y<<",actual_kappa,"<<k<<",analytic_graph_kappa,"<<orientation*exact<<",relative_error,"<<error
        <<",frozen_relative_gate,.025,actual_source_force_rows,"<<loaded<<",raw_stocks_planes_clock_unchanged,1\n";
    ++scenarios;}
    check(scenarios==2,"mixed-height phase-swap coverage incomplete");
    std::cout<<"NATIVE_SUPPORT_MIXED_HEIGHT_POINTS_PASS,phase_swapped_scenes,2,frozen_relative_gate,.025\n";
}
void coincidentNativeFace()
{
    namespace I=LAMP::Interface;constexpr double h=1./64,radius=.125,sigma=.7;unsigned scenes=0,checkedFaces=0,rejectedLayers=0;
    for(unsigned complement:{0u,1u}){
        auto mesh=heightMesh(64);FaceSegmentGraph graph;graph.build(*mesh);const auto domain=NativeLeafDomain::build(*mesh,graph,false);
        const auto stadiumArea=[](Vector2D lo,Vector2D hi){
            const auto x0=std::max(lo.x,.25),x1=std::min(hi.x,.75),y0=std::max(lo.y,.375),y1=std::min(hi.y,.625);
            long double area=x1>x0&&y1>y0?static_cast<long double>(x1-x0)*(y1-y0):0;
            if(lo.x<.25)area+=circleCell(lo,{std::min(hi.x,.25),hi.y},{.25,.5},.125,false);
            if(hi.x>.75)area+=circleCell({std::max(lo.x,.75),lo.y},hi,{.75,.5},.125,false);return area;};
        std::map<std::pair<int,int>,int> grid;std::vector<std::array<long double,2>> stocks(mesh->cellIndexNum());long double analyticArea=0;
        for(int ci=0;ci<mesh->cellIndexNum();++ci){const auto& cell=mesh->cells(mesh->cellIndexToLabel[ci]);
            const auto lo=mesh->verts(cell.cellToVert[0]).coord(),hi=mesh->verts(cell.cellToVert[2]).coord();
            const auto capacity=static_cast<long double>(hi.x-lo.x)*(hi.y-lo.y),inside=stadiumArea(lo,hi);
            check(inside>=0&&inside<=capacity,"coincident-face analytic closed stadium escaped its actual native capacity");
            stocks[ci][complement]=inside;stocks[ci][1-complement]=capacity-inside;analyticArea+=inside;
            grid[{static_cast<int>(std::llround(lo.x/h)),static_cast<int>(std::llround(lo.y/h))}]=ci;}
        close(analyticArea,.5L*.25L+pi*radius*radius,1e-12L,"coincident-face independent complete stadium area");
        const int pureLiquid=grid.at({32,complement?40:39});const long double capacity=h*h,microGas=1e-20L*capacity;
        check(stocks[pureLiquid][0]==capacity&&stocks[pureLiquid][1]==0,"coincident-face micro carrier is not analytically pure liquid");
        stocks[pureLiquid]={capacity,microGas};
        const I::GeometricVofTimeLayer clock{52,.002,.0022,LAMP::TimeLayer::Position::Midpoint};I::GeometricVofConfig config;
        config.localResolution=I::GeometricVofLocalResolution::NativeInterfaceSupport;
        const auto legacy=I::GeometricVofPublication::prepareLocalPlicInventory(*mesh,domain,clock,stocks,{});
        const auto publication=I::GeometricVofPublication::prepareLocalPlicInventory(*mesh,domain,clock,stocks,config);
        const auto id=publication.identity(),legacyId=legacy.identity();const auto saved=isolatedFragmentSnapshot(publication);const NativeMeshIdentity meshIdentity(*mesh);
        check(saved==isolatedFragmentSnapshot(legacy)&&!publication.pinnedCorner(),"coincident-face capability changed raw stocks/planes or acquired a pin");
        const auto& carrier=publication.plicCells()[pureLiquid];
        check(carrier.inventory[0]==capacity&&carrier.inventory[1]==microGas&&carrier.measured[0]==capacity&&carrier.measured[1]==0&&!carrier.resolved,
            "coincident-face fixture lost its retained raw gas tail and unresolved measured-pure liquid geometry");
        const auto metrics=HydroPhysicalDivergence::build(*mesh,HydroIntegralDivergence::build(*mesh,HydroDofLayout::build(*mesh,publication.domain())));
        const auto phase=HydroPhaseGeometry::build(*mesh,metrics,publication);const auto pressure=HydroPressureLoad::build(metrics.divergence());const auto& layout=phase.layout();
        const auto load=HydroBalancedCsfLoad::localPlicHeightFunctionCapillary(*mesh,pressure,phase,sigma,heightConfig());
        load.validateVofFor(*mesh,pressure,id,clock);const Reference incidence(*mesh,phase);std::vector<double> color(layout.pressureCount());
        for(int p=0;p<layout.pressureCount();++p){const auto& q=phase.cellVolumes()[p];color[p]=q.liquid==0?0:q.gas==0?1:static_cast<double>(q.liquid/metrics.cellVolumes()[p]);}
        const auto independentC=incidence.gradient(color);const std::vector<double> zero(layout.velocityCount());
        const auto authoritativeC=pressure.apply(*mesh,HydroNativeValues(layout,zero,color));
        const auto measuredBox=[&](Vector2D lo,Vector2D hi){long double liquid=0,covered=0;
            for(const auto& q:publication.plicCells()){const Vector2D a{std::max(lo.x,q.rectangle.lower.x),std::max(lo.y,q.rectangle.lower.y)},b{std::min(hi.x,q.rectangle.upper.x),std::min(hi.y,q.rectangle.upper.y)};
                if(a.x>=b.x||a.y>=b.y)continue;covered+=static_cast<long double>(b.x-a.x)*(b.y-a.y);liquid+=LocalPlicEndpointReference::planeArea(q,a,b);}
            const auto want=static_cast<long double>(hi.x-lo.x)*(hi.y-lo.y);close(covered,want,1e-12L*want,"coincident-face independent query has complete actual rectangle coverage");return liquid;};
        for(int column:{32,16}){
            const int lower=grid.at({column,39}),upper=grid.at({column,40});int faceRow=-1;
            for(int u=0;u<layout.velocityCount();++u){const auto f=layout.domain().face(*mesh,layout.nativeFace(u));
                if((f.negativeCell==lower&&f.positiveCell==upper)||(f.negativeCell==upper&&f.positiveCell==lower))faceRow=u;}
            check(faceRow>=0,"coincident-face stadium lacks its named actual positive-area internal face");
            const auto face=layout.domain().face(*mesh,layout.nativeFace(faceRow));
            check(face.horizontal&&face.area>0&&face.center.y==.625,"coincident-face named frame is not the actual horizontal stadium top");
            for(int ci:{lower,upper}){const auto& q=publication.plicCells()[ci];const int p=layout.pressureDof(ci);
                check((q.measured[0]==0)!=(q.measured[1]==0),"coincident-face adjacent cell was fabricated as measured-mixed");
                check(!std::isfinite(load.heightFunctionCellCurvatures()[p])&&load.heightFunctionCellCurvatureMethods()[p]==HydroCsfCellCurvatureMethod::Unavailable,
                    "coincident-face assigned a fictitious own-cell curvature to a pure geometry");}
            // Keep the real tangential face position but place the normal
            // centre in its actually liquid-side cell. These nine strips then
            // align to the native rows, rather than mixing two pure branches
            // into artificial half-filled face-centred slices.
            const Vector2D frame{face.center.x,face.center.y+(publication.plicCells()[lower].measured[0]>0?-h/2:h/2)};
            std::array<long double,3> heights{};
            for(int transverse=-1;transverse<=1;++transverse){std::array<long double,9> fractions{};
                for(int line=-4;line<=4;++line){const Vector2D center{frame.x+transverse*h,frame.y+line*h};
                    fractions[line+4]=measuredBox(center-Vector2D{h/2,h/2},center+Vector2D{h/2,h/2})/(h*h);
                    heights[transverse+1]+=fractions[line+4]*h;}
                const int first=complement?8:0,last=complement?0:8;
                check(fractions[first]>=1-1e-6L&&fractions[last]<=1e-6L,"coincident-face independent complete height lacks the original opposite pure ends");
                long double previous=1;for(int k=0;k<9;++k){const auto value=fractions[complement?8-k:k];
                    check(value<=previous+1e-6L,"coincident-face independent complete height violates the original monotone topology gate");previous=value;}}
            check(heights[1]==5*h,"coincident-face independent central integral is not anchored to its actual native face");
            const auto expected=LocalPlicEndpointReference::heightCircumcircle(heights,h);const auto& report=load.supportReports()[faceRow];
            check(report.coincidentPhaseHeight&&report.curvatureDefined&&!report.boundary&&!report.underresolvedFragmentClosure&&report.supportCount==1&&report.totalWeight==1,
                "coincident-face has no independently valid face-owned height or masquerades as a masked/cell closure");
            const auto arithmetic=4096.L*std::numeric_limits<double>::epsilon()/h;
            close(load.curvatureAverages()[faceRow],expected,arithmetic,"coincident-face matches independent actual-plane height triangle curvature");
            check(independentC[faceRow]!=0,"coincident-face has no actual measured-color jump");
            close(load.colorGradient()[faceRow],independentC[faceRow],1e-13L*std::max(1.,std::abs(independentC[faceRow])),"coincident-face unchanged independent native incidence C");
            close(load.colorGradient()[faceRow],authoritativeC[faceRow],1e-13L*std::max(1.,std::abs(authoritativeC[faceRow])),"coincident-face unchanged authoritative pressure C");
            auto potential=color;for(auto& v:potential)v*=sigma*load.curvatureAverages()[faceRow];
            const auto constantPotential=pressure.apply(*mesh,HydroNativeValues(layout,zero,potential));
            close(load.momentumRhsLoad()[faceRow],constantPotential[faceRow],1e-12L*std::max(1e-30,std::abs(constantPotential[faceRow])),
                "coincident-face sigma kappa C equals the retained pressure-gradient potential");
            if(column==32)check(expected==0&&load.curvatureAverages()[faceRow]==0&&load.momentumRhsLoad()[faceRow]==0,
                "coincident-face straight plateau must have valid exact geometric zero curvature and force");
            else check((complement?-1.L:1.L)*expected>0&&std::abs(expected)>100*arithmetic&&load.momentumRhsLoad()[faceRow]!=0,
                "coincident-face round/straight join lacks its independent nonzero phase-signed height witness");
            std::cout<<"NATIVE_SUPPORT_COINCIDENT_FACE,complement,"<<complement<<",column,"<<column<<",velocity,"<<faceRow
                <<",height_points,"<<heights[0]<<'/'<<heights[1]<<'/'<<heights[2]<<",independent_triangle_kappa,"<<expected
                <<",actual_face_kappa,"<<load.curvatureAverages()[faceRow]<<",same_C_force,"<<load.momentumRhsLoad()[faceRow]<<",pure_own_kappa_NA,1\n";++checkedFaces;
        }
        check(load.curvatureReport().pureJumpHeightFaces>=2,"coincident-face work lacks original face-owned height provenance");
        const auto actualClock=load.timeLayer();check(actualClock.step==clock.step&&actualClock.begin==clock.begin&&actualClock.end==clock.end&&actualClock.position==clock.position,
            "coincident-face changed its complete actual source clock");meshIdentity.validate(*mesh);
        for(std::size_t ci=0;ci<stocks.size();++ci)check(publication.plicCells()[ci].inventory==stocks[ci],"coincident-face deleted or edited independent retained stocks");
        check(publication.identity()==id&&legacy.identity()==legacyId&&isolatedFragmentSnapshot(publication)==saved&&isolatedFragmentSnapshot(legacy)==saved,
            "coincident-face changed actual planes, measurements, resolution flags, raw stocks or owner");
        const auto legacyPhase=HydroPhaseGeometry::build(*mesh,metrics,legacy);bool legacyRejected=false;
        try{(void)HydroBalancedCsfLoad::localPlicHeightFunctionCapillary(*mesh,pressure,legacyPhase,sigma,heightConfig());}
        catch(const std::invalid_argument& e){legacyRejected=std::string(e.what()).find("active height-function CSF face without adjacent interface curvature")!=std::string::npos;}
        check(legacyRejected,"coincident-face silently enabled the new internal pure-interface capability for legacy VolumeOnly");++scenes;
        if(complement==0){
            auto badStocks=stocks;const int source=grid.at({32,39}),immediate=grid.at({32,38}),hole=grid.at({32,37});
            check(badStocks[source][0]==capacity&&badStocks[immediate][0]==capacity&&badStocks[hole][0]==capacity,
                "coincident-face thin-layer negative does not start in the actual complete pure-liquid column");
            badStocks[hole]={0,capacity};
            const auto badPublication=I::GeometricVofPublication::prepareLocalPlicInventory(*mesh,domain,clock,badStocks,config);
            const auto badId=badPublication.identity();const auto badSaved=isolatedFragmentSnapshot(badPublication);
            check(badPublication.plicCells()[source].measured[1]==0&&badPublication.plicCells()[immediate].measured[1]==0
                &&badPublication.plicCells()[hole].measured[0]==0&&badPublication.plicCells()[hole].measured[1]==capacity,
                "coincident-face negative lacks real L/G/L/G geometry with the immediate lower liquid carrier intact");
            const auto badPhase=HydroPhaseGeometry::build(*mesh,metrics,badPublication);bool layerRejected=false;
            try{(void)HydroBalancedCsfLoad::localPlicHeightFunctionCapillary(*mesh,pressure,badPhase,sigma,heightConfig());}
            catch(const std::invalid_argument& e){layerRejected=std::string(e.what()).find("coincident")!=std::string::npos;}
            check(layerRejected,"coincident-face native thin gas layer was averaged into a half-filled plateau, zero or unrelated remote graph");
            check(badPublication.identity()==badId&&isolatedFragmentSnapshot(badPublication)==badSaved,
                "coincident-face rejected layer changed its retained stocks, actual geometry or source owner");
            const auto t=badPublication.timeLayer();check(t.step==clock.step&&t.begin==clock.begin&&t.end==clock.end&&t.position==clock.position,
                "coincident-face rejected layer changed its complete source clock");++rejectedLayers;
        }
    }
    check(scenes==2&&checkedFaces==4&&rejectedLayers==1,"coincident-face phase/straight/join/layer-rejection coverage incomplete");
    std::cout<<"NATIVE_SUPPORT_COINCIDENT_FACE_PASS,scenes,2,independent_actual_faces,4,retained_unresolved_pure_liquid_gas_tail,1,native_thin_layer_rejected,1,legacy0_rejected,1,stocks_planes_owner_clock_unchanged,1\n";
}
void validFaceCurvature()
{
    namespace I=LAMP::Interface;constexpr long double h=1.L/32,epsilon=1.L/(1u<<24);constexpr double sigma=.7;
    auto mesh=heightMesh(32);FaceSegmentGraph graph;graph.build(*mesh);const auto domain=NativeLeafDomain::build(*mesh,graph,false);
    // A closed analytic capsule has a flat top at 21.25h. Its top cut is
    // below majority, whereas the complete liquid row below supports it.
    // The two added, retained gas-side tails therefore have no liquid-majority
    // witness; neither their effective zeros nor their normals are physical k.
    const auto capsuleArea=[](Vector2D lo,Vector2D hi){constexpr double a=.0625,b=.1640625;const Vector2D c{.5,.5};
        const double x0=std::max(lo.x,c.x-a),x1=std::min(hi.x,c.x+a),y0=std::max(lo.y,c.y-b),y1=std::min(hi.y,c.y+b);
        long double area=x1>x0&&y1>y0?static_cast<long double>(x1-x0)*(y1-y0):0;
        if(lo.x<c.x-a)area+=circleCell(lo,{std::min(hi.x,c.x-a),hi.y},{c.x-a,c.y},b,false);
        if(hi.x>c.x+a)area+=circleCell({std::max(lo.x,c.x+a),lo.y},hi,{c.x+a,c.y},b,false);return area;};
    std::map<std::pair<int,int>,int> grid;std::vector<std::array<long double,2>> stocks(mesh->cellIndexNum());
    for(int ci=0;ci<mesh->cellIndexNum();++ci){const auto& cell=mesh->cells(mesh->cellIndexToLabel[ci]);
        const auto lo=mesh->verts(cell.cellToVert[0]).coord(),hi=mesh->verts(cell.cellToVert[2]).coord();
        const auto capacity=static_cast<long double>(hi.x-lo.x)*(hi.y-lo.y),liquid=capsuleArea(lo,hi);
        check(liquid>=0&&liquid<=capacity,"valid-face analytic capsule stock escaped actual capacity");stocks[ci]={liquid,capacity-liquid};
        grid[{static_cast<int>(std::llround(lo.x/h)),static_cast<int>(std::llround(lo.y/h))}]=ci;}
    const int source=grid.at({16,21}),tail=grid.at({16,22}),outerTail=grid.at({16,23});
    check(stocks[source][0]==h*h/4&&stocks[tail][0]==0&&stocks[outerTail][0]==0,"valid-face independent capsule lacks its frozen quarter-cut and pure gas-side cells");
    stocks[tail]={epsilon*h*h,(1-epsilon)*h*h};stocks[outerTail]={epsilon*h*h/2,(1-epsilon/2)*h*h};
    const I::GeometricVofTimeLayer clock{48,.002,.0022,LAMP::TimeLayer::Position::Midpoint};I::GeometricVofConfig config;
    config.localResolution=I::GeometricVofLocalResolution::NativeInterfaceSupport;
    const auto baseline=I::GeometricVofPublication::prepareLocalPlicInventory(*mesh,domain,clock,stocks,{});
    const auto publication=I::GeometricVofPublication::prepareLocalPlicInventory(*mesh,domain,clock,stocks,config);
    const auto saved=isolatedFragmentSnapshot(publication);const auto id=publication.identity();const NativeMeshIdentity meshIdentity(*mesh);
    check(saved==isolatedFragmentSnapshot(baseline),"valid-face law selection changed raw inventories or actual PLIC geometry");
    const auto mask=publication.localPlicSubgridMask(*mesh);
    check(!mask[source]&&mask[tail]&&mask[outerTail],"valid-face fixture lacks its supported main source and two actual owner-masked tails");
    for(int ci:{source,tail,outerTail}){const auto& q=publication.plicCells()[ci];check(q.resolved&&q.measured[0]>0&&q.measured[1]>0,
        "valid-face fixture must use three real publicly reconstructed two-phase segments");}
    const auto metrics=HydroPhysicalDivergence::build(*mesh,HydroIntegralDivergence::build(*mesh,HydroDofLayout::build(*mesh,publication.domain())));
    const auto phase=HydroPhaseGeometry::build(*mesh,metrics,publication);const auto pressure=HydroPressureLoad::build(metrics.divergence());const auto& layout=phase.layout();
    const int sourceRow=layout.pressureDof(source),tailRow=layout.pressureDof(tail),outerTailRow=layout.pressureDof(outerTail);
    LocalPlicEndpointReference independent{phase};const auto column=independent.estimate(sourceRow,1);
    check(column.valid,"valid-face supported flat source lacks its independent complete original height chart");
    const auto load=HydroBalancedCsfLoad::localPlicHeightFunctionCapillary(*mesh,pressure,phase,sigma,heightConfig());
    load.validateVofFor(*mesh,pressure,id,clock);const auto& methods=load.heightFunctionCellCurvatureMethods();const auto& physical=load.heightFunctionCellCurvatures();
    const auto& effective=load.heightFunctionCellEffectiveCurvatures();const double k=physical[sourceRow];
    check(methods[sourceRow]==HydroCsfCellCurvatureMethod::HeightFunction&&std::isfinite(k)&&k!=0&&effective[sourceRow]==k,
        "valid-face main source must retain a nonzero physical actual HF curvature");
    close(k,independent.declaredCompleteCurvature(column,h),independent.arithmeticBound(sourceRow),
        "valid-face Native3 source circle curvature includes the complete actual-plane tail columns");
    for(int p:{tailRow,outerTailRow})check(methods[p]==HydroCsfCellCurvatureMethod::Unavailable&&!std::isfinite(physical[p])&&effective[p]==0,
        "valid-face actual subgrid tails must retain physicalNA and explicit effective zero");
    const Reference incidence(*mesh,phase);std::vector<double> color(layout.pressureCount());
    for(int p=0;p<layout.pressureCount();++p){const auto& q=phase.cellVolumes()[p];color[p]=q.liquid==0?0:q.gas==0?1:static_cast<double>(q.liquid/metrics.cellVolumes()[p]);}
    const auto independentGradient=incidence.gradient(color);const std::vector<double> zero(layout.velocityCount());
    const auto actualC=pressure.apply(*mesh,HydroNativeValues(layout,zero,color));auto potential=color;for(auto& value:potential)value*=sigma*k;
    const auto constantPotential=pressure.apply(*mesh,HydroNativeValues(layout,zero,potential));
    int soleFace=-1,maskedFace=-1;unsigned twoValid=0;double maximumTwoValidError=0;
    for(int u=0;u<layout.velocityCount();++u){const int a=incidence.minus[u],b=incidence.plus[u];
        close(load.colorGradient()[u],independentGradient[u],1e-13*std::max(1.,std::abs(independentGradient[u])),"valid-face unchanged independent physical incidence gradient");
        close(load.colorGradient()[u],actualC[u],1e-13*std::max(1.,std::abs(actualC[u])),"valid-face unchanged authoritative pressure C");
        if((a==sourceRow&&b==tailRow)||(a==tailRow&&b==sourceRow))soleFace=u;
        if((a==tailRow&&b==outerTailRow)||(a==outerTailRow&&b==tailRow))maskedFace=u;
        if(a<0||b<0||actualC[u]==0||!std::isfinite(physical[a])||!std::isfinite(physical[b]))continue;
        const auto expected=(static_cast<long double>(physical[a])+physical[b])/2;const auto& report=load.supportReports()[u];
        check(report.curvatureDefined&&report.supportCount==2&&report.totalWeight==2&&!report.underresolvedFragmentClosure,
            "valid-face two physical sides lost their unchanged positive arithmetic average support");
        close(load.curvatureAverages()[u],expected,1e-12L*std::max(1.L,std::abs(expected)),"valid-face two physical sides preserve their mean");
        const auto rhs=sigma*expected*actualC[u];close(load.momentumRhsLoad()[u],rhs,1e-12L*std::max(1e-30L,std::abs(rhs)),"valid-face two physical sides retain sigma mean-kappa C");
        maximumTwoValidError=std::max(maximumTwoValidError,static_cast<double>(std::abs(load.curvatureAverages()[u]-expected)));++twoValid;
    }
    check(soleFace>=0&&maskedFace>=0&&twoValid>0,"valid-face fixture lacks one of its three real native-face support cases");
    const auto& sole=load.supportReports()[soleFace];
    check(actualC[soleFace]!=0&&constantPotential[soleFace]!=0&&sole.curvatureDefined&&sole.underresolvedFragmentClosure&&sole.supportCount==1&&sole.totalWeight==1,
        "valid-face sole physical side lost its retained masked-neighbor provenance or nonzero force witness");
    close(load.curvatureAverages()[soleFace],k,1e-13*std::abs(k),"valid-face masked physicalNA must not halve the sole valid source curvature");
    close(load.momentumRhsLoad()[soleFace],constantPotential[soleFace],1e-12*std::abs(constantPotential[soleFace]),
        "valid-face sole curvature must preserve the same-C constant-curvature pressure potential");
    check(std::abs(load.momentumRhsLoad()[soleFace]-.5*constantPotential[soleFace])>.49*std::abs(constantPotential[soleFace]),
        "valid-face fixture cannot distinguish the rejected physical-kappa/2 face defect");
    const auto& masked=load.supportReports()[maskedFace];
    check(actualC[maskedFace]!=0&&masked.curvatureDefined&&masked.underresolvedFragmentClosure&&masked.supportCount==1&&masked.totalWeight==1
        &&load.curvatureAverages()[maskedFace]==0&&load.momentumRhsLoad()[maskedFace]==0,
        "valid-face both-masked nonzero color-gradient face lost its explicit typed zero closure");
    const auto actualClock=load.timeLayer();check(actualClock.step==clock.step&&actualClock.begin==clock.begin&&actualClock.end==clock.end&&actualClock.position==clock.position,
        "valid-face gather changed the complete actual force-owner clock");meshIdentity.validate(*mesh);
    for(std::size_t ci=0;ci<stocks.size();++ci)check(publication.plicCells()[ci].inventory==stocks[ci],"valid-face gather edited authoritative raw phase stocks");
    check(publication.identity()==id&&isolatedFragmentSnapshot(publication)==saved&&isolatedFragmentSnapshot(baseline)==saved,
        "valid-face gather changed actual geometry, measurement, inventory or owner identity");
    std::cout<<"NATIVE_SUPPORT_VALID_FACE_CURVATURE_PASS,source,"<<source<<",tail,"<<tail<<",outer_tail,"<<outerTail<<",source_physical_kappa,"<<k
        <<",sole_valid_face,"<<soleFace<<",sole_C,"<<actualC[soleFace]<<",actual_force,"<<load.momentumRhsLoad()[soleFace]<<",constant_potential_force,"<<constantPotential[soleFace]
        <<",two_valid_faces,"<<twoValid<<",maximum_two_valid_mean_error,"<<maximumTwoValidError<<",both_masked_face,"<<maskedFace
        <<",both_masked_C,"<<actualC[maskedFace]<<",both_masked_typed_force,0,raw_stocks_planes_clock_unchanged,1\n";
}
void maskedNormalSupport()
{
    namespace I=LAMP::Interface;constexpr long double h=1.L/32,epsilon=1.L/(1u<<24);constexpr int sx=16,sy=21;
    auto mesh=heightMesh(32);FaceSegmentGraph graph;graph.build(*mesh);const auto domain=NativeLeafDomain::build(*mesh,graph,false);
    std::map<std::pair<int,int>,int> grid;std::vector<std::array<long double,2>> stocks(mesh->cellIndexNum());
    for(int ci=0;ci<mesh->cellIndexNum();++ci){const auto& cell=mesh->cells(mesh->cellIndexToLabel[ci]);
        const auto lo=mesh->verts(cell.cellToVert[0]).coord(),hi=mesh->verts(cell.cellToVert[2]).coord();
        const auto capacity=static_cast<long double>(hi.x-lo.x)*(hi.y-lo.y),liquid=localPlicCapsuleArea(lo,hi);
        check(liquid>=0&&liquid<=capacity,"masked-normal independent closed-capsule stock escaped native capacity");
        grid[{static_cast<int>(std::llround(lo.x/h)),static_cast<int>(std::llround(lo.y/h))}]=ci;stocks[ci]={liquid,capacity-liquid};}
    const int source=grid.at({sx,sy}),lowerTail=grid.at({sx,sy+3}),upperTail=grid.at({sx,sy+4});
    for(const auto tail:{std::pair<int,long double>{lowerTail,epsilon},{upperTail,2*epsilon}}){
        check(stocks[tail.first][0]==0,"masked-normal retained tail fixture is not initially outside the analytic main interface");
        stocks[tail.first]={tail.second*h*h,(1-tail.second)*h*h};}
    const I::GeometricVofTimeLayer clock{47,.002,.0022,LAMP::TimeLayer::Position::Midpoint};
    const auto original=I::GeometricVofPublication::prepareLocalPlicInventory(*mesh,domain,clock,stocks,{});I::GeometricVofConfig config;
    config.localResolution=I::GeometricVofLocalResolution::NativeInterfaceSupport;
    const auto publication=I::GeometricVofPublication::prepareLocalPlicInventory(*mesh,domain,clock,stocks,config);
    const auto saved=isolatedFragmentSnapshot(publication);const auto id=publication.identity();const NativeMeshIdentity meshIdentity(*mesh);
    check(saved==isolatedFragmentSnapshot(original),"masked-normal law selection changed input stocks or actual PLIC planes");
    const auto mask=publication.localPlicSubgridMask(*mesh);const auto& main=publication.plicCells()[source];const auto& tail=publication.plicCells()[lowerTail];
    check(!mask[source]&&mask[lowerTail]&&mask[upperTail],"masked-normal fixture lacks a resolved main source and owner-classified subgrid tails");
    check(main.resolved&&main.measured[0]>0&&main.measured[1]>0&&main.normal.y>0
        &&tail.resolved&&tail.measured[0]>0&&tail.measured[1]>0&&tail.normal.y<0,
        "masked-normal fixture did not publicly reconstruct a real opposite-oriented two-phase tail");
    std::array<long double,3> rawHeight{};
    for(int column=-1;column<=1;++column){long double previous=1;
        for(int line=-4;line<=4;++line){const auto& q=stocks[grid.at({sx+column,sy+line})];const auto fraction=q[0]/(q[0]+q[1]);
            check(fraction<=previous+1e-6L,"masked-normal raw stock column violates the unchanged monotonic tolerance");previous=fraction;rawHeight[column+1]+=q[0]/h;
            if(line==-4)check(q[1]/(q[0]+q[1])<=1e-6L,"masked-normal raw stock column lacks its original liquid end");
            if(line==4)check(fraction<=1e-6L,"masked-normal raw stock column lacks its original gas end");}}
    const auto first=(rawHeight[2]-rawHeight[0])/(2*h),metric=std::hypot(1.L,first);
    const auto rawDifferentialMoment=-(rawHeight[2]-2*rawHeight[1]+rawHeight[0])/(h*h*metric*metric*metric);
    const auto rawCurvature=LocalPlicEndpointReference::heightCircumcircle(rawHeight,h);
    // The analytic flat main graph contributes no height difference. Both
    // retained tail inventories add to the central column, so their complete
    // differential HEIGHT moment is 2*(QL_lower+QL_upper)/h^3. Retain this
    // exact stock identity; Native3's declared physical curvature is instead
    // the independent three-position circumcircle above.
    const auto tailDifferentialMoment=2*(stocks[lowerTail][0]+stocks[upperTail][0])/(h*h*h);
    const auto metrics=HydroPhysicalDivergence::build(*mesh,HydroIntegralDivergence::build(*mesh,HydroDofLayout::build(*mesh,publication.domain())));
    const auto phase=HydroPhaseGeometry::build(*mesh,metrics,publication);const auto pressure=HydroPressureLoad::build(metrics.divergence());
    const int row=phase.layout().pressureDof(source);LocalPlicEndpointReference independent{phase};const auto prior=independent.estimate(row,1);
    check(prior.opposite&&!prior.valid,"masked-normal source does not independently expose the original opposite-normal veto");
    const auto arithmetic=independent.arithmeticBound(row);close(rawDifferentialMoment,tailDifferentialMoment,arithmetic,
        "masked-normal unchanged raw-stock differential height moment matches the independent retained-tail contribution");
    for(unsigned column=0;column<3;++column)close(prior.heights[column],rawHeight[column],9e-12L*h,
        "masked-normal actual full-plane column integral preserves every raw inventory within its original inverse-volume budget");
    const auto load=HydroBalancedCsfLoad::localPlicHeightFunctionCapillary(*mesh,pressure,phase,.7,heightConfig());
    load.validateVofFor(*mesh,pressure,id,clock);const auto k=load.heightFunctionCellCurvatures()[row];
    check(load.heightFunctionCellCurvatureMethods()[row]==HydroCsfCellCurvatureMethod::HeightFunction&&std::isfinite(k)
        &&load.heightFunctionCellEffectiveCurvatures()[row]==k,"masked-normal main source did not recover an actual physical HF estimate");
    close(k,independent.declaredCompleteCurvature(prior,h),arithmetic,
        "masked-normal Native3 recovered HF equals independent complete-plane height circumcircle");
    // Each of nine native cells keeps the original 1e-12 capacity inverse
    // budget. The three-column second difference has absolute weights 1,2,1;
    // this forward error bound introduces no new curvature magnitude limit.
    const auto rawBudget=36e-12L/h+arithmetic;check(tailDifferentialMoment>100*rawBudget,"masked-normal fixture cannot distinguish retained tails from deleted stocks");
    close(k,rawCurvature,rawBudget,"masked-normal recovered HF retains the full raw-stock curvature contribution");
    const auto actualClock=load.timeLayer();check(actualClock.step==clock.step&&actualClock.begin==clock.begin&&actualClock.end==clock.end&&actualClock.position==clock.position,
        "masked-normal HF changed the complete source clock");meshIdentity.validate(*mesh);
    for(std::size_t ci=0;ci<stocks.size();++ci)check(publication.plicCells()[ci].inventory==stocks[ci],"masked-normal HF deleted or altered original raw phase inventory");
    check(publication.identity()==id&&isolatedFragmentSnapshot(publication)==saved&&isolatedFragmentSnapshot(original)==saved,
        "masked-normal recovery changed either owner stock, actual geometry or identity");
    std::cout<<"NATIVE_SUPPORT_MASKED_NORMAL_PASS,source,"<<source<<",opposite_masked_tail,"<<lowerTail<<",retained_tail_pair,"<<lowerTail<<'/'<<upperTail
        <<",actual_method,1,full_raw_column_geometric_kappa,"<<rawCurvature<<",independent_full_plane_geometric_kappa,"<<independent.declaredCompleteCurvature(prior,h)<<",actual_HF_kappa,"<<k
        <<",unchanged_raw_differential_height_moment,"<<rawDifferentialMoment<<",retained_tail_analytic_differential_height_moment,"<<tailDifferentialMoment<<",original_volume_budget_forward_bound,"<<rawBudget
        <<",raw_stocks_planes_clock_unchanged,1,no_curvature_magnitude_gate,1\n";
}
void hfNeighborExtension()
{
    namespace I=LAMP::Interface;namespace Cut=LAMP::Geometry::PlanarCellCut;
    constexpr double h=1./64,radius=4.8*h,sigma=.7;const Vector2D center{.5,.5};
    auto mesh=heightMesh(64);FaceSegmentGraph graph;graph.build(*mesh);const auto domain=NativeLeafDomain::build(*mesh,graph,false);
    std::vector<std::array<long double,2>> stocks(mesh->cellIndexNum());int source=-1;
    for(int ci=0;ci<mesh->cellIndexNum();++ci){const auto& cell=mesh->cells(mesh->cellIndexToLabel[ci]);
        const auto lo=mesh->verts(cell.cellToVert[0]).coord(),hi=mesh->verts(cell.cellToVert[2]).coord();
        const auto capacity=static_cast<long double>(hi.x-lo.x)*(hi.y-lo.y),inside=circleCell(lo,hi,center,radius,false);
        stocks[ci]={inside,capacity-inside};const auto midpoint=(lo+hi)*.5;
        if(midpoint.x==center.x-3.5*h&&midpoint.y==center.y-3.5*h)source=ci;}
    check(source>=0,"HF neighbour extension lacks its independent diagonal circle source");
    I::GeometricVofConfig config;config.localResolution=I::GeometricVofLocalResolution::NativeInterfaceSupport;
    const I::GeometricVofTimeLayer clock{46,.002,.0022,LAMP::TimeLayer::Position::Midpoint};
    const auto publication=I::GeometricVofPublication::prepareLocalPlicInventory(*mesh,domain,clock,stocks,config);
    const auto saved=isolatedFragmentSnapshot(publication);const auto id=publication.identity();const NativeMeshIdentity meshIdentity(*mesh);
    check(!publication.localPlicSubgridMask(*mesh).at(source),"HF neighbour extension source was replaced by the subgrid force law");
    const auto metrics=HydroPhysicalDivergence::build(*mesh,HydroIntegralDivergence::build(*mesh,HydroDofLayout::build(*mesh,publication.domain())));
    const auto phase=HydroPhaseGeometry::build(*mesh,metrics,publication);const auto pressure=HydroPressureLoad::build(metrics.divergence());
    const int row=phase.layout().pressureDof(source);LocalPlicEndpointReference independent{phase};
    const auto x=independent.estimate(row,0),y=independent.estimate(row,1);
    check(x.opposite&&y.opposite&&!x.valid&&!y.valid,"HF neighbour extension source lacks independent opposite-branch failure of both native charts");
    const auto load=HydroBalancedCsfLoad::localPlicHeightFunctionCapillary(*mesh,pressure,phase,sigma,heightConfig());
    load.validateVofFor(*mesh,pressure,id,clock);const auto& methods=load.heightFunctionCellCurvatureMethods();const auto& curvatures=load.heightFunctionCellCurvatures();
    const double k=curvatures[row];check(methods[row]==HydroCsfCellCurvatureMethod::SegmentFit&&std::isfinite(k)&&k>0
        &&load.heightFunctionCellEffectiveCurvatures()[row]==k,"HF neighbour extension did not provide an actual physical method2 fallback");
    check(std::abs(k-1/radius)*radius<.35,"HF neighbour constant fallback exceeded the original independent circle accuracy gate");
    const auto segmentAt=[&](int native){const auto& q=publication.plicCells()[native];auto segment=Cut::interfaceSegment(q.rectangle,q.normal,q.offset);
        if(segment){const long double cx=static_cast<long double>(q.rectangle.lower.x)+(static_cast<long double>(q.rectangle.upper.x)-q.rectangle.lower.x)/2;
            const long double cy=static_cast<long double>(q.rectangle.lower.y)+(static_cast<long double>(q.rectangle.upper.y)-q.rectangle.lower.y)/2;
            for(auto& p:*segment){p.x+=cx;p.y+=cy;}}return segment;};
    const auto sourceSegment=segmentAt(source);check(bool(sourceSegment),"HF neighbour extension source has no actual public PLIC segment");
    const auto& sourceCell=publication.plicCells()[source];
    const long double ox=((*sourceSegment)[0].x+(*sourceSegment)[1].x)/2,oy=((*sourceSegment)[0].y+(*sourceSegment)[1].y)/2;
    const long double half=2.5L*std::max(sourceCell.rectangle.upper.x-sourceCell.rectangle.lower.x,sourceCell.rectangle.upper.y-sourceCell.rectangle.lower.y);
    // Independent line/box measure: enumerate the actual endpoints and all
    // four edge intersections, rather than replaying the production locator,
    // connected-component selection or interval-clipping implementation.
    const auto lengthInside=[&](int native){const auto segment=segmentAt(native);if(!segment)return 0.L;
        const auto a=(*segment)[0],b=(*segment)[1];const long double dx=b.x-a.x,dy=b.y-a.y;std::vector<long double> sites;
        const auto add=[&](long double t){if(t<0||t>1)return;const auto px=a.x+t*dx,py=a.y+t*dy;
            if(px>=ox-half&&px<=ox+half&&py>=oy-half&&py<=oy+half)sites.push_back(t);};
        const auto crossing=[&](long double t,bool vertical){if(t<0||t>1)return;
            // Its fixed coordinate is the exact chosen box edge; only the
            // independent tangential coordinate needs an interval check.
            const auto tangent=vertical?a.y+t*dy:a.x+t*dx,midpoint=vertical?oy:ox;
            if(tangent>=midpoint-half&&tangent<=midpoint+half)sites.push_back(t);};
        add(0);add(1);if(dx!=0){crossing((ox-half-a.x)/dx,true);crossing((ox+half-a.x)/dx,true);}
        if(dy!=0){crossing((oy-half-a.y)/dy,false);crossing((oy+half-a.y)/dy,false);}
        if(sites.size()<2)return 0.L;const auto ends=std::minmax_element(sites.begin(),sites.end());return std::hypot(dx,dy)*(*ends.second-*ends.first);};
    long double weight=0,residual=0,scale=0,maximumOperatorDifference=0;double minimum=std::numeric_limits<double>::infinity(),maximum=-minimum;
    unsigned donors=0,finiteCircleWitnesses=0;
    for(int donor=0;donor<phase.layout().pressureCount();++donor){if(methods[donor]!=HydroCsfCellCurvatureMethod::HeightFunction)continue;
        const int native=phase.layout().nativeCell(donor);const auto& q=publication.plicCells()[native];
        if(!q.resolved||q.measured[0]<=0||q.measured[1]<=0||static_cast<long double>(q.normal.x)*sourceCell.normal.x+static_cast<long double>(q.normal.y)*sourceCell.normal.y<=0)continue;
        const auto length=lengthInside(native);if(!(length>0))continue;const double value=curvatures[donor];check(std::isfinite(value),"qualified HF donor has no physical curvature");
        const int first=std::abs(q.normal.x)>=std::abs(q.normal.y)?0:1;const auto preferred=independent.estimate(donor,first);
        const auto alternate=preferred.valid?preferred:independent.estimate(donor,1-first);
        check(alternate.valid,"HF neighbour fixture donor has no independently qualified original native chart");
        const int selectedAxis=preferred.valid?first:1-first;const auto& rectangle=q.rectangle;
        const long double spacing=selectedAxis?rectangle.upper.x-rectangle.lower.x:rectangle.upper.y-rectangle.lower.y;
        const auto expected=independent.declaredCompleteCurvature(alternate,spacing),budget=independent.arithmeticBound(donor);
        close(value,expected,budget,"Native3 HF donor equals its independently integrated height triangle circumcircle");
        const auto difference=std::abs(expected-static_cast<long double>(alternate.curvature));maximumOperatorDifference=std::max(maximumOperatorDifference,difference);
        if(difference>16*budget){++finiteCircleWitnesses;
            std::cout<<"NATIVE_HEIGHT_CIRCLE_PUBLIC_WITNESS,native,"<<native<<",phase_law,3,independent_heights,"<<alternate.heights[0]<<'/'<<alternate.heights[1]<<'/'<<alternate.heights[2]
                <<",spacing,"<<spacing<<",legacy_differential,"<<alternate.curvature<<",independent_triangle_circle,"<<expected
                <<",actual_HF,"<<value<<",operator_difference,"<<difference<<",unchanged_arithmetic_bound,"<<budget<<'\n';}
        weight+=length;residual+=length*(static_cast<long double>(k)-value);scale+=length*(std::abs(k)+std::abs(value));
        minimum=std::min(minimum,value);maximum=std::max(maximum,value);++donors;
    }
    check(donors>=1&&weight>0&&scale>0,"HF neighbour extension lacks an independently located positive-length qualified donor");
    check(finiteCircleWitnesses>=1,"Native3 existing analytic circle lacks a public complete-HF witness distinguishing the retired differential operator");
    check(std::abs(residual)<=1e-12L*scale,"HF neighbour constant fallback violates the independent positive-length least-squares normal equation");
    check(k>=minimum-1e-12*std::max(1.,std::abs(minimum))&&k<=maximum+1e-12*std::max(1.,std::abs(maximum)),
        "HF neighbour constant fallback lies outside its immutable physical donor range");
    const auto actualClock=load.timeLayer();check(actualClock.step==clock.step&&actualClock.begin==clock.begin&&actualClock.end==clock.end&&actualClock.position==clock.position,
        "HF neighbour extension changed its complete source clock");meshIdentity.validate(*mesh);
    check(publication.identity()==id&&isolatedFragmentSnapshot(publication)==saved,"HF neighbour extension changed source stocks, actual planes or owner");
    std::cout<<"NATIVE_SUPPORT_HF_NEIGHBOR_EXTENSION_PASS,native,"<<source<<",qualified_HF_donors,"<<donors<<",positive_clip_length,"<<weight
        <<",constant_fit_kappa,"<<k<<",donor_min,"<<minimum<<",donor_max,"<<maximum<<",normal_equation_relative_residual,"<<std::abs(residual)/scale
        <<",original_circle_relative_error,"<<std::abs(k-1/radius)*radius<<",finite_height_circle_witnesses,"<<finiteCircleWitnesses
        <<",maximum_geometric_vs_differential_difference,"<<maximumOperatorDifference<<",raw_stocks_planes_clock_unchanged,1,method2_not_assumed_Kasa,1\n";
}
void circleFit()
{
    namespace I=LAMP::Interface;constexpr double h=1./64,radius=4.8*h,sigma=.7;const Vector2D center{.5,.5};
    auto mesh=heightMesh(64);FaceSegmentGraph graph;graph.build(*mesh);const auto domain=NativeLeafDomain::build(*mesh,graph,false);
    std::vector<std::array<long double,2>> stocks(mesh->cellIndexNum());int source=-1;
    for(int ci=0;ci<mesh->cellIndexNum();++ci){const auto& cell=mesh->cells(mesh->cellIndexToLabel[ci]);
        const auto lo=mesh->verts(cell.cellToVert[0]).coord(),hi=mesh->verts(cell.cellToVert[2]).coord();
        const auto capacity=static_cast<long double>(hi.x-lo.x)*(hi.y-lo.y),inside=circleCell(lo,hi,center,radius,false);
        check(inside>=0&&inside<=capacity,"Native2 narrow-circle independent stock outside actual native capacity");stocks[ci]={inside,capacity-inside};
        const auto midpoint=(lo+hi)*.5;if(midpoint.x==center.x-3.5*h&&midpoint.y==center.y-3.5*h)source=ci;}
    check(source>=0,"Native2 narrow-circle fixture lacks its named independent diagonal native source");
    I::GeometricVofConfig config;config.localResolution=I::GeometricVofLocalResolution::NativeInterfaceSupport;
    const I::GeometricVofTimeLayer clock{45,.002,.0022,LAMP::TimeLayer::Position::Midpoint};
    const auto publication=I::GeometricVofPublication::prepareLocalPlicInventory(*mesh,domain,clock,stocks,config);
    const auto saved=isolatedFragmentSnapshot(publication);const auto id=publication.identity();const NativeMeshIdentity meshIdentity(*mesh);
    const auto& cell=publication.plicCells()[source];
    check(cell.resolved&&cell.measured[0]>0&&cell.measured[1]>0&&cell.rectangle.upper.x-cell.rectangle.lower.x==h
        &&cell.rectangle.upper.y-cell.rectangle.lower.y==h&&!publication.localPlicSubgridMask(*mesh).at(source),
        "Native2 narrow-circle source is not an actual resolved fine PLIC cut eligible for physical fit fallback");
    const auto metrics=HydroPhysicalDivergence::build(*mesh,HydroIntegralDivergence::build(*mesh,HydroDofLayout::build(*mesh,publication.domain())));
    const auto phase=HydroPhaseGeometry::build(*mesh,metrics,publication);const auto pressure=HydroPressureLoad::build(metrics.divergence());
    const int row=phase.layout().pressureDof(source);LocalPlicEndpointReference independent{phase};
    const auto x=independent.estimate(row,0),y=independent.estimate(row,1);
    check(x.opposite&&y.opposite&&!x.valid&&!y.valid,"Native2 narrow-circle source does not independently reject both original HF charts on actual opposite branches");
    const auto load=HydroBalancedCsfLoad::localPlicHeightFunctionCapillary(*mesh,pressure,phase,sigma,heightConfig());
    load.validateVofFor(*mesh,pressure,id,clock);const auto actualClock=load.timeLayer();const auto k=load.heightFunctionCellCurvatures()[row];
    const double exact=1/radius,error=std::abs(k-exact)/exact;
    check(load.heightFunctionCellCurvatureMethods()[row]==HydroCsfCellCurvatureMethod::SegmentFit&&std::isfinite(k)&&k>0
        &&load.heightFunctionCellEffectiveCurvatures()[row]==k,"Native2 narrow-circle actual factory did not preserve a positive physical existing-fit curvature");
    check(error<.35,"Native2 narrow-circle actual fit exceeds the original frozen 35 percent independent analytic accuracy gate");
    check(load.curvatureReport().heightFunctionFitCells>0&&load.curvatureReport().heightFunctionFitSegments>=load.curvatureReport().heightFunctionFitCells,
        "Native2 narrow-circle fallback lacks its positive qualified-HF/geometric support work evidence");
    check(actualClock.step==clock.step&&actualClock.begin==clock.begin&&actualClock.end==clock.end&&actualClock.position==clock.position,
        "Native2 narrow-circle actual fit changed its complete source clock");meshIdentity.validate(*mesh);
    for(std::size_t ci=0;ci<stocks.size();++ci)check(publication.plicCells()[ci].inventory==stocks[ci],"Native2 narrow-circle fit changed independent raw phase stocks");
    check(publication.identity()==id&&isolatedFragmentSnapshot(publication)==saved,"Native2 narrow-circle fit changed actual source stocks, planes or owner");
    std::cout<<"NATIVE_SUPPORT_CIRCLE_FIT_PASS,n,64,native,"<<source<<",independent_opposite_HF_charts,2,actual_method,2,kappa,"<<k
        <<",fallback_kind,qualified_HF_constant_or_position_circle,analytic_kappa,"<<exact<<",relative_error,"<<error<<",frozen_relative_gate,.35,fit_cells,"<<load.curvatureReport().heightFunctionFitCells
        <<",fit_segments,"<<load.curvatureReport().heightFunctionFitSegments<<",raw_stocks_planes_clock_unchanged,1\n";
}
void refinementProbe()
{
    namespace I=LAMP::Interface;constexpr double radius=.32,sigma=.7;const Vector2D center{.487,.521};
    std::array<double,3> rms{},maximum{};unsigned grid=0;
    for(int n:{32,64,128}){
        auto mesh=heightMesh(n);FaceSegmentGraph graph;graph.build(*mesh);const auto domain=NativeLeafDomain::build(*mesh,graph,false);
        std::vector<std::array<long double,2>> stocks(mesh->cellIndexNum());
        for(int ci=0;ci<mesh->cellIndexNum();++ci){const auto& cell=mesh->cells(mesh->cellIndexToLabel[ci]);
            const auto lo=mesh->verts(cell.cellToVert[0]).coord(),hi=mesh->verts(cell.cellToVert[2]).coord();
            const auto capacity=static_cast<long double>(hi.x-lo.x)*(hi.y-lo.y),inside=circleCell(lo,hi,center,radius,false);
            stocks[ci]={inside,capacity-inside};}
        I::GeometricVofConfig config;config.localResolution=I::GeometricVofLocalResolution::NativeInterfaceSupport;
        const I::GeometricVofTimeLayer clock{43,.002,.0022,LAMP::TimeLayer::Position::Midpoint};
        const auto publication=I::GeometricVofPublication::prepareLocalPlicInventory(*mesh,domain,clock,stocks,config);
        const auto saved=isolatedFragmentSnapshot(publication);const auto id=publication.identity();
        const auto mask=publication.localPlicSubgridMask(*mesh);
        const auto metrics=HydroPhysicalDivergence::build(*mesh,HydroIntegralDivergence::build(*mesh,HydroDofLayout::build(*mesh,publication.domain())));
        const auto phase=HydroPhaseGeometry::build(*mesh,metrics,publication);const auto pressure=HydroPressureLoad::build(metrics.divergence());
        const auto load=HydroBalancedCsfLoad::localPlicHeightFunctionCapillary(*mesh,pressure,phase,sigma,heightConfig());
        load.validateVofFor(*mesh,pressure,id,clock);const auto& layout=phase.layout();
        struct Extreme {double error=-1,kappa=0;int row=-1,native=-1;};std::array<Extreme,3> extremes;
        long double squared=0;unsigned mixed=0,available=0,unavailable=0,hf=0,fit=0,nonpositive=0;
        for(int row=0;row<layout.pressureCount();++row){const int native=layout.nativeCell(row);const auto& q=publication.plicCells()[native];
            if(q.measured[0]==0||q.measured[1]==0)continue;++mixed;const auto k=load.heightFunctionCellCurvatures()[row];
            if(!std::isfinite(k)){++unavailable;continue;}const auto error=std::abs(k-1/radius)*radius;const auto method=load.heightFunctionCellCurvatureMethods()[row];
            ++available;nonpositive+=k<=0;squared+=error*error;const auto keep=[&](unsigned slot){if(error>extremes[slot].error)extremes[slot]={error,k,row,native};};keep(0);
            if(method==HydroCsfCellCurvatureMethod::HeightFunction){++hf;keep(1);}else if(method==HydroCsfCellCurvatureMethod::SegmentFit){++fit;keep(2);}}
        rms[grid]=available?std::sqrt(static_cast<double>(squared/available)):std::numeric_limits<double>::quiet_NaN();
        maximum[grid]=extremes[0].error;
        std::cout<<"NATIVE_SUPPORT_CIRCLE_REFINEMENT_PROBE,n,"<<n<<",complement,0,actual_mixed,"<<mixed<<",physical_available,"<<available
            <<",physical_unavailable,"<<unavailable<<",nonpositive_physical,"<<nonpositive<<",subgrid_masked,"<<std::count(mask.begin(),mask.end(),true)
            <<",HF_cells,"<<hf<<",fit_cells,"<<fit<<",relative_RMS,";if(available)std::cout<<rms[grid];else std::cout<<"NA";
        std::cout<<",relative_max,";if(available)std::cout<<maximum[grid];else std::cout<<"NA";
        std::cout<<",diagnostic_only_no_acceptance_gate,1\n";
        for(unsigned slot=0;slot<extremes.size();++slot){const auto& peak=extremes[slot];
            std::cout<<"NATIVE_SUPPORT_CIRCLE_REFINEMENT_EXTREMUM,n,"<<n<<",branch,"<<(slot==0?"all":slot==1?"HF":"fit")<<",available,"<<(peak.row>=0?1:0);
            if(peak.row>=0){const auto& q=publication.plicCells()[peak.native];std::cout<<",native,"<<peak.native<<",pressure_row,"<<peak.row
                <<",box,"<<q.rectangle.lower.x<<'/'<<q.rectangle.lower.y<<'/'<<q.rectangle.upper.x<<'/'<<q.rectangle.upper.y
                <<",normal,"<<q.normal.x<<'/'<<q.normal.y<<",plic_offset,"<<q.offset<<",kappa,"<<peak.kappa<<",analytic_kappa,"<<1/radius
                <<",relative_error,"<<peak.error<<",method,"<<static_cast<int>(load.heightFunctionCellCurvatureMethods()[peak.row]);}
            std::cout<<'\n';}
        check(publication.identity()==id&&isolatedFragmentSnapshot(publication)==saved,"circle refinement diagnostic mutated source inventories/planes/owner");
        ++grid;
    }
    for(unsigned i=1;i<grid;++i)std::cout<<"NATIVE_SUPPORT_CIRCLE_REFINEMENT_RATIO,coarse_n,"<<(i==1?32:64)<<",fine_n,"<<(i==1?64:128)
        <<",RMS_fine_over_coarse,"<<rms[i]/rms[i-1]<<",maximum_fine_over_coarse,"<<maximum[i]/maximum[i-1]
        <<",diagnostic_only_no_acceptance_gate,1\n";
}
void run()
{
    namespace I=LAMP::Interface;
    constexpr double radius=.32,sigma=.7;const Vector2D center{.487,.521};
    std::array<std::array<double,2>,2> maximum{},rms{};unsigned scenarios=0,curvatureCells=0,forceRows=0,fitCells=0;
    for(unsigned complement:{0u,1u})for(int n:{32,64}){
        const unsigned grid=n==32?0u:1u;const double orientation=complement?-1.:1.;auto mesh=heightMesh(n);
        FaceSegmentGraph graph;graph.build(*mesh);const auto domain=NativeLeafDomain::build(*mesh,graph,false);
        std::vector<std::array<long double,2>> stocks(mesh->cellIndexNum());long double circleVolume=0;
        for(int ci=0;ci<mesh->cellIndexNum();++ci){const auto& cell=mesh->cells(mesh->cellIndexToLabel[ci]);
            const auto lo=mesh->verts(cell.cellToVert[0]).coord(),hi=mesh->verts(cell.cellToVert[2]).coord();
            const auto capacity=static_cast<long double>(hi.x-lo.x)*(hi.y-lo.y),inside=circleCell(lo,hi,center,radius,false);
            check(inside>=0&&inside<=capacity,"Native2 independent circle stock outside its actual native capacity");
            stocks[ci][complement]=inside;stocks[ci][1-complement]=capacity-inside;circleVolume+=inside;}
        check(std::abs(circleVolume-pi*radius*radius)<=1e-12L,"Native2 independent complete analytic circle area");
        I::GeometricVofConfig config;config.localResolution=I::GeometricVofLocalResolution::NativeInterfaceSupport;
        const I::GeometricVofTimeLayer clock{43,.002,.0022,LAMP::TimeLayer::Position::Midpoint};
        const auto publication=I::GeometricVofPublication::prepareLocalPlicInventory(*mesh,domain,clock,stocks,config);
        const auto mask=publication.localPlicSubgridMask(*mesh);check(std::none_of(mask.begin(),mask.end(),[](bool q){return q;}),
            "Native2 resolved analytic circle was classified as a subgrid interface");
        const auto saved=isolatedFragmentSnapshot(publication);const auto id=publication.identity();const NativeMeshIdentity meshIdentity(*mesh);
        const auto metrics=HydroPhysicalDivergence::build(*mesh,HydroIntegralDivergence::build(*mesh,HydroDofLayout::build(*mesh,publication.domain())));
        const auto phase=HydroPhaseGeometry::build(*mesh,metrics,publication);const auto pressure=HydroPressureLoad::build(metrics.divergence());
        const auto load=HydroBalancedCsfLoad::localPlicHeightFunctionCapillary(*mesh,pressure,phase,sigma,heightConfig());
        load.validateVofFor(*mesh,pressure,id,clock);const auto actualClock=load.timeLayer();
        check(actualClock.step==clock.step&&actualClock.begin==clock.begin&&actualClock.end==clock.end&&actualClock.position==clock.position,
            "Native2 actual LocalPlic force changed its full source clock");
        const auto& layout=phase.layout();std::vector<double> measuredColor(layout.pressureCount());long double squared=0;unsigned mixed=0;
        for(int row=0;row<layout.pressureCount();++row){const int native=layout.nativeCell(row);const auto& q=publication.plicCells()[native];
            check(q.inventory==stocks[native],"Native2 circle estimator changed authoritative independent phase stocks");
            measuredColor[row]=q.measured[0]==0?0:q.measured[1]==0?1:static_cast<double>(q.measured[0]/metrics.cellVolumes()[row]);
            const double k=load.heightFunctionCellCurvatures()[row];if(q.measured[0]==0||q.measured[1]==0){
                check(!std::isfinite(k),"Native2 raw-pure circle leaf acquired a physical curvature");continue;}
            const auto method=load.heightFunctionCellCurvatureMethods()[row];
            check(q.resolved&&std::isfinite(k)&&orientation*k>0&&(method==HydroCsfCellCurvatureMethod::HeightFunction||method==HydroCsfCellCurvatureMethod::SegmentFit),
                "Native2 resolved circle lost physical curvature or was replaced by an unresolved closure");
            check(load.heightFunctionCellEffectiveCurvatures()[row]==k,"Native2 resolved circle effective coefficient differs from its physical estimate");
            const double error=std::abs(k-orientation/radius)*radius;maximum[complement][grid]=std::max(maximum[complement][grid],error);
            squared+=error*error;++mixed;fitCells+=method==HydroCsfCellCurvatureMethod::SegmentFit;}
        check(mixed>0&&load.curvatureReport().heightFunctionCells+load.curvatureReport().heightFunctionFitCells==mixed
            &&load.curvatureReport().unresolvedIsolatedFragmentCells==0,"Native2 analytic circle estimator accounting includes an unresolved closure");
        rms[complement][grid]=std::sqrt(static_cast<double>(squared/mixed));
        check(maximum[complement][grid]<(n==32?.08:.025),"Native2 independent analytic LocalPlic circle exceeds frozen maximum relative curvature error");
        const Reference incidence(*mesh,phase);const auto gradient=incidence.gradient(measuredColor);unsigned supported=0;
        for(int u=0;u<layout.velocityCount();++u){close(load.colorGradient()[u],gradient[u],1e-12*std::max(1.,std::abs(gradient[u])),
                "Native2 force retained independent actual native C measured-color algebra");
            const auto& support=load.supportReports()[u];long double sum=0;unsigned count=0;
            for(int p:{incidence.minus[u],incidence.plus[u]})if(p>=0&&std::isfinite(load.heightFunctionCellCurvatures()[p])){sum+=load.heightFunctionCellCurvatures()[p];++count;}
            if(support.curvatureDefined){check(count>0&&!support.underresolvedFragmentClosure,"Native2 resolved circle force borrowed an unresolved curvature");
                const auto expected=sigma*(sum/count)*gradient[u];close(load.curvatureAverages()[u],sum/count,1e-12*std::max(1.L,std::abs(sum/count)),
                    "Native2 actual adjacent physical curvature mean");close(load.momentumRhsLoad()[u],expected,1e-12*std::max(1.L,std::abs(expected)),
                    "Native2 unchanged actual sigma kappa C measured-color force");++supported;}
            else check(load.momentumRhsLoad()[u]==0,"Native2 circle used an unavailable force coefficient");}
        check(supported>0,"Native2 resolved circle has no actual force rows");meshIdentity.validate(*mesh);
        check(publication.identity()==id&&isolatedFragmentSnapshot(publication)==saved,"Native2 estimator changed actual source stocks, planes or owner");
        curvatureCells+=mixed;forceRows+=supported;++scenarios;
        std::cout<<"NATIVE_SUPPORT_HEIGHT_CHART_CIRCLE,n,"<<n<<",complement,"<<complement<<",actual_mixed,"<<mixed
            <<",HF_cells,"<<load.curvatureReport().heightFunctionCells<<",fit_cells,"<<load.curvatureReport().heightFunctionFitCells
            <<",relative_RMS,"<<rms[complement][grid]<<",relative_max,"<<maximum[complement][grid]<<",PASS\n";
    }
    for(unsigned complement:{0u,1u}){check(rms[complement][1]<.65*rms[complement][0]&&maximum[complement][1]<.65*maximum[complement][0],
        "Native2 independent actual LocalPlic circle curvature errors do not meet frozen refinement ratios");
        std::cout<<"NATIVE_SUPPORT_HEIGHT_CHART_REFINEMENT,complement,"<<complement<<",RMS_fine_over_coarse,"<<rms[complement][1]/rms[complement][0]
            <<",maximum_fine_over_coarse,"<<maximum[complement][1]/maximum[complement][0]<<",PASS\n";}
    check(scenarios==4&&curvatureCells>0&&forceRows>0,"Native2 actual LocalPlic circle phase/refinement coverage incomplete");
    std::cout<<"NATIVE_SUPPORT_HEIGHT_CHART_PASS,scenarios,"<<scenarios<<",physical_curvature_cells,"<<curvatureCells<<",actual_force_rows,"<<forceRows
        <<",adopted_fit_cells,"<<fitCells<<",independent_analytic_circle,1,raw_stocks_planes_clock_unchanged,1\n";
}
}

int main(int argc,char** argv)
{
    try{
        std::cout<<std::setprecision(17);
        if(argc==2&&std::string(argv[1])=="--native-support-coincident-face-only"){NativeSupportHeightChartContract::coincidentNativeFace();return 0;}
        if(argc==2&&std::string(argv[1])=="--native-support-mixed-height-chart-guard-only"){NativeSupportHeightChartContract::mixedHeightChartGuard();return 0;}
        if(argc==2&&std::string(argv[1])=="--native-support-mixed-height-points-only"){NativeSupportHeightChartContract::mixedHeightPoints();return 0;}
        if(argc==2&&std::string(argv[1])=="--native-support-valid-face-curvature-only"){NativeSupportHeightChartContract::validFaceCurvature();return 0;}
        if(argc==2&&std::string(argv[1])=="--native-support-masked-normal-only"){NativeSupportHeightChartContract::maskedNormalSupport();return 0;}
        if(argc==2&&std::string(argv[1])=="--native-support-hf-neighbor-extension-only"){NativeSupportHeightChartContract::hfNeighborExtension();return 0;}
        if(argc==2&&std::string(argv[1])=="--native-support-circle-fit-only"){NativeSupportHeightChartContract::circleFit();return 0;}
        if(argc==2&&std::string(argv[1])=="--native-support-circle-refinement-probe-only"){NativeSupportHeightChartContract::refinementProbe();return 0;}
        if(argc==2&&std::string(argv[1])=="--native-support-height-chart-regression-only")
            throw std::invalid_argument("NativeSupportHeightChartContract::regression is retired with law2; this option is unavailable");
        if(argc==2&&std::string(argv[1])=="--native-support-height-chart-only"){NativeSupportHeightChartContract::run();return 0;}
        if(argc==2&&std::string(argv[1])=="--local-plic-workspace-only"){HfWorkspaceContract::run();return 0;}
        if(argc==2&&std::string(argv[1])=="--hf-dyadic-width-only"){DyadicWidthContract::run();return 0;}
        if(argc==2&&std::string(argv[1])=="--pinned-corner-hf-only"){PinnedCornerHfContract::run();return 0;}
        if(argc==2&&std::string(argv[1])=="--south-even90-only"){F7SouthEven90::run();return 0;}
        if(argc==2&&std::string(argv[1])=="--local-plic-vertex-star-fragment-only"){F25VertexStarFragmentContract::run();return 0;}
        if(argc>1&&(std::string(argv[1])=="--normal-difference-only"||std::string(argv[1])=="--normal-local-pressure-only"))
            throw std::invalid_argument("NormalDifference capillary research is retired; this option is unavailable");
        if(argc==2&&std::string(argv[1])=="--local-plic-two-cell-fragment-only"){localPlicTwoCellFragmentClosure();return 0;}
        if(argc==2&&std::string(argv[1])=="--local-plic-isolated-fragment-only"){localPlicIsolatedFragmentClosure();return 0;}
        if(argc==2&&std::string(argv[1])=="--local-plic-rz-regrid-geometry-only"){rzRegridHeightGeometryProbe();return 0;}
        if(argc==2&&std::string(argv[1])=="--local-plic-rz-cylinder-only"){rzCylinderHeightScene();rzPhysicalWallNegatives();return 0;}
        if(argc==2&&std::string(argv[1])=="--local-plic-rz-only"){localPlicRzHeightFunction();return 0;}
        if(argc==2&&std::string(argv[1])=="--local-plic-only"){localPlicHeightFunction();localPlicGradedHeightFunction();return 0;}
        if(argc==2&&std::string(argv[1])=="--height-function-only"){heightFunctionOnly();return 0;}
        if(argc==2&&std::string(argv[1])=="--height-increment"){heightFunctionIncrementOnly();return 0;}
        if(argc==2&&std::string(argv[1])=="--initial-acceleration"){initialDropAcceleration();return 0;}
        if(argc==2&&std::string(argv[1])=="--geometric-only"){geometricOnly();return 0;}
        if(argc==2&&std::string(argv[1])=="--adjoint-only"){adjointOnly();return 0;}
        if(argc==2&&std::string(argv[1])=="--projected-adjoint-only"){projectedAdjointOnly();return 0;}
        if(argc==2&&std::string(argv[1])=="--projected-macro-only"){projectedAdjointOnly(true);return 0;}
        std::cout<<"BF02_CELL_HEADER,scene,shape,grid,selection,rho,mu,sigma,kappa,pressure,native_cell,component,volume,epsilon,color,sharp_alpha,p,p_expected,component_color_mean\n";
        std::cout<<"BF02_FACE_HEADER,scene,shape,grid,selection,rho,mu,sigma,kappa,velocity,negative,positive,area,free,axis,Cc,S,Cp,Cp_minus_S,Au,reaction,reaction_expected\n";
        const bool domainsOnly=argc==2&&std::string(argv[1])=="--domains";
        check(argc==1||domainsOnly,"unknown BF02 contract option");
        if(domainsOnly)domains();else{core();provenance();geometricOnly();adjointOnly();projectedAdjointOnly();projectedAdjointOnly(true);heightSmallContract();localPlicHeightFunction();}
        if(!domainsOnly)localPlicGradedHeightFunction();
        if(domainsOnly)check(wallFaces>0&&axisFaces>0,"selected/pure/axis witness missing");
        else check(coarseFaces>0&&axisFaces>0&&wallFaces>0&&negativeControls>0&&guards>0,"core/provenance witness missing");
        std::cout<<"BF02_SUMMARY,"<<(domainsOnly?"domains":"core")<<','<<scenes<<','<<solves<<','<<guards<<','<<negativeControls<<','<<coarseFaces<<','<<axisFaces<<','
            <<wallFaces<<','<<peakLoad<<','<<peakVelocity<<','<<peakPressure<<','<<peakBalance<<','<<peakWall<<','<<peakContinuity<<",PASS\n";return 0;
    }catch(const std::exception& e){std::cerr<<"BF02_CSF_FAIL,"<<e.what()<<'\n';return 1;}
}
