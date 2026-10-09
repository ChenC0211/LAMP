#include "Modules/Hydrodynamics/CchBalancedCsfLoad.h"
#include "Modules/Hydrodynamics/CchMomentumGeometry.h"
#include "Modules/Hydrodynamics/CchPhaseViscosity.h"
#include "Modules/Hydrodynamics/CchSinglePhaseSaddle.h"
#include "interface_first_jet_fixture.h"
#include <cstdint>
#include <map>
#include <numeric>

namespace
{
using namespace Quantum::Hydrodynamics;
using Matrix=CchSaddleVelocityBlock::Matrix;
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
CchPhysicalDivergence metricsFor(const QuadtreeMesh& mesh,bool rz,int selection=0)
{
    FaceSegmentGraph graph;graph.build(mesh,[selection](const QuadtreeMesh&,const Cell& cell){const auto p=cell.centerCoord();
        return selection==0 || (selection==1?!(p.x>.375&&p.x<.625&&p.y>.375&&p.y<.625):
            selection==2?p.x<.375||p.x>.625:p.x>.25);});
    return CchPhysicalDivergence::build(mesh,CchIntegralDivergence::build(mesh,
        CchDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,rz))));
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
CchRegularizedColor forceColor(const QuadtreeMesh& mesh,const CchPhaseGeometry& phase,double epsilon=.375)
{CchRegularizedColorConfig cfg;cfg.epsilon=epsilon;return CchRegularizedColor::build(mesh,phase,cfg);}

// Independent reference uses native incidence/physical edge coordinates,
// including every actual AMR subface. It never reads C/D coefficients.
struct Reference
{
    std::vector<int> minus,plus,component;std::vector<double> area,volume;int components=0;
    Reference(const QuadtreeMesh& mesh,const CchPhaseGeometry& phase)
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
void algebra(const QuadtreeMesh& mesh,const CchPressureLoad& pressure,const CchRegularizedColor& color,
    const Reference& ref,double sigma,double kappa)
{
    const auto load=CchBalancedCsfLoad::constantCurvature(mesh,pressure,color,sigma,kappa);
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
    const auto cp=pressure.apply(mesh,CchNativeValues(color.phaseGeometry().layout(),std::vector<double>(gradient.size()),potential));
    error=0;for(std::size_t u=0;u<gradient.size();++u)error=std::hypot(error,cp[u]-load.momentumRhsLoad()[u]);
    check(relative(error,scale)<=1e-12,"all-row nongauge constant-curvature identity");
    // Virtual work uses independent physical cell/face incidence.
    std::vector<long double> div(color.colors().size());long double work=0,workScale=0;
    for(std::size_t u=0;u<gradient.size();++u){const double velocity=std::sin(.31*(u+1));const auto term=velocity*load.momentumRhsLoad()[u];
        work+=term;workScale+=std::abs(term);if(ref.minus[u]>=0)div[ref.minus[u]]+=ref.area[u]*velocity;if(ref.plus[u]>=0)div[ref.plus[u]]-=ref.area[u]*velocity;}
    for(std::size_t p=0;p<potential.size();++p){const auto term=potential[p]*div[p];work+=term;workScale+=std::abs(term);}
    check(relative(static_cast<double>(work),static_cast<double>(workScale))<=1e-12,"C=-D_I transpose virtual work");
}
void solveStatic(const QuadtreeMesh& mesh,const CchPhaseGeometry& phase,const CchRegularizedColor& color,
    const CchPressureLoad& pressure,const CchMomentumGeometry& dual,const Matrix& viscosity,
    const Reference& ref,const char* tag,int grid,int selection,double rho,double mu,bool wrongSign=false,double positiveCurvature=5)
{
    ++scenes;const auto& l=phase.layout();const auto mass=CchMomentumMass::build(mesh,dual,phase,rho,1);
    Matrix a=viscosity;if(mu==0)check(a.nonZeros()==0,"zero viscosity has artificial floor");
    for(int u=0;u<l.velocityCount();++u)a.coeffRef(u,u)+=mass.masses()[u]/.125;
    const auto saddle=CchSinglePhaseSaddle::buildClosed(pressure,CchSaddleVelocityBlock::build(l,a));
    check(saddle.pressureGaugeCount()==ref.components,"independent connected-component gauge");
    std::vector<bool> free(l.velocityCount(),false);for(int u:saddle.freeVelocityDofs())free[u]=true;
    for(int u=0;u<l.velocityCount();++u)check(free[u]==(ref.minus[u]>=0&&ref.plus[u]>=0&&ref.area[u]>0),"native free-row classification");
    const auto means=ref.meanByCell(color.colors()),meanGradient=ref.gradient(means);
    for(const auto parameters:std::array<std::array<double,2>,3>{{{.7,-2.5},{2.5,0},{.5,positiveCurvature}}}){
        const double sigma=parameters[0],kappa=parameters[1],beta=sigma*kappa;
        algebra(mesh,pressure,color,ref,sigma,kappa);
        const auto load=CchBalancedCsfLoad::constantCurvature(mesh,pressure,color,sigma,kappa);
        const auto result=saddle.solve(mesh,CchSaddleRhs(l,load.momentumRhsLoad(),std::vector<double>(l.velocityCount())));++solves;
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
                const auto wrong=saddle.solve(mesh,CchSaddleRhs(l,badRhs,std::vector<double>(l.velocityCount())));++solves;
                double defect=0;for(int p=0;p<l.pressureCount();++p)defect=std::max(defect,std::abs(wrong.values().pressure()[p]-expected[p]));
                check(defect/std::abs(beta)>1e-10,"actual wrong-RHS sign negative control invisible");
                if(ref.component[i]==ref.component[j]&&*high-*low>1e-6)
                    check(wrong.values().pressure()[i]<wrong.values().pressure()[j],"actual wrong RHS failed to reverse liquid pressure");++negativeControls;}
        }
        std::cout<<"BF02_STATIC,"<<scenes<<','<<tag<<','<<rho<<','<<mu<<','<<sigma<<','<<kappa<<','<<velocity<<','<<pError<<','<<balance<<','<<wall<<','
            <<result.report().linearRefinements<<",PASS\n";
    }
}
void inspect(const QuadtreeMesh& mesh,const CchPhaseGeometry& phase,const CchRegularizedColor& color,
    const char* tag,int grid=0,int selection=0,bool densityMatrix=true,double positiveCurvature=5)
{
    check(phase.publication().state().primary().phaseVolumes.empty(),"CSF introduced conservative Interface inventory");
    const auto sharpBefore=phase.liquidFractions();const auto pressure=CchPressureLoad::build(phase.metrics().divergence());
    const Reference ref(mesh,phase);const auto dual=CchMomentumGeometry::build(mesh,phase.metrics());
    const auto quadrature=HermiteGalsPhaseQuadrature::build(mesh,phase.publication());
    algebra(mesh,pressure,color,ref,0,std::numeric_limits<double>::max());
    algebra(mesh,pressure,color,ref,std::numeric_limits<double>::max(),0);
    for(double mu:{0.,1e-3}){const auto viscosity=CchPhaseViscosity::build(mesh,phase,quadrature,mu,.01*mu);
        for(double rho:{1.,1e3,1e6}){if(!densityMatrix&&rho!=1e3)continue;
            solveStatic(mesh,phase,color,pressure,dual,viscosity.viscosity().matrix(),ref,tag,grid,selection,rho,mu,
                densityMatrix&&grid==0&&rho==1&&mu==0,positiveCurvature);}}
    check(phase.liquidFractions()==sharpBefore,"CSF changed sharp-material fractions");
}
void core()
{
    for(bool rz:{false,true})for(int amr:{0,1,2}){auto mesh=makeMesh(4,amr);Case shape;shape.pureFirstJet=true;shape.radial=rz;shape.radius=rz?.34:.19;
        const auto phase=CchPhaseGeometry::build(*mesh,metricsFor(*mesh,rz),published(*mesh,shape));const auto color=forceColor(*mesh,phase);
        // The constant-curvature drop control uses its independent analytic
        // physical curvature: 1/R in XY and 2/R for the RZ sphere.
        inspect(*mesh,phase,color,rz?"sphere":"circle",amr,0,true,(rz?2.:1.)/shape.radius);}
    std::uint32_t random=170120u;
    for(bool rz:{false,true})for(int amr:{0,1,2}){random=1664525u*random+1013904223u;
        const double offset=.31+.38*static_cast<double>(random)/std::numeric_limits<std::uint32_t>::max();
        auto mesh=makeMesh(4,amr);Case shape;shape.pureFirstJet=true;shape.radial=rz;shape.planar=true;shape.planeOffset=offset;
        const auto met=metricsFor(*mesh,rz);const auto phase=CchPhaseGeometry::build(*mesh,met,published(*mesh,shape));const auto color=forceColor(*mesh,phase);
        const auto scaledPhase=CchPhaseGeometry::build(*mesh,met,published(*mesh,shape,7));const auto scaled=forceColor(*mesh,scaledPhase);
        for(std::size_t p=0;p<color.colors().size();++p)close(color.colors()[p],scaled.colors()[p],5e-9,"positive phi scaling changed CSF color");
        const auto pressure=CchPressureLoad::build(met.divergence());const auto original=CchBalancedCsfLoad::constantCurvature(*mesh,pressure,color,.5,5);
        const auto altered=CchBalancedCsfLoad::constantCurvature(*mesh,pressure,scaled,.5,5);
        const Reference ref(*mesh,phase);for(std::size_t u=0;u<ref.area.size();++u)
            close(original.momentumRhsLoad()[u],altered.momentumRhsLoad()[u],2.5*1e-8*ref.area[u],"scaled true-area CSF gradient");
        inspect(*mesh,phase,color,rz?"axial-plane":"plane",amr,0,false);
        std::cout<<"BF02_PLANE_OFFSET,"<<rz<<','<<amr<<','<<offset<<",PASS\n";}
}
void domains()
{
    for(bool rz:{false,true})for(int selection:{1,2}){auto mesh=makeMesh(8);Case shape;shape.pureFirstJet=true;shape.radial=rz;shape.planar=true;shape.planeOffset=.43;
        if(rz&&selection==2){shape.planar=false;shape.radius=.34;}
        const auto phase=CchPhaseGeometry::build(*mesh,metricsFor(*mesh,rz,selection),published(*mesh,shape));
        const auto color=forceColor(*mesh,phase,.1875);inspect(*mesh,phase,color,selection==1?"hole":rz?"disconnected-sphere":"disconnected",0,selection,false,
            rz&&selection==2?2/shape.radius:5);
        const Reference ref(*mesh,phase);check(ref.components==(selection==1?1:2),"selected connected-component witness");
        if(selection==2){const auto mean=ref.meanByCell(color.colors());
            check(*std::max_element(mean.begin(),mean.end())-*std::min_element(mean.begin(),mean.end())>1e-4,"disconnected gauge lacked distinct component means");}}
    for(bool rz:{false,true})for(double sign:{-1.,1.}){auto mesh=makeMesh(4,1);
        const auto phase=CchPhaseGeometry::build(*mesh,metricsFor(*mesh,rz),purePublication(*mesh,rz,sign));const auto color=forceColor(*mesh,phase);
        inspect(*mesh,phase,color,sign<0?"pure-liquid":"pure-gas",1,0,false);
        const auto load=CchBalancedCsfLoad::constantCurvature(*mesh,CchPressureLoad::build(phase.metrics().divergence()),color,.5,5);
        const Reference ref(*mesh,phase);for(std::size_t u=0;u<ref.area.size();++u)if(ref.minus[u]>=0&&ref.plus[u]>=0)
            check(load.colorGradient()[u]==0&&load.momentumRhsLoad()[u]==0,"pure internal/free CSF not exact zero");}
    auto mesh=makeMesh(4);Case shape;shape.pureFirstJet=true;shape.radial=true;shape.planar=true;shape.planeOffset=.43;
    const auto phase=CchPhaseGeometry::build(*mesh,metricsFor(*mesh,true,3),published(*mesh,shape));
    inspect(*mesh,phase,forceColor(*mesh,phase),"annular",0,3,false);
}
void provenance()
{
    auto mesh=makeMesh(4);Case shape;shape.pureFirstJet=true;shape.planar=true;shape.planeOffset=.43;
    const auto phase=CchPhaseGeometry::build(*mesh,metricsFor(*mesh,false),published(*mesh,shape));
    const auto pressure=CchPressureLoad::build(phase.metrics().divergence());const auto color=forceColor(*mesh,phase);
    const auto load=CchBalancedCsfLoad::constantCurvature(*mesh,pressure,color,.5,5),copy=load;
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
    const auto otherPressure=CchPressureLoad::build(CchIntegralDivergence::build(*mesh,phase.layout()));
    rejection([&]{CchBalancedCsfLoad::constantCurvature(*mesh,otherPressure,color,.5,5);},"divergence owner");
    rejection([&]{load.validateFor(*mesh,otherPressure,color);},"pressure owner");
    const auto rebuiltPressure=CchPressureLoad::build(phase.metrics().divergence());
    const auto rebuilt=CchBalancedCsfLoad::constantCurvature(*mesh,rebuiltPressure,color,.5,5);
    rebuilt.validateFor(*mesh,rebuiltPressure,color);check(rebuilt.momentumRhsLoad()==load.momentumRhsLoad(),"same authoritative D changed CSF construction");
    rejection([&]{load.validateFor(*mesh,rebuiltPressure,color);},"pressure owner");
    const auto foreign=CchPressureLoad::build(metricsFor(*mesh,false).divergence());
    rejection([&]{CchBalancedCsfLoad::constantCurvature(*mesh,foreign,color,.5,5);},"layout");
    for(double sigma:{-1.,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()})
        rejection([&]{CchBalancedCsfLoad::constantCurvature(*mesh,pressure,color,sigma,5);},"surface tension");
    for(double kappa:{std::numeric_limits<double>::infinity(),-std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()})
        rejection([&]{CchBalancedCsfLoad::constantCurvature(*mesh,pressure,color,.5,kappa);},"curvature");
    rejection([&]{CchBalancedCsfLoad::constantCurvature(*mesh,pressure,color,std::numeric_limits<double>::max(),2);},"coefficient");
    rejection([&]{CchBalancedCsfLoad::constantCurvature(*mesh,pressure,color,std::numeric_limits<double>::denorm_min(),.5);},"coefficient");
    rejection([&]{CchBalancedCsfLoad::constantCurvature(*mesh,pressure,color,1e-200,1e-200);},"coefficient");
    auto target=makeMesh(4,1);
    rejection([&]{load.validateFor(*target,pressure,color);},"mesh");
    for(bool rz:{false,true}){shape.radial=rz;const auto original=published(*mesh,shape);HermiteGalsPublicationConfig cfg;cfg.volume.symmetryPlanes=shape.symmetry(*mesh);
        const auto transfer=transferHermiteGalsFirstJetMesh(*mesh,original,*target,HermiteGalsMeshTransferKind::ExactRefinement,cfg);
        check(transfer.report().refinedCells>0&&transfer.publication().state().primary().phaseVolumes.empty(),"actual FirstJet regrid transfer missing");
        const auto nextPhase=CchPhaseGeometry::build(*target,metricsFor(*target,rz),transfer.publication());const auto nextColor=forceColor(*target,nextPhase);
        const auto nextPressure=CchPressureLoad::build(nextPhase.metrics().divergence());
        const auto next=CchBalancedCsfLoad::constantCurvature(*target,nextPressure,nextColor,.5,5);next.validateFor(*target,nextPressure,nextColor);
        algebra(*target,nextPressure,nextColor,Reference(*target,nextPhase),.5,5);
        rejection([&]{next.validateFor(*target,nextPressure,original.state().identity(),nextColor.timeLayer());},"identity");}
    DictionaryMesh::RefineBox box;box.min={.25,.25};box.max={.75,.75};box.level=1;mesh->applyRefineBoxes({box});
    rejection([&]{load.validateFor(*mesh,pressure,color);},"mesh");check(saved==load.momentumRhsLoad(),"failed validation mutated retained load");
}
}
int main(int argc,char** argv)
{
    try{
        std::cout<<std::setprecision(17);
        std::cout<<"BF02_CELL_HEADER,scene,shape,grid,selection,rho,mu,sigma,kappa,pressure,native_cell,component,volume,epsilon,color,sharp_alpha,p,p_expected,component_color_mean\n";
        std::cout<<"BF02_FACE_HEADER,scene,shape,grid,selection,rho,mu,sigma,kappa,velocity,negative,positive,area,free,axis,Cc,S,Cp,Cp_minus_S,Au,reaction,reaction_expected\n";
        const bool domainsOnly=argc==2&&std::string(argv[1])=="--domains";
        check(argc==1||domainsOnly,"unknown BF02 contract option");
        if(domainsOnly)domains();else{core();provenance();}
        if(domainsOnly)check(wallFaces>0&&axisFaces>0,"selected/pure/axis witness missing");
        else check(coarseFaces>0&&axisFaces>0&&wallFaces>0&&negativeControls>0&&guards>0,"core/provenance witness missing");
        std::cout<<"BF02_SUMMARY,"<<(domainsOnly?"domains":"core")<<','<<scenes<<','<<solves<<','<<guards<<','<<negativeControls<<','<<coarseFaces<<','<<axisFaces<<','
            <<wallFaces<<','<<peakLoad<<','<<peakVelocity<<','<<peakPressure<<','<<peakBalance<<','<<peakWall<<','<<peakContinuity<<",PASS\n";return 0;
    }catch(const std::exception& e){std::cerr<<"BF02_CSF_FAIL,"<<e.what()<<'\n';return 1;}
}
