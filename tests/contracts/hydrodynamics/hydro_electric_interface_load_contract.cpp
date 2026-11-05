#include "Modules/Electric/ElectricPotential.h"
#include "Modules/Electric/ElectricMaxwellTrace.h"
#include "Modules/Electric/SurfaceChargeField.h"
#include "Modules/Hydrodynamics/HydroBalancedCsfLoad.h"
#include "Modules/Hydrodynamics/HydroSpatialAdvector.h"
#include "Modules/Hydrodynamics/HydroSinglePhaseSaddle.h"
#include "Modules/Interface/HermiteGalsFirstJetStep.h"
#include "Core/Mesh/UniformLeafGrid.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/IO/DictionaryMesh.h"
#include <algorithm>
#include <array>
#include <sstream>
#include <cmath>
#include <chrono>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
namespace I=LAMP::Interface;
namespace E=LAMP::Electric;
namespace H=LAMP::Hydrodynamics;
constexpr long double pi=3.1415926535897932384626433832795L;
constexpr double radius=.213;
using Clock=I::HermiteGalsTimeLayer;
void check(bool ok,const std::string& why){if(!ok)throw std::runtime_error(why);}
bool sameClock(Clock a,Clock b){return a.step==b.step&&a.begin==b.begin&&a.end==b.end&&a.position==b.position;}
bool sameBits(const std::vector<double>& a,const std::vector<double>& b){
    return a.size()==b.size()&&(a.empty()||std::memcmp(a.data(),b.data(),a.size()*sizeof(double))==0);}
double relative(long double error,long double scale){
    check(std::isfinite(error)&&std::isfinite(scale)&&scale>=0,"invalid independent scale");
    if(scale==0){check(error==0,"unscaled nonzero error");return 0;}return double(std::abs(error)/scale);}
template<class F> void rejected(F action,const char* label){
    try{action();}catch(const std::exception& e){std::cout<<"EL1_REJECT,"<<label<<','<<e.what()<<'\n';return;}
    throw std::runtime_error(std::string("missing rejection: ")+label);}
std::unique_ptr<QuadtreeMesh> meshFor(){
    DictionaryMesh d;d.Nx=d.Ny=8;d.maxInitialIterations=0;d.maxCellNum=512;
    return std::make_unique<QuadtreeMesh>(d);}
I::HermiteGalsPublication sphere(const QuadtreeMesh& mesh){
    return I::initializeHermiteGalsFirstJet(mesh,true,[](Vector2D p){
        const double z=p.y-.5;
        return I::HermiteGalsBoundaryFirstJet{p.x*p.x+z*z-radius*radius,{2*p.x,2*z}};});}
I::HermiteGalsPublication translated(const QuadtreeMesh& mesh,const I::HermiteGalsPublication& source,
    Clock target,double displacement){
    const UniformLeafGrid grid(mesh,CartesianGridScope::UniformLeaves);
    const LevelSetBoundaryTrace::Domain domain{grid.minCoord(),grid.maxCoord(),grid.spacing(),true};
    I::HermiteGalsCandidate authority(mesh,target,source.state().primary());
    const double h=target.time()-source.state().timeLayer().time();
    check(h>0,"real midpoint has no positive physical increment");
    const auto velocity=[speed=displacement/h](Vector2D){return LevelSetBoundaryTrace::VelocityJet{{0,speed},{0,0,0,0}};};
    const I::HermiteGalsVelocityView previous(source.state().identity(),source.state().timeLayer(),velocity);
    const I::HermiteGalsVelocityView arrival(authority.identity(),target,velocity);
    const I::HermiteGalsBoundaryFirstJetQuery boundary=[&mesh,target,displacement](Vector2D p,Clock layer,Boundary){
        check(sameClock(layer,target),"prescribed boundary lost actual fullclock");
        check(p.x==mesh.minCoord.x||p.x==mesh.maxCoord.x||p.y==mesh.minCoord.y||p.y==mesh.maxCoord.y,
            "prescribed boundary used for an interior first-jet reset");
        const double z=p.y-.5-displacement;
        return I::HermiteGalsBoundaryFirstJet{p.x*p.x+z*z-radius*radius,{2*p.x,2*z}};};
    const auto step=I::prepareHermiteGalsFirstJetStep(mesh,source,domain,previous,arrival,boundary);
    check(step.sourceIdentity()==source.state().identity()&&sameClock(step.publication().state().timeLayer(),target),
        "midpoint is not the actual FirstJet/RK2 arrival");
    return step.publication();}
H::HydroPhysicalDivergence metricsFor(const QuadtreeMesh& mesh){
    FaceSegmentGraph graph;graph.build(mesh);
    const auto layout=H::HydroDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,true));
    return H::HydroPhysicalDivergence::build(mesh,H::HydroIntegralDivergence::build(mesh,layout));}
std::vector<int> wallsFor(const QuadtreeMesh& mesh,const H::HydroDofLayout& layout){
    std::vector<int> walls;for(int u=0;u<layout.velocityCount();++u){
        const auto f=layout.domain().face(mesh,layout.nativeFace(u));
        if(f.negativeCell<0||f.positiveCell<0)walls.push_back(u);}return walls;}
// Reference uses physical native edge coordinates and incidence, never C/D coefficients.
struct Incidence {
    std::vector<int> minus,plus;
    std::vector<double> area,volume;
    Incidence(const QuadtreeMesh& mesh,const H::HydroPhaseGeometry& phase){
        const auto& l=phase.layout();
        for(int p=0;p<l.pressureCount();++p){
            const auto& c=mesh.cells(mesh.cellIndexToLabel[l.nativeCell(p)]);
            const auto a=mesh.verts(c.cellToVert[0]).coord(),b=mesh.verts(c.cellToVert[2]).coord();
            volume.push_back(double(pi)*(b.x*b.x-a.x*a.x)*(b.y-a.y));
            check(relative(volume.back()-phase.metrics().cellVolumes()[p],volume.back())<=1e-14,"independent RZ cell volume");}
        for(int u=0;u<l.velocityCount();++u){
            const auto f=l.domain().face(mesh,l.nativeFace(u));
            const auto& e=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
            const auto a=mesh.verts(e.edgeToVert[0]).coord(),b=mesh.verts(e.edgeToVert[1]).coord();
            const double length=f.horizontal?std::abs(b.x-a.x):std::abs(b.y-a.y);
            area.push_back(double(2*pi)*length*(f.horizontal?(a.x+b.x)*.5:a.x));
            minus.push_back(f.negativeCell<0?-1:l.pressureDof(f.negativeCell));
            plus.push_back(f.positiveCell<0?-1:l.pressureDof(f.positiveCell));
            check(relative(area.back()-phase.metrics().faceMeasures()[u],area.back())<=1e-14,"independent RZ face area");}
    }
    std::vector<double> gradient(const std::vector<double>& p)const{
        std::vector<double> result(area.size());
        for(std::size_t u=0;u<area.size();++u)
            result[u]=area[u]*((plus[u]<0?0:p[plus[u]])-(minus[u]<0?0:p[minus[u]]));
        return result;}
    double mean(const std::vector<double>& values)const{
        long double sum=0,total=0;for(std::size_t p=0;p<values.size();++p){sum+=volume[p]*values[p];total+=volume[p];}
        return double(sum/total);}
};
struct PowerReference {
    long double surface=0,flux=0,tangent=0,normalMoment=0,area=0,scale=0;
    // Positive arithmetic scales of five independent scalar integrals.
    std::array<long double,5> positiveScale{},estimatedError{};
    std::size_t evaluations=0;int maximumDepth=0;double tolerance=0;
};
std::array<long double,5> scalarValues(const PowerReference& r){
    return {{r.surface,r.flux,r.tangent,r.normalMoment,r.area}};}
void addReference(PowerReference& a,const PowerReference& b){
    a.surface+=b.surface;a.flux+=b.flux;a.tangent+=b.tangent;
    a.normalMoment+=b.normalMoment;a.area+=b.area;a.scale+=b.scale;
    for(int i=0;i<5;++i){a.positiveScale[i]+=b.positiveScale[i];a.estimatedError[i]+=b.estimatedError[i];}}
double referenceScore(const PowerReference& r){
    double score=0;for(int i=0;i<5;++i)score=std::max(score,relative(r.estimatedError[i],r.positiveScale[i]));
    return score;}
double referenceDifference(const PowerReference& a,const PowerReference& b){
    const auto av=scalarValues(a),bv=scalarValues(b);double error=0;
    for(int i=0;i<5;++i)error=std::max(error,relative(av[i]-bv[i],std::max(a.positiveScale[i],b.positiveScale[i])));
    return error;}
constexpr std::size_t referenceNodeCap=262144;
constexpr int referenceDepthCap=24;
PowerReference scalarSurfaceRule(const H::HydroPhaseGeometry& phase,const E::ElectricMaxwellTrace& maxwell,
    const H::HydroSpatialAdvector& field,double sigma,std::size_t c,std::size_t g,
    double begin,double end,std::size_t& evaluations){
    PowerReference r;const auto& publication=phase.publication();const auto& l=phase.layout();
    // Original positive GL8 subarc rule; scalar contractions use the actual
    // evaluated field, never production native weights, cuts or load reports.
    for(const auto& q:publication.quadrature(c,g,begin,end)){
        check(evaluations<referenceNodeCap,"independent scalar reference exceeded its point-query cap");++evaluations;
        const int p=l.pressureDof(q.sample.frame.nativeCellIndex);
        check(p>=0&&q.physicalWeight>0,"independent surface has no actual native cell/physical weight");
        const auto traction=maxwell.sample(c,g,q.sample);
        const auto velocity=field.cellVelocity(p,q.sample.frame.point);
        const auto n=q.sample.frame.normal;
        const long double nv=static_cast<long double>(n.x)*velocity.x+static_cast<long double>(n.y)*velocity.y;
        const long double fx=static_cast<long double>(traction.tractionJump.x)-sigma*q.sample.curvature*n.x,
            fy=static_cast<long double>(traction.tractionJump.y)-sigma*q.sample.curvature*n.y;
        const long double w=q.physicalWeight,work=w*(fx*velocity.x+fy*velocity.y);
        const long double tangent=w*traction.tangentialTraction
            *(static_cast<long double>(traction.tangent.x)*velocity.x+static_cast<long double>(traction.tangent.y)*velocity.y);
        const double j=sigma*q.sample.curvature-traction.normalTraction;
        r.surface+=work;r.flux+=w*nv;r.tangent+=tangent;
        r.normalMoment+=w*j;r.area+=w;r.scale+=std::abs(work)+std::abs(tangent);
        r.positiveScale[0]+=w*(std::abs(fx*velocity.x)+std::abs(fy*velocity.y));
        r.positiveScale[1]+=w*(std::abs(static_cast<long double>(n.x)*velocity.x)+std::abs(static_cast<long double>(n.y)*velocity.y));
        r.positiveScale[2]+=w*std::abs(traction.tangentialTraction)
            *(std::abs(static_cast<long double>(traction.tangent.x)*velocity.x)+std::abs(static_cast<long double>(traction.tangent.y)*velocity.y));
        r.positiveScale[3]+=w*(std::abs(sigma*q.sample.curvature)+std::abs(traction.normalTraction));
        r.positiveScale[4]+=w;
        check(relative(traction.tangentialTraction-traction.chargeDensity*traction.tangentialField,
            std::max(1.,std::abs(traction.tangentialTraction)))<=1e-12,"actual q_h Et traction identity");
    }
    check(r.area>0,"independent scalar rule lost a positive original subarc");
    return r;}
struct ScalarSurfacePanel {
    std::size_t component=0,graph=0;double begin=0,end=1;int depth=0;
    PowerReference left,right,fine;
};
ScalarSurfacePanel scalarSurfacePanel(const H::HydroPhaseGeometry& phase,const E::ElectricMaxwellTrace& maxwell,
    const H::HydroSpatialAdvector& field,double sigma,std::size_t c,std::size_t g,
    double begin,double end,int depth,std::size_t& evaluations,const PowerReference* reusedCoarse=nullptr){
    ScalarSurfacePanel p;p.component=c;p.graph=g;p.begin=begin;p.end=end;p.depth=depth;
    const double middle=begin+(end-begin)*.5;
    check(middle>begin&&middle<end,"independent scalar reference interval cannot be subdivided");
    const auto coarse=reusedCoarse?*reusedCoarse:scalarSurfaceRule(phase,maxwell,field,sigma,c,g,begin,end,evaluations);
    p.left=scalarSurfaceRule(phase,maxwell,field,sigma,c,g,begin,middle,evaluations);
    p.right=scalarSurfaceRule(phase,maxwell,field,sigma,c,g,middle,end,evaluations);
    addReference(p.fine,p.left);addReference(p.fine,p.right);
    const auto cv=scalarValues(coarse),fv=scalarValues(p.fine);
    for(int i=0;i<5;++i)p.fine.estimatedError[i]=std::abs(fv[i]-cv[i]);
    return p;}
PowerReference surfacePower(const H::HydroPhaseGeometry& phase,const E::ElectricMaxwellTrace& maxwell,
    const H::HydroSpatialAdvector& field,double sigma,double tolerance){
    check(tolerance>0&&tolerance<=1e-12,"invalid independent scalar reference tolerance");
    const auto& publication=phase.publication();std::size_t evaluations=0;
    std::vector<ScalarSurfacePanel> panels;
    for(std::size_t c=0;c<publication.topology().curves.size();++c)
        for(std::size_t g=0;g<publication.topology().curves[c].graphs().size();++g)
            panels.push_back(scalarSurfacePanel(phase,maxwell,field,sigma,c,g,0,1,0,evaluations));
    check(!panels.empty(),"independent scalar reference has no actual interface");
    for(;;){
        PowerReference total;
        for(const auto& p:panels){addReference(total,p.fine);total.maximumDepth=std::max(total.maximumDepth,p.depth);}
        const double score=referenceScore(total);
        if(score<=tolerance){total.evaluations=evaluations;total.tolerance=tolerance;return total;}
        std::size_t worst=0;double worstScore=-1;
        for(std::size_t k=0;k<panels.size();++k){
            double local=0;for(int i=0;i<5;++i)
                local=std::max(local,relative(panels[k].fine.estimatedError[i],total.positiveScale[i]));
            if(local>worstScore){worstScore=local;worst=k;}}
        auto old=std::move(panels[worst]);
        if(old.depth>=referenceDepthCap){
            std::ostringstream why;why<<std::setprecision(17)
                <<"independent scalar reference depth cap, component="<<old.component<<", graph="<<old.graph
                <<", interval=("<<old.begin<<','<<old.end<<"), global_score="<<score<<", requested="<<tolerance;
            throw std::runtime_error(why.str());}
        const double middle=old.begin+(old.end-old.begin)*.5;
        // Retain each child coarse rule; only new quarter rules are evaluated.
        panels[worst]=scalarSurfacePanel(phase,maxwell,field,sigma,old.component,old.graph,
            old.begin,middle,old.depth+1,evaluations,&old.left);
        panels.push_back(scalarSurfacePanel(phase,maxwell,field,sigma,old.component,old.graph,
            middle,old.end,old.depth+1,evaluations,&old.right));
    }
}
void midpointLoad(){
    auto mesh=meshFor();const auto initial=sphere(*mesh);
    const auto electricGeometry=E::ElectricPhaseGeometry::build(*mesh,initial);
    E::SurfaceChargeField initialField(*mesh,initial.topology().curves.front(),0.);
    auto q0=initialField.transport().panelCharge();bool positive=false,negative=false;
    for(std::size_t i=0;i<q0.size();++i){
        q0[i]=initialField.transport().panelMeasure()[i]*(std::sin(1.7*(i+1))+.2+.01*i);
        positive=positive||q0[i]>0;negative=negative||q0[i]<0;}
    check(positive&&negative,"native nonuniform charge lacks both signs");initialField.assignIntegratedCharge(q0);
    const auto moments=E::ElectricSurfaceChargeMomentOperator::prepare(electricGeometry,{&initialField},electricGeometry.timeLayer());
    E::ElectricPhaseProblem problem;problem.epsilonByPhase={{5.,1.}};
    for(int side=0;side<4;++side){
        problem.boundaries[side].type=side==0?E::ElectricCartesianBoundary::Type::Neumann:E::ElectricCartesianBoundary::Type::Dirichlet;
        problem.boundaries[side].value=side==0?std::function<double(Vector2D)>([](Vector2D){return 0.;}):
            std::function<double(Vector2D)>([](Vector2D p){return .05*p.y;});}
    const auto initialPotential=E::solveElectricPhaseElectrostatics(electricGeometry,problem,moments.apply(q0));
    const auto solver=E::ElectricFixedGeometryChargeSolver::prepare(initialPotential,moments,{{.7,.03}});
    auto state=E::ElectricSurfaceChargeTimeState::initialize(solver,{&initialField},electricGeometry.timeLayer());
    const auto acceptedClock=state.timeLayer(),acceptedGeometryClock=state.geometrySourceTimeLayer();
    const auto* acceptedField=&state.field(0);const auto* acceptedSolver=&state.solver();
    const auto unchanged=[&]{
        check(sameBits(state.integratedCharge(),q0)&&sameBits(initialField.transport().panelCharge(),q0)
            &&sameClock(state.timeLayer(),acceptedClock)&&sameClock(state.geometrySourceTimeLayer(),acceptedGeometryClock)
            &&!state.previousTimeLayer()&&!state.previousIntegratedCharge()&&&state.field(0)==acceptedField
            &&&state.solver()==acceptedSolver,"readonly load/preparation/rejection mutated accepted electric state");};
    const Clock midpoint{1,0,.06,Clock::Position::Midpoint};constexpr double displacement=.004;
    const auto target=translated(*mesh,initial,midpoint,displacement);
    const auto targetElectric=E::ElectricPhaseGeometry::build(*mesh,target);
    std::size_t materialQueries=0;
    const auto transfer=initialField.prepareMaterialTransfer(*mesh,LevelSetInterfacePanels(target.topology().curves.front()),
        [&](Vector2D p){++materialQueries;return p+Vector2D{0,displacement};},1,0,midpoint.time());
    // This is a real E-C1 candidate, not a retagged End solve.
    const auto stage=state.prepareMovingStage(*mesh,targetElectric,{transfer},midpoint,problem,{{.7,.03}});
    check(stage.step().report().dt==midpoint.time()&&materialQueries==initialField.geometry().panels().size()+1,
        "midpoint used wrong physical dt or retraced original material endpoints");
    unchanged();const auto maxwell=E::ElectricMaxwellTrace::prepare(stage.step().potential());
    const auto view=maxwell.surfaceVectorView();
    const auto phase=H::HydroPhaseGeometry::build(*mesh,metricsFor(*mesh),target);
    const auto& l=phase.layout();const auto pressure=H::HydroPressureLoad::build(phase.metrics().divergence());
    const auto geometry=H::HydroMomentumGeometry::build(*mesh,phase.metrics());
    const auto walls=wallsFor(*mesh,l);const auto boundary=H::HydroAdvectorBoundary::impermeable(*mesh,geometry,walls);
    const auto reconstruction=H::HydroAdvectorReconstruction::prepare(*mesh,geometry,boundary);
    constexpr double sigma=.7;
    const auto load=H::HydroBalancedCsfLoad::adjointInterfaceLoad(*mesh,pressure,phase,sigma,view,reconstruction,&boundary);
    load.validateFor(*mesh,pressure,phase.geometryIdentity(),midpoint,view);
    load.validateFor(*mesh,pressure,phase.geometryIdentity(),midpoint,reconstruction,&boundary);
    check(load.hasExternalSurfaceVector()&&load.externalSurfaceVector().sameIdentity(view)
        &&&load.phaseGeometry().publication().state().geometry().snapshot()==&target.state().geometry().snapshot(),
        "load did not retain actual electric/geometry owner");
    std::vector<double> velocity(l.velocityCount());
    for(int u=0;u<l.velocityCount();++u)velocity[u]=std::sin(.731*(u+1))+.3*std::cos(.17*(u+1)*(u+3));
    for(int u:walls)velocity[u]=0;
    const auto field=reconstruction.evaluate(*mesh,H::HydroNativeValues(l,velocity,std::vector<double>(l.pressureCount())),
        phase.geometryIdentity(),midpoint,boundary);
    const Incidence incidence(*mesh,phase);
    // Control each raw scalar integral independently: normal flux and surface
    // work can otherwise conceal reference error through the J0 flux remainder.
    const auto coarse=surfacePower(phase,maxwell,field,sigma,1e-12),fine=surfacePower(phase,maxwell,field,sigma,1e-13);
    check(fine.area>0,"actual electric surface area missing");
    const double j0=double(fine.normalMoment/fine.area);
    std::vector<double> potential(l.pressureCount());
    for(int p=0;p<l.pressureCount();++p){
        check(relative(load.referenceNormalStresses()[p]-j0,std::max(1.,std::abs(j0)))<=1e-11,
            "reference stress not the actual positive physical surface mean");
        potential[p]=j0*phase.liquidFractions()[p];}
    const auto base=incidence.gradient(potential);
    long double native=0,basePower=0,tangentPower=0,scale=0,tangentNorm=0;
    for(int u=0;u<l.velocityCount();++u){
        check(relative(load.referencePressureLoad()[u]-base[u],std::max(1.,std::abs(base[u])))<=1e-12,
            "actual common normal reference differs from independent C incidence");
        native+=velocity[u]*load.momentumRhsLoad()[u];basePower+=velocity[u]*base[u];
        tangentPower+=velocity[u]*load.tangentialRhsLoad()[u];
        tangentNorm+=std::abs(load.tangentialRhsLoad()[u]);scale+=std::abs(velocity[u]*load.momentumRhsLoad()[u]);}
    const long double expected=fine.surface+basePower+j0*fine.flux;
    const double error=relative(native-expected,std::max(1.L,scale+fine.scale+std::abs(basePower)));
    const double tangentError=relative(tangentPower-fine.tangent,std::max(1.L,fine.scale));
    const double quadratureError=referenceDifference(coarse,fine);
    std::cout<<"ELECTRIC_HHO_EL1_POWER_CHECK,native="<<double(native)<<",expected="<<double(expected)
        <<",surface="<<double(fine.surface)<<",base="<<double(basePower)<<",flux="<<double(fine.flux)
        <<",j0="<<j0<<",native_scale="<<double(scale)<<",surface_scale="<<double(fine.scale)
        <<",base_abs="<<double(std::abs(basePower))<<",fullvector_relative="<<error
        <<",tangent_native="<<double(tangentPower)<<",tangent_surface="<<double(fine.tangent)
        <<",tangent_relative="<<tangentError<<",coarse_surface="<<double(coarse.surface)
        <<",fine_surface="<<double(fine.surface)<<",coarse_tolerance="<<coarse.tolerance
        <<",fine_tolerance="<<fine.tolerance<<",coarse_estimated_relative="<<referenceScore(coarse)
        <<",fine_estimated_relative="<<referenceScore(fine)<<",coarse_evaluations="<<coarse.evaluations
        <<",fine_evaluations="<<fine.evaluations<<",coarse_depth="<<coarse.maximumDepth
        <<",fine_depth="<<fine.maximumDepth<<",five_scalar_reference_relative="<<quadratureError<<std::endl;
    check(error<=1e-11&&tangentError<=1e-11&&quadratureError<=1e-11,"independent full-vector virtual work");
    check(tangentNorm>1e-10&&std::abs(fine.tangent)>1e-10,"actual q_h Et was omitted from native momentum load");
    const auto electricOnly=H::HydroBalancedCsfLoad::adjointInterfaceLoad(*mesh,pressure,phase,0,view,reconstruction,&boundary);
    long double electricNorm=0,electricWork=0,electricScale=0;
    for(int u=0;u<l.velocityCount();++u){electricNorm+=std::abs(electricOnly.momentumRhsLoad()[u]);
        electricWork+=velocity[u]*electricOnly.momentumRhsLoad()[u];electricScale+=std::abs(velocity[u]*electricOnly.momentumRhsLoad()[u]);}
    const auto electrical=surfacePower(phase,maxwell,field,0,1e-13);
    const double electricalJ0=double(electrical.normalMoment/electrical.area);
    for(int p=0;p<l.pressureCount();++p)potential[p]=electricalJ0*phase.liquidFractions()[p];
    const auto electricalBase=incidence.gradient(potential);long double electricalBaseWork=0;
    for(int u=0;u<l.velocityCount();++u)electricalBaseWork+=velocity[u]*electricalBase[u];
    const double zeroSigmaError=relative(electricWork-electrical.surface-electricalBaseWork-electricalJ0*electrical.flux,
        std::max(1.L,electricScale+electrical.scale+std::abs(electricalBaseWork)));
    check(electricNorm>1e-10&&zeroSigmaError<=1e-11,"sigma_capillary=0 discarded actual nonzero electric force");
    unchanged();

    int guards=0;const auto reject=[&](auto action,const char* label){rejected(action,label);++guards;unchanged();};
    Clock samePhysicalTime=midpoint;samePhysicalTime.begin=.01;samePhysicalTime.end=.05;
    check(samePhysicalTime.time()==midpoint.time(),"fullclock negative control lost equal physical time");
    reject([&]{load.validateFor(*mesh,pressure,phase.geometryIdentity(),samePhysicalTime,view);},"equal-time foreign fullclock");
    const auto oldView=E::ElectricMaxwellTrace::prepare(initialPotential).surfaceVectorView();
    reject([&]{H::HydroBalancedCsfLoad::adjointInterfaceLoad(*mesh,pressure,phase,sigma,oldView,reconstruction,&boundary);},
        "Initial electric view on actual Midpoint");
    const auto foreignPublication=translated(*mesh,initial,midpoint,displacement);
    reject([&]{load.validateFor(*mesh,pressure,foreignPublication.state().identity(),midpoint,view);},"foreign same-shape geometry identity");
    const auto otherPressure=H::HydroPressureLoad::build(H::HydroIntegralDivergence::build(*mesh,l));
    reject([&]{H::HydroBalancedCsfLoad::adjointInterfaceLoad(*mesh,otherPressure,phase,sigma,view,reconstruction,&boundary);},
        "foreign native divergence/pressure owner");
    const auto freshReconstruction=H::HydroAdvectorReconstruction::prepare(*mesh,geometry,boundary);
    reject([&]{load.validateFor(*mesh,pressure,phase.geometryIdentity(),midpoint,freshReconstruction,&boundary);},
        "fresh reconstruction owner");
    const auto foreignBoundary=H::HydroAdvectorBoundary::impermeable(*mesh,geometry,walls);
    reject([&]{load.validateFor(*mesh,pressure,phase.geometryIdentity(),midpoint,reconstruction,&foreignBoundary);},
        "foreign physical boundary owner");
    const auto newView=E::ElectricMaxwellTrace::prepare(stage.step().potential()).surfaceVectorView();
    reject([&]{load.validateFor(*mesh,pressure,phase.geometryIdentity(),midpoint,newView);},"different retained electric view owner");
    reject([&]{state.commit(stage);},"readonly Midpoint remains uncommittable after load preparation");
    const auto retainedView=[&]{
        const auto temporaryProducer=E::ElectricMaxwellTrace::prepare(stage.step().potential());
        return temporaryProducer.surfaceVectorView();}();
    const auto q=target.sample(0,0,.37);
    const auto retained=retainedView.sample(0,0,q),original=maxwell.sample(0,0,q).tractionJump;
    check(retained.x==original.x&&retained.y==original.y,"retained view lost original tensor after producer wrapper release");
    unchanged();
    std::cout<<"ELECTRIC_HHO_EL1_MIDPOINT,RZ,cells="<<mesh->cellIndexNum()<<",h="<<stage.step().report().dt
        <<",epsilon=5:1,conductivity=0.7:0.03,signed_Q=retained,material_queries="<<materialQueries
        <<",fullvector_power_relative="<<error<<",tangent_power_relative="<<tangentError
        <<",independent_quadrature_relative="<<quadratureError<<",tangent_native_l1="<<double(tangentNorm)
        <<",sigma_zero_electric_l1="<<double(electricNorm)<<",sigma_zero_power_relative="<<zeroSigmaError
        <<",atomic_guards="<<guards<<",accepted_history=unchanged,producer_wrapper_lifetime=retained,PASS\n";
}
void constantElectricStress(){
    const auto planeStart=std::chrono::steady_clock::now();auto planeStamp=planeStart;
    const auto mark=[&](const char* stage){
        const auto now=std::chrono::steady_clock::now();
        std::cout<<"ELECTRIC_HHO_EL1_PLANE_STAGE,"<<stage
            <<",elapsed_s="<<std::chrono::duration<double>(now-planeStart).count()
            <<",section_s="<<std::chrono::duration<double>(now-planeStamp).count()<<std::endl;
        planeStamp=now;};
    mark("publication_begin");
    auto mesh=meshFor();constexpr double z0=.483;
    // The actual plane is even in r at the outer section r=1. Declare its
    // sliding scalar endpoint there; the native RZ owner already supplies r=0.
    I::HermiteGalsPublicationConfig publicationConfig;
    publicationConfig.volume.symmetryPlanes.push_back({0,mesh->maxCoord.x});
    const auto publication=I::initializeHermiteGalsFirstJet(*mesh,true,[z0](Vector2D p){
        return I::HermiteGalsBoundaryFirstJet{p.y-z0,{0,1}};},publicationConfig);
    mark("publication_ready");
    check(publication.topology().curves.size()==1,"true RZ plane is not one native interface component");
    mark("electric_geometry_begin");
    const auto electricGeometry=E::ElectricPhaseGeometry::build(*mesh,publication);
    mark("electric_geometry_ready");
    mark("charge_moments_begin");
    E::SurfaceChargeField charge(*mesh,publication.topology().curves.front(),2.);
    const auto accepted=charge.transport().panelCharge();
    const auto moments=E::ElectricSurfaceChargeMomentOperator::prepare(electricGeometry,{&charge},publication.state().timeLayer());
    mark("charge_moments_ready");
    E::ElectricPhaseProblem problem;problem.epsilonByPhase={{5.,1.}};
    const auto exactPotential=[z0](Vector2D p){return (p.y<z0?1.:3.)*(p.y-z0);};
    for(int side=0;side<4;++side){
        problem.boundaries[side].type=side==0?E::ElectricCartesianBoundary::Type::Neumann:E::ElectricCartesianBoundary::Type::Dirichlet;
        problem.boundaries[side].value=side==0?std::function<double(Vector2D)>([](Vector2D){return 0.;}):
            std::function<double(Vector2D)>(exactPotential);}
    mark("electric_solve_begin");
    const auto potential=E::solveElectricPhaseElectrostatics(electricGeometry,problem,moments.apply(accepted));
    mark("electric_solve_ready");
    mark("maxwell_begin");
    const auto maxwell=E::ElectricMaxwellTrace::prepare(potential);const auto view=maxwell.surfaceVectorView();
    mark("maxwell_ready");
    double tensorError=0;std::size_t samples=0;
    for(auto face:maxwell.interfaceFaces())for(const auto& q:maxwell.quadrature(face)){
        check(std::abs(q.normal.x)<=1e-12&&std::abs(q.normal.y-1)<=1e-12,"plane normal is not liquid to gas");
        tensorError=std::max({tensorError,std::abs(q.normalDisplacement[0]+5),std::abs(q.normalDisplacement[1]+3),
            std::abs(q.chargeDensity-2),std::abs(q.normalTraction-2),std::abs(q.tangentialTraction),
            std::abs(q.tractionJump.x),std::abs(q.tractionJump.y-2)});
        ++samples;}
    check(samples>0&&tensorError<=1e-10,"actual piecewise-linear electric plane did not reproduce constant Maxwell stress");
    mark("hydro_phase_reconstruction_begin");
    const auto phase=H::HydroPhaseGeometry::build(*mesh,metricsFor(*mesh),publication);const auto& l=phase.layout();
    const auto pressure=H::HydroPressureLoad::build(phase.metrics().divergence());
    const auto geometry=H::HydroMomentumGeometry::build(*mesh,phase.metrics());
    const auto walls=wallsFor(*mesh,l);const auto boundary=H::HydroAdvectorBoundary::impermeable(*mesh,geometry,walls);
    const auto reconstruction=H::HydroAdvectorReconstruction::prepare(*mesh,geometry,boundary);
    mark("hydro_phase_reconstruction_ready");
    mark("common_load_begin");
    const auto load=H::HydroBalancedCsfLoad::adjointInterfaceLoad(*mesh,pressure,phase,0,view,reconstruction,&boundary);
    load.validateFor(*mesh,pressure,phase.geometryIdentity(),phase.timeLayer(),view);
    mark("common_load_ready");
    const Incidence incidence(*mesh,phase);constexpr double jump=-2;
    std::vector<double> pressurePotential(l.pressureCount());
    for(int p=0;p<l.pressureCount();++p){
        check(std::abs(load.referenceNormalStresses()[p]-jump)<=1e-10,"J=sigma*kappa-tau_n electric sign");
        check(std::abs(load.referenceCurvatures()[p])<=1e-12,"same-scalar plane curvature");
        pressurePotential[p]=jump*phase.liquidFractions()[p];}
    const auto expected=incidence.gradient(pressurePotential);
    long double loadError=0,loadScale=0,residualNorm=0,tangentNorm=0,boundaryReference=0;
    for(int u=0;u<l.velocityCount();++u){
        loadError+=std::abs(load.momentumRhsLoad()[u]-expected[u]);loadScale+=std::abs(expected[u]);
        residualNorm+=std::abs(load.adjointResidualLoad()[u]);tangentNorm+=std::abs(load.tangentialRhsLoad()[u]);
        if(incidence.minus[u]<0||incidence.plus[u]<0)boundaryReference+=std::abs(expected[u]);}
    check(boundaryReference>0,"constant J omitted nonzero full native boundary rows");
    const double fullRows=relative(loadError,std::max(1.L,loadScale));
    check(fullRows<=1e-11&&relative(residualNorm,std::max(1.L,loadScale))<=1e-11
        &&relative(tangentNorm,std::max(1.L,loadScale))<=1e-11,"actual constant electric J lost exact native pressure pairing");
    mark("saddle_begin");
    const auto mass=H::HydroMomentumMass::build(*mesh,geometry,phase,100,1);
    H::HydroSaddleVelocityBlock::Matrix block(l.velocityCount(),l.velocityCount());
    for(int u=0;u<l.velocityCount();++u)block.coeffRef(u,u)=mass.masses()[u]/.125;
    block.makeCompressed();
    const auto saddle=H::HydroSinglePhaseSaddle::buildClosed(pressure,H::HydroSaddleVelocityBlock::build(l,block));
    const auto solved=saddle.solve(*mesh,H::HydroSaddleRhs(l,load.momentumRhsLoad(),std::vector<double>(l.velocityCount())));
    mark("saddle_ready");
    check(saddle.pressureGaugeCount()==1,"plane selected domain acquired a false pressure component");
    const double mean=incidence.mean(phase.liquidFractions());
    std::vector<double> physicalPressure(l.pressureCount()),constantPressure(l.pressureCount(),-jump*mean);
    double pressureError=0,umax=0,uscale=0;
    for(int p=0;p<l.pressureCount();++p){
        physicalPressure[p]=jump*(phase.liquidFractions()[p]-mean);
        pressureError=std::max(pressureError,std::abs(solved.values().pressure()[p]-physicalPressure[p]));}
    const auto actualPressureLoad=incidence.gradient(solved.values().pressure()),gaugeReaction=incidence.gradient(constantPressure);
    std::vector<bool> free(l.velocityCount(),false);for(int u:saddle.freeVelocityDofs())free[u]=true;
    long double balanceError=0,balanceScale=0,reactionError=0,reactionScale=0;
    for(int u=0;u<l.velocityCount();++u){
        umax=std::max(umax,std::abs(solved.values().normalVelocity()[u]));
        uscale=std::max(uscale,std::abs(load.momentumRhsLoad()[u])/block.coeff(u,u));
        if(free[u]){
            balanceError+=std::abs(actualPressureLoad[u]-expected[u]);
            balanceScale+=std::abs(actualPressureLoad[u])+std::abs(expected[u]);
        }else{
            check(solved.values().normalVelocity()[u]==0,"prescribed physical wall velocity changed");
            reactionError+=std::abs(actualPressureLoad[u]-expected[u]-gaugeReaction[u]);
            reactionScale+=std::abs(actualPressureLoad[u])+std::abs(expected[u])+std::abs(gaugeReaction[u]);}}
    const double velocityError=relative(umax,uscale),pError=relative(pressureError,std::abs(jump));
    const double balance=relative(balanceError,std::max(1.L,balanceScale));
    const double reaction=relative(reactionError,std::max(1.L,reactionScale));
    check(velocityError<=1e-11&&pError<=1e-10&&balance<=1e-11&&reaction<=1e-11,
        "actual constant electric stress failed inviscid native static balance/physical boundary reactions");
    check(sameBits(charge.transport().panelCharge(),accepted),"readonly plane electric load mutated accepted charge");
    mark("verification_ready");
    std::cout<<"ELECTRIC_HHO_EL1_CONSTANT_STRESS,RZ,cells="<<mesh->cellIndexNum()
        <<",true_linear_V=1:3,Dn=-5:-3,q_h=2,tau_n=2,J=-2,rho=100:1,mu=0"
        <<",tensor_error="<<tensorError<<",full_native_rows_relative="<<fullRows
        <<",velocity_relative="<<velocityError<<",pressure_relative="<<pError
        <<",free_balance_relative="<<balance<<",boundary_reaction_relative="<<reaction<<",PASS\n";
}
}
int main()try{
    std::cout<<std::setprecision(17);
    midpointLoad();constantElectricStress();
    std::cout<<"ELECTRIC_HHO_EL1_SUMMARY,actual_native_midpoint_and_plane=2,fullvector_and_constant_J=PASS\n";
    return 0;
}catch(const std::exception& e){std::cerr<<"ELECTRIC_HHO_EL1_FAIL,"<<e.what()<<'\n';return 1;}
