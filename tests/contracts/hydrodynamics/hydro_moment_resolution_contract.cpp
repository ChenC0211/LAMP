#include "Modules/Interface/Compatibility/InterfaceResolutionPlan.h"
#include "Modules/Interface/Compatibility/InterfaceParametricMotion.h"
#include "Modules/Interface/Compatibility/InterfacePrescribedData.h"
#include "Modules/Hydrodynamics/HydroCompatibleMomentFlow.h"
#include "Modules/Hydrodynamics/HydroStep.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include <iostream>
#include <iomanip>
#include <numeric>

namespace {
using Flow=Quantum::Hydrodynamics::HydroCompatibleMomentFlow;
using Space=Quantum::Hydrodynamics::HydroCompatibleMomentSpace;
namespace IP=InterfaceParametricProjection;
namespace RP=InterfaceResolutionPlan;
void require(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
DictionaryMesh input(int n,int maximum){DictionaryMesh p;p.Nx=p.Ny=n;p.minCoord={0,0};p.maxCoord={1,1};p.maxCellNum=10000;p.maxInitialIterations=0;p.minLevel_neg=p.minLevel_pos={0};p.maxLevel_local={maximum};return p;}
std::unique_ptr<NativeLeafDomain> domain(const QuadtreeMesh& mesh,bool rz){FaceSegmentGraph graph;graph.build(mesh);return std::make_unique<NativeLeafDomain>(NativeLeafDomain::build(mesh,graph,rz));}
double liquid(const InterfaceState& state){const auto& v=state.parametricGeometry().inventory();return std::accumulate(v.begin(),v.end(),0.);}
void localTargets(){
    auto mesh=std::make_unique<QuadtreeMesh>(input(8,3));auto d=domain(*mesh,false);
    auto source=InterfaceParametricMotion::initialize(*mesh,*d,std::vector<IP::Curve>{
        InterfacePrescribedData::parametricEllipse({.22,.22},{.045,.04},24),
        InterfacePrescribedData::parametricEllipse({.7,.7},{.18,.16},24)});
    RP::Options policy;policy.maximumLevel=3;policy.maxCurvatureCells=.5;policy.minimumOpposedDistanceCells=2;
    const auto plan=RP::evaluate(*mesh,*d,source,policy,4096);
    const auto big=mesh->cells(mesh->findCellLabel({.8,.8})).index(),small=mesh->cells(mesh->findCellLabel({.22,.22})).index();
    std::cout<<"LOCAL_RESOLUTION cells="<<mesh->cellIndexNum()<<" planned="<<plan.plannedLeafCount<<" small_level="<<plan.targetLevels[small]<<" large_level="<<plan.targetLevels[big]<<std::endl;
    require(plan.localInterfaceLevels&&!plan.saturated&&plan.targetLevels[small]==3&&plan.targetLevels[big]<3&&plan.plannedLeafCount<4096,"local curvature request became a global finest interface level");
    auto previous=InterfaceParametricMotion::initialize(*mesh,*d,InterfacePrescribedData::parametricEllipse({.7,.7},{.18,.16},24));
    const auto old=RP::evaluate(*mesh,*d,previous,policy,4096),both=RP::evaluate(*mesh,*d,source,previous,policy,4096);
    for(int c=0;c<mesh->cellIndexNum();c++)require(both.targetLevels[c]==std::max(plan.targetLevels[c],old.targetLevels[c]),"history union changed local physical levels");
    const double before=liquid(source);const auto curve=source.parametricGeometry().boundaryHandle();
    auto target=RP::refine(*mesh,plan,4096);auto td=domain(*target,false);
    auto transfer=InterfaceParametricMotion::prepareMeshTransfer(*mesh,*d,source,*target,*td,{});
    require(transfer.state.parametricGeometry().boundaryHandle()==curve&&std::abs(liquid(transfer.state)-before)<2e-15,"local refinement changed geometry or independent liquid");
    const auto stationary=RP::evaluate(*target,*td,transfer.state,policy,4096);
    require(!stationary.saturated&&!stationary.changesMesh(),"stationary local mesh keeps spreading its fine band");
    const auto repeated=RP::evaluate(*target,*td,transfer.state,transfer.state,policy,4096);
    require(!repeated.changesMesh(),"stationary history union keeps spreading its fine band");
    bool budget=false;try{RP::refine(*mesh,plan,int(plan.requiredLeafCount)-1);}catch(const std::length_error&){budget=true;}
    require(budget&&source.parametricGeometry().boundaryHandle()==curve,"hard cell budget was relaxed or source changed");
    std::cout<<std::setprecision(17)<<"LOCAL_RESOLUTION_FIXED cells="<<target->cellIndexNum()<<" volume="<<liquid(transfer.state)-before<<" fixed=true history_local=true budget_strict=true"<<std::endl;
}
void terminalScope(){
    auto mesh=std::make_unique<QuadtreeMesh>(input(8,3));auto d=domain(*mesh,true);
    // In the local neck r(z)=r0+O(z^2), cap feasibility depends on r0/h^2,
    // not only r0/h. Retain the thin-neck scaling of the existing cap control.
    const double finest=std::ldexp(.125,-3),neck=.128*finest*finest;
    auto source=InterfaceParametricMotion::initialize(*mesh,*d,InterfacePrescribedData::parametricNeckedSpheroid(.5,{.35,.3},neck/.35,32));
    RP::Options policy;policy.maximumLevel=3;policy.maxCurvatureCells=.5;policy.minimumOpposedDistanceCells=2;
    const auto strict=RP::evaluate(*mesh,*d,source,policy,4096);require(strict.saturated&&!strict.unresolved.empty(),"underresolved neck lost its rejected geometric witness");
    for(const auto& w:strict.unresolved)if(w.kind==RP::Result::UnresolvedFeature::Kind::Curvature)
        require(w.observed>w.allowed&&w.component==1&&w.point.x>=0,"curvature saturation witness does not explain its rejected request");
    policy.parametricNeckRadiusCells=.1;auto plan=RP::evaluate(*mesh,*d,source,policy,4096);
    std::cout<<std::setprecision(17)<<"TERMINAL_RESOLUTION regions="<<plan.terminalRanges.size()<<" saturated="<<plan.saturated<<" cells="<<plan.plannedLeafCount<<" radius="<<plan.terminalRadius<<" volume="<<plan.terminalGeometricVolume<<std::endl;
    require(!plan.saturated&&plan.terminalRanges.size()==1&&!plan.terminalRanges[0].pole&&plan.terminalCurvatureRequests>0,"bounded terminal neck did not enter one explicit finite-scale region");
    const auto& range=plan.terminalRanges.front();
    Quantum::Geometry::Parametric::PhaseMeasure measure(*source.parametricGeometry().boundaryHandle()->components().front().curve);
    const double measured=measure.rectanglePhaseMoments(0,plan.terminalRadius,range.lowerZ,range.upperZ,true).values[0];
    require(std::abs(measured-plan.terminalGeometricVolume)<2e-15&&measured<liquid(source),"terminal geometric measure lost its physical parameter scale");
    const auto sites=IP::inspectNecks(source.parametricGeometry(),.1);require(!sites.empty(),"missing initial terminal neck");
    bool early=false;try{RP::admitsNeck(*mesh,source,plan,3,sites.front());}catch(const std::logic_error&){early=true;}
    require(early,"terminal model allowed an early coarse-cell split");
    auto target=RP::refine(*mesh,plan,4096);auto td=domain(*target,true);
    auto transfer=InterfaceParametricMotion::prepareMeshTransfer(*mesh,*d,source,*target,*td,{});
    source=std::move(transfer.state);mesh=std::move(target);d=std::move(td);
    plan=RP::evaluate(*mesh,*d,source,policy,4096);require(!plan.saturated&&!plan.changesMesh(),"terminal support is not a fixed local plan");
    const auto fineSites=IP::inspectNecks(source.parametricGeometry(),.1);
    require(fineSites.size()==1&&RP::admitsNeck(*mesh,source,plan,3,fineSites.front()),"finest terminal support did not admit its neck");
    auto pending=InterfaceParametricMotion::initialize(*mesh,*d,InterfacePrescribedData::parametricNeckedSpheroid(.5,{.35,.3},.001/.35,32));
    const auto pendingSites=IP::inspectNecks(pending.parametricGeometry(),.1);
    require(pendingSites.size()==1&&pendingSites.front().radius<.1*finest,"missing subthreshold pending neck");
    const auto pendingGeometry=IP::assessNeckSplit(pending.parametricGeometry(),pendingSites.front().options);
    require(pendingGeometry.readiness==IP::NeckReadiness::CapPending&&!pendingGeometry.split,"radius threshold alone forced a geometrically infeasible split");
    std::cout<<"TERMINAL_PENDING radius="<<pendingSites.front().radius<<" finest="<<finest<<" cap_pending=true retained_parent=true"<<std::endl;
    const double before=liquid(source);auto space=Space::build(*mesh,*d);const Eigen::VectorXd zero=Eigen::VectorXd::Zero(space.mass().rows());
    auto flow=Flow::createTwoPhase(*mesh,*d,std::move(space),std::move(source),{5,1,.02,.01},zero,.1);
    auto split=flow.assessNeckSplit(*mesh,*d,fineSites.front().options);
    std::cout<<"TERMINAL_MOMENT feasibility="<<int(split.readiness)<<" kinetic="<<split.report.momentum.kineticAfter<<" minimum="<<split.report.momentum.minimumKinetic<<" constraint_rank="<<split.report.momentum.independentConstraints<<std::endl;
    require(bool(split.candidate),"finest canonical cap did not enter the joint feasible set");flow.acceptNeckSplit(*mesh,*d,std::move(*split.candidate));
    plan=RP::evaluate(*mesh,*d,flow.interfaceState(),policy,4096);
    std::size_t poles=0;for(const auto& r:plan.terminalRanges)if(r.pole&&r.birthParent==1)++poles;
    std::cout<<"TERMINAL_BIRTH regions="<<plan.terminalRanges.size()<<" poles="<<poles<<" saturated="<<plan.saturated<<" gap_requests="<<plan.terminalGapRequests<<" volume="<<plan.terminalGeometricVolume<<std::endl;
    require(poles==2&&!plan.saturated&&!plan.changesMesh()&&plan.terminalGapRequests>0,"event-born paired poles lost their declared model scope");
    std::vector<IP::Curve> untagged;for(const auto& c:flow.interfaceState().parametricGeometry().boundaryHandle()->components())untagged.push_back(*c.curve);
    auto withoutEvent=InterfaceParametricMotion::initialize(*mesh,*d,std::move(untagged));
    require(RP::evaluate(*mesh,*d,withoutEvent,policy,4096).saturated,"unrecorded sharp poles received a model exception");
    auto closeDrops=InterfaceParametricMotion::initialize(*mesh,*d,std::vector<IP::Curve>{
        InterfacePrescribedData::parametricSpheroid(.39,{.1,.1},24),InterfacePrescribedData::parametricSpheroid(.61,{.1,.1},24)});
    const auto unrelated=RP::evaluate(*mesh,*d,closeDrops,policy,4096);
    require(unrelated.saturated&&unrelated.gapRequests>0&&unrelated.terminalGapRequests==0,"unrelated drop approach was mistaken for a modeled split gap");
    bool gapWitness=false;for(const auto& w:unrelated.unresolved)if(w.kind==RP::Result::UnresolvedFeature::Kind::OpposedGap){
        const auto delta=w.point-w.otherPoint;require(w.observed<w.allowed&&std::abs(std::hypot(delta.x,delta.y)-w.observed)<1e-14,"opposed-gap witness lost its actual endpoints");gapWitness=true;}
    require(gapWitness,"unrelated close drops did not preserve their unresolved gap");
    Flow::Options options;double maxEquation=0,maxVolume=0;
    for(int step=1;step<=4;step++){
        const HydroStepTime time{step,flow.state().time,.00025,HydroTimePolicy::Fixed};auto load=flow.makeLoad(*mesh,*d,Eigen::VectorXd::Zero(flow.space().mass().rows()),time.begin+.5*time.dt);
        auto candidate=flow.prepareNavierStokesStep(*mesh,*d,time,load,options);const auto report=candidate.report();flow.accept(*mesh,*d,std::move(candidate));
        maxVolume=std::max(maxVolume,std::abs(liquid(flow.interfaceState())-before));maxEquation=std::max(maxEquation,report.equationResidual);
        const auto next=RP::evaluate(*mesh,*d,flow.interfaceState(),policy,4096);
        std::cout<<"TERMINAL_STEP step="<<step<<" time="<<flow.state().time<<" saturated="<<next.saturated<<" changes_mesh="<<next.changesMesh()<<" equation="<<report.equationResidual<<std::endl;
        require(!next.saturated&&!next.changesMesh(),"post-event motion escaped the predeclared terminal/support control");
    }
    auto tight=policy;tight.maximumGeometryFractionDifference=1e-16;
    const auto defect=RP::evaluate(*mesh,*d,flow.interfaceState(),tight,4096);
    require(defect.inventoryGeometryRequests>0&&defect.saturated,"terminal scope hid an independent inventory defect");
    require(maxVolume<2e-15&&maxEquation<1e-10,"terminal local AMR continuation lost mass or momentum solve");
    std::cout<<"TERMINAL_CONTINUATION steps=4 time="<<flow.state().time<<" cells="<<mesh->cellIndexNum()<<" volume="<<maxVolume<<" equation="<<maxEquation<<" inventory_strict=true unrelated_gap_strict=true"<<std::endl;
}
}
void hydroMomentLocalResolutionContract(){localTargets();terminalScope();}
