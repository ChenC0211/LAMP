#include "IO/TwoPhaseInputComponents.h"
#include "IO/DictionaryParallel.h"
#include "Core/Mesh/DyadicGridCoordinates.h"
#include "Core/IO/DictionaryGeometry.h"
#include "Core/Units/UnitsContext.h"
#include "Core/Geometry/FittedDomainGeometry.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <fstream>
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>

namespace
{
using Node=DictionaryInput::Node;
using Trace=DictionaryInput::ConfigurationTrace;
using Expression=LAMP::ScalarExpression;
using PinWallTraceLaw=LAMP::Interface::GeometricVofPublication::PinCertificate::WallTraceLaw;
const char* pinWallTraceName(PinWallTraceLaw law)
{
    switch(law){
        case PinWallTraceLaw::StrictWholeCell:return "strictWholeCell";
        case PinWallTraceLaw::PrescribedWetDry:return "prescribedWetDry";
    }
    throw std::invalid_argument("unknown typed CSG pin wall trace law");
}
using namespace LAMP::Hydrodynamics;

std::ostringstream stream(){std::ostringstream out;out<<std::setprecision(17);return out;}
void vec(std::ostream& out,Vector2D v){out<<'('<<v.x<<' '<<v.y<<')';}
Vector2D vector(const Node& node,Trace& trace,const char* units,const char* consumer)
{
    if(node.isBlock || node.values.size()!=4 || node.values.front()!="(" || node.values.back()!=")")
        node.fail("a finite two-component vector (first second)");
    auto item=node;item.values={node.values[1]};const auto first=item.number();
    item.values={node.values[2]};const auto second=item.number();
    auto out=stream();vec(out,{first,second});trace.record(node,out.str(),units,consumer);
    return {first,second};
}
std::array<int,2> resolution(const Node& node,Trace& trace,const char* consumer="Mesh.uniformGrid")
{
    if(node.isBlock || node.values.size()!=4 || node.values.front()!="(" || node.values.back()!=")")
        node.fail("two explicit integer cell counts (nx ny)");
    auto item=node;item.values={node.values[1]};const int nx=item.integer();
    item.values={node.values[2]};const int ny=item.integer();
    if(nx<4 || ny<4 || static_cast<long long>(nx)*ny>1000000)
        node.fail("nx,ny >= 4 and at most 1000000 fixed uniform cells");
    trace.record(node,std::vector<int>{nx,ny},"cells",consumer);return {nx,ny};
}
std::vector<std::string> sides(const Node& node,Trace& trace)
{
    if(node.isBlock || node.values.size()<3 || node.values.front()!="(" || node.values.back()!=")")
        node.fail("a nonempty list of distinct domain sides (left right bottom top)");
    std::vector<std::string> result;
    for(std::size_t i=1;i+1<node.values.size();++i){
        const auto& name=node.values[i];
        if((name!="left" && name!="right" && name!="bottom" && name!="top") ||
           std::find(result.begin(),result.end(),name)!=result.end())
            node.fail("a nonempty list of distinct domain sides (left right bottom top)");
        result.push_back(name);
    }
    trace.record(node,node.raw(),"dimensionless","Interface.FirstJet.symmetryPlanes");
    return result;
}
double bounded(const Node& node,Trace& trace,double lower,bool strictLower,const char* consumer,const char* units="code")
{
    const double value=node.number();
    if(strictLower?!(value>lower):value<lower)node.fail(strictLower?"a finite positive value":"a finite nonnegative value");
    return trace.record(node,value,units,consumer);
}
double option(const Node& block,const char* key,double value,Trace& trace,const char* units,const char* consumer)
{
    const auto* node=block.find(key);
    return node?trace.record(*node,node->number(),units,consumer):trace.defaultValue(block,key,value,units,consumer);
}

// The same canonical level-0 coordinates used by BackgroundMesh, not a
// tolerance or a second a+h*k expression with a different rounding history.
void readVofOpenPatches(HydrodynamicsInput& input,const Node& boundary,Trace& trace)
{
    using Patch=HydrodynamicsInput::VofBoundaryPatch;
    using Mode=LAMP::Interface::GeometricVofAdvectorSlab::Boundary::Mode;
    const auto& patches=boundary.block("Patches");
    if(patches.entries.empty())patches.fail("nonempty named physical boundary patches");
    const LAMP::Mesh::DyadicGridCoordinates coordinates(input.domainMin,input.domainMax,input.nx,input.ny);
    const auto aligned=[&](double x,bool horizontal){int left=0,right=horizontal?input.nx:input.ny;
        while(left<=right){const int k=left+(right-left)/2;const auto point=coordinates.point(horizontal?k:0,horizontal?0:k,0);
            const auto value=horizontal?point.x:point.y;if(value==x)return true;if(value<x)left=k+1;else right=k-1;}
        return false;};
    bool incoming=false,outgoing=false;std::set<double> reservoirPressures;
    for(const auto& node:patches.entries){
        if(!node.isBlock)node.fail("a named physical patch block");
        Patch patch;patch.name=node.name;
        const auto type=trace.record(node.get("type"),node.get("type").choice({"prescribedVelocity","pressureOutlet"}),
            "dimensionless","Hydrodynamics.VOF.physicalPatchKind");
        patch.kind=type=="prescribedVelocity"?Patch::Kind::PrescribedVelocity:Patch::Kind::PressureOutlet;
        if(input.vofPressureDrivenReservoirs&&patch.kind!=Patch::Kind::PressureOutlet)
            node.get("type").fail("only pressureOutlet patches for pressureDrivenReservoirs; no prescribedVelocity or fixed feed");
        if(const auto* id=node.find("boundaryId")){
            if(!input.csgGeometry||!input.vofPressureDrivenReservoirs)
                id->fail("boundaryId only for actual A CSG pressureDrivenReservoirs");
            node.allow({"type","boundaryId","phase","pressure"});
            const int value=trace.record(*id,id->integer(),"1-based Core boundary ID","Hydrodynamics.VOF.physicalPatchBoundaryId");
            const auto found=input.csgGeometry->boundaryById.find(value);
            if(value<1||found==input.csgGeometry->boundaryById.end()||found->second.ignored)
                id->fail("one declared nonignored Core CSG boundary ID");
            for(const auto& previous:input.vofBoundaryPatches)if(previous.boundaryId&&*previous.boundaryId==value)
                id->fail("distinct Core boundary ID across physical flow patches");
            patch.boundaryId=value;
        }else{
            if(input.csgGeometry)node.fail("explicit boundaryId for each actual CSG pressure reservoir; side/interval are exclusive");
            if(patch.kind==Patch::Kind::PrescribedVelocity)node.allow({"type","side","interval","phase","velocity"});
            else node.allow({"type","side","interval","phase","pressure"});
            const auto side=trace.record(node.get("side"),node.get("side").choice({"W","S","E","N"}),
                "dimensionless","Hydrodynamics.VOF.physicalPatchSide");
            patch.side=side=="W"?Patch::Side::W:side=="S"?Patch::Side::S:side=="E"?Patch::Side::E:Patch::Side::N;
            const bool horizontal=patch.side==Patch::Side::S||patch.side==Patch::Side::N;
            if(input.radial&&input.domainMin.x==0&&patch.side==Patch::Side::W)
                node.get("side").fail("Closed true r=0 axis; no physical flow patch");
            patch.interval=vector(node.get("interval"),trace,"code_length","Hydrodynamics.VOF.physicalPatchInterval");
            const auto lo=horizontal?input.domainMin.x:input.domainMin.y,hi=horizontal?input.domainMax.x:input.domainMax.y;
            if(!(patch.interval.x>=lo&&patch.interval.x<patch.interval.y&&patch.interval.y<=hi
                &&aligned(patch.interval.x,horizontal)&&aligned(patch.interval.y,horizontal)))
                node.get("interval").fail("ordered in-domain physical tangent endpoints exactly on actual base-grid coordinates");
            for(const auto& previous:input.vofBoundaryPatches)if(!previous.boundaryId&&previous.side==patch.side&&
                std::max(previous.interval.x,patch.interval.x)<std::min(previous.interval.y,patch.interval.y))
                node.get("interval").fail("nonoverlapping physical patches on each side (shared endpoints are allowed)");
        }
        const auto phase=trace.record(node.get("phase"),node.get("phase").choice({"liquid","gas"}),
            "dimensionless","Interface.VOF.reservoirPhase");
        patch.reservoir=phase=="liquid"?Mode::PureLiquidReservoir:Mode::PureGasReservoir;
        if(patch.kind==Patch::Kind::PrescribedVelocity){
            patch.velocity=vector(node.get("velocity"),trace,"code_velocity","Hydrodynamics.VOF.steadyPhysicalPatchVector");
            const bool horizontal=patch.side==Patch::Side::S||patch.side==Patch::Side::N;
            const auto normal=horizontal?patch.velocity.y:patch.velocity.x;
            if(normal==0)node.get("velocity").fail("nonzero steady normal component of the complete physical vector");
            const auto outward=(patch.side==Patch::Side::W||patch.side==Patch::Side::S)?-normal:normal;
            incoming=incoming||outward<0;outgoing=outgoing||outward>0;
        }else{patch.pressure=trace.record(node.get("pressure"),node.get("pressure").number(),
            "code_pressure","Hydrodynamics.VOF.physicalPressureOutlet");outgoing=true;if(input.vofPressureDrivenReservoirs)reservoirPressures.insert(patch.pressure);}
        input.vofBoundaryPatches.push_back(std::move(patch));
    }
    if(input.vofPressureDrivenReservoirs){
        if(input.vofBoundaryPatches.size()<2||reservoirPressures.size()<2)
            patches.fail("at least two nonoverlapping pressureOutlet reservoir ports with distinct physical pressures");
    }else if(!incoming||!outgoing)patches.fail("at least one fixed incoming patch and a fixed outgoing or free pressure outlet");
    // Native physical flux compatibility and outlet direction are checked by
    // the actual initial/next Saddle. The reader never adjusts any g.
}

void readVofBackflow(HydrodynamicsInput& input,const Node& boundary,Trace& trace)
{
    const auto* selected=boundary.find("Backflow");if(!selected)return;
    const auto& block=boundary.block("Backflow");
    if(std::none_of(input.vofBoundaryPatches.begin(),input.vofBoundaryPatches.end(),[](const auto& patch){
        return patch.kind==HydrodynamicsInput::VofBoundaryPatch::Kind::PressureOutlet;}))
        block.fail("explicit backflow requires an actual pressureOutlet physical patch");
    using Model=HydroVofPressureOutletMomentumModel;
    HydrodynamicsInput::VofPressureOutletBackflow config;
    const auto model=trace.record(block.get("model"),block.get("model").choice({"ambient","interiorContinuationZeroShear"}),
        "dimensionless","HydroVofStep.prescribedOpenMomentum");
    if(model=="ambient"){
        block.allow({"model","velocity"});config.model=Model::Ambient;
        config.velocity=vector(block.get("velocity"),trace,"code_velocity","HydroVofStep.pressureOutletBackflowMomentum");
    }else{
        block.allow({"model"});config.model=Model::InteriorContinuationZeroShear;
    }
    input.vofPressureOutletBackflow=std::move(config);
}

void kernelDefaults(const Node& hydro,Trace& trace,const HydroTwoPhaseStepConfig& step)
{
    const HydroSaddleConfig saddle;
    trace.defaultValue(hydro,"KernelDefaults.linear.method",std::string("EigenSparseLU"),"dimensionless","HydroSinglePhaseSaddle.quotientFactor");
    trace.defaultValue(hydro,"KernelDefaults.linear.algebraicRelativeTolerance",saddle.algebraicRelativeTolerance,"dimensionless","HydroSinglePhaseSaddle.solve");
    trace.defaultValue(hydro,"KernelDefaults.linear.momentumRelativeTolerance",saddle.momentumRelativeTolerance,"dimensionless","HydroSinglePhaseSaddle.solve");
    trace.defaultValue(hydro,"KernelDefaults.linear.continuityRelativeTolerance",saddle.continuityRelativeTolerance,"dimensionless","HydroSinglePhaseSaddle.solve");
    trace.defaultValue(hydro,"KernelDefaults.linear.pressureMeanRelativeTolerance",saddle.pressureMeanRelativeTolerance,"dimensionless","HydroSinglePhaseSaddle.solve");
    trace.defaultValue(hydro,"KernelDefaults.linear.boundaryFluxRelativeTolerance",saddle.boundaryFluxRelativeTolerance,"dimensionless","HydroSinglePhaseSaddle.solve");
    trace.defaultValue(hydro,"KernelDefaults.momentumPowerTolerance",step.momentumPowerTolerance,"dimensionless","HydroTwoPhaseStep.acceptance");
    trace.defaultValue(hydro,"KernelDefaults.material.integrationRule",std::string("embeddedSelectiveBoole"),"dimensionless","HydroMaterialTransport.prepareSpaceTime");
    trace.defaultValue(hydro,"KernelDefaults.material.relativeTolerance",step.spaceTime.relativeTolerance,"dimensionless","HydroMaterialTransport.prepareSpaceTime");
    trace.defaultValue(hydro,"KernelDefaults.material.maxDepth",step.spaceTime.maxDepth,"levels","HydroMaterialTransport.prepareSpaceTime");
    trace.defaultValue(hydro,"KernelDefaults.material.maxGeometrySamples",step.spaceTime.maxGeometrySamples,"samples","HydroMaterialTransport.prepareSpaceTime");
    trace.defaultValue(hydro,"KernelDefaults.material.firstJetTimeScheme",std::string("characteristicPullback"),"dimensionless","HydroMaterialTransport.FirstJet");
    trace.defaultValue(hydro,"KernelDefaults.material.pointParallel",step.spaceTime.enablePointParallel,"dimensionless","HydroMaterialTransport.FirstJet");
    trace.defaultValue(hydro,"KernelDefaults.material.numericalTransportSampling",step.spaceTime.enableNumericalTransportSampling,"dimensionless","HydroMaterialTransport.FirstJet");
    trace.defaultValue(hydro,"KernelDefaults.material.atomLocalTimeIntegration",step.spaceTime.enableAtomLocalTimeIntegration,"dimensionless","HydroMaterialTransport.prepareSpaceTime");
    trace.defaultValue(hydro,"KernelDefaults.material.transportValueAbsoluteTolerance",step.spaceTime.transportValueAbsoluteTolerance,"code_phi","HydroMaterialTransport.FirstJet");
    trace.defaultValue(hydro,"KernelDefaults.material.transportGradientAbsoluteTolerance",step.spaceTime.transportGradientAbsoluteTolerance,"code_phi/code_length","HydroMaterialTransport.FirstJet");
    trace.defaultValue(hydro,"KernelDefaults.interface.relativeInterfaceResidualTolerance",step.publication.relativeInterfaceResidualTolerance,"dimensionless","Interface.HermiteGalsPublication");
    trace.defaultValue(hydro,"KernelDefaults.interface.unitNormalTolerance",step.publication.unitNormalTolerance,"dimensionless","Interface.HermiteGalsPublication");
}
}
namespace AppIO
{
std::array<std::string,3> readVofMethodInput(HydrodynamicsInput&,const Node&,const Node&,const Node&);
void recordVofInputDefaults(const Node&,const HydrodynamicsInput&,Trace&);
TwoPhaseInputComponents readTwoPhaseInputComponents(HydrodynamicsInput& input,
    const CaseDocuments& documents,ConstantsUsage& usage,std::optional<int> threadsOverride,int checkpointDefault)
{
    TwoPhaseInputComponents components;auto& trace=input.trace;
    documents.root("interface.dict").allow({"Interface"});
    const auto& interface=documents.root("interface.dict").block("Interface");
    const auto method=interface.get("method").choice({"hermiteGals","vof"});
    const bool vof=method=="vof",nativeGeometry=vof;
    const bool ehdFamily=documents.app()=="VofElectrohydrodynamics"||documents.app()=="Electrospray";
    if(vof&&documents.app()!="Hydrodynamics"&&!ehdFamily)interface.get("method").fail("vof only for Hydrodynamics or the VOF-EHD production family");
    input.interfaceMethod=vof?HydrodynamicsInput::InterfaceMethod::Vof:HydrodynamicsInput::InterfaceMethod::HermiteGals;
    trace.record(interface.get("method"),method,"dimensionless",vof?"Interface.VOF.dispatch":"Interface.FirstJet.dispatch");
    const bool csgInitial=documents.root("geometry.dict").block("Geometry").block("Initial").find("csgFile")!=nullptr;
    if(csgInitial&&!vof)documents.root("geometry.dict").block("Geometry").block("Initial").fail("CSG initial source is supported only for VOF");
    // Choose the capability default before tracing; never silently replace 10.
    if(csgInitial)checkpointDefault=0;
    const auto& run=documents.root("run.dict");
    if(vof)run.allow({"Run","Time","Output","Log","Parallel","Checkpoint","Restart"});
    else run.allow({"Run","Time","Output","Log","Parallel"});
    const char* const app=documents.app().c_str();
    input.parallelThreads=AppRuntime::readParallelThreads(run,&trace,threadsOverride);
    const auto& dispatch=run.block("Run");dispatch.allow({"app"});
    trace.record(dispatch.get("app"),dispatch.get("app").choice({app}),"dimensionless","Runtime.App.dispatch");
    if(vof){
        const char* consumer="Runtime.VOF.acceptedCheckpoint";
        if(checkpointDefault<0||checkpointDefault>1000000)
            throw std::logic_error("invalid caller-owned VOF checkpoint default");
        const auto* checkpoint=run.find("Checkpoint");input.run.checkpointExplicit=checkpoint!=nullptr;
        if(checkpoint){
            checkpoint->allow({"everySteps"});const auto& interval=checkpoint->get("everySteps");
            input.run.checkpointEverySteps=trace.record(interval,interval.integer(),"accepted steps",consumer);
            if(input.run.checkpointEverySteps<0||input.run.checkpointEverySteps>1000000)
                interval.fail("0 (disabled) or an integer accepted-step interval in [1,1000000]");
        }else{
            input.run.checkpointEverySteps=checkpointDefault;
            trace.entries.push_back({"Checkpoint.everySteps",run.file,"<omitted>","accepted steps",consumer,
                run.line,checkpointDefault,true,true,{}});
        }
        if(const auto* restart=run.find("Restart")){
            restart->allow({"directory"});const auto& directory=restart->get("directory");
            const auto value=directory.word();const auto path=std::filesystem::u8path(value);
            if(value.empty()||value.find_first_of("\"\r\n")!=std::string::npos||value.find('\0')!=std::string::npos||!path.is_absolute())
                directory.fail("a nonempty absolute checkpoint directory without newlines or NUL");
            input.run.restartDirectory=trace.record(directory,value,"absolute path",consumer);
        }
    }else input.run.checkpointEverySteps=0;
    if(csgInitial){
        if(input.run.checkpointEverySteps>0)run.block("Checkpoint").get("everySteps").fail("unsupported CSG initial checkpoint; use everySteps 0");
        if(!input.run.restartDirectory.empty())run.block("Restart").fail("unsupported CSG initial restart until selected-domain checkpoint support");
    }
    const auto& time=run.block("Time");
    const bool allowAdaptive=vof&&(documents.app()=="Hydrodynamics"||ehdFamily);
    if(allowAdaptive)time.allow({"units","dt","endTime","StepRetry","Adaptive"});
    else time.allow({"units","dt","endTime","StepRetry"});
    trace.record(time.get("units"),time.get("units").choice({"code"}),"dimensionless","Runtime.clock.units");
    input.run.dt=bounded(time.get("dt"),trace,0,true,"Runtime.clock","code_time");
    input.run.endTime=bounded(time.get("endTime"),trace,0,true,"Runtime.clock","code_time");
    if(!time.find("Adaptive")){
        const auto count=input.run.fixedStepCount();
        if(!count || *count>1000000)time.get("endTime").fail("an integer number of fixed dt steps in [1,1000000]");
        input.steps=*count;
    }
    if(const auto* retry=time.find("StepRetry")){
        retry->allow({"minDt","maxRetries","recoverAfter"});
        RunControls::StepRetry policy;
        policy.minDt=bounded(retry->get("minDt"),trace,0,true,"Runtime.transactionalStepRetry","code_time");
        policy.maxRetries=trace.record(retry->get("maxRetries"),retry->get("maxRetries").integer(),"attempts","Runtime.transactionalStepRetry");
        policy.recoverAfter=retry->find("recoverAfter")
            ?trace.record(retry->get("recoverAfter"),retry->get("recoverAfter").integer(),"accepted steps","Runtime.transactionalStepRetry")
            :trace.defaultValue(*retry,"recoverAfter",0,"accepted steps","Runtime.transactionalStepRetry");
        if(policy.minDt>input.run.dt)retry->get("minDt").fail("positive minDt no greater than Time.dt");
        if(policy.maxRetries<1 || policy.maxRetries>30)retry->get("maxRetries").fail("an integer retry budget in [1,30]");
        if(policy.recoverAfter<0 || policy.recoverAfter>1000000)retry->get("recoverAfter").fail("0 (disabled) or an integer in [1,1000000]");
        if(!std::isfinite(input.run.endTime/policy.minDt) || input.run.endTime/policy.minDt>=std::numeric_limits<int>::max()-1.)
            retry->get("minDt").fail("fewer than INT_MAX accepted steps to endTime");
        input.run.stepRetry=policy;
    }
    if(const auto* adaptive=time.find("Adaptive")){
        adaptive->allow({"minDt","maxDt","targetCourant","safety","maxGrowth"});
        RunControls::Adaptive policy;
        const auto positive=[&](const char* key,const char* units){const auto& node=adaptive->get(key);const auto value=node.number();
            if(!std::isfinite(value)||!(value>0))node.fail("a finite positive adaptive parameter");
            return trace.record(node,value,units,"Runtime.VOF.Adaptive");};
        policy.minDt=positive("minDt","code_time");policy.maxDt=positive("maxDt","code_time");
        policy.targetCourant=positive("targetCourant","dimensionless");
        policy.safety=adaptive->find("safety")?positive("safety","dimensionless"):
            trace.defaultValue(*adaptive,"safety",policy.safety,"dimensionless","Runtime.VOF.Adaptive");
        policy.maxGrowth=adaptive->find("maxGrowth")?positive("maxGrowth","dimensionless"):
            trace.defaultValue(*adaptive,"maxGrowth",policy.maxGrowth,"dimensionless","Runtime.VOF.Adaptive");
        if(policy.minDt>input.run.dt)adaptive->get("minDt").fail("Adaptive.minDt <= Time.dt (first-step seed)");
        if(input.run.dt>policy.maxDt)adaptive->get("maxDt").fail("Time.dt (first-step seed) <= Adaptive.maxDt");
        if(!(policy.safety<=1))adaptive->get("safety").fail("Adaptive.safety in (0,1]");
        if(!(policy.maxGrowth>=1))adaptive->get("maxGrowth").fail("Adaptive.maxGrowth >= 1");
        if(!std::isfinite(input.run.endTime/policy.minDt)||input.run.endTime/policy.minDt>=std::numeric_limits<int>::max()-1.)
            adaptive->get("minDt").fail("fewer than INT_MAX accepted steps to endTime");
        if(input.run.stepRetry){const auto& retry=*input.run.stepRetry;
            if(retry.minDt!=policy.minDt)adaptive->get("minDt").fail("Adaptive.minDt must exactly match Time.StepRetry.minDt");
            if(retry.recoverAfter!=0)time.block("StepRetry").get("recoverAfter").fail("recoverAfter 0 with Time.Adaptive; adaptive policy owns growth");
        }
        input.run.adaptive=policy;
        // Nominal seed metadata only. The runtime uses the actual accepted
        // clock/endTime, never this count, for completion or output clocks.
        input.steps=static_cast<int>(std::max(1.,std::ceil(input.run.endTime/input.run.dt)));
    }
    const auto& output=run.block("Output");output.allow({"Fields","Series"});
    std::string seriesProfile="physical";
    if(const auto* series=output.find("Series")){
        series->allow({"profile"});
        seriesProfile=series->find("profile")
            ?trace.record(series->get("profile"),series->get("profile").choice({"physical","diagnostic"}),"dimensionless","Runtime.output.seriesProfile")
            :trace.defaultValue(*series,"profile",std::string("physical"),"dimensionless","Runtime.output.seriesProfile");
    }else trace.defaultValue(output,"Series.profile",std::string("physical"),"dimensionless","Runtime.output.seriesProfile");
    input.run.diagnosticSeries=seriesProfile=="diagnostic";
    const auto& outFields=output.block("Fields");outFields.allow({"enabled","writeEverySteps"});
    input.run.writeFields=trace.record(outFields.get("enabled"),outFields.get("enabled").boolean(),"dimensionless","Runtime.output.fields");
    input.run.writeInitial=true;input.run.writeFirst=false;input.run.writeFinal=true;
    if(input.run.writeFields){
        input.run.writeEverySteps=outFields.find("writeEverySteps")
            ?trace.record(outFields.get("writeEverySteps"),outFields.get("writeEverySteps").integer(),"steps","Runtime.output.schedule")
            :trace.defaultValue(outFields,"writeEverySteps",1,"steps","Runtime.output.schedule");
        if(input.run.writeEverySteps<1 || input.run.writeEverySteps>1000000)
            outFields.get("writeEverySteps").fail("an integer interval in [1,1000000]");
    }else if(outFields.find("writeEverySteps"))outFields.get("writeEverySteps").fail("omitted when field output is disabled");
    std::string logLevel="summary";
    if(const auto* log=run.find("Log")){
        log->allow({"everySteps","level"});
        input.run.logEverySteps=log->find("everySteps")?trace.record(log->get("everySteps"),log->get("everySteps").integer(),"steps","Runtime.console.schedule"):
            trace.defaultValue(*log,"everySteps",1,"steps","Runtime.console.schedule");
        logLevel=log->find("level")?trace.record(log->get("level"),log->get("level").choice({"summary","detail"}),"dimensionless","Runtime.console.verbosity"):
            trace.defaultValue(*log,"level",std::string("summary"),"dimensionless","Runtime.console.verbosity");
        if(input.run.logEverySteps<1 || input.run.logEverySteps>1000000)log->get("everySteps").fail("an integer interval in [1,1000000]");
    }else{
        input.run.logEverySteps=trace.defaultValue(run,"Log.everySteps",1,"steps","Runtime.console.schedule");
        trace.defaultValue(run,"Log.level",std::string("summary"),"dimensionless","Runtime.console.verbosity");
    }
    input.run.detailedLog=logLevel=="detail";
    documents.root("geometry.dict").allow({"Geometry"});const auto& geometry=documents.root("geometry.dict").block("Geometry");
    geometry.allow({"coordinates","Domain","Phases","Initial","Constants"});
    input.radial=trace.record(geometry.get("coordinates"),geometry.get("coordinates").choice({"xy","rz"}),"dimensionless","Geometry.measure")=="rz";
    const auto& domain=geometry.block("Domain");domain.allow({"min","max"});
    input.domainMin=vector(domain.get("min"),trace,"code_length","Geometry.domain");
    input.domainMax=vector(domain.get("max"),trace,"code_length","Geometry.domain");
    const auto width=input.domainMax.x-input.domainMin.x,height=input.domainMax.y-input.domainMin.y;
    if(!(std::isfinite(width) && std::isfinite(height) && width>0 && height>0) || (input.radial && input.domainMin.x!=0))
        domain.fail("finite ordered bounds; rz starts at radial minimum zero");
    const auto& phases=geometry.block("Phases");phases.allow({"negative","positive"});
    input.negativePhase=trace.record(phases.get("negative"),phases.get("negative").integer(),"1-based ID","Geometry.negativePhase");
    input.positivePhase=trace.record(phases.get("positive"),phases.get("positive").integer(),"1-based ID","Geometry.positivePhase");
    if(input.negativePhase<1 || input.positivePhase<1 || input.negativePhase==input.positivePhase)
        phases.fail("two distinct positive phase IDs; negative phi binds HYDRO liquid material");
    const auto& initial=geometry.block("Initial");
    if(csgInitial){
        initial.allow({"csgFile","lengthScaleSI","excludedPhases","Pin"});
        input.csgFile=trace.record(initial.get("csgFile"),initial.get("csgFile").word(),"local basename","Core.InitialPhaseSampler.source");
        input.csgLengthScaleSI=bounded(initial.get("lengthScaleSI"),trace,0,true,"Core.DictionaryGeometry.lengthScale","m");
        input.csgExcludedPhases=trace.record(initial.get("excludedPhases"),initial.get("excludedPhases").ids(),"1-based IDs","Core.selectFixedSolids");
        if(input.csgExcludedPhases.empty())initial.get("excludedPhases").fail("nonempty distinct excluded solid phase IDs");
        UnitsScales scales;scales.L0_m=input.csgLengthScaleSI;const auto units=UnitsContext::MakeSI(scales);
        auto definition=std::make_shared<DictionaryGeometry>();
        const auto path=documents.path(input.csgFile);const auto& source=documents.sourceText(input.csgFile);
        const auto read=[&]{std::ifstream in(path,std::ios::binary);if(!in)throw std::runtime_error("cannot read declared CSG source: "+path);
            std::ostringstream text;text<<in.rdbuf();if(in.bad())throw std::runtime_error("cannot read declared CSG bytes: "+path);return text.str();};
        if(read()!=source)initial.get("csgFile").fail("CSG source bytes unchanged since case capture");
        definition->readDict(path,&units);
        if(read()!=source)initial.get("csgFile").fail("CSG source bytes unchanged through Core parse");
        if(!definition->hasDomain||definition->domainMin.x!=input.domainMin.x||definition->domainMin.y!=input.domainMin.y||
            definition->domainMax.x!=input.domainMax.x||definition->domainMax.y!=input.domainMax.y)
            initial.get("lengthScaleSI").fail("CSG physical Domain converted by lengthScaleSI must equal the declared code Domain exactly");
        const auto phases=definition->phaseName.size();std::set<int> selected;
        if(input.negativePhase>static_cast<int>(phases)||input.positivePhase>static_cast<int>(phases))
            initial.get("excludedPhases").fail("liquid/gas phase IDs present in the Core CSG definition");
        for(int id:input.csgExcludedPhases)if(id<1||id>static_cast<int>(phases)||id==input.negativePhase||id==input.positivePhase||!selected.insert(id).second)
            initial.get("excludedPhases").fail("distinct valid solid IDs disjoint from the two liquid/gas IDs");
        if(selected.size()+2!=phases)initial.get("excludedPhases").fail("every Core CSG phase belongs to liquid, gas or explicit excluded solids");
        if(const auto* pin=initial.find("Pin")){
            pin->allow({"pointId","wetBoundaryId","dryBoundaryId","wallTrace"});HydrodynamicsInput::CsgPin value;
            value.pointId=trace.record(pin->get("pointId"),pin->get("pointId").integer(),"1-based ID","Interface.PinnedCorner.initialPhysicalIDs");
            value.wetBoundaryId=trace.record(pin->get("wetBoundaryId"),pin->get("wetBoundaryId").integer(),"1-based ID","Interface.PinnedCorner.initialPhysicalIDs");
            value.dryBoundaryId=trace.record(pin->get("dryBoundaryId"),pin->get("dryBoundaryId").integer(),"1-based ID","Interface.PinnedCorner.initialPhysicalIDs");
            const auto wallTrace=pin->find("wallTrace")
                ?trace.record(pin->get("wallTrace"),pin->get("wallTrace").choice({"strictWholeCell","prescribedWetDry"}),"dimensionless","Interface.PinnedCorner.wallTrace")
                :trace.defaultValue(*pin,"wallTrace",std::string("strictWholeCell"),"dimensionless","Interface.PinnedCorner.wallTrace");
            value.wallTrace=wallTrace=="strictWholeCell"?PinWallTraceLaw::StrictWholeCell:PinWallTraceLaw::PrescribedWetDry;
            if(value.pointId<1||value.pointId>static_cast<int>(definition->points.size())||value.wetBoundaryId==value.dryBoundaryId||
                !definition->boundaryById.count(value.wetBoundaryId)||!definition->boundaryById.count(value.dryBoundaryId))
                pin->fail("actual distinct Core CSG point and wet/dry boundary IDs");
            input.csgPin=value;
        }
        input.csgGeometry=std::move(definition);
        input.documents.push_back({path,source,source});
    }else{
        initial.allow({"phi"});
        input.initialPhi=trace.record(initial.get("phi"),initial.get("phi").word(),"code_phi",vof?"Interface.GeometricVofInitialization.initialData":"Interface.FirstJet.initialData");
    }
    documents.root("mesh.dict").allow({"Mesh"});const auto& mesh=documents.root("mesh.dict").block("Mesh");if(nativeGeometry)mesh.allow({"resolution","initialMesh","Refinement","Regrid"});else mesh.allow({"resolution","initialMesh"});
    const auto cells=resolution(mesh.get("resolution"),trace,nativeGeometry?"Mesh.baseGrid":"Mesh.uniformGrid");input.nx=cells[0];input.ny=cells[1];
    if(!nativeGeometry){
    if(mesh.find("initialMesh"))trace.record(mesh.get("initialMesh"),mesh.get("initialMesh").choice({"uniform"}),"dimensionless","Mesh.initialGrid");
    else trace.defaultValue(mesh,"initialMesh",std::string("uniform"),"dimensionless","Mesh.initialGrid");
    }
    if(geometry.find("Constants"))for(const auto& node:geometry.block("Constants").entries){
        if(std::set<std::string>{"x","y","r","z","t","s","pi","sin","cos","tan","exp","log","sqrt"}.count(node.name))
            node.fail("a nonreserved expression constant name");
        const auto value=node.number();
        try{Expression test("0",{{node.name,value}},Expression::Coordinates::Time);}catch(const std::exception& e){node.fail(e.what());}
        input.constants.emplace(node.name,trace.record(node,value,"code","Geometry.expressionParameters"));
    }
    if(!csgInitial)try{
        input.expression_=std::make_shared<const Expression>(input.initialPhi,input.constants,input.radial?Expression::Coordinates::Meridional:Expression::Coordinates::Planar);
        for(const auto x:{input.domainMin.x,input.domainMax.x})for(const auto y:{input.domainMin.y,input.domainMax.y})input.initialJet({x,y});
    }catch(const std::exception& e){initial.get("phi").fail(std::string(vof?"a valid spatial VOF initial expression without t: ":"a valid spatial FirstJet expression without t: ")+e.what());}
    if(!csgInitial)usage.include(input.initialPhi);
    std::string scalarRepresentation,gradientUpdate,boundaryType;std::vector<std::string> symmetry;const Node* gradient=nullptr;
    if(!nativeGeometry){
    interface.allow({"method","scalarRepresentation","gradientUpdate","Boundary"});

    const auto* scalar=interface.find("scalarRepresentation");
    scalarRepresentation=scalar
        ?trace.record(*scalar,scalar->choice({"singleQuartic","c2MacroQuartic"}),"dimensionless","Interface.FirstJet.scalarRepresentation")
        :trace.defaultValue(interface,"scalarRepresentation",std::string("singleQuartic"),"dimensionless","Interface.FirstJet.scalarRepresentation");
    input.step.publication.scalarRepresentation=scalarRepresentation=="c2MacroQuartic"
        ?LAMP::Interface::HermiteGalsScalarRepresentation::C2MacroQuartic
        :LAMP::Interface::HermiteGalsScalarRepresentation::SingleQuartic;
    if(scalarRepresentation=="c2MacroQuartic"&&input.radial)
        scalar->fail("c2MacroQuartic requires the complete fixed uniform XY grid; RZ is unsupported");
    gradient=interface.find("gradientUpdate");
    gradientUpdate=gradient
        ?trace.record(*gradient,gradient->choice({"characteristicPullback","compatibleIncrement"}),"dimensionless","Interface.FirstJet.gradientUpdate")
        :trace.defaultValue(interface,"gradientUpdate",std::string("characteristicPullback"),"dimensionless","Interface.FirstJet.gradientUpdate");
    input.step.publication.gradientUpdate=gradientUpdate=="compatibleIncrement"
        ?LAMP::Interface::HermiteGalsGradientUpdate::CompatibleIncrement:LAMP::Interface::HermiteGalsGradientUpdate::CharacteristicPullback;
    const auto& interfaceBoundary=interface.block("Boundary");
    boundaryType=trace.record(interfaceBoundary.get("type"),interfaceBoundary.get("type").choice({"none","symmetry"}),"dimensionless","Interface.FirstJet.boundary");

    if(boundaryType=="symmetry"){
        interfaceBoundary.allow({"type","sides"});symmetry=sides(interfaceBoundary.get("sides"),trace);
        for(const auto& side:symmetry){
            const int axis=(side=="left" || side=="right")?0:1;
            const auto coordinate=side=="left"?input.domainMin.x:side=="right"?input.domainMax.x:side=="bottom"?input.domainMin.y:input.domainMax.y;
            input.step.publication.volume.symmetryPlanes.push_back({axis,coordinate});
        }
        input.step.spaceTime.firstJetBoundaryTreatment=LAMP::Interface::HermiteGalsBoundaryTreatment::SymmetryFirstJet;
    }else{
        interfaceBoundary.allow({"type"});
        input.step.spaceTime.firstJetBoundaryTreatment=LAMP::Interface::HermiteGalsBoundaryTreatment::InflowHistory;
    }
    }
    documents.root("materials.dict").allow({"Materials"});const auto& materials=documents.root("materials.dict").block("Materials");materials.allow({"Phases","Interfaces"});
    const auto& materialPhases=materials.block("Phases");
    if(materialPhases.entries.size()!=2)materialPhases.fail("exactly the two Geometry.Phases materials");
    std::set<int> materialIds;
    for(const auto& phase:materialPhases.entries){
        if(!phase.isBlock)phase.fail("a named material block");
        const auto id=trace.record(phase.get("phaseId"),phase.get("phaseId").integer(),"1-based ID","Materials.geometryBinding");
        if((id!=input.negativePhase && id!=input.positivePhase) || !materialIds.insert(id).second)
            phase.get("phaseId").fail("one distinct Geometry.Phases ID");
        const auto rho=bounded(phase.get("density"),trace,0,true,"Hydrodynamics.geometryMass","code_density");
        const auto mu=bounded(phase.get("dynamicViscosity"),trace,0,false,"Hydrodynamics.phaseViscosity","code_dynamic_viscosity");
        if(id==input.negativePhase){if(vof){input.vof.liquidDensity=rho;input.vof.liquidViscosity=mu;}else{input.step.liquidDensity=rho;input.step.liquidViscosity=mu;}}
        else{if(vof){input.vof.gasDensity=rho;input.vof.gasViscosity=mu;}else{input.step.gasDensity=rho;input.step.gasViscosity=mu;}}
        auto& material=components.materials[materialIds.size()-1];material.phaseId=id;
        auto prefix=stream();prefix<<phase.name<<" { phaseId "<<id<<"; density "<<rho
            <<"; dynamicViscosity "<<mu<<';';material.prefix=prefix.str();

    }
    const auto& pairs=materials.block("Interfaces");
    if(pairs.entries.size()!=1)pairs.fail("exactly one interface between the two Geometry.Phases");
    const auto& pair=pairs.entries.front();if(!pair.isBlock)pair.fail("a named material interface block");pair.allow({"phases","surfaceTension"});
    const auto ids=trace.record(pair.get("phases"),pair.get("phases").ids(),"1-based IDs","Materials.interfaceBinding");
    if(ids.size()!=2 || !materialIds.count(ids[0]) || !materialIds.count(ids[1]))pair.get("phases").fail("the two distinct Geometry.Phases IDs");
    const auto surfaceTension=bounded(pair.get("surfaceTension"),trace,0,false,"Hydrodynamics.capillaryLoad","code_surface_tension");if(vof)input.vof.surfaceTension=surfaceTension;else input.step.surfaceTension=surfaceTension;
    documents.root("hydro.dict").allow({"Hydrodynamics"});const auto& hydro=documents.root("hydro.dict").block("Hydrodynamics");
    if(vof)hydro.allow({"method","capillaryModel","densityRepresentation","Initial","Boundary","BodyForce","TimeIntegration","Numerics"});else hydro.allow({"method","capillaryModel","densityRepresentation","Initial","Boundary","BodyForce","TimeIntegration","Numerics","CSF"});
    trace.record(hydro.get("method"),hydro.get("method").choice({"staggered"}),"dimensionless","Hydrodynamics.dispatch");
    const auto& hydroInitial=hydro.block("Initial");hydroInitial.allow({"velocity"});
    input.initialVelocity=vector(hydroInitial.get("velocity"),trace,"code_velocity","Hydrodynamics.nativeInitialVelocity");
    const auto& hydroBoundary=hydro.block("Boundary");
    const bool allowOpen=vof&&ehdFamily;
    const auto requestedBoundary=hydroBoundary.get("type").word();
    const bool requestedOpen=allowOpen&&(requestedBoundary=="prescribedFlow"||requestedBoundary=="pressureDrivenReservoirs");
    if(hydroBoundary.find("Backflow")&&!requestedOpen)
        hydroBoundary.block("Backflow").fail("Backflow is consumed only by VOF-EHD prescribedFlow or pressureDrivenReservoirs");
    // The old zero initial field remains the closed initial state. For open
    // flow it is only a seed for Run's actual constrained initialization.
    if(input.initialVelocity.x!=0 || input.initialVelocity.y!=0)hydroInitial.get("velocity").fail(requestedOpen?
        "constant (0 0) seed for actual native-mass open initialization":"constant (0 0) for the declared fully impermeable domain");
    const auto physicalBoundary=trace.record(hydroBoundary.get("type"),allowOpen?
        hydroBoundary.get("type").choice({"impermeable","prescribedFlow","pressureDrivenReservoirs"}):hydroBoundary.get("type").choice({"impermeable"}),
        "dimensionless","Hydrodynamics.nativeBoundary");
    input.vofPressureDrivenReservoirs=physicalBoundary=="pressureDrivenReservoirs";
    if(physicalBoundary=="prescribedFlow"||input.vofPressureDrivenReservoirs){
        hydroBoundary.allow({"type","Patches","Backflow","StationaryNoSlip"});
        readVofOpenPatches(input,hydroBoundary,trace);
        if(input.vofPressureDrivenReservoirs&&!hydroBoundary.find("Backflow"))
            hydroBoundary.fail("explicit Backflow model required for pressureDrivenReservoirs");
        readVofBackflow(input,hydroBoundary,trace);
        if(const auto* walls=hydroBoundary.find("StationaryNoSlip")){
            walls->allow({"boundaryIds"});
            if(!csgInitial||!input.vofPressureDrivenReservoirs)
                walls->fail("StationaryNoSlip only for actual A CSG pressureDrivenReservoirs");
            const auto& ids=walls->get("boundaryIds");
            input.vofStationaryNoSlipBoundaryIds=trace.record(ids,ids.ids(),"1-based Core boundary IDs","HydroPhaseViscosity.stationaryNoSlipPhysicalWalls");
            if(input.vofStationaryNoSlipBoundaryIds.empty())ids.fail("nonempty distinct declared physical wall IDs");
            std::set<int> selected;
            for(const int id:input.vofStationaryNoSlipBoundaryIds){const auto found=input.csgGeometry->boundaryById.find(id);
                if(id<1||found==input.csgGeometry->boundaryById.end()||found->second.ignored||!selected.insert(id).second)
                    ids.fail("distinct declared nonignored Core CSG wall boundary IDs");
                for(const auto& port:input.vofBoundaryPatches)if(port.boundaryId&&*port.boundaryId==id)
                    ids.fail("disjoint physical reservoir and stationary no-slip boundary IDs");
            }
            if(input.vof.liquidViscosity==0&&input.vof.gasViscosity==0)
                walls->fail("stationary full-trace no-slip requires positive physical viscosity");
        }
    }else hydroBoundary.allow({"type"});
    const auto& force=hydro.block("BodyForce");force.allow({"acceleration"});
    const auto acceleration=vector(force.get("acceleration"),trace,"code_length/code_time^2","Hydrodynamics.midpointMassBodyForce");if(vof)input.vof.bodyAcceleration=acceleration;else input.step.bodyAcceleration=acceleration;
    if(input.radial && acceleration.x!=0)force.get("acceleration").fail("zero radial component for constant axisymmetric acceleration");
    std::array<std::string,3> nativeEffective;std::string densityRepresentation="sharpGeometry",capillaryModel="balancedCsf",csfForceMethod="motionAdjoint";
    if(vof){
        nativeEffective=readVofMethodInput(input,mesh,interface,hydro);
        if(input.initialInterfaceBand&&input.initialBandScale==HydrodynamicsInput::InitialBandScale::TargetCell&&!csgInitial)
            mesh.get("Refinement").fail("bandScale targetCell requires the explicit immutable CSG initial source");
        if(csgInitial){
            if(input.initialInterfaceBand&&input.initialBandScale!=HydrodynamicsInput::InitialBandScale::TargetCell)
                mesh.get("Refinement").fail("CSG initialInterfaceBand requires explicit bandScale targetCell");
            if(input.vofRegrid){const auto& policy=*input.vofRegrid;
                if(policy.liquidMinLevel.value_or(policy.minLevel)!=policy.minLevel||policy.gasMinLevel.value_or(policy.minLevel)!=policy.minLevel)
                    mesh.block("Regrid").fail("CSG initialization currently requires phase floors equal to the common minLevel");
            }
            if(!input.vofBoundaryPatches.empty()){
                if(!input.vofPressureDrivenReservoirs)hydroBoundary.get("type").fail("CSG open flow requires pressureDrivenReservoirs without fixed velocity feed");
                if(input.vofStationaryNoSlipBoundaryIds.empty())hydroBoundary.fail("explicit StationaryNoSlip boundaryIds for actual CSG fixed-solid walls");
            }
            if(surfaceTension>0&&!input.csgPin)pair.get("surfaceTension").fail("positive CSG surfaceTension requires its explicit stock-derived corner Pin");
            if(!input.vofSymmetryPlanes.empty()||input.vofSouthEven90Interval)interface.block("Boundary").fail("unsupported additional CSG interface boundary provider");
            if(ehdFamily){
                const auto& electric=documents.root("electric.dict").block("Electric");
                if(electric.get("method").word()!="vofFiniteVolume")electric.get("method").fail("CSG initialization requires the actual A finite-volume charge route");
                if(const auto* type=electric.block("Initial").find("type"))if(type->word()!="cellDensity")type->fail("CSG initialization supports cellDensity charge only");
            }
            // Pure geometry configuration: the sampler owns its original
            // callbacks/AST after this short-lived base probe is destroyed.
            DictionaryMesh probeInput;probeInput.minCoord=input.domainMin;probeInput.maxCoord=input.domainMax;
            probeInput.Nx=input.nx;probeInput.Ny=input.ny;probeInput.maxInitialIterations=0;probeInput.levelSetCount=0;
            probeInput.maxCellNum=std::max(16384,input.nx*input.ny);QuadtreeMesh probe(probeInput);
            input.csgSampler=std::make_shared<const FittedDomainGeometry::InitialPhaseSampler>(probe,*input.csgGeometry);
            (void)FittedDomainGeometry::selectFixedSolids(probe,*input.csgGeometry,*input.csgSampler,input.csgExcludedPhases);
        }
        if(input.vofSouthEven90Interval){
            const auto& declaration=interface.block("Boundary");
            if(!ehdFamily)declaration.get("type").fail("southEven90 only for the actual VOF-EHD stage provider");
            const LAMP::Mesh::DyadicGridCoordinates coordinates(input.domainMin,input.domainMax,input.nx,input.ny);
            const auto aligned=[&](double x){int left=0,right=input.nx;
                while(left<=right){const int k=left+(right-left)/2;const auto actual=coordinates.point(k,0,0).x;
                    if(actual==x)return true;if(actual<x)left=k+1;else right=k-1;}return false;};
            const auto interval=*input.vofSouthEven90Interval;
            if(!aligned(interval.x)||!aligned(interval.y))declaration.get("interval").fail("physical S endpoints exactly on actual base-grid coordinates");
        }
        if(requestedOpen&&surfaceTension>0&&!input.vofSouthEven90Interval&&!(csgInitial&&input.csgPin))
            pair.get("surfaceTension").fail("positive surfaceTension on open VOF-EHD requires explicit southEven90 or stock-derived CSG Pin geometry");

        if(input.run.adaptive&&input.run.adaptive->targetCourant>input.vof.transport.maximumLocalCourant)
            time.block("Adaptive").get("targetCourant").fail("Adaptive.targetCourant <= Hydrodynamics.Numerics.maximumLocalCourant");
    }else{

    if(const auto* density=hydro.find("densityRepresentation")){
        densityRepresentation=density->choice({"sharpGeometry","cellAverageProjection"});
        trace.record(*density,densityRepresentation,"dimensionless","HydroTwoPhaseStep.momentumDensityRepresentation");
    }else trace.defaultValue(hydro,"densityRepresentation",densityRepresentation,"dimensionless","HydroTwoPhaseStep.momentumDensityRepresentation");
    input.step.momentumDensityRepresentation=densityRepresentation=="sharpGeometry"
        ?HydroMomentumMass::DensityRepresentation::SharpGeometry:HydroMomentumMass::DensityRepresentation::CellAverageProjection;

    if(const auto* model=hydro.find("capillaryModel")){
        capillaryModel=model->choice({"balancedCsf","sharpNormalStress"});
        trace.record(*model,capillaryModel,"dimensionless","HydroTwoPhaseStep.capillaryModel");
    }else trace.defaultValue(hydro,"capillaryModel",capillaryModel,"dimensionless","HydroTwoPhaseStep.capillaryModel");
    input.step.capillaryModel=capillaryModel=="balancedCsf"?HydroCapillaryModel::BalancedCsf:HydroCapillaryModel::SharpNormalStress;
    const auto* csf=hydro.find("CSF");
    if(input.step.capillaryModel==HydroCapillaryModel::SharpNormalStress&&csf)
        csf->fail("no CSF block for capillaryModel=sharpNormalStress; its parameters would be unused");

    if(input.step.capillaryModel==HydroCapillaryModel::BalancedCsf){
    if(csf){
        csf->allow({"forceMethod","curvatureTolerance"});
        if(const auto* method=csf->find("forceMethod")){
            csfForceMethod=method->choice({"motionAdjoint","heightFunction","heightFunctionIncrement"});
            trace.record(*method,csfForceMethod,"dimensionless","Hydrodynamics.CSF.forceMethod");
        }else trace.defaultValue(*csf,"forceMethod",csfForceMethod,"dimensionless","Hydrodynamics.CSF.forceMethod");
    }else trace.defaultValue(hydro,"CSF.forceMethod",csfForceMethod,"dimensionless","Hydrodynamics.CSF.forceMethod");
    input.step.csfCurvature.forceMethod=csfForceMethod=="heightFunctionIncrement"
        ?HydroCsfForceMethod::HeightFunctionIncrement:csfForceMethod=="heightFunction"
        ?HydroCsfForceMethod::HeightFunction:HydroCsfForceMethod::MotionAdjoint;
    if(input.step.csfCurvature.forceMethod!=HydroCsfForceMethod::MotionAdjoint){
        const auto& method=csf->get("forceMethod");
        if(csfForceMethod=="heightFunctionIncrement"&&gradientUpdate!="compatibleIncrement")
            method.fail("heightFunctionIncrement requires Interface.gradientUpdate=compatibleIncrement");
        if(input.radial)method.fail("heightFunction for complete uniform XY; axisymmetric support is not implemented");
        if(input.step.publication.scalarRepresentation!=LAMP::Interface::HermiteGalsScalarRepresentation::SingleQuartic)
            method.fail("heightFunction with Interface.scalarRepresentation=singleQuartic");
        if(const auto* unused=csf->find("curvatureTolerance"))
            unused->fail("no curvatureTolerance for heightFunction; this surface quadrature setting would be unused");
    }
    if(csf && input.step.csfCurvature.forceMethod==HydroCsfForceMethod::MotionAdjoint){
        input.step.csfCurvature.relativeQuadratureTolerance=option(*csf,"curvatureTolerance",
            input.step.csfCurvature.relativeQuadratureTolerance,trace,"dimensionless","HydroTwoPhaseStep.csfCurvature");
    }else if(input.step.csfCurvature.forceMethod==HydroCsfForceMethod::MotionAdjoint)
        trace.defaultValue(hydro,"CSF.curvatureTolerance",input.step.csfCurvature.relativeQuadratureTolerance,
            "dimensionless","HydroTwoPhaseStep.csfCurvature");
    if(!(std::isfinite(input.step.csfCurvature.relativeQuadratureTolerance)
        &&input.step.csfCurvature.relativeQuadratureTolerance>0&&input.step.csfCurvature.relativeQuadratureTolerance<1)){
        if(csf&&csf->find("curvatureTolerance"))csf->get("curvatureTolerance").fail("a finite CSF relative tolerance strictly between zero and one");
        (csf?*csf:hydro).fail("a finite CSF relative tolerance strictly between zero and one");
    }
    }
    if(gradientUpdate=="compatibleIncrement"){
        if(input.radial||scalarRepresentation!="singleQuartic"
            ||input.step.capillaryModel!=HydroCapillaryModel::BalancedCsf
            ||(input.step.csfCurvature.forceMethod!=HydroCsfForceMethod::HeightFunction
                &&input.step.csfCurvature.forceMethod!=HydroCsfForceMethod::HeightFunctionIncrement))
            gradient->fail("compatibleIncrement requires standalone uniform XY / singleQuartic / balancedCsf / heightFunction or heightFunctionIncrement");
    }
    const auto& integration=hydro.block("TimeIntegration");
    integration.allow({"method","interfaceMismatchTolerance","maximumCorrectionRelative","correctionReferenceVelocity",
        "capillaryCorrectionCourant","maximumCapillaryCourant"});
    trace.record(integration.get("method"),integration.get("method").choice({"adaptiveLinear"}),"dimensionless","Hydrodynamics.AC.timeIntegration");
    input.step.timeScheme=HydroTwoPhaseTimeScheme::AdaptiveLinearMidpoint;
    input.step.interfaceMismatchTolerance=option(integration,"interfaceMismatchTolerance",input.step.interfaceMismatchTolerance,trace,"dimensionless","Hydrodynamics.AC.correctionTrigger");
    input.step.maximumCorrectionRelative=option(integration,"maximumCorrectionRelative",input.step.maximumCorrectionRelative,trace,"dimensionless","Hydrodynamics.AC.correctionGate");
    input.step.correctionReferenceVelocity=option(integration,"correctionReferenceVelocity",input.step.correctionReferenceVelocity,trace,"code_velocity","Hydrodynamics.AC.correctionScale");
    input.step.capillaryCorrectionCourant=option(integration,"capillaryCorrectionCourant",input.step.capillaryCorrectionCourant,trace,"dimensionless","Hydrodynamics.AC.capillaryCorrectionTrigger");
    input.step.maximumCapillaryCourant=option(integration,"maximumCapillaryCourant",input.step.maximumCapillaryCourant,trace,"dimensionless","Hydrodynamics.AC.capillaryRetryBudget");
    for(const char* name:{"capillaryCorrectionCourant","maximumCapillaryCourant"}){
        const double value=std::string(name)=="capillaryCorrectionCourant"?input.step.capillaryCorrectionCourant:input.step.maximumCapillaryCourant;
        if(!(value>0))integration.get(name).fail("a finite positive capillary Courant parameter");
    }
    if(input.step.capillaryCorrectionCourant>input.step.maximumCapillaryCourant)
        integration.fail("capillaryCorrectionCourant <= maximumCapillaryCourant");
    for(const char* name:{"interfaceMismatchTolerance","maximumCorrectionRelative","correctionReferenceVelocity"}){
        const auto value=std::string(name)=="interfaceMismatchTolerance"?input.step.interfaceMismatchTolerance:
            std::string(name)=="maximumCorrectionRelative"?input.step.maximumCorrectionRelative:input.step.correctionReferenceVelocity;
        if(!(value>0))integration.get(name).fail("a finite positive AC parameter");
    }
    if(const auto* numerics=hydro.find("Numerics")){
        numerics->allow({"characteristicSubsteps","materialInitialPanels"});
        input.characteristicSubsteps=numerics->find("characteristicSubsteps")
            ?trace.record(numerics->get("characteristicSubsteps"),numerics->get("characteristicSubsteps").integer(),"substeps","Hydrodynamics.App.actualTraceBudget")
            :trace.defaultValue(*numerics,"characteristicSubsteps",1,"substeps","Hydrodynamics.App.actualTraceBudget");
        input.step.spaceTime.initialPanels=numerics->find("materialInitialPanels")
            ?trace.record(numerics->get("materialInitialPanels"),numerics->get("materialInitialPanels").integer(),"panels","HydroMaterialTransport.prepareSpaceTime")
            :trace.defaultValue(*numerics,"materialInitialPanels",input.step.spaceTime.initialPanels,"panels","HydroMaterialTransport.prepareSpaceTime");
        if(input.characteristicSubsteps<1 || input.characteristicSubsteps>1024)
            numerics->get("characteristicSubsteps").fail("an integer in [1,1024] for the actual attempt dt/substeps trace budget");
        if(input.step.spaceTime.initialPanels<1 || input.step.spaceTime.initialPanels>1024)
            numerics->get("materialInitialPanels").fail("an integer initial quadrature count in [1,1024]");
    }else{
        input.characteristicSubsteps=trace.defaultValue(hydro,"Numerics.characteristicSubsteps",1,"substeps","Hydrodynamics.App.actualTraceBudget");
        trace.defaultValue(hydro,"Numerics.materialInitialPanels",input.step.spaceTime.initialPanels,"panels","HydroMaterialTransport.prepareSpaceTime");
    }
    if(!(input.run.dt/input.characteristicSubsteps>0))
        time.get("dt").fail("a positive resolved nominal characteristic dt/substeps budget");
    if(input.run.stepRetry && !(input.run.stepRetry->minDt/input.characteristicSubsteps>0))
        time.block("StepRetry").get("minDt").fail("a positive resolved minimum characteristic dt/substeps budget");
    }
    auto runText=stream(),geoText=stream(),meshText=stream(),interfaceText=stream(),materialText=stream(),hydroText=stream();
    runText<<"Run { app "<<app<<"; }\nParallel { threads "<<input.parallelThreads<<"; }\nTime { units code; dt "<<input.run.dt<<"; endTime "<<input.run.endTime<<';';
    if(input.run.stepRetry){const auto& retry=*input.run.stepRetry;runText<<" StepRetry { minDt "<<retry.minDt<<"; maxRetries "<<retry.maxRetries<<"; recoverAfter "<<retry.recoverAfter<<"; }";}
    if(input.run.adaptive){const auto& adaptive=*input.run.adaptive;
        runText<<" Adaptive { minDt "<<adaptive.minDt<<"; maxDt "<<adaptive.maxDt<<"; targetCourant "<<adaptive.targetCourant
            <<"; safety "<<adaptive.safety<<"; maxGrowth "<<adaptive.maxGrowth<<"; }";}
    runText<<" }\nOutput { Fields { enabled "<<(input.run.writeFields?"true":"false")<<';';
    if(input.run.writeFields)runText<<" writeEverySteps "<<input.run.writeEverySteps<<';';
    runText<<" } Series { profile "<<seriesProfile<<"; } }\nLog { everySteps "<<input.run.logEverySteps<<"; level "<<logLevel<<"; }\n";
    if(vof){
        runText<<"Checkpoint { everySteps "<<input.run.checkpointEverySteps<<"; }\n";
        if(!input.run.restartDirectory.empty())runText<<"Restart { directory \""<<input.run.restartDirectory<<"\"; }\n";
    }
    geoText<<"Geometry {\ncoordinates "<<(input.radial?"rz":"xy")<<";\nDomain { min ";vec(geoText,input.domainMin);geoText<<"; max ";vec(geoText,input.domainMax);
    geoText<<"; }\nPhases { negative "<<input.negativePhase<<"; positive "<<input.positivePhase<<"; }\n";
    if(!input.constants.empty()){geoText<<"Constants {\n";for(const auto& item:input.constants)geoText<<item.first<<' '<<item.second<<";\n";geoText<<"}\n";}
    if(csgInitial){
        geoText<<"Initial { csgFile "<<std::quoted(input.csgFile)<<"; lengthScaleSI "<<input.csgLengthScaleSI<<"; excludedPhases (";
        for(int id:input.csgExcludedPhases)geoText<<id<<' ';geoText<<");";
        if(input.csgPin)geoText<<" Pin { pointId "<<input.csgPin->pointId<<"; wetBoundaryId "<<input.csgPin->wetBoundaryId
            <<"; dryBoundaryId "<<input.csgPin->dryBoundaryId<<"; wallTrace "<<pinWallTraceName(input.csgPin->wallTrace)<<"; }";
        geoText<<" }\n}\n";
    }else geoText<<"Initial { phi "<<std::quoted(input.initialPhi)<<"; }\n}\n";
    if(nativeGeometry)meshText<<nativeEffective[0];else meshText<<"Mesh { resolution ("<<input.nx<<' '<<input.ny<<"); initialMesh uniform; }\n";
    if(nativeGeometry)interfaceText<<nativeEffective[1];else{
    interfaceText<<"Interface { method hermiteGals; scalarRepresentation "<<scalarRepresentation<<"; gradientUpdate "<<gradientUpdate<<"; Boundary { type "<<boundaryType<<';';
    if(!symmetry.empty()){interfaceText<<" sides (";for(const auto& side:symmetry)interfaceText<<side<<' ';interfaceText<<");";}
    interfaceText<<" } }\n";
    }
    materialText<<"Materials {\nPhases {\n";
    for(const auto& material:components.materials)materialText<<material.prefix<<" }\n";
    auto interfacesText=stream();interfacesText<<"}\nInterfaces { "<<pair.name<<" { phases ("<<ids[0]<<' '<<ids[1]
        <<"); surfaceTension "<<surfaceTension<<"; } }\n}\n";
    components.materialInterfaces=interfacesText.str();materialText<<components.materialInterfaces;
    if(nativeGeometry)hydroText<<nativeEffective[2];else{
    hydroText<<"Hydrodynamics {\nmethod staggered;\ncapillaryModel "<<capillaryModel<<";\ndensityRepresentation "<<densityRepresentation<<";\nInitial { velocity ";vec(hydroText,input.initialVelocity);
    hydroText<<"; }\nBoundary { type impermeable; }\nBodyForce { acceleration ";vec(hydroText,input.step.bodyAcceleration);
    hydroText<<"; }\nTimeIntegration { method adaptiveLinear; interfaceMismatchTolerance "<<input.step.interfaceMismatchTolerance
        <<"; maximumCorrectionRelative "<<input.step.maximumCorrectionRelative<<"; correctionReferenceVelocity "<<input.step.correctionReferenceVelocity
        <<"; capillaryCorrectionCourant "<<input.step.capillaryCorrectionCourant<<"; maximumCapillaryCourant "<<input.step.maximumCapillaryCourant<<"; }\n"
        ;
    if(input.step.capillaryModel==HydroCapillaryModel::BalancedCsf){
        hydroText<<"CSF { forceMethod "<<csfForceMethod<<';';
        if(input.step.csfCurvature.forceMethod==HydroCsfForceMethod::MotionAdjoint)
            hydroText<<" curvatureTolerance "<<input.step.csfCurvature.relativeQuadratureTolerance<<';';
        hydroText<<" }\n";
    }
    hydroText<<"Numerics { characteristicSubsteps "<<input.characteristicSubsteps<<"; materialInitialPanels "<<input.step.spaceTime.initialPanels<<"; }\n}\n";
    }
    components.effective={runText.str(),geoText.str(),meshText.str(),interfaceText.str(),materialText.str(),hydroText.str()};
    return components;
}
void recordTwoPhaseInputDefaults(const CaseDocuments& documents,const HydrodynamicsInput& input,
    DictionaryInput::ConfigurationTrace& trace)
{const auto& hydro=documents.root("hydro.dict").block("Hydrodynamics");if(input.interfaceMethod==HydrodynamicsInput::InterfaceMethod::Vof)recordVofInputDefaults(hydro,input,trace);else kernelDefaults(hydro,trace,input.step);}
double boundedCaseValue(const Node& node,Trace& trace,double lower,bool strictLower,const char* consumer,const char* units)
{return bounded(node,trace,lower,strictLower,consumer,units);}
}



