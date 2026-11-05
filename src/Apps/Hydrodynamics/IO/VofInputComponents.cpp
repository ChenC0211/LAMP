#include "IO/TwoPhaseInputComponents.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace AppIO
{
namespace
{
using Node=DictionaryInput::Node;using Trace=DictionaryInput::ConfigurationTrace;
using namespace LAMP::Hydrodynamics;
std::ostringstream stream(){std::ostringstream out;out<<std::setprecision(17);return out;}
void vectorText(std::ostream& out,Vector2D p){out<<'('<<p.x<<' '<<p.y<<')';}
double number(const Node& parent,const Node* block,const char* blockName,const char* key,double value,
    Trace& trace,const char* units,const char* consumer)
{
    if(block){if(const auto* node=block->find(key))return trace.record(*node,node->number(),units,consumer);
        return trace.defaultValue(*block,key,value,units,consumer);}
    return trace.defaultValue(parent,std::string(blockName)+'.'+key,value,units,consumer);
}
int integer(const Node& parent,const Node* block,const char* blockName,const char* key,int value,
    Trace& trace,const char* units,const char* consumer)
{
    if(block){if(const auto* node=block->find(key))return trace.record(*node,node->integer(),units,consumer);
        return trace.defaultValue(*block,key,value,units,consumer);}
    return trace.defaultValue(parent,std::string(blockName)+'.'+key,value,units,consumer);
}
const Node& where(const Node& parent,const Node* block,const char* key)
{return block?(block->find(key)?block->get(key):*block):parent;}
}
std::array<std::string,3> readVofMethodInput(HydrodynamicsInput& input,const Node& mesh,const Node& interface,const Node& hydro)
{
    auto& trace=input.trace;auto& config=input.vof;
    interface.allow({"method","Boundary","Reconstruction","Initialization"});
    const auto& boundary=interface.block("Boundary");
    const auto boundaryType=boundary.get("type").choice({"none","symmetry","southEven90"});
    trace.record(boundary.get("type"),boundaryType,"dimensionless",boundaryType=="southEven90"
        ?"HydroVofStep.southEven90Boundary.actualStageGeometry":"Interface.VofInitialization.boundary");
    std::vector<std::string> symmetry;
    if(boundaryType=="symmetry"){
        boundary.allow({"type","sides"});const auto& sides=boundary.get("sides");
        if(sides.isBlock||sides.values.size()<3||sides.values.front()!="("||sides.values.back()!=")")
            sides.fail("a nonempty list of distinct domain sides (left right bottom top)");
        for(std::size_t i=1;i+1<sides.values.size();++i){const auto& side=sides.values[i];
            if((side!="left"&&side!="right"&&side!="bottom"&&side!="top")||
                std::find(symmetry.begin(),symmetry.end(),side)!=symmetry.end())
                sides.fail("a nonempty list of distinct domain sides (left right bottom top)");
            symmetry.push_back(side);const int axis=(side=="left"||side=="right")?0:1;
            const auto coordinate=side=="left"?input.domainMin.x:side=="right"?input.domainMax.x:
                side=="bottom"?input.domainMin.y:input.domainMax.y;
            input.vofSymmetryPlanes.push_back({axis,coordinate});
        }
        trace.record(sides,sides.raw(),"dimensionless","Interface.VofInitialization.symmetryPlanes");
        if(config.surfaceTension>0)boundary.fail("VOF symmetry contact currently requires zero surface tension; height-function wall contact is unsupported");
    }else if(boundaryType=="southEven90"){
        boundary.allow({"type","interval"});const auto& node=boundary.get("interval");
        if(node.isBlock||node.values.size()!=4||node.values.front()!="("||node.values.back()!=")")
            node.fail("one finite physical S tangent interval (lower upper)");
        auto item=node;item.values={node.values[1]};const double lower=item.number();item.values={node.values[2]};const double upper=item.number();
        if(!(lower>=input.domainMin.x&&lower<upper&&upper<=input.domainMax.x))
            node.fail("ordered in-domain physical S tangent endpoints");
        input.vofSouthEven90Interval=Vector2D{lower,upper};auto value=stream();vectorText(value,*input.vofSouthEven90Interval);
        trace.record(node,value.str(),"code_length","HydroVofStep.southEven90Boundary.actualStageGeometry");
    }else boundary.allow({"type"});
    if(const auto* regrid=mesh.find("Regrid")){
        regrid->allow({"minLevel","maxLevel","interfaceBandLayers","vorticityTolerance","coarsenRatio","liquidMinLevel","gasMinLevel","interval"});
        HydrodynamicsVofRegridPolicy policy;
        const char* consumer="Hydrodynamics.VOF.interfaceFlowRegrid";
        policy.minLevel=integer(mesh,regrid,"Regrid","minLevel",policy.minLevel,trace,"levels",consumer);
        policy.maxLevel=integer(mesh,regrid,"Regrid","maxLevel",policy.maxLevel,trace,"levels",consumer);
        policy.interfaceBandLayers=integer(mesh,regrid,"Regrid","interfaceBandLayers",policy.interfaceBandLayers,trace,"native face layers",consumer);
        policy.vorticityTolerance=number(mesh,regrid,"Regrid","vorticityTolerance",policy.vorticityTolerance,trace,"code_velocity",consumer);
        policy.coarsenRatio=number(mesh,regrid,"Regrid","coarsenRatio",policy.coarsenRatio,trace,"dimensionless",consumer);
        if(const auto* interval=regrid->find("interval")){
            const int value=interval->integer();
            if(value<=0)interval->fail("a positive integer accepted-step interval");
            policy.interval=trace.record(*interval,value,"accepted physical steps",consumer);
        }
        if(policy.minLevel<0||policy.minLevel>policy.maxLevel||policy.maxLevel<1||policy.maxLevel>6)
            regrid->fail("0 <= minLevel <= maxLevel, maxLevel in [1,6]");
        const auto phaseMinimum=[&](const char* key)->std::optional<int>{
            // Omitted phase bounds inherit the already traced global bound.
            // Keep the legacy five-field Regrid trace unchanged in that case.
            if(!regrid->find(key))return std::nullopt;
            const int value=integer(mesh,regrid,"Regrid",key,policy.minLevel,trace,"levels",consumer);
            if(value<policy.minLevel||value>policy.maxLevel)
                where(mesh,regrid,key).fail("a phase minimum in [minLevel,maxLevel]");
            return regrid->find(key)?std::optional<int>{value}:std::nullopt;
        };
        policy.liquidMinLevel=phaseMinimum("liquidMinLevel");policy.gasMinLevel=phaseMinimum("gasMinLevel");
        if(policy.interfaceBandLayers<0||policy.interfaceBandLayers>16)
            regrid->fail("interfaceBandLayers in [0,16] actual face-neighbor layers");
        if(!std::isfinite(policy.vorticityTolerance)||!(policy.vorticityTolerance>0))
            regrid->fail("finite positive vorticityTolerance in code velocity");
        if(!std::isfinite(policy.coarsenRatio)||!(policy.coarsenRatio>0&&policy.coarsenRatio<.5))
            regrid->fail("finite coarsenRatio in (0,.5)");
        if(static_cast<long long>(input.nx)*input.ny*(1LL<<(2*policy.maxLevel))>1000000)
            regrid->fail("base resolution times 4^maxLevel no greater than 1000000 native leaves");
        input.vofRegrid=policy;
    }
    const auto* reconstruction=interface.find("Reconstruction");
    if(reconstruction)reconstruction->allow({"relativeVolumeTolerance","maximumVolumeIterations","localResolution"});
    const auto localResolution=reconstruction&&reconstruction->find("localResolution")?
        trace.record(reconstruction->get("localResolution"),reconstruction->get("localResolution").choice({"volumeOnly","nativeCurvatureSupport"}),
            "dimensionless","Interface.GeometricVofPublication.localResolution"):
        trace.defaultValue(interface,"Reconstruction.localResolution",std::string("volumeOnly"),
            "dimensionless","Interface.GeometricVofPublication.localResolution");
    config.transport.geometry.localResolution=localResolution=="nativeCurvatureSupport"?
        LAMP::Interface::GeometricVofLocalResolution::NativeInterfaceSupport:LAMP::Interface::GeometricVofLocalResolution::VolumeOnly;
    if(localResolution=="nativeCurvatureSupport"&&(input.radial||input.csgPin||!input.vofBoundaryPatches.empty()))
        where(interface,reconstruction,"localResolution").fail("nativeCurvatureSupport requires closed unpinned XY local VOF");
    config.transport.geometry.relativeVolumeTolerance=number(interface,reconstruction,"Reconstruction","relativeVolumeTolerance",
        static_cast<double>(config.transport.geometry.relativeVolumeTolerance),trace,"dimensionless","Interface.GeometricVofPublication.localPlicVolume");
    if(!(config.transport.geometry.relativeVolumeTolerance>0&&config.transport.geometry.relativeVolumeTolerance<=1e-2))
        where(interface,reconstruction,"relativeVolumeTolerance").fail("local PLIC relative volume tolerance in (0,1e-2]");
    const int iterations=integer(interface,reconstruction,"Reconstruction","maximumVolumeIterations",
        static_cast<int>(config.transport.geometry.maximumVolumeIterations),trace,"iterations","Interface.GeometricVofPublication.localPlicVolume");
    if(iterations<1||iterations>128)where(interface,reconstruction,"maximumVolumeIterations").fail("local PLIC volume iteration budget in [1,128]");
    if(input.vofRegrid&&config.transport.geometry.relativeVolumeTolerance>1e-12L)
        where(interface,reconstruction,"relativeVolumeTolerance").fail("Mesh.Regrid requires reconstruction relativeVolumeTolerance <= 1e-12 for independent phase transfer");
    config.transport.geometry.maximumVolumeIterations=static_cast<unsigned>(iterations);
    const auto* initialization=interface.find("Initialization");
    if(initialization)initialization->allow({"subcellsPerAxis"});
    const int subcells=integer(interface,initialization,"Initialization","subcellsPerAxis",
        static_cast<int>(input.vofInitialization.subcellsPerAxis),trace,"subcells per axis",
        "Interface.GeometricVofInitialization.fixedSubcells");
    if(subcells<1||subcells>32)where(interface,initialization,"subcellsPerAxis").fail("initial subcellsPerAxis in [1,32]");
    input.vofInitialization.subcellsPerAxis=static_cast<unsigned>(subcells);
    input.vofInitialization.geometry=config.transport.geometry;
    const auto initialMesh=mesh.find("initialMesh")?trace.record(mesh.get("initialMesh"),mesh.get("initialMesh").choice({"uniform","interfaceBand"}),"dimensionless","Mesh.initialGrid"):
        trace.defaultValue(mesh,"initialMesh",std::string("uniform"),"dimensionless","Mesh.initialGrid");
    input.initialInterfaceBand=initialMesh=="interfaceBand";
    if(input.initialInterfaceBand){const auto& refinement=mesh.block("Refinement");refinement.allow({"level","bandCells","bandScale"});
        input.initialRefinementLevel=trace.record(refinement.get("level"),refinement.get("level").integer(),"levels","Hydrodynamics.VOF.initialInterfaceBand");
        const auto bandScale=refinement.find("bandScale")
            ?trace.record(refinement.get("bandScale"),refinement.get("bandScale").choice({"baseCell","currentCell","targetCell"}),"dimensionless","Hydrodynamics.VOF.initialInterfaceBand")
            :trace.defaultValue(refinement,"bandScale",std::string("baseCell"),"dimensionless","Hydrodynamics.VOF.initialInterfaceBand");
        input.initialBandScale=bandScale=="targetCell"?HydrodynamicsInput::InitialBandScale::TargetCell:
            bandScale=="currentCell"?HydrodynamicsInput::InitialBandScale::CurrentCell:HydrodynamicsInput::InitialBandScale::BaseCell;
        const char* bandUnits=bandScale=="targetCell"?"target cells":bandScale=="currentCell"?"current cells":"base cells";
        input.initialInterfaceBandCells=refinement.find("bandCells")?trace.record(refinement.get("bandCells"),refinement.get("bandCells").number(),bandUnits,"Hydrodynamics.VOF.initialInterfaceBand"):
            trace.defaultValue(refinement,"bandCells",2.,bandUnits,"Hydrodynamics.VOF.initialInterfaceBand");
        if(input.initialRefinementLevel<1||input.initialRefinementLevel>6)refinement.get("level").fail("initial interface-band refinement level in [1,6], within the native leaf budget");
        if(!(input.initialInterfaceBandCells>0))refinement.fail("finite positive interface-band width in its declared bandScale");
        if(static_cast<long long>(input.nx)*input.ny*(1LL<<(2*input.initialRefinementLevel))>1000000)
            refinement.get("level").fail("base resolution times 4^level no greater than 1000000 native leaves");
    }else if(mesh.find("Refinement"))mesh.get("Refinement").fail("omitted for initialMesh=uniform; its fields would be unused");
    const auto densityRepresentation=hydro.find("densityRepresentation")
        ?trace.record(hydro.get("densityRepresentation"),hydro.get("densityRepresentation").choice({"cellAverageProjection"}),"dimensionless","HydroVofStep.momentumDensityRepresentation")
        :trace.defaultValue(hydro,"densityRepresentation",std::string("cellAverageProjection"),"dimensionless","HydroVofStep.momentumDensityRepresentation");
    config.momentumDensityRepresentation=HydroMomentumMass::DensityRepresentation::CellAverageProjection;
    if(const auto* capillary=hydro.find("capillaryModel"))trace.record(*capillary,capillary->choice({"heightFunction"}),"dimensionless","HydroBalancedCsfLoad.localPlicHeightFunctionCapillary");
    else trace.defaultValue(hydro,"capillaryModel",std::string("heightFunction"),"dimensionless","HydroBalancedCsfLoad.localPlicHeightFunctionCapillary");
    config.curvature.forceMethod=HydroCsfForceMethod::HeightFunction;
    const auto& integration=hydro.block("TimeIntegration");
    integration.allow({"method","interfaceMismatchTolerance","maximumCorrectionRelative","correctionReferenceVelocity","capillaryCorrectionCourant","maximumCapillaryCourant","resetVelocityHistory"});
    const auto scheme=trace.record(integration.get("method"),integration.get("method").choice({"adaptiveLinear","extrapolatedMidpoint"}),"dimensionless","HydroVofStep.AC");
    config.scheme=scheme=="adaptiveLinear"?HydroTwoPhaseTimeScheme::AdaptiveLinearMidpoint:HydroTwoPhaseTimeScheme::ExtrapolatedMidpoint;
    const auto ac=[&](const char* key,double value,const char* units){const auto result=number(hydro,&integration,"TimeIntegration",key,value,trace,units,"HydroVofStep.AC");
        if(!(result>0))where(hydro,&integration,key).fail("finite positive VOF AC parameter");return result;};
    config.maximumCorrectionRelative=ac("maximumCorrectionRelative",config.maximumCorrectionRelative,"dimensionless");
    config.correctionReferenceVelocity=ac("correctionReferenceVelocity",config.correctionReferenceVelocity,"code_velocity");
    if(config.scheme==HydroTwoPhaseTimeScheme::AdaptiveLinearMidpoint){
        config.interfaceMismatchTolerance=ac("interfaceMismatchTolerance",config.interfaceMismatchTolerance,"dimensionless");
        config.capillaryCorrectionCourant=ac("capillaryCorrectionCourant",config.capillaryCorrectionCourant,"dimensionless");
        config.maximumCapillaryCourant=ac("maximumCapillaryCourant",config.maximumCapillaryCourant,"dimensionless");
        if(config.capillaryCorrectionCourant>config.maximumCapillaryCourant)integration.fail("capillaryCorrectionCourant <= maximumCapillaryCourant");
    }else for(const char* key:{"interfaceMismatchTolerance","capillaryCorrectionCourant","maximumCapillaryCourant"})
        if(const auto* unused=integration.find(key))unused->fail("omitted for extrapolatedMidpoint; adaptive-only trigger would be unused");
    config.resetVelocityHistory=integration.find("resetVelocityHistory")?trace.record(integration.get("resetVelocityHistory"),integration.get("resetVelocityHistory").boolean(),"dimensionless","HydroVofStep.velocityHistory"):
        trace.defaultValue(integration,"resetVelocityHistory",config.resetVelocityHistory,"dimensionless","HydroVofStep.velocityHistory");
    const auto* numerics=hydro.find("Numerics");if(numerics)numerics->allow({"maximumLocalCourant","relativeConservationTolerance"});
    config.transport.maximumLocalCourant=number(hydro,numerics,"Numerics","maximumLocalCourant",static_cast<double>(config.transport.maximumLocalCourant),trace,"dimensionless","Interface.GeometricVofTransport.localPlicSweep");
    config.transport.relativeConservationTolerance=number(hydro,numerics,"Numerics","relativeConservationTolerance",static_cast<double>(config.transport.relativeConservationTolerance),trace,"dimensionless","Interface.GeometricVofTransport.localPlicSweep");
    if(!(config.transport.maximumLocalCourant>0&&config.transport.maximumLocalCourant<=.5))where(hydro,numerics,"maximumLocalCourant").fail("local donor-column Courant in (0,.5]");
    if(!(config.transport.relativeConservationTolerance>0&&config.transport.relativeConservationTolerance<=1e-2))where(hydro,numerics,"relativeConservationTolerance").fail("relative conservation tolerance in (0,1e-2]");
    auto meshText=stream(),interfaceText=stream(),hydroText=stream();
    meshText<<"Mesh { resolution ("<<input.nx<<' '<<input.ny<<"); initialMesh "<<initialMesh<<';';
    if(input.initialInterfaceBand){meshText<<" Refinement { level "<<input.initialRefinementLevel<<"; bandCells "<<input.initialInterfaceBandCells<<';';
        if(input.initialBandScale==HydrodynamicsInput::InitialBandScale::CurrentCell)meshText<<" bandScale currentCell;";
        else if(input.initialBandScale==HydrodynamicsInput::InitialBandScale::TargetCell)meshText<<" bandScale targetCell;";
        meshText<<" }";}
    if(input.vofRegrid){const auto& policy=*input.vofRegrid;
        meshText<<" Regrid { minLevel "<<policy.minLevel<<"; maxLevel "<<policy.maxLevel
            <<"; interfaceBandLayers "<<policy.interfaceBandLayers<<"; vorticityTolerance "<<policy.vorticityTolerance
            <<"; coarsenRatio "<<policy.coarsenRatio<<';';
        if(policy.liquidMinLevel)meshText<<" liquidMinLevel "<<*policy.liquidMinLevel<<';';
        if(policy.gasMinLevel)meshText<<" gasMinLevel "<<*policy.gasMinLevel<<';';
        if(policy.interval!=1)meshText<<" interval "<<policy.interval<<';';
        meshText<<" }";
    }
    meshText<<" }\n";
    interfaceText<<"Interface { method vof; Boundary { type "<<boundaryType<<';';
    if(!symmetry.empty()){interfaceText<<" sides (";for(const auto& side:symmetry)interfaceText<<side<<' ';interfaceText<<");";}
    if(input.vofSouthEven90Interval){interfaceText<<" interval ";vectorText(interfaceText,*input.vofSouthEven90Interval);interfaceText<<';';}
    interfaceText<<" } Initialization { subcellsPerAxis "<<input.vofInitialization.subcellsPerAxis
        <<"; } Reconstruction { relativeVolumeTolerance "<<static_cast<double>(config.transport.geometry.relativeVolumeTolerance)
        <<"; maximumVolumeIterations "<<config.transport.geometry.maximumVolumeIterations
        <<"; localResolution "<<localResolution<<"; } }\n";
    hydroText<<"Hydrodynamics { method staggered; capillaryModel heightFunction; densityRepresentation "<<densityRepresentation<<"; Initial { velocity ";vectorText(hydroText,input.initialVelocity);
    hydroText<<"; } Boundary { type ";
    if(input.vofBoundaryPatches.empty())hydroText<<"impermeable; }";
    else{hydroText<<(input.vofPressureDrivenReservoirs?"pressureDrivenReservoirs":"prescribedFlow")<<"; Patches {";
        constexpr const char* sides[]{"W","S","E","N"};
        for(const auto& patch:input.vofBoundaryPatches){using Patch=HydrodynamicsInput::VofBoundaryPatch;
            hydroText<<' '<<std::quoted(patch.name)<<" { ";
            if(patch.boundaryId)hydroText<<"boundaryId "<<*patch.boundaryId;
            else{hydroText<<"side "<<sides[static_cast<int>(patch.side)]<<"; interval ";vectorText(hydroText,patch.interval);}
            hydroText<<"; type "<<(patch.kind==Patch::Kind::PrescribedVelocity?"prescribedVelocity":"pressureOutlet")<<"; phase "
                <<(patch.reservoir==LAMP::Interface::GeometricVofAdvectorSlab::Boundary::Mode::PureLiquidReservoir?"liquid":"gas")<<';';
            if(patch.kind==Patch::Kind::PrescribedVelocity){hydroText<<" velocity ";vectorText(hydroText,patch.velocity);}
            else hydroText<<" pressure "<<patch.pressure;
            hydroText<<"; }";}
        hydroText<<" }";
        if(input.vofPressureOutletBackflow){const auto& backflow=*input.vofPressureOutletBackflow;
            const bool ambient=backflow.model==HydroVofPressureOutletMomentumModel::Ambient;
            hydroText<<" Backflow { model "<<(ambient?"ambient":"interiorContinuationZeroShear")<<';';
            if(ambient){hydroText<<" velocity ";vectorText(hydroText,*backflow.velocity);hydroText<<';';}
            hydroText<<" }";
        }
        if(!input.vofStationaryNoSlipBoundaryIds.empty()){
            hydroText<<" StationaryNoSlip { boundaryIds (";
            for(const int id:input.vofStationaryNoSlipBoundaryIds)hydroText<<id<<' ';
            hydroText<<"); }";
        }
        hydroText<<" }";}
    hydroText<<" BodyForce { acceleration ";vectorText(hydroText,config.bodyAcceleration);hydroText<<"; } TimeIntegration { method "<<scheme
        <<"; maximumCorrectionRelative "<<config.maximumCorrectionRelative<<"; correctionReferenceVelocity "<<config.correctionReferenceVelocity<<';';
    if(config.scheme==HydroTwoPhaseTimeScheme::AdaptiveLinearMidpoint)hydroText<<" interfaceMismatchTolerance "<<config.interfaceMismatchTolerance<<"; capillaryCorrectionCourant "<<config.capillaryCorrectionCourant<<"; maximumCapillaryCourant "<<config.maximumCapillaryCourant<<';';
    hydroText<<" resetVelocityHistory "<<(config.resetVelocityHistory?"true":"false")<<"; } Numerics { maximumLocalCourant "<<static_cast<double>(config.transport.maximumLocalCourant)
        <<"; relativeConservationTolerance "<<static_cast<double>(config.transport.relativeConservationTolerance)<<"; } }\n";
    return {meshText.str(),interfaceText.str(),hydroText.str()};
}
void recordVofInputDefaults(const Node& hydro,const HydrodynamicsInput& input,Trace& trace)
{
    const auto& c=input.vof;const HydroSaddleConfig saddle;
    trace.defaultValue(hydro,"KernelDefaults.linear.method",std::string("EigenSparseLU"),"dimensionless","HydroSinglePhaseSaddle.quotientFactor");
    trace.defaultValue(hydro,"KernelDefaults.linear.algebraicRelativeTolerance",saddle.algebraicRelativeTolerance,"dimensionless","HydroSinglePhaseSaddle.solve");
    trace.defaultValue(hydro,"KernelDefaults.linear.momentumRelativeTolerance",saddle.momentumRelativeTolerance,"dimensionless","HydroSinglePhaseSaddle.solve");
    trace.defaultValue(hydro,"KernelDefaults.linear.continuityRelativeTolerance",saddle.continuityRelativeTolerance,"dimensionless","HydroSinglePhaseSaddle.solve");
    trace.defaultValue(hydro,"KernelDefaults.linear.pressureMeanRelativeTolerance",saddle.pressureMeanRelativeTolerance,"dimensionless","HydroSinglePhaseSaddle.solve");
    trace.defaultValue(hydro,"KernelDefaults.linear.boundaryFluxRelativeTolerance",saddle.boundaryFluxRelativeTolerance,"dimensionless","HydroSinglePhaseSaddle.solve");
    trace.defaultValue(hydro,"KernelDefaults.momentumPowerTolerance",c.momentumPowerTolerance,"dimensionless","HydroVofStep.acceptance");
    trace.defaultValue(hydro,"KernelDefaults.material.method",std::string("localPlicSplit"),"dimensionless","Interface.GeometricVofTransport.localPlicSweep");
    trace.defaultValue(hydro,"KernelDefaults.viscosity.material",std::string(c.liquidViscosity==0&&c.gasViscosity==0?
        "exactZero":"localPlicInventoryCellAverage"),"dimensionless","HydroPhaseViscosity");
    trace.defaultValue(hydro,"KernelDefaults.interface.maximumBracketExpansions",static_cast<int>(c.transport.geometry.maximumBracketExpansions),"iterations","Interface.GeometricVofPublication.localPlicVolume");
}
}
