#include "Electrohydrodynamics/IO/VofElectrohydrodynamicsInput.h"
#include "IO/TwoPhaseInputComponents.h"
#include "Core/IO/DictionaryGeometry.h"
#include <set>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace
{
using Node=DictionaryInput::Node;
using Trace=DictionaryInput::ConfigurationTrace;
constexpr std::array<const char*,4> sides{{"W","S","E","N"}};
std::ostringstream stream(){std::ostringstream out;out<<std::setprecision(17);return out;}
double solverTolerance(const Node& electric,const Node* solver,const char* key,double value,double maximum,
    Trace& trace,const char* consumer)
{
    const auto* node=solver?solver->find(key):nullptr;
    if(node){value=node->number();if(!(value>0&&value<=maximum))node->fail("a positive tolerance no greater than the original Electric bound");
        return trace.record(*node,value,"dimensionless",consumer);}
    return solver?trace.defaultValue(*solver,key,value,"dimensionless",consumer)
        :trace.defaultValue(electric,std::string("Solver.")+key,value,"dimensionless",consumer);
}
}
double VofElectrohydrodynamicsInput::initialChargeDensity(Vector2D point) const
{
    if(method!=Method::VofFiniteVolume||initialChargeType!=InitialChargeType::CellDensity||!chargeExpression_)
        throw std::logic_error("VOF-EHD initial volume charge has not been validated for this method and Initial.type");
    return chargeExpression_->sample(point.x,point.y).value;
}
double VofElectrohydrodynamicsInput::initialSurfaceChargeDensity(Vector2D point) const
{
    if((method!=Method::VofSharpSurface&&!(method==Method::VofFiniteVolume&&initialChargeType==InitialChargeType::SurfaceChargeProjection))||!chargeExpression_)
        throw std::logic_error("VOF-EHD initial surface charge has not been validated for this method and Initial.type");
    return chargeExpression_->sample(point.x,point.y).value;
}
const char* vofElectrohydrodynamicsApplicationName(VofElectrohydrodynamicsInput::Application application)
{
    switch(application){
    case VofElectrohydrodynamicsInput::Application::General:return "VofElectrohydrodynamics";
    case VofElectrohydrodynamicsInput::Application::Electrospray:return "Electrospray";
    }
    throw std::invalid_argument("unknown VOF-EHD production application profile");
}
VofElectrohydrodynamicsInput readVofElectrohydrodynamicsInput(const std::string& systemDir,std::optional<int> threadsOverride,
    VofElectrohydrodynamicsInput::Application application)
{
    const auto* applicationName=vofElectrohydrodynamicsApplicationName(application);
    const AppIO::CaseDocuments documents(systemDir,
        {"run.dict","geometry.dict","mesh.dict","interface.dict","materials.dict","hydro.dict","electric.dict"},
        applicationName,"VOF-EHD seven-document input");
    const auto& method=documents.root("interface.dict").block("Interface").get("method");
    if(method.word()!="vof")method.fail("vof for VofElectrohydrodynamics production");
    for(const auto& phase:documents.root("materials.dict").block("Materials").block("Phases").entries)
        phase.allow({"phaseId","density","dynamicViscosity","permittivity","conductivity"});
    VofElectrohydrodynamicsInput input;input.application=application;auto& fluid=input.fluid;auto& trace=fluid.trace;AppIO::ConstantsUsage usage;
    documents.root("electric.dict").allow({"Electric"});const auto& electric=documents.root("electric.dict").block("Electric");
    const auto selected=electric.get("method").choice({"vofFiniteVolume","vofSharpSurface"});
    const bool sharp=selected=="vofSharpSurface";
    if(application==VofElectrohydrodynamicsInput::Application::Electrospray&&sharp)
        electric.get("method").fail("vofFiniteVolume for the implemented Electrospray production profile");
    const auto components=AppIO::readTwoPhaseInputComponents(fluid,documents,usage,threadsOverride,sharp?0:10);
    if(application==VofElectrohydrodynamicsInput::Application::Electrospray){
        const auto& geometryInitial=documents.root("geometry.dict").block("Geometry").block("Initial");
        if(!fluid.csgGeometry||!fluid.csgSampler)
            geometryInitial.fail("Electrospray requires its declared immutable Core CSG source and excluded fixed solids");
        if(!fluid.csgPin)geometryInitial.fail("Electrospray requires an explicit stock-derived CSG Pin declaration");
        const auto& hydroBoundary=documents.root("hydro.dict").block("Hydrodynamics").block("Boundary");
        if(!fluid.vofPressureDrivenReservoirs)
            hydroBoundary.get("type").fail("pressureDrivenReservoirs for Electrospray without prescribed feed velocity");
        if(fluid.vofStationaryNoSlipBoundaryIds.empty())
            hydroBoundary.fail("Electrospray requires explicit StationaryNoSlip fixed-solid boundary IDs");
    }
    trace.record(electric.get("method"),selected,"dimensionless",sharp?"ElectricVofSharpField.dispatch":"ElectricVofField.dispatch");
    input.method=sharp?VofElectrohydrodynamicsInput::Method::VofSharpSurface:VofElectrohydrodynamicsInput::Method::VofFiniteVolume;
    if(sharp){
        const auto& run=documents.root("run.dict");
        if(!fluid.run.restartDirectory.empty())run.block("Restart").fail("B1 restart is not implemented; only Hydro/A VOF checkpoints are supported");
        if(fluid.run.checkpointExplicit&&fluid.run.checkpointEverySteps!=0)
            run.block("Checkpoint").get("everySteps").fail("0 for B1; surface-charge checkpoints are not implemented");
    }
    if(sharp){electric.allow({"method","Initial","Boundary","Solver","SmoothGeometry","SurfaceTransport","Coupling"});
        if(fluid.radial)electric.get("method").fail("XY for the current B1 material surface-charge transport");
        if(fluid.vofRegrid||fluid.initialInterfaceBand)electric.get("method").fail("B1 fixed uniform mesh without initial refinement or Mesh.Regrid until actual surface-charge remap is implemented");
        if(!fluid.vofBoundaryPatches.empty()||!fluid.vofSymmetryPlanes.empty()||fluid.vofSouthEven90Interval||fluid.vof.southEven90Boundary)
            electric.get("method").fail("B1 closed interfaces without prescribedFlow patches, symmetry contact or southEven90 boundary capability");
    }else electric.allow({"method","Initial","Boundary","Solver","Advection","Coupling","Regrid"});
    if(const auto* regrid=electric.find("Regrid")){
        regrid->allow({"chargeDensityJumpTolerance"});
        if(!fluid.vofRegrid)regrid->fail("Electric.Regrid only with an actual Mesh.Regrid policy");
        input.chargeDensityJumpTolerance=AppIO::boundedCaseValue(regrid->get("chargeDensityJumpTolerance"),trace,0,true,
            "ElectricVofRefinement.signedDensityJump","code_charge/code_volume");
    }
    const char* materialConsumer=sharp?"ElectricVofSharpField.material":"ElectricVofField.material";
    const auto& phases=documents.root("materials.dict").block("Materials").block("Phases").entries;
    for(std::size_t i=0;i<phases.size();++i){const auto& phase=phases[i];const bool liquid=components.materials[i].phaseId==fluid.negativePhase;
        const auto epsilon=AppIO::boundedCaseValue(phase.get("permittivity"),trace,0,true,materialConsumer,"code_permittivity");
        const auto conductivity=AppIO::boundedCaseValue(phase.get("conductivity"),trace,0,false,materialConsumer,"code_conductivity");
        if(liquid){input.materials.epsilonLiquid=epsilon;input.materials.conductivityLiquid=conductivity;}
        else{input.materials.epsilonGas=epsilon;input.materials.conductivityGas=conductivity;}}
    const auto& initial=electric.block("Initial");const auto* initialType=sharp?nullptr:initial.find("type");
    std::string initialTypeName="cellDensity";
    if(!sharp){
        initialTypeName=initialType?trace.record(*initialType,initialType->choice({"cellDensity","surfaceChargeProjection"}),
            "dimensionless","ElectricVofCharge.initialization"):trace.defaultValue(initial,"type",initialTypeName,
            "dimensionless","ElectricVofCharge.initialization");
        input.initialChargeType=initialTypeName=="surfaceChargeProjection"?VofElectrohydrodynamicsInput::InitialChargeType::SurfaceChargeProjection
            :VofElectrohydrodynamicsInput::InitialChargeType::CellDensity;
    }
    const bool surfaceProjection=!sharp&&input.initialChargeType==VofElectrohydrodynamicsInput::InitialChargeType::SurfaceChargeProjection;
    if(sharp)initial.allow({"surfaceChargeDensity"});
    else if(surfaceProjection)initial.allow({"type","surfaceChargeDensity"});
    else initial.allow({"type","chargeDensity"});
    if(surfaceProjection&&(fluid.radial||!fluid.vofBoundaryPatches.empty()||!fluid.vofSymmetryPlanes.empty()
        ||fluid.vofSouthEven90Interval||fluid.vof.southEven90Boundary))
        initial.get("type").fail("surfaceChargeProjection requires closed XY Initial geometry without contact boundary data");
    const char* densityKey=sharp||surfaceProjection?"surfaceChargeDensity":"chargeDensity";
    const auto& charge=initial.get(densityKey);
    const auto formula=trace.record(charge,charge.word(),sharp||surfaceProjection?"code_charge/code_interface_measure":"code_charge/code_volume",
        sharp?"ElectricVofSurfaceCharge.initialMoments":surfaceProjection?"ElectricVofCharge.initialSurfaceMomentProjection":"ElectricVofCharge.initialIntegratedCellCharge");
    usage.include(formula);
    try{
        input.chargeExpression_=std::make_shared<const LAMP::ScalarExpression>(formula,fluid.constants,
            fluid.radial?LAMP::ScalarExpression::Coordinates::Meridional:LAMP::ScalarExpression::Coordinates::Planar);
        const auto sample=[&](Vector2D p){return sharp||surfaceProjection?input.initialSurfaceChargeDensity(p):input.initialChargeDensity(p);};
        for(const auto x:{fluid.domainMin.x,fluid.domainMax.x})for(const auto y:{fluid.domainMin.y,fluid.domainMax.y})sample({x,y});
        sample({fluid.domainMin.x+(fluid.domainMax.x-fluid.domainMin.x)/2,fluid.domainMin.y+(fluid.domainMax.y-fluid.domainMin.y)/2});
    }catch(const std::exception& error){charge.fail(std::string("a finite signed spatial charge expression without t: ")+error.what());}
    const auto& boundary=electric.block("Boundary");bool anchored=false;
    std::array<std::string,4> types;
    const bool csgBoundary=boundary.find("Patches")!=nullptr;
    if(fluid.csgGeometry&&!csgBoundary)
        boundary.fail("explicit Patches with Core boundaryId for the selected CSG fluid domain");
    if(csgBoundary){
        boundary.allow({"Patches"});
        if(sharp||!fluid.csgGeometry||!fluid.csgSampler)
            boundary.get("Patches").fail("A-route CSG geometry with its validated immutable sampler");
        const auto& patches=boundary.block("Patches");
        if(patches.entries.empty())patches.fail("nonempty named Core boundary-ID conditions");
        std::set<int> ids;
        for(const auto& entry:patches.entries){
            if(!entry.isBlock)entry.fail("a named electric boundary patch block");
            entry.allow({"boundaryId","type","value"});VofElectrohydrodynamicsInput::CsgElectricBoundary patch;
            patch.name=entry.name;patch.boundaryId=trace.record(entry.get("boundaryId"),entry.get("boundaryId").integer(),
                "1-based Core boundary ID","ElectricVofField.fittedCsgBoundary");
            const auto found=fluid.csgGeometry->boundaryById.find(patch.boundaryId);
            if(patch.boundaryId<1||found==fluid.csgGeometry->boundaryById.end()||found->second.ignored||!ids.insert(patch.boundaryId).second)
                entry.get("boundaryId").fail("one distinct declared nonignored Core CSG boundary ID");
            const auto type=trace.record(entry.get("type"),entry.get("type").choice({"Dirichlet","Neumann"}),
                "dimensionless","ElectricVofField.fittedCsgBoundary");
            patch.condition.type=type=="Dirichlet"?LAMP::Electric::ElectricVofBoundary::Type::Dirichlet:
                LAMP::Electric::ElectricVofBoundary::Type::HomogeneousNeumann;
            patch.condition.value=trace.record(entry.get("value"),entry.get("value").number(),
                type=="Dirichlet"?"code_potential":"code_displacement_flux","ElectricVofField.fittedCsgBoundary");
            if(type=="Neumann"&&patch.condition.value!=0)entry.get("value").fail("zero for homogeneous Neumann");
            anchored=anchored||type=="Dirichlet";input.csgElectricBoundaries.push_back(std::move(patch));
        }
    }else{
    boundary.allow({"W","S","E","N"});
    for(std::size_t side=0;side<sides.size();++side){const auto& data=boundary.block(sides[side]);data.allow({"type","value"});
        types[side]=trace.record(data.get("type"),data.get("type").choice({"Dirichlet","Neumann"}),"dimensionless","ElectricVofField.staticBoundary");
        const auto value=trace.record(data.get("value"),data.get("value").number(),
            types[side]=="Dirichlet"?"code_potential":"code_displacement_flux","ElectricVofField.staticBoundary");
        auto& result=input.boundaries[side];result.type=types[side]=="Dirichlet"
            ?LAMP::Electric::ElectricVofBoundary::Type::Dirichlet:LAMP::Electric::ElectricVofBoundary::Type::HomogeneousNeumann;
        result.value=value;
        if(result.type==LAMP::Electric::ElectricVofBoundary::Type::HomogeneousNeumann&&value!=0)
            data.get("value").fail("zero for the implemented homogeneous Neumann boundary");
        if(fluid.radial&&fluid.domainMin.x==0&&side==0&&result.type!=LAMP::Electric::ElectricVofBoundary::Type::HomogeneousNeumann)
            data.get("type").fail("homogeneous Neumann with zero value on the true RZ axis");
        anchored=anchored||result.type==LAMP::Electric::ElectricVofBoundary::Type::Dirichlet;
    }

    }
    if(!anchored)boundary.fail(csgBoundary?"at least one Dirichlet boundary ID to anchor the electric potential":
        "at least one Dirichlet side to anchor the electric potential");
    const auto* advection=electric.find("Advection");const Node* inlets=nullptr;
    std::string advectionScheme="donorCell";
    if(advection){
        advection->allow({"scheme","Inlets"});
        if(advection->find("Inlets")){
            if(fluid.vofBoundaryPatches.empty())advection->get("Inlets").fail("Electric.Advection.Inlets only with actual prescribedFlow Hydro patches");
            inlets=&advection->block("Inlets");}}
    if(!sharp){
        const auto* method=advection?advection->find("scheme"):nullptr;
        advectionScheme=method?trace.record(*method,method->choice({"donorCell","limitedLinearFct"}),"dimensionless","ElectricVofChargeTransport.sharedFluxAdvection")
            :advection?trace.defaultValue(*advection,"scheme",std::string("limitedLinearFct"),"dimensionless","ElectricVofChargeTransport.sharedFluxAdvection")
            :trace.defaultValue(electric,"Advection.scheme",std::string("limitedLinearFct"),"dimensionless","ElectricVofChargeTransport.sharedFluxAdvection");
        input.controls.chargeAdvection=advectionScheme=="limitedLinearFct"?LAMP::Electric::ElectricVofChargeAdvectionScheme::LimitedLinearFct
            :LAMP::Electric::ElectricVofChargeAdvectionScheme::DonorCell;
    }
    input.inletChargeDensityByPatch.resize(fluid.vofBoundaryPatches.size());
    const auto incoming=[&fluid](const HydrodynamicsInput::VofBoundaryPatch& patch){using Patch=HydrodynamicsInput::VofBoundaryPatch;
        // This physical charge reservoir is independent of the weak momentum
        // model/vector and is consumed only on actual incoming gross volume.
        if(patch.kind==Patch::Kind::PressureOutlet)return fluid.vofPressureOutletBackflow.has_value();
        if(patch.kind!=Patch::Kind::PrescribedVelocity)return false;
        const bool horizontal=patch.side==Patch::Side::S||patch.side==Patch::Side::N;
        const double normal=horizontal?patch.velocity.y:patch.velocity.x;
        return (patch.side==Patch::Side::W||patch.side==Patch::Side::S)?normal>0:normal<0;};
    if(inlets)for(const auto& node:inlets->entries){
        node.allow({"chargeDensity"});const auto it=std::find_if(fluid.vofBoundaryPatches.begin(),fluid.vofBoundaryPatches.end(),
            [&](const auto& patch){return patch.name==node.name;});
        if(it==fluid.vofBoundaryPatches.end()||!incoming(*it))node.fail("a real incoming prescribedVelocity or explicitly opted-in pressureOutlet Hydro patch name");}
    for(std::size_t k=0;k<fluid.vofBoundaryPatches.size();++k){const auto& patch=fluid.vofBoundaryPatches[k];if(!incoming(patch))continue;
        const auto* entry=inlets?inlets->find(patch.name):nullptr;const auto* density=entry?entry->find("chargeDensity"):nullptr;
        input.inletChargeDensityByPatch[k]=density?trace.record(*density,density->number(),"code_charge/code_volume","ElectricVofChargeBoundary.signedInflowDensity"):
            entry?trace.defaultValue(*entry,"chargeDensity",0.,"code_charge/code_volume","ElectricVofChargeBoundary.signedInflowDensity"):
            inlets?trace.defaultValue(*inlets,patch.name+".chargeDensity",0.,"code_charge/code_volume","ElectricVofChargeBoundary.signedInflowDensity"):
            trace.defaultValue(electric,"Advection.Inlets."+patch.name+".chargeDensity",0.,"code_charge/code_volume","ElectricVofChargeBoundary.signedInflowDensity");}
    const auto* solver=electric.find("Solver");
    if(solver){if(sharp)solver->allow({"relativeLinearTolerance","relativeConservationTolerance"});
        else solver->allow({"backend","relativeLinearTolerance","relativeConservationTolerance"});}
    std::string backend="gmres";
    if(!sharp){const auto* selectedBackend=solver?solver->find("backend"):nullptr;
        backend=selectedBackend?trace.record(*selectedBackend,selectedBackend->choice({"gmres","superlu"}),"dimensionless","ElectricVofField.strictIntegratedBackend")
            :solver?trace.defaultValue(*solver,"backend",backend,"dimensionless","ElectricVofField.strictIntegratedBackend")
            :trace.defaultValue(electric,"Solver.backend",backend,"dimensionless","ElectricVofField.strictIntegratedBackend");
        input.controls.backend=backend=="superlu"?LAMP::Electric::ElectricVofControls::Backend::SuperLU:LAMP::Electric::ElectricVofControls::Backend::Gmres;}
    input.controls.relativeLinearTolerance=solverTolerance(electric,solver,"relativeLinearTolerance",input.controls.relativeLinearTolerance,1e-11,
        trace,"ElectricVofField.strictLinearResidual");
    input.controls.relativeConservationTolerance=solverTolerance(electric,solver,"relativeConservationTolerance",input.controls.relativeConservationTolerance,1e-10,
        trace,"ElectricVofField.physicalGaussAndChargeBalance");
    const auto readBudget=[&](const char* block,const char* key,long double value,double maximum){
        const auto* parent=electric.find(block);const auto* node=parent?parent->find(key):nullptr;
        if(node){const auto x=node->number();if(!(x>0&&x<=maximum))node->fail("a positive value no larger than the original declared budget");
            return trace.record(*node,x,"dimensionless",block);}
        return parent?trace.defaultValue(*parent,key,static_cast<double>(value),"dimensionless",block)
            :trace.defaultValue(electric,std::string(block)+"."+key,static_cast<double>(value),"dimensionless",block);};
    if(sharp){
        if(const auto* b=electric.find("SmoothGeometry"))b->allow({"relativeVolumeTolerance"});
        if(const auto* b=electric.find("SurfaceTransport"))b->allow({"maximumDisplacementInCells","maximumProjectionDistanceInCells","relativeChargeTolerance"});
        input.smoothGeometry.relativeVolumeTolerance=readBudget("SmoothGeometry","relativeVolumeTolerance",input.smoothGeometry.relativeVolumeTolerance,1e-12);
        input.surfaceTransport.maximumDisplacementInCells=readBudget("SurfaceTransport","maximumDisplacementInCells",input.surfaceTransport.maximumDisplacementInCells,.45);
        input.surfaceTransport.maximumProjectionDistanceInCells=readBudget("SurfaceTransport","maximumProjectionDistanceInCells",input.surfaceTransport.maximumProjectionDistanceInCells,.2);
        input.surfaceTransport.relativeChargeTolerance=readBudget("SurfaceTransport","relativeChargeTolerance",input.surfaceTransport.relativeChargeTolerance,1e-10);
    }
    if(sharp||electric.find("Coupling")){
        if(const auto* b=electric.find("Coupling"))b->allow({"maximumForceChange","maximumWorkChange"});
        fluid.vof.requireMaxwellCorrection=true;
        fluid.vof.maximumMaxwellLoadChange=static_cast<double>(readBudget("Coupling","maximumForceChange",fluid.vof.maximumMaxwellLoadChange,.25));
        fluid.vof.maximumMaxwellWorkChange=static_cast<double>(readBudget("Coupling","maximumWorkChange",fluid.vof.maximumMaxwellWorkChange,.25));
    }
    usage.validate(documents.root("geometry.dict").block("Geometry"),fluid.constants,
        "a constant used by Geometry.Initial.phi or the selected Electric.Initial charge expression");
    AppIO::recordTwoPhaseInputDefaults(documents,fluid,trace);
    auto materialText=stream();materialText<<"Materials {\nPhases {\n";
    for(const auto& material:components.materials){const bool liquid=material.phaseId==fluid.negativePhase;
        materialText<<material.prefix<<" permittivity "<<(liquid?input.materials.epsilonLiquid:input.materials.epsilonGas)
            <<"; conductivity "<<(liquid?input.materials.conductivityLiquid:input.materials.conductivityGas)<<"; }\n";}
    materialText<<components.materialInterfaces;
    auto electricText=stream();electricText<<"Electric { method "<<selected<<"; Initial { ";
    if(!sharp&&initialType)electricText<<"type "<<initialTypeName<<"; ";
    electricText<<densityKey<<' '<<std::quoted(formula)<<"; }\nBoundary {\n";
    if(csgBoundary){electricText<<"Patches {\n";
        for(const auto& patch:input.csgElectricBoundaries)electricText<<std::quoted(patch.name)<<" { boundaryId "<<patch.boundaryId
            <<"; type "<<(patch.condition.type==LAMP::Electric::ElectricVofBoundary::Type::Dirichlet?"Dirichlet":"Neumann")
            <<"; value "<<patch.condition.value<<"; }\n";
        electricText<<"}\n";
    }else for(std::size_t side=0;side<sides.size();++side)electricText<<sides[side]<<" { type "<<types[side]<<"; value "<<input.boundaries[side].value<<"; }\n";
    electricText<<"}\n";
    if(!sharp){electricText<<"Advection { scheme "<<advectionScheme<<';';
        if(!fluid.vofBoundaryPatches.empty()){electricText<<" Inlets {";
            for(std::size_t k=0;k<fluid.vofBoundaryPatches.size();++k)if(incoming(fluid.vofBoundaryPatches[k]))
                electricText<<' '<<std::quoted(fluid.vofBoundaryPatches[k].name)<<" { chargeDensity "<<input.inletChargeDensityByPatch[k]<<"; }";
            electricText<<" }";}
        electricText<<" }\n";}
    if(input.chargeDensityJumpTolerance)electricText<<"Regrid { chargeDensityJumpTolerance "<<*input.chargeDensityJumpTolerance<<"; }\n";
    electricText<<"Solver { ";if(!sharp)electricText<<"backend "<<backend<<"; ";
    electricText<<"relativeLinearTolerance "<<input.controls.relativeLinearTolerance
        <<"; relativeConservationTolerance "<<input.controls.relativeConservationTolerance<<"; }\n";
    if(sharp)electricText<<"SmoothGeometry { relativeVolumeTolerance "<<input.smoothGeometry.relativeVolumeTolerance<<"; }\nSurfaceTransport { maximumDisplacementInCells "
        <<input.surfaceTransport.maximumDisplacementInCells<<"; maximumProjectionDistanceInCells "<<input.surfaceTransport.maximumProjectionDistanceInCells
        <<"; relativeChargeTolerance "<<input.surfaceTransport.relativeChargeTolerance<<"; }\n";
    if(sharp||fluid.vof.requireMaxwellCorrection)electricText<<"Coupling { maximumForceChange "<<fluid.vof.maximumMaxwellLoadChange
        <<"; maximumWorkChange "<<fluid.vof.maximumMaxwellWorkChange<<"; }\n";
    electricText<<"}\n";
    const std::array<const char*,6> files{{"run.dict","geometry.dict","mesh.dict","interface.dict","materials.dict","hydro.dict"}};
    for(std::size_t i=0;i<files.size();++i)fluid.documents.push_back({documents.path(files[i]),documents.sourceText(files[i]),i==4?materialText.str():components.effective[i]});
    fluid.documents.push_back({documents.path("electric.dict"),documents.sourceText("electric.dict"),electricText.str()});
    return input;
}
