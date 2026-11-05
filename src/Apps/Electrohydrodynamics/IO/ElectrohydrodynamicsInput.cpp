#include "Electrohydrodynamics/IO/ElectrohydrodynamicsInput.h"
#include "IO/TwoPhaseInputComponents.h"
#include <array>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace
{
using Node=DictionaryInput::Node;
using Trace=DictionaryInput::ConfigurationTrace;
using Expression=LAMP::ScalarExpression;
using AppIO::boundedCaseValue;
constexpr std::array<const char*,4> electricSides{"W","S","E","N"};
std::ostringstream stream(){std::ostringstream out;out<<std::setprecision(17);return out;}
void electricDefaults(const Node& electric,Trace& trace,const LAMP::Coupling::ElectrohydrodynamicStepConfig& config)
{
    const auto record=[&](const char* key,double value,const char* consumer){
        trace.defaultValue(electric,key,value,"dimensionless",consumer);
    };
    record("KernelDefaults.linear.relativeLinearTolerance",config.electric.relativeLinearTolerance,"Electric.electrostatic.linear");
    record("KernelDefaults.linear.relativeSymmetryTolerance",config.electric.relativeSymmetryTolerance,"Electric.electrostatic.symmetry");
    record("KernelDefaults.linear.relativeConservationTolerance",config.electric.relativeConservationTolerance,"Electric.electrostatic.Gauss");
    record("KernelDefaults.charge.relativeSchurTolerance",config.charge.relativeSchurTolerance,"Electric.Q.Schur");
    record("KernelDefaults.charge.relativeChargeBalanceTolerance",config.charge.relativeChargeBalanceTolerance,"Electric.Q.chargeBalance");
    record("KernelDefaults.charge.maximumStepRatio",config.charge.maximumStepRatio,"Electric.Q.variableBdf2");
    trace.defaultValue(electric,"KernelDefaults.charge.maximumIterations",config.charge.maximumIterations,"iterations","Electric.Q.Schur");
    trace.defaultValue(electric,"KernelDefaults.charge.restartLength",config.charge.restartLength,"iterations","Electric.Q.Schur");
    record("KernelDefaults.geometry.relativeMatrixTolerance",config.geometry.relativeMatrixTolerance,"Electric.electrostatic.phaseQuadrature");
    record("KernelDefaults.geometry.minimumRelativePhaseArea",config.geometry.minimumRelativePhaseArea,"Electric.electrostatic.smallCutRejection");
    trace.defaultValue(electric,"KernelDefaults.geometry.maximumRefinements",config.geometry.maximumRefinements,"levels","Electric.electrostatic.phaseQuadrature");
    trace.defaultValue(electric,"KernelDefaults.geometry.maximumPoints",std::to_string(config.geometry.maximumPoints),"points","Electric.electrostatic.phaseQuadrature");
    trace.defaultValue(electric,"KernelDefaults.geometry.maximumSliceQueriesPerCell",std::to_string(config.geometry.maximumSliceQueriesPerCell),"queries","Electric.electrostatic.phaseQuadrature");
    trace.defaultValue(electric,"KernelDefaults.geometry.maximumDepth",config.geometry.maximumDepth,"levels","Electric.electrostatic.phaseQuadrature");
    record("KernelDefaults.moments.relativeQuadratureTolerance",config.moments.relativeQuadratureTolerance,"Electric.P_Q");
    trace.defaultValue(electric,"KernelDefaults.moments.maximumAdditionalRefinements",config.moments.maximumAdditionalRefinements,"levels","Electric.P_Q");
    const LAMP::Electric::ElectricMaxwellControls maxwell;
    record("KernelDefaults.maxwell.relativeProjectionTolerance",maxwell.relativeProjectionTolerance,"Electric.Maxwell.projection");
    record("KernelDefaults.maxwell.relativeChargeJumpTolerance",maxwell.relativeChargeJumpTolerance,"Electric.Maxwell.jump");
    trace.defaultValue(electric,"KernelDefaults.bulkSource",std::string("zero"),"code_charge/code_volume","Electric.electrostatic.quasineutralBulk");
}
void validateCapabilities(const AppIO::CaseDocuments& documents)
{
    documents.root("interface.dict").block("Interface").get("method").choice({"hermiteGals"});
    const auto& geometry=documents.root("geometry.dict").block("Geometry");
    if(geometry.get("coordinates").choice({"xy","rz"})!="rz")
        geometry.get("coordinates").fail("rz for the implemented native moving sharp-interface electrostatics Hydro coupling");
    const auto& hydro=documents.root("hydro.dict").block("Hydrodynamics");
    if(const auto* model=hydro.find("capillaryModel"))
        if(model->choice({"balancedCsf","sharpNormalStress"})!="balancedCsf")
            model->fail("balancedCsf for the implemented sharp-interface electrostatics external traction; sharpNormalStress coupling is unavailable");
    if(const auto* csf=hydro.find("CSF"))if(const auto* method=csf->find("forceMethod")){
        const auto force=method->choice({"motionAdjoint","heightFunction","heightFunctionIncrement"});
        if(force!="motionAdjoint"){
            const auto* gradient=documents.root("interface.dict").block("Interface").find("gradientUpdate");
            if(force=="heightFunctionIncrement"&&(!gradient||gradient->choice({"characteristicPullback","compatibleIncrement"})!="compatibleIncrement"))
                method->fail("heightFunctionIncrement requires Interface.gradientUpdate=compatibleIncrement");
            method->fail("heightFunction for standalone Hydrodynamics; electric traction is not implemented");
        }
    }
}
}
double ElectrohydrodynamicsInput::initialSurfaceChargeDensity(Vector2D p) const
{
    if(!chargeExpression_)throw std::logic_error("coupled Hydro initial surface charge has not been validated");
    return chargeExpression_->sample(p.x,p.y).value;
}
ElectrohydrodynamicsInput readElectrohydrodynamicsInput(const std::string& systemDir,std::optional<int> threadsOverride)
{
    const AppIO::CaseDocuments documents(systemDir,
        {"run.dict","geometry.dict","mesh.dict","fields.dict","interface.dict","materials.dict","hydro.dict","electric.dict"},
        "Electrohydrodynamics");
    for(const auto& phase:documents.root("materials.dict").block("Materials").block("Phases").entries)
        phase.allow({"phaseId","density","dynamicViscosity","permittivity","conductivity"});
    validateCapabilities(documents);
    ElectrohydrodynamicsInput input;auto& fluid=input.fluid;auto& trace=fluid.trace;
    AppIO::ConstantsUsage usage;
    const auto components=AppIO::readTwoPhaseInputComponents(fluid,documents,usage,threadsOverride);
    const auto& phases=documents.root("materials.dict").block("Materials").block("Phases").entries;
    for(std::size_t i=0;i<phases.size();++i){
        const auto& phase=phases[i];const int index=components.materials[i].phaseId==fluid.negativePhase?0:1;
        input.electric.electricProblem.epsilonByPhase[index]=boundedCaseValue(phase.get("permittivity"),trace,0,true,
            "Electric.electrostatic.material","code_permittivity");
        input.electric.conductivityByPhase[index]=boundedCaseValue(phase.get("conductivity"),trace,0,false,
            "Electric.OhmicCurrent.material","code_conductivity");
    }
    std::string chargeFormula;std::array<std::string,4> electricTypes,electricValues;

    documents.root("electric.dict").allow({"Electric"});
    const auto& electric=documents.root("electric.dict").block("Electric");
    electric.allow({"method","interfaceMotion","Initial","Boundary","TimeIntegration"});
    trace.record(electric.get("method"),electric.get("method").choice({"HHO"}),
        "dimensionless","Electric.electrostatic.dispatch");
    trace.record(electric.get("interfaceMotion"),electric.get("interfaceMotion").choice({"coupled"}),
        "dimensionless","Coupling.ElectrohydrodynamicStep.actualHydroMaterialTransport");
    const auto parseExpression=[&](const Node& node,const char* units,const char* consumer){
        const auto formula=trace.record(node,node.word(),units,consumer);
        usage.include(formula);
        try{
            return std::make_shared<const Expression>(formula,fluid.constants,Expression::Coordinates::Meridional);
        }catch(const std::exception& e){
            node.fail(std::string("a spatial r,z ScalarExpression using Fields.Constants, without t: ")+e.what());
        }
    };
    const auto sampleExpression=[&](const Node& node,const Expression& expression,Vector2D point){
        try{expression.sample(point.x,point.y);}
        catch(const std::exception& e){node.fail(std::string("finite spatial value/first partials at the input check points: ")+e.what());}
    };
    const auto& chargeInitial=electric.block("Initial");
    chargeInitial.allow({"surfaceChargeDensity"});
    const auto& chargeNode=chargeInitial.get("surfaceChargeDensity");
    input.chargeExpression_=parseExpression(chargeNode,"code_charge/code_surface_measure","Electric.Q.initialPanelIntegral");
    chargeFormula=chargeNode.word();
    for(double r:{fluid.domainMin.x,fluid.domainMax.x})for(double z:{fluid.domainMin.y,fluid.domainMax.y})
        sampleExpression(chargeNode,*input.chargeExpression_,{r,z});
    const auto& boundary=electric.block("Boundary");boundary.allow({"W","S","E","N"});
    bool anchor=false;
    for(std::size_t side=0;side<electricSides.size();++side){
        const auto& data=boundary.block(electricSides[side]);data.allow({"type","value"});
        electricTypes[side]=trace.record(data.get("type"),data.get("type").choice({"Dirichlet","Neumann"}),
            "dimensionless","Electric.electrostatic.outerBoundary");
        const auto type=electricTypes[side]=="Dirichlet"
            ?LAMP::Electric::ElectricCartesianBoundary::Type::Dirichlet:LAMP::Electric::ElectricCartesianBoundary::Type::Neumann;
        anchor=anchor || type==LAMP::Electric::ElectricCartesianBoundary::Type::Dirichlet;
        const auto formula=parseExpression(data.get("value"),electricTypes[side]=="Dirichlet"
            ?"code_potential":"code_permittivity*code_potential/code_length","Electric.electrostatic.outerBoundary");
        electricValues[side]=data.get("value").word();
        const Vector2D first=side==0 || side==1?fluid.domainMin:
            side==2?Vector2D{fluid.domainMax.x,fluid.domainMin.y}:Vector2D{fluid.domainMin.x,fluid.domainMax.y};
        const Vector2D second=side==0?Vector2D{fluid.domainMin.x,fluid.domainMax.y}:
            side==1?Vector2D{fluid.domainMax.x,fluid.domainMin.y}:fluid.domainMax;
        sampleExpression(data.get("value"),*formula,first);sampleExpression(data.get("value"),*formula,second);
        if(side==0){
            if(type!=LAMP::Electric::ElectricCartesianBoundary::Type::Neumann)
                data.get("type").fail("Neumann on the zero-area rz axis");
            if(data.get("value").number()!=0)
                data.get("value").fail("literal zero for homogeneous rz-axis Neumann flux");
        }
        input.electric.electricProblem.boundaries[side].type=type;
        input.electric.electricProblem.boundaries[side].value=[formula](Vector2D point){
            return formula->sample(point.x,point.y).value;
        };
    }
    if(!anchor)boundary.fail("at least one Dirichlet side to anchor the electric potential; no unanchored fallback");
    const auto& electricIntegration=electric.block("TimeIntegration");electricIntegration.allow({"method"});
    trace.record(electricIntegration.get("method"),electricIntegration.get("method").choice({"variableBdf2"}),
        "dimensionless","Electric.Q.BE startup then variable-step BDF2");
    input.electric.charge.method=LAMP::Electric::ElectricChargeTimeMethod::VariableBdf2;
    input.electric.hydro=fluid.step;
    electricDefaults(electric,trace,input.electric);

    usage.validate(documents.root("fields.dict").block("Fields"),fluid.constants,
        "a constant used by Geometry.Initial.phi, Electric.Initial.surfaceChargeDensity or an electric Boundary.value");
    AppIO::recordTwoPhaseInputDefaults(documents,fluid,trace);
    auto materialText=stream();materialText<<"Materials {\nPhases {\n";
    for(const auto& material:components.materials){
        const int index=material.phaseId==fluid.negativePhase?0:1;
        materialText<<material.prefix<<" permittivity "<<input.electric.electricProblem.epsilonByPhase[index]
            <<"; conductivity "<<input.electric.conductivityByPhase[index]<<"; }\n";
    }
    materialText<<components.materialInterfaces;
    auto electricText=stream();
    electricText<<"Electric { method HHO; interfaceMotion coupled; Initial { surfaceChargeDensity "
        <<std::quoted(chargeFormula)<<"; } Boundary {\n";
    for(std::size_t side=0;side<electricSides.size();++side)
        electricText<<electricSides[side]<<" { type "<<electricTypes[side]<<"; value "<<std::quoted(electricValues[side])<<"; }\n";
    electricText<<"} TimeIntegration { method variableBdf2; } }\n";
    const std::array<const char*,7> files{"run.dict","geometry.dict","mesh.dict","fields.dict","interface.dict","materials.dict","hydro.dict"};
    for(std::size_t i=0;i<files.size();++i)
        input.documents.push_back({documents.path(files[i]),documents.sourceText(files[i]),i==5?materialText.str():components.effective[i]});
    input.documents.push_back({documents.path("electric.dict"),documents.sourceText("electric.dict"),electricText.str()});
    input.trace=std::move(trace);
    return input;
}

