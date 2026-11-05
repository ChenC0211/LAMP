#include "Electric/IO/ElectricCaseInput.h"
#include "Core/IO/StrictDictionary.h"
#include "Core/Units/UnitsContext.h"
#include "IO/DictionaryParallel.h"
#include "IO/DictionaryRunControls.h"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>

namespace
{
using Node = DictionaryInput::Node;
using Trace = DictionaryInput::ConfigurationTrace;
using Expression = LAMP::ScalarExpression;
constexpr std::array<const char*,5> files{
    "run.dict","geometry.dict","mesh.dict","materials.dict","electric.dict"};
constexpr std::array<const char*,4> sides{"W","S","E","N"};
std::ostringstream stream() { std::ostringstream out; out << std::setprecision(17); return out; }
std::string text(const std::filesystem::path& path)
{
    std::ifstream in(path,std::ios::binary);
    if(!in) throw std::runtime_error("cannot open Electric input: "+path.string());
    std::ostringstream out; out << in.rdbuf();
    if(in.bad()) throw std::runtime_error("cannot read Electric input: "+path.string());
    return out.str();
}
void vectorText(std::ostream& out,Vector2D point) { out << '(' << point.x << ' ' << point.y << ')'; }
Vector2D vector(const Node& node,Trace& trace)
{
    if(node.isBlock || node.values.size()!=4 || node.values.front()!="(" || node.values.back()!=")")
        node.fail("a finite two-component vector (first second)");
    auto item=node; item.values={node.values[1]}; const double first=item.number();
    item.values={node.values[2]}; const double second=item.number();
    auto out=stream(); vectorText(out,{first,second});
    trace.record(node,out.str(),"code_length","Geometry.domain");
    return {first,second};
}
std::array<int,2> resolution(const Node& node,Trace& trace)
{
    if(node.isBlock || node.values.size()!=4 || node.values.front()!="(" || node.values.back()!=")")
        node.fail("two explicit integer cell counts (nx ny)");
    auto item=node; item.values={node.values[1]}; const int nx=item.integer();
    item.values={node.values[2]}; const int ny=item.integer();
    if(nx<4 || ny<4 || static_cast<long long>(nx)*ny>1000000)
        node.fail("nx,ny >= 4 and at most 1000000 fixed uniform cells; no AMR");
    trace.record(node,std::vector<int>{nx,ny},"cells","Mesh.uniformGrid");
    return {nx,ny};
}
double material(const Node& node,Trace& trace,bool positive,const char* units,const char* consumer)
{
    const double value=node.number();
    if(positive?!(value>0):value<0)
        node.fail(positive?"a finite positive material coefficient":"a finite nonnegative material coefficient (zero supported)");
    return trace.record(node,value,units,consumer);
}
std::shared_ptr<const Expression> expression(const Node& node,bool radial,Trace& trace,const char* units,const char* consumer)
{
    const auto formula=trace.record(node,node.word(),units,consumer);
    try {
        return std::make_shared<const Expression>(formula,Expression::Constants{},
            radial?Expression::Coordinates::Meridional:Expression::Coordinates::Planar);
    } catch(const std::exception& error) {
        node.fail(std::string("a spatial ScalarExpression using ")+(radial?"r,z":"x,y")
            +", numeric literals and pi, without t or unknown constants: "+error.what());
    }
}
void sample(const Node& node,const Expression& expression,Vector2D point)
{
    try { expression.sample(point.x,point.y); }
    catch(const std::exception& error) { node.fail(std::string("finite spatial value/first partials at the input check points: ")+error.what()); }
}
void kernelDefaults(const Node& electric,Trace& trace,const ElectricCaseInput& input)
{
    const auto record=[&](const char* key,double value,const char* consumer) {
        trace.defaultValue(electric,key,value,"dimensionless",consumer);
    };
    record("KernelDefaults.linear.relativeLinearTolerance",input.linear.relativeLinearTolerance,"Electric.electrostatic.linear");
    record("KernelDefaults.linear.relativeSymmetryTolerance",input.linear.relativeSymmetryTolerance,"Electric.electrostatic.symmetry");
    record("KernelDefaults.linear.relativeConservationTolerance",input.linear.relativeConservationTolerance,"Electric.electrostatic.Gauss");
    record("KernelDefaults.charge.relativeSchurTolerance",input.charge.relativeSchurTolerance,"Electric.Q.Schur");
    record("KernelDefaults.charge.relativeChargeBalanceTolerance",input.charge.relativeChargeBalanceTolerance,"Electric.Q.chargeBalance");
    record("KernelDefaults.charge.maximumStepRatio",input.charge.maximumStepRatio,"Electric.Q.variableBdf2");
    trace.defaultValue(electric,"KernelDefaults.charge.maximumIterations",input.charge.maximumIterations,"iterations","Electric.Q.Schur");
    trace.defaultValue(electric,"KernelDefaults.charge.restartLength",input.charge.restartLength,"iterations","Electric.Q.Schur");
    record("KernelDefaults.geometry.relativeMatrixTolerance",input.geometry.relativeMatrixTolerance,"Electric.electrostatic.phaseQuadrature");
    record("KernelDefaults.geometry.minimumRelativePhaseArea",input.geometry.minimumRelativePhaseArea,"Electric.electrostatic.smallCutRejection");
    trace.defaultValue(electric,"KernelDefaults.geometry.maximumRefinements",input.geometry.maximumRefinements,"levels","Electric.electrostatic.phaseQuadrature");
    trace.defaultValue(electric,"KernelDefaults.geometry.maximumPoints",std::to_string(input.geometry.maximumPoints),"points","Electric.electrostatic.phaseQuadrature");
    trace.defaultValue(electric,"KernelDefaults.geometry.maximumSliceQueriesPerCell",std::to_string(input.geometry.maximumSliceQueriesPerCell),"queries","Electric.electrostatic.phaseQuadrature");
    trace.defaultValue(electric,"KernelDefaults.geometry.maximumDepth",input.geometry.maximumDepth,"levels","Electric.electrostatic.phaseQuadrature");
    record("KernelDefaults.moments.relativeQuadratureTolerance",input.moments.relativeQuadratureTolerance,"Electric.P_Q");
    trace.defaultValue(electric,"KernelDefaults.moments.maximumAdditionalRefinements",input.moments.maximumAdditionalRefinements,"levels","Electric.P_Q");
    record("KernelDefaults.maxwell.relativeProjectionTolerance",input.maxwell.relativeProjectionTolerance,"Electric.Maxwell.projection");
    record("KernelDefaults.maxwell.relativeChargeJumpTolerance",input.maxwell.relativeChargeJumpTolerance,"Electric.Maxwell.jump");
    record("KernelDefaults.firstJet.relativeInterfaceResidualTolerance",input.publication.relativeInterfaceResidualTolerance,"Interface.FirstJet.publication");
    record("KernelDefaults.firstJet.unitNormalTolerance",input.publication.unitNormalTolerance,"Interface.FirstJet.publication");
    record("KernelDefaults.firstJet.relativeLedgerTolerance",input.publication.volume.relativeLedgerTolerance,"Interface.FirstJet.phaseMeasure");
    record("KernelDefaults.firstJet.relativeMomentEstimateTolerance",input.publication.volume.relativeMomentEstimateTolerance,"Interface.FirstJet.phaseMeasure");
    trace.defaultValue(electric,"KernelDefaults.bulkSource",std::string("zero"),"code_charge/code_volume","Electric.electrostatic.quasineutralBulk");
}
}

LAMP::ScalarExpression::Jet ElectricCaseInput::initialJet(Vector2D point) const
{
    if(!initialPhi_) throw std::logic_error("Electric initial geometry has not been validated");
    return initialPhi_->sample(point.x,point.y);
}
double ElectricCaseInput::initialChargeDensity(Vector2D point) const
{
    if(!initialCharge_) throw std::logic_error("Electric initial surface charge has not been validated");
    return initialCharge_->sample(point.x,point.y).value;
}

ElectricCaseInput readElectricCaseInput(std::optional<int> threadsOverride)
{
    namespace fs=std::filesystem;
    const fs::path directory("system");
    const std::set<std::string> names(files.begin(),files.end());
    if(!fs::is_directory(directory)) throw std::invalid_argument("missing formal Electric system directory: system");
    for(const auto& entry:fs::directory_iterator(directory))
        if(!entry.is_regular_file() || entry.is_symlink() || !names.count(entry.path().filename().string()))
            throw std::invalid_argument("unused formal Electric input: "+entry.path().string()
                +"; expected only run/geometry/mesh/materials/electric.dict");
    for(const auto* name:files) if(!fs::is_regular_file(directory/name))
        throw std::invalid_argument("missing formal Electric input: "+(directory/name).string());
    const auto path=[&](const char* name) { return (directory/name).generic_string(); };
    DictionaryInput::Document runDoc(path(files[0])),geoDoc(path(files[1])),meshDoc(path(files[2])),
        materialDoc(path(files[3])),electricDoc(path(files[4]));
    ElectricCaseInput input; auto& trace=input.trace;
    const auto& run=runDoc.root(); run.allow({"Run","Time","Output","Log","Restart","Parallel"});
    const auto& dispatch=run.block("Run"); dispatch.allow({"app"});
    trace.record(dispatch.get("app"),dispatch.get("app").choice({"Electric"}),"dimensionless","Runtime.App.dispatch");
    run.block("Time").get("units").choice({"code"});
    input.parallelThreads=AppRuntime::readParallelThreads(run,&trace,threadsOverride);
    input.run=readRunControls(run,UnitsContext::MakeCodeUnits(),trace);
    if(input.run.fixedStepCount().value()>=std::numeric_limits<int>::max())
        run.block("Time").get("endTime").fail("at most INT_MAX-1 fixed steps plus the initial field evaluation");
    if(input.run.diagnosticFields)
        run.block("Output").block("Fields").get("profile").fail("physical; this electric writer has no additional diagnostic profile");

    geoDoc.root().allow({"Geometry"}); const auto& geometry=geoDoc.root().block("Geometry");
    geometry.allow({"coordinates","Domain","Phases","Initial"});
    input.radial=trace.record(geometry.get("coordinates"),geometry.get("coordinates").choice({"xy","rz"}),
        "dimensionless","Geometry.measure")=="rz";
    const auto& domain=geometry.block("Domain"); domain.allow({"min","max"});
    input.domainMin=vector(domain.get("min"),trace); input.domainMax=vector(domain.get("max"),trace);
    const double width=input.domainMax.x-input.domainMin.x,height=input.domainMax.y-input.domainMin.y;
    if(!(std::isfinite(width) && std::isfinite(height) && width>0 && height>0)
        || (input.radial && input.domainMin.x!=0))
        domain.fail("finite ordered rectangular bounds; rz starts at radial minimum zero");
    const auto& phases=geometry.block("Phases"); phases.allow({"negative","positive"});
    input.negativePhase=trace.record(phases.get("negative"),phases.get("negative").integer(),"1-based ID","Geometry.negativePhase");
    input.positivePhase=trace.record(phases.get("positive"),phases.get("positive").integer(),"1-based ID","Geometry.positivePhase");
    if(input.negativePhase<1 || input.positivePhase<1 || input.negativePhase==input.positivePhase)
        phases.fail("two distinct positive phase IDs; negative/positive phi bind electric phase indices 0/1");
    const auto& initial=geometry.block("Initial"); initial.allow({"phi"});
    input.initialPhi_=expression(initial.get("phi"),input.radial,trace,"code_phi","Interface.FirstJet.initialData");
    for(double first:{input.domainMin.x,input.domainMax.x})
        for(double second:{input.domainMin.y,input.domainMax.y}) sample(initial.get("phi"),*input.initialPhi_,{first,second});
    if(input.radial) {
        // This is a small preflight regularity check, not a mesh/topology proof.
        // The actual publication remains responsible for all native checks.
        const double probe=width*.01;
        for(double z:{input.domainMin.y,(input.domainMin.y+input.domainMax.y)*.5,input.domainMax.y}) {
            try {
                const auto axis=input.initialJet({0,z}),plus=input.initialJet({probe,z}),minus=input.initialJet({-probe,z});
                const double valueScale=std::max({1.,std::abs(plus.value),std::abs(minus.value)});
                const double radialScale=std::max({1.,std::abs(plus.derivative),std::abs(minus.derivative)});
                const double axialScale=std::max({1.,std::abs(plus.partialSecond),std::abs(minus.partialSecond)});
                if(std::abs(axis.derivative)>1e-12*radialScale
                    || std::abs(plus.value-minus.value)>1e-12*valueScale
                    || std::abs(plus.derivative+minus.derivative)>1e-12*radialScale
                    || std::abs(plus.partialSecond-minus.partialSecond)>1e-12*axialScale)
                    initial.get("phi").fail("an even radial FirstJet with zero radial derivative on rz axis");
            } catch(const std::exception& error) {
                initial.get("phi").fail(std::string("an even radial FirstJet evaluable across rz axis: ")+error.what());
            }
        }
        // The RZ axis already has its own native symmetry identity.
    }
    meshDoc.root().allow({"Mesh"}); const auto& mesh=meshDoc.root().block("Mesh");
    mesh.allow({"resolution","initialMesh"}); const auto cells=resolution(mesh.get("resolution"),trace);
    input.nx=cells[0]; input.ny=cells[1];
    if(const auto* kind=mesh.find("initialMesh"))
        trace.record(*kind,kind->choice({"uniform"}),"dimensionless","Mesh.initialGrid");
    else trace.defaultValue(mesh,"initialMesh",std::string("uniform"),"dimensionless","Mesh.initialGrid");

    materialDoc.root().allow({"Materials"}); const auto& materials=materialDoc.root().block("Materials");
    materials.allow({"Phases"}); const auto& materialPhases=materials.block("Phases");
    if(materialPhases.entries.size()!=2) materialPhases.fail("exactly the two Geometry.Phases material blocks");
    std::set<int> materialIds;
    for(const auto& phase:materialPhases.entries) {
        if(!phase.isBlock) phase.fail("a named material block");
        phase.allow({"phaseId","permittivity","conductivity"});
        const int id=trace.record(phase.get("phaseId"),phase.get("phaseId").integer(),"1-based ID","Materials.geometryBinding");
        if((id!=input.negativePhase && id!=input.positivePhase) || !materialIds.insert(id).second)
            phase.get("phaseId").fail("one distinct Geometry.Phases ID");
        const int index=id==input.negativePhase?0:1;
        input.electric.epsilonByPhase[index]=material(phase.get("permittivity"),trace,true,"code_permittivity","Electric.electrostatic.material");
        input.conductivity[index]=material(phase.get("conductivity"),trace,false,"code_conductivity","Electric.OhmicCurrent.material");
    }
    electricDoc.root().allow({"Electric"}); const auto& electric=electricDoc.root().block("Electric");
    electric.allow({"method","interfaceMotion","Initial","Boundary","TimeIntegration"});
    trace.record(electric.get("method"),electric.get("method").choice({"HHO"}),"dimensionless","Electric.electrostatic.dispatch");
    trace.record(electric.get("interfaceMotion"),electric.get("interfaceMotion").choice({"fixed"}),
        "dimensionless","Electric.fixedGeometry; moving geometry requires the Hydro coupling owner");
    const auto& chargeInitial=electric.block("Initial"); chargeInitial.allow({"surfaceChargeDensity"});
    input.initialCharge_=expression(chargeInitial.get("surfaceChargeDensity"),input.radial,trace,
        "code_charge/code_surface_measure","Electric.Q.initialPanelIntegral");
    for(double first:{input.domainMin.x,input.domainMax.x})
        for(double second:{input.domainMin.y,input.domainMax.y})
            sample(chargeInitial.get("surfaceChargeDensity"),*input.initialCharge_,{first,second});
    const auto& boundary=electric.block("Boundary"); boundary.allow({"W","S","E","N"});
    bool anchor=false; std::array<std::string,4> types,values;
    for(std::size_t side=0;side<sides.size();++side) {
        const auto& data=boundary.block(sides[side]); data.allow({"type","value"});
        types[side]=trace.record(data.get("type"),data.get("type").choice({"Dirichlet","Neumann"}),
            "dimensionless","Electric.electrostatic.outerBoundary");
        input.electric.boundaries[side].type=types[side]=="Dirichlet"
            ?LAMP::Electric::ElectricCartesianBoundary::Type::Dirichlet:LAMP::Electric::ElectricCartesianBoundary::Type::Neumann;
        anchor=anchor || types[side]=="Dirichlet";
        const auto formula=expression(data.get("value"),input.radial,trace,
            types[side]=="Dirichlet"?"code_potential":"code_permittivity*code_potential/code_length",
            "Electric.electrostatic.outerBoundary");
        values[side]=data.get("value").word();
        const Vector2D first=side==0?input.domainMin:side==1?input.domainMin:
            side==2?Vector2D{input.domainMax.x,input.domainMin.y}:Vector2D{input.domainMin.x,input.domainMax.y};
        const Vector2D second=side==0?Vector2D{input.domainMin.x,input.domainMax.y}:
            side==1?Vector2D{input.domainMax.x,input.domainMin.y}:input.domainMax;
        sample(data.get("value"),*formula,first); sample(data.get("value"),*formula,second);
        if(input.radial && side==0) {
            if(types[side]!="Neumann") data.get("type").fail("Neumann on the zero-area rz axis");
            // A literal guarantees homogeneous flux everywhere on the axis;
            // finite sampling would not certify a general expression is zero.
            if(data.get("value").number()!=0) data.get("value").fail("literal zero for homogeneous rz-axis Neumann flux");
        }
        input.electric.boundaries[side].value=[formula](Vector2D point) { return formula->sample(point.x,point.y).value; };
    }
    if(!anchor) boundary.fail("at least one Dirichlet side to anchor the electric potential; no unanchored fallback");
    const auto& integration=electric.block("TimeIntegration"); integration.allow({"method"});
    trace.record(integration.get("method"),integration.get("method").choice({"variableBdf2"}),
        "dimensionless","Electric.Q.BE startup then variable-step BDF2");
    input.charge.method=LAMP::Electric::ElectricChargeTimeMethod::VariableBdf2;
    kernelDefaults(electric,trace,input);

    auto runText=stream(),geoText=stream(),meshText=stream(),materialsText=stream(),electricText=stream();
    const auto boolean=[](bool value) { return value?"true":"false"; };
    runText << "Run { app Electric; }\nParallel { threads " << input.parallelThreads << "; }\n"
        << "Time { units code; dt " << input.run.dt << "; endTime " << input.run.endTime << "; }\n"
        << "Output { prefix " << std::quoted(input.run.outputPrefix) << "; Fields { enabled " << boolean(input.run.writeFields)
        << "; format tecplot; profile physical; writeControl " << (input.run.writeEverySteps?"stepInterval":"timeInterval")
        << "; writeInterval ";
    if(input.run.writeEverySteps) runText << input.run.writeEverySteps; else runText << input.run.writeInterval;
    runText << "; writeInitial " << boolean(input.run.writeInitial) << "; writeFirst " << boolean(input.run.writeFirst)
        << "; writeFinal " << boolean(input.run.writeFinal) << "; } Series { profile "
        << (input.run.diagnosticSeries?"diagnostic":"physical") << "; } }\nLog { everySteps " << input.run.logEverySteps
        << "; level " << (input.run.detailedLog?"detail":"summary") << "; }\nRestart { model none; }\n";
    geoText << "Geometry { coordinates " << (input.radial?"rz":"xy") << "; Domain { min ";
    vectorText(geoText,input.domainMin); geoText << "; max "; vectorText(geoText,input.domainMax);
    geoText << "; } Phases { negative " << input.negativePhase << "; positive " << input.positivePhase
        << "; } Initial { phi " << std::quoted(initial.get("phi").word()) << "; } }\n";
    meshText << "Mesh { resolution (" << input.nx << ' ' << input.ny << "); initialMesh uniform; }\n";
    materialsText << "Materials { Phases {\n";
    for(const auto& phase:materialPhases.entries) {
        const int id=phase.get("phaseId").integer(),index=id==input.negativePhase?0:1;
        materialsText << phase.name << " { phaseId " << id << "; permittivity " << input.electric.epsilonByPhase[index]
            << "; conductivity " << input.conductivity[index] << "; }\n";
    }
    materialsText << "} }\n";
    electricText << "Electric { method HHO; interfaceMotion fixed; Initial { surfaceChargeDensity "
        << std::quoted(chargeInitial.get("surfaceChargeDensity").word()) << "; } Boundary {\n";
    for(std::size_t side=0;side<sides.size();++side)
        electricText << sides[side] << " { type " << types[side] << "; value " << std::quoted(values[side]) << "; }\n";
    electricText << "} TimeIntegration { method variableBdf2; } }\n";
    const std::array<std::string,5> effective{runText.str(),geoText.str(),meshText.str(),materialsText.str(),electricText.str()};
    for(std::size_t file=0;file<files.size();++file)
        input.documents.push_back({path(files[file]),text(directory/files[file]),effective[file]});
    return input;
}
