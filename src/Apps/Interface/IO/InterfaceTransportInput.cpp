#include "Interface/InterfaceTransport.h"
#include "Runtime/RunControls.h"
#include "Core/IO/StrictDictionary.h"
#include <filesystem>
#include <iomanip>
#include <regex>
#include <set>
#include <sstream>

namespace
{
using Node=DictionaryInput::Node;
using Trace=DictionaryInput::ConfigurationTrace;
using Config=InterfaceTransportConfig;
using Expression=Quantum::ScalarExpression;
Vector2D vector(const Node& node,Trace& trace,const char* units,const char* consumer)
{
    if(node.isBlock || node.values.size()!=4 || node.values.front()!="(" || node.values.back()!=")")
        node.fail("a finite two-component vector (first second)");
    Node scalar=node;scalar.values={node.values[1]};const double first=scalar.number();
    scalar.values={node.values[2]};const double second=scalar.number();
    std::ostringstream effective;effective<<std::setprecision(17)<<'('<<first<<' '<<second<<')';
    trace.record(node,effective.str(),units,consumer);return {first,second};
}
void expression(const Node& node,const Config& c,Expression::Coordinates coordinates,bool boundary=false)
{
    try{
        const Expression compiled(node.word(),c.constants,coordinates);
        if(coordinates==Expression::Coordinates::Time){compiled.sample(0);compiled.sample(c.endTime);return;}
        for(double time:{0.,boundary?c.endTime:0.})
            for(double x:{c.domainMin.x,c.domainMax.x})for(double y:{c.domainMin.y,c.domainMax.y})compiled.sample(x,y,time);
    }catch(const std::exception& error){node.fail(std::string("a valid analytic expression for Interface.FirstJet: ")+error.what());}
}
void inputFiles()
{
    std::set<std::string> names;for(const char* file:interfaceTransportInputFiles)names.insert(std::filesystem::path(file).filename().string());
    for(const auto& entry:std::filesystem::directory_iterator("system"))
        if(!entry.is_regular_file() || !names.count(entry.path().filename().string()))
            throw std::invalid_argument("unused standalone Interface input: "+entry.path().string());
}
void readRun(const Node& root,Config& c,Trace& trace)
{
    root.allow({"Run","Time","Output"});const auto& run=root.block("Run");run.allow({"app"});
    trace.record(run.get("app"),run.get("app").choice({"Interface"}),"dimensionless","Runtime.App.dispatch");
    const auto& time=root.block("Time");time.allow({"units","dt","endTime"});
    trace.record(time.get("units"),time.get("units").choice({"code"}),"dimensionless","Runtime.clock.units");
    c.dt=trace.record(time.get("dt"),time.get("dt").number(),"code_time","Runtime.clock");
    c.endTime=trace.record(time.get("endTime"),time.get("endTime").number(),"code_time","Runtime.clock");
    if(!(c.dt>0))time.get("dt").fail("finite and positive dt for Interface.FirstJet");
    if(!(c.endTime>0))time.get("endTime").fail("finite and positive endTime for Interface.FirstJet");
    RunControls controls;controls.dt=c.dt;controls.endTime=c.endTime;
    const auto steps=controls.fixedStepCount();
    if(!steps || *steps>1000000)time.get("endTime").fail("an integer number of fixed dt steps in [1,1000000] for Interface.FirstJet");
    c.steps=*steps;
    const auto& output=root.block("Output");output.allow({"Fields"});
    const auto& fields=output.block("Fields");fields.allow({"enabled","writeEverySteps"});
    c.writeFields=trace.record(fields.get("enabled"),fields.get("enabled").boolean(),"dimensionless","Runtime.output.fields");
    if(c.writeFields)c.writeEverySteps=fields.find("writeEverySteps")?trace.record(fields.get("writeEverySteps"),fields.get("writeEverySteps").integer(),"steps","Runtime.output.schedule"):
        trace.defaultValue(fields,"writeEverySteps",1,"steps","Runtime.output.schedule");
    else if(fields.find("writeEverySteps"))fields.get("writeEverySteps").fail("omitted when enabled is false (unused schedule)");
}
void readGeometry(const Node& root,Config& c,Trace& trace)
{
    root.allow({"coordinates","Domain","Initial"});
    c.radial=trace.record(root.get("coordinates"),root.get("coordinates").choice({"xy","rz"}),"dimensionless","Geometry.measure")=="rz";
    const auto& domain=root.block("Domain");domain.allow({"min","max"});
    c.domainMin=vector(domain.get("min"),trace,"code_length","Geometry.domain");
    c.domainMax=vector(domain.get("max"),trace,"code_length","Geometry.domain");
    const auto& initial=root.block("Initial");initial.allow({"phi"});
    c.initialPhi=trace.record(initial.get("phi"),initial.get("phi").word(),"code_phi","Interface.FirstJet.initialData");
}
void readMesh(const Node& root,Config& c,Trace& trace)
{
    root.allow({"resolution","initialMesh","Refinement","Regrid"});
    c.resolution=trace.record(root.get("resolution"),root.get("resolution").integer(),"dimensionless","Mesh.baseGrid");
    c.initialAmr=(root.find("initialMesh")?trace.record(root.get("initialMesh"),root.get("initialMesh").choice({"uniform","amr"}),"dimensionless","Mesh.initialGrid"):
        trace.defaultValue(root,"initialMesh",std::string("uniform"),"dimensionless","Mesh.initialGrid"))=="amr";
    if(root.find("Regrid")){
        const auto& block=root.block("Regrid");block.allow({"atSteps"});
        c.regridSteps=trace.record(block.get("atSteps"),block.get("atSteps").ids(),"steps","Mesh.transactions");
    }else trace.defaultValue(root,"Regrid.atSteps",std::vector<int>{},"steps","Mesh.transactions");
    if(c.initialAmr || !c.regridSteps.empty()){
        const auto& block=root.block("Refinement");block.allow({"min","max","level"});
        c.refineMin=vector(block.get("min"),trace,"code_length","Mesh.refinement");
        c.refineMax=vector(block.get("max"),trace,"code_length","Mesh.refinement");
        c.refinementLevel=trace.record(block.get("level"),block.get("level").integer(),"dimensionless","Mesh.refinement");
    }else if(root.find("Refinement"))root.get("Refinement").fail("omitted when initialMesh=uniform and no Regrid (unused refinement)");
}
void readFields(const Node& root,Config& c,Trace& trace,std::vector<std::string>& formulas)
{
    root.allow({"Constants","Velocity"});
    if(root.find("Constants"))for(const auto& node:root.block("Constants").entries){
        if(std::set<std::string>{"x","y","r","z","t","s","pi","sin","cos","tan","exp","log","sqrt"}.count(node.name))node.fail("a nonreserved constant name");
        const double value=node.number();
        try{Expression test("0",{{node.name,value}},Expression::Coordinates::Time);}catch(const std::exception& e){node.fail(e.what());}
        c.constants.emplace(node.name,trace.record(node,value,"code","Fields.expressionParameters"));
    }
    const auto& velocity=root.block("Velocity");
    c.velocity=trace.record(velocity.get("type"),velocity.get("type").choice({"constant","polynomial"}),"dimensionless","Fields.prescribedVelocity")=="constant"?Config::Velocity::Constant:Config::Velocity::Polynomial;
    if(c.velocity==Config::Velocity::Constant){velocity.allow({"type","value"});c.constantVelocity=vector(velocity.get("value"),trace,"code_velocity","Fields.prescribedVelocity");}
    else{
        velocity.allow({"type","origin","scale","Coefficients"});
        c.velocityOrigin=vector(velocity.get("origin"),trace,"code_length","Fields.prescribedVelocity");
        c.velocityScale=vector(velocity.get("scale"),trace,"code_length","Fields.prescribedVelocity");
        if(!(c.velocityScale.x>0 && c.velocityScale.y>0))velocity.get("scale").fail("positive polynomial coordinate scales");
        const auto& coefficients=velocity.block("Coefficients");
        if(coefficients.entries.empty())coefficients.fail("at least one coefficient c00..c33; use type constant for a zero field");
        for(const auto& node:coefficients.entries){
            if(node.name.size()!=3 || node.name[0]!='c' || node.name[1]<'0' || node.name[1]>'3' || node.name[2]<'0' || node.name[2]>'3')node.fail("coefficient block c00..c33, index ij means xi^i eta^j");
            if(!node.isBlock)node.fail("a coefficient block with x and y time expressions");node.allow({"x","y"});
            auto& item=c.coefficients[4*(node.name[1]-'0')+node.name[2]-'0'];item.specified=true;
            item.x=trace.record(node.get("x"),node.get("x").word(),"code_velocity","Fields.prescribedVelocity");
            item.y=trace.record(node.get("y"),node.get("y").word(),"code_velocity","Fields.prescribedVelocity");
            for(const char* component:{"x","y"}){expression(node.get(component),c,Expression::Coordinates::Time);formulas.push_back(node.get(component).word());}
        }
        for(int k=0;k<16;++k)if(!c.coefficients[k].specified){const auto name="c"+std::to_string(k/4)+std::to_string(k%4);
            trace.defaultValue(coefficients,name+".x",std::string("0"),"code_velocity","Fields.prescribedVelocity");
            trace.defaultValue(coefficients,name+".y",std::string("0"),"code_velocity","Fields.prescribedVelocity");}
    }
}
}
InterfaceTransportInput readInterfaceTransportInput()
{
    DictionaryInput::Document interfaceDocument("system/interface.dict");interfaceDocument.root().allow({"Interface"});
    const auto& interface=interfaceDocument.root().block("Interface");
    if(interface.find("Benchmark"))interface.get("Benchmark").fail("five-file Interface input; see docs/api/G04_Interface_Dictionary_Schema.md");
    interface.allow({"method","Boundary"});
    InterfaceTransportInput input;auto& c=input.config;auto& trace=input.trace;
    trace.record(interface.get("method"),interface.get("method").choice({"hermiteGals"}),"dimensionless","Interface.FirstJet.dispatch");
    DictionaryInput::Document runDocument("system/run.dict"),geometryDocument("system/geometry.dict"),meshDocument("system/mesh.dict"),fieldsDocument("system/fields.dict");
    geometryDocument.root().allow({"Geometry"});meshDocument.root().allow({"Mesh"});fieldsDocument.root().allow({"Fields"});
    const auto& geometry=geometryDocument.root().block("Geometry");const auto& mesh=meshDocument.root().block("Mesh");const auto& fields=fieldsDocument.root().block("Fields");
    readRun(runDocument.root(),c,trace);readGeometry(geometry,c,trace);readMesh(mesh,c,trace);
    try{c.validate();}catch(const InterfaceTransportConfigError& error){
        const auto& node=[&]() -> const Node& {using F=Config::Field;switch(error.field){
            case F::Resolution:return mesh.get("resolution");case F::Steps:case F::EndTime:return runDocument.root().block("Time").get("endTime");
            case F::Dt:return runDocument.root().block("Time").get("dt");case F::WriteEverySteps:return runDocument.root().block("Output").block("Fields").get("writeEverySteps");
            case F::RegridSteps:return mesh.block("Regrid").get("atSteps");case F::Domain:return geometry.block("Domain");
            case F::Refinement:return mesh.get("Refinement");}throw std::logic_error("unmapped transport input field");}();
        node.fail(std::string("for method=hermiteGals / Interface.FirstJet: ")+error.what());
    }
    std::vector<std::string> formulas{c.initialPhi};readFields(fields,c,trace,formulas);
    expression(geometry.block("Initial").get("phi"),c,c.radial?Expression::Coordinates::Meridional:Expression::Coordinates::Planar);
    const auto& boundary=interface.block("Boundary");
    c.prescribedBoundary=trace.record(boundary.get("type"),boundary.get("type").choice({"none","prescribed"}),"dimensionless","Interface.FirstJet.boundaryData")=="prescribed";
    if(c.prescribedBoundary){boundary.allow({"type","phi"});
        c.boundaryPhi=trace.record(boundary.get("phi"),boundary.get("phi").word(),"code_phi","Interface.FirstJet.boundaryData");formulas.push_back(c.boundaryPhi);
        expression(boundary.get("phi"),c,c.radial?Expression::Coordinates::MeridionalTime:Expression::Coordinates::PlanarTime,true);
    }else boundary.allow({"type"});
    std::set<std::string> identifiers;const std::regex token("(?:^|[^A-Za-z_0-9.])([A-Za-z_][A-Za-z_0-9]*)");
    for(const auto& formula:formulas)for(auto it=std::sregex_iterator(formula.begin(),formula.end(),token);it!=std::sregex_iterator();++it)identifiers.insert((*it)[1].str());
    for(const auto& item:c.constants)if(!identifiers.count(item.first))fields.block("Constants").get(item.first).fail("a constant used by Initial, Boundary or Velocity (unused constant)");
    inputFiles();return input;
}
