#include "Modules/Interface/Compatibility/IO/DictionaryInterface.h"
#include "Core/IO/ConfigurationTrace.h"

#include <set>

InterfaceAdaptationControls readInterfaceAdaptation(const DictionaryInput::Node& a,
    double distanceBandCells,DictionaryInput::ConfigurationTrace& trace,bool allowCoarsening)
{
    a.allow({"maximumLevel","cellBudget","maxCurvatureCells","minimumGapCells","maximumGeometryFractionDifference",
        "remapScalarTolerance","remapRelativeGeometryTolerance","reserveBandCells","Coarsening"});
    InterfaceAdaptationControls result;auto& options=result.policy;
    options.maximumLevel=trace.record(a.get("maximumLevel"),a.get("maximumLevel").integer(),"level","Interface.ResolutionPlan");
    result.cellBudget=trace.record(a.get("cellBudget"),a.get("cellBudget").integer(),"cells","Interface.ResolutionPlan.refine");
    const auto number=[&](const char* key,double value,const char* consumer){return a.find(key)?
        trace.record(a.get(key),a.get(key).number(),"dimensionless",consumer):trace.defaultValue(a,key,value,"dimensionless",consumer);};
    options.maxCurvatureCells=number("maxCurvatureCells",options.maxCurvatureCells,"Interface.ResolutionPlan");
    options.minimumOpposedDistanceCells=number("minimumGapCells",options.minimumOpposedDistanceCells,"Interface.ResolutionPlan");
    options.maximumGeometryFractionDifference=number("maximumGeometryFractionDifference",options.maximumGeometryFractionDifference,"Interface.ResolutionPlan");
    options.distanceBandCells=distanceBandCells;
    options.reserveBandCells=number("reserveBandCells",options.reserveBandCells,"Interface.ResolutionPlan");
    result.remapScalarTolerance=number("remapScalarTolerance",result.remapScalarTolerance,"Interface.MeshTransfer");
    result.remapRelativeGeometryTolerance=number("remapRelativeGeometryTolerance",result.remapRelativeGeometryTolerance,"Interface.MeshTransfer");
    if(options.maximumLevel<1 || options.maximumLevel>8 || result.cellBudget<=0 || result.cellBudget>2000000
        || !(options.maxCurvatureCells>0 && options.maxCurvatureCells<=.5)
        || !(options.minimumOpposedDistanceCells>=2) || !(options.maximumGeometryFractionDifference>0 && options.maximumGeometryFractionDifference<=1)
        || !(options.reserveBandCells>=0 && options.reserveBandCells<=64) || !(distanceBandCells>0)
        || !(result.remapScalarTolerance>=0) || !(result.remapRelativeGeometryTolerance>0 && result.remapRelativeGeometryTolerance<=1))
        throw std::invalid_argument("invalid Interface adaptation levels, budget or numerical limits");
    if(a.find("Coarsening")){
        if(!allowCoarsening)throw std::invalid_argument("Coarsening requires the complete CCH parametric Hydro consumer");
        const auto& block=a.block("Coarsening");block.allow({"hysteresis","relativeVelocityError"});
        InterfaceAdaptationControls::CoarseningControls coarsening;
        const auto value=[&](const char* key,double fallback,const char* consumer){return block.find(key)?
            trace.record(block.get(key),block.get(key).number(),"dimensionless",consumer):trace.defaultValue(block,key,fallback,"dimensionless",consumer);};
        coarsening.hysteresis=value("hysteresis",coarsening.hysteresis,"Interface.ResolutionPlan.coarsening");
        coarsening.relativeVelocityError=value("relativeVelocityError",coarsening.relativeVelocityError,"Hydro.CompleteVelocityTransfer");
        if(!(coarsening.hysteresis>0&&coarsening.hysteresis<1)
            ||!(coarsening.relativeVelocityError>0&&coarsening.relativeVelocityError<1))
            throw std::invalid_argument("coarsening hysteresis and relativeVelocityError must be in (0,1)");
        result.coarsening=coarsening;
    }
    return result;
}

namespace
{
InterfaceNodalConfig::ScalarBoundary readInterfaceBoundary(DictionaryInput::Reader& input)
{
    InterfaceNodalConfig::ScalarBoundary result;
    input.expect("{");
    std::set<std::string> keys;
    while (input.peek() != "}")
    {
        const auto key = input.next();
        if (!keys.insert(key).second) input.fail("duplicate interfaceBoundary key '" + key + "'");
        if (key == "model") result.model = input.next();
        else if (key == "fixedContacts")
        {
            input.expect("(");
            while (input.peek() != ")")
            {
                input.expect("(");
                const double x = input.number("fixedContacts.x");
                if (input.peek() == ",") input.next();
                const double y = input.number("fixedContacts.y");
                input.expect(")");
                for (const auto& old : result.fixedContacts)
                    if (old.x == x && old.y == y) input.fail("duplicate fixed contact");
                result.fixedContacts.emplace_back(x, y);
                if (input.peek() == ",") input.next();
            }
            input.expect(")");
        }
        else if (key == "inflow")
        {
            input.expect("{");
            std::set<std::string> names;
            while (input.peek() != "}")
            {
                InterfaceNodalConfig::Inflow inflow;
                inflow.patch = input.next();
                if (!names.insert(inflow.patch).second) input.fail("duplicate interface inflow patch '" + inflow.patch + "'");
                input.expect("{");
                std::set<std::string> fields;
                while (input.peek() != "}")
                {
                    const auto field = input.next();
                    if (!fields.insert(field).second) input.fail("duplicate interface inflow key '" + field + "'");
                    if (field == "offset") inflow.offset = input.number(field);
                    else if (field == "slope") inflow.slope = input.number(field);
                    else if (field == "origin") inflow.origin = input.number(field);
                    else input.fail("unknown interface inflow key '" + field + "'");
                }
                input.expect("}");
                for (const auto* field : {"offset", "slope", "origin"})
                    if (!fields.count(field)) input.fail("interface inflow requires '" + std::string(field) + "'");
                result.inflow.push_back(std::move(inflow));
            }
            input.expect("}");
        }
        else input.fail("unknown interfaceBoundary key '" + key + "'");
    }
    input.expect("}");
    for (const auto* key : {"model", "fixedContacts", "inflow"})
        if (!keys.count(key)) input.fail("interfaceBoundary requires '" + std::string(key) + "'");
    return result;
}

}

bool readNodalInterfaceEntry(const std::string& key, DictionaryInput::Reader& input,
    InterfaceNodalConfig& config)
{
    if (key == "advectInterface") config.advect = input.boolean(key);
    else if(key=="transport") {
        config.transport=input.next();
        if(config.transport!="nodal" && config.transport!="conservative" && config.transport!="parametricGals") input.fail("transport must be nodal, conservative or parametricGals");
    }
    else if(key=="portPhases") {
        input.expect("{");std::set<std::string> names;
        while(input.peek()!="}") {
            const auto name=input.next(),phase=input.next();
            if(!names.insert(name).second || (phase!="liquid" && phase!="gas")) input.fail("portPhases requires unique patch names followed by liquid or gas");
            config.portPhases.push_back({name,phase=="liquid"?-1:1});
        }
        input.expect("}");
    }
    else if (key == "reinitializeEvery") config.reinitializeEvery = input.integer(key);
    else if(key=="neckBreakup") {
        InterfaceNodalConfig::NeckBreakup control;std::set<std::string> seen;
        input.expect("{");
        while(input.peek()!="}") {
            const auto name=input.next();if(!seen.insert(name).second)input.fail("duplicate neckBreakup key");
            if(name=="maximumLevel")control.maximumLevel=input.integer(name);
            else if(name=="radiusCells")control.radiusCells=input.number(name);
            else input.fail("unknown neckBreakup key '"+name+"'");
        }
        input.expect("}");
        if(seen.size()!=2 || control.maximumLevel<0 || control.maximumLevel>24
            || !(control.radiusCells>0 && control.radiusCells<1))
            input.fail("neckBreakup requires maximumLevel in [0,24] and radiusCells in (0,1)");
        config.neckBreakup=control;
    }
    else if (key == "interfaceBoundary") config.boundary = readInterfaceBoundary(input);
    else return false;
    return true;
}

InterfaceNodalConfig readNodalInterfaceCase(const std::string& path)
{
    DictionaryInput::Reader input(path);
    input.expect("Interface");
    input.expect("{");
    InterfaceNodalConfig result;
    std::set<std::string> seen;
    while (input.peek() != "}")
    {
        const auto key = input.next();
        if (!seen.insert(key).second) input.fail("duplicate Interface key '" + key + "'");
        if (!readNodalInterfaceEntry(key, input, result)) input.fail("unknown nodal Interface key '" + key + "'");
    }
    input.expect("}");
    if (!input.done()) input.fail("unexpected content after Interface block");
    if (!seen.count("advectInterface")) input.fail("Interface requires advectInterface");
    if(result.neckBreakup && (!result.advect || result.transport!="conservative"))
        input.fail("neckBreakup requires advectInterface true and transport conservative");
    if (result.reinitializeEvery < 0) input.fail("reinitializeEvery must be nonnegative (0 disables redistancing)");
    if (!result.advect && result.reinitializeEvery != 0) input.fail("reinitializeEvery requires advectInterface true");
    if (seen.count("interfaceBoundary"))
    {
        if (result.boundary.model != "prescribedScalar") input.fail("interfaceBoundary.model must be prescribedScalar");
        if (result.boundary.fixedContacts.empty()) input.fail("interfaceBoundary requires nonempty fixedContacts");
    }
    return result;
}

FittedInterfaceCaseConfig readFittedInterfaceCase(const std::string& path,
    DictionaryInput::ConfigurationTrace& trace)
{
    DictionaryInput::Document document(path);
    document.root().allow({"Interface"});
    return readFittedInterfaceCase(document.root().block("Interface"), trace);
}

FittedInterfaceCaseConfig readFittedInterfaceCase(const DictionaryInput::Node& interface,
    DictionaryInput::ConfigurationTrace& trace)
{
    FittedInterfaceCaseConfig result;
    interface.allow({"levelSet", "movingBoundary", "transport", "ContactLine"});
    const auto& levelSet = interface.get("levelSet");
    result.fluidInterfaceLevelSet = trace.record(levelSet, levelSet.word(), "name", "shared moving level set");
    if (result.fluidInterfaceLevelSet.empty()) levelSet.fail("a nonempty geometry level-set name");
    const auto& moving = interface.get("movingBoundary");
    const int boundaryId = moving.integer();
    if (boundaryId <= 0) moving.fail("a positive 1-based ID");
    result.movingBoundaryId = trace.record(moving, boundaryId, "1-based ID", "moving free-surface boundary; fixed walls stay fixed");
    const auto& transport = interface.get("transport");
    trace.record(transport, transport.choice({"activeHermiteMaterial"}), "model", "existing active-scalar material trajectory and volume correction");
    const auto& contact = interface.block("ContactLine");
    contact.allow({"model", "points"});
    const auto& model = contact.get("model");
    result.contactLineModel = trace.record(model, model.choice({"pinned"}), "model", "fixed contact position; contact angle not prescribed");
    const auto& points = contact.get("points");
    result.pinnedContactPointIds = trace.record(points, points.ids(), "1-based IDs", "pinned geometry point positions");
    return result;
}

InterfaceCaseConfig readInterfaceCase(const std::string& path,
    DictionaryInput::ConfigurationTrace& trace)
{
    DictionaryInput::Document document(path);
    document.root().allow({"Interface"});
    const auto& root=document.root().block("Interface");
    root.allow({"units","excludedPhase","Initial","Transport","Boundaries"});
    trace.record(root.get("units"),root.get("units").choice({"code"}),"","prescribed scalar/velocity units");
    InterfaceCaseConfig result;
    result.excludedPhase=trace.record(root.get("excludedPhase"),root.get("excludedPhase").integer(),"1-based ID","fixed-domain active phase selection");
    if(result.excludedPhase<1)root.get("excludedPhase").fail("positive solid phase ID");
    const auto& initial=root.block("Initial");
    initial.allow({"model","pinnedPoint","height"});
    trace.record(initial.get("model"),initial.get("model").choice({"parabolicCap"}),"","InterfacePrescribedData::initializeParabolicCap");
    result.pinnedPoint=trace.record(initial.get("pinnedPoint"),initial.get("pinnedPoint").integer(),"1-based ID","pinned coordinate and column radius");
    result.initialHeight=trace.record(initial.get("height"),initial.get("height").number(),"code length","initial scalar cap height");
    if(result.pinnedPoint<1)initial.get("pinnedPoint").fail("positive geometry point ID");
    if(result.initialHeight<=0)initial.get("height").fail("positive cap height");
    const auto& transport=root.block("Transport");
    transport.allow({"model","halfWidthCells","Velocity"});
    const auto method=trace.record(transport.get("model"),transport.get("model").choice({"pinnedToAxis","conservative"}),"","Interface motion service selection");
    result.transport=method=="conservative"?InterfaceCaseConfig::Transport::Conservative:InterfaceCaseConfig::Transport::PinnedToAxis;
    result.requiredHalfWidthCells=trace.record(transport.get("halfWidthCells"),transport.get("halfWidthCells").number(),"cells","physical support supplied to complete motion service");
    if(result.requiredHalfWidthCells<=0)transport.get("halfWidthCells").fail("positive support width");
    const auto& velocity=transport.block("Velocity");
    velocity.allow({"model","peakSpeed"});
    trace.record(velocity.get("model"),velocity.get("model").choice({"axialColumn"}),"","InterfacePrescribedData::AxialColumn velocity and exact native flux");
    result.peakSpeed=trace.record(velocity.get("peakSpeed"),velocity.get("peakSpeed").number(),"code velocity","prescribed axial velocity, not a solved flow");
    const auto& boundaries=root.block("Boundaries");
    boundaries.allow({"stationaryWall","symmetry","singlePhasePort"});
    result.stationaryWalls=trace.record(boundaries.get("stationaryWall"),boundaries.get("stationaryWall").ids(),"1-based IDs","fixed interface wall constraints");
    result.symmetry=trace.record(boundaries.get("symmetry"),boundaries.get("symmetry").ids(),"1-based IDs","interface symmetry boundary constraints");
    result.singlePhasePorts=trace.record(boundaries.get("singlePhasePort"),boundaries.get("singlePhasePort").ids(),"1-based IDs","strict single-phase port separation");
    std::set<int> used;
    for(const auto& list:{result.stationaryWalls,result.symmetry,result.singlePhasePorts})
        for(int id:list)if(!used.insert(id).second)boundaries.fail("each boundary ID in exactly one constraint group");
    return result;
}
