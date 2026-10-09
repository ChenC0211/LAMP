#include "Modules/Hydrodynamics/IO/DictionaryHydrodynamics.h"
#include "Modules/Hydrodynamics/IO/HydroInputParser.h"
#include "Modules/Interface/Compatibility/IO/DictionaryInterface.h"

#include <algorithm>
#include <fstream>
#include <limits>
#include <set>
#include <utility>

namespace
{
std::vector<HydroBoundaryPatch> readBoundaryPatches(HydroInput::Reader& input)
{
    std::vector<HydroBoundaryPatch> patches;
    input.expect("{");
    std::set<std::string> names;
    while (input.peek() != "}")
    {
        HydroBoundaryPatch patch;
        patch.name = input.next();
        if (!names.insert(patch.name).second) input.fail("duplicate boundary patch '" + patch.name + "'");
        input.expect("{");
        std::set<std::string> keys;
        while (input.peek() != "}")
        {
            const auto key = input.next();
            if (!keys.insert(key).second) input.fail("duplicate boundary patch key '" + key + "'");
            if (key == "side")
            {
                const auto side = input.next();
                if (side == "west") patch.side = HydroBoundarySide::West;
                else if (side == "south") patch.side = HydroBoundarySide::South;
                else if (side == "east") patch.side = HydroBoundarySide::East;
                else if (side == "north") patch.side = HydroBoundarySide::North;
                else input.fail("boundary side must be west, south, east or north");
            }
            else if (key == "type")
            {
                const auto type = input.next();
                if (type == "normalVelocity") patch.type = HydroBoundaryType::NormalVelocity;
                else if (type == "pressureOutlet") patch.type = HydroBoundaryType::PressureOutlet;
                else if (type == "freeSlip") patch.type = HydroBoundaryType::FreeSlip;
                else input.fail("boundary type must be normalVelocity, pressureOutlet or freeSlip");
            }
            else if (key == "start") patch.start = input.number(key);
            else if (key == "end") patch.end = input.number(key);
            else if (key == "value") patch.value = input.number(key);
            else if (key == "incomingVelocity") {
                input.expect("(");const double x=input.number(key);if(input.peek()==",")input.next();const double y=input.number(key);input.expect(")");patch.incomingVelocity=Vector2D{x,y};
            }
            else input.fail("unknown boundary patch key '" + key + "'");
        }
        input.expect("}");
        for (const auto* key : {"side", "type", "start", "end"})
            if (!keys.count(key)) input.fail("boundary patch '" + patch.name + "' requires '" + key + "'");
        if (patch.type != HydroBoundaryType::FreeSlip && !keys.count("value"))
            input.fail("boundary patch '" + patch.name + "' requires 'value'");
        if (patch.type == HydroBoundaryType::FreeSlip && (patch.value != 0 || patch.incomingVelocity))
            input.fail("freeSlip requires zero normal velocity and no incomingVelocity");
        patches.push_back(std::move(patch));
    }
    input.expect("}");
    if (patches.empty()) input.fail("boundaryPatches cannot be empty; omit the block for a closed domain");
    return patches;
}
HydrodynamicsConfig readHydrodynamics(const std::string& path,
    DictionaryHydrodynamics::RunInput* run, InterfaceNodalConfig* legacyInterface,
    const InterfaceNodalConfig& interface_)
{
    HydrodynamicsConfig cfg_;
    HydroInput::Reader input(path);
    input.expect("Hydrodynamics");
    input.expect("{");
    std::set<std::string> seen;
    std::set<std::string> conservationKeys;
    while (input.peek() != "}")
    {
        const std::string key = input.next();
        if (!seen.insert(key).second) input.fail("duplicate Hydrodynamics key '" + key + "'");
        if (!run && (key == "steps" || key == "dt" || key == "outputName" ||
            key == "writeFields" || key == "writeInitialFields" || key == "writeEvery" || key == "logEvery"))
            input.fail("run settings belong in run.dict/Time, Output and Log; mixed sources are not supported");
        if (!legacyInterface && (key == "advectInterface" || key == "reinitializeEvery" || key == "interfaceBoundary"))
            input.fail("interface settings belong in interface.dict; mixed sources are not supported");
        if (key == "steps") run->steps = input.integer(key);
        else if (key == "dt") run->dt = input.number(key);
        else if (key == "coordinateSystem") cfg_.coordinateSystem = input.next();
        else if (key == "timeScheme") cfg_.timeScheme = input.next();
        else if (key == "Moment") {
            input.expect("{");std::set<std::string> keys;
            while(input.peek()!="}"){const auto setting=input.next();
                if(!keys.insert(setting).second)input.fail("duplicate Moment key '"+setting+"'");
                if(setting=="maximumIterations")cfg_.momentMaximumIterations=input.integer(setting);
                else if(setting=="convection")cfg_.momentConvection=input.next();
                else if(setting=="neckRadiusCells")cfg_.momentNeckRadiusCells=input.number(setting);
                else input.fail("unknown Moment key '"+setting+"'");}
            input.expect("}");
        }
        else if (key == "Midpoint") {
            input.expect("{");std::set<std::string> keys;
            while(input.peek()!="}") {const auto setting=input.next();
                if(!keys.insert(setting).second)input.fail("duplicate Midpoint key '"+setting+"'");
                if(setting=="maximumIterations")cfg_.midpoint.maximumIterations=input.integer(setting);
                else if(setting=="relativeVelocityTolerance")cfg_.midpoint.relativeVelocityTolerance=input.number(setting);
                else if(setting=="absoluteVelocityTolerance")cfg_.midpoint.absoluteVelocityTolerance=input.number(setting);
                else input.fail("unknown Midpoint key '"+setting+"'");}
            input.expect("}");
        }
        else if (key == "density") cfg_.density = input.number(key);
        else if (key == "viscosity") cfg_.viscosity = input.number(key);
        else if (key == "materialModel") cfg_.materialModel = input.next();
        else if (key == "interfaceLevelSet") cfg_.interfaceLevelSet = input.next();
        else if (key == "densityNegative") cfg_.densityNegative = input.number(key);
        else if (key == "densityPositive") cfg_.densityPositive = input.number(key);
        else if (key == "viscosityNegative") cfg_.viscosityNegative = input.number(key);
        else if (key == "viscosityPositive") cfg_.viscosityPositive = input.number(key);
        else if (key == "interfaceThickness") cfg_.interfaceThickness = input.number(key);
        else if (legacyInterface && readNodalInterfaceEntry(key, input, *legacyInterface)) {}
        else if (key == "surfaceTension") cfg_.surfaceTension = input.number(key);
        else if (key == "interfaceForceGeometry") cfg_.interfaceForceGeometry = input.next();
        else if (key == "interfaceForceLayer") cfg_.interfaceForceLayer = input.next();
        else if (key == "interfaceNormalTraction") cfg_.interfaceNormalTraction = input.number(key);
        else if (key == "interfaceTangentialTraction") cfg_.interfaceTangentialTraction = input.number(key);
        else if (key == "accelerationX") cfg_.accelerationX = input.number(key);
        else if (key == "accelerationY") cfg_.accelerationY = input.number(key);
        else if (key == "boundaryPatches") cfg_.boundaryPatches = readBoundaryPatches(input);
        else if (key == "initialVelocity") cfg_.initialVelocity = input.next();
        else if (key == "initialPressureJump") cfg_.initialPressureJump = input.number(key);
        else if (key == "velocityAmplitude") cfg_.velocityAmplitude = input.number(key);
        else if (key == "linearResidualTolerance") cfg_.linearResidualTolerance = input.number(key);
        else if (key == "Conservation")
        {
            if (run) input.fail("Conservation requires separate Hydro inputs with Time/Output/Log in run.dict");
            input.expect("{");
            while (input.peek() != "}")
            {
                const auto setting = input.next();
                if (!conservationKeys.insert(setting).second) input.fail("duplicate Conservation key '" + setting + "'");
                if (setting == "momentumBalance") cfg_.conservation.momentumBalance = input.number(setting);
                else if (setting == "initialContinuity") cfg_.conservation.initialContinuity = input.number(setting);
                else if (setting == "outletAdmissibility") cfg_.conservation.outletAdmissibility = input.number(setting);
                else if (setting == "boundaryFluxBalance") cfg_.conservation.boundaryFluxBalance = input.number(setting);
                else input.fail("unknown Conservation key '" + setting + "'");
            }
            input.expect("}");
        }
        else if (key == "outputName") run->outputName = input.next();
        else if (key == "writeFields") run->writeFields = input.boolean(key);
        else if (key == "writeInitialFields") run->writeInitialFields = input.boolean(key);
        else if (key == "writeEvery") run->writeEvery = input.integer(key);
        else if (key == "logEvery") run->logEvery = input.integer(key);
        else input.fail("unknown or unsupported Hydrodynamics key '" + key +
                        "'; supported scope is fixed-grid flow with optional two-material level-set transport and interface traction");
    }
    input.expect("}");
    if (!input.done()) input.fail("unexpected content after Hydrodynamics block");
    for (const auto* key : {"steps", "dt", "coordinateSystem"})
        if ((run || std::string(key) == "coordinateSystem") && !seen.count(key)) input.fail(std::string("missing required Hydrodynamics key '") + key + "'");
    if (run && (run->steps < 1 || run->steps == std::numeric_limits<int>::max()))
        input.fail("steps must be a positive integer below INT_MAX");
    if (run && (!(run->dt > 0.0) || !std::isfinite(run->dt * run->steps))) input.fail("dt must be positive with a finite end time");
    if (cfg_.coordinateSystem != "planar2D" && cfg_.coordinateSystem != "axisymmetricRZ")
        input.fail("coordinateSystem must be planar2D or axisymmetricRZ");
    const bool conservative=interface_.transport=="conservative";
    const bool parametric=interface_.transport=="parametricGals";
    if(cfg_.momentNeckRadiusCells!=0&&(!parametric||cfg_.coordinateSystem!="axisymmetricRZ"||!(cfg_.momentNeckRadiusCells>0&&cfg_.momentNeckRadiusCells<=.5)))
        input.fail("Moment neckRadiusCells requires parametric RZ and a value in (0,0.5]; omit or use zero to disable");
    if(parametric && (cfg_.timeScheme!="cchMomentMidpoint" || cfg_.materialModel!="twoPhase"
        || !interface_.advect || interface_.reinitializeEvery!=0 || !interface_.boundary.model.empty() || interface_.neckBreakup
        || cfg_.interfaceForceGeometry!="parametric" || cfg_.interfaceForceLayer!="sharp" || cfg_.momentConvection!="central"
        || cfg_.interfaceNormalTraction!=0 || cfg_.interfaceTangentialTraction!=0))
        input.fail("parametricGals requires XY/RZ cchMomentMidpoint twoPhase, parametric/sharp traction, central convection and no redistance/moving contacts/prescribed traction");
    for(const auto& patch:cfg_.boundaryPatches)if(patch.incomingVelocity&&cfg_.timeScheme!="cchMomentMidpoint")input.fail("incomingVelocity requires cchMomentMidpoint");
    for(const auto& patch:cfg_.boundaryPatches)if(patch.type==HydroBoundaryType::FreeSlip&&cfg_.timeScheme!="cchMomentMidpoint")input.fail("freeSlip requires cchMomentMidpoint");
    if(cfg_.timeScheme!="cnAb2" && cfg_.timeScheme!="cchMidpoint" && cfg_.timeScheme!="cchMomentMidpoint")input.fail("timeScheme must be cnAb2, cchMidpoint or cchMomentMidpoint");
    if(seen.count("Moment") && cfg_.timeScheme!="cchMomentMidpoint")input.fail("Moment controls require timeScheme cchMomentMidpoint");
    if(cfg_.timeScheme=="cchMomentMidpoint") {
        if(run || (!parametric && cfg_.materialModel!="constant")
            || conservative || (!parametric&&interface_.advect) || cfg_.initialVelocity=="smooth" || seen.count("Conservation"))
            input.fail("cchMomentMidpoint requires separate inputs, XY/RZ constant material or parametric interface, zero/vortex initial velocity and no legacy Interface/Conservation controls");
        if(cfg_.momentMaximumIterations<1 || cfg_.momentMaximumIterations>100 || cfg_.momentConvection!="central")
            input.fail("Moment requires 1..100 iterations and convection central; traceResidual belongs to the retired low-order time update");
    }
    if(seen.count("Midpoint") && cfg_.timeScheme!="cchMidpoint")input.fail("Midpoint controls require timeScheme cchMidpoint");
    if(cfg_.timeScheme=="cchMidpoint") {
        if(cfg_.materialModel!="constant" || !cfg_.boundaryPatches.empty() || conservative)
            input.fail("cchMidpoint case currently requires constant single-phase material and closed boundaries");
        if(cfg_.midpoint.maximumIterations<1 || cfg_.midpoint.maximumIterations>100
            || !(cfg_.midpoint.relativeVelocityTolerance>0) || !(cfg_.midpoint.absoluteVelocityTolerance>0))
            input.fail("Midpoint requires 1..100 iterations and positive velocity tolerances");
    }
    if (seen.count("initialPressureJump") && (!conservative || !cfg_.boundaryPatches.empty()))
        input.fail("initialPressureJump requires conservative interface transport and a closed domain without boundaryPatches");
    if (cfg_.interfaceForceGeometry != "legacyContour" && cfg_.interfaceForceGeometry != "hermite" && cfg_.interfaceForceGeometry!="conservative" && !(parametric&&cfg_.interfaceForceGeometry=="parametric"))
        input.fail("interfaceForceGeometry must be legacyContour, hermite or conservative");
    if (cfg_.interfaceForceLayer != "geometricPL" && cfg_.interfaceForceLayer != "regularized" && cfg_.interfaceForceLayer!="inventory" && !(parametric&&cfg_.interfaceForceLayer=="sharp"))
        input.fail("interfaceForceLayer must be geometricPL, regularized or inventory");
    if(conservative) {
        if(cfg_.materialModel!="twoPhase"
            || cfg_.interfaceForceGeometry!="conservative" || cfg_.interfaceForceLayer!="inventory"
            || !interface_.advect || interface_.reinitializeEvery!=0 || !interface_.boundary.model.empty()
            || cfg_.interfaceNormalTraction!=0 || cfg_.interfaceTangentialTraction!=0)
            input.fail("conservative transport requires twoPhase, conservative/inventory capillarity, advection, and no redistance/contact/traction");
    }
    if(conservative||parametric){
        for(const auto& phase:interface_.portPhases) {
            const auto patch=std::find_if(cfg_.boundaryPatches.begin(),cfg_.boundaryPatches.end(),[&](const auto& p){return p.name==phase.patch;});
            if(patch==cfg_.boundaryPatches.end() || patch->type==HydroBoundaryType::FreeSlip || (patch->type==HydroBoundaryType::NormalVelocity && patch->value==0))
                input.fail("portPhases must name actual nonzero normal-velocity or pressure-outlet patches");
            if(parametric)patch->incomingPhase=phase.sign;
        }
        for(const auto& patch:cfg_.boundaryPatches) if(patch.type==HydroBoundaryType::PressureOutlet || patch.value!=0)
            if(std::none_of(interface_.portPhases.begin(),interface_.portPhases.end(),[&](const auto& p){return p.patch==patch.name;}))
                input.fail("every two-phase inlet/outlet requires an explicit portPhases declaration");
        if(parametric)for(const auto& patch:cfg_.boundaryPatches)if(patch.type!=HydroBoundaryType::FreeSlip&&patch.incomingPhase==0)input.fail("every parametric inlet/outlet requires an exterior portPhases declaration");
    } else if(cfg_.interfaceForceGeometry=="conservative" || cfg_.interfaceForceLayer=="inventory" || !interface_.portPhases.empty())
        input.fail("conservative force geometry, inventory layer and portPhases require conservative interface transport");
    if (cfg_.interfaceForceLayer == "regularized" && cfg_.interfaceForceGeometry != "hermite")
        input.fail("interfaceForceLayer regularized requires interfaceForceGeometry hermite");
    if (cfg_.interfaceForceGeometry == "hermite")
    {
        if (cfg_.materialModel != "twoPhase")
            input.fail("interfaceForceGeometry hermite requires materialModel twoPhase");
        if (cfg_.coordinateSystem == "axisymmetricRZ" && cfg_.interfaceForceLayer != "regularized")
            input.fail("axisymmetricRZ Hermite loads currently require interfaceForceLayer regularized");
        if (cfg_.interfaceNormalTraction != 0.0 || cfg_.interfaceTangentialTraction != 0.0)
            input.fail("interfaceForceGeometry hermite does not support prescribed normal/tangential traction");
    }
    if (cfg_.materialModel == "constant")
    {
        for (const auto* key : {"density", "viscosity"})
            if (!seen.count(key)) input.fail(std::string("missing required Hydrodynamics key '") + key + "'");
        if (!(cfg_.density > 0.0)) input.fail("density must be positive");
        if (cfg_.viscosity < 0.0) input.fail("viscosity must be nonnegative");
        for (const auto* key : {"interfaceLevelSet", "densityNegative", "densityPositive",
                "viscosityNegative", "viscosityPositive", "interfaceThickness",
                "advectInterface", "reinitializeEvery", "surfaceTension",
                "interfaceNormalTraction", "interfaceTangentialTraction"})
            if (seen.count(key)) input.fail(std::string(key) + " requires materialModel twoPhase");
    }
    else if (cfg_.materialModel == "twoPhase")
    {
        for (const auto* key : {"interfaceLevelSet", "densityNegative", "densityPositive",
                "viscosityNegative", "viscosityPositive", "interfaceThickness"})
            if (!seen.count(key) && !(parametric&&std::string(key)=="interfaceThickness")) input.fail(std::string("materialModel twoPhase requires '") + key + "'");
        if (seen.count("density") || seen.count("viscosity"))
            input.fail("twoPhase uses densityNegative/Positive and viscosityNegative/Positive, not density/viscosity");
        if (cfg_.interfaceLevelSet.empty()) input.fail("interfaceLevelSet must be a nonempty geometry LevelSet name");
        if (!(cfg_.densityNegative > 0.0) || !(cfg_.densityPositive > 0.0))
            input.fail("densityNegative and densityPositive must be positive");
        if (cfg_.viscosityNegative < 0.0 || cfg_.viscosityPositive < 0.0)
            input.fail("viscosityNegative and viscosityPositive must be nonnegative");
        if(parametric&&seen.count("interfaceThickness"))input.fail("sharp parametric geometry has no interfaceThickness");
        if (!parametric&&!(cfg_.interfaceThickness > 0.0)) input.fail("interfaceThickness must be positive");
        if (cfg_.surfaceTension < 0.0) input.fail("surfaceTension must be nonnegative");
        if (interface_.reinitializeEvery < 0) input.fail("reinitializeEvery must be nonnegative (0 disables redistancing)");
        if (!interface_.advect && interface_.reinitializeEvery != 0)
            input.fail("reinitializeEvery requires advectInterface true");
    }
    else input.fail("materialModel must be constant or twoPhase");
    const bool hasInterfaceBoundary = legacyInterface ? seen.count("interfaceBoundary") != 0 : !interface_.boundary.model.empty();
    if (hasInterfaceBoundary)
    {
        const auto& boundary = interface_.boundary;
        if (boundary.model != "prescribedScalar") input.fail("interfaceBoundary.model must be prescribedScalar");
        if (cfg_.materialModel != "twoPhase" || cfg_.interfaceForceGeometry != "hermite"
            || cfg_.interfaceForceLayer != "regularized" || interface_.reinitializeEvery != 0
            || cfg_.boundaryPatches.empty() || boundary.fixedContacts.empty())
            input.fail("interfaceBoundary requires twoPhase, hermite+regularized, reinitializeEvery 0, explicit boundaryPatches and nonempty fixedContacts");
        for (const auto& inflow : boundary.inflow)
        {
            const auto patch = std::find_if(cfg_.boundaryPatches.begin(), cfg_.boundaryPatches.end(),
                [&](const HydroBoundaryPatch& candidate) { return candidate.name == inflow.patch; });
            if (patch == cfg_.boundaryPatches.end() || patch->type != HydroBoundaryType::NormalVelocity || !(patch->value < 0.0))
                input.fail("interface inflow '" + inflow.patch + "' must name an actual negative normalVelocity patch");
            for (const double s : {patch->start, patch->end})
                if (!std::isfinite(inflow.valueAt(s)))
                    input.fail("interface inflow scalar is unrepresentable on its patch");
        }
        bool hasInflow = false;
        for (const auto& patch : cfg_.boundaryPatches)
            if (patch.type == HydroBoundaryType::NormalVelocity && patch.value < 0.0)
            {
                hasInflow = true;
                if (std::none_of(boundary.inflow.begin(), boundary.inflow.end(),
                    [&](const InterfaceNodalConfig::Inflow& value) { return value.patch == patch.name; }))
                    input.fail("every negative normalVelocity patch requires explicit interface inflow data: '" + patch.name + "'");
            }
        if (!hasInflow) input.fail("interfaceBoundary requires at least one actual negative normalVelocity inlet");
    }
    if (!cfg_.boundaryPatches.empty() && cfg_.materialModel != "constant" && !hasInterfaceBoundary && !conservative&&!parametric)
        input.fail("boundaryPatches currently require constant material; interface inflow/contact and charge boundary transport are not implemented");
    if (run && (run->writeEvery < 1 || run->logEvery < 1)) input.fail("writeEvery and logEvery must be positive integers");
    if (cfg_.initialVelocity != "zero" && cfg_.initialVelocity != "smooth" && cfg_.initialVelocity!="vortex")
        input.fail("initialVelocity must be zero, smooth or vortex");
    if (cfg_.velocityAmplitude < 0.0) input.fail("velocityAmplitude must be nonnegative");
    if (!(cfg_.linearResidualTolerance > 0.0) || cfg_.linearResidualTolerance > 1e-2)
        input.fail("linearResidualTolerance must be in (0, 0.01]");
    if (run) cfg_.conservation = HydroCompleteConservation(cfg_.linearResidualTolerance);
    else
    {
        try { cfg_.conservation.validate(); }
        catch (const std::invalid_argument& error) { input.fail(error.what()); }
        if (cfg_.boundaryPatches.empty())
            for (const auto* key : {"initialContinuity", "outletAdmissibility", "boundaryFluxBalance"})
                if (conservationKeys.count(key)) input.fail(std::string(key) + " requires boundaryPatches");
        if (conservationKeys.count("outletAdmissibility") && std::none_of(cfg_.boundaryPatches.begin(),
            cfg_.boundaryPatches.end(), [](const HydroBoundaryPatch& patch) { return patch.type == HydroBoundaryType::PressureOutlet; }))
            input.fail("outletAdmissibility requires a pressureOutlet patch");
    }
    if (run && (run->outputName.empty() || run->outputName.find_first_not_of(
            "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-") != std::string::npos))
        input.fail("outputName must be a nonempty filename stem containing letters, digits, '_' or '-'");
    return cfg_;
}
}

bool DictionaryHydrodynamics::readDict(const std::string& path)
{
    if (!std::ifstream(path)) return false;
    run_ = RunInput{};
    interface_ = InterfaceNodalConfig{};
    cfg_ = readHydrodynamics(path, &run_, &interface_, interface_);
    return true;
}

HydrodynamicsConfig readNodalHydrodynamicsConfig(const std::string& path,
    const InterfaceNodalConfig& interfaceConfig)
{
    return readHydrodynamics(path, nullptr, nullptr, interfaceConfig);
}
