#include "Electrospray/IO/ElectrosprayResolvedReport.h"
#include "Electrospray/Compatibility/ElectrosprayInputAdapters.h"
#include "Electrospray/IO/ElectrosprayVersionedInput.h"
#include <algorithm>
#include "Core/Units/UnitsContext.h"
#include "Core/IO/OutputLayout.h"
#include "Core/IO/DictionaryMesh.h"
#include <fstream>
#include <cmath>
#include <stdexcept>
#include <iomanip>
#include <type_traits>

namespace
{
void quoted(std::ostream& out, const std::string& value)
{
    out << '"';
    for (unsigned char c : value)
    {
        if (c == '"' || c == '\\') out << '\\' << c;
        else if (c < 32) out << "\\u00" << "0123456789abcdef"[c >> 4] << "0123456789abcdef"[c & 15];
        else out << c;
    }
    out << '"';
}
void jsonValue(std::ostream& out, const DictionaryInput::ConfigurationValue& value)
{
    std::visit([&](const auto& item) {
        using T = std::decay_t<decltype(item)>;
        if constexpr (std::is_same_v<T, std::string>) quoted(out, item);
        else if constexpr (std::is_same_v<T, bool>) out << (item ? "true" : "false");
        else if constexpr (std::is_same_v<T, std::vector<int>>)
        { out << '['; bool first = true; for (int id : item) { if (!first) out << ','; first = false; out << id; } out << ']'; }
        else
        {
            if constexpr (std::is_same_v<T, double>)
                if (!std::isfinite(item)) throw std::runtime_error("non-finite resolved configuration value");
            out << item;
        }
    }, value);
}
}

static std::map<std::string, DictionaryInput::ConfigurationValue> resolvedValues(
    const ElectrosprayModuleInputs& input, const UnitsContext& units, ElectrosprayProblem problem,
    int initialAdaptIterations = -1, bool initializeLevelSets = true, int restartLoop = -1)
{
    const auto& flow = input.flow;
    const auto& interfaceConfig = input.interfaceConfig;
    const auto& electric = input.electric;
    const auto& run = input.controls;
    const auto& liquid = input.materials.phases.at(flow.liquidPhase);
    const auto& gas = input.materials.phases.at(flow.gasPhase);
    std::map<std::string, DictionaryInput::ConfigurationValue> values;
    const bool charged = problem == ElectrosprayProblem::Conducting;
    values["workflow"] = std::string(electrosprayProblemLabel(problem));
    values["electric.model.actual"] = std::string(charged ? "leakyDielectric" : "perfectDielectric");
    values["electric.interfaceDiscretization"] = electric.interfaceDiscretization;
    values["electric.initialSurfaceCharge_C_m2"] = electric.initialSurfaceChargeSI;
    values["electric.conductorPhase"] = electric.conductorPhaseId;
    values["electric.conductorPotential_V"] = electric.conductorPotentialSI;
    for (const auto* phase : {&liquid, &gas})
    {
        const auto prefix = "materials.phase" + std::to_string(phase->phaseId) + '.';
        values[prefix + "density_kg_m3"] = phase->densitySI;
        values[prefix + "dynamicViscosity_Pa_s"] = phase->dynamicViscositySI;
        values[prefix + "relativePermittivity"] = phase->relativePermittivity;
        if (charged) values[prefix + "conductivity_S_m"] = phase->conductivitySI;
    }
    values["materials.surfaceTension_N_m"] = input.materials.surfaceTensionSI.at(std::minmax(flow.liquidPhase, flow.gasPhase));
    values["hydro.liquidPhase"] = liquid.phaseId; values["hydro.gasPhase"] = gas.phaseId;
    values["hydro.supply.id"] = flow.supply.boundaryId; values["hydro.supply.pressure_Pa"] = flow.supply.pressureSI;
    values["hydro.outlet.id"] = flow.outlet.boundaryId; values["hydro.outlet.pressure_Pa"] = flow.outlet.pressureSI;
    values["hydro.outlet.incomingPhase"] = flow.outletIncomingPhaseId;
    for (const auto& bc : flow.ambient) values["hydro.ambient" + std::to_string(bc.boundaryId) + ".pressure_Pa"] = bc.pressureSI;
    for (const auto& bc : electric.domainBoundaries)
    {
        const auto prefix = "electric.boundary" + std::to_string(bc.boundaryId) + '.';
        values[prefix + "type"] = bc.type;
        if (bc.type == "dirichlet") values[prefix + "potential_V"] = bc.potentialSI;
    }
    values["hydro.incomingMomentumModel"] = flow.incomingMomentumModel;
    values["hydro.solidWallModel"] = flow.solidWallModel;
    values["interface.levelSet"] = interfaceConfig.fluidInterfaceLevelSet; values["interface.movingBoundary"] = interfaceConfig.movingBoundaryId;
    values["interface.contactLineModel"] = interfaceConfig.contactLineModel; values["interface.pinnedPoints"] = interfaceConfig.pinnedContactPointIds;
    values["interface.forceHalfWidthCells"] = flow.capillaryHalfWidthCells;
    values["numerics.momentumLinear"] = input.systems.at("Momentum").relativeTolerance; values["numerics.pressureLinear"] = input.systems.at("Pressure").relativeTolerance;
    values["numerics.coupledLinear"] = input.systems.at("CoupledVelocityPressure").relativeTolerance; values["numerics.momentumBalance"] = flow.momentumBalance;
    values["numerics.localContinuity"] = flow.localContinuity;
    if (flow.outletIncomingPhaseId == -1) values["numerics.outletAdmissibility"] = flow.outletAdmissibility;
    values["numerics.electricPotentialLinear"] = input.systems.at("ElectricPotential").relativeTolerance; values["numerics.potentialBounds"] = electric.potentialBounds;
    values["numerics.volumeCorrectionLinear"] = input.systems.at("InterfaceVolumeCorrection").relativeTolerance; values["numerics.phaseFluxBalance"] = flow.phaseFluxBalance;
    if (charged)
    {
        values["numerics.reconstructedPoisson"] = electric.fluxChecks.reconstructedPoisson;
        values["numerics.gaussClosure"] = electric.fluxChecks.gaussClosure; values["numerics.chargeClosure"] = electric.chargeClosure;
    }
    else values["numerics.polarizationConsistency"] = electric.polarizationConsistency;
    values["linearBackend"] = std::string("SuperLU");
    values["time.dt_code"] = run.dt; values["time.endTime_code"] = run.endTime;
    values["output.prefix"] = run.outputPrefix; values["output.fields"] = run.writeFields;
    values["output.diagnosticFields"] = run.diagnosticFields;
    values["output.initial"] = run.writeInitial; values["output.first"] = run.writeFirst;
    values["output.final"] = run.writeFinal; values["output.everySteps"] = run.writeEverySteps;
    if (!run.writeEverySteps) values["output.interval_code"] = run.writeInterval;
    values["output.detailedFields"] = run.detailedLog;
    values["initial.refinementIterationsOverride"] = initialAdaptIterations;
    values["initial.levelSets"] = initializeLevelSets; values["initial.restartLoop"] = restartLoop;
    values["units.L0_m"] = units.scales().L0_m; values["units.T0_s"] = units.scales().T0_s;
    values["units.Phi0_V"] = units.scales().Phi0_V; values["units.rho0_C_m3"] = units.scales().rho0_Cpm3;
    values["units.v0_m_s"] = units.scales().v0_mps; values["units.E0_V_m"] = units.scales().E0_Vpm;
    return values;
}

static void writeResolvedReport(const std::map<std::string, DictionaryInput::ConfigurationValue>& values,
    const std::string& requestedModel, const std::string& legacyToleranceSource,
    int logEverySteps, bool detailedIterationLog, int schemaVersion,
    const DictionaryInput::ConfigurationTrace* trace, const DictionaryMesh* mesh)
{
    const auto path = OutputLayout::meta() / (std::get<std::string>(values.at("output.prefix")) + "_resolved.json");
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write resolved configuration: " + path.string());
    out << std::setprecision(17) << "{\n  \"schemaVersion\": " << schemaVersion
        << ",\n  \"executionStatus\": \"pre-initialization configuration snapshot; not an accepted step\",\n  \"requestedElectricModel\": ";
    quoted(out, requestedModel);
    out << ",\n  \"legacyToleranceSource\": "; quoted(out, legacyToleranceSource);
    out << ",\n  \"logEverySteps\": " << logEverySteps << ",\n  \"detailedIterationLog\": " << (detailedIterationLog ? "true" : "false")
        << ",\n  \"normalized\": {";
    bool first = true;
    for (const auto& item : values)
    {
        if (!first) out << ','; first = false;
        out << "\n    "; quoted(out, item.first); out << ": "; jsonValue(out, item.second);
    }
    out << "\n  },";
    if (mesh && mesh->runtimeAdaptation.present)
    {
        const auto& amr = mesh->runtimeAdaptation;
        out << "\n  \"runtimeAdaptation\": {\"source\":\"system/mesh.dict/RuntimeAdaptation\","
            "\"consumer\":\"ElectrosprayTimeLoop / QuadtreeMeshRefiner / QuadtreeMesh::adapt\","
            "\"enabled\":" << (amr.enabled ? "true" : "false");
        if (amr.enabled) out << ",\"everySteps\":" << amr.everySteps << ",\"maxSweeps\":" << amr.maxSweeps;
        out << "},";
    }
    if (mesh && mesh->physicsRefinement.enabled)
    {
        const auto& flow = mesh->physicsRefinement;
        out << "\n  \"flowRefinement\": {\"source\":\"system/mesh.dict/PhysicsRefinement\","
            "\"consumer\":\"HydroFittedRefinement / QuadtreeMeshRefiner\","
            "\"firstEvaluation\":\"after first accepted step\",\"hysteresis\":" << flow.hysteresis
            << ",\"indicators\":[";
        for (std::size_t i = 0; i < flow.indicators.size(); ++i)
        {
            if (i) out << ',';
            const bool speed = flow.indicators[i] == "Velocity";
            out << "{\"name\":"; quoted(out, flow.indicators[i]);
            out << ",\"thresholdSI\":" << (speed ? flow.velocityThresholdSI : flow.velocityGradientThresholdSI)
                << ",\"unit\":\"" << (speed ? "m/s" : "1/s") << "\",\"maxLevel\":"
                << flow.maxLevelCaps.at(speed ? "velocity" : "velocitygradient") << '}';
        }
        out << "]},";
    }
    out << "\n  \"sources\": [";
    first = true;
    if (trace) for (const auto& entry : trace->entries)
    {
        if (!first) out << ','; first = false;
        out << "\n    {\"field\":"; quoted(out, entry.field); out << ",\"file\":"; quoted(out, entry.file);
        out << ",\"line\":" << entry.line << ",\"raw\":"; quoted(out, entry.rawValue);
        out << ",\"value\":"; jsonValue(out, entry.value); out << ",\"units\":"; quoted(out, entry.units);
        out << ",\"defaulted\":" << (entry.defaulted ? "true" : "false") << ",\"used\":" << (entry.used ? "true" : "false");
        out << ",\"consumer\":"; quoted(out, entry.consumer); out << ",\"unusedReason\":"; quoted(out, entry.unusedReason); out << '}';
    }
    out << "\n  ]\n}\n";
    out.close();
    if (!out) throw std::runtime_error("failed to complete resolved configuration: " + path.string());
}

std::map<std::string, DictionaryInput::ConfigurationValue> electrosprayResolvedValues(
    const ElectrosprayVersionedInput& input)
{
    return resolvedValues(input.modules, *input.units, input.problem);
}

void writeElectrosprayResolvedReport(const ElectrosprayVersionedInput& input, const DictionaryMesh* mesh)
{
    writeResolvedReport(electrosprayResolvedValues(input), input.modules.electric.model, "",
        input.modules.controls.logEverySteps, input.modules.controls.detailedLog, 2, &input.trace, mesh);
}

// Compatibility overloads retain the old requested model/mode and independent
// initialization/log flags; they share the same owner-based value serialization.
std::map<std::string, DictionaryInput::ConfigurationValue> electrosprayResolvedValues(
    const ElectrosprayPhysicalConfig& physical, const ElectrosprayConfig& run, const UnitsContext& units)
{
    const auto workflow = adaptLegacyElectrosprayWorkflow(run);
    return resolvedValues(adaptLegacyElectrosprayInputs(physical, run, workflow), units, workflow.problem,
        run.initialAdaptIterations, run.initializeLevelSets, run.restartLoop);
}

void writeElectrosprayResolvedReport(const ElectrosprayPhysicalConfig& physical, const ElectrosprayConfig& run,
    const UnitsContext& units, int schemaVersion, const DictionaryInput::ConfigurationTrace* trace,
    const DictionaryMesh* mesh)
{
    writeResolvedReport(electrosprayResolvedValues(physical, run, units), physical.model, physical.legacyToleranceSource,
        run.logEvery, run.detailedIterationLog, schemaVersion, trace, mesh);
}
