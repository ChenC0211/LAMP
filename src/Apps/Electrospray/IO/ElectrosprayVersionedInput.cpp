#include "Electrospray/IO/ElectrosprayVersionedInput.h"
#include "Runtime/IO/DictionaryRunControls.h"
#include "Runtime/IO/DictionaryManualUnits.h"
#include "Core/IO/DictionaryMaterials.h"
#include "Core/Units/NormalizationContext.h"
#include "Core/Units/UnitsContext.h"
#include "Modules/Hydrodynamics/IO/DictionaryFittedFlow.h"
#include "Modules/Interface/Compatibility/IO/DictionaryInterface.h"
#include "Modules/Electricity/IO/DictionaryFittedElectric.h"
#include "Solvers/IO/DictionaryLinearSystems.h"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <set>
#include <utility>

ElectrosprayVersionedInput readElectrosprayVersionedInput(const std::string& directory)
{
    const std::set<std::string> files{"run.dict", "geometry.dict", "mesh.dict", "materials.dict",
        "hydro.dict", "interface.dict", "electricity.dict", "solvers.dict", "units.dict", "normalization.dict"};
    for (const auto& entry : std::filesystem::directory_iterator(directory))
        if (!entry.is_regular_file() || !files.count(entry.path().filename().string()))
            throw std::invalid_argument("unsupported schema 2 Electrospray input: " + entry.path().string());
    ElectrosprayVersionedInput input;
    auto& trace = input.trace;
    DictionaryInput::Document document(directory + "/run.dict");
    const auto& root = document.root();
    root.allow({"Run", "Time", "Output", "Log", "Restart", "Electrospray"});
    const auto& run = root.block("Run");
    run.allow({"schemaVersion", "app"});
    if (trace.record(run.get("schemaVersion"), run.get("schemaVersion").integer(), "", "versioned input router") != 2)
        run.get("schemaVersion").fail("2");
    trace.record(run.get("app"), run.get("app").choice({"Electrospray"}), "", "selected workflow");
    input.units = readManualSIUnits(directory + "/units.dict", trace);
    DictionaryInput::Document normalization(directory + "/normalization.dict");
    normalization.root().allow({"Normalization"});
    const auto& norm = normalization.root().block("Normalization");
    norm.allow({"dimension", "strategy"});
    NormalizationConfig normConfig;
    normConfig.dimension = trace.record(norm.get("dimension"), norm.get("dimension").integer(), "", "RZ full-volume closure metadata");
    if (normConfig.dimension != 3) norm.get("dimension").fail("3 for full RZ revolution volumes (not a Cartesian 3D mesh)");
    trace.record(norm.get("strategy"), norm.get("strategy").choice({"manual"}), "", "retain explicit SI references");
    input.normalization = std::make_shared<NormalizationContext>(NormalizationContext::MakeFromUnits(input.units, normConfig));
    auto& modules = input.modules;
    modules.controls = readRunControls(root, *input.units, trace);
    const auto& workflow = root.block("Electrospray");
    workflow.allow({"Coupling", "Verification"});
    const auto* verification = workflow.find("Verification");
    if ((verification != nullptr) == (workflow.find("Coupling") != nullptr))
        workflow.fail("exactly one Coupling or Verification block");
    if (verification)
    {
        verification->allow({"problem"});
        const auto& node = verification->get("problem");
        const auto problem = trace.record(node,
            node.choice({"pressureDriven", "fixedInterfaceCapillary", "fixedInterfaceDielectric"}),
            "", "fixed-material/interface verification using the existing NS loop; no interface or charge advance");
        input.problem = problem == "pressureDriven" ? ElectrosprayProblem::PressureDriven
            : problem == "fixedInterfaceCapillary" ? ElectrosprayProblem::FixedInterfaceCapillary
            : ElectrosprayProblem::FixedInterfaceDielectric;
    }
    else
    {
        const auto& coupling = workflow.block("Coupling");
        coupling.allow({"electricUpdate"});
        trace.record(coupling.get("electricUpdate"), coupling.get("electricUpdate").choice({"everyStep"}), "", "NS/material-interface/selected-charge-model/electric accepted-step sequence");
    }
    const auto& materials = modules.materials = readMaterialTable(directory + "/materials.dict", trace);
    DictionaryInput::Document hydroDocument(directory + "/hydro.dict");
    hydroDocument.root().allow({"Hydrodynamics"});
    const auto& hydro = hydroDocument.root().block("Hydrodynamics");
    const auto* inlineInterface = hydro.find("Interface");
    const auto interfacePath = directory + "/interface.dict";
    if (std::filesystem::exists(interfacePath))
    {
        if (inlineInterface)
            inlineInterface->fail("one Interface input source: remove this inline block when interface.dict is present");
        modules.flow = readFittedFlowCase(hydro, trace);
        modules.interfaceConfig = readFittedInterfaceCase(interfacePath, trace);
    }
    else
    {
        if (!inlineInterface) hydro.fail("interface.dict or the legacy inline Interface block");
        // Compatibility only: preserve original node paths/line numbers, and
        // send each module only its own fields. No old all-physics config.
        auto legacyFlow = hydro;
        auto& entries = legacyFlow.entries;
        entries.erase(std::remove_if(entries.begin(), entries.end(),
            [](const auto& node) { return node.name == "Interface"; }), entries.end());
        modules.flow = readFittedFlowCase(legacyFlow, trace);
        modules.interfaceConfig = readFittedInterfaceCase(hydro.block("Interface"), trace);
    }
    const auto& flow = modules.flow;
    const auto& electric = modules.electric = readFittedElectricCase(directory + "/electricity.dict", trace);
    if (verification)
    {
        if (electric.model != "perfectDielectric" || electric.interfaceDiscretization != "cellHarmonic")
            verification->get("problem").fail("perfectDielectric/cellHarmonic in electricity.dict for fixed verification; no conducting model is substituted");
    }
    else input.problem = electrosprayProblemForElectricModel(electric.model);
    if (electric.relaxationSI)
    {
        input.chargeRelaxation=electric.relaxationSI;
        const double densityScale=8.8541878128e-12*input.units->scales().Phi0_V/input.units->scales().L0_m;
        auto& tolerance=input.chargeRelaxation->densityChangeTolerance;
        tolerance/=densityScale;
        if (!std::isfinite(densityScale) || densityScale<=0.0 || !std::isfinite(tolerance) || tolerance<=0.0)
            throw std::invalid_argument(directory+"/electricity.dict: surface-charge density tolerance cannot be represented in code units");
    }
    modules.systems = Solvers::readDirectLinearSystems(directory + "/solvers.dict",
        {"ElectricPotential", "Momentum", "Pressure", "CoupledVelocityPressure", "InterfaceVolumeCorrection"}, trace);
    if (materials.phases.size() != 2 || !materials.phases.count(flow.liquidPhase) || !materials.phases.count(flow.gasPhase))
        throw std::invalid_argument(directory + "/materials.dict: Materials.Phases; expected exactly the two Hydro liquid/gas phase IDs");
    const auto pair = std::minmax(flow.liquidPhase, flow.gasPhase);
    if (materials.surfaceTensionSI.size() != 1 || !materials.surfaceTensionSI.count(pair))
        throw std::invalid_argument(directory + "/materials.dict: Materials.Interfaces; expected exactly the Hydro liquid/gas pair");
    if (electric.conductorPhaseId == flow.liquidPhase || electric.conductorPhaseId == flow.gasPhase)
        throw std::invalid_argument(directory + "/electricity.dict: Electricity.BoundaryConditions.Conductor.phaseId; expected a phase distinct from Hydro liquid and gas");
    for (auto& entry : trace.entries)
    {
        if (verification)
        {
            const bool pressureOnly = input.problem == ElectrosprayProblem::PressureDriven;
            const bool polarization = input.problem == ElectrosprayProblem::FixedInterfaceDielectric;
            if (entry.field == "Interface.transport" || entry.field == "Hydrodynamics.Interface.transport"
                || entry.field == "Hydrodynamics.Numerics.Conservation.phaseFluxBalance"
                || entry.field.find("Solvers.InterfaceVolumeCorrection.") == 0
                || (!polarization && entry.field == "Electricity.Numerics.polarizationConsistency")
                || (pressureOnly && (entry.field == "Hydrodynamics.Numerics.InterfaceForce.halfWidthCells"
                    || entry.field.find("Materials.Interfaces.") == 0)))
            {
                entry.used = false;
                entry.unusedReason = "not consumed by the selected fixed verification problem";
            }
        }
        if (entry.field == "Units.n0" || entry.field == "Units.rho0"
            || (electric.model == "perfectDielectric" && entry.field.size() >= 13
                && entry.field.substr(entry.field.size() - 13) == ".conductivity")
            || (flow.outletIncomingPhaseId != -1 && entry.field == "Hydrodynamics.Numerics.Conservation.outletAdmissibility"))
        {
            entry.used = false;
            entry.unusedReason = entry.field == "Units.n0" ? "no particles in this workflow"
                : entry.field == "Units.rho0" ? "retained volume-charge reference; surface charge uses epsilon0*Phi0/L0"
                : entry.field.find("conductivity") != std::string::npos ? "perfect dielectric: no conduction current is advanced"
                : "all configured pressure ports allow an explicit incoming reservoir state";
        }
    }
    return input;
}
