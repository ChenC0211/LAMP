#include "Electrospray/Compatibility/ElectrosprayInputAdapters.h"
#include "Electrospray/Setup/ElectrosprayModuleInputs.h"
#include <algorithm>

ElectrosprayWorkflow adaptLegacyElectrosprayWorkflow(const ElectrosprayConfig& run)
{
    ElectrosprayWorkflow result;
    if (run.physicsMode == "coupled") result.problem = ElectrosprayProblem::Conducting;
    else if (run.physicsMode == "movingDielectric") result.problem = ElectrosprayProblem::MovingDielectric;
    else if (run.physicsMode == "frozenMaterialHydro") result.problem = ElectrosprayProblem::PressureDriven;
    else if (run.physicsMode == "frozenInterfaceCapillary") result.problem = ElectrosprayProblem::FixedInterfaceCapillary;
    else if (run.physicsMode == "frozenInterfaceEhd") result.problem = ElectrosprayProblem::FixedInterfaceDielectric;
    else throw std::runtime_error("Electrospray physicsMode must be coupled, frozenMaterialHydro, frozenInterfaceCapillary, frozenInterfaceEhd or movingDielectric");
    result.restartLoop = run.restartLoop;
    result.initialAdaptIterations = run.initialAdaptIterations;
    result.initializeLevelSets = run.initializeLevelSets;
    result.legacyIterationLog = run.detailedIterationLog;
    return result;
}

namespace
{
FittedFlowCaseConfig adaptLegacyElectrosprayFlow(const ElectrosprayPhysicalConfig& p)
{
    FittedFlowCaseConfig f;
    f.liquidPhase = p.liquid.phaseId; f.gasPhase = p.gas.phaseId;
    f.supply = p.supply; f.outlet = p.outlet; f.ambient = p.ambient;
    f.outletIncomingPhaseId = p.outletIncomingPhaseId;
    f.incomingMomentumModel = p.incomingMomentumModel; f.solidWallModel = p.solidWallModel;
    f.momentumBalance = p.numerics.flow.momentumBalance;
    f.localContinuity = p.numerics.flow.localContinuity;
    f.outletAdmissibility = p.numerics.flow.outletAdmissibility;
    f.phaseFluxBalance = p.numerics.phaseFluxBalance;
    return f;
}

} // namespace

ElectrosprayModuleInputs adaptLegacyElectrosprayInputs(
    const ElectrosprayPhysicalConfig& p, const ElectrosprayConfig& run,
    const ElectrosprayWorkflow& workflow)
{
    ElectrosprayModuleInputs input;
    input.materials.phases.emplace(p.liquid.phaseId, p.liquid);
    input.materials.phases.emplace(p.gas.phaseId, p.gas);
    input.materials.surfaceTensionSI.emplace(std::minmax(p.liquid.phaseId, p.gas.phaseId), p.surfaceTensionSI);
    input.flow = adaptLegacyElectrosprayFlow(p);
    input.flow.capillaryHalfWidthCells = run.capillaryHalfWidthCells;
    auto& interface = input.interfaceConfig;
    interface.fluidInterfaceLevelSet = p.fluidInterfaceLevelSet;
    interface.movingBoundaryId = p.movingBoundaryId;
    interface.contactLineModel = p.contactLineModel;
    interface.pinnedContactPointIds = p.pinnedContactPointIds;
    auto& e = input.electric;
    e.model = p.model; e.interfaceDiscretization = p.interfaceDiscretization;
    e.surfaceChargeModel = workflow.conducting() ? "explicitConductionMaterialTransport" : "none";
    e.conductorPhaseId = p.conductorPhaseId; e.conductorPotentialSI = p.conductorPotentialSI;
    e.initialSurfaceChargeSI = p.initialSurfaceChargeSI; e.domainBoundaries = p.domainElectricBoundaries;
    e.potentialBounds = p.numerics.potentialBounds; e.fluxChecks = p.numerics.electricFlux;
    e.polarizationConsistency = p.numerics.polarizationConsistency; e.chargeClosure = p.numerics.chargeClosure;
    input.systems = {
        {"Momentum", {"SuperLU", p.numerics.flow.momentumLinear}},
        {"Pressure", {"SuperLU", p.numerics.flow.pressureLinear}},
        {"CoupledVelocityPressure", {"SuperLU", p.numerics.flow.coupledLinear}},
        {"ElectricPotential", {"SuperLU", p.numerics.electricPotentialLinear}},
        {"InterfaceVolumeCorrection", {"SuperLU", p.numerics.volumeCorrectionLinear}}};
    auto& c = input.controls;
    c.dt = run.dt; c.endTime = run.finalTime; c.writeInterval = run.writeDt;
    c.writeEverySteps = run.writeEverySteps; c.logEverySteps = run.logEvery;
    c.outputPrefix = run.outputName; c.writeFields = run.writeFields;
    c.diagnosticFields = run.diagnosticFields; c.writeInitial = run.writeInitial;
    c.writeFirst = run.writeFirst.value_or(workflow.movesInterface());
    c.writeFinal = run.writeFinal;
    // Legacy field detail and iteration detail may differ; the legacy report
    // wrapper preserves these separate flags instead of inferring one from the other.
    c.detailedLog = run.verboseDiagnostics;
    return input;
}
