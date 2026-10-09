#include "Electrospray/Coupling/ElectrosprayConducting.h"
#include "Electrospray/Setup/ElectrosprayWorkflow.h"
#include "Electrospray/Setup/ElectrosprayModuleInputs.h"
#include "Electrospray/IO/ElectrosprayForceDiagnostics.h"
#include "Electrospray/IO/ElectrosprayChargeDiagnostics.h"
#include "Electrospray/Setup/ElectrosprayElectrostatics.h"
#include "Electrospray/Coupling/ElectrosprayInterfaceMotion.h"
#include "Modules/Hydrodynamics/HydroFittedDomain.h"
#include "Core/Units/UnitsContext.h"
#include <algorithm>
#include <chrono>
#include <limits>
#include <stdexcept>

namespace {
Quantum::Electricity::FittedChargeParameters parameters(
    const ElectrosprayModuleInputs& inputs, const UnitsContext& units)
{
    return bindElectrosprayConductingParameters(inputs.electric,
        inputs.materials.phases.at(inputs.flow.liquidPhase),
        inputs.materials.phases.at(inputs.flow.gasPhase),
        inputs.systems.at("ElectricPotential").relativeTolerance,units);
}
}
ElectrosprayConductingState::ElectrosprayConductingState(const QuadtreeMesh& mesh) : forces(mesh) {}
ElectrosprayConductingCandidate::ElectrosprayConductingCandidate(const QuadtreeMesh& mesh)
    : state(std::make_unique<ElectrosprayConductingState>(mesh)) {}

ElectrosprayConductingStepReport preflightElectrosprayConduction(
    const ElectrosprayModuleInputs& inputs, const ElectrosprayWorkflow& workflow,
    const UnitsContext& units)
{
    try
    {
        if (workflow.problem != ElectrosprayProblem::Conducting)
            throw std::runtime_error("conducting evolution requires the coupled workflow");
        ElectrosprayConductingStepReport report;
        static_cast<Quantum::Electricity::FittedChargeStepReport&>(report) =
            Quantum::Electricity::preflightFittedChargeStep(parameters(inputs,units),inputs.controls.dt);
        return report;
    }
    catch (const std::exception& error)
    {
        throw std::runtime_error(std::string("[ElectrosprayConducting] preflight: ") + error.what());
    }
}

ElectrosprayConductingCandidate prepareElectrosprayConductingStep(
    const QuadtreeMesh& mesh, const HydroFittedDomain& domain,
    const ElectrosprayElectricState& initialBinding,
    const ElectrosprayModuleInputs& inputs, const ElectrosprayWorkflow& workflow,
    const UnitsContext& units, const ElectrosprayConductingState& oldState,
    const ElectrosprayInterfaceMotionCandidate& motion, double targetTime,
    const Quantum::Electricity::SurfaceChargeRelaxationOptions* relaxation)
{
    const char* stage = "time-level-binding";
    try
    {
        preflightElectrosprayConduction(inputs,workflow,units);
        domain.validateMeshIdentity(mesh);
        if (!domain.axisymmetric() || !oldState.forces.total.source()
            || oldState.electrical.completedStep()==std::numeric_limits<int>::max()
            || oldState.forces.total.source()->sourceTime()!=oldState.electrical.time()
            || motion.motion.report.dt != inputs.controls.dt)
            throw std::runtime_error("conducting step requires matching RZ flow, force time and material trajectory");
        const Quantum::Electricity::FittedChargeProblem problem{
            initialBinding.dielectric.graph(),initialBinding.dielectric.boundaries(),parameters(inputs,units)};
        const Quantum::Electricity::FittedChargeMotion materialMotion{
            motion.state.surface.geometry(),motion.state.surface.panels(),
            motion.state.surface.pinnedCoordinates(),motion.motion.materialMapFor(
                oldState.electrical.electric().potential.source(),motion.state.surface.geometry(),
                {oldState.electrical.completedStep()+1,oldState.electrical.time(),targetTime,inputs.controls.dt}),
            motion.motion.report.dt,oldState.electrical.time(),targetTime};
        stage = "electric-charge-candidate";
        auto electrical = Quantum::Electricity::prepareFittedChargeStep(
            mesh,problem,oldState.electrical,materialMotion,relaxation,reportElectrosprayChargeIteration);
        ElectrosprayConductingCandidate result(mesh);
        result.state->electrical = std::move(electrical.state);
        result.chargeStep = std::move(electrical.chargeStep);
        static_cast<Quantum::Electricity::FittedChargeStepReport&>(result.report) = std::move(electrical.report);
        stage = "candidate-shared-interface-loads";
        const auto& flowConfig = inputs.flow;
        const auto& electricConfig = inputs.electric;
        const auto& liquid = inputs.materials.phases.at(flowConfig.liquidPhase);
        const auto& gas = inputs.materials.phases.at(flowConfig.gasPhase);
        const double surfaceTensionSI = inputs.materials.surfaceTensionSI.at(std::minmax(flowConfig.liquidPhase,flowConfig.gasPhase));
        const auto loadsStarted = std::chrono::steady_clock::now();
        // The source retains this immutable electrical candidate across the
        // synchronous load consumer and any subsequent handle movement.
        const auto source = result.state->electrical.source();
        result.state->forces = buildElectrosprayConductingForces(mesh,domain,motion.state.surface.geometry(),
            result.state->electrical.charge()->geometry(),source,result.state->electrical.electric(),
            liquid,gas,surfaceTensionSI,electricConfig,flowConfig.capillaryHalfWidthCells,units,targetTime);
        result.report.interfaceLoads_s = std::chrono::duration<double>(
            std::chrono::steady_clock::now()-loadsStarted).count();
        reportElectrosprayConductingForces(result.state->forces,units);
        return result;
    }
    catch (const std::exception& error)
    {
        throw std::runtime_error(std::string("[ElectrosprayConducting] ") + stage + ": " + error.what());
    }
}
