#pragma once

struct ElectrosprayWorkflow;

struct ElectrosprayModuleInputs;

#include <memory>
#include <optional>
#include <vector>
#include "Electrospray/Coupling/ElectrosprayForces.h"
#include "Modules/Electricity/FittedChargeEvolution.h"

struct ElectrosprayInterfaceMotionCandidate;

// One accepted or candidate conducting time level. Charge owns the retained
// scalar/panels; electric contains that level's potential and phase currents;
// forces contains that level's shared capillary/Maxwell native-face loads.
// No borrowed SurfaceChargeSource callback is stored in this state.
struct ElectrosprayConductingState
{
    Quantum::Electricity::FittedChargeState electrical;
    ElectrosprayConductingForces forces;
    explicit ElectrosprayConductingState(const QuadtreeMesh& mesh);
};

// Extend the module's numerical report with App force-assembly timing only.
struct ElectrosprayConductingStepReport : Quantum::Electricity::FittedChargeStepReport
{
    double interfaceLoads_s = 0.0;
};

struct ElectrosprayConductingCandidate
{
    // Commit exchanges the owner without copying field arrays. Afterwards this
    // candidate owns the retired time level until publication finishes.
    std::unique_ptr<ElectrosprayConductingState> state;
    Quantum::Electricity::SurfaceChargeFieldStep chargeStep;
    ElectrosprayConductingStepReport report;
    explicit ElectrosprayConductingCandidate(const QuadtreeMesh& mesh);
};

// Run before the first NS solve and again at candidate preparation. The
// existing explicit-charge restriction is dt*max(sigma/epsilon)<=0.25;
// this is a step restriction, not a new physical-accuracy acceptance gate.
// Closure fields are zero until an actual charge candidate is evaluated.
ElectrosprayConductingStepReport preflightElectrosprayConduction(
    const ElectrosprayModuleInputs& inputs, const ElectrosprayWorkflow& workflow,
    const UnitsContext& units);

// App time-level binding only. Default explicit branch: add OLD current once, remap through the
// motion candidate's retained material map onto its FINAL volume-corrected
// scalar, then solve the charged potential/current and shared loads there.
// The old charged potential is the explicit warm-start input; initialBinding
// contributes only the fixed active electrical graph and actual boundary IDs.
// No extra tangential advection, NS solve, scalar transport or state commit.
// Initial state assembly belongs to the App: move its already solved charge,
// electric result and forces into ElectrosprayConductingState, without a solve.
// On failure oldState and motion are unchanged; no candidate is published.
// With explicit relaxation controls, remap current/older material histories
// together, iterate the existing P3 electric subproblem with the module time
// kernel, then use its coeval current in electrode accounting and Maxwell loads.
// No density-expansion or extra tangential-convection term is added.
// targetTime is supplied by the App's sole clock, not accumulated from metadata.
ElectrosprayConductingCandidate prepareElectrosprayConductingStep(
    const QuadtreeMesh& mesh, const HydroFittedDomain& domain,
    const ElectrosprayElectricState& initialBinding,
    const ElectrosprayModuleInputs& inputs, const ElectrosprayWorkflow& workflow,
    const UnitsContext& units, const ElectrosprayConductingState& oldState,
    const ElectrosprayInterfaceMotionCandidate& motion, double targetTime,
    const Quantum::Electricity::SurfaceChargeRelaxationOptions* relaxation = nullptr);
