#pragma once

struct ElectrosprayWorkflow;

struct ElectrosprayModuleInputs;

#include <memory>
#include <functional>
class HydroFittedDomain;
class HydroForceSource;

class QuadtreeMesh;
class UnitsContext;
template<typename T> class VEdgeField;
template<typename T> class HEdgeField;
struct ElectrosprayElectricState;
struct ElectrosprayInterfaceState;
struct ElectrosprayConductingState;
struct ElectrosprayMeshAdaptation;
struct ElectrosprayAcceptedStateView;
// Called after physical acceptance and successful publication, before the
// next step/AMR. May inspect the live state synchronously; no numerical input
// is supplied by an observer. An exception stops the already accepted run.
using ElectrosprayAcceptedObserver = std::function<void(int, const ElectrosprayAcceptedStateView&)>;

// Explicit development mode; never invoked by default/full coupled mode.
// Uses the real pressure/solid boundary data and shared Hydrodynamics class,
// but freezes materials and omits interface advection/capillary/electric force.
// Success means only this declared pressure-driven NS interval completed.
void runElectrosprayFrozenMaterialHydro(QuadtreeMesh& mesh, const HydroFittedDomain& domain,
    const ElectrosprayElectricState& initial, const ElectrosprayModuleInputs& inputs,
    const ElectrosprayWorkflow& workflow, const UnitsContext& units);

// Explicit frozenInterfaceCapillary mode: the SAME fixed-material NS loop
// consumes a retained interface's native-face force densities in code units.
// Geometry/materials and these forces stay fixed; pressure data remain real.
// No electric force, charge transport or interface motion is advanced here.
void runElectrosprayFrozenCapillaryHydro(QuadtreeMesh& mesh, const HydroFittedDomain& domain,
    const ElectrosprayElectricState& initial, const ElectrosprayModuleInputs& inputs,
    const ElectrosprayWorkflow& workflow, const UnitsContext& units,
    const VEdgeField<double>& forceX, const HEdgeField<double>& forceY,
    std::shared_ptr<const HydroForceSource> source);

// Explicit frozenInterfaceEhd mode, using the SAME fixed-material NS loop.
// Caller supplies the sum of capillary and q=0 cellHarmonic polarization
// native-face force densities in code units. Initial interface, material and
// electric field remain frozen: no charge advance or electrospray trajectory.
void runElectrosprayFrozenEhdHydro(QuadtreeMesh& mesh, const HydroFittedDomain& domain,
    const ElectrosprayElectricState& initial, const ElectrosprayModuleInputs& inputs,
    const ElectrosprayWorkflow& workflow, const UnitsContext& units,
    const VEdgeField<double>& totalForceX, const HEdgeField<double>& totalForceY,
    std::shared_ptr<const HydroForceSource> source);

// Explicit movingDielectric development mode: the same NS plus free-interface,
// material and q=0 electrostatic updates. Conductive charge is NOT advanced;
// this is never silently substituted for the requested leaky-dielectric mode.
void runElectrosprayMovingDielectricHydro(QuadtreeMesh& mesh, const HydroFittedDomain& domain,
    const ElectrosprayElectricState& initial, const ElectrosprayModuleInputs& inputs,
    const ElectrosprayWorkflow& workflow, const UnitsContext& units, ElectrosprayInterfaceState& interfaceState);

// Same NS/interface/volume loop, with one conservative material charge update
// and a new charged P3 field/full Maxwell load before each state is accepted.
// Takes ownership of the already solved initial charge/electric/load bundle.
// A post-NS failure stops the run; it does not roll back Hydro's internal state.
namespace Quantum::Electricity { struct SurfaceChargeRelaxationOptions; }
void runElectrosprayConductingHydro(QuadtreeMesh& mesh, const HydroFittedDomain& domain,
    const ElectrosprayElectricState& binding, const ElectrosprayModuleInputs& inputs,
    const ElectrosprayWorkflow& workflow, const UnitsContext& units, ElectrosprayInterfaceState& interfaceState,
    std::unique_ptr<ElectrosprayConductingState> initialState,
    const Quantum::Electricity::SurfaceChargeRelaxationOptions* relaxation = nullptr,
    ElectrosprayMeshAdaptation* adaptation = nullptr,
    const ElectrosprayAcceptedObserver& acceptedObserver = {});
