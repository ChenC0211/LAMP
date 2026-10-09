#pragma once

#include "Modules/Interface/Compatibility/InterfaceTransportInput.h"
#include "Modules/Interface/Compatibility/InterfaceParametricTransport.h"
#include "Modules/Hydrodynamics/HydroBoundary.h"
#include "Modules/Hydrodynamics/HydroCompatibleMomentSpline.h"

class HydroFittedDomain;
class HydroFittedVelocitySampler;
class HydroSplineVelocity;
namespace Quantum::Hydrodynamics {class HydroCompatibleMomentSpline;}
class NativeLeafDomain;
class QuadtreeMesh;
namespace LevelSetHermiteGeometry { class Snapshot; }
template<typename T> class VEdgeField;
template<typename T> class HEdgeField;

// Exports existing Hydro boundary/normal-flux data to Interface. These are
// explicit borrowed views, not another reconstruction, scalar transport or
// volume solver. Their Hydro objects must remain alive and immutable during
// the synchronous Interface call; no candidate retains either view.
namespace HydroInterfaceTransport
{
InterfaceParametricTransport::Flow frozenParametricFlow(const Quantum::Hydrodynamics::HydroCompatibleMomentSpline::NativeVelocity&,int blocks=8);
InterfaceParametricTransport::Flow timeLinearParametricFlow(const Quantum::Hydrodynamics::HydroCompatibleMomentSpline::NativeVelocity& previous,
    const Quantum::Hydrodynamics::HydroCompatibleMomentSpline::NativeVelocity& current,double stepDuration,int blocks=8);
InterfaceTransportDomain domainView(const HydroFittedDomain& domain);
InterfaceNormalFluxView normalFlux(const HydroFittedVelocitySampler& velocity);
InterfaceNormalFluxView normalFlux(const HydroSplineVelocity& velocity);
InterfaceNormalFluxView normalFlux(const Quantum::Hydrodynamics::HydroCompatibleMomentSpline&);
InterfaceParametricTransport::Flow frozenParametricFlow(const Quantum::Hydrodynamics::HydroCompatibleMomentSpline&,int blocks=8);
InterfaceParametricTransport::Flow timeLinearParametricFlow(const Quantum::Hydrodynamics::HydroCompatibleMomentSpline& previous,
    const Quantum::Hydrodynamics::HydroCompatibleMomentSpline& current,double stepDuration,int blocks=8);
// Explicit C1 flux-preserving reconstruction; the whole domain is one C1
// trajectory chart, so native spline knots do not trigger time restarts.
InterfaceParametricTransport::Flow frozenParametricFlow(const HydroSplineVelocity& velocity,int blocks=8);
// Linear physical-time interpolation of TWO immutable snapshots. Flow time
// is elapsed from the previous layer and lies in [0,stepDuration]. Pair with
// InterfaceParametricMotion::prepareTimeLinear and both snapshots' flux views.
InterfaceParametricTransport::Flow timeLinearParametricFlow(const HydroSplineVelocity& previous,
    const HydroSplineVelocity& current,double stepDuration,int blocks=8);
// Autonomous first-jet map for the current parametric XY producer. Supply the
// completeVelocity snapshot for this same uniform mesh/domain, and pair with
// normalFlux of that snapshot. Owns the native smooth-patch index; borrows the
// immutable velocity. It neither blends time layers nor extends physical walls.
InterfaceParametricTransport::Flow frozenParametricFlow(const QuadtreeMesh& mesh,
    const NativeLeafDomain& domain, const HydroFittedVelocitySampler& velocity, int blocks=8);
// Complete XY/RZ transport using Hydro zero-shear boundary traces. Mixed
// native leaves currently support the closed box; root grids also accept ports.
// Empty patches select the closed zero-normal/free-slip box. RZ alone uses
// the odd radial/even axial trace at r=0; x=0 has no special meaning in XY.
// Shared vertex traces are reconstructed from native NS fields; axial jets
// use three actual vertex coordinates on root grids; mixed leaves reuse the
// native least-squares value/derivative traces. Both trajectories and fluxes must
// then use the returned immutable sampler. Discontinuous prescribed normal
// traces at a shared vertex reject rather than being smoothed across an inlet lip.
HydroFittedVelocitySampler completeVelocity(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    const VEdgeField<double>& ux, const HEdgeField<double>& uy,
    const std::vector<HydroBoundaryPatch>& patches);
// Borrows mesh/domain and owns a copy of the physical patches. Zero normal
// boundaries are symmetry/slip, not stationary no-slip contacts. Inlet/outlet
// phase history must be supplied separately by the Interface caller.
InterfaceTransportDomain domainView(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    const std::vector<HydroBoundaryPatch>& patches);

// Original native-average pressure-port diagnostic. Keep its multiply/order
// and strict pure-phase certificate; do not replace it with reconstructed
// face integration just because both return a volume rate.
long double outwardPhaseFlux(const QuadtreeMesh& mesh, const HydroFittedDomain& domain,
    const LevelSetHermiteGeometry::Snapshot& snapshot,
    const VEdgeField<double>& ux, const HEdgeField<double>& uy);
}
