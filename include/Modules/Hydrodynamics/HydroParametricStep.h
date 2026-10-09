#pragma once
#include "Modules/Hydrodynamics/HydroStep.h"
#include "Modules/Hydrodynamics/HydroUniformLoad.h"
#include "Modules/Interface/Compatibility/InterfaceParametricMotion.h"
#include <functional>

namespace Quantum::Hydrodynamics {
struct HydroParametricStepOptions {
    double liquidDensity=1,gasDensity=1,liquidViscosity=0,gasViscosity=0,surfaceTension=0;
    int maximumIterations=16;
    double relativeVelocityTolerance=1e-9,absoluteVelocityTolerance=1e-11;
    InterfaceParametricProjection::Options geometry;
    HydroSolveOptions solver;
};
struct HydroParametricStepResult {
    InterfaceState interface;
    HydroSolveReport hydro;
    ParametricCapillaryReport midpointSurface;
    int iterations=0;
    double midpointVelocityDefect=0,inventoryBalance=0,kineticEnergy=0;
};
// Optional synchronous observation of the actual trial midpoint velocity.
// References are read-only and valid for this call only. The caller owns IO;
// observation happens before half/full interface preparation or fluid mutation.
using HydroParametricStepObserver=std::function<void(int,const VEdgeField<double>&,const HEdgeField<double>&)>;
// One physical step, with the caller's time and mesh. The sole coupled unknown
// is u_mid=(u_old+u_new)/2. Each Picard iterate regenerates half/full geometry,
// independent FV fluxes, mass and force from that same velocity. Only a
// converged iterate is accepted; failure restores Hydro and retains old
// Interface. No dt selection, retry, output or App loop belongs here.
HydroParametricStepResult advanceParametricCapillaryStep(QuadtreeMesh&,const NativeLeafDomain&,
    ::Hydrodynamics&,const InterfaceState&,const InterfaceStepTime&,const HydroParametricStepOptions&,
    const HydroParametricStepObserver& observer={});
}
