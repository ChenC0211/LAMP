#pragma once
#include "Modules/Hydrodynamics/HydroStep.h"

namespace Quantum::Hydrodynamics {
struct HydroMidpointIterationOptions {
    int maximumIterations=16;
    double relativeVelocityTolerance=1e-9,absoluteVelocityTolerance=1e-11;
};
struct HydroMidpointStepResult {
    HydroSolveReport hydro;
    int iterations=0;
    double midpointVelocityDefect=0;
};
// One physical incompressible single-phase step, with constant density and
// frozen owned viscosity (which may vary in space). Uses the same compatible
// midpoint momentum equation and Picard acceptance as the parametric two-phase
// step. Uniform closed XY/RZ; no swirl. Force declares its actual midpoint
// time/source. No clock, dt selection, retry or output is owned here.
HydroMidpointStepResult advanceSinglePhaseMidpointStep(QuadtreeMesh&,const NativeLeafDomain&,
    ::Hydrodynamics&,const HydroStepTime&,const HydroSolveOptions&,const HydroForceDensity&,
    const HydroMidpointIterationOptions& = {});
}
