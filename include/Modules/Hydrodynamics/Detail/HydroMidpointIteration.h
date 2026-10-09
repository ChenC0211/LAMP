#pragma once
#include "Modules/Hydrodynamics/HydroMidpointStep.h"
#include <functional>

namespace Quantum::Hydrodynamics::Detail {
// Common nonlinear solve/acceptance. A trial assembles one physical equation
// from the supplied arithmetic midpoint; it returns the real Hydro candidate.
// Publish updates any endpoint material/interface owner only after convergence.
// Each unsuccessful iterate or exception restores the original Hydro state.
using MidpointTrial=std::function<HydroStepCandidate(int,const VEdgeField<double>&,const HEdgeField<double>&)>;
HydroMidpointStepResult solveMidpointIterations(QuadtreeMesh&,::Hydrodynamics&,
    const HydroMidpointIterationOptions&,const MidpointTrial&,const std::function<void()>& publish);
}
