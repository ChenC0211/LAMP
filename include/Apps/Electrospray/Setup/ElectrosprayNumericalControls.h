#pragma once

#include "Modules/Hydrodynamics/HydroFittedNumerics.h"
#include "Modules/Electricity/DielectricFluxChecks.h"

// App-owned composition of independently consumed module settings. No reader,
// mesh, mutable state or physical material data lives here. The Compatibility
// adapter expands the old single tolerance once; numerical consumers must not
// derive physical acceptance gates from a linear solver tolerance again.
struct ElectrosprayNumericalControls
{
    HydroFittedNumerics flow;
    DielectricFluxChecks electricFlux;
    double electricPotentialLinear = 0.0;
    double volumeCorrectionLinear = 0.0;
    double phaseFluxBalance = 0.0;
    double potentialBounds = 0.0;
    double polarizationConsistency = 0.0;
    double chargeClosure = 0.0;
};
