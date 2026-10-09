#pragma once

#include "Core/Materials/MaterialTable.h"
#include "Modules/Hydrodynamics/FittedFlowCaseConfig.h"
#include "Modules/Interface/Compatibility/InterfaceCaseConfig.h"
#include "Modules/Electricity/FittedElectricCaseConfig.h"
#include "Solvers/DirectLinearSystemInput.h"
#include "Modules/Hydrodynamics/HydroFittedNumerics.h"
#include "Runtime/RunControls.h"
#include <map>
#include <string>

// Composition of owner-defined input values, not another flattened physics
// configuration. Setup borrows these values; no field/mesh/clock state lives here.
struct ElectrosprayModuleInputs
{
    MaterialTable materials;
    FittedFlowCaseConfig flow;
    FittedInterfaceCaseConfig interfaceConfig;
    FittedElectricCaseConfig electric;
    std::map<std::string, Solvers::DirectLinearSystemInput> systems;
    RunControls controls;
};

// Bind the existing Hydro controls once per solve/transfer setup. Linear
// residuals stay separate from the module's physical acceptance tolerances.
inline HydroFittedNumerics electrosprayHydroNumerics(const ElectrosprayModuleInputs& input)
{
    HydroFittedNumerics result;
    result.momentumLinear = input.systems.at("Momentum").relativeTolerance;
    result.pressureLinear = input.systems.at("Pressure").relativeTolerance;
    result.coupledLinear = input.systems.at("CoupledVelocityPressure").relativeTolerance;
    result.momentumBalance = input.flow.momentumBalance;
    result.localContinuity = input.flow.localContinuity;
    result.outletAdmissibility = input.flow.outletAdmissibility;
    return result;
}
