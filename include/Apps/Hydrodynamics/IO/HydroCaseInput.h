#pragma once

#include "Hydrodynamics/HydroCase.h"
#include "Modules/Hydrodynamics/HydrodynamicsConfig.h"
#include "Modules/Interface/Compatibility/InterfaceNodalConfig.h"

struct HydroCaseInput
{
    HydrodynamicsConfig flow;
    InterfaceNodalConfig interface;
    std::shared_ptr<const HydroCaseSchedule> schedule;
};

bool usesSeparateHydroInputs();
HydroCaseInput readHydroCaseInput(bool coupled);
