#pragma once

#include "Modules/Hydrodynamics/HydroDofLayout.h"
#include "Modules/Interface/GeometricVofTransport.h"

namespace LAMP::Hydrodynamics
{
// Copies the actual Hydro native normal means into one VOF temporal input.
// The two native values must share their layout and the old publication's
// mesh binding. Values have no physical clock; nextEnd explicitly declares
// the endpoint guess interval. No Hydro or Interface state is accepted here.
Interface::GeometricVofAdvectorSlab prepareHydroVofAdvectorSlab(
    const QuadtreeMesh&,const Interface::GeometricVofPublication& oldVof,
    Interface::GeometricVofTimeLayer nextEnd,
    const HydroNativeValues& beginValues,const HydroNativeValues& endGuessValues);
}