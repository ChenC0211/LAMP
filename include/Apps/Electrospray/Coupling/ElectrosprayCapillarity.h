#pragma once

#include "Modules/Hydrodynamics/HydroFittedCapillary.h"

class InterfaceState;
class UnitsContext;

// Unit/case adapter only; no file or console output. Geometry and force
// arithmetic stay in Core / Hydro; no source-circle curvature enters the load.
HydroFittedCapillaryState assembleElectrosprayCapillarity(const QuadtreeMesh& mesh,
    const HydroFittedDomain& domain, const InterfaceState& interfaceState,
    double referenceDensitySI, double surfaceTensionSI, double capillaryHalfWidthCells,
    const UnitsContext& units, double sourceTime, HydroForceEvolution evolution);
