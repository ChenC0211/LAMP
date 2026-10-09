#pragma once

#include <memory>
#include "Hydrodynamics/Coupling/HydroCaseCoupling.h"
struct DielectricConfig;

// Owns electric/charge time levels, never the mesh, Hydro state or clock.
// The App-owned callbacks preserve the existing explicit force/charge ordering.
std::unique_ptr<HydroCaseCoupling> makeDielectricCoupling(DielectricConfig config);
// Explicit research case: charged reflection column and grounded radial wall.
// Axis-cap closure remains a separately reported, unvalidated event model.
std::unique_ptr<HydroCaseCoupling> makeConservativeColumnCoupling(DielectricConfig config,bool axisCapModel);
