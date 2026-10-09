#pragma once

#include "Modules/Hydrodynamics/HydroInterfaceLoadGeometry.h"
#include "Core/Geometry/LevelSetBoundaryConditions.h"
#include "Modules/Hydrodynamics/HydroBoundary.h"
#include "Modules/Interface/Compatibility/InterfaceNodalConfig.h"

template<typename T> class VEdgeField;
template<typename T> class HEdgeField;

// Adapter from explicit Hydro boundary data to the shared scalar/geometry
// interfaces. No case shape, contact-angle law, mesh owner or clock lives here.
namespace HydroInterfaceBoundarySupport
{
LevelSetBoundaryConditions scalarConditions(const InterfaceNodalConfig::ScalarBoundary& boundary,
    const std::vector<HydroBoundaryPatch>& patches);
void configureLoads(const std::vector<Vector2D>& contacts, const std::vector<HydroBoundaryPatch>& patches,
    Vector2D lower, Vector2D upper,
    HydroInterfaceLoadGeometry::Options& options);

// Actual t^n phase geometry with the supplied native velocity. A diagnostic
// flux, not a conservative density/VOF transport discretization.
double negativePhaseOutflow(const QuadtreeMesh& mesh,
    const LevelSetHermiteGeometry::Snapshot& scalar, const VEdgeField<double>& ux,
    const HEdgeField<double>& uy, bool axisymmetric);
}
