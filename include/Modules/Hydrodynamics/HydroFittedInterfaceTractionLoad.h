#pragma once

#include "Modules/Hydrodynamics/HydroFittedInterfaceLoadGeometry.h"
#include "Core/Field/VEdgeField.h"
#include "Core/Field/HEdgeField.h"

namespace Quantum::Hydrodynamics
{
// Result is in the supplied traction units per length. The caller may apply
// one mechanical normalization after loading, preserving its original order.
struct HydroFittedInterfaceTractionLoad
{
    explicit HydroFittedInterfaceTractionLoad(const QuadtreeMesh& mesh);
    VEdgeField<double> forceX;
    HEdgeField<double> forceY;
    double maxForce = 0.0, totalForceX = 0.0, totalForceY = 0.0;
};

// One normal/tangential pair per actual ordered sampling point. Uses the full
// immutable primary/transverse support prepared by the receiving Hydro domain:
// tau_n*G(H) + tau_t*R(G(H)). No division by interface-normal components and no
// second curvature, pressure, electric or charge computation. All validation
// and finite/underflow checks finish before the candidate is returned.
HydroFittedInterfaceTractionLoad buildFittedInterfaceTractionLoad(const QuadtreeMesh& mesh,
    const HydroFittedDomain& domain, const HydroFittedInterfaceLoadGeometry& geometry,
    const std::shared_ptr<const InterfaceTraceGeometry>& traceGeometry,
    const std::vector<Vector2D>& normalTangentialTractions);
}
