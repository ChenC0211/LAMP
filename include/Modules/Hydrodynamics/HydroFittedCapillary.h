#pragma once

#include <vector>
#include "Core/Field/VEdgeField.h"
#include "Core/Field/HEdgeField.h"
#include "Modules/Hydrodynamics/HydroInterfaceLoadGeometry.h"
#include "Modules/Hydrodynamics/HydroForceSource.h"

class QuadtreeMesh;
class HydroFittedDomain;
class LevelSetInterfacePanels;
namespace Quantum::Hydrodynamics { struct HydroFittedInterfaceLoadGeometry; }

// A single immutable-interface force assembly, not an interface/charge time
// integrator. All returned force densities use the supplied fitted NS duals.
struct HydroFittedCapillaryState
{
    VEdgeField<double> forceX;
    HEdgeField<double> forceY;
    std::shared_ptr<const HydroForceSource> source;
    std::vector<double> indicator; // Native cells; inactive payload is NaN.
    std::vector<HydroInterfaceLoadGeometry::Face> faces;
    double halfWidth = 0.0, sigma = 0.0;
    double maxAbsComponent = 0.0, integralX = 0.0, integralY = 0.0;
    double maxGradientIdentityError = 0.0;
    int endpointProjections = 0;
    explicit HydroFittedCapillaryState(const QuadtreeMesh& mesh);
};

// Builds the same compact C2 H and -sigma*kappa*G(H) as HydroBalancedForce.
// G(H) is the ACTUAL fitted pressure-load matrix C applied to H, divided by
// actual dualVolume, not a new h or full-rectangle pressure stencil. Pressure
// ports must be separated from the band; constrained solid/axis rows get no
// force. Wall reaction/contact-line support are not part of the reported sum.
HydroFittedCapillaryState buildFittedCapillaryForce(const QuadtreeMesh& mesh,
    const HydroFittedDomain& domain, const LevelSetHermiteGeometry::Snapshot& snapshot,
    const std::vector<Vector2D>& boundaryContacts, double halfWidth, double sigma);

// Same capillary arithmetic, consuming an already validated fitted load
// geometry (uniform or native). Fills an unaccepted workspace without another
// force-field copy. On failure discard the workspace; source provenance is
// cleared and must be bound by the caller only after successful assembly.
void assembleFittedCapillaryForce(const QuadtreeMesh& mesh, const HydroFittedDomain& domain,
    const Quantum::Hydrodynamics::HydroFittedInterfaceLoadGeometry& geometry,
    double sigma, HydroFittedCapillaryState& destination);

// Explicit native-leaf entry. Reuses the same source/panels and native load
// geometry as the charged path; never queries uniform logical-grid metadata.
// Transverse geometry is retained, but capillarity consumes only the primary
// pressure gradient. No electric traction, wall reaction or contact law.
HydroFittedCapillaryState buildNativeFittedCapillaryForce(const QuadtreeMesh& mesh,
    const HydroFittedDomain& domain, std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> snapshot,
    const LevelSetInterfacePanels& panels, double halfWidth, double sigma);
