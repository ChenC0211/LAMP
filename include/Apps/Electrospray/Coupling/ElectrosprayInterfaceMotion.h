#pragma once

#include <memory>
#include "Electrospray/Setup/ElectrosprayInterface.h"
#include "Modules/Interface/Compatibility/InterfaceMotion.h"

class Hydrodynamics;

// App binding plus the module's one result. No duplicate transport report,
// phase ledger, geometry, correction or material-trajectory implementation.
struct ElectrosprayInterfaceMotionCandidate
{
    ElectrosprayInterfaceState state;
    InterfaceMotionResult motion;
    explicit ElectrosprayInterfaceMotionCandidate(const QuadtreeMesh& mesh) : state(mesh) {}
};

// Validates the pending Hydro step and extracts its immutable current/previous
// velocity snapshots for the public interface service. Copies
// only nozzle IDs/initial-contact diagnostics. The App's coupled time loop,
// not this adapter or Interface, decides whether the result is accepted.
ElectrosprayInterfaceMotionCandidate prepareElectrosprayInterfaceMotion(
    const QuadtreeMesh& mesh, const HydroFittedDomain& domain,
    const ElectrosprayInterfaceState& previous, const Hydrodynamics& hydro,
    const InterfaceStepTime& time, double capillaryHalfWidthCells,
    double linearTolerance, const InterfaceVolumeProgress& progress = {},
    const InterfaceGeometryObserver& observeGeometry = {});
