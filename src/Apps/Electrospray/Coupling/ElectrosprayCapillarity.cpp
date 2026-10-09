#include "Electrospray/Coupling/ElectrosprayCapillarity.h"
#include "Modules/Interface/Compatibility/InterfaceState.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Units/UnitsContext.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

HydroFittedCapillaryState assembleElectrosprayCapillarity(const QuadtreeMesh& mesh,
    const HydroFittedDomain& domain, const InterfaceState& interfaceState,
    double referenceDensitySI, double surfaceTensionSI, double capillaryHalfWidthCells,
    const UnitsContext& units, double sourceTime, HydroForceEvolution evolution)
{
    if (!interfaceState.geometry() || !interfaceState.panels() || units.mode() != UnitsMode::kSI
        || !std::isfinite(capillaryHalfWidthCells) || capillaryHalfWidthCells <= 0.0)
        throw std::invalid_argument("spray capillarity requires active geometry, SI scales and positive capillaryHalfWidthCells");
    const double length = units.scales().L0_m;
    const double velocity = length / units.scales().T0_s;
    const double pressure = referenceDensitySI * velocity * velocity;
    if (!std::isfinite(pressure) || pressure <= 0.0 || !std::isfinite(length) || length <= 0.0)
        throw std::invalid_argument("invalid spray capillary reference scales");
    const bool native = interfaceState.geometry()->usesNativeLeaves();
    const auto spacing = native ? interfaceState.panels()->nativeChain().minimumActiveSpacing()
        : interfaceState.geometry()->cellSize();
    const double halfWidth = capillaryHalfWidthCells * std::max(spacing.x, spacing.y);
    const double sigma = surfaceTensionSI / (pressure * length);
    auto state = native
        ? buildNativeFittedCapillaryForce(mesh, domain, interfaceState.geometry(), *interfaceState.panels(), halfWidth, sigma)
        : buildFittedCapillaryForce(mesh, domain, *interfaceState.geometry(),
            interfaceState.pinnedCoordinates(), halfWidth, sigma);
    state.source = HydroForceSource::fitted(domain, sourceTime, evolution, interfaceState.geometry());
    return state;
}
