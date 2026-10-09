#include "Electrospray/Coupling/ElectrosprayInterfaceMotion.h"
#include "Modules/Hydrodynamics/HydroInterfaceTransport.h"
#include "Modules/Hydrodynamics/HydroFittedVelocity.h"
#include "Modules/Hydrodynamics/Hydrodynamics.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

ElectrosprayInterfaceMotionCandidate prepareElectrosprayInterfaceMotion(
    const QuadtreeMesh& mesh, const HydroFittedDomain& domain, const ElectrosprayInterfaceState& previous,
    const Hydrodynamics& hydro, const InterfaceStepTime& time,
    double capillaryHalfWidthCells, double linearTolerance,
    const InterfaceVolumeProgress& progress, const InterfaceGeometryObserver& observeGeometry)
{
    const auto flow = hydro.state().progress;
    const double allowance = 64.0*std::numeric_limits<double>::epsilon()
        * std::max({std::abs(time.begin),std::abs(time.end),time.dt});
    if (flow.stopped || !flow.pendingAcceptance || flow.completedSteps != time.index
        || flow.acceptedSteps != previous.surface.progress().completedSteps
        || !std::isfinite(time.begin) || !std::isfinite(time.end)
        || std::abs(flow.acceptedEndTime-time.begin) > allowance
        || std::abs(flow.endTime-time.end) > allowance)
        throw std::invalid_argument("[ElectrosprayInterfaceMotion] velocity layers require the matching pending Hydro step");
    if (previous.surface.initialContactJets().size() != 1)
        throw std::invalid_argument("[ElectrosprayInterfaceMotion] this moving binding requires one retained pinned contact");
    const auto newVelocity = std::make_shared<const HydroFittedVelocitySampler>(
        hydro.fittedVelocity(mesh, domain, HydroVelocityTimeLevel::Current,
            HydroFittedVelocitySampler::Reconstruction::RZFluxCompatible));
    const auto oldVelocity = std::make_shared<const HydroFittedVelocitySampler>(
        hydro.fittedVelocity(mesh, domain, HydroVelocityTimeLevel::Previous,
            HydroFittedVelocitySampler::Reconstruction::RZFluxCompatible));
    InterfaceMotionVelocityLayers velocities;
    velocities.current=[newVelocity](const Vector2D& point){return newVelocity->sample(point);};
    velocities.previous=[oldVelocity](const Vector2D& point){return oldVelocity->sample(point);};
    velocities.independentQueries=true;
    ElectrosprayInterfaceMotionCandidate result(mesh);
    result.motion=InterfaceMotion::preparePinnedToAxis(mesh,HydroInterfaceTransport::domainView(domain),
        previous.surface,time,velocities,HydroInterfaceTransport::normalFlux(*newVelocity),
        capillaryHalfWidthCells,linearTolerance,result.state.surface,progress,observeGeometry);
    result.state.freeLineIds=previous.freeLineIds;
    result.state.fixedLineIds=previous.fixedLineIds;
    result.state.sourceLevelSetIndex=previous.sourceLevelSetIndex;
    return result;
}
