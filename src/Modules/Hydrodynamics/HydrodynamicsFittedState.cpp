#include "Modules/Hydrodynamics/Hydrodynamics.h"
#include "Modules/Hydrodynamics/HydroInterfaceTransport.h"

#include <stdexcept>

HydroFittedVelocitySampler Hydrodynamics::fittedVelocity(const QuadtreeMesh& mesh,
    const HydroFittedDomain& domain, HydroVelocityTimeLevel timeLevel,
    HydroFittedVelocitySampler::Reconstruction reconstruction) const
{
    if (stepProgress_.stopped)
        throw std::logic_error("cannot sample a stopped fitted Hydro state");
    if (timeLevel != HydroVelocityTimeLevel::Current && timeLevel != HydroVelocityTimeLevel::Previous)
        throw std::invalid_argument("unknown Hydro velocity time layer");
    if (timeLevel == HydroVelocityTimeLevel::Previous && fixedDt_ <= 0.0)
        throw std::invalid_argument("previous Hydro velocity requires an accepted fixed-dt step");
    const bool previous = timeLevel == HydroVelocityTimeLevel::Previous;
    return HydroFittedVelocitySampler(mesh, domain, previous ? uxPre : ux, previous ? uyPre : uy, reconstruction);
}

long double Hydrodynamics::fittedPhaseOutwardFlux(const QuadtreeMesh& mesh, const HydroFittedDomain& domain,
    const LevelSetHermiteGeometry::Snapshot& interfaceGeometry) const
{
    if (stepProgress_.stopped)
        throw std::logic_error("cannot query phase flux from a stopped fitted Hydro state");
    return HydroInterfaceTransport::outwardPhaseFlux(mesh, domain, interfaceGeometry, ux, uy);
}
