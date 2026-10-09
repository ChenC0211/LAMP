#include "Modules/Hydrodynamics/Hydrodynamics.h"
#include "Modules/Hydrodynamics/HydroBoundaryVelocity.h"
#include "Modules/Hydrodynamics/HydroInterfaceBoundary.h"
#include "Core/Mesh/NativeMeshIdentity.h"

#include <stdexcept>

void Hydrodynamics::validateCompleteVelocityQuery(const QuadtreeMesh& mesh,
    HydroVelocityTimeLevel timeLevel) const
{
    if (!initialStateSet_ || stepProgress_.stopped || !stepProgress_.geometry
        || *stepProgress_.geometry == HydroDomainKind::FittedRZ || !stepMeshIdentity_)
        throw std::logic_error("complete velocity query requires initialized, non-stopped complete-domain Hydro");
    if (timeLevel != HydroVelocityTimeLevel::Current && timeLevel != HydroVelocityTimeLevel::Previous)
        throw std::invalid_argument("unknown Hydro velocity time layer");
    if (timeLevel == HydroVelocityTimeLevel::Previous && stepProgress_.completedSteps == 0)
        throw std::logic_error("previous Hydro velocity requires a completed substep");
    stepMeshIdentity_->validate(mesh);
}

std::pair<std::unique_ptr<VertField<double>>, std::unique_ptr<VertField<double>>>
Hydrodynamics::vertexVelocity(const QuadtreeMesh& mesh, HydroVelocityTimeLevel timeLevel,
    bool explicitBoundaryPatches, bool radialAxisSymmetry) const
{
    validateCompleteVelocityQuery(mesh, timeLevel);
    const bool previous = timeLevel == HydroVelocityTimeLevel::Previous;
    return hydroVertexVelocity(mesh, previous ? uxPre : ux, previous ? uyPre : uy,
        explicitBoundaryPatches, radialAxisSymmetry);
}

std::function<Vector2D(Vector2D)> Hydrodynamics::boundaryAdvectionVelocity(const QuadtreeMesh& mesh,
    HydroVelocityTimeLevel timeLevel, const std::vector<HydroBoundaryPatch>& patches,
    bool radialAxisSymmetry) const
{
    validateCompleteVelocityQuery(mesh, timeLevel);
    const bool previous = timeLevel == HydroVelocityTimeLevel::Previous;
    return hydroBoundaryAdvectionVelocity(mesh, previous ? uxPre : ux, previous ? uyPre : uy,
        patches, radialAxisSymmetry);
}

StaggeredVelocitySampler Hydrodynamics::averagedNativeVelocity(const QuadtreeMesh& mesh) const
{
    validateCompleteVelocityQuery(mesh, HydroVelocityTimeLevel::Previous);
    return StaggeredVelocitySampler(mesh, uxPre, uyPre, ux, uy);
}

HydroSplineVelocity Hydrodynamics::splineVelocity(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    HydroVelocityTimeLevel timeLevel) const
{
    validateCompleteVelocityQuery(mesh,timeLevel);
    if(*stepProgress_.geometry!=HydroDomainKind::CompleteXY || !frozenBoundaryPatches_.empty())
        throw std::invalid_argument("C1 spline snapshot currently requires closed complete XY Hydro");
    const bool previous=timeLevel==HydroVelocityTimeLevel::Previous;
    return HydroSplineVelocity(mesh,domain,previous?uxPre:ux,previous?uyPre:uy);
}

double Hydrodynamics::completePhaseOutwardFlux(const QuadtreeMesh& mesh,
    const LevelSetHermiteGeometry::Snapshot& interfaceGeometry, bool axisymmetric) const
{
    validateCompleteVelocityQuery(mesh, HydroVelocityTimeLevel::Current);
    if (axisymmetric != (*stepProgress_.geometry == HydroDomainKind::CompleteRZ))
        throw std::invalid_argument("phase flux coordinate measure differs from Hydro state");
    return HydroInterfaceBoundarySupport::negativePhaseOutflow(mesh, interfaceGeometry, ux, uy, axisymmetric);
}
