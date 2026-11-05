#pragma once
#include "Modules/Electric/ElectricVofSharpField.h"
#include "Modules/Interface/GeometricVofTransport.h"
#include <memory>
#include <stdexcept>

namespace LAMP::Electric
{
enum class ElectricSurfaceChargeTransferFailure {
    UnsupportedMesh, UnsupportedTopology, ForeignSource, InvalidClock,
    UnsupportedBoundary, ExcessiveDisplacement, ProjectionDistance,
    AmbiguousProjection, NonconvexProjection, InvalidMoments, Conservation
};
class ElectricSurfaceChargeTransferError : public std::runtime_error
{
public:
    ElectricSurfaceChargeTransferError(ElectricSurfaceChargeTransferFailure,const char*);
    ElectricSurfaceChargeTransferFailure reason() const noexcept { return reason_; }
private:
    ElectricSurfaceChargeTransferFailure reason_;
};
struct ElectricSurfaceChargeTransferConfig
{
    long double maximumDisplacementInCells=.45L;
    long double maximumProjectionDistanceInCells=.2L;
    long double relativeChargeTolerance=1e-10L;
};
enum class ElectricSurfaceChargeTransferMode { MaterialProjection, ExactPartitionIdentity };
struct ElectricSurfaceChargeTransferReport
{
    ElectricSurfaceChargeTransferMode mode=ElectricSurfaceChargeTransferMode::MaterialProjection;
    // Exact identity copies moments without measuring a density, charge L1,
    // spatial moments or source quadrature residual. Those fields are NaN.
    bool integralDiagnosticsAvailable=true;
    std::size_t identityArcTransfers=0,sourceDensityReconstructions=0;
    std::size_t sourceArcs=0,targetArcs=0,materialPanels=0,materialSamples=0;
    std::size_t parameterSplits=0;
    std::size_t nativeArcProjections=0,exactArcProjections=0,localChainCertificates=0;
    std::size_t rootIsolationCalls=0; // Completed Core calls; excludes exact endpoint identity.
    std::size_t endpointProjectionRequests=0,endpointProjectionMemoHits=0;
    // Passive classifications of completed Core projections only.
    std::size_t coreHessianFallbacks=0,coreEndpointSignFallbacks=0,coreNonfiniteFallbacks=0;
    std::size_t coreRootWidthFallbacks=0,coreForcedDistanceFallbacks=0;
    std::size_t coreEndpointProjections=0,coreInteriorProjections=0,geometryBasePreparations=0;
    // Passive nested wall times. Integration includes point projection;
    // point projection includes all arc paths. Never add these as disjoint costs.
    double prepareSeconds=0,materialIntegrationSeconds=0,pointProjectionSeconds=0;
    double nativeArcProjectionSeconds=0,coreArcProjectionSeconds=0,identityArcProjectionSeconds=0;
    long double sourceCharge=0,targetCharge=0,sourceAbsoluteCharge=0;
    long double chargeResidual=0,relativeChargeResidual=0;
    std::array<long double,2> sourcePositionMoment{},advectedPositionMoment{},projectedPositionMoment{};
    long double maximumProjectionDistance=0,maximumDisplacement=0;
    long double maximumTerminalParameterWidth=0,sourceMomentQuadratureResidual=0;
};
// Immutable material receipt. Every attempt starts from acceptedSource;
// no iteration can advance an earlier candidate. Exact parameter/partition
// identities copy all three signed moments without projection or quadrature;
// a physical stage also requires every actual native velocity sample zero
// (LocalPlic retains native normal means, not reconstructed vectors).
// This narrow direct path does not establish a general moving map.
// XY, same uniform mesh and
// one nonbranching component only. Unsupported transfers throw before publish.
class ElectricSurfaceChargeTransport
{
public:
    static ElectricSurfaceChargeTransport prepareStage(const QuadtreeMesh&,
        const ElectricVofSurfaceCharge& acceptedSource,
        const Interface::GeometricVofTransport&,
        const Interface::VofSmoothGeometry& target,
        const ElectricSurfaceChargeTransferConfig& = {});
    // Same actual VOF publication and complete End/Initial clock; only the
    // derived smooth owner is rebuilt. No AMR or topology inference.
    static ElectricSurfaceChargeTransport prepareSameEnd(const QuadtreeMesh&,
        const ElectricVofSurfaceCharge& acceptedSource,
        const Interface::VofSmoothGeometry& target,
        const ElectricSurfaceChargeTransferConfig& = {});
    const ElectricVofSurfaceCharge& sourceCharge() const;
    const ElectricVofSurfaceCharge& chargeCandidate() const;
    const Interface::GeometricVofTransport* interfaceTransport() const;
    const ElectricSurfaceChargeTransferReport& report() const;
    bool sameInstance(const ElectricSurfaceChargeTransport&) const noexcept;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
    explicit ElectricSurfaceChargeTransport(std::shared_ptr<const Data>);
    static ElectricSurfaceChargeTransport prepare(const QuadtreeMesh&,
        const ElectricVofSurfaceCharge&,const Interface::GeometricVofTransport*,
        const Interface::VofSmoothGeometry&,const ElectricSurfaceChargeTransferConfig&);
};
}


