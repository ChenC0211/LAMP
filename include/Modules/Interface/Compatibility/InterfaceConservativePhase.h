#pragma once
#include "Modules/Interface/Compatibility/InterfaceState.h"
#include "Core/Geometry/QuadraticCellCut.h"
#include <memory>
#include <functional>
#include <cstddef>
#include <limits>
#include <optional>
#include <vector>

class NativeLeafDomain;
class NativeMeshIdentity;

// Immutable, independently transported liquid inventory. This is not the
// volume obtained by integrating the latest level-set zero surface. Creation
// and transport produce values; neither operation changes the input state.
// Initial geometry seeding remains uniform-leaf; explicit zero-time inventory
// remapping supports native parent/child leaves in XY/RZ. Geometric sweeping,
// scalar synchronization and full InterfaceState publication are separate.
class InterfaceConservativePhase
{
public:
    struct DirectionalFlux
    {
        // Nonnegative GROSS volumes: simultaneous spatial counterflow and
        // accumulated temporal reversal must not disappear into signed net F.
        long double totalPositive=0,totalNegative=0,liquidPositive=0,liquidNegative=0;
    };
    struct FaceFlux
    {
        long double totalVolume=0, highOrderLiquidVolume=0;
        // Required for any GROSS boundary inflow, including zero/outward net
        // volume. Unused values may be NaN.
        long double boundaryLiquidFraction=std::numeric_limits<long double>::quiet_NaN();
        // Optional spatial counterflow within this interval. Net fields must
        // equal these directional differences; fromDirectional establishes
        // that relation. Each liquid direction is bounded by its own total.
        std::optional<DirectionalFlux> directional;
        static FaceFlux fromDirectional(DirectionalFlux flux,
            long double inflowFraction=std::numeric_limits<long double>::quiet_NaN())
        {
            FaceFlux result;result.totalVolume=flux.totalPositive-flux.totalNegative;
            result.highOrderLiquidVolume=flux.liquidPositive-flux.liquidNegative;
            result.boundaryLiquidFraction=inflowFraction;result.directional=flux;return result;
        }
    };
    struct Report
    {
        long double liquidBefore=0, liquidAfter=0, boundaryOutward=0, closure=0;
        long double maxDivergence=0, maxOutflowCfl=0;
        std::size_t limitedFaces=0, intervals=0, limiterSweeps=0, limiterFallbacks=0;
        // Roundoff-only repair of the shared low-order interior phase flux.
        // No cell mass is clipped; boundary phase fluxes remain prescribed.
        std::size_t baselineRepairSweeps=0,baselineRepairedFaces=0;
        long double maxBaselineFluxCorrection=0;
    };
    struct Interval { double begin=0,end=0,dt=0; };
    using FluxProvider=std::function<std::vector<FaceFlux>(std::size_t,const Interval&,
        const InterfaceConservativePhase&)>;
    // Native cell indexing, inactive entries must be NaN. Initial volumes
    // must already be bounded; no geometric repair or hidden mass clipping.
    static InterfaceConservativePhase initialize(const QuadtreeMesh& mesh,
        const NativeLeafDomain& domain, const std::vector<long double>& liquidVolumes);
    struct MeshTransferReport
    {
        std::size_t refinedSourceCells=0,coarsenedTargetCells=0,uniformTailSplits=0;
        std::size_t geometricSourceCells=0,geometricFallbackCells=0;
        long double geometricCorrectionL1=0;
        long double liquidBefore=0,liquidAfter=0,closure=0;
        long double maxPartitionCorrection=0,maxSourceClosure=0;
        long double pureCapacityAdjustment=0; // Signed native-measure roundoff on coarsening.
    };
    // Zero-time native parent/child remap on the SAME background lattice and
    // fixed active domain. Resolved split cells use source-C-matched PLIC or
    // quadratic donors; unsplit cells retain their inventory exactly. Pure
    // children use native capacities, unresolved minorities retain their cell
    // average. Only roundoff-size partition residuals may be redistributed
    // within child bounds; any remaining measure-roundoff closure is reported.
    // Normals/shapes belong to this immutable source inventory.
    // Optional target geometric volumes guide ONLY new subdivision degrees
    // of freedom: a bounded weighted projection retains EACH source parent's
    // original C. It never reseeds unsplit/pure cells or overwrites old C.
    // Geometrically pure children are fixed even during roundoff closure;
    // residuals without an admissible mixed recipient remain reported.
    // If mixed children cannot carry
    // the parent's inventory, retain/report the existing source donor fallback.
    // This transfers only inventory, not geometry, momentum, charge, or App
    // acceptance. Mixed-leaf transport/synchronization remain separate gates.
    static InterfaceConservativePhase transferMesh(const QuadtreeMesh& sourceMesh,
        const NativeLeafDomain& sourceDomain,const InterfaceConservativePhase& source,
        const QuadtreeMesh& targetMesh,const NativeLeafDomain& targetDomain,
        const std::vector<Vector2D>& sourceNormals,
        const std::vector<std::optional<Quantum::Geometry::QuadraticCellCut::Shape>>& sourceCurves={},
        const std::vector<long double>& targetGeometricVolumes={});
    bool transferredFrom(const InterfaceConservativePhase& source) const noexcept;
    const std::optional<MeshTransferReport>& meshTransferReport() const noexcept {return meshTransferReport_;}
    // A single time interval: face orientation is +x/+y, shared by
    // both owners; inputs are time INTEGRATED volumes, not end-time rates.
    // Stage-wise temporal reversal reconstructs geometry at each stage. A
    // resolved full-step space-time flux may instead retain ONE old material
    // representation and supply both gross directions in FaceFlux::directional;
    // the same gross-outflow bound still applies. Spatial counterflow uses it too. High order
    // fluxes may be omitted (empty) for zero motion,
    // otherwise every native face has one entry. No wall physics is inferred.
    static InterfaceConservativePhase advance(const QuadtreeMesh& mesh,
        const NativeLeafDomain& domain, const InterfaceConservativePhase& previous,
        const InterfaceStepTime& time, const std::vector<FaceFlux>& fluxes);
    // One physical clock update over strictly increasing normalized interval
    // ends in (0,1], the last exactly 1. Reconstruct each interval's geometry
    // from the supplied intermediate inventory. Its progress remains the old
    // physical clock until the WHOLE candidate succeeds. Callbacks are borrowed
    // synchronously; no partial result is published on callback/interval failure.
    // Each interval still requires monotone flow and its own outflow CFL bound.
    static InterfaceConservativePhase advanceSubsteps(const QuadtreeMesh& mesh,
        const NativeLeafDomain& domain,const InterfaceConservativePhase& previous,
        const InterfaceStepTime& time,const std::vector<double>& intervalEnds,
        const FluxProvider& provider);
    void validate(const QuadtreeMesh& mesh, const NativeLeafDomain& domain) const;
    bool follows(const InterfaceConservativePhase& previous,const InterfaceStepTime& time) const noexcept;
    const std::vector<long double>& liquidVolumes() const noexcept { return liquid_; }
    const std::vector<long double>& capacities() const noexcept { return capacity_; }
    const std::vector<long double>& totalFluxes() const noexcept { return totalFlux_; }
    const std::vector<long double>& liquidFluxes() const noexcept { return liquidFlux_; }
    const std::vector<long double>& limiter() const noexcept { return limiter_; }
    // For substeps limiter() is the minimum theta per face; it is not an
    // equivalent limiter for the net flux. limitedFaces counts face-intervals.
    const std::vector<DirectionalFlux>& directionalFluxes() const noexcept { return directionalFlux_; }
    InterfaceStepProgress progress() const noexcept { return progress_; }
    const Report& report() const noexcept { return report_; }
private:
    InterfaceConservativePhase() = default;
    struct Revision {};
    std::shared_ptr<const Revision> revision_,parentRevision_;
    std::shared_ptr<const Revision> meshSourceRevision_;
    std::optional<MeshTransferReport> meshTransferReport_;
    static InterfaceConservativePhase initializeNative(const QuadtreeMesh&,const NativeLeafDomain&,
        const std::vector<long double>&);
    static InterfaceConservativePhase advanceInterval(const QuadtreeMesh& mesh,
        const NativeLeafDomain& domain,const InterfaceConservativePhase& previous,
        const std::vector<FaceFlux>& fluxes);
    std::shared_ptr<const NativeMeshIdentity> identity_;
    bool axisymmetric_=false;
    std::vector<int> cells_;
    std::vector<long double> capacity_,liquid_,totalFlux_,liquidFlux_,limiter_;
    std::vector<DirectionalFlux> directionalFlux_;
    InterfaceStepProgress progress_;
    InterfaceStepTime stepTime_;
    Report report_;
};
