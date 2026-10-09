#pragma once
#include "Modules/Interface/Compatibility/InterfaceConservativePhase.h"
#include "Modules/Interface/Compatibility/InterfaceCurvature.h"
#include "Modules/Interface/Compatibility/InterfacePlicTransport.h"
#include "Core/Base/Enumeration.h"

struct InterfaceMotionVelocityLayers;
class InterfaceSurfaceTransfer;
struct InterfaceInteriorMaterialTransfer;

struct InterfaceConservativeGeometryReport
{
    // The predictor and published scalar are different roles. This retained
    // source identifies the numerical synchronization; it is NOT a second
    // physical interface. An explicitly requested material plan binds the
    // physical trajectory to the final synchronized geometry.
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> predictor;
    std::shared_ptr<const InterfaceSurfaceTransfer> surfaceTransfer;
    std::shared_ptr<const InterfaceInteriorMaterialTransfer> interiorMaterialTransfer;
    std::vector<long double> geometricVolumeByCell;
    long double inventoryVolume=0,geometricVolume=0,signedDifference=0,l1Difference=0;
    long double maxFractionDifference=0;
    double distanceBandWidth=0,maxNodalScalarChange=0;
    // Zero-time-only, bidirectional endpoint/midpoint projection observation.
    // Each sample must stay within h/4 and retain normal orientation. This
    // finite native search is not a continuous Hausdorff certificate.
    double maxZeroTimeInterfaceShift=0;
    std::size_t plicSegments=0,ambiguousDistanceNodes=0;
    std::size_t fixedWallTraceVertices=0;
    // Section-detachment model handoff: once two local axial cell widths
    // adjacent to the old section are geometrically pure gas according to
    // independently transported C, ordinary inventory reconstruction can no
    // longer reconnect that section. No inventory is erased or reseeded.
    bool terminalInventoryHandoff=false;
    std::size_t terminalGasStripCells=0;
    bool interiorNeckSplit=false,interiorGapProtected=false;
    // Post-interior-event local C1 scalar correction. Match already cut
    // cells to the independent inventory while preserving pure-cell phase,
    // admissible monotone charts and the two-child topology. Geometric gas
    // inventories outside those cuts remain separate, unchanged diagnostics.
    int eventInventoryIterations=0,eventInventoryBacktracks=0;
    std::size_t eventInventoryCutCells=0,eventInventorySupportCells=0;
    double eventInventoryInitialCutError=0,eventInventoryFinalCutError=0;
    double eventInventorySupportError=0,eventInventoryPureTailError=0;
    // Implicit bounded concentration transport on accepted directional liquid
    // fluxes. Residual is the largest absolute global child mass change.
    int componentTransportIterations=0;
    std::size_t componentImplicitCells=0;
    long double componentMassResidual=0,componentRowResidual=0;
    // Counts in the distance band, not interface-measure coverage. Height jets
    // match column moments and carry consistent mixed derivatives; missing
    // supports retain the older PLIC-distance baseline, whose accuracy remains
    // visible in the independent geometric/inventory discrepancy.
    std::size_t heightReconstructedNodes=0,distanceFallbackNodes=0;
    std::size_t blendedHeightNodes=0;
    double maxHeightChartSpread=0; // Accepted chart scalar spread divided by h.
    // Connected PLIC segment-moment graph fit when volume columns are unavailable.
    // This is a pre-publication geometry fit, distinct from curvature queries.
    std::size_t localGraphNodes=0;
    // Nodes whose compact blend includes a cached local PLIC graph. This is
    // a subset of localGraphNodes; fit diagnostics cover these accepted blends.
    // Fit residuals use max(hx,hy); chart spread uses min(hx,hy). RMS is the
    // length-weighted segment-moment residual, also after local cell-volume
    // refinement; it is not the nonlinear volume objective. The maximum
    // midpoint gap remains separate even when a tiny facet is an outlier.
    std::size_t blendedLocalGraphNodes=0;
    std::size_t cubicLocalGraphNodes=0; // Blends containing a volume-refined cubic chart.
    double maxLocalGraphRms=0,maxLocalGraphResidual=0,maxLocalGraphChartSpread=0;
    // Band vertices with conflicting incident-cell PLIC phase traces, and
    // proposed height/local fits conflicting with incident pure C (including
    // those subsequently recovered by the constrained scalar projection).
    // No inventory fraction is clipped or corrected.
    std::size_t localPhaseDisagreementNodes=0,purePhaseRejectedFits=0;
    // Phase-conflicting fits recovered by the nearest pure-cell Bernstein
    // scalar bound, keeping their derivatives. The corrected scalar cannot
    // be farther into the pure phase than the existing PLIC-distance fallback.
    std::size_t phaseCorrectedFits=0;
    double maxPhaseScalarCorrection=0;
    // Shared derivative scaling makes pure-cell Bernstein controls retain
    // their inventory phase. Nodal values and conservative inventory stay fixed.
    std::size_t pureLimitedVertices=0;
    // Mixed inventories below the existing PLIC inverse resolution use the
    // same geometric phase in column fits, node signs and pure-cell controls.
    // Their full inventory and geometric discrepancy remain in the reports.
    std::size_t unresolvedPhaseConstrainedCells=0;
    double minPureJetScale=1;
    // Sparse axis trace diagnostics before/after the two shared-jet limiters.
    // Values stay fixed. This separates column/blend position error from
    // root movement caused by limiting the axial Hermite derivative.
    struct AxisJet {int vertex=-1;double value=0,dyBefore=0,dyMixed=0,dyPure=0;};
    std::vector<AxisJet> axisJets;
    // Band-node provenance without changing reconstruction/limiter decisions.
    // source: 0 ridge, 1 height blend, 2 local-graph blend, 3 single height,
    // 4 single local graph, 5 PLIC distance, 6 final event-volume correction.
    // For source 6, value/pure contain the final published nodal data.
    // Rejected signs are reported in
    // attempt order (blend, height, local graph); NaN means no phase rejection.
    struct NodalJet {
        int vertex=-1,source=0;
        double value=0;
        std::array<double,3> before{},mixed{},pure{};
        std::array<double,3> rejectedValues{};
        double scalarCorrection=0;
    };
    std::vector<NodalJet> nodalJets;
    // Single-graph mixed cells with monotone corner data constrain the cubic
    // Bernstein derivative along their PLIC normal axis. Shared jet scaling
    // preserves values/inventory; cells without monotone corner data are
    // reported and still go through the original geometry validation.
    std::size_t monotoneLimitedVertices=0,nonmonotoneNodalCells=0;
    double minMonotoneJetScale=1;
    struct MonotoneCell {int cell=-1,axis=0,sign=1;};
    std::vector<MonotoneCell> monotoneCells;
    // Mixed cells whose PLIC directions come from the same admissible column
    // moments. Remaining cells retain predictor directions; offsets are always
    // re-solved from unchanged C. Normal change is the norm of the unit-vector difference.
    std::size_t heightPlicCells=0;
    struct PlicDirection {int cell=-1;Vector2D normal;};
    // Resolved synchronization planes after normal selection and exact C
    // offset inversion. Physical donors use the same reconstruction service
    // on the old state; these records belong to the synchronized target.
    std::vector<PlicDirection> plicDirections;
    // Donors used in the physical step, not a second published surface.
    // Curves preserve each original cell volume through their solved offset.
    std::size_t curvedTransportCells=0;
    // TransportDonors::NormalFitStatus order; zero for the older models.
    std::array<std::size_t,6> normalFrameDonorStatus{};
    long double maxNormalFrameDonorRms=0;
    std::size_t regularizedDonorCells=0;
    double maxPlicNormalChange=0;
    // Mixed inventories with minority volume below the PLIC inverse's own
    // absolute volume precision. They do not determine a reliable geometric
    // facet, but are NOT removed or rounded in C. No degenerate
    // segment is used in the node fit; the associated raw volume defect is
    // explicit and the final C_Gamma-C still accounts for every cell. An
    // entirely unresolved interface remains unsupported, never "disappeared".
    std::vector<int> unresolvedPlicCells;
    long double unresolvedPlicVolumeDefect=0;
    // Native cell of the selected closest PLIC piece; -1 outside the band.
    // This is a local donor label, NOT a persistent physical branch ID.
    std::vector<int> closestPlicCellByVertex;
    struct ComponentDifference
    {
        std::uint64_t liquidId=0;
        long double inventoryVolume=0,geometricVolume=0,signedDifference=0,l1Difference=0;
    };
    std::vector<ComponentDifference> components;
    // Component sums use the published geometric regions. Positive inventory
    // in cells not represented by any such region remains explicit here; it
    // must not be dropped or arbitrarily assigned to a nearby liquid ID.
    long double unassignedInventoryVolume=0;
    struct SurfaceCurvature
    {
        std::uint64_t surfaceId=0;
        std::vector<InterfaceCurvature::Sample> graphMidpoints,axisEnds,symmetryEnds;
    };
    // Same published geometry, with primary-height status retained even when
    // the same-arc fit succeeds. Counts include axis and symmetry endpoints.
    // Unavailable values stay nonfinite. Actual force faces use extend(), not
    // interpolation of these midpoint diagnostics; contact/coupled accuracy
    // still needs its own acceptance evidence.
    std::vector<SurfaceCurvature> curvature;
    std::size_t validHeightSamples=0,invalidHeightSamples=0;
    std::size_t validLocalFitSamples=0,unavailableCurvatureSamples=0;
};

class InterfaceConservativeGeometry
{
public:
    // Midpoint/endpoint report samples only; actual force providers always
    // construct and validate their own curvature. Omit when the consumer
    // replaces the provider and would immediately discard these samples.
    enum class CurvatureDiagnostics { Record, Omit };
    struct BoundaryInput
    {
        const InterfaceTransportDomain& domain;
        // Immutable whole-step phase: -1 liquid, +1 gas on SinglePhasePort;
        // zero on all other faces. A pure old/new boundary layer is required.
        InterfacePlicTransport::SinglePhasePorts portPhase;
        // Scalar at the resolved inflow event, using ABSOLUTE physical time.
        // Its finite nonzero sign must match the same declared phase history.
        std::function<double(Boundary,Vector2D,double)> inflowScalar;
    };
    struct AdvanceReport
    {
        InterfacePlicTransport::TraceReport traces;
        InterfacePlicTransport::SweepReport sweeps;
        InterfacePlicTransport::SpaceTimeReport spaceTime;
        bool usedSpaceTime=false;
        std::size_t heightPlicCells=0; // Transport donor directions, before this step.
        int predictedNodes=0;
        double maxTraceRatio=0;
        // Wall-clock observations, collected only when an AdvanceReport is
        // requested: donors, phase fluxes, inventory, predictor, synchronization,
        // optional material transfer. Never used in a numerical decision.
        std::array<double,6> stageSeconds{};
    };
    enum class DonorModel { QuadraticMoments, ProjectedMoments, NormalFrameVolumeFit };
    static constexpr DonorModel defaultDonorModel=DonorModel::NormalFrameVolumeFit;
    struct TransportDonors
    {
        enum class NormalFitStatus { Fitted, MissingBranch, BoundarySupport,
            AmbiguousBranch, CurvatureUnavailable, OrientationBoundary };
        struct NormalFit {
            int cell=-1;
            NormalFitStatus status=NormalFitStatus::MissingBranch;
            double meridional=0,angleChange=0;
            long double initialRms=0,finalRms=0,inverseResidual=0;
            bool regularizedCurvature=false;
        };
        // Default normal-frame parabolas use valid branch meridional curvature
        // and fit direction to the actual native 3x3 phase fractions. The
        // central original C fixes the offset in every trial. A minimum at
        // the +/-30 degree search boundary is rejected. Missing/ambiguous
        // branch, curvature or physical-boundary support retains the previous
        // projected column-moment graph (then three moments/PLIC); axis parity
        // therefore retains its existing fit. Mixed-level gaps still request
        // the existing same-level support refinement before this fit.
        // Publication/curvature keep their higher-order stencils; no donor
        // replaces the one published interface or its conservative inventory.
        std::vector<Vector2D> normals;
        InterfacePlicTransport::CurvedDonors curved;
        std::size_t heightCells=0;
        // Populated only for NormalFrameVolumeFit, one record per resolved
        // mixed donor. Unsupported support explicitly retains ProjectedMoments.
        std::vector<NormalFit> normalFits;
    };
    // Immutable reconstruction inputs for the existing PLIC flux service.
    // Offsets are still inverted from the source inventory by that service;
    // these are donor shapes, not a second published physical interface.
    // Mixed meshes require the union of both five-column/radius-six stencils
    // around every resolved mixed donor to consist of actual same-level leaves.
    // An AMR gap raises InterfaceGeometryResolution::MixedResolution instead
    // of silently reducing reconstruction order. Physical boundaries retain
    // their existing laws; selected/contact domains still require their own
    // boundary capability. C is never virtually refined for a volume column.
    static TransportDonors reconstructTransportDonors(const QuadtreeMesh&,const NativeLeafDomain&,
        const InterfaceState&,const BoundaryInput* boundaries=nullptr,
        DonorModel model=defaultDonorModel);
    // Attempt one physical step using the existing immutable velocity and flux
    // layers: space-time PLIC/quadratic donor flux -> bounded inventory -> native jet
    // predictor -> synchronization. Fixed native fluid cells, including a
    // mixed mesh with the common-level interface support and scalar plateau
    // required by prepare(). Mixed meshes require the space-time flux path.
    // Selected domains require explicit boundary laws and space-time flux.
    // Mixed leaves admit the radial axis, separated outer slip layers and
    // pure outer ports, fixed wall/contact traces and retained wet corner
    // sectors with a locally same-level corner; common-level band rules apply.
    // Boundaries admit fixed no-slip contacts and explicit single-phase ports.
    // Explicit outer reflection sections admit moving orthogonal endpoints
    // and even column stencils (RZ: axial sections only). Other slip faces
    // away from the radial axis require an unchanged
    // geometrically pure adjacent phase layer; tangential motion is free,
    // contact events reject. Below-inverse-resolution minority is retained
    // under the same geometric classification as the volume reconstruction;
    // the impermeable boundary still has exactly zero gross phase/total flux.
    // Ports continue to require an exactly pure inventory and declared history.
    // Retained flat wet-corner sectors use explicit cap/wet chart binding.
    // A supplied NeckEventPolicy admits the source-bound RZ four-cell
    // interior split before identity matching. Its material closure is not
    // yet connected, so this path requires retainSurfaceTransfer=false.
    // Other topology/contact-event reconstruction remains unsupported.
    // Ports require the space-time path and matching scalar inflow history.
    // Unresolved/folded sweeps still reject. Both paths use physical
    // two-layer RK2, with no scalar-amplitude taper. A successful state is a
    // module candidate, not NS/EHD acceptance. The optional material plan is
    // constructed eagerly from this exact predictor's forward map; it retains
    // no velocity callbacks. Failure rejects the entire candidate. Consumers
    // still own and advance their integrated surface quantities.
    // Polynomial views and shared spatial views with separable amplitudes use
    // space-time flux by default. Explicit cuts or other views retain the
    // earlier straight-sweep contract.
    // Failed candidates/reports do
    // not publish partial physical data; the old state remains untouched.
    // parallelFlux promises concurrent-safe normal-flux queries; the velocity
    // layers must also explicitly promise independentQueries. Opt-in only.
    struct NeckEventPolicy {double radiusCells=0;int maximumLevel=0;bool retainMaterial=false;};
    static void advance(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
        const InterfaceState& previous,const InterfaceStepTime& time,
        const InterfaceMotionVelocityLayers& velocities,const InterfaceNormalFluxView& previousFlux,
        const InterfaceNormalFluxView& currentFlux,InterfaceState& destination,
        double distanceBandCells=4,const std::vector<std::vector<double>>& partitions={},AdvanceReport* report=nullptr,
        bool retainSurfaceTransfer=false,const BoundaryInput* boundaries=nullptr,bool parallelFlux=false,
        CurvatureDiagnostics curvatureDiagnostics=CurvatureDiagnostics::Record,
        const NeckEventPolicy* neckEvents=nullptr);
    // CLSVOF synchronization on fixed native fluid cells,
    // with uniform leaves or one common same-level interface band. For mixed
    // leaves all resolved donors need their complete same-level support, and
    // every coarse/fine transition must have a constant signed scalar plateau
    // outside the reconstruction band. Native hanging payload remains dependent.
    // Resolved components (RZ axis permitted); an explicit NeckEventPolicy
    // also permits the conservative interior-neck event. Contact AMR
    // additionally requires the boundary-aware InterfaceMeshTransfer contract.
    // An explicit
    // boundary input (required for selected domains) admits stationary impermeable walls: retain
    // the previous cubic wall trace, including declared endpoint locations,
    // while its normal derivative/contact angle remains free. Single-phase
    // ports must stay separated from the interface by their pure donor layer.
    // Predictor supplies cell normals; C supplies
    // offsets. Admissible column-moment graph jets, with PLIC-distance fallback,
    // form a shared nodal cubic Hermite scalar
    // -> existing topology/queries, with C_Gamma-C measured separately.
    // Finite-distance PLIC endpoints and unconnected face traces are not used
    // directly as downstream surfaces. The once-published scalar is authority.
    // The complete candidate is published only after geometry and diagnostics
    // succeed. Inventory lineage and physical time must match previous exactly.
    static void prepare(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
        const InterfaceState& previous,InterfaceConservativePhase inventory,
        std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> predictor,
        const InterfaceStepTime& time,InterfaceState& destination,double distanceBandCells=4,
        const BoundaryInput* boundaries=nullptr,CurvatureDiagnostics curvatureDiagnostics=CurvatureDiagnostics::Record,
        const NeckEventPolicy* neckEvents=nullptr);
    // Zero-time inventory-driven reconstruction, including after AMR. Copies
    // the exact immutable inventory (and its remap provenance), keeps the
    // physical clock, IDs, event history and ancestry. Reuses the physical
    // step's scalar reconstruction and pure-phase/quality gates. A changed
    // topology or ambiguous overlap rejects the whole fresh destination.
    // This creates no material/charge remap and no artificial zero-flux step.
    // Fixed contact endpoints require explicit boundary laws; pure phase is
    // checked without inventing a physical-step face-flux ledger after remap.
    static void resynchronize(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
        const InterfaceState& previous,InterfaceState& destination,double distanceBandCells=4,
        const BoundaryInput* boundaries=nullptr,CurvatureDiagnostics curvatureDiagnostics=CurvatureDiagnostics::Record);
    // End-to-end constant-velocity step for closed interfaces in a full uniform
    // rectangle, with pure-gas exterior. Exact translation jet pullback feeds
    // the same conservative donor flux/synchronization services; no manufactured end geometry
    // or geometric-volume overwrite. RZ velocity must be axial. This entry does
    // not stand in for the general velocity predictor. Optional material
    // transfer uses this exact translation and the published target geometry.
    static void translate(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
        const InterfaceState& previous,const InterfaceStepTime& time,Vector2D velocity,
        InterfaceState& destination,double distanceBandCells=4,bool retainSurfaceTransfer=false);
private:
    static void correctEventInventory(const QuadtreeMesh&,const NativeLeafDomain&,
        const InterfaceState&,const InterfaceConservativePhase&,
        Quantum::Geometry::ActiveHermitePatch&,InterfaceConservativeGeometryReport&);
    static void prepareImpl(const QuadtreeMesh&,const NativeLeafDomain&,const InterfaceState&,
        InterfaceConservativePhase,std::shared_ptr<const LevelSetHermiteGeometry::Snapshot>,
        const InterfaceStepTime*,InterfaceState&,double,const BoundaryInput*,CurvatureDiagnostics,const NeckEventPolicy*);
};
