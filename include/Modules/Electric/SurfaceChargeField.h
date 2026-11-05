#pragma once

#include <array>
#include <cstdint>
#include <functional>
#include <map>
#include <memory>
#include <mutex>
#include <utility>
#include <vector>
#include "Core/Geometry/LevelSetInterfacePanels.h"
#include "Modules/Electric/SurfaceChargeSource.h"
#include "Modules/Electric/SurfaceChargeTransport.h"

class StaggeredVelocitySampler;

namespace LAMP::Electric
{
class SurfaceChargeField;
// Synchronous borrow of a completed geometric transfer, not a trajectory.
// The orchestration layer validates its geometry-owner provenance. Electric
// independently binds source panels, target panels and the physical clock.
struct SurfaceChargePreparedTransfer
{
    const LevelSetInterfacePanels* source=nullptr;
    const LevelSetInterfacePanels* target=nullptr;
    std::vector<std::vector<SurfaceChargeTransport::PanelWeight>> weights;
    int step=0;
    double sourceTime=0,targetTime=0,dt=0;
    double maxMaterialDisplacement=0,maxProjectionDistance=0;
};
// Separate zero-time closure after an ordinary physical material transfer.
// Only one upper RZ reflection endpoint may become an axis pole. The geometry
// owner supplies a conservative local map; remote native panels must be exact
// identities. No normal current, advection or physical clock is advanced here.
struct SurfaceChargeSectionDetachment
{
    const LevelSetInterfacePanels* source=nullptr;
    const LevelSetInterfacePanels* target=nullptr;
    std::vector<std::vector<SurfaceChargeTransport::PanelWeight>> weights;
    std::array<int,2> cells{{-1,-1}};
    int step=0;
    double time=0,beginZ=0,upperSection=0;
};
struct SurfaceChargeGeometryOptions
{
    bool axisymmetric = false;
    std::vector<Vector2D> boundaryContacts;
    CartesianGridScope scope = CartesianGridScope::RootOnly;
};

struct SurfaceChargeFieldStep
{
    SurfaceChargeStepReport transport;
    double unassignedBulkCurrent = 0.0;
    double absoluteBulkCurrent = 0.0;
    bool remapped = false;
    std::size_t previousPanelCount = 0;
    double remapConservationResidual = 0.0;
    double maxMaterialDisplacement = 0.0;
    double maxProjectionDistance = 0.0;
    bool sectionDetached = false;
    double sectionRemapConservationResidual = 0.0;
};

// Current charge and additional historical integral arrays remapped together.
// Every input history uses the old Field's panel order; every output history
// uses field->geometry() order. No charge source, time formula or stretching
// term is applied. Geometry/trajectories are prepared exactly once.
struct SurfaceChargeMaterialRemap
{
    std::unique_ptr<SurfaceChargeField> field;
    SurfaceChargeFieldStep report;
    std::vector<std::vector<double>> history;
    std::vector<SurfaceChargeStepReport> historyReports;
};

// Completed native material geometry/endpoint-overlap preparation. Only Field's
// real trajectory/projection factory can construct it. Owns both immutable
// panel publications and local physical coordinates; no Q or time history.
class SurfaceChargeMaterialTransfer
{
public:
    const LevelSetInterfacePanels& sourcePanels() const;
    const LevelSetInterfacePanels& targetPanels() const;
    const std::vector<SurfaceChargeTransport::PanelCoordinate>& materialBoundaries() const;
    const SurfaceChargeFieldStep& geometryReport() const;
    int step() const;
    double sourceTime() const;
    double targetTime() const;
    double dt() const;
private:
    struct Data;
    explicit SurfaceChargeMaterialTransfer(std::shared_ptr<const Data>);
    const Data& checked() const;
    std::shared_ptr<const Data> data_;
    friend class SurfaceChargeField;
};

// A zero-time mesh rebind with unchanged physical panels. The bijection copies
// INTEGRATED charge and every supplied time-history entry without arithmetic.
// No current, characteristic, overlap diffusion or time step is introduced.
struct SurfaceChargeMeshTransfer
{
    std::unique_ptr<SurfaceChargeField> field;
    std::vector<std::size_t> sourceToTargetPanel;
    std::vector<std::vector<double>> history;
};

// Read-only value description of the existing conservative density reconstruction.
// Columns correspond to neighbor Q integrals; rows are polynomial coefficients
// a,b,c in meridional u=(s-origin)/scale. Ordinary panels use a+u*(b+u*c);
// an axis-pole panel uses a+u^2*(b+u^2*c). Physical charge measures are already
// included in the weights: callers do not divide Q by a second area.
struct SurfaceChargeDensityStencil
{
    std::array<std::size_t,3> neighbor{};
    std::array<double,9> polynomialWeights{}; // Row-major 3 x 3.
    double origin=0.0,scale=0.0,panelMeasure=0.0; // Meridional ds coordinates.
    bool evenAtAxis=false;
    // Closed-panel query, without Frame/point inversion or extrapolation.
    std::array<double,3> densityWeightsAtArclength(double s) const;
};

// One retained closed interface or explicitly supported open chain: common Hermite integration panels + conservative
// integrated charge state. No shape/material/electrode parameters or clock.
// Piecewise quadratic reconstruction matches three neighboring panel averages;
// it is not a globally C2 charge interpolant. Panel integrals and local jets
// come from this same reconstruction, not two independently prescribed fields.
// In RZ, charge and panel averages use physical surface area 2*pi*r*ds;
// polynomial coordinates and both charge derivatives still use meridional ds.
// Open-chain pole panels use an even quartic in ds from their AXIS endpoint,
// matching three one-sided panel averages and enforcing zero meridional derivative there. Declared
// non-axis endpoints instead use an unconstrained one-sided quadratic.
class SurfaceChargeField
{
public:
    SurfaceChargeField(const QuadtreeMesh& mesh,
        const LevelSetHermiteGeometry::Snapshot& snapshot, double initialDensity,
        bool axisymmetric = false, std::vector<Vector2D> boundaryContacts = {});
    // Explicit cell-local query mode on the retained immutable scalar. The
    // layout may be RootOnly or explicitly UniformLeaves. source() requires
    // a Frame carrying the matching cell-local chart; no global fallback is
    // used. Old non-sector Panel frames may lack that provenance: callers
    // must query the retained scalar/graph explicitly, not just reuse point.
    // Their saved physical GL8 points are adapted only inside the known cell,
    // with exact physical round-trip checks; the local jet can differ in last
    // bits from the old domain-normalized jet. Sector provenance is never
    // recovered from physical display coordinates.
    // Non-sector endpoint frames retain the original graph edge-root parameter
    // instead: its local-map display point can differ from the old GridEdge
    // display point by roundoff. No endpoint parameter, Panel or saved GL8
    // integral is changed to force those two display conventions to agree.
    SurfaceChargeField(const QuadtreeMesh& mesh,
        std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> snapshot, double initialDensity,
        SurfaceChargeGeometryOptions options = {});
    // Explicit fixed mixed-leaf layout. Reuses the same charge moments,
    // integrated state and material-overlap algorithm on native curve panels.
    SurfaceChargeField(const QuadtreeMesh& mesh,
        LevelSetInterfacePanels::NativeGraphChain chain, double initialDensity);
    // Already completed native material geometry: own a value copy/move of
    // the same immutable panels instead of reintegrating their scalar. Mesh
    // identity is still validated; non-native panels are not accepted here.
    SurfaceChargeField(const QuadtreeMesh& mesh,
        LevelSetInterfacePanels panels, double initialDensity);
    SurfaceChargeField(const SurfaceChargeField&) = delete;
    SurfaceChargeField& operator=(const SurfaceChargeField&) = delete;
    SurfaceChargeField(SurfaceChargeField&&) = delete;
    SurfaceChargeField& operator=(SurfaceChargeField&&) = delete;

    // Returned callback borrows this object's immutable panel geometry, but
    // retains its own reconstructed charge coefficients for this time level.
    // It must not outlive this object. The legacy constructor retains no
    // Snapshot; the explicit shared constructor retains its immutable source.
    // Native queries reuse only a previously validated frame's geometric
    // location on this Field's immutable source. Trial polynomials are never
    // cached. The bounded per-Field location cache is synchronized for
    // concurrent read-only callbacks; publishing charge still needs exclusive
    // caller ownership, as before. Every query checks the Core FP environment.
    SurfaceChargeSource source() const;
    // Reconstruct a trial time/iteration level on this SAME immutable geometry,
    // without changing the stored charge. The callback owns its trial polynomial
    // coefficients and borrows this Field exactly as source() does.
    SurfaceChargeSource source(const std::vector<double>& integratedCharge) const;
    // Same three-neighbor reconstruction as source(), expressed as a linear
    // stencil on this immutable panel geometry. Does not read or change Q,
    // create a density state, or attach a physical clock.
    SurfaceChargeDensityStencil densityStencil(std::size_t panel) const;
    // Publish caller-accepted integrals only after validating their size,
    // finiteness and representable panel densities. Failure leaves charge intact.
    void assignIntegratedCharge(const std::vector<double>& integratedCharge);
    // Read-only native-cell -> material-panel current binding. Uses the SAME
    // bulk-current rejection as explicit advance; bulk charge is not supported
    // and is never redistributed. Only the two bulk-current report fields are
    // replaced on success. A failed mapping leaves report and Field unchanged.
    std::vector<double> mapNormalCurrent(const std::vector<double>& cellCurrent,
        double currentScale, double relativeBulkTolerance, SurfaceChargeFieldStep& report) const;
    // Fixed retained panels: face i points from panel i into (i+1)%N.
    // The caller supplies velocity relative to these panel boundaries.
    // An open chain instead has N+1 boundaries, including two exact zero fluxes.
    SurfaceChargeFieldStep advance(const std::vector<double>& cellCurrent,
        const std::vector<double>& relativeAdvectiveMeasureFlux,
        double currentScale, double relativeBulkTolerance, double dt);
    // Fixed open chain with explicit external inflow densities. Non-axis end
    // fluxes may be nonzero; any Axis end still requires exactly zero flux.
    // Boundary current is reported separately from integrated normal current.
    SurfaceChargeFieldStep advance(const std::vector<double>& cellCurrent,
        const std::vector<double>& relativeAdvectiveMeasureFlux,
        double currentScale, double relativeBulkTolerance, double dt,
        const SurfaceChargeBoundaryState& boundary);
    SurfaceChargeFieldStep advanceNormalCurrent(const std::vector<double>& cellCurrent,
        double currentScale, double relativeBulkTolerance, double dt);
    // Sample the saved shared endpoints, not newly projected geometry. For
    // these fixed panels the advective measure flux is u.t (planar) or
    // 2*pi*r*(u.t) (RZ), with t=(-ny,nx). Declared boundary chains reject this
    // sampler: wall-aware velocities and caller-supplied explicit fluxes are
    // required; no extrapolated or zero wall velocity is silently supplied.
    std::vector<double> advectiveMeasureFlux(const StaggeredVelocitySampler& velocity) const;

    // Material endpoints follow the SAME characteristic map as the level set.
    // Build a new-time field from the actual new scalar, add old-time current,
    // then conservatively overlap the advected old panels with new cut panels.
    // Endpoint motion already includes tangential advection: no extra u.t flux.
    // Overlap coordinates use cumulative ds in planar mode and 2*pi*r*ds in RZ.
    // Local projection is limited to h/4; topology must remain unchanged.
    // Open RZ chains trace all N+1 endpoints. Each pole trajectory must remain
    // exactly on r=0 and is matched to the same oriented new-time axis endpoint;
    // their physical-area overlap is nonperiodic. No radius clipping is used.
    // Declared boundary contacts are not material endpoints and are rejected;
    // moving boundary-contact remapping is not implemented by this method.
    // The explicit cell-local mode is also rejected: this legacy remap must
    // not silently rebuild it in global-coordinate/RootOnly mode.
    // This object and report are unchanged if any part fails.
    std::unique_ptr<SurfaceChargeField> advanceMoving(const QuadtreeMesh& mesh,
        const LevelSetHermiteGeometry::Snapshot& nextSnapshot,
        const std::function<Vector2D(Vector2D)>& forwardMaterialPoint,
        const std::vector<double>& cellCurrent, double currentScale,
        double relativeBulkTolerance, double dt, SurfaceChargeFieldStep& report) const;

    // Explicit retained-chart OPEN-chain material remap on the same native
    // grid and stored scope. fixedMaterialContacts must exactly identify all
    // declared endpoints; forwardMaterialPoint must leave each of them fixed.
    // Other endpoints follow that map once. The old local root plus its
    // physical displacement is the projection query, never a global inverse
    // of the moved display point. Axis endpoints keep their oriented identity
    // and remain at r=0. Remapping uses per-panel local physical measures, so
    // a short panel is not lost in a whole-chain cumulative coordinate.
    // Old motion/projection/residual/conservation gates remain in force.
    // This does not infer a contact law or add a second tangential-advection
    // update. The original field and caller's report are unchanged on failure.
    // Optional preparedPanels is a synchronous borrow of completed NATIVE
    // geometry on exactly nextSnapshot, with the same symmetry and contacts.
    // The new Field owns a copy; no borrow escapes. Invalid supplied geometry
    // rejects, never silently rebuilds. Without it the existing Core factory
    // builds those same panels. Scalar/grid and remap gates are unchanged.
    std::unique_ptr<SurfaceChargeField> advanceMoving(const QuadtreeMesh& mesh,
        std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> nextSnapshot,
        const std::vector<Vector2D>& fixedMaterialContacts,
        const std::function<Vector2D(Vector2D)>& forwardMaterialPoint,
        const std::vector<double>& cellCurrent, double currentScale,
        double relativeBulkTolerance, double dt, SurfaceChargeFieldStep& report,
        const LevelSetInterfacePanels* preparedPanels = nullptr) const;

    // Pure material remap on the retained OPEN-chain path described above.
    // Used before an implicit current update. Additional histories follow the
    // exact SAME local overlap coordinates, not separately retraced endpoints.
    // No physical time step is invented to invoke a zero-current advance.
    // Input Field/history are unchanged on failure; no partial result escapes.
    // Trajectory callbacks remain sequential in old-endpoint order. Native
    // new-interface projections may run independently in parallel; overlap
    // deposition and diagnostic accumulation retain their original order.
    SurfaceChargeMaterialRemap remapMoving(const QuadtreeMesh& mesh,
        std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> nextSnapshot,
        const std::vector<Vector2D>& fixedMaterialContacts,
        const std::function<Vector2D(Vector2D)>& forwardMaterialPoint,
        const std::vector<std::vector<double>>& history = {},
        const LevelSetInterfacePanels* preparedPanels = nullptr) const;

    // Prepare the existing retained-chart OPEN native material map once.
    // targetPanels must be completed on the actual target immutable scalar.
    // Declared contacts are unsupported here; no law or fake endpoint is added.
    // This wrapper preserves the original donor-integral physical-overlap
    // remap, including local chart/projection/motion/coverage rejection.
    SurfaceChargeMaterialTransfer prepareMaterialTransfer(const QuadtreeMesh&,
        const LevelSetInterfacePanels& targetPanels,
        const std::function<Vector2D(Vector2D)>& forwardMaterialPoint,
        int step,double sourceTime,double targetTime) const;
    // Consume the same opaque geometry map for Q and all supplied histories.
    // No endpoint retrace or target geometry reintegration; source partition
    // and actual native source are checked. Failure leaves this Field intact.
    SurfaceChargeMaterialRemap remapPrepared(const QuadtreeMesh&,
        const SurfaceChargeMaterialTransfer&,
        const std::vector<std::vector<double>>& history={}) const;

    // Consume already completed native geometry and overlap weights. Copies
    // the target panels; no projection callback or reintegration is performed.
    // Applies the identical map to current Q and every supplied history.
    SurfaceChargeMaterialRemap remapPrepared(const QuadtreeMesh&,
        const SurfaceChargePreparedTransfer&,const std::vector<std::vector<double>>& history={}) const;
    SurfaceChargeMaterialRemap remapSectionDetachment(const QuadtreeMesh&,
        const SurfaceChargeSectionDetachment&,const std::vector<std::vector<double>>& history={}) const;
    std::unique_ptr<SurfaceChargeField> advancePrepared(const QuadtreeMesh&,
        const SurfaceChargePreparedTransfer&,const std::vector<double>& cellCurrent,
        double currentScale,double relativeBulkTolerance,double dt,SurfaceChargeFieldStep&) const;

    // AMR between accepted physical steps, restricted to a native Q4 curve
    // whose cut leaves retain their exact polynomial and physical bounds.
    // sourceToTargetCell uses SOURCE native indices (e.g. Core Q4 transfer's
    // verified map); it is checked against BOTH immutable scalar authorities,
    // never trusted as proof by itself. The old panel partition/orientation,
    // endpoint charts, ds/RZ areas and quadrature moments must be unchanged.
    // Rebuilds mesh bindings and charge reconstruction on the new snapshot;
    // old geometry/callbacks remain valid while their Field remains alive.
    // Cannot split/merge cut panels or replace a physical material remap.
    // All input histories use the old panel order; no input mutates on failure.
    // Optional prebuilt Core panels must own the exact target scalar. They
    // avoid repeated geometry discovery; all original transfer checks remain.
    SurfaceChargeMeshTransfer rebindMesh(const QuadtreeMesh& targetMesh,
        std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> targetSnapshot,
        const std::vector<int>& sourceToTargetCell,
        const std::vector<std::vector<double>>& history = {},
        const LevelSetInterfacePanels* transferredPanels = nullptr) const;

    const LevelSetInterfacePanels& geometry() const noexcept { return panels_; }
    const SurfaceChargeTransport& transport() const noexcept { return charge_; }
    struct MeasureInterpolationReport {
        std::size_t rationalArcIntervals=0,rationalAreaIntervals=0;
    };
    // Cubic intervals whose positive endpoint metrics require a monotone
    // rational interpolant. This reports coverage, not an integration error.
    const MeasureInterpolationReport& measureInterpolation() const noexcept {return measureInterpolation_;}

private:
    // Fixed-geometry charge transaction stages every component before publishing Q.
    friend class ElectricSurfaceChargeTimeState;
    // Only one live fixed-geometry time coordinator may publish this Field.
    std::weak_ptr<const void> chargeTimeOwner_;
    struct MaterialRemapPlan
    {
        std::unique_ptr<SurfaceChargeField> field;
        std::vector<SurfaceChargeTransport::PanelCoordinate> boundaries;
        SurfaceChargeFieldStep report;
    };
    MaterialRemapPlan prepareMaterialRemap(const QuadtreeMesh& mesh,
        std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> nextSnapshot,
        const std::vector<Vector2D>& fixedMaterialContacts,
        const std::function<Vector2D(Vector2D)>& forwardMaterialPoint,
        const LevelSetInterfacePanels* preparedPanels) const;
    SurfaceChargeSource source(const SurfaceChargeTransport& charge) const;
    struct MeasureNode { double coordinate, measure, derivative; bool rationalToNext=false; };
    struct Mapping
    {
        Vector2D low, high;
        double coordinateBegin = 0.0, coordinateSpan = 0.0, scale = 0.0;
        double polynomialOrigin = 0.0;
        bool evenAtAxis = false;
        std::array<std::size_t, 3> neighbor{};
        std::vector<MeasureNode> arc;
        std::vector<MeasureNode> physical;
        std::array<double, 9> inverseMoments{};
    };
    struct Polynomial { double a, b, c, scale, origin; bool evenAtAxis; };
    Polynomial reconstructPolynomial(std::size_t panel,double panelDensity,
        const std::array<double,3>& neighborDensity) const;
    struct Location { std::size_t panel; double arclength, localPhysicalMeasure; };
    using NativeLocationKey = std::array<std::uint64_t, 11>;
    mutable std::mutex nativeLocationMutex_;
    mutable std::map<NativeLocationKey, Location> nativeLocations_;
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source_;
    LevelSetInterfacePanels panels_;
    SurfaceChargeTransport charge_;
    std::vector<Mapping> mappings_;
    MeasureInterpolationReport measureInterpolation_;
    std::map<std::pair<int, int>, std::size_t> panelByCell_;
    CartesianGridScope scope_ = CartesianGridScope::RootOnly;
    std::vector<int> retainedCellLabels_, retainedVertexLabels_;
    double coordinateTolerance_ = 0.0;

    void initialize(const QuadtreeMesh& mesh);
    void validateNativeRemapGrid(const QuadtreeMesh& mesh,
        const LevelSetHermiteGeometry::Snapshot& next) const;
    double nativeFrameFraction(std::size_t panel,const LevelSetHermiteGeometry::Frame& frame) const;
    LevelSetHermiteGeometry::Frame retainedFrame(std::size_t panel,
        const LevelSetHermiteGeometry::Frame& frame, int endpoint) const;
    void validateLocalFrame(std::size_t panel, const LevelSetHermiteGeometry::Frame& frame) const;
    static std::size_t prepareMapping(std::vector<MeasureNode>& nodes,
        const char* orderReason, const char* monotonicityReason,
        bool zeroBeginDerivative = false, bool zeroEndDerivative = false);
    static double interpolateMapping(const std::vector<MeasureNode>& nodes, double coordinate);
    Location locate(Vector2D point) const;
    Location locate(const LevelSetHermiteGeometry::Frame& frame) const;
    Location locateNativeCached(const LevelSetHermiteGeometry::Frame& frame) const;
    Location locateFraction(std::size_t panel,double coordinate) const;
    std::vector<double> panelCurrent(const std::vector<double>& cellCurrent,
        double currentScale, double relativeBulkTolerance, SurfaceChargeFieldStep& report) const;
    SurfaceChargeJet evaluate(const Location& location, const std::vector<Polynomial>& polynomial) const;
};
}
