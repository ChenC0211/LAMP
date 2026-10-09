#pragma once
#include "Modules/Electricity/FittedSharpDielectricElectrostatics.h"
#include "Modules/Electricity/SurfaceChargeField.h"
#include "Modules/Electricity/SurfaceChargeRelaxation.h"
#include <functional>
#include <cstdint>
#include <limits>

namespace Quantum::Electricity {
// Normalized numerical inputs only. SI conversion and boundary IDs belong to
// the caller. The existing active RZ P3 capability and CFL <= 0.25 apply to
// both time policies; this is not a new large-step implicit method.
struct FittedChargeParameters
{
    DielectricConfig field;
    DielectricConductivity conductivity;
    double closureTolerance = 0.0;
};
struct FittedChargeProblem
{
    const FaceSegmentGraph& graph;
    const std::vector<PolarizationPoissonFVSolverCellCentered::FittedBoundary>& boundaries;
    FittedChargeParameters parameters;
};

// One complete electrical time level, published only by this module. Copies
// retain the same immutable charge, field and two-level history without a
// deep copy. A successful candidate still needs the caller's coupled commit.
class FittedChargeState
{
public:
    FittedChargeState() = default;
    explicit operator bool() const noexcept { return bool(data_); }
    std::shared_ptr<const SurfaceChargeField> charge() const;
    // Stable Interface surface identity, independent of vector order. The
    // no-index accessors reject a collection instead of selecting a child.
    std::size_t componentCount() const;
    std::uint64_t surfaceId(std::size_t component) const;
    std::shared_ptr<const SurfaceChargeField> charge(std::size_t component) const;
    const std::vector<double>& previousMaterialCharge(std::size_t component) const;
    const FittedSharpDielectricResult& electric() const;
    const std::vector<double>& previousMaterialCharge() const;
    double previousChargeDt() const;
    int completedStep() const;
    double time() const;
    void validateMesh(const QuadtreeMesh&) const;
    // Unlike SurfaceChargeField::source(), this callback owns its immutable
    // time level and remains valid when the App replaces/destroys this handle.
    SurfaceChargeSource source() const;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
    friend struct FittedChargeAccess;
    friend class FittedChargeMeshTransfer;
};

// Same-time charge/history rebind, before the corresponding new-grid solve.
// This is deliberately not a solved state. The caller may inspect panels to
// bind Interface/Hydro, but cannot edit charge/history or mark it accepted.
class FittedChargeMeshTransfer
{
public:
    std::shared_ptr<const SurfaceChargeField> charge() const;
private:
    std::shared_ptr<const FittedChargeState::Data> data_;
    friend struct FittedChargeAccess;
};
struct FittedChargeStepReport
{
    double dt = 0.0; // Code time, matching the interface motion candidate.
    double relaxationRate = 0.0; // max(code conductivity/permittivity).
    double relaxationCfl = 0.0;
    // Independent whole-domain charge/electrode accounting in electric code
    // units: Q_after-Q_before+dt*sum(J_old outward). This is distinct from
    // SurfaceChargeFieldStep's panel transport/remap conservation checks.
    double chargeClosureResidual = 0.0;
    double chargeClosureScale = 0.0;
    double relativeChargeClosure = 0.0;
    // Present only for the implicit branch. Closure then uses the actual BDF
    // derivative and NEW electrode current; explicit accounting is unchanged.
    std::optional<SurfaceChargeRelaxationReport> relaxation;
    int electricFactorReuses = 0, electricGeometryReuses = 0;
    int electricPreparedRhsSolves = 0; // Console performance observation, not an acceptance gate.
    FittedSharpDielectricResult::Timings electricTimings;
    double electricMatrix_s = 0.0, electricFactor_s = 0.0, electricBacksolve_s = 0.0;
    double materialTransfer_s = 0.0;
};

// Final, volume-corrected geometry and its material map, supplied by the
// caller. panels is borrowed for this synchronous call only. Neither scalar
// transport nor force assembly is repeated here.
struct FittedChargeMotion
{
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> geometry;
    const LevelSetInterfacePanels* panels = nullptr;
    std::vector<Vector2D> fixedContacts;
    std::function<Vector2D(Vector2D)> forwardMaterialPoint;
    double dt = 0.0;
    double sourceTime = std::numeric_limits<double>::quiet_NaN();
    double targetTime = std::numeric_limits<double>::quiet_NaN();
    // Exclusive alternative to forwardMaterialPoint, borrowed only during
    // prepare. Geometry/time bindings are checked against this motion and state.
    const SurfaceChargePreparedTransfer* preparedTransfer = nullptr;
    // Optional local closure after the physical leg, at targetTime. That leg
    // ends on sectionDetachment->source; the event ends on panels/geometry.
    // Old normal current is still applied exactly once.
    const SurfaceChargeSectionDetachment* sectionDetachment = nullptr;
    // Optional Interface identity for the final single-surface geometry.
    // A changed identity requires the explicit section event above.
    std::uint64_t targetSurfaceId = 0;
};
struct FittedChargeCandidate
{
    FittedChargeState state;
    SurfaceChargeFieldStep chargeStep;
    FittedChargeStepReport report;
};
struct FittedChargeComponentInitial
{
    std::uint64_t surfaceId = 0;
    const LevelSetInterfacePanels* panels = nullptr;
    double initialDensity = 0.0;
};
struct FittedChargeComponentMotion
{
    std::uint64_t surfaceId = 0;
    FittedChargeMotion motion;
};
struct FittedChargeCollectionCandidate
{
    FittedChargeState state;
    std::vector<SurfaceChargeFieldStep> chargeSteps; // Output surface order.
    FittedChargeStepReport report; // One global electrode balance.
    double unassignedBulkCurrent = 0.0, absoluteBulkCurrent = 0.0;
    bool interiorSplit = false;
    double eventRemapResidual = 0.0;
};
struct FittedChargeInteriorSplit
{
    struct Target {std::uint64_t surfaceId=0;const LevelSetInterfacePanels* panels=nullptr;};
    struct Weight {std::size_t targetSurface=0,targetPanel=0;long double share=0;};
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> geometry;
    std::vector<Target> targets;
    std::vector<std::vector<Weight>> weights;
    std::array<int,4> cells{};
    double section=0;
};
// Observes completed field evaluations; never supplies a numerical current.
// Arguments: evaluation count, iteration limit, geometry reuse, factor reuse.
using FittedChargeObserver = std::function<void(int, int, bool, bool)>;

FittedChargeStepReport preflightFittedChargeStep(const FittedChargeParameters&, double dt);
FittedSharpDielectricResult solveFittedChargeField(
    const QuadtreeMesh&, const FittedChargeProblem&,
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot>,
    const LevelSetInterfacePanels&, const SurfaceChargeSource&,
    const FittedPotentialInput& previousPotential, FittedSharpDielectricSolver* repeatedSolver = nullptr);

FittedChargeState initializeFittedChargeState(const QuadtreeMesh&, const FittedChargeProblem&,
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot>, const LevelSetInterfacePanels&,
    const std::vector<Vector2D>& fixedContacts, const FittedPotentialInput& previousPotential,
    double initialDensity = 0.0);
// Native collection: one scalar, complete surfaces, unique nonzero identities
// and distinct cut-cell ownership. Two surfaces in one cut cell need a
// resolved phase-current partition and are explicitly unsupported here.
FittedChargeState initializeFittedChargeState(const QuadtreeMesh&, const FittedChargeProblem&,
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot>,
    const std::vector<FittedChargeComponentInitial>&, const FittedPotentialInput&);
// Explicit same-topology material continuation. Applies the old global
// current once, remaps each child, then solves one new global electric field.
// All candidates are discarded together on failure. Split-event creation and
// implicit collection histories are separate operations, never inferred here.
FittedChargeCollectionCandidate prepareFittedChargeCollectionStep(
    const QuadtreeMesh&, const FittedChargeProblem&, const FittedChargeState&,
    const std::vector<FittedChargeComponentMotion>&);
// Physical current/material leg to a connected target, then a zero-time
// local 1:2 split. Remote panels must retain their exact native identities;
// local rows cannot transfer charge across the declared axial section.
FittedChargeCollectionCandidate prepareFittedChargeInteriorSplitStep(
    const QuadtreeMesh&,const FittedChargeProblem&,const FittedChargeState&,
    const FittedChargeMotion& physical,const FittedChargeInteriorSplit&);
FittedChargeMeshTransfer prepareFittedChargeMeshTransfer(const QuadtreeMesh& sourceMesh,
    const FittedChargeState&, const QuadtreeMesh& targetMesh,
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot>, const std::vector<int>& sourceToTargetCell,
    const LevelSetInterfacePanels* transferredPanels = nullptr);
FittedChargeState solveFittedChargeMeshTransfer(const QuadtreeMesh&, const FittedChargeProblem&,
    const FittedChargeMeshTransfer&, std::shared_ptr<const LevelSetHermiteGeometry::Snapshot>,
    const FittedPotentialInput& previousPotential);

// Explicit: old current once, conservative material remap, new P3 solve.
// Implicit: remap both histories, original BE/two-level relaxation, electrode
// closure from the SAME accepted charge/current iterate. Old state and motion
// remain unchanged even when a solve, remap, closure or observer throws;
// only a fully completed candidate is returned. Commit belongs to the caller.
FittedChargeCandidate prepareFittedChargeStep(
    const QuadtreeMesh&, const FittedChargeProblem&, const FittedChargeState&,
    const FittedChargeMotion&, const SurfaceChargeRelaxationOptions* relaxation = nullptr,
    const FittedChargeObserver& observer = {});
}
