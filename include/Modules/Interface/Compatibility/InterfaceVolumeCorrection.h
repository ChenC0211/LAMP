#pragma once

#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <vector>
#include "Core/Geometry/ActiveHermitePatch.h"
#include "Core/Geometry/ActiveQuarticHermitePatch.h"
#include "Core/Geometry/LevelSetInterfacePanels.h"
#include "Modules/Interface/InterfaceCellMeasure.h"
#include "Modules/Interface/Compatibility/InterfacePhaseFlux.h"
#include "Modules/Interface/Compatibility/InterfaceState.h"
#include "Modules/Interface/Compatibility/InterfaceTemporalPhaseFlux.h"



namespace InterfaceVolumeCorrection
{
struct GradientDof
{
    int vertexIndex = -1; // Actual native vertex index, independent of fixedValue.
    unsigned component = 0; // 1: dx, 2: dy. Mixed derivatives are not unknowns.
    double scale = 1.0; // deltaGradient = phi-equivalent unknown * scale.
};
struct LinearReport
{
    double linearResidual = 0.0;
    double normalEquationResidual = 0.0; // Max relative KKT stationarity/constraint block residual.
    int rows = 0, columns = 0;
    bool usedSuperLU = false;
    // Diagnostics for the largest |defect / free-row norm|. These do not
    // suppress small rows or alter the full set of linear constraints.
    int largestScaledDefectCell = -1;
    long double largestScaledDefect = 0.0L, correspondingRowNorm = 0.0L;
};
struct Result : LinearReport
{
    // Full native vertex indexing supplied by fixedValue. Fixed and unused
    // values are exactly zero. Explicit appended gradients have separate IDs.
    std::vector<double> deltaValue;
    std::vector<double> deltaGradient; // Physical dx/dy increment, appended-list order.
};

struct QuarticDof
{
    using Kind=Quantum::Geometry::ActiveQuarticHermiteSensitivity::Kind;
    Kind kind=Kind::Vertex;
    int index=-1,component=0; // Same independent physical identity as the scalar.
    double scale=1.0; // deltaData = phi-equivalent unknown * scale.
};
struct QuarticResult : LinearReport
{
    std::vector<double> deltaData; // Physical increments, explicit freeData order.
};
// Rows of R in ||delta||^2 + ||R delta||^2. Coefficients act on physical
// independent data; the kernel applies each declared phi-equivalent scale.
using QuarticRegularization = std::vector<std::vector<Quantum::Geometry::ActiveQuarticHermiteSensitivity>>;

// Original VALUES-ONLY Newton LINEAR step, without changing the patch or
// rebuilding geometry. This interface does not append gradient unknowns.
// Only native cells containing free arcs contribute volume constraints.
// Arc sensitivity[4*corner] identifies a shared VALUE DOF through the arc's
// source vertexIndices, not through its actual control-volume owner.
// Explicit native-source arcs instead consume independentSensitivity, already
// pulled back through the scalar's hanging constraints. Their formal corner
// arrays are never solve columns. The caller still fixes pin/axis values and
// selects admissible REGULAR gradient unknowns; sector variations are absent.
// fixedValue is the caller's full native-vertex binary mask; targets use the
// same full native-cell indexing as Measure::negativeVolumeByCell. Targets
// outside the constrained rows are not consumed or claimed to be satisfied.
//
// Row-2-norm scaling leaves the equality-constrained Euclidean minimum-norm
// problem unchanged: solve [I A^T; A 0] [delta;lambda]=[0;b] with SuperLU.
// Zero/free-empty rows, failed factorization, nonfinite/unrepresentable data,
// and complete normalized residuals above linearTolerance throw. There is
// no pin, damping, pseudoinverse, rank truncation, or alternative backend.
// No Gram matrix is formed; a successful linear step is NOT a
// nonlinear conservation certificate, geometric acceptance, or volume fix.
Result solve(const InterfaceCellMeasure::Measure& measure,
    const std::vector<long double>& targetNegativeVolumeByCell,
    const std::vector<std::uint8_t>& fixedValue, double linearTolerance);

// Same linear system as solve(), but consumes the desired volume increment
// directly. This avoids forming a large target volume only to subtract the
// measured volume again. Rows without free arcs remain outside this LINEAR
// system; correct() below still includes every active cell in its residual.
// Optional gradient columns use sensitivity[4*corner+component]*scale. Their
// unknowns have phi-equivalent units in the same minimum-norm problem; their
// returned increments are multiplied by scale. fixedValue does not fix them.
Result solveDefect(const InterfaceCellMeasure::Measure& measure,
    const std::vector<long double>& deltaVolumeByCell,
    const std::vector<std::uint8_t>& fixedValue, double linearTolerance,
    const std::vector<GradientDof>& freeGradients = {});

// Same minimum-norm/SuperLU kernel, with residual derivative dV+halfDt*dQ.
// The caller must supply the measure and flux variation of the SAME current
// scalar/cuts/frozen velocity. Flux-only rows are included, not disguised as
// geometric arcs. At exact face endpoints dQ is a symmetric generalized slope.
Result solveBalanceDefect(const InterfaceCellMeasure::Measure& measure,
    const InterfacePhaseFlux::FluxSensitivity& fluxSensitivity, long double halfDt,
    const std::vector<long double>& negativeResidualByCell,
    const std::vector<std::uint8_t>& fixedValue, double linearTolerance,
    const std::vector<GradientDof>& freeGradients = {});

// Same row normalization, minimum-norm/SuperLU solve and complete residual
// gates, with typed native quartic dV+halfDt*dQ. Only explicitly listed shared
// data can change; omitted data (pins, walls, fixed sectors) remain fixed.
// This is one LINEAR step, not geometry reconstruction or nonlinear acceptance.
QuarticResult solveQuarticBalanceDefect(const InterfaceCellMeasure::QuarticMeasure& measure,
    const InterfacePhaseFlux::QuarticFluxSensitivity& fluxSensitivity, long double halfDt,
    const std::vector<long double>& negativeResidualByCell,
    const std::vector<QuarticDof>& freeData, double linearTolerance,
    const QuarticRegularization& regularization = {});

struct Progress;
struct IterationOptions
{
    double linearTolerance = 1e-10;
    // Explicit nonlinear-iteration gates, not historical physical validation
    // tolerances. The denominator uses per-CV net outward phase fluxes.
    long double relativeBalanceTolerance = 1e-3L;
    // Optional dimensional L1 volume budget for zero/near-zero transport.
    // |R|_1 <= absoluteVolumeTolerance + relativeBalanceTolerance*|dV_flux|_1.
    // Zero preserves the legacy gate. Never added to a diagnostic denominator.
    long double absoluteVolumeTolerance = 0.0L;
    long double maxCellRelativeTolerance = 1e-10L;
    int maxIterations = 12, maxBacktracks = 8;
    std::size_t maxRootRefinements = 256;
    // Synchronous observer; the kernel itself does not print or write files.
    // Observer exceptions abort the correction, not a geometric backtrack.
    std::function<void(const Progress&)> progress;
};
struct Options : IterationOptions
{
    std::vector<GradientDof> freeGradients; // Caller-declared first derivatives only.
};
struct QuarticOptions : IterationOptions
{
    // Explicit independent values and first derivatives; mixed derivatives
    // stay fixed, as in the original iteration. Wall/pin values cannot be free.
    std::vector<QuarticDof> freeData;
    // Complete closed-component path only. Zero retains Euclidean correction.
    double geometricRegularization = 0;
};
struct Residual
{
    long double sumAbsolute = 0.0L, transportedMagnitude = 0.0L;
    long double relativeBalance = 0.0L, maxRelativeCell = 0.0L, signedSum = 0.0L;
    // Worst |R_K|/capacity_K, retaining the signed local balance ingredients.
    int worstCell = -1;
    long double worstResidual = 0.0L, worstOldVolume = 0.0L, worstVolume = 0.0L;
    long double worstTransport = 0.0L, worstCapacity = 0.0L;
    // Complete trajectory path: raw target violation and numerical budget.
    // Neither the target nor a shared-face flux is clipped to enforce bounds.
    long double maxTargetViolation = 0.0L, maxTargetUncertainty = 0.0L;
};
struct Report
{
    // Observational wall times, accumulated across completed candidate stages.
    // Not solver tolerances or physical diagnostics; do not use for acceptance.
    struct Timing
    {
        double geometry = 0.0, facePartitions = 0.0, volume = 0.0;
        double phaseFlux = 0.0, balanceChecks = 0.0, linearSolve = 0.0;
        int evaluations = 0;
    } timing;
    int iterations = 0, backtracks = 0, rows = 0, columns = 0;
    bool usedSuperLU = false;
    double maxCorrection = 0.0; // Actual max value increment (quartic: also edge/center), code phi.
    double maxLinearCorrection = 0.0;
    double maxGradientCorrection = 0.0, maxLinearGradientCorrection = 0.0; // Physical first-derivative units.
    int largestScaledDefectCell = -1;
    long double largestScaledDefect = 0.0L, correspondingRowNorm = 0.0L;
    Residual initial, final;
};
struct Progress
{
    enum class Stage { PredictorReady, LinearStepReady, TrialStarted, TrialRejected, TrialAccepted, Converged };
    Stage stage;
    Report report; // final is the latest accepted iterate, never a rejected trial.
    int backtrack = 0;
    double damping = 1.0;
    std::string detail;
};
template<class PatchType,class MeasureType> struct InterfaceCorrection
{
    PatchType patch;
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> geometry;
    std::unique_ptr<LevelSetInterfacePanels> panels;
    MeasureType measure;
    InterfacePhaseFlux::FluxBudget flux;
    Report report;
    // Exactly one geometry representation is returned: correct() owns legacy
    // panels; correctNative() owns the complete same-source native chain.
    // A native result does not supply the unconnected legacy motion/load API.
    std::unique_ptr<LevelSetInterfacePanels::NativeGraphChain> nativeChain;
    // Complete XY/RZ correction owns all components; the legacy single-chain
    // entries leave this empty. It has geometry only, no committed time step.
    std::unique_ptr<InterfaceState> completeGeometry;
    std::unique_ptr<InterfaceTemporalPhaseFlux::Integral> stepFlux;
};
using Correction=InterfaceCorrection<Quantum::Geometry::ActiveHermitePatch,InterfaceCellMeasure::Measure>;
using QuarticCorrection=InterfaceCorrection<Quantum::Geometry::ActiveQuarticHermitePatch,InterfaceCellMeasure::QuarticMeasure>;

// Explicit old-time authority for pinned chains. Borrowed only during correct;
// the accepted result owns its integrated ledger. No second physical clock.
template<class Patch> struct TemporalInput
{
    const Patch& previousPatch;
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> previousGeometry;
    const InterfaceNormalFluxView& previousVelocity;
    const InterfaceTemporalPhaseFlux::PhaseTrajectory<Patch>* trajectory=nullptr;
    InterfaceTemporalPhaseFlux::Options options;
};
using CubicTemporalInput=TemporalInput<Quantum::Geometry::ActiveHermitePatch>;
using QuarticTemporalInput=TemporalInput<Quantum::Geometry::ActiveQuarticHermitePatch>;

// Bounded, damped volume-defect iteration for the already computed old/new
// velocity pair. u^(n+1) is frozen; its phase flux is NOT: every trial rebuilds
// one Snapshot -> Panels -> refined mapped partitions -> volume AND flux.
// R_K = (V_K - V_K^n) + dt/2*(Q_K^n + Q_K(phi,u^(n+1))). solveDefect uses dV
// only, so this is a quasi-Newton/Picard correction, not the full flux Jacobian.
// All active cells, including pure-phase/no-arc cells, enter both final gates.
// Only unfixed native nodal VALUES and explicitly declared first derivatives
// change. All other derivatives, solid masks, corner sectors and pinned values
// remain unchanged; no derivative refit is performed. Snapshot must certify
// compatibility with unchanged sectors/axis data after every trial. A caller
// can free the contact's wall-tangential derivative without moving its zero.
// Returned volumes and trapezoidal targets lie in [0, actual CV capacity]
// without clipping. Internal iterates always have bounded ACTUAL volumes.
// Before the balance gates pass, require a strictly decreasing complete L1;
// once they pass, retain those gates and strictly decrease the sum of target
// range violations instead. Finally acceptable candidates take precedence.
// Internal flux-derived targets may temporarily be outside the bounds, so
// newly cut cells can enter the next linearization. This is not acceptance of
// a physical step; final target bounds remain mandatory. The merit phases do
// not introduce another tolerance or change the final residual thresholds.
// Failure returns no partially accepted geometry/fields.
// Supplying temporal replaces the trapezoidal target by the common within-step
// integral and uses dV+d(integral Q) in the linearization. Its accepted ledger
// is returned in stepFlux. Every candidate still passes the same nonlinear gates.
// Caller must still apply its transport-band, material, E/force and physical
// validation gates before committing the returned interface to an App state.
Correction correct(const QuadtreeMesh& mesh, const InterfaceTransportDomain& domain,
    const InterfaceCellMeasure::Measure& oldMeasure,
    const InterfacePhaseFlux::FluxBudget& oldFlux,
    const InterfaceNormalFluxView& currentVelocity, double dt,
    const Quantum::Geometry::ActiveHermitePatch& advectedPatch,
    const std::vector<Vector2D>& pinnedCoordinates,
    const std::vector<std::uint8_t>& fixedValue, const Options& options = {},
    const CubicTemporalInput* temporal = nullptr);

// Same nonlinear iteration, balance/target gates, root refinement and SuperLU
// kernel on actual rectangular leaves, using dV+dt/2*dQ (the uniform entry
// retains its original volume-only linearization). Each trial reconstructs the constrained
// native scalar and discovers its COMPLETE current chain; old graph IDs are
// not carried over after a nodal change. Hanging jets are dependent inputs:
// their unused patch payload is neither read nor modified, and appended
// gradient unknowns/pinned zeros must identify regular active vertices.
// An active hanging fixedValue entry must be zero: prescribing its value
// requires constraints on its independent support, not an ignored mask bit.
// Fixed no-slip values, declared pinned contacts and pure pressure-port cells
// keep the same contracts. This is a correction of a supplied predictor,
// not a native transport/Q4 reconstruction or an App AMR switch.
Correction correctNative(const QuadtreeMesh& mesh, const InterfaceTransportDomain& domain,
    const InterfaceCellMeasure::Measure& oldMeasure,
    const InterfacePhaseFlux::FluxBudget& oldFlux,
    const InterfaceNormalFluxView& currentVelocity, double dt,
    const Quantum::Geometry::ActiveHermitePatch& advectedPatch,
    const std::vector<Vector2D>& pinnedCoordinates,
    const std::vector<std::uint8_t>& fixedValue, const Options& options = {});

// Same damped iteration and ALL original balance, raw-volume, target-range,
// pin/sector and pressure-phase gates. Native quartic values/first derivatives
// are updated only through the explicit independent-data list; no Q4 refit.
// Dependent hanging/fine-side data are not unknowns. This corrects a supplied
// predictor, not the transport PDE or an App state. Failure returns no candidate.
QuarticCorrection correctQuarticNative(const QuadtreeMesh& mesh, const InterfaceTransportDomain& domain,
    const InterfaceCellMeasure::QuarticMeasure& oldMeasure,
    const InterfacePhaseFlux::FluxBudget& oldFlux,
    const InterfaceNormalFluxView& currentVelocity, double dt,
    const Quantum::Geometry::ActiveQuarticHermitePatch& advectedPatch,
    const std::vector<Vector2D>& pinnedCoordinates,
    const std::vector<std::uint8_t>& fixedValue, const QuarticOptions& options,
    const QuarticTemporalInput* temporal = nullptr);

// Same native Q4 balance iteration with the complete component inventory.
// Rebuilds all arcs and cell phase certificates at every trial, including XY.
// This low-level result has no material lineage or accepted physical clock.
QuarticCorrection correctQuarticComponents(const QuadtreeMesh& mesh, const InterfaceTransportDomain& domain,
    const InterfaceState& previous,
    const InterfaceNormalFluxView& previousVelocity, const InterfaceNormalFluxView& currentVelocity, double dt,
    const Quantum::Geometry::ActiveQuarticHermitePatch& advectedPatch,
    const std::vector<Vector2D>& boundaryContacts,
    const std::vector<std::uint8_t>& fixedValue, const QuarticOptions& options,
    const InterfaceTemporalPhaseFlux::Options& timeFluxOptions = {},
    const InterfaceTemporalPhaseFlux::Trajectory* trajectory = nullptr);
}
