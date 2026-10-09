#pragma once

#include <cstddef>
#include <functional>
#include <optional>
#include <variant>
#include "Modules/Interface/Compatibility/InterfaceState.h"
#include "Modules/Interface/Compatibility/InterfaceVolumeCorrection.h"
#include "Core/Geometry/ActiveLevelSetTransport.h"
#include "Core/Geometry/LevelSetBoundaryTrace.h"

// Retains the complete typed measure needed by the next correction. No
// quartic-to-cubic sensitivity conversion and no parallel old/new ledgers.
using InterfaceVolumeMeasure = std::variant<InterfaceCellMeasure::Measure,
    InterfaceCellMeasure::QuarticMeasure>;

// Stored inside the owning InterfaceState and exposed only as a const borrow.
// Both arrays describe that state's one geometry; callers cannot publish an
// independently constructed pair through the motion/AMR interfaces.
struct InterfacePhaseState
{
    InterfaceVolumeMeasure volumeMeasure;
    InterfacePhaseFlux::FluxBudget phaseBudget;
    // Accepted step-volume flux, distinct from the instantaneous end flux.
    // Absent for initialization and legacy endpoint-trapezoid workflows.
    struct Transport
    {
        InterfaceStepTime time;
        std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source,target;
        InterfaceTemporalPhaseFlux::Integral integral;
    };
    std::optional<Transport> transport;
};

struct InterfaceMotionReport
{
    // changed/stationary and maxValueIncrement describe the final stored phi;
    // taperedNodes and maxTraceRatio remain predictor/backtrace diagnostics.
    int activeNodes = 0, changedValues = 0, stationaryNodes = 0, taperedNodes = 0;
    int supportCells = 0, pressureCells = 0;
    double dt = 0.0, fullWidth = 0.0, zeroWidth = 0.0;
    double maxValueIncrement = 0.0, maxTraceRatio = 0.0;
    double transportSeconds = 0.0; // Wall time only; never a physical acceptance criterion.
    ActiveLevelSetTransport::QuarticTimings quarticTimings;
    // Outer workflow scopes, excluding transportSeconds above; correction
    // includes its DOF preparation. Finalization excludes return destruction.
    double preparationSeconds=0, correctionPreparationSeconds=0;
    double correctionSeconds=0, finalizationSeconds=0;
    std::size_t transportStencilRequests = 0, transportStencilBuilds = 0;
    std::size_t transportScalarRequests = 0, transportScalarEvaluations = 0;
    double areaBefore = 0.0, areaAfter = 0.0; // RZ free-surface area, code length squared.
    Vector2D pinNormalBefore{0.0, 0.0}, pinNormalAfter{0.0, 0.0};
    Vector2D poleBefore{0.0, 0.0}, poleAfter{0.0, 0.0};
};

// Motion diagnostics/material map. The fresh destination owns the new scalar,
// its volume measure and phase budget, with no second field allocation.
struct InterfaceMotionResult
{
    InterfaceMotionReport report;
    InterfaceVolumeCorrection::Report conservation;
    // Extract the SAME predictor trajectory only for its retained source,
    // final corrected target and actual step times. Diagnostics are not authority for
    // this binding. The returned callback owns the original velocity captures;
    // correction does not advect material again or project points to the target.
    std::function<Vector2D(Vector2D)> materialMapFor(
        const std::shared_ptr<const LevelSetHermiteGeometry::Snapshot>& source,
        const std::shared_ptr<const LevelSetHermiteGeometry::Snapshot>& target, const InterfaceStepTime& time) const;
private:
    friend struct InterfaceMotionAccess;
    friend class InterfaceComponentMotion;
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> materialSource_, materialTarget_;
    std::function<Vector2D(Vector2D)> forwardMaterialPoint_;
    InterfaceStepTime materialTime_;
};

// Optional synchronous, read-only iteration observer. It must not retain the
// borrowed Progress or change candidate state. The driver supplies the writer.
using InterfaceVolumeProgress = std::function<void(const InterfaceVolumeCorrection::Progress&)>;

// Optional synchronous, read-only observation of EXISTING stage geometry.
// No extra reconstruction or candidate acceptance. References expire when the
// callback returns; an observer failure aborts rather than retries the step.
enum class InterfaceMotionStage { Previous, Transported, VolumeCorrected };
using InterfaceGeometryObserver = std::function<void(
    InterfaceMotionStage, const LevelSetInterfacePanels&)>;

struct InterfaceMotionVelocityLayers
{
    // Immutable repeatable fields at the current and previous substep times.
    // Copies are retained by the material trajectory: captures must own their
    // data, or the caller must guarantee that data outlive the material map.
    LevelSetBoundaryTrace::Velocity current, previous;
    // Explicit concurrency promise for both callbacks and all captured data;
    // arbitrary prescribed fields remain sequential unless opted in.
    bool independentQueries = false;
};

namespace InterfaceMotion
{
// Existing RZ pinned-contact-to-axis capability: one wet stationary-wall
// corner, cubic uniform leaves or native quartic leaves. Other topologies
// are not silently substituted. No prescribed contact angle or redistance.
void validatePinnedToAxisChain(const LevelSetInterfacePanels& panels, Vector2D pin);

// Complete existing workflow: support/port separation, characteristic and
// constrained increment reconstruction, local volume correction, final
// geometry/flux and the SAME retained material trajectory. Uses Core geometry
// and supplied frozen velocity/normal flux; no NS or electric algorithm.
// requiredHalfWidthCells is the caller's required physical-load support; the
// original numerical halo/taper is added without changing that load width.
// velocities.previous describes time.begin and current describes time.end.
// currentFlux must describe velocities.current on this SAME domain/time layer.
// It and domain are synchronous borrows, never retained in the result.
// previousFlux, when supplied, must describe velocities.previous. It selects
// the common within-step phase integral and its full correction derivative;
// callers without an old flux trace retain the explicitly legacy trapezoid.
//
// destination must be a fresh, distinct unaccepted workspace on this mesh.
// It is filled only during finalization; on ANY exception discard it. This
// avoids copying a full display field through an App-to-module state wrapper.
// Return is NOT coupled-step acceptance: the App checks the remaining coupled
// physics and controls commit/output publication. The App chooses time;
// success records its index/end without a second running clock.
InterfaceMotionResult preparePinnedToAxis(const QuadtreeMesh& mesh,
    const InterfaceTransportDomain& domain, const InterfaceState& previous, const InterfaceStepTime& time,
    const InterfaceMotionVelocityLayers& velocities, const InterfaceNormalFluxView& currentFlux,
    double requiredHalfWidthCells, double linearTolerance,
    InterfaceState& destination, const InterfaceVolumeProgress& progress = {},
    const InterfaceGeometryObserver& observeGeometry = {},
    const InterfaceNormalFluxView* previousFlux = nullptr);
}
