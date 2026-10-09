#pragma once

#include <optional>
#include "Core/Field/CellField.h"
#include "Core/Field/HEdgeField.h"
#include "Core/Field/VEdgeField.h"
#include "Core/Field/VertField.h"
#include "Modules/Hydrodynamics/HydroSolveReport.h"
#include "Modules/Hydrodynamics/HydroMaterialFields.h"
#include "Modules/Hydrodynamics/HydroForceSource.h"

class HydroFittedDomain;
class Hydrodynamics;
class QuadtreeMesh;
class NativeLeafDomain;
class InterfaceState;
namespace Quantum::Hydrodynamics {class HydroConservativeTransport;}

enum class HydroDomainKind { CompleteXY, CompleteRZ, FittedRZ };
enum class HydroVelocityTimeLevel { Current, Previous };
enum class HydroForceTimeLevel { StepBegin, StepEnd, StepMidpoint };
enum class HydroTimePolicy { Fixed, ConservativeOneStep };

// Borrowed for one synchronous call. Physical boundary data are either the
// complete-domain patches in options or the fitted domain, never both.
struct HydroStepGeometry
{
    HydroDomainKind kind;
    const HydroFittedDomain* fitted = nullptr;
    static HydroStepGeometry complete(bool axisymmetric)
    { return {axisymmetric ? HydroDomainKind::CompleteRZ : HydroDomainKind::CompleteXY, nullptr}; }
    static HydroStepGeometry fittedRZ(const HydroFittedDomain& domain)
    { return {HydroDomainKind::FittedRZ, &domain}; }
};

// Supplied by the App's sole clock, in the SAME code units as the solver.
// The module records identity/continuity; it does not choose dt or run a loop.
struct HydroStepTime
{
    int index = 1;
    double begin = 0.0;
    double dt = 0.0;
    // Conservative donor and coupled midpoint paths are one-step methods.
    // They use the current interval; ordinary CN/AB2 remains fixed.
    HydroTimePolicy policy = HydroTimePolicy::Fixed;
};

// Frozen force PER VOLUME, not traction or integrated force. Cell and native
// staggered contributions are additive, so a physical source must not be
// supplied twice. RZ uses x=r, y=z; the existing kernel supplies all measures.
// Each pair is either present or absent. Values, labels and the source must
// remain unchanged for advance(). No field is retained or copied by this port.
struct HydroForceDensity
{
    const CellField<double>* cellX = nullptr;
    const CellField<double>* cellY = nullptr;
    const VEdgeField<double>* faceX = nullptr;
    const HEdgeField<double>* faceY = nullptr;
    double evaluationTime = 0.0;
    // Required for any supplied field pair. Set by its assembly/publication
    // stage, not reconstructed from the requested step at consumption time.
    std::shared_ptr<const HydroForceSource> source;
    // Default preserves the original frozen-source step. StepEnd is explicit
    // for a prepared conservative interface: drift with u^n, then assemble
    // capillarity on its new C/geometry before the momentum kick. It requires
    // a time-dependent retained source and the same conservative phase step.
    HydroForceTimeLevel timeLevel = HydroForceTimeLevel::StepBegin;
};

// Frozen data for ONE iterate of the compatible midpoint equation. The same
// native midpoint velocity transports the half/full interface and advects
// momentum. With an interface, viscosity is integrated on its actual half-
// step geometry. Without one, the mass source must be single-phase and its
// velocity must match this iterate exactly; owned material fields are frozen.
struct HydroMidpointData
{
    const NativeLeafDomain* domain = nullptr;
    const InterfaceState* interface = nullptr; // Absent for frozen single-phase material.
    const VEdgeField<double>* velocityX = nullptr;
    const HEdgeField<double>* velocityY = nullptr;
    double viscosityNegative = 0.0, viscosityPositive = 0.0;
};

struct HydroStepInput
{
    HydroStepGeometry geometry;
    HydroStepTime time;
    const HydroSolveOptions& options;
    HydroForceDensity force;
    // Current immutable interface authority, borrowed for this call. Null for
    // sources without a retained scalar. Complete-domain mutable phi remains
    // checked by the existing HydroInterfaceLoadGeometry source-validation contract.
    const LevelSetHermiteGeometry::Snapshot* forceInterface = nullptr;
    // Optional conservative temporal mass and first-order implicit upwind
    // advection from this physical Interface step. Closed native leaves or
    // root-grid boundary patches;
    // viscosity remains at step begin; force declares its own quadrature layer.
    const Quantum::Hydrodynamics::HydroConservativeTransport* conservativeTransport = nullptr;
    // Alternative unique geometry authority. It must be the exact beginning
    // or end candidate retained by conservativeTransport at the force layer.
    // A scalar snapshot and parameter geometry cannot be supplied together.
    const InterfaceState* forceParametricInterface = nullptr;
    const HydroMidpointData* midpoint = nullptr;
};

// Identity of successfully completed HYDRO substeps, not acceptance of a
// coupled physical step. AMR transfers this record at the same physical time.
// Geometry is bound by initialization/material input, before completedSteps=1.
struct HydroStepProgress
{
    int completedSteps = 0;
    double endTime = 0.0;
    std::optional<HydroDomainKind> geometry;
    bool stopped = false;
    // prepareStep has advanced the local fields, but the caller has not yet
    // accepted all of its coupled work. advance accepts immediately.
    int acceptedSteps = 0;
    double acceptedEndTime = 0.0;
    bool pendingAcceptance = false;
};

// Move-only acceptance of one locally completed Hydro substep. No field or
// mesh snapshot: read-only Hydro queries expose the candidate time layer.
// Accept after the caller's remaining physical checks, before publishing IO.
// Destruction without accept stops the owner; it does NOT restore old fields.
// A rejected input has no pending candidate and retains the old failure rules.
// The result/report survives owner destruction; accepting it then is invalid.
class HydroStepCandidate
{
public:
    ~HydroStepCandidate();
    HydroStepCandidate(HydroStepCandidate&&) noexcept;
    HydroStepCandidate& operator=(HydroStepCandidate&&) noexcept;
    HydroStepCandidate(const HydroStepCandidate&) = delete;
    HydroStepCandidate& operator=(const HydroStepCandidate&) = delete;
    const HydroSolveReport& report() const;
    bool pending() const noexcept;
    bool accepted() const noexcept;
    void accept(const QuadtreeMesh& mesh);
private:
    friend class Hydrodynamics;
    struct Control;
    explicit HydroStepCandidate(std::shared_ptr<Control> control) noexcept;
    void abandon() noexcept;
    std::shared_ptr<Control> control_;
};

// Optional rollback of ONE managed step on an unchanged mesh. Captures Hydro
// fields, materials and time history before any candidate mutation. A pending
// HydroStepCandidate still exposes the trial fields; commit only after it and
// the caller's remaining physical preparations have succeeded. Destruction
// without commit restores the snapshot. This is not a disk checkpoint, an AMR
// transfer, an automatic retry policy, or a reset of fixed-dt CN/AB2 history.
class HydroStepTransaction
{
public:
    ~HydroStepTransaction();
    HydroStepTransaction(HydroStepTransaction&&) noexcept;
    HydroStepTransaction& operator=(HydroStepTransaction&&) noexcept;
    HydroStepTransaction(const HydroStepTransaction&) = delete;
    HydroStepTransaction& operator=(const HydroStepTransaction&) = delete;
    bool pending() const noexcept;
    void rollback();
    void commit();
private:
    friend class Hydrodynamics;
    struct Control;
    explicit HydroStepTransaction(std::shared_ptr<Control> control) noexcept;
    void abandon() noexcept;
    std::shared_ptr<Control> control_;
};

// Narrow scalar queries; raw Field storage remains private to Hydro.
// No array ownership/copy. Borrow only until the next material update, step,
// mesh transfer, cache restore or destruction of the owning Hydro instance.
template<class FieldType>
class HydroReadOnlyField
{
public:
    explicit HydroReadOnlyField(const FieldType& field) : field_(&field) {}
    int size() const { return field_->size(); }
    bool hasStorage() const { return static_cast<bool>(field_->value); }
    int label(int index) const { return field_->label(index); }
    double operator[](int index) const { return field_->value[index]; }
    // Explicit independent storage for serialization or an existing consumer
    // that needs an owning input. Mutating it never mutates Hydro. Ordinary
    // queries above remain allocation-free synchronous borrows.
    FieldType copy() const { return *field_; }
    template<class Mesh>
    std::unique_ptr<VertField<double>> toVertField(const Mesh& mesh) const
    { return field_->toVertField(mesh); }
private:
    const FieldType* field_;
};

struct HydroStateView
{
    HydroReadOnlyField<CellField<double>> physicalPressure, density, viscosity;
    HydroReadOnlyField<VEdgeField<double>> velocityX, previousVelocityX, predictorX;
    HydroReadOnlyField<HEdgeField<double>> velocityY, previousVelocityY, predictorY;
    HydroReadOnlyField<CellField<double>> divergence, speed;
    HydroReadOnlyField<VertField<double>> vertexDensity, vertexViscosity;
    // Actual assembled cell force used by the latest solve, for its residual
    // diagnostics. These are not the caller's next force input.
    HydroReadOnlyField<CellField<double>> externalForceX, externalForceY;
    HydroStepProgress progress;
    bool initialStateSet = false;
    std::optional<HydroMaterialState> materials;
    // At end of step n: current=U^n, previous=U^(n-1), predictor is scratch
    // from step n, NEVER U^(n-2). Physical p is not the projection multiplier.
};
