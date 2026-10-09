#pragma once

#include "Modules/Electricity/DielectricElectrostatics.h"
#include "Modules/Electricity/SurfaceChargeField.h"

class StaggeredVelocitySampler;

namespace Quantum::Electricity
{
struct DielectricStepTime
{
    int step;
    double forceTime, dt;
};

// Synchronous borrowed inputs from the actual Interface substep. The service
// never advances phi or retains the map/snapshot past completeCharge().
struct DielectricChargeMotion
{
    const QuadtreeMesh& mesh;
    const LevelSetHermiteGeometry::Snapshot& snapshot;
    const std::function<Vector2D(Vector2D)>& forwardMaterialPoint;
};

// Public charge data owns its geometry and integrals. No borrowed evaluator
// or mutable SurfaceChargeField escapes the service. time=completedStep*dt;
// panel measure is physical (ds in XY, revolved area in RZ).
struct DielectricChargeSnapshot
{
    NativeMeshIdentity meshIdentity;
    int completedStep = 0;
    double time = 0.0;
    // Aliasing shared ownership retains the Field's immutable geometry only;
    // no full panel/quad copy per charge update and no mutable charge access.
    std::shared_ptr<const LevelSetInterfacePanels> geometry;
    std::vector<double> integratedCharge, measure;
    explicit DielectricChargeSnapshot(const QuadtreeMesh& mesh) : meshIdentity(mesh) {}
    std::vector<double> density() const;
};
struct DielectricEvaluation
{
    DielectricStepTime time;
    DielectricSolveResult electric;
    // Charge that actually generated this field/current; remains at t^n when
    // the latest charge publication advances to t^(n+1).
    std::shared_ptr<const DielectricChargeSnapshot> charge;
};

struct DielectricChargeUpdate
{
    std::shared_ptr<const DielectricConductionReport> current;
    SurfaceChargeFieldStep transport;
    double electrodeResidual = 0.0;
    DielectricStepTime sourceTime{0, 0.0, 0.0};
    std::shared_ptr<const DielectricChargeSnapshot> charge;
};

// Complete fixed-uniform-domain electric/charge service, using the existing
// diffuse or sharp FV/P3 implementation. No App/Hydro types or file IO.
//
// Preserves ordinary EHD's split: evaluate Q^n on phi^n for force at t^n,
// then complete the explicit charge update using that SAME staged current,
// optionally on phi^(n+1) through the Interface material map. The retained
// electric fields still belong to t^n after completion; they are not a new-time
// electric solve. Fitted/nozzle new-geometry implicit ordering is separate.
// The caller owns the mesh lifetime and sole clock. No mesh clone is made.
// Preflight rejection leaves state usable; once a numerical operation starts,
// failure stops this object permanently. Successful completion is an electric
// substep, not acceptance of the coupled step. No partial state can be read.
class DielectricEvolution
{
public:
    DielectricEvolution(const QuadtreeMesh& mesh, DielectricConfig config, double thickness);
    std::shared_ptr<const DielectricEvaluation> evaluate(const QuadtreeMesh& mesh, const VertField<double>& scalar,
        DielectricStepTime time, const InterfaceGeometry* geometry = nullptr,
        DielectricSharpAssemblyDiagnostics* diagnostics = nullptr,
        std::shared_ptr<const InterfaceTraceGeometry> traceQueries = nullptr);
    DielectricChargeUpdate completeCharge(DielectricStepTime time,
        const StaggeredVelocitySampler* velocity = nullptr,
        const DielectricChargeMotion* motion = nullptr);

    // Immutable, owning results. Previously retained successful publications
    // survive later advances, solver failure and service destruction.
    std::shared_ptr<const DielectricEvaluation> evaluated() const;
    std::shared_ptr<const DielectricChargeSnapshot> charge() const;
    bool failed() const noexcept { return failed_; }
    static double explicitRelaxationCfl(const DielectricConfig& config, double dt);

private:
    void requireUsable() const;
    const QuadtreeMesh* mesh_;
    DielectricConfig config_;
    double thickness_;
    DielectricElectrostatics electric_;
    std::shared_ptr<SurfaceChargeField> charge_;
    std::shared_ptr<const DielectricConductionReport> pendingCurrent_;
    DielectricStepTime evaluatedTime_{0, 0.0, 0.0};
    std::shared_ptr<const DielectricEvaluation> publication_;
    std::shared_ptr<const DielectricChargeSnapshot> chargePublication_;
    bool failed_ = false;
};
}
