#pragma once
// S08_FROZEN_BASELINE_BEGIN
// Frozen experimental baseline: CompatibleMoment/cchMomentMidpoint flow.
// Retained only for explicit A/B comparisons and legacy case regression.
// Do not expand its global high-order velocity state or add V2 features here.
// Native staggered CCH V2 must own separate face/cell state; no implicit fallback.
// Scope and source identities: docs/task/S08_Frozen_Baselines.json.
// S08_FROZEN_BASELINE_END
#include "Modules/Hydrodynamics/HydroCompatibleMomentViscosity.h"
#include "Modules/Hydrodynamics/HydroCompatibleMomentTransport.h"
#include "Modules/Hydrodynamics/HydroCompatibleMomentMaterial.h"
#include "Modules/Hydrodynamics/HydroCompatibleMomentSpline.h"
#include "Modules/Interface/Compatibility/InterfaceParametricProjection.h"
#include "Modules/Interface/Compatibility/InterfaceParametricNeck.h"
#include "Modules/Hydrodynamics/HydroBoundary.h"
#include <cstdint>
#include <optional>
#include <functional>

struct HydroStepTime;

namespace Quantum::Hydrodynamics {
class HydroCompatibleMomentBoundary;
// Owner of the new compatible-moment velocity history, distinct from the old
// single-number-per-face Hydro storage. Constant-density NS or frozen-material
// variable-density Stokes, with one owned stress,
// XY/RZ leaves (uniform or 2:1), impermeable no-slip/free-slip walls or native
// normal-velocity / full-stress ports for constant material; the RZ axis has
// zero radial trace and no fictitious face DOF.
// Stokes, constant-density NS and moving two-phase NS use the same saddle
// solver and acceptance. Moving material uses XY/RZ native leaves, with
// closed walls or declared material/momentum reservoirs at physical ports.
class HydroCompatibleMomentFlow {
    struct PhaseData;
    struct Impl;
public:
    struct State {
        // Every flow regime evolves these complete fields. The native
        // moment vectors below are immutable derived publication caches.
        std::optional<HydroCompatibleMomentSpline::NativeVelocity> physicalVelocity,previousPhysicalVelocity,midpointPhysicalVelocity;
        Eigen::VectorXd velocity,previousVelocity,midpointVelocity;
        // Pressure of the last midpoint. The retained viscousTrace member is
        // empty for complete C1 velocity (there is no separate trace state).
        Eigen::VectorXd pressure,viscousTrace;
        int acceptedSteps=0;
        double time=0,previousTime=0,lastDt=0;
        std::optional<double> stageTime;
        // Zero-time remapping preserves the physical stage time, but these
        // published stage fields have not been solved on the new mesh.
        bool remappedStage=false;
        // The two endpoints describe continuous time evolution. An
        // instantaneous topology/material event invalidates extrapolation,
        // while retaining the solved historical fields for diagnostics.
        bool continuousHistory=false;
    };
    struct Report {
        double equationResidual=0,divergenceResidual=0,viscousTraceResidual=0;
        double kineticBefore=0,kineticAfter=0,viscousDissipation=0;
        double externalWork=0,pressureWork=0,energyDefect=0;
        double boundaryTractionWork=0,boundaryReactionWork=0;
        double boundaryKineticFlux=0,boundaryInflowDissipation=0;
        double boundaryVolumeRate=0;
        double boundaryInflowRate=0,boundaryOutflowRate=0;
        double centralConvectionWork=0,transportDissipation=0,localMomentumDefect=0;
        // Arithmetic-midpoint time product rule and cut-geometry/FV mismatch
        // are published separately, never called physical dissipation.
        double temporalMassWork=0,materialMetricWork=0;
        double geometricMassWork=0,boundaryConstraintResidual=0,complementResidual=0;
        double phaseMassClosure=0,geometricMassClosure=0,massMetricDefect=0;
        double capillaryWork=0,surfaceEnergyBefore=0,surfaceEnergyAfter=0,surfaceWorkDefect=0;
        // Independent trajectory integrals, already multiplied by sigma.
        // Positive surfaceAdvection is outward; surfaceEndpointWork is the
        // conormal contribution in weak curvature traction, not a new load.
        double surfaceAdvection=0,surfaceEndpointWork=0,surfaceFluxError=0,surfaceFluxBudget=0,surfaceCriticalError=0;
        int nonlinearIterations=0;
        int activeTangentialDofs=0,inactiveTangentialDofs=0;
        bool reusedFactorization=false;
    };
    struct Options {
        double relativeResidualTolerance=1e-10;
        int maximumIterations=24;
        // Limited-memory acceleration of the SAME midpoint fixed point.
        // Zero retains plain Picard. Acceptance always uses the full equation.
        int accelerationDepth=4;
        HydroCompatibleMomentTransport::Flux convectionFlux=HydroCompatibleMomentTransport::Flux::Central;
        InterfaceParametricProjection::Options interfaceProjection;
    };
    enum class VelocityLevel {Endpoint,PreviousEndpoint,Midpoint};
    struct NormalPolynomial {
        double mean=0,linear=0,tangentCenter=0,tangentScale=0,area=0,time=0;
        double value(double tangentCoordinate) const{return mean+linear*(tangentCoordinate-tangentCenter)/tangentScale;}
        double volumeRate() const{return area*mean;}
    };
    class Load {
    public:
        double time() const noexcept{return time_;}
    private:
        friend class HydroCompatibleMomentFlow;
        std::shared_ptr<const int> owner_;
        Eigen::VectorXd coefficients_;
        Eigen::VectorXd physicalCoefficients_; // Complete [q,eta] virtual work when supplied.
        std::optional<Vector2D> acceleration_; // Reintegrated with the same midpoint phase density.
        double time_=0;
    };
    class Candidate {
    public:
        Candidate(Candidate&&) noexcept=default;
        Candidate& operator=(Candidate&&) noexcept=default;
        Candidate(const Candidate&)=delete;
        Candidate& operator=(const Candidate&)=delete;
        const State& state() const;
        const InterfaceState& interfaceState() const;
        const Report& report() const noexcept{return report_;}
    private:
        friend class HydroCompatibleMomentFlow;
        Candidate()=default;
        std::shared_ptr<const int> owner_;
        std::uint64_t revision_=0;
        std::shared_ptr<const State> state_;
        std::shared_ptr<const PhaseData> phase_;
        Report report_;
        bool fixedPolicy_=false;
    };
    struct TransferLevel {
        double time=0,kineticBefore=0,kineticAfter=0;
        double minimumKinetic=0,unconstrainedKinetic=0,retainedNullspace=1;
        double divergenceResidual=0,momentumResidual=0,correctionNorm=0;
        // Conservative upper estimate of ||u_new-u_old||_rho / ||u_old||_rho,
        // including cancellation roundoff in the quadratic distance.
        double relativeProjectionError=0;
        double constraintResidual=0,kineticCeiling=0;
        int patches=0,independentConstraints=0; // -1: exact zero minimizer, no rank factorization needed.
        bool energyLimited=false;
        bool exactFieldRestriction=false;
    };
    struct MeshTransferReport {
        InterfaceParametricProjection::MeshTransferReport interface;
        std::vector<TransferLevel> levels; // endpoint, previous endpoint, midpoint
    };
    struct NeckReport {
        InterfaceParametricProjection::NeckSplitReport geometry;
        TransferLevel momentum;
        double surfaceBefore=0,surfaceAfter=0,kineticChange=0,totalEnergyChange=0,modelDissipation=0;
        double relativeMassMatrixChange=0;
    };
    class NeckCandidate {
    public:
        NeckCandidate(NeckCandidate&&) noexcept=default;
        NeckCandidate& operator=(NeckCandidate&&) noexcept=default;
        NeckCandidate(const NeckCandidate&)=delete;
        NeckCandidate& operator=(const NeckCandidate&)=delete;
        const NeckReport& report() const noexcept{return report_;}
        const State& state() const;
        const InterfaceState& interfaceState() const;
    private:
        friend class HydroCompatibleMomentFlow;
        NeckCandidate()=default;
        std::shared_ptr<const int> owner_;
        std::uint64_t revision_=0;
        std::shared_ptr<const State> state_;
        std::shared_ptr<const PhaseData> phase_;
        NeckReport report_;
    };
    enum class NeckReadiness {ScalePending,CapPending,MomentumPending,EnergyPending,Ready};
    struct NeckAssessment {
        NeckReadiness readiness=NeckReadiness::ScalePending;
        NeckReport report;
        std::optional<NeckCandidate> candidate;
    };
    class MeshTransferCandidate {
    public:
        ~MeshTransferCandidate();
        MeshTransferCandidate(MeshTransferCandidate&&) noexcept;
        MeshTransferCandidate& operator=(MeshTransferCandidate&&) noexcept;
        MeshTransferCandidate(const MeshTransferCandidate&)=delete;
        MeshTransferCandidate& operator=(const MeshTransferCandidate&)=delete;
        const MeshTransferReport& report() const noexcept{return report_;}
    private:
        friend class HydroCompatibleMomentFlow;
        MeshTransferCandidate();
        std::unique_ptr<Impl> prepared_;
        std::shared_ptr<const int> sourceOwner_;
        std::uint64_t sourceRevision_=0;
        MeshTransferReport report_;
    };
    enum class MeshTransferReadiness {MomentumPending,EnergyPending,AccuracyPending,Ready};
    struct MeshTransferAssessment {
        MeshTransferReadiness readiness=MeshTransferReadiness::MomentumPending;
        int rejectedLevel=-1; // endpoint=0, previous=1, midpoint=2
        MeshTransferReport report;
        std::optional<MeshTransferCandidate> candidate;
    };
    // Cellwise-constant input adapter; stress must retain its physical source
    // from buildCellwise(). General six-moment stress is not a material law:
    // supply createFrozenMaterial() with explicit physical properties instead.
    static HydroCompatibleMomentFlow create(const QuadtreeMesh&,const NativeLeafDomain&,
        HydroCompatibleMomentSpace,HydroCompatibleMomentViscosity,double density,
        Eigen::VectorXd initialVelocity,double initialTime=0,const std::vector<HydroBoundaryPatch>& patches={});
    // With explicit ports, initialization projects the declared initial field
    // in the low-order moment mass norm onto D u=0 and prescribed normal
    // traces, before initializing the complete field. This is not a minimum
    // kinetic-energy projection in the complete velocity space.
    // This zero-time initialization is explicit; subsequent steps never apply
    // a pressure projection or change the pressure/divergence adjoint pair.
    double initialBoundaryCorrection() const;
    // Fixed material field, unsteady Stokes only for nonuniform density.
    // The SAME snapshot builds both inertia and stress. This does not accept
    // moving geometry or advertise a variable-density momentum flux.
    static HydroCompatibleMomentFlow createFrozenMaterial(const QuadtreeMesh&,const NativeLeafDomain&,
        HydroCompatibleMomentSpace,HydroCompatibleMomentMaterial,Eigen::VectorXd initialVelocity,double initialTime=0,
        const std::vector<HydroBoundaryPatch>& patches={});
    // XY/RZ native leaves: one accepted velocity/interface owner. Half/full
    // material and phase flux are recomputed from every Picard midpoint.
    // Central transport and sharp capillarity share the midpoint geometry.
    // External loads are evaluated at the midpoint; electric feedback is not
    // supplied by this owner yet.
    static HydroCompatibleMomentFlow createTwoPhase(const QuadtreeMesh&,const NativeLeafDomain&,
        HydroCompatibleMomentSpace,InterfaceState,HydroCompatibleMomentMaterial::Properties,Eigen::VectorXd initialVelocity,double surfaceTension=0,
        const std::vector<HydroBoundaryPatch>& patches={});
    ~HydroCompatibleMomentFlow();
    HydroCompatibleMomentFlow(HydroCompatibleMomentFlow&&) noexcept;
    HydroCompatibleMomentFlow& operator=(HydroCompatibleMomentFlow&&) noexcept;
    HydroCompatibleMomentFlow(const HydroCompatibleMomentFlow&)=delete;
    HydroCompatibleMomentFlow& operator=(const HydroCompatibleMomentFlow&)=delete;
    const State& state() const;
    const InterfaceState& interfaceState() const;
    const HydroCompatibleMomentMaterial& material() const;
    const HydroCompatibleMomentSpace& space() const;
    const HydroCompatibleMomentViscosity& viscosity() const;
    // The same immutable wall/port partition used by native elimination and
    // boundary work; shared with the common C1 characteristic reconstruction.
    const HydroCompatibleMomentBoundary& boundary() const;
    // Empty for the complete C1 velocity; retained for old observation clients.
    const std::vector<int>& activeTangentialRows() const;
    double density() const;
    // Low-order metric for moment initialization/diagnostics. It is NOT the
    // kinetic metric of the evolved complete field: use kineticEnergy() or
    // physicalCellBalance() to observe physical energy.
    const HydroCompatibleMomentSpace::Matrix& kineticMass() const;
    double kineticEnergy() const;
    // Read-only physical cell integrals of the owned COMPLETE field. Uses
    // the accepted material forms of the requested time layer, without
    // reintegration or replacing them by the old low-order moment metric.
    HydroCompatibleMomentSpline::NativeCellBalance physicalCellBalance(VelocityLevel=VelocityLevel::Endpoint) const;
    struct PhaseIntegral {
        double volume=0,time=0;
        // XY: centroid. RZ: volume-weighted mean radius and axial centroid;
        // mean radius is not the Cartesian centroid of the revolved body.
        Vector2D coordinateMean,velocityIntegral,meanVelocity;
        std::array<double,3> centralSecondMoment{}; // xx, xy, yy (RZ: rr, rz, zz)
        bool present=false;
    };
    // Read only: integrate the actual complete velocity in each phase using
    // its owned cut-material forms. No centroid differencing, density-jump
    // subtraction, new quadrature pass or alteration of the FV inventory.
    std::array<PhaseIntegral,2> phaseIntegrals(const QuadtreeMesh&,const NativeLeafDomain&,
        VelocityLevel=VelocityLevel::Endpoint) const;
    // Actual owned normal polynomial at one declared time layer; volumeRate
    // is the SAME mean used by FV divergence. No auxiliary interpolation or
    // time integral is manufactured by this read-only publication.
    NormalPolynomial normalPolynomial(const QuadtreeMesh&,const NativeLeafDomain&,int nativeFace,
        VelocityLevel=VelocityLevel::Endpoint) const;
    // The immutable load is already integrated against the BULK moment test
    // functions, in force units; no tangential/inertial trace forcing exists.
    Load makeLoad(const QuadtreeMesh&,const NativeLeafDomain&,Eigen::VectorXd integratedForce,double evaluationTime) const;
    // Smooth force per unit PHYSICAL volume, integrated against the actual
    // velocity and pressure-test complement on each native support piece.
    // Uses positive GL12 volume quadrature; discontinuous sharp tractions
    // remain the interface quadrature owner's responsibility.
    Load makeBodyForceLoad(const QuadtreeMesh&,const NativeLeafDomain&,
        const std::function<Vector2D(Vector2D)>& force,double evaluationTime) const;
    // Already integrated physical work at quadrature points, including the
    // physical measure (2*pi*r in RZ). The caller supplies a rule resolving
    // its source discontinuities; both q and pressure tests use the SAME
    // points. No projection to low-order moments or force-density fitting.
    Load makeQuadratureLoad(const QuadtreeMesh&,const NativeLeafDomain&,
        const std::vector<HydroCompatibleMomentSpline::PointForce>&,double evaluationTime) const;
    // Integral of rho*g on the actual midpoint material, against complete
    // velocity and pressure tests. Publishes the total physical pressure.
    Load makeUniformAccelerationLoad(const QuadtreeMesh&,const NativeLeafDomain&,Vector2D acceleration,double evaluationTime) const;
    Candidate prepareStokesStep(const QuadtreeMesh&,const NativeLeafDomain&,
        const HydroStepTime&,const Load&,const Options&);
    Candidate prepareStokesStep(const QuadtreeMesh&,const NativeLeafDomain&,
        const HydroStepTime&,const Load&);
    Candidate prepareNavierStokesStep(const QuadtreeMesh&,const NativeLeafDomain&,
        const HydroStepTime&,const Load&,const Options&);
    Candidate prepareNavierStokesStep(const QuadtreeMesh&,const NativeLeafDomain&,
        const HydroStepTime&,const Load&);
    // Preparation never mutates the accepted fields/clock. Discard is rollback;
    // accept verifies mesh and source revision before an atomic pointer swap.
    void accept(const QuadtreeMesh&,const NativeLeafDomain&,Candidate&&);
    // Zero-time, explicitly modelled topology event. Preserve every native
    // cell's phase inventory and momentum, divergence and prescribed boundary
    // values. Closest velocity in the new kinetic metric with nonincreasing
    // kinetic+surface energy. Historical time layers remain historical.
    NeckCandidate prepareNeckSplit(const QuadtreeMesh&,const NativeLeafDomain&,
        const InterfaceParametricProjection::NeckSplitOptions&) const;
    NeckAssessment assessNeckSplit(const QuadtreeMesh&,const NativeLeafDomain&,
        const InterfaceParametricProjection::NeckSplitOptions&) const;
    void acceptNeckSplit(const QuadtreeMesh&,const NativeLeafDomain&,NeckCandidate&&);
    // Zero-time joint remap retains physical wall/port definitions, all three
    // velocity layers and the absolute pressure stage. It conserves overlap-
    // patch momentum and cannot create kinetic energy; infeasible targets
    // are rejected before acceptance, including inadmissible coarsening.
    MeshTransferCandidate prepareMeshTransfer(const QuadtreeMesh& sourceMesh,const NativeLeafDomain& sourceDomain,
        const QuadtreeMesh& targetMesh,const NativeLeafDomain& targetDomain,const Options&) const;
    // Optional adaptation is declined explicitly for physical infeasibility
    // or excessive information loss. Numerical failures remain errors.
    MeshTransferAssessment assessMeshTransfer(const QuadtreeMesh& sourceMesh,const NativeLeafDomain& sourceDomain,
        const QuadtreeMesh& targetMesh,const NativeLeafDomain& targetDomain,const Options&,double relativeVelocityError) const;
    // All source/revision/target checks precede ownership exchange. Caller
    // keeps the target mesh/domain alive and commits their owning bundle.
    void acceptMeshTransfer(const QuadtreeMesh& sourceMesh,const NativeLeafDomain& sourceDomain,
        const QuadtreeMesh& targetMesh,const NativeLeafDomain& targetDomain,MeshTransferCandidate&&);
private:
    static HydroCompatibleMomentFlow createWithMass(const QuadtreeMesh&,const NativeLeafDomain&,
        HydroCompatibleMomentSpace,HydroCompatibleMomentViscosity,double,HydroCompatibleMomentSpace::Matrix,
        Eigen::VectorXd,double,const std::vector<HydroBoundaryPatch>&,bool projectOpenInitialVelocity=true);
    Candidate prepareStep(const QuadtreeMesh&,const NativeLeafDomain&,
        const HydroStepTime&,const Load&,const Options&,bool convection);
    Candidate prepareNativeStep(const QuadtreeMesh&,const NativeLeafDomain&,
        const HydroStepTime&,const Load&,const Options&,bool convection);
    MeshTransferAssessment prepareNativeMeshTransfer(const QuadtreeMesh&,const NativeLeafDomain&,
        const QuadtreeMesh&,const NativeLeafDomain&,const Options&,double relativeVelocityError) const;
    explicit HydroCompatibleMomentFlow(std::unique_ptr<Impl>);
    std::unique_ptr<Impl> impl_;
};
}
