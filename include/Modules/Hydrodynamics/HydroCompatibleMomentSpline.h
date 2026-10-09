#pragma once
// S08_FROZEN_BASELINE_BEGIN
// Frozen experimental baseline: CompatibleMoment global C1/Q8 velocity field.
// Retained only for explicit A/B comparisons and legacy case regression.
// This persistent global high-order representation must not become CCH V2 state.
// V2 AMR accuracy belongs to local compatible Hodge/reconstruction; no fallback.
// Scope and source identities: docs/task/S08_Frozen_Baselines.json.
// S08_FROZEN_BASELINE_END
#include "Modules/Hydrodynamics/HydroCompatibleMomentSpace.h"
#include "Modules/Interface/Compatibility/InterfaceTransportInput.h"
#include "Modules/Interface/Compatibility/InterfaceParametricTransport.h"

namespace Quantum::Hydrodynamics {
class HydroCompatibleMomentBoundary;
class HydroCompatibleMomentMaterial;
class HydroCompatibleMomentComplement;
class HydroCompatibleMomentMassCurrent;
// One C1 physical velocity from the existing solenoidal moment state.
// XY/RZ with the owned wall/port partition. Closed uniform XY uses a tensor
// C2 quartic streamfunction.
// Closed uniform RZ uses Psi=r^2*chi and one minimum-kinetic-energy C2 tensor
// quintic space enriched by the complete cell-interior moment space;
// all six RZ cell moments and both face moments are retained.
// The axis has a finite even radial limit, not a no-slip wall. The sparse
// geometric Riesz factor is reused for velocities and the exact work adjoint.
// Mixed leaves and open domains use the same constrained kinetic definition on actual leaf
// Hermite patches; hanging jets are restrictions of their coarse edge.
// Open RZ uses a C2 physical streamfunction Psi; its exact axis factor r^2
// is divided analytically. This also represents nonzero flux constants on
// horizontal walls between ports. Open RZ currently includes the axis.
// No velocity filtering, second evolving field or change to D is performed.
// Its actual weighted normal trace is cubic (uniform XY), up to degree six
// (closed uniform RZ), degree seven/nine on closed mixed XY/RZ leaves, and
// degree seven on open native leaves in either measure.
// It is NOT the old P1 representative. Consumers
// of products of traces must integrate this trace, not just its first moments.
class HydroCompatibleMomentSpline {
    struct Layout;
    struct NativeLayout;
public:
    struct FirstJet {Vector2D value,dx,dy;int nativeCell=-1;};
    // Bilinear forms of the COMPLETE native Q8 field. Global matrices use
    // the scaled Hermite coordinates of this layout, NOT the old moment
    // state. They must not be substituted into its existing NS equations.
    // Positive cut-phase quadrature gives M=V^T W_rho V and K=E^T W_mu E;
    // pressure/divergence and boundary elimination remain separate owners.
    struct NativeMaterialForms {
        struct Cell {
            int index=-1;
            Eigen::MatrixXd mass,viscosity,momentum;
            // Integral of each velocity test over the actual phase region,
            // without density. Index 0 = negative/liquid, 1 = positive/gas.
            std::array<Eigen::MatrixXd,2> phaseVelocityIntegral;
            std::array<double,6> liquid{},gas{};
        };
        std::vector<Cell> cells;
        HydroCompatibleMomentSpace::Matrix mass,viscosity;
        Eigen::MatrixXd momentum;
        std::array<Eigen::MatrixXd,2> phaseVelocityIntegral;
        int solenoidalCoordinates=0,complementCoordinates=0;
        double kinetic=0,dissipation=0;
        Vector2D physicalMomentum{0,0};
        double maximumPhaseMomentError=0,quadratureEstimate=0;
        std::size_t quadraturePoints=0;
        int reusedPureCells=0,preparedNormCells=0;
    };
    // Integrated work on value, dx and dy of the SAME first jet.
    struct PointForce {Vector2D point,measure,dx{0,0},dy{0,0};};
    struct NativeFaceProducts {
        InterfaceParametricTransport::FaceVolumes::Products x{},y{},speedSquared{};
    };
    struct NativeFaceTrace {
        // Each column is a physical test velocity on one native face:
        // numerator in degree-10 Bernstein / d, with d=r/r_min only on
        // off-axis horizontal RZ faces. Local coordinates, no fitting.
        Eigen::Matrix<double,11,Eigen::Dynamic> x,y;
    };
    struct NativeCellBalance {Eigen::MatrixXd momentum;Eigen::VectorXd kinetic;};
    struct NativeProjectionForms {
        NativeMaterialForms material;
        Eigen::VectorXd load; // integral rho_target * v_target * u_source
        double sourceNormSquared=0; // SAME target material, for projection distance
    };
    struct NativeVolumeTransport {
        // TIME-INTEGRATED volume term -integral j_dt.grad(v_test) v_trial.
        // Native face momentum and prescribed incoming data remain separate.
        HydroCompatibleMomentSpace::Matrix volume;
        // Observable weak-continuity defect from the geometric constant
        // mass mismatch. NEVER added to inertia or momentum. Its coefficient
        // is fixed by the constant-mode functional of the scalar space.
        HydroCompatibleMomentSpace::Matrix geometricMetric;
        int solenoidalCoordinates=0,complementCoordinates=0;
    };
    struct NativeMomentumTransport : NativeVolumeTransport {
        // Time-integrated central volume + exterior outflow operator.
        // The source for prescribed incoming momentum is on the RHS.
        HydroCompatibleMomentSpace::Matrix central,boundaryMetric;
        Eigen::VectorXd boundaryLoad;
        double boundaryKineticFlux=0,boundaryInflowDissipation=0,boundaryVolume=0;
        double traceDefect=0;
        int pressureBackflowFaces=0,bidirectionalFaces=0;
    };
    // Complete kinematic velocity: coordinates are the sole independent
    // data; local patches and native moments are derived caches. This is a
    // building block for the SAME Hydro owner, not a separate time solver.
    // Wall/axis conditions are checked here. Nonzero prescribed inlet data
    // and momentum/pressure acceptance remain the physical owner's job.
    class NativeVelocity {
    public:
        void validate(const QuadtreeMesh&,const NativeLeafDomain&,const HydroCompatibleMomentSpace&) const;
        const Eigen::VectorXd& coordinates() const noexcept{return coordinates_;}
        const Eigen::VectorXd& nativeMoments() const noexcept{return moments_;}
        const HydroCompatibleMomentSpace::Matrix& momentMap() const;
        const HydroCompatibleMomentSpace::Matrix& kinematicConstraints() const;
        // Owned cell-interior Hermite coordinates. Exact local elimination
        // may use these blocks; the recovered complete q remains the state.
        const std::vector<std::array<int,9>>& cellInteriorCoordinates() const;
        HydroCompatibleMomentSpace::Matrix normalBoundaryConstraints(const HydroCompatibleMomentBoundary&) const;
        NativeVelocity withCoordinates(const Eigen::VectorXd&,double time) const;
        FirstJet firstJet(Vector2D) const;
        struct CellSpeedSamples {int cell;double coarseMaximum=0,fineMaximum=0;};
        // Read-only nested 5x5/9x9 samples of the COMPLETE field on every
        // native leaf (including edges). Sample maxima, not L-infinity bounds.
        std::vector<CellSpeedSamples> cellSpeedSamples() const;
        // Kinematic continuation for the SAME material interface. Outside a
        // declared pressure outlet, extend the boundary streamfunction by
        // its quadratic normal Taylor polynomial, retaining velocity/C1
        // matching and divergence structure. It owns no exterior NS state.
        FirstJet transportJet(Vector2D) const;
        Vector2D transportLower() const;
        Vector2D transportUpper() const;
        InterfaceNormalFluxView normalFlux() const;
        // Exact algebraic trace coefficients for work against the SAME
        // retained Bernstein/rational phase-flux functionals. No fitting.
        NativeFaceProducts faceProducts(int nativeFace) const;
        NativeVolumeTransport nativeVolumeTransport(const QuadtreeMesh&,const NativeLeafDomain&,
            const HydroCompatibleMomentSpace&,const HydroCompatibleMomentMassCurrent&,
            const HydroCompatibleMomentComplement* complement=nullptr) const;
        NativeMomentumTransport nativeMomentumTransport(const QuadtreeMesh&,const NativeLeafDomain&,
            const HydroCompatibleMomentSpace&,const HydroCompatibleMomentMassCurrent&,
            const HydroCompatibleMomentBoundary&,const HydroCompatibleMomentComplement* complement=nullptr) const;
        // Work-conjugate load in COMPLETE coordinates, not old native U.
        Eigen::VectorXd transposePointForces(const std::vector<PointForce>&) const;
        // Gram of the individual work covectors in the SAME complete q space.
        // Positive quadrature factors are included as square roots in each
        // covector. This is a solver metric, not a new physical velocity.
        HydroCompatibleMomentSpace::Matrix pointWorkGram(const std::vector<PointForce>&) const;
        NativeMaterialForms nativeMaterialForms(const QuadtreeMesh&,const NativeLeafDomain&,
            const HydroCompatibleMomentSpace&,const HydroCompatibleMomentMaterial&,double tolerance=2e-12) const;
        // The same positive material quadrature, including ALL cross terms
        // with the smooth pressure-test complement. Coordinates are [q,eta].
        NativeMaterialForms mixedMaterialForms(const QuadtreeMesh&,const NativeLeafDomain&,
            const HydroCompatibleMomentSpace&,const HydroCompatibleMomentMaterial&,
            const HydroCompatibleMomentComplement&,double tolerance=2e-12) const;
        double time() const noexcept{return time_;}
        Vector2D lower() const;
        Vector2D upper() const;
        bool sameLayout(const NativeVelocity&) const;
        bool openBoundary() const;
        bool usesGeneralTraceSpace() const;
        NativeVelocity refine(const QuadtreeMesh&,const NativeLeafDomain&,const HydroCompatibleMomentSpace&,
            const HydroCompatibleMomentBoundary* boundary=nullptr) const;
        // Reuse this target layout for another historical source field.
        NativeVelocity restrictFrom(const NativeVelocity& source) const;
        NativeCellBalance cellBalance(const NativeMaterialForms&) const;
        HydroCompatibleMomentSpace::Matrix cellMomentumMap(const NativeMaterialForms&) const;
        // Overlay-resolved physical L2 pairing, using the same positive
        // cut-material quadrature as the momentum equation. No moment lift.
        NativeProjectionForms projectionForms(const QuadtreeMesh&,const NativeLeafDomain&,
            const HydroCompatibleMomentSpace&,const HydroCompatibleMomentMaterial&,
            const NativeVelocity& source,const HydroCompatibleMomentComplement* complement=nullptr,double tolerance=2e-12) const;
    private:
        friend class HydroCompatibleMomentSpline;
        NativeVelocity()=default;
        std::shared_ptr<const NativeLayout> layout_;
        Eigen::VectorXd coordinates_,moments_;
        Eigen::MatrixXd coefficients_;
        double time_=0;
    };
    // Transient exact image of the accepted Q8 field on refined native
    // leaves. It carries the retained polynomial detail, not just its native
    // moments. It is NOT an accepted Hydro state and deliberately has no
    // withVelocity or target-moment work transpose: those would discard detail.
    class RefinementSource {
    public:
        FirstJet firstJet(Vector2D) const;
        InterfaceNormalFluxView normalFlux() const;
        NativeMaterialForms nativeMaterialForms(const QuadtreeMesh&,const NativeLeafDomain&,
            const HydroCompatibleMomentSpace&,const HydroCompatibleMomentMaterial&,double tolerance=2e-12) const;
        const Eigen::VectorXd& nativeMoments() const noexcept{return moments_;}
        double time() const noexcept{return time_;}
        NativeVelocity nativeVelocity() const;
    private:
        friend class HydroCompatibleMomentSpline;
        RefinementSource()=default;
        std::shared_ptr<const NativeLayout> layout_;
        Eigen::MatrixXd coefficients_;
        Eigen::VectorXd moments_;
        double time_=0;
    };
    static HydroCompatibleMomentSpline build(const QuadtreeMesh&,const NativeLeafDomain&,
        const HydroCompatibleMomentSpace&,const Eigen::VectorXd&,double time,
        const HydroCompatibleMomentBoundary* boundary=nullptr);
    // Initial moment-to-field choice only. Subsequent evolution keeps q;
    // uniform closed meshes use this same complete native space as AMR.
    static NativeVelocity initializeNative(const QuadtreeMesh&,const NativeLeafDomain&,
        const HydroCompatibleMomentSpace&,const Eigen::VectorXd&,double time,
        const HydroCompatibleMomentBoundary* boundary=nullptr);
    // The immutable partition fixes the reconstruction SPACE, not its input
    // coefficients. withVelocity and the work transpose retain this space
    // even when an inlet moment is zero or changes sign. Port values remain
    // the native solver's responsibility; all supplied moments are retained.
    HydroCompatibleMomentSpline withVelocity(const QuadtreeMesh&,const NativeLeafDomain&,
        const HydroCompatibleMomentSpace&,const Eigen::VectorXd&,double time) const;
    // Complete native Q8 source layouts (open XY/RZ or mixed leaves).
    // Every target leaf must lie in one source leaf; coarsening is rejected.
    // This is algebraic restriction, with no fit or kinetic minimization.
    RefinementSource prepareNativeRefinement(const QuadtreeMesh&,const NativeLeafDomain&,
        const HydroCompatibleMomentSpace&,const HydroCompatibleMomentBoundary* boundary=nullptr) const;
    NativeMaterialForms nativeMaterialForms(const QuadtreeMesh&,const NativeLeafDomain&,
        const HydroCompatibleMomentSpace&,const HydroCompatibleMomentMaterial&,double tolerance=2e-12) const;
    NativeVelocity nativeVelocity() const;
    FirstJet firstJet(Vector2D) const;
    InterfaceNormalFluxView normalFlux() const;
    // Exact transpose of the same linear reconstruction, integrated point
    // forces in physical units. It defines work on the solenoidal state;
    // its extension off null(D) is an algebraic pressure-gauge choice.
    Eigen::VectorXd transposePointForces(const std::vector<PointForce>&) const;
    Vector2D lower() const;
    Vector2D upper() const;
    bool sameLayout(const HydroCompatibleMomentSpline&) const;
    bool openBoundary() const;
    bool usesGeneralTraceSpace() const;
    double time() const noexcept{return time_;}
private:
    std::shared_ptr<const Layout> layout_;
    Eigen::MatrixXd coefficients_;
    Eigen::VectorXd interior_;
    double time_=0;
    static HydroCompatibleMomentSpline bind(std::shared_ptr<const Layout>,const HydroCompatibleMomentSpace&,
        const Eigen::VectorXd&,double);
    static NativeVelocity bindNative(std::shared_ptr<const NativeLayout>,const Eigen::VectorXd&,double);
};
}
