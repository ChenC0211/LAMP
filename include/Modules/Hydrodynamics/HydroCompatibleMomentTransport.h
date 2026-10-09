#pragma once
#include "Modules/Hydrodynamics/HydroCompatibleMomentSpace.h"
#include "Modules/Interface/Compatibility/InterfaceParametricTransport.h"
#include <functional>
class InterfaceState;
struct InterfaceStepTime;

namespace Quantum::Hydrodynamics {
class HydroCompatibleMomentMaterial;
class HydroCompatibleMomentBoundary;
// XY/RZ, constant density, native leaves with at most 2:1 split faces.
// Optional port data define incoming momentum on the same native mass flux.
// All reconstructions derive from the ONE owned normal/cell moment state.
// Local RT1 pieces give a normal-continuous field with the exact P1 divergence
// and cell moments. This is H(div), NOT a globally C1 characteristic field.
// RZ uses regular solenoidal streamfunction pieces and the same scaled cell
// moments as mass and stress; its physical reconstruction is defined on null(D).
class HydroCompatibleMomentTransport {
public:
    using Matrix=HydroCompatibleMomentSpace::Matrix;
    enum class Flux {Central,TraceUpwind};
    struct Piece {Vector2D center,size;Eigen::MatrixXd coefficients;
        int degree=1;bool radialWeighted=false,axisRegular=false;double radialScale=1;};
    struct Cell {std::vector<Piece> pieces;double constraintDefect=0,constraintCondition=0,coefficientAmplification=0;
        // RZ: physical L2 Gram of the complete local minimum-energy lift.
        // Supplies the native kinetic Hodge, including its moment complement.
        Eigen::MatrixXd kineticMass;
    };
    struct Assembly {
        Matrix central,dissipation;
        // Budgets of the nonlinear load C(w)w + S(w)w. These are integrated
        // momentum rates (not a separate cell velocity update).
        std::vector<Eigen::Vector2d> cellMomentum,cellBoundaryMomentum;
        double maximumLocalMomentumDefect=0;
        double phaseMassClosure=0,geometricMassClosure=0,massMetricDefect=0;
        Eigen::VectorXd boundaryLoad;
        double boundaryKineticFlux=0,boundaryInflowDissipation=0,boundaryVolumeRate=0;
    };
    struct PhaseMassFlux {
        // Time-integrated mass flux in the SAME native +normal orientation:
        // integrals against (1,theta,theta^2), theta=(s-center)/tangentScale.
        std::vector<Eigen::Vector3d> face;
        std::vector<Eigen::Vector3d> exteriorPositive;
        std::vector<double> massBefore,massAfter;
        // Optional complete velocity-product functionals in the SAME
        // bounded Bernstein/rational face basis as Interface's ledger.
        std::vector<Eigen::VectorXd> products,exteriorPositiveProducts;
        double productConstantDefect=0;
        double maximumRelativeClosure=0;
    };
    static PhaseMassFlux phaseMassFlux(const QuadtreeMesh&,const NativeLeafDomain&,const HydroCompatibleMomentSpace&,
        const InterfaceState& previous,const InterfaceState& next,const InterfaceStepTime&,double liquidDensity,double gasDensity);
    // Native-row volumes from InterfaceParametricMotion::backgroundFaceVolumes.
    // Uses the SAME mass-functional conversion and cell-closure check as the
    // moving two-phase ledger; density is a constant physical material.
    static PhaseMassFlux backgroundMassFlux(const QuadtreeMesh&,const NativeLeafDomain&,const HydroCompatibleMomentSpace&,
        const std::vector<InterfaceParametricTransport::FaceVolumes::Face>&,double density);
    static HydroCompatibleMomentTransport build(const QuadtreeMesh&,const NativeLeafDomain&,const HydroCompatibleMomentSpace&);
    void validate(const QuadtreeMesh&,const NativeLeafDomain&,const HydroCompatibleMomentSpace&) const;
    const std::vector<Cell>& cells() const noexcept{return cells_;}
    static Eigen::Matrix<double,2,12> rtBasis(Vector2D local);
    // Physical velocity from the SAME stored reconstruction. In RZ the
    // reconstructed flux is (r/radialScale)*u; axis pieces factor r explicitly.
    // RZ maps are defined on the local solenoidal moment subspace; assemble
    // rejects non-solenoidal global input before consumption.
    static Eigen::MatrixXd physicalBasis(const Piece&,Vector2D local);
    Assembly assemble(const QuadtreeMesh&,const NativeLeafDomain&,const HydroCompatibleMomentSpace&,
        const Eigen::VectorXd& advectingVelocity,double density,Flux=Flux::Central,
        const HydroCompatibleMomentBoundary* boundary=nullptr) const;
    // Central moving-material transport on native XY/RZ. The actual retained
    // phase flux supplies every face product. Interior flux moments satisfy
    // weak continuity of the same old/new cut-density Gram; the remaining
    // rotational moment uses the midpoint material and velocity projection.
    // No face-density interpolation or post-assembly skew symmetrization.
    Assembly assemblePhase(const QuadtreeMesh&,const NativeLeafDomain&,const HydroCompatibleMomentSpace&,
        const Eigen::VectorXd& midpointVelocity,const InterfaceState& previous,const InterfaceState& next,
        const InterfaceStepTime&,const HydroCompatibleMomentMaterial& before,
        const HydroCompatibleMomentMaterial& midpoint,const HydroCompatibleMomentMaterial& after,
        const HydroCompatibleMomentBoundary* boundary=nullptr) const;
private:
    static PhaseMassFlux massFlux(const QuadtreeMesh&,const NativeLeafDomain&,const HydroCompatibleMomentSpace&,
        const std::function<const InterfaceParametricTransport::FaceVolumes::Face&(int)>& faceVolumes,
        const std::vector<double>& previousLiquid,const std::vector<double>& nextLiquid,
        double liquidDensity,double gasDensity,bool velocityProducts,bool exteriorDirectionMoments);
    std::shared_ptr<const NativeMeshIdentity> identity_;
    int velocityCount_=0;
    bool axisymmetric_=false;
    std::vector<Cell> cells_;
    static HydroCompatibleMomentTransport buildAxisymmetric(const QuadtreeMesh&,const NativeLeafDomain&,const HydroCompatibleMomentSpace&);
};
}
