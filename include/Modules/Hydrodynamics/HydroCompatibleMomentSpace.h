#pragma once
#include "Core/Base/Vector2D.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <memory>
#include <vector>
#include <functional>
#include <array>

class QuadtreeMesh;
class NativeLeafDomain;
class NativeMeshIdentity;

namespace Quantum::Hydrodynamics {
// Compatible weighted moment space on the EXISTING complete Cartesian leaves.
// This explicit research/assembly port is not the old native-face NS state.
// Each positive-area face has the mean normal velocity and its coefficient
// against (t-weightedCentroid)/weightedStdDev. Each cell has weighted mean
// ux,uy and mean(u dot (-Y,X)); X,Y are local coordinates divided by cell size.
// RZ additionally owns three scaled streamfunction moments. With q=(r/rc)u,
// div(q)=d0+dX*X+dY*Y, subtract qp=(integral_{r_left}^r div(q) dr,0).
// Its solenoidal remainder is (psi_z,-psi_r), psi(bottom left)=0. Extra DOFs
// are mean_A(psi*(X*s,X^2*s^2,Y*s^2))/sqrt(hx*hy), s=rc/r. They replace
// nearly dependent flat moments without changing the retained information.
// Pressure has coefficients in (1,X,Y). D integrates div(w*u) against these
// three tests, with w=1 (XY) or 2*pi*r (RZ); C=-D^T by construction.
// The local virtual space retains div(w*u) in P1. Its mass inner product
// reproduces ALL affine physical velocities, including the rotational mode.
// RZ uses the physical L2 Gram of its constrained local minimum-energy lift,
// including the divergence particular solution; raw stream-moment coordinate
// sizes do not determine kinetic penalties. This is the same native state.
// Zero-area axis faces have no flux unknown. Axis regularity of an NS state,
// material transport, stress, nonlinear stepping and AMR transfer are separate
// consumers and are not enabled merely by constructing this metric.
class HydroCompatibleMomentSpace {
public:
    using Matrix=Eigen::SparseMatrix<double>;
    struct Face {
        int first=-1; // mean, then centered linear coefficient; -1 on RZ axis
        double area=0,tangentCenter=0,tangentScale=0;
    };
    struct CellBlock {
        std::vector<int> rows;
        Eigen::MatrixXd mass;
        // Local velocity coefficients -> coefficients of an affine physical
        // velocity in (grad(h0*X),grad(h0*Y),grad(h0*X^2/2),
        // grad(h0*X*Y),grad(h0*Y^2/2),(-Y,X)), h0=sqrt(hx*hy).
        Eigen::MatrixXd affineProjection;
        Eigen::MatrixXd affineSamples; // Same six physical fields -> owned moments.
        // RZ only: flat means (ur,X*ur,Y*ur,uz,X*uz,Y*uz) from local DOFs.
        // A derived kinematic functional, not an additional evolving field.
        Eigen::MatrixXd flatVelocityMoments;
    };
    struct Report {
        double maximumConsistencyDefect=0,maximumSymmetryDefect=0;
        // XY: original baseline scaling; RZ: unit diagonal physical Hodge.
        // Different moment coordinates need not have the same spectrum.
        double minimumScaledEigenvalue=0,maximumScaledEigenvalue=0;
        int maximumLocalVelocityDofs=0;
    };
    // Algebraic complement to a complete solenoidal velocity. R maps its
    // coordinates to native moments and D*R=G. G has full column rank;
    // Mean-divergence coordinates use local spanning-tree edges (and one
    // outlet for an open domain), rather than root-to-cell path sums.
    // Closed domains exclude the one net-flux direction. The pressure rows
    // remain the original native P1 rows, with load -G^T. This is NOT a
    // physical reconstruction, time projection or final NS operator.
    struct DivergenceComplement {
        Matrix moments,divergence,pressureLoad;
        int rootCell=-1,outletFace=-1;
        double commutingResidual=0;
    };
    static HydroCompatibleMomentSpace build(const QuadtreeMesh&,const NativeLeafDomain&);
    // Same affine basis for physical reconstruction, loads and stress consumers.
    // local=(x-center.x)/hx,(y-center.y)/hy.
    static Eigen::Matrix<double,2,6> affineVelocityBasis(Vector2D local,Vector2D cellSize);
    void validate(const QuadtreeMesh&,const NativeLeafDomain&) const;
    DivergenceComplement divergenceComplement(const QuadtreeMesh&,const NativeLeafDomain&,
        const std::vector<bool>& fixedVelocity) const;
    // Physical input ports, integrated with eight-point Gauss rules. They do
    // not project divergence or alter boundary data; the owner validates the
    // resulting initial state. Body force is force per physical volume.
    using PhysicalVector=std::function<Eigen::Vector2d(Vector2D)>;
    // Shared DIRECT evaluation of stable coordinates for one or several
    // physical fields (columns). The particular flux is zero on the left
    // boundary; the canonical streamfunction is therefore obtained directly
    // from the physical field by a left-edge then horizontal path integral.
    using PhysicalFields=std::function<Eigen::MatrixXd(Vector2D)>;
    static Eigen::MatrixXd rzExtraMoments(Vector2D center,Vector2D size,
        const PhysicalFields&);
    using Partition=std::array<std::vector<double>,2>;
    static Eigen::MatrixXd rzExtraMoments(Vector2D center,Vector2D size,
        const PhysicalFields&,const Partition&);
    Eigen::VectorXd interpolateVelocity(const QuadtreeMesh&,const NativeLeafDomain&,const PhysicalVector&) const;
    // Same moment functional, split at the donor's actual polynomial planes.
    // Required when coarsening over multiple old patches; no sampled remesh.
    Eigen::VectorXd interpolateVelocity(const QuadtreeMesh&,const NativeLeafDomain&,const PhysicalVector&,const Partition&) const;
    Eigen::VectorXd integrateBodyForce(const QuadtreeMesh&,const NativeLeafDomain&,const PhysicalVector&) const;
    const Matrix& mass() const noexcept{return mass_;}
    const Matrix& divergence() const noexcept{return divergence_;}
    const Matrix& pressureLoad() const noexcept{return pressureLoad_;}
    const std::vector<Face>& faces() const noexcept{return faces_;}
    const std::vector<CellBlock>& cells() const noexcept{return cells_;}
    int cellVelocityCount() const noexcept{return axisymmetric_?6:3;}
    int firstCellVelocity(int cell) const noexcept{return cellVelocityBegin_+cellVelocityCount()*cell;}
    const Report& report() const noexcept{return report_;}
private:
    std::shared_ptr<const NativeMeshIdentity> identity_;
    bool axisymmetric_=false;
    int cellVelocityBegin_=0;
    Matrix mass_,divergence_,pressureLoad_;
    std::vector<Face> faces_;
    std::vector<CellBlock> cells_;
    Report report_;
    Eigen::MatrixXd rzFlatMoments(const QuadtreeMesh&,const NativeLeafDomain&,int cell,
        const std::vector<int>& nativeFaces,const Eigen::MatrixXd& localDivergence) const;
};
}
