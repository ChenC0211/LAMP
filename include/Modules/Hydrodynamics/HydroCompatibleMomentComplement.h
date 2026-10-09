#pragma once
#include "Modules/Hydrodynamics/HydroCompatibleMomentSpline.h"

namespace Quantum::Hydrodynamics {
// Continuous physical test directions completing the native solenoidal
// field. Smooth face patches supply mean divergence; cell bubbles supply
// its two linear tests. This owns no evolving velocity or pressure state.
class HydroCompatibleMomentComplement {
public:
    using Matrix=HydroCompatibleMomentSpace::Matrix;
    struct Cell {Vector2D lower,size;std::vector<int> basis;HydroCompatibleMomentSpace::Partition partition;};
    static HydroCompatibleMomentComplement build(const QuadtreeMesh&,const NativeLeafDomain&,
        const HydroCompatibleMomentSpace&,const std::vector<bool>& fixedVelocity);
    void validate(const QuadtreeMesh&,const NativeLeafDomain&,const HydroCompatibleMomentSpace&) const;
    const std::vector<Cell>& cells() const;
    // Local smooth basis -> (u,du/dx,du/dy), six rows in native jet order.
    Eigen::MatrixXd cellJets(int cell,Vector2D) const;
    HydroCompatibleMomentSpline::NativeFaceTrace faceTrace(int cell,
        const InterfaceNormalFluxView::FaceChart&,bool horizontal) const;
    const Matrix& basisMap() const;      // physical local-basis coefficients = E eta
    const Matrix& nativeMoments() const;// R=J_B E; NOT the old arbitrary algebraic R
    const Matrix& divergence() const;   // native D R = G
    const Matrix& pressureLoad() const; // -G^T
    double commutingResidual() const;
    int rootCell() const;
    int outletFace() const;
    HydroCompatibleMomentSpline::FirstJet firstJet(Vector2D,const Eigen::VectorXd& eta) const;
    Eigen::VectorXd transposePointForces(const std::vector<HydroCompatibleMomentSpline::PointForce>&) const;
private:
    struct Impl;
    std::shared_ptr<const Impl> impl_;
};
}
