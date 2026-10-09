#pragma once

#include "Modules/Hydrodynamics/Detail/HydroFittedMomentumDuals.h"
#include <Eigen/Sparse>
#include <memory>
#include <vector>

class NativeLeafDomain;
class NativeMeshIdentity;
class QuadtreeMesh;

namespace Quantum::Hydrodynamics
{
// Integrated +grad(p) on the physical staggered duals used by conservative
// mass/momentum transport. Complete XY/RZ native leaves, including 2:1 faces.
// A shared pressure trace supplies opposite loads to neighbouring duals.
// RZ radial loads include -integral(p/r dV), without division at the axis.
// This is an assembly service, not a pressure projection or an NS integrator.
class HydroConservativePressure
{
public:
    using Matrix = Eigen::SparseMatrix<double,Eigen::RowMajor>;
    static HydroConservativePressure build(const QuadtreeMesh&,const NativeLeafDomain&);
    // Explicit CCH pressure port on complete uniform XY/RZ leaves. D uses
    // physical face areas and C=-D^T exactly. Closed boundary normal rows are
    // eliminated by the momentum owner. Mixed-leaf accuracy requires a
    // compatible mass/Hodge construction and is not enabled by this entry.
    static HydroConservativePressure buildDivergenceAdjoint(const QuadtreeMesh&,const NativeLeafDomain&);
    void validate(const QuadtreeMesh&,const NativeLeafDomain&) const;
    const Matrix& gradientLoad() const noexcept { return gradientLoad_; }
    const Matrix& integratedDivergence() const noexcept { return divergence_; }
    bool divergenceAdjoint() const noexcept { return divergenceAdjoint_; }
    const HydroFittedDetail::MomentumDuals& duals() const noexcept { return duals_; }
    // Boundary pressure is an area average, at the weighted point for affine
    // fields. Fixed-velocity rows remain present; the NS caller eliminates them.
    const std::vector<double>& boundaryCoefficient() const noexcept { return boundaryCoefficient_; }
    const std::vector<Vector2D>& boundaryPoints() const noexcept { return boundaryPoints_; }
    // Same matrix, evaluated with scalar differences for an exactly zero
    // constant mode. No thresholding of nonzero scalar gradients.
    std::vector<double> apply(const QuadtreeMesh&,const NativeLeafDomain&,
        const std::vector<double>& scalar,const std::vector<double>& boundaryScalar) const;
private:
    static HydroConservativePressure topology(const QuadtreeMesh&,const NativeLeafDomain&);
    std::shared_ptr<const NativeMeshIdentity> identity_;
    bool axisymmetric_=false;
    bool divergenceAdjoint_=false;
    Matrix gradientLoad_,divergence_;
    HydroFittedDetail::MomentumDuals duals_;
    std::vector<double> boundaryCoefficient_;
    std::vector<Vector2D> boundaryPoints_;
    std::vector<int> anchors_;
};
}
