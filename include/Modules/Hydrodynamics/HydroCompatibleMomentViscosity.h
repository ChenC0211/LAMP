#pragma once
#include "Modules/Hydrodynamics/HydroCompatibleMomentSpace.h"

namespace Quantum::Hydrodynamics {
// Weak symmetric strain on the existing compatible moment velocity space.
// The bulk/normal unknowns retain their original indices. Two shared tangential
// coefficients per native face follow them. These are algebraic traces, with
// NO artificial face mass or second physical velocity field.
// XY/RZ complete Cartesian leaves, including 2:1 split sides. RZ uses the
// owned streamfunction moments and a radial reconstruction vanishing on-axis.
// The matrix is the Hessian of one nonnegative dissipation; no pressure terms
// or boundary conditions are imposed here. The material owner supplies
// integral_K mu*(1,X,Y)*(1,X,Y)^T: six moments of the SAME viscosity measure.
class HydroCompatibleMomentViscosity {
public:
    using Matrix=HydroCompatibleMomentSpace::Matrix;
    struct CellBlock {
        std::vector<int> rows;
        // XY: (exx, eyy, sqrt(2)*exy); RZ: (err,ezz,ett,sqrt(2)*erz).
        // Each component is in local (1,X,Y).
        Eigen::MatrixXd strain;
        // Quadratic viscous reconstruction (ux then uy, each in
        // 1,X,Y,X^2,XY,Y^2). Same strain and the three owned cell moments;
        // no independent velocity state is introduced.
        Eigen::MatrixXd reconstruction;
        Eigen::MatrixXd dissipation;
    };
    static HydroCompatibleMomentViscosity build(const QuadtreeMesh&,const NativeLeafDomain&,
        const HydroCompatibleMomentSpace&,const std::vector<Eigen::Matrix3d>& viscosityMoments);
    // Explicit piecewise constant material, not a high-order approximation to
    // a varying viscosity. It uses the same assembly and preserves positivity.
    static HydroCompatibleMomentViscosity buildCellwise(const QuadtreeMesh&,const NativeLeafDomain&,
        const HydroCompatibleMomentSpace&,const std::vector<double>& cellViscosity);
    void validate(const QuadtreeMesh&,const NativeLeafDomain&,const HydroCompatibleMomentSpace&) const;
    const Matrix& matrix() const noexcept{return matrix_;}
    // Available only when the caller supplied a physical cellwise law.
    // Six general viscosity moments cannot define a complete velocity form.
    const std::shared_ptr<const std::vector<double>>& cellViscosity() const noexcept{return cellViscosity_;}
    const std::vector<CellBlock>& cells() const noexcept{return cells_;}
    int firstTangent(int face) const noexcept{return tangentRows_[face];} // -1 on axis
    int bulkSize() const noexcept{return bulkSize_;}
private:
    std::shared_ptr<const NativeMeshIdentity> identity_;
    std::shared_ptr<const std::vector<double>> cellViscosity_;
    int bulkSize_=0;
    std::vector<int> tangentRows_;
    Matrix matrix_;
    std::vector<CellBlock> cells_;
    static HydroCompatibleMomentViscosity buildAxisymmetric(const QuadtreeMesh&,const NativeLeafDomain&,
        const HydroCompatibleMomentSpace&,const std::vector<Eigen::Matrix3d>&);
};
}
