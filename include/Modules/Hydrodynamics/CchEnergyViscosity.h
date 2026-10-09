#pragma once

#include "Modules/Hydrodynamics/CchPhysicalDivergence.h"
#include <Eigen/SparseCore>
#include <memory>
#include <vector>

namespace Quantum::Hydrodynamics
{
// A volume quadrature sample, not a material/phase classifier. Physical weight
// includes 2*pi*r in RZ. H07 must supply cut samples from the unique geometry.
struct CchStrainPoint
{
    int cellDof = -1;
    Vector2D point{0,0};
    double physicalWeight = 0;
};

// Derived local P2 coefficients are eliminated into this native-velocity map.
// Rows per point: XY (xx,yy,sqrt(2)xy), RZ (rr,zz,hoop,sqrt(2)rz).
class CchStrainOperator
{
public:
    using Matrix = Eigen::SparseMatrix<double,Eigen::RowMajor,Eigen::Index>;
    static std::vector<CchStrainPoint> volumeQuadrature(const QuadtreeMesh& mesh,
        const CchPhysicalDivergence& metrics);
    static CchStrainOperator build(const QuadtreeMesh& mesh, CchPhysicalDivergence metrics,
        const std::vector<CchStrainPoint>& points);
    const CchPhysicalDivergence& metrics() const noexcept;
    const CchDofLayout& layout() const noexcept;
    const std::vector<CchStrainPoint>& points() const noexcept;
    int components() const noexcept;
    const Matrix& matrix() const noexcept;
    std::vector<double> apply(const QuadtreeMesh& mesh, const CchNativeValues& values) const;
private:
    struct Data;
    explicit CchStrainOperator(std::shared_ptr<const Data> data);
    std::shared_ptr<const Data> data_;
};

// K_mu = E^T diag(2*mu_q*dV_q) E; momentum uses -K_mu*u.
// No minimum viscosity, boundary penalty, phase averaging or new primary DOF.
class CchEnergyViscosity
{
public:
    using Matrix = CchStrainOperator::Matrix;
    static CchEnergyViscosity build(CchStrainOperator strain, const std::vector<double>& viscosity);
    const CchStrainOperator& strain() const noexcept;
    const std::vector<double>& viscosity() const noexcept;
    const std::vector<double>& rowWeights() const noexcept;
    const Matrix& matrix() const noexcept;
    std::vector<double> apply(const QuadtreeMesh& mesh, const CchNativeValues& values) const;
    double dissipation(const QuadtreeMesh& mesh, const CchNativeValues& values) const;
private:
    struct Data;
    explicit CchEnergyViscosity(std::shared_ptr<const Data> data);
    std::shared_ptr<const Data> data_;
};
}
