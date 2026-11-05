#pragma once

#include "Modules/Hydrodynamics/HydroDofLayout.h"
#include <Eigen/SparseCore>
#include <memory>
#include <vector>

namespace LAMP::Hydrodynamics
{
// Prepared D_I = B A_f on the H01 native layout, before boundary elimination.
// Integral net outflow, not cell-volume-normalized divergence or a time flux.
// Owns only immutable sparse coefficients and a retained layout handle.
class HydroIntegralDivergence
{
public:
    using Matrix = Eigen::SparseMatrix<double, Eigen::RowMajor, Eigen::Index>;
    static HydroIntegralDivergence build(const QuadtreeMesh& mesh, HydroDofLayout layout);
    HydroIntegralDivergence(const HydroIntegralDivergence&) = default;
    HydroIntegralDivergence& operator=(const HydroIntegralDivergence&) = default;

    const HydroDofLayout& layout() const noexcept;
    const Matrix& matrix() const noexcept; // The same D_I consumed by H03.
    std::vector<double> apply(const QuadtreeMesh& mesh, const HydroNativeValues& values) const;

private:
    friend class HydroTwoPhaseState;
    std::vector<double> applyInBatch(const QuadtreeMesh&,const HydroNativeValues&,
        Detail::HydroStateMeshValidationBatch&) const;
    std::vector<double> applyChecked(const HydroNativeValues&) const;
    struct Data;
    explicit HydroIntegralDivergence(std::shared_ptr<const Data> data);
    std::shared_ptr<const Data> data_;
};
}
