#pragma once

#include "Modules/Hydrodynamics/CchDofLayout.h"
#include <Eigen/SparseCore>
#include <memory>
#include <vector>

namespace Quantum::Hydrodynamics
{
// Prepared D_I = B A_f on the H01 native layout, before boundary elimination.
// Integral net outflow, not cell-volume-normalized divergence or a time flux.
// Owns only immutable sparse coefficients and a retained layout handle.
class CchIntegralDivergence
{
public:
    using Matrix = Eigen::SparseMatrix<double, Eigen::RowMajor, Eigen::Index>;
    static CchIntegralDivergence build(const QuadtreeMesh& mesh, CchDofLayout layout);
    CchIntegralDivergence(const CchIntegralDivergence&) = default;
    CchIntegralDivergence& operator=(const CchIntegralDivergence&) = default;

    const CchDofLayout& layout() const noexcept;
    const Matrix& matrix() const noexcept; // The same D_I consumed by H03.
    std::vector<double> apply(const QuadtreeMesh& mesh, const CchNativeValues& values) const;

private:
    struct Data;
    explicit CchIntegralDivergence(std::shared_ptr<const Data> data);
    std::shared_ptr<const Data> data_;
};
}
