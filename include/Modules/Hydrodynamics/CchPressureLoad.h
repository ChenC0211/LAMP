#pragma once

#include "Modules/Hydrodynamics/CchIntegralDivergence.h"
#include <memory>
#include <vector>

namespace Quantum::Hydrodynamics
{
// Integral pressure load C = -D_I^T, before any boundary elimination/Hodge.
// Retains the authoritative divergence; never assembles a gradient stencil.
class CchPressureLoad
{
public:
    using Matrix = CchIntegralDivergence::Matrix;
    static CchPressureLoad build(CchIntegralDivergence divergence);
    CchPressureLoad(const CchPressureLoad&) = default;
    CchPressureLoad& operator=(const CchPressureLoad&) = default;

    const CchIntegralDivergence& divergence() const noexcept;
    const CchDofLayout& layout() const noexcept;
    const Matrix& matrix() const noexcept;
    std::vector<double> apply(const QuadtreeMesh& mesh, const CchNativeValues& values) const;

private:
    struct Data;
    explicit CchPressureLoad(std::shared_ptr<const Data> data);
    std::shared_ptr<const Data> data_;
};
}
