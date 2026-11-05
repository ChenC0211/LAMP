#pragma once

#include "Modules/Hydrodynamics/HydroIntegralDivergence.h"
#include <memory>
#include <vector>

namespace LAMP::Hydrodynamics
{
// Integral pressure load C = -D_I^T, before any boundary elimination/Hodge.
// Retains the authoritative divergence; never assembles a gradient stencil.
class HydroPressureLoad
{
public:
    using Matrix = HydroIntegralDivergence::Matrix;
    static HydroPressureLoad build(HydroIntegralDivergence divergence);
    HydroPressureLoad(const HydroPressureLoad&) = default;
    HydroPressureLoad& operator=(const HydroPressureLoad&) = default;

    const HydroIntegralDivergence& divergence() const noexcept;
    const HydroDofLayout& layout() const noexcept;
    const Matrix& matrix() const noexcept;
    std::vector<double> apply(const QuadtreeMesh& mesh, const HydroNativeValues& values) const;

private:
    struct Data;
    explicit HydroPressureLoad(std::shared_ptr<const Data> data);
    std::shared_ptr<const Data> data_;
};
}
