#pragma once

#include "Modules/Hydrodynamics/CchIntegralDivergence.h"
#include <memory>
#include <vector>

namespace Quantum::Hydrodynamics
{
// Physical cell-average divergence from the SAME integral D_I and true V_c.
// Axis reflections are derived boundary values, never ghost primary DOFs.
class CchPhysicalDivergence
{
public:
    static CchPhysicalDivergence build(const QuadtreeMesh& mesh, CchIntegralDivergence divergence);
    const CchIntegralDivergence& divergence() const noexcept;
    const CchDofLayout& layout() const noexcept;
    const std::vector<double>& cellVolumes() const noexcept;
    const std::vector<double>& faceMeasures() const noexcept;
    const std::vector<int>& axisVelocityDofs() const noexcept;
    std::vector<double> apply(const QuadtreeMesh& mesh, const CchNativeValues& values) const;

    // Caller supplies a sample at positive radius; return its negative-radius
    // partner. Only a selected RZ domain retaining the axis supports this BC.
    double reflectAxisNormalVelocity(CchDofLayout::NormalAxis component, double positiveRadiusValue) const;
    double reflectAxisPressure(double positiveRadiusValue) const;
private:
    struct Data;
    explicit CchPhysicalDivergence(std::shared_ptr<const Data> data);
    std::shared_ptr<const Data> data_;
};
}
