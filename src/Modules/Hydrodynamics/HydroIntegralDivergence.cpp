#include "Modules/Hydrodynamics/HydroIntegralDivergence.h"

#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>

namespace LAMP::Hydrodynamics
{
namespace
{
void require(bool ok, const char* reason)
{
    if (!ok) throw std::invalid_argument(std::string("[HydroIntegralDivergence] ") + reason);
}
}

struct HydroIntegralDivergence::Data
{
    explicit Data(HydroDofLayout input) : layout(std::move(input)) {}
    HydroDofLayout layout;
    Matrix integral;
};

HydroIntegralDivergence::HydroIntegralDivergence(std::shared_ptr<const Data> data) : data_(std::move(data)) {}

HydroIntegralDivergence HydroIntegralDivergence::build(const QuadtreeMesh& mesh, HydroDofLayout layout)
{
    layout.validateMesh(mesh); // One exact validation before the face-query batch.
    auto data = std::make_shared<Data>(std::move(layout));
    const auto& retained = data->layout;
    std::vector<Eigen::Triplet<double, Eigen::Index>> entries;
    entries.reserve(2 * static_cast<std::size_t>(retained.velocityCount()));
    for (int u = 0; u < retained.velocityCount(); ++u)
    {
        const auto face = retained.domain().face(mesh, retained.nativeFace(u));
        const int negative = face.negativeCell < 0 ? -1 : retained.pressureDof(face.negativeCell);
        const int positive = face.positiveCell < 0 ? -1 : retained.pressureDof(face.positiveCell);
        require((negative >= 0 || positive >= 0) && negative != positive, "invalid selected face incidence");
        require(std::isfinite(face.area) && face.area >= 0, "invalid native face measure");
        if (face.area == 0)
        {
            require(retained.domain().axisymmetric() && !face.horizontal && face.center.x == 0,
                "zero measure is only valid on the RZ axis");
            continue; // Keep the velocity column, with no artificial axis coefficient.
        }
        if (negative >= 0) entries.emplace_back(negative, u, face.area);
        if (positive >= 0) entries.emplace_back(positive, u, -face.area);
    }
    data->integral.resize(retained.pressureCount(), retained.velocityCount());
    data->integral.setFromTriplets(entries.begin(), entries.end());
    data->integral.makeCompressed();
    return HydroIntegralDivergence(std::move(data));
}

const HydroDofLayout& HydroIntegralDivergence::layout() const noexcept { return data_->layout; }
const HydroIntegralDivergence::Matrix& HydroIntegralDivergence::matrix() const noexcept { return data_->integral; }

std::vector<double> HydroIntegralDivergence::apply(const QuadtreeMesh& mesh, const HydroNativeValues& values) const
{
    values.validateFor(data_->layout, mesh);
    return applyChecked(values);
}
std::vector<double> HydroIntegralDivergence::applyInBatch(const QuadtreeMesh& mesh,
    const HydroNativeValues& values,Detail::HydroStateMeshValidationBatch& batch) const
{
    batch.validateValues(mesh,values,data_->layout);
    return applyChecked(values);
}
std::vector<double> HydroIntegralDivergence::applyChecked(const HydroNativeValues& values) const
{
    std::vector<double> result(static_cast<std::size_t>(data_->layout.pressureCount()));
    Eigen::Map<const Eigen::VectorXd> velocity(values.normalVelocity().data(), data_->layout.velocityCount());
    Eigen::Map<Eigen::VectorXd> integral(result.data(), data_->layout.pressureCount());
    integral.noalias() = data_->integral * velocity; // No second divergence loop/stencil.
    for (double value : result) require(std::isfinite(value), "nonfinite integrated divergence");
    return result;
}
}
