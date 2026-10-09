#include "Modules/Hydrodynamics/CchIntegralDivergence.h"

#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>

namespace Quantum::Hydrodynamics
{
namespace
{
void require(bool ok, const char* reason)
{
    if (!ok) throw std::invalid_argument(std::string("[CchIntegralDivergence] ") + reason);
}
}

struct CchIntegralDivergence::Data
{
    explicit Data(CchDofLayout input) : layout(std::move(input)) {}
    CchDofLayout layout;
    Matrix integral;
};

CchIntegralDivergence::CchIntegralDivergence(std::shared_ptr<const Data> data) : data_(std::move(data)) {}

CchIntegralDivergence CchIntegralDivergence::build(const QuadtreeMesh& mesh, CchDofLayout layout)
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
    return CchIntegralDivergence(std::move(data));
}

const CchDofLayout& CchIntegralDivergence::layout() const noexcept { return data_->layout; }
const CchIntegralDivergence::Matrix& CchIntegralDivergence::matrix() const noexcept { return data_->integral; }

std::vector<double> CchIntegralDivergence::apply(const QuadtreeMesh& mesh, const CchNativeValues& values) const
{
    values.validateFor(data_->layout, mesh);
    std::vector<double> result(static_cast<std::size_t>(data_->layout.pressureCount()));
    Eigen::Map<const Eigen::VectorXd> velocity(values.normalVelocity().data(), data_->layout.velocityCount());
    Eigen::Map<Eigen::VectorXd> integral(result.data(), data_->layout.pressureCount());
    integral.noalias() = data_->integral * velocity; // No second divergence loop/stencil.
    for (double value : result) require(std::isfinite(value), "nonfinite integrated divergence");
    return result;
}
}
