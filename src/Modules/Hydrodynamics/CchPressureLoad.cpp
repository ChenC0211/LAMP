#include "Modules/Hydrodynamics/CchPressureLoad.h"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace Quantum::Hydrodynamics
{
struct CchPressureLoad::Data
{
    explicit Data(CchIntegralDivergence input) : divergence(std::move(input)) {}
    CchIntegralDivergence divergence;
    Matrix load;
};

CchPressureLoad::CchPressureLoad(std::shared_ptr<const Data> data) : data_(std::move(data)) {}

CchPressureLoad CchPressureLoad::build(CchIntegralDivergence divergence)
{
    auto data = std::make_shared<Data>(std::move(divergence));
    data->load = -data->divergence.matrix().transpose();
    data->load.makeCompressed();
    return CchPressureLoad(std::move(data));
}

const CchIntegralDivergence& CchPressureLoad::divergence() const noexcept { return data_->divergence; }
const CchDofLayout& CchPressureLoad::layout() const noexcept { return data_->divergence.layout(); }
const CchPressureLoad::Matrix& CchPressureLoad::matrix() const noexcept { return data_->load; }

std::vector<double> CchPressureLoad::apply(const QuadtreeMesh& mesh, const CchNativeValues& values) const
{
    values.validateFor(layout(), mesh);
    std::vector<double> result(static_cast<std::size_t>(layout().velocityCount()));
    Eigen::Map<const Eigen::VectorXd> pressure(values.pressure().data(), layout().pressureCount());
    Eigen::Map<Eigen::VectorXd> load(result.data(), layout().velocityCount());
    load.noalias() = data_->load * pressure;
    for (double value : result)
        if (!std::isfinite(value)) throw std::invalid_argument("[CchPressureLoad] nonfinite integrated pressure load");
    return result;
}
}
