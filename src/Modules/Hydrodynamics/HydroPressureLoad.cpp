#include "Modules/Hydrodynamics/HydroPressureLoad.h"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace LAMP::Hydrodynamics
{
struct HydroPressureLoad::Data
{
    explicit Data(HydroIntegralDivergence input) : divergence(std::move(input)) {}
    HydroIntegralDivergence divergence;
    Matrix load;
};

HydroPressureLoad::HydroPressureLoad(std::shared_ptr<const Data> data) : data_(std::move(data)) {}

HydroPressureLoad HydroPressureLoad::build(HydroIntegralDivergence divergence)
{
    auto data = std::make_shared<Data>(std::move(divergence));
    data->load = -data->divergence.matrix().transpose();
    data->load.makeCompressed();
    return HydroPressureLoad(std::move(data));
}

const HydroIntegralDivergence& HydroPressureLoad::divergence() const noexcept { return data_->divergence; }
const HydroDofLayout& HydroPressureLoad::layout() const noexcept { return data_->divergence.layout(); }
const HydroPressureLoad::Matrix& HydroPressureLoad::matrix() const noexcept { return data_->load; }

std::vector<double> HydroPressureLoad::apply(const QuadtreeMesh& mesh, const HydroNativeValues& values) const
{
    values.validateFor(layout(), mesh);
    std::vector<double> result(static_cast<std::size_t>(layout().velocityCount()));
    Eigen::Map<const Eigen::VectorXd> pressure(values.pressure().data(), layout().pressureCount());
    Eigen::Map<Eigen::VectorXd> load(result.data(), layout().velocityCount());
    load.noalias() = data_->load * pressure;
    for (double value : result)
        if (!std::isfinite(value)) throw std::invalid_argument("[HydroPressureLoad] nonfinite integrated pressure load");
    return result;
}
}
