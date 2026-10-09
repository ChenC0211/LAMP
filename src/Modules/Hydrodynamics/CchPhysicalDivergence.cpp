#include "Modules/Hydrodynamics/CchPhysicalDivergence.h"

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
    if (!ok) throw std::invalid_argument(std::string("[CchPhysicalDivergence] ")+reason);
}
}
struct CchPhysicalDivergence::Data
{
    explicit Data(CchIntegralDivergence input) : divergence(std::move(input)) {}
    CchIntegralDivergence divergence;
    std::vector<double> faceMeasures;
    std::vector<int> axis;
};
CchPhysicalDivergence::CchPhysicalDivergence(std::shared_ptr<const Data> data) : data_(std::move(data)) {}
CchPhysicalDivergence CchPhysicalDivergence::build(const QuadtreeMesh& mesh, CchIntegralDivergence divergence)
{
    divergence.layout().validateMesh(mesh);
    auto data=std::make_shared<Data>(std::move(divergence));
    const auto& layout=data->divergence.layout();
    for (double volume : layout.domain().cellVolumes())
        require(std::isfinite(volume) && volume>0,"invalid true cell volume");
    data->faceMeasures.reserve(layout.velocityCount());
    for (int u=0;u<layout.velocityCount();++u) {
        const auto face=layout.domain().face(mesh,layout.nativeFace(u));
        require(std::isfinite(face.area) && face.area>=0,"invalid true face measure");
        data->faceMeasures.push_back(face.area);
        if (face.area==0) {
            require(layout.domain().axisymmetric() && !face.horizontal && face.center.x==0,
                "zero measure outside the radial axis");
            data->axis.push_back(u);
        }
    }
    return CchPhysicalDivergence(std::move(data));
}
const CchIntegralDivergence& CchPhysicalDivergence::divergence() const noexcept {return data_->divergence;}
const CchDofLayout& CchPhysicalDivergence::layout() const noexcept {return data_->divergence.layout();}
const std::vector<double>& CchPhysicalDivergence::cellVolumes() const noexcept {return layout().domain().cellVolumes();}
const std::vector<double>& CchPhysicalDivergence::faceMeasures() const noexcept {return data_->faceMeasures;}
const std::vector<int>& CchPhysicalDivergence::axisVelocityDofs() const noexcept {return data_->axis;}
std::vector<double> CchPhysicalDivergence::apply(const QuadtreeMesh& mesh, const CchNativeValues& values) const
{
    require(values.layout().sameIdentity(layout()),"different native layout identity");
    for (int u : data_->axis)
        require(values.normalVelocity()[u]==0,"nonzero radial axis velocity");
    auto result=data_->divergence.apply(mesh,values); // One actual-mesh validation and the unique D_I.
    const auto& volumes=cellVolumes();
    for (std::size_t p=0;p<result.size();++p) {
        result[p]/=volumes[p];
        require(std::isfinite(result[p]),"nonfinite volume-normalized divergence");
    }
    return result;
}
double CchPhysicalDivergence::reflectAxisNormalVelocity(CchDofLayout::NormalAxis component, double positiveRadiusValue) const
{
    require(!data_->axis.empty(),"axis reflection requires a selected RZ axis");
    require(std::isfinite(positiveRadiusValue),"nonfinite reflected velocity source");
    using Axis=CchDofLayout::NormalAxis;
    require(component==Axis::First || component==Axis::Second,"invalid reflected velocity component");
    return component==Axis::First ? -positiveRadiusValue : positiveRadiusValue;
}
double CchPhysicalDivergence::reflectAxisPressure(double positiveRadiusValue) const
{
    require(!data_->axis.empty(),"axis reflection requires a selected RZ axis");
    require(std::isfinite(positiveRadiusValue),"nonfinite reflected pressure source");
    return positiveRadiusValue;
}
}
