#pragma once
#include "Core/Geometry/ParametricHermiteCurve.h"
#include <cstdint>
#include <memory>

namespace Quantum::Geometry::Parametric {
// One liquid region, expressed by disjoint material components. Each curve
// is an effective oriented boundary (including its owned domain closure).
// IDs survive smooth transport and remap; a topology operation must explicitly
// assign descendants. Overlap/contact is not silently replaced by a union.
// An empty collection is the absent negative phase, with no free surface.
class Boundary {
public:
    struct Component {std::uint64_t id;std::shared_ptr<const Curve> curve;};
    explicit Boundary(std::vector<Component>);
    const std::vector<Component>& components() const noexcept{return components_;}
    bool matches(const Boundary& b) const noexcept {
        if(components_.size()!=b.components_.size())return false;
        for(std::size_t i=0;i<components_.size();i++)if(components_[i].id!=b.components_[i].id||components_[i].curve!=b.components_[i].curve)return false;return true;
    }
    const std::shared_ptr<const Curve>& single() const {
        check(components_.size()==1,"single-curve consumer requires exactly one active component");return components_.front().curve;
    }
private:
    std::vector<Component> components_;
};
}
