#pragma once
#include <functional>
#include <map>
#include <string>

namespace Quantum {
// Immutable arithmetic expression for analytic input data.
// No scripting, assignments, files or runtime dictionary access. Derivatives
// use the chain rule on the same expression as the value.
class ScalarExpression {
public:
    struct Jet {double value=0,derivative=0,partialSecond=0;}; // First partials, not a Hessian.
    enum class Coordinates {Parameter,Planar,Meridional,Time,PlanarTime,MeridionalTime};
    using Constants=std::map<std::string,double>;
    ScalarExpression(const std::string& formula,const Constants& constants={},Coordinates coordinates=Coordinates::Parameter);
    Jet sample(double s) const;
    Jet sample(double first,double second) const;
    Jet sample(double first,double second,double time) const; // Spatial partials at fixed time.
private:
    Coordinates coordinates_;
    std::function<Jet(double,double,double)> evaluate_;
};
}
