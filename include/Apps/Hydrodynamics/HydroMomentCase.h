#pragma once
#include "Modules/Hydrodynamics/HydroCompatibleMomentFlow.h"
#include "Modules/Hydrodynamics/IO/HydroGeometryLoader.h"
#include <memory>
#include <string>
struct HydrodynamicsConfig;
struct InterfaceAdaptationControls;

// Case adapter only. Numerical fields, clock, midpoint and acceptance remain
// in the module owner; no copies in the older staggered Hydro state.
class HydroMomentCase {
public:
    HydroMomentCase(std::unique_ptr<QuadtreeMesh>&,const HydrodynamicsConfig&,
        const std::optional<HydroGeometryInput::ParametricLoop>& = {},const InterfaceAdaptationControls* = nullptr);
    ~HydroMomentCase();
    const Quantum::Hydrodynamics::HydroCompatibleMomentFlow::State& state() const;
    Quantum::Hydrodynamics::HydroCompatibleMomentFlow::Report advance(std::unique_ptr<QuadtreeMesh>&,int step,double begin,double dt);
    void write(const QuadtreeMesh&,const std::string& prefix,bool fields,
        const Quantum::Hydrodynamics::HydroCompatibleMomentFlow::Report* report=nullptr);
private:
    // Zero-time joint transfer, before the step and after each accepted event.
    void adapt(std::unique_ptr<QuadtreeMesh>&,const char* stage);
    struct Impl;
    std::unique_ptr<Impl> impl_;
};
