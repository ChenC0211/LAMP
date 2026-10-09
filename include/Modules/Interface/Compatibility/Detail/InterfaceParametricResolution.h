#pragma once
#include "Modules/Interface/Compatibility/InterfaceResolutionPlan.h"
#include <functional>

namespace InterfaceResolutionPlan {
// Exact restrictions of the unique accepted curve provide observations and
// explicit terminal ranges to the common local target/band allocator.
void parametricIndicators(const QuadtreeMesh&,const NativeLeafDomain&,const InterfaceState&,
    const Options&,Result&,std::vector<unsigned char>&,const std::function<void(int,double)>&,double hysteresis=1);
}
