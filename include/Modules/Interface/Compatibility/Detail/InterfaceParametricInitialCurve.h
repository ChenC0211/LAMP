#pragma once
#include "Modules/Interface/Compatibility/InterfaceParametricProjection.h"
#include <functional>

namespace InterfacePrescribedData::detail {
using InitialJet=std::array<Quantum::Geometry::Parametric::Point,2>;
using InitialMap=std::function<InitialJet(double)>;
// Shared endpoint ownership and quartic construction for formula initial data.
InterfaceParametricProjection::Curve fromInitialJet(const InitialMap&,int segments,bool axisymmetric);
}
