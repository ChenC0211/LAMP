#pragma once

#include "Core/Lagrangian/InterfaceNode.h"
#include "Core/Field/Field.h"

// Adds normal/tangential stresses to existing Cartesian nodal loads.
// Labels and fields must describe the same active nodes, in the same order.
namespace LegacyNodalMaxwell
{
void add(const std::vector<int>& labels, const std::vector<InterfaceNode>& nodes, const Field<double>& nMaxwellStress, const Field<double>& tMaxwellStress, Field<double>& nodeFx, Field<double>& nodeFy);
}
