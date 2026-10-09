#pragma once

#include "Core/Lagrangian/InterfaceNode.h"
#include <iosfwd>

// SIC2's native node record predates the geometry/physical-state split.
// Keep that wire record in Hydro IO; Core nodes contain no evaporation state.
namespace LegacyNodalCache
{
void writeNodes(std::ostream& stream, const std::vector<InterfaceNode>& nodes,
    const std::vector<double>& evaporationVelocity, int count);
void readNodes(std::istream& stream, std::vector<InterfaceNode>& nodes,
    std::vector<double>& evaporationVelocity, int count);
}
