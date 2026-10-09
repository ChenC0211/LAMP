#pragma once

struct FittedInterfaceCaseConfig;

#include <string>
#include "Electrospray/Setup/ElectrosprayInterface.h"
class UnitsContext;

void reportElectrosprayInterfaceInitialization(const QuadtreeMesh& mesh,
    const ElectrosprayInterfaceState& state, const FittedInterfaceCaseConfig& interfaceConfig);
void writeElectrosprayInitialInterface(const QuadtreeMesh& mesh,
    const ElectrosprayInterfaceState& state, const UnitsContext& units, const std::string& prefix);

// Captures mesh/state/units by reference and the prefix by value. Use only
// synchronously while those objects live, during initial geometry construction.
// An empty prefix suppresses files, not the construction progress summaries.
InterfaceInitialGeometryObserver makeElectrosprayInitialGeometryObserver(const QuadtreeMesh& mesh,
    const ElectrosprayInterfaceState& state, const UnitsContext& units, const std::string& outputPrefix);
