#pragma once

#include <string>
class QuadtreeMesh;
class UnitsContext;
struct ElectrosprayElectricState;

// Called after the initialization solve/checks, even when spatial output is off.
void reportElectrosprayElectricInitialization(const QuadtreeMesh& mesh,
    const ElectrosprayElectricState& state, double conductorPotentialSI,
    const UnitsContext& units);

// Authoritative cell potentials and oriented FV face fluxes in SI. No interpolated
// nodal field is presented as an additional numerical solution.
void writeElectrosprayElectrostatics(const QuadtreeMesh& mesh,
    const ElectrosprayElectricState& state, const UnitsContext& units,
    const std::string& outputPrefix, bool initialReport = true, bool verboseDiagnostics = false);
