#pragma once

#include "Modules/Hydrodynamics/HydroFittedDomain.h"

class QuadtreeMesh;
class DictionaryGeometry;
class UnitsContext;
struct ElectrosprayElectricState;
struct FittedFlowCaseConfig;
struct PhaseMaterial;

// Bind dictionary IDs to the SAME actual active-fluid graph already certified
// against the solid geometry. No nozzle coordinates or material constants are
// embedded in the Hydro module.
HydroFittedDomain buildElectrosprayHydroDomain(const QuadtreeMesh& mesh,
    const ElectrosprayElectricState& initial, const DictionaryGeometry& definition,
    const FittedFlowCaseConfig& flow, const PhaseMaterial& liquid, const PhaseMaterial& gas,
    int conductorPhaseId, const UnitsContext& units);
