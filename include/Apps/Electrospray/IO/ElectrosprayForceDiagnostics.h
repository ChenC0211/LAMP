#pragma once

struct RunControls;

class QuadtreeMesh;
class HydroFittedDomain;
class UnitsContext;
struct HydroFittedCapillaryState;
struct ElectrosprayDielectricForces;
struct ElectrosprayConductingForces;

// Diagnostics consume completed computations. No force assembly or state
// mutation occurs here; the caller determines the candidate/accepted stage.
void writeElectrosprayCapillaryDiagnostics(const QuadtreeMesh& mesh,
    const HydroFittedCapillaryState& state, double referenceDensitySI, double capillaryHalfWidthCells,
    const RunControls& run, const UnitsContext& units);
void writeElectrosprayDielectricForceDiagnostics(const QuadtreeMesh& mesh,
    const HydroFittedDomain& domain, const ElectrosprayDielectricForces& forces,
    const HydroFittedCapillaryState& capillary, double referenceDensitySI,
    const RunControls& run, const UnitsContext& units);
void reportElectrosprayConductingForces(const ElectrosprayConductingForces& forces,
    const UnitsContext& units);
