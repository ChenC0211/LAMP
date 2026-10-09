#pragma once

struct FittedElectricCaseConfig;

struct PhaseMaterial;

#include "Modules/Electricity/FittedDielectricPolarization.h"
#include "Modules/Electricity/FittedSharpDielectricElectrostatics.h"
#include "Modules/Hydrodynamics/HydroFittedInterfaceLoadGeometry.h"
#include "Modules/Hydrodynamics/HydroFittedCapillary.h"
#include "Modules/Hydrodynamics/HydroFittedScalarGradientLoad.h"
#include "Modules/Hydrodynamics/HydroFittedInterfaceTractionLoad.h"
#include "Modules/Hydrodynamics/HydroFittedForce.h"

struct HydroFittedCapillaryState;
struct ElectrosprayElectricState;
class UnitsContext;

struct ElectrosprayDielectricForces
{
    FittedDielectricPolarization electric;
    HydroFittedScalarGradientLoad polarization;
    HydroFittedForce total;
    explicit ElectrosprayDielectricForces(const QuadtreeMesh& mesh);
};

// App unit binding without IO: Electricity provides the solved coefficient;
// Hydro applies its scalar-gradient load. Only the sum of native force
// densities enters momentum.
ElectrosprayDielectricForces assembleElectrosprayDielectricForces(const QuadtreeMesh& mesh,
    const HydroFittedDomain& domain, const ElectrosprayElectricState& initial,
    const HydroFittedCapillaryState& capillary, double referenceDensitySI, const FittedElectricCaseConfig& electricConfig,
    const UnitsContext& units);

// Full charged Maxwell and capillary loads consume one retained local frame
// and the same fitted pressure gradient. Electric fields/tractions remain in
// electric code units; capillary and total native densities use NS code units.
// The explicit active-side wall-gradient closure is recorded in geometry.
struct ElectrosprayConductingForces
{
    Quantum::Hydrodynamics::HydroFittedInterfaceLoadGeometry geometry;
    HydroFittedCapillaryState capillary;
    FittedSharpDielectricTractions electric;
    Quantum::Hydrodynamics::HydroFittedInterfaceTractionLoad maxwell;
    HydroFittedForce total;
    double electricToHydro = 0.0;
    double geometrySeconds = 0.0, maxwellSeconds = 0.0; // Wall-clock observation only.
    explicit ElectrosprayConductingForces(const QuadtreeMesh& mesh);
};

ElectrosprayConductingForces buildElectrosprayConductingForces(
    const QuadtreeMesh& mesh, const HydroFittedDomain& domain,
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> snapshot,
    const LevelSetInterfacePanels& panels, const SurfaceChargeSource& charge,
    const FittedSharpDielectricResult& electric, const PhaseMaterial& liquid, const PhaseMaterial& gas,
    double surfaceTensionSI, const FittedElectricCaseConfig& electricConfig,
    double capillaryHalfWidthCells, const UnitsContext& units, double sourceTime);
