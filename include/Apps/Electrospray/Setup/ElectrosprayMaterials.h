#pragma once

struct FittedElectricCaseConfig;

struct FittedFlowCaseConfig;

struct PhaseMaterial;

#include <vector>
#include "Modules/Hydrodynamics/HydroMaterialFields.h"

class QuadtreeMesh;
class HydroFittedDomain;
class UnitsContext;
class Hydrodynamics;
template<typename T> class CellField;
template<typename T> class VertField;
namespace LevelSetHermiteGeometry { class Snapshot; }
struct ElectrosprayElectricState;

// App-level material classification on the fixed active-fluid mask. Native
// cell indices and dictionary phase IDs are retained; conductor entries are
// copied unchanged. Cell epsilon is a phase property, not a face coefficient
// or an alternative reconstruction of the sharp electric operator.
struct ElectrosprayMaterialState
{
    std::vector<int> phaseByCell;
    std::vector<double> epsilonByCell;
};

// Select negative phi as liquid and nonnegative phi as gas on active cells.
// Validate the original fixed conductor voltage/mask, active scalar/graph
// identity, and pressure-port phase, pressure and incoming-material data.
// Explicit units are required for those SI-to-code boundary checks.
// No field is changed, no face harmonic coefficient is formed, and no
// Poisson/charge/flow solve occurs. No q=0 or electric-discretization choice
// is made here. Failure leaves the binding and scalar unchanged.
ElectrosprayMaterialState prepareElectrosprayMaterials(
    const QuadtreeMesh& mesh, const HydroFittedDomain& domain,
    const LevelSetHermiteGeometry::Snapshot& currentSnapshot,
    const FittedFlowCaseConfig& flowConfig, const FittedElectricCaseConfig& electricConfig,
    const PhaseMaterial& liquid, const PhaseMaterial& gas, const UnitsContext& units,
    const ElectrosprayElectricState& binding);

// Existing SI mechanical normalization, shared by boundary/material binding,
// the time driver and its diagnostics. Charge normalization remains separate.
struct ElectrosprayHydroScales
{
    double density, velocity, pressure, viscosity, length, time;
    ElectrosprayHydroScales(double referenceDensitySI, const UnitsContext& units);
};

// Populate cell and active-only vertex rho/mu on the same native fluid mask.
// Does not reclassify phases, modify the scalar, solve NS or write diagnostics.
void updateElectrosprayHydroMaterials(QuadtreeMesh& mesh, const HydroFittedDomain& domain,
    Hydrodynamics& hydro, const std::vector<int>& phaseByCell,
    const PhaseMaterial& liquid, const PhaseMaterial& gas, const ElectrosprayHydroScales& scale,
    double sourceTime, HydroMaterialEvolution evolution);

// Same closure into staged module-owned fields, without constructing another
// Hydrodynamics solver just to carry material arrays during AMR.
void updateElectrosprayHydroMaterials(const QuadtreeMesh& mesh, const HydroFittedDomain& domain,
    HydroMaterialFields fields, const std::vector<int>& phaseByCell,
    const PhaseMaterial& liquid, const PhaseMaterial& gas, const ElectrosprayHydroScales& scale);
