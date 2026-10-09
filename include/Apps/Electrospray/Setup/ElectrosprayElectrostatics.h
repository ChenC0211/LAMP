#pragma once

struct FittedElectricCaseConfig;

struct FittedFlowCaseConfig;
struct FittedInterfaceCaseConfig;

struct PhaseMaterial;

#include <string>
#include <vector>
#include "Core/Field/CellField.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Solvers/Poisson/PolarizationPoissonFVSolverCellCentered.h"
#include "Modules/Electricity/FittedChargeEvolution.h"
#include "Modules/Electricity/FittedCellHarmonic.h"

class QuadtreeMesh;
class DictionaryGeometry;
class UnitsContext;
class HydroFittedDomain;
namespace LevelSetHermiteGeometry { class Snapshot; }
template<class Mesh> class ImmersedGeometry;

// App binding, not another Poisson solver or mesh. The initial liquid/gas
// interface uses cell-material harmonic epsilon; solid boundaries must actually
// coincide with mesh faces. Conductivity and interface motion are NOT advanced.
struct ElectrosprayElectricState
{
    Quantum::Electricity::FittedCellHarmonicState dielectric;
    std::vector<int> phaseByCell; // App's dictionary material IDs, not numerical unknowns.
    std::vector<int> boundaryIdBySegment; // -1 for active-active segments.
};

// Rebuild the ORIGINAL fitted graph/electrode binding on explicit cell phases.
// By default this prepares an unsolved seed for AMR/charged initialization.
// Bootstrap controls additionally ask Electricity to solve before publication.
// Optional guess uses target mesh labels; conductor voltage comes from the
// unchanged physical input. No charge update or free-interface reconstruction.
ElectrosprayElectricState bindElectrosprayElectrostatics(
    const QuadtreeMesh& mesh, const DictionaryGeometry& definition,
    const FittedElectricCaseConfig& electricConfig,
    const PhaseMaterial& liquid, const PhaseMaterial& gas, const UnitsContext& units,
    const std::vector<int>& phaseByCell, const CellField<double>* potentialGuess = nullptr,
    const Quantum::Electricity::FittedCellHarmonicControls* solveControls = nullptr);

ElectrosprayElectricState initializeElectrosprayElectrostatics(
    const QuadtreeMesh& mesh, const ImmersedGeometry<QuadtreeMesh>& geometry,
    const DictionaryGeometry& definition, const FittedFlowCaseConfig& flowConfig,
    const FittedInterfaceCaseConfig& interfaceConfig, const FittedElectricCaseConfig& electricConfig,
    const PhaseMaterial& liquid, const PhaseMaterial& gas, double electricLinearTolerance,
    const UnitsContext& units);

// Re-solve the same fitted q=0 electric problem with cell materials selected by
// the current active free-interface scalar (negative = liquid). The fixed solid
// mask, electrical graph and voltage data are retained; pressure-port owner
// cells must remain in their configured phase. No conductivity/charge equation
// is advanced. Failure leaves previousState and the retained scalar unchanged.
ElectrosprayElectricState updateElectrosprayDielectricState(
    const QuadtreeMesh& mesh, const HydroFittedDomain& domain,
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> currentSnapshot,
    const FittedFlowCaseConfig& flowConfig, const FittedElectricCaseConfig& electricConfig,
    const PhaseMaterial& liquid, const PhaseMaterial& gas, double electricLinearTolerance, const UnitsContext& units,
    const ElectrosprayElectricState& previousState, double time);

// Bind SI inputs and initialize one complete module-owned electrical state.
Quantum::Electricity::FittedChargeState initializeElectrosprayConductingState(
    const QuadtreeMesh& mesh, std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> snapshot,
    const LevelSetInterfacePanels& panels, const std::vector<Vector2D>& fixedContacts,
    const FittedElectricCaseConfig& electricConfig,
    const PhaseMaterial& liquid, const PhaseMaterial& gas, double electricLinearTolerance, const UnitsContext& units,
    const ElectrosprayElectricState& binding);

// Convert the selected material/SI inputs without selecting a numerical backend.
Quantum::Electricity::FittedChargeParameters bindElectrosprayConductingParameters(
    const FittedElectricCaseConfig&, const PhaseMaterial& liquid, const PhaseMaterial& gas,
    double electricLinearTolerance, const UnitsContext&);
