#include "Electrospray/Setup/ElectrosprayElectrostatics.h"
#include "Modules/Electricity/FittedElectricCaseConfig.h"
#include "Modules/Hydrodynamics/FittedFlowCaseConfig.h"
#include "Modules/Interface/Compatibility/InterfaceCaseConfig.h"
#include "Core/Materials/MaterialTable.h"
#include "Electrospray/Setup/ElectrosprayMaterials.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>
#include "Core/Geometry/ImmersedGeometry.h"
#include "Core/Geometry/FittedDomainGeometry.h"
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/IO/DictionaryGeometry.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Units/UnitsContext.h"
#include "Modules/Hydrodynamics/HydroFittedDomain.h"

namespace
{
using Solver = PolarizationPoissonFVSolverCellCentered;
using BCType = Solver::FittedBoundaryType;
constexpr double vacuumPermittivitySI = 8.8541878128e-12;
Quantum::Electricity::FittedCellHarmonicControls cellHarmonicControls(
    const FittedElectricCaseConfig& config,const UnitsContext& units,double linearTolerance)
{
    double lower=config.conductorPotentialSI,upper=lower;
    for(const auto& bc:config.domainBoundaries) if(bc.type=="dirichlet") {
        lower=std::min(lower,bc.potentialSI);upper=std::max(upper,bc.potentialSI);
    }
    return {linearTolerance,config.potentialBounds,units.scales().Phi0_V,lower,upper};
}
}

ElectrosprayElectricState bindElectrosprayElectrostatics(
    const QuadtreeMesh& mesh, const DictionaryGeometry& definition,
    const FittedElectricCaseConfig& electricConfig,
    const PhaseMaterial& liquid, const PhaseMaterial& gas, const UnitsContext& units,
    const std::vector<int>& phaseByCell, const CellField<double>* potentialGuess,
    const Quantum::Electricity::FittedCellHarmonicControls* solveControls)
{
    ElectrosprayElectricState state;
    state.phaseByCell = phaseByCell;
    if (phaseByCell.size() != static_cast<std::size_t>(mesh.cellIndexNum()))
        throw std::invalid_argument("electric binding phase count differs from mesh");
    std::vector<double> epsilonByCell(mesh.cellIndexNum(),0.0);
    FaceSegmentGraph graph;
    std::vector<Solver::FittedBoundary> boundaries;
    std::map<int, int> phaseCounts;
    for (int i = 0; i < mesh.cellIndexNum(); ++i)
    {
        const int phase = state.phaseByCell[i];
        ++phaseCounts[phase];
        if (phase == electricConfig.conductorPhaseId) continue;
        else if (phase == liquid.phaseId) epsilonByCell[i] = liquid.relativePermittivity;
        else if (phase == gas.phaseId) epsilonByCell[i] = gas.relativePermittivity;
        else throw std::runtime_error("cell has no configured physical phase");
    }
    if (phaseCounts.size() != 3) throw std::runtime_error("initial spray requires represented conductor, liquid and gas");
    graph.build(mesh, [&](const QuadtreeMesh&, const Cell& cell)
        { return state.phaseByCell[cell.index()] != electricConfig.conductorPhaseId; });
    const auto& segments = graph.segments();
    boundaries.resize(segments.size());
    state.boundaryIdBySegment.assign(segments.size(), -1);
    std::map<int, FittedElectricBoundarySI> domainBCs;
    for (const auto& bc : electricConfig.domainBoundaries) domainBCs.emplace(bc.boundaryId, bc);
    std::set<int> usedDomainBCs;
    std::map<int, double> coveredLength;
    const double tolerance = 64.0 * std::numeric_limits<double>::epsilon()
        * std::max({1.0, std::abs(definition.domainMin.x), std::abs(definition.domainMin.y),
            std::abs(definition.domainMax.x), std::abs(definition.domainMax.y)});
    for (std::size_t j = 0; j < segments.size(); ++j)
    {
        const auto& segment = segments[j];
        const int k = mesh.cells(segment.K).index();
        const int l = segment.L < 0 ? -1 : mesh.cells(segment.L).index();
        if (l >= 0 && state.phaseByCell[l] != electricConfig.conductorPhaseId) continue;
        const auto& boundary = FittedDomainGeometry::bindFace(definition, segment, state.phaseByCell[k], electricConfig.conductorPhaseId, tolerance);
        state.boundaryIdBySegment[j] = boundary.id;
        coveredLength[boundary.id] += segment.Af;
        auto& fitted = boundaries[j];
        if (l >= 0)
        {
            fitted.type = BCType::Dirichlet;
            fitted.value = units.PotentialSIToCode(electricConfig.conductorPotentialSI);
        }
        else
        {
            const auto bc = domainBCs.find(boundary.id);
            if (bc == domainBCs.end()) throw std::runtime_error("missing electric data for domain boundary " + std::to_string(boundary.id));
            usedDomainBCs.insert(boundary.id);
            fitted.type = bc->second.type == "dirichlet" ? BCType::Dirichlet : BCType::HomogeneousNeumann;
            fitted.value = units.PotentialSIToCode(bc->second.potentialSI);
        }
    }
    if (usedDomainBCs.size() != domainBCs.size()) throw std::runtime_error("unused or inactive electric domain boundary configuration");
    // Certify whole declared solid boundary coverage, not only selected face
    // midpoints. No curved/slanted solid is silently converted to a stair-step.
    for (const auto& boundary : definition.boundaries)
    {
        if (boundary.ignored || boundary.kind == GeometricBoundary::domain
            || (boundary.phase + 1 != electricConfig.conductorPhaseId && boundary.opposedPhase + 1 != electricConfig.conductorPhaseId)) continue;
        if (boundary.kind != GeometricBoundary::fixed) throw std::runtime_error("moving conductor unsupported");
        double declaredLength = 0.0;
        for (int lineId : boundary.lineIds)
        {
            const auto& line = definition.lines.at(lineId);
            if (line.type() != LineType::straight) throw std::runtime_error("solid must have face-fitted straight boundaries");
            declaredLength += std::hypot(line.points[1].x - line.points[0].x, line.points[1].y - line.points[0].y);
        }
        if (!(declaredLength > 0.0) || std::abs(coveredLength[boundary.id] - declaredLength) > tolerance * segments.size())
            throw std::runtime_error("incomplete fitted conductor boundary " + std::to_string(boundary.id));
    }
    std::optional<FittedPotentialInput> seed;
    if(potentialGuess) seed.emplace(*potentialGuess);
    state.dielectric=Quantum::Electricity::prepareFittedCellHarmonic(mesh,std::move(graph),std::move(boundaries),
        std::move(epsilonByCell),true,units.PotentialSIToCode(electricConfig.conductorPotentialSI),
        seed?&*seed:nullptr,solveControls);
    return state;
}

ElectrosprayElectricState initializeElectrosprayElectrostatics(
    const QuadtreeMesh& mesh, const ImmersedGeometry<QuadtreeMesh>& geometry,
    const DictionaryGeometry& definition, const FittedFlowCaseConfig& flowConfig,
    const FittedInterfaceCaseConfig& interfaceConfig, const FittedElectricCaseConfig& electricConfig,
    const PhaseMaterial& liquid, const PhaseMaterial& gas, double electricLinearTolerance,
    const UnitsContext& units)
{
    if (units.mode() != UnitsMode::kSI || (electricConfig.model != "leakyDielectric" && electricConfig.model != "perfectDielectric")
        || (electricConfig.interfaceDiscretization != "cellHarmonic" && electricConfig.interfaceDiscretization != "sharpHermiteP3")
        || electricConfig.initialSurfaceChargeSI != 0.0)
        throw std::runtime_error("initial spray bootstrap requires an SI dielectric model and zero free surface charge");
    // In explicit P3 mode this solve only initializes the fixed electrode
    // binding and a starting potential. It is replaced after the shared free
    // surface and charge panels exist; it is not a P3/conducting time result.
    if (!definition.periodicPairs.empty()) throw std::runtime_error("periodic spray electric boundaries are not supported");
    const auto moving = definition.boundaryById.find(interfaceConfig.movingBoundaryId);
    const auto ls = std::find(definition.levelSetName.begin(), definition.levelSetName.end(), interfaceConfig.fluidInterfaceLevelSet);
    if (ls == definition.levelSetName.end() || moving == definition.boundaryById.end()
        || moving->second.kind != GeometricBoundary::moving || moving->second.ignored
        || moving->second.levelSetIndex != std::distance(definition.levelSetName.begin(), ls)
        || !((moving->second.phase + 1 == liquid.phaseId && moving->second.opposedPhase + 1 == gas.phaseId)
            || (moving->second.opposedPhase + 1 == liquid.phaseId && moving->second.phase + 1 == gas.phaseId)))
        throw std::runtime_error("moving electric interface must select the liquid/gas boundary, not a solid wall");
    for (const auto& pressure : {flowConfig.supply, flowConfig.outlet})
    {
        const auto boundary = definition.boundaryById.find(pressure.boundaryId);
        if (boundary == definition.boundaryById.end() || boundary->second.kind != GeometricBoundary::domain
            || boundary->second.ignored || (boundary->second.adjacentPhase + 1 != liquid.phaseId
                && boundary->second.adjacentPhase + 1 != gas.phaseId))
            throw std::runtime_error("pressure data must refer to an active fluid domain boundary");
    }

    const auto controls=cellHarmonicControls(electricConfig,units,electricLinearTolerance);
    return bindElectrosprayElectrostatics(mesh,definition,electricConfig,liquid,gas,units,
        FittedDomainGeometry::cellPhases(mesh,geometry),nullptr,&controls);
}

Quantum::Electricity::FittedChargeParameters bindElectrosprayConductingParameters(
    const FittedElectricCaseConfig& electricConfig, const PhaseMaterial& liquid,
    const PhaseMaterial& gas, double electricLinearTolerance, const UnitsContext& units)
{
    if(units.mode()!=UnitsMode::kSI || electricConfig.model!="leakyDielectric"
        || electricConfig.interfaceDiscretization!="sharpHermiteP3")
        throw std::runtime_error("conducting spray requires explicit SI leakyDielectric/sharpHermiteP3 inputs");
    Quantum::Electricity::FittedChargeParameters result;
    auto& electric = result.field;
    electric.model=electricConfig.model; electric.interfaceDiscretization=electricConfig.interfaceDiscretization;
    electric.permittivityNegative=liquid.relativePermittivity;
    electric.permittivityPositive=gas.relativePermittivity;
    electric.axisymmetric=true; electric.linearResidualTolerance=electricLinearTolerance;
    electric.sharpFluxChecks=electricConfig.fluxChecks;
    // kappa=epsilon_r, x*=x/L0, phi*=phi/Phi0 and t*=t/T0 imply
    // sigma*=sigma_SI*T0/epsilon0. Field charge/source already uses the same
    // q0=epsilon0*Phi0/L0 and RZ Q0=epsilon0*Phi0*L0 normalization.
    const double timeScale=units.scales().T0_s;
    const DielectricConductivity conductivity{
        liquid.conductivitySI*timeScale/vacuumPermittivitySI,
        gas.conductivitySI*timeScale/vacuumPermittivitySI};
    if (!std::isfinite(timeScale) || timeScale <= 0.0
        || !std::isfinite(conductivity.negative) || !std::isfinite(conductivity.positive)
        || (liquid.conductivitySI > 0.0 && conductivity.negative <= 0.0)
        || (gas.conductivitySI > 0.0 && conductivity.positive <= 0.0))
        throw std::runtime_error("conductivity SI-to-code conversion overflows or loses a nonzero conductivity");
    result.conductivity = conductivity;
    result.closureTolerance = electricConfig.chargeClosure;
    return result;
}

Quantum::Electricity::FittedChargeState initializeElectrosprayConductingState(
    const QuadtreeMesh& mesh, std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> snapshot,
    const LevelSetInterfacePanels& panels, const std::vector<Vector2D>& fixedContacts,
    const FittedElectricCaseConfig& electricConfig,
    const PhaseMaterial& liquid, const PhaseMaterial& gas, double electricLinearTolerance, const UnitsContext& units,
    const ElectrosprayElectricState& binding)
{
    const Quantum::Electricity::FittedChargeProblem problem{binding.dielectric.graph(),binding.dielectric.boundaries(),
        bindElectrosprayConductingParameters(electricConfig,liquid,gas,electricLinearTolerance,units)};
    return Quantum::Electricity::initializeFittedChargeState(mesh,problem,
        std::move(snapshot),panels,fixedContacts,binding.dielectric.potential());
}

ElectrosprayElectricState updateElectrosprayDielectricState(
    const QuadtreeMesh& mesh, const HydroFittedDomain& domain,
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> currentSnapshot,
    const FittedFlowCaseConfig& flowConfig, const FittedElectricCaseConfig& electricConfig,
    const PhaseMaterial& liquid, const PhaseMaterial& gas, double electricLinearTolerance, const UnitsContext& units,
    const ElectrosprayElectricState& previousState, double time)
{
    const auto require = [](bool condition, const char* reason) {
        if (!condition) throw std::runtime_error(std::string("electric material update: ") + reason);
    };
    require(currentSnapshot && units.mode() == UnitsMode::kSI && (electricConfig.model == "leakyDielectric" || electricConfig.model == "perfectDielectric")
        && electricConfig.interfaceDiscretization == "cellHarmonic" && electricConfig.initialSurfaceChargeSI == 0.0,
        "requires SI leakyDielectric/cellHarmonic with zero free surface charge; no charge fallback exists");
    const int count = mesh.cellIndexNum();
    const auto& segments = previousState.dielectric.graph().segments();
    require(previousState.dielectric.epsilonBySegment().size() == segments.size()
        && previousState.dielectric.boundaries().size() == segments.size()
        && previousState.boundaryIdBySegment.size() == segments.size(), "previous electric state storage mismatch");
    auto materials = prepareElectrosprayMaterials(mesh, domain, *currentSnapshot, flowConfig, electricConfig, liquid, gas, units, previousState);

    ElectrosprayElectricState state;
    state.boundaryIdBySegment=previousState.boundaryIdBySegment;
    state.phaseByCell = std::move(materials.phaseByCell);
    const double conductorVoltage = units.PotentialSIToCode(electricConfig.conductorPotentialSI);

    std::map<int, FittedElectricBoundarySI> electricPorts;
    for (const auto& bc : electricConfig.domainBoundaries)
        require(electricPorts.emplace(bc.boundaryId, bc).second, "duplicate electric domain boundary");
    std::set<int> usedElectricPorts;
    for (std::size_t j = 0; j < segments.size(); ++j)
    {
        const auto& segment = segments[j];
        const int k = mesh.cells(segment.K).index();
        const int l = segment.L < 0 ? -1 : mesh.cells(segment.L).index();
        const auto& boundary = previousState.dielectric.boundaries()[j];
        require(domain.cellRowByIndex()[k] >= 0, "retained electric segment owner is not active");
        if (l >= 0 && domain.cellRowByIndex()[l] >= 0)
        {
            require(boundary.type == BCType::Interior && state.boundaryIdBySegment[j] == -1,
                "retained active-active electric boundary changed");

        }
        else
        {
            if (l >= 0)
                require(state.phaseByCell[l] == electricConfig.conductorPhaseId
                    && boundary.type == BCType::Dirichlet && boundary.value == conductorVoltage,
                    "fixed conductor face or voltage changed");
            else
            {
                const auto found = electricPorts.find(state.boundaryIdBySegment[j]);
                require(found != electricPorts.end(), "retained electric domain boundary is unconfigured");
                const auto& bc = found->second;
                require((bc.type == "dirichlet" || bc.type == "homogeneousNeumann")
                    && boundary.type == (bc.type == "dirichlet" ? BCType::Dirichlet : BCType::HomogeneousNeumann)
                    && boundary.value == units.PotentialSIToCode(bc.potentialSI), "fixed electric domain data changed");
                usedElectricPorts.insert(bc.boundaryId);
            }
        }
    }
    require(usedElectricPorts.size() == electricPorts.size(), "unused electric domain boundary data");
    state.dielectric=Quantum::Electricity::updateFittedCellHarmonic(mesh,previousState.dielectric,
        std::move(materials.epsilonByCell),cellHarmonicControls(electricConfig,units,electricLinearTolerance),
        time,std::move(currentSnapshot));
    return state;
}
