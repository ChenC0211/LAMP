#include "Electrospray/Setup/ElectrosprayMaterials.h"
#include "Modules/Electricity/FittedElectricCaseConfig.h"
#include "Modules/Hydrodynamics/FittedFlowCaseConfig.h"
#include "Core/Materials/MaterialTable.h"
#include "Electrospray/Setup/ElectrosprayElectrostatics.h"
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/UniformLeafGrid.h"
#include "Core/Units/UnitsContext.h"
#include "Modules/Hydrodynamics/HydroFittedDomain.h"
#include "Modules/Hydrodynamics/Hydrodynamics.h"

#include <cmath>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <optional>

ElectrosprayMaterialState prepareElectrosprayMaterials(
    const QuadtreeMesh& mesh, const HydroFittedDomain& domain,
    const LevelSetHermiteGeometry::Snapshot& currentSnapshot,
    const FittedFlowCaseConfig& flowConfig, const FittedElectricCaseConfig& electricConfig,
    const PhaseMaterial& liquid, const PhaseMaterial& gas, const UnitsContext& units,
    const ElectrosprayElectricState& binding)
{
    const auto require = [](bool condition, const char* reason) {
        if (!condition) throw std::runtime_error(std::string("spray material update: ") + reason);
    };
    require(units.mode() == UnitsMode::kSI, "physical spray materials require explicit SI scales");
    domain.validateMeshIdentity(mesh);
    const bool native = !domain.uniformLeaves();
    std::optional<UniformLeafGrid> grid;
    if (!native) grid.emplace(mesh,CartesianGridScope::UniformLeaves);
    require(domain.axisymmetric() && currentSnapshot.radialAxisSymmetry() && currentSnapshot.restrictedDomain()
        && currentSnapshot.usesNativeLeaves() == native
        && currentSnapshot.minCoord().x == mesh.minCoord.x && currentSnapshot.minCoord().y == mesh.minCoord.y
        && currentSnapshot.maxCoord().x == mesh.maxCoord.x && currentSnapshot.maxCoord().y == mesh.maxCoord.y
        && (native ? currentSnapshot.nativeCellCount() == static_cast<std::size_t>(mesh.cellIndexNum())
            : currentSnapshot.nx() == grid->nx() && currentSnapshot.ny() == grid->ny()
                && currentSnapshot.cellSize().x == grid->spacing().x && currentSnapshot.cellSize().y == grid->spacing().y),
        "current active scalar and fixed RZ leaf layout differ");
    const int count = mesh.cellIndexNum();
    const auto& segments = binding.dielectric.graph().segments();
    const auto& offsets = binding.dielectric.graph().cellOffsets();
    require(binding.dielectric.potential().size() == count
        && binding.phaseByCell.size() == static_cast<std::size_t>(count)
        && binding.dielectric.epsilonByCell().size() == static_cast<std::size_t>(count)
        && offsets.size() == static_cast<std::size_t>(count + 1)
        && offsets.front() == 0 && offsets.back() == segments.size(), "retained material binding storage mismatch");
    require(std::isfinite(liquid.relativePermittivity) && liquid.relativePermittivity > 0.0
        && std::isfinite(gas.relativePermittivity) && gas.relativePermittivity > 0.0,
        "liquid/gas permittivity must be finite and positive");

    // Only the two material arrays are copied. In particular, this helper
    // neither copies a potential solve into a fake new-time electric state nor
    // changes the binding's fixed graph, voltages or face discretization.
    ElectrosprayMaterialState state{binding.phaseByCell, binding.dielectric.epsilonByCell()};
    const double conductorVoltage = units.PotentialSIToCode(electricConfig.conductorPotentialSI);
    for (int position = 0; position < count; ++position)
    {
        const auto& cell = mesh.cells(native ? mesh.cellIndexToLabel[position] : grid->cells()[position]);
        const int k = cell.index();
        const bool active = domain.cellRowByIndex()[k] >= 0;
        const bool scalarActive = native ? currentSnapshot.nativeCellActive(k)
            : currentSnapshot.cellActive(position % grid->nx(),position / grid->nx());
        if (native)
        {
            const auto bounds = currentSnapshot.nativeCellBounds(k);
            const auto lo = mesh.verts(cell.cellToVert[0]).coord(), hi = mesh.verts(cell.cellToVert[2]).coord();
            require(bounds[0].x == lo.x && bounds[0].y == lo.y && bounds[1].x == hi.x && bounds[1].y == hi.y,
                "native scalar leaf bounds differ from actual material control volume");
        }
        require(offsets[k] <= offsets[k + 1] && offsets[k + 1] <= segments.size()
            && (offsets[k] != offsets[k + 1]) == active && scalarActive == active
            && (binding.phaseByCell[k] != electricConfig.conductorPhaseId) == active,
            "fixed conductor, active scalar and retained electric graph disagree");
        if (!active)
        {
            require(binding.dielectric.potential()[k] == conductorVoltage,
                "fixed conductor voltage differs from the retained state");
            continue; // No scalar or material continuation is sampled in solid.
        }
        require(binding.phaseByCell[k] == liquid.phaseId || binding.phaseByCell[k] == gas.phaseId,
            "previous active cell is not a configured fluid phase");
        const double phi = currentSnapshot.value(cell.centerCoord());
        require(std::isfinite(phi), "nonfinite current active scalar");
        const auto& material = phi < 0.0 ? liquid : gas;
        state.phaseByCell[k] = material.phaseId;
        state.epsilonByCell[k] = material.relativePermittivity;
    }

    // Preserve the original separated-interface pressure-port contract.
    // Incoming phase is explicit reservoir data; it is never inferred from
    // a velocity sign or reassigned to the moving cell-center classification.
    using Domain = HydroFittedDomain;
    std::map<int, std::pair<double, int>> pressurePorts;
    require(pressurePorts.emplace(flowConfig.supply.boundaryId,
        std::make_pair(flowConfig.supply.pressureSI, liquid.phaseId)).second, "duplicate pressure port");
    for (const auto& port : flowConfig.ambient)
        require(pressurePorts.emplace(port.boundaryId,
            std::make_pair(port.pressureSI, gas.phaseId)).second, "duplicate pressure port");
    require(pressurePorts.emplace(flowConfig.outlet.boundaryId,
        std::make_pair(flowConfig.outlet.pressureSI, flowConfig.outletIncomingPhaseId)).second, "duplicate pressure port");
    const double velocityScale = units.scales().L0_m / units.scales().T0_s;
    const double pressureScale = liquid.densitySI * velocityScale * velocityScale;
    const double viscosityScale = liquid.densitySI * velocityScale * units.scales().L0_m;
    require(std::isfinite(pressureScale) && pressureScale > 0.0
        && std::isfinite(viscosityScale) && viscosityScale > 0.0, "invalid fixed mechanical scales");
    std::set<int> usedPorts;
    for (const auto& face : domain.faces())
    {
        if (face.boundary.type != Domain::BoundaryType::PressureReservoir
            && face.boundary.type != Domain::BoundaryType::PressureOutlet) continue;
        const auto found = pressurePorts.find(face.boundary.id);
        require(found != pressurePorts.end(), "fixed domain has an unconfigured pressure port");
        const int k = face.negativeCell >= 0 ? face.negativeCell : face.positiveCell;
        require(k >= 0 && k < count && (face.negativeCell < 0) != (face.positiveCell < 0),
            "pressure face must have exactly one active owner");
        const int incoming = found->second.second;
        const int expectedPhase = incoming == -1 ? binding.phaseByCell[k] : incoming;
        require(expectedPhase == liquid.phaseId || expectedPhase == gas.phaseId,
            "pressure port has an unsupported incoming phase");
        require(face.boundary.pressure == found->second.first / pressureScale,
            "fixed pressure boundary value changed");
        if (incoming == -1)
            require(face.boundary.type == Domain::BoundaryType::PressureOutlet,
                "a reservoir cannot omit incoming material");
        else
        {
            const auto& material = incoming == liquid.phaseId ? liquid : gas;
            require(face.boundary.type == Domain::BoundaryType::PressureReservoir
                && face.boundary.incomingDensity == material.densitySI / liquid.densitySI
                && face.boundary.incomingViscosity == material.dynamicViscositySI / viscosityScale,
                "fixed reservoir material differs from its explicit incoming phase");
        }
        require(binding.phaseByCell[k] == expectedPhase && state.phaseByCell[k] == expectedPhase
            && currentSnapshot.closedCellSign(mesh.cells(mesh.cellIndexToLabel[k]).centerCoord())
                == (expectedPhase == liquid.phaseId ? -1 : 1),
            "free interface reached or changed a pressure-port owner cell");
        usedPorts.insert(face.boundary.id);
    }
    require(usedPorts.size() == pressurePorts.size(), "configured pressure port is absent from the fixed domain");
    return state;
}

ElectrosprayHydroScales::ElectrosprayHydroScales(double referenceDensitySI, const UnitsContext& units)
{
    if (units.mode() != UnitsMode::kSI) throw std::invalid_argument("physical spray Hydro requires SI scales");
    length = units.scales().L0_m;
    time = units.scales().T0_s;
    density = referenceDensitySI;
    velocity = length / time;
    pressure = density * velocity * velocity;
    viscosity = density * velocity * length;
    for (double value : {length, time, density, velocity, pressure, viscosity})
        if (!std::isfinite(value) || !(value > 0.0)) throw std::invalid_argument("invalid spray Hydro reference scale");
}

void updateElectrosprayHydroMaterials(QuadtreeMesh& mesh, const HydroFittedDomain& domain, Hydrodynamics& hydro,
    const std::vector<int>& phaseByCell, const PhaseMaterial& liquid, const PhaseMaterial& gas, const ElectrosprayHydroScales& scale,
    double sourceTime, HydroMaterialEvolution evolution)
{
    hydro.updateMaterials(mesh, {HydroStepGeometry::fittedRZ(domain), sourceTime, evolution, HydroClassifiedMaterial{&phaseByCell,
        {liquid.phaseId, liquid.densitySI / scale.density, liquid.dynamicViscositySI / scale.viscosity},
        {gas.phaseId, gas.densitySI / scale.density, gas.dynamicViscositySI / scale.viscosity}}});
}

void updateElectrosprayHydroMaterials(const QuadtreeMesh& mesh, const HydroFittedDomain& domain, HydroMaterialFields fields,
    const std::vector<int>& phaseByCell, const PhaseMaterial& liquid, const PhaseMaterial& gas, const ElectrosprayHydroScales& scale)
{
    assignHydroFittedMaterials(mesh, domain, fields, phaseByCell,
        {liquid.phaseId, liquid.densitySI / scale.density, liquid.dynamicViscositySI / scale.viscosity},
        {gas.phaseId, gas.densitySI / scale.density, gas.dynamicViscositySI / scale.viscosity});
}
