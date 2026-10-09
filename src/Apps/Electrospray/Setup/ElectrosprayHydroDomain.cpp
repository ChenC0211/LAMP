#include "Electrospray/Setup/ElectrosprayHydroDomain.h"
#include "Electrospray/Setup/ElectrosprayElectrostatics.h"
#include "Electrospray/Setup/ElectrosprayMaterials.h"
#include "Modules/Hydrodynamics/HydroFittedMomentum.h"
#include "Core/Materials/MaterialTable.h"
#include "Modules/Hydrodynamics/FittedFlowCaseConfig.h"
#include "Core/IO/DictionaryGeometry.h"
#include "Core/Mesh/QuadtreeMesh.h"

#include <map>
#include <stdexcept>

namespace { using Domain = HydroFittedDomain; }

HydroFittedDomain buildElectrosprayHydroDomain(const QuadtreeMesh& mesh,
    const ElectrosprayElectricState& initial, const DictionaryGeometry& definition,
    const FittedFlowCaseConfig& flow, const PhaseMaterial& liquid, const PhaseMaterial& gas,
    int conductorPhaseId, const UnitsContext& units)
{
    if (flow.solidWallModel != "noSlip" || flow.incomingMomentumModel != "interiorContinuationZeroShear")
        throw std::invalid_argument("unsupported spray mechanical wall/incoming momentum model");
    if (liquid.phaseId != flow.liquidPhase || gas.phaseId != flow.gasPhase)
        throw std::invalid_argument("Hydro material identities differ from configured flow phases");
    const ElectrosprayHydroScales scale(liquid.densitySI, units);
    std::map<int, Domain::Boundary> pressureBoundaries;
    const auto reservoir = [&](int id, double pressure, const PhaseMaterial& material)
    {
        const auto found = definition.boundaryById.find(id);
        if (found == definition.boundaryById.end() || found->second.ignored
            || found->second.kind != GeometricBoundary::domain || found->second.adjacentPhase + 1 != material.phaseId)
            throw std::invalid_argument("pressure reservoir has missing/wrong initial incoming material boundary");
        Domain::Boundary boundary;
        boundary.type = Domain::BoundaryType::PressureReservoir;
        boundary.id = id;
        boundary.pressure = pressure / scale.pressure;
        boundary.incomingDensity = material.densitySI / scale.density;
        boundary.incomingViscosity = material.dynamicViscositySI / scale.viscosity;
        if (!pressureBoundaries.emplace(id, boundary).second) throw std::invalid_argument("duplicate Hydro pressure boundary");
    };
    reservoir(flow.supply.boundaryId, flow.supply.pressureSI, liquid);
    for (const auto& boundary : flow.ambient) reservoir(boundary.boundaryId, boundary.pressureSI, gas);
    if (flow.outletIncomingPhaseId != -1)
    {
        if (flow.outletIncomingPhaseId != liquid.phaseId
            && flow.outletIncomingPhaseId != gas.phaseId)
            throw std::invalid_argument("outlet incoming phase must be configured liquid or gas");
        const auto& material = flow.outletIncomingPhaseId == liquid.phaseId ? liquid : gas;
        reservoir(flow.outlet.boundaryId, flow.outlet.pressureSI, material);
    }
    else
    {
        Domain::Boundary outlet;
        outlet.type = Domain::BoundaryType::PressureOutlet;
        outlet.id = flow.outlet.boundaryId;
        outlet.pressure = flow.outlet.pressureSI / scale.pressure;
        if (!pressureBoundaries.emplace(outlet.id, outlet).second) throw std::invalid_argument("duplicate Hydro outlet boundary");
    }
    std::vector<Domain::Boundary> boundaries(initial.dielectric.graph().segments().size());
    std::map<int, int> covered;
    for (std::size_t i = 0; i < boundaries.size(); ++i)
    {
        const auto& segment = initial.dielectric.graph().segments()[i];
        const int boundaryId = initial.boundaryIdBySegment[i];
        if (boundaryId < 0) continue;
        auto& boundary = boundaries[i];
        boundary.id = boundaryId;
        if (segment.L >= 0)
        {
            if (initial.phaseByCell[mesh.cells(segment.L).index()] != conductorPhaseId)
                throw std::invalid_argument("internal Hydro boundary neighbor is not configured conductor");
            boundary.type = Domain::BoundaryType::NoSlipWall;
        }
        else if (const auto pressure = pressureBoundaries.find(boundaryId); pressure != pressureBoundaries.end())
        {
            boundary = pressure->second;
            ++covered[boundaryId];
        }
        else if (segment.faceCenter.x == 0.0 && segment.n.x == -1.0 && mesh.minCoord.x == 0.0)
            boundary.type = Domain::BoundaryType::Symmetry;
        else throw std::invalid_argument("unconfigured mechanical fluid domain boundary " + std::to_string(boundaryId));
    }
    if (covered.size() != pressureBoundaries.size()) throw std::invalid_argument("unused Hydro pressure boundary data");
    auto result = Domain::build(mesh, initial.dielectric.graph(), boundaries, true);
    // Reject an unsupported boundary halo before costly interface/electric loads.
    if (!result.uniformLeaves()) HydroFittedMomentumLayout::build(mesh, result);
    return result;
}
