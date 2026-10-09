#include "Modules/Hydrodynamics/HydroInterfaceBoundary.h"

#include "Modules/Interface/Compatibility/InterfaceBoundaryGeometry.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Field/HEdgeField.h"
#include "Core/Field/VEdgeField.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace HydroInterfaceBoundarySupport
{
namespace
{
Boundary sideOf(HydroBoundarySide side)
{
    switch (side)
    {
    case HydroBoundarySide::West: return Boundary::west;
    case HydroBoundarySide::South: return Boundary::south;
    case HydroBoundarySide::East: return Boundary::east;
    case HydroBoundarySide::North: return Boundary::north;
    }
    throw std::invalid_argument("invalid Hydro interface boundary side");
}

}

LevelSetBoundaryConditions scalarConditions(const InterfaceNodalConfig::ScalarBoundary& boundary,
    const std::vector<HydroBoundaryPatch>& patches)
{
    LevelSetBoundaryConditions result;
    result.fixedZeroVertices = boundary.fixedContacts;
    for (const auto& profile : boundary.inflow)
    {
        const auto patch = std::find_if(patches.begin(), patches.end(),
            [&](const auto& value) { return value.name == profile.patch; });
        if (patch == patches.end() || patch->type != HydroBoundaryType::NormalVelocity
            || !(patch->value < 0.0))
            throw std::invalid_argument("interface inflow must name a prescribed inward velocity patch");
        const bool horizontal = patch->side == HydroBoundarySide::South || patch->side == HydroBoundarySide::North;
        result.inflow.push_back({sideOf(patch->side), patch->start, patch->end,
            [profile, horizontal](Vector2D point, double) {
                const double coordinate = horizontal ? point.x : point.y;
                return profile.valueAt(coordinate);
            }});
    }
    return result;
}

void configureLoads(const std::vector<Vector2D>& contacts, const std::vector<HydroBoundaryPatch>& patches,
    Vector2D lower, Vector2D upper,
    HydroInterfaceLoadGeometry::Options& options)
{
    options.boundaryContacts = contacts;
    // Every unpatched physical normal row is also constrained to zero by
    // Hydro. Exclude precisely the pressure-outlet intervals, not whole sides.
    // This lets the regularized layer straddle a supply/lip/wall patch junction.
    for (const auto side : {HydroBoundarySide::West, HydroBoundarySide::South,
                           HydroBoundarySide::East, HydroBoundarySide::North})
    {
        if (options.axisymmetric && lower.x == 0.0 && side == HydroBoundarySide::West) continue;
        const bool horizontal = side == HydroBoundarySide::South || side == HydroBoundarySide::North;
        double start = horizontal ? lower.x : lower.y;
        const double end = horizontal ? upper.x : upper.y;
        std::vector<std::pair<double, double>> outlets;
        for (const auto& patch : patches)
            if (patch.side == side && patch.type == HydroBoundaryType::PressureOutlet)
                outlets.emplace_back(patch.start, patch.end);
        std::sort(outlets.begin(), outlets.end());
        for (const auto& outlet : outlets)
        {
            if (start < outlet.first) options.normalVelocityBoundaries.push_back({sideOf(side), start, outlet.first});
            start = outlet.second;
        }
        if (start < end) options.normalVelocityBoundaries.push_back({sideOf(side), start, end});
    }
}


double negativePhaseOutflow(const QuadtreeMesh& mesh,
    const LevelSetHermiteGeometry::Snapshot& scalar, const VEdgeField<double>& ux,
    const HEdgeField<double>& uy, bool axisymmetric)
{
    long double result = 0.0L;
    const auto add = [&](const Edge& edge, double velocity) {
        if (edge.masterCell != -1 && edge.slaveCell != -1) return;
        const double sign = edge.masterCell == -1 ? 1.0 : -1.0;
        result += sign*velocity*InterfaceBoundaryGeometry::negativeFaceMeasure(mesh, edge, scalar, axisymmetric);
    };
    for (int i = 0; i < mesh.vEdgeIndexNum(); ++i)
        add(mesh.vEdges(mesh.vEdgeIndexToLabel[i]), ux.value[i]);
    for (int i = 0; i < mesh.hEdgeIndexNum(); ++i)
        add(mesh.hEdges(mesh.hEdgeIndexToLabel[i]), uy.value[i]);
    if (!std::isfinite(result)) throw std::runtime_error("nonfinite negative-phase boundary flux");
    return static_cast<double>(result);
}
}
