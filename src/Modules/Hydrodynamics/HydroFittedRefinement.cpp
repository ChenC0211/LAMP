#include "Modules/Hydrodynamics/HydroFittedRefinement.h"

#include "Core/Mesh/QuadtreeMesh.h"
#include "Modules/Hydrodynamics/HydroFittedDomain.h"
#include "Modules/Hydrodynamics/HydroFittedVelocity.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

void protectHydroFittedBoundarySupport(const QuadtreeMesh& mesh,
    const HydroFittedDomain& domain, std::vector<std::uint8_t>& protectedCells)
{
    domain.validateMeshIdentity(mesh);
    if (protectedCells.size() != static_cast<std::size_t>(mesh.cellIndexNum()))
        throw std::invalid_argument("Hydro mesh protection requires a native cell mask");
    for (const auto& f : domain.faces())
    {
        if ((f.negativeCell < 0 && f.positiveCell < 0)
            || f.boundary.type == HydroFittedDomain::BoundaryType::Interior) continue;
        const auto& edge = f.horizontal ? mesh.hEdges(f.edgeLabel) : mesh.vEdges(f.edgeLabel);
        for (int label : {edge.masterCell, edge.slaveCell})
            if (label >= 0) protectedCells[mesh.cells(label).index()] = 1;
    }
}

HydroFittedRefinementFields evaluateHydroFittedRefinement(const QuadtreeMesh& mesh,
    const HydroFittedDomain& domain, const HydroFittedVelocitySampler& velocity)
{
    domain.validateMeshIdentity(mesh);
    HydroFittedRefinementFields result;
    result.speed.assign(mesh.cellIndexNum(), 0);
    result.velocityGradient.assign(mesh.cellIndexNum(), 0);
    for (int i : domain.cellIndexByRow())
    {
        const auto& cell = mesh.cells(mesh.cellIndexToLabel[i]);
        const auto center = cell.centerCoord(), size = cell.size();
        const auto v = velocity.sample(center);
        // Interior secants of the existing, wall-aware velocity polynomial.
        // A quarter-cell offset avoids crossing a phase/solid mesh boundary.
        const auto west = velocity.sample(center - Vector2D(.25 * size.x, 0));
        const auto east = velocity.sample(center + Vector2D(.25 * size.x, 0));
        const auto south = velocity.sample(center - Vector2D(0, .25 * size.y));
        const auto north = velocity.sample(center + Vector2D(0, .25 * size.y));
        const auto dr = (east - west) * (2 / size.x);
        const auto dz = (north - south) * (2 / size.y);
        result.speed[i] = std::max({v.mod(), west.mod(), east.mod(), south.mod(), north.mod()});
        // RZ without swirl: include the azimuthal u_r/r component of grad(u).
        // Cell centers lie strictly off axis; the shared axis trace is regular.
        const double hoop = domain.axisymmetric() ? v.x / center.x : 0;
        result.velocityGradient[i] = std::hypot(std::hypot(dr.mod(), dz.mod()), hoop);
        if (!std::isfinite(result.speed[i]) || !std::isfinite(result.velocityGradient[i]))
            throw std::runtime_error("Hydro AMR: nonfinite accepted velocity indicator");
    }
    return result;
}
