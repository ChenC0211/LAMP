#include "Modules/Hydrodynamics/HydroMaterialFields.h"
#include "Modules/Hydrodynamics/HydroFittedDomain.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Field/CellField.h"
#include "Core/Field/VertField.h"

#include <cmath>
#include <stdexcept>

void assignHydroFittedMaterials(const QuadtreeMesh& mesh, const HydroFittedDomain& domain,
    HydroMaterialFields hydro, const std::vector<int>& phaseByCell,
    const HydroPhaseProperties& first, const HydroPhaseProperties& second)
{
    domain.validateMeshIdentity(mesh);
    if (first.phaseId == second.phaseId
        || !std::isfinite(first.density) || first.density <= 0.0
        || !std::isfinite(second.density) || second.density <= 0.0
        || !std::isfinite(first.viscosity) || first.viscosity < 0.0
        || !std::isfinite(second.viscosity) || second.viscosity < 0.0)
        throw std::invalid_argument("Hydro requires distinct phases, positive densities and nonnegative viscosities");
    for (const auto* field : {&hydro.rho, &hydro.mu})
    {
        if (field->size() != mesh.cellIndexNum() || !field->value)
            throw std::invalid_argument("Hydro material destination does not match mesh");
        for (int i = 0; i < field->size(); ++i)
            if (field->label(i) != mesh.cellIndexToLabel[i])
                throw std::invalid_argument("Hydro material destination has stale native labels");
    }
    if (phaseByCell.size() != static_cast<std::size_t>(mesh.cellIndexNum()))
        throw std::invalid_argument("Hydro material field does not cover the actual mesh");
    hydro.vertRho.initialize(mesh, 0.0); hydro.vertMu.initialize(mesh, 0.0);
    for (int i : domain.cellIndexByRow())
    {
        if (phaseByCell[i] != first.phaseId && phaseByCell[i] != second.phaseId)
            throw std::invalid_argument("active fluid cell has an unbound material phase");
        const auto& material = phaseByCell[i] == first.phaseId ? first : second;
        hydro.rho.value[i] = material.density;
        hydro.mu.value[i] = material.viscosity;
    }
    // Same active-only material trace used by the original frozen loop.
    for (int i = 0; i < mesh.vertIndexNum(); ++i)
    {
        const auto& vertex = mesh.verts(mesh.vertIndexToLabel[i]);
        int count = 0;
        for (int label : vertex.vertToCell)
        {
            if (label < 0) continue;
            const int cell = mesh.cells(label).index();
            if (domain.cellRowByIndex()[cell] < 0) continue;
            hydro.vertRho.value[i] += hydro.rho.value[cell];
            hydro.vertMu.value[i] += hydro.mu.value[cell];
            ++count;
        }
        if (count > 0) { hydro.vertRho.value[i] /= count; hydro.vertMu.value[i] /= count; }
    }
}
