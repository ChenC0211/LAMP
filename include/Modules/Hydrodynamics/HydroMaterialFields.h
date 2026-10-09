#pragma once

#include <vector>
#include <cstdint>

class QuadtreeMesh;
class HydroFittedDomain;
template<class T> class CellField;
template<class T> class VertField;

struct HydroPhaseProperties
{
    int phaseId;
    double density;
    double viscosity;
};

enum class HydroMaterialEvolution { Frozen, TimeDependent };

// Module-owned provenance of a successfully completed material update.
// It is not an interface topology version or coupled-step acceptance token.
struct HydroMaterialState
{
    std::uint64_t revision = 0;
    double sourceTime = 0.0;
    HydroMaterialEvolution evolution = HydroMaterialEvolution::Frozen;
};

// Writable material destination owned by Hydro or its explicit AMR staging
// object. This is not permission to mutate velocity/pressure/history fields.
struct HydroMaterialFields
{
    CellField<double>& rho;
    CellField<double>& mu;
    VertField<double>& vertRho;
    VertField<double>& vertMu;
};

// Reuses the original active-cell assignment and active-only arithmetic
// vertex trace. Values are already in solver code units; phase classification
// and SI conversion belong to the caller. Inactive cell values are retained.
// No NS initialization, mesh construction, interface classification or IO.
void assignHydroFittedMaterials(const QuadtreeMesh& mesh, const HydroFittedDomain& domain,
    HydroMaterialFields fields, const std::vector<int>& phaseByCell,
    const HydroPhaseProperties& first, const HydroPhaseProperties& second);
