#pragma once

#include <map>
#include <string>
#include <utility>

// Shared SI material data; no geometry, solver, transport or App ownership.
struct PhaseMaterial
{
    int phaseId = -1; // Geometry's positive 1-based identity.
    double relativePermittivity = 0.0;
    double conductivitySI = 0.0;
    double densitySI = 0.0;
    double dynamicViscositySI = 0.0;
};

struct MaterialTable
{
    std::map<int, PhaseMaterial> phases;
    std::map<std::pair<int, int>, double> surfaceTensionSI;
};
