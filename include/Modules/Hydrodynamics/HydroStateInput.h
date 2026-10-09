#pragma once

#include "Modules/Hydrodynamics/HydroStep.h"
#include "Modules/Hydrodynamics/HydroMaterialFields.h"
#include <variant>

class LevelSet;
class NativeLeafDomain;
class InterfaceConservativePhase;
class InterfaceState;

// Prescribed complete-domain profiles, in code velocity units. Smooth is the
// original non-solenoidal trigonometric seed; it is not a solved flow state.
enum class HydroInitialVelocity { Zero, Smooth, Vortex };

// Native fields sampled by input/setup code, in solver code units. Null means
// zero; the velocity components must be supplied together. No time, history,
// derived diagnostic or analytical pressure recovery is supplied by the App.
struct HydroInitialFields
{
    const CellField<double>* pressure = nullptr;
    const VEdgeField<double>* velocityX = nullptr;
    const HEdgeField<double>* velocityY = nullptr;
};

struct HydroConstantMaterial
{
    double density, viscosity;
};

// Prescribed spatial coefficients on the existing native locations. These
// are effective code-unit material fields, not a second mixture law. Borrowed
// only during updateMaterials; the module copies them before publication.
struct HydroNativeMaterial
{
    const CellField<double>* density;
    const CellField<double>* viscosity;
    const VertField<double>* vertexDensity;
    const VertField<double>* vertexViscosity;
};

// The existing vertex tanh interpolation and vertex-to-cell average, not a
// new mixture law. The borrowed scalar is the actual source interface.
struct HydroLevelSetMaterial
{
    const LevelSet* levelSet;
    double densityNegative, densityPositive, viscosityNegative, viscosityPositive, thickness;
};

// The existing sharp cell classification and active-only vertex trace.
// Classification and SI conversion stay with their actual providers.
struct HydroClassifiedMaterial
{
    const std::vector<int>* phaseByCell;
    HydroPhaseProperties first, second;
};

// Complete native XY/RZ: cell density is mass/capacity from the authoritative
// liquid inventory C, with arithmetic phase viscosity. Vertex coefficients
// are traces of these same cell mixtures, not a second level-set mixture.
// This call updates coefficients only. HydroStepInput::conservativeTransport
// supplies the matching dual masses/fluxes separately. Domain/inventory are borrowed synchronously.
struct HydroConservativeMaterial
{
    const NativeLeafDomain* domain;
    const InterfaceConservativePhase* inventory;
    double densityNegative, densityPositive, viscosityNegative, viscosityPositive;
};
// Same conservative mixture, read directly from the parameter interface's
// independent FV inventory. No legacy scalar/phase state is reconstructed.
struct HydroParametricMaterial
{
    const NativeLeafDomain* domain;
    const InterfaceState* inventory;
    double densityNegative, densityPositive, viscosityNegative, viscosityPositive;
};

struct HydroMaterialUpdate
{
    HydroStepGeometry geometry;
    double sourceTime; // current completed Hydro substep time, not a future step
    HydroMaterialEvolution evolution;
    std::variant<HydroConstantMaterial, HydroNativeMaterial, HydroLevelSetMaterial, HydroClassifiedMaterial,
        HydroConservativeMaterial, HydroParametricMaterial> model;
};
