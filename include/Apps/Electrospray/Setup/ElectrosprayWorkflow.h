#pragma once

#include <optional>
#include <stdexcept>
#include <string_view>

// Closed App compositions of the same module services. The fixed-interface
// and pressure-only problems are verification runs, not separate NS models.
enum class ElectrosprayProblem
{
    Conducting, MovingDielectric, PressureDriven,
    FixedInterfaceCapillary, FixedInterfaceDielectric
};

// Stable legacy labels remain in existing output files and compatibility IO.
inline std::string_view electrosprayProblemLabel(ElectrosprayProblem problem)
{
    switch (problem)
    {
    case ElectrosprayProblem::Conducting: return "coupled";
    case ElectrosprayProblem::MovingDielectric: return "movingDielectric";
    case ElectrosprayProblem::PressureDriven: return "frozenMaterialHydro";
    case ElectrosprayProblem::FixedInterfaceCapillary: return "frozenInterfaceCapillary";
    case ElectrosprayProblem::FixedInterfaceDielectric: return "frozenInterfaceEhd";
    }
    throw std::invalid_argument("invalid Electrospray problem");
}

// Formal input has already validated these electric models. Keep its actual
// problem selection shared with resolved configuration output.
inline ElectrosprayProblem electrosprayProblemForElectricModel(std::string_view model)
{
    if (model == "leakyDielectric") return ElectrosprayProblem::Conducting;
    if (model == "perfectDielectric") return ElectrosprayProblem::MovingDielectric;
    throw std::invalid_argument("unsupported Electrospray electric model");
}

// App composition and initial-geometry choices only. Time/output/log controls
// are owned by RunControls; phase properties, force width and solver controls
// remain in their respective module inputs. No numerical state is stored here.
struct ElectrosprayWorkflow
{
    ElectrosprayProblem problem = ElectrosprayProblem::Conducting;
    bool conducting() const noexcept { return problem == ElectrosprayProblem::Conducting; }
    bool movesInterface() const noexcept { return conducting() || problem == ElectrosprayProblem::MovingDielectric; }
    bool usesDielectricPolarization() const noexcept {
        return problem == ElectrosprayProblem::MovingDielectric || problem == ElectrosprayProblem::FixedInterfaceDielectric;
    }
    int restartLoop = -1;
    int initialAdaptIterations = -1;
    bool initializeLevelSets = true;
    // Old inputs could select iteration progress independently of geometry
    // detail. Only the one-way legacy adapter sets this override; current
    // input uses RunControls::detailedLog for both, with no new dictionary key.
    std::optional<bool> legacyIterationLog;
};
