#pragma once

#include <memory>
#include <string>

namespace Solvers { struct SolversConfig; }

class UnitsContext;
class NormalizationContext;
struct PlasmaCollisionInputs;

// CLI startup intent. The selected App owns support, confirmation and cleanup.
enum class RunStartMode { Automatic, Ask, New, Resume };

// ----------------------------------------------------------------------------
// RunConfig
//
// Runtime configuration parsed from system/run.dict.
// 发布级规范：run.dict 只负责选择“要运行的 App”。
//
// In addition, runCase loads ./system/solvers.dict once and attaches the
// parsed solver-layer configuration here, so that apps/modules share a single
// consistent view of solver knobs.
// ----------------------------------------------------------------------------

struct RunConfig
{
    std::string app;
    int schemaVersion = 1; // Omitted version retains the legacy input path.

    // Validate the selected Hydro/EHD case without constructing a mesh or advancing.
    bool checkCase = false;
    RunStartMode startMode = RunStartMode::Automatic;

    // Shared solver-layer configuration (may be null if solvers.dict is absent).
    std::shared_ptr<const Solvers::SolversConfig> solvers;
    // Shared units context (never null; defaults to CodeUnits).
    std::shared_ptr<const UnitsContext> units;

    // Optional normalization context (never null; defaults to legacy).
    std::shared_ptr<const NormalizationContext> normalization;


    // Plasma owns collision parsing, validation and immutable services.
    std::shared_ptr<const PlasmaCollisionInputs> plasmaCollisions;

    // If true, selected outputs may write both code- and SI-scaled values.
    // The SI scaling is based on UnitsContext reference scales.
    bool writeSIOutputs = false;
};
