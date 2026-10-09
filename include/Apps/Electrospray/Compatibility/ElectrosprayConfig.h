#pragma once

#include <string>
#include <optional>

// ----------------------------------------------------------------------------
// ElectrosprayConfig
//
// App-level runtime configuration. Physical numbers and material/boundary
// choices are read separately by ElectrosprayPhysicalConfig.
//
// Runtime and numerical discretization controls only. Physical/model
// parameters remain in their models.
// ----------------------------------------------------------------------------

struct ElectrosprayConfig
{
    // Full coupled spray remains the default. Frozen development modes retain
    // initial materials/interface. movingDielectric updates geometry/material/E
    // at q=0 but omits conductivity; none replaces conducting electrospray.
    std::string physicsMode = "coupled";
    // Numerical compact C2 force-band half-width / max(leaf dr, leaf dz).
    // This is neither a physical diffuse interface nor a thin-limit claim.
    double capillaryHalfWidthCells = 3.0;
    // Format-independent layout: output/fields, output/series, output/meta.
    std::string outputName = "electrospray_init";

    // Write Tecplot vertex fields for quick verification.
    bool writeFields = true;
    // Legacy inputs retain their original output inventory. Schema 2 binds
    // Output.Fields.profile independently of the console Log.level.
    bool diagnosticFields = true;
    bool writeInitial = true;
    bool writeFinal = true;
    std::optional<bool> writeFirst; // absent retains legacy moving/frozen behavior
    int writeEverySteps = 0; // 0: legacy code-time interval
    bool detailedIterationLog = true; // legacy volume progress; schema 2 uses Log.level
    int logEvery = 1; // Always also report the first and final accepted step.
    bool verboseDiagnostics = false; // Detailed geometry/force diagnostics.

    // If >= 0, restart from ./cache/cache.<restartLoop>/
    int restartLoop = -1;

    // Driver timestep and times (kept for compatibility with the legacy driver).
    double dt = 5e-6;
    double finalTime = 20.0;
    double writeDt = 0.1;

    // Override initial adapt iterations (if <0, use mesh.maxInitialIterations()).
    int initialAdaptIterations = -1;

    // Initialize from the case geometry; no hard-coded replacement shape.
    bool initializeLevelSets = true;
};
