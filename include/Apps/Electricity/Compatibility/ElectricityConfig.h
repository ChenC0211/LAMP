#pragma once

#include <string>

// ----------------------------------------------------------------------------
// ElectricityConfig
//
// Historical runtime configuration for SingleModule/Electricity only.
//
// Notes:
// - This config is intentionally small. Algorithmic and physical settings
//   belong in the corresponding models; module dicts control *run* behavior.
// - Defaults reproduce the current Phase-2 module behavior.
// ----------------------------------------------------------------------------

struct ElectricityConfig
{
    // Prefix relative to output/. Explicit legacy prefixes retain their paths.
    std::string outputName = "fields/electricity_init";

    // Whether to write Tecplot vertex fields for quick verification.
    bool writeFields = true;

    // Deprecated compatibility input, parsed with the original aliases/defaults.
    // The former temporary Hydro owner discarded its computed material fields;
    // neither these values nor the selected level set affected electricity.
    // Retaining parsing preserves old malformed-integer rejection as well.
    bool computeHydroProperties = true;

    // No longer selects a discarded material calculation. No mesh/field effect.
    int propertyLevelSetIndex = 1;
};
