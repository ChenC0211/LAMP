#pragma once

#include <string>

#include "Core/IO/Dictionary.h"
#include "Modules/Hydrodynamics/HydrodynamicsConfig.h"
#include "Modules/Interface/Compatibility/InterfaceNodalConfig.h"

// Separate-input path: hydro.dict contains flow physics only. The resolved
// Interface input is borrowed solely to validate supported combinations.
HydrodynamicsConfig readNodalHydrodynamicsConfig(const std::string& path,
    const InterfaceNodalConfig& interfaceConfig);

// Reads the existing Hydro case configuration for the App-owned runtime adapter.
// Default file: ./system/hydro.dict
class DictionaryHydrodynamics : public Dictionary
{
public:
    // Transient legacy hydro.dict run input. The App adapts this to its
    // shared RunControls; numerical consumers receive config() only.
    struct RunInput
    {
        int steps = 1;
        double dt = 0.001;
        std::string outputName = "hydro";
        bool writeFields = true, writeInitialFields = false;
        int writeEvery = 1, logEvery = 1;
    };

    DictionaryHydrodynamics() = default;
    ~DictionaryHydrodynamics() = default;

    // Return false if file cannot be opened; throw with path/line on invalid
    // schema, duplicate/unknown keys, unsupported features or invalid values.
    bool readDict(const std::string& path = "./system/hydro.dict");

    const HydrodynamicsConfig& config() const { return cfg_; }
    const RunInput& runInput() const { return run_; }
    const InterfaceNodalConfig& interfaceConfig() const { return interface_; }

private:
    HydrodynamicsConfig cfg_{};
    RunInput run_{};
    InterfaceNodalConfig interface_{};

};
