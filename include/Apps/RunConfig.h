#pragma once

#include <optional>
#include <string>

// Process-level request handed to every App; each App reads its own dictionaries.
struct RunConfig
{
    std::string app;

    // Validate supported App inputs without constructing a mesh or advancing.
    bool checkCase = false;

    // Explicit command-line resource override, consumed by AppRuntime.
    std::optional<int> threadsOverride;
};
