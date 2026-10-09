#pragma once

#include <filesystem>

// Case-relative output categories, matching the Plasma layout. File format
// does not decide the directory: DAT/CSV/TSV spatial frames all belong in
// fields, time histories in series, and manifests/descriptions in meta.
// Only creates the requested directory; never cleans or relocates old output.
namespace OutputLayout
{
inline std::filesystem::path directory(const char* category)
{
    const auto path = std::filesystem::path("./output") / category;
    std::filesystem::create_directories(path);
    return path;
}

inline std::filesystem::path fields() { return directory("fields"); }
inline std::filesystem::path series() { return directory("series"); }
inline std::filesystem::path meta() { return directory("meta"); }
}
