#pragma once
#include <filesystem>
#include <fstream>
#include <string>

// Materialize an already extracted block for legacy file-based readers.
// Selection, dictionary aliases and cache naming belong to the caller.
inline bool writeDictionaryBlockFile(const std::filesystem::path& path,
                                     const std::string& blockName,
                                     const std::string& blockText)
{
    std::error_code error;
    std::filesystem::create_directories(path.parent_path(), error);
    std::ofstream out(path, std::ios::out | std::ios::trunc);
    if (!out) return false;
    out << blockName << " {\n" << blockText << "\n}\n";
    out.close();
    return !out.fail();
}
