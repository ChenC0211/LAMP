#pragma once
#include <filesystem>
#include <fstream>

namespace FileIO {
std::ofstream openOutputFile(const std::filesystem::path& path, std::ios::openmode mode);
bool closeOutputFile(std::ofstream& stream);
void finishOutputFile(std::ofstream& stream, const std::filesystem::path& path);
}
