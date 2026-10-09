#include "Core/IO/FileManagement.h"
#include <stdexcept>

namespace FileIO {
std::ofstream openOutputFile(const std::filesystem::path& path, std::ios::openmode mode)
{
    std::ofstream stream(path, mode);
    if (!stream) throw std::runtime_error("Cannot open output file: " + path.string());
    return stream;
}

bool closeOutputFile(std::ofstream& stream)
{
    stream.flush();
    const bool written = static_cast<bool>(stream);
    stream.close();
    return written && static_cast<bool>(stream);
}

void finishOutputFile(std::ofstream& stream, const std::filesystem::path& path)
{
    if (!closeOutputFile(stream))
        throw std::runtime_error("Failed to write or close output file: " + path.string());
}
}
