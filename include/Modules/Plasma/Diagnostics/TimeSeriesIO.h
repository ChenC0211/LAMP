#pragma once

#include <cctype>
#include <cerrno>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>
#include <system_error>
#include <vector>

namespace TimeSeriesIO
{

inline std::string sanitizeCsvToken(const std::string& raw)
{
    std::string out;
    out.reserve(raw.size());
    bool prevUnderscore = false;
    for (const unsigned char ch : raw)
    {
        if (std::isalnum(ch))
        {
            out.push_back(static_cast<char>(std::tolower(ch)));
            prevUnderscore = false;
        }
        else if (!prevUnderscore)
        {
            out.push_back('_');
            prevUnderscore = true;
        }
    }
    while (!out.empty() && out.front() == '_')
        out.erase(out.begin());
    while (!out.empty() && out.back() == '_')
        out.pop_back();
    if (out.empty())
        out = "species";
    return out;
}

inline std::string trimCopy(const std::string& value)
{
    std::size_t first = 0;
    while (first < value.size() && std::isspace(static_cast<unsigned char>(value[first])) != 0)
        ++first;
    std::size_t last = value.size();
    while (last > first && std::isspace(static_cast<unsigned char>(value[last - 1])) != 0)
        --last;
    return value.substr(first, last - first);
}

inline bool tryParseLeadingStep(const std::string& line, int& stepOut)
{
    const std::size_t comma = line.find(',');
    const std::string token = trimCopy(line.substr(0, comma));
    if (token.empty())
        return false;

    char* end = nullptr;
    errno = 0;
    const long parsed = std::strtol(token.c_str(), &end, 10);
    if (errno != 0 || end == nullptr || *end != '\0')
        return false;
    stepOut = static_cast<int>(parsed);
    return true;
}

inline bool prepareStepCsvFile(const std::filesystem::path& filePath, int startStep)
{
    std::error_code ec;
    if (filePath.has_parent_path())
        std::filesystem::create_directories(filePath.parent_path(), ec);

    if (startStep <= 0)
    {
        std::filesystem::remove(filePath, ec);
        return false;
    }

    if (!std::filesystem::exists(filePath, ec))
        return false;

    std::ifstream in(filePath, std::ios::in);
    if (!in)
        return false;

    std::vector<std::string> keptLines;
    std::string line;
    while (std::getline(in, line))
    {
        int rowStep = -1;
        if (tryParseLeadingStep(line, rowStep) && rowStep >= startStep)
            continue;
        keptLines.push_back(line);
    }
    in.close();

    if (keptLines.empty())
    {
        std::filesystem::remove(filePath, ec);
        return false;
    }

    std::ofstream out(filePath, std::ios::out | std::ios::trunc);
    if (!out)
        return false;
    for (const auto& kept : keptLines)
        out << kept << '\n';
    out.close();
    return true;
}

inline bool openPreparedStepCsv(std::ofstream& out,
                                const std::filesystem::path& filePath,
                                int startStep)
{
    const bool hasExistingContent = prepareStepCsvFile(filePath, startStep);
    out.open(filePath.string(), std::ios::out | std::ios::app);
    return hasExistingContent && static_cast<bool>(out);
}

inline bool tryParseTrailingStepToken(const std::filesystem::path& filePath, int& stepOut)
{
    const std::string stem = filePath.stem().string();
    if (stem.empty())
        return false;

    const std::size_t us = stem.rfind('_');
    const std::string token = (us == std::string::npos) ? stem : stem.substr(us + 1);
    if (token.empty())
        return false;

    char* end = nullptr;
    errno = 0;
    const long parsed = std::strtol(token.c_str(), &end, 10);
    if (errno != 0 || end == nullptr || *end != '\0')
        return false;
    stepOut = static_cast<int>(parsed);
    return stepOut >= 0;
}

inline void prepareStepTaggedDirectory(const std::filesystem::path& dirPath, int startStep)
{
    std::error_code ec;
    if (startStep <= 0)
    {
        std::filesystem::remove_all(dirPath, ec);
        std::filesystem::create_directories(dirPath, ec);
        return;
    }

    if (!std::filesystem::exists(dirPath, ec))
    {
        std::filesystem::create_directories(dirPath, ec);
        return;
    }

    for (const auto& entry : std::filesystem::recursive_directory_iterator(dirPath, ec))
    {
        if (ec)
            break;
        if (!entry.is_regular_file())
            continue;

        int fileStep = -1;
        if (!tryParseTrailingStepToken(entry.path(), fileStep))
            continue;
        if (fileStep >= startStep)
            std::filesystem::remove(entry.path(), ec);
    }
}

} // namespace TimeSeriesIO
