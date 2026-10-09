#pragma once

#include <string>

#include "Core/IO/Dictionary.h"
#include "Electricity/Compatibility/ElectricityConfig.h"

// Reads module-level runtime settings for ElectricityModule.
// Default file: ./system/electricity.dict
class DictionaryElectricity : public Dictionary
{
public:
    DictionaryElectricity() = default;
    ~DictionaryElectricity() = default;

    // Return false if file cannot be opened.
    bool readDict(const std::string& path = "./system/electricity.dict");

    const ElectricityConfig& config() const { return cfg_; }

private:
    ElectricityConfig cfg_{};

    static std::string toLowerCopy(std::string s);
    static bool parseBoolToken(const std::string& s, bool& out);
};
