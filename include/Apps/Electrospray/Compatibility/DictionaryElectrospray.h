#pragma once

#include <string>

#include "Core/IO/Dictionary.h"
#include "Electrospray/Compatibility/ElectrosprayConfig.h"

// Reads app-level runtime settings for ElectrosprayApp.
// 发布级规范：Electrospray 的驱动配置应与 app 选择合并到 ./system/run.dict。
// Default file: ./system/run.dict (Electrospray{ ... }).
class DictionaryElectrospray : public Dictionary
{
public:
    DictionaryElectrospray() = default;
    ~DictionaryElectrospray() = default;

    // Return false if file cannot be opened, or if Electrospray block is not found.
    bool readDict(const std::string& path = "./system/run.dict");

    const ElectrosprayConfig& config() const { return cfg_; }

private:
    ElectrosprayConfig cfg_{};

    static std::string toLowerCopy(std::string s);
    static bool parseBoolToken(const std::string& s, bool& out);
};
