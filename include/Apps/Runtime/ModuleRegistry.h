#pragma once

#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "Runtime/IModule.h"
#include "Runtime/RunConfig.h"

// ----------------------------------------------------------------------------
// ModuleRegistry
//
// Lightweight factory registry to create modules by name.
// Keys are case-insensitive.
// ----------------------------------------------------------------------------

class ModuleRegistry
{
public:
    using Factory = std::function<std::unique_ptr<IModule>(const RunConfig&)>;

    void registerModule(const std::string& name, Factory factory);

    std::unique_ptr<IModule> create(const std::string& name, const RunConfig& cfg) const;

    std::vector<std::string> registeredNames() const;

    static std::string normalizeName(std::string name);

private:
    std::unordered_map<std::string, Factory> factories_;
};
