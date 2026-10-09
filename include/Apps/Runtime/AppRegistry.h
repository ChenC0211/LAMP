#pragma once

#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "Runtime/IApp.h"
#include "Runtime/RunConfig.h"

// ----------------------------------------------------------------------------
// AppRegistry
//
// Lightweight factory registry to create apps by name.
// Keys are case-insensitive.
// ----------------------------------------------------------------------------

class AppRegistry
{
public:
    using Factory = std::function<std::unique_ptr<IApp>(const RunConfig&)>;

    void registerApp(const std::string& name, Factory factory);

    std::unique_ptr<IApp> create(const std::string& name, const RunConfig& cfg) const;

    // Execute one registered factory with the caller's prepared configuration.
    // Empty selection returns 1, unknown/null factory result returns 2;
    // otherwise preserve the App's return code and exception behavior.
    int run(const RunConfig& cfg) const;

    std::vector<std::string> registeredNames() const;

    static std::string normalizeName(std::string name);

private:
    std::unordered_map<std::string, Factory> factories_;
};
