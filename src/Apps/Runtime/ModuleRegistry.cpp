#include "Runtime/ModuleRegistry.h"

#include <algorithm>
#include <cctype>

using namespace std;

void ModuleRegistry::registerModule(const std::string& name, Factory factory)
{
    factories_[normalizeName(name)] = std::move(factory);
}

std::unique_ptr<IModule> ModuleRegistry::create(const std::string& name, const RunConfig& cfg) const
{
    const auto key = normalizeName(name);
    auto it = factories_.find(key);
    if (it == factories_.end())
        return nullptr;
    return (it->second)(cfg);
}

std::vector<std::string> ModuleRegistry::registeredNames() const
{
    std::vector<std::string> keys;
    keys.reserve(factories_.size());
    for (const auto& kv : factories_)
        keys.push_back(kv.first);
    std::sort(keys.begin(), keys.end());
    return keys;
}

std::string ModuleRegistry::normalizeName(std::string name)
{
    std::transform(name.begin(), name.end(), name.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return name;
}
