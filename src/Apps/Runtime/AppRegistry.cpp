#include "Runtime/AppRegistry.h"

#include <algorithm>
#include <cctype>
#include <iostream>

using namespace std;

void AppRegistry::registerApp(const std::string& name, Factory factory)
{
    factories_[normalizeName(name)] = std::move(factory);
}

std::unique_ptr<IApp> AppRegistry::create(const std::string& name, const RunConfig& cfg) const
{
    const auto key = normalizeName(name);
    auto it = factories_.find(key);
    if (it == factories_.end())
        return nullptr;
    return (it->second)(cfg);
}

std::vector<std::string> AppRegistry::registeredNames() const
{
    std::vector<std::string> keys;
    keys.reserve(factories_.size());
    for (const auto& kv : factories_)
        keys.push_back(kv.first);
    std::sort(keys.begin(), keys.end());
    return keys;
}

std::string AppRegistry::normalizeName(std::string name)
{
    std::transform(name.begin(), name.end(), name.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return name;
}

int AppRegistry::run(const RunConfig& cfg) const
{
    if (cfg.app.empty())
        return 1;

    unique_ptr<IApp> app = create(cfg.app, cfg);
    if (!app)
    {
        cerr << "[RunController] ERROR: Unknown app: '" << cfg.app << "'.";
        const auto names = registeredNames();
        if (!names.empty())
        {
            cerr << " Available:";
            for (const auto& k : names)
                cerr << " " << k;
        }
        cerr << endl;
        return 2;
    }

    cout << "[RunController] Running app: " << app->name() << endl;
    const int rc = app->run();
    if (rc != 0)
        cerr << "[RunController] App '" << app->name() << "' returned rc=" << rc << endl;
    return rc;
}

