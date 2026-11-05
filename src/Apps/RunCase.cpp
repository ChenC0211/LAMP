#include "RunCase.h"

#include <iostream>
#include <memory>

#include "Core/IO/StrictDictionary.h"
#include "AppRegistry.h"
#include "IO/DictionaryParallel.h"
#include "IO/DictionaryRun.h"

#include "Electric/ElectricApp.h"
#include "Hydrodynamics/HydrodynamicsApp.h"
#include "Electrohydrodynamics/VofElectrohydrodynamicsApp.h"
#include "Pic/PicApp.h"

namespace
{
AppRegistry builtInApps()
{
    AppRegistry apps;
    apps.registerApp("Electric", [](const RunConfig& c) { return std::make_unique<ElectricApp>(c); });
    apps.registerApp("Hydrodynamics", [](const RunConfig& c) { return std::make_unique<HydrodynamicsApp>(c); });
    apps.registerApp("VofElectrohydrodynamics", [](const RunConfig& c) { return std::make_unique<VofElectrohydrodynamicsApp>(c,
        VofElectrohydrodynamicsInput::Application::General); });
    apps.registerApp("Electrospray", [](const RunConfig& c) { return std::make_unique<VofElectrohydrodynamicsApp>(c,
        VofElectrohydrodynamicsInput::Application::Electrospray); });
    apps.registerApp("Pic", [](const RunConfig& c) { return std::make_unique<PicApp>(c); });
    return apps;
}
} // namespace

int runCase(bool checkCase, std::optional<int> threadsOverride)
{
    DictionaryRun dictRun;
    if (!dictRun.readDict("./system/run.dict"))
    {
        std::cerr << "[LAMP] ERROR: Cannot read ./system/run.dict. Expected:\n"
                  << "        Run { app <AppName> }\n";
        return 2;
    }

    RunConfig cfg = dictRun.config();
    if (cfg.app.empty())
    {
        std::cerr << "[LAMP] ERROR: Empty app selection in ./system/run.dict. Expected:\n"
                  << "        Run { app <AppName> }\n";
        return 2;
    }
    cfg.checkCase = checkCase;
    cfg.threadsOverride = threadsOverride;

    // Process resources are configured once by the shared entry before dispatch.
    const DictionaryInput::Document run("./system/run.dict");
    const int threads = AppRuntime::readParallelThreads(run.root(), nullptr, threadsOverride);
    AppRuntime::configureParallelThreads(threads, threadsOverride ? "command_line"
        : run.root().find("Parallel") ? "run.dict" : "default");

    return builtInApps().run(cfg);
}
