#include "Hydrodynamics/IO/HydroCaseInput.h"
#include "Modules/Hydrodynamics/IO/DictionaryHydrodynamics.h"
#include "Modules/Interface/Compatibility/IO/DictionaryInterface.h"
#include "Runtime/IO/DictionaryRunControls.h"
#include "Core/IO/DictionaryTokenReader.h"
#include "Core/IO/ConfigurationTrace.h"
#include "Core/Units/UnitsContext.h"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <limits>
#include <set>
#include <stdexcept>
#include <utility>

bool usesSeparateHydroInputs()
{
    DictionaryInput::Reader input("./system/run.dict");
    int depth = 0;
    while (!input.done())
    {
        const auto token = input.next();
        if (depth == 0 && (token == "Time" || token == "Output" || token == "Log" || token == "Restart"))
            return true; // Malformed blocks must reach the strict reader too.
        if (token == "{") ++depth;
        else if (token == "}") --depth;
    }
    return false;
}

HydroCaseInput readHydroCaseInput(bool coupled)
{
    HydroCaseInput result;
    auto schedule = std::make_shared<HydroCaseSchedule>();
    if (!usesSeparateHydroInputs())
    {
        if (std::filesystem::exists("./system/interface.dict"))
            throw std::runtime_error("interface.dict requires Time/Output/Log in run.dict; do not mix legacy and separate inputs");
        DictionaryHydrodynamics legacy;
        if (!legacy.readDict()) throw std::runtime_error("cannot open ./system/hydro.dict");
        result.flow = legacy.config();
        result.interface = legacy.interfaceConfig();
        const auto& input = legacy.runInput();
        schedule->steps = input.steps;
        auto& run = schedule->controls;
        run.dt = input.dt;
        run.endTime = input.steps * input.dt;
        run.writeEverySteps = input.writeEvery;
        run.logEverySteps = input.logEvery;
        run.outputPrefix = input.outputName;
        run.writeFields = input.writeFields;
        run.writeInitial = input.writeInitialFields;
        run.writeFirst = false;
        run.writeFinal = true;
    }
    else
    {
        std::set<std::string> allowed{"run.dict", "hydro.dict", "interface.dict", "geometry.dict", "mesh.dict"};
        if (coupled) allowed.insert("electricity.dict");
        for (const auto& entry : std::filesystem::directory_iterator("./system"))
            if (!entry.is_regular_file() || !allowed.count(entry.path().filename().string()))
                throw std::runtime_error("unsupported separate Hydro/EHD input: " + entry.path().string());
        DictionaryInput::Document document("./system/run.dict");
        const auto& root = document.root();
        if (coupled) root.allow({"Run", "Time", "Output", "Log", "Restart"});
        else root.allow({"Run", "SingleModule", "Time", "Output", "Log", "Restart"});
        const auto& header = root.block("Run");
        header.allow({"app"});
        auto app = header.get("app").word();
        std::transform(app.begin(), app.end(), app.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        if (coupled)
        {
            if (app != "electrohydrodynamics" && app != "electrohydrodynamicsapp" && app != "ehd")
                header.get("app").fail("Electrohydrodynamics, ElectrohydrodynamicsApp or EHD");
        }
        else
        {
            if (app != "singlemodule" && app != "singlemoduleapp") header.get("app").fail("SingleModule");
            const auto& module = root.block("SingleModule");
            module.allow({"modules"});
            const auto& names = module.get("modules");
            if (names.values.size() != 3 || names.values.front() != "(" || names.values.back() != ")")
                names.fail("(Hydrodynamics)");
            auto name = names.values[1];
            std::transform(name.begin(), name.end(), name.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
            if (name != "hydrodynamics" && name != "hydro") names.fail("(Hydrodynamics)");
        }
        root.block("Time").get("units").choice({"code"});
        DictionaryInput::ConfigurationTrace trace;
        schedule->controls = readRunControls(root, UnitsContext::MakeCodeUnits(), trace,true);
        schedule->steps = schedule->controls.stepRetry?static_cast<int>(std::ceil(schedule->controls.endTime/schedule->controls.dt))
            :schedule->controls.fixedStepCount().value();
        if (schedule->steps == std::numeric_limits<int>::max())
            root.block("Time").get("endTime").fail("fewer than INT_MAX steps");
        if (schedule->controls.diagnosticFields)
            root.block("Output").block("Fields").get("profile").fail("physical; no additional field profile in ordinary Hydro/EHD");
        const bool hasInterface = std::filesystem::exists("./system/interface.dict");
        if (hasInterface) result.interface = readNodalInterfaceCase("./system/interface.dict");
        if(schedule->controls.stepRetry && result.interface.transport!="conservative")
            root.block("Time").block("StepRetry").fail("uncoupled conservative Hydro transport");
        if(schedule->controls.columnAmplitudeLimit && result.interface.transport!="conservative")
            root.block("Time").block("Stop").fail("uncoupled conservative Hydro column transport");
        result.flow = readNodalHydrodynamicsConfig("./system/hydro.dict", result.interface);
        if (schedule->controls.detailedLog && result.interface.transport != "conservative")
            root.block("Log").get("level").fail("summary; detail timing requires uncoupled conservative Hydro");
        if ((result.flow.materialModel == "twoPhase") != hasInterface)
            throw std::runtime_error("separate Hydro/EHD requires interface.dict exactly for materialModel twoPhase");
    }
    result.schedule = std::move(schedule);
    return result;
}
