#include "Electrohydrodynamics/ElectrohydrodynamicsApp.h"
#include "Electrohydrodynamics/Coupling/DielectricCoupling.h"
#include "Modules/Electricity/DielectricElectrostatics.h"
#include "Modules/Electricity/IO/DictionaryDielectric.h"
#include "Hydrodynamics/HydroCase.h"
#include "Hydrodynamics/IO/HydroCaseInput.h"
#include "Modules/Hydrodynamics/IO/HydroInputParser.h"
#include "Core/IO/StrictDictionary.h"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <iostream>
#include <set>
#include <stdexcept>

namespace
{
bool validateAppSources()
{
    // The small code-unit App must not silently consume a legacy spray block,
    // PIC configuration, restart setting or ignored material/solver file.
    const bool separate = usesSeparateHydroInputs();
    std::set<std::string> allowed{
        "run.dict", "hydro.dict", "electricity.dict", "geometry.dict", "mesh.dict"};
    if (separate) allowed.insert("interface.dict");
    for (const auto& entry : std::filesystem::directory_iterator("./system"))
        if (!entry.is_regular_file() || !allowed.count(entry.path().filename().string()))
            throw std::runtime_error("Electrohydrodynamics supports run/hydro/electricity/geometry/mesh.dict and separate interface.dict; unsupported entry: " + entry.path().string());
    if (separate) return true; // Shared case input validates the complete run schema.
    HydroInput::Reader input("./system/run.dict");
    input.expect("Run");
    input.expect("{");
    input.expect("app");
    auto app = input.next();
    std::transform(app.begin(), app.end(), app.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (app != "electrohydrodynamics" && app != "electrohydrodynamicsapp" && app != "ehd")
        input.fail("expected Electrohydrodynamics, ElectrohydrodynamicsApp or EHD");
    input.expect("}");
    if (!input.done()) input.fail("EHD run.dict only selects the App; use hydro.dict and electricity.dict for physics, no legacy Electrospray block");
    return false;
}
}

int ElectrohydrodynamicsApp::run()
{
    try
    {
        const bool separate = validateAppSources();
        std::unique_ptr<HydroCaseCoupling> coupling;
        if(separate) {
            DictionaryInput::Document electric("./system/electricity.dict");
            electric.root().allow({"Electricity","ConservativeColumn"});
            auto config=readDielectricConfig(electric.root().block("Electricity"));
            if(const auto* column=electric.root().find("ConservativeColumn")) {
                column->allow({"boundary","axisCapTraction"});
                column->get("boundary").choice({"groundedRadialWall"});
                const bool model=column->get("axisCapTraction").choice({"strictTwoPhase","regularHarmonic"})=="regularHarmonic";
                coupling=makeConservativeColumnCoupling(std::move(config),model);
            } else coupling=makeDielectricCoupling(std::move(config));
        } else coupling=makeDielectricCoupling(readDielectricConfig());
        int result = 0;
        try
        {
            HydroCase flow(coupling.get());
            if (cfg_.checkCase)
            {
                std::cout << "[Hydrodynamics] Case check passed\n";
                return 0;
            }
            result = flow.initialize();
            if (result != 0) return result;
            while (flow.hasNextStep())
            {
                result = flow.advanceNext();
                if (result != 0) return result;
            }
            flow.finish();
        }
        catch (const std::exception& error)
        {
            // Preserve the original flow-stage failure category/exit code.
            std::cerr << "[Hydrodynamics] ERROR: " << error.what() << '\n';
            return 2;
        }
        if (result == 0 && !cfg_.checkCase)
            std::cout << "[EHD] 电场/流动耦合时间循环已完成；逐步诊断：./output/series/，场数据：./output/fields/。\n"
                      << "[EHD] 短程计算完成不等于完整电喷雾物理验收。\n";
        return result;
    }
    catch (const std::exception& error)
    {
        std::cerr << "[EHD] ERROR: " << error.what() << '\n';
        return 2;
    }
}
