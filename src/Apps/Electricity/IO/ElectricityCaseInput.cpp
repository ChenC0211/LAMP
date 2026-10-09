#include "Electricity/IO/ElectricityCaseInput.h"
#include "Core/IO/ConfigurationTrace.h"
#include "Core/Units/UnitsContext.h"
#include "Runtime/IO/DictionaryRunControls.h"
#include <filesystem>
#include <limits>
#include <utility>
#include <set>
#include <stdexcept>

ElectricityCaseInput readElectricityCaseInput()
{
    const std::set<std::string> allowed{"run.dict","electricity.dict","geometry.dict","mesh.dict"};
    for (const auto& entry : std::filesystem::directory_iterator("system"))
        if (!entry.is_regular_file() || !allowed.count(entry.path().filename().string()))
            throw std::invalid_argument("standalone Electricity supports only run/electricity/geometry/mesh.dict: "+entry.path().string());
    for (const auto& name : allowed)
        if (!std::filesystem::is_regular_file(std::filesystem::path("system")/name))
            throw std::invalid_argument("missing Electricity input: "+name);
    DictionaryInput::Document run("system/run.dict");
    run.root().allow({"Run","Time","Output","Log","Restart"});
    run.root().block("Run").allow({"app"});
    run.root().block("Run").get("app").choice({"Electricity"});
    run.root().block("Time").get("units").choice({"code"});
    const auto units=UnitsContext::MakeCodeUnits();
    DictionaryInput::ConfigurationTrace trace;
    auto controls=readRunControls(run.root(),units,trace);
    if (controls.fixedStepCount().value()>=std::numeric_limits<int>::max())
        run.root().block("Time").get("endTime").fail("at most INT_MAX-1 steps plus the initial field evaluation");
    if (controls.diagnosticFields)
        run.root().block("Output").block("Fields").get("profile").fail("physical; no additional diagnostic profile in this runner");

    const double dt=controls.dt;
    return {std::move(controls),readDielectricDomainInput(dt)};
}
