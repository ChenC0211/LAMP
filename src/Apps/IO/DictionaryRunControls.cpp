#include "IO/DictionaryRunControls.h"
#include "Core/IO/ConfigurationTrace.h"
#include "Core/Units/UnitsContext.h"
#include <cmath>

RunControls readRunControls(const DictionaryInput::Node& root,
    const UnitsContext& units, DictionaryInput::ConfigurationTrace& trace,bool allowStepRetry)
{
    using DictionaryInput::Node;
    RunControls result;
    const auto& time = root.block("Time");
    if(allowStepRetry) time.allow({"units", "dt", "endTime", "StepRetry"});
    else time.allow({"units", "dt", "endTime"});
    const auto timeUnits = trace.record(time.get("units"), time.get("units").choice({"code", "SI"}), "", "main clock/input conversion");
    const auto timeValue = [&](const Node& node) {
        const double value = timeUnits == "code" ? node.number()
            : units.TimeSIToCode(node.quantity(UnitParser::Dimension::Time, "s"));
        if (!(value > 0.0) || !std::isfinite(value)) node.fail("positive finite time after conversion to code units");
        return trace.record(node, value, "code time", "main clock/output schedule");
    };
    result.dt = timeValue(time.get("dt"));
    result.endTime = timeValue(time.get("endTime"));
    if(const auto* retry=time.find("StepRetry")) {
        retry->allow({"minDt","maxRetries","recoverAfter"});
        RunControls::StepRetry policy;
        policy.minDt=timeValue(retry->get("minDt"));
        policy.maxRetries=trace.record(retry->get("maxRetries"),retry->get("maxRetries").integer(),"attempts","whole-step halving limit");
        const auto* recovery=retry->find("recoverAfter");
        policy.recoverAfter=recovery
            ?trace.record(*recovery,recovery->integer(),"accepted steps","retry-free steps before a 1.25x dt proposal, capped by Time.dt")
            :trace.defaultValue(*retry,"recoverAfter",0,"accepted steps","disabled; retain the smaller accepted dt");
        if(policy.recoverAfter<0 || policy.recoverAfter>1000000)
            recovery->fail("an integer from 0 (disabled) to 1000000");
        if(policy.minDt>result.dt) retry->get("minDt").fail("minDt no greater than dt");
        if(policy.maxRetries<1 || policy.maxRetries>30) retry->get("maxRetries").fail("an integer from 1 to 30");
        if(!std::isfinite(result.endTime/policy.minDt) || result.endTime/policy.minDt>=std::numeric_limits<int>::max()-1.)
            retry->get("minDt").fail("fewer than INT_MAX accepted steps to endTime");
        result.stepRetry=policy;
    }
    if (!result.stepRetry && !result.fixedStepCount())
        time.get("endTime").fail("an integer number of fixed dt steps (1..INT_MAX)");
    const auto& output = root.block("Output");
    output.allow({"prefix", "Fields", "Series"});
    std::string seriesProfile="physical";
    if(const auto* series=output.find("Series")) {
        series->allow({"profile"});
        seriesProfile=series->find("profile")
            ? trace.record(series->get("profile"),series->get("profile").choice({"physical","diagnostic"}),"dimensionless","Runtime.output.seriesProfile")
            : trace.defaultValue(*series,"profile",std::string("physical"),"dimensionless","Runtime.output.seriesProfile");
    } else trace.defaultValue(output,"Series.profile",std::string("physical"),"dimensionless","Runtime.output.seriesProfile");
    result.diagnosticSeries=seriesProfile=="diagnostic";
    result.outputPrefix = trace.record(output.get("prefix"), output.get("prefix").word(), "", "field/history file prefix");
    if (result.outputPrefix.empty() || result.outputPrefix.find_first_of("/\\:") != std::string::npos
        || result.outputPrefix == "." || result.outputPrefix == "..") output.get("prefix").fail("a nonempty filename prefix, not a path");
    const auto& fields = output.block("Fields");
    fields.allow({"enabled", "format", "profile", "writeControl", "writeInterval", "writeInitial", "writeFirst", "writeFinal"});
    const auto* profile = fields.find("profile");
    const std::string fieldProfile = profile
        ? trace.record(*profile, profile->choice({"physical", "diagnostic"}), "", "physical fields or additional native/duplicate diagnostic dumps")
        : trace.defaultValue(fields, "profile", std::string("physical"), "", "physical fields without native/duplicate diagnostic dumps");
    result.diagnosticFields = fieldProfile == "diagnostic";
    result.writeFields = trace.record(fields.get("enabled"), fields.get("enabled").boolean(), "", "field writers only; histories remain enabled");
    trace.record(fields.get("format"), fields.get("format").choice({"tecplot"}), "", "Tecplot DAT field writers");
    const auto schedule = trace.record(fields.get("writeControl"), fields.get("writeControl").choice({"timeInterval", "stepInterval"}), "", "field schedule");
    if (schedule == "timeInterval") result.writeInterval = timeValue(fields.get("writeInterval"));
    else
    {
        result.writeEverySteps = trace.record(fields.get("writeInterval"), fields.get("writeInterval").integer(), "steps", "field schedule");
        if (result.writeEverySteps < 1) fields.get("writeInterval").fail("a positive integer step interval");
        result.writeInterval = result.dt; // positive adapter value, unused by step scheduling
    }
    result.writeInitial = trace.record(fields.get("writeInitial"), fields.get("writeInitial").boolean(), "", "initial field writers");
    result.writeFirst = trace.record(fields.get("writeFirst"), fields.get("writeFirst").boolean(), "", "force an accepted step-1 field frame");
    result.writeFinal = trace.record(fields.get("writeFinal"), fields.get("writeFinal").boolean(), "", "force an accepted final field frame");
    const auto& log = root.block("Log");
    log.allow({"everySteps", "level"});
    result.logEverySteps = trace.record(log.get("everySteps"), log.get("everySteps").integer(), "steps", "accepted-step console summary");
    if (result.logEverySteps < 1) log.get("everySteps").fail("a positive integer");
    result.detailedLog = trace.record(log.get("level"), log.get("level").choice({"summary", "detail"}), "", "detailed diagnostics/volume-correction progress") == "detail";
    if (const auto* restart = root.find("Restart"))
    {
        restart->allow({"model"});
        trace.record(restart->get("model"), restart->get("model").choice({"none"}), "", "fresh initialization; restart is not yet supported");
    }
    else trace.defaultValue(root, "Restart.model", std::string("none"), "", "fresh initialization");
    return result;
}
