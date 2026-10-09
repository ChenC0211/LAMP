#include "Runtime/IO/DictionaryManualUnits.h"
#include "Core/IO/ConfigurationTrace.h"
#include "Core/Units/UnitsContext.h"
#include <cmath>

std::shared_ptr<const UnitsContext> readManualSIUnits(const std::string& path,
    DictionaryInput::ConfigurationTrace& trace)
{
    DictionaryInput::Document document(path);
    document.root().allow({"Units"});
    const auto& block = document.root().block("Units");
    block.allow({"mode", "L0", "T0", "Phi0", "rho0", "n0", "writeSIOutputs"});
    trace.record(block.get("mode"), block.get("mode").choice({"SI"}), "", "common UnitsContext and physical field output");
    const auto quantity = [&](const char* name, UnitParser::Dimension dimension, const char* unit) {
        const auto& entry = block.get(name);
        const double value = entry.quantity(dimension, unit);
        if (!(value > 0.0)) entry.fail("a positive finite SI reference scale");
        return trace.record(entry, value, unit, "common SI/code conversion reference");
    };
    UnitsScales scales;
    scales.L0_m = quantity("L0", UnitParser::Dimension::Length, "m");
    scales.T0_s = quantity("T0", UnitParser::Dimension::Time, "s");
    scales.Phi0_V = quantity("Phi0", UnitParser::Dimension::Voltage, "V");
    scales.rho0_Cpm3 = quantity("rho0", UnitParser::Dimension::VolumeChargeDensity, "C/m^3");
    scales.n0_m3 = quantity("n0", UnitParser::Dimension::NumberDensity, "1/m^3");
    const bool si = trace.record(block.get("writeSIOutputs"), block.get("writeSIOutputs").boolean(), "", "SI Tecplot contract");
    if (!si) block.get("writeSIOutputs").fail("true; this workflow publishes SI physical fields");
    auto result = std::make_shared<UnitsContext>(UnitsContext::MakeSI(scales));
    if (!std::isfinite(result->scales().v0_mps) || !std::isfinite(result->scales().E0_Vpm)
        || result->scales().v0_mps <= 0 || result->scales().E0_Vpm <= 0)
        block.fail("finite positive derived velocity and electric-field reference scales");
    return result;
}
