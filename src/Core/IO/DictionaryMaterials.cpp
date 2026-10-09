#include "Core/IO/DictionaryMaterials.h"
#include "Core/IO/ConfigurationTrace.h"

#include <algorithm>
#include <utility>

MaterialTable readMaterialTable(const std::string& path, DictionaryInput::ConfigurationTrace& trace)
{
    using D = UnitParser::Dimension;
    DictionaryInput::Document document(path);
    document.root().allow({"Materials"});
    const auto& materials = document.root().block("Materials");
    materials.allow({"Phases", "Interfaces"});
    MaterialTable result;
    const auto& phases = materials.block("Phases");
    if (phases.entries.empty()) phases.fail("at least one named phase material");
    for (const auto& block : phases.entries)
    {
        block.allow({"phaseId", "density", "dynamicViscosity", "relativePermittivity", "conductivity"});
        PhaseMaterial phase;
        const auto& id = block.get("phaseId");
        phase.phaseId = trace.record(id, id.integer(), "1-based ID", "geometry phase binding");
        if (phase.phaseId <= 0 || result.phases.count(phase.phaseId)) id.fail("a distinct positive phase ID");
        const auto quantity = [&](const char* key, D dimension, const char* unit, const char* consumer, bool positive) {
            const auto& entry = block.get(key);
            const double value = entry.quantity(dimension, unit);
            if (positive ? value <= 0.0 : value < 0.0) entry.fail(positive ? "a positive SI material value" : "a nonnegative SI material value");
            return trace.record(entry, value, unit, consumer);
        };
        phase.densitySI = quantity("density", D::MassDensity, "kg/m^3", "Hydro cell/vertex density and pressure normalization", true);
        phase.dynamicViscositySI = quantity("dynamicViscosity", D::DynamicViscosity, "Pa*s", "Hydro cell/vertex dynamic viscosity", true);
        phase.conductivitySI = quantity("conductivity", D::Conductivity, "S/m", "Electricity phase conduction current", false);
        const auto& epsilon = block.get("relativePermittivity");
        phase.relativePermittivity = trace.record(epsilon, epsilon.number(), "dimensionless", "Electricity permittivity and Maxwell stress");
        if (phase.relativePermittivity <= 0.0) epsilon.fail("positive relative permittivity");
        result.phases.emplace(phase.phaseId, phase);
    }
    const auto& interfaces = materials.block("Interfaces");
    for (const auto& block : interfaces.entries)
    {
        block.allow({"phases", "surfaceTension"});
        const auto& ids = block.get("phases");
        auto pair = trace.record(ids, ids.ids(), "1-based IDs", "unordered material interface pair");
        if (pair.size() != 2 || !result.phases.count(pair[0]) || !result.phases.count(pair[1]))
            ids.fail("exactly two defined, distinct material phases");
        if (pair[0] > pair[1]) std::swap(pair[0], pair[1]);
        const auto key = std::make_pair(pair[0], pair[1]);
        if (result.surfaceTensionSI.count(key)) ids.fail("a unique unordered material interface pair");
        const auto& sigma = block.get("surfaceTension");
        const double value = trace.record(sigma, sigma.quantity(D::SurfaceTension, "N/m"), "N/m", "Hydro capillary coefficient");
        if (value < 0.0) sigma.fail("nonnegative surface tension");
        result.surfaceTensionSI.emplace(key, value);
    }
    return result;
}
