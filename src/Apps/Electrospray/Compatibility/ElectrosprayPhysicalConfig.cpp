#include "Electrospray/Compatibility/ElectrosprayPhysicalConfig.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace
{
std::string clean(std::string line)
{
    const auto comment = std::min(line.find('#'), line.find("//"));
    if (comment != std::string::npos) line.resize(comment);
    const auto first = line.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return {};
    return line.substr(first, line.find_last_not_of(" \t\r\n") - first + 1);
}

double number(const std::string& token)
{
    std::size_t consumed = 0;
    const double value = std::stod(token, &consumed);
    if (consumed != token.size() || !std::isfinite(value))
        throw std::invalid_argument("expected a finite SI number, without a unit suffix: " + token);
    return value;
}

int id(const std::string& token)
{
    std::size_t consumed = 0;
    const int value = std::stoi(token, &consumed);
    if (consumed != token.size() || value <= 0)
        throw std::invalid_argument("expected a positive 1-based ID: " + token);
    return value;
}
}

ElectrosprayPhysicalConfig ElectrosprayPhysicalConfig::read(const std::string& path)
{
    std::ifstream input(path);
    if (!input) throw std::runtime_error("[ElectrosprayPhysicalConfig] cannot open " + path);
    ElectrosprayPhysicalConfig config;
    double legacyTolerance = 0.0;
    const std::map<std::string, double*> scalars{
        {"liquidRelativePermittivity", &config.liquid.relativePermittivity},
        {"liquidConductivity_SI", &config.liquid.conductivitySI},
        {"liquidDensity_SI", &config.liquid.densitySI},
        {"liquidDynamicViscosity_SI", &config.liquid.dynamicViscositySI},
        {"gasRelativePermittivity", &config.gas.relativePermittivity},
        {"gasConductivity_SI", &config.gas.conductivitySI},
        {"gasDensity_SI", &config.gas.densitySI},
        {"gasDynamicViscosity_SI", &config.gas.dynamicViscositySI},
        {"conductorPotential_SI", &config.conductorPotentialSI},
        {"surfaceTension_SI", &config.surfaceTensionSI},
        {"initialSurfaceCharge_SI", &config.initialSurfaceChargeSI},
        {"linearResidualTolerance", &legacyTolerance}};
    const std::map<std::string, int*> ids{
        {"liquidPhase", &config.liquid.phaseId}, {"gasPhase", &config.gas.phaseId},
        {"conductorPhase", &config.conductorPhaseId}, {"movingBoundary", &config.movingBoundaryId}};
    const std::map<std::string, std::string*> strings{
        {"model", &config.model}, {"interfaceDiscretization", &config.interfaceDiscretization},
        {"fluidInterfaceLevelSet", &config.fluidInterfaceLevelSet},
        {"incomingMomentumModel", &config.incomingMomentumModel}, {"solidWallModel", &config.solidWallModel}};
    std::set<int> contactPointIds;
    std::set<std::string> seen;
    std::set<int> electricBoundaryIds;
    std::set<int> ambientBoundaryIds;
    int state = 0; // title, opening brace, entries, complete
    int lineNumber = 0;
    try
    {
        std::string line;
        while (std::getline(input, line))
        {
            ++lineNumber;
            line = clean(std::move(line));
            if (line.empty()) continue;
            if (state == 0 && line == "ElectrosprayPhysical") { state = 1; continue; }
            if (state == 0 && line == "ElectrosprayPhysical {") { state = 2; continue; }
            if (state == 1 && line == "{") { state = 2; continue; }
            if (state == 2 && line == "}") { state = 3; continue; }
            if (state != 2 || line.find_first_of("{}") != std::string::npos)
                throw std::invalid_argument("expected exactly one flat ElectrosprayPhysical block");
            std::istringstream stream(line);
            std::vector<std::string> tokens;
            for (std::string token; stream >> token;) tokens.push_back(std::move(token));
            const auto& key = tokens.front();
            const auto count = [&](std::size_t expected) {
                if (tokens.size() != expected) throw std::invalid_argument("wrong argument count for " + key);
            };
            if (key == "pinnedContactPoint")
            {
                count(2);
                const int point = id(tokens[1]);
                if (!contactPointIds.insert(point).second) throw std::invalid_argument("duplicate pinned contact point");
                config.pinnedContactPointIds.push_back(point);
                continue;
            }
            if (key == "domainBoundaryPotential_SI" || key == "domainBoundaryHomogeneousNeumann")
            {
                const bool dirichlet = key == "domainBoundaryPotential_SI";
                count(dirichlet ? 3 : 2);
                DomainElectricBoundary boundary;
                boundary.boundaryId = id(tokens[1]);
                boundary.type = dirichlet ? "dirichlet" : "homogeneousNeumann";
                boundary.potentialSI = dirichlet ? number(tokens[2]) : 0.0;
                if (!electricBoundaryIds.insert(boundary.boundaryId).second)
                    throw std::invalid_argument("duplicate electric domain boundary ID");
                config.domainElectricBoundaries.push_back(std::move(boundary));
                continue;
            }
            if (key == "ambientPressure_SI")
            {
                count(3);
                PressureBoundary boundary;
                boundary.boundaryId = id(tokens[1]);
                boundary.pressureSI = number(tokens[2]);
                if (!ambientBoundaryIds.insert(boundary.boundaryId).second)
                    throw std::invalid_argument("duplicate ambient pressure boundary ID");
                config.ambient.push_back(boundary);
                continue;
            }
            if (!seen.insert(key).second) throw std::invalid_argument("duplicate key " + key);
            if (key == "contactLineModel")
            {
                count(2);
                config.contactLineModel = tokens[1];
            }
            else if (key == "outletIncomingPhaseId")
            {
                count(2);
                config.outletIncomingPhaseId = id(tokens[1]);
            }
            else if (key == "supplyPressure_SI" || key == "outletPressure_SI")
            {
                count(3);
                auto& boundary = key == "supplyPressure_SI" ? config.supply : config.outlet;
                boundary.boundaryId = id(tokens[1]);
                boundary.pressureSI = number(tokens[2]);
            }
            else if (const auto found = scalars.find(key); found != scalars.end())
            {
                count(2);
                *found->second = number(tokens[1]);
            }
            else if (const auto found = ids.find(key); found != ids.end())
            {
                count(2);
                *found->second = id(tokens[1]);
            }
            else if (const auto found = strings.find(key); found != strings.end())
            {
                count(2);
                std::string value = tokens[1];
                if (value.size() >= 2 && value.front() == '"' && value.back() == '"')
                    value = value.substr(1, value.size() - 2);
                if (value.empty() || value.find('"') != std::string::npos)
                    throw std::invalid_argument("invalid string for " + key);
                *found->second = std::move(value);
            }
            else throw std::invalid_argument("unknown key " + key);
        }
        if (state != 3) throw std::invalid_argument("missing or unterminated ElectrosprayPhysical block");
        const auto required = [&](const std::string& key) {
            if (!seen.count(key)) throw std::invalid_argument("missing required key " + key);
        };
        for (const auto& entry : scalars) required(entry.first);
        for (const auto& entry : ids) required(entry.first);
        for (const auto& entry : strings) required(entry.first);
        required("supplyPressure_SI"); required("outletPressure_SI");
        if (config.model != "leakyDielectric") throw std::invalid_argument("model must be leakyDielectric");
        if (config.interfaceDiscretization != "cellHarmonic" && config.interfaceDiscretization != "sharpHermiteP3")
            throw std::invalid_argument("interfaceDiscretization must be cellHarmonic or sharpHermiteP3");
        if (config.incomingMomentumModel != "interiorContinuationZeroShear" || config.solidWallModel != "noSlip")
            throw std::invalid_argument("pressure inflow requires explicit interiorContinuationZeroShear and solidWallModel noSlip");
        // Old pressure-only snapshots need no contact-line law. A coupled
        // interface initialization below requires an explicit complete choice.
        if ((!config.contactLineModel.empty() && config.contactLineModel != "pinned")
            || (config.contactLineModel.empty() != config.pinnedContactPointIds.empty()))
            throw std::invalid_argument("contact line requires model pinned and explicit pinnedContactPoint IDs together");
        for (const auto* phase : {&config.liquid, &config.gas})
            if (!(phase->relativePermittivity > 0.0) || phase->conductivitySI < 0.0
                || !(phase->densitySI > 0.0) || !(phase->dynamicViscositySI > 0.0))
                throw std::invalid_argument("phase epsilon/density/viscosity must be positive and conductivity nonnegative");
        if (config.liquid.phaseId == config.gas.phaseId || config.liquid.phaseId == config.conductorPhaseId
            || config.gas.phaseId == config.conductorPhaseId)
            throw std::invalid_argument("liquid, gas and conductor phases must be distinct");
        if (config.outletIncomingPhaseId != -1 && config.outletIncomingPhaseId != config.liquid.phaseId
            && config.outletIncomingPhaseId != config.gas.phaseId)
            throw std::invalid_argument("outletIncomingPhaseId must name the configured liquid or gas phase");
        if (config.surfaceTensionSI < 0.0 || !(legacyTolerance > 0.0)
            || !(legacyTolerance < 1.0))
            throw std::invalid_argument("surface tension must be nonnegative and residual tolerance in (0,1)");
        if (config.supply.boundaryId == config.outlet.boundaryId
            || !electricBoundaryIds.count(config.supply.boundaryId)
            || !electricBoundaryIds.count(config.outlet.boundaryId))
            throw std::invalid_argument("distinct supply/outlet IDs require explicit electric domain conditions");
        if (config.domainElectricBoundaries.empty())
            throw std::invalid_argument("explicit electric domain boundary conditions are required");
        for (const auto& boundary : config.ambient)
            if (boundary.boundaryId == config.supply.boundaryId || boundary.boundaryId == config.outlet.boundaryId
                || !electricBoundaryIds.count(boundary.boundaryId))
                throw std::invalid_argument("ambient pressure boundary conflicts with supply/outlet or lacks electric data");
        // This is a model configuration, not a solver capability gate. The
        // current t=0 cellHarmonic App must separately reject nonzero free q.
    }
    catch (const std::exception& error)
    {
        throw std::runtime_error("[ElectrosprayPhysicalConfig] " + path + ":" + std::to_string(lineNumber)
            + ": " + error.what());
    }
    // Expand only here. Future versioned readers must supply each setting
    // independently, never infer conservation acceptance from solver accuracy.
    auto& numerical = config.numerics;
    numerical.flow = HydroFittedNumerics(legacyTolerance);
    numerical.electricFlux = DielectricFluxChecks(legacyTolerance);
    numerical.electricPotentialLinear = legacyTolerance;
    numerical.volumeCorrectionLinear = legacyTolerance;
    numerical.phaseFluxBalance = legacyTolerance;
    numerical.potentialBounds = legacyTolerance;
    numerical.polarizationConsistency = legacyTolerance;
    numerical.chargeClosure = 16.0 * legacyTolerance;
    config.legacyToleranceSource = path + "/ElectrosprayPhysical.linearResidualTolerance";
    return config;
}
