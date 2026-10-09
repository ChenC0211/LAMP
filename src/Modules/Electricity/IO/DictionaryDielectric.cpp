#include "Modules/Electricity/IO/DictionaryDielectric.h"
#include "Core/IO/StrictDictionary.h"
#include "Core/Units/UnitsContext.h"
#include "Modules/Electricity/DielectricEvolution.h"
#include <algorithm>
#include <iostream>
#include <cctype>
#include <cmath>
#include <fstream>
#include <set>
#include <stdexcept>

namespace {
class DictionaryReader
{
public:
    explicit DictionaryReader(const std::string& path) : path_(path)
    {
        std::ifstream input(path);
        if (!input) fail("cannot open dictionary");
        std::string line;
        int number = 0;
        while (std::getline(input, line))
        {
            ++number;
            for (std::size_t i = 0; i < line.size();)
            {
                const char c = line[i];
                if (std::isspace(static_cast<unsigned char>(c)) || c == ';') { ++i; continue; }
                if (c == '#' || (c == '/' && i + 1 < line.size() && line[i + 1] == '/')) break;
                if (c == '{' || c == '}')
                {
                    tokens_.push_back({std::string(1, c), number});
                    ++i;
                    continue;
                }
                const auto start = i;
                while (i < line.size() && !std::isspace(static_cast<unsigned char>(line[i])) &&
                       line[i] != ';' && line[i] != '{' && line[i] != '}' && line[i] != '#' &&
                       !(line[i] == '/' && i + 1 < line.size() && line[i + 1] == '/')) ++i;
                tokens_.push_back({line.substr(start, i - start), number});
            }
        }
    }
    // The standalone case uses a strict document with separate binding data.
    // Feed its physical block to the SAME parser as the historical EHD file.
    explicit DictionaryReader(const DictionaryInput::Node& block) : path_(block.file)
    {
        if (!block.isBlock || block.name != "Electricity") block.fail("an Electricity block");
        tokens_.push_back({"Electricity",block.line});
        tokens_.push_back({"{",block.line});
        const auto append = [&](const auto& self, const DictionaryInput::Node& entry) -> void {
            tokens_.push_back({entry.name,entry.line});
            if (entry.isBlock) {
                tokens_.push_back({"{",entry.line});
                for (const auto& child : entry.entries) self(self,child);
                tokens_.push_back({"}",entry.line});
            } else tokens_.push_back({entry.word(),entry.line});
        };
        for (const auto& entry : block.entries) append(append,entry);
        tokens_.push_back({"}",block.line});
    }
    bool done() const { return position_ == tokens_.size(); }
    std::string peek() const { return done() ? std::string() : tokens_[position_].text; }
    std::string next()
    {
        if (done()) fail("unexpected end of dictionary");
        return tokens_[position_++].text;
    }
    void expect(const std::string& value)
    {
        if (next() != value) fail("expected '" + value + "'");
    }
    double number(const std::string& key)
    {
        const auto token = next();
        std::size_t used = 0;
        double value = 0.0;
        try { value = std::stod(token, &used); }
        catch (...) { fail(key + " must be a finite number in code units"); }
        if (used != token.size() || !std::isfinite(value)) fail(key + " must be a finite number in code units");
        return value;
    }
    [[noreturn]] void fail(const std::string& message) const
    {
        const int line = done() ? (tokens_.empty() ? 1 : tokens_.back().line) : tokens_[position_].line;
        throw std::runtime_error(path_ + ":" + std::to_string(line) + ": " + message);
    }
private:
    struct Token { std::string text; int line; };
    std::string path_;
    std::vector<Token> tokens_;
    std::size_t position_ = 0;
};

DielectricConfig parseDielectricConfig(DictionaryReader& input, bool legacy)
{
    DielectricConfig config;
    input.expect("Electricity"); input.expect("{");
    std::set<std::string> seen;
    std::set<std::string> conservation;
    while (input.peek() != "}")
    {
        const auto key = input.next();
        if (!seen.insert(key).second) input.fail("duplicate Electricity key '" + key + "'");
        if (key == "model") config.model = input.next();
        else if (key == "interfaceDiscretization") config.interfaceDiscretization = input.next();
        else if (key == "permittivityNegative") config.permittivityNegative = input.number(key);
        else if (key == "permittivityPositive") config.permittivityPositive = input.number(key);
        else if (key == "electrodeAxis") config.electrodeAxis = input.next();
        else if (key == "lowerPotential") config.lowerPotential = input.number(key);
        else if (key == "upperPotential") config.upperPotential = input.number(key);
        else if (key == "linearResidualTolerance") config.linearResidualTolerance = input.number(key);
        else if (key == "conductivityNegative") config.conductivityNegative = input.number(key);
        else if (key == "conductivityPositive") config.conductivityPositive = input.number(key);
        else if (key == "initialSurfaceChargeDensity") config.initialSurfaceChargeDensity = input.number(key);
        else if (key == "chargeConvection")
        {
            const auto value = input.next();
            if (value != "true" && value != "false") input.fail("chargeConvection must be true or false");
            config.chargeConvection = value == "true";
        }
        else if (key == "Conservation")
        {
            if (legacy) input.fail("Conservation requires formal Electricity or separate EHD inputs");
            input.expect("{");
            while (input.peek() != "}")
            {
                const auto setting = input.next();
                if (!conservation.insert(setting).second) input.fail("duplicate Conservation key '" + setting + "'");
                if (setting == "boundaryFlux") config.boundaryFluxTolerance = input.number(setting);
                else if (setting == "chargeClosure") config.chargeClosureTolerance = input.number(setting);
                else if (setting == "reconstructedPoisson" || setting == "gaussClosure")
                {
                    if (!config.sharpFluxChecks) config.sharpFluxChecks.emplace();
                    const double value = input.number(setting);
                    if (setting == "reconstructedPoisson") config.sharpFluxChecks->reconstructedPoisson = value;
                    else config.sharpFluxChecks->gaussClosure = value;
                }
                else input.fail("unknown Conservation key '" + setting + "'");
            }
            input.expect("}");
        }
        else input.fail("unknown or unsupported Electricity key '" + key + "'");
    }
    input.expect("}");
    if (!input.done()) input.fail("unexpected content after Electricity block");
    for (const auto* key : {"model", "permittivityNegative", "permittivityPositive", "electrodeAxis", "lowerPotential", "upperPotential"})
        if (!seen.count(key)) input.fail(std::string("missing required Electricity key '") + key + "'");
    for (const auto* key : {"conductivityNegative", "conductivityPositive", "initialSurfaceChargeDensity", "chargeConvection"})
    {
        if (config.model == "leakyDielectric" && !seen.count(key)) input.fail(std::string("missing leakyDielectric key '") + key + "'");
        if (config.model != "leakyDielectric" && seen.count(key)) input.fail(std::string(key) + " is only valid for leakyDielectric");
    }
    if (legacy)
    {
        if (config.interfaceDiscretization == "diffuseTanh") config.boundaryFluxTolerance = config.linearResidualTolerance;
        else if (config.interfaceDiscretization == "sharpHermiteP3") config.sharpFluxChecks = DielectricFluxChecks(config.linearResidualTolerance);
        if (config.model == "leakyDielectric") config.chargeClosureTolerance = 16.0 * config.linearResidualTolerance;
    }
    try { validateDielectricConfig(config); }
    catch (const std::exception& error) { input.fail(error.what()); }
    return config;
}

}
DielectricConfig readDielectricConfig(const std::string& path)
{
    DictionaryReader input(path);
    return parseDielectricConfig(input, true);
}
DielectricConfig readDielectricConfig(const DictionaryInput::Node& block)
{
    DictionaryReader input(block);
    return parseDielectricConfig(input, false);
}

DielectricConfig readModularDielectricConfig(const std::string& path)
{
    DictionaryInput::Document document(path);
    document.root().allow({"Electricity"});
    return readDielectricConfig(document.root().block("Electricity"));
}

DielectricDomainInput readDielectricDomainInput(double dt, const std::string& electricPath,
    const std::string& geometryPath, const std::string& meshPath)
{
    DielectricDomainInput input;
    const auto units = UnitsContext::MakeCodeUnits();
    DictionaryInput::Document electric(electricPath);
    electric.root().allow({"Electricity","Domain"});
    input.electric=readDielectricConfig(electric.root().block("Electricity"));
    const auto& domain=electric.root().block("Domain");
    domain.allow({"coordinateSystem","interfaceLevelSet","interfaceHalfWidth","motion"});
    domain.get("motion").choice({"fixed"});
    input.electric.axisymmetric=domain.get("coordinateSystem").choice({"cartesianXY","axisymmetricRZ"})=="axisymmetricRZ";
    if(input.electric.interfaceDiscretization=="diffuseTanh") {
        input.diffuseThickness=domain.get("interfaceHalfWidth").number();
        if(!(input.diffuseThickness>0.0) || !std::isfinite(2*input.diffuseThickness))
            domain.get("interfaceHalfWidth").fail("a positive finite diffuse material half-width");
    } else if(const auto* legacyWidth=domain.find("interfaceHalfWidth")) {
        const double value=legacyWidth->number();
        if(!(value>0.0) || !std::isfinite(2*value))
            legacyWidth->fail("a positive finite legacy half-width (deprecated for sharpHermiteP3)");
        input.ignoredSharpHalfWidth=value;
        std::cout<<"[Electricity] Domain.interfaceHalfWidth is deprecated and ignored for sharpHermiteP3; remove it.\n";
    }
    validateDielectricConfig(input.electric);
    if (input.electric.chargeConvection)
        throw std::invalid_argument("fixed-interface Electricity has no velocity provider; chargeConvection must be false");
    if (input.electric.model=="leakyDielectric")
        Quantum::Electricity::DielectricEvolution::explicitRelaxationCfl(input.electric,dt);
    DictionaryInput::Document mesh(meshPath);
    mesh.root().allow({"adaptiveMesh","BackgroundMesh"});
    mesh.root().get("adaptiveMesh").choice({"false"});
    mesh.root().block("BackgroundMesh").allow({"divide"});
    // Geometry and mesh retain their public legacy readers. No discarded
    // auxiliary physics, refinement, periodic or moving-solid settings.
    input.geometry.readDict(geometryPath, &units);
    input.mesh.readDict(meshPath, geometryPath, &units);
    const auto& g=input.geometry;
    if (!g.hasDomain || !(g.domainMax.x>g.domainMin.x) || !(g.domainMax.y>g.domainMin.y)
        || !std::isfinite(g.domainMin.x) || !std::isfinite(g.domainMin.y)
        || !std::isfinite(g.domainMax.x) || !std::isfinite(g.domainMax.y)
        || !g.periodicPairs.empty() || g.extensionEnabled("MovingEB") || g.extensionEnabled("MotionSchedule")
        || !g.phaseName.empty() || !g.boundaries.empty() || g.hasSlabThickness
        || (g.hasDimension && g.dimension!=2) || g.levelSetName.size()!=1
        || g.levelSetName.front()!=domain.get("interfaceLevelSet").word()
        || std::any_of(g.moving.begin(),g.moving.end(),[](bool v){return v;}))
        throw std::invalid_argument("Electricity requires one named fixed scalar in a full nonperiodic 2D rectangle; no solid/phase/Boundary definitions");
    if (input.mesh.Nx<8 || input.mesh.Ny<8 || !input.mesh.refineBoxes.empty()
        || input.mesh.geometryRefinement.enabled || input.mesh.physicsRefinement.enabled
        || input.mesh.runtimeAdaptation.enabled)
        throw std::invalid_argument("Electricity requires a fixed uniform grid with Nx/Ny >= 8");
    if (input.electric.axisymmetric && (g.domainMin.x<0
        || (g.domainMin.x==0 && input.electric.electrodeAxis!="y")))
        throw std::invalid_argument("RZ requires nonnegative radius and an unelectrified symmetry axis");
    return input;
}
