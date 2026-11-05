#include "Pic/IO/PicInput.h"
#include "Core/IO/StrictDictionary.h"
#include "Core/IO/DictionaryBlockFile.h"
#include "Core/IO/DictionaryBlockText.h"
#include "Core/IO/DictionaryGeometry.h"
#include "Core/IO/DictionaryUnits.h"
#include "Modules/Plasma/IO/DictionaryNormalization.h"
#include "Modules/Plasma/IO/PlasmaCollisionLoader.h"
#include "Modules/Plasma/IO/PicCaseInputLoader.h"
#include "Solvers/IO/DictionarySolvers.h"
#include <algorithm>
#include <cctype>
#include <filesystem>
#include <iostream>
#include <set>
#include <utility>
#include <vector>

namespace
{
using NormalizationLoadResult=PicNormalizationInput;
std::string lowerName(std::string value)
{
    std::transform(value.begin(), value.end(), value.begin(),
        [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
    return value;
}

bool validatePicSystemFiles(const std::string& systemDir)
{
    const std::set<std::string> allowed = {
        "run.dict",
        "pic.dict",
        "geometry.dict",
        "mcc.dict",
        "mesh.dict",
        "solvers.dict",
        "units.dict",
        "normalization.dict",
        "gas.dict",
        "crossSections.dict",
        "reactions.dict",
        "collisions.dict",
    };

    namespace fs = std::filesystem;
    if (!fs::exists(systemDir) || !fs::is_directory(systemDir))
    {
        std::cerr << "[Pic] ERROR: missing system directory: " << systemDir << std::endl;
        return false;
    }

    std::vector<std::string> extras;
    for (const auto& e : fs::directory_iterator(systemDir))
    {
        if (!e.is_regular_file())
            continue;
        const std::string name = e.path().filename().string();
        if (allowed.find(name) == allowed.end())
            extras.push_back(name);
    }

    if (!extras.empty())
    {
        std::cerr << "[Pic] ERROR: system/ contains non-whitelisted file(s):";
        for (const auto& f : extras)
            std::cerr << " " << f;
        std::cerr << std::endl;
        std::cerr << "[Pic] Allowed files are:";
        for (const auto& f : allowed)
            std::cerr << " " << f;
        std::cerr << std::endl;
        return false;
    }

    return true;
}

bool validateSharedConfigSources(const DictionaryBlockText& parser,
                                 const std::string& picText)
{
    // A shared file and the embedded PIC block must not silently compete.
    for (const auto& name : {"Solvers", "Units", "Normalization"})
    {
        const std::string lower = lowerName(name);
        const std::string path = "./system/" + lower + ".dict";
        std::string block;
        if (std::filesystem::exists(path) &&
            parser.extractFirstBlockText(picText, {name, lower}, block, DictionaryBlockText::Scope::TopLevel))
        {
            std::cerr << "[Pic] ERROR: duplicate shared configuration: "
                      << path << " and pic.dict/" << name
                      << ". Keep exactly one source.\n";
            return false;
        }
    }
    return true;
}

std::shared_ptr<const Solvers::SolversConfig> loadSolversConfigOrDefault(const std::string& path)
{
    Solvers::DictionarySolvers dict;
    auto cfg = std::make_shared<Solvers::SolversConfig>();
    if (path.empty() || !dict.readDict(path))
        return cfg;

    *cfg = dict.config();
    return cfg;
}

struct UnitsLoadResult
{
    std::shared_ptr<const UnitsContext> ctx;
    bool writeSIOutputs = false;
};

UnitsLoadResult loadUnitsContextOrDefault(const std::string& path)
{
    DictionaryUnits dict;
    UnitsLoadResult r;
    if (!dict.readDict(path))
    {
        r.ctx = std::make_shared<UnitsContext>(UnitsContext::MakeCodeUnits());
        return r;
    }

    r.ctx = std::make_shared<UnitsContext>(dict.context());
    r.writeSIOutputs = dict.writeSIOutputs();
    return r;
}

NormalizationLoadResult loadNormalizationConfigOrDefault(const std::string& path)
{
    DictionaryNormalization dict;
    NormalizationLoadResult r;
    if (!dict.readDict(path))
    {
        // Default: 2D slab, thickness = 1 m.
        r.config = NormalizationConfig{};
        r.config.strategy = NormalizationStrategy::kManual;
        r.config.dimension = 2;
        r.config.slabThickness_m = 1.0;
        return r;
    }

    r.config = dict.config();
    r.hasDict = true;
    r.hasDimension = dict.hasDimension();
    r.hasSlabThickness = dict.hasSlabThickness();
    r.hasAny = dict.hasAny();
    return r;
}

bool extractPicBlockToCache(const DictionaryBlockText& parser,
                            const std::string& picText,
                            const std::vector<std::string>& names,
                            const std::string& label,
                            const std::string& blockName,
                            std::string& outPath,
                            bool& found)
{
    found = false;
    outPath.clear();
    std::string blockText;
    if (!parser.extractFirstBlockText(picText, names, blockText, DictionaryBlockText::Scope::TopLevel))
        return true;

    found = true;
    outPath = (std::filesystem::path("./cache") / ("pic_" + label + ".dict")).string();
    return writeDictionaryBlockFile(outPath, blockName, blockText);
}
} // namespace

PicAppConfig readPicAppConfig(const std::string& runDictPath)
{
    DictionaryInput::Document document(runDictPath);
    document.root().allow({"Run", "Pic", "Parallel"});
    const auto& run = document.root().block("Run");
    run.allow({"app"});
    run.get("app").choice({"Pic"});

    const auto* block = document.root().find("Pic");
    if (!block)
        return {};

    block->allow({"integrator"});
    const auto integrator = block->get("integrator").choice({"explicit", "implicit"});
    return {integrator == "explicit" ? PicIntegrator::Explicit : PicIntegrator::Implicit};
}

std::optional<PicSharedInput> readPicSharedInput()
{
    if (!validatePicSystemFiles("./system"))
        return std::nullopt;

    std::filesystem::create_directories("./output");
    std::filesystem::create_directories("./cache");

    PicSharedInput input;

    // pic.dict is the canonical shared PIC dictionary.
    const std::string primaryPicPath = std::filesystem::exists("./system/pic.dict")
        ? "./system/pic.dict" : std::string();
    DictionaryBlockText picParser;
    std::string picText;
    if (!primaryPicPath.empty() && !picParser.loadFile(primaryPicPath, picText))
    {
        std::cerr << "[Pic] ERROR: cannot read " << primaryPicPath << std::endl;
        return std::nullopt;
    }
    if (!validateSharedConfigSources(picParser, picText))
        return std::nullopt;

    const std::string picSource = primaryPicPath.empty() ? std::string("<missing pic dict>") : primaryPicPath;

    bool hasSolversBlock = false;
    std::string solversPath;
    if (!extractPicBlockToCache(picParser, picText, {"Solvers", "solvers"},
                                "solvers", "Solvers", solversPath, hasSolversBlock))
    {
        std::cerr << "[Pic] ERROR: failed to write Solvers block from " << picSource << std::endl;
        return std::nullopt;
    }
    input.solvers = loadSolversConfigOrDefault(hasSolversBlock ? solversPath : "./system/solvers.dict");

    bool hasUnitsBlock = false;
    std::string unitsPath;
    if (!extractPicBlockToCache(picParser, picText, {"Units", "units"},
                                "units", "Units", unitsPath, hasUnitsBlock))
    {
        std::cerr << "[Pic] ERROR: failed to write Units block from " << picSource << std::endl;
        return std::nullopt;
    }
    const auto unitsLoad = loadUnitsContextOrDefault(hasUnitsBlock ? unitsPath : "./system/units.dict");
    input.units = unitsLoad.ctx;
    input.writeSIOutputs = unitsLoad.writeSIOutputs;

    bool hasNormBlock = false;
    std::string normPath;
    if (!extractPicBlockToCache(picParser, picText, {"Normalization", "normalization"},
                                "normalization", "Normalization", normPath, hasNormBlock))
    {
        std::cerr << "[Pic] ERROR: failed to write Normalization block from " << picSource << std::endl;
        return std::nullopt;
    }
    auto normLoad = loadNormalizationConfigOrDefault(hasNormBlock ? normPath : "./system/normalization.dict");

    input.normalization = std::move(normLoad);
    input.normalizationSource = input.normalization.hasDict
        ? (hasNormBlock ? "pic.dict/Normalization" : "system/normalization.dict") : "default";
    input.primaryPicPath = primaryPicPath;
    return input;
}

PicGeometryClosure readPicGeometryClosure(const PicSharedInput& input,const UnitsContext& units)
{
    PicGeometryClosure closure;
    if (!input.normalization.hasDimension || !input.normalization.hasSlabThickness)
    {
        DictionaryGeometry geom;
        geom.readDict(&units);
        if (!input.normalization.hasDimension && geom.hasDimension)
        {
            closure.dimension = geom.dimension;
            closure.hasDimension = true;
        }
        if (!input.normalization.hasSlabThickness && geom.hasSlabThickness)
        {
            closure.slabThickness_m = geom.slabThickness_m;
            closure.hasSlabThickness = true;
        }
    }
    return closure;
}
std::shared_ptr<const PlasmaCollisionInputs> readPicRunCollisions(const PicSharedInput& input)
{
    return loadPlasmaCollisionInputs(input.primaryPicPath);
}
std::optional<PicCaseInput> readPicIntegratorInput(PicIntegrator integrator,
    const UnitsContext* units,const Solvers::SolversConfig* solvers)
{
    return readPicCaseInput(integrator,units,solvers);
}
