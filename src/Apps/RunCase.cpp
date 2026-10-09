#include "RunCase.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "Runtime/IO/DictionaryRun.h"
#include "Core/IO/DictionaryUnits.h"
#include "Modules/Plasma/IO/DictionaryNormalization.h"
#include "Runtime/AppRegistry.h"

#include "Core/IO/DictionaryGeometry.h"
#include "Core/IO/DictionaryBlockFile.h"
#include "Modules/Plasma/IO/PlasmaCollisionLoader.h"

#include "Solvers/IO/DictionarySolvers.h"
#include "Core/IO/DictionaryBlockText.h"

// Built-in apps
#include "Electrospray/ElectrosprayApp.h"
#include "Electrohydrodynamics/ElectrohydrodynamicsApp.h"
#include "SingleModule/SingleModuleApp.h"
#include "Interface/InterfaceApp.h"
#include "Electricity/ElectricityApp.h"
#include "Verification/Pic/PicVerificationApp.h"
#include "Verification/Mesh/MeshVerificationApp.h"

using namespace std;

namespace
{
std::string lowerName(std::string value)
{
    std::transform(value.begin(), value.end(), value.begin(),
        [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
    return value;
}

AppRegistry builtInApps()
{
    AppRegistry appRegistry;
    appRegistry.registerApp("MeshVerification", [](const RunConfig& c) { return make_unique<MeshVerificationApp>(c); });
    appRegistry.registerApp("Interface", [](const RunConfig& c) { return make_unique<InterfaceApp>(c); });
    appRegistry.registerApp("Electricity", [](const RunConfig& c) { return make_unique<ElectricityApp>(c); });
    // Built-in apps (case-insensitive keys handled by AppRegistry)
    appRegistry.registerApp("SingleModule", [](const RunConfig& c) { return make_unique<SingleModuleApp>(c); });
    appRegistry.registerApp("SingleModuleApp", [](const RunConfig& c) { return make_unique<SingleModuleApp>(c); });

    appRegistry.registerApp("Electrospray", [](const RunConfig& c) { return make_unique<ElectrosprayApp>(c); });
    appRegistry.registerApp("ElectrosprayApp", [](const RunConfig& c) { return make_unique<ElectrosprayApp>(c); });
    appRegistry.registerApp("Spray", [](const RunConfig& c) { return make_unique<ElectrosprayApp>(c); });

    appRegistry.registerApp("Electrohydrodynamics", [](const RunConfig& c) { return make_unique<ElectrohydrodynamicsApp>(c); });
    appRegistry.registerApp("ElectrohydrodynamicsApp", [](const RunConfig& c) { return make_unique<ElectrohydrodynamicsApp>(c); });
    appRegistry.registerApp("EHD", [](const RunConfig& c) { return make_unique<ElectrohydrodynamicsApp>(c); });

    appRegistry.registerApp("PicVerification", [](const RunConfig& c) { return make_unique<PicVerificationApp>(c); });
    appRegistry.registerApp("PicVerificationApp", [](const RunConfig& c) { return make_unique<PicVerificationApp>(c); });
    return appRegistry;
}

bool validateLegacySystemWhitelist(const std::string& systemDir, bool nodalFlow)
{
    // Legacy shared inputs. Formal Apps enforce their own file contracts.
    std::set<std::string> allowed = {
        "run.dict",
        // Canonical PIC dictionary.
        "pic.dict",
        "geometry.dict",
        "mcc.dict",
        "mesh.dict",
        // Reusable physics modules and the shared dictionaries loaded below.
        "hydro.dict",
        "electricity.dict",
        "solvers.dict",
        "units.dict",
        "normalization.dict",
        "gas.dict",
        "crossSections.dict",
        "reactions.dict",
        "collisions.dict",
    };
    if (nodalFlow) allowed.insert("interface.dict");

    namespace fs = std::filesystem;
    if (!fs::exists(systemDir) || !fs::is_directory(systemDir))
    {
        std::cerr << "[RunController] ERROR: missing system directory: " << systemDir << std::endl;
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
        std::cerr << "[RunController] ERROR: system/ contains non-whitelisted file(s):";
        for (const auto& f : extras)
            std::cerr << " " << f;
        std::cerr << std::endl;
        std::cerr << "[RunController] Allowed files are:";
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
    // A shared file and the compatibility PIC block must not silently compete.
    // Existing PIC-only cases keep their original embedded configuration path.
    for (const auto& name : {"Solvers", "Units", "Normalization"})
    {
        const std::string lower = lowerName(name);
        const std::string path = "./system/" + lower + ".dict";
        std::string block;
        if (std::filesystem::exists(path) &&
            parser.extractFirstBlockText(picText, {name, lower}, block, DictionaryBlockText::Scope::TopLevel))
        {
            std::cerr << "[RunController] ERROR: duplicate shared configuration: "
                      << path << " and pic.dict/" << name
                      << ". Keep exactly one source.\n";
            return false;
        }
    }
    return true;
}

std::vector<std::string> selectedSingleModules(const RunConfig& cfg)
{
    const auto app = lowerName(cfg.app);
    if (app != "singlemodule" && app != "singlemoduleapp")
        return {};
    return SingleModuleApp::readModulesFromRunDict("./system/run.dict");
}

bool isHydroName(const std::string& module)
{
    const auto name = lowerName(module);
    return name == "hydrodynamics" || name == "hydro";
}

bool validateStandaloneHydroSources(const DictionaryBlockText& parser,
                                    const std::string& picText)
{
    // This deliberately small entry point consumes raw code units and the
    // checked SparseLU options from hydro.dict. Do not silently ignore knobs
    // or emit SI/normalization reports for physics that does not consume them.
    for (const auto& name : {"Solvers", "Units", "Normalization"})
    {
        const std::string lower = lowerName(name);
        const std::string path = "./system/" + lower + ".dict";
        std::string block;
        if (std::filesystem::exists(path) ||
            parser.extractFirstBlockText(picText, {name, lower}, block, DictionaryBlockText::Scope::TopLevel))
        {
            std::cerr << "[RunController] ERROR: standalone Hydrodynamics does not support shared "
                      << name << " configuration (" << path << " or pic.dict/" << name
                      << "). Use numeric code units and solver options in hydro.dict.\n";
            return false;
        }
    }
    return true;
}

std::shared_ptr<const Solvers::SolversConfig> loadSolversConfigOrDefault(const std::string& path)
{
    Solvers::DictionarySolvers dict;
    if (path.empty() || !dict.readDict(path))
    {
        auto cfg = std::make_shared<Solvers::SolversConfig>();
        return cfg;
    }

    auto cfg = std::make_shared<Solvers::SolversConfig>();
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
    if (!dict.readDict(path))
    {
        UnitsLoadResult r;
        r.ctx = std::make_shared<UnitsContext>(UnitsContext::MakeCodeUnits());
        r.writeSIOutputs = false;
        return r;
    }

    UnitsLoadResult r;
    r.ctx = std::make_shared<UnitsContext>(dict.context());
    r.writeSIOutputs = dict.writeSIOutputs();
    return r;
}

void writeUnitsReportOrWarn(const std::string& outDir,
                            const UnitsContext& units,
                            bool writeSIOutputs)
{
    namespace fs = std::filesystem;
    std::error_code ec;
    fs::create_directories(outDir, ec);

    const std::string filePath = (fs::path(outDir) / "units_report.txt").string();
    std::ofstream out(filePath, std::ios::out);
    if (!out)
    {
        std::cerr << "[RunController] WARNING: failed to write units report: " << filePath << std::endl;
        return;
    }

    out << units.TableString();
    out << "\n";
    out << "OutputOptions\n";
    out << "  writeSIOutputs : " << (writeSIOutputs ? "true" : "false") << "\n";
    out << "\n";
    out << "Conversions\n";
    out << "  1 code length   = " << std::scientific << std::setprecision(8) << units.scales().L0_m << " m\n";
    out << "  1 code time     = " << std::scientific << std::setprecision(8) << units.scales().T0_s << " s\n";
    out << "  1 code potential= " << std::scientific << std::setprecision(8) << units.scales().Phi0_V << " V\n";
    out << "  1 code density  = " << std::scientific << std::setprecision(8) << units.scales().n0_m3 << " 1/m^3\n";
    out << "  1 code velocity = " << std::scientific << std::setprecision(8) << units.scales().v0_mps << " m/s\n";
    out << "  1 code field    = " << std::scientific << std::setprecision(8) << units.scales().E0_Vpm << " V/m\n";
    out << "  1 code rho      = " << std::scientific << std::setprecision(8) << units.scales().rho0_Cpm3 << " C/m^3\n";
}


struct NormalizationLoadResult
{
    NormalizationConfig config;
    bool hasDict = false;
    bool hasDimension = false;
    bool hasSlabThickness = false;
    bool hasAny = false;
};

NormalizationLoadResult loadNormalizationConfigOrDefault(const std::string& path)
{
    DictionaryNormalization dict;
    if (!dict.readDict(path))
    {
        NormalizationLoadResult r;
        // Default: 2D slab, thickness = 1 m.
        r.config = NormalizationConfig{};
        r.config.strategy = NormalizationStrategy::kManual;
        r.config.dimension = 2;
        r.config.slabThickness_m = 1.0;
        r.hasDict = false;
        r.hasDimension = false;
        r.hasSlabThickness = false;
        r.hasAny = false;
        return r;
    }

    NormalizationLoadResult r;
    r.config = dict.config();
    r.hasDict = true;
    r.hasDimension = dict.hasDimension();
    r.hasSlabThickness = dict.hasSlabThickness();
    r.hasAny = dict.hasAny();
    return r;
}

void writeNormalizationReportOrWarn(const std::string& outDir,
                                   const NormalizationContext& norm,
                                   bool explicitConfig,
                                   const std::string& source,
                                   const std::string& dimensionSource,
                                   const std::string& slabSource)
{
    namespace fs = std::filesystem;
    std::error_code ec;
    fs::create_directories(outDir, ec);

    const std::string filePath = (fs::path(outDir) / "normalization_report.txt").string();
    std::ofstream out(filePath, std::ios::out);
    if (!out)
    {
        std::cerr << "[RunController] WARNING: failed to write normalization report: " << filePath << std::endl;
        return;
    }

    out << norm.TableString();
    out << "\n";
    out << "NormalizationConfig\n";
    out << "  explicit : " << (explicitConfig ? "true" : "false") << "\n";
    out << "  source   : " << (source.empty() ? "(unknown)" : source) << "\n";
    out << "  dimension_source : " << (dimensionSource.empty() ? "(unknown)" : dimensionSource) << "\n";
    out << "  slabThickness_source : " << (slabSource.empty() ? "(unknown)" : slabSource) << "\n";
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
    if (!writeDictionaryBlockFile(outPath, blockName, blockText))
        return false;
    return true;
}

} // namespace

int runCase(bool checkCase, RunStartMode startMode)
{
    namespace fs = std::filesystem;

    DictionaryRun dictRun;
    RunConfig cfg;
    const bool hasRunDict = dictRun.readDict("./system/run.dict");
    if (!hasRunDict && fs::exists("./system/run.dict"))
    {
        std::cerr << "[RunController] ERROR: failed to read system/run.dict" << std::endl;
        return 2;
    }
    if (hasRunDict)
    {
        cfg = dictRun.config();
        if (cfg.app.empty())
        {
            std::cerr << "[RunController] ERROR: missing Run.app in system/run.dict" << std::endl;
            return 2;
        }
    }
    else
    {
        cfg.app = "SingleModule";
    }

    cfg.checkCase = checkCase;
    cfg.startMode = startMode;
    if (startMode != RunStartMode::Automatic && (cfg.schemaVersion != 2 || cfg.app != "Electrospray"))
    {
        std::cerr << "[LAMP] --start 当前仅支持新版 Electrospray 算例，其他入口保持原启动方式。\n";
        return 2;
    }
    const bool interfaceOnly = cfg.app == "Interface";
    const auto selectedModules = selectedSingleModules(cfg);
    const bool hydroOnly = selectedModules.size() == 1 && isHydroName(selectedModules.front());
    const auto appName = lowerName(cfg.app);
    const bool ehd = appName == "electrohydrodynamics" || appName == "electrohydrodynamicsapp" || appName == "ehd";
    if (interfaceOnly || cfg.app == "Electricity" || cfg.app == "MeshVerification")
    {
        // Prescribed inputs own their readers; no Plasma caches or reports.
        // Existing App schemas and startup policies remain unchanged.
        return builtInApps().run(cfg);
    }
    if (cfg.schemaVersion == 2)
    {
        // The strict App loader owns its module/material binding. Do not load
        // Plasma dictionaries, extract caches or write reports on this path.
        if (cfg.app != "Electrospray") return 2;
        return builtInApps().run(cfg);
    }
    if (!validateLegacySystemWhitelist("./system", hydroOnly || ehd)) return 2;
    if (cfg.checkCase && !hydroOnly && !ehd)
    {
        std::cerr << "[RunController] ERROR: --check-case currently supports only "
                     "SingleModule { modules ( Hydrodynamics ) } or Electrohydrodynamics; other modules "
                     "have no preflight contract yet.\n";
        return 2;
    }
    if (!hydroOnly && std::any_of(selectedModules.begin(), selectedModules.end(), isHydroName))
    {
        std::cerr << "[RunController] ERROR: standalone Hydrodynamics must run alone "
                     "in SingleModule.modules; multi-module coupling is not implemented.\n";
        return 2;
    }

    auto choosePrimaryPicDictPath = [&]() -> std::string {
        // pic.dict is the canonical shared PIC dictionary.
        namespace fs2 = std::filesystem;
        if (fs2::exists("./system/pic.dict"))
            return "./system/pic.dict";
        return std::string();
    };

    const std::string primaryPicPath = choosePrimaryPicDictPath();
    DictionaryBlockText picParser;
    std::string picText;
    if (!primaryPicPath.empty() && !picParser.loadFile(primaryPicPath, picText))
    {
        std::cerr << "[RunController] ERROR: cannot read " << primaryPicPath << std::endl;
        return 2;
    }
    if (!validateSharedConfigSources(picParser, picText))
        return 2;

    if (hydroOnly || ehd)
    {
        if (!validateStandaloneHydroSources(picParser, picText))
            return 2;
        cfg.units = std::make_shared<UnitsContext>(UnitsContext::MakeCodeUnits());
        cfg.normalization = std::make_shared<NormalizationContext>(
            NormalizationContext::MakeFromUnits(cfg.units, NormalizationConfig{}));
        cfg.solvers = std::make_shared<Solvers::SolversConfig>();
        // Strict Hydro preflight owns Domain/mesh validation. In particular,
        // --check-case must not enter the legacy geometry parser, extract PIC
        // blocks into cache, construct a mesh, or write runtime reports.
        return builtInApps().run(cfg);
    }

    bool hasSolversBlock = false;
    std::string solversPath;
    if (!extractPicBlockToCache(picParser, picText,
                                {"Solvers", "solvers"},
                                "solvers", "Solvers",
                                solversPath, hasSolversBlock))
    {
        std::cerr << "[RunController] ERROR: failed to write Solvers block from "
                  << (primaryPicPath.empty() ? std::string("<missing pic dict>") : primaryPicPath)
                  << std::endl;
        return 2;
    }

    if (hasSolversBlock)
    {
        cfg.solvers = loadSolversConfigOrDefault(solversPath);
    }
    else
    {
        cfg.solvers = loadSolversConfigOrDefault("./system/solvers.dict");
    }

    bool hasUnitsBlock = false;
    std::string unitsPath;
    if (!extractPicBlockToCache(picParser, picText,
                                {"Units", "units"},
                                "units", "Units",
                                unitsPath, hasUnitsBlock))
    {
        std::cerr << "[RunController] ERROR: failed to write Units block from "
                  << (primaryPicPath.empty() ? std::string("<missing pic dict>") : primaryPicPath)
                  << std::endl;
        return 2;
    }

    const std::string unitsPathToUse = hasUnitsBlock ? unitsPath : "./system/units.dict";
    const auto unitsLoad = loadUnitsContextOrDefault(unitsPathToUse);
    cfg.units = unitsLoad.ctx ? unitsLoad.ctx : std::make_shared<UnitsContext>(UnitsContext::MakeCodeUnits());
    cfg.writeSIOutputs = unitsLoad.writeSIOutputs;

    bool hasNormBlock = false;
    std::string normPath;
    if (!extractPicBlockToCache(picParser, picText,
                                {"Normalization", "normalization"},
                                "normalization", "Normalization",
                                normPath, hasNormBlock))
    {
        std::cerr << "[RunController] ERROR: failed to write Normalization block from "
                  << (primaryPicPath.empty() ? std::string("<missing pic dict>") : primaryPicPath)
                  << std::endl;
        return 2;
    }
    const std::string normPathToUse = hasNormBlock ? normPath : "./system/normalization.dict";
    auto normLoad = loadNormalizationConfigOrDefault(normPathToUse);

    // Track where normalization settings came from (for user-facing reports).
    bool normExplicit = normLoad.hasDict;
    const std::string normSourceConfig = normLoad.hasDict
                                             ? (hasNormBlock ? "pic.dict/Normalization" : "system/normalization.dict")
                                             : "default";
    std::string normSource = normSourceConfig;
    std::string normDimSource = normLoad.hasDimension ? normSourceConfig : "";
    std::string normSlabSource = normLoad.hasSlabThickness ? normSourceConfig : "";

    std::string normWarn;
    const UnitsScales effectiveScales = resolvePlasmaReferenceScales(*cfg.units, normLoad.config, normWarn);
    if (!normWarn.empty())
        std::cerr << "[RunController] WARNING: " << normWarn << std::endl;

    if (cfg.units->mode() == UnitsMode::kSI)
        cfg.units = std::make_shared<UnitsContext>(UnitsContext::MakeSI(effectiveScales));
    else
        cfg.units = std::make_shared<UnitsContext>(UnitsContext::MakeCodeUnits());

    // Allow Geometry.Domain to supply dimensional closure metadata (dimension / slabThickness)
    // when not explicitly set in Normalization.
    bool geometryProvidedClosure = false;
    if (!normLoad.hasDimension || !normLoad.hasSlabThickness)
    {
        DictionaryGeometry geom;
        geom.readDict(cfg.units.get());
        if (!normLoad.hasDimension && geom.hasDimension)
        {
            normLoad.config.dimension = geom.dimension;
            normDimSource = "geometry.dict/Domain";
            geometryProvidedClosure = true;
        }
        if (!normLoad.hasSlabThickness && geom.hasSlabThickness)
        {
            normLoad.config.slabThickness_m = geom.slabThickness_m;
            normSlabSource = "geometry.dict/Domain";
            geometryProvidedClosure = true;
        }
    }
    if (geometryProvidedClosure)
    {
        normExplicit = true;
        if (!normLoad.hasDict)
            normSource = "geometry.dict/Domain";
    }

    if (normDimSource.empty())
        normDimSource = "default";
    if (normLoad.config.dimension == 2)
    {
        if (normSlabSource.empty())
            normSlabSource = "default";
    }
    else
    {
        if (normSlabSource.empty())
            normSlabSource = "n/a";
    }

    cfg.normalization = std::make_shared<NormalizationContext>(NormalizationContext::MakeFromUnits(cfg.units, normLoad.config));

    std::cout << cfg.units->TableString();
    writeUnitsReportOrWarn("./output/meta", *cfg.units, cfg.writeSIOutputs);

    std::cout << cfg.normalization->TableString();
    {
        auto fmtSci = [](double v, int p = 6) {
            std::ostringstream oss;
            oss << std::scientific << std::setprecision(p) << v;
            return oss.str();
        };
        std::cout << "[RunController] Normalization source: " << normSource << "\n";
        std::cout << "[RunController] Normalization closure: dimension=" << normLoad.config.dimension
                  << " (" << normDimSource << ")";
        if (normLoad.config.dimension == 2)
        {
            std::cout << " slabThickness=" << fmtSci(normLoad.config.slabThickness_m) << " m"
                      << " (" << normSlabSource << ")";
        }
        std::cout << "\n";
    }
    if (!normExplicit)
    {
        std::cerr << "[RunController] WARNING: Normalization not specified; using defaults (dimension=2, slabThickness=1 m).\n";
    }
    if (normLoad.config.dimension == 2 && normSlabSource == "default")
    {
        std::cerr << "[RunController] WARNING: slabThickness not specified for 2D; defaulting to 1 m. "
                     "Set Normalization.slabThickness or Geometry.Domain.slabThickness.\n";
    }
    writeNormalizationReportOrWarn("./output/meta", *cfg.normalization, normExplicit, normSource, normDimSource, normSlabSource);

    cfg.plasmaCollisions = loadPlasmaCollisionInputs(primaryPicPath);
    if (!cfg.plasmaCollisions) return 2;

    const int rc = builtInApps().run(cfg);
    return rc;
}
