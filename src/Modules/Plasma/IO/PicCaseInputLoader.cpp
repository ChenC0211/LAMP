#include "Modules/Plasma/IO/PicCaseInputLoader.h"
#include "Modules/Plasma/Particles/PicSpeciesState.h"
#include <stdexcept>
#include "Modules/Plasma/ExplicitPIC/IO/DictionaryPIC.h"
#include "Modules/Plasma/IO/DictionaryMovingEB.h"
#include "Modules/Plasma/IO/DictionaryMotionSchedule.h"
#include "Solvers/IO/DictionarySolvers.h"
#include "Solvers/SolversConfig.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>

namespace {
static std::string toLowerCopy(std::string s)
{
    for (auto& c : s)
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}
}

ExplicitPicAmrSetup::ExplicitPicAmrSetup(const PicCaseInput& input)
    : config_(input.pic.amr), refinerMesh_(input.mesh)
{
    const auto& dictMesh = input.mesh;
    auto& amrCfg = config_;
    auto& meshRefinerDict = refinerMesh_;
    const bool meshBlocksPresent = dictMesh.geometryRefinement.present || dictMesh.physicsRefinement.present;
    const bool meshPhysEnabled = meshBlocksPresent ? dictMesh.physicsRefinement.enabled : false;
    meshPhysicsEnabled_ = meshPhysEnabled;
    // If mesh.dict declares PhysicsRefinement, map its policy knobs onto PIC AMR config.
    if (dictMesh.physicsRefinement.present)
    {
        if (dictMesh.physicsRefinement.adaptEverySteps > 0)
            amrCfg.adaptEverySteps = dictMesh.physicsRefinement.adaptEverySteps;
        if (dictMesh.physicsRefinement.indicatorEverySteps > 0)
            amrCfg.indicatorEverySteps = dictMesh.physicsRefinement.indicatorEverySteps;
        amrCfg.report.every = std::max(0, dictMesh.physicsRefinement.reportEvery);
        amrCfg.report.enabled = (amrCfg.report.every > 0);
        if (dictMesh.physicsRefinement.persistenceSteps > 0)
            amrCfg.persistenceSteps = dictMesh.physicsRefinement.persistenceSteps;
        if (dictMesh.physicsRefinement.hysteresis >= 0.0)
            amrCfg.hysteresis = dictMesh.physicsRefinement.hysteresis;

        if (!dictMesh.physicsRefinement.indicators.empty())
        {
            bool useDebye = false;
            bool useEField = false;
            for (const auto& name : dictMesh.physicsRefinement.indicators)
            {
                const std::string n = toLowerCopy(name);
                if (n == "debyelength")
                    useDebye = true;
                else if (n == "electricfield")
                    useEField = true;
            }
            amrCfg.debye.enabled = amrCfg.debye.enabled && useDebye;
            amrCfg.indicator.useEField = amrCfg.indicator.useEField && useEField;
        }

        if (!meshPhysEnabled)
        {
            amrCfg.debye.enabled = false;
            amrCfg.indicator.useEField = false;
        }
    }

    if (!meshBlocksPresent)
    {
        meshRefinerDict.physicsRefinement.enabled = amrCfg.enabled;
        meshRefinerDict.physicsRefinement.indicators.clear();
        if (amrCfg.debye.enabled)
            meshRefinerDict.physicsRefinement.indicators.push_back("DebyeLength");
        if (amrCfg.indicator.useEField)
            meshRefinerDict.physicsRefinement.indicators.push_back("ElectricField");
        meshRefinerDict.physicsRefinement.hysteresis = amrCfg.hysteresis;
    }

}

void ExplicitPicAmrSetup::activateForMesh(int maxLevel, bool geometryRefinementActive, std::ostream& errors)
{
    if (meshActivated_) throw std::logic_error("Explicit PIC AMR mesh activation may only occur once");
    meshActivated_ = true;
    geometryRefinementActive_ = geometryRefinementActive;
    auto& amrEnabled = enabled_;
    amrEnabled = meshPhysicsEnabled_ || geometryRefinementActive_;
    if (amrEnabled && maxLevel <= 0)
    {
        errors << "[ExplicitPIC][AMR] WARNING: AMR enabled but mesh.maxLevel=0. "
                  << "Check mesh.dict adaptive settings (maxLevel_local). Disabling AMR.\n";
        amrEnabled = false;
    }

}

void ExplicitPicAmrSetup::selectElectron(const PicSpeciesState& particles, std::ostream& errors)
{
    if (!meshActivated_ || electronSelected_)
        throw std::logic_error("Explicit PIC AMR electron selection requires its mesh and may only occur once");
    electronSelected_ = true;
    auto& amrCfg = config_;
    auto& amrEnabled = enabled_;
    auto& amrElectronIndex = electronIndex_;
    const bool geometryRefinementActive = geometryRefinementActive_;
    // AMR: choose an electron species index (for Debye length estimation).
    if (amrEnabled && amrCfg.debye.enabled)
    {
        amrElectronIndex = particles.electronIndex();

        if (amrElectronIndex < 0)
        {
            errors << "[ExplicitPIC][AMR] WARNING: no electron species found; disabling Debye indicator.\n";
            amrCfg.debye.enabled = false;
        }
    }
    if (amrEnabled && !amrCfg.debye.enabled && !amrCfg.indicator.useEField && !geometryRefinementActive)
    {
        errors << "[ExplicitPIC][AMR] WARNING: no active refinement indicators; disabling AMR.\n";
        amrEnabled = false;
    }
}

double ExplicitPicAmrSetup::temperatureEV(const std::vector<SpeciesConfig>& species) const
{
    if (!electronSelected_) throw std::logic_error("Explicit PIC AMR temperature requires electron selection");
    const auto& amrCfg = config_;
    const bool amrEnabled = enabled_;
    const int amrElectronIndex = electronIndex_;
    double amrTe_eV = amrCfg.debye.Te_eV;
    if (amrEnabled && amrCfg.debye.enabled && amrTe_eV <= 0.0 && amrElectronIndex >= 0 &&
        static_cast<std::size_t>(amrElectronIndex) < species.size())
    {
        amrTe_eV = species[static_cast<std::size_t>(amrElectronIndex)].Te_eV;
    }
    if (amrEnabled && amrCfg.debye.enabled && amrTe_eV <= 0.0)
        amrTe_eV = 2.0; // fallback

    return amrTe_eV;
}

std::optional<PicCaseInput> readPicCaseInput(PicIntegrator integrator,
    const UnitsContext* units, const Solvers::SolversConfig* solvers)
{
    const bool implicit = integrator == PicIntegrator::Implicit;
    const char* name = implicit ? "ImplicitPIC" : "ExplicitPIC";
    DictionaryPIC dictionary;
    const std::string path = "./system/pic.dict";
    if (!dictionary.readDict(path, units))
    {
        std::cerr << '[' << name << "] ERROR: failed to read " << path << std::endl;
        return std::nullopt;
    }

    PicCaseInput input;
    input.source = path;
    input.pic = dictionary.config();
    input.geometryFilePresent = std::filesystem::exists("./system/geometry.dict");
    if (input.geometryFilePresent)
    {
        const auto geometry = [&]() {
            input.geometry.readDict(units);
            input.movingBoundary = readMovingEBConfig(input.geometry);
            input.gridMotion = readMotionSchedule(input.geometry);
            input.haveGeometry = true;
        };
        if (implicit)
        {
            try { geometry(); }
            catch (const std::exception& error)
            {
                std::cerr << "[ImplicitPIC][Geometry] WARNING: failed to read geometry.dict: "
                          << error.what() << '\n';
                input.haveGeometry = false;
            }
        }
        else geometry();
    }
    input.mesh.readDict(units);
    if (input.haveGeometry)
        input.mesh.ensureLevelSetCount(static_cast<int>(input.geometry.levelSetName.size()));

    if (solvers) input.poisson = solvers->poissonCellCentered;
    else if (!implicit)
    {
        // Preserve the explicit driver's old optional fallback. The implicit
        // driver uses default options when no solver settings are supplied.
        Solvers::DictionarySolvers dictionary;
        if (dictionary.readDict("./system/solvers.dict"))
            input.poisson = dictionary.config().poissonCellCentered;
    }
    if (implicit) input.implicitOverrides = readImplicitPicOverrides();
    return input;
}

namespace {
std::string lower(std::string value)
{
    for (auto& c : value) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return value;
}
std::optional<std::string> stringValue(const char* key)
{
    if (const auto* value = std::getenv(key)) return std::string(value);
    return std::nullopt;
}
std::optional<std::string> nonemptyString(const char* key)
{
    auto value = stringValue(key);
    return value && !value->empty() ? value : std::nullopt;
}
std::optional<bool> boolean(const char* key)
{
    if (const auto value = stringValue(key))
    {
        const auto token = lower(*value);
        if (token == "1" || token == "true" || token == "yes" || token == "on") return true;
        if (token == "0" || token == "false" || token == "no" || token == "off") return false;
    }
    return std::nullopt;
}
std::optional<bool> debugEnabled(const char* key)
{
    if (const auto value = stringValue(key))
    {
        const auto token = lower(*value);
        return !(token.empty() || token == "0" || token == "false" || token == "off" || token == "no");
    }
    return std::nullopt;
}
std::optional<int> integer(const char* key)
{
    if (const auto* value = std::getenv(key))
        try { return std::stoi(value); } catch (...) {}
    return std::nullopt;
}
std::optional<int> positiveInteger(const char* key)
{
    auto value = integer(key);
    if (value) *value = std::max(1, *value);
    return value;
}
std::optional<double> finiteReal(const char* key)
{
    if (const auto* value = std::getenv(key))
        try { const auto v = std::stod(value); if (std::isfinite(v)) return v; } catch (...) {}
    return std::nullopt;
}
// The old gates and temporal omega deliberately differ from finiteReal:
// infinity survives, NaN becomes zero, and malformed omega becomes 0.8.
std::optional<double> nonnegativeReal(const char* key, double invalid = 0.)
{
    if (const auto* value = std::getenv(key))
    {
        try { return std::max(0., std::stod(value)); } catch (...) { return invalid; }
    }
    return std::nullopt;
}
} // namespace

ImplicitPicAmrSetup resolveImplicitPicAmr(const PicCaseInput& input,
    int electronIndex, bool geometryRefinementActive)
{
    const auto& mesh = input.mesh;
    const auto& cfg = input.pic;
    const auto& overrides = input.implicitOverrides;
    const bool physicsEnabled = (mesh.geometryRefinement.present || mesh.physicsRefinement.present)
        && mesh.physicsRefinement.enabled;
    ImplicitPicAmrSetup result;
    result.adaptEvery = mesh.physicsRefinement.adaptEverySteps > 0 ? mesh.physicsRefinement.adaptEverySteps : 0;
    result.fastRemap = overrides.amrFastRemap.value_or(false);
    result.resetTransients = overrides.resetTransientsOnAmr.value_or(false);
    auto& indicators = result.indicators;
    indicators.enabled = physicsEnabled || geometryRefinementActive;
    indicators.indicatorEvery = mesh.physicsRefinement.indicatorEverySteps > 0
        ? mesh.physicsRefinement.indicatorEverySteps : std::max(1, result.adaptEvery);
    indicators.persistenceSteps = mesh.physicsRefinement.persistenceSteps > 0 ? mesh.physicsRefinement.persistenceSteps : 1;
    indicators.electronIndex = electronIndex;
    indicators.maxChanges = std::max(0, overrides.amrMaxChanges.value_or(0));
    indicators.temperatureEV = cfg.amr.debye.Te_eV > 0.0 ? cfg.amr.debye.Te_eV
        : ((electronIndex >= 0 && electronIndex < static_cast<int>(cfg.species.size()))
            ? std::max(1e-6, cfg.species[static_cast<std::size_t>(electronIndex)].Te_eV) : 2.0);
    indicators.debyeScale = cfg.amr.debye.CD > 0.0 ? cfg.amr.debye.CD : 0.5;
    indicators.electricPercentile = std::isfinite(cfg.amr.indicator.eFieldPercentile) && cfg.amr.indicator.eFieldPercentile >= 0.0
        ? cfg.amr.indicator.eFieldPercentile : 95.0;
    indicators.gateElectricRelative = overrides.amrGateEFieldRel.value_or(0.0);
    indicators.gateDebyeRelative = overrides.amrGateDebyeRel.value_or(0.0);
    if (physicsEnabled)
        for (const auto& name : mesh.physicsRefinement.indicators)
        {
            const auto key = lower(name);
            if (key == "debyelength") indicators.wantDebye = true;
            else if (key == "electricfield") indicators.wantEField = true;
        }

    const auto& meshRefresh = mesh.geometryRefresh;
    const auto& picRefresh = cfg.amr.geometryRefresh;
    result.meshRefreshPresent = meshRefresh.present;
    result.picRefreshPresent = picRefresh.present;
    const auto normalized = [](std::string token) {
        token = lower(token);
        token.erase(std::remove_if(token.begin(), token.end(), [](unsigned char c) { return std::isspace(c) != 0; }), token.end());
        return token;
    };
    result.refreshConflict = meshRefresh.present && picRefresh.present &&
        (normalized(meshRefresh.mode) != normalized(picRefresh.mode) ||
         std::abs(meshRefresh.changeRel - picRefresh.changeRel) > 1e-12 ||
         meshRefresh.verbose != picRefresh.verbose || meshRefresh.allowWithMovingEb != picRefresh.allowWithMovingEb ||
         normalized(meshRefresh.movingEbRefreshMode) != normalized(picRefresh.movingEbRefreshMode));
    auto& geometry = result.geometry;
    // Conflict comparison strips whitespace; runtime selection historically only
    // lowercases. Preserve that distinction for prepared/programmatic inputs.
    const auto mode = lower(overrides.amrGeomRefreshMode.value_or(meshRefresh.present ? meshRefresh.mode : picRefresh.mode));
    geometry.mode = mode == "skip" ? PicAmrGeometry::RefreshMode::Skip
        : (mode == "adaptive" ? PicAmrGeometry::RefreshMode::Adaptive : PicAmrGeometry::RefreshMode::Full);
    geometry.amrGeomRefreshChangeRelTh = std::max(0.0, overrides.amrGeomRefreshChangeRel.value_or(meshRefresh.present ? meshRefresh.changeRel : picRefresh.changeRel));
    geometry.amrGeomRefreshVerbose = overrides.amrGeomRefreshVerbose.value_or(meshRefresh.present ? meshRefresh.verbose : picRefresh.verbose);
    geometry.amrGeomRefreshAllowWithMovingEb = overrides.amrGeomRefreshAllowWithMovingEb.value_or(meshRefresh.present ? meshRefresh.allowWithMovingEb : picRefresh.allowWithMovingEb);
    result.exactMovingRequested = lower(overrides.amrMovingEbRefreshMode.value_or(meshRefresh.present ? meshRefresh.movingEbRefreshMode : picRefresh.movingEbRefreshMode)) == "exact_refresh";
    result.exactMovingAllowed = overrides.amrAllowExactMovingEbRefresh.value_or(false);
    geometry.amrReinitActiveMovingOnly = overrides.amrReinitActiveMovingOnly.value_or(false);
    geometry.amrReinitHybridSensitiveInactive = overrides.amrReinitHybridSensitiveInactive.value_or(false);
    geometry.amrReinitHybridUseDelta = overrides.amrReinitHybridUseDelta.value_or(true);
    geometry.amrReinitHybridNearZeroMin = std::max(0, overrides.amrReinitHybridNearZeroMin.value_or(1));
    geometry.amrReinitHybridNarrowBandMin = std::max(0, overrides.amrReinitHybridNarrowBandMin.value_or(1));
    geometry.amrReinitHybridNearZeroDeltaMin = std::max(0, overrides.amrReinitHybridNearZeroDeltaMin.value_or(1));
    geometry.amrReinitHybridNarrowBandDeltaMin = std::max(0, overrides.amrReinitHybridNarrowBandDeltaMin.value_or(1));
    geometry.amrReinitHybridIncludeMaxIter = overrides.amrReinitHybridIncludeMaxIter.value_or(true);
    geometry.amrReinitHybridMaxIterThreshold = std::max(1, overrides.amrReinitHybridMaxIterThreshold.value_or(500));
    geometry.amrReinitActiveMovingMaxStreak = std::max(0, overrides.amrReinitActiveMovingMaxStreak.value_or(0));
    return result;
}

void ImplicitPicAmrSetup::writeReport(std::ostream& output, std::ostream& errors, bool consoleInfo) const
{
    if (fastRemap)
        output << "[ImplicitPIC][AMR] fast remap enabled (LAMP_IMPLICIT_AMR_FAST_REMAP=1)" << std::endl;
    if (consoleInfo && picRefreshPresent)
    {
        if (meshRefreshPresent)
            output << "[ImplicitPIC][AMR] note: pic.dict/AMR/GeometryRefresh is deprecated; "
                      "runtime GeometryRefresh now follows mesh.dict/GeometryRefresh." << std::endl;
        else
            output << "[ImplicitPIC][AMR] note: pic.dict/AMR/GeometryRefresh is deprecated; "
                      "move this block to mesh.dict/GeometryRefresh." << std::endl;
    }
    if (refreshConflict)
        errors << "[ImplicitPIC][AMR] WARNING: GeometryRefresh is defined in both mesh.dict and pic.dict "
                  "with different values; runtime is using mesh.dict/GeometryRefresh.\n";
    if (exactMovingRequested && !exactMovingAllowed)
        errors << "[ImplicitPIC][AMR] WARNING: movingEbRefreshMode=exact_refresh is experimental and "
                  "now requires LAMP_IMPLICIT_AMR_ALLOW_EXACT_MOVINGEB_REFRESH=1. Falling back to full.\n";
    if (consoleInfo)
        output << "[ImplicitPIC][AMR] geometry refresh mode="
               << (geometry.mode == PicAmrGeometry::RefreshMode::Skip ? "skip" : geometry.mode == PicAmrGeometry::RefreshMode::Adaptive ? "adaptive" : "full")
               << " relTh=" << geometry.amrGeomRefreshChangeRelTh
               << " allowWithMovingEB=" << (geometry.amrGeomRefreshAllowWithMovingEb ? 1 : 0)
               << " movingEbRefreshMode=" << (exactMovingRequested && exactMovingAllowed ? "exact_refresh" : "full")
               << " reinitActiveMovingOnly=" << (geometry.amrReinitActiveMovingOnly ? 1 : 0)
               << " hybridSensitiveInactive=" << (geometry.amrReinitHybridSensitiveInactive ? 1 : 0)
               << " hybridUseDelta=" << (geometry.amrReinitHybridUseDelta ? 1 : 0)
               << " reinitActiveMovingMaxStreak=" << geometry.amrReinitActiveMovingMaxStreak
               << " (dict: " << (meshRefreshPresent ? "mesh.dict/GeometryRefresh" : picRefreshPresent ? "pic.dict/AMR/GeometryRefresh" : "built-in defaults via pic.dict/AMR/GeometryRefresh")
               << ", env override supported)" << std::endl;
}

ImplicitPicOverrides readImplicitPicOverrides()
{
    ImplicitPicOverrides result;
    result.initGeomSnapshotPath = nonemptyString("LAMP_IMPLICIT_INIT_GEOM_SNAPSHOT_PATH");
    result.seeDebug = debugEnabled("LAMP_IMPLICIT_SEE_DEBUG");
    result.packedFlushEvery = integer("LAMP_IMPLICIT_PACKED_FLUSH_EVERY");
    result.amrGateEFieldRel = nonnegativeReal("LAMP_IMPLICIT_AMR_GATE_EFIELD_REL");
    result.amrGateDebyeRel = nonnegativeReal("LAMP_IMPLICIT_AMR_GATE_DEBYE_REL");
    result.amrMaxChanges = integer("LAMP_IMPLICIT_AMR_MAX_CHANGES");
    result.initGeomSnapshot = boolean("LAMP_IMPLICIT_INIT_GEOM_SNAPSHOT");
    result.periodicDiagDebug = boolean("LAMP_IMPLICIT_PERIODIC_DIAG_DEBUG");
    result.consoleFallbackFracWarn = finiteReal("LAMP_IMPLICIT_CONSOLE_FALLBACK_FRAC_WARN");
    result.consoleFallbackDebounceSteps = integer("LAMP_IMPLICIT_CONSOLE_FALLBACK_DEBOUNCE_STEPS");
    result.consoleFallbackPrintRecover = boolean("LAMP_IMPLICIT_CONSOLE_FALLBACK_PRINT_RECOVER");
    result.consoleDqrelWarn = finiteReal("LAMP_IMPLICIT_CONSOLE_DQREL_WARN");
    result.consoleDqrelRecover = finiteReal("LAMP_IMPLICIT_CONSOLE_DQREL_RECOVER");
    result.consoleDqrelDebounceSteps = integer("LAMP_IMPLICIT_CONSOLE_DQREL_DEBOUNCE_STEPS");
    result.consoleDqrelEventCooldownSteps = integer("LAMP_IMPLICIT_CONSOLE_DQREL_EVENT_COOLDOWN_STEPS");
    result.consoleDqrelPrintRecover = boolean("LAMP_IMPLICIT_CONSOLE_DQREL_PRINT_RECOVER");
    result.consoleRgL2Warn = finiteReal("LAMP_IMPLICIT_CONSOLE_RG_L2_WARN");
    result.consoleAmrSkipStreakWarn = positiveInteger("LAMP_IMPLICIT_CONSOLE_AMR_SKIP_STREAK_WARN");
    result.consolePhysicsEnabled = boolean("LAMP_IMPLICIT_CONSOLE_PHYSICS_ENABLED");
    result.consolePhysicsEvery = integer("LAMP_IMPLICIT_CONSOLE_PHYSICS_EVERY");
    result.consolePhysicsShowCharge = boolean("LAMP_IMPLICIT_CONSOLE_PHYSICS_SHOW_CHARGE");
    result.consolePhysicsShowCurrent = boolean("LAMP_IMPLICIT_CONSOLE_PHYSICS_SHOW_CURRENT");
    result.consolePhysicsShowWall = boolean("LAMP_IMPLICIT_CONSOLE_PHYSICS_SHOW_WALL");
    result.consolePhysicsShowMcc = boolean("LAMP_IMPLICIT_CONSOLE_PHYSICS_SHOW_MCC");
    result.consolePhysicsShowPerf = boolean("LAMP_IMPLICIT_CONSOLE_PHYSICS_SHOW_PERF");
    result.depositFillOmpMinParticles = integer("LAMP_IMPLICIT_DEPOSIT_FILL_OMP_MIN_PARTICLES");
    result.orbitAvgFused = boolean("LAMP_IMPLICIT_ORBITAVG_FUSED");
    result.resetTransientsOnAmr = boolean("LAMP_IMPLICIT_RESET_TRANSIENTS_ON_AMR");
    result.cutcellStateDumpStep = integer("LAMP_IMPLICIT_CUTCELL_STATE_DUMP_STEP");
    result.cutcellStateDumpPath = stringValue("LAMP_IMPLICIT_CUTCELL_STATE_DUMP_PATH");
    result.cutcellStateDumpOnlyActive = boolean("LAMP_IMPLICIT_CUTCELL_STATE_DUMP_ONLY_ACTIVE");
    result.movingTraceStep = integer("LAMP_IMPLICIT_MOVING_TRACE_STEP");
    result.movingTraceAll = boolean("LAMP_IMPLICIT_MOVING_TRACE_ALL");
    result.movingTracePath = stringValue("LAMP_IMPLICIT_MOVING_TRACE_PATH");
    result.outputFastDensity = boolean("LAMP_IMPLICIT_OUTPUT_FAST_DENSITY");
    result.outputCellConstant = boolean("LAMP_IMPLICIT_OUTPUT_CELL_CONSTANT");
    result.outputDisableCellConstant = boolean("LAMP_IMPLICIT_OUTPUT_DISABLE_CELL_CONSTANT");
    result.outputMaskDiag = boolean("LAMP_IMPLICIT_OUTPUT_MASK_DIAG");
    result.packedBatchWithAmr = boolean("LAMP_IMPLICIT_PACKED_BATCH_WITH_AMR");
    result.amrFastRemap = boolean("LAMP_IMPLICIT_AMR_FAST_REMAP");
    result.amrGeomRefreshMode = stringValue("LAMP_IMPLICIT_AMR_GEOM_REFRESH_MODE");
    result.amrGeomRefreshChangeRel = finiteReal("LAMP_IMPLICIT_AMR_GEOM_REFRESH_CHANGE_REL");
    result.amrGeomRefreshVerbose = boolean("LAMP_IMPLICIT_AMR_GEOM_REFRESH_VERBOSE");
    result.amrGeomRefreshAllowWithMovingEb = boolean("LAMP_IMPLICIT_AMR_GEOM_REFRESH_ALLOW_WITH_MOVINGEB");
    result.amrAllowExactMovingEbRefresh = boolean("LAMP_IMPLICIT_AMR_ALLOW_EXACT_MOVINGEB_REFRESH");
    result.amrMovingEbRefreshMode = stringValue("LAMP_IMPLICIT_AMR_MOVINGEB_REFRESH_MODE");
    result.amrReinitActiveMovingOnly = boolean("LAMP_IMPLICIT_AMR_REINIT_ACTIVE_MOVING_ONLY");
    result.amrReinitHybridSensitiveInactive = boolean("LAMP_IMPLICIT_AMR_REINIT_HYBRID_SENSITIVE_INACTIVE");
    result.amrReinitHybridUseDelta = boolean("LAMP_IMPLICIT_AMR_REINIT_HYBRID_USE_DELTA");
    result.amrReinitHybridNearZeroMin = positiveInteger("LAMP_IMPLICIT_AMR_REINIT_HYBRID_NEARZERO_MIN");
    result.amrReinitHybridNarrowBandMin = positiveInteger("LAMP_IMPLICIT_AMR_REINIT_HYBRID_NARROWBAND_MIN");
    result.amrReinitHybridNearZeroDeltaMin = positiveInteger("LAMP_IMPLICIT_AMR_REINIT_HYBRID_NEARZERO_DELTA_MIN");
    result.amrReinitHybridNarrowBandDeltaMin = positiveInteger("LAMP_IMPLICIT_AMR_REINIT_HYBRID_NARROWBAND_DELTA_MIN");
    result.amrReinitHybridIncludeMaxIter = boolean("LAMP_IMPLICIT_AMR_REINIT_HYBRID_INCLUDE_MAXITER");
    result.amrReinitHybridMaxIterThreshold = positiveInteger("LAMP_IMPLICIT_AMR_REINIT_HYBRID_MAXITER_THRESHOLD");
    result.amrReinitActiveMovingMaxStreak = positiveInteger("LAMP_IMPLICIT_AMR_REINIT_ACTIVE_MOVING_MAX_STREAK");
    result.orbitAvgSpeciesPar = boolean("LAMP_IMPLICIT_ORBITAVG_SPECIES_PAR");
    result.segmenterMaxSegments = integer("LAMP_IMPLICIT_SEGMENTER_MAXSEG");
    result.segmenterMaxLevel = integer("LAMP_IMPLICIT_SEGMENTER_MAXLEVEL");
    result.polarization.mode = stringValue("LAMP_IMPLICIT_POL_MODE");
    result.polarization.temporalGuess = boolean("LAMP_IMPLICIT_POL_CG_TEMPORAL_GUESS");
    result.polarization.reuseIterativePattern = boolean("LAMP_IMPLICIT_POL_CG_REUSE_PATTERN");
    result.polarization.temporalOmega = nonnegativeReal("LAMP_IMPLICIT_POL_CG_TEMPORAL_OMEGA", 0.8);
    // Presence, including an empty value or "0", enabled these old debug flags.
    result.segmenterTrace = std::getenv("LAMP_SEGMENTER_TRACE") != nullptr;
    result.segmenterIgnoreHostHints = std::getenv("LAMP_SEGMENTER_NOHINT") != nullptr;
    return result;
}
