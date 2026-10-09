#pragma once

#include <cstddef>
#include "Modules/Plasma/Physics/PicMeshRefinement.h"
#include "Modules/Plasma/Boundary/PicAmrGeometry.h"
#include <deque>
#include "Modules/Plasma/Physics/PicSpawnManagement.h"
#include "Modules/Plasma/Physics/PicSecondaryEmission.h"
#include <filesystem>
#include <fstream>
#include <string_view>
#include <string>
#include <vector>
#include "Modules/Plasma/ExplicitPIC/Diagnostics/CurrentWaveformWriter.h"
#include "Modules/Plasma/Diagnostics/PicStepDiagnostics.h"
#include "Modules/Plasma/ExplicitPIC/Diagnostics/PicDiagnostics.h"
#include "Modules/Plasma/Diagnostics/PhysicsDiagnostics.h"

class UnitsContext;
struct ParticleManagementReport;
struct PicConfig;
struct ImplicitPicOverrides;
struct LegacyGridMotionState;

// Owns the PIC diagnostic CSV streams, preserving each legacy layout.
// Callers supply prepared values at the original diagnostic/cadence stage.
class PicStepSeries
{
public:
    enum class Mode { Explicit, Implicit };
    // Prescribed-motion IO. The caller binds the original startup/step stages
    // and forces initial/final/report rows; changes and cadence use the sample.
    void openGeometryMetrics(const std::filesystem::path& seriesDir,
        std::filesystem::path configuredPath, int resumeStep, int every);
    void writeGeometryMetrics(int step, double timeS, const LegacyGridMotionState& state, bool force);
    void openGeometryTrace(const std::filesystem::path& path, int selectedStep, bool allSteps);
    void writeGeometryTrace(int step, const std::string& message);
    // Explicit legacy reports. SEE starts a new text report even on restart;
    // impact and RF rows use the caller's original CSV resume threshold.
    void openSecondaryEmission(const std::filesystem::path& path);
    bool secondaryEmissionReady() const { return secondaryEmission_.is_open() && secondaryEmission_.good(); }
    void writeSecondaryEmission(int step, double timeS, const PicSecondaryEmission::Stats& row,
        const PicEmissionEnergy& energy, bool flush);
    void openImpact(const std::filesystem::path& path, int resumeStep);
    bool impactReady() const { return impact_.is_open() && impact_.good(); }
    void writeImpact(int step, double timeS, const PicBoundaryLosses::Totals& losses);
    void openRfPhase(const std::filesystem::path& path, int resumeStep, bool writeSI);
    bool rfPhaseReady() const { return rfPhase_.is_open() && rfPhase_.good(); }
    void writeRfPhase(int step, double time, const std::vector<PicRfPhase>& rows);
    void flushRfPhase() { if (rfPhaseReady()) rfPhase_.flush(); }
    // Retains precision/flush behavior; caller chooses the original debug cadence.
    // Console charge totals are the deposition totals, independent of Gauss sums.
    static void writeExplicitDiagnostic(std::ostream& stream, int step,
        double totalGridCharge, double totalMacroCharge,
        const ExplicitPIC::Diagnostics::PicDiagnostics::ExplicitSample& sample,
        std::string_view poissonMessage);

    // AMR IO consumes prepared reports at the caller's original sampling stage.
    // Explicit indicators retain default float formatting; implicit tables use
    // fixed precision. Each written row/batch flushes, including restart rows.
    void openAmrIndicators(const std::filesystem::path& path, int resumeStep);
    bool amrIndicatorsReady() const { return amrIndicators_.is_open() && amrIndicators_.good(); }
    void closeAmrIndicators() { if (amrIndicators_.is_open()) amrIndicators_.close(); }
    void writeAmrIndicators(int step, double timeS, const ExplicitPicRefinement::Report& report,
        bool pending, bool adapted);
    void openAmrAdaptation(const std::filesystem::path& path, int resumeStep);
    void openAmrLevelSets(const std::filesystem::path& path, int resumeStep);
    bool amrAdaptationReady() const { return amrAdaptation_.is_open() && amrAdaptation_.good(); }
    bool amrLevelSetsReady() const { return amrLevelSets_.is_open() && amrLevelSets_.good(); }
    void writeAmrAdaptation(int step, const ImplicitPicRefinement::AdaptationReport& adaptation,
        const PicAmrGeometry::Report& geometry);
    void writeAmrLevelSets(int step, const PicAmrGeometry::Geometry* geometry,
        const PicAmrGeometry::Report& report);
    static void writeAmrDebyeWarning(std::ostream& stream, const ExplicitPicRefinement::Report& report);
    static void writeAmrConsole(std::ostream& stream, int step, int cellCount, double minimumLeafWidth,
        const ExplicitPicRefinement::AdaptationReport& adaptation);
    static void writeAmrConsole(std::ostream& stream, int step,
        const ImplicitPicRefinement::AdaptationReport& adaptation, const PicAmrGeometry::Report& geometry);

    template<class MacroCount> struct OverviewValues
    {
        double fieldEnergyJ, kineticElectronJ, kineticIonJ, kineticNeutralJ, kineticTotalJ;
        double totalEnergyJ, conductionCurrent, displacementCurrent, totalCurrent;
        MacroCount macroTotal;
        double physicalTotal, chargeBalanceRelative;
    };
    struct ExplicitOverview
    {
        double residualL2, residualLinf, residualL2PerActive;
        double residualL2RelativeRhs, residualLinfRelativeRhs, embeddedRelativeRhs;
        OverviewValues<double> values;
    };
    struct ImplicitOverview
    {
        double residualL2, residualLinf;
        OverviewValues<long long> values;
    };
    struct Process
    {
        double driveEnergyJ, injectionEnergyJ, reactionEnergyLossJ;
        double domainAbsorptionJ, embeddedAbsorptionJ, seeEmissionEnergyJ, ionImpactEnergyJ;
        double displacementCurrentFraction;
        std::string_view dominantSource;
        double dominantSourceFraction;
        std::string_view dominantBoundary;
        double dominantBoundaryPowerW, dominantBoundaryFraction;
        std::string_view dominantReactionChannel;
        double dominantReactionPowerW, dominantReactionFraction;
    };
    struct Geometry
    {
        int updated;
        std::size_t activeCutCellCount;
        double embeddedSurfaceAreaM2;
        std::size_t movingRegionCount, movingLevelsetCount;
        double movingRegionTotalAreaM2, movingRegionTotalFlux;
        double movingRegionMeanEnergyEV, movingRegionMeanCos;
        double movingRegionMeanAbsSpeedMps, movingRegionMaxAbsSpeedMps;
        double geometryGsaMm, geometryGadMm, geometryDaMm;
        double fieldEnergyJ, domainAbsorptionJ, embeddedAbsorptionJ, reactionEnergyLossJ;
    };
    struct Power
    {
        double residualL2, residualLinf, residualL2PerActive, residualL2RelativeRhs, residualLinfRelativeRhs;
        double coupling, embedded, embeddedRelativeRhs;
        std::array<double, 3> couplingPercentiles, embeddedPercentiles;
        double fieldEnergy, kineticElectron, kineticIon, kineticNeutral, kineticTotal;
        PicEnergyBalance balance;
        PicEnergyTransfers transfers;
        double potentialCurrent, potentialEnergy;
    };

    struct ExplicitPerformance
    {
        PicExplicitStepTiming field;
        std::array<double, 4> poissonSeconds; // assemble, factor, solve, total
        double amrAdaptSeconds, managementSeconds, gatherPushSeconds, outputSeconds, totalSeconds;
        std::uint64_t particles;
        int poissonN;
        long long poissonNnz;
        bool usedSuperLU;
    };
    struct ImplicitPerformance
    {
        std::uint64_t particles;
        bool outputStep, amrAdaptStep, amrSkippedByGate;
        int poissonN;
        long long poissonNnz;
        std::string_view solverMode;
    };
    // Explicit time is code time; implicit time is seconds. Explicit restart
    // migrates the known legacy header/short initial row before trimming history.
    void openPerformance(const std::filesystem::path& path, int resumeStep);
    bool performanceReady() const { return performance_.is_open() && performance_.good(); }
    void writePerformance(int step, double timeCode, const ExplicitPerformance& row);
    void writePerformance(int step, double timeS, const PicImplicitStepTiming& timing,
        const PicImplicitTimingSummary& summary, const ImplicitPerformance& row);
    void openCouplingFallback(const std::filesystem::path& path, int resumeStep);
    bool couplingFallbackReady() const { return fallback_.is_open() && fallback_.good(); }
    void writeCouplingFallback(int step, double timeS, const PicImplicitStepTiming& timing,
        const PicImplicitTimingSummary& summary);
    void flushPerformance(); // Implicit cadence remains owned by the integrator.

    // Implicit diagnostic streams. The integrator keeps the original stage,
    // cadence and final-step gates; the writer owns headers/restart/buffering.
    void openPhysics(const std::filesystem::path& path, int resumeStep);
    bool physicsReady() const { return physics_.is_open() && physics_.good(); }
    void writePhysics(int step, double timeS, const ImplicitOverview& row,
        const PicStepParticleBalance::Result& charge);
    void openAmount(const std::filesystem::path& path, int resumeStep);
    void openAmountSpecies(const std::filesystem::path& path, int resumeStep);
    bool amountReady() const { return amount_.is_open() && amount_.good(); }
    bool amountSpeciesReady() const { return amountSpecies_.is_open() && amountSpecies_.good(); }
    void writeAmount(int step, double timeS, const PicStepParticleBalance::Result& charge);
    void writeAmountSpecies(int step, double timeS, std::string_view species, const PicSpeciesAmount& row);
    // Mode preserves the distinct band schemas, formatting and open failures.
    // Both integrators retain their original sampling and flush stages.
    void openEmbeddedBand(const std::filesystem::path& path, int resumeStep);
    bool embeddedBandReady() const { return embeddedBand_.is_open() && embeddedBand_.good(); }
    void writeEmbeddedBand(int step, double timeS, const PicEmbeddedBandField& row);
    void flushPhysics() { if (physicsReady()) physics_.flush(); }
    void flushAmounts();
    void flushEmbeddedBand() { if (embeddedBandReady()) embeddedBand_.flush(); }

    explicit PicStepSeries(Mode mode) : mode_(mode) {}
    PicStepSeries(const PicStepSeries&) = delete;
    PicStepSeries& operator=(const PicStepSeries&) = delete;
    // Explicit detail output opens lazily at the first diagnostic row. One file,
    // resume step and unit mode per session; units are borrowed only for a write.
    // Existing headers must match before history is trimmed. Later rows retain
    // the initial species order, filling absent species/values with zero.
    void writeDetail(const std::filesystem::path& path, int step, double time,
        const ExplicitPIC::Diagnostics::PicDiagnostics::Result& result,
        const UnitsContext* units = nullptr, bool writeSI = false, int resumeStep = 0);
    void openParticleManagement(const std::filesystem::path& path, int resumeStep);
    bool particleManagementReady() const { return particleManagement_.is_open() && particleManagement_.good(); }
    void writeSpawnManagement(int step, double timeS, std::string_view name,
        const PicSpawnManagement::Stats& stats, double keepProbability);
    void writeParticleManagement(int step, double timeS, std::string_view name,
        const ParticleManagementReport& report);
    void flushParticleManagement() { particleManagement_.flush(); }
    void openOverview(const std::filesystem::path& path, int resumeStep);
    void openProcess(const std::filesystem::path& path, int resumeStep);
    void openGeometry(const std::filesystem::path& path, int resumeStep);
    void openCurrentBalance(const std::filesystem::path& path, int resumeStep);
    bool currentBalanceReady() const { return currentBalance_.is_open() && static_cast<bool>(currentBalance_); }
    void writeCurrentBalance(int step, double timeS, const PicExplicitCurrentBalance& row);
    void writeCurrentBalance(int step, double timeS, const PicStepParticleBalance::Result& row);
    void flushCurrentBalance() { currentBalance_.flush(); }
    void openPower(const std::filesystem::path& path, int resumeStep, bool writeSI = true);
    bool powerReady() const { return power_.is_open() && static_cast<bool>(power_); }
    void writePower(int step, double time, const Power& row);
    void flushPower() { power_.flush(); }
    void openPeriodic(const std::filesystem::path& path, int resumeStep);
    void openWaveform(const std::filesystem::path& path, int resumeStep, int every);
    bool waveformDue(int step) const { return waveform_.isOpen() && step % waveformEvery_ == 0; }
    void writeWaveform(int step, double timeS, double voltageV, double conduction, double displacement);
    bool periodicReady() const { return periodic_.is_open() && static_cast<bool>(periodic_); }
    bool overviewReady() const { return overview_.is_open() && static_cast<bool>(overview_); }
    bool processReady() const { return process_.is_open() && static_cast<bool>(process_); }
    bool geometryReady() const { return geometry_.is_open() && static_cast<bool>(geometry_); }

    // The overview layouts remain distinct; implicit rows retain their buffering.
    void writeOverview(int step, double time_s, const ExplicitOverview& row);
    void writeOverview(int step, double time_s, const ImplicitOverview& row);
    // String views are borrowed only during this synchronous write, not retained.
    void writeProcess(int step, double time_s, const Process& row);
    void writeGeometry(int step, double time_s, const Geometry& row);

    // Explicit rows always flush; implicit rows flush on the caller's report
    // step, with a separate flush at the original end-of-step IO checkpoint.
    void writePeriodic(int step, double time_s, double jumpX_V, double jumpY_V,
        bool periodicX, bool periodicY, bool reportStep = false);
    void flushPeriodic() { periodic_.flush(); }

private:
    void open(std::ofstream& stream, const std::filesystem::path& path, int resumeStep,
        const char* header, const char* role);
    template<class MacroCount> void writeOverviewValues(const OverviewValues<MacroCount>& values);
    const Mode mode_;
    std::ofstream overview_, process_, geometry_, periodic_, power_, currentBalance_;
    std::ofstream detail_, particleManagement_, performance_, fallback_;
    std::ofstream physics_, amount_, amountSpecies_, embeddedBand_;
    std::ofstream amrIndicators_, amrAdaptation_, amrLevelSets_;
    std::ofstream secondaryEmission_, impact_, rfPhase_;
    std::ofstream geometryMetrics_, geometryTrace_;
    int geometryMetricsEvery_ = 1, geometryTraceStep_ = -1;
    bool geometryTraceAll_ = false;
    std::filesystem::path detailPath_;
    std::vector<std::string> detailSpeciesNames_;
    bool detailWriteSI_ = false;
    int detailRowsSinceFlush_ = 0;
    bool powerWriteSI_ = true;
    ExplicitPIC::Diagnostics::CurrentWaveformWriter waveform_;
    int waveformEvery_ = 1;
};

// Per-run console configuration, sampled diagnostics and warning/ETA history.
// Step inputs are borrowed synchronously; reporting never advances particles,
// geometry or the physical clock. Restart begins a fresh console history.
class PicImplicitConsole
{
public:
    struct Walls
    {
        const std::array<double, 4>& electronCurrent;
        const std::array<double, 4>& ionCurrent;
        const std::array<double, 4>& electronCharge;
        const std::array<double, 4>& ionCharge;
    };
    struct Events
    {
        bool amrIndicator, amrSkippedByGate, amrAdapted, embeddedReady, collisionsEnabled;
        long long sourceSpawnCount, collisionSpawnCount;
    };
    struct Sample
    {
        const QuadtreeMesh& mesh;
        const CellField<double>& potential;
        const UnitsContext* units;
        const std::vector<ParticleType>& species;
        const std::vector<double>& defaultWeights;
        const PicStepParticleBalance::Result& charge;
        const PicBoundaryLosses::Totals& losses;
        Walls walls;
        const PhysicsDiagnostics::ElectrodeCurrentGroups& electrodes;
        const PicElectrodeCurrent::Sample& current;
        const MCCCollisionOperator::StepStats& collisions;
        Events events;
    };
    PicImplicitConsole(const PicConfig& config, const ImplicitPicOverrides& overrides,
        int firstStep, int endStep, int fieldEvery, bool infoEnabled);
    PicImplicitConsole(const PicImplicitConsole&) = delete;
    PicImplicitConsole& operator=(const PicImplicitConsole&) = delete;
    void observeDiagnostics(const ExplicitPIC::Diagnostics::PicDiagnostics::Result& diagnostic)
    { diagLast = diagnostic; haveDiagLast = true; }
    // step is the integrator's zero-based completed step. Returns the current
    // macro-particle count also needed by the separate performance CSV.
    std::uint64_t writeStep(int step, double timeS, double dtS, bool defaultSummaryStep,
        const PicImplicitStepTiming& timing, const PicImplicitTimingSummary& timingSummary,
        const Sample& sample, std::ostream& output);
private:
    int startStep, steps, fieldsEvery;
    bool consoleInfo;
    double fallbackFracWarnTh, dQrelWarnTh, dQrelRecoverTh, rGL2PerActiveWarnTh;
    int fallbackFracDebounceSteps, dQrelDebounceSteps, dQrelEventCooldownSteps, amrSkipWarnStreakTh;
    bool fallbackPrintRecover, dQrelPrintRecover, physicsSummaryEnabled;
    int physicsSummaryEvery;
    bool physicsShowCharge, physicsShowCurrent, physicsShowWall, physicsShowMCC, physicsShowPerf;
    std::string physicsDetail;
    bool summaryOneLine;
    ExplicitPIC::Diagnostics::PicDiagnostics::Result diagLast{};
    bool haveDiagLast = false, dQrelWarnActive = false, rGL2WarnActive = false, fallbackWarnActive = false;
    int dQrelWarnConsec = 0, dQrelRecoverConsec = 0, lastDqrelEventStep = -1;
    int fallbackWarnConsec = 0, fallbackRecoverConsec = 0, amrSkipGateStreak = 0;
    std::deque<double> recentStepCost_s;
    double recentStepCostSum_s = 0.0;
};

// Explicit progress formatting and per-run configuration. The shared latest
// diagnostic remains with its sampling/CSV caller; current physical inputs and
// that optional diagnostic are borrowed only for this write, never retained.
class PicExplicitConsole
{
public:
    struct Emission
    {
        long long impactMacro, emittedMacro;
        double impactPhysical, emittedPhysical, meanEV, stdEV;
    };
    struct Events
    {
        long long sourceSpawn;
        double sourceCharge, domainCharge, embeddedCharge;
        bool seeEnabled, embeddedReady, collisionsEnabled, adapted, outputWritten;
        std::uint64_t embeddedResidual;
    };
    struct Sample
    {
        PicExplicitConsoleInput physical;
        const PhysicsDiagnostics::ElectrodeCurrentGroups& electrodes;
        const ExplicitPIC::Diagnostics::PicDiagnostics::Result* diagnostic;
        Emission emission;
        Events events;
    };
    PicExplicitConsole(const PicConfig& config, int firstStep, int endStep);
    // Original one-based completed step; the caller supplies its existing ETA
    // and measured timings. Disabled reports do not perform physical reductions.
    void writeStep(int step, double timeCode, double dtCode, double etaSeconds,
        const PicStepSeries::ExplicitPerformance& timing, const Sample& sample, std::ostream& output) const;
private:
    int startStep, steps, physicsSummaryEvery;
    bool physicsSummaryEnabled, physicsShowCharge, physicsShowCurrent, physicsShowWall, physicsShowMCC, physicsShowPerf;
    std::string physicsDetail;
    bool summaryOneLine;
};
