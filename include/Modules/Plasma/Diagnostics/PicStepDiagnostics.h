#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <optional>
#include <string>
#include <string_view>
#include <limits>
#include <cmath>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "Modules/Plasma/Physics/MCCCollisionOperator.h"
#include "Solvers/Poisson/CutCellMetrics.h"
#include "Solvers/Poisson/PoissonBC.h"

class MovingEBController;
class ParticleType;
class PicSpeciesState;
class PicFieldBoundaryState;
class Particle;
class ParticleSpawnBuffer;
class UnitsContext;
class Vector2D;
namespace Solvers { namespace Poisson { struct DomainBC; } }
class QuadtreeMesh;
class FaceSegmentGraph;
namespace Solvers { namespace Poisson { class CutCellMetrics; } }
struct FaceSegment;
enum class Boundary;
template<class T> class CellField;
class NormalizationContext;
class LevelSet;

// Explicit report samples own their values. No particle reordering, host repair,
// RNG use or retained mesh/particle/boundary references.
struct PicParticleWeightStatistics
{
    std::size_t macro = 0;
    double physical = 0.0, minimum = 0.0, mean = 0.0;
    std::array<double, 3> percentiles{};
    double maximum = 0.0;
    int cellMacroMaximum = 0;
    double cellPhysicalMaximum = 0.0;
    int cellMacroIndex = -1, cellPhysicalIndex = -1;
};
PicParticleWeightStatistics summarizePicParticleWeights(const QuadtreeMesh& mesh,
    const ParticleType& species, double defaultWeight);

struct PicEmissionEnergy { double meanEV, stdEV; };
PicEmissionEnergy summarizePicEmissionEnergy(double physical, double energyJ, double energyJ2);

struct PicRfPhase
{
    std::string target, type;
    double frequencyHz, radians, degrees, potential;
};
// Evaluate at SI time; a null unit context leaves the potential in code units.
// Report the original domain/default embedded waveforms, in the original order.
std::vector<PicRfPhase> summarizePicRfPhases(const PicDomainBC& domain,
    const PicEmbeddedBC& embedded, double timeSI, const UnitsContext* outputUnits);

// Synchronous, read-only cut-cell diagnostic. All inputs and sampled boundary
// pointers are borrowed; neither the view nor its rows may outlive them or AMR.
// Only the embedded normal-field vector is prepared for the whole mesh.
class PicCutcellStateView
{
public:
    struct Input {
        const QuadtreeMesh& mesh;
        const CellField<double>& density;
        const std::vector<double>& kappa;
        const CellField<double>& potential;
        const CellField<Vector2D>& electric;
        const Solvers::Poisson::CutCellMetrics* metrics = nullptr;
        const PicEmbeddedBC* embedded = nullptr;
        const PicEmbeddedBCMap* embeddedById = nullptr;
        const std::vector<LevelSet>* levelSets = nullptr;
        const UnitsContext* units = nullptr;
        const NormalizationContext* normalization = nullptr;
        bool skipMissingEmbeddedBC = false;
    };
    struct Row {
        int ci, cLabel, level, boundary;
        Vector2D xc;
        double x_m, y_m, areaCode, area_m2;
        bool active;
        Solvers::Poisson::CutCellMetrics::CellMetric cm;
        double rhoCode, rhoSI, kappa, phiCode, phiV;
        Vector2D eCode;
        double exVpm, eyVpm, enGammaCode, enGammaVpm;
        const PicEmbeddedBC* eb;
        const char* ebSource;
        double aGamma_m, dGamma_m;
        std::vector<double> psi;
    };
    explicit PicCutcellStateView(const Input& input);
    int size() const;
    int levelSetCount() const;
    bool active(int cellIndex) const;
    Row sample(int cellIndex) const;
private:
    const Input input_;
    const std::vector<double> normalField_;
    const bool useEmbeddedMap_, useEmbeddedFallback_;
};

// Legacy implicit stage-amount comparison, in physical particle counts.
// A missing or differently sized amount field contributes zero.
struct PicSpeciesAmount
{
    double stage, physical, relative;
};
PicSpeciesAmount summarizePicSpeciesAmount(const CellField<double>* amount,
    int cellCount, double physicalCount);

// Read-only near-wall field reports. No geometry or field state is retained.
struct PicEmbeddedBandField
{
    int layers;
    long long cells;
    double rmsVpm, meanVpm;
    int cutCells = 0; // Explicit report only; implicit seeds use alpha instead.
};
// Explicit: isCut seeds, interior faces, no nonfinite filtering; SI conversion
// follows the code-unit reduction. Missing/stale metrics or field yield zeros.
PicEmbeddedBandField summarizeExplicitPicEmbeddedBandField(const QuadtreeMesh& mesh,
    const FaceSegmentGraph& faces, const Solvers::Poisson::CutCellMetrics* metrics,
    int layers, const CellField<Vector2D>& field, double fieldScale);
// Implicit: finite alpha < 1 - 1e-12 seeds. Missing/nonfinite fields contribute
// zero to the numerator but remain in the selected-band denominator.
PicEmbeddedBandField summarizePicEmbeddedBandField(const QuadtreeMesh& mesh,
    const FaceSegmentGraph& faces, bool useCutMetrics, const std::vector<double>& alpha,
    int layers, const CellField<Vector2D>& field, double fieldScale);

// One macro-step's absorbed-particle accounting. No particle mutation, RNG,
// geometry lookup or IO. The integrator selects the original impact phase.
class PicBoundaryLosses
{
public:
    enum class Mode { Explicit, Implicit };
    struct EmbeddedCounts
    {
        std::uint64_t inside = 0, crossing = 0, unknown = 0;
        std::unordered_map<int, std::uint64_t> byId;
    };
    struct Totals
    {
        std::array<std::uint64_t, 4> domainMacroE{}, domainMacroI{}, periodicLeakE{}, periodicLeakI{};
        std::uint64_t unknownE = 0, unknownI = 0;
        std::array<double, 4> domainChargeE{}, domainChargeI{};
        std::vector<std::array<std::uint64_t, 4>> domainMacro;
        std::vector<std::array<double, 4>> domainPhysical, domainEnergy;
        std::vector<long long> embeddedMacro;
        std::vector<double> embeddedPhysical, embeddedEnergy;
        double domainEnergyJ = 0.0, embeddedEnergyJ = 0.0;
        std::unordered_map<int, double> embeddedChargeById;
        EmbeddedCounts embeddedCounts;
    };
    struct Sample
    {
        bool accepted;
        double weight, chargeC, energyPerParticleJ;
    };
    struct ImplicitWall
    {
        std::array<double, 4> areaM2{}, electronCharge{}, ionCharge{}, electronCurrent{}, ionCurrent{};
    };
    PicBoundaryLosses(std::size_t speciesCount, double velocityScale,
        bool periodicX, bool periodicY, Mode mode);
    const Totals& totals() const { return totals_; }
    // side: left/right/bottom/top = 0/1/2/3; other values are unknown.
    // Explicit counts periodic leaks; implicit ignores them. Both exclude
    // these impacts from physical outflow, energy and downstream observers.
    Sample recordDomain(std::size_t speciesIndex, const ParticleType& species,
        const Particle& particle, double defaultWeight, int side);
    Sample recordEmbedded(std::size_t speciesIndex, const ParticleType& species,
        const Particle& particle, double defaultWeight);
    // Resolve geometry outside this object, after recording the particle.
    // Explicit inside/crossing/ID counts retain the console sampling gate;
    // implicit ID counts remain unconditional for positive boundary IDs.
    void recordEmbeddedBoundary(int boundaryId, const Sample& sample,
        bool crossing, bool explicitClassification);
    // Sample after all impacts: implicit charge is reduced by species, unlike
    // explicit event-wise charge. Own values only; preserve SI area/dt fallbacks.
    ImplicitWall summarizeImplicitWall(const std::vector<ParticleType>& species,
        const QuadtreeMesh& mesh, const UnitsContext* units,
        const NormalizationContext* normalization, double dtSI) const;
private:
    Sample sample(const ParticleType& species, const Particle& particle, double defaultWeight) const;
    bool countEnergy(const ParticleType& species, double weightedEnergy) const;
    Totals totals_;
    double velocityScale_;
    bool periodicX_, periodicY_;
    Mode mode_;
};

// Preserve the existing MCC report's planar electron energy proxy: vx/vy,
// negative charges, threshold 11.59 eV. This is not a full 3D EEDF reduction.
struct PicCollisionElectronEnergy
{
    double meanEV, tailFraction;
    long long macroCount;
    double physicalCount;
};
PicCollisionElectronEnergy summarizePicCollisionElectronEnergy(
    const std::vector<ParticleType>& species, const std::vector<double>& weights, double velocityScale);

// Read-only values for the legacy PIC progress reports. Their electron
// energy uses all three velocity components; its channel-name classification
// remains distinct from the collision model's physical channel dispatch.
struct PicConsoleValues
{
    struct ProjectileGroups
    {
        const std::unordered_set<std::string>& electrons;
        const std::unordered_set<std::string>& ions;
    };
    static constexpr double tailThresholdEV = 11.59;
    struct MccSplit
    {
        std::uint64_t att = 0, coll = 0, nul = 0, spawn = 0;
        std::uint64_t el = 0, ex = 0, ion = 0, at = 0, cx = 0, other = 0;
        double maxNuDt = 0.0;
    };
    double phiMinV, phiMaxV, meanElectronEV, temperatureEV, electronTailFraction;
    std::uint64_t elastic, excitation, ionization, attachment, chargeExchange, other;
    MccSplit electrons, ions, neutrals;
};
// With no configured groups, classify by the first current species of that
// name (implicit). Explicit uses its configured sets, with electrons first.
PicConsoleValues summarizePicConsole(const QuadtreeMesh& mesh,
    const CellField<double>& potential, const UnitsContext* units,
    const std::vector<ParticleType>& species, const std::vector<double>& defaultWeights,
    const MCCCollisionOperator::StepStats& collisions,
    const PicConsoleValues::ProjectileGroups* configuredGroups = nullptr);

struct PicExplicitConsoleInput
{
    const QuadtreeMesh& mesh;
    const CellField<double>& potential;
    const UnitsContext* units;
    const PicSpeciesState& particles;
    const MCCCollisionOperator::StepStats& collisions;
    const PicBoundaryLosses::Totals& losses;
    double gridCharge, particleCharge, widthCode, heightCode, slabMeters;
    bool embeddedStatistics;
};
struct PicExplicitConsoleValues
{
    struct Population { std::uint64_t macro; double physical; };
    PicConsoleValues state;
    std::vector<Population> species;
    std::uint64_t macroCount = 0, wallLossMacro = 0, embeddedAbsorbed = 0;
    std::array<double, 4> electronCurrent{}, ionCurrent{};
    double totalCurrent = 0.0, chargeDifference = 0.0, chargeRelative = 0.0;
};
// Current population, deposition charge and explicit side-current rules;
// no diagnostic sampling, particle mutation, persistent references or IO.
PicExplicitConsoleValues summarizeExplicitPicConsole(const PicExplicitConsoleInput& input, double dtSI);

struct PicNeutralEnergy
{
    std::size_t macroCount;
    double physicalCount, totalJ, meanEV;
};
// Neutral report includes all three velocity components and the existing
// per-particle/default weight fallback. Both reductions are read-only.
PicNeutralEnergy summarizePicNeutralEnergy(const ParticleType& species, double defaultWeight, double velocityScale);

// Borrow the original per-species and per-channel accumulators without copying
// or reordering them. Arrays must contain at least speciesCount entries.
struct PicProcessSummaryInput
{
    double conductionCurrent, displacementCurrent;
    double driveEnergyJ, injectionEnergyJ, domainAbsorptionJ, embeddedAbsorptionJ;
    double reactionEnergyLossJ, timeWindowS;
    std::size_t speciesCount;
    const std::vector<std::array<double, 4>>& domainAbsorptionBySpecies;
    const std::vector<double>& embeddedAbsorptionBySpecies;
    const std::unordered_map<std::string, double>& reactionEnergyByChannel;
};

struct PicProcessSummary
{
    double displacementFraction;
    std::string dominantSource;
    double dominantSourceFraction;
    const char* dominantBoundary; // Static legacy side name, never a borrowed temporary.
    double dominantBoundaryPowerW, dominantBoundaryFraction;
    std::string dominantReactionChannel;
    double dominantReactionPowerW, dominantReactionFraction;
};

struct PicGeometrySummary
{
    std::size_t activeCutCellCount, movingRegionCount, movingLevelsetCount;
    double movingRegionTotalAreaM2, movingRegionTotalFlux;
    double movingRegionMeanEnergyEV, movingRegionMeanCos;
    double movingRegionMeanAbsSpeedMps, movingRegionMaxAbsSpeedMps;
    double geometryGsaMm, geometryGadMm, geometryDaMm;
};

// Pure diagnostics: no state mutation, RNG, unit conversion or file IO.
// Callers select the original step/window and retain output readiness/cadence.
PicProcessSummary summarizePicProcesses(const PicProcessSummaryInput& input);
PicGeometrySummary summarizePicGeometry(bool cutMetricsReady,
    const std::vector<double>& cutAlphaByCellIndex, int cellCount,
    const MovingEBController* movingEb, int step);

// Code-unit, area-weighted normal fields. Values and segments are paired up to
// their common length; nonpositive/NaN areas are skipped. No field sanitizing.
double averagePicBoundaryNormalField(const std::vector<FaceSegment>& segments,
    const std::vector<double>& normalFields, Boundary boundary);
// Fixed order: west, east, south, north (distinct from FaceSide's enum order).
std::array<double, 4> averagePicBoundaryNormalFields(const std::vector<FaceSegment>& segments,
    const std::vector<double>& normalFields);

// Explicit PIC accepts any located mate; implicit PIC requires tangentially
// aligned cell centers. Preserve these policies when measuring the same mesh.
enum class PicPeriodicMatePolicy { AnyMate, CleanMate };
struct PicPeriodicJumps
{
    double maxJumpX = 0.0, maxJumpY = 0.0;
    std::string debugX, debugY; // Owned text; no borrowed mesh or field state.
};
// Unique cell pairs per axis, fully open faces and uncut cells when metrics
// match this graph/mesh. Debug ties retain the last unique pair in graph order.
PicPeriodicJumps measurePicPeriodicJumps(const QuadtreeMesh& mesh,
    const FaceSegmentGraph& faceGraph, const CellField<double>& phi,
    bool periodicX, bool periodicY, const Solvers::Poisson::CutCellMetrics* cutMetrics,
    PicPeriodicMatePolicy matePolicy, bool debug = false);

// Owns the selected electrode's sampled fields, conduction current and previous
// normal field. Independent of CSV enablement/cadence; no mesh/particle aliases
// or checkpoint state are retained. A resumed run starts with no previous field.
class PicElectrodeCurrent
{
public:
    struct Sample
    {
        bool ready = false;
        double voltageV = 0.0;
        double conduction = std::numeric_limits<double>::quiet_NaN();
        double displacement = std::numeric_limits<double>::quiet_NaN();
        double total = std::numeric_limits<double>::quiet_NaN();
    };
    PicElectrodeCurrent(std::string_view side, const Solvers::Poisson::DomainBC& boundaries);
    Boundary boundary() const { return boundary_; }
    int sideIndex() const { return sideIndex_; } // west/east/south/north; -1 if unresolved
    void sampleNormalFieldCode(const FaceSegmentGraph& graph, const std::vector<double>& normalFields);
    double normalFieldCode() const { return normalCode_; }
    void sampleFieldSI(double normalFieldVpm, double voltageV);
    // Use the previously sampled normal field and evaluate this electrode's
    // voltage at the supplied SI time. Without units, explicit passes its E0
    // fallback; absence retains the implicit raw field (no extra multiply).
    void sampleFieldSI(const PicDomainBC& boundaries, double timeSI, const UnitsContext* units,
        std::optional<double> fieldScaleWithoutUnits = std::nullopt);
    void sampleConductionSI(const std::vector<ParticleType>& species,
        const std::vector<std::array<double, 4>>& physicalOutflow, double areaM2, double dtS);
    // Complete every original bookkeeping step, including steps without output.
    // Current/time-window validity follows the legacy rules, not CSV readiness.
    void finishStep(double dtS);
    const Sample& current() const { return current_; }
private:
    const Boundary boundary_;
    const int sideIndex_;
    double normalCode_ = 0.0, normalSI_ = 0.0, previousNormalSI_ = 0.0;
    double voltageSI_ = 0.0, conductionSI_ = 0.0;
    bool havePreviousNormal_ = false;
    Sample current_;
};

// Per-step electrode work, separate from the selected-electrode waveform.
// Owns only SI samples, four-side displacement history and scalar results.
// Callers retain voltage evaluation time, enabled gate and reset/finish points.
class PicBoundaryWork
{
public:
    struct Result { double energyJ = 0.0, powerW = 0.0, potentialCurrentW = 0.0; };
    struct Injection
    {
        const std::array<double, 4>& sourceChargeC;
        const std::array<double, 4>& secondaryChargeC;
    };
    explicit PicBoundaryWork(bool trackPotentialCurrentSum) : trackPotentialCurrentSum_(trackPotentialCurrentSum) {}
    void sampleNormalFieldSI(std::size_t side, double fieldVpm) { normalSI_[side] = fieldVpm; }
    void sampleDirichlet(const std::array<bool, 4>& selected) { dirichlet_ = selected; }
    void sampleVoltagesSI(const std::array<double, 4>& voltageV) { voltageSI_ = voltageV; }
    // Separate stages preserve explicit/implicit sampling order and gates.
    // An absent time consumes already evaluated implicit boundary values.
    void sampleNormalFieldsSI(const FaceSegmentGraph& graph, const std::vector<double>& normalFields,
        const UnitsContext& units);
    void sampleVoltagesSI(const PicDomainBC& boundaries, const UnitsContext& units,
        std::optional<double> evaluateAtTimeSI = std::nullopt);
    void resetStep(double dtS) { dtS_ = dtS; result_ = Result{}; }
    // Borrow in species order. Providing Injection selects explicit net outflow;
    // absence preserves implicit raw outflow without a synthetic subtraction.
    void addDomain(const std::array<double, 4>& areaM2,
        const std::vector<ParticleType>& species,
        const std::vector<std::array<double, 4>>& physicalOutflow,
        const Injection* injection = nullptr);
    // Synchronous borrowed resolver; preserve the map's original order and skip
    // zero/NaN charges before resolving a boundary voltage. No callback retained.
    template<class ResolveVoltage>
    void addEmbedded(const std::unordered_map<int, double>& chargeOutC, ResolveVoltage&& resolve)
    {
        if (!(dtS_ > 0.0)) return;
        for (const auto& entry : chargeOutC)
        {
            if (!(std::abs(entry.second) > 0.0)) continue;
            double voltageV = 0.0;
            if (!resolve(entry.first, voltageV)) continue;
            addCurrent(voltageV, entry.second / dtS_);
        }
    }
    // Explicit evaluates only charged, selected IDs at this time; implicit
    // reads the existing boundary state. Both reuse the same ordered work sum.
    void addEmbedded(const std::unordered_map<int, double>& chargeOutC,
        const PicEmbeddedBC& fallback, const PicEmbeddedBCMap& byId,
        const UnitsContext& units, double timeSI);
    void addEmbedded(const std::unordered_map<int, double>& chargeOutC,
        const PicFieldBoundaryState& boundaries, const UnitsContext& units);
    // Only completed positive-dt work steps advance history; even skipped sides
    // acquire the current field reference. No history is restored from CSV/cache.
    void finishStep();
    const Result& current() const { return result_; }
private:
    void addCurrent(double voltageV, double currentA);
    const bool trackPotentialCurrentSum_;
    std::array<double, 4> normalSI_{}, previousNormalSI_{}, voltageSI_{};
    std::array<bool, 4> havePreviousNormal_{}, dirichlet_{};
    double dtS_ = 0.0;
    Result result_;
};

// A diagnostic count pair; macro arithmetic stays signed, physical arithmetic
// keeps the legacy addition order (including nonfinite inputs).
struct PicSpeciesCount
{
    long long macro;
    double physical;
    // Read the live container, not ParticleType::number. Positive particle
    // weights override the supplied default; preserve the particle sum order.
    static PicSpeciesCount sample(const ParticleType& species, double defaultWeight);
};
std::uint64_t countPicMacroParticles(const std::vector<ParticleType>& species);

// Read-only boundary rates shared by long/wide reports. Callers select the
// original counts, units and time window; integer counts keep their precision.
template<class MacroCount> struct PicBoundaryFlux
{
    MacroCount macro;
    double physical, macroFlux, physicalFlux, current;
    static PicBoundaryFlux sample(MacroCount macro, double physical,
        double charge, double area, double dt)
    {
        const bool valid = area > 0.0 && dt > 0.0;
        return {macro, physical,
            valid ? static_cast<double>(macro) / (area * dt) : 0.0,
            valid ? physical / (area * dt) : 0.0,
            valid ? (charge * physical) / (area * dt) : 0.0};
    }
};
struct PicBoundaryEnergy
{
    PicSpeciesCount out;
    double absorbedEnergyJ, absorbedPowerW, meanEnergyEV, energyFluxWpm2, surfaceAreaM2;
    // dtSI is the caller's positive, clamped step/report interval. A missing
    // surface area retains NaN energy flux, distinct from zero particle flux.
    static PicBoundaryEnergy sample(PicSpeciesCount out, double energyJ, double areaM2, double dtSI);
};
struct PicSpeciesBudget
{
    PicSpeciesCount begin, source, mccSpawn, seeSource, outDomain, outEmbedded;
    PicSpeciesCount end, particleManagement, reactionNet, residual;
};
struct PicSpeciesBudgetInput
{
    PicSpeciesCount begin, source, mccSpawn, seeSource;
    const std::array<std::uint64_t, 4>& domainMacro;
    const std::array<double, 4>& domainPhysical;
    PicSpeciesCount outEmbedded, end, particleManagement;
};
PicSpeciesBudget summarizePicSpeciesBudget(const PicSpeciesBudgetInput& input);

// SI transfer samples shared by the independent series and power intervals.
struct PicEnergyTransfers
{
    double driveJ = 0.0, injectionJ = 0.0, domainAbsorptionJ = 0.0, embeddedAbsorptionJ = 0.0;
    double reactionLossJ = 0.0, seeEmissionJ = 0.0, ionImpactJ = 0.0;
};

// Owns an explicit PIC series interval, independent of the field/power cadence.
// Inputs are borrowed synchronously; results remain owned until the next begin.
// Implicit PIC keeps its original per-step samples and uses the same budget kernel.
class PicDiagnosticWindow
{
public:
    using Energy = PicEnergyTransfers;
    struct Totals
    {
        explicit Totals(std::size_t count);
        Energy energy;
        std::vector<long long> speciesMacroBegin;
        std::vector<double> speciesPhysBegin;
        std::vector<long long> sourceMacroBySpecies;
        std::vector<double> sourcePhysBySpecies;
        std::vector<long long> mccMacroBySpecies;
        std::vector<double> mccPhysBySpecies;
        std::vector<long long> seeMacroBySpecies;
        std::vector<double> seePhysBySpecies;
        std::vector<long long> particleMgmtMacroLoss;
        std::vector<double> particleMgmtPhysLoss;
        std::vector<std::array<std::uint64_t, 4>> wallOutMacroBySpecies;
        std::vector<std::array<double, 4>> wallOutPhysBySpecies;
        std::vector<std::array<double, 4>> wallAbsEnergyJBySpecies;
        std::vector<double> ebOutPhysBySpecies;
        std::vector<long long> ebOutMacroBySpecies;
        std::vector<double> ebAbsEnergyJBySpecies;
        std::unordered_map<std::string, double> reactionEnergyJByChannel;
    };
    struct BoundarySample
    {
        const std::vector<std::array<std::uint64_t, 4>>& domainMacro;
        const std::vector<std::array<double, 4>>& domainPhysical;
        const std::vector<std::array<double, 4>>& domainEnergyJ;
        const std::vector<double>& embeddedPhysical;
        const std::vector<long long>& embeddedMacro;
        const std::vector<double>& embeddedEnergyJ;
        const std::vector<long long>& managementMacro;
        const std::vector<double>& managementPhysical;
    };
    enum class Spawn { Source, Collision };
    explicit PicDiagnosticWindow(std::size_t count) : speciesCount_(count), data_(count) {}
    // A repeated begin while active preserves the first snapshot and time.
    void begin(double timeCode, const std::vector<long long>& macroBegin,
        const std::vector<double>& physicalBegin);
    // Capture live counts only when a new report interval begins. An active
    // window retains its original snapshot without rescanning the particles.
    void beginFromParticles(double timeCode, const std::vector<ParticleType>& species,
        const std::vector<double>& weights);
    void addEnergy(const Energy& energy);
    void addBoundary(const BoundarySample& boundary);
    void addSpawns(Spawn kind, const ParticleSpawnBuffer& pending,
        const std::unordered_map<std::string, std::size_t>& speciesIndexByName,
        const std::vector<double>& macroWeight);
    void addSecondary(std::size_t seeEmittedIndex, const std::vector<Particle>& seeSpawn,
        const std::vector<double>& macroWeight);
    void addReactionEnergy(const std::unordered_map<std::string, double>& energy);
    PicSpeciesBudget budget(std::size_t species, PicSpeciesCount end) const;
    // Closing a report marks the window inactive; clearing waits for next begin.
    void finishReport() { active_ = false; }
    double startTimeCode() const { return startTimeCode_; }
    const Totals& totals() const { return data_; }
private:
    const std::size_t speciesCount_;
    Totals data_;
    bool active_ = false;
    double startTimeCode_ = 0.0;
};

struct PicKineticEnergy
{
    double electronJ = 0.0, ionJ = 0.0, neutralJ = 0.0;
    double total() const { return electronJ + ionJ + neutralJ; }
    // Implicit samples remain inside the caller's charge/energy-observer loop:
    // positive mass/weight, finite energy, and abs(charge)<1e-30 as neutral.
    void addImplicit(const Particle& particle, double mass, double charge, double weight, double velocityScale);
};
struct PicImplicitParticleSummary
{
    double macroCharge = 0.0;
    PicKineticEnergy kinetic;
};
// Synchronous observer at the original point between each particle's charge
// and kinetic-energy accumulation. Arguments are prepared bulk-energy samples;
// no IO dependency, particle mutation, callback retention or RNG consumption.
using PicBulkEnergyObserver = std::function<void(const std::string&, double, double, double)>;
// Explicit bulk histograms resolve fallback weights by the restored name
// index, unlike the implicit summary's positional weights. Borrow everything
// synchronously; keep species/particle order and propagate observer failures.
void observeExplicitPicBulkEnergy(const std::vector<ParticleType>& species,
    const std::unordered_map<std::string, std::size_t>& index,
    const std::vector<double>& weights, double velocityScale, const PicBulkEnergyObserver& observer);
PicImplicitParticleSummary summarizeImplicitPicParticles(const std::vector<ParticleType>& species,
    const std::vector<double>& macroWeight, double velocityScale,
    const PicBulkEnergyObserver& observeEnergy = {});
// Explicit reports sum each species first, classify by charge sign, and keep
// the original weight fallback/nonfinite behavior. No particle state is retained.
PicKineticEnergy summarizeExplicitPicKineticEnergy(const std::vector<ParticleType>& species,
    const std::vector<double>& macroWeight, double velocityScale);

struct PicEnergyBalance
{
    double totalJ = 0.0, absorbedJ = 0.0, residualJ = 0.0;
    double fieldRate = 0.0, kineticRate = 0.0, totalRate = 0.0;
};

// Explicit power history advances only when its original report is sampled.
// The first report has zero derivative/residual; finishReport clears only the
// transfer interval. Code-unit field/time samples retain their legacy meaning.
class PicPowerWindow
{
public:
    struct Report { PicEnergyBalance balance; double potentialCurrentW, potentialEnergyJ; };
    void add(const PicEnergyTransfers& sample);
    Report sample(double time, double fieldEnergy, double kineticEnergy);
    const PicEnergyTransfers& transfers() const { return transfers_; }
    void finishReport() { transfers_ = {}; }
private:
    PicEnergyTransfers transfers_;
    double previousField_ = 0.0, previousTime_ = 0.0, previousTotal_ = 0.0;
    bool havePrevious_ = false;
};

// Implicit energy history advances at every diagnostic step, independently of
// output readiness/cadence. Finite previous field AND kinetic samples are needed;
// both derivatives divide by the current physical dt, not the report interval.
class PicStepEnergyHistory
{
public:
    PicEnergyBalance sample(double dt, double fieldEnergy, double kineticEnergy,
        double driveEnergy, double injectionEnergy, double domainAbsorption, double embeddedAbsorption);
private:
    double previousField_ = std::numeric_limits<double>::quiet_NaN();
    double previousKinetic_ = std::numeric_limits<double>::quiet_NaN();
};

// Explicit current reports use source spawns and species-aggregated outflow at
// report time. Signed charge totals also feed the legacy console at that gate.
struct PicExplicitCurrentBalance
{
    double sourceC = 0.0, domainC = 0.0, embeddedC = 0.0;
    double ionInjectionA = 0.0, electronInjectionA = 0.0, injectionMagnitudeA = 0.0;
    double domainMagnitudeA = 0.0, embeddedMagnitudeA = 0.0;
    double magnitudeResidualA = 0.0, magnitudeRelative = 0.0;
    double injectionSignedA = 0.0, domainSignedA = 0.0, embeddedSignedA = 0.0;
    double signedResidualA = 0.0, signedRelative = 0.0;
};
PicExplicitCurrentBalance summarizeExplicitPicCurrentBalance(
    const ParticleSpawnBuffer& sourceSpawnPending,
    const std::unordered_map<std::string, std::size_t>& speciesIndexByName,
    const std::vector<ParticleType>& species, const std::vector<double>& macroWeight,
    const std::vector<std::array<double, 4>>& wallOutPhysBySpecies,
    const std::vector<double>& ebOutPhysBySpecies, double dtSI);

// Implicit per-step particle/source counts, injection energy and charge closure,
// independent of CSV cadence. Stage and physical charge closure remain distinct.
class PicStepParticleBalance
{
public:
    struct Result
    {
        double macroBeginC = 0.0, macroEndC = std::numeric_limits<double>::quiet_NaN();
        double stageMacroC = 0.0, stageGridC = std::numeric_limits<double>::quiet_NaN();
        double stageResidualC = std::numeric_limits<double>::quiet_NaN();
        double stageRelative = std::numeric_limits<double>::quiet_NaN();
        double sourceC = 0.0, secondaryC = 0.0, domainC = 0.0, embeddedC = 0.0;
        double expectedEndC = std::numeric_limits<double>::quiet_NaN();
        double residualC = std::numeric_limits<double>::quiet_NaN();
        double relative = std::numeric_limits<double>::quiet_NaN();
    };
    explicit PicStepParticleBalance(double beginC, std::size_t speciesCount = 0);
    PicStepParticleBalance(const std::vector<ParticleType>& species, const std::vector<double>& weights);
    const PicSpeciesCount& begin(std::size_t species) const { return begin_.at(species); }
    enum class Spawn { Collision, Secondary, Source };
    struct SpawnTotals { std::vector<PicSpeciesCount> species; double energyJ = 0.; };
    const SpawnTotals& spawns(Spawn kind) const { return spawns_.at(static_cast<std::size_t>(kind)); }
    // Count each retained buffer at its original integration stage. Two MCC
    // halves add directly to the same per-species sums; MCC contributes neither
    // physical source charge nor injected kinetic energy here.
    void addSpawns(Spawn kind, const ParticleSpawnBuffer& pending,
        const std::vector<ParticleType>& species, const std::vector<double>& weights,
        const std::unordered_map<std::string,std::size_t>& index, double velocityScale);
    double injectionEnergyJ() const { return spawns(Spawn::Source).energyJ + spawns(Spawn::Secondary).energyJ; }
    // Keep the caller's particle order and OpenMP reduction grouping unchanged.
    void addStageMacro(double chargeC) { result_.stageMacroC += chargeC; }
    void finishStage(double gridC);
    void addDomain(double chargeC) { result_.domainC += chargeC; }
    void addEmbedded(double chargeC) { result_.embeddedC += chargeC; }
    void addSecondary(double chargeC) { result_.secondaryC += chargeC; }
    // Add the SEE subtotal once, after its original merge, before source spawns.
    void includeSecondarySource() { result_.sourceC += result_.secondaryC; }
    void addSource(double chargeC) { result_.sourceC += chargeC; }
    void finish(double endC);
    const Result& current() const { return result_; }
private:
    Result result_;
    std::vector<PicSpeciesCount> begin_;
    std::array<SpawnTotals,3> spawns_;
};
double sumPicMacroCharge(const std::vector<ParticleType>& species,
    const std::vector<double>& macroWeight);

// Timers are sampled by the integrators at their original stage boundaries.
struct PicExplicitStepTiming
{
    double buildWeights_s = 0.0;
    double deposit_s = 0.0;
    double rhs_s = 0.0;
    double evalE_s = 0.0;
    double amr_indicator_s = 0.0;
};

struct PicImplicitStepTiming
{
    double mcc_pre_s = 0.0;
    double buildWeightsN_s = 0.0;
    double gatherEn_s = 0.0;
    double predictorPrep_s = 0.0;
    double predictorSegment_s = 0.0;
    std::uint64_t predictor_fast_exit_hits = 0u;
    std::uint64_t predictor_bisection_fallback_hits = 0u;
    std::uint64_t predictor_coarsen_hits = 0u;
    std::uint64_t predictor_coarsen_parent_climb = 0u;
    double orbitAvgWeights_s = 0.0;
    int orbitAvg_fused_mode = 0;
    double orbitAvg_fused_eval_s = 0.0;
    double orbitAvg_fused_accum_s = 0.0;
    double orbitAvg_sampleFill_s = 0.0;
    double orbitAvg_sampleCacheBuild_s = 0.0;
    double orbitAvg_localReduce_s = 0.0;
    double orbitAvg_concat_s = 0.0;
    std::uint64_t buildWeightsN_cwFallbackCount = 0u;
    std::uint64_t buildWeightsN_particleCount = 0u;
    std::uint64_t predictor_segments_total = 0u;
    std::uint64_t orbitAvg_sampleCwFallbackCount = 0u;
    std::uint64_t orbitAvg_sampleCount = 0u;
    std::uint64_t orbitAvg_sampleCache_nnz = 0u;
    std::uint64_t orbitAvg_out_nnz = 0u;
    std::uint64_t orbitAvg_adaptAcceptCount = 0u;
    std::uint64_t orbitAvg_adaptRefineCount = 0u;
    std::uint64_t orbitAvg_adaptHitMaxDepthCount = 0u;
    std::uint64_t orbitAvg_adaptAcceptedDepthSum = 0u;
    double depositOrbit_s = 0.0;
    std::uint64_t depositOrbit_qp_size_total = 0u;
    std::uint64_t depositOrbit_nnz_total = 0u;
    double assemblePolarization_s = 0.0;
    double polPoisson_assemble_s = 0.0;
    double polPoisson_factor_s = 0.0;
    double polPoisson_solve_s = 0.0;
    double polPoisson_total_s = 0.0;
    int pol_cg_iterations = 0;
    double pol_cg_residual = 0.0;
    double evalE_s = 0.0;
    double waveformPeriodicSample_s = 0.0;
    double profiles_s = 0.0;
    double gatherPush_s = 0.0;
    double fluxWaveform_s = 0.0;
    double mcc_post_s = 0.0;
    double source_s = 0.0;
    double pm_s = 0.0;
    double diagnostics_s = 0.0;
    double packedFieldCompute_s = 0.0;
    double packedFieldIo_s = 0.0;
    double packedFieldOutput_s = 0.0;
    double amr_indicator_s = 0.0;
    double amr_adapt_s = 0.0;
    double amr_buildTargetLevels_s = 0.0;
    double amr_faceGraph_s = 0.0;
    double amr_cutMetrics_s = 0.0;
    double amr_reinitLevelset_s = 0.0;
    double amr_geom_initializeLevelSets_s = 0.0;
    double amr_geom_pdeReinitLevelset_s = 0.0;
    double amr_geom_refreshDerived_s = 0.0;
    double amr_geom_checkCutInfo_s = 0.0;
    double amr_geom_rebuildComposite_s = 0.0;
    double amr_geom_markPhases_s = 0.0;
    double amr_fieldRebuild_s = 0.0;
    double amr_s = 0.0;
    double total_s = 0.0;
};

struct PicImplicitTimingSummary
{
    struct Stage { const char* name; double t; };
    double buildFallbackFraction, orbitFallbackFraction, totalFallbackFraction;
    std::uint64_t totalFallbackCount;
    std::array<Stage, 3> top;
};
PicImplicitTimingSummary summarizePicImplicitTiming(const PicImplicitStepTiming& timing);
