#pragma once
#include "Modules/Plasma/Diagnostics/PhysicsDiagnostics.h"
#include "Modules/Plasma/ExplicitPIC/PicConfig.h"
#include "Modules/Plasma/Physics/MCCCollisionOperator.h"
#include <filesystem>
#include <fstream>

namespace PhysicsDiagnostics
{
ElectrodeCurrentGroups makeElectrodeCurrentGroups(
    const std::vector<OutputConfig::DiagnosticsConfig::ElectrodeGroup>& config);

// Own output streams and sampled rows. isOpen requires an actually open, healthy
// file; default/disabled streams cannot activate diagnostic sampling or lazy IO.
class ElectrodeCurrentWriter
{
public:
    bool open(const std::filesystem::path& rootDir,
              const std::vector<ElectrodeCurrentGroup>& groups,
              bool writePerBoundary,
              int startStep = 0);

    bool isOpen() const;

    void beginStep();

    void record(int boundaryId, double chargeC, double neutralWeight);

    void writeStep(int step, double timeS, double dtS);

private:
    bool writePerBoundary_ = false;
    std::ofstream groupOut_;
    std::ofstream boundaryOut_;
    ElectrodeCurrentGroups groups_;
    std::unordered_map<int, ElectrodeCurrentBucket> boundaryBuckets_;
};

class ReactionDiagnosticsWriter
{
public:
    bool open(const std::filesystem::path& rootDir,
        const OutputConfig::DiagnosticsConfig::ReactionStatsConfig& config,
        const ReactionDatabase* reactions, const MCCCollisionOperator* operation, int startStep = 0);
    struct ChannelColumn
    {
        std::string key;
        std::string token;
    };

    bool open(const std::filesystem::path& rootDir,
              bool summaryEnabled,
              bool perChannelEnabled,
              bool perProjectileEnabled,
              bool spawnedSpeciesEnabled,
              const std::vector<std::string>& projectileNames,
              const std::vector<std::string>& channelKeys,
              const std::vector<std::string>& spawnedSpeciesNames,
              int startStep = 0);

    bool isOpen() const;

    void writeStep(int step, double timeS, double dtS, const MCCCollisionOperator::StepStats& stats);

private:
    std::ofstream summaryOut_;
    std::ofstream channelOut_;
    std::ofstream energyOut_;
    std::ofstream projectileOut_;
    std::ofstream spawnOut_;
    std::vector<std::string> projectileNames_;
    std::vector<ChannelColumn> channelColumns_;
    std::vector<std::string> spawnedSpeciesNames_;
};

class EnergyHistogramWriter
{
public:
    bool open(const std::filesystem::path& rootDir,
        const OutputConfig::DiagnosticsConfig::EnergyHistogramConfig& config,
        const std::vector<SpeciesConfig>& species, int startStep = 0);
    bool open(const std::filesystem::path& rootDir,
              int bins,
              double maxEnergyEV,
              const std::vector<EnergyHistogramSample>& samples,
              int startStep = 0);

    bool isOpen() const;

    void beginStep();

    void observe(const std::string& speciesName,
                 double chargeC,
                 const std::string& location,
                 double energyEV,
                 double weight);

    void writeStep(int step, double timeS);

private:
    static bool matchesSpecies(const std::string& selectorRaw,
                               const std::string& speciesLower,
                               double chargeC);

    static bool matchesLocation(const std::string& selectorRaw, const std::string& locationLower);

    int bins_ = 0;
    double maxEnergyEV_ = 0.0;
    std::ofstream out_;
    std::vector<EnergyHistogramSample> samples_;
    std::vector<std::vector<double>> hist_;
};

// Own the configured diagnostic streams, cadence and per-step sampling gates.
// Console electrode buckets are updated every step, independently of file IO.
// Inputs are copied; reaction/operator pointers are borrowed only during open.
class PicDiagnosticsOutput
{
public:
    enum class Mode { Explicit, Implicit };
    PicDiagnosticsOutput(Mode mode, const OutputConfig::DiagnosticsConfig& config,
        const std::filesystem::path& root, int baseEvery, int resumeStep,
        const std::vector<SpeciesConfig>& species, const ReactionDatabase* reactions = nullptr,
        const MCCCollisionOperator* operation = nullptr);
    PicDiagnosticsOutput(const PicDiagnosticsOutput&) = delete;
    PicDiagnosticsOutput& operator=(const PicDiagnosticsOutput&) = delete;

    // Explicit initialization/retry stays at its original caller phase. Implicit
    // reaction output is attempted once during construction, even without MCC.
    void prepareExplicitReactions(const ReactionDatabase* reactions, const MCCCollisionOperator* operation);
    void beginStep(int completedStep, const ReactionDatabase* reactions = nullptr,
        const MCCCollisionOperator* operation = nullptr);
    bool sampleElectrodes() const { return sampleElectrodes_; }
    bool sampleEnergy() const { return sampleEnergy_; }
    bool includesDomainElectrodes() const { return config_.electrodeCurrent.includeDomain; }
    const ElectrodeCurrentGroups& consoleElectrodes() const { return consoleElectrodes_; }
    void recordElectrode(int boundaryId, double chargeC, double neutralWeight);
    void observeEnergy(const std::string& species, double chargeC, const std::string& location,
        double energyEV, double weight);
    // Separate publication points retain reaction -> bulk sampling -> energy ->
    // electrode order in explicit PIC, and already-sampled bulk in implicit PIC.
    void writeReactions(int step, double timeS, double dtS, const MCCCollisionOperator::StepStats& stats);
    void writeEnergy(int step, double timeS);
    void writeElectrodes(int step, double timeS, double dtS);
private:
    void openReactions(const ReactionDatabase* reactions, const MCCCollisionOperator* operation);
    const Mode mode_;
    const OutputConfig::DiagnosticsConfig config_;
    const std::filesystem::path root_;
    const int resumeStep_, electrodeEvery_, reactionEvery_, energyEvery_;
    ElectrodeCurrentGroups consoleElectrodes_;
    ElectrodeCurrentWriter electrodes_;
    ReactionDiagnosticsWriter reactions_;
    EnergyHistogramWriter energy_;
    bool sampleElectrodes_ = false, sampleReactions_ = false, sampleEnergy_ = false;
};

} // namespace PhysicsDiagnostics
