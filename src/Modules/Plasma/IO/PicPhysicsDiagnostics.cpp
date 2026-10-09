#include "Modules/Plasma/IO/PicPhysicsDiagnostics.h"
#include "Modules/Plasma/Diagnostics/TimeSeriesIO.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <unordered_set>

namespace PhysicsDiagnostics
{
namespace
{
bool ready(const std::ofstream& stream) { return stream.is_open() && static_cast<bool>(stream); }
}
bool ElectrodeCurrentWriter::open(const std::filesystem::path& rootDir,
          const std::vector<ElectrodeCurrentGroup>& groups,
          bool writePerBoundary,
          int startStep)
{
    groups_ = ElectrodeCurrentGroups(groups);

    const std::filesystem::path dir = rootDir / "electrodes";
    std::filesystem::create_directories(dir);

    const bool groupHasExisting = TimeSeriesIO::openPreparedStepCsv(groupOut_, dir / "current_groups_step.csv", startStep);
    if (!ready(groupOut_))
        return false;
    if (!groupHasExisting)
        groupOut_ << "step,time_s,group,boundary_count,charge_C,current_A,electron_charge_C,electron_current_A,ion_charge_C,ion_current_A,neutral_count_phys\n";
    groupOut_ << std::setprecision(16) << std::fixed;

    writePerBoundary_ = writePerBoundary;
    if (writePerBoundary_)
    {
        const bool boundaryHasExisting = TimeSeriesIO::openPreparedStepCsv(boundaryOut_, dir / "current_boundaries_step.csv", startStep);
        if (!ready(boundaryOut_))
            return false;
        if (!boundaryHasExisting)
            boundaryOut_ << "step,time_s,boundary_id,group,charge_C,current_A,electron_charge_C,electron_current_A,ion_charge_C,ion_current_A,neutral_count_phys\n";
        boundaryOut_ << std::setprecision(16) << std::fixed;
    }
    return true;
}

bool ElectrodeCurrentWriter::isOpen() const
{
    return ready(groupOut_);
}

void ElectrodeCurrentWriter::beginStep()
{
    groups_.beginStep();
    boundaryBuckets_.clear();
}

void ElectrodeCurrentWriter::record(int boundaryId, double chargeC, double neutralWeight)
{
    if (boundaryId <= 0)
        return;
    groups_.record(boundaryId, chargeC, neutralWeight);
    if (writePerBoundary_)
        boundaryBuckets_[boundaryId].add(chargeC, neutralWeight);
}

void ElectrodeCurrentWriter::writeStep(int step, double timeS, double dtS)
{
    if (!ready(groupOut_))
        return;
    const double invDt = (dtS > 0.0) ? (1.0 / dtS) : 0.0;
    for (std::size_t gi = 0; gi < groups_.groups().size(); ++gi)
    {
        const auto& bucket = groups_.buckets()[gi];
        groupOut_ << step << ','
                  << timeS << ','
                  << groups_.groups()[gi].name << ','
                  << groups_.groups()[gi].boundaryIds.size() << ','
                  << bucket.chargeC << ','
                  << (bucket.chargeC * invDt) << ','
                  << bucket.electronChargeC << ','
                  << (bucket.electronChargeC * invDt) << ','
                  << bucket.ionChargeC << ','
                  << (bucket.ionChargeC * invDt) << ','
                  << bucket.neutralCountPhys << '\n';
    }
    groupOut_.flush();

    if (!writePerBoundary_ || !ready(boundaryOut_))
        return;

    for (const auto& kv : boundaryBuckets_)
    {
        const int boundaryId = kv.first;
        const auto& bucket = kv.second;
        std::string groupName = "ungrouped";
        const auto index = groups_.indexOf(boundaryId);
        if (index < groups_.groups().size()) groupName = groups_.groups()[index].name;
        boundaryOut_ << step << ','
                     << timeS << ','
                     << boundaryId << ','
                     << groupName << ','
                     << bucket.chargeC << ','
                     << (bucket.chargeC * invDt) << ','
                     << bucket.electronChargeC << ','
                     << (bucket.electronChargeC * invDt) << ','
                     << bucket.ionChargeC << ','
                     << (bucket.ionChargeC * invDt) << ','
                     << bucket.neutralCountPhys << '\n';
    }
    boundaryOut_.flush();
}

bool ReactionDiagnosticsWriter::open(const std::filesystem::path& rootDir,
          bool summaryEnabled,
          bool perChannelEnabled,
          bool perProjectileEnabled,
          bool spawnedSpeciesEnabled,
          const std::vector<std::string>& projectileNames,
          const std::vector<std::string>& channelKeys,
          const std::vector<std::string>& spawnedSpeciesNames,
          int startStep)
{
    projectileNames_ = projectileNames;
    std::sort(projectileNames_.begin(), projectileNames_.end());

    channelColumns_.clear();
    channelColumns_.reserve(channelKeys.size());
    for (const auto& key : channelKeys)
    {
        ChannelColumn column;
        column.key = key;
        const std::size_t pos = key.rfind("::");
        if (pos == std::string::npos)
            column.token = sanitizeToken(key);
        else
            column.token = sanitizeToken(key.substr(0, pos)) + "_" + sanitizeToken(key.substr(pos + 2));
        channelColumns_.push_back(std::move(column));
    }
    std::sort(channelColumns_.begin(), channelColumns_.end(), [](const ChannelColumn& a, const ChannelColumn& b) {
        return a.key < b.key;
    });

    spawnedSpeciesNames_ = spawnedSpeciesNames;
    std::sort(spawnedSpeciesNames_.begin(), spawnedSpeciesNames_.end());

    const std::filesystem::path dir = rootDir / "reactions";
    std::filesystem::create_directories(dir);

    auto openCsvWithHeader = [&](std::ofstream& out,
                                 const std::filesystem::path& path,
                                 const std::string& header) -> bool {
        bool compatibleExisting = false;
        if (startStep > 0)
        {
            std::error_code ec;
            if (std::filesystem::exists(path, ec))
            {
                std::ifstream in(path, std::ios::in);
                std::string firstLine;
                if (in && std::getline(in, firstLine) && TimeSeriesIO::trimCopy(firstLine) == header)
                    compatibleExisting = true;
            }
        }

        const bool hasExisting = TimeSeriesIO::openPreparedStepCsv(out, path, compatibleExisting ? startStep : 0);
        if (!out)
            return false;
        if (!hasExisting)
            out << header << '\n';
        out << std::setprecision(16) << std::fixed;
        out.flush();
        return true;
    };

    if (summaryEnabled)
    {
        const std::string header = "step,time_s,macro_particles,phys_particles,macro_attempts,phys_attempts,macro_null,phys_null,macro_collisions,phys_collisions,collision_frequency_1ps,macro_spawned,phys_spawned,spawn_frequency_1ps,energy_loss_J_total,max_nu_dt_majorant";
        if (!openCsvWithHeader(summaryOut_, dir / "reaction_step.csv", header))
            return false;
    }

    if (perChannelEnabled)
    {
        std::string channelHeader = "step,time_s";
        for (const auto& column : channelColumns_)
        {
            channelHeader += "," + column.token + "_macro_collisions";
            channelHeader += "," + column.token + "_phys_collisions";
            channelHeader += "," + column.token + "_collision_frequency_1ps";
        }
        if (!openCsvWithHeader(channelOut_, dir / "reaction_channel_step.csv", channelHeader))
            return false;

        std::string energyHeader = "step,time_s";
        for (const auto& column : channelColumns_)
        {
            energyHeader += "," + column.token + "_macro_collisions";
            energyHeader += "," + column.token + "_phys_collisions";
            energyHeader += "," + column.token + "_projectile_energy_loss_J_total";
            energyHeader += "," + column.token + "_projectile_energy_loss_power_W";
            energyHeader += "," + column.token + "_projectile_energy_loss_fraction_step";
        }
        if (!openCsvWithHeader(energyOut_, dir / "reaction_energy_step.csv", energyHeader))
            return false;
    }

    if (perProjectileEnabled)
    {
        std::string projectileHeader = "step,time_s";
        for (const auto& projectile : projectileNames_)
        {
            const std::string token = sanitizeToken(projectile);
            projectileHeader += "," + token + "_macro_particles";
            projectileHeader += "," + token + "_phys_particles";
            projectileHeader += "," + token + "_macro_attempts";
            projectileHeader += "," + token + "_phys_attempts";
            projectileHeader += "," + token + "_macro_null";
            projectileHeader += "," + token + "_phys_null";
            projectileHeader += "," + token + "_macro_collisions";
            projectileHeader += "," + token + "_phys_collisions";
            projectileHeader += "," + token + "_collision_frequency_1ps";
            projectileHeader += "," + token + "_macro_spawned";
            projectileHeader += "," + token + "_phys_spawned";
            projectileHeader += "," + token + "_spawn_frequency_1ps";
            projectileHeader += "," + token + "_energy_loss_J_total";
            projectileHeader += "," + token + "_max_nu_dt_majorant";
        }
        if (!openCsvWithHeader(projectileOut_, dir / "reaction_projectile_step.csv", projectileHeader))
            return false;
    }

    if (spawnedSpeciesEnabled)
    {
        std::string spawnHeader = "step,time_s";
        for (const auto& speciesName : spawnedSpeciesNames_)
        {
            const std::string token = sanitizeToken(speciesName);
            spawnHeader += "," + token + "_macro_spawned";
            spawnHeader += "," + token + "_phys_spawned";
            spawnHeader += "," + token + "_spawn_frequency_1ps";
        }
        if (!openCsvWithHeader(spawnOut_, dir / "reaction_spawn_step.csv", spawnHeader))
            return false;
    }

    return true;
}

bool ReactionDiagnosticsWriter::isOpen() const
{
    return ready(summaryOut_) || ready(channelOut_) || ready(energyOut_) || ready(projectileOut_) || ready(spawnOut_);
}

void ReactionDiagnosticsWriter::writeStep(int step, double timeS, double dtS, const MCCCollisionOperator::StepStats& stats)
{
    const double invDt = (dtS > 0.0) ? (1.0 / dtS) : 0.0;
    if (ready(summaryOut_))
    {
        summaryOut_ << step << ','
                    << timeS << ','
                    << stats.n_particles << ','
                    << stats.n_particles_weighted << ','
                    << stats.n_attempts << ','
                    << stats.n_attempts_weighted << ','
                    << stats.n_null << ','
                    << stats.n_null_weighted << ','
                    << stats.n_collisions << ','
                    << stats.n_collisions_weighted << ','
                    << (stats.n_collisions_weighted * invDt) << ','
                    << stats.n_spawned << ','
                    << stats.n_spawned_weighted << ','
                    << (stats.n_spawned_weighted * invDt) << ','
                    << stats.energy_loss_J_total << ','
                    << stats.max_nu_dt_majorant << '\n';
        summaryOut_.flush();
    }

    if (ready(projectileOut_))
    {
        projectileOut_ << step << ',' << timeS;
        for (const auto& projectile : projectileNames_)
        {
            const auto it = stats.per_projectile.find(projectile);
            const MCCCollisionOperator::StepStats::ProjectileStats* p =
                (it != stats.per_projectile.end()) ? &it->second : nullptr;
            projectileOut_ << ','
                           << (p ? p->n_particles : 0ull) << ','
                           << (p ? p->n_particles_weighted : 0.0) << ','
                           << (p ? p->n_attempts : 0ull) << ','
                           << (p ? p->n_attempts_weighted : 0.0) << ','
                           << (p ? p->n_null : 0ull) << ','
                           << (p ? p->n_null_weighted : 0.0) << ','
                           << (p ? p->n_collisions : 0ull) << ','
                           << (p ? p->n_collisions_weighted : 0.0) << ','
                           << ((p ? p->n_collisions_weighted : 0.0) * invDt) << ','
                           << (p ? p->n_spawned : 0ull) << ','
                           << (p ? p->n_spawned_weighted : 0.0) << ','
                           << ((p ? p->n_spawned_weighted : 0.0) * invDt) << ','
                           << (p ? p->energy_loss_J_total : 0.0) << ','
                           << (p ? p->max_nu_dt_majorant : 0.0);
        }
        projectileOut_ << '\n';
        projectileOut_.flush();
    }

    if (ready(channelOut_))
    {
        channelOut_ << step << ',' << timeS;
        for (const auto& column : channelColumns_)
        {
            const auto itMacro = stats.per_channel.find(column.key);
            const auto itPhys = stats.per_channel_weighted.find(column.key);
            const std::uint64_t macro = (itMacro != stats.per_channel.end()) ? itMacro->second : 0ull;
            const double phys = (itPhys != stats.per_channel_weighted.end()) ? itPhys->second : 0.0;
            channelOut_ << ','
                        << macro << ','
                        << phys << ','
                        << (phys * invDt);
        }
        channelOut_ << '\n';
        channelOut_.flush();
    }

    if (ready(energyOut_))
    {
        const double totalEnergyLossJ = std::max(0.0, stats.energy_loss_J_total);
        energyOut_ << step << ',' << timeS;
        for (const auto& column : channelColumns_)
        {
            const auto itMacro = stats.per_channel.find(column.key);
            const auto itPhys = stats.per_channel_weighted.find(column.key);
            const auto itEnergy = stats.per_channel_energy_loss_J.find(column.key);
            const std::uint64_t macroCollisions = (itMacro != stats.per_channel.end()) ? itMacro->second : 0ull;
            const double physCollisions = (itPhys != stats.per_channel_weighted.end()) ? itPhys->second : 0.0;
            const double energyLossJ = (itEnergy != stats.per_channel_energy_loss_J.end()) ? itEnergy->second : 0.0;
            const double powerW = energyLossJ * invDt;
            const double fraction = (totalEnergyLossJ > 0.0) ? (energyLossJ / totalEnergyLossJ) : 0.0;
            energyOut_ << ','
                       << macroCollisions << ','
                       << physCollisions << ','
                       << energyLossJ << ','
                       << powerW << ','
                       << fraction;
        }
        energyOut_ << '\n';
        energyOut_.flush();
    }

    if (ready(spawnOut_))
    {
        spawnOut_ << step << ',' << timeS;
        for (const auto& speciesName : spawnedSpeciesNames_)
        {
            const auto itMacro = stats.spawn_per_species.find(speciesName);
            const auto itPhys = stats.spawn_per_species_weighted.find(speciesName);
            const std::uint64_t macro = (itMacro != stats.spawn_per_species.end()) ? itMacro->second : 0ull;
            const double phys = (itPhys != stats.spawn_per_species_weighted.end()) ? itPhys->second : 0.0;
            spawnOut_ << ','
                      << macro << ','
                      << phys << ','
                      << (phys * invDt);
        }
        spawnOut_ << '\n';
        spawnOut_.flush();
    }
}

bool EnergyHistogramWriter::open(const std::filesystem::path& rootDir,
          int bins,
          double maxEnergyEV,
          const std::vector<EnergyHistogramSample>& samples,
          int startStep)
{
    bins_ = std::max(1, bins);
    maxEnergyEV_ = (maxEnergyEV > 0.0) ? maxEnergyEV : 200.0;
    samples_ = samples;
    hist_.assign(samples_.size(), std::vector<double>(static_cast<std::size_t>(bins_), 0.0));

    const std::filesystem::path dir = rootDir / "energy";
    std::filesystem::create_directories(dir);
    const bool hasExisting = TimeSeriesIO::openPreparedStepCsv(out_, dir / "histogram_step.csv", startStep);
    if (!ready(out_))
        return false;
    if (!hasExisting)
        out_ << "step,time_s,sample,species_selector,location,bin_index,bin_lo_eV,bin_hi_eV,phys_count\n";
    out_ << std::setprecision(16) << std::fixed;
    return true;
}

bool EnergyHistogramWriter::isOpen() const
{
    return ready(out_);
}

void EnergyHistogramWriter::beginStep()
{
    for (auto& row : hist_)
        std::fill(row.begin(), row.end(), 0.0);
}

void EnergyHistogramWriter::observe(const std::string& speciesName,
             double chargeC,
             const std::string& location,
             double energyEV,
             double weight)
{
    if (!ready(out_) || !(energyEV >= 0.0) || !(weight > 0.0) || bins_ <= 0)
        return;
    const std::string speciesLower = lowerCopy(speciesName);
    const std::string locationLower = lowerCopy(location);
    int bin = static_cast<int>(std::floor((energyEV / maxEnergyEV_) * static_cast<double>(bins_)));
    if (bin < 0)
        bin = 0;
    if (bin >= bins_)
        bin = bins_ - 1;

    for (std::size_t si = 0; si < samples_.size(); ++si)
    {
        if (!matchesSpecies(samples_[si].species, speciesLower, chargeC))
            continue;
        if (!matchesLocation(samples_[si].location, locationLower))
            continue;
        hist_[si][static_cast<std::size_t>(bin)] += weight;
    }
}

void EnergyHistogramWriter::writeStep(int step, double timeS)
{
    if (!ready(out_))
        return;
    const double dE = maxEnergyEV_ / static_cast<double>(bins_);
    for (std::size_t si = 0; si < samples_.size(); ++si)
    {
        for (int bi = 0; bi < bins_; ++bi)
        {
            const double lo = dE * static_cast<double>(bi);
            const double hi = lo + dE;
            out_ << step << ','
                 << timeS << ','
                 << samples_[si].name << ','
                 << samples_[si].species << ','
                 << samples_[si].location << ','
                 << bi << ','
                 << lo << ','
                 << hi << ','
                 << hist_[si][static_cast<std::size_t>(bi)] << '\n';
        }
    }
    out_.flush();
}

bool EnergyHistogramWriter::matchesSpecies(const std::string& selectorRaw,
                           const std::string& speciesLower,
                           double chargeC)
{
    const std::string selector = lowerCopy(selectorRaw);
    if (selector.empty() || selector == "all" || selector == "*")
        return true;
    if (selector == "electron" || selector == "electrons")
        return chargeC < 0.0;
    if (selector == "ion" || selector == "ions")
        return chargeC > 0.0;
    if (selector == "neutral" || selector == "neutrals")
        return chargeC == 0.0;
    return selector == speciesLower;
}

bool EnergyHistogramWriter::matchesLocation(const std::string& selectorRaw, const std::string& locationLower)
{
    const std::string selector = lowerCopy(selectorRaw);
    if (selector.empty() || selector == "all")
        return true;
    if (selector == "wall")
        return locationLower == "domain" || locationLower == "embedded" || locationLower == "wall";
    return selector == locationLower;
}

ElectrodeCurrentGroups makeElectrodeCurrentGroups(
    const std::vector<OutputConfig::DiagnosticsConfig::ElectrodeGroup>& config)
{
    std::vector<ElectrodeCurrentGroup> electrodeGroups;
    electrodeGroups.reserve(config.size());
    for (const auto& groupCfg : config)
    {
        ElectrodeCurrentGroup group;
        group.name = groupCfg.name;
        group.boundaryIds = groupCfg.boundaryIds;
        electrodeGroups.push_back(std::move(group));
    }
    return ElectrodeCurrentGroups(std::move(electrodeGroups));
}

bool ReactionDiagnosticsWriter::open(const std::filesystem::path& rootDir,
    const OutputConfig::DiagnosticsConfig::ReactionStatsConfig& config,
    const ReactionDatabase* reactions, const MCCCollisionOperator* operation, int startStep)
{
    std::vector<std::string> projectiles, channelKeys, spawnedSpecies;
    std::unordered_set<std::string> projectileSet;
    std::unordered_set<std::string> channelKeySet;
    std::unordered_set<std::string> spawnedSpeciesSet;

    if (reactions)
    {
        for (const auto& key : reactions->Keys())
        {
            const ReactionSet* reactionSet = reactions->Get(key);
            if (!reactionSet)
                continue;
            if (!reactionSet->projectile.empty())
                projectileSet.insert(reactionSet->projectile);
            for (const auto& channel : reactionSet->channels)
            {
                if (!reactionSet->projectile.empty() && !channel.name.empty())
                    channelKeySet.insert(reactionSet->projectile + "::" + channel.name);
                for (const auto& speciesName : channel.products)
                {
                    if (!speciesName.empty())
                        spawnedSpeciesSet.insert(speciesName);
                }
            }
        }
    }

    if (projectileSet.empty() && operation)
    {
        for (const auto& name : operation->ProjectileNames())
            projectileSet.insert(name);
        for (const auto& key : operation->ChannelKeys())
            channelKeySet.insert(key);
        for (const auto& name : operation->SpawnSpeciesNames())
            spawnedSpeciesSet.insert(name);
    }

    projectiles.assign(projectileSet.begin(), projectileSet.end());
    channelKeys.assign(channelKeySet.begin(), channelKeySet.end());
    spawnedSpecies.assign(spawnedSpeciesSet.begin(), spawnedSpeciesSet.end());
    return open(rootDir, config.summary, config.perChannel, config.perProjectile,
        config.spawnedSpecies, projectiles, channelKeys, spawnedSpecies, startStep);
}

bool EnergyHistogramWriter::open(const std::filesystem::path& rootDir,
    const OutputConfig::DiagnosticsConfig::EnergyHistogramConfig& config,
    const std::vector<SpeciesConfig>& species, int startStep)
{
    std::vector<EnergyHistogramSample> energySamples;
    energySamples.reserve(config.samples.size());
    for (const auto& sampleCfg : config.samples)
    {
        EnergyHistogramSample sample;
        sample.name = sampleCfg.name;
        sample.species = sampleCfg.species;
        sample.location = sampleCfg.location;
        energySamples.push_back(std::move(sample));
    }
    if (energySamples.empty())
    {
        std::vector<std::string> trackedSpeciesNames;
        trackedSpeciesNames.reserve(species.size());
        for (const auto& spCfg : species)
            trackedSpeciesNames.push_back(spCfg.name);
        energySamples = defaultBulkSamples(trackedSpeciesNames);
    }
    return open(rootDir, config.bins, config.maxEnergyEV, energySamples, startStep);
}
PicDiagnosticsOutput::PicDiagnosticsOutput(Mode mode, const OutputConfig::DiagnosticsConfig& config,
    const std::filesystem::path& root, int baseEvery, int resumeStep,
    const std::vector<SpeciesConfig>& species, const ReactionDatabase* reactions,
    const MCCCollisionOperator* operation)
    : mode_(mode), config_(config), root_(root), resumeStep_(resumeStep),
      electrodeEvery_(resolveEvery(baseEvery, config.electrodeCurrent.every)),
      reactionEvery_(resolveEvery(baseEvery, config.reactions.every)),
      energyEvery_(resolveEvery(baseEvery, config.energyHistogram.every)),
      consoleElectrodes_(makeElectrodeCurrentGroups(config.electrodeCurrent.groups))
{
    const char* label = mode_ == Mode::Explicit ? "ExplicitPIC" : "ImplicitPIC";
    if (config_.electrodeCurrent.enabled &&
        !electrodes_.open(root_, consoleElectrodes_.groups(), config_.electrodeCurrent.perBoundary, resumeStep_))
        std::cerr << '[' << label << "] WARNING: failed to create electrode diagnostics under "
                  << root_.string() << std::endl;
    if (mode_ == Mode::Implicit && config_.reactions.enabled) openReactions(reactions, operation);
    if (config_.energyHistogram.enabled && !energy_.open(root_, config_.energyHistogram, species, resumeStep_))
        std::cerr << '[' << label << "] WARNING: failed to create energy histogram diagnostics under "
                  << root_.string() << std::endl;
}

void PicDiagnosticsOutput::openReactions(const ReactionDatabase* reactions, const MCCCollisionOperator* operation)
{
    if (!reactions_.open(root_, config_.reactions, reactions, operation, resumeStep_))
        std::cerr << '[' << (mode_ == Mode::Explicit ? "ExplicitPIC" : "ImplicitPIC")
                  << "] WARNING: failed to create reaction diagnostics under " << root_.string() << std::endl;
}

void PicDiagnosticsOutput::prepareExplicitReactions(const ReactionDatabase* reactions, const MCCCollisionOperator* operation)
{
    if (mode_ == Mode::Explicit && config_.reactions.enabled && !reactions_.isOpen() && operation)
        openReactions(reactions, operation);
}

void PicDiagnosticsOutput::beginStep(int completedStep, const ReactionDatabase* reactions,
    const MCCCollisionOperator* operation)
{
    prepareExplicitReactions(reactions, operation);
    sampleElectrodes_ = electrodes_.isOpen() && electrodeEvery_ > 0 && completedStep % electrodeEvery_ == 0;
    sampleReactions_ = reactions_.isOpen() && reactionEvery_ > 0 && completedStep % reactionEvery_ == 0;
    sampleEnergy_ = energy_.isOpen() && energyEvery_ > 0 && completedStep % energyEvery_ == 0;
    if (sampleElectrodes_) electrodes_.beginStep();
    consoleElectrodes_.beginStep();
    if (sampleEnergy_) energy_.beginStep();
}

void PicDiagnosticsOutput::recordElectrode(int boundaryId, double chargeC, double neutralWeight)
{
    if (sampleElectrodes_) electrodes_.record(boundaryId, chargeC, neutralWeight);
    consoleElectrodes_.record(boundaryId, chargeC, neutralWeight);
}

void PicDiagnosticsOutput::observeEnergy(const std::string& species, double chargeC,
    const std::string& location, double energyEV, double weight)
{
    if (sampleEnergy_) energy_.observe(species, chargeC, location, energyEV, weight);
}

void PicDiagnosticsOutput::writeReactions(int step, double timeS, double dtS, const MCCCollisionOperator::StepStats& stats)
{
    if (sampleReactions_) reactions_.writeStep(step, timeS, dtS, stats);
}

void PicDiagnosticsOutput::writeEnergy(int step, double timeS)
{
    if (sampleEnergy_) energy_.writeStep(step, timeS);
}

void PicDiagnosticsOutput::writeElectrodes(int step, double timeS, double dtS)
{
    if (sampleElectrodes_) electrodes_.writeStep(step, timeS, dtS);
}
} // namespace PhysicsDiagnostics
