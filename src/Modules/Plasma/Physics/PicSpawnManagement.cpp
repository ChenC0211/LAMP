#include "Modules/Plasma/Physics/PicSpawnManagement.h"
#include "Modules/Plasma/Physics/ParticleSpawnBuffer.h"
#include "Modules/Plasma/Physics/PicSecondaryEmission.h"
#include "Modules/Plasma/Physics/PicParticleSources.h"
#include "Modules/Plasma/Particles/ParticleType.h"
#include "Modules/Plasma/Diagnostics/PicStepDiagnostics.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <utility>

namespace {
template<class Items, class ParticleOf, class WeightOf>
PicSpawnManagement::Stats thinItems(Items& items, bool enabled, double probability,
    std::mt19937_64& rng, std::uniform_real_distribution<double>& uniform,
    ParticleOf particleOf, WeightOf weightOf, bool implicit, std::vector<int>* tags = nullptr)
{
    if (tags && tags->size() != items.size())
        throw std::invalid_argument("Spawn thinning requires one tag per particle");
    PicSpawnManagement::Stats stats;
    stats.in = stats.kept = static_cast<long long>(items.size());
    if (stats.in > 0) {
        stats.keepProbWeightedMacro += (enabled ? probability : 1.0) * static_cast<double>(stats.in);
        stats.keepProbMacroCount += stats.in;
    }
    if (items.empty()) return stats;
    for (const auto& item : items) {
        const auto& particle = particleOf(item);
        stats.physIn += particle.w > 0.0 ? particle.w : weightOf(item);
    }
    stats.physKept = stats.physIn;
    if (!enabled || !(probability < 1.0)) return stats;
    std::size_t write = 0;
    double physicalKept = 0.0;
    for (std::size_t j = 0; j < items.size(); ++j) {
        auto& particle = particleOf(items[j]);
        const double weight = particle.w > 0.0 ? particle.w : weightOf(items[j]);
        if (uniform(rng) <= probability) {
            particle.w = weight / probability;
            physicalKept += particle.w;
            if (write != j) {
                items[write] = std::move(items[j]);
                if (tags) (*tags)[write] = (*tags)[j];
            }
            ++write;
        } else ++stats.drop;
    }
    items.resize(write);
    if (tags) tags->resize(write);
    stats.kept = static_cast<long long>(write);
    // Retain each driver's tally arithmetic even for nonfinite legacy weights.
    if (implicit) { stats.physKept -= stats.physIn; stats.physKept += physicalKept; }
    else stats.physKept = physicalKept;
    return stats;
}
}

PicSpawnManagement::PicSpawnManagement(Mode mode, ParticleManagementConfig config,
    const std::unordered_map<std::string,std::size_t>& speciesIndex,
    const std::vector<double>& defaultWeights)
    : mode_(mode), config_(std::move(config)), speciesIndex_(speciesIndex), weights_(defaultWeights) {}

void PicSpawnManagement::Stats::accumulate(const Stats& other)
{
    in += other.in; kept += other.kept; drop += other.drop;
    physIn += other.physIn; physKept += other.physKept;
    keepProbWeightedMacro += other.keepProbWeightedMacro;
    keepProbMacroCount += other.keepProbMacroCount;
}

double PicSpawnManagement::defaultWeight(const std::string& name) const
{
    auto it = speciesIndex_.find(name);
    if (it == speciesIndex_.end()) return 1.0;
    return it->second < weights_.size() ? weights_[it->second] : 1.0;
}

double PicSpawnManagement::keepProbability(const std::vector<ParticleType>& species,
    long long spawnCount, bool enabled) const
{
    if (!enabled)
        return 1.0;
    if (!(config_.softCapTotal > 0 || config_.hardCapTotal > 0))
        return 1.0;

    long long totalMacro = 0;
    for (const auto& sp : species) totalMacro += static_cast<long long>(sp.particle.size());
    const double totalAfter = static_cast<double>(totalMacro + std::max(0ll, spawnCount));
    if (!(totalAfter > 0.0))
        return 1.0;

    double keepProb = 1.0;
    if (config_.softCapTotal > 0 && totalAfter > static_cast<double>(config_.softCapTotal))
        keepProb = std::min(keepProb, static_cast<double>(config_.softCapTotal) / totalAfter);
    if (config_.hardCapTotal > 0 && totalAfter > static_cast<double>(config_.hardCapTotal))
        keepProb = std::min(keepProb, static_cast<double>(config_.hardCapTotal) / totalAfter);
    keepProb = std::min(1.0, keepProb);
    if (keepProb < 1.0)
        keepProb = std::max(keepProb, std::max(0.0, config_.minKeepProb));
    return keepProb;
}

PicSpawnManagement::Stats PicSpawnManagement::thin(ParticleSpawnBuffer& buffer,
    bool enabled, double probability, std::mt19937_64& rng,
    std::uniform_real_distribution<double>& uniform, std::vector<int>* tags) const
{
    return thinItems(buffer.items_mut(), enabled, probability, rng, uniform,
        [](auto& item) -> auto& { return item.particle; },
        [&](const auto& item) { return defaultWeight(item.species_name); }, mode_ == Mode::Implicit, tags);
}

PicSpawnManagement::Stats PicSpawnManagement::thin(std::vector<Particle>& particles,
    std::size_t speciesIndex, bool enabled, double probability, std::mt19937_64& rng,
    std::uniform_real_distribution<double>& uniform, std::vector<int>* tags) const
{
    const double weight = speciesIndex < weights_.size() ? weights_[speciesIndex] : 1.0;
    return thinItems(particles, enabled, probability, rng, uniform,
        [](auto& particle) -> auto& { return particle; },
        [&](const auto&) { return weight; }, mode_ == Mode::Implicit, tags);
}

void PicSpawnManagement::commit(const ParticleSpawnBuffer& buffer, const QuadtreeMesh& mesh,
    std::vector<ParticleType>& species, std::vector<unsigned char>* touched, bool mcc) const
{
    if (buffer.empty()) return;
    if (mode_ == Mode::Explicit)
    {
        const auto& pendingItems = buffer.items();
        std::vector<int> targetIndexByItem(pendingItems.size(), -1);
        std::vector<std::size_t> perTarget(species.size(), 0u);

        for (std::size_t i = 0; i < pendingItems.size(); ++i)
        {
            const auto& it = pendingItems[i];
            auto mit = speciesIndex_.find(it.species_name);
            if (mit == speciesIndex_.end())
                continue;
            const std::size_t dstIdx = mit->second;
            targetIndexByItem[i] = static_cast<int>(dstIdx);
            perTarget[dstIdx] += 1;
        }

        for (std::size_t si = 0; si < species.size(); ++si)
        {
            const std::size_t addCount = perTarget[si];
            if (addCount == 0u)
                continue;
            ParticleType& dst = species[si];
            if (touched) (*touched)[si] = 1u;
            dst.particle.reserve(dst.particle.size() + addCount);
        }

        for (std::size_t i = 0; i < pendingItems.size(); ++i)
        {
            const int dstI = targetIndexByItem[i];
            if (dstI < 0)
                continue;
            species[static_cast<std::size_t>(dstI)].particle.push_back(pendingItems[i].particle);
        }
        for (std::size_t si = 0; si < species.size(); ++si)
        {
            if (perTarget[si] == 0u)
                continue;
            ParticleType& dst = species[si];
            dst.number = static_cast<int>(dst.particle.size());
        }
        return;
    }
    for (const auto& it : buffer.items())
    {
        auto mit = speciesIndex_.find(it.species_name);
        if (mit == speciesIndex_.end())
        {
            if (mcc) std::cerr << "[ImplicitPIC][MCC] WARNING: spawned species '" << it.species_name
                      << "' is not defined in dict; dropping particle" << std::endl;
            continue;
        }

        const std::size_t dstIdx = mit->second;
        Particle np = it.particle;
        if (np.w <= 0.0 && dstIdx < weights_.size())
            np.w = weights_[dstIdx];

        int inCellL = -1;
        if (np.inCellL >= 0)
            inCellL = mesh.findCellLabelFromHint(np.inCellL, np.position);
        if (inCellL < 0)
            inCellL = mesh.findCellLabel(np.position);
        if (inCellL < 0)
            continue;
        np.inCellL = inCellL;
        np.inBackgroundCellL = mesh.findBackgroundCellLabel(np.position);

        species[dstIdx].particle.push_back(np);
    }

    for (auto& sp : species)
        sp.number = static_cast<int>(sp.particle.size());

    // Light sanity checks: keep particle container metadata consistent and
    // ensure cell labels are valid after spawn merge.
    if (mcc) for (const auto& sp : species)
    {
        if (sp.number != static_cast<int>(sp.particle.size()))
            throw std::runtime_error("[ImplicitPIC][MCC] internal error: species.number mismatch after spawn merge");
        for (const auto& p : sp.particle)
        {
            if (p.inCellL < 0 || p.inBackgroundCellL < 0)
                throw std::runtime_error("[ImplicitPIC][MCC] internal error: invalid cell label after spawn merge");
        }
    }
}

PicSpawnManagement::ExplicitResult PicSpawnManagement::acceptExplicit(const ExplicitStage& stage,
    ParticleSpawnBuffer& collision, PicSecondaryEmission& emission, std::vector<Particle>& secondary,
    PicParticleSources& sources, ParticleSpawnBuffer& source) const
{
    if (mode_ != Mode::Explicit)
        throw std::logic_error("Explicit spawn stages require explicit management");
    ExplicitResult result;
    const long long spawnCount = static_cast<long long>(collision.size()) +
        static_cast<long long>(secondary.size()) + static_cast<long long>(source.size());
    const bool thinning = stage.managementActive && config_.thinningEnabled;
    result.keepProbability = keepProbability(stage.species, spawnCount, thinning);
    result.collision = thin(collision, thinning && config_.thinningMccSpawn, result.keepProbability, stage.rng, stage.uniform);
    result.secondary = emission.thinExplicit(*this, secondary, thinning && config_.thinningSeeSpawn, result.keepProbability, stage.rng, stage.uniform);
    result.source = sources.thinExplicit(*this, source, thinning && config_.thinningMccSpawn, result.keepProbability, stage.rng, stage.uniform);
    emission.accountRetained(secondary, stage.species, weights_);
    const auto retained = sources.accountRetained(source, stage.species, weights_, speciesIndex_, stage.velocityScale, stage.domainBoundaryIds);
    result.injectionEnergyJ = emission.stats().emittedEnergyJ + retained.energyJ;
    result.sourceChargeC = retained.domainChargeC;
    result.touchedSpecies.assign(stage.species.size(), 0u);
    commit(collision, stage.mesh, stage.species, &result.touchedSpecies);
    commit(source, stage.mesh, stage.species, &result.touchedSpecies);
    if (emission.enabled() && emission.emittedIndex() != static_cast<std::size_t>(-1) && !secondary.empty())
        commitExplicitSee(secondary, emission.emittedIndex(), stage.species, result.touchedSpecies);
    return result;
}

void PicSpawnManagement::acceptImplicit(const ImplicitStage& stage, ParticleSpawnBuffer& buffer,
    ImplicitKind kind, Stats* statistics) const
{
    if (mode_ != Mode::Implicit)
        throw std::logic_error("Implicit spawn stages require implicit management");
    if (kind == ImplicitKind::Secondary && buffer.empty()) return;
    const bool thinning = stage.managementActive && config_.thinningEnabled &&
        (kind == ImplicitKind::Secondary ? config_.thinningSeeSpawn : config_.thinningMccSpawn);
    const double probability = keepProbability(stage.species, static_cast<long long>(buffer.size()), thinning);
    const auto local = thin(buffer, thinning, probability, stage.rng, stage.uniform);
    if (statistics)
    {
        if (kind == ImplicitKind::Collision) statistics->accumulate(local);
        else *statistics = local;
    }
    using Spawn = PicStepParticleBalance::Spawn;
    const Spawn spawn = kind == ImplicitKind::Collision ? Spawn::Collision :
        kind == ImplicitKind::Secondary ? Spawn::Secondary : Spawn::Source;
    if (kind != ImplicitKind::Source)
        stage.balance.addSpawns(spawn, buffer, stage.species, weights_, speciesIndex_, stage.velocityScale);
    commit(buffer, stage.mesh, stage.species, nullptr, kind == ImplicitKind::Collision);
    if (kind == ImplicitKind::Source)
        stage.balance.addSpawns(spawn, buffer, stage.species, weights_, speciesIndex_, stage.velocityScale);
    if (kind == ImplicitKind::Secondary) stage.balance.includeSecondarySource();
}

void PicSpawnManagement::commitExplicitSee(const std::vector<Particle>& particles,
    std::size_t speciesIndex, std::vector<ParticleType>& species, std::vector<unsigned char>& touched) const
{
    if (particles.empty()) return;
    ParticleType& dst = species[speciesIndex];
    touched[speciesIndex] = 1u;
    dst.particle.reserve(dst.particle.size() + particles.size());
    for (const auto& p : particles) dst.particle.push_back(p);
    dst.number = static_cast<int>(dst.particle.size());
}
