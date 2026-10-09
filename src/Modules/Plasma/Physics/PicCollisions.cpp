#include "Modules/Plasma/Physics/PicCollisions.h"
#include "Modules/Plasma/ExplicitPIC/PicConfig.h"
#include "Modules/Plasma/Physics/PhysicalConstants.h"
#include "Modules/Plasma/IO/PicCheckpointReader.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/IO/ConsoleControl.h"
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <utility>

PicCollisions::PicCollisions(PlasmaCollisionInputs inputs, double velocityScale)
    : inputs_(std::move(inputs)), velocityScale_(velocityScale),
      inverseVelocityScale_((velocityScale > 0.0) ? (1.0 / velocityScale) : 1.0)
{}

std::unique_ptr<PicCollisions> PicCollisions::create(PlasmaCollisionInputs inputs,
    const std::vector<SpeciesConfig>& species, double velocityScale, std::string& error,
    const PicCheckpointReader* restart, bool restoreRng)
{
    error.clear();
    if (!inputs.collisions || !inputs.collisions->enabled) return nullptr;
    if (!inputs.gas || !inputs.crossSections || !inputs.reactions)
    {
        error = "collisions enabled but gas/crossSections/reactions are not configured";
        return nullptr;
    }
    auto result = std::unique_ptr<PicCollisions>(new PicCollisions(std::move(inputs), velocityScale));
    for (const auto& sp : species)
    {
        Species s;
        s.name = sp.name;
        s.mass_kg = sp.m;
        s.charge_C = sp.q;
        s.macro_weight = sp.weight;
        if (std::abs(sp.q) < 1e-30)
            s.role = SpeciesRole::kNeutral;
        else if (sp.q < 0.0)
            s.role = SpeciesRole::kElectron;
        else
            s.role = SpeciesRole::kIon;
        s.init.kT_J = std::max(0.0, sp.Te_eV) * PhysicalConstants::eV_J;
        s.init.drift_mps = sp.drift;
        if (!result->species_.Add(s, &error))
        {
            error = "SpeciesDatabase add failed: " + error;
            return nullptr;
        }
    }
    result->operation_ = std::make_unique<MCCCollisionOperator>(result->species_,
        *result->inputs_.gas, *result->inputs_.crossSections, *result->inputs_.reactions);
    if (!result->operation_->Initialize(*result->inputs_.collisions, &error))
    {
        error = "MCCCollisionOperator init failed: " + error;
        return nullptr;
    }
    if (restart) restart->readCollisions(*result, restoreRng);
    return result;
}

bool PicCollisions::apply(ParticleType& projectile, double macroWeight, double dtSI,
    std::mt19937_64& rng, ParticleSpawnBuffer& spawned,
    MCCCollisionOperator::StepStats& stats, std::string* error) const
{
    return operation_->ApplyCodeUnits(projectile, macroWeight, inputs_.gas->number_density_m3,
        dtSI, rng, velocityScale_, inverseVelocityScale_, spawned, stats, error);
}

bool PicCollisions::Stage::applyExplicit(const QuadtreeMesh& mesh, ParticleType& sp,
    double macroWeight, double dtSI, std::mt19937_64& rng, int step)
{
    if (!collisions_ || !collisions_->operation().HasProjectile(sp.name())) return true;
    spawned_.Clear();
    std::string err;

    // MCC uses SI cross-sections but particle velocities are stored in code units.
    // Use ApplyCodeUnits() to avoid O(Np) pre/post scaling loops.
    if (!collisions_->apply(sp, macroWeight, dtSI, rng,
                               spawned_, batchStats_, &err))
    {
        std::cerr << "[ExplicitPIC] ERROR: MCC Apply failed at step " << step
                  << " projectile='" << sp.name() << "': " << err << std::endl;
        return false;
    }

    // Defer merge of spawned particles until *after* all species pushes.
    // Otherwise, a projectile (e.g. electron) can spawn particles into another
    // species (e.g. Ar_plus) during this loop, causing the gather cache for the
    // later species to be built for the old particle count.
    std::size_t kept = 0;
    if (!spawned_.empty())
    {
        std::unordered_map<std::string, std::size_t> counts;
        counts.reserve(8);
        for (const auto& it : spawned_.items())
            counts[it.species_name] += 1;

        for (const auto& it : spawned_.items())
        {
            auto mit = speciesIndex_.find(it.species_name);
            if (mit == speciesIndex_.end())
            {
                std::cerr << "[ExplicitPIC][MCC] WARNING: spawned species '" << it.species_name
                          << "' is not defined in pic.dict; dropping particle" << std::endl;
                continue;
            }

            Particle np = it.particle;
            if (np.w <= 0.0 && mit->second < weights_.size())
                np.w = weights_[mit->second];
            int inCellL = -1;
            if (np.inCellL >= 0)
                inCellL = mesh.findCellLabelFromHint(np.inCellL, np.position);
            if (inCellL < 0)
                inCellL = mesh.findCellLabel(np.position);
            if (inCellL < 0)
                continue;
            np.inCellL = inCellL;
            np.inBackgroundCellL = mesh.findBackgroundCellLabel(np.position);

            pending_.Push(it.species_name, np);
            kept += 1;
        }

        if (ConsoleControl::enabled(ConsoleControl::Level::kTrace))
        {
            std::cout << "[ExplicitPIC][MCC] step " << step
                      << " projectile='" << sp.name() << "' spawned(raw)=" << spawned_.size()
                      << " kept=" << kept;
            if (!counts.empty())
            {
                std::cout << " (";
                bool first = true;
                for (const auto& kv : counts)
                {
                    if (!first)
                        std::cout << ", ";
                    first = false;
                    std::cout << kv.first << "=" << kv.second;
                }
                std::cout << ")";
            }
            std::cout << std::endl;
            std::cout << spawned_.SummaryString(8) << std::endl;
        }
    }
    return true;
}

bool PicCollisions::Stage::applyImplicitHalf(std::vector<ParticleType>& species,
    double dtSI, std::mt19937_64& rng)
{
    if (!collisions_)
        return false;
    if (!(dtSI > 0.0) || !std::isfinite(dtSI))
        return false;

    beginBatch();

    for (std::size_t si = 0; si < species.size(); ++si)
    {
        auto& sp = species[si];
        if (!collisions_->operation().HasProjectile(sp.name()))
            continue;

        spawned_.Clear();
        std::string err;
        if (!collisions_->apply(sp,
                                   (si < weights_.size()) ? weights_[si] : 1.0,
                                   dtSI,
                                   rng,
                                   spawned_,
                                   batchStats_,
                                   &err))
        {
            std::cerr << "[ImplicitPIC] ERROR: MCC Apply failed projectile='" << sp.name() << "': " << err
                      << std::endl;
            throw std::runtime_error("MCC ApplyCodeUnits failed");
        }

        if (!spawned_.empty())
        {
            for (const auto& it : spawned_.items())
            {
                auto mit = speciesIndex_.find(it.species_name);
                if (mit == speciesIndex_.end())
                    continue;

                Particle np = it.particle;
                if (np.w <= 0.0 && mit->second < weights_.size())
                    np.w = weights_[mit->second];
                pending_.Push(it.species_name, np);
            }
        }
    }

    return true;
}

bool PicCollisions::Stage::advanceImplicitHalf(const PicSpawnManagement& management,
    const PicSpawnManagement::ImplicitStage& stage, double dtSI,
    PicSpawnManagement::Stats* statistics)
{
    if (!applyImplicitHalf(stage.species, dtSI, stage.rng)) return false;
    management.acceptImplicit(stage, pending_, PicSpawnManagement::ImplicitKind::Collision, statistics);
    return true;
}

PicSpawnManagement::ExplicitResult PicCollisions::Stage::acceptExplicit(const PicSpawnManagement& management,
    const PicSpawnManagement::ExplicitStage& stage, PicSecondaryEmission& emission,
    std::vector<Particle>& secondary, PicParticleSources& sources, ParticleSpawnBuffer& source)
{
    return management.acceptExplicit(stage, pending_, emission, secondary, sources, source);
}

PicSpawnManagement::Stats PicCollisions::Stage::thinSpawns(const PicSpawnManagement& management,
    bool enabled, double probability, std::mt19937_64& rng,
    std::uniform_real_distribution<double>& uniform)
{
    return management.thin(pending_, enabled, probability, rng, uniform);
}
