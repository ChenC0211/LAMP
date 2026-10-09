#include "Modules/Plasma/Physics/MCCCollisionOperator.h"

#include <algorithm>
#include <atomic>
#include <cctype>
#include <chrono>
#include <cmath>
#include <memory>
#include <ostream>
#include <sstream>
#include <limits>
#include <iterator>
#include <set>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "Modules/Plasma/Physics/PhysicalConstants.h"
#include "Core/IO/OMPPlacementDiag.h"
#include "Core/Random/FastRng.h"

using namespace std;

void MCCCollisionOperator::restoreRngState(const RngState& state)
{
    if (state.mode != cfg_.rng_mode || state.seed != cfg_.rng_seed)
        throw std::runtime_error("MCC RNG checkpoint mode/seed does not match collision configuration");
    apply_counter_ = state.applyCount;
}

namespace
{
std::string normalizeEnergyFrame(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

struct RelativeState
{
    PostCollisionModel::TargetState target;
    double v_rel_mps = 0.0;
    double lookup_E_eV = 0.0;
};

bool resolveLookupMassKg(const GasState& gas,
                         const ReactionSet& reactionSet,
                         double projectile_mass_kg,
                         double* lookup_mass_kg,
                         std::string* err)
{
    if (!lookup_mass_kg)
        return false;

    *lookup_mass_kg = projectile_mass_kg;
    const std::string energyFrame = normalizeEnergyFrame(reactionSet.energy_frame);
    if (energyFrame == "centerofmass")
    {
        const double m_tgt = gas.particle_mass_kg;
        if (!(m_tgt > 0.0) || !std::isfinite(m_tgt))
        {
            if (err)
                *err = "ReactionSet '" + reactionSet.name + "' uses energyFrame=centerOfMass but gas particle mass is invalid";
            return false;
        }
        const double mu = (projectile_mass_kg * m_tgt) / (projectile_mass_kg + m_tgt);
        if (!(mu > 0.0) || !std::isfinite(mu))
        {
            if (err)
                *err = "ReactionSet '" + reactionSet.name + "' has invalid reduced mass for energyFrame=centerOfMass";
            return false;
        }
        *lookup_mass_kg = mu;
        return true;
    }

    if (energyFrame != "lab")
    {
        if (err)
            *err = "ReactionSet '" + reactionSet.name + "' has unsupported energyFrame='" + reactionSet.energy_frame + "'";
        return false;
    }
    return true;
}

template <class Rng>
bool sampleGasTargetState(const GasState& gas,
                          bool use_gas_thermal_motion,
                          Rng& rng,
                          PostCollisionModel::TargetState* target,
                          std::string* err)
{
    if (!target)
        return false;

    *target = PostCollisionModel::TargetState{};
    if (!use_gas_thermal_motion)
        return true;

    if (!(gas.particle_mass_kg > 0.0) || !std::isfinite(gas.particle_mass_kg))
    {
        if (err)
            *err = "useGasThermalMotion=true requires a valid background gas particle mass";
        return false;
    }
    if (!(gas.temperature_K > 0.0) || !std::isfinite(gas.temperature_K))
    {
        if (err)
            *err = "useGasThermalMotion=true requires a valid positive background gas temperature";
        return false;
    }

    const double sigma = std::sqrt(PhysicalConstants::kB_JpK * gas.temperature_K / gas.particle_mass_kg);
    std::normal_distribution<double> nd(0.0, sigma);
    target->velocity = Vector2D(nd(rng), nd(rng));
    target->vz = nd(rng);
    return true;
}

template <class Rng>
bool buildRelativeState(const Particle& projectile_si,
                        const GasState& gas,
                        bool use_gas_thermal_motion,
                        double lookup_mass_kg,
                        Rng& rng,
                        RelativeState* rel,
                        std::string* err)
{
    if (!rel)
        return false;
    if (!sampleGasTargetState(gas, use_gas_thermal_motion, rng, &rel->target, err))
        return false;

    const double dvx = projectile_si.velocity.x - rel->target.velocity.x;
    const double dvy = projectile_si.velocity.y - rel->target.velocity.y;
    const double dvz = projectile_si.vz - rel->target.vz;
    const double v2_rel = dvx * dvx + dvy * dvy + dvz * dvz;
    rel->v_rel_mps = std::sqrt(std::max(0.0, v2_rel));

    // The lookup energy is evaluated in the configured frame using the current
    // relative speed. The null-collision rate still uses n_gas * sigma(E) * v_rel.
    rel->lookup_E_eV = 0.5 * lookup_mass_kg * v2_rel / PhysicalConstants::eV_J;
    return true;
}

template <class Rng>
bool applyParticleEventsSI(Particle& particle_si,
                           double particle_weight,
                           double n_gas_m3,
                           double dt_s,
                           const GasState& gas,
                           bool use_gas_thermal_motion,
                           double lookup_mass_kg,
                           const NullCollisionOperator& null_op,
                           const ReactionSet& reactionSet,
                           const std::vector<const ReactionChannel*>& enabled_channels,
                           const PostCollisionModel& post,
                           Rng& rng,
                           ParticleSpawnBuffer& spawnBuffer,
                           std::vector<std::uint64_t>* channel_counts,
                           std::vector<double>* channel_counts_weighted,
                           std::vector<double>* channel_energy_loss_J,
                           std::uint64_t* n_attempts,
                           double* n_attempts_weighted,
                           std::uint64_t* n_null,
                           double* n_null_weighted,
                           std::uint64_t* n_collisions,
                           double* n_collisions_weighted,
                           double* energy_loss_eV_macro,
                           double* energy_loss_J_total,
                           double* max_nu_dt_majorant,
                           double* t_sample_s,
                           double* t_post_s,
                           const std::string& projectile_name,
                           std::string* err)
{
    using Clock = std::chrono::steady_clock;
    auto secondsBetween = [](const Clock::time_point& a, const Clock::time_point& b) -> double {
        return std::chrono::duration_cast<std::chrono::duration<double>>(b - a).count();
    };

    if (!(dt_s > 0.0) || !(n_gas_m3 > 0.0))
        return true;

    double remaining_dt = dt_s;
    while (remaining_dt > 0.0)
    {
        RelativeState rel;
        if (!buildRelativeState(particle_si, gas, use_gas_thermal_motion, lookup_mass_kg, rng, &rel, err))
            return false;
        if (!(rel.v_rel_mps > 0.0))
            break;

        const double sigmaM = null_op.sigma_majorant_m2();
        const double nu_dt = n_gas_m3 * rel.v_rel_mps * sigmaM * remaining_dt;
        if (max_nu_dt_majorant && std::isfinite(nu_dt))
            *max_nu_dt_majorant = std::max(*max_nu_dt_majorant, nu_dt);

        NullCollisionOperator::NextEvent event;
        NullCollisionCounts counts;
        counts.Clear();
        {
            const auto t0 = Clock::now();
            null_op.SampleNextEvent(n_gas_m3, rel.v_rel_mps, rel.lookup_E_eV, remaining_dt, rng, &event, &counts);
            if (t_sample_s)
                *t_sample_s += secondsBetween(t0, Clock::now());
        }

        if (n_attempts)
            *n_attempts += counts.n_attempts;
        if (n_attempts_weighted)
            *n_attempts_weighted += static_cast<double>(counts.n_attempts) * std::max(0.0, particle_weight);
        if (n_null)
            *n_null += counts.n_null;
        if (n_null_weighted)
            *n_null_weighted += static_cast<double>(counts.n_null) * std::max(0.0, particle_weight);
        if (n_collisions)
            *n_collisions += counts.n_collisions;
        if (n_collisions_weighted)
            *n_collisions_weighted += static_cast<double>(counts.n_collisions) * std::max(0.0, particle_weight);

        const double consumed_dt = std::min(remaining_dt, std::max(0.0, event.elapsed_s));
        remaining_dt -= consumed_dt;

        if (event.kind == NullCollisionOperator::NextEvent::Kind::kNone)
            break;
        if (event.kind == NullCollisionOperator::NextEvent::Kind::kNull)
            continue;
        if (event.channel_index < 0 || static_cast<std::size_t>(event.channel_index) >= enabled_channels.size())
            continue;

        const ReactionChannel* ch = enabled_channels[static_cast<std::size_t>(event.channel_index)];
        if (!ch)
            continue;

        PostCollisionModel::Outcome out;
        std::string localErr;
        {
            const auto t0 = Clock::now();
            const bool ok = post.Apply(reactionSet, *ch, particle_si, rng, &rel.target, &spawnBuffer, &out, &localErr);
            if (t_post_s)
                *t_post_s += secondsBetween(t0, Clock::now());
            if (!ok)
            {
                if (err)
                    *err = "PostCollisionModel failed for projectile '" + projectile_name + "' channel '" + ch->name + "': " + localErr;
                return false;
            }
        }

        if (channel_counts && static_cast<std::size_t>(event.channel_index) < channel_counts->size())
            (*channel_counts)[static_cast<std::size_t>(event.channel_index)] += 1;
        if (channel_counts_weighted && static_cast<std::size_t>(event.channel_index) < channel_counts_weighted->size())
            (*channel_counts_weighted)[static_cast<std::size_t>(event.channel_index)] += std::max(0.0, particle_weight);

        const double dE = std::max(0.0, out.E_before_eV - out.E_after_eV);
        const double dE_J = dE * PhysicalConstants::eV_J * std::max(0.0, particle_weight);
        if (channel_energy_loss_J && static_cast<std::size_t>(event.channel_index) < channel_energy_loss_J->size())
            (*channel_energy_loss_J)[static_cast<std::size_t>(event.channel_index)] += dE_J;
        if (energy_loss_eV_macro)
            *energy_loss_eV_macro += dE;
        if (energy_loss_J_total)
            *energy_loss_J_total += dE_J;
    }

    return true;
}
}

MCCCollisionOperator::MCCCollisionOperator(const SpeciesDatabase& species_db,
                                         const GasState& gas,
                                         const CrossSectionDB& xs,
                                         const ReactionDatabase& reactions)
    : species_db_(species_db)
    , gas_(gas)
    , xs_(xs)
    , reactions_(reactions)
{
}

bool MCCCollisionOperator::Initialize(const CollisionsConfig& cfg, std::string* err)
{
    cfg_ = cfg;
    ops_.clear();
    summary_ = Summary{};

    std::string verr;
    if (!cfg_.Validate(&verr))
    {
        if (err) *err = verr;
        return false;
    }
    if (!cfg_.enabled)
        return true;

    for (const auto& m : cfg_.mappings)
    {
        const auto* set = reactions_.Get(m.reaction_set);
        if (!set)
        {
            if (err) *err = "ReactionSet not found: '" + m.reaction_set + "'";
            return false;
        }

        ProjectileOps po;
        po.projectile = m.projectile;
        po.reaction_set_key = m.reaction_set;
        PostCollisionModel::Options popt;
        popt.gas_species_name = gas_.name;
        popt.gas_temperature_K = gas_.temperature_K;
        popt.gas_particle_mass_kg = gas_.particle_mass_kg;
        // Ion-neutral ELASTIC+CX slow-down without fast-neutral tracking.
        popt.ion_neutral_elastic_thermalize = true;
        popt.ion_neutral_cx_thermalize = true;
        popt.track_fast_neutrals = cfg_.track_fast_neutrals;
        po.post = std::make_unique<PostCollisionModel>(species_db_, popt);

        std::string initErr;
        if (!po.null_op.Initialize(*set, xs_, cfg_.null_collision, &initErr))
        {
            if (err) *err = "NullCollisionOperator init failed for projectile '" + m.projectile + "': " + initErr;
            return false;
        }

        // Build enabled channel pointer list matching sampling indices.
        for (const auto& ch : po.null_op.reaction_set().channels)
        {
            if (!ch.enabled)
                continue;
            po.enabled_channels.push_back(&ch);
        }

        summary_.n_projectiles += 1;
        summary_.n_channels_total += static_cast<uint32_t>(po.enabled_channels.size());

        ops_.emplace(po.projectile, std::move(po));
    }

    return true;
}

bool MCCCollisionOperator::HasProjectile(const std::string& projectile) const
{
    return ops_.find(projectile) != ops_.end();
}

void MCCCollisionOperator::PrintSummary(std::ostream& os) const
{
    os << "MCCCollisionOperator\n";
    os << "  enabled         : " << (cfg_.enabled ? "true" : "false") << "\n";
    os << "  gas             : " << gas_.name << "\n";
    os << "  n_projectiles   : " << summary_.n_projectiles << "\n";
    os << "  n_channels_total: " << summary_.n_channels_total << "\n";
    for (const auto& kv : ops_)
    {
        const auto& po = kv.second;
        os << "  - projectile=" << po.projectile
           << " reactionSet=" << po.reaction_set_key
           << " sigmaM=" << po.null_op.sigma_majorant_m2() << " m^2"
           << " channels=" << po.enabled_channels.size() << "\n";
    }
}

std::string MCCCollisionOperator::SummaryString() const
{
    std::ostringstream ss;
    PrintSummary(ss);
    return ss.str();
}

std::vector<std::string> MCCCollisionOperator::ProjectileNames() const
{
    std::vector<std::string> out;
    out.reserve(ops_.size());
    for (const auto& kv : ops_)
        out.push_back(kv.first);
    std::sort(out.begin(), out.end());
    return out;
}

std::vector<std::string> MCCCollisionOperator::ChannelKeys() const
{
    std::vector<std::string> out;
    for (const auto& kv : ops_)
    {
        const auto& projectile = kv.first;
        const auto& po = kv.second;
        for (const ReactionChannel* ch : po.enabled_channels)
        {
            if (!ch)
                continue;
            out.push_back(projectile + "::" + ch->name);
        }
    }
    std::sort(out.begin(), out.end());
    return out;
}

std::string MCCCollisionOperator::ChannelKind(const std::string& projectile, const std::string& channel) const
{
    const auto it = ops_.find(projectile);
    if (it == ops_.end()) return {};
    const auto& po = it->second;
    for (const auto* ch : po.enabled_channels)
        if (ch && ch->name == channel)
            return PostCollisionModel::InferKind(po.null_op.reaction_set(), *ch);
    return {};
}

std::vector<std::string> MCCCollisionOperator::SpawnSpeciesNames() const
{
    std::set<std::string> names;
    for (const auto& kv : ops_)
    {
        const auto& po = kv.second;
        for (const ReactionChannel* ch : po.enabled_channels)
        {
            if (!ch)
                continue;
            for (const auto& speciesName : ch->products)
            {
                if (!speciesName.empty())
                    names.insert(speciesName);
            }
        }
    }
    return std::vector<std::string>(names.begin(), names.end());
}

std::vector<MCCCollisionOperator::MajorantInfo> MCCCollisionOperator::MajorantTable() const
{
    std::vector<MajorantInfo> out;
    out.reserve(ops_.size());
    for (const auto& kv : ops_)
    {
        MajorantInfo info;
        info.projectile = kv.first;
        info.sigma_majorant_m2 = kv.second.null_op.sigma_majorant_m2();
        out.push_back(std::move(info));
    }
    return out;
}

void MCCCollisionOperator::StepStats::Clear()
{
    // totals
    n_particles = 0;
    n_attempts = 0;
    n_null = 0;
    n_collisions = 0;
    n_spawned = 0;
    n_particles_weighted = 0.0;
    n_attempts_weighted = 0.0;
    n_null_weighted = 0.0;
    n_collisions_weighted = 0.0;
    n_spawned_weighted = 0.0;
    energy_loss_eV_macro = 0.0;
    energy_loss_J_total = 0.0;
    max_nu_dt_majorant = 0.0;
    per_channel.clear();
    per_channel_weighted.clear();
    per_channel_energy_loss_J.clear();
    spawn_per_species.clear();
    spawn_per_species_weighted.clear();
    per_projectile.clear();
    timing.Clear();
}

void MCCCollisionOperator::StepStats::ProjectileStats::Clear()
{
    n_particles = 0;
    n_attempts = 0;
    n_null = 0;
    n_collisions = 0;
    n_spawned = 0;
    n_particles_weighted = 0.0;
    n_attempts_weighted = 0.0;
    n_null_weighted = 0.0;
    n_collisions_weighted = 0.0;
    n_spawned_weighted = 0.0;
    energy_loss_eV_macro = 0.0;
    energy_loss_J_total = 0.0;
    max_nu_dt_majorant = 0.0;
    per_channel.clear();
    per_channel_weighted.clear();
    per_channel_energy_loss_J.clear();
    spawn_per_species.clear();
    spawn_per_species_weighted.clear();
    timing.Clear();
}

void MCCCollisionOperator::StepStats::ProjectileStats::Accumulate(const ProjectileStats& other)
{
    n_particles += other.n_particles;
    n_attempts += other.n_attempts;
    n_null += other.n_null;
    n_collisions += other.n_collisions;
    n_spawned += other.n_spawned;
    n_particles_weighted += other.n_particles_weighted;
    n_attempts_weighted += other.n_attempts_weighted;
    n_null_weighted += other.n_null_weighted;
    n_collisions_weighted += other.n_collisions_weighted;
    n_spawned_weighted += other.n_spawned_weighted;
    energy_loss_eV_macro += other.energy_loss_eV_macro;
    energy_loss_J_total += other.energy_loss_J_total;
    max_nu_dt_majorant = std::max(max_nu_dt_majorant, other.max_nu_dt_majorant);
    for (const auto& kv : other.per_channel)
        per_channel[kv.first] += kv.second;
    for (const auto& kv : other.per_channel_weighted)
        per_channel_weighted[kv.first] += kv.second;
    for (const auto& kv : other.per_channel_energy_loss_J)
        per_channel_energy_loss_J[kv.first] += kv.second;
    for (const auto& kv : other.spawn_per_species)
        spawn_per_species[kv.first] += kv.second;
    for (const auto& kv : other.spawn_per_species_weighted)
        spawn_per_species_weighted[kv.first] += kv.second;
    timing.Accumulate(other.timing);
}

void MCCCollisionOperator::StepStats::Accumulate(const StepStats& other)
{
    n_particles += other.n_particles;
    n_attempts += other.n_attempts;
    n_null += other.n_null;
    n_collisions += other.n_collisions;
    n_spawned += other.n_spawned;
    n_particles_weighted += other.n_particles_weighted;
    n_attempts_weighted += other.n_attempts_weighted;
    n_null_weighted += other.n_null_weighted;
    n_collisions_weighted += other.n_collisions_weighted;
    n_spawned_weighted += other.n_spawned_weighted;
    energy_loss_eV_macro += other.energy_loss_eV_macro;
    energy_loss_J_total += other.energy_loss_J_total;
    max_nu_dt_majorant = std::max(max_nu_dt_majorant, other.max_nu_dt_majorant);
    for (const auto& kv : other.per_channel)
        per_channel[kv.first] += kv.second;
    for (const auto& kv : other.per_channel_weighted)
        per_channel_weighted[kv.first] += kv.second;
    for (const auto& kv : other.per_channel_energy_loss_J)
        per_channel_energy_loss_J[kv.first] += kv.second;
    for (const auto& kv : other.spawn_per_species)
        spawn_per_species[kv.first] += kv.second;
    for (const auto& kv : other.spawn_per_species_weighted)
        spawn_per_species_weighted[kv.first] += kv.second;
    for (const auto& kv : other.per_projectile)
        per_projectile[kv.first].Accumulate(kv.second);
    timing.Accumulate(other.timing);
}

void MCCCollisionOperator::StepStats::PrintSummary(std::ostream& os) const
{
    os << "MCC step stats: particles=" << n_particles
         << " particles_phys=" << n_particles_weighted
       << " attempts=" << n_attempts
         << " attempts_phys=" << n_attempts_weighted
       << " null=" << n_null
         << " null_phys=" << n_null_weighted
       << " collisions=" << n_collisions
         << " collisions_phys=" << n_collisions_weighted
       << " spawned=" << n_spawned
         << " spawned_phys=" << n_spawned_weighted
       << " dE_macro_eV=" << energy_loss_eV_macro
       << " dE_J=" << energy_loss_J_total
       << " max_nu_dtM=" << max_nu_dt_majorant;
    if (!per_channel.empty())
    {
        os << "\n  per_channel:";
        for (const auto& kv : per_channel)
            os << " " << kv.first << "=" << kv.second;
    }
    if (!spawn_per_species.empty())
    {
        os << "\n  spawn_per_species:";
        for (const auto& kv : spawn_per_species)
            os << " " << kv.first << "=" << kv.second;
    }
    os << "\n";
}

std::string MCCCollisionOperator::StepStats::SummaryString() const
{
    std::ostringstream ss;
    PrintSummary(ss);
    return ss.str();
}

bool MCCCollisionOperator::Apply(ParticleType& projectile,
                                 double macro_weight,
                                 double n_gas_m3,
                                 double dt_s,
                                 std::mt19937_64& rng,
                                 ParticleSpawnBuffer& spawnBuffer,
                                 StepStats& stepStats,
                                 std::string* err) const
{
    if (!cfg_.enabled)
        return true;

    using Clock = std::chrono::steady_clock;
    auto secondsBetween = [](const Clock::time_point& a, const Clock::time_point& b) -> double {
        return std::chrono::duration_cast<std::chrono::duration<double>>(b - a).count();
    };
    const auto tApply0 = Clock::now();
    double t_sample_s = 0.0;
    double t_post_s = 0.0;
    double t_spawn_s = 0.0;

    const auto it = ops_.find(projectile.name());
    if (it == ops_.end())
        return true; // not configured

    if (!(dt_s >= 0.0))
    {
        if (err)
            *err = "Invalid dt_s";
        return false;
    }

    const Species* sp = species_db_.Find(projectile.name());
    if (!sp)
    {
        if (err)
            *err = "Unknown projectile species: " + projectile.name();
        return false;
    }
    const double m = sp->mass_kg;
    if (!(m > 0.0) || !std::isfinite(m))
    {
        if (err)
            *err = "Invalid projectile mass for species: " + projectile.name();
        return false;
    }

    const ProjectileOps& po = it->second;
    const ReactionSet& reactionSet = po.null_op.reaction_set();
    double lookup_mass_kg = m;
    if (!resolveLookupMassKg(gas_, reactionSet, m, &lookup_mass_kg, err))
        return false;
    const std::string& pname = projectile.name();
    StepStats::ProjectileStats& pstats = stepStats.per_projectile[pname];

    const CollisionsRngMode rng_mode = cfg_.rng_mode;
    const std::uint64_t call_counter = ++apply_counter_;
    auto ensureRngState = [&](Particle& p, std::uint64_t key) {
        if (p.rngState == 0)
        {
            // Seed is deterministic given cfg_.rng_seed and a stable key.
            p.rngState = FastRng::SeedFrom(cfg_.rng_seed, key);
        }
    };

    const std::size_t n0 = spawnBuffer.size();

#ifdef _OPENMP
    const bool can_parallel = cfg_.parallel &&
                              ((rng_mode == CollisionsRngMode::kPerParticle) || (rng_mode == CollisionsRngMode::kCounterBased));
#else
    const bool can_parallel = false;
#endif

#ifdef _OPENMP
    // Parallel particle loop (thread-safe RNG modes only). Stats/diagnostics are kept minimal
    // here; detailed per-channel and per-species merges are handled later.
    if (can_parallel)
    {
        struct LocalTotals
        {
            std::uint64_t n_particles = 0;
            std::uint64_t n_attempts = 0;
            std::uint64_t n_null = 0;
            std::uint64_t n_collisions = 0;
            double n_particles_weighted = 0.0;
            double n_attempts_weighted = 0.0;
            double n_null_weighted = 0.0;
            double n_collisions_weighted = 0.0;
            double energy_loss_eV_macro = 0.0;
            double energy_loss_J_total = 0.0;
            double max_nu_dt_majorant = 0.0;
            double t_sample_s = 0.0;
            double t_post_s = 0.0;
        };

        const int nthreads = std::max(1, omp_get_max_threads());
        std::vector<LocalTotals> tls_totals(static_cast<std::size_t>(nthreads));
        std::vector<ParticleSpawnBuffer> tls_spawn(static_cast<std::size_t>(nthreads));

        const std::size_t nch = po.enabled_channels.size();
        std::vector<std::vector<std::uint64_t>> tls_channel_counts(static_cast<std::size_t>(nthreads), std::vector<std::uint64_t>(nch, 0));
        std::vector<std::vector<double>> tls_channel_counts_weighted(static_cast<std::size_t>(nthreads), std::vector<double>(nch, 0.0));
        std::vector<std::vector<double>> tls_channel_energy_loss_J(static_cast<std::size_t>(nthreads), std::vector<double>(nch, 0.0));

        std::atomic<bool> had_error(false);
        std::string first_error;

#pragma omp parallel
        {
            const int tid = omp_get_thread_num();
            LocalTotals& lt = tls_totals[static_cast<std::size_t>(tid)];
            ParticleSpawnBuffer& localSpawn = tls_spawn[static_cast<std::size_t>(tid)];
            std::vector<int> accepted;
            accepted.reserve(4);

            auto& chCounts = tls_channel_counts[static_cast<std::size_t>(tid)];
            auto& chCountsWeighted = tls_channel_counts_weighted[static_cast<std::size_t>(tid)];
            auto& chEnergyLossJ = tls_channel_energy_loss_J[static_cast<std::size_t>(tid)];

#pragma omp for schedule(static)
            for (std::size_t pi = 0; pi < projectile.particle.size(); ++pi)
            {
                auto& p = projectile.particle[pi];

                if (had_error.load())
                    continue;
                const double w = (p.w > 0.0) ? p.w : macro_weight;
                lt.n_particles += 1;
                lt.n_particles_weighted += std::max(0.0, w);

                NullCollisionCounts nc;
                nc.Clear();
                FastRng::SplitMix64Rng prng(0);
                FastRng::SplitMix64Rng* frng = nullptr;
                if (rng_mode == CollisionsRngMode::kPerParticle)
                {
                    ensureRngState(p, static_cast<std::uint64_t>(pi + 1));
                    prng = FastRng::SplitMix64Rng(p.rngState);
                    frng = &prng;
                }
                else if (rng_mode == CollisionsRngMode::kCounterBased)
                {
                    ensureRngState(p, static_cast<std::uint64_t>(pi + 1));
                    const std::uint64_t seed = FastRng::Mix64(p.rngState ^ FastRng::Mix64(call_counter));
                    prng = FastRng::SplitMix64Rng(seed);
                    frng = &prng;
                }
                std::string localErr;
                const bool ok = applyParticleEventsSI(p,
                                                      w,
                                                      n_gas_m3,
                                                      dt_s,
                                                      gas_,
                                                      cfg_.use_gas_thermal_motion,
                                                      lookup_mass_kg,
                                                      po.null_op,
                                                      reactionSet,
                                                      po.enabled_channels,
                                                      *po.post,
                                                      *frng,
                                                      localSpawn,
                                                      &chCounts,
                                                      &chCountsWeighted,
                                                      &chEnergyLossJ,
                                                      &nc.n_attempts,
                                                      &lt.n_attempts_weighted,
                                                      &nc.n_null,
                                                      &lt.n_null_weighted,
                                                      &nc.n_collisions,
                                                      &lt.n_collisions_weighted,
                                                      &lt.energy_loss_eV_macro,
                                                      &lt.energy_loss_J_total,
                                                      &lt.max_nu_dt_majorant,
                                                      &lt.t_sample_s,
                                                      &lt.t_post_s,
                                                      pname,
                                                      &localErr);
                lt.n_attempts += nc.n_attempts;
                lt.n_null += nc.n_null;
                lt.n_collisions += nc.n_collisions;
                if (!ok)
                {
#pragma omp critical
                    {
                        if (!had_error.load())
                        {
                            had_error.store(true);
                            first_error = localErr;
                        }
                    }
                    continue;
                }

                // Persist per-particle RNG advancement.
                if (rng_mode == CollisionsRngMode::kPerParticle && frng)
                {
                    p.rngState = frng->state();
                }
            }
        } // omp parallel

        if (had_error.load())
        {
            if (err)
                *err = first_error.empty() ? std::string("MCC collision failed in parallel region") : first_error;
            return false;
        }

        // Merge basic totals (per-channel merges handled below; spawn-per-species is counted after merge).
        for (const auto& lt : tls_totals)
        {
            stepStats.n_particles += lt.n_particles;
            pstats.n_particles += lt.n_particles;

            stepStats.n_attempts += lt.n_attempts;
            stepStats.n_null += lt.n_null;
            stepStats.n_collisions += lt.n_collisions;
            stepStats.n_particles_weighted += lt.n_particles_weighted;
            stepStats.n_attempts_weighted += lt.n_attempts_weighted;
            stepStats.n_null_weighted += lt.n_null_weighted;
            stepStats.n_collisions_weighted += lt.n_collisions_weighted;

            pstats.n_attempts += lt.n_attempts;
            pstats.n_null += lt.n_null;
            pstats.n_collisions += lt.n_collisions;
            pstats.n_particles_weighted += lt.n_particles_weighted;
            pstats.n_attempts_weighted += lt.n_attempts_weighted;
            pstats.n_null_weighted += lt.n_null_weighted;
            pstats.n_collisions_weighted += lt.n_collisions_weighted;

            stepStats.energy_loss_eV_macro += lt.energy_loss_eV_macro;
            stepStats.energy_loss_J_total += lt.energy_loss_J_total;
            pstats.energy_loss_eV_macro += lt.energy_loss_eV_macro;
            pstats.energy_loss_J_total += lt.energy_loss_J_total;

            stepStats.max_nu_dt_majorant = std::max(stepStats.max_nu_dt_majorant, lt.max_nu_dt_majorant);
            pstats.max_nu_dt_majorant = std::max(pstats.max_nu_dt_majorant, lt.max_nu_dt_majorant);

            t_sample_s += lt.t_sample_s;
            t_post_s += lt.t_post_s;
        }


        // Merge per-channel hit counts into step/projectile maps.
        if (nch > 0)
        {
            for (std::size_t ci = 0; ci < nch; ++ci)
            {
                std::uint64_t csum = 0;
                double csumWeighted = 0.0;
                double csumEnergyLossJ = 0.0;
                for (const auto& v : tls_channel_counts)
                    csum += v[ci];
                for (const auto& v : tls_channel_counts_weighted)
                    csumWeighted += v[ci];
                for (const auto& v : tls_channel_energy_loss_J)
                    csumEnergyLossJ += v[ci];
                if (csum == 0)
                    continue;
                const ReactionChannel* ch = po.enabled_channels[ci];
                if (!ch)
                    continue;
                stepStats.per_channel[pname + "::" + ch->name] += csum;
                stepStats.per_channel_weighted[pname + "::" + ch->name] += csumWeighted;
                stepStats.per_channel_energy_loss_J[pname + "::" + ch->name] += csumEnergyLossJ;
                pstats.per_channel[ch->name] += csum;
                pstats.per_channel_weighted[ch->name] += csumWeighted;
                pstats.per_channel_energy_loss_J[ch->name] += csumEnergyLossJ;
            }
        }

        // Merge spawned particles deterministically by thread id.
        if (!tls_spawn.empty())
        {
            auto& dst = spawnBuffer.items_mut();
            for (auto& sb : tls_spawn)
            {
                auto& src = sb.items_mut();
                dst.insert(dst.end(), std::make_move_iterator(src.begin()), std::make_move_iterator(src.end()));
                src.clear();
            }
        }
    }
    else
#endif
    {
        std::vector<int> accepted;
        accepted.reserve(4);

        const std::size_t nch = po.enabled_channels.size();
        // Accumulate per-channel counts in a dense vector and merge once per Apply()
        // to avoid per-collision unordered_map/string overhead.
        std::vector<std::uint64_t> channel_counts(nch, 0);
        std::vector<double> channel_counts_weighted(nch, 0.0);
        std::vector<double> channel_energy_loss_J(nch, 0.0);

        for (std::size_t pi = 0; pi < projectile.particle.size(); ++pi)
        {
            auto& p = projectile.particle[pi];
            const double w = (p.w > 0.0) ? p.w : macro_weight;
            stepStats.n_particles += 1;
            pstats.n_particles += 1;
            stepStats.n_particles_weighted += std::max(0.0, w);
            pstats.n_particles_weighted += std::max(0.0, w);

            NullCollisionCounts nc;
            nc.Clear();
            FastRng::SplitMix64Rng prng(0);
            FastRng::SplitMix64Rng* frng = nullptr;
            if (rng_mode == CollisionsRngMode::kPerParticle)
            {
                // Stable key within this projectile container.
                ensureRngState(p, static_cast<std::uint64_t>(pi + 1));
                prng = FastRng::SplitMix64Rng(p.rngState);
                frng = &prng;
            }
            else if (rng_mode == CollisionsRngMode::kCounterBased)
            {
                // Use rngState as a persistent key, but do not advance it.
                ensureRngState(p, static_cast<std::uint64_t>(pi + 1));
                const std::uint64_t seed = FastRng::Mix64(p.rngState ^ FastRng::Mix64(call_counter));
                prng = FastRng::SplitMix64Rng(seed);
                frng = &prng;
            }
            double particle_energy_loss_eV_macro = 0.0;
            double particle_energy_loss_J_total = 0.0;
            double particle_max_nu_dt_majorant = 0.0;
            double particle_attempts_weighted = 0.0;
            double particle_null_weighted = 0.0;
            double particle_collisions_weighted = 0.0;
            std::string localErr;
            const bool ok = frng
                                ? applyParticleEventsSI(p,
                                                        w,
                                                        n_gas_m3,
                                                        dt_s,
                                                        gas_,
                                                        cfg_.use_gas_thermal_motion,
                                                        lookup_mass_kg,
                                                        po.null_op,
                                                        reactionSet,
                                                        po.enabled_channels,
                                                        *po.post,
                                                        *frng,
                                                        spawnBuffer,
                                                        &channel_counts,
                                                        &channel_counts_weighted,
                                                        &channel_energy_loss_J,
                                                        &nc.n_attempts,
                                                        &particle_attempts_weighted,
                                                        &nc.n_null,
                                                        &particle_null_weighted,
                                                        &nc.n_collisions,
                                                        &particle_collisions_weighted,
                                                        &particle_energy_loss_eV_macro,
                                                        &particle_energy_loss_J_total,
                                                        &particle_max_nu_dt_majorant,
                                                        &t_sample_s,
                                                        &t_post_s,
                                                        pname,
                                                        &localErr)
                                : applyParticleEventsSI(p,
                                                        w,
                                                        n_gas_m3,
                                                        dt_s,
                                                        gas_,
                                                        cfg_.use_gas_thermal_motion,
                                                        lookup_mass_kg,
                                                        po.null_op,
                                                        reactionSet,
                                                        po.enabled_channels,
                                                        *po.post,
                                                        rng,
                                                        spawnBuffer,
                                                        &channel_counts,
                                                        &channel_counts_weighted,
                                                        &channel_energy_loss_J,
                                                        &nc.n_attempts,
                                                        &particle_attempts_weighted,
                                                        &nc.n_null,
                                                        &particle_null_weighted,
                                                        &nc.n_collisions,
                                                        &particle_collisions_weighted,
                                                        &particle_energy_loss_eV_macro,
                                                        &particle_energy_loss_J_total,
                                                        &particle_max_nu_dt_majorant,
                                                        &t_sample_s,
                                                        &t_post_s,
                                                        pname,
                                                        &localErr);
            if (!ok)
            {
                if (err)
                    *err = localErr;
                return false;
            }

            stepStats.n_attempts += nc.n_attempts;
            stepStats.n_null += nc.n_null;
            stepStats.n_collisions += nc.n_collisions;
            stepStats.n_attempts_weighted += particle_attempts_weighted;
            stepStats.n_null_weighted += particle_null_weighted;
            stepStats.n_collisions_weighted += particle_collisions_weighted;

            pstats.n_attempts += nc.n_attempts;
            pstats.n_null += nc.n_null;
            pstats.n_collisions += nc.n_collisions;
            pstats.n_attempts_weighted += particle_attempts_weighted;
            pstats.n_null_weighted += particle_null_weighted;
            pstats.n_collisions_weighted += particle_collisions_weighted;
            stepStats.energy_loss_eV_macro += particle_energy_loss_eV_macro;
            stepStats.energy_loss_J_total += particle_energy_loss_J_total;
            pstats.energy_loss_eV_macro += particle_energy_loss_eV_macro;
            pstats.energy_loss_J_total += particle_energy_loss_J_total;
            stepStats.max_nu_dt_majorant = std::max(stepStats.max_nu_dt_majorant, particle_max_nu_dt_majorant);
            pstats.max_nu_dt_majorant = std::max(pstats.max_nu_dt_majorant, particle_max_nu_dt_majorant);

            // Persist per-particle RNG advancement.
            if (rng_mode == CollisionsRngMode::kPerParticle && frng)
            {
                p.rngState = frng->state();
            }
        }

        // Merge per-channel hit counts into step/projectile maps (diagnostics compatibility with legacy output).
        if (nch > 0)
        {
            for (std::size_t ci = 0; ci < nch; ++ci)
            {
                const std::uint64_t csum = channel_counts[ci];
                if (csum == 0)
                    continue;
                const ReactionChannel* ch = po.enabled_channels[ci];
                if (!ch)
                    continue;
                stepStats.per_channel[pname + "::" + ch->name] += csum;
                stepStats.per_channel_weighted[pname + "::" + ch->name] += channel_counts_weighted[ci];
                stepStats.per_channel_energy_loss_J[pname + "::" + ch->name] += channel_energy_loss_J[ci];
                pstats.per_channel[ch->name] += csum;
                pstats.per_channel_weighted[ch->name] += channel_counts_weighted[ci];
                pstats.per_channel_energy_loss_J[ch->name] += channel_energy_loss_J[ci];
            }
        }
    }

    // Update spawn counters for newly appended items.
    const std::size_t n1 = spawnBuffer.size();
    if (n1 > n0)
    {
        const auto t0 = Clock::now();
        stepStats.n_spawned += static_cast<std::uint64_t>(n1 - n0);
        pstats.n_spawned += static_cast<std::uint64_t>(n1 - n0);
        const auto& items = spawnBuffer.items();
        for (std::size_t i = n0; i < n1; ++i)
        {
            stepStats.spawn_per_species[items[i].species_name] += 1;
            pstats.spawn_per_species[items[i].species_name] += 1;
            const double w = (items[i].particle.w > 0.0) ? items[i].particle.w : macro_weight;
            stepStats.n_spawned_weighted += std::max(0.0, w);
            pstats.n_spawned_weighted += std::max(0.0, w);
            stepStats.spawn_per_species_weighted[items[i].species_name] += std::max(0.0, w);
            pstats.spawn_per_species_weighted[items[i].species_name] += std::max(0.0, w);
        }
        t_spawn_s += secondsBetween(t0, Clock::now());
    }

    // Accumulate timing into step stats.
    const double t_total_s = secondsBetween(tApply0, Clock::now());
    const double t_other_s = std::max(0.0, t_total_s - t_sample_s - t_post_s - t_spawn_s);
    StepStats::Timing t;
    t.total_s = t_total_s;
    t.sample_s = t_sample_s;
    t.post_s = t_post_s;
    t.spawn_count_s = t_spawn_s;
    t.other_s = t_other_s;
    stepStats.timing.Accumulate(t);
    pstats.timing.Accumulate(t);
    return true;
}

bool MCCCollisionOperator::ApplyCodeUnits(ParticleType& projectile,
                                         double macro_weight,
                                         double n_gas_m3,
                                         double dt_si,
                                         std::mt19937_64& rng,
                                         double v0_mps,
                                         double invV0,
                                         ParticleSpawnBuffer& spawnBuffer,
                                         StepStats& stepStats,
                                         std::string* err) const
{
    if (!cfg_.enabled)
        return true;

    // Guard against invalid scales (treat as identity).
    if (!(v0_mps > 0.0) || !std::isfinite(v0_mps))
    {
        v0_mps = 1.0;
        invV0 = 1.0;
    }
    if (!std::isfinite(invV0) || invV0 <= 0.0)
        invV0 = 1.0 / v0_mps;

    using Clock = std::chrono::steady_clock;
    auto secondsBetween = [](const Clock::time_point& a, const Clock::time_point& b) -> double {
        return std::chrono::duration_cast<std::chrono::duration<double>>(b - a).count();
    };
    const auto tApply0 = Clock::now();
    double t_sample_s = 0.0;
    double t_post_s = 0.0;
    double t_spawn_s = 0.0;

    const auto it = ops_.find(projectile.name());
    if (it == ops_.end())
        return true; // not configured

    if (!(dt_si >= 0.0))
    {
        if (err)
            *err = "Invalid dt_si";
        return false;
    }

    const Species* sp = species_db_.Find(projectile.name());
    if (!sp)
    {
        if (err)
            *err = "Unknown projectile species: " + projectile.name();
        return false;
    }
    const double m = sp->mass_kg;
    if (!(m > 0.0) || !std::isfinite(m))
    {
        if (err)
            *err = "Invalid projectile mass for species: " + projectile.name();
        return false;
    }

    const ProjectileOps& po = it->second;
    const ReactionSet& reactionSet = po.null_op.reaction_set();
    double lookup_mass_kg = m;
    if (!resolveLookupMassKg(gas_, reactionSet, m, &lookup_mass_kg, err))
        return false;
    const std::string& pname = projectile.name();
    StepStats::ProjectileStats& pstats = stepStats.per_projectile[pname];

    const CollisionsRngMode rng_mode = cfg_.rng_mode;
    const std::uint64_t call_counter = ++apply_counter_;
    auto ensureRngState = [&](Particle& p, std::uint64_t key) {
        if (p.rngState == 0)
        {
            p.rngState = FastRng::SeedFrom(cfg_.rng_seed, key);
        }
    };

    const std::size_t n0 = spawnBuffer.size();

#ifdef _OPENMP
    const bool can_parallel = cfg_.parallel &&
                              ((rng_mode == CollisionsRngMode::kPerParticle) || (rng_mode == CollisionsRngMode::kCounterBased));
#else
    const bool can_parallel = false;
#endif

#ifdef _OPENMP
    // Parallel particle loop (thread-safe RNG modes only).
    if (can_parallel)
    {
        PrintOMPPlacementOnce("mcc");

        struct LocalTotals
        {
            std::uint64_t n_particles = 0;
            std::uint64_t n_attempts = 0;
            std::uint64_t n_null = 0;
            std::uint64_t n_collisions = 0;
            double n_particles_weighted = 0.0;
            double n_attempts_weighted = 0.0;
            double n_null_weighted = 0.0;
            double n_collisions_weighted = 0.0;
            double energy_loss_eV_macro = 0.0;
            double energy_loss_J_total = 0.0;
            double max_nu_dt_majorant = 0.0;
            double t_sample_s = 0.0;
            double t_post_s = 0.0;
        };

        const int nthreads = std::max(1, omp_get_max_threads());
        std::vector<LocalTotals> tls_totals(static_cast<std::size_t>(nthreads));
        std::vector<ParticleSpawnBuffer> tls_spawn(static_cast<std::size_t>(nthreads));

        const std::size_t nch = po.enabled_channels.size();
        std::vector<std::vector<std::uint64_t>> tls_channel_counts(static_cast<std::size_t>(nthreads), std::vector<std::uint64_t>(nch, 0));
        std::vector<std::vector<double>> tls_channel_counts_weighted(static_cast<std::size_t>(nthreads), std::vector<double>(nch, 0.0));
        std::vector<std::vector<double>> tls_channel_energy_loss_J(static_cast<std::size_t>(nthreads), std::vector<double>(nch, 0.0));

        std::atomic<bool> had_error(false);
        std::string first_error;

#pragma omp parallel
        {
            const int tid = omp_get_thread_num();
            LocalTotals& lt = tls_totals[static_cast<std::size_t>(tid)];
            ParticleSpawnBuffer& localSpawn = tls_spawn[static_cast<std::size_t>(tid)];
            std::vector<int> accepted;
            accepted.reserve(4);

            auto& chCounts = tls_channel_counts[static_cast<std::size_t>(tid)];
            auto& chCountsWeighted = tls_channel_counts_weighted[static_cast<std::size_t>(tid)];
            auto& chEnergyLossJ = tls_channel_energy_loss_J[static_cast<std::size_t>(tid)];

#pragma omp for schedule(static)
            for (std::size_t pi = 0; pi < projectile.particle.size(); ++pi)
            {
                auto& p = projectile.particle[pi];
                const double w = (p.w > 0.0) ? p.w : macro_weight;
                if (had_error.load())
                    continue;
                lt.n_particles += 1;
                lt.n_particles_weighted += std::max(0.0, w);

                NullCollisionCounts nc;
                nc.Clear();
                FastRng::SplitMix64Rng prng(0);
                FastRng::SplitMix64Rng* frng = nullptr;
                if (rng_mode == CollisionsRngMode::kPerParticle)
                {
                    ensureRngState(p, static_cast<std::uint64_t>(pi + 1));
                    prng = FastRng::SplitMix64Rng(p.rngState);
                    frng = &prng;
                }
                else if (rng_mode == CollisionsRngMode::kCounterBased)
                {
                    ensureRngState(p, static_cast<std::uint64_t>(pi + 1));
                    const std::uint64_t seed = FastRng::Mix64(p.rngState ^ FastRng::Mix64(call_counter));
                    prng = FastRng::SplitMix64Rng(seed);
                    frng = &prng;
                }
                p.velocity *= v0_mps;
                p.vz *= v0_mps;
                std::string localErr;
                const bool ok = applyParticleEventsSI(p,
                                                      w,
                                                      n_gas_m3,
                                                      dt_si,
                                                      gas_,
                                                      cfg_.use_gas_thermal_motion,
                                                      lookup_mass_kg,
                                                      po.null_op,
                                                      reactionSet,
                                                      po.enabled_channels,
                                                      *po.post,
                                                      *frng,
                                                      localSpawn,
                                                      &chCounts,
                                                      &chCountsWeighted,
                                                      &chEnergyLossJ,
                                                      &nc.n_attempts,
                                                      &lt.n_attempts_weighted,
                                                      &nc.n_null,
                                                      &lt.n_null_weighted,
                                                      &nc.n_collisions,
                                                      &lt.n_collisions_weighted,
                                                      &lt.energy_loss_eV_macro,
                                                      &lt.energy_loss_J_total,
                                                      &lt.max_nu_dt_majorant,
                                                      &lt.t_sample_s,
                                                      &lt.t_post_s,
                                                      pname,
                                                      &localErr);
                p.velocity *= invV0;
                p.vz *= invV0;
                lt.n_attempts += nc.n_attempts;
                lt.n_null += nc.n_null;
                lt.n_collisions += nc.n_collisions;
                if (!ok)
                {
#pragma omp critical
                    {
                        if (!had_error.load())
                        {
                            had_error.store(true);
                            first_error = localErr;
                        }
                    }
                    continue;
                }

                // Persist per-particle RNG advancement.
                if (rng_mode == CollisionsRngMode::kPerParticle && frng)
                {
                    p.rngState = frng->state();
                }
            }
        } // omp parallel

        if (had_error.load())
        {
            if (err)
                *err = first_error.empty() ? std::string("MCC collision failed in parallel region") : first_error;
            return false;
        }

        // Convert spawned velocities back to code units (thread-local buffers).
        if (v0_mps != 1.0)
        {
            for (auto& sb : tls_spawn)
            {
                for (auto& it : sb.items_mut())
                {
                    it.particle.velocity *= invV0;
                    it.particle.vz *= invV0;
                }
            }
        }

        for (const auto& lt : tls_totals)
        {
            stepStats.n_particles += lt.n_particles;
            pstats.n_particles += lt.n_particles;
            stepStats.n_particles_weighted += lt.n_particles_weighted;
            pstats.n_particles_weighted += lt.n_particles_weighted;

            stepStats.n_attempts += lt.n_attempts;
            stepStats.n_null += lt.n_null;
            stepStats.n_collisions += lt.n_collisions;
            stepStats.n_attempts_weighted += lt.n_attempts_weighted;
            stepStats.n_null_weighted += lt.n_null_weighted;
            stepStats.n_collisions_weighted += lt.n_collisions_weighted;

            pstats.n_attempts += lt.n_attempts;
            pstats.n_null += lt.n_null;
            pstats.n_collisions += lt.n_collisions;
            pstats.n_attempts_weighted += lt.n_attempts_weighted;
            pstats.n_null_weighted += lt.n_null_weighted;
            pstats.n_collisions_weighted += lt.n_collisions_weighted;

            stepStats.energy_loss_eV_macro += lt.energy_loss_eV_macro;
            stepStats.energy_loss_J_total += lt.energy_loss_J_total;
            pstats.energy_loss_eV_macro += lt.energy_loss_eV_macro;
            pstats.energy_loss_J_total += lt.energy_loss_J_total;

            stepStats.max_nu_dt_majorant = std::max(stepStats.max_nu_dt_majorant, lt.max_nu_dt_majorant);
            pstats.max_nu_dt_majorant = std::max(pstats.max_nu_dt_majorant, lt.max_nu_dt_majorant);

            t_sample_s += lt.t_sample_s;
            t_post_s += lt.t_post_s;
        }

        if (nch > 0)
        {
            for (std::size_t ci = 0; ci < nch; ++ci)
            {
                std::uint64_t csum = 0;
                double csumWeighted = 0.0;
                double csumEnergyLossJ = 0.0;
                for (const auto& v : tls_channel_counts)
                    csum += v[ci];
                for (const auto& v : tls_channel_counts_weighted)
                    csumWeighted += v[ci];
                for (const auto& v : tls_channel_energy_loss_J)
                    csumEnergyLossJ += v[ci];
                if (csum == 0)
                    continue;
                const ReactionChannel* ch = po.enabled_channels[ci];
                if (!ch)
                    continue;
                stepStats.per_channel[pname + "::" + ch->name] += csum;
                stepStats.per_channel_weighted[pname + "::" + ch->name] += csumWeighted;
                stepStats.per_channel_energy_loss_J[pname + "::" + ch->name] += csumEnergyLossJ;
                pstats.per_channel[ch->name] += csum;
                pstats.per_channel_weighted[ch->name] += csumWeighted;
                pstats.per_channel_energy_loss_J[ch->name] += csumEnergyLossJ;
            }
        }

        // Merge spawned particles deterministically by thread id.
        if (!tls_spawn.empty())
        {
            auto& dst = spawnBuffer.items_mut();
            for (auto& sb : tls_spawn)
            {
                auto& src = sb.items_mut();
                dst.insert(dst.end(), std::make_move_iterator(src.begin()), std::make_move_iterator(src.end()));
                src.clear();
            }
        }
    }
    else
#endif
    {
        std::vector<int> accepted;
        accepted.reserve(4);

        const std::size_t nch = po.enabled_channels.size();
        std::vector<std::uint64_t> channel_counts(nch, 0);
        std::vector<double> channel_counts_weighted(nch, 0.0);
        std::vector<double> channel_energy_loss_J(nch, 0.0);

        for (std::size_t pi = 0; pi < projectile.particle.size(); ++pi)
        {
            auto& p = projectile.particle[pi];
            const double w = (p.w > 0.0) ? p.w : macro_weight;
            stepStats.n_particles += 1;
            pstats.n_particles += 1;
            stepStats.n_particles_weighted += std::max(0.0, w);
            pstats.n_particles_weighted += std::max(0.0, w);

            NullCollisionCounts nc;
            nc.Clear();
            FastRng::SplitMix64Rng prng(0);
            FastRng::SplitMix64Rng* frng = nullptr;

            if (rng_mode == CollisionsRngMode::kPerParticle)
            {
                ensureRngState(p, static_cast<std::uint64_t>(pi + 1));
                prng = FastRng::SplitMix64Rng(p.rngState);
                frng = &prng;
            }
            else if (rng_mode == CollisionsRngMode::kCounterBased)
            {
                ensureRngState(p, static_cast<std::uint64_t>(pi + 1));
                const std::uint64_t seed = FastRng::Mix64(p.rngState ^ FastRng::Mix64(call_counter));
                prng = FastRng::SplitMix64Rng(seed);
                frng = &prng;
            }
            double particle_energy_loss_eV_macro = 0.0;
            double particle_energy_loss_J_total = 0.0;
            double particle_max_nu_dt_majorant = 0.0;
            double particle_attempts_weighted = 0.0;
            double particle_null_weighted = 0.0;
            double particle_collisions_weighted = 0.0;
            p.velocity *= v0_mps;
            p.vz *= v0_mps;
            std::string localErr;
            const bool ok = frng
                                ? applyParticleEventsSI(p,
                                                        w,
                                                        n_gas_m3,
                                                        dt_si,
                                                        gas_,
                                                        cfg_.use_gas_thermal_motion,
                                                        lookup_mass_kg,
                                                        po.null_op,
                                                        reactionSet,
                                                        po.enabled_channels,
                                                        *po.post,
                                                        *frng,
                                                        spawnBuffer,
                                                        &channel_counts,
                                                        &channel_counts_weighted,
                                                        &channel_energy_loss_J,
                                                        &nc.n_attempts,
                                                        &particle_attempts_weighted,
                                                        &nc.n_null,
                                                        &particle_null_weighted,
                                                        &nc.n_collisions,
                                                        &particle_collisions_weighted,
                                                        &particle_energy_loss_eV_macro,
                                                        &particle_energy_loss_J_total,
                                                        &particle_max_nu_dt_majorant,
                                                        &t_sample_s,
                                                        &t_post_s,
                                                        pname,
                                                        &localErr)
                                : applyParticleEventsSI(p,
                                                        w,
                                                        n_gas_m3,
                                                        dt_si,
                                                        gas_,
                                                        cfg_.use_gas_thermal_motion,
                                                        lookup_mass_kg,
                                                        po.null_op,
                                                        reactionSet,
                                                        po.enabled_channels,
                                                        *po.post,
                                                        rng,
                                                        spawnBuffer,
                                                        &channel_counts,
                                                        &channel_counts_weighted,
                                                        &channel_energy_loss_J,
                                                        &nc.n_attempts,
                                                        &particle_attempts_weighted,
                                                        &nc.n_null,
                                                        &particle_null_weighted,
                                                        &nc.n_collisions,
                                                        &particle_collisions_weighted,
                                                        &particle_energy_loss_eV_macro,
                                                        &particle_energy_loss_J_total,
                                                        &particle_max_nu_dt_majorant,
                                                        &t_sample_s,
                                                        &t_post_s,
                                                        pname,
                                                        &localErr);
            p.velocity *= invV0;
            p.vz *= invV0;
            if (!ok)
            {
                if (err)
                    *err = localErr;
                return false;
            }

            stepStats.n_attempts += nc.n_attempts;
            stepStats.n_null += nc.n_null;
            stepStats.n_collisions += nc.n_collisions;
            stepStats.n_attempts_weighted += particle_attempts_weighted;
            stepStats.n_null_weighted += particle_null_weighted;
            stepStats.n_collisions_weighted += particle_collisions_weighted;
            pstats.n_attempts += nc.n_attempts;
            pstats.n_null += nc.n_null;
            pstats.n_collisions += nc.n_collisions;
            pstats.n_attempts_weighted += particle_attempts_weighted;
            pstats.n_null_weighted += particle_null_weighted;
            pstats.n_collisions_weighted += particle_collisions_weighted;
            stepStats.energy_loss_eV_macro += particle_energy_loss_eV_macro;
            stepStats.energy_loss_J_total += particle_energy_loss_J_total;
            pstats.energy_loss_eV_macro += particle_energy_loss_eV_macro;
            pstats.energy_loss_J_total += particle_energy_loss_J_total;
            stepStats.max_nu_dt_majorant = std::max(stepStats.max_nu_dt_majorant, particle_max_nu_dt_majorant);
            pstats.max_nu_dt_majorant = std::max(pstats.max_nu_dt_majorant, particle_max_nu_dt_majorant);

            if (rng_mode == CollisionsRngMode::kPerParticle && frng)
            {
                p.rngState = frng->state();
            }
        }

        // Merge per-channel hit counts into step/projectile maps (diagnostics compatibility with legacy output).
        if (nch > 0)
        {
            for (std::size_t ci = 0; ci < nch; ++ci)
            {
                const std::uint64_t csum = channel_counts[ci];
                if (csum == 0)
                    continue;
                const ReactionChannel* ch = po.enabled_channels[ci];
                if (!ch)
                    continue;
                stepStats.per_channel[pname + "::" + ch->name] += csum;
                stepStats.per_channel_weighted[pname + "::" + ch->name] += channel_counts_weighted[ci];
                stepStats.per_channel_energy_loss_J[pname + "::" + ch->name] += channel_energy_loss_J[ci];
                pstats.per_channel[ch->name] += csum;
                pstats.per_channel_weighted[ch->name] += channel_counts_weighted[ci];
                pstats.per_channel_energy_loss_J[ch->name] += channel_energy_loss_J[ci];
            }
        }

        // Convert spawned velocities back to code units for newly appended items.
        if (v0_mps != 1.0)
        {
            const std::size_t n1_pre = spawnBuffer.size();
            if (n1_pre > n0)
            {
                for (std::size_t i = n0; i < n1_pre; ++i)
                {
                    auto& it = spawnBuffer.items_mut()[i];
                    it.particle.velocity *= invV0;
                    it.particle.vz *= invV0;
                }
            }
        }
    }

    // Update spawn counters for newly appended items.
    const std::size_t n1 = spawnBuffer.size();
    if (n1 > n0)
    {
        const auto t0 = Clock::now();
        stepStats.n_spawned += static_cast<std::uint64_t>(n1 - n0);
        pstats.n_spawned += static_cast<std::uint64_t>(n1 - n0);
        const auto& items = spawnBuffer.items();
        for (std::size_t i = n0; i < n1; ++i)
        {
            stepStats.spawn_per_species[items[i].species_name] += 1;
            pstats.spawn_per_species[items[i].species_name] += 1;
            const double w = (items[i].particle.w > 0.0) ? items[i].particle.w : macro_weight;
            stepStats.n_spawned_weighted += std::max(0.0, w);
            pstats.n_spawned_weighted += std::max(0.0, w);
            stepStats.spawn_per_species_weighted[items[i].species_name] += std::max(0.0, w);
            pstats.spawn_per_species_weighted[items[i].species_name] += std::max(0.0, w);
        }
        t_spawn_s += secondsBetween(t0, Clock::now());
    }

    const double t_total_s = secondsBetween(tApply0, Clock::now());
    const double t_other_s = std::max(0.0, t_total_s - t_sample_s - t_post_s - t_spawn_s);
    StepStats::Timing t;
    t.total_s = t_total_s;
    t.sample_s = t_sample_s;
    t.post_s = t_post_s;
    t.spawn_count_s = t_spawn_s;
    t.other_s = t_other_s;
    stepStats.timing.Accumulate(t);
    pstats.timing.Accumulate(t);
    return true;
}
