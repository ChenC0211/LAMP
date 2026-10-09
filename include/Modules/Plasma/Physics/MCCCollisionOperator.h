#pragma once

#include <cstdint>
#include <iosfwd>
#include <memory>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

#include "Modules/Plasma/Particles/ParticleType.h"
#include "Modules/Plasma/Physics/BackgroundGas.h"
#include "Modules/Plasma/Physics/CollisionsConfig.h"
#include "Modules/Plasma/Physics/CrossSectionDB.h"
#include "Modules/Plasma/Physics/NullCollisionOperator.h"
#include "Modules/Plasma/Physics/PostCollisionModel.h"
#include "Modules/Plasma/Physics/ReactionDatabase.h"
#include "Modules/Plasma/Physics/SpeciesDatabase.h"
#include "Modules/Plasma/Physics/ParticleSpawnBuffer.h"

// -----------------------------------------------------------------------------
// MCCCollisionOperator
//
// Bridge object that wires:
//   - ReactionDatabase + CrossSectionDB
//   - NullCollisionOperator (event generation)
//   - PostCollisionModel (kinematics + secondary generation)
// into a per-projectile collision operator.
//
// This file provides the scaffolding and initialization summary.
// The actual coupling into PIC loops is done at module/app layer.
// -----------------------------------------------------------------------------

class MCCCollisionOperator
{
public:
    struct Summary
    {
        std::uint32_t n_projectiles = 0;
        std::uint32_t n_channels_total = 0;
    };

    MCCCollisionOperator(const SpeciesDatabase& species_db,
                         const GasState& gas,
                         const CrossSectionDB& xs,
                         const ReactionDatabase& reactions);

    bool Initialize(const CollisionsConfig& cfg, std::string* err = nullptr);

    const CollisionsConfig& config() const { return cfg_; }
    const Summary& summary() const { return summary_; }

    bool HasProjectile(const std::string& projectile) const;

    void PrintSummary(std::ostream& os) const;
    std::string SummaryString() const;
    std::vector<std::string> ProjectileNames() const;
    std::vector<std::string> ChannelKeys() const;
    // Query the initialized, enabled channel actually used by Apply. Empty for
    // an unmapped projectile or unknown/disabled channel; consumes no RNG.
    std::string ChannelKind(const std::string& projectile, const std::string& channel) const;
    std::vector<std::string> SpawnSpeciesNames() const;

    struct MajorantInfo
    {
        std::string projectile;
        double sigma_majorant_m2 = 0.0;
    };
    std::vector<MajorantInfo> MajorantTable() const;

    // Per-step collision statistics aggregated across Apply() calls.
    struct StepStats
    {
        // Wall-clock breakdown for profiling.
        // Notes:
        // - These timings are intended for relative profiling only.
        // - "other" captures everything not explicitly attributed.
        // - When collisions are disabled, all fields remain zero.
        struct Timing
        {
            double total_s = 0.0;      // total time spent inside Apply()
            double sample_s = 0.0;     // NullCollisionOperator sampling (ApplyOneCollect)
            double post_s = 0.0;       // PostCollisionModel + per-channel bookkeeping
            double spawn_count_s = 0.0; // spawnBuffer scan (by-species counters)
            double other_s = 0.0;      // residual (=total - sample - post - spawn_count)

            void Clear()
            {
                total_s = 0.0;
                sample_s = 0.0;
                post_s = 0.0;
                spawn_count_s = 0.0;
                other_s = 0.0;
            }
            void Accumulate(const Timing& o)
            {
                total_s += o.total_s;
                sample_s += o.sample_s;
                post_s += o.post_s;
                spawn_count_s += o.spawn_count_s;
                other_s += o.other_s;
            }
        };

        struct ProjectileStats
        {
            std::uint64_t n_particles = 0;
            std::uint64_t n_attempts = 0;
            std::uint64_t n_null = 0;
            std::uint64_t n_collisions = 0;
            std::uint64_t n_spawned = 0;
            double n_particles_weighted = 0.0;
            double n_attempts_weighted = 0.0;
            double n_null_weighted = 0.0;
            double n_collisions_weighted = 0.0;
            double n_spawned_weighted = 0.0;

            // Energy loss accumulated over accepted collisions.
            double energy_loss_eV_macro = 0.0;
            double energy_loss_J_total = 0.0;

            // Diagnostic: maximum majorant probability (n_gas * v_rel * sigma_M * dt) over particles in this call.
            // This is dimensionless and should typically stay << 1 for accuracy.
            double max_nu_dt_majorant = 0.0;

            // Channel hit counts (by channel name, unique within the projectile's ReactionSet).
            std::unordered_map<std::string, std::uint64_t> per_channel;
            std::unordered_map<std::string, double> per_channel_weighted;
            std::unordered_map<std::string, double> per_channel_energy_loss_J;

            // Spawned particle counts by target species.
            std::unordered_map<std::string, std::uint64_t> spawn_per_species;
            std::unordered_map<std::string, double> spawn_per_species_weighted;

            // Profiling: wall-clock breakdown for this projectile.
            Timing timing;

            void Clear();
            void Accumulate(const ProjectileStats& other);
        };

        std::uint64_t n_particles = 0;
        std::uint64_t n_attempts = 0;
        std::uint64_t n_null = 0;
        std::uint64_t n_collisions = 0;
        std::uint64_t n_spawned = 0;
        double n_particles_weighted = 0.0;
        double n_attempts_weighted = 0.0;
        double n_null_weighted = 0.0;
        double n_collisions_weighted = 0.0;
        double n_spawned_weighted = 0.0;

        // Energy loss accumulated over accepted collisions.
        // - energy_loss_eV_macro: sum over macro-particles (eV per real particle).
        // - energy_loss_J_total  : scaled by per-particle weight (fallback: macroWeight) and converted to Joules.
        double energy_loss_eV_macro = 0.0;
        double energy_loss_J_total = 0.0;

        // Diagnostic: maximum majorant probability (n_gas * v_rel * sigma_M * dt) across all Apply() calls this step.
        double max_nu_dt_majorant = 0.0;

        // Channel hit counts (by unique key "<projectile>::<channel>").
        std::unordered_map<std::string, std::uint64_t> per_channel;
        std::unordered_map<std::string, double> per_channel_weighted;
        std::unordered_map<std::string, double> per_channel_energy_loss_J;

        // Spawned particle counts by target species.
        std::unordered_map<std::string, std::uint64_t> spawn_per_species;
        std::unordered_map<std::string, double> spawn_per_species_weighted;

        // Profiling: wall-clock breakdown aggregated across all projectiles.
        Timing timing;

        // Per-projectile breakdown.
        std::unordered_map<std::string, ProjectileStats> per_projectile;

        void Clear();
        void Accumulate(const StepStats& other);
        void PrintSummary(std::ostream& os) const;
        std::string SummaryString() const;
    };

    // Runtime state used together with the caller's global RNG and particles.
    // Calls count actual successful entries into the sampling phase, not steps.
    struct RngState
    {
        CollisionsRngMode mode;
        std::uint64_t seed;
        std::uint64_t applyCount;
    };
    RngState rngState() const { return {cfg_.rng_mode, cfg_.rng_seed, apply_counter_}; }
    // Validate configuration before changing live state; no random draws.
    void restoreRngState(const RngState& state);

    // Apply MCC collisions to a projectile species for one time step.
    //
    // This function:
    // - runs a null-collision majorant process for each particle,
    // - selects reaction channels via accept-reject,
    // - updates projectile velocity and produces secondaries via PostCollisionModel,
    // - appends secondaries to spawnBuffer (caller merges them into PIC containers),
    // - accumulates basic step statistics.
    bool Apply(ParticleType& projectile,
               double macro_weight,
               double n_gas_m3,
               double dt_s,
               std::mt19937_64& rng,
               ParticleSpawnBuffer& spawnBuffer,
               StepStats& stepStats,
               std::string* err = nullptr) const;

    // Variant for codes that store particle velocities in "code units".
    //
    // This avoids an O(Np) pre/post scaling loop by converting only the
    // particles that actually undergo collisions and by converting spawned
    // products before returning.
    //
    // Conventions:
    //   - projectile velocities are in code units on entry and on return
    //   - dt_si is SI seconds
    //   - v0_mps is the velocity scale (code * v0_mps = m/s)
    //   - invV0 is 1.0 / v0_mps
    bool ApplyCodeUnits(ParticleType& projectile,
                        double macro_weight,
                        double n_gas_m3,
                        double dt_si,
                        std::mt19937_64& rng,
                        double v0_mps,
                        double invV0,
                        ParticleSpawnBuffer& spawnBuffer,
                        StepStats& stepStats,
                        std::string* err = nullptr) const;

private:
    struct ProjectileOps
    {
        std::string projectile;
        std::string reaction_set_key;
        NullCollisionOperator null_op;
        std::unique_ptr<PostCollisionModel> post;
        // Enabled channel definitions in the same order as null_op sampling indices.
        std::vector<const ReactionChannel*> enabled_channels;
    };

    const SpeciesDatabase& species_db_;
    const GasState& gas_;
    const CrossSectionDB& xs_;
    const ReactionDatabase& reactions_;

    CollisionsConfig cfg_;
    Summary summary_;

    std::unordered_map<std::string, ProjectileOps> ops_;

    // Monotonic counter for counter-based RNG mode (incremented per Apply call).
    // Used to decorrelate streams across time steps while remaining reproducible
    // for a fixed call order.
    mutable std::uint64_t apply_counter_ = 0;
};
