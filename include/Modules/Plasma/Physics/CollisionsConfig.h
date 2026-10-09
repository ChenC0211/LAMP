#pragma once

#include <cstdint>
#include <iosfwd>
#include <string>
#include <vector>

#include "Modules/Plasma/Physics/NullCollisionOperator.h"

// -----------------------------------------------------------------------------
// CollisionsConfig
//
// Runtime configuration for the MCC collision pipeline.
//
// Design goals:
// - Core-level config object (usable by modules and solvers).
// - Minimal: defines which projectile species use which ReactionSet.
// - Reuses NullCollisionOperator::Options for the majorant settings.
// -----------------------------------------------------------------------------

struct CollisionMapping
{
    // Projectile species name (must exist in SpeciesDatabase when MCC runs).
    std::string projectile;

    // ReactionSet key registered in ReactionDatabase.
    std::string reaction_set;
};

enum class CollisionsRngMode
{
    kGlobal = 0,
    kPerParticle = 1,
    kCounterBased = 2,
};

const char* ToString(CollisionsRngMode m);
bool ParseCollisionsRngMode(const std::string& s, CollisionsRngMode* out);

struct CollisionsConfig
{
    bool enabled = false;

    // If true and built with OpenMP, MCC may parallelize the per-particle loop
    // when RNG mode is thread-safe (per-particle or counter-based).
    bool parallel = true;

    // RNG mode for stochastic operators (collisions/reactions).
    // - kGlobal      : use a shared RNG stream (legacy behavior).
    // - kPerParticle : store RNG state per particle for reproducibility under threading.
    // - kCounterBased: particle key + persistent operator application count.
    CollisionsRngMode rng_mode = CollisionsRngMode::kGlobal;

    // Global seed used to initialize per-particle RNG states when rng_mode is kPerParticle.
    std::uint64_t rng_seed = 12345;

    // Stats/diagnostics output cadence (0 disables periodic outputs; callers may still
    // force a final write at the end of a run).
    int stats_every = 1;


    // Enable optional fast-neutral tracking in PostCollisionModel for
    // ion-neutral charge-exchange channels when reactions.dict provides neutral
    // products (e.g., "Ar_fast").
    bool track_fast_neutrals = false;

    // If true, sample a thermal background-gas velocity for each MCC event and
    // use the resulting relative velocity for the null-collision process and
    // lookup energy. Default is false to preserve legacy behavior.
    bool use_gas_thermal_motion = false;

    // Currently only null-collision majorant is supported.
    NullCollisionOperator::Options null_collision;

    // Projectile -> reaction set mappings.
    std::vector<CollisionMapping> mappings;

    // Basic validation for downstream wiring.
    bool Validate(std::string* err = nullptr) const;

    void PrintTable(std::ostream& os) const;
    std::string TableString() const;
};
