#pragma once

#include <cstdint>
#include <random>
#include <string>

#include "Core/Random/FastRng.h"
#include "Modules/Plasma/Particles/Particle.h"
#include "Modules/Plasma/Physics/ReactionSet.h"
#include "Modules/Plasma/Physics/SpeciesDatabase.h"
#include "Modules/Plasma/Physics/ParticleSpawnBuffer.h"

// -----------------------------------------------------------------------------
// PostCollisionModel
//
// Minimal post-collision kinematics for common electrostatic plasma MCC
// channels:
//   - elastic: isotropic or backscatter (backward-hemisphere) scattering with constant speed
//   - excitation: isotropic scattering with energy loss = threshold
//   - ionization: creates secondary electrons and ions (minimal energy partition)
//   - charge-exchange (CX): spawns a fast neutral and resets the ion velocity
//
// Notes:
// - This is a *minimal viable* model to validate the MCC pipeline end-to-end.
// - It assumes a cold background gas (stationary target); thermal and anisotropic
//   scattering are planned extensions.
// - All kinematics are in SI; thresholds are specified in eV.
// -----------------------------------------------------------------------------

class PostCollisionModel
{
public:
    // Same legacy inference used by Apply; diagnostics must not dispatch by
    // literal channel spelling independently of the physical model.
    static std::string InferKind(const ReactionSet& set, const ReactionChannel& channel);

    struct Options
    {
        // Default scatter model used when ReactionChannel.scatter_model is empty.
        std::string default_scatter_model = "isotropic";

        // Default energy-sharing model for ionization when ReactionChannel.energy_share_model is empty.
        std::string default_ionization_share_model = "equal";

        // Background gas hints for ion-neutral slow-down models.
        //
        // If gas_species_name matches ReactionSet.target (case-insensitive) and
        // the projectile is an ion, we may apply simplified ion-neutral
        // momentum/energy transfer models without explicitly sampling or
        // tracking target neutrals.
        std::string gas_species_name = "";  // e.g. "Ar"
        double gas_temperature_K = 0.0;     // background temperature [K]
        double gas_particle_mass_kg = 0.0;  // background particle mass [kg]

        // For ion-neutral ELASTIC and CHARGE_EXCHANGE, replace the projectile
        // velocity with a sample from a Maxwellian at gas_temperature_K.
        // This provides a minimal "ion slow-down" model when fast neutrals are
        // not tracked.
        bool ion_neutral_elastic_thermalize = true;
        bool ion_neutral_cx_thermalize = true;

        // For CX: whether to spawn a fast neutral (carrying pre-CX velocity)
        // when neutral products are present.
        bool track_fast_neutrals = false;

        // Floor for energies after inelastic losses (avoid negative due to rounding).
        double min_energy_eV = 0.0;
    };

    struct TargetState
    {
        Vector2D velocity = Vector2D(0.0, 0.0);
        double vz = 0.0;
    };

    struct Outcome
    {
        std::string channel_name;
        std::string kind; // human-readable inferred kind
        std::uint32_t n_spawned = 0;
        double E_before_eV = 0.0;
        double E_after_eV = 0.0;
    };

    PostCollisionModel(const SpeciesDatabase& species_db, Options opt = Options());

    // Apply post-collision kinematics to a projectile particle.
    //
    // - 'set' provides projectile/target names (used for mass lookup and kind inference).
    // - 'ch' provides threshold/products/scatter model hints.
    // - 'projectile' is modified in-place.
    // - spawned secondaries are appended into 'spawn' (if non-null).
    //
    // Returns true on success; false on invalid inputs.
    bool Apply(const ReactionSet& set, const ReactionChannel& ch,
               Particle& projectile,
               std::mt19937_64& rng,
               const TargetState* target_state,
               ParticleSpawnBuffer* spawn,
               Outcome* out = nullptr,
               std::string* err = nullptr) const;

    // FastRng overload (for reproducible per-particle / counter-based RNG streams).
    bool Apply(const ReactionSet& set, const ReactionChannel& ch,
               Particle& projectile,
               FastRng::SplitMix64Rng& rng,
               const TargetState* target_state,
               ParticleSpawnBuffer* spawn,
               Outcome* out = nullptr,
               std::string* err = nullptr) const;

private:
    const SpeciesDatabase& species_db_;
    Options opt_;
};
