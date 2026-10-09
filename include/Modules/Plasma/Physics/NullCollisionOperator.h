#pragma once

#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

#include "Core/Random/FastRng.h"
#include "Modules/Plasma/Physics/CrossSectionDB.h"
#include "Modules/Plasma/Physics/ReactionSet.h"

class TabulatedCrossSection;

// Null-collision (majorant) Monte Carlo collision operator.
//
// Conventions (SI):
//   - n_gas in 1/m^3
//   - v_rel in m/s
//   - dt in seconds
//   - Cross sections in m^2 (CrossSectionDB convention)
//   - Channel sampling uses sigma_i(E)/sigma_tot(E), where E is in eV.

struct NullCollisionStats
{
    std::uint64_t n_particles = 0;
    std::uint64_t n_steps = 0;
    std::uint64_t n_attempts = 0;
    std::uint64_t n_null = 0;
    std::uint64_t n_collisions = 0;
    std::unordered_map<std::string, std::uint64_t> per_channel;

    void Clear();
    void Print(std::ostream& os) const;
    std::string ToString() const;
};

// Lightweight counters used in hot loops (avoids unordered_map updates).
struct NullCollisionCounts
{
    std::uint64_t n_particles = 0;
    std::uint64_t n_steps = 0;
    std::uint64_t n_attempts = 0;
    std::uint64_t n_null = 0;
    std::uint64_t n_collisions = 0;

    void Clear()
    {
        n_particles = 0;
        n_steps = 0;
        n_attempts = 0;
        n_null = 0;
        n_collisions = 0;
    }

    void Accumulate(const NullCollisionCounts& o)
    {
        n_particles += o.n_particles;
        n_steps += o.n_steps;
        n_attempts += o.n_attempts;
        n_null += o.n_null;
        n_collisions += o.n_collisions;
    }
};


class NullCollisionOperator
{
public:
    struct NextEvent
    {
        enum class Kind
        {
            kNone = 0,
            kNull = 1,
            kCollision = 2,
        };

        Kind kind = Kind::kNone;
        double elapsed_s = 0.0;
        int channel_index = -1;
    };

    struct Options
    {
        // Energy range used to compute a safe majorant sigma.
        // Particle energies above this range will be clamped (if enabled).
        double majorant_Emax_eV = 1000.0;

        // Multiplier applied to the sampled max sigma_tot to provide a safety margin.
        double safety_factor = 1.2;

        // Number of energy samples used to estimate max sigma_tot on [0, Emax].
        std::size_t energy_samples = 2048;

        // If true, energies passed to ApplyOne() above majorant_Emax_eV are clamped.
        bool clamp_energy = true;

        // Use fast lookup hints for tabulated cross sections in hot loops.
        //
        // When enabled, each TabulatedCrossSection channel builds a uniform-energy hint table
        // that approximates the lower_bound index, allowing SigmaM2WithHint() to avoid the
        // per-call binary search. This preserves the exact interpolation logic and should
        // produce identical results to the legacy path.
        bool use_fast_lookup = true;

        // Number of uniform energy bins used for the hint table on [0, majorant_Emax_eV].
        // Larger values reduce the local correction walk at the cost of a small precompute.
        std::size_t lookup_bins = 4096;
    };

    NullCollisionOperator() = default;

    // Initialize from a reaction set and cross section database.
    // Stores pointers to enabled channel cross-section models.
    bool Initialize(const ReactionSet& set, const CrossSectionDB& xs_db,
                    const Options& opt, std::string* err = nullptr);

    const ReactionSet& reaction_set() const { return set_; }
    double sigma_majorant_m2() const { return sigma_majorant_m2_; }
    const Options& options() const { return opt_; }

    // Apply collisions for a single particle over a time interval dt.
    // Returns the number of accepted (real) collisions.
    //
    // E_eV should represent the relative kinetic energy of the projectile.
    // v_rel is the projectile-target relative speed.
    std::uint32_t ApplyOne(double n_gas_m3, double v_rel_mps, double E_eV, double dt_s,
                           std::mt19937_64& rng, NullCollisionStats* stats = nullptr) const;

    // FastRng overloads for reproducible per-particle / counter-based RNG streams.
    std::uint32_t ApplyOne(double n_gas_m3, double v_rel_mps, double E_eV, double dt_s,
                           FastRng::SplitMix64Rng& rng, NullCollisionStats* stats = nullptr) const;

    // Variant that also records the chosen channel indices for each accepted collision event.
    // The indices correspond to the enabled channels in reaction_set().channels (after filtering).
    std::uint32_t ApplyOneCollect(double n_gas_m3, double v_rel_mps, double E_eV, double dt_s,
                                  std::mt19937_64& rng,
                                  std::vector<int>* accepted_channel_indices,
                                  NullCollisionStats* stats = nullptr) const;

    std::uint32_t ApplyOneCollect(double n_gas_m3, double v_rel_mps, double E_eV, double dt_s,
                                  FastRng::SplitMix64Rng& rng,
                                  std::vector<int>* accepted_channel_indices,
                                  NullCollisionStats* stats = nullptr) const;


// Hot-loop variant: collect accepted channel indices and lightweight counters (no per-channel map).
std::uint32_t ApplyOneCollectCounts(double n_gas_m3, double v_rel_mps, double E_eV, double dt_s,
                                    std::mt19937_64& rng,
                                    std::vector<int>* accepted_channel_indices,
                                    NullCollisionCounts* counts = nullptr) const;

std::uint32_t ApplyOneCollectCounts(double n_gas_m3, double v_rel_mps, double E_eV, double dt_s,
                                    FastRng::SplitMix64Rng& rng,
                                    std::vector<int>* accepted_channel_indices,
                                    NullCollisionCounts* counts = nullptr) const;

    // Sample at most one majorant-process event within dt_s.
    //
    // This is the event-by-event helper used by MCCCollisionOperator when the
    // projectile state must be updated after every accepted collision.
    // Unlike ApplyOne*(), this helper does not increment n_particles/n_steps in
    // NullCollisionCounts because it is intended to be called repeatedly within
    // one particle step.
    void SampleNextEvent(double n_gas_m3, double v_rel_mps, double E_eV, double dt_s,
                         std::mt19937_64& rng,
                         NextEvent* event,
                         NullCollisionCounts* counts = nullptr) const;

    void SampleNextEvent(double n_gas_m3, double v_rel_mps, double E_eV, double dt_s,
                         FastRng::SplitMix64Rng& rng,
                         NextEvent* event,
                         NullCollisionCounts* counts = nullptr) const;

    // Helper: compute sigma_tot(E) in m^2.
    double SigmaTotalM2(double E_eV) const;

private:
    struct ChannelRuntime
    {
        std::string name;
        double threshold_eV = 0.0;
        CrossSectionModelPtr model;

        // Optional fast-path for TabulatedCrossSection: a precomputed hint table that maps
        // uniform energy bins to an initial lower_bound index guess.
        const TabulatedCrossSection* tabulated = nullptr;
        std::vector<std::uint32_t> hint_i1;  // lower_bound index i1 for each energy bin
        double hint_inv_dE = 0.0;
        std::size_t hint_bins = 0;
        double hint_Emax = 0.0;

        inline bool HasHintTable() const
        {
            return (tabulated != nullptr) && !hint_i1.empty() && (hint_bins > 0) && (hint_inv_dE > 0.0);
        }

        inline std::size_t HintI1ForEnergy(double E_eV) const
        {
            if (!HasHintTable())
                return 0;
            if (E_eV <= 0.0)
                return static_cast<std::size_t>(hint_i1.front());
            if (E_eV >= hint_Emax)
                return static_cast<std::size_t>(hint_i1.back());
            const std::size_t b = std::min(hint_bins, static_cast<std::size_t>(E_eV * hint_inv_dE));
            return static_cast<std::size_t>(hint_i1[b]);
        }
    };

    double ClampEnergy(double E_eV) const;
    bool RecomputeMajorant(std::string* err);
    int SampleChannelIndex(double E_eV, std::mt19937_64& rng) const;
    int SampleChannelIndex(double E_eV, FastRng::SplitMix64Rng& rng) const;

    ReactionSet set_;
    Options opt_;
    std::vector<ChannelRuntime> channels_;
    double sigma_majorant_m2_ = 0.0;
};
