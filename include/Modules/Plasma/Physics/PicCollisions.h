#pragma once
#include "Modules/Plasma/PlasmaCollisionInputs.h"
#include "Modules/Plasma/Physics/MCCCollisionOperator.h"
#include "Modules/Plasma/Physics/PicSpawnManagement.h"

struct SpeciesConfig;
class PicCheckpointReader;

// One initialized collision operator per PIC run. Own every object referenced
// by its majorants and post-collision models. Particle/RNG state and the explicit
// whole-step or implicit half-step schedule remain with the caller.
class PicCollisions
{
public:
    class Stage;
    // Disabled inputs return null with an empty error. Failed enabled inputs
    // return null with a diagnostic; never publish a partially built operator.
    // An optional selected checkpoint is borrowed only during construction.
    // Restore the collision RNG before returning the owner; restore failures
    // propagate from checkpoint IO. Disabled/failed creation never reads it.
    static std::unique_ptr<PicCollisions> create(PlasmaCollisionInputs inputs,
        const std::vector<SpeciesConfig>& species, double velocityScale,
        std::string& error, const PicCheckpointReader* restart = nullptr,
        bool restoreRng = true);

    PicCollisions(const PicCollisions&) = delete;
    PicCollisions& operator=(const PicCollisions&) = delete;
    PicCollisions(PicCollisions&&) = delete;
    PicCollisions& operator=(PicCollisions&&) = delete;

    // Borrowed for existing diagnostics/projectile selection. Initialization
    // cannot be repeated through this const view (including its RNG counter).
    const MCCCollisionOperator& operation() const { return *operation_; }

    using RngState = MCCCollisionOperator::RngState;
    RngState rngState() const { return operation_->rngState(); }

    bool apply(ParticleType& projectile, double macroWeight, double dtSI,
        std::mt19937_64& rng, ParticleSpawnBuffer& spawned,
        MCCCollisionOperator::StepStats& stats, std::string* error) const;

private:
    friend class PicCheckpointReader;
    PicCollisions(PlasmaCollisionInputs inputs, double velocityScale);
    // Declaration order makes the operator die before its borrowed inputs.
    const PlasmaCollisionInputs inputs_;
    SpeciesDatabase species_;
    std::unique_ptr<MCCCollisionOperator> operation_;
    const double velocityScale_;
    const double inverseVelocityScale_;
};

// Reused per run. Own collision output and statistics, while borrowing the
// initialized operator, species lookup and weights for this stage's lifetime.
// Explicit batches span all species pushes; implicit batches span one half-step.
// The caller schedules explicit acceptance, culling and batch accumulation.
// advanceImplicitHalf includes the original spawn acceptance.
class PicCollisions::Stage
{
public:
    Stage(const PicCollisions* collisions,
        const std::unordered_map<std::string, std::size_t>& speciesIndex,
        const std::vector<double>& weights)
        : collisions_(collisions), speciesIndex_(speciesIndex), weights_(weights) {}
    Stage(const Stage&) = delete;
    Stage& operator=(const Stage&) = delete;

    void beginBatch() { batchStats_.Clear(); pending_.Clear(); }
    void clearAccumulated() { accumulatedStats_.Clear(); }
    void accumulateBatch() { accumulatedStats_.Accumulate(batchStats_); }
    const MCCCollisionOperator::StepStats& batchStats() const { return batchStats_; }
    const MCCCollisionOperator::StepStats& accumulatedStats() const { return accumulatedStats_; }
    const ParticleSpawnBuffer& pending() const { return pending_; }

    // Preserve explicit unknown-species warnings and host repair before queuing.
    // False retains the old failure diagnostic; already applied work is not rolled back.
    bool applyExplicit(const QuadtreeMesh& mesh, ParticleType& projectile,
        double macroWeight, double dtSI, std::mt19937_64& rng, int step);
    // Reset only for an enabled, positive finite half-step. Visit all projectiles
    // before publishing pending spawns to the caller; false means no work, and
    // errors retain the old throw.
    bool applyImplicitHalf(std::vector<ParticleType>& species, double dtSI, std::mt19937_64& rng);
    // Complete one implicit collision half-step: process all projectiles, then
    // thin/account/commit their retained spawns. Disabled/invalid dt is a no-op,
    // preserving the previous batch and caller's statistics. No batch totals are
    // accumulated here; the integrator retains those original publication points.
    bool advanceImplicitHalf(const PicSpawnManagement& management,
        const PicSpawnManagement::ImplicitStage& stage, double dtSI,
        PicSpawnManagement::Stats* statistics = nullptr);
    // Join the private collision batch to the explicit deferred spawn stage.
    // Does not clear/accumulate collision statistics or expose mutable pending.
    PicSpawnManagement::ExplicitResult acceptExplicit(const PicSpawnManagement& management,
        const PicSpawnManagement::ExplicitStage& stage, PicSecondaryEmission& emission,
        std::vector<Particle>& secondary, PicParticleSources& sources, ParticleSpawnBuffer& source);
    PicSpawnManagement::Stats thinSpawns(const PicSpawnManagement& management,
        bool enabled, double probability, std::mt19937_64& rng,
        std::uniform_real_distribution<double>& uniform);
private:
    const PicCollisions* collisions_;
    const std::unordered_map<std::string, std::size_t>& speciesIndex_;
    const std::vector<double>& weights_;
    ParticleSpawnBuffer spawned_, pending_;
    MCCCollisionOperator::StepStats batchStats_, accumulatedStats_;
};
