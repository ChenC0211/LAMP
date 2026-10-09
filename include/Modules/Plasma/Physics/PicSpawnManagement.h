#pragma once
#include "Modules/Plasma/Particles/ParticleManagementConfig.h"
#include <array>
#include <random>
#include <unordered_map>

class Particle;
class ParticleType;
class ParticleSpawnBuffer;
class QuadtreeMesh;
class PicStepParticleBalance;
class PicSecondaryEmission;
class PicParticleSources;

// Shared population control for newly spawned particles. Settings are owned;
// species lookup/default weights must outlive this object. Buffers, current
// particles, mesh and RNG are borrowed only during calls. No implicit clear,
// culling or change to the caller's commit/accounting stage.
class PicSpawnManagement
{
public:
    enum class Mode { Explicit, Implicit };
    struct Stats {
        long long in = 0, kept = 0, drop = 0;
        double physIn = 0., physKept = 0., keepProbWeightedMacro = 0.;
        long long keepProbMacroCount = 0;
        void accumulate(const Stats& other);
        double reportKeepProb() const {
            return keepProbMacroCount > 0 ? keepProbWeightedMacro / static_cast<double>(keepProbMacroCount) : 1.;
        }
    };
    PicSpawnManagement(Mode mode, ParticleManagementConfig config,
        const std::unordered_map<std::string,std::size_t>& speciesIndex,
        const std::vector<double>& defaultWeights);
    double keepProbability(const std::vector<ParticleType>& species, long long spawnCount, bool enabled) const;
    Stats thin(ParticleSpawnBuffer& buffer, bool enabled, double probability,
        std::mt19937_64& rng, std::uniform_real_distribution<double>& uniform,
        std::vector<int>* tags = nullptr) const;
    // Optional per-particle tags are compacted with the same stable selection.
    // Size mismatches fail before modifying particles or consuming RNG.
    Stats thin(std::vector<Particle>& particles, std::size_t speciesIndex, bool enabled, double probability,
        std::mt19937_64& rng, std::uniform_real_distribution<double>& uniform,
        std::vector<int>* tags = nullptr) const;
    // Explicit appends prepared particles and marks touched species. Implicit
    // repairs hosts/weights and retains the MCC-only warning/sanity checks.
    // Neither path clears the buffer used by subsequent diagnostic stages.
    void commit(const ParticleSpawnBuffer& buffer, const QuadtreeMesh& mesh,
        std::vector<ParticleType>& species, std::vector<unsigned char>* touched = nullptr, bool mcc = false) const;
    void commitExplicitSee(const std::vector<Particle>& particles, std::size_t speciesIndex,
        std::vector<ParticleType>& species, std::vector<unsigned char>& touched) const;
    struct ExplicitStage {
        const QuadtreeMesh& mesh;
        std::vector<ParticleType>& species;
        std::mt19937_64& rng;
        std::uniform_real_distribution<double>& uniform;
        double velocityScale;
        std::array<int,4> domainBoundaryIds;
        bool managementActive;
    };
    struct ExplicitResult {
        Stats collision, secondary, source;
        double keepProbability = 1., injectionEnergyJ = 0.;
        std::array<double,4> sourceChargeC{};
        std::vector<unsigned char> touchedSpecies;
    };
    // One probability from all pending groups before any thinning or commit.
    // Thin MCC/SEE/source (including provenance), account SEE/source, then
    // commit MCC/source/SEE. Buffers remain available for diagnostics; the
    // caller culls touched species afterward. No borrowed inputs are retained.
    ExplicitResult acceptExplicit(const ExplicitStage& stage, ParticleSpawnBuffer& collision,
        PicSecondaryEmission& emission, std::vector<Particle>& secondary,
        PicParticleSources& sources, ParticleSpawnBuffer& source) const;
    enum class ImplicitKind { Collision, Secondary, Source };
    // Borrow the live mesh/container and current macro-step accounting only.
    // No references are retained by acceptImplicit or collision half-steps.
    struct ImplicitStage {
        const QuadtreeMesh& mesh;
        std::vector<ParticleType>& species;
        std::mt19937_64& rng;
        std::uniform_real_distribution<double>& uniform;
        PicStepParticleBalance& balance;
        double velocityScale;
        bool managementActive;
    };
    // Implicit-only thinning/accounting/commit. Collision stats accumulate;
    // Secondary/Source stats replace. Empty Secondary is a complete no-op.
    // Collision/Secondary account before commit; Source accounts afterward.
    // Secondary enters source balance only after commit. Exceptions retain all
    // already completed mutations/accounting. The retained buffer stays intact.
    void acceptImplicit(const ImplicitStage& stage, ParticleSpawnBuffer& buffer,
        ImplicitKind kind, Stats* statistics = nullptr) const;
private:
    double defaultWeight(const std::string& name) const;
    const Mode mode_;
    const ParticleManagementConfig config_;
    const std::unordered_map<std::string,std::size_t>& speciesIndex_;
    const std::vector<double>& weights_;
};
