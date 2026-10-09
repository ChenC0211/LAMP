#pragma once
#include "Modules/Plasma/ExplicitPIC/PicConfig.h"
#include "Modules/Plasma/Physics/PicSpawnManagement.h"
#include <array>
#include <optional>
#include "Modules/Plasma/Particles/Particle.h"
#include <random>
#include <unordered_map>
#include <unordered_set>

namespace ExplicitPIC { namespace Pusher { enum class DomainSide; struct DomainImpact; struct EmbeddedImpact; class DomainImpactBuffer; class EmbeddedImpactBuffer; } }

class QuadtreeMesh;
class PicEmbeddedGeometry;
class DictionaryGeometry;
class ParticleType;
class ParticleSpawnBuffer;
template<class Mesh> class ImmersedGeometry;

// SEE configuration, impact/emission stage and Gaussian sampling state.
// The caller owns RNG and pending buffers, and schedules thinning/commit.
class PicSecondaryEmission
{
public:
    enum class Mode { Explicit, Implicit };
    using Engine = std::mt19937_64;
    using Uniform = std::uniform_real_distribution<double>;
    PicSecondaryEmission(SEEConfig config, Mode mode);
    PicSecondaryEmission(const PicSecondaryEmission&) = delete;
    PicSecondaryEmission& operator=(const PicSecondaryEmission&) = delete;
    // Explicit sampling historically starts a fresh normal distribution each
    // step. Implicit sampling retains its cached variate across steps.
    void beginStep(int step = 0);
    double yield(double ionEnergyEV) const;
    int sampleCount(double expectedMacro, Engine& rng, Uniform& uniform) const;
    double sampleEnergyJ(Engine& rng, Uniform& uniform);
    std::array<double,3> sampleDirection(const std::array<double,3>& normal, Engine& rng, Uniform& uniform) const;
    std::array<double,3> sampleDirection(ExplicitPIC::Pusher::DomainSide side, Engine& rng, Uniform& uniform) const;

    struct Kinematics { double massKg, weight, velocityScale, epsilon; };
    struct DomainBox { Vector2D min, max, margin; };
    struct Spawn { Particle particle; double energyJ = 0.; };
    static bool embeddedNormal(const QuadtreeMesh& mesh, ImmersedGeometry<QuadtreeMesh>* geometry,
        const DictionaryGeometry& dictionary, const PicEmbeddedGeometry& embedded,
        int boundaryId, const ExplicitPIC::Pusher::EmbeddedImpact& impact,
        Vector2D& normal, double& gradientMagnitude);
    // Energy/direction draws occur before host lookup, including rejected
    // placements. Inputs are borrowed only for this call; no buffer is changed.
    std::optional<Spawn> spawnDomain(const QuadtreeMesh& mesh,
        const ExplicitPIC::Pusher::DomainImpact& impact, const DomainBox& box,
        const Kinematics& settings, Engine& rng, Uniform& uniform);
    std::optional<Spawn> spawnEmbedded(const QuadtreeMesh& mesh,
        const ExplicitPIC::Pusher::EmbeddedImpact& impact, const Vector2D& normal,
        const Kinematics& settings, Engine& rng, Uniform& uniform);

    // Resolve once after particle restore/index rebuilding. Length/velocity
    // scales are the original mode-specific driver values. Setup draws no RNG.
    void initialize(const std::unordered_set<int>& boundaryIds,
        const std::vector<ParticleType>& species, const std::vector<double>& weights,
        const std::unordered_map<std::string,std::size_t>& index,
        double lengthScale, double velocityScale, bool debug = false);
    bool enabled() const { return enabled_; }
    std::size_t emittedIndex() const { return emittedIndex_; }
    struct Stats {
        long long ionImpactMacroCathode = 0, emittedMacro = 0;
        double ionImpactPhysCathode = 0., ionImpactEnergyJCathode = 0.;
        double emittedPhys = 0., emittedEnergyJ = 0., emittedEnergyJ2 = 0.;
        std::array<double,4> injectedCharge{}; // Retained domain SEE charge, SI C.
    };
    const Stats& stats() const { return stats_; }
    struct BoundaryGeometry {
        const QuadtreeMesh& mesh;
        ImmersedGeometry<QuadtreeMesh>* geometry;
        const DictionaryGeometry& dictionary;
        const PicEmbeddedGeometry& embedded;
    };
    void processExplicit(const BoundaryGeometry& geometry,
        const std::vector<ParticleType>& species, std::size_t incidentIndex, double defaultWeight,
        const ExplicitPIC::Pusher::DomainImpactBuffer& domainImpacts,
        const ExplicitPIC::Pusher::EmbeddedImpactBuffer* embeddedImpacts,
        const DomainBox& box, const std::vector<int>& lineToBoundaryId,
        std::vector<Particle>& pending, Engine& rng, Uniform& uniform);
    // Implicit wall diagnostics have already resolved/sanitized weight and
    // energy. Keep their order and pass those exact values to the SEE stage.
    void processImplicitEmbedded(const BoundaryGeometry& geometry,
        const ExplicitPIC::Pusher::EmbeddedImpact& impact, int boundaryId,
        const std::vector<ParticleType>& species, std::size_t incidentIndex,
        double weight, double energyJ, ParticleSpawnBuffer& pending, Engine& rng, Uniform& uniform);
    // Preserve boundary provenance through shared population control. Do not
    // resize/reorder the explicit pending buffer outside this stage.
    PicSpawnManagement::Stats thinExplicit(const PicSpawnManagement& manager,
        std::vector<Particle>& pending, bool enabled, double probability, Engine& rng, Uniform& uniform);
    // Replace emission totals with actual retained weights/energies/charges;
    // an empty buffer contributes zero. Incident-ion totals are unaffected.
    void accountRetained(const std::vector<Particle>& pending,
        const std::vector<ParticleType>& species, const std::vector<double>& weights);

private:
    struct Burst { int count = 0; double weight = 0.; };
    Burst sampleBurst(double weight, double energyJ, int boundaryId, double gradientMagnitude,
        Engine& rng, Uniform& uniform);
    double impactEnergyJ(const Particle& particle, double mass) const;
    void recordImpact(double weight, double energyJ);
    void appendExplicit(const Spawn& spawned, int side, std::vector<Particle>& pending);
    void requireMode(Mode mode) const;
    const SEEConfig config_;
    const Mode mode_;
    const std::string yieldModel_, energyModel_, angular_;
    const double monoJ_, thermalJ_, meanJ_;
    std::normal_distribution<double> normal_{0.,1.};
    bool initialized_ = false, enabled_ = false, debug_ = false;
    std::unordered_set<int> boundaryIds_;
    std::array<bool,4> sideMask_{};
    std::size_t emittedIndex_ = static_cast<std::size_t>(-1);
    double macroWeight_ = 1., mass_ = 0., epsilon_ = 0., velocityScale_ = 0.;
    int step_ = 0, debugBudget_ = 24;
    Stats stats_;
    std::vector<int> pendingSides_; // Domain side 0..3; -1 for embedded emission.
};
