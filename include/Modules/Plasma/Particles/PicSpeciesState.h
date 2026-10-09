#pragma once

#include <functional>
#include <unordered_map>
#include <unordered_set>
#include "Modules/Plasma/Particles/ParticleType.h"
#include "Modules/Plasma/ExplicitPIC/Pusher/ParticlePusherLeapfrog.h"
#include <random>
#include <optional>

class QuadtreeMesh;
class PicDomain;
class PicCheckpointReader;
class Vector2D;
struct SpeciesConfig;
struct NormalizationScales;
class CouplingWeightsCache;
template<class T> class CellField;

using PicInitialParticleExclusion = std::function<bool(const Vector2D&, int)>;

// Fill a fresh ParticleType constructed with this configuration's name/mass/
// charge. The caller retains ownership and diagnostic/species ordering.
// Draw position and all velocity components before the original exclusion
// predicate; rejected samples are not retried. Checkpoint initialization only
// prepares capacity/kappa and consumes no RNG or exclusion calls.
void initializePicSpecies(ParticleType& species, const SpeciesConfig& config,
    const NormalizationScales& scales, const QuadtreeMesh& mesh, std::mt19937_64& rng,
    bool fromCheckpoint, const PicInitialParticleExclusion& exclude = {});

// Owns one PIC session's species, configured weight/classification snapshots,
// master RNG and restored-name index. References are tied to this object;
// particle kernels and checkpoint IO borrow the live species/RNG synchronously.
// Initialization is one-shot; an exception stops the session without rollback.
class PicSpeciesState
{
public:
    using Engine = std::mt19937_64;
    using Observer = std::function<void(const SpeciesConfig&, const ParticleType&)>;
    PicSpeciesState(std::size_t capacity, unsigned int randomSeed);
    PicSpeciesState(const PicSpeciesState&) = delete;
    PicSpeciesState& operator=(const PicSpeciesState&) = delete;
    void initialize(const std::vector<SpeciesConfig>& config, const NormalizationScales& scales,
        const QuadtreeMesh& mesh, bool fromCheckpoint,
        const PicInitialParticleExclusion& exclude = {}, const Observer& observer = {});
    // Prepare fresh species or restore the checkpoint's RNG/particles, then
    // cull on the current domain using its explicit/implicit policy. Returns
    // the first loop step only on success; failure leaves partial state and
    // stops the session. No time advance, ID assignment or backward half-kick.
    // Rebuild the name index at the caller's original name-consumer stage.
    std::optional<int> initializeForRun(const std::vector<SpeciesConfig>& config,
        const NormalizationScales& scales, const PicDomain& domain,
        const PicCheckpointReader* restart, bool restoreRng);
    // Stable removal from live species, shared by startup and runtime callers.
    // Implicit survivors refresh hosts; explicit survivors keep existing hosts.
    // Optional arrays accumulate losses; no RNG or ID assignment occurs here.
    std::size_t cullEmbedded(const PicDomain& domain, const char* tag,
        const std::vector<unsigned char>* speciesMask = nullptr,
        std::vector<long long>* macroLoss = nullptr, std::vector<double>* physicalLoss = nullptr);
    // Rebuild after checkpoint metadata restoration, before name-based consumers.
    // Duplicate names retain their first index. Configured weights/classes stay
    // fixed; changes to runtime names/order require rebuilding this index again.
    void rebuildIndex();
    int electronIndex() const;
    std::vector<ParticleType>& species() { return species_; }
    const std::vector<ParticleType>& species() const { return species_; }
    Engine& rng() { return rng_; }
    const Engine& rng() const { return rng_; }
    const std::vector<double>& weights() const { return weights_; }
    const std::unordered_map<std::string,std::size_t>& index() const { return index_; }
    const std::unordered_set<std::string>& ions() const { return ions_; }
    const std::unordered_set<std::string>& electrons() const { return electrons_; }
    const std::unordered_set<std::string>& neutrals() const { return neutrals_; }

    struct PushResult {
        ExplicitPIC::Pusher::DomainImpactBuffer domain;
        ExplicitPIC::Pusher::EmbeddedImpactBuffer embedded;
    };
    using PushOptions = ExplicitPIC::Pusher::ParticlePusherLeapfrog::Options;
    // Own one immutable policy snapshot for startup and all later pushes.
    // Configuration has no particle/RNG effects and cannot be repeated. The
    // domain overload retains each integrator's original embedded-wall rules;
    // no domain, geometry or options argument is retained by reference.
    void configurePush(PushOptions options);
    void configurePush(const PicDomain& domain, bool useEmbeddedGeometry,
        bool periodicX, bool periodicY, bool absorbEmbeddedCrossing = false);
    const PushOptions& pushOptions() const;
    // Fresh v^0 -> v^-1/2 using current-position weights and E^0. This calls
    // the original half-kick with -dt; it may absorb particles but never drifts
    // positions or draws RNG. Restored, already staggered or pushed states
    // cannot be initialized again. Invalid shapes reject before any kick;
    // once kicking starts, failure stops the session without retry/rollback.
    void initializeBackwardHalfStep(QuadtreeMesh& mesh,
        const std::vector<CouplingWeightsCache>& weights, const CellField<Vector2D>& field,
        double dt, ImmersedGeometry<QuadtreeMesh>* geometry);
    struct RelabelOptions {
        bool useHint = true, useFallback = true;
        static RelabelOptions fromExplicitConfig(std::string mode, std::string fallback);
    };
    // After mesh adaptation, repair particle hosts and counts in place. Without
    // hints, the legacy fallback searches only negative labels; positive labels
    // are deliberately retained. No drift, cull, weight rebuild or RNG draw.
    void relabelAfterMeshChange(const QuadtreeMesh& mesh, const RelabelOptions& options);
    // Clear the previous step's impacts, including now-empty species. Explicit
    // PIC then interleaves each species push with its original wall/SEE work.
    void beginPush();
    const PushResult& pushSpecies(std::size_t index, QuadtreeMesh& mesh,
        const CouplingWeightsCache& weights, const CellField<Vector2D>& field,
        double dt, ImmersedGeometry<QuadtreeMesh>* geometry,
        bool recordEmbedded);
    // Implicit PIC pushes all species before consuming any wall/SEE records.
    void pushAll(QuadtreeMesh& mesh, const std::vector<CouplingWeightsCache>& weights,
        const CellField<Vector2D>& field, double dt,
        ImmersedGeometry<QuadtreeMesh>* geometry);
    const PushResult& pushResult(std::size_t index) const { return pushResults_.at(index); }
    // Results expire on beginPush(), another push of that species, or resizing
    // the species collection. A failed push is not rolled back; caller stops.
private:
    friend class ImplicitPicRefinement;
    friend class ExplicitPicRefinement;
    void preparePushStorage();
    std::vector<ParticleType> species_;
    std::vector<double> weights_;
    std::unordered_map<std::string,std::size_t> index_;
    std::unordered_set<std::string> ions_, electrons_, neutrals_;
    Engine rng_;
    bool initialized_ = false;
    bool halfStepInitializationAllowed_ = false;
    std::optional<PushOptions> pushOptions_;
    std::vector<std::vector<Vector2D>> pushGather_;
    std::vector<PushResult> pushResults_;
};
