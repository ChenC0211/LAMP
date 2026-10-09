#pragma once

#include "Modules/Plasma/IO/PicCheckpointStore.h"
#include <optional>

class ParticleType;
class QuadtreeMesh;
class PicCollisions;
struct CollisionsConfig;

// One selected legacy checkpoint. Owns only its location/metadata; restores
// into the caller's live state at the original initialization stages.
class PicCheckpointReader
{
public:
    enum class MeshCache { Optional, Required };
    class PreparedExplicitRestart;

    static std::optional<PicCheckpointReader> select(const PicCheckpointStore& store,
        int step, double& dt, std::string moduleName);
    // Select once before explicit output preparation. Keep the permissive output
    // metadata reader and collision validation here; read full restart metadata
    // and apply its dt override later, from this same directory.
    static std::optional<PreparedExplicitRestart> prepareExplicitRestart(const PicCheckpointStore& store,
        int step, const CollisionsConfig* collisions, bool restoreRng);
    const std::filesystem::path& directory() const noexcept { return directory_; }
    const PicCheckpointStore::Metadata& metadata() const noexcept { return metadata_; }
    int nextStep() const noexcept;

    bool readMesh(QuadtreeMesh& mesh, MeshCache policy) const;
    // Preserve RNG -> species caches -> relabel order. Missing later caches
    // leave earlier restores in place; the caller must stop on failure.
    bool readParticles(const QuadtreeMesh& mesh, std::vector<ParticleType>& species,
        std::mt19937_64& rng, bool restoreRng) const;
    void report() const;
    // Explicit output opens before the collision owner; validate first so an
    // unsupported restart cannot truncate existing output history.
    static void validateCollisions(const std::filesystem::path& directory,
        const CollisionsConfig& config, bool restoreRng);

private:
    friend class PicCollisions;
    // Used only by the collision factory before publishing its owner. Old
    // global/per-particle checkpoints need no counter; counterBased requires it.
    void readCollisions(PicCollisions& collisions, bool restoreRng) const;
    static std::optional<PicCheckpointReader> readSelected(std::filesystem::path directory,
        double& dt, std::string moduleName);
    PicCheckpointReader(std::filesystem::path directory, PicCheckpointStore::Metadata metadata,
        std::string moduleName);
    std::filesystem::path directory_;
    PicCheckpointStore::Metadata metadata_;
    std::string moduleName_;
};

// Pins directory identity across output preparation and physical startup. It
// neither keeps the store alive nor locks/copies checkpoint contents. A missing
// or newly incomplete selected checkpoint fails; never fall back to a newer one.
class PicCheckpointReader::PreparedExplicitRestart
{
public:
    int outputResumeStep() const noexcept { return outputResumeStep_; }
    const std::filesystem::path& directory() const noexcept { return directory_; }
    std::optional<PicCheckpointReader> read(double& dt) const;
private:
    friend class PicCheckpointReader;
    PreparedExplicitRestart(std::filesystem::path directory, int outputResumeStep);
    std::filesystem::path directory_;
    int outputResumeStep_;
};
