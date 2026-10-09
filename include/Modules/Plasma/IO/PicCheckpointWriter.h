#pragma once

#include <filesystem>
#include <functional>
#include <random>
#include <string>
#include <vector>

class PicCheckpointStore;
class PicCollisions;
class PicDomain;
class PicMovingGeometry;
class ParticleType;
class ImplicitPicFieldHistory;

// Shared legacy write order. Physical components retain their own encodings;
// the caller decides when a checkpoint is requested. Store/RNG are borrowed
// from the live PIC state, so no stale RNG or field snapshot is introduced.
// Failed writes throw and remain marked incomplete. Retention runs only after
// every required component has been written and the completion marker removed.
class PicCheckpointWriter
{
public:
    using Components = std::function<void(const std::filesystem::path&, const std::string&)>;
    PicCheckpointWriter(const PicCheckpointStore& store, const std::mt19937_64& rng,
        bool saveRng, int keep, std::string moduleName);
    void write(int nextStep, double time, double dt, const Components& components) const;
    struct LiveState
    {
        const PicDomain& domain;
        const PicMovingGeometry& moving;
        const std::vector<ParticleType>& species;
        const PicCollisions* collisions = nullptr;
        const ImplicitPicFieldHistory* implicitFields = nullptr;
    };
    // Borrow the current owners only for this call, including the late-created
    // collision state. Preserve mesh -> geometry -> wall -> MCC -> species ->
    // optional implicit history order inside the existing write transaction.
    void write(int nextStep, double time, double dt, const LiveState& state) const;
    // Called with the other live components; shares the saveRng gate.
    void writeCollisions(const std::filesystem::path& directory, const PicCollisions& collisions) const;

private:
    const PicCheckpointStore& store_;
    const std::mt19937_64& rng_;
    const bool saveRng_;
    const int keep_;
    const std::string moduleName_;
};
