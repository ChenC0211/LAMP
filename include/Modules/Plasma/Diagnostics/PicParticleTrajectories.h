#pragma once

#include "Modules/Plasma/Diagnostics/PicParticleTrajectoryConfig.h"
#include <filesystem>
#include <memory>
#include <random>
#include <string>
#include <vector>

class ParticleType;

// One run's trajectory batches and runtime ID sequence. Particle state and RNG
// are borrowed per call; sampling uses the simulation RNG in electron/ion order.
// Sampling cadence is based on the supplied completed-step number, not call count.
class PicParticleTrajectories
{
public:
    PicParticleTrajectories();
    ~PicParticleTrajectories(); // Releases samples; never writes or retries IO.
    PicParticleTrajectories(const PicParticleTrajectories&) = delete;
    PicParticleTrajectories& operator=(const PicParticleTrajectories&) = delete;

    // Configure once. No filesystem access or RNG use here.
    void configure(const PicParticleTrajectoryConfig& config, int seriesEvery,
        const std::filesystem::path& seriesDir, const std::string& moduleTag);
    // Also used when diagnostics are disabled; preserves existing nonzero IDs.
    // May precede configuration to retain the implicit driver's startup order.
    void assignMissingIds(std::vector<ParticleType>& species);
    // Original startup policy: clear the configured output root only if enabled,
    // including on restart. The caller supplies the intended run output folder.
    void resetOutput();
    void sample(std::vector<ParticleType>& species, std::mt19937_64& rng,
        int step, double time_s, double lengthScale_m, double velocityScale_m);
    // Explicit normal-completion output of any unfinished batch. No checkpoint
    // serialization or continuation of trajectory histories across restart.
    void finish();

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};
