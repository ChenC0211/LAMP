#pragma once

#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

#include "Modules/Plasma/Particles/Particle.h"

// -----------------------------------------------------------------------------
// ParticleSpawnBuffer
//
// A small buffer that accumulates secondary particles produced by MCC / reaction
// models during one time step.
//
// Design notes:
// - Core-level, no dependency on specific PIC container types.
// - Stores the target species name and a Particle payload.
// - Deterministic insertion order (useful for debugging and self-tests).
// -----------------------------------------------------------------------------

struct SpawnedParticle
{
    std::string species_name;
    Particle particle;
};

class ParticleSpawnBuffer
{
public:
    ParticleSpawnBuffer() = default;

    void Clear() { items_.clear(); }
    void Reserve(std::size_t n) { items_.reserve(n); }

    std::size_t size() const { return items_.size(); }
    bool empty() const { return items_.empty(); }

    const std::vector<SpawnedParticle>& items() const { return items_; }
    std::vector<SpawnedParticle>& items_mut() { return items_; }

    // Append a spawned particle.
    void Push(const std::string& species_name, const Particle& p);

    // Pretty-print a brief summary.
    void PrintSummary(std::ostream& os, std::size_t max_items = 8) const;
    std::string SummaryString(std::size_t max_items = 8) const;

private:
    std::vector<SpawnedParticle> items_;
};
