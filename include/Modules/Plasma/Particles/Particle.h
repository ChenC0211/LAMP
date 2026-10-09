#pragma once

#include <cstdint>
#include <iostream>
#include <vector>
#include "Core/Base/Base.h"
#include "Core/Mesh/Cell.h"
#include "Core/Base/Vector2D.h"
#include "Core/Mesh/Vert.h"

class Particle
{
public:
    Particle();
    Particle(double& px, double& py, double& vx, double& vy, int inCellL);
    Particle(const Particle& particle);
    ~Particle();

    Vector2D position;

    Vector2D velocity;

    // Out-of-plane velocity component for 2D geometry / 3D velocity-space.
    // In the current PIC mode, the field and pusher are 2D, so vz is unchanged by E.
    double vz = 0.0;

    int inCellL = -1;
    int inBackgroundCellL = -1;
    int discreteVertL = -1;

    // Physical macro-particle weight (number of real particles represented).
    double w = 1.0;

    double weight[4] = {0.0, 0.0, 0.0, 0.0};

    double emitChance = 0.0;

    std::vector<double> levelSetValues;

    // Per-particle RNG state for reproducible stochastic operators (e.g., MCC collisions).
    // 0 means uninitialized (caller may initialize deterministically from a global seed + particle index).
    std::uint64_t rngState = 0;

    // Stable runtime identifier used by diagnostics that need to follow the same
    // macro-particle across vector compaction and wall-loss filtering.
    std::uint64_t trackId = 0;

    template<typename Mesh>
    void updateWeight(Mesh& mesh, Cell& interpolationCell);
    template<typename Mesh, typename Field>
    void getDiscreteVert(Mesh& mesh, Field& LSField, const int levelSetI);

    void ensureLevelSetSize(std::size_t n)
    {
        if (levelSetValues.size() != n)
            levelSetValues.assign(n, 0.0);
    }
};

#include "Modules/Plasma/Particles/Particle_impl.h"
