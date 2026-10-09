#pragma once

#include <string>

struct PicParticleTrajectoryConfig
{
    bool enabled = false;
    int every = 0;
    int electronEvery = 0;
    int ionEvery = 0;
    int electronSampleEvery = 1;
    int ionSampleEvery = 10;
    int electronCount = 0;
    int ionCount = 0;
    int electronMaxSamples = 0;
    int ionMaxSamples = 0;
    std::string baseDir = "diagnostics/particle_trajectories";
};
