#pragma once

#include <string>
#include <vector>

// ----------------------------------------------------------------------------
// ParticleManagementConfig
//
// Configuration for particle population control (merge/split/thinning).
// Parsed from PIC dictionary (ParticleManagement block).
// ----------------------------------------------------------------------------
struct ParticleManagementConfig
{
    bool enabled = false;
    int every = 1;

    // Global caps (all species combined). 0 disables the cap.
    long long softCapTotal = 0;
    long long hardCapTotal = 0;

    // Per-cell caps (defaults). 0 disables.
    int perCellMax = 120;
    int perCellMin = 20;

    // Merge controls
    bool mergeEnabled = true;
    double mergeMaxDeltaE_eV = 1.0;
    // Optional per-species overrides (by Species order in pic.dict).
    // If empty or shorter than species count, fallback to mergeMaxDeltaE_eV.
    std::vector<double> mergeMaxDeltaE_bySpecies_eV;
    int mergeMaxPairsPerCell = 64;

    // Split controls
    bool splitEnabled = true;
    bool splitByCountEnabled = false;
    double splitWMax = 4.0;
    int splitMaxPerCell = 8;
    double splitJitterFrac = 0.25;
    bool splitEbSafe = true;
    int splitMaxRetries = 8;

    // Spawn thinning controls (applied to mccSpawnPending / seeSpawn)
    bool thinningEnabled = true;
    bool thinningMccSpawn = true;
    bool thinningSeeSpawn = true;
    double minKeepProb = 0.05;

    // Diagnostics / report
    bool reportEnabled = true;
    int reportEvery = 1;
};
