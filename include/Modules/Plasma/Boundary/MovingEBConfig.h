#pragma once

#include <string>
#include <vector>

// Legacy geometry.dict/MovingEB physical inputs; lengths are stored in SI.
struct MovingEBImpactModel
{
    bool ionsOnly = true;
    int minSamples = 1;
    std::string smoothing = "none";
    double smoothingLength_m = 0.0;
};

struct MovingEBSelector
{
    bool wholeBoundary = false;
    std::vector<int> lineIds; // 0-based
    bool hasXRange = false;
    bool hasYRange = false;
    double xMin_m = 0.0;
    double xMax_m = 0.0;
    double yMin_m = 0.0;
    double yMax_m = 0.0;
};

struct MovingEBRegion
{
    std::string name;
    int boundaryId = -1;      // 1-based
    int levelSetIndex = -1;   // 0-based
    MovingEBSelector selector;

    std::string motionMode = "fixed";   // fixed|rigidShift|normalErosion
    std::string model = "prescribed";   // prescribed|ionErosion

    int startStep = 0;
    int endStep = -1;         // inclusive; <0 means forever

    double velocityX_mps = 0.0;
    double velocityY_mps = 0.0;
    double normalSpeed_mps = 0.0;

    double C0 = 0.0;
    double alpha = 1.0;
    std::string energyModel = "linear";
    std::string angleModel = "cosine";
    double energyRef_eV = 100.0;
    double angleExponent = 1.0;
    double smoothingLength_m = 0.0;
    double vmax_mps = 0.0;
};

struct MovingEBConfig
{
    bool enabled = false;
    int updateEvery = 1;
    int outputEvery = 200;
    bool restart = true;
    std::string metricsPath = "output/series/geometry_metrics.csv";
    std::string consistencyPath;
    std::string erosionMetricsPath = "output/series/erosion_region_metrics.csv";
    std::string impactStatsPath = "output/series/impact_wall_stats.csv";
    std::string restartPath = "moving_eb_cache.bin";
    MovingEBImpactModel impactModel;
    std::vector<MovingEBRegion> regions;
};

