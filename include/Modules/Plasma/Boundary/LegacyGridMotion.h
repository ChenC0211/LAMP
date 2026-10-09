#pragma once
#include "Core/Geometry/Line.h"
#include <string>
#include <vector>

// Retains the legacy three-grid scenario and its SI-valued input convention.
struct LegacyGridMotionConfig
{
    bool enabled = false;

    int updateEvery = 1;
    int outputEvery = 200;
    std::string metricsPath = "output/series/geometry_step.csv";

    int startupHoldSteps = 0;
    int stage1Steps = 0;
    int stage2Steps = 0;

    double gsaInitial_m = 0.0;
    double gsaHot_m = 0.0;
    double gadInitial_m = 0.0;
    double gadHot_m = 0.0;
    double daInitial_m = 0.0;
    double daFinal_m = 0.0;

    double accelShiftStage1_m = 0.0;
    double decelShiftStage1_m = 0.0;
    double accelRecessionStage2_m = 0.0;

    double accelRegionXMin_m = 0.0;
    double accelRegionXMax_m = 0.0;
    double decelRegionXMin_m = 0.0;
    double decelRegionXMax_m = 0.0;
    double apertureHalfGapInitial_m = 0.0;
    double accelDownstreamFaceX_m = 0.0;
    double recessionBandXMin_m = 0.0;
    double recessionBandXMax_m = 0.0;
    double recessionBandYMin_m = 0.0;
    double recessionBandYMax_m = 0.0;
};


struct LegacyGridMotionState
{
    std::string stageName = "static";
    double stage1Alpha = 0.0;
    double stage2Alpha = 0.0;
    double accelShift_m = 0.0;
    double decelShift_m = 0.0;
    double accelRecession_m = 0.0;
    double gsa_m = 0.0;
    double gad_m = 0.0;
    double da_m = 0.0;
    double screenDownstreamX_m = 0.0;
    double accelUpstreamX_m = 0.0;
    double accelDownstreamLocalX_m = 0.0;
    double decelUpstreamX_m = 0.0;
    bool changed = false;
};

class LegacyGridMotion
{
public:
    struct Bindings
    {
        std::vector<int> screenUpperLineIds;
        std::vector<int> screenLowerLineIds;
        std::vector<int> accelUpperLineIds;
        std::vector<int> accelLowerLineIds;
        std::vector<int> decelUpperLineIds;
        std::vector<int> decelLowerLineIds;
    };

    LegacyGridMotion(LegacyGridMotionConfig config, Bindings bindings);
    LegacyGridMotionState state(int step, bool geometryMotionEnabled) const;
    // Caller restores the reference polylines before applying absolute motion.
    void deform(std::vector<Line>& lines, const LegacyGridMotionState& state) const;
    void measure(const std::vector<Line>& lines, LegacyGridMotionState& state) const;
private:
    const LegacyGridMotionConfig config_;
    const std::vector<int> screenUpperLineIds_;
    const std::vector<int> screenLowerLineIds_;
    const std::vector<int> accelUpperLineIds_;
    const std::vector<int> accelLowerLineIds_;
    const std::vector<int> decelUpperLineIds_;
    const std::vector<int> decelLowerLineIds_;
};
