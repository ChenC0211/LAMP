#include "Modules/Plasma/Boundary/LegacyGridMotionEvolution.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

LegacyGridMotionEvolution::LegacyGridMotionEvolution(
    LegacyGridMotionConfig config, LegacyGridMotion::Bindings lines,
    Reference reference, InnerPointBindings innerPoints)
    : config_(config), motion_(std::move(config), std::move(lines)),
      reference_(std::move(reference)), innerPoints_(std::move(innerPoints))
{}

LegacyGridMotionState LegacyGridMotionEvolution::prepare(int step, bool force) const
{
    auto state = motion_.state(step, true);
    auto approxSame = [](double a, double b) -> bool {
        const double scale = std::max(1.0, std::max(std::abs(a), std::abs(b)));
        return std::abs(a - b) <= 1e-12 * scale;
    };
    const bool stateChanged =
        !haveApplied_ ||
        !approxSame(state.stage1Alpha, applied_.stage1Alpha) ||
        !approxSame(state.stage2Alpha, applied_.stage2Alpha) ||
        !approxSame(state.accelShift_m, applied_.accelShift_m) ||
        !approxSame(state.decelShift_m, applied_.decelShift_m) ||
        !approxSame(state.accelRecession_m, applied_.accelRecession_m);
    const int every = std::max(1, config_.updateEvery);
    const bool cadenceHit = (step == 0) || ((step % every) == 0);
    const bool stageChanged = haveApplied_ && state.stageName != applied_.stageName;
    state.changed = stateChanged && (force || cadenceHit || stageChanged);
    if (!force && (!stateChanged || (!cadenceHit && !stageChanged)) && haveApplied_)
    {
        auto retained = applied_;
        retained.changed = false;
        return retained;
    }
    return state;
}

void LegacyGridMotionEvolution::deformReference(
    std::vector<Line>& lines, const LegacyGridMotionState& state) const
{
    lines = reference_.lines;
    motion_.deform(lines, state);
}

Vector2D LegacyGridMotionEvolution::innerPoint(
    std::size_t index, const LegacyGridMotionState& state) const
{
    auto point = reference_.innerPoints.at(index);
    const int id = static_cast<int>(index);
    if (std::find(innerPoints_.accel.begin(), innerPoints_.accel.end(), id) != innerPoints_.accel.end())
        point = point + Vector2D(state.accelShift_m, 0.0);
    if (std::find(innerPoints_.decel.begin(), innerPoints_.decel.end(), id) != innerPoints_.decel.end())
        point = point + Vector2D(state.decelShift_m, 0.0);
    return point;
}

void LegacyGridMotionEvolution::accept(LegacyGridMotionState state, const std::vector<Line>& lines)
{
    motion_.measure(lines, state);
    applied_ = std::move(state);
    haveApplied_ = true;
}

void LegacyGridMotionEvolution::restore(int nextStep, const std::vector<Line>& checkpointLines)
{
    if (nextStep < 0)
        throw std::invalid_argument("Legacy grid checkpoint has a negative completed step count");
    const int last = std::max(0, nextStep - 1);
    const int every = std::max(1, config_.updateEvery);
    int appliedStep = (last / every) * every;
    // A stage-name change forces an update even between cadence hits.
    const long long hold = std::max(0, config_.startupHoldSteps);
    const long long gapEnd = hold + std::max(0, config_.stage1Steps);
    for (long long boundary : {hold, gapEnd})
    {
        if (boundary > appliedStep && boundary <= last &&
            motion_.state(static_cast<int>(boundary - 1), true).stageName !=
            motion_.state(static_cast<int>(boundary), true).stageName)
            appliedStep = static_cast<int>(boundary);
    }
    auto restored = motion_.state(appliedStep, true);
    restored.changed = false;
    std::vector<Line> expected;
    deformReference(expected, restored);
    auto sameCoordinate = [](double a, double b) {
        return std::isfinite(a) && std::isfinite(b) &&
               std::abs(a - b) <= 1e-12 * std::max(1.0, std::max(std::abs(a), std::abs(b)));
    };
    if (expected.size() != checkpointLines.size())
        throw std::runtime_error("Legacy grid checkpoint/reference line counts differ");
    for (std::size_t i = 0; i < expected.size(); ++i)
        for (std::size_t j = 0; j < 2; ++j)
            if (!sameCoordinate(expected[i].points[j].x, checkpointLines[i].points[j].x) ||
                !sameCoordinate(expected[i].points[j].y, checkpointLines[i].points[j].y))
                throw std::runtime_error("Legacy grid checkpoint does not match the input reference and completed motion");
    accept(std::move(restored), checkpointLines);
}
