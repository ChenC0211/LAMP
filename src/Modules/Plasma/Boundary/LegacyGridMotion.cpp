#include "Modules/Plasma/Boundary/LegacyGridMotion.h"
#include <algorithm>
#include <cmath>
#include <utility>

LegacyGridMotion::LegacyGridMotion(LegacyGridMotionConfig config, Bindings bindings)
    : config_(std::move(config)),
      screenUpperLineIds_(std::move(bindings.screenUpperLineIds)),
      screenLowerLineIds_(std::move(bindings.screenLowerLineIds)),
      accelUpperLineIds_(std::move(bindings.accelUpperLineIds)),
      accelLowerLineIds_(std::move(bindings.accelLowerLineIds)),
      decelUpperLineIds_(std::move(bindings.decelUpperLineIds)),
      decelLowerLineIds_(std::move(bindings.decelLowerLineIds))
{}

LegacyGridMotionState LegacyGridMotion::state(int step, bool geometryMotionEnabled) const
{
    LegacyGridMotionState state;
    if (!geometryMotionEnabled)
        return state;

    const auto& ms = config_;
    const int hold = std::max(0, ms.startupHoldSteps);
    const int stage1 = std::max(0, ms.stage1Steps);
    const int stage2 = std::max(0, ms.stage2Steps);
    auto progress01 = [](int numer, int denom) -> double {
        if (denom <= 0)
            return 1.0;
        return std::clamp(static_cast<double>(numer) / static_cast<double>(denom), 0.0, 1.0);
    };

    if (step < hold)
    {
        state.stageName = "stage0_hold";
    }
    else if (step < hold + stage1)
    {
        state.stageName = "stage1_hot_gap";
        state.stage1Alpha = progress01(step - hold + 1, stage1);
    }
    else
    {
        state.stageName = (stage2 > 0) ? "stage2_recession" : "stage1_hot_gap";
        state.stage1Alpha = 1.0;
        state.stage2Alpha = progress01(step - hold - stage1 + 1, stage2);
    }

    state.accelShift_m = ms.accelShiftStage1_m * state.stage1Alpha;
    state.decelShift_m = ms.decelShiftStage1_m * state.stage1Alpha;
    state.accelRecession_m = ms.accelRecessionStage2_m * state.stage2Alpha;

    state.gsa_m = ms.gsaInitial_m + (ms.gsaHot_m - ms.gsaInitial_m) * state.stage1Alpha;
    state.gad_m = ms.gadInitial_m + (ms.gadHot_m - ms.gadInitial_m) * state.stage1Alpha;
    state.da_m = ms.daInitial_m + (ms.daFinal_m - ms.daInitial_m) * state.stage2Alpha;
    return state;

}

void LegacyGridMotion::deform(std::vector<Line>& lines, const LegacyGridMotionState& state) const
{
    auto translateBoundary = [&](const std::vector<int>& lineIds, double dx, double dy) {
        for (int lineId : lineIds)
        {
            if (lineId < 0 || lineId >= static_cast<int>(lines.size()))
                continue;
            auto& ln = lines[static_cast<std::size_t>(lineId)];
            for (auto& pt : ln.points)
            {
                pt.x += dx;
                pt.y += dy;
            }
        }
    };

    translateBoundary(accelUpperLineIds_, state.accelShift_m, 0.0);
    translateBoundary(accelLowerLineIds_, state.accelShift_m, 0.0);
    translateBoundary(decelUpperLineIds_, state.decelShift_m, 0.0);
    translateBoundary(decelLowerLineIds_, state.decelShift_m, 0.0);


    auto applyAccelRecession = [&](const std::vector<int>& lineIds) {
        if (!(state.accelRecession_m > 0.0))
            return;
        const auto& ms = config_;
        const double bandXMin = ms.recessionBandXMin_m + state.accelShift_m;
        const double bandXMax = ms.recessionBandXMax_m + state.accelShift_m;
        const double faceX = ms.accelDownstreamFaceX_m + state.accelShift_m;
        const double bandYMin = std::max(0.0, ms.recessionBandYMin_m);
        const double bandYMax = std::max(bandYMin, ms.recessionBandYMax_m);
        const double spanX = std::max(1e-15, bandXMax - bandXMin);
        const double faceTol = std::max(1e-12, 1e-6 * std::abs(faceX));

        for (int lineId : lineIds)
        {
            if (lineId < 0 || lineId >= static_cast<int>(lines.size()))
                continue;
            auto& ln = lines[static_cast<std::size_t>(lineId)];
            for (auto& pt : ln.points)
            {
                const double ay = std::abs(pt.y);
                if ((ay + 1e-15) < bandYMin || (ay - 1e-15) > bandYMax)
                    continue;
                const double xi = std::clamp((pt.x - bandXMin) / spanX, 0.0, 1.0);
                if (!(xi > 0.0))
                    continue;
                const double ramp = std::sqrt(xi);
                const double signY = (pt.y >= 0.0) ? 1.0 : -1.0;
                pt.y += signY * state.accelRecession_m * ramp;
                if (std::abs(pt.x - faceX) <= faceTol)
                    pt.x += state.accelRecession_m * ramp;
            }
        }
    };

    applyAccelRecession(accelUpperLineIds_);
    applyAccelRecession(accelLowerLineIds_);


}

void LegacyGridMotion::measure(const std::vector<Line>& lines, LegacyGridMotionState& state) const
{
    auto collectEndpoints = [&](const std::vector<int>& lineIds) {
        std::vector<Vector2D> pts;
        for (int lineId : lineIds)
        {
            if (lineId < 0 || lineId >= static_cast<int>(lines.size()))
                continue;
            const auto& ln = lines[static_cast<std::size_t>(lineId)];
            pts.push_back(ln.points[0]);
            pts.push_back(ln.points[1]);
        }
        return pts;
    };
    auto minX = [](const std::vector<Vector2D>& pts, double fallback) {
        double out = fallback;
        bool first = true;
        for (const auto& p : pts)
        {
            if (first || p.x < out)
            {
                out = p.x;
                first = false;
            }
        }
        return out;
    };
    auto maxX = [](const std::vector<Vector2D>& pts, double fallback) {
        double out = fallback;
        bool first = true;
        for (const auto& p : pts)
        {
            if (first || p.x > out)
            {
                out = p.x;
                first = false;
            }
        }
        return out;
    };
    auto maxXWithinAbsY = [](const std::vector<Vector2D>& pts, double yMax, double fallback) {
        double out = fallback;
        bool first = true;
        for (const auto& p : pts)
        {
            if (std::abs(p.y) > yMax)
                continue;
            if (first || p.x > out)
            {
                out = p.x;
                first = false;
            }
        }
        return out;
    };
    auto minXWithinAbsY = [](const std::vector<Vector2D>& pts, double yMax, double fallback) {
        double out = fallback;
        bool first = true;
        for (const auto& p : pts)
        {
            if (std::abs(p.y) > yMax)
                continue;
            if (first || p.x < out)
            {
                out = p.x;
                first = false;
            }
        }
        return out;
    };
    auto minAbsY = [](const std::vector<Vector2D>& pts, double yCap, double fallback) {
        double out = fallback;
        bool first = true;
        for (const auto& p : pts)
        {
            const double ay = std::abs(p.y);
            if (ay > yCap)
                continue;
            if (first || ay < out)
            {
                out = ay;
                first = false;
            }
        }
        return out;
    };

    std::vector<Vector2D> screenPts = collectEndpoints(screenUpperLineIds_);
    {
        auto lower = collectEndpoints(screenLowerLineIds_);
        screenPts.insert(screenPts.end(), lower.begin(), lower.end());
    }
    std::vector<Vector2D> accelPts = collectEndpoints(accelUpperLineIds_);
    {
        auto lower = collectEndpoints(accelLowerLineIds_);
        accelPts.insert(accelPts.end(), lower.begin(), lower.end());
    }
    std::vector<Vector2D> decelPts = collectEndpoints(decelUpperLineIds_);
    {
        auto lower = collectEndpoints(decelLowerLineIds_);
        decelPts.insert(decelPts.end(), lower.begin(), lower.end());
    }

    const double screenX = maxX(screenPts, config_.gsaInitial_m);
    const double accelUpX = minX(accelPts, screenX + config_.gsaInitial_m);
    const double accelDownX = maxXWithinAbsY(
        accelPts,
        std::max(0.0, config_.recessionBandYMax_m) + 1e-12,
        maxX(accelPts, accelUpX));
    const double decelUpX = minXWithinAbsY(
        decelPts,
        std::max(0.0, config_.recessionBandYMax_m) + 1e-12,
        minX(decelPts, accelDownX + config_.gadInitial_m));
    const double accelHalfGap = minAbsY(accelPts, 0.8 * 0.001, 0.5 * config_.daInitial_m);

    state.screenDownstreamX_m = screenX;
    state.accelUpstreamX_m = accelUpX;
    state.accelDownstreamLocalX_m = accelDownX;
    state.decelUpstreamX_m = decelUpX;
    state.gsa_m = accelUpX - screenX;
    state.gad_m = decelUpX - accelDownX;
    state.da_m = 2.0 * accelHalfGap;

}
