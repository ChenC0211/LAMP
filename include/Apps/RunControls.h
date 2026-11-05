#pragma once

#include <string>
#include <cmath>
#include <limits>
#include <optional>
#include <algorithm>
#include <stdexcept>

// One main clock and one field-output schedule; no module physics is owned here.
struct RunControls
{
    struct StepRetry {
        double minDt=0.0; int maxRetries=0;
        // Zero preserves the monotone halving policy. Otherwise, after this
        // many consecutive steps without retries, propose 1.25 times the
        // accepted dt, capped by Time.dt. All physical gates still apply.
        int recoverAfter=0;
    };
    // Parsed only for an App that explicitly opts into transactional retry.
    std::optional<StepRetry> stepRetry;
    // Explicit Time.Adaptive block. This is the only owner of the selected
    // step-size policy; absence preserves the original fixed/retry behavior.
    struct Adaptive {
        double minDt=0.0,maxDt=0.0,targetCourant=0.0;
        double safety=0.9,maxGrowth=1.25;
    };
    std::optional<Adaptive> adaptive;
    // Accepted VOF End checkpoints; zero disables new publications. Restart
    // input is an explicit immutable checkpoint directory, never a field DAT.
    int checkpointEverySteps=10;
    bool checkpointExplicit=false;
    std::string restartDirectory;
    double dt = 0.0, endTime = 0.0, writeInterval = 0.0; // code time; adaptive dt is the first-step seed
    int writeEverySteps = 0, logEverySteps = 1;
    std::string outputPrefix;
    bool writeFields = true, writeInitial = true, writeFirst = true, writeFinal = true;
    bool diagnosticFields = false; // physical profile omits duplicate/native field dumps
    bool diagnosticSeries = false; // physical profile keeps process histories and essential diagnostics
    bool detailedLog = false;

    // The accepted clock is accumulated in binary floating point. Include
    // the rounding budget of prior additions/subtractions when selecting its
    // last endpoint, so a nominal full window cannot create an ulp-size step.
    // This only selects the actual next endpoint; physics validates that dt.
    double stepEnd(double begin, double proposedDt, int acceptedSteps) const
    {
        if (!std::isfinite(begin) || !std::isfinite(proposedDt) || proposedDt <= 0.
            || !std::isfinite(endTime) || !(endTime > begin) || acceptedSteps < 0)
            throw std::invalid_argument("invalid next-step clock");
        double end = std::min(endTime, begin + proposedDt);
        const double scale = std::max(std::abs(begin), std::abs(endTime));
        const double rounding = 8. * std::numeric_limits<double>::epsilon()
            * (static_cast<double>(acceptedSteps) + 1.) * scale;
        if (endTime - end <= rounding) end = endTime;
        if (!(end > begin)) throw std::runtime_error("next time is not representable");
        return end;
    }

    // Shared validation for endTime/dt input. The App selects this count once
    // for its loop and passes the final step to IO; writers do not infer it
    // from rounded physical times. Legacy integer-step inputs keep that integer.
    std::optional<int> fixedStepCount() const noexcept
    {
        if (!std::isfinite(dt) || dt <= 0.0 || !std::isfinite(endTime) || endTime <= 0.0)
            return std::nullopt;
        const double count = std::round(endTime / dt);
        if (!std::isfinite(count) || count < 1.0 || count > std::numeric_limits<int>::max()
            || std::abs(count * dt - endTime) > 1e-12 * endTime)
            return std::nullopt;
        return static_cast<int>(count);
    }

    // Select fields from the App's accepted step and explicit final step.
    // This does not advance a clock or infer completion from endTime. Time
    // intervals are anchored at zero; forced frames never shift the cadence.
    bool writesFieldsAtStep(int step, int finalStep) const noexcept
    {
        if (!writeFields || finalStep < 1 || step < 0 || step > finalStep) return false;
        if (step == 0) return writeInitial;
        if ((step == 1 && writeFirst) || (step == finalStep && writeFinal)) return true;
        if (writeEverySteps != 0)
            return writeEverySteps > 0 && step % writeEverySteps == 0;
        if (!(dt > 0.0) || !std::isfinite(dt) || !(writeInterval > 0.0) || !std::isfinite(writeInterval))
            return false;
        // At most one frame per accepted step. Avoid overflowing the interval
        // index when the requested interval is much smaller than dt.
        if (writeInterval <= dt) return true;
        // With interval > dt and an int step, these indices remain bounded by
        // INT_MAX. Match a represented interval boundary to 1e-12 intervals.
        return std::floor(step * dt / writeInterval + 1e-12)
            > std::floor((step - 1) * dt / writeInterval + 1e-12);
    }
};
