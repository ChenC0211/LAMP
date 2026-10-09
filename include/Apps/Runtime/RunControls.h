#pragma once

#include <string>
#include <cmath>
#include <limits>
#include <optional>

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
    // Optional standalone conservative Hydro stopping quantity. The App
    // binds this value to its actual column geometry and accepted state.
    std::optional<double> columnAmplitudeLimit;
    double dt = 0.0, endTime = 0.0, writeInterval = 0.0; // code time
    int writeEverySteps = 0, logEverySteps = 1;
    std::string outputPrefix;
    bool writeFields = true, writeInitial = true, writeFirst = true, writeFinal = true;
    bool diagnosticFields = false; // physical profile omits duplicate/native field dumps
    bool detailedLog = false;

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
