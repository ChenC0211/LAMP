#include "Modules/Plasma/PicStepSession.h"
#include <cmath>
#include <stdexcept>

PicStepResult PicStepSession::advance(bool writeCheckpoint)
{
    if (failed_ || advancing_)
        throw std::logic_error("PIC step session is stopped or already advancing");
    if (nextStep_ >= endStep_)
        throw std::out_of_range("PIC step session has no remaining macro-step");
    if (writeCheckpoint && !completion_)
        throw std::invalid_argument("PIC session has no checkpoint completion operation");
    lastStep_ = PicStepResult{{nextStep_, time(), static_cast<double>(nextStep_ + 1) * dt_, dt_},
        PicStepResult::Stage::MacroStep, 0, writeCheckpoint, false};
    advancing_ = true;
    int result;
    try
    {
        result = step_(lastStep_->time);
        if (result == 0 && completion_)
        {
            lastStep_->stage = PicStepResult::Stage::Completion;
            result = completion_(lastStep_->time, writeCheckpoint);
        }
    }
    catch (...)
    {
        advancing_ = false;
        failed_ = true;
        lastStep_->threw = true;
        throw;
    }
    advancing_ = false;
    lastStep_->status = result;
    if (result != 0)
    {
        failed_ = true;
        failureCode_ = result;
        return *lastStep_;
    }
    ++nextStep_;
    lastStep_->stage = PicStepResult::Stage::Accepted;
    return *lastStep_;
}

int drivePicStepSession(int firstStep, int endStep, double dt, const PicStepOperation& step,
    const PicSessionDriver& driver, const PicStepCompletion& completion)
{
    if (!step || !driver)
        throw std::invalid_argument("PIC step operation and session driver are required");
    if (firstStep < 0 || endStep < 0 || !std::isfinite(dt) || dt <= 0.0
        || !std::isfinite(static_cast<double>(firstStep) * dt)
        || !std::isfinite(static_cast<double>(endStep) * dt))
        throw std::invalid_argument("PIC session requires nonnegative steps and finite positive code-time range");
    PicStepSession session(firstStep, endStep, dt, step, completion);
    const int result = driver(session);
    if (session.failureCode_ != 0) return session.failureCode_;
    if (session.failed_)
        throw std::logic_error("PIC driver suppressed a failed macro-step");
    if (result != 0) return result;
    if (session.nextStep() < session.endStep())
        throw std::logic_error("PIC driver returned success before completing its macro-steps");
    return 0;
}
