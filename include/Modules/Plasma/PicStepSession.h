#pragma once

#include <functional>
#include <optional>

class PicStepSession;

// One macro-step in code units, formed from the restored step index and dt.
// begin/end use index multiplication, never accumulated floating-point time.
struct PicStepTime
{
    int index;
    double begin, end, dt;
};

struct PicStepResult
{
    enum class Stage { MacroStep, Completion, Accepted };
    PicStepTime time;
    Stage stage = Stage::MacroStep;
    int status = 0;
    bool checkpointRequested = false;
    bool threw = false;
    bool accepted() const noexcept { return stage == Stage::Accepted; }
};

using PicStepOperation = std::function<int(const PicStepTime&)>;
using PicStepCompletion = std::function<int(const PicStepTime&, bool)>;
using PicSessionDriver = std::function<int(PicStepSession&)>;

// Synchronous, non-owning access to one initialized PIC state. The caller
// controls the macro-step loop; the module retains particle substeps and the
// existing per-step output/checkpoint order. Do not retain this reference.
class PicStepSession
{
public:
    PicStepSession(const PicStepSession&) = delete;
    PicStepSession& operator=(const PicStepSession&) = delete;

    int nextStep() const noexcept { return nextStep_; }
    int endStep() const noexcept { return endStep_; }
    double dt() const noexcept { return dt_; }
    double time() const noexcept { return static_cast<double>(nextStep_) * dt_; }
    bool stopped() const noexcept { return failed_; }
    // A value copy remains valid across later steps. An exception records its
    // attempted phase with threw=true; status alone never certifies success.
    std::optional<PicStepResult> lastStep() const { return lastStep_; }
    // Execute exactly one complete original macro-step. A nonzero result or
    // exception stops the session; a partially changed state is never retried.
    // The caller supplies the checkpoint decision. Completion runs after the
    // physical step, preserving the old checkpoint/log order before acceptance.
    PicStepResult advance(bool writeCheckpoint = false);

private:
    friend int drivePicStepSession(int, int, double, const PicStepOperation&, const PicSessionDriver&,
        const PicStepCompletion&);
    PicStepSession(int firstStep, int endStep, double dt, const PicStepOperation& step,
        const PicStepCompletion& completion)
        : nextStep_(firstStep), endStep_(endStep), dt_(dt), step_(step), completion_(completion) {}

    int nextStep_, endStep_;
    double dt_;
    const PicStepOperation& step_;
    const PicStepCompletion& completion_;
    bool advancing_ = false, failed_ = false;
    int failureCode_ = 0;
    std::optional<PicStepResult> lastStep_;
};

// Binds the module's live state for the duration of the caller's driver.
// Success requires the configured range to finish. Nonzero driver results
// stop the run; original module finalization is performed only after success.
// dt is the effective value after checkpoint restore. A restart already past
// endStep remains an empty run. MacroStep includes the engine's existing
// interleaved IO; neither failure phase promises physical rollback.
int drivePicStepSession(int firstStep, int endStep, double dt, const PicStepOperation& step,
    const PicSessionDriver& driver, const PicStepCompletion& completion = {});
