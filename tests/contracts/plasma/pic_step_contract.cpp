#include "Modules/Plasma/PicStepSession.h"
#include <iostream>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace {
void check(bool good, const char* why) { if (!good) throw std::runtime_error(why); }
template<class Action> void rejects(Action action, const char* why)
{
    bool rejected = false;
    try { action(); } catch (const std::logic_error&) { rejected = true; }
    check(rejected, why);
}
}

int main() try
{
    static_assert(!std::is_copy_constructible_v<PicStepSession>);
    static_assert(!std::is_move_constructible_v<PicStepSession>);
    std::vector<int> executed;
    const auto action = [&](const PicStepTime& time) { executed.push_back(time.index); return 0; };
    check(drivePicStepSession(2, 5, .125, action, [&](PicStepSession& session) {
        check(!session.lastStep() && session.time() == .25 && session.dt() == .125,
            "restored session did not start at its effective time");
        for (int expected = 2; expected < 5; ++expected)
        {
            check(session.nextStep() == expected && session.endStep() == 5, "restart range changed");
            const auto result = session.advance();
            check(result.accepted() && !result.threw && result.status == 0
                && result.time.index == expected && result.time.begin == expected * .125
                && result.time.end == (expected + 1) * .125 && result.time.dt == .125
                && !result.checkpointRequested && session.time() == result.time.end
                && session.lastStep()->accepted() && session.nextStep() == expected + 1,
                "single-step completion and published interval disagree");
        }
        rejects([&] { session.advance(); }, "advanced beyond configured range");
        return 0;
    }) == 0 && executed == std::vector<int>({2, 3, 4}), "exactly one call per requested step");

    rejects([&] { drivePicStepSession(0, 2, .125, action, [](PicStepSession& session) {
        session.advance(); return 0;
    }); }, "early return claimed a complete run");
    const auto before = executed.size();
    check(drivePicStepSession(0, 2, .125, action, [](PicStepSession&) { return 7; }) == 7
        && executed.size() == before, "driver cancellation advanced state");

    int failedCalls = 0;
    check(drivePicStepSession(0, 2, .125, [&](const PicStepTime&) { ++failedCalls; return -4; }, [](PicStepSession& session) {
        check(session.advance().status == -4 && session.nextStep() == 0, "failed step was accepted");
        check(session.stopped() && session.time() == 0 && !session.lastStep()->accepted()
            && session.lastStep()->stage == PicStepResult::Stage::MacroStep
            && !session.lastStep()->threw, "physical return failure lost its stage");
        rejects([&] { session.advance(); }, "failed step was retried");
        return 0; // A caller cannot erase the original physical failure code.
    }) == -4 && failedCalls == 1, "physical failure was hidden");

    rejects([&] { drivePicStepSession(0, 1, .125, [](const PicStepTime&) -> int { throw std::runtime_error("partial step"); },
        [](PicStepSession& session) {
            try { session.advance(); } catch (const std::runtime_error&) {}
            check(session.lastStep()->threw && session.stopped()
                && session.lastStep()->stage == PicStepResult::Stage::MacroStep,
                "physical exception lost its stage");
            rejects([&] { session.advance(); }, "exception allowed another step");
            return 0;
        }); }, "suppressed physical exception was reported as success");

    PicStepSession* active = nullptr;
    check(drivePicStepSession(0, 1, .125, [&](const PicStepTime&) {
        rejects([&] { active->advance(); }, "recursive step entered the physical kernel");
        return 0;
    }, [&](PicStepSession& session) { active = &session; return session.advance().status; }) == 0,
        "reentrancy guard altered the outer step");
    check(drivePicStepSession(5, 3, .125, action, [](PicStepSession&) { return 0; }) == 0
        && executed.size() == before, "already-past-end restart executed a step");

    std::vector<int> phases;
    check(drivePicStepSession(1, 4, .125, [&](const PicStepTime& time) { phases.push_back(time.index); return 0; },
        [](PicStepSession& session) {
            check(session.advance(true).status == 0, "first requested checkpoint");
            check(session.advance(false).status == 0, "step without checkpoint");
            return session.advance(true).status;
        }, [&](const PicStepTime& time, bool checkpoint) {
            phases.push_back(checkpoint ? 10 + time.index : 20 + time.index); return 0;
        }) == 0 && phases == std::vector<int>({1, 11, 2, 22, 3, 13}),
        "checkpoint decision or physical/completion order changed");
    int completions = 0;
    check(drivePicStepSession(0, 1, .125, [](const PicStepTime&) { return -2; }, [](PicStepSession& session) {
        return session.advance(true).status;
    }, [&](const PicStepTime&, bool) { ++completions; return 0; }) == -2 && completions == 0,
        "checkpoint ran after a failed physical step");
    check(drivePicStepSession(0, 2, .125, [](const PicStepTime&) { return 0; }, [](PicStepSession& session) {
        check(session.advance(true).status == -3 && session.nextStep() == 0,
            "failed completion advanced the accepted step");
        check(session.lastStep()->stage == PicStepResult::Stage::Completion
            && session.lastStep()->checkpointRequested && !session.lastStep()->accepted()
            && session.time() == 0 && session.stopped(), "completion failure lost its phase or accepted time");
        rejects([&] { session.advance(); }, "failed checkpoint completion was retried");
        return 0;
    }, [](const PicStepTime&, bool) { return -3; }) == -3, "checkpoint failure was hidden");
    rejects([&] { drivePicStepSession(0, 1, .125, [](const PicStepTime&) { return 0; }, [](PicStepSession& session) {
        try { session.advance(true); } catch (const std::runtime_error&) {}
        check(session.lastStep()->stage == PicStepResult::Stage::Completion && session.lastStep()->threw,
            "completion exception lost its stage");
        return 0;
    }, [](const PicStepTime&, bool) -> int { throw std::runtime_error("checkpoint failure"); }); },
        "completion exception was reported as success");
    int unsavedSteps = 0;
    check(drivePicStepSession(0, 1, .125, [&](const PicStepTime&) { ++unsavedSteps; return 0; }, [](PicStepSession& session) {
        rejects([&] { session.advance(true); }, "missing checkpoint operation was ignored");
        check(session.nextStep() == 0 && !session.lastStep() && !session.stopped(),
            "invalid checkpoint request advanced state");
        return session.advance().status;
    }) == 0 && unsavedSteps == 1, "checkpoint preflight altered physical execution");
    const PicSessionDriver noSteps = [](PicStepSession&) { return 0; };
    for (double dt : {0., -1., std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::max()})
        rejects([&] { drivePicStepSession(0, 2, dt, action, noSteps); }, "invalid time range reached driver");
    rejects([&] { drivePicStepSession(-1, 2, .125, action, noSteps); }, "negative first step reached driver");
    rejects([&] { drivePicStepSession(0, -1, .125, action, noSteps); }, "negative end step reached driver");
    PicStepTime observed{};
    const int final = std::numeric_limits<int>::max();
    check(drivePicStepSession(final - 1, final, .1,
        [&](const PicStepTime& time) { observed = time; return 0; },
        [&](PicStepSession& session) {
            const auto saved = session.advance(true);
            rejects([&] { session.advance(); }, "last representable step overflowed");
            check(saved.accepted() && saved.time.end == static_cast<double>(final) * .1,
                "final time accumulated instead of multiplying the step index");
            return saved.status;
        }, [&](const PicStepTime& time, bool checkpoint) {
            check(time.index == observed.index && time.begin == observed.begin && time.end == observed.end
                && time.dt == observed.dt && checkpoint, "completion used another macro-step interval");
            return 0;
        }) == 0, "last representable step did not finish");
    std::cout << "PIC caller-controlled intervals/results, restart range and failure isolation passed\n";
    return 0;
}
catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
