#pragma once
#include "Runtime/IModule.h"
#include "Runtime/RunConfig.h"
#include "Modules/Plasma/PlasmaRunInputs.h"
#include "Modules/Plasma/PlasmaCollisionInputs.h"
#include "Modules/Plasma/IO/PicCaseInputLoader.h"
#include <iostream>

// Runtime protocol belongs to the App. Preserve the same immutable service
// owners and PIC selection aliases; do not reread or reinterpret their inputs.
template<class Pic>
class PicModuleAdapter final : public IModule
{
public:
    explicit PicModuleAdapter(const RunConfig& config)
        : units_(config.units), solvers_(config.solvers), pic_(bind(config)) {}
    std::string name() const override { return pic_.name(); }
    int run() override
    {
        auto input = readPicCaseInput(Pic::integrator, units_.get(), solvers_.get());
        if (!input) return -1;
        const int checkpointEvery = input->pic.checkpoint.enabled ? input->pic.checkpoint.every : 0;
        return pic_.withSession(std::move(*input), [checkpointEvery](PicStepSession& session) {
            while (session.nextStep() < session.endStep())
            {
                const int completedStep = session.nextStep() + 1;
                const bool checkpoint = checkpointEvery > 0
                    && (completedStep % checkpointEvery == 0 || completedStep == session.endStep());
                const auto result = session.advance(checkpoint);
                if (!result.accepted())
                {
                    std::cerr << "[PIC] stopped in "
                        << (result.stage == PicStepResult::Stage::MacroStep ? "macro-step" : "completion")
                        << " at step=" << result.time.index + 1 << ", time=" << result.time.end
                        << " (code units), status=" << result.status << '\n';
                    return result.status;
                }
            }
            return 0;
        });
    }
private:
    static PlasmaRunInputs bind(const RunConfig& config)
    {
        PlasmaRunInputs inputs;
        inputs.units = config.units;
        inputs.normalization = config.normalization;
        if (config.plasmaCollisions)
        {
            inputs.gas = config.plasmaCollisions->gas;
            inputs.crossSections = config.plasmaCollisions->crossSections;
            inputs.reactions = config.plasmaCollisions->reactions;
            inputs.collisions = config.plasmaCollisions->collisions;
        }
        inputs.writeSIOutputs = config.writeSIOutputs;
        return inputs;
    }
    std::shared_ptr<const UnitsContext> units_;
    std::shared_ptr<const Solvers::SolversConfig> solvers_;
    Pic pic_;
};
