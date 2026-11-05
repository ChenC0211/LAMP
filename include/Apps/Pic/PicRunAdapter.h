#pragma once

#include <exception>
#include <iostream>
#include <memory>
#include <utility>

#include "Pic/Setup/PicRunSetup.h"
#include "Pic/IO/PicInput.h"
#include "Core/IO/TecplotWriter.h"
#include "Modules/Plasma/PlasmaRunInputs.h"
#include "Modules/Plasma/PlasmaCollisionInputs.h"

template<class Pic>
class PicRunAdapter final
{
public:
    explicit PicRunAdapter(const PicRunSetup& setup)
        : units_(setup.units), solvers_(setup.solvers), pic_(bindPicRunInputs(setup)) {}

    int run()
    {
        auto input = readPicIntegratorInput(Pic::integrator, units_.get(), solvers_.get());
        if (!input) return -1;
        const int checkpointEvery = input->pic.checkpoint.enabled ? input->pic.checkpoint.every : 0;
        bool ownsSeries=false;
        const auto exportAfterFailure=[&ownsSeries](){
            if(!ownsSeries)return;
            try{TecplotWriter::exportCsvHistories("output/series");}
            catch(const std::exception& e){std::cerr<<"[PIC] series Tecplot export failed: "<<e.what()<<'\n';}
            catch(...){std::cerr<<"[PIC] series Tecplot export failed with an unknown error\n";}
        };
        int status;
        try{
            status=pic_.withSession(std::move(*input), [checkpointEvery,&ownsSeries](PicStepSession& session) {
                ownsSeries=true;
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
        }catch(...){
            // Module-local CSV streams have already unwound and closed.
            exportAfterFailure();throw;
        }
        if(status==0){
            if(ownsSeries)TecplotWriter::exportCsvHistories("output/series");
        }
        else exportAfterFailure();
        return status;
    }

private:
    std::shared_ptr<const UnitsContext> units_;
    std::shared_ptr<const Solvers::SolversConfig> solvers_;
    Pic pic_;
};