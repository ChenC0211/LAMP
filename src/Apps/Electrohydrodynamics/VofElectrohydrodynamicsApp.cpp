#include "Electrohydrodynamics/VofElectrohydrodynamicsApp.h"
#include "Electrohydrodynamics/VofElectrohydrodynamicsRun.h"
#include "Electrohydrodynamics/IO/VofElectrohydrodynamicsOutput.h"

int VofElectrohydrodynamicsApp::run()
{
    const auto input=readVofElectrohydrodynamicsInput("system",config_.threadsOverride,application_);
    VofElectrohydrodynamicsOutput output(input);
    VofElectrohydrodynamicsTransaction transaction(input,output);
    const HydrodynamicsRunIO io{
        [&]{output.begin();},
        [&](auto attempted,auto accepted,double time,double oldDt,double newDt,const std::string& reason){
            output.retry(attempted,accepted,time,oldDt,newDt,reason);},
        [&](auto attempted,auto accepted,double time,const char* stage,const std::string& reason){
            output.failure(attempted,accepted,time,stage,reason);},
        [&]{output.complete();},
        [&](const auto& report,const auto& acceptedClock,int retries){
            output.timeStep(report,acceptedClock,retries);},
        [&](const auto& acceptedClock,double wallSeconds){output.stepCompleted(acceptedClock,wallSeconds);}};
    return runHydrodynamics(input.fluid,config_,transaction,io);
}
