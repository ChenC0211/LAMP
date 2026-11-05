#include "Electrohydrodynamics/ElectrohydrodynamicsApp.h"
#include "Electrohydrodynamics/ElectrohydrodynamicsRun.h"

int ElectrohydrodynamicsApp::run(){
    const auto input=readElectrohydrodynamicsInput("system",config_.threadsOverride);
    return runElectrohydrodynamics(input,config_);
}
