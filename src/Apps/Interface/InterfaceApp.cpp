#include "Interface/InterfaceApp.h"
#include "Interface/InterfaceTransport.h"
#include <iostream>
#include <string>

int InterfaceApp::run()
{
    const char* stage="input";
    try {
        const auto input=readInterfaceTransportInput();
        if(!config_.checkCase) stage="hermite-gals-first-jet";
        return runInterfaceTransport(input,config_.checkCase);
    } catch(const std::exception& error) {
        std::cerr<<"[LAMP Interface] stage="<<stage<<"; "<<error.what()<<'\n';
        return std::string(stage)=="input"?2:3;
    }
}
