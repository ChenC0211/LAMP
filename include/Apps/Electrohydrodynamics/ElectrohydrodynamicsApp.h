#pragma once
#include "IApp.h"
#include "RunConfig.h"
// Formal native HYDRO + sharp-interface electrostatics App; owns the joint accepted transaction.
class ElectrohydrodynamicsApp final:public IApp {
public:
    explicit ElectrohydrodynamicsApp(const RunConfig& c):config_(c){}
    std::string name()const override{return "Electrohydrodynamics";}
    int run()override;
private:
    RunConfig config_;
};
