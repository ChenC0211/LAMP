#pragma once
#include "IApp.h"
#include "RunConfig.h"
// Formal native HYDRO / FirstJet App; clock, transactions and IO only.
class HydrodynamicsApp final:public IApp {
public:
    explicit HydrodynamicsApp(const RunConfig& c):config_(c){}
    std::string name()const override{return "Hydrodynamics";}
    int run()override;
private:
    RunConfig config_;
};
