#pragma once

#include "Runtime/IApp.h"
#include "Runtime/RunConfig.h"

// Configured planar/RZ dielectric or conducting workflow. The shared Hydro
// case state owns the only mesh/interface; this App owns the main loop. Its adapter supplies
// native face loads and advances charge at the original callback time levels.
// This entry does not bind nozzle boundaries or implement evaporation/emission.
class ElectrohydrodynamicsApp : public IApp
{
public:
    explicit ElectrohydrodynamicsApp(const RunConfig& cfg) : cfg_(cfg) {}
    std::string name() const override { return "Electrohydrodynamics"; }
    int run() override;

private:
    RunConfig cfg_;
};
