#pragma once

#include "Runtime/IModule.h"
#include "Runtime/RunConfig.h"
#include <string>

// Standalone Hydro driver registered through SingleModule for CLI compatibility.
// Owns the main loop; no Runtime dependency remains in the numerical Hydro module.
class HydrodynamicsApp : public IModule
{
public:
    explicit HydrodynamicsApp(const RunConfig& cfg) : cfg_(cfg) {}
    std::string name() const override { return "Hydrodynamics"; }
    int run() override;

private:
    RunConfig cfg_;
};
