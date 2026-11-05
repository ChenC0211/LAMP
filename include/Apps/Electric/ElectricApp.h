#pragma once
#include "IApp.h"
#include "RunConfig.h"

// Production standalone electrostatics/charging on a fixed supplied interface.
// Uses the same Electric service as EHD, without constructing Hydro.
class ElectricApp final : public IApp {
public:
    explicit ElectricApp(const RunConfig& config) : config_(config) {}
    std::string name() const override { return "Electric"; }
    int run() override;
private:
    RunConfig config_;
};
