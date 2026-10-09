#pragma once
#include "Runtime/IApp.h"
#include "Runtime/RunConfig.h"

// Production standalone electrostatics/charging on a fixed supplied interface.
// Uses the same Electricity service as EHD, without constructing Hydro.
class ElectricityApp final : public IApp {
public:
    explicit ElectricityApp(const RunConfig& config) : config_(config) {}
    std::string name() const override { return "Electricity"; }
    int run() override;
private:
    RunConfig config_;
};
