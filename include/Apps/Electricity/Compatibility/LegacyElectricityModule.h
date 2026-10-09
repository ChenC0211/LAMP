#pragma once

#include <string>

#include "Runtime/IModule.h"
#include "Runtime/RunConfig.h"

// ----------------------------------------------------------------------------
// ElectricityModule
//
// Historical SingleModule/Electricity runner, preserved without numerical
// changes. The production standalone App is ElectricityApp. This adapter is
// deliberately outside the numerical module; it still uses the old model.
// Its SparseLU vertex initializer and old input semantics are not advertised
// as verification of the cell-centered SuperLU production service.
// ----------------------------------------------------------------------------

class ElectricityModule : public IModule
{
public:
    explicit ElectricityModule(const RunConfig& cfg);
    ~ElectricityModule() override = default;

    std::string name() const override { return "Electricity"; }

    int run() override;

private:
    RunConfig cfg_;
};
