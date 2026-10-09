#pragma once

#include <string>

#include "Runtime/IApp.h"
#include "Runtime/RunConfig.h"

// Application-level driver for electrospray workflows.
class ElectrosprayApp : public IApp
{
public:
    explicit ElectrosprayApp(const RunConfig& cfg);
    ~ElectrosprayApp() override = default;

    std::string name() const override { return "Electrospray"; }

    int run() override;

private:
    RunConfig cfg_;
};
