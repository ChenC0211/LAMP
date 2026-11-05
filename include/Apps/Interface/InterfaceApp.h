#pragma once

#include "IApp.h"
#include "RunConfig.h"

// Independent Interface runner. The same numerical services are used by
// Electrospray; prescribed data replace the external Hydro input provider.
class InterfaceApp final : public IApp
{
public:
    explicit InterfaceApp(const RunConfig& config) : config_(config) {}
    std::string name() const override { return "Interface"; }
    int run() override;
private:
    RunConfig config_;
};
