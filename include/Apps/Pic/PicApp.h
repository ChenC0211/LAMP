#pragma once

#include <string>

#include "IApp.h"
#include "RunConfig.h"
#include "Pic/IO/PicInput.h"

class PicApp final : public IApp
{
public:
    explicit PicApp(const RunConfig& config);

    std::string name() const override;
    int run() override;

    static PicAppConfig readConfig(const std::string& runDictPath = "./system/run.dict");

private:
    RunConfig config_;
};