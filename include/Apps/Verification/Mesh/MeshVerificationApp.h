#pragma once
#include "Runtime/IApp.h"
#include "Runtime/RunConfig.h"

// Independent Core mesh workflow. Owns no physical model or simulation clock.
class MeshVerificationApp : public IApp
{
public:
    explicit MeshVerificationApp(const RunConfig& config) : checkOnly_(config.checkCase) {}
    std::string name() const override { return "MeshVerification"; }
    int run() override;
private:
    bool checkOnly_;
};
