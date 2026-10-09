#pragma once

#include <string>

#include "Runtime/IApp.h"
#include "Runtime/RunConfig.h"

// -----------------------------------------------------------------------------
// PicVerificationApp
//
// App-level driver for verification / regression cases of the Explicit-PIC
// toolchain (cut-cell metrics and FV Poisson regressions).
//
// The app reads ./system/pic.dict and runs the enabled regression branches.
// -----------------------------------------------------------------------------

class PicVerificationApp final : public IApp
{
public:
    explicit PicVerificationApp(const RunConfig& cfg);

    std::string name() const override { return "PicVerification"; }

    int run() override;

private:
    RunConfig cfg_;
};
