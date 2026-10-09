#pragma once

#include <string>
#include <vector>

#include "Runtime/IApp.h"
#include "Runtime/RunConfig.h"

// ----------------------------------------------------------------------------
// SingleModuleApp
//
// Runs one or more Modules selected in ./system/run.dict.
//
// Expected syntax:
//
//   Run { app SingleModule }
//
//   SingleModule
//   {
//       modules ( Electricity ExplicitPIC )
//   }
//
// 发布级规范：system/ 仅允许 mesh.dict/geometry.dict/run.dict 及各模块 dict。
// ----------------------------------------------------------------------------

class SingleModuleApp final : public IApp
{
public:
    explicit SingleModuleApp(const RunConfig& cfg);

    std::string name() const override;
    int run() override;

    // Shared by Runtime preflight and execution. Returns no modules for a missing,
    // empty, duplicate, or malformed block/list; never returns a partial list.
    static std::vector<std::string> readModulesFromRunDict(const std::string& path);

private:
    RunConfig cfg_;
};
