#pragma once

#include <string>

#include "Core/IO/Dictionary.h"
#include "Runtime/RunConfig.h"

// ----------------------------------------------------------------------------
// DictionaryRun
//
// 发布级规范：run.dict 的 Run{...} 块只负责选择“要运行的 App”。
// 其它与 App 相关的运行参数应放在与 App 同名的块中（例如 SingleModule{...}、Electrospray{...}），
// 并由对应的 App 自己解析；DictionaryRun 只解析 Run{app ...}。
//
// 允许语法（建议使用花括号块）：
//
// Run
// {
//     app PicVerification
// }
// ----------------------------------------------------------------------------

class DictionaryRun : public Dictionary
{
public:
    DictionaryRun() = default;
    ~DictionaryRun() = default;

    // Return false if file does not exist / cannot be opened.
    bool readDict(const std::string& runDictPath = "./system/run.dict");

    const RunConfig& config() const { return config_; }

private:
    RunConfig config_{};
    static std::string toLowerCopy(std::string s);
};
