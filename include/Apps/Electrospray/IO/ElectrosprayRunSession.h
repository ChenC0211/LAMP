#pragma once

#include "Runtime/RunConfig.h"
#include <filesystem>
#include <iosfwd>
#include <memory>

// Owns the case-local process lock until the App exits. Startup decisions and
// output cleanup happen before any field/history writer or mesh initialization.
// Only output/ may be cleared; system/, cache/ and other case files are retained.
class ElectrosprayRunSession
{
public:
    explicit ElectrosprayRunSession(const std::filesystem::path& caseDirectory);
    ~ElectrosprayRunSession();
    ElectrosprayRunSession(const ElectrosprayRunSession&) = delete;
    ElectrosprayRunSession& operator=(const ElectrosprayRunSession&) = delete;

    // false means user cancellation, with old output unchanged. Resume is not
    // implemented and fails explicitly rather than appending to diagnostic files.
    bool prepare(RunStartMode mode, bool interactive, std::istream& input, std::ostream& output);
    static bool terminalInput();

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};
