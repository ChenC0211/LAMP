#pragma once

#include <map>
#include <string>

#include "Verification/Pic/DictionaryPicVerification.h"

// -----------------------------------------------------------------------------
// RegressionHarness
//
// Minimal baseline record/check utility for App-level regressions.
//
// The harness operates on a set of scalar metrics (name -> value) produced by
// a verification App. Depending on mode, it either:
//   - records (writes) a baseline file
//   - checks (compares) against an existing baseline file
//   - is disabled
//
// Baseline files are case-local (recommended under ./cache/).
// -----------------------------------------------------------------------------

class RegressionHarness
{
public:
    using Metrics = std::map<std::string, double>;

    explicit RegressionHarness(const PicVerificationConfig::RegressionHarnessConfig& cfg);

    // Execute record/check logic.
    // Return 0 on success.
    // Return non-zero on missing baseline or metric mismatch.
    int process(const Metrics& metrics) const;

private:
    enum class Mode
    {
        Off,
        Record,
        Check
    };

    Mode mode_ = Mode::Off;
    std::string baselinePath_;
    std::string reportPath_;
    double absTol_ = 1e-12;
    double relTol_ = 1e-10;

    static Mode parseMode(const std::string& s);

    static bool readBaseline(const std::string& path, Metrics& out);
    static bool writeBaseline(const std::string& path, const Metrics& metrics);
    bool writeReport(const std::string& path,
                     const Metrics& baseline,
                     const Metrics& current,
                     const std::string& summary) const;

    bool nearlyEqual(double a, double b) const;
};
