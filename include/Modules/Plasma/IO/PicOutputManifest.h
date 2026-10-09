#pragma once

#include <filesystem>
#include <string>

struct PicConfig;
class UnitsContext;
struct NormalizationConfig;

// Owns output paths, without creating files during construction. The integrator
// invokes each filesystem stage at its original position relative to restore,
// trajectory setup and report opening. No step clock or checkpoint selection.
class PicOutputDirectory
{
public:
    enum class Mode { Explicit, Implicit };
    PicOutputDirectory(Mode mode, const std::filesystem::path& diagnosticsBase,
        std::filesystem::path outputRoot = "./output");
    void create() const;
    // Fresh runs retain the legacy full fields-tree reset; positive thresholds
    // prune only files whose trailing step token is at or above that threshold.
    void prepareFields(int resumeStep) const;
    // Remove only each mode's retired names, non-recursively and without throwing.
    void removeObsoleteReports() const;
    const std::filesystem::path root, fields, series, meta, diagnostics;
private:
    const Mode mode_;
};

// Immutable file-output cadence shared by the run and its manifest. It never
// advances a clock or samples physics. Explicit periodic fields have a separate
// final-output stage; implicit fields remain periodic without a forced end frame.
class PicOutputCadence
{
public:
    using Mode = PicOutputDirectory::Mode;
    PicOutputCadence(Mode mode, const PicConfig& config, int finalStep);
    Mode mode() const noexcept { return mode_; }
    int fieldsEvery() const noexcept { return fieldsEvery_; }
    int seriesEvery() const noexcept { return seriesEvery_; }
    int diagnosticsEvery() const noexcept { return diagnosticsEvery_; }
    int amountEvery() const noexcept { return amountEvery_; }
    bool fields(int step) const noexcept;
    bool series(int step) const noexcept;
    bool periodicSeries(int step) const noexcept;
    bool amounts(int step) const noexcept;
    bool explicitDetailFlush(int step) const noexcept;
    bool explicitFinalFields() const noexcept;
private:
    Mode mode_;
    int finalStep_, fieldsEvery_, seriesEvery_, diagnosticsEvery_, amountEvery_;
    bool fieldsEnabled_;
};

namespace PicOutputManifest
{
// Synchronous, borrowed input. The caller creates output directories and chooses
// the original startup stage; this writer neither changes inputs nor opens steps.
struct Input
{
    const PicConfig& config;
    const std::filesystem::path& outputRoot;
    const std::string& dictionaryPath;
    const UnitsContext* units;
    const NormalizationConfig* normalization;
};
// Retain each integrator's declared catalog, schema order, optional sections and
// warning-only open failures. Declarations are not a list of enabled files.
void write(const Input& input, const PicOutputCadence& cadence);
}
