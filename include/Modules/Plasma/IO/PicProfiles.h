#pragma once

#include "Modules/Plasma/IO/PicFieldSelection.h"
#include "Modules/Plasma/ExplicitPIC/Diagnostics/Profile1DWriter.h"
#include "Modules/Plasma/ExplicitPIC/Diagnostics/Profile1DAvgWriter.h"
#include <filesystem>

// Owns the instantaneous stream and rolling-average window for one PIC run.
// The driver supplies the original sampling time/fields/caches at each call;
// no mesh, physical state, coupling cache or RNG is retained here. Restart
// prepares the step CSV, but starts a fresh average window as before.
class PicProfiles
{
public:
    struct Options
    {
        bool instantaneous = false;
        int bins = 128;
        bool averaged = false;
        int windowSteps = 0;
    };
    struct Sample
    {
        const QuadtreeMesh& mesh;
        const UnitsContext* units;
        const NormalizationContext* norm;
        const CellField<double>& potential;
        const CellField<Vector2D>& electric;
        const std::vector<ParticleType>& species;
        const std::vector<double>& weights;
        const std::vector<double>* alpha;
    };
    using BuildCache = ExplicitPIC::Diagnostics::Profile1DWriter::BuildCacheFn;

    PicProfiles(PicFieldMode mode, const std::filesystem::path& fieldsDirectory,
        int resumeStep, int fieldsEvery, const Options& options);
    bool instantaneousDue(int step) const;
    void writeInstantaneous(int step, double timeCode, double timeSI, const Sample& sample,
        const CellField<double>& charge, const BuildCache& buildCache, bool writeSI);
    // Accumulate every supplied valid step, even when output is not due. The
    // supplied solve caches are borrowed synchronously. Count mismatch skips
    // both accumulation and writing, retaining the previous window.
    void accumulate(int step, const Sample& sample, const std::vector<CouplingWeightsCache>& caches);

private:
    const PicFieldMode mode_;
    const int every_;
    const bool instantaneousEnabled_;
    ExplicitPIC::Diagnostics::Profile1DWriter instantaneous_;
    ExplicitPIC::Diagnostics::Profile1DAvgWriter average_;
};
