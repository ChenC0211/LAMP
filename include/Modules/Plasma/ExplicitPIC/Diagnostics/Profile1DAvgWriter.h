#pragma once

#include <cstddef>
#include <functional>
#include <string>
#include <vector>

#include "Core/Base/Vector2D.h"
#include "Core/Field/CellField.h"
#include "Modules/Plasma/Particles/ParticleType.h"

#include "Modules/Plasma/ExplicitPIC/Diagnostics/Profile1DRollingAccumulator.h"

class QuadtreeMesh;
class UnitsContext;
class NormalizationContext;
class CouplingWeightsCache;

namespace ExplicitPIC
{
namespace Diagnostics
{

// Profile1DAvgWriter
// ---------------------------------------------------------------------------
// Computes and writes time-averaged 1D (x-binned, y-averaged) profiles for
// Turner-style CCP benchmarks.
//
// Output: a single CSV snapshot (overwritten on each write) with columns
// designed to be consumed directly by postprocess/compare_turner.py:
//   bin,x_m,
//   phi_mean_V,phi_std_V,
//   E_mean_Vpm,E_std_Vpm,
//   ne_mean_m-3,ne_std_m-3,
//   ni_mean_m-3,ni_std_m-3,
//   Te_mean_eV,Te_std_eV
//
// The writer maintains a fixed-size sliding window in step units through a
// rolling accumulator. For minimal steady-periodic diagnostics, it is typical
// to accumulate only the last N steps.
// ---------------------------------------------------------------------------

struct Profile1DAvgWriterConfig
{
    bool enabled = false;
    int every = 50;             // write cadence (steps). If <=0, caller decides.
    int avgWindowSteps = 0;     // window size in steps; 0 disables
    std::string path = "./output/profiles_avg.csv";
};

class Profile1DAvgWriter
{
public:
    using CacheRefFn = std::function<const CouplingWeightsCache&(std::size_t speciesIndex)>;

    enum Field
    {
        NE = 0,
        NI = 1,
        TE = 2,
        PHI = 3,
        E = 4,
        NUM_FIELDS = 5
    };

    Profile1DAvgWriter() = default;

    void setConfig(const Profile1DAvgWriterConfig& cfg) { cfg_ = cfg; }
    const Profile1DAvgWriterConfig& config() const { return cfg_; }

    bool enabled() const { return cfg_.enabled && cfg_.avgWindowSteps > 0; }

    // Configure bin count and window size; resets accumulator.
    void reset(std::size_t nBinsX);

    // Push one step sample into the rolling window.
    void accumulate(const QuadtreeMesh& mesh,
                    const UnitsContext* units,
                    const NormalizationContext* norm,
                    const CellField<double>& phi,
                    const CellField<Vector2D>& Ecell,
                    const std::vector<ParticleType>& species,
                    const std::vector<double>& macroWeight,
                    const CacheRefFn& cacheForSpecies,
                    const std::vector<double>* alphaByCellIndex);

    // Write current mean/std snapshot to CSV (overwrites file each time).
    void writeSnapshotCSV(const QuadtreeMesh& mesh, const UnitsContext* units) const;

    bool windowReady() const { return enabled() && acc_.count() >= acc_.windowSteps(); }
    std::size_t count() const { return acc_.count(); }

private:
    Profile1DAvgWriterConfig cfg_;
    std::size_t nBinsX_ = 0;
    Profile1DRollingAccumulator acc_;

    static void ensureParentDir(const std::string& filePath);
};

} // namespace Diagnostics
} // namespace ExplicitPIC
