#include "Modules/Plasma/IO/PicProfiles.h"
#include "Modules/Plasma/ExplicitPIC/Deposit/CouplingWeightsCache.h"
#include <algorithm>
#include <iostream>

PicProfiles::PicProfiles(PicFieldMode mode, const std::filesystem::path& fieldsDirectory,
    int resumeStep, int fieldsEvery, const Options& options)
    : mode_(mode), every_(std::max(1, fieldsEvery)), instantaneousEnabled_(options.instantaneous)
{
    const auto instantPath = fieldsDirectory / "profile1d.csv";
    ExplicitPIC::Diagnostics::Profile1DWriterConfig instant;
    instant.enabled = options.instantaneous;
    instant.every = every_;
    instant.nbinsX = mode == PicFieldMode::Implicit ? std::max(1, options.bins) : options.bins;
    instant.path = instantPath.string();
    instantaneous_.setConfig(instant);
    if (options.instantaneous)
    {
        instantaneous_.open(instantPath.string(), resumeStep);
        if (!instantaneous_.isOpen())
            std::cerr << '[' << (mode == PicFieldMode::Explicit ? "ExplicitPIC" : "ImplicitPIC")
                      << "] WARNING: failed to create profile1d csv: " << instantPath.string() << std::endl;
    }

    ExplicitPIC::Diagnostics::Profile1DAvgWriterConfig average;
    average.enabled = options.averaged && options.windowSteps > 0;
    average.every = every_;
    average.avgWindowSteps = std::max(0, options.windowSteps);
    average.path = (fieldsDirectory / "profiles_avg.csv").string();
    average_.setConfig(average);
    const std::size_t bins = options.bins > 0 ? static_cast<std::size_t>(options.bins) : static_cast<std::size_t>(128);
    if (average_.enabled())
        average_.reset(bins);
}

bool PicProfiles::instantaneousDue(int step) const
{
    return instantaneousEnabled_ && instantaneous_.isOpen() &&
        ((mode_ == PicFieldMode::Explicit && step == 0) || step % every_ == 0);
}

void PicProfiles::writeInstantaneous(int step, double timeCode, double timeSI, const Sample& sample,
    const CellField<double>& charge, const BuildCache& buildCache, bool writeSI)
{
    if (!instantaneousDue(step))
        return;
    instantaneous_.write(step, timeCode, timeSI, sample.mesh, sample.units, sample.norm,
        sample.potential, sample.electric, charge, sample.species, sample.weights, buildCache, sample.alpha, writeSI);
}

void PicProfiles::accumulate(int step, const Sample& sample, const std::vector<CouplingWeightsCache>& caches)
{
    if (!average_.enabled() || caches.size() != sample.species.size())
        return;
    const auto cacheRef = [&](std::size_t si) -> const CouplingWeightsCache& { return caches[si]; };
    average_.accumulate(sample.mesh, sample.units, sample.norm, sample.potential, sample.electric,
        sample.species, sample.weights, cacheRef, sample.alpha);
    if (step % every_ == 0 && average_.windowReady())
        average_.writeSnapshotCSV(sample.mesh, sample.units);
}
