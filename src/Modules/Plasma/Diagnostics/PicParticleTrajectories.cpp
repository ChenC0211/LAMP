#include "Modules/Plasma/Diagnostics/PicParticleTrajectories.h"
#include "Modules/Plasma/Diagnostics/PhysicsDiagnostics.h"
#include "Modules/Plasma/Diagnostics/TimeSeriesIO.h"
#include "Modules/Plasma/Particles/ParticleType.h"
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace
{
struct ParticleTrajectoryTrack
{
    std::uint64_t trackId = 0u;
    int trackSlot = 0;
    int speciesCode = 0;
    int initialIndex = -1;
    int samplesWritten = 0;
    int maxSamples = 0;
    bool finished = false;
    std::size_t speciesIndex = 0u;
    std::string speciesName;
    struct Sample
    {
        int step = 0;
        double time_s = 0.0;
        int alive = 0;
        double x_m = std::numeric_limits<double>::quiet_NaN();
        double y_m = std::numeric_limits<double>::quiet_NaN();
        double vx_mps = std::numeric_limits<double>::quiet_NaN();
        double vy_mps = std::numeric_limits<double>::quiet_NaN();
        double vz_mps = std::numeric_limits<double>::quiet_NaN();
        int inCellL = -1;
        double weight = std::numeric_limits<double>::quiet_NaN();
    };
    std::vector<Sample> samples;
};

struct ParticleTrajectoryGroup
{
    std::string prefix;
    int targetCount = 0;
    int maxSamples = 0;
    int sampleEvery = 1;
    std::vector<ParticleTrajectoryTrack> tracks;
};

static void clearParticleTrajectoryRoot(const std::filesystem::path& root)
{
    std::error_code ec;
    std::filesystem::create_directories(root, ec);
    for (const auto& entry : std::filesystem::directory_iterator(root, ec))
    {
        if (ec)
            break;
        std::filesystem::remove_all(entry.path(), ec);
        if (ec)
            break;
    }
}

static bool speciesMatchesTrajectoryRole(const ParticleType& species, bool wantElectrons)
{
    return wantElectrons ? (species.charge() < 0.0) : (species.charge() > 0.0);
}

static void assignMissingParticleTrackIds(std::vector<ParticleType>& species, std::uint64_t& nextTrackId)
{
    if (nextTrackId == 0u)
        nextTrackId = 1u;
    for (auto& sp : species)
    {
        for (auto& particle : sp.particle)
        {
            if (particle.trackId == 0u)
                particle.trackId = nextTrackId++;
        }
        sp.number = static_cast<int>(sp.particle.size());
    }
}

static void initParticleTrajectoryGroup(ParticleTrajectoryGroup& group,
                                        const std::vector<ParticleType>& species,
                                        bool wantElectrons,
                                        const std::filesystem::path& root,
                                        const std::string& moduleTag,
                                        std::mt19937_64& rng)
{
    group.tracks.clear();
    if (group.targetCount <= 0 || group.maxSamples <= 0)
        return;

    struct Candidate
    {
        std::size_t speciesIndex = 0u;
        int particleIndex = -1;
        std::uint64_t trackId = 0u;
        std::string speciesName;
    };

    std::vector<Candidate> candidates;
    for (std::size_t si = 0; si < species.size(); ++si)
    {
        const auto& sp = species[si];
        if (!speciesMatchesTrajectoryRole(sp, wantElectrons))
            continue;
        for (std::size_t pi = 0; pi < sp.particle.size(); ++pi)
        {
            const auto& particle = sp.particle[pi];
            if (particle.trackId == 0u)
                continue;
            Candidate cand;
            cand.speciesIndex = si;
            cand.particleIndex = static_cast<int>(pi);
            cand.trackId = particle.trackId;
            cand.speciesName = sp.name();
            candidates.push_back(std::move(cand));
        }
    }

    if (candidates.empty())
        return;

    std::filesystem::create_directories(root);
    const int wanted = std::min(group.targetCount, static_cast<int>(candidates.size()));
    std::shuffle(candidates.begin(), candidates.end(), rng);

    group.tracks.reserve(static_cast<std::size_t>(wanted));
    for (int k = 0; k < wanted; ++k)
    {
        const auto& cand = candidates[static_cast<std::size_t>(k)];
        ParticleTrajectoryTrack track;
        track.trackId = cand.trackId;
        track.trackSlot = k + 1;
        track.speciesCode = wantElectrons ? 0 : 1;
        track.initialIndex = cand.particleIndex;
        track.maxSamples = group.maxSamples;
        track.speciesIndex = cand.speciesIndex;
        track.speciesName = cand.speciesName;
        track.samples.reserve(static_cast<std::size_t>(std::max(0, group.maxSamples)));
        group.tracks.push_back(std::move(track));
    }
}

static void dumpParticleTrajectoryGroup(ParticleTrajectoryGroup& group,
                                        const std::vector<ParticleType>& species,
                                        int step,
                                        double time_s,
                                        double lengthScale_m,
                                        double velocityScale_m)
{
    const int sampleEvery = std::max(1, group.sampleEvery);
    if ((step % sampleEvery) != 0)
        return;

    for (std::size_t trackSlot = 0; trackSlot < group.tracks.size(); ++trackSlot)
    {
        auto& track = group.tracks[trackSlot];
        if (track.finished)
            continue;
        if (track.samplesWritten >= track.maxSamples)
        {
            track.finished = true;
            continue;
        }

        const Particle* found = nullptr;
        if (track.speciesIndex < species.size())
        {
            for (const auto& particle : species[track.speciesIndex].particle)
            {
                if (particle.trackId == track.trackId)
                {
                    found = &particle;
                    break;
                }
            }
        }

        if (found != nullptr)
        {
            ParticleTrajectoryTrack::Sample sample;
            sample.step = step;
            sample.time_s = time_s;
            sample.alive = 1;
            sample.x_m = found->position.x * lengthScale_m;
            sample.y_m = found->position.y * lengthScale_m;
            sample.vx_mps = found->velocity.x * velocityScale_m;
            sample.vy_mps = found->velocity.y * velocityScale_m;
            sample.vz_mps = found->vz * velocityScale_m;
            sample.inCellL = found->inCellL;
            sample.weight = found->w;
            track.samples.push_back(sample);
            ++track.samplesWritten;
            if (track.samplesWritten >= track.maxSamples)
                track.finished = true;
        }
        else
        {
            track.finished = true;
        }
    }
}

static void writeParticleTrajectoryTecplot(const std::filesystem::path& outPath,
                                           const ParticleTrajectoryGroup& electronGroup,
                                           const ParticleTrajectoryGroup& ionGroup,
                                           const std::string& moduleTag)
{
    std::filesystem::create_directories(outPath.parent_path());

    std::ofstream out(outPath.string(), std::ios::out | std::ios::trunc);
    if (!out)
    {
        std::cerr << '[' << moduleTag << "][trajectory] WARNING: failed to open "
                  << outPath.string() << std::endl;
        return;
    }

    out << "Title = \"Particle Trajectories\"\n";
    out << "Variables = \"step\" \"time_s\" \"track_slot\" \"track_id\" \"species_code\" \"alive\" \"x_m\" \"y_m\" \"vx_mps\" \"vy_mps\" \"vz_mps\" \"inCellL\" \"weight\" \"initial_index\"\n";
    out << std::setprecision(16) << std::fixed;

    auto writeGroup = [&](const ParticleTrajectoryGroup& group) {
        for (const auto& track : group.tracks)
        {
            if (track.samples.empty())
                continue;

            out << "ZONE T=\""
                << group.prefix << "_track_" << track.trackSlot << '_' << TimeSeriesIO::sanitizeCsvToken(track.speciesName)
                << "_id" << track.trackId << "\", I=" << track.samples.size() << ", F=POINT\n";
            for (const auto& sample : track.samples)
            {
                out << sample.step << ' '
                    << sample.time_s << ' '
                    << track.trackSlot << ' '
                    << track.trackId << ' '
                    << track.speciesCode << ' '
                    << sample.alive << ' '
                    << sample.x_m << ' '
                    << sample.y_m << ' '
                    << sample.vx_mps << ' '
                    << sample.vy_mps << ' '
                    << sample.vz_mps << ' '
                    << sample.inCellL << ' '
                    << sample.weight << ' '
                    << track.initialIndex << '\n';
            }
        }
    };

    writeGroup(electronGroup);
    writeGroup(ionGroup);
    out.flush();
}

static std::filesystem::path makeParticleTrajectorySnapshotPath(const std::filesystem::path& root,
                                                               const std::string& prefix,
                                                               int step)
{
    std::ostringstream name;
    name << "particle_trajectories_" << TimeSeriesIO::sanitizeCsvToken(prefix) << "_"
         << std::setw(7) << std::setfill('0') << std::max(0, step) << ".dat";
    return root / name.str();
}

static void writeSingleParticleTrajectoryTecplot(const std::filesystem::path& outPath,
                                                 const ParticleTrajectoryGroup& group,
                                                 bool isElectronGroup,
                                                 const char* moduleTag)
{
    const ParticleTrajectoryGroup emptyGroup{};
    if (isElectronGroup)
    {
        writeParticleTrajectoryTecplot(outPath, group, emptyGroup, moduleTag);
    }
    else
    {
        writeParticleTrajectoryTecplot(outPath, emptyGroup, group, moduleTag);
    }
}

static bool hasActiveParticleTrajectoryTracks(const ParticleTrajectoryGroup& group)
{
    for (const auto& track : group.tracks)
    {
        if (!track.finished && track.samplesWritten < track.maxSamples)
            return true;
    }
    return false;
}
}

struct PicParticleTrajectories::Impl
{
    bool configured = false;
    bool particleTrajectoryEnabled = false;
    int electronTrajectoryEvery = 0;
    int ionTrajectoryEvery = 0;
    std::filesystem::path particleTrajectoryRoot;
    ParticleTrajectoryGroup electronTrajectoryGroup;
    ParticleTrajectoryGroup ionTrajectoryGroup;
    std::filesystem::path activeElectronTrajectoryPath;
    std::filesystem::path activeIonTrajectoryPath;
    bool electronTrajectoryBatchActive = false;
    bool ionTrajectoryBatchActive = false;
    std::uint64_t nextParticleTrackId = 1u;
    std::string moduleTag;

    void requireConfigured() const
    {
        if (!configured) throw std::logic_error("PIC trajectories require configuration");
    }

    void configure(const PicParticleTrajectoryConfig& config, int seriesEvery,
                   const std::filesystem::path& seriesDir, const std::string& tag)
    {
        if (configured) throw std::logic_error("PIC trajectories cannot be reconfigured");
        const bool electronTrajectoryEnabled = config.enabled &&
            config.electronCount > 0 &&
            config.electronMaxSamples > 0;
        const bool ionTrajectoryEnabled = config.enabled &&
            config.ionCount > 0 &&
            config.ionMaxSamples > 0;
        particleTrajectoryEnabled = electronTrajectoryEnabled || ionTrajectoryEnabled;
        const int particleTrajectoryEvery = PhysicsDiagnostics::resolveEvery(seriesEvery, config.every);
        electronTrajectoryEvery = electronTrajectoryEnabled
            ? PhysicsDiagnostics::resolveEvery(particleTrajectoryEvery, config.electronEvery)
            : 0;
        ionTrajectoryEvery = ionTrajectoryEnabled
            ? PhysicsDiagnostics::resolveEvery(particleTrajectoryEvery, config.ionEvery)
            : 0;
        particleTrajectoryRoot = seriesDir / config.baseDir;
        electronTrajectoryGroup = ParticleTrajectoryGroup{"electron",
                                                       config.electronCount,
                                                       config.electronMaxSamples,
                                                       std::max(1, config.electronSampleEvery)};
        ionTrajectoryGroup = ParticleTrajectoryGroup{"ion",
                                                   config.ionCount,
                                                   config.ionMaxSamples,
                                                   std::max(1, config.ionSampleEvery)};
        activeElectronTrajectoryPath.clear();
        activeIonTrajectoryPath.clear();
        electronTrajectoryBatchActive = false;
        ionTrajectoryBatchActive = false;
        moduleTag = tag;
        configured = true;
    }

    void sample(std::vector<ParticleType>& species, std::mt19937_64& rng,
                int trajectoryStep, double time_s, double lengthScale_m, double velocityScale_m)
    {
        requireConfigured();
        if (particleTrajectoryEnabled)
        {
            const bool startElectronTrajectoryBatch = electronTrajectoryEvery > 0 && ((trajectoryStep % electronTrajectoryEvery) == 0);
            const bool startIonTrajectoryBatch = ionTrajectoryEvery > 0 && ((trajectoryStep % ionTrajectoryEvery) == 0);
            if (startElectronTrajectoryBatch || startIonTrajectoryBatch)
            {
                assignMissingParticleTrackIds(species, nextParticleTrackId);
            }

            if (startElectronTrajectoryBatch)
            {
                if (electronTrajectoryBatchActive)
                {
                    writeSingleParticleTrajectoryTecplot(activeElectronTrajectoryPath,
                                                         electronTrajectoryGroup,
                                                         true,
                                                         moduleTag.c_str());
                }
                initParticleTrajectoryGroup(electronTrajectoryGroup, species, true, particleTrajectoryRoot, moduleTag, rng);
                activeElectronTrajectoryPath = makeParticleTrajectorySnapshotPath(particleTrajectoryRoot,
                                                                                  electronTrajectoryGroup.prefix,
                                                                                  trajectoryStep);
                electronTrajectoryBatchActive = true;
            }

            if (startIonTrajectoryBatch)
            {
                if (ionTrajectoryBatchActive)
                {
                    writeSingleParticleTrajectoryTecplot(activeIonTrajectoryPath,
                                                         ionTrajectoryGroup,
                                                         false,
                                                         moduleTag.c_str());
                }
                initParticleTrajectoryGroup(ionTrajectoryGroup, species, false, particleTrajectoryRoot, moduleTag, rng);
                activeIonTrajectoryPath = makeParticleTrajectorySnapshotPath(particleTrajectoryRoot,
                                                                             ionTrajectoryGroup.prefix,
                                                                             trajectoryStep);
                ionTrajectoryBatchActive = true;
            }

            if (electronTrajectoryBatchActive)
            {
                dumpParticleTrajectoryGroup(electronTrajectoryGroup, species, trajectoryStep, time_s, lengthScale_m, velocityScale_m);
                if (!hasActiveParticleTrajectoryTracks(electronTrajectoryGroup))
                {
                    writeSingleParticleTrajectoryTecplot(activeElectronTrajectoryPath,
                                                         electronTrajectoryGroup,
                                                         true,
                                                         moduleTag.c_str());
                    electronTrajectoryBatchActive = false;
                }
            }

            if (ionTrajectoryBatchActive)
            {
                dumpParticleTrajectoryGroup(ionTrajectoryGroup, species, trajectoryStep, time_s, lengthScale_m, velocityScale_m);
                if (!hasActiveParticleTrajectoryTracks(ionTrajectoryGroup))
                {
                    writeSingleParticleTrajectoryTecplot(activeIonTrajectoryPath,
                                                         ionTrajectoryGroup,
                                                         false,
                                                         moduleTag.c_str());
                    ionTrajectoryBatchActive = false;
                }
            }
        }
    }

    void finish()
    {
        requireConfigured();
        if (electronTrajectoryBatchActive)
        {
            writeSingleParticleTrajectoryTecplot(activeElectronTrajectoryPath,
                                                 electronTrajectoryGroup,
                                                 true,
                                                 moduleTag.c_str());
        }
        if (ionTrajectoryBatchActive)
        {
            writeSingleParticleTrajectoryTecplot(activeIonTrajectoryPath,
                                                 ionTrajectoryGroup,
                                                 false,
                                                 moduleTag.c_str());
        }
    }
};

PicParticleTrajectories::PicParticleTrajectories() : impl_(std::make_unique<Impl>()) {}
PicParticleTrajectories::~PicParticleTrajectories() = default;

void PicParticleTrajectories::configure(const PicParticleTrajectoryConfig& config, int seriesEvery,
    const std::filesystem::path& seriesDir, const std::string& moduleTag)
{
    impl_->configure(config, seriesEvery, seriesDir, moduleTag);
}

void PicParticleTrajectories::assignMissingIds(std::vector<ParticleType>& species)
{
    assignMissingParticleTrackIds(species, impl_->nextParticleTrackId);
}

void PicParticleTrajectories::resetOutput()
{
    impl_->requireConfigured();
    if (impl_->particleTrajectoryEnabled) clearParticleTrajectoryRoot(impl_->particleTrajectoryRoot);
}

void PicParticleTrajectories::sample(std::vector<ParticleType>& species, std::mt19937_64& rng,
    int step, double time_s, double lengthScale_m, double velocityScale_m)
{
    impl_->sample(species, rng, step, time_s, lengthScale_m, velocityScale_m);
}

void PicParticleTrajectories::finish() { impl_->finish(); }
