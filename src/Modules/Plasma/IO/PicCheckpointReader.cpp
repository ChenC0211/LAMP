#include "Modules/Plasma/IO/PicCheckpointReader.h"
#include "Modules/Plasma/Physics/PicCollisions.h"
#include "Modules/Plasma/Particles/ParticleType.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <utility>
#include <charconv>
#include <fstream>
#include <stdexcept>

namespace {
std::optional<PicCollisions::RngState> readCollisionState(const std::filesystem::path& directory,
    const PicCollisions::RngState& expected)
{
    const auto path = directory / PicCheckpointStore::collisionName;
    auto state = expected;
    std::error_code ec;
    if (!std::filesystem::exists(path, ec) && !ec && state.mode != CollisionsRngMode::kCounterBased)
        return std::nullopt; // Old global/per-particle streams do not use this counter.
    auto invalid = [&]() {
        return std::runtime_error("[PIC][Checkpoint] Missing or invalid MCC RNG state: "
            + path.string() + " (counterBased restart requires the saved application count)");
    };
    std::ifstream in(path);
    std::string magic, version, mode, seed, count, extra;
    if (!(in >> magic >> version >> mode >> seed >> count) || (magic != "LAMP_MCC_RNG" && magic != "QUANTUM_MCC_RNG") || version != "1"
        || (in >> extra) || in.bad()) throw invalid();
    auto number = [&](const std::string& token) {
        std::uint64_t value = 0;
        const auto parsed = std::from_chars(token.data(), token.data() + token.size(), value);
        if (parsed.ec != std::errc{} || parsed.ptr != token.data() + token.size()) throw invalid();
        return value;
    };
    const auto modeValue = number(mode);
    if (modeValue > static_cast<std::uint64_t>(CollisionsRngMode::kCounterBased)) throw invalid();
    state.mode = static_cast<CollisionsRngMode>(modeValue);
    state.seed = number(seed);
    state.applyCount = number(count);
    if (state.mode != expected.mode || state.seed != expected.seed)
        throw std::runtime_error("[PIC][Checkpoint] MCC RNG mode/seed mismatch: " + path.string());
    return state;
}
}

void PicCheckpointReader::readCollisions(PicCollisions& collisions, bool restoreRng) const
{
    if (!restoreRng) return;
    if (const auto state = readCollisionState(directory_, collisions.rngState()))
        collisions.operation_->restoreRngState(*state);
}

void PicCheckpointReader::validateCollisions(const std::filesystem::path& directory,
    const CollisionsConfig& config, bool restoreRng)
{
    if (restoreRng && config.enabled)
        (void)readCollisionState(directory, {config.rng_mode, config.rng_seed, 0});
}

PicCheckpointReader::PicCheckpointReader(std::filesystem::path directory,
    PicCheckpointStore::Metadata metadata, std::string moduleName)
    : directory_(std::move(directory)), metadata_(metadata), moduleName_(std::move(moduleName)) {}

std::optional<PicCheckpointReader> PicCheckpointReader::select(const PicCheckpointStore& store,
    int step, double& dt, std::string moduleName)
{
    return readSelected(store.select(step), dt, std::move(moduleName));
}

std::optional<PicCheckpointReader> PicCheckpointReader::readSelected(std::filesystem::path directory,
    double& dt, std::string moduleName)
{
    PicCheckpointStore::Metadata metadata;
    std::error_code error;
    if (directory.empty() || std::filesystem::exists(directory / PicCheckpointStore::incompleteName, error)
        || error || !PicCheckpointStore::readMetadata(directory, metadata))
    {
        std::cerr << '[' << moduleName << "][Checkpoint] ERROR: restart requested but no checkpoint found.\n";
        return std::nullopt;
    }
    if (metadata.dt_code > 0.0 && std::abs(metadata.dt_code - dt) > 1e-16)
    {
        std::cout << '[' << moduleName << "][Checkpoint] NOTE: override dt from pic.dict ("
                  << dt << ") to checkpoint dt (" << metadata.dt_code << ").\n";
        dt = metadata.dt_code;
    }
    return PicCheckpointReader(std::move(directory), metadata, std::move(moduleName));
}

int PicCheckpointReader::nextStep() const noexcept
{
    return std::max(0, metadata_.nextStep);
}

std::optional<PicCheckpointReader::PreparedExplicitRestart> PicCheckpointReader::prepareExplicitRestart(const PicCheckpointStore& store,
    int step, const CollisionsConfig* collisions, bool restoreRng)
{
    const auto directory = store.select(step);
    PicCheckpointStore::Metadata metadata;
    if (directory.empty() || !PicCheckpointStore::readOutputMetadata(directory, metadata))
    {
        std::cerr << "[ExplicitPIC][Checkpoint] ERROR: no compatible checkpoint found; check required components"
                     " (counterBased with saveRng requires pic_collisions_state.txt).\n";
        return std::nullopt;
    }
    if (collisions) validateCollisions(directory, *collisions, restoreRng);
    return PreparedExplicitRestart(directory, std::max(0, metadata.nextStep));
}

PicCheckpointReader::PreparedExplicitRestart::PreparedExplicitRestart(
    std::filesystem::path directory, int outputResumeStep)
    : directory_(std::move(directory)), outputResumeStep_(outputResumeStep) {}

std::optional<PicCheckpointReader> PicCheckpointReader::PreparedExplicitRestart::read(double& dt) const
{
    return PicCheckpointReader::readSelected(directory_, dt, "ExplicitPIC");
}

bool PicCheckpointReader::readMesh(QuadtreeMesh& mesh, MeshCache policy) const
{
    const std::string prefix = PicCheckpointStore::pathPrefix(directory_);
    const std::filesystem::path meshCachePath = directory_ / "mesh_cache.bin";
    std::error_code ec;
    if (!std::filesystem::exists(meshCachePath, ec))
    {
        if (policy == MeshCache::Optional) return true;
        std::cerr << '[' << moduleName_ << "][Checkpoint] ERROR: missing mesh cache: "
                  << meshCachePath.string() << std::endl;
        return false;
    }
    mesh.readCache(prefix);
    return true;
}

bool PicCheckpointReader::readParticles(const QuadtreeMesh& mesh,
    std::vector<ParticleType>& species, std::mt19937_64& rng, bool restoreRng) const
{
    if (restoreRng)
    {
        const std::filesystem::path rngPath = directory_ / PicCheckpointStore::rngName;
        if (!PicCheckpointStore::readRng(directory_, rng))
            std::cerr << '[' << moduleName_ << "][Checkpoint] WARNING: RNG state not found: "
                      << rngPath.string() << std::endl;
    }
    std::error_code ec;
    if (directory_.empty() || !std::filesystem::exists(directory_, ec)) return false;
    const std::string prefix = PicCheckpointStore::pathPrefix(directory_);
    for (auto& sp : species)
    {
        const std::filesystem::path filePath = directory_ / ("particle_" + sp.name() + "_cache.bin");
        if (!std::filesystem::exists(filePath, ec))
        {
            std::cerr << '[' << moduleName_ << "][Checkpoint] ERROR: missing particle cache: "
                      << filePath.string() << std::endl;
            return false;
        }
        sp.readCache(prefix);
        sp.number = static_cast<int>(sp.particle.size());
    }
    for (auto& sp : species)
    {
        for (auto& p : sp.particle)
        {
            int cellL = p.inCellL;
            cellL = mesh.findCellLabelFromHint(cellL, p.position);
            if (cellL < 0)
                cellL = mesh.findCellLabel(p.position);
            p.inCellL = cellL;
            p.inBackgroundCellL = (cellL >= 0) ? mesh.findBackgroundCellLabel(p.position) : -1;
        }
        sp.number = static_cast<int>(sp.particle.size());
    }
    return true;
}

void PicCheckpointReader::report() const
{
    std::string directory = directory_.string();
    for (char& c : directory)
        if (c == '\\') c = '/';
    std::cout << '[' << moduleName_ << "][Checkpoint] Restart from step=" << nextStep()
              << " (" << directory << ")" << std::endl;
}
