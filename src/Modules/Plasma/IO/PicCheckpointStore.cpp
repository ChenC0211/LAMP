#include "Modules/Plasma/IO/PicCheckpointStore.h"
#include "Modules/Plasma/PicCaseInput.h"
#include "Modules/Plasma/Physics/CollisionsConfig.h"
#include "Core/IO/FileManagement.h"
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <utility>

namespace {
bool nonemptyFile(const std::filesystem::path& path)
{
    std::error_code ec;
    if (!std::filesystem::is_regular_file(path, ec) || ec) return false;
    const auto size = std::filesystem::file_size(path, ec);
    return !ec && size > 0;
}
}

PicCheckpointStore::PicCheckpointStore(std::filesystem::path base, Requirements requirements)
    : base_(base.empty() ? std::filesystem::path("../cache") : std::move(base)),
      requirements_(std::move(requirements))
{}

PicCheckpointStore PicCheckpointStore::forRun(PicIntegrator integrator,
    const PicCaseInput& input, const CollisionsConfig* collisions)
{
    const auto& mesh = input.mesh;
    const bool meshBlocksPresent = mesh.geometryRefinement.present || mesh.physicsRefinement.present;
    const bool meshGeoEnabled = meshBlocksPresent ? mesh.geometryRefinement.enabled : false;
    const bool embedded = input.haveGeometry && !input.geometry.immersedBoundaryIds().empty();
    const bool refinement = meshGeoEnabled &&
        (!mesh.geometryRefinement.boundaryRules.empty() ||
         !mesh.geometryRefinement.phaseRules.empty() ||
         !mesh.geometryRefinement.regionRules.empty());
    Requirements requirements;
    for (const auto& species : input.pic.species) requirements.species.push_back(species.name);
    requirements.geometry = input.haveGeometry && (embedded || refinement)
        && (integrator == PicIntegrator::Implicit || !input.geometry.levelSetName.empty());
    requirements.implicitFields = integrator == PicIntegrator::Implicit;
    requirements.movingBoundary = requirements.geometry && input.movingBoundary.enabled && input.movingBoundary.restart;
    requirements.movingBoundaryCache = input.movingBoundary.restartPath;
    requirements.collisionState = input.pic.checkpoint.saveRng && collisions && collisions->enabled
        && collisions->rng_mode == CollisionsRngMode::kCounterBased;
    return PicCheckpointStore(input.pic.checkpoint.path, std::move(requirements));
}

std::filesystem::path PicCheckpointStore::directory(int step) const
{
    return base_ / ("pic." + std::to_string(step));
}

bool PicCheckpointStore::parseStep(const std::string& name, int& step)
{
    const std::string prefix = "pic.";
    if (name.size() <= prefix.size() || name.rfind(prefix, 0) != 0)
        return false;
    const std::string tail = name.substr(prefix.size());
    char* end = nullptr;
    const long value = std::strtol(tail.c_str(), &end, 10);
    if (end == nullptr || *end != '\0')
        return false;
    step = static_cast<int>(value);
    return step >= 0;
}

bool PicCheckpointStore::componentsPresent(const std::filesystem::path& dir) const
{
    std::error_code ec;
    if (dir.empty() || !std::filesystem::is_directory(dir, ec) || ec) return false;
    if (!nonemptyFile(dir / metaName)) return false;
    if (!nonemptyFile(dir / "mesh_cache.bin")) return false;
    for (const auto& name : requirements_.species)
        if (!nonemptyFile(dir / ("particle_" + name + "_cache.bin"))) return false;
    if (requirements_.implicitFields && !nonemptyFile(dir / fieldName)) return false;
    if (requirements_.collisionState && !nonemptyFile(dir / collisionName)) return false;
    if (requirements_.geometry && !nonemptyFile(dir / "iGeometry_cache.bin")) return false;
    if (requirements_.movingBoundary)
    {
        const auto& configured = requirements_.movingBoundaryCache;
        const auto path = configured.empty() ? dir / "moving_eb_cache.bin"
            : configured.is_absolute() ? configured : dir / configured;
        if (!nonemptyFile(path)) return false;
    }
    return true;
}

bool PicCheckpointStore::usable(const std::filesystem::path& dir) const
{
    std::error_code ec;
    if (std::filesystem::exists(dir / incompleteName, ec) || ec) return false;
    return componentsPresent(dir);
}

std::filesystem::path PicCheckpointStore::latest() const
{
    std::error_code ec;
    if (base_.empty() || !std::filesystem::exists(base_, ec)) return {};
    int bestStep = -1;
    std::filesystem::path best;
    for (const auto& entry : std::filesystem::directory_iterator(base_, ec))
    {
        if (ec || !entry.is_directory()) continue;
        int step = -1;
        if (!parseStep(entry.path().filename().string(), step) || !usable(entry.path())) continue;
        if (step > bestStep) { bestStep = step; best = entry.path(); }
    }
    return best;
}

std::filesystem::path PicCheckpointStore::select(int step) const
{
    if (step < 0) return latest();
    const auto dir = directory(step);
    std::error_code ec;
    if (std::filesystem::exists(dir / incompleteName, ec) || ec) return {};
    return dir;
}

bool PicCheckpointStore::createDirectory(int step, std::filesystem::path& dir) const
{
    std::error_code ec;
    if (!base_.empty()) std::filesystem::create_directories(base_, ec);
    ec.clear();
    dir = directory(step);
    std::filesystem::create_directories(dir, ec);
    return !ec;
}

void PicCheckpointStore::beginWrite(const std::filesystem::path& dir) const
{
    const auto marker = dir / incompleteName;
    auto out = FileIO::openOutputFile(marker, std::ios::out | std::ios::trunc);
    out << "incomplete\n";
    FileIO::finishOutputFile(out, marker);
}

void PicCheckpointStore::completeWrite(const std::filesystem::path& dir, bool requireRng) const
{
    if (!componentsPresent(dir) || (requireRng && !nonemptyFile(dir / rngName)))
        throw std::runtime_error("Checkpoint components are missing or empty: " + dir.string());
    std::error_code ec;
    if (!std::filesystem::remove(dir / incompleteName, ec) || ec)
        throw std::runtime_error("Cannot complete checkpoint: " + dir.string());
}

void PicCheckpointStore::keepLatest(int count) const
{
    if (count <= 0) return;
    std::error_code ec;
    if (base_.empty() || !std::filesystem::exists(base_, ec)) return;
    std::vector<int> steps;
    for (const auto& entry : std::filesystem::directory_iterator(base_, ec))
    {
        if (ec || !entry.is_directory()) continue;
        int step = -1;
        if (parseStep(entry.path().filename().string(), step) && usable(entry.path())) steps.push_back(step);
    }
    if (steps.size() <= static_cast<std::size_t>(count)) return;
    std::sort(steps.begin(), steps.end());
    const auto removeCount = steps.size() - static_cast<std::size_t>(count);
    for (std::size_t i = 0; i < removeCount; ++i)
        std::filesystem::remove_all(directory(steps[i]), ec);
}

std::string PicCheckpointStore::pathPrefix(const std::filesystem::path& dir)
{
    std::string prefix = dir.string();
    if (!prefix.empty() && prefix.back() != '/' && prefix.back() != '\\')
        prefix.push_back(std::filesystem::path::preferred_separator);
    return prefix;
}

bool PicCheckpointStore::readMetadataImpl(const std::filesystem::path& dir, Metadata& meta, bool readVersion)
{
    std::ifstream in(dir / metaName, std::ios::in);
    if (!in) return false;
    std::string key;
    while (in >> key)
    {
        if (key == "version" && readVersion) in >> meta.version;
        else if (key == "nextStep" || key == "step") in >> meta.nextStep;
        else if (key == "t_code") in >> meta.t_code;
        else if (key == "dt_code") in >> meta.dt_code;
        else { std::string rest; std::getline(in, rest); }
    }
    return true;
}

bool PicCheckpointStore::readMetadata(const std::filesystem::path& dir, Metadata& meta)
{ return readMetadataImpl(dir, meta, true); }
bool PicCheckpointStore::readOutputMetadata(const std::filesystem::path& dir, Metadata& meta)
{ return readMetadataImpl(dir, meta, false); }

bool PicCheckpointStore::writeMetadata(const std::filesystem::path& dir, int nextStep, double time, double dt)
{
    std::ofstream out(dir / metaName, std::ios::out | std::ios::trunc);
    if (!out) return false;
    out << "version 1\n";
    out << "nextStep " << nextStep << "\n";
    out << std::setprecision(16) << "t_code " << time << "\n";
    out << std::setprecision(16) << "dt_code " << dt << "\n";
    return FileIO::closeOutputFile(out);
}

bool PicCheckpointStore::readRng(const std::filesystem::path& dir, std::mt19937_64& rng)
{
    std::ifstream in(dir / rngName, std::ios::in);
    if (!in) return false;
    in >> rng;
    return true;
}

bool PicCheckpointStore::writeRng(const std::filesystem::path& dir, const std::mt19937_64& rng)
{
    std::ofstream out(dir / rngName, std::ios::out | std::ios::trunc);
    if (!out) return false;
    out << rng;
    return FileIO::closeOutputFile(out);
}
