#include "Modules/Plasma/IO/PicCheckpointWriter.h"
#include "Modules/Plasma/IO/PicCheckpointStore.h"
#include "Modules/Plasma/Physics/PicCollisions.h"
#include "Modules/Plasma/Physics/PicDomain.h"
#include "Modules/Plasma/Physics/PicFieldState.h"
#include "Modules/Plasma/Boundary/PicMovingGeometry.h"
#include "Modules/Plasma/Particles/ParticleType.h"
#include "Core/IO/FileManagement.h"
#include <stdexcept>
#include <utility>

void PicCheckpointWriter::writeCollisions(const std::filesystem::path& directory,
    const PicCollisions& collisions) const
{
    if (!saveRng_) return;
    const auto state = collisions.rngState();
    const auto path = directory / PicCheckpointStore::collisionName;
    auto out = FileIO::openOutputFile(path, std::ios::out | std::ios::trunc);
    out << "LAMP_MCC_RNG 1\n" << static_cast<int>(state.mode) << ' ' << state.seed << ' ' << state.applyCount << '\n';
    FileIO::finishOutputFile(out, path);
}

PicCheckpointWriter::PicCheckpointWriter(const PicCheckpointStore& store,
    const std::mt19937_64& rng, bool saveRng, int keep, std::string moduleName)
    : store_(store), rng_(rng), saveRng_(saveRng), keep_(keep), moduleName_(std::move(moduleName)) {}

void PicCheckpointWriter::write(int nextStep, double time, double dt, const Components& components) const
{
    if (!components) throw std::invalid_argument("PIC checkpoint components are required");
    std::filesystem::path dir;
    if (!store_.createDirectory(nextStep, dir))
        throw std::runtime_error("[" + moduleName_ + "][Checkpoint] Cannot create directory: " + dir.string());
    store_.beginWrite(dir);
    // Reusing a step directory must not retain a prior collision RNG component
    // when this run disables collisions or saveRng. Fail with the marker intact.
    std::error_code ec;
    std::filesystem::remove(dir / PicCheckpointStore::collisionName, ec);
    if (ec) throw std::runtime_error("Cannot invalidate old MCC RNG state: " + dir.string());
    if (!PicCheckpointStore::writeMetadata(dir, nextStep, time, dt))
        throw std::runtime_error("[" + moduleName_ + "][Checkpoint] Cannot write metadata: " + dir.string());
    if (saveRng_ && !PicCheckpointStore::writeRng(dir, rng_))
        throw std::runtime_error("[" + moduleName_ + "][Checkpoint] Cannot write RNG state: " + dir.string());
    const std::string prefix = PicCheckpointStore::pathPrefix(dir);
    components(dir, prefix);
    store_.completeWrite(dir, saveRng_);
    store_.keepLatest(keep_);
}

void PicCheckpointWriter::write(int nextStep, double time, double dt, const LiveState& state) const
{
    write(nextStep, time, dt, [&](const std::filesystem::path& dir, const std::string& prefix) {
        state.domain.mesh().writeCache(prefix);
        state.domain.writeGeometryCache(prefix);
        state.moving.writeCache(dir);
        if (state.collisions) writeCollisions(dir, *state.collisions);
        for (const auto& species : state.species) species.writeCache(prefix);
        if (state.implicitFields) state.implicitFields->writeCheckpoint(dir);
    });
}
