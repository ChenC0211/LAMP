#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <random>
#include <vector>

#include "Modules/Plasma/Particles/ParticleType.h"
#include "Core/Mesh/BackgroundMesh.h"
#include "Modules/Plasma/Particles/ParticleManagementConfig.h"

// ----------------------------------------------------------------------------
// ParticleManager
//
// Simple per-cell merge/split controls for PIC particle populations.
// Serial implementation intended to preserve existing OpenMP hot paths.
// ----------------------------------------------------------------------------

struct ParticleManagementStats
{
    long long mergedPairs = 0;
    long long splitCount = 0;
    long long splitEmergencyCount = 0;
    long long splitCountBasedCount = 0;
};

class QuadtreeMesh;
class FaceSegmentGraph;
class PicEmbeddedGeometry;
struct SpeciesConfig;
template<class Mesh> class ImmersedGeometry;

struct ParticleWeightTailStats
{
    double p50 = 0., p90 = 0., p99 = 0., max = 0., top1PhysFrac = 0.;
};
struct ParticleManagementReport
{
    long long macroBefore = 0, macroAfter = 0;
    double physBefore = 0., physAfter = 0.;
    ParticleManagementStats counts;
    ParticleWeightTailStats weights;
};

class ParticleManager
{
public:
    class Step;
    struct CullReport { std::size_t removed = 0, kept = 0; };
    // Stable removal on the current geometry. Optional arrays accumulate into
    // existing step totals; only implicit callers refresh surviving hosts.
    static CullReport cullEmbedded(const QuadtreeMesh& mesh,
        const ImmersedGeometry<QuadtreeMesh>* geometry, const PicEmbeddedGeometry& embedded,
        std::vector<ParticleType>& species, const std::vector<double>& defaultWeights,
        bool refreshHosts, const std::vector<unsigned char>* speciesMask = nullptr,
        std::vector<long long>* macroLoss = nullptr, std::vector<double>* physicalLoss = nullptr);
    static double kineticEnergyEV(const Particle& particle, double massKg, double velocityScale);
    static std::size_t pruneBelowEnergyCutoff(ParticleType& species, double cutoffEV, double velocityScale);
    static std::size_t pruneByEnergyCutoff(std::vector<ParticleType>& species,
        const std::vector<SpeciesConfig>& config, const std::vector<double>& defaultWeights,
        double velocityScale, std::vector<long long>* macroLoss = nullptr,
        std::vector<double>* physicalLoss = nullptr);
    template <class Mesh>
    static void Apply(const ParticleManagementConfig& cfg,
                      const Mesh& mesh,
                      ParticleType& sp,
                      int speciesIndex,
                      double mass_kg,
                      double v0_mps,
                      std::mt19937_64& rng,
                      ParticleManagementStats& stats,
                      bool refreshCellLabels = true,
                      const std::vector<unsigned char>* skipMergeCellMask = nullptr)
    {
        if (!cfg.enabled)
            return;
        if (!cfg.mergeEnabled && !cfg.splitEnabled)
            return;
        if (sp.particle.empty())
            return;

        const int nCells = mesh.cellIndexNum();
        std::vector<std::vector<int>> bins(static_cast<std::size_t>(nCells));

        // Bin particles by cached cell index. When refreshCellLabels is false,
        // trust the caller-maintained labels and repair only invalid ones.
        for (std::size_t i = 0; i < sp.particle.size(); ++i)
        {
            Particle& p = sp.particle[i];
            int cellL = p.inCellL;
            if (refreshCellLabels)
                cellL = mesh.findCellLabelFromHint(cellL, p.position);
            if (cellL < 0)
                cellL = mesh.findCellLabel(p.position);
            p.inCellL = cellL;
            if (cellL < 0)
                continue;
            const int cellI = mesh.cells(cellL).index();
            if (cellI < 0 || cellI >= nCells)
                continue;
            bins[static_cast<std::size_t>(cellI)].push_back(static_cast<int>(i));
        }

        std::vector<Particle> newParticles;
        newParticles.reserve(sp.particle.size());

        const double v0 = (v0_mps > 0.0) ? v0_mps : 1.0;
        const double invEV = 1.0 / std::max(1e-300, 1.602176634e-19); // J -> eV

        auto mix64 = [](std::uint64_t a, std::uint64_t b) -> std::uint64_t {
            std::uint64_t x = a ^ (b + 0x9e3779b97f4a7c15ULL + (a << 6) + (a >> 2));
            x ^= x >> 33;
            x *= 0xff51afd7ed558ccdULL;
            x ^= x >> 33;
            x *= 0xc4ceb9fe1a85ec53ULL;
            x ^= x >> 33;
            return x;
        };

        auto particleWeight = [](const Particle& p) -> double {
            return (p.w > 0.0) ? p.w : 1.0;
        };

        auto velocityNorm2 = [](const Particle& p) -> double {
            return p.velocity.x * p.velocity.x + p.velocity.y * p.velocity.y + p.vz * p.vz;
        };

        auto cellNeedsSplitWork = [&](const std::vector<int>& idx) -> bool {
            if (!cfg.splitEnabled || idx.empty())
                return false;
            if (cfg.splitByCountEnabled && cfg.perCellMin > 0 && static_cast<int>(idx.size()) < cfg.perCellMin)
                return true;
            if (cfg.splitWMax > 0.0)
            {
                for (const int pi : idx)
                {
                    if (particleWeight(sp.particle[static_cast<std::size_t>(pi)]) > cfg.splitWMax)
                        return true;
                }
            }
            return false;
        };

        auto dominantVelocityAxis = [&](const std::vector<Particle>& particles) -> int {
            double sumAbsVx = 0.0;
            double sumAbsVy = 0.0;
            double sumAbsVz = 0.0;
            for (const auto& p : particles)
            {
                sumAbsVx += std::abs(p.velocity.x);
                sumAbsVy += std::abs(p.velocity.y);
                sumAbsVz += std::abs(p.vz);
            }

            int axis = 0;
            double best = sumAbsVx;
            if (sumAbsVy > best)
            {
                axis = 1;
                best = sumAbsVy;
            }
            if (sumAbsVz > best)
                axis = 2;
            return axis;
        };

        auto velocityComponent = [](const Particle& p, int axis) -> double {
            if (axis == 1)
                return p.velocity.y;
            if (axis == 2)
                return p.vz;
            return p.velocity.x;
        };

        auto mergePair = [&](Particle left, Particle right) -> Particle {
            const double w1 = particleWeight(left);
            const double w2 = particleWeight(right);
            const double wsum = w1 + w2;
            if (wsum > 0.0)
            {
                left.velocity = (left.velocity * w1 + right.velocity * w2) * (1.0 / wsum);
                left.vz = (left.vz * w1 + right.vz * w2) / wsum;
                left.w = wsum;
            }

            const std::uint64_t a = left.rngState;
            const std::uint64_t b = right.rngState;
            const std::uint64_t r = static_cast<std::uint64_t>(rng());
            left.rngState = mix64(a == 0 ? r : a, b == 0 ? r : b);
            return left;
        };

        auto shouldSkipPmCell = [&](const int cellI) -> bool {
            if (skipMergeCellMask == nullptr)
                return false;
            if (cellI < 0 || cellI >= nCells)
                return false;
            if (skipMergeCellMask->size() != static_cast<std::size_t>(nCells))
                return false;
            return (*skipMergeCellMask)[static_cast<std::size_t>(cellI)] != 0u;
        };

        // Helper for splitting positions inside a cell box.
        auto placeInsideCell = [&](Particle& p, const Cell& cell, const Vector2D& base, const Vector2D& offset) {
            const Vector2D center = cell.centerCoord();
            const Vector2D size = cell.size();
            const double hx = 0.5 * size.x;
            const double hy = 0.5 * size.y;
            Vector2D pos(base.x + offset.x, base.y + offset.y);
            pos.x = std::min(std::max(pos.x, center.x - hx), center.x + hx);
            pos.y = std::min(std::max(pos.y, center.y - hy), center.y + hy);
            p.position = pos;
            p.inCellL = cell.label();
        };

        for (int cellI = 0; cellI < nCells; ++cellI)
        {
            const int cellL = mesh.cellIndexToLabel[cellI];
            if (cellL < 0)
                continue;
            const Cell& cell = mesh.cells(cellL);
            const auto& idx = bins[static_cast<std::size_t>(cellI)];
            if (idx.empty())
                continue;

            const bool skipPmCell = shouldSkipPmCell(cellI);
            const bool needsMergeWork = !skipPmCell && cfg.mergeEnabled && cfg.perCellMax > 0 &&
                                        static_cast<int>(idx.size()) > cfg.perCellMax;
            const bool needsSplitWork = !skipPmCell && cellNeedsSplitWork(idx);

            if (!needsMergeWork && !needsSplitWork)
            {
                for (const int pi : idx)
                    newParticles.push_back(std::move(sp.particle[static_cast<std::size_t>(pi)]));
                continue;
            }

            std::vector<Particle> cellParticles;
            cellParticles.reserve(idx.size());
            for (const int pi : idx)
                cellParticles.push_back(std::move(sp.particle[static_cast<std::size_t>(pi)]));

            // Merge: reduce counts above perCellMax.
            if (needsMergeWork)
            {
                int merges = 0;
                double maxDeltaE = cfg.mergeMaxDeltaE_eV;
                if (speciesIndex >= 0 &&
                    static_cast<std::size_t>(speciesIndex) < cfg.mergeMaxDeltaE_bySpecies_eV.size())
                {
                    const double v = cfg.mergeMaxDeltaE_bySpecies_eV[static_cast<std::size_t>(speciesIndex)];
                    if (v >= 0.0)
                        maxDeltaE = v;
                }
                maxDeltaE = std::max(0.0, maxDeltaE);
                while (static_cast<int>(cellParticles.size()) > cfg.perCellMax &&
                       merges < cfg.mergeMaxPairsPerCell &&
                       cellParticles.size() >= 2)
                {
                    const int dominantAxis = dominantVelocityAxis(cellParticles);
                    std::array<std::vector<Particle>, 3> groups;
                    for (const auto& p : cellParticles)
                    {
                        const double vc = velocityComponent(p, dominantAxis);
                        const int group = (vc > 0.0) ? 2 : ((vc < 0.0) ? 0 : 1);
                        groups[static_cast<std::size_t>(group)].push_back(p);
                    }

                    for (auto& group : groups)
                    {
                        std::stable_sort(group.begin(), group.end(), [&](const Particle& lhs, const Particle& rhs) {
                            const double lhsVc = velocityComponent(lhs, dominantAxis);
                            const double rhsVc = velocityComponent(rhs, dominantAxis);
                            if (lhsVc != rhsVc)
                                return lhsVc < rhsVc;

                            const double lhsSpeed2 = velocityNorm2(lhs);
                            const double rhsSpeed2 = velocityNorm2(rhs);
                            if (lhsSpeed2 != rhsSpeed2)
                                return lhsSpeed2 < rhsSpeed2;

                            return particleWeight(lhs) < particleWeight(rhs);
                        });
                    }

                    std::vector<Particle> mergedParticles;
                    mergedParticles.reserve(cellParticles.size());
                    int mergesThisRound = 0;

                    for (auto& group : groups)
                    {
                        for (std::size_t i = 0; i < group.size();)
                        {
                            if (static_cast<int>(mergedParticles.size() + (group.size() - i)) <= cfg.perCellMax ||
                                merges >= cfg.mergeMaxPairsPerCell)
                            {
                                for (; i < group.size(); ++i)
                                    mergedParticles.push_back(group[i]);
                                break;
                            }

                            if (i + 1 >= group.size())
                            {
                                mergedParticles.push_back(group[i]);
                                break;
                            }

                            const Particle& p1 = group[i];
                            const Particle& p2 = group[i + 1];
                            const double dvx = (p1.velocity.x - p2.velocity.x) * v0;
                            const double dvy = (p1.velocity.y - p2.velocity.y) * v0;
                            const double dvz = (p1.vz - p2.vz) * v0;
                            const double dE_eV = 0.5 * mass_kg * (dvx * dvx + dvy * dvy + dvz * dvz) * invEV;

                            if (dE_eV <= maxDeltaE)
                            {
                                mergedParticles.push_back(mergePair(p1, p2));
                                merges += 1;
                                mergesThisRound += 1;
                                i += 2;
                            }
                            else
                            {
                                mergedParticles.push_back(p1);
                                i += 1;
                            }
                        }
                    }

                    if (mergesThisRound == 0)
                        break;

                    cellParticles.swap(mergedParticles);
                }
                stats.mergedPairs += merges;
            }

            // Split: emergency-only fallback for overweight particles.
            if (needsSplitWork)
            {
                int splits = 0;
                int emergencySplits = 0;
                int countBasedSplits = 0;
                const double jitterFrac = std::max(0.0, cfg.splitJitterFrac);
                const Vector2D baseSize = cell.size();
                const double jx = jitterFrac * baseSize.x;
                const double jy = jitterFrac * baseSize.y;

                for (std::size_t i = 0; i < cellParticles.size() && splits < cfg.splitMaxPerCell; ++i)
                {
                    const bool needByCount = cfg.splitByCountEnabled && (cfg.perCellMin > 0) &&
                                             (static_cast<int>(cellParticles.size()) < cfg.perCellMin);
                    const bool needByWeight = (cfg.splitWMax > 0.0) &&
                                              (particleWeight(cellParticles[i]) > cfg.splitWMax);
                    if (!needByCount && !needByWeight)
                        continue;

                    if (!needByCount && static_cast<int>(cellParticles.size()) < 2)
                        continue;

                    Particle base = cellParticles[i];
                    const double w = particleWeight(base);
                    const double w2 = w * 0.5;
                    if (!(w2 > 0.0))
                        continue;

                    const int dominantAxis = dominantVelocityAxis(cellParticles);
                    Vector2D offsetA(0.0, 0.0);
                    Vector2D offsetB(0.0, 0.0);
                    if (dominantAxis == 0)
                    {
                        const double dy = (jy > 0.0) ? jy : 0.25 * baseSize.y;
                        offsetA = Vector2D(0.0, +dy);
                        offsetB = Vector2D(0.0, -dy);
                    }
                    else
                    {
                        const double dx = (jx > 0.0) ? jx : 0.25 * baseSize.x;
                        offsetA = Vector2D(+dx, 0.0);
                        offsetB = Vector2D(-dx, 0.0);
                    }

                    const std::uint64_t seed = (base.rngState != 0u)
                        ? base.rngState
                        : static_cast<std::uint64_t>(rng());

                    base.w = w2;
                    base.rngState = mix64(seed, 0x9e3779b97f4a7c15ULL);
                    placeInsideCell(base, cell, base.position, offsetA);
                    cellParticles[i] = base;

                    Particle np = base;
                    np.w = w2;
                    np.rngState = mix64(seed, 0xc4ceb9fe1a85ec53ULL);
                    placeInsideCell(np, cell, base.position, offsetB);
                    cellParticles.push_back(np);

                    splits += 1;
                    if (needByWeight)
                        emergencySplits += 1;
                    else if (needByCount)
                        countBasedSplits += 1;
                }
                stats.splitCount += splits;
                stats.splitEmergencyCount += emergencySplits;
                stats.splitCountBasedCount += countBasedSplits;
            }

            for (auto& p : cellParticles)
                newParticles.push_back(std::move(p));
        }

        sp.particle.swap(newParticles);
        sp.number = static_cast<int>(sp.particle.size());
    }
};

// One resampling stage on the current mesh. Owns its protected-cell mask;
// mesh/configuration must outlive this scope. Recreate after any mesh/geometry
// change. Particle mutation still uses the original Apply kernel above.
class ParticleManager::Step
{
public:
    enum class Mode { Explicit, Implicit };
    Step(const ParticleManagementConfig& config, const QuadtreeMesh& mesh,
        const FaceSegmentGraph& faces, const PicEmbeddedGeometry& embedded,
        Mode mode, bool useCutMetrics = false);
    Step(const Step&) = delete;
    Step& operator=(const Step&) = delete;
    ParticleManagementReport apply(ParticleType& species, std::size_t speciesIndex,
        double defaultWeight, double velocityScale, std::mt19937_64& rng, bool reportWeights) const;
    using Observer = std::function<void(const ParticleType&, const ParticleManagementReport&)>;
    // Process species in storage order. Implicit mode accumulates signed
    // before-after population changes in the supplied, possibly shorter arrays;
    // explicit mode leaves them untouched. Missing default weights use 1.
    // Observe each species immediately after mutation and accounting, including
    // empty species. The synchronous observer must not mutate the collection;
    // no references are retained. An exception preserves completed mutations.
    ParticleManagementStats applyAll(std::vector<ParticleType>& species,
        const std::vector<double>& defaultWeights, double velocityScale,
        std::mt19937_64& rng, bool reportWeights, const Observer& observer = {},
        std::vector<long long>* macroLoss = nullptr, std::vector<double>* physicalLoss = nullptr) const;
    const std::vector<unsigned char>* skipMergeCells() const { return maskEnabled_ ? &skipMask_ : nullptr; }
private:
    const ParticleManagementConfig& config_;
    const QuadtreeMesh& mesh_;
    const Mode mode_;
    bool maskEnabled_ = false;
    std::vector<unsigned char> skipMask_;
};
