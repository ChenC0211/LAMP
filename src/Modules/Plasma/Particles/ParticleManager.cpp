#include "Modules/Plasma/Particles/ParticleManager.h"
#include "Modules/Plasma/ExplicitPIC/PicConfig.h"
#include "Modules/Plasma/Physics/PhysicalConstants.h"
#include "Modules/Plasma/Boundary/PicEmbeddedGeometry.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/FaceSegmentGraph.h"

namespace {
double sumPhys(const ParticleType& sp, double w_default)
{
    double sum = 0.0;
    for (const auto& p : sp.particle)
        sum += (p.w > 0.0) ? p.w : w_default;
    return sum;
}
ParticleWeightTailStats weightTail(const ParticleType& sp, double w_default)
{
    ParticleWeightTailStats stats;
    if (sp.particle.empty())
        return stats;

    std::vector<double> weights;
    weights.reserve(sp.particle.size());
    double totalPhys = 0.0;
    for (const auto& p : sp.particle)
    {
        const double w = (p.w > 0.0) ? p.w : w_default;
        weights.push_back(w);
        totalPhys += w;
    }

    std::sort(weights.begin(), weights.end());
    auto quantile = [&](double q) -> double {
        if (weights.empty())
            return 0.0;
        const double qq = std::min(1.0, std::max(0.0, q));
        const std::size_t idx = static_cast<std::size_t>(qq * static_cast<double>(weights.size() - 1));
        return weights[idx];
    };

    stats.p50 = quantile(0.50);
    stats.p90 = quantile(0.90);
    stats.p99 = quantile(0.99);
    stats.max = weights.back();

    const std::size_t topCount = std::max<std::size_t>(1u, static_cast<std::size_t>(std::ceil(0.01 * static_cast<double>(weights.size()))));
    double topPhys = 0.0;
    for (std::size_t i = weights.size() - topCount; i < weights.size(); ++i)
        topPhys += weights[i];
    if (totalPhys > 0.0)
        stats.top1PhysFrac = topPhys / totalPhys;
    return stats;
}
}

ParticleManager::Step::Step(const ParticleManagementConfig& config, const QuadtreeMesh& mesh,
    const FaceSegmentGraph& faceGraph, const PicEmbeddedGeometry& embeddedGeometry,
    Mode mode, bool useCutMetrics)
    : config_(config), mesh_(mesh), mode_(mode)
{
    if (mode_ == Mode::Explicit)
    {
        const auto& cutMetrics = embeddedGeometry.metrics();
        if (embeddedGeometry.metricsReady() && cutMetrics.cellMetricsByIndex().size() == static_cast<std::size_t>(mesh.cellIndexNum()))
        {
            const auto& cm = cutMetrics.cellMetricsByIndex();
            skipMask_.assign(static_cast<std::size_t>(mesh.cellIndexNum()), 0u);

            auto markCellLabel = [&](const int cellLabel) {
                if (cellLabel < 0)
                    return;
                const int neighborIndex = mesh.cells(cellLabel).index();
                if (neighborIndex < 0 || neighborIndex >= mesh.cellIndexNum())
                    return;
                skipMask_[static_cast<std::size_t>(neighborIndex)] = 1u;
            };

            for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
            {
                if (!cm[static_cast<std::size_t>(ci)].isCut)
                    continue;

                skipMask_[static_cast<std::size_t>(ci)] = 1u;
                const int cellLabel = mesh.cellIndexToLabel[ci];
                const Cell& cell = mesh.cells(cellLabel);
                for (int face = 0; face < 4; ++face)
                {
                    for (int branch = 0; branch < 2; ++branch)
                        markCellLabel(cell.neiCell[face][branch]);
                }
            }

            maskEnabled_ = true;
        }
    }
    else
    {
        const auto& cutAlphaByCellIndex = embeddedGeometry.alphaByCell();
        if (useCutMetrics && (cutAlphaByCellIndex.size() == static_cast<std::size_t>(mesh.cellIndexNum())))
        {
            skipMask_.assign(static_cast<std::size_t>(mesh.cellIndexNum()), 0u);

            for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
            {
                const double a = cutAlphaByCellIndex[static_cast<std::size_t>(ci)];
                if (!(std::isfinite(a) && (a < 1.0 - 1e-12)))
                    continue;

                skipMask_[static_cast<std::size_t>(ci)] = 1u;
            }

            for (const auto& seg : faceGraph.segments())
            {
                int kIndex = -1;
                int lIndex = -1;
                if (seg.K >= 0)
                    kIndex = mesh.cells(seg.K).index();
                if (seg.L >= 0)
                    lIndex = mesh.cells(seg.L).index();

                const bool kCut = (kIndex >= 0) && (kIndex < mesh.cellIndexNum()) &&
                                  (skipMask_[static_cast<std::size_t>(kIndex)] != 0u);
                const bool lCut = (lIndex >= 0) && (lIndex < mesh.cellIndexNum()) &&
                                  (skipMask_[static_cast<std::size_t>(lIndex)] != 0u);

                if (kCut && lIndex >= 0 && lIndex < mesh.cellIndexNum())
                    skipMask_[static_cast<std::size_t>(lIndex)] = 1u;
                if (lCut && kIndex >= 0 && kIndex < mesh.cellIndexNum())
                    skipMask_[static_cast<std::size_t>(kIndex)] = 1u;
            }

            maskEnabled_ = true;
        }
    }
}

ParticleManagementReport ParticleManager::Step::apply(ParticleType& species, std::size_t speciesIndex,
    double defaultWeight, double velocityScale, std::mt19937_64& rng, bool reportWeights) const
{
    ParticleManagementReport report;
    report.macroBefore = static_cast<long long>(species.particle.size());
    report.physBefore = sumPhys(species, defaultWeight);
    ParticleManager::Apply(config_, mesh_, species, static_cast<int>(speciesIndex),
        species.mass(), velocityScale, rng, report.counts, false, skipMergeCells());
    report.macroAfter = static_cast<long long>(species.particle.size());
    report.physAfter = sumPhys(species, defaultWeight);
    if (mode_ == Mode::Explicit || reportWeights) report.weights = weightTail(species, defaultWeight);
    return report;
}

ParticleManagementStats ParticleManager::Step::applyAll(std::vector<ParticleType>& species,
    const std::vector<double>& defaultWeights, double velocityScale,
    std::mt19937_64& rng, bool reportWeights, const Observer& observer,
    std::vector<long long>* macroLoss, std::vector<double>* physicalLoss) const
{
    ParticleManagementStats totals;
    for (std::size_t si = 0; si < species.size(); ++si)
    {
        auto& sp = species[si];
        const double weight = si < defaultWeights.size() ? defaultWeights[si] : 1.0;
        const auto report = apply(sp, si, weight, velocityScale, rng, reportWeights);
        if (mode_ == Mode::Implicit)
        {
            if (macroLoss && si < macroLoss->size())
                (*macroLoss)[si] += report.macroBefore - report.macroAfter;
            if (physicalLoss && si < physicalLoss->size())
                (*physicalLoss)[si] += report.physBefore - report.physAfter;
        }
        totals.mergedPairs += report.counts.mergedPairs;
        totals.splitCount += report.counts.splitCount;
        totals.splitEmergencyCount += report.counts.splitEmergencyCount;
        totals.splitCountBasedCount += report.counts.splitCountBasedCount;
        if (observer) observer(sp, report);
    }
    return totals;
}

double ParticleManager::kineticEnergyEV(const Particle& particle, double mass_kg, double v0_mps)
{
    if (!(mass_kg > 0.0) || !(v0_mps > 0.0) || !(PhysicalConstants::eV_J > 0.0))
        return 0.0;

    const double vx = particle.velocity.x * v0_mps;
    const double vy = particle.velocity.y * v0_mps;
    const double vz = particle.vz * v0_mps;
    const double v2 = vx * vx + vy * vy + vz * vz;
    return 0.5 * mass_kg * v2 / PhysicalConstants::eV_J;
}

std::size_t ParticleManager::pruneBelowEnergyCutoff(ParticleType& species, double cutoff_eV, double v0_mps)
{
    if (!(cutoff_eV > 0.0) || !(species.mass() > 0.0) || !(v0_mps > 0.0))
    {
        species.number = static_cast<int>(species.particle.size());
        return 0;
    }

    std::size_t write = 0;
    const std::size_t oldSize = species.particle.size();
    for (std::size_t read = 0; read < oldSize; ++read)
    {
        const Particle& particle = species.particle[read];
        if (!(kineticEnergyEV(particle, species.mass(), v0_mps) >= cutoff_eV))
            continue;
        if (write != read)
            species.particle[write] = particle;
        ++write;
    }

    species.particle.resize(write);
    species.number = static_cast<int>(write);
    return oldSize - write;
}

ParticleManager::CullReport ParticleManager::cullEmbedded(const QuadtreeMesh& mesh,
    const ImmersedGeometry<QuadtreeMesh>* iGeometry, const PicEmbeddedGeometry& embeddedGeometry,
    std::vector<ParticleType>& species, const std::vector<double>& macroWeight, bool refreshHosts,
    const std::vector<unsigned char>* speciesMask, std::vector<long long>* macroLossBySpecies,
    std::vector<double>* physLossBySpecies)
{
    if (!embeddedGeometry.ready() || !embeddedGeometry.absorbing() || !iGeometry)
        return {};

    std::size_t removedTotal = 0;
    std::size_t keptTotal = 0;
    for (std::size_t si = 0; si < species.size(); ++si)
    {
        if (speciesMask && (si >= speciesMask->size() || (*speciesMask)[si] == 0u)) continue;
        auto& sp = species[si];
        auto& vec = sp.particle;
        const std::size_t before = vec.size();
        std::size_t write = 0;
        for (std::size_t i = 0; i < before; ++i)
        {
            Particle& p = vec[i];
            if (embeddedGeometry.containsPoint(mesh, iGeometry, p.position, p.inCellL))
            {
                if (macroLossBySpecies != nullptr && si < macroLossBySpecies->size())
                    (*macroLossBySpecies)[si] += 1;
                if (physLossBySpecies != nullptr && si < physLossBySpecies->size())
                {
                    const double wDefault = (si < macroWeight.size()) ? macroWeight[si] : 1.0;
                    (*physLossBySpecies)[si] += (p.w > 0.0) ? p.w : wDefault;
                }
                continue;
            }
            if (refreshHosts)
            {
                int cellL = p.inCellL;
                cellL = mesh.findCellLabelFromHint(cellL, p.position);
                if (cellL < 0)
                    cellL = mesh.findCellLabel(p.position);
                p.inCellL = cellL;
                p.inBackgroundCellL = (cellL >= 0) ? mesh.findBackgroundCellLabel(p.position) : -1;
            }
            vec[write++] = p;
        }
        vec.resize(write);
        sp.number = static_cast<int>(write);
        removedTotal += (before - write);
        keptTotal += write;
    }
    return {removedTotal, keptTotal};
}

std::size_t ParticleManager::pruneByEnergyCutoff(std::vector<ParticleType>& species,
    const std::vector<SpeciesConfig>& settings, const std::vector<double>& macroWeight,
    double v0_mps_energyCutoff, std::vector<long long>* macroLossBySpecies,
    std::vector<double>* physLossBySpecies)
{
    std::size_t removedTotal = 0;
    for (std::size_t si = 0; si < species.size() && si < settings.size(); ++si)
    {
        const double wDefault = (si < macroWeight.size()) ? macroWeight[si] : 1.0;
        const long long macroBefore = static_cast<long long>(species[si].particle.size());
        const double physBefore = sumPhys(species[si], wDefault);
        removedTotal += pruneBelowEnergyCutoff(species[si], settings[si].energyCutoff_eV,
                                                        v0_mps_energyCutoff);
        const long long macroAfter = static_cast<long long>(species[si].particle.size());
        const double physAfter = sumPhys(species[si], wDefault);
        if (macroLossBySpecies != nullptr && si < macroLossBySpecies->size())
            (*macroLossBySpecies)[si] += std::max(0ll, macroBefore - macroAfter);
        if (physLossBySpecies != nullptr && si < physLossBySpecies->size())
            (*physLossBySpecies)[si] += std::max(0.0, physBefore - physAfter);
    }
    return removedTotal;
}
