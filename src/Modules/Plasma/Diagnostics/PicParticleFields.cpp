#include "Modules/Plasma/Diagnostics/PicParticleFields.h"
#include "Modules/Plasma/ExplicitPIC/Deposit/CouplingWeightsCache.h"
#include "Modules/Plasma/Particles/ParticleType.h"
#include "Modules/Plasma/Physics/PhysicalConstants.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Field/CellField.h"
#include "Core/Units/NormalizationContext.h"
#include <algorithm>
#include <cmath>

namespace
{
void depositCSR(const std::vector<double>& qp, const CouplingWeightsCache& cache,
    CellField<double>& outCell)
{
    const auto& idx = cache.cellIndices();
    const auto& w = cache.weights();
    const auto& off = cache.offsets();
    const std::size_t np = qp.size();
    for (std::size_t p = 0; p < np; ++p)
    {
        const double q = qp[p];
        const std::size_t a = static_cast<std::size_t>(off[p]);
        const std::size_t b = static_cast<std::size_t>(off[p + 1]);
        for (std::size_t k = a; k < b; ++k)
        {
            const int ci = idx[k];
            outCell.value[ci] += q * w[k];
        }
    }
}

template<bool FinitePositive, class Emit>
std::size_t sampleElectronEnergy(const ParticleType& species, double w_default,
    const NormalizationContext& norm, const Emit& emit)
{
    const int n = static_cast<int>(species.particle.size());
    const double v0_mps = std::max(1e-300, norm.scales().v0_mps);
    const double m_kg = species.mass();
    std::size_t nBad = 0;
    for (int p = 0; p < n; ++p)
    {
        const auto& part = species.particle[static_cast<std::size_t>(p)];
        double w;
        if constexpr (!FinitePositive)
        {
            w = (part.w > 0.0) ? part.w : w_default;
            if (!std::isfinite(w) || !(w > 0.0))
                w = w_default;
        }
        const double vx = part.velocity.x * v0_mps;
        const double vy = part.velocity.y * v0_mps;
        const double vz = part.vz * v0_mps;
        if (!(std::isfinite(vx) && std::isfinite(vy) && std::isfinite(vz) && std::isfinite(m_kg) && (m_kg > 0.0)))
        {
            if constexpr (FinitePositive) nBad += 1;
            continue;
        }
        const double v2 = vx * vx + vy * vy + vz * vz;
        const double ke_J = 0.5 * m_kg * v2;
        const double ke_eV = ke_J / PhysicalConstants::eV_J;
        if constexpr (FinitePositive)
        {
            if (!std::isfinite(ke_eV))
            {
                nBad += 1;
                continue;
            }
            w = (part.w > 0.0) ? part.w : w_default;
            if (!std::isfinite(w) || !(w > 0.0))
                w = w_default;
            if (!std::isfinite(w) || !(w > 0.0))
            {
                nBad += 1;
                continue;
            }
        }
        else if (!std::isfinite(ke_eV) || !(w > 0.0))
            continue;
        emit(p, w, ke_eV);
    }
    return nBad;
}

void finishElectronEnergy(const QuadtreeMesh& mesh, const CellField<double>& Namount,
    const CellField<double>& Eamount_eV, CellField<double>& Emean_eV, CellField<double>& Te_eV)
{
    for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
    {
        const double denom = Namount.value[ci];
        if (denom > 0.0)
        {
            const double mean_eV = Eamount_eV.value[ci] / denom;
            if (std::isfinite(mean_eV))
            {
                Emean_eV.value[ci] = mean_eV;
                Te_eV.value[ci] = (2.0 / 3.0) * mean_eV;
            }
        }
    }
}
}

namespace PicParticleFields
{
void depositNumberAmount(const ParticleType& species, double w_default,
    const CouplingWeightsCache& cache, CellField<double>& amount)
{
    const int n = static_cast<int>(species.particle.size());
    std::vector<double> qp(static_cast<std::size_t>(n), 0.0);
    for (int p = 0; p < n; ++p)
    {
        const auto& part = species.particle[static_cast<std::size_t>(p)];
        qp[static_cast<std::size_t>(p)] = (part.w > 0.0) ? part.w : w_default;
    }
    depositCSR(qp, cache, amount);
}

void numberDensitySI(const QuadtreeMesh& mesh, const NormalizationContext& norm,
    const CellField<double>& Namount, const std::vector<double>* alphaByCellIndex,
    CellField<double>& nSI)
{
    const bool hasAlpha = alphaByCellIndex && alphaByCellIndex->size() == static_cast<std::size_t>(mesh.cellIndexNum());
    auto alphaForCell = [&](int ci) -> double {
        if (!hasAlpha) return 1.0;
        double alpha = (*alphaByCellIndex)[static_cast<std::size_t>(ci)];
        if (!std::isfinite(alpha)) alpha = 1.0;
        if (alpha < 0.0) alpha = 0.0;
        else if (alpha > 1.0) alpha = 1.0;
        return alpha;
    };
    for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
    {
        const int cLbl = mesh.cellIndexToLabel[ci];
        const Cell& c = mesh.cells(cLbl);
        const double area_code = std::max(1e-300, c.area());
        const double alpha = alphaForCell(ci);
        if (!(alpha > 0.0))
        {
            nSI.value[ci] = 0.0;
            continue;
        }
        const double vol_SI = std::max(1e-300, norm.VolumeCodeToSI(area_code * alpha));
        nSI.value[ci] = Namount.value[ci] / vol_SI;
    }
}

std::size_t electronEnergyCSR(const QuadtreeMesh& mesh, const ParticleType& species,
    double defaultWeight, const NormalizationContext& norm, const CouplingWeightsCache& cache,
    ElectronWeightPolicy policy, CellField<double>& meanEnergyEV, CellField<double>& temperatureEV)
{
    CellField<double> Namount(mesh, 0.0);
    CellField<double> Eamount_eV(mesh, 0.0);
    const int n = static_cast<int>(species.particle.size());
    std::vector<double> qpN(static_cast<std::size_t>(n), 0.0);
    std::vector<double> qpE(static_cast<std::size_t>(n), 0.0);
    auto emit = [&](int p, double w, double ke_eV) {
        qpN[static_cast<std::size_t>(p)] = w;
        qpE[static_cast<std::size_t>(p)] = w * ke_eV;
    };
    const std::size_t skipped = policy == ElectronWeightPolicy::FinitePositive
        ? sampleElectronEnergy<true>(species, defaultWeight, norm, emit)
        : sampleElectronEnergy<false>(species, defaultWeight, norm, emit);
    depositCSR(qpN, cache, Namount);
    depositCSR(qpE, cache, Eamount_eV);
    finishElectronEnergy(mesh, Namount, Eamount_eV, meanEnergyEV, temperatureEV);
    return skipped;
}

void electronEnergyHostCells(const QuadtreeMesh& mesh, const ParticleType& species,
    double defaultWeight, const NormalizationContext& norm,
    CellField<double>& meanEnergyEV, CellField<double>& temperatureEV)
{
    CellField<double> Namount(mesh, 0.0);
    CellField<double> Eamount_eV(mesh, 0.0);
    sampleElectronEnergy<false>(species, defaultWeight, norm, [&](int p, double w, double ke_eV) {
        const auto& part = species.particle[static_cast<std::size_t>(p)];
        int hostLabel = part.inCellL;
        if (hostLabel < 0)
            hostLabel = mesh.findCellLabel(part.position);
        if (hostLabel < 0)
            return;
        const int ci = mesh.cells(hostLabel).index();
        if (ci < 0 || ci >= mesh.cellIndexNum())
            return;
        Namount.value[ci] += w;
        Eamount_eV.value[ci] += w * ke_eV;
    });
    finishElectronEnergy(mesh, Namount, Eamount_eV, meanEnergyEV, temperatureEV);
}

void particleWeights(const QuadtreeMesh& mesh, const std::vector<ParticleType>& species,
    const std::vector<double>& macroWeight, CellField<double>& Wsum,
    CellField<double>& Wmean, CellField<double>& Wmax, CellField<double>& Nmacro)
{
    const int cellCount = mesh.cellIndexNum();
    std::vector<double> wsum(cellCount, 0.0);
    std::vector<double> wmax(cellCount, 0.0);
    std::vector<double> wcount(cellCount, 0.0);
    for (std::size_t si = 0; si < species.size(); ++si)
    {
        const auto& sp = species[si];
        const double w_default = (si < macroWeight.size()) ? macroWeight[si] : 1.0;
        for (const auto& p : sp.particle)
        {
            if (p.inCellL < 0)
                continue;
            const int ci = mesh.cells(p.inCellL).index();
            if (ci < 0 || ci >= cellCount)
                continue;
            const double w = (p.w > 0.0) ? p.w : w_default;
            wsum[ci] += w;
            wcount[ci] += 1.0;
            if (w > wmax[ci])
                wmax[ci] = w;
        }
    }
    for (int ci = 0; ci < cellCount; ++ci)
    {
        Wsum.value[ci] = wsum[ci];
        Wmax.value[ci] = wmax[ci];
        Nmacro.value[ci] = wcount[ci];
        Wmean.value[ci] = (wcount[ci] > 0.0) ? (wsum[ci] / wcount[ci]) : 0.0;
    }
}
}
