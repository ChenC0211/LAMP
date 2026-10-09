#include "Modules/Plasma/Physics/PicChargeDeposition.h"
#include "Modules/Plasma/Physics/PicFieldState.h"
#include "Modules/Plasma/Boundary/PicEmbeddedGeometry.h"
#include "Modules/Plasma/Particles/ParticleType.h"
#include "Modules/Plasma/ExplicitPIC/Deposit/CouplingWeightsCache.h"
#include "Modules/Plasma/ExplicitPIC/Deposit/DepositQK.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Units/NormalizationContext.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace {
void requireInputs(const QuadtreeMesh& mesh, const CellField<double>& q, const CellField<double>& rho,
    std::size_t speciesCount, std::size_t cacheCount)
{
    if (q.size() != mesh.cellIndexNum() || rho.size() != mesh.cellIndexNum())
        throw std::logic_error("PIC charge fields must be reset after mesh adaptation");
    if (cacheCount < speciesCount)
        throw std::invalid_argument("PIC charge deposition species/cache count mismatch");
}
}

ExplicitPicChargeDeposition::ExplicitPicChargeDeposition(PicFieldState& fields)
    : fields_(fields), amrNeInstant(fields.mesh_, 0.0) {}

void ExplicitPicChargeDeposition::resetAfterMeshChange()
{
    amrNeInstant.initialize(fields_.mesh_, 0.0);
}

void ExplicitPicChargeDeposition::deposit(const std::vector<ParticleType>& species,
    const std::vector<double>& macroWeight, const std::vector<CouplingWeightsCache>& caches,
    const NormalizationContext& norm, const PicEmbeddedGeometry& embeddedGeometry,
    const ImmersedGeometry<QuadtreeMesh>* geometry, bool amrEnabled, int amrElectronIndex)
{
    const auto& mesh = fields_.mesh_;
    auto& Q = fields_.q_;
    auto& rho = fields_.rho_;
    requireInputs(mesh, Q, rho, species.size(), caches.size());
    if (amrNeInstant.size() != mesh.cellIndexNum())
        throw std::logic_error("Explicit PIC electron amount must be reset after mesh adaptation");
    depositQpScratch.resize(species.size());
    depositWpScratch.resize(species.size());
    const DepositQK depositor;
    const auto& cutAlphaByCellIndex = embeddedGeometry.alphaByCell();
    for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
    {
        Q.value[ci] = 0.0;
        rho.value[ci] = 0.0;
        if (amrEnabled)
            amrNeInstant.value[ci] = 0.0;
    }

    double totalMacroCharge = 0.0;
    for (std::size_t si = 0; si < species.size(); ++si)
    {
        const auto& sp = species[si];
        const auto& cache = caches[si];

        const int n = static_cast<int>(sp.particle.size());
        if (n <= 0)
            continue;

        auto& qp = depositQpScratch[si];
        qp.resize(n);
        const double q = sp.charge();
        const double w_default = (si < macroWeight.size()) ? macroWeight[si] : 1.0;
        if (amrEnabled && static_cast<int>(si) == amrElectronIndex)
        {
            auto& wp = depositWpScratch[si];
            wp.resize(n);
            for (int p = 0; p < n; ++p)
            {
                const double w = (sp.particle[static_cast<std::size_t>(p)].w > 0.0)
                                     ? sp.particle[static_cast<std::size_t>(p)].w
                                     : w_default;
                const double qmacro = q * w;
                qp[p] = qmacro;
                wp[p] = w;
                totalMacroCharge += qmacro;
            }
            (void)depositor.depositDual(qp, wp, cache, &Q, &amrNeInstant);
        }
        else
        {
            for (int p = 0; p < n; ++p)
            {
                const double w = (sp.particle[static_cast<std::size_t>(p)].w > 0.0)
                                     ? sp.particle[static_cast<std::size_t>(p)].w
                                     : w_default;
                const double qmacro = q * w;
                qp[p] = qmacro;
                totalMacroCharge += qmacro;
            }
            (void)depositor.deposit(qp, cache, &Q);
        }
    }

    // Convert deposited charge amount to *code charge density* (rho_code = rho_SI / rho0).
    // Required closure:
    //   V_SI   = VolumeCodeToSI(area_code * alpha_K)   (2D slab closure is inside VolumeCodeToSI)
    //   rho_SI = Q_SI / V_SI
    //   rho_code = rho_SI / rho0

    const bool hasAlpha = embeddedGeometry.metricsReady() && cutAlphaByCellIndex.size() == static_cast<size_t>(mesh.cellIndexNum());
    auto alphaForCell = [&](int ci) -> double {
        if (!hasAlpha)
            return 1.0;
        double alpha = cutAlphaByCellIndex[static_cast<size_t>(ci)];
        if (!std::isfinite(alpha))
            alpha = 1.0;
        if (alpha < 0.0)
            alpha = 0.0;
        else if (alpha > 1.0)
            alpha = 1.0;
        return alpha;
    };
    auto isSolidCell = [&](int ci) -> bool {
        if (hasAlpha)
            return alphaForCell(ci) <= 0.0;
        if (embeddedGeometry.ready())
        {
            const int cLbl = mesh.cellIndexToLabel[ci];
            const Cell& c = mesh.cells(cLbl);
            if (embeddedGeometry.containsPoint(mesh, geometry, c.centerCoord(), cLbl))
                return true;
        }
        return false;
    };

    double totalGridCharge = 0.0;
    for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
    {
        const double alpha = alphaForCell(ci);
        if (isSolidCell(ci) || !(alpha > 0.0))
        {
            Q.value[ci] = 0.0;
            rho.value[ci] = 0.0;
            if (amrEnabled)
                amrNeInstant.value[ci] = 0.0;
            continue;
        }
        const int cLbl = mesh.cellIndexToLabel[ci];
        const Cell& c = mesh.cells(cLbl);
        const double area_code = std::max(1e-300, c.area());

        const double vol_SI = std::max(1e-300, norm.VolumeCodeToSI(area_code * alpha));
        const double rho_SI = Q.value[ci] / vol_SI; // [C/m^3]
        rho.value[ci] = norm.RhoSIToTilde(rho_SI);  // rho_code = rho_SI/rho0

        totalGridCharge += Q.value[ci];
    }

    // Cache totals for console summaries
    totals_.grid = totalGridCharge;
    totals_.particles = totalMacroCharge;

}

void PicFieldState::depositFullCellCharge(const std::vector<ParticleType>& species,
    const std::vector<double>& macroWeight, const std::vector<CouplingWeightsCache>& caches,
    const NormalizationContext& norm)
{
    const auto& mesh = mesh_;
    auto& Qout = q_;
    auto& rhoOut = rho_;
    requireInputs(mesh, Qout, rhoOut, species.size(), caches.size());
    const DepositQK depositor;
    for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
    {
        Qout.value[ci] = 0.0;
        rhoOut.value[ci] = 0.0;
    }

    for (std::size_t si = 0; si < species.size(); ++si)
    {
        const auto& sp = species[si];
        const auto& cache = caches[si];
        const int n = static_cast<int>(sp.particle.size());
        if (n <= 0)
            continue;

        std::vector<double> qp;
        qp.resize(static_cast<std::size_t>(n));
        const double q = sp.charge();
        const double w_default = (si < macroWeight.size()) ? macroWeight[si] : 1.0;
        for (int p = 0; p < n; ++p)
        {
            const double wp = (sp.particle[static_cast<std::size_t>(p)].w > 0.0)
                                  ? sp.particle[static_cast<std::size_t>(p)].w
                                  : w_default;
            qp[static_cast<std::size_t>(p)] = q * wp;
        }
        (void)depositor.deposit(qp, cache, &Qout);
    }

    for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
    {
        const int cLbl = mesh.cellIndexToLabel[ci];
        const Cell& c = mesh.cells(cLbl);
        if (c.son[0] != -1)
        {
            rhoOut.value[ci] = 0.0;
            continue;
        }
        const double area_code = std::max(1e-300, c.area());
        const double vol_SI = std::max(1e-300, norm.VolumeCodeToSI(area_code));
        const double rho_SI = Qout.value[ci] / vol_SI;
        rhoOut.value[ci] = norm.RhoSIToTilde(rho_SI);
    }
}
