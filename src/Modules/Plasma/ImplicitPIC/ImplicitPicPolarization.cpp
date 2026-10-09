#include "Modules/Plasma/ImplicitPIC/ImplicitPicPolarization.h"
#include "Modules/Plasma/ExplicitPIC/Deposit/CouplingWeightsCache.h"
#include "Modules/Plasma/Particles/ParticleType.h"
#include "Modules/Plasma/Physics/PhysicalConstants.h"
#include "Modules/Plasma/Diagnostics/PicStepDiagnostics.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Units/NormalizationContext.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

ImplicitPicPolarization::ImplicitPicPolarization(QuadtreeMesh& mesh,
    std::size_t speciesCount, std::size_t segmentCount)
    : mesh(mesh), Qstar(mesh, 0.0), rhoStar(mesh, 0.0),
      depositQpScratch(speciesCount), depositWpScratch(speciesCount), depositDualWorkspaces(speciesCount)
{
    Nstar.reserve(speciesCount);
    for (std::size_t si = 0; si < speciesCount; ++si)
        Nstar.emplace_back(mesh, 0.0);
    resetCoefficients(segmentCount);
}

void ImplicitPicPolarization::resetDepositedFields()
{
    Qstar.initialize(mesh, 0.0);
    rhoStar.initialize(mesh, 0.0);
    for (auto& nf : Nstar)
        nf.initialize(mesh, 0.0);
}

void ImplicitPicPolarization::resetCoefficients(std::size_t segmentCount)
{
    muCell_Fpm.assign(static_cast<std::size_t>(mesh.cellIndexNum()), 0.0);
    kappaByCell.assign(static_cast<std::size_t>(mesh.cellIndexNum()), 1.0);
    kappaBySeg.assign(segmentCount, 1.0);
}

void ImplicitPicPolarization::resetScratch()
{
    for (auto& buf : depositQpScratch)
        buf.clear();
    for (auto& buf : depositWpScratch)
        buf.clear();
    depositDualWorkspaces.assign(Nstar.size(), {});
}

void ImplicitPicPolarization::requireCurrentSizes() const
{
    const int count = mesh.cellIndexNum();
    if (Qstar.size() != count || rhoStar.size() != count
        || muCell_Fpm.size() != static_cast<std::size_t>(count)
        || kappaByCell.size() != static_cast<std::size_t>(count)
        || std::any_of(Nstar.begin(), Nstar.end(), [&](const auto& field) { return field.size() != count; }))
        throw std::logic_error("Implicit PIC polarization must be rebuilt after mesh adaptation");
}

ImplicitPicPolarization::DepositStats ImplicitPicPolarization::deposit(
    const std::vector<ParticleType>& species, const std::vector<double>& macroWeight,
    const std::vector<CouplingWeightsCache>& cachesBar, int depositFillOmpMinParticles,
    PicStepParticleBalance& chargeBalance)
{
    requireCurrentSizes();
    if (species.size() != Nstar.size() || cachesBar.size() < species.size())
        throw std::invalid_argument("Implicit PIC deposition species/cache count mismatch");
    DepositStats result;
    const DepositQK depositor;
    for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
    {
        Qstar.value[ci] = 0.0;
        rhoStar.value[ci] = 0.0;
        muCell_Fpm[static_cast<std::size_t>(ci)] = 0.0;
        kappaByCell[static_cast<std::size_t>(ci)] = 1.0;
        for (std::size_t si = 0; si < Nstar.size(); ++si)
            Nstar[si].value[ci] = 0.0;
    }

    for (std::size_t si = 0; si < species.size(); ++si)
    {
        const auto& sp = species[si];
        const auto& cache = cachesBar[si];
        const int n = static_cast<int>(sp.particle.size());
        if (n <= 0)
            continue;

        auto& qp = depositQpScratch[si];
        auto& wp = depositWpScratch[si];
        qp.resize(static_cast<std::size_t>(n));
        wp.resize(static_cast<std::size_t>(n));
        const double q = sp.charge();
        const double w_default = (si < macroWeight.size()) ? macroWeight[si] : 1.0;
#ifdef _OPENMP
        if (n >= depositFillOmpMinParticles)
        {
            double qMacroLocal = 0.0;
            #pragma omp parallel for reduction(+:qMacroLocal) schedule(static)
            for (int p = 0; p < n; ++p)
            {
                const double wphys = (sp.particle[static_cast<std::size_t>(p)].w > 0.0)
                                         ? sp.particle[static_cast<std::size_t>(p)].w
                                         : w_default;
                const double qpv = q * wphys;
                qp[static_cast<std::size_t>(p)] = qpv;
                wp[static_cast<std::size_t>(p)] = wphys;
                qMacroLocal += qpv;
            }
            chargeBalance.addStageMacro(qMacroLocal);
        }
        else
#endif
        {
            for (int p = 0; p < n; ++p)
            {
                const double wphys = (sp.particle[static_cast<std::size_t>(p)].w > 0.0)
                                         ? sp.particle[static_cast<std::size_t>(p)].w
                                         : w_default;
                qp[static_cast<std::size_t>(p)] = q * wphys;
                wp[static_cast<std::size_t>(p)] = wphys;
                chargeBalance.addStageMacro(qp[static_cast<std::size_t>(p)]);
            }
        }

        result.particles += static_cast<std::uint64_t>(n);
        result.entries += static_cast<std::uint64_t>(cache.nnz());
        (void)depositor.depositDual(qp, wp, cache, &Qstar, &Nstar[si], &depositDualWorkspaces[si]);
    }
    return result;
}

void ImplicitPicPolarization::assembleCells(const std::vector<ParticleType>& species,
    const NormalizationContext& norm, double dt_si)
{
    requireCurrentSizes();
    if (species.size() != Nstar.size())
        throw std::invalid_argument("Implicit PIC polarization species count mismatch");
    for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
    {
        const int cLbl = mesh.cellIndexToLabel[ci];
        const Cell& c = mesh.cells(cLbl);
        if (c.son[0] != -1)
        {
            rhoStar.value[ci] = 0.0;
            kappaByCell[static_cast<std::size_t>(ci)] = 1.0;
            continue;
        }
        const double area_code = std::max(1e-300, c.area());
        const double vol_SI = std::max(1e-300, norm.VolumeCodeToSI(area_code));
        const double rho_SI = Qstar.value[ci] / vol_SI;
        rhoStar.value[ci] = norm.RhoSIToTilde(rho_SI);

        double mu = 0.0;
        for (std::size_t si = 0; si < species.size(); ++si)
        {
            const double m = species[si].mass();
            const double q = species[si].charge();
            if (!(m > 0.0) || q == 0.0)
                continue;
            const double n_SI = Nstar[si].value[ci] / vol_SI;
            mu += (q * q / m) * n_SI;
        }
        mu *= (dt_si * dt_si) * 0.25;
        muCell_Fpm[static_cast<std::size_t>(ci)] = mu;
        const double kappa = 1.0 + mu / PhysicalConstants::epsilon0_Fpm;
        kappaByCell[static_cast<std::size_t>(ci)] = (kappa > 0.0) ? kappa : 1.0;
    }
}

ImplicitPicPolarization::NeutralityCorrection ImplicitPicPolarization::removeMeanCharge(
    const NormalizationContext& norm)
{
    requireCurrentSizes();
    double totalQ = 0.0;
    double totalVol_SI = 0.0;
    for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
    {
        const int cLbl = mesh.cellIndexToLabel[ci];
        const Cell& c = mesh.cells(cLbl);
        if (c.son[0] != -1)
            continue;
        const double area_code = std::max(1e-300, c.area());
        const double vol_SI = std::max(1e-300, norm.VolumeCodeToSI(area_code));
        totalQ += Qstar.value[ci];
        totalVol_SI += vol_SI;
    }
    if (std::isfinite(totalQ) && std::isfinite(totalVol_SI) && totalVol_SI > 0.0)
    {
        const double meanRho_SI = totalQ / totalVol_SI;
        const double meanRho_tilde = norm.RhoSIToTilde(meanRho_SI);
        for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
        {
            const int cLbl = mesh.cellIndexToLabel[ci];
            const Cell& c = mesh.cells(cLbl);
            if (c.son[0] != -1)
                continue;
            const double area_code = std::max(1e-300, c.area());
            const double vol_SI = std::max(1e-300, norm.VolumeCodeToSI(area_code));
            Qstar.value[ci] -= meanRho_SI * vol_SI;
            rhoStar.value[ci] -= meanRho_tilde;
        }

        return {true, totalQ, meanRho_SI};
    }
    return {};
}

double ImplicitPicPolarization::totalCharge() const
{
    requireCurrentSizes();
    double qGrid = 0.0;
    for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
        qGrid += Qstar.value[ci];
    return qGrid;
}

void ImplicitPicPolarization::assembleFaces(const FaceSegmentGraph& faceGraph)
{
    requireCurrentSizes();
    const auto& segs = faceGraph.segments();
    if (kappaBySeg.size() != segs.size())
        kappaBySeg.resize(segs.size(), 1.0);
    for (std::size_t sid = 0; sid < segs.size(); ++sid)
    {
        const FaceSegment& s = segs[sid];
        const int kIdx = (s.K >= 0) ? mesh.cells(s.K).index() : -1;
        const int lIdx = (s.L >= 0) ? mesh.cells(s.L).index() : -1;
        double kK = 1.0;
        double kL = 1.0;
        if (kIdx >= 0 && static_cast<std::size_t>(kIdx) < kappaByCell.size())
            kK = kappaByCell[static_cast<std::size_t>(kIdx)];
        if (lIdx >= 0 && static_cast<std::size_t>(lIdx) < kappaByCell.size())
            kL = kappaByCell[static_cast<std::size_t>(lIdx)];
        double kf = kK;
        if (lIdx >= 0)
        {
            const double denom = std::max(1e-300, kK + kL);
            kf = 2.0 * kK * kL / denom;
        }
        if (!(kf > 0.0) || !std::isfinite(kf))
            kf = 1.0;
        kappaBySeg[sid] = kf;
    }
}
