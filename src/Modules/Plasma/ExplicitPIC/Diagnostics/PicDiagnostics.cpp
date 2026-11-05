#include "Modules/Plasma/ExplicitPIC/Diagnostics/PicDiagnostics.h"

#include <algorithm>
#include <cctype>
#include <utility>
#include "Core/Units/NormalizationContext.h"
#include "Modules/Plasma/Particles/ParticleType.h"
#include "Modules/Plasma/Physics/PhysicalConstants.h"
#include <cmath>
#include <limits>
#include <vector>

#include "Core/Mesh/PeriodicFaceUtils.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/Cell.h"
#include "Solvers/Poisson/CutCellMetrics.h"

namespace ExplicitPIC
{
namespace Diagnostics
{
namespace
{
inline const PicBoundarySide& bcSide(const PicDomainBC& bc, const Boundary b)
{
    switch (b)
    {
    case Boundary::west:  return bc.left;
    case Boundary::east:  return bc.right;
    case Boundary::south: return bc.bottom;
    case Boundary::north: return bc.top;
    default:              return bc.left;
    }
}

inline bool isDirichletLike(PicBCType t)
{
    return (t == PicBCType::Dirichlet) || (t == PicBCType::Robin);
}

inline double safeCellArea(const QuadtreeMesh& mesh, int cellIndex)
{
    if (cellIndex < 0 || cellIndex >= mesh.cellIndexNum())
        return 0.0;
    const int lab = mesh.cellIndexToLabel[cellIndex];
    return mesh.cells(lab).area();
}

// Coarse-fine composite band predicate, kept consistent with CouplingWeightsCache.cpp.
inline bool isCompositeBandCell(const QuadtreeMesh& mesh, const Cell& c)
{
    // Non-leaf treated as in-band.
    if (!(c.son[0] == -1 && c.son[1] == -1 && c.son[2] == -1 && c.son[3] == -1))
        return true;

    const int lvl = c.level();
    for (int dir = 0; dir < 4; ++dir)
    {
        const int nb0 = c.neiCell[dir][0];
        const int nb1 = c.neiCell[dir][1];
        if (nb1 != -1)
            return true;
        if (nb0 != -1)
        {
            const Cell& nb = mesh.cells(nb0);
            const bool nbLeaf = (nb.son[0] == -1 && nb.son[1] == -1 && nb.son[2] == -1 && nb.son[3] == -1);
            if (!nbLeaf || nb.level() != lvl)
                return true;
        }
    }
    return false;
}
} // namespace

namespace
{
void appendParticleSummaries(PicDiagnostics::Result& diag, const QuadtreeMesh& mesh,
    const std::vector<ParticleType>& species, const std::vector<double>& macroWeight,
    const NormalizationContext& norm, const Solvers::Poisson::CutCellMetrics* cutMetricsPtr,
    bool allChildren)
{
    double Vdomain_SI = 0.0;
    for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
    {
        const int cLbl = mesh.cellIndexToLabel[ci];
        const Cell& c = mesh.cells(cLbl);
        if (c.son[0] != -1 || (allChildren && (c.son[1] != -1 || c.son[2] != -1 || c.son[3] != -1)))
            continue;
        const double area_code = std::max(1e-300, c.area());
        double alpha = 1.0;
        if (cutMetricsPtr)
        {
            const auto& cm = cutMetricsPtr->cellMetricsByIndex();
            if (static_cast<size_t>(ci) < cm.size())
            {
                alpha = cm[static_cast<size_t>(ci)].alpha;
                if (!std::isfinite(alpha))
                    alpha = 1.0;
                if (alpha < 0.0)
                    alpha = 0.0;
                else if (alpha > 1.0)
                    alpha = 1.0;
            }
        }
        Vdomain_SI += norm.VolumeCodeToSI(area_code * alpha);
    }
    Vdomain_SI = std::max(1e-300, Vdomain_SI);

    diag.speciesNames.clear();
    diag.speciesNmean_1pm3.clear();
    diag.speciesMacroCount.clear();
    diag.speciesPhysCount.clear();
    diag.macro_total = 0;
    diag.phys_total = 0.0;
    diag.speciesNames.reserve(species.size());
    diag.speciesNmean_1pm3.reserve(species.size());
    diag.speciesMacroCount.reserve(species.size());
    diag.speciesPhysCount.reserve(species.size());
    for (std::size_t si = 0; si < species.size(); ++si)
    {
        diag.speciesNames.push_back(species[si].name());
        double Nphys = 0.0;
        const long long Nmacro = static_cast<long long>(species[si].particle.size());
        const double w_default = (si < macroWeight.size()) ? macroWeight[si] : 1.0;
        for (const auto& p : species[si].particle)
            Nphys += (p.w > 0.0) ? p.w : w_default;
        diag.speciesNmean_1pm3.push_back(Nphys / Vdomain_SI);
        diag.speciesMacroCount.push_back(Nmacro);
        diag.speciesPhysCount.push_back(Nphys);
        diag.macro_total += Nmacro;
        diag.phys_total += Nphys;
    }

    // Select electrons: prefer name "e"/"electron", else any negative charge.
    int eIndex = -1;
    for (std::size_t si = 0; si < species.size(); ++si)
    {
        auto nm = species[si].name();
        for (auto& ch : nm)
            ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
        if (nm == "e" || nm == "electron" || nm.find("electron") != std::string::npos)
        {
            eIndex = static_cast<int>(si);
            break;
        }
    }
    if (eIndex < 0)
    {
        for (std::size_t si = 0; si < species.size(); ++si)
        {
            if (species[si].charge() < 0.0)
            {
                eIndex = static_cast<int>(si);
                break;
            }
        }
    }

    if (eIndex >= 0)
    {
        const std::size_t si = static_cast<std::size_t>(eIndex);
        const double m_kg = species[si].mass();
        const double v0_mps = std::max(1e-300, norm.scales().v0_mps);

        const double w_default = (si < macroWeight.size()) ? macroWeight[si] : 1.0;
        double Nphys = 0.0;
        double Etot_eV = 0.0;
        for (const auto& p : species[si].particle)
        {
            const double vx = p.velocity.x * v0_mps;
            const double vy = p.velocity.y * v0_mps;
            const double vz = p.vz * v0_mps;
            const double v2 = vx * vx + vy * vy + vz * vz;
            const double ke_J = 0.5 * m_kg * v2;
            const double w = (p.w > 0.0) ? p.w : w_default;
            Nphys += w;
            Etot_eV += (ke_J / PhysicalConstants::eV_J) * w;
        }
        diag.electron_mean_energy_eV = (Nphys > 0.0) ? (Etot_eV / Nphys) : 0.0;
        diag.Te_e_mean_eV = (2.0 / 3.0) * diag.electron_mean_energy_eV;
    }
}
}

PicDiagnostics::ExplicitSample PicDiagnostics::sampleExplicit(const QuadtreeMesh& mesh,
    const FaceSegmentGraph& graph, const PicDomainBC& domainBC,
    const CellField<double>& rhs, const CellField<double>& phi, const CellField<double>& Q,
    double totalMacroCharge, const std::vector<ParticleType>& species,
    const std::vector<double>& macroWeight, const NormalizationContext& norm,
    const Solvers::Poisson::CutCellMetrics* cutMetricsPtr,
    const PicEmbeddedBC& embeddedBC, const PicEmbeddedBCMap& embeddedBCById)
{
    double phiMin = std::numeric_limits<double>::infinity();
    double phiMax = -std::numeric_limits<double>::infinity();
    for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
    {
        phiMin = std::min(phiMin, phi.value[ci]);
        phiMax = std::max(phiMax, phi.value[ci]);
    }

    auto diag = compute(mesh, graph, domainBC, rhs, phi, &Q, totalMacroCharge,
        nullptr, cutMetricsPtr, &embeddedBC, &embeddedBCById);

    appendParticleSummaries(diag, mesh, species, macroWeight, norm, cutMetricsPtr, false);
    return {std::move(diag), phiMin, phiMax};
}

PicDiagnostics::Result PicDiagnostics::sampleImplicit(const QuadtreeMesh& mesh,
    const FaceSegmentGraph& graph, const PicDomainBC& domainBC,
    const CellField<double>& rhs, const CellField<double>& phi, const CellField<double>& Q,
    double totalMacroCharge, const std::vector<ParticleType>& species,
    const std::vector<double>& macroWeight, const NormalizationContext& norm,
    const Solvers::Poisson::CutCellMetrics* cutMetrics,
    const PicEmbeddedBC* embeddedBC, const PicEmbeddedBCMap* embeddedBCById,
    const std::vector<double>& kappaBySegment, const std::vector<double>& kappaByCellIndex)
{
    auto diag = compute(mesh, graph, domainBC, rhs, phi, &Q, totalMacroCharge,
        nullptr, cutMetrics, embeddedBC, embeddedBCById, &kappaBySegment, &kappaByCellIndex);
    appendParticleSummaries(diag, mesh, species, macroWeight, norm, cutMetrics, true);
    return diag;
}

PicDiagnostics::Result PicDiagnostics::compute(const QuadtreeMesh& mesh,
                                              const FaceSegmentGraph& graph,
                                              const PicDomainBC& domainBC,
                                              const CellField<double>& rhs,
                                              const CellField<double>& phi,
                                              const CellField<double>* Q,
                                              double Qmacro,
                                              CellField<double>* rGout,
                                              const Solvers::Poisson::CutCellMetrics* cutMetrics,
                                              const PicEmbeddedBC* embeddedBC,
                                              const PicEmbeddedBCMap* embeddedBCById,
                                              const std::vector<double>* kappaBySegment,
                                              const std::vector<double>* kappaByCellIndex)
{
    Result out;
    out.Qmacro = Qmacro;

    // Active unknown set: cells with at least one segment.
    const int nCellIndex = mesh.cellIndexNum();
    std::vector<int> rowOfCellIndex(static_cast<size_t>(nCellIndex), -1);
    std::vector<int> cellIndexOfRow;
    cellIndexOfRow.reserve(static_cast<size_t>(nCellIndex));

    for (int ci = 0; ci < nCellIndex; ++ci)
    {
        const auto range = graph.segmentRangeForCellIndex(ci);
        if (range.first < range.second)
        {
            rowOfCellIndex[static_cast<size_t>(ci)] = static_cast<int>(cellIndexOfRow.size());
            cellIndexOfRow.push_back(ci);
        }
    }

    const int N = static_cast<int>(cellIndexOfRow.size());
    out.nActive = N;

    if (rGout)
    {
        if (rGout->size() != mesh.cellIndexNum() || !rGout->value)
            rGout->initialize(mesh, 0.0);
        for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
            rGout->value[ci] = 0.0;
    }

    if (N == 0)
        return out;

    const bool useCutMetrics =
        (cutMetrics != nullptr) &&
        (cutMetrics->betaBySegment().size() == graph.segments().size()) &&
        (cutMetrics->cellMetricsByIndex().size() == static_cast<size_t>(mesh.cellIndexNum()));

    const bool useEmbeddedBCMap =
        useCutMetrics && (embeddedBCById != nullptr) && (!embeddedBCById->empty());
    const bool useEmbeddedBCFallback =
        useCutMetrics && (embeddedBC != nullptr) && embeddedBC->enabled;
    const bool useEmbeddedBC = useEmbeddedBCMap || useEmbeddedBCFallback;

    const bool useKappaBySegment =
        (kappaBySegment != nullptr) && (kappaBySegment->size() == graph.segments().size());
    const bool useKappaByCell =
        (kappaByCellIndex != nullptr) && (kappaByCellIndex->size() == static_cast<size_t>(mesh.cellIndexNum()));

    auto clampKappa = [](double k) -> double {
        if (!std::isfinite(k) || !(k > 0.0))
            return 1.0;
        // Avoid extreme coefficients blowing up diagnostics.
        constexpr double kMin = 1e-300;
        constexpr double kMax = 1e+300;
        if (k < kMin)
            return kMin;
        if (k > kMax)
            return kMax;
        return k;
    };

    auto kappaForSegment = [&](size_t sid) -> double {
        if (!useKappaBySegment)
            return 1.0;
        return clampKappa((*kappaBySegment)[sid]);
    };

    auto kappaForCellIndex = [&](int cellIndex) -> double {
        if (!useKappaByCell || cellIndex < 0 || cellIndex >= mesh.cellIndexNum())
            return 1.0;
        return clampKappa((*kappaByCellIndex)[static_cast<size_t>(cellIndex)]);
    };

    // Build per-row volumes (alpha*area) for rhs integration and compatibility correction.
    std::vector<double> volOfRow(static_cast<size_t>(N), 0.0);
    for (int i = 0; i < N; ++i)
    {
        const int cellIndex = cellIndexOfRow[i];
        double alpha = 1.0;
        if (useCutMetrics)
        {
            alpha = cutMetrics->cellMetricsByIndex()[static_cast<size_t>(cellIndex)].alpha;
            if (!std::isfinite(alpha))
                alpha = 1.0;
            alpha = std::min(1.0, std::max(0.0, alpha));
        }
        volOfRow[static_cast<size_t>(i)] = safeCellArea(mesh, cellIndex) * alpha;
    }

    // Assemble A*phi (matrix-vector) and b (rhs) without explicitly building A.
    std::vector<double> Aphi(static_cast<size_t>(N), 0.0);
    std::vector<double> b(static_cast<size_t>(N), 0.0);

    const bool rhsOk = (rhs.size() == mesh.cellIndexNum()) && (rhs.value != nullptr);
    for (int i = 0; i < N; ++i)
    {
        const int cellIndex = cellIndexOfRow[i];
        const double f = rhsOk ? rhs.value[cellIndex] : 0.0;
        b[static_cast<size_t>(i)] = f * volOfRow[static_cast<size_t>(i)];
    }

    const auto& segs = graph.segments();
    for (size_t sid = 0; sid < segs.size(); ++sid)
    {
        const auto& seg = segs[sid];
        const int Kindex = mesh.cells(seg.K).index();
        if (Kindex < 0)
            continue;
        const int i = rowOfCellIndex[static_cast<size_t>(Kindex)];
        if (i < 0)
            continue;

        double beta = 1.0;
        if (useCutMetrics)
        {
            beta = cutMetrics->betaBySegment()[sid];
            if (!std::isfinite(beta))
                beta = 1.0;
            beta = std::min(1.0, std::max(0.0, beta));
        }

        const double AfOpen = beta * seg.Af;
        if (!(AfOpen > 0.0) || !std::isfinite(AfOpen))
            continue;
        const double phiK = phi.value[Kindex];
        const double kappa_f = kappaForSegment(sid);

        if (seg.boundary == Boundary::non)
        {
            if (!(seg.df > 0.0) || !std::isfinite(seg.df))
                continue;
            const double c = kappa_f * AfOpen / seg.df;

            if (seg.L < 0)
                continue;
            const int Lindex = mesh.cells(seg.L).index();
            if (Lindex < 0)
                continue;
            const int j = rowOfCellIndex[static_cast<size_t>(Lindex)];
            if (j < 0)
                continue; // inactive neighbor is skipped in the solver

            const double phiL = phi.value[Lindex];
            Aphi[static_cast<size_t>(i)] += c * phiK - c * phiL;
        }
        else
        {
            const auto& s = bcSide(domainBC, seg.boundary);
            if (s.type == PicBCType::Dirichlet || s.type == PicBCType::Robin)
            {
                if (!(seg.df > 0.0) || !std::isfinite(seg.df))
                    continue;
                const double c = kappa_f * AfOpen / seg.df;
                // Robin on domain boundary is treated as Dirichlet-like in the solver.
                Aphi[static_cast<size_t>(i)] += c * phiK;
                b[static_cast<size_t>(i)] += c * s.value;
            }
            else if (s.type == PicBCType::Neumann)
            {
                // b += g * AfOpen
                b[static_cast<size_t>(i)] += (s.value * kappa_f) * AfOpen;
            }
            else if (s.type == PicBCType::Periodic)
            {
                const Vector2D shift = PeriodicFaceUtils::shiftToMate(mesh, seg.boundary);
                int Llabel = -1;
                if (!PeriodicFaceUtils::findMateCellLabel(mesh, seg, shift, Llabel))
                    continue;

                const int Lindex = mesh.cells(Llabel).index();
                if (Lindex < 0)
                    continue;

                const int j = rowOfCellIndex[static_cast<size_t>(Lindex)];
                if (j < 0)
                    continue;

                const double df = PeriodicFaceUtils::centerDistance(mesh, seg, Llabel, shift);
                if (!(df > 0.0) || !std::isfinite(df))
                    continue;

                const double c = kappa_f * AfOpen / df;
                const double phiL = phi.value[Lindex];
                Aphi[static_cast<size_t>(i)] += c * phiK - c * phiL;
            }
        }
    }

    // Embedded-boundary closure: add EB Robin/Dirichlet/Neumann terms.
    size_t nEBCellsConstrained = 0u;
    if (useEmbeddedBC)
    {
        for (int i = 0; i < N; ++i)
        {
            const int cellIndex = cellIndexOfRow[i];
            const auto& cm = cutMetrics->cellMetricsByIndex()[static_cast<size_t>(cellIndex)];
            if (!cm.isCut || !(cm.A_gamma > 0.0))
                continue;
            if (!(volOfRow[static_cast<size_t>(i)] > 0.0))
                continue;

            const PicEmbeddedBC* eb = nullptr;
            if (useEmbeddedBCMap && cm.gammaBoundaryId > 0)
            {
                auto it = embeddedBCById->find(cm.gammaBoundaryId);
                if (it != embeddedBCById->end() && it->second.enabled)
                    eb = &it->second;
            }
            if (!eb && useEmbeddedBCFallback)
                eb = embeddedBC;
            if (!eb)
                continue;

            const double phiK = phi.value[cellIndex];
            const double A_gamma = cm.A_gamma;
            const double kappa_g = kappaForCellIndex(cellIndex);
            const double a = eb->a;
            const double bcoef = eb->b;
            const double g = eb->g;
            const PicBCType t = eb->type;

            if (t == PicBCType::Robin)
            {
                if (std::abs(bcoef) > 0.0)
                {
                    const double addDiag = kappa_g * (a / bcoef) * A_gamma;
                    const double addRhs = kappa_g * (g / bcoef) * A_gamma;
                    if (std::isfinite(addDiag))
                        Aphi[static_cast<size_t>(i)] += addDiag * phiK;
                    if (std::isfinite(addRhs))
                        b[static_cast<size_t>(i)] += addRhs;
                    if (a != 0.0)
                        ++nEBCellsConstrained;
                }
                else
                {
                    // Dirichlet limit: a*phi=g
                    if (std::abs(a) <= 0.0)
                        continue;
                    const double phiD = g / a;
                    double d = cm.d_gamma;
                    if (!(d > 0.0) || !std::isfinite(d))
                        d = 1e-12;
                    const double c = kappa_g * (A_gamma / d);
                    if (!std::isfinite(c) || !(c > 0.0))
                        continue;
                    Aphi[static_cast<size_t>(i)] += c * phiK;
                    b[static_cast<size_t>(i)] += c * phiD;
                    ++nEBCellsConstrained;
                }
            }
            else if (t == PicBCType::Dirichlet)
            {
                const double phiD = g;
                double d = cm.d_gamma;
                if (!(d > 0.0) || !std::isfinite(d))
                    d = 1e-12;
                const double c = kappa_g * (A_gamma / d);
                if (!std::isfinite(c) || !(c > 0.0))
                    continue;
                Aphi[static_cast<size_t>(i)] += c * phiK;
                b[static_cast<size_t>(i)] += c * phiD;
                ++nEBCellsConstrained;
            }
            else if (t == PicBCType::Neumann)
            {
                const double addRhs = (g * kappa_g) * A_gamma;
                if (std::isfinite(addRhs))
                    b[static_cast<size_t>(i)] += addRhs;
            }
        }
    }

    // Pure-Neumann compatibility correction + gauge fix, mirroring the solver.
    const bool domainPureNeumann =
        (!isDirichletLike(domainBC.left.type)) &&
        (!isDirichletLike(domainBC.right.type)) &&
        (!isDirichletLike(domainBC.bottom.type)) &&
        (!isDirichletLike(domainBC.top.type));

    const bool embeddedAnchors = (useEmbeddedBC && (nEBCellsConstrained > 0u));
    const bool pureNeumann = domainPureNeumann && (!embeddedAnchors);

    if (pureNeumann)
    {
        double totalB = 0.0;
        double totalV = 0.0;
        for (int i = 0; i < N; ++i)
        {
            totalB += b[static_cast<size_t>(i)];
            totalV += volOfRow[static_cast<size_t>(i)];
        }
        if (totalV > 0.0)
        {
            const double meanBperV = totalB / totalV;
            for (int i = 0; i < N; ++i)
                b[static_cast<size_t>(i)] -= meanBperV * volOfRow[static_cast<size_t>(i)];
        }

        // Gauge: pin one DOF (phi_pin=0) by adding +1 to the diagonal.
        int pinRow = 0;
        for (int i = 0; i < N; ++i)
        {
            if (volOfRow[static_cast<size_t>(i)] > 0.0)
            {
                pinRow = i;
                break;
            }
        }
        const int pinCellIndex = cellIndexOfRow[pinRow];
        const double phiPin = phi.value[pinCellIndex];
        Aphi[static_cast<size_t>(pinRow)] += 1.0 * phiPin;
        b[static_cast<size_t>(pinRow)] = 0.0;
    }

    // Final residuals and norms.
    double sum2 = 0.0;
    double bsum2 = 0.0;
    double blinf = 0.0;
    double linf = 0.0;
    double scf = 0.0;
    double seb = 0.0;
    std::vector<double> cfAbs;
    std::vector<double> ebAbs;
    cfAbs.reserve(static_cast<size_t>(N));
    ebAbs.reserve(static_cast<size_t>(N));
    for (int i = 0; i < N; ++i)
    {
        const double r = Aphi[static_cast<size_t>(i)] - b[static_cast<size_t>(i)];
        const double bi = b[static_cast<size_t>(i)];
        const double ar = std::abs(r);
        const double abi = std::abs(bi);
        sum2 += r * r;
        bsum2 += bi * bi;
        blinf = std::max(blinf, abi);
        linf = std::max(linf, ar);
        if (rGout)
        {
            const int cellIndex = cellIndexOfRow[i];
            rGout->value[cellIndex] = r;
        }

        const int cellIndex = cellIndexOfRow[i];
        if (cellIndex >= 0 && cellIndex < mesh.cellIndexNum())
        {
            const int cLbl = mesh.cellIndexToLabel[cellIndex];
            const Cell& c = mesh.cells(cLbl);
            if (isCompositeBandCell(mesh, c))
            {
                scf = std::max(scf, ar);
                cfAbs.push_back(ar);
            }

            if (useCutMetrics)
            {
                const auto& cmv = cutMetrics->cellMetricsByIndex();
                if (static_cast<size_t>(cellIndex) < cmv.size())
                {
                    const auto& cm = cmv[static_cast<size_t>(cellIndex)];
                    if (cm.isCut)
                    {
                        seb = std::max(seb, ar);
                        ebAbs.push_back(ar);
                    }
                }
            }
        }
    }
    out.rG_L2 = std::sqrt(sum2);
    out.rG_Linf = linf;
    out.rG_L2_per_active = (N > 0) ? std::sqrt(sum2 / static_cast<double>(N)) : 0.0;
    out.rG_L2_rel_rhs = (bsum2 > 0.0) ? (out.rG_L2 / std::sqrt(bsum2)) : 0.0;
    out.rG_Linf_rel_rhs = (blinf > 0.0) ? (out.rG_Linf / blinf) : 0.0;
    out.S_cf = scf;
    out.S_eb = seb;
    out.S_eb_rel_rhs = (blinf > 0.0) ? (out.S_eb / blinf) : 0.0;

    auto quantileInPlace = [](std::vector<double>& v, double q) -> double {
        if (v.empty())
            return 0.0;
        if (!std::isfinite(q))
            q = 0.0;
        q = std::min(1.0, std::max(0.0, q));
        const std::size_t n = v.size();
        const std::size_t k = static_cast<std::size_t>(std::llround(q * static_cast<double>(n - 1)));
        auto it = v.begin() + static_cast<std::vector<double>::difference_type>(k);
        std::nth_element(v.begin(), it, v.end());
        return *it;
    };

    {
        std::vector<double> tmp = cfAbs;
        out.S_cf_p50 = quantileInPlace(tmp, 0.50);
        tmp = cfAbs;
        out.S_cf_p90 = quantileInPlace(tmp, 0.90);
        tmp = cfAbs;
        out.S_cf_p99 = quantileInPlace(tmp, 0.99);

        tmp = ebAbs;
        out.S_eb_p50 = quantileInPlace(tmp, 0.50);
        tmp = ebAbs;
        out.S_eb_p90 = quantileInPlace(tmp, 0.90);
        tmp = ebAbs;
        out.S_eb_p99 = quantileInPlace(tmp, 0.99);
    }

    // Energy and charge totals.
    double W = 0.0;
    double Qgrid = 0.0;
    if (Q && Q->value && Q->size() == mesh.cellIndexNum())
    {
        for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
        {
            const double Qi = Q->value[ci];
            Qgrid += Qi;
            W += 0.5 * Qi * phi.value[ci];
        }
    }
    else
    {
        // Fall back: integrate 0.5 * rho * phi * area, with rho = rhs * epsilon0.
        // NOTE: rhs here is the Poisson RHS (rho/epsilon0), so rho = rhs*epsilon0.
        constexpr double epsilon0 = 8.854187817e-12;
        for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
        {
            const int cLbl = mesh.cellIndexToLabel[ci];
            const double area = mesh.cells(cLbl).area();
            const double rhoi = rhs.value[ci] * epsilon0;
            Qgrid += rhoi * area;
            W += 0.5 * rhoi * phi.value[ci] * area;
        }
    }
    out.energyW = W;
    out.Qgrid = Qgrid;
    out.dQ = Qgrid - Qmacro;
    out.dQrel = (std::abs(Qmacro) > 0.0) ? (out.dQ / Qmacro) : 0.0;

    return out;
}


} // namespace Diagnostics
} // namespace ExplicitPIC
