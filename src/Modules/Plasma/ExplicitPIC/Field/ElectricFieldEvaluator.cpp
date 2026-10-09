#include "Modules/Plasma/ExplicitPIC/Field/ElectricFieldEvaluator.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>

#include "Core/IO/ConsoleControl.h"
#include "Core/Mesh/PeriodicFaceUtils.h"

namespace ExplicitPIC
{
namespace Field
{

namespace
{
static inline double safeDiv(const double a, const double b, const double eps)
{
    return a / (std::abs(b) > eps ? b : (b >= 0.0 ? eps : -eps));
}

static inline const PicBoundarySide& pickDomainSide(
    const PicDomainBC& bc,
    const FaceSide side)
{
    // FaceSide ordering follows mesh convention: west/south/east/north.
    switch (side)
    {
    case FaceSide::West:  return bc.left;
    case FaceSide::East:  return bc.right;
    case FaceSide::South: return bc.bottom;
    case FaceSide::North: return bc.top;
    default:              return bc.left;
    }
}

static inline bool isFinite2(const Vector2D& v)
{
    return std::isfinite(v.x) && std::isfinite(v.y);
}

static inline bool mapHasEnabledBC(const PicEmbeddedBCMap* map)
{
    if (!map)
        return false;
    for (const auto& kv : *map)
    {
        if (kv.second.enabled)
            return true;
    }
    return false;
}

static inline double stabilizedEmbeddedDistance(const double dGamma, const double cellArea, const double eps)
{
    const double h = (std::isfinite(cellArea) && cellArea > 0.0) ? std::sqrt(cellArea) : 1.0;
    const double floor = std::max(std::max(1e-12, eps), 2.0e-2 * h);
    if (!std::isfinite(dGamma))
        return floor;
    return std::max(dGamma, floor);
}

} // namespace

std::vector<double> ElectricFieldEvaluator::computeFaceNormalField(
    const QuadtreeMesh& mesh,
    const FaceSegmentGraph& faceGraph,
    const PicDomainBC& domainBC,
    const CellField<double>& phi,
    const Options& opt)
{
    const auto& segs = faceGraph.segments();
    std::vector<double> En(segs.size(), 0.0);

    const double eps = opt.eps;

    for (size_t s = 0; s < segs.size(); ++s)
    {
        const FaceSegment& seg = segs[s];
        const int kLabel = seg.K;
        if (kLabel < 0) continue;

        const int kIndex = mesh.cells(kLabel).index();
        if (kIndex < 0 || kIndex >= phi.size()) continue;
        const double phiK = phi.value[kIndex];

        const double df = std::max(seg.df, eps);

        if (seg.boundary != Boundary::non || seg.L < 0)
        {
            const PicBoundarySide& bc = pickDomainSide(domainBC, seg.side);
            if (bc.type == PicBCType::Periodic)
            {
                const Vector2D shift = PeriodicFaceUtils::shiftToMate(mesh, seg.side);
                const int lLabel = PeriodicFaceUtils::findMateCellLabel(mesh, seg, shift);
                if (lLabel >= 0)
                {
                    const int lIndex = mesh.cells(lLabel).index();
                    if (lIndex >= 0 && lIndex < phi.size())
                    {
                        const double phiL = phi.value[lIndex];
                        const double df_eff =
                            std::max(PeriodicFaceUtils::centerDistance(mesh, seg, lLabel, shift), eps);
                        En[s] = -(phiL - phiK) / df_eff;
                    }
                }
            }
            else if (bc.type == PicBCType::Neumann)
            {
                // bc.value is dphi/dn (outward). E_n = -dphi/dn.
                En[s] = -bc.value;
            }
            else
            {
                // Dirichlet or (domain) Robin treated as Dirichlet: bc.value is phi_b.
                const double phiB = bc.value;
                En[s] = -(phiB - phiK) / df;
            }
        }
        else
        {
            const int lLabel = seg.L;
            const int lIndex = mesh.cells(lLabel).index();
            if (lIndex < 0 || lIndex >= phi.size()) continue;
            const double phiL = phi.value[lIndex];
            En[s] = -(phiL - phiK) / df;
        }

        if (!std::isfinite(En[s])) En[s] = 0.0;
    }

    return En;
}

std::vector<double> ElectricFieldEvaluator::computeEmbeddedNormalField(
    const QuadtreeMesh& mesh,
    const Solvers::Poisson::CutCellMetrics* metrics,
    const PicEmbeddedBC* embeddedBC,
    const PicEmbeddedBCMap* embeddedBCById,
    const CellField<double>& phi,
    const Options& opt)
{
    std::vector<double> EnGamma(mesh.cellIndexNum(), 0.0);

    if (!metrics) return EnGamma;

    const double eps = opt.eps;
    const auto& cm = metrics->cellMetricsByIndex();
    if ((int)cm.size() != mesh.cellIndexNum()) return EnGamma;

    const bool useMap = mapHasEnabledBC(embeddedBCById);
    const bool useFallback = (embeddedBC != nullptr) && embeddedBC->enabled;
    if (!useMap && !useFallback)
        return EnGamma;

    size_t nApplied = 0u;
    size_t nMissing = 0u;
    size_t nFallback = 0u;
    std::map<int, size_t> byId;

    for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
    {
        if (!metrics->isActiveCellIndex(ci))
            continue;

        const auto& m = cm[ci];
        if (!(m.A_gamma > 0.0) || !(m.d_gamma > eps))
            continue;

        const PicEmbeddedBC* eb = nullptr;
        if (useMap && m.gammaBoundaryId > 0)
        {
            auto it = embeddedBCById->find(m.gammaBoundaryId);
            if (it != embeddedBCById->end() && it->second.enabled)
            {
                eb = &it->second;
            }
        }
        if (!eb && useFallback)
        {
            eb = embeddedBC;
            ++nFallback;
        }
        if (!eb)
        {
            ++nMissing;
            if (opt.embeddedBCMissingPolicy == Options::EmbeddedBCMissingPolicy::Error)
            {
                throw std::runtime_error("Missing embedded BC for boundaryId " +
                                         std::to_string(m.gammaBoundaryId) +
                                         " (no map entry and no fallback).");
            }
            continue;
        }

        const double phiK = (ci < phi.size() ? phi.value[ci] : 0.0);
        const int cLbl = mesh.cellIndexToLabel[ci];
        const double d = stabilizedEmbeddedDistance(m.d_gamma, mesh.cells(cLbl).area(), eps);

        double En = 0.0;
        if (eb->type == PicBCType::Dirichlet)
        {
            const double phiG = eb->g;
            En = -(phiG - phiK) / d;
        }
        else if (eb->type == PicBCType::Neumann)
        {
            // g is dphi/dn.
            En = -eb->g;
        }
        else
        {
            // Robin: a*phi + b*dphi/dn = g
            const double denom = eb->a * d + eb->b;
            const double dphi_dn = safeDiv(eb->g - eb->a * phiK, denom, eps);
            En = -dphi_dn;
        }

        if (!std::isfinite(En)) En = 0.0;
        EnGamma[ci] = En;
        ++nApplied;
        if (m.gammaBoundaryId > 0)
            byId[m.gammaBoundaryId] += 1;
    }

    static long long s_evalCalls = 0;
    const long long callIndex = s_evalCalls++;
    if (ConsoleControl::enabled(ConsoleControl::Level::kDebug) &&
        ConsoleControl::shouldPrintPeriodic(static_cast<int>(callIndex), ConsoleControl::g_config.every))
    {
        std::cout << "[ElectricFieldEvaluator] EB En_gamma: applied=" << nApplied
                  << ", missing=" << nMissing;
        if (nFallback > 0)
            std::cout << ", fallback=" << nFallback;
        if (!byId.empty())
        {
            std::cout << ", byBoundaryId={";
            bool first = true;
            for (const auto& kv : byId)
            {
                if (!first) std::cout << ", ";
                first = false;
                std::cout << kv.first << ":" << kv.second;
            }
            std::cout << "}";
        }
        std::cout << "\n";
    }

    return EnGamma;
}

ElectricFieldEvaluator::NormalFields ElectricFieldEvaluator::evaluate(
    const QuadtreeMesh& mesh,
    const FaceSegmentGraph& faceGraph,
    const Solvers::Poisson::CutCellMetrics* metrics,
    const PicDomainBC& domainBC,
    const PicEmbeddedBC* embeddedBC,
    const PicEmbeddedBCMap* embeddedBCById,
    const CellField<double>& phi,
    const Options& opt)
{
    NormalFields out;
    out.EnBySegment = computeFaceNormalField(mesh, faceGraph, domainBC, phi, opt);
    out.EnGammaByCellIndex = computeEmbeddedNormalField(mesh, metrics, embeddedBC, embeddedBCById, phi, opt);
    return out;
}

void ElectricFieldEvaluator::reconstructCellCenteredField(
    const QuadtreeMesh& mesh,
    const FaceSegmentGraph& faceGraph,
    const Solvers::Poisson::CutCellMetrics* metrics,
    const std::vector<double>& EnBySegment,
    const std::vector<double>* EnGammaByCellIndex,
    const PicEmbeddedBC* embeddedBC,
    const PicEmbeddedBCMap* embeddedBCById,
    CellField<Vector2D>& Ecell,
    const Options& opt)
{
    Ecell.initialize(mesh, Vector2D(0.0, 0.0));

    // Weighted least-squares reconstruction:
    //   minimize sum_i omega_i (n_i · E - En_i)^2
    // which yields the normal equations:
    //   A = sum_i omega_i (n_i n_i^T),  b = sum_i omega_i (n_i En_i),  A E = b.
    //
    // For axis-aligned faces only, A is diagonal and the reconstruction reduces
    // to the familiar split by x/y faces. We optionally include the EB normal
    // constraint (when provided), which introduces off-diagonal coupling.
    //
    // Weight model: omega = beta_f * A_f for face segments, and omega = A_gamma
    // for the embedded segment.
    //
    // beta_f priority:
    //   (1) cutMetrics->betaBySegment(sid)
    //   (2) fallback (when metrics are absent):
    //         - if the underlying edge is not cut by EB: beta = 1
    //         - if the edge is cut: beta = edgeCutLength / edgeLength (clamped)

    const auto& segs = faceGraph.segments();
    const double eps = opt.eps;

    const std::vector<double>* betaBySeg = nullptr;
    const std::vector<Solvers::Poisson::CutCellMetrics::CellMetric>* cm = nullptr;
    if (metrics)
    {
        betaBySeg = &metrics->betaBySegment();
        cm = &metrics->cellMetricsByIndex();
    }

    const bool useEmbedded = opt.includeEmbedded && (EnGammaByCellIndex != nullptr) && metrics &&
        (mapHasEnabledBC(embeddedBCById) || (embeddedBC != nullptr && embeddedBC->enabled));

    for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
    {
        if (metrics && !metrics->isActiveCellIndex(ci))
        {
            Ecell.value[ci] = Vector2D(0.0, 0.0);
            continue;
        }
        if (cm && (*cm)[ci].alpha <= 0.0)
        {
            Ecell.value[ci] = Vector2D(0.0, 0.0);
            continue;
        }

        const auto range = faceGraph.segmentRangeForCellIndex(ci);

        double A11 = 0.0;
        double A12 = 0.0;
        double A22 = 0.0;
        double b1 = 0.0;
        double b2 = 0.0;

        for (size_t s = range.first; s < range.second; ++s)
        {
            if (s >= segs.size() || s >= EnBySegment.size()) break;
            const FaceSegment& seg = segs[s];

            // beta priority: metrics->betaBySegment else edgeCutLength/edgeLength fallback.
            double betaS = 1.0;
            if (betaBySeg && s < betaBySeg->size())
            {
                betaS = (*betaBySeg)[s];
            }
            else
            {
                // Fallback when CutCellMetrics are not available.
                // IMPORTANT: for a non-cut edge, Edge::cutLength() is typically 0.
                // Using cutLength/length naively would incorrectly set beta=0 on regular
                // edges and zero out the reconstructed field (Ex/Ey become identically 0).
                if (!seg.edgeCut)
                {
                    betaS = 1.0;
                }
                else if (seg.edgeLength > eps)
                {
                    betaS = seg.edgeCutLength / seg.edgeLength;
                    betaS = std::max(0.0, std::min(1.0, betaS));
                }
                else
                {
                    // Degenerate edge length; keep a safe default.
                    betaS = (seg.edgeCutLength > 0.0 ? 1.0 : 0.0);
                }
            }

            if (!(betaS > 0.0)) continue;

            const double omega = betaS * seg.Af;
            if (!(omega > 0.0)) continue;

            const double En = EnBySegment[s];
            const Vector2D& n = seg.n;
            A11 += omega * n.x * n.x;
            A12 += omega * n.x * n.y;
            A22 += omega * n.y * n.y;
            b1 += omega * n.x * En;
            b2 += omega * n.y * En;
        }

        // Embedded-boundary constraint (one segment per cut cell).
        if (useEmbedded && cm && (ci < static_cast<int>(cm->size())) &&
            (ci < static_cast<int>(EnGammaByCellIndex->size())))
        {
            const auto& m = (*cm)[ci];
            if (m.A_gamma > 0.0 && std::isfinite(m.n_gamma.x) && std::isfinite(m.n_gamma.y))
            {
                const double omega = m.A_gamma;
                if (omega > 0.0)
                {
                    const double En = (*EnGammaByCellIndex)[ci];
                    const Vector2D& n = m.n_gamma;
                    A11 += omega * n.x * n.x;
                    A12 += omega * n.x * n.y;
                    A22 += omega * n.y * n.y;
                    b1 += omega * n.x * En;
                    b2 += omega * n.y * En;
                }
            }
        }

        Vector2D E(0.0, 0.0);
        const double offDiag = std::abs(A12);
        const double diagScale = std::max(1.0, A11 + A22);
        if (opt.axisAlignedSimplification && offDiag <= eps * diagScale)
        {
            if (A11 > eps) E.x = b1 / A11;
            if (A22 > eps) E.y = b2 / A22;
        }
        else
        {
            const double det = A11 * A22 - A12 * A12;
            if (std::abs(det) > eps * diagScale)
            {
                E.x = (b1 * A22 - b2 * A12) / det;
                E.y = (-b1 * A12 + b2 * A11) / det;
            }
            else
            {
                if (A11 > eps) E.x = b1 / A11;
                if (A22 > eps) E.y = b2 / A22;
            }
        }

        // requirement (3): no corresponding faces/constraints -> component remains 0.
        if (!isFinite2(E)) E = Vector2D(0.0, 0.0);
        Ecell.value[ci] = E;
    }
}

void ElectricFieldEvaluator::evaluateCellCenteredField(
    const QuadtreeMesh& mesh,
    const FaceSegmentGraph& faceGraph,
    const Solvers::Poisson::CutCellMetrics* metrics,
    const PicDomainBC& domainBC,
    const PicEmbeddedBC* embeddedBC,
    const PicEmbeddedBCMap* embeddedBCById,
    const CellField<double>& phi,
    CellField<Vector2D>& Ecell,
    std::vector<double>* EnBySegmentOut,
    std::vector<double>* EnGammaByCellIndexOut,
    const Options& opt)
{
    std::vector<double> EnBySegment = computeFaceNormalField(mesh, faceGraph, domainBC, phi, opt);

    std::vector<double> EnGamma;
    if (mapHasEnabledBC(embeddedBCById) || (embeddedBC != nullptr && embeddedBC->enabled))
        EnGamma = computeEmbeddedNormalField(mesh, metrics, embeddedBC, embeddedBCById, phi, opt);

    const std::vector<double>* EnGammaPtr = (!EnGamma.empty()) ? &EnGamma : nullptr;

    reconstructCellCenteredField(mesh, faceGraph, metrics, EnBySegment, EnGammaPtr,
                                 embeddedBC, embeddedBCById, Ecell, opt);

    if (EnBySegmentOut) *EnBySegmentOut = std::move(EnBySegment);
    if (EnGammaByCellIndexOut) *EnGammaByCellIndexOut = std::move(EnGamma);
}

} // namespace Field
} // namespace ExplicitPIC
