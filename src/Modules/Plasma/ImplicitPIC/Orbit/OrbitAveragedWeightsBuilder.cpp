#include "Modules/Plasma/ImplicitPIC/Orbit/OrbitAveragedWeightsBuilder.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <exception>
#include <limits>
#include <stdexcept>
#include <utility>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "Core/Mesh/BackgroundMesh.h"
#include "Modules/Plasma/ExplicitPIC/Deposit/CouplingWeightsCache.h"

#include "Modules/Plasma/ImplicitPIC/Orbit/PredictorOrbitSegmenter.h"

namespace ImplicitPIC
{
namespace Orbit
{
namespace
{
#if defined(LAMP_ORBITAVG_FINE_TIMING)
constexpr bool kOrbitAvgFineTimingEnabled = (LAMP_ORBITAVG_FINE_TIMING != 0);
#else
constexpr bool kOrbitAvgFineTimingEnabled = false;
#endif

inline Vector2D lerp(const Vector2D& a, const Vector2D& b, double s)
{
    return a * (1.0 - s) + b * s;
}

inline double wrap01(double a)
{
    // Wrap into [0,1).
    a = std::fmod(a, 1.0);
    if (a < 0.0)
        a += 1.0;
    // Avoid landing exactly at 1.0 due to roundoff.
    if (a >= 1.0)
        a = std::nextafter(1.0, 0.0);
    return a;
}

inline Vector2D wrapToDomain(const BackgroundMesh& mesh, const Vector2D& x, bool px, bool py, double eps)
{
    if (!(px || py))
        return x;

    const double xmin = mesh.minCoord.x;
    const double xmax = mesh.maxCoord.x;
    const double ymin = mesh.minCoord.y;
    const double ymax = mesh.maxCoord.y;

    const bool xInside = (!px) || ((x.x >= xmin) && (x.x < xmax));
    const bool yInside = (!py) || ((x.y >= ymin) && (x.y < ymax));
    if (xInside && yInside)
        return x;

    Vector2D out = x;

    if (px)
    {
        const double L = xmax - xmin;
        if (L > 0.0)
        {
            const double s = (out.x - xmin) / L;
            out.x = xmin + wrap01(s) * L;
            // Keep away from maxCoord to satisfy half-open convention.
            if (out.x >= xmax)
                out.x = xmax - std::abs(eps);
        }
    }
    if (py)
    {
        const double L = ymax - ymin;
        if (L > 0.0)
        {
            const double s = (out.y - ymin) / L;
            out.y = ymin + wrap01(s) * L;
            if (out.y >= ymax)
                out.y = ymax - std::abs(eps);
        }
    }

    return out;
}

inline Vector2D unwrapToNearestImage(const BackgroundMesh& mesh, const Vector2D& xRef, const Vector2D& x, bool px, bool py)
{
    if (!(px || py))
        return x;

    Vector2D out = x;
    if (px)
    {
        const double xmin = mesh.minCoord.x;
        const double xmax = mesh.maxCoord.x;
        const double L = xmax - xmin;
        if (L > 0.0)
        {
            const double dx = out.x - xRef.x;
            out.x = out.x - std::round(dx / L) * L;
        }
    }
    if (py)
    {
        const double ymin = mesh.minCoord.y;
        const double ymax = mesh.maxCoord.y;
        const double L = ymax - ymin;
        if (L > 0.0)
        {
            const double dy = out.y - xRef.y;
            out.y = out.y - std::round(dy / L) * L;
        }
    }

    return out;
}

inline void reduceTripsToCSR(std::vector<std::pair<int, double>>& trips,
                             std::vector<int>& outCell,
                             std::vector<double>& outW)
{
    outCell.clear();
    outW.clear();
    if (trips.empty())
        return;

    std::sort(trips.begin(), trips.end(), [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
        return a.first < b.first;
    });

    outCell.reserve(trips.size());
    outW.reserve(trips.size());

    int currCell = trips.front().first;
    double currW = trips.front().second;
    for (std::size_t k = 1; k < trips.size(); ++k)
    {
        const int cell = trips[k].first;
        const double wv = trips[k].second;
        if (cell == currCell)
        {
            currW += wv;
        }
        else
        {
            if (currW != 0.0)
            {
                outCell.push_back(currCell);
                outW.push_back(currW);
            }
            currCell = cell;
            currW = wv;
        }
    }
    if (currW != 0.0)
    {
        outCell.push_back(currCell);
        outW.push_back(currW);
    }
}

inline double l1DiffCSR(const std::vector<int>& aCell,
                        const std::vector<double>& aW,
                        const std::vector<int>& bCell,
                        const std::vector<double>& bW)
{
    // Both inputs are expected sorted by cell.
    double s = 0.0;
    std::size_t i = 0u;
    std::size_t j = 0u;
    while (i < aCell.size() || j < bCell.size())
    {
        if (j >= bCell.size() || (i < aCell.size() && aCell[i] < bCell[j]))
        {
            s += std::abs(aW[i]);
            ++i;
        }
        else if (i >= aCell.size() || (j < bCell.size() && bCell[j] < aCell[i]))
        {
            s += std::abs(bW[j]);
            ++j;
        }
        else
        {
            s += std::abs(aW[i] - bW[j]);
            ++i;
            ++j;
        }
    }
    return s;
}

inline void reduceTripsToCSRStamp(const std::vector<std::pair<int, double>>& trips,
                                  std::vector<int>& outCell,
                                  std::vector<double>& outW,
                                  std::vector<double>& accumByCell,
                                  std::vector<std::uint32_t>& stampByCell,
                                  std::vector<int>& touchedCells,
                                  std::uint32_t& stampToken,
                                  bool deterministicOrdering)
{
    outCell.clear();
    outW.clear();
    if (trips.empty())
        return;

    if (accumByCell.empty() || stampByCell.size() != accumByCell.size())
    {
        auto tmp = trips;
        reduceTripsToCSR(tmp, outCell, outW);
        return;
    }

    ++stampToken;
    if (stampToken == 0u)
    {
        std::fill(stampByCell.begin(), stampByCell.end(), 0u);
        stampToken = 1u;
    }

    touchedCells.clear();
    touchedCells.reserve(std::min<std::size_t>(trips.size(), accumByCell.size()));

    const int nCell = static_cast<int>(accumByCell.size());
    for (const auto& tw : trips)
    {
        const int c = tw.first;
        if (c < 0 || c >= nCell)
            continue;
        const double wv = tw.second;
        if (wv == 0.0)
            continue;

        if (stampByCell[static_cast<std::size_t>(c)] != stampToken)
        {
            stampByCell[static_cast<std::size_t>(c)] = stampToken;
            accumByCell[static_cast<std::size_t>(c)] = wv;
            touchedCells.push_back(c);
        }
        else
        {
            accumByCell[static_cast<std::size_t>(c)] += wv;
        }
    }

    if (touchedCells.empty())
        return;

    if (deterministicOrdering)
        std::sort(touchedCells.begin(), touchedCells.end());
    outCell.reserve(touchedCells.size());
    outW.reserve(touchedCells.size());

    for (int c : touchedCells)
    {
        const double wv = accumByCell[static_cast<std::size_t>(c)];
        if (wv != 0.0)
        {
            outCell.push_back(c);
            outW.push_back(wv);
        }
    }
}
} // namespace

bool buildOrbitAveragedCacheAdaptiveGauss24(const BackgroundMesh& mesh,
                                           const std::vector<double>& hp,
                                           const CouplingWeightsCache::BuildOptions& buildOpt,
                                           const std::size_t numParticles,
                                           const std::vector<Vector2D>& x0,
                                           const std::vector<Vector2D>& x1,
                                           const std::vector<int>& hostCellLabelByParticle,
                                           const double dt,
                                           CouplingWeightsCache& out,
                                           OrbitAveragedWeightsBuilderWorkspace* workspace,
                                           const OrbitAveragedWeightsBuilderOptions& opt,
                                           OrbitAveragedWeightsBuilderTiming* timing)
{
    using Clock = std::chrono::steady_clock;
    auto secondsBetween = [](const Clock::time_point& a, const Clock::time_point& b) -> double {
        return std::chrono::duration_cast<std::chrono::duration<double>>(b - a).count();
    };
    const auto tTotal0 = Clock::now();
    if (timing)
        *timing = OrbitAveragedWeightsBuilderTiming{};

    OrbitAveragedWeightsBuilderWorkspace localWorkspace;
    OrbitAveragedWeightsBuilderWorkspace& ws = workspace ? *workspace : localWorkspace;

    if (dt <= 0.0)
        return false;
    if (hp.empty())
        return false;
    if (x0.size() != numParticles || x1.size() != numParticles)
        return false;
    if (hostCellLabelByParticle.size() != numParticles)
        return false;

    // Only gauss2_4 currently.
    if (!opt.adaptiveRule.empty() && opt.adaptiveRule != "gauss2_4")
        return false;

    if (timing)
        timing->fusedMode = 2; // adaptive fused

    const bool hpIsConst = (hp.size() == 1u);
    const double hpConst = hpIsConst ? hp.front() : 0.0;
    auto hpForParticle = [&](std::size_t p) -> double {
        if (hpIsConst)
            return hpConst;
        if (p < hp.size())
            return hp[p];
        return hp.back();
    };

    auto& particleNnz = ws.particleNnz;
    particleNnz.assign(numParticles, 0);
    auto& outOffsets = ws.outOffsets;
    outOffsets.assign(numParticles + 1u, 0);
    std::vector<int> outCell;
    std::vector<double> outW;
    std::size_t outNnz = 0u;

    // Gauss-Legendre nodes/weights.
    constexpr double gl2_x[2] = {-0.5773502691896257645, 0.5773502691896257645};
    constexpr double gl2_w[2] = {1.0, 1.0};
    constexpr double gl4_x[4] = {-0.8611363115940525752,
                                 -0.3399810435848562648,
                                 0.3399810435848562648,
                                 0.8611363115940525752};
    constexpr double gl4_w[4] = {0.3478548451374538574,
                                 0.6521451548625461426,
                                 0.6521451548625461426,
                                 0.3478548451374538574};

    const bool periodicX = buildOpt.periodicX;
    const bool periodicY = buildOpt.periodicY;
    constexpr double epsWrap = 1e-12;

    auto hostLabelForSample = [&](const Vector2D& x, int fallbackLabel) -> int {
        int lbl = mesh.findCellLabel(x);
        if (lbl < 0 && (periodicX || periodicY))
        {
            const Vector2D xw = wrapToDomain(mesh, x, periodicX, periodicY, epsWrap);
            lbl = mesh.findCellLabel(xw);
        }
        if (lbl >= 0)
            return lbl;
        if (fallbackLabel >= 0)
            return fallbackLabel;
        return 0;
    };

    auto buildParticleAdaptiveFused = [&](std::size_t p,
                                          CouplingWeightsCache& evalCache,
                                          std::vector<Vector2D>& evalPos,
                                          std::vector<int>& evalHostLabel,
                                          std::vector<double>& hpSingle,
                                          CouplingWeightsCache::BuildOptions& evalOpt,
                                          std::vector<std::pair<int, double>>& scratchTrips,
                                          std::vector<std::pair<int, double>>& trips2,
                                          std::vector<std::pair<int, double>>& trips4,
                                          std::vector<int>& tmpCell2,
                                          std::vector<double>& tmpW2,
                                          std::vector<int>& tmpCell4,
                                          std::vector<double>& tmpW4,
                                          std::vector<int>& outCellLocal,
                                          std::vector<double>& outWLocal,
                                          std::vector<double>& accumByCell,
                                          std::vector<std::uint32_t>& stampByCell,
                                          std::vector<int>& touchedCells,
                                          std::uint32_t& stampToken,
                                          std::size_t& sampleCountOut,
                                          std::size_t& sampleFallbackOut,
                                          std::size_t& sampleNnzOut,
                                          std::size_t& acceptCountOut,
                                          std::size_t& refineCountOut,
                                          std::size_t& hitMaxDepthCountOut,
                                          std::size_t& acceptedDepthSumOut,
                                          double& evalTimeOut,
                                          double& accumTimeOut)
    {
        const auto tParticle0 = kOrbitAvgFineTimingEnabled ? Clock::now() : Clock::time_point{};
        double evalTimeLocal = 0.0;
        std::size_t localSampleCount = 0u;
        std::size_t localSampleFallback = 0u;
        std::size_t localSampleNnz = 0u;

        // Batch accepted Gauss2 samples to reduce tiny evalCache.build calls.
        std::vector<double> pendingTVals;
        std::vector<double> pendingCoeffVals;
        std::vector<int> pendingHostVals;
        pendingTVals.reserve(128);
        pendingCoeffVals.reserve(128);
        pendingHostVals.reserve(128);

        scratchTrips.clear();
        outCellLocal.clear();
        outWLocal.clear();

        evalHostLabel[0] = hostCellLabelByParticle[p];
        if (!hpIsConst)
            hpSingle[0] = hpForParticle(p);

        const Vector2D x0p = x0[p];
        const Vector2D x1p = unwrapToNearestImage(mesh, x0p, x1[p], periodicX, periodicY);

        auto posAtTime = [&](double t) -> Vector2D {
            const double s = (dt > 0.0) ? (t / dt) : 0.0;
            const Vector2D xr = lerp(x0p, x1p, std::clamp(s, 0.0, 1.0));
            return wrapToDomain(mesh, xr, periodicX, periodicY, epsWrap);
        };

        auto probeHostLabelAt = [&](double t) -> int {
            const Vector2D xw = posAtTime(t);
            int lbl = mesh.findCellLabel(xw);
            if (lbl >= 0)
                return lbl;
            const int fallbackLabel = hostCellLabelByParticle[p];
            if (fallbackLabel >= 0)
                return fallbackLabel;
            return 0;
        };

        auto isLabelCutCell = [&](int lbl) -> bool {
            if (lbl < 0)
                return false;
            const auto& c = mesh.cells(lbl);
            if (buildOpt.levelSetI >= 0 && static_cast<std::size_t>(buildOpt.levelSetI) < c.isCut.size())
                return c.isCut[static_cast<std::size_t>(buildOpt.levelSetI)] != 0u;
            if (buildOpt.cutByAlpha && buildOpt.alphaByCellIndex)
            {
                const int ci = c.index();
                if (ci >= 0 && static_cast<std::size_t>(ci) < buildOpt.alphaByCellIndex->size())
                {
                    const double a = (*buildOpt.alphaByCellIndex)[static_cast<std::size_t>(ci)];
                    return (a > buildOpt.cutAlphaEps) && (a < (1.0 - buildOpt.cutAlphaEps));
                }
            }
            return false;
        };

        const double sheathPad = std::max(0.0, opt.adaptiveCrossSplitSheathPad);
        auto nearDomainWall = [&](const Vector2D& xw) -> bool {
            if (!(sheathPad > 0.0))
                return false;
            const double dxMin = std::abs(xw.x - mesh.minCoord.x);
            const double dxMax = std::abs(mesh.maxCoord.x - xw.x);
            const double dyMin = std::abs(xw.y - mesh.minCoord.y);
            const double dyMax = std::abs(mesh.maxCoord.y - xw.y);
            return (dxMin <= sheathPad) || (dxMax <= sheathPad) || (dyMin <= sheathPad) || (dyMax <= sheathPad);
        };

        auto evalBatch = [&](const double* tVals,
                             const double* coeffVals,
                             const int* knownHostLabels,
                             int n,
                             std::vector<std::pair<int, double>>& dest) {
            if (n <= 0)
                return;

            evalPos.resize(static_cast<std::size_t>(n));
            evalHostLabel.resize(static_cast<std::size_t>(n));
            if (!hpIsConst)
                hpSingle.resize(static_cast<std::size_t>(n));

            for (int i = 0; i < n; ++i)
            {
                evalPos[static_cast<std::size_t>(i)] = posAtTime(tVals[i]);
                const int knownHost = knownHostLabels ? knownHostLabels[i] : -1;
                evalHostLabel[static_cast<std::size_t>(i)] =
                    (knownHost >= 0) ? knownHost : hostLabelForSample(evalPos[static_cast<std::size_t>(i)], hostCellLabelByParticle[p]);
                if (!hpIsConst)
                    hpSingle[static_cast<std::size_t>(i)] = hpForParticle(p);
            }

            const auto tEval0 = kOrbitAvgFineTimingEnabled ? Clock::now() : Clock::time_point{};
            const bool prevParallelBuild = evalOpt.parallelBuild;
            // Fused orbit-averaged path uses thread-local scratch and reduction-heavy flows;
            // keep inner cache build serial here to avoid unstable nested/combined parallel behavior.
            evalOpt.parallelBuild = false;
            if (hpIsConst)
                evalCache.build(mesh, evalPos, hp, evalOpt);
            else
                evalCache.build(mesh, evalPos, hpSingle, evalOpt);
            evalOpt.parallelBuild = prevParallelBuild;
            if (kOrbitAvgFineTimingEnabled)
                evalTimeLocal += secondsBetween(tEval0, Clock::now());

            localSampleCount += static_cast<std::size_t>(n);
            localSampleFallback += evalCache.lastBuildFallbackCount();
            localSampleNnz += evalCache.nnz();

            const auto& offEval = evalCache.offsets();
            const auto& ciEval = evalCache.cellIndices();
            const auto& wEval = evalCache.weights();
            if (ciEval.size() != wEval.size() || offEval.size() != static_cast<std::size_t>(n + 1))
                throw std::runtime_error("[OrbitAveragedWeightsBuilder] eval cache size mismatch");

            for (int i = 0; i < n; ++i)
            {
                const int e0 = offEval[static_cast<std::size_t>(i)];
                const int e1 = offEval[static_cast<std::size_t>(i) + 1u];
                if (e0 < 0 || e1 < e0)
                    throw std::runtime_error("[OrbitAveragedWeightsBuilder] eval cache invalid offsets");

                const double coeff = coeffVals[i];
                for (int e = e0; e < e1; ++e)
                {
                    const std::size_t ee = static_cast<std::size_t>(e);
                    if (ee >= ciEval.size())
                        throw std::runtime_error("[OrbitAveragedWeightsBuilder] eval cache nnz out of range");
                    const double wv = coeff * wEval[ee];
                    if (wv == 0.0)
                        continue;
                    dest.emplace_back(ciEval[ee], wv);
                }
            }
        };

        auto flushPendingGauss2 = [&]() {
            const int nPending = static_cast<int>(pendingTVals.size());
            if (nPending <= 0)
                return;
            evalBatch(pendingTVals.data(),
                      pendingCoeffVals.data(),
                      pendingHostVals.data(),
                      nPending,
                      scratchTrips);
            pendingTVals.clear();
            pendingCoeffVals.clear();
            pendingHostVals.clear();
        };

        // Adaptive recursion on [a,b].
        struct Interval
        {
            double a;
            double b;
            int depth;
        };

        std::vector<Interval> stack;
        stack.reserve(static_cast<std::size_t>(std::max(8, opt.adaptiveMaxRefineDepth + 2)));
        stack.push_back(Interval{0.0, dt, 0});

        const double tol = std::max(0.0, opt.adaptiveTol);
        const int maxDepth = std::max(0, opt.adaptiveMaxRefineDepth);
        const bool splitAtMaxDepth = opt.adaptiveCrossSplitAtMaxDepth;
        const int splitMaxIters = std::max(1, opt.adaptiveCrossSplitMaxIters);
        int splitBudget = std::max(0, opt.adaptiveCrossSplitMaxSplits);

        while (!stack.empty())
        {
            const Interval it = stack.back();
            stack.pop_back();
            const double a = it.a;
            const double b = it.b;
            const double dur = b - a;
            if (!(dur > 0.0))
                continue;

            const double mid = 0.5 * (a + b);
            const double half = 0.5 * dur;
            const double scale = half / dt;

            if (opt.adaptiveCheapRefine)
            {
                // Cheap "补刀" criterion: if host-cell labels differ at probe
                // points inside [a,b], bisect. Otherwise accept and integrate
                // with Gauss2 only.
                const double tL = mid + half * gl2_x[0];
                const double tR = mid + half * gl2_x[1];
                const int lblL = probeHostLabelAt(tL);
                const int lblM = probeHostLabelAt(mid);
                int lblR = -1;
                bool sameCell = false;
                if (lblL == lblM)
                {
                    lblR = probeHostLabelAt(tR);
                    sameCell = (lblM == lblR);
                }

                const bool cutSensitiveSameCell =
                    sameCell && (isLabelCutCell(lblL) || isLabelCutCell(lblM) || isLabelCutCell(lblR));

                if (sameCell && !cutSensitiveSameCell)
                {
                    acceptCountOut += 1u;
                    acceptedDepthSumOut += static_cast<std::size_t>(std::max(0, it.depth));
                    const double tVals[2] = {tL, tR};
                    const double coeffVals[2] = {scale * gl2_w[0], scale * gl2_w[1]};
                    const int knownHosts[2] = {lblL, lblR};
                    pendingTVals.push_back(tVals[0]);
                    pendingTVals.push_back(tVals[1]);
                    pendingCoeffVals.push_back(coeffVals[0]);
                    pendingCoeffVals.push_back(coeffVals[1]);
                    pendingHostVals.push_back(knownHosts[0]);
                    pendingHostVals.push_back(knownHosts[1]);
                    if (pendingTVals.size() >= 128u)
                        flushPendingGauss2();
                    continue;
                }

                if (it.depth >= maxDepth)
                {
                    bool didSplit = false;
                    if (splitAtMaxDepth && splitBudget > 0)
                    {
                        const Vector2D xa = posAtTime(a);
                        const Vector2D xb = posAtTime(b);
                        const Vector2D xm = posAtTime(mid);
                        const int lblA = probeHostLabelAt(a);
                        const int lblB = probeHostLabelAt(b);
                        const bool sheathGate = nearDomainWall(xa) || nearDomainWall(xm) || nearDomainWall(xb) ||
                                                isLabelCutCell(lblA) || isLabelCutCell(lblM) || isLabelCutCell(lblB);

                        if (sheathGate)
                        {
                            double lo = a;
                            double hi = b;
                            int lblLo = lblA;
                            int lblHi = lblB;
                            bool haveBracket = false;

                            if (lblL != lblM)
                            {
                                lo = tL;
                                hi = mid;
                                lblLo = lblL;
                                lblHi = lblM;
                                haveBracket = true;
                            }
                            else
                            {
                                if (lblR < 0)
                                    lblR = probeHostLabelAt(tR);
                            }

                            if (lblR < 0)
                                lblR = probeHostLabelAt(tR);

                            if (lblM != lblR)
                            {
                                lo = mid;
                                hi = tR;
                                lblLo = lblM;
                                lblHi = lblR;
                                haveBracket = true;
                            }
                            else if (lblL != lblR)
                            {
                                lo = tL;
                                hi = tR;
                                lblLo = lblL;
                                lblHi = lblR;
                                haveBracket = true;
                            }
                            else if (lblA != lblB)
                            {
                                haveBracket = true;
                            }

                            if (haveBracket && lblLo != lblHi)
                            {
                                for (int itBis = 0; itBis < splitMaxIters; ++itBis)
                                {
                                    const double m = 0.5 * (lo + hi);
                                    const int lblM2 = probeHostLabelAt(m);
                                    if (lblM2 == lblLo)
                                        lo = m;
                                    else
                                        hi = m;
                                }

                                const double tStar = 0.5 * (lo + hi);
                                const double epsT = std::max(1e-14 * std::max(1.0, dt), 1e-15);
                                if (tStar > a + epsT && tStar < b - epsT)
                                {
                                    refineCountOut += 1u;
                                    splitBudget -= 1;
                                    stack.push_back(Interval{tStar, b, it.depth});
                                    stack.push_back(Interval{a, tStar, it.depth});
                                    didSplit = true;
                                }
                            }
                        }
                    }

                    if (didSplit)
                        continue;

                    hitMaxDepthCountOut += 1u;
                    acceptedDepthSumOut += static_cast<std::size_t>(std::max(0, it.depth));
                    if (lblR < 0)
                        lblR = probeHostLabelAt(tR);
                    const double tVals[2] = {tL, tR};
                    const double coeffVals[2] = {scale * gl2_w[0], scale * gl2_w[1]};
                    const int knownHosts[2] = {lblL, lblR};
                    pendingTVals.push_back(tVals[0]);
                    pendingTVals.push_back(tVals[1]);
                    pendingCoeffVals.push_back(coeffVals[0]);
                    pendingCoeffVals.push_back(coeffVals[1]);
                    pendingHostVals.push_back(knownHosts[0]);
                    pendingHostVals.push_back(knownHosts[1]);
                    if (pendingTVals.size() >= 128u)
                        flushPendingGauss2();
                    continue;
                }

                refineCountOut += 1u;
                const double m = mid;
                stack.push_back(Interval{m, b, it.depth + 1});
                stack.push_back(Interval{a, m, it.depth + 1});
                continue;
            }

            trips2.clear();
            trips4.clear();

            // I2
            {
                const double tVals[2] = {mid + half * gl2_x[0], mid + half * gl2_x[1]};
                const double coeffVals[2] = {scale * gl2_w[0], scale * gl2_w[1]};
                evalBatch(tVals, coeffVals, nullptr, 2, trips2);
            }

            // I4
            {
                const double tVals[4] = {
                    mid + half * gl4_x[0],
                    mid + half * gl4_x[1],
                    mid + half * gl4_x[2],
                    mid + half * gl4_x[3]
                };
                const double coeffVals[4] = {
                    scale * gl4_w[0],
                    scale * gl4_w[1],
                    scale * gl4_w[2],
                    scale * gl4_w[3]
                };
                evalBatch(tVals, coeffVals, nullptr, 4, trips4);
            }

            reduceTripsToCSR(trips2, tmpCell2, tmpW2);
            reduceTripsToCSR(trips4, tmpCell4, tmpW4);
            const double err = l1DiffCSR(tmpCell2, tmpW2, tmpCell4, tmpW4);

            if (err <= tol)
            {
                // Accept higher-order estimate.
                acceptCountOut += 1u;
                acceptedDepthSumOut += static_cast<std::size_t>(std::max(0, it.depth));
                scratchTrips.insert(scratchTrips.end(), trips4.begin(), trips4.end());
                continue;
            }

            if (it.depth >= maxDepth)
            {
                hitMaxDepthCountOut += 1u;
                acceptedDepthSumOut += static_cast<std::size_t>(std::max(0, it.depth));
                scratchTrips.insert(scratchTrips.end(), trips4.begin(), trips4.end());
                continue;
            }

            refineCountOut += 1u;
            // Bisect; push right then left to keep left processed first (stack LIFO).
            const double m = mid;
            stack.push_back(Interval{m, b, it.depth + 1});
            stack.push_back(Interval{a, m, it.depth + 1});
        }

        flushPendingGauss2();

        reduceTripsToCSRStamp(scratchTrips,
                      outCellLocal,
                      outWLocal,
                      accumByCell,
                      stampByCell,
                      touchedCells,
                      stampToken,
                      opt.deterministic);

        if (opt.renormalize)
        {
            double sW = 0.0;
            for (double wv : outWLocal)
                sW += wv;
            if (sW > 0.0)
            {
                const double inv = 1.0 / sW;
                for (double& wv : outWLocal)
                    wv *= inv;
            }
            else
            {
                int cellLabel = hostCellLabelByParticle[p];
                int idxFallback = 0;
                if (cellLabel >= 0)
                {
                    const int ciFallback = mesh.cells(cellLabel).index();
                    if (ciFallback >= 0)
                        idxFallback = ciFallback;
                }
                outCellLocal.push_back(idxFallback);
                outWLocal.push_back(1.0);
            }
        }

        sampleCountOut += localSampleCount;
        sampleFallbackOut += localSampleFallback;
        sampleNnzOut += localSampleNnz;

        if (kOrbitAvgFineTimingEnabled)
        {
            evalTimeOut += evalTimeLocal;
            const double particleTotal = secondsBetween(tParticle0, Clock::now());
            accumTimeOut += std::max(0.0, particleTotal - evalTimeLocal);
        }
    };

    std::size_t sampleCountTotal = 0u;
    std::size_t sampleFallbackTotal = 0u;
    std::size_t sampleNnzTotal = 0u;
    std::size_t acceptCountTotal = 0u;
    std::size_t refineCountTotal = 0u;
    std::size_t hitMaxDepthCountTotal = 0u;
    std::size_t acceptedDepthSumTotal = 0u;
    double fusedEvalTotal_s = 0.0;
    double fusedAccumTotal_s = 0.0;

#ifdef _OPENMP
    if (numParticles >= 64u)
    {
        const int maxThreads = std::max(1, omp_get_max_threads());
        const int activeThreads = std::max(1, std::min<int>(maxThreads, static_cast<int>(numParticles)));
        std::vector<std::size_t> sampleCountByThread(static_cast<std::size_t>(activeThreads), 0u);
        std::vector<std::size_t> sampleFallbackByThread(static_cast<std::size_t>(activeThreads), 0u);
        std::vector<std::size_t> sampleNnzByThread(static_cast<std::size_t>(activeThreads), 0u);
        std::vector<std::size_t> acceptByThread(static_cast<std::size_t>(activeThreads), 0u);
        std::vector<std::size_t> refineByThread(static_cast<std::size_t>(activeThreads), 0u);
        std::vector<std::size_t> hitMaxByThread(static_cast<std::size_t>(activeThreads), 0u);
        std::vector<std::size_t> depthSumByThread(static_cast<std::size_t>(activeThreads), 0u);
        std::vector<double> evalTimeByThread(static_cast<std::size_t>(activeThreads), 0.0);
        std::vector<double> accumTimeByThread(static_cast<std::size_t>(activeThreads), 0.0);

        auto& reducedCellByThread = ws.reducedCellByThread;
        auto& reducedWByThread = ws.reducedWByThread;
        auto& reducedOffsetsByThread = ws.reducedOffsetsByThread;
        auto& reducedThreadBegin = ws.reducedThreadBegin;
        auto& reducedThreadEnd = ws.reducedThreadEnd;
        auto& reducedChunkNnzHint = ws.reducedChunkNnzHint;
        reducedCellByThread.resize(activeThreads);
        reducedWByThread.resize(activeThreads);
        reducedOffsetsByThread.resize(activeThreads);
        reducedThreadBegin.resize(activeThreads);
        reducedThreadEnd.resize(activeThreads);
        if (reducedChunkNnzHint.size() != static_cast<std::size_t>(activeThreads))
            reducedChunkNnzHint.assign(static_cast<std::size_t>(activeThreads), 0u);

        std::vector<std::vector<int>>().swap(ws.reducedCellByParticle);
        std::vector<std::vector<double>>().swap(ws.reducedWByParticle);

        std::exception_ptr firstError;
        // Schedule every logical chunk even if the runtime forms a smaller team.
        #pragma omp parallel for num_threads(activeThreads) schedule(static)
        for (int tid = 0; tid < activeThreads; ++tid)
        {
            const std::size_t pBeg = (numParticles * static_cast<std::size_t>(tid)) / static_cast<std::size_t>(activeThreads);
            const std::size_t pEnd = (numParticles * static_cast<std::size_t>(tid + 1)) / static_cast<std::size_t>(activeThreads);
            reducedThreadBegin[static_cast<std::size_t>(tid)] = pBeg;
            reducedThreadEnd[static_cast<std::size_t>(tid)] = pEnd;

            std::size_t localSampleCount = 0u;
            std::size_t localSampleFallback = 0u;
            std::size_t localSampleNnz = 0u;
            std::size_t localAccept = 0u;
            std::size_t localRefine = 0u;
            std::size_t localHitMax = 0u;
            std::size_t localDepthSum = 0u;
            double localEvalTime = 0.0;
            double localAccumTime = 0.0;

            CouplingWeightsCache evalCache;
            std::vector<Vector2D> evalPos(1u);
            std::vector<int> evalHostLabel(1u, -1);
            std::vector<double> hpSingle(1u, hpConst);
            CouplingWeightsCache::BuildOptions evalOpt = buildOpt;
            evalOpt.hostCellLabels = &evalHostLabel;
            evalOpt.parallelBuild = false;
            evalOpt.recordHpStats = false;
            evalOpt.recordPerParticleHp = false;
            evalOpt.recordPerParticleCompositeBand = false;
            evalOpt.fallbackToHostOnZeroWeight = false;

            std::vector<std::pair<int, double>> scratchTrips;
            scratchTrips.reserve(256);
            std::vector<std::pair<int, double>> trips2;
            std::vector<std::pair<int, double>> trips4;
            trips2.reserve(128);
            trips4.reserve(256);
            std::vector<int> tmpCell2, tmpCell4;
            std::vector<double> tmpW2, tmpW4;
            tmpCell2.reserve(64);
            tmpW2.reserve(64);
            tmpCell4.reserve(64);
            tmpW4.reserve(64);

            std::vector<int> outCellLocal;
            std::vector<double> outWLocal;
            outCellLocal.reserve(64);
            outWLocal.reserve(64);

            std::vector<double> accumByCell(static_cast<std::size_t>(mesh.cellIndexNum()), 0.0);
            std::vector<std::uint32_t> stampByCell(static_cast<std::size_t>(mesh.cellIndexNum()), 0u);
            std::vector<int> touchedCells;
            touchedCells.reserve(256);
            std::uint32_t stampToken = 0u;

            auto& cellFlat = reducedCellByThread[static_cast<std::size_t>(tid)];
            auto& wFlat = reducedWByThread[static_cast<std::size_t>(tid)];
            auto& localOffsets = reducedOffsetsByThread[static_cast<std::size_t>(tid)];
            cellFlat.clear();
            wFlat.clear();
            const std::size_t chunkHint = reducedChunkNnzHint[static_cast<std::size_t>(tid)];
            if (chunkHint > 0u)
            {
                cellFlat.reserve(chunkHint);
                wFlat.reserve(chunkHint);
            }
            localOffsets.clear();
            localOffsets.reserve((pEnd - pBeg) + 1u);
            localOffsets.push_back(0);

            try
            {
                for (std::size_t p = pBeg; p < pEnd; ++p)
                {
                    buildParticleAdaptiveFused(p,
                                               evalCache,
                                               evalPos,
                                               evalHostLabel,
                                               hpSingle,
                                               evalOpt,
                                               scratchTrips,
                                               trips2,
                                               trips4,
                                               tmpCell2,
                                               tmpW2,
                                               tmpCell4,
                                               tmpW4,
                                               outCellLocal,
                                               outWLocal,
                                               accumByCell,
                                               stampByCell,
                                               touchedCells,
                                               stampToken,
                                               localSampleCount,
                                               localSampleFallback,
                                               localSampleNnz,
                                               localAccept,
                                               localRefine,
                                               localHitMax,
                                               localDepthSum,
                                               localEvalTime,
                                               localAccumTime);
                    particleNnz[p] = static_cast<int>(outCellLocal.size());
                    if (!outCellLocal.empty())
                    {
                        const std::size_t oldN = cellFlat.size();
                        const std::size_t nLocal = outCellLocal.size();
                        cellFlat.resize(oldN + nLocal);
                        wFlat.resize(oldN + nLocal);
                        std::memcpy(cellFlat.data() + oldN, outCellLocal.data(), nLocal * sizeof(int));
                        std::memcpy(wFlat.data() + oldN, outWLocal.data(), nLocal * sizeof(double));
                    }
                    localOffsets.push_back(static_cast<int>(cellFlat.size()));
                }

                reducedChunkNnzHint[static_cast<std::size_t>(tid)] = std::max<std::size_t>(cellFlat.size(), 32u);
                sampleCountByThread[static_cast<std::size_t>(tid)] = localSampleCount;
                sampleFallbackByThread[static_cast<std::size_t>(tid)] = localSampleFallback;
                sampleNnzByThread[static_cast<std::size_t>(tid)] = localSampleNnz;
                acceptByThread[static_cast<std::size_t>(tid)] = localAccept;
                refineByThread[static_cast<std::size_t>(tid)] = localRefine;
                hitMaxByThread[static_cast<std::size_t>(tid)] = localHitMax;
                depthSumByThread[static_cast<std::size_t>(tid)] = localDepthSum;
                evalTimeByThread[static_cast<std::size_t>(tid)] = localEvalTime;
                accumTimeByThread[static_cast<std::size_t>(tid)] = localAccumTime;
            }
            catch (...)
            {
                #pragma omp critical(lamp_orbitavg_adapt_error)
                {
                    if (!firstError)
                        firstError = std::current_exception();
                }
            }
        }
        if (firstError)
            std::rethrow_exception(firstError);

        for (int t = 0; t < activeThreads; ++t)
        {
            const std::size_t tt = static_cast<std::size_t>(t);
            sampleCountTotal += sampleCountByThread[tt];
            sampleFallbackTotal += sampleFallbackByThread[tt];
            sampleNnzTotal += sampleNnzByThread[tt];
            acceptCountTotal += acceptByThread[tt];
            refineCountTotal += refineByThread[tt];
            hitMaxDepthCountTotal += hitMaxByThread[tt];
            acceptedDepthSumTotal += depthSumByThread[tt];
            fusedEvalTotal_s = std::max(fusedEvalTotal_s, evalTimeByThread[tt]);
            fusedAccumTotal_s = std::max(fusedAccumTotal_s, accumTimeByThread[tt]);
        }
    }
    else
#endif
    {
        auto& reducedCellByParticle = ws.reducedCellByParticle;
        auto& reducedWByParticle = ws.reducedWByParticle;
        reducedCellByParticle.resize(numParticles);
        reducedWByParticle.resize(numParticles);

        CouplingWeightsCache evalCache;
        std::vector<Vector2D> evalPos(1u);
        std::vector<int> evalHostLabel(1u, -1);
        std::vector<double> hpSingle(1u, hpConst);
        CouplingWeightsCache::BuildOptions evalOpt = buildOpt;
        evalOpt.hostCellLabels = &evalHostLabel;
        evalOpt.parallelBuild = false;
        evalOpt.recordHpStats = false;
        evalOpt.recordPerParticleHp = false;
        evalOpt.recordPerParticleCompositeBand = false;
        evalOpt.fallbackToHostOnZeroWeight = false;

        std::vector<std::pair<int, double>> scratchTrips;
        scratchTrips.reserve(256);
        std::vector<std::pair<int, double>> trips2;
        std::vector<std::pair<int, double>> trips4;
        trips2.reserve(128);
        trips4.reserve(256);
        std::vector<int> tmpCell2, tmpCell4;
        std::vector<double> tmpW2, tmpW4;
        tmpCell2.reserve(64);
        tmpW2.reserve(64);
        tmpCell4.reserve(64);
        tmpW4.reserve(64);

        std::vector<int> outCellLocal;
        std::vector<double> outWLocal;
        outCellLocal.reserve(64);
        outWLocal.reserve(64);

        std::vector<double> accumByCell(static_cast<std::size_t>(mesh.cellIndexNum()), 0.0);
        std::vector<std::uint32_t> stampByCell(static_cast<std::size_t>(mesh.cellIndexNum()), 0u);
        std::vector<int> touchedCells;
        touchedCells.reserve(256);
        std::uint32_t stampToken = 0u;

        for (std::size_t p = 0; p < numParticles; ++p)
        {
            auto& outCellP = reducedCellByParticle[p];
            auto& outWP = reducedWByParticle[p];
            buildParticleAdaptiveFused(p,
                                       evalCache,
                                       evalPos,
                                       evalHostLabel,
                                       hpSingle,
                                       evalOpt,
                                       scratchTrips,
                                       trips2,
                                       trips4,
                                       tmpCell2,
                                       tmpW2,
                                       tmpCell4,
                                       tmpW4,
                                       outCellLocal,
                                       outWLocal,
                                      accumByCell,
                                      stampByCell,
                                      touchedCells,
                                      stampToken,
                                       sampleCountTotal,
                                       sampleFallbackTotal,
                                       sampleNnzTotal,
                                       acceptCountTotal,
                                       refineCountTotal,
                                       hitMaxDepthCountTotal,
                                       acceptedDepthSumTotal,
                                       fusedEvalTotal_s,
                                       fusedAccumTotal_s);
            outCellP = outCellLocal;
            outWP = outWLocal;
            particleNnz[p] = static_cast<int>(outCellP.size());
        }
    }

#ifdef _OPENMP
    if (numParticles >= 64u)
    {
        const int activeThreads = static_cast<int>(ws.reducedCellByThread.size());
        if (activeThreads <= 0)
            throw std::runtime_error("[OrbitAveragedWeightsBuilder] invalid adaptive fused thread chunks");

        const auto& reducedCellByThread = ws.reducedCellByThread;
        const auto& reducedWByThread = ws.reducedWByThread;
        const auto& reducedOffsetsByThread = ws.reducedOffsetsByThread;
        const auto& reducedThreadBegin = ws.reducedThreadBegin;
        const auto& reducedThreadEnd = ws.reducedThreadEnd;

        std::vector<std::size_t> threadBase(static_cast<std::size_t>(activeThreads) + 1u, 0u);
        for (int t = 0; t < activeThreads; ++t)
        {
            threadBase[static_cast<std::size_t>(t) + 1u] =
                threadBase[static_cast<std::size_t>(t)] + reducedCellByThread[static_cast<std::size_t>(t)].size();
        }
        outNnz = threadBase[static_cast<std::size_t>(activeThreads)];

        for (int t = 0; t < activeThreads; ++t)
        {
            const std::size_t pBeg = reducedThreadBegin[static_cast<std::size_t>(t)];
            const std::size_t pEnd = reducedThreadEnd[static_cast<std::size_t>(t)];
            const auto& localOffsets = reducedOffsetsByThread[static_cast<std::size_t>(t)];
            if (localOffsets.size() != (pEnd - pBeg + 1u))
                throw std::runtime_error("[OrbitAveragedWeightsBuilder] inconsistent adaptive fused chunk offsets");

            const int base = static_cast<int>(threadBase[static_cast<std::size_t>(t)]);
            for (std::size_t i = 0; i < localOffsets.size(); ++i)
                outOffsets[pBeg + i] = base + localOffsets[i];
        }

        outCell.resize(outNnz);
        outW.resize(outNnz, 0.0);

        std::exception_ptr firstError;
        #pragma omp parallel for num_threads(activeThreads) schedule(static)
        for (int t = 0; t < activeThreads; ++t)
        {
            try
            {
                const std::size_t base = threadBase[static_cast<std::size_t>(t)];
                const auto& cellFlat = reducedCellByThread[static_cast<std::size_t>(t)];
                const auto& wFlat = reducedWByThread[static_cast<std::size_t>(t)];
                if (cellFlat.size() != wFlat.size())
                    throw std::runtime_error("[OrbitAveragedWeightsBuilder] inconsistent adaptive fused chunk data");
                if (!cellFlat.empty())
                {
                    std::memcpy(outCell.data() + base, cellFlat.data(), cellFlat.size() * sizeof(int));
                    std::memcpy(outW.data() + base, wFlat.data(), wFlat.size() * sizeof(double));
                }
            }
            catch (...)
            {
                #pragma omp critical(lamp_orbitavg_adapt_copy_error)
                {
                    if (!firstError)
                        firstError = std::current_exception();
                }
            }
        }
        if (firstError)
            std::rethrow_exception(firstError);
    }
    else
#endif
    {
        const auto& reducedCellByParticle = ws.reducedCellByParticle;
        const auto& reducedWByParticle = ws.reducedWByParticle;

        for (std::size_t p = 0; p < numParticles; ++p)
            outOffsets[p + 1u] = outOffsets[p] + particleNnz[p];

        outNnz = static_cast<std::size_t>(outOffsets.back());
        outCell.resize(outNnz);
        outW.resize(outNnz, 0.0);
        for (std::size_t p = 0; p < numParticles; ++p)
        {
            const auto& localCell = reducedCellByParticle[p];
            const auto& localW = reducedWByParticle[p];
            if (localCell.size() != static_cast<std::size_t>(particleNnz[p]) || localCell.size() != localW.size())
                throw std::runtime_error("[OrbitAveragedWeightsBuilder] adaptive fused CSR count mismatch");
            const std::size_t writeBeg = static_cast<std::size_t>(outOffsets[p]);
            const std::size_t nLocal = localCell.size();
            if (nLocal > 0u)
            {
                std::memcpy(outCell.data() + writeBeg, localCell.data(), nLocal * sizeof(int));
                std::memcpy(outW.data() + writeBeg, localW.data(), nLocal * sizeof(double));
            }
        }
    }

    if (timing)
    {
        timing->fusedEval_s += fusedEvalTotal_s;
        timing->fusedAccum_s += fusedAccumTotal_s;
        timing->sampleCount += sampleCountTotal;
        timing->sampleCacheFallbackCount += sampleFallbackTotal;
        timing->sampleCacheNnz += sampleNnzTotal;
        timing->outNnz += outNnz;
        timing->adaptAcceptCount += acceptCountTotal;
        timing->adaptRefineCount += refineCountTotal;
        timing->adaptHitMaxDepthCount += hitMaxDepthCountTotal;
        timing->adaptAcceptedDepthSum += acceptedDepthSumTotal;
        timing->total_s += secondsBetween(tTotal0, Clock::now());
    }

    out.assignCSR(numParticles, std::move(outOffsets), std::move(outCell), std::move(outW));
    return true;
}

bool buildOrbitAveragedCacheAdaptiveGauss24(const BackgroundMesh& mesh,
                                           const std::vector<double>& hp,
                                           const CouplingWeightsCache::BuildOptions& buildOpt,
                                           std::size_t numParticles,
                                           const std::vector<Vector2D>& x0,
                                           const std::vector<Vector2D>& x1,
                                           const std::vector<int>& hostCellLabelByParticle,
                                           double dt,
                                           CouplingWeightsCache& out,
                                           const OrbitAveragedWeightsBuilderOptions& opt,
                                           OrbitAveragedWeightsBuilderTiming* timing)
{
    return buildOrbitAveragedCacheAdaptiveGauss24(mesh,
                                                 hp,
                                                 buildOpt,
                                                 numParticles,
                                                 x0,
                                                 x1,
                                                 hostCellLabelByParticle,
                                                 dt,
                                                 out,
                                                 nullptr,
                                                 opt,
                                                 timing);
}

bool buildOrbitAveragedCacheGauss2(const BackgroundMesh& mesh,
                                  const std::vector<double>& hp,
                                  const CouplingWeightsCache::BuildOptions& buildOpt,
                                  const std::size_t numParticles,
                                  const std::vector<int>& segmentOffsets,
                                  const std::vector<OrbitSegment>& segments,
                                  const std::vector<int>& hostCellLabelByParticle,
                                  const double dt,
                                  CouplingWeightsCache& out,
                                  OrbitAveragedWeightsBuilderWorkspace* workspace,
                                  const OrbitAveragedWeightsBuilderOptions& opt,
                                  OrbitAveragedWeightsBuilderTiming* timing)
{
    using Clock = std::chrono::steady_clock;
    auto secondsBetween = [](const Clock::time_point& a, const Clock::time_point& b) -> double {
        return std::chrono::duration_cast<std::chrono::duration<double>>(b - a).count();
    };
    const auto tTotal0 = Clock::now();
    if (timing)
        *timing = OrbitAveragedWeightsBuilderTiming{};

    OrbitAveragedWeightsBuilderWorkspace localWorkspace;
    OrbitAveragedWeightsBuilderWorkspace& ws = workspace ? *workspace : localWorkspace;

    if (dt <= 0.0)
        return false;
    if (segmentOffsets.size() != numParticles + 1u)
        return false;
    if (hostCellLabelByParticle.size() != numParticles)
        return false;
    if (hp.empty())
        return false;

    const bool periodicX = buildOpt.periodicX;
    const bool periodicY = buildOpt.periodicY;
    constexpr double epsWrap = 1e-12;

    // Patch mode uses fused-style per-particle accumulation regardless of opt.useFusedPath,
    // since the number of samples per segment becomes data-dependent.
    if (opt.useFusedPath || opt.segmentPatchEnabled)
    {
        if (timing)
            timing->fusedMode = opt.segmentPatchEnabled ? 4 : 1;

        const bool hpIsConst = (hp.size() == 1u);
        const double hpConst = hpIsConst ? hp.front() : 0.0;
        auto hpForParticle = [&](std::size_t p) -> double {
            if (hpIsConst)
                return hpConst;
            if (p < hp.size())
                return hp[p];
            return hp.back();
        };

        auto& particleNnz = ws.particleNnz;
        particleNnz.assign(numParticles, 0);
        auto& outOffsets = ws.outOffsets;
        outOffsets.assign(numParticles + 1u, 0);
        std::vector<int> outCell;
        std::vector<double> outW;
        std::size_t outNnz = 0u;

        // Gauss-Legendre nodes/weights (mapped by mid+half*x).
        constexpr double gl2_x[2] = {-0.5773502691896257645, 0.5773502691896257645};
        constexpr double gl2_w[2] = {1.0, 1.0};
        constexpr double gl4_x[4] = {-0.8611363115940525752,
                                     -0.3399810435848562648,
                                     0.3399810435848562648,
                                     0.8611363115940525752};
        constexpr double gl4_w[4] = {0.3478548451374538574,
                                     0.6521451548625461426,
                                     0.6521451548625461426,
                                     0.3478548451374538574};

        auto hostLabelForSample = [&](const Vector2D& x, int fallbackLabel) -> int {
            int lbl = mesh.findCellLabel(x);
            if (lbl < 0 && (periodicX || periodicY))
            {
                const Vector2D xw = wrapToDomain(mesh, x, periodicX, periodicY, epsWrap);
                lbl = mesh.findCellLabel(xw);
            }
            if (lbl >= 0)
                return lbl;
            if (fallbackLabel >= 0)
                return fallbackLabel;
            return 0;
        };

        auto buildParticleFused = [&](std::size_t p,
                                      CouplingWeightsCache& evalCache,
                                      std::vector<Vector2D>& evalPos,
                                      std::vector<int>& evalHostLabel,
                                      std::vector<double>& hpSingle,
                                      CouplingWeightsCache::BuildOptions& evalOpt,
                                      std::vector<std::pair<int, double>>& scratchTrips,
                                      std::vector<int>& outCellLocal,
                                      std::vector<double>& outWLocal,
                                      std::vector<double>& accumByCell,
                                      std::vector<std::uint32_t>& stampByCell,
                                      std::vector<int>& touchedCells,
                                      std::uint32_t& stampToken,
                                      std::size_t& sampleCountOut,
                                      std::size_t& sampleFallbackOut,
                                      std::size_t& sampleNnzOut,
                                      std::size_t& acceptCountOut,
                                      std::size_t& refineCountOut,
                                      std::size_t& hitMaxDepthCountOut,
                                      std::size_t& acceptedDepthSumOut,
                                      double& evalTimeOut,
                                      double& accumTimeOut)
        {
            const auto tParticle0 = kOrbitAvgFineTimingEnabled ? Clock::now() : Clock::time_point{};
            double evalTimeLocal = 0.0;
            std::size_t localSampleCount = 0u;
            std::size_t localSampleFallback = 0u;
            std::size_t localSampleNnz = 0u;

            scratchTrips.clear();
            outCellLocal.clear();
            outWLocal.clear();

            const int s0 = segmentOffsets[p];
            const int s1 = segmentOffsets[p + 1u];
            if (s0 < 0 || s1 < s0)
                throw std::runtime_error("[OrbitAveragedWeightsBuilder] invalid segment offsets");

            for (int si = s0; si < s1; ++si)
            {
                if (static_cast<std::size_t>(si) >= segments.size())
                    throw std::runtime_error("[OrbitAveragedWeightsBuilder] segment index out of range");
                const OrbitSegment& seg = segments[static_cast<std::size_t>(si)];
                const double a = seg.t0;
                const double b = seg.t1;
                const double dur = b - a;
                if (!(dur > 0.0))
                    continue;

                const Vector2D x0s = seg.x0;
                const Vector2D x1s = unwrapToNearestImage(mesh, x0s, seg.x1, periodicX, periodicY);

                const double mid = 0.5 * (a + b);
                const double half = 0.5 * dur;

                auto posAt = [&](double t) -> Vector2D {
                    const double s = (dur > 0.0) ? ((t - a) / dur) : 0.0;
                    const Vector2D xr = lerp(x0s, x1s, std::clamp(s, 0.0, 1.0));
                    return wrapToDomain(mesh, xr, periodicX, periodicY, epsWrap);
                };

                auto evalBatch = [&](const double* tVals,
                                     const double* coeffVals,
                                     int n,
                                     std::vector<std::pair<int, double>>& dest,
                                     int* outHostLabels = nullptr) {
                    if (n <= 0)
                        return;

                    evalPos.resize(static_cast<std::size_t>(n));
                    evalHostLabel.resize(static_cast<std::size_t>(n));
                    if (!hpIsConst)
                        hpSingle.resize(static_cast<std::size_t>(n));

                    for (int i = 0; i < n; ++i)
                    {
                        evalPos[static_cast<std::size_t>(i)] = posAt(tVals[i]);

                        int host = -1;
                        if (opt.segmentPatchEnabled)
                        {
                            // Patch mode: host label follows the sample position.
                            host = hostLabelForSample(evalPos[static_cast<std::size_t>(i)], hostCellLabelByParticle[p]);
                        }
                        else
                        {
                            // Baseline: keep the cheap host-cell hint from segmenter.
                            host = seg.cellLabel;
                            if (host < 0)
                                host = hostCellLabelByParticle[p];
                            if (host < 0)
                                host = 0;
                        }

                        evalHostLabel[static_cast<std::size_t>(i)] = host;
                        if (outHostLabels)
                            outHostLabels[i] = host;
                        if (!hpIsConst)
                            hpSingle[static_cast<std::size_t>(i)] = hpForParticle(p);
                    }

                    const auto tEval0 = kOrbitAvgFineTimingEnabled ? Clock::now() : Clock::time_point{};
                    const bool prevParallelBuild = evalOpt.parallelBuild;
                    // Fused orbit-averaged path uses thread-local scratch and reduction-heavy flows;
                    // keep inner cache build serial here to avoid unstable nested/combined parallel behavior.
                    evalOpt.parallelBuild = false;
                    if (hpIsConst)
                        evalCache.build(mesh, evalPos, hp, evalOpt);
                    else
                        evalCache.build(mesh, evalPos, hpSingle, evalOpt);
                    evalOpt.parallelBuild = prevParallelBuild;
                    if (kOrbitAvgFineTimingEnabled)
                        evalTimeLocal += secondsBetween(tEval0, Clock::now());

                    localSampleCount += static_cast<std::size_t>(n);
                    localSampleFallback += evalCache.lastBuildFallbackCount();
                    localSampleNnz += evalCache.nnz();

                    const auto& offEval = evalCache.offsets();
                    const auto& ciEval = evalCache.cellIndices();
                    const auto& wEval = evalCache.weights();
                    if (ciEval.size() != wEval.size() || offEval.size() != static_cast<std::size_t>(n + 1))
                        throw std::runtime_error("[OrbitAveragedWeightsBuilder] eval cache size mismatch");

                    for (int i = 0; i < n; ++i)
                    {
                        const int e0 = offEval[static_cast<std::size_t>(i)];
                        const int e1 = offEval[static_cast<std::size_t>(i) + 1u];
                        if (e0 < 0 || e1 < e0)
                            throw std::runtime_error("[OrbitAveragedWeightsBuilder] eval cache invalid offsets");

                        const double coeff = coeffVals[i];
                        for (int e = e0; e < e1; ++e)
                        {
                            const std::size_t ee = static_cast<std::size_t>(e);
                            if (ee >= ciEval.size())
                                throw std::runtime_error("[OrbitAveragedWeightsBuilder] eval cache nnz out of range");
                            const double wv = coeff * wEval[ee];
                            if (wv == 0.0)
                                continue;
                            dest.emplace_back(ciEval[ee], wv);
                        }
                    }
                };

                if (!opt.segmentPatchEnabled)
                {
                    // Baseline: fixed Gauss2 per segment.
                    const double scale = half / dt;
                    const double tVals[2] = {mid + half * gl2_x[0], mid + half * gl2_x[1]};
                    const double coeffVals[2] = {scale * gl2_w[0], scale * gl2_w[1]};
                    evalBatch(tVals, coeffVals, 2, scratchTrips);
                    continue;
                }

                // Hybrid patching: only refine when the segment appears to cross leaf cells.
                struct Interval
                {
                    double a;
                    double b;
                    int depth;
                };

                const double tol = std::max(0.0, opt.segmentPatchTol);
                const int maxDepth = std::max(0, opt.segmentPatchMaxRefineDepth);

                std::vector<Interval> stack;
                stack.reserve(static_cast<std::size_t>(std::max(8, maxDepth + 2)));
                stack.push_back(Interval{a, b, 0});

                std::vector<std::pair<int, double>> trips2;
                std::vector<std::pair<int, double>> trips4;
                std::vector<int> tmpCell2, tmpCell4;
                std::vector<double> tmpW2, tmpW4;
                trips2.reserve(128);
                trips4.reserve(256);
                tmpCell2.reserve(64);
                tmpW2.reserve(64);
                tmpCell4.reserve(64);
                tmpW4.reserve(64);

                while (!stack.empty())
                {
                    const Interval it = stack.back();
                    stack.pop_back();
                    const double ia = it.a;
                    const double ib = it.b;
                    const double idur = ib - ia;
                    if (!(idur > 0.0))
                        continue;

                    const double imid = 0.5 * (ia + ib);
                    const double ihalf = 0.5 * idur;
                    const double scale = (ihalf / dt);

                    trips2.clear();
                    trips4.clear();

                    int lbl0 = -1;
                    int lbl1 = -2;
                    {
                        const double tVals[2] = {imid + ihalf * gl2_x[0], imid + ihalf * gl2_x[1]};
                        const double coeffVals[2] = {scale * gl2_w[0], scale * gl2_w[1]};
                        int hosts[2] = {-1, -1};
                        evalBatch(tVals, coeffVals, 2, trips2, hosts);
                        lbl0 = hosts[0];
                        lbl1 = hosts[1];
                    }

                    // If the two Gauss2 nodes land in the same leaf cell, accept I2.
                    // This makes patching effectively "rare": extra work only happens
                    // when the interval likely crosses cell boundaries.
                    if (lbl0 == lbl1)
                    {
                        acceptCountOut += 1u;
                        acceptedDepthSumOut += static_cast<std::size_t>(std::max(0, it.depth));
                        scratchTrips.insert(scratchTrips.end(), trips2.begin(), trips2.end());
                        continue;
                    }

                    {
                        const double tVals[4] = {
                            imid + ihalf * gl4_x[0],
                            imid + ihalf * gl4_x[1],
                            imid + ihalf * gl4_x[2],
                            imid + ihalf * gl4_x[3]
                        };
                        const double coeffVals[4] = {
                            scale * gl4_w[0],
                            scale * gl4_w[1],
                            scale * gl4_w[2],
                            scale * gl4_w[3]
                        };
                        evalBatch(tVals, coeffVals, 4, trips4);
                    }

                    reduceTripsToCSR(trips2, tmpCell2, tmpW2);
                    reduceTripsToCSR(trips4, tmpCell4, tmpW4);
                    const double err = l1DiffCSR(tmpCell2, tmpW2, tmpCell4, tmpW4);

                    if (err <= tol)
                    {
                        acceptCountOut += 1u;
                        acceptedDepthSumOut += static_cast<std::size_t>(std::max(0, it.depth));
                        scratchTrips.insert(scratchTrips.end(), trips4.begin(), trips4.end());
                        continue;
                    }

                    if (it.depth >= maxDepth)
                    {
                        hitMaxDepthCountOut += 1u;
                        acceptedDepthSumOut += static_cast<std::size_t>(std::max(0, it.depth));
                        scratchTrips.insert(scratchTrips.end(), trips4.begin(), trips4.end());
                        continue;
                    }

                    refineCountOut += 1u;
                    const double m = imid;
                    stack.push_back(Interval{m, ib, it.depth + 1});
                    stack.push_back(Interval{ia, m, it.depth + 1});
                }
            }

            reduceTripsToCSRStamp(scratchTrips,
                                  outCellLocal,
                                  outWLocal,
                                  accumByCell,
                                  stampByCell,
                                  touchedCells,
                          stampToken,
                          opt.deterministic);

            if (opt.renormalize)
            {
                double sW = 0.0;
                for (double wv : outWLocal)
                    sW += wv;
                if (sW > 0.0)
                {
                    const double inv = 1.0 / sW;
                    for (double& wv : outWLocal)
                        wv *= inv;
                }
                else
                {
                    int cellLabel = hostCellLabelByParticle[p];
                    int idxFallback = 0;
                    if (cellLabel >= 0)
                    {
                        const int ciFallback = mesh.cells(cellLabel).index();
                        if (ciFallback >= 0)
                            idxFallback = ciFallback;
                    }
                    outCellLocal.push_back(idxFallback);
                    outWLocal.push_back(1.0);
                }
            }

            sampleCountOut += localSampleCount;
            sampleFallbackOut += localSampleFallback;
            sampleNnzOut += localSampleNnz;

            if (kOrbitAvgFineTimingEnabled)
            {
                evalTimeOut += evalTimeLocal;
                const double particleTotal = secondsBetween(tParticle0, Clock::now());
                accumTimeOut += std::max(0.0, particleTotal - evalTimeLocal);
            }
        };

        std::size_t sampleCountTotal = 0u;
        std::size_t sampleFallbackTotal = 0u;
        std::size_t sampleNnzTotal = 0u;
        std::size_t acceptCountTotal = 0u;
        std::size_t refineCountTotal = 0u;
        std::size_t hitMaxDepthCountTotal = 0u;
        std::size_t acceptedDepthSumTotal = 0u;
        double fusedEvalTotal_s = 0.0;
        double fusedAccumTotal_s = 0.0;

#ifdef _OPENMP
        if (numParticles >= 64u)
        {
            const int maxThreads = std::max(1, omp_get_max_threads());
            const int activeThreads = std::max(1, std::min<int>(maxThreads, static_cast<int>(numParticles)));
            std::vector<std::size_t> sampleCountByThread(static_cast<std::size_t>(activeThreads), 0u);
            std::vector<std::size_t> sampleFallbackByThread(static_cast<std::size_t>(activeThreads), 0u);
            std::vector<std::size_t> sampleNnzByThread(static_cast<std::size_t>(activeThreads), 0u);
            std::vector<std::size_t> acceptByThread(static_cast<std::size_t>(activeThreads), 0u);
            std::vector<std::size_t> refineByThread(static_cast<std::size_t>(activeThreads), 0u);
            std::vector<std::size_t> hitMaxByThread(static_cast<std::size_t>(activeThreads), 0u);
            std::vector<std::size_t> depthSumByThread(static_cast<std::size_t>(activeThreads), 0u);
            std::vector<double> evalTimeByThread(static_cast<std::size_t>(activeThreads), 0.0);
            std::vector<double> accumTimeByThread(static_cast<std::size_t>(activeThreads), 0.0);

            auto& reducedCellByThread = ws.reducedCellByThread;
            auto& reducedWByThread = ws.reducedWByThread;
            auto& reducedOffsetsByThread = ws.reducedOffsetsByThread;
            auto& reducedThreadBegin = ws.reducedThreadBegin;
            auto& reducedThreadEnd = ws.reducedThreadEnd;
            auto& reducedChunkNnzHint = ws.reducedChunkNnzHint;
            reducedCellByThread.resize(activeThreads);
            reducedWByThread.resize(activeThreads);
            reducedOffsetsByThread.resize(activeThreads);
            reducedThreadBegin.resize(activeThreads);
            reducedThreadEnd.resize(activeThreads);
            if (reducedChunkNnzHint.size() != static_cast<std::size_t>(activeThreads))
                reducedChunkNnzHint.assign(static_cast<std::size_t>(activeThreads), 0u);

            std::vector<std::vector<int>>().swap(ws.reducedCellByParticle);
            std::vector<std::vector<double>>().swap(ws.reducedWByParticle);

            std::exception_ptr firstError;
            // Schedule every logical chunk even if the runtime forms a smaller team.
            #pragma omp parallel for num_threads(activeThreads) schedule(static)
            for (int tid = 0; tid < activeThreads; ++tid)
            {
                const std::size_t pBeg = (numParticles * static_cast<std::size_t>(tid)) / static_cast<std::size_t>(activeThreads);
                const std::size_t pEnd = (numParticles * static_cast<std::size_t>(tid + 1)) / static_cast<std::size_t>(activeThreads);
                reducedThreadBegin[static_cast<std::size_t>(tid)] = pBeg;
                reducedThreadEnd[static_cast<std::size_t>(tid)] = pEnd;

                std::size_t localSampleCount = 0u;
                std::size_t localSampleFallback = 0u;
                std::size_t localSampleNnz = 0u;
                std::size_t localAccept = 0u;
                std::size_t localRefine = 0u;
                std::size_t localHitMax = 0u;
                std::size_t localDepthSum = 0u;
                double localEvalTime = 0.0;
                double localAccumTime = 0.0;

                CouplingWeightsCache evalCache;
                std::vector<Vector2D> evalPos(1u);
                std::vector<int> evalHostLabel(1u, -1);
                std::vector<double> hpSingle(1u, hpConst);
                CouplingWeightsCache::BuildOptions evalOpt = buildOpt;
                evalOpt.hostCellLabels = &evalHostLabel;
                evalOpt.parallelBuild = false;
                evalOpt.recordHpStats = false;
                evalOpt.recordPerParticleHp = false;
                evalOpt.recordPerParticleCompositeBand = false;
                evalOpt.fallbackToHostOnZeroWeight = false;

                std::vector<std::pair<int, double>> scratchTrips;
                scratchTrips.reserve(128);
                std::vector<int> outCellLocal;
                std::vector<double> outWLocal;
                outCellLocal.reserve(64);
                outWLocal.reserve(64);
                std::vector<double> accumByCell(static_cast<std::size_t>(mesh.cellIndexNum()), 0.0);
                std::vector<std::uint32_t> stampByCell(static_cast<std::size_t>(mesh.cellIndexNum()), 0u);
                std::vector<int> touchedCells;
                touchedCells.reserve(128);
                std::uint32_t stampToken = 0u;

                auto& cellFlat = reducedCellByThread[static_cast<std::size_t>(tid)];
                auto& wFlat = reducedWByThread[static_cast<std::size_t>(tid)];
                auto& localOffsets = reducedOffsetsByThread[static_cast<std::size_t>(tid)];
                cellFlat.clear();
                wFlat.clear();
                const std::size_t chunkHint = reducedChunkNnzHint[static_cast<std::size_t>(tid)];
                if (chunkHint > 0u)
                {
                    cellFlat.reserve(chunkHint);
                    wFlat.reserve(chunkHint);
                }
                localOffsets.clear();
                localOffsets.reserve((pEnd - pBeg) + 1u);
                localOffsets.push_back(0);

                try
                {
                    for (std::size_t p = pBeg; p < pEnd; ++p)
                    {
                        buildParticleFused(p,
                                           evalCache,
                                           evalPos,
                                           evalHostLabel,
                                           hpSingle,
                                           evalOpt,
                                           scratchTrips,
                                           outCellLocal,
                                           outWLocal,
                                           accumByCell,
                                           stampByCell,
                                           touchedCells,
                                           stampToken,
                                           localSampleCount,
                                           localSampleFallback,
                                           localSampleNnz,
                                           localAccept,
                                           localRefine,
                                           localHitMax,
                                           localDepthSum,
                                           localEvalTime,
                                           localAccumTime);
                        particleNnz[p] = static_cast<int>(outCellLocal.size());
                        if (!outCellLocal.empty())
                        {
                            const std::size_t oldN = cellFlat.size();
                            const std::size_t nLocal = outCellLocal.size();
                            cellFlat.resize(oldN + nLocal);
                            wFlat.resize(oldN + nLocal);
                            std::memcpy(cellFlat.data() + oldN, outCellLocal.data(), nLocal * sizeof(int));
                            std::memcpy(wFlat.data() + oldN, outWLocal.data(), nLocal * sizeof(double));
                        }
                        localOffsets.push_back(static_cast<int>(cellFlat.size()));
                    }

                    reducedChunkNnzHint[static_cast<std::size_t>(tid)] = std::max<std::size_t>(cellFlat.size(), 32u);
                    sampleCountByThread[static_cast<std::size_t>(tid)] = localSampleCount;
                    sampleFallbackByThread[static_cast<std::size_t>(tid)] = localSampleFallback;
                    sampleNnzByThread[static_cast<std::size_t>(tid)] = localSampleNnz;
                    acceptByThread[static_cast<std::size_t>(tid)] = localAccept;
                    refineByThread[static_cast<std::size_t>(tid)] = localRefine;
                    hitMaxByThread[static_cast<std::size_t>(tid)] = localHitMax;
                    depthSumByThread[static_cast<std::size_t>(tid)] = localDepthSum;
                    evalTimeByThread[static_cast<std::size_t>(tid)] = localEvalTime;
                    accumTimeByThread[static_cast<std::size_t>(tid)] = localAccumTime;
                }
                catch (...)
                {
                    #pragma omp critical(lamp_orbitavg_fused_error)
                    {
                        if (!firstError)
                            firstError = std::current_exception();
                    }
                }
            }
            if (firstError)
                std::rethrow_exception(firstError);

            for (int t = 0; t < activeThreads; ++t)
            {
                const std::size_t tt = static_cast<std::size_t>(t);
                sampleCountTotal += sampleCountByThread[tt];
                sampleFallbackTotal += sampleFallbackByThread[tt];
                sampleNnzTotal += sampleNnzByThread[tt];
                acceptCountTotal += acceptByThread[tt];
                refineCountTotal += refineByThread[tt];
                hitMaxDepthCountTotal += hitMaxByThread[tt];
                acceptedDepthSumTotal += depthSumByThread[tt];
                fusedEvalTotal_s = std::max(fusedEvalTotal_s, evalTimeByThread[tt]);
                fusedAccumTotal_s = std::max(fusedAccumTotal_s, accumTimeByThread[tt]);
            }
        }
        else
#endif
        {
            auto& reducedCellByParticle = ws.reducedCellByParticle;
            auto& reducedWByParticle = ws.reducedWByParticle;
            reducedCellByParticle.resize(numParticles);
            reducedWByParticle.resize(numParticles);

            CouplingWeightsCache evalCache;
            std::vector<Vector2D> evalPos(1u);
            std::vector<int> evalHostLabel(1u, -1);
            std::vector<double> hpSingle(1u, hpConst);
            CouplingWeightsCache::BuildOptions evalOpt = buildOpt;
            evalOpt.hostCellLabels = &evalHostLabel;
            evalOpt.parallelBuild = false;
            evalOpt.recordHpStats = false;
            evalOpt.recordPerParticleHp = false;
            evalOpt.recordPerParticleCompositeBand = false;
            evalOpt.fallbackToHostOnZeroWeight = false;

            std::vector<std::pair<int, double>> scratchTrips;
            scratchTrips.reserve(128);
            std::vector<int> outCellLocal;
            std::vector<double> outWLocal;
            outCellLocal.reserve(64);
            outWLocal.reserve(64);
            std::vector<double> accumByCell(static_cast<std::size_t>(mesh.cellIndexNum()), 0.0);
            std::vector<std::uint32_t> stampByCell(static_cast<std::size_t>(mesh.cellIndexNum()), 0u);
            std::vector<int> touchedCells;
            touchedCells.reserve(128);
            std::uint32_t stampToken = 0u;

            for (std::size_t p = 0; p < numParticles; ++p)
            {
                auto& outCellP = reducedCellByParticle[p];
                auto& outWP = reducedWByParticle[p];
                buildParticleFused(p,
                                   evalCache,
                                   evalPos,
                                   evalHostLabel,
                                   hpSingle,
                                   evalOpt,
                                   scratchTrips,
                                   outCellLocal,
                                   outWLocal,
                                   accumByCell,
                                   stampByCell,
                                   touchedCells,
                                   stampToken,
                                   sampleCountTotal,
                                   sampleFallbackTotal,
                                   sampleNnzTotal,
                                   acceptCountTotal,
                                   refineCountTotal,
                                   hitMaxDepthCountTotal,
                                   acceptedDepthSumTotal,
                                   fusedEvalTotal_s,
                                   fusedAccumTotal_s);
                outCellP = outCellLocal;
                outWP = outWLocal;
                particleNnz[p] = static_cast<int>(outCellP.size());
            }
        }

#ifdef _OPENMP
        if (numParticles >= 64u)
        {
            const int activeThreads = static_cast<int>(ws.reducedCellByThread.size());
            if (activeThreads <= 0)
                throw std::runtime_error("[OrbitAveragedWeightsBuilder] invalid fused thread chunks");

            const auto& reducedCellByThread = ws.reducedCellByThread;
            const auto& reducedWByThread = ws.reducedWByThread;
            const auto& reducedOffsetsByThread = ws.reducedOffsetsByThread;
            const auto& reducedThreadBegin = ws.reducedThreadBegin;
            const auto& reducedThreadEnd = ws.reducedThreadEnd;

            std::vector<std::size_t> threadBase(static_cast<std::size_t>(activeThreads) + 1u, 0u);
            for (int t = 0; t < activeThreads; ++t)
            {
                threadBase[static_cast<std::size_t>(t) + 1u] =
                    threadBase[static_cast<std::size_t>(t)] + reducedCellByThread[static_cast<std::size_t>(t)].size();
            }
            outNnz = threadBase[static_cast<std::size_t>(activeThreads)];

            for (int t = 0; t < activeThreads; ++t)
            {
                const std::size_t pBeg = reducedThreadBegin[static_cast<std::size_t>(t)];
                const std::size_t pEnd = reducedThreadEnd[static_cast<std::size_t>(t)];
                const auto& localOffsets = reducedOffsetsByThread[static_cast<std::size_t>(t)];
                if (localOffsets.size() != (pEnd - pBeg + 1u))
                    throw std::runtime_error("[OrbitAveragedWeightsBuilder] inconsistent fused chunk offsets");

                const int base = static_cast<int>(threadBase[static_cast<std::size_t>(t)]);
                for (std::size_t i = 0; i < localOffsets.size(); ++i)
                    outOffsets[pBeg + i] = base + localOffsets[i];
            }

            outCell.resize(outNnz);
            outW.resize(outNnz, 0.0);

            std::exception_ptr firstError;
            #pragma omp parallel for num_threads(activeThreads) schedule(static)
            for (int t = 0; t < activeThreads; ++t)
            {
                try
                {
                    const std::size_t base = threadBase[static_cast<std::size_t>(t)];
                    const auto& cellFlat = reducedCellByThread[static_cast<std::size_t>(t)];
                    const auto& wFlat = reducedWByThread[static_cast<std::size_t>(t)];
                    if (cellFlat.size() != wFlat.size())
                        throw std::runtime_error("[OrbitAveragedWeightsBuilder] inconsistent fused chunk data");
                    if (!cellFlat.empty())
                    {
                        std::memcpy(outCell.data() + base, cellFlat.data(), cellFlat.size() * sizeof(int));
                        std::memcpy(outW.data() + base, wFlat.data(), wFlat.size() * sizeof(double));
                    }
                }
                catch (...)
                {
                    #pragma omp critical(lamp_orbitavg_fused_copy_error)
                    {
                        if (!firstError)
                            firstError = std::current_exception();
                    }
                }
            }
            if (firstError)
                std::rethrow_exception(firstError);
        }
        else
#endif
        {
            const auto& reducedCellByParticle = ws.reducedCellByParticle;
            const auto& reducedWByParticle = ws.reducedWByParticle;

            for (std::size_t p = 0; p < numParticles; ++p)
                outOffsets[p + 1u] = outOffsets[p] + particleNnz[p];

            outNnz = static_cast<std::size_t>(outOffsets.back());
            outCell.resize(outNnz);
            outW.resize(outNnz, 0.0);
            for (std::size_t p = 0; p < numParticles; ++p)
            {
                const auto& localCell = reducedCellByParticle[p];
                const auto& localW = reducedWByParticle[p];
                if (localCell.size() != static_cast<std::size_t>(particleNnz[p]) || localCell.size() != localW.size())
                    throw std::runtime_error("[OrbitAveragedWeightsBuilder] fused path CSR count mismatch");
                const std::size_t writeBeg = static_cast<std::size_t>(outOffsets[p]);
                const std::size_t nLocal = localCell.size();
                if (nLocal > 0u)
                {
                    std::memcpy(outCell.data() + writeBeg, localCell.data(), nLocal * sizeof(int));
                    std::memcpy(outW.data() + writeBeg, localW.data(), nLocal * sizeof(double));
                }
            }
        }

        if (timing)
        {
            timing->fusedEval_s += fusedEvalTotal_s;
            timing->fusedAccum_s += fusedAccumTotal_s;
            timing->sampleCount += sampleCountTotal;
            timing->sampleCacheFallbackCount += sampleFallbackTotal;
            timing->sampleCacheNnz += sampleNnzTotal;
            timing->outNnz += outNnz;
            timing->adaptAcceptCount += acceptCountTotal;
            timing->adaptRefineCount += refineCountTotal;
            timing->adaptHitMaxDepthCount += hitMaxDepthCountTotal;
            timing->adaptAcceptedDepthSum += acceptedDepthSumTotal;
            timing->total_s += secondsBetween(tTotal0, Clock::now());
        }

        out.assignCSR(numParticles, std::move(outOffsets), std::move(outCell), std::move(outW));
        return true;
    }

    // ------------------------------------------------------------------
    // Build pseudo-particles for quadrature samples (particle-major order).
    // Keep samples for each particle contiguous so later reduction can be done
    // per-particle without a global sort.
    // ------------------------------------------------------------------
    auto& sampleOffsetsByParticle = ws.sampleOffsetsByParticle;
    sampleOffsetsByParticle.assign(numParticles + 1u, 0u);

    const double invSqrt3 = 1.0 / std::sqrt(3.0);

    for (std::size_t p = 0; p < numParticles; ++p)
    {
        const int s0 = segmentOffsets[p];
        const int s1 = segmentOffsets[p + 1u];
        if (s0 < 0 || s1 < s0)
            return false;
        for (int si = s0; si < s1; ++si)
        {
            if (static_cast<std::size_t>(si) >= segments.size())
                return false;
            const OrbitSegment& seg = segments[static_cast<std::size_t>(si)];
            const double a = seg.t0;
            const double b = seg.t1;
            const double dur = b - a;
            if (!(dur > 0.0))
                continue;
            sampleOffsetsByParticle[p + 1u] += 2u;
        }
    }

    for (std::size_t p = 1; p < sampleOffsetsByParticle.size(); ++p)
        sampleOffsetsByParticle[p] += sampleOffsetsByParticle[p - 1u];

    auto& samplePos = ws.samplePos;
    auto& sampleHostLabel = ws.sampleHostLabel;
    auto& sampleCoeff = ws.sampleCoeff;
    samplePos.resize(sampleOffsetsByParticle.back());
    sampleHostLabel.assign(sampleOffsetsByParticle.back(), -1);
    sampleCoeff.assign(sampleOffsetsByParticle.back(), 0.0);

    auto fillSamplesForParticle = [&](std::size_t p) {
        std::size_t write = sampleOffsetsByParticle[p];
        const int s0 = segmentOffsets[p];
        const int s1 = segmentOffsets[p + 1u];
        for (int si = s0; si < s1; ++si)
        {
            const OrbitSegment& seg = segments[static_cast<std::size_t>(si)];
            const double a = seg.t0;
            const double b = seg.t1;
            const double dur = b - a;
            if (!(dur > 0.0))
                continue;

            const Vector2D x0s = seg.x0;
            const Vector2D x1s = unwrapToNearestImage(mesh, x0s, seg.x1, periodicX, periodicY);

            const double mid = 0.5 * (a + b);
            const double half = 0.5 * dur;
            const double tA = mid - half * invSqrt3;
            const double tB = mid + half * invSqrt3;

            auto posAt = [&](double t) -> Vector2D {
                const double s = (dur > 0.0) ? ((t - a) / dur) : 0.0;
                const Vector2D xr = lerp(x0s, x1s, std::clamp(s, 0.0, 1.0));
                return wrapToDomain(mesh, xr, periodicX, periodicY, epsWrap);
            };

            const double coeff = half / dt; // (dur/2)/dt

            samplePos[write] = posAt(tA);
            sampleHostLabel[write] = seg.cellLabel;
            sampleCoeff[write] = coeff;
            ++write;

            samplePos[write] = posAt(tB);
            sampleHostLabel[write] = seg.cellLabel;
            sampleCoeff[write] = coeff;
            ++write;
        }
    };

    const auto tSampleFill0 = Clock::now();
#ifdef _OPENMP
    if (numParticles >= 64u)
    {
        std::exception_ptr firstError;
        #pragma omp parallel for schedule(guided, 32)
        for (int pi = 0; pi < static_cast<int>(numParticles); ++pi)
        {
            try
            {
                fillSamplesForParticle(static_cast<std::size_t>(pi));
            }
            catch (...)
            {
                #pragma omp critical(lamp_orbitavg_fill_error)
                {
                    if (!firstError)
                        firstError = std::current_exception();
                }
            }
        }
        if (firstError)
            std::rethrow_exception(firstError);
    }
    else
#endif
    {
        for (std::size_t p = 0; p < numParticles; ++p)
            fillSamplesForParticle(p);
    }
    if (timing)
        timing->sampleFill_s += secondsBetween(tSampleFill0, Clock::now());

    const std::size_t nSamples = samplePos.size();
    if (timing)
        timing->sampleCount += nSamples;
    if (nSamples == 0u)
    {
        // Fallback: no segments provided. Build a trivial cache from host cells.
        std::vector<int> offsets(numParticles + 1u, 0);
        std::vector<int> cellIdx;
        std::vector<double> w;
        cellIdx.reserve(numParticles);
        w.reserve(numParticles);
        for (std::size_t p = 0; p < numParticles; ++p)
        {
            int cellLabel = hostCellLabelByParticle[p];
            int ci = -1;
            if (cellLabel >= 0)
                ci = mesh.cells(cellLabel).index();
            if (ci < 0)
                ci = 0;
            cellIdx.push_back(ci);
            w.push_back(1.0);
            offsets[p + 1u] = static_cast<int>(cellIdx.size());
        }
        if (timing)
            timing->outNnz += cellIdx.size();
        out.assignCSR(numParticles, std::move(offsets), std::move(cellIdx), std::move(w));
        return true;
    }

    CouplingWeightsCache& sampleCache = ws.sampleCache;
    CouplingWeightsCache::BuildOptions wopt = buildOpt;
    wopt.hostCellLabels = &sampleHostLabel;
    wopt.fallbackToHostOnZeroWeight = false;
    const auto tSampleCacheBuild0 = Clock::now();
    sampleCache.build(mesh, samplePos, hp, wopt);
    if (timing)
    {
        timing->sampleCacheBuild_s += secondsBetween(tSampleCacheBuild0, Clock::now());
        timing->sampleCacheFallbackCount += sampleCache.lastBuildFallbackCount();
        timing->sampleCacheNnz += sampleCache.nnz();
    }

    const auto& off = sampleCache.offsets();
    const auto& ci = sampleCache.cellIndices();
    const auto& ww = sampleCache.weights();

    if (off.size() != nSamples + 1u)
        return false;

    auto reduceParticle = [&](std::size_t p,
                              std::vector<int>& outCellLocal,
                              std::vector<double>& outWLocal,
                              std::vector<std::pair<int, double>>& scratchTrips) {
        const std::size_t sBeg = sampleOffsetsByParticle[p];
        const std::size_t sEnd = sampleOffsetsByParticle[p + 1u];

        scratchTrips.clear();

        for (std::size_t s = sBeg; s < sEnd; ++s)
        {
            const double c = sampleCoeff[s];
            const int e0 = off[s];
            const int e1 = off[s + 1u];
            if (e0 < 0 || e1 < e0)
                throw std::runtime_error("[OrbitAveragedWeightsBuilder] invalid sample CSR offsets");
            for (int e = e0; e < e1; ++e)
            {
                const std::size_t ee = static_cast<std::size_t>(e);
                if (ee >= ci.size() || ee >= ww.size())
                    throw std::runtime_error("[OrbitAveragedWeightsBuilder] sample CSR nnz out of range");
                const int cellIndex = ci[ee];
                const double wv = c * ww[ee];
                if (wv == 0.0)
                    continue;
                scratchTrips.emplace_back(cellIndex, wv);
            }
        }

        if (scratchTrips.empty())
        {
            outCellLocal.clear();
            outWLocal.clear();
        }
        else
        {
            std::sort(scratchTrips.begin(), scratchTrips.end(),
                      [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
                          return a.first < b.first;
                      });

            outCellLocal.clear();
            outWLocal.clear();
            outCellLocal.reserve(scratchTrips.size());
            outWLocal.reserve(scratchTrips.size());

            int currCell = scratchTrips.front().first;
            double currW = scratchTrips.front().second;
            for (std::size_t k = 1; k < scratchTrips.size(); ++k)
            {
                const int cell = scratchTrips[k].first;
                const double wv = scratchTrips[k].second;
                if (cell == currCell)
                {
                    currW += wv;
                }
                else
                {
                    if (currW != 0.0)
                    {
                        outCellLocal.push_back(currCell);
                        outWLocal.push_back(currW);
                    }
                    currCell = cell;
                    currW = wv;
                }
            }
            if (currW != 0.0)
            {
                outCellLocal.push_back(currCell);
                outWLocal.push_back(currW);
            }
        }

        if (opt.renormalize)
        {
            double sW = 0.0;
            for (double wv : outWLocal)
                sW += wv;
            if (sW > 0.0)
            {
                const double inv = 1.0 / sW;
                for (double& wv : outWLocal)
                    wv *= inv;
            }
            else
            {
                int cellLabel = hostCellLabelByParticle[p];
                int idxFallback = 0;
                if (cellLabel >= 0)
                {
                    const int ciFallback = mesh.cells(cellLabel).index();
                    if (ciFallback >= 0)
                        idxFallback = ciFallback;
                }
                outCellLocal.push_back(idxFallback);
                outWLocal.push_back(1.0);
            }
        }
    };

    auto& particleNnz = ws.particleNnz;
    particleNnz.assign(numParticles, 0);
    auto& outOffsets = ws.outOffsets;
    outOffsets.assign(numParticles + 1u, 0);
    std::vector<int> outCell;
    std::vector<double> outW;
    std::size_t outNnz = 0u;

    const auto tLocalReduce0 = Clock::now();
#ifdef _OPENMP
    if (numParticles >= 64u)
    {
        const int maxThreads = std::max(1, omp_get_max_threads());
        const int activeThreads = std::max(1, std::min<int>(maxThreads, static_cast<int>(numParticles)));

        auto& reducedCellByThread = ws.reducedCellByThread;
        auto& reducedWByThread = ws.reducedWByThread;
        auto& reducedOffsetsByThread = ws.reducedOffsetsByThread;
        auto& reducedThreadBegin = ws.reducedThreadBegin;
        auto& reducedThreadEnd = ws.reducedThreadEnd;
        auto& reducedChunkNnzHint = ws.reducedChunkNnzHint;
        reducedCellByThread.resize(activeThreads);
        reducedWByThread.resize(activeThreads);
        reducedOffsetsByThread.resize(activeThreads);
        reducedThreadBegin.resize(activeThreads);
        reducedThreadEnd.resize(activeThreads);
        if (reducedChunkNnzHint.size() != static_cast<std::size_t>(activeThreads))
            reducedChunkNnzHint.assign(static_cast<std::size_t>(activeThreads), 0u);

        std::exception_ptr firstError;
        // Schedule every logical chunk even if the runtime forms a smaller team.
        #pragma omp parallel for num_threads(activeThreads) schedule(static)
        for (int tid = 0; tid < activeThreads; ++tid)
        {
            const std::size_t pBeg = (numParticles * static_cast<std::size_t>(tid)) / static_cast<std::size_t>(activeThreads);
            const std::size_t pEnd = (numParticles * static_cast<std::size_t>(tid + 1)) / static_cast<std::size_t>(activeThreads);
            reducedThreadBegin[tid] = pBeg;
            reducedThreadEnd[tid] = pEnd;

            auto& cellFlat = reducedCellByThread[tid];
            auto& wFlat = reducedWByThread[tid];
            auto& localOffsets = reducedOffsetsByThread[tid];
            cellFlat.clear();
            wFlat.clear();
            const std::size_t chunkHint = reducedChunkNnzHint[static_cast<std::size_t>(tid)];
            if (chunkHint > 0u)
            {
                cellFlat.reserve(chunkHint);
                wFlat.reserve(chunkHint);
            }

            localOffsets.clear();
            localOffsets.reserve((pEnd - pBeg) + 1u);
            localOffsets.push_back(0);

            std::vector<std::pair<int, double>> scratchTrips;
            scratchTrips.reserve(128);
            std::vector<int> localCell;
            std::vector<double> localW;
            localCell.reserve(32);
            localW.reserve(32);

            try
            {
                for (std::size_t p = pBeg; p < pEnd; ++p)
                {
                    reduceParticle(p, localCell, localW, scratchTrips);
                    particleNnz[p] = static_cast<int>(localCell.size());
                    if (!localCell.empty())
                    {
                        const std::size_t oldN = cellFlat.size();
                        const std::size_t nLocal = localCell.size();
                        cellFlat.resize(oldN + nLocal);
                        wFlat.resize(oldN + nLocal);
                        std::memcpy(cellFlat.data() + oldN, localCell.data(), nLocal * sizeof(int));
                        std::memcpy(wFlat.data() + oldN, localW.data(), nLocal * sizeof(double));
                    }
                    localOffsets.push_back(static_cast<int>(cellFlat.size()));
                }
                reducedChunkNnzHint[static_cast<std::size_t>(tid)] = std::max<std::size_t>(cellFlat.size(), 32u);
            }
            catch (...)
            {
                #pragma omp critical(lamp_orbitavg_reduce_chunk_error)
                {
                    if (!firstError)
                        firstError = std::current_exception();
                }
            }
        }
        if (firstError)
            std::rethrow_exception(firstError);
    }
    else
#endif
    {
        auto& reducedCellByParticle = ws.reducedCellByParticle;
        auto& reducedWByParticle = ws.reducedWByParticle;
        reducedCellByParticle.resize(numParticles);
        reducedWByParticle.resize(numParticles);

        std::vector<std::pair<int, double>> scratchTrips;
        scratchTrips.reserve(128);
        for (std::size_t p = 0; p < numParticles; ++p)
        {
            auto& localCell = reducedCellByParticle[p];
            auto& localW = reducedWByParticle[p];
            reduceParticle(p, localCell, localW, scratchTrips);
            particleNnz[p] = static_cast<int>(localCell.size());
        }
    }
    if (timing)
        timing->localReduce_s += secondsBetween(tLocalReduce0, Clock::now());

    const auto tConcat0 = Clock::now();
#ifdef _OPENMP
    if (numParticles >= 64u)
    {
        const int activeThreads = static_cast<int>(ws.reducedCellByThread.size());
        if (activeThreads <= 0)
            throw std::runtime_error("[OrbitAveragedWeightsBuilder] invalid active thread chunks");
        const auto& reducedCellByThread = ws.reducedCellByThread;
        const auto& reducedWByThread = ws.reducedWByThread;
        const auto& reducedOffsetsByThread = ws.reducedOffsetsByThread;
        const auto& reducedThreadBegin = ws.reducedThreadBegin;
        const auto& reducedThreadEnd = ws.reducedThreadEnd;

        std::vector<std::size_t> threadBase(static_cast<std::size_t>(activeThreads) + 1u, 0u);
        for (int t = 0; t < activeThreads; ++t)
        {
            threadBase[static_cast<std::size_t>(t) + 1u] =
                threadBase[static_cast<std::size_t>(t)] + reducedCellByThread[t].size();
        }
        outNnz = threadBase[static_cast<std::size_t>(activeThreads)];

        for (int t = 0; t < activeThreads; ++t)
        {
            const std::size_t pBeg = reducedThreadBegin[t];
            const std::size_t pEnd = reducedThreadEnd[t];
            const auto& localOffsets = reducedOffsetsByThread[t];
            if (localOffsets.size() != (pEnd - pBeg + 1u))
                throw std::runtime_error("[OrbitAveragedWeightsBuilder] inconsistent chunk local offsets");
            const int base = static_cast<int>(threadBase[static_cast<std::size_t>(t)]);
            for (std::size_t i = 0; i < localOffsets.size(); ++i)
                outOffsets[pBeg + i] = base + localOffsets[i];
        }

        outCell.resize(outNnz);
        outW.resize(outNnz, 0.0);

        std::exception_ptr firstError;
        #pragma omp parallel for num_threads(activeThreads) schedule(static)
        for (int t = 0; t < activeThreads; ++t)
        {
            try
            {
                const std::size_t base = threadBase[static_cast<std::size_t>(t)];
                const auto& cellFlat = reducedCellByThread[t];
                const auto& wFlat = reducedWByThread[t];
                if (cellFlat.size() != wFlat.size())
                    throw std::runtime_error("[OrbitAveragedWeightsBuilder] inconsistent chunk flat data");
                if (!cellFlat.empty())
                {
                    std::memcpy(outCell.data() + base, cellFlat.data(), cellFlat.size() * sizeof(int));
                    std::memcpy(outW.data() + base, wFlat.data(), wFlat.size() * sizeof(double));
                }
            }
            catch (...)
            {
                #pragma omp critical(lamp_orbitavg_chunk_copy_error)
                {
                    if (!firstError)
                        firstError = std::current_exception();
                }
            }
        }
        if (firstError)
            std::rethrow_exception(firstError);
    }
    else
#endif
    {
        auto& reducedCellByParticle = ws.reducedCellByParticle;
        auto& reducedWByParticle = ws.reducedWByParticle;
        for (std::size_t p = 0; p < numParticles; ++p)
            outOffsets[p + 1u] = outOffsets[p] + particleNnz[p];

        outNnz = static_cast<std::size_t>(outOffsets.back());
        outCell.resize(outNnz);
        outW.resize(outNnz, 0.0);

        for (std::size_t p = 0; p < numParticles; ++p)
        {
            const auto& localCell = reducedCellByParticle[p];
            const auto& localW = reducedWByParticle[p];
            const std::size_t writeBeg = static_cast<std::size_t>(outOffsets[p]);
            if (localCell.size() != static_cast<std::size_t>(particleNnz[p]))
                throw std::runtime_error("[OrbitAveragedWeightsBuilder] inconsistent CSR count in fill pass");
            const std::size_t nLocal = localCell.size();
            if (nLocal > 0u)
            {
                std::memcpy(outCell.data() + writeBeg, localCell.data(), nLocal * sizeof(int));
                std::memcpy(outW.data() + writeBeg, localW.data(), nLocal * sizeof(double));
            }
        }
    }

    if (timing)
    {
        timing->concat_s += secondsBetween(tConcat0, Clock::now());
        timing->outNnz += outNnz;
        timing->total_s += secondsBetween(tTotal0, Clock::now());
    }

    out.assignCSR(numParticles, std::move(outOffsets), std::move(outCell), std::move(outW));
    return true;
}

bool buildOrbitAveragedCacheGauss2(const BackgroundMesh& mesh,
                                  const std::vector<double>& hp,
                                  const CouplingWeightsCache::BuildOptions& buildOpt,
                                  const std::size_t numParticles,
                                  const std::vector<int>& segmentOffsets,
                                  const std::vector<OrbitSegment>& segments,
                                  const std::vector<int>& hostCellLabelByParticle,
                                  const double dt,
                                  CouplingWeightsCache& out,
                                  const OrbitAveragedWeightsBuilderOptions& opt,
                                  OrbitAveragedWeightsBuilderTiming* timing)
{
    return buildOrbitAveragedCacheGauss2(mesh,
                                         hp,
                                         buildOpt,
                                         numParticles,
                                         segmentOffsets,
                                         segments,
                                         hostCellLabelByParticle,
                                         dt,
                                         out,
                                         /*workspace=*/nullptr,
                                         opt,
                                         timing);
}

} // namespace Orbit
} // namespace ImplicitPIC
