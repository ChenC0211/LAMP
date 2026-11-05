#include "Modules/Plasma/ExplicitPIC/Deposit/CouplingWeightsCache.h"

#include "Modules/Plasma/ExplicitPIC/Deposit/CICAnalyticIntegral.h"
#include "Modules/Plasma/ExplicitPIC/Deposit/CICMomentExactCut.h"
#include "Modules/Plasma/ExplicitPIC/Deposit/ParticleStencilQuery.h"

#include "Core/Base/Base.h"
#include "Core/Mesh/BackgroundMesh.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <limits>
#include <stdexcept>
#include "Core/IO/ConsoleControl.h"

#include <fstream>

#if defined(_OPENMP)
#include <omp.h>
#endif

namespace {

// Simple Kahan summation for doubles.
inline double kahanSum(const std::vector<double>& v)
{
    double sum = 0.0;
    double c = 0.0;
    for (double x : v)
    {
        const double y = x - c;
        const double t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
    return sum;
}

inline double clamp01(const double a)
{
    return std::max(0.0, std::min(1.0, a));
}

// Coarse–fine composite band predicate.
// A leaf cell is considered in the coarse–fine neighborhood set B if it touches
// a coarse–fine interface across any face. We detect this using the neighbor
// connectivity:
//   - two neighbors in a direction (neiCell[dir][1] != -1), or
//   - a neighbor at a different refinement level, or
//   - a neighbor that is not a leaf.
inline bool isCompositeBandCell(const BackgroundMesh& mesh, const Cell& c)
{
    // We only care about leaf cells; for non-leaf we conservatively treat as true.
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
            // If neighbor is not leaf or at different level, this is a coarse–fine interface.
            const bool nbLeaf = (nb.son[0] == -1 && nb.son[1] == -1 && nb.son[2] == -1 && nb.son[3] == -1);
            if (!nbLeaf || nb.level() != lvl)
                return true;
        }
    }
    return false;
}

inline double minCellSize(const Cell& c)
{
    const Vector2D s = c.size();
    return std::min(s.x, s.y);
}

inline bool isLeafCellSimple(const Cell& c)
{
    return (c.son[0] == -1 && c.son[1] == -1 && c.son[2] == -1 && c.son[3] == -1);
}

inline bool nearlyEqualRel(const double a, const double b, const double relTol = 1e-12)
{
    const double s = std::max({1.0, std::abs(a), std::abs(b)});
    return std::abs(a - b) <= relTol * s;
}

// Stencils are typically small (< O(10)), so a stable insertion sort avoids
// per-particle temporary allocations while keeping deterministic ordering.
inline void sortCellWeightsByCellIndex(std::vector<int>& cellIndex,
                                       std::vector<double>& weights)
{
    const std::size_t n = cellIndex.size();
    if (n <= 1 || weights.size() != n)
        return;

    for (std::size_t i = 1; i < n; ++i)
    {
        const int keyIdx = cellIndex[i];
        const double keyW = weights[i];
        std::size_t j = i;
        while (j > 0 && cellIndex[j - 1] > keyIdx)
        {
            cellIndex[j] = cellIndex[j - 1];
            weights[j] = weights[j - 1];
            --j;
        }
        cellIndex[j] = keyIdx;
        weights[j] = keyW;
    }
}

} // namespace

void CouplingWeightsCache::clear()
{
    numParticles_ = 0;
    offsets_.clear();
    cellIndices_.clear();
    weights_.clear();
    rawWeights_.clear();
    effectiveHp_.clear();
    inCompositeBand_.clear();
    lastBuildFallbackCount_ = 0;
    lastBuildParticleCount_ = 0;
}

void CouplingWeightsCache::assignCSR(std::size_t numParticles,
                                    std::vector<int> offsets,
                                    std::vector<int> cellIndices,
                                    std::vector<double> weights)
{
    // Validate offsets.
    if (offsets.size() != (numParticles + 1u))
        throw std::runtime_error("CouplingWeightsCache::assignCSR: offsets.size() must be numParticles+1");

    if (offsets.empty() || offsets.front() != 0)
        throw std::runtime_error("CouplingWeightsCache::assignCSR: offsets[0] must be 0");

    for (std::size_t i = 1; i < offsets.size(); ++i)
    {
        if (offsets[i] < offsets[i - 1])
            throw std::runtime_error("CouplingWeightsCache::assignCSR: offsets must be non-decreasing");
    }

    const std::size_t nnz = static_cast<std::size_t>(offsets.back());
    if (cellIndices.size() != nnz || weights.size() != nnz)
        throw std::runtime_error("CouplingWeightsCache::assignCSR: cellIndices/weights size mismatch with offsets.back()");

    // Reset auxiliary state that is only meaningful for build().
    numParticles_ = numParticles;
    offsets_ = std::move(offsets);
    cellIndices_ = std::move(cellIndices);
    weights_ = std::move(weights);

    // Raw weights are not supplied via this path.
    rawWeights_.clear();

    // hp stats and other per-particle diagnostics are unknown.
    hasHpStats_ = false;
    lastHpStats_ = HpStats{};
    effectiveHp_.clear();
    inCompositeBand_.clear();
    lastBuildFallbackCount_ = 0;
    lastBuildParticleCount_ = numParticles_;
}

void CouplingWeightsCache::buildFromCSR(std::size_t numParticles,
                                       std::vector<int> offsets,
                                       std::vector<int> cellIndices,
                                       std::vector<double> weights)
{
    assignCSR(numParticles, std::move(offsets), std::move(cellIndices), std::move(weights));
}

void CouplingWeightsCache::reserve(std::size_t numParticles, std::size_t nnzGuess)
{
    // CSR offsets has length (numParticles + 1).
    offsets_.reserve(numParticles + 1);
    cellIndices_.reserve(nnzGuess);
    weights_.reserve(nnzGuess);
    if (keepRawWeights_)
        rawWeights_.reserve(nnzGuess);
    effectiveHp_.reserve(numParticles);
    inCompositeBand_.reserve(numParticles);
}

void CouplingWeightsCache::enableRawWeights(const bool enabled)
{
    keepRawWeights_ = enabled;
    if (!keepRawWeights_)
    {
        // Keep memory simple and deterministic: drop raw weights when disabled.
        rawWeights_.clear();
        rawWeights_.shrink_to_fit();
    }
}

void CouplingWeightsCache::buildPlaceholder(const std::size_t numParticles)
{
    // Placeholder: no coupling stencil yet.
    clear();
    numParticles_ = numParticles;
    offsets_.assign(numParticles_ + 1, 0);
    lastBuildFallbackCount_ = 0;
    lastBuildParticleCount_ = numParticles_;
}

void CouplingWeightsCache::build(const BackgroundMesh& mesh,
                               const std::vector<Vector2D>& positions,
                               const std::vector<double>& hp,
                               const BuildOptions& opt)
{
    clear();
    hasHpStats_ = false;
    lastHpStats_ = HpStats{};

    numParticles_ = positions.size();
    lastBuildParticleCount_ = numParticles_;
    offsets_.assign(numParticles_ + 1, 0);
    if (opt.recordPerParticleHp)
        effectiveHp_.assign(numParticles_, 0.0);
    else
        effectiveHp_.clear();
    if (opt.recordPerParticleCompositeBand)
        inCompositeBand_.assign(numParticles_, 0);
    else
        inCompositeBand_.clear();

    if (numParticles_ == 0)
        return;

    const bool hpIsConst = (hp.size() == 1);
    const double hpConst = hpIsConst ? hp[0] : 0.0;

    // Debug/dump controls
    const bool wantDump = (opt.debugParticleId != std::numeric_limits<std::size_t>::max())
        && !opt.dumpWeightsPath.empty();
    const bool willDumpThisBuild = wantDump && (!opt.dumpWeightsOnce || !dumpedWeightsOnce_);
    const bool storeRawThisBuild = keepRawWeights_ || willDumpThisBuild;

    // If debugPrintStencil is enabled, the stencil query prints to stdout.
    // Disable OpenMP in this mode to avoid interleaved output.
    const bool allowParallel = opt.parallelBuild && !opt.debugPrintStencil;

#if defined(_OPENMP)
    const bool useParallel = allowParallel && (static_cast<int>(numParticles_) >= opt.ompMinParticles);
#else
    const bool useParallel = false;
#endif

    auto fillWeightsForStencil = [&](double xp,
                                     double yp,
                                     double hp_eff,
                                     bool inBand,
                                     const std::vector<int>& stencilLabels,
                                     std::vector<int>& outCellIndex,
                                     std::vector<double>& outWraw)
    {
        outCellIndex.clear();
        outWraw.clear();

        for (int cellLabel : stencilLabels)
        {
            const Cell& cell = mesh.cells(cellLabel);
            const int cellIndex = cell.index();

            bool isCut = false;
            if (opt.levelSetI >= 0 && opt.levelSetI < static_cast<int>(cell.isCut.size()))
                isCut = (cell.isCut[opt.levelSetI] != 0);
            else if (opt.cutByAlpha && opt.alphaByCellIndex &&
                     cellIndex < static_cast<int>(opt.alphaByCellIndex->size()))
            {
                double alphaCut = (*opt.alphaByCellIndex)[cellIndex];
                if (!std::isfinite(alphaCut))
                    alphaCut = 1.0;
                if (alphaCut < 0.0)
                    alphaCut = 0.0;
                if (alphaCut > 1.0)
                    alphaCut = 1.0;
                isCut = (alphaCut > opt.cutAlphaEps) && (alphaCut < (1.0 - opt.cutAlphaEps));
            }

            const Vector2D cc = cell.centerCoord();
            const Vector2D sz = cell.size();
            const double xL = cc.x - 0.5 * sz.x;
            const double xR = cc.x + 0.5 * sz.x;
            const double yL = cc.y - 0.5 * sz.y;
            const double yR = cc.y + 0.5 * sz.y;

            const double hK = minCellSize(cell);
            const bool doVirtualPartition =
                opt.enableCompositeBand && inBand && (hK > hp_eff * (1.0 + 1e-12));

            // Uncut overlap weight (with optional virtual partition).
            double wUncut = CICAnalytic::rectWeight2D(xL, xR, yL, yR, xp, yp, hp_eff);
            if (doVirtualPartition)
            {
                const double xM = 0.5 * (xL + xR);
                const double yM = 0.5 * (yL + yR);
                const double w00 = CICAnalytic::rectWeight2D(xL, xM, yL, yM, xp, yp, hp_eff);
                const double w10 = CICAnalytic::rectWeight2D(xM, xR, yL, yM, xp, yp, hp_eff);
                const double w01 = CICAnalytic::rectWeight2D(xL, xM, yM, yR, xp, yp, hp_eff);
                const double w11 = CICAnalytic::rectWeight2D(xM, xR, yM, yR, xp, yp, hp_eff);
                wUncut = w00 + w10 + w01 + w11;
            }

            double w = wUncut;
            if (!(w > 0.0))
                continue;

            double alpha = 1.0;
            if (opt.alphaByCellIndex && cellIndex < static_cast<int>(opt.alphaByCellIndex->size()))
                alpha = (*opt.alphaByCellIndex)[cellIndex];
            alpha = clamp01(alpha);
            if (!(alpha > 0.0))
                continue;

            if (isCut)
            {
                bool usedMomentExact = false;

                if (opt.enableMomentExactCut && opt.levelSetVert)
                {
                    if (doVirtualPartition)
                    {
                        CICMomentExactCut::Plane pl;
                        if (CICMomentExactCut::buildCellPlaneFromVertexLevelSet(mesh, cell, *opt.levelSetVert, &pl))
                        {
                            const double xM = 0.5 * (xL + xR);
                            const double yM = 0.5 * (yL + yR);
                            double w00 = 0.0, w10 = 0.0, w01 = 0.0, w11 = 0.0;
                            const bool ok00 = CICMomentExactCut::cutRectWeight2DFromPlane(xL, xM, yL, yM, pl, xp, yp, hp_eff, opt.insideIsNegative, &w00);
                            const bool ok10 = CICMomentExactCut::cutRectWeight2DFromPlane(xM, xR, yL, yM, pl, xp, yp, hp_eff, opt.insideIsNegative, &w10);
                            const bool ok01 = CICMomentExactCut::cutRectWeight2DFromPlane(xL, xM, yM, yR, pl, xp, yp, hp_eff, opt.insideIsNegative, &w01);
                            const bool ok11 = CICMomentExactCut::cutRectWeight2DFromPlane(xM, xR, yM, yR, pl, xp, yp, hp_eff, opt.insideIsNegative, &w11);
                            if (ok00 && ok10 && ok01 && ok11 &&
                                std::isfinite(w00) && std::isfinite(w10) && std::isfinite(w01) && std::isfinite(w11))
                            {
                                const double wCut = std::max(0.0, w00) + std::max(0.0, w10) +
                                                    std::max(0.0, w01) + std::max(0.0, w11);
                                w = std::min(wUncut, wCut);
                                usedMomentExact = true;
                            }
                        }
                    }
                    else
                    {
                        double wCut = 0.0;
                        if (CICMomentExactCut::cutCellWeight2D(mesh,
                                                              cellLabel,
                                                              *opt.levelSetVert,
                                                              xp,
                                                              yp,
                                                              hp_eff,
                                                              opt.insideIsNegative,
                                                              &wCut))
                        {
                            if (std::isfinite(wCut) && wCut >= 0.0)
                            {
                                w = std::min(wUncut, std::max(0.0, wCut));
                                usedMomentExact = true;
                            }
                        }
                    }
                }

                if (!usedMomentExact)
                {
                    w = wUncut * alpha;
                }
            }
            else if (alpha < 1.0)
            {
                // Cut flag may miss partially-filled cells; scale by alpha as a safety.
                w = wUncut * alpha;
            }

            if (w > 0.0)
            {
                outCellIndex.push_back(cellIndex);
                outWraw.push_back(w);
            }
        }
    };

    // Domain extents (used by periodic ghost images)
    const double xMin = mesh.minCoord.x;
    const double xMax = mesh.maxCoord.x;
    const double yMin = mesh.minCoord.y;
    const double yMax = mesh.maxCoord.y;
    const double Lx = xMax - xMin;
    const double Ly = yMax - yMin;
    const double eps = 1e-12 * std::max(1.0, std::max(std::abs(Lx), std::abs(Ly)));

    auto clampToDomain = [&](double x, const double a, const double b) -> double
    {
        // Clamp into the open interval (a,b) to avoid boundary ambiguity.
        return std::min(std::max(x, a + eps), b - eps);
    };

    auto appendCellWeights = [&](std::vector<int>& idx,
                                 std::vector<double>& w,
                                 const std::vector<int>& addIdx,
                                 const std::vector<double>& addW)
    {
        if (addIdx.empty())
            return;
        const std::size_t oldN = idx.size();
        const std::size_t addN = addIdx.size();
        idx.resize(oldN + addN);
        w.resize(oldN + addN);
        std::memcpy(idx.data() + oldN, addIdx.data(), addN * sizeof(int));
        std::memcpy(w.data() + oldN, addW.data(), addN * sizeof(double));
    };

    auto mergeSortedCellWeightsInPlace = [&](std::vector<int>& idx, std::vector<double>& w)
    {
        if (idx.empty())
            return;
        sortCellWeightsByCellIndex(idx, w);
        std::size_t out = 0u;
        int currCell = idx[0];
        double currW = w[0];
        for (std::size_t i = 1; i < idx.size(); ++i)
        {
            if (idx[i] == currCell)
            {
                currW += w[i];
            }
            else
            {
                idx[out] = currCell;
                w[out] = currW;
                ++out;
                currCell = idx[i];
                currW = w[i];
            }
        }
        idx[out] = currCell;
        w[out] = currW;
        ++out;
        idx.resize(out);
        w.resize(out);
    };

    // Fast stencil path for common non-EB / non-periodic / small-hp cases:
    // host + face neighbors + one extra ring, then AABB-filter candidates.
    auto tryBuildFastHostRingStencil = [&](const double xp,
                                           const double yp,
                                           const double hp_eff,
                                           const int hostCellLabel,
                                           std::vector<int>* stencilOut) -> bool
    {
        if (!stencilOut || hostCellLabel < 0)
            return false;
        if (opt.periodicX || opt.periodicY)
            return false;
        if (opt.levelSetVert != nullptr || opt.cutByAlpha)
            return false;

        const Cell& host = mesh.cells(hostCellLabel);
        if (!isLeafCellSimple(host))
            return false;

        // Keep this path conservative near coarse-fine interfaces.
        if (opt.enableCompositeBand && isCompositeBandCell(mesh, host))
            return false;

        const double hHost = minCellSize(host);
        if (!(std::isfinite(hHost) && hHost > 0.0))
            return false;
        if (hp_eff > hHost * (1.0 + 1e-12))
            return false;

        std::vector<int>& out = *stencilOut;
        out.clear();
        out.reserve(24);

        const auto box = ParticleStencilQuery::supportAABB(xp, yp, hp_eff);
        auto tryAppend = [&](const int lbl)
        {
            if (lbl < 0)
                return;
            if (std::find(out.begin(), out.end(), lbl) != out.end())
                return;
            const Cell& c = mesh.cells(lbl);
            if (!isLeafCellSimple(c))
                return;
            if (ParticleStencilQuery::cellIntersectsAABB(c, box))
                out.push_back(lbl);
        };

        tryAppend(hostCellLabel);

        std::vector<int> faceNbrs;
        faceNbrs.reserve(8);
        for (int dir = 0; dir < 4; ++dir)
        {
            for (int slot = 0; slot < 2; ++slot)
            {
                const int nb = host.neiCell[dir][slot];
                if (nb != -1)
                {
                    faceNbrs.push_back(nb);
                    tryAppend(nb);
                }
            }
        }

        for (const int nb : faceNbrs)
        {
            if (nb < 0)
                continue;
            const Cell& c = mesh.cells(nb);
            if (!isLeafCellSimple(c))
                continue;
            for (int dir = 0; dir < 4; ++dir)
            {
                for (int slot = 0; slot < 2; ++slot)
                {
                    tryAppend(c.neiCell[dir][slot]);
                }
            }
        }

        if (opt.deterministicOrdering)
            std::sort(out.begin(), out.end());

        return !out.empty();
    };

    // Reuse cache for repeated single-sample queries (hot in orbit-avg evalAt).
    int cachedStencilHost = -1;
    double cachedStencilHp = 0.0;
    bool cachedStencilValid = false;
    std::vector<int> cachedStencil;
    cachedStencil.reserve(64);

    int fallbackCountTotal = 0;

    if (!useParallel)
    {
        // ------------------------------
        // Legacy serial build (kept for small particle counts and debug output)
        // ------------------------------
        // Conservative guess: CIC support on quadtree is typically < 16 cells.
        reserve(numParticles_, numParticles_ * 9);

        double hpMin = std::numeric_limits<double>::infinity();
        double hpMax = 0.0;
        double hpSum = 0.0;
        std::size_t hpCount = 0;

        // Build periodic image shifts for a particle position.
        // Each returned shift (sx,sy) will be added to (xp,yp) and used as the
        // CIC kernel center. The shifted center may be outside the domain;
        // stencil enumeration will seed from a clamped point inside the domain.
        auto buildPeriodicImageShifts = [&](const double xp, const double yp, const double hp_eff) -> std::vector<Vector2D>
        {
            std::vector<double> sxList{0.0};
            std::vector<double> syList{0.0};

            if (opt.periodicX && Lx > 0.0)
            {
                if (xp + hp_eff > xMax)
                    sxList.push_back(-Lx);
                if (xp - hp_eff < xMin)
                    sxList.push_back(+Lx);
            }
            if (opt.periodicY && Ly > 0.0)
            {
                if (yp + hp_eff > yMax)
                    syList.push_back(-Ly);
                if (yp - hp_eff < yMin)
                    syList.push_back(+Ly);
            }

            std::vector<Vector2D> shifts;
            shifts.reserve(sxList.size() * syList.size());
            for (double sx : sxList)
                for (double sy : syList)
                    shifts.push_back(Vector2D(sx, sy));
            return shifts;
        };

        for (std::size_t p = 0; p < numParticles_; ++p)
        {
            offsets_[p] = static_cast<int>(cellIndices_.size());

            const double xp = positions[p].x;
            const double yp = positions[p].y;

            const double hp_in = hpIsConst ? hpConst : hp[p];
            double hp_eff = hp_in;

            int hostCellLabel = -1;
            if (opt.hostCellLabels && p < opt.hostCellLabels->size())
                hostCellLabel = (*opt.hostCellLabels)[p];
            if (hostCellLabel < 0)
                hostCellLabel = mesh.findCellLabel(Vector2D(xp, yp));
            if (hostCellLabel < 0)
            {
                // Outside domain: emit an empty stencil for this particle.
                offsets_[p + 1] = static_cast<int>(cellIndices_.size());
                continue;
            }

            const double hHost = minCellSize(mesh.cells(hostCellLabel));
            // Keep a positive user-provided hp; otherwise fall back to host-cell size.
            if (!(std::isfinite(hp_eff) && hp_eff > 0.0) &&
                std::isfinite(hHost) && hHost > 0.0)
                hp_eff = hHost;

            // Build base + periodic ghost stencils (if needed)
            const std::vector<Vector2D> shifts = buildPeriodicImageShifts(xp, yp, hp_eff);
            const bool singleShift = (shifts.size() == 1u);
            std::vector<std::vector<int>> stencils(shifts.size());

            // First pass stencil enumeration (used for hp shrink + composite band detection)
            double minHK = 1.0e300;
            double maxHK = 0.0;
            if (singleShift)
            {
                bool usedFastStencil = tryBuildFastHostRingStencil(xp, yp, hp_eff, hostCellLabel, &stencils[0]);
                if (!usedFastStencil)
                {
                    const bool canReuse = cachedStencilValid
                        && cachedStencilHost == hostCellLabel
                        && nearlyEqualRel(cachedStencilHp, hp_eff);
                    if (canReuse)
                    {
                        stencils[0] = cachedStencil;
                    }
                    else
                    {
                        const bool dbg = opt.debugPrintStencil && (p == opt.debugParticleId);
                        ParticleStencilQuery::enumerateLeafCellsIntersectParticleSupport(
                            mesh, xp, yp, hp_eff, &stencils[0], dbg, hostCellLabel, opt.deterministicOrdering);
                        cachedStencilHost = hostCellLabel;
                        cachedStencilHp = hp_eff;
                        cachedStencil = stencils[0];
                        cachedStencilValid = true;
                    }
                }
                else
                {
                    cachedStencilHost = hostCellLabel;
                    cachedStencilHp = hp_eff;
                    cachedStencil = stencils[0];
                    cachedStencilValid = true;
                }
                for (int cellLabel : stencils[0])
                {
                    const Cell& cell = mesh.cells(cellLabel);
                    const double hK = minCellSize(cell);
                    minHK = std::min(minHK, hK);
                    maxHK = std::max(maxHK, hK);
                }
            }
            else
            {
                for (std::size_t si = 0; si < shifts.size(); ++si)
                {
                    const double xg = xp + shifts[si].x;
                    const double yg = yp + shifts[si].y;

                    // Seed from a clamped in-domain point to ensure host cell exists.
                    const double xs = clampToDomain(xg, xMin, xMax);
                    const double ys = clampToDomain(yg, yMin, yMax);
                    const int hostSeed = mesh.findCellLabel(Vector2D(xs, ys));

                    const bool dbg = opt.debugPrintStencil && (p == opt.debugParticleId)
                        && (std::abs(shifts[si].x) < 0.5 * eps) && (std::abs(shifts[si].y) < 0.5 * eps);

                    ParticleStencilQuery::enumerateLeafCellsIntersectParticleSupport(
                        mesh, xg, yg, hp_eff, &stencils[si], dbg, hostSeed, opt.deterministicOrdering);

                    for (int cellLabel : stencils[si])
                    {
                        const Cell& cell = mesh.cells(cellLabel);
                        const double hK = minCellSize(cell);
                        minHK = std::min(minHK, hK);
                        maxHK = std::max(maxHK, hK);
                    }
                }
            }

            const bool hasLevelMismatch = (maxHK > minHK * (1.0 + 1e-12));
            bool doShrink = false;
            if (opt.enableCompositeBand && hasLevelMismatch && minHK < 1.0e299)
            {
                const double hp_star = minHK;
                if (hp_star < hp_eff)
                {
                    hp_eff = hp_star;
                    doShrink = true;
                }
            }

            if (doShrink)
            {
                // Re-enumerate stencils at reduced hp.
                if (singleShift)
                {
                    bool usedFastStencil = tryBuildFastHostRingStencil(xp, yp, hp_eff, hostCellLabel, &stencils[0]);
                    if (!usedFastStencil)
                    {
                        const bool canReuse = cachedStencilValid
                            && cachedStencilHost == hostCellLabel
                            && nearlyEqualRel(cachedStencilHp, hp_eff);
                        if (canReuse)
                        {
                            stencils[0] = cachedStencil;
                        }
                        else
                        {
                            const bool dbg = opt.debugPrintStencil && (p == opt.debugParticleId);
                            ParticleStencilQuery::enumerateLeafCellsIntersectParticleSupport(
                                mesh, xp, yp, hp_eff, &stencils[0], dbg, hostCellLabel, opt.deterministicOrdering);
                            cachedStencilHost = hostCellLabel;
                            cachedStencilHp = hp_eff;
                            cachedStencil = stencils[0];
                            cachedStencilValid = true;
                        }
                    }
                    else
                    {
                        cachedStencilHost = hostCellLabel;
                        cachedStencilHp = hp_eff;
                        cachedStencil = stencils[0];
                        cachedStencilValid = true;
                    }
                }
                else
                {
                    for (std::size_t si = 0; si < shifts.size(); ++si)
                    {
                        const double xg = xp + shifts[si].x;
                        const double yg = yp + shifts[si].y;
                        const double xs = clampToDomain(xg, xMin, xMax);
                        const double ys = clampToDomain(yg, yMin, yMax);
                        const int hostSeed = mesh.findCellLabel(Vector2D(xs, ys));
                        const bool dbg = opt.debugPrintStencil && (p == opt.debugParticleId)
                            && (std::abs(shifts[si].x) < 0.5 * eps) && (std::abs(shifts[si].y) < 0.5 * eps);
                        ParticleStencilQuery::enumerateLeafCellsIntersectParticleSupport(
                            mesh, xg, yg, hp_eff, &stencils[si], dbg, hostSeed, opt.deterministicOrdering);
                    }
                }
            }

            if (opt.recordHpStats && std::isfinite(hp_eff) && hp_eff > 0.0)
            {
                hpMin = std::min(hpMin, hp_eff);
                hpMax = std::max(hpMax, hp_eff);
                hpSum += hp_eff;
                hpCount += 1;
            }

            const bool inCompositeBand = opt.enableCompositeBand ? hasLevelMismatch : false;
            if (opt.recordPerParticleHp)
                effectiveHp_[p] = hp_eff;
            if (opt.recordPerParticleCompositeBand)
                inCompositeBand_[p] = static_cast<unsigned char>(inCompositeBand ? 1 : 0);

            // Build and merge raw weights across periodic images.
            std::vector<int> tmpCellIndex;
            std::vector<double> tmpWraw;
            tmpCellIndex.reserve(64);
            tmpWraw.reserve(64);

            std::vector<int> mergedCellIndex;
            std::vector<double> mergedWraw;
            mergedCellIndex.reserve(64);
            mergedWraw.reserve(64);

            if (singleShift)
            {
                fillWeightsForStencil(xp, yp, hp_eff, inCompositeBand, stencils[0], mergedCellIndex, mergedWraw);
            }
            else
            {
                for (std::size_t si = 0; si < shifts.size(); ++si)
                {
                    const double xg = xp + shifts[si].x;
                    const double yg = yp + shifts[si].y;
                    tmpCellIndex.clear();
                    tmpWraw.clear();
                    fillWeightsForStencil(xg, yg, hp_eff, inCompositeBand, stencils[si], tmpCellIndex, tmpWraw);
                    appendCellWeights(mergedCellIndex, mergedWraw, tmpCellIndex, tmpWraw);
                }

                // Deterministic ordering + duplicate merge by cellIndex.
                mergeSortedCellWeightsInPlace(mergedCellIndex, mergedWraw);
            }

            double sumW = kahanSum(mergedWraw);

            if (!(sumW > 0.0) || mergedCellIndex.empty())
            {
                if (opt.fallbackToHostOnZeroWeight)
                {
                    // Fallback to host cell.
                    const int hostIndex = mesh.cells(hostCellLabel).index();
                    cellIndices_.push_back(hostIndex);
                    weights_.push_back(1.0);
                    if (storeRawThisBuild)
                        rawWeights_.push_back(0.0);
                    fallbackCountTotal += 1;
                }
                offsets_[p + 1] = static_cast<int>(cellIndices_.size());
                continue;
            }

            // Renormalize and write entries
            int maxLocal = -1;
            double maxVal = -1.0;

            for (std::size_t k = 0; k < mergedCellIndex.size(); ++k)
            {
                const double wtilde = mergedWraw[k] / sumW;
                cellIndices_.push_back(mergedCellIndex[k]);
                weights_.push_back(wtilde);
                if (storeRawThisBuild)
                    rawWeights_.push_back(mergedWraw[k]);
                if (wtilde > maxVal)
                {
                    maxVal = wtilde;
                    maxLocal = static_cast<int>(k);
                }
            }

            // Residual correction to enforce sum=1 exactly.
            {
                const std::size_t a = static_cast<std::size_t>(offsets_[p]);
                const std::size_t b = weights_.size();
                double sumTilde = 0.0;
                for (std::size_t i = a; i < b; ++i)
                    sumTilde += weights_[i];
                const double r = 1.0 - sumTilde;
                if (maxLocal >= 0)
                    weights_[a + static_cast<std::size_t>(maxLocal)] += r;
            }

            offsets_[p + 1] = static_cast<int>(cellIndices_.size());
        }

        if (opt.recordHpStats && hpCount > 0)
        {
            lastHpStats_.min = hpMin;
            lastHpStats_.max = hpMax;
            lastHpStats_.mean = hpSum / static_cast<double>(hpCount);
            lastHpStats_.count = hpCount;
            hasHpStats_ = true;
        }
    }
    else
    {
        // ------------------------------
        // OpenMP parallel build (single-pass by thread-chunk):
        //   - each thread builds a contiguous particle chunk once
        //   - local flat CSR data are concatenated after join
        // This avoids rebuilding stencils/weights twice in the old two-pass path.
        // ------------------------------
        std::vector<int> hostCellLabel(numParticles_, -1);
        std::vector<double> hpEff(numParticles_, 0.0);
        std::vector<unsigned char> inCompositeBand(numParticles_, 0);

        // Helper: build periodic image shifts without allocating.
        auto buildShiftLists = [&](const double xp, const double yp, const double hp_eff,
                                   double* sx, int& nx, double* sy, int& ny)
        {
            nx = 1;
            ny = 1;
            sx[0] = 0.0;
            sy[0] = 0.0;

            if (opt.periodicX && Lx > 0.0)
            {
                if (xp + hp_eff > xMax)
                    sx[nx++] = -Lx;
                if (xp - hp_eff < xMin)
                    sx[nx++] = +Lx;
            }
            if (opt.periodicY && Ly > 0.0)
            {
                if (yp + hp_eff > yMax)
                    sy[ny++] = -Ly;
                if (yp - hp_eff < yMin)
                    sy[ny++] = +Ly;
            }
        };

#if defined(_OPENMP)
        const int maxThreads = std::max(1, omp_get_max_threads());
        const int activeThreads = std::max(1, std::min<int>(maxThreads, static_cast<int>(numParticles_)));
#else
        const int activeThreads = 1;
#endif
        const bool nonPeriodic = (!opt.periodicX && !opt.periodicY);

        std::vector<std::vector<int>> cellByThread(static_cast<std::size_t>(activeThreads));
        std::vector<std::vector<double>> weightByThread(static_cast<std::size_t>(activeThreads));
        std::vector<std::vector<double>> rawByThread;
        if (storeRawThisBuild)
            rawByThread.resize(static_cast<std::size_t>(activeThreads));
        std::vector<std::vector<int>> offsetsByThread(static_cast<std::size_t>(activeThreads));
        std::vector<std::size_t> threadBegin(static_cast<std::size_t>(activeThreads), 0u);
        std::vector<std::size_t> threadEnd(static_cast<std::size_t>(activeThreads), 0u);
        std::vector<int> fallbackByThread(static_cast<std::size_t>(activeThreads), 0);

        std::exception_ptr firstError;
        // Process every logical chunk even when a nested team is serialized
        // or the runtime supplies fewer workers than requested.
        #pragma omp parallel for num_threads(activeThreads) schedule(static)
        for (int tid = 0; tid < activeThreads; ++tid)
        {
            try
            {
                const std::size_t pBeg =
                    (numParticles_ * static_cast<std::size_t>(tid)) / static_cast<std::size_t>(activeThreads);
                const std::size_t pEnd =
                    (numParticles_ * static_cast<std::size_t>(tid + 1)) / static_cast<std::size_t>(activeThreads);
                threadBegin[static_cast<std::size_t>(tid)] = pBeg;
                threadEnd[static_cast<std::size_t>(tid)] = pEnd;

                auto& cellFlat = cellByThread[static_cast<std::size_t>(tid)];
                auto& wFlat = weightByThread[static_cast<std::size_t>(tid)];
                cellFlat.clear();
                wFlat.clear();
                const std::size_t roughNnz = (pEnd - pBeg) * 8u;
                cellFlat.reserve(roughNnz);
                wFlat.reserve(roughNnz);
                std::vector<double>* rawFlat = nullptr;
                if (storeRawThisBuild)
                {
                    rawFlat = &rawByThread[static_cast<std::size_t>(tid)];
                    rawFlat->clear();
                    rawFlat->reserve(roughNnz);
                }

                auto& localOffsets = offsetsByThread[static_cast<std::size_t>(tid)];
                localOffsets.clear();
                localOffsets.reserve((pEnd - pBeg) + 1u);
                localOffsets.push_back(0);

                std::array<std::vector<int>, 9> stencils;
                std::vector<int> tmpCellIndex;
                std::vector<double> tmpWraw;
                std::vector<int> mergedCellIndex;
                std::vector<double> mergedWraw;
                tmpCellIndex.reserve(64);
                tmpWraw.reserve(64);
                mergedCellIndex.reserve(64);
                mergedWraw.reserve(64);
                for (auto& s : stencils)
                    s.reserve(64);

                int localFallback = 0;

                for (std::size_t p = pBeg; p < pEnd; ++p)
                {
                    const double xp = positions[p].x;
                    const double yp = positions[p].y;

                    const double hp_in = hpIsConst ? hpConst : hp[p];
                    double hp_eff = hp_in;

                    int host = -1;
                    if (opt.hostCellLabels && p < opt.hostCellLabels->size())
                        host = (*opt.hostCellLabels)[p];
                    if (host < 0)
                        host = mesh.findCellLabel(Vector2D(xp, yp));
                    hostCellLabel[p] = host;

                    if (host < 0)
                    {
                        hpEff[p] = hp_eff;
                        inCompositeBand[p] = 0;
                        localOffsets.push_back(static_cast<int>(cellFlat.size()));
                        continue;
                    }

                    const double hHost = minCellSize(mesh.cells(host));
                    if (std::isfinite(hHost) && hHost > 0.0)
                        hp_eff = hHost;

                    double sx[3];
                    double sy[3];
                    int nx = 1;
                    int ny = 1;
                    if (nonPeriodic)
                    {
                        sx[0] = 0.0;
                        sy[0] = 0.0;
                    }
                    else
                    {
                        buildShiftLists(xp, yp, hp_eff, sx, nx, sy, ny);
                    }

                    const int nShifts = nx * ny;
                    const bool singleShift = (nShifts == 1);
                    for (int si = 0; si < nShifts; ++si)
                        stencils[static_cast<std::size_t>(si)].clear();

                    double minHK = 1.0e300;
                    double maxHK = 0.0;
                    if (singleShift)
                    {
                        ParticleStencilQuery::enumerateLeafCellsIntersectParticleSupport(
                            mesh, xp, yp, hp_eff, &stencils[0],
                            /*debugPrint=*/false, host, opt.deterministicOrdering);

                        for (int cellLabel : stencils[0])
                        {
                            const Cell& cell = mesh.cells(cellLabel);
                            const double hK = minCellSize(cell);
                            minHK = std::min(minHK, hK);
                            maxHK = std::max(maxHK, hK);
                        }
                    }
                    else
                    {
                        int si = 0;
                        for (int ix = 0; ix < nx; ++ix)
                        {
                            for (int iy = 0; iy < ny; ++iy)
                            {
                                const double xg = xp + sx[ix];
                                const double yg = yp + sy[iy];

                                int hostSeed = host;
                                if (!(sx[ix] == 0.0 && sy[iy] == 0.0))
                                {
                                    const double xs = clampToDomain(xg, xMin, xMax);
                                    const double ys = clampToDomain(yg, yMin, yMax);
                                    hostSeed = mesh.findCellLabel(Vector2D(xs, ys));
                                }

                                ParticleStencilQuery::enumerateLeafCellsIntersectParticleSupport(
                                    mesh, xg, yg, hp_eff, &stencils[static_cast<std::size_t>(si)],
                                    /*debugPrint=*/false, hostSeed, opt.deterministicOrdering);

                                for (int cellLabel : stencils[static_cast<std::size_t>(si)])
                                {
                                    const Cell& cell = mesh.cells(cellLabel);
                                    const double hK = minCellSize(cell);
                                    minHK = std::min(minHK, hK);
                                    maxHK = std::max(maxHK, hK);
                                }
                                ++si;
                            }
                        }
                    }

                    const bool hasLevelMismatch = (maxHK > minHK * (1.0 + 1e-12));
                    if (opt.enableCompositeBand && hasLevelMismatch && minHK < 1.0e299)
                    {
                        const double hp_star = minHK;
                        if (hp_star < hp_eff)
                        {
                            hp_eff = hp_star;
                            if (singleShift)
                            {
                                ParticleStencilQuery::enumerateLeafCellsIntersectParticleSupport(
                                    mesh, xp, yp, hp_eff, &stencils[0],
                                    /*debugPrint=*/false, host, opt.deterministicOrdering);
                            }
                            else
                            {
                                int si = 0;
                                for (int ix = 0; ix < nx; ++ix)
                                {
                                    for (int iy = 0; iy < ny; ++iy)
                                    {
                                        const double xg = xp + sx[ix];
                                        const double yg = yp + sy[iy];

                                        int hostSeed = host;
                                        if (!(sx[ix] == 0.0 && sy[iy] == 0.0))
                                        {
                                            const double xs = clampToDomain(xg, xMin, xMax);
                                            const double ys = clampToDomain(yg, yMin, yMax);
                                            hostSeed = mesh.findCellLabel(Vector2D(xs, ys));
                                        }

                                        ParticleStencilQuery::enumerateLeafCellsIntersectParticleSupport(
                                            mesh, xg, yg, hp_eff, &stencils[static_cast<std::size_t>(si)],
                                            /*debugPrint=*/false, hostSeed, opt.deterministicOrdering);
                                        ++si;
                                    }
                                }
                            }
                        }
                    }

                    const bool inBand = opt.enableCompositeBand ? hasLevelMismatch : false;
                    hpEff[p] = hp_eff;
                    inCompositeBand[p] = static_cast<unsigned char>(inBand ? 1 : 0);

                    mergedCellIndex.clear();
                    mergedWraw.clear();
                    if (singleShift)
                    {
                        fillWeightsForStencil(xp, yp, hp_eff, inBand, stencils[0], mergedCellIndex, mergedWraw);
                    }
                    else
                    {
                        int si = 0;
                        for (int ix = 0; ix < nx; ++ix)
                        {
                            for (int iy = 0; iy < ny; ++iy)
                            {
                                const double xg = xp + sx[ix];
                                const double yg = yp + sy[iy];
                                tmpCellIndex.clear();
                                tmpWraw.clear();
                                fillWeightsForStencil(xg, yg, hp_eff, inBand, stencils[static_cast<std::size_t>(si)],
                                                      tmpCellIndex, tmpWraw);
                                appendCellWeights(mergedCellIndex, mergedWraw, tmpCellIndex, tmpWraw);
                                ++si;
                            }
                        }
                        mergeSortedCellWeightsInPlace(mergedCellIndex, mergedWraw);
                    }

                    const double sumW = kahanSum(mergedWraw);
                    if (!(sumW > 0.0) || mergedCellIndex.empty())
                    {
                        if (opt.fallbackToHostOnZeroWeight)
                        {
                            const int hostIndex = mesh.cells(host).index();
                            cellFlat.push_back(hostIndex);
                            wFlat.push_back(1.0);
                            if (rawFlat)
                                rawFlat->push_back(0.0);
                            ++localFallback;
                        }
                        localOffsets.push_back(static_cast<int>(cellFlat.size()));
                        continue;
                    }

                    const std::size_t writeBeg = cellFlat.size();
                    int maxLocal = -1;
                    double maxVal = -1.0;
                    double sumTilde = 0.0;
                    for (std::size_t k = 0; k < mergedCellIndex.size(); ++k)
                    {
                        const double wtilde = mergedWraw[k] / sumW;
                        cellFlat.push_back(mergedCellIndex[k]);
                        wFlat.push_back(wtilde);
                        if (rawFlat)
                            rawFlat->push_back(mergedWraw[k]);
                        if (wtilde > maxVal)
                        {
                            maxVal = wtilde;
                            maxLocal = static_cast<int>(k);
                        }
                        sumTilde += wtilde;
                    }
                    const double r = 1.0 - sumTilde;
                    if (maxLocal >= 0)
                        wFlat[writeBeg + static_cast<std::size_t>(maxLocal)] += r;

                    localOffsets.push_back(static_cast<int>(cellFlat.size()));
                }

                fallbackByThread[static_cast<std::size_t>(tid)] = localFallback;
            }
            catch (...)
            {
                #pragma omp critical(lamp_cwcache_build_error)
                {
                    if (!firstError)
                        firstError = std::current_exception();
                }
            }
        }

        if (firstError)
            std::rethrow_exception(firstError);

        {
            double hpMin = std::numeric_limits<double>::infinity();
            double hpMax = 0.0;
            double hpSum = 0.0;
            std::size_t hpCount = 0;
            for (std::size_t p = 0; p < numParticles_; ++p)
            {
                if (hostCellLabel[p] < 0)
                    continue;
                const double hp_eff = hpEff[p];
                if (!std::isfinite(hp_eff) || hp_eff <= 0.0)
                    continue;
                hpMin = std::min(hpMin, hp_eff);
                hpMax = std::max(hpMax, hp_eff);
                hpSum += hp_eff;
                hpCount += 1;
            }
            if (opt.recordHpStats && hpCount > 0)
            {
                lastHpStats_.min = hpMin;
                lastHpStats_.max = hpMax;
                lastHpStats_.mean = hpSum / static_cast<double>(hpCount);
                lastHpStats_.count = hpCount;
                hasHpStats_ = true;
            }
        }

        if (opt.recordPerParticleHp)
            effectiveHp_ = hpEff;
        else
            effectiveHp_.clear();
        if (opt.recordPerParticleCompositeBand)
            inCompositeBand_ = inCompositeBand;
        else
            inCompositeBand_.clear();

        std::vector<std::size_t> threadBase(static_cast<std::size_t>(activeThreads) + 1u, 0u);
        for (int t = 0; t < activeThreads; ++t)
        {
            threadBase[static_cast<std::size_t>(t) + 1u] =
                threadBase[static_cast<std::size_t>(t)] + cellByThread[static_cast<std::size_t>(t)].size();
        }

        offsets_.assign(numParticles_ + 1u, 0);
        for (int t = 0; t < activeThreads; ++t)
        {
            const std::size_t pBeg = threadBegin[static_cast<std::size_t>(t)];
            const std::size_t pEnd = threadEnd[static_cast<std::size_t>(t)];
            const auto& localOffsets = offsetsByThread[static_cast<std::size_t>(t)];
            if (localOffsets.size() != (pEnd - pBeg + 1u))
                throw std::runtime_error("CouplingWeightsCache::build: inconsistent local offsets in parallel build");
            const int base = static_cast<int>(threadBase[static_cast<std::size_t>(t)]);
            for (std::size_t i = 0; i < localOffsets.size(); ++i)
                offsets_[pBeg + i] = base + localOffsets[i];
        }

        const std::size_t totalNnz = threadBase.back();
        cellIndices_.assign(totalNnz, 0);
        weights_.assign(totalNnz, 0.0);
        if (storeRawThisBuild)
            rawWeights_.assign(totalNnz, 0.0);

        std::exception_ptr copyError;
        #pragma omp parallel for num_threads(activeThreads) schedule(static)
        for (int t = 0; t < activeThreads; ++t)
        {
            try
            {
                const std::size_t base = threadBase[static_cast<std::size_t>(t)];
                const auto& srcCell = cellByThread[static_cast<std::size_t>(t)];
                const auto& srcW = weightByThread[static_cast<std::size_t>(t)];
                if (srcCell.size() != srcW.size())
                    throw std::runtime_error("CouplingWeightsCache::build: inconsistent thread CSR data");
                if (!srcCell.empty())
                {
                    std::memcpy(cellIndices_.data() + base, srcCell.data(), srcCell.size() * sizeof(int));
                    std::memcpy(weights_.data() + base, srcW.data(), srcW.size() * sizeof(double));
                    if (storeRawThisBuild)
                    {
                        const auto& srcRaw = rawByThread[static_cast<std::size_t>(t)];
                        if (srcRaw.size() != srcW.size())
                            throw std::runtime_error("CouplingWeightsCache::build: inconsistent thread raw CSR data");
                        std::memcpy(rawWeights_.data() + base, srcRaw.data(), srcRaw.size() * sizeof(double));
                    }
                }
            }
            catch (...)
            {
                #pragma omp critical(lamp_cwcache_copy_error)
                {
                    if (!copyError)
                        copyError = std::current_exception();
                }
            }
        }
        if (copyError)
            std::rethrow_exception(copyError);

        fallbackCountTotal = 0;
        for (int v : fallbackByThread)
            fallbackCountTotal += v;

    }

    lastBuildFallbackCount_ = static_cast<std::size_t>(std::max(0, fallbackCountTotal));
    if (fallbackCountTotal > 0)
    {
        static std::size_t warnCount = 0;
        ++warnCount;
        const bool verbose = ConsoleControl::enabled(ConsoleControl::Level::kDebug);
        const bool shouldPrint = verbose || (warnCount <= 3);
        if (shouldPrint)
        {
            std::cerr << "[CouplingWeightsCache] WARNING: " << fallbackCountTotal
                      << " particles had zero-weight stencils; fell back to host cell.";
            if (!verbose && warnCount >= 3)
                std::cerr << " (suppressed; enable debug for full logs)";
            std::cerr << std::endl;
        }
    }

    // Optional debug dump (done after building; safe in both serial and parallel builds).
    if (willDumpThisBuild)
    {
        if (dumpParticleWeightsDetailedCSV(mesh, opt.levelSetI, opt.alphaByCellIndex, opt.dumpWeightsPath, opt.debugParticleId))
        {
            dumpedWeightsOnce_ = true;
            std::cout << "[CouplingWeightsCache] dumped particle weights to: " << opt.dumpWeightsPath << std::endl;
        }
        else
        {
            std::cerr << "[CouplingWeightsCache] WARNING: failed to dump particle weights to: " << opt.dumpWeightsPath << std::endl;
        }
    }

    // Drop raw weights unless the user requested keeping them.
    if (!keepRawWeights_)
        rawWeights_.clear();
}

bool CouplingWeightsCache::dumpParticleWeightsCSV(const std::string& filePath,
                                                  const std::size_t particleId) const
{
    std::ofstream out(filePath);
    if (!out)
        return false;

    out << "particleId,entry,cellIndex,weight";
    if (rawWeights_.size() == weights_.size())
        out << ",rawWeight";
    out << "\n";

    if (particleId >= numParticles_ || offsets_.empty())
        return true;

    const int begin = offsets_.at(particleId);
    const int end = offsets_.at(particleId + 1);
    for (int e = begin; e < end; ++e)
    {
        out << particleId << ',' << (e - begin) << ',' << cellIndices_.at(e) << ',' << weights_.at(e);
        if (rawWeights_.size() == weights_.size() && (std::size_t)e < rawWeights_.size())
            out << ',' << rawWeights_.at(e);
        out << "\n";
    }
    return true;
}

bool CouplingWeightsCache::dumpParticleWeightsDetailedCSV(const BackgroundMesh& mesh,
                                                          const int levelSetI,
                                                          const std::vector<double>* alphaByCellIndex,
                                                          const std::string& filePath,
                                                          const std::size_t particleId) const
{
    std::ofstream out(filePath);
    if (!out)
        return false;

    // Columns in the detailed dump:
    //   (cell, center, alpha, Wraw, Wtilde)
    // We also include cellIndex/isCut for convenience.
    out << "particleId,entry,cellLabel,cellIndex,centerX,centerY,isCut,alpha,Wraw,Wtilde\n";

    if (particleId >= numParticles_ || offsets_.empty())
        return true;

    const bool hasRaw = (rawWeights_.size() == weights_.size());

    const int begin = offsets_.at(particleId);
    const int end = offsets_.at(particleId + 1);
    for (int e = begin; e < end; ++e)
    {
        const int cellIndex = cellIndices_.at(e);
        int cellLabel = -1;
        if (mesh.cellIndexToLabel)
            cellLabel = mesh.cellIndexToLabel[cellIndex];

        Vector2D cc(0.0, 0.0);
        bool isCut = false;
        if (cellLabel >= 0)
        {
            const Cell& c = mesh.cells(cellLabel);
            cc = c.centerCoord();
            if (levelSetI >= 0 && levelSetI < static_cast<int>(c.isCut.size()))
                isCut = (c.isCut[levelSetI] != 0);
        }

        double alpha = 1.0;
        if (alphaByCellIndex && (std::size_t)cellIndex < alphaByCellIndex->size())
            alpha = (*alphaByCellIndex)[cellIndex];
        alpha = clamp01(alpha);

        if (!isCut)
        {
            const bool outOfRange = (levelSetI < 0);
            const bool unknownCut = (cellLabel < 0) || (levelSetI >= static_cast<int>(mesh.cells(cellLabel).isCut.size()));
            if (outOfRange || unknownCut)
                isCut = (alpha > 1e-12) && (alpha < (1.0 - 1e-12));
        }

        const double wtilde = weights_.at(e);
        const double wraw = hasRaw ? rawWeights_.at(e) : std::numeric_limits<double>::quiet_NaN();

        out << particleId << ',' << (e - begin) << ',' << cellLabel << ',' << cellIndex << ','
            << cc.x << ',' << cc.y << ',' << (isCut ? 1 : 0) << ',' << alpha << ',' << wraw << ',' << wtilde << "\n";
    }

    return true;
}
