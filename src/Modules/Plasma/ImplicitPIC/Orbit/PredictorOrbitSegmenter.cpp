#include "Modules/Plasma/ImplicitPIC/Orbit/PredictorOrbitSegmenter.h"

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "Core/Mesh/BackgroundMesh.h"

#include "Modules/Plasma/ExplicitPIC/Deposit/CICMomentExactCut.h"

namespace ImplicitPIC
{
namespace Orbit
{
namespace
{
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

    Vector2D out = x;

    if (px)
    {
        const double xmin = mesh.minCoord.x;
        const double xmax = mesh.maxCoord.x;
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
        const double ymin = mesh.minCoord.y;
        const double ymax = mesh.maxCoord.y;
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

inline double clamp01Closed(double s)
{
    if (s < 0.0)
        return 0.0;
    if (s > 1.0)
        return 1.0;
    return s;
}

inline bool isCutByAlpha(double alpha, double eps)
{
    alpha = std::clamp(alpha, 0.0, 1.0);
    return (alpha > eps) && (alpha < (1.0 - eps));
}

// Return true if the segment from xA to xB (within a single leaf cell) crosses the EB
// half-plane from fluid into solid, and output the hit time fraction sHit in (0,1].
inline bool firstEmbeddedHitFractionInLeaf(const BackgroundMesh& mesh,
                                           int cellLabel,
                                           const Vector2D& xA,
                                           const Vector2D& xB,
                                           const PredictorOrbitSegmenter::Options& opt,
                                           double* outSHit)
{
    if (outSHit)
        *outSHit = 0.0;
    if (!opt.enableEmbeddedIntercept)
        return false;
    if (!opt.embeddedLevelSetVert)
        return false;
    if (!opt.alphaByCellIndex)
        return false;

    if (cellLabel < 0 || cellLabel >= mesh.cellLabelCapacity())
        return false;

    const int cellIndex = mesh.cells(cellLabel).index();
    if (cellIndex < 0 || static_cast<std::size_t>(cellIndex) >= opt.alphaByCellIndex->size())
        return false;

    const double alpha = (*opt.alphaByCellIndex)[static_cast<std::size_t>(cellIndex)];
    if (!isCutByAlpha(alpha, opt.cutAlphaEps))
        return false;

    const Cell& cell = mesh.cells(cellLabel);
    CICMomentExactCut::Plane pl;
    if (!CICMomentExactCut::buildCellPlaneFromVertexLevelSet(mesh, cell, *opt.embeddedLevelSetVert, &pl))
        return false;

    auto phi = [&](const Vector2D& x) -> double {
        return pl.a * x.x + pl.b * x.y + pl.c;
    };

    const double phiA = phi(xA);
    const double phiB = phi(xB);

    // Fluid-side condition matches the convention used by EB Poisson + deposition:
    // insideIsNegative==true => solid: psi<0, fluid: psi>=0.
    auto isFluid = [&](double psi) -> bool {
        return opt.embeddedInsideIsNegative ? (psi >= 0.0) : (psi <= 0.0);
    };

    if (!isFluid(phiA))
    {
        // Already in solid (numerical or bad hint). Treat as immediate termination.
        // Caller will emit a zero-length segment (ignored) and stop.
        if (outSHit)
            *outSHit = 0.0;
        return true;
    }
    if (isFluid(phiB))
        return false;

    const double denom = (phiB - phiA);
    if (!std::isfinite(denom) || denom == 0.0)
        return false;

    // Solve phiA + s*(phiB-phiA) = 0.
    double sHit = -phiA / denom;
    if (!std::isfinite(sHit))
        return false;
    sHit = clamp01Closed(sHit);
    if (!(sHit > 0.0))
        return false;

    if (outSHit)
        *outSHit = sHit;
    return true;
}

// Return first analytic exit time from an axis-aligned leaf rectangle.
// The trajectory is x(t) = xStart + (t - tStart) * velocity.
inline bool firstAxisAlignedExitTimeInLeaf(const BackgroundMesh& mesh,
                                           int cellLabel,
                                           const Vector2D& xStart,
                                           const Vector2D& velocity,
                                           double tStart,
                                           double tEnd,
                                           double inCellEps,
                                           double* outTExit)
{
    if (outTExit)
        *outTExit = tStart;
    if (cellLabel < 0 || cellLabel >= mesh.cellLabelCapacity())
        return false;
    if (!(tEnd > tStart))
        return false;

    const Cell& c = mesh.cells(cellLabel);
    const Vector2D cc = c.centerCoord();
    const Vector2D sz = c.size();
    const double eps = std::abs(inCellEps);
    const double xL = cc.x - 0.5 * sz.x - eps;
    const double xR = cc.x + 0.5 * sz.x + eps;
    const double yL = cc.y - 0.5 * sz.y - eps;
    const double yR = cc.y + 0.5 * sz.y + eps;

    double tExit = std::numeric_limits<double>::infinity();
    auto update = [&](double tCand) {
        if (std::isfinite(tCand) && (tCand > tStart) && (tCand < tExit))
            tExit = tCand;
    };

    if (velocity.x > 0.0)
        update(tStart + (xR - xStart.x) / velocity.x);
    else if (velocity.x < 0.0)
        update(tStart + (xL - xStart.x) / velocity.x);

    if (velocity.y > 0.0)
        update(tStart + (yR - xStart.y) / velocity.y);
    else if (velocity.y < 0.0)
        update(tStart + (yL - xStart.y) / velocity.y);

    if (!std::isfinite(tExit) || !(tExit < tEnd))
        return false;

    if (outTExit)
        *outTExit = tExit;
    return true;
}
} // namespace

void PredictorOrbitSegmenter::segmentLinearOrbits(const BackgroundMesh& mesh,
                                                  const std::vector<Vector2D>& x0,
                                                  const std::vector<Vector2D>& x1,
                                                  const std::vector<int>* hostCellHints,
                                                  const double dt,
                                                  std::vector<int>& segmentOffsets,
                                                  std::vector<OrbitSegment>& segments,
                                                  const Options& opt,
                                                  Stats* stats)
{
    const bool trace = opt.trace;
    const bool noHint = opt.ignoreHostHints;

    auto clampCellToMaxLevel = [&](int cellLabel) -> int {
        if (cellLabel < 0 || cellLabel >= mesh.cellLabelCapacity())
            return -1;
        const int maxLevel = opt.segmentCellMaxLevel;
        if (maxLevel < 0)
            return cellLabel;
        // EB interception splitting assumes leaf/cut-leaf semantics.
        if (opt.enableEmbeddedIntercept)
            return cellLabel;

        const int original = cellLabel;
        int climbed = 0;
        int c = cellLabel;
        while (c >= 0)
        {
            if (c >= mesh.cellLabelCapacity())
                return -1;
            const Cell& cell = mesh.cells(c);
            if (cell.label() == -1)
                return -1;
            if (cell.level() <= maxLevel)
            {
                if (stats && c >= 0 && c != original)
                {
                    stats->coarsenHits += 1u;
                    stats->coarsenParentClimb += static_cast<std::uint64_t>(std::max(0, climbed));
                }
                return c;
            }
            c = cell.parent;
            ++climbed;
        }
        return -1;
    };

    const std::size_t nP = x0.size();

    segmentOffsets.assign(nP + 1u, 0);
    segments.clear();

    if (nP == 0u)
        return;
    if (x1.size() != nP)
        throw std::runtime_error("[PredictorOrbitSegmenter] x0/x1 size mismatch");
    if (dt <= 0.0)
        throw std::runtime_error("[PredictorOrbitSegmenter] dt must be positive");
    if (hostCellHints && hostCellHints->size() != nP)
        throw std::runtime_error("[PredictorOrbitSegmenter] hostCellHints size mismatch");

    const bool px = opt.periodicX;
    const bool py = opt.periodicY;

    auto inDomain = [&](const Vector2D& x) -> bool {
        if (!std::isfinite(x.x) || !std::isfinite(x.y))
            return false;
        if (!px)
        {
            if (x.x < mesh.minCoord.x || x.x > mesh.maxCoord.x)
                return false;
        }
        if (!py)
        {
            if (x.y < mesh.minCoord.y || x.y > mesh.maxCoord.y)
                return false;
        }
        return true;
    };

    auto posAt = [&](const Vector2D& a, const Vector2D& b, double t) -> Vector2D {
        const double s = std::clamp(t / dt, 0.0, 1.0);
        const Vector2D x = a * (1.0 - s) + b * s;
        return wrapToDomain(mesh, x, px, py, opt.wrapEps);
    };

    auto isSolidCell = [&](int cellLabel) -> bool {
        if (!opt.terminateOnSolidCell || !opt.alphaByCellIndex)
            return false;
        if (cellLabel < 0 || cellLabel >= mesh.cellLabelCapacity())
            return false;
        const int cellIndex = mesh.cells(cellLabel).index();
        if (cellIndex < 0 || static_cast<std::size_t>(cellIndex) >= opt.alphaByCellIndex->size())
            return false;
        double alpha = (*opt.alphaByCellIndex)[static_cast<std::size_t>(cellIndex)];
        if (!std::isfinite(alpha))
            alpha = 0.0;
        alpha = std::clamp(alpha, 0.0, 1.0);
        return alpha <= opt.solidAlphaEps;
    };

    auto cellFromHint = [&](int hint, const Vector2D& x) -> int {
        // Non-periodic boundaries: once outside the domain, stop segmenting.
        if (!inDomain(x))
            return -1;

        int c = noHint ? -1 : hint;
        if (!noHint)
        {
            if (c >= 0 && c < mesh.cellLabelCapacity())
                c = mesh.findCellLabelFromHint(c, x, opt.inCellEps);
        }
        if (c < 0 || c >= mesh.cellLabelCapacity())
            c = mesh.findCellLabel(x);
        if (c < 0 || c >= mesh.cellLabelCapacity())
            return -1;
        return clampCellToMaxLevel(c);
    };

    const double tTol = std::max(opt.timeTolAbs, opt.timeTolRel * dt);

    auto segmentOneParticle = [&](std::size_t p,
                                  std::vector<OrbitSegment>* outVec,
                                  std::uint64_t* outFastExitHits,
                                  std::uint64_t* outBisectionFallbackHits) -> int
    {
        int written = 0;
        std::uint64_t fastExitHitsLocal = 0u;
        std::uint64_t bisectionFallbackHitsLocal = 0u;
        auto emitSeg = [&](const OrbitSegment& seg) {
            if (!(seg.t1 > seg.t0))
                return;
            if (outVec)
            {
                outVec->push_back(seg);
            }
            ++written;
        };

        const Vector2D a = wrapToDomain(mesh, x0[p], px, py, opt.wrapEps);
        const Vector2D b = wrapToDomain(mesh, x1[p], px, py, opt.wrapEps);
        const Vector2D velocity = (b - a) * (1.0 / dt);
        const bool allowAnalyticExit =
            (!px && !py) && (std::isfinite(velocity.x) && std::isfinite(velocity.y));

        const int hint0 = hostCellHints ? (*hostCellHints)[p] : -1;

        if (trace)
        {
            std::cout << "[SegTrace] p=" << p
                      << " hint=" << hint0
                      << " x0=" << a
                      << " x1=" << b
                      << " dt=" << dt
                      << " noHint=" << (noHint ? 1 : 0)
                      << std::endl;
        }

        int cell0 = cellFromHint(hint0, a);
        if (cell0 < 0)
        {
            return written;
        }
        if (isSolidCell(cell0))
        {
            // In solid at predictor start: emit no orbit-averaged segment.
            return written;
        }

        double tStart = 0.0;
        Vector2D xStart = a;
        int segmentsBuilt = 0;

        while (tStart < dt - tTol)
        {
            if (++segmentsBuilt > opt.maxSegmentsPerParticle)
                throw std::runtime_error("[PredictorOrbitSegmenter] exceeded maxSegmentsPerParticle");

            if (trace)
            {
                std::cout << "[SegTrace] p=" << p
                          << " seg=" << segmentsBuilt
                          << " tStart=" << tStart
                          << " cell0=" << cell0
                          << " xStart=" << xStart
                          << std::endl;
            }

            const Vector2D xEnd = b;
            const double tEnd = dt;

            if (cell0 < 0 || cell0 >= mesh.cellLabelCapacity())
            {
                cell0 = cellFromHint(cell0, xStart);
                if (cell0 < 0)
                    break;
            }
            if (isSolidCell(cell0))
            {
                // Once entering a fully solid cell, stop segmenting this particle.
                tStart = dt;
                break;
            }

            // If the endpoint is still inside the current cell, close the segment.
            const bool insideEnd = mesh.isInsideCell(cell0, xEnd, opt.inCellEps);
            if (insideEnd)
            {
                // Optional: split/terminate at first EB interception inside this cut leaf.
                double sHit = 0.0;
                const bool hit = firstEmbeddedHitFractionInLeaf(mesh, cell0, xStart, xEnd, opt, &sHit);
                if (hit)
                {
                    const double tHit = tStart + sHit * (tEnd - tStart);
                    const Vector2D xHit = posAt(a, b, tHit);

                    OrbitSegment seg;
                    seg.cellLabel = cell0;
                    seg.t0 = tStart;
                    seg.t1 = tHit;
                    seg.x0 = xStart;
                    seg.x1 = xHit;
                    emitSeg(seg);

                    // Absorbing EB: stop further segments.
                    tStart = dt;
                    break;
                }

                OrbitSegment seg;
                seg.cellLabel = cell0;
                seg.t0 = tStart;
                seg.t1 = tEnd;
                seg.x0 = xStart;
                seg.x1 = xEnd;
                emitSeg(seg);
                tStart = dt;
                break;
            }

            const double tNudge = std::max(4.0 * tTol, 1e-15);

            // Otherwise, find first exit time from cell0.
            double lo = tStart;
            double hi = dt;
            Vector2D xLo = xStart;

            // Ensure bracketing: inside(lo)=true, inside(hi)=false.
            if (!mesh.isInsideCell(cell0, xLo, opt.inCellEps))
            {
                // Degenerate (numerical): re-seed cell0 from position.
                cell0 = cellFromHint(cell0, xLo);
                if (cell0 < 0)
                    break;
            }

            bool usedAnalyticExit = false;
            if (allowAnalyticExit)
            {
                double tExit = tStart;
                if (firstAxisAlignedExitTimeInLeaf(mesh, cell0, xStart, velocity, tStart, dt, opt.inCellEps, &tExit))
                {
                    double loCand = tExit;
                    Vector2D xLoCand = posAt(a, b, loCand);
                    if (!mesh.isInsideCell(cell0, xLoCand, opt.inCellEps))
                    {
                        loCand = std::max(tStart, tExit - tNudge);
                        xLoCand = posAt(a, b, loCand);
                    }

                    const double hiCand = std::min(dt, tExit + tNudge);
                    if (hiCand > loCand)
                    {
                        const Vector2D xHiCand = posAt(a, b, hiCand);
                        if (mesh.isInsideCell(cell0, xLoCand, opt.inCellEps) &&
                            !mesh.isInsideCell(cell0, xHiCand, opt.inCellEps))
                        {
                            lo = loCand;
                            hi = hiCand;
                            xLo = xLoCand;
                            usedAnalyticExit = true;
                            ++fastExitHitsLocal;
                        }
                    }
                }
            }

            if (!usedAnalyticExit)
            {
                ++bisectionFallbackHitsLocal;
                for (int it = 0; it < opt.maxBisectionIter; ++it)
                {
                    if ((hi - lo) <= tTol)
                        break;
                    const double mid = 0.5 * (lo + hi);
                    const Vector2D xMid = posAt(a, b, mid);
                    const bool insideMid = mesh.isInsideCell(cell0, xMid, opt.inCellEps);
                    if (insideMid)
                    {
                        lo = mid;
                        xLo = xMid;
                    }
                    else
                    {
                        hi = mid;
                    }
                }
            }

            // Emit segment [tStart, lo] (still inside cell0).
            {
                // Optional: split/terminate at first EB interception before the leaf exit.
                double sHit = 0.0;
                const bool hit = firstEmbeddedHitFractionInLeaf(mesh, cell0, xStart, xLo, opt, &sHit);
                if (hit)
                {
                    const double tHit = tStart + sHit * (lo - tStart);
                    const Vector2D xHit = posAt(a, b, tHit);

                    OrbitSegment seg;
                    seg.cellLabel = cell0;
                    seg.t0 = tStart;
                    seg.t1 = tHit;
                    seg.x0 = xStart;
                    seg.x1 = xHit;
                    emitSeg(seg);

                    // Absorbing EB: stop further segments.
                    tStart = dt;
                    break;
                }

                OrbitSegment seg;
                seg.cellLabel = cell0;
                seg.t0 = tStart;
                seg.t1 = lo;
                seg.x0 = xStart;
                seg.x1 = xLo;
                emitSeg(seg);
            }

            // Move slightly past the crossing to determine the next cell.
            // NOTE: Using only a relative tolerance (tTol) can be too small in double
            // precision, causing xStart to remain numerically on the boundary and the
            // loop to make no progress (eventually hitting maxSegmentsPerParticle).
            // Nudge time by a small absolute floor to ensure forward progress.
            const double tAfter = std::min(dt, hi + tNudge);
            if (!(tAfter > tStart))
            {
                // Degenerate: cannot advance (numerical). Stop this particle.
                tStart = dt;
                break;
            }

            tStart = tAfter;
            xStart = posAt(a, b, tStart);
            cell0 = cellFromHint(cell0, xStart);
            if (cell0 < 0)
            {
                // Could be outside domain in non-periodic cases; stop.
                tStart = dt;
                break;
            }
            if (isSolidCell(cell0))
            {
                tStart = dt;
                break;
            }
        }

        if (outFastExitHits)
            *outFastExitHits = fastExitHitsLocal;
        if (outBisectionFallbackHits)
            *outBisectionFallbackHits = bisectionFallbackHitsLocal;
        return written;
    };

    std::vector<int> segCountByParticle(nP, 0);
    std::vector<std::uint64_t> fastExitHitsByParticle(nP, 0u);
    std::vector<std::uint64_t> bisectionFallbackHitsByParticle(nP, 0u);

#ifdef _OPENMP
    const bool canParallel = (!trace) && (nP >= 64u);
    if (canParallel)
    {
        const int maxThreads = std::max(1, omp_get_max_threads());
        const int activeThreads = std::max(1, std::min<int>(maxThreads, static_cast<int>(nP)));
        std::vector<std::vector<OrbitSegment>> segByThread(static_cast<std::size_t>(activeThreads));
        std::vector<std::vector<int>> localOffsetsByThread(static_cast<std::size_t>(activeThreads));
        std::vector<std::size_t> threadBegin(static_cast<std::size_t>(activeThreads), 0u);
        std::vector<std::size_t> threadEnd(static_cast<std::size_t>(activeThreads), 0u);

        std::exception_ptr firstError;
        // Schedule every logical chunk even if the runtime forms a smaller team.
        #pragma omp parallel for num_threads(activeThreads) schedule(static)
        for (int tid = 0; tid < activeThreads; ++tid)
        {
            const std::size_t pBeg = (nP * static_cast<std::size_t>(tid)) / static_cast<std::size_t>(activeThreads);
            const std::size_t pEnd = (nP * static_cast<std::size_t>(tid + 1)) / static_cast<std::size_t>(activeThreads);
            threadBegin[static_cast<std::size_t>(tid)] = pBeg;
            threadEnd[static_cast<std::size_t>(tid)] = pEnd;

            auto& segFlat = segByThread[static_cast<std::size_t>(tid)];
            auto& localOffsets = localOffsetsByThread[static_cast<std::size_t>(tid)];
            segFlat.clear();
            localOffsets.clear();
            localOffsets.reserve((pEnd - pBeg) + 1u);
            localOffsets.push_back(0);

            try
            {
                for (std::size_t p = pBeg; p < pEnd; ++p)
                {
                    std::uint64_t fastHits = 0u;
                    std::uint64_t fallbackHits = 0u;
                    const int produced = segmentOneParticle(p, &segFlat, &fastHits, &fallbackHits);
                    segCountByParticle[p] = produced;
                    fastExitHitsByParticle[p] = fastHits;
                    bisectionFallbackHitsByParticle[p] = fallbackHits;
                    localOffsets.push_back(static_cast<int>(segFlat.size()));
                }
            }
            catch (...)
            {
                #pragma omp critical(lamp_segmenter_error)
                {
                    if (!firstError)
                        firstError = std::current_exception();
                }
            }
        }
        if (firstError)
            std::rethrow_exception(firstError);

        std::vector<std::size_t> threadBase(static_cast<std::size_t>(activeThreads) + 1u, 0u);
        for (int t = 0; t < activeThreads; ++t)
        {
            threadBase[static_cast<std::size_t>(t) + 1u] =
                threadBase[static_cast<std::size_t>(t)] + segByThread[static_cast<std::size_t>(t)].size();
        }
        const std::size_t totalSegs = threadBase.back();
        segments.resize(totalSegs);

        for (int t = 0; t < activeThreads; ++t)
        {
            const std::size_t pBeg = threadBegin[static_cast<std::size_t>(t)];
            const std::size_t pEnd = threadEnd[static_cast<std::size_t>(t)];
            const auto& localOffsets = localOffsetsByThread[static_cast<std::size_t>(t)];
            if (localOffsets.size() != (pEnd - pBeg + 1u))
                throw std::runtime_error("[PredictorOrbitSegmenter] inconsistent local offsets");

            const int base = static_cast<int>(threadBase[static_cast<std::size_t>(t)]);
            for (std::size_t i = 0; i < localOffsets.size(); ++i)
                segmentOffsets[pBeg + i] = base + localOffsets[i];
        }

        std::exception_ptr copyError;
        #pragma omp parallel for num_threads(activeThreads) schedule(static)
        for (int t = 0; t < activeThreads; ++t)
        {
            try
            {
                const std::size_t base = threadBase[static_cast<std::size_t>(t)];
                const auto& src = segByThread[static_cast<std::size_t>(t)];
                if (!src.empty())
                {
                    std::copy(src.begin(), src.end(), segments.data() + base);
                }
            }
            catch (...)
            {
                #pragma omp critical(lamp_segmenter_fill_error)
                {
                    if (!copyError)
                        copyError = std::current_exception();
                }
            }
        }
        if (copyError)
            std::rethrow_exception(copyError);
    }
    else
#endif
    {
        segments.clear();
        segmentOffsets[0] = 0;
        for (std::size_t p = 0; p < nP; ++p)
        {
            std::uint64_t fastHits = 0u;
            std::uint64_t fallbackHits = 0u;
            const int produced = segmentOneParticle(p, &segments, &fastHits, &fallbackHits);
            segCountByParticle[p] = produced;
            fastExitHitsByParticle[p] = fastHits;
            bisectionFallbackHitsByParticle[p] = fallbackHits;
            segmentOffsets[p + 1u] = static_cast<int>(segments.size());
        }
    }

    if (stats)
    {
        stats->fastExitHits = 0u;
        stats->bisectionFallbackHits = 0u;
        for (std::size_t p = 0; p < nP; ++p)
        {
            stats->fastExitHits += fastExitHitsByParticle[p];
            stats->bisectionFallbackHits += bisectionFallbackHitsByParticle[p];
        }
    }
}

} // namespace Orbit
} // namespace ImplicitPIC
