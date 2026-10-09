#include "Solvers/Poisson/CutCellMetrics.h"

#include <fstream>
#include <iostream>
#include <cmath>
#include <limits>

#include <algorithm>
#include <map>
#include <unordered_set>

#include "Core/IO/ConsoleControl.h"

#include "Core/Mesh/Edge.h"

#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/Cell.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Geometry/ImmersedGeometry.h"

#include "Core/Base/Interpolation.h"

namespace Solvers
{
namespace Poisson
{

namespace
{
static inline double sqr(const double x) { return x * x; }

static inline double clamp01(const double t)
{
    if (t < 0.0) return 0.0;
    if (t > 1.0) return 1.0;
    return t;
}

static inline double embeddedDistanceFloor(const double h)
{
    // Fixed absolute floors are too weak once MovingEB creates near-tangent
    // cut cells on a refined mesh. Keep the EB center-to-interface distance
    // proportional to the local cell size to avoid pathological A_gamma/d_gamma.
    const double hSafe = (std::isfinite(h) && h > 0.0) ? h : 1.0;
    return std::max(1e-12, 2.0e-2 * hSafe);
}

static inline double pointSegmentDist2(const Vector2D& p, const Vector2D& a, const Vector2D& b)
{
    const Vector2D ab = b - a;
    const Vector2D ap = p - a;
    const double denom = ab.x * ab.x + ab.y * ab.y;
    if (!(denom > 0.0))
        return sqr(p.x - a.x) + sqr(p.y - a.y);

    const double t = clamp01((ap.x * ab.x + ap.y * ab.y) / denom);
    const Vector2D q = a + ab * t;
    return sqr(p.x - q.x) + sqr(p.y - q.y);
}

static std::vector<Vector2D> discretizeLineBySign(const std::vector<Line>& lines, const SignInt& line)
{
    std::vector<Vector2D> discretePoints;

    int lineIdx = line.value;
    if (lineIdx < 0 || lineIdx >= static_cast<int>(lines.size()))
    {
        if (lineIdx > 0 && (lineIdx - 1) >= 0 && (lineIdx - 1) < static_cast<int>(lines.size()))
            lineIdx -= 1;
    }
    if (lineIdx < 0 || lineIdx >= static_cast<int>(lines.size()))
        return discretePoints;

    const int numPoints = 50;
    const Line& L = lines[static_cast<size_t>(lineIdx)];

    switch (L.type())
    {
        case LineType::straight:
        {
            for (int i = 0; i <= numPoints; ++i)
            {
                const double t = static_cast<double>(i) / numPoints;
                Vector2D p = line.isNegative
                    ? L.points[1] + (L.points[0] - L.points[1]) * t
                    : L.points[0] + (L.points[1] - L.points[0]) * t;
                discretePoints.push_back(p);
            }
        }
        break;

        case LineType::arc:
        {
            Vector2D start = L.points[0];
            const double radius = L.radius();
            const double radian = L.radian();
            Vector2D center = L.arcCenter();
            const double startAngle = std::atan2(start.y - center.y, start.x - center.x);
            const double endAngle = startAngle + radian;

            for (int i = 0; i <= numPoints; ++i)
            {
                const double t = static_cast<double>(i) / numPoints;
                const double angle = line.isNegative ? endAngle - t * radian * Pi
                                                     : startAngle + t * radian * Pi;
                Vector2D p = center + Vector2D(radius * std::cos(angle), radius * std::sin(angle));
                discretePoints.push_back(p);
            }
        }
        break;

        case LineType::spline:
        {
            BSplineCurve bspline(L.featurePoints(), 3);
            for (int i = 0; i <= numPoints; ++i)
            {
                const double t = static_cast<double>(i) / numPoints;
                Vector2D p = line.isNegative ? bspline.evaluate(1 - t) : bspline.evaluate(t);
                discretePoints.push_back(p);
            }
        }
        break;

        default:
            break;
    }

    return discretePoints;
}

static inline double pointPolylineDist2(const Vector2D& p, const std::vector<Vector2D>& poly)
{
    if (poly.empty())
        return std::numeric_limits<double>::infinity();
    if (poly.size() == 1)
        return sqr(p.x - poly[0].x) + sqr(p.y - poly[0].y);
    double best = std::numeric_limits<double>::infinity();
    for (size_t i = 1; i < poly.size(); ++i)
        best = std::min(best, pointSegmentDist2(p, poly[i - 1], poly[i]));
    return best;
}

// Add a point to list if it is not duplicated within eps.
static void addUniquePoint(std::vector<Vector2D>& pts, const Vector2D& p, const double eps2)
{
    for (const auto& q : pts)
    {
        const double d2 = sqr(p.x - q.x) + sqr(p.y - q.y);
        if (d2 <= eps2) return;
    }
    pts.push_back(p);
}

// Compute the (approximate) two intersection points of psi=0 with a line segment [a,b],
// using linear interpolation on endpoint values. Add any found points to `out`.
static void addZeroCrossingsOnEdge(std::vector<Vector2D>& out,
                                  const Vector2D& a, const Vector2D& b,
                                  const double psiA, const double psiB,
                                  const double eps2)
{
    // Endpoint hits
    if (psiA == 0.0) addUniquePoint(out, a, eps2);
    if (psiB == 0.0) addUniquePoint(out, b, eps2);

    // Sign change -> one intersection
    if (psiA * psiB < 0.0)
    {
        const double t = psiA / (psiA - psiB); // in (0,1)
        const Vector2D p = a + (b - a) * t;
        addUniquePoint(out, p, eps2);
    }
}

// Pick two points that are farthest apart. Returns false if <2.
static bool farthestPair(const std::vector<Vector2D>& pts, Vector2D& p0, Vector2D& p1)
{
    if (pts.size() < 2) return false;
    double best = -1.0;
    size_t iBest = 0, jBest = 1;
    for (size_t i = 0; i < pts.size(); ++i)
    {
        for (size_t j = i + 1; j < pts.size(); ++j)
        {
            const double d2 = sqr(pts[i].x - pts[j].x) + sqr(pts[i].y - pts[j].y);
            if (d2 > best)
            {
                best = d2;
                iBest = i;
                jBest = j;
            }
        }
    }
    p0 = pts[iBest];
    p1 = pts[jBest];
    return true;
}
} // namespace

void CutCellMetrics::clear()
{
    cell_.clear();
    beta_.clear();
    active_.clear();
    levelSetI_ = 0;
    ++revision_;
}

bool CutCellMetrics::isActiveCellIndex(const int cellIndex) const
{
    if (cellIndex < 0 || cellIndex >= static_cast<int>(active_.size())) return false;
    return active_[cellIndex] != 0;
}

void CutCellMetrics::build(
    const QuadtreeMesh& mesh,
    ImmersedGeometry<QuadtreeMesh>& iGeometry,
    const FaceSegmentGraph& faceGraph,
    const BuildOptions& opt)
{
    // Build per-cell metrics:
    //   - default values for (A_gamma, d_gamma, n_gamma)
    //   - beta_f from Edge::cutLength()/length()
    //   - alpha_K from existing cut info if available, otherwise sampling

    levelSetI_ = opt.levelSetI;
    ++revision_;
    const int alphaN = std::max(1, opt.alphaSampleN);
    const bool insideNeg = opt.insideIsNegative;

    const bool useOverride = (opt.psiOverride != nullptr);
    if (useOverride && (!opt.psiDxOverride || !opt.psiDyOverride))
    {
        throw std::runtime_error("[CutCellMetrics] psiOverride requires psiDxOverride and psiDyOverride.");
    }

    const auto& psiField = useOverride ? *opt.psiOverride : iGeometry.levelSets.at(levelSetI_);
    const auto& psiDxField = useOverride ? *opt.psiDxOverride : iGeometry.LS_dx(levelSetI_);
    const auto& psiDyField = useOverride ? *opt.psiDyOverride : iGeometry.LS_dy(levelSetI_);

    const int nCellIndex = mesh.cellIndexNum();
    cell_.assign(nCellIndex, CellMetric{});
    active_.assign(nCellIndex, 0u);

    auto cellIsCut = [&](const Cell& c) -> bool {
        if (levelSetI_ < 0)
            return false;
        if (levelSetI_ >= static_cast<int>(iGeometry.levelSets.size()))
            return false;
        return c.isCut[levelSetI_];
    };

    // Active selection:
    // - If user provides predicate, mark those leaf cells.
    // - Else infer from FaceSegmentGraph: a cell is active if it has at least one segment.
    if (opt.isActive)
    {
        for (int ci = 0; ci < nCellIndex; ++ci)
        {
            const int cLabel = mesh.cellIndexToLabel[ci];
            const Cell& c = mesh.cells(cLabel);
            if (c.son[0] != -1) continue; // leaf only
            if (opt.isActive(mesh, c))
            {
                active_[ci] = 1u;
                // Defaults already set: alpha=1, etc.
                cell_[ci].isCut = cellIsCut(c);
            }
        }
    }
    else
    {
        for (int ci = 0; ci < nCellIndex; ++ci)
        {
            const auto range = faceGraph.segmentRangeForCellIndex(ci);
            if (range.second > range.first)
            {
                active_[ci] = 1u;
                const int cLabel = mesh.cellIndexToLabel[ci];
                const Cell& c = mesh.cells(cLabel);
                cell_[ci].isCut = cellIsCut(c);
            }
        }
    }

    // --------------------------------------------------------------------
    // alpha_K
    // --------------------------------------------------------------------
    // Prefer to reuse existing cut-cell utilities if present. In the current
    // QUANTUM infrastructure, we do not have a direct "cell volume fraction"
    // API exposed by ImmersedGeometry. We therefore approximate alpha_K by
    // sub-cell sampling of the level-set inside each active leaf cell.
    //
    // This is intentionally non-intrusive (no mesh/geometry modification) and
    // works for both uniform and quadtree meshes, including diagonal cut
    // patterns.

    auto isInside = [insideNeg](double psi) {
        return insideNeg ? (psi <= 0.0) : (psi >= 0.0);
    };

    auto sampleAlpha = [&](const Cell& c) -> double {
        const Vector2D xc = c.centerCoord();
        const Vector2D h = c.size();
        const double x0 = xc.x - 0.5 * h.x;
        const double y0 = xc.y - 0.5 * h.y;
        const double dx = h.x / static_cast<double>(alphaN);
        const double dy = h.y / static_cast<double>(alphaN);

        int insideCount = 0;
        for (int j = 0; j < alphaN; ++j)
        {
            for (int i = 0; i < alphaN; ++i)
            {
                const double xs = x0 + (static_cast<double>(i) + 0.5) * dx;
                const double ys = y0 + (static_cast<double>(j) + 0.5) * dy;
                const Vector2D x(xs, ys);

                // Bilinear interpolation on the *same* cell for speed.
                const double psi = Interpolation::biLinearInterpolation(mesh, psiField, c, x);
                if (isInside(psi)) ++insideCount;
            }
        }
        return static_cast<double>(insideCount) / static_cast<double>(alphaN * alphaN);
    };

    double alphaMin = 1.0;
    double alphaMax = 0.0;
    size_t nAlphaSampled = 0;
    size_t nActiveLeaf = 0;

    for (int ci = 0; ci < nCellIndex; ++ci)
    {
        if (!isActiveCellIndex(ci)) continue;

        const int cLabel = mesh.cellIndexToLabel[ci];
        const Cell& c = mesh.cells(cLabel);
        if (c.son[0] != -1) continue; // leaf only
        ++nActiveLeaf;

        // Fast path: if all four vertices are on the same side, alpha is 0/1.
        // Otherwise, fall back to sampling.
        const Vert& v0 = mesh.verts(c.cellToVert[0]);
        const Vert& v1 = mesh.verts(c.cellToVert[1]);
        const Vert& v2 = mesh.verts(c.cellToVert[2]);
        const Vert& v3 = mesh.verts(c.cellToVert[3]);
        const bool in0 = isInside(psiField.value[v0.index()]);
        const bool in1 = isInside(psiField.value[v1.index()]);
        const bool in2 = isInside(psiField.value[v2.index()]);
        const bool in3 = isInside(psiField.value[v3.index()]);

        double alpha = 1.0;
        if (in0 && in1 && in2 && in3)
        {
            alpha = 1.0;
        }
        else if ((!in0) && (!in1) && (!in2) && (!in3))
        {
            alpha = 0.0;
        }
        else
        {
            alpha = sampleAlpha(c);
            ++nAlphaSampled;
        }

        // Safety clamp
        if (!std::isfinite(alpha)) alpha = 0.0;
        if (alpha < 0.0) alpha = 0.0;
        if (alpha > 1.0) alpha = 1.0;

        cell_[ci].alpha = alpha;
        alphaMin = std::min(alphaMin, alpha);
        alphaMax = std::max(alphaMax, alpha);
    }

    // --------------------------------------------------------------------
    // A_gamma[K], d_gamma[K], n_gamma[K]
    // --------------------------------------------------------------------
    // A_gamma: interface segment length inside the cell (2D). We locate zero
    // crossings of psi along the *cell boundary leaf edges* (face segments)
    // and connect the two intersection points.
    //
    // d_gamma: distance from cell center to the interface along (approximate)
    // normal direction. We use |psi(xc)|/|grad psi(xc)| and clamp to a
    // reasonable range to avoid degeneracy.
    //
    // n_gamma: unit normal approximated from grad psi at cell center. We
    // orient it as outward-from-inside normal using insideIsNegative.

    // Reuse pre-bound references declared at the top of build().
    const auto& dxField = psiDxField;
    const auto& dyField = psiDyField;

    size_t nCutCellGamma = 0;
    size_t nRejectedFarCut = 0;
    double AgammaMin = std::numeric_limits<double>::infinity();
    double AgammaMax = 0.0;
    double AgammaSum = 0.0;
    double dgammaMin = std::numeric_limits<double>::infinity();
    double dgammaMax = 0.0;

    struct CachedLinePolyline
    {
        int globalLine = -1;
        std::vector<Vector2D> poly;
    };
    std::vector<CachedLinePolyline> cachedLinePolylines;
    if (opt.lineToBoundaryId && !opt.lineToBoundaryId->empty())
    {
        const std::vector<SignInt>* lsLinesPtr = nullptr;
        std::vector<SignInt> levelSetLines;
        const bool useLinesOverride = (opt.linesOverride && !opt.linesOverride->empty());
        if (useLinesOverride)
        {
            lsLinesPtr = opt.linesOverride;
        }
        else if (levelSetI_ >= 0 && levelSetI_ < static_cast<int>(iGeometry.levelSets.size()))
        {
            // lines() returns a value; keep it alive throughout tagging.
            levelSetLines = iGeometry.levelSets[levelSetI_].lines();
            lsLinesPtr = &levelSetLines;
        }

        if (lsLinesPtr)
        {
            const auto& lsLines = *lsLinesPtr;
            cachedLinePolylines.reserve(lsLines.size());
            std::unordered_set<int> seenGlobalLines;
            seenGlobalLines.reserve(lsLines.size());
            for (int localLineI = 0; localLineI < static_cast<int>(lsLines.size()); ++localLineI)
            {
                const SignInt si = lsLines[static_cast<std::size_t>(localLineI)];
                int globalLine = si.value;

                // Accept either 0-based or 1-based.
                if (globalLine < 0 || globalLine >= static_cast<int>(iGeometry.lines.size()))
                {
                    if (globalLine > 0 && (globalLine - 1) >= 0 && (globalLine - 1) < static_cast<int>(iGeometry.lines.size()))
                        globalLine -= 1;
                }
                if (globalLine < 0 || globalLine >= static_cast<int>(iGeometry.lines.size()))
                    continue;
                if (!seenGlobalLines.insert(globalLine).second)
                    continue;

                std::vector<Vector2D> poly = useLinesOverride
                    ? discretizeLineBySign(iGeometry.lines, si)
                    : iGeometry.discreteLines(levelSetI_, localLineI);
                if (poly.size() < 2)
                    continue;
                if (opt.transformLinePoint)
                {
                    for (auto& p : poly)
                        p = opt.transformLinePoint(globalLine, p);
                }

                CachedLinePolyline cached;
                cached.globalLine = globalLine;
                cached.poly = std::move(poly);
                cachedLinePolylines.push_back(std::move(cached));
            }
        }
    }


    for (int ci = 0; ci < nCellIndex; ++ci)
    {
        if (!isActiveCellIndex(ci)) continue;

        const int cLabel = mesh.cellIndexToLabel[ci];
        const Cell& c = mesh.cells(cLabel);
        if (c.son[0] != -1) continue; // leaf only

        // Collect intersection points on cell boundary by scanning its face segments.
        std::vector<Vector2D> xInts;
        xInts.reserve(4);

        const auto range = faceGraph.segmentRangeForCellIndex(ci);
        const Vector2D h = c.size();
        const double hmax = std::max(h.x, h.y);
        const double eps2 = sqr(1e-12 * (hmax > 0.0 ? hmax : 1.0));

        for (size_t s = range.first; s < range.second; ++s)
        {
            const FaceSegment& seg = faceGraph.segments().at(s);
            if (seg.edgeLabel < 0) continue;

            const Edge& e = seg.isHorizontal ? mesh.hEdges(seg.edgeLabel) : mesh.vEdges(seg.edgeLabel);
            // Use edge cut flag as a filter (still robust to sign-change below).
            const bool edgeCut = (levelSetI_ >= 0 ? e.isCut[levelSetI_] : false) || e.cut();
            if (!edgeCut)
            {
                // Still allow sign-change detection; skip only if both endpoints on same side.
            }

            const int vA = e.edgeToVert[0];
            const int vB = e.edgeToVert[1];
            if (vA < 0 || vB < 0) continue;
            const Vert& VA = mesh.verts(vA);
            const Vert& VB = mesh.verts(vB);
            const Vector2D a = VA.coord();
            const Vector2D b = VB.coord();

            const double psiA = psiField.value[VA.index()];
            const double psiB = psiField.value[VB.index()];

            // Quick skip: if no sign change and no zero endpoints and not flagged cut.
            if (!edgeCut && psiA != 0.0 && psiB != 0.0 && psiA * psiB > 0.0)
            {
                continue;
            }

            addZeroCrossingsOnEdge(xInts, a, b, psiA, psiB, eps2);
        }

        Vector2D p0, p1;
        double A_gamma = 0.0;
        if (farthestPair(xInts, p0, p1))
        {
            const double dx = p0.x - p1.x;
            const double dy = p0.y - p1.y;
            A_gamma = std::sqrt(dx * dx + dy * dy);
        }

        // Compute grad psi at cell center and derived normal/distance.
        const Vector2D xc = c.centerCoord();
        const double psiC = Interpolation::biLinearInterpolation(mesh, psiField, c, xc);
        const double gx = Interpolation::biLinearInterpolation(mesh, dxField, c, xc);
        const double gy = Interpolation::biLinearInterpolation(mesh, dyField, c, xc);
        const double gmag = std::sqrt(gx * gx + gy * gy);

        Vector2D n_gamma(0.0, 0.0);
        double d_gamma = 0.0;
        if (gmag > 0.0 && std::isfinite(gmag))
        {
            const double sgn = (opt.insideIsNegative ? 1.0 : -1.0);
            n_gamma = Vector2D(sgn * gx / gmag, sgn * gy / gmag);

            // Signed-distance approximation from level set.
            double d = std::fabs(psiC) / gmag;
            // Clamp to avoid tiny/huge values (robust for near-flat grad in coarse cells).
            const double dMin = embeddedDistanceFloor(hmax);
            const double dMax = 2.0 * (hmax > 0.0 ? hmax : 1.0);
            if (!std::isfinite(d)) d = dMax;
            if (d < dMin) d = dMin;
            if (d > dMax) d = dMax;
            d_gamma = d;
        }

        double nearestExplicitLineDist2 = std::numeric_limits<double>::infinity();
        if (!cachedLinePolylines.empty())
        {
            for (const auto& cached : cachedLinePolylines)
                nearestExplicitLineDist2 = std::min(nearestExplicitLineDist2, pointPolylineDist2(xc, cached.poly));
        }
        const bool alphaLooksFull = (cell_[ci].alpha >= (1.0 - 1.0e-12));
        const bool alphaLooksEmpty = (cell_[ci].alpha <= 1.0e-12);
        const bool hasExplicitGeometryDistance = std::isfinite(nearestExplicitLineDist2);
        const double nearestExplicitLineDist = hasExplicitGeometryDistance ? std::sqrt(nearestExplicitLineDist2) : 0.0;
        const double farCutRejectTol = 2.5 * (hmax > 0.0 ? hmax : 1.0);
        const bool suspiciousFarCut =
            hasExplicitGeometryDistance &&
            (nearestExplicitLineDist > farCutRejectTol) &&
            (alphaLooksFull || alphaLooksEmpty) &&
            ((A_gamma > 0.0) || cell_[ci].isCut);
        if (suspiciousFarCut)
        {
            A_gamma = 0.0;
            d_gamma = 0.0;
            n_gamma = Vector2D(0.0, 0.0);
            cell_[ci].isCut = false;
            cell_[ci].gammaBoundaryId = 0;
            cell_[ci].gammaLineIndex = -1;
            ++nRejectedFarCut;
        }

        // Only tag gamma metrics when the cell is cut (alpha in (0,1)) or we found a nonzero A_gamma.
        const bool isCutByAlpha = (cell_[ci].alpha > 0.0 && cell_[ci].alpha < 1.0);
        const bool isCutByFlag = cellIsCut(c);
        const bool hasGamma = (A_gamma > 0.0) || isCutByAlpha || isCutByFlag;

        if (hasGamma)
        {
            cell_[ci].A_gamma = A_gamma;
            cell_[ci].d_gamma = d_gamma;
            cell_[ci].n_gamma = n_gamma;
            cell_[ci].isCut = cell_[ci].isCut || (A_gamma > 0.0);

            // ----------------------------------------------------------------
            // Tag the cut-cell with a geometry boundaryId (gammaBoundaryId)
            // by finding the closest level-set interface line to a reference
            // point on the local EB segment.
            // ----------------------------------------------------------------
            if (opt.lineToBoundaryId && !opt.lineToBoundaryId->empty())
            {
                // Reference point on the local EB segment.
                Vector2D xRef = xc;
                if (A_gamma > 0.0)
                {
                    xRef = (p0 + p1) * 0.5;
                }
                else if (!xInts.empty())
                {
                    xRef = xInts.front();
                }
                double bestD2 = std::numeric_limits<double>::infinity();
                int bestGlobalLine = -1;

                for (const auto& cached : cachedLinePolylines)
                {
                    const double d2 = pointPolylineDist2(xRef, cached.poly);
                    if (d2 < bestD2)
                    {
                        bestD2 = d2;
                        bestGlobalLine = cached.globalLine;
                    }
                }

                if (bestGlobalLine >= 0 && bestGlobalLine < static_cast<int>(opt.lineToBoundaryId->size()))
                {
                    cell_[ci].gammaLineIndex = bestGlobalLine;
                    cell_[ci].gammaBoundaryId = (*opt.lineToBoundaryId)[static_cast<size_t>(bestGlobalLine)];
                }

            }

            ++nCutCellGamma;
            AgammaMin = std::min(AgammaMin, A_gamma);
            AgammaMax = std::max(AgammaMax, A_gamma);
            AgammaSum += A_gamma;
            dgammaMin = std::min(dgammaMin, d_gamma);
            dgammaMax = std::max(dgammaMax, d_gamma);
        }
        else
        {
            // keep defaults
        }
    }

    // --------------------------------------------------------------
    // Debug: cut-cell stats by embedded boundaryId
    // --------------------------------------------------------------
    if (opt.lineToBoundaryId && ConsoleControl::enabled(ConsoleControl::Level::kDebug))
    {
        std::map<int, size_t> byId;
        size_t nTagged = 0;
        size_t nUnknown = 0;
        for (int ci = 0; ci < nCellIndex; ++ci)
        {
            if (!isActiveCellIndex(ci))
                continue;
            const int cLabel = mesh.cellIndexToLabel[ci];
            const Cell& c = mesh.cells(cLabel);
            if (c.son[0] != -1)
                continue;

            const bool isCutByAlpha = (cell_[ci].alpha > 0.0 && cell_[ci].alpha < 1.0);
            const bool hasGamma = (cell_[ci].A_gamma > 0.0) || isCutByAlpha || cell_[ci].isCut;
            if (!hasGamma)
                continue;

            const int bid = cell_[ci].gammaBoundaryId;
            if (bid > 0)
            {
                byId[bid] += 1;
                ++nTagged;
            }
            else
            {
                ++nUnknown;
            }
        }

        std::cout << "[CutCellMetrics] gammaBoundaryId stats: tagged=" << nTagged
                  << ", unknown=" << nUnknown;
        if (!byId.empty())
        {
            std::cout << ", buckets={";
            bool first = true;
            for (const auto& kv : byId)
            {
                if (!first) std::cout << ", ";
                first = false;
                std::cout << kv.first << ":" << kv.second;
            }
            std::cout << "}";
        }
        std::cout << std::endl;
    }
    if (ConsoleControl::enabled(ConsoleControl::Level::kDebug) && nRejectedFarCut > 0)
    {
        std::cout << "[CutCellMetrics] rejected far-field cut tags=" << nRejectedFarCut << std::endl;
    }

    // Per-segment beta_f defaults to 1.0, aligned with faceGraph.segments().
    beta_.assign(faceGraph.segments().size(), 1.0);


    // Compute beta_f for each face segment.
    //
    // beta_f is the *open fraction* of the face segment with respect to the active
    // region ("inside" side of the level-set, controlled by Options::insideIsNegative).
    //
    // Important: Edge::cutLength() is measured from edgeToVert[0] to the intersection.
    // Therefore on a cut edge, the open length is:
    //   - cutLength       if vert0 is inside and vert1 is outside
    //   - (len-cutLength) if vert0 is outside and vert1 is inside
    //
    // On a 2:1 balanced quadtree, coarse-fine interfaces are represented by splitting
    // a coarse face into leaf edges (Edge::son[]). FaceSegmentGraph stores segments over
    // these leaf edges ("leaf edges"), so Edge::cutLength()/length() already corresponds
    // to the correct sub-segment length.
    const auto& segs = faceGraph.segments();
    double betaMin = 1.0;
    double betaMax = 1.0;
    size_t nCutSeg = 0;
    size_t nClamped = 0;
    size_t nDegenerate = 0;

    for (size_t s = 0; s < segs.size(); ++s)
    {
        const FaceSegment& seg = segs[s];
        if (seg.edgeLabel < 0)
        {
            continue;
        }

        const Edge& e = seg.isHorizontal ? mesh.hEdges(seg.edgeLabel) : mesh.vEdges(seg.edgeLabel);

        // Determine which side of this edge is "inside".
        const int vLbl0 = e.edgeToVert[0];
        const int vLbl1 = e.edgeToVert[1];
        const Vert& v0 = mesh.verts(vLbl0);
        const Vert& v1 = mesh.verts(vLbl1);

        const double psi0 = psiField.value[v0.index()];
        const double psi1 = psiField.value[v1.index()];
        const bool in0 = isInside(psi0);
        const bool in1 = isInside(psi1);

        // If the whole edge is inside/outside, beta is trivially 1/0.
        if (in0 && in1)
        {
            beta_[s] = 1.0;
        }
        else if (!in0 && !in1)
        {
            beta_[s] = 0.0;
        }
        else
        {
            // Mixed: must be cut (or nearly cut).
            const double len = e.length();
            if (!(len > 0.0))
            {
                beta_[s] = 0.0;
                ++nDegenerate;
            }
            else
            {
                double cutLen = e.cutLength();
                if (!std::isfinite(cutLen) || cutLen <= 0.0 || cutLen >= len)
                {
                    // Fallback: linear root from psi values.
                    const double denom = (psi0 - psi1);
                    double t = 0.5;
                    if (std::abs(denom) > 1e-30) t = psi0 / denom;
                    if (t < 0.0) t = 0.0;
                    if (t > 1.0) t = 1.0;
                    cutLen = len * t;
                    ++nClamped;
                }

                // Clamp cutLen to [0, len].
                if (cutLen < 0.0)
                {
                    cutLen = 0.0;
                    ++nClamped;
                }
                if (cutLen > len)
                {
                    cutLen = len;
                    ++nClamped;
                }

                const double insideLen = in0 ? cutLen : (len - cutLen);
                beta_[s] = insideLen / len;
                ++nCutSeg;
            }
        }

        // Final safety clamp.
        if (!std::isfinite(beta_[s]))
        {
            beta_[s] = 0.0;
            ++nDegenerate;
        }
        if (beta_[s] < 0.0)
        {
            beta_[s] = 0.0;
            ++nClamped;
        }
        if (beta_[s] > 1.0)
        {
            beta_[s] = 1.0;
            ++nClamped;
        }

        betaMin = std::min(betaMin, beta_[s]);
        betaMax = std::max(betaMax, beta_[s]);
    }

    // Summary
    size_t nActive = 0;
    for (unsigned char m : active_) if (m) ++nActive;

    if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
    {
        std::cout << "[CutCellMetrics] build levelSetI=" << levelSetI_
                  << ", activeCells=" << nActive
                  << ", activeLeaf=" << nActiveLeaf
                  << ", alphaN=" << alphaN
                  << ", alphaMin=" << alphaMin
                  << ", alphaMax=" << alphaMax
                  << ", alphaSampled=" << nAlphaSampled
                  << ", gammaCutCells=" << nCutCellGamma
                  << ", AgammaSum=" << AgammaSum
                  << ", AgammaMin=" << (std::isfinite(AgammaMin) ? AgammaMin : 0.0)
                  << ", AgammaMax=" << AgammaMax
                  << ", dGammaMin=" << (std::isfinite(dgammaMin) ? dgammaMin : 0.0)
                  << ", dGammaMax=" << dgammaMax
                  << ", segments=" << beta_.size()
                  << ", betaMin=" << betaMin
                  << ", betaMax=" << betaMax
                  << ", cutSeg=" << nCutSeg
                  << ", clamped=" << nClamped
                  << ", degenerate=" << nDegenerate
                  << "\n";
    }
}

void CutCellMetrics::dumpCellsCSV(
    const QuadtreeMesh& mesh,
    const std::string& filePath,
    const bool onlyActive) const
{
    std::ofstream os(filePath);
    if (!os)
    {
        std::cerr << "[CutCellMetrics] Failed to open file for writing: " << filePath << "\n";
        return;
    }

    os << "cellIndex,cellLabel,level,boundary,center_x,center_y,area,alpha,A_gamma,d_gamma,nx,ny,isCut,active\n";

    const int nCellIndex = mesh.cellIndexNum();
    for (int ci = 0; ci < nCellIndex; ++ci)
    {
        const bool act = isActiveCellIndex(ci);
        if (onlyActive && !act) continue;

        const int cLabel = mesh.cellIndexToLabel[ci];
        const Cell& c = mesh.cells(cLabel);
        const Vector2D xc = c.centerCoord();
        const double area = c.area();

        const CellMetric& m = cell_.at(ci);

        os << ci << "," << cLabel << "," << c.level() << "," << static_cast<int>(c.boundary()) << ","
           << xc.x << "," << xc.y << "," << area << ","
           << m.alpha << "," << m.A_gamma << "," << m.d_gamma << ","
           << m.n_gamma.x << "," << m.n_gamma.y << ","
           << (m.isCut ? 1 : 0) << "," << (act ? 1 : 0)
           << "\n";
    }

    if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
    if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
        std::cout << "[CutCellMetrics] dumpCellsCSV -> " << filePath << "\n";
}

void CutCellMetrics::dumpSegmentsCSV(
    const QuadtreeMesh& mesh,
    const FaceSegmentGraph& faceGraph,
    const std::string& filePath,
    const bool onlyBoundary) const
{
    (void)mesh;

    std::ofstream os(filePath);
    if (!os)
    {
        std::cerr << "[CutCellMetrics] Failed to open file for writing: " << filePath << "\n";
        return;
    }

    os << "segId,K,L,boundary,side,isHorizontal,edgeLabel,Af,df,nx,ny,faceCenter_x,faceCenter_y,edgeCut,edgeLen,edgeCutLen,beta,openAf\n";

    const auto& segs = faceGraph.segments();
    const size_t nSeg = segs.size();

    for (size_t s = 0; s < nSeg; ++s)
    {
        const auto& seg = segs[s];
        if (onlyBoundary && seg.boundary == Boundary::non) continue;

        const double beta = (s < beta_.size() ? beta_[s] : 1.0);

        const double openAf = beta * seg.Af;

        os << s << "," << seg.K << "," << seg.L << "," << static_cast<int>(seg.boundary) << ","
           << static_cast<int>(seg.side) << "," << (seg.isHorizontal ? 1 : 0) << "," << seg.edgeLabel << ","
           << seg.Af << "," << seg.df << "," << seg.n.x << "," << seg.n.y << ","
           << seg.faceCenter.x << "," << seg.faceCenter.y << ","
           << (seg.edgeCut ? 1 : 0) << "," << seg.edgeLength << "," << seg.edgeCutLength << ","
           << beta << "," << openAf
           << "\n";
    }

    if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
        std::cout << "[CutCellMetrics] dumpSegmentsCSV -> " << filePath << "\n";
}

} // namespace Poisson
} // namespace Solvers
