#include "Modules/Plasma/Boundary/PicEmbeddedGeometry.h"
#include "Modules/Plasma/Boundary/MovingEBController.h"
#include "Core/Base/Interpolation.h"
#include "Core/Geometry/ImmersedGeometry.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Units/UnitsContext.h"
#include "Core/Units/NormalizationContext.h"
#include <cmath>
#include <limits>
#include <unordered_set>

using Solvers::Poisson::CutCellMetrics;
namespace
{
static double pointSegmentDist2(const Vector2D& p, const Vector2D& a, const Vector2D& b)
{
    const double vx = b.x - a.x;
    const double vy = b.y - a.y;
    const double wx = p.x - a.x;
    const double wy = p.y - a.y;
    const double c1 = vx * wx + vy * wy;
    if (c1 <= 0.0)
    {
        const double dx = p.x - a.x;
        const double dy = p.y - a.y;
        return dx * dx + dy * dy;
    }
    const double c2 = vx * vx + vy * vy;
    if (c2 <= c1)
    {
        const double dx = p.x - b.x;
        const double dy = p.y - b.y;
        return dx * dx + dy * dy;
    }
    const double t = c1 / c2;
    const double projx = a.x + t * vx;
    const double projy = a.y + t * vy;
    const double dx = p.x - projx;
    const double dy = p.y - projy;
    return dx * dx + dy * dy;
}

static double pointPolylineDist2(const Vector2D& p, const std::vector<Vector2D>& poly)
{
    if (poly.size() < 2)
        return std::numeric_limits<double>::infinity();
    double best = std::numeric_limits<double>::infinity();
    for (std::size_t i = 1; i < poly.size(); ++i)
    {
        const double d2 = pointSegmentDist2(p, poly[i - 1], poly[i]);
        if (d2 < best)
            best = d2;
    }
    return best;
}

}

PicEmbeddedGeometry::PicEmbeddedGeometry(Mode mode) : mode_(mode) {}
PicEmbeddedGeometry::~PicEmbeddedGeometry() = default;

void PicEmbeddedGeometry::resolveSelection(const DictionaryGeometry& dictGeometry, int nLS)
{
    clearComposite();
    cutMetrics = Metrics{};
    cutMetricsReady = false;
    cutAlphaByCellIndex.clear();
    embeddedBoundaryIdByCellIndex.clear();
    embeddedLevelSetIndices.clear();
    embeddedBoundaryIds.clear();
    embeddedLevelSetIndex = -1;
    embeddedGeomReady = embeddedAbsorb = false;
    embeddedInsideIsNegative = true;
    const char* prefix = mode_ == Mode::Explicit ? "[ExplicitPIC]" : "[ImplicitPIC]";
    auto appendUnique = [](std::vector<int>& out, int v)
    {
        if (v < 0)
            return;
        if (std::find(out.begin(), out.end(), v) == out.end())
            out.push_back(v);
    };

    auto addLevelSet = [&](int idx, const std::string& reason)
    {
        if (idx < 0 || idx >= nLS)
        {
            if (mode_ == Mode::Explicit)
            std::cerr << "[ExplicitPIC][EB] WARNING: " << reason << " levelSetIndex="
                      << idx << " out of range (0.." << (nLS - 1) << "); skipping.\n";
            return;
        }
        appendUnique(embeddedLevelSetIndices, idx);
    };

    embeddedBoundaryIds = dictGeometry.immersedBoundaryIds();
    for (int bid : embeddedBoundaryIds)
    {
        auto it = dictGeometry.boundaryById.find(bid);
        if (it == dictGeometry.boundaryById.end())
            continue;
        const auto& ent = it->second;
        if (ent.levelSetIndex >= 0)
            addLevelSet(ent.levelSetIndex, "geometry boundary");
    }

    if (!embeddedBoundaryIds.empty() && !embeddedLevelSetIndices.empty())
    {
        std::vector<int> filteredBoundaryIds;
        filteredBoundaryIds.reserve(embeddedBoundaryIds.size());
        for (int bid : embeddedBoundaryIds)
        {
            auto it = dictGeometry.boundaryById.find(bid);
            if (it == dictGeometry.boundaryById.end())
            {
                continue;
            }
            const auto& ent = it->second;
            if (ent.ignored)
                continue;
            if (std::find(embeddedLevelSetIndices.begin(), embeddedLevelSetIndices.end(), ent.levelSetIndex) ==
                embeddedLevelSetIndices.end())
            {
                continue;
            }
            filteredBoundaryIds.push_back(bid);
        }
        embeddedBoundaryIds = std::move(filteredBoundaryIds);
    }

    if (embeddedLevelSetIndices.empty())
    {
        if (embeddedBoundaryIds.empty())
        {
            embeddedLevelSetIndex = 0;
            addLevelSet(embeddedLevelSetIndex, "fallback levelSetIndex");
        }
    }
    else
    {
        embeddedLevelSetIndex = embeddedLevelSetIndices.front();
    }

    if (embeddedLevelSetIndices.empty())
    {
        if (mode_ == Mode::Explicit)
        std::cerr << "[ExplicitPIC][EB] WARNING: no valid embedded level sets resolved; disabling EB.\n";
        embeddedGeomReady = false;
    }
    else
    {
        embeddedGeomReady = true;
        embeddedAbsorb = true;

        bool inferredInsideIsNegative = true;
        bool inferredSign = false;
        if (!embeddedBoundaryIds.empty())
            inferredSign = dictGeometry.inferEmbeddedInsideIsNegative(inferredInsideIsNegative);
        if (!inferredSign)
        {
            for (int lsI : embeddedLevelSetIndices)
            {
                bool localInsideIsNegative = true;
                if (!dictGeometry.inferLevelSetInsideIsNegative(lsI, localInsideIsNegative))
                    continue;
                if (!inferredSign)
                {
                    inferredInsideIsNegative = localInsideIsNegative;
                    inferredSign = true;
                    continue;
                }
                if (inferredInsideIsNegative != localInsideIsNegative)
                {
                    inferredSign = false;
                    break;
                }
            }
        }
        if (!inferredSign)
        {
            std::cerr << prefix << "[EB] WARNING: failed to infer insideIsNegative from geometry.dict. "
                      << "Using runtime default true.\n";
        }
        else
        {
            embeddedInsideIsNegative = inferredInsideIsNegative;
        }
    }
}

void PicEmbeddedGeometry::clearComposite()
{
    embeddedUseComposite = false;
    embeddedLevelSetVert = nullptr;
    embeddedPsiComposite.reset();
    embeddedDxComposite.reset();
    embeddedDyComposite.reset();
    embeddedCompositeLines.clear();

}

void PicEmbeddedGeometry::rebuildComposite(const QuadtreeMesh& mesh, const Geometry* iGeometry)
{
    clearComposite();

    if (!embeddedGeomReady || !iGeometry)
        return;
    if (embeddedLevelSetIndices.empty())
        return;

    if (embeddedLevelSetIndices.size() == 1)
    {
        embeddedLevelSetIndex = embeddedLevelSetIndices.front();
        if (embeddedLevelSetIndex >= 0 &&
            embeddedLevelSetIndex < static_cast<int>(iGeometry->levelSets.size()))
        {
            embeddedLevelSetVert = &iGeometry->levelSets[embeddedLevelSetIndex];
        }
        return;
    }

    embeddedUseComposite = true;

    embeddedPsiComposite = std::make_unique<VertField<double>>(mesh, 0.0);
    embeddedDxComposite = std::make_unique<VertField<double>>(mesh, 0.0);
    embeddedDyComposite = std::make_unique<VertField<double>>(mesh, 0.0);

    // Union of interface lines (used for gammaBoundaryId mapping).
    std::unordered_set<int> seenLines;
    if (mode_ == Mode::Explicit) seenLines.reserve(embeddedLevelSetIndices.size() * 64);
    for (int lsI : embeddedLevelSetIndices)
    {
        if (lsI < 0 || lsI >= static_cast<int>(iGeometry->levelSets.size()))
            continue;
        for (const auto& si : iGeometry->levelSets[lsI].lines())
        {
            int globalLine = si.value;
            if (globalLine < 0 || globalLine >= static_cast<int>(iGeometry->lines.size()))
            {
                if (globalLine > 0 && (globalLine - 1) >= 0 &&
                    (globalLine - 1) < static_cast<int>(iGeometry->lines.size()))
                    globalLine -= 1;
            }
            if (globalLine < 0 || globalLine >= static_cast<int>(iGeometry->lines.size()))
                continue;
            if (seenLines.insert(globalLine).second)
                embeddedCompositeLines.push_back(si);
        }
    }

    const int nVert = mesh.vertIndexNum();
    for (int vi = 0; vi < nVert; ++vi)
    {
        int bestIdx = -1;
        double bestPsi = 0.0;
        for (int lsI : embeddedLevelSetIndices)
        {
            if (lsI < 0 || lsI >= static_cast<int>(iGeometry->levelSets.size()))
                continue;
            const double psi = iGeometry->levelSets[lsI].value[vi];
            if (!std::isfinite(psi))
                continue;
            if (bestIdx < 0 ||
                (embeddedInsideIsNegative ? (psi < bestPsi) : (psi > bestPsi)))
            {
                bestPsi = psi;
                bestIdx = lsI;
            }
        }

        if (bestIdx >= 0)
        {
            embeddedPsiComposite->value[vi] = bestPsi;
            embeddedDxComposite->value[vi] = iGeometry->LS_dx(bestIdx).value[vi];
            embeddedDyComposite->value[vi] = iGeometry->LS_dy(bestIdx).value[vi];
        }
    }

    embeddedLevelSetVert = embeddedPsiComposite.get();
}

double PicEmbeddedGeometry::surfaceAreaM2(const Geometry* iGeometry, const DictionaryGeometry& dictGeometry,
    const std::vector<int>& lineToBoundaryId, const UnitsContext* units, const NormalizationContext* normalization) const
{
    if (!embeddedGeomReady || embeddedBoundaryIds.empty() || normalization == nullptr)
        return 0.0;

    const double slabThickness_m = normalization->config().slabThickness_m;
    if (!(slabThickness_m > 0.0))
        return 0.0;

    const std::vector<Line>* lines = nullptr;
    if (iGeometry && !iGeometry->lines.empty())
        lines = &iGeometry->lines;
    else if (!dictGeometry.lines.empty())
        lines = &dictGeometry.lines;
    if (lines == nullptr || lineToBoundaryId.empty())
        return 0.0;

    const std::size_t nLines = std::min(lines->size(), lineToBoundaryId.size());
    double totalLength_m = 0.0;
    for (std::size_t li = 0; li < nLines; ++li)
    {
        const int boundaryId = lineToBoundaryId[li];
        if (boundaryId <= 0)
            continue;
        if (std::find(embeddedBoundaryIds.begin(), embeddedBoundaryIds.end(), boundaryId) == embeddedBoundaryIds.end())
            continue;

        const Line& line = (*lines)[li];
        const double dx = line.points[1].x - line.points[0].x;
        const double dy = line.points[1].y - line.points[0].y;
        const double lenCode = std::hypot(dx, dy);
        totalLength_m += (units != nullptr) ? units->LengthCodeToSI(lenCode) : lenCode;
    }

    return std::max(0.0, totalLength_m * slabThickness_m);
}

bool PicEmbeddedGeometry::activeCell(const QuadtreeMesh& m, const Geometry* iGeometry, const Cell& c) const
{
    if (!embeddedGeomReady || !iGeometry)
        return true;
    if (embeddedLevelSetIndices.empty())
        return true;
    const double hMax = std::max(c.size().x, c.size().y);
    // psi is a signed-distance-like field in code units; require a clear
    // inside margin to avoid classifying cells by near-zero sign noise.
    const double solidMargin = std::max(1e-12, 0.25 * hMax);
    for (int lsI : embeddedLevelSetIndices)
    {
        if (lsI < 0 || lsI >= static_cast<int>(iGeometry->levelSets.size()))
            continue;
        bool allInside = true;
        double maxPsi = -std::numeric_limits<double>::infinity();
        for (int i = 0; i < 4; ++i)
        {
            const Vert& v = m.verts(c.cellToVert[i]);
            const double psi = iGeometry->levelSets[lsI].value[v.index()];
            maxPsi = std::max(maxPsi, psi);
            const bool inside = embeddedInsideIsNegative ? (psi < 0.0) : (psi > 0.0);
            if (!inside)
            {
                allInside = false;
                break;
            }
        }
        if (allInside)
        {
            const bool confidentlyInside = embeddedInsideIsNegative
                                              ? (maxPsi <= -solidMargin)
                                              : (maxPsi >= solidMargin);
            if (mode_ == Mode::Explicit || confidentlyInside)
                return false;
        }
    }
    return true;
}

void PicEmbeddedGeometry::rebuildMetrics(const QuadtreeMesh& mesh, Geometry* iGeometry,
    const FaceSegmentGraph& faceGraph, const std::vector<int>& lineToBoundaryId, const MovingEBController* movingEb)
{
    cutMetricsReady = false;
    cutAlphaByCellIndex.clear();
    if (mode_ == Mode::Implicit) embeddedBoundaryIdByCellIndex.clear();

    if (!embeddedGeomReady || !iGeometry)
        return;

    CutCellMetrics::BuildOptions opt;
    opt.alphaSampleN = alphaSampleCount;
    // Cut metrics should represent the plasma (fluid) region.
    // If solid is negative, fluid is positive => insideIsNegative=false.
    opt.insideIsNegative = !embeddedInsideIsNegative;

    if (!lineToBoundaryId.empty())
        opt.lineToBoundaryId = &lineToBoundaryId;

    if (mode_ == Mode::Implicit && movingEb && movingEb->enabled() && !lineToBoundaryId.empty() &&
        !movingEb->linesInCurrentGeometryFrame())
    {
        opt.transformLinePoint =
            [movingEbPtr = movingEb, &lineToBoundaryId](int globalLine, const Vector2D& p) -> Vector2D {
                if (movingEbPtr == nullptr || globalLine < 0 || globalLine >= static_cast<int>(lineToBoundaryId.size()))
                    return p;
                const int boundaryId = lineToBoundaryId[static_cast<std::size_t>(globalLine)];
                if (boundaryId <= 0)
                    return p;
                return p + movingEbPtr->boundaryRigidShiftCode(boundaryId);
            };
    }

    if (embeddedUseComposite)
    {
        opt.levelSetI = -1;
        opt.psiOverride = embeddedPsiComposite.get();
        opt.psiDxOverride = embeddedDxComposite.get();
        opt.psiDyOverride = embeddedDyComposite.get();
        if (!embeddedCompositeLines.empty())
            opt.linesOverride = &embeddedCompositeLines;
    }
    else
    {
        opt.levelSetI = embeddedLevelSetIndex;
    }

    cutMetrics.build(mesh, *iGeometry, faceGraph, opt);
    cutMetricsReady = true;

    const auto& cm = cutMetrics.cellMetricsByIndex();
    if (cm.size() == static_cast<size_t>(mesh.cellIndexNum()))
    {
        cutAlphaByCellIndex.resize(mesh.cellIndexNum(), 1.0);
        if (mode_ == Mode::Explicit) embeddedBoundaryIdByCellIndex.assign(mesh.cellIndexNum(), 0);
        else embeddedBoundaryIdByCellIndex.resize(mesh.cellIndexNum(), 0);
        for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
        {
            double alpha = cm[static_cast<size_t>(ci)].alpha;
            if (!cutMetrics.isActiveCellIndex(ci))
                alpha = 0.0;
            cutAlphaByCellIndex[ci] = alpha;
            embeddedBoundaryIdByCellIndex[ci] = cm[static_cast<size_t>(ci)].gammaBoundaryId;
        }
    }
    else
    {
        embeddedBoundaryIdByCellIndex.clear();
    }
}

int PicEmbeddedGeometry::boundaryIdAt(const QuadtreeMesh& mesh, Geometry* iGeometry,
    const std::vector<int>& lineToBoundaryId, const Vector2D& pos, int cellLabelHint) const
{
    auto posInCell = [&](int cellLabel) -> bool {
        if (cellLabel < 0)
            return false;
        const Cell& c = mesh.cells(cellLabel);
        const Vector2D cc = c.centerCoord();
        const Vector2D s = c.size();
        const double h = std::min(s.x, s.y);
        const double tol = 1e-12 + 1e-9 * h;
        if (pos.x < cc.x - 0.5 * s.x - tol || pos.x > cc.x + 0.5 * s.x + tol)
            return false;
        if (pos.y < cc.y - 0.5 * s.y - tol || pos.y > cc.y + 0.5 * s.y + tol)
            return false;
        return true;
    };

    if (cutMetricsReady && !embeddedBoundaryIdByCellIndex.empty())
    {
        int cellLabel = cellLabelHint;
        if (cellLabel < 0 || !posInCell(cellLabel))
        {
            cellLabel = mesh.findCellLabelFromHint(cellLabelHint, pos);
            if (cellLabel < 0 || !posInCell(cellLabel))
                cellLabel = mesh.findCellLabel(pos);
        }
        if (cellLabel >= 0)
        {
            const int cellIndex = mesh.cells(cellLabel).index();
            if (cellIndex >= 0 && cellIndex < static_cast<int>(embeddedBoundaryIdByCellIndex.size()))
            {
                const int bid = embeddedBoundaryIdByCellIndex[static_cast<size_t>(cellIndex)];
                if (bid > 0)
                    return bid;
            }
        }
    }

    if (!iGeometry || lineToBoundaryId.empty())
        return 0;

    double bestD2 = std::numeric_limits<double>::infinity();
    int bestGlobalLine = -1;
    std::vector<int> lsList = embeddedLevelSetIndices;
    if (lsList.empty())
        lsList.push_back(embeddedLevelSetIndex);
    for (int lsI : lsList)
    {
        if (lsI < 0 || lsI >= static_cast<int>(iGeometry->levelSets.size()))
            continue;
        const auto& lsLines = iGeometry->levelSets[lsI].lines();
        for (int localLineI = 0; localLineI < static_cast<int>(lsLines.size()); ++localLineI)
        {
            const SignInt si = lsLines[static_cast<size_t>(localLineI)];
            int globalLine = si.value;

            // Accept either 0-based or 1-based.
            if (globalLine < 0 || globalLine >= static_cast<int>(iGeometry->lines.size()))
            {
                if (globalLine > 0 && (globalLine - 1) >= 0 &&
                    (globalLine - 1) < static_cast<int>(iGeometry->lines.size()))
                    globalLine -= 1;
            }
            if (globalLine < 0 || globalLine >= static_cast<int>(iGeometry->lines.size()))
                continue;

            const std::vector<Vector2D> poly = iGeometry->discreteLines(lsI, localLineI);
            if (poly.size() < 2)
                continue;
            const double d2 = pointPolylineDist2(pos, poly);
            if (d2 < bestD2)
            {
                bestD2 = d2;
                bestGlobalLine = globalLine;
            }
        }
    }

    if (bestGlobalLine >= 0 && bestGlobalLine < static_cast<int>(lineToBoundaryId.size()))
        return lineToBoundaryId[static_cast<size_t>(bestGlobalLine)];
    return 0;
}

bool PicEmbeddedGeometry::containsPoint(const QuadtreeMesh& mesh, const Geometry* iGeometry,
    const Vector2D& pos, int cellLabelHint) const
{
    if (mode_ == Mode::Implicit)
    {
        if (!embeddedGeomReady || !embeddedAbsorb || !iGeometry)
            return false;
        constexpr double epsPsi = 1e-12;
        for (const int lsI : embeddedLevelSetIndices)
        {
            if (lsI < 0 || lsI >= static_cast<int>(iGeometry->levelSets.size()))
                continue;
            int cellLabel = cellLabelHint;
            if (cellLabel >= 0)
                cellLabel = mesh.findCellLabelFromHint(cellLabel, pos);
            if (cellLabel < 0)
                cellLabel = mesh.findCellLabel(pos);
            if (cellLabel < 0)
                continue;
            const Cell& c = mesh.cells(cellLabel);
            const double psi = Interpolation::biLinearInterpolation(mesh, iGeometry->levelSets[lsI], c, pos);
            if (embeddedInsideIsNegative ? (psi <= epsPsi) : (psi >= -epsPsi))
                return true;
        }
        return false;
    }
    if (!embeddedGeomReady || !iGeometry)
        return false;

    int cLabel = cellLabelHint;
    if (cLabel >= 0)
        cLabel = mesh.findCellLabelFromHint(cLabel, pos);
    if (cLabel < 0)
        cLabel = mesh.findCellLabel(pos);
    if (cLabel < 0)
        return false;

    const Cell& c = mesh.cells(cLabel);
    const int cIndex = c.index();
    const bool hasCutAlpha = cutMetricsReady &&
                             (cutAlphaByCellIndex.size() == static_cast<std::size_t>(mesh.cellIndexNum()));
    if (hasCutAlpha && cIndex >= 0)
    {
        double alpha = cutAlphaByCellIndex[static_cast<std::size_t>(cIndex)];
        if (!std::isfinite(alpha))
            alpha = 1.0;
        if (alpha < 0.0)
            alpha = 0.0;
        else if (alpha > 1.0)
            alpha = 1.0;

        constexpr double kAlphaEps = 1e-12;
        if (alpha <= kAlphaEps)
            return true;   // fully solid
        if (alpha >= (1.0 - kAlphaEps))
            return false;  // fully fluid
    }

    const double h = std::min(c.size().x, c.size().y);
    const double epsPsi = std::max(1e-12, 1e-3 * h);

    if (!embeddedLevelSetIndices.empty())
    {
        for (int lsI : embeddedLevelSetIndices)
        {
            if (lsI < 0 || lsI >= static_cast<int>(iGeometry->levelSets.size()))
                continue;
            const double psi = Interpolation::biLinearInterpolation(mesh, iGeometry->levelSets[lsI], c, pos);
            const bool inside = embeddedInsideIsNegative ? (psi <= epsPsi) : (psi >= -epsPsi);
            if (inside)
                return true;
        }
        return false;
    }

    const int lsI = embeddedLevelSetIndex;
    if (lsI < 0 || lsI >= static_cast<int>(iGeometry->levelSets.size()))
        return false;
    const double psi = Interpolation::biLinearInterpolation(mesh, iGeometry->levelSets[lsI], c, pos);
    return embeddedInsideIsNegative ? (psi <= epsPsi) : (psi >= -epsPsi);
}
