#include "Solvers/Poisson/PoissonFVSolverCellCentered.h"

#include <cmath>
#include <chrono>
#include <algorithm>
#include <iostream>
#include <limits>
#include <queue>
#include <unordered_map>
#include <vector>

#include "Core/IO/ConsoleControl.h"

#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>

#ifdef EIGEN_SUPERLU_SUPPORT
#include <Eigen/SuperLUSupport>
#endif

#include "Core/Mesh/Cell.h"
#include "Core/Mesh/PeriodicFaceUtils.h"
#include "Core/Mesh/QuadtreeMesh.h"

#include "Solvers/Poisson/CutCellMetrics.h"

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
    default:              return bc.left; // fallback (should not happen)
    }
}

inline const char* bcTypeName(PicBCType t)
{
    switch (t)
    {
    case PicBCType::Dirichlet: return "Dirichlet";
    case PicBCType::Neumann:   return "Neumann";
    case PicBCType::Robin:     return "Robin";
    case PicBCType::Periodic:  return "Periodic";
    default:                   return "Unknown";
    }
}

inline bool isDirichletLike(PicBCType t)
{
    return (t == PicBCType::Dirichlet) || (t == PicBCType::Robin);
}

inline double safeCellArea(const QuadtreeMesh& mesh, const int cellIndex)
{
    if (cellIndex < 0 || cellIndex >= mesh.cellIndexNum())
        return 0.0;
    const int lab = mesh.cellIndexToLabel[cellIndex];
    return mesh.cells(lab).area();
}

inline double stabilizedEmbeddedDistance(const double dGamma, const double cellArea)
{
    const double h = (std::isfinite(cellArea) && cellArea > 0.0) ? std::sqrt(cellArea) : 1.0;
    const double floor = std::max(1e-12, 2.0e-2 * h);
    if (!std::isfinite(dGamma))
        return floor;
    return std::max(dGamma, floor);
}

inline bool sameBoundarySide(const PicBoundarySide& a, const PicBoundarySide& b)
{
    return (a.type == b.type) && (a.value == b.value);
}

inline bool sameBoundarySideStructure(const PicBoundarySide& a, const PicBoundarySide& b)
{
    return (a.type == b.type);
}

inline bool sameDomainBC(const PicDomainBC& a, const PicDomainBC& b)
{
    return sameBoundarySide(a.left, b.left) &&
           sameBoundarySide(a.right, b.right) &&
           sameBoundarySide(a.bottom, b.bottom) &&
           sameBoundarySide(a.top, b.top);
}

inline bool sameDomainBCStructure(const PicDomainBC& a, const PicDomainBC& b)
{
    return sameBoundarySideStructure(a.left, b.left) &&
           sameBoundarySideStructure(a.right, b.right) &&
           sameBoundarySideStructure(a.bottom, b.bottom) &&
           sameBoundarySideStructure(a.top, b.top);
}

inline const char* boundaryToString(const Boundary b)
{
    switch (b)
    {
    case Boundary::west:  return "west";
    case Boundary::south: return "south";
    case Boundary::east:  return "east";
    case Boundary::north: return "north";
    default: return "non";
    }
}

size_t stabilizeFloatingComponents(const std::vector<std::vector<int>>& neighbors,
                                   std::vector<unsigned char>& anchoredRow,
                                   const std::vector<double>& volOfRow,
                                   Eigen::VectorXd& bVec,
                                   std::vector<Eigen::Triplet<double>>& trips,
                                   std::vector<double>& diag,
                                   std::vector<int>* pinnedRows,
                                   std::vector<std::vector<int>>* floatingComponents)
{
    const int N = static_cast<int>(neighbors.size());
    std::vector<unsigned char> visited(static_cast<size_t>(N), 0u);
    std::vector<int> componentRows;
    componentRows.reserve(static_cast<size_t>(N));
    size_t floatingComponentCount = 0u;

    if (floatingComponents)
        floatingComponents->clear();

    for (int seed = 0; seed < N; ++seed)
    {
        if (visited[static_cast<size_t>(seed)] != 0u)
            continue;

        componentRows.clear();
        bool hasAnchor = false;
        std::queue<int> q;
        q.push(seed);
        visited[static_cast<size_t>(seed)] = 1u;

        while (!q.empty())
        {
            const int row = q.front();
            q.pop();
            componentRows.push_back(row);
            hasAnchor = hasAnchor || (anchoredRow[static_cast<size_t>(row)] != 0u);

            for (const int nbr : neighbors[static_cast<size_t>(row)])
            {
                if (nbr < 0 || nbr >= N)
                    continue;
                if (visited[static_cast<size_t>(nbr)] != 0u)
                    continue;
                visited[static_cast<size_t>(nbr)] = 1u;
                q.push(nbr);
            }
        }

        if (hasAnchor)
            continue;

        ++floatingComponentCount;

        double totalB = 0.0;
        double totalV = 0.0;
        for (const int row : componentRows)
        {
            totalB += bVec[row];
            totalV += volOfRow[static_cast<size_t>(row)];
        }
        if (totalV > 0.0)
        {
            const double meanBperV = totalB / totalV;
            for (const int row : componentRows)
                bVec[row] -= meanBperV * volOfRow[static_cast<size_t>(row)];
        }

        int pinRow = componentRows.front();
        for (const int row : componentRows)
        {
            if (volOfRow[static_cast<size_t>(row)] > 0.0)
            {
                pinRow = row;
                break;
            }
        }

        trips.emplace_back(pinRow, pinRow, 1.0);
        diag[static_cast<size_t>(pinRow)] += 1.0;
        bVec[pinRow] = 0.0;
        anchoredRow[static_cast<size_t>(pinRow)] = 1u;
        if (pinnedRows)
            pinnedRows->push_back(pinRow);
        if (floatingComponents)
            floatingComponents->push_back(componentRows);
    }

    return floatingComponentCount;
}

}

inline bool sameEmbeddedBC(const PicEmbeddedBC& a, const PicEmbeddedBC& b)
{
    return (a.enabled == b.enabled) &&
           (a.type == b.type) &&
           (a.a == b.a) &&
           (a.b == b.b) &&
           (a.g == b.g);
}

inline bool sameEmbeddedBCStructure(const PicEmbeddedBC& a, const PicEmbeddedBC& b)
{
    return (a.enabled == b.enabled) &&
           (a.type == b.type) &&
           (a.a == b.a) &&
           (a.b == b.b);
}

inline bool sameEmbeddedBCMapStructure(const PicEmbeddedBCMap& a, const PicEmbeddedBCMap& b)
{
    if (a.size() != b.size())
        return false;
    for (const auto& kv : a)
    {
        auto it = b.find(kv.first);
        if (it == b.end())
            return false;
        if (!sameEmbeddedBCStructure(kv.second, it->second))
            return false;
    }
    return true;
}

inline bool sameEmbeddedBCMap(const PicEmbeddedBCMap& a, const PicEmbeddedBCMap& b)
{
    if (a.size() != b.size())
        return false;
    for (const auto& kv : a)
    {
        auto it = b.find(kv.first);
        if (it == b.end())
            return false;
        if (!sameEmbeddedBC(kv.second, it->second))
            return false;
    }
    return true;
}

bool rebuildBCConst(const QuadtreeMesh& mesh,
                    const FaceSegmentGraph& graph,
                    const PicDomainBC& bc,
                    const Solvers::Poisson::CutCellMetrics* cutMetrics,
                    const PicEmbeddedBC* embeddedBC,
                    const PicEmbeddedBCMap* embeddedBCById,
                    const PoissonFVSolverCellCentered::Options& opt,
                    const std::vector<int>& rowOfCellIndex,
                    const std::vector<int>& cellIndexOfRow,
                    const std::vector<double>& volOfRow,
                    Eigen::VectorXd& bConst,
                    std::string* errMsg)
{
    const bool useCutMetrics =
        (cutMetrics != nullptr) &&
        (cutMetrics->betaBySegment().size() == graph.segments().size()) &&
        (cutMetrics->cellMetricsByIndex().size() == static_cast<size_t>(mesh.cellIndexNum()));

    const bool useEmbeddedBCMap =
        useCutMetrics &&
        (embeddedBCById != nullptr) &&
        (!embeddedBCById->empty());
    const bool useEmbeddedBCFallback =
        useCutMetrics &&
        (embeddedBC != nullptr) &&
        embeddedBC->enabled;
    const bool useEmbeddedBC = useEmbeddedBCMap || useEmbeddedBCFallback;

    const int N = static_cast<int>(cellIndexOfRow.size());
    bConst.resize(N);
    bConst.setZero();

    const auto& segs = graph.segments();
    for (size_t sid = 0; sid < segs.size(); ++sid)
    {
        const auto& seg = segs[sid];
        if (seg.boundary == Boundary::non)
            continue;

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
            if (beta < 0.0)
                beta = 0.0;
            else if (beta > 1.0)
                beta = 1.0;
        }

        const double AfOpen = beta * seg.Af;
        if (!(seg.df > 0.0) || !(AfOpen >= 0.0))
            continue;
        if (AfOpen == 0.0)
            continue;

        const PicBoundarySide& s = bcSide(bc, seg.boundary);
        PicBCType t = s.type;
        if (t == PicBCType::Robin)
            t = PicBCType::Dirichlet;

        if (t == PicBCType::Periodic)
            continue;

        if (t == PicBCType::Dirichlet)
        {
            const double c = AfOpen / seg.df;
            bConst[i] += c * s.value;
        }
        else if (t == PicBCType::Neumann)
        {
            bConst[i] += s.value * AfOpen;
        }
    }

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
            {
                if (opt.embeddedBCMissingPolicy == PoissonFVSolverCellCentered::Options::EmbeddedBCMissingPolicy::Error)
                {
                    if (errMsg)
                    {
                        *errMsg = "Missing embedded BC for boundaryId " + std::to_string(cm.gammaBoundaryId) +
                                  " (no map entry and no fallback).";
                    }
                    return false;
                }
                continue;
            }

            const double a = eb->a;
            const double b = eb->b;
            const double g = eb->g;
            const PicBCType t = eb->type;
            const double A_gamma = cm.A_gamma;

            if (t == PicBCType::Robin)
            {
                if (std::abs(b) > 0.0)
                {
                    const double addRhs = (g / b) * A_gamma;
                    if (!std::isfinite(addRhs))
                        continue;
                    bConst[i] += addRhs;
                }
                else
                {
                    if (std::abs(a) <= 0.0)
                        continue;
                    const double phiD = g / a;
                    const double d = stabilizedEmbeddedDistance(cm.d_gamma,
                                                                volOfRow[static_cast<size_t>(i)]);

                    const double c = A_gamma / d;
                    if (!std::isfinite(c) || c <= 0.0)
                        continue;
                    bConst[i] += c * phiD;
                }
            }
            else if (t == PicBCType::Dirichlet)
            {
                const double phiD = g;
                const double d = stabilizedEmbeddedDistance(cm.d_gamma,
                                                            volOfRow[static_cast<size_t>(i)]);
                const double c = A_gamma / d;
                if (!std::isfinite(c) || c <= 0.0)
                    continue;
                bConst[i] += c * phiD;
            }
            else if (t == PicBCType::Neumann)
            {
                const double addRhs = g * A_gamma;
                if (!std::isfinite(addRhs))
                    continue;
                bConst[i] += addRhs;
            }
        }
    }
    return true;
}


// ============================================================================
// Two-stage API + factorization cache (setup + solveWithCache)
// ============================================================================

struct PoissonFVSolverCellCentered::FactorCache
{
    // Cache signature (conservative) used to decide cache hit/miss.
    int meshCellIndexNum = 0;
    size_t graphSegmentCount = 0;
    PicDomainBC bc;
    bool requestSuperLU = true;
    bool useCutMetrics = false;
    bool useEmbeddedBC = false;
    bool useEmbeddedBCMap = false;
    bool useEmbeddedBCFallback = false;
    const Solvers::Poisson::CutCellMetrics* cutMetricsPtr = nullptr;
    unsigned long long cutMetricsRevision = 0;
    const PicEmbeddedBCMap* embeddedBCByIdPtr = nullptr;
    PicEmbeddedBCMap embeddedBCByIdValue;
    bool hasEmbeddedBCByIdValue = false;
    const PicEmbeddedBC* embeddedBCPtr = nullptr;

    PicEmbeddedBC embeddedBCValue;
    bool hasEmbeddedBCValue = false;
    PoissonFVSolverCellCentered::Options::EmbeddedBCMissingPolicy embeddedBCMissingPolicy =
        PoissonFVSolverCellCentered::Options::EmbeddedBCMissingPolicy::Error;

    // Unknown maps
    std::vector<int> rowOfCellIndex;   // size = mesh.cellIndexNum(), -1 for inactive
    std::vector<int> cellIndexOfRow;   // size = N

    // Geometry weights
    std::vector<double> volOfRow;      // size = N, plasma volume weight (alpha*area)

    // Cached RHS constant part from BC/EB contributions (independent of rhs field).
    Eigen::VectorXd bConst;            // size = N

    // Pure-Neumann info
    bool pureNeumann = false;
    int pinRow = 0;
    std::vector<int> pinnedRows;
    std::vector<std::vector<int>> floatingComponents;

    // Matrix and factorization
    Eigen::SparseMatrix<double> A;
    long long nnz = 0;

#ifdef EIGEN_SUPERLU_SUPPORT
    Eigen::SuperLU<Eigen::SparseMatrix<double>> superlu;
    bool superluReady = false;
#endif
    Eigen::SparseLU<Eigen::SparseMatrix<double>> sparselu;
    bool sparseluReady = false;

    bool usedSuperLU = false;
};

bool PoissonFVSolverCellCentered::cacheMatches(const FactorCache& c,
                                               const QuadtreeMesh& mesh,
                                               const FaceSegmentGraph& graph,
                                               const PicDomainBC& bc,
                                               const Options& opt,
                                               const Solvers::Poisson::CutCellMetrics* cutMetrics,
                                               const PicEmbeddedBC* embeddedBC,
                                               const PicEmbeddedBCMap* embeddedBCById)
{
    if (c.meshCellIndexNum != mesh.cellIndexNum())
        return false;
    if (c.graphSegmentCount != graph.segments().size())
        return false;
    if (!sameDomainBCStructure(c.bc, bc))
        return false;
    if (c.requestSuperLU != opt.useSuperLU)
        return false;

    const bool useCutMetricsNow =
        (cutMetrics != nullptr) &&
        (cutMetrics->betaBySegment().size() == graph.segments().size()) &&
        (cutMetrics->cellMetricsByIndex().size() == static_cast<size_t>(mesh.cellIndexNum()));

    const bool useEmbeddedBCMapNow =
        useCutMetricsNow &&
        (embeddedBCById != nullptr) &&
        (!embeddedBCById->empty());
    const bool useEmbeddedBCFallbackNow =
        useCutMetricsNow &&
        (embeddedBC != nullptr) &&
        embeddedBC->enabled;
    const bool useEmbeddedBCNow = useEmbeddedBCMapNow || useEmbeddedBCFallbackNow;

    if (c.useCutMetrics != useCutMetricsNow)
        return false;
    if (useCutMetricsNow && c.cutMetricsPtr != cutMetrics)
        return false;
    if (useCutMetricsNow && c.cutMetricsRevision != cutMetrics->revision())
        return false;

    if (c.useEmbeddedBC != useEmbeddedBCNow)
        return false;
    if (c.useEmbeddedBCMap != useEmbeddedBCMapNow)
        return false;
    if (useEmbeddedBCMapNow)
    {
        if (!c.hasEmbeddedBCByIdValue)
            return false;
        if (!sameEmbeddedBCMapStructure(c.embeddedBCByIdValue, *embeddedBCById))
            return false;
    }
    if (c.useEmbeddedBCFallback != useEmbeddedBCFallbackNow)
        return false;
    if (useEmbeddedBCFallbackNow)
    {
        if (!c.hasEmbeddedBCValue)
            return false;
        if (!sameEmbeddedBCStructure(c.embeddedBCValue, *embeddedBC))
            return false;
    }

    if (c.embeddedBCMissingPolicy != opt.embeddedBCMissingPolicy)
        return false;

    return true;
}

PoissonFVSolverCellCentered::PoissonFVSolverCellCentered() = default;
PoissonFVSolverCellCentered::~PoissonFVSolverCellCentered() = default;

PoissonFVSolverCellCentered::SolveReport
PoissonFVSolverCellCentered::setup(const QuadtreeMesh& mesh,
                                  const FaceSegmentGraph& graph,
                                  const PicDomainBC& bc,
                                  const Options& opt,
                                  const Solvers::Poisson::CutCellMetrics* cutMetrics,
                                  const PicEmbeddedBC* embeddedBC,
                                  const PicEmbeddedBCMap* embeddedBCById) const
{
    using Clock = std::chrono::steady_clock;
    const auto tBegin = Clock::now();
    auto secondsSince = [&](const Clock::time_point& t0) -> double {
        return std::chrono::duration_cast<std::chrono::duration<double>>(Clock::now() - t0).count();
    };

    SolveReport rep;

    clearCache();
    cache_ = std::make_unique<FactorCache>();
    cache_->meshCellIndexNum = mesh.cellIndexNum();
    cache_->graphSegmentCount = graph.segments().size();
    cache_->bc = bc;
    cache_->requestSuperLU = opt.useSuperLU;
    cache_->cutMetricsPtr = cutMetrics;
    cache_->cutMetricsRevision = cutMetrics ? cutMetrics->revision() : 0;
    cache_->embeddedBCPtr = embeddedBC;
    cache_->embeddedBCByIdPtr = embeddedBCById;
    cache_->embeddedBCMissingPolicy = opt.embeddedBCMissingPolicy;

    if (embeddedBC)
    {
        cache_->embeddedBCValue = *embeddedBC;
        cache_->hasEmbeddedBCValue = true;
    }
    else
    {
        cache_->hasEmbeddedBCValue = false;
    }
    if (embeddedBCById)
    {
        cache_->embeddedBCByIdValue = *embeddedBCById;
        cache_->hasEmbeddedBCByIdValue = true;
    }
    else
    {
        cache_->embeddedBCByIdValue.clear();
        cache_->hasEmbeddedBCByIdValue = false;
    }

    const int nCellIndex = mesh.cellIndexNum();



    // ---------------------------------------------------------------------
    // Build active unknown set from FaceSegmentGraph offsets.
    // ---------------------------------------------------------------------
    cache_->rowOfCellIndex.assign(static_cast<size_t>(nCellIndex), -1);
    cache_->cellIndexOfRow.clear();
    cache_->cellIndexOfRow.reserve(static_cast<size_t>(nCellIndex));

    for (int ci = 0; ci < nCellIndex; ++ci)
    {
        const auto range = graph.segmentRangeForCellIndex(ci);
        if (range.first < range.second)
        {
            cache_->rowOfCellIndex[static_cast<size_t>(ci)] = static_cast<int>(cache_->cellIndexOfRow.size());
            cache_->cellIndexOfRow.push_back(ci);
        }
    }

    const int N = static_cast<int>(cache_->cellIndexOfRow.size());
    rep.N = N;

    if (N == 0)
    {
        rep.success = false;
        rep.message = "No active cells in FaceSegmentGraph (N=0).";
        rep.timeTotal_s = secondsSince(tBegin);
        rep.timeAssemble_s = rep.timeTotal_s;
        rep.timeFactor_s = 0.0;
        rep.timeSolve_s = 0.0;
        rep.nnz = 0;
        rep.usedSuperLU = false;
        clearCache();
        return rep;
    }

    // Optional cut-cell metrics sanity checks.
    const bool useCutMetrics =
        (cutMetrics != nullptr) &&
        (cutMetrics->betaBySegment().size() == graph.segments().size()) &&
        (cutMetrics->cellMetricsByIndex().size() == static_cast<size_t>(mesh.cellIndexNum()));

    const bool useEmbeddedBCMap =
        useCutMetrics &&
        (embeddedBCById != nullptr) &&
        (!embeddedBCById->empty());
    const bool useEmbeddedBCFallback =
        useCutMetrics &&
        (embeddedBC != nullptr) &&
        embeddedBC->enabled;
    const bool useEmbeddedBC = useEmbeddedBCMap || useEmbeddedBCFallback;

    cache_->useCutMetrics = useCutMetrics;
    cache_->useEmbeddedBC = useEmbeddedBC;
    cache_->useEmbeddedBCMap = useEmbeddedBCMap;
    cache_->useEmbeddedBCFallback = useEmbeddedBCFallback;

    // Cache per-row plasma volume weight (alpha*area).
    cache_->volOfRow.assign(static_cast<size_t>(N), 0.0);

    double alphaMin = std::numeric_limits<double>::infinity();
    double alphaMax = 0.0;
    size_t nAlphaClamped = 0u;
    size_t nAlphaNonFinite = 0u;

    for (int i = 0; i < N; ++i)
    {
        const int cellIndex = cache_->cellIndexOfRow[i];
        const double area = safeCellArea(mesh, cellIndex);

        double alpha = 1.0;
        if (useCutMetrics)
        {
            alpha = cutMetrics->cellMetricsByIndex()[static_cast<size_t>(cellIndex)].alpha;
            if (!std::isfinite(alpha))
            {
                alpha = 1.0;
                ++nAlphaNonFinite;
            }
            if (alpha < 0.0)
            {
                alpha = 0.0;
                ++nAlphaClamped;
            }
            else if (alpha > 1.0)
            {
                alpha = 1.0;
                ++nAlphaClamped;
            }
        }

        cache_->volOfRow[static_cast<size_t>(i)] = area * alpha;
        alphaMin = std::min(alphaMin, alpha);
        alphaMax = std::max(alphaMax, alpha);
    }

    if (!useCutMetrics)
    {
        alphaMin = 1.0;
        alphaMax = 1.0;
    }

    // RHS constant part from BC/EB contributions.
    cache_->bConst.resize(N);
    cache_->bConst.setZero();

    // ---------------------------------------------------------------------
    // Identify cells with degenerate volume fraction for extrapolation.
    // Instead of assembling a near-singular flux equation, these cells
    // will receive a Laplacian averaging stencil: phi_i = mean(phi_nbrs).
    // ---------------------------------------------------------------------
    std::vector<bool> extrapolateRow(static_cast<size_t>(N), false);
    size_t nExtrapolate = 0u;
    if (useCutMetrics && opt.alphaMinThreshold > 0.0)
    {
        for (int i = 0; i < N; ++i)
        {
            const int cellIndex = cache_->cellIndexOfRow[static_cast<size_t>(i)];
            const double alpha = cutMetrics->cellMetricsByIndex()[static_cast<size_t>(cellIndex)].alpha;
            if (std::isfinite(alpha) && alpha < opt.alphaMinThreshold)
            {
                extrapolateRow[static_cast<size_t>(i)] = true;
                cache_->volOfRow[static_cast<size_t>(i)] = 0.0;
                ++nExtrapolate;
            }
        }
    }

    // ---------------------------------------------------------------------
    // Assemble sparse matrix (triplets) and diagonal diagnostics.
    // This largely mirrors the legacy solve() assembly but excludes rhs term.
    // ---------------------------------------------------------------------
    std::vector<Eigen::Triplet<double>> trips;
    trips.reserve(graph.segments().size() * 2u + 8u);

    std::vector<double> diag(static_cast<size_t>(N), 0.0);
    std::vector<std::vector<int>> neighbors(static_cast<size_t>(N));
    std::vector<unsigned char> anchoredRow(static_cast<size_t>(N), 0u);

    size_t nSegTotal = graph.segments().size();
    size_t nInteriorUsed = 0u;
    size_t nInteriorSkippedMissing = 0u;
    size_t nInteriorSkippedDegenerate = 0u;

    size_t nBoundaryDirichlet = 0u;
    size_t nBoundaryNeumann = 0u;
    size_t nBoundaryPeriodic = 0u;
    size_t nBoundarySkippedMissing = 0u;
    size_t nBoundarySkippedDegenerate = 0u;

    bool warnedRobin = false;
    bool warnedPeriodicMismatch = false;

    const auto& segs = graph.segments();

    // Track beta statistics if cut-metrics are used.
    double betaMin = std::numeric_limits<double>::infinity();
    double betaMax = 0.0;
    size_t nBetaClamped = 0u;
    size_t nBetaNonFinite = 0u;
    size_t nZeroOpenArea = 0u;

    for (size_t sid = 0; sid < segs.size(); ++sid)
    {
        const auto& seg = segs[sid];

        const int Kindex = mesh.cells(seg.K).index();
        if (Kindex < 0)
            continue;

        const int i = cache_->rowOfCellIndex[static_cast<size_t>(Kindex)];
        if (i < 0)
            continue;

        // Face open fraction beta_f (default=1).
        double beta = 1.0;
        if (useCutMetrics)
        {
            beta = cutMetrics->betaBySegment()[sid];
            if (!std::isfinite(beta))
            {
                beta = 1.0;
                ++nBetaNonFinite;
            }
            if (beta < 0.0)
            {
                beta = 0.0;
                ++nBetaClamped;
            }
            else if (beta > 1.0)
            {
                beta = 1.0;
                ++nBetaClamped;
            }
            betaMin = std::min(betaMin, beta);
            betaMax = std::max(betaMax, beta);
        }

        const double AfOpen = beta * seg.Af;

        if (!(seg.df > 0.0) || !(AfOpen >= 0.0))
        {
            if (seg.boundary == Boundary::non)
                ++nInteriorSkippedDegenerate;
            else
                ++nBoundarySkippedDegenerate;
            continue;
        }

        if (AfOpen == 0.0)
        {
            ++nZeroOpenArea;
            continue;
        }

        if (seg.boundary == Boundary::non)
        {
            if (seg.L == -1)
            {
                ++nInteriorSkippedMissing;
                continue;
            }

            const int Lindex = mesh.cells(seg.L).index();
            if (Lindex < 0)
            {
                ++nInteriorSkippedMissing;
                continue;
            }

            const int j = cache_->rowOfCellIndex[static_cast<size_t>(Lindex)];
            if (j < 0)
            {
                ++nInteriorSkippedMissing;
                continue;
            }

            const double c = AfOpen / seg.df;
            trips.emplace_back(i, i, c);
            trips.emplace_back(i, j, -c);
            diag[static_cast<size_t>(i)] += c;
            neighbors[static_cast<size_t>(i)].push_back(j);
            neighbors[static_cast<size_t>(j)].push_back(i);
            ++nInteriorUsed;
        }
        else
        {
            // Domain boundary segment
            const PicBoundarySide& s = bcSide(bc, seg.boundary);
            PicBCType t = s.type;

            if (t == PicBCType::Robin)
            {
                // Only Dirichlet/Neumann are supported on the domain boundary.
                // Treat Robin as Dirichlet(value) for now.
                if (!warnedRobin)
                {
                    warnedRobin = true;
                    if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
                    {
                        std::cout << "[PoissonFVSolverCellCentered] Warning: Robin BC encountered on domain boundary; "
                                  << "treating as Dirichlet(value).\n";
                    }
                }
                t = PicBCType::Dirichlet;
            }

            if (t == PicBCType::Dirichlet)
            {
                if (!extrapolateRow[static_cast<size_t>(i)])
                {
                    const double c = AfOpen / seg.df;
                    trips.emplace_back(i, i, c);
                    diag[static_cast<size_t>(i)] += c;
                    cache_->bConst[i] += c * s.value;
                    anchoredRow[static_cast<size_t>(i)] = 1u;
                }
                ++nBoundaryDirichlet;
            }
            else if (t == PicBCType::Neumann)
            {
                if (!extrapolateRow[static_cast<size_t>(i)])
                    cache_->bConst[i] += s.value * AfOpen;
                ++nBoundaryNeumann;
            }
            else if (t == PicBCType::Periodic)
            {
                const Boundary mateB = PeriodicFaceUtils::mateBoundary(seg.boundary);
                if (!warnedPeriodicMismatch)
                {
                    if (mateB == Boundary::non || bcSide(bc, mateB).type != PicBCType::Periodic)
                    {
                        warnedPeriodicMismatch = true;
                        if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
                        {
                            std::cout << "[PoissonFVSolverCellCentered] Warning: Periodic BC on '"
                                      << boundaryToString(seg.boundary)
                                      << "' but mate side is not Periodic. Assembly will still wrap.\n";
                        }
                    }
                }

                const Vector2D shift = PeriodicFaceUtils::shiftToMate(mesh, seg.boundary);
                int Llabel = -1;
                if (!PeriodicFaceUtils::findMateCellLabel(mesh, seg, shift, Llabel))
                {
                    ++nBoundarySkippedMissing;
                    continue;
                }

                const int Lindex = mesh.cells(Llabel).index();
                if (Lindex < 0)
                {
                    ++nBoundarySkippedMissing;
                    continue;
                }

                const int j = cache_->rowOfCellIndex[static_cast<size_t>(Lindex)];
                if (j < 0)
                {
                    ++nBoundarySkippedMissing;
                    continue;
                }

                const double df = PeriodicFaceUtils::centerDistance(mesh, seg, Llabel, shift);
                if (!(df > 0.0))
                {
                    ++nBoundarySkippedDegenerate;
                    continue;
                }

                const double c = AfOpen / df;
                if (!extrapolateRow[static_cast<size_t>(i)])
                {
                    trips.emplace_back(i, i, c);
                    trips.emplace_back(i, j, -c);
                    diag[static_cast<size_t>(i)] += c;
                }
                neighbors[static_cast<size_t>(i)].push_back(j);
                neighbors[static_cast<size_t>(j)].push_back(i);
                ++nBoundaryPeriodic;
            }
        }
    }

    if (!useCutMetrics)
    {
        betaMin = 1.0;
        betaMax = 1.0;
    }

    // ---------------------------------------------------------------------
    // Embedded-boundary closure on cut cells (adds diagonal and/or RHS const).
    // ---------------------------------------------------------------------
    size_t nEBCellsTouched = 0u;
    size_t nEBCellsConstrained = 0u;
    size_t nEBRobin = 0u;
    size_t nEBDirichlet = 0u;
    size_t nEBNeumann = 0u;
    size_t nEBSkipped = 0u;
    size_t nEBDegenerate = 0u;

    double ebAddedDiagMin = std::numeric_limits<double>::infinity();
    double ebAddedDiagMax = 0.0;

    if (useEmbeddedBC)
    {
        size_t nEBMapUsed = 0u;
        size_t nEBFallbackUsed = 0u;
        size_t nEBMissing = 0u;
        bool ebMissingFatal = false;
        std::string ebMissingMsg;

        for (int i = 0; i < N; ++i)
        {
            const int cellIndex = cache_->cellIndexOfRow[i];
            const auto& cm = cutMetrics->cellMetricsByIndex()[static_cast<size_t>(cellIndex)];

            if (!cm.isCut || !(cm.A_gamma > 0.0))
                continue;

            if (extrapolateRow[static_cast<size_t>(i)])
                continue;

            if (!(cache_->volOfRow[static_cast<size_t>(i)] > 0.0))
                continue;

            ++nEBCellsTouched;

            const PicEmbeddedBC* eb = nullptr;
            if (useEmbeddedBCMap && cm.gammaBoundaryId > 0)
            {
                auto it = embeddedBCById->find(cm.gammaBoundaryId);
                if (it != embeddedBCById->end() && it->second.enabled)
                {
                    eb = &it->second;
                    ++nEBMapUsed;
                }
            }
            if (!eb && useEmbeddedBCFallback)
            {
                eb = embeddedBC;
                ++nEBFallbackUsed;
            }
            if (!eb)
            {
                ++nEBMissing;
                if (opt.embeddedBCMissingPolicy == Options::EmbeddedBCMissingPolicy::Error)
                {
                    ebMissingFatal = true;
                    ebMissingMsg = "Missing embedded BC for boundaryId " + std::to_string(cm.gammaBoundaryId) +
                                   " (no map entry and no fallback).";
                    break;
                }
                continue;
            }

            const double a = eb->a;
            const double b = eb->b;
            const double g = eb->g;
            const PicBCType t = eb->type;
            const double A_gamma = cm.A_gamma;

            if (t == PicBCType::Robin)
            {
                if (std::abs(b) > 0.0)
                {
                    const double addDiag = (a / b) * A_gamma;
                    const double addRhs = (g / b) * A_gamma;

                    if (!std::isfinite(addDiag) || !std::isfinite(addRhs))
                    {
                        ++nEBDegenerate;
                        continue;
                    }

                    if (addDiag != 0.0)
                    {
                        trips.emplace_back(i, i, addDiag);
                        diag[static_cast<size_t>(i)] += addDiag;
                        ebAddedDiagMin = std::min(ebAddedDiagMin, addDiag);
                        ebAddedDiagMax = std::max(ebAddedDiagMax, addDiag);
                        if (a != 0.0)
                        {
                            ++nEBCellsConstrained;
                            anchoredRow[static_cast<size_t>(i)] = 1u;
                        }
                    }

                    cache_->bConst[i] += addRhs;
                    ++nEBRobin;
                }
                else
                {
                    // Dirichlet limit of Robin (b=0): a*phi=g.
                    if (std::abs(a) <= 0.0)
                    {
                        ++nEBDegenerate;
                        continue;
                    }
                    const double phiD = g / a;
                    const double d = stabilizedEmbeddedDistance(cm.d_gamma,
                                                                cache_->volOfRow[static_cast<size_t>(i)]);

                    const double c = A_gamma / d;
                    if (!std::isfinite(c) || c <= 0.0)
                    {
                        ++nEBDegenerate;
                        continue;
                    }

                    trips.emplace_back(i, i, c);
                    diag[static_cast<size_t>(i)] += c;
                    cache_->bConst[i] += c * phiD;

                    ebAddedDiagMin = std::min(ebAddedDiagMin, c);
                    ebAddedDiagMax = std::max(ebAddedDiagMax, c);
                    ++nEBCellsConstrained;
                    anchoredRow[static_cast<size_t>(i)] = 1u;
                    ++nEBDirichlet;
                }
            }
            else if (t == PicBCType::Dirichlet)
            {
                const double phiD = g;
                const double d = stabilizedEmbeddedDistance(cm.d_gamma,
                                                            cache_->volOfRow[static_cast<size_t>(i)]);

                const double c = A_gamma / d;
                if (!std::isfinite(c) || c <= 0.0)
                {
                    ++nEBDegenerate;
                    continue;
                }

                trips.emplace_back(i, i, c);
                diag[static_cast<size_t>(i)] += c;
                cache_->bConst[i] += c * phiD;

                ebAddedDiagMin = std::min(ebAddedDiagMin, c);
                ebAddedDiagMax = std::max(ebAddedDiagMax, c);
                ++nEBCellsConstrained;
                anchoredRow[static_cast<size_t>(i)] = 1u;
                ++nEBDirichlet;
            }
            else if (t == PicBCType::Neumann)
            {
                const double addRhs = g * A_gamma;
                if (!std::isfinite(addRhs))
                {
                    ++nEBDegenerate;
                    continue;
                }
                cache_->bConst[i] += addRhs;
                ++nEBNeumann;
            }
            else
            {
                ++nEBSkipped;
            }
        }

        if (ebMissingFatal)
        {
            rep.success = false;
            rep.message = ebMissingMsg;
            rep.timeTotal_s = secondsSince(tBegin);
            rep.timeAssemble_s = rep.timeTotal_s;
            rep.timeFactor_s = 0.0;
            rep.timeSolve_s = 0.0;
            rep.nnz = 0;
            rep.usedSuperLU = false;
            clearCache();
            return rep;
        }
        (void)nEBMapUsed;
        (void)nEBFallbackUsed;
        (void)nEBMissing;
    }

    // ---------------------------------------------------------------------
    // Extrapolation stencil for small-alpha cells.
    // Replace degenerate flux equations with phi_i = mean(phi_neighbors).
    // ---------------------------------------------------------------------
    if (nExtrapolate > 0u)
    {
        for (int i = 0; i < N; ++i)
        {
            if (!extrapolateRow[static_cast<size_t>(i)])
                continue;

            // Deduplicate neighbor list and remove self-references.
            auto& nbrs = neighbors[static_cast<size_t>(i)];
            std::sort(nbrs.begin(), nbrs.end());
            nbrs.erase(std::unique(nbrs.begin(), nbrs.end()), nbrs.end());
            nbrs.erase(std::remove(nbrs.begin(), nbrs.end(), i), nbrs.end());

            const int count = static_cast<int>(nbrs.size());
            if (count > 0)
            {
                trips.emplace_back(i, i, static_cast<double>(count));
                for (int j : nbrs)
                    trips.emplace_back(i, j, -1.0);
                diag[static_cast<size_t>(i)] = static_cast<double>(count);
            }
            else
            {
                trips.emplace_back(i, i, 1.0);
                diag[static_cast<size_t>(i)] = 1.0;
            }
            cache_->bConst[i] = 0.0;
            anchoredRow[static_cast<size_t>(i)] = 1u;
        }

        if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
        {
            std::cout << "[PoissonFVSolverCellCentered] setup(): "
                      << nExtrapolate
                      << " cell(s) with alpha < " << opt.alphaMinThreshold
                      << " replaced with Laplacian averaging stencil.\n";
        }
    }

    // ---------------------------------------------------------------------
    // Detect pure Neumann and apply gauge fix (matrix-only).
    // Compatibility correction depends on rhs and will be applied in the RHS pass.
    // ---------------------------------------------------------------------
    const bool domainPureNeumann =
        (!isDirichletLike(bc.left.type)) &&
        (!isDirichletLike(bc.right.type)) &&
        (!isDirichletLike(bc.bottom.type)) &&
        (!isDirichletLike(bc.top.type));

    const bool embeddedAnchors = (useEmbeddedBC && (nEBCellsConstrained > 0u));
    const bool pureNeumann = domainPureNeumann && (!embeddedAnchors);

    cache_->pureNeumann = pureNeumann;
    cache_->pinRow = 0;
    cache_->pinnedRows.clear();
    cache_->floatingComponents.clear();
    size_t nFloatingComponents = 0u;

    if (pureNeumann)
    {
        int pinRow = 0;
        for (int i = 0; i < N; ++i)
        {
            if (cache_->volOfRow[static_cast<size_t>(i)] > 0.0)
            {
                pinRow = i;
                break;
            }
        }
        cache_->pinRow = pinRow;

        trips.emplace_back(pinRow, pinRow, 1.0);
        diag[static_cast<size_t>(pinRow)] += 1.0;

        if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
            std::cout << "[PoissonFVSolverCellCentered] Pure Neumann detected in setup(): applied gauge fix (pin row "
                      << pinRow << ")\n";
    }
    else
    {
        nFloatingComponents = stabilizeFloatingComponents(neighbors,
                                                          anchoredRow,
                                                          cache_->volOfRow,
                                                          cache_->bConst,
                                                          trips,
                                                          diag,
                                                          &cache_->pinnedRows,
                                                          &cache_->floatingComponents);

        if ((nFloatingComponents > 0u) && ConsoleControl::enabled(ConsoleControl::Level::kDebug))
        {
            std::cout << "[PoissonFVSolverCellCentered] setup(): stabilized "
                      << nFloatingComponents
                      << " floating component(s) with compatibility correction + gauge pin.\n";
        }
    }

    // ---------------------------------------------------------------------
    // Safety: pin rows with zero/invalid diagonal to avoid singular factorization.
    // This can occur if an active cell ends up fully blocked (no open faces).
    // ---------------------------------------------------------------------
    size_t nPinnedZeroDiag = 0u;
    for (int i = 0; i < N; ++i)
    {
        const double d = diag[static_cast<size_t>(i)];
        if (!std::isfinite(d) || !(d > 0.0))
        {
            trips.emplace_back(i, i, 1.0);
            diag[static_cast<size_t>(i)] += 1.0;
            cache_->pinnedRows.push_back(i);
            ++nPinnedZeroDiag;
        }
    }

    // Finalize sparse matrix.
    cache_->A.resize(N, N);
    cache_->A.setFromTriplets(trips.begin(), trips.end());
    cache_->A.makeCompressed();
    cache_->nnz = static_cast<long long>(cache_->A.nonZeros());

    rep.timeAssemble_s = secondsSince(tBegin);
    rep.nnz = cache_->nnz;

    // ---------------------------------------------------------------------
    // Factorize using SuperLU if available/requested; otherwise SparseLU.
    // ---------------------------------------------------------------------
    cache_->usedSuperLU = false;

#ifdef EIGEN_SUPERLU_SUPPORT
    if (opt.useSuperLU)
    {
        cache_->superlu.options().ColPerm = COLAMD;
        cache_->superlu.options().DiagPivotThresh = 0.1;

        const auto tFac0 = Clock::now();
        cache_->superlu.analyzePattern(cache_->A);
        cache_->superlu.factorize(cache_->A);
        rep.timeFactor_s += std::chrono::duration_cast<std::chrono::duration<double>>(Clock::now() - tFac0).count();

        if (cache_->superlu.info() == Eigen::Success)
        {
            cache_->usedSuperLU = true;
            cache_->superluReady = true;
        }
        else
        {
            if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
                std::cout << "[PoissonFVSolverCellCentered] SuperLU factorization failed in setup(); falling back to SparseLU.\n";
            cache_->superluReady = false;
        }
    }
#endif

    if (!cache_->usedSuperLU)
    {
        const auto tFac0 = Clock::now();
        cache_->sparselu.analyzePattern(cache_->A);
        cache_->sparselu.factorize(cache_->A);
        rep.timeFactor_s += std::chrono::duration_cast<std::chrono::duration<double>>(Clock::now() - tFac0).count();

        if (cache_->sparselu.info() != Eigen::Success)
        {
            rep.success = false;
            rep.message = "SparseLU factorization failed in setup().";
            rep.timeTotal_s = secondsSince(tBegin);
            rep.usedSuperLU = false;
            clearCache();
            return rep;
        }

        cache_->sparseluReady = true;
        cache_->usedSuperLU = false;
    }

    rep.success = true;
    rep.message = cache_->usedSuperLU ? "setup() completed: factorized with SuperLU." : "setup() completed: factorized with SparseLU.";
    rep.usedSuperLU = cache_->usedSuperLU;
    rep.timeSolve_s = 0.0;
    rep.timeTotal_s = secondsSince(tBegin);

    // Info-level summary (debug only to avoid spam).
    if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
    {
        std::cout << "[PoissonFVSolverCellCentered] setup(): N=" << N << ", nnz=" << cache_->A.nonZeros()
                  << ", solver=" << (cache_->usedSuperLU ? "SuperLU" : "SparseLU")
                  << ", cutMetrics=" << (useCutMetrics ? "on" : "off")
                  << ", embeddedBC=" << (useEmbeddedBC ? "on" : "off")
                  << "\n";

        if (useCutMetrics)
        {
            std::cout << "  alphaMin=" << alphaMin << ", alphaMax=" << alphaMax
                      << ", betaMin=" << (std::isfinite(betaMin) ? betaMin : 1.0)
                      << ", betaMax=" << betaMax
                      << ", betaClamped=" << nBetaClamped
                      << ", betaNonFinite=" << nBetaNonFinite
                      << ", alphaClamped=" << nAlphaClamped
                      << ", alphaNonFinite=" << nAlphaNonFinite
                      << ", zeroOpenSeg=" << nZeroOpenArea
                      << "\n";
        }

        if (useEmbeddedBC)
        {
            std::cout << "  EB: touched=" << nEBCellsTouched
                      << ", constrained=" << nEBCellsConstrained
                      << ", robin=" << nEBRobin
                      << ", dir=" << nEBDirichlet
                      << ", neu=" << nEBNeumann
                      << ", skipped=" << nEBSkipped
                      << ", degenerate=" << nEBDegenerate;

            if (std::isfinite(ebAddedDiagMin))
                std::cout << ", ebDiagMin=" << ebAddedDiagMin << ", ebDiagMax=" << ebAddedDiagMax;
            else
                std::cout << ", ebDiagMin=0, ebDiagMax=0";

            std::cout << ", embeddedAnchors=" << (embeddedAnchors ? "true" : "false")
                      << "\n";
        }

        if (nPinnedZeroDiag > 0u)
        {
            std::cout << "  pinned zero-diag rows=" << nPinnedZeroDiag << "\n";
        }
        if (nFloatingComponents > 0u)
        {
            std::cout << "  floating components stabilized=" << nFloatingComponents << "\n";
        }

        if (ConsoleControl::enabled(ConsoleControl::Level::kTrace))
        {
            std::cout << "  segments total                = " << nSegTotal << "\n";
            std::cout << "  interior used                 = " << nInteriorUsed << "\n";
            std::cout << "  interior skipped (L missing/inactive) = " << nInteriorSkippedMissing << "\n";
            std::cout << "  interior skipped (degenerate df/Af)   = " << nInteriorSkippedDegenerate << "\n";
            std::cout << "  boundary Dirichlet segments    = " << nBoundaryDirichlet << "\n";
            std::cout << "  boundary Neumann segments      = " << nBoundaryNeumann << "\n";
            std::cout << "  boundary Periodic segments     = " << nBoundaryPeriodic << "\n";
            std::cout << "  boundary skipped (mate missing/inactive) = " << nBoundarySkippedMissing << "\n";
            std::cout << "  boundary skipped (degenerate df/Af)   = " << nBoundarySkippedDegenerate << "\n";
            std::cout << "  BC types: left=" << bcTypeName(bc.left.type)
                      << ", right=" << bcTypeName(bc.right.type)
                      << ", bottom=" << bcTypeName(bc.bottom.type)
                      << ", top=" << bcTypeName(bc.top.type) << "\n";
        }
    }

    return rep;
}

bool PoissonFVSolverCellCentered::solveWithCache(const QuadtreeMesh& mesh,
                                                const FaceSegmentGraph& graph,
                                                const PicDomainBC& bc,
                                                const CellField<double>& rhs,
                                                CellField<double>& phi,
                                                SolveReport* report,
                                                const Options& opt,
                                                const Solvers::Poisson::CutCellMetrics* cutMetrics,
                                                const PicEmbeddedBC* embeddedBC,
                                                const PicEmbeddedBCMap* embeddedBCById) const
{
    using Clock = std::chrono::steady_clock;
    const auto tBegin = Clock::now();
    auto secondsSince = [&](const Clock::time_point& t0) -> double {
        return std::chrono::duration_cast<std::chrono::duration<double>>(Clock::now() - t0).count();
    };

    // Ensure phi has correct storage.
    if (phi.size() != mesh.cellIndexNum() || !phi.value)
    {
        phi.initialize(mesh, 0.0);
    }

    // ---------------------------------------------------------------------
    // Cache control + logging
    // - cacheEnabled=false: solve still works, but factorization is not reused
    // - cacheRebuildEvery=N: force periodic rebuild every N solve calls
    // - cacheHit/miss logs:
    //     verbose=true: per-call log
    //     else cacheLogEvery>0: periodic summary
    //     else cacheLogOnMiss: log misses only
    // ---------------------------------------------------------------------
    const long long callIndex = cacheSolveCalls_++;
    const bool forceRebuild = (opt.cacheRebuildEvery > 0) && (callIndex > 0) &&
                              (callIndex % static_cast<long long>(opt.cacheRebuildEvery) == 0);

    auto logLine = [&](const std::string& s) {
        std::cout << "[PoissonFVSolverCellCentered] " << s << "\n";
    };
    const bool doSummaryLog = opt.cacheEnabled &&
                              (!opt.verbose) && (opt.cacheLogEvery > 0) &&
                              ((callIndex + 1) % static_cast<long long>(opt.cacheLogEvery) == 0);

    bool rebuiltThisCall = false;
    SolveReport repSetup;

    std::string missReason;
    const bool haveCache = static_cast<bool>(cache_);
    const bool sigOk = haveCache && cacheMatches(*cache_, mesh, graph, bc, opt, cutMetrics, embeddedBC, embeddedBCById);
    const bool needRebuild = (!opt.cacheEnabled) || (!haveCache) || (!sigOk) || forceRebuild;

    if (needRebuild)
    {
        ++cacheMisses_;
        if (!opt.cacheEnabled)
            missReason = "disabled";
        else if (!haveCache)
            missReason = "empty";
        else if (forceRebuild)
            missReason = "rebuildEvery";
        else
            missReason = "signatureMismatch";

        if (opt.cacheEnabled && (opt.verbose || opt.cacheLogOnMiss))
            logLine("cache MISS (" + missReason + ") -> rebuilding factorization");

        // Drop any stale cache then rebuild.
        clearCache();
        repSetup = setup(mesh, graph, bc, opt, cutMetrics, embeddedBC, embeddedBCById);
        if (!repSetup.success)
        {
            if (report)
                *report = repSetup;
            return false;
        }
        rebuiltThisCall = true;
    }
    else
    {
        ++cacheHits_;
        if (opt.verbose)
            logLine("cache HIT -> reuse factorization");
    }

    if (doSummaryLog)
    {
        logLine("cache summary: calls=" + std::to_string(cacheSolveCalls_) +
                " hits=" + std::to_string(cacheHits_) +
                " misses=" + std::to_string(cacheMisses_));
    }

    // Cache must now exist and match.
    if (!cache_)
    {
        if (report)
        {
            report->success = false;
            report->message = "solveWithCache(): internal error, cache is null after setup.";
            report->timeTotal_s = secondsSince(tBegin);
        }
        return false;
    }
    if (!cacheMatches(*cache_, mesh, graph, bc, opt, cutMetrics, embeddedBC, embeddedBCById))
    {
        if (report)
        {
            report->success = false;
            report->message = "solveWithCache(): internal error, cache signature mismatch after setup.";
            report->timeTotal_s = secondsSince(tBegin);
        }
        return false;
    }

    const bool useCutMetricsNow =
        (cutMetrics != nullptr) &&
        (cutMetrics->betaBySegment().size() == graph.segments().size()) &&
        (cutMetrics->cellMetricsByIndex().size() == static_cast<size_t>(mesh.cellIndexNum()));

    const bool useEmbeddedBCMapNow =
        useCutMetricsNow &&
        (embeddedBCById != nullptr) &&
        (!embeddedBCById->empty());
    const bool useEmbeddedBCFallbackNow =
        useCutMetricsNow &&
        (embeddedBC != nullptr) &&
        embeddedBC->enabled;
    const bool useEmbeddedBCNow = useEmbeddedBCMapNow || useEmbeddedBCFallbackNow;

    bool bcChanged = !sameDomainBC(cache_->bc, bc);
    bool embeddedChanged = false;
    if (useEmbeddedBCMapNow)
    {
        if (!cache_->hasEmbeddedBCByIdValue)
            embeddedChanged = true;
        else
            embeddedChanged = !sameEmbeddedBCMap(cache_->embeddedBCByIdValue, *embeddedBCById);
    }
    else if (cache_->hasEmbeddedBCByIdValue)
    {
        embeddedChanged = true;
    }

    if (useEmbeddedBCFallbackNow)
    {
        if (!cache_->hasEmbeddedBCValue)
            embeddedChanged = true;
        else
            embeddedChanged = embeddedChanged || !sameEmbeddedBC(cache_->embeddedBCValue, *embeddedBC);
    }
    else if (cache_->hasEmbeddedBCValue)
    {
        embeddedChanged = true;
    }

    if (bcChanged || embeddedChanged)
    {
        std::string errMsg;
        if (!rebuildBCConst(mesh, graph, bc, cutMetrics, embeddedBC, embeddedBCById, opt,
                            cache_->rowOfCellIndex, cache_->cellIndexOfRow, cache_->volOfRow, cache_->bConst, &errMsg))
        {
            if (report)
            {
                report->success = false;
                report->message = errMsg.empty() ? "solveWithCache(): embedded BC rebuild failed." : errMsg;
                report->timeTotal_s = secondsSince(tBegin);
            }
            return false;
        }

        cache_->bc = bc;
        if (useEmbeddedBCMapNow && embeddedBCById)
        {
            cache_->embeddedBCByIdValue = *embeddedBCById;
            cache_->hasEmbeddedBCByIdValue = true;
        }
        else
        {
            cache_->embeddedBCByIdValue.clear();
            cache_->hasEmbeddedBCByIdValue = false;
        }
        if (useEmbeddedBCFallbackNow && embeddedBC)
        {
            cache_->embeddedBCValue = *embeddedBC;
            cache_->hasEmbeddedBCValue = true;
        }
        else
        {
            cache_->hasEmbeddedBCValue = false;
        }
    }

    // RHS assembly: b = (rhs * plasmaVolume) + bConst.
    const int N = static_cast<int>(cache_->cellIndexOfRow.size());
    if (N <= 0)
    {
        if (report)
        {
            report->success = false;
            report->iterations = 0;
            report->finalResidual = 0.0;
            report->message = "solveWithCache(): cache has N=0.";
            report->timeTotal_s = secondsSince(tBegin);
            report->timeAssemble_s = report->timeTotal_s;
            report->timeFactor_s = 0.0;
            report->timeSolve_s = 0.0;
            report->N = 0;
            report->nnz = 0;
            report->usedSuperLU = false;
        }
        return false;
    }

    const auto tRhs0 = Clock::now();
    Eigen::VectorXd bVec = cache_->bConst;
    for (int i = 0; i < N; ++i)
    {
        const int cellIndex = cache_->cellIndexOfRow[i];
        double rhsVal = 0.0;
        if (cellIndex >= 0 && cellIndex < rhs.size() && rhs.value)
            rhsVal = rhs.value[cellIndex];
        if (!std::isfinite(rhsVal))
            rhsVal = 0.0;
        bVec[i] += rhsVal * cache_->volOfRow[static_cast<size_t>(i)];
    }

    // Pure Neumann: compatibility correction uses the same convention as legacy solve().
    if (cache_->pureNeumann)
    {
        double totalB = 0.0;
        double totalV = 0.0;
        for (int i = 0; i < N; ++i)
        {
            totalB += bVec[i];
            totalV += cache_->volOfRow[static_cast<size_t>(i)];
        }
        if (totalV > 0.0)
        {
            const double meanBperV = totalB / totalV;
            for (int i = 0; i < N; ++i)
            {
                bVec[i] -= meanBperV * cache_->volOfRow[static_cast<size_t>(i)];
            }
        }

        // Gauge pin row (matrix already pinned in setup). Ensure RHS respects it.
        const int pinRow = cache_->pinRow;
        if (pinRow >= 0 && pinRow < N)
            bVec[pinRow] = 0.0;
    }
    else if (!cache_->floatingComponents.empty())
    {
        for (const auto& componentRows : cache_->floatingComponents)
        {
            double totalB = 0.0;
            double totalV = 0.0;
            for (const int row : componentRows)
            {
                if (row < 0 || row >= N)
                    continue;
                totalB += bVec[row];
                totalV += cache_->volOfRow[static_cast<size_t>(row)];
            }
            if (totalV <= 0.0)
                continue;
            const double meanBperV = totalB / totalV;
            for (const int row : componentRows)
            {
                if (row < 0 || row >= N)
                    continue;
                bVec[row] -= meanBperV * cache_->volOfRow[static_cast<size_t>(row)];
            }
        }
    }

    if (!cache_->pinnedRows.empty())
    {
        for (const int row : cache_->pinnedRows)
        {
            if (row >= 0 && row < N)
                bVec[row] = 0.0;
        }
    }

    const double tAssemble_s = std::chrono::duration_cast<std::chrono::duration<double>>(Clock::now() - tRhs0).count();

    // Solve using cached factorization.
    const auto tSol0 = Clock::now();
    Eigen::VectorXd x;
    bool usedSuperLU = cache_->usedSuperLU;
    bool ok = false;

#ifdef EIGEN_SUPERLU_SUPPORT
    if (usedSuperLU && cache_->superluReady)
    {
        x = cache_->superlu.solve(bVec);
        ok = (cache_->superlu.info() == Eigen::Success);
    }
#endif

    if (!ok)
    {
        // Fallback to SparseLU (either because SuperLU not ready, not supported,
        // or SuperLU solve failed).
        if (!cache_->sparseluReady)
        {
            if (report)
            {
                report->success = false;
                report->iterations = 0;
                report->finalResidual = 0.0;
                report->message = "solveWithCache(): cache has no ready factorization.";
                report->timeAssemble_s = tAssemble_s;
                report->timeFactor_s = 0.0;
                report->timeSolve_s = secondsSince(tSol0);
                report->timeTotal_s = secondsSince(tBegin);
                report->N = N;
                report->nnz = cache_->nnz;
                report->usedSuperLU = false;
            }
            return false;
        }

        x = cache_->sparselu.solve(bVec);
        ok = (cache_->sparselu.info() == Eigen::Success);
        usedSuperLU = false;
    }

    const double tSolve_s = std::chrono::duration_cast<std::chrono::duration<double>>(Clock::now() - tSol0).count();

    if (!ok)
    {
        if (report)
        {
            report->success = false;
            report->iterations = 0;
            report->finalResidual = 0.0;
            report->message = usedSuperLU ? "solveWithCache(): SuperLU solve failed." : "solveWithCache(): SparseLU solve failed.";
            report->timeAssemble_s = tAssemble_s;
            report->timeFactor_s = 0.0;
            report->timeSolve_s = tSolve_s;
            report->timeTotal_s = secondsSince(tBegin);
            report->N = N;
            report->nnz = cache_->nnz;
            report->usedSuperLU = usedSuperLU;
        }
        return false;
    }

    // Scatter solution back to phi (cell-index aligned).
    for (int i = 0; i < N; ++i)
    {
        const int cellIndex = cache_->cellIndexOfRow[i];
        if (cellIndex >= 0 && cellIndex < mesh.cellIndexNum())
            phi.value[cellIndex] = x[i];
    }

    if (report)
    {
        report->success = true;
        report->iterations = 0;
        report->finalResidual = 0.0;
        const std::string solverName = usedSuperLU ? "SuperLU" : "SparseLU";
        if (needRebuild)
            report->message = "cache miss (" + missReason + "): rebuilt + solved using " + solverName + ".";
        else
            report->message = "cache hit: solved using cached " + solverName + ".";

        report->timeAssemble_s = tAssemble_s + (rebuiltThisCall ? repSetup.timeAssemble_s : 0.0);
        report->timeFactor_s = rebuiltThisCall ? repSetup.timeFactor_s : 0.0;
        report->timeSolve_s = tSolve_s;
        report->timeTotal_s = secondsSince(tBegin);
        report->N = N;
        report->nnz = cache_->nnz;
        report->usedSuperLU = usedSuperLU;
    }

    const bool clearAfter = !opt.cacheEnabled;
    if (clearAfter)
        clearCache();

    return true;
}

bool PoissonFVSolverCellCentered::hasCache() const
{
    return static_cast<bool>(cache_);
}

void PoissonFVSolverCellCentered::clearCache() const
{
    cache_.reset();
}

bool PoissonFVSolverCellCentered::solve(const QuadtreeMesh& mesh,
                                       const FaceSegmentGraph& graph,
                                       const PicDomainBC& bc,
                                       const CellField<double>& rhs,
                                       CellField<double>& phi,
                                       SolveReport* report,
                                       const Options& opt,
                                       const Solvers::Poisson::CutCellMetrics* cutMetrics,
                                       const PicEmbeddedBC* embeddedBC,
                                       const PicEmbeddedBCMap* embeddedBCById) const
{
    using Clock = std::chrono::steady_clock;
    const auto tBegin = Clock::now();
    auto secondsSince = [&](const Clock::time_point& t0) -> double {
        return std::chrono::duration_cast<std::chrono::duration<double>>(Clock::now() - t0).count();
    };

    double tAssemble_s = 0.0;
    double tFactor_s = 0.0;
    double tSolve_s = 0.0;

    // Assembly pipeline:
    // - Assemble interior contributions
    // - Apply domain BC (Dirichlet / Neumann) from pic.dict
    // - Detect pure-Neumann and apply compatibility + gauge fix
    // - Solve using Eigen::SparseLU (robust baseline)

    // Ensure phi has correct storage.
    if (phi.size() != mesh.cellIndexNum() || !phi.value)
    {
        phi.initialize(mesh, 0.0);
    }

    // ---------------------------------------------------------------------
    // The legacy solve() is now a cache-aware wrapper.
    // solveWithCache() owns the full cache policy
    // (hit/miss, rebuildEvery, disabled), so solve() simply delegates to it
    // and only falls back to the legacy one-shot path if the cached path fails.
    // ---------------------------------------------------------------------
    {
        SolveReport repSolve;
        const bool ok = solveWithCache(mesh, graph, bc, rhs, phi, &repSolve, opt, cutMetrics, embeddedBC, embeddedBCById);
        if (ok)
        {
            if (report)
            {
                *report = repSolve;
                // Include wrapper overhead in the reported total time.
                report->timeTotal_s = secondsSince(tBegin);
            }
            return true;
        }

        // Cached path failed unexpectedly: drop cache and fall back to legacy one-shot solve below.
        clearCache();
    }

    // ---------------------------------------------------------------------
    // Build active unknown set from FaceSegmentGraph offsets.
    // A cell index is considered active if it has at least one segment.
    // ---------------------------------------------------------------------
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

    if (N == 0)
    {
        if (report)
        {
            report->success = false;
            report->iterations = 0;
            report->finalResidual = 0.0;
            report->message = "No active cells in FaceSegmentGraph (N=0).";
            report->timeTotal_s = secondsSince(tBegin);
            report->timeAssemble_s = report->timeTotal_s;
            report->timeFactor_s = 0.0;
            report->timeSolve_s = 0.0;
            report->N = 0;
            report->nnz = 0;
            report->usedSuperLU = false;
        }
        if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
            std::cout << "[PoissonFVSolverCellCentered] N=0, nothing to solve.\n";
        return false;
    }

    // Optional cut-cell metrics:
    // - beta_f multiplies face conductance (beta*Af/df)
    // - alpha_K multiplies cell volume integration of rhs (rhs*(alpha*V))
    const bool useCutMetrics =
        (cutMetrics != nullptr) &&
        (cutMetrics->betaBySegment().size() == graph.segments().size()) &&
        (cutMetrics->cellMetricsByIndex().size() == static_cast<size_t>(mesh.cellIndexNum()));

    // Embedded-boundary (EB) closure requires both cut metrics and
    // an enabled embedded BC (parsed from pic.dict: BC/embedded).
    const bool useEmbeddedBCMap =
        useCutMetrics &&
        (embeddedBCById != nullptr) &&
        (!embeddedBCById->empty());
    const bool useEmbeddedBCFallback =
        useCutMetrics &&
        (embeddedBC != nullptr) &&
        embeddedBC->enabled;
    const bool useEmbeddedBC = useEmbeddedBCMap || useEmbeddedBCFallback;

    if (cutMetrics && !useCutMetrics && ConsoleControl::enabled(ConsoleControl::Level::kDebug))
        std::cout << "[PoissonFVSolverCellCentered] cutMetrics provided but size mismatch; ignoring.\n";

    // Cache per-row cell geometric area and plasma volume (alpha*area).
    std::vector<double> areaOfRow(static_cast<size_t>(N), 0.0);
    std::vector<double> volOfRow(static_cast<size_t>(N), 0.0);

    double alphaMin = std::numeric_limits<double>::infinity();
    double alphaMax = 0.0;
    size_t nAlphaClamped = 0u;
    size_t nAlphaNonFinite = 0u;

    for (int i = 0; i < N; ++i)
    {
        const int cellIndex = cellIndexOfRow[i];
        const double area = safeCellArea(mesh, cellIndex);
        areaOfRow[static_cast<size_t>(i)] = area;

        double alpha = 1.0;
        if (useCutMetrics)
        {
            alpha = cutMetrics->cellMetricsByIndex()[static_cast<size_t>(cellIndex)].alpha;
            if (!std::isfinite(alpha))
            {
                alpha = 1.0;
                ++nAlphaNonFinite;
            }
            if (alpha < 0.0)
            {
                alpha = 0.0;
                ++nAlphaClamped;
            }
            else if (alpha > 1.0)
            {
                alpha = 1.0;
                ++nAlphaClamped;
            }
        }

        volOfRow[static_cast<size_t>(i)] = area * alpha;
        alphaMin = std::min(alphaMin, alpha);
        alphaMax = std::max(alphaMax, alpha);
    }

    if (!useCutMetrics)
    {
        alphaMin = 1.0;
        alphaMax = 1.0;
    }

    // ---------------------------------------------------------------------
    // Identify cells with degenerate volume fraction for extrapolation.
    // Instead of assembling a near-singular flux equation, these cells
    // will receive a Laplacian averaging stencil: phi_i = mean(phi_nbrs).
    // ---------------------------------------------------------------------
    std::vector<bool> extrapolateRow(static_cast<size_t>(N), false);
    size_t nExtrapolate = 0u;
    if (useCutMetrics && opt.alphaMinThreshold > 0.0)
    {
        for (int i = 0; i < N; ++i)
        {
            const int cellIndex = cellIndexOfRow[static_cast<size_t>(i)];
            const double alpha = cutMetrics->cellMetricsByIndex()[static_cast<size_t>(cellIndex)].alpha;
            if (std::isfinite(alpha) && alpha < opt.alphaMinThreshold)
            {
                extrapolateRow[static_cast<size_t>(i)] = true;
                volOfRow[static_cast<size_t>(i)] = 0.0;
                ++nExtrapolate;
            }
        }
    }

    // Initialize RHS vector b with cell-integrated rhs: b_i = rhs * (alpha*area).
    Eigen::VectorXd bVec(N);
    bVec.setZero();
    const bool rhsOk = (rhs.size() == mesh.cellIndexNum()) && (rhs.value != nullptr);
    for (int i = 0; i < N; ++i)
    {
        const int cellIndex = cellIndexOfRow[i];
        const double f = rhsOk ? rhs.value[cellIndex] : 0.0;
        bVec[i] = f * volOfRow[static_cast<size_t>(i)];
    }

    // ---------------------------------------------------------------------
    // Assemble sparse matrix (triplets) and diagonal diagnostics.
    // For each oriented face-segment owned by K:
    //   interior: c = Af/df, A(i,i)+=c, A(i,j)-=c
    //   boundary Dirichlet: A(i,i)+=c, b(i)+=c*g
    //   boundary Neumann: b(i)+=g*Af
    // ---------------------------------------------------------------------
    std::vector<Eigen::Triplet<double>> trips;
    trips.reserve(graph.segments().size() * 2u + 8u);

    std::vector<double> diag(static_cast<size_t>(N), 0.0);
    std::vector<std::vector<int>> neighbors(static_cast<size_t>(N));
    std::vector<unsigned char> anchoredRow(static_cast<size_t>(N), 0u);

    size_t nSegTotal = graph.segments().size();
    size_t nInteriorUsed = 0u;
    size_t nInteriorSkippedMissing = 0u;
    size_t nInteriorSkippedDegenerate = 0u;

    size_t nBoundaryDirichlet = 0u;
    size_t nBoundaryNeumann = 0u;
    size_t nBoundaryPeriodic = 0u;
    size_t nBoundarySkippedMissing = 0u;
    size_t nBoundarySkippedDegenerate = 0u;

    bool warnedRobin = false;
    bool warnedPeriodicMismatch = false;

    const auto& segs = graph.segments();

    // Track beta statistics if cut-metrics are used.
    double betaMin = std::numeric_limits<double>::infinity();
    double betaMax = 0.0;
    size_t nBetaClamped = 0u;
    size_t nBetaNonFinite = 0u;
    size_t nZeroOpenArea = 0u;

    for (size_t sid = 0; sid < segs.size(); ++sid)
    {
        const auto& seg = segs[sid];

        const int Kindex = mesh.cells(seg.K).index();
        if (Kindex < 0)
        {
            // Owner cell is not in index set.
            continue;
        }
        const int i = rowOfCellIndex[static_cast<size_t>(Kindex)];
        if (i < 0)
        {
            continue;
        }

        // Face open fraction beta_f (default=1).
        double beta = 1.0;
        if (useCutMetrics)
        {
            beta = cutMetrics->betaBySegment()[sid];
            if (!std::isfinite(beta))
            {
                beta = 1.0;
                ++nBetaNonFinite;
            }
            if (beta < 0.0)
            {
                beta = 0.0;
                ++nBetaClamped;
            }
            else if (beta > 1.0)
            {
                beta = 1.0;
                ++nBetaClamped;
            }
            betaMin = std::min(betaMin, beta);
            betaMax = std::max(betaMax, beta);
        }

        const double AfOpen = beta * seg.Af;

        if (!(seg.df > 0.0) || !(AfOpen >= 0.0))
        {
            // Degenerate metric.
            if (seg.boundary == Boundary::non)
                ++nInteriorSkippedDegenerate;
            else
                ++nBoundarySkippedDegenerate;
            continue;
        }

        if (AfOpen == 0.0)
        {
            // Fully blocked face segment (by EB).
            ++nZeroOpenArea;
            continue;
        }

        if (seg.boundary == Boundary::non)
        {
            // Interior segment
            if (seg.L == -1)
            {
                ++nInteriorSkippedMissing;
                continue;
            }

            const int Lindex = mesh.cells(seg.L).index();
            if (Lindex < 0)
            {
                ++nInteriorSkippedMissing;
                continue;
            }

            const int j = rowOfCellIndex[static_cast<size_t>(Lindex)];
            if (j < 0)
            {
                // Neighbor not in the active unknown set (future: masked regions)
                ++nInteriorSkippedMissing;
                continue;
            }

            const double c = AfOpen / seg.df;
            if (!extrapolateRow[static_cast<size_t>(i)])
            {
                trips.emplace_back(i, i, c);
                trips.emplace_back(i, j, -c);
                diag[static_cast<size_t>(i)] += c;
            }
            neighbors[static_cast<size_t>(i)].push_back(j);
            neighbors[static_cast<size_t>(j)].push_back(i);
            ++nInteriorUsed;
        }
        else
        {
            // Domain boundary segment
            const PicBoundarySide& s = bcSide(bc, seg.boundary);
            PicBCType t = s.type;

            if (t == PicBCType::Robin)
            {
                // Only Dirichlet/Neumann are supported on the domain boundary.
                // Treat Robin as Dirichlet(value) for now (EB Robin comes later).
                if (!warnedRobin)
                {
                    warnedRobin = true;
                    if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
                    {
                        std::cout << "[PoissonFVSolverCellCentered] Warning: Robin BC encountered on domain boundary; "
                                  << "treating as Dirichlet(value).\n";
                    }
                }
                t = PicBCType::Dirichlet;
            }

            if (t == PicBCType::Dirichlet)
            {
                if (!extrapolateRow[static_cast<size_t>(i)])
                {
                    const double c = AfOpen / seg.df;
                    trips.emplace_back(i, i, c);
                    diag[static_cast<size_t>(i)] += c;
                    bVec[i] += c * s.value;
                    anchoredRow[static_cast<size_t>(i)] = 1u;
                }
                ++nBoundaryDirichlet;
            }
            else if (t == PicBCType::Neumann)
            {
                // Neumann specifies dphi/dn = g along the outward normal of the domain.
                // In the FV balance for -Delta(phi)=rhs, this contributes: b += g * (beta*Af).
                if (!extrapolateRow[static_cast<size_t>(i)])
                    bVec[i] += s.value * AfOpen;
                ++nBoundaryNeumann;
            }
            else if (t == PicBCType::Periodic)
            {
                const Boundary mateB = PeriodicFaceUtils::mateBoundary(seg.boundary);
                if (!warnedPeriodicMismatch)
                {
                    if (mateB == Boundary::non || bcSide(bc, mateB).type != PicBCType::Periodic)
                    {
                        warnedPeriodicMismatch = true;
                        if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
                        {
                            std::cout << "[PoissonFVSolverCellCentered] Warning: Periodic BC on '"
                                      << boundaryToString(seg.boundary)
                                      << "' but mate side is not Periodic. Assembly will still wrap.\n";
                        }
                    }
                }

                const Vector2D shift = PeriodicFaceUtils::shiftToMate(mesh, seg.boundary);
                int Llabel = -1;
                if (!PeriodicFaceUtils::findMateCellLabel(mesh, seg, shift, Llabel))
                {
                    ++nBoundarySkippedMissing;
                    continue;
                }

                const int Lindex = mesh.cells(Llabel).index();
                if (Lindex < 0)
                {
                    ++nBoundarySkippedMissing;
                    continue;
                }

                const int j = rowOfCellIndex[static_cast<size_t>(Lindex)];
                if (j < 0)
                {
                    ++nBoundarySkippedMissing;
                    continue;
                }

                const double df = PeriodicFaceUtils::centerDistance(mesh, seg, Llabel, shift);
                if (!(df > 0.0))
                {
                    ++nBoundarySkippedDegenerate;
                    continue;
                }

                const double c = AfOpen / df;
                if (!extrapolateRow[static_cast<size_t>(i)])
                {
                    trips.emplace_back(i, i, c);
                    trips.emplace_back(i, j, -c);
                    diag[static_cast<size_t>(i)] += c;
                }
                neighbors[static_cast<size_t>(i)].push_back(j);
                neighbors[static_cast<size_t>(j)].push_back(i);
                ++nBoundaryPeriodic;
            }
            else
            {
                // Unknown type: skip.
            }
        }
    }

    if (!useCutMetrics)
    {
        betaMin = 1.0;
        betaMax = 1.0;
    }

    // ---------------------------------------------------------------------
    // Embedded-boundary Robin/Dirichlet/Neumann closure on cut cells.
    // We use CutCellMetrics-provided per-cell metrics:
    //   A_gamma[K], d_gamma[K], n_gamma[K] (at least A_gamma)
    // and pic.dict embedded BC parameters (a,b,g).
    //
    // Robin form on EB: a*phi + b*dphi/dn = g
    // FV contribution for -Delta(phi)=rhs:
    //   A_kk += (a/b) * A_gamma
    //   b_k  += (g/b) * A_gamma
    // Dirichlet limit (b=0): a*phi=g -> phiD=g/a, approximated as
    //   flux ~ (phi_K - phiD)/d_gamma
    //   A_kk += (A_gamma/d_gamma)
    //   b_k  += (A_gamma/d_gamma) * phiD
    // Neumann on EB: dphi/dn = g -> b_k += g * A_gamma
    // ---------------------------------------------------------------------
    size_t nEBCellsTouched = 0u;
    size_t nEBCellsConstrained = 0u; // Dirichlet-like (adds diagonal penalty)
    size_t nEBRobin = 0u;
    size_t nEBDirichlet = 0u;
    size_t nEBNeumann = 0u;
    size_t nEBSkipped = 0u;
    size_t nEBDegenerate = 0u;

    double ebAddedDiagMin = std::numeric_limits<double>::infinity();
    double ebAddedDiagMax = 0.0;

    if (useEmbeddedBC)
    {
        size_t nEBMapUsed = 0u;
        size_t nEBFallbackUsed = 0u;
        size_t nEBMissing = 0u;
        bool ebMissingFatal = false;
        std::string ebMissingMsg;

        for (int i = 0; i < N; ++i)
        {
            const int cellIndex = cellIndexOfRow[i];
            const auto& cm = cutMetrics->cellMetricsByIndex()[static_cast<size_t>(cellIndex)];

            // Only apply EB closure on cut cells with positive EB measure.
            if (!cm.isCut || !(cm.A_gamma > 0.0))
                continue;

            if (extrapolateRow[static_cast<size_t>(i)])
                continue;

            // Skip fully-solid / inactive-volume cells.
            if (!(volOfRow[static_cast<size_t>(i)] > 0.0))
                continue;

            ++nEBCellsTouched;

            const PicEmbeddedBC* eb = nullptr;
            if (useEmbeddedBCMap && cm.gammaBoundaryId > 0)
            {
                auto it = embeddedBCById->find(cm.gammaBoundaryId);
                if (it != embeddedBCById->end() && it->second.enabled)
                {
                    eb = &it->second;
                    ++nEBMapUsed;
                }
            }
            if (!eb && useEmbeddedBCFallback)
            {
                eb = embeddedBC;
                ++nEBFallbackUsed;
            }
            if (!eb)
            {
                ++nEBMissing;
                if (opt.embeddedBCMissingPolicy == Options::EmbeddedBCMissingPolicy::Error)
                {
                    ebMissingFatal = true;
                    ebMissingMsg = "Missing embedded BC for boundaryId " + std::to_string(cm.gammaBoundaryId) +
                                   " (no map entry and no fallback).";
                    break;
                }
                continue;
            }

            const double a = eb->a;
            const double b = eb->b;
            const double g = eb->g;
            const PicBCType t = eb->type;
            const double A_gamma = cm.A_gamma;

            if (t == PicBCType::Robin)
            {
                if (std::abs(b) > 0.0)
                {
                    const double addDiag = (a / b) * A_gamma;
                    const double addRhs = (g / b) * A_gamma;

                    if (!std::isfinite(addDiag) || !std::isfinite(addRhs))
                    {
                        ++nEBDegenerate;
                        continue;
                    }

                    if (addDiag != 0.0)
                    {
                        trips.emplace_back(i, i, addDiag);
                        diag[static_cast<size_t>(i)] += addDiag;
                        ebAddedDiagMin = std::min(ebAddedDiagMin, addDiag);
                        ebAddedDiagMax = std::max(ebAddedDiagMax, addDiag);
                        if (a != 0.0)
                        {
                            ++nEBCellsConstrained;
                            anchoredRow[static_cast<size_t>(i)] = 1u;
                        }
                    }
                    bVec[i] += addRhs;
                    ++nEBRobin;
                }
                else
                {
                    // Dirichlet limit of Robin: b=0 => a*phi=g.
                    if (std::abs(a) <= 0.0)
                    {
                        // Ill-posed (a=b=0). Skip with warning.
                        ++nEBDegenerate;
                        continue;
                    }
                    const double phiD = g / a;
                    const double d = stabilizedEmbeddedDistance(cm.d_gamma,
                                                                volOfRow[static_cast<size_t>(i)]);

                    const double c = A_gamma / d;
                    if (!std::isfinite(c) || c <= 0.0)
                    {
                        ++nEBDegenerate;
                        continue;
                    }

                    trips.emplace_back(i, i, c);
                    diag[static_cast<size_t>(i)] += c;
                    bVec[i] += c * phiD;

                    ebAddedDiagMin = std::min(ebAddedDiagMin, c);
                    ebAddedDiagMax = std::max(ebAddedDiagMax, c);
                    ++nEBCellsConstrained;
                    anchoredRow[static_cast<size_t>(i)] = 1u;
                    ++nEBDirichlet;
                }
            }
            else if (t == PicBCType::Dirichlet)
            {
                // Treat g as phi_D.
                const double phiD = g;
                const double d = stabilizedEmbeddedDistance(cm.d_gamma,
                                                            volOfRow[static_cast<size_t>(i)]);

                const double c = A_gamma / d;
                if (!std::isfinite(c) || c <= 0.0)
                {
                    ++nEBDegenerate;
                    continue;
                }

                trips.emplace_back(i, i, c);
                diag[static_cast<size_t>(i)] += c;
                bVec[i] += c * phiD;

                ebAddedDiagMin = std::min(ebAddedDiagMin, c);
                ebAddedDiagMax = std::max(ebAddedDiagMax, c);
                ++nEBCellsConstrained;
                anchoredRow[static_cast<size_t>(i)] = 1u;
                ++nEBDirichlet;
            }
            else if (t == PicBCType::Neumann)
            {
                // Treat g as dphi/dn on EB (outward normal from plasma).
                const double addRhs = g * A_gamma;
                if (!std::isfinite(addRhs))
                {
                    ++nEBDegenerate;
                    continue;
                }
                bVec[i] += addRhs;
                ++nEBNeumann;
            }
            else
            {
                ++nEBSkipped;
            }
        }

        if (ebMissingFatal)
        {
            if (report)
            {
                report->success = false;
                report->iterations = 0;
                report->finalResidual = 0.0;
                report->message = ebMissingMsg;
                report->timeTotal_s = secondsSince(tBegin);
                report->timeAssemble_s = report->timeTotal_s;
                report->timeFactor_s = 0.0;
                report->timeSolve_s = 0.0;
                report->N = N;
                report->nnz = 0;
                report->usedSuperLU = false;
            }
            return false;
        }
        (void)nEBMapUsed;
        (void)nEBFallbackUsed;
        (void)nEBMissing;
    }

    // ---------------------------------------------------------------------
    // Extrapolation stencil for small-alpha cells (legacy solve path).
    // Replace degenerate flux equations with phi_i = mean(phi_neighbors).
    // ---------------------------------------------------------------------
    if (nExtrapolate > 0u)
    {
        for (int i = 0; i < N; ++i)
        {
            if (!extrapolateRow[static_cast<size_t>(i)])
                continue;

            auto& nbrs = neighbors[static_cast<size_t>(i)];
            std::sort(nbrs.begin(), nbrs.end());
            nbrs.erase(std::unique(nbrs.begin(), nbrs.end()), nbrs.end());
            nbrs.erase(std::remove(nbrs.begin(), nbrs.end(), i), nbrs.end());

            const int count = static_cast<int>(nbrs.size());
            if (count > 0)
            {
                trips.emplace_back(i, i, static_cast<double>(count));
                for (int j : nbrs)
                    trips.emplace_back(i, j, -1.0);
                diag[static_cast<size_t>(i)] = static_cast<double>(count);
            }
            else
            {
                trips.emplace_back(i, i, 1.0);
                diag[static_cast<size_t>(i)] = 1.0;
            }
            bVec[i] = 0.0;
            anchoredRow[static_cast<size_t>(i)] = 1u;
        }

        if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
        {
            std::cout << "[PoissonFVSolverCellCentered] solve(): "
                      << nExtrapolate
                      << " cell(s) with alpha < " << opt.alphaMinThreshold
                      << " replaced with Laplacian averaging stencil.\n";
        }
    }

    // ---------------------------------------------------------------------
    // Detect pure Neumann and enforce compatibility + gauge fix.
    // "Pure Neumann" here means: (i) domain boundary has no Dirichlet-like
    // constraint AND (ii) embedded boundary closure does not add any
    // Dirichlet-like diagonal penalty (i.e., no constrained EB cells).
    // ---------------------------------------------------------------------
    const bool domainPureNeumann =
        (!isDirichletLike(bc.left.type)) &&
        (!isDirichletLike(bc.right.type)) &&
        (!isDirichletLike(bc.bottom.type)) &&
        (!isDirichletLike(bc.top.type));

    const bool embeddedAnchors = (useEmbeddedBC && (nEBCellsConstrained > 0u));

    const bool pureNeumann = domainPureNeumann && (!embeddedAnchors);
    std::vector<int> pinnedRows;
    size_t nFloatingComponents = 0u;

    if (pureNeumann)
    {
        // Enforce compatibility: sum(b) = 0 (since interior fluxes cancel globally).
        // With cut metrics, b already contains cell-integrated rhs over plasma volume (alpha*V).
        double totalB = 0.0;
        double totalV = 0.0;
        for (int i = 0; i < N; ++i)
        {
            totalB += bVec[i];
            totalV += volOfRow[static_cast<size_t>(i)];
        }

        if (totalV > 0.0)
        {
            const double meanBperV = totalB / totalV;
            for (int i = 0; i < N; ++i)
            {
                bVec[i] -= meanBperV * volOfRow[static_cast<size_t>(i)];
            }
        }

        // Gauge fix: pin one DOF (phi_pin = 0) by adding a unit diagonal.
        // Prefer a row with positive plasma volume to avoid pinning a fully-solid cell.
        int pinRow = 0;
        for (int i = 0; i < N; ++i)
        {
            if (volOfRow[static_cast<size_t>(i)] > 0.0)
            {
                pinRow = i;
                break;
            }
        }

        trips.emplace_back(pinRow, pinRow, 1.0);
        diag[static_cast<size_t>(pinRow)] += 1.0;
        bVec[pinRow] = 0.0;
        pinnedRows.push_back(pinRow);

        if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
            std::cout << "[PoissonFVSolverCellCentered] Pure Neumann detected: applied compatibility correction and gauge fix (pin row "
                      << pinRow << ")\n";
    }
    else
    {
        nFloatingComponents = stabilizeFloatingComponents(neighbors,
                                                          anchoredRow,
                                                          volOfRow,
                                                          bVec,
                                                          trips,
                                                          diag,
                                                          &pinnedRows,
                                                          nullptr);

        if ((nFloatingComponents > 0u) && ConsoleControl::enabled(ConsoleControl::Level::kDebug))
        {
            std::cout << "[PoissonFVSolverCellCentered] stabilized "
                      << nFloatingComponents
                      << " floating component(s) with compatibility correction + gauge pin.\n";
        }
    }

    // Safety: pin rows with zero/invalid diagonal to avoid singular factorization.
    size_t nPinnedZeroDiag = 0u;
    for (int i = 0; i < N; ++i)
    {
        const double d = diag[static_cast<size_t>(i)];
        if (!std::isfinite(d) || !(d > 0.0))
        {
            trips.emplace_back(i, i, 1.0);
            diag[static_cast<size_t>(i)] += 1.0;
            bVec[i] = 0.0;
            pinnedRows.push_back(i);
            ++nPinnedZeroDiag;
        }
    }

    // Finalize sparse matrix.
    Eigen::SparseMatrix<double> A(N, N);
    A.setFromTriplets(trips.begin(), trips.end());
    A.makeCompressed();

    tAssemble_s = secondsSince(tBegin);
    if (report)
    {
        report->timeAssemble_s = tAssemble_s;
        report->N = N;
        report->nnz = static_cast<long long>(A.nonZeros());
    }

    // Diagonal diagnostics.
    double minDiag = std::numeric_limits<double>::infinity();
    double maxDiag = 0.0;
    for (double d : diag)
    {
        minDiag = std::min(minDiag, d);
        maxDiag = std::max(maxDiag, d);
    }

    const bool consoleDebug = ConsoleControl::enabled(ConsoleControl::Level::kDebug);
    const bool consoleTrace = ConsoleControl::enabled(ConsoleControl::Level::kTrace);

    if (consoleDebug)
    {
        std::cout << "[PoissonFVSolverCellCentered] Assembled (interior+domainBC"
                  << (useEmbeddedBC ? "+EB" : "")
                  << "): "
                  << "N=" << N << ", nnz=" << A.nonZeros()
                  << ", diagMin=" << minDiag << ", diagMax=" << maxDiag
                  << (pureNeumann ? " (pure Neumann + gauge)" : "")
                  << "\n";

        if (useCutMetrics)
        {
            std::cout << "  cutMetrics=on"
                      << ", alphaMin=" << alphaMin << ", alphaMax=" << alphaMax
                      << ", betaMin=" << betaMin << ", betaMax=" << betaMax
                      << ", zeroOpenSeg=" << nZeroOpenArea
                      << ", betaClamped=" << nBetaClamped
                      << ", betaNonFinite=" << nBetaNonFinite
                      << ", alphaClamped=" << nAlphaClamped
                      << ", alphaNonFinite=" << nAlphaNonFinite
                      << "\n";
        }

        if (useEmbeddedBC)
        {
            std::cout << "  embeddedBC=on"
                      << ", type=" << bcTypeName(embeddedBC->type)
                      << ", a=" << embeddedBC->a
                      << ", b=" << embeddedBC->b
                      << ", g=" << embeddedBC->g
                      << ", ebTouched=" << nEBCellsTouched
                      << ", ebConstrained=" << nEBCellsConstrained
                      << ", ebRobin=" << nEBRobin
                      << ", ebDir=" << nEBDirichlet
                      << ", ebNeu=" << nEBNeumann
                      << ", ebSkipped=" << nEBSkipped
                      << ", ebDegenerate=" << nEBDegenerate;

            if (std::isfinite(ebAddedDiagMin))
                std::cout << ", ebDiagMin=" << ebAddedDiagMin << ", ebDiagMax=" << ebAddedDiagMax;
            else
                std::cout << ", ebDiagMin=0, ebDiagMax=0";

            std::cout << ", embeddedAnchors=" << (embeddedAnchors ? "true" : "false")
                      << "\n";
        }

        if (nPinnedZeroDiag > 0u)
        {
            std::cout << "  pinned zero-diag rows=" << nPinnedZeroDiag << "\n";
        }
        if (nFloatingComponents > 0u)
        {
            std::cout << "  floating components stabilized=" << nFloatingComponents << "\n";
        }

        if (consoleTrace)
        {
            std::cout << "  segments total                = " << nSegTotal << "\n";
            std::cout << "  interior used                 = " << nInteriorUsed << "\n";
            std::cout << "  interior skipped (L missing/inactive) = " << nInteriorSkippedMissing << "\n";
            std::cout << "  interior skipped (degenerate df/Af)   = " << nInteriorSkippedDegenerate << "\n";
            std::cout << "  boundary Dirichlet segments    = " << nBoundaryDirichlet << "\n";
            std::cout << "  boundary Neumann segments      = " << nBoundaryNeumann << "\n";
            std::cout << "  boundary Periodic segments     = " << nBoundaryPeriodic << "\n";
            std::cout << "  boundary skipped (mate missing/inactive) = " << nBoundarySkippedMissing << "\n";
            std::cout << "  boundary skipped (degenerate df/Af)   = " << nBoundarySkippedDegenerate << "\n";
            std::cout << "  BC types: left=" << bcTypeName(bc.left.type)
                      << ", right=" << bcTypeName(bc.right.type)
                      << ", bottom=" << bcTypeName(bc.bottom.type)
                      << ", top=" << bcTypeName(bc.top.type) << "\n";
        }
    }

    // ---------------------------------------------------------------------
    // Solve using SuperLU (direct) if available; otherwise fall back to SparseLU.
    // Notes:
    // - maxIter/relTol are kept in Options for future CG-based solvers.
    // - For direct solvers, iterations/residual are not tracked.
    // ---------------------------------------------------------------------
    Eigen::VectorXd x;
    bool usedSuperLU = false;
    bool usedLDLT = false;
    bool usedBiCGSTAB = false;
    int iterativeIterations = 0;
    double iterativeResidual = 0.0;

    auto dumpMatrixDiagnostics = [&](bool force = false) {
        if (!force && !ConsoleControl::enabled(ConsoleControl::Level::kDebug))
            return;
        // Basic structural diagnostics helpful when SuperLU reports a singular factorization.
        // We use the pre-assembly diag[] as a proxy for missing couplings.
        size_t nZeroDiag = 0u;
        size_t nNonFiniteDiag = 0u;
        size_t nTinyDiag1e12 = 0u;
        size_t nTinyDiag1e10 = 0u;
        size_t nTinyDiag1e8 = 0u;
        double minPositiveDiag = std::numeric_limits<double>::infinity();
        for (double d : diag)
        {
            if (!std::isfinite(d))
            {
                ++nNonFiniteDiag;
                continue;
            }
            const double ad = std::abs(d);
            if (ad == 0.0)
            {
                ++nZeroDiag;
                continue;
            }
            minPositiveDiag = std::min(minPositiveDiag, ad);
            if (ad < 1e-12)
                ++nTinyDiag1e12;
            if (ad < 1e-10)
                ++nTinyDiag1e10;
            if (ad < 1e-8)
                ++nTinyDiag1e8;
        }

        size_t nNonFiniteB = 0u;
        for (int i = 0; i < bVec.size(); ++i)
        {
            if (!std::isfinite(bVec[i]))
                ++nNonFiniteB;
        }

        std::cout << "[PoissonFVSolverCellCentered] Diagnostics:"
                  << " zeroDiagRows=" << nZeroDiag
                  << ", nonFiniteDiagRows=" << nNonFiniteDiag
                  << ", tinyDiag<1e-12=" << nTinyDiag1e12
                  << ", tinyDiag<1e-10=" << nTinyDiag1e10
                  << ", tinyDiag<1e-8=" << nTinyDiag1e8
                  << ", minPositiveDiag=" << (std::isfinite(minPositiveDiag) ? minPositiveDiag : 0.0)
                  << ", nonFiniteBEntries=" << nNonFiniteB
                  << ", pureNeumann=" << (pureNeumann ? "true" : "false")
                  << ", floatingComponents=" << nFloatingComponents
                  << ", pinnedZeroDiagRows=" << nPinnedZeroDiag
                  << ", ebConstrainedRows=" << nEBCellsConstrained
                  << ", embeddedAnchors=" << (embeddedAnchors ? "true" : "false")
                  << "\n";
        if (nZeroDiag > 0u && ConsoleControl::enabled(ConsoleControl::Level::kTrace))
            std::cout << "  Hint: zero diagonal often means an active cell has no interior couplings and no Dirichlet penalty terms (e.g. isolated active cell or degenerate df/Af).\n";
    };

    auto tryLDLTFallback = [&]() -> bool {
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
        const auto tFac0 = Clock::now();
        solver.compute(A);
        tFactor_s += std::chrono::duration_cast<std::chrono::duration<double>>(Clock::now() - tFac0).count();
        if (solver.info() != Eigen::Success)
        {
            if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
                std::cout << "[PoissonFVSolverCellCentered] SimplicialLDLT factorization failed.\n";
            return false;
        }

        const auto tSol0 = Clock::now();
        x = solver.solve(bVec);
        tSolve_s += std::chrono::duration_cast<std::chrono::duration<double>>(Clock::now() - tSol0).count();
        if (solver.info() != Eigen::Success)
        {
            if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
                std::cout << "[PoissonFVSolverCellCentered] SimplicialLDLT solve failed.\n";
            return false;
        }

        for (int i = 0; i < x.size(); ++i)
        {
            if (!std::isfinite(x[i]))
                return false;
        }

        usedLDLT = true;
        return true;
    };

    auto tryBiCGSTABFallback = [&]() -> bool {
        Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>> solver;
        solver.setMaxIterations(std::max(1000, opt.maxIter));
        solver.setTolerance((opt.relTol > 0.0) ? opt.relTol : 1e-10);
        solver.preconditioner().setDroptol(1e-3);
        solver.preconditioner().setFillfactor(10);

        const auto tFac0 = Clock::now();
        solver.compute(A);
        tFactor_s += std::chrono::duration_cast<std::chrono::duration<double>>(Clock::now() - tFac0).count();
        if (solver.info() != Eigen::Success)
        {
            if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
                std::cout << "[PoissonFVSolverCellCentered] BiCGSTAB preconditioner setup failed.\n";
            return false;
        }

        const auto tSol0 = Clock::now();
        x = solver.solve(bVec);
        tSolve_s += std::chrono::duration_cast<std::chrono::duration<double>>(Clock::now() - tSol0).count();
        iterativeIterations = static_cast<int>(solver.iterations());
        iterativeResidual = solver.error();
        if (solver.info() != Eigen::Success)
        {
            if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
                std::cout << "[PoissonFVSolverCellCentered] BiCGSTAB solve failed after "
                          << iterativeIterations << " iterations, error=" << iterativeResidual << ".\n";
            return false;
        }

        for (int i = 0; i < x.size(); ++i)
        {
            if (!std::isfinite(x[i]))
                return false;
        }

        usedBiCGSTAB = true;
        return true;
    };

#ifdef EIGEN_SUPERLU_SUPPORT
    if (opt.useSuperLU)
    {
        Eigen::SuperLU<Eigen::SparseMatrix<double>> solver;

        // SuperLU can be sensitive to ordering on ill-conditioned / nearly singular systems.
        // COLAMD is a robust default for general sparse matrices.
        solver.options().ColPerm = COLAMD;
        solver.options().DiagPivotThresh = 0.1;

        const auto tFac0 = Clock::now();
        solver.analyzePattern(A);
        solver.factorize(A);
        tFactor_s += std::chrono::duration_cast<std::chrono::duration<double>>(Clock::now() - tFac0).count();

        if (solver.info() != Eigen::Success)
        {
            if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
                std::cout << "[PoissonFVSolverCellCentered] SuperLU factorization failed; falling back to SparseLU.\n";
            dumpMatrixDiagnostics(/*force=*/true);
        }
        else
        {
            const auto tSol0 = Clock::now();
            x = solver.solve(bVec);
            tSolve_s += std::chrono::duration_cast<std::chrono::duration<double>>(Clock::now() - tSol0).count();

            if (solver.info() != Eigen::Success)
            {
                if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
                    std::cout << "[PoissonFVSolverCellCentered] SuperLU solve failed; falling back to SparseLU.\n";
                dumpMatrixDiagnostics(/*force=*/true);
            }
            else
            {
                usedSuperLU = true;
            }
        }
    }
#else
    if (opt.useSuperLU)
    {
        if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
            std::cout << "[PoissonFVSolverCellCentered] SuperLU requested but EIGEN_SUPERLU_SUPPORT is not enabled at compile time; "
                      << "falling back to SparseLU. (Install/link SuperLU and reconfigure.)\n";
    }
#endif

    if (!usedSuperLU)
    {
        if (opt.useSuperLU && ConsoleControl::enabled(ConsoleControl::Level::kDebug))
            std::cout << "[PoissonFVSolverCellCentered] Using SparseLU fallback.\n";

        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;

        const auto tFac0 = Clock::now();
        solver.analyzePattern(A);
        solver.factorize(A);
        tFactor_s += std::chrono::duration_cast<std::chrono::duration<double>>(Clock::now() - tFac0).count();

        if (solver.info() != Eigen::Success)
        {
            dumpMatrixDiagnostics(/*force=*/true);
            if (!tryLDLTFallback() && !tryBiCGSTABFallback())
            {
                std::cerr << "[PoissonFVSolverCellCentered] SparseLU factorization failed.\n";
                if (report)
                {
                    report->success = false;
                    report->iterations = 0;
                    report->finalResidual = 0.0;
                    report->message = "SparseLU factorization failed.";
                    report->timeAssemble_s = tAssemble_s;
                    report->timeFactor_s = tFactor_s;
                    report->timeSolve_s = tSolve_s;
                    report->timeTotal_s = secondsSince(tBegin);
                    report->N = N;
                    report->nnz = static_cast<long long>(A.nonZeros());
                    report->usedSuperLU = false;
                }
                return false;
            }
            std::cout << "[PoissonFVSolverCellCentered] SparseLU factorization failed; recovered with "
                      << (usedLDLT ? "SimplicialLDLT" : "BiCGSTAB")
                      << " fallback.\n";
        }
        else
        {
            const auto tSol0 = Clock::now();
            x = solver.solve(bVec);
            tSolve_s += std::chrono::duration_cast<std::chrono::duration<double>>(Clock::now() - tSol0).count();

            if (solver.info() != Eigen::Success)
            {
                dumpMatrixDiagnostics(/*force=*/true);
                if (!tryLDLTFallback() && !tryBiCGSTABFallback())
                {
                    std::cerr << "[PoissonFVSolverCellCentered] SparseLU solve failed.\n";
                    if (report)
                    {
                        report->success = false;
                        report->iterations = 0;
                        report->finalResidual = 0.0;
                        report->message = "SparseLU solve failed.";
                        report->timeAssemble_s = tAssemble_s;
                        report->timeFactor_s = tFactor_s;
                        report->timeSolve_s = tSolve_s;
                        report->timeTotal_s = secondsSince(tBegin);
                        report->N = N;
                        report->nnz = static_cast<long long>(A.nonZeros());
                        report->usedSuperLU = false;
                    }
                    return false;
                }
                std::cout << "[PoissonFVSolverCellCentered] SparseLU solve failed; recovered with "
                          << (usedLDLT ? "SimplicialLDLT" : "BiCGSTAB")
                          << " fallback.\n";
            }
        }
    }

    if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
        std::cout << "[PoissonFVSolverCellCentered] Solved using "
                  << (usedSuperLU ? "SuperLU" : (usedLDLT ? "SimplicialLDLT" : (usedBiCGSTAB ? "BiCGSTAB" : "SparseLU")))
                  << ".\n";

    // Scatter solution back to phi (cell-index aligned).

    for (int i = 0; i < N; ++i)
    {
        const int cellIndex = cellIndexOfRow[i];
        phi.value[cellIndex] = x[i];
    }

    if (report)
    {
        report->success = true;
        report->iterations = usedBiCGSTAB ? iterativeIterations : 0;
        report->finalResidual = usedBiCGSTAB ? iterativeResidual : 0.0;
        report->message = usedSuperLU ? "Solved using SuperLU." : (usedLDLT ? "Solved using SimplicialLDLT fallback." : (usedBiCGSTAB ? "Solved using BiCGSTAB fallback." : "Solved using SparseLU."));

        report->timeAssemble_s = tAssemble_s;
        report->timeFactor_s = tFactor_s;
        report->timeSolve_s = tSolve_s;
        report->timeTotal_s = secondsSince(tBegin);
        report->N = N;
        report->nnz = static_cast<long long>(A.nonZeros());
        report->usedSuperLU = usedSuperLU;
    }

    return true;
}

PoissonFVSolverCellCentered::ErrorNorms
PoissonFVSolverCellCentered::regressionUnitBoxSinSin(const QuadtreeMesh& mesh,
                                                     const FaceSegmentGraph& graph,
                                                     CellField<double>& phiOut,
                                                     const Options& opt) const
{
    // Dirichlet (phi=0) on the entire domain boundary.
    PicDomainBC bc;
    bc.left.type = PicBCType::Dirichlet;
    bc.right.type = PicBCType::Dirichlet;
    bc.bottom.type = PicBCType::Dirichlet;
    bc.top.type = PicBCType::Dirichlet;

    bc.left.value = 0.0;
    bc.right.value = 0.0;
    bc.bottom.value = 0.0;
    bc.top.value = 0.0;

    const double Lx = mesh.maxCoord.x - mesh.minCoord.x;
    const double Ly = mesh.maxCoord.y - mesh.minCoord.y;
    if (Lx <= 0.0 || Ly <= 0.0)
    {
        if (ConsoleControl::enabled(ConsoleControl::Level::kInfo))
            std::cout << "[PoissonFVSolverCellCentered] Regression skipped: invalid domain size." << std::endl;
        return ErrorNorms{};
    }

    const double pi = std::acos(-1.0);
    const double ax = pi / Lx;
    const double ay = pi / Ly;
    const double coeff = ax * ax + ay * ay;

    // Build RHS on active cells.
    CellField<double> rhs;
    rhs.initialize(mesh, 0.0);
    phiOut.initialize(mesh, 0.0);

    const int nCellIndex = mesh.cellIndexNum();
    size_t nActive = 0;
    for (int ci = 0; ci < nCellIndex; ++ci)
    {
        const auto range = graph.segmentRangeForCellIndex(ci);
        if (range.first == range.second)
            continue;

        const int K = mesh.cellIndexToLabel[ci];
        const Cell& cell = mesh.cells(K);
        const Vector2D c = cell.centerCoord();

        const double xhat = (c.x - mesh.minCoord.x) / Lx;
        const double yhat = (c.y - mesh.minCoord.y) / Ly;
        const double phiExact = std::sin(pi * xhat) * std::sin(pi * yhat);

        rhs.value[ci] = coeff * phiExact;
        ++nActive;
    }

    // Solve
    SolveReport rep;
    CellField<double> phi;
    phi.initialize(mesh, 0.0);

    const bool ok = solve(mesh, graph, bc, rhs, phi, &rep, opt);
    if (!ok)
    {
        if (ConsoleControl::enabled(ConsoleControl::Level::kInfo))
        if (ConsoleControl::enabled(ConsoleControl::Level::kInfo))
            std::cout << "[PoissonFVSolverCellCentered] Regression solve failed: " << rep.message << std::endl;
        phiOut = phi;
        return ErrorNorms{};
    }

    // Error norms on active cells.
    double sumE2V = 0.0;
    double sumV = 0.0;
    double linf = 0.0;

    for (int ci = 0; ci < nCellIndex; ++ci)
    {
        const auto range = graph.segmentRangeForCellIndex(ci);
        if (range.first == range.second)
            continue;

        const int K = mesh.cellIndexToLabel[ci];
        const Cell& cell = mesh.cells(K);
        const Vector2D c = cell.centerCoord();

        const double xhat = (c.x - mesh.minCoord.x) / Lx;
        const double yhat = (c.y - mesh.minCoord.y) / Ly;
        const double phiExact = std::sin(pi * xhat) * std::sin(pi * yhat);

        const double e = phi.value[ci] - phiExact;
        const double V = cell.area();

        sumE2V += e * e * V;
        sumV += V;
        linf = std::max(linf, std::abs(e));
    }

    ErrorNorms norms;
    norms.L2 = (sumV > 0.0) ? std::sqrt(sumE2V / sumV) : 0.0;
    norms.Linf = linf;

    if (ConsoleControl::enabled(ConsoleControl::Level::kInfo))
    {
        std::cout << "[PoissonFVSolverCellCentered] Unit-box sin-sin regression:" << std::endl;
        std::cout << "  activeCells = " << nActive << std::endl;
        std::cout << "  L2   = " << norms.L2 << std::endl;
        std::cout << "  Linf = " << norms.Linf << std::endl;
    }

    // Return solution
    phiOut = phi;
    return norms;
}


PoissonFVSolverCellCentered::ErrorNorms PoissonFVSolverCellCentered::regressionUnitBoxSinCosPeriodicY(
    const QuadtreeMesh& mesh,
    const FaceSegmentGraph& graph,
    CellField<double>& phiOut,
    const Options& opt) const
{
    // Manufactured solution on a rectangular unit box with periodic-Y (bottom/top) and Dirichlet on left/right.
    // u(x,y) = sin(pi xhat) cos(2 pi yhat) is periodic in y and vanishes on xhat=0,1.

    const Vector2D x0 = mesh.minCoord;
    const Vector2D x1 = mesh.maxCoord;
    const double Lx = x1.x - x0.x;
    const double Ly = x1.y - x0.y;

    PicDomainBC bc;
    bc.left.type = PicBCType::Dirichlet;
    bc.right.type = PicBCType::Dirichlet;
    bc.left.value = 0.0;
    bc.right.value = 0.0;
    bc.bottom.type = PicBCType::Periodic;
    bc.top.type = PicBCType::Periodic;
    bc.bottom.value = 0.0;
    bc.top.value = 0.0;

    // Build RHS and exact solution on active cells.
    CellField<double> rhs;
    rhs.initialize(mesh, 0.0);

    CellField<double> phiExact;
    phiExact.initialize(mesh, 0.0);

    const double pi = std::acos(-1.0);
    const double ax = pi / Lx;
    const double ay = 2.0 * pi / Ly;
    const double k2 = ax * ax + ay * ay;

    for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
    {
        const int cLbl = mesh.cellIndexToLabel[ci];
        const Cell& c = mesh.cells(cLbl);
        const Vector2D xc = c.centerCoord();
        const double xhat = (xc.x - x0.x) / Lx;
        const double yhat = (xc.y - x0.y) / Ly;
        const double u = std::sin(pi * xhat) * std::cos(2.0 * pi * yhat);
        phiExact.value[ci] = u;
        rhs.value[ci] = k2 * u; // -Laplace(u) = k^2 u
    }

    // Solve.
    solve(mesh, graph, bc, rhs, phiOut, nullptr, opt);

    // Compute error norms.
    ErrorNorms norms;
    double sum2 = 0.0;
    double maxAbs = 0.0;
    double volSum = 0.0;

    for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
    {
        const int cLbl = mesh.cellIndexToLabel[ci];
        const Cell& c = mesh.cells(cLbl);
        const double e = phiOut.value[ci] - phiExact.value[ci];
        const double w = c.area();
        sum2 += e * e * w;
        volSum += w;
        maxAbs = std::max(maxAbs, std::abs(e));
    }

    if (volSum > 0.0)
        norms.L2 = std::sqrt(sum2 / volSum);
    norms.Linf = maxAbs;
    return norms;
}
