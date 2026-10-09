#pragma once

#include <memory>
#include <string>

#include <Eigen/Sparse>

#include "Core/Field/CellField.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Solvers/Poisson/PoissonBC.h"

namespace Solvers { namespace Poisson { class CutCellMetrics; } }

class QuadtreeMesh;

// ============================================================================
// Cell-centered FV Poisson solver on quadtree leaf-face segments
//
// This solver targets the explicit PIC module (cell-centered FV discretization).
//
// Governing form (no EB / cut terms yet):
//   -Delta(phi) = rhs
// with domain Dirichlet/Neumann/Robin BC.
//
// BC source:
//   The physics module supplies domain BC explicitly.
// ============================================================================
class PoissonFVSolverCellCentered
{
public:
    struct Options
    {
        // Linear solver controls
        // Prefer SuperLU if it is available at compile/link time; otherwise
        // fall back to Eigen::SparseLU.
        bool useSuperLU = true;

        // Iterative-solver knobs are kept for future CG-based solvers.
        // (SuperLU/SparseLU are direct solvers and will ignore these.)
        bool useCG = true;
        bool useIC = true;
        int maxIter = 2000;
        double relTol = 1e-10;

        // Assembly / runtime
        bool verbose = false;

        // -----------------------------------------------------------------
        // Factorization cache controls
        //
        // When enabled, this solver caches the assembled sparse matrix and
        // its direct-solver factorization between calls.
        //
        // Typical PIC usage: mesh/graph/BC are fixed, only RHS changes.
        // In that case, caching avoids repeated analyzePattern+factorize.
        // -----------------------------------------------------------------
        bool cacheEnabled = true;

        // Force periodic cache rebuild every N solve calls (0 disables).
        // This is mainly for debugging or when coefficients change slowly
        // without changing the cache signature.
        int cacheRebuildEvery = 0;

        // Cache logging:
        // - If verbose=true, log cache hit/miss every call.
        // - Else if cacheLogEvery>0, log a short hit/miss summary every N calls.
        int cacheLogEvery = 0;

        // Additionally log every cache miss event even when verbose=false.
        bool cacheLogOnMiss = true;

        // -----------------------------------------------------------------
        // Embedded-boundary (EB) BC lookup behavior
        // -----------------------------------------------------------------
        enum class EmbeddedBCMissingPolicy
        {
            Error, // abort solve/setup if a cut-cell has no matching EB BC
            Skip   // skip EB closure for that cut-cell
        };
        EmbeddedBCMissingPolicy embeddedBCMissingPolicy = EmbeddedBCMissingPolicy::Error;

        // -----------------------------------------------------------------
        // Small cut-cell treatment (alpha-threshold extrapolation)
        //
        // Cells with volume fraction alpha < alphaMinThreshold are treated
        // as "extrapolation cells": instead of assembling the degenerate
        // flux equation (which yields near-zero diagonal and phi~0), the
        // solver replaces the equation with a Laplacian averaging stencil:
        //     phi_i = mean(phi_neighbors)
        // This produces a smooth, physically reasonable field extension
        // through near-solid cells without excluding them from the matrix.
        //
        // Set to 0.0 (default) to disable.
        // -----------------------------------------------------------------
        double alphaMinThreshold = 0.0;
    };

    struct SolveReport
    {
        bool success = false;
        int iterations = 0;
        double finalResidual = 0.0;
        std::string message;

        // Timing diagnostics (wall-clock seconds).
        // For direct solvers, iterations/finalResidual may remain 0.
        double timeAssemble_s = 0.0; // build unknown map, rhs, triplets, A.setFromTriplets
        double timeFactor_s = 0.0;   // analyzePattern + factorize
        double timeSolve_s = 0.0;    // triangular solves (solver.solve)
        double timeTotal_s = 0.0;    // total wall time in solve()

        // Matrix diagnostics
        int N = 0;
        long long nnz = 0;
        bool usedSuperLU = false;
    };

    struct ErrorNorms
    {
        double L2 = 0.0;
        double Linf = 0.0;
    };

    PoissonFVSolverCellCentered();
    ~PoissonFVSolverCellCentered();

    PoissonFVSolverCellCentered(const PoissonFVSolverCellCentered&) = delete;
    PoissonFVSolverCellCentered& operator=(const PoissonFVSolverCellCentered&) = delete;

    PoissonFVSolverCellCentered(PoissonFVSolverCellCentered&&) noexcept = default;
    PoissonFVSolverCellCentered& operator=(PoissonFVSolverCellCentered&&) noexcept = default;

    // ---------------------------------------------------------------------
    // One-shot solve (legacy API).
    // The wrapper uses setup()+solveWithCache() and reuses cached factorization
    // when the cache signature matches (mesh/graph/BC/cut-metrics/EB config).
    // ---------------------------------------------------------------------
    bool solve(const QuadtreeMesh& mesh,
               const FaceSegmentGraph& graph,
               const PicDomainBC& bc,
               const CellField<double>& rhs,
               CellField<double>& phi,
               SolveReport* report = nullptr,
               const Options& opt = Options(),
               const Solvers::Poisson::CutCellMetrics* cutMetrics = nullptr,
               const PicEmbeddedBC* embeddedBC = nullptr,
               const PicEmbeddedBCMap* embeddedBCById = nullptr) const;

    // ---------------------------------------------------------------------
    // Two-stage API (cache skeleton).
    // ---------------------------------------------------------------------
    SolveReport setup(const QuadtreeMesh& mesh,
                      const FaceSegmentGraph& graph,
                      const PicDomainBC& bc,
                      const Options& opt = Options(),
                      const Solvers::Poisson::CutCellMetrics* cutMetrics = nullptr,
                      const PicEmbeddedBC* embeddedBC = nullptr,
                      const PicEmbeddedBCMap* embeddedBCById = nullptr) const;

    bool solveWithCache(const QuadtreeMesh& mesh,
                        const FaceSegmentGraph& graph,
                        const PicDomainBC& bc,
                        const CellField<double>& rhs,
                        CellField<double>& phi,
                        SolveReport* report = nullptr,
                        const Options& opt = Options(),
                        const Solvers::Poisson::CutCellMetrics* cutMetrics = nullptr,
                        const PicEmbeddedBC* embeddedBC = nullptr,
                        const PicEmbeddedBCMap* embeddedBCById = nullptr) const;

    bool hasCache() const;
    void clearCache() const;

    // Unit-box manufactured-solution regression.
    ErrorNorms regressionUnitBoxSinSin(const QuadtreeMesh& mesh,
                                       const FaceSegmentGraph& graph,
                                       CellField<double>& phiOut,
                                       const Options& opt = Options()) const;

    // Unit-box manufactured-solution regression with y-periodic wrapping (bottom/top Periodic).
    ErrorNorms regressionUnitBoxSinCosPeriodicY(const QuadtreeMesh& mesh,
                                               const FaceSegmentGraph& graph,
                                               CellField<double>& phiOut,
                                               const Options& opt = Options()) const;

private:
    struct FactorCache;
    mutable std::unique_ptr<FactorCache> cache_;

    static bool cacheMatches(const FactorCache& c,
                             const QuadtreeMesh& mesh,
                             const FaceSegmentGraph& graph,
                             const PicDomainBC& bc,
                             const Options& opt,
                             const Solvers::Poisson::CutCellMetrics* cutMetrics,
                             const PicEmbeddedBC* embeddedBC,
                             const PicEmbeddedBCMap* embeddedBCById);

    // Cache statistics / control. These are best-effort diagnostics.
    mutable long long cacheSolveCalls_ = 0;
    mutable long long cacheHits_ = 0;
    mutable long long cacheMisses_ = 0;
};
