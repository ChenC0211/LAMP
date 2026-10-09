#pragma once

#include <string>
#include <vector>
#include <cstdint>
#include <memory>

#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

#ifdef EIGEN_SUPERLU_SUPPORT
#include <Eigen/SuperLUSupport>
#endif

#include "Core/Field/CellField.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Solvers/Poisson/PoissonBC.h"

namespace Solvers { namespace Poisson { class CutCellMetrics; } }

class QuadtreeMesh;

// ============================================================================
// Cell-centered FV variable-coefficient Poisson solver (polarization closure)
//
// Governing form (dimensionless):
//   -div( kappa * grad(phi) ) = rhs
// where
//   kappa = 1 + mu/epsilon0  (mu in [F/m])
//
// Discretization mirrors PoissonFVSolverCellCentered (same mesh/segment graph,
// same cut-cell metrics + EB closure), but scales face conductance by kappa_f.
//
// Inputs:
//   - kappaBySegment: per face-segment coefficient (size==graph.segments().size())
//   - kappaByCellIndex: optional per-cell coefficient used for EB closure terms
//     (size==mesh.cellIndexNum()). If null, EB terms assume kappa=1.
// ============================================================================
class PolarizationPoissonFVSolverCellCentered
{
public:
    struct FluxCoefficient
    {
        int cellIndex = -1; // mesh cell INDEX, never a cell label
        double value = 0.0;
    };

    // Additional integrated outward flux of this oriented graph segment:
    // q_add = sum(coefficients[j].value * phi[cellIndex]) + constant.
    // Its owner row receives +coefficients in A and -constant in b.
    // In axisymmetric mode this must ALREADY include the physical 2*pi*r
    // surface measure; the solver does not multiply corrections by radius.
    struct SegmentFluxCorrection
    {
        std::vector<FluxCoefficient> coefficients;
        double constant = 0.0;
    };

    enum class FittedBoundaryType { Interior, Dirichlet, HomogeneousNeumann };
    struct FittedBoundary
    {
        FittedBoundaryType type = FittedBoundaryType::Interior;
        double value = 0.0; // static potential; exactly zero for Interior/Neumann
    };

    struct Options
    {
        bool useSuperLU = true;
        // Strict caller contract: require SuperLU at build/runtime and never
        // fall back to SparseLU or CG. Defaults preserve existing PIC callers.
        bool requireSuperLU = false;

        // Iterative-solver knobs are kept for future CG-based solvers.
        bool useCG = true;
        bool useIC = true;
        int maxIter = 2000;
        double relTol = 1e-10;
        bool reuseIterativePattern = true;
        bool reuseDirectPattern = true;
        // Reuse direct LU factors only when the complete assembled matrix
        // (indices and numeric bytes) is unchanged; otherwise refactorize.
        bool freezeDirectMatrixUntilPatternChange = false;
        // Explicit fitted-boundary port only. Reuse SuperLU factors only when
        // the complete assembled matrix (indices AND numeric bytes) and native
        // row identities are unchanged. All input/full-residual gates still run.
        // This is independent of legacy pattern/freezing options and environment.
        bool reuseFittedFactorization = false;
        bool useTemporalExtrapolatedGuess = false;
        double temporalExtrapolationOmega = 1.0;

        bool verbose = false;

        enum class EmbeddedBCMissingPolicy
        {
            Error,
            Skip
        };
        EmbeddedBCMissingPolicy embeddedBCMissingPolicy = EmbeddedBCMissingPolicy::Error;

        // Strict fluxCorrections or fitted-boundary port: x=r, y=z, r>=0. TPFA face
        // measure is 2*pi*r_face*Af; cell RHS measure is 2*pi*r_cell*area.
        // The legacy nullptr port rejects true and retains its planar method.
        bool axisymmetric = false;
    };

    struct SolveReport
    {
        bool success = false;
        int iterations = 0;
        // With fluxCorrections or fitted boundaries: ||A*x-b||_2 / ||b||_2
        // (absolute if b is zero),
        // evaluated on the complete corrected operator and gated by relTol.
        // Legacy nullptr callers retain the previous report semantics.
        double finalResidual = 0.0;
        std::string message;

        double timeAssemble_s = 0.0;
        double timeFactor_s = 0.0;
        double timeSolve_s = 0.0;
        double timeTotal_s = 0.0;

        int N = 0;
        long long nnz = 0;
        bool usedSuperLU = false;
        bool usedCG = false;
        bool cgAttempted = false;
        bool cgFallbackToDirect = false;
        int cgIterations = 0;
        double cgResidual = 0.0;
        bool cgPatternReused = false;
        bool cgTemporalGuessUsed = false;
        bool directPatternReused = false;
        bool directMatrixFrozen = false;
        bool fittedFactorizationReused = false;
        bool fittedRhsOnlySolve = false; // Explicit owned-operator RHS entry, not a full reassembly.
    };

    PolarizationPoissonFVSolverCellCentered();
    ~PolarizationPoissonFVSolverCellCentered();

    PolarizationPoissonFVSolverCellCentered(const PolarizationPoissonFVSolverCellCentered&) = delete;
    PolarizationPoissonFVSolverCellCentered& operator=(const PolarizationPoissonFVSolverCellCentered&) = delete;

    // Optional fluxCorrections has one entry per oriented graph segment.
    // This port supports complete uniform Cartesian meshes, planar by default
    // or axisymmetric with Options::axisymmetric=true (x=r>=0, y=z), and static
    // finite Dirichlet or exactly homogeneous Neumann domain BCs, with at least
    // one positive-measure Dirichlet side. At r=0 the face measure is exactly
    // zero and only homogeneous Neumann with zero corrections is accepted.
    // Nonzero Neumann, all-Neumann, time-dependent, Robin,
    // and periodic BCs are rejected. Every Neumann segment's correction constant
    // and coefficient values must each equal zero (empty or explicit signed-zero
    // coefficients are legal, still subject to the index/finite checks below).
    // No phase/interface assumption is made by this generic solver.
    // No cut metrics/embedded BCs; positive explicitly supplied segment kappa
    // and strict SuperLU options are required. Coefficients
    // must be sorted/unique active cell indices; reverse interior entries must
    // be exact negatives with identical geometry/kappa. Zero entries are legal.
    // The complete matrix must retain a positive finite diagonal; no diagonal
    // repair, compatibility adjustment, or gauge pin is applied by this port.
    // It disables all pattern/factor reuse, validates before changing phi, and
    // commits phi only after successful finite full-operator residual checking.
    // A default-constructed phi is supported. The nullptr path is unchanged.
    bool solve(const QuadtreeMesh& mesh,
               const FaceSegmentGraph& graph,
               const PicDomainBC& bc,
               const CellField<double>& rhs,
               CellField<double>& phi,
               const std::vector<double>* kappaBySegment,
               const std::vector<double>* kappaByCellIndex,
               SolveReport* report = nullptr,
               const Options& opt = Options(),
               const Solvers::Poisson::CutCellMetrics* cutMetrics = nullptr,
               const PicEmbeddedBC* embeddedBC = nullptr,
               const PicEmbeddedBCMap* embeddedBCById = nullptr,
               const std::vector<SegmentFluxCorrection>* fluxCorrections = nullptr) const;

    // Explicit active-leaf TPFA operator on the actual quadtree graph (including
    // refined identities and split faces). Arrays follow oriented segment order.
    // Interior means an active neighbor with an identical reverse kappa. Every
    // domain/inactive-neighbor face must instead specify static Dirichlet or
    // exactly homogeneous Neumann. Dirichlet distance is owner center to the
    // actual face, not to the inactive cell center. No cut-cell/solid fraction
    // is inferred: the caller is responsible for a genuinely face-fitted mask.
    // Uses full cell/face measures, revolved by 2*pi*r in axisymmetric mode;
    // zero-area axis faces must be homogeneous Neumann. Every active connected
    // component needs a positive-area Dirichlet anchor. Strict SuperLU options
    // and positive finite per-segment kappa are mandatory; no gauge/repair or
    // iterative fallback. Complete residual uses the unchanged relTol gate.
    // phi must already be mesh-initialized. Success preserves inactive values,
    // field labels and storage; failure leaves all phi values unchanged.
    // Optional reuseFittedFactorization retains exact-matrix factors on this
    // solver instance, not a global cache. Changes refactorize; any failure or
    // clearDirectCache invalidates them. Concurrent calls on one instance are
    // not supported. Default false preserves one-shot factorization.
    bool solveFittedBoundary(const QuadtreeMesh& mesh, const FaceSegmentGraph& graph,
               const CellField<double>& rhs, CellField<double>& phi,
               const std::vector<double>& kappaBySegment,
               const std::vector<FittedBoundary>& boundaryBySegment,
               const Options& opt, SolveReport* report = nullptr) const;

    // The SAME fitted operator with explicit additional integrated outward
    // flux on every oriented segment. This is not a different discretization:
    // Dirichlet TPFA still uses owner-center to actual-face distance, while
    // correction coefficients enter A and their constants enter -b exactly as
    // documented above. In RZ corrections already include physical 2*pi*r.
    // All fitted-boundary contracts remain in force. Correction coefficients
    // must have sorted, unique ACTIVE mesh cell indices (not labels or packed
    // solver rows), with finite values/constants. Interior reverse entries
    // must have identical indices and exactly opposite values/constants.
    // Every homogeneous-Neumann/zero-area-axis correction must be identically
    // zero, including EACH coefficient; signed zeros remain legal but their
    // cell indices must still be active. Dirichlet corrections may be nonzero.
    // No phase, surface-charge law or interface-to-face mapping is inferred.
    // Uses the same strict SuperLU, complete residual and no-repair gates;
    // all validation and report construction precede the final value commit.
    // Inactive values, field labels and storage are preserved on success, and
    // every phi value is unchanged on failure. The zero-correction overload
    // above and the legacy solve entry keep their existing behavior.
    bool solveFittedBoundary(const QuadtreeMesh& mesh, const FaceSegmentGraph& graph,
               const CellField<double>& rhs, CellField<double>& phi,
               const std::vector<double>& kappaBySegment,
               const std::vector<FittedBoundary>& boundaryBySegment,
               const std::vector<SegmentFluxCorrection>& fluxCorrections,
               const Options& opt, SolveReport* report = nullptr) const;

    // Solve the EXACT owned operator from the latest successful fitted solve
    // with reuseFittedFactorization=true. Only volumetric RHS and correction
    // constants may change. No mesh/graph/material/boundary inputs are borrowed
    // or reinterpreted here: this explicitly solves the prepared snapshot, NOT
    // a possibly changed live mesh. Reprepare via the full fitted entry when
    // geometry, coefficients, static boundary values or native identity change.
    // constants must be null for the TPFA overload, or cover every original
    // oriented segment for the corrected overload. Finite/zero-Neumann/reverse
    // constant, native field identity, strict SuperLU and full residual gates
    // still run. Preserves all phi data/storage on failure and inactive payload
    // on success. Any failure or clearDirectCache destroys the prepared state.
    // No concurrent calls on one solver; no pattern-freezing mode is involved.
    bool solvePreparedFittedRhs(const CellField<double>& rhs, CellField<double>& phi,
               const std::vector<double>* constants, double relTol,
               SolveReport* report = nullptr) const;

    // Clear iterative/direct internal solver caches (e.g. after AMR topology change).
    void clearDirectCache() const;

private:
    bool solveImpl(const QuadtreeMesh& mesh, const FaceSegmentGraph& graph,
                   const PicDomainBC& bc, const CellField<double>& rhs,
                   CellField<double>& phi, const std::vector<double>* kappaBySegment,
                   const std::vector<double>* kappaByCellIndex, SolveReport* report,
                   const Options& opt, const Solvers::Poisson::CutCellMetrics* cutMetrics,
                   const PicEmbeddedBC* embeddedBC, const PicEmbeddedBCMap* embeddedBCById,
                   const std::vector<SegmentFluxCorrection>* fluxCorrections,
                   const std::vector<FittedBoundary>* fittedBoundaries = nullptr) const;

    using CgNoPrecond = Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper>;
    using CgWithIC = Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper,
                                              Eigen::IncompleteCholesky<double>>;

    struct IterativeCache
    {
        int icN = -1;
        long long icNnz = -1;
        std::uint64_t icPatternSig = 0;
        bool icAnalyzed = false;

        int cgN = -1;
        long long cgNnz = -1;
        std::uint64_t cgPatternSig = 0;
        bool cgAnalyzed = false;

        CgWithIC cgIC;
        CgNoPrecond cg;

        int xN = -1;
        std::uint64_t xPatternSig = 0;
        bool hasLastX = false;
        bool hasPrevX = false;
        Eigen::VectorXd lastX;
        Eigen::VectorXd prevX;
    };

    struct DirectCache
    {
#ifdef EIGEN_SUPERLU_SUPPORT
        Eigen::SparseMatrix<double> superMatrix;
        bool superAnalyzed = false;
        bool superFactorized = false;
        Eigen::SuperLU<Eigen::SparseMatrix<double>> superLU;
#endif

        Eigen::SparseMatrix<double> sparseMatrix;
        bool sparseAnalyzed = false;
        bool sparseFactorized = false;
        Eigen::SparseLU<Eigen::SparseMatrix<double>> sparseLU;
    };

    mutable IterativeCache iterCache_{};
    mutable DirectCache directCache_{};
    struct FittedFactorCache;
    mutable std::unique_ptr<FittedFactorCache> fittedFactorCache_;
};
