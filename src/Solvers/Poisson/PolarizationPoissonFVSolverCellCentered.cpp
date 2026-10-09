#include "Solvers/Poisson/PolarizationPoissonFVSolverCellCentered.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <new>
#include <queue>
#include <sstream>

#include <Eigen/SparseLU>
#include <Eigen/IterativeLinearSolvers>

#ifdef EIGEN_SUPERLU_SUPPORT
#include <Eigen/SuperLUSupport>
#endif

#include "Core/IO/ConsoleControl.h"
#include "Core/Mesh/PeriodicFaceUtils.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/CartesianMeasure.h"
#include "Solvers/Poisson/CutCellMetrics.h"
#include "Solvers/Poisson/detail/PolarizationFluxCorrectionValidation.h"
#include "Solvers/Poisson/detail/PolarizationFittedBoundaryValidation.h"
#include "Solvers/Poisson/detail/PolarizationFittedRhsPlan.h"

// Like the fitted Hydro cache, this compares actual compressed coefficients,
// not a pattern hash. It never reuses the legacy matrix-freezing path.
struct PolarizationPoissonFVSolverCellCentered::FittedFactorCache
{
    Eigen::SparseMatrix<double> matrix;
    const QuadtreeMesh* mesh = nullptr; // identity only; never dereferenced
    std::vector<int> cellIndices, cellLabels;
    std::unique_ptr<const PolarizationFittedBoundaryDetail::RhsPlan> rhsPlan;
#ifdef EIGEN_SUPERLU_SUPPORT
    Eigen::SuperLU<Eigen::SparseMatrix<double>> solver;
#endif

    bool matches(const Eigen::SparseMatrix<double>& candidate,
        const QuadtreeMesh& candidateMesh, const std::vector<int>& rows) const
    {
        if (mesh != &candidateMesh || cellIndices != rows ||
            matrix.rows() != candidate.rows() || matrix.cols() != candidate.cols() ||
            !matrix.isCompressed() || !candidate.isCompressed() ||
            matrix.nonZeros() != candidate.nonZeros()) return false;
        for (std::size_t i = 0; i < rows.size(); ++i)
            if (cellLabels[i] != candidateMesh.cellIndexToLabel[rows[i]]) return false;
        using Index = Eigen::SparseMatrix<double>::StorageIndex;
        return std::memcmp(matrix.outerIndexPtr(), candidate.outerIndexPtr(),
                   (candidate.outerSize() + 1) * sizeof(Index)) == 0 &&
               std::memcmp(matrix.innerIndexPtr(), candidate.innerIndexPtr(),
                   candidate.nonZeros() * sizeof(Index)) == 0 &&
               std::memcmp(matrix.valuePtr(), candidate.valuePtr(),
                   candidate.nonZeros() * sizeof(double)) == 0;
    }
};

PolarizationPoissonFVSolverCellCentered::PolarizationPoissonFVSolverCellCentered() = default;
PolarizationPoissonFVSolverCellCentered::~PolarizationPoissonFVSolverCellCentered() = default;

namespace
{
using Clock = std::chrono::steady_clock;

auto secondsBetween(const Clock::time_point& a, const Clock::time_point& b) -> double
{
    return std::chrono::duration_cast<std::chrono::duration<double>>(b - a).count();
}

inline bool isDirichletLike(const PicBCType t)
{
    return (t == PicBCType::Dirichlet) || (t == PicBCType::Robin);
}

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

inline const char* boundaryToString(const Boundary b)
{
    switch (b)
    {
    case Boundary::west:  return "west";
    case Boundary::south: return "south";
    case Boundary::east:  return "east";
    case Boundary::north: return "north";
    default:              return "non";
    }
}

inline double safeCellArea(const QuadtreeMesh& mesh, int cellIndex)
{
    if (cellIndex < 0 || cellIndex >= mesh.cellIndexNum())
        return 0.0;
    const int cLbl = mesh.cellIndexToLabel[cellIndex];
    const Cell& c = mesh.cells(cLbl);
    const double a = c.area();
    if (!std::isfinite(a) || !(a > 0.0))
        return 0.0;
    return a;
}

inline double stabilizedEmbeddedDistance(const double dGamma, const double cellArea)
{
    const double h = (std::isfinite(cellArea) && cellArea > 0.0) ? std::sqrt(cellArea) : 1.0;
    const double floor = std::max(1e-12, 2.0e-2 * h);
    if (!std::isfinite(dGamma))
        return floor;
    return std::max(dGamma, floor);
}

inline double clampKappa(double v)
{
    if (!std::isfinite(v))
        return 1.0;
    if (v < 0.0)
        return 0.0;
    return v;
}

inline std::uint64_t sparsePatternSignature(const Eigen::SparseMatrix<double>& A)
{
    constexpr std::uint64_t kOffset = 1469598103934665603ull;
    constexpr std::uint64_t kPrime = 1099511628211ull;
    auto mix = [kPrime](std::uint64_t h, std::uint64_t v) -> std::uint64_t {
        h ^= v;
        h *= kPrime;
        return h;
    };

    std::uint64_t h = kOffset;
    h = mix(h, static_cast<std::uint64_t>(A.rows()));
    h = mix(h, static_cast<std::uint64_t>(A.cols()));
    h = mix(h, static_cast<std::uint64_t>(A.nonZeros()));

    const Eigen::Index outerSize = A.outerSize();
    if (outerSize > 0)
    {
        const int* outer = A.outerIndexPtr();
        const Eigen::Index stride = std::max<Eigen::Index>(1, outerSize / 64);
        for (Eigen::Index i = 0; i < outerSize; i += stride)
            h = mix(h, static_cast<std::uint64_t>(outer[i]));
        h = mix(h, static_cast<std::uint64_t>(outer[outerSize]));
    }

    const Eigen::Index nnz = A.nonZeros();
    if (nnz > 0)
    {
        const int* inner = A.innerIndexPtr();
        const Eigen::Index stride = std::max<Eigen::Index>(1, nnz / 256);
        for (Eigen::Index k = 0; k < nnz; k += stride)
            h = mix(h, static_cast<std::uint64_t>(inner[k]));
        h = mix(h, static_cast<std::uint64_t>(inner[nnz - 1]));
    }

    return h;
}

inline bool sameSparsePattern(const Eigen::SparseMatrix<double>& a,
                              const Eigen::SparseMatrix<double>& b)
{
    if (a.rows() != b.rows() || a.cols() != b.cols() ||
        !a.isCompressed() || !b.isCompressed() || a.nonZeros() != b.nonZeros())
        return false;
    using Index = Eigen::SparseMatrix<double>::StorageIndex;
    return std::memcmp(a.outerIndexPtr(), b.outerIndexPtr(),
               (a.outerSize() + 1) * sizeof(Index)) == 0 &&
           std::memcmp(a.innerIndexPtr(), b.innerIndexPtr(),
               a.nonZeros() * sizeof(Index)) == 0;
}

inline bool sameSparseValues(const Eigen::SparseMatrix<double>& a,
                             const Eigen::SparseMatrix<double>& b)
{
    // Called only after exact pattern equality has been established.
    return std::memcmp(a.valuePtr(), b.valuePtr(), a.nonZeros() * sizeof(double)) == 0;
}

inline bool envBool(const char* name, bool fallback)
{
    if (!name)
        return fallback;
    std::string v;
#if defined(_MSC_VER)
    char* s = nullptr;
    size_t len = 0;
    if (_dupenv_s(&s, &len, name) != 0 || !s)
        return fallback;
    v.assign(s);
    std::free(s);
#else
    const char* s = std::getenv(name);
    if (!s)
        return fallback;
    v.assign(s);
#endif
    std::transform(v.begin(), v.end(), v.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (v == "1" || v == "true" || v == "yes" || v == "on")
        return true;
    if (v == "0" || v == "false" || v == "no" || v == "off")
        return false;
    return fallback;
}

// Shared by the full fitted/corrected solve and the owned-operator RHS solve.
// Keep the original complete residual expression and tolerance semantics.
bool fullResidualAccepted(const Eigen::SparseMatrix<double>& matrix, const Eigen::VectorXd& x,
    const Eigen::VectorXd& b, bool usedSuperLU, double relTol,
    double& relativeResidual, std::string& reason)
{
    const Eigen::VectorXd residual = matrix * x - b;
    const double rhsNorm = b.norm();
    const double residualNorm = residual.norm();
    if (!usedSuperLU || !x.allFinite() || !residual.allFinite() ||
        !std::isfinite(rhsNorm) || !std::isfinite(residualNorm))
    {
        reason = "fluxCorrections: non-finite solution/full-operator residual norm or required SuperLU was not used";
        return false;
    }
    relativeResidual = rhsNorm > 0.0 ? residualNorm / rhsNorm : residualNorm;
    if (!std::isfinite(relativeResidual) || relativeResidual > relTol)
    {
        std::ostringstream message;
        message << "fluxCorrections: full-operator residual " << relativeResidual
                << " exceeds relTol=" << relTol;
        reason = message.str();
        return false;
    }
    return true;
}

// Called ONLY after the fitted validator has checked the complete actual
// active-leaf graph, fields, coefficients, boundary types and reverse geometry.
// The legacy flux validator deliberately remains uniform/full-domain-only;
// this supplement adds its algebraic contract without substituting its mesh.
void validateFittedFluxCorrections(const QuadtreeMesh& mesh, const FaceSegmentGraph& graph,
    const std::vector<PolarizationPoissonFVSolverCellCentered::FittedBoundary>& boundaries,
    const std::vector<PolarizationPoissonFVSolverCellCentered::SegmentFluxCorrection>& corrections)
{
    using Type = PolarizationPoissonFVSolverCellCentered::FittedBoundaryType;
    const auto require = [](bool condition, const char* message)
    {
        if (!condition) throw std::runtime_error(std::string("fitted fluxCorrections: ") + message);
    };
    const auto& segments = graph.segments();
    const auto& offsets = graph.cellOffsets();
    const int count = mesh.cellIndexNum();
    require(corrections.size() == segments.size(), "correction count must equal oriented segment count");
    for (std::size_t sid = 0; sid < segments.size(); ++sid)
    {
        const auto& correction = corrections[sid];
        require(std::isfinite(correction.constant), "non-finite correction constant");
        int previous = -1;
        for (const auto& coefficient : correction.coefficients)
        {
            const int index = coefficient.cellIndex;
            require(index > previous && index < count && std::isfinite(coefficient.value),
                "coefficients must have sorted unique active cell indices and finite values");
            require(offsets[index] < offsets[index + 1], "correction coefficient references an inactive cell index");
            previous = index;
        }
        if (boundaries[sid].type == Type::HomogeneousNeumann)
        {
            // The fitted geometry validator requires this same type on the
            // zero-area axis. Neither cancellation nor a small-value tolerance
            // can turn a prescribed zero-flux boundary into another operator.
            require(correction.constant == 0.0,
                "homogeneous Neumann boundary correction constant must be exactly zero");
            for (const auto& coefficient : correction.coefficients)
                require(coefficient.value == 0.0,
                    "homogeneous Neumann boundary correction coefficients must each be exactly zero");
        }
    }
    for (std::size_t sid = 0; sid < segments.size(); ++sid)
    {
        if (boundaries[sid].type != Type::Interior) continue;
        const auto& face = segments[sid];
        const int neighbor = mesh.cells(face.L).index(); // Already validated as active.
        std::size_t reverseId = segments.size();
        for (std::size_t candidate = offsets[neighbor]; candidate < offsets[neighbor + 1]; ++candidate)
        {
            const auto& reverse = segments[candidate];
            if (reverse.isHorizontal != face.isHorizontal || reverse.edgeLabel != face.edgeLabel) continue;
            require(reverseId == segments.size() && reverse.K == face.L && reverse.L == face.K,
                "interior correction requires one reverse incidence of the same actual face");
            reverseId = candidate;
        }
        require(reverseId < segments.size(), "interior correction has no reverse actual face");
        if (sid > reverseId) continue; // Check every physical pair once, without coalescing split faces.
        const auto& forward = corrections[sid];
        const auto& reverse = corrections[reverseId];
        require(reverse.constant == -forward.constant && reverse.coefficients.size() == forward.coefficients.size(),
            "interior correction pair must be exactly antisymmetric");
        for (std::size_t j = 0; j < forward.coefficients.size(); ++j)
            require(reverse.coefficients[j].cellIndex == forward.coefficients[j].cellIndex
                && reverse.coefficients[j].value == -forward.coefficients[j].value,
                "interior correction coefficients must be exactly antisymmetric");
    }
}

} // namespace

bool PolarizationPoissonFVSolverCellCentered::solve(const QuadtreeMesh& mesh,
                                                   const FaceSegmentGraph& graph,
                                                   const PicDomainBC& bc,
                                                   const CellField<double>& rhs,
                                                   CellField<double>& phi,
                                                   const std::vector<double>* kappaBySegment,
                                                   const std::vector<double>* kappaByCellIndex,
                                                   SolveReport* report,
                                                   const Options& opt,
                                                   const Solvers::Poisson::CutCellMetrics* cutMetrics,
                                                   const PicEmbeddedBC* embeddedBC,
                                                   const PicEmbeddedBCMap* embeddedBCById,
                                                   const std::vector<SegmentFluxCorrection>* fluxCorrections) const
{
    // This option has a consumer only on the explicit fitted-boundary port.
    fittedFactorCache_.reset();
    if (opt.reuseFittedFactorization)
    {
        if (report)
        {
            *report = {};
            report->message = "reuseFittedFactorization requires solveFittedBoundary";
        }
        return false;
    }
    if (opt.axisymmetric && !fluxCorrections)
    {
        if (report)
        {
            *report = {};
            report->message = "axisymmetric mode requires the strict integrated fluxCorrections port";
        }
        return false;
    }
    if (!fluxCorrections)
        return solveImpl(mesh, graph, bc, rhs, phi, kappaBySegment, kappaByCellIndex,
                         report, opt, cutMetrics, embeddedBC, embeddedBCById, nullptr);

    const auto begin = Clock::now();
    // The corrected operator has no cache identity contract. Never consult old
    // pattern/value hashes (including environment-enabled freeze modes).
    clearDirectCache();
    try
    {
        PolarizationFluxCorrectionValidation::validate(mesh, graph, bc, rhs, phi,
            kappaBySegment, kappaByCellIndex, opt, cutMetrics, embeddedBC,
            embeddedBCById, *fluxCorrections);
        Options strictOptions = opt;
        strictOptions.reuseDirectPattern = false;
        strictOptions.freezeDirectMatrixUntilPatternChange = false;
        strictOptions.reuseIterativePattern = false;
        strictOptions.useTemporalExtrapolatedGuess = false;
        const bool success = solveImpl(mesh, graph, bc, rhs, phi, kappaBySegment,
            kappaByCellIndex, report, strictOptions, nullptr, nullptr, nullptr, fluxCorrections);
        if (!success) clearDirectCache();
        return success;
    }
    catch (const std::exception& error)
    {
        clearDirectCache();
        if (report)
        {
            *report = {};
            report->message = std::string("Integrated flux-correction solve failed: ") + error.what();
            report->timeTotal_s = secondsBetween(begin, Clock::now());
        }
        return false;
    }
}

bool PolarizationPoissonFVSolverCellCentered::solveFittedBoundary(
    const QuadtreeMesh& mesh, const FaceSegmentGraph& graph,
    const CellField<double>& rhs, CellField<double>& phi,
    const std::vector<double>& kappaBySegment,
    const std::vector<FittedBoundary>& boundaryBySegment,
    const Options& opt, SolveReport* report) const
{
    const auto begin = Clock::now();
    if (!opt.reuseFittedFactorization) clearDirectCache();
    try
    {
        PolarizationFittedBoundaryValidation::validate(mesh, graph, rhs, phi,
            kappaBySegment, boundaryBySegment, opt);
        // The shared solver publishes into this owned stage. Inactive values
        // are retained, and final publication needs no allocation or relabeling.
        CellField<double> nextPhi(phi);
        SolveReport nextReport;
        Options strictOptions = opt;
        strictOptions.reuseDirectPattern = false;
        strictOptions.freezeDirectMatrixUntilPatternChange = false;
        strictOptions.reuseIterativePattern = false;
        strictOptions.useTemporalExtrapolatedGuess = false;
        const PicDomainBC unusedDomainBC{};
        const bool success = solveImpl(mesh, graph, unusedDomainBC, rhs, nextPhi,
            &kappaBySegment, nullptr, &nextReport, strictOptions, nullptr, nullptr,
            nullptr, nullptr, &boundaryBySegment);
        if (!success) clearDirectCache();
        // Report construction/copy also precedes the no-throw value commit.
        if (report) *report = std::move(nextReport);
        if (!success) return false;
        std::copy_n(nextPhi.value.get(), phi.size(), phi.value.get());
        return true;
    }
    catch (const std::exception& error)
    {
        clearDirectCache();
        if (report)
        {
            *report = {};
            report->message = std::string("Fitted-boundary solve failed: ") + error.what();
            report->timeTotal_s = secondsBetween(begin, Clock::now());
        }
        return false;
    }
}

bool PolarizationPoissonFVSolverCellCentered::solveFittedBoundary(
    const QuadtreeMesh& mesh, const FaceSegmentGraph& graph,
    const CellField<double>& rhs, CellField<double>& phi,
    const std::vector<double>& kappaBySegment,
    const std::vector<FittedBoundary>& boundaryBySegment,
    const std::vector<SegmentFluxCorrection>& fluxCorrections,
    const Options& opt, SolveReport* report) const
{
    const auto begin = Clock::now();
    if (!opt.reuseFittedFactorization) clearDirectCache();
    try
    {
        PolarizationFittedBoundaryValidation::validate(mesh, graph, rhs, phi,
            kappaBySegment, boundaryBySegment, opt);
        validateFittedFluxCorrections(mesh, graph, boundaryBySegment, fluxCorrections);
        // Keep the old fitted entry untouched. Both explicit inputs use the
        // existing combined assembly and publish only this owned staged field.
        CellField<double> nextPhi(phi);
        SolveReport nextReport;
        Options strictOptions = opt;
        strictOptions.reuseDirectPattern = false;
        strictOptions.freezeDirectMatrixUntilPatternChange = false;
        strictOptions.reuseIterativePattern = false;
        strictOptions.useTemporalExtrapolatedGuess = false;
        const PicDomainBC unusedDomainBC{};
        const bool success = solveImpl(mesh, graph, unusedDomainBC, rhs, nextPhi,
            &kappaBySegment, nullptr, &nextReport, strictOptions, nullptr, nullptr,
            nullptr, &fluxCorrections, &boundaryBySegment);
        if (!success) clearDirectCache();
        // No observable phi write precedes report construction, complete
        // corrected-operator residual verification, or a possible exception.
        if (report) *report = std::move(nextReport);
        if (!success) return false;
        std::copy_n(nextPhi.value.get(), phi.size(), phi.value.get());
        return true;
    }
    catch (const std::exception& error)
    {
        clearDirectCache();
        if (report)
        {
            *report = {};
            report->message = std::string("Fitted-boundary integrated flux-correction solve failed: ") + error.what();
            report->timeTotal_s = secondsBetween(begin, Clock::now());
        }
        return false;
    }
}

bool PolarizationPoissonFVSolverCellCentered::solvePreparedFittedRhs(
    const CellField<double>& rhs, CellField<double>& phi,
    const std::vector<double>* constants, double relTol, SolveReport* report) const
{
    const auto begin = Clock::now();
    try
    {
        if (!fittedFactorCache_ || !fittedFactorCache_->rhsPlan)
            throw std::runtime_error("no successfully prepared fitted operator");
        if (!std::isfinite(relTol) || relTol <= 0.0 || relTol >= 1.0)
            throw std::runtime_error("relTol must be finite and in (0,1)");
        auto& prepared = *fittedFactorCache_;
        const auto b = prepared.rhsPlan->assemble(rhs, phi, constants);
        SolveReport next;
        next.N = static_cast<int>(prepared.matrix.rows());
        next.nnz = static_cast<long long>(prepared.matrix.nonZeros());
        next.fittedFactorizationReused = next.fittedRhsOnlySolve = true;
        next.timeAssemble_s = secondsBetween(begin, Clock::now());
        Eigen::VectorXd x;
#ifdef EIGEN_SUPERLU_SUPPORT
        const auto solveBegin = Clock::now();
        x = prepared.solver.solve(b);
        next.timeSolve_s = secondsBetween(solveBegin, Clock::now());
        next.usedSuperLU = prepared.solver.info() == Eigen::Success;
#endif
        std::string reason = "prepared fitted SuperLU backsolve failed; no fallback is allowed";
        next.success = next.usedSuperLU && fullResidualAccepted(prepared.matrix, x, b,
            next.usedSuperLU, relTol, next.finalResidual, reason);
        next.timeTotal_s = secondsBetween(begin, Clock::now());
        next.message = next.success
            ? "Solved owned face-fitted operator with updated RHS; full-operator residual verified."
            : reason;
        const bool accepted = next.success;
        // Report construction and all fallible work precede phi publication.
        if (report) *report = std::move(next);
        if (!accepted)
        {
            clearDirectCache();
            return false;
        }
        for (std::size_t row = 0; row < prepared.cellIndices.size(); ++row)
            phi.value[prepared.cellIndices[row]] = x[row];
        return true;
    }
    catch (const std::exception& error)
    {
        clearDirectCache();
        if (report)
        {
            *report = {};
            report->message = std::string("Prepared fitted RHS solve failed: ") + error.what();
            report->timeTotal_s = secondsBetween(begin, Clock::now());
        }
        return false;
    }
    catch (...)
    {
        clearDirectCache();
        throw;
    }
}

bool PolarizationPoissonFVSolverCellCentered::solveImpl(const QuadtreeMesh& mesh,
                                                   const FaceSegmentGraph& graph,
                                                   const PicDomainBC& bc,
                                                   const CellField<double>& rhs,
                                                   CellField<double>& phi,
                                                   const std::vector<double>* kappaBySegment,
                                                   const std::vector<double>* kappaByCellIndex,
                                                   SolveReport* report,
                                                   const Options& opt,
                                                   const Solvers::Poisson::CutCellMetrics* cutMetrics,
                                                   const PicEmbeddedBC* embeddedBC,
                                                   const PicEmbeddedBCMap* embeddedBCById,
                                                   const std::vector<SegmentFluxCorrection>* fluxCorrections,
                                                   const std::vector<FittedBoundary>* fittedBoundaries) const
{
    const auto tBegin = Clock::now();
    const bool correctedOperator = fluxCorrections != nullptr;
    const bool fittedOperator = fittedBoundaries != nullptr;
    const bool checkedOperator = correctedOperator || fittedOperator;

    // Validate the strict backend contract before changing phi or solver state.
    auto requiredBackendFailure = [&](const std::string& message) {
        if (report)
        {
            *report = {};
            report->message = message;
        }
        return false;
    };
    if (opt.requireSuperLU && (!opt.useSuperLU || opt.useCG || opt.useIC))
        return requiredBackendFailure("requireSuperLU requires useSuperLU=true, useCG=false and useIC=false; no alternative backend is allowed");
#ifndef EIGEN_SUPERLU_SUPPORT
    if (opt.requireSuperLU)
        return requiredBackendFailure("requireSuperLU requested but SuperLU support is unavailable in this build; SparseLU fallback is disabled");
#endif

    auto secondsSince = [&](const Clock::time_point& t0) -> double {
        return secondsBetween(t0, Clock::now());
    };

    if (!checkedOperator && (phi.size() != mesh.cellIndexNum() || !phi.value))
        phi.initialize(mesh, 0.0);

    // Active unknown set: cell indices with at least one face segment.
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
            *report = {};
            report->success = false;
            report->message = "No active cells in FaceSegmentGraph (N=0).";
            report->timeTotal_s = secondsSince(tBegin);
            report->timeAssemble_s = report->timeTotal_s;
        }
        return false;
    }

    const auto& segs = graph.segments();
    const bool useKappaBySegment =
        (kappaBySegment != nullptr) && (kappaBySegment->size() == segs.size());
    const bool useKappaByCell =
        (kappaByCellIndex != nullptr) && (kappaByCellIndex->size() == static_cast<size_t>(mesh.cellIndexNum()));

    const bool useCutMetrics =
        (cutMetrics != nullptr) &&
        (cutMetrics->betaBySegment().size() == segs.size()) &&
        (cutMetrics->cellMetricsByIndex().size() == static_cast<size_t>(mesh.cellIndexNum()));

    const bool useEmbeddedBCMap =
        useCutMetrics && (embeddedBCById != nullptr) && (!embeddedBCById->empty());
    const bool useEmbeddedBCFallback =
        useCutMetrics && (embeddedBC != nullptr) && embeddedBC->enabled;
    const bool useEmbeddedBC = useEmbeddedBCMap || useEmbeddedBCFallback;

    // Per-row volume weight: legacy alpha*area, or revolved full-cell volume
    // for the validated strict axisymmetric port (which forbids cut metrics).
    std::vector<double> volOfRow(static_cast<size_t>(N), 0.0);
    for (int i = 0; i < N; ++i)
    {
        const int cellIndex = cellIndexOfRow[i];
        const double area = opt.axisymmetric
            ? CartesianMeasure::cell(mesh.cells(mesh.cellIndexToLabel[cellIndex]).area(),
                mesh.cells(mesh.cellIndexToLabel[cellIndex]).centerCoord().x, true)
            : safeCellArea(mesh, cellIndex);

        double alpha = 1.0;
        if (useCutMetrics)
        {
            alpha = cutMetrics->cellMetricsByIndex()[static_cast<size_t>(cellIndex)].alpha;
            if (!std::isfinite(alpha))
                alpha = 1.0;
            if (alpha < 0.0)
                alpha = 0.0;
            else if (alpha > 1.0)
                alpha = 1.0;
        }
        volOfRow[static_cast<size_t>(i)] = area * alpha;
    }

    std::unique_ptr<PolarizationFittedBoundaryDetail::RhsPlan> preparedRhs;
    if (fittedOperator && opt.reuseFittedFactorization)
        preparedRhs = std::make_unique<PolarizationFittedBoundaryDetail::RhsPlan>(mesh, graph,
            *fittedBoundaries, cellIndexOfRow, rowOfCellIndex, volOfRow, correctedOperator);

    // b_i = volumetric rhs * physical cell measure.
    Eigen::VectorXd bVec(N);
    bVec.setZero();
    const bool rhsOk = (rhs.size() == mesh.cellIndexNum()) && (rhs.value != nullptr);
    for (int i = 0; i < N; ++i)
    {
        const int cellIndex = cellIndexOfRow[i];
        const double f = rhsOk ? rhs.value[cellIndex] : 0.0;
        bVec[i] = f * volOfRow[static_cast<size_t>(i)];
    }

    const auto tAsm0 = Clock::now();

    std::vector<Eigen::Triplet<double>> trips;
    trips.reserve(segs.size() * 2u + static_cast<std::size_t>(N) + 32u);

    std::vector<double> diag(static_cast<size_t>(N), 0.0);
    std::vector<std::vector<int>> neighbors(static_cast<size_t>(N));
    std::vector<unsigned char> anchoredRow(static_cast<size_t>(N), 0u);

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

    bool warnedRobin = false;
    bool warnedPeriodicMismatch = false;

    for (size_t sid = 0; sid < segs.size(); ++sid)
    {
        const auto& seg = segs[sid];

        const int Kindex = mesh.cells(seg.K).index();
        if (Kindex < 0)
            continue;
        const int i = rowOfCellIndex[static_cast<size_t>(Kindex)];
        if (i < 0)
            continue;

        if (correctedOperator)
        {
            // Corrections already represent integrated physical flux. Their
            // RZ metric belongs to the caller and must not be applied twice.
            const auto& correction = (*fluxCorrections)[sid];
            for (const auto& coefficient : correction.coefficients)
            {
                const int j = rowOfCellIndex[static_cast<std::size_t>(coefficient.cellIndex)];
                trips.emplace_back(i, j, coefficient.value);
                if (j == i) diag[static_cast<std::size_t>(i)] += coefficient.value;
            }
            bVec[i] -= correction.constant;
        }

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

        const double AfOpen = opt.axisymmetric
            ? CartesianMeasure::face(seg.Af, seg.faceCenter.x, true)
            : beta * seg.Af;
        if (!(seg.df > 0.0) || !(AfOpen >= 0.0))
            continue;
        if (AfOpen == 0.0)
            continue;

        const double kappa_f = kappaForSegment(sid);

        if (fittedOperator)
        {
            const auto& boundary = (*fittedBoundaries)[sid];
            if (boundary.type == FittedBoundaryType::HomogeneousNeumann)
                continue; // Prescribed integrated flux is exactly zero.
            if (boundary.type == FittedBoundaryType::Dirichlet)
            {
                // An inactive neighbor is not an unknown at its center. The
                // supplied potential is on this actual physical face instead.
                const double distance = (seg.faceCenter - mesh.cells(seg.K).centerCoord()).dotProduct(seg.n);
                const double c = kappa_f * AfOpen / distance;
                if (preparedRhs) preparedRhs->recordDirichlet(sid, c, boundary.value);
                trips.emplace_back(i, i, c);
                diag[static_cast<std::size_t>(i)] += c;
                bVec[i] += c * boundary.value;
                anchoredRow[static_cast<std::size_t>(i)] = 1u;
                continue;
            }
        }

        if (seg.boundary == Boundary::non)
        {
            if (seg.L == -1)
                continue;

            const int Lindex = mesh.cells(seg.L).index();
            if (Lindex < 0)
                continue;

            const int j = rowOfCellIndex[static_cast<size_t>(Lindex)];
            if (j < 0)
                continue;

            const double c = kappa_f * AfOpen / seg.df;
            trips.emplace_back(i, i, c);
            trips.emplace_back(i, j, -c);
            diag[static_cast<size_t>(i)] += c;
            neighbors[static_cast<size_t>(i)].push_back(j);
            neighbors[static_cast<size_t>(j)].push_back(i);
        }
        else
        {
            const PicBoundarySide& s = bcSide(bc, seg.boundary);
            PicBCType t = s.type;

            if (t == PicBCType::Robin)
            {
                if (!warnedRobin)
                {
                    warnedRobin = true;
                    if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
                    {
                        std::cout << "[PolarizationPoissonFVSolverCellCentered] Warning: Robin BC encountered on domain boundary; "
                                  << "treating as Dirichlet(value).\n";
                    }
                }
                t = PicBCType::Dirichlet;
            }

            if (t == PicBCType::Dirichlet)
            {
                const double c = kappa_f * AfOpen / seg.df;
                trips.emplace_back(i, i, c);
                diag[static_cast<size_t>(i)] += c;
                bVec[i] += c * s.value;
                anchoredRow[static_cast<size_t>(i)] = 1u;
            }
            else if (t == PicBCType::Neumann)
            {
                // Interpret Neumann value as dphi/dn = g.
                // For -div(kappa grad phi)=rhs, boundary flux contribution scales by kappa.
                bVec[i] += (s.value * kappa_f) * AfOpen;
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
                            std::cout << "[PolarizationPoissonFVSolverCellCentered] Warning: Periodic BC on '"
                                      << boundaryToString(seg.boundary)
                                      << "' but mate side is not Periodic. Assembly will still wrap.\n";
                        }
                    }
                }

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
                if (!(df > 0.0))
                    continue;

                const double c = kappa_f * AfOpen / df;
                trips.emplace_back(i, i, c);
                trips.emplace_back(i, j, -c);
                diag[static_cast<size_t>(i)] += c;
                neighbors[static_cast<size_t>(i)].push_back(j);
                neighbors[static_cast<size_t>(j)].push_back(i);
            }
        }
    }

    // Embedded-boundary closure on cut cells.
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
                if (opt.embeddedBCMissingPolicy == Options::EmbeddedBCMissingPolicy::Error)
                {
                    if (report)
                    {
                        *report = {};
                        report->success = false;
                        report->message = "Missing embedded BC for boundaryId " + std::to_string(cm.gammaBoundaryId) +
                                          " (no map entry and no fallback).";
                        report->timeTotal_s = secondsSince(tBegin);
                        report->timeAssemble_s = report->timeTotal_s;
                        report->N = N;
                    }
                    return false;
                }
                continue;
            }

            const double kappa_g = kappaForCellIndex(cellIndex);
            const double a = eb->a;
            const double b = eb->b;
            const double g = eb->g;
            const PicBCType t = eb->type;
            const double A_gamma = cm.A_gamma;

            if (t == PicBCType::Robin)
            {
                if (std::abs(b) > 0.0)
                {
                    const double addDiag = kappa_g * (a / b) * A_gamma;
                    const double addRhs = kappa_g * (g / b) * A_gamma;
                    if (!std::isfinite(addDiag) || !std::isfinite(addRhs))
                        continue;
                    if (addDiag != 0.0)
                    {
                        trips.emplace_back(i, i, addDiag);
                        diag[static_cast<size_t>(i)] += addDiag;
                        if (a != 0.0)
                        {
                            anchoredRow[static_cast<size_t>(i)] = 1u;
                        }
                    }
                    bVec[i] += addRhs;
                }
                else
                {
                    if (std::abs(a) <= 0.0)
                        continue;
                    const double phiD = g / a;
                    const double d = stabilizedEmbeddedDistance(cm.d_gamma,
                                                                volOfRow[static_cast<size_t>(i)]);

                    const double c = kappa_g * (A_gamma / d);
                    if (!std::isfinite(c) || c <= 0.0)
                        continue;

                    trips.emplace_back(i, i, c);
                    diag[static_cast<size_t>(i)] += c;
                    bVec[i] += c * phiD;
                    anchoredRow[static_cast<size_t>(i)] = 1u;
                }
            }
            else if (t == PicBCType::Dirichlet)
            {
                const double phiD = g;
                const double d = stabilizedEmbeddedDistance(cm.d_gamma,
                                                            volOfRow[static_cast<size_t>(i)]);

                const double c = kappa_g * (A_gamma / d);
                if (!std::isfinite(c) || c <= 0.0)
                    continue;

                trips.emplace_back(i, i, c);
                diag[static_cast<size_t>(i)] += c;
                bVec[i] += c * phiD;
                anchoredRow[static_cast<size_t>(i)] = 1u;
            }
            else if (t == PicBCType::Neumann)
            {
                const double addRhs = (g * kappa_g) * A_gamma;
                if (!std::isfinite(addRhs))
                    continue;
                bVec[i] += addRhs;
            }
        }
    }

    // Stabilize any floating connected component. This covers fully-Neumann
    // problems and mixed-BC cases where moving EB temporarily disconnects a
    // plasma island from all Dirichlet/Robin anchors.
    std::vector<unsigned char> visited(static_cast<size_t>(N), 0u);
    std::vector<int> componentRows;
    componentRows.reserve(static_cast<size_t>(N));
    size_t floatingComponentCount = 0u;

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

        if (checkedOperator)
            throw std::runtime_error("fluxCorrections: unanchored component; compatibility correction/gauge pin is forbidden");

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
    }

    if ((floatingComponentCount > 0u) && (opt.verbose || ConsoleControl::enabled(ConsoleControl::Level::kDebug)))
    {
        std::cout << "[PolarizationPoissonFVSolverCellCentered] stabilized "
                  << floatingComponentCount
                  << " floating component(s) with compatibility correction + gauge pin.\n";
    }

    // Safety pin for zero/invalid diagonal.
    for (int i = 0; i < N; ++i)
    {
        const double d = diag[static_cast<size_t>(i)];
        if (!std::isfinite(d) || !(d > 0.0))
        {
            if (checkedOperator)
                throw std::runtime_error("fluxCorrections: non-positive/non-finite assembled diagonal; safety pin is forbidden");
            trips.emplace_back(i, i, 1.0);
            diag[static_cast<size_t>(i)] += 1.0;
            bVec[i] = 0.0;
        }
    }

    Eigen::SparseMatrix<double> A(N, N);
    A.setFromTriplets(trips.begin(), trips.end());
    A.makeCompressed();
    if (checkedOperator)
    {
        for (Eigen::Index entry = 0; entry < A.nonZeros(); ++entry)
            if (!std::isfinite(A.valuePtr()[entry]))
                throw std::runtime_error("fluxCorrections: assembled matrix contains non-finite entries");
        if (!bVec.allFinite())
            throw std::runtime_error("fluxCorrections: assembled RHS contains non-finite entries");
        for (int row = 0; row < N; ++row)
            if (!(A.coeff(row, row) > 0.0))
                throw std::runtime_error("fluxCorrections: complete assembled diagonal must be positive; no safety pin is allowed");
    }
    const bool useIterativeCG = opt.useCG || opt.useIC;
    const bool needPatternSig = !checkedOperator && (
        (useIterativeCG && (opt.reuseIterativePattern || opt.useTemporalExtrapolatedGuess)) ||
        ((!useIterativeCG) && (opt.reuseDirectPattern || opt.freezeDirectMatrixUntilPatternChange || opt.useSuperLU)));
    const std::uint64_t patternSig = needPatternSig ? sparsePatternSignature(A) : 0ull;

    const auto tAsm1 = Clock::now();
    const double tAssemble_s = secondsBetween(tAsm0, tAsm1);

    double tFactor_s = 0.0;
    double tSolve_s = 0.0;
    bool usedSuperLU = false;
    bool usedCG = false;
    const bool cgAttempted = useIterativeCG;
    bool cgFallbackToDirect = false;
    int cgIterations = 0;
    double cgResidual = 0.0;
    bool cgPatternReused = false;
    bool cgTemporalGuessUsed = false;
    bool directPatternReused = false;
    bool directMatrixFrozen = false;
    bool fittedFactorizationReused = false;
    std::string directFailReason;
    bool ok = false;
    Eigen::VectorXd x;

    if (useIterativeCG)
    {
        Eigen::VectorXd x0(N);
        const bool canUseTemporalGuess =
            opt.useTemporalExtrapolatedGuess &&
            iterCache_.hasLastX &&
            (iterCache_.xN == N) &&
            (iterCache_.xPatternSig == patternSig) &&
            (iterCache_.lastX.size() == N);

        if (canUseTemporalGuess)
        {
            x0 = iterCache_.lastX;
            if (iterCache_.hasPrevX && (iterCache_.prevX.size() == N))
            {
                const double w = std::max(0.0, opt.temporalExtrapolationOmega);
                if (w > 0.0)
                {
                    x0.noalias() += w * (iterCache_.lastX - iterCache_.prevX);
                    cgTemporalGuessUsed = true;
                }
            }
        }
        else
        {
            for (int i = 0; i < N; ++i)
            {
                const int cellIndex = cellIndexOfRow[i];
                x0[i] = (cellIndex >= 0 && cellIndex < mesh.cellIndexNum()) ? phi.value[cellIndex] : 0.0;
            }
        }

        if (opt.useIC)
        {
            auto& cg = iterCache_.cgIC;
            cg.setMaxIterations(std::max(1, opt.maxIter));
            cg.setTolerance(std::max(1e-16, opt.relTol));

            const bool canReuse =
                opt.reuseIterativePattern &&
                iterCache_.icAnalyzed &&
                iterCache_.icN == N &&
                iterCache_.icNnz == static_cast<long long>(A.nonZeros()) &&
                iterCache_.icPatternSig == patternSig;

            const auto tCgFact0 = Clock::now();
            bool haveFactorization = false;
            if (canReuse)
            {
                cg.factorize(A);
                if (cg.info() == Eigen::Success)
                {
                    haveFactorization = true;
                    cgPatternReused = true;
                }
            }
            if (!haveFactorization)
            {
                cg.compute(A);
                haveFactorization = (cg.info() == Eigen::Success);
                iterCache_.icAnalyzed = haveFactorization;
                if (haveFactorization)
                {
                    iterCache_.icN = N;
                    iterCache_.icNnz = static_cast<long long>(A.nonZeros());
                    iterCache_.icPatternSig = patternSig;
                }
            }
            const auto tCgFact1 = Clock::now();
            const auto tCgSol0 = tCgFact1;
            if (haveFactorization)
            {
                x = cg.solveWithGuess(bVec, x0);
                ok = (cg.info() == Eigen::Success);
                cgIterations = static_cast<int>(cg.iterations());
                cgResidual = cg.error();
            }
            const auto tCgSol1 = Clock::now();
            tFactor_s += secondsBetween(tCgFact0, tCgFact1);
            tSolve_s += secondsBetween(tCgSol0, tCgSol1);
        }
        else
        {
            auto& cg = iterCache_.cg;
            cg.setMaxIterations(std::max(1, opt.maxIter));
            cg.setTolerance(std::max(1e-16, opt.relTol));

            const bool canReuse =
                opt.reuseIterativePattern &&
                iterCache_.cgAnalyzed &&
                iterCache_.cgN == N &&
                iterCache_.cgNnz == static_cast<long long>(A.nonZeros()) &&
                iterCache_.cgPatternSig == patternSig;

            const auto tCgFact0 = Clock::now();
            bool haveFactorization = false;
            if (canReuse)
            {
                cg.factorize(A);
                if (cg.info() == Eigen::Success)
                {
                    haveFactorization = true;
                    cgPatternReused = true;
                }
            }
            if (!haveFactorization)
            {
                cg.compute(A);
                haveFactorization = (cg.info() == Eigen::Success);
                iterCache_.cgAnalyzed = haveFactorization;
                if (haveFactorization)
                {
                    iterCache_.cgN = N;
                    iterCache_.cgNnz = static_cast<long long>(A.nonZeros());
                    iterCache_.cgPatternSig = patternSig;
                }
            }
            const auto tCgFact1 = Clock::now();
            const auto tCgSol0 = tCgFact1;
            if (haveFactorization)
            {
                x = cg.solveWithGuess(bVec, x0);
                ok = (cg.info() == Eigen::Success);
                cgIterations = static_cast<int>(cg.iterations());
                cgResidual = cg.error();
            }
            const auto tCgSol1 = Clock::now();
            tFactor_s += secondsBetween(tCgFact0, tCgFact1);
            tSolve_s += secondsBetween(tCgSol0, tCgSol1);
        }

        if (ok)
        {
            for (int i = 0; i < N; ++i)
            {
                if (!std::isfinite(x[i]))
                {
                    ok = false;
                    break;
                }
            }
        }
        if (ok)
            usedCG = true;
    }

    if (!ok)
    {
        if (cgAttempted)
            cgFallbackToDirect = true;
        const auto tDirect0 = Clock::now();
        double directSolveWall_s = 0.0;
        const bool canReuseDirectPattern = opt.reuseDirectPattern && needPatternSig;
        const bool canFreezeDirect = opt.freezeDirectMatrixUntilPatternChange && needPatternSig;
        static int s_directDebugCached = -1;
        if (s_directDebugCached < 0)
            s_directDebugCached = envBool("QUANTUM_POL_DIRECT_FREEZE_DEBUG", false) ? 1 : 0;
        const bool directDebug = (s_directDebugCached == 1);
        static std::uint64_t s_directCallId = 0;
        const std::uint64_t directCallId = directDebug ? (++s_directCallId) : 0;
#ifdef EIGEN_SUPERLU_SUPPORT
        if (opt.useSuperLU)
        {
            // Default to a per-call SuperLU instance unless direct matrix freezing is enabled.
            // This avoids long-lived internal factorization workspaces accumulating across
            // many macro-steps in the common "refactorize every step" path.
            const bool superluLocal = checkedOperator || envBool("QUANTUM_POL_SUPERLU_LOCAL", !canFreezeDirect);

            auto configureSuperlu = [](auto& solver) {
                solver.options().ColPerm = COLAMD;
                solver.options().DiagPivotThresh = 0.1;
            };

            if (fittedOperator && opt.reuseFittedFactorization)
            {
                fittedFactorizationReused = fittedFactorCache_ &&
                    fittedFactorCache_->matches(A, mesh, cellIndexOfRow);
                if (!fittedFactorizationReused)
                {
                    // Destroy old factors before allocating/factorizing a new
                    // operator. A failed factorization never becomes reusable.
                    fittedFactorCache_.reset();
                    auto candidate = std::make_unique<FittedFactorCache>();
                    candidate->matrix = A;
                    candidate->mesh = &mesh;
                    candidate->cellIndices = cellIndexOfRow;
                    candidate->cellLabels.reserve(cellIndexOfRow.size());
                    for (int cell : cellIndexOfRow)
                        candidate->cellLabels.push_back(mesh.cellIndexToLabel[cell]);
                    configureSuperlu(candidate->solver);
                    candidate->solver.analyzePattern(candidate->matrix);
                    candidate->solver.factorize(candidate->matrix);
                    if (candidate->solver.info() == Eigen::Success)
                        fittedFactorCache_ = std::move(candidate);
                    else
                        directFailReason = "fitted SuperLU factorize failed, info=" +
                            std::to_string(static_cast<int>(candidate->solver.info()));
                }
                if (fittedFactorCache_)
                {
                    const auto tSolve0 = Clock::now();
                    x = fittedFactorCache_->solver.solve(bVec);
                    ok = fittedFactorCache_->solver.info() == Eigen::Success;
                    if (!ok)
                        directFailReason = "fitted SuperLU solve failed, info=" +
                            std::to_string(static_cast<int>(fittedFactorCache_->solver.info()));
                    directSolveWall_s += secondsBetween(tSolve0, Clock::now());
                    usedSuperLU = ok;
                }
            }
            else if (superluLocal)
            {
                Eigen::SuperLU<Eigen::SparseMatrix<double>> superlu;
                configureSuperlu(superlu);
                if (directDebug)
                    std::cout << "[PolPoisson][direct][SuperLU] call=" << directCallId << " LOCAL solver mode" << std::endl;

                if (directDebug)
                    std::cout << "[PolPoisson][direct][SuperLU] call=" << directCallId << " analyzePattern..." << std::endl;
                superlu.analyzePattern(A);
                if (directDebug)
                    std::cout << "[PolPoisson][direct][SuperLU] call=" << directCallId << " analyzePattern done (info="
                              << static_cast<int>(superlu.info()) << ")" << std::endl;

                if (directDebug)
                    std::cout << "[PolPoisson][direct][SuperLU] call=" << directCallId << " factorize..." << std::endl;
                superlu.factorize(A);
                if (directDebug)
                    std::cout << "[PolPoisson][direct][SuperLU] call=" << directCallId << " factorize done (info="
                              << static_cast<int>(superlu.info()) << ")" << std::endl;

                const bool factored = (superlu.info() == Eigen::Success);
                if (!factored)
                    directFailReason = "SuperLU factorize failed, info=" + std::to_string(static_cast<int>(superlu.info()));
                if (factored)
                {
                    const auto tSolve0 = Clock::now();
                    if (directDebug)
                        std::cout << "[PolPoisson][direct][SuperLU] call=" << directCallId << " solve..." << std::endl;
                    x = superlu.solve(bVec);
                    if (directDebug)
                        std::cout << "[PolPoisson][direct][SuperLU] call=" << directCallId << " solve done (info="
                                  << static_cast<int>(superlu.info()) << ")" << std::endl;
                    ok = (superlu.info() == Eigen::Success);
                    if (!ok)
                        directFailReason = "SuperLU solve failed, info=" + std::to_string(static_cast<int>(superlu.info()));
                    directSolveWall_s += secondsBetween(tSolve0, Clock::now());
                }
                if (ok)
                {
                    usedSuperLU = true;
                }
            }
            else
            {
                auto& superlu = directCache_.superLU;
                configureSuperlu(superlu);
            const bool samePattern =
                directCache_.superAnalyzed &&
                sameSparsePattern(directCache_.superMatrix, A);

            const bool reusePatternNow = samePattern && canReuseDirectPattern;
            const bool freezeNow =
                samePattern && canFreezeDirect && directCache_.superFactorized &&
                sameSparseValues(directCache_.superMatrix, A);
            bool factored = directCache_.superFactorized;

            if (directDebug)
            {
                std::cout << "[PolPoisson][direct][SuperLU]"
                          << " call=" << directCallId
                          << " N=" << N
                          << " nnz=" << static_cast<long long>(A.nonZeros())
                          << " samePattern=" << (samePattern ? 1 : 0)
                          << " reusePatternNow=" << (reusePatternNow ? 1 : 0)
                          << " hadFactor=" << (directCache_.superFactorized ? 1 : 0)
                          << " canFreeze=" << (canFreezeDirect ? 1 : 0)
                          << " freezeNow=" << (freezeNow ? 1 : 0)
                          << std::endl;
            }

            if (!freezeNow)
            {
                bool analyzed = reusePatternNow;
                directCache_.superAnalyzed = false;
                directCache_.superFactorized = false;
                if (!analyzed)
                {
                    if (directDebug)
                        std::cout << "[PolPoisson][direct][SuperLU] call=" << directCallId << " analyzePattern..." << std::endl;
                    superlu.analyzePattern(A);
                    if (directDebug)
                        std::cout << "[PolPoisson][direct][SuperLU] call=" << directCallId << " analyzePattern done (info="
                                  << static_cast<int>(superlu.info()) << ")" << std::endl;
                    analyzed = true;
                }
                if (analyzed)
                {
                    if (directDebug)
                    {
                        std::size_t nonFinite = 0u;
                        const double* vals = A.valuePtr();
                        const Eigen::Index nnz = A.nonZeros();
                        for (Eigen::Index k = 0; k < nnz; ++k)
                        {
                            if (!std::isfinite(vals[k]))
                                ++nonFinite;
                        }
                        if (nonFinite > 0u)
                        {
                            std::cout << "[PolPoisson][direct][SuperLU] call=" << directCallId
                                      << " WARNING: matrix has non-finite values, skipping SuperLU" << std::endl;
                            analyzed = false;
                        }
                    }
                    if (analyzed && directDebug)
                        std::cout << "[PolPoisson][direct][SuperLU] call=" << directCallId << " factorize..." << std::endl;
                    superlu.factorize(A);
                    if (directDebug)
                        std::cout << "[PolPoisson][direct][SuperLU] call=" << directCallId << " factorize done (info="
                                  << static_cast<int>(superlu.info()) << ")" << std::endl;
                    factored = (superlu.info() == Eigen::Success);
                    if (!factored && directFailReason.empty())
                    {
                        std::ostringstream oss;
                        oss << "SuperLU factorize failed, info=" << static_cast<int>(superlu.info());
                        directFailReason = oss.str();
                    }
                    if (factored)
                        directCache_.superMatrix = A;
                    directCache_.superAnalyzed = factored;
                    directCache_.superFactorized = factored;
                    if (factored && reusePatternNow)
                        directPatternReused = true;
                }
                else
                {
                    factored = false;
                    directCache_.superFactorized = false;
                }
            }
            else
            {
                directPatternReused = true;
                directMatrixFrozen = true;
            }

            if (factored)
            {
                const auto tSolve0 = Clock::now();
                if (directDebug)
                    std::cout << "[PolPoisson][direct][SuperLU] call=" << directCallId << " solve..." << std::endl;
                x = superlu.solve(bVec);
                if (directDebug)
                    std::cout << "[PolPoisson][direct][SuperLU] call=" << directCallId << " solve done (info="
                              << static_cast<int>(superlu.info()) << ")" << std::endl;
                ok = (superlu.info() == Eigen::Success);
                directSolveWall_s += secondsBetween(tSolve0, Clock::now());
                if (!ok && directFailReason.empty())
                {
                    std::ostringstream oss;
                    oss << "SuperLU solve failed, info=" << static_cast<int>(superlu.info());
                    directFailReason = oss.str();
                }
            }
            if (ok)
            {
                usedSuperLU = true;
            }
            else
            {
                usedSuperLU = false;
                directCache_.superFactorized = false;
            }
            }
        }
#endif

        if (opt.requireSuperLU && ok && (!usedSuperLU || !x.allFinite()))
        {
            ok = false;
            usedSuperLU = false;
            directFailReason = "required SuperLU solve returned non-finite values or did not use SuperLU";
        }

        if (!ok && !opt.requireSuperLU)
        {
            auto& sparselu = directCache_.sparseLU;
            const bool samePattern =
                directCache_.sparseAnalyzed &&
                sameSparsePattern(directCache_.sparseMatrix, A);

            const bool reusePatternNow = samePattern && canReuseDirectPattern;
            const bool freezeNow =
                samePattern && canFreezeDirect && directCache_.sparseFactorized &&
                sameSparseValues(directCache_.sparseMatrix, A);
            bool factored = directCache_.sparseFactorized;

            if (directDebug)
            {
                std::cout << "[PolPoisson][direct][SparseLU]"
                          << " call=" << directCallId
                          << " N=" << N
                          << " nnz=" << static_cast<long long>(A.nonZeros())
                          << " samePattern=" << (samePattern ? 1 : 0)
                          << " reusePatternNow=" << (reusePatternNow ? 1 : 0)
                          << " hadFactor=" << (directCache_.sparseFactorized ? 1 : 0)
                          << " canFreeze=" << (canFreezeDirect ? 1 : 0)
                          << " freezeNow=" << (freezeNow ? 1 : 0)
                          << std::endl;
            }

            if (!freezeNow)
            {
                bool analyzed = reusePatternNow;
                directCache_.sparseAnalyzed = false;
                directCache_.sparseFactorized = false;
                if (!analyzed)
                {
                    sparselu.analyzePattern(A);
                    analyzed = true;
                }
                if (analyzed)
                {
                    sparselu.factorize(A);
                    factored = (sparselu.info() == Eigen::Success);
                    if (!factored && directFailReason.empty())
                    {
                        std::ostringstream oss;
                        oss << "SparseLU factorize failed, info=" << static_cast<int>(sparselu.info());
                        directFailReason = oss.str();
                    }
                    if (factored)
                        directCache_.sparseMatrix = A;
                    directCache_.sparseAnalyzed = factored;
                    directCache_.sparseFactorized = factored;
                    if (factored && reusePatternNow)
                        directPatternReused = true;
                }
                else
                {
                    factored = false;
                    directCache_.sparseFactorized = false;
                }
            }
            else
            {
                directPatternReused = true;
                directMatrixFrozen = true;
            }

            if (factored)
            {
                const auto tSolve0 = Clock::now();
                x = sparselu.solve(bVec);
                ok = (sparselu.info() == Eigen::Success);
                directSolveWall_s += secondsBetween(tSolve0, Clock::now());
                if (!ok && directFailReason.empty())
                {
                    std::ostringstream oss;
                    oss << "SparseLU solve failed, info=" << static_cast<int>(sparselu.info());
                    directFailReason = oss.str();
                }
            }
            if (!ok)
            {
                directCache_.sparseFactorized = false;
            }
        }
        const double directTotal_s = secondsBetween(tDirect0, Clock::now());
        tSolve_s += directSolveWall_s;
        tFactor_s += std::max(0.0, directTotal_s - directSolveWall_s);
        usedCG = false;
        if (!cgAttempted)
        {
            cgIterations = 0;
            cgResidual = 0.0;
        }
    }

    double correctedResidual = 0.0;
    if (checkedOperator && ok)
    {
        // Check the operator actually factorized, including every integrated
        // flux correction and boundary constant. Do not use CG's estimate or a
        // residual of only the underlying two-point diffusion stencil.
        ok = fullResidualAccepted(A, x, bVec, usedSuperLU, opt.relTol,
            correctedResidual, directFailReason);
    }

    if (!ok)
    {
        if (report)
        {
            *report = {};
            report->success = false;
            report->finalResidual = correctedResidual;
            if (opt.requireSuperLU)
                report->message = "requireSuperLU solve failed; SparseLU fallback is disabled: " +
                    (directFailReason.empty() ? std::string("SuperLU factorization/solve failed") : directFailReason);
            else if (!directFailReason.empty())
                report->message = "Polarization Poisson solve failed (CG and direct fallback): " + directFailReason;
            else
                report->message = "Polarization Poisson solve failed (CG and direct fallback).";
            report->timeAssemble_s = tAssemble_s;
            report->timeFactor_s = tFactor_s;
            report->timeSolve_s = tSolve_s;
            report->timeTotal_s = secondsSince(tBegin);
            report->N = N;
            report->nnz = static_cast<long long>(A.nonZeros());
            report->usedSuperLU = usedSuperLU;
            report->usedCG = usedCG;
            report->cgAttempted = cgAttempted;
            report->cgFallbackToDirect = cgFallbackToDirect;
            report->cgIterations = cgIterations;
            report->cgResidual = cgResidual;
            report->cgPatternReused = cgPatternReused;
            report->cgTemporalGuessUsed = cgTemporalGuessUsed;
            report->directPatternReused = directPatternReused;
            report->directMatrixFrozen = directMatrixFrozen;
            report->fittedFactorizationReused = fittedFactorizationReused;
        }
        return false;
    }

    if (useIterativeCG && opt.useTemporalExtrapolatedGuess)
    {
        const bool haveCompatibleLast =
            iterCache_.hasLastX &&
            (iterCache_.xN == N) &&
            (iterCache_.xPatternSig == patternSig) &&
            (iterCache_.lastX.size() == N);
        if (haveCompatibleLast)
        {
            iterCache_.prevX = iterCache_.lastX;
            iterCache_.hasPrevX = true;
        }
        else
        {
            iterCache_.prevX.resize(0);
            iterCache_.hasPrevX = false;
        }
        iterCache_.lastX = x;
        iterCache_.xN = N;
        iterCache_.xPatternSig = patternSig;
        iterCache_.hasLastX = true;
    }

    if (preparedRhs)
    {
        if (!fittedFactorCache_)
            throw std::runtime_error("fitted RHS plan has no successful SuperLU factorization");
        fittedFactorCache_->rhsPlan = std::move(preparedRhs);
    }

    // New-port validation, factorization, solve and residual failures have not
    // touched phi. Allocate a missing output only at the successful commit.
    if (correctedOperator && (phi.size() != mesh.cellIndexNum() || !phi.value))
        phi.initialize(mesh, 0.0);

    for (int i = 0; i < N; ++i)
    {
        const int cellIndex = cellIndexOfRow[i];
        if (cellIndex >= 0 && cellIndex < mesh.cellIndexNum())
            phi.value[cellIndex] = x[i];
    }

    if (report)
    {
        *report = {};
        report->success = true;
        report->finalResidual = correctedResidual;
        report->message = fittedOperator && correctedOperator ?
            "Solved face-fitted active-leaf Poisson with integrated flux corrections; full-operator residual verified." :
            correctedOperator ?
            "Solved variable-coefficient Poisson with integrated flux corrections; full-operator residual verified." :
            fittedOperator ? "Solved face-fitted active-leaf Poisson; full-operator residual verified." :
            "Solved variable-coefficient Poisson (one-shot).";
        report->timeAssemble_s = tAssemble_s;
        report->timeFactor_s = tFactor_s;
        report->timeSolve_s = tSolve_s;
        report->timeTotal_s = secondsSince(tBegin);
        report->N = N;
        report->nnz = static_cast<long long>(A.nonZeros());
        report->usedSuperLU = usedSuperLU;
        report->usedCG = usedCG;
        report->cgAttempted = cgAttempted;
        report->cgFallbackToDirect = cgFallbackToDirect;
        report->cgIterations = cgIterations;
        report->cgResidual = cgResidual;
        report->cgPatternReused = cgPatternReused;
        report->cgTemporalGuessUsed = cgTemporalGuessUsed;
        report->directPatternReused = directPatternReused;
        report->directMatrixFrozen = directMatrixFrozen;
        report->fittedFactorizationReused = fittedFactorizationReused;
    }

    return true;
}

void PolarizationPoissonFVSolverCellCentered::clearDirectCache() const
{
    fittedFactorCache_.reset();
    iterCache_.icN = -1;
    iterCache_.icNnz = -1;
    iterCache_.icPatternSig = 0;
    iterCache_.icAnalyzed = false;
    iterCache_.cgN = -1;
    iterCache_.cgNnz = -1;
    iterCache_.cgPatternSig = 0;
    iterCache_.cgAnalyzed = false;
    iterCache_.xN = -1;
    iterCache_.xPatternSig = 0;
    iterCache_.hasLastX = false;
    iterCache_.hasPrevX = false;
    iterCache_.lastX.resize(0);
    iterCache_.prevX.resize(0);

#ifdef EIGEN_SUPERLU_SUPPORT
    directCache_.superLU.~SuperLU();
    new (&directCache_.superLU) Eigen::SuperLU<Eigen::SparseMatrix<double>>();
    directCache_.superMatrix.resize(0, 0);
    directCache_.superAnalyzed = false;
    directCache_.superFactorized = false;
#endif
    directCache_.sparseLU.~SparseLU();
    new (&directCache_.sparseLU) Eigen::SparseLU<Eigen::SparseMatrix<double>>();
    directCache_.sparseMatrix.resize(0, 0);
    directCache_.sparseAnalyzed = false;
    directCache_.sparseFactorized = false;
}
