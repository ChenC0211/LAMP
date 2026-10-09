#include "Modules/Hydrodynamics/Hydrodynamics.h"
#include <cstring>
#include <sstream>
#ifdef EIGEN_SUPERLU_SUPPORT
#include <Eigen/SuperLUSupport>
#endif

using namespace std;

struct Hydrodynamics::LinearSolveCache
{
#ifdef EIGEN_SUPERLU_SUPPORT
    // Own the exact compressed matrix as well as its factorization. No factor
    // can survive a changed coefficient, topology, dimension or failed compute.
    Eigen::SparseMatrix<double> matrix;
    Eigen::SuperLU<Eigen::SparseMatrix<double>> solver;

    bool matches(const Eigen::SparseMatrix<double>& candidate) const
    {
        if (!candidate.isCompressed() || matrix.rows() != candidate.rows()
            || matrix.cols() != candidate.cols() || matrix.nonZeros() != candidate.nonZeros())
            return false;
        using StorageIndex = Eigen::SparseMatrix<double>::StorageIndex;
        const auto outerBytes = static_cast<std::size_t>(matrix.outerSize() + 1) * sizeof(StorageIndex);
        if (std::memcmp(matrix.outerIndexPtr(), candidate.outerIndexPtr(), outerBytes) != 0)
            return false;
        if (matrix.nonZeros() == 0) return true;
        const auto entries = static_cast<std::size_t>(matrix.nonZeros());
        return std::memcmp(matrix.innerIndexPtr(), candidate.innerIndexPtr(), entries * sizeof(StorageIndex)) == 0
            && std::memcmp(matrix.valuePtr(), candidate.valuePtr(), entries * sizeof(double)) == 0;
    }
#endif
};

Hydrodynamics::~Hydrodynamics() { discardStepTransaction(); discardPendingStep(); }

void Hydrodynamics::requireLegacyMaterialEntry() const
{
    if (initialStateSet_ || materialState_ || stepProgress_.stopped)
        throw std::logic_error("managed Hydro materials require updateMaterials with source time and geometry");
}

std::optional<HydroSolveReport> Hydrodynamics::claimStepEntry(StepEntry entry)
{
    if (stepEntry_ != StepEntry::Unset && stepEntry_ != entry)
    {
        HydroSolveReport rejected;
        rejected.stage = "step-entry";
        rejected.reason = "cannot mix legacy solve/cache and managed advance on one Hydro instance";
        rejected.requiresReinitialization = stepProgress_.stopped;
        return rejected;
    }
    stepEntry_ = entry;
    return std::nullopt;
}

Eigen::VectorXd Hydrodynamics::checkedLinearSolve(const Eigen::SparseMatrix<double>& matrix,
    const Eigen::VectorXd& rhs, const std::string& stage, double& residual,
    std::optional<double> relativeTolerance,int maximumRefinements)
{
    const double tolerance = relativeTolerance.value_or(activeOptions_.linearResidualTolerance);
    if (!std::isfinite(tolerance) || tolerance <= 0.0)
        throw std::invalid_argument("linear relative tolerance must be finite and positive");
    if(maximumRefinements<0 || maximumRefinements>4)
        throw std::invalid_argument("linear residual refinement count must be in [0,4]");
    activeReport_.stage = stage + "-assembly";
    if (matrix.rows() <= 0 || matrix.rows() != matrix.cols() || rhs.size() != matrix.rows())
        throw std::runtime_error("linear system must be nonempty, square and match the right-hand side");
    if (!rhs.allFinite())
        throw std::runtime_error("right-hand side contains non-finite values");
    for (int column = 0; column < matrix.outerSize(); ++column)
        for (Eigen::SparseMatrix<double>::InnerIterator entry(matrix, column); entry; ++entry)
            if (!std::isfinite(entry.value()))
                throw std::runtime_error("matrix contains non-finite values");

    activeReport_.stage = stage + "-factorization";
#ifndef EIGEN_SUPERLU_SUPPORT
    throw std::runtime_error("SuperLU backend is required for Hydro but was not enabled in this build; SparseLU fallback is disabled");
#else
    Eigen::SparseMatrix<double> compressed;
    const Eigen::SparseMatrix<double>* candidate = &matrix;
    if (!matrix.isCompressed())
    {
        compressed = matrix;
        compressed.makeCompressed();
        candidate = &compressed;
    }
    // Unknown internal stages retain one-shot behavior instead of sharing a
    // factor with any of the three independently cached physical operators.
    std::shared_ptr<LinearSolveCache> transient;
    auto& cache = stage == "momentum" ? momentumLinearCache_
        : stage == "pressure" ? pressureLinearCache_
        : stage == "coupled" ? coupledLinearCache_ : transient;
    const bool reused = cache && cache->matches(*candidate);
    if (!reused)
    {
        // Clear first: a failed allocation/factorization must not make an old
        // operator available as a fallback for this different matrix.
        cache.reset();
        auto next = std::make_shared<LinearSolveCache>();
        next->matrix = *candidate;
        next->matrix.makeCompressed();
        next->solver.compute(next->matrix);
        if (next->solver.info() != Eigen::Success)
            throw std::runtime_error("SuperLU factorization failed, info="
                + std::to_string(static_cast<int>(next->solver.info())));
        cache = std::move(next);
    }

    activeReport_.stage = stage + "-solve";
    Eigen::VectorXd result = cache->solver.solve(rhs);
    if (cache->solver.info() != Eigen::Success || !result.allFinite())
    {
        const int info = static_cast<int>(cache->solver.info());
        cache.reset();
        throw std::runtime_error("SuperLU solve failed or returned non-finite values, info=" + std::to_string(info));
    }
    const double rhsNorm = rhs.stableNorm();
    Eigen::VectorXd defect=rhs-matrix*result;
    double errorNorm = defect.stableNorm();
    // Opt-in for conservative AMR projection only. Reuse its checked factor
    // to solve the actual residual; the original acceptance tolerance stays
    // unchanged. Ordinary physical solves retain their one-solve behavior.
    for(int pass=0;pass<maximumRefinements && std::isfinite(errorNorm)
        && errorNorm>tolerance*rhsNorm;++pass) {
        const Eigen::VectorXd delta=cache->solver.solve(defect);
        if(cache->solver.info()!=Eigen::Success || !delta.allFinite()) {
            cache.reset();throw std::runtime_error("SuperLU residual correction failed or returned non-finite values");
        }
        result+=delta;defect=rhs-matrix*result;errorNorm=defect.stableNorm();
    }
    if (!std::isfinite(rhsNorm) || !std::isfinite(errorNorm))
    {
        residual = std::numeric_limits<double>::quiet_NaN();
        cache.reset();
        throw std::runtime_error("linear right-hand-side or residual norm is non-finite");
    }
    residual = rhsNorm > 0.0 ? errorNorm / rhsNorm : errorNorm;
    if (!std::isfinite(residual) || residual > tolerance)
    {
        cache.reset();
        std::ostringstream message;
        message << std::scientific << "relative linear residual " << residual
            << " exceeds tolerance " << tolerance;
        throw std::runtime_error(message.str());
    }
    if (stage == "momentum")
    {
        activeReport_.momentumUsedSuperLU = true;
        activeReport_.momentumFactorReused = reused;
    }
    if (stage == "pressure")
    {
        activeReport_.pressureUsedSuperLU = true;
        activeReport_.pressureFactorReused = reused;
    }
    if (stage == "coupled")
    {
        activeReport_.coupledUsedSuperLU = true;
        activeReport_.coupledFactorReused = reused;
    }
    return result;
#endif
}
