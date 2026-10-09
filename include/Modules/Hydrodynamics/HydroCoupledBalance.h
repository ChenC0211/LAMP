#pragma once

#include <Eigen/Sparse>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

// Shared closure of the already assembled CN momentum and FV continuity
// equations. Geometry/physical boundary strategies provide C, D and their
// evaluated loads; this service owns the coupled block layout and residual.
// It does not select a time scheme, pressure boundary, solver or tolerance.
// Inline definitions preserve the caller's existing floating-point boundary.
namespace HydroCoupledBalance
{
class Assembly;

// One frozen coupled operator for one assembled momentum equation. Owns the
// block matrix, borrows the unchanged momentum matrix for residual evaluation.
// Both native and sparse strategies use this same assembly/RHS lifetime; a
// fitted balance refinement changes only the RHS, never the operator/gauge.
// The momentum owner (and C/D for a sparse strategy) must outlive this object
// and keep those operators unchanged. No physical field is retained.
class Operator
{
public:
    Operator(const Operator&) = delete;
    Operator& operator=(const Operator&) = delete;
    Operator(Operator&&) = default;
    const Eigen::SparseMatrix<double>& matrix() const noexcept { return matrix_; }

    // pressureLoad = C*pi; continuityDefect = D*u - b_boundary.
    // Preserve native accumulation order, including (V*dp)/h vs V*(dp/h).
    // Keep sparse C*pi as an expression until assignment: materializing it
    // changes Eigen's b-A*u-C*pi evaluation, including signed zero.
    template<class PressureLoad>
    Eigen::VectorXd rightHandSide(const Eigen::VectorXd& momentumRhs, const Eigen::VectorXd& velocity,
        const Eigen::MatrixBase<PressureLoad>& pressureLoad, const Eigen::VectorXd& continuityDefect) const
    {
        if (matrix_.rows() != velocities_ + cells_ || momentum_.rows() != velocities_
            || momentum_.cols() != velocities_ || momentumRhs.size() != velocities_ || velocity.size() != velocities_
            || pressureLoad.size() != velocities_ || continuityDefect.size() != cells_
            || !momentumRhs.allFinite() || !velocity.allFinite() || !continuityDefect.allFinite())
            throw std::invalid_argument("invalid Hydro coupled correction loads");
        Eigen::VectorXd rhs(velocities_ + cells_);
        rhs.head(velocities_) = momentumRhs - momentum_ * velocity - pressureLoad.derived();
        rhs.tail(cells_) = continuityDefect;
        if (reference_) rhs[velocities_ + *reference_] = 0.0;
        if (!rhs.allFinite()) throw std::invalid_argument("non-finite Hydro coupled correction rhs");
        return rhs;
    }

    // Sparse strategies retain their actual C/D, so refinements cannot bind
    // a different geometry's loads to this step's block matrix. Validate pi
    // itself as before, including components in a nullspace of C.
    Eigen::VectorXd rightHandSide(const Eigen::VectorXd& momentumRhs, const Eigen::VectorXd& velocity,
        const Eigen::VectorXd& pressureIncrement) const
    {
        if (!gradient_ || !divergence_ || velocity.size() != velocities_
            || pressureIncrement.size() != cells_ || !pressureIncrement.allFinite() || !velocity.allFinite())
            throw std::invalid_argument("invalid Hydro coupled sparse operands");
        const Eigen::VectorXd continuityDefect = *divergence_ * velocity;
        return rightHandSide(momentumRhs, velocity, *gradient_ * pressureIncrement, continuityDefect);
    }

private:
    friend class Assembly;
    friend Operator assemble(const Eigen::SparseMatrix<double>&,
        const Eigen::SparseMatrix<double>&, const Eigen::SparseMatrix<double>&);
    Operator(const Eigen::SparseMatrix<double>& momentum, int cells, std::optional<int> reference,
        Eigen::SparseMatrix<double>&& matrix)
        : momentum_(momentum), cells_(cells), velocities_(static_cast<int>(momentum.rows())),
          reference_(reference), matrix_(std::move(matrix)) {}
    const Eigen::SparseMatrix<double>& momentum_;
    int cells_, velocities_;
    std::optional<int> reference_;
    Eigen::SparseMatrix<double> matrix_;
    const Eigen::SparseMatrix<double>* gradient_ = nullptr;
    const Eigen::SparseMatrix<double>* divergence_ = nullptr;
};

class Assembly
{
public:
    // Short-lived borrow: momentum must outlive this assembly and remain
    // unchanged. No geometry/field copy or persistent matrix cache is created.
    // nullopt retains every continuity row (physical pressure anchoring);
    // otherwise replace exactly that redundant row with delta_pi[cell] = 0.
    // eliminateReferenceColumn additionally removes the known-zero pressure
    // column, retaining the free-block adjoint relation for CCH assembly.
    Assembly(const Eigen::SparseMatrix<double>& momentum, int cells,
        std::size_t couplingEntriesHint, std::optional<int> pressureReference = std::nullopt,
        bool eliminateReferenceColumn = false)
        : momentum_(momentum), cells_(cells), reference_(pressureReference), eliminateReferenceColumn_(eliminateReferenceColumn)
    {
        require(momentum.rows() > 0 && momentum.rows() == momentum.cols()
            && cells > 0 && momentum.rows() <= std::numeric_limits<int>::max() - cells,
            "invalid Hydro coupled operator dimensions");
        velocities_ = static_cast<int>(momentum.rows());
        require(!reference_ || (*reference_ >= 0 && *reference_ < cells_),
            "invalid Hydro coupled pressure-reference row");
        const auto momentumEntries = static_cast<std::size_t>(momentum.nonZeros());
        require(momentumEntries < entries_.max_size()
            && couplingEntriesHint <= entries_.max_size() - momentumEntries - 1,
            "Hydro coupled entry count overflow");
        entries_.reserve(momentumEntries + couplingEntriesHint + 1);
        for (int column = 0; column < momentum.outerSize(); ++column)
            for (Eigen::SparseMatrix<double>::InnerIterator entry(momentum, column); entry; ++entry)
            {
                require(std::isfinite(entry.value()), "non-finite Hydro momentum coefficient");
                entries_.emplace_back(static_cast<int>(entry.row()), static_cast<int>(entry.col()), entry.value());
            }
    }
    Assembly(Eigen::SparseMatrix<double>&&, int, std::size_t, std::optional<int> = std::nullopt, bool = false) = delete;
    Assembly(const Assembly&) = delete;
    Assembly& operator=(const Assembly&) = delete;

    void addPressureGradient(int velocityRow, int cell, double coefficient)
    {
        checkEntry(velocityRow, cell, coefficient);
        // delta_p[reference]=0. Removing its column as well as the redundant
        // continuity row preserves adjoint coupling on the free unknowns.
        if(eliminateReferenceColumn_ && reference_ && cell==*reference_)return;
        entries_.emplace_back(velocityRow, velocities_ + cell, coefficient);
    }

    // coefficient belongs to D, not -D. The shared lower block is -D.
    void addDivergence(int cell, int velocityRow, double coefficient)
    {
        checkEntry(velocityRow, cell, coefficient);
        if (!reference_ || cell != *reference_)
            entries_.emplace_back(velocities_ + cell, velocityRow, -coefficient);
    }

    // Finish once. No further rows or second finalization are allowed, even
    // after a failed allocation. Returned operator owns the assembled block.
    Operator finish()
    {
        // Keep the gauge insertion after native row assembly, as in the
        // complete-domain kernel. finish is a single assembly operation.
        require(!finished_, "Hydro coupled assembly already finished");
        finished_ = true;
        if (reference_) entries_.emplace_back(velocities_ + *reference_, velocities_ + *reference_, 1.0);
        Eigen::SparseMatrix<double> matrix(velocities_ + cells_, velocities_ + cells_);
        matrix.setFromTriplets(entries_.begin(), entries_.end());
        matrix.makeCompressed();
        return Operator(momentum_, cells_, reference_, std::move(matrix));
    }

private:
    static void require(bool condition, const char* message)
    {
        if (!condition) throw std::invalid_argument(message);
    }
    void checkEntry(int velocityRow, int cell, double coefficient) const
    {
        require(!finished_, "Hydro coupled assembly already finished");
        require(velocityRow >= 0 && velocityRow < velocities_ && cell >= 0 && cell < cells_
            && std::isfinite(coefficient), "invalid Hydro coupled gradient/divergence entry");
    }
    const Eigen::SparseMatrix<double>& momentum_;
    int cells_ = 0, velocities_ = 0;
    std::optional<int> reference_;
    std::vector<Eigen::Triplet<double>> entries_;
    bool finished_ = false;
    bool eliminateReferenceColumn_ = false;
};

// Sparse native operators are a row strategy, not a second coupled solver.
// Used by the fitted domain, including its non-adjoint coarse/fine C and D.
inline Operator assemble(const Eigen::SparseMatrix<double>& momentum,
    const Eigen::SparseMatrix<double>& gradient, const Eigen::SparseMatrix<double>& divergence)
{
    if (divergence.rows() > std::numeric_limits<int>::max()
        || gradient.rows() != momentum.rows() || gradient.cols() != divergence.rows()
        || divergence.cols() != momentum.rows())
        throw std::invalid_argument("invalid Hydro coupled sparse operands");
    Assembly assembly(momentum, static_cast<int>(divergence.rows()),
        static_cast<std::size_t>(gradient.nonZeros()) + static_cast<std::size_t>(divergence.nonZeros()));
    for (int column = 0; column < gradient.outerSize(); ++column)
        for (Eigen::SparseMatrix<double>::InnerIterator entry(gradient, column); entry; ++entry)
            assembly.addPressureGradient(static_cast<int>(entry.row()), static_cast<int>(entry.col()), entry.value());
    for (int column = 0; column < divergence.outerSize(); ++column)
        for (Eigen::SparseMatrix<double>::InnerIterator entry(divergence, column); entry; ++entry)
            assembly.addDivergence(static_cast<int>(entry.row()), static_cast<int>(entry.col()), entry.value());
    auto equations = assembly.finish();
    equations.gradient_ = &gradient;
    equations.divergence_ = &divergence;
    return equations;
}
inline Operator assemble(Eigen::SparseMatrix<double>&&,
    const Eigen::SparseMatrix<double>&, const Eigen::SparseMatrix<double>&) = delete;
inline Operator assemble(const Eigen::SparseMatrix<double>&,
    Eigen::SparseMatrix<double>&&, const Eigen::SparseMatrix<double>&) = delete;
inline Operator assemble(const Eigen::SparseMatrix<double>&,
    const Eigen::SparseMatrix<double>&, Eigen::SparseMatrix<double>&&) = delete;

struct MomentumResidual
{
    double relative = 0.0;
    double maxForceDensity = 0.0;
};

// The original unsplit momentum norm, shared by complete and fitted domains.
// A zero dual volume marks a prescribed-velocity constraint: its algebraic
// defect still enters the relative norm, but is not a physical force density.
// A borrowed volume query avoids allocating another per-face geometry array.
template<class DualVolumeAt>
MomentumResidual momentumResidual(const Eigen::SparseMatrix<double>& momentum,
    const Eigen::VectorXd& velocity, const Eigen::VectorXd& pressureLoad,
    const Eigen::VectorXd& momentumRhs, DualVolumeAt dualVolumeAt)
{
    if (momentum.rows() <= 0 || momentum.rows() != momentum.cols()
        || velocity.size() != momentum.rows() || pressureLoad.size() != momentum.rows()
        || momentumRhs.size() != momentum.rows() || !velocity.allFinite()
        || !pressureLoad.allFinite() || !momentumRhs.allFinite())
        throw std::invalid_argument("invalid Hydro momentum residual operands");
    const Eigen::VectorXd momentumLoad = momentum * velocity;
    const Eigen::VectorXd defect = momentumLoad + pressureLoad - momentumRhs;
    const double scale = std::max(momentumRhs.stableNorm(), momentumLoad.stableNorm() + pressureLoad.stableNorm());
    MomentumResidual result;
    result.relative = scale > 0.0 ? defect.stableNorm() / scale : defect.stableNorm();
    for (int i = 0; i < velocity.size(); ++i)
    {
        const double volume = dualVolumeAt(i);
        if (!std::isfinite(volume) || volume < 0.0)
            throw std::invalid_argument("invalid Hydro momentum residual dual volume");
        if (volume > 0.0)
            result.maxForceDensity = std::max(result.maxForceDensity, std::abs(defect[i] / volume));
    }
    if (!defect.allFinite() || !std::isfinite(scale) || !std::isfinite(result.relative)
        || !std::isfinite(result.maxForceDensity))
        throw std::runtime_error("non-finite Hydro unsplit momentum residual");
    return result;
}
}
