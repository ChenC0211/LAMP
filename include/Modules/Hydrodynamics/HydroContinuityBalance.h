#pragma once

#include <Eigen/Sparse>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>
#include "Core/Field/CellField.h"
#include "Core/Field/HEdgeField.h"
#include "Core/Field/VEdgeField.h"
#include "Core/Mesh/CartesianMeasure.h"

// FV volume-flux diagnostics of a supplied velocity, not a pressure solver or
// an acceptance policy. SI/code units follow the supplied velocity and actual
// XY/RZ measures. No App config, time loop, mesh copy or persistent cache.
namespace HydroContinuityBalance
{
// One immutable evaluation of the original compact divergence operator.
// Re-evaluate after EVERY velocity correction; no borrowed velocity/operator
// survives. Prescribed boundary flux is already represented by the supplied
// domain's unknowns. Do not reuse this zero-source evaluation for a different
// continuity equation with an omitted nonzero boundary/source RHS.
class LocalFlux
{
public:
    LocalFlux(const Eigen::SparseMatrix<double>& divergence, const Eigen::VectorXd& velocity)
    {
        if (divergence.rows() <= 0 || divergence.cols() != velocity.size()
            || velocity.size() <= 0 || !velocity.allFinite())
            throw std::invalid_argument("invalid Hydro continuity operator/velocity");
        flux_ = divergence * velocity;
        maxFaceSpeed_ = velocity.cwiseAbs().maxCoeff();
        scale_ = divergence.cwiseAbs()
            * Eigen::VectorXd::Constant(velocity.size(), maxFaceSpeed_);
        if (!flux_.allFinite() || !scale_.allFinite())
            throw std::runtime_error("non-finite Hydro continuity flux or scale");
    }

    const Eigen::VectorXd& netFlux() const noexcept { return flux_; }
    const Eigen::VectorXd& referenceFlux() const noexcept { return scale_; }
    double maxFaceSpeed() const noexcept { return maxFaceSpeed_; }

    // Original local ratio: no absolute floor, predictor scale or density.
    // A zero scale permits only exactly zero flux.
    double relativeResidual() const
    {
        double result = 0.0;
        for (int i = 0; i < flux_.size(); ++i)
        {
            const double relative = scale_[i] > 0 ? std::abs(flux_[i]) / scale_[i]
                : flux_[i] == 0 ? 0 : std::numeric_limits<double>::infinity();
            result = std::max(result, relative);
        }
        return result;
    }

    // Preserve the NS gate abs(q) <= tolerance*scale. Comparing the rounded
    // ratio to tolerance is NOT interchangeable at a threshold. AMR retains
    // its separately selected relativeResidual() gate and correction count.
    bool satisfiesScaledFlux(double tolerance) const
    {
        if (!std::isfinite(tolerance) || tolerance < 0.0)
            throw std::invalid_argument("invalid Hydro local-continuity tolerance");
        bool accepted = true;
        for (int i = 0; i < flux_.size(); ++i)
            accepted = accepted && std::abs(flux_[i]) <= tolerance * scale_[i];
        return accepted;
    }

    double maxDivergence(const std::vector<double>& volumes) const
    { return maximumDivergence(flux_, volumes); }

    // Also usable for a predictor when only dimensional divergence is needed;
    // it does not build or impose a local residual scale.
    static double maximumDivergence(const Eigen::VectorXd& flux, const std::vector<double>& volumes)
    {
        if (flux.size() <= 0 || flux.size() != static_cast<Eigen::Index>(volumes.size()) || !flux.allFinite())
            throw std::invalid_argument("invalid Hydro divergence flux/volume layout");
        double result = 0.0;
        for (int i = 0; i < flux.size(); ++i)
        {
            if (!std::isfinite(volumes[i]) || volumes[i] <= 0.0)
                throw std::invalid_argument("invalid Hydro continuity cell volume");
            result = std::max(result, std::abs(flux[i] / volumes[i]));
        }
        if (!std::isfinite(result)) throw std::runtime_error("non-finite Hydro fluid divergence");
        return result;
    }

private:
    Eigen::VectorXd flux_, scale_;
    double maxFaceSpeed_ = 0.0;
};

// Original zero-time, complete-domain projection check. Its global norm and
// roundoff allowance are deliberately distinct from the fitted local gate.
// flux includes every physical boundary contribution, including prescribed
// velocities. The caller supplies its existing tolerance; no policy is chosen
// from a solver backend or dictionary here.
inline bool compatibleProjection(const Eigen::VectorXd& flux, const Eigen::VectorXd& rhs,
    const Eigen::VectorXd& absoluteFlux, double tolerance)
{
    if (flux.size() <= 0 || rhs.size() != flux.size() || absoluteFlux.size() != flux.size()
        || !std::isfinite(tolerance) || tolerance <= 0.0)
        throw std::invalid_argument("invalid Hydro projection continuity operands");
    const double rhsNorm = rhs.stableNorm();
    const double fluxNorm = flux.stableNorm();
    const double allowance = tolerance * rhsNorm
        + 512.0 * std::numeric_limits<double>::epsilon() * (rhsNorm + absoluteFlux.stableNorm());
    return std::isfinite(fluxNorm) && std::isfinite(allowance) && fluxNorm <= allowance;
}

// Original complete-leaf row strategy: visit each cell's real half-faces in
// side/part order, including boundary velocities. The length-weighted center
// speed is not the fitted domain's separate compact/native reconstruction.
// Optional output fields must be provided together. In the original parallel
// path nonfinite values are reduced and reported after workers join; no worker
// throws. This low-level diagnostic does not promise rollback of output fields.
template<class Mesh>
std::pair<double, double> nativeVelocity(Mesh& mesh,
    const VEdgeField<double>& vx, const HEdgeField<double>& vy, bool axisymmetric,
    int threads, CellField<double>* storedDivergence = nullptr, CellField<double>* storedMagnitude = nullptr)
{
    if ((storedDivergence == nullptr) != (storedMagnitude == nullptr) || threads <= 0)
        throw std::invalid_argument("invalid Hydro native diagnostic outputs/thread count");
    if (vx.size() != mesh.vEdgeIndexNum() || vy.size() != mesh.hEdgeIndexNum() || !vx.value || !vy.value)
        throw std::invalid_argument("Hydro native diagnostic velocity fields do not match the mesh");
    const bool storeFields = storedDivergence != nullptr;
    if (storeFields)
    {
        storedDivergence->initialize(mesh, 0.0);
        storedMagnitude->initialize(mesh, 0.0);
    }
    double maxDivergence = 0.0;
    double maxSpeed = 0.0;
    int invalidValues = 0;
    int invalidLabels = 0;
#pragma omp parallel for num_threads(threads) reduction(max: maxDivergence, maxSpeed) reduction(+: invalidValues, invalidLabels)
    for (int i = 0; i < mesh.cellIndexNum(); ++i)
    {
        const auto& cell = mesh.cells(mesh.cellIndexToLabel[i]);
        double flux = 0.0;
        double normalX[2] = {0.0, 0.0}, normalY[2] = {0.0, 0.0};
        for (int side = 0; side < 2; ++side)
        {
            for (int part = 0; part < 2; ++part)
            {
                const int vLabel = cell.cellToVEdge[side][part];
                if (vLabel != -1)
                {
                    const auto& edge = mesh.vEdges(vLabel);
                    const double velocity = vx.value[edge.index()];
                    const double area = edge.length() * (axisymmetric ? CartesianMeasure::twoPi * edge.midpointCoord().x : 1.0);
                    flux += (side == 0 ? -1.0 : 1.0) * velocity * area;
                    normalX[side] += velocity * edge.length() / cell.size().y;
                    maxSpeed = std::max(maxSpeed, std::abs(velocity));
                    invalidValues += !std::isfinite(velocity);
                    invalidLabels += vx.label(edge.index()) != vLabel;
                }
                const int hLabel = cell.cellToHEdge[side][part];
                if (hLabel != -1)
                {
                    const auto& edge = mesh.hEdges(hLabel);
                    const double velocity = vy.value[edge.index()];
                    const double area = edge.length() * (axisymmetric ? CartesianMeasure::twoPi * edge.midpointCoord().x : 1.0);
                    flux += (side == 0 ? -1.0 : 1.0) * velocity * area;
                    normalY[side] += velocity * edge.length() / cell.size().x;
                    maxSpeed = std::max(maxSpeed, std::abs(velocity));
                    invalidValues += !std::isfinite(velocity);
                    invalidLabels += vy.label(edge.index()) != hLabel;
                }
            }
        }
        const double volume = cell.area() * (axisymmetric ? CartesianMeasure::twoPi * cell.centerCoord().x : 1.0);
        const double divergence = flux / volume;
        const double speed = std::hypot(0.5 * (normalX[0] + normalX[1]), 0.5 * (normalY[0] + normalY[1]));
        invalidValues += !std::isfinite(divergence) || !std::isfinite(speed);
        maxDivergence = std::max(maxDivergence, std::abs(divergence));
        maxSpeed = std::max(maxSpeed, speed);
        if (storeFields)
        {
            storedDivergence->value[i] = divergence;
            storedMagnitude->value[i] = speed;
        }
    }
    if (invalidValues) throw std::runtime_error("velocity or flux diagnostics contain non-finite values");
    if (invalidLabels) throw std::invalid_argument("Hydro native diagnostic velocity fields have stale labels");
    return {maxDivergence, maxSpeed};
}
}
