#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <Eigen/QR>

#include "Core/Base/Vector2D.h"

// Shared implementation-only local trace reconstruction, not a global
// sharp-interface electric-field solver. It builds no global matrix, finds
// no interface, partitions no mesh face, and establishes neither conservation
// nor stability of a global flux stencil.
// Implementation-only: the anonymous ancestor gives every definition and
// type TU-local linkage, so different floating-point modes cannot share these
// engines through COMDAT. All includes stay outside these namespaces; no
// Detail type is a cross-TU public API.
namespace Quantum::Electricity::Detail
{
namespace
{
namespace DielectricInterfaceTrace
{
enum class Phase { Negative, Positive };

struct SamplePoint
{
    Vector2D point;
    Phase phase = Phase::Negative;
};

struct Result
{
    double interfacePotential = 0.0;
    double tangentialElectricField = 0.0;
    double normalDisplacement = 0.0;
};

namespace Detail
{
inline void require(bool condition, const char* reason)
{
    if (!condition) throw std::runtime_error(std::string("[DielectricInterfaceTraceCandidate] ") + reason);
}

inline double checked(long double value, const char* reason)
{
    const double result = static_cast<double>(value);
    require(std::isfinite(value) && std::isfinite(result), reason);
    return result;
}
}

// Given n directed negative -> positive and t=(-ny,nx), fit with equal weights:
//   V^s(p) = V_Gamma - Et * t.(p-x_Gamma) - Dn/epsilon_s * n.(p-x_Gamma).
// Thus E^s = Et*t + (Dn/epsilon_s)*n. The shared unknowns enforce local
// potential/tangential-field continuity and zero-free-charge normal-D
// continuity. They do NOT model an independently prescribed surface charge.
//
// Phase is supplied by the caller's actual geometry. Its sign is deliberately
// not inferred from the local tangent plane: a legitimate sample near a curved
// interface can be on the other side of that plane. A single sample exactly at
// x_Gamma is allowed; duplicate samples and a rank-deficient complete design
// are not. Both phase tags must occur.
class Reconstruction
{
public:
    Reconstruction(const Vector2D& interfacePoint, const Vector2D& unitNormal,
        double epsilonNegative, double epsilonPositive,
        const std::vector<SamplePoint>& samples, double lengthScale)
        : lengthScale_(lengthScale)
    {
        using Detail::require;
        require(samples.size() >= 3
            && samples.size() <= static_cast<std::size_t>(std::numeric_limits<Eigen::Index>::max()),
            "at least three samples within the supported index range are required");
        require(std::isfinite(interfacePoint.x) && std::isfinite(interfacePoint.y),
            "interface point must be finite");
        require(std::isfinite(unitNormal.x) && std::isfinite(unitNormal.y), "normal must be finite");
        const double normalLength = std::hypot(unitNormal.x, unitNormal.y);
        const double tolerance = 128.0 * std::numeric_limits<double>::epsilon();
        require(std::isfinite(normalLength) && std::abs(normalLength - 1.0) <= tolerance,
            "normal must have unit length");
        require(std::isfinite(lengthScale_) && lengthScale_ > 0.0,
            "local length scale must be finite and positive");
        require(std::isfinite(epsilonNegative) && epsilonNegative > 0.0
            && std::isfinite(epsilonPositive) && epsilonPositive > 0.0,
            "both permittivities must be finite and positive");

        // Harmonic epsilon is only a change of unknown scale, not a mixture
        // law or interface model: epsilon_s remains phase-specific in A.
        const double smaller = std::min(epsilonNegative, epsilonPositive);
        const double larger = std::max(epsilonNegative, epsilonPositive);
        permittivityScale_ = smaller * (2.0 / (1.0 + smaller / larger));
        require(std::isfinite(permittivityScale_) && permittivityScale_ > 0.0,
            "unrepresentable permittivity normalization");
        const double negativeRatio = permittivityScale_ / epsilonNegative;
        const double positiveRatio = permittivityScale_ / epsilonPositive;
        require(std::isfinite(negativeRatio) && negativeRatio > 0.0
            && std::isfinite(positiveRatio) && positiveRatio > 0.0,
            "permittivity contrast is not representable in the normalized design");

        const Eigen::Index count = static_cast<Eigen::Index>(samples.size());
        Eigen::MatrixXd design(count, 3);
        std::vector<std::array<double, 2>> coordinates;
        coordinates.reserve(samples.size());
        std::size_t negativeCount = 0, positiveCount = 0;
        for (std::size_t i = 0; i < samples.size(); ++i)
        {
            const auto& sample = samples[i];
            require(sample.phase == Phase::Negative || sample.phase == Phase::Positive, "invalid sample phase tag");
            negativeCount += sample.phase == Phase::Negative;
            positiveCount += sample.phase == Phase::Positive;
            require(std::isfinite(sample.point.x) && std::isfinite(sample.point.y), "sample point must be finite");
            const double x = Detail::checked((static_cast<long double>(sample.point.x) - interfacePoint.x) / lengthScale_,
                "unrepresentable normalized sample x coordinate");
            const double y = Detail::checked((static_cast<long double>(sample.point.y) - interfacePoint.y) / lengthScale_,
                "unrepresentable normalized sample y coordinate");
            for (const auto& other : coordinates)
            {
                const double separation = std::hypot(x - other[0], y - other[1]);
                const double scale = std::max({1.0, std::abs(x), std::abs(y), std::abs(other[0]), std::abs(other[1])});
                require(std::isfinite(separation) && separation > tolerance * scale,
                    "duplicate, coincident or roundoff-unresolved sample points");
            }
            coordinates.push_back({{x, y}});
            const double tangential = Detail::checked(-static_cast<long double>(unitNormal.y) * x
                + static_cast<long double>(unitNormal.x) * y, "non-finite tangential sample coordinate");
            const double normal = Detail::checked(static_cast<long double>(unitNormal.x) * x
                + static_cast<long double>(unitNormal.y) * y, "non-finite normal sample coordinate");
            design(static_cast<Eigen::Index>(i), 0) = 1.0;
            design(static_cast<Eigen::Index>(i), 1) = -tangential;
            design(static_cast<Eigen::Index>(i), 2) = -normal
                * (sample.phase == Phase::Negative ? negativeRatio : positiveRatio);
        }
        require(negativeCount > 0 && positiveCount > 0, "both phases must have at least one sample");
        require(design.allFinite(), "normalized design contains non-finite entries");

        // Dimensionless unknowns are [V_Gamma, ell*Et, ell*Dn/epsilon_ref].
        // Thin Q avoids an M-by-M identity/pseudoinverse construction: only
        // O(M) storage is needed for this three-unknown local fit.
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(design);
        qr.setThreshold(tolerance);
        rank_ = static_cast<int>(qr.rank());
        require(rank_ == 3, "rank-deficient two-phase trace design");
        const Eigen::MatrixXd qThin = qr.householderQ() * Eigen::MatrixXd::Identity(count, 3);
        const Eigen::Matrix3d upper = qr.matrixR().topLeftCorner<3, 3>().triangularView<Eigen::Upper>();
        const Eigen::MatrixXd unpermuted = upper.triangularView<Eigen::Upper>().solve(qThin.transpose());
        const Eigen::MatrixXd projection = qr.colsPermutation() * unpermuted;
        require(projection.allFinite(), "non-finite QR trace weights");
        condition_ = design.cwiseAbs().rowwise().sum().maxCoeff()
            * projection.cwiseAbs().rowwise().sum().maxCoeff();
        require(std::isfinite(condition_) && condition_ >= 1.0 && condition_ <= 1e8,
            "ill-conditioned normalized trace design");
        reproductionError_ = (projection * design - Eigen::Matrix3d::Identity()).cwiseAbs().maxCoeff();
        require(std::isfinite(reproductionError_)
            && reproductionError_ <= 256.0 * std::numeric_limits<double>::epsilon() * condition_,
            "QR weights fail the linear patch reproduction check");

        interfaceWeights_.resize(samples.size());
        tangentialWeights_.resize(samples.size());
        displacementWeights_.resize(samples.size());
        long double interfaceTail = 0.0L, tangentialTail = 0.0L, displacementTail = 0.0L;
        for (std::size_t i = 1; i < samples.size(); ++i)
        {
            interfaceWeights_[i] = projection(0, static_cast<Eigen::Index>(i));
            tangentialWeights_[i] = Detail::checked(static_cast<long double>(projection(1, static_cast<Eigen::Index>(i))) / lengthScale_,
                "non-finite physical tangential-field weights");
            displacementWeights_[i] = Detail::checked(static_cast<long double>(projection(2, static_cast<Eigen::Index>(i)))
                * permittivityScale_ / lengthScale_, "non-finite physical normal-displacement weights");
            interfaceTail += interfaceWeights_[i];
            tangentialTail += tangentialWeights_[i];
            displacementTail += displacementWeights_[i];
        }
        // Define the exposed raw linear map to be exactly the anchored map
        // used below (in exact arithmetic), including its constant response.
        interfaceWeights_[0] = Detail::checked(1.0L - interfaceTail, "non-finite anchored potential weight");
        tangentialWeights_[0] = Detail::checked(-tangentialTail, "non-finite anchored tangential-field weight");
        displacementWeights_[0] = Detail::checked(-displacementTail, "non-finite anchored normal-displacement weight");

        // Verify the map actually exposed to consumers as well: physical
        // scaling must not silently underflow an entire derivative functional.
        Eigen::MatrixXd anchored(3, count);
        for (std::size_t i = 0; i < samples.size(); ++i)
        {
            anchored(0, static_cast<Eigen::Index>(i)) = interfaceWeights_[i];
            anchored(1, static_cast<Eigen::Index>(i)) = Detail::checked(
                static_cast<long double>(tangentialWeights_[i]) * lengthScale_,
                "unrepresentable normalized tangential-field weights");
            anchored(2, static_cast<Eigen::Index>(i)) = Detail::checked(
                static_cast<long double>(displacementWeights_[i]) * lengthScale_ / permittivityScale_,
                "unrepresentable normalized normal-displacement weights");
        }
        const double anchoredError = (anchored * design - Eigen::Matrix3d::Identity()).cwiseAbs().maxCoeff();
        require(std::isfinite(anchoredError)
            && anchoredError <= 256.0 * std::numeric_limits<double>::epsilon() * condition_,
            "physical anchored weights fail the linear patch reproduction check");
        reproductionError_ = std::max(reproductionError_, anchoredError);
    }

    std::size_t sampleCount() const { return interfaceWeights_.size(); }
    int rank() const { return rank_; }
    double lengthScale() const { return lengthScale_; }
    double permittivityScale() const { return permittivityScale_; }
    // ||A||_inf ||A^+||_inf for the rectangular dimensionless design, not the
    // singular-value condition number or a global operator condition number.
    double conditionEstimate() const { return condition_; }
    // Maximum entry error of P*A-I, checking both the QR map and the exposed
    // anchored physical weights transformed back to dimensionless units.
    double normalizedReproductionError() const { return reproductionError_; }

    // Physical units: 1, 1/length, epsilon/length, respectively. Downstream
    // consumers must also remove the common potential offset when applying
    // these weights. A naive dot(weights,V) can lose accuracy for a large V
    // offset despite the algebraically correct constant-nullspace constraint.
    // No method can recover potential differences already lost when its input
    // doubles were rounded at a large common offset.
    const std::vector<double>& interfacePotentialWeights() const { return interfaceWeights_; }
    const std::vector<double>& tangentialElectricFieldWeights() const { return tangentialWeights_; }
    const std::vector<double>& normalDisplacementWeights() const { return displacementWeights_; }

    Result evaluate(const std::vector<double>& potentials) const
    {
        using Detail::require;
        require(potentials.size() == sampleCount(), "potential count must match the reconstruction samples");
        for (double value : potentials) require(std::isfinite(value), "sample potential must be finite");
        const long double reference = potentials[0];
        long double potentialCorrection = 0.0L, tangential = 0.0L, displacement = 0.0L;
        for (std::size_t i = 1; i < potentials.size(); ++i)
        {
            const long double difference = static_cast<long double>(potentials[i]) - reference;
            require(std::isfinite(difference), "unrepresentable potential difference");
            potentialCorrection += interfaceWeights_[i] * difference;
            tangential += tangentialWeights_[i] * difference;
            displacement += displacementWeights_[i] * difference;
        }
        Result result;
        result.interfacePotential = Detail::checked(reference + potentialCorrection, "non-finite reconstructed interface potential");
        result.tangentialElectricField = Detail::checked(tangential, "non-finite reconstructed tangential field");
        result.normalDisplacement = Detail::checked(displacement, "non-finite reconstructed normal displacement");
        return result;
    }

private:
    double lengthScale_, permittivityScale_ = 0.0, condition_ = 0.0, reproductionError_ = 0.0;
    int rank_ = 0;
    std::vector<double> interfaceWeights_, tangentialWeights_, displacementWeights_;
};
}
} // TU-local engine
} // namespace Quantum::Electricity::Detail
