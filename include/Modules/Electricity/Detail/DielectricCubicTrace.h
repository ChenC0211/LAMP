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

#include "Modules/Electricity/Detail/DielectricInterfaceTrace.h"
#include "Modules/Electricity/Detail/DielectricTraceFaceGeometry.h"
#include "Core/Mesh/CartesianMeasure.h"

// Shared implementation-only local, curvature-aware harmonic Taylor
// reconstruction, not a global dielectric solver. Geometry and phase
// labels are caller supplied.
// Implementation-only: the anonymous ancestor gives every definition and
// type TU-local linkage, so different floating-point modes cannot share these
// engines through COMDAT. All includes stay outside these namespaces; no
// Detail type is a cross-TU public API.
namespace Quantum::Electricity::Detail
{
namespace
{
namespace DielectricCubicTrace
{
using Phase = DielectricInterfaceTrace::Phase;
using SamplePoint = DielectricInterfaceTrace::SamplePoint;
using SegmentPiece = DielectricTraceFaceGeometry::SegmentPiece;
using LinkGeometry = DielectricTraceFaceGeometry::LinkGeometry;

struct State
{
    double interfacePotential = 0.0;
    double tangentialElectricField = 0.0;
    double normalDisplacement = 0.0;
    double tangentialFieldDerivative = 0.0;
    double normalDisplacementDerivative = 0.0;
    double tangentialFieldSecondDerivative = 0.0;
    double normalDisplacementSecondDerivative = 0.0;
};

// Known jump [Dn]=Dn_positive-Dn_negative and its arclength derivatives.
// With a jump, State's displacement entries always describe the negative side.
struct Jump
{
    double value = 0.0;
    double firstDerivative = 0.0;
    double secondDerivative = 0.0;
};

struct AffineCorrection
{
    std::vector<double> weights;
    double constant = 0.0;
};

namespace Detail
{
inline void require(bool condition, const char* reason)
{
    if (!condition) throw std::runtime_error(std::string("[DielectricCubicTraceCandidate] ") + reason);
}

inline double checked(long double value, const char* reason)
{
    const double result = static_cast<double>(value);
    require(std::isfinite(value) && std::isfinite(result), reason);
    return result;
}

inline void phase(Phase value)
{
    require(value == Phase::Negative || value == Phase::Positive, "invalid phase tag");
}

inline void unit(const Vector2D& value, const char* reason)
{
    const double length = std::hypot(value.x, value.y);
    require(std::isfinite(value.x) && std::isfinite(value.y) && std::isfinite(length)
        && std::abs(length - 1.0) <= 128.0 * std::numeric_limits<double>::epsilon(), reason);
}
}

// n points negative -> positive; t=(-ny,nx); kappa=div(n), ks=d_s kappa.
// Along the interface, t_s=-kappa*n and n_s=kappa*t. Each phase has
// Htt=-dsEt-kappa*Dn/epsilon, Htn=kappa*Et-dsDn/epsilon, Hnn=-Htt.
// Differentiating these moving-frame contractions gives Cartesian derivatives:
// Tttt=-dssEt+2*kappa^2*Et-ks*Dn/epsilon-3*kappa*dsDn/epsilon,
// Tttn=ks*Et+3*kappa*dsEt+2*kappa^2*Dn/epsilon-dssDn/epsilon,
// Ttnn=-Tttt, Tnnn=-Tttn. Thus both differentiated Laplacian constraints hold.
// The cubic Taylor contribution is
// Tttt*(s^3-3*s*z^2)/6 + Tttn*(3*s^2*z-z^3)/6.
// These are local transmission jets, NOT exact continuity on a finite curved
// interface away from xGamma. In particular, ks must not silently be set to zero.
// Optional RZ uses x=r>0 at the anchor and MERIDIONAL kappa,ks. With
// Vr=t_r*Vt+n_r*Vn, Hnn=-Htt-Vr/r; Tttt is unchanged, while Tttn gains
// -kappa*Vr/r. Differentiating the RZ PDE gives
// Ttnn=-Tttt-Htr/r+Vr*t_r/r^2 and Tnnn=-Tttn-Hnr/r+Vr*n_r/r^2.
// Thus the PDE and its first derivatives hold at the anchor. This is not a
// claim that a general cubic solves the variable-coefficient PDE everywhere.
//
// Equal-weight QR fits [V,h*Et,h*Dn/epsRef,h^2*dsEt,h^2*dsDn/epsRef,
// h^3*dssEt,h^3*dssDn/epsRef]. Actual phase labels are caller supplied,
// including tangent-plane and interface samples; both phases must occur.
// This object owns an immutable design/weight snapshot. No shape, global
// electric solver, flux assembly, or Hydro dependency is present.
// Finite/representable normalized coordinates, coefficients and physical
// functionals are required; this is floating-point QR, not exact arithmetic.
class Reconstruction
{
public:
    Reconstruction(const Vector2D& interfacePoint, const Vector2D& unitNormal,
        double curvature, double curvatureTangentialDerivative,
        double epsilonNegative, double epsilonPositive,
        const std::vector<SamplePoint>& samples, double lengthScale, bool axisymmetric = false)
        : interfacePoint_(interfacePoint), normal_(unitNormal), lengthScale_(lengthScale),
          epsilonNegative_(epsilonNegative), epsilonPositive_(epsilonPositive), axisymmetric_(axisymmetric)
    {
        using Detail::require;
        require(samples.size() >= 7
            && samples.size() <= static_cast<std::size_t>(std::numeric_limits<Eigen::Index>::max()),
            "at least seven samples within the supported index range are required");
        require(std::isfinite(interfacePoint.x) && std::isfinite(interfacePoint.y), "interface point must be finite");
        require(!axisymmetric_ || interfacePoint.x > 0.0, "RZ trace anchor must have strictly positive radius");
        Detail::unit(unitNormal, "interface normal must have finite unit length");
        require(std::isfinite(curvature), "curvature must be finite");
        require(std::isfinite(curvatureTangentialDerivative), "tangential curvature derivative must be finite");
        require(std::isfinite(lengthScale_) && lengthScale_ > 0.0, "length scale must be finite and positive");
        require(std::isfinite(epsilonNegative_) && epsilonNegative_ > 0.0
            && std::isfinite(epsilonPositive_) && epsilonPositive_ > 0.0,
            "both permittivities must be finite and positive");
        scaledCurvature_ = Detail::checked(static_cast<long double>(curvature) * lengthScale_,
            "unrepresentable normalized curvature");
        scaledCurvatureDerivative_ = Detail::checked(static_cast<long double>(curvatureTangentialDerivative)
            * lengthScale_ * lengthScale_, "unrepresentable normalized curvature derivative");
        require((curvature == 0.0 || scaledCurvature_ != 0.0)
            && (curvatureTangentialDerivative == 0.0 || scaledCurvatureDerivative_ != 0.0),
            "nonzero normalized curvature data underflow");
        const double smaller = std::min(epsilonNegative_, epsilonPositive_);
        const double larger = std::max(epsilonNegative_, epsilonPositive_);
        permittivityScale_ = smaller * (2.0 / (1.0 + smaller / larger));
        require(std::isfinite(permittivityScale_) && permittivityScale_ > 0.0,
            "unrepresentable permittivity normalization");
        negativeRatio_ = permittivityScale_ / epsilonNegative_;
        positiveRatio_ = permittivityScale_ / epsilonPositive_;
        require(std::isfinite(negativeRatio_) && negativeRatio_ > 0.0
            && std::isfinite(positiveRatio_) && positiveRatio_ > 0.0,
            "permittivity contrast is not representable in the normalized design");
        if (axisymmetric_) prepareAxisymmetricJets();

        const double tolerance = 128.0 * std::numeric_limits<double>::epsilon();
        const Eigen::Index count = static_cast<Eigen::Index>(samples.size());
        Eigen::MatrixXd design(count, 7);
        std::vector<std::array<double, 2>> coordinates;
        coordinates.reserve(samples.size());
        sampleJumpBasis_.reserve(samples.size());
        std::size_t negativeCount = 0, positiveCount = 0;
        for (std::size_t i = 0; i < samples.size(); ++i)
        {
            Detail::phase(samples[i].phase);
            negativeCount += samples[i].phase == Phase::Negative;
            positiveCount += samples[i].phase == Phase::Positive;
            const auto coordinate = normalized(samples[i].point);
            for (const auto& other : coordinates)
            {
                const double separation = std::hypot(coordinate[0] - other[0], coordinate[1] - other[1]);
                const double scale = std::max({1.0, std::abs(coordinate[0]), std::abs(coordinate[1]),
                    std::abs(other[0]), std::abs(other[1])});
                require(std::isfinite(separation) && separation > tolerance * scale,
                    "duplicate, coincident or roundoff-unresolved sample points");
            }
            coordinates.push_back(coordinate);
            const Basis values = potentialBasis(coordinate, samples[i].phase);
            for (int j = 0; j < 7; ++j) design(static_cast<Eigen::Index>(i), j) = values[j];
            sampleJumpBasis_.push_back(samples[i].phase == Phase::Positive
                ? std::array<double, 3>{{values[2], values[4], values[6]}}
                : std::array<double, 3>{{0.0, 0.0, 0.0}});
        }
        require(negativeCount > 0 && positiveCount > 0, "both phases must have at least one sample");
        require(design.allFinite(), "normalized design contains non-finite entries");

        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(design);
        qr.setThreshold(tolerance);
        rank_ = static_cast<int>(qr.rank());
        require(rank_ == 7, "rank-deficient two-phase cubic trace design");
        const Eigen::MatrixXd qThin = qr.householderQ() * Eigen::MatrixXd::Identity(count, 7);
        const Eigen::Matrix<double, 7, 7> upper = qr.matrixR().topLeftCorner<7, 7>().triangularView<Eigen::Upper>();
        const Eigen::MatrixXd unpermuted = upper.triangularView<Eigen::Upper>().solve(qThin.transpose());
        projection_ = qr.colsPermutation() * unpermuted;
        require(projection_.allFinite(), "non-finite QR trace weights");
        condition_ = design.cwiseAbs().rowwise().sum().maxCoeff()
            * projection_.cwiseAbs().rowwise().sum().maxCoeff();
        require(std::isfinite(condition_) && condition_ >= 1.0 && condition_ <= 1e8,
            "ill-conditioned normalized cubic trace design");
        reproductionError_ = (projection_ * design - Eigen::Matrix<double, 7, 7>::Identity()).cwiseAbs().maxCoeff();
        require(std::isfinite(reproductionError_)
            && reproductionError_ <= 256.0 * std::numeric_limits<double>::epsilon() * condition_,
            "QR weights fail the cubic patch reproduction check");

        std::array<double, 7> physicalScales{};
        for (int j = 0; j < 7; ++j)
        {
            Basis functional{};
            long double scale = j == 2 || j == 4 || j == 6 ? permittivityScale_ : 1.0L;
            if (j != 0) scale /= lengthScale_;
            if (j >= 3) scale /= lengthScale_;
            if (j >= 5) scale /= lengthScale_;
            functional[j] = Detail::checked(scale, "unrepresentable physical trace functional");
            require(functional[j] > 0.0, "physical trace functional scale underflow");
            physicalScales[j] = functional[j];
            stateWeights_[j] = project(functional);
        }
        Eigen::MatrixXd anchored(7, count);
        for (int j = 0; j < 7; ++j)
            for (std::size_t i = 0; i < samples.size(); ++i)
            {
                const long double value = static_cast<long double>(stateWeights_[j][i]) / physicalScales[j];
                anchored(j, static_cast<Eigen::Index>(i)) = Detail::checked(value,
                    "unrepresentable normalized physical trace weights");
            }
        const double anchoredError = (anchored * design - Eigen::Matrix<double, 7, 7>::Identity()).cwiseAbs().maxCoeff();
        require(std::isfinite(anchoredError)
            && anchoredError <= 256.0 * std::numeric_limits<double>::epsilon() * condition_,
            "physical anchored weights fail the cubic patch reproduction check");
        reproductionError_ = std::max(reproductionError_, anchoredError);
    }

    std::size_t sampleCount() const { return static_cast<std::size_t>(projection_.cols()); }
    int rank() const { return rank_; }
    double lengthScale() const { return lengthScale_; }
    double permittivityScale() const { return permittivityScale_; }
    // Rectangular infinity-norm design/left-inverse condition, not cond_2.
    double conditionEstimate() const { return condition_; }
    double normalizedReproductionError() const { return reproductionError_; }

    State evaluate(const std::vector<double>& potentials) const
    {
        Detail::require(potentials.size() == sampleCount(), "potential count must match the reconstruction samples");
        for (double value : potentials) Detail::require(std::isfinite(value), "sample potential must be finite");
        std::array<long double, 7> values{{static_cast<long double>(potentials[0]), 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L}};
        for (std::size_t i = 1; i < potentials.size(); ++i)
        {
            const long double difference = static_cast<long double>(potentials[i]) - potentials[0];
            Detail::require(std::isfinite(difference), "unrepresentable potential difference");
            for (int j = 0; j < 7; ++j) values[j] += stateWeights_[j][i] * difference;
        }
        State result;
        result.interfacePotential = Detail::checked(values[0], "non-finite reconstructed interface potential");
        result.tangentialElectricField = Detail::checked(values[1], "non-finite reconstructed tangential field");
        result.normalDisplacement = Detail::checked(values[2], "non-finite reconstructed normal displacement");
        result.tangentialFieldDerivative = Detail::checked(values[3], "non-finite reconstructed tangential-field derivative");
        result.normalDisplacementDerivative = Detail::checked(values[4], "non-finite reconstructed displacement derivative");
        result.tangentialFieldSecondDerivative = Detail::checked(values[5], "non-finite reconstructed second tangential-field derivative");
        result.normalDisplacementSecondDerivative = Detail::checked(values[6], "non-finite reconstructed second displacement derivative");
        return result;
    }

    // The positive-side particular polynomial uses the very same Dn columns
    // as the homogeneous design. Subtract its sampled values before applying
    // the existing QR map: no new unknowns, fit or factorization are needed.
    // Difference V and the particular separately to retain gauge cancellation.
    State evaluate(const std::vector<double>& potentials, const Jump& jump) const
    {
        const auto coefficients = normalizedJump(jump);
        if (zero(jump)) return evaluate(potentials);
        Detail::require(potentials.size() == sampleCount(), "potential count must match the reconstruction samples");
        for (double value : potentials) Detail::require(std::isfinite(value), "sample potential must be finite");
        const auto particular = sampleParticular(coefficients);
        std::array<long double, 7> values{{static_cast<long double>(potentials[0]) - particular[0],
            0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L}};
        for (std::size_t i = 1; i < potentials.size(); ++i)
        {
            const long double difference = (static_cast<long double>(potentials[i]) - potentials[0])
                - (static_cast<long double>(particular[i]) - particular[0]);
            Detail::require(std::isfinite(difference), "unrepresentable affine potential difference");
            for (int j = 0; j < 7; ++j) values[j] += stateWeights_[j][i] * difference;
        }
        State result;
        result.interfacePotential = Detail::checked(values[0], "non-finite reconstructed interface potential");
        result.tangentialElectricField = Detail::checked(values[1], "non-finite reconstructed tangential field");
        result.normalDisplacement = Detail::checked(values[2], "non-finite reconstructed normal displacement");
        result.tangentialFieldDerivative = Detail::checked(values[3], "non-finite reconstructed tangential-field derivative");
        result.normalDisplacementDerivative = Detail::checked(values[4], "non-finite reconstructed displacement derivative");
        result.tangentialFieldSecondDerivative = Detail::checked(values[5], "non-finite reconstructed second tangential-field derivative");
        result.normalDisplacementSecondDerivative = Detail::checked(values[6], "non-finite reconstructed second displacement derivative");
        return result;
    }

    double particularPotential(const Vector2D& point, Phase phase, const Jump& jump) const
    {
        const auto coefficients = normalizedJump(jump);
        return particular(potentialBasis(normalized(point), phase), phase, coefficients);
    }

    double particularDisplacementNormal(const Vector2D& point, Phase phase,
        const Vector2D& unitDirection, const Jump& jump) const
    {
        const auto coefficients = normalizedJump(jump);
        Detail::unit(unitDirection, "displacement direction must have finite unit length");
        return particular(displacementBasis(normalized(point), phase, unitDirection), phase, coefficients);
    }

    // Weights use the input sample order. Evaluate potential functionals as
    // V0 + sum_{j>0} w_j*(Vj-V0), and all field/flux functionals as the sum
    // without V0. Raw dot products can lose accuracy at large gauge offsets.
    std::vector<double> potentialWeights(const Vector2D& point, Phase phase) const
    {
        return project(potentialBasis(normalized(point), phase));
    }

    std::vector<double> displacementNormalWeights(const Vector2D &point, Phase phase,
                                                  const Vector2D &unitDirection) const
    {
        Detail::unit(unitDirection, "displacement direction must have finite unit length");
        return project(displacementBasis(normalized(point), phase, unitDirection));
    }

    std::vector<double> correctionWeights(const Vector2D &owner, Phase ownerPhase, const Vector2D &neighbor,
                                          Phase neighborPhase, const Vector2D &ownerFaceNormal,
                                          double faceLength, double linkDistance, double twoPointPermittivity,
                                          const std::vector<SegmentPiece> &pieces,
                                          LinkGeometry geometry = LinkGeometry::Orthogonal) const;

    AffineCorrection affineCorrection(const Vector2D &owner, Phase ownerPhase, const Vector2D &neighbor,
                                      Phase neighborPhase, const Vector2D &ownerFaceNormal, double faceLength,
                                      double linkDistance, double twoPointPermittivity,
                                      const std::vector<SegmentPiece> &pieces, const Jump &jump,
                                      LinkGeometry geometry = LinkGeometry::Orthogonal) const;

  private:
    using Basis = std::array<double, 7>;

    struct PolynomialJet
    {
        long double vs = 0.0L, vn = 0.0L, htt = 0.0L, htn = 0.0L, hnn = 0.0L;
        long double ttt = 0.0L, ttn = 0.0L, tnn = 0.0L, nnn = 0.0L;
    };

    void prepareAxisymmetricJets()
    {
        const long double inverseRadius = static_cast<long double>(lengthScale_) / interfacePoint_.x;
        const long double inverseSquared = inverseRadius * inverseRadius;
        Detail::require(std::isfinite(inverseRadius) && inverseRadius > 0.0L
            && std::isfinite(inverseSquared) && inverseSquared > 0.0L,
            "unrepresentable normalized RZ inverse radius");
        const long double tr = -static_cast<long double>(normal_.y), nr = normal_.x;
        const long double c = scaledCurvature_, d = scaledCurvatureDerivative_;
        for (int j = 1; j < 7; ++j)
        {
            const long double e = j == 1, dn = j == 2, es = j == 3, dns = j == 4;
            const long double ess = j == 5, dnss = j == 6;
            auto& a = axisymmetricJets_[j];
            a.vs = -e; a.vn = -dn;
            a.htt = -es - c * dn;
            a.htn = c * e - dns;
            const long double vr = tr * a.vs + nr * a.vn;
            a.hnn = -a.htt - inverseRadius * vr;
            a.ttt = -ess + 2.0L * c * c * e - d * dn - 3.0L * c * dns;
            a.ttn = d * e + 3.0L * c * es + 2.0L * c * c * dn - dnss - c * inverseRadius * vr;
            const long double htr = tr * a.htt + nr * a.htn;
            const long double hnr = tr * a.htn + nr * a.hnn;
            a.tnn = -a.ttt - inverseRadius * htr + inverseSquared * vr * tr;
            a.nnn = -a.ttn - inverseRadius * hnr + inverseSquared * vr * nr;
            for (long double value : {a.vs, a.vn, a.htt, a.htn, a.hnn, a.ttt, a.ttn, a.tnn, a.nnn})
                Detail::require(std::isfinite(value), "non-finite RZ Taylor coefficient");
        }
    }

    Basis axisymmetricPotentialBasis(const std::array<double, 2>& point, Phase phase) const
    {
        const auto coordinate = local(point);
        const long double s = coordinate[0], z = coordinate[1];
        const long double ratio = phase == Phase::Negative ? negativeRatio_ : positiveRatio_;
        Basis result{}; result[0] = 1.0;
        for (int j = 1; j < 7; ++j)
        {
            const auto& a = axisymmetricJets_[j];
            const long double value = a.vs * s + a.vn * z
                + (a.htt * s * s + 2.0L * a.htn * s * z + a.hnn * z * z) / 2.0L
                + (a.ttt * s * s * s + 3.0L * a.ttn * s * s * z
                    + 3.0L * a.tnn * s * z * z + a.nnn * z * z * z) / 6.0L;
            result[j] = Detail::checked((j == 2 || j == 4 || j == 6 ? ratio : 1.0L) * value,
                "non-finite RZ potential basis");
        }
        return result;
    }

    Basis axisymmetricDisplacementBasis(const std::array<double, 2>& point, Phase phase,
        const Vector2D& direction) const
    {
        const auto coordinate = local(point);
        const long double s = coordinate[0], z = coordinate[1];
        const long double et = -static_cast<long double>(normal_.y) * direction.x
            + static_cast<long double>(normal_.x) * direction.y;
        const long double en = static_cast<long double>(normal_.x) * direction.x
            + static_cast<long double>(normal_.y) * direction.y;
        const long double epsilon = phase == Phase::Negative ? epsilonNegative_ : epsilonPositive_;
        Basis result{};
        for (int j = 1; j < 7; ++j)
        {
            const auto& a = axisymmetricJets_[j];
            const long double ds = a.vs + a.htt * s + a.htn * z
                + (a.ttt * s * s + 2.0L * a.ttn * s * z + a.tnn * z * z) / 2.0L;
            const long double dz = a.vn + a.htn * s + a.hnn * z
                + (a.ttn * s * s + 2.0L * a.tnn * s * z + a.nnn * z * z) / 2.0L;
            const long double scale = j == 2 || j == 4 || j == 6 ? permittivityScale_ : epsilon;
            result[j] = Detail::checked(-scale * (ds * et + dz * en) / lengthScale_,
                "non-finite RZ displacement functional");
        }
        return result;
    }

    double radialWeightAt(const std::array<double, 2>& point) const
    {
        const double radius = Detail::checked(static_cast<long double>(interfacePoint_.x)
            + static_cast<long double>(lengthScale_) * point[0], "unrepresentable RZ quadrature radius");
        return CartesianMeasure::radialWeight(radius);
    }

    double axisymmetricFaceMeasure(double faceLength, const std::vector<SegmentPiece>& pieces) const
    {
        const double radius = Detail::checked(0.5L * pieces.front().begin.x + 0.5L * pieces.back().end.x,
            "unrepresentable RZ face midpoint radius");
        return CartesianMeasure::face(faceLength, radius, true);
    }

    static bool zero(const Jump& jump)
    { return jump.value == 0.0 && jump.firstDerivative == 0.0 && jump.secondDerivative == 0.0; }

    std::array<double, 3> normalizedJump(const Jump& jump) const
    {
        const std::array<double, 3> input{{jump.value, jump.firstDerivative, jump.secondDerivative}};
        std::array<double, 3> result{};
        for (int j = 0; j < 3; ++j)
        {
            Detail::require(std::isfinite(input[j]), "known displacement jump and derivatives must be finite");
            long double value = static_cast<long double>(input[j]) / permittivityScale_;
            for (int k = 0; k <= j; ++k) value *= lengthScale_;
            result[j] = Detail::checked(value, "unrepresentable normalized displacement jump");
            Detail::require(input[j] == 0.0 || result[j] != 0.0, "nonzero normalized displacement jump underflow");
        }
        return result;
    }

    static double particular(const Basis& basis, Phase phase, const std::array<double, 3>& coefficients)
    {
        if (phase == Phase::Negative) return 0.0;
        long double value = 0.0L;
        for (int j = 0; j < 3; ++j) value += static_cast<long double>(basis[2 + 2 * j]) * coefficients[j];
        return Detail::checked(value, "non-finite known-jump particular functional");
    }

    std::vector<double> sampleParticular(const std::array<double, 3>& coefficients) const
    {
        std::vector<double> result(sampleCount());
        for (std::size_t i = 0; i < result.size(); ++i)
        {
            long double value = 0.0L;
            for (int j = 0; j < 3; ++j) value += static_cast<long double>(sampleJumpBasis_[i][j]) * coefficients[j];
            result[i] = Detail::checked(value, "non-finite sampled known-jump particular potential");
        }
        return result;
    }

    std::array<double, 2> normalized(const Vector2D& point) const
    {
        Detail::require(std::isfinite(point.x) && std::isfinite(point.y), "point must be finite");
        Detail::require(!axisymmetric_ || point.x >= 0.0, "RZ sample/query radius must be nonnegative");
        return {{Detail::checked((static_cast<long double>(point.x) - interfacePoint_.x) / lengthScale_,
            "unrepresentable normalized x coordinate"),
            Detail::checked((static_cast<long double>(point.y) - interfacePoint_.y) / lengthScale_,
                "unrepresentable normalized y coordinate")}};
    }

    std::array<double, 2> local(const std::array<double, 2>& point) const
    {
        return {{Detail::checked(-static_cast<long double>(normal_.y) * point[0] + static_cast<long double>(normal_.x) * point[1],
            "unrepresentable local tangential coordinate"),
            Detail::checked(static_cast<long double>(normal_.x) * point[0] + static_cast<long double>(normal_.y) * point[1],
                "unrepresentable local normal coordinate")}};
    }

    Basis potentialBasis(const std::array<double, 2>& point, Phase phase) const
    {
        Detail::phase(phase);
        if (axisymmetric_) return axisymmetricPotentialBasis(point, phase);
        const auto coordinate = local(point);
        const long double s = coordinate[0], z = coordinate[1], c = scaledCurvature_;
        const long double d = scaledCurvatureDerivative_;
        const long double ratio = phase == Phase::Negative ? negativeRatio_ : positiveRatio_;
        const long double q = 0.5L * (z * z - s * s);
        const long double a = (s * s * s - 3.0L * s * z * z) / 6.0L;
        const long double b = (3.0L * s * s * z - z * z * z) / 6.0L;
        return {{1.0,
            Detail::checked(-s + c * s * z + 2.0L * c * c * a + d * b,
                "non-finite cubic tangential basis"),
            Detail::checked(ratio * (-z + c * q - d * a + 2.0L * c * c * b),
                "non-finite cubic displacement basis"),
            Detail::checked(q + 3.0L * c * b, "non-finite tangential derivative basis"),
            Detail::checked(ratio * (-s * z - 3.0L * c * a), "non-finite displacement derivative basis"),
            Detail::checked(-a, "non-finite second tangential derivative basis"),
            Detail::checked(-ratio * b, "non-finite second displacement derivative basis")}};
    }

    Basis displacementBasis(const std::array<double, 2>& point, Phase phase, const Vector2D& direction) const
    {
        Detail::phase(phase);
        if (axisymmetric_) return axisymmetricDisplacementBasis(point, phase, direction);
        const auto coordinate = local(point);
        const long double s = coordinate[0], z = coordinate[1], c = scaledCurvature_;
        const long double d = scaledCurvatureDerivative_;
        const long double et = -static_cast<long double>(normal_.y) * direction.x
            + static_cast<long double>(normal_.x) * direction.y;
        const long double en = static_cast<long double>(normal_.x) * direction.x
            + static_cast<long double>(normal_.y) * direction.y;
        const long double epsilon = phase == Phase::Negative ? epsilonNegative_ : epsilonPositive_;
        const long double q = 0.5L * (z * z - s * s);
        // a_s=-q, a_z=-s*z, b_s=s*z, b_z=-q for the two harmonic cubics.
        const std::array<long double, 7> ds{{0.0L,
            -1.0L + c * z - 2.0L * c * c * q + d * s * z,
            -c * s + d * q + 2.0L * c * c * s * z,
            -s + 3.0L * c * s * z,
            -z + 3.0L * c * q,
            q, -s * z}};
        const std::array<long double, 7> dz{{0.0L,
            c * s - 2.0L * c * c * s * z - d * q,
            -1.0L + c * z + d * s * z - 2.0L * c * c * q,
            z - 3.0L * c * q,
            -s + 3.0L * c * s * z,
            s * z, q}};
        Basis result{};
        for (int j = 1; j < 7; ++j)
        {
            // epsilon*(epsRef/epsilon) is algebraically epsRef. Keep that
            // cancellation before arithmetic, as in P2, to avoid losing a
            // displacement functional through an intermediate tiny ratio.
            const long double scale = j == 2 || j == 4 || j == 6 ? permittivityScale_ : epsilon;
            result[j] = Detail::checked(-scale * (ds[j] * et + dz[j] * en) / lengthScale_,
                "non-finite cubic displacement functional");
        }
        return result;
    }

    std::vector<double> project(const Basis& functional) const
    {
        std::vector<double> weights(sampleCount(), 0.0);
        long double tail = 0.0L;
        bool nonzeroTail = false;
        for (std::size_t i = 1; i < weights.size(); ++i)
        {
            long double value = 0.0L;
            for (int j = 0; j < 7; ++j)
                value += static_cast<long double>(functional[j]) * projection_(j, static_cast<Eigen::Index>(i));
            weights[i] = Detail::checked(value, "unrepresentable projected functional weight");
            nonzeroTail = nonzeroTail || weights[i] != 0.0;
            tail += weights[i];
        }
        const bool nonzeroFunctional = std::any_of(functional.begin(), functional.end(),
            [](double value) { return value != 0.0; });
        Detail::require(functional[0] != 0.0 || !nonzeroFunctional || nonzeroTail,
            "nonzero constant-null functional weights underflow");
        weights[0] = Detail::checked(static_cast<long double>(functional[0]) - tail, "unrepresentable anchored functional weight");
        return weights;
    }

    Vector2D interfacePoint_, normal_;
    double lengthScale_, epsilonNegative_, epsilonPositive_;
    bool axisymmetric_ = false;
    std::array<PolynomialJet, 7> axisymmetricJets_{};
    double permittivityScale_ = 0.0, negativeRatio_ = 0.0, positiveRatio_ = 0.0, scaledCurvature_ = 0.0;
    double scaledCurvatureDerivative_ = 0.0;
    double condition_ = 0.0, reproductionError_ = 0.0;
    int rank_ = 0;
    Eigen::MatrixXd projection_;
    std::array<std::vector<double>, 7> stateWeights_;
    std::vector<std::array<double, 3>> sampleJumpBasis_;
};

// ADDITIVE correction to twoPointPermittivity*L/d*(Vowner-Vneighbor), with
// supplied positive two-point permittivity, not geometric curvature. On each
// actual-phase straight piece, two-point Gauss integrates the quadratic
// displacement of the SAME local cubic fit (exact for that polynomial in
// real arithmetic). Neither frame nor curvature is re-evaluated at the Gauss
// points. Subtract in the seven-dimensional basis and project only once;
// the constant column stays exactly zero before projection/anchoring.
// A global caller assembles each physical face once and uses its exact
// negative for the other owner; independently repeated fits are not a proof
// of global conservation. No boundary-potential constant is added here.
inline std::vector<double> Reconstruction::correctionWeights(
    const Vector2D &owner, Phase ownerPhase, const Vector2D &neighbor, Phase neighborPhase,
    const Vector2D &ownerFaceNormal, double faceLength, double linkDistance, double twoPointPermittivity,
    const std::vector<SegmentPiece> &pieces, LinkGeometry geometry) const
{
    using Detail::require;
    std::array<long double, 7> integral{};
    const auto link = DielectricTraceFaceGeometry::visitValidatedFace(
        owner, ownerPhase, neighbor, neighborPhase, ownerFaceNormal, faceLength, linkDistance,
        twoPointPermittivity, pieces, lengthScale_,
        [this](const Vector2D &point) { return normalized(point); },
        [&](const std::array<double, 2> &begin, const std::array<double, 2> &end,
            const std::array<double, 2> &midpoint, double length, Phase phase,
            const std::array<long double, 2> &span) {
            const long double offsetScale = 0.5L / std::sqrt(3.0L);
            std::array<double, 2> minus{}, plus{};
            for (int axis = 0; axis < 2; ++axis)
            {
                const long double offset = span[axis] * offsetScale;
                minus[axis] = Detail::checked(static_cast<long double>(midpoint[axis]) - offset,
                                              "non-finite normalized negative Gauss point");
                plus[axis] = Detail::checked(static_cast<long double>(midpoint[axis]) + offset,
                                             "non-finite normalized positive Gauss point");
            }
            const Basis displacementMinus = displacementBasis(minus, phase, ownerFaceNormal);
            const Basis displacementPlus = displacementBasis(plus, phase, ownerFaceNormal);
            const long double physicalLength = static_cast<long double>(length) * lengthScale_;
            require(std::isfinite(physicalLength) && physicalLength > 0.0L,
                    "unrepresentable physical piece length");
            const long double halfLength = 0.5L * physicalLength;
            require(halfLength > 0.0L, "physical Gauss weight underflow");
            if (axisymmetric_)
            {
                const double radialMinus = radialWeightAt(minus), radialPlus = radialWeightAt(plus);
                const long double weightMinus = halfLength * radialMinus,
                                  weightPlus = halfLength * radialPlus;
                require(std::isfinite(weightMinus) && std::isfinite(weightPlus) &&
                            (radialMinus == 0.0 || weightMinus > 0.0L) &&
                            (radialPlus == 0.0 || weightPlus > 0.0L),
                        "RZ Gauss weight overflow or underflow");
                for (int j = 1; j < 7; ++j)
                    integral[j] += weightMinus * displacementMinus[j] + weightPlus * displacementPlus[j];
            }
            else
                for (int j = 1; j < 7; ++j)
                    integral[j] +=
                        halfLength * (static_cast<long double>(displacementMinus[j]) + displacementPlus[j]);
        },
        "[DielectricCubicTraceCandidate] ", geometry);
    const Basis ownerValues = potentialBasis(link.owner, ownerPhase);
    const Basis neighborValues = potentialBasis(link.neighbor, neighborPhase);
    long double conductance = static_cast<long double>(twoPointPermittivity) * faceLength / linkDistance;
    double faceMeasure = faceLength;
    if (axisymmetric_)
    {
        faceMeasure = axisymmetricFaceMeasure(faceLength, pieces);
        conductance = static_cast<long double>(twoPointPermittivity) * faceMeasure / linkDistance;
    }
    require(std::isfinite(conductance) &&
                (axisymmetric_ && faceMeasure == 0.0 ? conductance == 0.0L : conductance > 0.0L),
            "unrepresentable two-point conductance");
    Basis correction{};
    for (int j = 1; j < 7; ++j)
        correction[j] = Detail::checked(
            integral[j] - conductance * (static_cast<long double>(ownerValues[j]) - neighborValues[j]),
            "non-finite integrated cubic flux correction");
    return project(correction);
}

// For actual potentials V, this returns Q(Vfit)-TP(Vfit)=weights*V+constant.
// The weights are the unchanged charge-free map; the constant removes its
// response to sampled particular potentials and adds the physical particular
// flux. This is ONLY a face correction. A global charged Poisson equation
// also needs its cell-integrated interface-charge source; do not omit it.
inline AffineCorrection Reconstruction::affineCorrection(const Vector2D &owner, Phase ownerPhase,
                                                         const Vector2D &neighbor, Phase neighborPhase,
                                                         const Vector2D &ownerFaceNormal, double faceLength,
                                                         double linkDistance, double twoPointPermittivity,
                                                         const std::vector<SegmentPiece> &pieces,
                                                         const Jump &jump, LinkGeometry geometry) const
{
    const auto coefficients = normalizedJump(jump);
    AffineCorrection result{correctionWeights(owner, ownerPhase, neighbor, neighborPhase, ownerFaceNormal,
                                              faceLength, linkDistance, twoPointPermittivity, pieces,
                                              geometry),
                            0.0};
    if (zero(jump))
        return result;
    long double integral = 0.0L;
    const auto link = DielectricTraceFaceGeometry::visitValidatedFace(
        owner, ownerPhase, neighbor, neighborPhase, ownerFaceNormal, faceLength, linkDistance,
        twoPointPermittivity, pieces, lengthScale_,
        [this](const Vector2D &point) { return normalized(point); },
        [&](const std::array<double, 2> &begin, const std::array<double, 2> &end,
            const std::array<double, 2> &midpoint, double length, Phase phase,
            const std::array<long double, 2> &span) {
            const long double offsetScale = 0.5L / std::sqrt(3.0L);
            std::array<double, 2> minus{}, plus{};
            for (int axis = 0; axis < 2; ++axis)
            {
                const long double offset = span[axis] * offsetScale;
                minus[axis] = Detail::checked(static_cast<long double>(midpoint[axis]) - offset,
                                              "non-finite normalized negative Gauss point");
                plus[axis] = Detail::checked(static_cast<long double>(midpoint[axis]) + offset,
                                             "non-finite normalized positive Gauss point");
            }
            const double a =
                particular(displacementBasis(minus, phase, ownerFaceNormal), phase, coefficients);
            const double b = particular(displacementBasis(plus, phase, ownerFaceNormal), phase, coefficients);
            const long double physicalLength = static_cast<long double>(length) * lengthScale_;
            Detail::require(std::isfinite(physicalLength) && physicalLength > 0.0L,
                            "unrepresentable physical piece length");
            const long double halfLength = 0.5L * physicalLength;
            Detail::require(halfLength > 0.0L, "physical Gauss weight underflow");
            if (axisymmetric_)
            {
                const double radialMinus = radialWeightAt(minus), radialPlus = radialWeightAt(plus);
                const long double weightMinus = halfLength * radialMinus,
                                  weightPlus = halfLength * radialPlus;
                Detail::require(std::isfinite(weightMinus) && std::isfinite(weightPlus) &&
                                    (radialMinus == 0.0 || weightMinus > 0.0L) &&
                                    (radialPlus == 0.0 || weightPlus > 0.0L),
                                "RZ Gauss weight overflow or underflow");
                integral += weightMinus * a + weightPlus * b;
            }
            else
                integral += halfLength * (static_cast<long double>(a) + b);
        },
        "[DielectricCubicTraceCandidate] ", geometry);
    const double ownerValue = particular(potentialBasis(link.owner, ownerPhase), ownerPhase, coefficients);
    const double neighborValue =
        particular(potentialBasis(link.neighbor, neighborPhase), neighborPhase, coefficients);
    long double conductance = static_cast<long double>(twoPointPermittivity) * faceLength / linkDistance;
    double faceMeasure = faceLength;
    if (axisymmetric_)
    {
        faceMeasure = axisymmetricFaceMeasure(faceLength, pieces);
        conductance = static_cast<long double>(twoPointPermittivity) * faceMeasure / linkDistance;
    }
    Detail::require(std::isfinite(conductance) &&
                        (axisymmetric_ && faceMeasure == 0.0 ? conductance == 0.0L : conductance > 0.0L),
                    "unrepresentable two-point conductance");
    long double constant = integral - conductance * (static_cast<long double>(ownerValue) - neighborValue);
    const auto sampled = sampleParticular(coefficients);
    for (std::size_t i = 1; i < sampled.size(); ++i)
        constant -= result.weights[i] * (static_cast<long double>(sampled[i]) - sampled[0]);
    result.constant = Detail::checked(constant, "non-finite affine integrated cubic flux correction");
    return result;
}
} // namespace DielectricCubicTrace
} // namespace
} // namespace Quantum::Electricity::Detail
