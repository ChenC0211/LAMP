#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "Modules/Electricity/Detail/DielectricInterfaceTrace.h"

// Shared implementation-only correction to an orthogonal two-point flux,
// not a global electric solver. Planar affine exactness does not establish
// curved-interface accuracy, positive definiteness, monotonicity or global
// stencil stability.
// Implementation-only: the anonymous ancestor gives every definition and
// type TU-local linkage, so different floating-point modes cannot share these
// engines through COMDAT. All includes stay outside these namespaces; no
// Detail type is a cross-TU public API.
namespace Quantum::Electricity::Detail
{
namespace
{
namespace DielectricInterfaceFlux
{
struct Geometry
{
    double kappa = 0.0;
    double tangentMultiplier = 0.0;
};

struct Result
{
    double kappa = 0.0;
    double tangentMultiplier = 0.0;
    // Same order as the supplied trace samples/cell indexes. These are only
    // the additive correction, NOT the full face-flux stencil.
    std::vector<double> correctionWeights;
};

namespace Detail
{
inline void require(bool condition, const char* reason)
{
    if (!condition) throw std::runtime_error(std::string("[DielectricInterfaceFluxCandidate] ") + reason);
}

inline double checked(long double value, const char* reason)
{
    const double result = static_cast<double>(value);
    require(std::isfinite(value) && std::isfinite(result), reason);
    return result;
}

inline void checkNormal(const Vector2D& normal, const char* reason)
{
    const double tolerance = 128.0 * std::numeric_limits<double>::epsilon();
    const double length = std::hypot(normal.x, normal.y);
    require(std::isfinite(normal.x) && std::isfinite(normal.y)
        && std::isfinite(length) && std::abs(length - 1.0) <= tolerance, reason);
}

inline void checkGeometry(const Geometry& geometry)
{
    require(std::isfinite(geometry.kappa) && geometry.kappa > 0.0,
        "harmonic permittivity must be finite and positive");
    require(std::isfinite(geometry.tangentMultiplier), "tangential multiplier must be finite");
}
}

// e is the owner-outward face normal AND the center-to-center link direction;
// d is that link's length and L is the physical face length (2D per unit depth).
// n points negative -> positive; t=(-ny,nx). Fractions come from the caller's
// actual geometry, not from the trace's local tangent-plane phase signs.
//
// For ONE planar interface with constant Et,Dn and phase permittivities:
//   Delta V = -Et*(t.e)*d - Dn*(n.e)*d*(f_d/eps- + (1-f_d)/eps+),
//   q_TP    = -(kappa*L/d)*Delta V,
//   q_exact = L*(eps_a*Et*(t.e) + Dn*(n.e)),
// so q_exact-q_TP = L*(eps_a-kappa)*(t.e)*Et.
// Link and face fractions need not agree. A same-phase pair of centers does
// not imply an uncut physical face or a vanishing correction. On curved
// geometry these are local planar approximations, not exact identities.
// The caller must use the same n/eps-/eps+ when constructing the trace.
inline Geometry evaluateGeometry(const Vector2D& ownerFaceNormal,
    double faceLength, double linkDistance, double negativeLinkFraction,
    double negativeFaceFraction, double epsilonNegative, double epsilonPositive,
    const Vector2D& interfaceNormal)
{
    using Detail::require;
    Detail::checkNormal(ownerFaceNormal, "owner face normal/link direction must have finite unit length");
    Detail::checkNormal(interfaceNormal, "interface normal must have finite unit length");
    require(std::isfinite(faceLength) && faceLength > 0.0
        && std::isfinite(linkDistance) && linkDistance > 0.0,
        "face length and link distance must be finite and positive");
    require(std::isfinite(negativeLinkFraction) && negativeLinkFraction >= 0.0 && negativeLinkFraction <= 1.0
        && std::isfinite(negativeFaceFraction) && negativeFaceFraction >= 0.0 && negativeFaceFraction <= 1.0,
        "negative phase fractions must be finite and in [0,1]");
    require(std::isfinite(epsilonNegative) && epsilonNegative > 0.0
        && std::isfinite(epsilonPositive) && epsilonPositive > 0.0,
        "both permittivities must be finite and positive");

    Geometry result;
    // Preserve the exact homogeneous-material limit, without constructing a
    // trace or introducing a roundoff-sized correction at arbitrary fractions.
    if (epsilonNegative == epsilonPositive)
    {
        result.kappa = epsilonNegative;
        return result;
    }
    if (negativeLinkFraction == 0.0) result.kappa = epsilonPositive;
    else if (negativeLinkFraction == 1.0) result.kappa = epsilonNegative;
    else
    {
        // Normalize by the smaller epsilon to avoid reciprocal overflow.
        // A numerically unrepresentable active ratio is rejected, not silently
        // interpreted as a conductor or silently replaced by a fallback flux.
        const double scale = std::min(epsilonNegative, epsilonPositive);
        const double negativeRatio = scale / epsilonNegative;
        const double positiveRatio = scale / epsilonPositive;
        require(negativeRatio > 0.0 && positiveRatio > 0.0,
            "permittivity contrast is not representable in the harmonic link");
        const long double denominator = static_cast<long double>(negativeLinkFraction) * negativeRatio
            + (1.0L - negativeLinkFraction) * positiveRatio;
        require(std::isfinite(denominator) && denominator > 0.0L,
            "harmonic link denominator must be finite and positive");
        result.kappa = Detail::checked(static_cast<long double>(scale) / denominator,
            "unrepresentable harmonic link permittivity");
    }
    const double arithmetic = negativeFaceFraction == 0.0 ? epsilonPositive
        : negativeFaceFraction == 1.0 ? epsilonNegative
        : Detail::checked(static_cast<long double>(negativeFaceFraction) * epsilonNegative
            + (1.0L - negativeFaceFraction) * epsilonPositive,
            "unrepresentable phase-integrated face permittivity");
    require(arithmetic > 0.0, "phase-integrated face permittivity must be positive");
    const double tangentDotFace = Detail::checked(-static_cast<long double>(interfaceNormal.y) * ownerFaceNormal.x
        + static_cast<long double>(interfaceNormal.x) * ownerFaceNormal.y,
        "non-finite tangent/face-normal dot product");
    const double contrast = arithmetic - result.kappa;
    if (contrast != 0.0 && tangentDotFace != 0.0)
    {
        result.tangentMultiplier = Detail::checked(static_cast<long double>(faceLength) * tangentDotFace * contrast,
            "unrepresentable tangential flux multiplier");
        require(result.tangentMultiplier != 0.0, "tangential flux multiplier underflow");
    }
    Detail::checkGeometry(result);
    return result;
}

// Strictly zero geometry needs no two-phase sample neighborhood or trace QR.
// The caller still assembles q_TP=-(kappa*L/d)*(V_neighbor-V_owner).
inline Result reconstruct(const Geometry& geometry)
{
    Detail::checkGeometry(geometry);
    Detail::require(geometry.tangentMultiplier == 0.0,
        "nonzero tangential correction requires a trace reconstruction");
    Result result;
    result.kappa = geometry.kappa;
    result.tangentMultiplier = geometry.tangentMultiplier;
    return result;
}

inline Result reconstruct(const Geometry& geometry,
    const DielectricInterfaceTrace::Reconstruction& trace,
    const std::vector<int>& sampleCellIndexes)
{
    using Detail::require;
    Detail::checkGeometry(geometry);
    const auto& tangentialWeights = trace.tangentialElectricFieldWeights();
    require(sampleCellIndexes.size() == trace.sampleCount()
        && tangentialWeights.size() == sampleCellIndexes.size() && sampleCellIndexes.size() >= 3,
        "sample cell index count must match the trace samples");
    std::vector<int> sortedIndexes = sampleCellIndexes;
    for (int index : sortedIndexes) require(index >= 0, "sample cell indexes must be nonnegative");
    std::sort(sortedIndexes.begin(), sortedIndexes.end());
    require(std::adjacent_find(sortedIndexes.begin(), sortedIndexes.end()) == sortedIndexes.end(),
        "sample cell indexes must be unique");
    for (double weight : tangentialWeights) require(std::isfinite(weight), "trace tangential weight must be finite");

    Result result;
    result.kappa = geometry.kappa;
    result.tangentMultiplier = geometry.tangentMultiplier;
    result.correctionWeights.assign(tangentialWeights.size(), 0.0);
    if (geometry.tangentMultiplier == 0.0) return result;

    long double tail = 0.0L;
    bool nonzeroTail = false;
    for (std::size_t i = 1; i < tangentialWeights.size(); ++i)
    {
        result.correctionWeights[i] = Detail::checked(static_cast<long double>(geometry.tangentMultiplier)
            * tangentialWeights[i], "unrepresentable tangential flux correction weight");
        nonzeroTail = nonzeroTail || result.correctionWeights[i] != 0.0;
        tail += result.correctionWeights[i];
    }
    require(nonzeroTail, "tangential flux correction weights underflow");
    result.correctionWeights[0] = Detail::checked(-tail, "unrepresentable anchored flux correction weight");
    // Apply as sum_{j>0} w_j*(V_j-V_0), or assemble those exact same differences.
    // This annihilates a common potential exactly when evaluated in anchored
    // form. A naive raw dot product can still suffer common-offset roundoff.
    // Build each full physical-face stencil once, then use its exact opposite
    // for the neighbor. Independent owner/neighbor fits do not imply global
    // conservation even though the planar flux formula itself is consistent.
    return result;
}
}
} // TU-local engine
} // namespace Quantum::Electricity::Detail
