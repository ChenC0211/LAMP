#pragma once

#include "Modules/Hydrodynamics/CchMomentumGeometry.h"
#include "Core/Geometry/LevelSetBoundaryTrace.h"
#include "Core/Geometry/LevelSetMappedSegment.h"
#include <cstdint>

namespace Quantum::Hydrodynamics
{
namespace Detail { class CchRegularAdvectorGeometry; }
class CchSpatialAdvector;
class CchAdvectorIdentity
{
public:
    CchAdvectorIdentity() = default;
    std::uint64_t value() const noexcept { return value_; }
    friend bool operator==(CchAdvectorIdentity a, CchAdvectorIdentity b) noexcept
    { return a.value_ == b.value_; }
    friend bool operator!=(CchAdvectorIdentity a, CchAdvectorIdentity b) noexcept { return !(a==b); }
private:
    explicit CchAdvectorIdentity(std::uint64_t v) : value_(v) {}
    std::uint64_t value_ = 0;
    friend class CchAdvectorReconstruction;
};
struct CchAdvectorAtom
{
    int normalAxis = 0, primalVelocity = -1, ownerPressure = -1;
    Vector2D begin{0,0}, end{0,0}; // canonical positive tangent orientation
    long double area = 0;
};
struct CchAdvectorPreparationReport
{
    // Legacy tangentUnidentified counts unidentifiable priority-space rows.
    std::size_t svdCount = 0, rankDeficient = 0, tangentUnidentified = 0;
    std::size_t weights = 0, cellPieces = 0, atoms = 0, viewEntries = 0;
    double maxTangentDeficit = 0;
};
// One polynomial support piece of the same C1 physical velocity field.
// Its closed bounds are metadata; values and first derivatives agree at knots.
struct CchAdvectorCellJet
{
    LevelSetBoundaryTrace::VelocityJet jet;
    Vector2D lower{0,0}, upper{0,0};
};
// Actual cell/knot rectangle; its index belongs to this reconstruction only.
struct CchAdvectorPiece
{
    int pressure=-1;
    Vector2D lower{0,0}, upper{0,0};
};
// Independent internal streamfunction mode of this same prepared source.
// Coefficients are velocity state, not fitted from native face means.
struct CchTangentialMode
{
    int pressure=-1,nativeCell=-1;
    Vector2D lower{0,0},upper{0,0};
    double gradientScale=0;
};
// Physical normal velocity in an endpoint-relative distance chart. This is
// metric data, not phase/root authority. Coefficients at both endpoints avoid
// reconstructing an unrepresentable global endpoint for a tiny phase range.
struct CchNormalPolynomial
{
    std::array<long double,5> beginCoefficients{}, endCoefficients{}; // ascending powers through quartic
    long double length = 0, radiusBegin = 0, radiusEnd = 0;
    bool radialMeasure = false, radialTangent = false;
    long double value(bool fromEnd, long double offset) const;
    long double integral(bool fromEnd, long double offset, long double width, int normalSign=1) const;
};
// Explicit homogeneous normal law on caller-selected native boundary faces.
// Static geometry binding only: no velocity, phase, clock or accepted state.
class CchAdvectorBoundary
{
public:
    static CchAdvectorBoundary impermeable(const QuadtreeMesh&, CchMomentumGeometry,
        const std::vector<int>& velocityDofs);
    const CchMomentumGeometry& momentumGeometry() const noexcept;
    const std::vector<int>& velocityDofs() const noexcept;
    bool sameIdentity(const CchAdvectorBoundary&) const noexcept;
    void validateFor(const QuadtreeMesh&, const CchMomentumGeometry&) const;
private:
    struct Data;
    explicit CchAdvectorBoundary(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
};
// Base prepare binds mesh geometry only. withTangentialModes additionally
// retains the authoritative phase publication/clock; it never owns velocity.
// Copies retain the same momentum geometry and immutable source weights.
class CchAdvectorReconstruction
{
public:
    static CchAdvectorReconstruction prepare(const QuadtreeMesh&, CchMomentumGeometry);
    CchAdvectorReconstruction withTangentialModes(const QuadtreeMesh&,CchPhaseGeometry) const;
    const CchPhaseGeometry* phaseGeometry() const noexcept;
    const std::vector<CchTangentialMode>& tangentialModes() const noexcept;
    int tangentialModeForPressure(int pressure) const;
    // Physical one-sided phase branch of the same retained Q4 curl.
    LevelSetBoundaryTrace::VelocityJet tangentialModeJet(int mode,Vector2D point,bool liquid) const;
    long double tangentialModeStreamfunction(int mode,Vector2D point,bool liquid) const;
    const CchMomentumGeometry& momentumGeometry() const noexcept;
    const std::vector<CchAdvectorAtom>& atoms() const noexcept;
    const std::vector<CchAdvectorPiece>& pieces() const noexcept;
    const std::vector<std::vector<int>>& primalAtoms() const noexcept;
    const std::vector<std::vector<int>>& dualAtoms() const noexcept;
    const CchAdvectorPreparationReport& report() const noexcept;
    // Same immutable preparation used by every actual value/J/atom query.
    const Detail::CchRegularAdvectorGeometry& regularGeometry() const noexcept;
    bool sameIdentity(const CchAdvectorReconstruction&) const noexcept;
    CchSpatialAdvector evaluate(const QuadtreeMesh&, const CchNativeValues&,
        Interface::HermiteGalsIdentity velocityGeometry, Interface::HermiteGalsTimeLayer fullClock) const;
    CchSpatialAdvector evaluate(const QuadtreeMesh&, const CchNativeValues&,
        Interface::HermiteGalsIdentity velocityGeometry, Interface::HermiteGalsTimeLayer fullClock,
        const CchAdvectorBoundary&) const;
    // Explicit modal initial/current state; omitting q on a modal source rejects.
    CchSpatialAdvector evaluate(const QuadtreeMesh&,const CchNativeValues&,
        const std::vector<double>& tangentialVelocity,Interface::HermiteGalsIdentity,
        Interface::HermiteGalsTimeLayer,const CchAdvectorBoundary* =nullptr) const;
    // Dense sample of the SAME endpoint coefficients, bound to a derived
    // geometry at clock.time(). No new fit or native-to-field application.
    CchSpatialAdvector interpolate(const QuadtreeMesh&,const CchSpatialAdvector& previous,
        const CchSpatialAdvector& arrival,Interface::HermiteGalsIdentity geometry,
        Interface::HermiteGalsTimeLayer clock) const;
private:
    CchSpatialAdvector evaluateBound(const QuadtreeMesh&, const CchNativeValues&,
        Interface::HermiteGalsIdentity, Interface::HermiteGalsTimeLayer, const CchAdvectorBoundary*,
        const std::vector<double>*) const;
    struct Data;
    explicit CchAdvectorReconstruction(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
    friend class CchSpatialAdvector;
};
// Caller-declared provenance is checked exactly; native values themselves
// have no physical clock. This owner does not certify an H14 coupled commit.
class CchSpatialAdvector
{
public:
    const CchAdvectorReconstruction& reconstruction() const noexcept;
    const CchMomentumGeometry& momentumGeometry() const noexcept;
    const CchNativeValues& values() const noexcept;
    const std::vector<double>& tangentialVelocity() const noexcept;
    Interface::HermiteGalsIdentity geometryIdentity() const noexcept;
    Interface::HermiteGalsTimeLayer timeLayer() const noexcept;
    CchAdvectorIdentity identity() const noexcept;
    const CchAdvectorBoundary* boundary() const noexcept; // null means no prescribed normal law
    // Nonzero q rejects this legacy single-polynomial port. Combined phase
    // flux uses regularAtomPolynomials + tangentialAtomRate.
    const std::vector<CchNormalPolynomial>& atomPolynomials() const;
    // Explicit base contribution for the combined physical phase-rate owner.
    const std::vector<CchNormalPolynomial>& regularAtomPolynomials() const noexcept;
    // Exact curl endpoint formula on the retained Q4 trace. Gas contribution;
    // liquid is zero. Native boundary face atoms are identically zero.
    long double tangentialAtomRate(int atom,const LevelSetMappedSegment::NativeQuarticTrace&) const;
    // Original D_I cell average; not a pointwise divergence constraint.
    const std::vector<long double>& cellDivergence() const noexcept;
    Vector2D cellVelocity(int pressure, Vector2D point) const;
    CchAdvectorCellJet cellVelocityJet(int pressure, Vector2D point) const;
    // Physical phase view. cellVelocityJet remains the C1 geometry extension;
    // modal interface normal and native face means are exactly zero.
    CchAdvectorCellJet cellPhysicalVelocityJet(int pressure,Vector2D point,bool liquid) const;
    // Explicit one-sided branch, including BOTH endpoints of its closure.
    CchAdvectorCellJet pieceVelocityJet(int piece, Vector2D point) const;
    // Same physical field at the actual face center, distinct from its mean.
    Vector2D nativeVelocity(int velocity) const;
    void validateFor(const QuadtreeMesh&, const CchMomentumGeometry& expected,
        Interface::HermiteGalsIdentity, Interface::HermiteGalsTimeLayer) const;
private:
    struct Data;
    explicit CchSpatialAdvector(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
    friend class CchAdvectorReconstruction;
};
}
