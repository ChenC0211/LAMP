#pragma once

#include "Modules/Hydrodynamics/HydroMomentumGeometry.h"
#include "Core/Geometry/LevelSetBoundaryTrace.h"
#include "Core/Geometry/LevelSetMappedSegment.h"
#include <cstdint>
#include <memory>
#include <utility>
#include <vector>

namespace LAMP::Hydrodynamics
{
namespace Detail { class HydroRegularAdvectorGeometry; }
class HydroSpatialAdvector;
class HydroAdvectorIdentity
{
public:
    HydroAdvectorIdentity() = default;
    std::uint64_t value() const noexcept { return value_; }
    friend bool operator==(HydroAdvectorIdentity a, HydroAdvectorIdentity b) noexcept
    { return a.value_ == b.value_; }
    friend bool operator!=(HydroAdvectorIdentity a, HydroAdvectorIdentity b) noexcept { return !(a==b); }
private:
    explicit HydroAdvectorIdentity(std::uint64_t v) : value_(v) {}
    std::uint64_t value_ = 0;
    friend class HydroAdvectorReconstruction;
};
struct HydroAdvectorAtom
{
    int normalAxis = 0, primalVelocity = -1, ownerPressure = -1;
    Vector2D begin{0,0}, end{0,0}; // canonical positive tangent orientation
    long double area = 0;
};
struct HydroAdvectorPreparationReport
{
    // Retained preparation structure, not work repeated by a cached factory.
    // Legacy tangentUnidentified counts unidentifiable priority-space rows.
    std::size_t svdCount = 0, rankDeficient = 0, tangentUnidentified = 0;
    std::size_t weights = 0, cellPieces = 0, atoms = 0, viewEntries = 0;
    double maxTangentDeficit = 0;
};
// One polynomial support piece of the same C1 physical velocity field.
// Its closed bounds are metadata; values and first derivatives agree at knots.
struct HydroAdvectorCellJet
{
    LevelSetBoundaryTrace::VelocityJet jet;
    Vector2D lower{0,0}, upper{0,0};
};
// Both original endpoint fields on one actual shared reconstruction piece.
struct HydroAdvectorCellJetPair
{
    HydroAdvectorCellJet previous,arrival;
};
// Actual cell/knot rectangle; its index belongs to this reconstruction only.
struct HydroAdvectorPiece
{
    int pressure=-1;
    Vector2D lower{0,0}, upper{0,0};
};
// Independent internal streamfunction mode of this same prepared source.
// Coefficients are velocity state, not fitted from native face means.
struct HydroTangentialMode
{
    int pressure=-1,nativeCell=-1;
    Vector2D lower{0,0},upper{0,0};
    double gradientScale=0;
};
// Physical normal velocity in an endpoint-relative distance chart. This is
// metric data, not phase/root authority. Coefficients at both endpoints avoid
// reconstructing an unrepresentable global endpoint for a tiny phase range.
struct HydroNormalPolynomial
{
    std::array<long double,5> beginCoefficients{}, endCoefficients{}; // ascending powers through quartic
    long double length = 0, radiusBegin = 0, radiusEnd = 0;
    bool radialMeasure = false, radialTangent = false;
    long double value(bool fromEnd, long double offset) const;
    long double integral(bool fromEnd, long double offset, long double width, int normalSign=1) const;
};
// Integrals in the original curve parameter theta, including span width.
// K rows are velocity components, columns arrival coordinates: xx,xy,yx,yy.
struct HydroAdvectorCurveIntegral
{
    std::array<long double,2> V0{},V1{}; // integral u, theta*u
    std::array<long double,4> K0{},K1{}; // integral Du*P, theta*Du*P
};
// Value-only integrals of two original endpoint fields on one actual curve.
// Neither field, geometry nor clock is owned by this temporary result.
struct HydroAdvectorCurveIntegralPair
{
    HydroAdvectorCurveIntegral previous,arrival;
};
// Explicit homogeneous normal law on caller-selected native boundary faces.
// Static geometry binding only: no velocity, phase, clock or accepted state.
class HydroAdvectorBoundary
{
public:
    static HydroAdvectorBoundary impermeable(const QuadtreeMesh&, HydroMomentumGeometry,
        const std::vector<int>& velocityDofs);
    const HydroMomentumGeometry& momentumGeometry() const noexcept;
    const std::vector<int>& velocityDofs() const noexcept;
    bool completeBoundaryLines() const noexcept; // immutable coplanar native-boundary mask, prepared once
    bool sameIdentity(const HydroAdvectorBoundary&) const noexcept;
    void validateFor(const QuadtreeMesh&, const HydroMomentumGeometry&) const;
private:
    struct Data;
    explicit HydroAdvectorBoundary(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
    friend class HydroAdvectorReconstruction;
};
// Base prepare binds mesh geometry only. withTangentialModes additionally
// retains the authoritative phase publication/clock; it never owns velocity.
// Copies retain the same momentum geometry and immutable source weights.
class HydroAdvectorReconstruction
{
public:
    static HydroAdvectorReconstruction prepare(const QuadtreeMesh&, HydroMomentumGeometry);
    // Select the one prepared representation from immutable physical boundary
    // geometry. Partial boundary lines use the existing constrained P2 source.
    static HydroAdvectorReconstruction prepare(const QuadtreeMesh&, HydroMomentumGeometry,
        const HydroAdvectorBoundary&);
    // Same bound geometry with a boundary-owner retained bare preparation.
    // Mesh/boundary are validated on every call; every wrapper has a new
    // identity. report() describes retained SVD/weights/atoms, not new work.
    static HydroAdvectorReconstruction prepareCached(const QuadtreeMesh&, HydroMomentumGeometry,
        const HydroAdvectorBoundary&);
    // Geometry-only/null-boundary preparation. The boundary handle is only
    // its validated lifetime/cache owner; no wall law enters the source.
    static HydroAdvectorReconstruction prepareUnboundCached(const QuadtreeMesh&, HydroMomentumGeometry,
        const HydroAdvectorBoundary& cacheOwner);
    HydroAdvectorReconstruction withTangentialModes(const QuadtreeMesh&,HydroPhaseGeometry) const;
    const HydroPhaseGeometry* phaseGeometry() const noexcept;
    const std::vector<HydroTangentialMode>& tangentialModes() const noexcept;
    int tangentialModeForPressure(int pressure) const;
    // Physical one-sided phase branch of the same retained Q4 curl.
    LevelSetBoundaryTrace::VelocityJet tangentialModeJet(int mode,Vector2D point,bool liquid) const;
    long double tangentialModeStreamfunction(int mode,Vector2D point,bool liquid) const;
    const HydroMomentumGeometry& momentumGeometry() const noexcept;
    const std::vector<HydroAdvectorAtom>& atoms() const noexcept;
    const std::vector<HydroAdvectorPiece>& pieces() const noexcept;
    const std::vector<std::vector<int>>& primalAtoms() const noexcept;
    const std::vector<std::vector<int>>& dualAtoms() const noexcept;
    const HydroAdvectorPreparationReport& report() const noexcept;
    // Same immutable preparation used by every actual value/J/atom query.
    const Detail::HydroRegularAdvectorGeometry& regularGeometry() const noexcept;
    bool sameIdentity(const HydroAdvectorReconstruction&) const noexcept;
    HydroSpatialAdvector evaluate(const QuadtreeMesh&, const HydroNativeValues&,
        Interface::HermiteGalsIdentity velocityGeometry, Interface::HermiteGalsTimeLayer fullClock) const;
    HydroSpatialAdvector evaluate(const QuadtreeMesh&, const HydroNativeValues&,
        Interface::HermiteGalsIdentity velocityGeometry, Interface::HermiteGalsTimeLayer fullClock,
        const HydroAdvectorBoundary&) const;
    // Explicit modal initial/current state; omitting q on a modal source rejects.
    HydroSpatialAdvector evaluate(const QuadtreeMesh&,const HydroNativeValues&,
        const std::vector<double>& tangentialVelocity,Interface::HermiteGalsIdentity,
        Interface::HermiteGalsTimeLayer,const HydroAdvectorBoundary* =nullptr) const;
    // Dense sample of the SAME endpoint coefficients, bound to a derived
    // geometry at clock.time(). No new fit or native-to-field application.
    HydroSpatialAdvector interpolate(const QuadtreeMesh&,const HydroSpatialAdvector& previous,
        const HydroSpatialAdvector& arrival,Interface::HermiteGalsIdentity geometry,
        Interface::HermiteGalsTimeLayer clock) const;
private:
    static HydroAdvectorReconstruction prepareBound(const QuadtreeMesh&,HydroMomentumGeometry,
        const HydroAdvectorBoundary* boundary,const HydroAdvectorBoundary* cacheOwner);
    HydroSpatialAdvector evaluateBound(const QuadtreeMesh&, const HydroNativeValues&,
        Interface::HermiteGalsIdentity, Interface::HermiteGalsTimeLayer, const HydroAdvectorBoundary*,
        const std::vector<double>*) const;
    struct Preparation;
    struct Data;
    explicit HydroAdvectorReconstruction(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
    friend class HydroSpatialAdvector;
    friend class HydroAdvectorBoundary;
};
// Caller-declared provenance is checked exactly; native values themselves
// have no physical clock. This owner does not certify an H14 coupled commit.
class HydroSpatialAdvector
{
public:
    const HydroAdvectorReconstruction& reconstruction() const noexcept;
    const HydroMomentumGeometry& momentumGeometry() const noexcept;
    const HydroNativeValues& values() const noexcept;
    const std::vector<double>& tangentialVelocity() const noexcept;
    Interface::HermiteGalsIdentity geometryIdentity() const noexcept;
    Interface::HermiteGalsTimeLayer timeLayer() const noexcept;
    HydroAdvectorIdentity identity() const noexcept;
    const HydroAdvectorBoundary* boundary() const noexcept; // null means no prescribed normal law
    // Nonzero q rejects this legacy single-polynomial port. Combined phase
    // flux uses regularAtomPolynomials + tangentialAtomRate.
    const std::vector<HydroNormalPolynomial>& atomPolynomials() const;
    // Explicit base contribution for the combined physical phase-rate owner.
    const std::vector<HydroNormalPolynomial>& regularAtomPolynomials() const noexcept;
    // Exact curl endpoint formula on the retained Q4 trace. Gas contribution;
    // liquid is zero. Native boundary face atoms are identically zero.
    long double tangentialAtomRate(int atom,const LevelSetMappedSegment::NativeQuarticTrace&) const;
    // Original D_I cell average; not a pointwise divergence constraint.
    const std::vector<long double>& cellDivergence() const noexcept;
    Vector2D cellVelocity(int pressure, Vector2D point) const;
    HydroAdvectorCellJet cellVelocityJet(int pressure, Vector2D point) const;
    // Physical phase view. cellVelocityJet remains the C1 geometry extension;
    // modal interface normal and native face means are exactly zero.
    HydroAdvectorCellJet cellPhysicalVelocityJet(int pressure,Vector2D point,bool liquid) const;
    // Explicit one-sided branch, including BOTH endpoints of its closure.
    HydroAdvectorCellJet pieceVelocityJet(int piece, Vector2D point) const;
    // Share geometric factors, retaining each endpoint's original field and
    // the selected piece's closed bounds. Requires the same reconstruction.
    HydroAdvectorCellJetPair pieceVelocityJetPair(const HydroSpatialAdvector& arrival,
        int piece,Vector2D point) const;
    // Same actual piece polynomial along r(theta)=r0+r1*theta+r2*theta^2.
    // P(theta) is the fixed-theta arrival derivative of this same curve.
    // Endpoint powers are retained lazily once by this immutable field owner.
    HydroAdvectorCurveIntegral pieceCurveIntegral(int piece,Vector2D r0,Vector2D r1,Vector2D r2,
        LevelSetBoundaryTrace::Jacobian P0,LevelSetBoundaryTrace::Jacobian P1,
        LevelSetBoundaryTrace::Jacobian P2,double lo,double hi) const;
    // Share only this span's geometric moments; retain both immutable fields.
    HydroAdvectorCurveIntegralPair pieceCurveIntegralPair(const HydroSpatialAdvector& arrival,
        int piece,Vector2D r0,Vector2D r1,Vector2D r2,
        LevelSetBoundaryTrace::Jacobian P0,LevelSetBoundaryTrace::Jacobian P1,
        LevelSetBoundaryTrace::Jacobian P2,double lo,double hi) const;
    // Same physical field at the actual face center, distinct from its mean.
    Vector2D nativeVelocity(int velocity) const;
    void validateFor(const QuadtreeMesh&, const HydroMomentumGeometry& expected,
        Interface::HermiteGalsIdentity, Interface::HermiteGalsTimeLayer) const;
private:
    struct Data;
    explicit HydroSpatialAdvector(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
    friend class HydroAdvectorReconstruction;
};
// Closed spatial support of one ordered atom group, bound to its actual slab.
// Callers cannot fabricate the owner, IDs or cell coverage.
class HydroAdvectorAtomSupport
{
public:
    const std::vector<int>& atomIds() const noexcept { return atomIds_; }
    const std::vector<int>& nativeCells() const noexcept { return nativeCells_; }
private:
    HydroAdvectorAtomSupport(std::shared_ptr<const void> owner,std::vector<int> ids,std::vector<int> cells)
        :owner_(std::move(owner)),atomIds_(std::move(ids)),nativeCells_(std::move(cells)){}
    std::shared_ptr<const void> owner_;
    std::vector<int> atomIds_,nativeCells_;
    friend class HydroAdvectorTimeSlab;
};
// One authenticated immutable endpoint slab for requested local phase atoms.
// No field fit, full dense sample, new velocity ID or mutable mesh is retained.
// Construction binds the actual original publication and complete endpoint
// clocks; query work scales with the requested atoms and their closed support.
class HydroAdvectorTimeSlab
{
public:
    static HydroAdvectorTimeSlab prepare(const QuadtreeMesh&,Interface::HermiteGalsPublication,
        HydroAdvectorReconstruction,const HydroSpatialAdvector& previous,const HydroSpatialAdvector& arrival);
    const HydroAdvectorReconstruction& reconstruction() const noexcept;
    Interface::HermiteGalsIdentity sourceIdentity() const noexcept;
    Interface::HermiteGalsIdentity arrivalGeometryIdentity() const noexcept;
    HydroAdvectorIdentity previousIdentity() const noexcept;
    HydroAdvectorIdentity arrivalIdentity() const noexcept;
    Interface::HermiteGalsTimeLayer endpointTimeLayer() const noexcept;
    void validateSample(Interface::HermiteGalsIdentity geometry,Interface::HermiteGalsIdentity source,
        Interface::HermiteGalsIdentity arrival,Interface::HermiteGalsTimeLayer) const;
    HydroNormalPolynomial atomPolynomial(int atom,Interface::HermiteGalsTimeLayer) const;
    HydroAdvectorAtomSupport prepareAtomSupport(const std::vector<int>& atomIds) const;
    void validateAtomSupport(const HydroAdvectorAtomSupport&) const;
    std::vector<int> nativeCells(const std::vector<int>& atomIds) const;
    std::size_t nativeCellCount() const noexcept;
    std::array<Vector2D,2> nativeCellBounds(int cell) const;
private:
    struct Data;
    explicit HydroAdvectorTimeSlab(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
};
}
