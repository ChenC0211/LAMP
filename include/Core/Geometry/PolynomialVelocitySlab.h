#pragma once

#include "Core/Geometry/LevelSetBoundaryTrace.h"
#include <array>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace LAMP::Geometry
{
// Immutable reconstruction input after copying into a slab, not global DOFs.
// coefficient[4*i+j] multiplies xi^i eta^j; each entry contains (u_x,u_y).
struct PolynomialVelocityData
{
    Vector2D origin{0,0},scale{1,1};
    std::array<Vector2D,16> previous{},arrival{};
};
struct PolynomialVelocityRange {double lower=0,upper=0;};
struct PolynomialVelocityBounds
{
    std::array<PolynomialVelocityRange,2> value;
    std::array<PolynomialVelocityRange,4> derivative; // xx,xy,yx,yy
};
struct PolynomialRk2Bounds
{
    Vector2D arrivalLower,arrivalUpper,midpointLower,midpointUpper,departureLower,departureUpper;
    std::array<PolynomialVelocityRange,4> jacobian; // xx,xy,yx,yy
    double maximumL1Travel=0; // conservative sufficient CFL bound, not an ODE error estimate
};
struct PolynomialReachableSourceBox
{
    Vector2D lower,upper;
    std::array<bool,4> boundary; // west,south,east,north; conservative possible sources
};
// Finite-speed/rounding source support ONLY. Does not certify completion,
// uniqueness of a boundary event, B or J. Used for constant-content algebra.
class PolynomialRk2Reachability
{
public:
    PolynomialReachableSourceBox cover(Vector2D lower,Vector2D upper) const;
private:
    friend class PolynomialVelocitySlab;
    PolynomialRk2Reachability(Vector2D lo,Vector2D hi,Vector2D radius)
        :lower_(lo),upper_(hi),radius_(radius){}
    Vector2D lower_,upper_,radius_;
};
struct PolynomialAffineCoordinate
{
    double center=0,tangent=0,time=0;
    PolynomialVelocityRange remainder;
    PolynomialVelocityRange enclose(PolynomialVelocityRange tangentParameter,
        PolynomialVelocityRange timeParameter) const;
};
// Strict composition of ONE tensor Q4 polynomial with two shared-parameter
// physical forms. Coefficients [5*y+x] enclose the actual chart powers, not
// fitted data. Evaluates the polynomial extension; owns no donor/phase/mesh.
PolynomialAffineCoordinate encloseTensorQuarticAffine(
    const std::array<std::array<double,2>,25>& chartPowerCoefficients,Vector2D cellLower,Vector2D cellUpper,
    const std::array<PolynomialAffineCoordinate,2>& physicalCoordinates);
// q_s/q_t are the exact real affine coordinates of the retained double face
// and time endpoints, in [-1,1]. Both components share the SAME parameters.
struct PolynomialAffineRk2FaceBounds
{
    Vector2D faceBegin,faceEnd;
    double timeBegin=0,timeEnd=0;
    PolynomialRk2Bounds branch;
    std::array<PolynomialAffineCoordinate,2> departure;
    std::array<PolynomialVelocityRange,2> encloseDeparture(Vector2D point,double time) const;
};
// A valid box whose complete original RK2 branch has not been proved. A
// caller may subdivide, but may not discard this condition on an error budget.
class PolynomialRk2EnclosureUnresolved : public std::invalid_argument
{
public:
    explicit PolynomialRk2EnclosureUnresolved(const std::string& reason);
};
struct PolynomialFaceIntegralBounds
{
    // Integrated total+, total-, weighted-liquid+, weighted-liquid-.
    std::array<PolynomialVelocityRange,4> volumes{};
    PolynomialFaceIntegralBounds added(const PolynomialFaceIntegralBounds&) const;
    PolynomialFaceIntegralBounds intersected(const PolynomialFaceIntegralBounds&) const;
    PolynomialFaceIntegralBounds withContentRange(long double lower,long double upper) const;
    double errorBound() const; // outward sum of all four widths
    std::array<double,4> representative() const; // inside bounds, liquid <= total
};
// Strict subdivision helper; nullopt means no represented interior midpoint.
std::optional<double> polynomialFaceMidpoint(double lower,double upper);

// Tensor-cubic space, linear endpoint time reconstruction. Strict out-of-line
// evaluation and outward bounds of the SAME data/first derivatives. Bounds
// cover closed boxes, including degenerate faces/points. They are NOT a
// certificate of a characteristic, physical boundary history or phase flux.
class PolynomialVelocitySlab
{
public:
    PolynomialVelocitySlab(double begin,double end,const PolynomialVelocityData&);
    double timeBegin() const noexcept {return begin_;}
    double timeEnd() const noexcept {return end_;}
    const PolynomialVelocityData& data() const noexcept {return data_;}
    LevelSetBoundaryTrace::VelocityJet evaluate(Vector2D point,double time) const;
    PolynomialVelocityBounds enclose(Vector2D lower,Vector2D upper,double begin,double end) const;
    // Same Horner/normalization/mixing expressions, additionally intersected
    // with exact sign facts preserved by their rounded operations. This is
    // not a small-value tolerance or a characteristic/phase certificate.
    PolynomialVelocityBounds encloseSignPreserving(Vector2D lower,Vector2D upper,double begin,double end) const;
    // Source support over all old->arrival durations. Sufficient global CFL;
    // no opaque coverage; RZ restricted to identically zero radial velocity.
    // Does not run or replace the original characteristic/Jacobian.
    PolynomialRk2Reachability encloseRk2Reachability(const LevelSetBoundaryTrace::Domain&,
        double arrivalEnd) const;
    // Succeeds only when the existing tracer's complete, unshortened no-hit
    // branch is covered for EVERY input point and arrival time. Unresolved
    // domain/axis/CFL/opaque-coverage cases reject; no clipping or fallback.
    PolynomialRk2Bounds encloseInteriorRk2(const LevelSetBoundaryTrace::Domain&,
        Vector2D lower,Vector2D upper,double arrivalBegin,double arrivalEnd) const;
    // Same branch/expression, additionally intersects time mixing with an
    // endpoint hull plus an outward rounding envelope. Original query stays
    // numerically unchanged for its existing consumers.
    PolynomialRk2Bounds encloseInteriorRk2TightTime(const LevelSetBoundaryTrace::Domain&,
        Vector2D lower,Vector2D upper,double arrivalBegin,double arrivalEnd) const;
    // Every point in this closed box, at every POSITIVE duration in the
    // interval, takes the original tracer's unique immediate boundary event.
    // If arrivalBegin==timeBegin(), the zero-duration identity time slice is
    // excluded from this certificate (zero space-time measure). No later
    // event, exterior velocity, scalar or Jacobian is certified here.
    Boundary encloseInitialBoundaryHit(const LevelSetBoundaryTrace::Domain&,
        Vector2D lower,Vector2D upper,double arrivalBegin,double arrivalEnd) const;
    // Correlated pointwise affine form, with a rigorously bounded nonlinear
    // and original floating-operation remainder. Same complete branch guard.
    PolynomialAffineRk2FaceBounds encloseAffineFaceRk2(const LevelSetBoundaryTrace::Domain&,
        Vector2D faceBegin,Vector2D faceEnd,double arrivalBegin,double arrivalEnd) const;
    // Same reconstructed polynomial, integrated over an ascending native
    // axis-aligned face and time interval. Content bounds must be supplied by
    // an independent phase/transport proof; this Core method owns no phase.
    PolynomialFaceIntegralBounds encloseNormalFaceIntegral(Vector2D beginPoint,Vector2D endPoint,
        double timeBegin,double timeEnd,bool radial,long double contentLower=0,long double contentUpper=1) const;
    // All forms use this face/time parameterization and together enclose
    // every possible donor scalar. The caller proves complete source/branch
    // coverage. Exact half-plane subset measures refine the same gross bounds.
    PolynomialFaceIntegralBounds encloseNormalFaceIntegralWithAffinePhase(Vector2D beginPoint,Vector2D endPoint,
        double timeBegin,double timeEnd,bool radial,const std::vector<PolynomialAffineCoordinate>& phaseForms,
        long double liquidLower=1,long double liquidUpper=1,long double gasLower=0,long double gasUpper=0) const;
private:
    PolynomialVelocityBounds encloseTightTime(Vector2D lower,Vector2D upper,double begin,double end) const;
    PolynomialRk2Bounds encloseInteriorRk2Impl(const LevelSetBoundaryTrace::Domain&,
        Vector2D lower,Vector2D upper,double arrivalBegin,double arrivalEnd,bool tightTime) const;
    double begin_,end_;
    PolynomialVelocityData data_;
};
}
