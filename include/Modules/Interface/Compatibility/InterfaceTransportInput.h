#pragma once

#include <functional>
#include <memory>
#include "Core/Base/Vector2D.h"
#include "Core/Geometry/UnitPolynomialRoots.h"
#include <stdexcept>

class NativeLeafDomain;
class LevelSetMappedSegment;

// Borrowed selected domain and scalar boundary constraints. Pressure, density,
// viscosity and momentum unknowns are deliberately absent. The owner of the
// domain and boundary callback must outlive all synchronous consumers and keep
// them immutable. This view does not create another geometry or face array.
class InterfaceTransportDomain
{
public:
    enum class Boundary { Interior, StationaryWall, Symmetry, SinglePhasePort, PrescribedScalar };
    using BoundaryQuery = std::function<Boundary(int nativeFace)>;
    InterfaceTransportDomain(const NativeLeafDomain& geometry, BoundaryQuery boundary);
    const NativeLeafDomain& geometry() const noexcept { return *geometry_; }
    Boundary boundary(int nativeFace) const;
private:
    const NativeLeafDomain* geometry_;
    BoundaryQuery boundary_;
};

// Frozen normal-flux input, not a velocity solver. Callbacks may borrow an
// existing immutable reconstruction; they are not retained by phase-volume
// correction or its result. The provider must keep its data alive/unchanged for
// the call. No callback here advances time, refits a field or accepts a step.
//
// All required queries and the optional root query describe ONE native trace: +r/+z orientation, actual two
// half-face charts and the domain's XY or full RZ measure. Integrals and density must come from
// that same primitive. Phase cutting checks exact chart/source correspondence;
// this contract does not invent a time/version certificate for an arbitrary
// user callback. The driving step supplies the correct frozen velocity layer.
class InterfaceNormalFluxView
{
public:
    struct FaceChart { Vector2D begin, midpoint, end; };
    struct Queries
    {
        std::function<FaceChart(int)> chart;
        // t=anchor0+(anchor1-anchor0)*s. Preserve factored small intervals;
        // neither compose/round their endpoints nor multiply a second metric.
        std::function<long double(int,int,double,double,double,double)> integral;
        std::function<long double(int,int,const LevelSetMappedSegment&)> phaseIntegral;
        // Derivative of the SAME integral w.r.t. the half-face parameter.
        std::function<long double(int,int,bool,long double)> density;
        // Optional exact-representation capability, not reconstructed from
        // samples of an arbitrary callback. Bounds describe that SAME density.
        std::function<Quantum::Geometry::UnitPolynomialRoots(int,int)> roots;
        // Same retained basis, for interpolation at actual quadrature times.
        std::function<std::vector<Quantum::Geometry::UnitPolynomialTerm>(int,int)> polynomial;
        // Optional factored range primitive for non-polynomial traces. Native
        // half coordinate t=a0+(a1-a0)*s, s=offset+q or 1-offset-q.
        std::function<long double(int,int,double,double,bool,long double,long double)> rangeIntegral;
        // Optional analytic bound of |density| on the entire native half.
        // Polynomial views can derive it from their retained coefficients.
        std::function<long double(int,int)> densityMagnitudeBound;
        // Required signed face-moment degree for the consuming momentum
        // space (0 or 2), in the whole-face unit tangential coordinate.
        int retainedMomentDegree=0;
        // Complete Q8 velocity products use a bounded Bernstein basis.
        // This only requests extra integrals of this SAME phase flux.
        bool retainedVelocityProducts=false;
    };
    explicit InterfaceNormalFluxView(Queries queries);
    // Exact separable capability q(face,s,t)=a(t)*q0(face,s). The factory
    // scales all primitives/densities from ONE retained spatial view. Merely
    // similar callbacks do not establish identity between time layers.
    static InterfaceNormalFluxView scaled(std::shared_ptr<const InterfaceNormalFluxView> spatial,double amplitude);
    bool canInterpolateRoots(const InterfaceNormalFluxView&) const noexcept;
    Quantum::Geometry::UnitPolynomialRoots interpolatedRoots(const InterfaceNormalFluxView& next,
        int row,int half,double fraction) const;
    long double interpolatedIntegral(const InterfaceNormalFluxView& next,int row,int half,double fraction,
        double lo,double hi,double partLo=0,double partHi=1) const;
    bool hasRangeIntegral() const noexcept {return bool(queries_.rangeIntegral) || hasNormalFluxPolynomial();}
    long double integratedRange(int row,int half,double a0,double a1,bool fromEnd,long double offset,long double length) const;
    long double interpolatedRange(const InterfaceNormalFluxView& next,int row,int half,double fraction,
        double a0,double a1,bool fromEnd,long double offset,long double length) const;
    // Shared scalar amplitude has at most one temporal zero. Polynomial
    // isolation retains its enclosure; general polynomial traces use spatial
    // isolation at each quadrature time and supply no global temporal cut.
    Quantum::Geometry::UnitPolynomialRoots interpolationTimeRoots(const InterfaceNormalFluxView& next) const;
    FaceChart faceChart(int row) const { return queries_.chart(row); }
    long double integratedNormalFlux(int row, int half, double begin, double end) const
    { return integratedNormalFlux(row, half, 0.0, 1.0, begin, end); }
    long double integratedNormalFlux(int row, int half, double anchor0, double anchor1,
        double begin, double end) const { return queries_.integral(row,half,anchor0,anchor1,begin,end); }
    long double integratedNormalFlux(int row, int half, const LevelSetMappedSegment& segment) const
    { return queries_.phaseIntegral(row,half,segment); }
    long double normalFluxDensity(int row, int half, bool fromEnd, long double offset) const
    { return queries_.density(row,half,fromEnd,offset); }
    bool hasNormalFluxRoots() const noexcept { return static_cast<bool>(queries_.roots); }
    bool hasNormalFluxPolynomial() const noexcept { return static_cast<bool>(queries_.polynomial); }
    int retainedMomentDegree() const noexcept {return queries_.retainedMomentDegree;}
    bool retainedVelocityProducts() const noexcept {return queries_.retainedVelocityProducts;}
    long double densityMagnitudeBound(int row,int half) const;
    std::vector<Quantum::Geometry::UnitPolynomialTerm> normalFluxPolynomial(int row,int half) const
    {
        if (!queries_.polynomial) throw std::logic_error("normal flux view has no retained polynomial");
        return queries_.polynomial(row,half);
    }
    Quantum::Geometry::UnitPolynomialRoots normalFluxRoots(int row,int half) const
    {
        if (!queries_.roots) throw std::logic_error("normal flux view has no retained root capability");
        return queries_.roots(row,half);
    }
private:
    Queries queries_;
    std::shared_ptr<const InterfaceNormalFluxView> spatial_;
    double amplitude_=1;
};
