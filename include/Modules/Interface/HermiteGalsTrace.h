#pragma once

#include "Modules/Interface/HermiteGalsState.h"
#include "Core/Geometry/LevelSetBoundaryTrace.h"
#include <optional>

namespace LAMP::Interface
{
class HermiteGalsVelocityView;
class HermiteGalsTraceResult;

// Fixed explicit midpoint RK2, with tangent differentiation of the very same
// stages/substeps. No finite-difference Jacobian or alternate-integrator fallback.
// The two velocity samples bracket [previous.time(), arrival.time()]; for a
// Midpoint layer, arrival must contain velocity AT that physical midpoint.
HermiteGalsTraceResult traceHermiteGalsWithJacobian(
    const LevelSetBoundaryTrace::Domain& domain, Vector2D arrivalPoint,
    const HermiteGalsVelocityView& previous, const HermiteGalsVelocityView& arrival);

// Immutable endpoint velocity contract. Value and Du come from ONE evaluation
// of the same spatial interpolant/branch. Only first derivatives are required.
// Owns the callable, not borrowed capture objects. Providers must keep captures
// alive, immutable and repeatable for every trace (including parallel callers).
class HermiteGalsVelocityView
{
public:
    using Query = LevelSetBoundaryTrace::DifferentiableVelocity;
    HermiteGalsVelocityView(HermiteGalsIdentity geometry, HermiteGalsTimeLayer time, Query query);
    HermiteGalsIdentity geometryIdentity() const noexcept { return geometry_; }
    HermiteGalsTimeLayer timeLayer() const noexcept { return time_; }
private:
    HermiteGalsIdentity geometry_;
    HermiteGalsTimeLayer time_;
    Query query_;
    friend HermiteGalsTraceResult traceHermiteGalsWithJacobian(
        const LevelSetBoundaryTrace::Domain&, Vector2D,
        const HermiteGalsVelocityView&, const HermiteGalsVelocityView&);
};

// Read-only result with both endpoint identities/time layers. The old and new
// velocity may legitimately belong to DIFFERENT geometry revisions. The caller
// binds these views to its coupling candidate; this operator never changes it.
class HermiteGalsTraceResult
{
public:
    Vector2D arrivalPoint() const { return arrivalPoint_; }
    Vector2D departurePoint() const { return trace_.trajectory.point; }
    double elapsed() const noexcept { return trace_.trajectory.elapsed; }
    double departureTime() const;
    bool boundaryHit() const noexcept { return trace_.trajectory.boundaryHit; }
    std::optional<Boundary> boundarySide() const { return trace_.trajectory.side; }
    int acceptedSubsteps() const noexcept { return trace_.acceptedSubsteps; }
    const std::optional<LevelSetBoundaryTrace::Jacobian>& jacobian() const noexcept
    { return trace_.derivative; }
    // Boundary-event time/history has NOT been differentiated. This accessor
    // rejects a hit instead of supplying an invalid interior pullback matrix.
    LevelSetBoundaryTrace::Jacobian interiorJacobian() const;
    HermiteGalsIdentity previousGeometryIdentity() const noexcept { return previousGeometry_; }
    HermiteGalsIdentity arrivalGeometryIdentity() const noexcept { return arrivalGeometry_; }
    HermiteGalsTimeLayer previousTimeLayer() const noexcept { return previousTime_; }
    HermiteGalsTimeLayer arrivalTimeLayer() const noexcept { return arrivalTime_; }
private:
    HermiteGalsTraceResult(Vector2D point, const HermiteGalsVelocityView& previous,
        const HermiteGalsVelocityView& arrival, LevelSetBoundaryTrace::JetResult trace);
    Vector2D arrivalPoint_;
    HermiteGalsIdentity previousGeometry_, arrivalGeometry_;
    HermiteGalsTimeLayer previousTime_, arrivalTime_;
    LevelSetBoundaryTrace::JetResult trace_;
    friend HermiteGalsTraceResult traceHermiteGalsWithJacobian(
        const LevelSetBoundaryTrace::Domain&, Vector2D,
        const HermiteGalsVelocityView&, const HermiteGalsVelocityView&);
};
} // namespace LAMP::Interface
