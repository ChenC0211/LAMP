#pragma once

#include "Core/Geometry/ActiveQuarticHermitePatch.h"
#include "Core/Geometry/LevelSetBoundaryTrace.h"
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include <cstddef>
#include <functional>
#include <memory>

class QuadtreeMesh;

namespace ActiveQuarticHermiteGals
{
enum class Integrator { Midpoint, ClosedRk4 };
struct Options
{
    bool radialAxisSymmetry=false;
    // Both callbacks and all captured fields must allow concurrent immutable
    // reads. False retains serial callback order; no parallelism is inferred.
    bool independentVelocityQueries=false;
    Integrator integrator=Integrator::ClosedRk4;
};
struct Report
{
    std::size_t tracedPoints=0, changedFirstJetPoints=0, acceptedSubsteps=0;
    double maxTraceRatio=0.0;
    // Difference between Dx(phi_y) and Dy(phi_x), scaled by hx*hy.
    // Diagnostic only; never used to relax a geometry/transport acceptance.
    double maxScaledMixedDisagreement=0.0;
};
struct Candidate
{
    LAMP::Geometry::ActiveQuarticHermitePatch patch;
    std::unique_ptr<LevelSetHermiteGeometry::Snapshot> geometry;
    // Owns the callbacks/domain, not objects borrowed by callback captures.
    // Exactly the numerical map used for this candidate's first-jet pullback.
    std::function<LevelSetBoundaryTrace::JetResult(Vector2D)> backward;
    // Departure from normalized physical time theta to the OLD time, with
    // the same velocity interpolation and selected integrator as backward.
    std::function<LevelSetBoundaryTrace::JetResult(double,Vector2D)> backwardAt;
    bool independentVelocityQueries=false;
    Report report;
    // Conditional bound for every successful trace under the declared
    // half-cell CFL. Allows a phase query to exclude interface support
    // without evaluating trajectories in a certified pure swept region.
    double maximumDisplacement=0;
};

// Complete uniform native leaf domain, at least two cells per direction.
// Transports independent corner (phi,dx,dy), shared edge (phi,dn), and a center
// scalar prediction. Mixed derivatives close from the transported skeleton's
// two five-point one-dimensional derivatives, at actual physical abscissae.
// No two-dimensional scalar/map finite differences or velocity Hessians.
// RZ axis parity is checked, never repaired after constructing a candidate.
// Zero first-jet change retains the old mixed derivatives exactly; nonzero
// steps close ALL regular vertices, including stationary points next to flow.
//
// This is the GALS PREDICTOR, not a conservative time-step acceptance. The
// caller still supplies FV phase transport and admissible interior-bubble
// closure before publishing InterfaceState. No physical inflow, contacts,
// sectors or mixed-level AMR are enabled by this initial uniform entry.
// Boundary hits reject explicitly; old inputs are never modified.
Candidate predictUniform(const QuadtreeMesh& mesh,
    const LAMP::Geometry::ActiveQuarticHermitePatch& previous,double dt,
    const LevelSetBoundaryTrace::DifferentiableVelocity& current,
    const LevelSetBoundaryTrace::DifferentiableVelocity& old,const Options& options={});
}
