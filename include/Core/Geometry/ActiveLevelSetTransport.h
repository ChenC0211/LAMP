#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <vector>
#include "Core/Geometry/ActiveHermitePatch.h"
#include "Core/Geometry/ActiveQuarticHermitePatch.h"
#include "Core/Geometry/ActiveQuarticHermiteReconstruction.h"
#include "Core/Geometry/LevelSetBoundaryTrace.h"
#include "Core/Geometry/LevelSetHermiteGeometry.h"

class QuadtreeMesh;

namespace ActiveLevelSetTransport
{
struct FixedWallCornerExtension
{
    int cellIndex = -1, vertexIndex = -1, phaseSign = 0;
    bool wallHorizontal = false;
};
// Value at a resolved physical inflow event. Time is elapsed backwards from
// this advance's new time; the caller owns the absolute boundary history.
using InflowScalar=std::function<double(Boundary,Vector2D,double)>;
struct Options
{
    // Explicit scalar-amplitude velocity extension, not a distance claim:
    // chi=1 inside fullWidth, C2 taper to 0 at zeroWidth. Both time levels
    // use the SAME old scalar. Physical momentum velocities are not changed.
    double fullWidth = 0.0, zeroWidth = 0.0;
    bool radialAxisSymmetry = false;
    // Native jet port only. No exterior velocity or scalar sampling. Scalar
    // data at boundary events enter only the derivative supports they affect.
    InflowScalar inflowScalar;
    std::vector<int> fixedZeroVertices; // Actual native IDs; fixes phi, not angle.
    // Optional one-sided derivative/data support. The cubic entry uses 5x5;
    // the quartic entry uses the existing 7/8-point recursive extension.
    // Nodes admitted here use this mask; no extrapolation.
    std::vector<std::uint8_t> oneSidedSamples;
    // Optional complementary smooth extension chart. Nodes on the shared
    // seam use oneSidedSamples (the cap/base jet), not this side's derivative.
    std::vector<std::uint8_t> extensionSideSamples;
    // Explicit pure-phase corner extension: copy the evolved cap/base jet,
    // but give this cell a flat wall-tangential endpoint derivative. This
    // does NOT modify the base jet/contact angle. The scalar constructor
    // must certify C0 edge traces and phaseSign except the fixed zero corner.
    std::vector<FixedWallCornerExtension> fixedWallCornerExtensions;
    // Optional nonphysical, sequential workspace. Q4 support selection only;
    // exact topology/activity/charts are checked by the Core generator.
    std::shared_ptr<ActiveQuarticHermiteReconstruction::StencilCache> quarticStencils;
    // Native quartic only. Explicit execution contract, not a physical option:
    // both velocity callbacks and optional inflowScalar must allow independent concurrent read-only
    // queries (including any captured state). Default keeps callbacks on the
    // caller in first-use order. Per-trajectory RK2 callback order is unchanged.
    bool independentVelocityQueries = false;
    // Use the supplied physical velocity everywhere in the active domain.
    // Explicit opt-in for conservative transport sharing the same trajectory;
    // fullWidth/zeroWidth are unused. Existing callers retain scalar tapering.
    bool physicalVelocity = false;
};

template<class PatchType> struct TransportCandidate
{
    using Patch = PatchType;
    Patch patch;
    std::unique_ptr<LevelSetHermiteGeometry::Snapshot> geometry;
    // Inverse of the SAME retained numerical backtrace used for patch. Owns
    // the old scalar, dt, widths and copies of both velocity callbacks. Any
    // objects borrowed by those callbacks must still outlive their use; an
    // owning immutable velocity snapshot is recommended. Geometry correction
    // after advance does not redefine this material trajectory.
    // Inflow characteristics have no old interior material identity and are
    // rejected by this inverse; scalar boundary history does not invent one.
    std::function<Vector2D(Vector2D)> forwardMaterialPoint;
    // Native entry counts independent regular nodes, not dependent hanging data.
    int activeNodes = 0, changedValues = 0, stationaryNodes = 0, taperedNodes = 0;
    double maxValueIncrement = 0.0, maxTraceRatio = 0.0;
    // Observational counts for native Q4 samples only; not acceptance criteria
    // or default file output. Nodal value and material-map queries are separate.
    std::size_t nativeStencilSampleRequests = 0, nativeStencilSampleEvaluations = 0;
};
struct Candidate : TransportCandidate<LAMP::Geometry::ActiveHermitePatch> {};
// Wall-clock observations only, with nonoverlapping scopes inside the native
// quartic advance. Never used by a timestep, acceptance gate or field writer.
struct QuarticTimings
{
    double sourceGeometry=0, stencilPreparation=0, scalarSampling=0;
    double dataReconstruction=0, dataUpdate=0, candidateGeometry=0;
};
struct QuarticCandidate : TransportCandidate<LAMP::Geometry::ActiveQuarticHermitePatch>
{
    // Values below supplement the nodal diagnostics, not acceptance gates.
    int changedEdgeValues = 0, changedCenterValues = 0;
    std::size_t scalarIncrementEvaluations = 0;
    std::size_t scalarIncrementRequests = 0;
    ActiveQuarticHermiteReconstruction::Report reconstruction;
    QuarticTimings timings;
};

// Retains old DOFs and adds Q4(delta_phi) derivatives after the existing
// two-time-level midpoint backtrace. Thus zero displacement does not refit
// an already compatible contact jet. Inactive payload is never sampled.
// The new snapshot is constructed before returning; caller must additionally
// validate its topology/contacts, transport-band coverage and port separation
// before committing. No material, contact-angle, inflow or redistance law here.
Candidate advance(const QuadtreeMesh& mesh, const LAMP::Geometry::ActiveHermitePatch& previous,
    double dt, const LevelSetBoundaryTrace::Velocity& newVelocity,
    const LevelSetBoundaryTrace::Velocity& oldVelocity, const Options& options);
// Fixed mixed-leaf counterpart: same characteristic, scalar-increment Q4
// kernel and contact rules; only independent regular jets are advanced.
// 25 local query samples use the finest incident active-leaf spacing, not
// interpolated hanging payload or a replacement mesh. A sample-mask chart
// must admit every corner of every leaf overlapping the stencil interior.
// Both velocity callbacks must represent fixed, repeatable time-level fields
// throughout this call and subsequent use of forwardMaterialPoint. Within
// this call only, bitwise-identical Q4 sample coordinates reuse their scalar
// increment. No cache is retained by the candidate or material trajectory.
// Trajectories retain the conservative GLOBAL finest active-cell half-step
// bound. Runtime remeshing, material projection and physical inflow are not
// enabled by this entry; the caller still validates/corrects the candidate.
Candidate advanceNative(const QuadtreeMesh& mesh, const LAMP::Geometry::ActiveHermitePatch& previous,
    double dt, const LevelSetBoundaryTrace::Velocity& newVelocity,
    const LevelSetBoundaryTrace::Velocity& oldVelocity, const Options& options);
// Cubic counterpart of advanceQuarticNativeJet. Shares the retained RK2 map
// and old native jet; 5x5 physical supports differentiate its displacement.
// Axis map derivatives are one-sided (no even extension of radial motion).
// Supports hitting explicit scalar inflow reconstruct transported derivatives;
// untouched supports use the jet chain rule. Explicit fixed zeros require
// actual no-slip velocity in both layers. Retained flat wet-corner extensions
// require both complete active half-chart masks about their actual wall seam;
// the base jet is pulled back from the actual cap
// cell across the seam, then the existing pure wet-sector rule is applied.
// The base derivatives remain free; no contact-angle value is prescribed.
Candidate advanceNativeJet(const QuadtreeMesh& mesh, const LAMP::Geometry::ActiveHermitePatch& previous,
    double dt, const LevelSetBoundaryTrace::Velocity& newVelocity,
    const LevelSetBoundaryTrace::Velocity& oldVelocity, const Options& options);
// Same retained characteristic and GLOBAL finest-cell trace bound. Add the
// pointwise delta_phi at each independent regular vertex, maximal-edge
// midpoint and center; Q4 recursive extension supplies derivative increments.
// Thus a stationary point retains its value even if derivative supports cross
// moving data. Never refit the full scalar or prescribe a wall value here.
// The complete 7/8-point supports and one-sided charts are owned by the shared
// ActiveQuarticHermiteReconstruction generator. Hanging/fine-side payload is
// not advanced. Independent edge identity/order from previous is preserved.
// Missing midpoint/center scalar increments form a second ordered batch using
// the SAME backtrace and explicit read-only concurrency contract. Its starting
// scalar is reused from that backtrace, never queried again or approximated.
// A valid geometry is returned, NOT an accepted physical timestep: topology,
// band/port checks, local volume/phase-flux correction and material projection
// remain caller responsibilities. No runtime remeshing or new contact law.
QuarticCandidate advanceQuarticNative(const QuadtreeMesh& mesh,
    const LAMP::Geometry::ActiveQuarticHermitePatch& previous, double dt,
    const LevelSetBoundaryTrace::Velocity& newVelocity,
    const LevelSetBoundaryTrace::Velocity& oldVelocity, const Options& options);
// Closed-component migration port: shares the characteristic, pointwise value
// updates and candidate construction above. Derivatives instead pull back the
// old native scalar jet by the same numerical backtrace. The shared Q4 support
// differentiates backtrace DISPLACEMENT, not scalar increments. This is a
// numerical map-derivative approximation, not exact automatic differentiation.
// Radial map components use one-sided physical supports (radial displacement
// is odd, unlike the scalar); scalar axis parity is imposed on the final jet.
// Contact sectors/pinned contacts are not enabled by this port.
// With explicit inflow data, supports touching a physical boundary event
// reconstruct derivatives of the transported scalar (including its boundary
// history). Unaffected supports retain the interior jet pullback.
QuarticCandidate advanceQuarticNativeJet(const QuadtreeMesh& mesh,
    const LAMP::Geometry::ActiveQuarticHermitePatch& previous, double dt,
    const LevelSetBoundaryTrace::Velocity& newVelocity,
    const LevelSetBoundaryTrace::Velocity& oldVelocity, const Options& options);
}
