#pragma once
#include "Modules/Interface/Compatibility/InterfaceConservativePhase.h"
#include "Core/Geometry/PlanarCellCut.h"
#include "Core/Geometry/QuadraticCellCut.h"
#include "Core/Geometry/UnitPolynomialRoots.h"
#include <functional>
class InterfaceNormalFluxView;
struct InterfaceMotionVelocityLayers;

namespace InterfacePlicTransport
{
using CurvedDonors=std::vector<std::optional<Quantum::Geometry::QuadraticCellCut::Shape>>;
// Optional native-cell graph shapes replace only the corresponding resolved
// donor PLIC. Their positions are solved from this call's immutable inventory;
// an empty entry retains PLIC. Shapes are auxiliary flux geometry, not a second
// published physical interface. All candidates use the same shared-face bounds.
struct CellData { long double fraction=0; Vector2D normal{1,0};
    std::optional<Quantum::Geometry::QuadraticCellCut::Shape> curved; };
// A minority fraction below Core's PLIC inverse resolution has no reliable
// facet. Both polygon and space-time paths transport its unmodified cell
// average over each donor intersection. It needs no invented normal and is
// neither dropped nor rounded to a pure phase; resolved donors remain PLIC.
// Boundary data for actual exterior Cartesian ghost rectangles. The returned
// fraction/normal define a PLIC continuation, not a sampled point phase. A
// callback is required only when a nonzero swept piece reaches the exterior;
// it is synchronous, not retained, and must return the same immutable data for
// repeated queries in this interval. Wall/contact boundary laws are not inferred.
using ExteriorCell=std::function<CellData(const Quantum::Geometry::PlanarCellCut::Rectangle&)>;
// Explicit phase history on a whole native boundary face for this entire
// physical step: -1 liquid, +1 gas, 0 no port. The adjacent old donor cell
// must have exactly that pure inventory, keeping the interface away from the
// inlet layer. No time-varying mixture or undeclared exterior state is inferred.
using SinglePhasePorts=std::function<int(int nativeFace)>;
struct SweepReport
{
    std::size_t sweptFaces=0,correctedFaces=0,sweptPieces=0;
    long double maxRelativeVolumeCorrection=0,maxNormalScaleChange=0;
};
struct SpaceTimeOptions
{
    double relativeFluxTolerance=1e-5, mapToleranceCells=1e-4;
    int maxSpaceDepth=20,maxTimeDepth=10;
    // Opt-in promise: both normal-flux views support independent concurrent
    // queries. Velocity callbacks must also declare independentQueries.
    // Each native face keeps its serial quadrature and directional ledger.
    bool parallelFaces=false;
    // Resolve localized grazing intervals without expanding a whole binary
    // tree. The leaf budget equals the former depth-12 worst case; error and
    // map tolerances are unchanged. Both bounds apply per signed slice piece.
    std::size_t maxSpaceLeaves=4096;
};
struct SpaceTimeReport
{
    std::size_t tracedPoints=0,timeSlices=0,spaceLeaves=0,timeLeaves=0,rootEnclosures=0;
    std::size_t maxAcceptedSpaceLeaves=0;
    int maxAcceptedSpaceDepth=0;
    std::size_t constantPhaseChords=0; // Exact 0/1 throughout the CFL donor halo; no phase clipping.
    // Observations of one immutable half-face/time slice, not state or tolerances.
    std::size_t primitiveRequests=0,primitiveEvaluations=0;
    // Old donor boundary vertices are immutable for this entire spaceTime
    // call. Lazy reuse never spans a new inventory, mesh or physical step.
    std::size_t arrangementVertexRequests=0,arrangementVertexBuilds=0,arrangementEdgeRootQueries=0;
    std::size_t temporalRootEnclosures=0;
    std::size_t temporalGeometryCuts=0; // Characteristic/donor crossing seeds; not root certificates.
    std::size_t indivisibleTimeIntervals=0;
    long double indivisibleTimeVolumeBound=0; // Full gross-flow bound charged to temporal error.
    std::size_t temporalVertexSeeds=0,temporalTangencySeeds=0; // Before shared-time deduplication.
    double maxTemporalRootEnclosureWidth=0;
    long double estimatedTemporalError=0,estimatedSpatialError=0,maxNetFluxDifference=0;
    double maxAcceptedMapDeviationCells=0,maxAcceptedSpatialRelativeError=0,maxRootEnclosureWidth=0;
    std::size_t inflowTraceEvents=0;
    long double boundaryInflowVolume=0,boundaryInflowLiquid=0;
};
// Integrate the old reconstructed PLIC phase at the departure of each actual
// face/time point, against the SAME retained normal-flux primitive. Piecewise
// linear departure traces are clipped into disjoint donor intervals; spatial
// refinement compares phase integrals and midpoint map curvature. Positive
// time quadrature and retained roots of interpolated polynomial/separable
// densities retain gross flows. Separable amplitudes split at temporal zeros.
// Donor-plane/curve and cell crossings along half-face endpoint/midpoint
// characteristics also split time, preventing endpoint phase aliasing for
// resolved traces. Bilinear trace patches also seed arrangement vertices and
// interior quadratic tangencies, detecting phase pockets absent on those three
// trajectories. Simple events are polished on the same RK trace. Map refinement
// and nonlinear-event error control remain observational, not certified.
// Error estimates are refinement observations, not interval certificates for
// arbitrary velocity callbacks. Selected native fluid leaves, including mixed
// levels; characteristic CFL and map resolution use minimum selected spacing.
// Actual native subfaces share one flux; departure donor intervals use leaf
// bounds and positive-side ownership at aligned coarse/fine boundaries. Inactive
// and internal solid faces require identically zero density in both layers.
// RK stages and departure chords must stay in their actual closed union.
// Optional outer single-phase ports close resolved boundary events without
// sampling exterior velocity. Internal inflow/contact physics is not inferred.
std::vector<InterfaceConservativePhase::FaceFlux> spaceTime(const QuadtreeMesh& mesh,
    const NativeLeafDomain& domain,const InterfaceConservativePhase& inventory,
    const std::vector<Vector2D>& normals,const InterfaceStepTime& time,
    const InterfaceMotionVelocityLayers& velocities,const InterfaceNormalFluxView& previousFlux,
    const InterfaceNormalFluxView& currentFlux,const SpaceTimeOptions& options={},SpaceTimeReport* report=nullptr,
    const CurvedDonors& curved={},const SinglePhasePorts& ports={});
struct FaceSweep
{
    // Positive tangent coordinate on the actual native face (x on horizontal,
    // y on vertical). Pieces must exactly partition [0,1], without overlaps.
    // Adjacent inputs must share their endpoint displacement exactly. Normal
    // volume correction is per piece, not a certified common material map.
    double begin=0,end=1;
    long double totalVolume=0;
    Vector2D beginDisplacement,endDisplacement;
};
struct TraceReport
{
    struct RootUse
    {
        int nativeFace=-1,half=-1;
        bool current=false;
        Quantum::Geometry::UnitPolynomialRoots::Root root; // Original half coordinate and enclosure.
    };
    std::vector<RootUse> roots;
    std::size_t zeroTraceHalves=0,automaticInteriorCuts=0;
};
// Bridge from the existing frozen velocity/flux views. Exact half-face
// primitives are trapezoid-integrated between the two external time layers;
// point trajectories reuse Core's two-layer RK2 and exact-coordinate cache.
// Without explicit cuts, views with polynomial-root capability contribute all
// old/current roots in addition to {0,1/2,1}. Root enclosures remain in report;
// roots that cannot be represented in the shared face chart reject. Views
// without that capability retain the explicit two-half baseline (no claim of
// automatic sign isolation). Optional caller cuts must contain 0,1/2,1.
// Temporal reversal at sampled
// traces and physical inflow events reject and require the caller's substeps
// or boundary trajectory closure. No outside velocity is fabricated.
std::vector<std::vector<FaceSweep>> traceSweeps(const QuadtreeMesh& mesh,
    const NativeLeafDomain& domain,const InterfaceStepTime& time,
    const InterfaceMotionVelocityLayers& velocities,const InterfaceNormalFluxView& previousFlux,
    const InterfaceNormalFluxView& currentFlux,const std::vector<std::vector<double>>& partitions = {},
    TraceReport* report=nullptr);
// Spatially partitioned counterpart: each piece has one direction; opposite
// pieces on the SAME shared face retain their separate total/liquid volumes.
// No zero-net cancellation before upwind, CFL, limiting or momentum consumers.
// Empty per-face pieces mean zero transport. Temporal reversal still requires
// advanceSubsteps with a reconstruction from each intermediate inventory.
std::vector<InterfaceConservativePhase::FaceFlux> sweptSubfaces(
    const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const InterfaceConservativePhase& inventory,const std::vector<Vector2D>& normals,
    const std::vector<std::vector<FaceSweep>>& faceSweeps,
    const ExteriorCell& exterior = {},SweepReport* report=nullptr,const CurvedDonors& curved={});
// Shared, time-integrated native face volumes from the velocity layer, and
// native vertex displacements x-X_old(x) of its backward characteristics. Both inputs
// must describe ONE monotone interval. Endpoint trajectories are joined by
// straight segments; normal depth is scaled to match each prescribed face
// volume in the actual XY/RZ measure. A correction above 25% of |F| rejects.
// Full-vector corner donor transport is retained; RZ includes radial motion.
// Each endpoint must satisfy |dx|/hx+|dy|/hy<=1 after correction. Opposite normal
// endpoint motions, folded polygons and unresolved volume constraints reject;
// spatial counterflow needs face subdivision, temporal reversal needs substeps.
// This port does not infer velocity time layers or invent a divergence-free F:
// InterfaceConservativePhase::advance still checks GCL/CFL and shared bounds.
std::vector<InterfaceConservativePhase::FaceFlux> swept(
    const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const InterfaceConservativePhase& inventory,const std::vector<Vector2D>& normals,
    const std::vector<long double>& faceVolumes,const std::vector<Vector2D>& vertexDisplacements,
    const ExteriorCell& exterior = {},SweepReport* report=nullptr,const CurvedDonors& curved={});
// Complete fixed uniform rectangle. One spatially constant displacement during
// a monotone interval; XY can be diagonal, RZ permits axial displacement only.
// CFL |dx|/hx+|dy|/hy <= 1. Exact face preimages use the FULL displacement,
// partitioned by donor cells, so diagonal corner transport is not double spent.
// Normals use native cell indexing; pure-cell normals are not consumed.
// Exact constant-velocity specialization, using the same donor geometry.
std::vector<InterfaceConservativePhase::FaceFlux> translation(
    const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const InterfaceConservativePhase& inventory,const std::vector<Vector2D>& normals,
    Vector2D displacement,const ExteriorCell& exterior = {},const CurvedDonors& curved={});
}
