#pragma once

#include <array>
#include <cstddef>
#include <memory>
#include <vector>
#include "Core/Base/Vector2D.h"
#include "Core/Geometry/ActiveHermitePatch.h"
#include "Core/Geometry/ActiveQuarticHermitePatch.h"
#include "Core/Geometry/ActiveHermiteSensitivity.h"
#include "Core/Geometry/CellChartPoint.h"
#include "Core/Mesh/UniformLeafGrid.h"

class QuadtreeMesh;
class LevelSetMappedSegment;
class LevelSetMappedInterface;
class LevelSetInterfacePanels;
template<typename T> class VertField;
namespace Quantum::Geometry {struct PolynomialAffineCoordinate;}

// Shared geometry, independent of material, force, solver and reference shape.
// Result structures contain no implementation types; Vector2D is not promised
// to be trivially copyable, so callers must not compare structure padding.
namespace LevelSetHermiteGeometry
{
// Selects the retained scalar itself, not a separate curvature estimator.
enum class Reconstruction { Cubic, EnrichedQuartic };

// Physical constant normal-speed transfer measure on an ascending axis-aligned
// face over [begin,end]. Exact input algebra in the strict geometry TU, then
// one floating representation (and pi for RZ). Independent of phase/BC data;
// does not evaluate or extrapolate a scalar. XY area / true RZ volume.
long double constantNormalFaceTransferVolume(Vector2D a,Vector2D b,double normalSpeed,
    double begin,double end,bool axisymmetric);

struct Piece
{
    Vector2D begin, end;
    bool negative = false;
};

struct ParametricPiece
{
    // Actual endpoints of one Cartesian subsegment, with a retained chart
    // orientation chosen for endpoint resolution. They define
    // p(t)=anchorBegin+t*(anchorEnd-anchorBegin); do not materialize global
    // root points before taking differences or integrating a small phase.
    Vector2D anchorBegin, anchorEnd;
    // Traversal follows the requested segment; tBegin may exceed tEnd.
    double tBegin = 0.0, tEnd = 1.0;
    bool negative = false;
};

struct SegmentRootIncidence
{
    int cellI = -1, cellJ = -1;
    // The original represented canonical subsegment, not rounded root points.
    // The root is anchorBegin+parameter*(anchorEnd-anchorBegin).
    Vector2D anchorBegin, anchorEnd;
    // Original cell-local restriction and its floating isolated-root point.
    Vector2D localBegin, localEnd, localPoint;
    double parameter = 0.0;
    // A cell meeting only this subsegment endpoint has no positive-length
    // restriction: localBegin=localEnd=localPoint, parameter is still 0 or 1
    // in the ORIGINAL physical anchor chart above.
    bool pointOnly = false;
};

enum class SegmentPartitionPolicy
{
    RequireRepresentableGridPieces,
    // In the phase-piece API, a grid-cut sliver lacking an interior point is kept
    // only when enclosing Bernstein bounds certify its ENTIRE closed phase.
    // An unresolved sliver containing/touching a zero still rejects.
    // A root-inventory query may instead isolate the complete polynomial in
    // an active closed cell containing both distinct physical/local anchors;
    // it reports endpoint incidences without claiming a physical phase sample.
    CertifyUnresolvedGridPieces
};

struct GridEdgeZeros
{
    std::vector<Vector2D> points; // All isolated roots, ordered from begin to end.
    bool coincident = false; // Entire edge is zero; points is then empty.
    std::vector<double> parameters; // Same roots in the oriented edge's [0,1] chart.
};

struct Frame
{
    Vector2D point, normal;
    double curvature = 0.0; // +div(n), negative-to-positive normal.
    double levelSetResidual = 0.0; // |phi|/|grad phi|, physical length.
    double stationarityResidual = 0.0; // Tangential projection residual, length.
    int iterations = 0, attemptedSeeds = 0, convergedSeeds = 0;
    // Explicit local-chart queries retain their actual represented position
    // here. point is its rounded global display coordinate, which may lose a
    // sub-ULP offset at a contact. Never recover this chart by subtracting point.
    Vector2D cellLocalPoint;
    int cellI = -1, cellJ = -1;
    bool hasCellLocalPoint = false;
    // Explicit endpoint-relative queries use this alternative provenance.
    // The two flags are mutually exclusive; cellLocalPoint is not a rounded
    // replacement for a chart point whose upper-end tail would be lost.
    Quantum::Geometry::CellChartPoint cellChartPoint;
    bool hasCellChartPoint = false;
    // Native provenance uses this actual leaf index, with i/j both -1.
    // Its chart is evaluated using that leaf's actual dimensions.
    int nativeCellIndex = -1;
};

// Read either explicit provenance representation, rejecting absent or
// conflicting provenance. Does not reconstruct a chart from frame.point.
Quantum::Geometry::CellChartPoint cellChartPoint(const Frame& frame);

struct BoundaryProjection
{
    Frame frame;
    bool endpoint = false;
    // Endpoint: (frame.point-query).t_in, where t_in points into the domain
    // along the curve. Interior candidate: zero (no one-sided direction).
    double inwardDerivative = 0.0;
    // Endpoint: max(0,-inwardDerivative); interior: stationarityResidual.
    // Both have units of length. frame.stationarityResidual is NOT replaced.
    double kktResidual = 0.0;
};

struct PatchEndpoint
{
    Vector2D point;
    int cellI = -1, cellJ = -1; // Incident active Cartesian cell containing the curve.
    int inwardTangentSign = 0; // Exactly +/-1 multiplying t=(-ny,nx) into that cell.
};

// Owns a retained snapshot: subsequent mesh/phi changes do not update it.
// The legacy constructors support complete uniform Cartesian layouts with
// regular-resolved active cells. The explicit native-leaf factory has the
// narrower query support documented below. The
// shared Hermite scalar is C1, not C2/C3, except explicitly declared pure-phase
// corner-sector seams, which retain C0 traces and one-sided jets. Ambiguous
// topology is rejected; such seams may not hide an interior zero.
// Geometric computations require nearest rounding and gradual underflow;
// unsupported environments are rejected, never changed on the caller's behalf.
// Enabled FP traps are outside the supported contract. The implementation is
// out-of-line so callers do not instantiate its exact-arithmetic machinery.
class Snapshot
{
public:
    struct ScalarJet2
    {
        double value = 0.0;
        double dx = 0.0, dy = 0.0;
        double dxx = 0.0, dxy = 0.0, dyy = 0.0;
    };

    // RZ opt-in imposes an even scalar extension at the actual r=0 boundary.
    // Nodal values are unchanged; boundary derivative stencils use symmetry.
    // RootOnly preserves the original complete root-grid restriction. Explicit
    // UniformLeaves accepts a complete equal-level leaf grid with its actual
    // labels/ancestry, not hanging-node AMR or active-fluid scalar support.
    Snapshot(const QuadtreeMesh& mesh, const VertField<double>& phi, bool radialAxisSymmetry = false,
        CartesianGridScope scope = CartesianGridScope::RootOnly,
        Reconstruction reconstruction = Reconstruction::Cubic);
    // EnrichedQuartic reuses the Q4 data extension and the common geometry
    // engine with 25 shared cell DOFs. Requires complete 7/8-node data stencils
    // (at least 7 cells per direction), never falls back to Cubic. No active
    // patches, hanging-node constraints, shared-data sensitivity pullback or
    // separate reconstructed-D3 field are supplied by this selection.
    // Mapped segments and the explicitly typed 25-column
    // LevelSetMappedInterface::integrateQuartic use this same scalar.
    Reconstruction reconstruction() const;
    // Retained active-cell patch with caller-supplied shared nodal derivatives.
    // Inactive payload is not geometry. Queries require an active cell; at an
    // EXACT Cartesian edge/vertex an incident active cell may supply the jet.
    // Nearby points are never moved across an inactive boundary. No solid
    // extension, contact-angle law or internal endpoint projection is inferred.
    Snapshot(const QuadtreeMesh& mesh, const Quantum::Geometry::ActiveHermitePatch& patch,
        bool radialAxisSymmetry = false, CartesianGridScope scope = CartesianGridScope::RootOnly);
    // Explicit native rectangular 2:1 layout. Uses the same Hermite scalar
    // with coarse-edge constraints at hanging vertices and declared sectors.
    // Currently supports value/negative/closed-cell sign, scalarJetInNativeCell,
    // native metadata/activity coverage, LevelSetMappedSegment (including root refinement), and
    // LevelSetMappedInterface, automatic NativeGraphChain discovery and the
    // explicit native-chart differential queries below. Logical i/j/global
    // spacing, global jet selection, projection and legacy full panels reject.
    // This factory does not select a new level-set PDE or enable App AMR.
    static Snapshot fromNativeLeaves(const QuadtreeMesh& mesh,
        const Quantum::Geometry::ActiveHermitePatch& patch, bool radialAxisSymmetry = false);
    // Same native geometry/chain/projection consumers, retaining the explicit
    // shared 25-DOF scalar instead. This factory does NOT generate Q4 data or
    // corner-sector extensions (declared replacements are validated).
    // Reconstruction reports EnrichedQuartic for this scalar family; the
    // caller, not the factory, supplies every datum. Use integrateQuartic,
    // nativeQuarticCellSensitivity and nativeQuarticNegativePhaseVariation.
    // Old 16-column APIs reject this source instead of converting.
    static Snapshot fromNativeLeaves(const QuadtreeMesh& mesh,
        const Quantum::Geometry::ActiveQuarticHermitePatch& patch, bool radialAxisSymmetry = false);
    // Pure refinement of the SAME retained Q4 scalar, including cut cells.
    // Full active coverage, unchanged domain, no sectors or coarsening. Every
    // target leaf must lie inside one source leaf; all 25 restricted Hermite
    // data are checked by exact identity. Owns its data after source destruction.
    // No M transfer, Interface state/publication, transport or accepted regrid.
    static Snapshot refineNativeQuartic(const QuadtreeMesh& targetMesh,const Snapshot& source);
    // Temporary local-query source. Includes requested cells and recursive
    // coarse constraint owners, preserving the ORIGINAL vertex/edge/center
    // independent IDs. Unselected leaves are inactive; this is not a complete
    // transported state. No corner-sector subsets or extrapolation.
    static Snapshot fromNativeQuarticQueryCells(const QuadtreeMesh& mesh,
        const Quantum::Geometry::ActiveQuarticHermitePatch& patch,const std::vector<int>& cells,
        bool radialAxisSymmetry = false);
    bool usesNativeLeaves() const;
    std::size_t nativeCellCount() const;
    bool nativeCellActive(int index) const;
    std::array<Vector2D, 2> nativeCellBounds(int index) const;
    // Existing outward-Bernstein certificate: one strict sign of the selected
    // scalar partial on the complete zero set in this active native leaf.
    // Axis 0/1 denotes x/y. False is inconclusive, including a turning point;
    // no nodal/curve samples replace the complete scalar predicate. Optional
    // local bounds restrict it to one retained native chart's closed subbox,
    // preserving the original-cell derivative and resolution scale.
    bool nativeZeroSetPartialResolved(int index, int axis,
        Vector2D localLower={0,0},Vector2D localUpper={1,1}) const;
    // Diagnostic only: the 25 shared tensor Q4 coefficients in one active
    // native quartic leaf's own [0,1]^2 chart (tensor index order). Native
    // quartic sources only; rounded display data, never an exact replacement.
    std::array<double, 25> nativeQuarticLocalDofs(int index) const;
    // Same authoritative homogeneous native Q4, composed with correlated
    // physical forms. Outside this leaf it bounds only that polynomial's
    // mathematical extension: the caller must prove actual donor membership.
    // Full coordinate type/API: Core/Geometry/PolynomialVelocitySlab.h.
    Quantum::Geometry::PolynomialAffineCoordinate nativeQuarticAffineEnclosure(int index,
        const std::array<Quantum::Geometry::PolynomialAffineCoordinate,2>& physicalCoordinates) const;
    struct AffineFaceTimeMeasure
    {
        long double total=0, liquidFraction=0;
    };
    // Integral on an ascending axis-aligned segment x [begin,end] of the
    // SAME native Q4 translated by constant velocity from referenceTime.
    // Exact polynomial identities certify one affine scalar on the complete
    // swept active box. Exact half-plane predicates/events and rational
    // polygon moments precede the final floating representation. Liquid < 0.
    // Includes true RZ metric; no exterior history or tolerance plane fit.
    AffineFaceTimeMeasure integrateNativeAffineTranslatedFace(Vector2D a, Vector2D b,
        Vector2D velocity, double referenceTime, double begin, double end) const;
    struct NormalTranslationSweep
    {
        Vector2D innerLower,innerUpper,outerLower,outerUpper;
        long double totalVolume=0,coordinateVolumeBand=0;
    };
    // Exact expansion of a - speed*(time-referenceTime) in a strict TU.
    // Inner/outer DOUBLE rectangles bracket the exact normal-translation
    // sweep; their volume difference bounds coordinate uncertainty only.
    // Nonzero normal speed, ascending axis-aligned face, full active native
    // Q4 support, and a resolved inner rectangle are required. No exterior
    // extension. RZ is an explicit metric and requires a horizontal face.
    NormalTranslationSweep nativeNormalTranslationSweep(Vector2D a,Vector2D b,
        double normalSpeed,double referenceTime,double begin,double end,bool axisymmetric) const;
    struct QuarticCenterInterval
    {
        double lower=0,upper=0; // Absolute center scalar values, not increments.
        int partialAxis=-1,partialSign=0;
        double minimumPhysicalDerivative=0;
        double maximumDisplacementPerCenter=0; // max psi/|partial phi| on the possible zero-set tube.
    };
    // Conservative admissible interval for varying ONLY the independent cell
    // center. Exact shared-data Bernstein differences certify one signed
    // Cartesian partial wherever a zero can occur, throughout the interval.
    // Outward subdivision covers every intermediate center; pure boxes need
    // no derivative bound. All face values/first traces stay fixed.
    // |delta center| max(psi/|partial phi|) bounds the
    // displacement along that direction by maxDisplacementCells*min(hx,hy).
    // This excludes critical zeros/new interior loops; it is not a curvature
    // error bound or a maximal feasibility interval. Cells
    // lacking this sufficient certificate explicitly request refinement.
    enum class QuarticCenterDirection { Both, Decrease, Increase };
    // Both preserves the symmetric search. Directed searches keep the
    // predicted center as one endpoint and use the SAME signed-partial and
    // displacement proof over the entire one-sided interval.
    QuarticCenterInterval nativeQuarticCenterInterval(int index,double maxDisplacementCells,
        double retainedDerivativeFraction=0.25,QuarticCenterDirection direction=QuarticCenterDirection::Both) const;
    struct QuarticChangeCertificate
    {
        bool resolved=false;
        double maxDisplacementCells=0;
        double minimumPhysicalPartial=0;
        int visitedBoxes=0;
    };
    // Certify the straight scalar homotopy between two same-cell native Q4
    // fields. Outward Bernstein boxes cover every intermediate zero; at least
    // one physical partial is nonzero there and |delta phi|/|grad phi| bounds
    // the normal homotopy motion. Shared DOF changes are allowed. This is a
    // sufficient regularity/displacement certificate, not a curvature bound.
    // A false result is inconclusive. Pure boxes need no gradient bound.
    QuarticChangeCertificate nativeQuarticChangeCertificate(int index,const Snapshot& other,
        double maxDisplacementCells) const;
    // Strict sign on a physical closed subbox of this native Q4 leaf, using
    // outward chart conversion and authoritative Bernstein controls. Zero
    // means inconclusive; degenerate point/line boxes are allowed.
    int nativeQuarticBoxSign(int index,Vector2D lower,Vector2D upper) const;
    // Exact constrained polynomial AND physical-box equality for two active
    // native Q4 leaves (including sectors/hanging data). No sampled tolerance.
    // A true result permits reusing that leaf's local chart after reindexing;
    // it says nothing about the remainder of either mesh or the physical state.
    // Other reconstruction families and inactive/invalid indices reject.
    bool identicalNativeQuarticCell(int index, const Snapshot& other, int otherIndex) const;
    // Owned leaf-search index; returns active AND inactive closed-box candidates.
    // Callers must check actual overlap and activity, not sample solid payload.
    std::vector<int> nativeCellCandidates(Vector2D lower, Vector2D upper) const;
    ScalarJet2 scalarJetInNativeCell(int index, Vector2D point) const;
    // Same deterministic active-cell selection as value(); no projection or
    // averaging across a cell boundary. Native scalar derivatives, including
    // at non-interface points where a normalized geometric frame is undefined.
    ScalarJet2 scalarJetNative(Vector2D point) const;
    // Floating cell-cache jet; homogeneous value/phase/root queries remain
    // the scalar authority, just as in the original uniform chart evaluator.
    ScalarJet2 scalarJetAtNativeCellChart(int index,const Quantum::Geometry::CellChartPoint& point) const;
    // Same scaled cell evaluation and original differential formulas as the
    // uniform engine. Curvature is meridional div(n), not RZ pole completion.
    // At a shared face D2/D3 need not agree: retain the caller's actual leaf.
    // This is a jet at a given point, NOT projection or a curvature fit.
    Frame jetAtNativeCellChart(int index,const Quantum::Geometry::CellChartPoint& point) const;
    double curvatureTangentialDerivativeAtNativeCellChart(int index,
        const Quantum::Geometry::CellChartPoint& point) const;
    // C^T times the 16 formal physical corner-jet sensitivities (SW/SE/NE/NW).
    // Returns unique regular native vertices in ascending order, all four
    // physical components retained. Pin/axis/free-gradient constraints remain
    // the caller's responsibility. Corner-sector cells reject by default;
    // HoldFixed omits ONLY the substituted corner jet, not the whole cell.
    // Sector data must not masquerade as shared regular-vertex unknowns.
    std::vector<Quantum::Geometry::ActiveHermiteNodalSensitivity> nativeCellSensitivity(
        int index, const std::array<long double, 16>& coefficients,
        Quantum::Geometry::CornerSectorVariation sectors =
            Quantum::Geometry::CornerSectorVariation::RequireAbsent) const;
    std::vector<Quantum::Geometry::ActiveQuarticHermiteSensitivity> nativeQuarticCellSensitivity(
        int index,const std::array<long double,25>& coefficients,
        Quantum::Geometry::CornerSectorVariation sectors =
            Quantum::Geometry::CornerSectorVariation::RequireAbsent) const;
    // Independent-data derivative of w.[phi,phi_x,phi_y,phi_xx,phi_xy,phi_yy]
    // at a retained native chart. Physical derivatives and the same Q4 basis;
    // no fitting, curvature model or regularization is introduced here.
    std::vector<Quantum::Geometry::ActiveQuarticHermiteSensitivity> nativeQuarticDifferentialSensitivity(
        int index,const Quantum::Geometry::CellChartPoint& point,
        const std::array<long double,6>& weights) const;
    ~Snapshot();
    Snapshot(Snapshot&&) noexcept;
    Snapshot& operator=(Snapshot&&) noexcept;
    Snapshot(const Snapshot&) = delete;
    Snapshot& operator=(const Snapshot&) = delete;

    // Metadata copies do not require a floating-point environment check.
    // All queries, including metadata, reject a moved-from snapshot.
    int nx() const;
    int ny() const;
    bool radialAxisSymmetry() const;
    bool restrictedDomain() const;
    bool hasCornerSectors() const;
    bool cellActive(int i, int j) const; // Cartesian cell coordinates, not native indices; bounds checked.
    Vector2D minCoord() const;
    Vector2D maxCoord() const;
    Vector2D cellSize() const;
    std::size_t candidateCellCount() const;
    std::size_t seedCount() const;
    std::size_t vertexCount() const;

    double value(Vector2D point) const;
    bool negative(Vector2D point) const;
    // Sufficient strict sign certificate on the entire CLOSED selected cell
    // containing point: -1 or +1 from the existing authoritative Bernstein
    // enclosure. Zero means NOT certified, not "phi is zero"; do not infer
    // either absence or presence of roots from zero. This includes cell edges
    // and corners and can certify wall/half-link separation without ignoring
    // same-phase tangencies. Uses the same FP/query contract as negative().
    int closedCellSign(Vector2D point) const;
    // Sufficient strict sign on the closed canonical cell EXCEPT explicitly
    // allowed zero corners. Bits 0..3 are SW, SE, NE, NW. Uses the signs of
    // exact shared-DOF Bernstein controls, not sampled phases. Zero remains
    // inconclusive. In particular no zero edge or hidden zero is exempted.
    int closedCellSignExceptVertices(Vector2D point, unsigned allowedCornerMask) const;
    // Exact signs of the authoritative shared hx*phi_x, hy*phi_y nodal DOFs.
    // Requires one exact Cartesian vertex with exactly zero nodal phi; no
    // near-zero normal or coordinate snapping defines topology here. A zero
    // derivative stays zero (the caller must not infer a transverse crossing).
    std::array<int, 2> zeroVertexGradientSigns(Vector2D point) const;
    // Incident-cell version for declared corner sectors; uses this cell's
    // authoritative DOFs, not the base cap jet from another incident sector.
    std::array<int, 2> zeroVertexGradientSignsInCell(int i, int j, Vector2D point) const;
    // Sufficient certificate that the entire CLOSED canonical cell is
    // strictly outside [-halfWidth,+halfWidth] in phi units. False means
    // uncertified, not necessarily that the band actually intersects it.
    bool closedCellOutsideBand(Vector2D point, double halfWidth) const;
    // Sufficient convex-hull certificate that the entire selected CLOSED
    // active cell is in [-halfWidth,+halfWidth]. False is inconclusive.
    // Used to verify a transport extension equals physical velocity around
    // the whole zero-curve reconstruction support, not just sampled nodes.
    bool closedCellInsideBand(Vector2D point, double halfWidth) const;
    enum class BandRequirement { EntireCellInside, EntireCellOutside };
    // Native-cell counterpart of an ordered short-circuit band check. Returns
    // the first failing POSITION, or cellIndices.size() if all satisfy the
    // original closed-cell predicate. Earlier query errors propagate; errors
    // after the first violation are not visited by the logical ordered check.
    // Implementation may query immutable cells concurrently, with independent
    // results and ordered failure reduction. No sampled-point substitute.
    std::size_t firstNativeBandViolation(const std::vector<int>& cellIndices,
        double halfWidth, BandRequirement requirement) const;
    Frame jet(Vector2D point) const;
    // Same shared scalar, using the specified active closed cell's own
    // gradient/Hessian. No projection, Q4 refit or cross-cell averaging.
    Frame jetInCell(int i, int j, Vector2D point) const;
    // Multistart local stationary projection, NOT a certified globally
    // closest point or distance oracle. Does not clip queries to the domain.
    Frame projectedFrame(Vector2D point) const;
    // Explicit local-provenance version for a represented PHYSICAL query,
    // e.g. an electric face center. Reuses the same seeds, Newton iterations,
    // backtracking and residual limits. At convergence it retains the actual
    // selected cell/sample parameters and produces its jetAtCellLocal Frame.
    // That local chart's rounded display point may differ from the legacy
    // global iterate by roundoff; residuals are rechecked before acceptance.
    // This is not a way to recover an old sub-ULP root from its display point.
    // No declared boundary/axis endpoint candidates or root snapping are
    // added; callers needing endpoint projection must supply a separate
    // authoritative endpoint path. The legacy projectedFrame is unchanged.
    Frame projectedCellLocalFrame(Vector2D point) const;
    // Adds declared noncorner, nontangential physical-wall zero endpoints to
    // the SAME local multistart stationary search. Endpoints must be exact
    // scalar zeros or exact canonical intersectGridEdge root representatives;
    // a small sampled phi alone is never accepted. A one-sided KKT condition
    // is required at the endpoint. No endpoint clipping/contact-angle law is
    // inferred; this remains a local search, NOT a global nearest certificate.
    // An empty endpoint list preserves projectedFrame's result and failures.
    // Declared endpoints remain limited to the outer rectangular domain;
    // internal active/inactive contacts are not supported by this API.
    BoundaryProjection boundaryProjectedFrame(Vector2D query,
        const std::vector<Vector2D>& boundaryEndpoints) const;
    // Restricted-patch counterpart: endpoints must be authoritative zeros on
    // actual active-domain boundary edges (including axis/corner incidence).
    // Their explicit incident cell supplies the one-sided jet and +/- tangent
    // direction; a tangent along a cell side is allowed, an outward one is not.
    // Reuses the same local stationary search and endpoint KKT selection.
    // Does not infer a contact law, global nearest point or RZ pole curvature.
    BoundaryProjection patchBoundaryProjectedFrame(Vector2D query,
        const std::vector<PatchEndpoint>& boundaryEndpoints) const;
    // t.grad(+div(n)), t=(-ny,nx), from the same canonical cell as jet.
    // Evaluates the supplied point without projection or cross-cell averaging.
    double curvatureTangentialDerivative(Vector2D point) const;
    // Same RAW Hermite curvature derivative, in physical length^-2, using
    // exactly the active closed cell/local point supplied to jetAtCellLocal.
    // D1/D2/D3 all come from that cell's retained polynomial, including its
    // declared one-sided corner sector. Pass a retained Frame's cellI/J and
    // cellLocalPoint directly; no global-display round trip is performed.
    // The C1 scalar is not C2/C3 (a sector seam is only C0), so this derivative
    // can jump across cells/sectors. It is NOT the five-point reconstructed-D3
    // discretization below. No neighboring/solid values, projection or fitting
    // are introduced. This is the meridional +div(n) derivative, not an added
    // RZ azimuthal-curvature term. The original local jet gradient, activity,
    // [0,1]^2, finite-result, FP and moved-from guards remain in force.
    double curvatureTangentialDerivativeAtCellLocal(int i, int j, Vector2D localPoint) const;
    // Same raw, one-sided CellData D3 as the local API, evaluated with
    // endpoint-reflected Hermite bases. This is not reconstructed C3 data.
    double curvatureTangentialDerivativeAtCellChart(int i, int j,
        const Quantum::Geometry::CellChartPoint& point) const;
    // Same-point gradient/Hessian as jet, with D3(phi) reconstructed from
    // five-point nodal differences and bilinear interpolation. Avoids the
    // one-sided C1 polynomial's third-derivative jump; second-order consistent
    // for smooth resolved phi with nonzero gradient. Not the exact derivative
    // of the piecewise Hermite curvature, nor a globally C1 curvature promise.
    // Leaves scalar values, phases, projection, normals and curvature intact.
    // EnrichedQuartic rejects: use its same-cell D3, not another fitted field.
    // ActiveHermitePatch also rejects: the surrounding
    // five-point derivative stencil has no implicit inactive-domain extension.
    double reconstructedCurvatureTangentialDerivative(Vector2D point) const;
    // Pure retained-domain predicate, independent of phi, phase roots and
    // representable subsegment midpoints. True exactly when the entire CLOSED
    // segment is covered by active CLOSED canonical cells. A segment along a
    // grid edge needs an active incident cell throughout; endpoint activity
    // alone is insufficient. Zero length checks the point's activity.
    // Finite/outer-rectangle and FP query guards still apply. Unsupported exact
    // arithmetic scales throw rather than dropping a small uncovered interval.
    bool segmentInsideActiveDomain(Vector2D a, Vector2D b) const;
    // Sufficient strict sign on the entire CLOSED segment, including endpoints.
    // Uses exact active-cell interval coverage and authoritative Bernstein
    // bounds on the enclosing cell/segment boxes. No roots or subsegment midpoints are
    // required. Returns +/-1 only with complete same-sign coverage; 0 means
    // inconclusive (including an activity gap), never "the segment is zero".
    int certifiedSegmentSign(Vector2D a, Vector2D b) const;
    // Every positive-length subsegment must lie in an active closed cell.
    // Traversal through inactive cells is rejected, never skipped or clipped.
    std::vector<Piece> partitionSegment(Vector2D a, Vector2D b) const;
    // Same exact restrictions and phase certificates, with explicit
    // endpoint-relative root isolation. Retains every positive parameter
    // interval even if its global image has no representable interior point.
    // Records follow input traversal. Reversal preserves each anchor chart,
    // reverses records and swaps tBegin/tEnd (which may therefore decrease),
    // never computes 1-t. The original grid-cut representation/activity gates
    // still apply by default; the explicit certified-grid-piece policy only
    // exempts midpoint representability after a whole-piece strict sign proof.
    // Collapsed grid-cut parameters/anchors and unresolved zeros still reject.
    // No arbitrary near-endpoint/root accuracy is promised.
    std::vector<ParametricPiece> partitionSegmentParametric(Vector2D a, Vector2D b,
        SegmentPartitionPolicy policy = SegmentPartitionPolicy::RequireRepresentableGridPieces) const;
    // All isolated zeros of the represented physical segment, retaining the
    // canonical subsegment and SAME-cell local restriction/root provenance.
    // Includes active cells on both sides of a gridline and cells meeting
    // only a segment/grid-cut endpoint. Repeated incidences are NOT merged.
    // Records follow input traversal; reversal only reverses records, leaving
    // anchor charts, local charts and parameters unchanged (no 1-parameter).
    // Uses endpoint-relative isolation. Coincident positive-length zero
    // restrictions, collapsed local roots and inactive traversal reject.
    // The explicit certified-grid-piece policy may use a complete cell-local
    // root isolation for a grid sliver without a physical midpoint. Both
    // anchors must remain distinct and lie in the same active closed cell;
    // strict policy retains the physical midpoint requirement.
    // Root representatives and the existing rounded canonical embedding are
    // floating geometry, NOT an exact/watertight intersection certificate.
    std::vector<SegmentRootIncidence> intersectSegmentWithProvenance(Vector2D a, Vector2D b,
        SegmentPartitionPolicy policy = SegmentPartitionPolicy::CertifyUnresolvedGridPieces) const;
    // Full zero inventory on ONE complete native Cartesian edge. Endpoints
    // must be exact adjacent grid vertices. Unlike phase partitioning, this
    // retains endpoint roots and reports coincident edges. Unresolved roots
    // still throw; no near-zero snapping or root deletion is performed.
    GridEdgeZeros intersectGridEdge(Vector2D a, Vector2D b) const;
    // Isolated zeros on a nonzero segment in ONE cell's [0,1]^2 chart.
    // No round-trip through global coordinates; coincident segments reject.
    std::vector<Vector2D> intersectCellSegment(int i, int j, Vector2D localA, Vector2D localB) const;
    Frame jetAtCellLocal(int i, int j, Vector2D localPoint) const;
    Frame jetAtCellChart(int i, int j, const Quantum::Geometry::CellChartPoint& point) const;
    // All roots of a complete Cartesian cell section, ordered along the
    // unfixed coordinate. Each half is isolated from its nearest endpoint;
    // only the identical exact midpoint root is shared between halves.
    // fixedCoordinate is retained verbatim. Coincident/unresolved sections
    // reject; no positive tail or root identity is rounded away.
    std::vector<Quantum::Geometry::CellChartPoint> intersectCellSection(int i, int j,
        bool coordinateIsX, Quantum::Geometry::CellChartCoordinate fixedCoordinate) const;
    // Raw value and physical derivatives from the SAME retained cell
    // evaluation as jetAtCellLocal, including its one-sided corner sector.
    // No projection, global-coordinate round trip, Q4 refit or averaging.
    // Unlike a normal/curvature Frame, this permits zero gradients. Cell and
    // local [0,1]^2 bounds, activity, finite results, FP and moved-from guards
    // are retained; this floating evaluation is not an exact sign certificate.
    ScalarJet2 scalarJetAtCellLocal(int i, int j, Vector2D localPoint) const;
    ScalarJet2 scalarJetAtCellChart(int i, int j,
        const Quantum::Geometry::CellChartPoint& point) const;
    double negativeFraction(Vector2D a, Vector2D b) const;

private:
    friend class ::LevelSetMappedSegment;
    friend class ::LevelSetMappedInterface;
    friend class ::LevelSetInterfacePanels;
    // Stronger than polynomial equality: only identical terminal/cache data
    // and homogeneous expansion sequences permit reusing computed GL8 values.
    bool identicalNativeQuarticRepresentation(int index, const Snapshot& other, int otherIndex) const;
    // Internal cooperation with the retained panel graph: it must filter
    // actual curve incidences before selecting a stationary/endpoint result.
    // This is the original bounded multistart inventory, not a new search.
    struct LocalStationaryCandidates
    {
        std::vector<Frame> frames;
        int attemptedSeeds = 0, convergedSeeds = 0;
    };
    LocalStationaryCandidates stationaryCellLocalFrames(Vector2D point) const;
    struct Impl;
    explicit Snapshot(std::unique_ptr<const Impl> impl);
    std::unique_ptr<const Impl> impl_;
    const Impl& checked() const;
};
}
