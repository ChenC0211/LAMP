#pragma once

#include <array>
#include <memory>
#include <optional>
#include <vector>
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Geometry/LevelSetMappedInterface.h"

class QuadtreeMesh;

// Retained, material-independent integration geometry for ONE meridional
// interface on a complete fixed Cartesian mesh. One regular graph segment per
// cut cell is supported.
// RootOnly is the default layout contract. Explicit UniformLeaves accepts the
// actual complete same-level leaf lattice matching the supplied Snapshot;
// this does not select an active-fluid subset or change endpoint rules.
// By default wall/vertex contact, along-edge zeros,
// multiple crossings on an edge, branches and multiple components are rejected. RZ with
// an explicitly radial-symmetric Snapshot also accepts one open chain whose
// only boundary contacts are two crossings of the r=0 axis.
// Explicit boundaryContacts opt into a separate open-chain path: declared
// noncorner wall endpoints may join an axis endpoint or another declaration.
// Exact grid vertices share one graph node. Strict same-sign inward nodal
// partials can exclude only that corner's local cell incidence; an otherwise
// empty cell must still be certified single-phase away from its zero corners.
// This opt-in also retains multiple resolved roots on a shared edge and
// interior vertex crossings, but still only one regular graph per cell and
// never along-edge zeros/branches. No global root is deleted by incidence tests.
// Declaration points must exactly equal retained roots; no snapping is used.
// Declarations identify geometry only; they impose no pinning/contact-angle law.
// With an explicit active-patch Snapshot, declarations may instead lie on an
// actual active/inactive native boundary (including a fitted corner). Only
// active cells participate; all roots/incidences on their boundary remain
// checked. The no-declaration closed-chain path retains its full-domain scope.
// An explicit corner-sector Snapshot uses original cell-local edge-root
// parameters for every declared-chain panel, including cross-edge panels.
// Frames, orientation and GL8 integration stay in that local chart; a rounded
// physical display point is not inverted to recover the curve coordinate.
//
// The original const-reference construction owns its results and retains no
// mesh/snapshot references. The explicit shared-source overload additionally
// owns the immutable source and cell-local graphs for new subarc queries.
// Arc lengths and weights are fixed GL8 approximations on the actual Hermite
// zero curve, NOT exact arc lengths and NOT lengths of a replacement polygon.
// In axisymmetric mode the curve interior must stay at r=x>0; physical measures
// use 2*pi*r*ds. The original arc-length fields and closed-loop behavior are unchanged.
class LevelSetInterfacePanels
{
public:
    enum class EndpointKind { None, Axis, DeclaredBoundary, SymmetryPlane };
    // Fixed outer Cartesian section, with a sliding interface endpoint.
    // The represented scalar must have identically zero normal derivative
    // on its complete active boundary trace (C2 even reflection). This does
    // not assert higher odd-derivative parity, a contact law, or periodicity.
    struct SymmetryPlane
    {
        int normalAxis = 0; // 0=x/r, 1=y/z
        double coordinate = 0;
        bool operator==(const SymmetryPlane& other) const
        { return normalAxis==other.normalAxis && coordinate==other.coordinate; }
    };
    struct ArcSample;
    struct LocalProjection;
    struct NativeProjectionQuery
    {
        int cellIndex = -1; // Actual unchanged leaf identity of the old anchor.
        Quantum::Geometry::CellChartPoint cellChartPoint;
        Vector2D displacement; // Added to the retained chart, not a rounded global inverse.
    };
    // Native topology stage of the same interface construction. This owns
    // ordered original edge-root identities and its immutable scalar source,
    // Explicit frameAt/quadrature queries below use the same scalar's raw
    // normal/meridional curvature and the original GL8 rule. They are not
    // full material panels or a reconstructed-curvature model. Bounded local
    // projection is explicit below. Each object is one complete closed or
    // declared/axis open chain. The complete-discovery entry additionally
    // permits several disjoint regular graph intervals in a native leaf.
    // Real split faces share vertex/edge identity; no uniform resampling or
    // geometric root snapping. Empty cells need exact phase certificates.
    class NativeGraphChain
    {
    public:
        const std::vector<LevelSetMappedInterface::Graph>& graphs() const { return graphs_; }
        // Compatibility view: throws if a cell has several arcs in this chain.
        const std::vector<int>& cellToGraph() const;
        const std::vector<std::vector<int>>& cellToGraphs() const { return cellGraphs_; }
        bool closed() const { return closed_; }
        bool axisymmetric() const { return axisymmetric_; }
        EndpointKind beginEndpointKind() const { return beginKind_; }
        EndpointKind endEndpointKind() const { return endKind_; }
        // Index into the negative regions of the complete native discovery.
        // Several surface components (e.g. a shell and its cavity) may share it.
        int negativeRegion() const { return negativeRegion_; }
        bool usesSource(const LevelSetHermiteGeometry::Snapshot& source) const
        { return source_.get()==&source; }
        // Fraction of the ORIGINAL chosen-coordinate extent, not arclength.
        // Endpoints retain the selected side/root; terminal scalar evaluation
        // uses a nearest-end chart approximation, never a global-point inverse.
        // Distinct ends and interior samples must remain distinguishable in
        // that chart. Unresolved samples reject; mapped volume root identities
        // are not merged or replaced to make a floating frame possible.
        LevelSetHermiteGeometry::Frame frameAt(int graphIndex,double fraction) const;
        // Conservative partial certificate on the connecting subarc. Bound
        // its retained monotone chart coordinates by a dyadic slab, then use
        // the original scalar's outward zero-set predicate in that slab.
        // False remains inconclusive; no sampled normals replace the proof.
        // The containing graph may extend beyond this interval to an axis pole.
        bool zeroSetPartialResolved(int graphIndex,int axis,double firstFraction,double lastFraction) const;
        std::vector<ArcSample> quadrature(int graphIndex) const;
        // Positive oriented subarc in original graph fractions. GL8 samples
        // retain the exact affine interval until terminal chart evaluation.
        std::vector<ArcSample> quadrature(int graphIndex,double firstFraction,double lastFraction) const;
        double physicalCoordinateSpan(int graphIndex) const;
        double graphFraction(int graphIndex,const Quantum::Geometry::CellChartPoint& point) const;
        // Find the retained interval inside an actual leaf, including leaves
        // partitioned into several native graphs. Uses exact chart bounds and
        // endpoint identities, never a rounded physical-coordinate inverse.
        // Adjacent shared endpoints use the first graph in chain order; all
        // other multiple incidences reject. Consumers still validate the full
        // source Frame/jet, rather than treating this lookup as a projection.
        int graphAt(int nativeCellIndex,const Quantum::Geometry::CellChartPoint& point) const;
        const std::shared_ptr<const LevelSetHermiteGeometry::Snapshot>& source() const { return source_; }
        const std::vector<Vector2D>& boundaryContacts() const { return boundaryContacts_; }
        const std::vector<SymmetryPlane>& symmetryPlanes() const { return symmetryPlanes_; }
        const SymmetryPlane* endpointSymmetryPlane(bool atEnd) const
        {
            const int id=endpointSymmetryPlanes_[atEnd?1:0];
            return id<0?nullptr:&symmetryPlanes_.at(id);
        }
        // Same bounded stationary/endpoint search as projectLocal below.
        // Retains actual leaf/root identity and endpoint-relative coordinates.
        // Uses global minimum active-leaf spacing for the original h/2 motion,
        // h/4 projection and 1e-9*h residual gates, consistent with advanceNative.
        // The caller must establish unchanged old/new topology and anchor provenance.
        // This is not a global nearest-point certificate or a load-band query.
        LocalProjection projectLocal(const NativeProjectionQuery& query) const;
        LocalProjection projectLocalEndpoint(const NativeProjectionQuery& query,
            int graphIndex,bool atGraphBegin) const;
        // Physical load/electric query on the retained curve, not a material
        // displacement. Reuses the original retained-graph stationary search
        // and open-end KKT conditions without imposing the h/4 remap radius.
        // The query must be in the active domain; no solid extension is made.
        // Successful identical physical queries can reuse a bounded result
        // inventory owned by this immutable chain/source, shared by copies.
        // Full coordinate bits are keys; current FP checks remain mandatory.
        LocalProjection projectedFrame(Vector2D query) const;
        // Finite retained-graph stationary inventory used by load projection,
        // including admissible original open endpoints. No nearest candidate
        // is discarded here: consumers can detect competing branch distances.
        // Duplicates at shared graph endpoints may remain. Empty is permitted;
        // this is a sampled search, not a globally complete root certificate.
        // Complete inventories can be reused for identical coordinate bits on
        // this immutable source/chain. Copies returned to callers are independent;
        // branch/force acceptance is still the consumer's responsibility.
        std::vector<LocalProjection> loadProjectionCandidates(Vector2D query) const;
        Vector2D minimumActiveSpacing() const { return minimumSpacing_; }
    private:
        friend class LevelSetInterfacePanels;
        NativeGraphChain() = default;
        LevelSetHermiteGeometry::Frame frameAtInterval(int graphIndex,double begin,double end,double fraction) const;
        LocalProjection projectedFrameUncached(Vector2D query) const;
        std::vector<LocalProjection> loadProjectionCandidatesUncached(Vector2D query) const;
        std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source_;
        std::vector<LevelSetMappedInterface::Graph> graphs_;
        std::vector<int> mapping_;
        std::vector<std::vector<int>> cellGraphs_;
        // IDs from the authoritative actual-subedge discovery, never display positions.
        std::vector<std::array<int,2>> rootNodes_;
        // Filled only by completed material-panel construction, from its
        // already evaluated endpoints and original GL8 frames. Owned by this
        // immutable chain/source; no global or cross-time mutable cache.
        // Slots 0/1 are endpoints, slots 2..9 the original GL8 fractions.
        std::vector<std::array<LevelSetHermiteGeometry::Frame,10>> fixedProjectionFrames_;
        // Exact physical coordinates of the same fixed samples and graph
        // boxes. Prepared once after panel construction, immutable and shared
        // by chain copies. Rebuilt when native cells are rebound by AMR.
        struct ProjectionCoordinates;
        std::shared_ptr<const ProjectionCoordinates> projectionCoordinates_;
        void prepareProjectionCoordinates();
        // Geometry-only endpoint-root data prepared on this source. Immutable,
        // owning through the scalar snapshot; copying/moving a chain cannot
        // leave graph/vector borrows behind. Arbitrary fractions still use
        // the original transverse-root and projection algorithms.
        struct SamplingData;
        std::shared_ptr<const SamplingData> sampling_;
        Vector2D minimumSpacing_;
        std::vector<Vector2D> boundaryContacts_;
        std::vector<SymmetryPlane> symmetryPlanes_;
        std::array<int,2> endpointSymmetryPlanes_{{-1,-1}};
        bool closed_=true,axisymmetric_=false;
        EndpointKind beginKind_=EndpointKind::None,endKind_=EndpointKind::None;
        int negativeRegion_ = -1;
    };
    static NativeGraphChain discoverNativeGraphs(const QuadtreeMesh& mesh,
        std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source,
        bool axisymmetric=false,const std::vector<Vector2D>& boundaryContacts={});
    // Complete native zero set. Components share one immutable scalar and
    // exact shared-edge root inventory. No component is selected or dropped.
    // Multiple arcs per chart require a certified monotone transverse partial.
    // Quartic complete discovery can partition a native cell into bounded
    // dyadic geometry charts while retaining its original scalar and DOFs;
    // unresolved cells throw InterfaceGeometryResolution for refinement/retry.
    // Native Q4 isolated contacts retain one root with zero/one/two inward
    // graph incidences, certified by original boundary phases and nonzero
    // normal derivative on its enclosure. Negative fragments join across
    // positive-length intervals, not through an isolated touching point.
    // This does not admit singular zeros or positive-length zero traces.
    // Axis-to-axis open meridians represent detached closed RZ surfaces.
    // Explicit symmetryPlanes additionally admit sliding endpoints on outer
    // sections. They remain actual domain boundaries in the volume ledger;
    // no mirror copies or artificial end disks enter the interface inventory.
    struct NativeNegativeRegion
    {
        // Actual cells intersecting this region; a cut cell may occur in
        // several regions. Pure cells provide unambiguous interior witnesses.
        std::vector<int> cells, pureCells;
        bool touchesNonAxisBoundary = false;
        long double fixedBoundaryVolume = 0; // Green contribution of actual negative boundary intervals.
    };
    struct NativeComponents
    {
        std::vector<NativeGraphChain> curves;
        std::vector<NativeNegativeRegion> negativeRegions;
        // +/-1 certifies the interior phase of an active cell with no interior
        // interface arc (an explicit side-attached arc may bound that cell);
        // 0 denotes a cut or inactive cell. Includes the monotone-boundary
        // certificate used where a Bernstein control-sign test is insufficient.
        std::vector<int> cellPhase;
    };
    static NativeComponents discoverNativeComponents(const QuadtreeMesh& mesh,
        std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source,
        bool axisymmetric=false,const std::vector<Vector2D>& boundaryContacts={},
        const std::vector<SymmetryPlane>& symmetryPlanes={});
    // Explicit native Q4 zero-edge topology. Adds canonical side-attached
    // support with actual incident phase/vertex/subedge identity; the default
    // discovery and legacy sensitivity still reject positive-length zeros.
    static NativeComponents discoverNativeTraceComponents(const QuadtreeMesh& mesh,
        std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source,
        bool axisymmetric=false,const std::vector<Vector2D>& boundaryContacts={},
        const std::vector<SymmetryPlane>& symmetryPlanes={});
    struct QuadratureNode
    {
        LevelSetHermiteGeometry::Frame frame;
        double weight = 0.0; // Positive arc-length quadrature weight.
        double localArclength = 0.0; // Cumulative GL8 integral from panel begin; composite for native panels.
        double physicalWeight = 0.0; // weight in planar mode, 2*pi*r*weight in RZ.
        double localPhysicalMeasure = 0.0; // Corresponding cumulative physical integral.
        double nativeGraphFraction = 0.0; // Native material panels only.
    };
    struct Panel
    {
        int cellIndex = -1; // Solver index, not label or Cartesian index.
        Vector2D begin, end; // Consecutive endpoints, oriented t=(-ny,nx).
        bool coordinateIsX = true; // Monotone coordinate used by quadrature.
        LevelSetHermiteGeometry::Frame beginFrame, endFrame;
        // Chord-midpoint projection by default. Declared same-edge panels and
        // all explicit sector declared-chain panels use the actual graph-
        // coordinate midpoint, retaining its true stationarity residual.
        LevelSetHermiteGeometry::Frame frame;
        double measure = 0.0;
        double physicalMeasure = 0.0; // Sum of physical quadrature weights.
        std::vector<QuadratureNode> quadratureNodes; // Increasing local arclength.
        // Native positive composite GL8 diagnostics. Relative error ESTIMATES
        // compare each accepted rule with its halves; not rigorous bounds.
        // Zero intervals means the legacy, uncontrolled non-native path.
        std::size_t measureQuadratureIntervals = 0;
        std::array<double,2> measureQuadratureRelativeError{{0,0}}; // ds, physical area.
    };
    struct PanelGraph
    {
        struct Endpoint
        {
            // Exact topology identity in the original complete grid-edge
            // root inventory. Sides W/S/E/N = 0/1/2/3; roots ordered in the
            // positive coordinate direction, including exact corner roots.
            // Not a rounded coordinate or an ordinal among graph-only roots.
            int side = -1, rootOrdinal = -1;
        };
        int cellI = -1, cellJ = -1;
        Vector2D localBegin, localEnd; // Original cell-local edge-root parameters.
        Vector2D cellSize;
        bool coordinateIsX = true;
        Endpoint beginRoot, endRoot;
    };
    struct ArcSample
    {
        LevelSetHermiteGeometry::Frame frame; // Always carries cell-local provenance.
        double localGraphCoordinate = 0.0;
        double weight = 0.0; // Positive ds quadrature weight.
        double physicalWeight = 0.0; // ds or 2*pi*r*ds, according to axisymmetric().
        double gradientMagnitude = 0.0; // SAME raw cell jet, never assumed to be one.
        // Native quadrature carries this chart instead of a potentially
        // collapsed localGraphCoordinate (which is NaN for that path).
        Quantum::Geometry::CellChartCoordinate localGraphChart;
        bool hasLocalGraphChart=false;
    };
    struct LocalProjectionQuery
    {
        int cellI = -1, cellJ = -1;
        Vector2D cellLocalPoint; // Closed-cell anchor on the unchanged source lattice.
        Vector2D displacement; // Physical displacement, NOT a new rounded global point.
    };
    struct LocalProjection
    {
        LevelSetHermiteGeometry::Frame frame; // Actual source jet in its retained local/leaf chart.
        int panelIndex = -1;
        bool endpoint = false, atPanelBegin = false;
        PanelGraph::Endpoint root; // Original root identity; valid only at an endpoint.
        double distance = 0.0;
        double inwardDerivative = 0.0; // Endpoint: (frame-query).t_in; interior: zero.
        double kktResidual = 0.0; // Endpoint: max(0,-inwardDerivative); interior: |d.t|.
        // Native entry only: original graph fraction that generated frame.
        // Retained for subarc queries, never recovered from its display point.
        double nativeGraphFraction = 0.0;
    };

    LevelSetInterfacePanels(const QuadtreeMesh& mesh,
        const LevelSetHermiteGeometry::Snapshot& snapshot, bool axisymmetric = false,
        std::vector<Vector2D> boundaryContacts = {},
        CartesianGridScope scope = CartesianGridScope::RootOnly);
    // Explicit retained-query mode. The referent must remain immutable (also
    // through any other alias); ownership prevents a dangling borrowed source.
    // Both overloads use the same Panel construction path; source retention
    // additionally enables cell-local graph and subarc queries.
    LevelSetInterfacePanels(const QuadtreeMesh& mesh,
        std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> snapshot, bool axisymmetric = false,
        std::vector<Vector2D> boundaryContacts = {},
        CartesianGridScope scope = CartesianGridScope::RootOnly);

    // Explicit native material-panel adapter on the SAME retained graph chain.
    // No Cartesian i/j or uniform replacement grid is manufactured. Frames and
    // GL8 measures come from that chain; legacy constructors remain unchanged.
    explicit LevelSetInterfacePanels(NativeGraphChain chain);
    // Zero-time mesh rebinding, not material motion: target is a freshly
    // discovered complete Q4 chain. Every retained leaf must have identical
    // bounds/polynomial, orientation, endpoint roots and ordered neighbors.
    // Reuse owned GL8 frames/measures only if the numerical representation is
    // also identical; otherwise use the original full panel constructor.
    // Target owns its NEW scalar/sampling data. No old projection cache moves.
    LevelSetInterfacePanels rebindNative(NativeGraphChain target,
        const std::vector<int>& sourceToTargetCell) const;
    bool nativeMeasuresReused() const noexcept { return nativeMeasuresReused_; } // Observational only.
    bool usesNativeLeaves() const noexcept { return nativeChain_.has_value(); }
    const NativeGraphChain& nativeChain() const;

    const std::vector<Panel>& panels() const noexcept { return panels_; }
    const std::vector<int>& cellToPanel() const;
    bool axisymmetric() const noexcept { return axisymmetric_; }
    bool closed() const noexcept { return closed_; }
    EndpointKind beginEndpointKind() const noexcept { return beginEndpointKind_; }
    EndpointKind endEndpointKind() const noexcept { return endEndpointKind_; }
    const std::vector<Vector2D>& boundaryContacts() const noexcept { return boundaryContacts_; }
    bool hasRetainedSource() const noexcept { return static_cast<bool>(source_); }
    bool usesSource(const LevelSetHermiteGeometry::Snapshot& snapshot) const noexcept
    { return source_.get() == &snapshot; }
    const PanelGraph& panelGraph(int panelIndex) const;
    // Endpoints are the retained graph's chosen CELL-LOCAL coordinate, NOT
    // global coordinates or arclength. Either order is accepted; samples follow
    // that order, weights stay positive, and the negative-to-positive normal
    // is unchanged. Equal endpoints return no samples. A nonzero interval must
    // resolve every GL8 interior coordinate in its local chart, or it rejects.
    // Split GL8 sums are approximations, not an exact additive primitive.
    std::vector<ArcSample> subarcQuadrature(int panelIndex,
        double localCoordinateBegin, double localCoordinateEnd) const;

    // Bounded local projection on the retained original graphs, using their
    // actual transverse roots and local jets. Requires shared-source mode.
    // The caller verifies the old anchor's provenance and unchanged lattice;
    // this object cannot verify which old scalar supplied the query.
    // Physical differences retain the Snapshot's actual cell origins, nominal
    // h*local offsets and exact local==1 boundary convention in expansion
    // arithmetic. Only terminal distances/dot products are rounded; no global
    // display coordinate is inverted and no endpoint is snapped or merged.
    // With h=min(cellSize), displacement <= h/2, projection <= h/4 and the
    // original 1e-9*h level-set/stationarity (endpoint KKT) gates are required.
    // This is a finite, sampled, local stationary search with the original
    // nonadjacent whole-cell-box guard, NOT a unique global nearest certificate.
    LocalProjection projectLocal(const LocalProjectionQuery& query) const;
    // Evaluate the specified ORIGINAL endpoint without searching/replacing it.
    // The caller owns its material pin/axis identity. A one-sided endpoint KKT
    // residual is different from interior orthogonality; frame still reports
    // the actual absolute tangential residual in both cases.
    LocalProjection projectLocalEndpoint(const LocalProjectionQuery& query,
        int panelIndex, bool atPanelBegin) const;

    // Load/trace projection of a represented PHYSICAL query, independent of
    // the material-motion bounds above. The caller owns its load-band extent
    // (e.g. 3h); no h/2 motion, h/4 projection or nonadjacent material-box
    // policy is changed or reused here. Requires the retained shared source.
    // Reuses the Snapshot's original eight-seed stationary search, but accepts
    // only candidates in an actual retained curve-panel incidence. A pure
    // wet-sector corner is not a free-interface stationary candidate.
    // Open-chain ends are evaluated at their ORIGINAL local edge roots, with
    // one-sided KKT and the original 2048*eps*coordinate-scale boundary gate.
    // Here endpoint=true denotes a constrained OPEN-CHAIN end (and retains
    // its root identity); a stationary internal panel crossing is not a wall
    // endpoint. No root snapping/global endpoint inverse or closest-point
    // certificate is introduced. If that inventory has no admissible result,
    // complete the search on retained panel graphs using their original GL8
    // section seeds and the shared 128-iteration stationary bracket. Section
    // roots retain endpoint-relative coordinates, including sub-ULP upper-side
    // offsets; scalar/stationarity gates stay unchanged. This remains a finite
    // sampled stationary search, not a certified global nearest-point solver.
    // Native panels call the same retained-graph search directly (there is no
    // uniform Newton seed lattice); they keep real leaf/chart provenance and
    // the same boundary gate. Old material projection bounds are unchanged.
    LocalProjection projectedFrame(Vector2D query) const;

private:
    static NativeComponents discoverNativeComponentsImpl(const QuadtreeMesh& mesh,
        std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source, bool axisymmetric,
        const std::vector<Vector2D>& boundaryContacts, bool legacy,
        const std::vector<SymmetryPlane>& symmetryPlanes={},bool nativeTrace=false);
    struct NativeMetadata {};
    LevelSetInterfacePanels(NativeGraphChain chain, NativeMetadata);
    bool nativeMeasuresReused_ = false;
    std::vector<Panel> panels_;
    std::vector<int> cellToPanel_;
    bool axisymmetric_ = false;
    bool closed_ = true;
    EndpointKind beginEndpointKind_ = EndpointKind::None;
    EndpointKind endEndpointKind_ = EndpointKind::None;
    std::vector<Vector2D> boundaryContacts_;
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source_;
    std::vector<PanelGraph> graphs_;
    std::optional<NativeGraphChain> nativeChain_;
    void initialize(const QuadtreeMesh& mesh, const LevelSetHermiteGeometry::Snapshot& snapshot,
        CartesianGridScope scope);
};
