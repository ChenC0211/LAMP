#pragma once

#include <array>
#include <cstddef>
#include <vector>
#include <stdexcept>
#include <string>
#include <utility>
#include <type_traits>
#include <memory>
#include "Core/Base/Vector2D.h"
#include "Core/Geometry/ActiveHermiteSensitivity.h"
#include "Core/Geometry/LevelSetMappedSegment.h"

class LevelSetMappedSegment;
namespace LevelSetHermiteGeometry { class Snapshot; }

// Same-source free-curve integration in the shared-canonical-bound embedding.
// Does not change the legacy nominal-h point, panel or force geometry.
class LevelSetMappedInterface
{
public:
    // Prepared by explicit same-source trace discovery only. These endpoints
    // are actual subedge vertices, not finite roots of a zero polynomial.
    class NativeSideSupport
    {
    public:
        struct Side {int nativeCellIndex=-1,phase=0;bool available()const{return nativeCellIndex>=0;}};
        bool usesSource(const LevelSetHermiteGeometry::Snapshot& s)const{return source_.get()==&s;}
        bool usesNativeScalar(const void* s)const{return scalarIdentity_==s;}
        int ownerCell()const{return owner_;}
        int side()const{return side_;} // Actual owner W/S/E/N = 0/1/2/3.
        int normalAxis()const{return side_==0 || side_==2?0:1;}
        bool horizontal()const{return horizontal_;}
        int edgeLabel()const{return edgeLabel_;}
        Vector2D begin()const{return begin_;}
        Vector2D end()const{return end_;}
        int beginVertex()const{return beginVertex_;}
        int endVertex()const{return endVertex_;}
        const Side& minus()const{return minus_;}
        const Side& plus()const{return plus_;}
    private:
        friend class LevelSetInterfacePanels;
        NativeSideSupport(std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source,const void* scalar,
            int owner,int side,bool horizontal,int edge,Vector2D begin,Vector2D end,int beginVertex,int endVertex,Side minus,Side plus)
            :source_(std::move(source)),scalarIdentity_(scalar),owner_(owner),side_(side),horizontal_(horizontal),edgeLabel_(edge),
            begin_(begin),end_(end),beginVertex_(beginVertex),endVertex_(endVertex),minus_(minus),plus_(plus){}
        std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source_;
        const void* scalarIdentity_;
        int owner_,side_;bool horizontal_;int edgeLabel_;
        Vector2D begin_,end_;int beginVertex_,endVertex_;Side minus_,plus_;
    };
    struct RootResolutionRequest
    {
        std::size_t segmentIndex=0; // Index in the integrate() faceCuts argument.
        std::size_t rootIndex=0; // Index in that segment's public roots() order.
    };
    class RootResolutionError : public std::runtime_error
    {
    public:
        std::vector<RootResolutionRequest> requests;
        RootResolutionError(std::string message,std::vector<RootResolutionRequest> requested)
            : std::runtime_error(std::move(message)),requests(std::move(requested)) {}
    };
    struct Endpoint
    {
        int side = -1; // W,S,E,N = 0,1,2,3.
        int rootOrdinal = -1; // Complete edge roots in increasing local coordinate, including endpoints.
    };
    struct Graph
    {
        int cellI = -1, cellJ = -1;
        bool coordinateIsX = true;
        Endpoint begin, end;
        // Native graph: i/j stay -1 and this is an actual leaf index.
        // Uniform graph: this stays -1. Mixed addressing is rejected.
        int nativeCellIndex = -1;
        // A graph may cover a dyadic subbox of the ORIGINAL native cell.
        // These are geometry charts only: cell identity, scalar coefficients,
        // physical Jacobian and sensitivity columns remain those of that cell.
        // Endpoint side/ordinal refers to the complete subbox-edge inventory.
        Vector2D localLower{0.0,0.0}, localUpper{1.0,1.0};
        // Complete Q4 components use error-controlled moment/Jacobian
        // quadrature. Legacy single-chain graphs retain their GL8 rule.
        bool controlledIntegration = false;
        // Present only on a proved regular native zero subedge. Ordinary
        // begin/end root ordinals remain unset for this distinct support.
        std::shared_ptr<const NativeSideSupport> nativeSide;
    };
    struct Region
    {
        int regionId = -1;
        std::array<Vector2D,4> vertices; // Actual convex quadrilateral, cyclic in either direction.
    };
    template<std::size_t Count> struct CellArc
    {
        static_assert(Count==16 || Count==25,"unsupported mapped cell sensitivity layout");
        int regionId = -1;
        // Chosen dimensionless CELL coordinate q=offset+s or 1-offset-s,
        // 0<=s<=width. These independent terminal approximations are not two
        // rounded endpoints and are not interval bounds. All integrals use
        // positive arclength and the phi<0 -> phi>0 normal, not graph order.
        bool fromEnd = false;
        long double offset = 0.0L, width = 0.0L;
        long double axialMoment = 0.0L; // integral y*n_y ds (XY), or 2*pi*r*z*n_z ds (RZ).
        long double radialMoment = 0.0L; // integral x*n_x ds, or pi*r*r*n_r ds.
        long double momentErrorEstimate = 0.0L; // GL8 versus subdivided GL8; not an interval bound.
        // PHYSICAL data columns. Count=16: original SW SE NE NW, each
        // value,dx,dy,dxy. Count=25: tensor ordering documented below.
        // Uniform sources use their original nominal hx,hy; native sources
        // use the actual leaf widths. The shared-bound physical Jacobian
        // scales integration and gradient consistently in both cases.
        std::array<long double,Count> sensitivity{};
        long double minTransverseDerivative = 0.0L;
        long double maxScalarResidual = 0.0L;
        std::size_t quadraturePoints = 0;
        // sensitivity[] remains the formal cell-data diagnostic. Cubic native
        // integration also applies the retained scalar's hanging constraints
        // to produce independentSensitivity. Use THAT list for a global
        // correction, never the formal columns. Pin/axis constraints and the
        // choice of free derivatives still belong to the caller.
        int nativeCellIndex = -1;
        bool requiresNativeDofConstraints = false;
        using IndependentSensitivity=std::conditional_t<Count==16,
            LAMP::Geometry::ActiveHermiteNodalSensitivity,LAMP::Geometry::ActiveQuarticHermiteSensitivity>;
        std::vector<IndependentSensitivity> independentSensitivity;
    };
    using Arc = CellArc<16>;
    // Different result type: these 25 columns must never be passed to an
    // old 16-corner correction. Tensor order is 5*y+x; each direction uses
    // value(0), value(1), derivative(0), derivative(1), midpoint value.
    // Derivative columns are with respect to PHYSICAL cell data, with the
    // original nominal hx/hy scaling. They are formal independent cell
    // variations, not a pullback through the Q4 nodal-data reconstruction.
    // Equal-level input leaves independentSensitivity empty. The native data
    // adapter applies its shared constraints and returns typed vertex/edge/
    // center identities; they are not old four-component nodal columns.
    // The explicit native quartic Snapshot exposes that adapter; these data
    // columns still do not claim a global correction or a Q4 nodal pullback.
    using QuarticArc = CellArc<25>;
    template<class ArcType> struct IntegrationResult
    {
        std::vector<ArcType> arcs;
        std::size_t sourceRootIncidences = 0, retainedCuts = 0;
    };
    using Result = IntegrationResult<Arc>;
    using QuarticResult = IntegrationResult<QuarticArc>;

    // Synchronous query: source and face cuts must remain alive/immutable
    // for this call; Result retains no references. Every cut must use exactly
    // this Snapshot. Caller supplies the complete regular panel graph (one
    // source cell, either logical uniform or explicit native-leaf ID) and all
    // external region-boundary cuts. A native side may be covered by multiple
    // actual collinear half faces with exactly gap-free endpoint coverage.
    // Their retained root inventories remain unchanged. Convex pieces sharing a regionId form one
    // region; their internal seams need not cut the arc. Different owners may
    // not overlap. A coincident, unresolved, ambiguous, uncovered, or unordered
    // graph rejects; no root/short arc is merged using a tolerance.
    // When a retained face-root enclosure prevents a quadrature-side test,
    // RootResolutionError identifies ONLY the input roots needing further
    // resolution. The caller must refine shared face partitions and restart
    // all dependent measures; this query never edits an input face root.
    // Legacy graphs use fixed GL8. Complete native Q4 graphs request bounded
    // GL8 subdivision controlling both Green moments and all sensitivity
    // columns. These are floating quadratures and error estimates, NOT
    // certified enclosures or nonlinear volume-conservation acceptance.
    static Result integrate(const LevelSetHermiteGeometry::Snapshot& source,
        const Graph& graph, const std::vector<const LevelSetMappedSegment*>& faceCuts,
        const std::vector<Region>& regions, bool axisymmetric = true);
    // Same roots, regions and graph-selected quadrature, with the selected
    // quartic scalar and all 25 formal cell sensitivities. The legacy entry
    // rejects quartic sources rather than returning empty/wrong columns.
    static QuarticResult integrateQuartic(const LevelSetHermiteGeometry::Snapshot& source,
        const Graph& graph, const std::vector<const LevelSetMappedSegment*>& faceCuts,
        const std::vector<Region>& regions, bool axisymmetric = true);
    // Same integration with a constant physical Green-potential origin.
    // Axial: (y-y0), or 2*pi*r*(z-z0). Radial: (x-x0), or
    // pi*(r-r0)*(r+r0). These retain the original divergence/volume law.
    // Fixed-boundary moments must use THIS SAME origin. Sensitivities remain
    // volume variations, not derivatives of an isolated free-boundary moment.
    // The finite origin may lie outside the domain; it is not a query point.
    // Zero origin preserves the original arithmetic path. Differences use
    // retained rational locations before rounding the translated coordinate.
    static QuarticResult integrateQuartic(const LevelSetHermiteGeometry::Snapshot& source,
        const Graph& graph, const std::vector<const LevelSetMappedSegment*>& faceCuts,
        const std::vector<Region>& regions, bool axisymmetric, Vector2D greenOrigin);
    // Optional caller-owned cumulative worker attribution. Each actual
    // batch/group owns its profile; no global counters or atomics are used.
    // Time fields are exclusive: isolation/refinement exclude restrictions.
    // restrictionConstruction times only the actual layout.polynomial backend;
    // coordinate and endpoint preparation belong to remainingIntegration.
    // remainingIntegration includes graph/cut preparation, ownership, jets,
    // quadrature, recovery and other work not assigned to a named segment.
    struct NativeMomentWorkProfile
    {
        double restrictionConstruction=0,initialRootCacheLookup=0,initialRootCacheRemember=0;
        double initialRootIsolation=0,rootMappingRefinement=0,fixedAndSetup=0,remainingIntegration=0,total=0;
        std::size_t restrictionQueries=0,initialRootCacheQueries=0,initialRootCacheHits=0,initialRootCacheMisses=0;
        std::size_t initialRootCacheRememberQueries=0,initialRootIsolations=0,rootMappings=0,rootRefinements=0;
    };
    struct GraphMoments { long double radial=0, axial=0, errorEstimate=0; };
    // Geometry-only Green moments of a complete native cubic/Q4 graph. Uses the
    // original scalar and controlled rule, integrating the graph-aligned
    // moment and recovering the other by parts with retained endpoint roots.
    // Controls moments only; no correction Jacobian is computed or exposed.
    // Error estimates remain floating quadrature diagnostics, not bounds.
    static GraphMoments nativeMoments(const LevelSetHermiteGeometry::Snapshot& source,
        const Graph& graph, bool axisymmetric);

    struct RegionMoments
    {
        int regionId=-1;
        long double capacity=0;
        long double axialLiquid=0,radialLiquid=0,axialGas=0,radialGas=0;
        long double fixedAxialLiquid=0,fixedRadialLiquid=0;
        long double freeAxialLiquid=0,freeRadialLiquid=0;
        long double errorEstimate=0;
        std::size_t freeArcs=0;
    };
    struct NativeRegionMoments
    {
        std::vector<RegionMoments> regions;
        // Geometric clip ledger, not an antisymmetric material flux ledger.
        long double exteriorFixedAxialLiquid=0,exteriorFixedRadialLiquid=0;
        long double coincidentAxialLiquid=0,coincidentRadialLiquid=0;
        std::size_t coincidentArcs=0,sourceRootIncidences=0,localRootRefinements=0;
    };
    // Each target is an independent complete rectangular partition and keeps
    // its own arithmetic Green gauge. Targets may overlap one another; their
    // regions must individually cover the same native cell exactly once.
    struct NativeMomentTarget
    {
        std::vector<Region> regions;
        Vector2D greenOrigin{0,0};
    };
    // Only fixed/capacity fields exist for unselected regions. No missing
    // phase/free integral is represented as a zero volume.
    struct FixedRegionMoments
    {
        int regionId=-1;
        long double capacity=0,fixedAxialLiquid=0,fixedRadialLiquid=0;
    };
    struct NativeSelectedMomentTarget
    {
        std::vector<Region> regions; // Complete original native-cell partition.
        Vector2D greenOrigin{0,0};
        std::vector<int> selectedRegionIds;
    };
    struct NativeSelectedRegionMoments
    {
        std::vector<RegionMoments> selectedRegions;
        std::vector<FixedRegionMoments> fixedRegions;
        long double exteriorFixedAxialLiquid=0,exteriorFixedRadialLiquid=0;
        long double coincidentAxialLiquid=0,coincidentRadialLiquid=0;
        std::size_t coincidentArcs=0,sourceRootIncidences=0,localRootRefinements=0;
    };
    // Same complete fixed/coincident partition and original per-target gauge.
    // Free moments and controlled quadrature run only on selected owners.
    static std::vector<NativeSelectedRegionMoments> integrateNativeQuarticSelectedRegionMomentsBatch(
        const LevelSetHermiteGeometry::Snapshot& source,int nativeCellIndex,
        const std::vector<Graph>& graphs,
        const std::vector<const LevelSetMappedSegment::NativeQuarticTrace*>& faceCuts,
        const std::vector<NativeSelectedMomentTarget>& targets,bool axisymmetric,
        NativeMomentWorkProfile* work = nullptr);
    // Shared native graph cuts/root/jet traversal, with direct moment channels
    // for each target. faceCuts are common immutable input inventories; every
    // RootResolutionError uses this one input's segment/event indices.
    // Optional global moments consume each complete input graph at zero
    // origin in the SAME traversal, independent of target clip ownership.
    // Output is in input graph order and assigned only after full success.
    static std::vector<NativeRegionMoments> integrateNativeQuarticRegionMomentsBatch(
        const LevelSetHermiteGeometry::Snapshot& source,int nativeCellIndex,
        const std::vector<Graph>& graphs,
        const std::vector<const LevelSetMappedSegment::NativeQuarticTrace*>& faceCuts,
        const std::vector<NativeMomentTarget>& targets,bool axisymmetric,
        std::vector<GraphMoments>* nativeGraphMoments = nullptr,
        NativeMomentWorkProfile* work = nullptr);
    // Geometry only: actual axis-aligned rectangular pieces inside ONE active
    // native Q4 cell, complete regular graphs and all actual typed boundary
    // cuts from this identical Snapshot. Each fixed side takes its inside
    // phase trace. Certified free/fixed coincidence is counted only once.
    // No Jacobian, phase repair or volume inventory is produced. Legacy
    // sensitivity and native-edge publication entry points are unchanged.
    // RootResolutionError indices refer to faceCuts/events(), not legacy roots.
    static NativeRegionMoments integrateNativeQuarticRegionMoments(
        const LevelSetHermiteGeometry::Snapshot& source,int nativeCellIndex,
        const std::vector<Graph>& graphs,
        const std::vector<const LevelSetMappedSegment::NativeQuarticTrace*>& faceCuts,
        const std::vector<Region>& regions,bool axisymmetric,Vector2D greenOrigin);
private:
    template<bool Selected>
    static std::vector<std::conditional_t<Selected,NativeSelectedRegionMoments,NativeRegionMoments>>
        nativeQuarticRegionMomentsBatchImpl(const LevelSetHermiteGeometry::Snapshot&,int,
        const std::vector<Graph>&,
        const std::vector<const LevelSetMappedSegment::NativeQuarticTrace*>&,
        const std::vector<NativeMomentTarget>&,bool,std::vector<GraphMoments>*,
        const std::vector<std::vector<int>>*,NativeMomentWorkProfile*);

};
