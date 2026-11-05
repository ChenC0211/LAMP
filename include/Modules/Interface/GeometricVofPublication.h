#pragma once

#include "Core/Geometry/PlanarCellCut.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Base/TimeLayer.h"
#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace LAMP::Geometry::Detail::FiniteVolume { class Dyadic; }
namespace LAMP::Mesh { class NativeLeafOverlap; }

namespace LAMP::Interface
{
class GeometricVofAdvectorSlab;
class GeometricVofTransport;
// Clock value only: no Hermite state/identity is used by this VOF method.
using GeometricVofTimeLayer=LAMP::TimeLayer;
class GeometricVofIdentity
{
public:
    GeometricVofIdentity()=default;
    std::uint64_t value() const noexcept { return value_; }
    friend bool operator==(GeometricVofIdentity a,GeometricVofIdentity b) noexcept { return a.value_==b.value_; }
    friend bool operator!=(GeometricVofIdentity a,GeometricVofIdentity b) noexcept { return !(a==b); }
private:
    explicit GeometricVofIdentity(std::uint64_t value):value_(value){}
    std::uint64_t value_=0;
    friend class GeometricVofPublication;
};
enum class GeometricVofGeometryKind { CompatibleFan, LocalPlic };
// One actual native rectangle and a cell-centred PLIC plane. Both inventories
// are authoritative; measured[] describes the actual plane, not those stocks.
// Unresolved cells retain both stocks and explicitly require the declared
// donor-fraction approximation instead of claiming resolved phase shapes.
struct GeometricVofPlicCell
{
    Geometry::PlanarCellCut::Rectangle rectangle;
    Vector2D normal{1,0};
    long double offset=0;
    std::array<long double,2> inventory{},measured{}; // Liquid, gas.
    bool resolved=false;
};
// Explicit finite-resolution model; it is not an interface-isolation proof.
// Values 1/2 are retired experimental laws, never aliases of the new model.
enum class GeometricVofLocalResolution {VolumeOnly=0,NativeInterfaceSupport=3};
// Dense native results of one owner query: final curvature availability and
// the frozen base/opposed-orientation eligibility of actual geometric data.
struct GeometricVofLocalResolutionMasks
{
    std::vector<bool> capillarySubgrid;
    std::vector<bool> geometryUnsupported;
};
struct GeometricVofConfig
{
    long double relativeVolumeTolerance=1e-12L;
    unsigned maximumVolumeIterations=80,maximumBracketExpansions=32;
    GeometricVofLocalResolution localResolution=GeometricVofLocalResolution::VolumeOnly;
};
struct GeometricVofTriangle
{
    using Point=Geometry::PlanarCellCut::Point;
    using AnchoredPoint=Geometry::PlanarCellCut::AnchoredPoint;
    std::array<AnchoredPoint,3> anchoredVertices{}; // Authoritative actual sites in the ONE common frame.
    Geometry::PlanarCellCut::AnchoredTriangleCut anchoredCut;
    std::array<Point,3> vertices{}; // Same-owner lossy diagnostic projections only.
    std::array<long double,3> scalar{};
    std::array<std::size_t,3> vertexIds{};
    Geometry::PlanarCellCut::LinearTriangleCut cut; // Projection of anchoredCut; never recut or authoritative.
};
struct GeometricVofCell
{
    int nativeCell=-1;
    long double capacity=0,inventory=0,measured=0;
    long double centerScalar=0;
    std::vector<std::size_t> boundaryVertices; // Actual CCW, including hanging nodes.
    std::vector<GeometricVofTriangle> triangles; // Mixed cells only; pure cells homogeneous.
};
struct GeometricVofPhaseRange
{
    long double begin=0,end=1; // Rounded diagnostic parameters; actual geometry is vertices below.
    bool minusLiquid=false,plusLiquid=false;
    std::array<std::size_t,2> vertices{}; // Same-owner actual AP segment endpoints, native tangential order.
};
struct GeometricVofFace
{
    int nativeRow=-1,negativeCell=-1,positiveCell=-1;
    std::array<std::size_t,2> vertices{}; // In increasing native tangential coordinate.
    bool coincident=false;
    std::array<GeometricVofPhaseRange,2> ranges{};
    unsigned rangeCount=0;
    long double liquidMinus=0,liquidPlus=0; // Actual physical side traces, not transport.
};
struct GeometricVofFacet
{
    enum class Kind { CellFan, NativeCoincident };
    Kind kind=Kind::CellFan;
    std::array<std::size_t,2> vertices{}; // Canonical IDs; negative phase is on left.
    Geometry::PlanarCellCut::Point unitNormal{};
    int liquidCell=-1,gasCell=-1,nativeRow=-1;
    long double length=0,physicalMeasure=0;
};
struct GeometricVofTopology
{
    std::vector<std::size_t> facetComponent;
    std::vector<std::size_t> junctionVertices; // True degree !=2, never duplicated to hide it.
    std::vector<std::size_t> axisEndpoints,exteriorEndpoints;
    std::size_t components=0;
};
// Receipt of the most recent actual full LocalPlic build. Pin attachment copies
// that receipt; checkpoint restore has no ephemeral build receipt. These counts
// never replace physical reports/gates, and bytes include native support once.
struct GeometricVofLocalPlicReuseReport
{
    std::size_t geometryCells=0,vertexSets=0,faceSets=0,nodeRings=0;
    std::size_t mixedReused=0,mixedReconstructed=0,mixedDependencyChanged=0,pinnedReconstructed=0;
    std::size_t retainedBytes=0;
};
struct GeometricVofReport
{
    std::size_t mixedCells=0,pureLiquidCells=0,pureGasCells=0,zeroCollars=0,coincidentFaces=0;
    unsigned maximumVolumeIterations=0;
    long double maximumRelativeVolumeResidual=0,maximumRelativeCapacityResidual=0;
    long double inventory=0,measured=0;
    GeometricVofLocalPlicReuseReport localPlicReuse;
};
// ONE immutable alpha owner and its actual derived geometry. Preparation does
// not accept an App step, advance time, transport alpha, or mutate the mesh.
// Input alpha uses native cell indices. LocalPlic can select complete leaves;
// inactive native slots are explicit empty sentinels, never a third phase.
// Typed same-mesh LocalPlic descendants may share bounded private geometry.
// Every publication owns fresh stocks/clock and repeats all global physical
// gates. Only dependency-identical resolved planes with an identical complete
// quiet native FP/errno entry may retain the predecessor's actual geometry.
class GeometricVofPublication
{
    class MeshValidationBatch;
public:
    // Actual selected-domain reentrant native vertex, retained as physical
    // point and wall-leg directions. No inventory, alpha or clock is changed.
    // Four actual corner leaves (three fluid, one solid). Default law keeps
    // exact whole-neighbour wet-liquid/dry-gas purity. Explicit prescribed
    // wall data may coexist with unresolved majority-neighbour bulk stocks,
    // only for a truly resolved mixed contact endpoint. It certifies no
    // internal shared trace, complete continuous star or minor geometry.
    // Unresolved positive contact microphases and ambiguous topology reject.
    class PinCertificate
    {
    public:
        enum class WallTraceLaw { StrictWholeCell, PrescribedWetDry };
        enum class Kind { ResolvedCorner, NativeCoincident };
        GeometricVofIdentity sourceIdentity() const;
        GeometricVofTimeLayer sourceTime() const;
        Vector2D physicalPoint() const;
        int nativeVertex() const;
        int nativeCell() const; // Fluid cell opposite the excluded solid.
        Geometry::PlanarCellCut::RectangleCorner corner() const;
        int wetNativeRow() const;
        int dryNativeRow() const;
        Kind kind() const;
        WallTraceLaw wallTraceLaw() const;
        // Successful current-layer exterior phase declarations, not measured
        // minority support. False for strict/default or NativeCoincident.
        bool prescribedWetWallTrace() const;
        bool prescribedDryWallTrace() const;
        void validateFor(const QuadtreeMesh&,const GeometricVofPublication&) const;
    private:
        struct Data;
        std::shared_ptr<const Data> data_;
        explicit PinCertificate(std::shared_ptr<const Data>);
        void validateForInBatch(const QuadtreeMesh&,const GeometricVofPublication&,
            const MeshValidationBatch&) const;
        friend class GeometricVofPublication;
        friend class GeometricVofTransport;
    };
    static GeometricVofPublication withPinnedCorner(const QuadtreeMesh&,
        const GeometricVofPublication&,Vector2D physicalPoint,
        int actualWetNativeRow,int actualDryNativeRow,
        PinCertificate::WallTraceLaw=PinCertificate::WallTraceLaw::StrictWholeCell);
    std::optional<PinCertificate> pinnedCorner() const;
    // A narrow XY geometric qualification for a complete, isolated, closed
    // single-cell phase polygon smaller than the actual source leaf Delta.
    // It is NOT a curvature value, HF success, or a resolved-droplet claim.
    class IsolatedFragment
    {
    public:
        GeometricVofIdentity sourceIdentity() const;
        GeometricVofTimeLayer sourceTime() const;
        int nativeCell() const;
        unsigned phase() const; // 0 liquid, 1 gas; inventories remain unchanged.
        long double diameter() const; // Certified upper bound, strictly < delta().
        long double delta() const; // Actual min(dx,dy), never a stock threshold.
        std::size_t coveredNativeFaces() const;
        std::size_t checkedCornerCells() const;
        void validateFor(const QuadtreeMesh&,const GeometricVofPublication&) const;
    private:
        struct Data;
        std::shared_ptr<const Data> data_;
        explicit IsolatedFragment(std::shared_ptr<const Data>);
        friend class GeometricVofPublication;
    };
    // All actual split faces and their endpoint/corner incidence are checked.
    // Walls, RZ, unresolved source normals, uncertain/contact topology and
    // insufficient diameter margin return no qualification. Invalid native
    // index or a foreign mesh is an error. No geometry or stock is changed.
    std::optional<IsolatedFragment> isolatedLocalPlicFragment(const QuadtreeMesh&,int nativeCell) const;
    // A different sufficient condition: a complete closed component formed
    // by exactly two same-level, same-size, face-adjacent XY cells. Its AABB
    // fits two bins, each with L1 diameter < Delta, so no three component
    // positions can be pairwise Delta-separated. NOT a diameter certificate.
    class TwoCellPositionDeficit
    {
    public:
        GeometricVofIdentity sourceIdentity() const;
        GeometricVofTimeLayer sourceTime() const;
        int nativeCell() const; // Requested source; nativeCells() is canonical.
        const std::array<int,2>& nativeCells() const; // Increasing native index.
        unsigned phase() const;
        // Conservative AABB in the cell-centred frame of nativeCells()[0].
        std::array<long double,2> bboxLower() const;
        std::array<long double,2> bboxUpper() const;
        long double delta() const;
        long double twoBinL1Bound() const; // Upward bound, strictly < delta().
        unsigned splitAxis() const; // 0 splits bbox width, 1 splits height.
        unsigned maximumIndependentPositions() const; // At most two at Delta.
        std::size_t coveredNativeFaces() const;
        std::size_t checkedCornerCells() const;
        void validateFor(const QuadtreeMesh&,const GeometricVofPublication&) const;
    private:
        struct Data;
        std::shared_ptr<const Data> data_;
        explicit TwoCellPositionDeficit(std::shared_ptr<const Data>);
        friend class GeometricVofPublication;
    };
    // Full shared phase trace, complete actual external split-face coverage
    // and all corner incidences are required. Walls/RZ/uncertain topology
    // return nullopt. Stocks, PLIC planes and the singleton contract are kept.
    std::optional<TwoCellPositionDeficit> twoCellLocalPlicPositionDeficit(const QuadtreeMesh&,int nativeCell) const;
    // Complete published-geometry phase component carried by an actual
    // same-level interior vertex star (at most four equal-size leaves).
    // Strict whole-polygon diameter < Delta, not a physical droplet/stock claim.
    class VertexStarFragment
    {
    public:
        GeometricVofIdentity sourceIdentity() const;
        GeometricVofTimeLayer sourceTime() const;
        int nativeCell() const; // Requested source; nativeCells is canonical.
        const std::vector<int>& nativeCells() const;
        Geometry::PlanarCellCut::Point starVertex() const;
        unsigned phase() const;
        // Conservative whole-component AABB in nativeCells()[0]'s centre frame.
        std::array<long double,2> bboxLower() const;
        std::array<long double,2> bboxUpper() const;
        long double diameter() const; // Upward L1 bound, strictly < delta().
        long double delta() const;
        std::size_t coveredNativeFaces() const;
        std::size_t checkedCornerCells() const;
        void validateFor(const QuadtreeMesh&,const GeometricVofPublication&) const;
    private:
        struct Data;
        std::shared_ptr<const Data> data_;
        explicit VertexStarFragment(std::shared_ptr<const Data>);
        friend class GeometricVofPublication;
    };
    // Existing measured Empty/Full semantics retain all raw P0 stocks.
    // All actual split faces/corners and positive trace connectivity are proved;
    // wall/open/RZ/unknown/contact/external-connection cases return nullopt.
    std::optional<VertexStarFragment> vertexStarLocalPlicFragment(const QuadtreeMesh&,int nativeCell) const;
    static GeometricVofPublication prepare(const QuadtreeMesh&,NativeLeafDomain,
        GeometricVofTimeLayer,const std::vector<double>& alpha,const GeometricVofConfig& = {});
    // LocalPlic supports Core-selected whole leaves. Arrays remain native-
    // indexed at mesh.cellIndexNum(); domain.cellRowByIndex is the sole active
    // authority. Inactive entries have nativeCell=-1, zero capacity/stocks/
    // alpha/closure defect, no resolved geometry and no interface. Their zero
    // diagnostic alpha is an inactive sentinel, never a pure-gas assignment.
    // Nonzero inactive inputs are rejected. CompatibleFan remains complete-only.
    static GeometricVofPublication prepareLocalPlic(const QuadtreeMesh&,NativeLeafDomain,
        GeometricVofTimeLayer,const std::vector<double>& alpha,const GeometricVofConfig& = {});
    // Explicit independent liquid/gas stocks, native-cell indexed. Neither
    // phase is recovered from a rounded alpha or from its scalar complement.
    // Inactive native entries must be exactly {0,0}; only active leaves count.
    static GeometricVofPublication prepareLocalPlicInventory(const QuadtreeMesh&,NativeLeafDomain,
        GeometricVofTimeLayer,const std::vector<std::array<long double,2>>&,
        const GeometricVofConfig&);
    // Sealed native-ABI LocalPlic Initial/End image, issued only by this owner.
    // Full rectangular XY/RZ domains only. Restores raw independent stocks,
    // actual planes/measures and their exact capacity certificate unchanged;
    // does not reconstruct normals/geometry or manufacture a raw-stock root.
    // The target domain must have the identical actual native mesh/indexing.
    // Restored publications receive a fresh identity and retain the full clock.
    std::string checkpoint() const;
    static GeometricVofPublication restoreCheckpoint(const QuadtreeMesh&,NativeLeafDomain,
        const std::string&);
    GeometricVofGeometryKind geometryKind() const;
    GeometricVofLocalResolution localPlicResolution() const;
    // Dense native arrays: test the domain map before interpreting an entry.
    const std::vector<GeometricVofPlicCell>& plicCells() const;
    // Dense native no-resolved-capillary classification under the owner's
    // explicit NativeInterfaceSupport model. The fixed support is the union of
    // leaves incident to ALL actual source-boundary vertices (including self).
    // Each such vertex must have complete four-quadrant fluid coverage. A
    // strictly mixed source first lacks support if either phase has no raw
    // Qp>=Qother witness; ties witness both. A second, non-propagating pass
    // also marks opposite (native normal dot<=0) measured-mixed PLIC neighbours,
    // only when BOTH endpoints have complete support and a false first-pass
    // mask. First-pass tails never contribute their normals. A third,
    // simultaneous pass requires each phase to reach an independent raw
    // majority witness inside that same star through positive native-face
    // trace overlap. Measured-pure carriers retain their full phase; mixed
    // carriers must be eligible/resolved and outside the frozen second-pass
    // mask. Third-pass marks never propagate. Incomplete/wall support retains
    // the original method; pin/RZ reject the opt-in. Stocks,
    // actual PLIC planes, measured volumes and resolved flags are unchanged.
    // Both roles share the same support construction and classification.
    GeometricVofLocalResolutionMasks localPlicResolutionMasks(const QuadtreeMesh&) const;
    // Compatibility projection of that one query's final capillary mask.
    std::vector<bool> localPlicSubgridMask(const QuadtreeMesh&) const;
    // Exact QL+QG-nativeC, issued by a strict root or a typed predecessor
    // operation. Bound to this immutable mesh, identity and full clock.
    const std::vector<LAMP::Geometry::Detail::FiniteVolume::Dyadic>& capacityDefects() const;
    GeometricVofIdentity identity() const;
    GeometricVofTimeLayer timeLayer() const;
    const NativeLeafDomain& domain() const;
    const std::vector<double>& alpha() const;
    Geometry::PlanarCellCut::Point frameOrigin() const;
    const std::vector<Geometry::PlanarCellCut::AnchoredPoint>& anchoredVertices() const;
    const std::vector<Geometry::PlanarCellCut::Point>& vertices() const; // Lossy same-owner diagnostic projection.
    const std::vector<long double>& scalarValues() const; // Derived only, not transported.
    const std::vector<GeometricVofCell>& cells() const;
    // LocalPlic native rows with two excluded cells have no active incidence
    // or interface ranges; their nonzero mesh metric is not fluid coverage.
    const std::vector<GeometricVofFace>& faces() const;
    const std::vector<GeometricVofFacet>& facets() const;
    const GeometricVofTopology& topology() const;
    const GeometricVofReport& report() const;
    // Same authoritative actual fan; pure cells return the real native boundary.
    std::vector<std::vector<Geometry::PlanarCellCut::AnchoredPoint>> anchoredPhasePolygons(int nativeCell,bool liquid) const;
    // Lossy diagnostic projection of anchoredPhasePolygons; never a geometry predicate/measure.
    std::vector<std::vector<Geometry::PlanarCellCut::Point>> phasePolygons(int nativeCell,bool liquid) const;
    void validateMesh(const QuadtreeMesh&) const;
    void validateFor(const QuadtreeMesh&,GeometricVofIdentity,GeometricVofTimeLayer) const;
private:
    // Stack-only proof for one callback-free Boundary/Slab validation or
    // Transport validation/LocalPlic prepare, after the original complete
    // source validation. Never retained in a publication/slab/boundary.
    class MeshValidationBatch
    {
        friend class GeometricVofPublication;
        friend class GeometricVofTransport;
        friend class GeometricVofAdvectorSlab;
        MeshValidationBatch(const QuadtreeMesh& mesh,
            const std::shared_ptr<const NativeMeshIdentity>& identity)
            :mesh_(&mesh),identity_(identity){}
        MeshValidationBatch(const MeshValidationBatch&)=delete;
        MeshValidationBatch& operator=(const MeshValidationBatch&)=delete;
        const QuadtreeMesh* const mesh_;
        const std::shared_ptr<const NativeMeshIdentity> identity_;
    };
    void validateForInBatch(const QuadtreeMesh&,GeometricVofIdentity,
        GeometricVofTimeLayer,const MeshValidationBatch&) const;
    static void validateDomainInBatch(const QuadtreeMesh&,const NativeLeafDomain&,
        const MeshValidationBatch&);
    // Stack-only receipt for one LocalPlic sweep after its ORIGINAL complete
    // slab validation. Match both the original immutable parent and its binding;
    // this never authorizes a later public call or a different mesh/parent.
    class ConstructionMeshReceipt
    {
        friend class GeometricVofPublication;
        friend class GeometricVofTransport;
        ConstructionMeshReceipt(const QuadtreeMesh&,const GeometricVofPublication&);
        ConstructionMeshReceipt(const ConstructionMeshReceipt&)=delete;
        ConstructionMeshReceipt& operator=(const ConstructionMeshReceipt&)=delete;
        const QuadtreeMesh* const mesh_;
        const void* const parentData_;
        const std::shared_ptr<const NativeMeshIdentity> identity_;
    };
    // Callback-free capability for ONE actual mesh-transfer publication.
    // Only the real transfer can construct it from its immutable overlap;
    // no public raw-stock caller can self-sign finite-capacity defects.
    class MeshTransferReceipt
    {
        friend class GeometricVofPublication;
        friend class GeometricVofMeshTransfer;
        MeshTransferReceipt(const QuadtreeMesh&,const GeometricVofPublication&,
            const QuadtreeMesh&,const Mesh::NativeLeafOverlap&);
        MeshTransferReceipt(const MeshTransferReceipt&)=delete;
        MeshTransferReceipt& operator=(const MeshTransferReceipt&)=delete;
        const QuadtreeMesh* const sourceMesh_;
        const QuadtreeMesh* const targetMesh_;
        const void* const parentData_;
        const Mesh::NativeLeafOverlap* const coverage_;
    };
    enum class CertifiedOperation { MeshTransfer,DirectionalScratch,HalfTransport,FullTransport };
    // Only the real transfer/transport owners can derive a publication.
    // Same-mesh derivations retain the exact selected map. Selected-domain
    // MeshTransfer additionally requires its private actual-overlap receipt.
    // predictedDefects comes from predecessor/operation algebra; computing
    // QL+QG-C here validates that certificate and never authorizes its origin.
    static GeometricVofPublication prepareCertifiedLocalPlicInventory(
        const QuadtreeMesh&,NativeLeafDomain,GeometricVofTimeLayer,
        const GeometricVofPublication&,CertifiedOperation,
        const std::vector<std::array<long double,2>>&,
        const std::vector<LAMP::Geometry::Detail::FiniteVolume::Dyadic>&,
        const GeometricVofConfig&,const ConstructionMeshReceipt* =nullptr,
        const MeshTransferReceipt* =nullptr);
    struct Data;
    std::shared_ptr<const Data> data_;
    explicit GeometricVofPublication(std::shared_ptr<const Data>);
    const Data& checked() const;
    friend class GeometricVofTransport;
    friend class GeometricVofAdvectorSlab;
    friend class GeometricVofMeshTransfer;
};
}
