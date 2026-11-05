#pragma once

#include "Modules/Interface/HermiteGalsTrace.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Geometry/PolynomialVelocitySlab.h"
#include <array>
#include <memory>
#include <optional>
#include <vector>

namespace LAMP::Interface
{
class HermiteGalsIntegratedInventoryResult;
class HermiteGalsCertifiedInventoryResult;
class HermiteGalsConstantInflowSource;
class HermiteGalsBoundaryPhaseSource;
class HermiteGalsInventoryRegionSource;
struct HermiteGalsIntegratedInventoryConfig;
// Source identity is the retained immutable owner, as in the existing native
// field/phase contracts. Copies share identity; independently bound callbacks
// never become equal merely because geometry/time or sampled values agree.
// Captures have the SAME immutable, owning/repeatable-provider prerequisite as
// HermiteGalsVelocityView; this API cannot certify arbitrary mutable closures.
class HermiteGalsTemporalVelocitySource
{
public:
    // Bind queries through the existing validated endpoint-view constructor;
    // do not adopt externally writable handles or moved-from velocity views.
    HermiteGalsTemporalVelocitySource(const HermiteGalsState& previousAccepted,
        HermiteGalsIdentity arrivalGeometry, HermiteGalsTimeLayer arrivalTime,
        HermiteGalsVelocityView::Query previous, HermiteGalsVelocityView::Query arrival);
    // Own the value AND generate the actual two endpoint queries. An opaque
    // query cannot acquire this capability by agreement of sampled values.
    static HermiteGalsTemporalVelocitySource constant(const HermiteGalsState& previousAccepted,
        HermiteGalsIdentity arrivalGeometry, HermiteGalsTimeLayer arrivalTime, Vector2D value);
    static HermiteGalsTemporalVelocitySource polynomial(const HermiteGalsState& previousAccepted,
        HermiteGalsIdentity arrivalGeometry,HermiteGalsTimeLayer arrivalTime,
        const LAMP::Geometry::PolynomialVelocityData& data);
    // Null for opaque/constant sources. The returned immutable object is the
    // SAME owner captured by the actual endpoint/substage callbacks.
    const LAMP::Geometry::PolynomialVelocitySlab* polynomialVelocity() const;
    HermiteGalsVelocityView polynomialArrivalView(HermiteGalsIdentity actualCandidate,
        HermiteGalsTimeLayer actualArrival) const;
    std::optional<Vector2D> constantValue() const;
    // Actual substage view of the SAME owned constant provider. The new
    // candidate/time is explicit; opaque callbacks cannot claim this ability.
    // Keeps this immutable source alive without manufacturing a second source.
    HermiteGalsVelocityView constantArrivalView(HermiteGalsIdentity actualCandidate,
        HermiteGalsTimeLayer actualArrival) const;
    const HermiteGalsVelocityView& previous() const;
    const HermiteGalsVelocityView& arrival() const;
    bool sameSource(const HermiteGalsTemporalVelocitySource&) const noexcept;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
    const Data& checked() const;
};

// Owns the actual domain/coverage callback used by the EXISTING fixed RK2
// tracer, plus its exact native mesh identity. No alternate tracer or numeric
// source token. Construct once and share across this trajectory's intervals.
class HermiteGalsTemporalCharacteristicSource
{
public:
    enum class Method { CoherentExplicitMidpointRK2 };
    HermiteGalsTemporalCharacteristicSource(const QuadtreeMesh& mesh,
        const LevelSetBoundaryTrace::Domain& domain);
    const LevelSetBoundaryTrace::Domain& domain() const;
    const NativeMeshIdentity& meshIdentity() const;
    Method method() const;
    void validateMesh(const QuadtreeMesh& mesh) const;
    bool sameSource(const HermiteGalsTemporalCharacteristicSource&) const noexcept;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
    const Data& checked() const;
};

// Explicit [begin,end] extent; no fabricated Midpoint geometry ID. The App
// supplies its accepted Initial/End handle. That immutable authority is kept
// alive, not copied/rebuilt. A subinterval may start later than the old state,
// but remains inside the retained velocity pair's full physical time span.
class HermiteGalsPhaseFluxIntervalSource
{
public:
    HermiteGalsPhaseFluxIntervalSource(const QuadtreeMesh& mesh,
        const HermiteGalsState& previousAccepted,
        const HermiteGalsTemporalVelocitySource& velocity,
        const HermiteGalsTemporalCharacteristicSource& characteristic,
        double timeBegin, double timeEnd);
    HermiteGalsPhaseFluxIntervalSource(const QuadtreeMesh& mesh,
        const HermiteGalsState& previousAccepted,
        const HermiteGalsTemporalVelocitySource& velocity,
        const HermiteGalsTemporalCharacteristicSource& characteristic,
        double timeBegin,double timeEnd,const HermiteGalsConstantInflowSource& inflow);
    HermiteGalsPhaseFluxIntervalSource(const QuadtreeMesh& mesh,
        const HermiteGalsState& previousAccepted,
        const HermiteGalsTemporalVelocitySource& velocity,
        const HermiteGalsTemporalCharacteristicSource& characteristic,
        double timeBegin,double timeEnd,const HermiteGalsBoundaryPhaseSource& boundary);
    const HermiteGalsState& previousAccepted() const;
    HermiteGalsIdentity previousGeometryIdentity() const;
    const HermiteGalsTemporalVelocitySource& velocitySource() const;
    const HermiteGalsTemporalCharacteristicSource& characteristicSource() const;
    const HermiteGalsConstantInflowSource* inflowSource() const;
    const HermiteGalsBoundaryPhaseSource* boundaryPhaseSource() const;
    const NativeMeshIdentity& meshIdentity() const;
    HermiteGalsMeasure measure() const;
    double timeBegin() const;
    double timeEnd() const;
    double dt() const;
    // Equality of the complete provenance/extent, not a generated numeric ID.
    bool sameInterval(const HermiteGalsPhaseFluxIntervalSource&) const;
    void validateMesh(const QuadtreeMesh& mesh) const;
    NativeLeafFace face(const QuadtreeMesh& mesh, int nativeFace) const;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
    const Data& checked() const;
};

// Read-only proof of a complete no-hit RK2 branch over this source's arrival
// time interval and spatial box. It retains provenance, not a fake time-node
// geometry. No phase classification, integrated flux or commit is implied.
class HermiteGalsCharacteristicEnclosure
{
public:
    const HermiteGalsPhaseFluxIntervalSource& interval() const;
    const LAMP::Geometry::PolynomialRk2Bounds& bounds() const;
    void validateMesh(const QuadtreeMesh&) const;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
    explicit HermiteGalsCharacteristicEnclosure(std::shared_ptr<const Data>);
    const Data& checked() const;
    friend HermiteGalsCharacteristicEnclosure encloseHermiteGalsCharacteristic(const QuadtreeMesh&,
        const HermiteGalsPhaseFluxIntervalSource&,Vector2D,Vector2D);
};
HermiteGalsCharacteristicEnclosure encloseHermiteGalsCharacteristic(const QuadtreeMesh&,
    const HermiteGalsPhaseFluxIntervalSource&,Vector2D lower,Vector2D upper);

// Strict sign on the entire closed region. Unresolved proves neither an
// interface intersection nor its absence; it must never be treated as gas.
enum class HermiteGalsCertifiedPhase { Liquid=-1, Unresolved=0, Gas=1 };
struct HermiteGalsCharacteristicPhaseCell
{
    int cell=-1;
    Vector2D lower,upper; // closed intersection with this OLD native leaf
    HermiteGalsCertifiedPhase phase=HermiteGalsCertifiedPhase::Unresolved;
    long double minimumLiquidContent=0,maximumLiquidContent=1;
};
// Complete closed departure-box coverage on the same old Q4 snapshot.
// Optional r18 inventory weights describe bounded content on that geometry;
// a pure geometric phase is not a claim that the independent M is 0 or C.
class HermiteGalsCharacteristicPhaseCoverage
{
public:
    const HermiteGalsCharacteristicEnclosure& characteristic() const;
    const std::vector<HermiteGalsCharacteristicPhaseCell>& cells() const;
    HermiteGalsCertifiedPhase phase() const;
    long double minimumLiquidContent() const;
    long double maximumLiquidContent() const;
    const HermiteGalsInventoryRegionSource* inventorySource() const;
    void validateMesh(const QuadtreeMesh&) const;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
    explicit HermiteGalsCharacteristicPhaseCoverage(std::shared_ptr<const Data>);
    const Data& checked() const;
    friend HermiteGalsCharacteristicPhaseCoverage coverHermiteGalsCharacteristicPhase(const QuadtreeMesh&,
        const HermiteGalsCharacteristicEnclosure&);
    friend HermiteGalsCharacteristicPhaseCoverage coverHermiteGalsCharacteristicPhase(const QuadtreeMesh&,
        const HermiteGalsCharacteristicEnclosure&,const HermiteGalsInventoryRegionSource&);
};
HermiteGalsCharacteristicPhaseCoverage coverHermiteGalsCharacteristicPhase(const QuadtreeMesh&,
    const HermiteGalsCharacteristicEnclosure&);
HermiteGalsCharacteristicPhaseCoverage coverHermiteGalsCharacteristicPhase(const QuadtreeMesh&,
    const HermiteGalsCharacteristicEnclosure&,const HermiteGalsInventoryRegionSource&);

enum class HermiteGalsInflowPhase { Gas, Liquid };
// Explicit whole-side, time-independent physical phase history for an owned
// polynomial velocity. Array order west/south/east/north. The unbound base
// interval keeps ownership acyclic. This data alone certifies NO trace branch.
class HermiteGalsBoundaryPhaseSource
{
public:
    HermiteGalsBoundaryPhaseSource(const QuadtreeMesh&,const HermiteGalsPhaseFluxIntervalSource& baseExtent,
        const std::array<HermiteGalsInflowPhase,4>& phases);
    const HermiteGalsPhaseFluxIntervalSource& baseExtent() const;
    HermiteGalsInflowPhase phaseFor(Boundary side) const;
    void validateInterval(const QuadtreeMesh&,const HermiteGalsPhaseFluxIntervalSource&) const;
    bool sameSource(const HermiteGalsBoundaryPhaseSource&) const noexcept;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
    const Data& checked() const;
};
struct HermiteGalsInflowFacePhase
{
    int nativeFace=-1;
    HermiteGalsInflowPhase phase=HermiteGalsInflowPhase::Gas;
};
// Owned time-independent, whole-face inflow data; never inferred from old phi.
// Base extent MUST have no BC attached (acyclic ownership). Entries name ALL
// actual inflow faces in ascending order for this owned single-axis constant
// velocity. XY or RZ axial, full active domain, no opaque coverage. No spatial
// phase partition, time interpolation or extrapolated ghost scalar is implied.
class HermiteGalsConstantInflowSource
{
public:
    HermiteGalsConstantInflowSource(const QuadtreeMesh&,const HermiteGalsPhaseFluxIntervalSource& baseExtent,
        const std::vector<HermiteGalsInflowFacePhase>& entries);
    const HermiteGalsPhaseFluxIntervalSource& baseExtent() const;
    const std::vector<HermiteGalsInflowFacePhase>& entries() const;
    HermiteGalsInflowPhase phaseFor(int nativeFace) const;
    void validateInterval(const QuadtreeMesh&,const HermiteGalsPhaseFluxIntervalSource&) const;
    bool sameSource(const HermiteGalsConstantInflowSource&) const noexcept;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
    const Data& checked() const;
};

struct HermiteGalsIntegratedFaceVolumes
{
    // XY: area per unit depth; RZ: true physical volume (including 2*pi*r).
    // Inputs are already integrated in space AND time. No rate is stored.
    long double integratedTotalPositiveVolume = 0;
    long double integratedTotalNegativeVolume = 0;
    long double integratedLiquidPositiveVolume = 0;
    long double integratedLiquidNegativeVolume = 0;
};

enum class HermiteGalsFaceTransferOrientation { PositiveX, PositiveY };

// Pure algebraic contract, not a certificate that any quadrature has run.
// Canonical native +x/+y orientation is shared by both adjacent cells.
class HermiteGalsIntegratedFaceTransfer
{
public:
    HermiteGalsIntegratedFaceTransfer(const QuadtreeMesh& mesh,
        const HermiteGalsPhaseFluxIntervalSource& interval, int nativeFace,
        HermiteGalsFaceTransferOrientation orientation,
        const HermiteGalsIntegratedFaceVolumes& volumes);
    const HermiteGalsPhaseFluxIntervalSource& interval() const;
    const NativeLeafFace& face() const;
    HermiteGalsFaceTransferOrientation orientation() const;
    const HermiteGalsIntegratedFaceVolumes& integratedVolumes() const;
    long double signedTotalVolume() const;
    long double signedLiquidVolume() const;
    long double averageTotalFlux() const;
    long double averageLiquidFlux() const;
    // Chronological, contiguous composition only. Does not sample, integrate,
    // update M, or mutate either operand, even when validation/addition fails.
    HermiteGalsIntegratedFaceTransfer appended(const QuadtreeMesh& mesh,
        const HermiteGalsIntegratedFaceTransfer& next) const;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
    const Data& checked() const;
    friend HermiteGalsIntegratedInventoryResult assembleHermiteGalsIntegratedInventory(const QuadtreeMesh&,
        const HermiteGalsPhaseFluxIntervalSource&, const std::vector<HermiteGalsIntegratedFaceTransfer>&,
        const HermiteGalsIntegratedInventoryConfig&);
    friend class HermiteGalsIntegratedInventoryResult;
};

struct HermiteGalsSpaceTimeIntegralConfig
{
    double absoluteVolumeTolerance=1e-8; // four-width L1 budget, not an M repair tolerance
    std::size_t maximumPanels=4096;
    int maximumDepth=40;
};
struct HermiteGalsSpaceTimeIntegralReport
{
    double absoluteErrorBound=0;
    std::size_t acceptedPanels=0,evaluatedPanels=0,phaseQueries=0,branchRefusals=0,unresolvedPhasePanels=0;
    std::size_t initialBoundaryQueries=0; // certified trial panels, including discarded trials
    std::size_t uniformReachableQueries=0; // constant source content, NOT a B/J completion certificate
    int maximumDepthReached=0;
};
class HermiteGalsBoundedFaceIntegral
{
public:
    const HermiteGalsIntegratedFaceTransfer& transfer() const;
    const HermiteGalsSpaceTimeIntegralConfig& config() const;
    const LAMP::Geometry::PolynomialFaceIntegralBounds& bounds() const;
    const HermiteGalsSpaceTimeIntegralReport& report() const;
    const HermiteGalsInventoryRegionSource* inventorySource() const;
    void validateMesh(const QuadtreeMesh&) const;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
    explicit HermiteGalsBoundedFaceIntegral(std::shared_ptr<const Data>);
    const Data& checked() const;
    friend HermiteGalsBoundedFaceIntegral integrateHermiteGalsPolynomialFace(const QuadtreeMesh&,
        const HermiteGalsPhaseFluxIntervalSource&,int,const HermiteGalsSpaceTimeIntegralConfig&);
    friend HermiteGalsBoundedFaceIntegral integrateHermiteGalsPolynomialFace(const QuadtreeMesh&,
        const HermiteGalsPhaseFluxIntervalSource&,int,const HermiteGalsInventoryRegionSource&,const HermiteGalsSpaceTimeIntegralConfig&);
};
HermiteGalsBoundedFaceIntegral integrateHermiteGalsPolynomialFace(const QuadtreeMesh&,
    const HermiteGalsPhaseFluxIntervalSource&,int nativeFace,const HermiteGalsSpaceTimeIntegralConfig& = {});
HermiteGalsBoundedFaceIntegral integrateHermiteGalsPolynomialFace(const QuadtreeMesh&,
    const HermiteGalsPhaseFluxIntervalSource&,int nativeFace,const HermiteGalsInventoryRegionSource&,
    const HermiteGalsSpaceTimeIntegralConfig& = {});

// Narrow exact-event operator: affine retained Q4 throughout the swept box,
// constant owned velocity, complete interior coverage. RZ requires axial
// translation. No candidate geometry, inventory update or midpoint fallback.
HermiteGalsIntegratedFaceTransfer integrateHermiteGalsPlanarFace(
    const QuadtreeMesh& mesh, const HermiteGalsPhaseFluxIntervalSource& interval, int nativeFace);

class HermiteGalsPhaseRegionSource;
class HermiteGalsInventoryRegionSource;
struct HermiteGalsSweptFaceConfig
{
    double relativeCoordinateTolerance=1e-10;
    double relativeMomentEstimateTolerance=1e-10;
};
struct HermiteGalsSweptFaceReport
{
    long double coordinateVolumeBand=0,momentErrorEstimate=0;
    long double innerLiquid=0,outerLiquid=0,innerGas=0,outerGas=0;
    std::size_t rectangleQueries=0,integratedArcs=0,rootRefinements=0;
    bool inventoryWeighted=false;
    // Same rectangle queries; geometricLiquidVolume is an unsigned reference,
    // not a second production face flux. Correction has global face orientation.
    long double geometricLiquidVolume=0,signedInventoryCorrection=0,l1InventoryCorrection=0;
};
struct HermiteGalsSweptFaceResult
{
    HermiteGalsIntegratedFaceTransfer transfer;
    HermiteGalsSweptFaceReport report;
};
// Complete old-Q4 geometry in the normal-translation swept rectangle.
// Owned constant velocity, zero tangential velocity, full active support and
// interior characteristics only; RZ allows axial motion. Inner/outer spatial
// queries cover coordinate rounding; their floating moments are not certified
// bounds. No candidate geometry, BC extrapolation or inventory modification.
HermiteGalsSweptFaceResult integrateHermiteGalsSweptFace(const QuadtreeMesh& mesh,
    const HermiteGalsPhaseFluxIntervalSource& interval,int nativeFace,
    const HermiteGalsPhaseRegionSource& regions,const HermiteGalsSweptFaceConfig& config = {});
// Same sweep/coverage contract, integrating the old source's bounded inventory
// content. The caller explicitly chooses this producer; no automatic fallback.
HermiteGalsSweptFaceResult integrateHermiteGalsSweptFace(const QuadtreeMesh& mesh,
    const HermiteGalsPhaseFluxIntervalSource& interval,int nativeFace,
    const HermiteGalsInventoryRegionSource& regions,const HermiteGalsSweptFaceConfig& config = {});

struct HermiteGalsTranslationFaceResult
{
    HermiteGalsIntegratedFaceTransfer transfer;
    bool prescribedInflow=false;
    std::optional<HermiteGalsSweptFaceReport> swept;
    bool inventoryWeighted=false; // also identifies prescribed-inflow/zero-velocity results
};
// Complete single-axis constant-translation face consumer: actual inflow
// boundaries require the explicit source attached to interval; other faces
// use the existing swept query. Interior traces leaving the domain still
// reject. No inventory update or new geometry is performed here.
HermiteGalsTranslationFaceResult integrateHermiteGalsTranslationFace(const QuadtreeMesh& mesh,
    const HermiteGalsPhaseFluxIntervalSource& interval,int nativeFace,
    const HermiteGalsPhaseRegionSource& regions,const HermiteGalsSweptFaceConfig& config = {});
HermiteGalsTranslationFaceResult integrateHermiteGalsTranslationFace(const QuadtreeMesh& mesh,
    const HermiteGalsPhaseFluxIntervalSource& interval,int nativeFace,
    const HermiteGalsInventoryRegionSource& regions,const HermiteGalsSweptFaceConfig& config = {});

struct HermiteGalsIntegratedInventoryConfig
{
    double relativeDivergenceTolerance=1e-12;
    double relativeConservationTolerance=1e-12;
};
struct HermiteGalsInventoryAssemblyReport
{
    std::size_t activeFaces=0,internalFaces=0,counterflowFaces=0;
    std::size_t roundoffRepairSweeps=0,repairedFaces=0;
    long double liquidBefore=0,liquidAfter=0,boundaryOutward=0,closure=0;
    long double maxOutflowCfl=0,maxRelativeDivergence=0,maxRelativeFluxCorrection=0;
};
class HermiteGalsIntegratedInventoryResult
{
public:
    const HermiteGalsPhaseFluxIntervalSource& interval() const;
    // Raw integrated transfers and final shared transfers retain the SAME
    // interval source. Only the existing bounded interior roundoff repair
    // may alter final liquid volumes; total volumes never change.
    const std::vector<HermiteGalsIntegratedFaceTransfer>& rawFaces() const;
    const std::vector<HermiteGalsIntegratedFaceTransfer>& faces() const;
    const std::vector<HermiteGalsCellVolume>& phaseVolumes() const;
    const HermiteGalsInventoryAssemblyReport& report() const;
    void validateMesh(const QuadtreeMesh& mesh) const;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
    explicit HermiteGalsIntegratedInventoryResult(std::shared_ptr<const Data>);
    const Data& checked() const;
    static HermiteGalsIntegratedInventoryResult assemble(const QuadtreeMesh&,
        const HermiteGalsPhaseFluxIntervalSource&,const std::vector<HermiteGalsIntegratedFaceTransfer>&,
        const HermiteGalsIntegratedInventoryConfig&,const std::vector<HermiteGalsBoundedFaceIntegral>*);
    friend HermiteGalsIntegratedInventoryResult assembleHermiteGalsIntegratedInventory(const QuadtreeMesh&,
        const HermiteGalsPhaseFluxIntervalSource&, const std::vector<HermiteGalsIntegratedFaceTransfer>&,
        const HermiteGalsIntegratedInventoryConfig&);
    friend HermiteGalsCertifiedInventoryResult assembleHermiteGalsCertifiedInventory(const QuadtreeMesh&,
        const HermiteGalsPhaseFluxIntervalSource&,const std::vector<HermiteGalsBoundedFaceIntegral>&,
        const HermiteGalsIntegratedInventoryConfig&);
};
// Algebraic assembly, not certification of arbitrary externally supplied
// quadrature or boundary history. Complete, ordered, compatible native faces
// are required. No rate->volume reconstruction, no new geometry, no commit.
HermiteGalsIntegratedInventoryResult assembleHermiteGalsIntegratedInventory(const QuadtreeMesh& mesh,
    const HermiteGalsPhaseFluxIntervalSource& interval,
    const std::vector<HermiteGalsIntegratedFaceTransfer>& faces,
    const HermiteGalsIntegratedInventoryConfig& config = {});

struct HermiteGalsCertifiedInventoryReport
{
    bool selected=false;
    std::size_t changedFaces=0,augmentations=0;
    long double maximumRelativeSelectionChange=0;
};
// One batch selection inside the original certified liquid intervals, followed
// by the shared inventory assembler, with certificate-constrained roundoff
// repair inside the original budget. All original certificates stay owned.
class HermiteGalsCertifiedInventoryResult
{
public:
    const HermiteGalsIntegratedInventoryResult& inventory() const;
    const std::vector<HermiteGalsBoundedFaceIntegral>& certificates() const;
    const HermiteGalsCertifiedInventoryReport& report() const;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
    explicit HermiteGalsCertifiedInventoryResult(std::shared_ptr<const Data>);
    const Data& checked() const;
    friend HermiteGalsCertifiedInventoryResult assembleHermiteGalsCertifiedInventory(const QuadtreeMesh&,
        const HermiteGalsPhaseFluxIntervalSource&,const std::vector<HermiteGalsBoundedFaceIntegral>&,
        const HermiteGalsIntegratedInventoryConfig&);
};
HermiteGalsCertifiedInventoryResult assembleHermiteGalsCertifiedInventory(const QuadtreeMesh&,
    const HermiteGalsPhaseFluxIntervalSource&,const std::vector<HermiteGalsBoundedFaceIntegral>&,
    const HermiteGalsIntegratedInventoryConfig& = {});
}
