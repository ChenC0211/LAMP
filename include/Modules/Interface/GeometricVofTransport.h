#pragma once
#include "Modules/Interface/GeometricVofPublication.h"
#include <array>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace LAMP::Interface
{
class VofSmoothGeometry;
struct VofSmoothGeometryConfig;
// Actual LocalPlic donor-column or complete-cell Courant excess only.
// Configuration, owner, inventory and geometry errors remain hard failures.
class GeometricVofCourantExceeded : public std::invalid_argument
{
public:
    explicit GeometricVofCourantExceeded(const std::string& reason)
        :std::invalid_argument("[Interface/GeometricVofTransport] "+reason){}
};
// Immutable native-face temporal INPUT snapshot, not an accepted velocity
// state. It is bound to the actual old VOF owner and one complete next clock.
// Normal means are in canonical +x/r or +y/z direction, native-row indexed.
class GeometricVofAdvectorSlab
{
public:
    // Immutable actual-source boundary rules; no exterior ghost inventories.
    // Selected domains retain dense native rows: excluded/excluded faces are
    // Closed with exactly zero normals and flux, never pure-gas reservoirs.
    class Boundary
    {
    public:
        enum class Mode {Interior,Closed,PureLiquidReservoir,PureGasReservoir};
        static Boundary prepare(const QuadtreeMesh&,const GeometricVofPublication&,
            const std::vector<Mode>& rulesByNativeRow);
        GeometricVofIdentity sourceIdentity() const;
        GeometricVofTimeLayer sourceTime() const;
        const NativeLeafDomain& domain() const;
        const std::vector<Mode>& rules() const;
        bool allClosed() const;
        void validateFor(const QuadtreeMesh&,const GeometricVofPublication&) const;
        bool sameInstance(const Boundary&) const noexcept;
    private:
        void validateForInBatch(const QuadtreeMesh&,const GeometricVofPublication&,
            const GeometricVofPublication::MeshValidationBatch&) const;
        friend class GeometricVofAdvectorSlab;
        struct Data;
        std::shared_ptr<const Data> data_;
        explicit Boundary(std::shared_ptr<const Data>);
    };
    // Explicit LocalPlic reservoir slab. Boundary rules are source-bound before
    // dt selection; nextEnd still authenticates the complete next interval.
    // Exterior endpoint signs may agree or include zero, never reverse in-slab.
    // All-Closed delegates to prepare and retains its nullptr boundary.
    static GeometricVofAdvectorSlab prepareBoundary(const QuadtreeMesh&,
        const GeometricVofPublication&,const Boundary&,GeometricVofTimeLayer nextEnd,
        const std::vector<double>& beginNormal,const std::vector<double>& endNormal);
    const Boundary* boundary() const;
    static GeometricVofAdvectorSlab prepare(const QuadtreeMesh&,
        const GeometricVofPublication&,GeometricVofTimeLayer nextEnd,
        const std::vector<double>& beginNormal,const std::vector<double>& endNormal);
    GeometricVofIdentity sourceIdentity() const;
    GeometricVofTimeLayer sourceTime() const;
    GeometricVofTimeLayer endTime() const;
    const NativeLeafDomain& domain() const;
    const std::vector<double>& beginNormal() const;
    const std::vector<double>& endNormal() const;
    const std::vector<Geometry::PlanarCellCut::Point>& beginVelocity() const;
    const std::vector<Geometry::PlanarCellCut::Point>& endVelocity() const;
    void validateFor(const QuadtreeMesh&,const GeometricVofPublication&) const;
    bool sameInstance(const GeometricVofAdvectorSlab&) const noexcept;
private:
    void validateForInBatch(const QuadtreeMesh&,const GeometricVofPublication&,
        const GeometricVofPublication::MeshValidationBatch&) const;
    friend class GeometricVofTransport;
    static GeometricVofAdvectorSlab prepareImpl(const QuadtreeMesh&,
        const GeometricVofPublication&,GeometricVofTimeLayer,
        const std::vector<double>&,const std::vector<double>&,const Boundary*);
    struct Data;
    std::shared_ptr<const Data> data_;
    explicit GeometricVofAdvectorSlab(std::shared_ptr<const Data>);
};
struct GeometricVofTransportConfig
{
    long double maximumOutflowCfl=.9L; // CompatibleFan gross limit; LocalPlic uses maximumLocalCourant.
    long double relativeConservationTolerance=1e-12L;
    unsigned maximumTemporalEvents=512;
    GeometricVofConfig geometry;
    long double maximumLocalCourant=.5L; // LocalPlic: incoming, outgoing, donor-column and sealed-pin peak-depth bounds.
};
// Instantaneous complete LocalPlic native-face rates, in inverse code time.
// These are dimensional rates, not accepted full/half-interval Courant values.
struct GeometricVofLocalCourantRates
{
    long double maximumLocalColumn=0,maximumIncoming=0,maximumOutflow=0;
};
struct GeometricVofIntegratedFlux
{
    long double total=0,positiveTotal=0,negativeTotal=0; // Negative magnitude >=0.
    // CompatibleFan signed phases are once-rounded internal exact gross differences;
    // subtracting its exported gross diagnostics can round differently.
    // LocalPlic exports the declared frozen-mean directional swept integrals.
    long double rawLiquid=0,liquid=0,gas=0,mass=0;
    long double positiveLiquid=0,negativeLiquid=0;
    long double positiveGas=0,negativeGas=0; // Independent minor phase is retained.
    long double limiter=1; // Rounded diagnostic of the internally shared exact ratio.
    // Native accepted-face correction from a strict whole-phase sweep certificate,
    // distinct from raw geometry and limiter. Numerical subface traces do not
    // carry this geometry certificate; their physical accepted phases remain authoritative.
    long double exhaustionLiquidCorrection=0,exhaustionGasCorrection=0;
};
struct GeometricVofTransportFace
{
    int nativeRow=-1,negativeCell=-1,positiveCell=-1;
    GeometricVofIntegratedFlux half,full;
};
// Boundary physical phase balance is independent of the WY/RN algebraic
// receipt. All entries use the actual combined half/full directional result.
struct GeometricVofBoundaryLedger
{
    bool enabled=false;
    GeometricVofIdentity sourceIdentity;
    GeometricVofTimeLayer sourceTime,targetTime;
    std::array<long double,2> source{},target{},inflow{},outflow{};
    std::array<long double,2> compression{},roundingError{};
    // Physical source from canonical net = RN(positive-negative):
    // physicalResidual = compression + roundingError + netFluxRounding
    //                  + algebraicResidual. No inventory is repaired.
    std::array<long double,2> netFluxRounding{};
    std::array<long double,2> physicalResidual{},algebraicResidual{};
    long double physicalRelativeConservation=0,algebraicRelativeConservation=0;
    long double massResidual=0,massRelativeConservation=0;
};
// One prepareLocalPlicSweep invocation, including half/full and both directions.
// Pure execution diagnostics; no retained field, geometry or State owner.
struct GeometricVofLocalFaceExecutionReport
{
    bool requested=false,eligible=false;
    bool faceValuesParallelUsed=false,cellGatherParallelUsed=false,fallback=false;
    std::size_t nativeFaces=0,incidences=0,directionsUsed=0,faceEvaluations=0,cellGathers=0;
    std::size_t scratchBytes=0,scratchLimitBytes=64u*1024u*1024u;
    unsigned maximumWorkers=0;
};
struct GeometricVofTransportReport
{
    std::size_t donorQueries=0,polygonIntegrals=0,temporalIntervals=0,limitedFaces=0;
    std::size_t temporalEvents=0;
    std::size_t localDirectionalSweeps=0,localStripClips=0,localUnresolvedDonors=0;
    std::size_t localPinnedSweeps=0; // Actual nonzero sealed-corner triangular donor sweeps.
    std::size_t localWholePhaseExhaustions=0; // Strict actual-geometry certificates, before shared ledger.
    // Actual returned serial phaseUpdate calls. Production uses original Core.
    // Fallback counts every returned Core call, including its internal Noop;
    // Fast/Promoted/Noop are retained legacy filtered receipts and remain zero.
    // Core preserves its original old-stock signed-zero early return.
    std::size_t localPhaseUpdateCalls=0,localPhaseUpdateFast=0,localPhaseUpdatePromoted=0;
    std::size_t localPhaseUpdateFallback=0,localPhaseUpdateNoop=0;
    std::size_t localNativeP0Restrictions=0,localInteriorP0Extensions=0;
    long double maximumLocalColumnCfl=0,maximumIncomingCfl=0;
    long double maximumPinnedPeakCfl=0; // Actual peak normal depth / donor-column width.
    long double maximumWholePhaseExhaustionRelativeCorrection=0; // Scaled by the consumed phase stock.
    long double maximumCompressionResidual=0,maximumLocalStripResidual=0;
    unsigned maximumTemporalRefinements=0;
    long double maximumFaceTimeRelativeUncertainty=0; // Relative to a proven positive branch, never total capacity.
    long double maximumOutflowCfl=0,maximumTotalDivergence=0;
    long double maximumCoverageResidual=0,maximumFluxCorrection=0;
    long double maximumCapacityReconciliation=0,maximumInventoryResidual=0;
    long double halfRelativeConservation=0,fullRelativeConservation=0;
    GeometricVofBoundaryLedger boundaryHalf,boundaryFull;
    GeometricVofLocalFaceExecutionReport localFaceExecution;
};
// Actual world-coordinate subinterval of one native face. Endpoint order
// does not change canonical flux signs. Zero length is an exact empty query.
struct GeometricVofNativeSubface
{
    int nativeRow=-1;
    Geometry::PlanarCellCut::Point begin{},end{};
};
// Actual world interval strictly inside one original native pressure cell.
// It has no native parent face or pressure-inventory limiter demand.
struct GeometricVofInteriorSubface
{
    int nativeCell=-1; // Actual mesh native-cell index, not a compact row.
    Geometry::PlanarCellCut::Point begin{},end{};
};
struct GeometricVofInteriorFlux
{
    int nativeCell=-1;
    Geometry::PlanarCellCut::Point beginVelocity{},endVelocity{};
    GeometricVofIntegratedFlux half,full;
};
// LocalPlic requested-full results carry no unrequested half flux or endpoint
// velocity samples. Native incidence remains the actual accepted face's.
struct GeometricVofNativeFullFlux
{
    int nativeRow=-1,negativeCell=-1,positiveCell=-1;
    GeometricVofIntegratedFlux full;
};
struct GeometricVofInteriorFullFlux
{
    int nativeCell=-1;
    GeometricVofIntegratedFlux full;
};
// CompatibleFan prepare shares one simultaneous update. LocalPlic uses the
// explicitly named split sweep. Both methods start half/full from the SAME
// old owner; a C stage never advances the A candidate.
// This prepares candidates only, without an App commit or Hydro solve.
class GeometricVofTransport
{
public:
    // Read-only actual native-row means. No reconstruction, time slab,
    // transport candidate or divergence projection is prepared. Nonzero
    // exterior volume flow is unsupported, as in the closed local sweep.
    static GeometricVofLocalCourantRates localPlicCourantRates(const QuadtreeMesh&,
        const GeometricVofPublication&,const std::vector<double>& nativeNormal);
    // Instantaneous source-bound open rates; legal before selecting next dt.
    static GeometricVofLocalCourantRates localPlicBoundaryCourantRates(const QuadtreeMesh&,
        const GeometricVofPublication&,const GeometricVofAdvectorSlab::Boundary&,
        const std::vector<double>& nativeNormal);
    static GeometricVofTransport prepare(const QuadtreeMesh&,
        const GeometricVofPublication&,const GeometricVofAdvectorSlab&,
        long double liquidDensity,long double gasDensity,
        const GeometricVofTransportConfig& = {});
    // Explicit local split candidate. Half/full each run both directions from
    // the same old inventory; frozen old cc, with scratch PLIC reconstruction
    // between directions. Odd next step x/y, even y/x. This does not establish
    // second-order A+C compatibility. Material consumers use its retained
    // accepted flux P0 traces, not an additional sharp PLIC transport.
    static GeometricVofTransport prepareLocalPlicSweep(const QuadtreeMesh&,
        const GeometricVofPublication&,const GeometricVofAdvectorSlab&,
        long double liquidDensity,long double gasDensity,
        const GeometricVofTransportConfig& = {});
    // Explicit B1 source geometry. Each directional donor strip consumes its
    // actual smooth curve; the second direction rebuilds real scratch stocks.
    static GeometricVofTransport prepareLocalSmoothSweep(const QuadtreeMesh&,
        const GeometricVofPublication&,const VofSmoothGeometry&,
        const GeometricVofAdvectorSlab&,long double liquidDensity,long double gasDensity,
        const GeometricVofTransportConfig&,const VofSmoothGeometryConfig&);
    // Explicit open Interface transport only. Uses the slab retained owner;
    // pure reservoirs govern inflow, actual inside donor stocks govern outflow.
    // Temporal counterflow admits strict scratch-pure cells matching the
    // reservoir (no other-phase directional inflow), or producer-unresolved
    // double-positive stocks using the existing frozen-directional P0 law.
    // Both gross branches are retained. Resolved mixed chronology is unsupported.
    // No Hydro/material/KEP open-boundary capability is implied.
    static GeometricVofTransport prepareLocalPlicBoundarySweep(const QuadtreeMesh&,
        const GeometricVofPublication&,const GeometricVofAdvectorSlab&,
        long double liquidDensity,long double gasDensity,
        const GeometricVofTransportConfig& = {});
    GeometricVofIdentity sourceIdentity() const;
    GeometricVofTimeLayer sourceTime() const;
    GeometricVofTimeLayer endTime() const;
    const GeometricVofPublication& midpoint() const;
    const GeometricVofPublication& endpoint() const;
    const std::vector<GeometricVofTransportFace>& faces() const;
    const GeometricVofTransportReport& report() const;
    const GeometricVofPublication& source() const;
    const GeometricVofAdvectorSlab& advector() const;
    // Null for the original PLIC/fan paths; exact immutable owners for B1.
    const VofSmoothGeometry* smoothSource() const;
    const VofSmoothGeometry* smoothMidpoint() const;
    const VofSmoothGeometry* smoothEndpoint() const;
    long double liquidDensity() const;
    long double gasDensity() const;
    bool sameInstance(const GeometricVofTransport&) const noexcept;
    void validateFor(const QuadtreeMesh&,const GeometricVofPublication&,
        const GeometricVofAdvectorSlab&,long double liquidDensity,long double gasDensity) const;
    // LocalPlic: numerical P0 restriction of the SAME accepted split half/full
    // phase flux, using actual XY/RZ subface measure. Whole faces are cached.
    // CompatibleFan: same retained exact acceptance operator on old-fan subintervals.
    // Parent-normalized high phases define independent positive measures;
    // their sum defines the numerical total of a proper subinterval. Each
    // result retains original native incidence and both complete time extents.
    // Whole faces return cached accepted values exactly. Proper partitions
    // are additive within the original own-phase integration budget, not
    // necessarily bitwise after scalar conversion. No last-piece correction.
    std::vector<GeometricVofTransportFace> restrictNativeFaces(const QuadtreeMesh&,
        const std::vector<GeometricVofNativeSubface>&,
        GeometricVofTransportReport* diagnostics=nullptr) const;
    // LocalPlic: opposite actual side P0 accepted traces are integrated over
    // the query tangent interval (including hanging subfaces) and linearly
    // interpolated in the normal cell coordinate. RZ radial interpolation uses
    // the measure-compatible r^2 coordinate for r*u_r. This is a numerical
    // accepted-flux extension, not another phase transport or geometric stock.
    // Its begin/end velocity samples contain the atom normal component only.
    // CompatibleFan: same original AP source/slab with a fixed per-cell P1 fit of the
    // original native endpoint vectors. Evaluate at the physical centroid,
    // then freeze space along this atom and keep the full linear time law.
    // Interior acceptance is exactly unit: it never updates primal stock.
    // This does not replace native acceptance or construct a global C1 field.
    std::vector<GeometricVofInteriorFlux> measureInteriorFaces(const QuadtreeMesh&,
        const std::vector<GeometricVofInteriorSubface>&,
        GeometricVofTransportReport* diagnostics=nullptr) const;
    // LocalPlic only: the same actual-query validation and full arithmetic as
    // the paired APIs above. Unrequested half flux and endpoint samples are
    // not evaluated. Full limiter still includes the retained half limiter.
    std::vector<GeometricVofNativeFullFlux> restrictNativeFacesFull(const QuadtreeMesh&,
        const std::vector<GeometricVofNativeSubface>&,
        GeometricVofTransportReport* diagnostics=nullptr) const;
    std::vector<GeometricVofInteriorFullFlux> measureInteriorFacesFull(const QuadtreeMesh&,
        const std::vector<GeometricVofInteriorSubface>&,
        GeometricVofTransportReport* diagnostics=nullptr) const;
private:
    static GeometricVofLocalCourantRates localPlicCourantRatesImpl(const QuadtreeMesh&,
        const GeometricVofPublication&,const std::vector<double>&,const GeometricVofAdvectorSlab::Boundary*);
    static GeometricVofTransport prepareLocalPlicSweepImpl(const QuadtreeMesh&,
        const GeometricVofPublication&,const GeometricVofAdvectorSlab&,long double,long double,
        const GeometricVofTransportConfig&,bool open,
        const VofSmoothGeometry* = nullptr,const VofSmoothGeometryConfig* = nullptr,
        const GeometricVofPublication::MeshValidationBatch* = nullptr);
    struct Data;
    std::shared_ptr<const Data> data_;
    explicit GeometricVofTransport(std::shared_ptr<const Data>);
};
}
