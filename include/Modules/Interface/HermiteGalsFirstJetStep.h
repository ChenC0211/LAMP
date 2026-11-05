#pragma once

#include "Modules/Interface/HermiteGalsPublication.h"
#include "Modules/Interface/HermiteGalsTrace.h"
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Geometry/LevelSetMappedSegment.h"
#include <cstddef>
#include <functional>
#include <memory>
#include <vector>
#include <variant>

namespace LAMP::Interface
{
struct HermiteGalsBoundaryFirstJet { double phi = 0; Vector2D gradient; };
// Neutral consumer of an already resolved map and its SAME derivative.
// Owns its immutable callable and explicit provenance; no Hydro dependency.
class HermiteGalsResolvedMotion
{
public:
    using Query=std::function<LevelSetBoundaryTrace::JetResult(Vector2D)>;
    HermiteGalsResolvedMotion(HermiteGalsIdentity source,HermiteGalsTimeLayer sourceTime,
        HermiteGalsIdentity velocity,HermiteGalsTimeLayer targetTime,Query);
    HermiteGalsIdentity sourceIdentity() const noexcept {return source_;}
    HermiteGalsIdentity velocityIdentity() const noexcept {return velocity_;}
    HermiteGalsTimeLayer sourceTime() const noexcept {return sourceTime_;}
    HermiteGalsTimeLayer targetTime() const noexcept {return targetTime_;}
    LevelSetBoundaryTrace::JetResult trace(Vector2D) const;
private:
    HermiteGalsIdentity source_,velocity_;
    HermiteGalsTimeLayer sourceTime_,targetTime_;
    Query query_;
};
// Actual prescribed data at the OUTER arrival point/time, not a source-field
// extension or an interior reset. Captures must be immutable and repeatable.
using HermiteGalsBoundaryFirstJetQuery = std::function<HermiteGalsBoundaryFirstJet(
    Vector2D, HermiteGalsTimeLayer, Boundary)>;
// Generic expression/data initialization, using the same master/hanging Q4
// owner and FirstJet publication closure as transport. No IO or case identity.
using HermiteGalsInitialFirstJetQuery=std::function<HermiteGalsBoundaryFirstJet(Vector2D)>;
HermiteGalsPublication initializeHermiteGalsFirstJet(const QuadtreeMesh&,bool radial,
    const HermiteGalsInitialFirstJetQuery&,const HermiteGalsPublicationConfig& = {});
// SymmetryFirstJet explicitly imposes homogeneous scalar Neumann data on the
// declared outer config.volume.symmetryPlanes. The real pulled value and
// tangential gradient remain unchanged; this is not velocity parity.
enum class HermiteGalsBoundaryTreatment { InflowHistory, PrescribedOuterFirstJet, SymmetryFirstJet };
// An explicit temporal projection rule, independent of case identifiers.
// The compensated rule requires two interior queries of the same closed map.
enum class HermiteGalsTimeProjection { CharacteristicPullback, ClosedProjectionCompensated };
// Parallel explicitly permits concurrent queries of the same immutable resolved motion.
// Unknown external callables retain the default Serial execution.
enum class HermiteGalsPointExecution { Serial, Parallel };
enum class HermiteGalsPointSampling { Exact, Numerical };
struct HermiteGalsPointSamplingConfig
{
    HermiteGalsPointSampling mode=HermiteGalsPointSampling::Exact;
    double valueAbsoluteTolerance=0,gradientAbsoluteTolerance=0; // Explicit physical units for Numerical.
};
// Passive steady-clock diagnostics. Never inputs to numerical acceptance.
// The seven phases partition preparation; total also includes timer/report overhead.
struct HermiteGalsFirstJetStepTimings
{
    double input=0,boundary=0,vertices=0,edges=0,centers=0,mixed=0,publication=0,total=0;
};
struct HermiteGalsFirstJetStepReport
{
    std::size_t transportedVertices = 0, transportedEdges = 0, tracedCenters = 0;
    std::size_t prescribedVertices = 0, prescribedEdges = 0;
    std::size_t sourceCells = 0, sourceNativeFaces = 0, sourceCutCells = 0, sourceGraphs = 0;
    // Local kernel diagnostics only; the full producer leaves these at zero.
    std::size_t requestedNativeCells=0, stagedNativeCells=0, evaluatedNativeCells=0;
    // Actual snapshot point-query calls, including staged/second-map queries; not motion calls.
    std::size_t evaluatedPointQueries=0;
    HermiteGalsPointExecution pointExecution=HermiteGalsPointExecution::Serial;
    std::size_t pointJobs=0,pointTraceQueries=0; // Actual distinct-point jobs and motion.trace calls.
    int pointThreads=1; // Largest actual point-team size, never a requested thread count.
    double pointWorkerSeconds=0; // Sum of active worker elapsed intervals; excludes final barrier, not wall time.
    HermiteGalsPointSampling pointSampling=HermiteGalsPointSampling::Exact;
    std::size_t numericalPointQueries=0,numericalValueFallbacks=0,numericalGradientFallbacks=0;
    double maxAcceptedValueErrorBound=0,maxAcceptedGradientErrorBound=0; // Actual returned query enclosures.
    HermiteGalsFirstJetStepTimings timings;
};
class HermiteGalsFirstJetStep
{
public:
    HermiteGalsIdentity sourceIdentity() const;
    HermiteGalsIdentity arrivalVelocityIdentity() const;
    const HermiteGalsPublication& publication() const;
    const HermiteGalsFirstJetStepReport& report() const;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
    explicit HermiteGalsFirstJetStep(std::shared_ptr<const Data>);
    const Data& checked() const;
    friend HermiteGalsFirstJetStep prepareHermiteGalsFirstJetStep(const QuadtreeMesh&,
        const HermiteGalsPublication&, const LevelSetBoundaryTrace::Domain&,
        const HermiteGalsVelocityView&, const HermiteGalsVelocityView&,
        const HermiteGalsBoundaryFirstJetQuery&, const HermiteGalsPublicationConfig&,HermiteGalsPointExecution,const HermiteGalsPointSamplingConfig&);
    friend HermiteGalsFirstJetStep prepareHermiteGalsFirstJetStep(const QuadtreeMesh&,
        const HermiteGalsPublication&,const HermiteGalsResolvedMotion&,
        const HermiteGalsBoundaryFirstJetQuery&,const HermiteGalsPublicationConfig&,HermiteGalsBoundaryTreatment,HermiteGalsTimeProjection,
        HermiteGalsPointExecution,const HermiteGalsPointSamplingConfig&);
};
// A private temporal scalar sample for consumers that need only the actual
// transported Q4 and clock. This is not an I11 publication or an accepted
// state: no inventory, topology or surface-normal/curvature certificate is
// claimed. Copying its geometry view retains the same unique scalar owner.
class HermiteGalsFirstJetSample
{
public:
    HermiteGalsIdentity identity() const;
    HermiteGalsTimeLayer timeLayer() const;
    HermiteGalsMeasure measure() const;
    HermiteGalsGeometryView geometry() const;
    HermiteGalsIdentity sourceIdentity() const;
    HermiteGalsIdentity arrivalVelocityIdentity() const;
    const HermiteGalsFirstJetStepReport& report() const;
    void validateMesh(const QuadtreeMesh&) const;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
    explicit HermiteGalsFirstJetSample(std::shared_ptr<const Data>);
    const Data& checked() const;
    friend HermiteGalsFirstJetSample prepareHermiteGalsFirstJetSample(const QuadtreeMesh&,
        const HermiteGalsPublication&,const HermiteGalsResolvedMotion&,
        const HermiteGalsBoundaryFirstJetQuery&,const HermiteGalsPublicationConfig&,
        HermiteGalsBoundaryTreatment,HermiteGalsTimeProjection);
};
HermiteGalsFirstJetSample prepareHermiteGalsFirstJetSample(const QuadtreeMesh& mesh,
    const HermiteGalsPublication& source,const HermiteGalsResolvedMotion&,
    const HermiteGalsBoundaryFirstJetQuery& boundary = {},const HermiteGalsPublicationConfig& config = {},
    HermiteGalsBoundaryTreatment treatment=HermiteGalsBoundaryTreatment::InflowHistory,
    HermiteGalsTimeProjection projection=HermiteGalsTimeProjection::CharacteristicPullback);
// Restricted uniform-XY scalar support. Only nativeCells() have actual Q4
// data; no full-domain state/publication or value for missing support is exposed.
// Mixed five-point halos and compensated departure cells are evaluated from
// the same original source and resolved motion, not filled with zero data.
class HermiteGalsFirstJetLocalSample
{
public:
    HermiteGalsIdentity identity() const;
    HermiteGalsTimeLayer timeLayer() const;
    HermiteGalsMeasure measure() const;
    HermiteGalsGeometryView geometry() const;
    const std::vector<int>& nativeCells() const;
    HermiteGalsIdentity sourceIdentity() const;
    HermiteGalsIdentity arrivalVelocityIdentity() const;
    const HermiteGalsFirstJetStepReport& report() const;
    void validateMesh(const QuadtreeMesh&) const;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
    explicit HermiteGalsFirstJetLocalSample(std::shared_ptr<const Data>);
    const Data& checked() const;
    friend class HermiteGalsFirstJetLocalPreparation;
    friend HermiteGalsFirstJetLocalSample prepareHermiteGalsFirstJetLocalSample(const QuadtreeMesh&,
        const HermiteGalsPublication&,const HermiteGalsResolvedMotion&,const std::vector<int>&,
        const HermiteGalsBoundaryFirstJetQuery&,const HermiteGalsPublicationConfig&,
        HermiteGalsBoundaryTreatment,HermiteGalsTimeProjection,HermiteGalsPointExecution,const HermiteGalsPointSamplingConfig&);
};
bool supportsHermiteGalsFirstJetLocalSample(const QuadtreeMesh&,const HermiteGalsPublication&);
HermiteGalsFirstJetLocalSample prepareHermiteGalsFirstJetLocalSample(const QuadtreeMesh& mesh,
    const HermiteGalsPublication& source,const HermiteGalsResolvedMotion&,
    const std::vector<int>& requestedNativeCells,
    const HermiteGalsBoundaryFirstJetQuery& boundary = {},const HermiteGalsPublicationConfig& config = {},
    HermiteGalsBoundaryTreatment treatment=HermiteGalsBoundaryTreatment::InflowHistory,
    HermiteGalsTimeProjection projection=HermiteGalsTimeProjection::CharacteristicPullback,
    HermiteGalsPointExecution execution=HermiteGalsPointExecution::Serial,
    const HermiteGalsPointSamplingConfig& sampling={});
// Read-only time sampling of the original source for actual phase-rate traces.
// This owns cardinal data or an authentic complete local Q4 fallback. It has
// no new primary/state identity and exposes no geometry or publication view.
class HermiteGalsFirstJetTraceSample
{
public:
    using Trace=std::variant<LevelSetMappedSegment::NativeQuarticTrace,LevelSetMappedSegment::CardinalQuarticTrace>;
    struct Segment {Vector2D begin{0,0},end{0,0};int normalAxis=0;};
    HermiteGalsIdentity sourceIdentity() const;
    HermiteGalsIdentity arrivalVelocityIdentity() const;
    HermiteGalsTimeLayer sourceTimeLayer() const;
    HermiteGalsTimeLayer timeLayer() const;
    HermiteGalsMeasure measure() const;
    const std::vector<int>& nativeCells() const;
    std::size_t nativeCellCount() const;
    std::array<Vector2D,2> nativeCellBounds(int) const;
    const HermiteGalsFirstJetStepReport& report() const;
    Trace trace(std::size_t segment) const;
    bool ownsTrace(const Trace&) const;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
    explicit HermiteGalsFirstJetTraceSample(std::shared_ptr<const Data>);
    const Data& checked() const;
    friend class HermiteGalsFirstJetLocalPreparation;
};
// One actual accepted publication owns the frozen native topology used by a
// synchronous transport batch. Sampling reads no live mesh and accepts no
// boundary callback; the batch caller still validates its live mesh at its
// entrance and successful return.
class HermiteGalsFirstJetLocalPreparation
{
public:
    static HermiteGalsFirstJetLocalPreparation prepare(const QuadtreeMesh&,const HermiteGalsPublication&);
    HermiteGalsIdentity sourceIdentity() const;
    void validateMesh(const QuadtreeMesh&) const;
    HermiteGalsFirstJetLocalSample sample(const HermiteGalsResolvedMotion&,const std::vector<int>& requestedNativeCells,
        const HermiteGalsPublicationConfig& config={},
        HermiteGalsTimeProjection projection=HermiteGalsTimeProjection::CharacteristicPullback,
        HermiteGalsPointExecution execution=HermiteGalsPointExecution::Serial,
        const HermiteGalsPointSamplingConfig& sampling={},
        HermiteGalsBoundaryTreatment treatment=HermiteGalsBoundaryTreatment::InflowHistory) const;
    // Only actual cardinal-row datum points are transported. The raw staged
    // compensation source remains a complete Q4; unavailable cardinal/zero
    // traces retain the original complete local scalar/Taylor path.
    HermiteGalsFirstJetTraceSample sampleTraces(const HermiteGalsResolvedMotion&,
        const std::vector<int>& requestedNativeCells,const std::vector<HermiteGalsFirstJetTraceSample::Segment>&,
        const HermiteGalsPublicationConfig& config={},
        HermiteGalsTimeProjection projection=HermiteGalsTimeProjection::CharacteristicPullback,
        HermiteGalsPointExecution execution=HermiteGalsPointExecution::Serial,
        const HermiteGalsPointSamplingConfig& sampling={},
        HermiteGalsBoundaryTreatment treatment=HermiteGalsBoundaryTreatment::InflowHistory) const;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
    explicit HermiteGalsFirstJetLocalPreparation(std::shared_ptr<const Data>);
    const Data& checked() const;
    HermiteGalsState bindLocal(HermiteGalsTimeLayer,const HermiteGalsPrimaryState&,const HermiteGalsClosure&,
        const LevelSetHermiteGeometry::Snapshot::NativeLeafSubset&) const;
};
// Prepare one Midpoint/End candidate from an accepted Initial/End publication.
// Arrival velocity time is the explicit target clock; its geometry identity
// identifies the velocity provenance, not the newly generated publication.
// No accepted commit, regrid, inventory correction, output or fallback.
HermiteGalsFirstJetStep prepareHermiteGalsFirstJetStep(const QuadtreeMesh& mesh,
    const HermiteGalsPublication& source, const LevelSetBoundaryTrace::Domain& domain,
    const HermiteGalsVelocityView& previousVelocity, const HermiteGalsVelocityView& arrivalVelocity,
    const HermiteGalsBoundaryFirstJetQuery& boundary = {}, const HermiteGalsPublicationConfig& config = {},
    HermiteGalsPointExecution execution=HermiteGalsPointExecution::Serial,
    const HermiteGalsPointSamplingConfig& sampling={});
HermiteGalsFirstJetStep prepareHermiteGalsFirstJetStep(const QuadtreeMesh& mesh,
    const HermiteGalsPublication& source,const HermiteGalsResolvedMotion&,
    const HermiteGalsBoundaryFirstJetQuery& boundary = {},const HermiteGalsPublicationConfig& config = {},
    HermiteGalsBoundaryTreatment treatment=HermiteGalsBoundaryTreatment::InflowHistory,
    HermiteGalsTimeProjection projection=HermiteGalsTimeProjection::CharacteristicPullback,
    HermiteGalsPointExecution execution=HermiteGalsPointExecution::Serial,
    const HermiteGalsPointSamplingConfig& sampling={});
}
