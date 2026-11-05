#pragma once

#include "Modules/Hydrodynamics/HydroCharacteristic.h"
#include "Modules/Interface/HermiteGalsFirstJetStep.h"
#include "Modules/Interface/GeometricVofTransport.h"
#include <optional>
#include <array>
#include <chrono>
#include <functional>
#include <memory>
#include <vector>

namespace LAMP::Hydrodynamics
{
struct HydroMaterialSources
{
    Interface::HermiteGalsIdentity begin, midpoint, end;
    HydroAdvectorIdentity advector; // Actual FirstJet spatial owner.
    HydroPhaseGeometry::Model model=HydroPhaseGeometry::Model::FirstJet;
    std::array<Interface::GeometricVofIdentity,3> vofIds{};
    std::optional<Interface::GeometricVofTransport> vofStamp; // Same immutable accepted operator.
};
struct HydroMaterialCell
{
    long double massBegin=0, massMidpoint=0, massEnd=0;
    long double massChange=0, outwardIntegral=0, integratedDefect=0, defectRate=0;
    long double relativeDefect=0; // |defect|/(|massChange|+sum |face integral|).
};
struct HydroMaterialReport
{
    // True when the material report uses retained primary-alpha inventory.
    bool usesPrimaryAlphaInventory=false;
    std::array<long double,3> geometryMass{{0,0,0}};
    long double massChange=0, boundaryIntegral=0, integratedDefect=0, l1IntegratedDefect=0;
    long double maxAbsDefectRate=0, maxRelativeCellDefect=0, balanceIdentity=0;
};
// Select once from the declared boundary/model, never after a failed trace.
// The compensated rule is currently certified only for complete uniform XY
// base-native sources with a full homogeneous impermeable normal law.
enum class HydroFirstJetTimeScheme
{ BoundaryLaw, CharacteristicPullback, ClosedProjectionCompensated };
// All embedded rules retain the same per-atom phase error/work budget.
enum class HydroSpaceTimeIntegrationRule
{ EmbeddedGaussLegendre35, EmbeddedSimpsonBoole, EmbeddedSelectiveBoole };
// Actual completed rate work within one prepare. A local sample is one
// (time, missing-atom group) evaluation. Prepared native cells count genuine
// new Q4 final and compensated-stage channels; the retained begin is excluded.
// Cache hits count no work. Full rate samples include begin/midpoint/end.
struct HydroSpaceTimeSamplingCost
{
    bool atomLocal=false;
    std::size_t localSamples=0,fullSamples=0,localPreparedNativeCells=0,fullPreparedNativeCells=0;
    std::size_t localAtomRateEvaluations=0,fullAtomRateEvaluations=0,localScalarPointQueries=0,localPointJobs=0,localPointTraceQueries=0,localNumericalPointQueries=0,localExactValueFallbacks=0,localExactGradientFallbacks=0;
    int pointThreads=1;
    double maxAcceptedTransportValueError=0,maxAcceptedTransportGradientError=0;
    Interface::HermiteGalsIdentity sourceIdentity;
    HydroAdvectorIdentity previousAdvectorIdentity,arrivalAdvectorIdentity;
};
// Caller-thread diagnostics only. The error is the latest global estimate
// (zero before its first evaluation), not an acceptance status. During a
// pending split, panels counts the live partition before both children land.
struct HydroSpaceTimeProgress
{
    std::size_t geometrySamples=0,panels=0;
    int deepestPanel=0;
    long double relativePhaseErrorEstimate=0;
    double elapsedSeconds=0;
    LAMP::TimeLayer timeLayer; // Actual complete Midpoint interval.
    HydroSpaceTimeSamplingCost samplingCost;
    HydroSpaceTimeIntegrationRule integrationRule=HydroSpaceTimeIntegrationRule::EmbeddedGaussLegendre35;
};
struct HydroSpaceTimeConfig
{
    HydroCharacteristicConfig characteristic;
    double relativeTolerance=1e-7;
    int initialPanels=4,maxDepth=24,maxGeometrySamples=8192;
    HydroFirstJetTimeScheme firstJetTimeScheme=HydroFirstJetTimeScheme::BoundaryLaw;
    // With a callback, select outer arrival data on all sides or only actual inflow hits.
    // No callback always uses InflowHistory; interior boundary events still reject.
    Interface::HermiteGalsBoundaryTreatment firstJetBoundaryTreatment=
        Interface::HermiteGalsBoundaryTreatment::PrescribedOuterFirstJet;
    HydroSpaceTimeIntegrationRule integrationRule=HydroSpaceTimeIntegrationRule::EmbeddedGaussLegendre35;
    bool enablePointParallel=true; // Typed immutable/no-callback local samples; disables outer time teams.
    double transportValueAbsoluteTolerance=1e-13,transportGradientAbsoluteTolerance=1e-11;
    bool enableNumericalTransportSampling=true; // Separate bounded Q4 point value/gradient contract; exact fallback.
    bool enableAtomLocalTimeIntegration=true; // Preselect supported uniform XY; never a failure fallback.
    std::function<void(const HydroSpaceTimeProgress&)> progress; // Caller only; exceptions propagate.
};
// Passive elapsed diagnostics for ONE prepared transport. FirstJet and rate
// fields sum individual task elapsed times; internal tasks can overlap and
// these sums are not wall-clock partitions. FirstJet total includes its child
// stages. fullGeometrySeconds includes the existing initial geometry build
// plus midpoint/end builds; fullSamples counts only the latter two samples.
// Both rate fields include creation of their actual interpolated advector.
// fullRateSeconds also includes begin/end rates on the retained full geometry
// used by the endpoint support guard; these do not add geometry samples.
// Characteristic sums are nested in FirstJet stages, not additional work.
struct HydroSpaceTimeTimings
{
    Interface::HermiteGalsFirstJetStepTimings fullFirstJet{},internalFirstJet{};
    double fullGeometrySeconds=0,fullRateSeconds=0,internalRateSeconds=0,totalWallSeconds=0;
    double fullCharacteristicSeconds=0,internalCharacteristicSeconds=0;
    // Non-overlapping caller wall stages: input/preparation, initial geometry,
    // full midpoint/end jobs (including their overlap), adaptive quadrature,
    // final source/partition checks and flux, material assembly.
    std::array<double,6> workflowWallSeconds{};
    std::size_t fullSamples=0,internalSamples=0,fullInternalBatchTasks=0;
};
struct HydroSpaceTimeReport
{
    std::size_t geometrySamples=0,acceptedPanels=0,initialPlanPanels=0;
    std::size_t endpointSupportPanels=0; // Evaluated panels with an unseen-phase endpoint guard.
    bool reusedTimePlan=false;
    HydroSpaceTimeIntegrationRule integrationRule=HydroSpaceTimeIntegrationRule::EmbeddedGaussLegendre35;
    int deepestPanel=0;
    long double absolutePhaseErrorEstimate=0,maxRelativePhaseErrorEstimate=0;
    Interface::HermiteGalsTimeProjection firstJetTimeProjection=
        Interface::HermiteGalsTimeProjection::CharacteristicPullback;
    // Effective law actually consumed by every full/local FirstJet sample.
    Interface::HermiteGalsBoundaryTreatment firstJetBoundaryTreatment=
        Interface::HermiteGalsBoundaryTreatment::InflowHistory;
    HydroSpaceTimeTimings timings;
    HydroSpaceTimeSamplingCost samplingCost; // acceptedPanels counts nonempty interval/atom groups on this route.
};
// UNIQUE instantaneous density-weighted spatial flux. Valid at any declared
// Initial/Midpoint/End layer; no dt, inventory, mass update or accepted state.
class HydroMaterialFlux
{
public:
    static HydroMaterialFlux prepare(const QuadtreeMesh&, HydroPhaseRate,
        double liquidDensity, double gasDensity);
    const HydroPhaseRate& phaseRate() const noexcept;
    const HydroMomentumGeometry& momentumGeometry() const;
    const HydroDofLayout& layout() const noexcept;
    HydroAdvectorIdentity advectorIdentity() const;
    Interface::HermiteGalsIdentity geometryIdentity() const;
    LAMP::TimeLayer timeLayer() const;
    double liquidDensity() const noexcept;
    double gasDensity() const noexcept;
    const std::vector<long double>& atomRates() const noexcept;
    const std::vector<long double>& faceRates() const noexcept;
    const std::vector<long double>& dualFaceRates() const noexcept;
    // One-sided point view of this SAME density-weighted physical flux.
    // Phase is supplied by the retained publication/cubature.
    Vector2D physicalFlux(int pressure,Vector2D point,bool liquid) const;
    void validateFor(const QuadtreeMesh&, const HydroMomentumGeometry& expected,
        Interface::HermiteGalsIdentity, HydroAdvectorIdentity,
        LAMP::TimeLayer) const;
private:
    struct Data;
    explicit HydroMaterialFlux(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
};
// One fixed-layout interval: measured geometry masses + UNIQUE material
// integral + explicit nonconservative defect. No update/commit or repair.
// prepare keeps the original midpoint rule; prepareSpaceTime integrates the
// actual moving geometry and endpoint source over the whole interval.
class HydroMaterialTransport
{
public:
    static HydroMaterialTransport prepare(const QuadtreeMesh&, HydroPhaseGeometry begin,
        HydroPhaseRate midpoint, HydroPhaseGeometry end, double liquidDensity, double gasDensity);
    static HydroMaterialTransport prepare(const QuadtreeMesh&, HydroPhaseGeometry begin,
        HydroMaterialFlux midpoint, HydroPhaseGeometry end);
    // Full interval of the actual accepted VOF operator. Geometry-only atom
    // topology is shared by primal and dual views; no FirstJet rate is invented.
    static HydroMaterialTransport prepareGeometricVof(const QuadtreeMesh&, HydroMomentumGeometry,
        HydroPhaseGeometry begin, HydroPhaseGeometry midpoint, HydroPhaseGeometry end,
        Interface::GeometricVofTransport);
    // Same numerical unbound reconstruction/atom partition. This boundary is
    // only the validated lifetime owner of geometry-only preparation; all
    // stage phase integrals, masses and R remain freshly assembled.
    static HydroMaterialTransport prepareGeometricVof(const QuadtreeMesh&, HydroMomentumGeometry,
        HydroPhaseGeometry begin, HydroPhaseGeometry midpoint, HydroPhaseGeometry end,
        Interface::GeometricVofTransport,const HydroAdvectorBoundary& cacheOwner);
    HydroPhaseGeometry::Model model() const noexcept;
    const HydroPhaseGeometry& midpointGeometry() const noexcept;
    const Interface::GeometricVofTransport& vofTransport() const;
    const HydroAdvectorReconstruction& fluxReconstruction() const;
    // True spacetime integration along ONE original FirstJet/endpoint source.
    // Every sampled geometry uses that original state, not interpolated snapshots.
    // previousPlan reuses ONLY an authenticated same-source/interval partition.
    // All new geometry/rates/errors are recomputed; invalid plans reject.
    // The callback supplies outer ARRIVAL FirstJet data, not an exterior scalar
    // extension. Config chooses all outer points (legacy default) or only actual
    // inflow hits; tangent/outflow points then use the original scalar pullback.
    // An interior boundary hit remains unsupported and rejects explicitly.
    static HydroMaterialTransport prepareSpaceTime(const QuadtreeMesh&,
        Interface::HermiteGalsPublication original,HydroCharacteristicSource,
        double liquidDensity,double gasDensity,HydroSpaceTimeConfig,
        const Interface::HermiteGalsBoundaryFirstJetQuery& boundary={},
        const Interface::HermiteGalsPublicationConfig& publication={},
        const HydroMaterialTransport* previousPlan=nullptr);
    bool isSpaceTime() const noexcept;
    const HydroCharacteristicSource& characteristicSource() const;
    const HydroSpaceTimeReport& spaceTimeReport() const;
    const std::vector<HydroPhaseMeasure>& phaseAtomIntegrals() const;
    const std::vector<long double>& atomIntegrals() const noexcept;
    const std::vector<long double>& dualFaceIntegrals() const noexcept;
    const HydroPhaseGeometry& beginGeometry() const noexcept;
    const HydroPhaseRate& midpointRate() const; // FirstJet only; VOF rejects.
    const HydroPhaseGeometry& endGeometry() const noexcept;
    const HydroDofLayout& layout() const noexcept;
    HydroMaterialSources sources() const;
    LAMP::TimeLayer timeLayer() const; // Full Midpoint clock.
    double deltaTime() const noexcept;
    double liquidDensity() const noexcept;
    double gasDensity() const noexcept;
    bool isSpatial() const noexcept;
    // Original spatial geometry and single density-weighted atom owner.
    // These ports reject temporary legacy preparations.
    const HydroMomentumGeometry& momentumGeometry() const;
    const HydroMaterialFlux& midpointFlux() const;
    const std::vector<long double>& atomRates() const;
    const std::vector<long double>& dualFaceRates() const;
    const std::vector<long double>& faceRates() const noexcept; // interval average for spacetime; midpoint otherwise.
    const std::vector<long double>& faceIntegrals() const noexcept;
    const std::vector<HydroMaterialCell>& cells() const noexcept; // pressure DOF order.
    const HydroMaterialReport& report() const noexcept;
    void validateFor(const QuadtreeMesh&, HydroMaterialSources,
        LAMP::TimeLayer midpointClock) const;
    void validateFor(const QuadtreeMesh&, const HydroMomentumGeometry& expected,
        HydroMaterialSources, LAMP::TimeLayer midpointClock) const;
private:
    static HydroMaterialTransport prepareGeometricVofBound(const QuadtreeMesh&,HydroMomentumGeometry,
        HydroPhaseGeometry begin,HydroPhaseGeometry midpoint,HydroPhaseGeometry end,
        Interface::GeometricVofTransport,const HydroAdvectorBoundary* cacheOwner);
    struct TimePanel { double lo,hi;int depth;std::vector<int> atomIds; }; // Sealed intervals/IDs only; empty IDs are the full route.
    static HydroMaterialTransport assemble(const QuadtreeMesh&, HydroPhaseGeometry,
        HydroPhaseRate, HydroPhaseGeometry, double, double, const HydroMaterialFlux*,
        const std::vector<HydroPhaseMeasure>* integrated=nullptr,const HydroCharacteristicSource* motion=nullptr,
        const HydroSpaceTimeReport* report=nullptr,
        const std::chrono::steady_clock::time_point* profilingStart=nullptr,
        const std::vector<TimePanel>* acceptedTimePanels=nullptr);
    struct Data;
    static HydroMaterialTransport assembleIntegrated(std::shared_ptr<Data>,
        std::vector<Interface::GeometricVofIntegratedFlux>);
    explicit HydroMaterialTransport(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
};
}

