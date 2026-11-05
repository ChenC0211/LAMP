#pragma once

#include "Modules/Hydrodynamics/HydroKepConvection.h"
#include "Modules/Hydrodynamics/HydroBalancedCsfLoad.h"
#include "Modules/Hydrodynamics/HydroInterfaceNormalLoad.h"
#include "Modules/Hydrodynamics/HydroPhaseViscosity.h"
#include "Modules/Hydrodynamics/HydroSinglePhaseSaddle.h"
#include <array>
#include <functional>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace LAMP::Hydrodynamics
{
// Steady native face-mean inlet/prescribed outlet and typed pressure outlets.
// All port lists are disjoint. Pressure outlets are free; their normalVelocity
// payload is exactly zero. All unlisted boundary normals are impermeable.
// Callables supply immutable, repeatable physical inlet data at the given clock.
// Inlet FirstJet is ARRIVAL data (not a hit-time/exterior extension).
struct HydroTwoPhasePrescribedFlow
{
    std::vector<int> inletVelocityDofs,outletVelocityDofs;
    std::vector<HydroSaddlePressureOutlet> pressureOutlets;
    std::vector<double> normalVelocity;
    using InletVelocityQuery=std::function<Vector2D(Vector2D,Interface::HermiteGalsTimeLayer,Boundary)>;
    InletVelocityQuery inletVelocity;
    Interface::HermiteGalsBoundaryFirstJetQuery inletFirstJet;
};
// A value-only steady native-row query. Construction seals the real layout,
// coverage and clock role; copies share immutable values, never user callbacks.
class HydroVofSteadyVelocity
{
public:
    enum class ClockRole { Endpoint, FullMidpoint };
    static HydroVofSteadyVelocity prepare(const QuadtreeMesh&,const HydroDofLayout&,
        std::vector<Vector2D> byNativeRow,std::vector<bool> coveredByNativeRow,ClockRole);
    Vector2D operator()(int nativeRow,LAMP::TimeLayer) const;
    bool matches(const HydroDofLayout&,ClockRole) const noexcept;
    void validateFor(const QuadtreeMesh&,const HydroDofLayout&,ClockRole) const;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
    explicit HydroVofSteadyVelocity(std::shared_ptr<const Data>);
    HydroVofSteadyVelocity reboundFor(const QuadtreeMesh&,const HydroDofLayout&,ClockRole) const;
    friend class HydroTwoPhaseState;
};
// Separate VOF physical boundary specification, with no FirstJet or exterior
// geometry owner. By default native normal means are STEADY; complete exterior
// vectors at accepted and next-End clocks provide fixed-inlet F2 components.
// Explicit pressure-driven reservoirs instead retain free pressure normals.
struct HydroVofPrescribedFlow
{
    std::vector<Interface::GeometricVofAdvectorSlab::Boundary::Mode> rulesByNativeRow;
    std::vector<int> prescribedVelocityDofs;
    std::vector<double> normalVelocity; // Hydro layout order; other entries 0.
    std::vector<HydroSaddlePressureOutlet> pressureOutlets;
    using EndpointVelocityQuery=std::function<Vector2D(int nativeRow,LAMP::TimeLayer endpoint)>;
    EndpointVelocityQuery endpointVelocity; // Repeatable declared face-P0 vector.
    // Opt-in steady face-P0 weak incoming momentum on actual free pressure
    // ports. Frozen at the full actual stage midpoint; never a normal constraint.
    // This generic ambient model is not interiorContinuationZeroShear.
    using PressureOutletBackflowMomentumQuery=std::function<Vector2D(int nativeRow,LAMP::TimeLayer fullMid)>;
    PressureOutletBackflowMomentumQuery pressureOutletBackflowMomentum;
    // Empty retains the original one-way model. The continuation token covers
    // momentum only; actual zero-shear viscous traction is a separate contract.
    std::optional<HydroVofPressureOutletMomentumModel> pressureOutletBackflowModel;
    // Opt-in pressure supply: every open reservoir is an actual free pressure
    // port, with at least two distinct physical tractions. Requires an explicit
    // backflow momentum model, no prescribedVelocityDofs or endpointVelocity,
    // and an exactly zero normalVelocity payload. Initial U=0 is valid; actual
    // later net/gross transport determines inflow and outflow. The existing
    // Interface pure-phase reservoir rules remain the phase donor authority.
    bool pressureDrivenReservoirs=false;
    // Explicit stationary full-trace no-slip for positive physical viscosity.
    // Sorted unique Hydro rows, a subset of the genuine positive-area Closed
    // walls; pressure/fixed ports, interiors and the RZ axis are excluded.
    // Empty preserves the original normal-impermeable strain boundary model.
    std::vector<int> stationaryNoSlipVelocityDofs;
};
// One actual accepted predecessor value, with separate model-specific source
// identities. No state-handle chain or fabricated FirstJet ID is retained.
struct HydroTwoPhaseVelocityHistory
{
    HydroNativeValues values;
    Interface::HermiteGalsTimeLayer clock;
    Interface::HermiteGalsIdentity source;
    HydroPhaseGeometry::Model model=HydroPhaseGeometry::Model::FirstJet;
    Interface::GeometricVofIdentity vofSource;
};
// FirstJet retains its complete uniform XY/RZ scope. VOF accepted state
// additionally supports a complete native AMR domain through createVof.
// Closed, balanced prescribed flow, or inlet with outward total-traction outlets.
// RZ uses true revolved measures and eliminates the zero-area radial axis velocity.
// All arrays and prepared geometry are retained through immutable handles.
class HydroTwoPhaseState
{
public:
    static HydroTwoPhaseState create(const QuadtreeMesh&, HydroPhaseGeometry,
        const HydroNativeValues&);
    static HydroTwoPhaseState createVof(const QuadtreeMesh&,HydroPhaseGeometry,const HydroNativeValues&);
    // Same accepted End on a changed native tree. Reuse the genuine target
    // M/geometry and explicitly discard incompatible velocity history.
    static HydroTwoPhaseState createVofRemapped(const QuadtreeMesh&, HydroPhaseGeometry,
        const HydroNativeValues&, HydroMomentumMass retainedMass);
    static HydroTwoPhaseState createVofAccepted(const QuadtreeMesh&,HydroPhaseGeometry,
        const HydroNativeValues&,const HydroTwoPhaseState& previous,
        std::optional<HydroMomentumMass> retainedMass=std::nullopt);
    static HydroTwoPhaseState createPrescribedFlow(const QuadtreeMesh&, HydroPhaseGeometry,
        const HydroNativeValues&, HydroTwoPhasePrescribedFlow);
    const HydroTwoPhasePrescribedFlow* prescribedFlow() const noexcept;
    // Explicit fixed-native-mesh steady-normal VOF flow; walls are exactly
    // the Closed F1 rows. Legacy pressure outlets are strictly outward unless a
    // backflow model is explicit; pressureDrivenReservoirs permits resting ports.
    static HydroTwoPhaseState createVofPrescribedFlow(const QuadtreeMesh&,HydroPhaseGeometry,
        const HydroNativeValues&,HydroVofPrescribedFlow);
    // Changed native mesh, same accepted End: preserve actual target M and
    // physical ports, with no incompatible old-layout velocity history.
    static HydroTwoPhaseState createVofPrescribedFlowRemapped(const QuadtreeMesh&,HydroPhaseGeometry,
        const HydroNativeValues&,HydroMomentumMass retainedMass,HydroVofPrescribedFlow);
    const HydroVofPrescribedFlow* vofPrescribedFlow() const noexcept;
    // Present only for fully value-owned flow; bound to this actual publication
    // and complete Initial/End clock, never an earlier state's source.
    const Interface::GeometricVofAdvectorSlab::Boundary* vofTransportBoundary() const noexcept;
    const HydroPhaseGeometry& phaseGeometry() const noexcept;
    HydroPhaseGeometry::Model model() const noexcept;
    const Interface::HermiteGalsPublication& publication() const;
    const HydroTwoPhaseVelocityHistory* velocityHistory() const noexcept;
    const HydroMomentumMass* retainedMass() const noexcept;
    const HydroNativeValues& values() const noexcept;
    // Passive samples of the retained native reconstruction in pressure order.
    // Same actual geometry/velocity owner; no field fit or invented identity.
    std::vector<Vector2D> sampleCellCenterVelocities(const QuadtreeMesh&) const;
    const HydroMomentumGeometry& momentumGeometry() const noexcept;
    const HydroAdvectorReconstruction& reconstruction() const noexcept;
    const HydroAdvectorBoundary& boundary() const noexcept;
    const HydroCharacteristicGeometry& characteristicGeometry() const;
    Interface::HermiteGalsIdentity geometryIdentity() const;
    Interface::HermiteGalsTimeLayer timeLayer() const;
    bool sameIdentity(const HydroTwoPhaseState&) const noexcept;
    void validateFor(const QuadtreeMesh&) const;
    void validateVofFor(const QuadtreeMesh&,Interface::GeometricVofIdentity,Interface::GeometricVofTimeLayer) const;
    // Versioned native checkpoint, issued only by a validated VOF owner.
    // Existing velocity history requires its genuine still-live predecessor;
    // a restored state retains that one predecessor for a subsequent capture.
    // These methods perform no file IO, initialization or mesh transfer.
    std::string checkpoint(const QuadtreeMesh&,const HydroTwoPhaseState* previous=nullptr) const;
    static HydroTwoPhaseState restoreCheckpoint(const QuadtreeMesh&,const std::string&,
        std::optional<HydroVofPrescribedFlow> configuredFlow=std::nullopt);
private:
    void validateForInBatch(const QuadtreeMesh&,Detail::HydroStateMeshValidationBatch&) const;
    static void validateVofPhaseInBatch(const QuadtreeMesh&,const HydroPhaseGeometry&,
        Interface::GeometricVofIdentity,LAMP::TimeLayer,Detail::HydroStateMeshValidationBatch&);
    static void supportedVofGeometryInBatch(const QuadtreeMesh&,const HydroPhaseGeometry&,
        Detail::HydroStateMeshValidationBatch&);
    static void acceptedValuesInBatch(const QuadtreeMesh&,const HydroPhaseGeometry&,
        const HydroNativeValues&,const std::vector<int>& walls,bool pressureAnchored,
        Detail::HydroStateMeshValidationBatch&);
    struct Data;
    explicit HydroTwoPhaseState(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
    friend class HydroTwoPhaseStep;
};

enum class HydroTwoPhaseTimeScheme
{ FullyCoupledMidpoint, ExtrapolatedMidpoint, AdaptiveLinearMidpoint };
// Explicit force candidates share the actual midpoint geometry/native momentum.
enum class HydroCapillaryModel { BalancedCsf, SharpNormalStress };
enum class HydroTwoPhaseRetryReason { Geometry, Correction, VelocityBlock, Capillary, Transport };
// Only recoverable numerical failures invite the App to reduce its own dt.
class HydroTwoPhaseStepRetry : public std::runtime_error
{
public:
    HydroTwoPhaseStepRetry(HydroTwoPhaseRetryReason reason,const std::string& message)
        :std::runtime_error(message),reason_(reason){}
    HydroTwoPhaseRetryReason reason() const noexcept{return reason_;}
private:
    HydroTwoPhaseRetryReason reason_;
};
struct HydroTwoPhaseIterationReport;
struct HydroTwoPhaseStepConfig
{
    double liquidDensity=1, gasDensity=1;
    HydroMomentumMass::DensityRepresentation momentumDensityRepresentation=HydroMomentumMass::DensityRepresentation::SharpGeometry;
    double liquidViscosity=0, gasViscosity=0;
    double surfaceTension=0; // Finite physical sigma >= 0; zero skips only a pure capillary load.
    // Synchronous App/Coupling wiring: receives this actual stage transport and
    // returns an authenticated coeval vector view. No electric solver dependency
    // or arbitrary point-force callback enters Hydro. Each A/C stage calls anew.
    std::function<Interface::HermiteGalsSurfaceVectorView(
        const QuadtreeMesh&,const HydroMaterialTransport&)> interfaceLoad;
    HydroCapillaryModel capillaryModel=HydroCapillaryModel::BalancedCsf;
    // BalancedCsf only: curvature uses adjacent interface cells.
    HydroCsfCurvatureConfig csfCurvature=[] {HydroCsfCurvatureConfig c;
        c.relativeQuadratureTolerance=1e-7;return c;}();
    // Physical constant acceleration (XY components; RZ permits only axial y).
    Vector2D bodyAcceleration{0,0}; // Exact zero retains the original momentum arithmetic.
    // Adopted AC route; legacy Picard remains an explicit scheme selection.
    HydroTwoPhaseTimeScheme timeScheme=HydroTwoPhaseTimeScheme::AdaptiveLinearMidpoint;
    double interfaceMismatchTolerance=.01; // Kinematic trigger, not an accuracy certificate.
    double maximumCorrectionRelative=.25,correctionReferenceVelocity=1e-8;
    // Adaptive AC only: a continuous-wave budget, not a discrete-spectrum certificate.
    double capillaryCorrectionCourant=.1,maximumCapillaryCourant=.5;
    bool resetVelocityHistory=false; // Caller reports a topology/reset event; create() also starts without history.
    int maximumIterations=12; // Explicit Picard budget 1..12; AC uses at most two stages.
    // Explicit Picard update targets; AC retains these gaps as diagnostics.
    double velocityUpdateTolerance=1e-11,geometryUpdateTolerance=1e-11;
    double momentumPowerTolerance=1e-11; // Actual retained-owner momentum power residual gate.
    // Reuse actual begin/midpoint/end and quarter rates; only atoms whose
    // lower-order indicator is over budget receive four extra eighth nodes.
    // Independent cumulative phase and endpoint-support gates remain.
    // Use one characteristic map consistently at every material time node.
    // Explicit spaceTime selection retains the compensated projection.
    HydroSpaceTimeConfig spaceTime=[] {HydroSpaceTimeConfig c;
        c.firstJetTimeScheme=HydroFirstJetTimeScheme::CharacteristicPullback;
        c.integrationRule=HydroSpaceTimeIntegrationRule::EmbeddedSelectiveBoole;return c;}();
    Interface::HermiteGalsPublicationConfig publication;
    std::optional<HydroConvectionDissipationConfig> dissipation;
    // Caller-thread passive diagnostics; an observer failure rejects this preparation.
    std::function<void(const HydroTwoPhaseIterationReport&)> onIteration;
};

struct HydroTwoPhaseIterationReport
{
    int iteration=0;
    HydroTwoPhaseTimeScheme timeScheme=HydroTwoPhaseTimeScheme::AdaptiveLinearMidpoint;
    double interfaceMismatch=0,linearMomentumDefect=0;
    bool geometryCompared=false;
    double velocityUpdateRelative=0, geometryUpdateRelative=0;
    double midpointGeometryUpdateRelative=0, endGeometryUpdateRelative=0;
    double frozenOwnerMomentumDefect=0; // Actual reassessment on this iteration retained F/M family.
    double momentumPowerRelativeDefect=0;
    HydroSaddleReport saddle;
    // Passive wall timings: endpoint, transport, dual mass, convection, capillary, viscosity, saddle, reassessment.
    std::array<double,8> stageSeconds{};
    double csfColorSeconds=0,csfProjectionSeconds=0; // Passive split of the retained capillary stage.
    std::size_t geometrySamples=0,acceptedPanels=0,initialPlanPanels=0,endpointSupportPanels=0;
    int deepestPanel=0;
    bool reusedTimePlan=false;
    Interface::HermiteGalsTimeLayer timeLayer;
    HydroSpaceTimeSamplingCost samplingCost; // Actual completed transport; no progress observer required.
};

struct HydroTwoPhaseStepReport
{
    int iterations=0; // Actual one/two AC stages, or explicitly selected Picard iterations.
    HydroTwoPhaseTimeScheme timeScheme=HydroTwoPhaseTimeScheme::AdaptiveLinearMidpoint;
    bool usedVelocityHistory=false,startupCorrection=false,corrected=false;
    double capillaryCourant=0; // Zero when adaptive capillary control is inactive.
    bool capillaryForcedCorrection=false;
    double interfaceMismatch=0,finalInterfaceMismatch=0,correctionRelative=0,linearMomentumDefect=0;
    double deltaTime=0, velocityUpdateRelative=0, geometryUpdateRelative=0;
    double advectorGuessSolveRelative=0;
    // The final nonlinear reassessment keeps the retained F/M/geometry family.
    double frozenOwnerMomentumDefect=0;
    long double kineticBegin=0, kineticEnd=0, kineticChange=0;
    long double inertialImpulseWork=0, temporalMassCrossTerm=0, temporalIdentityDefect=0;
    long double inertialPower=0, pressurePower=0, viscousPower=0;
    long double convectionPower=0, consistencyPower=0, numericalDissipationPower=0;
    long double sourceNumericalDissipationPower=0, wallReactionPower=0, momentumPowerDefect=0;
    // Actual prescribed normal reaction work; a free pressure outlet is excluded.
    long double pressureOutletPower=0; // Solved midpoint U dot original Saddle T, added once.
    // wallReactionPower includes actual work of all prescribed normal ports.
    std::array<long double,2> boundaryMomentumRate{};
    long double boundaryKineticPower=0,boundaryPhysicalKineticPower=0,inletDissipationPower=0;
    double momentumPowerRelativeDefect=0;
    double surfaceTension=0;
    HydroMomentumMass::DensityRepresentation momentumDensityRepresentation=HydroMomentumMass::DensityRepresentation::SharpGeometry;
    HydroCapillaryModel capillaryModel=HydroCapillaryModel::BalancedCsf;
    long double interfacePower=0; // Complete retained RHS work, including external traction.
    bool hasExternalInterfaceLoad=false,capillaryPowerAvailable=true;
    // Mixed-load separation is unavailable: these two diagnostics are NaN,
    // rather than labeling electric work as capillary work. Pure path is unchanged.
    long double capillaryPower=0;
    long double bodyPower=0; // Sum solved-midpoint U * retained midpoint-mass body load; external input work.
    long double interfaceMeasureBegin=0, interfaceMeasureMidpoint=0, interfaceMeasureEnd=0;
    long double surfaceEnergyBegin=0, surfaceEnergyEnd=0, surfaceEnergyChange=0;
    // Observed Delta(sigma*area)+dt*capillaryPower; not an acceptance certificate.
    long double capillaryEnergyDefect=0;
    std::vector<HydroTwoPhaseIterationReport> history;
};

// Adopted AC one/two-stage preparation, or explicitly selected legacy Picard.
// Uses the caller's explicit next End clock.
// Throws on any failure without modifying the accepted handle or its arrays.
// Success retains the actual characteristic/FirstJet/F/three-M/R/KEP/K/saddle
// family and the distinct endpoint guess and solve output. No IO or dt choice.
class HydroTwoPhaseStep
{
public:
    using Matrix=HydroSaddleVelocityBlock::Matrix;
    static HydroTwoPhaseStep prepare(const QuadtreeMesh&, const HydroTwoPhaseState& accepted,
        Interface::HermiteGalsTimeLayer nextEndClock, const HydroTwoPhaseStepConfig&);
    const HydroTwoPhaseState& beginState() const noexcept;
    const HydroTwoPhaseState& candidateState() const noexcept;
    const HydroNativeValues& endpointGuess() const noexcept;
    const HydroNativeValues& solvedMidpointValues() const noexcept;
    const HydroCharacteristicSource& characteristicSource() const noexcept;
    const HydroMaterialTransport& materialTransport() const noexcept;
    const HydroDualMaterial& dualMaterial() const noexcept;
    const HydroKepConvection& convection() const noexcept;
    const HydroKepConvectionResult& convectionResult() const noexcept;
    const HydroKepConvectionResult& reassessedConvectionResult() const noexcept;
    bool hasCapillaryLoad() const noexcept;
    HydroCapillaryModel capillaryModel() const noexcept;
    // Physical RHS for either selected model; sharp balance loads are negated once.
    const std::vector<double>& capillaryMomentumRhsLoad() const; // Historical getter: actual selected total RHS.
    const HydroBalancedCsfLoad& capillaryLoad() const; // BalancedCsf only; other model rejects.
    const HydroInterfaceNormalLoad& sharpNormalStressLoad() const;
    const HydroPhaseViscosity& viscosity() const noexcept;
    const HydroSinglePhaseSaddle& saddle() const noexcept;
    const HydroSaddleResult& solveResult() const noexcept;
    const Matrix& velocityMatrix() const noexcept;
    // Bulk RHS excludes the original saddle's separately retained outlet T.
    const std::vector<double>& integralMomentumRhs() const noexcept;
    // Actual native load M_mid[u] * bodyAcceleration[normalAxis(u)], including fixed rows.
    const std::vector<long double>& bodyLoads() const noexcept;
    const HydroTwoPhaseStepReport& report() const noexcept;
    // Mesh and expected begin handle are checked before one pointer assignment.
    // An obsolete candidate cannot overwrite a newer accepted state.
    void commit(const QuadtreeMesh&, HydroTwoPhaseState& accepted) const;
private:
    struct Data;
    explicit HydroTwoPhaseStep(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
};
}

