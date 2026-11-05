#pragma once

#include "Modules/Hydrodynamics/HydroDualMaterial.h"
#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <functional>
#include <utility>
#include <vector>

namespace LAMP::Hydrodynamics
{
struct HydroConvectionWeight { int velocity=-1; double weight=0; };
// Prepared final-operator correction: J_ij=coefficient, J_ji=-coefficient.
// Not another material face flux; both indices belong to one local patch.
struct HydroConvectionCoupling { int first=-1, second=-1; long double coefficient=0; };
// Explicit exterior segment relation: second translated into first's image.
// Geometry/port data only; no extra velocity unknown or material flux.
struct HydroConvectionPeriodicPair
{ int firstSegment=-1, secondSegment=-1; Vector2D translation{}; };
struct HydroConvectionMomentReport
{
    // iterations counts typed GI add/drop pivots; other paths retain their original solver iterations.
    int interiorRows=0, constraints=0, iterations=0;
    // Uniform typed closed path: actual physical F/2 pairs, no moment QP.
    int physicalPairing=0, physicalSegments=0, physicalPairs=0;
    long double constantError=0, affineError=0, quadraticError=0;
    int boundedRows=0;
    // Typed impermeable coefficientRatio bounds FREE complete rows only.
    long double coefficientRatio=0, quadraticBudgetRatio=0, maximumSupportRatio=0;
    long double maxReactionConstantRawError=0; // Unnormalized diagnostic; fixed C uses actual coefficient arithmetic.
    long double maxReactionCoefficientL1=0; // Actual fixed-normal row norm; finite diagnostic, not a stability bound.
    int projectionCalls=0;
    long double quadraticBoundError=0;
    // Explicit impermeable scope: only zero-gross FREE rows require zero K.
    int fixedRows=0, freeRows=0, zeroGrossRows=0, zeroGrossReactionRows=0;
    int meanP1Rows=0; // Typed free rows using the actual native face-mean functional.
    int meanP2Rows=0; // The same preparation also integrates xx/xy/yy advective responses.
    long double maxQuadraticTargetRatio=0; // Typed |P2 center| / original positive slab halfwidth.
    long double maxMeanPointDifference=0; // Absolute vector-component difference, in velocity units.
    long double maxNativeMeanDefect=0; // Normal defect divided by positive sample gross plus |native U|.
    int physicalCertificateAfterWorkLimit=0; // Historical diagnostic; typed GI rejects its pivot limit and leaves this zero.
    // New typed scope: raw sum(face rates) minus canonical sum(integrals)/dt.
    // Retained diagnostics only; no material R or flux is corrected.
    long double rateSumMinusIntegralRateSigned=0, rateSumMinusIntegralRateL1=0;
    long double maxAbsRateSumMinusIntegralRate=0;
};
struct HydroKepComponentReport
{
    long double momentumRate=0, boundaryMomentumRate=0, convectivePower=0;
    // materialKineticPower is the signed residual -1/2*u^2*R on certified
    // conservative VOF, +1/2*u^2*R on the retained corrected reference route.
    long double boundaryKineticPower=0, materialKineticPower=0, kineticRate=0;
    // Prescribed-flow boundary work: physical outward energy plus inlet upwind dissipation.
    long double boundaryPhysicalKineticPower=0, inletDissipationPower=0;
    // Explicit bidirectional boundary only: signed canonical/gross representation work.
    long double boundaryFluxRoundingPower=0, maximumBoundaryCanonicalGrossRateDefect=0;
    // Actual frozen exterior diagonal coefficient contracted with current U squared.
    long double boundaryDiagonalPower=0;
    long double momentumIdentity=0, kineticIdentity=0;
    long double dissipationPower=0, maximumSensor=0;
    int dissipationConnections=0, activeDissipationConnections=0;
};
// Only the six boundary contractions consumed by the production step.
// Absence is explicit: a loads-only request returns no budget.
struct HydroKepOpenBoundaryBudget
{
    struct Component
    {
        long double boundaryMomentumRate=0,boundaryKineticPower=0;
        long double boundaryPhysicalKineticPower=0,inletDissipationPower=0;
        long double boundaryFluxRoundingPower=0,boundaryDiagonalPower=0;
    };
    std::array<Component,2> components{};
};
struct HydroKepConvectionReport { std::array<HydroKepComponentReport,2> components{}; };
// Passive receipt of one ACTUAL open apply; not a cached/frozen-U evaluation.
// Children are disjoint intervals within totalSeconds; checks match the explicit request.
struct HydroKepOpenApplyReport
{
    bool available=false,parallelRequested=false,parallelEligible=false,parallelUsed=false,parallelFallback=false;
    unsigned contributionWorkers=0,gatherWorkers=0;
    std::size_t nativeRows=0,segments=0,pairs=0,scratchBytes=0;
    std::array<std::uint64_t,3> vofSources{};
    LAMP::TimeLayer timeLayer;
    double validationSeconds=0,preparationSeconds=0,segmentSeconds=0,diagonalSeconds=0;
    double pairSeconds=0,rowSeconds=0,ledgerSeconds=0,fallbackSeconds=0,totalSeconds=0;
};

// Explicit local stabilization. Density contrast alone never activates it.
// The default build overloads remain the unchanged center-only operator.
struct HydroConvectionDissipationConfig
{
    double strength=1, activation=0.2, saturation=0.5;
};
struct HydroConvectionDissipationSample
{
    int first=-1, second=-1, firstSegment=-1, secondSegment=-1;
    long double sensor=0, chi=0, halfJump=0, momentumRate=0, power=0;
};

enum class HydroVofPressureOutletMomentumModel { Ambient,InteriorContinuationZeroShear };

// Same actual F1 native subface gross, in outward physical orientation.
// Rates are passive rounded diagnostics; all derived momentum/work contractions
// consume the independent phase integrals and original material densities/dt.
struct HydroVofBoundaryMomentumFlux
{
    int nativeRow=-1,velocityDof=-1;
    HydroVofPressureOutletMomentumModel incomingModel=HydroVofPressureOutletMomentumModel::Ambient;
    std::array<long double,2> outgoingPhaseIntegral{},incomingPhaseIntegral{};
    long double canonicalOutwardMassRate=0,outgoingMassRate=0,incomingMassRate=0;
    double incomingMomentum=0;
};

// Explicit transported-component traces on external dual segments. Internal
// entries must be zero. Bound to the actual H09 owner, identities and clock;
// no default boundary law, hidden ghosts or accepted state.
class HydroConvectionBoundary
{
public:
    static HydroConvectionBoundary bind(const QuadtreeMesh&, HydroDualMaterial,
        const std::vector<double>& segmentTraces, HydroMaterialSources,
        LAMP::TimeLayer);
    // Actual F1 inflow velocity data, not a normal-only inferred vector.
    // Unused internal/zero/outflow entries are zero. Incoming endpoint normal
    // components match the retained source slab; tangential values are explicit.
    static HydroConvectionBoundary bindVofEndpointVelocity(const QuadtreeMesh&,HydroDualMaterial,
        const std::vector<Vector2D>& beginVelocity,const std::vector<Vector2D>& endVelocity,
        HydroMaterialSources,LAMP::TimeLayer fullMidpointClock);
    // Opt-in weak incoming-momentum model. Begin/end arrays are dual-segment
    // order for non-free prescribed inflow only, and must be steady. Free-port
    // arrays are zero; their finite ambient vector is native-row ordered and
    // is not a prescribed normal velocity. All unused ambient rows are zero.
    // InteriorContinuationZeroShear uses the CURRENT unknown component U as
    // weak P0 incoming momentum, not a frozen value; ambient rows must be zero.
    // This momentum token does not certify a viscous zero-shear discretization.
    static HydroConvectionBoundary bindVofOpenMomentum(const QuadtreeMesh&,HydroDualMaterial,
        const std::vector<Vector2D>& prescribedBeginBySegment,
        const std::vector<Vector2D>& prescribedEndBySegment,
        const std::vector<int>& freePressureOutletVelocityDofs,
        const std::vector<Vector2D>& ambientByNativeRow,
        HydroMaterialSources,LAMP::TimeLayer fullMidpointClock,
        HydroVofPressureOutletMomentumModel=HydroVofPressureOutletMomentumModel::Ambient);
    // Dedicated production producer: one full accepted-flux query batch also
    // certifies the original finer free-port atom view. Endpoint callbacks run
    // only for actual nonfree incoming native rows, source then End once per row.
    static HydroConvectionBoundary bindVofOpenProductionMomentum(const QuadtreeMesh&,HydroDualMaterial,
        const std::vector<int>& freePressureOutletVelocityDofs,
        const std::vector<Vector2D>& ambientByNativeRow,
        const std::function<Vector2D(int,LAMP::TimeLayer)>& endpointVelocity,
        HydroMaterialSources,LAMP::TimeLayer fullMidpointClock,
        HydroVofPressureOutletMomentumModel=HydroVofPressureOutletMomentumModel::Ambient);
    bool hasVofOpenMomentum() const noexcept;
    const std::vector<HydroVofBoundaryMomentumFlux>& vofOpenMomentumFluxes() const;
    const std::vector<int>& freePressureOutletVelocityDofs() const;
    bool hasVofEndpointVelocity() const noexcept;
    const std::vector<Vector2D>& beginVelocity() const;
    const std::vector<Vector2D>& endVelocity() const;
    const HydroDualMaterial& material() const noexcept;
    const std::vector<double>& traces() const noexcept;
    void validateFor(const QuadtreeMesh&, const HydroDualMaterial&, HydroMaterialSources,
        LAMP::TimeLayer) const;
private:
    friend class HydroKepConvection;
    struct Data;
    explicit HydroConvectionBoundary(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
};

class HydroKepConvectionResult
{
public:
    const HydroDualMaterial& material() const noexcept;
    const HydroConvectionBoundary& boundary() const noexcept;
    // Raw F*t^T*u diagnostic, NOT the final compatible native scatter.
    // Explicit omission is unavailable and returns an empty vector.
    bool hasCenterMomentumRates() const noexcept;
    const std::vector<long double>& centerMomentumRates() const noexcept;
    // Positive outward convection N, in native velocity order.
    const std::vector<long double>& momentumRates() const noexcept;
    // Certified LocalPlic inventory-P0 VOF: -N, with no uR source.
    // Retained reference representations: -N + u_mid*R_dual.
    // Excludes pressure/viscosity/external loads.
    const std::vector<long double>& momentumLoads() const noexcept;
    // Actual native stabilization L, including periodic mass redistribution.
    const std::vector<long double>& dissipationRates() const noexcept;
    const std::vector<HydroConvectionDissipationSample>& dissipationSamples() const noexcept;
    // The retained MomentumLoads request defines momentum/eta and known-zero
    // stabilization; its nine power fields are unavailable quiet NaNs.
    bool hasPowerReport() const noexcept;
    // New production requests omit component diagnostics: all thirteen
    // non-stabilization fields are quiet NaNs; known-zero stabilization stays zero.
    bool hasComponentReport() const noexcept;
    const HydroKepOpenBoundaryBudget* openBoundaryBudget() const noexcept;
    const HydroKepConvectionReport& report() const noexcept;
    const HydroKepOpenApplyReport& openApplyReport() const noexcept;
private:
    friend class HydroKepConvection;
    struct Data;
    explicit HydroKepConvectionResult(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
};

// Compatible native transport: physical two-point typed closed pairing;
// historical opaque/periodic moment construction retains its declared scope.
// Consumes the unique H09 flux; optional explicit local nonnegative dissipation.
// No mass/velocity update, IO or clock ownership.
class HydroKepConvection
{
public:
    // Raw phase-volume Galerkin A_ij = integral B_i dot(F dot grad B_j).
    // Uses the unique flux's physical view and the same full [U,q] source.
    // No forced skew projection or R/boundary cancellation. General temporal,
    // boundary and continuity pairing remains the caller's obligation.
    static HydroKineticHodge::Matrix modalBulkTransportMatrix(const QuadtreeMesh&,
        const HydroKineticHodge&,const HydroMaterialFlux&);
    static HydroKepConvection build(const QuadtreeMesh&, HydroDualMaterial);
    static HydroKepConvection build(const QuadtreeMesh&, HydroDualMaterial,
        const std::vector<HydroConvectionPeriodicPair>&);
    static HydroKepConvection buildControlled(const QuadtreeMesh&, HydroDualMaterial,
        HydroConvectionDissipationConfig);
    static HydroKepConvection buildControlled(const QuadtreeMesh&, HydroDualMaterial,
        const std::vector<HydroConvectionPeriodicPair>&, HydroConvectionDissipationConfig);
    // Complete uniform FirstJet or actual 2:1 AMR VOF XY/RZ domain, genuine
    // typed impermeable source and zero
    // exterior F. Physical neighboring pairs consume the unique F/2; canonical
    // constant/skew and the original free-row L1 bound are certified.
    // Fixed normal rows carry wall reaction. P1/P2 moment-QP selection is
    // replaced by this physical operator; VOF authenticates its native slab
    // directly, with actual space accuracy checked
    // independently; no periodic aliases or opaque-trace fallback.
    static HydroKepConvection buildImpermeable(const QuadtreeMesh&, HydroDualMaterial,
        const HydroAdvectorBoundary&);
    static HydroKepConvection buildImpermeableControlled(const QuadtreeMesh&, HydroDualMaterial,
        const HydroAdvectorBoundary&, HydroConvectionDissipationConfig);
    // Complete FirstJet uniform XY/RZ source with typed
    // zero-normal walls and explicit
    // inlet transported-component traces. Outflow uses the native owner value;
    // inlet data and unique F are retained at the same full midpoint clock.
    // The caller certifies atom flow direction throughout the full interval.
    static HydroKepConvection buildPrescribedFlow(const QuadtreeMesh&, HydroDualMaterial,
        const HydroAdvectorBoundary& walls,const HydroConvectionBoundary& prescribedInlet);
    // Explicit physical pressure-outlet normals remain unknown in momentum.
    // All other exterior normals retain the four-argument prescribed-flow law.
    // Free ports must carry no incoming unique material flux; their actual
    // source atom integrals must retain the original native normal mean.
    static HydroKepConvection buildPrescribedFlow(const QuadtreeMesh&, HydroDualMaterial,
        const HydroAdvectorBoundary& walls,const HydroConvectionBoundary& prescribedInlet,
        const std::vector<int>& freeOutflowVelocityDofs);
    // Explicit LocalPlic F1 boundary source. Full-vector inlet/backflow
    // data must come from bindVofEndpointVelocity; pressure outlets are
    // one-way and all other outer normals remain prescribed. Reuses the
    // same physical pairing/upwind boundary-work kernel as FirstJet.
    static HydroKepConvection buildVofPrescribedFlow(const QuadtreeMesh&,HydroDualMaterial,
        const HydroAdvectorBoundary& closedWalls,const HydroConvectionBoundary& endpointInlet,
        const std::vector<int>& freeOutflowVelocityDofs = {});
    // Actual bidirectional gross and explicit steady ambient momentum. The
    // free native normal remains an unknown; this factory never imposes b.n.
    static HydroKepConvection buildVofOpenFlow(const QuadtreeMesh&,HydroDualMaterial,
        const HydroAdvectorBoundary& closedWalls,const HydroConvectionBoundary& openMomentum);
    // Production preparation omits points()/weights() center-diagnostic
    // tables (both empty). Explicit open Full apply still evaluates its center
    // diagnostic directly from the genuine two-point segment incidence.
    // The three unrequested raw-rate/integral difference diagnostics are quiet NaNs.
    static HydroKepConvection buildVofOpenProductionFlow(const QuadtreeMesh&,HydroDualMaterial,
        const HydroAdvectorBoundary& closedWalls,const HydroConvectionBoundary& openMomentum);
    bool isVofOpenFlow() const noexcept;
    // Complete result, including the raw center diagnostic.
    HydroKepConvectionResult applyVofOpen(const QuadtreeMesh&,const HydroNativeValues&,
        HydroMaterialSources,LAMP::TimeLayer fullMidpointClock) const;
    // Same actual loads, component reports and identity gates; omits only the
    // unconsumed raw center diagnostic and its diagnostic-only failures/flags.
    HydroKepConvectionResult applyVofOpenWithoutCenterDiagnostic(const QuadtreeMesh&,const HydroNativeValues&,
        HydroMaterialSources,LAMP::TimeLayer fullMidpointClock) const;
    // Predictor demand: exact rates/loads, momentum balance and eta only.
    // Omits center and energy-only checks/errors; solved-U power is a separate request.
    HydroKepConvectionResult applyVofOpenMomentumLoads(const QuadtreeMesh&,const HydroNativeValues&,
        HydroMaterialSources,LAMP::TimeLayer fullMidpointClock) const;
    // Production row contraction. Retains actual rates/loads and source/boundary
    // checks; omits diagnostic component identities and their errors/FP flags.
    HydroKepConvectionResult applyVofOpenProductionLoads(const QuadtreeMesh&,const HydroNativeValues&,
        HydroMaterialSources,LAMP::TimeLayer fullMidpointClock) const;
    // Same row contraction, plus the six actual exterior-only budget fields.
    HydroKepConvectionResult applyVofOpenBoundaryBudget(const QuadtreeMesh&,const HydroNativeValues&,
        HydroMaterialSources,LAMP::TimeLayer fullMidpointClock) const;
    bool isImpermeable() const noexcept;
    bool isPrescribedFlow() const noexcept;
    const std::vector<int>& fixedVelocityDofs() const noexcept;
    const HydroDualMaterial& material() const noexcept;
    const HydroDofLayout& layout() const noexcept;
    const std::vector<Vector2D>& points() const noexcept;
    const std::vector<std::vector<HydroConvectionWeight>>& weights() const noexcept;
    const std::vector<HydroConvectionCoupling>& couplings() const noexcept;
    const HydroConvectionMomentReport& momentReport() const noexcept;
    const std::vector<HydroConvectionPeriodicPair>& periodicPairs() const noexcept;
    const std::vector<int>& representatives() const noexcept;
    // Frozen typed physical, sensor-free operator. At the same owner/clock,
    // apply gives momentumLoads = -H*v+B. Closed B is exactly zero; prescribed
    // flow retains the inlet affine load and the outward upwind diagonal.
    // Uses final physical pairs, without probing or a second material flux.
    HydroKineticHodge::Matrix frozenMomentumOperator(const QuadtreeMesh&,
        HydroMaterialSources actualSources,
        LAMP::TimeLayer actualMidpointClock) const;
    std::vector<long double> frozenBoundaryMomentumLoad(const QuadtreeMesh&,
        HydroMaterialSources actualSources,
        LAMP::TimeLayer actualMidpointClock) const;
    HydroKepConvectionResult apply(const QuadtreeMesh&, const HydroNativeValues&,
        const HydroConvectionBoundary&, HydroMaterialSources actualSources,
        LAMP::TimeLayer actualMidpointClock) const;
    HydroKepConvectionResult applyPaired(const QuadtreeMesh&, const HydroNativeValues&,
        HydroMaterialSources actualSources, LAMP::TimeLayer actualMidpointClock) const;
private:
    friend class HydroFrozenMidpoint;
    // Only Frozen calls this after validating the same actual material,
    // convection and boundary in its callback-free preparation chain.
    std::pair<HydroKineticHodge::Matrix,std::vector<long double>> prepareFrozenMomentum(
        const QuadtreeMesh&,const HydroDualMaterial& validatedMaterial,
        HydroMaterialSources,LAMP::TimeLayer) const;
    HydroKineticHodge::Matrix frozenMomentumOperatorAfterValidation() const;
    struct Data;
    HydroKepConvectionResult applyVofOpenImpl(const QuadtreeMesh&,const HydroNativeValues&,
        HydroMaterialSources,LAMP::TimeLayer,bool includeCenterDiagnostic,bool includePowerReport) const;
    HydroKepConvectionResult applyVofOpenProductionImpl(const QuadtreeMesh&,const HydroNativeValues&,
        HydroMaterialSources,LAMP::TimeLayer,bool includeBoundaryBudget) const;
    explicit HydroKepConvection(std::shared_ptr<const Data>);
    static HydroKepConvection buildVofOpenImpl(const QuadtreeMesh&,HydroDualMaterial,
        const HydroAdvectorBoundary&,const HydroConvectionBoundary&,bool includeCenterTables);
    static HydroKepConvection buildPrescribedFlowImpl(const QuadtreeMesh&,HydroDualMaterial,
        const HydroAdvectorBoundary&,const HydroConvectionBoundary&,const std::vector<int>&,bool vof);
    static HydroKepConvection buildWithPorts(const QuadtreeMesh&, HydroDualMaterial,
        const std::vector<HydroConvectionPeriodicPair>*, const HydroConvectionDissipationConfig* = nullptr,
        const HydroAdvectorBoundary* = nullptr);
    std::shared_ptr<const Data> data_;
};
}
