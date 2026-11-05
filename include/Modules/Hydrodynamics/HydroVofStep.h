#pragma once
#include "Modules/Hydrodynamics/HydroFrozenMidpoint.h"
#include "Modules/Hydrodynamics/HydroTwoPhaseStep.h"
#include "Modules/Hydrodynamics/HydroVofAdvector.h"
#include <functional>
#include <cstdint>
#include <cstddef>

namespace LAMP::Hydrodynamics {
struct HydroVofAcConfig {
    double liquidDensity=1,gasDensity=1,liquidViscosity=0,gasViscosity=0,surfaceTension=0;
    HydroMomentumMass::DensityRepresentation momentumDensityRepresentation=HydroMomentumMass::DensityRepresentation::CellAverageProjection;
    Vector2D bodyAcceleration{0,0};
    HydroTwoPhaseTimeScheme scheme=HydroTwoPhaseTimeScheme::AdaptiveLinearMidpoint;
    double interfaceMismatchTolerance=.01,maximumCorrectionRelative=.25,correctionReferenceVelocity=1e-8;
    double capillaryCorrectionCourant=.1,maximumCapillaryCourant=.5,momentumPowerTolerance=1e-11;
    bool resetVelocityHistory=false;
    HydroSaddleConfig saddle;
    Interface::GeometricVofTransportConfig transport;
    // Explicit closed-flow transport preparation, called at EACH actual AC
    // stage. It must retain the passed accepted source and actual slab owner.
    // Empty preserves the existing PLIC path; open entry points reject it.
    using TransportProvider=std::function<Interface::GeometricVofTransport(const QuadtreeMesh&,
        const Interface::GeometricVofPublication&,const Interface::GeometricVofAdvectorSlab&,
        long double liquidDensity,long double gasDensity,const Interface::GeometricVofTransportConfig&)>;
    TransportProvider transportProvider;
    HydroCsfCurvatureConfig curvature;
    // App-owned physical preparation, called anew for EACH actual AC stage.
    // The returned typed load must retain this exact Midpoint/pressure owner.
    using MaxwellLoadProvider=std::function<HydroMaxwellLoad(const QuadtreeMesh&,
        const Interface::GeometricVofTransport&,const HydroPhaseGeometry&,const HydroPressureLoad&)>;
    MaxwellLoadProvider maxwellLoad;
    bool requireMaxwellCorrection=false;
    double maximumMaxwellLoadChange=.25,maximumMaxwellWorkChange=.25;
    // Called anew at EACH actual positive-sigma stage; the immutable result
    // must retain this exact Midpoint publication and complete clock.
    using CapillaryBoundaryProvider=std::function<HydroVofSouthEven90Boundary(const QuadtreeMesh&,
        const Interface::GeometricVofPublication&)>;
    CapillaryBoundaryProvider southEven90Boundary;

};
// Read-only stability scales of one accepted closed LocalPlic state.
struct HydroVofStabilityRates {
    Interface::GeometricVofLocalCourantRates transport;
    long double flowRate=0,capillaryRate=0; // Inverse code time.
    bool capillaryAvailable=false; // False for sigma=0 or no actual mixed leaf.
    double interfaceMinimumSpacing=0; // Zero when capillaryAvailable is false.
};
struct HydroVofStageReport {
    double interfaceMismatch=0,momentumPowerRelativeDefect=0;
    long double maximumAbsDualDefect=0,consistencyPower=0,kineticChange=0,capillaryWork=0,maxwellWork=0;
    // Actual outward F2 kinetic flux and original Saddle boundary work.
    std::array<long double,2> boundaryMomentumRate{};
    long double boundaryKineticWork=0,boundaryPhysicalKineticWork=0,inletDissipationWork=0;
    long double pressureOutletWork=0,prescribedReactionWork=0;
    // Signed F16 canonical/gross representation work; already included in
    // boundaryKineticWork, so this passive split is never added again.
    long double boundaryFluxRoundingWork=0;
    // Passive endpoint-KE ledger of this stage's actual M/F/R/u and loads.
    // halfRWork and its split are signed -1/2*u_mid^2*r on conservative VOF.
    // consistencyPower is the actually applied uR power (zero on this route).
    // Only the final accepted stage is a step work budget; AC stages are not summed.
    struct EnergyBudget {
        long double halfRWork=0,mappedPrimalHalfRWork=0,geometryTransferHalfRWork=0,spatialExchangeHalfRWork=0;
        long double timeMassWork=0,roundedPhaseMassWork=0,operatorRateRoundoffWork=0;
        long double rawMinusCanonicalOutflowHalfWork=0,halfRDecompositionResidual=0;
        long double bodyWork=0,viscousWork=0,pressureWork=0,momentumEquationWork=0;
        long double kineticBudgetResidual=0,closureResidual=0,closureRelativeResidual=0;
        // Stable pure gas: actual represented liquid volume is exactly zero
        // at begin/mid/end. Ratios use actual midpoint mass, not a mass floor.
        int stableGasDofs=0,peakGasNativeFace=-1;
        long double stableGasHalfRWork=0,stableGasSpatialExchangeHalfRWork=0;
        long double stableGasMaxAbsIntegratedDefectOverMidMass=0,peakGasAbsEndNormal=0;
        long double peakGasIntegratedDefectOverMidMass=0,peakGasMassChangeOverMidMass=0,peakGasOutwardIntegralOverMidMass=0;
        long double peakGasGeometryTransferOverMidMass=0,peakGasSpatialExchangeOverMidMass=0;
    } energy;
    bool maxwellLoadReportAvailable=false;
    HydroMaxwellLoadReport maxwellLoadReport;
    bool viscosityProductReportAvailable=false;
    HydroViscosityProductReport viscosityProduct;
    std::uint64_t viscosityProductPhaseIdentity=0;
    LAMP::TimeLayer viscosityProductTime;
    HydroSaddleReport saddle;
    // Two actual different-U open evaluations in this exact AC stage.
    // The solved midpoint reassessment remains part of diagnostics.
    struct OpenConvection {
        bool available=false;
        double boundaryBindingSeconds=0,operatorSeconds=0,advectorBoundarySeconds=0;
        double advectorApplySeconds=0,solvedMidpointApplySeconds=0;
        HydroKepOpenApplyReport advector,solvedMidpoint;
    } openConvection;
    // Passive producer receipt for this actual stage's optional rate rows.
    // rows is published only after every native rate is FP/errno-neutral;
    // scratchBytes is the conservative peak budget, including integer reserve.
    struct MomentumRateRows {
        bool requested=false,eligible=false,used=false,fallback=false;
        std::size_t rows=0,scratchBytes=0;
        int workers=0;
        double seconds=0;
    } momentumRateRows;
    struct Timing {
        double transport=0,material=0,capillary=0,maxwell=0,frozenNs=0,diagnostics=0,total=0;
        bool concurrentViscosityCapillary=false;
        // Passive subdivision of material; the original total remains authoritative.
        double phaseGeometry=0,materialTransport=0,dualMass=0,convection=0,viscosity=0;
    } timing;
};
struct HydroVofStepReport {
    int linearSolves=0;
    bool usedVelocityHistory=false,startupCorrection=false,corrected=false,capillaryForcedCorrection=false;
    double deltaTime=0,capillaryCourant=0,interfaceMismatch=0,finalInterfaceMismatch=0,correctionRelative=0;
    bool maxwellCorrectionChecked=false;
    double maxwellLoadChange=0,maxwellWorkChange=0;
    std::vector<HydroVofStageReport> stages;
};
// Pure alpha-only local PLIC orchestration over the same original native NS.
// No material-curve transport, surface-area inventory/inverse or extra solve.
class HydroVofStep {
public:
    // No next clock, solve, reconstruction, candidate or state mutation.
    // Actual accepted native normal means are mapped into VOF native rows.
    static HydroVofStabilityRates stabilityRates(const QuadtreeMesh&,
        const HydroTwoPhaseState&,const HydroVofAcConfig& = {});
    static HydroVofStep prepare(const QuadtreeMesh&,const HydroTwoPhaseState&,
        Interface::GeometricVofTimeLayer,const HydroVofAcConfig& = {});
    // Explicit actual F1/F2 steady-normal reservoir path, same original A/C
    // and native NS. Empty backflow model retains one-way pressure ports;
    // an explicit F16 model consumes actual gross/unknown-U momentum at each
    // stage. Free normal velocities remain original Saddle unknowns.
    static HydroVofStabilityRates stabilityRatesPrescribedFlow(const QuadtreeMesh&,
        const HydroTwoPhaseState&,const HydroVofAcConfig& = {});
    static HydroVofStep preparePrescribedFlow(const QuadtreeMesh&,const HydroTwoPhaseState&,
        Interface::GeometricVofTimeLayer,const HydroVofAcConfig& = {});
    const HydroTwoPhaseState& beginState() const noexcept;
    const HydroTwoPhaseState& candidateState() const noexcept;
    const HydroFrozenMidpoint& frozenMidpoint() const noexcept;
    const HydroNativeValues& solvedMidpointValues() const noexcept;
    const Interface::GeometricVofTransport& transport() const noexcept;
    bool hasCapillaryLoad() const noexcept;
    const HydroBalancedCsfLoad& capillaryLoad() const;
    const HydroVofStepReport& report() const noexcept;
    void commit(const QuadtreeMesh&,HydroTwoPhaseState&) const;
private:
    static HydroVofStabilityRates stabilityRatesImpl(const QuadtreeMesh&,
        const HydroTwoPhaseState&,const HydroVofAcConfig&,bool prescribed);
    static HydroVofStep prepareImpl(const QuadtreeMesh&,const HydroTwoPhaseState&,
        Interface::GeometricVofTimeLayer,const HydroVofAcConfig&,bool prescribed);
    struct Data;
    explicit HydroVofStep(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
};
}
