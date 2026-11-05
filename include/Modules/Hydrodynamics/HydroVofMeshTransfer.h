#pragma once
#include "Modules/Hydrodynamics/HydroTwoPhaseStep.h"
#include "Modules/Interface/GeometricVofMeshTransfer.h"
#include <array>
#include <cstddef>
#include <memory>

namespace LAMP::Hydrodynamics
{
enum class HydroVofMeshVelocityTransfer { RegularFluxThenNativeMassProjection };
struct HydroVofMeshTransferConfig
{
    HydroVofMeshVelocityTransfer velocityTransfer=HydroVofMeshVelocityTransfer::RegularFluxThenNativeMassProjection;
    HydroMomentumMass::DensityRepresentation momentumDensityRepresentation=HydroMomentumMass::DensityRepresentation::SharpGeometry;
    HydroSaddleConfig saddle;
};
struct HydroVofRemapNativeBudget
{
    std::array<long double,2> mass{},momentum{};
    long double kinetic=0;
};
// Disjoint children of the original fluid-remap total; no State/mesh owner.
struct HydroVofRemapPreparationTimings
{
    bool available=false;
    double inputValidationSeconds=0,targetLayoutMetricSeconds=0,targetMomentumPhaseSeconds=0;
    double targetMassSeconds=0,sourceMassSeconds=0,postProjectionValidationSeconds=0;
    double targetStateSeconds=0,budgetValidationSeconds=0,unattributedSeconds=0;
};
struct HydroVofMeshTransferReport
{
    HydroVofMeshVelocityTransfer velocityTransfer=HydroVofMeshVelocityTransfer::RegularFluxThenNativeMassProjection;
    HydroMomentumMass::DensityRepresentation momentumDensityRepresentation=HydroMomentumMass::DensityRepresentation::SharpGeometry;
    std::size_t sourceLeafCandidates=0,faceSpans=0,polynomialIntegrals=0,wallNormals=0,axisNormals=0;
    int remapSolves=0;
    bool historyReset=false;
    bool prescribedFlow=false,pressureGaugeApplied=true;
    std::size_t fixedNormals=0,pressureOutletNormals=0,boundaryCoverageSpans=0;
    // Same-End projection constraint work, never physical timestep work.
    long double projectionBoundaryWork=0,projectionRhsRoundingWork=0;
    long double projectionFreeResidualWork=0,projectionContinuityWork=0;
    long double projectionEnergyIdentityScale=0,projectionEnergyIdentityRelativeResidual=0;
    double maximumRawRelativeDivergence=0,maximumFinalRelativeDivergence=0;
    long double rawDivergenceL1=0,finalDivergenceL1=0;
    long double pressureIntegralBeforeGauge=0,pressureGaugeShift=0,pressureGaugeRelativeResidual=0;
    HydroVofRemapNativeBudget source,rawTarget,target;
    std::array<long double,2> momentumChange{},projectionMomentumChange{};
    long double kineticChange=0,projectionKineticChange=0,projectionCorrectionEnergy=0,projectionEnergyIdentityResidual=0;
    double rawFluxSeconds=0,projectionSeconds=0,totalSeconds=0;
    HydroSaddleReport projectionSaddle;
    HydroVofRemapPreparationTimings preparation;
};
// Same accepted End clock, two immutable native meshes. Interface alone
// transfers both independent VOF stocks. Actual old regular-source face integrals seed one named
// native-mass constrained remap using the ORIGINAL saddle solver. Its multiplier
// is a remap diagnostic, never physical pressure or an AC time-stage solve.
// Physical p is overlap-volume averaged, gauged only without a physical
// pressure outlet; history is reset. Open projection uses lambda_outlet=0.
// Changes of diagonal native momentum/KE are reported, never external work.
// This law is M-orthogonal, not strict conservation of total native momentum.
// A no-op tree must bypass this changed-tree service.
class HydroVofMeshTransfer
{
public:
    static HydroVofMeshTransfer prepare(const QuadtreeMesh& sourceMesh,const HydroTwoPhaseState&,
        const QuadtreeMesh& targetMesh,Interface::GeometricVofMeshTransfer,double liquidDensity,double gasDensity,
        const HydroVofMeshTransferConfig& = {});
    // Source polynomial integrals in NATIVE-ROW order (interior entries 0).
    // Caller maps only fixed port entries into target Hydro layout g; free
    // pressure outlet payload remains 0. This service changes no flow/state.
    static std::vector<double> boundaryNormalMeans(const QuadtreeMesh& sourceMesh,
        const HydroTwoPhaseState& sourceEnd,const QuadtreeMesh& targetMesh,const NativeLeafDomain& targetDomain);
    static HydroVofMeshTransfer preparePrescribedFlow(const QuadtreeMesh& sourceMesh,const HydroTwoPhaseState&,
        const QuadtreeMesh& targetMesh,Interface::GeometricVofMeshTransfer,double liquidDensity,double gasDensity,
        HydroVofPrescribedFlow targetFlow,const HydroVofMeshTransferConfig& = {});
    // Actual stored RN(M*v) supplied to the original saddle; readonly receipt.
    const std::vector<double>& projectionMomentumRhs() const noexcept;
    const HydroTwoPhaseState& sourceState() const noexcept;
    const Interface::GeometricVofMeshTransfer& interfaceTransfer() const noexcept;
    const HydroNativeValues& rawValues() const noexcept;
    const std::vector<double>& projectionMultiplier() const noexcept;
    const HydroTwoPhaseState& candidateState() const noexcept;
    const HydroVofMeshTransferReport& report() const noexcept;
private:
    static HydroVofMeshTransfer prepareImpl(const QuadtreeMesh&,const HydroTwoPhaseState&,
        const QuadtreeMesh&,Interface::GeometricVofMeshTransfer,double,double,
        const HydroVofMeshTransferConfig&,const HydroVofPrescribedFlow*);
    struct Data;
    explicit HydroVofMeshTransfer(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
};
}
