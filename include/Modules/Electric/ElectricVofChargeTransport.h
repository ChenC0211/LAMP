#pragma once
#include "Modules/Electric/ElectricVofField.h"
#include "Modules/Interface/GeometricVofTransport.h"
#include "Modules/Interface/GeometricVofMeshTransfer.h"
#include <memory>
#include <vector>

namespace LAMP::Electric
{
struct ElectricVofChargeFace
{
    int nativeRow=-1,negativeCell=-1,positiveCell=-1;
    long double totalVolume=0; // Accepted canonical VOF T, never phase flux.
    double integratedCharge=0; // Canonical signed Q_t; boundary donors follow both gross T branches, never Q sign.
    double lowOrderIntegratedCharge=0,highOrderIntegratedCharge=0,limiter=1;
};
struct ElectricVofChargeAdvectionReport
{
    bool boundaryEnabled=false;
    std::size_t boundaryReservoirFaces=0;
    // Signed charge carried along the physical flow, not positive magnitudes.
    double boundaryIncomingCharge=0,boundaryOutgoingCharge=0,boundaryNetCharge=0;
    double boundaryProductRoundingResidual=0,maximumBoundaryProductRoundingError=0;
    // Canonical shared Qt minus the exact sum of its represented outgoing and
    // incoming branches, with physical outward orientation; not product RNE.
    double boundaryNetRoundingResidual=0,maximumBoundaryNetRoundingError=0;
    double boundaryModelChargeResidual=0,relativeBoundaryModelChargeResidual=0;
    std::size_t cells=0,nativeFaces=0,donors=0,nonzeroFaceCharges=0,quantizedZeroPieces=0;
    double maximumOutflowCfl=0,sourceTotalCharge=0,candidateTotalCharge=0;
    double chargeResidual=0,relativeChargeResidual=0,roundingResidual=0,maximumCellRoundingError=0,totalSeconds=0;
    // Per actual half/full advance. No PDE, phase or charge owner is cached.
    bool donorParallelUsed=false,targetParallelUsed=false,parallelFallback=false;
    std::size_t parallelWorkers=0,parallelScratchBytes=0;
    bool limitedLinearFct=false;
    std::size_t reconstructedCells=0,highOrderFaces=0,limitedFaces=0;
    double minimumFaceLimiter=1,maximumDensityBoundDefect=0;
};
// Source/transport-bound constant-per-face signed reservoir charge density.
// Actual F1 gross volume selects boundary inflow/outflow, including temporal
// counterflow. The fluid reservoir phase never sets Q. Interior donor flow
// uses the selected shared-flux scheme; exterior gross branches keep P0.
// Half/full are independently advanced from the same accepted source.
class ElectricVofChargeBoundary
{
public:
    static ElectricVofChargeBoundary prepare(const QuadtreeMesh&,
        const ElectricVofCharge& acceptedSource,const Interface::GeometricVofTransport&,
        const std::vector<double>& inflowDensityByNativeFace);
    const ElectricVofCharge& sourceCharge() const;
    const Interface::GeometricVofTransport& transport() const;
    const std::vector<double>& inflowChargeDensities() const;
    bool allClosed() const;
    void validateFor(const QuadtreeMesh&,const ElectricVofCharge&,
        const Interface::GeometricVofTransport&) const;
    bool sameInstance(const ElectricVofChargeBoundary&) const noexcept;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
    explicit ElectricVofChargeBoundary(std::shared_ptr<const Data>);
    friend class ElectricVofAdvectedCharge;
    void validateForInBatch(const QuadtreeMesh&,const ElectricVofCharge&,
        const Interface::GeometricVofTransport&,Detail::ElectricVofValidationBatch&) const;
};
// Typed actual-stage receipt; there is no public vector/clock constructor.
class ElectricVofAdvectedCharge
{
public:
    const ElectricVofCharge& sourceCharge() const;
    const Interface::GeometricVofTransport& transport() const;
    const Interface::GeometricVofPublication& geometry() const;
    LAMP::TimeLayer timeLayer() const;
    const std::vector<double>& integratedCharge() const;
    const std::vector<ElectricVofChargeFace>& faces() const;
    const ElectricVofChargeAdvectionReport& report() const;
    const ElectricVofChargeBoundary* boundary() const; // nullptr on the original closed/delegate path.
    void validateFor(const QuadtreeMesh&,const ElectricVofCharge&,
        const Interface::GeometricVofTransport&,LAMP::TimeLayer) const;
    bool sameInstance(const ElectricVofAdvectedCharge&) const noexcept;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
    explicit ElectricVofAdvectedCharge(std::shared_ptr<const Data>);
    friend class ElectricVofChargeTransport;
    friend class ElectricVofField;
    void validateForInBatch(const QuadtreeMesh&,const ElectricVofCharge&,
        const Interface::GeometricVofTransport&,LAMP::TimeLayer,Detail::ElectricVofValidationBatch&) const;
};
class ElectricVofChargeTransport
{
public:
    static ElectricVofChargeTransport prepareBoundary(const QuadtreeMesh&,
        const ElectricVofCharge& acceptedSource,const Interface::GeometricVofTransport&,
        const ElectricVofChargeBoundary&,const ElectricVofControls& = {});
    static ElectricVofChargeTransport prepareClosed(const QuadtreeMesh&,
        const ElectricVofCharge& acceptedSource,const Interface::GeometricVofTransport&,
        const ElectricVofControls& = {});
    // Only the actual source->Midpoint extent is advanced. The receipt owns
    // its original controls/law so completion cannot substitute new inputs.
    static ElectricVofAdvectedCharge prepareMidpointClosed(const QuadtreeMesh&,
        const ElectricVofCharge& acceptedSource,const Interface::GeometricVofTransport&,
        const ElectricVofControls& = {});
    static ElectricVofAdvectedCharge prepareMidpointBoundary(const QuadtreeMesh&,
        const ElectricVofCharge& acceptedSource,const Interface::GeometricVofTransport&,
        const ElectricVofChargeBoundary&,const ElectricVofControls& = {});
    // Keep the SAME half receipt; independently advance accepted Q with the
    // actual full VOF flux. End input is rejected, never relabeled Midpoint.
    static ElectricVofChargeTransport completeEndpoint(const QuadtreeMesh&,const ElectricVofAdvectedCharge& midpoint);
    const ElectricVofCharge& sourceCharge() const;
    const Interface::GeometricVofTransport& transport() const;
    const ElectricVofAdvectedCharge& midpoint() const noexcept { return midpoint_; }
    const ElectricVofAdvectedCharge& endpoint() const noexcept { return endpoint_; }
private:
    ElectricVofAdvectedCharge midpoint_,endpoint_;
    ElectricVofChargeTransport(ElectricVofAdvectedCharge,ElectricVofAdvectedCharge);
};
struct ElectricVofChargeRemapReport
{
    std::size_t pieces=0,refinedDonors=0,unchangedTargets=0,coarsenedTargets=0,quantizedZeroPieces=0;
    std::size_t highOrderRefinedDonors=0,limitedRefinedDonors=0,finiteEligibilityP0Donors=0;
    std::size_t zeroNetQuantizedDonors=0;
    bool limitedLinearRemap=false;
    double sourceTotalCharge=0,candidateTotalCharge=0,chargeResidual=0,relativeChargeResidual=0;
    double roundingResidual=0,maximumCellRoundingError=0,totalSeconds=0;
    double prolongationRoundingResidual=0,restrictionRoundingResidual=0;
    double maximumDonorRelativeChargeResidual=0,minimumDonorLimiter=1;
    double maximumZeroNetQuantizationAdjustment=0;
    double maximumLinearQuantizationAdjustment=0,linearQuantizationL1Adjustment=0;
};
// Actual same-End AMR receipt. Old electric potential/face arrays are not
// published on the target layout; a new field solve is required there.
class ElectricVofChargeMeshTransfer
{
public:
    static ElectricVofChargeMeshTransfer prepare(const QuadtreeMesh& sourceMesh,
        const Interface::GeometricVofPublication& sourceGeometry,const ElectricVofCharge& source,
        const QuadtreeMesh& targetMesh,const Interface::GeometricVofMeshTransfer&,
        const ElectricVofControls& = {});
    const Interface::GeometricVofPublication& sourceGeometry() const;
    const ElectricVofCharge& sourceCharge() const;
    const Interface::GeometricVofMeshTransfer& interfaceTransfer() const;
    const ElectricVofCharge& chargeCandidate() const;
    const ElectricVofChargeRemapReport& report() const;
    void validateFor(const QuadtreeMesh& sourceMesh,
        const Interface::GeometricVofPublication& sourceGeometry,const ElectricVofCharge& source,
        const QuadtreeMesh& targetMesh,const Interface::GeometricVofMeshTransfer&) const;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
    explicit ElectricVofChargeMeshTransfer(std::shared_ptr<const Data>);
};
}
