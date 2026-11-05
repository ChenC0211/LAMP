#pragma once

#include "Modules/Hydrodynamics/HydroMaterialTransport.h"
#include "Modules/Hydrodynamics/HydroMomentumGeometry.h"
#include <array>
#include <memory>
#include <vector>

namespace LAMP::Hydrodynamics
{
struct HydroDualMaterialCell
{
    long double massChange=0, outwardIntegral=0, integratedDefect=0, defectRate=0, relativeDefect=0;
    long double mappedPrimalMassChange=0, geometryTransferChange=0, mappedPrimalIntegratedDefect=0;
    // Canonical interval averages from the same retained integrals:
    // I_dual/dt, (T I_primal)/dt, (I_dual-T I_primal)/dt.
    long double outwardRate=0, mappedPrimalOutwardRate=0, spatialExchangeRate=0;
    // Original sums of individually averaged face rates, never a second F
    // or a replacement for the canonical defectRate/consistency impulse.
    long double rawOutwardRate=0, rawMappedPrimalOutwardRate=0, rawSpatialExchangeRate=0;
};
struct HydroDualMaterialBalance
{
    std::array<long double,3> geometryMass{};
    long double massChange=0, boundaryIntegral=0, integratedDefect=0, l1IntegratedDefect=0;
    long double maxAbsDefectRate=0, maxRelativeCellDefect=0, balanceIdentity=0;
};
struct HydroDualMaterialReport
{
    std::array<HydroDualMaterialBalance,2> components{}; // x/r, y/z
    // Endpoint-mass/gross-flux scaled continuity diagnostics, not divided by dt.
    // The conservative LocalPlic P0 route requires both <=256 binary64 eps.
    long double maxRelativeContinuityDefect=0, maxRelativeStoredMassContinuityDefect=0;
    long double maxRelativePrimalContinuityDefect=0, maxRelativeMassTransferDefect=0, maxRelativeFluxCommutatorDefect=0;
    // Actual E_sp diagnostics; reference representations need not commute.
    long double maxRelativeSpatialExchange=0, maxAbsSpatialExchange=0, maxDefectDecomposition=0;
};
struct HydroModalMaterialEnergyReport
{
    HydroKineticHodgeEnergyChange kinetic;
    std::array<long double,2> legacyDiagonalMaterialChange{};
};
// Fixed-layout consumer of the UNIQUE H08 spatial mass flux and the SAME
// actual momentumGeometry. Retains three true H09a masses; continuity is diagnostic,
// never repaired. No velocity update, pressure solve, clock or commit.
class HydroDualMaterial
{
public:
    static HydroDualMaterial prepare(const QuadtreeMesh&, HydroMaterialTransport,
        HydroMomentumMass begin, HydroMomentumMass midpoint, HydroMomentumMass end);
    const HydroMaterialTransport& material() const noexcept;
    const HydroMomentumMass& beginMass() const noexcept;
    const HydroMomentumMass& midpointMass() const noexcept;
    const HydroMomentumMass& endMass() const noexcept;
    const HydroMomentumGeometry& geometry() const noexcept;
    const HydroDofLayout& layout() const noexcept;
    HydroMaterialSources sources() const;
    LAMP::TimeLayer timeLayer() const;
    double deltaTime() const noexcept;
    const std::vector<long double>& faceRates() const noexcept; // Original H08 dual view, H09a segment order.
    const std::vector<long double>& faceIntegrals() const noexcept; // Original H08 time integral; never reintegrated here.
    const std::vector<HydroDualMaterialCell>& cells() const noexcept; // velocity DOFs
    const HydroDualMaterialReport& report() const noexcept;
    // Fixed model selection: LocalPlic + inventory P0, with canonical AND
    // actual stored endpoint mass continuity certified during prepare.
    // A failed certificate rejects; it never selects a different operator.
    bool usesConservativeVofMomentum() const noexcept;
    void validateFor(const QuadtreeMesh&, HydroMaterialSources,
        LAMP::TimeLayer actualMidpointClock) const;
    // Diagnostic integral u_mid*r_dual, in native velocity order. The
    // conservative VOF operator does not apply this as a momentum load.
    // Caller must retain this prepared owner and validate actual provenance;
    // a returned load array is not an accepted/candidate state.
    std::vector<long double> consistencyImpulse(const QuadtreeMesh&, const HydroNativeValues&,
        HydroMaterialSources actualSources, LAMP::TimeLayer actualClock) const;
    // Integral material kinetic term 1/2*u_mid^2*r_dual, separately by component.
    // This is not the complete H10/H14 energy balance or exact endpoint KE change.
    std::array<long double,2> materialKineticChange(const QuadtreeMesh&, const HydroNativeValues&,
        HydroMaterialSources actualSources, LAMP::TimeLayer actualClock) const;
    // Full [U,q] endpoint metric diagnostic bound to these EXACT three masses
    // and the unique midpoint F/R owner. Keeps legacy diagonal u^2 R separate;
    // this does not fabricate the still-required augmented consistency load.
    HydroModalMaterialEnergyReport modalEnergyChange(const QuadtreeMesh&,
        const HydroKineticHodge& begin,const HydroKineticHodge& midpoint,const HydroKineticHodge& end,
        const HydroSpatialAdvector& before,const HydroSpatialAdvector& after) const;
private:
    friend class HydroKepConvection;
    // Same-call consumer only, after the original public source and the first
    // independent values validation. No mesh receipt is stored in this owner.
    std::vector<long double> consistencyImpulseAfterValidation(const HydroNativeValues&) const;
    struct Data;
    explicit HydroDualMaterial(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
};
}
