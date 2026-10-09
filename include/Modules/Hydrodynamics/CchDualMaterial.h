#pragma once

#include "Modules/Hydrodynamics/CchMaterialTransport.h"
#include "Modules/Hydrodynamics/CchMomentumGeometry.h"
#include <array>
#include <memory>
#include <vector>

namespace Quantum::Hydrodynamics
{
struct CchDualMaterialCell
{
    long double massChange=0, outwardIntegral=0, integratedDefect=0, defectRate=0, relativeDefect=0;
    long double mappedPrimalMassChange=0, geometryTransferChange=0, mappedPrimalIntegratedDefect=0;
    long double outwardRate=0, mappedPrimalOutwardRate=0, spatialExchangeRate=0;
};
struct CchDualMaterialBalance
{
    std::array<long double,3> geometryMass{};
    long double massChange=0, boundaryIntegral=0, integratedDefect=0, l1IntegratedDefect=0;
    long double maxAbsDefectRate=0, maxRelativeCellDefect=0, balanceIdentity=0;
};
struct CchDualMaterialReport
{
    std::array<CchDualMaterialBalance,2> components{}; // x/r, y/z
    // Actual E_sp diagnostics, never required to vanish.
    long double maxRelativeSpatialExchange=0, maxAbsSpatialExchange=0, maxDefectDecomposition=0;
};
struct CchModalMaterialEnergyReport
{
    CchKineticHodgeEnergyChange kinetic;
    std::array<long double,2> legacyDiagonalMaterialChange{};
};
// Fixed-layout consumer of the UNIQUE H08 spatial mass flux and the SAME
// actual momentumGeometry. Retains three true H09a masses; continuity is diagnostic,
// never repaired. No velocity update, pressure solve, clock or commit.
class CchDualMaterial
{
public:
    static CchDualMaterial prepare(const QuadtreeMesh&, CchMaterialTransport,
        CchMomentumMass begin, CchMomentumMass midpoint, CchMomentumMass end);
    const CchMaterialTransport& material() const noexcept;
    const CchMomentumMass& beginMass() const noexcept;
    const CchMomentumMass& midpointMass() const noexcept;
    const CchMomentumMass& endMass() const noexcept;
    const CchMomentumGeometry& geometry() const noexcept;
    const CchDofLayout& layout() const noexcept;
    CchMaterialSources sources() const;
    Interface::HermiteGalsTimeLayer timeLayer() const;
    double deltaTime() const noexcept;
    const std::vector<long double>& faceRates() const noexcept; // Original H08 dual view, H09a segment order.
    const std::vector<long double>& faceIntegrals() const noexcept; // Original H08 time integral; never reintegrated here.
    const std::vector<CchDualMaterialCell>& cells() const noexcept; // velocity DOFs
    const CchDualMaterialReport& report() const noexcept;
    void validateFor(const QuadtreeMesh&, CchMaterialSources,
        Interface::HermiteGalsTimeLayer actualMidpointClock) const;
    // Derived integral momentum load u_mid*r_dual, in native velocity order.
    // Caller must retain this prepared owner and validate actual provenance;
    // a returned load array is not an accepted/candidate state.
    std::vector<long double> consistencyImpulse(const QuadtreeMesh&, const CchNativeValues&,
        CchMaterialSources actualSources, Interface::HermiteGalsTimeLayer actualClock) const;
    // Integral material kinetic term 1/2*u_mid^2*r_dual, separately by component.
    // This is not the complete H10/H14 energy balance or exact endpoint KE change.
    std::array<long double,2> materialKineticChange(const QuadtreeMesh&, const CchNativeValues&,
        CchMaterialSources actualSources, Interface::HermiteGalsTimeLayer actualClock) const;
    // Full [U,q] endpoint metric diagnostic bound to these EXACT three masses
    // and the unique midpoint F/R owner. Keeps legacy diagonal u^2 R separate;
    // this does not fabricate the still-required augmented consistency load.
    CchModalMaterialEnergyReport modalEnergyChange(const QuadtreeMesh&,
        const CchKineticHodge& begin,const CchKineticHodge& midpoint,const CchKineticHodge& end,
        const CchSpatialAdvector& before,const CchSpatialAdvector& after) const;
private:
    struct Data;
    explicit CchDualMaterial(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
};
}
