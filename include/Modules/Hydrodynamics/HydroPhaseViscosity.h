#pragma once

#include "Modules/Hydrodynamics/HydroEnergyViscosity.h"
#include "Modules/Hydrodynamics/HydroPhaseGeometry.h"
#include "Modules/Interface/HermiteGalsPhaseQuadrature.h"
#include <cstdint>
#include <memory>

namespace LAMP::Hydrodynamics
{
class HydroKineticHodge;
class HydroMomentumGeometry;
// Retained phase source -> existing native P2 strain -> nonnegative mu.
// Hermite uses phase samples; VOF explicitly selects cell-average material
// coefficients. Neither creates a second strain or material inventory.
class HydroPhaseViscosity
{
public:
    enum class MaterialModel { HermitePhaseQuadrature, VofCellAverage=2, ExactZero };
    // Same physical source and complete phase-volume partition, including U-q
    // cross terms. This bulk form does not supply the finite-mu jump trace.
    static HydroEnergyViscosity::Matrix modalBulkMatrix(const QuadtreeMesh&,const HydroKineticHodge&,
        double liquidViscosity,double gasViscosity);
    // Exact zero phase coefficients need no volume/strain sampling. The same
    // phase geometry and full clock remain authoritative for every use.
    static HydroPhaseViscosity zero(const QuadtreeMesh&, HydroPhaseGeometry);
    static HydroPhaseViscosity build(const QuadtreeMesh&, HydroPhaseGeometry,
        Interface::HermiteGalsPhaseQuadrature,double liquidViscosity,double gasViscosity);
    // Genuine LocalPlic independent L/G stocks divided by actual cell volume.
    // No rounded alpha/complement, recut geometry or Hermite sampling source.
    static HydroPhaseViscosity buildVofCellAverage(const QuadtreeMesh&, HydroPhaseGeometry,
        double liquidViscosity,double gasViscosity);
    // Reuse only the canonical raw strain of a genuine previous VofCellAverage
    // material on the same actual mesh/layout/D owner. Current mu and K are fresh.
    static HydroPhaseViscosity buildVofCellAverage(const QuadtreeMesh&, HydroPhaseGeometry,
        double liquidViscosity,double gasViscosity,const HydroPhaseViscosity& previous);
    // Consume only the genuine raw view of the builtin complete stress factory
    // on this actual layout/D. Current phase stocks, mu and K remain fresh.
    static HydroPhaseViscosity buildVofCellAverage(const QuadtreeMesh&, HydroPhaseGeometry,
        double liquidViscosity,double gasViscosity,const HydroStrainOperator& canonicalStressGeometry);
    // Explicit stationary wall raw strain. Only identical wall rows and the
    // same actual mesh/layout/D owner permit previous geometry reuse; current
    // inventories, phase coefficients and K are always prepared afresh.
    static HydroPhaseViscosity buildVofCellAverage(const QuadtreeMesh&, HydroPhaseGeometry,
        double liquidViscosity,double gasViscosity,const std::vector<int>& stationaryNoSlipVelocityDofs);
    static HydroPhaseViscosity buildVofCellAverage(const QuadtreeMesh&, HydroPhaseGeometry,
        double liquidViscosity,double gasViscosity,const std::vector<int>& stationaryNoSlipVelocityDofs,
        const HydroPhaseViscosity& previous);
    // Optional geometry-only stock on this actual Momentum owner. Original
    // phase/wall/current-mu validation and K assembly still run every call.
    static HydroPhaseViscosity buildVofCellAverage(const QuadtreeMesh&, HydroPhaseGeometry,
        double liquidViscosity,double gasViscosity,const std::vector<int>& stationaryNoSlipVelocityDofs,
        const HydroMomentumGeometry& actualMomentumGeometry);
    MaterialModel materialModel() const noexcept;
    const HydroPhaseGeometry& geometry() const noexcept;
    // Only an actual Hermite phase-quadrature source sets this flag.
    bool hasQuadrature() const noexcept;
    // Cell-average and exact zero forms have no Hermite source; then rejects.
    const Interface::HermiteGalsPhaseQuadrature& quadrature() const;
    const HydroEnergyViscosity& viscosity() const noexcept;
    double liquidViscosity() const noexcept;
    double gasViscosity() const noexcept;
    void validateFor(const QuadtreeMesh&,Interface::HermiteGalsIdentity,LAMP::TimeLayer) const;
    std::vector<double> apply(const QuadtreeMesh&,const HydroNativeValues&,
        Interface::HermiteGalsIdentity,LAMP::TimeLayer) const;
    double dissipation(const QuadtreeMesh&,const HydroNativeValues&,
        Interface::HermiteGalsIdentity,LAMP::TimeLayer) const;
private:
    struct Data;
    static HydroPhaseViscosity buildVofCellAverageImpl(const QuadtreeMesh&,HydroPhaseGeometry,
        double liquidViscosity,double gasViscosity,const HydroPhaseViscosity* previous,
        const HydroStrainOperator* canonicalStressGeometry=nullptr,
        const std::vector<int>* stationaryNoSlipVelocityDofs=nullptr,
        const HydroMomentumGeometry* rawGeometryOwner=nullptr);
    explicit HydroPhaseViscosity(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
};
}
