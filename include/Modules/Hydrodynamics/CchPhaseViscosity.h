#pragma once

#include "Modules/Hydrodynamics/CchEnergyViscosity.h"
#include "Modules/Hydrodynamics/CchPhaseGeometry.h"
#include "Modules/Interface/HermiteGalsPhaseQuadrature.h"
#include <memory>

namespace Quantum::Hydrodynamics
{
class CchKineticHodge;
// Same-geometry phase quadrature -> native strain -> nonnegative phase mu.
// No blended viscosity, second strain implementation or material inventory.
class CchPhaseViscosity
{
public:
    // Same physical source and complete phase-volume partition, including U-q
    // cross terms. This bulk form does not supply the finite-mu jump trace.
    static CchEnergyViscosity::Matrix modalBulkMatrix(const QuadtreeMesh&,const CchKineticHodge&,
        double liquidViscosity,double gasViscosity);
    static CchPhaseViscosity build(const QuadtreeMesh&, CchPhaseGeometry,
        Interface::HermiteGalsPhaseQuadrature,double liquidViscosity,double gasViscosity);
    const CchPhaseGeometry& geometry() const noexcept;
    const Interface::HermiteGalsPhaseQuadrature& quadrature() const noexcept;
    const CchEnergyViscosity& viscosity() const noexcept;
    double liquidViscosity() const noexcept;
    double gasViscosity() const noexcept;
    void validateFor(const QuadtreeMesh&,Interface::HermiteGalsIdentity,Interface::HermiteGalsTimeLayer) const;
    std::vector<double> apply(const QuadtreeMesh&,const CchNativeValues&,
        Interface::HermiteGalsIdentity,Interface::HermiteGalsTimeLayer) const;
    double dissipation(const QuadtreeMesh&,const CchNativeValues&,
        Interface::HermiteGalsIdentity,Interface::HermiteGalsTimeLayer) const;
private:
    struct Data;
    explicit CchPhaseViscosity(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
};
}
