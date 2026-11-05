#pragma once

#include "Modules/Hydrodynamics/HydroPhaseGeometry.h"
#include "Modules/Hydrodynamics/HydroPressureLoad.h"
#include <memory>
#include <vector>

namespace LAMP::Hydrodynamics
{
class HydroAdvectorReconstruction;
class HydroAdvectorBoundary;
// Native pressure-compatible normal stress. The constant factory consumes
// measured FirstJet alpha. Capillarity uses same-scalar pressure-center phase
// and physical curvature at the real interface on each native center link.
// This capillary pressure is a cell-center phase value, not a cell-volume
// average. This normal-jump discretization does not add a viscous GFM route.
class HydroInterfaceNormalLoad
{
public:
    static HydroInterfaceNormalLoad constant(const QuadtreeMesh&, HydroPressureLoad,
        HydroPhaseGeometry, double jump);
    // Physical sigma >= 0. FirstJet uses an outward normal, whereas the
    // frozen HYDRO stress coefficient uses the inward normal: J=-sigma*kappa.
    // Native center links use the SAME accepted scalar for the phase jump
    // and the actual zero-crossing curvature; the stress matrix is diagonal.
    static HydroInterfaceNormalLoad capillary(const QuadtreeMesh&, HydroPressureLoad,
        HydroPhaseGeometry, double surfaceTension,
        const HydroAdvectorReconstruction&, const HydroAdvectorBoundary* =nullptr);
    const HydroPressureLoad& pressureLoad() const noexcept;
    const HydroPhaseGeometry& geometry() const noexcept;
    const HydroDofLayout& layout() const noexcept;
    bool isCapillary() const noexcept;
    double constantJump() const; // rejects for capillary loads
    double surfaceTension() const; // capillary only
    const std::vector<double>& phaseIndicators() const; // pressure-center phi<0:1, phi>0:0, phi==0:.5
    const std::vector<double>& linkCurvatures() const; // same-scalar curvature at active native link zero
    const std::vector<double>& linkLengths() const; // actual active pressure-center link length
    const HydroPressureLoad::Matrix& stressMatrix() const noexcept; // diagonal J_f; inactive capillary rows zero
    const HydroAdvectorReconstruction& reconstruction() const; // capillary only, retained SAME preparation
    const HydroAdvectorBoundary* boundary() const; // explicit virtual normal law; capillary only
    const std::vector<double>& constantReactions() const; // capillary internal links have no prescribed-wall load
    // Constant: F=C(J0*alpha). Capillary: F=diag(-sigma*kappa_link)*C*chi_center.
    // Both are balance terms in Cp + F = 0.
    // For A*u + Cp = rhs, the caller subtracts F from rhs.
    const std::vector<double>& balanceLoads() const noexcept;
    // Independent curl modes have zero interface normal trace; their normal
    // stress balance rows are exactly zero, including the constant channel.
    const std::vector<double>& tangentialBalanceLoads() const;
    void validateFor(const QuadtreeMesh&, const HydroPressureLoad& actualPressure,
        Interface::HermiteGalsIdentity actualGeometry,
        Interface::HermiteGalsTimeLayer actualClock) const;
    void validateFor(const QuadtreeMesh&, const HydroPressureLoad&,
        Interface::HermiteGalsIdentity, Interface::HermiteGalsTimeLayer,
        const HydroAdvectorReconstruction&, const HydroAdvectorBoundary*) const;
private:
    struct Data;
    explicit HydroInterfaceNormalLoad(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
};
}
