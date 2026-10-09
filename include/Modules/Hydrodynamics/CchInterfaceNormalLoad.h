#pragma once

#include "Modules/Hydrodynamics/CchPhaseGeometry.h"
#include "Modules/Hydrodynamics/CchPressureLoad.h"
#include <memory>
#include <vector>

namespace Quantum::Hydrodynamics
{
class CchAdvectorReconstruction;
class CchAdvectorBoundary;
namespace Detail
{
// Derived algebra of this same publication and prepared source. No new field
// or state owner. Samples follow publication curve/graph/quadrature order,
// including scalar support outside the selected Hydro domain.
struct CchInterfaceStressAssembly
{
    CchPressureLoad::Matrix matrix;
    std::vector<double> load, averages, surfaceMeasures, anchorMeasures;
    std::vector<double> constantReactions;
    double constantChannel=0;
};
CchInterfaceStressAssembly assembleCchInterfaceStress(const QuadtreeMesh&,
    const CchPressureLoad&, const CchPhaseGeometry&,
    const CchAdvectorReconstruction&, const CchAdvectorBoundary*,
    const std::vector<double>& normalStress);
}
// One pressure-compatible normal-stress deposition owner. Both factories
// consume the same FirstJet alpha and C. No material inventory, second
// gradient, clock or accepted-state update.
class CchInterfaceNormalLoad
{
public:
    static CchInterfaceNormalLoad constant(const QuadtreeMesh&, CchPressureLoad,
        CchPhaseGeometry, double jump);
    // Physical sigma >= 0. FirstJet uses an outward normal, whereas the
    // frozen CCH stress coefficient uses the inward normal: J=-sigma*kappa.
    // Complete same-publication stress is tested by the SAME source's local
    // normal basis. Positive local anchors preserve the constant channel.
    static CchInterfaceNormalLoad capillary(const QuadtreeMesh&, CchPressureLoad,
        CchPhaseGeometry, double surfaceTension,
        const CchAdvectorReconstruction&, const CchAdvectorBoundary* =nullptr);
    const CchPressureLoad& pressureLoad() const noexcept;
    const CchPhaseGeometry& geometry() const noexcept;
    const CchDofLayout& layout() const noexcept;
    bool isCapillary() const noexcept;
    double constantJump() const; // rejects for capillary loads
    double surfaceTension() const; // capillary only
    const std::vector<double>& curvatureAverages() const; // local positive anchors; zero without real support
    const std::vector<double>& supportMeasures() const; // ds or 2*pi*r ds on actual source footprint
    const CchPressureLoad::Matrix& stressMatrix() const noexcept; // complete native J_f, including inactive rows
    const CchAdvectorReconstruction& reconstruction() const; // capillary only, retained SAME preparation
    const CchAdvectorBoundary* boundary() const; // explicit virtual normal law; capillary only
    const std::vector<double>& constantReactions() const; // known prescribed-wall constant channel
    // Constant: F=C(J0*alpha). Capillary: local anchored weak stress, F=J_f*C*alpha.
    // Both are balance terms in Cp + F = 0.
    // For A*u + Cp = rhs, the caller subtracts F from rhs.
    const std::vector<double>& balanceLoads() const noexcept;
    // Independent curl modes have zero interface normal trace; their normal
    // stress balance rows are exactly zero, including the constant channel.
    const std::vector<double>& tangentialBalanceLoads() const;
    void validateFor(const QuadtreeMesh&, const CchPressureLoad& actualPressure,
        Interface::HermiteGalsIdentity actualGeometry,
        Interface::HermiteGalsTimeLayer actualClock) const;
    void validateFor(const QuadtreeMesh&, const CchPressureLoad&,
        Interface::HermiteGalsIdentity, Interface::HermiteGalsTimeLayer,
        const CchAdvectorReconstruction&, const CchAdvectorBoundary*) const;
private:
    struct Data;
    explicit CchInterfaceNormalLoad(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
};
}
