#pragma once

#include "Core/Base/Vector2D.h"
#include "Modules/Electric/ElectricPhaseGeometry.h"
#include <array>
#include <cstddef>
#include <memory>
#include <vector>

class QuadtreeMesh;
namespace LAMP::Interface { class HermiteGalsPublication; }

namespace LAMP::Electric
{
class ElectricPhaseGeometry;
class ElectricPhaseResult;

struct ElectricMaxwellControls
{
    double relativeProjectionTolerance = 1e-11;
    double relativeChargeJumpTolerance = 1e-10;
};
struct ElectricMaxwellReport
{
    std::size_t interfaceFaces = 0, quadratureNodes = 0;
    // Maximum of the cubic shared-P2 constraint and weighted q_h mass solve.
    double maxRelativeProjectionResidual = 0;
    double maxRelativeOptimalityResidual = 0;
    // Original electric recovery's three-moment jump receipt and arithmetic scale.
    double maxRelativeChargeJumpResidual = 0;
};
struct ElectricMaxwellSample
{
    Vector2D point, normal, tangent, tractionJump;
    double coordinate = 0, physicalWeight = 0;
    std::array<double,2> normalDisplacement{{0,0}};
    double tangentialField = 0, chargeDensity = 0;
    double normalTraction = 0, tangentialTraction = 0;
};
// Read-only (T_plus-T_minus)*n on the original shared Gamma skeleton.
// Normal displacement is the equilibrated electric flux; the ONE tangential field
// comes from the shared P2 trace plus a symmetric, weighted cubic completion.
// It neither solves a second electric field nor deposits a Hydro load. This
// first recovery is a continuous tensor evaluation, not an energy shape derivative.
class ElectricMaxwellTrace
{
public:
    static ElectricMaxwellTrace prepare(const ElectricPhaseResult&,const ElectricMaxwellControls& = {});
    const ElectricPhaseGeometry& geometry() const;
    // The physical electric clock and immutable geometry source clock may
    // differ in an explicitly fixed-geometry P4 experiment.
    Interface::HermiteGalsTimeLayer timeLayer() const;
    Interface::HermiteGalsTimeLayer geometrySourceTimeLayer() const;
    const ElectricMaxwellReport& report() const;
    const std::vector<std::size_t>& interfaceFaces() const;
    // Every face argument is geometry().faces()'s original ID. Other faces reject.
    // Coefficients use oriented-arclength Legendre P0,P1,P2,P3 / P0,P1,P2.
    std::array<double,4> potentialPolynomialP3(std::size_t face) const;
    std::array<double,3> tangentialFieldPolynomialP2(std::size_t face) const;
    const std::vector<ElectricMaxwellSample>& quadrature(std::size_t face) const;
    // Same P2 Dn/q and shared cubic Et at an original certified graph sample.
    // A point query has physicalWeight=0; it is not an integration rule.
    ElectricMaxwellSample sample(std::size_t component,std::size_t graph,
        const Interface::HermiteGalsSurfaceSample&) const;
    // Coeval read-only interface view of the original physical traction jump.
    // A fixed-geometry field with a different physical clock cannot supply it.
    Interface::HermiteGalsSurfaceVectorView surfaceVectorView() const;
    // Integrates the same sample vector with ds (XY) or 2*pi*r*ds (RZ).
    // The RZ radial entry is meridional; it is not a Cartesian 3D net force.
    Vector2D integratedTractionJump(std::size_t face) const;
    void validateSource(const QuadtreeMesh&,const Interface::HermiteGalsPublication&,
        const Interface::HermiteGalsTimeLayer& electricLayer) const;
private:
    struct Data;
    explicit ElectricMaxwellTrace(std::shared_ptr<const Data>);
    const Data& checked() const;
    std::shared_ptr<const Data> data_;
};
}
