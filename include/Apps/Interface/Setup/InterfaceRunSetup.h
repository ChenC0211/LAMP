#pragma once

#include "Interface/IO/InterfaceTransportInput.h"
#include "Core/Base/ScalarExpression.h"
#include "Core/Geometry/PolynomialVelocitySlab.h"
#include "Core/Geometry/LevelSetBoundaryTrace.h"
#include <array>
#include <memory>

class QuadtreeMesh;

// Expressions own immutable compiled functions. Time coefficients feed the
// existing polynomial endpoint provider, including its gradient/CFL contract.
class InterfacePrescribedData
{
public:
    explicit InterfacePrescribedData(const InterfaceTransportConfig&);
    LAMP::ScalarExpression::Jet initial(Vector2D) const;
    LAMP::ScalarExpression::Jet boundary(Vector2D,double time) const;
    LAMP::Geometry::PolynomialVelocityData polynomial(double previous,double arrival) const;
private:
    LAMP::ScalarExpression initial_;
    std::shared_ptr<const LAMP::ScalarExpression> boundary_;
    Vector2D origin_,scale_;
    std::array<std::array<std::shared_ptr<const LAMP::ScalarExpression>,2>,16> coefficients_;
};
// Typed preparation of the original mesh and boundary sampling support.
// The App owns the accepted publication, main clock and mesh-transfer commit.
std::unique_ptr<QuadtreeMesh> makeInterfaceTransportMesh(const InterfaceTransportConfig&,bool amr);
LevelSetBoundaryTrace::Domain interfaceTransportTraceDomain(const QuadtreeMesh&,bool radial);
