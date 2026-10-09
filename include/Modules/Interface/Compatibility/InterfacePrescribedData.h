#pragma once

#include "Modules/Interface/Compatibility/InterfaceTransportInput.h"
#include "Core/Geometry/LevelSetBoundaryTrace.h"
#include "Modules/Interface/Compatibility/InterfaceParametricTransport.h"
#include "Core/Base/ScalarExpression.h"

class QuadtreeMesh;
class NativeLeafDomain;
class InterfaceState;

// Analytic input data for independent interface calculations, not another
// advection/volume/NS solver. Code length/time units; r=0 is the symmetry axis.
namespace InterfacePrescribedData
{
// Reversible unit-square incompressible vortices. Analytic data only; the
// retained spatial flux authority is shared by all time layers. Mesh/domain
// must outlive this object and every returned normal-flux view.
class ReversibleVortex
{
public:
    enum class Profile {ClosedWall,SlipReversal};
    ReversibleVortex(const QuadtreeMesh&,const NativeLeafDomain&,double period,Profile=Profile::ClosedWall);
    Vector2D sample(Vector2D point,double time) const;
    LevelSetBoundaryTrace::VelocityJet sampleJet(Vector2D point,double time) const;
    InterfaceNormalFluxView normalFlux(double time) const;
    InterfaceNormalFluxView spatialNormalFlux() const {return *spatial_;}
    InterfaceParametricTransport::Flow parametricSpatialFlow() const;
    double parameterIncrement(double begin,double end) const;
private:
    double amplitude(double time) const;
    double period_;
    Profile profile_;
    std::shared_ptr<const InterfaceNormalFluxView> spatial_;
};
// Analytic initial data: a closed Hermite ellipse with its analytic arc areas.
// Only initialization seeds an inventory from this geometry.
InterfaceParametricProjection::Curve parametricEllipse(Vector2D center,Vector2D radii,int segments);
// Analytic initial coordinates on s in [0,1]. XY is a C1-periodic loop;
// RZ is a regular meridian from the lower axis pole to the upper one.
// Formula values/jets and midpoint values define one quartic Hermite curve.
InterfaceParametricProjection::Curve parametricAnalytic(const Quantum::ScalarExpression& first,
    const Quantum::ScalarExpression& second,int segments,bool axisymmetric);
// A seed selects one regular bounded F<0 component. Trace tolerance controls
// initial zero-set continuation, independently of the final Hermite segments.
// XY traces a closed loop; RZ starts at its lower regular axis pole.
InterfaceParametricProjection::Curve parametricImplicit(const Quantum::ScalarExpression& field,
    const Quantum::Geometry::Parametric::DomainRectangle& domain,Vector2D seed,
    int segments,bool axisymmetric,double traceTolerance);
// Physical RZ meridian with two regular axis poles; no material axis segment.
InterfaceParametricProjection::Curve parametricSpheroid(double centerZ,Vector2D radii,int segments);
// Connected analytic necked meridian, used only to initialize the same
// authoritative Hermite representation. ratio is r_neck / radial scale.
InterfaceParametricProjection::Curve parametricNeckedSpheroid(double centerZ,Vector2D radii,double ratio,int segments);
// Free paraboloidal cap from the fixed south inlet rim to its regular axis
// pole. Axis and inlet disk close phase measure only, never surface energy.
InterfaceParametricProjection::Curve parametricSupplyCap(double radius,double height,int segments,double bottom=0);
class AxialColumn
{
public:
    AxialColumn(double radius, double peakSpeed);
    // u=(0,A(1-r^2/R^2)^5) inside R, zero outside. Divergence-free,
    // regular on the axis, with a C4 stationary continuation at the wall.
    Vector2D sample(const Vector2D& point) const;
    // SAME analytic velocity, integrated with the full 2*pi*r measure.
    // The returned flux borrows mesh/domain, while velocity copies own data.
    InterfaceNormalFluxView normalFlux(const QuadtreeMesh&, const NativeLeafDomain&) const;
private:
    double radius_, speed_;
};

// phi=z-pin.z-height*(1-r^2/pin.r^2), with the existing fixed wet-corner
// sector. Uses Core's cubic source and panel checks, not a curve fit. This
// initial-data capability is uniform RZ with one vertical wet wall at pin.
// Other mesh/initial-field combinations must use their supported initializer.
void initializeParabolicCap(const QuadtreeMesh&, const InterfaceTransportDomain&,
    Vector2D pin, double height, InterfaceState& freshDestination);
}
