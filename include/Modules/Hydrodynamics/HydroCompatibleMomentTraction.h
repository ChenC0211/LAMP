#pragma once
#include "Modules/Hydrodynamics/HydroCompatibleMomentMaterial.h"
#include "Modules/Hydrodynamics/HydroCompatibleMomentSpline.h"
#include <functional>

namespace Quantum::Hydrodynamics {
// One sharp traction, integrated on the same material geometry/time as Hydro.
// Positive normal stress acts inward; tangent follows the CCW curve. A future
// electro-capillary consumer supplies the already combined normal stress here.
class HydroCompatibleMomentTraction {
public:
    struct Components {double normal=0,tangential=0;};
    using Law=std::function<Components(const Geometry::Parametric::SurfaceSample&)>;
    struct Report {double length=0,meanCurvature=0,curvatureSquaredIntegral=0,meanNormal=0,quadratureError=0;int pieces=0,samples=0;bool clipped=false;};
    struct Result {
        Eigen::VectorXd variation,pressurePotential;
        Report report;
    };
    static Result build(const QuadtreeMesh&,const NativeLeafDomain&,const HydroCompatibleMomentSpace&,
        const HydroCompatibleMomentSpline&,const HydroCompatibleMomentMaterial&,const Law&,double constantSurfaceTension=0);
    // Variation in [q,eta], using the same sharp quadrature and pressure lift.
    static Result build(const QuadtreeMesh&,const NativeLeafDomain&,const HydroCompatibleMomentSpace&,
        const HydroCompatibleMomentSpline::NativeVelocity&,const HydroCompatibleMomentComplement&,
        const HydroCompatibleMomentMaterial&,const Law&,double constantSurfaceTension=0);
    // Law supplies the TOTAL normal stress. Its sigma*kappa part, when sigma
    // is provided, is evaluated by -sigma integral(t . grad(v) t) ds; this is
    // the identical sharp traction in weak form and only uses first jets.
    // RZ uses -2*pi*sigma integral(r*t.grad(v).t+v_r) ds on the physical
    // meridian; the virtual axis closure has no surface energy or load.
    // Report::length denotes physical interface measure: length in XY,
    // surface area in RZ. Curvature is respectively planar/total curvature.
    static Report measure(const Geometry::Parametric::Curve&,bool axisymmetric=false);
    static Report measure(const Geometry::Parametric::Boundary&,bool axisymmetric=false);
    // Physical control-volume measure. Geometry may extend beyond the box;
    // only its original polynomial intervals inside the domain contribute.
    // No free surface in the domain gives zero integral/mean conventions.
    static Report measure(const Geometry::Parametric::Boundary&,const QuadtreeMesh&,const NativeLeafDomain&);
    // Positive principal shape response: sigma*integral grad_Gamma(u.n) .
    // grad_Gamma(v.n), on the accepted geometry for one nonlinear solve.
    // Tangential reparametrization has zero normal-displacement derivative.
    // RZ uses the same scalar normal speed and physical 2*pi*r measure.
    // Used on BOTH sides of Picard; never enters physical residual/energy.
    static HydroCompatibleMomentSpace::Matrix iterationMetric(const QuadtreeMesh&,const NativeLeafDomain&,
        const HydroCompatibleMomentSpace&,const HydroCompatibleMomentSpline::NativeVelocity&,
        const HydroCompatibleMomentMaterial&,double surfaceTension);
private:
    template<class Velocity> static Result buildComponent(const QuadtreeMesh&,const NativeLeafDomain&,const HydroCompatibleMomentSpace&,
        const Velocity&,const HydroCompatibleMomentMaterial&,std::size_t,const Law&,double,bool,
        const HydroCompatibleMomentComplement*);
};
}
