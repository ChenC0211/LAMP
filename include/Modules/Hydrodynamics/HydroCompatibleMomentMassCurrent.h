#pragma once
#include "Modules/Hydrodynamics/HydroCompatibleMomentSpline.h"
#include "Modules/Hydrodynamics/HydroCompatibleMomentMaterial.h"
#include "Modules/Hydrodynamics/HydroCompatibleMomentTransport.h"

class InterfaceState;
struct InterfaceStepTime;
namespace Quantum::Hydrodynamics {
// Time-integrated VOLUME mass current belonging to the same phase ledger.
// This is a local integration object, not another velocity or time solver.
// Its gradient moments satisfy weak continuity for the complete native
// velocity-product space. The orthogonal remainder is the midpoint rho*u.
class HydroCompatibleMomentMassCurrent {
public:
    struct MaterialPoint {Vector2D point;double massWeight=0,density=0;};
    struct CurrentPoint {Vector2D point,weight;};
    struct FacePhase {
        double total=0,totalPositive=0,totalNegative=0,liquidPositive=0,liquidNegative=0;
    };
    struct Cell {
        int index=-1,radialDegree=0,axialDegree=16;
        Vector2D lower,size;
        bool rational=false;
        HydroCompatibleMomentSpace::Partition partition;
        // Tensor shifted Legendre space; in off-axis RZ, radial modes
        // k>0 are divided by (r/r_min)^2; mode zero is exactly constant.
        Eigen::VectorXd beforeMoments,afterMoments,boundaryMoments,referenceMoments,potential;
        std::vector<MaterialPoint> before,midpoint,after;
        // Integrate j_dt dot grad(f) by summing weight dot grad(f(point)).
        std::vector<CurrentPoint> current,correction;
        double constantDefect=0,gradientDefect=0,referenceDefect=0;
        double correctionNorm=0,referenceNorm=0,quadratureEstimate=0;
    };
    static HydroCompatibleMomentMassCurrent build(const QuadtreeMesh&,const NativeLeafDomain&,
        const HydroCompatibleMomentSpace&,const HydroCompatibleMomentSpline::NativeVelocity& midpointVelocity,
        const InterfaceState& previous,const InterfaceState& next,const InterfaceStepTime&,
        const HydroCompatibleMomentMaterial& midpointMaterial,double tolerance=2e-12,
        const HydroCompatibleMomentComplement* complement=nullptr);
    // Constant density with no moving interface. A fixed viscosity interface
    // may remain in the material, but it does not change the mass current.
    // Full velocity-product face moments
    // and the volume current use exactly the same integration/continuity
    // construction as build(); no fictitious phase state is initialized.
    static HydroCompatibleMomentMassCurrent buildConstantDensity(const QuadtreeMesh&,const NativeLeafDomain&,
        const HydroCompatibleMomentSpace&,const HydroCompatibleMomentSpline::NativeVelocity& midpointVelocity,
        const InterfaceStepTime&,const HydroCompatibleMomentMaterial& midpointMaterial,double tolerance=2e-12,
        const HydroCompatibleMomentComplement* complement=nullptr);
    void validate(const QuadtreeMesh&,const NativeLeafDomain&,const HydroCompatibleMomentSpace&) const;
    double time() const noexcept{return time_;}
    const std::vector<Cell>& cells() const noexcept{return cells_;}
    const HydroCompatibleMomentTransport::PhaseMassFlux& massFlux() const noexcept{return mass_;}
    const std::vector<FacePhase>& facePhase() const noexcept{return facePhase_;}
    bool tracksInterface() const noexcept{return tracksInterface_;}
    double inventoryClosure() const noexcept{return inventoryClosure_;}
private:
    static HydroCompatibleMomentMassCurrent assemble(const QuadtreeMesh&,const NativeLeafDomain&,
        const HydroCompatibleMomentSpace&,const HydroCompatibleMomentSpline::NativeVelocity&,
        const Geometry::Parametric::Boundary& previous,const Geometry::Parametric::Boundary& next,
        const InterfaceStepTime&,const HydroCompatibleMomentMaterial&,double,const HydroCompatibleMomentComplement*,
        HydroCompatibleMomentTransport::PhaseMassFlux,std::vector<FacePhase>,bool tracksInterface);
    std::vector<Cell> cells_;
    double inventoryClosure_=0;
    double time_=0;
    bool tracksInterface_=true;
    std::shared_ptr<const NativeMeshIdentity> identity_;
    HydroCompatibleMomentTransport::PhaseMassFlux mass_;
    std::vector<FacePhase> facePhase_;
};
}
