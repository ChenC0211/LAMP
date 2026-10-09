#pragma once
#include "Modules/Hydrodynamics/HydroCompatibleMomentViscosity.h"
#include "Modules/Hydrodynamics/HydroBoundary.h"
#include <stdexcept>

namespace Quantum::Hydrodynamics {
// One immutable partition of native exterior faces. NormalVelocity prescribes
// outward flux with zero tangential traction; PressureOutlet supplies the FULL
// normal stress -p_out n and zero tangential stress. FreeSlip has zero normal
// velocity and natural zero shear traction. Unlisted walls are no-slip.
// The same face identity is consumed by elimination, momentum flux and work.
class HydroCompatibleMomentBoundary {
public:
    enum class Kind {Wall,NormalVelocity,PressureOutlet,FreeSlip};
    struct Face {Kind kind=Kind::Wall;double outwardVelocity=0,pressure=0;
        int incomingPhase=0;std::optional<Vector2D> incomingVelocity;};
    struct NormalReference {bool horizontal=false;double plane=0,lo=0,hi=0,area=0,center=0,scale=0,value=0;Kind kind=Kind::Wall;};
    static HydroCompatibleMomentBoundary build(const QuadtreeMesh&,const NativeLeafDomain&,
        const HydroCompatibleMomentSpace&,const HydroCompatibleMomentViscosity&,
        const std::vector<HydroBoundaryPatch>&);
    // Preserve the ORIGINAL weak normal boundary data on its original
    // physical intervals. Refined child means are not fresh independent BCs.
    HydroCompatibleMomentBoundary inheritNormalConstraints(const QuadtreeMesh&,const NativeLeafDomain&,
        const HydroCompatibleMomentSpace&,const HydroCompatibleMomentBoundary& source) const;
    void validate(const QuadtreeMesh&,const NativeLeafDomain&,const HydroCompatibleMomentSpace&) const;
    const std::vector<Face>& faces() const noexcept{return faces_;}
    // Native fixed-face mask for the legacy lift and homogeneous pressure
    // complement. Complete q uses normalConstraints(), including after AMR.
    const std::vector<bool>& fixed() const noexcept{return fixed_;}
    const Eigen::VectorXd& values() const noexcept{return values_;}
    const Eigen::VectorXd& traction() const noexcept{return traction_;}
    bool open() const noexcept{return open_;}
    bool hasPressureOutlet() const noexcept{return pressureOutlet_;}
    // Free-slip walls need unconstrained tangential traces without opening
    // the domain or fixing its pressure reference.
    bool usesGeneralTraceSpace() const noexcept{return generalTrace_;}
    static bool impermeable(Kind kind) noexcept{return kind==Kind::Wall||kind==Kind::FreeSlip;}
    const HydroCompatibleMomentSpace::Matrix& normalConstraints() const {
        if(!nestedNormal_)throw std::invalid_argument("non-nested normal data require complete trace quadrature");return normalConstraints_;
    }
    const std::vector<NormalReference>& normalReferences() const noexcept{return normalReference_;}
    bool nestedNormalConstraints() const noexcept{return nestedNormal_;}
    const Eigen::VectorXd& normalConstraintValues() const noexcept{return normalValues_;}
private:
    std::vector<NormalReference> normalReference_;
    HydroCompatibleMomentSpace::Matrix normalConstraints_;
    Eigen::VectorXd normalValues_;
    void assembleNormalConstraints(const QuadtreeMesh&,const NativeLeafDomain&,const HydroCompatibleMomentSpace&);
    std::shared_ptr<const NativeMeshIdentity> identity_;
    bool open_=false;
    bool pressureOutlet_=false,generalTrace_=false;
    bool nestedNormal_=true;
    std::vector<Face> faces_;
    std::vector<bool> fixed_;
    Eigen::VectorXd values_,traction_;
};
}
