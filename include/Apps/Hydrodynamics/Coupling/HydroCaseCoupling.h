#pragma once

#include <memory>
#include <stdexcept>
#include "Modules/Interface/Compatibility/InterfaceNeckCap.h"
struct HydroCaseSchedule;

class DictionaryGeometry;
class DictionaryMesh;
class QuadtreeMesh;
class HydroInterfaceLoadGeometry;
class InterfaceNodalMaterialMap;
class Hydrodynamics;
struct HydrodynamicsConfig;
struct InterfaceNodalConfig;
template<typename Mesh> class ImmersedGeometry;
template<typename T> class CellField;
template<typename T> class VEdgeField;
template<typename T> class HEdgeField;
class InterfaceSurfaceTransfer;
class NativeLeafDomain;
namespace Quantum::Hydrodynamics {class HydroConservativePressure;struct ConservativeTractionGeometry;struct HydroUniformLoad;}

// Owned unpublished coupling state. Destruction is rollback; publication is
// a nonthrowing handle swap, after all flow/geometry acceptance checks.
class HydroConservativeCouplingCandidate
{
public:
    virtual ~HydroConservativeCouplingCandidate()=default;
    virtual const Quantum::Hydrodynamics::HydroUniformLoad& load() const=0;
    virtual void commit() noexcept=0;
};

// App-owned coupling protocol for the shared fixed-grid Hydro case stages.
// It borrows the same mesh and geometry; it must not build another mesh or
// advance a separate clock. No electric model is part of this interface.
// The provider outlives HydroCase; field/geometry borrows are call-scoped.
class HydroCaseCoupling
{
public:
    virtual ~HydroCaseCoupling() = default;

    // Geometry capability request, not an electric model in the flow layer.
    virtual bool needsTangentialLoadGeometry() const { return false; }
    // Moving charge consumes the SAME Interface transport's material mapping.
    virtual bool needsInterfaceMaterialMap() const { return false; }
    virtual bool supportsConservative() const {return false;}
    // Preserve the already accepted material source before an event. No
    // physical/history advancement; the App uses its pre-event IO stage.
    virtual void writeConservativeEventSource(int,double) {}
    virtual void initializeConservative(const QuadtreeMesh&,const NativeLeafDomain&,const InterfaceState&)
    {throw std::logic_error("coupling does not support conservative initialization");}
    virtual std::unique_ptr<HydroConservativeCouplingCandidate> prepareConservative(const QuadtreeMesh&,
        const NativeLeafDomain&,const InterfaceState&,const InterfaceState&,const InterfaceSurfaceTransfer&,
        const InterfaceNeckCap::MaterialTransfer*,const InterfaceStepTime&,
        const Quantum::Hydrodynamics::HydroConservativePressure*,
        const Quantum::Hydrodynamics::ConservativeTractionGeometry&)
    {throw std::logic_error("coupling does not support joint conservative candidates");}

    // Complete physical coupling after Hydro and Interface, before acceptance
    // and new-time output. Consume staged t^n currents here; do not write
    // completed-step fields/history. Retain only owned results for publication.
    // Exceptions stop the case without accepting the full physical step.
    // The actual Hydro substep and post-advection geometry/material map.
    // App adapters request an owning named-time velocity sampler only when
    // needed; raw mutable Hydro storage never crosses this coupling port.
    // References are borrowed for this call only. This is not joint acceptance.
    virtual void completeStepWithVelocity(const QuadtreeMesh&,
        const ImmersedGeometry<QuadtreeMesh>*, int, const InterfaceNodalMaterialMap*,
        const Hydrodynamics&,
        int, double, double) {}

    // Called once after the App accepts all physical substeps. No borrowed
    // geometry/velocity/map is available here. Emit the owned completed result;
    // an IO failure stops the case but does not undo physical acceptance.
    virtual void publishCompletedStep(int, double, double) {}

    // Called before mesh construction and before creating any output. This
    // must remain read-only, including when --check-case is requested.
    // The immutable schedule may be retained by shared ownership; do
    // not copy run settings into physical config or start another clock.
    virtual void preflight(const HydrodynamicsConfig& config, const InterfaceNodalConfig& interface,
        const DictionaryGeometry& geometry, const DictionaryMesh& mesh,
        int interfaceIndex, std::shared_ptr<const HydroCaseSchedule> schedule) = 0;

    // Called after the shared level set, material fields and initial velocity
    // are ready. Geometry may be null for a constant-material Hydro case;
    // providers requiring an interface must reject it during preflight.
    virtual void initialize(const QuadtreeMesh& mesh,
        const ImmersedGeometry<QuadtreeMesh>* geometry, int interfaceIndex,
        const HydrodynamicsConfig& config) = 0;

    // Called at t^n before the momentum solve and before interface advection.
    // fx/fy are mesh-sized, zero-initialized force densities (force per volume).
    // Hydro's load assembly checks and adds them to the capillary force.
    // Providers may write their own diagnostics and independent field frames
    // here using the shared forceTime. Step 1, forceTime 0 supplies initial
    // electric output; never label these fields with the later Hydro time.
    // Exceptions are recorded as coupling-force failures by the App-owned case step.
    virtual void computeForce(const QuadtreeMesh& mesh,
        const ImmersedGeometry<QuadtreeMesh>* geometry, int interfaceIndex,
        int step, double forceTime, double dt,
        CellField<double>& fx, CellField<double>& fy) = 0;

    // Optional staggered force densities, read immediately AFTER computeForce
    // for the same step/time. This must reuse that solve, never advance or solve
    // a second time. Outputs are zero-initialized x/y face fields; return true
    // when supplying them. A face provider normally leaves its cell outputs
    // zero to avoid applying the same load twice. Existing cell-only providers
    // inherit false unchanged. No density or dual-volume factor belongs here.
    virtual bool computeFaceForce(const QuadtreeMesh&,
        const ImmersedGeometry<QuadtreeMesh>*, int,
        int, double, double, VEdgeField<double>&, HEdgeField<double>&)
    {
        return false;
    }

    // Opt-in immutable load geometry, owned by the App's case state for this
    // force/solve phase only. The context and all borrowed views must NOT be
    // retained past the call/step or used after interface advection. Existing
    // providers retain their original paths through these default forwards.
    virtual void computeForceWithGeometry(const QuadtreeMesh& mesh,
        const ImmersedGeometry<QuadtreeMesh>* geometry, int interfaceIndex,
        int step, double forceTime, double dt, const HydroInterfaceLoadGeometry&,
        CellField<double>& fx, CellField<double>& fy)
    {
        computeForce(mesh, geometry, interfaceIndex, step, forceTime, dt, fx, fy);
    }

    // Reads the SAME electric solve as computeForceWithGeometry. Supplies
    // electric load only; capillarity is already provided by HydroInterfaceForce.
    virtual bool computeFaceForceWithGeometry(const QuadtreeMesh& mesh,
        const ImmersedGeometry<QuadtreeMesh>* geometry, int interfaceIndex,
        int step, double forceTime, double dt, const HydroInterfaceLoadGeometry&,
        VEdgeField<double>& fx, HEdgeField<double>& fy)
    {
        return computeFaceForce(mesh, geometry, interfaceIndex, step, forceTime, dt, fx, fy);
    }
};
