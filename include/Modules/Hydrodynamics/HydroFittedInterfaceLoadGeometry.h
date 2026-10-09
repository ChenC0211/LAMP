#pragma once

#include <memory>
#include <vector>
#include "Modules/Hydrodynamics/HydroInterfaceLoadGeometry.h"
#include "Core/Geometry/InterfaceTraceGeometry.h"

class QuadtreeMesh;
class HydroFittedDomain;
class LevelSetInterfacePanels;
class NativeMeshIdentity;

namespace Quantum::Hydrodynamics
{
// An explicit scalar-gradient closure, NOT a velocity boundary condition.
// Missing transverse source gradients on actual no-slip walls use the unique
// active-side Hermite jet: H'(phi) times the physical Cartesian derivative.
// Interior sources retain the fitted C(H)/dualVolume difference. This is a
// new boundary discretization, not a claim of exact near-wall force balance.
enum class FittedTransverseGradientClosure { ActiveHermiteWallGradient };

struct HydroFittedInterfaceLoadGeometry
{
    std::vector<double> indicator; // Native cell indices; inactive values are NaN.
    // Native free faces with nonzero primary OR transverse gradient. The same
    // retained local Frame is intended for capillarity and electric traces.
    const std::vector<HydroInterfaceLoadGeometry::Face>& faces() const noexcept { return faces_; }
    const std::shared_ptr<const InterfaceTraceGeometry>& traceGeometry() const noexcept { return traceGeometry_; }
    // A completed builder owns the full pressure-gradient union. Only its
    // identical ordered trace collection can be loaded on its receiving mesh.
    void validate(const QuadtreeMesh&, const HydroFittedDomain&,
        const std::shared_ptr<const InterfaceTraceGeometry>& traceGeometry) const;
    double halfWidth = 0.0;
    double maxGradientIdentityError = 0.0;
    double maxAbsWallGradient = 0.0;
    double maxAbsTransverseGradient = 0.0;
    // Distinct source faces actually consumed by the transverse interpolation.
    int wallGradientSources = 0;
    int nonzeroWallGradientSources = 0;
    int axisGradientSources = 0;
    int pressurePlateauSources = 0;
    int transverseOnlyFaces = 0;
    int endpointProjections = 0;
    int nativeTransverseFits = 0;
    double maxTransverseAffineError = 0.0;
private:
    std::vector<HydroInterfaceLoadGeometry::Face> faces_;
    std::shared_ptr<const InterfaceTraceGeometry> traceGeometry_;
    std::shared_ptr<const NativeMeshIdentity> meshIdentity_;
    static HydroFittedInterfaceLoadGeometry build(const QuadtreeMesh&, const HydroFittedDomain&,
        std::shared_ptr<const LevelSetHermiteGeometry::Snapshot>, const LevelSetInterfacePanels&,
        double, FittedTransverseGradientClosure, bool nativeLeaves);
    friend HydroFittedInterfaceLoadGeometry buildFittedInterfaceLoadGeometry(
        const QuadtreeMesh&, const HydroFittedDomain&, std::shared_ptr<const LevelSetHermiteGeometry::Snapshot>,
        const LevelSetInterfacePanels&, double, FittedTransverseGradientClosure);
    friend HydroFittedInterfaceLoadGeometry buildNativeFittedInterfaceLoadGeometry(
        const QuadtreeMesh&, const HydroFittedDomain&, std::shared_ptr<const LevelSetHermiteGeometry::Snapshot>,
        const LevelSetInterfacePanels&, double, FittedTransverseGradientClosure);
};

// Explicit active RZ, same-level leaves, retained shared-source open panels.
// No sigma, epsilon, charge or electric field enters this geometric assembly.
// The PRIMARY gradient retains the old fitted capillary C(H)/dualVolume and
// its independent 512-epsilon identity gate. The cross component is the old
// four-face average wherever all four sources are interior. At a true r=0
// symmetry face its radial source is zero by the even scalar, not by u_r=0.
// Pressure-port owner cells must remain entirely outside the load band.
// Constrained velocity rows receive no load; their scalar-gradient traces
// remain distinct data and may contribute to neighboring free-row loads.
//
// This does not change buildFittedCapillaryForce's legacy normal-only contract
// or its projection arithmetic. Wall reactions/contact-line forces are not
// supplied. All queries and checks finish before a result is returned.
HydroFittedInterfaceLoadGeometry buildFittedInterfaceLoadGeometry(
    const QuadtreeMesh& mesh, const HydroFittedDomain& domain,
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> snapshot,
    const LevelSetInterfacePanels& retainedPanels, double halfWidth,
    FittedTransverseGradientClosure transverseClosure);

// Explicit fixed 2:1 leaf adapter. Primary loads use the SAME fitted C(H),
// including the old virtual-child pressure stencil and its actual dual volume.
// Transverse values use the old linear least-squares interpolation on all real
// opposite-component subfaces of the two incident active leaves. Its fixed
// incidence stencil must reproduce constants/affine fields; it is not enlarged
// or replaced by a velocity boundary condition. This adaptation does not claim
// exact tangential stress balance or enable runtime refinement/coarsening.
// Independent native-face frame queries may run in parallel on immutable
// geometry. Face ordering, scalar-gradient arithmetic and diagnostic sums
// remain deterministic; worker failures propagate before any result returns.
HydroFittedInterfaceLoadGeometry buildNativeFittedInterfaceLoadGeometry(
    const QuadtreeMesh& mesh, const HydroFittedDomain& domain,
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> snapshot,
    const LevelSetInterfacePanels& retainedPanels, double halfWidth,
    FittedTransverseGradientClosure transverseClosure);
}
