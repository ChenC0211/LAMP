#pragma once

#include "Core/Field/CellField.h"
#include "Core/Field/VEdgeField.h"
#include "Core/Field/HEdgeField.h"
#include "Core/Geometry/InterfaceTraceGeometry.h"
#include <memory>
#include <vector>

class HydroInterfaceLoadGeometry;
class NativeMeshIdentity;
class NativeLeafDomain;
class InterfaceState;

namespace Quantum::Hydrodynamics
{
class HydroConservativePressure;
// Face values are the authoritative momentum load;
// cell averages and integrals are diagnostics, never a second momentum source.
struct HydroUniformLoad
{
    CellField<double> forceX, forceY;
    VEdgeField<double> faceForceX;
    HEdgeField<double> faceForceY;
    double maxForce = 0.0, totalForceX = 0.0, totalForceY = 0.0;
    explicit HydroUniformLoad(const QuadtreeMesh& mesh);
};

// Planar complete-domain pressure gradient. Native coefficients use x rows
// then y rows; boundary coefficients must be zero. Source identity is checked
// before any scalar/coefficient is read. Both arrays describe that source.
HydroUniformLoad buildUniformScalarGradientLoad(const QuadtreeMesh& mesh,
    const NativeMeshIdentity& source, const std::vector<double>& scalarByCell,
    const std::vector<double>& coefficientByNative);

// CCH normal load: coefficient_f * (C alpha)_f divided by the SAME physical
// momentum dual volume. alpha is liquid inventory/cell volume. With outward
// liquid->gas n, kappa=div(n)>0 for a drop, coefficient=sigma*kappa gives
// p=p0+sigma*kappa*alpha at rest. No second geometry or density factor enters.
// Closed normal boundary coefficients must be zero; viscosity and materials
// remain the caller's step inputs. The supplied pressure operator is retained
// only for this synchronous assembly and must be the adjoint variant.
HydroUniformLoad buildCompatibleScalarGradientLoad(const QuadtreeMesh& mesh,
    const NativeLeafDomain& domain,const HydroConservativePressure& pressure,
    const std::vector<double>& scalarByCell,const std::vector<double>& coefficientByNative);

struct ParametricCapillaryReport {
    double length=0,meanCurvature=0,curvatureSquaredIntegral=0;
    double quadratureError=0;
    std::size_t surfacePieces=0,surfaceSamples=0;
};
// One globally defined balanced surface integral on the unique parameter
// curve: sigma*k0*C(alpha) - sigma*integral W_f*(kappa-k0)*n_i ds,
// where k0=int(kappa ds)/length and W is the SAME C1 spline reconstruction
// which advects the interface. Its exact algebraic transpose deposits force.
// All curvature variation and its normal correlation remain in the integral;
// k0 is a pressure mode, not a curvature filter. Numerical FV/geometry tails
// need neither a nearest-point extension nor a different force provider.
// For a closed curve and a discretely divergence-free native velocity, the
// constant mode does no work. The remaining surface work is paired with the
// actual reconstructed velocity (up to quadrature), not a different RT field.
// Uniform XY unit-square states only; no contact force, RZ or AMR claim.
HydroUniformLoad buildParametricCapillaryLoad(const QuadtreeMesh&,const NativeLeafDomain&,
    const HydroConservativePressure&,const InterfaceState&,double surfaceTension,
    ParametricCapillaryReport* report=nullptr);

// Normal/tangential tractions in the exact retained query order. The same
// query collection must come back from the physical module, not a separately
// reconstructed/reordered geometry. Geometry validates current mesh/phi.
// includeTangential=false retains the original normal-only arithmetic.
HydroUniformLoad buildUniformInterfaceTractionLoad(const QuadtreeMesh& mesh,
    const VertField<double>& scalar, const HydroInterfaceLoadGeometry& geometry,
    const std::shared_ptr<const InterfaceTraceGeometry>& traceGeometry,
    const std::vector<Vector2D>& traction, bool includeTangential);

struct ConservativeCapillaryReport
{
    // Read-only samples of the ACTUAL receiving-face contributions. These
    // belong to the supplied source geometry. In the conservative App this
    // is the prepared end-step interface, before joint physical acceptance.
    struct Face
    {
        int nativeRow=0;
        Vector2D center,interfacePoint;
        double gradient=0,curvature=0,meridional=0,azimuthal=0,force=0;
        // These three flags describe the first contributing knot. A mixed
        // interpolation is counted separately and retains both method IDs.
        bool height=false,evenAxis=false,implicitAxis=false;
        bool interpolated=false;
        int primaryMethod=0,secondaryMethod=0;
        double secondaryWeight=0;
        bool secondaryEvenAxis=false,secondaryImplicitAxis=false;
        bool unresolvedInventoryTail=false;
        double projectionDistance=0,maxMinorityFraction=0;
        std::size_t sameChartTailCandidates=0;
        bool regularizedCurvature=false;
        double secondaryCurvatureWeight=0;
        int regularizationStatus=0;
        std::size_t eventChartCandidates=0;
        double eventCurvatureSpread=0;
        bool boundedTailRidge=false;
        double tailCurvatureSpread=0,scaledTailForceUncertainty=0;
        std::size_t tailAlternatives=0;
        std::uint64_t liquidId=0; // Zero: unsplit/all-surface or unassigned-tail contribution.
        // Integrated faces report weighted mean curvature for diagnostics;
        // force includes its correlation with each signed interface normal.
        bool surfaceIntegrated=false;
        std::size_t surfaceSamples=0;
        double normalMomentDefect=0; // -G(C_k)-integral(W*n_i), units 1/length.
    };
    std::size_t heightFaces=0,localFitFaces=0,mixedFaces=0,eventModelFaces=0;
    double maxEventModelForce=0;
    double maxEventProjectionForceSpread=0;
    std::size_t unresolvedTailFaces=0;
    double maxUnresolvedTailForce=0;
    double maxProjectionDistance=0;
    std::size_t boundedTailRidgeFaces=0;
    double maxScaledTailForceUncertainty=0,maxTailForceUncertainty=0;
    std::size_t surfaceIntegratedFaces=0;
    double maxNormalMomentDefect=0;
    std::vector<Face> faces;
};
// Receiving geometry for Maxwell traction, sharing capillarity's inventory
// moment constraints and pressure support. Includes transverse-only rows; the physical module
// returns the retained trace collection unchanged. C is the liquid fraction.
struct ConservativeTractionGeometry
{
    struct Face {int row;bool horizontal;double gradient,transverse;};
    std::shared_ptr<const NativeMeshIdentity> meshIdentity;
    std::shared_ptr<const InterfaceTraceGeometry> trace;
    std::vector<Face> faces;
    // Primary trace queries are first, in faces order. Additional actual
    // interface queries retain the competing projections of certified tails;
    // Maxwell must evaluate them and pass its OWN traction-error check.
    struct TailAlternative {std::size_t faceIndex=0,queryIndex=0;double resolutionScale=0;};
    std::vector<TailAlternative> tailAlternatives;
    // Event-cell owner supports retain ALL their actual interface quadrature
    // points. Queries follow ordinary faces and tail witnesses; nodes shared
    // by several receiving faces occur only once in the physical trace solve.
    struct Moment {std::size_t queryIndex;long double normal,tangent;};
    struct IntegratedFace {Face face;std::vector<Moment> moments;};
    std::vector<IntegratedFace> integratedFaces;
    std::size_t integrationQueryCount=0;
    std::vector<long double> liquidVolumes;
    std::uint64_t liquidId=0;
    // Split children retain separate G(C_k) rows and query sets. Their loads
    // add even when several children contribute to the same pressure face.
    std::vector<std::shared_ptr<const ConservativeTractionGeometry>> components;
    bool nativePressure=false;
};
struct ConservativeTractionReport
{
    std::size_t tailAlternatives=0;
    double maxScaledTailForceUncertainty=0,maxTailForceUncertainty=0;
};
HydroUniformLoad buildConservativeInterfaceTractionLoad(const QuadtreeMesh&,
    const NativeLeafDomain&,const InterfaceState&,const HydroConservativePressure*,
    const ConservativeTractionGeometry&,const std::shared_ptr<const InterfaceTraceGeometry>&,
    const std::vector<Vector2D>& traction,ConservativeTractionReport* report=nullptr);
// Closed or symmetry-section conservative interfaces, full XY/RZ root mesh.
// Section-normal boundary loads vanish by the declared even C/curvature law;
// the physical caller must bind those sections to impermeable slip boundaries.
// Assemble
// +sigma*kappa*G(C) on the SAME native pressure faces: C is liquid inventory,
// kappa is +div(n_liquid->gas) extended from the one published interface.
// Explicit interior-split children use sum sigma*kappa_k*G(C_k), with unique
// geometric cell ownership and a separate, retained unassigned-tail channel.
// A pressure row may have several reported contributions (liquidId); their
// gradients sum to the original G(C) and their forces sum to the actual load.
// Requested Maxwell geometry retains each component's own inventory support.
// On uniform event-cell owner supports, integrate all retained native surface
// nodes with the unit-flux face basis. Correct its zeroth oriented moments to
// -G(C_k) and the original transverse gradient. Constant traction preserves
// pressure balance; variable traction keeps the signed contributions of folds.
// This is an explicit event-model discretization, not a curvature/grid accuracy
// certificate. Non-event and mixed-mesh faces retain the projection policy.
// Fixed contact/wetting forces and AMR are not provided by this entry.
// An ambiguous/unavailable resolved load query rejects the whole candidate.
// A pure-phase unresolved inventory tail may keep the actual nearest provider
// only when |G(C)|*spread(kappa)*h^2 <= the existing PLIC volume resolution,
// using the complete pressure support. All competing providers must exist;
// the whole spread is recorded. No C clipping or gradient change. Constant
// curvature therefore retains the same pressure balance. This bounds branch
// selection sensitivity, not the physical error of an unresolved interface.
// report is assigned only on success. No independent scalar reconstruction.
HydroUniformLoad buildUniformConservativeCapillaryLoad(const QuadtreeMesh& mesh,
    const NativeLeafDomain& domain,const InterfaceState& state,double surfaceTension,
    double projectionHalfWidth,ConservativeCapillaryReport* report=nullptr,
    std::shared_ptr<const ConservativeTractionGeometry>* tractionGeometry=nullptr);

// Complete mixed XY/RZ native leaves. Uses the supplied physical pressure
// operator and its actual dual volumes for both gradient and diagnostics.
// The same curvature/boundary acceptance applies; no fallback across an
// unsupported coarse/fine curvature stencil. This alone does not enable AMR NS.
HydroUniformLoad buildNativeConservativeCapillaryLoad(const QuadtreeMesh& mesh,
    const NativeLeafDomain& domain,const InterfaceState& state,
    const HydroConservativePressure& pressure,double surfaceTension,
    double projectionHalfWidth,ConservativeCapillaryReport* report=nullptr,
    std::shared_ptr<const ConservativeTractionGeometry>* tractionGeometry=nullptr);
}
