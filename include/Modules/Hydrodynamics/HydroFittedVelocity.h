#pragma once

#include <memory>
#include <vector>
#include "Core/Base/Vector2D.h"
#include "Modules/Hydrodynamics/HydroFittedFaceFlux.h"

class QuadtreeMesh;
class HydroFittedDomain;
class NativeLeafDomain;
class LevelSetMappedSegment;
template<typename T> class VEdgeField;
template<typename T> class HEdgeField;

// Immutable velocity snapshot on the actual active leaf domain. The original
// constructor retains its same-level behavior described here; RZFluxCompatible
// also reconstructs shared vertex data from native NS fields on rectangular
// 2:1 leaves. The explicit VertexTrace overload accepts caller-supplied data.
// Solid-only field payload is never read. Full-fluid native stencils retain
// staggered bilinear interpolation; missing wall strips use the fitted Hydro
// no-slip, symmetry or pressure zero-shear trace, not ghost-cell values.
// Every cell center, native midpoint and vertex coordinate is retained; actual
// leaf coordinates need not form a bitwise Cartesian product. Shared sides
// pass through their own native midpoint, and inverse bilinear coordinates use
// those actual local points. Strictly rectangular native stencils preserve the
// original bilinear expression. Only the closed union of active cells is
// supported; there is no tolerance-based crossing of a wall or point snapping.
class HydroFittedVelocitySampler
{
public:
    enum class Reconstruction { StaggeredBilinear, RZFluxCompatible };
    using FaceChart = HydroFittedFaceFlux::Chart;
    HydroFittedVelocitySampler(const QuadtreeMesh& mesh, const HydroFittedDomain& domain,
        const VEdgeField<double>& ux, const HEdgeField<double>& uy,
        Reconstruction reconstruction = Reconstruction::StaggeredBilinear);

    struct VertexTrace
    {
        Vector2D velocity;
        double axialDerivative=0.0; // d(u_z)/dz in caller units, shared per vertex.
    };
    // Explicit trace-driven RZ reconstruction on actual rectangular 2:1 leaves.
    // The caller supplies full native vertex indexing, with finite shared jets
    // on active face endpoints. Wall/axis speeds must obey the domain exactly.
    // Uses the existing face primitive and piecewise Hermite streamfunction
    // extension; does not infer a pressure-boundary/solid/contact trace model.
    // All inputs are retained by value. Same-level original behavior is unchanged.
    HydroFittedVelocitySampler(const QuadtreeMesh& mesh, const HydroFittedDomain& domain,
        const VEdgeField<double>& ux, const HEdgeField<double>& uy,
        const std::vector<VertexTrace>& vertexTraces);
    // Complete native XY/RZ domain with caller-owned physical boundary traces.
    // Reuses the same face primitive/streamfunction extension, without imposing
    // fitted pressure/no-slip BCs on a complete-domain normal-velocity inlet.
    // Native face averages remain authoritative. No implicit wall law is added.
    HydroFittedVelocitySampler(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
        const VEdgeField<double>& ux, const HEdgeField<double>& uy,
        const std::vector<VertexTrace>& vertexTraces);
    ~HydroFittedVelocitySampler();
    HydroFittedVelocitySampler(HydroFittedVelocitySampler&&) noexcept;
    HydroFittedVelocitySampler& operator=(HydroFittedVelocitySampler&&) noexcept;
    HydroFittedVelocitySampler(const HydroFittedVelocitySampler&) = delete;
    HydroFittedVelocitySampler& operator=(const HydroFittedVelocitySampler&) = delete;

    // Fitted-domain constructors keep true no-slip vertices exactly stationary.
    // The original constructor supplies the following inferred boundary traces:
    // At a true no-slip vertex (including the fixed lip) both components are
    // exactly zero. RZ axis u_r is exactly zero and the missing u_z strip is
    // even in r, using the first two actual native axial samples. A pressure
    // corner lacking the fitted inward zero-shear stencil is rejected.
    // Read-only after construction, including native leaf lookup and flux
    // reconstruction. Independent concurrent queries share no mutable work.
    Vector2D sample(Vector2D point) const;

    struct FirstJet
    {
        Vector2D value, dx, dy; // dx/dy differentiate both velocity components.
        int nativeCell = -1;
        bool oneSided = false; // Cell side or piecewise face-trace seam.
    };
    // Analytic derivative of the SAME piecewise streamfunction reconstruction
    // used by the explicit VertexTrace constructors, the mixed-leaf flux
    // constructor and completeVelocity().
    // At a seam, returns the same owning branch as sample(), marked oneSided;
    // does not claim a globally C1 velocity. RZ includes the regular axis limit.
    // Legacy bilinear/same-level Piola constructors explicitly reject for now.
    FirstJet firstJet(Vector2D point) const;
    // Same source in volume coordinates: XY is unchanged; RZ uses s=r^2/2,
    // value=(r*u_r,u_z), dx=d/ds and dy=d/dz. Input remains physical (r,z).
    // The axis limit is analytic, without dividing a sampled derivative by r.
    FirstJet volumeFirstJet(Vector2D physicalPoint) const;

    // In the original constructor, RZFluxCompatible preserves each native flux using a
    // weighted Piola map of the same actual four subquadrilaterals. A local
    // streamfunction correction supplies shared vertex traces and whole-wall
    // no-slip, without changing any face integral or erasing the FV residual.
    // Native values are face averages here, not prescribed midpoint samples.
    // Two perpendicular no-slip walls in one cell are currently unsupported.
    // This does not make pointwise level-set advection volume-conservative.
    // The explicit overload uses the same face fluxes on piecewise cell sides
    // and a Hermite streamfunction extension; it also retains the FV residual.
    // Analytic physical RZ divergence of either flux reconstruction (bilinear rejects).
    double divergence(Vector2D point) const;

    // RZFluxCompatible only. Native rows use the fitted domain's full V+H
    // ordering; a solid-only row is not a valid face. These retained actual
    // endpoints and native midpoint define the two distinct affine charts.
    FaceChart faceChart(int nativeRow) const;
    // half=0: begin -> midpoint; half=1: midpoint -> end. The two parameters
    // belong to THAT half's [0,1] chart, may decrease, and are never converted
    // to rounded global root points. Integrates the existing reconstruction's
    // normal flux, with positive normal toward +r (vertical) or +z (horizontal),
    // including the actual half-edge slope. Reversing parameters negates the
    // result. Units are caller code velocity*RZ area. No extra measure factor
    // is needed. Two full halves sum algebraically to the retained native Q.
    long double integratedNormalFlux(int nativeRow, int half, double tBegin, double tEnd) const;
    // Anchored variant: half parameter t=anchor0+(anchor1-anchor0)*s. Both
    // anchors and both s endpoints lie in [0,1]; either pair may reverse.
    // Retains the factored interval even when its composed t endpoints round
    // to the same number. The caller must establish that these anchors denote
    // its actual geometric subsegment; no rounded root remapping occurs here.
    long double integratedNormalFlux(int nativeRow, int half, double anchor0, double anchor1,
        double sBegin, double sEnd) const;

    // Negative-phase flux of a mapped segment whose original begin/end must
    // EXACTLY match this actual half-face's positive-direction chart. Near-end
    // offsets and positive widths are integrated separately; fromEnd does not
    // reverse the +r/+z normal. No extra RZ/embedding Jacobian is applied.
    // Floating polynomial integration, not an interval enclosure.
    long double integratedNormalFlux(int nativeRow, int half, const LevelSetMappedSegment& segment) const;

    // Same retained polynomial as integratedNormalFlux; includes the actual
    // half-face/RZ metric, not normal speed and not a second reconstruction.
    long double normalFluxDensity(int nativeRow, int half, bool fromEnd, long double offset) const;
    Quantum::Geometry::UnitPolynomialRoots normalFluxRoots(int nativeRow,int half) const;
    std::vector<Quantum::Geometry::UnitPolynomialTerm> normalFluxPolynomial(int nativeRow,int half) const;

private:
    struct Impl;
    std::unique_ptr<const Impl> impl_;
};
