#pragma once

#include "Modules/Hydrodynamics/HydroPressureLoad.h"
#include "Modules/Hydrodynamics/HydroRegularizedColor.h"
#include "Modules/Interface/HermiteGalsAmrClosure.h"
#include "Modules/Interface/HermiteGalsPublication.h"
#include <memory>
#include <array>
#include <utility>
#include <optional>
#include <optional>

namespace LAMP::Hydrodynamics
{
class HydroAdvectorReconstruction;
class HydroAdvectorBoundary;
class HydroSpatialAdvector;
// Declared geometric contact law only: H(z) is even across the real south
// plane (90 degrees, freely moving contact line). Native rows are WHOLE actual
// exterior subfaces. No pinned radius, ghost inventory or transport BC.
class HydroVofSouthEven90Boundary
{
public:
    static HydroVofSouthEven90Boundary prepare(const QuadtreeMesh&,
        const Interface::GeometricVofPublication&,const std::vector<int>& actualSouthNativeRows);
    const Interface::GeometricVofPublication& source() const;
    const std::vector<int>& nativeRows() const;
    const std::vector<std::array<double,2>>& intervals() const;
    double coordinate() const;
    void validateFor(const QuadtreeMesh&,const Interface::GeometricVofPublication&) const;
    bool sameIdentity(const HydroVofSouthEven90Boundary&) const noexcept;
private:
    struct Data;
    explicit HydroVofSouthEven90Boundary(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
};
enum class HydroCsfForceMethod {MotionAdjoint,HeightFunction,HeightFunctionIncrement};
struct HydroCsfCurvatureConfig
{
    HydroCsfForceMethod forceMethod=HydroCsfForceMethod::MotionAdjoint;
    std::size_t maxCandidateQueriesPerFace=262144;
    std::size_t maxVisibilityCellQueriesPerFace=262144;
    std::size_t maxSurfaceQueriesPerFace=262144,maxGeometryQueriesPerFace=32768;
    std::size_t maxShadowWallQueriesPerFace=262144;
    int maxQuadratureDepth=24;
    double relativeQuadratureTolerance=1e-13;
    double relativeGradientRoundoff=256*std::numeric_limits<double>::epsilon(); // Audit only; never filters S.
};
enum class HydroCsfDiscretization {RegularizedColor,GeometricPhaseFraction,GeometricAdjoint,GeometricProjectedAdjoint,GeometricHeightFunction,GeometricHeightFunctionIncrement};
enum class HydroCsfSupportStatus {Projected,ZeroGradient,PureBoundaryReaction,SigmaZero};
// The unresolved value is a declared force closure, never a geometric/HF
// curvature estimate. Its independent phase inventories remain unchanged.
enum class HydroCsfCellCurvatureMethod {Unavailable,HeightFunction,SegmentFit,UnderresolvedIsolatedFragmentZero,MixedHeightFit=4,TwoSegmentCircleFit=5};
struct HydroCsfFaceSupportReport
{
    Vector2D center{0,0};double radius=0;
    int pressureComponent=-1;
    std::size_t candidateQueries=0,supportCount=0,visibilityRejected=0,visibilityCellQueries=0;
    long double totalWeight=0; // Local mean: one/two cell weights; adjoint: integral abs(n.B_u or Q_Gamma,u).
    double curvatureMin=0,curvatureMax=0,maxSupportDistance=0;
    double gradientRoundoffScale=0,colorGradientErrorEstimate=0;
    std::size_t surfaceQueries=0,geometryQueries=0,shadowWallQueries=0,acceptedPanels=0;
    std::size_t trialPositiveGraphsWithoutFineSupport=0; // Audit only, never a force filter.
    int deepestPanel=0;
    long double weightErrorEstimate=0,anchoredMomentErrorEstimate=0;
    long double normalMomentErrorEstimate=0,tractionMomentErrorEstimate=0; // Both adjoint kinds: actual native moments
    double curvatureErrorEstimate=0,curvatureScale=0,quadratureAnchor=0;
    bool quadratureAnchorDefined=false;
    bool curvatureErrorEstimateAvailable=true; // Height-function truncation error is not estimated by a quadrature difference.
    bool curvatureDefined=false,boundary=false,nearBoundary=false,nearAxis=false,crossesAmr=false,subResolutionGradient=false;
    bool coincidentPhaseHeight=false; // Curvature evaluated at a true pure-L/G shared face.
    bool underresolvedFragmentClosure=false; // At least one actual adjacent force coefficient is the qualified zero closure.
    HydroCsfSupportStatus status=HydroCsfSupportStatus::ZeroGradient;
};
struct HydroCsfQuadraturePanel
{
    std::size_t component=0,graph=0;
    double begin=0,end=1;
    int depth=0;
    std::size_t coarseSubdivisions=1; // Actual accepted rule uses 2*coarseSubdivisions GL8 subarcs.
    bool wholeGraphEstimate=false; // Zero tails carry no independent absence certificate.
    long double weightErrorEstimate=0,anchoredMomentErrorEstimate=0;
};
// A successful read-only receipt from one actual bounded real-column call.
// Only masked columns use the exact canonical finite-inventory RN law;
// unmasked physical heights retain the original Kahan calculation.
struct HydroCsfFiniteHeightCertificate
{
    int heightAxis=0,outward=1;
    Vector2D center{0,0};
    double normalSpacing=0,tangentSpacing=0;
    std::array<double,3> tangentCenters{};
    long double origin=0,nativeInverseSpacing=0;
    std::array<long double,3> physicalHeights{};
    std::array<bool,3> certified{};
    std::size_t queryCount=0;
};
struct HydroCsfCurvatureReport
{
    std::size_t sourceNodes=0,selectedNodes=0,candidateQueries=0,visibilityRejected=0;
    std::size_t projectedFaces=0,pureBoundaryRows=0,nearBoundaryFaces=0,nearAxisFaces=0,amrFaces=0;
    std::size_t heightFunctionCells=0,heightFunctionColumns=0,heightFunctionCellQueries=0,heightFunctionAlternateDirections=0;
    std::size_t heightFunctionXCells=0,heightFunctionYCells=0,heightFunctionPartialCertificates=0;
    std::size_t heightFunctionExtendedCells=0,heightFunctionExtendedSlices=0;
    std::size_t heightFunctionRefinedCells=0,heightFunctionRefinedSlices=0;
    // HF Cells/X/Y count successful heights only; LocalPlic mixed cells = HF + Fit + Unresolved.
    std::size_t heightFunctionFitCells=0,heightFunctionFitSegments=0;
    std::size_t unresolvedIsolatedFragmentCells=0; // Separate from successful HF/fit estimates.
    std::size_t unresolvedVertexStarFragmentCells=0; // Subset of unresolved total above, never add twice.
    std::size_t heightFunctionCoarseAlignedCharts=0,heightFunctionCoarseAlignedCells=0;
    std::size_t southEven90BoundaryCells=0,southEven90ReflectedSlices=0,pureJumpHeightFaces=0;
    std::size_t pinnedCornerHeightCells=0,pinnedCornerHeightFaces=0,pinnedCornerHeightColumns=0;
    // Unique final-native-height certificates only. Source minor geometry
    // remains unresolved and both independent inventories remain unchanged.
    std::size_t finiteInventoryHeightColumns=0,finiteInventoryHeightQueries=0;
    std::vector<HydroCsfFiniteHeightCertificate> finiteInventoryHeightCertificates;
    std::size_t incrementNodes=0,incrementRatePoints=0,incrementLoadDofs=0;
    bool projectionSkippedSigmaZero=false;
};
// Qualified column jet: H is positive liquid thickness along outward;
// the physical graph is origin+outward*H. Derivatives use the physical
// tangent coordinate. A true r=0 axial-height even jet alone has axisRegular.
struct HydroCsfHeightJet
{
    long double height=0,first=0,second=0,origin=0,tangentCoordinate=0;
    int heightAxis=0,outward=1;
    bool axisRegular=false,axisReflected=false;
    // Actual selected RZ chart, possibly an ancestor-centred material frame.
    // Its curvature is not a fine-cell-centre second-order certificate.
    long double normalSpacing=0,tangentSpacing=0;
    bool coarseAligned=false;
    bool southEven90Reflected=false; // Declared S reflection, distinct from the true RZ axis.
    bool pinnedCornerConstrained=false; // Same-source Interface certificate plus three real bins.
    bool finiteInventoryHeightSupport=false; // Height interval certificate, never a recovered minor plane.
    int pinnedCornerBin=-1; // Actual first/second tangent bin, otherwise unavailable.
};
// Read-only balanced CSF load on the retained authoritative C=-D_I^T.
// The output is the ACTUAL momentum RHS of A*u+C*p=rhs: add it directly.
// Native boundary rows are retained; the existing saddle handles wall reactions.
class HydroBalancedCsfLoad
{
public:
    static HydroBalancedCsfLoad constantCurvature(const QuadtreeMesh&,HydroPressureLoad,
        HydroRegularizedColor,double surfaceTension,double curvature);
    static HydroBalancedCsfLoad capillary(const QuadtreeMesh&,HydroPressureLoad,
        HydroRegularizedColor,double surfaceTension,HydroCsfCurvatureConfig={});
    // C*alpha of this SAME geometric material owner. Curvature is a true
    // cell-local positive surface average, then the arithmetic mean of the
    // one/two adjacent interface cells. No force band or wider extension.
    static HydroBalancedCsfLoad geometricCapillary(const QuadtreeMesh&,HydroPressureLoad,
        HydroPhaseGeometry,double surfaceTension,HydroCsfCurvatureConfig={});
    // Three-column second-order height curvature of the SAME measured phase
    // volumes. Complete uniform XY / SingleQuartic only; missing single-branch
    // pure-phase support rejects. S_u=sigma*kappa_HF,f*(C*alpha)_u. No surface
    // quadrature, raw curvature fallback or exact motion-adjoint power claim.
    static HydroBalancedCsfLoad heightFunctionCapillary(const QuadtreeMesh&,HydroPressureLoad,
        HydroPhaseGeometry,double surfaceTension,HydroCsfCurvatureConfig={});
    // Actual LocalPlic owner; exact selected XY/RZ native-face 2:1 AMR.
    // A same-publication Interface pin certificate alone permits a bounded
    // one-sided cubic chart from three actual equal-width phase-area bins
    // when the ordinary first/second-bin support meets its real solid wall.
    // Pure coincident pin branches retain unavailable mixed-cell curvature.
    // RZ heights use independent ordinary meridional phase areas; force color
    // and C retain physical 2*pi*r volumes/measures. The SAME qualified jet
    // supplies meridional+hoop curvature. Only the true r=0 axis permits
    // read-only even reflection/truncation; other physical walls still reject.
    // Original 9x3 equal-width columns integrate partial actual PLIC leaf boxes.
    // Only after both original axes fail, missing pure ends may add at most
    // four common rings/end (17x3 maximum), batched over the same actual owner.
    // After both original/extended axes fail, coarse mixed cells alone may
    // use their actual PLIC segment midpoint and actual finest-leaf hx/hy for
    // a virtual 9x3 frame with the same four-ring/end bound and branch guards.
    // Reversed/unresolved branches, physical boundaries and missing ends reject;
    // After all height frames fail, a fixed 5h box may fit a quadratic to
    // the source native-face-connected, same-oriented actual PLIC segments.
    // At least three independent two-sided segment midpoints and rank three
    // are required; no neighbourhood search, curvature clipping or new owner.
    // Only after all estimates fail, an Interface-issued XY closed single-cell
    // component certificate with actual diameter<min(dx,dy) permits the named
    // unresolved zero force closure (Popinet 2009 Algorithm 7(4b)). It supplies
    // no physical curvature estimate. Connected/wall/ambiguous/RZ failures reject.
    // No FirstJet, changed face averaging, clipped heights or A owner.
    // RZ only: if native tangent strips would subdivide a coarser mixed PLIC
    // leaf, select once the actual source ancestor at the observed coarse
    // level. All27/extended strips retain that aligned frame; any further
    // coarser partial mixed support rejects instead of recursively widening.
    // Position and effective spacings are retained in the successful jet.
    static HydroBalancedCsfLoad localPlicHeightFunctionCapillary(const QuadtreeMesh&,HydroPressureLoad,
        HydroPhaseGeometry,double surfaceTension,HydroCsfCurvatureConfig={});
    // Explicit geometric boundary overload. The view names this exact phase
    // publication/full clock. Other physical walls and uncovered S intervals
    // remain unsupported. All phase inventories and force C are unchanged.
    static HydroBalancedCsfLoad localPlicHeightFunctionCapillary(const QuadtreeMesh&,HydroPressureLoad,
        HydroPhaseGeometry,double surfaceTension,const HydroVofSouthEven90Boundary&,HydroCsfCurvatureConfig={});
    const HydroVofSouthEven90Boundary* southEven90Boundary() const noexcept;
    // Pure jet arithmetic also consumed by the actual RZ HF estimator.
    // Returns {meridional,hoop}; no radius floor or curvature clipping.
    // At actual r=0 only an axial even regular jet (first==0) is accepted.
    static std::array<long double,2> axisymmetricHeightCurvatures(const HydroCsfHeightJet&);
    // Native pressure order. No value for pure cells/sigma-zero or a rotated
    // segment fit; every value is the SAME successful measured column jet.
    const std::vector<std::optional<HydroCsfHeightJet>>& heightFunctionCellJets() const;
    // HF base plus raw-minus-HF residual paired with the compatibleIncrement
    // scalar generator. Cached publication GL8 scatter; no recursive traction
    // certificate or full surface-energy identity is implied.
    static HydroBalancedCsfLoad heightFunctionIncrementCapillary(const QuadtreeMesh&,HydroPressureLoad,
        HydroPhaseGeometry,double surfaceTension,const HydroSpatialAdvector&,
        const Interface::HermiteGalsAmrClosureSensitivity* =nullptr,HydroCsfCurvatureConfig={});
    const std::vector<double>& heightFunctionIncrementReferenceLoad() const;
    const std::vector<double>& heightFunctionIncrementResidualLoad() const;
    const std::vector<double>& heightFunctionIncrementMeanCurvatureDifference() const;
    std::vector<std::pair<int,long double>> incrementNormalRateAtSample(
        const Interface::HermiteGalsSurfaceSample&) const;
    // Weak geometric CSF on the SAME native-to-motion reconstruction:
    // S = C*(sigma*kappa0*alpha) - sigma*integral((kappa-kappa0)*n.B_u ds).
    // Each native pressure component uses its actual positive surface mean
    // kappa0. All pointwise curvature variations are retained in the residual.
    // Independent tangential modes are not supported by this native-only load.
    static HydroBalancedCsfLoad adjointCapillary(const QuadtreeMesh&,HydroPressureLoad,
        HydroPhaseGeometry,double surfaceTension,const HydroAdvectorReconstruction&,
        const HydroAdvectorBoundary* =nullptr,HydroCsfCurvatureConfig={});
    // Frozen zero-step FirstJet/scalar generator, including the SAME compact
    // mixed closure and actual backward donor selected by the midpoint field.
    // Single Q4 and C2 macro use their actual complete source sensitivities;
    // macro completion introduces no transported Hessian or extra source ID.
    // Q_Gamma,u=-phiDot_u/|grad(phi)|; this is not a finite GALS Jacobian.
    // Initial scope: complete uniform nonhanging XY, homogeneous closed law.
    static HydroBalancedCsfLoad projectedMotionCapillary(const QuadtreeMesh&,HydroPressureLoad,
        HydroPhaseGeometry,double surfaceTension,const HydroSpatialAdvector&,
        const Interface::HermiteGalsAmrClosureSensitivity* =nullptr,HydroCsfCurvatureConfig={});
    // The external vector is (T_plus-T_minus)*n. Its normal contribution
    // joins capillary stress BEFORE the same pressure-compatible deposition:
    // J = sigma*kappa - n.external; physical tangent keeps its positive sign.
    // The source is an authenticated coeval trace, never a caller callback.
    static HydroBalancedCsfLoad adjointInterfaceLoad(const QuadtreeMesh&,HydroPressureLoad,
        HydroPhaseGeometry,double surfaceTension,const Interface::HermiteGalsSurfaceVectorView&,
        const HydroAdvectorReconstruction&,const HydroAdvectorBoundary* =nullptr,HydroCsfCurvatureConfig={});
    static HydroBalancedCsfLoad projectedMotionInterfaceLoad(const QuadtreeMesh&,HydroPressureLoad,
        HydroPhaseGeometry,double surfaceTension,const Interface::HermiteGalsSurfaceVectorView&,
        const HydroSpatialAdvector&,const Interface::HermiteGalsAmrClosureSensitivity* =nullptr,
        HydroCsfCurvatureConfig={});
    bool hasExternalSurfaceVector() const noexcept;
    const Interface::HermiteGalsSurfaceVectorView& externalSurfaceVector() const;
    const std::vector<double>& referenceNormalStresses() const; // same-component actual J0
    const std::vector<double>& tangentialRhsLoad() const; // actual +external tangent dot R_u
    const HydroSpatialAdvector& projectedMotionAdvector() const;
    const Interface::HermiteGalsAmrClosureSensitivity& mixedSensitivity() const;
    const HydroAdvectorReconstruction& reconstruction() const; // Both adjoint kinds
    const HydroAdvectorBoundary* boundary() const; // Both adjoint kinds
    // Same-component actual kappa0 at each pressure DOF. NaN denotes no
    // interface in that component (or sigma-zero, when no samples are taken).
    const std::vector<double>& referenceCurvatures() const;
    const std::vector<double>& referencePressureLoad() const; // actual C*(J0*alpha); J0=sigma*kappa0 without external
    const std::vector<double>& adjointResidualLoad() const; // actual residual traction, not a face multiplier
    // Pure positive-average arithmetic, also used by capillary(). This is NOT
    // a geometry/source override port; capillary reads only its publication.
    static double positiveCurvatureAverage(const std::vector<double>& weights,const std::vector<double>& curvatures);
    // Exact closed-segment/rectangle-union predicate. The factory supplies
    // only its actual selected leaves; this is not a support override port.
    static bool selectedSegmentVisible(Vector2D begin,Vector2D end,
        const std::vector<std::array<Vector2D,2>>& selectedBoxes);
    const HydroPressureLoad& pressureLoad() const noexcept;
    const HydroRegularizedColor& color() const; // Only RegularizedColor; other kinds reject.
    const HydroPhaseGeometry& phaseGeometry() const noexcept;
    HydroCsfDiscretization discretization() const noexcept;
    const HydroDofLayout& layout() const noexcept;
    double surfaceTension() const noexcept;
    double constantCurvature() const; // rejects projected loads, including skipped sigma=0
    bool isCapillary() const noexcept;
    const HydroCsfCurvatureConfig& curvatureConfig() const;
    const std::vector<double>& curvatureAverages() const; // Actual face multiplier; adjoint and increment kinds reject.
    const std::vector<double>& heightFunctionFaceCurvatures() const; // HF base only, including the increment kind.
    const std::vector<double>& heightFunctionCellCurvatures() const; // Physical HF/fit estimates; NaN for unresolved closure/non-interface/sigma=0.
    // LocalPlic only: actual force coefficients and their distinct methods.
    // A qualified unresolved fragment has effective 0, physical curvature NaN.
    const std::vector<double>& heightFunctionCellEffectiveCurvatures() const;
    const std::vector<HydroCsfCellCurvatureMethod>& heightFunctionCellCurvatureMethods() const;
    const std::vector<HydroCsfFaceSupportReport>& supportReports() const;
    const HydroCsfCurvatureReport& curvatureReport() const;
    const std::vector<HydroCsfQuadraturePanel>& quadraturePanels(std::size_t velocity) const;
    const std::vector<double>& colorGradient() const noexcept; // C*c or C*alpha by discretization(), integral native face load
    const std::vector<double>& momentumRhsLoad() const noexcept; // Actual positive physical force; formula depends on discretization().
    Interface::HermiteGalsIdentity geometryIdentity() const;
    Interface::GeometricVofIdentity vofIdentity() const;
    void validateVofFor(const QuadtreeMesh&,const HydroPressureLoad& actualPressure,
        Interface::GeometricVofIdentity actualGeometry,Interface::GeometricVofTimeLayer actualClock) const;
    LAMP::TimeLayer timeLayer() const;
    void validateFor(const QuadtreeMesh&,const HydroPressureLoad& actualPressure,
        Interface::HermiteGalsIdentity actualGeometry,LAMP::TimeLayer actualClock) const;
    // Adjoint load must retain the SAME preparation and boundary law.
    void validateFor(const QuadtreeMesh&,const HydroPressureLoad&,
        Interface::HermiteGalsIdentity,LAMP::TimeLayer,
        const HydroAdvectorReconstruction&,const HydroAdvectorBoundary*) const;
    // Projected load additionally binds the actual direction field/full clock
    // and the authenticated mixed-closure constraints.
    void validateFor(const QuadtreeMesh&,const HydroPressureLoad&,
        Interface::HermiteGalsIdentity,LAMP::TimeLayer,
        const HydroSpatialAdvector&,const Interface::HermiteGalsAmrClosureSensitivity* =nullptr) const;
    // External owner validation supplements the matching reconstruction or
    // projected-motion overload; a new equal-valued trace is a different owner.
    void validateFor(const QuadtreeMesh&,const HydroPressureLoad&,
        Interface::HermiteGalsIdentity,LAMP::TimeLayer,
        const Interface::HermiteGalsSurfaceVectorView&) const;
    // Also requires the SAME immutable color owner, including its epsilon.
    void validateFor(const QuadtreeMesh&,const HydroPressureLoad&,const HydroRegularizedColor& actualColor) const;
private:
    static HydroBalancedCsfLoad localPlicHeightFunctionCapillaryImpl(const QuadtreeMesh&,HydroPressureLoad,
        HydroPhaseGeometry,double,HydroCsfCurvatureConfig,const HydroVofSouthEven90Boundary*);
    static HydroBalancedCsfLoad weakCapillary(const QuadtreeMesh&,HydroPressureLoad,
        HydroPhaseGeometry,double,const HydroAdvectorReconstruction&,const HydroAdvectorBoundary*,
        HydroCsfCurvatureConfig,const HydroSpatialAdvector*,const Interface::HermiteGalsAmrClosureSensitivity*,
        const Interface::HermiteGalsSurfaceVectorView* =nullptr);
    struct Data;
    explicit HydroBalancedCsfLoad(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
};
}
