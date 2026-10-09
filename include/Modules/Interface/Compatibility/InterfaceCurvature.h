#pragma once
#include "Modules/Interface/Compatibility/InterfaceHeightCurvature.h"
#include "Core/Geometry/LevelSetInterfacePanels.h"
#include <cstdint>
#include <string>

// Curvature of the one published conservative interface. Height functions are
// primary; a local implicit isotropic quadratic with cubic/quartic tangent
// corrections samples only ordered neighbors of the queried graph. Both
// available sides are retained and clipped native arcs use positive ds weights.
// This borrows the same immutable mesh/domain/state as InterfaceHeightCurvature.
// Neither fit moves the zero set or changes the conservative phase inventory.
// Mixed meshes use the queried leaf's scale and require same-level selected
// arc support. A cross-level fallback fit requests refinement, not a silently
// rescaled derivative from remote finest leaves.
// Fixed native graph midpoints and physical endpoints supply pointwise knots;
// all other queries interpolate by ds within one connected native arc. Normal
// and the two curvature parts are interpolated quantities, not the jet of a
// single blended polynomial. RZ axis endpoints keep their original even limit.
// At an event-model boundary the ordinary graph uses its own pointwise endpoint
// as a fixed knot; interpolation never borrows the model's curvature.
// Scalar roughness is regularized along a resolved same-scale branch. Open
// declared-boundary endpoints keep the original one-sided curvature values;
// actual RZ axis and declared reflection-section ends receive even reflection.
// The azimuthal axis limit remains specific to r=0. No boundary normal or
// contact-angle condition is inferred from scalar regularization.
class InterfaceCurvature
{
public:
    enum class Regularization { None, ArcThirdDerivative };
    static constexpr Regularization defaultRegularization=Regularization::ArcThirdDerivative;
    enum class RegularizationStatus { NotRequested, Applied, UnsupportedGeometry,
        IncompleteCurvature, MixedScale, UnresolvedTurn, EventPartition };
    enum class Method { Unavailable, HeightFunction, LocalFit, Mixed, EventModel };
    enum class FitStatus { NotAttempted, Valid, InsufficientArc, PoorOrientation,
        RefinementRequired, SingularFit, GeometryMismatch, UnsupportedEndpoint };
    enum class RefinementCause { None, OpposedArc, CurvatureRadius, MixedResolution };
    enum class OrientationCause { None, SupportTurn, FittedFold, NormalMismatch };
    enum class LocalModel { IsotropicQuadratic, TangentQuartic, ArcMeasureQuartic,
        SymmetricArcQuartic, ClippedArcQuartic };
    static constexpr LocalModel defaultLocalModel=LocalModel::ClippedArcQuartic;
    struct Fit
    {
        FitStatus status=FitStatus::NotAttempted;
        OrientationCause orientationCause=OrientationCause::None;
        Vector2D point,normal;
        double meridional=std::numeric_limits<double>::quiet_NaN();
        double azimuthal=std::numeric_limits<double>::quiet_NaN();
        double curvature=std::numeric_limits<double>::quiet_NaN(); // +div(n)
        double scaledRankRatio=0,scaledResidual=0,scaledPositionError=0;
        double negativeSupport=0,positiveSupport=0;
        double maximumArcCells=3.5; // Requested one-sided sampled arc bound, in local h.
        // Preserve the evidence behind a refinement request without publishing
        // the rejected curvature as a valid force. The opposed distance is the
        // triggering endpoint/midpoint sample, NOT a minimum-distance bound.
        RefinementCause refinementCause=RefinementCause::None;
        double resolutionScale=0;
        double scaledMeridionalMagnitude=std::numeric_limits<double>::quiet_NaN();
        double opposedSampleDistance=std::numeric_limits<double>::quiet_NaN();
        std::size_t samples=0;
        std::vector<int> branchGraphs;
        // A failed shorter-support attempt must remain visible alongside the
        // original wide-fit rejection. Diagnostic only; never used to admit
        // a curvature or to select a different physical stencil.
        std::string rejectedShorterFit;
        bool valid() const noexcept { return status==FitStatus::Valid; }
    };
    struct Source
    {
        Method method=Method::Unavailable;
        Vector2D point,normal;
        double meridional=std::numeric_limits<double>::quiet_NaN();
        double azimuthal=std::numeric_limits<double>::quiet_NaN();
        double curvature=std::numeric_limits<double>::quiet_NaN();
        InterfaceHeightCurvature::Sample height;
        Fit fit;
        bool regularized=false; // height/fit retain the original raw evidence.
        // EventModel uses derivatives of the CURRENT published scalar within
        // producer-declared event cells, including the even RZ pole limit.
        // It is never advertised as a resolved/high-order curvature estimate.
        bool valid() const noexcept { return method!=Method::Unavailable; }
    };
    struct Sample : Source
    {
        // Fixed native graph midpoints (and physical endpoints) are curvature
        // knots. Between them interpolate the two curvatures by native arc
        // length, not derivatives of a blended height polynomial. The inherited
        // height/fit record belongs to the first knot; secondary to the second.
        // Arc regularization changes curvature before interpolation;
        // the underlying height/fit evidence and normal remain unmodified.
        Method primaryMethod=Method::Unavailable;
        Vector2D primaryPoint;
        Source secondary;
        double secondaryWeight=0;
        double secondaryCurvatureWeight=0; // Even interpolation at an RZ pole.
        bool interpolated=false;
    };
    enum class ExtensionStatus { Valid, NoProjection, OutsideBand, AmbiguousBranch, CurvatureUnavailable };
    struct Extension
    {
        ExtensionStatus status=ExtensionStatus::NoProjection;
        Vector2D query;
        std::size_t surface=std::numeric_limits<std::size_t>::max();
        std::uint64_t surfaceId=0;
        LevelSetInterfacePanels::LocalProjection projection; // Raw canonical geometry, not fitted curvature.
        Sample sample;
        double distance=std::numeric_limits<double>::quiet_NaN();
        double signedDistance=std::numeric_limits<double>::quiet_NaN();
        double competingDistance=std::numeric_limits<double>::infinity();
        std::size_t competingSurface=std::numeric_limits<std::size_t>::max();
        LevelSetInterfacePanels::LocalProjection competingProjection;
        // Only extendConservativeFace can admit an out-of-band inventory tail.
        // Both owners must be the same geometric pure phase and their retained
        // minority inventories below the existing PLIC inverse resolution.
        bool unresolvedInventoryTail=false;
        double maxMinorityFraction=0;
        std::size_t sameChartTailCandidates=0;
        // Terminal model only: near-tied stationary points on one certified
        // native radial graph use the nearest point. Keep their curvature
        // spread as model sensitivity. Only a certified unresolved inventory
        // tail may span contiguous radial charts of a terminal-model state,
        // including the resolved/model boundary. Each retains its own
        // available provider; folds and distinct surfaces still reject.
        std::size_t eventChartCandidates=0;
        double eventCurvatureSpread=0;
        // An ambiguous pure-phase inventory tail remains AmbiguousBranch.
        // Supply ALL separated candidates in the original near-tie layer to
        // a force consumer, which may certify its actual pressure-row error.
        // Empty means no complete available-provider witness. This is not a
        // new curvature/branch acceptance and never authorizes material remap.
        struct TailAlternative {
            std::size_t surface=0;
            LevelSetInterfacePanels::LocalProjection projection;
            double curvature=0;
        };
        std::vector<TailAlternative> tailAlternatives;
        double tailCurvatureSpread=0;
        bool valid() const noexcept { return status==ExtensionStatus::Valid; }
    };
    InterfaceCurvature(const QuadtreeMesh&,const NativeLeafDomain&,const InterfaceState&,
        Regularization=defaultRegularization);
    ~InterfaceCurvature();
    InterfaceCurvature(const InterfaceCurvature&)=delete;
    InterfaceCurvature& operator=(const InterfaceCurvature&)=delete;
    Sample evaluate(std::size_t surface,int graph,double fraction=.5) const;
    // Unregularized provider at an actual knot, without interpolation. The
    // terminal-region producer uses this to validate the ordinary endpoints
    // and midpoint it will publish at a model boundary. Event tags still apply;
    // ordinary-region selection must first clear inherited model tags.
    Source evaluatePointwise(std::size_t surface,int graph,double fraction) const;
    RegularizationStatus regularizationStatus(std::size_t surface) const;
    // Same geometry/support/quality gates, with either three isotropic-circle
    // coefficients or five including tangent cubic/quartic corrections.
    // Production is height-first, then ClippedArcQuartic: reintegrate both
    // available sides without truncating the longer side to match the shorter.
    // A geometric mismatch or a support turn discovered during reintegration
    // of the ordinary 3.5h fit may retry at 2h in the production provider.
    // Fitted folds/normal mismatches still reject. Both retain >=h on each side and every
    // original branch, rank, residual, orientation and resolution gate.
    // maximumArcCells=0 selects that production policy; [2,3.5] requests one
    // fixed support for a read-only control. Other models keep their 3.5h
    // default and retain explicit read-only controls of degree, measure
    // and support. None changes geometry or conservative inventory.
    Fit localFit(std::size_t surface,int graph,double fraction=.5,
        LocalModel model=defaultLocalModel,double maximumArcCells=0) const;
    // Physical load point in the active domain; width is a physical distance.
    // Reuses retained stationary candidates on ALL surfaces and rejects near
    // ties between separated points (also within one folded surface). Curvature
    // is evaluated at the selected actual zero point, never at the offset level.
    // Search completeness is limited by Core's finite stationary inventory.
    Extension extend(Vector2D query,double halfWidth) const;
    // Preserve capillary G(C) for sub-resolution inventory transported by the
    // uniform-donor rule. This does not create a facet, delete C, or bypass
    // distinct-branch ambiguity/curvature availability. For a tail only,
    // near minima inside the SAME admitted height stencil or contiguous
    // terminal radial model chart form one chart;
    // retain its nearest point and report this continuation explicitly.
    // Tail classification comes from inventory and complete pressure support,
    // independently of whether the point is inside or outside the search band.
    // Resolved faces keep the
    // strict physical band. Uses this immutable state's actual native owners.
    Extension extendConservativeFace(int nativeFace,double halfWidth) const;
    // For a wider physical pressure row, the caller supplies ALL nonzero
    // coefficient cells plus both face owners. Tail permission requires every
    // cell in that actual support to be the same geometrically pure phase and
    // its retained minority to remain below the existing resolution threshold.
    // No coefficients/inventory are discarded and branch/height gates remain.
    // A nonzero liquidId restricts projection to that actual component. The
    // force caller must first partition C by unique geometric cell ownership
    // and apply the original pressure gradient to each component separately;
    // this is not permission to select a surface for the unsplit total G(C).
    Extension extendConservativeFace(int nativeFace,double halfWidth,const std::vector<int>& pressureSupport,
        std::uint64_t liquidId=0) const;
private:
    Extension extendImpl(Vector2D query,double halfWidth,bool permitTail,std::uint64_t liquidId=0) const;
    struct Impl;
    std::unique_ptr<Impl> impl_;
};
