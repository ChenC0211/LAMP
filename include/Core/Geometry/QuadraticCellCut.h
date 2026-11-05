#pragma once
#include "Core/Geometry/PlanarCellCut.h"
#include <optional>
#include <memory>

namespace LAMP::Geometry::Detail {struct QuadraticGraphExact;struct QuadraticParametricData;}

namespace LAMP::Geometry::QuadraticCellCut
{
using Rectangle=PlanarCellCut::Rectangle;
using Point=PlanarCellCut::Point;
// Cell-centred graph n = offset + slope*t + second*t*t/2.
// axis is the normal coordinate (0=x/r, 1=y/z); sign*(n-H)<=0 is liquid.
// Optional frameNormal replaces the coordinate axes by its unit normal n and
// t=(-n.y,n.x). The same sign*(n-H)<=0 convention applies. Coefficients remain
// physical derivatives about the cell centre; RZ still measures 2*pi*r dA.
struct Shape {int axis=1,sign=1;long double slope=0,second=0;
    std::optional<Vector2D> frameNormal;};
struct Basis {Point tangent,normal;};
Basis basis(const Shape&);
struct Reconstruction {Shape shape;long double offset=0,volume=0,residual=0;
    std::shared_ptr<const Detail::QuadraticGraphExact> exact;};
// Exact same-quadratic numerical encoding. Frame/side and the ORIGINAL shared
// physical ports form c0/c1 before any native rounding; c2=eta/max(cell widths).
Reconstruction fromFixedPorts(const Rectangle&,Shape frameAndPhase,Point first,Point last,long double eta);
// One exact polynomial through THREE actual open-face ports in the declared
// fixed frame. This is a geometric chart, not a stock/one-arc candidate.
// All three original port events are retained; no tolerance merges tangent
// parameters and no rounded eta is used to reconstruct the polynomial.
struct ThreePortGraphResult {
    Reconstruction graph;
    // [input port][physical x/y coordinate][normalized q0/h,q1,q2*h].
    // Each derivative has units 1/length; h=max(actual cell widths).
    std::array<std::array<std::array<long double,3>,2>,3> normalizedCoefficientJacobian{};
};
ThreePortGraphResult fromThreePorts(const Rectangle&,Shape frameAndPhase,const std::array<Point,3>& ports);
// Divided incidence on one actual open physical face. knownRoot MUST be
// exactly on the encoded curve; fourthNode need only be on that same face.
// Equal input points are allowed: the residual has NO coincident-root factor.
struct SameFaceDividedIncidence {
    long double residual=0;Point companion{};int face=-1; // 0L,1R,2B,3T
    // wrt normalized (q0/h,q1,q2*h), dimensionless.
    std::array<long double,3> normalizedCoefficientJacobian{};
    // [knownRoot/fourthNode][physical x/y], units 1/length.
    std::array<std::array<long double,2>,2> pointJacobian{};
};
// residual=(fourth tangential coordinate - Vieta companion coordinate)/h.
// The companion is a same-curve numerical view on the supporting face line,
// possibly outside its finite segment; this does not certify an arc or stock.
// Jacobians represent the extension valid on the KNOWN-ROOT identity: total
// variations must retain that identity and both points' original physical
// face. In particular, zero coefficient-0 derivative is NOT an independent
// fixed-point root derivative when changing q0. No new geometry owner exists.
SameFaceDividedIncidence sameFaceDividedIncidence(const Rectangle&,const Reconstruction& currentExactGraph,
    Point knownRoot,Point fourthNode);
// One original fixed-port chart, with volume eliminated analytically. The
// native eta is a diagnostic; graph retains the unrounded algebraic solution.
struct FixedPortVolumeResult {
    Reconstruction graph;
    long double eta=0,chordPhaseVolume=0,volumeResidual=0,ownPhaseResidual=0,volumeDerivativeEta=0;
    bool quadratic=false;unsigned realRoots=0,branchCandidates=0;
    // Same selected-phase volume as graph.volume, differentiated with respect
    // to normalized (c0/h,c1,c2*h); reused from the final exactPolygon measure.
    std::array<long double,3> selectedPhaseDerivative{};
};
// targetPhaseVolume belongs to frameAndPhase.sign. branchEta declares the
// original monotone V(eta) half-branch, never a nearest-root/fallback policy.
// Covers ordinary fixed ports only; the regular-axis chart is separate.
FixedPortVolumeResult fixedPortVolumeGraph(const Rectangle&,Shape frameAndPhase,
    Point first,Point last,long double targetPhaseVolume,bool axisymmetric,long double branchEta);
// Second derivatives, NOT coefficients of lambda squared. Ports move as
// p_i(lambda)=p_i+lambda*v_i along their actual original cell faces.
struct Jet2 {long double value=0,first=0,second=0;};
// One physical-coordinate stationary point, eligible only while its OTHER
// coordinate lies in the actual face's tangent interval. side[] indexes
// lower/upper, with physical side ids 0Left/1Right/2Bottom/3Top.
struct StationaryFaceBoundary {
    long double parameter=0;Point position{};bool inOriginalArc=false;
    std::array<int,2> side{},branchSign{};
    std::array<long double,2> slack{}; // fixed branchSign * inwardGap / h
    std::array<std::array<long double,3>,2> normalizedCoeffGradient{}; // wrt(c0/h,q1,c2*h)
};
// Pure current-geometry events on any nonempty set of actual arcs.
// inOriginalArc means membership in their union of closed parameter intervals.
// No target-volume requirement: the existing joint phase-I can consume these
// same branch bounds before volume feasibility. The fixed-volume jet below
// independently retains its original ONE-arc requirement.
std::array<std::optional<StationaryFaceBoundary>,2> stationaryFaceBoundaries(const Rectangle&,const Reconstruction& currentExactGraph);
struct StationaryFaceBoundaryJet {StationaryFaceBoundary current;std::array<Jet2,2> slack{};};
struct FixedPortVolumeBoundaryJet {
    Jet2 eta,slope,span; // c2*h, q1, abs(t1-t0)/h
    std::array<Jet2,4> cornerImplicit{}; // LL,LR,UR,UL; selected signed implicit / h
    std::array<std::array<Jet2,2>,2> portTangent{}; // [original input port][physical x/y] dX/dt
    std::array<std::optional<StationaryFaceBoundaryJet>,2> stationaryFace{}; // physical coordinate x/y
};
// A numerical second-order model on the CURRENT exact fixed-volume branch.
// No candidate/owner/acceptance: every actual trial retains the original gates.
// Physical-axis ports require the independent axis chart, not this free-port jet.
FixedPortVolumeBoundaryJet fixedPortVolumeBoundaryJet(const Rectangle&,const Reconstruction& currentExactGraph,
    Point first,Point last,Point firstVelocity,Point lastVelocity,long double targetPhaseVolume,bool axisymmetric);
Reconstruction fromAxisPorts(const Rectangle&,int liquidSign,Point axis,Point other);
// RZ regular-axis chart, with a fixed actual r=0/open axis point and ONE
// selected actual companion subface. phaseSign=+1 measures below z(r), -1
// measures above, exactly as Shape.sign. The two segment endpoints can be
// reversed; companionFraction follows that original direction in open (0,1).
// Only the right vertical or top/bottom horizontal cell face is supported.
// Solves a rational k or r^2 analytically, retains it in the same exact graph,
// and returns the companion from that graph's original side-root event. It
// never rounds a sqrt endpoint then rebuilds fromAxisPorts, clips an inventory,
// runs a volume Newton, or adds mesh/Interface ownership. Irrational endpoints
// use the existing bounded higher-precision Core root model, not symbolic roots.
struct AxisPortVolumeResult {
    Reconstruction graph;
    Point companion;
    long double companionFraction=0,phaseVolume=0,volumeResidual=0,ownPhaseResidual=0;
    bool horizontalCompanion=false;
};
AxisPortVolumeResult axisPortVolumeGraph(const Rectangle&,int phaseSign,Point axis,
    Point faceBegin,Point faceEnd,long double targetPhaseVolume);
Reconstruction recenter(const Rectangle& source,const Rectangle& target,const Reconstruction&);
std::shared_ptr<const Detail::QuadraticParametricData> materialCurveData(const Rectangle&,
    const Reconstruction&,long double lower,long double upper);
// Same centred phase polynomial used by interval and temporal-event queries.
long double implicitValue(const Rectangle&,const Reconstruction&,Point);
// Rotated, scaled graph n=c[0]+c[1]*t+c[2]*t^2, liquid below the graph.
// (t,n) are coordinates in the tangent=(-ny,nx)/normal frame, divided
// by scale. The supplied normal is normalized internally. The derivatives
// are with respect to these three dimensionless coefficients, including
// the actual XY/RZ cell measure and moving edge intersections.
struct GraphMeasure {long double volume=0;std::array<long double,3> derivative{};};
GraphMeasure graphMeasure(const Rectangle&,Point origin,Vector2D normal,long double scale,
    const std::array<long double,3>& coefficients,bool axisymmetric);
// Actual sign-side volume and derivatives wrt normalized (c0/h,c1,c2*h).
GraphMeasure graphMeasure(const Rectangle&,const Reconstruction&,bool axisymmetric);
// Closed polynomial moments after floating quadratic root splitting; this is
// not an interval certificate. The inverse enforces the original PLIC absolute
// volume tolerance and never changes the requested fraction.
Reconstruction reconstruct(const Rectangle&,Shape,long double fraction,bool axisymmetric);
// Project an entirely interior single graph onto a parabola preserving its
// two endpoint heights and the supplied XY/RZ liquid volume. Heights are
// absolute physical normal coordinates at the tangent bounds. Returns no
// graph when these data admit no contained regular parabola; no clipping.
std::optional<Reconstruction> projectEndpointVolume(const Rectangle&,int axis,int sign,
    long double lowerHeight,long double upperHeight,long double liquidVolume,bool axisymmetric);
struct StripPhaseMeasure {long double total=0,liquid=0,gas=0;
    // Only an explicitly requested linear-sweep donor-phase comparison sets
    // these bits (1=liquid, 2=gas). Zero-measure donor phases do not set a bit.
    unsigned coveredPhases=0;};
// A donor-interior axis-aligned strip, retaining near plane and
// positive physical depth separately. Its far plane is near +/- depth as an
// exact coordinate sum; never construct and subtract rounded far coordinates.
// The original axis graph (no frameNormal) is integrated in a unit strip chart.
// Both phase measures use true shared event intervals and independent positive
// bounds; total is the separate factored XY/RZ rectangle metric, never a phase
// normalization or an inventory value. Graph/root approximation scope matches
// the existing closed-polynomial measure; no alpha, clipping or owner is added.
StripPhaseMeasure axisStripMeasure(const Rectangle&,const Reconstruction&,int normalAxis,
    Point nearBegin,Point nearEnd,long double depth,bool positive,bool axisymmetric);
// Original frame graph counterpart. Exact near/depth corners and original
// frame coefficients retain the actual convex strip in (t,n), with reduced
// exact integer ratios per slab and true physical Jacobian/radius. The private
// integer backend has no binary64 expansion exponent window or public ABI.
// Shares the original floating quadratic event-root scope, not interval roots.
StripPhaseMeasure framedStripMeasure(const Rectangle&,const Reconstruction&,int normalAxis,
    Point nearBegin,Point nearEnd,long double depth,bool positive,bool axisymmetric);
// Affine normal sweep of the ORIGINAL near-face segment. Its signed
// displacement is d(s)=meanDepth+halfVariationDepth*(2*s-1), s in [0,1].
// Inputs remain separate exact binary values: neither endpoint depth nor a
// spatial zero/back-edge is rounded in world coordinates. positive selects
// d>0, negative selects d<0; in either case the donor back-edge is near-d.
// Optional actual query points restrict this ORIGINAL law by an exact tangent
// coordinate ratio; they must lie on the same near segment, in either order.
// Each selected part must stay inside the donor. Exact graph only; no new
// geometry, phase inventory, velocity or time owner. XY/RZ use the same true
// polygon measure. In RZ meanDepth is the PARAMETER-midpoint value, not a
// general radius-weighted face mean. Caller owns that trace conversion.
// Coverage is opt-in: it costs one additional whole-donor exact measure.
StripPhaseMeasure linearSignedNormalSweepMeasure(const Rectangle&,const Reconstruction&,int normalAxis,
    Point nearBegin,Point nearEnd,long double meanDepth,long double halfVariationDepth,
    bool positive,bool axisymmetric,std::optional<std::array<Point,2>> query=std::nullopt,
    bool identifyCoveredPhases=false);
// Nonnegative endpoint-depth form. Negative depths reject: spatial sign
// changes belong to the signed, factored interface above. Zero endpoints and
// zero total sweep are valid, without a fake phase or rounded zero event.
StripPhaseMeasure linearNormalSweepMeasure(const Rectangle&,const Reconstruction&,int normalAxis,
    Point nearBegin,Point nearEnd,long double depthBegin,long double depthEnd,
    bool positive,bool axisymmetric,std::optional<std::array<Point,2>> query=std::nullopt,
    bool identifyCoveredPhases=false);
PlanarCellCut::Overlap overlapInFrame(const Rectangle&,const std::vector<Point>& polygon,Point origin,
    const Reconstruction&,bool axisymmetric);
// Liquid portions of a retained chord interval, in its original [0,1]
// parameter. Caller has already clipped that interval to the donor rectangle.
std::vector<std::array<long double,2>> intervals(const Rectangle&,const Reconstruction&,
    Point begin,Point end,long double lower=0,long double upper=1);
}

