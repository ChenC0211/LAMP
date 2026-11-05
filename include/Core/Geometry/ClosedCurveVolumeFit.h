#pragma once
#include "Core/Geometry/PlanarCellCut.h"
#include <vector>

namespace LAMP::Geometry::ClosedCurveVolumeFit
{
using Point=PlanarCellCut::Point;
using Rectangle=PlanarCellCut::Rectangle;
struct Constraint {Rectangle box;long double insideVolume=0;};
struct Result {
    bool converged=false;
    int controls=0,rank=0,iterations=0;
    int initialSupportRepairs=0;
    double maxInitialDisplacement=0;
    double maxNormalDisplacement=0;
    long double maxFractionResidual=0;
    std::vector<Point> polygon;
    // Preserve these vertices when observing polygon quadrature refinement:
    // fixed physical endpoints and all corners of a supplied boundary closure.
    std::vector<std::size_t> fixedPolygonVertices;
};
// Auxiliary reconstruction of ONE simple closed curve, sampled in cyclic
// order (either orientation). Uniform arc-parameter cubic B-splines have four
// local controls. Only fixed reference-normal displacements are unknown;
// tangential parameterization is not an optimization variable.
// If initial uniform resampling misses an originally mixed cell, an interior
// trace witness initializes nearby normal displacement DOFs to restore support.
// Reference curve/normals, inventory and constraints stay fixed, reflected
// parity is retained and the total normal displacement stays within the original
// quarter-cell bound. No witness or unresolved support returns failure.
// Minimize third parameter differences of the resulting vector curve under native rectangle
// volume equations. This is a branch-sized constrained fit, not a distance
// field or a topology discovery algorithm. The caller owns branch separation,
// pure-phase admission and preservation of the independent inventory.
// XY area or RZ integral 2*pi*r dA; this entry supports closed off-axis curves.
// Open axis/contact representations are supplied by the separate entries below.
// The integration polygon rejects intersections and unresolved nonadjacent
// contacts. This does not certify the continuous spline between its samples.
// Compact local derivatives, sparse interpolation and an augmented equality
// system avoid a dense branch nullspace. Sparse QR checks full constraint rank.
// Polygon quadrature and Newton/linear residuals are numerical observations,
// not interval certificates. Failure returns converged=false; no input changes.
Result fit(const std::vector<Point>& cyclicTrace,const std::vector<Constraint>&,
    double nativeSpacing,bool axisymmetric,int samplesPerControl=48);
// A single ordered RZ meridian with two distinct axis endpoints. All interior
// trace points have r>0. Reflection ties opposite normal displacements and
// imposes odd r/even z spline parity at both poles, without fixing pole z.
// Constraints are the ORIGINAL positive-radius native RZ inventories only.
// Returns the full reflected CCW auxiliary polygon; measure it only in native
// r>=0 rectangles. The mirrored half is not a physical interface component.
// Declared contact/port endpoints are not supported by this representation.
Result fitAxisToAxis(const std::vector<Point>& meridianTrace,const std::vector<Constraint>&,
    double nativeSpacing,int samplesPerControl=48);
// Open physical curve, with both endpoint coordinates fixed and no prescribed
// endpoint slope. Roughness integrates the squared third derivative of the
// clamped cubic B-spline, including its repeated endpoint knots. The fixed
// closure is NOT smoothed or assigned fit unknowns.
// fixedClosure lists intermediate boundary vertices from trace.back() to
// trace.front(), excluding those endpoints. It must close a simple polygon.
// The caller supplies the actual domain boundary and the appropriate enclosed
// phase inventory. An empty closure means a straight boundary segment.
Result fitOpen(const std::vector<Point>& trace,const std::vector<Point>& fixedClosure,
    const std::vector<Constraint>&,double nativeSpacing,bool axisymmetric,int samplesPerControl=48);
// RZ curve from one axis endpoint to one fixed contact (either input order).
// The closure runs back along the physical boundary and ends along the axis.
// Reflection produces an OPEN, symmetric curve between two fixed contacts;
// axis height is free, with odd r/even z parity. Only original r>=0 parent
// inventories are constrained. The returned polygon includes the reflected
// auxiliary half and must only be measured in native positive-r rectangles.
Result fitAxisToContact(const std::vector<Point>& meridianTrace,const std::vector<Point>& fixedClosure,
    const std::vector<Constraint>&,double nativeSpacing,int samplesPerControl=48);
// Positive (counterclockwise) polygon orientation, as returned by fit.
long double measure(const std::vector<Point>& polygon,const Rectangle&,bool axisymmetric);
// Reuses one ordered-arc bounding hierarchy for all rectangle integrations.
std::vector<long double> measureMany(const std::vector<Point>& polygon,
    const std::vector<Rectangle>&,bool axisymmetric);
}
