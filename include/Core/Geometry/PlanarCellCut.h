#pragma once
#include "Core/Base/Vector2D.h"
#include <vector>
#include <cstddef>
#include <array>
#include <optional>
#include <limits>

namespace LAMP::Geometry::PlanarCellCut
{
// Absolute volume precision of the inverse, normalized by cell capacity.
inline constexpr long double relativeVolumeResolution=128*std::numeric_limits<long double>::epsilon();
// Actual axis-aligned physical rectangle. The plane is
// normal dot (x - rectangle centre) <= offset (negative/liquid side).
// Normal must be finite and nonzero; it is normalized internally. Offset has
// length units relative to that UNIT normal. No nominal mesh spacing is used.
struct Rectangle { Vector2D lower, upper; };
struct Reconstruction
{
    Vector2D normal;
    long double offset = 0, volume = 0, residual = 0;
};
long double capacity(const Rectangle& cell, bool axisymmetric);
long double volume(const Rectangle& cell, Vector2D normal, long double offset,
    bool axisymmetric);
// Solves one monotone PLIC offset. Fractions outside [0,1] reject; no clipping.
// RZ integrates 2*pi*r over the clipped polygon, not midpoint-radius area.
Reconstruction reconstruct(const Rectangle& cell, Vector2D normal,
    long double fraction, bool axisymmetric);
struct Point { long double x=0, y=0; };
enum class RectangleCorner { LowerLeft=0, LowerRight=1, UpperRight=2, UpperLeft=3 };
struct CornerReconstruction
{
    Reconstruction plane;
    Point localCorner;
    RectangleCorner corner;
};
// Solves one monotone normal parameter while the plane passes through the
// specified actual rectangle corner. True makes the adjacent horizontal edge
// negative and the adjacent vertical edge positive, except their common corner;
// false reverses the phase sides. Fractions outside [0,1] reject, without clipping.
// The returned near-unit double normal is the actual argument to volume/interfaceSegment;
// offset uses their same internal normalization and corner expression, so the
// exact local corner remains a segment endpoint for representable mixed cuts.
// Pure fractions may give a coincident edge, not a physical-interface claim.
// Unrepresentable positive mixed geometry rejects; no phase stock is created.
CornerReconstruction reconstructThroughCorner(const Rectangle&,RectangleCorner,
    long double fraction,bool axisymmetric,bool negativeAlongHorizontalEdge=true);
struct Overlap { long double total=0, liquid=0; };
// Validated convex polygon intersection with an actual rectangle. Coordinates
// remain in the supplied common frame; used by planar and curved donor cuts.
std::vector<Point> clipPolygonInFrame(const Rectangle& cell,const std::vector<Point>& polygon,Point origin);
// Convex swept polygon, with coordinates RELATIVE TO this rectangle's centre.
// Clip to the actual cell, then to its PLIC plane; both moments use the same
// XY/RZ measure. Either cyclic orientation is accepted, concavity rejects.
// Degenerate polygons have zero measure. No union/overlap repair is performed.
Overlap overlap(const Rectangle& cell, const std::vector<Point>& polygon,
    Vector2D normal, long double offset, bool axisymmetric);
// Same kernel, with one explicit GLOBAL frame origin shared by all donors.
// Polygon coordinates stay relative to this origin; do not translate a thin
// swept polygon separately into each cell's frame. The PLIC offset still uses
// the rectangle centre, as in reconstruct().
Overlap overlapInFrame(const Rectangle& cell,const std::vector<Point>& polygon,Point origin,
    Vector2D normal,long double offset,bool axisymmetric);
// Plane/rectangle segment in cell-centred local coordinates. No result when
// fewer than two distinct boundary intersections are representable. Invalid
// geometry still throws; no fraction threshold or coordinate snapping is used.
std::optional<std::array<Point,2>> interfaceSegment(const Rectangle& cell,Vector2D normal,long double offset);

// Ordinary area and first moments in the SUPPLIED coordinate frame. Either
// cyclic orientation is accepted; concave/crossed/nonfinite polygons reject.
// Repeated/collinear vertices and lower-dimensional polygons have zero area.
struct PolygonMoments { long double area=0, firstX=0, firstY=0; };
PolygonMoments polygonMoments(const std::vector<Point>& polygon);
// RZ uses the actual r=origin.x+p.x, integrating 2*pi*r over the polygon.
// XY ignores the translation. Negative radius is rejected, never clipped.
long double polygonMeasureInFrame(const std::vector<Point>& polygon,Point origin,bool axisymmetric);
// n dot p <= offset in the supplied frame; n is NOT normalized. The result
// stays in that frame and is CCW. A zero normal or unresolved crossing rejects.
std::vector<Point> clipConvexPolygon(const std::vector<Point>& polygon,Point normal,long double offset);

// Strictly opposite finite scalar signs on a nondegenerate segment. Uses the
// same robust crossing as polygon/triangle cuts; exact zero endpoints are sites.
Point linearZeroCrossing(Point a,Point b,long double va,long double vb);

enum class TriangleRegion { Negative, Positive, Cut, AllZeroScalar };
// One affine scalar on an actual nondegenerate triangle. No transported jet,
// mesh/state identity or geometric inventory is created. The output has fixed
// storage; accurate internal expansions may spill beyond their inline storage.
struct LinearTriangleCut
{
    TriangleRegion kind=TriangleRegion::AllZeroScalar;
    std::array<Point,4> negativeVertices{};
    std::size_t negativeCount=0;
    // Bits refer to INPUT vertex i / edge (i,(i+1)%3), preserving exact zeros.
    unsigned zeroVertexMask=0, zeroEdgeMask=0;
    // Present for a positive-length zero line, including a coincident edge.
    // Ordered with negative phase on the left; normal points negative->positive.
    std::optional<std::array<Point,2>> interface;
    Point unitNormal{}; // Meaningful only when interface is present.
};
LinearTriangleCut cutTriangle(const std::array<Point,3>& triangle,const std::array<long double,3>& scalar);
// AllZeroScalar has positive geometric area but no defined phase, so this
// query throws. Ordinary vertex/edge zeros have zero area and retain tags.
long double negativeMeasure(const LinearTriangleCut&,Point origin,bool axisymmetric);
// Same cut/measure generator without segment/normal preparation, for repeated
// local volume queries. Fixed output storage; internal exact arithmetic may
// allocate for expansion tails. AllZeroScalar still rejects.
long double negativeMeasure(const std::array<Point,3>& triangle,const std::array<long double,3>& scalar,
    Point origin,bool axisymmetric);

// Explicit actual sites: each coordinate is the exact sum of these two
// publicly retained fields in the caller's common frame. No conversion to
// Point is authoritative, and no hidden ideal intersection is measured.
struct AnchoredPoint {Point anchor{},offset{};};
Point approximateAnchoredPoint(AnchoredPoint); // Display/diagnostic only.
Point anchoredDifference(AnchoredPoint a,AnchoredPoint b); // Once-rounded a-b.
// Complete exact q-frame-center, rounded directly once to binary64 RNE.
// Derived coordinates may round to zero; nonfinite inputs/final overflow reject.
Vector2D anchoredDifferenceDouble(AnchoredPoint query,AnchoredPoint frameOrigin,
    Point localCenter={});
bool sameAnchoredPoint(AnchoredPoint,AnchoredPoint); // Exact coordinate sums.
AnchoredPoint anchoredZeroCrossing(AnchoredPoint a,AnchoredPoint b,
    long double valueA,long double valueB);
struct AnchoredTriangleCut
{
    TriangleRegion kind=TriangleRegion::AllZeroScalar;
    std::array<AnchoredPoint,4> negativeVertices{};
    std::size_t negativeCount=0;
    unsigned zeroVertexMask=0,zeroEdgeMask=0;
    std::optional<std::array<AnchoredPoint,2>> interface;
    Point unitNormal{}; // Metric metadata, negative phase -> positive phase.
};
AnchoredTriangleCut cutAnchoredTriangle(const std::array<AnchoredPoint,3>&,
    const std::array<long double,3>&);
long double negativeAnchoredMeasure(const AnchoredTriangleCut&,Point origin,bool axisymmetric);
long double negativeAnchoredMeasure(const std::array<AnchoredPoint,3>&,
    const std::array<long double,3>&,Point origin,bool axisymmetric);
PolygonMoments anchoredPolygonMoments(const std::vector<AnchoredPoint>&);
long double anchoredPolygonMeasureInFrame(const std::vector<AnchoredPoint>&,
    Point origin,bool axisymmetric);
std::vector<AnchoredPoint> clipAnchoredPolygonInFrame(const Rectangle&,
    const std::vector<AnchoredPoint>&,Point origin);
std::vector<AnchoredPoint> clipAnchoredConvexPolygon(const std::vector<AnchoredPoint>&,
    Point normal,long double offset);
bool anchoredPolygonInsideBox(const std::vector<AnchoredPoint>&,
    Point origin,const Rectangle&);

}
