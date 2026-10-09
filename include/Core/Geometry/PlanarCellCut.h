#pragma once
#include "Core/Base/Vector2D.h"
#include <vector>
#include <array>
#include <optional>
#include <limits>

namespace Quantum::Geometry::PlanarCellCut
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
}
