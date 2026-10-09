#pragma once
#include "Core/Geometry/LevelSetHermiteGeometry.h"
class LevelSetMappedSegment;

namespace CartesianBoundaryIntegral
{
// Flux of (y-reference)e_y through the right co-normal of the ORIGINAL
// segment. XY is per unit depth; RZ includes 2*pi*r. A closed CCW boundary
// yields its area/volume. Negative mapped ranges retain offset/width charts;
// no global root endpoint is formed. The caller owns phase selection for a
// ParametricPiece. These are floating polynomial integrals, not enclosures.
long double yVolume(const LevelSetHermiteGeometry::ParametricPiece& segment,
    bool axisymmetric, double reference = 0);
long double yVolume(const LevelSetMappedSegment& segment,
    bool axisymmetric, double reference = 0);
}
