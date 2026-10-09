#pragma once

#include "Core/Geometry/LevelSetHermiteGeometry.h"

class LevelSetMappedSegment;

// Boundary moments of represented straight segments revolved around r=x=0.
// These are polynomial integrals, not Cartesian midpoint approximations.
namespace AxisymmetricBoundaryIntegral
{
// Integral of (z-axialReference)e_z through the segment, including 2*pi*r.
// The right-hand co-normal is (dz,-dr): CCW traversal therefore has the
// outward normal of its enclosed region. Summing a CLOSED boundary gives
// its volume; a single segment contribution can have either sign.
//
// Uses the retained anchor/root parameter chart directly, including reversed
// tBegin/tEnd. The phase flag is NOT consumed; the caller selects wet pieces.
// Never materializes/subtracts global root points or cumulative primitives.
// Zero span and dr=0 return zero; invalid radii/parameters/nonfinite arithmetic
// reject. This is a floating polynomial evaluation, not an interval enclosure.
long double axialVolume(const LevelSetHermiteGeometry::ParametricPiece& segment,
    double axialReference = 0.0);

// Negative-phase moment of the mapped segment's ORIGINAL physical input
// direction. Each range retains a separate near-end offset and positive width;
// fromEnd changes only the integration chart, not the right-hand co-normal.
// No rounded root endpoints or cumulative primitive differences are formed.
// This is a floating polynomial integral, not an interval enclosure.
long double axialVolume(const LevelSetMappedSegment& segment,
    double axialReference = 0.0);
}
