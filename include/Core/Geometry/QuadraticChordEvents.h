#pragma once
#include "Core/Geometry/QuadraticCellCut.h"
#include "Core/Geometry/UnitPolynomialRoots.h"

namespace LAMP::Geometry::QuadraticChordEvents
{
using Point=PlanarCellCut::Point;
// Endpoints move affinely with tau in [0,1]; the chord is bilinear in
// (its own parameter,tau). Points use the same physical Cartesian frame.
struct Sweep {Point begin0,end0,begin1,end1;};
// All graph/chord intersections, including isolated tangential contacts on a
// gas-only cell edge. These are arrangement features, not liquid intervals.
UnitPolynomialRoots intersections(const PlanarCellCut::Rectangle&,
    const QuadraticCellCut::Reconstruction&,Point begin,Point end);
// Candidate times when the chord passes through a fixed arrangement vertex.
// Coefficient sums/products retain exact expansions of the input coordinates;
// roots retain the Core enclosures. Exact zero discriminants permit algebraic
// deflation of temporal multiple roots, without weakening spatial-chart rules.
// They do not certify a nonlinear characteristic approximated by this sweep.
UnitPolynomialRoots vertex(const Sweep&,Point);
// Times of an INTERIOR tangency to a cell-centred quadratic graph. Events on
// its cell boundary are handled by vertex(); callers own rectangle clipping.
// A graph with zero second derivative has no isolated interior tangency.
UnitPolynomialRoots tangencies(const PlanarCellCut::Rectangle&,
    const QuadraticCellCut::Reconstruction&,const Sweep&);
// Matching scalar event functions for polishing a seed on a nonlinear trace.
long double vertexValue(Point begin,Point end,Point vertex);
long double tangencyValue(const PlanarCellCut::Rectangle&,
    const QuadraticCellCut::Reconstruction&,Point begin,Point end);
}
