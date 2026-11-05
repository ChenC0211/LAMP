#pragma once

#include "Core/Geometry/QuadraticSurfaceMeasure.h"

namespace LAMP::Geometry::QuadraticSurfaceMeasure
{
// Explicit constrained chart: exact one arc with at least one true region
// corner endpoint. Corner parameters are fixed (dt=0); regular side-root
// endpoints retain the original moving-boundary sensitivity.
// The returned ambient covector is valid ONLY after pullback through a
// coefficient chain preserving delta F(corner)=0 at each fixed corner.
// It is not the derivative for three unconstrained coefficient variations.
Sensitivity differentiateWithFixedCornerEndpoints(const Rectangle& source,const Graph&,const Rectangle& region,
    bool axisymmetric,unsigned includedCoincidentSides=All);
}
