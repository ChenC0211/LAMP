#pragma once

#include "Core/Geometry/QuadraticCellCut.h"
#include <array>

namespace LAMP::Geometry::QuadraticCellCut
{
struct CornerBoundaryKind
{
    long double value=0; // Same selected signed implicit F/h.
    // d(F/h)/d(inward physical distance/h), dimensionless, physical x/y.
    std::array<long double,2> inwardDerivative{};
    bool isolatedClosed=false;
};
// Pure local classification, corner order LL/LR/UR/UL. Exact, strictly
// same-sign inward derivatives and F of that sign (or zero) imply that
// moving only the constant term to F=0 gives a locally isolated corner.
// This is sufficient, not necessary. No graph is moved or published, no
// exact-contact assertion is made about the current graph, and no stock,
// arc-count, global feasibility, shared C0 or acceptance certificate is given.
CornerBoundaryKind cornerBoundaryKind(const Rectangle&,const Reconstruction& currentExactGraph,unsigned corner);
}
