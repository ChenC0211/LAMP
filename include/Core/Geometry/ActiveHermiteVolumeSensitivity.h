#pragma once

#include <array>
#include <vector>
#include "Core/Geometry/LevelSetInterfacePanels.h"

namespace ActiveHermiteVolumeSensitivity
{
// Four cell corners [SW,SE,NE,NW], each [value,dx,dy,dxy]. Derivatives
// refer to the PHYSICAL nodal DOFs, not the already scaled Hermite cache.
using Coefficients = std::array<long double, 16>;

// First variation of negative-phase volume on the supplied resolved subarcs:
//   delta V = - integral (B * delta dof) / |grad phi| dA.
// physicalWeight already contains the caller's planar/RZ measure; no radial
// factor is applied again. dx/dy/dxy basis columns include hx/hy/hx*hy.
// Samples must retain the same active cell's local Frame and raw gradient.
// Their Snapshot provenance, cellSize and actual-control-volume clipping are
// the caller's responsibility; these cannot be inferred from copied samples.
// No global point subtraction, normalisation of phi, or assumed unit gradient.
// Empty samples return zeros. Invalid data or unrepresentable arithmetic
// throws before returning; the inputs are never changed.
//
// This is GL8-based local linear sensitivity, not an exact volume integral or
// a nonlinear conservation solve. Global shared/sector DOF identification,
// pin/axis constraints and any reconstruction chain rule belong to the caller.
Coefficients assemble(int cellI, int cellJ, Vector2D cellSize,
    const std::vector<LevelSetInterfacePanels::ArcSample>& samples);
}
