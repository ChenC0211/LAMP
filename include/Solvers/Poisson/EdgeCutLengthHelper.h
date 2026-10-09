#pragma once

#include "Core/Field/VertField.h"

class QuadtreeMesh;

namespace Solvers
{
namespace Poisson
{

// -----------------------------------------------------------------------------
// Helper:
// Fill Edge::cut() / Edge::cutLength() consistently from a given level-set field.
//
// - Uses a simple linear root model along each edge (robust enough for diagnostics/regression).
// - Writes both Edge::isCut[levelSetI] and Edge::cut()/cutLength().
// - Does NOT change any core mesh/geometry types.
// -----------------------------------------------------------------------------
struct EdgeCutLengthHelper
{
    static void computeLinear(QuadtreeMesh& mesh, const VertField<double>& psi, int levelSetI, double clampEps = 1e-12);
};

} // namespace Poisson
} // namespace Solvers
