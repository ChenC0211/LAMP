#pragma once

#include <vector>

#include "Core/Base/Vector2D.h"
#include "Core/Field/CellField.h"
#include "Modules/Plasma/ExplicitPIC/Deposit/CouplingWeightsCache.h"

namespace ExplicitPIC
{
namespace Gather
{

// ============================================================================
// MatchedGather 
// ----------------------------------------------------------------------------
// A gather operator that reuses the *same* coupling weights cached by CouplingWeightsCache
// (CouplingWeightsCache). This makes gather the adjoint of charge-amount
// deposition in the Euclidean inner product:
//
//   sum_p q_p * E(x_p) == sum_K Q_K * E_K
//
// when Q_K is obtained by DepositQK with the same cache and E_K is any
// cell-centered vector field.
// ============================================================================

class MatchedGather
{
public:
    MatchedGather() = default;
    ~MatchedGather() = default;

    // Gather a cell-centered vector field to particle positions.
    // EpOut will be resized to numParticles and overwritten.
    static void gatherVectorField(
        const CouplingWeightsCache& cache,
        const CellField<Vector2D>& Ecell,
        std::vector<Vector2D>& EpOut);
};

} // namespace Gather
} // namespace ExplicitPIC
