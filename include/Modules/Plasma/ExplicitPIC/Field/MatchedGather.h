#pragma once

#include <cstddef>
#include <vector>

#include "Core/Base/Vector2D.h"
#include "Core/Field/CellField.h"
#include "Modules/Plasma/ExplicitPIC/Deposit/CouplingWeightsCache.h"

namespace ExplicitPIC
{
namespace Field
{

// ============================================================================
// MatchedGather 
// ----------------------------------------------------------------------------
// Gather operators that reuse the *same* CSR coupling weights cached in CouplingWeightsCache
// (CouplingWeightsCache). This makes gather the discrete adjoint of charge-amount
// deposition (DepositQK) in the Euclidean inner product:
//
//   sum_p q_p * E(x_p)  ==  sum_K Q_K * E_K
//
// provided Q_K is deposited with the same cache and E(x_p) is gathered with the
// same cache.
// ============================================================================

class MatchedGather
{
public:
    struct Diagnostics
    {
        std::size_t numParticles = 0;
        std::size_t nnz = 0;
        int minNNZ = 0;
        int maxNNZ = 0;
        double avgNNZ = 0.0;
        bool validCSR = false;
    };

    // Gather a cell-centered Vector2D field to particle positions:
    //   Ep[p] = sum_e w[e] * Ecell[cellIndex[e]]
    // EpOut will be resized to cache.numParticles() and overwritten.
    static Diagnostics gatherVector2D(
        const CouplingWeightsCache& cache,
        const CellField<Vector2D>& Ecell,
        std::vector<Vector2D>& EpOut,
        bool printDiagnostics = true);

    // Optional: gather a cell-centered scalar field to particles.
    static Diagnostics gatherScalar(
        const CouplingWeightsCache& cache,
        const CellField<double>& fcell,
        std::vector<double>& fpOut,
        bool printDiagnostics = true);

private:
    static Diagnostics computeDiagnostics_(
        const CouplingWeightsCache& cache,
        bool csrLooksValid);

    static void printDiagnostics_(
        const char* tag,
        const Diagnostics& d);
};

} // namespace Field
} // namespace ExplicitPIC
