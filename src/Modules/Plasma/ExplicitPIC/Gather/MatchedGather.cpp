#include "Modules/Plasma/ExplicitPIC/Gather/MatchedGather.h"

#include <cassert>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace ExplicitPIC
{
namespace Gather
{

void MatchedGather::gatherVectorField(
    const CouplingWeightsCache& cache,
    const CellField<Vector2D>& Ecell,
    std::vector<Vector2D>& EpOut)
{
    const int nP = static_cast<int>(cache.numParticles());
    // Avoid unconditional zero-fill; the loop overwrites every entry.
    EpOut.resize(static_cast<size_t>(nP));

    const auto& offsets = cache.offsets();
    const auto& cellIndex = cache.cellIndices();
    const auto& w = cache.weights();

    // offsets has size nP+1.
    assert(offsets.size() == static_cast<size_t>(nP + 1));

    // Each particle gather is independent.
    #ifdef _OPENMP
    const bool useOmp = (nP >= 20000);
    #pragma omp parallel for if(useOmp) schedule(static)
    #endif
    for (int p = 0; p < nP; ++p)
    {
        Vector2D sum(0.0, 0.0);
        const size_t a = offsets[static_cast<size_t>(p)];
        const size_t b = offsets[static_cast<size_t>(p + 1)];
        for (size_t i = a; i < b; ++i)
        {
            const int ci = cellIndex[i];
            if (ci < 0 || ci >= Ecell.size()) continue;
            sum += Ecell.value[ci] * w[i];
        }
        EpOut[static_cast<size_t>(p)] = sum;
    }
}

} // namespace Gather
} // namespace ExplicitPIC
