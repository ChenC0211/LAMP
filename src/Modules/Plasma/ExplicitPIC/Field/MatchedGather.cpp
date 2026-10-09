#include "Modules/Plasma/ExplicitPIC/Field/MatchedGather.h"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "Core/IO/ConsoleControl.h"

namespace ExplicitPIC
{
namespace Field
{

namespace
{
inline void updateMinMax(int v, int& vmin, int& vmax)
{
    vmin = std::min(vmin, v);
    vmax = std::max(vmax, v);
}
} // namespace

MatchedGather::Diagnostics MatchedGather::computeDiagnostics_(
    const CouplingWeightsCache& cache,
    const bool csrLooksValid)
{
    Diagnostics d;
    d.numParticles = cache.numParticles();
    d.nnz = cache.nnz();
    d.validCSR = csrLooksValid;

    if (!csrLooksValid || d.numParticles == 0)
    {
        d.minNNZ = 0;
        d.maxNNZ = 0;
        d.avgNNZ = 0.0;
        return d;
    }

    const auto& off = cache.offsets();
    int minNNZ = std::numeric_limits<int>::max();
    int maxNNZ = 0;

    for (std::size_t p = 0; p < d.numParticles; ++p)
    {
        const int nnzRow = off[p + 1] - off[p];
        updateMinMax(nnzRow, minNNZ, maxNNZ);
    }

    d.minNNZ = (minNNZ == std::numeric_limits<int>::max()) ? 0 : minNNZ;
    d.maxNNZ = maxNNZ;
    d.avgNNZ = (d.numParticles > 0) ? static_cast<double>(d.nnz) / static_cast<double>(d.numParticles) : 0.0;

    return d;
}

void MatchedGather::printDiagnostics_(
    const char* tag,
    const Diagnostics& d)
{
    if (!ConsoleControl::enabled(ConsoleControl::Level::kDebug))
        return;
    std::cout << std::setprecision(6) << std::fixed
              << "[MatchedGather::" << tag << "]"
              << " np=" << d.numParticles
              << " nnz=" << d.nnz
              << " avgNNZ=" << d.avgNNZ
              << " minNNZ=" << d.minNNZ
              << " maxNNZ=" << d.maxNNZ
              << " csr=" << (d.validCSR ? "ok" : "bad")
              << std::defaultfloat << std::endl;
}

MatchedGather::Diagnostics MatchedGather::gatherVector2D(
    const CouplingWeightsCache& cache,
    const CellField<Vector2D>& Ecell,
    std::vector<Vector2D>& EpOut,
    const bool printDiagnostics)
{
    const std::size_t np = cache.numParticles();
    // Avoid an unconditional zero-fill (EpOut.assign) on the hot path.
    // The loop below overwrites every entry when CSR is valid.
    EpOut.resize(np);

    const auto& off = cache.offsets();
    const auto& idx = cache.cellIndices();
    const auto& w = cache.weights();

    const bool csrLooksValid = (!off.empty() && off.size() == np + 1 && idx.size() == w.size());
    Diagnostics diag;
    if (printDiagnostics)
    {
        diag = computeDiagnostics_(cache, csrLooksValid);
    }
    else
    {
        diag.numParticles = np;
        diag.nnz = cache.nnz();
        diag.validCSR = csrLooksValid;
        diag.minNNZ = 0;
        diag.maxNNZ = 0;
        diag.avgNNZ = 0.0;
    }

    if (!csrLooksValid || np == 0)
    {
        // Preserve legacy behavior: in invalid-CSR cases, EpOut is zeroed.
        std::fill(EpOut.begin(), EpOut.end(), Vector2D(0.0, 0.0));
        if (printDiagnostics)
            printDiagnostics_("gatherVector2D", diag);
        return diag;
    }

    // Defensive checks: offsets range must not exceed nnz.
    assert(static_cast<std::size_t>(off.back()) <= w.size());

    #ifdef _OPENMP
    const bool useOmp = (np >= 20000);
    #pragma omp parallel for if(useOmp) schedule(static)
    #endif
    for (std::int64_t pp = 0; pp < static_cast<std::int64_t>(np); ++pp)
    {
        const std::size_t p = static_cast<std::size_t>(pp);
        const int begin = off[p];
        const int end = off[p + 1];

        Vector2D sum(0.0, 0.0);
        for (int e = begin; e < end; ++e)
        {
            const int ci = idx[static_cast<std::size_t>(e)];
            if (ci < 0 || ci >= Ecell.size())
                continue;
            sum += Ecell.value[ci] * w[static_cast<std::size_t>(e)];
        }
        EpOut[p] = sum;
    }

    if (printDiagnostics)
        printDiagnostics_("gatherVector2D", diag);

    return diag;
}

MatchedGather::Diagnostics MatchedGather::gatherScalar(
    const CouplingWeightsCache& cache,
    const CellField<double>& fcell,
    std::vector<double>& fpOut,
    const bool printDiagnostics)
{
    const std::size_t np = cache.numParticles();
    // Avoid unconditional zero-fill on the valid-CSR hot path.
    fpOut.resize(np);

    const auto& off = cache.offsets();
    const auto& idx = cache.cellIndices();
    const auto& w = cache.weights();

    const bool csrLooksValid = (!off.empty() && off.size() == np + 1 && idx.size() == w.size());
    Diagnostics diag;
    if (printDiagnostics)
    {
        diag = computeDiagnostics_(cache, csrLooksValid);
    }
    else
    {
        diag.numParticles = np;
        diag.nnz = cache.nnz();
        diag.validCSR = csrLooksValid;
        diag.minNNZ = 0;
        diag.maxNNZ = 0;
        diag.avgNNZ = 0.0;
    }

    if (!csrLooksValid || np == 0)
    {
        // Preserve legacy behavior: in invalid-CSR cases, fpOut is zeroed.
        std::fill(fpOut.begin(), fpOut.end(), 0.0);
        if (printDiagnostics)
            printDiagnostics_("gatherScalar", diag);
        return diag;
    }

    assert(static_cast<std::size_t>(off.back()) <= w.size());

    #ifdef _OPENMP
    const bool useOmp = (np >= 20000);
    #pragma omp parallel for if(useOmp) schedule(static)
    #endif
    for (std::int64_t pp = 0; pp < static_cast<std::int64_t>(np); ++pp)
    {
        const std::size_t p = static_cast<std::size_t>(pp);
        const int begin = off[p];
        const int end = off[p + 1];

        double sum = 0.0;
        for (int e = begin; e < end; ++e)
        {
            const int ci = idx[static_cast<std::size_t>(e)];
            if (ci < 0 || ci >= fcell.size())
                continue;
            sum += fcell.value[ci] * w[static_cast<std::size_t>(e)];
        }
        fpOut[p] = sum;
    }

    if (printDiagnostics)
        printDiagnostics_("gatherScalar", diag);

    return diag;
}

} // namespace Field
} // namespace ExplicitPIC
