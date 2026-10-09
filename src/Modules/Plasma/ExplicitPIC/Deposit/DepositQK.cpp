#include "Modules/Plasma/ExplicitPIC/Deposit/DepositQK.h"

#include "Modules/Plasma/ExplicitPIC/Deposit/CouplingWeightsCache.h"
#include "Modules/Plasma/ExplicitPIC/Deposit/ParticleBinner.h"

// CellField is only needed in the implementation.
#include "Core/Field/CellField.h"
#include "Core/IO/ConsoleControl.h"

#include <algorithm>
#include <cstdint>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

// Kahan summation for improved roundoff behavior in global charge diagnostics.
inline void kahanAdd(const double x, double& sum, double& c)
{
    const double y = x - c;
    const double t = sum + y;
    c = (t - sum) - y;
    sum = t;
}

inline void printChargeDiagnostics(const DepositQK::Diagnostics& d,
                                   const std::size_t np,
                                   const std::size_t nnz)
{
    const double deltaQ = d.sumQgrid - d.sumQP;
    if (!ConsoleControl::enabled(ConsoleControl::Level::kTrace))
        return;
    std::cout << std::setprecision(17) << std::scientific
              << "[DepositQK] np=" << np
              << " nnz=" << nnz
              << " sumQP=" << d.sumQP
              << " sumQgrid=" << d.sumQgrid
              << " dQ=" << deltaQ
              << " absDiff=" << d.absDiff
              << " relDiff=" << d.relDiff
              << std::defaultfloat << std::endl;
}

} // namespace

void DepositQK::depositPlaceholder() const
{
    // Intentionally empty (placeholder hook).
}

DepositQK::Diagnostics DepositQK::deposit(const std::vector<double>& qp,
                                          const CouplingWeightsCache& cache,
                                          CellField<double>* Q) const
{
    // NOTE: cellIndices in the cache are *cell indices* compatible with CellField storage.
    Diagnostics diag;

    // Defensive: allow calling with null Q to only compute sums.
    const bool doWrite = (Q != nullptr);

    const std::size_t np = qp.size();
    diag.sumQP = 0.0;
    double cQP = 0.0;
    for (double q : qp)
        kahanAdd(q, diag.sumQP, cQP);

    // If cache is empty or particle counts mismatch, we still report sumQP.
    if (cache.numParticles() == 0 || cache.offsets().empty() || cache.numParticles() != np)
    {
        diag.sumQgrid = 0.0;
        diag.absDiff = std::abs(diag.sumQgrid - diag.sumQP);
        const double denom = std::max(1e-300, std::abs(diag.sumQP));
        diag.relDiff = diag.absDiff / denom;
        printChargeDiagnostics(diag, np, cache.nnz());
        return diag;
    }

    const auto& off = cache.offsets();
    const auto& idx = cache.cellIndices();
    const auto& w = cache.weights();

    diag.sumQgrid = 0.0;
    double cQgrid = 0.0;

    // CSR loop: for each particle p, entries in [off[p], off[p+1]).
    //
    // Reduce atomic contention by re-binning *entries* by target cell and
    // performing a per-cell reduction (each cell written once by a unique
    // thread). This turns scatter-add into a gather-by-cell.
#ifdef _OPENMP
    const int maxThreads = omp_get_max_threads();
    std::vector<double> sumThread(static_cast<std::size_t>(maxThreads), 0.0);

    // Map each entry e to its particle p (needed when reordering entries by cell).
    const std::size_t nnz = cache.nnz();
    const std::size_t invalidParticle = std::numeric_limits<std::size_t>::max();
    std::vector<std::size_t> particleOfEntry(nnz, invalidParticle);

    // Pass 1: fill particleOfEntry and compute sumQgrid (all entries, including invalid).
    // This loop is safe to parallelize since each particle owns a disjoint CSR range.
#pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        double sumLocal = 0.0;
        double cLocal = 0.0;

#pragma omp for schedule(static)
        for (std::int64_t p64 = 0; p64 < static_cast<std::int64_t>(np); ++p64)
        {
            const std::size_t p = static_cast<std::size_t>(p64);
            const double q = qp[p];
            const int begin = off[p];
            const int end = off[p + 1];
            for (int e = begin; e < end; ++e)
            {
                particleOfEntry[static_cast<std::size_t>(e)] = p;
                const double contrib = q * w[static_cast<std::size_t>(e)];
                kahanAdd(contrib, sumLocal, cLocal);
            }
        }
        sumThread[static_cast<std::size_t>(tid)] = sumLocal;
    }

    // Reduce thread-local sums (serial Kahan) for stable diagnostics.
    for (int t = 0; t < maxThreads; ++t)
        kahanAdd(sumThread[static_cast<std::size_t>(t)], diag.sumQgrid, cQgrid);

    // If we are only computing sums, we are done.
    if (!doWrite || Q == nullptr)
    {
        diag.absDiff = std::abs(diag.sumQgrid - diag.sumQP);
        const double denom = std::max(1e-300, std::abs(diag.sumQP));
        diag.relDiff = diag.absDiff / denom;
        printChargeDiagnostics(diag, np, cache.nnz());
        return diag;
    }

    // Pass 2: build entry bins by target cell and perform per-cell reductions.
    const int qSize = Q->size();
    ExplicitPIC::ParticleBinner binner;
    binner.build(idx, qSize);
    const auto& cellOffsets = binner.cellOffsets();
    const auto& perm = binner.perm();
    const auto& touched = binner.touchedCells();

    // Each touched cell is processed by exactly one thread, so the final write
    // is race-free (no atomics on Q).
#pragma omp parallel for schedule(static)
    for (std::int64_t it64 = 0; it64 < static_cast<std::int64_t>(touched.size()); ++it64)
    {
        const int cellIndex = touched[static_cast<std::size_t>(it64)];
        const int begin = cellOffsets[static_cast<std::size_t>(cellIndex)];
        const int end = cellOffsets[static_cast<std::size_t>(cellIndex + 1)];
        double sumCell = 0.0;
        for (int k = begin; k < end; ++k)
        {
            const int e = perm[static_cast<std::size_t>(k)];
            const std::size_t p = particleOfEntry[static_cast<std::size_t>(e)];
            if (p == invalidParticle)
                continue;
            sumCell += qp[p] * w[static_cast<std::size_t>(e)];
        }
        // Bounds are guaranteed by the binner (valid indices only).
        Q->value[cellIndex] += sumCell;
    }
#else
    for (std::size_t p = 0; p < np; ++p)
    {
        const double q = qp[p];
        const int begin = off[p];
        const int end = off[p + 1];
        for (int e = begin; e < end; ++e)
        {
            const int cellIndex = idx[e];
            const double contrib = q * w[e];
            kahanAdd(contrib, diag.sumQgrid, cQgrid);
            if (doWrite)
            {
                if (cellIndex >= 0 && cellIndex < Q->size())
                    Q->value[cellIndex] += contrib;
            }
        }
    }
#endif

    diag.absDiff = std::abs(diag.sumQgrid - diag.sumQP);
    const double denom = std::max(1e-300, std::abs(diag.sumQP));
    diag.relDiff = diag.absDiff / denom;

    // Print global charge conservation diagnostic (ΔQ) to roundoff.
    printChargeDiagnostics(diag, np, cache.nnz());
    return diag;
}

DepositQK::Diagnostics DepositQK::depositDual(const std::vector<double>& qp,
                                              const std::vector<double>& wp,
                                              const CouplingWeightsCache& cache,
                                              CellField<double>* Q,
                                              CellField<double>* N) const
{
    return depositDual(qp, wp, cache, Q, N, /*workspace=*/nullptr);
}

DepositQK::Diagnostics DepositQK::depositDual(const std::vector<double>& qp,
                                              const std::vector<double>& wp,
                                              const CouplingWeightsCache& cache,
                                              CellField<double>* Q,
                                              CellField<double>* N,
                                              Workspace* workspace) const
{
    Diagnostics diag;
    const bool doWriteQ = (Q != nullptr);
    const bool doWriteN = (N != nullptr);

    const std::size_t np = qp.size();
    diag.sumQP = 0.0;
    double cQP = 0.0;
    for (double q : qp)
        kahanAdd(q, diag.sumQP, cQP);

    if (cache.numParticles() == 0 || cache.offsets().empty() || cache.numParticles() != np || wp.size() != np)
    {
        diag.sumQgrid = 0.0;
        diag.absDiff = std::abs(diag.sumQgrid - diag.sumQP);
        const double denom = std::max(1e-300, std::abs(diag.sumQP));
        diag.relDiff = diag.absDiff / denom;
        printChargeDiagnostics(diag, np, cache.nnz());
        return diag;
    }

    const auto& off = cache.offsets();
    const auto& idx = cache.cellIndices();
    const auto& w = cache.weights();

    diag.sumQgrid = 0.0;
    double cQgrid = 0.0;

#ifdef _OPENMP
    Workspace localWorkspace;
    Workspace& ws = workspace ? *workspace : localWorkspace;

    const int maxThreads = omp_get_max_threads();
    ws.sumThread.assign(static_cast<std::size_t>(maxThreads), 0.0);

    const std::size_t nnz = cache.nnz();
    ws.particleOfEntry.resize(nnz);

#pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        double sumLocal = 0.0;
        double cLocal = 0.0;

#pragma omp for schedule(static)
        for (std::int64_t p64 = 0; p64 < static_cast<std::int64_t>(np); ++p64)
        {
            const std::size_t p = static_cast<std::size_t>(p64);
            const double q = qp[p];
            const int begin = off[p];
            const int end = off[p + 1];
            for (int e = begin; e < end; ++e)
            {
                ws.particleOfEntry[static_cast<std::size_t>(e)] = p;
                const double contrib = q * w[static_cast<std::size_t>(e)];
                kahanAdd(contrib, sumLocal, cLocal);
            }
        }
        ws.sumThread[static_cast<std::size_t>(tid)] = sumLocal;
    }

    for (int t = 0; t < maxThreads; ++t)
        kahanAdd(ws.sumThread[static_cast<std::size_t>(t)], diag.sumQgrid, cQgrid);

    if ((!doWriteQ && !doWriteN) || (Q == nullptr && N == nullptr))
    {
        diag.absDiff = std::abs(diag.sumQgrid - diag.sumQP);
        const double denom = std::max(1e-300, std::abs(diag.sumQP));
        diag.relDiff = diag.absDiff / denom;
        printChargeDiagnostics(diag, np, cache.nnz());
        return diag;
    }

    int qSize = 0;
    if (doWriteQ && Q)
        qSize = Q->size();
    else if (doWriteN && N)
        qSize = N->size();

    ws.binner.build(idx, qSize);
    const auto& cellOffsets = ws.binner.cellOffsets();
    const auto& perm = ws.binner.perm();
    const auto& touched = ws.binner.touchedCells();

#pragma omp parallel for schedule(static)
    for (std::int64_t it64 = 0; it64 < static_cast<std::int64_t>(touched.size()); ++it64)
    {
        const int cellIndex = touched[static_cast<std::size_t>(it64)];
        const int begin = cellOffsets[static_cast<std::size_t>(cellIndex)];
        const int end = cellOffsets[static_cast<std::size_t>(cellIndex + 1)];
        double sumCellQ = 0.0;
        double sumCellN = 0.0;
        for (int k = begin; k < end; ++k)
        {
            const int e = perm[static_cast<std::size_t>(k)];
            const std::size_t p = ws.particleOfEntry[static_cast<std::size_t>(e)];
            const double wpe = w[static_cast<std::size_t>(e)];
            sumCellQ += qp[p] * wpe;
            sumCellN += wp[p] * wpe;
        }
        if (doWriteQ && Q)
            Q->value[cellIndex] += sumCellQ;
        if (doWriteN && N)
            N->value[cellIndex] += sumCellN;
    }
#else
    for (std::size_t p = 0; p < np; ++p)
    {
        const double q = qp[p];
        const double wpv = wp[p];
        const int begin = off[p];
        const int end = off[p + 1];
        for (int e = begin; e < end; ++e)
        {
            const int cellIndex = idx[e];
            const double we = w[e];
            const double contribQ = q * we;
            kahanAdd(contribQ, diag.sumQgrid, cQgrid);
            if (doWriteQ && Q)
            {
                if (cellIndex >= 0 && cellIndex < Q->size())
                    Q->value[cellIndex] += contribQ;
            }
            if (doWriteN && N)
            {
                if (cellIndex >= 0 && cellIndex < N->size())
                    N->value[cellIndex] += wpv * we;
            }
        }
    }
#endif

    diag.absDiff = std::abs(diag.sumQgrid - diag.sumQP);
    const double denom = std::max(1e-300, std::abs(diag.sumQP));
    diag.relDiff = diag.absDiff / denom;
    printChargeDiagnostics(diag, np, cache.nnz());
    return diag;
}

bool DepositQK::dumpParticleWeightsCSV(const CouplingWeightsCache& cache,
                                       const std::string& filePath,
                                       const std::size_t particleId) const
{
    return cache.dumpParticleWeightsCSV(filePath, particleId);
}
