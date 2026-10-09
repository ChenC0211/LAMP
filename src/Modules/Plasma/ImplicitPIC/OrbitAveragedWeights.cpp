#include "Modules/Plasma/ImplicitPIC/OrbitAveragedWeights.h"

#include <algorithm>
#include <unordered_map>
#include <utility>
#include <vector>

#include "Modules/Plasma/ExplicitPIC/Deposit/CouplingWeightsCache.h"

namespace ImplicitPIC
{

bool buildTrapezoidAveragedCache(const CouplingWeightsCache& a,
                                const CouplingWeightsCache& b,
                                CouplingWeightsCache& out)
{
    if (a.numParticles() != b.numParticles())
        return false;

    const std::size_t nP = a.numParticles();
    const auto& offA = a.offsets();
    const auto& offB = b.offsets();

    if (offA.size() != nP + 1u || offB.size() != nP + 1u)
        return false;

    std::vector<int> outOffsets(nP + 1u, 0);
    std::vector<int> outCell;
    std::vector<double> outW;

    // A conservative reserve: union could be up to nnzA+nnzB.
    outCell.reserve(a.nnz() + b.nnz());
    outW.reserve(a.nnz() + b.nnz());

    std::unordered_map<int, double> acc;
    std::vector<std::pair<int, double>> entries;

    for (std::size_t p = 0; p < nP; ++p)
    {
        acc.clear();

        const int a0 = offA[p];
        const int a1 = offA[p + 1u];
        const int b0 = offB[p];
        const int b1 = offB[p + 1u];

        if (a0 < 0 || a1 < a0 || b0 < 0 || b1 < b0)
            return false;

        // Accumulate 0.5 * W_a
        {
            const auto& ci = a.cellIndices();
            const auto& w = a.weights();
            if (static_cast<std::size_t>(a1) > ci.size() || static_cast<std::size_t>(a1) > w.size())
                return false;
            for (int k = a0; k < a1; ++k)
                acc[ci[static_cast<std::size_t>(k)]] += 0.5 * w[static_cast<std::size_t>(k)];
        }

        // Accumulate 0.5 * W_b
        {
            const auto& ci = b.cellIndices();
            const auto& w = b.weights();
            if (static_cast<std::size_t>(b1) > ci.size() || static_cast<std::size_t>(b1) > w.size())
                return false;
            for (int k = b0; k < b1; ++k)
                acc[ci[static_cast<std::size_t>(k)]] += 0.5 * w[static_cast<std::size_t>(k)];
        }

        // Materialize in deterministic order.
        entries.clear();
        entries.reserve(acc.size());
        for (const auto& kv : acc)
        {
            if (kv.second != 0.0)
                entries.emplace_back(kv.first, kv.second);
        }
        std::sort(entries.begin(), entries.end(), [](const auto& x, const auto& y) {
            return x.first < y.first;
        });

        for (const auto& e : entries)
        {
            outCell.push_back(e.first);
            outW.push_back(e.second);
        }

        outOffsets[p + 1u] = static_cast<int>(outCell.size());
    }

    out.assignCSR(nP, std::move(outOffsets), std::move(outCell), std::move(outW));
    return true;
}

} // namespace ImplicitPIC
