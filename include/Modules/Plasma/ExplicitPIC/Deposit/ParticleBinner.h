#pragma once

#include <cstdint>
#include <vector>

namespace ExplicitPIC {

// A light-weight binning helper for deposit operations.
//
// Given a list of deposition entries (cellIndices) of length nnz, it builds a
// cell->entries mapping using counting sort:
//   - cellOffsets_[c] gives the starting index in perm_ for cell c.
//   - perm_[k] gives the entry id e (index into the original CSR arrays).
//
// Invalid indices (<0 or >=numCells) are ignored (not inserted into perm_).
//
// This helper is intentionally header-only to avoid touching build scripts for
// standalone tools that compile DepositQK.cpp explicitly.
class ParticleBinner
{
public:
    ParticleBinner() = default;

    void build(const std::vector<int>& cellIndices, const int numCells)
    {
        numCells_ = numCells;
        if (numCells_ <= 0)
        {
            counts_.clear();
            cellOffsets_.assign(1, 0);
            perm_.clear();
            touchedCells_.clear();
            return;
        }

        counts_.assign(static_cast<std::size_t>(numCells_), 0);
        touchedCells_.clear();
        touchedCells_.reserve(1024);

        // Count valid entries per cell and record touched cells.
        for (std::size_t e = 0; e < cellIndices.size(); ++e)
        {
            const int c = cellIndices[e];
            if (c < 0 || c >= numCells_)
                continue;
            if (counts_[static_cast<std::size_t>(c)] == 0)
                touchedCells_.push_back(c);
            ++counts_[static_cast<std::size_t>(c)];
        }

        // Prefix sum to offsets.
        cellOffsets_.assign(static_cast<std::size_t>(numCells_ + 1), 0);
        int running = 0;
        for (int c = 0; c < numCells_; ++c)
        {
            cellOffsets_[static_cast<std::size_t>(c)] = running;
            running += counts_[static_cast<std::size_t>(c)];
        }
        cellOffsets_[static_cast<std::size_t>(numCells_)] = running;

        // Counting-sort permutation (entry ids) by cell.
        perm_.assign(static_cast<std::size_t>(running), 0);
        // Reuse counts_ as write pointers.
        writePtr_.assign(static_cast<std::size_t>(numCells_), 0);
        for (int c = 0; c < numCells_; ++c)
            writePtr_[static_cast<std::size_t>(c)] = cellOffsets_[static_cast<std::size_t>(c)];

        for (std::size_t e = 0; e < cellIndices.size(); ++e)
        {
            const int c = cellIndices[e];
            if (c < 0 || c >= numCells_)
                continue;
            const int pos = writePtr_[static_cast<std::size_t>(c)]++;
            perm_[static_cast<std::size_t>(pos)] = static_cast<int>(e);
        }
    }

    int numCells() const { return numCells_; }
    int nnzValid() const
    {
        if (cellOffsets_.empty())
            return 0;
        return cellOffsets_.back();
    }

    const std::vector<int>& cellOffsets() const { return cellOffsets_; }
    const std::vector<int>& perm() const { return perm_; }
    const std::vector<int>& touchedCells() const { return touchedCells_; }

private:
    int numCells_ = 0;
    std::vector<int> counts_;
    std::vector<int> writePtr_;
    std::vector<int> cellOffsets_;
    std::vector<int> perm_;
    std::vector<int> touchedCells_;
};

} // namespace ExplicitPIC
