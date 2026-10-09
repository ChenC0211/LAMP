#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "Modules/Plasma/ExplicitPIC/Deposit/ParticleBinner.h"

// Forward declarations to keep this header lightweight.
template<typename T>
class CellField;

class CouplingWeightsCache;

// ----------------------------------------------------------------------------
// DepositQK
//
// Integrated-charge deposition on cell-centered unknowns.
//
// ----------------------------------------------------------------------------

class DepositQK
{
public:
    struct Workspace
    {
        std::vector<std::size_t> particleOfEntry;
        std::vector<double> sumThread;
        ExplicitPIC::ParticleBinner binner;
    };

    struct Diagnostics
    {
        // Sum over particle charges: sum_p q_p
        double sumQP = 0.0;

        // Sum over deposited charge increment: sum over all cache entries (qp[p] * w_{pK}).
        // If per-particle weights are normalized, this should match sumQP up to roundoff.
        double sumQgrid = 0.0;

        // |sumQgrid - sumQP|
        double absDiff = 0.0;

        // absDiff / max(|sumQP|, eps)
        double relDiff = 0.0;
    };

    DepositQK() = default;

    // Placeholder hook: does nothing.
    void depositPlaceholder() const;

    // Deposit per-particle charges qp into cell field Q using cached weights.
    //
    // qp[p] is the total charge represented by particle p (e.g. q_s * w_p).
    //
    // This routine performs CSR-based deposition directly by cell index.
    Diagnostics deposit(const std::vector<double>& qp,
                        const CouplingWeightsCache& cache,
                        CellField<double>* Q) const;

    // Deposit per-particle charges qp and number-amount wp into two fields
    // using a single cache traversal (avoids duplicate binning work).
    Diagnostics depositDual(const std::vector<double>& qp,
                            const std::vector<double>& wp,
                            const CouplingWeightsCache& cache,
                            CellField<double>* Q,
                            CellField<double>* N) const;

    Diagnostics depositDual(const std::vector<double>& qp,
                            const std::vector<double>& wp,
                            const CouplingWeightsCache& cache,
                            CellField<double>* Q,
                            CellField<double>* N,
                            Workspace* workspace) const;

    // Debug: dump a single particle's weights (delegates to cache).
    bool dumpParticleWeightsCSV(const CouplingWeightsCache& cache,
                                const std::string& filePath,
                                std::size_t particleId) const;
};
