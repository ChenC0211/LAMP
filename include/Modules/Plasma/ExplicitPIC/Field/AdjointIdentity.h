#pragma once

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "Core/Base/Vector2D.h"
#include "Core/Field/CellField.h"
#include "Modules/Plasma/ExplicitPIC/Deposit/CouplingWeightsCache.h"

namespace ExplicitPIC
{
namespace Field
{

// ============================================================================
// verifyAdjointIdentity 
// ----------------------------------------------------------------------------
// Check the discrete adjoint identity when deposit and gather share the same
// coupling weights \tilde W_{pK}:
//
//   LHS = \sum_p q_p E(x_p) ,   E(x_p) = \sum_K E_K \tilde W_{pK}
//   RHS = \sum_K Q_K E_K   ,   Q_K    = \sum_p q_p \tilde W_{pK}
//
// This header provides:
//   - AdjointIdentityResult (lhs/rhs/diff + abs/rel per component)
//   - verifyAdjointIdentity(...) with Kahan summation for numerical robustness.
// ============================================================================

struct AdjointIdentityResult
{
    Vector2D lhs{0.0, 0.0};
    Vector2D rhs{0.0, 0.0};
    Vector2D diff{0.0, 0.0};

    double absDiffX = 0.0;
    double absDiffY = 0.0;
    double relDiffX = 0.0;
    double relDiffY = 0.0;

    // vector norm relative error
    double relDiffVec = 0.0;

    // Diagnostics
    std::size_t numParticles = 0;
    std::size_t nnz = 0;
    bool sizesConsistent = true;
};

namespace detail
{
struct KahanSum
{
    double sum = 0.0;
    double c = 0.0;
    inline void add(double x)
    {
        const double y = x - c;
        const double t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
};
} // namespace detail

inline AdjointIdentityResult verifyAdjointIdentity(
    const std::vector<double>& qp,
    const CouplingWeightsCache& cache,
    const CellField<double>& Qcell,
    const CellField<Vector2D>& Ecell,
    const std::vector<Vector2D>& Ep,
    std::ostream& os = std::cout)
{
    using detail::KahanSum;

    AdjointIdentityResult r;
    r.numParticles = cache.numParticles();
    r.nnz = cache.nnz();

    const std::size_t np = cache.numParticles();
    if (qp.size() != np || Ep.size() != np)
        r.sizesConsistent = false;

    // LHS: sum_p q_p Ep[p]
    KahanSum sxL, syL;
    const std::size_t nLp = std::min({qp.size(), Ep.size(), np});
    for (std::size_t p = 0; p < nLp; ++p)
    {
        sxL.add(qp[p] * Ep[p].x);
        syL.add(qp[p] * Ep[p].y);
    }
    r.lhs = Vector2D(sxL.sum, syL.sum);

    // RHS: sum_K Q_K E_K
    KahanSum sxR, syR;
    const int nK = std::min(Qcell.size(), Ecell.size());
    for (int ci = 0; ci < nK; ++ci)
    {
        const double qk = Qcell.value[ci];
        sxR.add(qk * Ecell.value[ci].x);
        syR.add(qk * Ecell.value[ci].y);
    }
    r.rhs = Vector2D(sxR.sum, syR.sum);

    r.diff = Vector2D(r.lhs.x - r.rhs.x, r.lhs.y - r.rhs.y);
    r.absDiffX = std::abs(r.diff.x);
    r.absDiffY = std::abs(r.diff.y);

    const double denomX = std::max(std::abs(r.rhs.x), 1e-300);
    const double denomY = std::max(std::abs(r.rhs.y), 1e-300);
    r.relDiffX = r.absDiffX / denomX;
    r.relDiffY = r.absDiffY / denomY;

    const double denomV = std::max(r.rhs.mod(), 1e-300);
    r.relDiffVec = r.diff.mod() / denomV;

    os << std::setprecision(16);
    os << "[AdjointIdentity] np=" << r.numParticles
       << ", nnz=" << r.nnz
       << ", sizesConsistent=" << (r.sizesConsistent ? "true" : "false")
       << std::endl;
    os << "[AdjointIdentity] LHS=sum_p q_p E(x_p) = (" << r.lhs.x << ", " << r.lhs.y << ")" << std::endl;
    os << "[AdjointIdentity] RHS=sum_K Q_K E_K   = (" << r.rhs.x << ", " << r.rhs.y << ")" << std::endl;
    os << "[AdjointIdentity] diff=LHS-RHS        = (" << r.diff.x << ", " << r.diff.y << ")" << std::endl;
    os << "[AdjointIdentity] absDiff=(" << r.absDiffX << ", " << r.absDiffY << ")" << std::endl;
    os << "[AdjointIdentity] relDiffX=" << r.relDiffX << ", relDiffY=" << r.relDiffY
       << ", relDiffVec=" << r.relDiffVec << std::endl;

    return r;
}

} // namespace Field
} // namespace ExplicitPIC
