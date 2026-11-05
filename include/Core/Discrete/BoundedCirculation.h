#pragma once
#include <array>
#include <cstddef>
#include <vector>

namespace LAMP::Discrete
{
struct BoundedCirculationArc
{
    std::size_t from=0,to=0;
    double lower=0,upper=0,initial=0;
};
struct BoundedCirculationResult
{
    // Original arc order; every represented flow is inside its input bounds.
    // The internal exact circulation is feasible. Rounding may leave the
    // stated net-inflow residuals; consumers must enforce their own budget.
    std::vector<double> flows;
    std::vector<std::array<double,2>> roundedBalance;
    std::size_t augmentations=0;
};
// Deterministic shortest augmenting paths with exact binary64 sums/differences.
// Invalid data, unsupported FP modes and genuinely infeasible networks reject.
// No epsilon residual elimination, scalar physics, mesh or mutable callback.
BoundedCirculationResult solveBoundedCirculation(std::size_t nodes,
    const std::vector<BoundedCirculationArc>& arcs);
}
