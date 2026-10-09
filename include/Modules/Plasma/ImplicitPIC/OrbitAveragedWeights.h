#pragma once

#include <cstddef>

class CouplingWeightsCache;

namespace ImplicitPIC
{

// Build an orbit-averaged (trapezoid) coupling cache:
//   Wbar = 0.5 * (W(x0) + W(x1))
//
// This is a minimal Stage-1 approximation for orbit-averaged overlaps.
// Later stages may replace it with event-segment integration while keeping
// the same output cache format.
//
// Returns false if inputs are incompatible.
bool buildTrapezoidAveragedCache(const CouplingWeightsCache& a,
                                const CouplingWeightsCache& b,
                                CouplingWeightsCache& out);

} // namespace ImplicitPIC
