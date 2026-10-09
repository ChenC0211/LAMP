#pragma once

#include <vector>
#include <cstdint>

class QuadtreeMesh;
class HydroFittedDomain;
class HydroFittedVelocitySampler;

// Add both actual neighbors of every physical fluid boundary face, including
// inactive solid neighbors. Preserve other module marks; no halo is applied.
void protectHydroFittedBoundarySupport(const QuadtreeMesh& mesh,
    const HydroFittedDomain& domain, std::vector<std::uint8_t>& protectedCells);

// Resolution indicators only, not an alternative velocity/NS reconstruction.
// Read the accepted module velocity snapshot, never solid-only payloads.
// Values retain caller velocity/length units; inactive cells are exactly zero.
struct HydroFittedRefinementFields
{
    std::vector<double> speed;
    std::vector<double> velocityGradient;
};

HydroFittedRefinementFields evaluateHydroFittedRefinement(const QuadtreeMesh& mesh,
    const HydroFittedDomain& domain, const HydroFittedVelocitySampler& velocity);
