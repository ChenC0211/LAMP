#pragma once

#include "Core/Geometry/ActiveHermitePatch.h"
#include <array>
#include <vector>

namespace LAMP::Geometry
{
// Explicit shared data on actual rectangular leaves. This is not a fit,
// contact condition or AMR transfer algorithm. A snapshot owns its copy.
// Derivatives are PHYSICAL, before multiplication by each leaf's dimensions.
struct ActiveQuarticHermitePatch
{
    using NodalDof = ActiveHermitePatch::NodalDof;
    using CornerSector = ActiveHermitePatch::CornerSector;
    struct EdgeMidpoint
    {
        // Native vertex indices, ordered in increasing x (horizontal) or y
        // (vertical). Identifies a full leaf side, including a coarse side
        // consisting of two native fine edges. Not an inactive parent index.
        std::array<int,2> vertices{{-1,-1}};
        // Data at the EXACT parameter midpoint of this endpoint pair.
        double value=0.0, normalDerivative=0.0;
    };
    std::vector<std::uint8_t> activeCells; // full native cell-index length
    std::vector<NodalDof> nodes;          // full native vertex-index length
    std::vector<EdgeMidpoint> edgeMidpoints;
    std::vector<double> cellCenters;     // full native cell-index length
    // Same explicitly declared replacement corner jet as the cubic port.
    // Only its four physical corner data are substituted; all edge/center
    // values remain shared. Exact active-edge C0 and a full QUARTIC pure-phase
    // certificate are required, in the original locally same-level lip band.
    std::vector<CornerSector> cornerSectors;

    // Only regular used vertices, maximal active side midpoints and active
    // centers are independent. Hanging-node and fine-subedge payloads are
    // ignored; their values AND normal traces come from the same coarse side.
    // Unused payloads may be NaN. Duplicate/foreign side identities reject.
    // Cell centers and edge midpoints are exact local parameter locations,
    // not independently rounded physical sample positions. Corner sectors
    // are explicit input, not a generated wall/contact-angle condition.
};

} // namespace LAMP::Geometry
