#include "Modules/Hydrodynamics/IO/LegacyNodalCache.h"
#include <istream>
#include <ostream>

// The existing native SIC2 ABI is 504 bytes of node geometry/stencils followed
// by one double. It still includes native padding and is not a portable format.
// Fail compilation if Core's layout changes; never silently reinterpret a new
// node layout as the old file format. The cross-version contract checks actual
// old records, their fields and both directions of the 512-byte representation.
static_assert(sizeof(InterfaceNode) == 504 && alignof(InterfaceNode) == 8
    && sizeof(double) == 8 && sizeof(int) == 4, "SIC2 node ABI requires an explicit migration");

namespace LegacyNodalCache
{
void writeNodes(std::ostream& stream, const std::vector<InterfaceNode>& nodes,
    const std::vector<double>& evaporationVelocity, int count)
{
    for (int i = 0; i < count; ++i)
    {
        stream.write(reinterpret_cast<const char*>(&nodes[i]), sizeof(InterfaceNode));
        // Directly added legacy geometry nodes have the old zero default.
        const double value = i < evaporationVelocity.size() ? evaporationVelocity[i] : 0.0;
        stream.write(reinterpret_cast<const char*>(&value), sizeof(double));
    }
}

void readNodes(std::istream& stream, std::vector<InterfaceNode>& nodes,
    std::vector<double>& evaporationVelocity, int count)
{
    evaporationVelocity.resize(count, 0.0);
    for (int i = 0; i < count; ++i)
    {
        stream.read(reinterpret_cast<char*>(&nodes[i]), sizeof(InterfaceNode));
        stream.read(reinterpret_cast<char*>(&evaporationVelocity[i]), sizeof(double));
    }
}
}
