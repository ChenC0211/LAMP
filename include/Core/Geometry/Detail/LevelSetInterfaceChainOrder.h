#pragma once

#include <array>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace Quantum::Geometry::Detail
{
namespace
{
// Connectivity only: the caller owns root identity, local incidence tests,
// boundary declarations and the first oriented segment. Both layouts use
// this same one-component walk. No coordinates or tolerances enter it.
struct InterfaceChainVisit
{
    int cell, begin, end;
};
inline std::vector<InterfaceChainVisit> orderInterfaceChain(
    const std::vector<std::array<int, 2>> &connections, const std::vector<std::vector<int>> &nodeCells,
    int start, int firstCell, std::size_t expectedCells, bool closed)
{
    const auto require = [](bool valid, const char *reason) {
        if (!valid)
            throw std::runtime_error(std::string("[LevelSetInterfaceChain] ") + reason);
    };
    std::vector<InterfaceChainVisit> result;
    std::vector<bool> visited(connections.size(), false);
    int node = start, cell = firstCell;
    for (std::size_t k = 0; k < expectedCells; ++k)
    {
        require(cell >= 0 && static_cast<std::size_t>(cell) < connections.size() && !visited[cell],
                "chain repeats a cell or closes before covering all components");
        visited[cell] = true;
        const auto pair = connections[cell];
        require(pair[0] == node || pair[1] == node, "broken cell/root incidence");
        const int end = pair[0] == node ? pair[1] : pair[0];
        require(end >= 0 && static_cast<std::size_t>(end) < nodeCells.size() && end != node,
                "invalid chain endpoint");
        result.push_back({cell, node, end});
        node = end;
        const auto &incident = nodeCells[end];
        require(incident.size() == 1 || incident.size() == 2, "branch or missing endpoint incidence");
        require(incident[0] == cell || (incident.size() == 2 && incident[1] == cell),
                "endpoint does not own the current cell");
        if (incident.size() == 1)
        {
            require(!closed && k + 1 == expectedCells && node != start,
                    "chain ends early or leaves disconnected components");
            break;
        }
        cell = incident[0] == cell ? incident[1] : incident[0];
    }
    require(!result.empty() && result.size() == expectedCells &&
                (closed ? node == start && cell == firstCell : node != start && nodeCells[node].size() == 1),
            "interface did not form one complete chain");
    return result;
}
} // namespace
} // namespace Quantum::Geometry::Detail
