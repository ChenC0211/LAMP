#include "Core/Mesh/QuadtreeMesh.h"
#include <stdexcept>

void QuadtreeMesh::expandCellMask(std::vector<std::uint8_t>& mask, int layers) const
{
    if (mask.size() != static_cast<std::size_t>(cellIndexNum()) || layers < 0)
        throw std::invalid_argument("Mesh protection requires a native cell mask and nonnegative layers");
    for (int layer = 0; layer < layers; ++layer)
    {
        auto expanded = mask;
        for (int i = 0; i < cellIndexNum(); ++i)
            if (mask[i])
                for (const auto& side : cells(cellIndexToLabel[i]).neiCell)
                    for (int label : side)
                        if (label >= 0) expanded[cells(label).index()] = 1;
        mask.swap(expanded);
    }
}
