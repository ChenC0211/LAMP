#pragma once

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>

namespace LAMP::FieldTransfer
{
// Call before entering an interpolation OpenMP region. Saved source indices
// refer to the PREVIOUS field, not the new mesh's active count or capacity.
template<class Entity>
void validateSources(const Entity& entity, int previousSize, const char* location)
{
    const std::string context = std::string(location) + " mesh transfer, label=" + std::to_string(entity.label());
    if (entity.relatedPreIndex.empty() || entity.relatedPreIndex.size() != entity.relatedPreCoord.size())
        throw std::runtime_error(context + ": missing or inconsistent interpolation sources");
    for (std::size_t i = 0; i < entity.relatedPreIndex.size(); ++i)
    {
        const int index = entity.relatedPreIndex[i];
        if (index < 0 || index >= previousSize)
            throw std::runtime_error(context + ": previous source index=" + std::to_string(index)
                + " is outside [0," + std::to_string(previousSize) + ")");
        const auto point = entity.relatedPreCoord[i];
        if (!std::isfinite(point.x) || !std::isfinite(point.y))
            throw std::runtime_error(context + ": non-finite previous source coordinate");
    }
}
}
