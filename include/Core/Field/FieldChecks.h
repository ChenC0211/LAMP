#pragma once

#include <cmath>
#include <stdexcept>
#include <string>

// Storage/value checks shared by field consumers, independent of physics.
// These do not publish state, read configuration or decide physical acceptance.
namespace FieldChecks
{
template<typename FieldType>
void checkFinite(const FieldType& field, const std::string& name, int expectedSize,
                 bool positive = false, bool nonnegative = false)
{
    if (field.size() != expectedSize || !field.value)
        throw std::runtime_error(name + " storage does not match the fixed mesh");
    for (int i = 0; i < field.size(); ++i)
        if (!std::isfinite(field.value[i]) || (positive && !(field.value[i] > 0.0)) ||
            (nonnegative && field.value[i] < 0.0))
            throw std::runtime_error(name + " contains a non-finite or invalid value at index " + std::to_string(i));
}
}
