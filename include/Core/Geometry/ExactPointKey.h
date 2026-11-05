#pragma once

#include "Core/Base/Vector2D.h"
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>

namespace LAMP::Geometry
{
// Identity of a physical query, not geometric proximity. Preserve signed
// zeros and every mantissa bit; never merge adjacent representable points.
using ExactPointKey = std::array<std::uint64_t, 2>;

inline ExactPointKey exactPointKey(Vector2D point) noexcept
{
    static_assert(sizeof(double) == sizeof(std::uint64_t), "binary64 point keys required");
    ExactPointKey result;
    std::memcpy(&result[0], &point.x, sizeof(double));
    std::memcpy(&result[1], &point.y, sizeof(double));
    return result;
}

struct ExactPointKeyHash
{
    std::size_t operator()(const ExactPointKey& key) const noexcept
    {
        // Dyadic mesh coordinates have many identical LOW mantissa bits.
        // Avalanche each full word before combining: a shifts-only combiner
        // degenerates on implementations with power-of-two bucket counts.
        return static_cast<std::size_t>(mix(mix(key[0]) ^ (mix(key[1]) + 0x9e3779b97f4a7c15ULL)));
    }
  private:
    static std::uint64_t mix(std::uint64_t value) noexcept
    {
        value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
        value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
        return value ^ (value >> 31);
    }
};
}
