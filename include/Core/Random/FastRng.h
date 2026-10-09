#pragma once

#include <cstdint>
#include <limits>

// Fast RNG utilities (SplitMix64). Not cryptographic.

namespace FastRng
{

// SplitMix64 mix/hash function.
inline std::uint64_t Mix64(std::uint64_t x) noexcept
{
    x += 0x9e3779b97f4a7c15ull;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ull;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebull;
    x = (x ^ (x >> 31));
    return x;
}

// Advance SplitMix64 state and return next value.
inline std::uint64_t SplitMix64Next(std::uint64_t& state) noexcept
{
    state += 0x9e3779b97f4a7c15ull;
    std::uint64_t z = state;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ull;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebull;
    return z ^ (z >> 31);
}

// SplitMix64 generator satisfying UniformRandomBitGenerator.
class SplitMix64Rng
{
public:
    using result_type = std::uint64_t;

    explicit SplitMix64Rng(std::uint64_t seed = 0) noexcept : state_(seed) {}

    void seed(std::uint64_t s) noexcept { state_ = s; }

    std::uint64_t state() const noexcept { return state_; }

    static constexpr result_type min() { return std::numeric_limits<result_type>::min(); }
    static constexpr result_type max() { return std::numeric_limits<result_type>::max(); }

    result_type operator()() noexcept { return SplitMix64Next(state_); }

private:
    std::uint64_t state_ = 0;
};

// Deterministically derive a per-particle seed from a global seed and a key (e.g., particle index).
inline std::uint64_t SeedFrom(std::uint64_t global_seed, std::uint64_t key) noexcept
{
    const std::uint64_t s = Mix64(global_seed ^ Mix64(key));
    return (s == 0) ? 0x9e3779b97f4a7c15ull : s;
}

} // namespace FastRng
