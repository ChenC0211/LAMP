#pragma once

#include <cstdint>
#include <vector>

class QuadtreeMesh;

// Exact native geometry/connectivity identity extracted from Hydro's fitted
// domain. No physical boundary/material data or mutable interpolation caches.
// Capture once per mesh state; validate before reusing fields or operators.
class NativeMeshIdentity
{
public:
    explicit NativeMeshIdentity(const QuadtreeMesh& mesh);
    void validate(const QuadtreeMesh& mesh) const;
private:
    std::vector<std::uint64_t> words_;
};

