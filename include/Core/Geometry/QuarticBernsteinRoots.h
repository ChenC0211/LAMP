#pragma once
#include <array>
#include <vector>

namespace LAMP::Geometry
{
struct QuarticBernsteinRoots
{
    std::vector<double> roots;
    bool identicallyZero=false;
};
// Strict out-of-line isolation on [0,1], including exact endpoint roots.
// Unresolved interior multiple roots reject; no coefficient truncation.
QuarticBernsteinRoots isolateQuarticBernstein(const std::array<double,5>& controls);
}
