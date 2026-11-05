#pragma once
// Internal Electric detail. Only typed results and a declaration live here;
// the sole implementation stays in the strict production Q1 translation unit.
#include "Core/Geometry/PlanarCellCut.h"
#include "Core/Geometry/UnitPolynomialRoots.h"
#include <array>
#include <cstddef>
#include <optional>
#include <vector>

namespace LAMP::Electric::Detail::SurfaceChargeProjection
{
using Point=Geometry::PlanarCellCut::Point;
enum class MinimumCertificate { None,StrictEndpoint,StrictInterior };
struct Candidate
{
    Geometry::UnitPolynomialRoots::Root parameter;
    std::array<long double,2> squaredDistance{};
    int endpoint=-1;
};
struct Selection
{
    std::vector<Candidate> candidates;
    std::size_t roundedWinner=0,interiorRoots=0;
    std::optional<std::size_t> selected,overlap;
    bool uniqueInteriorCertificate=false;
    MinimumCertificate minimum=MinimumCertificate::None;
};
Selection select(const std::array<Point,5>& controls,Point query);
struct LocalMinimum
{
    std::size_t arc=0;
    std::array<std::size_t,2> nodes{};
    std::array<Point,2> endpoints{};
    std::array<long double,2> squaredDistance{};
    int endpoint=-1;
    MinimumCertificate minimum=MinimumCertificate::None;
};
// All candidates must use one query and one immutable target owner. The
// production adapter authenticates endpoint controls against its actual nodes.
// No result means that at least one comparison lacks a strict certificate.
std::optional<std::size_t> selectLocalMinimum(const std::vector<LocalMinimum>&);
}
