#pragma once

#include "Core/Geometry/QuadraticCellCut.h"
#include <array>
#include <vector>

namespace LAMP::Geometry::QuadraticSurfaceMeasure
{
using Point = PlanarCellCut::Point;
using Rectangle = PlanarCellCut::Rectangle;
using Graph = QuadraticCellCut::Reconstruction;
enum BoundaryMask : unsigned { Left=1, Right=2, Bottom=4, Top=8, All=15 };
struct Sample
{
    Point position{}, unitNormal{}; // Liquid -> gas; follows the graph's sign.
    long double meridionalCurvature=0, curvature=0;
};
struct Integral
{
    long double length=0, physicalMeasure=0;
    long double quadratureErrorEstimate=0; // Numerical estimate, not a certificate.
};
struct Patch
{
    long double lower=0, upper=0; // Original graph's physical tangent coordinate.
    Point begin{}, end{};
    unsigned beginBoundary=0, endBoundary=0;
    Integral integral;
};
struct Result
{
    std::vector<Patch> patches;
    Integral integral;
};
struct Sensitivity
{
    Result surface;
    // Derivatives w.r.t. original physical (offset,slope,second), including
    // moving intersections with the region's boundary. Rejects nonregular
    // endpoint/corner configurations; no frozen-endpoint approximation.
    std::array<long double,3> derivative{};
};
// Derived coefficients for ONE retained exact polynomial arc. Derivative
// arrays are indexed by physical (offset,slope,second), in the same units as
// Sensitivity. No quadrature is performed while preparing this plan.
struct MomentPatch
{
    long double lower=0,upper=0,m0=0,m1=0,scale=0;
    Point begin{},end{};unsigned beginBoundary=0,endBoundary=0;
    std::array<long double,3> radius{};
    std::array<long double,3> m0Derivative{},m1Derivative{},scaleDerivative{};
    std::array<std::array<long double,3>,3> radiusDerivative{}; // [radius power][physical variable]
};
struct MomentValues
{
    std::array<long double,3> moments{},firstDerivative{},lastDerivative{};
};
class MomentPlan
{
public:
    static MomentPlan prepare(const Rectangle& source,const Graph&,const Rectangle& region,
        bool axisymmetric,unsigned includedCoincidentSides=All);
    const std::vector<MomentPatch>& patches() const noexcept;
    bool axisymmetric() const noexcept;
    // Values and partials must belong to the SAME supplied moment model.
    // Error estimate zero means no quadrature HERE, not an interpolation/error
    // certificate. Final physical consumers must certify with evaluate/differentiate.
    Sensitivity apply(const std::vector<MomentValues>&) const;
private:
    struct Data;
    explicit MomentPlan(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
};
// Same original graph coefficients and centre as QuadraticCellCut. No refit,
// mesh/state ownership, phase inventory, curvature differencing or IO.
// RZ uses full 2*pi*r ds and div(n) including the regular axis limit.
Sample sample(const Rectangle& source, const Graph&, long double tangent, bool axisymmetric);
Integral integrate(const Rectangle& source, const Graph&, long double lower,
    long double upper, bool axisymmetric);
// Restrict the ORIGINAL curve to an actual subrectangle. Boundary masks only
// assign arcs identically coincident with a side; they never discard an
// isolated contact or move a root. AMR partitions can use half-open ownership
// (e.g. lower sides owned) to avoid counting an aligned interface twice.
Result evaluate(const Rectangle& source, const Graph&, const Rectangle& region,
    bool axisymmetric, unsigned includedCoincidentSides=All);
Sensitivity differentiate(const Rectangle& source, const Graph&, const Rectangle& region,
    bool axisymmetric, unsigned includedCoincidentSides=All);
}
