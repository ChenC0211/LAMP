#pragma once
#include "Core/Geometry/QuadraticCurveMeasure.h"
#include <array>
#include <cstddef>
#include <memory>
#include <vector>

namespace LAMP::Geometry::QuadraticMaterialVolume
{
namespace C=QuadraticCurveMeasure;
using Point=C::Point;
using Rectangle=C::Rectangle;
struct ArcMotion
{
    C::Curve source;
    // The retained physical displacements at s=0,1/2,1, never image-source.
    std::array<Point,3> displacement{};
    int liquidLeftOrientation=1; // +1/-1 relative to the original parameter.
};
struct Config
{
    long double relativeTolerance=1e-12L;
    unsigned maximumDepth=12;
    // Separate numerical accuracy of the local area/volume multiplier.
    // Third member preserves the meaning of existing {volumeTolerance,depth}.
    long double areaVolumeRelativeTolerance=1e-8L;
};
struct ArcCovector
{
    std::size_t arc=0;
    std::array<std::array<long double,2>,3> sample{};
};
struct CellMotion
{
    // The finite-gradient contraction retains tiny displacements even when independently
    // rounded endpoint volume views coincide. Neither value is an inventory.
    long double displacementVolumeChange=0,independentVolumeDifference=0;
    long double volumeIdentityResidual=0,quadratureErrorEstimate=0,nativeRoundoffEstimate=0;
    std::vector<ArcCovector> gradient;
    // Physical A/V derivative for fixed clipped endpoints, using the SAME
    // midpoint material curve and normal bubble. Not an optimizer multiplier.
    bool areaVolumeAvailable=false;
    long double areaDerivative=0,volumeDerivative=0,areaVolumeCovector=0;
    long double areaQuadratureErrorEstimate=0;
    std::size_t quadratureSamples=0;
};
// A read-only numerical functional of a genuinely closed oriented curve.
// Does not own a mesh, phase stock, accepted geometry or time layer.
// XY only. Quadrature errors are estimates, not interval certificates.
class Prepared
{
public:
    static Prepared prepare(const std::vector<ArcMotion>&);
    std::vector<CellMotion> evaluate(const std::vector<Rectangle>&,const Config& ={}) const;
    std::size_t arcCount() const noexcept;
private:
    struct Data;
    explicit Prepared(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
};
}
