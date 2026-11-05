#pragma once

#include "Core/Geometry/ParametricHermiteCurve.h"
#include "Core/Geometry/QuadraticSurfaceMeasure.h"
#include <array>
#include <optional>
#include <memory>
#include <vector>

namespace LAMP::Geometry::Detail {struct QuadraticDerivativeNativeCache;}

namespace LAMP::Geometry::QuadraticCurveMeasure
{
using Point=Parametric::Point;
using Rectangle=PlanarCellCut::Rectangle;
using Integral=QuadraticSurfaceMeasure::Integral;
struct Curve
{
    // Direct coefficient construction retains the original legacy power basis.
    // On fromSamples() curves these are lossy diagnostics, never authority.
    std::array<Point,3> coefficients{};
    struct EndpointBasis { Point first,last,bubble; };
    Curve()=default;
    Curve(std::array<Point,3> legacy):coefficients(legacy){}
    // X(s)=(1-s)*first+s*last+s*(1-s)*bubble. Endpoint identity is
    // represented directly, rather than recovered by power cancellation.
    static Curve fromSamples(Point first,Point middle,Point last);
    static Curve fromGraph(std::shared_ptr<const Detail::QuadraticParametricData>);
    bool usesExactBasis() const noexcept{return static_cast<bool>(exact_);}
    const std::shared_ptr<const Detail::QuadraticParametricData>& exactData() const noexcept{return exact_;}
    bool exactViewUnchanged() const noexcept;
    bool usesEndpointBasis() const noexcept{return endpoints_.has_value();}
    const EndpointBasis* endpointBasis() const noexcept{return endpoints_?&*endpoints_:nullptr;}
private:
    std::optional<EndpointBasis> endpoints_;
    std::shared_ptr<const Detail::QuadraticParametricData> exact_;
    // Read-only derived native evaluator/constants, never curve authority.
    // Incomplete here: public users do not acquire the private exact backend.
    std::shared_ptr<const Detail::QuadraticDerivativeNativeCache> derivativeNative_;
    friend struct Detail::QuadraticDerivativeNativeCache;
    std::array<Point,3> exactView_{};
};
// Same polynomial value/parameter, without reading lossy diagnostic views.
bool samePolynomial(const Curve&,const Curve&);
Curve affineMapped(const Curve&,Point sourceAnchor,Point imageAnchor,const std::array<double,4>& derivative);
Curve displacedSamples(const Curve&,Point firstDisplacement,Point middleDisplacement,Point lastDisplacement);
struct Bounds {Point lower{},upper{};};
Point position(const Curve&,double parameter);
// Same physical parameter and authoritative basis. Order >=3 is exact zero.
Point derivative(const Curve&,double parameter,int order=1);
Bounds bounds(const Curve&);
// Original parameter roots of coordinate(axis)=value; no fitted power trace.
std::vector<double> coordinateRoots(const Curve&,int axis,double value,double lower=0,double upper=1);
std::vector<std::array<double,2>> intervals(const Curve&,const Rectangle&,
    unsigned includedCoincidentSides=QuadraticSurfaceMeasure::All);
// Positive adaptive quadrature of ds / 2*pi*r ds on the retained parameter.
// Numerical estimate, no interval certificate or topology/phase assertion.
Integral integrate(const Curve&,double lower=0,double upper=1,bool axisymmetric=false,
    long double relativeTolerance=1e-13L,unsigned maximumDepth=16);
}