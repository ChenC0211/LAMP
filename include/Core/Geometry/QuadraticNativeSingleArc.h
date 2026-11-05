#pragma once
#include "Core/Geometry/QuadraticCellCut.h"
#include "Core/Geometry/QuadraticCornerBoundary.h"
#include "Core/Geometry/QuadraticSurfaceMeasure.h"
#include <array>
#include <memory>
#include <optional>
#include <string>

namespace LAMP::Geometry::QuadraticNativeSingleArc
{
namespace Q=QuadraticCellCut;
namespace S=QuadraticSurfaceMeasure;
using Point=Q::Point;
struct Request
{
    Q::Rectangle box{};
    Q::Shape shape{}; // Actual liquid phase and the same frozen physical frame.
    std::array<Point,2> ports{};
    long double targetPhaseVolume=0; // Liquid iff volumeLiquid, otherwise gas.
    bool volumeLiquid=true;
    long double branchEta=0,relativeVolumeTolerance=1e-12L;
    double minimumSlope=-8,maximumSlope=8;
};
struct SeedSample
{
    std::array<long double,3> value{}; // Actual-phase world tangent x/y, h*kappa.
    long double slope=0,tangentOverH=0,radius=0;
};
namespace Detail {struct Data;}
enum class Reason { Ready, Unsupported, UncertainRoundoff, NearEvent };
struct Preparation;
// Numerical trial/cache only. Never an exact graph, geometry/state owner or
// acceptance certificate. Same-input exact reconstruction remains mandatory.
struct Trial
{
    Request input;
    Q::Shape shape{};
    long double offset=0,eta=0,h=0,frameD=0,capacity=0,selectedVolume=0;
    std::array<long double,3> selectedDerivative{}; // wrt(q0/h,q1,q2*h).
    S::MomentPatch moment;
    std::array<SeedSample,2> seed{}; // Sorted positive tangent patch order.
    std::array<Q::CornerBoundaryKind,4> corners{}; // Actual liquid phase.
    std::array<std::optional<Q::StationaryFaceBoundary>,2> stationary{};
    long double relativeVolumeErrorBound=0;
private:
    std::shared_ptr<const Detail::Data> data_;
    friend Preparation prepare(const Request&);
    friend S::Sensitivity area(const Trial&,const S::MomentValues&);
    friend Q::FixedPortVolumeBoundaryJet direction(const Trial&,Point,Point);
    friend long double pairedAreaDifference(const Trial&,const S::MomentValues&,const Trial&,const S::MomentValues&);
};
struct Preparation
{
    std::optional<Trial> value;
    Reason reason=Reason::Unsupported;
    std::string detail;
};
// Invalid inputs throw. Unsupported/uncertain/near-event predicates select
// the original strict chart BEFORE use; no catch/failure-driven geometry swap.
Preparation prepare(const Request&);
// Same caller-supplied moment model and physical units as MomentPlan::apply.
// No quadrature and no claim that model/interval checks replace final Core.
S::Sensitivity area(const Trial&,const S::MomentValues&);
// Same target, fixed input port order; velocities stay on actual open faces.
// cornerImplicit uses the SELECTED volume phase, unlike Trial::corners.
Q::FixedPortVolumeBoundaryJet direction(const Trial&,Point firstVelocity,Point lastVelocity);
// Signed after-before XY area of this same supplied moment model.
long double pairedAreaDifference(const Trial&,const S::MomentValues&,const Trial&,const S::MomentValues&);
}
