#pragma once

namespace LAMP::Geometry::QuadraticCurveMeasure
{
struct Curve;
struct PhysicalMeasureDifference
{
    long double delta=0; // Signed after-before physical measure on each curve's [0,1].
    long double quadratureErrorEstimate=0; // Numerical estimate, not a rigorous bound.
};
// Paired integration retains differences of authoritative coefficients BEFORE
// native rounding. Each curve keeps its own normalized parameter on [0,1].
// Adaptive error uses the signed difference/absolute difference integrand,
// never the two large individual physical measures. No physical owner is made.
PhysicalMeasureDifference physicalMeasureDifference(const Curve& before,const Curve& after,
    bool axisymmetric=false,long double relativeTolerance=1e-13L,unsigned maximumDepth=16);
}
