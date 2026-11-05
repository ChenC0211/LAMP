#pragma once

#include "Modules/Interface/VofSmoothGeometry.h"

namespace LAMP::Interface::Detail
{
// Same represented-curve kernel used by volume fitting and final quadrature.
// Internal contract port: the caller supplies the actual cell and regular
// graph arc; fitting may use its signed polynomial extension. With output,
// the original positive physical quadrature-weight checks apply.
std::array<long double,2> measureSmoothArcPhases(const VofSmoothArc&,
    const GeometricVofPlicCell&,bool axisymmetric,
    std::array<std::vector<VofSmoothWeightedPoint>,2>* output=nullptr);
// The original native dyadic geometry certificate (same depth/bounds/gates).
bool certifySmoothArc(const VofSmoothArc&,const Geometry::PlanarCellCut::Rectangle&,
    unsigned& maximumDepth);
// Exact-degree XY Bernstein/Green integral of the actual world controls.
// Both phases have their own anchored-width integral and outer strips.
std::array<long double,2> analyticSmoothArcPhasesXY(const VofSmoothArc&,
    const GeometricVofPlicCell&);
struct SmoothDirectXYWork {
    std::size_t domainLeaves=0,candidateConstructions=0,geometryCertifications=0,
        phaseQuadratureBuilds=0,fitPhaseMeasurements=0,volumeNewtonSteps=0;
};
struct SmoothDirectXYResult {
    VofSmoothArc arc;
    std::array<long double,2> measured{};
    std::array<std::vector<VofSmoothWeightedPoint>,2> phaseQuadrature;
    unsigned maximumCertificateSubdivisionDepth=0;
};
// Existing quartic family and world-control authority only. Fixed outer
// controls/nodes/graph axis; full depth-8 sufficient rectangle, exact-minor
// target first, deterministic metric selection and ONE native candidate.
// One original certificate and one independent positive quadrature publish
// the result. No fit measurement, Newton, history, reflection or neighbor scan.
// Insufficient world-native representation remains a typed rejection.
SmoothDirectXYResult constructDirectSmoothArcXY(const VofSmoothArc&,
    const GeometricVofPlicCell&,long double relativeVolumeTolerance,
    SmoothDirectXYWork* work=nullptr);
}
