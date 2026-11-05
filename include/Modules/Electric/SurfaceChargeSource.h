#pragma once

#include <functional>
#include <vector>
#include "Core/Geometry/LevelSetHermiteGeometry.h"

// Prescribed free surface charge q=Dn_positive-Dn_negative. Derivatives are
// with respect to arclength along t=(-ny,nx), for n negative -> positive.
struct SurfaceChargeJet
{
    double density = 0.0;
    double firstDerivative = 0.0;
    double secondDerivative = 0.0;
};

// Two views of the SAME caller-owned surface charge field at one time:
// cellIntegratedCharge[k] is integral(q ds), or integral(q 2*pi*r ds) in
// axisymmetric mode, in solver cell-index order. Jets always differentiate
// with respect to meridional arclength, never cumulative revolved area;
// evaluate supplies its local jet on the requested actual interface frame.
// The caller must supply consistent quadrature/geometry, stable data and a
// deterministic, side-effect-free callback throughout solve(). This interface
// neither constructs charge panels nor advances/remaps their charge.
struct SurfaceChargeSource
{
    std::vector<double> cellIntegratedCharge;
    std::function<SurfaceChargeJet(const LevelSetHermiteGeometry::Frame&)> evaluate;
    bool axisymmetric = false; // Must match the electric operator and geometry.
};
