#pragma once

#include <string>
#include "Modules/Hydrodynamics/HydroBoundary.h"
#include "Modules/Hydrodynamics/HydroCompleteConservation.h"
#include "Modules/Hydrodynamics/HydroMidpointStep.h"

// Fixed-step flow on a rectangular domain (closed unless patches are supplied).
// A named level set may supply
// two material values and be transported by the returned flow velocity.
// All physical values are code units, independent of electrospray constants.
// Main-clock and output settings are not module physics; the IO adapter
// exposes legacy run and Interface input separately for the owning App.
struct HydrodynamicsConfig
{
    std::string coordinateSystem = "planar2D";
    std::string timeScheme = "cnAb2";
    Quantum::Hydrodynamics::HydroMidpointIterationOptions midpoint;
    int momentMaximumIterations=24;
    std::string momentConvection="central";
    // Explicit coupled unresolved-neck model, disabled unless selected.
    double momentNeckRadiusCells=0;
    double density = 1.0;
    double viscosity = 0.01;
    std::string materialModel = "constant";
    std::string interfaceLevelSet;
    double densityNegative = 1.0;
    double densityPositive = 1.0;
    double viscosityNegative = 0.01;
    double viscosityPositive = 0.01;
    double interfaceThickness = 0.0;
    // Traction on the named fluid interface, in code units. The normal points
    // from negative to positive phi; tangent = (-normal.y, normal.x).
    double surfaceTension = 0.0;
    // Explicit choice shared by standalone Hydro and coupled applications.
    // The default preserves the existing contour capillary discretization.
    std::string interfaceForceGeometry = "legacyContour";
    // H used by the shared Hermite force context. Regularized uses the existing
    // positive interfaceThickness as its half-width; no extra physical scale.
    std::string interfaceForceLayer = "geometricPL";
    double interfaceNormalTraction = 0.0;
    double interfaceTangentialTraction = 0.0;
    double accelerationX = 0.0;
    double accelerationY = 0.0;
    std::vector<HydroBoundaryPatch> boundaryPatches;
    std::string initialVelocity = "zero";
    // Prescribed p_negative - p_positive at t=0, in code pressure units.
    // Closed conservative cases use p_i = initialPressureJump * C_i.
    // This is initial data only, never a replacement for computed curvature.
    double initialPressureJump = 0.0;
    double velocityAmplitude = 0.05;
    double linearResidualTolerance = 1e-8;
    HydroCompleteConservation conservation;
};
