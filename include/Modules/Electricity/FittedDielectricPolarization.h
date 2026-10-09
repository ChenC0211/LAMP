#pragma once

#include <vector>
#include "Modules/Electricity/FittedElectricPotential.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Solvers/Poisson/PolarizationPoissonFVSolverCellCentered.h"

class FaceSegmentGraph;

// Saved quantities of the charge-free cell-harmonic FV model, NOT sharp
// two-sided interface traces. No contact-line/wall traction is represented.
struct FittedDielectricPolarization
{
    // Retain the geometry and material that generated the coefficients. These
    // travel together to a receiving discretization; no Hydro operator is used.
    NativeLeafDomain geometry;
    std::vector<double> permittivityByCell;
    std::vector<double> displacementByNative, tangentByNative, coefficientByNative;
    std::vector<double> fieldX, fieldY; // D-area reconstruction / actual cell epsilon.
    double chargeFreeResidual = 0.0, boundaryFluxRelativeClosure = 0.0;
    explicit FittedDielectricPolarization(const QuadtreeMesh& mesh);
};

// Electric-code-unit coefficient -0.5*(Et^2 + Dn^2/(epsLow*epsHigh)) for a
// receiving scalar-gradient load. Uses actual solved fitted BC fluxes. Native
// boundary and equal-permittivity faces have zero coefficients. The receiver
// owns the gradient, dual measure, constraints and mechanical normalization.
// q=0 only: checks actual cell Gauss residual, does not infer free charge or
// add qE. Conductive time advancement must supply its own charge closure.
FittedDielectricPolarization buildFittedDielectricPolarization(const QuadtreeMesh& mesh,
    const NativeLeafDomain& domain, const FaceSegmentGraph& graph,
    const FittedPotentialInput& potential, const std::vector<double>& epsilonCell,
    const std::vector<double>& epsilonSegment,
    const std::vector<PolarizationPoissonFVSolverCellCentered::FittedBoundary>& boundary,
    double relativeTolerance);
