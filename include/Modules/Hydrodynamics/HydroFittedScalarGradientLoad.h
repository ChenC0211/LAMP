#pragma once

#include "Core/Field/VEdgeField.h"
#include "Core/Field/HEdgeField.h"
#include <vector>

class HydroFittedDomain;
class NativeLeafDomain;
class QuadtreeMesh;

// Coefficient times the actual pressure-compatible scalar gradient on free
// active-active faces. This port supplies no boundary traction. It retains the
// two-cell identity on equal-level faces and the existing virtual-child
// pressure stencil at coarse/fine faces. No second gradient operator is used.
struct HydroFittedScalarGradientLoad
{
    struct Face
    {
        int nativeRow = -1, lowCell = -1, highCell = -1;
        double gradient = 0.0, dualVolume = 0.0, forceDensity = 0.0;
    };
    VEdgeField<double> forceX;
    HEdgeField<double> forceY;
    std::vector<Face> faces;
    double maxGradientIdentityError = 0.0;
    double maxAbsComponent = 0.0, integralX = 0.0, integralY = 0.0;
    explicit HydroFittedScalarGradientLoad(const QuadtreeMesh& mesh);
};

// sourceGeometry describes the native storage/material domain of both arrays.
// Validate exact mesh identity, selected-cell maps and coordinate measure
// before using the current Hydro pressure operator. scale converts the supplied
// coefficient*gradient to the receiving force-density units, in that order.
HydroFittedScalarGradientLoad buildFittedScalarGradientLoad(const QuadtreeMesh& mesh,
    const HydroFittedDomain& domain, const NativeLeafDomain& sourceGeometry,
    const std::vector<double>& scalarByCell, const std::vector<double>& coefficientByNative,
    double scale);
