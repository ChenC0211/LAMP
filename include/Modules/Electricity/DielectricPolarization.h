#pragma once

#include "Core/Mesh/NativeMeshIdentity.h"
#include <vector>

// Immutable publication of the diffuse planar solve's material and physical
// coefficient. Hydro owns the receiving scalar gradient and velocity rows.
// Native ordering is x edges followed by y edges; boundary coefficients zero.
struct DielectricPolarization
{
    NativeMeshIdentity meshIdentity;
    std::vector<double> permittivityByCell, coefficientByNative;
    explicit DielectricPolarization(const QuadtreeMesh& mesh) : meshIdentity(mesh) {}
};
