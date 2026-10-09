#pragma once

#include "Core/IO/TecplotWriter.h"
#include <vector>

class QuadtreeMesh;
class HydroFittedDomain;
class Hydrodynamics;

namespace HydroFittedOutput
{
struct Scales
{
    double length = 0.0, time = 0.0, density = 0.0;
    double viscosity = 0.0, pressure = 0.0, velocity = 0.0;
};

struct Frame
{
    TecplotWriter::Zone zone;
    // Zone element/node order -> original native indices. Apps may append
    // matching cell/nodal electric or scalar fields without changing topology.
    std::vector<int> nativeCellIndices, nativeVertexIndices;
};

// Active-fluid actual leaf quadrilaterals only; solid payload is not read.
// Physical rho/mu/p and stored FV divergence remain cell-centered. Velocity
// is the existing RZFluxCompatible reconstruction sampled at each ACTUAL
// cell center, not a face value or an independently divergence-free field.
// This explicit RZ output path inherits that sampler's supported topology.
// Coordinates/fields/time are SI; timeSI is already dimensional. Only builds
// an owned Zone, with no file I/O, dictionary access, field or mesh mutation.
Frame cells(const QuadtreeMesh& mesh, const HydroFittedDomain& domain,
    const Hydrodynamics& hydro, const Scales& scales, int step,
    double timeSI, int strandId = 1);
}
