#pragma once

#include <functional>
#include <vector>
#include "Core/Base/Enumeration.h"
#include "Core/Base/Vector2D.h"

// Explicit scalar boundary data, independent of Hydro/material/contact-angle
// models. Callbacks own their data and must be deterministic and side-effect
// free. The caller supplies physical time; this object owns no clock.
struct LevelSetBoundaryConditions
{
    struct InflowPatch
    {
        Boundary side = Boundary::non;
        double start = 0.0, end = 0.0; // Physical tangential coordinates, closed interval.
        std::function<double(Vector2D, double)> value; // Boundary point, physical time -> phi.
    };
    std::vector<InflowPatch> inflow;
    // Exact mesh vertices on a physical boundary. This imposes scalar zero,
    // not zero material velocity or a contact angle. No coordinate snapping.
    std::vector<Vector2D> fixedZeroVertices;

    bool empty() const noexcept { return inflow.empty() && fixedZeroVertices.empty(); }
};

// Optional complete, old/new spatial velocity interpolants for explicit
// scalar-boundary advancement. Callbacks own/capture stable input data for
// the duration of one update; geometry never retains them or chooses their
// physical wall law. They avoid resampling a discontinuous boundary trace
// into a single vertex value shared by two different boundary patches.
struct LevelSetAdvectionVelocity
{
    std::function<Vector2D(Vector2D)> current, previous;
};

// Optional retained interpolation of phi^n, indexed by the actual level set.
// The caller must bind the SAME old nodal scalar used by the geometry/loads,
// not a reference shape or a corrected phase volume. It is consumed only
// while staging an update and is never retained. Boundary hits still use
// their prescribed scalar at the actual hit time, independently of this data.
struct LevelSetAdvectionScalar
{
    std::function<double(int, Vector2D)> previous;
};
