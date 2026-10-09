#pragma once

#include "Core/Field/StaggeredVelocitySampler.h"
#include "Modules/Hydrodynamics/HydroBoundary.h"
#include <functional>

// Explicit fixed-uniform-grid Hydro wall closure. Reconstructs the missing
// tangential velocity nodes using the SAME zero-shear trace as the checked
// NormalVelocity/full-traction PressureOutlet momentum assembly. Native
// values (including solved outlet values) remain authoritative; no pressure
// value is converted to a velocity and no wall value is silently clipped.
// Corner derivatives use the two inward normal faces on the same real side.
// At r=0, retain exact zero radial velocity and the existing even axial trace.
StaggeredVelocityBoundaryTrace hydroZeroShearBoundaryTrace(const QuadtreeMesh& mesh,
    const VEdgeField<double>& ux, const HEdgeField<double>& uy,
    bool radialAxisSymmetry = false);

// With no explicit Hydro patches, preserve the original vertex conversion
// verbatim. The explicit path samples the retained native/wall interpolant
// at vertices; it does not use an independent least-squares wall fit.
std::pair<std::unique_ptr<VertField<double>>, std::unique_ptr<VertField<double>>>
hydroVertexVelocity(const QuadtreeMesh& mesh,
    const VEdgeField<double>& ux, const HEdgeField<double>& uy,
    bool explicitBoundaryPatches, bool radialAxisSymmetry = false);

// Owned full-domain velocity for advection with discontinuous normal patches.
// In the first native-normal grid layer, interpolate from the actual piecewise
// wall constraint to the same-tangent inner native sample; elsewhere preserve
// the native/wall interpolant. Pressure outlets retain their solved velocity.
// NormalVelocity intervals include endpoints: at an inlet/unpatched junction
// the endpoint uses the inlet trace, never an average with the wall. Two
// different explicit normal values at one endpoint are unresolved and throw.
// A fixed scalar contact is not thereby assigned a material trajectory.
std::function<Vector2D(Vector2D)> hydroBoundaryAdvectionVelocity(const QuadtreeMesh& mesh,
    const VEdgeField<double>& ux, const HEdgeField<double>& uy,
    const std::vector<HydroBoundaryPatch>& patches, bool radialAxisSymmetry = false);
