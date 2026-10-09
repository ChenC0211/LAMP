#pragma once

#include <memory>
#include <utility>
#include <vector>
#include "Core/Base/Vector2D.h"

class QuadtreeMesh;
template<typename T> class VEdgeField;
template<typename T> class HEdgeField;
template<typename T> class VertField;

// Caller-supplied tangential wall velocities, ordered by increasing physical
// grid coordinate. uxSouth/North have Nx+1 values on vertex x coordinates;
// uyWest/East have Ny+1 values on vertex y coordinates. Corner component
// values belong to these same arrays; Core supplies no slip/no-slip law or
// second, independently inferred corner value.
struct StaggeredVelocityBoundaryTrace
{
    std::vector<double> uxSouth, uxNorth;
    std::vector<double> uyWest, uyEast;
};

// Shared native-face -> vertex conversion for level-set and material tracing.
// Default preserves the legacy least-squares conversion. The explicit RZ mode
// enforces odd radial/even axial velocity at the actual r=0 boundary, after
// validating the native zero radial velocity; all non-axis values are unchanged.
std::pair<std::unique_ptr<VertField<double>>, std::unique_ptr<VertField<double>>>
staggeredToVertexVelocity(const QuadtreeMesh& mesh,
    const VEdgeField<double>& ux, const HEdgeField<double>& uy,
    bool radialAxisSymmetry = false);

// Explicit new path: sample the SAME retained native/boundary interpolant at
// each actual mesh vertex, instead of using the legacy least-squares fit.
// Both components cover the entire closed domain. RZ requires every native
// axis radial value and both radial corner traces to be exactly zero; the
// supplied uyWest is authoritative, not replaced by an inferred axis value.
std::pair<std::unique_ptr<VertField<double>>, std::unique_ptr<VertField<double>>>
staggeredToVertexVelocity(const QuadtreeMesh& mesh,
    const VEdgeField<double>& ux, const HEdgeField<double>& uy,
    const StaggeredVelocityBoundaryTrace& boundary, bool radialAxisSymmetry = false);

// Retained native-face velocity on a complete, fixed uniform Cartesian mesh.
// Each component is bilinearly sampled on its own staggered coordinate grid.
// Without explicit boundary traces, queries must lie in BOTH native grids'
// closed rectangles. Trace overloads add only the missing wall strips. There
// is no wall clipping, ghost extrapolation, phase model or inferred wall law.
class StaggeredVelocitySampler
{
public:
    StaggeredVelocitySampler(const QuadtreeMesh& mesh,
        const VEdgeField<double>& ux, const HEdgeField<double>& uy);
    // Copies 0.5*(old+new) on each native face before spatial interpolation.
    StaggeredVelocitySampler(const QuadtreeMesh& mesh,
        const VEdgeField<double>& uxOld, const HEdgeField<double>& uyOld,
        const VEdgeField<double>& uxNew, const HEdgeField<double>& uyNew);
    // Explicit full-closed-domain extension through real supplied boundary
    // traces. Only the missing half-cell strips are added; native-grid
    // interpolation keeps its original expression and actual coordinates.
    StaggeredVelocitySampler(const QuadtreeMesh& mesh,
        const VEdgeField<double>& ux, const HEdgeField<double>& uy,
        const StaggeredVelocityBoundaryTrace& boundary, bool radialAxisSymmetry = false);
    // Both the native faces and supplied traces are averaged before spatial
    // interpolation. RZ checks EACH time level before any cancellation can
    // hide a nonzero axis radial velocity. All trace arrays are owned copies.
    StaggeredVelocitySampler(const QuadtreeMesh& mesh,
        const VEdgeField<double>& uxOld, const HEdgeField<double>& uyOld,
        const VEdgeField<double>& uxNew, const HEdgeField<double>& uyNew,
        const StaggeredVelocityBoundaryTrace& boundaryOld,
        const StaggeredVelocityBoundaryTrace& boundaryNew, bool radialAxisSymmetry = false);
    ~StaggeredVelocitySampler();
    StaggeredVelocitySampler(StaggeredVelocitySampler&&) noexcept;
    StaggeredVelocitySampler& operator=(StaggeredVelocitySampler&&) noexcept;
    StaggeredVelocitySampler(const StaggeredVelocitySampler&) = delete;
    StaggeredVelocitySampler& operator=(const StaggeredVelocitySampler&) = delete;

    // Snapshot owns values and coordinate/index grids, not borrowed fields.
    // A moved-from sampler and a non-finite/out-of-domain point are rejected.
    Vector2D sample(Vector2D point) const;

private:
    struct Impl;
    std::unique_ptr<const Impl> impl_;
};
