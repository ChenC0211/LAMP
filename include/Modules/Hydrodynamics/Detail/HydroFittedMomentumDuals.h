#pragma once

#include <array>
#include <vector>
#include "Core/Base/Vector2D.h"

class QuadtreeMesh;
class HydroFittedDomain;
class NativeLeafDomain;

namespace HydroFittedDetail
{
// Geometry of the physical staggered momentum control volumes. This is an
// immutable connectivity view of existing native faces, NOT another mesh or
// another set of velocity unknowns. In particular a coarse/fine dual uses
// actual owner centres, unlike the legacy virtual-fine pressure metric.
// Do not substitute these volumes into only one term of that legacy operator.
class MomentumDuals
{
public:
    struct Volume
    {
        bool active = false;
        Vector2D lower{0.0, 0.0}, upper{0.0, 0.0};
        double measure = 0.0;
        // W, E, S, N; each common segment appears in both adjacent volumes.
        std::array<std::vector<int>, 4> sides;
    };
    struct Segment
    {
        // Positive-coordinate normal. -1 denotes physical exterior, never a
        // solid/ghost unknown. Constrained wall/axis native rows are retained.
        int negativeNative = -1, positiveNative = -1;
        int normalAxis = 0; // 0: r/x, 1: z/y
        Vector2D begin{0.0, 0.0}, end{0.0, 0.0};
        double measure = 0.0; // Area in RZ, length in planar geometry
        // Exactly one source: an existing primal face, or a plane within an
        // active primal cell. Exterior boundary data remain owned by domain.
        int primalNative = -1;
        int ownerCell = -1;
    };

    static MomentumDuals build(const QuadtreeMesh& mesh, const HydroFittedDomain& domain);
    // Same physical dual geometry for Interface-owned native domains; no
    // pressure boundary conditions or extra velocity unknowns are invented.
    static MomentumDuals build(const QuadtreeMesh& mesh, const NativeLeafDomain& domain);
    const std::vector<Volume>& volumes() const noexcept { return volumes_; }
    const std::vector<Segment>& segments() const noexcept { return segments_; }

private:
    std::vector<Volume> volumes_;
    std::vector<Segment> segments_;
};
}
