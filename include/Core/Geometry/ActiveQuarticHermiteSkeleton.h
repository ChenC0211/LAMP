#pragma once

#include "Core/Geometry/ActiveQuarticHermitePatch.h"
#include "Core/Geometry/ActiveHermiteSensitivity.h"
#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

class QuadtreeMesh;

namespace LAMP::Geometry
{
// Shared first/mixed jet inputs only. No cell centers, sectors, geometry,
// transport, phase inventory or physical clock. Interface owns those roles.
struct ActiveQuarticHermiteSkeletonInput
{
    std::vector<std::uint8_t> activeCells; // full native cell-index length
    std::vector<ActiveQuarticHermitePatch::NodalDof> nodes; // full vertex length
    std::vector<ActiveQuarticHermitePatch::EdgeMidpoint> edgeMidpoints;
    // Only regular used nodes and maximal active sides are independent.
    // Hanging/fine payloads are ignored (may be NaN); unused nodes are ignored.
    // All coarse traces supporting an active hanging vertex must be active.
};

// Immutable algebraic skeleton owner, using the SAME exact constraint engine
// as native Q4 geometry. It does not construct any cell polynomial/snapshot.
class ActiveQuarticHermiteSkeleton
{
public:
    using EdgeKey=std::array<int,2>;
    using Jet=std::array<double,4>; // physical phi, dx, dy, dxy
    struct Vertex
    {
        int index=-1;
        EdgeKey masterSide{{-1,-1}}; // invalid pair means independent vertex
    };
    struct Edge
    {
        EdgeKey vertices{{-1,-1}}; // actual active side, positive x/y tangent
        EdgeKey parentSide{{-1,-1}}; // invalid pair means independent midpoint
    };
    explicit ActiveQuarticHermiteSkeleton(const QuadtreeMesh& mesh,
        const ActiveQuarticHermiteSkeletonInput& input,bool radialAxisSymmetry=false);

    const std::vector<Vertex>& vertices() const;
    const std::vector<Edge>& edges() const;
    std::size_t independentDofCount() const;
    void validateMesh(const QuadtreeMesh& mesh) const;

    // Returned doubles are rounded query results, not new independent DOFs.
    Jet nodeJet(int nativeVertex) const;
    // Same retained constraint DAG, pulled back to source master inputs.
    // component: phi/dx/dy/dxy = 0..3; edge indices refer to input order.
    std::vector<ActiveQuarticHermiteSensitivity> nodeSensitivity(int nativeVertex,int component) const;
    std::array<double,2> edgeMidpoint(EdgeKey edge) const; // value, +normal
    // Physical tangential coordinate (x for horizontal / y for vertical).
    // Evaluates retained EXACT traces, including constrained fine sides;
    // never reconstructs from the rounded nodeJet/edgeMidpoint results.
    Jet edgeJet(EdgeKey edge,double coordinate) const;

private:
    struct Impl;
    std::shared_ptr<const Impl> impl_;
    const Impl& implementation() const;
};
} // namespace LAMP::Geometry
