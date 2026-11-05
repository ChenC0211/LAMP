#pragma once

#include "Modules/Interface/HermiteGalsTrace.h"
#include <vector>

namespace LAMP::Interface
{
class HermiteGalsVertexTransportResult;

// Transport a nonempty, strictly ascending selection of existing independent
// master vertices. Coordinates come from mesh; phi and grad(phi) come from the
// SAME retained source snapshot at B_h, followed by J_Bh^T grad(phi).
// This is a vertex-stage result, not a new closed geometry or accepted step.
HermiteGalsVertexTransportResult transportHermiteGalsVertices(
    const QuadtreeMesh& mesh, const HermiteGalsState& source,
    const std::vector<int>& masterVertices, const LevelSetBoundaryTrace::Domain& domain,
    const HermiteGalsVelocityView& previousVelocity, const HermiteGalsVelocityView& arrivalVelocity);

class HermiteGalsVertexTransportResult
{
public:
    const std::vector<HermiteGalsVertexJet>& vertices() const noexcept { return vertices_; }
    HermiteGalsIdentity sourceIdentity() const noexcept { return source_; }
    HermiteGalsTimeLayer sourceTimeLayer() const noexcept { return sourceTime_; }
    // Identity of the velocity's coupling geometry, NOT a newly published
    // Interface geometry. Closure and geometry publication follow in I05-I11.
    HermiteGalsIdentity arrivalVelocityGeometryIdentity() const noexcept { return arrival_; }
    HermiteGalsTimeLayer arrivalTimeLayer() const noexcept { return arrivalTime_; }
private:
    HermiteGalsVertexTransportResult(HermiteGalsIdentity source, HermiteGalsTimeLayer sourceTime,
        const HermiteGalsVelocityView& arrival, std::vector<HermiteGalsVertexJet> vertices);
    HermiteGalsIdentity source_, arrival_;
    HermiteGalsTimeLayer sourceTime_, arrivalTime_;
    std::vector<HermiteGalsVertexJet> vertices_;
    friend HermiteGalsVertexTransportResult transportHermiteGalsVertices(
        const QuadtreeMesh&, const HermiteGalsState&, const std::vector<int>&,
        const LevelSetBoundaryTrace::Domain&, const HermiteGalsVelocityView&, const HermiteGalsVelocityView&);
};
} // namespace LAMP::Interface
