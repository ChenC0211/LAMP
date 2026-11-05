#pragma once

#include "Modules/Interface/HermiteGalsTrace.h"
#include <array>
#include <vector>

namespace LAMP::Interface
{
class HermiteGalsEdgeTransportResult;

// Nonempty, lexicographically ascending native endpoint pairs. Each pair must
// be an existing source master and a full native leaf side, not a constrained
// subedge of a coarse side. Endpoints follow +x/+y; dn uses +y/+x respectively.
// The full old gradient is pulled back BEFORE its normal component is taken.
HermiteGalsEdgeTransportResult transportHermiteGalsEdges(
    const QuadtreeMesh& mesh,const HermiteGalsState& source,
    const std::vector<std::array<int,2>>& masterEdges,const LevelSetBoundaryTrace::Domain& domain,
    const HermiteGalsVelocityView& previousVelocity,const HermiteGalsVelocityView& arrivalVelocity);

class HermiteGalsEdgeTransportResult
{
public:
    const std::vector<HermiteGalsEdgeJet>& edges() const noexcept { return edges_; }
    HermiteGalsIdentity sourceIdentity() const noexcept { return source_; }
    HermiteGalsTimeLayer sourceTimeLayer() const noexcept { return sourceTime_; }
    // Coupling geometry of the supplied velocity, not a newly published state.
    HermiteGalsIdentity arrivalVelocityGeometryIdentity() const noexcept { return arrival_; }
    HermiteGalsTimeLayer arrivalTimeLayer() const noexcept { return arrivalTime_; }
private:
    HermiteGalsEdgeTransportResult(HermiteGalsIdentity source,HermiteGalsTimeLayer sourceTime,
        const HermiteGalsVelocityView& arrival,std::vector<HermiteGalsEdgeJet> edges);
    HermiteGalsIdentity source_,arrival_;
    HermiteGalsTimeLayer sourceTime_,arrivalTime_;
    std::vector<HermiteGalsEdgeJet> edges_;
    friend HermiteGalsEdgeTransportResult transportHermiteGalsEdges(const QuadtreeMesh&,
        const HermiteGalsState&,const std::vector<std::array<int,2>>&,const LevelSetBoundaryTrace::Domain&,
        const HermiteGalsVelocityView&,const HermiteGalsVelocityView&);
};
} // namespace LAMP::Interface
