#pragma once

#include "Modules/Interface/HermiteGalsState.h"
#include "Core/Geometry/ActiveQuarticHermiteSkeleton.h"
#include "Core/Geometry/LevelSetInterfacePanels.h"
#include <cstddef>
#include <vector>

namespace Quantum::Interface
{
struct HermiteGalsAmrClosureReport
{
    std::size_t masterVertices=0,hangingVertices=0,masterEdges=0,constrainedEdges=0;
    std::size_t lineComponents=0,compatibilityNonzeros=0;
    std::size_t symmetryConstrainedVertices=0;
    double algebraicResidual=0,compatibilityResidual=0,maxRelativeDirectionalDisagreement=0;
};
class HermiteGalsAmrClosureResult;

// Same-revision compact closure + exact hanging constraints. No geometry
// publication, cell interior, volume update or new independent transported jet.
// Declared outer reflection planes constrain derived mixed data, and require
// exactly zero primary normal first jets. The RZ axis remains automatic.
HermiteGalsAmrClosureResult closeHermiteGalsMixedAmr(
    const QuadtreeMesh& mesh,const HermiteGalsCandidate& candidate,
    const std::vector<LevelSetInterfacePanels::SymmetryPlane>& symmetryPlanes={});

class HermiteGalsAmrClosureResult
{
public:
    const std::vector<double>& vertexMixedDerivatives() const noexcept {return mixed_;}
    const Quantum::Geometry::ActiveQuarticHermiteSkeleton& skeleton() const noexcept {return skeleton_;}
    HermiteGalsIdentity primaryIdentity() const noexcept {return primary_;}
    HermiteGalsTimeLayer timeLayer() const noexcept {return time_;}
    const HermiteGalsAmrClosureReport& report() const noexcept {return report_;}
private:
    HermiteGalsAmrClosureResult(HermiteGalsIdentity primary,HermiteGalsTimeLayer time,std::vector<double> mixed,
        Quantum::Geometry::ActiveQuarticHermiteSkeleton skeleton,HermiteGalsAmrClosureReport report);
    HermiteGalsIdentity primary_;
    HermiteGalsTimeLayer time_;
    std::vector<double> mixed_;
    Quantum::Geometry::ActiveQuarticHermiteSkeleton skeleton_;
    HermiteGalsAmrClosureReport report_;
    friend HermiteGalsAmrClosureResult closeHermiteGalsMixedAmr(const QuadtreeMesh&,const HermiteGalsCandidate&,
        const std::vector<LevelSetInterfacePanels::SymmetryPlane>&);
};
} // namespace Quantum::Interface
