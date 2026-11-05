#pragma once

#include "Modules/Interface/HermiteGalsState.h"
#include "Core/Geometry/ActiveQuarticHermiteSkeleton.h"
#include "Core/Geometry/LevelSetInterfacePanels.h"
#include <cstddef>
#include <vector>

namespace LAMP::Interface
{
struct HermiteGalsAmrClosureReport
{
    std::size_t masterVertices=0,hangingVertices=0,masterEdges=0,constrainedEdges=0;
    std::size_t lineComponents=0,compatibilityNonzeros=0;
    std::size_t symmetryConstrainedVertices=0;
    double algebraicResidual=0,compatibilityResidual=0,maxRelativeDirectionalDisagreement=0;
};
class HermiteGalsAmrClosureResult;

// Sparse sensitivity of the actual compact mixed closure, including its
// independent master edge normal derivatives. Ordinals index the retained
// primary vertices (VertexDx/Dy) or edges (EdgeNormal), never mesh labels.
enum class HermiteGalsPrimaryDerivative { VertexDx,VertexDy,EdgeNormal };
struct HermiteGalsPrimaryDerivativeWeight
{
    HermiteGalsPrimaryDerivative component=HermiteGalsPrimaryDerivative::VertexDx;
    std::size_t primaryOrdinal=0;
    double weight=0;
};
class HermiteGalsAmrClosureSensitivity
{
public:
    // Complete uniform nonhanging XY FirstJet only. Same DirectDirection
    // five-point source and exact declared-symmetry constraints as closure.
    static HermiteGalsAmrClosureSensitivity buildUniformFirstJet(const QuadtreeMesh&,
        const HermiteGalsState&,
        const std::vector<LevelSetInterfacePanels::SymmetryPlane>& symmetryPlanes={});
    // Row slot is the retained primary vertex ordinal.
    const std::vector<std::vector<HermiteGalsPrimaryDerivativeWeight>>& rows() const noexcept {return rows_;}
    const std::vector<LevelSetInterfacePanels::SymmetryPlane>& symmetryPlanes() const noexcept {return symmetry_;}
    HermiteGalsIdentity primaryIdentity() const noexcept {return primary_;}
    HermiteGalsTimeLayer timeLayer() const noexcept {return time_;}
    void validateFor(const QuadtreeMesh&,const HermiteGalsState&) const;
    void validateFor(const QuadtreeMesh&,HermiteGalsIdentity,HermiteGalsTimeLayer) const;
private:
    HermiteGalsAmrClosureSensitivity(HermiteGalsState,
        std::vector<std::vector<HermiteGalsPrimaryDerivativeWeight>>,
        std::vector<LevelSetInterfacePanels::SymmetryPlane>);
    HermiteGalsState source_;
    HermiteGalsIdentity primary_;
    HermiteGalsTimeLayer time_;
    std::vector<std::vector<HermiteGalsPrimaryDerivativeWeight>> rows_;
    std::vector<LevelSetInterfacePanels::SymmetryPlane> symmetry_;
};

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
    const LAMP::Geometry::ActiveQuarticHermiteSkeleton& skeleton() const noexcept {return skeleton_;}
    HermiteGalsIdentity primaryIdentity() const noexcept {return primary_;}
    HermiteGalsTimeLayer timeLayer() const noexcept {return time_;}
    const HermiteGalsAmrClosureReport& report() const noexcept {return report_;}
private:
    HermiteGalsAmrClosureResult(HermiteGalsIdentity primary,HermiteGalsTimeLayer time,std::vector<double> mixed,
        LAMP::Geometry::ActiveQuarticHermiteSkeleton skeleton,HermiteGalsAmrClosureReport report);
    HermiteGalsIdentity primary_;
    HermiteGalsTimeLayer time_;
    std::vector<double> mixed_;
    LAMP::Geometry::ActiveQuarticHermiteSkeleton skeleton_;
    HermiteGalsAmrClosureReport report_;
    friend HermiteGalsAmrClosureResult closeHermiteGalsMixedAmr(const QuadtreeMesh&,const HermiteGalsCandidate&,
        const std::vector<LevelSetInterfacePanels::SymmetryPlane>&);
};
} // namespace LAMP::Interface
