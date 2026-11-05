#pragma once

#include "Modules/Interface/HermiteGalsPhaseVolume.h"
#include <memory>
#include <vector>

namespace LAMP::Interface
{
struct HermiteGalsRefinementParentBalance
{
    int parent = -1;
    std::size_t children = 0;
    long double capacity = 0, originalInventory = 0;
    long double childGeometry = 0, childInventory = 0, residual = 0;
};
struct HermiteGalsRefinementInventoryReport
{
    std::size_t refinedParents = 0, unchangedParents = 0, childCells = 0;
    long double originalInventory = 0, transferredInventory = 0;
    long double maxRelativeSourceGeometryDefect = 0;
    long double maxRelativeParentReferenceDifference = 0;
    long double maxRelativeParentResidual = 0, maxRelativeChildGeometryDefect = 0;
};

// Prepared, immutable inventory and exact scalar for a pure refinement.
// This does not produce an Interface state/publication or accept a time step.
class HermiteGalsRefinementInventory
{
public:
    HermiteGalsIdentity sourceIdentity() const;
    HermiteGalsTimeLayer timeLayer() const;
    const HermiteGalsState& source() const;
    const LevelSetHermiteGeometry::Snapshot& geometry() const;
    const std::vector<HermiteGalsCellVolume>& phaseVolumes() const;
    const std::vector<HermiteGalsPhaseVolumeCell>& cells() const;
    const std::vector<int>& parentCells() const;
    const std::vector<HermiteGalsRefinementParentBalance>& parentBalances() const;
    const HermiteGalsRefinementInventoryReport& report() const;
    void validateMeshes(const QuadtreeMesh& sourceMesh,const QuadtreeMesh& targetMesh) const;
private:
    struct Data;
    explicit HermiteGalsRefinementInventory(std::shared_ptr<const Data> data);
    const Data& checked() const;
    std::shared_ptr<const Data> data_;
    friend HermiteGalsRefinementInventory refineHermiteGalsInventory(const QuadtreeMesh&,
        const HermiteGalsState&,const QuadtreeMesh&,const HermiteGalsPhaseVolumeConfig&);
};

// M remains an independent input. A parent-local, bounded roundoff allocation
// preserves its original target; incompatible M/geometry or storage rejects.
// No coarsening, partial support, new closure, clipping or analytic reset.
HermiteGalsRefinementInventory refineHermiteGalsInventory(const QuadtreeMesh& sourceMesh,
    const HermiteGalsState& source,const QuadtreeMesh& targetMesh,
    const HermiteGalsPhaseVolumeConfig& config = {});
}
