#pragma once

#include "Modules/Interface/HermiteGalsState.h"
#include "Core/Geometry/LevelSetInterfacePanels.h"
#include <memory>
#include <vector>

namespace Quantum::Interface {
struct HermiteGalsPhaseRegionConfig {
    double relativeLedgerTolerance=1e-10;
    double relativeMomentEstimateTolerance=1e-10;
    std::vector<Vector2D> boundaryContacts;
    std::vector<LevelSetInterfacePanels::SymmetryPlane> symmetryPlanes;
};
struct HermiteGalsPhaseRegionLeafMeasure {
    int cell=-1;
    Vector2D lower,upper;
    long double capacity=0,liquid=0,gas=0; // same query, canonical smaller phase/complement
};
struct HermiteGalsPhaseRegionMeasure {
    HermiteGalsIdentity sourceIdentity;
    Vector2D lower,upper;
    long double capacity=0,liquid=0,gas=0;
    long double rawAxialLiquid=0,rawRadialLiquid=0;
    long double rawAxialGas=0,rawRadialGas=0;
    long double momentErrorEstimate=0,maxRelativeLedgerDefect=0;
    std::size_t intersectedLeaves=0,cutLeaves=0,integratedArcs=0,rootRefinements=0;
    std::vector<HermiteGalsPhaseRegionLeafMeasure> leaves;
};
// Prepared complete topology of the SAME immutable native Q4 scalar. Does not
// initialize or update M. A query clips actual source leaves; it never fits a
// replacement polynomial. Unsupported cuts/work limits throw transactionally.
// FirstJet uses complete native trace discovery: side-attached geometry has
// certified pure leaf interiors and is not passed to the regular arc consumer.
// ConservativeVolume retains the original legacy discovery contract.
class HermiteGalsPhaseRegionSource {
public:
    HermiteGalsPhaseRegionSource(const QuadtreeMesh&,const HermiteGalsState&,
        const HermiteGalsPhaseRegionConfig& = {});
    const HermiteGalsState& source() const;
    std::size_t components() const;
    std::size_t nativeSideGraphs() const;
    void validateMesh(const QuadtreeMesh&) const;
    HermiteGalsPhaseRegionMeasure measure(const QuadtreeMesh&,Vector2D lower,Vector2D upper) const;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
    const Data& checked() const;
    HermiteGalsPhaseRegionMeasure measureValidatedMesh(const QuadtreeMesh&,Vector2D lower,Vector2D upper) const;
    friend class HermiteGalsInventoryRegionSource; // one checked constructor batch; no public unchecked query
};
}
