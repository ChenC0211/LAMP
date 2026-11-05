#pragma once

#include "Modules/Interface/HermiteGalsState.h"
#include "Core/Geometry/LevelSetInterfacePanels.h"
#include "Core/Geometry/LevelSetMappedInterface.h"
#include <array>
#include <cstddef>
#include <exception>
#include <memory>
#include <vector>

namespace LAMP::Interface {
class HermiteGalsPublication;
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
// A checked batch keeps later query failures in their original ordered slot.
// Query zero is evaluated by the caller and its failure is thrown immediately.
struct HermiteGalsPhaseRegionBatchMeasure {
    HermiteGalsPhaseRegionMeasure measure;
    std::exception_ptr failure;
};
// Actual cut-integration work accumulated by this source. Worker durations
// are nested in query/batch wall time; pure phases, cache hits and waits add none.
struct HermiteGalsPhaseRegionWorkProfile {
    double pairWorkerSeconds=0,fallbackWorkerSeconds=0;
    std::size_t pairCalls=0,fallbackCalls=0;
    double pairRegionGaugeWorkerSeconds=0,pairCutPreparationWorkerSeconds=0,pairCoreWorkerSeconds=0,
        pairRetryWorkerSeconds=0,pairAssemblyWorkerSeconds=0;
    std::size_t pairCutConstructions=0,pairCoreAttempts=0,pairFaceCutRefinements=0;
    // Nested Core worker work; its children do not add to the pair/caller wall.
    LevelSetMappedInterface::NativeMomentWorkProfile pairCoreWork;
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
    // Retain only an actual immutable publication's complete certification.
    // Config declarations must match its curves' contact/plane union exactly;
    // the original numerical region tolerances still control every query.
    HermiteGalsPhaseRegionSource(const QuadtreeMesh&,const HermiteGalsPublication&,
        const HermiteGalsPhaseRegionConfig& = {});
    const HermiteGalsState& source() const;
    std::size_t components() const;
    std::size_t nativeSideGraphs() const;
    // Read after synchronous queries complete; sums actual native-cell slots
    // in native index order, including work done before an integration throws.
    HermiteGalsPhaseRegionWorkProfile workProfile() const;
    void validateMesh(const QuadtreeMesh&) const;
    HermiteGalsPhaseRegionMeasure measure(const QuadtreeMesh&,Vector2D lower,Vector2D upper) const;
    // One mesh-binding check for this immutable-source batch; each query
    // retains its actual-thread FP check and the original measure arithmetic.
    std::vector<HermiteGalsPhaseRegionBatchMeasure> measureBatch(const QuadtreeMesh&,
        const std::vector<std::array<Vector2D,2>>& orderedBoxes) const;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
    const Data& checked() const;
    HermiteGalsPhaseRegionSource(const QuadtreeMesh&,const HermiteGalsState&,
        const HermiteGalsPhaseRegionConfig&,const HermiteGalsPublication*);
    HermiteGalsPhaseRegionMeasure measureValidatedMesh(const QuadtreeMesh&,Vector2D lower,Vector2D upper) const;
    friend class HermiteGalsInventoryRegionSource; // one checked constructor batch; no public unchecked query
};
}
