#pragma once

#include "Modules/Interface/HermiteGalsState.h"
#include "Core/Geometry/LevelSetInterfacePanels.h"
#include <cstddef>
#include <memory>
#include <vector>
namespace InterfaceCellMeasure { struct NativeTraceLedger; }
namespace InterfacePhaseGeometry { class MappedFacePartitions; }
class NativeLeafDomain;

namespace Quantum::Interface
{
struct HermiteGalsPhaseVolumeConfig
{
    // Dimensionless, normalized by the relevant native cell / total capacity.
    double relativeLedgerTolerance = 1e-10;
    double relativeMomentEstimateTolerance = 1e-10;
    std::vector<Vector2D> boundaryContacts;
    std::vector<LevelSetInterfacePanels::SymmetryPlane> symmetryPlanes;
};
struct HermiteGalsPhaseVolumeCell
{
    int cell = -1, certifiedPhase = 0; // -1 liquid, +1 gas, 0 cut
    double capacity = 0.0;
    long double rawLiquid = 0.0L; // Green integral, never clipped
};
struct HermiteGalsPhaseVolumeReport
{
    std::size_t liquidCells = 0, gasCells = 0, cutCells = 0;
    std::size_t components = 0, graphs = 0, faceRootIncidences = 0;
    long double totalCapacity = 0, totalRawLiquid = 0, totalInitialLiquid = 0;
    long double relativeBoundaryClosure = 0, relativeMomentEstimate = 0;
    long double maxRelativeCapacityDefect = 0, maxRelativePureDefect = 0;
    long double totalSourceInventory = 0, geometricInventoryDifference = 0;
    long double maxCellRelativeInventoryDifference = 0;
    // The GL8/subdivision difference above is an estimate, not an error bound.
};
class HermiteGalsPhaseVolumeResult
{
public:
    HermiteGalsIdentity sourceIdentity() const;
    HermiteGalsTimeLayer timeLayer() const;
    HermiteGalsGeometryView geometry() const;
    const std::vector<HermiteGalsCellVolume>& phaseVolumes() const;
    const std::vector<HermiteGalsPhaseVolumeCell>& cells() const;
    const HermiteGalsPhaseVolumeReport& report() const;
    // FirstJet only; unavailable on the unchanged historical inventory port.
    const InterfaceCellMeasure::NativeTraceLedger& nativeTraceLedger() const;
    // Only a fresh step-0, time-0 Initial template may seed M. Returns a NEW
    // primary revision without geometry/closure; never updates this source.
    // I11 must subsequently publish geometry for that revision. This is not
    // restart, regrid, time advancement or coupled acceptance.
    HermiteGalsCandidate initializeCandidate(const QuadtreeMesh& mesh) const;
private:
    struct Data;
    explicit HermiteGalsPhaseVolumeResult(std::shared_ptr<const Data> data);
    std::shared_ptr<const Data> data_;
    const Data& checked() const;
    static HermiteGalsPhaseVolumeResult integrate(const QuadtreeMesh&, const HermiteGalsState&,
        const HermiteGalsPhaseVolumeConfig&, const NativeLeafDomain*,
        const InterfacePhaseGeometry::MappedFacePartitions*, const LevelSetInterfacePanels::NativeComponents*);
    friend HermiteGalsPhaseVolumeResult integrateHermiteGalsPhaseVolumes(
        const QuadtreeMesh&, const HermiteGalsState&, const HermiteGalsPhaseVolumeConfig&);
    friend HermiteGalsPhaseVolumeResult integrateHermiteGalsNativeTracePhaseVolumes(const QuadtreeMesh&,
        const HermiteGalsState&, const NativeLeafDomain&, const InterfacePhaseGeometry::MappedFacePartitions&,
        const LevelSetInterfacePanels::NativeComponents&, const HermiteGalsPhaseVolumeConfig&);
};
// Read-only measurement on ONE retained Q4 geometry at its declared time.
// All inventory cells must have scalar support. The source's M is not an
// input to the quadrature and is never overwritten, even on initial templates.
HermiteGalsPhaseVolumeResult integrateHermiteGalsPhaseVolumes(const QuadtreeMesh& mesh,
    const HermiteGalsState& source, const HermiteGalsPhaseVolumeConfig& config = {});
// Publication reuses its final immutable topology and actual half cache.
// Wrong source/layout/transport rejects before any result can escape.
HermiteGalsPhaseVolumeResult integrateHermiteGalsNativeTracePhaseVolumes(const QuadtreeMesh& mesh,
    const HermiteGalsState& source, const NativeLeafDomain& domain,
    const InterfacePhaseGeometry::MappedFacePartitions& partitions,
    const LevelSetInterfacePanels::NativeComponents& topology, const HermiteGalsPhaseVolumeConfig& config = {});
}
