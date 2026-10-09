#pragma once
#include "Modules/Interface/Compatibility/InterfaceState.h"
#include "Core/Geometry/ActiveQuarticHermiteTransfer.h"
#include "Modules/Interface/Compatibility/InterfaceConservativePhase.h"
#include "Modules/Interface/Compatibility/InterfaceConservativeGeometry.h"
#include <string>

class NativeLeafDomain;
namespace InterfaceMeshTransfer
{
// Add every active cell without a closed-cell sign certificate to the mask.
// The caller unions module support before applying its mesh-neighbor halo.
// Existing marks are retained; no interface or mesh state is changed.
void protectScalarSupport(const QuadtreeMesh& mesh, const InterfaceState& source,
    const NativeLeafDomain& domain, std::vector<std::uint8_t>& protectedCells);

struct Result
{
    std::vector<int> sourceToTargetCell, sourceToTargetVertex;
    ActiveQuarticHermiteTransfer::Report report;
};
// Existing zero-time native Q4 transfer for a pinned-to-axis RZ interface.
// Reuses Core's polynomial/phase certificates and material-panel rebinding;
// protected leaves retain the original scalar acceptance mask. No motion,
// volume correction, charge evolution or coupled AMR acceptance occurs.
// Source is unchanged; destination must be fresh and is discarded on failure.
Result preparePinnedToAxis(const QuadtreeMesh& sourceMesh, const InterfaceState& source,
    const QuadtreeMesh& target, const NativeLeafDomain& targetDomain,
    const std::vector<std::uint8_t>& protectedSourceCells, InterfaceState& destination);

struct ConservativeLimits
{
    // Caller-supplied absolute budgets; no amplitude/distance equivalence is
    // assumed. The second is sum |V_geometry - V_inventory| on the target.
    double scalarDifference = 0;
    long double geometryInventoryL1 = 0;
    // Candidate auxiliary curve fit at AMR only; retains the scalar unless
    // synchronizeGeometry is requested. Closed branches, RZ axis ends and
    // fixed contacts use their respective curve representations.
    // Unsupported/failed fits retain the original geometric allocation and
    // are reported below; invalid boundary declarations reject the candidate.
    bool sharedCurveReference = false;
    // Joint C/scalar candidate: new child phase follows the constrained curve,
    // then the existing inventory-driven reconstruction publishes one scalar.
    // Original pure parents and every parent inventory remain fixed. The
    // scalarDifference budget still applies to polynomial mesh transfer;
    // geometryInventoryL1 applies to the FINAL synchronized scalar versus C.
    // No momentum/charge transfer is supplied by this option.
    bool synchronizeGeometry = false;
    // Explicit permission to reextend a geometrically certified pure corner
    // sector during refinement. Its independent C must also be exactly pure.
    // Raw scalar defects stay in the report; the scalarDifference gate applies
    // unchanged to every other target leaf. Default preserves the strict gate.
    bool reextendPureSectors = false;
};
struct ConservativeReport
{
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> sourceGeometry;
    std::vector<double> scalarDifferenceBound;
    std::vector<long double> geometricVolumeByCell;
    double maxScalarDifferenceBound = 0;
    // Raw scalar defects above include certified pure wet-sector reextension.
    // Its target cell IDs are separate; independent source C must be exactly
    // the same pure phase. All other cells still obey scalarDifference.
    std::vector<int> pureSectorExtensionCells;
    double maxNonExtensionDifferenceBound = 0;
    bool geometrySynchronized = false;
    double maxSynchronizationScalarChange = 0; // Diagnostic, not a distance bound.
    double maxSynchronizationInterfaceShift = 0; // Sampled native projection, not a continuous bound.
    long double inventoryVolume = 0, geometricVolume = 0, geometryInventoryL1 = 0;
    long double maxFractionDifference = 0, unassignedInventoryVolume = 0;
    struct Component
    {
        std::uint64_t liquidId = 0;
        long double inventoryVolume = 0, geometricVolume = 0, geometryInventoryL1 = 0;
    };
    std::vector<Component> components;
    std::size_t retainedLiquids = 0, retainedSurfaces = 0;
    // Terminal model leaves may be renumbered as the surrounding support is
    // refined. Their physical boxes and independent C must remain identical;
    // refining/coarsening an already modeled leaf is not supported.
    std::size_t retainedEventCells = 0;
    struct SharedCurve {
        std::uint64_t surfaceId=0;
        bool applied=false;
        int controls=0,rank=0,iterations=0;
        int initialSupportRepairs=0;
        double maxInitialDisplacement=0;
        double maxNormalDisplacement=0;
        long double parentFractionResidual=0,maxIntegrationDifference=0;
        std::size_t changedPureChildren=0;
        std::string reason;
    };
    std::vector<SharedCurve> sharedCurves;
};
// Zero-time cubic geometry + independent C remap, with unchanged domain and
// one-to-one liquid/surface lineage. Retains IDs, event history, ancestry and
// physical clock; old face fluxes and synchronization reports are not reused.
// Pure-liquid overlap witnesses and native curve-cell overlap must resolve
// every component. Ambiguous/unresolved lineage rejects, never invents events.
// Supports closed XY/RZ interfaces, axis ends, retained outer reflection
// sections (with explicit target symmetry laws), and fixed wall
// contacts. Contacts require targetBoundaries bound to the actual target
// domain, with stationary wall laws at the retained endpoints. Pure wet-corner
// refinement may explicitly enable reextendPureSectors; its reentrant vertex,
// local same-level topology and exact pure C must survive. Legacy pinned-state
// metadata remains unsupported. Normals/curves are
// source-inventory donors, as in InterfaceConservativePhase::transferMesh.
// For refined mixed parents, retained native geometry first supplies child
// volumes. A bounded local projection preserves each parent's independent C;
// reference-pure children stay fixed. With synchronizeGeometry the new child
// reference comes from the constrained curve, and the final scalar is rebuilt
// from its C at zero physical time. Original pure parents stay fixed. No fit
// applied (e.g. only moving a pure support band): keep the transferred scalar.
// Infeasible geometric partitions use
// the supplied C-matched donors. Counts/correction are in meshTransferReport().
// The returned report measures geometry/C independently. The destination has
// no physical-step material report; optional zero-time synchronization has its
// own geometry/curvature diagnostics. Mixed-leaf transport and
// height stencils still enforce their existing capability checks. This does
// not transfer momentum/charge/material panels or accept an App AMR step.
// Fresh destination required; any preparation failure discards it, source unchanged.
ConservativeReport prepareConservative(const QuadtreeMesh& sourceMesh,
    const NativeLeafDomain& sourceDomain, const InterfaceState& source,
    const QuadtreeMesh& targetMesh, const NativeLeafDomain& targetDomain,
    const std::vector<Vector2D>& sourceNormals,
    const std::vector<std::optional<Quantum::Geometry::QuadraticCellCut::Shape>>& sourceCurves,
    const ConservativeLimits& limits, InterfaceState& destination,
    const InterfaceConservativeGeometry::BoundaryInput* targetBoundaries=nullptr);
}
