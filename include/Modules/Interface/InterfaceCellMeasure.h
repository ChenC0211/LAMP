#pragma once

#include <array>
#include <cstddef>
#include <memory>
#include <type_traits>
#include <vector>
#include "Core/Geometry/ActiveHermiteSensitivity.h"
#include "Core/Geometry/LevelSetMappedInterface.h"
#include "Core/Geometry/LevelSetInterfacePanels.h"
#include "Modules/Interface/InterfacePhaseGeometry.h"

class QuadtreeMesh;
class NativeLeafDomain;
class InterfaceState;
class LevelSetInterfacePanels;
namespace LevelSetHermiteGeometry { class Snapshot; }
namespace LAMP::Interface { class HermiteGalsPhaseVolumeResult; namespace Detail { class PublicationTracePermit; } }

namespace InterfaceCellMeasure
{
template<std::size_t Count> struct CellArc
{
    static_assert(Count==16 || Count==25,"unsupported interface geometry sensitivity layout");
    int panelIndex = -1, cellIndex = -1;
    double coordinateBegin = 0.0, coordinateEnd = 0.0;
    // Explicit shared-bound path: use independent offset/length; the legacy
    // begin/end above are display approximations only and may coincide.
    bool mappedCoordinate = false, fromEnd = false;
    long double coordinateOffset = 0.0L, coordinateLength = 0.0L;
    std::array<int, 4> vertexIndices{{-1,-1,-1,-1}}; // Uniform source, SW/SE/NE/NW.
    std::array<long double,Count> sensitivity{};
    long double freeBoundaryMoment = 0.0L;
    // Native leaf source: formal corner arrays above are not solve columns.
    // The retained scalar has applied its hanging-node chain rule here.
    int nativeSourceCell = -1;
    using IndependentSensitivity=std::conditional_t<Count==16,
        LAMP::Geometry::ActiveHermiteNodalSensitivity,LAMP::Geometry::ActiveQuarticHermiteSensitivity>;
    std::vector<IndependentSensitivity> independentSensitivity;
};
using Arc=CellArc<16>;
using QuarticArc=CellArc<25>; // Native only: no logical uniform/corner-only correction.
template<class ArcType> struct VolumeMeasure
{
    InterfacePhaseGeometry::BoundaryMoments fixedBoundary;
    std::vector<ArcType> arcs;
    // Full native cell indexing; inactive values remain NaN. Raw floating
    // boundary integrals, without clipping negative/overfilled roundoff.
    std::vector<long double> freeBoundaryMomentByCell, negativeVolumeByCell;
    long double summedNegativeVolume = 0.0L;
    int crossedPanels = 0;
    std::size_t rootRefinements = 0;
};
// Named structs retain lightweight forward declarations in IO consumers.
struct Measure : VolumeMeasure<Arc> {};
struct QuarticMeasure : VolumeMeasure<QuarticArc> {};

// Whole native CVs with a complete topology: no correction Jacobian is
// requested or returned. A distinct type prevents absent derivatives from
// being interpreted as zero sensitivities by correction consumers.
struct Volumes
{
    InterfacePhaseGeometry::BoundaryMoments fixedBoundary;
    // Global-origin free moments retain the shared-face ledger. Q4 discovery
    // uses the region consumer's local Green/smaller-phase volume in cut cells;
    // its rounded value need not equal global fixed + free. Legacy evaluation
    // retains the original global sum. No clipping or inventory correction.
    std::vector<long double> freeBoundaryMomentByCell, negativeVolumeByCell;
    long double summedNegativeVolume = 0.0L;
    long double momentErrorEstimate = 0.0L; // Sum of GL8/subdivided estimates, not a bound.
    std::vector<long double> momentErrorEstimateByCell; // Same estimates, native indexing; inactive NaN.
};


// Complete discovery and native phase volume integration, without a
// legacy InterfaceState. Owns no second scalar. Certificates cannot be supplied
// by the caller: they come from discovery on this exact immutable source.
// Native Q4 cut volumes/estimates share PhaseRegion's local integration path;
// global fixed/free ledger moments and pure certificates remain unchanged.
struct DiscoveredVolumes
{
    Volumes measure;
    std::vector<int> cellPhase;
    std::size_t components = 0, graphs = 0, faceRootIncidences = 0;
};
DiscoveredVolumes discoverVolumes(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source,
    const std::vector<Vector2D>& boundaryContacts = {},
    const std::vector<LevelSetInterfacePanels::SymmetryPlane>& symmetryPlanes = {});

// Explicit geometry-only trace port. Zero supports have no fixed half phase:
// their free moments belong only to the actual available liquid incident CV.
// The complete topology certifies its interior phase, independent of the
// canonical graph owner. Both Green ledgers retain the original global origin.
// Caller durations and summed independent cell-worker durations. Worker
// fields can overlap in wall time and must not be added to caller durations.
struct NativeTraceTimings
{
    double topologyValidation = 0, fixedBoundary = 0, ordinarySetup = 0;
    double ordinaryPrepare = 0, ordinaryMerge = 0, graphMomentsWorker = 0;
    double localRegionsWorker = 0, pureCellsWorker = 0, cutSetupWorker = 0;
    double sideMoments = 0, finalLedger = 0, total = 0;
    std::size_t activeCells = 0, cutCells = 0, pureCells = 0;
    std::size_t ordinaryGraphs = 0, faceHalves = 0;
    // Actual joint global/local Core work, merged in original native-cell order.
    LevelSetMappedInterface::NativeMomentWorkProfile coreWork;
};
struct NativeTraceLedger
{
    NativeTraceTimings timings;
    std::vector<long double> fixedRadialByCell, freeRadialByCell, axialVolumeByCell, radialVolumeByCell, radialCapacityByCell;
    long double exteriorFixedRadial = 0, fixedRadialClosure = 0;
    long double maxRelativeGreenDefect = 0, maxRelativeCapacityDefect = 0;
    std::size_t sideGraphs = 0, sideLiquidOwners = 0, zeroHalfFaces = 0;
};
struct DiscoveredNativeTraceVolumes : DiscoveredVolumes { NativeTraceLedger ledger; };
DiscoveredNativeTraceVolumes measureNativeTraceVolumes(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    const InterfacePhaseGeometry::MappedFacePartitions& partitions,
    const LevelSetInterfacePanels::NativeComponents& topology,
    const std::vector<Vector2D>& boundaryContacts = {},
    const std::vector<LevelSetInterfacePanels::SymmetryPlane>& symmetryPlanes = {});
namespace Detail
{
class PublicationTraceVolumeAccess
{
    friend class ::LAMP::Interface::HermiteGalsPhaseVolumeResult;
    static DiscoveredNativeTraceVolumes measure(const QuadtreeMesh&,const NativeLeafDomain&,
        const InterfacePhaseGeometry::MappedFacePartitions&,const LevelSetInterfacePanels::NativeComponents&,
        const std::vector<Vector2D>&,const std::vector<LevelSetInterfacePanels::SymmetryPlane>&,
        const ::LAMP::Interface::Detail::PublicationTracePermit&);
};
}
DiscoveredNativeTraceVolumes discoverNativeTraceVolumes(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source,
    const std::vector<Vector2D>& boundaryContacts = {},
    const std::vector<LevelSetInterfacePanels::SymmetryPlane>& symmetryPlanes = {});

// Same complete discovery and local moments, plus the existing Q4
// independent-center sensitivity integral. No legacy state or second scalar.
// The derivative uses the original mapped GL8 sensitivity, while volumes use
// the shared local region integral; safeguarded solvers check actual residuals.
struct DiscoveredBubbleVolumes : DiscoveredVolumes
{
    std::vector<long double> centerDerivative; // native indexing; inactive NaN
};
DiscoveredBubbleVolumes discoverBubbleVolumes(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source,
    const std::vector<Vector2D>& boundaryContacts = {},
    const std::vector<LevelSetInterfacePanels::SymmetryPlane>& symmetryPlanes = {});

// Partition retained scalar graphs at actual vertex-midpoint half-face roots.
// Each resolved subarc is assigned to the actual four-subquad CV using its
// cell-local position, not a rounded global Frame.point or panel.cellIndex.
// All GL8 nodes must agree on one active owner. Root inventories and source
// identity are checked; unresolved/coincident intervals reject, not merge.
//
// This is floating root/GL8 geometry, not an exact watertight overlay or an
// interval enclosure. The shared scalar's nominal h embedding is preserved;
// no (physical upper-lower) metric replacement is introduced. Independent
// volume/flux error checks remain necessary before accepting a time update.
Measure evaluate(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    const LevelSetHermiteGeometry::Snapshot& snapshot, const LevelSetInterfacePanels& panels);

// Shared-bound counterpart: fixed faces and free-curve cuts reuse exactly
// these retained partitions. All curve nodes, physical metric, gradient and
// sensitivity use that same embedding. Panel endpoint identities select the
// original regular graph, not a rounded global root. No old point-map query
// is used to restore a tiny interval. Same floating-volume caveats as above;
// this does not by itself correct or accept a transported interface.
Measure evaluate(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    const InterfacePhaseGeometry::MappedFacePartitions& partitions,
    const LevelSetInterfacePanels& panels);

// Native rectangular leaves: each scalar cell is the actual straight-sided
// CV, including split coarse/fine sides. Use ALL retained actual half faces,
// the same mapped curve integral and its independent-node pullback. The list
// must provide one regular graph for every non-pure active leaf; duplicates,
// hidden/omitted zero sets and mixed addressing reject. This entry does not
// discover/order the whole interface or supply its force/projection panels.
Measure evaluate(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    const InterfacePhaseGeometry::MappedFacePartitions& partitions,
    const std::vector<LevelSetMappedInterface::Graph>& nativeGraphs);

// Same actual CV/face-root/GL8 assembly, with the complete native quartic
// independent-data transpose. Distinct type prevents silent 25-to-16 loss.
QuarticMeasure evaluateQuartic(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    const InterfacePhaseGeometry::MappedFacePartitions& partitions,
    const std::vector<LevelSetMappedInterface::Graph>& nativeGraphs);

// Explicit preparation step: refine only face roots whose enclosures prevent
// curve ownership classification, then recompute the COMPLETE measure. The
// updated cache must also be used for subsequent fluxes. Previously computed
// moments/fluxes are invalid after any refinement. The work limit is not an
// accuracy gate: exhaustion rejects instead of accepting an unresolved arc.
// A later failure can leave earlier successful root refinements in the cache;
// no partial Measure is returned and no prior derived integral remains valid.
Measure evaluateRefined(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    InterfacePhaseGeometry::MappedFacePartitions& partitions,
    const LevelSetInterfacePanels& panels, std::size_t maxRootRefinements = 256);
Measure evaluateRefined(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    InterfacePhaseGeometry::MappedFacePartitions& partitions,
    const std::vector<LevelSetMappedInterface::Graph>& nativeGraphs, std::size_t maxRootRefinements = 256);
// Automatically connected native interface; additionally verifies the owned
// scalar identity and RZ measure contract before consuming its original graphs.
Measure evaluateRefined(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    InterfacePhaseGeometry::MappedFacePartitions& partitions,
    const LevelSetInterfacePanels::NativeGraphChain& chain, std::size_t maxRootRefinements = 256);
QuarticMeasure evaluateQuarticRefined(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    InterfacePhaseGeometry::MappedFacePartitions& partitions,
    const LevelSetInterfacePanels::NativeGraphChain& chain, std::size_t maxRootRefinements = 256);

// Complete module-owned topology in native XY or RZ. All components/arcs
// participate, including several graphs per cell. Graph-free cells consume
// the same discovery certificate; callers cannot supply an unchecked phase
// mask. State/source/domain identities must agree. No volume clipping occurs.


}
