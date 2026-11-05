#pragma once

#include <cstddef>
#include <memory>
#include <vector>
#include "Core/Geometry/LevelSetMappedSegment.h"

class QuadtreeMesh;
class NativeLeafDomain;
class LevelSetInterfacePanels;
class LevelSetMappedSegment;
namespace LevelSetHermiteGeometry { class Snapshot; }
namespace LAMP::Interface { class HermiteGalsPublication; }

// Phase geometry on the existing selected native leaves. No momentum boundary
// law, velocity reconstruction, time integration or Hydro solver is required.
namespace InterfacePhaseGeometry
{
// Original exact scalar/domain layout check, shared by geometric measures and
// the remaining Hydro flux adapter. Does not change either object.
void validateLayout(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    const LevelSetHermiteGeometry::Snapshot& snapshot);

// One retained partition of every actual active vertex-midpoint half-face,
// shared by fixed-boundary moments, phase fluxes and later free-arc cutting.
// Explicitly uses the shared-bound map, NOT the legacy nominal-h partition.
// Supports uniform and native-leaf snapshots, retaining the actual fine
// subfaces at 2:1 transitions. Native queries require these retained cuts;
// the legacy non-cached partition overloads explicitly reject that layout.
// The source and selected native domain must remain immutable through every alias.
class MappedFacePartitions
{
public:
    MappedFacePartitions(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
        std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source);
    // Explicit same-source native Q4 cache, including positive-length zero
    // support and its two normal Taylor phase traces. Default construction
    // retains the legacy zero rejection. No phase rate/volume is produced.
    static MappedFacePartitions prepareNativeTrace(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
        std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source);
    ~MappedFacePartitions();
    MappedFacePartitions(MappedFacePartitions&&) noexcept;
    MappedFacePartitions& operator=(MappedFacePartitions&&) noexcept;
    MappedFacePartitions(const MappedFacePartitions&) = delete;
    MappedFacePartitions& operator=(const MappedFacePartitions&) = delete;
    void validate(const QuadtreeMesh& mesh, const NativeLeafDomain& domain) const;
    const LevelSetHermiteGeometry::Snapshot& source() const;
    const LevelSetMappedSegment& halfFace(int nativeRow, int half) const;
    // Available only on prepareNativeTrace; ordinary halfFace uses this
    // query's legacy view and rejects any positive-length coincident support.
    const LevelSetMappedSegment::NativeQuarticTrace& traceHalfFace(int nativeRow,int half) const;
    // Typed mode counts value events only here, never normal Taylor events.
    std::size_t rootIncidenceCount() const;
    std::size_t traceEventIncidenceCount() const; // All typed value/normal events.
    // Refine roots()[rootIndex] on this actual half-face, without changing
    // the scalar, segment or root-incidence count. Replacement is committed
    // only after successful refinement and validation; failure leaves the
    // cache unchanged. A successful call invalidates references previously
    // returned by halfFace(nativeRow, half). All derived integrals must be
    // discarded and recomputed against the updated common partitions.
    void refineRoot(int nativeRow, int half, std::size_t rootIndex);
    // Same immutable event refinement; success invalidates both typed and
    // ordinary references on this half. Discard/recompute derived measures.
    // A failed refinement preserves the entire cache and its incidence counts.
    void refineTraceEvent(int nativeRow,int half,std::size_t eventIndex);
private:
    // Only the immutable publication owner may reuse the full layout proof
    // made by this cache's constructor. It still validates the current mesh
    // through its state before checking these actual retained bindings.
    friend class LAMP::Interface::HermiteGalsPublication;
    void validatePublicationBinding(const NativeLeafDomain&,
        const LevelSetHermiteGeometry::Snapshot&) const;
    struct Impl;
    MappedFacePartitions(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
        std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source,bool nativeTrace);
    std::unique_ptr<Impl> impl_;
    const Impl& checked() const;
};

// Negative-phi phase volume in a fixed fitted RZ domain. Panels must describe
// the supplied shared scalar's complete negative-to-positive free interface.
// Uses F=z e_z: free-surface GL8 plus the actual wetted horizontal fixed/domain
// boundary pieces. Not a cell-center phase-count volume or a cap-only integral.
// Returned units are code length cubed; positive finite volume is required.
long double volume(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    const LevelSetHermiteGeometry::Snapshot& snapshot, const LevelSetInterfacePanels& panels);

struct BoundaryMoments
{
    // Fixed native CONTROL-VOLUME boundary contributions, NOT complete phase
    // volumes. The free-interface contribution still needs its actual CV
    // overlay. All moments use F=z e_z and a common axial reference of zero,
    // so shared face moments cancel with exactly opposite signs.
    // Face orientation is outward from its negative-coordinate owner.
    std::vector<long double> negativePhaseByNativeFace, negativePhaseByCell;
    // Floating boundary-integral volume of the actual closed CV boundary;
    // not an interval enclosure or a certified positive lower bound.
    std::vector<long double> capacityByCell;
    long double exteriorNegativeMoment = 0.0L, summedCellNegativeMoment = 0.0L;
    long double assemblyClosure = 0.0L;
};

// Uses each actual vertex->midpoint->vertex half-face, retaining its small
// slopes. No horizontal-face midpoint shortcut, rectangular capacity, or
// canonical-panel-to-native-cell ownership assumption is used. Solid-only
// faces/cells are NaN. A phase boundary moment is not transported volume.
BoundaryMoments boundaryMoments(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    const LevelSetHermiteGeometry::Snapshot& snapshot);
BoundaryMoments boundaryMoments(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    const MappedFacePartitions& partitions);
}
