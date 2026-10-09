#pragma once

#include <cstddef>
#include <memory>
#include <vector>
#include "Core/Geometry/ActiveHermiteSensitivity.h"
#include "Modules/Interface/InterfacePhaseGeometry.h"
#include "Modules/Interface/Compatibility/InterfaceTransportInput.h"

class QuadtreeMesh;

class LevelSetInterfacePanels;

class LevelSetMappedSegment;
namespace LevelSetHermiteGeometry { class Snapshot; }
template<typename T> class VEdgeField;
template<typename T> class HEdgeField;

namespace InterfacePhaseFlux
{
struct FluxBudget
{
    // Native indices, not labels or logical positions. Solid-only faces and
    // inactive cells contain NaN and must not be interpreted as zero fluid.
    // Face flux is oriented +r/+z; a shared flux is added with opposite signs
    // to its two active owners. Units are code velocity times domain face measure.
    std::vector<long double> positiveFluxByNativeFace, outwardFluxByCell;
    long double boundaryOutward = 0.0L, summedCellOutward = 0.0L;
    long double absoluteFaceFlux = 0.0L, assemblyClosure = 0.0L;
    int activeFaces = 0, partitionedHalfFaces = 0;
};

// Negative-phase flux ledger on the SAME actual vertex-midpoint-vertex faces
// as the flux-compatible velocity. Phase roots retain nested local parameters;
// no global root subtraction or centre-phase replacement is used. No-slip and
// actual r=0 axis fluxes are exactly zero independently of phase. This is NOT a
// per-cell geometric volume: canonical Hermite panels have not been recut by
// the slightly non-Cartesian native control volumes. It supplies conservative
// shared face budgets, not a certificate for the current pointwise transport.
FluxBudget phaseFluxBudget(const QuadtreeMesh& mesh, const InterfaceTransportDomain& domain,
    const LevelSetHermiteGeometry::Snapshot& snapshot, const InterfaceNormalFluxView& velocity);
FluxBudget phaseFluxBudget(const QuadtreeMesh& mesh, const InterfaceTransportDomain& domain,
    const InterfacePhaseGeometry::MappedFacePartitions& partitions, const InterfaceNormalFluxView& velocity);

template<class IndependentSensitivity> struct PhaseFluxSensitivity
{
    // Native indexing; columns retain their scalar's independent-data type. Empty
    // entries mean no variation, including inactive cells; not a fluid value.
    // Each positive-face contribution is assembled ONCE and added with exact
    // opposite signs to the two owners. This is dQ, never geometric dV.
    std::vector<std::vector<IndependentSensitivity>>
        positiveFluxByNativeFace, outwardFluxByCell;
    int roots = 0, endpointRoots = 0;
};
using FluxSensitivity=PhaseFluxSensitivity<Quantum::Geometry::ActiveHermiteNodalSensitivity>;
using QuarticFluxSensitivity=PhaseFluxSensitivity<Quantum::Geometry::ActiveQuarticHermiteSensitivity>;
// Frozen velocity, current immutable native scalar and SAME refined cuts as
// volume/flux. Simple transverse roots only; exact half endpoints use the
// symmetric generalized derivative of the clipped phase interval (half
// weight). No-slip/axis contributions vanish identically. Does not refit jets.
// HoldFixed is valid only for updates that preserve declared corner-sector jets.
FluxSensitivity phaseFluxSensitivity(const QuadtreeMesh& mesh, const InterfaceTransportDomain& domain,
    const InterfacePhaseGeometry::MappedFacePartitions& partitions, const InterfaceNormalFluxView& velocity,
    Quantum::Geometry::CornerSectorVariation sectors =
        Quantum::Geometry::CornerSectorVariation::RequireAbsent);

// Same frozen-velocity/shared-face assembly; carries independent vertex,
// maximal-edge-midpoint and center columns instead of cubic nodal jets.
QuarticFluxSensitivity quarticPhaseFluxSensitivity(const QuadtreeMesh& mesh, const InterfaceTransportDomain& domain,
    const InterfacePhaseGeometry::MappedFacePartitions& partitions, const InterfaceNormalFluxView& velocity,
    Quantum::Geometry::CornerSectorVariation sectors =
        Quantum::Geometry::CornerSectorVariation::RequireAbsent);


}
