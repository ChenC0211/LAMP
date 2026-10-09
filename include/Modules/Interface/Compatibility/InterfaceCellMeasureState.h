#pragma once
#include "Modules/Interface/InterfaceCellMeasure.h"
class InterfaceState;
namespace InterfaceCellMeasure {
Volumes evaluateVolumes(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    const InterfacePhaseGeometry::MappedFacePartitions& partitions, const InterfaceState& state);
Measure evaluateRefined(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    InterfacePhaseGeometry::MappedFacePartitions& partitions,
    const InterfaceState& state, std::size_t maxRootRefinements = 256);
QuarticMeasure evaluateQuarticRefined(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    InterfacePhaseGeometry::MappedFacePartitions& partitions,
    const InterfaceState& state, std::size_t maxRootRefinements = 256);
}
