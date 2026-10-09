#pragma once
#include "Modules/Interface/InterfaceCellMeasure.h"

// Internal arithmetic port, not a second geometry/state owner. Callers must
// retain the source-bound discovery certificate for these exact graphs.
namespace InterfaceCellMeasure::Detail {
struct NativeInventory {
    std::vector<LevelSetMappedInterface::Graph> graphs;
    const std::vector<int>* phase = nullptr;
};
Measure evaluateRefined(const QuadtreeMesh&,const NativeLeafDomain&,
    InterfacePhaseGeometry::MappedFacePartitions&,const NativeInventory&,std::size_t);
QuarticMeasure evaluateQuarticRefined(const QuadtreeMesh&,const NativeLeafDomain&,
    InterfacePhaseGeometry::MappedFacePartitions&,const NativeInventory&,std::size_t);
Volumes evaluateVolumes(const QuadtreeMesh&,const NativeLeafDomain&,
    const InterfacePhaseGeometry::MappedFacePartitions&,const NativeInventory&);
}