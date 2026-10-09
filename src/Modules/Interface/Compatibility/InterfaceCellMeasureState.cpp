#include "Modules/Interface/Compatibility/InterfaceTopology.h"
#include "Modules/Interface/Detail/InterfaceNativeMeasure.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include <stdexcept>
namespace InterfaceCellMeasure {
namespace {
using CompleteNativeInventory=Detail::NativeInventory;
void require(bool ok,const char* reason){if(!ok)throw std::runtime_error(std::string("[InterfaceCellMeasure] ")+reason);}
CompleteNativeInventory completeInventory(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    const InterfacePhaseGeometry::MappedFacePartitions& partitions, const InterfaceState& state)
{
    state.validateMesh(mesh);
    require(state.geometry().get()==&partitions.source() && state.hasTopology(),
        "complete volume inventory requires the same module-owned scalar and topology");
    CompleteNativeInventory result;
    result.phase=&state.topology().cellPhase();
    for (const auto& surface:state.topology().surfaces()) {
        const auto& chain=surface.geometry.nativeChain();
        require(chain.usesSource(partitions.source()) && chain.axisymmetric()==domain.axisymmetric(),
            "complete volume inventory changed source or XY/RZ measure");
        result.graphs.insert(result.graphs.end(),chain.graphs().begin(),chain.graphs().end());
    }
    return result;
}
}
Measure evaluateRefined(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    InterfacePhaseGeometry::MappedFacePartitions& partitions, const InterfaceState& state,
    std::size_t maxRootRefinements)
{ return Detail::evaluateRefined(mesh,domain,partitions,completeInventory(mesh,domain,partitions,state),maxRootRefinements); }

QuarticMeasure evaluateQuarticRefined(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
    InterfacePhaseGeometry::MappedFacePartitions& partitions, const InterfaceState& state,
    std::size_t maxRootRefinements)
{ return Detail::evaluateQuarticRefined(mesh,domain,partitions,completeInventory(mesh,domain,partitions,state),maxRootRefinements); }

Volumes evaluateVolumes(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const InterfacePhaseGeometry::MappedFacePartitions& partitions,const InterfaceState& state)
{ return Detail::evaluateVolumes(mesh,domain,partitions,completeInventory(mesh,domain,partitions,state)); }
}
