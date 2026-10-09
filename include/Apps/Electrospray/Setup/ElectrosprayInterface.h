#pragma once

struct FittedFlowCaseConfig;
struct FittedInterfaceCaseConfig;

#include <functional>
#include <vector>
#include "Modules/Interface/Compatibility/InterfaceInitialization.h"

class QuadtreeMesh;
class DictionaryGeometry;
class HydroFittedDomain;
template<class Mesh> class ImmersedGeometry;
struct ElectrosprayElectricState;

// The App owns case binding only. Its surface
// is the module state, not an inherited/duplicated App scalar. The initial
// composite solid/phase level set remains separate and immutable.
struct ElectrosprayInterfaceState
{
    InterfaceState surface;
    std::vector<int> freeLineIds, fixedLineIds; // 0-based dictionary line IDs.
    int sourceLevelSetIndex = -1;
    explicit ElectrosprayInterfaceState(const QuadtreeMesh& mesh);
};

ElectrosprayInterfaceState initializeElectrosprayInterface(const QuadtreeMesh& mesh,
    ImmersedGeometry<QuadtreeMesh>& geometry, const DictionaryGeometry& definition,
    const HydroFittedDomain& fluid, const ElectrosprayElectricState& initial,
    const FittedFlowCaseConfig& flowConfig, const FittedInterfaceCaseConfig& interfaceConfig,
    int conductorPhaseId);

void buildElectrosprayActiveInterfaceGeometry(const QuadtreeMesh& mesh, const HydroFittedDomain& fluid,
    ElectrosprayInterfaceState& state, const InterfaceInitialGeometryObserver& observer = {});
