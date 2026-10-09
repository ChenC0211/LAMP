#include "Electrospray/Setup/ElectrosprayInterface.h"
#include "Modules/Hydrodynamics/FittedFlowCaseConfig.h"
#include "Modules/Interface/Compatibility/InterfaceCaseConfig.h"
#include "Modules/Hydrodynamics/HydroFittedDomain.h"
#include "Modules/Hydrodynamics/HydroInterfaceTransport.h"
#include "Electrospray/Setup/ElectrosprayElectrostatics.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Geometry/ImmersedGeometry.h"
#include <set>
#include <stdexcept>

namespace {
void require(bool valid, const char* reason) {
    if (!valid) throw std::invalid_argument(std::string("[ElectrosprayInterface] ")+reason);
}
bool same(Vector2D a, Vector2D b) { return a.x == b.x && a.y == b.y; }
}
ElectrosprayInterfaceState::ElectrosprayInterfaceState(const QuadtreeMesh& mesh) : surface(mesh) {}

ElectrosprayInterfaceState initializeElectrosprayInterface(const QuadtreeMesh &mesh,
                                                           ImmersedGeometry<QuadtreeMesh> &geometry,
                                                           const DictionaryGeometry &definition,
                                                           const HydroFittedDomain &fluid,
                                                           const ElectrosprayElectricState &initial,
                                                           const FittedFlowCaseConfig& flowConfig,
                                                           const FittedInterfaceCaseConfig& interfaceConfig,
                                                           int conductorPhaseId)
{
    require(interfaceConfig.contactLineModel == "pinned" && !interfaceConfig.pinnedContactPointIds.empty(),
            "moving interface requires explicit pinned contact-line data");
    fluid.validateMeshIdentity(mesh);
    require(fluid.axisymmetric() && mesh.minCoord.x == 0.0,
            "the current fitted interface binding requires the actual RZ axis");
    require(initial.phaseByCell.size() == static_cast<std::size_t>(mesh.cellIndexNum()),
            "initial material mapping does not match the fluid domain");
    const auto found = definition.boundaryById.find(interfaceConfig.movingBoundaryId);
    require(found != definition.boundaryById.end() && !found->second.ignored &&
                found->second.kind == GeometricBoundary::moving,
            "moving boundary is missing or not moving");
    const auto &moving = found->second;
    const int sourceIndex = moving.levelSetIndex;
    require(sourceIndex >= 0 && sourceIndex < static_cast<int>(geometry.levelSets.size()) &&
                definition.levelSetName[sourceIndex] == interfaceConfig.fluidInterfaceLevelSet,
            "moving boundary must identify the configured source level set");
    require(moving.phase + 1 == flowConfig.liquidPhase && moving.opposedPhase + 1 == flowConfig.gasPhase,
            "moving boundary must be oriented from configured liquid to gas");
    ElectrosprayInterfaceState state(mesh);
    state.sourceLevelSetIndex = sourceIndex;
    std::set<int> selected(moving.lineIds.begin(), moving.lineIds.end()), fixed;
    require(!selected.empty() && selected.size() == moving.lineIds.size(),
            "empty or duplicate moving line selection");
    for (const auto &boundary : definition.boundaries)
        if (!boundary.ignored && boundary.kind == GeometricBoundary::fixed &&
            boundary.levelSetIndex == sourceIndex)
        {
            require((boundary.phase + 1 == conductorPhaseId &&
                     boundary.opposedPhase + 1 == flowConfig.liquidPhase) ||
                        (boundary.opposedPhase + 1 == conductorPhaseId &&
                         boundary.phase + 1 == flowConfig.liquidPhase),
                    "fixed part of composite interface must separate liquid and configured conductor");
            for (int line : boundary.lineIds)
                require(!selected.count(line) && fixed.insert(line).second,
                        "overlapping fixed/moving line ownership");
        }
    std::vector<SignInt> freeReferences;
    std::set<int> covered;
    for (const auto reference : definition.linesOnInterface[sourceIndex])
    {
        require(reference.value >= 0 && reference.value < static_cast<int>(definition.lines.size()),
                "invalid physical line index");
        require(covered.insert(reference.value).second &&
                    (selected.count(reference.value) || fixed.count(reference.value)),
                "source physical lines must have unique explicit fixed/moving ownership");
        if (selected.count(reference.value))
            freeReferences.push_back(reference);
    }
    require(covered.size() == selected.size() + fixed.size() && !freeReferences.empty(),
            "uncovered fixed/moving physical line");
    state.freeLineIds.assign(selected.begin(), selected.end());
    state.fixedLineIds.assign(fixed.begin(), fixed.end());
    const auto endpoint = [&](SignInt reference, bool end) {
        return definition.lines[reference.value].points[(end != reference.isNegative) ? 1 : 0];
    };
    for (std::size_t i = 1; i < freeReferences.size(); ++i)
        require(same(endpoint(freeReferences[i - 1], true), endpoint(freeReferences[i], false)),
                "selected moving lines must form one directed open chain");
    const Vector2D begin = endpoint(freeReferences.front(), false),
                   end = endpoint(freeReferences.back(), true);
    require(!same(begin, end), "this fitted contact binding requires an open free-interface chain");
    VertField<double> unsignedDistance(mesh, 0.0);
    std::vector<std::array<Vector2D, 2>> inwardTangents;
    geometry.initializeUnsignedDistanceToLines(mesh, freeReferences, unsignedDistance, &inwardTangents);
    InterfaceInitialCurveInput input{unsignedDistance, geometry.levelSets[sourceIndex], initial.phaseByCell,
        flowConfig.liquidPhase, flowConfig.gasPhase, begin, end, {},
        {inwardTangents.front()[0], inwardTangents.back()[1]}, {}};
    for (int id : interfaceConfig.pinnedContactPointIds) {
        require(id > 0 && id <= static_cast<int>(definition.points.size()), "pinned Point ID is outside geometry");
        input.pinnedCoordinates.push_back(definition.points[id-1]);
    }
    if (!fluid.uniformLeaves()) {
        const auto primitive = geometry.initialCurveSampler(mesh, freeReferences, sourceIndex);
        input.primitive = [primitive](Vector2D point) {
            const auto value = primitive(point);
            return InterfaceInitialCurveSample{value.unsignedDistance, value.negative};
        };
    }
    InterfaceInitialization::prepareCurve(mesh, HydroInterfaceTransport::domainView(fluid), input, state.surface);
    return state;
}

void buildElectrosprayActiveInterfaceGeometry(const QuadtreeMesh& mesh, const HydroFittedDomain& fluid,
    ElectrosprayInterfaceState& state, const InterfaceInitialGeometryObserver& observer)
{
    InterfaceInitialization::publishCurve(mesh, HydroInterfaceTransport::domainView(fluid), state.surface, observer);
}
