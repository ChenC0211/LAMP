#pragma once
#include "Modules/Interface/Compatibility/InterfaceState.h"
#include "Modules/Interface/Compatibility/InterfaceTransportInput.h"
#include "Core/Field/VertField.h"
#include <array>
#include <functional>

struct InterfaceInitialCurveSample { double unsignedDistance; bool negative; };
// Bound geometry in code units. All inputs/callbacks are synchronous borrows,
// never retained. No dictionary, Hydro or electric model enters the module.
struct InterfaceInitialCurveInput
{
    const VertField<double>& unsignedDistance;
    const VertField<double>& compositeSign;
    const std::vector<int>& phaseByCell;
    int negativePhase, positivePhase;
    Vector2D begin, end;
    std::vector<Vector2D> pinnedCoordinates;
    std::array<Vector2D, 2> inwardEndpointTangents;
    std::function<InterfaceInitialCurveSample(Vector2D)> primitive;
};
// Read-only call-scoped observations of original stages, including rejected
// scalar DOFs. Observer exceptions abort and fail the candidate.
struct InterfaceInitialGeometryObserver
{
    std::function<void(int, double, double)> initialTangencyAdjusted;
    std::function<void(const InterfaceScalarPatch&)> patchPrepared;
    std::function<void(const LevelSetHermiteGeometry::Snapshot&)> scalarConstructed;
    std::function<void(const LevelSetInterfacePanels&)> panelsConstructed;
};
namespace InterfaceInitialization
{
// Original signed free distance, exact pins and one-sided jets, including Q4
// data. Produces Prepared, without a usable Hermite scalar/panel pair.
void prepareCurve(const QuadtreeMesh& mesh, const InterfaceTransportDomain& domain,
    const InterfaceInitialCurveInput& input, InterfaceState& freshDestination);
// Original cubic/native-Q4 scalar and panels with initial-only tangent
// compatibility. Publishes Ready, without any physical time advance.
void publishCurve(const QuadtreeMesh& mesh, const InterfaceTransportDomain& domain,
    InterfaceState& prepared, const InterfaceInitialGeometryObserver& observer = {});
}
