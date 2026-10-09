#include "Modules/Hydrodynamics/IO/DictionaryFittedFlow.h"
#include "Core/IO/ConfigurationTrace.h"

#include <set>
#include <string>

FittedFlowCaseConfig readFittedFlowCase(const DictionaryInput::Node& root, DictionaryInput::ConfigurationTrace& trace)
{
    using Node = DictionaryInput::Node;
    root.allow({"Model", "Initial", "BoundaryConditions", "Numerics"});
    FittedFlowCaseConfig result;
    const auto word = [&](const Node& node, const char* value, const char* consumer) {
        return trace.record(node, node.choice({value}), "model", consumer);
    };
    const auto id = [&](const Node& node, const char* consumer) {
        const int value = node.integer();
        if (value <= 0) node.fail("a positive 1-based ID");
        return trace.record(node, value, "1-based ID", consumer);
    };
    const auto& model = root.block("Model");
    model.allow({"formulation", "coordinateSystem", "liquidPhase", "gasPhase"});
    word(model.get("formulation"), "incompressibleTwoPhase", "fitted variable-property NS equations");
    word(model.get("coordinateSystem"), "axisymmetricRZ", "fitted RZ domain and physical measures");
    result.liquidPhase = id(model.get("liquidPhase"), "Hydro liquid material binding");
    result.gasPhase = id(model.get("gasPhase"), "Hydro gas/backflow material binding");
    if (result.liquidPhase == result.gasPhase) model.get("gasPhase").fail("distinct liquid and gas phase IDs");
    const auto& initial = root.block("Initial");
    initial.allow({"model"});
    word(initial.get("model"), "rest", "zero initial velocity and pressure; boundary pressure applied by NS");

    const auto& boundaries = root.block("BoundaryConditions");
    boundaries.allow({"Supply", "Outlet", "Ambient", "solidWallModel", "incomingMomentumModel"});
    result.solidWallModel = word(boundaries.get("solidWallModel"), "noSlip", "fitted fixed solid-wall velocity");
    result.incomingMomentumModel = word(boundaries.get("incomingMomentumModel"), "interiorContinuationZeroShear",
        "pressure-reservoir incoming momentum closure");
    std::set<int> boundaryIds;
    const auto pressure = [&](const Node& block, bool incomingField) {
        if (incomingField) block.allow({"boundaryId", "pressure", "incomingPhase"});
        else block.allow({"boundaryId", "pressure"});
        FittedPressureBoundarySI boundary;
        boundary.boundaryId = id(block.get("boundaryId"), "fitted pressure port geometry binding");
        if (!boundaryIds.insert(boundary.boundaryId).second) block.get("boundaryId").fail("distinct pressure boundary IDs");
        const auto& p = block.get("pressure");
        boundary.pressureSI = trace.record(p, p.quantity(UnitParser::Dimension::Pressure, "Pa"), "Pa",
            "fitted pressure gradient and incoming reservoir; does not prescribe flow rate");
        return boundary;
    };
    result.supply = pressure(boundaries.block("Supply"), false);
    const auto& outlet = boundaries.block("Outlet");
    result.outlet = pressure(outlet, true);
    const auto& incoming = outlet.get("incomingPhase");
    if (incoming.word() == "none")
    {
        trace.record(incoming, std::string("none"), "model", "pressure outlet rejects unconfigured backflow");
        result.outletIncomingPhaseId = -1;
    }
    else
    {
        result.outletIncomingPhaseId = id(incoming, "outlet reservoir incoming material");
        if (result.outletIncomingPhaseId != result.liquidPhase && result.outletIncomingPhaseId != result.gasPhase)
            incoming.fail("the liquid or gas phase ID, or none");
    }
    for (const auto& ambient : boundaries.block("Ambient").entries)
        result.ambient.push_back(pressure(ambient, false));

    const auto& numerics = root.block("Numerics");
    numerics.allow({"InterfaceForce", "Conservation"});
    const auto& force = numerics.block("InterfaceForce");
    force.allow({"halfWidthCells"});
    const auto& halfWidth = force.get("halfWidthCells");
    result.capillaryHalfWidthCells = trace.record(halfWidth, halfWidth.number(), "cells", "compact force-band half-width/max(dr,dz)");
    if (result.capillaryHalfWidthCells <= 0.0) halfWidth.fail("a positive force-band half-width in cells");
    const auto& conservation = numerics.block("Conservation");
    conservation.allow({"momentumBalance", "localContinuity", "outletAdmissibility", "phaseFluxBalance"});
    const auto tolerance = [&](const char* name, const char* consumer) {
        const auto& entry = conservation.get(name);
        const double value = entry.number();
        if (value <= 0.0 || value >= 1.0) entry.fail("a relative acceptance tolerance in (0,1)");
        return trace.record(entry, value, "dimensionless", consumer);
    };
    result.momentumBalance = tolerance("momentumBalance", "full fitted CN momentum balance");
    result.localContinuity = tolerance("localContinuity", "fitted per-cell fluid continuity");
    result.outletAdmissibility = tolerance("outletAdmissibility", "pressure outlet admissibility plus unchanged roundoff term");
    result.phaseFluxBalance = tolerance("phaseFluxBalance", "moving phase ledger and pressure-port consistency");
    return result;
}
