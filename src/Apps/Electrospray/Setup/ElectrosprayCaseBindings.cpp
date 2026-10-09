#include "Electrospray/Setup/ElectrosprayCaseBindings.h"
#include "Modules/Hydrodynamics/FittedFlowCaseConfig.h"
#include "Modules/Interface/Compatibility/InterfaceCaseConfig.h"
#include "Modules/Electricity/FittedElectricCaseConfig.h"
#include "Core/IO/DictionaryGeometry.h"
#include "Core/IO/DictionaryMesh.h"
#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>

void validateElectrosprayCaseBindings(const DictionaryGeometry& g, const DictionaryMesh& mesh,
    const FittedFlowCaseConfig& flow, const FittedInterfaceCaseConfig& interfaceConfig,
    const FittedElectricCaseConfig& electric)
{
    const auto require = [](bool valid, const std::string& message) {
        if (!valid) throw std::invalid_argument("[case binding] " + message);
    };
    require(g.hasDomain && g.domainMin.x == 0 && g.domainMax.x > 0 && g.domainMax.y > g.domainMin.y,
        "geometry.dict/Domain: expected a finite RZ domain starting at r=0");
    require(std::isfinite(g.domainMax.x) && std::isfinite(g.domainMin.y) && std::isfinite(g.domainMax.y), "non-finite Domain");
    require(g.phaseName.size() == 3, "geometry.dict/Phase: current fitted workflow requires conductor, liquid and gas only");
    for (int id : {flow.liquidPhase, flow.gasPhase, electric.conductorPhaseId})
        require(id > 0 && id <= static_cast<int>(g.phaseName.size()), "phase ID " + std::to_string(id) + " is outside geometry.dict/Phase");
    require(g.periodicPairs.empty() && !g.extensionEnabled("MovingEB") && !g.extensionEnabled("MotionSchedule"),
        "geometry.dict: periodic/prescribed solid motion is unsupported by this free-surface workflow");
    if (mesh.physicsRefinement.enabled)
    {
        const auto& flow = mesh.physicsRefinement;
        require(mesh.runtimeAdaptation.enabled && !flow.indicators.empty(),
            "mesh.dict/PhysicsRefinement requires runtime AMR and named Hydro indicators");
        std::set<std::string> caps;
        for (const auto& name : flow.indicators)
        {
            require(name == "Velocity" || name == "VelocityGradient", "unsupported Hydro refinement indicator: " + name);
            const bool speed = name == "Velocity";
            const auto key = speed ? "velocity" : "velocitygradient";
            const double reference = speed ? flow.velocityThresholdSI : flow.velocityGradientThresholdSI;
            require(std::isfinite(reference) && reference > 0 && flow.maxLevelCaps.count(key)
                && flow.maxLevelCaps.at(key) >= 0, "each flow indicator requires a positive SI threshold and explicit cap");
            caps.insert(key);
        }
        require(caps.size() == flow.maxLevelCaps.size(), "unused flow refinement cap");
        require((caps.count("velocity") != 0) == (flow.velocityThresholdSI > 0)
            && (caps.count("velocitygradient") != 0) == (flow.velocityGradientThresholdSI > 0),
            "unused flow refinement threshold");
    }
    require(mesh.Nx > 0 && mesh.Ny > 0 && mesh.maxInitialIterations >= 0,
        "mesh.dict: positive background divisions and nonnegative initial refinement iterations required");
    const auto entity = [&](int id) -> const DictionaryGeometry::BoundaryEntity& {
        const auto found = g.boundaryById.find(id);
        require(found != g.boundaryById.end(), "boundary ID " + std::to_string(id) + " does not exist in geometry.dict/Boundary");
        require(!found->second.ignored, "boundary ID " + std::to_string(id) + " is ignored in geometry.dict");
        return found->second;
    };
    const auto& moving = entity(interfaceConfig.movingBoundaryId);
    if (mesh.runtimeAdaptation.enabled)
    {
        require(electric.model == "leakyDielectric" && electric.interfaceDiscretization == "sharpHermiteP3"
            && flow.outletIncomingPhaseId != -1, "runtime AMR requires conducting P3 and explicit open-boundary incoming material");
        require(mesh.geometryRefinement.enabled, "runtime AMR requires GeometryRefinement");
        std::set<int> phases;
        for (const auto& rule : mesh.geometryRefinement.phaseRules)
            require(rule.phase > 0 && rule.phase <= static_cast<int>(g.phaseName.size())
                && rule.minLevel >= 0 && rule.maxLevel >= rule.minLevel && phases.insert(rule.phase).second,
                "refinePhase requires unique valid phase IDs and 0 <= minLevel <= maxLevel");
        int movingRules = 0;
        for (const auto& rule : mesh.geometryRefinement.boundaryRules)
        {
            const auto& boundary = entity(rule.boundaryId);
            require(boundary.kind != GeometricBoundary::domain && !boundary.lineIds.empty()
                && std::isfinite(rule.cGamma) && rule.cGamma > 0 && rule.maxLevel >= 0,
                "runtime AMR requires resolved boundary polylines and explicit refinement levels");
            if (rule.boundaryId == interfaceConfig.movingBoundaryId) ++movingRules;
        }
        require(movingRules == 1, "runtime AMR requires exactly one refineBoundary rule for the moving liquid surface");
    }
    require(moving.kind == GeometricBoundary::moving && moving.levelSetIndex >= 0
        && moving.levelSetIndex < static_cast<int>(g.levelSetName.size())
        && g.levelSetName[moving.levelSetIndex] == interfaceConfig.fluidInterfaceLevelSet
        && moving.phase + 1 == flow.liquidPhase && moving.opposedPhase + 1 == flow.gasPhase,
        "Interface: moving boundary must select the named level set oriented liquid-to-gas");
    std::set<int> electricIds;
    for (const auto& bc : electric.domainBoundaries)
    {
        const auto& b = entity(bc.boundaryId);
        require(b.kind == GeometricBoundary::domain && (b.adjacentPhase + 1 == flow.liquidPhase
            || b.adjacentPhase + 1 == flow.gasPhase), "electricity.dict: electric domain BC must adjoin active fluid");
        require(electricIds.insert(bc.boundaryId).second, "duplicate electric boundary ID");
    }
    std::set<int> pressureIds;
    const auto pressure = [&](const FittedPressureBoundarySI& port, int adjacent) {
        const auto& b = entity(port.boundaryId);
        require(b.kind == GeometricBoundary::domain && b.adjacentPhase + 1 == adjacent
            && electricIds.count(port.boundaryId) && pressureIds.insert(port.boundaryId).second,
            "hydro.dict: pressure port " + std::to_string(port.boundaryId) + " needs the proper adjacent phase and explicit electric BC");
    };
    pressure(flow.supply, flow.liquidPhase);
    const auto& outlet = entity(flow.outlet.boundaryId);
    require(outlet.adjacentPhase + 1 == flow.liquidPhase || outlet.adjacentPhase + 1 == flow.gasPhase,
        "hydro.dict: outlet must adjoin liquid or gas");
    pressure(flow.outlet, outlet.adjacentPhase + 1);
    for (const auto& ambient : flow.ambient) pressure(ambient, flow.gasPhase);
    std::set<int> freeLines(moving.lineIds.begin(), moving.lineIds.end()), fixedLines;
    require(!freeLines.empty() && freeLines.size() == moving.lineIds.size(), "empty or repeated free-interface line IDs");
    for (const auto& b : g.boundaries)
    {
        if (b.ignored) continue;
        for (int id : b.lineIds) require(id >= 0 && id < static_cast<int>(g.lines.size()), "invalid geometric line ID");
        if (b.kind == GeometricBoundary::moving) require(b.id == interfaceConfig.movingBoundaryId, "additional moving boundary has no consumer");
        if (b.kind == GeometricBoundary::domain && b.adjacentPhase + 1 != electric.conductorPhaseId)
            require(electricIds.count(b.id) != 0, "missing electric BC for active domain boundary " + std::to_string(b.id));
        if (b.kind != GeometricBoundary::fixed) continue;
        for (int id : b.lineIds)
        {
            const auto& line = g.lines[id];
            require(line.type() == LineType::straight && (line.points[0].x == line.points[1].x
                || line.points[0].y == line.points[1].y), "current fitted solid walls must be axis-aligned straight segments");
            if (b.levelSetIndex == moving.levelSetIndex)
            {
                require((b.phase + 1 == electric.conductorPhaseId && b.opposedPhase + 1 == flow.liquidPhase)
                    || (b.opposedPhase + 1 == electric.conductorPhaseId && b.phase + 1 == flow.liquidPhase), "composite fixed wall must separate conductor and liquid");
                require(!freeLines.count(id) && fixedLines.insert(id).second, "overlapping fixed/moving line ownership");
            }
        }
    }
    std::vector<SignInt> chain;
    std::set<int> covered;
    for (auto line : g.linesOnInterface[moving.levelSetIndex])
    {
        require(covered.insert(line.value).second && (freeLines.count(line.value) || fixedLines.count(line.value)), "source level set has unbound/repeated physical lines");
        if (freeLines.count(line.value)) chain.push_back(line);
    }
    require(!chain.empty() && covered.size() == freeLines.size() + fixedLines.size(), "incomplete free/fixed line coverage");
    const auto endpoint = [&](SignInt line, bool end) { return g.lines[line.value].points[end != line.isNegative ? 1 : 0]; };
    const auto same = [](Vector2D a, Vector2D b) { return a.x == b.x && a.y == b.y; };
    for (std::size_t i = 1; i < chain.size(); ++i)
        require(same(endpoint(chain[i-1], true), endpoint(chain[i], false)), "free lines must form one directed open chain");
    const auto begin = endpoint(chain.front(), false), end = endpoint(chain.back(), true);
    require(!same(begin, end), "closed free-interface chains are not supported by this contact binding");
    for (int id : interfaceConfig.pinnedContactPointIds)
    {
        require(id > 0 && id <= static_cast<int>(g.points.size()), "pinned Point ID " + std::to_string(id) + " is outside geometry");
        const auto point = g.points[id-1];
        require(point.x > 0 && (same(point, begin) || same(point, end)), "pin must be an off-axis free-chain endpoint");
        bool onWall = false;
        for (int line : fixedLines) onWall |= same(point, g.lines[line].points[0]) || same(point, g.lines[line].points[1]);
        require(onWall, "pin must also be a fixed-wall endpoint");
    }
    for (const auto point : {begin, end})
    {
        bool pinned = false;
        for (int id : interfaceConfig.pinnedContactPointIds) pinned |= same(point, g.points[id-1]);
        require(pinned || point.x == 0, "free-chain endpoints must be pinned contacts or the RZ axis");
    }
}
