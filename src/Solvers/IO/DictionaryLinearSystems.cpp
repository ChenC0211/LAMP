#include "Solvers/IO/DictionaryLinearSystems.h"
#include "Core/IO/ConfigurationTrace.h"

#include <algorithm>

namespace Solvers
{
std::map<std::string, DirectLinearSystemInput> readDirectLinearSystems(
    const std::string& path, const std::vector<std::string>& requiredSystems,
    DictionaryInput::ConfigurationTrace& trace)
{
    DictionaryInput::Document document(path);
    document.root().allow({"Solvers"});
    const auto& systems = document.root().block("Solvers");
    std::map<std::string, DirectLinearSystemInput> result;
    for (const auto& block : systems.entries)
    {
        if (std::find(requiredSystems.begin(), requiredSystems.end(), block.name) == requiredSystems.end())
            block.fail("a linear system actually used by the selected workflow");
        block.allow({"backend", "relativeTolerance"});
        DirectLinearSystemInput input;
        const auto& backend = block.get("backend");
        input.backend = trace.record(backend, backend.choice({"SuperLU"}), "backend", block.name + " direct linear solve");
        const auto& tolerance = block.get("relativeTolerance");
        input.relativeTolerance = trace.record(tolerance, tolerance.number(), "dimensionless", block.name + " actual algebraic residual");
        if (input.relativeTolerance <= 0.0 || input.relativeTolerance >= 1.0)
            tolerance.fail("a finite relative linear tolerance in (0,1)");
        result.emplace(block.name, input);
    }
    for (const auto& name : requiredSystems) systems.block(name);
    return result;
}
}
