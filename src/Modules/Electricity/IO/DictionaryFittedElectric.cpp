#include "Modules/Electricity/IO/DictionaryFittedElectric.h"
#include "Core/IO/ConfigurationTrace.h"

#include <set>

FittedElectricCaseConfig readFittedElectricCase(const std::string& path, DictionaryInput::ConfigurationTrace& trace)
{
    using Node = DictionaryInput::Node;
    DictionaryInput::Document document(path);
    document.root().allow({"Electricity"});
    const auto& root = document.root().block("Electricity");
    root.allow({"Model", "Initial", "BoundaryConditions", "SurfaceCharge", "Numerics"});
    FittedElectricCaseConfig result;
    const auto& model = root.block("Model");
    model.allow({"model", "interfaceDiscretization"});
    const auto& choice = model.get("model");
    result.model = trace.record(choice, choice.choice({"perfectDielectric", "leakyDielectric"}), "model", "electric constitutive model");
    const auto& discretization = model.get("interfaceDiscretization");
    result.interfaceDiscretization = trace.record(discretization,
        discretization.choice({"cellHarmonic", "sharpHermiteP3"}), "model", "actual fitted electric interface discretization");
    const bool conducting = result.model == "leakyDielectric";
    if (result.interfaceDiscretization != (conducting ? "sharpHermiteP3" : "cellHarmonic"))
        discretization.fail("sharpHermiteP3 for conducting motion, cellHarmonic for the restored q=0 motion");
    const auto& initial = root.block("Initial");
    initial.allow({"surfaceChargeDensity"});
    const auto& charge = initial.get("surfaceChargeDensity");
    result.initialSurfaceChargeSI = trace.record(charge, charge.quantity(UnitParser::Dimension::SurfaceChargeDensity, "C/m^2"),
        "C/m^2", "fitted bootstrap free surface charge");
    if (result.initialSurfaceChargeSI != 0.0) charge.fail("zero initial free surface charge for the current bootstrap capability");

    const auto& boundaries = root.block("BoundaryConditions");
    boundaries.allow({"Conductor", "Domain"});
    const auto& conductor = boundaries.block("Conductor");
    conductor.allow({"phaseId", "potential"});
    const auto& phase = conductor.get("phaseId");
    result.conductorPhaseId = trace.record(phase, phase.integer(), "1-based ID", "fixed conducting solid geometry phase");
    if (result.conductorPhaseId <= 0) phase.fail("a positive conductor phase ID");
    const auto potential = [&](const Node& node, const char* consumer) {
        return trace.record(node, node.quantity(UnitParser::Dimension::Voltage, "V"), "V", consumer);
    };
    result.conductorPotentialSI = potential(conductor.get("potential"), "fitted solid electrode Dirichlet potential");
    std::set<int> ids;
    const auto& domain = boundaries.block("Domain");
    if (domain.entries.empty()) domain.fail("explicit domain electrical boundary conditions");
    for (const auto& block : domain.entries)
    {
        const auto& type = block.get("type");
        FittedElectricBoundarySI boundary;
        boundary.type = trace.record(type, type.choice({"dirichlet", "homogeneousNeumann"}), "model", "electric boundary operator");
        if (boundary.type == "dirichlet") block.allow({"boundaryId", "type", "potential"});
        else block.allow({"boundaryId", "type"});
        const auto& id = block.get("boundaryId");
        boundary.boundaryId = trace.record(id, id.integer(), "1-based ID", "electric domain boundary binding");
        if (boundary.boundaryId <= 0 || !ids.insert(boundary.boundaryId).second) id.fail("distinct positive electric boundary IDs");
        if (boundary.type == "dirichlet") boundary.potentialSI = potential(block.get("potential"), "domain electrode Dirichlet potential");
        result.domainBoundaries.push_back(boundary);
    }
    const auto tolerance = [&](const Node& entry, const char* consumer) {
        const double value = entry.number();
        if (value <= 0.0 || value >= 1.0) entry.fail("a relative acceptance tolerance in (0,1)");
        return trace.record(entry, value, "dimensionless", consumer);
    };
    const auto& surface = root.block("SurfaceCharge");
    const auto& surfaceModel = surface.get("model");
    result.surfaceChargeModel = trace.record(surfaceModel,
        conducting ? surfaceModel.choice({"explicitConductionMaterialTransport", "implicitConductionMaterialTransport"})
                   : surfaceModel.choice({"none"}), "model", "charge time integration on the shared material trajectory");
    const bool implicit = result.surfaceChargeModel == "implicitConductionMaterialTransport";
    if (implicit) surface.allow({"model", "closureTolerance", "Iteration"});
    else if (conducting) surface.allow({"model", "closureTolerance"});
    else surface.allow({"model"});
    if (conducting) result.chargeClosure = tolerance(surface.get("closureTolerance"), "panel bulk-current and time-consistent electrode charge closure");
    if (implicit)
    {
        const auto& iteration = surface.block("Iteration");
        iteration.allow({"underRelaxation", "maxIterations", "densityChangeTolerance", "relativeEquationTolerance"});
        Quantum::Electricity::SurfaceChargeRelaxationOptions options;
        const auto& omega=iteration.get("underRelaxation");
        options.underRelaxation=trace.record(omega,omega.number(),"dimensionless","surface-charge fixed-point under-relaxation");
        if (options.underRelaxation<=0.0 || options.underRelaxation>1.0) omega.fail("under-relaxation in (0,1]");
        const auto& limit=iteration.get("maxIterations");
        options.maxIterations=trace.record(limit,limit.integer(),"iterations","finite electric/charge iteration limit");
        if (options.maxIterations<=0) limit.fail("a positive iteration limit");
        const auto& change=iteration.get("densityChangeTolerance");
        options.densityChangeTolerance=trace.record(change,
            change.quantity(UnitParser::Dimension::SurfaceChargeDensity,"C/m^2"),"C/m^2",
            "surface-charge adjacent-iterate density change; explicitly normalized by App");
        if (options.densityChangeTolerance<=0.0) change.fail("a positive surface-charge density tolerance");
        options.relativeEquationTolerance=tolerance(iteration.get("relativeEquationTolerance"),"full material charge time-equation residual");
        result.relaxationSI=options;
    }
    const auto& numerics = root.block("Numerics");
    if (conducting) numerics.allow({"potentialBounds", "reconstructedPoisson", "gaussClosure"});
    else numerics.allow({"potentialBounds", "polarizationConsistency"});
    result.potentialBounds = tolerance(numerics.get("potentialBounds"), "charge-free bootstrap/update potential bounds");
    if (conducting)
    {
        result.fluxChecks.reconstructedPoisson = tolerance(numerics.get("reconstructedPoisson"), "P3 reconstructed cell displacement balance");
        result.fluxChecks.gaussClosure = tolerance(numerics.get("gaussClosure"), "P3 global Gauss closure with original machine floor");
    }
    else result.polarizationConsistency = tolerance(numerics.get("polarizationConsistency"), "q=0 polarization field/flux consistency");
    return result;
}
