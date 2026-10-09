#include "Electrospray/IO/ElectrosprayNumericsReport.h"
#include "Electrospray/Compatibility/ElectrosprayConfig.h"
#include "Electrospray/Compatibility/ElectrosprayPhysicalConfig.h"
#include "Electrospray/Setup/ElectrosprayWorkflow.h"
#include "Core/IO/OutputLayout.h"
#include "Solvers/SolversConfig.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace
{
void quote(std::ostream& out, const std::string& value)
{
    out << '"';
    for (unsigned char c : value)
    {
        if (c == '"' || c == '\\') out << '\\' << static_cast<char>(c);
        else if (c < 32) out << "\\u00" << "0123456789abcdef"[c >> 4] << "0123456789abcdef"[c & 15];
        else out << static_cast<char>(c);
    }
    out << '"';
}
}

void writeElectrosprayNumericsReport(const ElectrosprayPhysicalConfig& physical,
    const ElectrosprayConfig& run, const ElectrosprayWorkflow& workflow,
    const Solvers::SolversConfig* legacySharedSolvers)
{
    const auto path = OutputLayout::meta() / (run.outputName + "_numerics.json");
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot open numerical configuration report: " + path.string());
    const auto& n = physical.numerics;
    const bool charged = workflow.conducting();
    const bool moving = workflow.movesInterface();
    const bool polarization = workflow.usesDielectricPolarization();
    out << std::setprecision(17) << std::boolalpha;
    out << "{\n  \"reportVersion\": 1,\n  \"inputSchema\": \"legacy\",\n"
        "  \"scope\": \"effective numerical controls only; not complete resolved physical configuration\",\n"
        "  \"executionStatus\": \"not started; actual solver flags and accepted steps remain in run diagnostics\",\n"
        "  \"physicsMode\": ";
    quote(out, run.physicsMode);
    out << ",\n  \"requiredLinearBackend\": \"SuperLU\",\n  \"parameters\": [\n";
    bool first = true;
    const auto entry = [&](const char* name, double value, const char* consumer, bool used, double multiplier = 1.0)
    {
        if (!std::isfinite(value)) throw std::runtime_error("nonfinite numerical configuration: " + std::string(name));
        if (!first) out << ",\n";
        first = false;
        out << "    {\"name\": "; quote(out, name);
        out << ", \"value\": " << value << ", \"units\": \"dimensionless\", \"source\": ";
        quote(out, physical.legacyToleranceSource);
        out << ", \"sourceMultiplier\": " << multiplier << ", \"defaulted\": false, \"usedBySelectedMode\": " << used
            << ", \"consumer\": "; quote(out, consumer); out << '}';
    };
    entry("electricPotentialLinear", n.electricPotentialLinear, "bootstrap cellHarmonic; q=0 update or conducting fitted P3 Poisson", true);
    entry("momentumLinear", n.flow.momentumLinear, "fitted NS momentum predictor SuperLU residual", true);
    entry("pressureLinear", n.flow.pressureLinear, "fitted NS pressure correction SuperLU residual", true);
    entry("coupledLinear", n.flow.coupledLinear, "fitted NS coupled velocity-pressure correction SuperLU residual", true);
    entry("volumeCorrectionLinear", n.volumeCorrectionLinear, "fitted interface-volume correction linear subproblem", moving);
    entry("momentumBalance", n.flow.momentumBalance, "fitted full CN momentum balance; original load-norm denominator", true);
    entry("localContinuity", n.flow.localContinuity, "fitted cell flux divided by original local absolute flux scale", true);
    const bool pressureOutlet = physical.outletIncomingPhaseId < 0;
    entry("outletAdmissibility", n.flow.outletAdmissibility, "pressure outlet without incoming material only; retains machine-roundoff impulse term", pressureOutlet);
    entry("phaseFluxBalance", n.phaseFluxBalance, "moving phase ledger assembly and pressure-port consistency", moving);
    entry("potentialBounds", n.potentialBounds, "charge-free bootstrap/update Dirichlet maximum principle", true);
    entry("reconstructedPoisson", n.electricFlux.reconstructedPoisson, "P3 independently reconstructed cell displacement balance", charged);
    entry("gaussClosure", n.electricFlux.gaussClosure, "P3 global Gauss closure; effective gate=max(value,128*epsilon*cellCount)", charged);
    entry("polarizationConsistency", n.polarizationConsistency, "q=0 polarization cell/face displacement consistency; original roundoff floor", polarization);
    entry("chargeClosure", n.chargeClosure, "old-current material remap and electrode charge closure; original absolute-Q/current scale", charged, 16.0);
    out << "\n  ],\n  \"unchangedNonlinearVolumeCriteria\": {\"relativeTransportBalance\": 0.001, \"cellCapacityNormalized\": 1e-10},\n"
        "  \"unchangedChargeRelaxationCflLimit\": 0.25,\n"
        "  \"maxwellReconstruction\": \"no linear solver tolerance consumed; retains existing local trace/traction arithmetic gates\",\n"
        "  \"legacySharedPoissonCellCentered\": {\"usedByThisApp\": false, \"loadedObjectPresent\": " << (legacySharedSolvers != nullptr);
    if (legacySharedSolvers)
    {
        const auto& requested = legacySharedSolvers->poissonCellCentered;
        out << ", \"useSuperLU\": " << requested.useSuperLU << ", \"useCG\": " << requested.useCG
            << ", \"useIC\": " << requested.useIC << ", \"maxIter\": " << requested.maxIter << ", \"relTol\": ";
        // The old shared reader is not strict. Do not emit NaN/Inf as invalid
        // JSON, nor promote an unused legacy value to a physical input gate.
        if (std::isfinite(requested.relTol)) out << requested.relTol;
        else out << "null";
        out << ", \"relTolFinite\": " << std::isfinite(requested.relTol) << ", \"verbose\": " << requested.verbose;
    }
    out << "}\n}\n";
    out.close();
    if (!out) throw std::runtime_error("failed to finish numerical configuration report: " + path.string());
    std::cout << "[电喷雾] 数值设置已展开：已启用线性系统的后端固定为 SuperLU；"
        "线性残差与动量、连续性、电通量、电荷收支分别检查。\n"
        "[电喷雾] 兼容输入提示：本入口不消费共享 solvers.dict/PoissonCellCentered 选项；"
        "当前各容限来自旧 ElectrosprayPhysical.linearResidualTolerance 的显式映射，不要通过修改共享 CG/IC 开关调节喷嘴。\n"
        "[电喷雾] 生效数值设置：" << path.string() << "；这是求解前配置，不表示算例已通过。\n";
}
