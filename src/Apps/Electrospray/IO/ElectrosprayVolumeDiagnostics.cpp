#include "Electrospray/IO/ElectrosprayVolumeDiagnostics.h"
#include "Modules/Hydrodynamics/HydroFittedDomain.h"
#include "Modules/Interface/Compatibility/InterfaceVolumeCorrection.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <initializer_list>
#include <iostream>
#include <locale>
#include <sstream>
#include <stdexcept>

namespace
{
void requireFinite(std::initializer_list<long double> values, const char* description)
{
    for (const auto value : values)
        if (!std::isfinite(value)) throw std::runtime_error(description);
}

std::ofstream openOutput(const std::string& path, std::ios_base::openmode mode)
{
    std::ofstream out(path, mode);
    if (!out) throw std::runtime_error("cannot open local liquid volume diagnostics: " + path);
    out.imbue(std::locale::classic());
    out << std::setprecision(17);
    return out;
}

void finishOutput(std::ofstream& out, const std::string& path)
{
    out.flush();
    if (!out) throw std::runtime_error("failed writing local liquid volume diagnostics: " + path);
    out.close();
    if (!out) throw std::runtime_error("failed closing local liquid volume diagnostics: " + path);
}
}

ElectrosprayVolumeDiagnostics::ElectrosprayVolumeDiagnostics(
    const std::string& seriesPrefix, double lengthSI)
    : seriesPath_(seriesPrefix + "_local_volume_diagnostics.dat"), lengthSI_(lengthSI),
      volumeFactor_(static_cast<long double>(lengthSI) * lengthSI * lengthSI)
{
    if (seriesPrefix.empty() || !std::isfinite(lengthSI_) || lengthSI_ <= 0.0
        || !std::isfinite(volumeFactor_) || volumeFactor_ <= 0.0L)
        throw std::invalid_argument("invalid local volume diagnostic prefix or SI length scale");
    history_.emplace(seriesPath_, "Local liquid volume correction",
        std::vector<std::string>{"step","time_s","iterations","backtracks","rows","columns","used_superlu",
            "predictor_relative_balance","corrected_relative_balance","max_cell_capacity_relative_residual",
            "sum_absolute_residual_m3","signed_residual_m3","max_phi_correction_m","max_gradient_correction",
            "predictor_relative_balance_defined"},
        std::map<std::string,std::string>{{"UndefinedPredictor", "defined=0 means positive defect / zero transport (legacy +inf); value 0 is ONLY a display placeholder; exclude it from ratio plots"},
            {"State", "volume correction diagnostic; not alone a completed coupled step or complete field frame"}});
}

void ElectrosprayVolumeDiagnostics::writeStep(
    int step, double timeSI, const InterfaceVolumeCorrection::Report& report)
{
    if (step < 1 || !std::isfinite(timeSI) || timeSI < 0.0)
        throw std::invalid_argument("invalid local volume diagnostic step or SI time");
    const long double sumAbsoluteSI = report.final.sumAbsolute * volumeFactor_;
    const long double signedResidualSI = report.final.signedSum * volumeFactor_;
    const double maxCorrectionSI = report.maxCorrection * lengthSI_;
    requireFinite({report.initial.sumAbsolute, report.initial.transportedMagnitude,
        report.final.relativeBalance,
        report.final.maxRelativeCell, sumAbsoluteSI, signedResidualSI, maxCorrectionSI,
        report.maxGradientCorrection},
        "nonfinite local volume step diagnostics");
    // A positive defect divided by zero initial transport has no finite
    // relative diagnostic. Preserve its meaning with an explicit validity
    // flag; Tecplot POINT data cannot contain the solver's +inf sentinel.
    const bool undefinedInitialRatio = std::isinf(report.initial.relativeBalance)
        && report.initial.relativeBalance > 0.0L
        && report.initial.transportedMagnitude == 0.0L && report.initial.sumAbsolute > 0.0L;
    if (!std::isfinite(report.initial.relativeBalance) && !undefinedInitialRatio)
        throw std::runtime_error("invalid initial local volume relative balance diagnostic");
    history_->append(TecplotWriter::row(step, timeSI, report.iterations,
        report.backtracks, report.rows, report.columns, report.usedSuperLU,
        undefinedInitialRatio ? 0.0L : report.initial.relativeBalance,
        report.final.relativeBalance, report.final.maxRelativeCell,
        sumAbsoluteSI, signedResidualSI, maxCorrectionSI, report.maxGradientCorrection, !undefinedInitialRatio));
}

void ElectrosprayVolumeDiagnostics::writeProgress(const InterfaceVolumeCorrection::Progress& progress)
{
    using Stage = InterfaceVolumeCorrection::Progress::Stage;
    const auto& report = progress.report;
    std::ostringstream message;
    message.imbue(std::locale::classic());
    message << std::scientific << std::setprecision(6) << "[电喷雾] 局部体积修正：";
    const auto writeResidual = [&] {
        message << "绝对残差总和（计算体积单位）=" << report.final.sumAbsolute
            << "，相对输运收支=" << report.final.relativeBalance
            << "，最大单元容量相对残差=" << report.final.maxRelativeCell;
    };
    const auto writeTrial = [&] {
        message << "迭代 " << report.iterations << "，回退 " << progress.backtrack
            << "，步长系数=" << progress.damping;
    };
    switch (progress.stage)
    {
    case Stage::PredictorReady:
        message << "预测液面检查完成；";
        writeResidual();
        break;
    case Stage::LinearStepReady:
        message << "迭代 " << report.iterations << "，线性修正完成；约束=" << report.rows
            << "，自由度=" << report.columns << "，SuperLU=" << (report.usedSuperLU ? "是" : "否")
            << "，累计最大线性φ修正（计算长度单位）=" << report.maxLinearCorrection
            << "，累计最大线性梯度修正=" << report.maxLinearGradientCorrection;
        break;
    case Stage::TrialStarted:
        writeTrial();
        message << "；开始重建候选液面几何，并计算体积与相通量。";
        break;
    case Stage::TrialRejected:
        writeTrial();
        // final may still describe the last accepted iterate, not this failed
        // geometry. Only the solver's failure detail describes this rejection.
        message << "；候选未接受，原因：" << progress.detail;
        break;
    case Stage::TrialAccepted:
        writeTrial();
        message << "；保留内迭代候选（尚非时间步接受），";
        writeResidual();
        if (!progress.detail.empty()) message << "；目标体积暂不可行，须继续修正：" << progress.detail;
        break;
    case Stage::Converged:
        message << "本步迭代达到既定准则，迭代=" << report.iterations
            << "，累计回退=" << report.backtracks << "；";
        writeResidual();
        message << "，最大φ修正（计算长度单位）=" << report.maxCorrection
            << "，最大梯度修正=" << report.maxGradientCorrection
            << "。不替代原全局体积诊断与物理验收。";
        break;
    }
    message << '\n';
    std::cout << message.str() << std::flush;
    if (!std::cout) throw std::runtime_error("failed writing local volume progress diagnostics");
}

void ElectrosprayVolumeDiagnostics::writeFields(const std::string& fieldPrefix,
    const HydroFittedDomain& domain, const InterfaceVolumeMeasure& oldStored,
    const InterfacePhaseFlux::FluxBudget& oldFlux,
    const InterfaceVolumeMeasure& newStored,
    const InterfacePhaseFlux::FluxBudget& newFlux, double dtCode)
{
    if (fieldPrefix.empty() || !std::isfinite(dtCode) || dtCode <= 0.0)
        throw std::invalid_argument("invalid local volume field prefix or code timestep");
    if (oldStored.index() != newStored.index())
        throw std::invalid_argument("local volume diagnostic scalar representation changed");
    std::visit([&](const auto& oldMeasure, const auto& newMeasure) {
        // Only active native cells are exported. Inactive entries may retain NaN;
        // checking their values would accidentally impose a different solid mask.
        for (int k : domain.cellIndexByRow())
        {
            if (k < 0 || static_cast<std::size_t>(k) >= oldMeasure.negativeVolumeByCell.size()
                || static_cast<std::size_t>(k) >= newMeasure.negativeVolumeByCell.size()
                || static_cast<std::size_t>(k) >= newMeasure.fixedBoundary.capacityByCell.size()
                || static_cast<std::size_t>(k) >= oldFlux.outwardFluxByCell.size()
                || static_cast<std::size_t>(k) >= newFlux.outwardFluxByCell.size())
                throw std::invalid_argument("local volume diagnostics lost native-cell indexing");
            requireFinite({oldMeasure.negativeVolumeByCell[k], newMeasure.negativeVolumeByCell[k],
                newMeasure.fixedBoundary.capacityByCell[k], oldFlux.outwardFluxByCell[k],
                newFlux.outwardFluxByCell[k]}, "nonfinite active-cell volume diagnostics");
        }
        const std::string path = fieldPrefix + "_volume_balance_cells.csv";
        auto out = openOutput(path, std::ios::out);
        out << "cell_index,capacity_m3,old_liquid_volume_m3,new_liquid_volume_m3,trapezoidal_outward_transport_m3,residual_m3\n";
        for (int k : domain.cellIndexByRow())
        {
            const auto oldVolume = oldMeasure.negativeVolumeByCell[k];
            const auto newVolume = newMeasure.negativeVolumeByCell[k];
            const auto transport = 0.5L * dtCode * (oldFlux.outwardFluxByCell[k] + newFlux.outwardFluxByCell[k]);
            const auto capacitySI = newMeasure.fixedBoundary.capacityByCell[k] * volumeFactor_;
            const auto oldVolumeSI = oldVolume * volumeFactor_;
            const auto newVolumeSI = newVolume * volumeFactor_;
            const auto transportSI = transport * volumeFactor_;
            const auto residualSI = ((newVolume - oldVolume) + transport) * volumeFactor_;
            requireFinite({capacitySI, oldVolumeSI, newVolumeSI, transportSI, residualSI},
                "local volume field diagnostic SI conversion overflow");
            out << k << ',' << capacitySI << ',' << oldVolumeSI << ',' << newVolumeSI
                << ',' << transportSI << ',' << residualSI << '\n';
        }
        finishOutput(out, path);
    }, oldStored, newStored);
}
