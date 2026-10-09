#include "Electrospray/IO/ElectrosprayMeshDiagnostics.h"
#include "Electrospray/Coupling/ElectrosprayMeshAdaptation.h"
#include "Core/Units/UnitsContext.h"

#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <locale>
#include <sstream>
#include <stdexcept>

ElectrosprayMeshDiagnostics::ElectrosprayMeshDiagnostics(const std::string& prefix,
    const ElectrosprayHydroScales& scale, const UnitsContext& units)
    : path_(prefix + "_amr.dat"), scale_(scale),
      chargeScale_(8.8541878128e-12 * units.scales().Phi0_V * scale.length)
{
    if (units.mode() != UnitsMode::kSI || prefix.empty() || std::filesystem::exists(path_))
        throw std::invalid_argument("AMR 收支输出要求 SI 单位和本次运行的新文件路径");
    history_.emplace(path_, "Accepted AMR events at unchanged physical time",
        std::vector<std::string>{"accepted_step","time_s","old_cells","new_cells","refined_source_cells","coarsened_source_cells","protected_cells","sweeps",
            "delta_volume_m3","volume_roundoff_bound_m3","delta_charge_C","delta_liquid_boundary_flux_m3_s",
            "current_max_div_1_s","previous_max_div_1_s","current_max_velocity_correction_m_s","previous_max_velocity_correction_m_s",
            "delta_axial_momentum_kg_m_s","delta_pressure_integral_Pa_m3","wall_seconds"});
}

void ElectrosprayMeshDiagnostics::progress(ElectrosprayMeshStage stage, int cells)
{
    const char* message = stage == ElectrosprayMeshStage::Grid ? "候选网格已生成，转移完整液面与电荷历史"
        : stage == ElectrosprayMeshStage::InterfaceAndCharge ? "液面与电荷转移完成，恢复流场连续性"
        : stage == ElectrosprayMeshStage::Flow ? "流场转移完成，重算带电电场及毛细/Maxwell 载荷"
        : stage == ElectrosprayMeshStage::ElectricAndForces ? "电场与受力重建完成，核对液体体积及相通量"
        : "候选状态准备完成，等待统一接受";
    std::cout << "[电喷雾 AMR] " << message << "；叶单元数=" << cells << "。\n" << std::flush;
}

void ElectrosprayMeshDiagnostics::flowPolicy(double speedThresholdSI, double gradientThresholdSI)
{
    std::ostringstream message;
    message << std::setprecision(6) << "[电喷雾 AMR] 流场判据从首个已接受步开始生效：";
    if (speedThresholdSI > 0) message << "速度门限=" << speedThresholdSI << " m/s；";
    if (gradientThresholdSI > 0) message << "速度梯度门限=" << gradientThresholdSI << " 1/s；";
    message << "初始网格只按界面、液体和固定边界细化。\n";
    std::cout << message.str();
}

void ElectrosprayMeshDiagnostics::write(int step, double time, const ElectrosprayMeshReport& r)
{
    const long double volume = static_cast<long double>(scale_.length) * scale_.length * scale_.length;
    const long double volumeChange = (r.volumeAfter - r.volumeBefore) * volume;
    const long double chargeChange = (r.chargeAfter - r.chargeBefore) * chargeScale_;
    const long double fluxChange = (r.liquidFluxAfter - r.liquidFluxBefore)
        * scale_.velocity * scale_.length * scale_.length;
    const long double momentumChange = (r.hydro.current.targetMomentumZ - r.hydro.current.sourceMomentumZ)
        * scale_.density * scale_.velocity * volume;
    const long double pressureChange = (r.hydro.targetPressureIntegral - r.hydro.sourcePressureIntegral)
        * scale_.pressure * volume;
    if (step < 1 || !std::isfinite(time) || time <= 0 || r.oldCells <= 0 || r.newCells <= 0)
        throw std::invalid_argument("AMR 记录缺少有效的接受步与网格数量");
    for (long double value : {volumeChange, chargeChange, fluxChange, momentumChange, pressureChange})
        if (!std::isfinite(value)) throw std::runtime_error("AMR 收支 SI 换算产生非有限值");
    history_->append(TecplotWriter::row(step, time * scale_.time, r.oldCells, r.newCells,
        r.refinedSourceCells, r.coarsenedSourceCells, r.protectedCells, r.sweeps,
        volumeChange, r.volumeRoundoffBound * volume, chargeChange, fluxChange,
        r.hydro.current.maxDivergenceAfter / scale_.time, r.hydro.previous.maxDivergenceAfter / scale_.time,
        r.hydro.current.maxCorrection * scale_.velocity, r.hydro.previous.maxCorrection * scale_.velocity,
        momentumChange, pressureChange, r.seconds));
    std::ostringstream message;
    message << std::scientific << std::setprecision(4) << "[电喷雾 AMR] 已接受：第 " << step
        << " 步之后，物理时间不变；单元数 " << r.oldCells << " → " << r.newCells
        << "，被细化/粗化的原单元=" << r.refinedSourceCells << '/' << r.coarsenedSourceCells
        << "；液体体积变化=" << volumeChange << " m^3，电荷变化=" << chargeChange << " C。\n";
    std::cout << message.str() << std::flush;
}
