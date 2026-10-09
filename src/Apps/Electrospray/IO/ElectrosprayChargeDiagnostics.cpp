#include "Electrospray/IO/ElectrosprayChargeDiagnostics.h"
#include "Electrospray/Coupling/ElectrosprayConducting.h"
#include "Core/Units/UnitsContext.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <locale>
#include <sstream>
#include <stdexcept>

namespace
{
double converted(double value, double factor)
{
    const long double product = static_cast<long double>(value) * factor;
    const double result = static_cast<double>(product);
    if (!std::isfinite(value) || !std::isfinite(product) || !std::isfinite(result)
        || (value != 0.0 && result == 0.0))
        throw std::runtime_error("unrepresentable conducting charge diagnostic SI value");
    return result;
}

}

ElectrosprayChargeDiagnostics::ElectrosprayChargeDiagnostics(const std::string& prefix, const UnitsContext& units, bool implicit)
    : path_(prefix + "_conducting_charge_diagnostics.dat"), implicit_(implicit)
{
    constexpr double epsilon0 = 8.8541878128e-12;
    if (units.mode() != UnitsMode::kSI)
        throw std::invalid_argument("导电电荷收支输出要求使用 SI 单位");
    if (prefix.empty())
        throw std::invalid_argument("导电电荷收支输出缺少文件前缀");
    if (std::filesystem::exists(path_))
        throw std::invalid_argument("电荷收支文件已存在，拒绝重复创建：" + path_
            + "；重新计算请在启动时确认覆盖旧输出，不能直接追加到旧结果。");
    length_ = units.scales().L0_m; time_ = units.scales().T0_s;
    charge_ = epsilon0 * units.scales().Phi0_V * length_;
    density_ = epsilon0 * units.scales().Phi0_V / length_;
    current_ = charge_ / time_;
    for (double scale : {length_,time_,charge_,density_,current_})
        if (!std::isfinite(scale) || scale <= 0.0) throw std::invalid_argument("invalid conducting charge SI scales");
    std::vector<std::string> columns{"step","time_s","current_source_time_s","dt_s","panel_count","total_charge_before_C","total_charge_after_C",
        "old_boundary_outward_current_A","new_boundary_outward_current_A","integrated_normal_current_A",
        "charge_boundary_closure_C","charge_closure_scale_C","relative_charge_closure","absolute_unassigned_bulk_current_A",
        "signed_unassigned_bulk_current_A","max_panel_mean_abs_charge_Cpm2","relaxation_CFL","remapped",
        "remap_charge_residual_C","max_material_displacement_m","max_projection_distance_m","new_electric_superlu","new_electric_residual"};
    if (implicit_) columns.insert(columns.end(), {"charge_temporal_order","charge_iterations","charge_alpha","charge_beta",
        "charge_max_density_change_Cpm2","charge_relative_equation_residual","electric_geometry_reuses","electric_factor_reuses"});
    history_.emplace(path_, "Accepted surface charge balance", columns,
        std::map<std::string,std::string>{{"charge_model", implicit_ ? "implicitConductionMaterialTransport" : "explicitConductionMaterialTransport"},
            {"CurrentMeaning", "boundary conduction current, not emission current; use current_source_time_s for its time layer"}});
}

void ElectrosprayChargeDiagnostics::writeStep(int step, double timeCode,
    const ElectrosprayConductingState& state, const DielectricConductionReport& previousCurrent,
    const Quantum::Electricity::SurfaceChargeFieldStep& update, const ElectrosprayConductingStepReport& report)
{
    if (step < 1 || !state.electrical || state.electrical.completedStep()!=step || state.electrical.time()!=timeCode
        || !state.electrical.electric().report.success || !state.electrical.electric().report.conductionCurrent
        || !update.transport.success || !std::isfinite(timeCode) || timeCode <= 0.0
        || !std::isfinite(report.dt) || report.dt <= 0.0 || timeCode < report.dt
        || implicit_!=report.relaxation.has_value() || (implicit_ && !report.relaxation->success))
        throw std::invalid_argument("conducting diagnostic requires a completed state and its actual step report");
    double peak = 0.0;
    for (double value : state.electrical.charge()->transport().density())
    {
        if (!std::isfinite(value)) throw std::runtime_error("nonfinite accepted panel mean charge density");
        peak = std::max(peak,std::abs(value));
    }
    const auto& q = update.transport;
    const double total = converted(q.totalChargeAfter,charge_);
    const double residual = converted(report.chargeClosureResidual,charge_);
    const double maximumDensity = converted(peak,density_);
    const double newCurrent = converted(state.electrical.electric().report.conductionCurrent->netBoundaryCurrent,current_);
    auto values = TecplotWriter::row(step, converted(timeCode,time_), converted(implicit_ ? timeCode : timeCode-report.dt,time_),
        converted(report.dt,time_), state.electrical.charge()->geometry().panels().size(),
        converted(q.totalChargeBefore,charge_), total,
        converted(previousCurrent.netBoundaryCurrent,current_), newCurrent,
        converted(q.integratedNormalCurrent,current_), residual,
        converted(report.chargeClosureScale,charge_), report.relativeChargeClosure,
        converted(update.absoluteBulkCurrent,current_), converted(update.unassignedBulkCurrent,current_),
        maximumDensity, report.relaxationCfl, update.remapped,
        converted(update.remapConservationResidual,charge_),
        converted(update.maxMaterialDisplacement,length_), converted(update.maxProjectionDistance,length_),
        state.electrical.electric().report.usedSuperLU, state.electrical.electric().report.relativeResidual);
    if (implicit_)
    {
        const auto& iteration=*report.relaxation;
        const auto iterationValues = TecplotWriter::row(iteration.temporalOrder, iteration.iterations,
            iteration.alpha, iteration.beta, converted(iteration.maxDensityChange,density_),
            iteration.relativeEquationResidual, report.electricGeometryReuses, report.electricFactorReuses);
        values.insert(values.end(), iterationValues.begin(), iterationValues.end());
    }
    history_->append(values);
    std::ostringstream message;
    message << std::scientific << std::setprecision(6)
        << "[电喷雾] 导电步 " << step << " 已接受：液面总电荷=" << total
        << " C，最大面均电荷密度=" << maximumDensity << " C/m^2，电极收支残差=" << residual << " C。\n";
    if (implicit_)
    {
        message << "  电荷—电场步内迭代收敛：" << report.relaxation->iterations
            << " 次，时间式阶数=" << report.relaxation->temporalOrder
            << "，电荷方程相对残差=" << report.relaxation->relativeEquationResidual
            << "；源电流与本步电荷同一时刻。边界总流入电流=" << -newCurrent << " A，不是发射电流。\n";
    }
    else
        message << "  本步源电流取自上一时刻；当前新边界总流入电流=" << -newCurrent
            << " A，供下一步使用；不是发射电流。\n";
    std::cout << message.str() << std::flush;
}

void reportElectrosprayChargeIteration(int iteration, int limit, bool geometryReused, bool factorReused)
{
    if (iteration==1 || iteration%25==0)
        std::cout << "[电喷雾] 电荷—电场步内迭代 " << iteration << '/' << limit
            << "：本次 P3 电场及导电电流求解完成，几何拟合" << (geometryReused ? "复用" : "新建")
            << "，SuperLU 分解" << (factorReused ? "复用" : "新建") << "。\n" << std::flush;
}
