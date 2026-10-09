#include "Electrospray/IO/ElectrosprayHydroDiagnostics.h"
#include "Modules/Hydrodynamics/FittedFlowCaseConfig.h"
#include "Electrospray/Setup/ElectrosprayMaterials.h"
#include "Electrospray/Coupling/ElectrosprayInterfaceMotion.h"
#include "Modules/Hydrodynamics/Hydrodynamics.h"
#include "Modules/Hydrodynamics/HydroFittedDomain.h"
#include "Core/Mesh/QuadtreeMesh.h"

#include <iomanip>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace { using Domain = HydroFittedDomain; }

void writeElectrosprayHydroFields(const QuadtreeMesh& mesh, const Domain& domain, const Hydrodynamics& hydro,
    const std::vector<int>& phaseByCell, const ElectrosprayHydroScales& scale, const std::string& prefix)
{
    const auto state = hydro.state();
    std::ofstream cells(prefix + "_hydro_cells.csv"), faces(prefix + "_hydro_faces.csv");
    if (!cells || !faces) throw std::runtime_error("cannot open fitted Hydro fields");
    cells << std::setprecision(17) << "cell_index,cell_label,level,r_m,z_m,phase_id,active,density_kgpm3,viscosity_Pas,pressure_Pa,speed_mps,divergence_per_s\n";
    faces << std::setprecision(17) << "native_row,edge_label,horizontal,velocity_row,boundary_id,boundary_type,r_m,z_m,velocity_mps,outward_volume_flux_m3ps\n";
    for (int i = 0; i < mesh.cellIndexNum(); ++i)
    {
        const auto& cell = mesh.cells(mesh.cellIndexToLabel[i]);
        const bool active = domain.cellRowByIndex()[i] >= 0;
        cells << i << ',' << cell.label() << ',' << cell.level() << ',' << cell.centerCoord().x * scale.length
            << ',' << cell.centerCoord().y * scale.length << ',' << phaseByCell[i] << ',' << active
            << ',' << state.density[i] * scale.density << ',' << state.viscosity[i] * scale.viscosity
            << ',' << state.physicalPressure[i] * scale.pressure << ',' << state.speed[i] * scale.velocity
            << ',' << state.divergence[i] / scale.time << '\n';
    }
    for (const auto& face : domain.faces())
    {
        const double velocity = face.horizontal ? state.velocityY[face.nativeRow - mesh.vEdgeIndexNum()]
            : state.velocityX[face.nativeRow];
        faces << face.nativeRow << ',' << face.edgeLabel << ',' << face.horizontal << ',' << face.velocityRow
            << ',' << face.boundary.id << ',' << static_cast<int>(face.boundary.type) << ',' << face.center.x * scale.length
            << ',' << face.center.y * scale.length << ',' << velocity * scale.velocity
            << ',' << face.outwardSign * velocity * face.area * scale.velocity * scale.length * scale.length << '\n';
    }
    cells.close(); faces.close();
    if (!cells || !faces) throw std::runtime_error("incomplete fitted Hydro field output");
}

void writeElectrosprayPhaseFluxBudget(const Domain& domain, const InterfacePhaseFlux::FluxBudget& budget,
    const ElectrosprayHydroScales& scale, const std::string& prefix)
{
    std::ofstream cells(prefix + "_phase_flux_cells.csv"), faces(prefix + "_phase_flux_faces.csv");
    if (!cells || !faces) throw std::runtime_error("cannot open fitted phase flux ledger");
    const long double factor = static_cast<long double>(scale.velocity)*scale.length*scale.length;
    cells << std::setprecision(17) << "cell_index,liquid_outward_Q_m3ps\n";
    faces << std::setprecision(17)
        << "native_row,negative_cell,positive_cell,horizontal,liquid_positive_Q_m3ps\n";
    for (int k : domain.cellIndexByRow()) cells << k << ',' << budget.outwardFluxByCell[k]*factor << '\n';
    for (const auto& face : domain.faces())
        if (face.negativeCell >= 0 || face.positiveCell >= 0)
            faces << face.nativeRow << ',' << face.negativeCell << ',' << face.positiveCell << ','
                << face.horizontal << ',' << budget.positiveFluxByNativeFace[face.nativeRow]*factor << '\n';
    cells.close(); faces.close();
    if (!cells || !faces) throw std::runtime_error("incomplete fitted phase flux ledger output");
}

void reportElectrosprayHydroDomain(const HydroFittedDomain& domain)
{
    std::cout << "[电喷雾] 轴对称流体区域完成：活动单元=" << domain.cellIndexByRow().size()
        << "，自由面速度未知量=" << domain.nativeByVelocityRow().size()
        << "；固壁无滑移，压力开口已指定回流相物性。\n";
}

ElectrosprayBoundaryFlowRatesSI evaluateElectrosprayBoundaryFlow(const QuadtreeMesh& mesh,
    const HydroFittedDomain& domain, const Hydrodynamics& hydro,
    const FittedFlowCaseConfig& flowConfig, const ElectrosprayHydroScales& scale)
{
    ElectrosprayBoundaryFlowRatesSI result;
    const auto state = hydro.state();
    for (const auto& face : domain.faces())
    {
        if (face.outwardSign == 0.0) continue;
        const double velocity = face.horizontal ? state.velocityY[face.nativeRow - mesh.vEdgeIndexNum()]
            : state.velocityX[face.nativeRow];
        const double flux = face.outwardSign * face.area * velocity * scale.length * scale.length * scale.velocity;
        result.netOutward += flux;
        if (face.boundary.id == flowConfig.supply.boundaryId) result.supplyOutward += flux;
        else if (face.boundary.id == flowConfig.outlet.boundaryId) result.outletOutward += flux;
        else if (face.boundary.type == Domain::BoundaryType::PressureReservoir) result.ambientOutward += flux;
    }
    return result;
}

ElectrosprayHydroDiagnostics::ElectrosprayHydroDiagnostics(const std::string& prefix,
    bool charged, bool moving, bool ehd, bool capillary)
{
    diagnostics.emplace(prefix + (charged ? "_conducting_hydro_diagnostics.dat" : moving ? "_moving_dielectric_diagnostics.dat" : ehd ? "_frozen_ehd_diagnostics.dat"
        : capillary ? "_frozen_capillary_diagnostics.dat" : "_frozen_hydro_diagnostics.dat"), "Accepted hydrodynamic steps",
        std::vector<std::string>{"step","time_code","time_s","momentum_superlu","pressure_superlu","coupled_superlu",
            "momentum_residual","projection_residual","coupled_residual","full_momentum_residual","max_divergence_code",
            "max_speed_mps","supply_outward_Q_m3ps","ambient_outward_Q_m3ps","outlet_outward_Q_m3ps","net_boundary_Q_m3ps",
            "regular_momentum_rows","boundary_momentum_rows"});
}

void ElectrosprayHydroDiagnostics::initializeMotion(const std::string& prefix)
{
    motionDiagnostics.emplace(prefix + "_interface_motion_diagnostics.dat", "Accepted interface motion",
        std::vector<std::string>{"step","time_s","changed_phi_nodes","max_phi_increment_m","max_trace_ratio","liquid_volume_m3",
            "liquid_outward_Q_m3ps","integrated_liquid_input_m3","volume_defect_m3","changed_material_cells","electric_superlu","electric_residual"});
    phaseFluxDiagnostics.emplace(prefix + "_phase_flux_diagnostics.dat", "Accepted liquid phase flux balance",
        std::vector<std::string>{"step","time_s","active_faces","partitioned_half_faces","boundary_liquid_Q_m3ps",
            "summed_cell_liquid_Q_m3ps","assembly_closure_m3ps","absolute_face_liquid_Q_m3ps","pressure_flux_difference_m3ps"});
}

void ElectrosprayHydroDiagnostics::writeStep(int step, double time, const HydroSolveReport& report,
    const ElectrosprayBoundaryFlowRatesSI& boundaryFlow, const ElectrosprayHydroScales& scale)
{
    diagnostics->append(TecplotWriter::row(step, time, time * scale.time, report.momentumUsedSuperLU,
        report.pressureUsedSuperLU, report.coupledUsedSuperLU,
        report.momentumRelativeResidual, report.pressureRelativeResidual,
        report.coupledRelativeResidual, report.momentumBalanceRelativeResidual,
        report.divergenceAfterProjection, report.maxSpeed * scale.velocity,
        boundaryFlow.supplyOutward, boundaryFlow.ambientOutward, boundaryFlow.outletOutward, boundaryFlow.netOutward,
        report.fittedRegularMomentumRows, report.fittedBoundaryMomentumRows));
}

void ElectrosprayHydroDiagnostics::writeAcceptedMotion(int step, double time,
    const InterfaceMotionReport& report, const InterfacePhaseFlux::FluxBudget& phaseBudget,
    const ElectrosprayMotionStepValues& values, const ElectrosprayHydroScales& scale)
{
    if (!motionDiagnostics || !phaseFluxDiagnostics) throw std::logic_error("motion histories were not initialized");
    motionDiagnostics->append(TecplotWriter::row(step, time * scale.time, report.changedValues,
        report.maxValueIncrement * scale.length, report.maxTraceRatio,
        values.volumeSI, values.liquidFluxSI, values.expectedVolumeChangeSI,
        values.volumeSI-values.initialVolumeSI-values.expectedVolumeChangeSI, values.changedMaterials,
        values.electricSuperLU, values.electricResidual));
    phaseFluxDiagnostics->append(TecplotWriter::row(step, time*scale.time, phaseBudget.activeFaces,
        phaseBudget.partitionedHalfFaces, values.ledgerBoundarySI,
        phaseBudget.summedCellOutward*values.fluxFactorSI, phaseBudget.assemblyClosure*values.fluxFactorSI,
        phaseBudget.absoluteFaceFlux*values.fluxFactorSI, values.pressureDifferenceSI));
}
