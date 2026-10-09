#include "Modules/Electricity/IO/SurfaceChargeDiagnostics.h"
#include "Modules/Electricity/DielectricElectrostatics.h"
#include "Modules/Electricity/DielectricEvolution.h"
#include "Core/IO/OutputLayout.h"

#include <cmath>
#include <filesystem>
#include <iomanip>
#include <locale>
#include <sstream>
#include <stdexcept>

namespace
{
std::ofstream newOutput(const std::filesystem::path& path)
{
    if (std::filesystem::exists(path))
        throw std::runtime_error("refusing to overwrite surface charge diagnostics: " + path.string());
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot create surface charge diagnostics: " + path.string());
    out.imbue(std::locale::classic());
    out << std::scientific << std::setprecision(17);
    return out;
}
void flush(std::ofstream& out, const char* description)
{
    out.flush();
    if (!out) throw std::runtime_error(description);
}
}

void SurfaceChargeDiagnostics::initialize()
{
    if (summary_.is_open()) throw std::runtime_error("surface charge diagnostics already initialized");
    summary_ = newOutput(OutputLayout::series() / "surface_charge_diagnostics.csv");
    summary_ << "step,force_time,charge_time,dt,total_charge_before,total_charge_after,max_abs_density,"
        "net_boundary_current,unassigned_bulk_current,absolute_bulk_current,charge_balance_residual,"
        "current_closure_residual,relative_current_closure_residual,integrated_normal_current,"
        "transport_conservation_residual,relative_transport_conservation_residual,"
        "lower_boundary_current,upper_boundary_current,total_cell_current,absolute_face_current,"
        "max_displacement_split_residual,panel_count,max_advective_cfl,absolute_advective_charge_flux,"
        "remapped,previous_panel_count,remap_conservation_residual,max_material_displacement,max_projection_distance\n";
    flush(summary_, "failed to write surface charge diagnostics header");
}

void SurfaceChargeDiagnostics::writeStep(const Quantum::Electricity::DielectricChargeSnapshot& charge,
    const DielectricConductionReport& current,
    const Quantum::Electricity::SurfaceChargeFieldStep& update,
    int step, double forceTime, double dt, bool writeFields)
{
    if (!summary_.is_open() || !summary_)
        throw std::runtime_error("surface charge diagnostics is not initialized or has failed");
    if (step < 1 || !std::isfinite(forceTime) || !std::isfinite(dt) || dt <= 0.0
        || !std::isfinite(forceTime + dt))
        throw std::runtime_error("invalid surface charge diagnostic time label");
    if (charge.completedStep != step || charge.time != step * dt)
        throw std::runtime_error("charge snapshot does not match the completed diagnostic time layer");
    if (!update.transport.success)
        throw std::runtime_error("surface charge diagnostics requires a completed charge update");
    const auto& panels = charge.geometry->panels();
    const auto& charges = charge.integratedCharge;
    const auto& measures = charge.measure;
    if (panels.size() != charges.size() || panels.size() != measures.size())
        throw std::runtime_error("surface charge diagnostic state/panel size mismatch");

    const double chargeTime = forceTime + dt;
    const auto& transport = update.transport;
    // Net boundary current is OUTWARD; cell current and the panel source are
    // positive for accumulation. Do not silently transfer bulk current onto
    // the interface or subtract it from this independently reported balance.
    const double chargeBalance = (transport.totalChargeAfter - transport.totalChargeBefore)
        + dt * current.netBoundaryCurrent;
    double absoluteAdvectiveFlux = 0.0;
    for (double value : transport.advectiveChargeFlux) absoluteAdvectiveFlux += std::abs(value);
    if (!std::isfinite(absoluteAdvectiveFlux)) throw std::runtime_error("advective charge diagnostic overflow");

    if (writeFields)
    {
        std::ostringstream name;
        name << "surface_charge_step" << std::setw(4) << std::setfill('0') << step << ".csv";
        auto out = newOutput(OutputLayout::fields() / name.str());
        out << "step,force_time,charge_time,dt,panel_index,cell_index,x,y,nx,ny,curvature,"
            "measure,integrated_charge,mean_surface_charge_density";
        if (charge.geometry->axisymmetric()) out << ",meridional_arclength,axisymmetric";
        if (!charge.geometry->closed()) out << ",begin_r,begin_z,end_r,end_z";
        out << '\n';
        for (std::size_t i = 0; i < panels.size(); ++i)
        {
            const auto& frame = panels[i].frame;
            // Current retained geometry and charge are at n+1. Q/measure is a
            // panel mean, not a new local trace or a pointwise reconstruction.
            out << step << ',' << forceTime << ',' << chargeTime << ',' << dt << ',' << i
                << ',' << panels[i].cellIndex << ',' << frame.point.x << ',' << frame.point.y
                << ',' << frame.normal.x << ',' << frame.normal.y << ',' << frame.curvature
                << ',' << measures[i] << ',' << charges[i] << ',' << charges[i] / measures[i];
            if (charge.geometry->axisymmetric()) out << ',' << panels[i].measure << ",1";
            if (!charge.geometry->closed())
                out << ',' << panels[i].begin.x << ',' << panels[i].begin.y
                    << ',' << panels[i].end.x << ',' << panels[i].end.y;
            out << '\n';
        }
        flush(out, "failed to write surface charge panel diagnostics");
        out.close();
        if (!out) throw std::runtime_error("failed to close surface charge panel diagnostics");
    }

    summary_ << step << ',' << forceTime << ',' << chargeTime << ',' << dt
        << ',' << transport.totalChargeBefore << ',' << transport.totalChargeAfter << ',' << transport.maxAbsDensity
        << ',' << current.netBoundaryCurrent << ',' << update.unassignedBulkCurrent << ',' << update.absoluteBulkCurrent
        << ',' << chargeBalance << ',' << current.closureResidual << ',' << current.relativeClosure
        << ',' << transport.integratedNormalCurrent << ',' << transport.conservationResidual
        << ',' << transport.relativeConservationResidual << ',' << current.lowerBoundaryCurrent
        << ',' << current.upperBoundaryCurrent << ',' << current.totalCellCurrent << ',' << current.absoluteFaceCurrent
        << ',' << current.maxDisplacementSplitResidual << ',' << panels.size()
        << ',' << transport.maxAdvectiveCfl << ',' << absoluteAdvectiveFlux
        << ',' << update.remapped << ',' << update.previousPanelCount << ',' << update.remapConservationResidual
        << ',' << update.maxMaterialDisplacement << ',' << update.maxProjectionDistance << '\n';
    flush(summary_, "failed to write surface charge step diagnostics");
}
