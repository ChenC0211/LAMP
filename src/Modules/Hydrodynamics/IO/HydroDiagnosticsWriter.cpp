#include "Modules/Hydrodynamics/IO/HydroDiagnosticsWriter.h"
#include "Modules/Hydrodynamics/Hydrodynamics.h"
#include "Modules/Hydrodynamics/HydrodynamicsConfig.h"
#include "Modules/Hydrodynamics/HydroBoundaryVelocity.h"
#include "Core/Field/FieldChecks.h"
#include "Core/IO/OutputLayout.h"
#include "Core/IO/TecplotMeshOutput.h"
#include "Core/Geometry/ImmersedGeometry.h"
#include "Core/Mesh/QuadtreeMesh.h"

#include <filesystem>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
std::string csvText(const std::string& text)
{
    std::string result = "\"";
    for (char c : text) { if (c == '"') result += '"'; result += c; }
    return result + "\"";
}
}
namespace HydroOutput
{
using FieldChecks::checkFinite;

DiagnosticsWriter::DiagnosticsWriter(const HydrodynamicsConfig& config, bool scalarBoundary)
    : config_(config), scalarBoundary_(scalarBoundary)
{
    const bool openBoundary = !config.boundaryPatches.empty();
    const auto diagnosticsPath = OutputLayout::series() / "hydro_diagnostics.csv";
    diagnostics.open(diagnosticsPath);
    if (!diagnostics) throw std::runtime_error("cannot open " + diagnosticsPath.string());
    diagnostics << "step,time,dt,success,momentum_relative_residual,pressure_relative_residual,"
                   "divergence_before_projection,divergence_after_projection,max_speed,"
                   "axis_max_radial_velocity,stage,reason,interface_max_change,negative_volume,"
                   "rho_min,rho_max,mu_min,mu_max,cut_cells,levelset_min,levelset_max,"
                   "surface_tension,interface_node_count,interface_surface_measure,max_abs_curvature,"
                   "force_time,interface_force_x,interface_force_y,spread_force_x,spread_force_y,"
                   "max_interface_force_density,interface_normal_traction,interface_tangential_traction,"
                   "active_capillary_face_force_x,active_capillary_face_force_y,max_active_capillary_face_force_density,"
                   "max_force_pressure_imbalance,l2_force_pressure_imbalance,"
                   "coupled_correction_used,coupled_relative_residual,momentum_balance_relative_residual,"
                   "momentum_balance_max_force_density,momentum_factor_reused,pressure_factor_reused,"
                   "coupled_factor_reused,coupled_superlu,legacy_reference_available";
    if (openBoundary) diagnostics << ",inflow_volume_flux,outflow_volume_flux,net_boundary_volume_flux,net_boundary_mass_flux";
    if (scalarBoundary_) diagnostics << ",geometric_negative_volume,initial_geometric_negative_volume,integrated_negative_phase_inflow,geometric_volume_budget_defect";
    diagnostics << '\n' << std::setprecision(17);
    if(config.timeScheme=="cchMidpoint"){
        energy.open(OutputLayout::series()/"hydro_energy.csv");
        if(!energy)throw std::runtime_error("cannot open Hydro midpoint energy history");
        energy<<"step,time,dt,pressure_time,mass,initial_mass,mass_change,kinetic_energy,initial_kinetic_energy,"
            "cumulative_viscous_dissipation,cumulative_external_work,energy_balance,pressure_power,convection_power,"
            "viscous_power,external_power,dual_mass_relative_closure,midpoint_iterations,midpoint_velocity_defect\n"<<std::setprecision(17);
    }
}
void DiagnosticsWriter::appendMidpoint(int step,double time,double dt,double pressureTime,const HydroSolveReport& report,
    const HydroFlowDiagnostics::NativeConservationDiagnostics& values,int iterations,double velocityDefect)
{
    if(!energy.is_open() || step!=energyStep_+1 || (step>0 && (!report.success || !report.midpointNativeMomentumUsed)))
        throw std::logic_error("midpoint energy output requires its next accepted physical state");
    if(step==0){initialMass_=values.mass;initialEnergy_=values.kineticEnergy;}
    else {dissipation_-=dt*report.physicalViscousPower;externalWork_+=dt*report.physicalExternalPower;}
    energy<<step<<','<<time<<','<<dt<<','<<pressureTime<<','<<values.mass<<','<<initialMass_<<','
        <<values.mass-initialMass_<<','<<values.kineticEnergy<<','<<initialEnergy_<<','<<dissipation_<<','<<externalWork_<<','
        <<values.kineticEnergy+dissipation_-externalWork_-initialEnergy_<<','
        <<(step?report.physicalPressurePower:0)<<','<<(step?report.physicalConvectionPower:0)<<','
        <<(step?report.physicalViscousPower:0)<<','<<(step?report.physicalExternalPower:0)<<','
        <<(step?report.conservativeMassRelativeClosure:0)<<','<<iterations<<','<<velocityDefect<<'\n';
    energy.flush();if(!energy)throw std::runtime_error("failed to write Hydro midpoint energy history");
    energyStep_=step;
}

void DiagnosticsWriter::appendReport(int step, double time, double dt, const HydroSolveReport& report, double axisRadialSpeed,
    const HydroFlowDiagnostics::GeometryMaterialDiagnostics& geometryDiagnostic,
    const HydroInterfaceForce::Report& forceDiagnostic,
    const HydroFlowDiagnostics::BoundaryFlowDiagnostics& flowDiagnostic,
    double geometricVolume, double initialGeometricVolume, double integratedNegativeInflow)
{
    const auto& config = config_;
    const bool openBoundary = !config.boundaryPatches.empty();
    diagnostics << step << ',' << time << ',' << dt << ',' << report.success
                << ',' << report.momentumRelativeResidual << ',' << report.pressureRelativeResidual
                << ',' << report.divergenceBeforeProjection << ',' << report.divergenceAfterProjection
                << ',' << report.maxSpeed << ',' << axisRadialSpeed
                << ',' << csvText(report.stage) << ',' << csvText(report.reason)
                << ',' << geometryDiagnostic.interfaceMaxChange << ',' << geometryDiagnostic.negativeVolume
                << ',' << geometryDiagnostic.rhoMin << ',' << geometryDiagnostic.rhoMax
                << ',' << geometryDiagnostic.muMin << ',' << geometryDiagnostic.muMax
                << ',' << geometryDiagnostic.cutCells << ',' << geometryDiagnostic.levelSetMin
                << ',' << geometryDiagnostic.levelSetMax
                << ',' << config.surfaceTension;
    const auto legacyValue = [&](const auto& value) {
        diagnostics << ',';
        if (forceDiagnostic.legacyReferenceAvailable) diagnostics << value;
    };
    legacyValue(forceDiagnostic.nodeCount);
    legacyValue(forceDiagnostic.surfaceMeasure);
    legacyValue(forceDiagnostic.maxAbsCurvature);
    diagnostics << ',' << forceDiagnostic.forceTime;
    legacyValue(forceDiagnostic.forceX);
    legacyValue(forceDiagnostic.forceY);
    legacyValue(forceDiagnostic.spreadForceX);
    legacyValue(forceDiagnostic.spreadForceY);
    legacyValue(forceDiagnostic.maxForceDensity);
    diagnostics << ',' << config.interfaceNormalTraction << ',' << config.interfaceTangentialTraction
                << ',' << forceDiagnostic.activeCapillaryFaceX << ',' << forceDiagnostic.activeCapillaryFaceY
                << ',' << forceDiagnostic.maxActiveCapillaryFaceDensity
                << ',' << forceDiagnostic.maxForcePressureImbalance << ',' << forceDiagnostic.l2ForcePressureImbalance
                << ',' << report.coupledCorrectionUsed << ',' << report.coupledRelativeResidual
                << ',' << report.momentumBalanceRelativeResidual << ',' << report.momentumBalanceMaxForceDensity
                << ',' << report.momentumFactorReused << ',' << report.pressureFactorReused
                << ',' << report.coupledFactorReused << ',' << report.coupledUsedSuperLU
                << ',' << forceDiagnostic.legacyReferenceAvailable;
    if (openBoundary)
    {
        if (report.success)
            diagnostics << ',' << flowDiagnostic.inflow << ',' << flowDiagnostic.outflow
                        << ',' << flowDiagnostic.net << ',' << flowDiagnostic.massNet;
        else diagnostics << ",,,,";
    }
    if (scalarBoundary_)
    {
        if (report.success)
            diagnostics << ',' << geometricVolume << ',' << initialGeometricVolume
                        << ',' << integratedNegativeInflow << ','
                        << geometricVolume-initialGeometricVolume-integratedNegativeInflow;
        else diagnostics << ",,,,";
    }
    diagnostics << '\n';
    diagnostics.flush();
    if (!diagnostics) throw std::runtime_error("failed to write Hydro diagnostics");
}

void writeFields(const QuadtreeMesh& mesh, const Hydrodynamics& hydro,
                 const ImmersedGeometry<QuadtreeMesh>* geometry, int interfaceIndex,
                 const HydrodynamicsConfig& config, int step, double time,
                 const std::string& outputName, bool radialAxisSymmetry, bool scalarBoundary,double pressureTime)
{
    std::ostringstream name;
    name << outputName << "_step" << std::setw(4) << std::setfill('0') << step;
    std::vector<std::string> names{"rho", "mu", "p", "ux", "uy", "divU"};
    const auto [vx, vy] = hydro.vertexVelocity(mesh, HydroVelocityTimeLevel::Current,
        !config.boundaryPatches.empty(), radialAxisSymmetry);
    if (scalarBoundary)
    {
        const auto velocity = hydro.boundaryAdvectionVelocity(mesh, HydroVelocityTimeLevel::Current,
            config.boundaryPatches, radialAxisSymmetry);
        for (int i = 0; i < mesh.vertIndexNum(); ++i)
        {
            const auto value = velocity(mesh.verts(mesh.vertIndexToLabel[i]).coord());
            vx->value[i] = value.x;
            vy->value[i] = value.y;
        }
    }
    const auto state = hydro.state();
    std::vector<VertField<double>> fields{state.vertexDensity.copy(), state.vertexViscosity.copy(),
        *state.physicalPressure.toVertField(mesh), *vx, *vy, *state.divergence.toVertField(mesh)};
    if (geometry)
    {
        names.push_back("levelSet");
        fields.push_back(geometry->levelSets[interfaceIndex]);
    }
    for (std::size_t i = 0; i < fields.size(); ++i)
        checkFinite(fields[i], "output " + names[i], mesh.vertIndexNum());
    const auto file = OutputLayout::fields() / (name.str() + ".dat");
    auto zone = TecplotMeshOutput::nodalZone(mesh, names, fields);
    zone.solutionTime = time;
    zone.strandId = 1;
    zone.auxiliary = {{"CoordinateUnits", "code_length"}, {"FieldUnits", "code_units"},
        {"TimeUnits", "code_time"}, {"TimeRole", "completed_step_time"},
        {"FieldRole", "hydrodynamic_nodal_diagnostics"},
        {"CoordinateSystem", config.coordinateSystem}};
    if(config.timeScheme=="cchMidpoint"){
        if(pressureTime<0 || !std::isfinite(pressureTime) || pressureTime>time)
            throw std::invalid_argument("midpoint field output requires its actual pressure time");
        std::ostringstream label;label<<std::setprecision(17)<<pressureTime;
        zone.auxiliary.emplace("PressureTime",label.str());
        zone.auxiliary.emplace("PressureTimeRole",step?"step_midpoint":"initial_condition");
        zone.auxiliary.emplace("VelocityTimeRole","completed_step_time");
        zone.auxiliary.emplace("TimeScheme","cchMidpoint");
    }
    TecplotWriter::write(file, zone);
    if (!std::filesystem::exists(file) || std::filesystem::file_size(file) == 0)
        throw std::runtime_error("failed to write Hydro field output: " + file.string());
}
}
