#include "Modules/Plasma/IO/PicStepSeries.h"
#include "Modules/Plasma/Boundary/LegacyGridMotion.h"
#include "Modules/Plasma/ExplicitPIC/PicConfig.h"
#include "Modules/Plasma/ImplicitPIC/ImplicitPicOverrides.h"
#include "Modules/Plasma/Particles/ParticleManager.h"
#include "Modules/Plasma/Particles/PicSpeciesState.h"
#include "Modules/Plasma/Diagnostics/TimeSeriesIO.h"
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <charconv>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include "Core/Units/UnitsContext.h"

namespace
{
void openExplicitReport(std::ofstream& stream, const std::filesystem::path& path,
    int resumeStep, const char* header, const char* role)
{
    const bool hasExisting = TimeSeriesIO::openPreparedStepCsv(stream, path, resumeStep);
    if (stream)
    {
        if (!hasExisting) stream << header;
        stream << std::setprecision(16);
    }
    else
        std::cerr << "[ExplicitPIC] WARNING: failed to create " << role << ": " << path.string() << std::endl;
}

std::string formatPicEta(double seconds)
{
    if (!std::isfinite(seconds) || seconds < 0.0)
        return "N/A";
    long long sec = static_cast<long long>(std::llround(seconds));
    const long long h = sec / 3600;
    sec %= 3600;
    const long long m = sec / 60;
    const long long s = sec % 60;
    std::ostringstream oss;
    oss << std::setfill('0') << std::setw(2) << h << ':'
        << std::setw(2) << m << ':'
        << std::setw(2) << s;
    return oss.str();
}
constexpr const char* amrAdaptationHeader = "step,cells_before,cells_after,use_fast_remap,"
                         "mesh_adapt_s,movingeb_onMeshChanged_s,geometry_refresh_s,geometry_reinit_skipped,"
                         "geom_initialize_levelsets_s,geom_pde_reinit_s,geom_refresh_derived_s,geom_check_cutinfo_s,"
                         "geom_rebuild_composite_s,geom_mark_phases_s,face_graph_s,embedded_metrics_s,"
                         "field_buffer_reset_s,relabel_particles_s,field_rebuild_s,coord_refresh_s,total_adapt_s\n";
constexpr const char* amrLevelSetsHeader = "step,level_set_index,moving,active_moving_levelset,selective_reinit_target,"
                                "hybrid_sensitive_inactive_target,"
                                "reinitialized,reinit_iterations,reinit_fine_residual,reinit_coarse_residual,"
                                "min_value,max_value,near_zero_count,nonfinite_count,narrow_band_count\n";
constexpr const char* explicitPerformanceHeader = "step,time_code,buildWeights_s,deposit_s,rhs_s,poisson_assemble_s,poisson_factor_s,poisson_solve_s,poisson_total_s,evalE_s,amr_indicator_s,amr_adapt_s,pm_s,gatherPush_s,output_s,total_s,particles_total,poisson_N,poisson_nnz,poisson_usedSuperLU\n";
constexpr const char* implicitPerformanceHeader = "step,time_s,np,total_s,mcc_pre_s,buildWeightsN_s,predictorSegment_s,orbitAvgWeights_s,polPoisson_total_s,pm_s,amr_s,amr_indicator_s,amr_adapt_s,amr_faceGraph_s,amr_cutMetrics_s,amr_reinitLevelset_s,amr_geom_initializeLevelSets_s,amr_geom_pdeReinitLevelset_s,amr_geom_refreshDerived_s,amr_geom_checkCutInfo_s,amr_geom_rebuildComposite_s,amr_geom_markPhases_s,amr_fieldRebuild_s,output_s,is_output_step,is_amr_adapt_step,amr_skipped_by_gate,pol_N,pol_nnz,pol_solver_mode,predictor_segments_per_particle,predictor_coarsen_hits_per_particle,predictor_coarsen_parentClimb_per_hit,orbitAvg_samples_per_particle,orbitAvg_sampleCache_nnz_per_sample,orbitAvg_adapt_accept_per_particle,orbitAvg_adapt_refine_per_particle,orbitAvg_adapt_hitMaxDepth_per_particle,orbitAvg_adapt_meanAcceptedDepth,top1_name,top1_s,top2_name,top2_s,top3_name,top3_s\n";
constexpr const char* fallbackHeader = "step,time_s,buildWeightsN_cwFallbackCount,buildWeightsN_particleCount,buildWeightsN_cwFallbackFrac,orbitAvg_sampleCwFallbackCount,orbitAvg_sampleCount,orbitAvg_sampleCwFallbackFrac,cwFallbackCount_total\n";
constexpr const char* explicitOverviewHeader = "step,time_s,rG_L2,rG_Linf,rG_L2_per_active,rG_L2_rel_rhs,rG_Linf_rel_rhs,S_eb_rel_rhs,field_energy_J,kinetic_e_J,kinetic_ion_J,kinetic_neu_J,kinetic_total_J,total_energy_J,J_cond_Apm2,J_disp_Apm2,J_total_Apm2,macro_total,phys_total,charge_balance_rel\n";
constexpr const char* implicitOverviewHeader = "step,time_s,rG_L2,rG_Linf,field_energy_J,kinetic_e_J,kinetic_ion_J,kinetic_neu_J,kinetic_total_J,total_energy_J,J_cond_Apm2,J_disp_Apm2,J_total_Apm2,macro_total,phys_total,charge_balance_rel\n";
constexpr const char* processHeader = "step,time_s,drive_energy_J,injection_energy_J,reaction_energy_loss_J,domain_absorption_J,embedded_absorption_J,see_emission_energy_J,ion_impact_energy_J,displacement_current_fraction,dominant_source,dominant_source_fraction,dominant_boundary,dominant_boundary_absorbed_power_W,dominant_boundary_absorbed_fraction,dominant_reaction_channel,dominant_reaction_power_W,dominant_reaction_fraction\n";
constexpr const char* periodicHeader = "step,time_s,phiJumpX_V,phiJumpY_V,periodicX,periodicY\n";
constexpr const char* geometryHeader = "step,time_s,geometry_updated,active_cut_cell_count,embedded_surface_area_m2,moving_region_count,moving_levelset_count,moving_region_total_area_m2,moving_region_total_flux_1pm2s,moving_region_mean_energy_eV,moving_region_mean_cos,moving_region_mean_abs_speed_mps,moving_region_max_abs_speed_mps,geometry_gsa_mm,geometry_gad_mm,geometry_da_mm,field_energy_J,domain_absorption_J,embedded_absorption_J,reaction_energy_loss_J\n";
}

void PicStepSeries::openGeometryMetrics(const std::filesystem::path& seriesDir,
    std::filesystem::path path, int resumeStep, int every)
{
    geometryMetricsEvery_ = std::max(1, every);
    if (path.empty())
        path = seriesDir / "geometry_step.csv";
    else if (path.is_relative())
    {
        const auto parent = path.parent_path();
        if (parent.empty() || parent == ".") path = seriesDir / path;
        else path = std::filesystem::path(".") / path;
    }
    const auto parent = path.parent_path();
    if (!parent.empty()) std::filesystem::create_directories(parent);
    open(geometryMetrics_, path, resumeStep,
        "step,time_s,stage,stage1_alpha,stage2_alpha,gsa_mm,gad_mm,da_mm,"
        "accel_shift_mm,decel_shift_mm,accel_recession_mm,"
        "screen_downstream_x_mm,accel_upstream_x_mm,accel_downstream_local_x_mm,decel_upstream_x_mm\n",
        "geometry metrics");
}
void PicStepSeries::writeGeometryMetrics(int step, double timeS, const LegacyGridMotionState& state, bool force)
{
    if (!geometryMetrics_.is_open() || !geometryMetrics_.good()) return;
    if (!force && !state.changed && step % geometryMetricsEvery_ != 0) return;
    auto toMm = [](double value_m) { return value_m * 1.0e3; };
    geometryMetrics_ << step << ',' << timeS << ',' << state.stageName << ','
        << state.stage1Alpha << ',' << state.stage2Alpha << ','
        << toMm(state.gsa_m) << ',' << toMm(state.gad_m) << ',' << toMm(state.da_m) << ','
        << toMm(state.accelShift_m) << ',' << toMm(state.decelShift_m) << ',' << toMm(state.accelRecession_m) << ','
        << toMm(state.screenDownstreamX_m) << ',' << toMm(state.accelUpstreamX_m) << ','
        << toMm(state.accelDownstreamLocalX_m) << ',' << toMm(state.decelUpstreamX_m) << "\n";
    geometryMetrics_.flush();
}
void PicStepSeries::openGeometryTrace(const std::filesystem::path& path, int selectedStep, bool allSteps)
{
    geometryTraceStep_ = selectedStep;
    geometryTraceAll_ = allSteps;
    if (!path.empty())
    {
        const auto parent = path.parent_path();
        if (!parent.empty()) std::filesystem::create_directories(parent);
        geometryTrace_.open(path, std::ios::out | std::ios::app);
    }
}
void PicStepSeries::writeGeometryTrace(int step, const std::string& message)
{
    if (!geometryTraceAll_ && !(geometryTraceStep_ >= 0 && step >= geometryTraceStep_ - 2 && step <= geometryTraceStep_ + 2)) return;
    const std::string line = "[moving-trace] step=" + std::to_string(step) + ' ' + message;
    std::cout << line << std::endl;
    if (geometryTrace_.is_open())
    {
        geometryTrace_ << line << '\n';
        geometryTrace_.flush();
    }
}

void PicStepSeries::openSecondaryEmission(const std::filesystem::path& path)
{
    secondaryEmission_.open(path, std::ios::out | std::ios::trunc);
    if (secondaryEmission_)
        secondaryEmission_ << "# see_report (ion impacts + SEE emission; 3D hemisphere distribution)\n"
            << "# note: ionImpact_* counts sum over configured SEE sides\n"
            << "# columns: step,time_s,"
            << "ionImpact_macro_cathode,ionImpact_phys_cathode,ionImpactEnergy_J_cathode,"
            << "emitted_macro,emitted_phys,emittedEnergy_J,emittedEnergy_mean_eV,emittedEnergy_std_eV\n";
    else
        std::cerr << "[ExplicitPIC] WARNING: failed to create see report: " << path.string() << std::endl;
}
void PicStepSeries::writeSecondaryEmission(int step, double timeS, const PicSecondaryEmission::Stats& row,
    const PicEmissionEnergy& energy, bool flush)
{
    secondaryEmission_ << step << ',' << std::setprecision(16) << timeS << ','
        << row.ionImpactMacroCathode << ',' << row.ionImpactPhysCathode << ',' << row.ionImpactEnergyJCathode << ','
        << row.emittedMacro << ',' << row.emittedPhys << ',' << row.emittedEnergyJ << ','
        << energy.meanEV << ',' << energy.stdEV << "\n";
    if (flush) secondaryEmission_.flush();
}
void PicStepSeries::openImpact(const std::filesystem::path& path, int resumeStep)
{
    openExplicitReport(impact_, path, resumeStep,
        "step,time_s,out_macro_e_L,out_macro_e_R,out_macro_e_B,out_macro_e_T,out_macro_e_unknown,"
        "out_macro_i_L,out_macro_i_R,out_macro_i_B,out_macro_i_T,out_macro_i_unknown\n", "impact report");
}
void PicStepSeries::writeImpact(int step, double timeS, const PicBoundaryLosses::Totals& losses)
{
    impact_ << step << ',' << timeS;
    for (const auto count : losses.domainMacroE) impact_ << ',' << count;
    impact_ << ',' << losses.unknownE;
    for (const auto count : losses.domainMacroI) impact_ << ',' << count;
    impact_ << ',' << losses.unknownI << "\n";
    impact_.flush();
}
void PicStepSeries::openRfPhase(const std::filesystem::path& path, int resumeStep, bool writeSI)
{
    openExplicitReport(rfPhase_, path, resumeStep,
        writeSI ? "step,time_s,target,type,freq_hz,phase_rad,phase_deg,phi_SI_V\n"
                : "step,time,target,type,freq_hz,phase_rad,phase_deg,phi\n", "RF phase report");
}
void PicStepSeries::writeRfPhase(int step, double time, const std::vector<PicRfPhase>& rows)
{
    for (const auto& row : rows)
        rfPhase_ << step << ',' << time << ',' << row.target << ',' << row.type << ',' << row.frequencyHz
            << ',' << row.radians << ',' << row.degrees << ',' << row.potential << "\n";
}

void PicStepSeries::openAmrIndicators(const std::filesystem::path& path, int resumeStep)
{
    const bool hasExisting = TimeSeriesIO::openPreparedStepCsv(amrIndicators_, path, resumeStep);
    if (amrIndicators_)
    {
        if (!hasExisting)
            amrIndicators_ << "step,time_s,lambdaD_min_m,h_min_m,ratio_max,ratio_mean,refine_cells,coarsen_cells,pending,adapted\n";
        amrIndicators_ << std::setprecision(16);
    }
    else
    {
        std::cerr << "[ExplicitPIC] WARNING: failed to create AMR report: " << path.string() << std::endl;
    }
}

void PicStepSeries::writeAmrIndicators(int stepForPrint, double t_si,
    const ExplicitPicRefinement::Report& amrStats, bool pending, bool adaptedThisStep)
{
    if (!amrIndicatorsReady()) return;
    amrIndicators_ << stepForPrint << ',' << std::setprecision(16) << t_si << ','
           << amrStats.lambdaD_min_m << ','
           << amrStats.h_min_m << ','
           << amrStats.ratio_max << ','
           << amrStats.ratio_mean << ','
           << amrStats.refine_cells << ','
           << amrStats.coarsen_cells << ','
           << (pending ? 1 : 0) << ','
           << (adaptedThisStep ? 1 : 0)
           << "\n";
    amrIndicators_.flush();
}

void PicStepSeries::writeAmrDebyeWarning(std::ostream& stream, const ExplicitPicRefinement::Report& report)
{
    const int debyeLowNeCount = report.debyeLowNeCount, debyeBadCount = report.debyeBadCount;
    const double neFloor = report.neFloor;
    if (debyeLowNeCount > 0 || debyeBadCount > 0)
    {
        stream << "[AMR][Debye] lowNe=" << debyeLowNeCount
                  << " badLambdaD=" << debyeBadCount
                  << " (neFloor=" << std::scientific << std::setprecision(3) << neFloor << ")"
                  << std::defaultfloat << std::setprecision(6)
                  << std::endl;
    }
}

void PicStepSeries::writeAmrConsole(std::ostream& stream, int stepForPrint, int cellsAfter,
    double minimumLeafWidth, const ExplicitPicRefinement::AdaptationReport& adaptation)
{
    const auto& levelsBefore = adaptation.levelsBefore;
    const auto& levelsAfter = adaptation.levelsAfter;
    std::size_t nLevel = std::max(levelsBefore.size(), levelsAfter.size());
    std::uint64_t added = 0;
    std::uint64_t removed = 0;
    for (std::size_t i = 0; i < nLevel; ++i)
    {
        const int b = (i < levelsBefore.size()) ? levelsBefore[i] : 0;
        const int a = (i < levelsAfter.size()) ? levelsAfter[i] : 0;
        const int d = a - b;
        if (d > 0)
            added += static_cast<std::uint64_t>(d);
        else if (d < 0)
            removed += static_cast<std::uint64_t>(-d);
    }

    std::ostringstream levels;
    levels << "{";
    bool first = true;
    for (std::size_t i = 0; i < levelsAfter.size(); ++i)
    {
        if (!first) levels << ',';
        first = false;
        levels << i << ":" << levelsAfter[i];
    }
    levels << "}";

    stream << "[AMR] step=" << stepForPrint
              << " cells=" << cellsAfter
              << " add=" << added
              << " del=" << removed
              << " hpMinLeaf=" << std::scientific << std::setprecision(6) << minimumLeafWidth
              << " levels=" << levels.str()
              << std::defaultfloat << std::setprecision(6)
              << std::endl;
}

void PicStepSeries::openAmrAdaptation(const std::filesystem::path& path, int resumeStep)
{ open(amrAdaptation_, path, resumeStep, amrAdaptationHeader, "amr adapt"); }

void PicStepSeries::openAmrLevelSets(const std::filesystem::path& path, int resumeStep)
{ open(amrLevelSets_, path, resumeStep, amrLevelSetsHeader, "amr levelset diag"); }

void PicStepSeries::writeAmrAdaptation(int stepForPrint,
    const ImplicitPicRefinement::AdaptationReport& adaptation, const PicAmrGeometry::Report& geometryRefresh)
{
    if (!amrAdaptationReady()) return;
    amrAdaptation_ << stepForPrint << ','
                << adaptation.cellsBefore << ','
                << adaptation.cellsAfter << ','
                << (adaptation.useFastRemap ? 1 : 0) << ','
                << adaptation.meshAdaptS << ','
                << adaptation.movingEbMeshChangedS << ','
                << geometryRefresh.geometryRefreshS << ','
                << (geometryRefresh.geometryReinitSkipped ? 1 : 0) << ','
                << geometryRefresh.geomInitializeLevelSetsS << ','
                << geometryRefresh.geomPdeReinitS << ','
                << geometryRefresh.geomRefreshDerivedS << ','
                << geometryRefresh.geomCheckCutInfoS << ','
                << geometryRefresh.geomRebuildCompositeS << ','
                << geometryRefresh.geomMarkPhasesS << ','
                << adaptation.faceGraphS << ','
                << adaptation.embeddedMetricsS << ','
                << adaptation.fieldResetS << ','
                << adaptation.relabelParticlesS << ','
                << adaptation.fieldRebuildS << ','
                << adaptation.coordRefreshS << ','
                << adaptation.seconds
                << "\n";
    amrAdaptation_.flush();
}

void PicStepSeries::writeAmrLevelSets(int stepForPrint, const PicAmrGeometry::Geometry* iGeometry,
    const PicAmrGeometry::Report& geometryRefresh)
{
    if (!amrLevelSetsReady() || !iGeometry) return;
    std::vector<unsigned char> activeMask(iGeometry->levelSets.size(), static_cast<unsigned char>(0));
    std::vector<unsigned char> selectiveMask(iGeometry->levelSets.size(), static_cast<unsigned char>(0));
    std::vector<unsigned char> hybridInactiveMask(iGeometry->levelSets.size(), static_cast<unsigned char>(0));
    for (const int idx : geometryRefresh.amrRawActiveMovingLevelSets)
    {
        if (idx >= 0 && idx < static_cast<int>(activeMask.size()))
            activeMask[static_cast<std::size_t>(idx)] = 1u;
    }
    for (const int idx : geometryRefresh.amrSelectiveMovingLevelSets)
    {
        if (idx >= 0 && idx < static_cast<int>(selectiveMask.size()))
            selectiveMask[static_cast<std::size_t>(idx)] = 1u;
    }
    for (const int idx : geometryRefresh.amrHybridSensitiveInactiveLevelSets)
    {
        if (idx >= 0 && idx < static_cast<int>(hybridInactiveMask.size()))
            hybridInactiveMask[static_cast<std::size_t>(idx)] = 1u;
    }
    const auto& lsDiags = iGeometry->levelSetUpdateDiagnostics();
    for (int lsI = 0; lsI < static_cast<int>(iGeometry->levelSets.size()) &&
                      lsI < static_cast<int>(lsDiags.size()); ++lsI)
    {
        const auto& ls = iGeometry->levelSets[lsI];
        const auto& diag = lsDiags[lsI];
        amrLevelSets_ << stepForPrint << ','
                           << (lsI + 1) << ','
                           << (ls.moving() ? 1 : 0) << ','
                           << (activeMask[static_cast<std::size_t>(lsI)] ? 1 : 0) << ','
                           << (selectiveMask[static_cast<std::size_t>(lsI)] ? 1 : 0) << ','
                           << (hybridInactiveMask[static_cast<std::size_t>(lsI)] ? 1 : 0) << ','
                           << (diag.reinitialized ? 1 : 0) << ','
                           << diag.reinitIterations << ','
                           << diag.reinitFineResidual << ','
                           << diag.reinitCoarseResidual << ','
                           << diag.minValue << ','
                           << diag.maxValue << ','
                           << diag.nearZeroCount << ','
                           << diag.nonFiniteCount << ','
                           << diag.narrowBandCount
                           << "\n";
    }
    amrLevelSets_.flush();
}

void PicStepSeries::writeAmrConsole(std::ostream& stream, int stepForPrint,
    const ImplicitPicRefinement::AdaptationReport& adaptation, const PicAmrGeometry::Report& geometryRefresh)
{
    stream << "[ImplicitPIC][Event][AMR] step=" << stepForPrint
              << " cells " << adaptation.cellsBefore << "->" << adaptation.cellsAfter
              << " adapt_s=" << std::fixed << std::setprecision(3) << adaptation.seconds
              << " fieldRebuild=" << (adaptation.useFastRemap ? "remap+evalE" : "deposit+solve")
              << " geom[pde=" << geometryRefresh.geomPdeReinitS
              << ",derived=" << geometryRefresh.geomRefreshDerivedS
              << ",cut=" << geometryRefresh.geomCheckCutInfoS
              << ",comp=" << geometryRefresh.geomRebuildCompositeS
              << ",mark=" << geometryRefresh.geomMarkPhasesS
              << "]"
              << std::defaultfloat << std::setprecision(6)
              << std::endl;
}

void PicStepSeries::writeExplicitDiagnostic(std::ostream& stream, int stepForPrint,
    double totalGridCharge, double totalMacroCharge,
    const ExplicitPIC::Diagnostics::PicDiagnostics::ExplicitSample& sample, std::string_view message)
{
    const auto& diag = sample.values;
    const double phiMin = sample.phiMin, phiMax = sample.phiMax;
    const double dq = totalGridCharge - totalMacroCharge;
    const double rel = (std::abs(totalMacroCharge) > 0.0) ? (dq / totalMacroCharge) : 0.0;
    stream << std::setprecision(16);
    stream << "[ExplicitPIC][diag] step=" << stepForPrint
              << " Qgrid=" << totalGridCharge
              << " Qmacro=" << totalMacroCharge
              << " dQ=" << dq
              << " rel=" << rel
              << " rG_L2=" << diag.rG_L2
              << " rG_Linf=" << diag.rG_Linf
              << " W=" << diag.energyW
              << " phi[min,max]=[" << phiMin << ", " << phiMax << "]";

    if (!message.empty())
        stream << " poisson: " << message;
    stream << std::endl;
}

void PicStepSeries::open(std::ofstream& stream, const std::filesystem::path& path,
    int resumeStep, const char* header, const char* role)
{
    const bool hasExisting = TimeSeriesIO::openPreparedStepCsv(stream, path, resumeStep);
    if (stream)
    {
        if (!hasExisting) stream << header;
        stream << std::setprecision(16) << std::fixed;
        if (mode_ == Mode::Implicit) stream.flush();
    }
    else if (role)
        std::cerr << '[' << (mode_ == Mode::Explicit ? "ExplicitPIC" : "ImplicitPIC")
                  << "] WARNING: failed to create " << role << " csv: " << path.string() << std::endl;
}

void PicStepSeries::openOverview(const std::filesystem::path& path, int resumeStep)
{
    open(overview_, path, resumeStep,
        mode_ == Mode::Explicit ? explicitOverviewHeader : implicitOverviewHeader, "overview primary");
}
void PicStepSeries::openProcess(const std::filesystem::path& path, int resumeStep)
{ open(process_, path, resumeStep, processHeader, "process summary"); }
void PicStepSeries::openGeometry(const std::filesystem::path& path, int resumeStep)
{ open(geometry_, path, resumeStep, geometryHeader, "geometry coupling"); }

void PicStepSeries::openPeriodic(const std::filesystem::path& path, int resumeStep)
{ open(periodic_, path, resumeStep, periodicHeader, "periodic diagnostics"); }

void PicStepSeries::openWaveform(const std::filesystem::path& path, int resumeStep, int every)
{
    waveformEvery_ = std::max(1, every);
    waveform_.open(path.string(), resumeStep);
    if (!waveform_.isOpen())
        std::cerr << '[' << (mode_ == Mode::Explicit ? "ExplicitPIC" : "ImplicitPIC")
                  << "] WARNING: failed to create current waveform csv: " << path.string() << std::endl;
}
void PicStepSeries::writeWaveform(int step, double timeS, double voltageV, double conduction, double displacement)
{ waveform_.write(step, timeS, voltageV, conduction, displacement); }

void PicStepSeries::writePeriodic(int step, double time_s, double jumpX_V, double jumpY_V,
    bool periodicX, bool periodicY, bool reportStep)
{
    periodic_ << step << ',' << std::setprecision(16) << time_s << ','
              << std::setprecision(16) << jumpX_V << ','
              << std::setprecision(16) << jumpY_V << ','
              << (periodicX ? 1 : 0) << ',' << (periodicY ? 1 : 0) << "\n";
    if (mode_ == Mode::Explicit || reportStep) periodic_.flush();
}

template<class MacroCount>
void PicStepSeries::writeOverviewValues(const OverviewValues<MacroCount>& row)
{
    overview_ << ',' << row.fieldEnergyJ << ',' << row.kineticElectronJ << ',' << row.kineticIonJ
              << ',' << row.kineticNeutralJ << ',' << row.kineticTotalJ << ',' << row.totalEnergyJ
              << ',' << row.conductionCurrent << ',' << row.displacementCurrent << ',' << row.totalCurrent
              << ',' << row.macroTotal << ',' << row.physicalTotal << ',' << row.chargeBalanceRelative
              << "\n";
}
void PicStepSeries::writeOverview(int step, double time_s, const ExplicitOverview& row)
{
    overview_ << step << ',' << time_s << ',' << row.residualL2 << ',' << row.residualLinf
              << ',' << row.residualL2PerActive << ',' << row.residualL2RelativeRhs
              << ',' << row.residualLinfRelativeRhs << ',' << row.embeddedRelativeRhs;
    writeOverviewValues(row.values);
    overview_.flush();
}
void PicStepSeries::writeOverview(int step, double time_s, const ImplicitOverview& row)
{
    overview_ << step << ',' << time_s << ',' << row.residualL2 << ',' << row.residualLinf;
    writeOverviewValues(row.values);
}
void PicStepSeries::writeProcess(int step, double time_s, const Process& row)
{
    process_ << step << ',' << time_s << ',' << row.driveEnergyJ << ',' << row.injectionEnergyJ
             << ',' << row.reactionEnergyLossJ << ',' << row.domainAbsorptionJ << ',' << row.embeddedAbsorptionJ
             << ',' << row.seeEmissionEnergyJ << ',' << row.ionImpactEnergyJ << ',' << row.displacementCurrentFraction
             << ',' << row.dominantSource << ',' << row.dominantSourceFraction << ',' << row.dominantBoundary
             << ',' << row.dominantBoundaryPowerW << ',' << row.dominantBoundaryFraction
             << ',' << row.dominantReactionChannel << ',' << row.dominantReactionPowerW << ',' << row.dominantReactionFraction
             << "\n";
    process_.flush();
}
void PicStepSeries::writeGeometry(int step, double time_s, const Geometry& row)
{
    geometry_ << step << ',' << time_s << ',' << row.updated << ',' << row.activeCutCellCount
              << ',' << row.embeddedSurfaceAreaM2 << ',' << row.movingRegionCount << ',' << row.movingLevelsetCount
              << ',' << row.movingRegionTotalAreaM2 << ',' << row.movingRegionTotalFlux
              << ',' << row.movingRegionMeanEnergyEV << ',' << row.movingRegionMeanCos
              << ',' << row.movingRegionMeanAbsSpeedMps << ',' << row.movingRegionMaxAbsSpeedMps
              << ',' << row.geometryGsaMm << ',' << row.geometryGadMm << ',' << row.geometryDaMm
              << ',' << row.fieldEnergyJ << ',' << row.domainAbsorptionJ << ',' << row.embeddedAbsorptionJ
              << ',' << row.reactionEnergyLossJ << "\n";
    geometry_.flush();
}

namespace {
constexpr const char* explicitPowerSI = "step,time_s,rG_L2,rG_Linf,rG_L2_per_active,rG_L2_rel_rhs,rG_Linf_rel_rhs,S_cf,S_eb,S_eb_rel_rhs,"
                        "S_cf_p50,S_cf_p90,S_cf_p99,S_eb_p50,S_eb_p90,S_eb_p99,"
                        "field_energy_J,field_dWdt_W,"
                        "kinetic_e_J,kinetic_ion_J,kinetic_neu_J,kinetic_total_J,"
                        "mcc_energy_loss_J,see_emitted_energy_J,ion_impact_energy_J,"
                        "dE_drive_J,dE_inj_J,dE_abs_domain_J,dE_abs_embedded_J,dE_abs_total_J,R_W_J,"
                        "sum_phiI_W,sum_phiI_dt_J\n";
constexpr const char* explicitPowerCode = "step,time,rG_L2,rG_Linf,rG_L2_per_active,rG_L2_rel_rhs,rG_Linf_rel_rhs,S_cf,S_eb,S_eb_rel_rhs,"
                        "S_cf_p50,S_cf_p90,S_cf_p99,S_eb_p50,S_eb_p90,S_eb_p99,"
                        "field_energy,field_dWdt,kinetic_e,kinetic_ion,kinetic_neu,kinetic_total,"
                        "mcc_energy_loss,see_emitted_energy,ion_impact_energy,"
                        "dE_drive,dE_inj,dE_abs_domain,dE_abs_embedded,dE_abs_total,R_W,"
                        "sum_phiI,sum_phiI_dt\n";
constexpr const char* implicitPowerSI = "step,time_s,rG_L2,rG_Linf,S_cf,S_eb,"
                     "S_cf_p50,S_cf_p90,S_cf_p99,S_eb_p50,S_eb_p90,S_eb_p99,"
                     "field_energy_J,field_dWdt_W,"
                             "kinetic_e_J,kinetic_ion_J,kinetic_neu_J,kinetic_total_J,kinetic_dWdt_W,"
                             "total_energy_J,total_dWdt_W,dE_drive_J,dE_inj_J,dE_abs_domain_J,dE_abs_embedded_J,dE_abs_total_J,R_W_J,"
                             "sum_phiI_W,sum_phiI_dt_J\n";
}

void PicStepSeries::openPower(const std::filesystem::path& path, int resumeStep, bool writeSI)
{
    powerWriteSI_ = writeSI;
    const bool hasExisting = TimeSeriesIO::openPreparedStepCsv(power_, path, resumeStep);
    if (power_)
    {
        if (!hasExisting)
            power_ << (mode_ == Mode::Implicit ? implicitPowerSI : writeSI ? explicitPowerSI : explicitPowerCode);
        power_ << std::setprecision(16);
    }
    else
        std::cerr << '[' << (mode_ == Mode::Explicit ? "ExplicitPIC" : "ImplicitPIC")
                  << "] WARNING: failed to create " << (mode_ == Mode::Explicit ? "power step report: " : "power report: ")
                  << path.string() << std::endl;
}

void PicStepSeries::writePower(int step, double time, const Power& row)
{
    power_ << step << ',' << time << ',' << row.residualL2 << ',' << row.residualLinf;
    if (mode_ == Mode::Explicit)
        power_ << ',' << row.residualL2PerActive << ',' << row.residualL2RelativeRhs << ',' << row.residualLinfRelativeRhs;
    power_ << ',' << row.coupling << ',' << row.embedded;
    if (mode_ == Mode::Explicit) power_ << ',' << row.embeddedRelativeRhs;
    power_ << ',' << row.couplingPercentiles[0] << ',' << row.couplingPercentiles[1] << ',' << row.couplingPercentiles[2]
           << ',' << row.embeddedPercentiles[0] << ',' << row.embeddedPercentiles[1] << ',' << row.embeddedPercentiles[2]
           << ',' << row.fieldEnergy << ',' << row.balance.fieldRate
           << ',' << row.kineticElectron << ',' << row.kineticIon << ',' << row.kineticNeutral << ',' << row.kineticTotal;
    if (mode_ == Mode::Explicit)
        power_ << ',' << row.transfers.reactionLossJ << ',' << row.transfers.seeEmissionJ << ',' << row.transfers.ionImpactJ;
    else
        power_ << ',' << row.balance.kineticRate << ',' << row.balance.totalJ << ',' << row.balance.totalRate;
    const bool writeDrive = mode_ == Mode::Implicit || powerWriteSI_;
    power_ << ',' << (writeDrive ? row.transfers.driveJ : 0.0) << ',' << row.transfers.injectionJ
           << ',' << row.transfers.domainAbsorptionJ << ',' << row.transfers.embeddedAbsorptionJ
           << ',' << row.balance.absorbedJ << ',' << row.balance.residualJ
           << ',' << (writeDrive ? row.potentialCurrent : 0.0) << ',' << (writeDrive ? row.potentialEnergy : 0.0) << "\n";
}

namespace { constexpr const char* explicitCurrentBalanceHeader = "step,time_s,"
    "Ii_inj_A,Ie_inj_A,I_inj_mag_A,"
    "I_out_domain_A,I_wall_eb_A,"
    "deltaI_mag_A,deltaI_mag_rel,"
    "I_inj_signed_A,I_out_domain_signed_A,I_wall_eb_signed_A,"
    "deltaI_signed_A,deltaI_signed_rel\n"; }

namespace { constexpr const char* implicitCurrentBalanceHeader = "step,time_s,Qsource_injected_C,Qsink_domain_C,Qsink_embedded_C,"
    "Qmacro_begin_C,Qmacro_end_C,dQmacro_step_C,dQphys_C,dQphys_rel\n"; }

void PicStepSeries::openCurrentBalance(const std::filesystem::path& path, int resumeStep)
{
    const bool hasExisting = TimeSeriesIO::openPreparedStepCsv(currentBalance_, path, resumeStep);
    if (currentBalance_)
    {
        if (!hasExisting) currentBalance_ << (mode_ == Mode::Explicit ? explicitCurrentBalanceHeader : implicitCurrentBalanceHeader);
        currentBalance_ << std::setprecision(16);
        if (mode_ == Mode::Implicit) { currentBalance_ << std::fixed; currentBalance_.flush(); }
    }
    else
        std::cerr << (mode_ == Mode::Explicit ? "[ExplicitPIC] WARNING: failed to create current-balance csv: "
                                             : "[ImplicitPIC] WARNING: failed to create current balance csv: ")
                  << path.string() << std::endl;
}
void PicStepSeries::writeCurrentBalance(int step, double timeS, const PicExplicitCurrentBalance& row)
{
    currentBalance_ << step << ','
                      << std::setprecision(16) << timeS << ','
                      << row.ionInjectionA << ',' << row.electronInjectionA << ',' << row.injectionMagnitudeA << ','
                      << row.domainMagnitudeA << ',' << row.embeddedMagnitudeA << ','
                      << row.magnitudeResidualA << ',' << row.magnitudeRelative << ','
                      << row.injectionSignedA << ',' << row.domainSignedA << ',' << row.embeddedSignedA << ','
                      << row.signedResidualA << ',' << row.signedRelative
                      << "\n";
    currentBalance_.flush();
}

void PicStepSeries::writeCurrentBalance(int step, double timeS, const PicStepParticleBalance::Result& row)
{
    currentBalance_ << step << ','
                      << timeS << ','
                      << row.sourceC << ','
                      << row.domainC << ','
                      << row.embeddedC << ','
                      << row.macroBeginC << ','
                      << row.macroEndC << ','
                      << (row.macroEndC - row.macroBeginC) << ','
                      << row.residualC << ','
                      << row.relative
                      << "\n";
}

void PicStepSeries::writeDetail(const std::filesystem::path& path, int step, double time,
    const ExplicitPIC::Diagnostics::PicDiagnostics::Result& r,
    const UnitsContext* units, bool writeSI, int resumeStep)
{
    using Result = ExplicitPIC::Diagnostics::PicDiagnostics::Result;
    const bool doSI = writeSI && units != nullptr;
    if (detail_.is_open() && (path != detailPath_ || doSI != detailWriteSI_))
        throw std::runtime_error("PIC detail CSV path or units changed within a session");
    if (!detail_.is_open())
    {
        // Stable, safe column keys for per-species diagnostics.
        //  - lower-case
        //  - map '+/-/*' to readable tokens
        //  - collapse underscores
        //  - ensure uniqueness by appending _2, _3, ... when needed
        auto safeKey = [](const std::string& raw) {
            std::string s;
            s.reserve(raw.size());
            for (char ch : raw)
            {
                if (ch >= 'A' && ch <= 'Z')
                    s.push_back(static_cast<char>(ch - 'A' + 'a'));
                else
                    s.push_back(ch);
            }
            std::string out;
            out.reserve(s.size() + 8);
            for (char ch : s)
            {
                if ((ch >= 'a' && ch <= 'z') || (ch >= '0' && ch <= '9'))
                {
                    out.push_back(ch);
                }
                else if (ch == '+')
                {
                    out += "plus";
                }
                else if (ch == '-')
                {
                    out += "minus";
                }
                else if (ch == '*')
                {
                    out += "star";
                }
                else
                {
                    out.push_back('_');
                }
            }

            // Collapse consecutive underscores.
            std::string collapsed;
            collapsed.reserve(out.size());
            bool prevUnderscore = false;
            for (char ch : out)
            {
                if (ch == '_')
                {
                    if (!prevUnderscore)
                        collapsed.push_back(ch);
                    prevUnderscore = true;
                }
                else
                {
                    collapsed.push_back(ch);
                    prevUnderscore = false;
                }
            }
            // Trim underscores.
            while (!collapsed.empty() && collapsed.front() == '_')
                collapsed.erase(collapsed.begin());
            while (!collapsed.empty() && collapsed.back() == '_')
                collapsed.pop_back();
            return collapsed.empty() ? std::string("unnamed") : collapsed;
        };

        std::ostringstream header;
        if (doSI)
        {
            header << "step,time_code,time_SI_s,nActive,rG_L2,rG_Linf,rG_L2_per_active,rG_L2_rel_rhs,rG_Linf_rel_rhs,S_cf,S_eb,S_eb_rel_rhs,"
                   "S_cf_p50,S_cf_p90,S_cf_p99,S_eb_p50,S_eb_p90,S_eb_p99,"
                   "energyW_code,energyW_SI_J,Qgrid,Qmacro,dQ,dQrel";
        }
        else
        {
            header << "step,time,nActive,rG_L2,rG_Linf,rG_L2_per_active,rG_L2_rel_rhs,rG_Linf_rel_rhs,S_cf,S_eb,S_eb_rel_rhs,"
                   "S_cf_p50,S_cf_p90,S_cf_p99,S_eb_p50,S_eb_p90,S_eb_p99,"
                   "energyW,Qgrid,Qmacro,dQ,dQrel";
        }

        // Optional per-species mean number density (volume-averaged) and electron energy summaries.
        const std::size_t ns = std::min(
            std::min(r.speciesNames.size(), r.speciesNmean_1pm3.size()),
            std::min(r.speciesMacroCount.size(), r.speciesPhysCount.size()));
        detailSpeciesNames_.assign(r.speciesNames.begin(), r.speciesNames.begin() + static_cast<std::ptrdiff_t>(ns));
        std::vector<std::string> keys;
        keys.reserve(ns);
        std::unordered_map<std::string, int> keyCount;
        keyCount.reserve(ns);
        for (std::size_t i = 0; i < ns; ++i)
        {
            const std::string base = safeKey(r.speciesNames[i]);
            int& c = keyCount[base];
            c += 1;
            const std::string key = (c == 1) ? base : (base + "_" + std::to_string(c));
            keys.push_back(key);
            header << ",n_" << key << "_mean_1pm3";
        }
        header << ",electron_mean_energy_eV,Te_e_mean_eV";
        header << ",macro_total,phys_total";
        for (const auto& key : keys)
            header << ",macro_" << key << ",phys_" << key;
        header << "\n";

        // Reconstruct the original species layout before trimming or appending.
        // Sanitized column names cannot be inverted to raw species identities;
        // require the same initial layout and unit mode instead of guessing.
        const std::string expected = header.str();
        if (resumeStep > 0)
        {
            std::ifstream existing(path);
            std::string first;
            if (std::getline(existing, first) && first + "\n" != expected)
                throw std::runtime_error("PIC detail CSV header mismatch: " + path.string());
        }
        const bool hasExisting = TimeSeriesIO::openPreparedStepCsv(detail_, path, resumeStep);
        if (!detail_) return; // Retain the optional diagnostic's open-failure policy.
        detailPath_ = path;
        detailWriteSI_ = doSI;
        if (!hasExisting) detail_ << expected;
    }
    auto& out = detail_;
    if (!out) return;
    // Helper: find species index by name in current Result.
    auto findSpeciesIndex = [](const Result& r, const std::string& nm) -> int {
        for (std::size_t i = 0; i < r.speciesNames.size(); ++i)
        {
            if (r.speciesNames[i] == nm)
                return static_cast<int>(i);
        }
        return -1;
    };

    out << std::setprecision(16);
    out << step << ',';
    if (doSI)
    {
        const double tSI = units->TimeCodeToSI(time);
        out << time << ',' << tSI << ',';
    }
    else
    {
        out << time << ',';
    }
    out << r.nActive << ',' << r.rG_L2 << ',' << r.rG_Linf << ','
        << r.rG_L2_per_active << ',' << r.rG_L2_rel_rhs << ',' << r.rG_Linf_rel_rhs << ','
        << r.S_cf << ',' << r.S_eb << ',' << r.S_eb_rel_rhs << ','
        << r.S_cf_p50 << ',' << r.S_cf_p90 << ',' << r.S_cf_p99 << ','
        << r.S_eb_p50 << ',' << r.S_eb_p90 << ',' << r.S_eb_p99 << ',';
    if (doSI)
    {
        // r.energyW is accumulated as sum(0.5 * Q[C] * phi_code). Convert to Joules by Phi0.
        out << r.energyW << ',' << (r.energyW * units->scales().Phi0_V) << ',';
    }
    else
    {
        out << r.energyW << ',';
    }

    out << r.Qgrid << ',' << r.Qmacro << ',' << r.dQ << ',' << r.dQrel;

    // IMPORTANT: keep row schema identical to header by writing only the species
    // set captured when the header was created.
    const std::size_t nsHeader = detailSpeciesNames_.size();
    for (std::size_t i = 0; i < nsHeader; ++i)
    {
        const int idx = findSpeciesIndex(r, detailSpeciesNames_[i]);
        const double v = (idx >= 0 && static_cast<std::size_t>(idx) < r.speciesNmean_1pm3.size())
                             ? r.speciesNmean_1pm3[static_cast<std::size_t>(idx)]
                             : 0.0;
        out << ',' << v;
    }
    out << ',' << r.electron_mean_energy_eV << ',' << r.Te_e_mean_eV;
    out << ',' << r.macro_total << ',' << r.phys_total;
    for (std::size_t i = 0; i < nsHeader; ++i)
    {
        const int idx = findSpeciesIndex(r, detailSpeciesNames_[i]);
        const long long macro = (idx >= 0 && static_cast<std::size_t>(idx) < r.speciesMacroCount.size())
                                    ? r.speciesMacroCount[static_cast<std::size_t>(idx)]
                                    : 0;
        const double phys = (idx >= 0 && static_cast<std::size_t>(idx) < r.speciesPhysCount.size())
                                 ? r.speciesPhysCount[static_cast<std::size_t>(idx)]
                                 : 0.0;
        out << ',' << macro << ',' << phys;
    }
    out << "\n";
    if (++detailRowsSinceFlush_ >= 64)
    {
        out.flush();
        detailRowsSinceFlush_ = 0;
    }
}

void PicStepSeries::openParticleManagement(const std::filesystem::path& path, int resumeStep)
{
    const bool existing = TimeSeriesIO::openPreparedStepCsv(particleManagement_, path, resumeStep);
    if (particleManagement_)
    {
        if (!existing) particleManagement_ << "step,time_s,kind,name,"
            "macro_before,phys_before,macro_after,phys_after,"
            "merged_pairs,split_count,split_emergency_count,split_count_based_count,"
            "spawn_in,spawn_kept,spawn_drop,keep_prob,"
            "weight_p50,weight_p90,weight_p99,weight_max,weight_top1_phys_frac\n";
        particleManagement_ << std::setprecision(16);
        if (mode_ == Mode::Implicit) { particleManagement_ << std::fixed; particleManagement_.flush(); }
    }
    else std::cerr << '[' << (mode_ == Mode::Explicit ? "ExplicitPIC" : "ImplicitPIC")
        << "] WARNING: failed to create particle management " << (mode_ == Mode::Explicit ? "report: " : "csv: ")
        << path.string() << std::endl;
}

void PicStepSeries::writeSpawnManagement(int step, double timeS, std::string_view name,
    const PicSpawnManagement::Stats& stats, double keepProbability)
{
    particleManagement_ << step << ',' << std::setprecision(16) << timeS << ",spawn," << name << ','
        << stats.in << ',' << stats.physIn << ',' << stats.kept << ',' << stats.physKept << ','
        << 0 << ',' << 0 << ',' << 0 << ',' << 0 << ','
        << stats.in << ',' << stats.kept << ',' << stats.drop << ',' << keepProbability << ','
        << 0.0 << ',' << 0.0 << ',' << 0.0 << ',' << 0.0 << ',' << 0.0 << "\n";
}

void PicStepSeries::writeParticleManagement(int step, double timeS, std::string_view name,
    const ParticleManagementReport& report)
{
    const auto& counts = report.counts;
    const auto& weights = report.weights;
    particleManagement_ << step << ',' << std::setprecision(16) << timeS << ",species," << name << ','
        << report.macroBefore << ',' << report.physBefore << ',' << report.macroAfter << ',' << report.physAfter << ','
        << counts.mergedPairs << ',' << counts.splitCount << ',' << counts.splitEmergencyCount << ',' << counts.splitCountBasedCount << ','
        << 0 << ',' << 0 << ',' << 0 << ',' << 1.0 << ','
        << weights.p50 << ',' << weights.p90 << ',' << weights.p99 << ',' << weights.max << ',' << weights.top1PhysFrac << "\n";
}

namespace
{
// Validate the entire history before changing either its schema or resume tail.
// Only the known explicit legacy layout is migrated; numeric tokens stay intact.
void migrateExplicitPerformance(const std::filesystem::path& path)
{
    if (!std::filesystem::exists(path)) return;
    const auto fail = [&]() {
        throw std::runtime_error("Incompatible or unreadable explicit performance CSV: " + path.string());
    };
    std::ifstream input(path);
    if (!input) fail();
    std::string line;
    if (!std::getline(input, line)) { if (!input.eof()) fail(); return; }
    std::string expected(explicitPerformanceHeader);
    expected.pop_back(); // trailing newline
    std::string legacy = expected;
    legacy.replace(5, 9, "time_s");
    const bool oldHeader = line == legacy;
    if (!oldHeader && line != expected) fail();
    std::vector<std::string> lines{expected};
    while (std::getline(input, line))
    {
        const auto comma = line.find(',');
        if (comma == std::string::npos) fail();
        int step = -1;
        const auto parsed = std::from_chars(line.data(), line.data() + comma, step);
        if (parsed.ec != std::errc{} || parsed.ptr != line.data() + comma || step < 0) fail();
        const auto columns = std::count(line.begin(), line.end(), ',') + 1;
        if (oldHeader && step == 0 && columns == 17)
        {
            std::size_t position = 0;
            for (int column = 0; column < 10; ++column) position = line.find(',', position) + 1;
            line.insert(position, "0,0,0,");
        }
        else if (columns != 20) fail();
        lines.push_back(line);
    }
    if (!input.eof()) fail();
    input.close();
    if (!oldHeader) return;
    std::ofstream output(path, std::ios::trunc);
    for (const auto& record : lines) output << record << '\n';
    output.close();
    if (!output) throw std::runtime_error("Failed to migrate explicit performance CSV: " + path.string());
}
}

void PicStepSeries::openPerformance(const std::filesystem::path& path, int resumeStep)
{
    if (mode_ == Mode::Explicit && resumeStep > 0) migrateExplicitPerformance(path);
    open(performance_, path, resumeStep,
        mode_ == Mode::Explicit ? explicitPerformanceHeader : implicitPerformanceHeader,
        mode_ == Mode::Explicit ? "timing" : "perf");
}
void PicStepSeries::openCouplingFallback(const std::filesystem::path& path, int resumeStep)
{ open(fallback_, path, resumeStep, fallbackHeader, "fallback"); }
void PicStepSeries::flushPerformance()
{
    if (performanceReady()) performance_.flush();
    if (couplingFallbackReady()) fallback_.flush();
}
void PicStepSeries::writePerformance(int step, double timeCode, const ExplicitPerformance& row)
{
    performance_ << step << ',' << timeCode
        << ',' << row.field.buildWeights_s << ',' << row.field.deposit_s << ',' << row.field.rhs_s;
    for (double seconds : row.poissonSeconds) performance_ << ',' << seconds;
    performance_ << ',' << row.field.evalE_s << ',' << row.field.amr_indicator_s
        << ',' << row.amrAdaptSeconds << ',' << row.managementSeconds
        << ',' << row.gatherPushSeconds << ',' << row.outputSeconds << ',' << row.totalSeconds
        << ',' << row.particles << ',' << row.poissonN << ',' << row.poissonNnz
        << ',' << (row.usedSuperLU ? 1 : 0) << "\n";
    performance_.flush();
}
void PicStepSeries::writeCouplingFallback(int step, double timeS, const PicImplicitStepTiming& stepTiming,
    const PicImplicitTimingSummary& summary)
{
    fallback_ << step << ',' << timeS
                << ',' << stepTiming.buildWeightsN_cwFallbackCount
                << ',' << stepTiming.buildWeightsN_particleCount
                << ',' << summary.buildFallbackFraction
                << ',' << stepTiming.orbitAvg_sampleCwFallbackCount
                << ',' << stepTiming.orbitAvg_sampleCount
                << ',' << summary.orbitFallbackFraction
                << ',' << summary.totalFallbackCount
                << "\n";
}
void PicStepSeries::writePerformance(int step, double timeS, const PicImplicitStepTiming& stepTiming,
    const PicImplicitTimingSummary& summary, const ImplicitPerformance& row)
{
    const double npAsDouble = static_cast<double>(std::max<std::uint64_t>(1u, row.particles));
    const double predictorSegPerParticle =
        static_cast<double>(stepTiming.predictor_segments_total) / npAsDouble;
    const double predictorCoarsenHitsPerParticle =
        static_cast<double>(stepTiming.predictor_coarsen_hits) / npAsDouble;
    const double predictorCoarsenParentClimbPerHit =
        (stepTiming.predictor_coarsen_hits > 0u)
            ? (static_cast<double>(stepTiming.predictor_coarsen_parent_climb) /
               static_cast<double>(stepTiming.predictor_coarsen_hits))
            : 0.0;
    const double orbitSamplesPerParticle =
        static_cast<double>(stepTiming.orbitAvg_sampleCount) / npAsDouble;
    const double orbitSampleCacheNnzPerSample =
        (stepTiming.orbitAvg_sampleCount > 0u)
            ? (static_cast<double>(stepTiming.orbitAvg_sampleCache_nnz) /
               static_cast<double>(stepTiming.orbitAvg_sampleCount))
            : 0.0;
    const double orbitAdaptAcceptPerParticle =
        static_cast<double>(stepTiming.orbitAvg_adaptAcceptCount) / npAsDouble;
    const double orbitAdaptRefinePerParticle =
        static_cast<double>(stepTiming.orbitAvg_adaptRefineCount) / npAsDouble;
    const double orbitAdaptHitMaxDepthPerParticle =
        static_cast<double>(stepTiming.orbitAvg_adaptHitMaxDepthCount) / npAsDouble;
    const double orbitAdaptMeanAcceptedDepth =
        (stepTiming.orbitAvg_adaptAcceptCount > 0u)
            ? (static_cast<double>(stepTiming.orbitAvg_adaptAcceptedDepthSum) /
               static_cast<double>(stepTiming.orbitAvg_adaptAcceptCount))
            : 0.0;

    performance_ << step << ','
            << timeS << ','
            << row.particles << ','
            << stepTiming.total_s << ','
            << stepTiming.mcc_pre_s << ','
            << stepTiming.buildWeightsN_s << ','
            << stepTiming.predictorSegment_s << ','
            << stepTiming.orbitAvgWeights_s << ','
            << stepTiming.polPoisson_total_s << ','
            << stepTiming.pm_s << ','
            << stepTiming.amr_s << ','
            << stepTiming.amr_indicator_s << ','
            << stepTiming.amr_adapt_s << ','
            << stepTiming.amr_faceGraph_s << ','
            << stepTiming.amr_cutMetrics_s << ','
            << stepTiming.amr_reinitLevelset_s << ','
            << stepTiming.amr_geom_initializeLevelSets_s << ','
            << stepTiming.amr_geom_pdeReinitLevelset_s << ','
            << stepTiming.amr_geom_refreshDerived_s << ','
            << stepTiming.amr_geom_checkCutInfo_s << ','
            << stepTiming.amr_geom_rebuildComposite_s << ','
            << stepTiming.amr_geom_markPhases_s << ','
            << stepTiming.amr_fieldRebuild_s << ','
            << stepTiming.packedFieldOutput_s << ','
            << (row.outputStep ? 1 : 0) << ','
            << (row.amrAdaptStep ? 1 : 0) << ','
            << (row.amrSkippedByGate ? 1 : 0) << ','
            << row.poissonN << ','
            << row.poissonNnz << ','
            << row.solverMode << ','
            << predictorSegPerParticle << ','
            << predictorCoarsenHitsPerParticle << ','
            << predictorCoarsenParentClimbPerHit << ','
            << orbitSamplesPerParticle << ','
            << orbitSampleCacheNnzPerSample << ','
            << orbitAdaptAcceptPerParticle << ','
            << orbitAdaptRefinePerParticle << ','
            << orbitAdaptHitMaxDepthPerParticle << ','
            << orbitAdaptMeanAcceptedDepth << ','
            << summary.top[0].name << ',' << summary.top[0].t << ','
            << summary.top[1].name << ',' << summary.top[1].t << ','
            << summary.top[2].name << ',' << summary.top[2].t
            << "\n";
}

PicImplicitConsole::PicImplicitConsole(const PicConfig& cfg, const ImplicitPicOverrides& overrides,
    int firstStep, int endStep, int fieldEvery, bool infoEnabled)
    : startStep(firstStep), steps(endStep), fieldsEvery(fieldEvery), consoleInfo(infoEnabled)
{
    const std::string outputV2ConsoleFormat = PhysicsDiagnostics::lowerCopy(cfg.output.v2.console.format);
    fallbackFracWarnTh = std::max(0.0, overrides.consoleFallbackFracWarn.value_or(1e-3));
    fallbackFracDebounceSteps =
        std::max(1, overrides.consoleFallbackDebounceSteps.value_or(3));
    fallbackPrintRecover =
        overrides.consoleFallbackPrintRecover.value_or(true);
    const double dQrelWarnCfg = std::max(0.0, cfg.output.console.dQrelWarn);
    dQrelWarnTh = std::max(0.0, overrides.consoleDqrelWarn.value_or(dQrelWarnCfg));
    const double dQrelRecoverCfg =
        (cfg.output.console.dQrelRecover >= 0.0) ? cfg.output.console.dQrelRecover : (0.7 * dQrelWarnTh);
    const double dQrelRecoverRaw =
        std::max(0.0, overrides.consoleDqrelRecover.value_or(dQrelRecoverCfg));
    dQrelRecoverTh =
        (dQrelWarnTh > 0.0) ? std::min(dQrelRecoverRaw, std::nextafter(dQrelWarnTh, 0.0)) : dQrelRecoverRaw;
    dQrelDebounceSteps =
        std::max(1, overrides.consoleDqrelDebounceSteps.value_or(cfg.output.console.dQrelDebounceSteps));
    dQrelEventCooldownSteps =
        std::max(0, overrides.consoleDqrelEventCooldownSteps.value_or(cfg.output.console.dQrelEventCooldownSteps));
    dQrelPrintRecover =
        overrides.consoleDqrelPrintRecover.value_or(cfg.output.console.dQrelPrintRecover);
    rGL2PerActiveWarnTh =
        std::max(0.0, overrides.consoleRgL2Warn.value_or(1e2));
    amrSkipWarnStreakTh = overrides.consoleAmrSkipStreakWarn.value_or(5);
    physicsSummaryEnabled =
        overrides.consolePhysicsEnabled.value_or(cfg.output.console.physicsEnabled);
    physicsSummaryEvery =
        std::max(0, overrides.consolePhysicsEvery.value_or(cfg.output.console.every));
    if (physicsSummaryEvery <= 0)
        physicsSummaryEvery = std::max(1, cfg.output.console.every);
    physicsShowCharge =
        overrides.consolePhysicsShowCharge.value_or(cfg.output.console.physicsShowCharge);
    physicsShowCurrent =
        overrides.consolePhysicsShowCurrent.value_or(cfg.output.console.physicsShowCurrent);
    physicsShowWall =
        overrides.consolePhysicsShowWall.value_or(cfg.output.console.physicsShowWall);
    physicsShowMCC =
        overrides.consolePhysicsShowMcc.value_or(cfg.output.console.physicsShowMCC);
    physicsShowPerf =
        overrides.consolePhysicsShowPerf.value_or(cfg.output.console.physicsShowPerf);
    physicsDetail = PhysicsDiagnostics::lowerCopy(cfg.output.console.physicsDetail);
    summaryOneLine = (outputV2ConsoleFormat != "multiline");
}

std::uint64_t PicImplicitConsole::writeStep(int step, double t_si_step, double dt_si, bool summaryStepDefault,
    const PicImplicitStepTiming& stepTiming, const PicImplicitTimingSummary& timingSummary,
    const Sample& sample, std::ostream& output)
{

    auto fmtCompact = [](double v, int p = 3) -> std::string {
        std::ostringstream oss;
        if (!std::isfinite(v))
            return "nan";
        oss << std::setprecision(p) << std::scientific << v;
        return oss.str();
    };

    const auto& stepCharge = sample.charge;
    const auto& species = sample.species;
    const auto& mccStepStats = sample.collisions;
    const bool mccEnabled = sample.events.collisionsEnabled;
    const bool isAmrIndicatorStep = sample.events.amrIndicator;
    const bool amrSkippedByGateStep = sample.events.amrSkippedByGate;
    const bool didAmrAdaptStep = sample.events.amrAdapted;
    const auto& electrodeGroups = sample.electrodes.groups();
    const auto& consoleElectrodeGroupBuckets = sample.electrodes.buckets();
    const auto& wallCurrentEBySideStep_Am2 = sample.walls.electronCurrent;
    const auto& wallCurrentIBySideStep_Am2 = sample.walls.ionCurrent;
    const auto& wallOutChargeEBySideStep_C = sample.walls.electronCharge;
    const auto& wallOutChargeIBySideStep_C = sample.walls.ionCharge;
    const auto& wallOutMacroEBySideStep = sample.losses.domainMacroE;
    const auto& wallOutMacroIBySideStep = sample.losses.domainMacroI;
    const auto& wallOutMacroEUnknownStep = sample.losses.unknownE;
    const auto& wallOutMacroIUnknownStep = sample.losses.unknownI;
    const auto& ebOutMacroBySpeciesStep = sample.losses.embeddedMacro;
    const auto& embeddedLossByBoundaryIdStep = sample.losses.embeddedCounts.byId;
    const double dE_abs_embedded_J = sample.losses.embeddedEnergyJ;
    const bool summaryStepPhysics =
        physicsSummaryEnabled
            ? ((((step + 1) % std::max(1, physicsSummaryEvery)) == 0) || ((step + 1) == steps))
            : summaryStepDefault;

    const double totalFallbackFrac = timingSummary.totalFallbackFraction;
    const std::string top1Name = timingSummary.top[0].name;
    const std::string top2Name = timingSummary.top[1].name;
    const std::string top3Name = timingSummary.top[2].name;
    const double top1S = timingSummary.top[0].t;
    const double top2S = timingSummary.top[1].t;
    const double top3S = timingSummary.top[2].t;

    if (isAmrIndicatorStep)
    {
        if (amrSkippedByGateStep)
        {
            ++amrSkipGateStreak;
            if (consoleInfo && amrSkipGateStreak == amrSkipWarnStreakTh)
            {
                output << "[ImplicitPIC][Event][WARN] AMR gate skipped for "
                          << amrSkipGateStreak << " consecutive indicator steps"
                          << " at step=" << (step + 1) << std::endl;
            }
        }
        else
        {
            amrSkipGateStreak = 0;
        }
    }

    if (haveDiagLast)
    {
        const double dQMonitor =
            std::isfinite(stepCharge.stageRelative) ? stepCharge.stageRelative : diagLast.dQrel;
        const double dQAbs = std::abs(dQMonitor);
        const bool dQFinite = std::isfinite(dQAbs);
        const bool dQAboveWarn = dQFinite && (dQAbs > dQrelWarnTh);
        const bool dQBelowRecover = dQFinite && (dQAbs < dQrelRecoverTh);
        const int stepOut = step + 1;
        const bool dqrelCooldownReady =
            (dQrelEventCooldownSteps <= 0) ||
            (lastDqrelEventStep < 0) ||
            ((stepOut - lastDqrelEventStep) >= dQrelEventCooldownSteps);

        if (!dQrelWarnActive)
        {
            dQrelRecoverConsec = 0;
            dQrelWarnConsec = dQAboveWarn ? (dQrelWarnConsec + 1) : 0;
            if (dQrelWarnConsec >= dQrelDebounceSteps)
            {
                dQrelWarnActive = true;
                dQrelWarnConsec = 0;
                if (consoleInfo && dqrelCooldownReady)
                {
                    const char* tag = std::isfinite(stepCharge.stageRelative) ? "|dQrel_stage|" : "|dQrel|";
                    output << "[ImplicitPIC][Event][WARN] " << tag << "="
                              << std::scientific << std::setprecision(3) << dQAbs
                              << " exceeded threshold=" << dQrelWarnTh
                              << " at step=" << stepOut
                              << " dQphys_rel=" << fmtCompact(stepCharge.relative, 3)
                              << std::defaultfloat << std::setprecision(6) << std::endl;
                    lastDqrelEventStep = stepOut;
                }
            }
        }
        else
        {
            dQrelWarnConsec = 0;
            dQrelRecoverConsec = dQBelowRecover ? (dQrelRecoverConsec + 1) : 0;
            if (dQrelRecoverConsec >= dQrelDebounceSteps)
            {
                dQrelWarnActive = false;
                dQrelRecoverConsec = 0;
                if (consoleInfo && dQrelPrintRecover && dqrelCooldownReady)
                {
                    output << "[ImplicitPIC][Event][RECOVER] dQrel recovered at step=" << stepOut << std::endl;
                    lastDqrelEventStep = stepOut;
                }
            }
        }

        const bool rGWarnNow = std::isfinite(diagLast.rG_L2_per_active) &&
                               (diagLast.rG_L2_per_active > rGL2PerActiveWarnTh);
        if (consoleInfo && rGWarnNow && !rGL2WarnActive)
        {
            output << "[ImplicitPIC][Event][WARN] rG_L2_per_active="
                      << std::scientific << std::setprecision(3) << diagLast.rG_L2_per_active
                      << " exceeded threshold=" << rGL2PerActiveWarnTh
                      << " at step=" << (step + 1)
                      << std::defaultfloat << std::setprecision(6) << std::endl;
        }
        else if (consoleInfo && !rGWarnNow && rGL2WarnActive)
        {
            output << "[ImplicitPIC][Event][RECOVER] rG_L2_per_active recovered at step="
                      << (step + 1) << std::endl;
        }
        rGL2WarnActive = rGWarnNow;
    }

    const bool fallbackWarnNow = (fallbackFracWarnTh > 0.0) && std::isfinite(totalFallbackFrac) &&
                                 (totalFallbackFrac > fallbackFracWarnTh);
    if (!fallbackWarnActive)
    {
        fallbackRecoverConsec = 0;
        fallbackWarnConsec = fallbackWarnNow ? (fallbackWarnConsec + 1) : 0;
        if (fallbackWarnConsec >= fallbackFracDebounceSteps)
        {
            fallbackWarnActive = true;
            fallbackWarnConsec = 0;
            if (consoleInfo)
            {
                output << "[ImplicitPIC][Event][WARN] coupling_fallback_frac="
                          << std::scientific << std::setprecision(3) << totalFallbackFrac
                          << " exceeded threshold=" << fallbackFracWarnTh
                          << " at step=" << (step + 1)
                          << std::defaultfloat << std::setprecision(6) << std::endl;
            }
        }
    }
    else
    {
        fallbackWarnConsec = 0;
        fallbackRecoverConsec = fallbackWarnNow ? 0 : (fallbackRecoverConsec + 1);
        if (fallbackRecoverConsec >= fallbackFracDebounceSteps)
        {
            fallbackWarnActive = false;
            fallbackRecoverConsec = 0;
            if (consoleInfo && fallbackPrintRecover)
            {
                output << "[ImplicitPIC][Event][RECOVER] coupling_fallback_frac recovered at step="
                          << (step + 1) << std::endl;
            }
        }
    }

    recentStepCost_s.push_back(stepTiming.total_s);
    recentStepCostSum_s += stepTiming.total_s;
    while (recentStepCost_s.size() > 20u)
    {
        recentStepCostSum_s -= recentStepCost_s.front();
        recentStepCost_s.pop_front();
    }

    std::uint64_t nP = 0u;
    for (const auto& sp : species)
        nP += static_cast<std::uint64_t>(sp.particle.size());
    if (consoleInfo && summaryStepPhysics && haveDiagLast)
    {
        const double avgStep_s = (recentStepCost_s.empty())
                                     ? stepTiming.total_s
                                     : (recentStepCostSum_s / static_cast<double>(recentStepCost_s.size()));
        const int stepsRemain = std::max(0, steps - (step + 1));
        const double eta_s = avgStep_s * static_cast<double>(stepsRemain);
        const char* amrFlag = didAmrAdaptStep ? "adapt" : (amrSkippedByGateStep ? "skip" : "none");
        const char* outFlag = (((step + 1) % fieldsEvery) == 0) ? "out" : "-";
        const double dQRelOut = std::isfinite(stepCharge.stageRelative) ? stepCharge.stageRelative : diagLast.dQrel;

        const auto values = summarizePicConsole(sample.mesh, sample.potential, sample.units,
            species, sample.defaultWeights, mccStepStats);
        const auto& [phiMinV, phiMaxV, eMeanE_eV, Te_eV, eTailFrac, nElastic, nExc, nIoniz, nAttach, nCx, nOther, mccE, mccI, mccN] = values;

        auto fmt4u = [](const std::array<std::uint64_t, 4>& a) {
            std::ostringstream oss;
            oss << a[0] << ',' << a[1] << ',' << a[2] << ',' << a[3];
            return oss.str();
        };
        auto fmt4d = [&](const std::array<double, 4>& a, int p) {
            std::ostringstream oss;
            oss << fmtCompact(a[0], p) << ','
                << fmtCompact(a[1], p) << ','
                << fmtCompact(a[2], p) << ','
                << fmtCompact(a[3], p);
            return oss.str();
        };

        auto fmtSides = [&](const std::array<double, 4>& a, int p) {
            std::ostringstream oss;
            oss << "L=" << fmtCompact(a[0], p)
                << " R=" << fmtCompact(a[1], p)
                << " B=" << fmtCompact(a[2], p)
                << " T=" << fmtCompact(a[3], p);
            return oss.str();
        };
        auto fmtEventTag = [&](const char* amrFlagLocal, const char* outFlagLocal) {
            std::string evt;
            if (amrFlagLocal != nullptr && std::string(amrFlagLocal) != "none")
            {
                evt += "amr-";
                evt += amrFlagLocal;
            }
            if (outFlagLocal != nullptr && std::string(outFlagLocal) == "out")
            {
                if (!evt.empty())
                    evt += ',';
                evt += "output";
            }
            if (evt.empty())
                evt = "advance";
            return evt;
        };
        auto classifyStage = [&](int stepNumber, const char* amrFlagLocal) {
            const int stepsSinceStart = std::max(0, stepNumber - startStep);
            if (stepsSinceStart <= 1)
                return std::string("field-init");
            if (stepsSinceStart <= 10)
                return std::string("startup");
            const bool quasiSteady =
                ((amrFlagLocal == nullptr) || std::string(amrFlagLocal) == "none") &&
                std::abs(dQRelOut) <= 1e-3 &&
                std::abs(stepCharge.relative) <= 1e-3;
            return quasiSteady ? std::string("quasi-steady") : std::string("transient");
        };
        auto fmtGroupCircuit = [&](double dtSeconds) {
            const bool highDetail = (physicsDetail == "high");
            std::ostringstream oss;
            oss << "Igrp_A{";
            for (std::size_t gi = 0; gi < electrodeGroups.size(); ++gi)
            {
                if (gi)
                    oss << ' ';
                const auto& bucket = consoleElectrodeGroupBuckets[gi];
                const double invDt = (dtSeconds > 0.0) ? (1.0 / dtSeconds) : 0.0;
                oss << electrodeGroups[gi].name << '=' << fmtCompact(bucket.chargeC * invDt, 3);
                if (highDetail)
                {
                    oss << "(e=" << fmtCompact(bucket.electronChargeC * invDt, 3)
                        << " i=" << fmtCompact(bucket.ionChargeC * invDt, 3) << ')';
                }
            }
            oss << '}';
            return oss.str();
        };

        std::ostringstream particleSummary;
        const std::size_t nSpeciesDiag =
            std::min(diagLast.speciesNames.size(),
                     std::min(diagLast.speciesMacroCount.size(), diagLast.speciesPhysCount.size()));
        for (std::size_t si = 0; si < nSpeciesDiag; ++si)
        {
            if (si > 0)
                particleSummary << ' ';
            particleSummary << diagLast.speciesNames[si]
                            << '=' << diagLast.speciesMacroCount[si]
                            << '(' << fmtCompact(diagLast.speciesPhysCount[si], 3) << ')';
        }

        const std::string detail =
            (physicsDetail == "minimal" || physicsDetail == "medium" || physicsDetail == "high")
                ? physicsDetail
                : std::string("medium");
        const bool detailMinimal = (detail == "minimal");
        const bool detailHigh = (detail == "high");
        const bool showParticles = !detailMinimal;
        const bool showCurrent = physicsShowCurrent;
        const bool showWall = !detailMinimal && physicsShowWall;
        const bool showEmbedded = detailHigh &&
                                  (sample.events.embeddedReady || std::abs(stepCharge.embeddedC) > 0.0);
        const bool showMccSplit = detailHigh && physicsShowMCC && mccEnabled;

        const double tailThr_eV = PicConsoleValues::tailThresholdEV;
        using MccSplit = PicConsoleValues::MccSplit;

        double jTot_Apm2 = 0.0;
        if (sample.current.ready && std::isfinite(sample.current.total))
        {
            jTot_Apm2 = sample.current.total;
        }
        else
        {
            for (int s = 0; s < 4; ++s)
                jTot_Apm2 += wallCurrentEBySideStep_Am2[static_cast<std::size_t>(s)] +
                             wallCurrentIBySideStep_Am2[static_cast<std::size_t>(s)];
        }

        std::uint64_t wallLossMacroTotal = 0u;
        for (int s = 0; s < 4; ++s)
            wallLossMacroTotal += wallOutMacroEBySideStep[static_cast<std::size_t>(s)] +
                                  wallOutMacroIBySideStep[static_cast<std::size_t>(s)];
        long long embeddedLossMacroTotal = 0ll;
        for (long long value : ebOutMacroBySpeciesStep)
            embeddedLossMacroTotal += value;
        const long long sourceSpawnMacroTotal = sample.events.sourceSpawnCount;
        const long long mccSpawnMacroTotal = sample.events.collisionSpawnCount;

        std::vector<std::pair<std::string, std::string>> sections;
        sections.emplace_back(
            "Progress",
            "step=" + std::to_string(step + 1) + "/" + std::to_string(steps) +
                " t_s=" + fmtCompact(t_si_step, 3) +
                " dt_s=" + fmtCompact(dt_si, 3) +
                " np=" + std::to_string(nP) +
                " stage=" + classifyStage(step + 1, amrFlag) +
                " evt=" + fmtEventTag(amrFlag, outFlag) +
                " ETA=" + formatPicEta(eta_s));

        if (showParticles)
            sections.emplace_back("Particles", particleSummary.str());

        {
            std::string stateLine;
            if (physicsShowCharge)
            {
                stateLine += "phi[V]={[" + fmtCompact(phiMinV, 3) + "," + fmtCompact(phiMaxV, 3) + "]}";
                stateLine += " ";
            }
            stateLine +=
                "eMeanE_eV=" + fmtCompact(eMeanE_eV, 3) +
                " Te_eV=" + fmtCompact(Te_eV, 3) +
                " tail_e_gt_" + fmtCompact(tailThr_eV, 3) + "eV=" + fmtCompact(eTailFrac, 3);
            sections.emplace_back("State", stateLine);
        }

        if (showCurrent)
        {
            std::string currentLine;
            if (!electrodeGroups.empty())
            {
                currentLine = fmtGroupCircuit(dt_si) +
                              " Jtot_Apm2=" + fmtCompact(jTot_Apm2, 3);
            }
            else
            {
                currentLine =
                    "Jtot_Apm2=" + fmtCompact(jTot_Apm2, 3) +
                    " Je_Apm2{" + fmtSides(wallCurrentEBySideStep_Am2, 2) + "}" +
                    " Ji_Apm2{" + fmtSides(wallCurrentIBySideStep_Am2, 2) + "}";
            }
            if (detailHigh && sample.current.ready)
            {
                currentLine += " Jcond_Apm2=" + fmtCompact(sample.current.conduction, 3) +
                               " Jdisp_Apm2=" + fmtCompact(sample.current.displacement, 3);
            }
            sections.emplace_back("Circuit", currentLine);
        }

        {
            auto fmtMccCompact = [&]() {
                if (!mccEnabled)
                    return std::string("mcc=off");
                const std::array<std::pair<const char*, std::uint64_t>, 6> counts = {{
                    {"el", nElastic},
                    {"ex", nExc},
                    {"ion", nIoniz},
                    {"att", nAttach},
                    {"cx", nCx},
                    {"other", nOther},
                }};
                std::ostringstream oss;
                oss << "mcc{";
                bool first = true;
                for (const auto& entry : counts)
                {
                    if (entry.second == 0)
                        continue;
                    if (!first)
                        oss << ' ';
                    oss << entry.first << '=' << entry.second;
                    first = false;
                }
                if (first)
                    oss << '0';
                oss << '}';
                return oss.str();
            };
            std::string sourceSinkLine =
                "N[s/m/d/e]=" + std::to_string(sourceSpawnMacroTotal) + '/' +
                std::to_string(mccSpawnMacroTotal) + '/' +
                std::to_string(wallLossMacroTotal) + '/' +
                std::to_string(embeddedLossMacroTotal) +
                " Q[s/d/e]=" + fmtCompact(stepCharge.sourceC, 3) + '/' +
                fmtCompact(stepCharge.domainC, 3) + '/' +
                fmtCompact(stepCharge.embeddedC, 3);
            if (physicsShowMCC)
                sourceSinkLine += ' ' + fmtMccCompact();
            sections.emplace_back("SourceSink", sourceSinkLine);
        }

        {
            std::string healthLine =
                "dQ_rel=" + fmtCompact(dQRelOut, 3) +
                " dQphys_rel=" + fmtCompact(stepCharge.relative, 3) +
                " rG_L2_per_active=" + fmtCompact(diagLast.rG_L2_per_active, 3) +
                " rG_Linf=" + fmtCompact(diagLast.rG_Linf, 3) +
                " S_cf=" + fmtCompact(diagLast.S_cf, 3) +
                " S_eb=" + fmtCompact(diagLast.S_eb, 3);
            sections.emplace_back("Health", healthLine);
        }

        if (showWall)
        {
            std::string wallLine =
                "out_e[L/R/B/T]=" + fmt4u(wallOutMacroEBySideStep) +
                " out_i[L/R/B/T]=" + fmt4u(wallOutMacroIBySideStep);
            if (detailHigh)
            {
                wallLine += " unk_e=" + std::to_string(wallOutMacroEUnknownStep) +
                            " unk_i=" + std::to_string(wallOutMacroIUnknownStep) +
                            " qout_e[C]=" + fmt4d(wallOutChargeEBySideStep_C, 2) +
                            " qout_i[C]=" + fmt4d(wallOutChargeIBySideStep_C, 2);
            }
            sections.emplace_back("Wall", wallLine);
        }

        if (detailHigh && physicsShowMCC)
        {
            if (mccEnabled)
            {
                sections.emplace_back(
                    "MCC",
                    "att=" + std::to_string(mccStepStats.n_attempts) +
                        " coll=" + std::to_string(mccStepStats.n_collisions) +
                        " null=" + std::to_string(mccStepStats.n_null) +
                        " spawn=" + std::to_string(mccStepStats.n_spawned) +
                        " byChan(el/ex/ion/att/cx/other)=" +
                        std::to_string(nElastic) + '/' + std::to_string(nExc) + '/' +
                        std::to_string(nIoniz) + '/' + std::to_string(nAttach) + '/' +
                        std::to_string(nCx) + '/' + std::to_string(nOther) +
                        " maxNuDt=" + fmtCompact(mccStepStats.max_nu_dt_majorant, 3));
            }
            else
            {
                sections.emplace_back("MCC", "off");
            }
        }

        if (showEmbedded)
        {
            std::ostringstream embeddedLine;
            embeddedLine << "absorb=" << embeddedLossMacroTotal
                         << " Qsink_emb[C]=" << fmtCompact(stepCharge.embeddedC, 3)
                         << " Eabs_emb[J]=" << fmtCompact(dE_abs_embedded_J, 3);
            if (!embeddedLossByBoundaryIdStep.empty())
            {
                std::vector<std::pair<int, std::uint64_t>> byId(
                    embeddedLossByBoundaryIdStep.begin(), embeddedLossByBoundaryIdStep.end());
                std::sort(byId.begin(), byId.end(), [](const auto& a, const auto& b) {
                    if (a.second != b.second)
                        return a.second > b.second;
                    return a.first < b.first;
                });
                embeddedLine << " byId(";
                for (std::size_t i = 0; i < byId.size(); ++i)
                {
                    if (i)
                        embeddedLine << ',';
                    embeddedLine << byId[i].first << ':' << byId[i].second;
                }
                embeddedLine << ')';
            }
            sections.emplace_back("Embedded", embeddedLine.str());
        }

        if (showMccSplit && mccEnabled && !mccStepStats.per_projectile.empty())
        {
            auto fmtSplit = [&](const char* name, const MccSplit& a) {
                std::ostringstream oss;
                oss << name
                    << "(att/coll/null/spawn)=" << a.att << '/' << a.coll << '/' << a.nul << '/' << a.spawn
                    << " byChan=" << a.el << '/' << a.ex << '/' << a.ion << '/' << a.at << '/' << a.cx << '/' << a.other
                    << " maxNuDt=" << fmtCompact(a.maxNuDt, 3);
                return oss.str();
            };

            std::ostringstream splitLine;
            splitLine << fmtSplit("e", mccE) << " | " << fmtSplit("ion", mccI) << " | " << fmtSplit("neu", mccN);
            if (!mccStepStats.spawn_per_species.empty())
            {
                std::vector<std::pair<std::string, std::uint64_t>> sp;
                sp.reserve(mccStepStats.spawn_per_species.size());
                for (const auto& kv : mccStepStats.spawn_per_species)
                    sp.push_back(kv);
                std::sort(sp.begin(), sp.end(), [](const auto& a, const auto& b) {
                    return a.second > b.second;
                });
                splitLine << " spawnBySpecies=";
                for (std::size_t i = 0; i < std::min<std::size_t>(4, sp.size()); ++i)
                {
                    if (i)
                        splitLine << ',';
                    splitLine << sp[i].first << ':' << sp[i].second;
                }
            }

            sections.emplace_back("MCCSplit", splitLine.str());
        }

        if (detailHigh && physicsShowPerf)
        {
            sections.emplace_back(
                "Perf",
                "total_s=" + fmtCompact(stepTiming.total_s, 3) +
                    " pm_s=" + fmtCompact(stepTiming.pm_s, 3) +
                    " top=" + top1Name + ":" + fmtCompact(top1S, 2) +
                    "|" + top2Name + ":" + fmtCompact(top2S, 2) +
                    "|" + top3Name + ":" + fmtCompact(top3S, 2) +
                    " rG=" + fmtCompact(diagLast.rG_L2_per_active, 3) +
                    " fallbackFrac=" + fmtCompact(totalFallbackFrac, 2) +
                    " amr=" + amrFlag +
                    " out=" + outFlag);
        }

        if (summaryOneLine)
        {
            output << "[ImplicitPIC] ";
            for (std::size_t i = 0; i < sections.size(); ++i)
            {
                if (i)
                    output << " | ";
                output << sections[i].first << ' ' << sections[i].second;
            }
            output << std::endl;
        }
        else
        {
            const std::string sep = "--------------------------------";
            output << "[ImplicitPIC] " << sep << std::endl;
            for (const auto& sec : sections)
                output << "[ImplicitPIC][" << sec.first << "] " << sec.second << std::endl;
            output << "[ImplicitPIC] " << sep << std::endl;
        }
    }

    return nP;
}

void PicStepSeries::openPhysics(const std::filesystem::path& path, int resumeStep)
{
    open(physics_, path, resumeStep,
        "step,time_s,rG_L2,rG_Linf,dQrel_stage,dQphys_rel,"
        "Qsource_injected_C,Qsink_domain_C,Qsink_embedded_C,Qmacro_begin_C,Qmacro_end_C,"
        "field_energy_J,kinetic_e_J,kinetic_ion_J,kinetic_neu_J,kinetic_total_J,total_energy_J,"
        "J_cond_Apm2,J_disp_Apm2,J_total_Apm2,macro_total,phys_total\n", "physics");
}

void PicStepSeries::writePhysics(int step, double timeS, const ImplicitOverview& row,
    const PicStepParticleBalance::Result& charge)
{
    const auto& v = row.values;
    physics_ << step << ',' << timeS << ',' << row.residualL2 << ',' << row.residualLinf << ','
             << charge.stageRelative << ',' << charge.relative << ',' << charge.sourceC << ','
             << charge.domainC << ',' << charge.embeddedC << ',' << charge.macroBeginC << ','
             << charge.macroEndC << ',' << v.fieldEnergyJ << ',' << v.kineticElectronJ << ','
             << v.kineticIonJ << ',' << v.kineticNeutralJ << ',' << v.kineticTotalJ << ','
             << v.totalEnergyJ << ',' << v.conductionCurrent << ',' << v.displacementCurrent << ','
             << v.totalCurrent << ',' << v.macroTotal << ',' << v.physicalTotal << "\n";
}

void PicStepSeries::openAmount(const std::filesystem::path& path, int resumeStep)
{
    open(amount_, path, resumeStep,
        "step,time_s,Qgrid_stage_C,Qmacro_stage_C,dQ_stage_C,dQrel_stage\n", "amount");
}
void PicStepSeries::openAmountSpecies(const std::filesystem::path& path, int resumeStep)
{
    open(amountSpecies_, path, resumeStep,
        "step,time_s,species,N_amount_stage,N_phys_particles,dN_rel_stage\n", "amount species");
}
void PicStepSeries::writeAmount(int step, double timeS, const PicStepParticleBalance::Result& charge)
{
    amount_ << step << ',' << timeS << ',' << charge.stageGridC << ',' << charge.stageMacroC << ','
            << charge.stageResidualC << ',' << charge.stageRelative << "\n";
}
void PicStepSeries::writeAmountSpecies(int step, double timeS, std::string_view species,
    const PicSpeciesAmount& row)
{
    amountSpecies_ << step << ',' << timeS << ',' << species << ',' << row.stage << ','
                   << row.physical << ',' << row.relative << "\n";
}
void PicStepSeries::openEmbeddedBand(const std::filesystem::path& path, int resumeStep)
{
    if (mode_ == Mode::Explicit)
    {
        open(embeddedBand_, path, resumeStep,
            "step,time_s,band_layers,cut_cells,band_cells,E_rms_Vpm,E_mean_Vpm\n", "EB-band RMS");
        embeddedBand_ << std::defaultfloat;
        return;
    }
    const auto parent = path.parent_path();
    if (!parent.empty()) std::filesystem::create_directories(parent);
    // This optional stream historically fails silently when it cannot open.
    open(embeddedBand_, path, resumeStep,
        "step,time_s,band_layers,n_cells,E_rms_Vpm,E_mean_Vpm\n", nullptr);
}
void PicStepSeries::writeEmbeddedBand(int step, double timeS, const PicEmbeddedBandField& row)
{
    embeddedBand_ << step << ',' << timeS << ',' << row.layers << ',';
    if (mode_ == Mode::Explicit) embeddedBand_ << row.cutCells << ',';
    embeddedBand_ << row.cells << ','
                  << row.rmsVpm << ',' << row.meanVpm << "\n";
}
void PicStepSeries::flushAmounts()
{
    if (amountReady()) amount_.flush();
    if (amountSpeciesReady()) amountSpecies_.flush();
}

PicExplicitConsole::PicExplicitConsole(const PicConfig& cfg, int firstStep, int endStep)
    : startStep(firstStep), steps(endStep)
{
    physicsSummaryEvery = std::max(1, cfg.output.console.every);
    physicsSummaryEnabled = cfg.output.console.physicsEnabled;
    physicsShowCharge = cfg.output.console.physicsShowCharge;
    physicsShowCurrent = cfg.output.console.physicsShowCurrent;
    physicsShowWall = cfg.output.console.physicsShowWall;
    physicsShowMCC = cfg.output.console.physicsShowMCC;
    physicsShowPerf = cfg.output.console.physicsShowPerf;
    physicsDetail = PhysicsDiagnostics::lowerCopy(cfg.output.console.physicsDetail);
    summaryOneLine = (PhysicsDiagnostics::lowerCopy(cfg.output.v2.console.format) != "multiline");
}

void PicExplicitConsole::writeStep(int stepForPrint, double time_code, double dt, double eta_s,
    const PicStepSeries::ExplicitPerformance& timing, const Sample& sample, std::ostream& output) const
{
    if (!ConsoleControl::enabled(ConsoleControl::Level::kInfo))
        return;
    if (!physicsSummaryEnabled)
        return;
    if (!ConsoleControl::shouldPrintPeriodic(stepForPrint, physicsSummaryEvery))
        return;

    const std::string detail =
        (physicsDetail == "minimal" || physicsDetail == "medium" || physicsDetail == "high")
            ? physicsDetail
            : std::string("medium");
    const bool detailMinimal = (detail == "minimal");
    const bool detailHigh = (detail == "high");
    const bool showParticles = !detailMinimal;
    const bool showWall = !detailMinimal && physicsShowWall;
    const bool showCurrent = physicsShowCurrent;
    const bool showSee = detailHigh && (sample.events.seeEnabled || sample.emission.impactMacro != 0 || sample.emission.emittedMacro != 0);
    const bool showEmbedded = detailHigh && (sample.events.embeddedReady ||
                                             sample.physical.embeddedStatistics || sample.events.embeddedResidual != 0);
    const bool showMccSplit = detailHigh && physicsShowMCC;

    const double t_si = (sample.physical.units != nullptr) ? sample.physical.units->TimeCodeToSI(time_code) : time_code;
    const double dt_si = (sample.physical.units != nullptr) ? sample.physical.units->TimeCodeToSI(dt) : dt;

    const auto values = summarizeExplicitPicConsole(sample.physical, dt_si);
    const auto& [phiMin_SI, phiMax_SI, eMeanE_eV, TeMean_eV, eTailFrac,
        n_elastic, n_exc, n_ioniz, n_attach, n_cx, n_other, mccE, mccI, mccN] = values.state;
    constexpr double tailThr_eV = PicConsoleValues::tailThresholdEV;
    using MccSplit = PicConsoleValues::MccSplit;
    const auto& species = sample.physical.particles.species();
    const auto& electrodeGroups = sample.electrodes.groups();
    const auto& consoleElectrodeGroupBuckets = sample.electrodes.buckets();
    const auto* ebStats = sample.physical.embeddedStatistics ? &sample.physical.losses.embeddedCounts : nullptr;

    auto fmtSci = [](double v, int p = 3) {
        std::ostringstream oss;
        oss.setf(std::ios::scientific);
        oss << std::setprecision(p) << v;
        return oss.str();
    };
    auto fmtFix = [](double v, int p = 3) {
        std::ostringstream oss;
        oss.setf(std::ios::fixed);
        oss << std::setprecision(p) << v;
        return oss.str();
    };
    auto fmtDef = [](double v, int p = 6) {
        std::ostringstream oss;
        oss << std::setprecision(p) << v;
        return oss.str();
    };
    auto fmtMaybeSci = [&](double v, int p = 3) {
        return std::isfinite(v) ? fmtSci(v, p) : std::string("nan");
    };

    std::ostringstream particleSummary;
    for (std::size_t si = 0; si < species.size(); ++si)
    {
        if (si) particleSummary << ' ';
        particleSummary << species[si].name() << '=' << values.species[si].macro
                        << '(' << fmtSci(values.species[si].physical) << ')';
    }

    auto fmt4 = [](const std::array<uint64_t, 4>& a) {
        std::ostringstream oss;
        oss << a[0] << ',' << a[1] << ',' << a[2] << ',' << a[3];
        return oss.str();
    };
    auto fmtSidesSci = [&](const std::array<double, 4>& a) {
        return std::string("L=") + fmtMaybeSci(a[0]) +
               " R=" + fmtMaybeSci(a[1]) +
               " B=" + fmtMaybeSci(a[2]) +
               " T=" + fmtMaybeSci(a[3]);
    };
    auto fmtEventTag = [&](bool amrAdapted, bool didOutput) {
        std::string evt;
        if (amrAdapted)
            evt += "amr-adapt";
        if (didOutput)
        {
            if (!evt.empty())
                evt += ',';
            evt += "output";
        }
        if (evt.empty())
            evt = "advance";
        return evt;
    };
    auto classifyStage = [&](int stepNumber, bool amrAdapted) {
        const int stepsSinceStart = std::max(0, stepNumber - startStep);
        if (stepsSinceStart <= 1)
            return std::string("field-init");
        if (stepsSinceStart <= 10)
            return std::string("startup");
        const bool quasiSteady =
            !amrAdapted && (sample.diagnostic != nullptr) &&
            std::abs(values.chargeRelative) <= 1e-3 &&
            std::abs(sample.diagnostic->rG_Linf_rel_rhs) <= 1e-3;
        return quasiSteady ? std::string("quasi-steady") : std::string("transient");
    };
    auto fmtGroupCircuit = [&](double dtSeconds) {
        std::ostringstream oss;
        oss << "Igrp_A{";
        for (std::size_t gi = 0; gi < electrodeGroups.size(); ++gi)
        {
            if (gi)
                oss << ' ';
            const auto& bucket = consoleElectrodeGroupBuckets[gi];
            const double invDt = (dtSeconds > 0.0) ? (1.0 / dtSeconds) : 0.0;
            oss << electrodeGroups[gi].name << '=' << fmtMaybeSci(bucket.chargeC * invDt);
            if (detailHigh)
            {
                oss << "(e=" << fmtMaybeSci(bucket.electronChargeC * invDt)
                    << " i=" << fmtMaybeSci(bucket.ionChargeC * invDt) << ')';
            }
        }
        oss << '}';
        return oss.str();
    };

    std::array<std::pair<std::string, double>, 8> perfStages = {{
        {"buildWeights", timing.field.buildWeights_s},
        {"deposit", timing.field.deposit_s},
        {"rhs", timing.field.rhs_s},
        {"poisson", timing.poissonSeconds[3]},
        {"evalE", timing.field.evalE_s},
        {"amr", timing.amrAdaptSeconds},
        {"pm", timing.managementSeconds},
        {"gatherPush", timing.gatherPushSeconds}
    }};
    std::sort(perfStages.begin(), perfStages.end(), [](const auto& a, const auto& b) {
        return a.second > b.second;
    });
    std::ostringstream perfTop;
    for (std::size_t i = 0; i < 3 && i < perfStages.size(); ++i)
    {
        if (i)
            perfTop << '|';
        perfTop << perfStages[i].first << ':' << fmtMaybeSci(perfStages[i].second);
    }

    std::vector<std::pair<std::string, std::string>> sections;
    sections.emplace_back(
        "Progress",
        "step=" + std::to_string(stepForPrint) + "/" + std::to_string(steps) +
            " t_s=" + fmtSci(t_si) +
            " dt_s=" + fmtSci(dt_si) +
            " np=" + std::to_string(values.macroCount) +
            " stage=" + classifyStage(stepForPrint, sample.events.adapted) +
            " evt=" + fmtEventTag(sample.events.adapted, sample.events.outputWritten) +
            " ETA=" + formatPicEta(eta_s));

    if (showParticles)
        sections.emplace_back("Particles", particleSummary.str());

    {
        std::string stateLine;
        if (physicsShowCharge)
        {
            stateLine += "phi[V]={[" + fmtDef(phiMin_SI, 6) + "," + fmtDef(phiMax_SI, 6) + "]}";
            stateLine += " ";
        }
        stateLine +=
            "eMeanE_eV=" + fmtFix(eMeanE_eV) +
            " Te_eV=" + fmtFix(TeMean_eV) +
            " tail_e_gt_" + fmtFix(tailThr_eV, 3) + "eV=" + fmtFix(eTailFrac);
        sections.emplace_back("State", stateLine);
    }

    if (showCurrent)
    {
        if (!electrodeGroups.empty())
        {
            sections.emplace_back(
                "Circuit",
                fmtGroupCircuit(dt_si) +
                    " Jtot_Apm2=" + fmtMaybeSci(values.totalCurrent));
        }
        else
        {
            sections.emplace_back(
                "Circuit",
                "Jtot_Apm2=" + fmtMaybeSci(values.totalCurrent) +
                    " Je_Apm2{" + fmtSidesSci(values.electronCurrent) + "}" +
                    " Ji_Apm2{" + fmtSidesSci(values.ionCurrent) + "}");
        }
    }

    {
        auto fmtMccCompact = [&]() {
            if (!sample.events.collisionsEnabled)
                return std::string("mcc=off");
            const std::array<std::pair<const char*, std::uint64_t>, 6> counts = {{
                {"el", n_elastic},
                {"ex", n_exc},
                {"ion", n_ioniz},
                {"att", n_attach},
                {"cx", n_cx},
                {"other", n_other},
            }};
            std::ostringstream oss;
            oss << "mcc{";
            bool first = true;
            for (const auto& entry : counts)
            {
                if (entry.second == 0)
                    continue;
                if (!first)
                    oss << ' ';
                oss << entry.first << '=' << entry.second;
                first = false;
            }
            if (first)
                oss << '0';
            oss << '}';
            return oss.str();
        };
        std::string sourceSinkLine =
            "N[s/m/d/e]=" + std::to_string(sample.events.sourceSpawn) + '/' +
            std::to_string(sample.physical.collisions.n_spawned) + '/' +
            std::to_string(values.wallLossMacro) + '/' +
            std::to_string(values.embeddedAbsorbed) +
            " Q[s/d/e]=" + fmtSci(sample.events.sourceCharge) + '/' +
            fmtSci(sample.events.domainCharge) + '/' +
            fmtSci(sample.events.embeddedCharge);
        if (physicsShowMCC)
            sourceSinkLine += ' ' + fmtMccCompact();
        if (sample.events.seeEnabled || sample.emission.emittedMacro != 0 || sample.emission.impactMacro != 0)
            sourceSinkLine += " seeEmit=" + std::to_string(sample.emission.emittedMacro) +
                              " ionImpCath=" + std::to_string(sample.emission.impactMacro);
        sections.emplace_back("SourceSink", sourceSinkLine);
    }

    {
        std::string healthLine =
            "dQ_C=" + fmtSci(values.chargeDifference) +
            " dQ_rel=" + fmtSci(values.chargeRelative);
        if ((sample.diagnostic != nullptr))
        {
            healthLine += " rG_L2_rel_rhs=" + fmtSci(sample.diagnostic->rG_L2_rel_rhs) +
                          " rG_Linf_rel_rhs=" + fmtSci(sample.diagnostic->rG_Linf_rel_rhs) +
                          " S_eb_rel_rhs=" + fmtSci(sample.diagnostic->S_eb_rel_rhs);
        }
        else
        {
            healthLine += " diag=pending";
        }
        sections.emplace_back("Health", healthLine);
    }

    if (showWall)
    {
        sections.emplace_back(
            "Wall",
            "out_e[L/R/B/T]=" + fmt4(sample.physical.losses.domainMacroE) +
                " out_i[L/R/B/T]=" + fmt4(sample.physical.losses.domainMacroI));
    }

    if (detailHigh && physicsShowMCC)
    {
        if (sample.events.collisionsEnabled)
        {
            sections.emplace_back(
                "MCC",
                "att=" + std::to_string(sample.physical.collisions.n_attempts) +
                    " coll=" + std::to_string(sample.physical.collisions.n_collisions) +
                    " null=" + std::to_string(sample.physical.collisions.n_null) +
                    " spawn=" + std::to_string(sample.physical.collisions.n_spawned) +
                    " byChan(el/ex/ion/att/cx/other)=" +
                    std::to_string(n_elastic) + '/' + std::to_string(n_exc) + '/' +
                    std::to_string(n_ioniz) + '/' + std::to_string(n_attach) + '/' +
                    std::to_string(n_cx) + '/' + std::to_string(n_other) +
                    " maxNuDt=" + fmtSci(sample.physical.collisions.max_nu_dt_majorant));
        }
        else
        {
            sections.emplace_back("MCC", "off");
        }
    }

    if (showSee)
    {
        sections.emplace_back(
            "SEE",
            "ionImp_cath=" + std::to_string(sample.emission.impactMacro) +
                "(" + fmtSci(sample.emission.impactPhysical) + ")" +
                " emit=" + std::to_string(sample.emission.emittedMacro) +
                "(" + fmtSci(sample.emission.emittedPhysical) + ")" +
                " eEmitMean_eV=" + fmtFix(sample.emission.meanEV) +
                " eEmitStd_eV=" + fmtFix(sample.emission.stdEV));
    }

    if (showEmbedded)
    {
        std::ostringstream embeddedLine;
        const std::uint64_t total = values.embeddedAbsorbed;
        embeddedLine << "absorb=" << total;
        if (ebStats != nullptr)
        {
            embeddedLine << " inside=" << ebStats->inside
                         << " crossing=" << ebStats->crossing
                         << " unk=" << ebStats->unknown;
            if (!ebStats->byId.empty())
            {
                embeddedLine << " byId(";
                bool first = true;
                for (const auto& kv : ebStats->byId)
                {
                    if (!first)
                        embeddedLine << ',';
                    first = false;
                    embeddedLine << kv.first << ':' << kv.second;
                }
                embeddedLine << ')';
            }
        }
        embeddedLine << " residual=" << sample.events.embeddedResidual;
        sections.emplace_back("Embedded", embeddedLine.str());
    }

    if (showMccSplit && sample.events.collisionsEnabled && !sample.physical.collisions.per_projectile.empty())
    {
        auto fmtSplit = [&](const char* name, const MccSplit& a) {
            std::ostringstream oss;
            oss << name
                << "(att/coll/null/spawn)=" << a.att << '/' << a.coll << '/' << a.nul << '/' << a.spawn
                << " byChan=" << a.el << '/' << a.ex << '/' << a.ion << '/' << a.at << '/' << a.cx << '/' << a.other
                << " maxNuDt=" << fmtSci(a.maxNuDt);
            return oss.str();
        };

        std::ostringstream splitLine;
        splitLine << fmtSplit("e", mccE)
                  << " | " << fmtSplit("ion", mccI)
                  << " | " << fmtSplit("neu", mccN);

        if (!sample.physical.collisions.spawn_per_species.empty())
        {
            std::vector<std::pair<std::string, uint64_t>> sp;
            sp.reserve(sample.physical.collisions.spawn_per_species.size());
            for (const auto& kv : sample.physical.collisions.spawn_per_species)
                sp.push_back(kv);
            std::sort(sp.begin(), sp.end(), [](const auto& a, const auto& b) { return a.second > b.second; });
            splitLine << " spawnBySpecies=";
            for (std::size_t i = 0; i < std::min<std::size_t>(4, sp.size()); ++i)
            {
                if (i)
                    splitLine << ',';
                splitLine << sp[i].first << ':' << sp[i].second;
            }
        }

        sections.emplace_back("MCCSplit", splitLine.str());
    }

    if (detailHigh && physicsShowPerf)
    {
        sections.emplace_back(
            "Perf",
            "total_s=" + fmtMaybeSci(timing.totalSeconds) +
                " top=" + perfTop.str() +
                " rG=off fallbackFrac=off" +
                " amr=" + std::string(sample.events.adapted ? "adapt" : "none") +
                " out=" + std::string(sample.events.outputWritten ? "out" : "-"));
    }

    if (summaryOneLine)
    {
        output << "[ExplicitPIC] ";
        for (std::size_t i = 0; i < sections.size(); ++i)
        {
            if (i)
                output << " | ";
            output << sections[i].first << ' ' << sections[i].second;
        }
        output << std::endl;
    }
    else
    {
        const std::string sep = "--------------------------------";
        output << "[ExplicitPIC] " << sep << std::endl;
        for (const auto& sec : sections)
            output << "[ExplicitPIC][" << sec.first << "] " << sec.second << std::endl;
        output << "[ExplicitPIC] " << sep << std::endl;
    }
}
