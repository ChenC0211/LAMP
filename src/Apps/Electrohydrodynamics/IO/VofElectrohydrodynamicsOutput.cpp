#include "Electrohydrodynamics/IO/VofElectrohydrodynamicsOutput.h"
#include "Electrohydrodynamics/VofElectrohydrodynamicsRun.h"
#include "IO/HydroOutputComponents.h"
#include "Hydrodynamics/HydrodynamicsRun.h"
#include "Hydrodynamics/HydrodynamicsVofRegrid.h"
#include "Hydrodynamics/HydrodynamicsTimeStep.h"
#include "Modules/Hydrodynamics/HydroVofMeshTransfer.h"
#include "Modules/Electric/ElectricVofChargeTransport.h"
#include "Core/IO/OutputProvenance.h"
#include "Core/IO/TecplotWriter.h"
#include "Core/IO/BinaryArchive.h"
#include "Core/Geometry/Detail/FiniteVolumeClosure.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cctype>
#include <utility>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <locale>
#include <map>
#include <sstream>
#include <stdexcept>
#include <system_error>
#include <thread>
#include <optional>
#include <limits>
#include <cfenv>
#include <cerrno>
#include <exception>
#if defined(_M_X64) || defined(__SSE2__)
#include <xmmintrin.h>
#endif
#ifdef _WIN32
#include <winerror.h>
#endif

namespace {
namespace fs=std::filesystem;
namespace H=LAMP::Hydrodynamics;
namespace E=LAMP::Electric;
namespace I=LAMP::Interface;
namespace Fluid=AppIO::HydroOutput;
using EhdClock=LAMP::TimeLayer;
using LAMP::IO::stream;
using LAMP::IO::quoted;
constexpr const char* app="VofElectrohydrodynamics";
struct PassivePerformanceEnvironment {
    std::fenv_t native{};int error=errno;
#if defined(_M_X64) || defined(__SSE2__)
    unsigned csr=_mm_getcsr();
#endif
    PassivePerformanceEnvironment(){if(std::fegetenv(&native)!=0){errno=error;throw std::runtime_error("performance serialization environment capture failed");}}
    ~PassivePerformanceEnvironment() noexcept {
        if(std::fesetenv(&native)!=0)std::terminate();
#if defined(_M_X64) || defined(__SSE2__)
        _mm_setcsr(csr);
#endif
        errno=error;
    }
};
void require(bool ok,const char* why){if(!ok)throw std::invalid_argument(std::string("[VofEHD output] ")+why);}
bool sameClock(EhdClock a,EhdClock b){return a.step==b.step&&a.begin==b.begin&&a.end==b.end&&a.position==b.position;}
template<class T> std::string text(T x){auto s=stream();s<<x;return s.str();}
void write(const fs::path& path,const std::string& contents){
    std::ofstream f;f.exceptions(std::ios::failbit|std::ios::badbit);
    f.open(path,std::ios::binary|std::ios::trunc);f<<contents;f.close();
}
void replace(const fs::path& path,const std::string& contents){
    const fs::path temporary=path.string()+".tmp";write(temporary,contents);
    std::error_code ec;fs::rename(temporary,path,ec);
    constexpr std::array<int,6> delays{{5,10,20,40,80,160}};
    for(std::size_t retry=0;ec&&retry<delays.size();++retry){
#ifdef _WIN32
        const bool sharing=ec.category()==std::system_category()&&(ec.value()==ERROR_ACCESS_DENIED||ec.value()==ERROR_SHARING_VIOLATION||ec.value()==ERROR_LOCK_VIOLATION);
#else
        const bool sharing=false;
#endif
        if(!sharing)break;std::this_thread::sleep_for(std::chrono::milliseconds(delays[retry]));fs::rename(temporary,path,ec);
    }
    if(ec)throw std::runtime_error("VofEHD output index "+path.generic_string()+": "+ec.message());
}
std::vector<TecplotWriter::TableColumn> schema(const std::string& names){
    std::istringstream in(names);std::vector<TecplotWriter::TableColumn> result;std::string name;
    while(std::getline(in,name,',')){while(!name.empty()&&(name.back()=='\n'||name.back()=='\r'))name.pop_back();
        const bool category=name=="stage"||name=="status"||name=="reason"||name=="event"||name=="field_kind"
            ||name=="transfer_rule"||name=="evaluation";
        result.push_back({name,category?TecplotWriter::TableColumnKind::Category:TecplotWriter::TableColumnKind::Numeric});}
    return result;
}
// Serialization schemas copied from src/Apps/Hydrodynamics/IO/HydrodynamicsOutput.cpp.
// These match the public HydroOutput vofRows/vofPerformance column order;
// this writer borrows actual rows and performs no fluid or electric calculation.
const std::string conservationHeader=R"SCHEMA(step,time,dt,geometry_identity,negative_phase_volume,initial_volume,relative_volume_drift,total_geometric_mass,divergence_max,dual_R_max,native_max_speed,kappa_max,phase_volume_closure,mass_closure,liquid_p0_moment_x,liquid_p0_moment_y,liquid_p0_moment_xx,liquid_p0_moment_yy,liquid_p0_moment_xy,vof_work_kinetic_change,vof_work_capillary,vof_work_body,vof_work_viscous,vof_work_pressure,vof_work_half_R,vof_work_half_R_mapped_primal,vof_work_half_R_geometry_transfer,vof_work_half_R_spatial_exchange,vof_work_time_mass,vof_work_rounded_phase_mass,vof_work_operator_rate_roundoff,vof_work_raw_minus_canonical_outflow_half,vof_work_half_R_decomposition_residual,vof_work_momentum_equation,vof_work_kinetic_budget_residual,vof_work_closure_residual,vof_work_closure_relative_residual,vof_stable_gas_duals,vof_stable_gas_half_R_work,vof_stable_gas_spatial_exchange_half_R_work,vof_stable_gas_max_abs_integrated_R_over_mid_mass,vof_peak_gas_native_face,vof_peak_gas_abs_end_normal_velocity,vof_peak_gas_integrated_R_over_mid_mass,vof_peak_gas_mass_change_over_mid_mass,vof_peak_gas_outward_integral_over_mid_mass,vof_peak_gas_geometry_transfer_over_mid_mass,vof_peak_gas_spatial_exchange_over_mid_mass)SCHEMA";
const std::string solverHeader=R"SCHEMA(step,time,stage,status,AC_stages,linear_iterations,linear_refinements,algebraic_residual,momentum_residual,reason
)SCHEMA";
const std::string energyHeader=R"SCHEMA(step,time,kinetic,surface,total,body_power,cumulative_body_work,cumulative_viscous_loss,work_corrected_total,capillary_energy_defect,interface_power,interface_power_available,capillary_power,capillary_power_available
)SCHEMA";
const std::string performanceHeader=R"SCHEMA(step,time,stage,wall_seconds,cells,native_velocity_dofs,pressure_dofs,concurrent_viscosity_capillary_stage_count
)SCHEMA";
const std::string capillaryHeader=R"SCHEMA(step,accepted_time,force_time,force_geometry_identity,force_kappa_min,force_kappa_max,loaded_dofs,height_cells,height_columns,height_column_cell_queries,height_alternate_directions,height_extended_cells,height_extended_slices,height_refined_cells,height_refined_slices,height_fit_cells,height_fit_segments,unresolved_isolated_fragment_cells
)SCHEMA";
const std::string reconstructionTemplatesHeader="step,time,entries,row_queries,row_hits,row_builds,port_queries,port_hits,port_builds,actual_svd_calls,retained_bytes";
std::string reconstructionTemplateRow(const H::HydroTwoPhaseState& fluid){
    const auto t=fluid.timeLayer();const auto templates=fluid.momentumGeometry().regularTemplateStatistics();
    auto row=stream();row<<t.step<<','<<t.time()<<','<<templates.entries<<','
        <<templates.rowQueries<<','<<templates.rowHits<<','<<templates.rowBuilds<<','
        <<templates.portQueries<<','<<templates.portHits<<','<<templates.portBuilds<<','
        <<templates.actualSvdCalls<<','<<templates.retainedBytes<<'\n';
    return row.str();
}
const std::string eventsHeader=R"SCHEMA(event,attempted_step,accepted_step,accepted_time,old_dt,new_dt,stage,reason
)SCHEMA";
// Only the retained final-AC transport invocation, never an all-AC ledger.
const std::string localPlicExecutionHeader="step,time,begin,end,position,physical_end_geometry_identity,accepted_geometry_identity,ac_stage,source_geometry_identity,source_step,source_begin,source_end,source_position,midpoint_geometry_identity,endpoint_geometry_identity,requested,eligible,face_values_parallel_used,cell_gather_parallel_used,fallback,native_faces,incidences,directions_used,face_evaluations,cell_gathers,scratch_bytes,scratch_limit_bytes,maximum_workers";
std::string localPlicExecutionRow(const H::HydroTwoPhaseState& accepted,const H::HydroVofStep* step){
    if(!step)return {};
    PassivePerformanceEnvironment environment;
    const auto t=accepted.timeLayer();const auto& actual=step->transport();
    const auto source=actual.sourceTime();const auto& r=actual.report().localFaceExecution;
    auto row=stream();row<<t.step<<','<<t.time()<<','<<t.begin<<','<<t.end<<','<<static_cast<int>(t.position)<<','
        <<step->candidateState().phaseGeometry().vofIdentity().value()<<','<<accepted.phaseGeometry().vofIdentity().value()<<','
        <<step->report().stages.size()<<','<<actual.sourceIdentity().value()<<','<<source.step<<','<<source.begin<<','<<source.end<<','
        <<static_cast<int>(source.position)<<','<<actual.midpoint().identity().value()<<','<<actual.endpoint().identity().value()<<','
        <<r.requested<<','<<r.eligible<<','<<r.faceValuesParallelUsed<<','<<r.cellGatherParallelUsed<<','<<r.fallback;
    if(r.eligible)row<<','<<r.nativeFaces<<','<<r.incidences<<','<<r.directionsUsed<<','<<r.faceEvaluations<<','<<r.cellGathers<<','<<r.scratchBytes;
    else row<<",NA,NA,NA,NA,NA,NA";
    row<<','<<r.scratchLimitBytes<<',';if(r.directionsUsed)row<<r.maximumWorkers;else row<<"NA";
    row<<'\n';return row.str();
}
const std::string timeStepHeader=R"SCHEMA(step,time,begin,actual_dt,proposed_dt,flow_rate,capillary_rate,interface_h,flow_dt_bound,capillary_dt_bound,growth_dt_bound,max_dt,reason,retries,endpoint_aligned,final_below_min,minimum_clock_rounded
)SCHEMA";
const std::string stressGeometryReuseHeader="step,time,begin,end,position,accepted_geometry_identity,requested,enabled,target_cells,reused_raw_cells,reused_stress_cells,recomputed_fit_cells,reused_faces,rebuilt_faces,stencil_rejected_cells,quadrature_rejected_cells,column_rejected_cells,direct_rebound_rows,literal_rebound_rows,csr_rebind_scratch_bytes,raw_direct_rebound_rows,volume_direct_rebound_rows";
std::string stressGeometryReuseRow(const H::HydroTwoPhaseState& accepted,const H::HydroVofStep* step){
    if(!step||!step->frozenMidpoint().forces().maxwellLoad)return {};
    PassivePerformanceEnvironment environment;
    const auto t=accepted.timeLayer();
    const auto& report=step->frozenMidpoint().forces().maxwellLoad->geometry().geometryReuseReport();
    auto row=stream();row<<t.step<<','<<t.time()<<','<<t.begin<<','<<t.end<<','<<static_cast<int>(t.position)<<','
        <<accepted.phaseGeometry().vofIdentity().value()<<','<<report.requested<<','<<report.enabled<<','
        <<report.targetCells<<','<<report.reusedRawCells<<','<<report.reusedStressCells<<','<<report.recomputedFitCells<<','
        <<report.reusedFaces<<','<<report.rebuiltFaces<<','<<report.stencilRejectedCells<<','
        <<report.quadratureRejectedCells<<','<<report.columnRejectedCells<<','
        <<report.directReboundRows<<','<<report.literalReboundRows<<','<<report.csrRebindScratchBytes<<','
        <<report.rawDirectReboundRows<<','<<report.volumeDirectReboundRows<<'\n';
    return row.str();
}
const std::string solverStagesHeader="step,time,begin,end,physical_end_position,physical_end_geometry_identity,accepted_geometry_identity,ac_stage,final_ac,used_block_gmres,used_direct_fallback,direct_fallback_reason,used_superlu,krylov_iterations,krylov_restarts,krylov_solves,krylov_stagnated,linear_refinements,velocity_block_seconds,preconditioner_seconds,factor_seconds,direct_factor_seconds,solve_seconds,algebraic_relative_residual,free_momentum_defect,max_cell_continuity_defect,max_pressure_mean_defect,max_boundary_flux_defect,max_abs_gauge_multiplier,krylov_relative_residual,krylov_scaled_relative_residual,factor_dimension,factor_nonzeros,free_velocity_count,pressure_gauge_count,free_tangential_count,velocity_spd_tier,velocity_spd_gram_attempted,velocity_spd_gram_bounds_available,velocity_spd_gram_bound_upper,velocity_spd_scaled_mass_lower,velocity_spd_certificate_seconds,velocity_spd_llt_seconds,frozen_preparation_timings_available,saddle_assembly_timings_available,frozen_input_validation_seconds,frozen_force_validation_seconds,frozen_matrix_action_seconds,frozen_rhs_seconds,frozen_velocity_prepare_seconds,saddle_partition_seconds,saddle_system_assembly_seconds,saddle_quotient_assembly_seconds,saddle_rhs_elimination_seconds,velocity_sparse_plan_requested,velocity_sparse_plan_eligible,velocity_sparse_plan_hit,velocity_sparse_exact_symmetry_used,velocity_sparse_exact_skew_used,velocity_sparse_nonzeros,velocity_sparse_retained_bytes,saddle_scatter_requested,saddle_scatter_eligible,saddle_scatter_plan_hit,saddle_system_scatter_used,saddle_quotient_scatter_used,saddle_system_nonzeros,saddle_quotient_nonzeros,saddle_scatter_retained_bytes,velocity_split_plan_requested,velocity_split_plan_eligible,velocity_split_plan_hit,velocity_split_symmetric_values_used,velocity_split_skew_values_used,velocity_split_input_nonzeros,velocity_split_union_nonzeros,velocity_split_retained_bytes,velocity_append_plan_used,velocity_append_scratch_bytes,ordered_spmv_requested,ordered_spmv_eligible,ordered_spmv_used,ordered_spmv_fallback,ordered_spmv_calls,ordered_spmv_parallel_calls,ordered_spmv_view_bytes,ordered_spmv_maximum_workers,direct_input_deferred,direct_input_materialized,parallel_pressure_block_count,parallel_pressure_maximum_workers,pressure_coarse_prepared,pressure_coarse_dimension,pressure_coarse_reason,pressure_coarse_applications,pressure_coarse_condition_upper,pressure_coarse_formation_norm_upper,pressure_global_requested,pressure_global_available,pressure_global_prepared,pressure_global_timings_available,pressure_global_rows,pressure_global_nonzeros,pressure_global_factor_nonzeros,pressure_global_storage_bytes,pressure_global_admission_bytes,pressure_global_byte_limit,pressure_global_symbolic_seconds,pressure_global_numeric_seconds,pressure_global_solve_calls,pressure_global_solve_seconds,velocity_split_parallel_used,velocity_split_parallel_fallback,velocity_split_maximum_workers,saddle_symbolic_built,saddle_symbolic_workers,saddle_scatter_parallel_used,saddle_scatter_parallel_fallback,saddle_scatter_maximum_workers,saddle_quotient_aliases_system,solve_breakdown_timings_available,krylov_precondition_apply_seconds,krylov_spmv_seconds,krylov_orthogonalization_seconds,krylov_true_residual_seconds,physical_assessment_seconds,orthogonalization_parallel_calls,orthogonalization_maximum_workers";
std::string solverStageRows(const H::HydroTwoPhaseState& accepted,const H::HydroVofStep* step){
    if(!step)return {};
    PassivePerformanceEnvironment environment;
    const auto& physical=step->candidateState();const auto t=physical.timeLayer();const auto& stages=step->report().stages;
    auto rows=stream();
    for(std::size_t ac=0;ac<stages.size();++ac){const auto& s=stages[ac].saddle;
        rows<<t.step<<','<<t.time()<<','<<t.begin<<','<<t.end<<','<<static_cast<int>(t.position)<<','
            <<physical.phaseGeometry().vofIdentity().value()<<','<<accepted.phaseGeometry().vofIdentity().value()<<','
            <<ac+1<<','<<(ac+1==stages.size()?1:0)<<','<<s.usedBlockGmres<<','<<s.usedDirectFallback<<','
            <<s.directFallbackReason<<','<<s.usedSuperLU<<','<<s.krylovIterations<<','<<s.krylovRestarts<<','
            <<s.krylovSolves<<','<<s.krylovStagnated<<','<<s.linearRefinements<<','<<s.velocityBlockSeconds<<','
            <<s.preconditionerSeconds<<','<<s.factorSeconds<<','<<s.directFactorSeconds<<','<<s.solveSeconds<<','
            <<s.algebraicRelativeResidual<<','<<s.freeMomentumDefect<<','<<s.maxCellContinuityDefect<<','
            <<s.maxPressureMeanDefect<<','<<s.maxBoundaryFluxDefect<<','<<s.maxAbsGaugeMultiplier<<','
            <<s.krylovRelativeResidual<<','<<s.krylovScaledRelativeResidual<<','<<s.factorDimension<<','
            <<s.factorNonZeros<<','<<s.freeVelocityCount<<','<<s.pressureGaugeCount<<','<<s.freeTangentialCount<<','
            <<s.velocitySpdTier<<','<<s.velocitySpdGramAttempted<<','<<s.velocitySpdGramBoundsAvailable<<',';
        if(s.velocitySpdGramBoundsAvailable)rows<<s.velocitySpdGramBoundUpper<<','<<s.velocitySpdScaledMassLower;
        else rows<<"NA,NA";
        rows<<','<<s.velocitySpdCertificateSeconds<<','<<s.velocitySpdLltSeconds<<','
            <<s.frozenPreparationTimingsAvailable<<','<<s.saddleAssemblyTimingsAvailable;
        const std::array<double,5> frozenTimes{{s.frozenInputValidationSeconds,s.frozenForceValidationSeconds,
            s.frozenMatrixActionSeconds,s.frozenRhsSeconds,s.frozenVelocityPrepareSeconds}};
        const std::array<double,4> saddleTimes{{s.saddlePartitionSeconds,s.saddleSystemAssemblySeconds,
            s.saddleQuotientAssemblySeconds,s.saddleRhsEliminationSeconds}};
        for(double value:frozenTimes){rows<<',';if(s.frozenPreparationTimingsAvailable&&std::isfinite(value))rows<<value;else rows<<"NA";}
        for(double value:saddleTimes){rows<<',';if(s.saddleAssemblyTimingsAvailable&&std::isfinite(value))rows<<value;else rows<<"NA";}
        const auto& p=s.velocitySparseValidation;
        rows<<','<<p.requested<<','<<p.eligible<<','<<p.planHit<<','<<p.exactSymmetryUsed<<','<<p.exactSkewUsed<<',';
        if(p.requested)rows<<p.nonZeros;else rows<<"NA";
        rows<<',';if(p.retainedBytes)rows<<p.retainedBytes;else rows<<"NA";
        const auto& a=s.sparseAssemblyPlan;
        rows<<','<<a.requested<<','<<a.eligible<<','<<a.planHit<<','<<a.systemScatterUsed<<','<<a.quotientScatterUsed<<',';
        if(a.requested)rows<<a.systemNonZeros<<','<<a.quotientNonZeros;else rows<<"NA,NA";
        rows<<',';if(a.retainedBytes)rows<<a.retainedBytes;else rows<<"NA";
        const auto& split=s.velocitySplitPlan;
        rows<<','<<split.requested<<','<<split.eligible<<','<<split.planHit<<','
            <<split.symmetricValuesUsed<<','<<split.skewValuesUsed<<',';
        if(split.requested)rows<<split.inputNonZeros<<','<<split.unionNonZeros;else rows<<"NA,NA";
        rows<<',';if(split.retainedBytes)rows<<split.retainedBytes;else rows<<"NA";
        rows<<','<<p.appendPlanUsed<<',';
        if(p.appendPlanUsed)rows<<p.appendScratchBytes;else rows<<"NA";
        rows<<','<<s.orderedSpmvRequested<<','<<s.orderedSpmvEligible<<','
            <<s.orderedSpmvUsed<<','<<s.orderedSpmvFallback<<','
            <<s.orderedSpmvCalls<<','<<s.orderedSpmvParallelCalls<<',';
        if(s.orderedSpmvEligible)rows<<s.orderedSpmvViewBytes;else rows<<"NA";
        rows<<',';if(s.orderedSpmvUsed)rows<<s.orderedSpmvMaximumWorkers;else rows<<"NA";
        rows<<','<<s.directInputDeferred<<','<<s.directInputMaterialized;
        rows<<',';if(s.parallelPressureBlockCount>0)rows<<s.parallelPressureBlockCount;else rows<<"NA";
        rows<<',';if(s.parallelPressureBlockCount>0)rows<<s.parallelPressureMaximumWorkers;else rows<<"NA";
        if(s.parallelPressureBlockCount>0){
            rows<<','<<s.pressureCoarsePrepared<<','<<s.pressureCoarseDimension<<','<<s.pressureCoarseReason<<','<<s.pressureCoarseApplications;
            if(s.pressureCoarsePrepared)rows<<','<<s.pressureCoarseConditionUpper<<','<<s.pressureCoarseFormationNormUpper;
            else rows<<",NA,NA";
        }else rows<<",NA,NA,NA,NA,NA,NA";
        rows<<','<<s.pressureGlobalRequested<<','<<s.pressureGlobalAvailable<<','<<s.pressureGlobalPrepared<<','<<s.pressureGlobalTimingsAvailable;
        if(s.pressureGlobalAvailable)rows<<','<<s.pressureGlobalRows<<','<<s.pressureGlobalNonZeros<<','<<s.pressureGlobalFactorNonZeros
            <<','<<s.pressureGlobalStorageBytes<<','<<s.pressureGlobalAdmissionBytes<<','<<s.pressureGlobalByteLimit;
        else rows<<",NA,NA,NA,NA,NA,NA";
        if(s.pressureGlobalTimingsAvailable)rows<<','<<s.pressureGlobalSymbolicSeconds<<','<<s.pressureGlobalNumericSeconds;
        else rows<<",NA,NA";
        rows<<',';if(s.pressureGlobalAvailable)rows<<s.pressureGlobalSolveCalls;else rows<<"NA";
        rows<<',';if(s.pressureGlobalTimingsAvailable)rows<<s.pressureGlobalSolveSeconds;else rows<<"NA";
        rows<<','<<split.parallelUsed<<','<<split.parallelFallback<<',';
        if(split.parallelUsed)rows<<split.maximumWorkers;else rows<<"NA";
        rows<<','<<a.symbolicBuilt<<',';
        if(a.symbolicBuilt&&a.symbolicWorkers>0)rows<<a.symbolicWorkers;else rows<<"NA";
        rows<<','<<a.scatterParallelUsed<<','<<a.scatterParallelFallback<<',';
        if(a.scatterParallelUsed)rows<<a.scatterMaximumWorkers;else rows<<"NA";
        rows<<','<<a.quotientAliasesSystem;
        rows<<','<<s.solveBreakdownTimingsAvailable;
        const std::array<double,5> solveTimes{{s.krylovPreconditionApplySeconds,s.krylovSpmvSeconds,
            s.krylovOrthogonalizationSeconds,s.krylovTrueResidualSeconds,s.physicalAssessmentSeconds}};
        for(double value:solveTimes){rows<<',';if(s.solveBreakdownTimingsAvailable&&std::isfinite(value))rows<<value;else rows<<"NA";}
        rows<<',';if(s.usedBlockGmres)rows<<s.orthogonalizationParallelCalls;else rows<<"NA";
        rows<<',';if(s.usedBlockGmres)rows<<s.orthogonalizationMaximumWorkers;else rows<<"NA";
        rows<<'\n';
    }
    return rows.str();
}
const std::string csgBoundaryPreparationHeader="step,time,physical_geometry_identity,accepted_geometry_identity,scope,requests,hits,rebuilds,retained_builds,selector_calls,complete_graph_builds,peak_retained_bytes,byte_limit,timing_available,total_seconds,rebuild_seconds";
std::string csgBoundaryPreparationRows(const H::HydroTwoPhaseState& accepted,const H::HydroVofStep* step,
    const VofEhdCsgBoundaryReport* report){
    if(!report)return {};
    PassivePerformanceEnvironment environment;
    const auto t=accepted.timeLayer();const auto acceptedIdentity=accepted.phaseGeometry().vofIdentity().value();
    const auto physicalIdentity=step?step->candidateState().phaseGeometry().vofIdentity().value():acceptedIdentity;
    auto rows=stream();
    for(std::size_t scope=0;scope<report->stages.size();++scope){const auto& r=report->stages[scope];
        if(r.requests==0)continue;
        rows<<t.step<<','<<t.time()<<','<<physicalIdentity<<','<<acceptedIdentity<<','<<scope<<','
            <<r.requests<<','<<r.hits<<','<<r.rebuilds<<','<<r.retainedBuilds<<','<<r.selectorCalls<<','<<r.completeGraphBuilds
            <<','<<r.peakRetainedBytes<<','<<report->byteLimit<<','<<r.timingAvailable;
        if(r.timingAvailable)rows<<','<<r.totalSeconds<<','<<r.rebuildSeconds;else rows<<",NA,NA";
        rows<<'\n';
    }
    return rows.str();
}
const std::string viscosityProductStagesHeader="step,time,begin,end,physical_end_position,physical_end_geometry_identity,accepted_geometry_identity,ac_stage,final_ac,weight_source_geometry_identity,weight_source_step,weight_source_begin,weight_source_end,weight_source_position,requested,eligible,source_available,source_fp_receipt,current_fp_qualified,local_used,full_assembly,memo_published,fallback_reason,weight_rows,changed_weight_rows,output_rows,recomputed_rows,reused_rows,symbolic_stored_values,retained_bytes,peak_optional_bytes,byte_limit,timing_available,symbolic_built,symbolic_workers,column_routing_seconds,symbolic_graph_seconds,numeric_assembly_seconds,no_slip_preparation_available,no_slip_raw_cache_hit,no_slip_cold_built,no_slip_parallel_used,no_slip_parallel_fallback,no_slip_fit_workers,no_slip_timing_available,no_slip_cold_seconds,no_slip_fit_requests,no_slip_fit_hits,no_slip_qr_factorizations,no_slip_svd_factorizations,no_slip_bank_peak_bytes,no_slip_bank_limit_bytes";
std::string viscosityProductStageRows(const H::HydroTwoPhaseState& accepted,const H::HydroVofStep* step){
    if(!step)return {};
    PassivePerformanceEnvironment environment;
    const auto& physical=step->candidateState();const auto t=physical.timeLayer();const auto& stages=step->report().stages;
    auto rows=stream();
    for(std::size_t ac=0;ac<stages.size();++ac){
        rows<<t.step<<','<<t.time()<<','<<t.begin<<','<<t.end<<','<<static_cast<int>(t.position)<<','
            <<physical.phaseGeometry().vofIdentity().value()<<','<<accepted.phaseGeometry().vofIdentity().value()<<','
            <<ac+1<<','<<(ac+1==stages.size()?1:0);
        if(!stages[ac].viscosityProductReportAvailable){for(unsigned k=0;k<29;++k)rows<<",NA";
            rows<<",0";for(unsigned k=0;k<13;++k)rows<<",NA";rows<<'\n';continue;}
        const auto& wt=stages[ac].viscosityProductTime;
        rows<<','<<stages[ac].viscosityProductPhaseIdentity<<','<<wt.step<<','<<wt.begin<<','<<wt.end<<','<<static_cast<int>(wt.position);
        const auto& r=stages[ac].viscosityProduct;
        rows<<','<<r.requested<<','<<r.eligible<<','<<r.sourceAvailable<<','<<r.sourceFpReceipt<<','<<r.currentFpQualified
            <<','<<r.localUsed<<','<<r.fullAssembly<<','<<r.memoPublished<<','<<r.fallbackReason;
        if(r.requested)rows<<','<<r.weightRows;else rows<<",NA";
        if(r.eligible)rows<<','<<r.changedWeightRows<<','<<r.outputRows<<','<<r.recomputedRows<<','<<r.reusedRows;
        else rows<<",NA,NA,NA,NA";
        if(r.requested)rows<<','<<r.symbolicStoredValues;else rows<<",NA";
        if(r.retainedBytes)rows<<','<<r.retainedBytes;else rows<<",NA";
        if(r.peakOptionalBytes)rows<<','<<r.peakOptionalBytes;else rows<<",NA";
        rows<<','<<r.byteLimit<<','<<r.timingAvailable<<','<<r.symbolicBuilt<<','<<r.symbolicWorkers;
        if(r.timingAvailable)rows<<','<<r.columnRoutingSeconds<<','<<r.symbolicGraphSeconds<<','<<r.numericAssemblySeconds;
        else rows<<",NA,NA,NA";
        const auto& ns=r.noSlip;rows<<','<<ns.available;
        if(ns.available){rows<<','<<ns.rawCacheHit<<','<<ns.coldBuilt<<','<<ns.parallelUsed<<','<<ns.parallelFallback
                <<','<<ns.maximumWorkers<<','<<ns.timingAvailable<<',';
            if(ns.timingAvailable)rows<<ns.coldSeconds;else rows<<"NA";
            rows<<','<<ns.fitRequests<<','<<ns.fitHits<<','<<ns.qrFactorizations<<','<<ns.svdFactorizations
                <<','<<ns.bankPeakBytes<<','<<ns.bankLimitBytes;
        }else for(unsigned k=0;k<13;++k)rows<<",NA";
        rows<<'\n';
    }
    return rows.str();
}
const std::string maxwellStagesHeader="step,time,begin,end,physical_end_position,physical_end_geometry_identity,accepted_geometry_identity,ac_stage,final_ac,compiled_contraction,compiled_geometry_hit,compiled_timing_available,compiled_volume_terms,compiled_boundary_terms,compiled_retained_bytes,validation_seconds,compiled_plan_seconds,sampling_and_contraction_seconds,final_rne_seconds";
std::string maxwellStageRows(const H::HydroTwoPhaseState& accepted,const H::HydroVofStep* step){
    if(!step)return {};
    PassivePerformanceEnvironment environment;
    const auto& physical=step->candidateState();const auto t=physical.timeLayer();const auto& stages=step->report().stages;
    auto rows=stream();
    for(std::size_t ac=0;ac<stages.size();++ac){
        if(!stages[ac].maxwellLoadReportAvailable)continue;
        const auto& r=stages[ac].maxwellLoadReport;
        rows<<t.step<<','<<t.time()<<','<<t.begin<<','<<t.end<<','<<static_cast<int>(t.position)<<','
            <<physical.phaseGeometry().vofIdentity().value()<<','<<accepted.phaseGeometry().vofIdentity().value()<<','
            <<ac+1<<','<<(ac+1==stages.size()?1:0)<<','<<r.compiledContraction<<','<<r.compiledGeometryHit<<','<<r.compiledTimingAvailable<<',';
        if(r.compiledContraction)rows<<r.compiledVolumeTerms<<','<<r.compiledBoundaryTerms<<','<<r.compiledRetainedBytes;
        else rows<<"NA,NA,NA";
        const std::array<double,4> times{{r.validationSeconds,r.compiledPlanSeconds,r.samplingAndContractionSeconds,r.finalRneSeconds}};
        for(double value:times){rows<<',';if(r.compiledTimingAvailable&&std::isfinite(value))rows<<value;else rows<<"NA";}
        rows<<'\n';
    }
    return rows.str();
}
// Every actual AC contributes BOTH evaluations. Times below nest inside the
// existing convection/diagnostics parents; never add them again to step cost.
const std::string openConvectionStagesHeader="step,time,begin,end,position,physical_end_geometry_identity,accepted_geometry_identity,ac_stage,final_accepted_ac,evaluation,material_begin_geometry_identity,material_midpoint_geometry_identity,material_end_geometry_identity,material_step,material_begin,material_end,material_position,native_rows,segments,pairs,boundary_binding_seconds,operator_seconds,advector_boundary_seconds,actual_call_seconds,validation_seconds,preparation_seconds,segment_seconds,diagonal_seconds,pair_seconds,row_seconds,ledger_seconds,failed_optional_seconds,total_seconds,parallel_requested,parallel_eligible,parallel_used,parallel_fallback,contribution_workers,gather_workers,scratch_bytes";
std::string openConvectionStageRows(const H::HydroTwoPhaseState& accepted,const H::HydroVofStep* step)
{
    if(!step)return {};
    PassivePerformanceEnvironment environment;
    const auto& candidate=step->candidateState();const auto t=candidate.timeLayer();const auto& stages=step->report().stages;auto rows=stream();
    for(std::size_t ac=0;ac<stages.size();++ac){const auto& c=stages[ac].openConvection;if(!c.available)continue;
        const auto append=[&](const char* evaluation,const H::HydroKepOpenApplyReport& r,double callSeconds,bool first){
            require(r.available&&r.timeLayer.step==t.step&&r.timeLayer.begin==t.begin&&r.timeLayer.end==t.end
                &&r.timeLayer.position==LAMP::TimeLayer::Position::Midpoint,"open convection actual apply full clock differs");
            rows<<t.step<<','<<t.time()<<','<<t.begin<<','<<t.end<<','<<static_cast<int>(t.position)<<','
                <<candidate.phaseGeometry().vofIdentity().value()<<','<<accepted.phaseGeometry().vofIdentity().value()<<','<<ac+1<<','<<(ac+1==stages.size())<<','<<evaluation<<','
                <<r.vofSources[0]<<','<<r.vofSources[1]<<','<<r.vofSources[2]<<','<<r.timeLayer.step<<','<<r.timeLayer.begin<<','<<r.timeLayer.end<<','<<static_cast<int>(r.timeLayer.position)<<','
                <<r.nativeRows<<','<<r.segments<<','<<r.pairs<<',';
            if(first)rows<<c.boundaryBindingSeconds<<','<<c.operatorSeconds<<','<<c.advectorBoundarySeconds;else rows<<"NA,NA,NA";
            const std::array<double,10> times{{callSeconds,r.validationSeconds,r.preparationSeconds,r.segmentSeconds,r.diagonalSeconds,r.pairSeconds,r.rowSeconds,r.ledgerSeconds,r.fallbackSeconds,r.totalSeconds}};
            for(double value:times){rows<<',';if(std::isfinite(value)&&value>=0)rows<<value;else rows<<"NA";}
            rows<<','<<r.parallelRequested<<','<<r.parallelEligible<<','<<r.parallelUsed<<','<<r.parallelFallback<<',';
            if(r.parallelUsed)rows<<r.contributionWorkers<<','<<r.gatherWorkers<<','<<r.scratchBytes;else rows<<"NA,NA,NA";
            rows<<'\n';
        };
        require(c.advector.vofSources==c.solvedMidpoint.vofSources&&sameClock(c.advector.timeLayer,c.solvedMidpoint.timeLayer),
            "open convection evaluations differ in actual source/clock");
        append("advector",c.advector,c.advectorApplySeconds,true);
        append("solved_midpoint",c.solvedMidpoint,c.solvedMidpointApplySeconds,false);
    }
    return rows.str();
}
const std::string electricPreparationHeader="step,time,accepted_geometry_identity,stage,field_begin,field_end,field_position,field_geometry_identity,source_charge_identity,candidate_charge_identity,geometry_seconds,assemble_seconds,factor_seconds,solve_seconds,total_seconds,advection_stage_seconds,charge_pair_midpoint_advance_seconds,charge_pair_endpoint_advance_seconds,used_superlu,electrostatic_snapshot,field_postprocess_parallel_used,field_postprocess_parallel_fallback,field_postprocess_workers,field_postprocess_scratch_bytes,used_gmres,linear_iterations,linear_restarts,linear_maximum_workers,field_postprocess_timing_available,field_postprocess_face_seconds,field_postprocess_cell_seconds,field_postprocess_finalize_seconds,field_postprocess_seconds,solver_expanded_residual_available,solver_expanded_residual_parallel_used,solver_expanded_residual_parallel_fallback,solver_expanded_residual_workers,solver_expanded_residual_rows,solver_expanded_residual_scratch_bytes,solver_expanded_residual_timing_available,solver_expanded_residual_route_seconds,solver_expanded_residual_seconds";
std::string electricPreparationRows(const H::HydroTwoPhaseState& accepted,const H::HydroVofStep* step,
    const E::ElectricVofChargeTransport* physicalTransport,const E::ElectricVofField* physicalEnd){
    if(!step)return {};
    PassivePerformanceEnvironment environment;
    const auto t=accepted.timeLayer();auto rows=stream();
    const auto append=[&](const char* stage,const E::ElectricVofField& field,const E::ElectricVofChargeTransport* pair){
        const auto c=field.candidateTime();const auto& r=field.report();
        const double unavailable=std::numeric_limits<double>::quiet_NaN();
        const auto* advection=field.advection();
        const std::array<double,8> times{{r.geometrySeconds,r.assembleSeconds,r.factorSeconds,r.solveSeconds,r.totalSeconds,
            !pair&&advection?advection->report().totalSeconds:unavailable,
            pair?pair->midpoint().report().totalSeconds:unavailable,
            pair?pair->endpoint().report().totalSeconds:unavailable}};
        rows<<t.step<<','<<t.time()<<','<<accepted.phaseGeometry().vofIdentity().value()<<','<<stage<<','
            <<c.begin<<','<<c.end<<','<<static_cast<int>(c.position)<<','<<field.geometry().identity().value()<<','
            <<field.sourceCharge().identity()<<','<<field.chargeCandidate().identity();
        for(double value:times){rows<<',';if(std::isfinite(value))rows<<value;else rows<<"NA";}
        rows<<','<<r.usedSuperLU<<','<<r.electrostaticSnapshot<<','
            <<r.postprocessParallelUsed<<','<<r.postprocessParallelFallback<<',';
        if(r.postprocessParallelUsed)rows<<r.postprocessWorkers<<','<<r.postprocessScratchBytes;else rows<<"NA,NA";
        rows<<','<<r.usedGMRES<<',';
        if(r.usedGMRES)rows<<r.linearIterations<<','<<r.linearRestarts<<','<<r.linearMaximumWorkers;else rows<<"NA,NA,NA";
        const auto tailTime=[&](bool available,double value){rows<<',';if(available&&std::isfinite(value))rows<<value;else rows<<"NA";};
        rows<<','<<r.postprocessTimingAvailable;
        for(double value:{r.postprocessFaceSeconds,r.postprocessCellSeconds,r.postprocessFinalizeSeconds,r.postprocessSeconds})
            tailTime(r.postprocessTimingAvailable,value);
        rows<<','<<r.integratedExpandedResidualAvailable;
        if(r.integratedExpandedResidualAvailable)rows<<','<<r.integratedExpandedResidualParallelUsed<<','<<r.integratedExpandedResidualParallelFallback<<','
            <<r.integratedExpandedResidualWorkers<<','<<r.integratedExpandedResidualRows<<','<<r.integratedExpandedResidualScratchBytes<<','<<r.integratedExpandedResidualTimingAvailable;
        else rows<<",NA,NA,NA,NA,NA,NA";
        tailTime(r.integratedExpandedResidualAvailable&&r.integratedExpandedResidualTimingAvailable,r.integratedExpandedResidualRouteSeconds);
        tailTime(r.integratedExpandedResidualAvailable&&r.integratedExpandedResidualTimingAvailable,r.integratedExpandedResidualSeconds);
        rows<<'\n';
    };
    const auto& force=step->frozenMidpoint().forces().maxwellLoad;
    if(force&&!force->isSharp())append("final_ac_midpoint",force->stress().field(),nullptr);
    if(physicalEnd)append("physical_end",*physicalEnd,physicalTransport);
    return rows.str();
}
// Actual authenticated target construction only. The VOF identities below
// describe this remap's physical/accepted provenance, not a Regular build ID.
const std::string frozenMatrixPreparationHeader="step,time,begin,end,physical_end_position,physical_end_geometry_identity,accepted_geometry_identity,stage,material_midpoint_geometry_identity,material_midpoint_step,material_midpoint_begin,material_midpoint_end,material_midpoint_position,actual_velocity_dofs,requested,eligible,physical_plan_hit,physical_values_used,union_plan_hit,union_values_used,viscosity_nonzeros,convection_nonzeros,physical_nonzeros,velocity_nonzeros,retained_bytes,physical_cold_compiled,physical_cold_used,union_cold_compiled,union_cold_used,values_parallel_used,values_parallel_fallback,values_maximum_workers,values_scratch_bytes";
std::string frozenMatrixPreparationRow(const H::HydroTwoPhaseState& accepted,const H::HydroVofStep* step){
    if(!step)return {};
    PassivePerformanceEnvironment environment;
    const auto& frozen=step->frozenMidpoint();const auto& material=frozen.material();
    const auto physical=step->candidateState().timeLayer(),actual=material.timeLayer();
    const auto& r=frozen.matrixAssemblyPlanReport();auto row=stream();
    row<<physical.step<<','<<physical.time()<<','<<physical.begin<<','<<physical.end<<','<<static_cast<int>(physical.position)<<','
        <<step->candidateState().phaseGeometry().vofIdentity().value()<<','<<accepted.phaseGeometry().vofIdentity().value()<<",final_ac_midpoint,"
        <<material.sources().vofIds[1].value()<<','<<actual.step<<','<<actual.begin<<','<<actual.end<<','<<static_cast<int>(actual.position)<<','
        <<material.layout().velocityCount()<<','<<r.requested<<','<<r.eligible<<','<<r.physicalPlanHit<<','<<r.physicalValuesUsed<<','
        <<r.unionPlanHit<<','<<r.unionValuesUsed<<','<<r.viscosityNonZeros<<','<<r.convectionNonZeros<<','<<r.physicalNonZeros<<','<<r.velocityNonZeros<<',';
    if(r.retainedBytes)row<<r.retainedBytes;else row<<"NA";
    row<<','<<r.physicalColdCompiled<<','<<r.physicalColdUsed<<','<<r.unionColdCompiled<<','<<r.unionColdUsed
        <<','<<r.valuesParallelUsed<<','<<r.valuesParallelFallback<<',';
    if(r.valuesParallelUsed)row<<r.valuesMaximumWorkers;else row<<"NA";
    row<<',';if(r.valuesParallelUsed||r.valuesParallelFallback)row<<r.valuesScratchBytes;else row<<"NA";
    row<<'\n';return row.str();
}
const std::string regularPreparationReuseHeader="step,time,begin,end,position,physical_source_geometry_identity,accepted_geometry_identity,source_cells,target_cells,requested,enabled,target_datums,reused_rows,reused_ports,recomputed_rows,recomputed_ports,stencil_rejected_rows,port_span_rejected_rows,column_rejected_rows,fp_rejected_rows,source_metadata_bytes";
std::string regularPreparationReuseRow(const H::HydroTwoPhaseState& accepted,
    const HydrodynamicsVofRemapReceipt* receipt){
    if(!receipt||!receipt->interfaceTransfer)return {};
    PassivePerformanceEnvironment environment;
    const auto r=accepted.momentumGeometry().regularPreparationReuseReport();
    if(!r.requested)return {};
    const auto t=accepted.timeLayer();const auto& tree=*receipt->decision;
    auto row=stream();row<<t.step<<','<<t.time()<<','<<t.begin<<','<<t.end<<','<<static_cast<int>(t.position)<<','
        <<receipt->physicalEnd->phaseGeometry().vofIdentity().value()<<','<<accepted.phaseGeometry().vofIdentity().value()<<','
        <<tree.sourceLeaves<<','<<tree.targetLeaves<<','<<r.requested<<','<<r.enabled<<','<<r.targetDatums<<','
        <<r.reusedRows<<','<<r.reusedPorts<<','<<r.recomputedRows<<','<<r.recomputedPorts<<','
        <<r.stencilRejectedRows<<','<<r.portSpanRejectedRows<<','<<r.columnRejectedRows<<','<<r.fpRejectedRows<<',';
    if(r.retainedMetadataBytes)row<<r.retainedMetadataBytes;else row<<"NA";
    row<<'\n';return row.str();
}
const std::string chargeTransportStagesHeader="step,time,begin,end,position,accepted_geometry_identity,stage,source_charge_identity,source_step,source_begin,source_end,source_position,transport_source_geometry_identity,stage_geometry_identity,stage_step,stage_begin,stage_end,stage_position,cells,faces,donors,nonzero_face_charges,total_seconds,donor_parallel_used,target_parallel_used,parallel_fallback,parallel_workers,parallel_scratch_bytes,limited_linear_fct,reconstructed_cells,high_order_faces,limited_faces,minimum_face_limiter,maximum_density_bound_defect";
std::string chargeTransportStageRows(const H::HydroTwoPhaseState& accepted,const H::HydroVofStep* step,
    const E::ElectricVofChargeTransport* physicalTransport){
    if(!step)return {};
    PassivePerformanceEnvironment environment;
    const auto t=step->candidateState().timeLayer();auto rows=stream();
    const auto append=[&](const char* stage,const E::ElectricVofAdvectedCharge& actual){
        const auto source=actual.sourceCharge().timeLayer(),candidate=actual.timeLayer();const auto& r=actual.report();
        rows<<t.step<<','<<t.time()<<','<<t.begin<<','<<t.end<<','<<static_cast<int>(t.position)<<','
            <<accepted.phaseGeometry().vofIdentity().value()<<','<<stage<<','<<actual.sourceCharge().identity()<<','
            <<source.step<<','<<source.begin<<','<<source.end<<','<<static_cast<int>(source.position)<<','
            <<actual.transport().source().identity().value()<<','<<actual.geometry().identity().value()<<','
            <<candidate.step<<','<<candidate.begin<<','<<candidate.end<<','<<static_cast<int>(candidate.position)<<','
            <<actual.integratedCharge().size()<<','<<actual.faces().size()<<','<<r.donors<<','<<r.nonzeroFaceCharges<<',';
        if(std::isfinite(r.totalSeconds))rows<<r.totalSeconds;else rows<<"NA";
        rows<<','<<r.donorParallelUsed<<','<<r.targetParallelUsed<<','<<r.parallelFallback<<',';
        if(r.donorParallelUsed||r.targetParallelUsed)rows<<r.parallelWorkers<<','<<r.parallelScratchBytes;
        else rows<<"NA,NA";
        rows<<','<<r.limitedLinearFct<<',';
        if(r.limitedLinearFct)rows<<r.reconstructedCells<<','<<r.highOrderFaces<<','<<r.limitedFaces<<','<<r.minimumFaceLimiter<<','<<r.maximumDensityBoundDefect;
        else rows<<"NA,NA,NA,NA,NA";
        rows<<'\n';
    };
    const auto& force=step->frozenMidpoint().forces().maxwellLoad;
    if(force&&!force->isSharp())if(const auto* actual=force->stress().field().advection())append("final_ac_midpoint",*actual);
    if(physicalTransport){append("physical_end_pair_midpoint",physicalTransport->midpoint());
        append("physical_end_pair_endpoint",physicalTransport->endpoint());}
    return rows.str();
}
const std::string stepCompleteHeader="step,time,begin,end,position,geometry_identity,actual_ac,wall_seconds,cells,native_velocity_dofs,pressure_dofs";
const std::string amrPerformanceHeader="step,time,begin,end,position,source_geometry_identity,target_geometry_identity,source_cells,target_cells,changed,tree_seconds,tree_indicator_seconds,tree_planning_seconds,tree_adaptation_seconds,inventory_seconds,fluid_remap_seconds,fluid_raw_flux_seconds,fluid_projection_seconds,charge_remap_seconds,target_snapshot_seconds,target_geometry_seconds,target_assemble_seconds,target_factor_seconds,target_solve_seconds,regrid_prepare_seconds,unattributed_seconds,target_stress_geometry_seconds,final_joint_validation_seconds,fluid_input_validation_seconds,fluid_target_layout_metric_seconds,fluid_target_momentum_phase_seconds,fluid_target_mass_seconds,fluid_source_mass_seconds,fluid_post_projection_validation_seconds,fluid_target_state_seconds,fluid_budget_validation_seconds,fluid_unattributed_seconds,tree_clone_seconds,tree_sweep_graph_seconds,tree_sweep_domain_seconds,tree_sweep_overlap_seconds,tree_native_adapt_seconds,tree_adaptation_unattributed_seconds,tree_clone_calls,tree_sweep_graph_calls,tree_sweep_domain_calls,tree_sweep_overlap_calls,tree_native_adapt_calls,tree_clone_logical_capacity_vertices,tree_clone_logical_capacity_horizontal_edges,tree_clone_logical_capacity_vertical_edges,tree_clone_logical_capacity_cells,tree_clone_source_constructed_vertices,tree_clone_source_constructed_horizontal_edges,tree_clone_source_constructed_vertical_edges,tree_clone_source_constructed_cells,tree_clone_target_constructed_vertices,tree_clone_target_constructed_horizontal_edges,tree_clone_target_constructed_vertical_edges,tree_clone_target_constructed_cells";
// Called only after the original full source/target/clock receipt validation.
// Child report timers are nested diagnostics, not additional physical work.
std::string amrPerformanceRow(const H::HydroTwoPhaseState& fluid,const E::ElectricVofField& field,
    const HydrodynamicsVofRemapReceipt* receipt,const E::ElectricVofChargeMeshTransfer* qremap){
    if(!receipt)return {};
    PassivePerformanceEnvironment environment;
    const auto& tree=*receipt->decision;const auto t=fluid.timeLayer();
    const bool changed=receipt->interfaceTransfer!=nullptr;
    const auto* hydro=changed?&receipt->fluidTransfer->report():nullptr;
    const auto* charge=changed?&qremap->report():nullptr;
    const auto* snapshot=changed?&field.report():nullptr;
    const double unavailable=std::numeric_limits<double>::quiet_NaN();
    const double inventory=changed?receipt->inventoryTransferSeconds:unavailable;
    double accounted=tree.totalSeconds;
    if(changed){accounted+=inventory;accounted+=hydro->totalSeconds;accounted+=charge->totalSeconds;accounted+=snapshot->totalSeconds;
        if(receipt->targetStressGeometrySeconds)accounted+=*receipt->targetStressGeometrySeconds;
        if(receipt->finalJointValidationSeconds)accounted+=*receipt->finalJointValidationSeconds;}
    const std::array<double,18> seconds{{tree.totalSeconds,tree.indicatorSeconds,tree.planningSeconds,tree.adaptationSeconds,
        inventory,hydro?hydro->totalSeconds:unavailable,hydro?hydro->rawFluxSeconds:unavailable,
        hydro?hydro->projectionSeconds:unavailable,charge?charge->totalSeconds:unavailable,
        snapshot?snapshot->totalSeconds:unavailable,snapshot?snapshot->geometrySeconds:unavailable,
        snapshot?snapshot->assembleSeconds:unavailable,snapshot?snapshot->factorSeconds:unavailable,
        snapshot?snapshot->solveSeconds:unavailable,receipt->regridSeconds,receipt->regridSeconds-accounted,
        changed?receipt->targetStressGeometrySeconds.value_or(unavailable):unavailable,
        changed?receipt->finalJointValidationSeconds.value_or(unavailable):unavailable}};
    const auto* fluidPreparation=hydro&&hydro->preparation.available?&hydro->preparation:nullptr;
    const std::array<double,9> fluidSeconds{{
        fluidPreparation?fluidPreparation->inputValidationSeconds:unavailable,
        fluidPreparation?fluidPreparation->targetLayoutMetricSeconds:unavailable,
        fluidPreparation?fluidPreparation->targetMomentumPhaseSeconds:unavailable,
        fluidPreparation?fluidPreparation->targetMassSeconds:unavailable,
        fluidPreparation?fluidPreparation->sourceMassSeconds:unavailable,
        fluidPreparation?fluidPreparation->postProjectionValidationSeconds:unavailable,
        fluidPreparation?fluidPreparation->targetStateSeconds:unavailable,
        fluidPreparation?fluidPreparation->budgetValidationSeconds:unavailable,
        fluidPreparation?fluidPreparation->unattributedSeconds:unavailable}};
    auto row=stream();row<<t.step<<','<<t.time()<<','<<t.begin<<','<<t.end<<','<<static_cast<int>(t.position)<<','
        <<receipt->physicalEnd->phaseGeometry().vofIdentity().value()<<','<<fluid.phaseGeometry().vofIdentity().value()<<','
        <<tree.sourceLeaves<<','<<tree.targetLeaves<<','<<(changed?1:0);
    for(double value:seconds){row<<',';if(std::isfinite(value))row<<value;else row<<"NA";}
    for(double value:fluidSeconds){row<<',';if(std::isfinite(value))row<<value;else row<<"NA";}
    const auto& treePreparation=tree.treePreparation;
    const std::array<double,6> treeSeconds{{treePreparation.cloneSeconds,treePreparation.graphSeconds,
        treePreparation.domainSeconds,treePreparation.overlapSeconds,treePreparation.nativeAdaptSeconds,treePreparation.unattributedSeconds}};
    const std::array<std::size_t,5> treeCalls{{treePreparation.cloneCalls,treePreparation.graphCalls,
        treePreparation.domainCalls,treePreparation.overlapCalls,treePreparation.nativeAdaptCalls}};
    for(double value:treeSeconds){row<<',';if(treePreparation.available&&std::isfinite(value))row<<value;else row<<"NA";}
    for(auto value:treeCalls){row<<',';if(treePreparation.available)row<<value;else row<<"NA";}
    for(const auto* stock:{&tree.cloneStorage.logicalCapacity,&tree.cloneStorage.sourceConstructed,&tree.cloneStorage.cloneConstructed})
        for(auto value:*stock){row<<',';if(tree.cloneStorage.available)row<<value;else row<<"NA";}
    row<<'\n';return row.str();
}
const I::GeometricVofMeshTransfer* validate(const QuadtreeMesh& mesh,const H::HydroTwoPhaseState& fluid,
    const E::ElectricVofCharge& charge,const E::ElectricVofField& field,const H::HydroVofStep* step,
    const HydrodynamicsVofRemapReceipt* receipt,const E::ElectricVofChargeMeshTransfer* qremap){
    fluid.validateFor(mesh);require(fluid.model()==H::HydroPhaseGeometry::Model::GeometricVof,"requires actual VOF state");
    const auto& geometry=fluid.phaseGeometry().vofPublication();const auto clock=fluid.timeLayer();
    require(clock.position==EhdClock::Position::Initial||clock.position==EhdClock::Position::End,"accepted field is not Initial/End");
    charge.validateFor(mesh,clock);
    field.validateFor(mesh,geometry,field.sourceCharge(),clock);
    require(field.geometry().identity()==geometry.identity()&&sameClock(field.geometryTime(),clock)&&sameClock(field.candidateTime(),clock)
        &&field.chargeCandidate().sameInstance(charge),"field/accepted geometry, charge owner or full clock differ");
    const bool changed=receipt&&receipt->interfaceTransfer;
    if(clock.position==EhdClock::Position::Initial||changed){
        require(field.isElectrostaticSnapshot()&&field.effectiveDt()==0&&field.sourceCharge().sameInstance(charge)
            &&sameClock(field.chargeSourceTime(),clock),"Initial/remapped fields require actual dt0 accepted snapshot");
    }else require(field.isElectrostaticSnapshot()||field.effectiveDt()>0,"End field lacks actual snapshot/Ohmic provenance");
    if(!receipt){require(!qremap,"charge remap lacks fluid receipt");
        require(!step||step->candidateState().sameIdentity(fluid),"physical step differs from accepted fluid");return nullptr;}
    require(step&&receipt->sourceMesh&&receipt->physicalEnd&&receipt->physicalStep==step&&receipt->decision,
        "incomplete original physical remap receipt");
    const auto& physical=*receipt->physicalEnd;physical.validateFor(*receipt->sourceMesh);
    require(physical.sameIdentity(step->candidateState())&&sameClock(physical.timeLayer(),clock)
        &&sameClock(receipt->decision->sourceClock,clock)&&receipt->decision->sourceLeaves==static_cast<std::size_t>(receipt->sourceMesh->cellIndexNum())
        &&receipt->decision->targetLeaves==static_cast<std::size_t>(mesh.cellIndexNum()),"physical source/target End receipt differs");
    require(std::isfinite(receipt->inventoryTransferSeconds)&&receipt->inventoryTransferSeconds>=0
        &&std::isfinite(receipt->regridSeconds)&&receipt->regridSeconds>=0,"invalid passive remap timings");
    require(bool(receipt->interfaceTransfer)==bool(receipt->fluidTransfer),"one fluid remap owner is missing");
    if(!changed){require(!receipt->decision->changed&&!qremap&&physical.sameIdentity(fluid)
            &&receipt->inventoryTransferSeconds==0,"no-op invented a remap");
        require(!receipt->targetStressGeometrySeconds&&!receipt->finalJointValidationSeconds,"no-op invented AMR child timings");return nullptr;}
    const auto& transfer=*receipt->interfaceTransfer;const auto& hydro=*receipt->fluidTransfer;
    require(receipt->decision->changed&&qremap&&hydro.sourceState().sameIdentity(physical)&&hydro.candidateState().sameIdentity(fluid)
        && &hydro.interfaceTransfer()==&transfer&&transfer.sourceIdentity()==physical.phaseGeometry().vofIdentity()
        &&transfer.candidate().identity()==geometry.identity()&&sameClock(transfer.candidate().timeLayer(),clock),
        "changed fluid source/target/transfer owner differs");
    qremap->validateFor(*receipt->sourceMesh,physical.phaseGeometry().vofPublication(),qremap->sourceCharge(),mesh,transfer);
    require(qremap->interfaceTransfer().sourceIdentity()==transfer.sourceIdentity()
        &&qremap->interfaceTransfer().candidate().identity()==transfer.candidate().identity()
        &&sameClock(qremap->interfaceTransfer().candidate().timeLayer(),clock)
        && &qremap->interfaceTransfer().nativeOverlaps()==&transfer.nativeOverlaps()
        &&qremap->chargeCandidate().sameInstance(charge)
        &&qremap->sourceGeometry().identity()==physical.phaseGeometry().vofIdentity()
        &&sameClock(qremap->sourceCharge().timeLayer(),clock),"charge remap physical-End/target owner differs");
    const auto validChild=[](const std::optional<double>& value){return !value||(std::isfinite(*value)&&*value>=0);};
    require(validChild(receipt->targetStressGeometrySeconds)&&validChild(receipt->finalJointValidationSeconds),"invalid passive AMR child timings");
    return &transfer;
}
// Output observations only: totals and dt times a rounded instantaneous J
// do not certify finite Ohmic transfer or become solver acceptance gates.
double chargeTotal(const std::vector<double>& values)
{
    PassivePerformanceEnvironment environment;
    LAMP::Geometry::Detail::FiniteVolume::Dyadic total;
    for(double value:values)total.add(value);
    const auto result=static_cast<double>(total.rounded());
    require(std::isfinite(result),"unrepresentable output charge total");return result;
}
double boundaryCurrent(const E::ElectricVofField& field)
{
    PassivePerformanceEnvironment environment;
    LAMP::Geometry::Detail::FiniteVolume::Dyadic total;
    for(const auto& face:field.faces())if((face.negativeCell>=0)!=(face.positiveCell>=0))
        total.add(face.current,face.negativeCell>=0?1:-1);
    const auto result=static_cast<double>(total.rounded());
    require(std::isfinite(result),"unrepresentable output boundary current");return result;
}
const std::string chargeBalanceHeader="step,time,begin,end,position,source_step,source_begin,source_end,source_position,source_geometry_identity,physical_end_geometry_identity,accepted_geometry_identity,source_charge_identity,physical_end_charge_identity,accepted_charge_identity,remapped,boundary_enabled,source_total_charge,advected_total_charge,physical_end_total_charge,accepted_total_charge,physical_delta_charge,advective_charge_in,advective_charge_out,advective_charge_net_out,advection_charge_residual,advection_relative_residual,advection_rounding_residual,inlet_product_rounding_residual,inlet_model_charge_residual,inlet_model_relative_residual,boundary_current,dt_boundary_current_from_report_total";
const std::string phaseBoundaryHeader="step,time,begin,end,position,source_geometry_identity,physical_end_geometry_identity,accepted_geometry_identity,source_liquid,target_liquid,inflow_liquid,outflow_liquid,WY_compression_liquid,rounding_liquid,physical_residual_liquid,algebraic_residual_liquid,source_gas,target_gas,inflow_gas,outflow_gas,WY_compression_gas,rounding_gas,physical_residual_gas,algebraic_residual_gas,physical_relative_residual,algebraic_relative_residual,mass_residual,mass_relative_residual";
const std::string boundaryWorkHeader="step,time,begin,end,position,source_geometry_identity,physical_end_geometry_identity,accepted_geometry_identity,open_flow,pressure_outlet_work,prescribed_reaction_work,boundary_kinetic_work,boundary_physical_kinetic_work,inlet_dissipation_work,boundary_momentum_rate_x_r,boundary_momentum_rate_y_z,kinetic_change,kinetic_budget_residual,closure_residual,closure_relative_residual,boundary_flux_rounding_work";
const std::string projectionWorkHeader="step,time,begin,end,position,source_geometry_identity,target_geometry_identity,projection_boundary_work,projection_rhs_rounding_work,projection_free_residual_work,projection_continuity_work,projection_kinetic_change,projection_correction_energy,projection_identity_residual,projection_identity_relative_residual,source_to_target_kinetic_change";
struct PhysicalRows {std::string charge,phase,work,projection;};
PhysicalRows physicalRows(const QuadtreeMesh& mesh,const H::HydroTwoPhaseState& accepted,
    const E::ElectricVofCharge& charge,const H::HydroVofStep* step,
    const HydrodynamicsVofRemapReceipt* remap,const E::ElectricVofChargeMeshTransfer* qremap,
    const E::ElectricVofChargeTransport* moved,const E::ElectricVofField* field,std::uint64_t previousCharge){
    require(bool(moved)==bool(field),"physical charge transport/End field receipts must be paired");
    const auto t=accepted.timeLayer();PhysicalRows rows;
    if(t.position==EhdClock::Position::Initial){require(!moved,"Initial output cannot contain a physical timestep receipt");return rows;}
    if(!moved)return rows;
    require(step&&t.position==EhdClock::Position::End,"physical receipt requires an actual End step");
    const auto& physicalMesh=remap?*remap->sourceMesh:mesh;
    const auto& physical=remap?*remap->physicalEnd:accepted;
    const auto& pub=physical.phaseGeometry().vofPublication();
    const auto& source=step->beginState();const auto s=source.timeLayer();
    require(physical.sameIdentity(step->candidateState())&&sameClock(physical.timeLayer(),t),"physical receipt has a different End state");
    require(moved->transport().sameInstance(step->transport())
        &&moved->transport().sourceIdentity()==source.phaseGeometry().vofIdentity()
        &&sameClock(moved->transport().sourceTime(),s)&&sameClock(moved->transport().endTime(),t),
        "physical charge transport is not the actual final fluid transport");
    moved->sourceCharge().validateFor(physicalMesh,s);
    require(moved->sourceCharge().identity()==previousCharge,"physical charge source differs from the previous accepted output owner");
    const auto& end=moved->endpoint();end.validateFor(physicalMesh,moved->sourceCharge(),step->transport(),t);
    require(end.geometry().identity()==pub.identity()&&sameClock(end.timeLayer(),t),"physical charge endpoint has a different geometry/full End clock");
    field->validateFor(physicalMesh,pub,moved->sourceCharge(),t);
    require(!field->isElectrostaticSnapshot()&&field->effectiveDt()==t.end-t.begin
        &&sameClock(field->geometryTime(),t)&&sameClock(field->candidateTime(),t)
        &&sameClock(field->chargeSourceTime(),s)&&field->sourceCharge().sameInstance(moved->sourceCharge())
        &&field->advection()&&field->advection()->sameInstance(end)
        &&field->implicitInputCharge()==end.integratedCharge(),"physical field lacks the actual final End advection/Ohmic owner");
    const bool changed=remap&&remap->interfaceTransfer;
    require(changed?(qremap&&field->chargeCandidate().sameInstance(qremap->sourceCharge())):
        field->chargeCandidate().sameInstance(charge),"physical End candidate is not the accepted/remap-source charge owner");
    const auto sourceId=source.phaseGeometry().vofIdentity().value();
    const auto physicalId=pub.identity().value(),acceptedId=accepted.phaseGeometry().vofIdentity().value();
    const auto& ar=end.report();
    require(field->report().usedAdvection,"physical Ohmic field is missing its typed advected input");
    const double sourceCharge=chargeTotal(field->sourceCharge().integratedCharge()),
        inputCharge=chargeTotal(field->implicitInputCharge()),
        physicalCharge=chargeTotal(field->chargeCandidate().integratedCharge()),current=boundaryCurrent(*field);
    auto row=stream();row<<t.step<<','<<t.time()<<','<<t.begin<<','<<t.end<<','<<static_cast<int>(t.position)<<','
        <<s.step<<','<<s.begin<<','<<s.end<<','<<static_cast<int>(s.position)<<','<<sourceId<<','<<physicalId<<','<<acceptedId<<','
        <<moved->sourceCharge().identity()<<','<<field->chargeCandidate().identity()<<','<<charge.identity()<<','<<changed<<','<<ar.boundaryEnabled<<','
        <<sourceCharge<<','<<inputCharge<<','<<physicalCharge<<',';
    // The changed target total is a separate remap observation;
    // source/input/physical-End totals above come from their actual arrays.
    row<<(changed?qremap->report().candidateTotalCharge:physicalCharge)<<','
        <<static_cast<long double>(physicalCharge)-sourceCharge<<','
        <<ar.boundaryIncomingCharge<<','<<ar.boundaryOutgoingCharge<<','<<ar.boundaryNetCharge<<','
        <<ar.chargeResidual<<','<<ar.relativeChargeResidual<<','<<ar.roundingResidual<<','
        <<ar.boundaryProductRoundingResidual<<','<<ar.boundaryModelChargeResidual<<','<<ar.relativeBoundaryModelChargeResidual<<','
        <<current<<','<<static_cast<long double>(field->effectiveDt())*current<<'\n';rows.charge=row.str();
    const auto& ledger=step->transport().report().boundaryFull;
    if(ledger.enabled){
        require(ledger.sourceIdentity==source.phaseGeometry().vofIdentity()&&sameClock(ledger.sourceTime,s)
            &&sameClock(ledger.targetTime,t),"physical phase boundary ledger has a different source/End clock");
        row=stream();row<<t.step<<','<<t.time()<<','<<t.begin<<','<<t.end<<','<<static_cast<int>(t.position)<<','
            <<sourceId<<','<<physicalId<<','<<acceptedId;
        for(int phase=0;phase<2;++phase)row<<','<<ledger.source[phase]<<','<<ledger.target[phase]<<','<<ledger.inflow[phase]<<','<<ledger.outflow[phase]
            <<','<<ledger.compression[phase]<<','<<ledger.roundingError[phase]<<','<<ledger.physicalResidual[phase]<<','<<ledger.algebraicResidual[phase];
        row<<','<<ledger.physicalRelativeConservation<<','<<ledger.algebraicRelativeConservation<<','<<ledger.massResidual<<','<<ledger.massRelativeConservation<<'\n';rows.phase=row.str();
    }
    require(!step->report().stages.empty(),"physical step has no final AC stage work report");
    const auto& work=step->report().stages.back();row=stream();
    row<<t.step<<','<<t.time()<<','<<t.begin<<','<<t.end<<','<<static_cast<int>(t.position)<<','<<sourceId<<','<<physicalId<<','<<acceptedId<<','
        <<bool(source.vofPrescribedFlow())<<','<<work.pressureOutletWork<<','<<work.prescribedReactionWork<<','
        <<work.boundaryKineticWork<<','<<work.boundaryPhysicalKineticWork<<','<<work.inletDissipationWork<<','
        <<work.boundaryMomentumRate[0]<<','<<work.boundaryMomentumRate[1]<<','<<work.kineticChange<<','
        <<work.energy.kineticBudgetResidual<<','<<work.energy.closureResidual<<','<<work.energy.closureRelativeResidual<<',';
    // This is the actual final-stage passive split, already inside BKE.
    // Legacy one-way/closed sources have no typed representation-work receipt.
    if(const auto* flow=source.vofPrescribedFlow();flow&&flow->pressureOutletBackflowModel)
        row<<work.boundaryFluxRoundingWork;
    else row<<"NA";
    row<<'\n';rows.work=row.str();
    if(changed&&remap->fluidTransfer->report().prescribedFlow){
        const auto& r=remap->fluidTransfer->report();row=stream();
        row<<t.step<<','<<t.time()<<','<<t.begin<<','<<t.end<<','<<static_cast<int>(t.position)<<','<<physicalId<<','<<acceptedId<<','
            <<r.projectionBoundaryWork<<','<<r.projectionRhsRoundingWork<<','<<r.projectionFreeResidualWork<<','<<r.projectionContinuityWork<<','
            <<r.projectionKineticChange<<','<<r.projectionCorrectionEnergy<<','<<r.projectionEnergyIdentityResidual<<','
            <<r.projectionEnergyIdentityRelativeResidual<<','<<r.kineticChange<<'\n';rows.projection=row.str();
    }
    return rows;
}
TecplotWriter::Zone electricFaces(const QuadtreeMesh& mesh,const E::ElectricVofField& field){
    TecplotWriter::Zone z;z.title="VOF-EHD native electric faces";z.element=TecplotWriter::Element::Line;
    z.coordinateNames=field.geometry().domain().axisymmetric()?std::array<std::string,2>{{"R","Z"}}:std::array<std::string,2>{{"X","Y"}};
    z.solutionTime=field.candidateTime().time();
    for(const char* name:{"native_face","negative_cell","positive_cell","physical_measure","H_integrated","D_integrated","J_integrated"})
        z.fields.push_back({name,TecplotWriter::Location::CellCentered,{}});
    for(const auto& f:field.faces()){
        const auto face=field.geometry().domain().face(mesh,f.nativeRow);
        const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
        const int first=static_cast<int>(z.points.size());z.points.push_back(mesh.verts(edge.edgeToVert[0]).coord());z.points.push_back(mesh.verts(edge.edgeToVert[1]).coord());
        z.connectivity.push_back({first,first+1});
        const std::array<double,7> values{{double(f.nativeRow),double(f.negativeCell),double(f.positiveCell),f.physicalMeasure,f.integratedField,f.displacement,f.current}};
        for(std::size_t k=0;k<values.size();++k)z.fields[k].values.push_back(values[k]);
    }
    z.auxiliary["CanonicalDirection"]="positive X/R or Y/Z; H is integrated negative potential gradient; D/J include physical measure once";
    return z;
}
void electricMetadata(TecplotWriter::Zone& z,const E::ElectricVofCharge& q,const E::ElectricVofField& field){
    const auto t=field.candidateTime(),s=field.chargeSourceTime();z.auxiliary["App"]=app;
    z.auxiliary["ChargeIdentity"]=text(q.identity());z.auxiliary["ElectricGeometryIdentity"]=text(field.geometry().identity().value());
    z.auxiliary["ElectricStep"]=text(t.step);z.auxiliary["ElectricBegin"]=text(t.begin);z.auxiliary["ElectricEnd"]=text(t.end);z.auxiliary["ElectricPosition"]=text(static_cast<int>(t.position));
    z.auxiliary["ElectricSourceStep"]=text(s.step);z.auxiliary["ElectricSourceBegin"]=text(s.begin);z.auxiliary["ElectricSourceEnd"]=text(s.end);z.auxiliary["ElectricSourcePosition"]=text(static_cast<int>(s.position));
    z.auxiliary["ElectricEffectiveDt"]=text(field.effectiveDt());z.auxiliary["ElectricFieldKind"]=field.isElectrostaticSnapshot()?"electrostaticSnapshot":"ohmicEnd";
    z.auxiliary["ElectricMeaning"]="accepted End/Initial diagnostic; never the applied Midpoint Maxwell field";
}
void sharpMetadata(TecplotWriter::Zone& z,const E::ElectricVofSharpField& field){
    const auto t=field.candidateTime(),s=field.sourceCharge().timeLayer();z.auxiliary["App"]=app;
    z.auxiliary["ChargeIdentity"]="NA";z.auxiliary["ElectricMethod"]="vofSharpSurface";
    z.auxiliary["ElectricGeometryIdentity"]=text(field.geometry().sourceIdentity().value());
    z.auxiliary["ElectricStep"]=text(t.step);z.auxiliary["ElectricBegin"]=text(t.begin);z.auxiliary["ElectricEnd"]=text(t.end);
    z.auxiliary["ElectricPosition"]=text(static_cast<int>(t.position));z.auxiliary["ElectricSourceStep"]=text(s.step);
    z.auxiliary["ElectricSourceBegin"]=text(s.begin);z.auxiliary["ElectricSourceEnd"]=text(s.end);
    z.auxiliary["ElectricSourcePosition"]=text(static_cast<int>(s.position));z.auxiliary["ElectricEffectiveDt"]=text(field.effectiveDt());
    z.auxiliary["ElectricMeaning"]="accepted Initial/End diagnostic; actual applied Midpoint Maxwell work is separate";
}
TecplotWriter::Zone sharpInterface(const E::ElectricVofSurfaceCharge& charge,const E::ElectricVofSharpField& field){
    TecplotWriter::Zone z;z.title="B1 actual smooth interface";z.element=TecplotWriter::Element::Line;
    z.coordinateNames={"X","Y"};z.solutionTime=field.candidateTime().time();
    for(const char* name:{"arc_id","native_cell","arc_parameter","surface_charge_density","En_liquid","En_gas","Et","arc_Q0","arc_Q1","arc_Q2"})
        z.fields.push_back({name,TecplotWriter::Location::CellCentered,{}});
    const auto& arcs=field.geometry().arcs();constexpr unsigned panels=8;
    for(std::size_t a=0;a<arcs.size();++a)for(unsigned k=0;k<panels;++k){
        const double t=(static_cast<double>(k)+.5)/panels;const auto s=field.interfaceField(a,t);
        const auto lo=arcs[a].sample(static_cast<long double>(k)/panels),hi=arcs[a].sample(static_cast<long double>(k+1)/panels);
        const int first=static_cast<int>(z.points.size());z.points.push_back({static_cast<double>(lo.x),static_cast<double>(lo.y)});
        z.points.push_back({static_cast<double>(hi.x),static_cast<double>(hi.y)});z.connectivity.push_back({first,first+1});
        const auto& q=charge.moments()[a];const std::array<double,10> values{{static_cast<double>(a),static_cast<double>(arcs[a].nativeCell),t,
            s.surfaceChargeDensity,s.normalElectricField[0],s.normalElectricField[1],s.tangentialElectricField,q[0],q[1],q[2]}};
        for(std::size_t j=0;j<values.size();++j)z.fields[j].values.push_back(values[j]);
    }
    z.auxiliary["ArcMomentMeaning"]="Q0/Q1/Q2 are repeated whole-arc annotations for plotting, not independent panel inventories; sum once per arc_id";
    z.auxiliary["Tessellation"]="8 display chords per actual quartic arc; field sampled on the actual curve";sharpMetadata(z,field);return z;
}
}
struct VofElectrohydrodynamicsOutput::Data{
    Data(const VofElectrohydrodynamicsInput& i,Context c):input(i),root(std::move(c.runDirectory)/"output"){}
    const VofElectrohydrodynamicsInput& input;fs::path root;bool begun=false,ioFailed=false,complete=false,finalSeen=false;
    std::int64_t step=-1;double time=0;std::uint64_t geometry=0,charge=0;std::string status="RUNNING",started,hashes;
    std::optional<EhdClock> acceptedClock,restoredClock;
    std::uint64_t frameOffset=0;
    Fluid::Totals totals;
    struct StepCompletion {
        EhdClock clock;std::uint64_t geometry=0;
        int actualAc=0,cells=0,velocityDofs=0,pressureDofs=0;
    };
    std::optional<StepCompletion> pendingCompletion;
    std::int64_t completedPerformanceStep=-1;
    void rememberCompletedOutput(const QuadtreeMesh& mesh,const H::HydroTwoPhaseState& fluid,const H::HydroVofStep* physicalStep){
        const auto t=fluid.timeLayer();
        if(t.position!=EhdClock::Position::End){pendingCompletion.reset();return;}
        const auto& layout=fluid.values().layout();
        pendingCompletion=StepCompletion{t,fluid.phaseGeometry().vofIdentity().value(),
            physicalStep?physicalStep->report().linearSolves:0,mesh.cellIndexNum(),layout.velocityCount(),layout.pressureCount()};
    }
    struct Frame{std::int64_t step;double time;std::uint64_t geometry,charge;std::string volume,surface,faces,kind;};
    std::vector<Frame> frames;
    std::map<std::string,std::unique_ptr<TecplotWriter::TableHistory>> histories;
    std::map<std::string,std::string> names;
    void table(const std::string& name,const std::string& row){
        if(row.empty())return;auto& p=histories[name];if(!p){const auto i=names.find(name);require(i!=names.end(),"series schema is missing");
            std::map<std::string,std::string> metadata{{"App",app},{"Interface","LocalPlicVOF"},{"Series",name}};
            if(name=="open_convection_stages"){
                metadata["StatisticsScope"]="two actual different-U evaluations for every AC stage; material advector and required solved-midpoint diagnostics";
                metadata["TimingScope"]="nested in existing convection and post diagnostics; children disjoint within actual apply, not additional step time";
                metadata["ParallelScope"]="actual consumed bounded chunk contributions and original-native-order gather; global ledger order retained; worker maxima are not CPU busy percentage";
                metadata["ScratchScope"]="conservative optional plan/chunk/bigint/worker estimate capped at64MiB; excludes original mandatory result/rates storage";
            }
            if(name=="reconstruction_templates"){
                metadata["StatisticsScope"]="actual owned regular-template bank snapshot at Initial/accepted End";
                metadata["AttemptCounters"]="row_queries,row_hits,row_builds,port_queries,port_hits,port_builds,actual_svd_calls are cumulative bank attempts, including failed preparations; not accepted-only";
                metadata["StorageFields"]="entries and retained_bytes are current owned bank storage; bytes are the accounted retained estimate";
                metadata["Units"]="counts and bytes; not elapsed time and not added to performance series";
            }
            if(name=="stress_geometry_reuse"){
                metadata["ReportSource"]="actual final-AC frozenMidpoint.forces.maxwellLoad.geometry.geometryReuseReport; no row without that immutable owner";
                metadata["SnapshotScope"]="snapshot of the geometry construction actually used; retained geometry reuse may repeat this snapshot at later accepted Ends";
                metadata["NotStepWork"]="not newly executed work in this step; do not sum these snapshots across steps or add counts to elapsed time";
                metadata["RowProvenance"]="step/time/full clock and accepted_geometry_identity describe current accepted End only; no geometry-construction identity or clock is inferred";
                metadata["Availability"]="Initial and accepted states without an actual Maxwell owner have no invented row; requested/enabled are actual report flags, not overall physical acceptance";
                metadata["CsrRebindScope"]="original direct/literal counters cover actual scalar source copies in shared face, cell integral and final Ehat rows only; direct requires ascending unique rebound columns, otherwise literal; counts are construction snapshots";
                metadata["RawVolumeCsrRebindScope"]="two appended counters are actual direct copied raw/volume scalar rows, including empty rows, after a successful fresh product assembly; nonascending/unqualified/resource paths remain literal; construction snapshots may repeat and are never summed across accepted steps";
                metadata["CsrScratchScope"]="accounted extra routing and temporary CSR construction peak <=64MiB, not RSS or cumulative; zero is an actual zero work-space report";
            }
            if(name=="local_plic_execution"){
                metadata["ReportSource"]="only actual retained final-AC step.transport.report.localFaceExecution; no Initial row, never an all-AC/retry ledger";
                metadata["ClockScope"]="actual source geometry/full previous-accepted clock, Midpoint and Endpoint identities; physical End and accepted-after-AMR identities are separate output provenance";
                metadata["CountScope"]="native_faces/incidences are one invocation route size; directions/evaluations/gathers aggregate successful half/full directional execution in that invocation, not unique geometry stock or newly allocated work at every output";
                metadata["ScratchScope"]="scratch_bytes is per-invocation extra optional execution peak <=64MiB, not RSS or retained stock; repeated reports are never summed";
                metadata["Availability"]="route/count/scratch fields require eligible; maximum_workers requires a successful direction; absent is NA, never fabricated zero";
                metadata["Flags"]="requested/eligible/face-used/gather-used/fallback describe actual execution; fallback may follow earlier successful directions, never a physical acceptance rule";
            }
            if(name=="performance_step_complete"){
                metadata["WallScope"]="controller start before stability/adaptive/retries through accepted commit, time_step, payload/series/index/manifest, releaseStep and recovery bookkeeping";
                metadata["Excluded"]="this timing row's own append/flush and final whole-App complete/close";
                metadata["AcceptanceScope"]="one matching actual accepted End, once only; alternative AC stages are not separate steps";
            }
            if(name=="performance_amr"){
                metadata["ReceiptScope"]="actual successful accepted same-End remap receipts; full-step time includes retries separately";
                metadata["NestedTimes"]="tree indicator/planning/adaptation belong to tree; fluid raw_flux/projection belong to fluid; target geometry/assemble/factor/solve belong to snapshot";
                metadata["Unattributed"]="regrid_prepare minus top-level tree, inventory, fluid, charge, snapshot and available target-stress/final-joint-validation intervals; unclamped diagnostic, not a conservation gate";
                metadata["AmrChildTimes"]="target_stress_geometry covers the actual requested Full/RT0-Green remapped or mode-switch cold target factory and is available only when its target owner is retained; final_joint_validation covers only the original final accepted target validate, excluding construction; neither is nested in fluid/snapshot reports";
                metadata["FluidPreparationChildren"]="eight disjoint fluid preparation intervals plus existing raw_flux/projection are nested in fluid_remap; do not add children again to parent/regrid/complete time";
                metadata["FluidTargetState"]="target_state is the complete actual createVof*Remapped call, including reconstruction and its full validation; not exclusive Regular time";
                metadata["FluidPreparationResidual"]="fluid_unattributed is the unclamped original fluid total minus ten children, an elapsed-time diagnostic rather than a physical residual; original total excludes final transfer Data packing/return/destruction";
                metadata["FluidPreparationAvailability"]="all nine new fluid columns require changed actual fluid receipt and preparation.available; otherwise NA; measured zero is valid";
                metadata["TreeCloneStorage"]="12 current-decision integer stocks: four logical capacities, source constructed slots before clone, and target constructed slots immediately after clone before adaptation; order vertices/horizontal-edges/vertical-edges/cells; NA if no actual clone; not RSS, cumulative work or final target storage";
                metadata["Availability"]="no-op retains tree/total; unexecuted operations, absent target stress owner or unsupported passive timing are NA; disabled/Initial has no invented row";
            }
            if(name=="solver_stages"){
                metadata["SaddleScatterPlan"]="actual Saddle report; requested/eligible/hit/system-used/quotient-used are independent; pure system/quotient index routing with fresh matrix values; cold symbolic build may be used in this same call with hit=false; unsupported/resource fallback keeps original packing";
                metadata["MatrixExecution"]="appended execution fields are actual current factory work; symbolic_built is not a cache hit; symbolic_workers observes successful integer preparation; value workers require successful parallel use; fallback is original arithmetic replay, never a physical-method fallback";
                metadata["IdentityQuotient"]="quotient_aliases_system means the actual no-gauge identity quotient factor borrows this same immutable system; quotient_scatter_used remains false; full augmented residual and physical gates still use original system";
                metadata["SaddleScatterAvailability"]="system/quotient nonzeros require requested; bytes require a participating retained plan and are stock, never summed; no numeric factor or prior State reuse";
                metadata["SparseValidationPlan"]="actual velocity-block report; requested/eligible/hit/exact symmetry/skew are separate flags; pure transpose-index plan, never cached matrix values";
                metadata["SparsePlanAvailability"]="nonzeros require requested; retained bytes require an actual participating plan and are stock, never summed; hit alone does not imply exact validation used; fallback retains original gates";
            }
            if(name=="frozen_matrix_preparation"){
                metadata["ReportSource"]="only actual retained final-AC frozenMidpoint.matrixAssemblyPlanReport; not an all-AC/retry ledger";
                metadata["ClockScope"]="physical-End and accepted-after-AMR provenance are separate from actual material Midpoint source geometry and full clock";
                metadata["PlanScope"]="physical K+diagonal then physical+freshH union; every physical coefficient remains fresh; flags distinguish pattern hit from numerical use";
                metadata["StorageScope"]="retained_bytes is current single pure-index slot stock bounded by 128MiB, never cumulative; absence is NA";
                metadata["ColdExecution"]="physical/union cold_compiled describe actual new pure-index routes; cold_used requires this call actually evaluated values with that new route, independently of prior-plan hit";
                metadata["ValueExecution"]="parallel_used and fallback may both occur in separate phases; maximum_workers requires a successful current numerical team; values_scratch_bytes is peak added optional worker storage, not matrix stock, RSS, or accumulated work; fields observe final AC only";
            }
            if(name=="regular_preparation_reuse"){
                metadata["ReportSource"]="actual one-shot authenticated target Regular construction; emitted only at changed remap; no row without requested receipt";
                metadata["ClockScope"]="physical-source and accepted VOF identities are remap output provenance, not invented Regular construction identities";
                metadata["StorageScope"]="source_metadata_bytes is the source dependency-certificate stock observed at construction, not current process/bank storage or cumulative memory";
            }
            if(name=="charge_transport_stages"){
                metadata["ReportSource"]="only actual retained final-AC Midpoint advection and physical-End pair Midpoint/Endpoint advances; not all AC/retries or a complete provider ledger";
                metadata["ClockScope"]="actual source Q identity/fullclock and actual transport source geometry plus stage geometry/fullclock; no candidate Q identity exists on this vector receipt";
                metadata["ParallelAvailability"]="workers/scratch only when actual donor or target multiworker segment used; scratch is the per-call optional budget observation, never cumulative process memory";
            }
            if(name=="electric_preparation"){
                metadata["ReportSource"]="actual final-AC A Maxwell stress.field report and retained physical End field/charge transport; not an all-AC ledger";
                metadata["LinearBackend"]="used_gmres is the actual explicit solver backend; iterations/restarts sum the initial solve and original at most two same-operator defect solves; maximum_workers is an observed maximum, not configured threads; counts are NA without GMRES";
                metadata["ClockScope"]="actual field candidate clock/geometry/source/candidate charge identities; step/time and accepted geometry describe output provenance only";
                metadata["FieldTimes"]="geometry/assemble/factor/solve are nested diagnostics in field total; total begins inside prepareImpl and excludes the outer typed advection-entry checks";
                metadata["ChargeTimes"]="advection_stage is only the actual retained Midpoint advance; physical End pair columns measure its separate midpoint/endpoint advance calls, not the full advectCharge provider";
                metadata["Unmeasured"]="outer transport/boundary validation, boundary vector preparation and stress preparation are not inferred from these reports; no sum is called full Maxwell provider wall";
                metadata["Availability"]="Initial or absent actual owners have no invented row; charge times without the corresponding real receipt are NA";
                metadata["PostprocessTimes"]="actual face/cell/finalize substages and enclosing postprocess time; finalize includes original global/receipt gates and requested cell gradients; setup and required serial replay remain in the parent, not omitted work; nested in Field total";
                metadata["ExpandedResidual"]="actual Solver complete-matrix three-term exact residual receipt, forwarded by this Field; flags/workers/rows/scratch describe completed or attempted optional work, not configured threads or a replacement linear solve";
                metadata["ExpandedResidualStorage"]="scratch bytes is the accounted route-preparation/final-residual-helper optional peak, not process RSS, retained history or measurement of prior defect-helper peaks; route time includes the necessary matrix/RHS/diagonal scan; residual time covers exact contraction/join or literal replay before the unchanged norm/relative gate; unavailable producer/timing fields are NA";
            }
            if(name=="csg_boundary_preparation"){
                metadata["ReportSource"]="actual App csgBoundaryIds direct-call branches only; absent report or zero requests emits no row";
                metadata["Scope"]="numeric0 Initial/restore;1 physical candidate all AC plus End;2 AMR target boundary consumers; per-output successful attempt scopes, not a cumulative global selector ledger";
                metadata["Excluded"]="separate AMR targetDomain provider and Core-internal graph work are excluded from selector_calls/complete_graph_builds";
                metadata["Provenance"]="step/time and physical_geometry_identity/accepted_geometry_identity describe current output owners; they are not cached geometry construction clocks or IDs";
                metadata["Counters"]="requests/hits/rebuilds/retained_builds/selector_calls/complete_graph_builds count actual helper branches; counters alone do not certify numerical equivalence";
                metadata["Timing"]="total_seconds encloses direct helper requests and rebuild_seconds is its nested original classification branch; do not sum these times or add them to parent performance buckets; unavailable times are NA";
                metadata["Storage"]="peak_retained_bytes is the observed single-slot retained payload peak under32MiB, not summed stock, process RSS or total preparation memory";
            }
            if(name=="viscosity_product_stages"){
                metadata["ReportSource"]="actual HydroVofStep AC-stage viscosity Energy product report; no cached physical State";
                metadata["ClockScope"]="physical End/accepted target provenance is separate from actual weight-source geometry and full Midpoint clock saved from the built PhaseViscosity; every executed AC stage, not accumulated work";
                metadata["MemoScope"]="same actual immutable raw E owner; fresh complete weights and weighted E; dirty actual E-column support rows recomputed in full original k order";
                metadata["FPReceipt"]="clean rows merge actual pruned K with complete symbolic zero slots; all-clean calls share the immutable matrix and replay genuine row plus canonical prune flags; source/current controls qualified; unsupported executes complete original assembly";
                metadata["FallbackReason"]="0 qualified row work; 1 first Eigen/no row receipt; 2 unsupported FP; 3 optional memory scope; 4 concurrent slot busy; 5 optional failure";
                metadata["Availability"]="Initial has no invented row; qualified-only counts are NA on fallback; retained bytes include actual cached K allocated values/indices/outer and small metadata; peak reserves max(old,new) one optional K plus simultaneous metadata/plan/workers, separately from mandatory new output; bounds include allocator headroom and are never summed across AC or RSS";
                metadata["NoInferredFraction"]="changed_weight_rows and recomputed/reused_rows are actual bit comparison and actual support closure, not unchanged-cell or AMR fractions";
                metadata["NoSlipReceipt"]="current Owned raw preparation only; public Energy build and previous/canonical raw owners do not repeat historical preparation";
                metadata["NoSlipAvailability"]="no_slip_preparation_available is always present;0 makes all13 remaining no-slip values unavailable;1 permits real zero work/counts/stock";
                metadata["NoSlipTiming"]="cold_seconds is complete actual cold raw construction nested in viscosity parent; missing timing is NA; never add the child to its parent";
                metadata["NoSlipWorkers"]="fit_workers is the observed peak attempted fit team; parallel_used and parallel_fallback are independent; raw-cache hit has no fit work/workers; requested threads alone do not prove parallel execution";
                metadata["NoSlipCounts"]="fit_requests/hits count actual complete coefficient-fit helper entries/hits; QR/SVD count actual Eigen factorization constructors including discarded attempts and original-order replay; do not infer factorizations from misses";
                metadata["NoSlipStorage"]="bank_peak_bytes and bank_limit_bytes are bounded logical optional mathematical bank stock, not RSS; zero with Available1 is actual absence of a current bank";
            }
            if(name=="maxwell_stages"){
                metadata["ReportSource"]="actual HydroVofStep.report.stages[ac_stage-1].maxwellLoadReport; no row without actual Maxwell owner";
                metadata["ClockScope"]="physical End interval and accepted target identities only; no force Midpoint identity is inferred";
                metadata["PlanScope"]="compiled geometry/index/exact-coefficient payload; each actual stress and tensor remains fresh";
                metadata["TimingScope"]="nested Maxwell owner intervals; sampling/contraction share the original worker team interval; no added barrier";
                metadata["Availability"]="term/byte counts require compiled_contraction; group times require compiled_timing_available, otherwise NA; flags do not assert physical accuracy";
            }
            if(name=="solver_stages"){
                metadata["ReportSource"]="actual HydroVofStep.report.stages[ac_stage-1].saddle; every executed AC of the accepted physical step";
                metadata["VelocityAppendScope"]="velocity_append_plan_used is fresh rounded S+=W plus original transpose-defect/norm using authenticated existing routing; scratch is a per-call scalar-buffer budget <=64MiB and is NA unless used";
                metadata["DirectInputStorage"]="current-factor stock: deferred means actual full-input certification omitted eager DirectSparse materialization; materialized means actual DirectSparse exists, including eager unsupported inputs or earlier fallback; flags are not per-solve Direct counts; lazy conversion belongs to direct_factor_seconds on actual fallback";
                metadata["OrderedSpmvScope"]="current-factor owned compact row view of fresh immutable quotient values; stable original CSC encounter order per row; only original GMRES main Q*z calls; no cross-factor numeric/cache reuse";
                metadata["OrderedSpmvCalls"]="calls/parallel_calls count actual main products within this public solve including outer linear refinements; used means successful parallel_calls>0; fallback means a main product used original Eigen and is not a direct or physical fallback";
                metadata["OrderedSpmvStorage"]="view_bytes is current-factor actual owned view stock, available only when eligible; total optional peak guard <=64MiB includes construction cursor and worker scratch; maximum_workers is maximum successful actual team, NA unless used; stocks and per-solve counts are not cumulative across rows";
                metadata["ParallelPressure"]="actual fixed block-Jacobi pressure preconditioner; block count and observed maximum_workers are NA without prepared blocks; these are current solve observations, not configured thread counts";
                metadata["PressureCoarse"]="current-factor same-S core256 Galerkin correction; prepared is admission, applications is actual solve usage including refinements; not prepared fixes original RAS before GMRES; bounds require prepared; no physical acceptance or speed inferred";
                metadata["PressureCoarseReason"]="0 not prepared,1 too few cores,2 environment,3 matrix,4 storage,5 rank,6 inverse certificate,7 condition or formation,8 workspace,9 prepared; all fields NA without actual pressure blocks";
                metadata["GlobalPressure"]="explicit GMRES-only same-S pressure LDLT, never full-K direct fallback; factor_nonzeros is stored strict L with implicit unit diagonal and separate D";
                metadata["GlobalPressureAvailability"]="requested is configuration; available requires an actual prepared factor; indexed storage/admission are payload models, not RSS; real symbolic/numeric preparation and solve_calls/time including refinements; absent data/timing and unused RAS/coarse are NA";
                metadata["ClockScope"]="row clock and identities are physical End interval/accepted target; no stage Midpoint geometry identity is inferred";
                metadata["Times"]="original report diagnostics; preconditioner/factor/direct-factor scopes may nest; never add these columns to complete-step wall";
                metadata["FallbackReason"]="0 none; 1 preconditioner qualification; 2 inner Krylov; 3 original outer residual gates";
                metadata["VelocitySpdTier"]="actual successful validation: 0 none, 1 positive diagonal, 2 strict diagonal dominance, 3 certified energy Gram, 4 original LLT";
                metadata["VelocitySpdBounds"]="upper residual/error row bound and lower scaled mass are available only when gram_bounds_available=1; unavailable/unsupported bounds are NA; computed nonqualifying bounds remain available";
                metadata["VelocitySpdTimes"]="certificate includes diagonal dominance and any Gram attempt; LLT covers only actual check.compute; both are nested in velocity_block_seconds";
                metadata["Acceptance"]="backend flags describe actual work, not success by themselves; original outer residual gates remain authoritative";
                metadata["FrozenSegments"]="actual input-binding, force validation, matrix/actions, boundary/RHS, whole velocity preparation; unknown groups are NA";
                metadata["SaddleSegments"]="actual partition, augmented system assembly, quotient assembly; RHS elimination is nested in solve_seconds";
                metadata["NestedSegments"]="Frozen velocity preparation includes SPD/certificate/LLT; component intervals must not be added to complete-step wall";
            }
            p=std::make_unique<TecplotWriter::TableHistory>(root/"series"/(name+".dat"),name,schema(i->second),std::move(metadata));}
        p->append(row);
    }
    void flush(){for(auto& h:histories)if(h.second)h.second->flush();}
    std::string index(const std::vector<Frame>& f)const{
        const bool sharp=input.method==VofElectrohydrodynamicsInput::Method::VofSharpSurface;
        auto s=stream();s<<"frame,step,time,geometry_identity,charge_identity,volume,interface,"<<(sharp?"electric_interface":"electric_faces")<<",field_kind\n";
        for(std::size_t k=0;k<f.size();++k){const auto& x=f[k];s<<frameOffset+k<<','<<x.step<<','<<x.time<<','<<x.geometry<<',';if(sharp)s<<"NA";else s<<x.charge;s<<",fields/"<<x.volume<<',';
            if(!x.surface.empty())s<<"fields/"<<x.surface;s<<",fields/"<<x.faces<<','<<x.kind<<'\n';}return s.str();
    }
    void manifest(){
        const bool sharp=input.method==VofElectrohydrodynamicsInput::Method::VofSharpSurface;
        auto s=stream();s<<"{\"app\":"<<quoted(app)<<",\"status\":"<<quoted(status)<<",\"started_at\":"<<quoted(started)
            <<",\"accepted_step\":"<<std::max<std::int64_t>(0,step)<<",\"accepted_time\":"<<time<<",\"geometry_identity\":"<<geometry<<",\"charge_identity\":";
        if(sharp)s<<"null";else s<<charge;
        s
            <<",\"geometry_identity_available\":"<<(geometry?"true":"false")<<",\"charge_identity_available\":"<<(charge?"true":"false")
            <<",\"complete_frames\":"<<frames.size()<<",\"frame_offset\":"<<frameOffset
            <<",\"cumulative_complete_frames\":"<<frameOffset+frames.size()
            <<",\"restart_checkpoint_directory\":"<<quoted(input.fluid.run.restartDirectory)
            <<",\"restored_checkpoint_step\":"<<(restoredClock?std::to_string(restoredClock->step):"null")
            <<",\"volume_format\":\"tecplotDat\",\"series_format\":\"tecplotDat\",\"frame_index\":\"meta/frames.csv\""
            <<",\"series_profile\":"<<quoted(input.fluid.run.diagnosticSeries?"diagnostic":"physical")
            <<(sharp?",\"electric_interface_format\":\"tecplotDat\"":",\"electric_faces_format\":\"tecplotDat\"")
            <<",\"electric_field\":\"accepted diagnostic; actual Midpoint Maxwell work is separate\""
            <<",\"electric_method\":"<<quoted(sharp?"vofSharpSurface":"vofFiniteVolume")
            <<(sharp?std::string{}:std::string(",\"charge_advection\":\"")
                +(input.controls.chargeAdvection==E::ElectricVofChargeAdvectionScheme::LimitedLinearFct?"limitedLinearFct":"donorCell")+"\"")
            <<",\"charge_representation\":"<<quoted(sharp?"independent signed surface moments; neutral phase interiors":"integrated signed Eulerian cell charge")
            <<",\"input_hashes\":"<<hashes<<"}\n";replace(root/"meta/manifest.json",s.str());
    }
};
VofElectrohydrodynamicsOutput::VofElectrohydrodynamicsOutput(const VofElectrohydrodynamicsInput& i,Context c):data_(std::make_unique<Data>(i,std::move(c))){}
VofElectrohydrodynamicsOutput::~VofElectrohydrodynamicsOutput()=default;
std::string VofElectrohydrodynamicsOutput::checkpoint() const{
    namespace Archive=LAMP::IO::Checkpoint;const auto& d=*data_;
    Archive::require(d.begun&&!d.complete&&!d.ioFailed&&d.status=="RUNNING"&&d.acceptedClock
        &&d.input.method==VofElectrohydrodynamicsInput::Method::VofFiniteVolume,
        "A output checkpoint requires a successful actual accepted output; B1 is unsupported");
    const auto t=*d.acceptedClock;(void)t.time();
    Archive::require(t.step==d.step&&t.time()==d.time&&t.position!=EhdClock::Position::Midpoint,
        "A output checkpoint has a different accepted full clock");
    Archive::require(std::isfinite(d.totals.initialVolume)&&d.totals.initialVolume>=0
        &&std::isfinite(d.totals.bodyWork)&&std::isfinite(d.totals.viscousLoss),"nonfinite A output reference totals");
    Archive::require(d.frames.size()<=std::numeric_limits<std::uint64_t>::max()-d.frameOffset,
        "A output cumulative frame counter overflow");
    Archive::Writer out("VofElectrohydrodynamicsOutput.A.v1");out.scalar(t.step);out.scalar(t.begin);out.scalar(t.end);out.scalar(t.position);
    out.scalar(d.totals.initialVolume);out.scalar(d.totals.bodyWork);out.scalar(d.totals.viscousLoss);
    out.scalar<std::uint64_t>(d.frameOffset+d.frames.size());return out.finish();
}
void VofElectrohydrodynamicsOutput::restoreCheckpoint(const QuadtreeMesh& mesh,const std::string& bytes,
    const H::HydroTwoPhaseState& restored,const E::ElectricVofCharge& charge){
    namespace Archive=LAMP::IO::Checkpoint;auto& d=*data_;
    Archive::require(d.begun&&!d.complete&&!d.ioFailed&&d.status=="RUNNING"&&d.step==-1&&!d.acceptedClock
        &&d.frames.empty()&&d.frameOffset==0&&d.input.method==VofElectrohydrodynamicsInput::Method::VofFiniteVolume,
        "A output restore requires a fresh begun A writer; B1 is unsupported");
    Archive::Reader in(bytes,"VofElectrohydrodynamicsOutput.A.v1");EhdClock t;
    t.step=in.scalar<std::int64_t>();t.begin=in.scalar<double>();t.end=in.scalar<double>();t.position=in.scalar<EhdClock::Position>();
    Fluid::Totals totals;totals.initialVolume=in.scalar<long double>();totals.bodyWork=in.scalar<long double>();totals.viscousLoss=in.scalar<long double>();
    const auto frames=in.scalar<std::uint64_t>();in.finish();(void)t.time();
    Archive::require(t.position!=EhdClock::Position::Midpoint&&frames<=static_cast<std::uint64_t>(t.step)+1,
        "A output checkpoint has an invalid accepted clock or frame count");
    Archive::require(std::isfinite(totals.initialVolume)&&totals.initialVolume>=0
        &&std::isfinite(totals.bodyWork)&&std::isfinite(totals.viscousLoss),"invalid A output reference totals");
    Archive::require(restored.model()==H::HydroPhaseGeometry::Model::GeometricVof&&sameClock(restored.timeLayer(),t),
        "A output checkpoint differs from the genuine restored full clock");
    const auto& publication=restored.phaseGeometry().vofPublication();
    restored.validateVofFor(mesh,publication.identity(),t);charge.validateFor(mesh,publication,t);
    d.step=t.step;d.time=t.time();d.geometry=publication.identity().value();d.charge=charge.identity();
    d.acceptedClock=t;d.restoredClock=t;d.totals=totals;d.frameOffset=frames;
    d.pendingCompletion.reset();d.completedPerformanceStep=t.step;d.finalSeen=false;d.manifest();
}
void VofElectrohydrodynamicsOutput::begin(){
    auto& d=*data_;require(!d.begun&&!fs::exists(d.root),"output requires a fresh run directory");
    require(d.input.fluid.interfaceMethod==HydrodynamicsInput::InterfaceMethod::Vof,"input is not VOF");
    const bool sharp=d.input.method==VofElectrohydrodynamicsInput::Method::VofSharpSurface;
    fs::create_directories(d.root);d.begun=true;for(const char* name:{"fields","series","meta/effective"})fs::create_directories(d.root/name);
    d.started=LAMP::IO::utc();auto hash=stream(),effective=stream();hash<<'{';unsigned count=0;
    for(const auto& doc:d.input.fluid.documents){if(count++)hash<<',';hash<<quoted(doc.file)<<':'<<quoted("sha256:"+LAMP::IO::sha256(doc.sourceText));
        write(d.root/"meta/effective"/fs::path(doc.file).filename(),doc.effectiveText);effective<<"// "<<fs::path(doc.file).filename().generic_string()<<'\n'<<doc.effectiveText<<'\n';}
    hash<<'}';d.hashes=hash.str();write(d.root/"meta/effective_config.dict",effective.str());
    auto trace=stream();trace<<"field,file,line,raw_value,effective_value,units,consumer,defaulted,used,unused_reason\n";
    for(const auto& e:d.input.fluid.trace.entries)trace<<LAMP::IO::csv(e.field)<<','<<LAMP::IO::csv(e.file)<<','<<e.line<<','<<LAMP::IO::csv(e.rawValue)<<','
        <<LAMP::IO::csv(LAMP::IO::value(e.value))<<','<<LAMP::IO::csv(e.units)<<','<<LAMP::IO::csv(e.consumer)<<','<<e.defaulted<<','<<e.used<<','<<LAMP::IO::csv(e.unusedReason)<<'\n';
    write(d.root/"meta/configuration_trace.csv",trace.str());
    write(d.root/"meta/build_info.txt",std::string("app=")+app+"\ncompiler="+LAMP::IO::compiler()+"\nbuild_type="+LAMP::IO::buildType()+
        "\nunits=code\n"+(sharp?std::string("primary=independent VOF stocks and signed surface moments\nsurface_charge_density=charge per actual interface measure; neutral bulk\ncharge_identity=NA\nelectric_interface=actual smooth curve, shared Et and equilibrated two-sided En\n"):
        std::string("primary=independent VOF stocks and integrated signed Eulerian charge\ncharge_density=integrated_charge/actual native XY or RZ capacity\nelectric_faces=actual native canonical integrated H/D/J; no cell interpolation\n"))+
        "field=accepted Initial/End diagnostic, snapshot or true End Ohmic\n"
        "maxwell_work=actual final AC stage only, original Midpoint owner; never summed over alternative stages\nremap=target snapshot, original physical work/charge remap separate\n"
        "physical_boundary_ledgers=actual final transport and old-tree End Ohmic receipts; absent receipts are unavailable\n"
        "dt_boundary_current_from_report_total=passive dt times rounded actual-face current at output; never exact Ohmic transfer\n"
        "ohmic_diagnostics=raw/finite/sigma/quantization ledger table and acceptance self-audits retired; actual transfer rule remains minimum-binary64-quantum RNE\n"
        "phase_boundary=physical Delta+out-in and WY/RN algebraic ledger are separate; relative_volume_drift is not an open balance\n"
        "vof_boundary_work=actual final AC only; regrid_projection_work=artificial same-End WB/WRHS, never physical timestep work\n"
        "surface_area=NA\nsurface_energy=NA\ntotal_energy=NA\nTecplot2017=readability and microphase precision are separate qualifications\n");
    d.names={{"conservation",conservationHeader},{"solver",solverHeader},{"energy",energyHeader},{"vof_csf",capillaryHeader},{"events",eventsHeader},{"time_step",timeStepHeader}};
    d.names["vof_native_velocity"]=Fluid::vofNativeVelocityCsvHeader();
    d.names["reconstruction_templates"]=reconstructionTemplatesHeader;
    d.names["performance_step_complete"]=stepCompleteHeader;
    d.names["performance_amr"]=amrPerformanceHeader;
    d.names["solver_stages"]=solverStagesHeader;
    d.names["maxwell_stages"]=maxwellStagesHeader;
    d.names["viscosity_product_stages"]=viscosityProductStagesHeader;
    d.names["csg_boundary_preparation"]=csgBoundaryPreparationHeader;
    d.names["electric_preparation"]=electricPreparationHeader;
    d.names["regular_preparation_reuse"]=regularPreparationReuseHeader;
    d.names["charge_transport_stages"]=chargeTransportStagesHeader;
    d.names["open_convection_stages"]=openConvectionStagesHeader;
    d.names["frozen_matrix_preparation"]=frozenMatrixPreparationHeader;
    d.names["stress_geometry_reuse"]=stressGeometryReuseHeader;
    d.names["local_plic_execution"]=localPlicExecutionHeader;
    d.names["electric"]="step,time,geometry_identity,charge_identity,field_kind,effective_dt,charge_source_step,charge_source_begin,charge_source_end,charge_source_position,total_charge,relative_linear_residual,relative_gauss_residual,instantaneous_boundary_current";
    if(sharp){d.names["electric"]="step,time,geometry_identity,charge_identity,field_kind,effective_dt,charge_source_step,charge_source_begin,charge_source_end,charge_source_position,total_charge,relative_linear_residual,relative_gauss_residual,instantaneous_boundary_current,global_gauss_relative,global_charge_update_relative,skeleton_dofs,nonzeros,assemble_seconds,factor_seconds,solve_seconds,recovery_seconds,mid_q1_prepare_seconds,mid_q1_integration_seconds,mid_q1_project_seconds,mid_q1_native_arc_seconds,mid_q1_core_arc_seconds,mid_q1_identity_arc_seconds,mid_q1_material_samples,mid_q1_endpoint_requests,mid_q1_endpoint_memo_hits,mid_q1_native_arcs,mid_q1_exact_arcs,mid_q1_core_roots,end_q1_prepare_seconds,end_q1_integration_seconds,end_q1_project_seconds,end_q1_native_arc_seconds,end_q1_core_arc_seconds,end_q1_identity_arc_seconds,end_q1_material_samples,end_q1_endpoint_requests,end_q1_endpoint_memo_hits,end_q1_native_arcs,end_q1_exact_arcs,end_q1_core_roots,mid_electric_assemble_seconds,mid_electric_factor_seconds,mid_electric_solve_seconds,mid_electric_recovery_seconds";
        d.names["electric_correction"]="step,time,checked,stages,force_change_relative,work_change_relative";
        d.names["surface_transport"]="step,time,source_charge,target_charge,relative_charge_residual,source_moment_quadrature_residual,maximum_displacement,maximum_projection_distance,material_samples,parameter_splits,identity_transport,integral_diagnostics_available,identity_arc_transfers,source_density_reconstructions";
        d.names["electric"]+=",mid_q1_core_hessian,mid_q1_core_endpoint_sign,mid_q1_core_nonfinite,mid_q1_core_root_width,mid_q1_core_forced_distance,mid_q1_core_endpoints,mid_q1_core_interiors,mid_q1_geometry_bases,end_q1_core_hessian,end_q1_core_endpoint_sign,end_q1_core_nonfinite,end_q1_core_root_width,end_q1_core_forced_distance,end_q1_core_endpoints,end_q1_core_interiors,end_q1_geometry_bases";}
    d.names["maxwell_work"]="step,accepted_time,force_time,force_geometry_identity,maxwell_work,physical_end_geometry_identity,accepted_geometry_identity,remapped";
    d.names["regrid"]="step,time,source_geometry_identity,target_geometry_identity,source_charge_identity,target_charge_identity,source_leaves,target_leaves,refined_donors,coarsened_targets,source_liquid,target_liquid,source_gas,target_gas,source_charge,target_charge,charge_residual,relative_charge_residual,delta_native_kinetic,history_reset,charge_linear_remap,charge_high_order_donors,charge_limited_donors,charge_finite_p0_donors,charge_max_linear_quantization,charge_linear_quantization_l1,charge_remap_seconds";
    d.names["regrid_decision"]="step,time,source_geometry_identity,target_geometry_identity,source_leaves,target_leaves,requested,changed,clones,sweeps,decision_seconds,inventory_transfer_seconds,regrid_seconds,charge_protection_enabled,additional_refine_cells,additional_hold_cells";
    d.names["charge_balance"]=chargeBalanceHeader;d.names["phase_boundary"]=phaseBoundaryHeader;
    d.names["vof_boundary_work"]=boundaryWorkHeader;d.names["regrid_projection_work"]=projectionWorkHeader;
    replace(d.root/"meta/frames.csv",d.index(d.frames));d.manifest();
}

void VofElectrohydrodynamicsOutput::accepted(const QuadtreeMesh& mesh,const H::HydroTwoPhaseState& fluid,
    const E::ElectricVofCharge& charge,const E::ElectricVofField& field,const H::HydroVofStep* step,
    Clock::time_point start,double compute,bool final,const HydrodynamicsVofRemapReceipt* receipt,const E::ElectricVofChargeMeshTransfer* qremap,
    const E::ElectricVofChargeTransport* physicalChargeTransport,const E::ElectricVofField* physicalEndField,
    const VofEhdCsgBoundaryReport* csgBoundaryPreparation){
    auto& d=*data_;require(d.begun&&!d.complete&&!d.ioFailed&&d.input.method==VofElectrohydrodynamicsInput::Method::VofFiniteVolume,"volume-charge writer is not active");
    const auto* transfer=validate(mesh,fluid,charge,field,step,receipt,qremap);const auto t=fluid.timeLayer();
    require(step?(t.step==d.step+1&&t.begin==d.time&&t.time()>d.time):(d.step==-1&&t.step==0&&t.time()==0),
        "output requires sequential actual accepted states");
    require(!final||t.time()==d.input.fluid.run.endTime,"final flag differs from actual accepted endTime");
    const auto rows=Fluid::vofRows(mesh,fluid,step,d.input.fluid,d.totals,compute,std::cout.flags(),std::cout.precision(),final,transfer);
    const bool diagnostics=d.input.fluid.run.diagnosticSeries;
    const auto reconstructionTemplates=diagnostics?reconstructionTemplateRow(fluid):std::string{};
    const auto physical=physicalRows(mesh,fluid,charge,step,receipt,qremap,physicalChargeTransport,physicalEndField,d.charge);
    const auto electricPreparation=diagnostics?electricPreparationRows(fluid,step,physicalChargeTransport,physicalEndField):std::string{};
    const auto regularPreparation=diagnostics?regularPreparationReuseRow(fluid,receipt):std::string{};
    const auto chargeTransportStages=diagnostics?chargeTransportStageRows(fluid,step,physicalChargeTransport):std::string{};
    const auto openConvectionStages=diagnostics?openConvectionStageRows(fluid,step):std::string{};
    const auto amrPerformance=diagnostics?amrPerformanceRow(fluid,field,receipt,qremap):std::string{};
    const auto solverStages=diagnostics?solverStageRows(fluid,step):std::string{};
    const auto localPlicExecution=diagnostics?localPlicExecutionRow(fluid,step):std::string{};
    const auto frozenMatrixPreparation=diagnostics?frozenMatrixPreparationRow(fluid,step):std::string{};
    const auto maxwellStages=diagnostics?maxwellStageRows(fluid,step):std::string{};
    const auto viscosityProductStages=diagnostics?viscosityProductStageRows(fluid,step):std::string{};
    const auto csgBoundaryRows=diagnostics?csgBoundaryPreparationRows(fluid,step,csgBoundaryPreparation):std::string{};
    const auto stressGeometryReuse=diagnostics?stressGeometryReuseRow(fluid,step):std::string{};
    const auto& q=charge.integratedCharge();const auto& pub=fluid.phaseGeometry().vofPublication();
    const double current=boundaryCurrent(field);
    require(q.size()==pub.cells().size()&&field.potential().size()==q.size()&&field.epsilonByCell().size()==q.size(),"electric native dimensions differ");
    const auto& run=d.input.fluid.run;
    const bool periodic=run.writeEverySteps>0?t.step%run.writeEverySteps==0:
        run.writeInterval>0&&std::floor(t.time()/run.writeInterval+1e-12)>std::floor(d.time/run.writeInterval+1e-12);
    const bool emit=run.writeFields&&(t.step==0?run.writeInitial:(final&&run.writeFinal)||(t.step==1&&run.writeFirst)||periodic);
    TecplotWriter::Zone volume,surface,faces;
    Data::Frame frame{t.step,t.time(),pub.identity().value(),charge.identity(),"","","",field.isElectrostaticSnapshot()?"electrostaticSnapshot":"ohmicEnd"};
    if(emit){
        volume=Fluid::vofVolumeZone(mesh,fluid,step,transfer);surface=Fluid::vofInterfaceZone(fluid,step,d.input.fluid.radial,transfer);
        faces=electricFaces(mesh,field);
        for(const char* name:{"integrated_charge","charge_density","electric_potential","permittivity"})volume.fields.push_back({name,TecplotWriter::Location::CellCentered,{}});
        const auto offset=volume.fields.size()-4;
        for(std::size_t row=0;row<volume.connectivity.size();++row){
            const int native=fluid.values().layout().nativeCell(static_cast<int>(row));
            require(native>=0&&static_cast<std::size_t>(native)<q.size()&&volume.fields[0].values[row]==native,"volume row is not actual native cell");
            const auto c=pub.cells()[native].capacity;require(c>0,"native charge capacity is not positive");
            const double density=static_cast<double>(static_cast<long double>(q[native])/c);
            require(std::isfinite(density),"charge density diagnostic is nonfinite");
            const std::array<double,4> values{{q[native],density,field.potential()[native],field.epsilonByCell()[native]}};
            for(std::size_t k=0;k<4;++k)volume.fields[offset+k].values.push_back(values[k]);
        }
        electricMetadata(volume,charge,field);electricMetadata(faces,charge,field);
        if(!surface.points.empty())electricMetadata(surface,charge,field);
        std::ostringstream number;number<<std::setw(6)<<std::setfill('0')<<t.step;const auto suffix=number.str();
        frame.volume="volume_"+suffix+".dat";frame.faces="electric_faces_"+suffix+".dat";
        if(!surface.points.empty())frame.surface="interface_"+suffix+".dat";
    }
    // Every owner/schema/derived value above is checked before file mutation.
    // On I/O failure the complete frame prefix remains authoritative.
    const auto ioStart=Clock::now();d.step=t.step;d.time=t.time();d.geometry=frame.geometry;d.charge=frame.charge;
    try{
        if(emit){TecplotWriter::write(d.root/"fields"/frame.volume,volume);if(!frame.surface.empty())TecplotWriter::write(d.root/"fields"/frame.surface,surface);
            TecplotWriter::write(d.root/"fields"/frame.faces,faces);}
        d.table("conservation",rows.conservation);d.table("solver",rows.solver);
        d.table("energy",rows.energyPrefix+"NA"+rows.energySuffix);
        if(diagnostics)d.table("vof_native_velocity",Fluid::vofNativeVelocityRow(mesh,fluid,step,transfer));
        d.table("reconstruction_templates",reconstructionTemplates);
        if(!rows.capillary.empty())d.table(fs::path(rows.capillaryFile).stem().string(),rows.capillary);
        long double sum=0;for(double value:q)sum+=value;
        const auto& er=field.report();const auto sourceClock=field.chargeSourceTime();auto row=stream();
        row<<t.step<<','<<t.time()<<','<<pub.identity().value()<<','<<charge.identity()<<','<<LAMP::IO::csv(frame.kind)<<','<<field.effectiveDt()<<','
            <<sourceClock.step<<','<<sourceClock.begin<<','<<sourceClock.end<<','<<static_cast<int>(sourceClock.position)<<','<<sum<<','
            <<er.relativeLinearResidual<<','<<er.maximumRelativeGaussResidual<<','<<current<<'\n';d.table("electric",row.str());
        row=stream();row<<t.step<<','<<t.time()<<',';
        if(step&&step->frozenMidpoint().hasMaxwellLoad()){
            const auto& load=*step->frozenMidpoint().forces().maxwellLoad;const auto& actual=step->report().stages.back();
            row<<load.phaseGeometry().timeLayer().time()<<','<<load.phaseGeometry().vofIdentity().value()<<','<<actual.maxwellWork;
        }else row<<"NA,NA,NA";
        row<<',';if(step)row<<step->candidateState().phaseGeometry().vofIdentity().value();else row<<"NA";
        row<<','<<pub.identity().value()<<','<<(transfer?1:0)<<'\n';d.table("maxwell_work",row.str());
        if(receipt){
            const auto& r=*receipt->decision;row=stream();row<<t.step<<','<<t.time()<<','<<receipt->physicalEnd->phaseGeometry().vofIdentity().value()<<','
                <<pub.identity().value()<<','<<r.sourceLeaves<<','<<r.targetLeaves<<','<<r.requested<<','<<r.changed<<','<<r.clones<<','<<r.sweeps<<','
                <<r.totalSeconds<<','<<receipt->inventoryTransferSeconds<<','<<receipt->regridSeconds<<','
                <<d.input.chargeDensityJumpTolerance.has_value()<<','<<r.additionalRefineCells<<','<<r.additionalHoldCells<<'\n';d.table("regrid_decision",row.str());
            if(transfer){const auto& r=transfer->report();const auto& qr=qremap->report();const auto& hr=receipt->fluidTransfer->report();
                row=stream();row<<t.step<<','<<t.time()<<','<<transfer->sourceIdentity().value()<<','<<pub.identity().value()<<','
                    <<qremap->sourceCharge().identity()<<','<<charge.identity()<<','<<transfer->nativeOverlaps().bySource().size()<<','
                    <<transfer->nativeOverlaps().byTarget().size()<<','<<r.refinedDonors<<','<<r.coarsenedTargets<<','
                    <<r.sourcePhaseInventory[0]<<','<<r.targetPhaseInventory[0]<<','<<r.sourcePhaseInventory[1]<<','<<r.targetPhaseInventory[1]<<','
                    <<qr.sourceTotalCharge<<','<<qr.candidateTotalCharge<<','<<qr.chargeResidual<<','<<qr.relativeChargeResidual<<','
                    <<hr.kineticChange<<','<<hr.historyReset<<','<<qr.limitedLinearRemap<<','<<qr.highOrderRefinedDonors<<','
                    <<qr.limitedRefinedDonors<<','<<qr.finiteEligibilityP0Donors<<','<<qr.maximumLinearQuantizationAdjustment<<','
                    <<qr.linearQuantizationL1Adjustment<<','<<qr.totalSeconds<<'\n';d.table("regrid",row.str());}
        }
        // All optional physical rows were authenticated and rendered before
        // the first payload mutation. They share the original flush/index gate.
        d.table("charge_balance",physical.charge);d.table("phase_boundary",physical.phase);
        d.table("vof_boundary_work",physical.work);d.table("regrid_projection_work",physical.projection);
        d.table("performance_amr",amrPerformance);
        d.table("solver_stages",solverStages);
        d.table("local_plic_execution",localPlicExecution);
        d.table("maxwell_stages",maxwellStages);
        d.table("viscosity_product_stages",viscosityProductStages);
        d.table("csg_boundary_preparation",csgBoundaryRows);
        d.table("electric_preparation",electricPreparation);
        d.table("regular_preparation_reuse",regularPreparation);
        d.table("charge_transport_stages",chargeTransportStages);
        d.table("open_convection_stages",openConvectionStages);
        d.table("frozen_matrix_preparation",frozenMatrixPreparation);
        d.table("stress_geometry_reuse",stressGeometryReuse);
        // Performance uses elapsed cost through the start of its own serialization.
        // Every current history must succeed before publishing a frame index.
        if(diagnostics){
            const auto io=std::chrono::duration<double>(Clock::now()-ioStart).count();
            std::istringstream performance(Fluid::vofPerformance(mesh,fluid,step,compute,io,std::chrono::duration<double>(Clock::now()-start).count(),transfer));
            unsigned concurrentStages=0;if(step)for(const auto& stage:step->report().stages)concurrentStages+=stage.timing.concurrentViscosityCapillary;
            std::string line;while(std::getline(performance,line)){if(line.empty())continue;const auto a=line.find(','),b=line.find(',',a+1),c=line.find(',',b+1);
                require(a!=std::string::npos&&b!=std::string::npos&&c!=std::string::npos,"performance row is malformed");
                const auto name=line.substr(b+1,c-b-1);require(!name.empty()&&std::all_of(name.begin(),name.end(),[](char x){return std::isalnum(static_cast<unsigned char>(x))||x=='_';}),"performance stage is not a fixed identifier");
                const auto file="performance_"+name;d.names[file]=performanceHeader;d.table(file,line+','+std::to_string(concurrentStages)+'\n');}
        }
        d.flush();
        if(emit){auto frames=d.frames;frames.push_back(frame);replace(d.root/"meta/frames.csv",d.index(frames));d.frames.swap(frames);}
        d.totals=rows.next;d.finalSeen=final;d.manifest();
        if(!rows.progress.empty())std::cout<<rows.progress<<std::endl;
        d.acceptedClock=t;
        d.rememberCompletedOutput(mesh,fluid,step);
    }catch(...){d.ioFailed=true;throw;}
}
void VofElectrohydrodynamicsOutput::accepted(const QuadtreeMesh& mesh,const H::HydroTwoPhaseState& fluid,
    const E::ElectricVofSurfaceCharge& charge,const E::ElectricVofSharpField& field,const H::HydroVofStep* step,
    Clock::time_point start,double compute,bool final){
    auto& d=*data_;require(d.begun&&!d.complete&&!d.ioFailed&&d.input.method==VofElectrohydrodynamicsInput::Method::VofSharpSurface,"B1 writer is not active");
    fluid.validateFor(mesh);const auto t=fluid.timeLayer();const auto& pub=fluid.phaseGeometry().vofPublication();
    const auto& g=field.geometry();g.validateFor(mesh,pub);charge.validateFor(mesh,g,t);
    require(sameClock(t,g.sourceTime())&&sameClock(t,field.candidateTime())&&field.chargeCandidate().sameInstance(charge),"B1 accepted owner/complete clock differs");
    require(t.position==EhdClock::Position::Initial||t.position==EhdClock::Position::End,"B1 output is not Initial/End");
    if(step){require(step->candidateState().sameIdentity(fluid)&&field.advection()&&field.effectiveDt()>0,"B1 End lacks actual final joint receipt");
        const auto* transport=field.advection()->interfaceTransport();require(transport&&transport->sameInstance(step->transport())
            &&field.sourceCharge().sameInstance(field.advection()->sourceCharge()),"B1 End field is not from the actual final transport");
        require(step->frozenMidpoint().hasMaxwellLoad()&&step->frozenMidpoint().forces().maxwellLoad->isSharp(),"B1 step lost its actual typed sharp force");
        require(step->frozenMidpoint().forces().maxwellLoad->sharpStress().field().sourceCharge().sameInstance(field.sourceCharge()),
            "B1 actual Midpoint force and End diagnostic have different accepted charge owners");
        require(step->report().maxwellCorrectionChecked&&step->report().stages.size()==2,"B1 step lacks actual two-stage electric correction");
    }else require(t.position==EhdClock::Position::Initial&&field.effectiveDt()==0&&field.sourceCharge().sameInstance(charge),"B1 Initial requires actual dt0 field");
    require(step?(t.step==d.step+1&&t.begin==d.time&&t.time()>d.time):(d.step==-1&&t.step==0&&t.time()==0),"B1 output requires sequential accepted states");
    require(!final||t.time()==d.input.fluid.run.endTime,"B1 final flag differs from actual endTime");
    const auto rows=Fluid::vofRows(mesh,fluid,step,d.input.fluid,d.totals,compute,std::cout.flags(),std::cout.precision(),final,nullptr);
    const bool diagnostics=d.input.fluid.run.diagnosticSeries;
    const auto reconstructionTemplates=diagnostics?reconstructionTemplateRow(fluid):std::string{};
    const auto solverStages=diagnostics?solverStageRows(fluid,step):std::string{};
    const auto localPlicExecution=diagnostics?localPlicExecutionRow(fluid,step):std::string{};
    const auto frozenMatrixPreparation=diagnostics?frozenMatrixPreparationRow(fluid,step):std::string{};
    const auto maxwellStages=diagnostics?maxwellStageRows(fluid,step):std::string{};
    const auto viscosityProductStages=diagnostics?viscosityProductStageRows(fluid,step):std::string{};
    const auto stressGeometryReuse=diagnostics?stressGeometryReuseRow(fluid,step):std::string{};
    const auto& run=d.input.fluid.run;const bool periodic=run.writeEverySteps>0?t.step%run.writeEverySteps==0:
        run.writeInterval>0&&std::floor(t.time()/run.writeInterval+1e-12)>std::floor(d.time/run.writeInterval+1e-12);
    const bool emit=run.writeFields&&(t.step==0?run.writeInitial:(final&&run.writeFinal)||(t.step==1&&run.writeFirst)||periodic);
    Data::Frame frame{t.step,t.time(),pub.identity().value(),0,"","","",field.effectiveDt()==0?"electrostaticSnapshot":"ohmicEnd"};
    TecplotWriter::Zone volume,surface,electric;
    if(emit){volume=Fluid::vofVolumeZone(mesh,fluid,step,nullptr);surface=Fluid::vofInterfaceZone(fluid,step,false,nullptr);
        electric=sharpInterface(charge,field);sharpMetadata(volume,field);if(!surface.points.empty())sharpMetadata(surface,field);
        std::ostringstream number;number<<std::setw(6)<<std::setfill('0')<<t.step;const auto suffix=number.str();
        frame.volume="volume_"+suffix+".dat";frame.faces="electric_interface_"+suffix+".dat";
        if(!surface.points.empty())frame.surface="interface_"+suffix+".dat";
    }
    const auto ioStart=Clock::now();d.step=t.step;d.time=t.time();d.geometry=frame.geometry;d.charge=0;
    try{
        if(emit){TecplotWriter::write(d.root/"fields"/frame.volume,volume);if(!frame.surface.empty())TecplotWriter::write(d.root/"fields"/frame.surface,surface);
            TecplotWriter::write(d.root/"fields"/frame.faces,electric);}
        d.table("conservation",rows.conservation);d.table("solver",rows.solver);d.table("energy",rows.energyPrefix+"NA"+rows.energySuffix);
        if(diagnostics)d.table("vof_native_velocity",Fluid::vofNativeVelocityRow(mesh,fluid,step,nullptr));
        d.table("reconstruction_templates",reconstructionTemplates);
        d.table("solver_stages",solverStages);
        d.table("local_plic_execution",localPlicExecution);
        d.table("maxwell_stages",maxwellStages);
        d.table("viscosity_product_stages",viscosityProductStages);
        d.table("frozen_matrix_preparation",frozenMatrixPreparation);
        d.table("stress_geometry_reuse",stressGeometryReuse);
        if(!rows.capillary.empty())d.table(fs::path(rows.capillaryFile).stem().string(),rows.capillary);
        long double sum=0;for(double value:charge.integratedCharge())sum+=value;
        const auto& er=field.report();const auto s=field.sourceCharge().timeLayer();auto row=stream();
        row<<t.step<<','<<t.time()<<','<<pub.identity().value()<<",NA,"<<LAMP::IO::csv(frame.kind)<<','<<field.effectiveDt()<<','
            <<s.step<<','<<s.begin<<','<<s.end<<','<<static_cast<int>(s.position)<<','<<sum<<','<<er.relativeLinearResidual<<','
            <<er.maximumRelativeGaussResidual<<','<<er.boundaryCurrent<<','<<er.relativeGlobalGaussResidual<<','<<er.relativeGlobalChargeUpdateResidual<<','
            <<er.skeletonDofs<<','<<er.nonzeros<<','<<er.assembleSeconds<<','<<er.factorSeconds<<','<<er.solveSeconds<<','<<er.recoverySeconds;
        // Only the retained final AC Midpoint and true End receipts exist
        // here. Nested Q1 timers are inclusive, not disjoint stage totals.
        const E::ElectricVofSharpField* midpoint=nullptr;
        if(step&&step->frozenMidpoint().hasMaxwellLoad())
            midpoint=&step->frozenMidpoint().forces().maxwellLoad->sharpStress().field();
        const auto appendQ1=[&](const E::ElectricSurfaceChargeTransport* receipt){
            if(!receipt){for(unsigned i=0;i<12;++i)row<<",NA";return;}
            const auto& r=receipt->report();row<<','<<r.prepareSeconds<<','<<r.materialIntegrationSeconds<<','<<r.pointProjectionSeconds
                <<','<<r.nativeArcProjectionSeconds<<','<<r.coreArcProjectionSeconds<<','<<r.identityArcProjectionSeconds
                <<','<<r.materialSamples<<','<<r.endpointProjectionRequests<<','<<r.endpointProjectionMemoHits
                <<','<<r.nativeArcProjections<<','<<r.exactArcProjections<<','<<r.rootIsolationCalls;};
        appendQ1(midpoint?midpoint->advection():nullptr);appendQ1(field.advection());
        if(midpoint){const auto& r=midpoint->report();row<<','<<r.assembleSeconds<<','<<r.factorSeconds<<','<<r.solveSeconds<<','<<r.recoverySeconds;}
        else for(unsigned i=0;i<4;++i)row<<",NA";
        const auto appendCoreReasons=[&](const E::ElectricSurfaceChargeTransport* receipt){
            if(!receipt){for(unsigned i=0;i<8;++i)row<<",NA";return;}
            const auto& r=receipt->report();row<<','<<r.coreHessianFallbacks<<','<<r.coreEndpointSignFallbacks
                <<','<<r.coreNonfiniteFallbacks<<','<<r.coreRootWidthFallbacks<<','<<r.coreForcedDistanceFallbacks
                <<','<<r.coreEndpointProjections<<','<<r.coreInteriorProjections<<','<<r.geometryBasePreparations;};
        appendCoreReasons(midpoint?midpoint->advection():nullptr);appendCoreReasons(field.advection());
        row<<'\n';d.table("electric",row.str());
        if(field.advection()){const auto& r=field.advection()->report();row=stream();row<<t.step<<','<<t.time()<<','<<r.sourceCharge<<','<<r.targetCharge<<','
            <<r.relativeChargeResidual<<',';
            if(r.integralDiagnosticsAvailable)row<<r.sourceMomentQuadratureResidual;else row<<"NA";
            row<<','<<r.maximumDisplacement<<','<<r.maximumProjectionDistance<<','<<r.materialSamples<<','<<r.parameterSplits<<','
                <<(r.mode==E::ElectricSurfaceChargeTransferMode::ExactPartitionIdentity)<<','<<r.integralDiagnosticsAvailable<<','
                <<r.identityArcTransfers<<','<<r.sourceDensityReconstructions<<'\n';d.table("surface_transport",row.str());}
        row=stream();row<<t.step<<','<<t.time()<<',';
        if(step&&step->frozenMidpoint().hasMaxwellLoad()){
            const auto& load=*step->frozenMidpoint().forces().maxwellLoad;require(load.isSharp(),"B1 step retained a different electric stress model");
            const auto& actual=step->report().stages.back();row<<load.phaseGeometry().timeLayer().time()<<','<<load.phaseGeometry().vofIdentity().value()<<','
                <<actual.maxwellWork;
        }else row<<"NA,NA,NA";
        row<<',';if(step)row<<step->candidateState().phaseGeometry().vofIdentity().value();else row<<"NA";
        row<<','<<pub.identity().value()<<",0\n";d.table("maxwell_work",row.str());
        if(step){row=stream();const auto& r=step->report();row<<t.step<<','<<t.time()<<','<<r.maxwellCorrectionChecked<<','<<r.stages.size()<<','
            <<r.maxwellLoadChange<<','<<r.maxwellWorkChange<<'\n';d.table("electric_correction",row.str());}
        if(diagnostics){
            const auto io=std::chrono::duration<double>(Clock::now()-ioStart).count();
            std::istringstream performance(Fluid::vofPerformance(mesh,fluid,step,compute,io,std::chrono::duration<double>(Clock::now()-start).count(),nullptr));
            unsigned concurrentStages=0;if(step)for(const auto& stage:step->report().stages)concurrentStages+=stage.timing.concurrentViscosityCapillary;
            std::string line;while(std::getline(performance,line)){if(line.empty())continue;const auto a=line.find(','),b=line.find(',',a+1),c=line.find(',',b+1);
                require(a!=std::string::npos&&b!=std::string::npos&&c!=std::string::npos,"B1 performance row malformed");
                const auto name=line.substr(b+1,c-b-1);require(!name.empty()&&std::all_of(name.begin(),name.end(),[](char x){return std::isalnum(static_cast<unsigned char>(x))||x=='_';}),"B1 performance stage invalid");
                d.names["performance_"+name]=performanceHeader;d.table("performance_"+name,line+','+std::to_string(concurrentStages)+'\n');}
        }
        d.flush();if(emit){auto frames=d.frames;frames.push_back(frame);replace(d.root/"meta/frames.csv",d.index(frames));d.frames.swap(frames);}
        d.totals=rows.next;d.finalSeen=final;d.manifest();if(!rows.progress.empty())std::cout<<rows.progress<<std::endl;
        d.rememberCompletedOutput(mesh,fluid,step);
    }catch(...){d.ioFailed=true;throw;}
}
void VofElectrohydrodynamicsOutput::stepCompleted(const EhdClock& t,double wallSeconds){
    PassivePerformanceEnvironment environment;
    auto& d=*data_;
    require(d.begun&&!d.complete&&!d.ioFailed&&t.position==EhdClock::Position::End&&d.pendingCompletion.has_value(),
        "completed-step performance lacks successful accepted output");
    const auto snapshot=*d.pendingCompletion;
    require(sameClock(t,snapshot.clock)&&t.step==d.step&&t.time()==d.time&&snapshot.geometry==d.geometry
        &&t.step>d.completedPerformanceStep,"completed-step performance has a stale or duplicate accepted owner/clock");
    require(std::isfinite(wallSeconds)&&wallSeconds>=0,"invalid completed-step duration");
    if(!d.input.fluid.run.diagnosticSeries){
        d.completedPerformanceStep=t.step;d.pendingCompletion.reset();return;
    }
    auto row=stream();row<<t.step<<','<<t.time()<<','<<t.begin<<','<<t.end<<','<<static_cast<int>(t.position)<<','
        <<snapshot.geometry<<','<<snapshot.actualAc<<','<<wallSeconds<<','<<snapshot.cells<<','
        <<snapshot.velocityDofs<<','<<snapshot.pressureDofs<<'\n';
    try{
        d.table("performance_step_complete",row.str());
        // All other histories already passed accepted()'s original flush gate.
        d.histories.at("performance_step_complete")->flush();
        d.completedPerformanceStep=t.step;d.pendingCompletion.reset();
    }catch(...){d.ioFailed=true;throw;}
}
void VofElectrohydrodynamicsOutput::timeStep(const HydrodynamicsTimeStepReport& r,const EhdClock& t,int retries){
    auto& d=*data_;require(d.begun&&!d.complete&&!d.ioFailed&&d.input.fluid.run.adaptive&&t.position==EhdClock::Position::End
        &&t.step==d.step+1&&t.begin==d.time&&retries>=0,"adaptive output lacks next committed clock");
    auto row=stream();row<<t.step<<','<<t.time()<<','<<t.begin<<','<<t.end-t.begin<<','<<r.proposed_dt<<',';
    const auto add=[&](long double x,bool available){if(available&&std::isfinite(x))row<<x;else row<<"NA";row<<',';};
    add(r.flow_rate,r.flow_rate_available);add(r.capillary_rate,r.capillary_rate_available);add(r.interface_minimum_spacing,r.interface_spacing_available);
    add(r.flow_bound,r.flow_bound_available);add(r.capillary_bound,r.capillary_bound_available);add(r.growth_bound,r.growth_bound_available);add(r.max_bound,r.max_bound_available);
    row<<LAMP::IO::csv(r.limiter)<<','<<retries<<','<<(r.endpoint_aligned&&t.end==r.end)<<','
        <<(t.end==d.input.fluid.run.endTime&&t.end-t.begin<d.input.fluid.run.adaptive->minDt)<<','<<(r.minimum_clock_rounded&&t.end==r.end)<<'\n';
    d.table("time_step",row.str());d.flush();
}
void VofElectrohydrodynamicsOutput::retry(std::int64_t attempt,std::int64_t accepted,double time,double oldDt,double nextDt,const std::string& reason){
    auto& d=*data_;require(d.begun&&!d.complete&&!d.ioFailed&&accepted==d.step&&time==d.time,"retry differs from accepted prefix");
    auto row=stream();row<<"retry,"<<attempt<<','<<accepted<<','<<time<<','<<oldDt<<','<<nextDt<<",VOF_EHD,"<<LAMP::IO::csv(reason)<<'\n';d.table("events",row.str());d.flush();
}
void VofElectrohydrodynamicsOutput::failure(std::int64_t attempt,std::int64_t accepted,double time,const std::string& stage,const std::string& reason){
    auto& d=*data_;if(!d.begun)return;d.status="FAILED";
    // A committed clock can advance before timeStep I/O fails. Old output IDs
    // then have no available witness for that new accepted numerical state.
    if(accepted!=d.step||time!=d.time){d.geometry=0;d.charge=0;}
    d.step=accepted;d.time=time;
    auto f=stream();f<<"{\"app\":"<<quoted(app)<<",\"status\":\"FAILED\",\"attempted_step\":"<<attempt<<",\"accepted_step\":"<<accepted
        <<",\"accepted_time\":"<<time<<",\"stage\":"<<quoted(stage)<<",\"reason\":"<<quoted(reason)
        <<",\"accepted_state_preserved\":true,\"numerical_state_rolled_back\":false}\n";
    write(d.root/"meta/failure.json",f.str());d.manifest();
    auto row=stream();row<<"failure,"<<attempt<<','<<accepted<<','<<time<<",NA,NA,"<<LAMP::IO::csv(stage)<<','<<LAMP::IO::csv(reason)<<'\n';d.table("events",row.str());d.flush();
}
void VofElectrohydrodynamicsOutput::complete(){
    auto& d=*data_;require(d.begun&&!d.complete&&!d.ioFailed&&d.status=="RUNNING"&&d.finalSeen&&d.time==d.input.fluid.run.endTime,
        "completion lacks true final accepted clock or successful output");
    for(auto& h:d.histories)if(h.second)h.second->close();d.status="COMPLETED";d.manifest();d.complete=true;
}
