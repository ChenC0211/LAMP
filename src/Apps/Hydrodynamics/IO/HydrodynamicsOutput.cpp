#include "Hydrodynamics/IO/HydrodynamicsOutput.h"
#include "Hydrodynamics/HydrodynamicsRun.h"
#include "Hydrodynamics/HydrodynamicsVofRegrid.h"
#include "Modules/Hydrodynamics/HydroVofMeshTransfer.h"
#include "Modules/Hydrodynamics/HydroVofStep.h"
#include "Hydrodynamics/HydrodynamicsTimeStep.h"
#include "IO/HydroOutputComponents.h"
#include "Core/IO/OutputProvenance.h"
#include "Core/IO/TecplotWriter.h"
#include "Core/IO/BinaryArchive.h"
#include "Core/Mesh/QuadtreeMesh.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <locale>
#include <limits>
#include <map>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <system_error>
#include <thread>
#ifdef _WIN32
#include <winerror.h>
#endif
namespace {
using namespace LAMP::Hydrodynamics;using namespace LAMP::Interface;

using LAMP::IO::stream;using LAMP::IO::quoted;
namespace fs=std::filesystem;namespace Fluid=AppIO::HydroOutput;
void write(const fs::path& path,const std::string& text,bool append=false){
    std::ofstream out;out.exceptions(std::ios::failbit|std::ios::badbit);
    try{out.open(path,std::ios::binary|(append?std::ios::app:std::ios::trunc));out<<text;out.close();}
    catch(const std::exception& e){throw std::runtime_error("Hydro output "+path.generic_string()+": "+e.what());}
}
bool temporaryIndexReplacementError(const std::error_code& error){
#ifdef _WIN32
    return error.category()==std::system_category() &&
        (error.value()==ERROR_ACCESS_DENIED || error.value()==ERROR_SHARING_VIOLATION || error.value()==ERROR_LOCK_VIOLATION);
#else
    (void)error;return false;
#endif
}
void replace(const fs::path& path,const std::string& text){
    const fs::path temp=path.string()+".tmp";write(temp,text);
    std::error_code error;fs::rename(temp,path,error);
    // A reader may briefly deny delete/rename sharing. Keep the complete temp
    // and old index intact; only these named Windows errors get bounded waits.
    constexpr std::array<int,6> delays{{5,10,20,40,80,160}};
    for(std::size_t retry=0;error&&temporaryIndexReplacementError(error)&&retry<delays.size();++retry){
        std::this_thread::sleep_for(std::chrono::milliseconds(delays[retry]));fs::rename(temp,path,error);
    }
    if(error)throw std::runtime_error("Hydro output index replacement "+path.generic_string()+": "+error.message()+" (code="+std::to_string(error.value())+")");
}
bool sameReceiptClock(const LAMP::TimeLayer& a,const LAMP::TimeLayer& b){
    return a.step==b.step&&a.begin==b.begin&&a.end==b.end&&a.position==b.position;
}
const GeometricVofMeshTransfer* validateReceipt(const QuadtreeMesh& mesh,const HydroTwoPhaseState& state,
    const HydroVofStep* step,const HydrodynamicsVofRemapReceipt* receipt){
    state.validateFor(mesh);if(!receipt)return nullptr;
    const auto invalid=[](){throw std::invalid_argument("Hydro VOF output receipt has incompatible physical source, target or accepted End clock");};
    if(!step||!receipt->sourceMesh||!receipt->physicalEnd||receipt->physicalStep!=step||!receipt->decision)invalid();
    const auto& physical=*receipt->physicalEnd;const auto& clock=state.timeLayer();const auto& decision=*receipt->decision;
    physical.validateFor(*receipt->sourceMesh);step->candidateState().validateFor(*receipt->sourceMesh);
    if(!physical.sameIdentity(step->candidateState())||!sameReceiptClock(clock,physical.timeLayer())||
        !sameReceiptClock(clock,decision.sourceClock)||clock.position!=LAMP::TimeLayer::Position::End||
        decision.sourceLeaves!=static_cast<std::size_t>(receipt->sourceMesh->cellIndexNum())||
        decision.targetLeaves!=static_cast<std::size_t>(mesh.cellIndexNum()))invalid();
    if(!std::isfinite(receipt->inventoryTransferSeconds)||receipt->inventoryTransferSeconds<0||
        !std::isfinite(receipt->regridSeconds)||receipt->regridSeconds<0)invalid();
    if(static_cast<bool>(receipt->fluidTransfer)!=static_cast<bool>(receipt->interfaceTransfer))invalid();
    if(!receipt->fluidTransfer){
        if(decision.changed||!physical.sameIdentity(state)||receipt->inventoryTransferSeconds!=0)invalid();
        return nullptr;
    }
    const auto& fluid=*receipt->fluidTransfer;const auto& interfaceTransfer=*receipt->interfaceTransfer;
    if(!decision.changed||!fluid.sourceState().sameIdentity(physical)||!fluid.candidateState().sameIdentity(state)||
        &fluid.interfaceTransfer()!=receipt->interfaceTransfer||!interfaceTransfer.report().localPlic||
        interfaceTransfer.sourceIdentity()!=physical.phaseGeometry().vofIdentity()||
        interfaceTransfer.candidate().identity()!=state.phaseGeometry().vofIdentity()||
        !sameReceiptClock(clock,interfaceTransfer.candidate().timeLayer()))invalid();
    return &interfaceTransfer;
}
struct DiagnosticRecord{std::string header,row;};
DiagnosticRecord decisionRecord(const HydroTwoPhaseState& state,const HydrodynamicsVofRemapReceipt& receipt){
    auto header=stream(),row=stream();bool first=true;
    const auto add=[&](const char* name,auto value){if(!first){header<<',';row<<',';}first=false;header<<name;row<<value;};
    const auto& d=*receipt.decision;const auto& t=state.timeLayer();
    add("step",t.step);add("time",t.time());add("begin",t.begin);add("end",t.end);add("same_end",1);
    add("source_geometry_identity",receipt.physicalEnd->phaseGeometry().vofIdentity().value());
    add("target_geometry_identity",state.phaseGeometry().vofIdentity().value());
    add("source_leaves",d.sourceLeaves);add("target_leaves",d.targetLeaves);add("mixed_cells",d.mixedCells);
    add("geometric_mixed_cells",d.geometricMixedCells);add("unresolved_mixed_cells",d.unresolvedMixedCells);
    add("pure_jump_faces",d.pureJumpFaces);add("interface_seed_cells",d.interfaceSeedCells);add("interface_band_cells",d.interfaceBandCells);
    add("flow_refine_cells",d.flowRefineCells);add("flow_coarsen_cells",d.flowCoarsenCells);add("maximum_h_abs_vorticity",d.maximumHAbsVorticity);
    add("additional_refine_cells",d.additionalRefineCells);add("additional_hold_cells",d.additionalHoldCells);
    add("requested",d.requested);add("changed",d.changed);add("clones",d.clones);add("sweeps",d.sweeps);
    add("indicator_seconds",d.indicatorSeconds);add("planning_seconds",d.planningSeconds);add("adaptation_seconds",d.adaptationSeconds);
    add("decision_seconds",d.totalSeconds);add("inventory_transfer_seconds",receipt.inventoryTransferSeconds);add("regrid_seconds",receipt.regridSeconds);
    header<<'\n';row<<'\n';return {header.str(),row.str()};
}
DiagnosticRecord regridRecord(const HydroTwoPhaseState& state,const HydrodynamicsVofRemapReceipt& receipt){
    auto header=stream(),row=stream();bool first=true;
    const auto add=[&](const char* name,auto value){if(!first){header<<',';row<<',';}first=false;header<<name;row<<value;};
    const auto& transfer=*receipt.interfaceTransfer;const auto& i=transfer.report();const auto& h=receipt.fluidTransfer->report();const auto& t=state.timeLayer();
    add("step",t.step);add("time",t.time());add("begin",t.begin);add("end",t.end);add("same_end",1);
    add("source_geometry_identity",transfer.sourceIdentity().value());add("target_geometry_identity",transfer.candidate().identity().value());
    add("source_leaves",transfer.nativeOverlaps().bySource().size());add("target_leaves",transfer.nativeOverlaps().byTarget().size());
    add("refined_donors",i.refinedDonors);add("refined_targets",i.refinedTargets);add("coarsened_targets",i.coarsenedTargets);add("unchanged_targets",i.unchangedTargets);
    add("resolved_refined_donors",i.resolvedRefinedDonors);add("unresolved_refined_donors",i.unresolvedRefinedDonors);add("p0_refined_pieces",i.p0RefinedPieces);
    add("corrected_refined_donors",i.correctedRefinedDonors);add("roundoff_closed_phases",i.roundoffClosedPhases);
    add("source_liquid",i.sourcePhaseInventory[0]);add("target_liquid",i.targetPhaseInventory[0]);
    add("source_gas",i.sourcePhaseInventory[1]);add("target_gas",i.targetPhaseInventory[1]);
    add("donor_relative_liquid",i.maximumRelativeDonorPhaseResidual[0]);add("donor_relative_gas",i.maximumRelativeDonorPhaseResidual[1]);
    add("global_relative_liquid",i.relativeGlobalPhaseResidual[0]);add("global_relative_gas",i.relativeGlobalPhaseResidual[1]);
    add("raw_inventory_relative",i.maximumRelativeRawInventoryResidual);add("raw_geometry_relative",i.maximumRelativeLocalGeometryResidual);
    add("target_phase_closure",i.maximumRelativeTargetPhaseClosure);add("history_reset",h.historyReset);add("remap_solves",h.remapSolves);
    add("source_native_P_x",h.source.momentum[0]);add("target_native_P_x",h.target.momentum[0]);
    add("source_native_P_y",h.source.momentum[1]);add("target_native_P_y",h.target.momentum[1]);
    add("delta_native_P_x",h.momentumChange[0]);add("delta_native_P_y",h.momentumChange[1]);
    add("source_native_kinetic",h.source.kinetic);add("target_native_kinetic",h.target.kinetic);add("delta_native_kinetic",h.kineticChange);
    add("projection_delta_P_x",h.projectionMomentumChange[0]);add("projection_delta_P_y",h.projectionMomentumChange[1]);
    add("projection_delta_kinetic",h.projectionKineticChange);add("projection_correction_energy",h.projectionCorrectionEnergy);
    add("projection_energy_identity_residual",h.projectionEnergyIdentityResidual);
    add("raw_relative_divergence",h.maximumRawRelativeDivergence);add("final_relative_divergence",h.maximumFinalRelativeDivergence);
    add("raw_divergence_L1",h.rawDivergenceL1);add("final_divergence_L1",h.finalDivergenceL1);
    add("pressure_gauge_shift",h.pressureGaugeShift);add("pressure_gauge_relative_residual",h.pressureGaugeRelativeResidual);
    add("source_leaf_candidates",h.sourceLeafCandidates);add("face_spans",h.faceSpans);add("polynomial_integrals",h.polynomialIntegrals);
    add("wall_normals",h.wallNormals);add("axis_normals",h.axisNormals);
    add("decision_seconds",receipt.decision->totalSeconds);add("inventory_transfer_seconds",receipt.inventoryTransferSeconds);
    add("raw_flux_seconds",h.rawFluxSeconds);add("projection_seconds",h.projectionSeconds);add("fluid_remap_seconds",h.totalSeconds);
    add("regrid_seconds",receipt.regridSeconds);header<<'\n';row<<'\n';return {header.str(),row.str()};
}

}
struct HydrodynamicsOutput::Data {
    explicit Data(const HydrodynamicsInput& i):input(i){}
    const HydrodynamicsInput& input;
    bool owns=false;std::int64_t step=0;double time=0;std::uint64_t geometry=0;
    std::optional<LAMP::TimeLayer> acceptedClock,restoredClock;
    std::uint64_t frameOffset=0;
    Fluid::Totals totals;Fluid::Observation observation;
    std::string status="RUNNING",startTime,hashes;
    struct Frame{std::int64_t step;double time;std::uint64_t geometry;std::string volume,surface;};
    std::vector<Frame> frames;
    // Reuse each scalar file handle; accepted boundaries still publish all
    // rows immediately, and explicit flush/close errors remain App failures.
    std::map<fs::path,std::ofstream> series;
    std::map<fs::path,std::vector<TecplotWriter::TableColumn>> tableSchemas;
    std::map<fs::path,std::unique_ptr<TecplotWriter::TableHistory>> tableSeries;
    bool vofSeries() const{return input.interfaceMethod==HydrodynamicsInput::InterfaceMethod::Vof;}
    void beginSeries(const fs::path& path,const std::string& header){
        if(!vofSeries()){write(path,header);return;}
        auto& schema=tableSchemas[path];if(!schema.empty())throw std::logic_error("duplicate Hydro series schema");
        std::istringstream columns(header);std::string name;
        while(std::getline(columns,name,',')){
            while(!name.empty()&&(name.back()=='\n'||name.back()=='\r'))name.pop_back();
            const bool category=name=="stage"||name=="status"||name=="reason"||name=="event";
            schema.push_back({name,category?TecplotWriter::TableColumnKind::Category:TecplotWriter::TableColumnKind::Numeric});
        }
    }
    void appendTable(const fs::path& logical,const fs::path& actual,const std::string& records){
        const auto schema=tableSchemas.find(logical);
        if(schema==tableSchemas.end())throw std::logic_error("Hydro series row has no declared schema");
        auto found=tableSeries.find(actual);
        if(found==tableSeries.end()){
            std::map<std::string,std::string> metadata{{"App","Hydrodynamics"},{"Interface","LocalPlicVOF"},{"Series",logical.stem().string()}};
            if(logical.filename()=="reconstruction_templates.csv"){
                metadata["StatisticsScope"]="actual owned regular-template bank snapshot at Initial/accepted End";
                metadata["AttemptCounters"]="row_queries,row_hits,row_builds,port_queries,port_hits,port_builds,actual_svd_calls are cumulative bank attempts, including failed preparations; not accepted-only";
                metadata["StorageFields"]="entries and retained_bytes are current owned bank storage; bytes are the accounted retained estimate";
                metadata["Units"]="counts and bytes; not elapsed time and not added to performance series";
            }
            auto history=std::make_unique<TecplotWriter::TableHistory>(actual,logical.stem().string(),schema->second,std::move(metadata));
            found=tableSeries.emplace(actual,std::move(history)).first;
        }
        found->second->append(records);
    }
    void appendVofSeries(const fs::path& path,const std::string& records){
        if(records.empty())return;
        if(path.filename()!="performance.csv"){
            auto actual=path;actual.replace_extension(".dat");appendTable(path,actual,records);return;
        }
        // These provider-owned stage names are fixed identifiers. Keep one
        // append-only stream per stage so every channel remains a real curve.
        std::istringstream lines(records);std::string line;
        while(std::getline(lines,line)){
            if(line.empty())continue;
            const auto first=line.find(','),second=first==std::string::npos?first:line.find(',',first+1);
            const auto third=second==std::string::npos?second:line.find(',',second+1);
            if(third==std::string::npos)throw std::logic_error("malformed Hydro performance row");
            const auto stage=line.substr(second+1,third-second-1);
            if(stage.empty()||std::any_of(stage.begin(),stage.end(),[](char c){
                return !((c>='a'&&c<='z')||(c>='A'&&c<='Z')||(c>='0'&&c<='9')||c=='_');}))
                throw std::logic_error("Hydro performance stage is not a fixed identifier");
            appendTable(path,path.parent_path()/("performance_"+stage+".dat"),line+'\n');
        }
    }
    void append(const fs::path& path,const std::string& text,bool flush=false){
        if(vofSeries()){appendVofSeries(path,text);if(flush)flushSeries();return;}
        auto& out=series[path];
        try{if(!out.is_open()){out.exceptions(std::ios::failbit|std::ios::badbit);out.open(path,std::ios::binary|std::ios::app);}
            out<<text;if(flush)out.flush();}
        catch(const std::exception& e){throw std::runtime_error("Hydro output "+path.generic_string()+": "+e.what());}
    }
    void diagnostic(const fs::path& path,const DiagnosticRecord& record){
        if(tableSchemas.find(path)==tableSchemas.end())beginSeries(path,record.header);
        append(path,record.row);
    }
    void flushSeries(){
        for(auto& file:series)try{file.second.flush();}
            catch(const std::exception& e){throw std::runtime_error("Hydro output flush "+file.first.generic_string()+": "+e.what());}
        for(auto& file:tableSeries)file.second->flush();
    }
    void closeSeries(){
        flushSeries();for(auto& file:series)try{file.second.close();}
            catch(const std::exception& e){throw std::runtime_error("Hydro output close "+file.first.generic_string()+": "+e.what());}
        for(auto& file:tableSeries)file.second->close();
    }
    void buildInfo() const {write("output/meta/build_info.txt",input.interfaceMethod==HydrodynamicsInput::InterfaceMethod::Vof?Fluid::vofBuildInfo(input,observation):Fluid::buildInfo(input,observation));}
    void manifest(){auto s=stream();s<<"{\"app\":"<<quoted("Hydrodynamics")<<",\"method\":\"staggered\",\"interface\":"<<quoted(input.interfaceMethod==HydrodynamicsInput::InterfaceMethod::Vof?"LocalPlicVOF":"FirstJet")<<",\"time_integration\":\"AC\","
        <<"\"status\":"<<quoted(status)<<",\"started_at\":"<<quoted(startTime)<<",\"accepted_step\":"<<step
        <<",\"accepted_time\":"<<time<<",\"geometry_identity\":"<<geometry<<(input.interfaceMethod==HydrodynamicsInput::InterfaceMethod::Vof?Fluid::vofManifestFields():Fluid::manifestFields(input,observation))
        <<(vofSeries()?",\"series_format\":\"tecplotDat\"":"")
        <<",\"series_profile\":"<<quoted(input.run.diagnosticSeries?"diagnostic":"physical")
        <<(vofSeries()?std::string(",\"initialization_method\":\"piecewiseLinearSubcellFan\",\"initialization_once\":true,\"initialization_subcells_per_axis\":")
            +std::to_string(input.vofInitialization.subcellsPerAxis):std::string{})
        <<",\"last_physical_frame_geometry_identity\":"<<(frames.empty()?0:frames.back().geometry)
        <<",\"complete_frames\":"<<frames.size()<<",\"frame_offset\":"<<frameOffset
        <<",\"cumulative_complete_frames\":"<<frameOffset+frames.size()
        <<",\"restart_checkpoint_directory\":"<<quoted(input.run.restartDirectory)
        <<",\"restored_checkpoint_step\":"<<(restoredClock?std::to_string(restoredClock->step):"null")
        <<",\"input_hashes\":"<<hashes<<"}\n";
        replace("output/meta/manifest.json",s.str());}
};
HydrodynamicsOutput::HydrodynamicsOutput(const HydrodynamicsInput& i):data_(std::make_unique<Data>(i)){}
HydrodynamicsOutput::~HydrodynamicsOutput()=default;
std::string HydrodynamicsOutput::checkpoint() const{
    namespace Archive=LAMP::IO::Checkpoint;const auto& d=*data_;
    Archive::require(d.owns&&d.vofSeries()&&d.status=="RUNNING"&&d.acceptedClock.has_value(),
        "Hydro output checkpoint requires a successful actual accepted output");
    const auto t=*d.acceptedClock;(void)t.time();
    Archive::require(t.step==d.step&&t.time()==d.time&&t.position!=LAMP::TimeLayer::Position::Midpoint,
        "Hydro output checkpoint has a different accepted full clock");
    Archive::require(std::isfinite(d.totals.initialVolume)&&d.totals.initialVolume>=0
        &&std::isfinite(d.totals.bodyWork)&&std::isfinite(d.totals.viscousLoss),"nonfinite Hydro output reference totals");
    Archive::require(d.frames.size()<=std::numeric_limits<std::uint64_t>::max()-d.frameOffset,
        "Hydro output cumulative frame counter overflow");
    Archive::Writer out("HydrodynamicsVofOutput.v1");out.scalar(t.step);out.scalar(t.begin);out.scalar(t.end);out.scalar(t.position);
    out.scalar(d.totals.initialVolume);out.scalar(d.totals.bodyWork);out.scalar(d.totals.viscousLoss);
    out.scalar<std::uint64_t>(d.frameOffset+d.frames.size());return out.finish();
}
void HydrodynamicsOutput::restoreCheckpoint(const QuadtreeMesh& mesh,const std::string& bytes,
    const HydroTwoPhaseState& restored){
    namespace Archive=LAMP::IO::Checkpoint;auto& d=*data_;
    Archive::require(d.owns&&d.vofSeries()&&d.status=="RUNNING"&&!d.acceptedClock&&d.frames.empty()
        &&d.frameOffset==0,"Hydro output restore requires a fresh begun VOF writer");
    Archive::Reader in(bytes,"HydrodynamicsVofOutput.v1");LAMP::TimeLayer t;
    t.step=in.scalar<std::int64_t>();t.begin=in.scalar<double>();t.end=in.scalar<double>();t.position=in.scalar<LAMP::TimeLayer::Position>();
    Fluid::Totals totals;totals.initialVolume=in.scalar<long double>();totals.bodyWork=in.scalar<long double>();totals.viscousLoss=in.scalar<long double>();
    const auto frames=in.scalar<std::uint64_t>();in.finish();(void)t.time();
    Archive::require(t.position!=LAMP::TimeLayer::Position::Midpoint&&frames<=static_cast<std::uint64_t>(t.step)+1,
        "Hydro output checkpoint has an invalid accepted clock or frame count");
    Archive::require(std::isfinite(totals.initialVolume)&&totals.initialVolume>=0
        &&std::isfinite(totals.bodyWork)&&std::isfinite(totals.viscousLoss),"invalid Hydro output reference totals");
    Archive::require(restored.model()==HydroPhaseGeometry::Model::GeometricVof&&sameReceiptClock(restored.timeLayer(),t),
        "Hydro output checkpoint differs from the genuine restored full clock");
    const auto& publication=restored.phaseGeometry().vofPublication();
    restored.validateVofFor(mesh,publication.identity(),t);
    const auto observation=Fluid::observeVof(restored,nullptr);
    d.step=t.step;d.time=t.time();d.geometry=publication.identity().value();d.acceptedClock=t;d.restoredClock=t;
    d.totals=totals;d.frameOffset=frames;d.observation=observation;d.buildInfo();d.manifest();
}
void HydrodynamicsOutput::begin(){

    auto& d=*data_;if(fs::exists("output"))throw std::runtime_error("output already exists; use a fresh run directory");
    fs::create_directories("output");d.owns=true;
    for(const char* p:{"fields","series","meta/effective"})fs::create_directories(fs::path("output")/p);
    if(d.input.interfaceMethod!=HydrodynamicsInput::InterfaceMethod::Vof)
        for(const char* p:{"fields/volume","fields/interface","fields/rejected"})fs::create_directories(fs::path("output")/p);
    d.startTime=LAMP::IO::utc();auto hashes=stream(),effective=stream();hashes<<'{';int count=0;
    for(const auto& doc:d.input.documents){
        if(count++)hashes<<',';hashes<<quoted(doc.file)<<':'<<quoted("sha256:"+LAMP::IO::sha256(doc.sourceText));
        const auto name=fs::path(doc.file).filename();write(fs::path("output/meta/effective")/name,doc.effectiveText);
        effective<<"// "<<name.generic_string()<<'\n'<<doc.effectiveText<<'\n';
    }
    hashes<<'}';d.hashes=hashes.str();write("output/meta/effective_config.dict",effective.str());
    auto trace=stream();trace<<"field,file,line,raw_value,effective_value,units,consumer,defaulted,used,unused_reason\n";
    for(const auto& e:d.input.trace.entries)trace<<LAMP::IO::csv(e.field)<<','<<LAMP::IO::csv(e.file)<<','<<e.line<<','
        <<LAMP::IO::csv(e.rawValue)<<','<<LAMP::IO::csv(LAMP::IO::value(e.value))<<','<<LAMP::IO::csv(e.units)<<','
        <<LAMP::IO::csv(e.consumer)<<','<<(e.defaulted?"true":"false")<<','<<(e.used?"true":"false")<<','<<LAMP::IO::csv(e.unusedReason)<<'\n';
    write("output/meta/configuration_trace.csv",trace.str());
    d.buildInfo();
    write("output/meta/frames.csv","frame,step,time,geometry_identity,volume,interface\n");
    std::string conservationHeader="step,time,dt,geometry_identity,negative_phase_volume,initial_volume,relative_volume_drift,total_geometric_mass,divergence_max,dual_R_max,native_max_speed,kappa_max,phase_volume_closure,mass_closure";
    if(d.input.interfaceMethod==HydrodynamicsInput::InterfaceMethod::Vof)conservationHeader+=",liquid_p0_moment_x,liquid_p0_moment_y,liquid_p0_moment_xx,liquid_p0_moment_yy,liquid_p0_moment_xy"
        ",vof_work_kinetic_change,vof_work_capillary,vof_work_body,vof_work_viscous,vof_work_pressure,vof_work_half_R,vof_work_half_R_mapped_primal,vof_work_half_R_geometry_transfer,vof_work_half_R_spatial_exchange,vof_work_time_mass,vof_work_rounded_phase_mass,vof_work_operator_rate_roundoff,vof_work_raw_minus_canonical_outflow_half,vof_work_half_R_decomposition_residual,vof_work_momentum_equation,vof_work_kinetic_budget_residual,vof_work_closure_residual,vof_work_closure_relative_residual,vof_stable_gas_duals,vof_stable_gas_half_R_work,vof_stable_gas_spatial_exchange_half_R_work,vof_stable_gas_max_abs_integrated_R_over_mid_mass,vof_peak_gas_native_face,vof_peak_gas_abs_end_normal_velocity,vof_peak_gas_integrated_R_over_mid_mass,vof_peak_gas_mass_change_over_mid_mass,vof_peak_gas_outward_integral_over_mid_mass,vof_peak_gas_geometry_transfer_over_mid_mass,vof_peak_gas_spatial_exchange_over_mid_mass";
    d.beginSeries("output/series/conservation.csv",conservationHeader+"\n");
    d.beginSeries("output/series/solver.csv","step,time,stage,status,AC_stages,linear_iterations,linear_refinements,algebraic_residual,momentum_residual,reason\n");
    d.beginSeries("output/series/energy.csv","step,time,kinetic,surface,total,body_power,cumulative_body_work,cumulative_viscous_loss,work_corrected_total,capillary_energy_defect,interface_power,interface_power_available,capillary_power,capillary_power_available\n");
    if(d.input.run.diagnosticSeries)
        d.beginSeries("output/series/performance.csv","step,time,stage,wall_seconds,cells,native_velocity_dofs,pressure_dofs\n");
    if(d.input.interfaceMethod==HydrodynamicsInput::InterfaceMethod::Vof){
        d.beginSeries("output/series/vof_csf.csv","step,accepted_time,force_time,force_geometry_identity,force_kappa_min,force_kappa_max,loaded_dofs,height_cells,height_columns,height_column_cell_queries,height_alternate_directions,height_extended_cells,height_extended_slices,height_refined_cells,height_refined_slices,height_fit_cells,height_fit_segments,unresolved_isolated_fragment_cells\n");
        if(d.input.run.diagnosticSeries){
            d.beginSeries("output/series/vof_native_velocity.csv",Fluid::vofNativeVelocityCsvHeader());
            d.beginSeries("output/series/reconstruction_templates.csv",
                "step,time,entries,row_queries,row_hits,row_builds,port_queries,port_hits,port_builds,actual_svd_calls,retained_bytes\n");
        }
        if(d.input.run.adaptive)d.beginSeries("output/series/time_step.csv",
            "step,time,begin,actual_dt,proposed_dt,flow_rate,capillary_rate,interface_h,flow_dt_bound,capillary_dt_bound,growth_dt_bound,max_dt,reason,retries,endpoint_aligned,final_below_min,minimum_clock_rounded\n");
    }
    else if(d.input.step.surfaceTension>0){
        if(d.input.step.capillaryModel==HydroCapillaryModel::BalancedCsf)d.beginSeries("output/series/csf.csv","step,accepted_time,force_time,force_geometry_identity,phase_fraction_min,phase_fraction_max,reference_kappa_min,reference_kappa_max,adjoint_residual_max,loaded_dofs,force_kappa_min,force_kappa_max,force_discretization,height_cells,height_columns,height_column_cell_queries,height_alternate_directions,increment_mean_difference_min,increment_mean_difference_max,increment_reference_load_max,increment_residual_load_max,increment_nodes,increment_rate_points,increment_load_dofs\n");
        else d.beginSeries("output/series/normal_stress.csv","step,accepted_time,force_time,force_geometry_identity,phase_indicator_min,phase_indicator_max,link_kappa_min,link_kappa_max,loaded_links\n");
    }
    d.beginSeries("output/series/events.csv","event,attempted_step,accepted_step,accepted_time,old_dt,new_dt,stage,reason\n");d.manifest();

}
void HydrodynamicsOutput::accepted(const QuadtreeMesh& mesh,const HydroTwoPhaseState& state,const HydroTwoPhaseStep* prepared,
    Clock::time_point start,double compute,bool final){acceptedImpl(mesh,state,prepared,nullptr,start,compute,final);}
void HydrodynamicsOutput::acceptedVof(const QuadtreeMesh& mesh,const HydroTwoPhaseState& state,const HydroVofStep* prepared,
    Clock::time_point start,double compute,bool final,const HydrodynamicsVofRemapReceipt* receipt){
    acceptedImpl(mesh,state,nullptr,prepared,start,compute,final,receipt);}
void HydrodynamicsOutput::acceptedImpl(const QuadtreeMesh& mesh,const HydroTwoPhaseState& state,const HydroTwoPhaseStep* prepared,
    const HydroVofStep* vofStep,Clock::time_point start,double compute,bool final,const HydrodynamicsVofRemapReceipt* receipt){
    auto& d=*data_;const auto& t=state.timeLayer();const bool r=prepared||vofStep;
    const bool vof=state.model()==HydroPhaseGeometry::Model::GeometricVof;
    if((vof&&prepared)||(!vof&&(vofStep||receipt)))throw std::logic_error("VOF output received another model step/receipt");
    if(!d.owns || (r?(t.step!=d.step+1 || t.time()<=d.time):(t.step!=0 || t.time()!=0)))
        throw std::logic_error("Hydro output requires sequential accepted states");
    const auto ioStart=Clock::now();
    
    const auto* transfer=vof?validateReceipt(mesh,state,vofStep,receipt):nullptr;
    const auto observed=vof?Fluid::observeVof(state,vofStep,transfer):Fluid::observe(state,prepared,d.observation);
    d.step=t.step;d.time=t.time();d.geometry=vof?state.phaseGeometry().vofIdentity().value():Fluid::geometryIdentityValue(state);
    if(observed.gradientUpdate!=d.observation.gradientUpdate){d.observation.gradientUpdate=observed.gradientUpdate;d.buildInfo();}
    if(observed.scalarRepresentation!=d.observation.scalarRepresentation){d.observation.scalarRepresentation=observed.scalarRepresentation;d.buildInfo();}
    if(observed.capillaryDiscretization!=d.observation.capillaryDiscretization){d.observation.capillaryDiscretization=observed.capillaryDiscretization;d.buildInfo();}
    d.observation=observed;
    const auto rows=vof?Fluid::vofRows(mesh,state,vofStep,d.input,d.totals,compute,std::cout.flags(),std::cout.precision(),final,transfer):Fluid::rows(mesh,state,prepared,d.input,d.totals,d.observation,compute,std::cout.flags(),std::cout.precision(),final);
    d.totals.initialVolume=rows.next.initialVolume;
    if(!rows.progress.empty()){
        std::cout<<rows.progress;
        std::cout<<std::endl;
    }
    d.append("output/series/conservation.csv",rows.conservation);d.append("output/series/solver.csv",rows.solver);
    d.totals.bodyWork=rows.next.bodyWork;d.totals.viscousLoss=rows.next.viscousLoss;
    auto energy=stream();energy<<rows.energyPrefix;if(vof)energy<<"NA";else energy<<rows.workCorrected;
    energy<<rows.energySuffix;d.append("output/series/energy.csv",energy.str());
    
    if(vof&&d.input.run.diagnosticSeries)
        d.append("output/series/vof_native_velocity.csv",Fluid::vofNativeVelocityRow(mesh,state,vofStep,transfer));
    if(vof&&d.input.run.diagnosticSeries){
        // The accepted target geometry owns this bank; counters also retain
        // unsuccessful preparation attempts, not only accepted-step work.
        const auto templates=state.momentumGeometry().regularTemplateStatistics();
        auto record=stream();record<<t.step<<','<<t.time()<<','<<templates.entries<<','
            <<templates.rowQueries<<','<<templates.rowHits<<','<<templates.rowBuilds<<','
            <<templates.portQueries<<','<<templates.portHits<<','<<templates.portBuilds<<','
            <<templates.actualSvdCalls<<','<<templates.retainedBytes<<'\n';
        d.append("output/series/reconstruction_templates.csv",record.str());
    }
    if(receipt){d.diagnostic("output/series/regrid_decision.csv",decisionRecord(state,*receipt));
        if(transfer)d.diagnostic("output/series/regrid.csv",regridRecord(state,*receipt));}
    if(!rows.capillary.empty())d.append(fs::path("output/series")/rows.capillaryFile,rows.capillary);
    if(d.input.run.writeFields && (t.step==0 || final || t.step%d.input.run.writeEverySteps==0)){
        std::ostringstream number;number<<std::setw(6)<<std::setfill('0')<<t.step;const auto suffix=number.str();
        Data::Frame f{t.step,t.time(),d.geometry,(vof?"volume_":"volume/volume_")+suffix+(vof?".dat":".vtu"),(vof?"interface_":"interface/interface_")+suffix+".dat"};
        // Publish the actual payload only after every same-owner file is complete.
        if(vof)TecplotWriter::write(fs::path("output/fields")/f.volume,Fluid::vofVolumeZone(mesh,state,vofStep,transfer));
        else write(fs::path("output/fields")/f.volume,Fluid::volumeXml(mesh,state));
        const auto surface=vof?Fluid::vofInterfaceZone(state,vofStep,d.input.radial,transfer):Fluid::interfaceZone(state,d.input.radial);
        if(surface.points.empty())f.surface.clear();else TecplotWriter::write(fs::path("output/fields")/f.surface,surface);
        auto frames=d.frames;frames.push_back(f);
        auto pvd=stream(),list=stream();
        if(!vof)pvd<<"<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\"><Collection>\n";
        list<<"frame,step,time,geometry_identity,volume,interface\n";
        for(std::size_t k=0;k<frames.size();++k){const auto& q=frames[k];
            if(!vof)pvd<<"<DataSet timestep=\""<<q.time<<"\" group=\"volume\" part=\"0\" file=\""<<q.volume<<"\"/>\n";
            list<<d.frameOffset+k<<','<<q.step<<','<<q.time<<','<<q.geometry<<",fields/"<<q.volume<<','<<(q.surface.empty()?"":"fields/"+q.surface)<<'\n';
        }
        if(!vof)pvd<<"</Collection></VTKFile>\n";
        replace("output/meta/frames.csv",list.str());
        if(!vof)replace("output/fields/solution.pvd",pvd.str());
        d.frames=std::move(frames);
        if(vof&&d.input.run.detailedLog&&vofStep&&vofStep->hasCapillaryLoad()
            &&vofStep->capillaryLoad().phaseGeometry().vofPublication().geometryKind()==LAMP::Interface::GeometricVofGeometryKind::LocalPlic){
            // A changed-tree receipt retains the physical step's original
            // mesh. The force owner is its Midpoint, never the accepted End.
            const auto& sourceMesh=receipt?*receipt->sourceMesh:mesh;
            const auto forceGeometry=Fluid::actualForceGeometryZone(sourceMesh,vofStep->capillaryLoad());
            fs::create_directories("output/diagnostics");
            TecplotWriter::write(fs::path("output/diagnostics")/("force_geometry_"+suffix+".dat"),forceGeometry);
        }
    }

    d.flushSeries();d.manifest();
    if(d.input.run.diagnosticSeries){
        const auto io=std::chrono::duration<double>(Clock::now()-ioStart).count();
        const auto elapsed=std::chrono::duration<double>(Clock::now()-start).count();
        d.append("output/series/performance.csv",vof?Fluid::vofPerformance(mesh,state,vofStep,compute,io,elapsed,transfer):Fluid::performance(mesh,state,prepared,d.input,d.observation,"HYDRO_AC",compute,io,elapsed),true);
    }
    d.acceptedClock=t;
}
void HydrodynamicsOutput::timeStep(const HydrodynamicsTimeStepReport& report,const LAMP::TimeLayer& t,int retries){
    auto& d=*data_;
    if(!d.owns||!d.vofSeries()||!d.input.run.adaptive||t.position!=LAMP::TimeLayer::Position::End
        ||t.step!=d.step+1||t.begin!=d.time||t.time()<=d.time||retries<0)
        throw std::logic_error("Hydro time step output requires a committed adaptive accepted clock");
    const auto actualDt=t.end-t.begin;
    auto row=stream();row<<t.step<<','<<t.end<<','<<t.begin<<','<<actualDt<<','<<report.proposed_dt<<',';
    const auto value=[&](long double number,bool available){if(available)row<<number;else row<<"NA";row<<',';};
    value(report.flow_rate,report.flow_rate_available);
    value(report.capillary_rate,report.capillary_rate_available);
    value(report.interface_minimum_spacing,report.interface_spacing_available);
    value(report.flow_bound,report.flow_bound_available);
    value(report.capillary_bound,report.capillary_bound_available);
    value(report.growth_bound,report.growth_bound_available);
    value(report.max_bound,report.max_bound_available);
    // reason is the original proposal limiter, even when a later typed retry
    // reduced its actual accepted interval. Only actual endpoint flags follow t.
    const bool endpointAligned=report.endpoint_aligned&&t.end==report.end;
    const bool finalBelowMin=t.end==d.input.run.endTime&&actualDt<d.input.run.adaptive->minDt;
    const bool minimumClockRounded=report.minimum_clock_rounded&&t.end==report.end;
    row<<LAMP::IO::csv(report.limiter)<<','<<retries<<','<<endpointAligned<<','<<finalBelowMin<<','<<minimumClockRounded<<'\n';
    d.append("output/series/time_step.csv",row.str(),true);
}
void HydrodynamicsOutput::retry(std::int64_t attempt,std::int64_t accepted,double time,double oldDt,double newDt,const std::string& reason){
    auto row=stream();row<<"retry,"<<attempt<<','<<accepted<<','<<time<<','<<oldDt<<','<<newDt<<",HYDRO_AC,"<<LAMP::IO::csv(reason)<<'\n';
    data_->append("output/series/events.csv",row.str(),true);
    std::cout<<"[RETRY] step="<<attempt<<" accepted_step="<<accepted<<" accepted_t="<<time
        <<" dt="<<oldDt<<" -> "<<newDt<<" reason="<<reason<<std::endl;
}
void HydrodynamicsOutput::failure(std::int64_t attempt,std::int64_t accepted,double time,const std::string& stage,const std::string& reason){
    auto& d=*data_;if(!d.owns)return;d.step=accepted;d.time=time;d.status="FAILED";
    auto f=stream();f<<"{\"status\":\"FAILED\",\"attempted_step\":"<<attempt<<",\"accepted_step\":"<<accepted<<",\"accepted_time\":"<<time
        <<",\"stage\":"<<quoted(stage)<<",\"reason\":"<<quoted(reason)<<",\"accepted_state_preserved\":true,\"numerical_state_rolled_back\":false}\n";
    write("output/meta/failure.json",f.str());
    // The physical terminal record must survive a later diagnostic I/O error.
    d.manifest();auto row=stream();row<<"failure,"<<attempt<<','<<accepted<<','<<time<<",NA,NA,"
        <<LAMP::IO::csv(stage)<<','<<LAMP::IO::csv(reason)<<'\n';d.append("output/series/events.csv",row.str(),true);d.flushSeries();d.manifest();
    try{if(!d.vofSeries())TecplotWriter::exportCsvHistories("output/series");}
    catch(const std::exception& e){std::cerr<<"[OUTPUT] Tecplot series export failed after preserved physical failure: "<<e.what()<<std::endl;}
}
void HydrodynamicsOutput::complete(){auto& d=*data_;d.closeSeries();if(!d.vofSeries())TecplotWriter::exportCsvHistories("output/series");d.status="COMPLETED";d.manifest();}


