#include "Electrohydrodynamics/IO/ElectrohydrodynamicsOutput.h"
#include "IO/HydroOutputComponents.h"
#include "Core/IO/OutputProvenance.h"
#include "Core/IO/TecplotWriter.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Modules/Electric/ElectricPotential.h"
#include "Modules/Electric/ElectricMaxwellTrace.h"
#include "Modules/Electric/SurfaceChargeField.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <locale>
#include <sstream>
#include <stdexcept>
namespace {
using namespace LAMP::Hydrodynamics;using namespace LAMP::Interface;
using namespace LAMP::Electric;
using LAMP::IO::stream;using LAMP::IO::quoted;
namespace fs=std::filesystem;namespace Fluid=AppIO::HydroOutput;
void write(const fs::path& path,const std::string& text,bool append=false){
    std::ofstream out;out.exceptions(std::ios::failbit|std::ios::badbit);
    try{out.open(path,std::ios::binary|(append?std::ios::app:std::ios::trunc));out<<text;out.close();}
    catch(const std::exception& e){throw std::runtime_error("Hydro output "+path.generic_string()+": "+e.what());}
}
void replace(const fs::path& path,const std::string& text){
    const fs::path temp=path.string()+".tmp";write(temp,text);
    std::error_code error;fs::rename(temp,path,error);
    if(error)throw std::runtime_error("Hydro output index replacement "+path.generic_string()+": "+error.message());
}
void requireElectric(bool condition,const char* reason){
    if(!condition)throw std::runtime_error(std::string("Hydro electric output: ")+reason);
}
bool sameClock(HermiteGalsTimeLayer a,HermiteGalsTimeLayer b){
    return a.step==b.step && a.begin==b.begin && a.end==b.end && a.position==b.position;
}
std::string number(double value){auto out=stream();out<<value;return out.str();}
std::array<double,3> basis(double x){return {{1,x,.5*(3*x*x-1)}};}
double polynomial(const std::array<double,3>& c,double x){
    const auto b=basis(x);return c[0]*b[0]+c[1]*b[1]+c[2]*b[2];
}
double polynomial(const std::array<double,4>& c,double x){
    return polynomial(std::array<double,3>{{c[0],c[1],c[2]}},x)+c[3]*.5*(5*x*x*x-3*x);
}
void sameArithmetic(long double defect,long double scale,double tolerance,const char* reason){
    requireElectric(std::isfinite(defect)&&std::isfinite(scale)
        &&(scale>0?defect<=tolerance*scale:defect==0),reason);
}
// Checks retained arithmetic receipts, never preparing a second field/trace.
void validateElectric(const QuadtreeMesh& mesh,const HydroTwoPhaseState& hydro,
    const ElectrohydrodynamicsOutputView& view,const ElectrohydrodynamicsInput& input){
    const auto& config=input.electric;const auto t=hydro.timeLayer();
    const auto& pub=hydro.publication();const auto& field=view.potential;
    const auto& current=view.current;const auto& maxwell=view.maxwell;
    requireElectric(t.position==(t.step==0?HermiteGalsTimeLayer::Position::Initial:HermiteGalsTimeLayer::Position::End),
        "requires an actual accepted Initial/End clock");
    requireElectric((t.step==0)==(view.report==nullptr),"time report must be present exactly for accepted End");
    requireElectric(!view.report||view.report->dt==t.end-t.begin,"time report interval differs from End clock");
    requireElectric(sameClock(view.accepted.timeLayer(),t)&&sameClock(view.accepted.geometrySourceTimeLayer(),t)
        &&sameClock(field.timeLayer(),t)&&sameClock(field.geometry().timeLayer(),t)
        &&sameClock(current.timeLayer(),t)&&sameClock(current.geometry().timeLayer(),t)
        &&sameClock(maxwell.timeLayer(),t)&&sameClock(maxwell.geometrySourceTimeLayer(),t),
        "accepted Q/V/current/Maxwell or geometry full clock differs from Hydro");
    field.validateSource(mesh,pub);current.geometry().validateSource(mesh,pub);
    current.momentOperator().geometry().validateSource(mesh,pub);
    view.accepted.solver().momentOperator().geometry().validateSource(mesh,pub);
    maxwell.validateSource(mesh,pub,t);
    requireElectric(field.axisymmetric()==input.fluid.radial&&input.fluid.radial
        &&field.epsilonByPhase()==config.electricProblem.epsilonByPhase
        &&current.epsilonByPhase()==config.electricProblem.epsilonByPhase
        &&current.conductivityByPhase()==config.conductivityByPhase,
        "coordinate measure or phase epsilon/conductivity differs from configured input");
    const auto* load=field.surfaceChargeLoad();
    requireElectric(load!=nullptr,"actual potential has no retained integrated-Q load");
    load->validateGeometry(field.geometry());
    requireElectric(sameClock(load->momentOperator().timeLayer(),t)
        &&sameClock(current.momentOperator().timeLayer(),t)
        &&sameClock(view.accepted.solver().momentOperator().timeLayer(),t),
        "panel moments do not belong to the accepted geometry clock");
    const auto& charge=view.accepted.integratedCharge();
    requireElectric(charge==load->panelCharge()&&charge.size()==current.panelCurrent().size()
        &&charge.size()==current.momentOperator().panelCount()
        &&charge.size()==view.accepted.solver().momentOperator().panelCount(),
        "accepted Q and potential/current panel ordering differ");
    requireElectric(view.accepted.componentCount()==pub.topology().curves.size(),
        "accepted component coverage differs from Hydro publication");
    std::size_t flat=0;
    for(std::size_t c=0;c<view.accepted.componentCount();++c){
        for(double q:view.accepted.field(c).transport().panelCharge()){
            requireElectric(flat<charge.size()&&charge[flat]==q,"accepted material Field Q was edited outside its owner");
            const auto& a=view.accepted.solver().momentOperator();const auto& b=current.momentOperator();
            const auto& d=load->momentOperator();
            requireElectric(a.component(flat)==static_cast<int>(c)&&a.component(flat)==b.component(flat)
                &&a.component(flat)==d.component(flat)&&a.graph(flat)==b.graph(flat)&&a.graph(flat)==d.graph(flat)
                &&a.interfaceFace(flat)==b.interfaceFace(flat)&&a.interfaceFace(flat)==d.interfaceFace(flat),
                "accepted Q/current/potential do not use the same native panel ordering");
            ++flat;
        }
    }
    requireElectric(flat==charge.size(),"accepted material Field coverage differs from flat Q");
    const auto applied=current.apply(field.skeletonSolution());
    requireElectric(applied.size()==charge.size(),"actual current affine receipt has different size");
    long double defect=0,scale=0;
    for(std::size_t i=0;i<charge.size();++i){
        defect+=std::abs(static_cast<long double>(applied[i])-current.panelCurrent()[i]);
        scale+=std::abs(current.offset()[i])+std::abs(current.panelCurrent()[i]);
    }
    for(int column=0;column<current.matrix().outerSize();++column)
        for(ElectricPanelCurrentOperator::Matrix::InnerIterator e(current.matrix(),column);e;++e)
            scale+=std::abs(static_cast<long double>(e.value())*field.skeletonSolution()[static_cast<std::size_t>(e.col())]);
    sameArithmetic(defect,scale,config.electric.relativeConservationTolerance,
        "actual potential and current affine receipts differ");
    requireElectric(maxwell.interfaceFaces().size()==field.report().interfaceFaces,
        "actual Maxwell interface coverage differs from the potential");
    long double displacementDefect=0,displacementScale=0;
    for(const auto face:maxwell.interfaceFaces()){
        const auto minus=field.normalDisplacementPolynomialP2(face,true);
        const auto plus=field.normalDisplacementPolynomialP2(face,false);
        const auto shared=field.skeletalFaceTrace(face),moments=field.interfaceChargeMoments(face);
        const auto completed=maxwell.potentialPolynomialP3(face);
        std::array<long double,3> qMoments{{0,0,0}},vDefect{{0,0,0}};
        long double qScale=0,vScale=0;
        for(const auto& node:maxwell.quadrature(face)){
            const auto b=basis(node.coordinate);
            const std::array<double,2> expected{{polynomial(minus,node.coordinate),polynomial(plus,node.coordinate)}};
            for(std::size_t phase=0;phase<2;++phase){
                displacementDefect+=node.physicalWeight*std::abs(static_cast<long double>(expected[phase])-node.normalDisplacement[phase]);
                displacementScale+=node.physicalWeight*(std::abs(static_cast<long double>(expected[phase]))+std::abs(node.normalDisplacement[phase]));
            }
            const auto v=polynomial(completed,node.coordinate),vShared=polynomial(shared,node.coordinate);
            for(std::size_t mode=0;mode<3;++mode){
                qMoments[mode]+=static_cast<long double>(node.physicalWeight)*b[mode]*node.chargeDensity;
                qScale+=std::abs(static_cast<long double>(node.physicalWeight)*b[mode]*node.chargeDensity);
                vDefect[mode]+=static_cast<long double>(node.physicalWeight)*b[mode]*(v-vShared);
                vScale+=node.physicalWeight*std::abs(b[mode])*(std::abs(static_cast<long double>(v))+std::abs(vShared));
            }
        }
        long double qDefect=0,traceDefect=0;
        for(std::size_t mode=0;mode<3;++mode){
            qDefect=std::max(qDefect,std::abs(qMoments[mode]-moments[mode]));
            qScale+=std::abs(moments[mode]);traceDefect=std::max(traceDefect,std::abs(vDefect[mode]));
        }
        sameArithmetic(qDefect,qScale,config.electric.relativeConservationTolerance,"actual Q and Maxwell charge moments differ");
        sameArithmetic(traceDefect,vScale,ElectricMaxwellControls{}.relativeProjectionTolerance,"actual potential and Maxwell shared trace differ");
    }
    sameArithmetic(displacementDefect,displacementScale,config.electric.relativeConservationTolerance,
        "actual potential and Maxwell normal displacement differ");
}
TecplotWriter::Zone electricZone(const HydroTwoPhaseState& state,const char* title,int strand){
    TecplotWriter::Zone zone;zone.title=title;zone.element=TecplotWriter::Element::OrderedPoints;
    zone.coordinateNames={{"R","Z"}};zone.solutionTime=state.timeLayer().time();zone.strandId=strand;
    const auto t=state.timeLayer();
    zone.auxiliary={{"Units","code"},{"Method","HHO"},{"InterfaceMotion","moving"},
        {"AcceptedStep",std::to_string(t.step)},{"GeometryIdentity",std::to_string(state.geometryIdentity().value())},
        {"ClockPosition",t.position==HermiteGalsTimeLayer::Position::Initial?"Initial":"End"},
        {"ClockBegin",number(t.begin)},{"ClockEnd",number(t.end)},
        {"ClockMeaning","same accepted Hydro/Interface/Q/V/current/Maxwell clock"}};
    return zone;
}
TecplotWriter::Zone electricPhaseZone(const HydroTwoPhaseState& state,
    const ElectrohydrodynamicsOutputView& view,const ElectrohydrodynamicsInput& input,bool liquid){
    auto zone=electricZone(state,liquid?"accepted liquid electric samples":"accepted gas electric samples",liquid?3:4);
    zone.auxiliary["PhaseId"]=std::to_string(liquid?input.fluid.negativePhase:input.fluid.positivePhase);
    zone.auxiliary["PhaseSign"]=liquid?"negative phi":"positive phi";
    zone.auxiliary["Epsilon"]=number(view.potential.epsilonByPhase()[liquid?0:1]);
    zone.auxiliary["Conductivity"]=number(input.electric.conductivityByPhase[liquid?0:1]);
    zone.auxiliary["Sampling"]="one retained positive-weight physical volume node nearest phase element center; one-sided P3 sample, not a phase average";
    zone.fields={{"V",TecplotWriter::Location::Nodal,{}},{"Er",TecplotWriter::Location::Nodal,{}},
        {"Ez",TecplotWriter::Location::Nodal,{}}};
    const auto& field=view.potential;
    for(std::size_t e=0;e<field.geometry().elements().size();++e){
        if(field.liquid(e)!=liquid)continue;
        const auto& volume=field.geometry().elements()[e].geometry.volume;
        requireElectric(!volume.empty(),"phase element has no retained volume quadrature");
        const auto center=field.elementCenter(e);
        const auto distance=[&](const auto& q){const auto delta=q.point-center;return delta.x*delta.x+delta.y*delta.y;};
        const auto node=std::min_element(volume.begin(),volume.end(),[&](const auto& a,const auto& b){return distance(a)<distance(b);});
        requireElectric(node->weight>0&&node->point.x>0,"phase sample is not a positive-weight physical volume node");
        const auto gradient=field.gradientInElement(e,node->point);
        zone.points.push_back(node->point);zone.fields[0].values.push_back(field.potentialInElement(e,node->point));
        zone.fields[1].values.push_back(-gradient.x);zone.fields[2].values.push_back(-gradient.y);
    }
    requireElectric(!zone.points.empty(),"phase has no actual electric samples");return zone;
}
TecplotWriter::Zone electricInterfaceZone(const HydroTwoPhaseState& state,
    const ElectrohydrodynamicsOutputView& view){
    auto zone=electricZone(state,"accepted interface electric and Maxwell samples",5);
    zone.auxiliary["Sampling"]="original Maxwell nodes on the shared interface skeleton; scatter view, no graph refit";
    zone.auxiliary["NormalOrientation"]="negative phi to positive phi";
    zone.auxiliary["TractionMeaning"]="(T_gas-T_liquid)n; tangent=(-nz,nr)";
    zone.auxiliary["ChargeMeaning"]="q_h derived from accepted integrated Q by shared weighted P2 load";
    zone.auxiliary["PotentialMeaning"]="shared P2 trace with retained weighted cubic completion";
    for(const char* name:{"NR","NZ","V_Gamma","q_h","Et","Dn_liquid","Dn_gas","Maxwell_normal","Maxwell_tangent","Maxwell_r","Maxwell_z"})
        zone.fields.push_back({name,TecplotWriter::Location::Nodal,{}});
    for(const auto face:view.maxwell.interfaceFaces()){
        const auto potential=view.maxwell.potentialPolynomialP3(face);
        for(const auto& q:view.maxwell.quadrature(face)){
            zone.points.push_back(q.point);
            const std::array<double,11> values{{q.normal.x,q.normal.y,polynomial(potential,q.coordinate),
                q.chargeDensity,q.tangentialField,q.normalDisplacement[0],q.normalDisplacement[1],
                q.normalTraction,q.tangentialTraction,q.tractionJump.x,q.tractionJump.y}};
            for(std::size_t k=0;k<values.size();++k)zone.fields[k].values.push_back(values[k]);
        }
    }
    requireElectric(!zone.points.empty(),"interface has no actual Maxwell samples");return zone;
}
}
struct ElectrohydrodynamicsOutput::Data {
    explicit Data(const ElectrohydrodynamicsInput& i):input(i.fluid),joint(i){}
    const HydrodynamicsInput& input;const ElectrohydrodynamicsInput& joint;
    bool owns=false;std::int64_t step=0;double time=0;std::uint64_t geometry=0;
    Fluid::Totals totals;Fluid::Observation observation;
    std::string status="RUNNING",startTime,hashes;
    struct Frame{std::int64_t step;double time;std::uint64_t geometry;std::string volume,surface,electricLiquid,electricGas,electricSurface;};
    std::vector<Frame> frames;
    void buildInfo() const {write("output/meta/build_info.txt",Fluid::buildInfo(input,observation)+std::string("electric=HHO\ninterface_motion=moving\naccepted_clock=one common Hydro/Interface/Q/V/current/Maxwell Initial or End; no stage frames\ncharge_integration=BE startup; variable BDF2 thereafter\nelectric_phase_field=one retained positive-weight volume node per phase element; one-sided P3 point sample, not an average\nelectric_interface_field=actual Maxwell nodes and retained weighted cubic trace; equilibrated Dn\nenergy=kinetic+surface is a partial mechanical diagnostic, not total electrohydrodynamic energy; complete retained interface RHS power is separate\n"));}
    void manifest(){auto s=stream();s<<"{\"app\":"<<quoted("Electrohydrodynamics")<<",\"method\":\"staggered\",\"interface\":\"FirstJet\",\"time_integration\":\"AC\","
        <<"\"status\":"<<quoted(status)<<",\"started_at\":"<<quoted(startTime)<<",\"accepted_step\":"<<step
        <<",\"accepted_time\":"<<time<<",\"geometry_identity\":"<<geometry<<Fluid::manifestFields(input,observation)
        <<",\"electric\":\"HHO\",\"electric_interface_motion\":\"moving\""
        <<",\"complete_frames\":"<<frames.size()<<",\"input_hashes\":"<<hashes<<"}\n";
        replace("output/meta/manifest.json",s.str());}
};
ElectrohydrodynamicsOutput::ElectrohydrodynamicsOutput(const ElectrohydrodynamicsInput& i):data_(std::make_unique<Data>(i)){}
ElectrohydrodynamicsOutput::~ElectrohydrodynamicsOutput()=default;
void ElectrohydrodynamicsOutput::begin(){

    auto& d=*data_;if(fs::exists("output"))throw std::runtime_error("output already exists; use a fresh run directory");
    fs::create_directories("output");d.owns=true;
    for(const char* p:{"fields/volume","fields/interface","fields/rejected","series","meta/effective"})fs::create_directories(fs::path("output")/p);
    d.startTime=LAMP::IO::utc();auto hashes=stream(),effective=stream();hashes<<'{';int count=0;
    for(const auto& doc:d.joint.documents){
        if(count++)hashes<<',';hashes<<quoted(doc.file)<<':'<<quoted("sha256:"+LAMP::IO::sha256(doc.sourceText));
        const auto name=fs::path(doc.file).filename();write(fs::path("output/meta/effective")/name,doc.effectiveText);
        effective<<"// "<<name.generic_string()<<'\n'<<doc.effectiveText<<'\n';
    }
    hashes<<'}';d.hashes=hashes.str();write("output/meta/effective_config.dict",effective.str());
    auto trace=stream();trace<<"field,file,line,raw_value,effective_value,units,consumer,defaulted,used,unused_reason\n";
    for(const auto& e:d.joint.trace.entries)trace<<LAMP::IO::csv(e.field)<<','<<LAMP::IO::csv(e.file)<<','<<e.line<<','
        <<LAMP::IO::csv(e.rawValue)<<','<<LAMP::IO::csv(LAMP::IO::value(e.value))<<','<<LAMP::IO::csv(e.units)<<','
        <<LAMP::IO::csv(e.consumer)<<','<<(e.defaulted?"true":"false")<<','<<(e.used?"true":"false")<<','<<LAMP::IO::csv(e.unusedReason)<<'\n';
    write("output/meta/configuration_trace.csv",trace.str());
    d.buildInfo();
    write("output/meta/frames.csv","frame,step,time,geometry_identity,volume,interface,electric_liquid,electric_gas,electric_interface\n");
    write("output/series/conservation.csv","step,time,dt,geometry_identity,negative_phase_volume,initial_volume,relative_volume_drift,total_geometric_mass,divergence_max,dual_R_max,native_max_speed,kappa_max,phase_volume_closure,mass_closure\n");
    write("output/series/solver.csv","step,time,stage,status,AC_stages,linear_iterations,linear_refinements,algebraic_residual,momentum_residual,reason\n");
    write("output/series/energy.csv","step,time,kinetic,surface,total,body_power,cumulative_body_work,cumulative_viscous_loss,work_corrected_total,capillary_energy_defect,interface_power,interface_power_available,capillary_power,capillary_power_available\n");
    write("output/series/performance.csv","step,time,stage,wall_seconds,cells,native_velocity_dofs,pressure_dofs\n");
    write("output/series/electric.csv","step,time,geometry_identity,total_Q,panel_current,exterior_current,bulk_current,current_balance,current_balance_relative,linear_residual,local_Gauss_relative,global_Gauss,global_Gauss_relative,time_update_present,temporal_order,dt,previous_dt,a0,a1,a2,Schur_iterations,Schur_relative,panel_equation_relative,charge_equation,charge_equation_relative,Maxwell_charge_jump_relative\n");
    { // Joint output retains capillary diagnostics even when surface tension is zero.
        if(d.input.step.capillaryModel==HydroCapillaryModel::BalancedCsf)write("output/series/csf.csv","step,accepted_time,force_time,force_geometry_identity,phase_fraction_min,phase_fraction_max,reference_kappa_min,reference_kappa_max,adjoint_residual_max,loaded_dofs,force_kappa_min,force_kappa_max,force_discretization,height_cells,height_columns,height_column_cell_queries,height_alternate_directions,increment_mean_difference_min,increment_mean_difference_max,increment_reference_load_max,increment_residual_load_max,increment_nodes,increment_rate_points,increment_load_dofs\n");
        else write("output/series/normal_stress.csv","step,accepted_time,force_time,force_geometry_identity,phase_indicator_min,phase_indicator_max,link_kappa_min,link_kappa_max,loaded_links\n");
    }
    write("output/series/events.csv","event,attempted_step,accepted_step,accepted_time,old_dt,new_dt,stage,reason\n");d.manifest();

}
void ElectrohydrodynamicsOutput::accepted(const QuadtreeMesh& mesh,const HydroTwoPhaseState& state,const HydroTwoPhaseStep* prepared,
    Clock::time_point start,double compute,bool final,const ElectrohydrodynamicsOutputView& electric){
    auto& d=*data_;const auto& t=state.timeLayer();const auto* r=prepared?&prepared->report():nullptr;
    if(!d.owns || (r?(t.step!=d.step+1 || t.time()<=d.time):(t.step!=0 || t.time()!=0)))
        throw std::logic_error("Hydro output requires sequential accepted states");
    const auto ioStart=Clock::now();
    validateElectric(mesh,state,electric,d.joint);
    d.step=t.step;d.time=t.time();d.geometry=state.geometryIdentity().value();
    const auto observed=Fluid::observe(state,prepared,d.observation);
    if(observed.gradientUpdate!=d.observation.gradientUpdate){d.observation.gradientUpdate=observed.gradientUpdate;d.buildInfo();}
    if(observed.scalarRepresentation!=d.observation.scalarRepresentation){d.observation.scalarRepresentation=observed.scalarRepresentation;d.buildInfo();}
    if(observed.capillaryDiscretization!=d.observation.capillaryDiscretization){d.observation.capillaryDiscretization=observed.capillaryDiscretization;d.buildInfo();}
    d.observation=observed;
    const auto rows=Fluid::rows(mesh,state,prepared,d.input,d.totals,d.observation,compute,std::cout.flags(),std::cout.precision(),final);
    d.totals.initialVolume=rows.next.initialVolume;
    if(!rows.progress.empty()){
        std::cout<<rows.progress;
        const auto& e=electric.potential.report();std::cout<<" Esym="<<e.skeleton.relativeSymmetryDefect<<" Gauss="<<e.relativeGlobalGaussResidual<<" Qeq=";
        if(electric.report)std::cout<<electric.report->relativeTotalEquationResidual;else std::cout<<"NA";
        std::cout<<std::endl;
    }
    d.manifest();
    write("output/series/conservation.csv",rows.conservation,true);write("output/series/solver.csv",rows.solver,true);
    d.totals.bodyWork=rows.next.bodyWork;d.totals.viscousLoss=rows.next.viscousLoss;
    auto energy=stream();energy<<rows.energyPrefix;energy<<"NA";
    energy<<rows.energySuffix;write("output/series/energy.csv",energy.str(),true);
        auto row=stream();

        const auto& e=electric.potential.report();const auto& c=electric.current.report();
        long double totalCharge=0;for(double q:electric.accepted.integratedCharge())totalCharge+=q;
        row=stream();row<<t.step<<','<<t.time()<<','<<d.geometry<<','<<totalCharge<<','
            <<c.totalPanelCurrent<<','<<c.exteriorCurrent<<','<<c.bulkCurrent<<','<<c.globalCurrentBalance<<','
            <<c.relativeGlobalCurrentBalance<<','<<e.skeleton.relativeLinearResidual<<','<<e.maxRelativeLocalGaussResidual<<','
            <<e.globalGaussResidual<<','<<e.relativeGlobalGaussResidual<<','<<(electric.report?"true":"false")<<',';
        if(const auto* charge=electric.report)
            row<<charge->temporalOrder<<','<<charge->dt<<','<<charge->previousDt<<','<<charge->a0<<','<<charge->a1<<','<<charge->a2<<','
                <<charge->iterations<<','<<charge->relativeSchurResidual<<','<<charge->maxRelativePanelEquationResidual<<','
                <<charge->integratedEquationResidual<<','<<charge->relativeTotalEquationResidual;
        else row<<"NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA";
        row<<','<<electric.maxwell.report().maxRelativeChargeJumpResidual<<'\n';
        write("output/series/electric.csv",row.str(),true);
    
    if(!rows.capillary.empty())write(fs::path("output/series")/rows.capillaryFile,rows.capillary,true);
    if(d.input.run.writeFields && (t.step==0 || final || t.step%d.input.run.writeEverySteps==0)){
        std::ostringstream number;number<<std::setw(6)<<std::setfill('0')<<t.step;const auto suffix=number.str();
        Data::Frame f{t.step,t.time(),d.geometry,"volume/volume_"+suffix+".vtu","interface/interface_"+suffix+".dat","","",""};
        // Complete all same-owner payloads before adding any PVD entry.
        write(fs::path("output/fields")/f.volume,Fluid::volumeXml(mesh,state));
        const auto surface=Fluid::interfaceZone(state,d.input.radial);
        if(surface.points.empty())f.surface.clear();else TecplotWriter::write(fs::path("output/fields")/f.surface,surface);
        {
            f.electricLiquid="volume/electric_liquid_"+suffix+".dat";f.electricGas="volume/electric_gas_"+suffix+".dat";
            f.electricSurface="interface/electric_interface_"+suffix+".dat";
            TecplotWriter::write(fs::path("output/fields")/f.electricLiquid,electricPhaseZone(state,electric,d.joint,true));
            TecplotWriter::write(fs::path("output/fields")/f.electricGas,electricPhaseZone(state,electric,d.joint,false));
            TecplotWriter::write(fs::path("output/fields")/f.electricSurface,electricInterfaceZone(state,electric));
        }
        auto frames=d.frames;frames.push_back(f);
        auto pvd=stream(),list=stream();pvd<<"<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\"><Collection>\n";
        list<<"frame,step,time,geometry_identity,volume,interface,electric_liquid,electric_gas,electric_interface\n";
        for(std::size_t k=0;k<frames.size();++k){const auto& q=frames[k];
            pvd<<"<DataSet timestep=\""<<q.time<<"\" group=\"volume\" part=\"0\" file=\""<<q.volume<<"\"/>\n";
            list<<k<<','<<q.step<<','<<q.time<<','<<q.geometry<<",fields/"<<q.volume<<','<<(q.surface.empty()?"":"fields/"+q.surface);
            list<<','<<(q.electricLiquid.empty()?"":"fields/"+q.electricLiquid)
                <<','<<(q.electricGas.empty()?"":"fields/"+q.electricGas)<<','<<(q.electricSurface.empty()?"":"fields/"+q.electricSurface);
            list<<'\n';
        }pvd<<"</Collection></VTKFile>\n";replace("output/meta/frames.csv",list.str());replace("output/fields/solution.pvd",pvd.str());d.frames=std::move(frames);
    }

    const auto io=std::chrono::duration<double>(Clock::now()-ioStart).count();
    write("output/series/performance.csv",Fluid::performance(mesh,state,prepared,d.input,d.observation,"ELECTROHYDRO_AC",compute,io,
        std::chrono::duration<double>(Clock::now()-start).count()),true);d.manifest();
}
void ElectrohydrodynamicsOutput::retry(std::int64_t attempt,std::int64_t accepted,double time,double oldDt,double newDt,const std::string& reason){
    auto row=stream();row<<"retry,"<<attempt<<','<<accepted<<','<<time<<','<<oldDt<<','<<newDt<<",HYDRO_AC,"<<LAMP::IO::csv(reason)<<'\n';
    write("output/series/events.csv",row.str(),true);
    std::cout<<"[RETRY] step="<<attempt<<" accepted_step="<<accepted<<" accepted_t="<<time
        <<" dt="<<oldDt<<" -> "<<newDt<<" reason="<<reason<<std::endl;
}
void ElectrohydrodynamicsOutput::failure(std::int64_t attempt,std::int64_t accepted,double time,const std::string& stage,const std::string& reason){
    auto& d=*data_;if(!d.owns)return;d.step=accepted;d.time=time;d.status="FAILED";
    auto f=stream();f<<"{\"status\":\"FAILED\",\"attempted_step\":"<<attempt<<",\"accepted_step\":"<<accepted<<",\"accepted_time\":"<<time
        <<",\"stage\":"<<quoted(stage)<<",\"reason\":"<<quoted(reason)<<",\"accepted_state_preserved\":true,\"numerical_state_rolled_back\":false}\n";
    write("output/meta/failure.json",f.str());auto row=stream();row<<"failure,"<<attempt<<','<<accepted<<','<<time<<",NA,NA,"
        <<LAMP::IO::csv(stage)<<','<<LAMP::IO::csv(reason)<<'\n';write("output/series/events.csv",row.str(),true);d.manifest();
    try{TecplotWriter::exportCsvHistories("output/series");}
    catch(const std::exception& e){std::cerr<<"[Electrohydrodynamics] series Tecplot export failed: "<<e.what()<<'\n';}
    catch(...){std::cerr<<"[Electrohydrodynamics] series Tecplot export failed with an unknown error\n";}
}
void ElectrohydrodynamicsOutput::complete(){
    auto& d=*data_;d.status="COMPLETED";d.manifest();
    TecplotWriter::exportCsvHistories("output/series");
}
