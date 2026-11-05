#include "Interface/IO/InterfaceFirstJetOutput.h"
#include "Interface/IO/InterfaceTransportInput.h"
#include "Modules/Interface/HermiteGalsMeshTransfer.h"
#include "Modules/Interface/HermiteGalsFirstJetStep.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/IO/OutputProvenance.h"
#include "Core/IO/TecplotWriter.h"
#include <array>
#include <cmath>
#include <cstdint>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <iostream>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <thread>
#include <vector>

namespace
{
using namespace LAMP::Interface;
namespace fs=std::filesystem;
std::ostringstream stream(){std::ostringstream s;s.imbue(std::locale::classic());s<<std::scientific<<std::setprecision(16);return s;}
void write(const fs::path& path,const std::string& text,bool append=false)
{
    try{std::ofstream file;file.exceptions(std::ios::failbit|std::ios::badbit);
        file.open(path,std::ios::binary|(append?std::ios::app:std::ios::trunc));file<<text;file.close();}
    catch(const std::ios_base::failure& e){throw std::runtime_error("FirstJet output write failed: "+path.generic_string()+": "+e.what());}
}
void replace(const fs::path& path,const std::string& text)
{
    // Preserve the previous complete manifest if serialization fails. Each
    // index is replaced separately; both only refer to complete accepted data.
    const fs::path temporary=path.string()+".tmp";write(temporary,text);
    for(unsigned attempt=0;;++attempt){std::error_code error;fs::rename(temporary,path,error);if(!error)return;
#ifdef _WIN32
        // File indexers/readers can briefly deny replacement on Windows.
        // Retry the same completed payload, never remove the old manifest.
        if(error==std::errc::permission_denied && attempt<6){
            std::cerr<<"[Interface FirstJet IO] atomic replace retry="<<attempt+1<<" file="<<path.generic_string()<<'\n';
            std::this_thread::sleep_for(std::chrono::milliseconds(10u<<attempt));continue;
        }
#endif
        throw std::runtime_error("FirstJet atomic replace failed: "+path.generic_string()+"; error_code="+
            std::to_string(error.value())+"; category="+error.category().name());
    }
}
std::string read(const fs::path& path)
{std::ifstream f(path,std::ios::binary);if(!f)throw std::runtime_error("cannot read FirstJet input provenance: "+path.string());
 std::string s{std::istreambuf_iterator<char>(f),{}};if(f.bad())throw std::runtime_error("FirstJet input provenance read failed");return s;}
using LAMP::IO::quoted;
using LAMP::IO::csv;
using LAMP::IO::value;
using LAMP::IO::sha256;
using LAMP::IO::utc;
using LAMP::IO::compiler;
using LAMP::IO::buildType;
void array(std::ostream& s,const char* name,const char* type,const std::string& values,int components=1)
{s<<"<DataArray type=\""<<type<<"\" Name=\""<<name<<"\" NumberOfComponents=\""<<components<<"\" format=\"ascii\">\n"<<values<<"\n</DataArray>\n";}
void fieldData(std::ostream& s,const HermiteGalsPublication& pub,int meshRevision)
{
    s<<"<FieldData>\n";
    s<<"<DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"ascii\">"<<pub.state().timeLayer().time()<<"</DataArray>\n";
    s<<"<DataArray type=\"UInt64\" Name=\"geometry_identity\" NumberOfTuples=\"1\" format=\"ascii\">"<<pub.state().identity().value()<<"</DataArray>\n";
    s<<"<DataArray type=\"Int32\" Name=\"accepted_step\" NumberOfTuples=\"1\" format=\"ascii\">"<<pub.state().timeLayer().step<<"</DataArray>\n";
    s<<"<DataArray type=\"Int32\" Name=\"mesh_revision\" NumberOfTuples=\"1\" format=\"ascii\">"<<meshRevision<<"</DataArray>\n</FieldData>\n";
}
std::string meshDescription(const QuadtreeMesh& mesh)
{
    auto s=stream();for(int v=0;v<mesh.vertIndexNum();++v){const auto p=mesh.verts(mesh.vertIndexToLabel[v]).coord();s<<p.x<<','<<p.y<<'\n';}
    for(int c=0;c<mesh.cellIndexNum();++c){const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);s<<cell.level();
        for(auto v:cell.cellToVert)s<<','<<mesh.verts(v).index();s<<'\n';}return s.str();
}
std::string volumeXml(const QuadtreeMesh& mesh,const HermiteGalsPublication& pub,int meshRevision)
{
    auto s=stream(),points=stream(),connectivity=stream(),offsets=stream(),types=stream(),level=stream(),ids=stream(),alpha=stream(),volume=stream(),capacity=stream(),phi=stream();
    s<<"<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n<UnstructuredGrid>\n";
    fieldData(s,pub,meshRevision);s<<"<Piece NumberOfPoints=\""<<mesh.vertIndexNum()<<"\" NumberOfCells=\""<<mesh.cellIndexNum()<<"\">\n<Points>\n";
    for(int v=0;v<mesh.vertIndexNum();++v){const auto p=mesh.verts(mesh.vertIndexToLabel[v]).coord();points<<p.x<<' '<<p.y<<" 0\n";}
    array(s,"Points","Float64",points.str(),3);s<<"</Points>\n<Cells>\n";
    std::vector<const HermiteGalsPhaseVolumeCell*> byCell(static_cast<std::size_t>(mesh.cellIndexNum()),nullptr);
    for(const auto& v:pub.volumes().cells())byCell.at(v.cell)=&v;
    for(int c=0;c<mesh.cellIndexNum();++c){const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto* v=byCell[c];
        if(!v || !(v->capacity>0))throw std::runtime_error("FirstJet output lacks same-publication cell measure");
        for(auto corner:cell.cellToVert)connectivity<<mesh.verts(corner).index()<<' ';connectivity<<'\n';offsets<<4*(c+1)<<' ';types<<"9 ";
        level<<cell.level()<<' ';ids<<c<<' ';alpha<<v->rawLiquid/v->capacity<<' ';volume<<v->rawLiquid<<' ';capacity<<v->capacity<<' ';
        phi<<pub.state().geometry().snapshot().scalarJetAtNativeCellChart(c,{{.5,false},{.5,false}}).value<<' ';
    }
    array(s,"connectivity","Int32",connectivity.str());array(s,"offsets","Int32",offsets.str());array(s,"types","UInt8",types.str());
    s<<"</Cells>\n<CellData Scalars=\"phase_volume_fraction\">\n";
    array(s,"level","Int32",level.str());array(s,"cell_id","Int32",ids.str());array(s,"phase_volume_fraction","Float64",alpha.str());
    array(s,"geometric_liquid_volume","Float64",volume.str());array(s,"cell_capacity","Float64",capacity.str());array(s,"phi_center","Float64",phi.str());
    s<<"</CellData>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";return s.str();
}
TecplotWriter::Zone interfaceZone(const HermiteGalsPublication& pub,bool radial,int meshRevision)
{
    TecplotWriter::Zone zone;zone.title="accepted interface";zone.element=TecplotWriter::Element::Line;
    zone.coordinateNames=radial?std::array<std::string,2>{"R","Z"}:std::array<std::string,2>{"X","Y"};
    zone.fields={{radial?"NR":"NX",TecplotWriter::Location::Nodal,{}},
                 {radial?"NZ":"NY",TecplotWriter::Location::Nodal,{}},
                 {"Kappa",TecplotWriter::Location::Nodal,{}}};
    zone.solutionTime=pub.state().timeLayer().time();zone.strandId=1;
    zone.auxiliary={{"AcceptedStep",std::to_string(pub.state().timeLayer().step)},
                    {"GeometryIdentity",std::to_string(pub.state().identity().value())},
                    {"MeshRevision",std::to_string(meshRevision)}};
    std::size_t graphs=0;
    for(std::size_t c=0;c<pub.topology().curves.size();++c)for(std::size_t g=0;g<pub.topology().curves[c].graphs().size();++g){
        const int base=static_cast<int>(zone.points.size());
        for(int i=0;i<9;++i){const auto q=pub.sample(c,g,double(i)/8);const auto& f=q.frame;
            if(q.geometryIdentity!=pub.state().identity() || !std::isfinite(q.curvature))
                throw std::runtime_error("FirstJet output sample identity/curvature invalid");
            zone.points.push_back(f.point);zone.fields[0].values.push_back(f.normal.x);
            zone.fields[1].values.push_back(f.normal.y);zone.fields[2].values.push_back(q.curvature);
            if(i)zone.connectivity.push_back({base+i-1,base+i});
        }++graphs;
    }
    zone.auxiliary["GraphCount"]=std::to_string(graphs);return zone;
}

}

struct InterfaceFirstJetOutput::Data
{
    InterfaceTransportInput input;bool ownsRoot=false,initial=false;std::int64_t step=0;int meshRevision=0;double time=0;
    long double initialVolume=0;std::uint64_t initialGeometry=0,geometry=0;std::size_t leafCells=0;
    std::string runId,startTime,casePath,inputHashes="{}",meshHash,status="running";
    struct Frame {std::int64_t step;double time;std::uint64_t geometry;int meshRevision;std::string volume,surface;};
    std::vector<Frame> frames;
    void manifest()
    {
        auto s=stream();s<<"{\n\"run_id\":"<<quoted(runId)<<",\n\"case_path\":"<<quoted(casePath)<<",\n\"start_time\":"<<quoted(startTime)
            <<",\n\"git_commit\":\"unavailable\",\n\"program_version\":\"LAMP G04\",\n\"build_type\":"<<quoted(buildType())<<",\n\"compiler\":"<<quoted(compiler())
            <<",\n\"selected_methods\":{\"Interface\":\"hermiteGals\",\"transport\":\"FirstJet\"},\n\"scope\":\"standalone prescribed velocity; no Hydro/Electric\","
            <<"\n\"volume_measure\":"<<quoted(input.config.radial?"2*pi*r*dr*dz":"dx*dy")
            <<",\n\"inventory\":null,\n\"phase_flux_ledger\":null,\n\"mesh_identity\":"<<(meshHash.empty()?"null":quoted("sha256:"+meshHash))
            <<",\n\"mesh_revision\":"<<meshRevision<<",\n\"leaf_cells\":"<<leafCells<<",\n\"initial_geometry_identity\":"<<(initial?std::to_string(initialGeometry):"null")
            <<",\n\"geometry_identity\":"<<(initial?std::to_string(geometry):"null")<<",\n\"input_file_hashes\":"<<inputHashes
            <<",\n\"accepted_steps\":"<<step<<",\n\"final_time\":"<<time<<",\n\"accepted_frames\":"<<frames.size()<<",\n\"exit_status\":"<<quoted(status)<<"\n}\n";
        replace("output/meta/manifest.json",s.str());
    }
};
InterfaceFirstJetOutput::InterfaceFirstJetOutput(const InterfaceTransportInput& input):data_(std::make_unique<Data>()){data_->input=input;}
InterfaceFirstJetOutput::~InterfaceFirstJetOutput()=default;
void InterfaceFirstJetOutput::begin()
{
    auto& d=*data_;if(d.ownsRoot)throw std::logic_error("FirstJet output already started");
    if(fs::exists("output") || fs::is_symlink(fs::symlink_status("output")))throw std::runtime_error("FirstJet output already exists; use a fresh run directory");
    if(!fs::create_directory("output"))throw std::runtime_error("cannot create fresh FirstJet output root");d.ownsRoot=true;
    for(const auto* path:{"fields/volume","fields/interface","fields/rejected","series","meta"})fs::create_directories(fs::path("output")/path);
    d.startTime=utc();d.casePath=fs::current_path().generic_u8string();d.runId=sha256(d.casePath+d.startTime+std::to_string(Clock::now().time_since_epoch().count()));
    auto hashes=stream();hashes<<'{';int n=0;for(const auto* path:interfaceTransportInputFiles){
        if(n++)hashes<<',';hashes<<quoted(path)<<":"<<quoted("sha256:"+sha256(read(path)));}hashes<<'}';d.inputHashes=hashes.str();
    const auto& c=d.input.config;
    const auto vec=[](std::ostringstream& out,Vector2D p){out<<'('<<p.x<<' '<<p.y<<')';};
    auto run=stream(),geometry=stream(),mesh=stream(),fields=stream(),interface=stream();
    run<<"Run { app Interface; }\nTime { units code; dt "<<c.dt<<"; endTime "<<c.endTime<<"; }\nOutput { Fields { enabled "<<(c.writeFields?"true":"false")<<';';
    if(c.writeFields)run<<" writeEverySteps "<<c.writeEverySteps<<';';run<<" } }\n";
    geometry<<"Geometry {\ncoordinates "<<(c.radial?"rz":"xy")<<";\nDomain { min ";vec(geometry,c.domainMin);geometry<<"; max ";vec(geometry,c.domainMax);
    geometry<<"; }\nInitial { phi "<<quoted(c.initialPhi)<<"; }\n}\n";
    mesh<<"Mesh {\nresolution "<<c.resolution<<";\ninitialMesh "<<(c.initialAmr?"amr":"uniform")<<";\n";
    if(!c.regridSteps.empty()){mesh<<"Regrid { atSteps (";for(int k:c.regridSteps)mesh<<k<<' ';mesh<<"); }\n";}
    if(c.initialAmr || !c.regridSteps.empty()){mesh<<"Refinement { min ";vec(mesh,c.refineMin);mesh<<"; max ";vec(mesh,c.refineMax);mesh<<"; level "<<c.refinementLevel<<"; }\n";}mesh<<"}\n";
    fields<<"Fields {\n";if(!c.constants.empty()){fields<<"Constants {\n";for(const auto& item:c.constants)fields<<item.first<<' '<<item.second<<";\n";fields<<"}\n";}
    fields<<"Velocity { type "<<(c.velocity==InterfaceTransportConfig::Velocity::Constant?"constant":"polynomial")<<";\n";
    if(c.velocity==InterfaceTransportConfig::Velocity::Constant){fields<<"value ";vec(fields,c.constantVelocity);fields<<";\n";}
    else{fields<<"origin ";vec(fields,c.velocityOrigin);fields<<"; scale ";vec(fields,c.velocityScale);fields<<";\nCoefficients {\n";
        for(int k=0;k<16;++k)fields<<'c'<<k/4<<k%4<<" { x "<<quoted(c.coefficients[k].x)<<"; y "<<quoted(c.coefficients[k].y)<<"; }\n";fields<<"}\n";}fields<<"}\n}\n";
    interface<<"Interface {\nmethod hermiteGals;\nBoundary { type "<<(c.prescribedBoundary?"prescribed":"none")<<';';
    if(c.prescribedBoundary)interface<<" phi "<<quoted(c.boundaryPhi)<<';';interface<<" }\n}\n";
    auto combined=stream();combined<<"// Effective configuration; sections retain their owning input file.\n";
    fs::create_directories("output/meta/effective");
    const std::array<std::string,5> bodies{run.str(),geometry.str(),mesh.str(),fields.str(),interface.str()};
    for(std::size_t i=0;i<bodies.size();++i){const auto name=fs::path(interfaceTransportInputFiles[i]).filename();
        write(fs::path("output/meta/effective")/name,bodies[i]);combined<<"// system/"<<name.generic_string()<<'\n'<<bodies[i]<<'\n';}
    write("output/meta/effective_config.dict",combined.str());
    auto trace=stream();trace<<"field,file,line,raw_value,effective_value,units,consumer,defaulted,used,unused_reason\n";
    for(const auto& e:d.input.trace.entries)trace<<csv(e.field)<<','<<csv(e.file)<<','<<e.line<<','<<csv(e.rawValue)<<','<<csv(value(e.value))<<','<<csv(e.units)<<','
        <<csv(e.consumer)<<','<<(e.defaulted?"true":"false")<<','<<(e.used?"true":"false")<<','<<csv(e.unusedReason)<<'\n';
    write("output/meta/configuration_trace.csv",trace.str());
    write("output/meta/build_info.txt","compiler="+compiler()+"\nbuild_type="+buildType()+"\ngit_commit=unavailable\ncompiled="+__DATE__+" "+__TIME__+"\nVTK=ASCII XML; graph polylines use 9 same-publication samples; no geometric refit\nunits=code\nmesh_identity=SHA256 of ordered vertex coordinates and leaf levels/connectivity\nstep_total excludes performance CSV and manifest serialization\n");
    write("output/meta/frames.csv","frame,step,time,geometry_identity,mesh_revision,volume,interface\n");
    write("output/series/conservation.csv","step,time,dt,geometry_identity,phase_volume,initial_phase_volume,relative_volume_drift,analytic_volume,absolute_volume_error,phase_volume_boundary_flux,phase_volume_closure,mass,mass_boundary_flux,mass_closure,divergence_norm,momentum_balance_norm,inventory,flux_ledger,cut_cells,interface_panels,kappa_max,kappa_rms,curvature_quadrature_samples,curvature_unavailable_samples,curvature_measure\n");
    write("output/series/solver.csv","step,time,module,stage,status,nonlinear_iterations,linear_iterations,linear_residual,physical_residual,rejected_reason\n");
    write("output/series/events.csv","event,step,time,geometry_identity,kind,refined_cells,coarsened_cells,unchanged_cells,relative_domain_volume_change,reason\n");
    write("output/series/performance.csv","step,time,stage,wall_seconds,leaf_cells,native_faces,cut_cells,interface_panels,primary_vertices,primary_edges,centers,hydro_dofs,electric_skeleton_dofs,nnz\n");
    d.manifest();
}
void InterfaceFirstJetOutput::accepted(const QuadtreeMesh& mesh,const HermiteGalsPublication& pub,
    Clock::time_point start,double stepSeconds,double transferSeconds,const HermiteGalsMeshTransferReport* transfer,const HermiteGalsFirstJetStepReport* stepReport)
{
    auto& d=*data_;if(!d.ownsRoot)throw std::logic_error("FirstJet output not started");pub.validateMesh(mesh);
    const auto ioStart=Clock::now();const auto& t=pub.state().timeLayer();
    if(pub.state().primary().transport!=HermiteGalsTransport::FirstJet || !pub.state().primary().phaseVolumes.empty() ||
       (d.initial?(t.step!=d.step+1 || !(t.time()>d.time)):(t.step!=0 || t.time()!=0)))throw std::logic_error("FirstJet output requires sequential App-accepted states");
    const auto volume=pub.volumes().report().totalRawLiquid;
    if(!d.initial){d.initialVolume=volume;d.initialGeometry=pub.state().identity().value();d.initial=true;}
    d.step=t.step;d.time=t.time();d.geometry=pub.state().identity().value();d.leafCells=mesh.cellIndexNum();
    const auto hash=sha256(meshDescription(mesh));if(hash!=d.meshHash){d.meshHash=hash;++d.meshRevision;}
    const auto drift=(volume-d.initialVolume)/d.initialVolume;
    // Same-publication physical quadrature, independent of field output frequency.
    // Scale by the already validated maximum before squaring the curvature.
    const double curvatureScale=pub.report().maxAbsCurvature;
    long double curvatureMeasure=0,curvatureSquare=0;
    std::size_t curvatureSamples=0,curvatureUnavailable=0;
    for(std::size_t c=0;c<pub.topology().curves.size();++c)
        for(std::size_t g=0;g<pub.topology().curves[c].graphs().size();++g)
            for(const auto& q:pub.quadrature(c,g)){
                ++curvatureSamples;
                if(q.sample.geometryIdentity!=pub.state().identity() || !std::isfinite(q.sample.curvature)
                    || !std::isfinite(q.physicalWeight) || !(q.physicalWeight>0)){
                    ++curvatureUnavailable;continue;
                }
                const long double scaled=curvatureScale>0?q.sample.curvature/curvatureScale:0;
                curvatureMeasure+=q.physicalWeight;curvatureSquare+=q.physicalWeight*scaled*scaled;
            }
    if(!std::isfinite(curvatureMeasure) || !std::isfinite(curvatureSquare))
        throw std::runtime_error("FirstJet output curvature diagnostic accumulation nonfinite");
    auto row=stream();row<<d.step<<','<<d.time<<','<<(d.step?d.input.config.dt:0)<<','<<d.geometry<<','<<volume<<','<<d.initialVolume<<','<<drift<<",NA,NA"<<",NA,NA,NA,NA,NA,NA,NA,NA,NA,"<<pub.volumes().report().cutCells<<','<<pub.report().graphs<<','<<pub.report().maxAbsCurvature<<',';
    if(curvatureUnavailable==0 && curvatureMeasure>0)row<<curvatureScale*std::sqrt(curvatureSquare/curvatureMeasure);else row<<"NA";
    row<<','<<curvatureSamples<<','<<curvatureUnavailable<<','<<curvatureMeasure<<'\n';
    write("output/series/conservation.csv",row.str(),true);
    row=stream();row<<d.step<<','<<d.time<<",Interface,"<<(d.step?"first_jet_step":"initialization")<<",accepted,NA,NA,NA,"<<pub.report().maxRelativeInterfaceResidual<<",\n";
    write("output/series/solver.csv",row.str(),true);
    if(transfer){row=stream();row<<"regrid,"<<d.step<<','<<d.time<<','<<d.geometry<<','<<(transfer->kind==HermiteGalsMeshTransferKind::ExactRefinement?"exact":"project")<<','
        <<transfer->refinedCells<<','<<transfer->coarsenedCells<<','<<transfer->unchangedCells<<','<<transfer->relativeDomainVolumeChange<<",\n";write("output/series/events.csv",row.str(),true);}
    if(d.input.config.writeFields && (d.step==0 || d.step==d.input.config.steps || d.step%d.input.config.writeEverySteps==0)){
        std::ostringstream number;number<<std::setw(6)<<std::setfill('0')<<d.step;const auto index=number.str();
        Data::Frame frame{d.step,d.time,d.geometry,d.meshRevision,"volume/volume_"+index+".vtu","interface/interface_"+index+".dat"};
        write(fs::path("output/fields")/frame.volume,volumeXml(mesh,pub,d.meshRevision));
        const auto surface=interfaceZone(pub,d.input.config.radial,d.meshRevision);
        if(surface.points.empty())frame.surface.clear();else TecplotWriter::write(fs::path("output/fields")/frame.surface,surface);
        // Payload files must be complete before either accepted-frame list.
        auto frames=d.frames;frames.push_back(frame);auto pvd=stream(),list=stream();
        pvd<<"<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n<Collection>\n";
        list<<"frame,step,time,geometry_identity,mesh_revision,volume,interface\n";
        for(std::size_t i=0;i<frames.size();++i){const auto& f=frames[i];
            pvd<<"<DataSet timestep=\""<<f.time<<"\" group=\"volume\" part=\"0\" file=\""<<f.volume<<"\"/>\n";
            list<<i<<','<<f.step<<','<<f.time<<','<<f.geometry<<','<<f.meshRevision<<",fields/"<<f.volume<<','<<(f.surface.empty()?"":"fields/"+f.surface)<<'\n';}
        pvd<<"</Collection>\n</VTKFile>\n";
        replace("output/meta/frames.csv",list.str());replace("output/fields/solution.pvd",pvd.str());d.frames=std::move(frames);
    }
    const auto ioSeconds=std::chrono::duration<double>(Clock::now()-ioStart).count(),totalSeconds=std::chrono::duration<double>(Clock::now()-start).count();
    auto performance=stream();const auto stage=[&](const char* name,double seconds,bool producer=false){performance<<d.step<<','<<d.time<<','<<name<<','<<seconds<<',';
        if(producer && stepReport)performance<<stepReport->sourceCells<<','<<stepReport->sourceNativeFaces<<','<<stepReport->sourceCutCells<<','<<stepReport->sourceGraphs<<','
            <<stepReport->transportedVertices+stepReport->prescribedVertices<<','<<stepReport->transportedEdges+stepReport->prescribedEdges<<','<<stepReport->tracedCenters;
        else performance<<mesh.cellIndexNum()<<','<<pub.domain().nativeFaceCount()<<','<<pub.volumes().report().cutCells<<','<<pub.report().graphs<<','<<pub.state().primary().vertices.size()<<','
            <<pub.state().primary().edges.size()<<','<<pub.state().primary().geometryCells.size();performance<<",NA,NA,NA\n";};
    stage(d.step?"interface_step":"interface_initialization",stepSeconds,true);stage("interface_mesh_transfer",transferSeconds);stage("io",ioSeconds);stage("step_total",totalSeconds);
    if(stepReport){const auto& t=stepReport->timings;
        stage("interface_input",t.input,true);stage("interface_boundary",t.boundary,true);stage("interface_vertices",t.vertices,true);stage("interface_edges",t.edges,true);
        stage("interface_centers",t.centers,true);stage("interface_mixed_closure",t.mixed,true);stage("interface_publication",t.publication,true);stage("interface_producer_total",t.total,true);}
    write("output/series/performance.csv",performance.str(),true);d.manifest();
}
void InterfaceFirstJetOutput::complete(long double peakDrift)
{
    auto& d=*data_;auto s=stream();s<<"{\"step\":"<<d.step<<",\"time\":"<<d.time
        <<",\"reference\":null,\"position\":null,\"gradient\":null,\"curvature\":null,\"absolute_volume_error\":null,\"peak_relative_volume_drift\":"<<peakDrift<<"}\n";
    replace("output/meta/accuracy.json",s.str());d.status="success";d.manifest();
    TecplotWriter::exportCsvHistories("output/series");
}
void InterfaceFirstJetOutput::rejected(int attempted,std::int64_t acceptedStep,double acceptedTime,const std::string& stage,const std::string& reason)
{
    auto& d=*data_;if(!d.ownsRoot)return;d.step=acceptedStep;d.time=acceptedTime;d.status="failed";
    auto s=stream();s<<"{\"module\":\"Interface\",\"stage\":"<<quoted(stage)<<",\"attempted_step\":"<<attempted
        <<",\"accepted_step\":"<<acceptedStep<<",\"accepted_time\":"<<acceptedTime<<",\"reason\":"<<quoted(reason)
        <<",\"diagnostics\":{\"accepted_frames\":"<<d.frames.size()<<",\"inventory\":null,\"flux_ledger\":null}}\n";
    replace("output/meta/failure.json",s.str());
    s=stream();s<<attempted<<','<<acceptedTime<<",Interface,"<<csv(stage)<<",rejected,NA,NA,NA,NA,"<<csv(reason)<<'\n';write("output/series/solver.csv",s.str(),true);
    s=stream();s<<"rejected,"<<attempted<<','<<acceptedTime<<','<<d.geometry<<",NA,NA,NA,NA,NA,"<<csv(reason)<<'\n';write("output/series/events.csv",s.str(),true);d.manifest();
    try{TecplotWriter::exportCsvHistories("output/series");}
    catch(const std::exception& e){std::cerr<<"[Interface FirstJet] series Tecplot export failed: "<<e.what()<<'\n';}
    catch(...){std::cerr<<"[Interface FirstJet] series Tecplot export failed with an unknown error\n";}
}
