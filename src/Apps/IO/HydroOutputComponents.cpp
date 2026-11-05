#include "IO/HydroOutputComponents.h"
#include "Core/IO/OutputProvenance.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <sstream>
#include <stdexcept>
namespace AppIO::HydroOutput {
using namespace LAMP::Hydrodynamics;using namespace LAMP::Interface;
using LAMP::IO::stream;using LAMP::IO::quoted;
const char* gradientUpdateName(const LAMP::Interface::HermiteGalsState& state){
    if(state.timeLayer().position==LAMP::TimeLayer::Position::Initial)return "initialData";
    switch(state.primary().gradientUpdate){
    case LAMP::Interface::HermiteGalsGradientUpdate::CharacteristicPullback:return "characteristicPullback";
    case LAMP::Interface::HermiteGalsGradientUpdate::CompatibleIncrement:return "compatibleIncrement";
    default:throw std::logic_error("Hydro output encountered an unknown actual gradient update");
    }
}
const char* scalarRepresentationName(const LevelSetHermiteGeometry::Snapshot& scalar){
    switch(scalar.reconstruction()){
    case LevelSetHermiteGeometry::Reconstruction::EnrichedQuartic:return "singleQuartic";
    case LevelSetHermiteGeometry::Reconstruction::C2MacroQuartic:return "c2MacroQuartic";
    default:throw std::logic_error("Hydro output encountered an unsupported actual scalar representation");
    }
}
void array(std::ostream& out,const char* name,const std::string& values,const char* type="Float64",int components=1){
    out<<"<DataArray type=\""<<type<<"\" Name=\""<<name<<"\" NumberOfComponents=\""<<components
       <<"\" format=\"ascii\">\n"<<values<<"\n</DataArray>\n";
}
std::uint64_t geometryIdentityValue(const HydroTwoPhaseState& state)
{return state.geometryIdentity().value();}
void clock(std::ostream& out,const HydroTwoPhaseState& state){
    const auto t=state.timeLayer();out<<"<FieldData>\n";
    out<<"<DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"ascii\">"<<t.time()<<"</DataArray>\n"
       <<"<DataArray type=\"UInt64\" Name=\"geometry_identity\" NumberOfTuples=\"1\" format=\"ascii\">"<<geometryIdentityValue(state)<<"</DataArray>\n"
       <<"<DataArray type=\"Int64\" Name=\"accepted_step\" NumberOfTuples=\"1\" format=\"ascii\">"<<t.step<<"</DataArray>\n"
       <<"<DataArray type=\"Float64\" Name=\"pressure_sample_time\" NumberOfTuples=\"1\" format=\"ascii\">"
       <<(t.step?t.begin+.5*(t.end-t.begin):0)<<"</DataArray>\n</FieldData>\n";
}
std::string volumeXml(const QuadtreeMesh& mesh,const HydroTwoPhaseState& state){
    const auto& layout=state.values().layout();const auto& phase=state.phaseGeometry();const auto scalar=state.publication().state().geometry();
    const auto velocity=state.reconstruction().evaluate(mesh,state.values(),state.geometryIdentity(),state.timeLayer(),state.boundary());
    auto xml=stream(),points=stream(),conn=stream(),offset=stream(),types=stream(),ids=stream(),cap=stream(),alpha=stream(),phi=stream(),grad=stream(),u=stream(),p=stream();
    for(int v=0;v<mesh.vertIndexNum();++v){const auto q=mesh.verts(mesh.vertIndexToLabel[v]).coord();points<<q.x<<' '<<q.y<<" 0\n";}
    for(int row=0;row<layout.pressureCount();++row){
        const int native=layout.nativeCell(row);const auto& cell=mesh.cells(mesh.cellIndexToLabel[native]);const auto x=cell.centerCoord();
        const auto f=scalar.snapshot().scalarJetInNativeCell(native,x);const auto v=velocity.cellVelocityJet(row,x).jet.value;
        for(auto corner:cell.cellToVert)conn<<mesh.verts(corner).index()<<' ';offset<<4*(row+1)<<' ';types<<"9 ";ids<<native<<' ';
        cap<<phase.metrics().cellVolumes()[row]<<' ';alpha<<phase.liquidFractions()[row]<<' ';phi<<f.value<<' ';grad<<f.dx<<' '<<f.dy<<" 0\n";
        u<<v.x<<' '<<v.y<<" 0\n";p<<state.values().pressure()[row]<<' ';
    }
    xml<<"<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\"><UnstructuredGrid>\n";clock(xml,state);
    xml<<"<Piece NumberOfPoints=\""<<mesh.vertIndexNum()<<"\" NumberOfCells=\""<<layout.pressureCount()<<"\"><Points>\n";
    array(xml,"Points",points.str(),"Float64",3);xml<<"</Points><Cells>\n";array(xml,"connectivity",conn.str(),"Int32");array(xml,"offsets",offset.str(),"Int32");
    array(xml,"types",types.str(),"UInt8");xml<<"</Cells><CellData>\n";array(xml,"cell_id",ids.str(),"Int32");
    array(xml,"cell_capacity",cap.str());array(xml,"phase_volume_fraction",alpha.str());array(xml,"phi",phi.str());
    array(xml,"grad_phi",grad.str(),"Float64",3);array(xml,"velocity",u.str(),"Float64",3);array(xml,"pressure",p.str());
    xml<<"</CellData></Piece></UnstructuredGrid></VTKFile>\n";return xml.str();
}
TecplotWriter::Zone interfaceZone(const HydroTwoPhaseState& state,bool radial){
    const auto& pub=state.publication();TecplotWriter::Zone zone;
    zone.title="accepted interface";zone.element=TecplotWriter::Element::Line;
    zone.coordinateNames=radial?std::array<std::string,2>{"R","Z"}:std::array<std::string,2>{"X","Y"};
    zone.fields={{radial?"NR":"NX",TecplotWriter::Location::Nodal,{}},
                 {radial?"NZ":"NY",TecplotWriter::Location::Nodal,{}},
                 {"Kappa",TecplotWriter::Location::Nodal,{}}};
    zone.solutionTime=state.timeLayer().time();zone.strandId=1;
    zone.auxiliary={{"AcceptedStep",std::to_string(state.timeLayer().step)},
                    {"GeometryIdentity",std::to_string(state.geometryIdentity().value())},
                    {"GradientUpdate",gradientUpdateName(pub.state())},
                    {"ScalarRepresentation",scalarRepresentationName(pub.state().geometry().snapshot())}};
    std::size_t graphs=0;
    for(std::size_t c=0;c<pub.topology().curves.size();++c)for(std::size_t g=0;g<pub.topology().curves[c].graphs().size();++g){
        const int base=static_cast<int>(zone.points.size());
        for(int k=0;k<9;++k){const auto q=pub.sample(c,g,double(k)/8);const auto& f=q.frame;
            if(q.geometryIdentity!=state.geometryIdentity())throw std::runtime_error("Hydro output interface identity mismatch");
            zone.points.push_back(f.point);zone.fields[0].values.push_back(f.normal.x);
            zone.fields[1].values.push_back(f.normal.y);zone.fields[2].values.push_back(q.curvature);
            if(k)zone.connectivity.push_back({base+k-1,base+k});
        }++graphs;
    }
    zone.auxiliary["GraphCount"]=std::to_string(graphs);return zone;
}
long double measure(const HermiteGalsPublication& pub){
    long double result=0;for(std::size_t c=0;c<pub.topology().curves.size();++c)
        for(std::size_t g=0;g<pub.topology().curves[c].graphs().size();++g)for(const auto& q:pub.quadrature(c,g))result+=q.physicalWeight;
    return result;
}
Observation observe(const HydroTwoPhaseState& state,const HydroTwoPhaseStep* prepared,const Observation& previous){
    auto out=previous;
    out.gradientUpdate=gradientUpdateName(state.publication().state());
    out.scalarRepresentation=scalarRepresentationName(state.publication().state().geometry().snapshot());
    out.heightFunctionForce=false;out.incrementForce=false;out.faceCurvatureDefined=false;out.faceCurvatureMin=0;out.faceCurvatureMax=0;
    if(prepared && prepared->hasCapillaryLoad() && prepared->capillaryModel()==HydroCapillaryModel::BalancedCsf){
        const auto kind=prepared->capillaryLoad().discretization();
        out.incrementForce=kind==HydroCsfDiscretization::GeometricHeightFunctionIncrement;
        out.heightFunctionForce=kind==HydroCsfDiscretization::GeometricHeightFunction||out.incrementForce;
        if(out.heightFunctionForce){
            const auto& load=prepared->capillaryLoad();const auto& values=load.heightFunctionFaceCurvatures();
            const auto& support=load.supportReports();
            for(std::size_t u=0;u<values.size();++u)if(support[u].curvatureDefined){
                if(!out.faceCurvatureDefined){out.faceCurvatureMin=out.faceCurvatureMax=values[u];out.faceCurvatureDefined=true;}
                else{out.faceCurvatureMin=std::min(out.faceCurvatureMin,values[u]);out.faceCurvatureMax=std::max(out.faceCurvatureMax,values[u]);}
            }
        }
        const std::string observed=kind==HydroCsfDiscretization::GeometricHeightFunctionIncrement?"geometricHeightFunctionIncrement":
            kind==HydroCsfDiscretization::GeometricHeightFunction?"geometricHeightFunction":
            kind==HydroCsfDiscretization::GeometricProjectedAdjoint?"geometricProjectedAdjoint":
            kind==HydroCsfDiscretization::GeometricAdjoint?"geometricAdjoint":
            kind==HydroCsfDiscretization::GeometricPhaseFraction?"geometricPhaseFraction":"regularizedColor";
        out.capillaryDiscretization=observed;
    }
    return out;
}
Rows rows(const QuadtreeMesh& mesh,const HydroTwoPhaseState& state,const HydroTwoPhaseStep* prepared,
    const HydrodynamicsInput& input,const Totals& totals,const Observation& observation,double compute,
    std::ios_base::fmtflags flags,std::streamsize precision,bool final){
    Rows out;out.next=totals;const auto& t=state.timeLayer();const auto* r=prepared?&prepared->report():nullptr;
    std::ostringstream progress;progress.flags(flags);progress.precision(precision);
    long double liquid=0,mass=0;for(const auto& c:state.phaseGeometry().cellVolumes()){liquid+=c.liquid;mass+=input.step.liquidDensity*c.liquid+input.step.gasDensity*c.gas;}
    if(!r)out.next.initialVolume=liquid;
    if(!(out.next.initialVolume>0))throw std::runtime_error("initial negative phase volume must be positive");
    const auto drift=(liquid-out.next.initialVolume)/out.next.initialVolume;
    double umax=0,divMax=0;for(double u:state.values().normalVelocity())umax=std::max(umax,std::abs(u));
    for(double v:state.phaseGeometry().metrics().apply(mesh,state.values()))divMax=std::max(divMax,std::abs(v));
    long double R=0;if(prepared)for(const auto& c:prepared->dualMaterial().report().components)R=std::max(R,c.maxAbsDefectRate);
    // Print and persist the real accepted prefix before any field serialization.
    if(t.step==0 || t.step==1 || final || t.step%input.run.logEverySteps==0){
        progress<<"[Step "<<t.step<<"] t="<<t.time()<<" dt="<<(r?r->deltaTime:0)<<" cells="<<mesh.cellIndexNum()
            <<" scalar="<<observation.scalarRepresentation<<" gradient="<<observation.gradientUpdate<<" AC="<<(r?r->iterations:0)<<" LIN=NA div="<<divMax<<" dV/V="<<drift<<" umax="<<umax
            <<" kappa="<<state.publication().report().maxAbsCurvature<<" wall="<<compute<<"s";
        if(observation.heightFunctionForce){
            progress<<" force="<<(observation.incrementForce?"heightFunctionIncrement":"heightFunction")
                <<(observation.incrementForce?" kappa_HF_base_max=":" kappa_force_max=");
            if(observation.faceCurvatureDefined)progress<<std::max(std::abs(observation.faceCurvatureMin),std::abs(observation.faceCurvatureMax));
            else progress<<"NA";
        }
        out.progress=progress.str();
    }
    auto row=stream();row<<t.step<<','<<t.time()<<','<<(r?r->deltaTime:0)<<','<<state.geometryIdentity().value()<<','<<liquid<<','<<out.next.initialVolume<<','<<drift
        <<','<<mass<<','<<divMax<<',';if(r)row<<R;else row<<"NA";
    row<<','<<umax<<','<<state.publication().report().maxAbsCurvature<<",NA,NA\n";out.conservation=row.str();
    row=stream();row<<t.step<<','<<t.time()<<','<<(r?"AC":"initialization")<<",accepted,"<<(r?r->iterations:0)<<",NA,";
    if(prepared){const auto& s=prepared->solveResult().report();row<<s.linearRefinements<<','<<s.algebraicRelativeResidual<<','<<r->frozenOwnerMomentumDefect;}else row<<"NA,NA,NA";
    row<<",\n";out.solver=row.str();
    const long double surface=input.step.surfaceTension==0?0:input.step.surfaceTension*(r?r->interfaceMeasureEnd:measure(state.publication()));
    const long double kinetic=r?r->kineticEnd:0;
    if(r){out.next.bodyWork+=r->deltaTime*r->bodyPower;out.next.viscousLoss+=r->deltaTime*r->viscousPower;}
    row=stream();row<<t.step<<','<<t.time()<<','<<kinetic<<','<<surface<<','<<kinetic+surface<<','<<(r?r->bodyPower:0)
        <<','<<out.next.bodyWork<<','<<out.next.viscousLoss<<',';
    out.energyPrefix=row.str();out.workCorrected=kinetic+surface-out.next.bodyWork+out.next.viscousLoss;row=stream();
    row<<',';if(r&&r->capillaryPowerAvailable)row<<r->capillaryEnergyDefect;else row<<"NA";
    row<<',';if(r)row<<r->interfacePower;else row<<"NA";
    row<<','<<(r?"true":"false")<<',';
    if(r&&r->capillaryPowerAvailable)row<<r->capillaryPower;else row<<"NA";
    row<<','<<(r&&r->capillaryPowerAvailable?"true":"false")<<'\n';
    out.energySuffix=row.str();
    if(prepared && prepared->hasCapillaryLoad()&&prepared->capillaryModel()==HydroCapillaryModel::BalancedCsf){
        const auto& csf=prepared->capillaryLoad();const auto& colors=csf.phaseGeometry().liquidFractions();
        const auto colorRange=std::minmax_element(colors.begin(),colors.end());
        double kmin=0,kmax=0;bool defined=false;
        const bool adjoint=csf.discretization()==HydroCsfDiscretization::GeometricAdjoint
            ||csf.discretization()==HydroCsfDiscretization::GeometricProjectedAdjoint;
        if(adjoint)for(double k:csf.referenceCurvatures())if(std::isfinite(k)){
            if(!defined){kmin=kmax=k;defined=true;}
            else{kmin=std::min(kmin,k);kmax=std::max(kmax,k);}
        }
        double residualMax=0;std::size_t loaded=0;
        if(adjoint)for(double f:csf.adjointResidualLoad())residualMax=std::max(residualMax,std::abs(f));
        for(double f:csf.momentumRhsLoad())if(f!=0)++loaded;
        row=stream();row<<t.step<<','<<t.time()<<','<<csf.timeLayer().time()<<','<<csf.geometryIdentity().value()
            <<','<<*colorRange.first<<','<<*colorRange.second<<',';
        if(defined)row<<kmin<<','<<kmax;else row<<"NA,NA";
        row<<',';if(adjoint)row<<residualMax;else row<<"NA";
        row<<','<<loaded<<',';
        if(observation.faceCurvatureDefined)row<<observation.faceCurvatureMin<<','<<observation.faceCurvatureMax;else row<<"NA,NA";
        row<<','<<observation.capillaryDiscretization<<',';
        if(observation.heightFunctionForce){const auto& h=csf.curvatureReport();
            row<<h.heightFunctionCells<<','<<h.heightFunctionColumns<<','<<h.heightFunctionCellQueries<<','<<h.heightFunctionAlternateDirections;
        }else row<<"NA,NA,NA,NA";
        row<<',';
        if(observation.incrementForce){
            double lo=0,hi=0,ref=0,res=0;bool found=false;
            for(double k:csf.heightFunctionIncrementMeanCurvatureDifference())if(std::isfinite(k)){
                if(!found){lo=hi=k;found=true;}else{lo=std::min(lo,k);hi=std::max(hi,k);}
            }
            for(double f:csf.heightFunctionIncrementReferenceLoad())ref=std::max(ref,std::abs(f));
            for(double f:csf.heightFunctionIncrementResidualLoad())res=std::max(res,std::abs(f));
            if(found)row<<lo<<','<<hi;else row<<"NA,NA";
            const auto& h=csf.curvatureReport();
            row<<','<<ref<<','<<res<<','<<h.incrementNodes<<','<<h.incrementRatePoints<<','<<h.incrementLoadDofs;
        }else row<<"NA,NA,NA,NA,NA,NA,NA";
        row<<'\n';out.capillaryFile="csf.csv";out.capillary=row.str();
    }else if(prepared&&prepared->hasCapillaryLoad()){
        const auto& sharp=prepared->sharpNormalStressLoad();
        const auto range=std::minmax_element(sharp.phaseIndicators().begin(),sharp.phaseIndicators().end());
        double kmin=0,kmax=0;int links=0;
        for(std::size_t u=0;u<sharp.linkLengths().size();++u)if(sharp.linkLengths()[u]>0){
            const double k=sharp.linkCurvatures()[u];if(!links)kmin=kmax=k;else{kmin=std::min(kmin,k);kmax=std::max(kmax,k);}++links;
        }
        row=stream();row<<t.step<<','<<t.time()<<','<<sharp.geometry().timeLayer().time()<<','<<sharp.geometry().geometryIdentity().value()
            <<','<<*range.first<<','<<*range.second<<',';if(links)row<<kmin<<','<<kmax;else row<<"NA,NA";
        row<<','<<links<<'\n';out.capillaryFile="normal_stress.csv";out.capillary=row.str();
    }
    return out;
}
std::string performance(const QuadtreeMesh& mesh,const HydroTwoPhaseState& state,const HydroTwoPhaseStep* prepared,
    const HydrodynamicsInput& input,const Observation& observation,const char* acceptedStage,double compute,double io,double elapsed){
    const auto& t=state.timeLayer();const auto* r=prepared?&prepared->report():nullptr;auto performance=stream();
    const auto stage=[&](const char* name,double seconds){performance<<t.step<<','<<t.time()<<','<<name<<','<<seconds<<','<<mesh.cellIndexNum()
        <<','<<state.values().layout().velocityCount()<<','<<state.values().layout().pressureCount()<<'\n';};
    stage(r?acceptedStage:"initialization",compute);stage("io",io);stage("step_total",elapsed);
    if(r && input.run.detailedLog){constexpr std::array<const char*,8> names{{"endpoint","transport","mass","convection","capillary","viscosity","saddle","reassessment"}};
        for(std::size_t k=0;k<names.size();++k){double seconds=0;for(const auto& h:r->history)seconds+=h.stageSeconds[k];stage(names[k],seconds);}}
    if(r && prepared->hasCapillaryLoad()){
        double forceSeconds=0;
        if(prepared->capillaryModel()==HydroCapillaryModel::BalancedCsf){
            for(const auto& h:r->history)forceSeconds+=h.csfProjectionSeconds;
            stage(prepared->capillaryLoad().discretization()==HydroCsfDiscretization::GeometricProjectedAdjoint
                ?"csf_projected_motion_load":observation.incrementForce?"csf_height_function_increment_load":observation.heightFunctionForce?"csf_height_function_load":"csf_adjoint_load",forceSeconds);
        }else{
            for(const auto& h:r->history)forceSeconds+=h.stageSeconds[4];
            stage("sharp_normal_stress",forceSeconds);
        }
    }
    return performance.str();
}
std::string buildInfo(const HydrodynamicsInput& input,const Observation& observation){
    const bool csfModel=input.step.capillaryModel==HydroCapillaryModel::BalancedCsf;
    const char* densityRepresentation=input.step.momentumDensityRepresentation==HydroMomentumMass::DensityRepresentation::SharpGeometry
        ?"sharpGeometry":"cellAverageProjection";
    return "compiler="+LAMP::IO::compiler()+"\nbuild_type="+LAMP::IO::buildType()+
        "\ncompiled="+__DATE__+" "+__TIME__+"\ngit_commit=unavailable\nunits=code\ncapillary_model="+
        (csfModel?"balancedCsf":"sharpNormalStress")+"\ncapillary_discretization="+
        (csfModel?observation.capillaryDiscretization:"nativePressureCenterJump")+"\ndensity_representation="+densityRepresentation+
        "\ngradient_update="+observation.gradientUpdate+
        "\nscalar_representation="+observation.scalarRepresentation+
        "\ncsf_force_setting="+(input.step.csfCurvature.forceMethod==HydroCsfForceMethod::HeightFunctionIncrement?"heightFunctionIncrement":
        input.step.csfCurvature.forceMethod==HydroCsfForceMethod::HeightFunction?"heightFunction":"motionAdjoint")+
        "\npressure=initial supplied zero; thereafter accepted step midpoint saddle cell pressure\nvelocity=VTU center samples use retained reconstruction of actual native face normal values\ninterface=9 raw Hermite samples of same accepted publication; no refit\ncapillary="+
        (csfModel?"same-publication pressure-compatible CSF; heightFunction uses phase-volume heights and face curvature; heightFunctionIncrement adds the compatibleIncrement residual reaction (HF curve columns describe only its base); motionAdjoint uses its retained motion residual; actual discretization is recorded above":"same-publication pressure-center phase and curvature at actual native link zero")+
        "; separate diagnostics use the retained midpoint force clock\n"+
        "";
}
std::string manifestFields(const HydrodynamicsInput& input,const Observation& observation){
    auto s=stream();s
          <<",\"density_representation\":"<<quoted(input.step.momentumDensityRepresentation==HydroMomentumMass::DensityRepresentation::SharpGeometry?"sharpGeometry":"cellAverageProjection")
          <<",\"gradient_update\":"<<quoted(observation.gradientUpdate)
          <<",\"scalar_representation\":"<<quoted(observation.scalarRepresentation)
          <<",\"capillary_method\":"<<quoted(input.step.capillaryModel==HydroCapillaryModel::BalancedCsf?"balancedCsf":"sharpNormalStress")
          <<",\"capillary_discretization\":"<<quoted(input.step.capillaryModel==HydroCapillaryModel::BalancedCsf?observation.capillaryDiscretization:"nativePressureCenterJump")
;return s.str();}
}

