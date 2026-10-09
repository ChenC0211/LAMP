#include "Modules/Hydrodynamics/HydroMidpointStep.h"
#include "Modules/Hydrodynamics/Hydrodynamics.h"
#include "Modules/Hydrodynamics/HydroConservativePressure.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>

namespace {
using namespace Quantum::Hydrodynamics;
void require(bool ok,const char* reason){if(!ok)throw std::runtime_error(reason);}
struct Run {std::vector<double> velocity,volume;double initialEnergy,finalEnergy,balance;};
Run evolve(bool rz,bool viscous,int steps,std::ostream* history){
    const int n=12;const double endTime=.4,dt=endTime/steps,pi=std::acos(-1.);
    DictionaryMesh input;input.minCoord={0,0};input.maxCoord={1,1};input.Nx=input.Ny=n;
    input.maxCellNum=2048;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,rz);
    const auto pressure=HydroConservativePressure::buildDivergenceAdjoint(mesh,domain);
    const auto kind=HydroStepGeometry::complete(rz);
    VEdgeField<double> x(mesh,0);HEdgeField<double> y(mesh,0);
    const auto stream=[&](Vector2D p){return rz?.2*2*pi*p.x*p.x*std::pow(1-p.x*p.x,3)*std::sin(pi*p.y)
        :.2/pi*(std::sin(pi*p.x)*std::sin(pi*p.y)+.17*std::sin(2*pi*p.x)*std::sin(3*pi*p.y));};
    for(const auto face:domain.faces(mesh)){
        const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
        const auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
        (face.horizontal?y.value:x.value)[edge.index()]=face.outwardSign?0:(face.horizontal?1:-1)*(stream(b)-stream(a))/face.area;
    }
    Hydrodynamics hydro(mesh);hydro.initializeState(mesh,kind,{nullptr,&x,&y});
    CellField<double> density(mesh,1),mu(mesh,0);VertField<double> vertexDensity(mesh,1),vertexMu(mesh,0);
    for(int c=0;c<mu.size();++c){const auto p=mesh.cells(mesh.cellIndexToLabel[c]).centerCoord();mu.value[c]=viscous?.02*(1+p.x*p.x+p.y):0;}
    for(int v=0;v<vertexMu.size();++v){const auto p=mesh.verts(mesh.vertIndexToLabel[v]).coord();vertexMu.value[v]=viscous?.02*(1+p.x*p.x+p.y):0;}
    hydro.updateMaterials(mesh,{kind,0,HydroMaterialEvolution::Frozen,HydroNativeMaterial{&density,&mu,&vertexDensity,&vertexMu}});
    HydroSolveOptions solver;solver.linearResidualTolerance=1e-10;
    HydroMidpointIterationOptions iteration;iteration.maximumIterations=30;
    iteration.absoluteVelocityTolerance=1e-13;iteration.relativeVelocityTolerance=1e-12;
    const auto kinetic=[&]{long double sum=0;const auto state=hydro.state();
        for(const auto face:domain.faces(mesh)){const double u=face.horizontal?state.velocityY[face.nativeRow-mesh.vEdgeIndexNum()]:state.velocityX[face.nativeRow];
            sum+=.5L*pressure.duals().volumes()[face.nativeRow].measure*u*u;}
        return static_cast<double>(sum);};
    const double initial=kinetic();double previous=initial,energy=initial,dissipation=0,maxIncrease=0,maxBalance=0,maxDivergence=0,
        maxPressurePower=0,maxConvectionPower=0,maxMass=0,maxDefect=0;int maxIterations=0,totalIterations=0;
    const auto writeHistory=[&](int step){if(history)*history<<std::setprecision(17)<<(rz?"RZ":"XY")<<','<<viscous<<','<<steps<<','
        <<step*dt<<','<<energy<<','<<dissipation<<','<<(energy+dissipation-initial)/initial<<'\n';};
    writeHistory(0);
    if(rz && !viscous && steps==20){
        auto fail=iteration;fail.maximumIterations=1;HydroForceDensity force;force.evaluationTime=.5*dt;force.timeLevel=HydroForceTimeLevel::StepMidpoint;
        bool rejected=false;try{advanceSinglePhaseMidpointStep(mesh,domain,hydro,{1,0,dt,HydroTimePolicy::ConservativeOneStep},solver,force,fail);}
        catch(const std::runtime_error&){rejected=true;}
        const auto restored=hydro.state();double difference=0;
        for(int i=0;i<x.size();++i)difference=std::max(difference,std::abs(restored.velocityX[i]-x.value[i]));
        for(int i=0;i<y.size();++i)difference=std::max(difference,std::abs(restored.velocityY[i]-y.value[i]));
        require(rejected && difference==0 && restored.progress.acceptedSteps==0 && !restored.progress.stopped,
            "unconverged single-phase midpoint did not restore its original physical state");
    }
    for(int step=1;step<=steps;++step){const double begin=(step-1)*dt;
        HydroForceDensity force;force.evaluationTime=begin+.5*dt;force.timeLevel=HydroForceTimeLevel::StepMidpoint;
        const auto result=advanceSinglePhaseMidpointStep(mesh,domain,hydro,
            {step,begin,dt,HydroTimePolicy::ConservativeOneStep},solver,force,iteration);
        const auto& report=result.hydro;
        require(report.midpointNativeMomentumUsed && report.compatibleNativePressureUsed && report.midpointMassEnergyDefect==0,
            "single-phase step bypassed shared midpoint momentum or changed its constant mass");
        energy=kinetic();dissipation-=dt*report.physicalViscousPower;
        maxIncrease=std::max(maxIncrease,(energy-previous)/initial);
        maxBalance=std::max(maxBalance,std::abs((energy+dissipation-initial)/initial));
        maxDivergence=std::max(maxDivergence,report.divergenceAfterProjection);
        maxPressurePower=std::max(maxPressurePower,std::abs(report.physicalPressurePower));
        maxConvectionPower=std::max(maxConvectionPower,std::abs(report.physicalConvectionPower));
        maxMass=std::max(maxMass,report.conservativeMassRelativeClosure);
        maxDefect=std::max(maxDefect,result.midpointVelocityDefect);
        maxIterations=std::max(maxIterations,result.iterations);totalIterations+=result.iterations;
        previous=energy;writeHistory(step);
    }
    const auto state=hydro.state();Run run;run.initialEnergy=initial;run.finalEnergy=energy;run.balance=maxBalance;
    for(const auto face:domain.faces(mesh)){run.velocity.push_back(face.horizontal?state.velocityY[face.nativeRow-mesh.vEdgeIndexNum()]:state.velocityX[face.nativeRow]);
        run.volume.push_back(pressure.duals().volumes()[face.nativeRow].measure);}
    require(state.progress.acceptedSteps==steps && std::abs(state.progress.acceptedEndTime-endTime)<1e-14 && maxBalance<2e-10
        && maxIncrease<2e-10 && maxDivergence<1e-10 && maxPressurePower<1e-12 && maxConvectionPower<1e-12,
        "nonlinear midpoint flow lost its accepted clock, incompressibility or energy identity");
    require(viscous?(energy<.98*initial && dissipation>0):(std::abs(energy-initial)<2e-10*initial && dissipation==0),
        "nonlinear midpoint flow has the wrong inviscid/viscous energy behavior");
    std::cout<<std::setprecision(17)<<"HYDRO_MIDPOINT_TIME metric="<<(rz?"RZ":"XY")<<" variable_viscosity="<<viscous<<" grid="<<n
        <<" steps="<<steps<<" dt="<<dt<<" end_time="<<endTime<<" initial_energy="<<initial<<" final_energy="<<energy
        <<" dissipation="<<dissipation<<" relative_energy_balance="<<maxBalance<<" max_relative_energy_increase="<<maxIncrease
        <<" divergence="<<maxDivergence<<" pressure_power="<<maxPressurePower<<" convection_power="<<maxConvectionPower
        <<" mass_closure="<<maxMass<<" midpoint_defect="<<maxDefect<<" max_iterations="<<maxIterations
        <<" total_iterations="<<totalIterations<<'\n'<<std::flush;
    return run;
}
double difference(const Run& a,const Run& b){
    require(a.velocity.size()==b.velocity.size() && a.volume==b.volume,"temporal comparison changed its spatial system");
    long double error=0,norm=0;for(std::size_t i=0;i<a.velocity.size();++i){const double d=a.velocity[i]-b.velocity[i];
        error+=a.volume[i]*d*d;norm+=a.volume[i]*b.velocity[i]*b.velocity[i];}
    return std::sqrt(static_cast<double>(error/norm));
}
}
void hydroMidpointTimeContract(){
    std::ofstream output;if(const char* path=std::getenv("QUANTUM_MIDPOINT_TIME_HISTORY")){
        output.open(path);require(bool(output),"cannot open midpoint time history");
        output<<"metric,variable_viscosity,steps,time,kinetic_energy,dissipation,relative_energy_balance\n";}
    for(bool rz:{false,true})for(bool viscous:{false,true}){
        const auto coarse=evolve(rz,viscous,20,output.is_open()?&output:nullptr),
            middle=evolve(rz,viscous,40,output.is_open()?&output:nullptr),
            fine=evolve(rz,viscous,80,output.is_open()?&output:nullptr);
        const double first=difference(coarse,middle),second=difference(middle,fine),ratio=first/second;
        std::cout<<std::setprecision(17)<<"HYDRO_MIDPOINT_TIME_ORDER metric="<<(rz?"RZ":"XY")<<" variable_viscosity="<<viscous
            <<" coarse_middle_relative_l2="<<first<<" middle_fine_relative_l2="<<second<<" refinement_ratio="<<ratio
            <<" observed_order="<<std::log2(ratio)<<" fixed_spatial_system=1\n"<<std::flush;
        require(second>1e-10 && ratio>3.6 && ratio<4.5,"nonlinear midpoint flow did not exhibit second-order temporal self-convergence");
    }
}
