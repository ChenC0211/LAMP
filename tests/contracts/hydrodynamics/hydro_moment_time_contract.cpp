#include "Modules/Hydrodynamics/HydroCompatibleMomentFlow.h"
#include "Modules/Hydrodynamics/HydroCompatibleMomentBoundary.h"
#include "Modules/Hydrodynamics/HydroStep.h"
#include "hydro_complete_field_observation.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Geometry/ParametricHermiteCurve.h"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>

void hydroMomentPortsContract(){
    using Space=Quantum::Hydrodynamics::HydroCompatibleMomentSpace;
    using Stress=Quantum::Hydrodynamics::HydroCompatibleMomentViscosity;
    using Flow=Quantum::Hydrodynamics::HydroCompatibleMomentFlow;
    const auto check=[](bool ok,const char* why){if(!ok)throw std::runtime_error(why);};
    std::ofstream history;if(const char* output=std::getenv("QUANTUM_MOMENT_PORTS_OUTPUT")){
        history.open(output);check(bool(history),"cannot open port control history");
        history<<std::setprecision(17)<<"rz,viscosity,step,time,kinetic,viscous,traction_work,reaction_work,kinetic_flux,inflow_dissipation,energy_balance,advective_identity,volume_rate,equation\n";
    }
    for(bool rz:{false,true})for(double mu:{0.,.02}){
        DictionaryMesh input;input.minCoord={0,0};input.maxCoord={1,1};input.Nx=input.Ny=4;input.maxCellNum=2048;input.maxInitialIterations=0;
        QuadtreeMesh mesh(input);DictionaryMesh::RefineBox fine;fine.min={0,.5};fine.max={.5,1};fine.level=1;mesh.applyRefineBoxes({fine});
        FaceSegmentGraph graph;graph.build(mesh);auto domain=NativeLeafDomain::build(mesh,graph,rz);auto space=Space::build(mesh,domain);
        std::vector<HydroBoundaryPatch> ports{{"supply",HydroBoundarySide::South,HydroBoundaryType::NormalVelocity,0,1,-.1},
            {"outlet",HydroBoundarySide::North,HydroBoundaryType::PressureOutlet,0,1,.2}};
        auto flow=Flow::create(mesh,domain,space,Stress::buildCellwise(mesh,domain,space,std::vector<double>(mesh.cellIndexNum(),mu)),1.,Eigen::VectorXd::Zero(space.mass().rows()),0,ports);
        const auto plug=space.interpolateVelocity(mesh,domain,[](Vector2D){return Eigen::Vector2d{0,.1};});
        const double initial=flow.kineticEnergy();
        check(flow.initialBoundaryCorrection()>0&&!flow.state().stageTime&&flow.state().time==0,"boundary initialization claimed a physical stage");
        const double expectedRate=.1*(rz?std::acos(-1.):1.);double budget=0,maxBudget=0,maxFlux=0,maxEquation=0,maxAdvective=0,maxPressure=0,maxPlug=0,dissipation=0;
        const int steps=mu==0?3:12;Flow::Options options;options.relativeResidualTolerance=1e-10;
        for(int step=1;step<=steps;step++){
            HydroStepTime time{step,flow.state().time,.01,HydroTimePolicy::Fixed};
            const auto load=flow.makeLoad(mesh,domain,Eigen::VectorXd::Zero(space.mass().rows()),time.begin+.5*time.dt);
            auto candidate=flow.prepareNavierStokesStep(mesh,domain,time,load,options);const auto report=candidate.report();
            check(flow.state().acceptedSteps==step-1,"port preparation advanced the source clock");
            flow.accept(mesh,domain,std::move(candidate));
            budget+=report.viscousDissipation+report.centralConvectionWork+report.transportDissipation
                -report.boundaryTractionWork-report.boundaryReactionWork-report.pressureWork;
            const double error=report.kineticAfter-initial+budget,advective=report.centralConvectionWork-report.boundaryKineticFlux-report.boundaryInflowDissipation;
            maxBudget=std::max(maxBudget,std::abs(error)/initial);maxEquation=std::max(maxEquation,report.equationResidual);
            maxAdvective=std::max(maxAdvective,std::abs(advective)/initial);dissipation+=report.viscousDissipation;
            maxFlux=std::max({maxFlux,std::abs(report.boundaryInflowRate-expectedRate)/expectedRate,std::abs(report.boundaryOutflowRate-expectedRate)/expectedRate});
            if(mu==0){const double error=HydroCompleteObservation::errorSquared(mesh,domain,*flow.state().physicalVelocity,[](Vector2D){return Vector2D{0,.1};});
                maxPlug=std::max(maxPlug,std::sqrt(error/(2*initial)));
                for(int c=0;c<mesh.cellIndexNum();c++){maxPressure=std::max(maxPressure,std::abs(flow.state().pressure[3*c]-.2));
                    maxPressure=std::max(maxPressure,flow.state().pressure.segment<2>(3*c+1).norm());}}
            check(report.boundaryInflowDissipation>=0&&report.transportDissipation==0,"central port flow introduced interior dissipation");
            if(history)history<<rz<<','<<mu<<','<<step<<','<<flow.state().time<<','<<report.kineticAfter<<','<<report.viscousDissipation<<','<<report.boundaryTractionWork<<','<<report.boundaryReactionWork<<','
                <<report.boundaryKineticFlux<<','<<report.boundaryInflowDissipation<<','<<error<<','<<advective<<','<<report.boundaryVolumeRate<<','<<report.equationResidual<<'\n';
        }
        std::cout<<std::setprecision(17)<<"MOMENT_PORTS metric="<<(rz?"RZ":"XY")<<" cells="<<mesh.cellIndexNum()<<" viscosity="<<mu<<" steps="<<steps<<" time="<<flow.state().time
            <<" flux_relative="<<maxFlux<<" energy_relative="<<maxBudget<<" advective_identity="<<maxAdvective<<" equation="<<maxEquation<<" plug_error="<<maxPlug<<" pressure_error="<<maxPressure<<" viscous="<<dissipation<<'\n';
        check(maxFlux<2e-11&&maxBudget<2e-9&&maxAdvective<2e-10&&maxEquation<=1e-10,"native open-boundary mass/energy identity failed");
        if(mu==0)check(maxPlug<2e-11&&maxPressure<2e-10&&dissipation==0,"inviscid through-flow or absolute outlet pressure failed");
        else check(dissipation>0,"viscous port flow did not dissipate");
        if(mu==0){
            auto backflowPorts=ports;backflowPorts[0].value=.1;
            auto reversed=Flow::create(mesh,domain,space,Stress::buildCellwise(mesh,domain,space,std::vector<double>(mesh.cellIndexNum(),0.)),1.,Eigen::VectorXd::Zero(space.mass().rows()),0,backflowPorts);
            const auto before=reversed.state().velocity;const auto load=reversed.makeLoad(mesh,domain,Eigen::VectorXd::Zero(space.mass().rows()),.005);bool rejected=false;
            try{reversed.prepareNavierStokesStep(mesh,domain,{1,0,.01,HydroTimePolicy::Fixed},load);}
            catch(const std::invalid_argument& error){rejected=std::string(error.what()).find("backflow")!=std::string::npos;}
            check(rejected&&reversed.state().time==0&&(reversed.state().velocity-before).norm()==0,"undefined pressure backflow did not reject transactionally");
            const auto transport=Quantum::Hydrodynamics::HydroCompatibleMomentTransport::build(mesh,domain,space);
            Quantum::Hydrodynamics::HydroCompatibleMomentBoundary unbound;rejected=false;
            try{transport.assemble(mesh,domain,space,plug,1.,Quantum::Hydrodynamics::HydroCompatibleMomentTransport::Flux::Central,&unbound);}
            catch(const std::invalid_argument&){rejected=true;}
            check(rejected,"unbound native port data accepted");
            std::cout<<"MOMENT_PORTS_REJECTION metric="<<(rz?"RZ":"XY")<<" backflow=1 source_unchanged=1 unbound_domain=1\n";
        }
    }
}

namespace {
using Space=Quantum::Hydrodynamics::HydroCompatibleMomentSpace;
using Stress=Quantum::Hydrodynamics::HydroCompatibleMomentViscosity;
using Flow=Quantum::Hydrodynamics::HydroCompatibleMomentFlow;
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
Eigen::Vector2d profile(Vector2D p){
    const double x=p.x,y=p.y,f=x*x*(1-x)*(1-x),g=y*y*(1-y)*(1-y),fx=2*x-6*x*x+4*x*x*x,gy=2*y-6*y*y+4*y*y*y;
    return Eigen::Vector2d{32*f*gy,-32*fx*g};
}
template<class F>void cellIntegral(Vector2D center,Vector2D h,F&& f){
    const auto& g=Quantum::Geometry::Parametric::gaussEight();
    for(int j=0;j<8;j++)for(int i=0;i<8;i++)f(Vector2D{center.x+.5*h.x*g.x[i],center.y+.5*h.y*g.x[j]},.25*h.x*h.y*g.w[i]*g.w[j]);
}
struct Fixture {
    QuadtreeMesh mesh;
    FaceSegmentGraph graph;
    NativeLeafDomain domain;
    Space space;
    Eigen::VectorXd initial,force;
    Fixture(int resolution=8):mesh(input(resolution)){
        DictionaryMesh::RefineBox outer;outer.min={.25,.25};outer.max={.75,.75};outer.level=1;
        DictionaryMesh::RefineBox inner;inner.min={.4375,.375};inner.max={.625,.5625};inner.level=2;mesh.applyRefineBoxes({outer,inner});
        graph.build(mesh);domain=NativeLeafDomain::build(mesh,graph,false);space=Space::build(mesh,domain);
        initial=Eigen::VectorXd::Zero(space.mass().rows());force=initial;
        const auto& g=Quantum::Geometry::Parametric::gaussEight();
        for(const auto face:domain.faces(mesh)){
            const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
            const auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();const auto& info=space.faces()[face.nativeRow];
            for(int k=0;k<8;k++){const double t=.5*(1+g.x[k]);const Vector2D p{a.x+t*(b.x-a.x),a.y+t*(b.y-a.y)};
                const double normal=profile(p)[face.horizontal?1:0],theta=((face.horizontal?p.x:p.y)-info.tangentCenter)/info.tangentScale;
                initial[info.first]+=.5*g.w[k]*normal;initial[info.first+1]+=.5*g.w[k]*theta*normal;}
        }
        for(int c=0;c<mesh.cellIndexNum();c++){
            const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto center=cell.centerCoord(),h=cell.size();const int row=space.firstCellVelocity(c);
            Eigen::Matrix<double,6,1> integrals=Eigen::Matrix<double,6,1>::Zero();
            cellIntegral(center,h,[&](Vector2D p,double da){const double x=(p.x-center.x)/h.x,y=(p.y-center.y)/h.y;const auto u=profile(p);
                initial[row]+=da/(h.x*h.y)*u[0];initial[row+1]+=da/(h.x*h.y)*u[1];initial[row+2]+=da/(h.x*h.y)*(-y*u[0]+x*u[1]);
                integrals+=da*Space::affineVelocityBasis({x,y},h).transpose()*u;});
            const auto& block=space.cells()[c];const Eigen::VectorXd local=block.affineProjection.transpose()*integrals;
            for(int j=0;j<local.size();j++)force[block.rows[j]]+=local[j];
        }
    }
    static DictionaryMesh input(int n){DictionaryMesh d;d.minCoord={0,0};d.maxCoord={1,1};d.Nx=n;d.Ny=5*n/4;d.maxCellNum=40000;d.maxInitialIterations=0;return d;}
    Flow owner(int mode){
        using Material=Quantum::Hydrodynamics::HydroCompatibleMomentMaterial;
        Material::ViscosityField mu(mesh.cellIndexNum(),Eigen::Matrix3d::Zero());
        for(int c=0;c<mesh.cellIndexNum();c++){
            const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto center=cell.centerCoord(),h=cell.size();
            for(int i=0;i<3;i++)for(int j=0;j<3;j++){
                const Vector2D p{center.x+.5*(i-1)*h.x,center.y+.5*(j-1)*h.y};
                mu[c](i,j)=mode==0?0:mode==3?1e-20:mode==2?(center.x<.5?.02:0):.02*(1+p.x+p.y);
            }
        }
        return Flow::createFrozenMaterial(mesh,domain,space,Material::buildSinglePhase(mesh,domain,space,1.3,mu),initial);
    }
};
struct Result {HydroCompleteObservation::Velocity velocity;double energyDefect=0;};
Result trajectory(Fixture& fixture,int mode,double dt,int steps,std::ofstream* history){
    auto flow=fixture.owner(mode);const double e0=flow.kineticEnergy();
    double cumulativeD=0,cumulativeW=0,cumulativeP=0,maxBalance=0,maxDivergence=0,maxEquation=0,maxTrace=0;int factors=0;
    int active=-1,inactive=-1;
    for(int i=0;i<steps;i++){
        HydroStepTime time;time.index=i+1;time.begin=flow.state().time;time.dt=dt;time.policy=HydroTimePolicy::ConservativeOneStep;
        const double midpoint=time.begin+.5*dt;
        const auto load=mode==0?flow.makeBodyForceLoad(fixture.mesh,fixture.domain,[&](Vector2D p){
            const Eigen::Vector2d force=std::cos(3*midpoint)*profile(p);return Vector2D{force[0],force[1]};},midpoint)
            :flow.makeLoad(fixture.mesh,fixture.domain,Eigen::VectorXd::Zero(fixture.initial.size()),midpoint);
        auto candidate=flow.prepareStokesStep(fixture.mesh,fixture.domain,time,load);const auto report=candidate.report();
        require(flow.state().acceptedSteps==i&&flow.state().time==time.begin,"preparation mutated accepted moment state");
        require(candidate.state().stageTime&&*candidate.state().stageTime==midpoint&&candidate.state().time==time.begin+dt,"moment state time layers differ");
        require((candidate.state().previousVelocity-flow.state().velocity).norm()==0,"moment previous velocity history differs");
        require((candidate.state().previousPhysicalVelocity->coordinates()-flow.state().physicalVelocity->coordinates()).norm()==0,"complete previous velocity history differs");
        flow.accept(fixture.mesh,fixture.domain,std::move(candidate));
        cumulativeD+=report.viscousDissipation;cumulativeW+=report.externalWork;cumulativeP+=report.pressureWork;
        const double balance=std::abs(report.kineticAfter-e0+cumulativeD-cumulativeW-cumulativeP)/e0;
        maxBalance=std::max(maxBalance,balance);maxDivergence=std::max(maxDivergence,report.divergenceResidual);
        maxEquation=std::max(maxEquation,report.equationResidual);maxTrace=std::max(maxTrace,report.viscousTraceResidual);
        if(!report.reusedFactorization)factors++;active=report.activeTangentialDofs;inactive=report.inactiveTangentialDofs;
        if(mode!=0)require(report.kineticAfter<=report.kineticBefore+2e-12*e0,"unforced moment Stokes created kinetic energy");
        if(history)*history<<mode<<','<<dt<<','<<i+1<<','<<flow.state().time<<','<<midpoint<<','<<report.kineticAfter<<','<<cumulativeD<<','<<cumulativeW<<','<<balance<<','<<report.divergenceResidual<<'\n';
    }
    require(factors==1&&maxBalance<3e-10&&maxDivergence<2e-11,"moment time energy/divergence/factor reuse failed");
    Eigen::VectorXd balance=Eigen::VectorXd::Zero(fixture.mesh.cellIndexNum());double faceIntegralDefect=0;
    for(const auto face:fixture.domain.faces(fixture.mesh)){
        const auto trace=flow.normalPolynomial(fixture.mesh,fixture.domain,face.nativeRow,Flow::VelocityLevel::Midpoint);
        const auto previous=flow.normalPolynomial(fixture.mesh,fixture.domain,face.nativeRow,Flow::VelocityLevel::PreviousEndpoint);
        require(trace.time==*flow.state().stageTime&&previous.time==flow.state().previousTime,"published face velocity has wrong time layer");
        const auto& edge=face.horizontal?fixture.mesh.hEdges(face.edgeLabel):fixture.mesh.vEdges(face.edgeLabel);
        const auto a=fixture.mesh.verts(edge.edgeToVert[0]).coord(),b=fixture.mesh.verts(edge.edgeToVert[1]).coord();
        const double t0=face.horizontal?a.x:a.y,t1=face.horizontal?b.x:b.y;const auto& g=Quantum::Geometry::Parametric::gaussEight();
        double integrated=0;for(int k=0;k<8;k++)integrated+=.5*face.area*g.w[k]*trace.value(.5*(t0+t1)+.5*(t1-t0)*g.x[k]);
        faceIntegralDefect=std::max(faceIntegralDefect,std::abs(integrated-trace.volumeRate()));
        if(face.negativeCell>=0)balance[face.negativeCell]+=integrated;if(face.positiveCell>=0)balance[face.positiveCell]-=integrated;
    }
    double inventoryDefect=0;for(int c=0;c<balance.size();c++)inventoryDefect=std::max(inventoryDefect,dt*std::abs(balance[c])/fixture.domain.cellVolumes()[fixture.domain.cellRowByIndex()[c]]);
    require(faceIntegralDefect<2e-14&&inventoryDefect<2e-11,"published velocity does not preserve actual cell flux balance");
    if(mode==0)require(active==0&&flow.state().viscousTrace.size()==0&&cumulativeD==0,"zero viscosity retained artificial trace inertia");
    if(mode==2)require(cumulativeD>0&&flow.state().viscousTrace.size()==0,"mixed viscosity lost physical dissipation or introduced trace state");
    if(mode==3)require(cumulativeD>0&&flow.state().viscousTrace.size()==0,"small positive viscosity was discarded or introduced trace state");
    std::cout<<std::setprecision(17)<<"MOMENT_TIME_TRAJECTORY mode="<<mode<<" dt="<<dt<<" steps="<<steps<<" end_time="<<flow.state().time
        <<" initial_energy="<<e0<<" final_energy="<<flow.kineticEnergy()
        <<" dissipation="<<cumulativeD<<" external_work="<<cumulativeW<<" pressure_work="<<cumulativeP<<" maximum_energy_defect="<<maxBalance
        <<" maximum_divergence="<<maxDivergence<<" maximum_equation_residual="<<maxEquation<<" maximum_trace_residual="<<maxTrace
        <<" active_traces="<<active<<" inactive_traces="<<inactive<<" factorizations="<<factors<<" face_integral_defect="<<faceIntegralDefect
        <<" final_step_inventory_defect="<<inventoryDefect<<" convection=0\n";
    return {*flow.state().physicalVelocity,maxBalance};
}
void acceptance(Fixture& f){
    auto owner=f.owner(1);const auto initial=owner.state().velocity;
    const auto initialField=*owner.state().physicalVelocity;
    const auto unchanged=[&]{return (owner.state().velocity-initial).norm()==0
        &&(owner.state().physicalVelocity->coordinates()-initialField.coordinates()).norm()==0;};
    HydroStepTime time;time.index=1;time.begin=0;time.dt=.005;
    const auto load=owner.makeLoad(f.mesh,f.domain,Eigen::VectorXd::Zero(initial.size()),.0025);
    {auto discarded=owner.prepareStokesStep(f.mesh,f.domain,time,load);}
    require(owner.state().acceptedSteps==0&&unchanged()&&!owner.state().stageTime,"discarded candidate changed accepted state");
    bool residualRejected=false;Flow::Options strict;strict.relativeResidualTolerance=1e-30;
    try{owner.prepareStokesStep(f.mesh,f.domain,time,load,strict);}catch(const std::invalid_argument&){residualRejected=true;}
    require(residualRejected&&owner.state().acceptedSteps==0&&unchanged(),"failed equation left a partially advanced state");
    const auto wrongTime=owner.makeLoad(f.mesh,f.domain,Eigen::VectorXd::Zero(initial.size()),0);bool rejected=false;
    try{owner.prepareStokesStep(f.mesh,f.domain,time,wrongTime);}catch(const std::invalid_argument&){rejected=true;}require(rejected,"wrong force time accepted");
    auto other=f.owner(1);const auto foreign=other.makeLoad(f.mesh,f.domain,Eigen::VectorXd::Zero(initial.size()),.0025);rejected=false;
    try{owner.prepareStokesStep(f.mesh,f.domain,time,foreign);}catch(const std::invalid_argument&){rejected=true;}require(rejected,"foreign force source accepted");
    auto a=owner.prepareStokesStep(f.mesh,f.domain,time,load);auto stale=owner.prepareStokesStep(f.mesh,f.domain,time,load);
    owner.accept(f.mesh,f.domain,std::move(a));const auto accepted=owner.state().velocity,acceptedCoordinates=owner.state().physicalVelocity->coordinates();rejected=false;
    try{owner.accept(f.mesh,f.domain,std::move(stale));}catch(const std::invalid_argument&){rejected=true;}
    require(rejected&&owner.state().acceptedSteps==1&&(owner.state().velocity-accepted).norm()==0
        &&(owner.state().physicalVelocity->coordinates()-acceptedCoordinates).norm()==0,"stale candidate overwrote accepted state");
    HydroStepTime next;next.index=2;next.begin=owner.state().time;next.dt=.0025;
    const auto nextLoad=owner.makeLoad(f.mesh,f.domain,Eigen::VectorXd::Zero(initial.size()),next.begin+.5*next.dt);rejected=false;
    try{owner.prepareStokesStep(f.mesh,f.domain,next,nextLoad);}catch(const std::invalid_argument&){rejected=true;}require(rejected,"fixed time policy changed dt");
    auto variable=f.owner(0);
    for(double dt:{.007,.003,.011}){HydroStepTime t;t.index=variable.state().acceptedSteps+1;t.begin=variable.state().time;t.dt=dt;t.policy=HydroTimePolicy::ConservativeOneStep;
        const auto l=variable.makeLoad(f.mesh,f.domain,Eigen::VectorXd::Zero(initial.size()),t.begin+.5*dt);
        auto trial=variable.prepareStokesStep(f.mesh,f.domain,t,l);variable.accept(f.mesh,f.domain,std::move(trial));}
    require(variable.state().acceptedSteps==3&&(variable.state().velocity-initial).norm()<2e-11*initial.norm(),"variable-step inviscid resting trajectory differs");
    auto badInitial=initial;badInitial[f.space.firstCellVelocity(0)]+=.1;bool initialRejected=false;
    try{Flow::create(f.mesh,f.domain,f.space,Stress::buildCellwise(f.mesh,f.domain,f.space,std::vector<double>(f.mesh.cellIndexNum(),.02)),1.3,badInitial);}catch(const std::invalid_argument&){initialRejected=true;}
    require(initialRejected,"non-solenoidal initial moments were silently projected");
    initialRejected=false;
    try{Flow::create(f.mesh,f.domain,f.space,Stress::buildCellwise(f.mesh,f.domain,f.space,std::vector<double>(f.mesh.cellIndexNum(),.02)),1.3,1e-14*badInitial);}catch(const std::invalid_argument&){initialRejected=true;}
    require(initialRejected,"velocity rescaling bypassed incompressible initialization");
    next.dt=.005;const auto meshLoad=owner.makeLoad(f.mesh,f.domain,Eigen::VectorXd::Zero(initial.size()),next.begin+.5*next.dt);
    auto meshCandidate=owner.prepareStokesStep(f.mesh,f.domain,next,meshLoad);
    DictionaryMesh::RefineBox change;change.min={.02,.02};change.max={.08,.08};change.level=1;f.mesh.applyRefineBoxes({change});rejected=false;
    try{owner.accept(f.mesh,f.domain,std::move(meshCandidate));}catch(const std::exception&){rejected=true;}
    require(rejected&&owner.state().acceptedSteps==1&&(owner.state().velocity-accepted).norm()==0
        &&(owner.state().physicalVelocity->coordinates()-acceptedCoordinates).norm()==0,"mesh change accepted stale moment candidate");
    std::cout<<"MOMENT_TIME_ACCEPTANCE discard=1 solve_failure=1 force_time=1 foreign_source=1 stale_candidate=1 fixed_dt=1 variable_dt=1 non_solenoidal_initial=1 scaled_initial_rejection=1 changed_mesh=1\n";
}
}
void hydroMomentTimeContract(){
    Fixture f;std::ofstream history;
    if(const char* prefix=std::getenv("QUANTUM_MOMENT_TIME_OUTPUT")){history.open(std::string(prefix)+"_history.csv");history<<std::setprecision(17)
        <<"mode,dt,step,time,pressure_time,kinetic_energy,cumulative_dissipation,cumulative_external_work,relative_energy_balance,divergence\n";}
    for(int mode:{0,1}){
        const double coarse=mode==0?.02:.005;std::vector<Result> runs;
        for(int refinement=0;refinement<3;refinement++){const double dt=coarse/std::pow(2.,refinement);runs.push_back(trajectory(f,mode,dt,int(std::lround(.4/dt)),history?&history:nullptr));}
        const double a=std::sqrt(HydroCompleteObservation::distanceSquared(f.mesh,f.domain,runs[0].velocity,runs[1].velocity));
        const double b=std::sqrt(HydroCompleteObservation::distanceSquared(f.mesh,f.domain,runs[1].velocity,runs[2].velocity)),order=std::log2(a/b);
        std::cout<<std::setprecision(17)<<"MOMENT_TIME_ORDER mode="<<mode<<" coarse_difference="<<a<<" fine_difference="<<b<<" order="<<order<<'\n';
        require(order>1.8&&order<2.2,"moment midpoint not approaching second-order time convergence");
    }
    trajectory(f,2,.005,20,history?&history:nullptr);trajectory(f,3,.005,1,nullptr);
    acceptance(f);
}

namespace {
using Transport=Quantum::Hydrodynamics::HydroCompatibleMomentTransport;
Eigen::VectorXd localVector(const Space::CellBlock& b,const Eigen::VectorXd& global){
    Eigen::VectorXd r(b.rows.size());for(int i=0;i<r.size();i++)r[i]=global[b.rows[i]];return r;
}
double reconstruction(Fixture& f,int n){
    const auto transport=Transport::build(f.mesh,f.domain,f.space);
    std::vector<std::vector<int>> faces(f.mesh.cellIndexNum());
    for(const auto face:f.domain.faces(f.mesh)){if(face.negativeCell>=0)faces[face.negativeCell].push_back(face.nativeRow);if(face.positiveCell>=0)faces[face.positiveCell].push_back(face.nativeRow);}
    double affineError=0,divergenceError=0,momentError=0,traceError=0,l2=0,norm=0,constraint=0;
    for(int c=0;c<f.mesh.cellIndexNum();c++){
        const auto& cell=f.mesh.cells(f.mesh.cellIndexToLabel[c]);const auto center=cell.centerCoord(),h=cell.size();
        const auto& bulk=f.space.cells()[c];const auto& reconstruction=transport.cells()[c];constraint=std::max(constraint,reconstruction.constraintDefect);
        const Eigen::VectorXd initial=localVector(bulk,f.initial);
        Eigen::VectorXd arbitrary(initial.size());for(int k=0;k<arbitrary.size();k++)arbitrary[k]=std::sin(1.2*(k+1));
        Eigen::Vector3d div=Eigen::Vector3d::Zero();for(int k=0;k<3;k++)for(int j=0;j<arbitrary.size();j++)div[k]+=f.space.divergence().coeff(3*c+k,bulk.rows[j])*arbitrary[j]*(k?12:1)/(h.x*h.y);
        Eigen::Vector3d moments=Eigen::Vector3d::Zero();Eigen::Matrix<double,6,1> projected=Eigen::Matrix<double,6,1>::Zero();
        Eigen::Matrix<double,6,6> gram=Eigen::Matrix<double,6,6>::Zero();
        for(const auto& piece:reconstruction.pieces){
            const auto evaluate=[&](Vector2D x){return Eigen::Vector2d(Transport::rtBasis(x)*piece.coefficients*arbitrary);};
            cellIntegral(piece.center,piece.size,[&](Vector2D p,double da){
                const Vector2D local{(p.x-piece.center.x)/piece.size.x,(p.y-piece.center.y)/piece.size.y},q{(p.x-center.x)/h.x,(p.y-center.y)/h.y};
                const auto rt=Transport::rtBasis(local);const auto basis=Space::affineVelocityBasis(q,h);
                affineError=std::max(affineError,(rt*piece.coefficients*bulk.affineSamples-basis).norm());
                const auto u=evaluate(local);moments+=da/(h.x*h.y)*Eigen::Vector3d(u[0],u[1],-q.y*u[0]+q.x*u[1]);
                projected+=da*basis.transpose()*u;gram+=da*basis.transpose()*basis;
                // RT normal components are quadratic in their own coordinate:
                // centered differences are exact polynomials, independently
                // of the derivative coefficients in the constraint builder.
                const double d=(evaluate({local.x+.125,local.y})[0]-evaluate({local.x-.125,local.y})[0])/(.25*piece.size.x)
                    +(evaluate({local.x,local.y+.125})[1]-evaluate({local.x,local.y-.125})[1])/(.25*piece.size.y);
                divergenceError=std::max(divergenceError,std::abs(d-div.dot(Eigen::Vector3d(1,q.x,q.y)))*std::min(h.x,h.y));
                l2+=da*(rt*piece.coefficients*initial-profile(p)).squaredNorm();norm+=da*profile(p).squaredNorm();
            });
            for(int axis=0;axis<2;axis++)for(int sign:{-1,1})for(double t:{-.37,0.,.29}){
                const Vector2D q=axis==0?Vector2D{.5*sign,t}:Vector2D{t,.5*sign};
                const Vector2D p{piece.center.x+piece.size.x*q.x,piece.center.y+piece.size.y*q.y};const double coordinate=axis==0?p.x:p.y;
                const bool outside=std::abs(coordinate-(axis==0?center.x:center.y)-.5*sign*(axis==0?h.x:h.y))<1e-12;
                double target=0;bool found=false;
                if(outside){for(int native:faces[c]){const auto face=f.domain.face(f.mesh,native);
                    if(face.horizontal!=(axis==1)||(face.negativeCell!=c&&face.positiveCell!=c)||(face.negativeCell==c?1:-1)!=sign)continue;
                    const auto& info=f.space.faces()[face.nativeRow];const double s=axis==0?p.y:p.x;
                    if(std::abs(s-info.tangentCenter)>.5*face.area+1e-12)continue;
                    int k=0;while(bulk.rows[k]!=info.first)k++;
                    target=arbitrary[k]+arbitrary[k+1]*(s-info.tangentCenter)/info.tangentScale;found=true;break;
                }}else for(const auto& other:reconstruction.pieces){
                    if(&other==&piece)continue;const Vector2D r{(p.x-other.center.x)/other.size.x,(p.y-other.center.y)/other.size.y};
                    if(std::abs((axis==0?r.x:r.y)+.5*sign)<1e-12&&std::abs(axis==0?r.y:r.x)<.5){target=(Transport::rtBasis(r)*other.coefficients*arbitrary)[axis];found=true;break;}
                }
                require(found,"RT trace test could not find face/piece");traceError=std::max(traceError,std::abs(evaluate(q)[axis]-target));
            }
        }
        momentError=std::max(momentError,(moments-arbitrary.tail(3)).norm());
        momentError=std::max(momentError,(projected-gram*bulk.affineProjection*arbitrary).norm()/(h.x*h.y));
    }
    const double error=std::sqrt(l2/norm);
    std::cout<<std::setprecision(17)<<"MOMENT_RECONSTRUCTION n="<<n<<" cells="<<f.mesh.cellIndexNum()<<" affine="<<affineError<<" divergence="<<divergenceError<<" moments="<<momentError<<" trace="<<traceError<<" constraints="<<constraint<<" relative_l2="<<error<<'\n';
    require(std::max({affineError,divergenceError,momentError,traceError})<3e-10,"RT polynomial/divergence/moment/trace contract failed");return error;
}
void convectionAlgebra(Fixture& f){
    auto transport=Transport::build(f.mesh,f.domain,f.space);auto a=transport.assemble(f.mesh,f.domain,f.space,f.initial,1.3,Transport::Flux::TraceUpwind);
    const auto baseline=transport.assemble(f.mesh,f.domain,f.space,f.initial,1.3);
    require(baseline.dissipation.norm()==0 && (baseline.central-a.central).norm()==0
        && Flow::Options{}.convectionFlux==Transport::Flux::Central,"central baseline contains global artificial dissipation");
    Eigen::VectorXd x=Eigen::VectorXd::Zero(f.initial.size()),y=x,affine=x;
    for(const auto face:f.domain.faces(f.mesh)){
        const auto& b=f.space.faces()[face.nativeRow];const auto p=face.center;const int normal=face.horizontal?1:0;
        (normal?y:x)[b.first]=1;
        affine[b.first]=normal?-.1+.5*p.x+.7*p.y:1+.2*p.x-.3*p.y;
        affine[b.first+1]=(normal?.5:-.3)*b.tangentScale;
    }
    for(int c=0;c<f.mesh.cellIndexNum();c++){
        const auto& cell=f.mesh.cells(f.mesh.cellIndexToLabel[c]);const auto p=cell.centerCoord(),h=cell.size();const int row=f.space.firstCellVelocity(c);
        x[row]=y[row+1]=1;affine.segment<3>(row)<<1+.2*p.x-.3*p.y,-.1+.5*p.x+.7*p.y,(.3*h.y+.5*h.x)/12;
    }
    const double skew=(a.central+Space::Matrix(a.central.transpose())).norm()/a.central.norm();
    const double constants=std::max({(a.central*x).norm(),(a.central*y).norm(),(a.central.transpose()*x).norm(),(a.central.transpose()*y).norm()})/a.central.norm();
    const double stabilization=(a.dissipation*affine).norm()/a.dissipation.norm();
    const double symmetry=(a.dissipation-Space::Matrix(a.dissipation.transpose())).norm()/a.dissipation.norm();
    double minimumWork=1e100;for(int k=0;k<8;k++){Eigen::VectorXd v(x.size());for(int i=0;i<v.size();i++)v[i]=std::sin((i+1)*(k+.31));minimumWork=std::min(minimumWork,v.dot(a.dissipation*v));}
    Eigen::Vector2d sum=Eigen::Vector2d::Zero();double scale=0;for(const auto& v:a.cellBoundaryMomentum){sum+=v;scale+=v.norm();}
    std::cout<<std::setprecision(17)<<"MOMENT_CONVECTION_ALGEBRA skew="<<skew<<" constant="<<constants<<" affine_stabilization="<<stabilization<<" symmetry="<<symmetry<<" minimum_sample_work="<<minimumWork<<" local_momentum="<<a.maximumLocalMomentumDefect<<" global_flux="<<sum.norm()/scale<<'\n';
    require(std::max({skew,constants,stabilization,symmetry})<2e-11&&minimumWork>=0&&a.maximumLocalMomentumDefect<1e-12&&sum.norm()<1e-12*scale,"conservative convection algebra failed");
}
Result nonlinearTrajectory(Fixture& f,int mode,double dt,int steps,std::ofstream* history,const char* prefix){
    auto flow=f.owner(mode==2?1:0);Flow::Options options;options.relativeResidualTolerance=1e-12;
    const auto initialField=*flow.state().physicalVelocity;
    options.convectionFlux=Transport::Flux::Central;
    const double e0=flow.kineticEnergy();double dissipated=0,viscous=0,central=0,maxBalance=0,maxEquation=0,maxDiv=0,maxLocal=0;int maximumIterations=0;
    for(int step=1;step<=steps;step++){
        HydroStepTime t;t.index=step;t.begin=flow.state().time;t.dt=dt;
        auto load=flow.makeLoad(f.mesh,f.domain,Eigen::VectorXd::Zero(f.initial.size()),t.begin+.5*dt);
        auto candidate=flow.prepareNavierStokesStep(f.mesh,f.domain,t,load,options);const auto r=candidate.report();
        flow.accept(f.mesh,f.domain,std::move(candidate));dissipated+=r.transportDissipation;viscous+=r.viscousDissipation;central+=r.centralConvectionWork;
        const double defect=std::abs(r.kineticAfter-e0+dissipated+viscous)/e0;maxBalance=std::max(maxBalance,defect);
        maxEquation=std::max(maxEquation,r.equationResidual);maxDiv=std::max(maxDiv,r.divergenceResidual);maxLocal=std::max(maxLocal,r.localMomentumDefect);maximumIterations=std::max(maximumIterations,r.nonlinearIterations);
        require(r.transportDissipation>=-1e-13*e0&&r.kineticAfter<=r.kineticBefore+2e-11*e0,"unforced nonlinear flow creates kinetic energy");
        if(history)*history<<mode<<','<<dt<<','<<step<<','<<flow.state().time<<','<<r.kineticAfter<<','<<dissipated<<','<<viscous<<','<<defect<<','<<r.equationResidual<<','<<r.nonlinearIterations<<'\n';
    }
    const double change=std::sqrt(HydroCompleteObservation::distanceSquared(f.mesh,f.domain,*flow.state().physicalVelocity,initialField)
        /HydroCompleteObservation::normSquared(f.mesh,f.domain,initialField));
    require(maxBalance<2e-10&&change>1e-3&&std::abs(central)<1e-12*e0,"nonlinear midpoint energy balance or actual evolution failed");
    std::cout<<std::setprecision(17)<<"MOMENT_NS_TRAJECTORY mode="<<mode<<" dt="<<dt<<" steps="<<steps<<" end_time="<<flow.state().time<<" initial_energy="<<e0<<" final_energy="<<flow.kineticEnergy()<<" transport_dissipation="<<dissipated<<" viscous_dissipation="<<viscous<<" central_work="<<central<<" maximum_balance="<<maxBalance<<" maximum_equation="<<maxEquation<<" maximum_divergence="<<maxDiv<<" maximum_local_momentum="<<maxLocal<<" maximum_iterations="<<maximumIterations<<" velocity_change="<<change<<'\n';
    if(prefix&&dt==.005){std::ofstream field(std::string(prefix)+"_field_"+std::to_string(mode)+".csv");field<<std::setprecision(17)<<"x,y,hx,hy,initial_ux,initial_uy,ux,uy,p\n";
        for(int c=0;c<f.mesh.cellIndexNum();c++){const auto& cell=f.mesh.cells(f.mesh.cellIndexToLabel[c]);const auto p=cell.centerCoord(),h=cell.size();
            const auto initial=initialField.firstJet(p).value,u=flow.state().physicalVelocity->firstJet(p).value;
            field<<p.x<<','<<p.y<<','<<h.x<<','<<h.y<<','<<initial.x<<','<<initial.y<<','<<u.x<<','<<u.y<<','<<flow.state().pressure[3*c]<<'\n';}}
    return {*flow.state().physicalVelocity,maxBalance};
}
}
void hydroMomentTransportContract(){
    Fixture f;double previous=0;
    for(int n:{8,16,32}){double error=0;if(n==8)error=reconstruction(f,n);else{Fixture fine(n);error=reconstruction(fine,n);}
        if(previous){const double order=std::log2(previous/error);std::cout<<"MOMENT_RECONSTRUCTION_ORDER n="<<n<<" order="<<order<<'\n';require(order>1.8,"RT reconstructed physical field lost second-order accuracy");}previous=error;}
    convectionAlgebra(f);std::ofstream history;const char* prefix=std::getenv("QUANTUM_MOMENT_TRANSPORT_OUTPUT");
    if(prefix){history.open(std::string(prefix)+"_history.csv");history<<std::setprecision(17)<<"mode,dt,step,time,kinetic_energy,transport_dissipation,viscous_dissipation,relative_balance,equation_residual,iterations\n";}
    for(int mode:{0,2}){std::vector<Result> runs;
        for(double dt:{.02,.01,.005})runs.push_back(nonlinearTrajectory(f,mode,dt,int(std::lround(.2/dt)),history?&history:nullptr,prefix));
        const double first=std::sqrt(HydroCompleteObservation::distanceSquared(f.mesh,f.domain,runs[0].velocity,runs[1].velocity));
        const double second=std::sqrt(HydroCompleteObservation::distanceSquared(f.mesh,f.domain,runs[1].velocity,runs[2].velocity)),order=std::log2(first/second);
        std::cout<<std::setprecision(17)<<"MOMENT_NS_ORDER mode="<<mode<<" coarse_difference="<<first<<" fine_difference="<<second<<" order="<<order<<'\n';require(order>1.8&&order<2.2,"nonlinear midpoint time convergence differs from second order");
    }
    auto flow=f.owner(0);HydroStepTime t;t.index=1;t.begin=0;t.dt=.02;
    const auto before=flow.state();
    const auto load=flow.makeLoad(f.mesh,f.domain,Eigen::VectorXd::Zero(f.initial.size()),.01);Flow::Options fail;fail.maximumIterations=1;fail.relativeResidualTolerance=1e-12;
    bool rejected=false;try{flow.prepareNavierStokesStep(f.mesh,f.domain,t,load,fail);}catch(const std::invalid_argument&){rejected=true;}
    require(rejected&&flow.state().time==0&&flow.state().acceptedSteps==0&&(flow.state().velocity-before.velocity).norm()==0
        &&(flow.state().physicalVelocity->coordinates()-before.physicalVelocity->coordinates()).norm()==0,"nonlinear rejection changed accepted physical state");
    auto stokes=flow.prepareStokesStep(f.mesh,f.domain,t,load);flow.accept(f.mesh,f.domain,std::move(stokes));
    require((flow.state().velocity-f.initial).norm()<1e-10*f.initial.norm(),"nonlinear factorization contaminated Stokes cache");
    std::cout<<"MOMENT_NS_ACCEPTANCE iteration_failure=1 old_state_retained=1 stokes_cache_isolated=1\n";
}
void hydroMomentPhysicalInputContract(){
    Fixture f;const auto sampled=f.space.interpolateVelocity(f.mesh,f.domain,profile),load=f.space.integrateBodyForce(f.mesh,f.domain,profile);
    const double velocity=(sampled-f.initial).norm()/f.initial.norm(),force=(load-f.force).norm()/f.force.norm();
    require(velocity<2e-14&&force<2e-14,"physical moment ports differ from independent integrals");
    bool rejected=false;try{f.space.interpolateVelocity(f.mesh,f.domain,[](Vector2D){return Eigen::Vector2d(std::numeric_limits<double>::quiet_NaN(),0);});}catch(const std::invalid_argument&){rejected=true;}
    require(rejected,"nonfinite initial physical vector was accepted");
    auto resting=Flow::create(f.mesh,f.domain,f.space,Stress::buildCellwise(f.mesh,f.domain,f.space,std::vector<double>(f.mesh.cellIndexNum(),.02)),1.3,Eigen::VectorXd::Zero(f.initial.size()));
    const Eigen::VectorXd gravity=f.space.integrateBodyForce(f.mesh,f.domain,[](Vector2D){return Eigen::Vector2d{1.3*.7,-1.3*.2};});double pressureError=0,restingVelocity=0,fullEquation=0;
    for(int step=1;step<=3;step++){
        HydroStepTime t;t.index=step;t.begin=resting.state().time;t.dt=.01;
        auto input=resting.makeUniformAccelerationLoad(f.mesh,f.domain,{.7,-.2},t.begin+.5*t.dt);
        auto candidate=resting.prepareNavierStokesStep(f.mesh,f.domain,t,input);fullEquation=std::max(fullEquation,candidate.report().equationResidual);
        resting.accept(f.mesh,f.domain,std::move(candidate));
        restingVelocity=std::max(restingVelocity,std::sqrt(HydroCompleteObservation::normSquared(f.mesh,f.domain,*resting.state().physicalVelocity)));
        require(restingVelocity<2e-13,"uniform gravity generated a false hydrostatic velocity");
        for(int c=0;c<f.mesh.cellIndexNum();c++){
            const auto& cell=f.mesh.cells(f.mesh.cellIndexToLabel[c]);const auto center=cell.centerCoord(),h=cell.size();
            cellIntegral(center,h,[&](Vector2D p,double){const Eigen::Vector3d q(1,(p.x-center.x)/h.x,(p.y-center.y)/h.y);
                pressureError=std::max(pressureError,std::abs(resting.state().pressure.segment<3>(3*c).dot(q)-1.3*(.7*(p.x-.5)-.2*(p.y-.5))));});
        }
    }
    // Dirichlet normal rows have wall reactions, not momentum equations.
    // Compare the physical force on every admissible velocity test function.
    Eigen::VectorXd difference=f.space.pressureLoad()*resting.state().pressure-gravity,admissibleForce=gravity;
    for(const auto face:f.domain.faces(f.mesh))if(face.outwardSign){const int row=f.space.faces()[face.nativeRow].first;difference.segment<2>(row).setZero();admissibleForce.segment<2>(row).setZero();}
    const double balance=difference.norm()/admissibleForce.norm();
    require(pressureError<1e-14&&balance<2e-13,"physical hydrostatic pressure/force does not match independent integrals");
    auto a=f.owner(1),b=f.owner(1);HydroStepTime t;t.index=1;t.begin=0;t.dt=.01;
    const auto none=a.makeLoad(f.mesh,f.domain,Eigen::VectorXd::Zero(f.initial.size()),.005),withGravity=b.makeUniformAccelerationLoad(f.mesh,f.domain,{.7,-.2},.005);
    auto ca=a.prepareNavierStokesStep(f.mesh,f.domain,t,none),cb=b.prepareNavierStokesStep(f.mesh,f.domain,t,withGravity);
    const double movingDifference=std::sqrt(HydroCompleteObservation::distanceSquared(f.mesh,f.domain,*ca.state().physicalVelocity,*cb.state().physicalVelocity));
    require(movingDifference<2e-13&&(cb.state().pressure-ca.state().pressure-resting.state().pressure).norm()<1e-13,"gradient force changed nonzero flow dynamics or physical pressure lift");
    std::cout<<std::setprecision(17)<<"MOMENT_PHYSICAL_INPUT velocity="<<velocity<<" force="<<force<<" nonfinite_rejected=1 hydrostatic_velocity="<<restingVelocity
        <<" pressure_error="<<pressureError<<" hydrostatic_force_balance="<<balance<<" full_equation="<<fullEquation<<" moving_gradient_difference="<<movingDifference<<'\n';
}
