#include "Modules/Hydrodynamics/Hydrodynamics.h"
#include "Modules/Hydrodynamics/HydroConservativePressure.h"
#include "Modules/Hydrodynamics/HydroConservativeTransport.h"
#include "Modules/Hydrodynamics/HydroUniformLoad.h"
#include "Modules/Interface/Compatibility/InterfaceConservativePhase.h"
#include "Modules/Interface/Compatibility/InterfaceParametricProjection.h"
#include "Modules/Interface/Compatibility/InterfacePrescribedData.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include <cmath>
#include <iostream>
#include <iomanip>

namespace {
void require(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}

void compatiblePressure(bool rz,double liquidDensity) {
    using namespace Quantum::Hydrodynamics;
    constexpr int n=16;constexpr double radius=.213,sigma=.7;
    DictionaryMesh input;input.minCoord={0,0};input.maxCoord={1,1};input.Nx=input.Ny=n;input.maxCellNum=4096;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,rz);
    const auto pressure=HydroConservativePressure::buildDivergenceAdjoint(mesh,domain);
    const auto& d=pressure.integratedDivergence();const auto& c=pressure.gradientLoad();
    const Eigen::SparseMatrix<double> adjoint=c+HydroConservativePressure::Matrix(d.transpose());
    require(adjoint.norm()==0,"CCH pressure is not the exact negative divergence transpose");
    Eigen::VectorXd velocity(domain.nativeFaceCount()),polynomial(mesh.cellIndexNum());
    std::vector<double> scalar(mesh.cellIndexNum()),boundary(domain.nativeFaceCount(),0);
    for(int i=0;i<mesh.cellIndexNum();++i){const auto point=mesh.cells(mesh.cellIndexToLabel[i]).centerCoord();
        polynomial[i]=scalar[i]=1+.3*point.x-.7*point.y;}
    for(const auto f:domain.faces(mesh)) {
        velocity[f.nativeRow]=f.outwardSign?0:std::sin(.17*(f.nativeRow+1));
        boundary[f.nativeRow]=1+.3*f.center.x-.7*f.center.y;
    }
    const Eigen::VectorXd gradient=c*polynomial,divergence=d*velocity;
    const double pressureWork=velocity.dot(gradient),divergenceWork=divergence.dot(polynomial);
    const double powerError=std::abs(pressureWork+divergenceWork)/std::max(1.,std::abs(pressureWork)+std::abs(divergenceWork));
    const auto loads=pressure.apply(mesh,domain,scalar,boundary);double affineError=0;
    for(const auto f:domain.faces(mesh))if(f.outwardSign==0)
        affineError=std::max(affineError,std::abs(loads[f.nativeRow]/pressure.duals().volumes()[f.nativeRow].measure-(f.horizontal?-.7:.3)));
    require(powerError<5e-14 && affineError<5e-12,"compatible pressure lost work identity or uniform-grid consistency");

    std::vector<long double> volumes(mesh.cellIndexNum());std::vector<double> alpha(mesh.cellIndexNum());
    if(!rz) {
        // One accepted parameter curve initializes the inventory. No nodal
        // level set, PLIC reconstruction, or second force geometry is built.
        const auto curve=InterfacePrescribedData::parametricEllipse({.47,.51},{radius,radius},32);
        const auto geometry=InterfaceParametricProjection::initialize(curve,{n,n});
        for(int cell=0;cell<mesh.cellIndexNum();++cell){const auto p=mesh.cells(mesh.cellIndexToLabel[cell]).centerCoord();
            const int i=static_cast<int>(p.x*n),j=static_cast<int>(p.y*n);
            volumes[cell]=geometry.inventory()[i*n+j];}
    } else {
        // Prescribed sphere-shaped phase profile, for the RZ FORCE algebra.
        // This does not claim that the new parameter transport supports RZ.
        for(int cell=0;cell<mesh.cellIndexNum();++cell){const auto p=mesh.cells(mesh.cellIndexToLabel[cell]).centerCoord();
            const double fraction=std::clamp(.5-(std::hypot(p.x,p.y-.5)-radius)*n,0.,1.);
            volumes[cell]=static_cast<long double>(fraction)*domain.cellVolumes()[domain.cellRowByIndex()[cell]];}
    }
    auto phase=InterfaceConservativePhase::initialize(mesh,domain,volumes);
    for(int cell=0;cell<mesh.cellIndexNum();++cell)alpha[cell]=static_cast<double>(phase.liquidVolumes()[cell]/phase.capacities()[cell]);
    const double jump=sigma*(rz?2:1)/radius;
    std::vector<double> coefficient(domain.nativeFaceCount());
    for(const auto f:domain.faces(mesh))coefficient[f.nativeRow]=f.outwardSign?0:jump;
    const auto load=buildCompatibleScalarGradientLoad(mesh,domain,pressure,alpha,coefficient);
    const auto kind=HydroStepGeometry::complete(rz);Hydrodynamics hydro(mesh);
    hydro.initializeState(mesh,kind);HydroSolveOptions options;options.compatibleNativePressure=true;options.linearResidualTolerance=1e-10;
    double time=0,maxVelocity=0,maxPressureError=0,maxPressurePower=0,maxGauge=0,maxResidual=0;
    long double meanFraction=0,totalVolume=0;
    for(int cell=0;cell<mesh.cellIndexNum();++cell){meanFraction+=phase.liquidVolumes()[cell];totalVolume+=phase.capacities()[cell];}
    meanFraction/=totalVolume;
    for(int step=1;step<=8;++step) {
        const double dt=step%2?.001:.002;
        const InterfaceStepTime interval{step,time,time+dt,dt};
        hydro.updateMaterials(mesh,{kind,time,HydroMaterialEvolution::TimeDependent,
            HydroConservativeMaterial{&domain,&phase,liquidDensity,1,.1,.001}});
        const auto next=InterfaceConservativePhase::advance(mesh,domain,phase,interval,{});
        const auto mass=HydroConservativeTransport::build(mesh,domain,phase,next,interval,liquidDensity,1);
        const auto provenance=hydro.captureForceSource(time,HydroForceEvolution::TimeDependent);
        const HydroForceDensity force{nullptr,nullptr,&load.faceForceX,&load.faceForceY,time,provenance};
        const HydroStepInput request{kind,{step,time,dt,HydroTimePolicy::ConservativeOneStep},options,force,nullptr,&mass};
        if(step==2) {
            auto wrongOptions=options;wrongOptions.compatibleNativePressure=false;
            const auto rejected=hydro.advance(mesh,{kind,request.time,wrongOptions,force,nullptr,&mass});
            require(!rejected.success && rejected.stateChange==HydroStateChange::Unchanged && hydro.state().progress.completedSteps==1,
                "pressure discretization changed inside an accepted history");
        }
        const auto report=hydro.advance(mesh,request);
        if(!report.success)throw std::runtime_error("compatible static load: "+report.stage+": "+report.reason);
        require(report.compatibleNativePressureUsed && report.physicalNativeMomentumUsed && report.conservativeMomentumUsed
            && report.coupledUsedSuperLU && !report.momentumUsedSuperLU,"compatible step bypassed the direct native coupled solve");
        const auto state=hydro.state();long double pressureIntegral=0;
        for(int i=0;i<state.velocityX.size();++i)maxVelocity=std::max(maxVelocity,std::abs(state.velocityX[i]));
        for(int i=0;i<state.velocityY.size();++i)maxVelocity=std::max(maxVelocity,std::abs(state.velocityY[i]));
        for(int cell=0;cell<mesh.cellIndexNum();++cell) {
            const double exact=jump*(alpha[cell]-static_cast<double>(meanFraction));
            maxPressureError=std::max(maxPressureError,std::abs(state.physicalPressure[cell]-exact)/jump);
            pressureIntegral+=phase.capacities()[cell]*state.physicalPressure[cell];
        }
        maxGauge=std::max(maxGauge,std::abs(static_cast<double>(pressureIntegral/totalVolume))/jump);
        maxPressurePower=std::max(maxPressurePower,std::abs(report.physicalPressurePower));
        maxResidual=std::max(maxResidual,report.momentumBalanceRelativeResidual);
        phase=next;time+=dt;
    }
    require(maxVelocity<1e-11 && maxPressureError<1e-10 && maxGauge<5e-14 && maxPressurePower<1e-12,
        "compatible analytic-curvature equilibrium lost Laplace pressure, rest state, pressure work or gauge");
    std::cout<<std::setprecision(17)<<"HYDRO_COMPATIBLE_PRESSURE metric="<<(rz?"RZ":"XY")<<" density_ratio="<<liquidDensity
        <<" steps=8 end_time="<<time<<" adjoint_defect="<<adjoint.norm()<<" power_identity="<<powerError<<" affine_error="<<affineError
        <<" max_native_velocity="<<maxVelocity<<" laplace_pressure_relative="<<maxPressureError<<" pressure_mean_relative="<<maxGauge
        <<" pressure_power="<<maxPressurePower<<" momentum_residual="<<maxResidual
        <<" analytic_curvature=1 interface_advection=0 amr_accuracy_claim=0\n";
}

void compatibleSolenoidalKick(bool rz) {
    using namespace Quantum::Hydrodynamics;
    DictionaryMesh input;input.minCoord={0,0};input.maxCoord={1,1};input.Nx=input.Ny=8;input.maxCellNum=1024;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,rz);
    const double pi=std::acos(-1.),dt=.003;
    // Integrating this streamfunction over actual faces gives an independent
    // nonzero divergence-free target, including RZ face areas and the axis.
    const auto psi=[&](Vector2D p){return .1*(rz?2*pi:1)*p.x*p.x*(1-p.x)*(1-p.x)*p.y*p.y*(1-p.y)*(1-p.y);};
    VEdgeField<double> forceX(mesh,0),targetX(mesh,0);HEdgeField<double> forceY(mesh,0),targetY(mesh,0);
    for(const auto face:domain.faces(mesh)) {
        const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
        const auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
        const double value=face.outwardSign?0:(face.horizontal?-1:1)*(psi(b)-psi(a))/face.area;
        if(face.horizontal){targetY.value[edge.index()]=value;forceY.value[edge.index()]=value/dt;}
        else{targetX.value[edge.index()]=value;forceX.value[edge.index()]=value/dt;}
    }
    const auto kind=HydroStepGeometry::complete(rz);
    const auto phase=InterfaceConservativePhase::initialize(mesh,domain,std::vector<long double>(mesh.cellIndexNum(),0));
    const InterfaceStepTime time{1,0,dt,dt};const auto next=InterfaceConservativePhase::advance(mesh,domain,phase,time,{});
    const auto mass=HydroConservativeTransport::build(mesh,domain,phase,next,time,1,1);
    Hydrodynamics hydro(mesh);hydro.initializeState(mesh,kind);
    hydro.updateMaterials(mesh,{kind,0,HydroMaterialEvolution::TimeDependent,HydroConservativeMaterial{&domain,&phase,1,1,0,0}});
    HydroSolveOptions options;options.compatibleNativePressure=true;options.linearResidualTolerance=1e-10;
    const HydroForceDensity force{nullptr,nullptr,&forceX,&forceY,0,hydro.captureForceSource(0,HydroForceEvolution::TimeDependent)};
    const auto report=hydro.advance(mesh,{kind,{1,0,dt,HydroTimePolicy::ConservativeOneStep},options,force,nullptr,&mass});
    if(!report.success)throw std::runtime_error("compatible solenoidal kick: "+report.reason);
    const auto state=hydro.state();double velocityError=0,targetMagnitude=0,pressureError=0;
    for(int i=0;i<targetX.size();++i){velocityError=std::max(velocityError,std::abs(state.velocityX[i]-targetX.value[i]));targetMagnitude=std::max(targetMagnitude,std::abs(targetX.value[i]));}
    for(int i=0;i<targetY.size();++i){velocityError=std::max(velocityError,std::abs(state.velocityY[i]-targetY.value[i]));targetMagnitude=std::max(targetMagnitude,std::abs(targetY.value[i]));}
    for(int i=0;i<state.physicalPressure.size();++i)pressureError=std::max(pressureError,std::abs(state.physicalPressure[i]));
    require(report.compatibleNativePressureUsed && targetMagnitude>1e-4 && velocityError<2e-12*targetMagnitude && pressureError<2e-12
        && report.divergenceAfterProjection<1e-12,"compatible saddle solve changed the exact solenoidal force response");
    std::cout<<std::setprecision(17)<<"HYDRO_COMPATIBLE_KICK metric="<<(rz?"RZ":"XY")<<" target_speed="<<targetMagnitude
        <<" velocity_error="<<velocityError<<" pressure_error="<<pressureError<<" divergence="<<report.divergenceAfterProjection<<'\n';
}
}

void hydroCompatiblePressureContract(){for(bool rz:{false,true}){for(double density:{1.,1000.})compatiblePressure(rz,density);compatibleSolenoidalKick(rz);}}
