#include "Modules/Hydrodynamics/HydroConservativeStress.h"
#include "Modules/Hydrodynamics/HydroConservativePressure.h"
#include "Modules/Hydrodynamics/Hydrodynamics.h"
#include "Modules/Interface/Compatibility/InterfacePrescribedData.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include <Eigen/SparseCholesky>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>

namespace {
using namespace Quantum::Hydrodynamics;
namespace P=Quantum::Geometry::Parametric;
void require(bool ok,const char* reason){if(!ok)throw std::runtime_error(reason);}
constexpr double mu0=.02;
const double pi=std::acos(-1.);
// Analytic, divergence-free meridional flow from
// psi=r^2(1-r^2)^3 sin(pi*z). It is regular on the axis and satisfies
// zero normal velocity and zero tangential traction on the entire box.
struct Exact {double radial,axial,rr,zz,hoop,shear,forceR,forceZ,mu;};
Exact exact(double r,double z,bool variable){
    const double r2=r*r,r4=r2*r2,r6=r4*r2;
    const double f=r*(1-3*r2+3*r4-r6),df=1-9*r2+15*r4-7*r6;
    const double g=2-12*r2+18*r4-8*r6,dg=r*(-24+72*r2-48*r4);
    const double c=std::cos(pi*z),s=std::sin(pi*z);
    Exact q{};q.radial=-pi*f*c;q.axial=g*s;q.rr=-pi*df*c;q.zz=pi*g*c;
    q.hoop=-pi*(1-3*r2+3*r4-r6)*c;q.shear=(pi*pi*f+dg)*s;
    q.mu=mu0*(1+(variable?r2+z:0));
    const double mur=variable?2*mu0*r:0,muz=variable?mu0:0;
    // Independent cylindrical div(2*mu*D), with analytic derivatives.
    q.forceR=-q.mu*pi*(dg-pi*pi*f)*c+2*mur*q.rr+muz*q.shear;
    q.forceZ=q.mu*(-48+288*r2-288*r4-pi*pi*g)*s+mur*q.shear+2*muz*q.zz;
    return q;
}
template<class F> double integrate(double r0,double r1,double z0,double z1,const F& value,int panels=1){
    long double sum=0;const auto& gauss=P::gaussEight();
    const double dr=(r1-r0)/panels,dz=(z1-z0)/panels;
    for(int j=0;j<panels;++j)for(int i=0;i<panels;++i)
        for(int b=0;b<8;++b)for(int a=0;a<8;++a){
            const double r=r0+(i+.5*(1+gauss.x[a]))*dr,z=z0+(j+.5*(1+gauss.x[b]))*dz;
            sum+=.25L*dr*dz*gauss.w[a]*gauss.w[b]*(2*pi*r)*value(r,z);}
    return static_cast<double>(sum);
}
struct Error {double force,energy,velocity;};
Error check(int n,bool variable){
    DictionaryMesh input;input.minCoord={0,0};input.maxCoord={1,1};input.Nx=input.Ny=n;
    input.maxCellNum=4096;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);FaceSegmentGraph graph;graph.build(mesh);
    const auto domain=NativeLeafDomain::build(mesh,graph,true);
    const auto pressure=HydroConservativePressure::buildDivergenceAdjoint(mesh,domain);
    const auto material=[&](double r0,double r1,double z0,double z1){
        const double r2=r1*r1-r0*r0,z=z1-z0;
        return mu0*2*pi*(.5*r2*z+(variable?.25*(r1*r1*r1*r1-r0*r0*r0*r0)*z+.25*r2*(z1*z1-z0*z0):0));};
    const auto k=HydroConservativeStress::assembleEnergy(mesh,domain,material);
    const Eigen::SparseMatrix<double> asymmetry=k-Eigen::SparseMatrix<double>(k.transpose());
    const auto zero=HydroConservativeStress::assembleEnergy(mesh,domain,[](double,double,double,double){return 0.;});
    Eigen::VectorXd u(domain.nativeFaceCount()),affine(u.size()),translation(u.size());
    for(const auto face:domain.faces(mesh)){
        const auto q=exact(face.center.x,face.center.y,variable);
        u[face.nativeRow]=face.outwardSign?0:(face.horizontal?q.axial:q.radial);
        affine[face.nativeRow]=face.horizontal?-.3*face.center.y:.7*face.center.x;
        translation[face.nativeRow]=face.horizontal?1:0;
    }
    // Nonhomogeneous affine traces are used only to check strain energy and
    // regularity, not as a solution of the closed-box boundary problem.
    const double affineReference=2*(2*.7*.7+.3*.3)*material(0,1,0,1);
    const double affineError=std::abs(affine.dot(k*affine)-affineReference);
    const double symmetry=asymmetry.norm()/k.norm(),translationError=(k*translation).norm();
    require(symmetry<1e-14 && zero.norm()==0 && affineError<1e-12 && translationError<1e-12,
        "RZ viscosity lost symmetry, affine hoop strain, translation nullspace or zero-mu limit");
    const Eigen::VectorXd force=-k*u;long double error2=0,reference2=0;double axisError=0,maxError=0;
    Vector2D worstPoint{};bool worstAxial=false;
    for(const auto face:domain.faces(mesh))if(!face.outwardSign){
        const auto& dual=pressure.duals().volumes()[face.nativeRow];
        const double reference=integrate(dual.lower.x,dual.upper.x,dual.lower.y,dual.upper.y,
            [&](double r,double z){const auto q=exact(r,z,variable);return face.horizontal?q.forceZ:q.forceR;});
        const double difference=force[face.nativeRow]-reference;
        error2+=difference*difference/dual.measure;reference2+=reference*reference/dual.measure;
        if(face.center.x<1.01/n)axisError=std::max(axisError,std::abs(difference)/dual.measure);
        if(std::abs(difference)/dual.measure>maxError){maxError=std::abs(difference)/dual.measure;
            worstPoint=face.center;worstAxial=face.horizontal;}
    }
    const auto density=[&](double r,double z){const auto q=exact(r,z,variable);
        return q.mu*(2*(q.rr*q.rr+q.zz*q.zz+q.hoop*q.hoop)+q.shear*q.shear);};
    const double reference=integrate(0,1,0,1,density,4),refined=integrate(0,1,0,1,density,8);
    require(std::abs(reference-refined)<1e-12*reference,"independent dissipation integral not resolved");
    const double energy=u.dot(k*u);require(energy>0,"nonzero smooth RZ flow lacks viscous dissipation");
    // Solve the manufactured vector Helmholtz problem u-div(2*mu*D(u))=b.
    // This isolates the viscous operator, including its natural slip boundary;
    // it is not an incompressible time step or a moving-interface validation.
    std::vector<int> freeRow(u.size(),-1);int freeCount=0;
    for(const auto face:domain.faces(mesh))if(!face.outwardSign)freeRow[face.nativeRow]=freeCount++;
    std::vector<Eigen::Triplet<double>> entries;
    for(int column=0;column<k.outerSize();++column)for(Eigen::SparseMatrix<double>::InnerIterator it(k,column);it;++it)
        if(freeRow[it.row()]>=0 && freeRow[it.col()]>=0)entries.emplace_back(freeRow[it.row()],freeRow[it.col()],it.value());
    Eigen::VectorXd rhs(freeCount);
    for(const auto face:domain.faces(mesh))if(!face.outwardSign){
        const auto& dual=pressure.duals().volumes()[face.nativeRow];const int row=freeRow[face.nativeRow];
        entries.emplace_back(row,row,dual.measure);
        rhs[row]=integrate(dual.lower.x,dual.upper.x,dual.lower.y,dual.upper.y,[&](double r,double z){
            const auto q=exact(r,z,variable);return face.horizontal?q.axial-q.forceZ:q.radial-q.forceR;});
    }
    Eigen::SparseMatrix<double> system(freeCount,freeCount);system.setFromTriplets(entries.begin(),entries.end());
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;solver.compute(system);
    require(solver.info()==Eigen::Success && solver.vectorD().minCoeff()>0,"viscous resolvent is not positive definite");
    const Eigen::VectorXd solved=solver.solve(rhs);const double residual=(system*solved-rhs).norm()/rhs.norm();
    require(solver.info()==Eigen::Success && residual<1e-11,"manufactured viscous solve failed");
    long double velocityError=0,velocityNorm=0;
    for(const auto face:domain.faces(mesh))if(!face.outwardSign){
        const double volume=pressure.duals().volumes()[face.nativeRow].measure,difference=solved[freeRow[face.nativeRow]]-u[face.nativeRow];
        velocityError+=volume*difference*difference;velocityNorm+=volume*u[face.nativeRow]*u[face.nativeRow];}
    const Error result{std::sqrt(static_cast<double>(error2/reference2)),std::abs(energy-reference)/reference,
        std::sqrt(static_cast<double>(velocityError/velocityNorm))};
    std::cout<<std::setprecision(17)<<"HYDRO_ENERGY_VISCOSITY grid="<<n<<" variable_mu="<<variable
        <<" force_relative_l2="<<result.force<<" energy_relative_error="<<result.energy
        <<" solved_velocity_relative_l2="<<result.velocity<<" solve_relative_residual="<<residual
        <<" axis_force_max_error="<<axisError<<" discrete_dissipation="<<energy<<" exact_dissipation="<<reference
        <<" symmetry_relative="<<symmetry<<" affine_energy_error="<<affineError<<" translation_error="<<translationError
        <<" force_max_error="<<maxError<<" worst_r="<<worstPoint.x<<" worst_z="<<worstPoint.y<<" worst_axial="<<worstAxial
        <<" zero_viscosity_norm="<<zero.norm()<<'\n'<<std::flush;
    return result;
}
struct CoupledError {double velocity,pressure;};
CoupledError checkCoupled(int n,bool variable){
    DictionaryMesh input;input.minCoord={0,0};input.maxCoord={1,1};input.Nx=input.Ny=n;
    input.maxCellNum=4096;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);FaceSegmentGraph graph;graph.build(mesh);
    const auto domain=NativeLeafDomain::build(mesh,graph,true);
    const auto pressure=HydroConservativePressure::buildDivergenceAdjoint(mesh,domain);
    const auto kind=HydroStepGeometry::complete(true);
    const double dt=.5;
    const HydroStepTime time{1,0,dt,HydroTimePolicy::ConservativeOneStep};
    // An actual managed CN Stokes start, with zero initial advecting velocity.
    // Manufactured final u is the analytic solenoidal field above; pressure
    // is nonconstant. The continuous force is integrated independently on
    // each physical dual, never obtained by applying the tested matrices.
    std::vector<double> advector(domain.nativeFaceCount(),0);
    const auto transport=HydroConservativeTransport::singlePhase(mesh,domain,time,1,advector);
    CellField<double> rho(mesh,1),mu(mesh,0);VertField<double> vrho(mesh,1),vmu(mesh,0);
    for(int c=0;c<mesh.cellIndexNum();++c){const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);
        const auto a=mesh.verts(cell.cellToVert[0]).coord(),b=mesh.verts(cell.cellToVert[2]).coord();
        mu.value[c]=mu0*(1+(variable?.5*(a.x*a.x+b.x*b.x+a.y+b.y):0));}
    for(int v=0;v<mesh.vertIndexNum();++v){const auto point=mesh.verts(mesh.vertIndexToLabel[v]).coord();
        vmu.value[v]=exact(point.x,point.y,variable).mu;}
    VEdgeField<double> forceR(mesh,0);HEdgeField<double> forceZ(mesh,0);
    for(const auto face:domain.faces(mesh))if(!face.outwardSign){
        const auto& dual=pressure.duals().volumes()[face.nativeRow];
        const double load=integrate(dual.lower.x,dual.upper.x,dual.lower.y,dual.upper.y,[&](double r,double z){
            const auto q=exact(r,z,variable);
            return face.horizontal?q.axial/dt-.5*q.forceZ-pi*(r*r-.5)*std::sin(pi*z)
                :q.radial/dt-.5*q.forceR+2*r*std::cos(pi*z);})/dual.measure;
        (face.horizontal?forceZ.value[face.nativeRow-mesh.vEdgeIndexNum()]:forceR.value[face.nativeRow])=load;
    }
    Hydrodynamics hydro(mesh);hydro.initializeState(mesh,kind);
    hydro.updateMaterials(mesh,{kind,0,HydroMaterialEvolution::Frozen,HydroNativeMaterial{&rho,&mu,&vrho,&vmu}});
    HydroSolveOptions options;options.compatibleNativePressure=true;options.linearResidualTolerance=1e-10;
    const HydroForceDensity force{nullptr,nullptr,&forceR,&forceZ,0,hydro.captureForceSource(0,HydroForceEvolution::Frozen)};
    const auto report=hydro.advance(mesh,{kind,time,options,force,nullptr,&transport});
    if(!report.success)throw std::runtime_error("RZ coupled manufactured start: "+report.stage+": "+report.reason);
    require(report.physicalNativeMomentumUsed && report.compatibleNativePressureUsed && report.coupledUsedSuperLU,
        "RZ coupled manufactured flow bypassed the managed production solver");
    const auto state=hydro.state();long double velocityError=0,velocityNorm=0,pError=0,pNorm=0,pMean=0,volume=0;
    std::vector<double> solved(domain.nativeFaceCount());
    for(const auto face:domain.faces(mesh)){
        solved[face.nativeRow]=face.horizontal?state.velocityY[face.nativeRow-mesh.vEdgeIndexNum()]:state.velocityX[face.nativeRow];
        const auto q=exact(face.center.x,face.center.y,variable);const double target=face.outwardSign?0:(face.horizontal?q.axial:q.radial);
        const double weight=pressure.duals().volumes()[face.nativeRow].measure,difference=solved[face.nativeRow]-target;
        velocityError+=weight*difference*difference;velocityNorm+=weight*target*target;
    }
    for(int c=0;c<mesh.cellIndexNum();++c){const auto point=mesh.cells(mesh.cellIndexToLabel[c]).centerCoord();
        const double target=(point.x*point.x-.5)*std::cos(pi*point.y),weight=domain.cellVolumes()[domain.cellRowByIndex()[c]],
            difference=state.physicalPressure[c]-target;
        pError+=weight*difference*difference;pNorm+=weight*target*target;pMean+=weight*state.physicalPressure[c];volume+=weight;}
    const auto k=HydroConservativeStress::assembleCellEnergy(mesh,domain,std::vector<double>(mu.value.get(),mu.value.get()+mu.size()));
    const Eigen::Map<const Eigen::VectorXd> vector(solved.data(),static_cast<int>(solved.size()));
    const double workIdentity=std::abs(report.physicalViscousPower+.25*vector.dot(k*vector));
    require(report.divergenceAfterProjection<2e-11 && std::abs(report.physicalPressurePower)<1e-12
        && report.physicalViscousPower<0 && workIdentity<1e-12 && std::abs(pMean/volume)<1e-13,
        "RZ coupled step lost incompressibility, pressure gauge/work or viscous energy");
    // Its nonzero solved velocity must now be usable as the next constant-
    // density mass source, without manufacturing an Interface inventory.
    const auto following=HydroConservativeTransport::singlePhase(mesh,domain,
        {2,dt,.001,HydroTimePolicy::ConservativeOneStep},1,solved);
    const auto uniformMomentum=following.upwind(std::vector<double>(solved.size(),1));
    double constantDefect=0;for(auto value:uniformMomentum.outward)constantDefect=std::max(constantDefect,std::abs(static_cast<double>(value)));
    require(constantDefect<1e-14,"single-phase dual mass flux does not preserve uniform transported velocity");
    if(n==8 && !variable){auto divergent=solved;
        for(const auto face:domain.faces(mesh))if(!face.outwardSign){divergent[face.nativeRow]+=.01;break;}
        bool rejected=false;try{(void)HydroConservativeTransport::singlePhase(mesh,domain,time,1,divergent);}
        catch(const std::invalid_argument&){rejected=true;}
        require(rejected,"single-phase transport accepted a non-solenoidal mass source");}
    if(n==32 && variable)if(const char* path=std::getenv("QUANTUM_RZ_COUPLED_FIELDS")){
        std::ofstream output(path);require(bool(output),"cannot write RZ coupled fields");
        output<<std::setprecision(17)<<"r,z,ur,uz,p,ur_exact,uz_exact,p_exact\n";
        for(int c=0;c<mesh.cellIndexNum();++c){const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto point=cell.centerCoord();
            const auto q=exact(point.x,point.y,variable);
            const double ur=.5*(state.velocityX[mesh.vEdges(cell.cellToVEdge[0][0]).index()]+state.velocityX[mesh.vEdges(cell.cellToVEdge[1][0]).index()]),
                uz=.5*(state.velocityY[mesh.hEdges(cell.cellToHEdge[0][0]).index()]+state.velocityY[mesh.hEdges(cell.cellToHEdge[1][0]).index()]);
            output<<point.x<<','<<point.y<<','<<ur<<','<<uz<<','<<state.physicalPressure[c]<<','<<q.radial<<','<<q.axial<<','
                <<(point.x*point.x-.5)*std::cos(pi*point.y)<<'\n';}
    }
    const CoupledError result{std::sqrt(static_cast<double>(velocityError/velocityNorm)),std::sqrt(static_cast<double>(pError/pNorm))};
    std::cout<<std::setprecision(17)<<"HYDRO_ENERGY_RZ_COUPLED grid="<<n<<" variable_mu="<<variable
        <<" velocity_relative_l2="<<result.velocity<<" pressure_relative_l2="<<result.pressure
        <<" divergence="<<report.divergenceAfterProjection<<" pressure_power="<<report.physicalPressurePower
        <<" viscous_power="<<report.physicalViscousPower<<" viscous_work_identity="<<workIdentity
        <<" momentum_residual="<<report.momentumBalanceRelativeResidual<<" single_phase_mass_closure="<<following.maxRelativeClosure()
        <<" uniform_momentum_defect="<<constantDefect<<" managed_steps="<<state.progress.acceptedSteps
        <<" time_order_claim=0 nonlinear_advection_claim=0\n"<<std::flush;
    return result;
}
}
void hydroEnergyViscosityContract(){
    bool converged=true;
    for(bool variable:{false,true}){
        auto previous=check(8,variable);
        for(int n:{16,32}){const auto current=check(n,variable);
            converged=converged && current.force<previous.force && previous.velocity/current.velocity>3.3
                && previous.energy/current.energy>3.3;
            previous=current;}
    }
    require(converged,"RZ manufactured velocity/energy failed second-order refinement");
    bool coupledConverged=true;
    for(bool variable:{false,true}){
        auto previous=checkCoupled(8,variable);
        for(int n:{16,32}){const auto current=checkCoupled(n,variable);
            coupledConverged=coupledConverged && previous.velocity/current.velocity>3.3 && previous.pressure/current.pressure>3.3;
            previous=current;}
    }
    require(coupledConverged,"RZ coupled velocity/pressure failed second-order refinement");
}
