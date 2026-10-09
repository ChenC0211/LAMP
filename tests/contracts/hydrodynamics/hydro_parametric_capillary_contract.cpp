#include "Modules/Hydrodynamics/Hydrodynamics.h"
#include "Modules/Hydrodynamics/HydroConservativePressure.h"
#include "Modules/Hydrodynamics/HydroConservativeStress.h"
#include "Modules/Hydrodynamics/HydroParametricStep.h"
#include "Modules/Hydrodynamics/HydroConservativeTransport.h"
#include "Modules/Hydrodynamics/HydroUniformLoad.h"
#include "Modules/Hydrodynamics/HydroSplineVelocity.h"
#include "Modules/Hydrodynamics/HydroInterfaceTransport.h"
#include "Modules/Interface/Compatibility/InterfaceParametricMotion.h"
#include "Modules/Interface/Compatibility/InterfacePrescribedData.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <fstream>
#include <cstdlib>

namespace {
using namespace Quantum::Hydrodynamics;
namespace P=Quantum::Geometry::Parametric;
void require(bool ok,const char* reason){if(!ok)throw std::runtime_error(reason);}
DictionaryMesh meshInput(int n){DictionaryMesh input;input.minCoord={0,0};input.maxCoord={1,1};
    input.Nx=input.Ny=n;input.maxCellNum=4096;input.maxInitialIterations=0;return input;}
std::array<double,2> diameters(const P::Curve& curve){std::array<double,2> result{};
    for(int axis=0;axis<2;++axis){double lo=1,hi=0;for(const auto& arc:curve.arcs){const auto bounds=P::range(arc.coordinate(axis));
        lo=std::min(lo,bounds[0]);hi=std::max(hi,bounds[1]);}result[axis]=hi-lo;}return result;}
HydroForceDensity forceInput(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,Hydrodynamics& hydro,
    const InterfaceState& state,const HydroUniformLoad& load,HydroForceTimeLevel level){
    HydroForceDensity force;force.faceX=&load.faceForceX;force.faceY=&load.faceForceY;
    force.evaluationTime=state.progress().endTime;force.timeLevel=level;
    force.source=hydro.captureParametricForceSource(mesh,domain,state);return force;
}
double maxVelocity(const Hydrodynamics& hydro){const auto state=hydro.state();double result=0;
    for(int row=0;row<state.velocityX.size();++row)result=std::max(result,std::abs(state.velocityX[row]));
    for(int row=0;row<state.velocityY.size();++row)result=std::max(result,std::abs(state.velocityY[row]));return result;}
void workPair(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const InterfaceState& state,
    const HydroConservativePressure& pressure,const HydroUniformLoad& load){
    VEdgeField<double> ux(mesh,0);HEdgeField<double> uy(mesh,0);const double pi=std::acos(-1.);
    const auto stream=[&](Vector2D p){return p.x==0 || p.x==1 || p.y==0 || p.y==1?0:
        .07*std::sin(2*pi*p.x)*std::sin(2*pi*p.y)/(2*pi)+.02*std::sin(pi*p.x)*std::sin(2*pi*p.y)/pi;};
    for(const auto face:domain.faces(mesh)){const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
        const auto lo=mesh.verts(edge.edgeToVert[0]).coord(),hi=mesh.verts(edge.edgeToVert[1]).coord();
        (face.horizontal?uy.value:ux.value)[edge.index()]=face.outwardSign?0:(face.horizontal?-1:1)*(stream(hi)-stream(lo))/face.area;}
    const HydroSplineVelocity velocity(mesh,domain,ux,uy);long double nativePower=0,lengthRate=0,normalPower=0,areaRate=0;
    for(const auto face:domain.faces(mesh)){const int row=face.nativeRow;const double u=face.horizontal?uy.value[row-mesh.vEdgeIndexNum()]:ux.value[row];
        const double f=face.horizontal?load.faceForceY.value[row-mesh.vEdgeIndexNum()]:load.faceForceX.value[row];
        nativePower+=static_cast<long double>(u)*f*pressure.duals().volumes()[row].measure;}
    const int n=state.parametricGeometry().grid().nx;
    for(const auto& arc:state.parametricGeometry().geometry().arcs){std::vector<double> cuts{0,1};
        for(int axis=0;axis<2;++axis)for(int line=1;line<n;++line){auto equation=arc.coordinate(axis);equation[0]-=double(line)/n;
            for(double t:P::roots(equation,0,1))cuts.push_back(t);}
        std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
        for(std::size_t part=1;part<cuts.size();++part)for(int sub=0;sub<2;++sub)for(int k=0;k<8;++k){
            const double width=(cuts[part]-cuts[part-1])*.5,lo=cuts[part-1]+sub*width;
            const auto q=P::sampleSurface(arc,lo+.5*width*(1+P::gaussEight().x[k]));const auto u=velocity.firstJet({q.point.x,q.point.y});
            const auto tangent=q.derivative/q.speed;const double weight=.5*width*P::gaussEight().w[k]*q.speed;
            const double normal=u.value.x*q.normal.x+u.value.y*q.normal.y;
            areaRate+=weight*normal;normalPower+=weight*q.curvature*normal;
            // Independent first variation of length: t.(grad u).t, using the
            // forward velocity's actual jet rather than the load's curvature.
            lengthRate+=weight*(tangent.x*(u.dx.x*tangent.x+u.dy.x*tangent.y)
                +tangent.y*(u.dx.y*tangent.x+u.dy.y*tangent.y));
        }
    }
    const double error=static_cast<double>(std::abs(nativePower+lengthRate));
    require(error<1e-10 && std::abs(lengthRate-normalPower)<1e-10 && std::abs(areaRate)<1e-12 && std::abs(lengthRate)>1e-5,
        "capillary work is not the negative first variation of length under the actual C1 velocity");
    // General (not divergence-free) native data exercise the reconstruction
    // transpose independently of the capillary pressure-mode identity.
    for(const auto face:domain.faces(mesh))(face.horizontal?uy.value:ux.value)[face.horizontal?face.nativeRow-mesh.vEdgeIndexNum():face.nativeRow]
        =face.outwardSign?0:.1*std::sin(.37*face.nativeRow);
    const HydroSplineVelocity arbitrary(mesh,domain,ux,uy);HydroSplineVelocityAdjoint transpose(mesh,domain);long double pointPower=0;
    for(int cell:domain.cellIndexByRow()){const auto& c=mesh.cells(mesh.cellIndexToLabel[cell]);
        const auto lower=mesh.verts(c.cellToVert[0]).coord(),upper=mesh.verts(c.cellToVert[2]).coord();
        const Vector2D local{.23,.61},point{lower.x+local.x*(upper.x-lower.x),lower.y+local.y*(upper.y-lower.y)},traction{std::sin(.31*cell),std::cos(.27*cell)};
        transpose.accumulate(cell,transpose.moment(cell,local,traction));const auto value=arbitrary.firstJet(point).value;
        pointPower+=traction.x*value.x+traction.y*value.y;}
    const auto native=transpose.integratedForce();long double transposePower=0;
    for(const auto face:domain.faces(mesh)){const double u=face.horizontal?uy.value[face.nativeRow-mesh.vEdgeIndexNum()]:ux.value[face.nativeRow];
        transposePower+=static_cast<long double>(u)*native[face.nativeRow];}
    const double transposeError=static_cast<double>(std::abs(pointPower-transposePower)/std::max(1.L,std::abs(pointPower)));
    require(transposeError<5e-13,"native load is not the transpose of the existing spline reconstruction");
    std::cout<<std::setprecision(17)<<"HYDRO_PARAMETRIC_CAPILLARY_WORK native_power="<<static_cast<double>(nativePower)
        <<" surface_length_rate="<<static_cast<double>(lengthRate)<<" curvature_normal_power="<<static_cast<double>(normalPower)
        <<" area_rate="<<static_cast<double>(areaRate)<<" virtual_work_error="<<error<<" arbitrary_transpose_relative_error="<<transposeError
        <<" full_discrete_energy_claim=0\n"<<std::flush;
}
struct StaticResult {double curvatureError=0,velocity=0;};
StaticResult circle(int arcs){
    constexpr int n=16;constexpr double radius=.19,dt=.001;
    auto input=meshInput(n);QuadtreeMesh mesh(input);FaceSegmentGraph graph;graph.build(mesh);
    const auto domain=NativeLeafDomain::build(mesh,graph,false);const auto pressure=HydroConservativePressure::buildDivergenceAdjoint(mesh,domain);
    auto state=InterfaceParametricMotion::initialize(mesh,domain,InterfacePrescribedData::parametricEllipse({.47,.53},{radius,radius},arcs));
    const auto kind=HydroStepGeometry::complete(false);Hydrodynamics hydro(mesh);hydro.initializeState(mesh,kind,{});
    hydro.updateMaterials(mesh,{kind,0,HydroMaterialEvolution::TimeDependent,HydroParametricMaterial{&domain,&state,1,.1,.01,.001}});
    HydroSolveOptions options;options.compatibleNativePressure=true;InterfaceParametricProjection::Options geometry;geometry.materialSegments=arcs;
    const auto velocity=hydro.splineVelocity(mesh,domain,HydroVelocityTimeLevel::Current);
    const auto flux=HydroInterfaceTransport::normalFlux(velocity);const auto flow=HydroInterfaceTransport::frozenParametricFlow(velocity);
    const InterfaceStepTime time{1,0,dt,dt};auto next=InterfaceParametricMotion::prepareTimeLinear(mesh,domain,state,time,flow,flux,flux,geometry);
    const auto mass=HydroConservativeTransport::build(mesh,domain,state,next.state,time,1,.1);
    ParametricCapillaryReport capillary;const auto load=buildParametricCapillaryLoad(mesh,domain,pressure,state,1,&capillary);
    const auto force=forceInput(mesh,domain,hydro,state,load,HydroForceTimeLevel::StepBegin);
    // Rejected candidate/source/time combinations must leave Hydro fields and
    // accepted physical time intact. This uses the same real mass step.
    int rejected=0;
    auto unrelated=InterfaceParametricMotion::initialize(mesh,domain,InterfacePrescribedData::parametricEllipse({.47,.53},{radius,radius},arcs));
    for(int wrong=0;wrong<3;++wrong){auto bad=force;const InterfaceState* context=&state;
        if(wrong==0)context=&unrelated;
        if(wrong==1)bad.source=hydro.captureParametricForceSource(mesh,domain,unrelated);
        if(wrong==2){bad.timeLevel=HydroForceTimeLevel::StepEnd;bad.evaluationTime=dt;context=&next.state;}
        auto transaction=hydro.beginStepTransaction(mesh);
        const auto report=hydro.advance(mesh,{kind,{1,0,dt,HydroTimePolicy::ConservativeOneStep},options,bad,nullptr,&mass,context});
        require(!report.success && report.stage=="step-input","parameter force accepted inconsistent geometry/time provenance");
        transaction.rollback();require(hydro.state().progress.completedSteps==0 && maxVelocity(hydro)==0,"rejected force changed accepted Hydro state");++rejected;
    }
    const auto solved=hydro.advance(mesh,{kind,{1,0,dt,HydroTimePolicy::ConservativeOneStep},options,force,nullptr,&mass,&state});
    if(!solved.success)throw std::runtime_error("parameter circle solve: "+solved.stage+": "+solved.reason);
    // Analytic radius is used ONLY by this independent error diagnostic. The
    // production load gets numerical curve coefficients, inventory and sigma.
    double curvatureIntegral=0;
    for(const auto& arc:state.parametricGeometry().geometry().arcs)for(int k=0;k<8;++k){const auto q=P::sampleSurface(arc,.5*(1+P::gaussEight().x[k]));
        curvatureIntegral+=.5*P::gaussEight().w[k]*q.speed*std::pow(radius*q.curvature-1,2);}
    StaticResult result{std::sqrt(curvatureIntegral/capillary.length),maxVelocity(hydro)};
    const double volume=std::accumulate(state.parametricGeometry().inventory().begin(),state.parametricGeometry().inventory().end(),0.);
    const auto& cells=InterfaceParametricMotion::cellIndices(state);double pressureError=0;
    for(std::size_t i=0;i<cells.size();++i){const int cell=cells[i];const double alpha=state.parametricGeometry().inventory()[i]*n*n;
        pressureError=std::max(pressureError,std::abs(hydro.state().physicalPressure[cell]-capillary.meanCurvature*(alpha-volume)));}
    require(solved.compatibleNativePressureUsed && solved.conservativeMomentumUsed && solved.coupledUsedSuperLU,
        "numerical curvature bypassed the actual compatible Hydro solve");
    require(std::abs(load.totalForceX)<1e-10 && std::abs(load.totalForceY)<1e-10,
        "parameter force lost the closed-curve resultant including wall reactions");
    std::cout<<std::setprecision(17)<<"HYDRO_PARAMETRIC_CAPILLARY_CIRCLE grid="<<n<<" arcs="<<arcs<<" dt="<<dt
        <<" curvature_relative_l2="<<result.curvatureError<<" mean_curvature="<<capillary.meanCurvature<<" max_native_velocity="<<result.velocity
        <<" pressure_constant_mode_error="<<pressureError<<" total_force_x="<<load.totalForceX<<" total_force_y="<<load.totalForceY
        <<" quadrature_error="<<capillary.quadratureError
        <<" rejected_provenance="<<rejected<<" analytic_curvature_input=0\n"<<std::flush;
    return result;
}
void movingEllipse(int steps){
    constexpr int n=16;constexpr double sigma=1;const double dt=.024/steps;
    auto input=meshInput(n);QuadtreeMesh mesh(input);FaceSegmentGraph graph;graph.build(mesh);
    const auto domain=NativeLeafDomain::build(mesh,graph,false);const auto pressure=HydroConservativePressure::buildDivergenceAdjoint(mesh,domain);
    const auto initial=InterfacePrescribedData::parametricEllipse({.5,.5},{.23,.17},32);
    auto state=InterfaceParametricMotion::initialize(mesh,domain,initial);
    const double initialVolume=std::accumulate(state.parametricGeometry().inventory().begin(),state.parametricGeometry().inventory().end(),0.);
    const auto initialDiameter=diameters(initial);ParametricCapillaryReport initialSurface;
    const auto initialLoad=buildParametricCapillaryLoad(mesh,domain,pressure,state,sigma,&initialSurface);
    workPair(mesh,domain,state,pressure,initialLoad);
    const auto kind=HydroStepGeometry::complete(false);Hydrodynamics hydro(mesh);hydro.initializeState(mesh,kind,{});
    HydroSolveOptions options;options.compatibleNativePressure=true;InterfaceParametricProjection::Options geometry;geometry.materialSegments=32;
    double maxFraction=0,maxVolume=0,maxMass=0,maxMomentum=0,maxForceResultant=0,maxEnergyIncrease=0;
    double kinetic=0,peakVelocity=0;ParametricCapillaryReport surface;
    for(int step=1;step<=steps;++step){const InterfaceStepTime time{step,(step-1)*dt,step*dt,dt};
        auto transaction=hydro.beginStepTransaction(mesh);
        hydro.updateMaterials(mesh,{kind,time.begin,HydroMaterialEvolution::TimeDependent,HydroParametricMaterial{&domain,&state,1,.1,.01,.001}});
        const auto velocity=hydro.splineVelocity(mesh,domain,HydroVelocityTimeLevel::Current);
        const auto flux=HydroInterfaceTransport::normalFlux(velocity);const auto flow=HydroInterfaceTransport::frozenParametricFlow(velocity);
        auto next=InterfaceParametricMotion::prepareTimeLinear(mesh,domain,state,time,flow,flux,flux,geometry);
        const auto mass=HydroConservativeTransport::build(mesh,domain,state,next.state,time,1,.1);
        const auto load=buildParametricCapillaryLoad(mesh,domain,pressure,next.state,sigma,&surface);
        const auto force=forceInput(mesh,domain,hydro,next.state,load,HydroForceTimeLevel::StepEnd);
        auto candidate=hydro.prepareStep(mesh,{kind,{step,time.begin,dt,HydroTimePolicy::ConservativeOneStep},options,force,nullptr,&mass,&next.state});
        const auto& solved=candidate.report();
        if(!solved.success)throw std::runtime_error("parameter capillary motion: "+solved.stage+": "+solved.reason);
        const auto fields=hydro.state();kinetic=0;
        for(int row=0;row<domain.nativeFaceCount();++row){const double u=row<mesh.vEdgeIndexNum()?fields.velocityX[row]:fields.velocityY[row-mesh.vEdgeIndexNum()];
            kinetic+=static_cast<double>(.5L*mass.massAfter()[row]*u*u);}
        const double total=std::accumulate(next.state.parametricGeometry().inventory().begin(),next.state.parametricGeometry().inventory().end(),0.);
        maxVolume=std::max(maxVolume,std::abs(total-initialVolume));maxMass=std::max(maxMass,mass.maxRelativeClosure());
        maxFraction=std::max(maxFraction,next.state.parametricGeometry().report().maximumInventoryFractionError);
        maxMomentum=std::max(maxMomentum,solved.momentumBalanceRelativeResidual);
        maxForceResultant=std::max(maxForceResultant,std::hypot(load.totalForceX,load.totalForceY));
        maxEnergyIncrease=std::max(maxEnergyIncrease,(kinetic+sigma*surface.length)/(sigma*initialSurface.length)-1);
        peakVelocity=std::max(peakVelocity,maxVelocity(hydro));
        candidate.accept(mesh);state=std::move(next.state);transaction.commit();
        const auto size=diameters(state.parametricGeometry().geometry());
        std::cout<<std::setprecision(17)<<"HYDRO_PARAMETRIC_CAPILLARY_STEP step="<<step<<" time="<<time.end
            <<" diameter_x="<<size[0]<<" diameter_y="<<size[1]<<" length="<<surface.length<<" kinetic_energy="<<kinetic
            <<" total_energy="<<kinetic+sigma*surface.length<<" max_velocity="<<maxVelocity(hydro)<<'\n'<<std::flush;
    }
    const auto finalDiameter=diameters(state.parametricGeometry().geometry());
    require(finalDiameter[0]<initialDiameter[0]-1e-5 && finalDiameter[1]>initialDiameter[1]+1e-5 && kinetic>1e-7,
        "surface tension did not drive the initially motionless ellipse toward a circle");
    require(maxVolume<1e-14 && maxFraction<3e-12 && state.progress().completedSteps==steps && hydro.state().progress.acceptedSteps==steps,
        "capillary motion did not preserve one accepted geometry/inventory/fluid history");
    std::cout<<std::setprecision(17)<<"HYDRO_PARAMETRIC_CAPILLARY_DYNAMIC grid="<<n<<" steps="<<steps<<" end_time="<<steps*dt
        <<" initial_diameter_x="<<initialDiameter[0]<<" final_diameter_x="<<finalDiameter[0]
        <<" initial_diameter_y="<<initialDiameter[1]<<" final_diameter_y="<<finalDiameter[1]
        <<" initial_length="<<initialSurface.length<<" final_length="<<surface.length<<" final_kinetic="<<kinetic
        <<" peak_velocity="<<peakVelocity<<" maximum_energy_increase_relative="<<maxEnergyIncrease
        <<" inventory_drift="<<maxVolume<<" max_fraction_error="<<maxFraction<<" dual_mass_relative_closure="<<maxMass
        <<" momentum_relative_residual="<<maxMomentum<<" max_force_resultant="<<maxForceResultant
        <<" capillary_feedback=1 analytic_curvature_input=0 midpoint_claim=0 full_oscillation_claim=0\n";
}
}
void hydroParametricCapillaryContract(int steps){
    if(steps>0){movingEllipse(steps);return;}
    const auto coarse=circle(16),medium=circle(32),fine=circle(64);
    require(medium.curvatureError<coarse.curvatureError && fine.curvatureError<medium.curvatureError
        && medium.velocity<coarse.velocity && fine.velocity<medium.velocity,
        "numerical circle curvature refinement did not reduce curvature and induced velocity together");
}

namespace {
double shapeMoment(const P::Curve& curve){long double xx=0,yy=0;
    for(const auto& arc:curve.arcs)for(int k=0;k<8;++k){const auto q=P::sampleSurface(arc,.5*(1+P::gaussEight().x[k]));
        xx+=P::gaussEight().w[k]*std::pow(q.point.x-.5,3)*q.derivative.y/6;
        yy-=P::gaussEight().w[k]*std::pow(q.point.y-.5,3)*q.derivative.x/6;}
    return static_cast<double>((xx-yy)/(xx+yy));
}
void midpointOperators(){
    auto input=meshInput(8);QuadtreeMesh mesh(input);FaceSegmentGraph graph;graph.build(mesh);
    const auto domain=NativeLeafDomain::build(mesh,graph,false);
    const auto initial=InterfacePrescribedData::parametricEllipse({.5,.5},{.23,.17},32);
    auto old=InterfaceParametricMotion::initialize(mesh,domain,initial);
    HydroConservativeStress::EnergyReport viscousReport;
    const auto k=HydroConservativeStress::assembleParametricEnergy(mesh,domain,old,.03,.004,&viscousReport);
    const Eigen::SparseMatrix<double> asymmetry=k-Eigen::SparseMatrix<double>(k.transpose());
    const auto zero=HydroConservativeStress::assembleParametricEnergy(mesh,domain,old,0,0);
    Eigen::VectorXd affine(domain.nativeFaceCount()),random(affine.size());
    for(const auto face:domain.faces(mesh)){affine[face.nativeRow]=face.horizontal?face.center.y:face.center.x;random[face.nativeRow]=std::sin(.23*face.nativeRow);}
    const double expected=4*(.03*initial.area()+.004*(1-initial.area()));
    const double affineError=std::abs(affine.dot(k*affine)-expected),power=random.dot(k*random);
    require(asymmetry.norm()<1e-13*k.norm() && zero.norm()==0 && affineError<1e-12 && power>0,
        "energy viscosity lost symmetry, positive work, zero-viscosity limit or phase-integrated affine strain");
    VEdgeField<double> ux(mesh,0);HEdgeField<double> uy(mesh,0);const double pi=std::acos(-1.);
    const auto psi=[&](Vector2D p){return p.x==0 || p.x==1 || p.y==0 || p.y==1?0:.08*std::sin(pi*p.x)*std::sin(pi*p.y)/pi;};
    for(const auto face:domain.faces(mesh)){const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
        const auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
        (face.horizontal?uy.value:ux.value)[edge.index()]=face.outwardSign?0:(face.horizontal?-1:1)*(psi(b)-psi(a))/face.area;}
    const HydroSplineVelocity velocity(mesh,domain,ux,uy);const auto flux=HydroInterfaceTransport::normalFlux(velocity);
    const auto flow=HydroInterfaceTransport::frozenParametricFlow(velocity);InterfaceParametricProjection::Options geometry;geometry.materialSegments=32;
    const InterfaceStepTime time{1,0,.004,.004},half{1,0,.002,.002};
    auto middle=InterfaceParametricMotion::prepareTimeLinear(mesh,domain,old,half,flow,flux,flux,geometry);
    auto next=InterfaceParametricMotion::prepareTimeLinear(mesh,domain,old,time,flow,flux,flux,geometry);
    const auto mass=HydroConservativeTransport::build(mesh,domain,old,next.state,time,1000,1,&middle.state);
    std::vector<double> u(domain.nativeFaceCount());for(std::size_t row=0;row<u.size();++row)u[row]=std::sin(.37*row);
    const auto momentum=mass.centered(u);long double work=0,scale=0;double constantError=0;
    const auto constant=mass.centered(std::vector<double>(u.size(),.23));
    for(std::size_t row=0;row<u.size();++row){const auto dm=mass.massAfter()[row]-mass.massBefore()[row];
        work+=u[row]*momentum.outward[row]+.5L*dm*u[row]*u[row];scale+=std::abs(u[row]*momentum.outward[row]);
        constantError=std::max(constantError,static_cast<double>(std::abs(constant.outward[row]+.23L*dm)));}
    const double workRelative=static_cast<double>(std::abs(work)/std::max(1e-30L,scale));
    bool wrongMiddle=false;try{HydroConservativeTransport::build(mesh,domain,old,next.state,time,1000,1,&next.state);}catch(const std::exception&){wrongMiddle=true;}
    require(workRelative<1e-11 && constantError<1e-11 && wrongMiddle,"centered flux violates mass/kinetic pairing or midpoint lineage");
    const auto kind=HydroStepGeometry::complete(false);Hydrodynamics hydro(mesh);hydro.initializeState(mesh,kind,{});
    HydroParametricStepOptions options;options.liquidDensity=1;options.gasDensity=.1;options.liquidViscosity=.01;options.gasViscosity=.001;
    options.surfaceTension=1;options.geometry=geometry;options.maximumIterations=1;
    bool rollback=false;try{advanceParametricCapillaryStep(mesh,domain,hydro,old,half,options);}
    catch(const std::runtime_error& e){rollback=std::string(e.what()).find("did not converge")!=std::string::npos;}
    require(rollback && hydro.state().progress.completedSteps==0 && !hydro.state().progress.stopped && maxVelocity(hydro)==0,
        "unconverged midpoint iterate changed the accepted state");
    options.maximumIterations=16;const auto advanced=advanceParametricCapillaryStep(mesh,domain,hydro,old,half,options);
    require(advanced.hydro.midpointNativeMomentumUsed && advanced.hydro.physicalViscousPower<=0 && advanced.iterations>1,
        "midpoint control did not converge the actual coupled native equation");
    std::cout<<std::setprecision(17)<<"HYDRO_PARAMETRIC_MIDPOINT_OPERATORS viscosity_symmetry="<<asymmetry.norm()
        <<" affine_energy_error="<<affineError<<" positive_viscous_energy="<<power<<" zero_viscosity_norm="<<zero.norm()
        <<" measure_fraction_adjustment="<<viscousReport.maximumMeasureFractionAdjustment<<" centered_work_relative="<<workRelative
        <<" constant_momentum_error="<<constantError<<" wrong_midpoint_rejected="<<wrongMiddle<<" failed_iterate_rolled_back="<<rollback
        <<" coupled_iterations="<<advanced.iterations<<" coupled_velocity_defect="<<advanced.midpointVelocityDefect<<'\n'<<std::flush;
}
void midpointOscillation(int steps){
    constexpr int n=16;constexpr double endTime=.4;const double dt=endTime/steps;
    auto input=meshInput(n);QuadtreeMesh mesh(input);FaceSegmentGraph graph;graph.build(mesh);
    const auto domain=NativeLeafDomain::build(mesh,graph,false);const auto pressure=HydroConservativePressure::buildDivergenceAdjoint(mesh,domain);
    const auto initial=InterfacePrescribedData::parametricEllipse({.5,.5},{.23,.17},32);
    auto state=InterfaceParametricMotion::initialize(mesh,domain,initial);const auto kind=HydroStepGeometry::complete(false);
    Hydrodynamics hydro(mesh);hydro.initializeState(mesh,kind,{});HydroParametricStepOptions options;
    options.liquidDensity=1;options.gasDensity=.1;options.liquidViscosity=.01;options.gasViscosity=.001;options.surfaceTension=1;options.geometry.materialSegments=32;
    ParametricCapillaryReport surface;buildParametricCapillaryLoad(mesh,domain,pressure,state,1,&surface);
    const double initialEnergy=surface.length,initialVolume=std::accumulate(state.parametricGeometry().inventory().begin(),state.parametricGeometry().inventory().end(),0.);
    std::ofstream curveOutput;
    if(const char* path=std::getenv("QUANTUM_MIDPOINT_CURVES")){
        curveOutput.open(path);require(bool(curveOutput),"cannot open midpoint curve output");
        curveOutput<<std::setprecision(17)<<"step,time,arc,local_parameter,x,y,curvature\n";}
    const auto writeCurve=[&](int step,double time){
        if(!curveOutput.is_open())return;int index=0;
        for(const auto& arc:state.parametricGeometry().geometry().arcs){
            for(int k=0;k<=8;++k){const double s=double(k)/8;const auto sample=P::sampleSurface(arc,s);
                curveOutput<<step<<','<<time<<','<<index<<','<<s<<','<<sample.point.x<<','<<sample.point.y<<','<<sample.curvature<<'\n';}
            ++index;}
        curveOutput.flush();require(bool(curveOutput),"cannot write accepted midpoint curve");};
    writeCurve(0,0);
    double maximumIncrease=0,maxFraction=0,maxVolume=0,maxMass=0,maxDefect=0,maxDissipation=0,dissipation=0,temporalMass=0;
    double oldShape=shapeMoment(initial),previousSlope=0,finalEnergy=0;int signCrossings=0,turningPoints=0,maxIterations=0,totalIterations=0;
    for(int step=1;step<=steps;++step){const InterfaceStepTime time{step,(step-1)*dt,step*dt,dt};
        HydroParametricStepObserver observer;
        if(const char* archive=std::getenv("QUANTUM_MIDPOINT_CAPTURE"))observer=[&,archive](int iteration,const VEdgeField<double>& x,const HEdgeField<double>& y){
            std::ofstream out(archive);out<<std::setprecision(17)<<"QUANTUM_MIDPOINT_FLUX_INPUT_V1\n"
                <<n<<' '<<time.index<<' '<<iteration<<' '<<time.begin<<' '<<time.dt<<'\n';
            const auto& geometry=state.parametricGeometry();out<<geometry.geometry().arcs.size()<<'\n';
            for(const auto& arc:geometry.geometry().arcs){out<<arc.a<<' '<<arc.b;for(const auto& p:arc.p)out<<' '<<p.x<<' '<<p.y;out<<'\n';}
            out<<geometry.inventory().size()<<'\n';for(double value:geometry.inventory())out<<value<<'\n';
            out<<x.size()<<' '<<y.size()<<'\n';for(int row=0;row<x.size();++row)out<<x.value[row]<<'\n';
            for(int row=0;row<y.size();++row)out<<y.value[row]<<'\n';out.close();require(bool(out),"cannot archive actual midpoint flux input");
        };
        auto result=advanceParametricCapillaryStep(mesh,domain,hydro,state,time,options,observer);state=std::move(result.interface);
        if(step%std::max(1,steps/8)==0 || step==steps)writeCurve(step,time.end);
        buildParametricCapillaryLoad(mesh,domain,pressure,state,1,&surface);
        const double shape=shapeMoment(state.parametricGeometry().geometry()),slope=shape-oldShape;
        if(shape*oldShape<0)++signCrossings;if(slope*previousSlope<0)++turningPoints;previousSlope=slope;oldShape=shape;
        dissipation-=time.dt*result.hydro.physicalViscousPower;temporalMass+=time.dt*result.hydro.midpointMassEnergyDefect;
        finalEnergy=result.kineticEnergy+surface.length;maximumIncrease=std::max(maximumIncrease,finalEnergy/initialEnergy-1);
        maxDissipation=std::max(maxDissipation,result.hydro.physicalViscousPower);maxMass=std::max(maxMass,result.hydro.conservativeMassRelativeClosure);
        maxFraction=std::max(maxFraction,state.parametricGeometry().report().maximumInventoryFractionError);
        const double volume=std::accumulate(state.parametricGeometry().inventory().begin(),state.parametricGeometry().inventory().end(),0.);
        maxVolume=std::max(maxVolume,std::abs(volume-initialVolume));maxDefect=std::max(maxDefect,result.midpointVelocityDefect);
        maxIterations=std::max(maxIterations,result.iterations);totalIterations+=result.iterations;
        std::cout<<std::setprecision(17)<<"HYDRO_PARAMETRIC_MIDPOINT_STEP step="<<step<<" time="<<time.end<<" shape_moment="<<shape
            <<" surface_length="<<surface.length<<" kinetic_energy="<<result.kineticEnergy<<" total_energy="<<finalEnergy
            <<" dissipated_energy="<<dissipation<<" mass_time_defect="<<temporalMass<<" energy_balance="<<finalEnergy+dissipation-initialEnergy-temporalMass
            <<" iterations="<<result.iterations<<" velocity_defect="<<result.midpointVelocityDefect<<" arcs="<<state.parametricGeometry().geometry().arcs.size()<<'\n'<<std::flush;
    }
    require(signCrossings>=2 && turningPoints>=2 && state.progress().completedSteps==steps && hydro.state().progress.acceptedSteps==steps,
        "requested capillary interval did not complete an oscillation and shared physical history");
    require(maxVolume<1e-13 && maxFraction<3e-12 && maxDissipation<=0,"midpoint motion violated inventory or viscous dissipation");
    std::cout<<std::setprecision(17)<<"HYDRO_PARAMETRIC_MIDPOINT_OSCILLATION grid="<<n<<" steps="<<steps<<" dt="<<dt<<" end_time="<<endTime
        <<" sign_crossings="<<signCrossings<<" turning_points="<<turningPoints<<" final_shape_moment="<<oldShape
        <<" initial_energy="<<initialEnergy<<" final_energy="<<finalEnergy<<" dissipated_energy="<<dissipation<<" mass_time_defect="<<temporalMass
        <<" energy_balance="<<finalEnergy+dissipation-initialEnergy-temporalMass<<" maximum_energy_increase_relative="<<maximumIncrease
        <<" inventory_drift="<<maxVolume<<" max_fraction_error="<<maxFraction<<" dual_mass_relative_closure="<<maxMass
        <<" max_iterations="<<maxIterations<<" total_iterations="<<totalIterations<<" max_midpoint_defect="<<maxDefect
        <<" full_discrete_energy_theorem=0 spatial_convergence_claim=0\n";
}
}
void hydroParametricMidpointContract(int steps){if(steps==0)midpointOperators();else midpointOscillation(steps);}

// Replays only the captured interface flux preparation. No reseeding or
// modification of the archived FV inventory, and no new fluid solve.
void hydroParametricFluxReplay(const char* path){
    std::ifstream in(path);std::string header;std::getline(in,header);
    require(header=="QUANTUM_MIDPOINT_FLUX_INPUT_V1","wrong midpoint flux archive");
    int n=0,step=0,iteration=0;double begin=0,dt=0;in>>n>>step>>iteration>>begin>>dt;
    require(n>=4 && n<=256 && dt>0,"invalid archived interval");
    std::size_t count=0;in>>count;require(count>0 && count<100000,"invalid archived curve size");P::Curve old;old.arcs.resize(count);
    for(auto& arc:old.arcs){in>>arc.a>>arc.b;for(auto& p:arc.p)in>>p.x>>p.y;}
    in>>count;require(count==std::size_t(n)*n,"invalid archived inventory size");std::vector<double> inventory(count);
    for(auto& value:inventory)in>>value;
    auto input=meshInput(n);input.maxCellNum=2*(n+1)*(n+1);QuadtreeMesh mesh(input);FaceSegmentGraph graph;graph.build(mesh);
    const auto domain=NativeLeafDomain::build(mesh,graph,false);VEdgeField<double> x(mesh,0);HEdgeField<double> y(mesh,0);
    int nx=0,ny=0;in>>nx>>ny;require(nx==x.size() && ny==y.size(),"archived native velocity layout differs");
    for(int row=0;row<nx;++row)in>>x.value[row];for(int row=0;row<ny;++row)in>>y.value[row];require(bool(in),"truncated midpoint archive");
    const HydroSplineVelocity velocity(mesh,domain,x,y);const auto view=HydroInterfaceTransport::normalFlux(velocity);
    const auto flow=HydroInterfaceTransport::frozenParametricFlow(velocity);
    std::array<std::vector<int>,2> rows;for(auto& list:rows)list.resize((n+1)*n,-1);
    for(const auto face:domain.faces(mesh)){const int line=int(std::lround((face.horizontal?face.center.y:face.center.x)*n));
        const int strip=int(std::floor((face.horizontal?face.center.x:face.center.y)*n));rows[face.horizontal][line*n+strip]=face.nativeRow;}
    InterfaceParametricTransport::NativeFlux flux;flux.grid.nx=flux.grid.ny=n;
    flux.integral=[&](bool horizontal,int line,int strip,double a,double b){long double value=0;const int row=rows[horizontal][line*n+strip];
        for(int half=0;half<2;++half){const double origin=(strip+.5*half)/n;
            value+=view.integratedNormalFlux(row,half,std::clamp((a-origin)*2*n,0.,1.),std::clamp((b-origin)*2*n,0.,1.));}return double(value);};
    flux.density=[&](bool horizontal,int line,int strip,double z){const int half=z<=(strip+.5)/n?0:1;
        return double(2*n*view.normalFluxDensity(rows[horizontal][line*n+strip],half,false,std::clamp((z-(strip+.5*half)/n)*2*n,0.,1.)));};
    flux.densityMagnitudeBound=[&](bool horizontal,int line,int strip){const int row=rows[horizontal][line*n+strip];
        return std::nextafter(double(2*n*std::max(view.densityMagnitudeBound(row,0),view.densityMagnitudeBound(row,1))),std::numeric_limits<double>::infinity());};
    flux.polynomial=[&](bool horizontal,int line,int strip,int half){return view.normalFluxPolynomial(rows[horizontal][line*n+strip],half);};
    for(double duration:{.5*dt,dt}){
        InterfaceParametricTransport::SourceReport sourceReport;InterfaceParametricTransport::FluxReport report;InterfaceParametricTransport::FaceVolumes faces;
        std::cout<<std::setprecision(17)<<"MIDPOINT_FLUX_REPLAY step="<<step<<" iteration="<<iteration<<" begin="<<begin<<" duration="<<duration<<'\n'<<std::flush;
        const auto source=InterfaceParametricTransport::transportSource(old,{n,n},duration,flow,sourceReport);
        const auto next=InterfaceParametricTransport::advanceInventory(old,source,flow,{flux,flux},duration,inventory,report,&faces);
        // Independent endpoint geometry, rather than the face-time charts,
        // checks the physical volume increment. Preserve inherited inventory
        // error by comparing increments, not by overwriting the old ledger.
        P::PhaseMeasure initialPhase(old);P::PolynomialPhaseMeasure<9> finalPhase(source);
        double incrementError=0;
        for(int i=0;i<n;++i)for(int j=0;j<n;++j){
            const double x0=double(i)/n,x1=double(i+1)/n,y0=double(j)/n,y1=double(j+1)/n;
            const double change=finalPhase.simpleRectangle(x0,x1,y0,y1)-initialPhase.simpleRectangle(x0,x1,y0,y1);
            incrementError=std::max(incrementError,std::abs(next[i*n+j]-inventory[i*n+j]-change)*n*n);}
        std::cout<<"MIDPOINT_FLUX_REPLAY_RESULT balance="<<std::accumulate(next.begin(),next.end(),0.)-std::accumulate(inventory.begin(),inventory.end(),0.)
            <<" crossing_residual="<<report.crossingResidual<<" quadrature_error="<<report.quadratureError
            <<" refinements="<<report.quadratureRefinements<<" endpoint_increment_fraction_error="<<incrementError<<'\n'<<std::flush;
        require(incrementError<3e-12,"frozen face flux differs from independent endpoint volume increment");
    }
}
