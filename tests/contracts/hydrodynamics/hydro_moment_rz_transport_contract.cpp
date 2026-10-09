#include "Modules/Hydrodynamics/HydroCompatibleMomentFlow.h"
#include "Modules/Hydrodynamics/HydroStep.h"
#include "hydro_complete_field_observation.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/CartesianMeasure.h"
#include "Core/Geometry/ParametricHermiteCurve.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>

namespace {
using Space=Quantum::Hydrodynamics::HydroCompatibleMomentSpace;
using Stress=Quantum::Hydrodynamics::HydroCompatibleMomentViscosity;
using Transport=Quantum::Hydrodynamics::HydroCompatibleMomentTransport;
using Flow=Quantum::Hydrodynamics::HydroCompatibleMomentFlow;
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
Eigen::Vector2d profile(Vector2D p){const double r=p.x,z=p.y,g=z*z*(1-z)*(1-z),gz=2*z-6*z*z+4*z*z*z;
    return {8*(r-2*r*r*r+std::pow(r,5))*gz,8*(-2+8*r*r-6*std::pow(r,4))*g};}
Eigen::VectorXd localValues(const Space::CellBlock& cell,const Eigen::VectorXd& global){Eigen::VectorXd r(cell.rows.size());for(int i=0;i<r.size();i++)r[i]=global[cell.rows[i]];return r;}
template<class F>void integrate(Vector2D center,Vector2D h,F&& f){
    const auto& g=Quantum::Geometry::Parametric::gaussEight();
    for(int panel=0;panel<4;panel++)for(int j=0;j<8;j++)for(int i=0;i<8;i++){
        const Vector2D q{-.5+(panel+.5+.5*g.x[i])/4,.5*g.x[j]};f(Vector2D{center.x+h.x*q.x,center.y+h.y*q.y},.25*h.x*h.y*g.w[i]*g.w[j]/4);}
}
struct Fixture {
    QuadtreeMesh mesh;FaceSegmentGraph graph;NativeLeafDomain domain;Space space;Eigen::VectorXd initial;
    static DictionaryMesh input(int n){DictionaryMesh d;d.minCoord={0,0};d.maxCoord={1,1};d.Nx=n;d.Ny=5*n/4;d.maxCellNum=40000;d.maxInitialIterations=0;return d;}
    Fixture(int n):mesh(input(n)){
        DictionaryMesh::RefineBox outer;outer.min={.25,.25};outer.max={.75,.75};outer.level=1;
        DictionaryMesh::RefineBox inner;inner.min={.4375,.375};inner.max={.625,.5625};inner.level=2;mesh.applyRefineBoxes({outer,inner});
        graph.build(mesh);domain=NativeLeafDomain::build(mesh,graph,true);space=Space::build(mesh,domain);initial=space.interpolateVelocity(mesh,domain,profile);
    }
    Flow owner(bool viscous){
        using Material=Quantum::Hydrodynamics::HydroCompatibleMomentMaterial;
        Material::ViscosityField material(mesh.cellIndexNum(),Eigen::Matrix3d::Zero());
        if(viscous)for(int c=0;c<mesh.cellIndexNum();c++){const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto center=cell.centerCoord(),h=cell.size();
            const double a=center.x-.5*h.x,b=a+h.x,r2[3]={a*a,a*b,b*b};
            for(int i=0;i<3;i++)for(int j=0;j<3;j++)material[c](i,j)=.02*(1+r2[i]+center.y+.5*(j-1)*h.y);}
        return Flow::createFrozenMaterial(mesh,domain,space,Material::buildSinglePhase(mesh,domain,space,1.3,material),initial);
    }
};
double reconstruction(Fixture& f,int resolution){
    const auto transport=Transport::build(f.mesh,f.domain,f.space);
    const auto affine=f.space.interpolateVelocity(f.mesh,f.domain,[](Vector2D p){return Eigen::Vector2d{.3*p.x,.2-.6*p.y+.1*p.x};});
    long double error=0,norm=0;double affineError=0,divError=0,momentError=0,axisError=0,constraint=0,projectionError=0,inputError=0,condition=0,amplification=0;
    for(int c=0;c<f.mesh.cellIndexNum();c++){
        const auto& bulk=f.space.cells()[c];const auto& cell=f.mesh.cells(f.mesh.cellIndexToLabel[c]);const auto center=cell.centerCoord(),h=cell.size();
        const auto initial=localValues(bulk,f.initial),af=localValues(bulk,affine);const auto& reconstructed=transport.cells()[c];constraint=std::max(constraint,reconstructed.constraintDefect);
        condition=std::max(condition,reconstructed.constraintCondition);amplification=std::max(amplification,reconstructed.coefficientAmplification);
        const auto psi=[](Vector2D p){return 8*p.x*p.x*std::pow(1-p.x*p.x,2)*p.y*p.y*std::pow(1-p.y,2);};
        const double gauge=psi({center.x-.5*h.x,center.y-.5*h.y});Eigen::Vector3d analytic=Eigen::Vector3d::Zero();
        integrate(center,h,[&](Vector2D p,double da){const double X=(p.x-center.x)/h.x,Y=(p.y-center.y)/h.y,s=center.x/p.x;
            analytic+=da/(h.x*h.y*std::sqrt(h.x*h.y))*(psi(p)-gauge)/center.x*Eigen::Vector3d(X*s,X*X*s*s,Y*s*s);});
        inputError=std::max(inputError,(analytic-initial.tail(3)).norm());
        Eigen::MatrixXd d(3,initial.size());for(int i=0;i<3;i++)for(int j=0;j<initial.size();j++)d(i,j)=f.space.divergence().coeff(3*c+i,bulk.rows[j]);
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(d,Eigen::ComputeFullV);const Eigen::MatrixXd nullspace=svd.matrixV().rightCols(initial.size()-3);
        Eigen::VectorXd seed(nullspace.cols());for(int i=0;i<seed.size();i++)seed[i]=std::sin(.31*(i+1));const Eigen::VectorXd arbitrary=nullspace*seed;
        Eigen::Matrix<double,6,1> flat=Eigen::Matrix<double,6,1>::Zero(),projected=flat;Eigen::Matrix<double,6,6> gram=Eigen::Matrix<double,6,6>::Zero();
        Eigen::Vector3d retained=Eigen::Vector3d::Zero();
        for(const auto& piece:reconstructed.pieces){
            const auto evaluate=[&](Vector2D q){return Eigen::Vector2d(Transport::physicalBasis(piece,q)*piece.coefficients*arbitrary);};
            const auto flux=[&](Vector2D q){return Eigen::Vector2d((piece.center.x+piece.size.x*q.x)/piece.radialScale*evaluate(q));};
            integrate(piece.center,piece.size,[&](Vector2D p,double da){const Vector2D q{(p.x-piece.center.x)/piece.size.x,(p.y-piece.center.y)/piece.size.y};const double X=(p.x-center.x)/h.x,Y=(p.y-center.y)/h.y;
                const auto basis=Transport::physicalBasis(piece,q);const Eigen::Vector2d u=basis*piece.coefficients*initial,v=evaluate(q),a=basis*piece.coefficients*af;
                affineError=std::max(affineError,(a-Eigen::Vector2d(.3*p.x,.2-.6*p.y+.1*p.x)).norm());
                constexpr double t=.03125;
                const double div=(-flux({q.x-3*t,q.y})[0]+9*flux({q.x-2*t,q.y})[0]-45*flux({q.x-t,q.y})[0]+45*flux({q.x+t,q.y})[0]-9*flux({q.x+2*t,q.y})[0]+flux({q.x+3*t,q.y})[0])/(60*t*piece.size.x)
                    +(-flux({q.x,q.y-3*t})[1]+9*flux({q.x,q.y-2*t})[1]-45*flux({q.x,q.y-t})[1]+45*flux({q.x,q.y+t})[1]-9*flux({q.x,q.y+2*t})[1]+flux({q.x,q.y+3*t})[1])/(60*t*piece.size.y);
                divError=std::max(divError,std::abs(div)*std::min(h.x,h.y));
                const double dw=da*p.x/(center.x*h.x*h.y),dv=CartesianMeasure::twoPi*p.x*da;
                Eigen::Matrix<double,6,1> flatValue;flatValue<<v[0],X*v[0],Y*v[0],v[1],X*v[1],Y*v[1];flat+=da/(h.x*h.y)*flatValue;
                retained+=dw*Eigen::Vector3d(v[0],v[1],-Y*v[0]+X*v[1]);
                const auto p1=Space::affineVelocityBasis({X,Y},h);projected+=dw*p1.transpose()*v;gram+=dw*p1.transpose()*p1;
                error+=dv*(u-profile(p)).squaredNorm();norm+=dv*profile(p).squaredNorm();
            });
            if(piece.axisRegular)for(double y:{-.49,0.,.49})axisError=std::max(axisError,std::abs(evaluate({-.5,y})[0]));
        }
        momentError=std::max({momentError,(retained-arbitrary.tail(6).head(3)).norm(),(flat-bulk.flatVelocityMoments*arbitrary).norm()});
        projectionError=std::max(projectionError,(gram.ldlt().solve(projected)-bulk.affineProjection*arbitrary).norm());
    }
    const double relative=std::sqrt(double(error/norm));
    std::cout<<std::setprecision(17)<<"MOMENT_RZ_RECONSTRUCTION n="<<resolution<<" cells="<<f.mesh.cellIndexNum()<<" velocity_l2="<<relative<<" affine="<<affineError
        <<" divergence="<<divError<<" moments="<<momentError<<" projection="<<projectionError<<" axis="<<axisError<<" constraints="<<constraint
        <<" input_moments="<<inputError<<" condition="<<condition<<" amplification="<<amplification<<'\n';
    require(affineError<2e-10&&divError<2e-9&&momentError<2e-10&&projectionError<2e-10&&axisError==0&&inputError<1e-12,"RZ reconstruction does not preserve one physical moment state");return relative;
}
void algebra(Fixture& f){
    const auto transport=Transport::build(f.mesh,f.domain,f.space);const auto center=transport.assemble(f.mesh,f.domain,f.space,f.initial,1.3,Transport::Flux::Central);
    const auto upwind=transport.assemble(f.mesh,f.domain,f.space,f.initial,1.3,Transport::Flux::TraceUpwind);
    const double skew=(center.central+Space::Matrix(center.central.transpose())).norm()/center.central.norm();
    const auto constant=f.space.interpolateVelocity(f.mesh,f.domain,[](Vector2D){return Eigen::Vector2d{.3,.7};});
    const auto affine=f.space.interpolateVelocity(f.mesh,f.domain,[](Vector2D p){return Eigen::Vector2d{.3+.1*p.x-.2*p.y,.7+.2*p.x+.3*p.y};});
    const double uniform=(center.central*constant).norm()/(center.central.norm()*constant.norm()),penalty=(upwind.dissipation*affine).norm()/(upwind.dissipation.norm()*affine.norm());
    const double energy=f.initial.dot(upwind.dissipation*f.initial),local=upwind.maximumLocalMomentumDefect;
    std::cout<<std::setprecision(17)<<"MOMENT_RZ_TRANSPORT_ALGEBRA skew="<<skew<<" constant="<<uniform<<" affine_penalty="<<penalty<<" dissipation="<<energy<<" local_momentum="<<local<<'\n';
    require(skew<1e-11&&uniform<1e-11&&penalty<1e-11&&energy>0&&local<1e-11,"RZ convection lost conservative energy relation");
}
struct Run {HydroCompleteObservation::Velocity velocity;};
Run trajectory(Fixture& f,int mode,double dt,double endTime,std::ofstream* history){
    auto flow=f.owner(mode==2);Flow::Options options;options.relativeResidualTolerance=1e-11;options.convectionFlux=Transport::Flux::Central;
    const auto initialField=*flow.state().physicalVelocity;
    const double e0=flow.kineticEnergy();double viscous=0,transport=0,maxBalance=0,maxDiv=0,maxResidual=0;int maxIterations=0;
    const int steps=int(std::lround(endTime/dt));
    for(int i=0;i<steps;i++){
        HydroStepTime t;t.index=i+1;t.begin=flow.state().time;t.dt=dt;
        const auto load=flow.makeLoad(f.mesh,f.domain,Eigen::VectorXd::Zero(f.initial.size()),t.begin+.5*dt);
        auto candidate=flow.prepareNavierStokesStep(f.mesh,f.domain,t,load,options);const auto r=candidate.report();flow.accept(f.mesh,f.domain,std::move(candidate));
        viscous+=r.viscousDissipation;transport+=r.transportDissipation;const double balance=std::abs(r.kineticAfter-e0+viscous+transport)/e0;
        maxBalance=std::max(maxBalance,balance);maxDiv=std::max(maxDiv,r.divergenceResidual);maxResidual=std::max(maxResidual,r.equationResidual);maxIterations=std::max(maxIterations,r.nonlinearIterations);
        if(history)*history<<mode<<','<<dt<<','<<i+1<<','<<flow.state().time<<','<<r.kineticAfter<<','<<viscous<<','<<transport<<','<<balance<<','<<r.divergenceResidual<<'\n';
    }
    const double evolution=std::sqrt(HydroCompleteObservation::distanceSquared(f.mesh,f.domain,*flow.state().physicalVelocity,initialField)
        /HydroCompleteObservation::normSquared(f.mesh,f.domain,initialField));
    const double energy=flow.kineticEnergy();
    for(const auto face:f.domain.faces(f.mesh))if(face.area==0){const auto trace=flow.normalPolynomial(f.mesh,f.domain,face.nativeRow);require(trace.value(face.center.y)==0&&trace.volumeRate()==0,"RZ axis publication is not zero");}
    std::cout<<std::setprecision(17)<<"MOMENT_RZ_NS_TRAJECTORY mode="<<mode<<" dt="<<dt<<" steps="<<steps<<" time="<<flow.state().time<<" initial_energy="<<e0<<" final_energy="<<energy
        <<" viscous="<<viscous<<" transport="<<transport<<" energy_defect="<<maxBalance<<" divergence="<<maxDiv<<" equation="<<maxResidual<<" iterations="<<maxIterations<<" velocity_change="<<evolution<<'\n';
    require(flow.state().acceptedSteps==steps&&maxBalance<1e-9&&maxDiv<1e-10&&evolution>1e-5,"RZ nonlinear trajectory failed energy or did not evolve");
    return {*flow.state().physicalVelocity};
}
}
void hydroMomentRZReconstructionContract(){
    double previous=0;for(int n:{8,16,32}){Fixture f(n);const double error=reconstruction(f,n);if(previous){const double order=std::log2(previous/error);std::cout<<std::setprecision(17)<<"MOMENT_RZ_RECONSTRUCTION_ORDER n="<<n<<" order="<<order<<'\n';require(order>1.7,"RZ velocity reconstruction not approaching second order");}previous=error;}
    Fixture f(8);algebra(f);
}
void hydroMomentRZTimeContract(){
    Fixture f(8);std::ofstream history;if(const char* prefix=std::getenv("QUANTUM_MOMENT_RZ_TRANSPORT_OUTPUT")){history.open(std::string(prefix)+"_history.csv");history<<std::setprecision(17)<<"mode,dt,step,time,energy,viscous,transport,energy_defect,divergence\n";}
    for(int mode:{0,2}){std::vector<Run> runs;for(double dt:{.02,.01,.005})runs.push_back(trajectory(f,mode,dt,.2,history?&history:nullptr));
        const double order=.5*std::log2(HydroCompleteObservation::distanceSquared(f.mesh,f.domain,runs[0].velocity,runs[1].velocity)
            /HydroCompleteObservation::distanceSquared(f.mesh,f.domain,runs[1].velocity,runs[2].velocity));
        std::cout<<std::setprecision(17)<<"MOMENT_RZ_NS_TIME_ORDER mode="<<mode<<" order="<<order<<'\n';require(order>1.8&&order<2.2,"RZ NS midpoint not approaching second order in time");}
}
void hydroMomentRZTransportContract(){hydroMomentRZReconstructionContract();hydroMomentRZTimeContract();}

// A bounded local conditioning study, not a completed flow or replacement for
// --moment-rz-transport. It compares equivalent radial streamfunction moment
// coordinates on the same three bubble modes, without changing the field.
void hydroMomentRZCoordinatesContract(){
    const auto& g=Quantum::Geometry::Parametric::gaussEight();
    for(double radius:{2.,4.,8.,16.,32.,64.,128.}){
        Eigen::Matrix3d raw=Eigen::Matrix3d::Zero(),stable=raw;
        for(int panel=0;panel<4;panel++)for(int i=0;i<8;i++){
            const double x=-.5+(panel+.5+.5*g.x[i])/4,w=.125*g.w[i],factor=radius/(radius+x),bubble=.25-x*x;
            const Eigen::Vector3d p(1,x,x*x),unscaled(1,factor,factor*factor),scaled(1,x*factor,x*x*factor*factor);
            raw+=w*bubble*unscaled*p.transpose();stable+=w*bubble*scaled*p.transpose();
        }
        Eigen::JacobiSVD<Eigen::Matrix3d> a(raw),b(stable);const double ca=a.singularValues()[0]/a.singularValues()[2],cb=b.singularValues()[0]/b.singularValues()[2];
        // (1,R/r,R^2/r^2) and (1,X*R/r,X^2*R^2/r^2) span the SAME tests.
        Eigen::Matrix3d transformation;transformation<<1,0,0,radius,-radius,0,radius*radius,-2*radius*radius,radius*radius;
        const double equivalence=(stable-transformation*raw).norm()/stable.norm();
        std::cout<<std::setprecision(17)<<"MOMENT_RZ_COORDINATES radius_over_h="<<radius<<" raw_condition="<<ca<<" stable_condition="<<cb<<" equivalence_defect="<<equivalence<<'\n';
        require(cb<500&&equivalence<1e-9,"stable radial coordinate identity failed");
    }
}
