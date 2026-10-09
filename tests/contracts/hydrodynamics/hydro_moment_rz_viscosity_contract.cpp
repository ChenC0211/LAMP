#include "Modules/Hydrodynamics/HydroCompatibleMomentViscosity.h"
#include "Modules/Hydrodynamics/HydroCompatibleMomentFlow.h"
#include "Modules/Hydrodynamics/HydroStep.h"
#include "hydro_complete_field_observation.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/CartesianMeasure.h"
#include "Core/Geometry/ParametricHermiteCurve.h"
#include <Eigen/SuperLUSupport>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <limits>

namespace {
using Space=Quantum::Hydrodynamics::HydroCompatibleMomentSpace;
using Stress=Quantum::Hydrodynamics::HydroCompatibleMomentViscosity;
const double twoPi=CartesianMeasure::twoPi;
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
template<class F>void cellIntegral(Vector2D center,Vector2D h,F&& f){
    const auto& g=Quantum::Geometry::Parametric::gaussEight();
    for(int j=0;j<8;j++)for(int i=0;i<8;i++)f(Vector2D{center.x+.5*h.x*g.x[i],center.y+.5*h.y*g.x[j]},.25*h.x*h.y*g.w[i]*g.w[j]);
}
template<class F>void faceIntegral(const QuadtreeMesh& mesh,const NativeLeafFace& face,F&& f){
    const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
    const auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
    const auto& g=Quantum::Geometry::Parametric::gaussEight();const double length=face.horizontal?b.x-a.x:b.y-a.y;
    for(int i=0;i<8;i++){const double t=.5*(1+g.x[i]);f(Vector2D{a.x+t*(b.x-a.x),a.y+t*(b.y-a.y)},.5*length*g.w[i]);}
}
struct Exact {Eigen::Vector2d u,f;Eigen::Vector4d e;double p,mu;};
Exact exact(Vector2D p,bool variable){
    const double r=p.x,z=p.y,r2=r*r,r3=r2*r,r4=r2*r2;
    const double q=r-2*r3+r4*r,qr=1-6*r2+5*r4,qt=1-2*r2+r4,b=-2+8*r2-6*r4,br=16*r-24*r3;
    const double g=z*z*(1-z)*(1-z),gz=2*z-6*z*z+4*z*z*z,gzz=2-12*z+12*z*z,gzzz=-12+24*z;
    const double amplitude=8,mu=.02*(variable?1+r2+z:1),mur=variable?.04*r:0,muz=variable?.02:0;
    const double gamma=amplitude*(q*gzz+br*g);Exact e;
    e.u<<amplitude*q*gz,amplitude*b*g;e.e<<amplitude*qr*gz,amplitude*b*gz,amplitude*qt*gz,gamma/std::sqrt(2.);
    e.p=.05*(r2-.5)*(z*z-1./3);e.mu=mu;
    e.f<<.1*r*(z*z-1./3)-mu*amplitude*((-16*r+24*r3)*gz+q*gzzz)-2*mur*e.e[0]-muz*gamma,
        .1*z*(r2-.5)-mu*amplitude*((32-96*r2)*g+b*gzz)-mur*gamma-2*muz*e.e[1];return e;
}
// Independent interpolation of all six independent cell functionals. Do not
// use the operator under test to manufacture its own polynomial input.
template<class F,class G>Eigen::VectorXd sample(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const Space& space,const Stress& stress,F&& velocity,G&& streamfunction){
    Eigen::VectorXd value=Eigen::VectorXd::Zero(stress.matrix().rows());
    for(const auto face:domain.faces(mesh)){
        const auto& info=space.faces()[face.nativeRow];if(info.first<0){require(stress.firstTangent(face.nativeRow)<0,"axis acquired fictitious face unknown");continue;}
        faceIntegral(mesh,face,[&](Vector2D p,double ds){const auto v=velocity(p);const double t=((face.horizontal?p.x:p.y)-info.tangentCenter)/info.tangentScale;
            for(int k=0;k<2;k++){const double q=k?t:1.,dw=twoPi*p.x*ds/face.area;
                value[info.first+k]+=dw*q*v[face.horizontal?1:0];value[stress.firstTangent(face.nativeRow)+k]+=dw*q*v[face.horizontal?0:1];}});
    }
    for(int c=0;c<mesh.cellIndexNum();c++){
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto center=cell.centerCoord(),h=cell.size();const int row=space.firstCellVelocity(c);
        const double volume=domain.cellVolumes()[domain.cellRowByIndex()[c]];
        cellIntegral(center,h,[&](Vector2D p,double da){const auto v=velocity(p);const double x=(p.x-center.x)/h.x,y=(p.y-center.y)/h.y,dv=twoPi*p.x*da/volume;
            value[row]+=dv*v[0];value[row+1]+=dv*v[1];value[row+2]+=dv*(-y*v[0]+x*v[1]);});
        // Independent analytic streamfunction, with the declared cell gauge.
        const double gauge=streamfunction({center.x-.5*h.x,center.y-.5*h.y});
        const auto& g=Quantum::Geometry::Parametric::gaussEight();
        for(int panel=0;panel<4;panel++)for(int j=0;j<8;j++)for(int i=0;i<8;i++){
            const double X=-.5+(panel+.5+.5*g.x[i])/4,Y=.5*g.x[j];const Vector2D p{center.x+h.x*X,center.y+h.y*Y};
            const double s=center.x/p.x,psi=(streamfunction(p)-gauge)/center.x,dw=.25*g.w[i]*g.w[j]/(4*std::sqrt(h.x*h.y));
            value.segment<3>(row+3)+=dw*psi*Eigen::Vector3d(X*s,X*X*s*s,Y*s*s);
        }
    }return value;
}
void patch(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const Space& space,const Stress& stress){
    // An admissible quadratic: div(r*u)=0 and every normal trace is P1.
    // A radial-square axial component has a quadratic horizontal normal
    // trace; on split sides its P1 interpolation is not this exact field.
    const auto velocity=[](Vector2D p){const double r=p.x,z=p.y;return Eigen::Vector2d{r*(.2+.3*r-.4*z),.7+.2*r-.4*z-.9*r*z+.4*z*z};};
    const auto q=sample(mesh,domain,space,stress,velocity,[](Vector2D p){return p.x*p.x*(.2*p.y-.2*p.y*p.y-.35)+p.x*p.x*p.x*(.3*p.y-1./15);});
    const auto axial=sample(mesh,domain,space,stress,[](Vector2D){return Eigen::Vector2d{0,1};},[](Vector2D p){return -.5*p.x*p.x;});
    const auto physical=space.interpolateVelocity(mesh,domain,velocity);
    double reproduction=0,strainError=0,commutation=0,spectral=1,axisDefect=0,momentDefect=0;long double energy=0;
    for(int c=0;c<mesh.cellIndexNum();c++){
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto center=cell.centerCoord(),h=cell.size();const auto& b=stress.cells()[c];
        const double volume=domain.cellVolumes()[domain.cellRowByIndex()[c]];
        Eigen::VectorXd local(b.rows.size());for(int i=0;i<local.size();i++)local[i]=q[b.rows[i]];
        const Eigen::VectorXd r=b.reconstruction*local,e=b.strain*local;
        Eigen::Matrix<double,3,12> owned=Eigen::Matrix<double,3,12>::Zero();
        const double ratio=h.x/center.x;Eigen::Matrix<double,1,6> mean;mean<<1,ratio/12,0,1./12,0,1./12;
        owned.block<1,6>(0,0)=mean;owned.block<1,6>(1,6)=mean;
        owned(2,2)=-1./12;owned(2,4)=-ratio/144;owned(2,6)=ratio/12;owned(2,7)=1./12;owned(2,9)=ratio/80;owned(2,11)=ratio/144;
        Eigen::MatrixXd retained=Eigen::MatrixXd::Zero(3,local.size());retained.block(0,int(space.cells()[c].rows.size())-6,3,3).setIdentity();
        momentDefect=std::max(momentDefect,(owned*b.reconstruction-retained).norm());
        Eigen::Matrix3d gram=Eigen::Matrix3d::Zero();gram.diagonal()<<volume,volume/12,volume/12;gram(0,1)=gram(1,0)=volume*h.x/(12*center.x);
        Eigen::MatrixXd d=Eigen::MatrixXd::Zero(3,local.size());for(int i=0;i<local.size();i++)if(b.rows[i]<space.mass().rows())for(int k=0;k<3;k++)d(k,i)=space.divergence().coeff(3*c+k,b.rows[i]);
        commutation=std::max(commutation,(gram*(b.strain.topRows(3)+b.strain.middleRows(3,3)+b.strain.middleRows(6,3))-d).norm()/d.norm());
        cellIntegral(center,h,[&](Vector2D p,double da){const double x=(p.x-center.x)/h.x,y=(p.y-center.y)/h.y;Eigen::Matrix<double,6,1> poly;poly<<1,x,y,x*x,x*y,y*y;
            const Eigen::Vector2d u(r.head<6>().dot(poly),r.tail<6>().dot(poly));const Eigen::Vector3d s(1,x,y);Eigen::Vector4d eh;
            for(int k=0;k<4;k++)eh[k]=e.segment<3>(3*k).dot(s);
            const Eigen::Vector4d expected(.2+.6*p.x-.4*p.y,-.4-.9*p.x+.8*p.y,.2+.3*p.x-.4*p.y,(.2-.4*p.x-.9*p.y)/std::sqrt(2.));
            reproduction=std::max(reproduction,(u-velocity(p)).norm());strainError=std::max(strainError,(eh-expected).norm());energy+=twoPi*p.x*da*2*.02*expected.squaredNorm();});
        if(center.x==.5*h.x){for(double y:{-.5,0.,.5}){Eigen::Matrix<double,1,12> trace=Eigen::Matrix<double,1,12>::Zero();trace.head<6>()<<1,-.5,y,.25,-.5*y,y*y;
            axisDefect=std::max(axisDefect,(trace*b.reconstruction).norm());}}
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(b.dissipation);require(eig.info()==Eigen::Success,"RZ stress eigensolve failed");
        const auto ev=eig.eigenvalues();require(std::abs(ev[0])<3e-12*ev.tail(1)[0]&&ev[1]>1e-10*ev.tail(1)[0],"RZ stress has modes beyond axial translation");spectral=std::min(spectral,ev[1]/ev.tail(1)[0]);
    }
    const double work=q.dot(stress.matrix()*q),roundoff=512*std::numeric_limits<double>::epsilon()*q.cwiseAbs().dot(stress.matrix().cwiseAbs()*q.cwiseAbs());
    const double input=(q.head(physical.size())-physical).norm()/physical.norm(),rigid=(stress.matrix()*axial).norm()/stress.matrix().norm();
    std::cout<<std::setprecision(17)<<"MOMENT_RZ_VISCOUS_PATCH reconstruction="<<reproduction<<" strain="<<strainError<<" divergence_commutation="<<commutation
        <<" axis_trace="<<axisDefect<<" retained_moments="<<momentDefect<<" positive_spectral_ratio="<<spectral<<" energy="<<work<<" exact_energy="<<double(energy)<<" axial_residual="<<rigid<<" physical_input="<<input<<'\n';
    require(reproduction<1e-10&&strainError<1e-10&&commutation<1e-12&&axisDefect<1e-10&&momentDefect<1e-11&&input<1e-13&&rigid<1e-12&&std::abs(work-double(energy))<roundoff,"RZ weighted quadratic strain/energy patch failed");
    const auto zero=Stress::buildCellwise(mesh,domain,space,std::vector<double>(mesh.cellIndexNum(),0));require(zero.matrix().nonZeros()==0,"zero RZ viscosity leaves residual force");
}
struct Result {double u,p,e;};
Result solve(int n,bool variable,bool nonlinear=false){
    DictionaryMesh input;input.minCoord={0,0};input.maxCoord={1,1};input.Nx=n;input.Ny=5*n/4;input.maxCellNum=40000;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);DictionaryMesh::RefineBox outer;outer.min={.25,.25};outer.max={.75,.75};outer.level=1;
    DictionaryMesh::RefineBox inner;inner.min={.4375,.375};inner.max={.625,.5625};inner.level=2;mesh.applyRefineBoxes({outer,inner});
    FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,true);const auto space=Space::build(mesh,domain);
    const int nc=mesh.cellIndexNum(),nb=int(space.mass().rows()),np=int(space.divergence().rows());std::vector<Eigen::Matrix3d> material(nc,Eigen::Matrix3d::Zero());
    for(int c=0;c<nc;c++){const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto center=cell.centerCoord(),h=cell.size();
        cellIntegral(center,h,[&](Vector2D p,double da){const Eigen::Vector3d q(1,(p.x-center.x)/h.x,(p.y-center.y)/h.y);material[c]+=twoPi*p.x*da*exact(p,variable).mu*q*q.transpose();});}
    const auto stress=Stress::build(mesh,domain,space,material);stress.validate(mesh,domain,space);const auto& k=stress.matrix();const auto& d=space.divergence();const int nv=int(k.rows());
    require((k-Space::Matrix(k.transpose())).norm()<1e-12*k.norm(),"RZ dissipation is not symmetric");
    if(n==8&&!variable&&!nonlinear){patch(mesh,domain,space,stress);const auto constant=Stress::buildCellwise(mesh,domain,space,std::vector<double>(nc,.02));require((constant.matrix()-k).norm()<1e-12*k.norm(),"RZ cellwise material omitted weighted cross moment");}
    Eigen::VectorXd velocity=Eigen::VectorXd::Zero(nv),strainVelocity=velocity,pressureState;
    std::optional<Quantum::Hydrodynamics::HydroCompatibleMomentSpline::NativeVelocity> endpoint,midpoint;
    double residual=0,divergence=0,power=0,dissipation=0,pressurePower=0;
    if(!nonlinear){
    std::vector<bool> fixed(nv,false);for(const auto face:domain.faces(mesh))if(face.outwardSign&&face.area>0)for(int j=0;j<2;j++){fixed[space.faces()[face.nativeRow].first+j]=true;fixed[stress.firstTangent(face.nativeRow)+j]=true;}
    std::vector<int> free(nv,-1);int count=0;for(int i=0;i<nv;i++)if(!fixed[i])free[i]=count++;
    std::vector<Eigen::Triplet<double>> entries;Eigen::VectorXd rhs=Eigen::VectorXd::Zero(count+np+1),force=Eigen::VectorXd::Zero(nv);
    force.head(nb)=space.integrateBodyForce(mesh,domain,[&](Vector2D p){return exact(p,variable).f;});
    for(int j=0;j<k.outerSize();j++)for(Space::Matrix::InnerIterator it(k,j);it;++it)if(free[it.row()]>=0&&free[it.col()]>=0)entries.emplace_back(free[it.row()],free[it.col()],it.value());
    for(int j=0;j<d.outerSize();j++)for(Space::Matrix::InnerIterator it(d,j);it;++it)if(free[it.col()]>=0){entries.emplace_back(count+it.row(),free[it.col()],-it.value());entries.emplace_back(free[it.col()],count+it.row(),-it.value());}
    for(int c=0;c<nc;c++){const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const double volume=domain.cellVolumes()[domain.cellRowByIndex()[c]];
        for(int m=0;m<2;m++){const double integral=m?volume*cell.size().x/(12*cell.centerCoord().x):volume;entries.emplace_back(count+3*c+m,count+np,integral);entries.emplace_back(count+np,count+3*c+m,integral);}}
    for(int i=0;i<nv;i++)if(free[i]>=0)rhs[free[i]]=force[i];
    Space::Matrix system(count+np+1,count+np+1);system.setFromTriplets(entries.begin(),entries.end());Eigen::SuperLU<Space::Matrix> solver;solver.compute(system);require(solver.info()==Eigen::Success,"RZ moment Stokes factorization failed");
    const Eigen::VectorXd solution=solver.solve(rhs);require(solver.info()==Eigen::Success&&solution.allFinite(),"RZ moment Stokes solve failed");
    residual=(system*solution-rhs).norm()/rhs.norm();for(int i=0;i<nv;i++)if(free[i]>=0)velocity[i]=solution[free[i]];
    divergence=(d*velocity.head(nb)).norm();power=velocity.dot(force);dissipation=velocity.dot(k*velocity);pressurePower=-solution.segment(count,np).dot(d*velocity.head(nb));
    require(residual<1e-9&&divergence<1e-9&&dissipation>0&&std::abs(dissipation-power)<1e-9*power&&std::abs(pressurePower)<1e-10,"RZ Stokes residual/energy balance failed");
    pressureState=solution.segment(count,np);strainVelocity=velocity;
    }else{
        using Flow=Quantum::Hydrodynamics::HydroCompatibleMomentFlow;
        constexpr double density=1.3,dt=.02;constexpr int steps=10;
        using Material=Quantum::Hydrodynamics::HydroCompatibleMomentMaterial;
        Material::ViscosityField law(nc);for(int c=0;c<nc;c++){
            const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto x=cell.centerCoord(),h=cell.size();
            const double a=x.x-.5*h.x,b=a+h.x,r2[3]={a*a,a*b,b*b};
            for(int i=0;i<3;i++)for(int j=0;j<3;j++)law[c](i,j)=.02*(variable?1+r2[i]+x.y+.5*(j-1)*h.y:1);}
        auto flow=Flow::createFrozenMaterial(mesh,domain,space,Material::buildSinglePhase(mesh,domain,space,density,law),
            space.interpolateVelocity(mesh,domain,[&](Vector2D p){return exact(p,variable).u;}));
        Flow::Options options;options.relativeResidualTolerance=1e-11;
        const auto force=[&](Vector2D p)->Vector2D{
            const auto e=exact(p,variable);const double r=p.x,z=p.y,g=z*z*(1-z)*(1-z),gzz=2-12*z+12*z*z;
            const double urz=8*(r-2*r*r*r+std::pow(r,5))*gzz,uzr=8*(16*r-24*r*r*r)*g;
            const Eigen::Vector2d f=e.f+density*Eigen::Vector2d(e.u[0]*e.e[0]+e.u[1]*urz,e.u[0]*uzr+e.u[1]*e.e[1]);return {f.x(),f.y()};};
        double balance=0;int iterations=0;
        for(int step=1;step<=steps;step++){
            HydroStepTime time{step,flow.state().time,dt,HydroTimePolicy::Fixed};
            const auto load=flow.makeBodyForceLoad(mesh,domain,force,time.begin+.5*dt);
            auto candidate=flow.prepareNavierStokesStep(mesh,domain,time,load,options);const auto report=candidate.report();
            residual=std::max(residual,report.equationResidual);divergence=std::max(divergence,report.divergenceResidual);iterations=std::max(iterations,report.nonlinearIterations);
            dissipation=report.viscousDissipation/dt;power=report.externalWork/dt;pressurePower=report.pressureWork/dt;
            const double error=report.kineticAfter-report.kineticBefore+report.viscousDissipation+report.transportDissipation+report.centralConvectionWork-report.externalWork-report.pressureWork;
            balance=std::max(balance,std::abs(error)/report.kineticBefore);flow.accept(mesh,domain,std::move(candidate));
        }
        require(flow.state().acceptedSteps==steps&&balance<1e-9,"RZ manufactured trajectory failed energy or final time");
        velocity.head(nb)=flow.state().velocity;strainVelocity.head(nb)=flow.state().midpointVelocity;pressureState=flow.state().pressure;
        for(int i=0;i<int(flow.activeTangentialRows().size());i++)strainVelocity[flow.activeTangentialRows()[i]]=flow.state().viscousTrace[i];
        endpoint=flow.state().physicalVelocity;midpoint=flow.state().midpointPhysicalVelocity;
        std::cout<<std::setprecision(17)<<"MOMENT_RZ_NS_MANUFACTURED_TRAJECTORY n="<<n<<" steps="<<steps<<" time="<<flow.state().time<<" dt="<<dt<<" energy_defect="<<balance<<" iterations="<<iterations<<'\n';
    }
    long double ue=0,un=0,pe=0,pn=0,se=0,sn=0;std::ofstream field;
    if(n==32)if(const char* prefix=std::getenv("QUANTUM_MOMENT_RZ_OUTPUT")){field.open(std::string(prefix)+(variable?"_variable.csv":"_constant.csv"));field<<std::setprecision(17)<<"x,y,hx,hy,ux,uy,exact_ux,exact_uy,pressure,exact_pressure\n";}
    for(int c=0;c<nc;c++){
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto center=cell.centerCoord(),h=cell.size();const auto& b=space.cells()[c];const auto& sb=stress.cells()[c];
        Eigen::VectorXd local(b.rows.size()),sl(sb.rows.size());for(int i=0;i<local.size();i++)local[i]=velocity[b.rows[i]];for(int i=0;i<sl.size();i++)sl[i]=strainVelocity[sb.rows[i]];
        const Eigen::VectorXd coeff=b.affineProjection*local,strain=sb.strain*sl;
        const auto observe=[&](Vector2D p,double da){const Vector2D localPoint{(p.x-center.x)/h.x,(p.y-center.y)/h.y};const auto e=exact(p,variable);const Eigen::Vector3d q(1,localPoint.x,localPoint.y);
            Eigen::Vector2d u=Space::affineVelocityBasis(localPoint,h)*coeff;const double ph=pressureState.segment<3>(3*c).dot(q),dv=twoPi*p.x*da;Eigen::Vector4d eh;for(int i=0;i<4;i++)eh[i]=strain.segment<3>(3*i).dot(q);
            if(endpoint){const auto v=endpoint->firstJet(p),m=midpoint->firstJet(p);u={v.value.x,v.value.y};
                eh<<m.dx.x,m.dy.y,(p.x?m.value.x/p.x:m.dx.x),(m.dx.y+m.dy.x)/std::sqrt(2.);}
            ue+=dv*(u-e.u).squaredNorm();un+=dv*e.u.squaredNorm();pe+=dv*(ph-e.p)*(ph-e.p);pn+=dv*e.p*e.p;se+=dv*(eh-e.e).squaredNorm();sn+=dv*e.e.squaredNorm();};
        if(endpoint){const auto lo=center-h*.5;for(int panel=0;panel<4;panel++)
            HydroCompleteObservation::rectangle({lo.x+panel*h.x/4,lo.y},{lo.x+(panel+1)*h.x/4,lo.y+h.y},observe);}
        else cellIntegral(center,h,observe);
        if(field){Eigen::Vector2d u=Space::affineVelocityBasis({0,0},h)*coeff;if(endpoint){const auto v=endpoint->firstJet(center).value;u={v.x,v.y};}
            const auto e=exact(center,variable);field<<center.x<<','<<center.y<<','<<h.x<<','<<h.y<<','<<u[0]<<','<<u[1]<<','<<e.u[0]<<','<<e.u[1]<<','<<pressureState[3*c]<<','<<e.p<<'\n';}
    }
    Result result{std::sqrt(double(ue/un)),std::sqrt(double(pe/pn)),std::sqrt(double(se/sn))};
    std::cout<<std::setprecision(17)<<(nonlinear?"MOMENT_RZ_NS_MANUFACTURED":"MOMENT_RZ_VISCOUS_STOKES")<<" n="<<n<<" variable="<<variable<<" cells="<<nc<<" bulk_dofs="<<nb<<" trace_dofs="<<nv-nb
        <<" velocity_l2="<<result.u<<" pressure_l2="<<result.p<<" strain_l2="<<result.e<<" residual="<<residual<<" divergence="<<divergence<<" dissipation="<<dissipation<<" external_power="<<power<<" pressure_power="<<pressurePower<<'\n';return result;
}
}
void hydroMomentRZViscosityContract(){
    bool accurate=true;
    for(bool variable:{false,true}){Result previous{};for(int n:{8,16,32}){const auto current=solve(n,variable);if(previous.u){
        const double u=std::log2(previous.u/current.u),p=std::log2(previous.p/current.p),e=std::log2(previous.e/current.e);
        std::cout<<std::setprecision(17)<<"MOMENT_RZ_VISCOUS_ORDER n="<<n<<" variable="<<variable<<" velocity="<<u<<" pressure="<<p<<" strain="<<e<<'\n';accurate=accurate&&u>1.7&&p>1.7;}previous=current;}}
    require(accurate,"RZ moment Stokes velocity/pressure not approaching second order");
}
void hydroMomentRZManufacturedContract(){
    Result previous{};for(int n:{8,16,32}){const auto current=solve(n,true,true);if(previous.u){
        const double u=std::log2(previous.u/current.u),p=std::log2(previous.p/current.p),e=std::log2(previous.e/current.e);
        std::cout<<std::setprecision(17)<<"MOMENT_RZ_NS_MANUFACTURED_ORDER n="<<n<<" velocity="<<u<<" pressure="<<p<<" strain="<<e<<'\n';
        require(u>1.7&&p>1.7,"RZ nonlinear velocity/pressure not approaching second order");}previous=current;}
}
