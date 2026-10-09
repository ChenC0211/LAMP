#include "Modules/Hydrodynamics/HydroCompatibleMomentViscosity.h"
#include "Modules/Hydrodynamics/HydroCompatibleMomentFlow.h"
#include "Modules/Hydrodynamics/HydroStep.h"
#include "hydro_complete_field_observation.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/NativeLeafDomain.h"
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
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
template<class F>void cellIntegral(Vector2D center,Vector2D h,F&& f){
    const auto& g=Quantum::Geometry::Parametric::gaussEight();
    for(int j=0;j<8;j++)for(int i=0;i<8;i++)f(Vector2D{center.x+.5*h.x*g.x[i],center.y+.5*h.y*g.x[j]},.25*h.x*h.y*g.w[i]*g.w[j]);
}
template<class F>void faceIntegral(const QuadtreeMesh& mesh,const NativeLeafFace& face,F&& f){
    const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
    const auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
    const auto& g=Quantum::Geometry::Parametric::gaussEight();
    for(int i=0;i<8;i++){const double t=.5*(1+g.x[i]);f(Vector2D{a.x+t*(b.x-a.x),a.y+t*(b.y-a.y)},.5*face.area*g.w[i]);}
}
struct Exact {Eigen::Vector2d velocity,force;double pressure,mu;Eigen::Vector3d strain;};
Exact exact(double x,double y,bool variable,bool nonlinear=false){
    const double f=x*x*(1-x)*(1-x),fx=2*x-6*x*x+4*x*x*x,fxx=2-12*x+12*x*x,fxxx=-12+24*x;
    const double g=y*y*(1-y)*(1-y),gy=2*y-6*y*y+4*y*y*y,gyy=2-12*y+12*y*y,gyyy=-12+24*y;
    const double amplitude=32,mu=.02*(variable?1+x+y:1),mux=variable?.02:0,muy=mux;
    const double exx=amplitude*fx*gy,eyy=-exx,gamma=amplitude*(f*gyy-fxx*g);
    const double lapx=amplitude*(fxx*gy+f*gyyy),lapy=-amplitude*(fxxx*g+fx*gyy);
    Exact e;e.velocity<<amplitude*f*gy,-amplitude*fx*g;e.strain<<exx,eyy,gamma/std::sqrt(2.);
    e.pressure=.05*(x*x-1./3)*(y*y-1./3);e.mu=mu;
    e.force<<.1*x*(y*y-1./3)-mu*lapx-2*mux*exx-muy*gamma,
        .1*y*(x*x-1./3)-mu*lapy-mux*gamma-2*muy*eyy;
    if(nonlinear)e.force+=1.3*Eigen::Vector2d(e.velocity[0]*exx+e.velocity[1]*amplitude*f*gyy,
        -e.velocity[0]*amplitude*fxx*g+e.velocity[1]*eyy);
    return e;
}
template<class F>Eigen::VectorXd sample(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const Space& space,const Stress& stress,F&& velocity){
    Eigen::VectorXd value=Eigen::VectorXd::Zero(stress.matrix().rows());
    for(const auto face:domain.faces(mesh)){
        const auto& info=space.faces()[face.nativeRow];
        faceIntegral(mesh,face,[&](Vector2D p,double ds){const auto v=velocity(p);const double t=((face.horizontal?p.x:p.y)-info.tangentCenter)/info.tangentScale;
            for(int k=0;k<2;k++){const double q=k?t:1.;value[info.first+k]+=ds/face.area*q*v[face.horizontal?1:0];
                value[stress.firstTangent(face.nativeRow)+k]+=ds/face.area*q*v[face.horizontal?0:1];}});
    }
    for(int c=0;c<mesh.cellIndexNum();c++){
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto center=cell.centerCoord(),h=cell.size();const int row=space.firstCellVelocity(c);
        cellIntegral(center,h,[&](Vector2D p,double da){const auto v=velocity(p);const double x=(p.x-center.x)/h.x,y=(p.y-center.y)/h.y;
            value[row]+=da/(h.x*h.y)*v[0];value[row+1]+=da/(h.x*h.y)*v[1];value[row+2]+=da/(h.x*h.y)*(-y*v[0]+x*v[1]);});
    }
    return value;
}
struct Result {double velocity,pressure,strain;};
Result solve(int n,bool variable){
    DictionaryMesh input;input.minCoord={0,0};input.maxCoord={1,1};input.Nx=n;input.Ny=5*n/4;input.maxCellNum=40000;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);DictionaryMesh::RefineBox outer;outer.min={.25,.25};outer.max={.75,.75};outer.level=1;
    DictionaryMesh::RefineBox inner;inner.min={.4375,.375};inner.max={.625,.5625};inner.level=2;mesh.applyRefineBoxes({outer,inner});
    FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,false);const auto space=Space::build(mesh,domain);
    const int nc=mesh.cellIndexNum(),nb=int(space.mass().rows()),np=int(space.divergence().rows());
    std::vector<Eigen::Matrix3d> viscosity(nc,Eigen::Matrix3d::Zero());
    for(int c=0;c<nc;c++){
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto center=cell.centerCoord(),h=cell.size();
        cellIntegral(center,h,[&](Vector2D p,double da){Eigen::Vector3d q;q<<1,(p.x-center.x)/h.x,(p.y-center.y)/h.y;
            viscosity[c]+=da*exact(p.x,p.y,variable).mu*q*q.transpose();});
    }
    const auto stress=Stress::build(mesh,domain,space,viscosity);stress.validate(mesh,domain,space);
    const auto& k=stress.matrix();const auto& d=space.divergence();const int nv=int(k.rows());
    require((k-Space::Matrix(k.transpose())).norm()<1e-13*k.norm(),"moment strain energy is not symmetric");
    if(n==8&&!variable){
        const auto affine=sample(mesh,domain,space,stress,[](Vector2D p){return Eigen::Vector2d{.2+.3*p.x-.4*p.y,.7+.2*p.x+.1*p.y};});
        const auto rigid=sample(mesh,domain,space,stress,[](Vector2D p){return Eigen::Vector2d{.2-.4*p.y,.7+.4*p.x};});
        const auto quadraticVelocity=[](Vector2D p){const double x=p.x,y=p.y;
            return Eigen::Vector2d{.2+.3*x-.4*y+.2*x*x-.1*x*y+.15*y*y,.7+.2*x+.1*y-.3*x*x+.12*x*y+.11*y*y};};
        const auto quadratic=sample(mesh,domain,space,stress,quadraticVelocity);
        double reconstructionDefect=0,strainDefect=0,divergenceCommutation=0;
        long double quadraticEnergy=0;
        for(int c=0;c<nc;c++){
            const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto center=cell.centerCoord(),h=cell.size();const auto& block=stress.cells()[c];
            Eigen::VectorXd local(block.rows.size());for(int j=0;j<local.size();j++)local[j]=quadratic[block.rows[j]];
            const Eigen::VectorXd r=block.reconstruction*local,e=block.strain*local;
            Eigen::MatrixXd localD=Eigen::MatrixXd::Zero(3,local.size());
            for(int j=0;j<local.size();j++)if(block.rows[j]<nb)for(int m=0;m<3;m++)localD(m,j)=d.coeff(3*c+m,block.rows[j]);
            Eigen::Vector3d gram;gram<<h.x*h.y,h.x*h.y/12,h.x*h.y/12;
            divergenceCommutation=std::max(divergenceCommutation,(gram.asDiagonal()*(block.strain.topRows(3)+block.strain.middleRows(3,3))-localD).norm());
            cellIntegral(center,h,[&](Vector2D p,double da){const double x=(p.x-center.x)/h.x,y=(p.y-center.y)/h.y;
                Eigen::Matrix<double,6,1> q;q<<1,x,y,x*x,x*y,y*y;Eigen::Vector3d s;s<<1,x,y;
                const Eigen::Vector2d rh{r.head<6>().dot(q),r.tail<6>().dot(q)};
                Eigen::Vector3d eh;for(int j=0;j<3;j++)eh[j]=e.segment<3>(3*j).dot(s);
                const Eigen::Vector3d expected{.3+.4*p.x-.1*p.y,.1+.12*p.x+.22*p.y,(-.2-.7*p.x+.42*p.y)/std::sqrt(2.)};
                reconstructionDefect=std::max(reconstructionDefect,(rh-quadraticVelocity(p)).norm());strainDefect=std::max(strainDefect,(eh-expected).norm());
                quadraticEnergy+=da*2*.02*expected.squaredNorm();});
        }
        const double quadraticWork=quadratic.dot(k*quadratic),quadraticRoundoff=256*std::numeric_limits<double>::epsilon()*quadratic.cwiseAbs().dot(k.cwiseAbs()*quadratic.cwiseAbs());
        require(reconstructionDefect<5e-12&&strainDefect<5e-12&&divergenceCommutation<1e-12
            &&std::abs(quadraticWork-double(quadraticEnergy))<quadraticRoundoff,"quadratic viscous reproduction/divergence failed");
        const double exactEnergy=2*.02*(.3*.3+.1*.1+.5*.2*.2),energy=affine.dot(k*affine),rigidResidual=(k*rigid).norm()/k.norm();
        std::cout<<std::setprecision(17)<<"MOMENT_VISCOUS_PATCH_DIAGNOSTIC energy="<<energy<<" exact="<<exactEnergy<<" rigid="<<rigidResidual<<'\n';
        // Rigid translations cancel in K. Bound arithmetic by the absolute
        // quadratic work, not by the much smaller residual strain energy.
        const double roundoff=256*std::numeric_limits<double>::epsilon()*affine.cwiseAbs().dot(k.cwiseAbs()*affine.cwiseAbs());
        require(std::abs(energy-exactEnergy)<roundoff&&rigidResidual<2e-13,"moment affine/rigid viscous patch failed");
        double smallest=1e100;int maximumDofs=0;
        for(const auto& block:stress.cells()){
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(block.dissipation);require(eig.info()==Eigen::Success,"local stress eigensolve failed");
            const auto ev=eig.eigenvalues();require(std::abs(ev[0])<2e-12*ev[ev.size()-1]&&std::abs(ev[2])<2e-12*ev[ev.size()-1]
                &&ev[3]>1e-7*ev[ev.size()-1],"local weak strain has modes beyond three rigid motions");
            smallest=std::min(smallest,ev[3]/ev[ev.size()-1]);maximumDofs=std::max(maximumDofs,int(ev.size()));
        }
        const auto zero=Stress::buildCellwise(mesh,domain,space,std::vector<double>(nc,0));require(zero.matrix().nonZeros()==0,"zero viscosity leaves a penalty");
        const auto constant=Stress::buildCellwise(mesh,domain,space,std::vector<double>(nc,.02));
        require((constant.matrix()-k).norm()<1e-13*k.norm(),"constant material differs between moment and cellwise ports");
        std::cout<<std::setprecision(17)<<"MOMENT_VISCOUS_PATCH energy="<<energy<<" exact="<<exactEnergy<<" rigid_residual="<<rigidResidual
            <<" positive_local_spectral_ratio="<<smallest<<" local_dofs="<<maximumDofs<<" quadratic_reconstruction="<<reconstructionDefect
            <<" quadratic_strain="<<strainDefect<<" divergence_commutation="<<divergenceCommutation<<" zero_viscosity=1\n";
    }
    // Strong homogeneous exterior velocity moments. Pressure has one integral
    // gauge, without dropping a cell divergence equation or pinning its slope.
    std::vector<bool> fixed(nv,false);for(const auto face:domain.faces(mesh))if(face.outwardSign){
        for(int j=0;j<2;j++){fixed[space.faces()[face.nativeRow].first+j]=true;fixed[stress.firstTangent(face.nativeRow)+j]=true;}}
    std::vector<int> free(nv,-1);int count=0;for(int i=0;i<nv;i++)if(!fixed[i])free[i]=count++;
    std::vector<Eigen::Triplet<double>> entries;Eigen::VectorXd rhs=Eigen::VectorXd::Zero(count+np+1),fullForce=Eigen::VectorXd::Zero(nv);
    for(int j=0;j<k.outerSize();j++)for(Space::Matrix::InnerIterator it(k,j);it;++it)if(free[it.row()]>=0&&free[it.col()]>=0)entries.emplace_back(free[it.row()],free[it.col()],it.value());
    for(int j=0;j<d.outerSize();j++)for(Space::Matrix::InnerIterator it(d,j);it;++it)if(free[it.col()]>=0){
        entries.emplace_back(count+it.row(),free[it.col()],-it.value());entries.emplace_back(free[it.col()],count+it.row(),-it.value());}
    for(int c=0;c<nc;c++){
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto center=cell.centerCoord(),h=cell.size();const auto& block=space.cells()[c];
        Eigen::Matrix<double,6,1> force=Eigen::Matrix<double,6,1>::Zero();
        cellIntegral(center,h,[&](Vector2D p,double da){const auto e=exact(p.x,p.y,variable);
            force+=da*Space::affineVelocityBasis({(p.x-center.x)/h.x,(p.y-center.y)/h.y},h).transpose()*e.force;});
        const Eigen::VectorXd localForce=block.affineProjection.transpose()*force;
        for(int j=0;j<localForce.size();j++)fullForce[block.rows[j]]+=localForce[j];
        entries.emplace_back(count+3*c,count+np,h.x*h.y);entries.emplace_back(count+np,count+3*c,h.x*h.y);
    }
    for(int i=0;i<nv;i++)if(free[i]>=0)rhs[free[i]]=fullForce[i];
    Space::Matrix system(count+np+1,count+np+1);system.setFromTriplets(entries.begin(),entries.end());
    Eigen::SuperLU<Space::Matrix> solver;solver.compute(system);require(solver.info()==Eigen::Success,"moment Stokes factorization failed");
    const Eigen::VectorXd solution=solver.solve(rhs);require(solver.info()==Eigen::Success&&solution.allFinite(),"moment Stokes solve failed");
    const double residual=(system*solution-rhs).norm()/std::max(1.,rhs.norm());require(residual<2e-10,"moment Stokes residual failed");
    Eigen::VectorXd velocity=Eigen::VectorXd::Zero(nv);for(int i=0;i<nv;i++)if(free[i]>=0)velocity[i]=solution[free[i]];
    const double divergence=(d*velocity.head(nb)).norm(),dissipation=velocity.dot(k*velocity),power=velocity.dot(fullForce),pressurePower=-solution.segment(count,np).dot(d*velocity.head(nb));
    require(divergence<2e-10&&dissipation>0&&std::abs(dissipation-power)<2e-10*power&&std::abs(pressurePower)<1e-12,"moment Stokes incompressibility or work balance failed");
    long double ue=0,un=0,pe=0,pn=0,se=0,sn=0;
    std::ofstream field;if(n==32)if(const char* prefix=std::getenv("QUANTUM_MOMENT_VISCOUS_FIELDS")){
        field.open(std::string(prefix)+(variable?"_variable.csv":"_constant.csv"));field<<std::setprecision(17)<<"x,y,hx,hy,ux,uy,exact_ux,exact_uy,pressure,exact_pressure\n";}
    for(int c=0;c<nc;c++){
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto center=cell.centerCoord(),h=cell.size();const auto& block=space.cells()[c];
        Eigen::VectorXd local(block.rows.size());for(int j=0;j<local.size();j++)local[j]=velocity[block.rows[j]];
        const Eigen::VectorXd coeff=block.affineProjection*local;const auto& sb=stress.cells()[c];Eigen::VectorXd sl(sb.rows.size());for(int j=0;j<sl.size();j++)sl[j]=velocity[sb.rows[j]];
        const Eigen::VectorXd strain=sb.strain*sl;
        cellIntegral(center,h,[&](Vector2D p,double da){const Vector2D q{(p.x-center.x)/h.x,(p.y-center.y)/h.y};const auto e=exact(p.x,p.y,variable);
            const Eigen::Vector2d uh=Space::affineVelocityBasis(q,h)*coeff;Eigen::Vector3d shape;shape<<1,q.x,q.y;
            const double ph=solution.segment<3>(count+3*c).dot(shape);Eigen::Vector3d sh;
            for(int j=0;j<3;j++)sh[j]=strain.segment<3>(3*j).dot(shape);
            ue+=da*(uh-e.velocity).squaredNorm();un+=da*e.velocity.squaredNorm();pe+=da*(ph-e.pressure)*(ph-e.pressure);pn+=da*e.pressure*e.pressure;
            se+=da*(sh-e.strain).squaredNorm();sn+=da*e.strain.squaredNorm();});
        if(field){const Eigen::Vector2d uh=Space::affineVelocityBasis({0,0},h)*coeff;const auto e=exact(center.x,center.y,variable);
            field<<center.x<<','<<center.y<<','<<h.x<<','<<h.y<<','<<uh[0]<<','<<uh[1]<<','<<e.velocity[0]<<','<<e.velocity[1]<<','<<solution[count+3*c]<<','<<e.pressure<<'\n';}
    }
    Result result{std::sqrt(double(ue/un)),std::sqrt(double(pe/pn)),std::sqrt(double(se/sn))};
    std::cout<<std::setprecision(17)<<"MOMENT_VISCOUS_STOKES n="<<n<<" ny="<<input.Ny<<" variable="<<variable<<" cells="<<nc<<" bulk_dofs="<<nb<<" trace_dofs="<<nv-nb
        <<" velocity_l2="<<result.velocity<<" pressure_l2="<<result.pressure<<" strain_l2="<<result.strain<<" residual="<<residual<<" divergence="<<divergence
        <<" dissipation="<<dissipation<<" external_power="<<power<<" pressure_power="<<pressurePower<<" full_ns=0\n";
    return result;
}
}
void hydroMomentViscosityContract(){
    for(bool variable:{false,true}){Result previous{};for(int n:{8,16,32}){const auto current=solve(n,variable);
        if(previous.velocity){const double u=std::log2(previous.velocity/current.velocity),p=std::log2(previous.pressure/current.pressure),s=std::log2(previous.strain/current.strain);
            std::cout<<std::setprecision(17)<<"MOMENT_VISCOUS_ORDER n="<<n<<" variable="<<variable<<" velocity="<<u<<" pressure="<<p<<" strain="<<s<<'\n';
            require(u>1.7&&p>1.7,"mixed-cell viscous velocity/pressure not approaching second order");}previous=current;}}
}

namespace {
Result manufacturedNS(int n){
    using Flow=Quantum::Hydrodynamics::HydroCompatibleMomentFlow;
    DictionaryMesh input;input.minCoord={0,0};input.maxCoord={1,1};input.Nx=n;input.Ny=5*n/4;input.maxCellNum=40000;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);DictionaryMesh::RefineBox outer;outer.min={.25,.25};outer.max={.75,.75};outer.level=1;
    DictionaryMesh::RefineBox inner;inner.min={.4375,.375};inner.max={.625,.5625};inner.level=2;mesh.applyRefineBoxes({outer,inner});
    FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,false);const auto space=Space::build(mesh,domain);
    const int nc=mesh.cellIndexNum(),nb=int(space.mass().rows());
    const Eigen::VectorXd initial=space.interpolateVelocity(mesh,domain,[](Vector2D p){return exact(p.x,p.y,true,true).velocity;});
    using Material=Quantum::Hydrodynamics::HydroCompatibleMomentMaterial;
    Material::ViscosityField law(nc);for(int c=0;c<nc;c++){
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto x=cell.centerCoord(),h=cell.size();
        for(int i=0;i<3;i++)for(int j=0;j<3;j++)law[c](i,j)=.02*(1+x.x+.5*(i-1)*h.x+x.y+.5*(j-1)*h.y);}
    auto flow=Flow::createFrozenMaterial(mesh,domain,space,Material::buildSinglePhase(mesh,domain,space,1.3,law),initial);Flow::Options options;options.relativeResidualTolerance=1e-11;
    double maxEquation=0,maxDivergence=0,maxBalance=0;int iterations=0;
    for(int step=1;step<=10;step++){
        HydroStepTime t;t.index=step;t.begin=flow.state().time;t.dt=.02;
        auto load=flow.makeBodyForceLoad(mesh,domain,[](Vector2D p){const auto f=exact(p.x,p.y,true,true).force;return Vector2D{f[0],f[1]};},t.begin+.5*t.dt);
        auto candidate=flow.prepareNavierStokesStep(mesh,domain,t,load,options);const auto r=candidate.report();
        maxEquation=std::max(maxEquation,r.equationResidual);maxDivergence=std::max(maxDivergence,r.divergenceResidual);maxBalance=std::max(maxBalance,std::abs(r.energyDefect)/r.kineticBefore);iterations+=r.nonlinearIterations;
        flow.accept(mesh,domain,std::move(candidate));
    }
    const auto& state=flow.state();
    long double ue=0,un=0,pe=0,pn=0,se=0,sn=0;
    std::ofstream field;if(n==32)if(const char* prefix=std::getenv("QUANTUM_MOMENT_NS_OUTPUT")){field.open(std::string(prefix)+"_field.csv");field<<std::setprecision(17)<<"x,y,hx,hy,ux,uy,exact_ux,exact_uy,pressure,exact_pressure\n";}
    for(int c=0;c<nc;c++){
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto center=cell.centerCoord(),h=cell.size();
        HydroCompleteObservation::rectangle(center-h*.5,center+h*.5,[&](Vector2D p,double da){const Vector2D q{(p.x-center.x)/h.x,(p.y-center.y)/h.y};const auto e=exact(p.x,p.y,true,true);
            const auto jet=state.physicalVelocity->firstJet(p),middle=state.midpointPhysicalVelocity->firstJet(p);
            const Eigen::Vector2d u(jet.value.x,jet.value.y);const Eigen::Vector3d shape(1,q.x,q.y);const double pressure=state.pressure.segment<3>(3*c).dot(shape);
            const Eigen::Vector3d s(middle.dx.x,middle.dy.y,(middle.dx.y+middle.dy.x)/std::sqrt(2.));
            ue+=da*(u-e.velocity).squaredNorm();un+=da*e.velocity.squaredNorm();pe+=da*(pressure-e.pressure)*(pressure-e.pressure);pn+=da*e.pressure*e.pressure;
            se+=da*(s-e.strain).squaredNorm();sn+=da*e.strain.squaredNorm();});
        if(field){const auto v=state.physicalVelocity->firstJet(center).value;const Eigen::Vector2d u(v.x,v.y);const auto e=exact(center.x,center.y,true,true);
            field<<center.x<<','<<center.y<<','<<h.x<<','<<h.y<<','<<u[0]<<','<<u[1]<<','<<e.velocity[0]<<','<<e.velocity[1]<<','<<state.pressure[3*c]<<','<<e.pressure<<'\n';}
    }
    Result result{std::sqrt(double(ue/un)),std::sqrt(double(pe/pn)),std::sqrt(double(se/sn))};
    std::cout<<std::setprecision(17)<<"MOMENT_NS_MANUFACTURED n="<<n<<" cells="<<nc<<" steps="<<state.acceptedSteps<<" end_time="<<state.time<<" dt=0.02"
        <<" velocity_l2="<<result.velocity<<" pressure_l2="<<result.pressure<<" strain_l2="<<result.strain<<" equation="<<maxEquation<<" divergence="<<maxDivergence<<" energy_defect="<<maxBalance<<" total_iterations="<<iterations<<'\n';
    require(state.acceptedSteps==10&&maxBalance<1e-9,"nonlinear manufactured trajectory did not complete accurately");return result;
}
}
void hydroMomentNSManufacturedContract(){
    Result previous{};for(int n:{8,16,32}){const auto current=manufacturedNS(n);
        if(previous.velocity){const double u=std::log2(previous.velocity/current.velocity),p=std::log2(previous.pressure/current.pressure),s=std::log2(previous.strain/current.strain);
            std::cout<<std::setprecision(17)<<"MOMENT_NS_SPATIAL_ORDER n="<<n<<" velocity="<<u<<" pressure="<<p<<" strain="<<s<<'\n';
            require(u>1.7&&p>1.7,"actual nonlinear velocity/pressure not approaching second-order space convergence");}previous=current;}
}
