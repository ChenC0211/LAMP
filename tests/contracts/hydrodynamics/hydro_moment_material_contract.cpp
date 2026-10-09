#include "Modules/Hydrodynamics/HydroCompatibleMomentFlow.h"
#include "Modules/Hydrodynamics/HydroStep.h"
#include "hydro_complete_field_observation.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>

namespace {
namespace P=Quantum::Geometry::Parametric;
using Space=Quantum::Hydrodynamics::HydroCompatibleMomentSpace;
using Material=Quantum::Hydrodynamics::HydroCompatibleMomentMaterial;
using Flow=Quantum::Hydrodynamics::HydroCompatibleMomentFlow;
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
void line(P::Curve& curve,P::Point a,P::Point b,int i,int n){P::Arc arc;arc.a=double(i)/n;arc.b=double(i+1)/n;arc.p[0]=a;arc.p[1]=b-a;curve.arcs.push_back(arc);}
P::Curve rectangle(double a,double b,double c,double d){P::Curve curve;line(curve,{a,c},{b,c},0,4);line(curve,{b,c},{b,d},1,4);line(curve,{b,d},{a,d},2,4);line(curve,{a,d},{a,c},3,4);return curve;}
double cap(double x){return .27+x*(.18+x*(.13+x*(.12+.2*x)));}
P::Curve curvedRegion(){
    P::Curve curve;line(curve,{0,0},{1,0},0,4);line(curve,{1,0},{1,cap(1)},1,4);
    P::Arc arc;arc.a=.5;arc.b=.75;arc.p[0]={1,cap(1)};arc.p[1]={-1,-1.6};arc.p[2]={0,1.69};arc.p[3]={0,-.92};arc.p[4]={0,.2};curve.arcs.push_back(arc);
    line(curve,{0,cap(0)},{0,0},3,4);return curve;
}
template<class F>void rectangleIntegral(double a,double b,double c,double d,F&& f){
    if(!(b>a&&d>c))return;const auto& g=P::gaussEight();
    for(int j=0;j<8;j++)for(int i=0;i<8;i++)f(Vector2D{a+.5*(b-a)*(1+g.x[i]),c+.5*(d-c)*(1+g.x[j])},.25*(b-a)*(d-c)*g.w[i]*g.w[j]);
}
// Independent vertical integration below a monotone quartic graph. Its x
// partition comes from monotonic bisection, not production arc clipping.
template<class F>void capIntegral(double a,double b,double c,double d,F&& f){
    std::vector<double> cuts{a,b};for(double height:{c,d})if(cap(a)<height&&cap(b)>height){
        double lo=a,hi=b;for(int i=0;i<60;i++){const double m=.5*(lo+hi);if(cap(m)<height)lo=m;else hi=m;}cuts.push_back(.5*(lo+hi));}
    std::sort(cuts.begin(),cuts.end());const auto& g=P::gaussEight();
    for(int k=1;k<int(cuts.size());k++)for(int i=0;i<8;i++){
        const double x=.5*(cuts[k]+cuts[k-1]+(cuts[k]-cuts[k-1])*g.x[i]),top=std::min(d,cap(x));if(!(top>c))continue;
        for(int j=0;j<8;j++)f(Vector2D{x,c+.5*(top-c)*(1+g.x[j])},.25*(cuts[k]-cuts[k-1])*(top-c)*g.w[i]*g.w[j]);
    }
}
DictionaryMesh input(int n){DictionaryMesh d;d.minCoord={0,0};d.maxCoord={1,1};d.Nx=n;d.Ny=5*n/4;d.maxCellNum=40000;d.maxInitialIterations=0;return d;}
void refine(QuadtreeMesh& mesh){DictionaryMesh::RefineBox a;a.min={.25,.25};a.max={.75,.75};a.level=1;DictionaryMesh::RefineBox b;b.min={.4375,.375};b.max={.625,.5625};b.level=2;mesh.applyRefineBoxes({a,b});}
void geometryAndMaterial(){
    double aligned=0,curved=0,area=0,additivity=0,affineEnergy=0,constantMass=0,viscousEnergy=0;
    for(bool rz:{false,true}){
        const auto shape=rectangle(0,1,0,1);P::PhaseMeasure measure(shape);
        for(const std::array<double,4> box:{std::array<double,4>{0,1,0,1},{0,1,-1,0},{0,1,1,2},{.2,.8,0,1},{0,1,.5,1},{0,1,-.5,.5}}){
            const auto m=measure.rectangleMoments(box[0],box[1],box[2],box[3],rz);std::array<double,6> exact{};
            rectangleIntegral(std::max(0.,box[0]),std::min(1.,box[1]),std::max(0.,box[2]),std::min(1.,box[3]),[&](Vector2D p,double da){
                const double X=(p.x-.5*(box[0]+box[1]))/(box[1]-box[0]),Y=(p.y-.5*(box[2]+box[3]))/(box[3]-box[2]);const double w=rz?2*std::acos(-1.)*p.x:1;
                const double q[6]={1,X,Y,X*X,X*Y,Y*Y};for(int k=0;k<6;k++)exact[k]+=da*w*q[k];});
            for(int k=0;k<6;k++)aligned=std::max(aligned,std::abs(m[k]-exact[k]));
            if(!rz)area=std::max(area,std::abs(measure.simpleRectangle(box[0],box[1],box[2],box[3])-exact[0]));
        }
        QuadtreeMesh mesh(input(8));refine(mesh);FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,rz);const auto space=Space::build(mesh,domain);
        auto curve=std::make_shared<const P::Curve>(curvedRegion());P::PhaseMeasure phase(*curve);
        Material::Properties properties;properties.liquidDensity=1000;properties.gasDensity=1;properties.liquidViscosity=.2;properties.gasViscosity=.02;
        const auto material=Material::buildFrozen(mesh,domain,space,curve,0,properties);
        double total=0,exactEnergy=0;const auto velocity=[&](Vector2D p)->Eigen::Vector2d{return rz?Eigen::Vector2d(.3*p.x,.2+.5*p.y):Eigen::Vector2d(.2+.3*p.x+.4*p.y,.1-.2*p.x+.5*p.y);};
        for(int c=0;c<mesh.cellIndexNum();c++){
            const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto h=cell.size(),x=cell.centerCoord();const double a=x.x-.5*h.x,b=x.x+.5*h.x,lo=x.y-.5*h.y,hi=x.y+.5*h.y;
            Eigen::Matrix3d reference=Eigen::Matrix3d::Zero();
            capIntegral(a,b,lo,hi,[&](Vector2D p,double da){const double w=rz?2*std::acos(-1.)*p.x:1;const Eigen::Vector3d v(1,(p.x-x.x)/h.x,(p.y-x.y)/h.y);
                reference+=da*w*v*v.transpose();exactEnergy+=999*da*w*velocity(p).squaredNorm();});
            rectangleIntegral(a,b,lo,hi,[&](Vector2D p,double da){exactEnergy+=da*(rz?2*std::acos(-1.)*p.x:1)*velocity(p).squaredNorm();});
            curved=std::max(curved,(material.cells()[c].liquid-reference).norm()/material.cells()[c].volume(0,0));total+=material.cells()[c].liquid(0,0);
            if(!rz)area=std::max(area,std::abs(phase.simpleRectangle(a,b,lo,hi)-reference(0,0)));
        }
        additivity=std::max(additivity,std::abs(total-phase.rectangleMoments(0,1,0,1,rz)[0]));
        const auto u=space.interpolateVelocity(mesh,domain,velocity);const double energy=u.dot(material.kineticMass()*u);affineEnergy=std::max(affineEnergy,std::abs(energy-exactEnergy)/exactEnergy);
        const auto stress=Quantum::Hydrodynamics::HydroCompatibleMomentViscosity::build(mesh,domain,space,material.viscosityMoments());
        require(stress.matrix().nonZeros()>0&&material.geometryHandle()==curve,"material consumers lost geometry identity or viscosity");
        Eigen::VectorXd strainInput=Eigen::VectorXd::Zero(stress.matrix().rows());strainInput.head(u.size())=u;
        const auto& rule=P::gaussEight();
        for(const auto face:domain.faces(mesh))if(face.area>0){
            const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
            const auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();const auto& info=space.faces()[face.nativeRow];
            const double length=face.horizontal?b.x-a.x:b.y-a.y;
            for(int j=0;j<8;j++){const double t=.5*(1+rule.x[j]);const Vector2D p{a.x+t*(b.x-a.x),a.y+t*(b.y-a.y)};
                const double dw=.5*length*rule.w[j]*(rz?2*std::acos(-1.)*p.x:1),theta=((face.horizontal?p.x:p.y)-info.tangentCenter)/info.tangentScale;
                strainInput[stress.firstTangent(face.nativeRow)]+=dw/face.area*velocity(p)[face.horizontal?0:1];
                strainInput[stress.firstTangent(face.nativeRow)+1]+=dw/face.area*theta*velocity(p)[face.horizontal?0:1];}
        }
        const double exactDissipation=(rz?.86:.72)*(.02*(rz?std::acos(-1.):1)+.18*total),actualDissipation=strainInput.dot(stress.matrix()*strainInput);
        viscousEnergy=std::max(viscousEnergy,std::abs(actualDissipation-exactDissipation)/exactDissipation);
        // A zero viscosity in either physical phase follows the same measure;
        // it must not be replaced by a small positive floor.
        for(bool zeroLiquid:{false,true}){
            auto zeroProperties=properties;if(zeroLiquid)zeroProperties.liquidViscosity=0;else zeroProperties.gasViscosity=0;
            const auto zeroMaterial=Material::buildFrozen(mesh,domain,space,curve,0,zeroProperties);
            Quantum::Hydrodynamics::HydroCompatibleMomentViscosity::build(mesh,domain,space,zeroMaterial.viscosityMoments());
        }
        properties.liquidDensity=properties.gasDensity=1.3;const auto uniform=Material::buildFrozen(mesh,domain,space,curve,0,properties);
        constantMass=std::max(constantMass,(uniform.kineticMass()-1.3*space.mass()).norm()/space.mass().norm());
        Eigen::SimplicialLDLT<Space::Matrix> factor;factor.compute(material.kineticMass());require(factor.info()==Eigen::Success&&factor.vectorD().minCoeff()>0,"material mass lost positivity");
        std::cout<<std::setprecision(17)<<"MOMENT_PHASE_MATERIAL rz="<<rz<<" cells="<<mesh.cellIndexNum()<<" liquid_volume="<<total<<" affine_kinetic_norm="<<energy<<" independent_norm="<<exactEnergy<<'\n';
    }
    std::cout<<std::setprecision(17)<<"MOMENT_PHASE_INTEGRALS aligned="<<aligned<<" curved_relative="<<curved<<" area="<<area<<" additive="<<additivity<<" affine_energy_relative="<<affineEnergy<<" viscous_energy_relative="<<viscousEnergy<<" constant_mass="<<constantMass<<'\n';
    require(aligned<2e-13&&curved<2e-11&&area<2e-13&&additivity<2e-12&&affineEnergy<2e-12&&viscousEnergy<2e-11&&constantMass<2e-14,"phase integration or weighted affine mass failed");
}
constexpr double interfaceY=.43;
struct Exact {Eigen::Vector2d velocity,laplacian;};
Exact exact(Vector2D p){
    const double x=p.x,y=p.y,d=y-interfaceY,f=x*x*(1-x)*(1-x),fx=2*x-6*x*x+4*x*x*x,fxx=2-12*x+12*x*x,fxxx=-12+24*x;
    const double a=y*y*(1-y)*(1-y),a1=2*y-6*y*y+4*y*y*y,a2=2-12*y+12*y*y,a3=-12+24*y;
    const double g=a*d*d*d,g1=a1*d*d*d+3*a*d*d,g2=a2*d*d*d+6*a1*d*d+6*a*d,g3=a3*d*d*d+9*a2*d*d+18*a1*d+6*a;
    return {512*Eigen::Vector2d(f*g1,-fx*g),512*Eigen::Vector2d(fxx*g1+f*g3,-fxxx*g-fx*g2)};
}
double pressure(Vector2D p){return .05*(p.x*p.x-1./3)*(p.y*p.y-1./3);}
template<class F>void splitIntegral(Vector2D center,Vector2D h,F&& f){
    // Resolve both the analytic material jump and the half-leaf supports at
    // balanced AMR faces. GL12 integrates the polynomial source/test products.
    const auto lo=center-h*.5;
    for(int i=0;i<2;i++)for(int j=0;j<2;j++){
        const double a=lo.x+.5*i*h.x,b=a+.5*h.x,c=lo.y+.5*j*h.y,d=c+.5*h.y;
        HydroCompleteObservation::rectangle({a,c},{b,std::min(d,interfaceY)},f);
        HydroCompleteObservation::rectangle({a,std::max(c,interfaceY)},{b,d},f);
    }
}
struct Fixture {
    QuadtreeMesh mesh;FaceSegmentGraph graph;NativeLeafDomain domain;Space space;Material material;
    using PointForce=Quantum::Hydrodynamics::HydroCompatibleMomentSpline::PointForce;
    std::vector<PointForce> acceleration,steady;
    explicit Fixture(int n):mesh(input(n)){
        refine(mesh);graph.build(mesh);domain=NativeLeafDomain::build(mesh,graph,false);space=Space::build(mesh,domain);
        Material::Properties properties;properties.liquidDensity=1000;properties.gasDensity=1;properties.liquidViscosity=.2;properties.gasViscosity=.02;
        material=Material::buildFrozen(mesh,domain,space,std::make_shared<const P::Curve>(rectangle(0,1,0,interfaceY)),0,properties);
        for(int c=0;c<mesh.cellIndexNum();c++){
            const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);
            splitIntegral(cell.centerCoord(),cell.size(),[&](Vector2D p,double da){const auto e=exact(p);const bool liquid=p.y<interfaceY;
                const Eigen::Vector2d a=da*(liquid?1000:1)*e.velocity;
                const Eigen::Vector2d b=da*(-(liquid?.2:.02)*e.laplacian+Eigen::Vector2d(.1*p.x*(p.y*p.y-1./3),.1*p.y*(p.x*p.x-1./3)));
                acceleration.push_back({p,{a[0],a[1]}});steady.push_back({p,{b[0],b[1]}});});
        }
    }
};
struct Result {HydroCompleteObservation::Velocity velocity;double u=0,p=0;};
Result trajectory(Fixture& f,int n,double dt,double end,bool field){
    auto flow=Flow::createFrozenMaterial(f.mesh,f.domain,f.space,f.material,Eigen::VectorXd::Zero(f.space.mass().rows()));
    double maxEquation=0,maxDiv=0,maxEnergy=0,dissipation=0,work=0;const int steps=int(std::lround(end/dt));
    std::ofstream history,output;const char* prefix=std::getenv("QUANTUM_MOMENT_MATERIAL_OUTPUT");
    if(prefix&&field){history.open(std::string(prefix)+"_history.csv");history<<std::setprecision(17)<<"step,time,energy,dissipation,external_work,balance\n";}
    for(int step=1;step<=steps;step++){
        HydroStepTime t;t.index=step;t.begin=flow.state().time;t.dt=dt;const double middle=t.begin+.5*dt;
        auto points=f.acceleration;for(std::size_t i=0;i<points.size();i++)
            points[i].measure=f.acceleration[i].measure*(2*std::cos(2*middle))+f.steady[i].measure*std::sin(2*middle);
        const auto load=flow.makeQuadratureLoad(f.mesh,f.domain,points,middle);
        if(step==1){bool rejected=false;try{flow.prepareNavierStokesStep(f.mesh,f.domain,t,load);}catch(const std::invalid_argument&){rejected=true;}
            require(rejected&&flow.state().acceptedSteps==0,"unimplemented variable-density convection was admitted");}
        auto candidate=flow.prepareStokesStep(f.mesh,f.domain,t,load);const auto r=candidate.report();flow.accept(f.mesh,f.domain,std::move(candidate));
        dissipation+=r.viscousDissipation;work+=r.externalWork;const double balance=std::abs(r.kineticAfter+dissipation-work)/std::max(r.kineticAfter,work);
        maxEnergy=std::max(maxEnergy,balance);maxEquation=std::max(maxEquation,r.equationResidual);maxDiv=std::max(maxDiv,r.divergenceResidual);
        if(history)history<<step<<','<<flow.state().time<<','<<r.kineticAfter<<','<<dissipation<<','<<work<<','<<balance<<'\n';
    }
    if(prefix&&field){output.open(std::string(prefix)+"_field.csv");output<<std::setprecision(17)<<"x,y,hx,hy,alpha,ux,uy,exact_ux,exact_uy,pressure,exact_pressure\n";}
    double ue=0,un=0,pe=0,pn=0;const auto& state=flow.state();
    for(int c=0;c<f.mesh.cellIndexNum();c++){
        const auto& cell=f.mesh.cells(f.mesh.cellIndexToLabel[c]);const auto x=cell.centerCoord(),h=cell.size();
        splitIntegral(x,h,[&](Vector2D p,double da){const Vector2D q{(p.x-x.x)/h.x,(p.y-x.y)/h.y};const auto value=state.physicalVelocity->firstJet(p).value;
            const Eigen::Vector2d uh(value.x,value.y),reference=std::sin(2*end)*exact(p).velocity;
            const double ph=state.pressure.segment<3>(3*c).dot(Eigen::Vector3d(1,q.x,q.y)),ref=std::sin(2**state.stageTime)*pressure(p);
            ue+=da*(uh-reference).squaredNorm();un+=da*reference.squaredNorm();pe+=da*(ph-ref)*(ph-ref);pn+=da*ref*ref;});
        if(output){const auto value=state.physicalVelocity->firstJet(x).value;const Eigen::Vector2d uh(value.x,value.y),reference=std::sin(2*end)*exact(x).velocity;
            output<<x.x<<','<<x.y<<','<<h.x<<','<<h.y<<','<<f.material.cells()[c].liquid(0,0)/f.material.cells()[c].volume(0,0)<<','<<uh[0]<<','<<uh[1]<<','<<reference[0]<<','<<reference[1]<<','<<state.pressure[3*c]<<','<<std::sin(2**state.stageTime)*pressure(x)<<'\n';}
    }
    Result result{*state.physicalVelocity,std::sqrt(ue/un),std::sqrt(pe/pn)};
    std::cout<<std::setprecision(17)<<"MOMENT_FROZEN_PHASE_STOKES n="<<n<<" cells="<<f.mesh.cellIndexNum()<<" density_ratio=1000 viscosity_ratio=10 dt="<<dt<<" steps="<<steps<<" end="<<state.time
        <<" velocity_l2="<<result.u<<" pressure_l2="<<result.p<<" equation="<<maxEquation<<" divergence="<<maxDiv<<" energy_defect="<<maxEnergy<<'\n';
    require(state.acceptedSteps==steps&&maxEnergy<1e-9,"frozen phase trajectory did not complete with energy balance");return result;
}
}
void hydroMomentMaterialContract(){geometryAndMaterial();}
void hydroMomentMaterialFlowContract(){
    std::optional<Result> previous;for(int n:{8,16,32}){Fixture f(n);const auto result=trajectory(f,n,.01,.2,n==32);
        if(previous){const double u=std::log2(previous->u/result.u),p=std::log2(previous->p/result.p);
            std::cout<<std::setprecision(17)<<"MOMENT_FROZEN_PHASE_SPACE_ORDER n="<<n<<" velocity="<<u<<" pressure="<<p<<'\n';
            require(u>1.5&&p>1.2,"frozen phase spatial velocity/pressure convergence insufficient");}previous=result;}
    Fixture f(8);std::vector<Result> runs;for(double dt:{.04,.02,.01})runs.push_back(trajectory(f,8,dt,.4,false));
    const auto distance=[&](const Result& a,const Result& b){long double sum=0;
        for(int c:f.domain.cellIndexByRow()){const auto& cell=f.mesh.cells(f.mesh.cellIndexToLabel[c]);
            splitIntegral(cell.centerCoord(),cell.size(),[&](Vector2D p,double da){const auto u=a.velocity.firstJet(p).value,v=b.velocity.firstJet(p).value;
                const long double x=u.x-v.x,y=u.y-v.y;sum+=da*(p.y<interfaceY?1000:1)*(x*x+y*y);});}return double(sum);};
    const double order=.5*std::log2(distance(runs[0],runs[1])/distance(runs[1],runs[2]));
    std::cout<<std::setprecision(17)<<"MOMENT_FROZEN_PHASE_TIME_ORDER order="<<order<<'\n';require(order>1.8&&order<2.2,"frozen phase midpoint time accuracy failed");
}
