#include "Modules/Hydrodynamics/HydroInterfaceNormalLoad.h"
#include "Modules/Hydrodynamics/HydroMomentumGeometry.h"
#include "Modules/Hydrodynamics/HydroSinglePhaseSaddle.h"
#include "Modules/Hydrodynamics/HydroSpatialAdvector.h"
#include "interface_first_jet_fixture.h"

namespace {
using namespace LAMP::Hydrodynamics;
using Matrix=HydroSaddleVelocityBlock::Matrix;
int scenes=0,guards=0,coarseFaces=0,negativeControls=0;
double peakLoad=0,peakMean=0,peakVelocity=0,peakPressure=0;
double relative(double e,double scale){if(scale==0){check(e==0,"nonzero unscaled error");return 0;}check(std::isfinite(e)&&std::isfinite(scale),"nonfinite independent error");return std::abs(e)/scale;}
HydroPhysicalDivergence metric(const QuadtreeMesh& mesh,bool rz,bool selected=false){FaceSegmentGraph graph;graph.build(mesh,[selected](const QuadtreeMesh&,const Cell& c){return !selected||c.centerCoord().x<.75;});return HydroPhysicalDivergence::build(mesh,HydroIntegralDivergence::build(mesh,HydroDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,rz))));}
HermiteGalsPublication publish(const QuadtreeMesh& mesh,const Profile& profile,bool rz,const std::vector<LevelSetInterfacePanels::SymmetryPlane>& planes={}){
    const auto seed=stateFor(mesh,profile,rz,{},[](int,Vector2D,Vector2D){return 0.;},planes,HermiteGalsTransport::FirstJet);
    HermiteGalsCandidate candidate(mesh,{},seed.primary());const auto id=candidate.identity();const auto mixed=closeHermiteGalsMixedAmr(mesh,candidate,planes);
    HermiteGalsPublicationConfig cfg;cfg.volume.symmetryPlanes=planes;return publishHermiteGalsFirstJetGeometry(mesh,std::move(candidate),mixed,{id,centers(mesh,profile)},cfg);
}
HydroSaddleVelocityBlock massBlock(const HydroMomentumMass& m,double dt){Matrix a(m.layout().velocityCount(),m.layout().velocityCount());for(int u=0;u<a.rows();++u)a.coeffRef(u,u)=m.masses()[u]/dt;return HydroSaddleVelocityBlock::build(m.layout(),a);}
std::vector<int> wallDofs(const QuadtreeMesh& mesh,const HydroMomentumGeometry& geometry){std::vector<int> ids;const auto& l=geometry.layout();
    for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(mesh,l.nativeFace(u));if((f.negativeCell<0||f.positiveCell<0)&&geometry.metrics().faceMeasures()[u]>0)ids.push_back(u);}return ids;}
struct SourceContext {
    HydroMomentumGeometry geometry;HydroAdvectorReconstruction source;HydroAdvectorBoundary law;
    SourceContext(const QuadtreeMesh& mesh,const HydroPhaseGeometry& phase):geometry(HydroMomentumGeometry::build(mesh,phase.metrics())),
        source(HydroAdvectorReconstruction::prepare(mesh,geometry)),law(HydroAdvectorBoundary::impermeable(mesh,geometry,wallDofs(mesh,geometry))){}
};
void inspect(const QuadtreeMesh& mesh,const HydroPhaseGeometry& phase,double exactCurvature=-1){
    ++scenes;const auto& l=phase.layout();const auto pressure=HydroPressureLoad::build(phase.metrics().divergence());
    const SourceContext context(mesh,phase);const auto& source=context.source;const auto& law=context.law;
    const auto cap=HydroInterfaceNormalLoad::capillary(mesh,pressure,phase,.7,source,&law);check(cap.isCapillary()&&cap.surfaceTension()==.7,"capillary kind");
    cap.validateFor(mesh,pressure,phase.geometryIdentity(),phase.timeLayer(),source,&law);
    const auto zero=HydroInterfaceNormalLoad::capillary(mesh,pressure,phase,0.,source,&law);for(double f:zero.balanceLoads())check(f==0,"zero sigma force");
    const auto& pub=phase.publication();std::vector<double> chi(l.pressureCount());
    for(int p=0;p<l.pressureCount();++p){
        const int native=l.nativeCell(p);const auto x=mesh.cells(mesh.cellIndexToLabel[native]).centerCoord();
        const double value=pub.state().geometry().snapshot().scalarJetInNativeCell(native,x).value;
        chi[p]=value<0?1:value>0?0:.5;
        check(cap.phaseIndicators()[p]==chi[p],"same-scalar center phase");
    }
    double error=0,scale=0,averageError=0,averageScale=0;
    for(int u=0;u<l.velocityCount();++u){
        const auto f=l.domain().face(mesh,l.nativeFace(u));const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
        const auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
        const double area=(f.horizontal?std::abs(b.x-a.x):std::abs(b.y-a.y))*(l.domain().axisymmetric()?double(2*pi)*(f.horizontal?(a.x+b.x)/2:a.x):1);
        const bool internal=f.negativeCell>=0&&f.positiveCell>=0&&area>0;
        const double difference=internal?chi[l.pressureDof(f.positiveCell)]-chi[l.pressureDof(f.negativeCell)]:0;
        if(difference==0){
            check(cap.balanceLoads()[u]==0&&cap.linkCurvatures()[u]==0&&cap.linkLengths()[u]==0,"same-phase or boundary capillary load");
        }else{
            const auto x=mesh.cells(mesh.cellIndexToLabel[f.negativeCell]).centerCoord(),y=mesh.cells(mesh.cellIndexToLabel[f.positiveCell]).centerCoord();
            close(cap.linkLengths()[u],std::hypot(y.x-x.x,y.y-x.y),1e-14,"actual pressure-center link");
            check(std::isfinite(cap.linkCurvatures()[u]),"finite native link curvature");
            if(exactCurvature>=0){
                close(cap.linkCurvatures()[u],exactCurvature,1e-10*std::max(1.,exactCurvature),"analytic circle/sphere/plane link curvature");
                const double expected=-.7*exactCurvature*area*difference;
                // Every fixture has a unit-length domain. A planar reference
                // has zero traction, so use its physical unit-curvature load
                // to normalize roundoff instead of dividing by analytic zero.
                const double referenceCurvature=exactCurvature>0?exactCurvature:1.;
                error=std::hypot(error,cap.balanceLoads()[u]-expected);
                scale=std::hypot(scale,.7*referenceCurvature*area*difference);
                averageError=std::hypot(averageError,cap.linkCurvatures()[u]-exactCurvature);
                averageScale=std::hypot(averageScale,referenceCurvature);
            }
        }
        if(internal&&mesh.cells(mesh.cellIndexToLabel[f.negativeCell]).level()!=mesh.cells(mesh.cellIndexToLabel[f.positiveCell]).level())++coarseFaces;
        check(cap.constantReactions()[u]==0,"capillary excludes prescribed-wall reaction");
    }
    const double e=relative(error,scale),av=relative(averageError,averageScale);
    check(e<=1e-12&&av<=1e-12,"analytic native traction");peakLoad=std::max(peakLoad,e);peakMean=std::max(peakMean,av);
    const auto incidence=pressure.apply(mesh,HydroNativeValues(l,std::vector<double>(l.velocityCount()),chi));
    const Eigen::Map<const Eigen::VectorXd> values(incidence.data(),incidence.size());const Eigen::VectorXd factored=cap.stressMatrix()*values;
    double factorError=0,factorScale=0;
    for(int u=0;u<l.velocityCount();++u){factorError=std::hypot(factorError,factored[u]-cap.balanceLoads()[u]);factorScale=std::hypot(factorScale,cap.balanceLoads()[u]);}
    check(relative(factorError,factorScale)<=1e-12,"native stress matrix pressure factorization");
    for(int outer=0;outer<cap.stressMatrix().outerSize();++outer)for(Matrix::InnerIterator it(cap.stressMatrix(),outer);it;++it)check(it.row()==it.col(),"capillary contains nonlocal signed stress");
    ++negativeControls;
    if(exactCurvature>0){const auto mass=HydroMomentumMass::build(mesh,context.geometry,phase,1000,1);
        const auto block=massBlock(mass,.125);const auto saddle=HydroSinglePhaseSaddle::buildClosed(pressure,block);auto rhs=cap.balanceLoads();for(auto& v:rhs)v=-v;
        const auto result=saddle.solve(mesh,HydroSaddleRhs(l,rhs,std::vector<double>(l.velocityCount())));double sum=0,volume=0,maxU=0,scaleU=0,maxP=0;
        for(int p=0;p<l.pressureCount();++p){sum+=phase.metrics().cellVolumes()[p]*chi[p];volume+=phase.metrics().cellVolumes()[p];}
        const double jump=.7*exactCurvature;for(int p=0;p<l.pressureCount();++p)maxP=std::max(maxP,std::abs(result.values().pressure()[p]-jump*(chi[p]-sum/volume)));
        for(int u=0;u<l.velocityCount();++u){maxU=std::max(maxU,std::abs(result.values().normalVelocity()[u]));scaleU=std::max(scaleU,std::abs(cap.balanceLoads()[u])/block.matrix().coeff(u,u));}
        const double ev=relative(maxU,scaleU),ep=maxP/jump;check(ev<=1e-11&&ep<=1e-10,"true Hermite sphere static balance");peakVelocity=std::max(peakVelocity,ev);peakPressure=std::max(peakPressure,ep);
        double pressureAlpha=0;for(int p=0;p<l.pressureCount();++p)pressureAlpha+=phase.metrics().cellVolumes()[p]*result.values().pressure()[p]*(chi[p]-sum/volume);
        check(pressureAlpha>0,"liquid interior must have positive pressure jump");++negativeControls;}
    check(phase.publication().state().primary().phaseVolumes.empty(),"capillary introduced inventory");
    std::cout<<"CAPILLARY,"<<scenes<<','<<l.pressureCount()<<','<<l.velocityCount()<<','<<e<<','<<av<<std::endl;
}
void wave(int n,bool amr){
    auto mesh=makeMesh(n,amr?1:0);const double amplitude=1e-4,k=double(2*pi);
    const Profile shape=[=](Vector2D p){return Jet{p.y-.5-amplitude*std::cos(k*p.x),(p.x==0||p.x==1)?0:amplitude*k*std::sin(k*p.x),1,0};};
    const auto phase=HydroPhaseGeometry::build(*mesh,metric(*mesh,false),publish(*mesh,shape,false,{{0,0},{0,1}}));
    const auto& l=phase.layout();const auto pressure=HydroPressureLoad::build(phase.metrics().divergence());const SourceContext context(*mesh,phase);
    const auto cap=HydroInterfaceNormalLoad::capillary(*mesh,pressure,phase,1,context.source,&context.law);
    const auto& geometry=context.geometry;const auto mass=HydroMomentumMass::build(*mesh,geometry,phase,1,1);
    const auto saddle=HydroSinglePhaseSaddle::buildClosed(pressure,massBlock(mass,1));auto rhs=cap.balanceLoads();for(auto& x:rhs)x=-x;
    const auto result=saddle.solve(*mesh,HydroSaddleRhs(l,rhs,std::vector<double>(l.velocityCount())));
    const auto& reconstruction=context.source;const auto field=reconstruction.evaluate(*mesh,result.values(),phase.geometryIdentity(),phase.timeLayer(),context.law);
    double modal=0,normalizer=0,power=0;const auto& pub=phase.publication();
    for(std::size_t c=0;c<pub.topology().curves.size();++c)for(std::size_t g=0;g<pub.topology().curves[c].graphs().size();++g)for(const auto& q:pub.quadrature(c,g)){
        const auto p=q.sample.frame.point;const auto a=field.cellVelocity(l.pressureDof(q.sample.frame.nativeCellIndex),p);const double slope=-q.sample.frame.normal.x/q.sample.frame.normal.y;
        const double dx=q.physicalWeight*q.sample.frame.normal.y,mode=std::cos(k*p.x);modal+=dx*mode*(a.y-slope*a.x);normalizer+=dx*mode*mode;}
    const auto psi=[&](Vector2D p){return -std::sin(k*p.x)*std::sinh(k*std::min(p.y,1-p.y))/(k*std::sinh(k*.5));};
    for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(*mesh,l.nativeFace(u));const auto& e=f.horizontal?mesh->hEdges(f.edgeLabel):mesh->vEdges(f.edgeLabel);auto a=mesh->verts(e.edgeToVert[0]).coord(),b=mesh->verts(e.edgeToVert[1]).coord();if(f.horizontal?a.x>b.x:a.y>b.y)std::swap(a,b);
        const double mode=f.horizontal?-(psi(b)-psi(a))/(b.x-a.x):(psi(b)-psi(a))/(b.y-a.y);power+=mode*rhs[u];}
    const double omegaSquared=-modal/(amplitude*normalizer),reference=k*k*k*std::tanh(k*.5)/2;
    check(omegaSquared>0&&power<0,"capillary restoring sign and virtual work");negativeControls+=2;
    std::cout<<"WAVE_ACCELERATION,"<<n<<','<<amr<<','<<omegaSquared<<','<<reference<<','<<std::abs(omegaSquared/reference-1)<<','<<power<<std::endl;
}
void ellipseLinks(){
    auto mesh=makeMesh(16);const double a=.168,b=.1523809523809524;
    const Profile ellipse=[=](Vector2D p){const double x=p.x-.5,y=p.y-.5;return Jet{x*x/(a*a)+y*y/(b*b)-1,2*x/(a*a),2*y/(b*b),0};};
    const auto phase=HydroPhaseGeometry::build(*mesh,metric(*mesh,false),publish(*mesh,ellipse,false));
    const auto pressure=HydroPressureLoad::build(phase.metrics().divergence());const SourceContext context(*mesh,phase);
    const auto cap=HydroInterfaceNormalLoad::capillary(*mesh,pressure,phase,341.642,context.source,&context.law);const auto& l=phase.layout();
    int active=0,pureGas=0;
    for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(*mesh,l.nativeFace(u));if(f.negativeCell<0||f.positiveCell<0)continue;
        if(cap.linkLengths()[u]==0)continue;++active;
        const auto first=mesh->cells(mesh->cellIndexToLabel[f.negativeCell]).centerCoord(),last=mesh->cells(mesh->cellIndexToLabel[f.positiveCell]).centerCoord(),d=last-first;
        const double x=first.x-.5,y=first.y-.5,A=d.x*d.x/(a*a)+d.y*d.y/(b*b),B=2*(x*d.x/(a*a)+y*d.y/(b*b)),C=x*x/(a*a)+y*y/(b*b)-1;
        const double discriminant=std::sqrt(B*B-4*A*C);double t=(-B-discriminant)/(2*A);if(t<0||t>1)t=(-B+discriminant)/(2*A);
        check(t>=0&&t<=1,"analytic ellipse pressure-link zero");const auto q=first+d*t;
        const double theta=std::atan2((q.y-.5)/b,(q.x-.5)/a),length=std::hypot(a*std::sin(theta),b*std::cos(theta)),curvature=a*b/(length*length*length);
        close(cap.linkCurvatures()[u],curvature,1e-10,"ellipse uses actual link crossing curvature");
    }
    const auto mass=HydroMomentumMass::build(*mesh,context.geometry,phase,100,1);
    for(int u=0;u<l.velocityCount();++u)if(mass.phaseVolumes()[u].liquid==0){check(cap.balanceLoads()[u]==0,"pure-gas dual receives no displaced capillary jump");++pureGas;}
    check(active>0&&pureGas>0,"ellipse critical force support");negativeControls+=2;
    // Also exercise phi==0 exactly at a pressure center, with the prescribed
    // half indicator and an endpoint zero on the same accepted scalar.
    const Profile centered=[](Vector2D p){return Jet{(p.x-.5)*(p.x-.5)+(p.y-.375)*(p.y-.375)-.125*.125,2*(p.x-.5),2*(p.y-.375),0};};
    auto endpointMesh=makeMesh(4);const auto endpointPhase=HydroPhaseGeometry::build(*endpointMesh,metric(*endpointMesh,false),publish(*endpointMesh,centered,false));
    inspect(*endpointMesh,endpointPhase,8);
}
void provenance(){auto mesh=makeMesh(4);Case s;s.pureFirstJet=true;const auto phase=HydroPhaseGeometry::build(*mesh,metric(*mesh,false),initial(*mesh,s));const auto pressure=HydroPressureLoad::build(phase.metrics().divergence());const SourceContext context(*mesh,phase);const auto cap=HydroInterfaceNormalLoad::capillary(*mesh,pressure,phase,1,context.source,&context.law),copy=cap;
    check(&cap.balanceLoads()==&copy.balanceLoads()&&&cap.geometry().publication().state().geometry().snapshot()==&phase.publication().state().geometry().snapshot(),"immutable shared source");
    const auto reject=[&](auto action,const char* message){rejects(action,message);++guards;};
    for(double sigma:{-1.,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()})reject([&]{HydroInterfaceNormalLoad::capillary(*mesh,pressure,phase,sigma,context.source,&context.law);},"surface tension");
    reject([&]{cap.constantJump();},"no prescribed");const auto fixed=HydroInterfaceNormalLoad::constant(*mesh,pressure,phase,-1);reject([&]{fixed.linkCurvatures();},"no curvature");
    reject([&]{cap.validateFor(*mesh,pressure,initial(*mesh,s).state().identity(),phase.timeLayer());},"identity");auto time=phase.timeLayer();time.step=1;reject([&]{cap.validateFor(*mesh,pressure,phase.geometryIdentity(),time);},"time");
    const auto other=HydroPressureLoad::build(phase.metrics().divergence());reject([&]{cap.validateFor(*mesh,other,phase.geometryIdentity(),phase.timeLayer());},"pressure owner");
    const auto different=HydroPressureLoad::build(HydroIntegralDivergence::build(*mesh,phase.layout()));reject([&]{HydroInterfaceNormalLoad::capillary(*mesh,different,phase,1,context.source,&context.law);},"divergence owner");
    const auto otherSource=HydroAdvectorReconstruction::prepare(*mesh,context.geometry);
    reject([&]{cap.validateFor(*mesh,pressure,phase.geometryIdentity(),phase.timeLayer(),otherSource,&context.law);},"source reconstruction owner");
    reject([&]{cap.validateFor(*mesh,pressure,phase.geometryIdentity(),phase.timeLayer(),context.source,nullptr);},"source boundary owner");
    const auto otherLaw=HydroAdvectorBoundary::impermeable(*mesh,context.geometry,wallDofs(*mesh,context.geometry));
    reject([&]{cap.validateFor(*mesh,pressure,phase.geometryIdentity(),phase.timeLayer(),context.source,&otherLaw);},"source boundary owner");
    auto target=makeMesh(4,1);reject([&]{cap.validateFor(*target,pressure,phase.geometryIdentity(),phase.timeLayer());},"mesh");
    const auto transfer=transferHermiteGalsFirstJetMesh(*mesh,phase.publication(),*target,HermiteGalsMeshTransferKind::ExactRefinement);inspect(*target,HydroPhaseGeometry::build(*target,metric(*target,false),transfer.publication()),1/s.radius);
}
}
int main(){try{std::cout<<std::setprecision(17);
    for(bool rz:{false,true})for(int amr:{0,1,2}){auto mesh=makeMesh(4,amr);Case s;s.pureFirstJet=true;s.radial=rz;
        inspect(*mesh,HydroPhaseGeometry::build(*mesh,metric(*mesh,rz),initial(*mesh,s)),(rz?2:1)/s.radius);
        s.nonpolynomial=true;inspect(*mesh,HydroPhaseGeometry::build(*mesh,metric(*mesh,rz),initial(*mesh,s)));
        s.planar=true;s.planeOffset=.43;inspect(*mesh,HydroPhaseGeometry::build(*mesh,metric(*mesh,rz),initial(*mesh,s)),0);}
    for(double phase:{-1.,1.}){auto mesh=makeMesh(4);const Profile pure=[=](Vector2D){return Jet{phase,0,0,0};};inspect(*mesh,HydroPhaseGeometry::build(*mesh,metric(*mesh,false),publish(*mesh,pure,false)));}
    {auto mesh=makeMesh(8);const Profile pair=[](Vector2D p){const double x=p.x-.3,y=p.y-.4,X=p.x-.7,Y=p.y-.6,q=x*x+y*y-.13*.13,Q=X*X+Y*Y-.13*.13;return Jet{q*Q,2*(x*Q+X*q),2*(y*Q+Y*q),4*(x*Y+X*y)};};
        const auto pub=publish(*mesh,pair,false);check(pub.topology().curves.size()==2,"two real interface components");inspect(*mesh,HydroPhaseGeometry::build(*mesh,metric(*mesh,false),pub));}
    {auto mesh=makeMesh(8,1);Case s;s.pureFirstJet=true;inspect(*mesh,HydroPhaseGeometry::build(*mesh,metric(*mesh,false,true),initial(*mesh,s)),1/s.radius);}
    ellipseLinks();provenance();for(int n:{8,16})for(bool amr:{false,true})wave(n,amr);
    std::cout<<"H13_SUMMARY,"<<scenes<<','<<guards<<','<<coarseFaces<<','<<negativeControls<<','<<peakLoad<<','<<peakMean<<','<<peakVelocity<<','<<peakPressure<<",PASS\n";return 0;
}catch(const std::exception& e){std::cerr<<"H13 FAIL: "<<e.what()<<'\n';return 1;}}
