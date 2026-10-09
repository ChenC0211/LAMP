#include "Modules/Hydrodynamics/CchInterfaceNormalLoad.h"
#include "Modules/Hydrodynamics/CchMomentumGeometry.h"
#include "Modules/Hydrodynamics/CchSinglePhaseSaddle.h"
#include "Modules/Hydrodynamics/CchSpatialAdvector.h"
#include "Modules/Hydrodynamics/Detail/CchRegularAdvectorGeometry.h"
#include "interface_first_jet_fixture.h"
#include <map>

namespace {
using namespace Quantum::Hydrodynamics;
using Matrix=CchSaddleVelocityBlock::Matrix;
int scenes=0,guards=0,coarseFaces=0,negativeControls=0;
double peakLoad=0,peakMean=0,peakVelocity=0,peakPressure=0;
double relative(double e,double scale){if(scale==0){check(e==0,"nonzero unscaled error");return 0;}check(std::isfinite(e)&&std::isfinite(scale),"nonfinite independent error");return std::abs(e)/scale;}
CchPhysicalDivergence metric(const QuadtreeMesh& mesh,bool rz,bool selected=false){FaceSegmentGraph graph;graph.build(mesh,[selected](const QuadtreeMesh&,const Cell& c){return !selected||c.centerCoord().x<.75;});return CchPhysicalDivergence::build(mesh,CchIntegralDivergence::build(mesh,CchDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,rz))));}
HermiteGalsPublication publish(const QuadtreeMesh& mesh,const Profile& profile,bool rz,const std::vector<LevelSetInterfacePanels::SymmetryPlane>& planes={}){
    const auto seed=stateFor(mesh,profile,rz,{},[](int,Vector2D,Vector2D){return 0.;},planes,HermiteGalsTransport::FirstJet);
    HermiteGalsCandidate candidate(mesh,{},seed.primary());const auto id=candidate.identity();const auto mixed=closeHermiteGalsMixedAmr(mesh,candidate,planes);
    HermiteGalsPublicationConfig cfg;cfg.volume.symmetryPlanes=planes;return publishHermiteGalsFirstJetGeometry(mesh,std::move(candidate),mixed,{id,centers(mesh,profile)},cfg);
}
CchSaddleVelocityBlock massBlock(const CchMomentumMass& m,double dt){Matrix a(m.layout().velocityCount(),m.layout().velocityCount());for(int u=0;u<a.rows();++u)a.coeffRef(u,u)=m.masses()[u]/dt;return CchSaddleVelocityBlock::build(m.layout(),a);}
std::vector<int> wallDofs(const QuadtreeMesh& mesh,const CchMomentumGeometry& geometry){std::vector<int> ids;const auto& l=geometry.layout();
    for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(mesh,l.nativeFace(u));if((f.negativeCell<0||f.positiveCell<0)&&geometry.metrics().faceMeasures()[u]>0)ids.push_back(u);}return ids;}
struct SourceContext {
    CchMomentumGeometry geometry;CchAdvectorReconstruction source;CchAdvectorBoundary law;
    SourceContext(const QuadtreeMesh& mesh,const CchPhaseGeometry& phase):geometry(CchMomentumGeometry::build(mesh,phase.metrics())),
        source(CchAdvectorReconstruction::prepare(mesh,geometry)),law(CchAdvectorBoundary::impermeable(mesh,geometry,wallDofs(mesh,geometry))){}
};
void inspect(const QuadtreeMesh& mesh,const CchPhaseGeometry& phase,double exactCurvature=-1){
    ++scenes;const auto& l=phase.layout();const auto pressure=CchPressureLoad::build(phase.metrics().divergence());
    const SourceContext context(mesh,phase);const auto& source=context.source;const auto& law=context.law;
    const auto cap=CchInterfaceNormalLoad::capillary(mesh,pressure,phase,.7,source,&law);check(cap.isCapillary()&&cap.surfaceTension()==.7,"capillary kind");
    cap.validateFor(mesh,pressure,phase.geometryIdentity(),phase.timeLayer(),source,&law);
    const auto zero=CchInterfaceNormalLoad::capillary(mesh,pressure,phase,0.,source,&law);for(double f:zero.balanceLoads())check(f==0,"zero sigma force");
    // Independently integrate actual unit-velocity virtual work. The source
    // primitive supplies footprint metadata; its values are checked against
    // the real field before the complete weak force is integrated anew.
    const auto& pub=phase.publication();std::vector<HermiteGalsSurfaceQuadrature> nodes;std::vector<double> allStress,p1,p2,jump;std::vector<std::size_t> selectedPoints;
    for(std::size_t c=0;c<pub.topology().curves.size();++c)for(std::size_t g=0;g<pub.topology().curves[c].graphs().size();++g)for(const auto& q:pub.quadrature(c,g)){
        const auto x=q.sample.frame.point;p1.push_back(x.x);p2.push_back(x.x*x.x+.37*x.y*x.y);jump.push_back(x.x*x.x+.2*g-.13*c);
        allStress.push_back(-.7*q.sample.curvature);if(l.pressureDof(q.sample.frame.nativeCellIndex)>=0){selectedPoints.push_back(allStress.size()-1);nodes.push_back(q);}
        check(q.sample.geometryIdentity==phase.geometryIdentity(),"surface source identity");}
    const double anchor=allStress.empty()?0:-allStress.front()/.7;
    std::vector<std::vector<double>> canonical(l.velocityCount(),std::vector<double>(nodes.size()));
    std::vector<std::vector<double>> actualNormal=canonical;
    std::map<int,Quantum::Hydrodynamics::Detail::CchRegularCellStencil> stencils;
    for(std::size_t q=0;q<nodes.size();++q){const auto& f=nodes[q].sample.frame;const int p=l.pressureDof(f.nativeCellIndex);
        auto found=stencils.find(p);if(found==stencils.end())found=stencils.emplace(p,source.regularGeometry().cellStencil(mesh,p,&law)).first;
        for(const auto& row:source.regularGeometry().cellJetWeights(found->second,f.point))canonical[row.velocityDof][q]=row.jet[0]*f.normal.x+row.jet[1]*f.normal.y;}
    double error=0,scale=0,averageError=0,averageScale=0;
    for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(mesh,l.nativeFace(u));const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);const auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
        const double area=(f.horizontal?std::abs(b.x-a.x):std::abs(b.y-a.y))*(l.domain().axisymmetric()?double(2*pi)*(f.horizontal?(a.x+b.x)/2:a.x):1);
        long double weight=0,positive=0,offset=0;double difference=0;
        for(int side:{-1,1}){const int cell=side<0?f.negativeCell:f.positiveCell;if(cell>=0)difference+=side*phase.liquidFractions()[l.pressureDof(cell)];}
        for(std::size_t q=0;q<nodes.size();++q){const double v=canonical[u][q];if(v!=0)weight+=nodes[q].physicalWeight;
            const long double w=nodes[q].physicalWeight*v*v;positive+=w;offset+=w*(nodes[q].sample.curvature-anchor);}
        const double k=positive>0?anchor+double(offset/positive):anchor;long double expected=-.7*k*area*difference;
        const bool prescribed=std::binary_search(law.velocityDofs().begin(),law.velocityDofs().end(),u)||area==0;
        if(!prescribed){std::vector<double> unit(l.velocityCount());unit[u]=1;
            const auto field=source.evaluate(mesh,CchNativeValues(l,unit,std::vector<double>(l.pressureCount())),phase.geometryIdentity(),phase.timeLayer(),law);
            for(std::size_t q=0;q<nodes.size();++q){const auto& frame=nodes[q].sample.frame;const auto v=field.cellVelocity(l.pressureDof(frame.nativeCellIndex),frame.point);const double normal=v.x*frame.normal.x+v.y*frame.normal.y;
                actualNormal[u][q]=normal;close(normal,canonical[u][q],1e-11*(1+std::abs(normal)),"canonical versus actual unit source");expected+=.7*nodes[q].physicalWeight*(nodes[q].sample.curvature-k)*normal;}}
        close(cap.supportMeasures()[u],weight,1e-12*std::max(double(weight),1e-30),"actual source support measure");
        const double observed=weight>0?cap.curvatureAverages()[u]:0,target=weight>0?k:0;averageError=std::hypot(averageError,observed-target);averageScale=std::hypot(averageScale,target);
        error=std::hypot(error,cap.balanceLoads()[u]-expected);scale=std::hypot(scale,expected);
        if(exactCurvature>=0&&weight>0)close(k,exactCurvature,1e-10*std::max(1.,exactCurvature),"analytic quadratic curvature");
        if(f.negativeCell>=0&&f.positiveCell>=0&&mesh.cells(mesh.cellIndexToLabel[f.negativeCell]).level()!=mesh.cells(mesh.cellIndexToLabel[f.positiveCell]).level())++coarseFaces;
        if(area==0)check(cap.balanceLoads()[u]==0,"axis force");}
    const double e=relative(error,scale),av=relative(averageError,averageScale);check(e<=1e-12&&av<=1e-12,"independent traction assembly");peakLoad=std::max(peakLoad,e);peakMean=std::max(peakMean,av);
    const auto incidence=pressure.apply(mesh,CchNativeValues(l,std::vector<double>(l.velocityCount()),phase.liquidFractions()));
    const Eigen::Map<const Eigen::VectorXd> sigma(incidence.data(),incidence.size());const Eigen::VectorXd factored=cap.stressMatrix()*sigma;
    double factorError=0;for(int u=0;u<l.velocityCount();++u)factorError=std::hypot(factorError,factored[u]-cap.balanceLoads()[u]);check(relative(factorError,scale)<=1e-12,"complete stress matrix pressure factorization");
    if(!allStress.empty())for(double j:{0.,1.,-3.,1e120,1e-280}){
        const auto assembled=Quantum::Hydrodynamics::Detail::assembleCchInterfaceStress(mesh,pressure,phase,source,&law,std::vector<double>(allStress.size(),j));
        check(assembled.matrix.rows()==l.velocityCount()&&assembled.matrix.cols()==l.velocityCount(),"full native constant matrix dimensions");
        for(int u=0;u<l.velocityCount();++u)check(assembled.matrix.coeff(u,u)==j,"full constant stress diagonal");
        for(int outer=0;outer<assembled.matrix.outerSize();++outer)for(CchPressureLoad::Matrix::InnerIterator it(assembled.matrix,outer);it;++it)check(it.row()==it.col()||it.value()==0,"constant stress off-diagonal");++negativeControls;
    }
    for(const auto* profile:{&p1,&p2,&jump}){
        const auto assembled=Quantum::Hydrodynamics::Detail::assembleCchInterfaceStress(mesh,pressure,phase,source,&law,*profile);
        const double base=profile->empty()?0:profile->front();double defect=0,norm=0;
        for(int u=0;u<l.velocityCount();++u){long double positive=0,offset=0;
            for(std::size_t q=0;q<nodes.size();++q){const long double weight=nodes[q].physicalWeight*canonical[u][q]*canonical[u][q];positive+=weight;offset+=weight*((*profile)[selectedPoints[q]]-base);}
            const long double mean=base+(positive>0?offset/positive:0);long double expected=mean*incidence[u];
            for(std::size_t q=0;q<nodes.size();++q)expected-=nodes[q].physicalWeight*((*profile)[selectedPoints[q]]-mean)*actualNormal[u][q];
            defect=std::hypot(defect,assembled.load[u]-expected);norm=std::hypot(norm,expected);}
        check(relative(defect,norm)<=1e-12,"complete P1/P2/jump actual virtual work");++negativeControls;
    }
    std::vector<double> sum(p1.size());for(std::size_t i=0;i<sum.size();++i)sum[i]=p1[i]+p2[i];
    const auto one=Quantum::Hydrodynamics::Detail::assembleCchInterfaceStress(mesh,pressure,phase,source,&law,p1);
    const auto two=Quantum::Hydrodynamics::Detail::assembleCchInterfaceStress(mesh,pressure,phase,source,&law,p2);
    const auto both=Quantum::Hydrodynamics::Detail::assembleCchInterfaceStress(mesh,pressure,phase,source,&law,sum);
    const CchPressureLoad::Matrix linear=both.matrix-one.matrix-two.matrix;
    check(linear.norm()<=1e-12*std::max(1.,both.matrix.norm()),"complete native stress matrix linearity");
    // Adding a spatially constant stress is a pressure jump. It must leave
    // the variable weak force unchanged and retain every prescribed reaction,
    // including native rows without surface support. This checks the full
    // matrix shift, not only its action on this particular C*alpha.
    if(!p2.empty()){
        constexpr double shift=2.75;auto shifted=p2;for(double& j:shifted)j+=shift;
        const auto changed=Quantum::Hydrodynamics::Detail::assembleCchInterfaceStress(mesh,pressure,phase,source,&law,shifted);
        CchPressureLoad::Matrix identity(l.velocityCount(),l.velocityCount());identity.setIdentity();
        const CchPressureLoad::Matrix defect=changed.matrix-two.matrix-shift*identity;
        check(defect.norm()<=1e-12*std::max({1.,changed.matrix.norm(),two.matrix.norm()}),"constant pressure-jump full matrix shift");
        double residual=0,scale=0;for(int u=0;u<l.velocityCount();++u){
            residual=std::hypot(residual,changed.load[u]-two.load[u]-shift*incidence[u]);
            scale=std::hypot(scale,shift*incidence[u]);
            close(changed.constantReactions[u]-two.constantReactions[u],
                std::binary_search(law.velocityDofs().begin(),law.velocityDofs().end(),u)?shift*incidence[u]:0,
                1e-12*(1+std::abs(shift*incidence[u])),"constant pressure-jump prescribed reaction");
        }
        check(residual<=1e-12*std::max(1.,scale),"constant pressure-jump load shift");++negativeControls;
    }
    if(exactCurvature>0){const auto mass=CchMomentumMass::build(mesh,context.geometry,phase,1000,1);
        const auto block=massBlock(mass,.125);const auto saddle=CchSinglePhaseSaddle::buildClosed(pressure,block);auto rhs=cap.balanceLoads();for(auto& v:rhs)v=-v;
        const auto result=saddle.solve(mesh,CchSaddleRhs(l,rhs,std::vector<double>(l.velocityCount())));double sum=0,volume=0,maxU=0,scaleU=0,maxP=0;
        for(int p=0;p<l.pressureCount();++p){sum+=phase.cellVolumes()[p].liquid;volume+=phase.metrics().cellVolumes()[p];}
        const double jump=.7*exactCurvature;for(int p=0;p<l.pressureCount();++p)maxP=std::max(maxP,std::abs(result.values().pressure()[p]-jump*(phase.liquidFractions()[p]-sum/volume)));
        for(int u=0;u<l.velocityCount();++u){maxU=std::max(maxU,std::abs(result.values().normalVelocity()[u]));scaleU=std::max(scaleU,std::abs(cap.balanceLoads()[u])/block.matrix().coeff(u,u));}
        const double ev=relative(maxU,scaleU),ep=maxP/jump;check(ev<=1e-11&&ep<=1e-10,"true Hermite sphere static balance");peakVelocity=std::max(peakVelocity,ev);peakPressure=std::max(peakPressure,ep);
        double pressureAlpha=0;for(int p=0;p<l.pressureCount();++p)pressureAlpha+=phase.metrics().cellVolumes()[p]*result.values().pressure()[p]*(phase.liquidFractions()[p]-sum/volume);
        check(pressureAlpha>0,"liquid interior must have positive pressure jump");++negativeControls;}
    check(phase.publication().state().primary().phaseVolumes.empty(),"capillary introduced inventory");
    std::cout<<"CAPILLARY,"<<scenes<<','<<l.pressureCount()<<','<<l.velocityCount()<<','<<e<<','<<av<<std::endl;
}
void wave(int n,bool amr){
    auto mesh=makeMesh(n,amr?1:0);const double amplitude=1e-4,k=double(2*pi);
    const Profile shape=[=](Vector2D p){return Jet{p.y-.5-amplitude*std::cos(k*p.x),(p.x==0||p.x==1)?0:amplitude*k*std::sin(k*p.x),1,0};};
    const auto phase=CchPhaseGeometry::build(*mesh,metric(*mesh,false),publish(*mesh,shape,false,{{0,0},{0,1}}));
    const auto& l=phase.layout();const auto pressure=CchPressureLoad::build(phase.metrics().divergence());const SourceContext context(*mesh,phase);
    const auto cap=CchInterfaceNormalLoad::capillary(*mesh,pressure,phase,1,context.source,&context.law);
    const auto& geometry=context.geometry;const auto mass=CchMomentumMass::build(*mesh,geometry,phase,1,1);
    const auto saddle=CchSinglePhaseSaddle::buildClosed(pressure,massBlock(mass,1));auto rhs=cap.balanceLoads();for(auto& x:rhs)x=-x;
    const auto result=saddle.solve(*mesh,CchSaddleRhs(l,rhs,std::vector<double>(l.velocityCount())));
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
void provenance(){auto mesh=makeMesh(4);Case s;s.pureFirstJet=true;const auto phase=CchPhaseGeometry::build(*mesh,metric(*mesh,false),initial(*mesh,s));const auto pressure=CchPressureLoad::build(phase.metrics().divergence());const SourceContext context(*mesh,phase);const auto cap=CchInterfaceNormalLoad::capillary(*mesh,pressure,phase,1,context.source,&context.law),copy=cap;
    check(&cap.balanceLoads()==&copy.balanceLoads()&&&cap.geometry().publication().state().geometry().snapshot()==&phase.publication().state().geometry().snapshot(),"immutable shared source");
    const auto reject=[&](auto action,const char* message){rejects(action,message);++guards;};
    for(double sigma:{-1.,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()})reject([&]{CchInterfaceNormalLoad::capillary(*mesh,pressure,phase,sigma,context.source,&context.law);},"surface tension");
    reject([&]{cap.constantJump();},"no prescribed");const auto fixed=CchInterfaceNormalLoad::constant(*mesh,pressure,phase,-1);reject([&]{fixed.curvatureAverages();},"no curvature");
    reject([&]{cap.validateFor(*mesh,pressure,initial(*mesh,s).state().identity(),phase.timeLayer());},"identity");auto time=phase.timeLayer();time.step=1;reject([&]{cap.validateFor(*mesh,pressure,phase.geometryIdentity(),time);},"time");
    const auto other=CchPressureLoad::build(phase.metrics().divergence());reject([&]{cap.validateFor(*mesh,other,phase.geometryIdentity(),phase.timeLayer());},"pressure owner");
    const auto different=CchPressureLoad::build(CchIntegralDivergence::build(*mesh,phase.layout()));reject([&]{CchInterfaceNormalLoad::capillary(*mesh,different,phase,1,context.source,&context.law);},"divergence owner");
    const auto otherSource=CchAdvectorReconstruction::prepare(*mesh,context.geometry);
    reject([&]{cap.validateFor(*mesh,pressure,phase.geometryIdentity(),phase.timeLayer(),otherSource,&context.law);},"source reconstruction owner");
    reject([&]{cap.validateFor(*mesh,pressure,phase.geometryIdentity(),phase.timeLayer(),context.source,nullptr);},"source boundary owner");
    const auto otherLaw=CchAdvectorBoundary::impermeable(*mesh,context.geometry,wallDofs(*mesh,context.geometry));
    reject([&]{cap.validateFor(*mesh,pressure,phase.geometryIdentity(),phase.timeLayer(),context.source,&otherLaw);},"source boundary owner");
    auto target=makeMesh(4,1);reject([&]{cap.validateFor(*target,pressure,phase.geometryIdentity(),phase.timeLayer());},"mesh");
    const auto transfer=transferHermiteGalsFirstJetMesh(*mesh,phase.publication(),*target,HermiteGalsMeshTransferKind::ExactRefinement);inspect(*target,CchPhaseGeometry::build(*target,metric(*target,false),transfer.publication()),1/s.radius);
}
}
int main(){try{std::cout<<std::setprecision(17);
    for(bool rz:{false,true})for(int amr:{0,1,2}){auto mesh=makeMesh(4,amr);Case s;s.pureFirstJet=true;s.radial=rz;
        inspect(*mesh,CchPhaseGeometry::build(*mesh,metric(*mesh,rz),initial(*mesh,s)),(rz?2:1)/s.radius);
        s.nonpolynomial=true;inspect(*mesh,CchPhaseGeometry::build(*mesh,metric(*mesh,rz),initial(*mesh,s)));
        s.planar=true;s.planeOffset=.43;inspect(*mesh,CchPhaseGeometry::build(*mesh,metric(*mesh,rz),initial(*mesh,s)),0);}
    for(double phase:{-1.,1.}){auto mesh=makeMesh(4);const Profile pure=[=](Vector2D){return Jet{phase,0,0,0};};inspect(*mesh,CchPhaseGeometry::build(*mesh,metric(*mesh,false),publish(*mesh,pure,false)));}
    {auto mesh=makeMesh(8);const Profile pair=[](Vector2D p){const double x=p.x-.3,y=p.y-.4,X=p.x-.7,Y=p.y-.6,q=x*x+y*y-.13*.13,Q=X*X+Y*Y-.13*.13;return Jet{q*Q,2*(x*Q+X*q),2*(y*Q+Y*q),4*(x*Y+X*y)};};
        const auto pub=publish(*mesh,pair,false);check(pub.topology().curves.size()==2,"two real interface components");inspect(*mesh,CchPhaseGeometry::build(*mesh,metric(*mesh,false),pub));}
    {auto mesh=makeMesh(8,1);Case s;s.pureFirstJet=true;inspect(*mesh,CchPhaseGeometry::build(*mesh,metric(*mesh,false,true),initial(*mesh,s)),1/s.radius);}
    provenance();for(int n:{8,16})for(bool amr:{false,true})wave(n,amr);
    std::cout<<"H13_SUMMARY,"<<scenes<<','<<guards<<','<<coarseFaces<<','<<negativeControls<<','<<peakLoad<<','<<peakMean<<','<<peakVelocity<<','<<peakPressure<<",PASS\n";return 0;
}catch(const std::exception& e){std::cerr<<"H13 FAIL: "<<e.what()<<'\n';return 1;}}
