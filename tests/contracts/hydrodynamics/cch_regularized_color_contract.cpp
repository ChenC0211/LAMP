#include "Modules/Hydrodynamics/CchRegularizedColor.h"
#include "interface_first_jet_fixture.h"
#include <chrono>

namespace
{
using namespace Quantum::Hydrodynamics;
long double referenceH(long double s,long double e)
{if(s<=-e)return 0;if(s>=e)return 1;return .5L+s/(2*e)+std::sin(pi*s/e)/(2*pi);}
long double primitive(long double s,long double e)
{if(s<=-e)return 0;if(s>=e)return s;return .5L*(s+e)+(s*s-e*e)/(4*e)-e*(std::cos(pi*s/e)+1)/(2*pi*pi);}
CchPhysicalDivergence metrics(const QuadtreeMesh& mesh,bool rz,bool selected=false)
{FaceSegmentGraph graph;graph.build(mesh,[selected](const QuadtreeMesh&,const Cell& c){return !selected||c.centerCoord().x<.5;});
    return CchPhysicalDivergence::build(mesh,CchIntegralDivergence::build(mesh,CchDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,rz))));}
HermiteGalsPublication publication(const QuadtreeMesh& mesh,const Case& shape,double scale=1)
{
    const auto base=shape.profile(0);const Profile f=[=](Vector2D p){auto q=base(p);q.value*=scale;q.dx*=scale;q.dy*=scale;q.dxy*=scale;return q;};
    const auto seed=stateFor(mesh,f,shape.radial,{},[](int,Vector2D,Vector2D){return 0.;},shape.symmetry(mesh),HermiteGalsTransport::FirstJet);
    HermiteGalsCandidate candidate(mesh,{},seed.primary());const auto mixed=closeHermiteGalsMixedAmr(mesh,candidate,shape.symmetry(mesh));
    HermiteGalsPublicationConfig cfg;cfg.volume.symmetryPlanes=shape.symmetry(mesh);const auto id=candidate.identity();
    return publishHermiteGalsFirstJetGeometry(mesh,std::move(candidate),mixed,{id,centers(mesh,f)},cfg);
}
HermiteGalsPublication profilePublication(const QuadtreeMesh& mesh,const Profile& original,double scale=1)
{
    const Profile f=[=](Vector2D p){auto q=original(p);q.value*=scale;q.dx*=scale;q.dy*=scale;q.dxy*=scale;return q;};
    const std::vector<LevelSetInterfacePanels::SymmetryPlane> symmetry{{1,mesh.minCoord.y},{1,mesh.maxCoord.y}};
    const auto seed=stateFor(mesh,f,false,{},[](int,Vector2D,Vector2D){return 0.;},symmetry,HermiteGalsTransport::FirstJet);
    HermiteGalsCandidate candidate(mesh,{},seed.primary());const auto mixed=closeHermiteGalsMixedAmr(mesh,candidate,symmetry);
    HermiteGalsPublicationConfig cfg;cfg.volume.symmetryPlanes=symmetry;
    const auto id=candidate.identity();return publishHermiteGalsFirstJetGeometry(mesh,std::move(candidate),mixed,{id,centers(mesh,f)},cfg);
}
void outputCell(const char* name,bool rz,int amr,double parameter,const CchRegularizedColor& color,int p)
{
    const auto& r=color.cellReports()[p];
    const auto cap=color.phaseGeometry().metrics().cellVolumes()[p];
    check(r.transitionVolumeLower>=0&&r.transitionVolumeLower<=r.transitionVolumeEstimate+1e-14*cap&&
        r.transitionVolumeEstimate<=r.transitionVolumeUpper+1e-14*cap&&r.transitionVolumeUpper<=cap*(1+1e-14),"invalid force-band volume enclosure");
    check(r.queries+r.enclosureQueries<=color.config().maxQueriesPerCell,"original-source work limit exceeded");
    check(r.minimumCertifiedGradientRatio>color.config().relativeGradientResolution,"nonpure accepted panel lacked gradient certificate");
    std::cout<<"BF01_COLOR_CELL,"<<name<<','<<rz<<','<<amr<<','<<parameter<<','<<color.config().epsilon<<','<<p<<','
        <<color.colors()[p]<<','<<color.phaseGeometry().liquidFractions()[p]<<','<<r.transitionVolumeEstimate<<','
        <<r.transitionVolumeLower<<','<<r.transitionVolumeUpper<<','<<r.gradientScale<<','<<r.minimumResolvedGradientRatio<<','<<r.minimumCertifiedGradientRatio<<','
        <<r.absoluteErrorEstimate<<','<<r.queries<<','<<r.enclosureQueries<<','<<r.acceptedPanels<<','<<r.deepestPanel<<','<<r.certifiedPure<<'\n';
}
template<class F> void rejected(F f,const char* name="input"){bool caught=false;try{f();}catch(const std::exception&){caught=true;}
    check(caught,std::string("invalid color input accepted: ")+name);}
void kernel()
{
    constexpr double e=.17;double last=-1;
    for(int k=0;k<=400;++k){const double s=e*(-1+2.*k/400),v=CchRegularizedColor::heaviside(s,e);
        check(v>=0&&v<=1&&v>=last,"kernel boundedness/monotonicity");last=v;
        close(v,referenceH(s,e),3e-15,"kernel independent formula");
        close(v+CchRegularizedColor::heaviside(-s,e),1,3e-15,"kernel symmetry");
        if(k>4&&k<396){const double h=1e-7*e;
            close((CchRegularizedColor::heaviside(s+h,e)-CchRegularizedColor::heaviside(s-h,e))/(2*h),
                CchRegularizedColor::delta(s,e),2e-8/e,"kernel derivative");}}
    for(double s:{-e,std::nextafter(-e,0.),std::nextafter(e,0.),e})
        check(CchRegularizedColor::heaviside(s,e)>=0&&CchRegularizedColor::heaviside(s,e)<=1,"endpoint kernel cancellation");
    check(CchRegularizedColor::heaviside(-e,e)==0&&CchRegularizedColor::heaviside(e,e)==1
        &&CchRegularizedColor::delta(e,e)==0,"exact kernel endpoints");
    rejected([]{CchRegularizedColor::heaviside(0,0);});
    long double mass=0;for(int k=0;k<4096;++k)mass+=CchRegularizedColor::delta(e*(-1+2*(k+.5)/4096),e)*(2*e/4096);
    close(mass,1,2e-14,"delta positive kernel normalization");
    const double huge=std::numeric_limits<double>::max();
    const double hugeDelta=CchRegularizedColor::delta(huge/2,huge);
    check(std::isfinite(hugeDelta)&&hugeDelta>0,"finite huge-width delta lost positive representable value");
    close(static_cast<long double>(hugeDelta)*huge,.5L,2e-14,"huge-width delta normalized scale");
    rejected([]{CchRegularizedColor::delta(0,std::numeric_limits<double>::denorm_min());},"nonfinite tiny-width delta");
}
void planes()
{
    int cells=0;long double peak=0;
    for(bool rz:{false,true})for(int amr:{0,1})for(double offset:{.375,.43,.49}){
        auto mesh=makeMesh(4,amr);Case shape;shape.pureFirstJet=true;shape.planar=true;shape.radial=rz;shape.planeOffset=offset;
        const auto met=metrics(*mesh,rz);const auto phase=CchPhaseGeometry::build(*mesh,met,publication(*mesh,shape));
        const auto oldFractions=phase.liquidFractions();CchRegularizedColorConfig cfg;cfg.epsilon=amr?.1875:.375;
        const auto color=CchRegularizedColor::build(*mesh,phase,cfg);const auto repeat=CchRegularizedColor::build(*mesh,phase,cfg);
        check(color.colors()==repeat.colors(),"repeat color construction not deterministic");
        check(phase.liquidFractions()==oldFractions&&phase.publication().state().primary().phaseVolumes.empty(),"color changed sharp material/inventory");
        color.validateFor(*mesh,phase);
        rejected([&]{color.validateFor(*mesh,color.geometryIdentity(),{1,0,.01,Position::Midpoint});},"clock");
        for(int p=0;p<phase.layout().pressureCount();++p){const auto bounds=phase.publication().state().geometry().snapshot().nativeCellBounds(phase.layout().nativeCell(p));
            const long double a=shape.axial()?bounds[0].y:bounds[0].x,b=shape.axial()?bounds[1].y:bounds[1].x;
            const auto exact=(primitive(offset-a,cfg.epsilon)-primitive(offset-b,cfg.epsilon))/(b-a);
            peak=std::max(peak,close(color.colors()[p],exact,4e-9,"physical plane cell average"));++cells;
            outputCell("plane",rz,amr,offset,color,p);
            if(b<=offset-cfg.epsilon)check(color.colors()[p]==1,"pure liquid color not exact");
            if(a>=offset+cfg.epsilon)check(color.colors()[p]==0,"pure gas color not exact");}
        const auto scaled=CchPhaseGeometry::build(*mesh,met,publication(*mesh,shape,7));
        const auto scaleColor=CchRegularizedColor::build(*mesh,scaled,cfg);
        for(int p=0;p<phase.layout().pressureCount();++p)close(color.colors()[p],scaleColor.colors()[p],5e-9,"phi scaling changed length proxy");
        rejected([&]{color.validateFor(*mesh,scaled);},"geometry");
        if(offset==.375)for(double ratio:{1.,2.}){auto alternative=cfg;alternative.epsilon=(amr?.125:.25)*ratio;
            const auto band=CchRegularizedColor::build(*mesh,phase,alternative);
            for(int p=0;p<phase.layout().pressureCount();++p){const auto box=phase.publication().state().geometry().snapshot().nativeCellBounds(phase.layout().nativeCell(p));
                const long double a=shape.axial()?box[0].y:box[0].x,b=shape.axial()?box[1].y:box[1].x;
                close(band.colors()[p],(primitive(offset-a,alternative.epsilon)-primitive(offset-b,alternative.epsilon))/(b-a),4e-9,"epsilon scan physical average");
                outputCell("plane_epsilon_scan",rz,amr,offset,band,p);}
            std::cout<<"BF01_EPSILON,"<<rz<<','<<amr<<','<<ratio<<','<<band.report().colorMinusSharpL1<<",PASS\n";}
        cfg.relativeGradientResolution=.99;rejected([&]{CchRegularizedColor::build(*mesh,phase,cfg);},"gradient resolution");
        cfg.relativeGradientResolution=CchRegularizedColorConfig{}.relativeGradientResolution;
        cfg.epsilon=.001;rejected([&]{CchRegularizedColor::build(*mesh,phase,cfg);},"underresolved epsilon");
    }
    std::cout<<"BF01_PLANE,"<<cells<<','<<peak<<",PASS\n";
}
// Analytic Q4 scalar/gradient in one physical dimension; the independent
// midpoint reference never queries the production polynomial or GL nodes.
void quartic()
{
    const auto analytic=[](long double x){const auto t=x-.43L;
        return std::array<long double,2>{{t+.35L*t*t+.2L*t*t*t+.1L*t*t*t*t,
            1+.7L*t+.6L*t*t+.4L*t*t*t}};};
    const Profile f=[=](Vector2D p){const auto q=analytic(p.x);return Jet{static_cast<double>(q[0]),static_cast<double>(q[1]),0,0};};
    const auto average=[&](long double a,long double b,double epsilon,int n){long double sum=0;
        for(int k=0;k<n;++k){const auto q=analytic(a+(b-a)*(k+.5L)/n);
            check(q[1]>0,"quartic reference not monotone");sum+=referenceH(-q[0]/q[1],epsilon);}return sum/n;};
    long double peak=0;int cells=0;
    for(int amr:{0,1}){auto mesh=makeMesh(4,amr);const auto met=metrics(*mesh,false);
        const auto phase=CchPhaseGeometry::build(*mesh,met,profilePublication(*mesh,f));
        const auto scaled=CchPhaseGeometry::build(*mesh,met,profilePublication(*mesh,f,3.7));
        CchRegularizedColorConfig cfg;cfg.epsilon=amr?.1875:.375;
        const auto color=CchRegularizedColor::build(*mesh,phase,cfg),other=CchRegularizedColor::build(*mesh,scaled,cfg);
        for(int p=0;p<phase.layout().pressureCount();++p){const auto b=phase.publication().state().geometry().snapshot().nativeCellBounds(phase.layout().nativeCell(p));
            const auto coarse=average(b[0].x,b[1].x,cfg.epsilon,1024),fine=average(b[0].x,b[1].x,cfg.epsilon,2048);
            peak=std::max(peak,close(color.colors()[p],(4*fine-coarse)/3,4e-9,"independent true-Q4 distance-proxy average"));
            close(color.colors()[p],other.colors()[p],5e-9,"true-Q4 positive scaling changed proxy");
            outputCell("quartic",false,amr,.43,color,p);++cells;}
    }
    std::cout<<"BF01_QUARTIC_REFERENCE,"<<cells<<','<<peak<<",PASS\n";
}
void unresolvedFaceGradient()
{
    auto mesh=makeMesh(4);constexpr double eta=1e-5;
    const Profile f=[=](Vector2D p){const double t=p.x-.5;return Jet{t*t*t+eta*t,3*t*t+eta,0,0};};
    const auto phase=CchPhaseGeometry::build(*mesh,metrics(*mesh,false),profilePublication(*mesh,f));
    // Source geometry remains resolved. The explicit force-color resolution
    // exceeds the minimum face gradient but is below all initial GL samples.
    CchRegularizedColorConfig cfg;cfg.epsilon=.375;cfg.relativeGradientResolution=2e-4;
    cfg.maxDepth=6;cfg.maxQueriesPerCell=4096;
    bool caught=false;
    try{CchRegularizedColor::build(*mesh,phase,cfg);}catch(const std::exception& e){const std::string reason=e.what();
        check(reason.find("gradient")!=std::string::npos &&
            (reason.find("unresolved")!=std::string::npos||reason.find("degenerate")!=std::string::npos),
            "face minimum gradient rejected for another reason: "+reason);caught=true;
        std::cout<<"BF01_GRADIENT_FACE_REJECT,"<<eta<<','<<reason<<'\n';}
    check(caught,"underresolved band gradient on native face silently accepted");
}
// Independent composite midpoint/Richardson reference; no production Gauss nodes
// or geometry/phase cubature are reused for the polynomial circle/sphere.
long double referenceVolume(Vector2D lo,Vector2D hi,const Case& shape,double e,int n)
{
    long double sum=0;const auto f=shape.profile(0);const long double hx=(hi.x-lo.x)/n,hy=(hi.y-lo.y)/n;
    for(int i=0;i<n;++i)for(int j=0;j<n;++j){const Vector2D p{static_cast<double>(lo.x+(i+.5L)*hx),static_cast<double>(lo.y+(j+.5L)*hy)};
        const auto q=f(p);const long double g=std::hypot(q.dx,q.dy);
        const auto c=g==0?(q.value<0?1.L:0.L):referenceH(-q.value/g,e);
        sum+=c*(shape.radial?2*pi*p.x:1);}
    return hx*hy*sum;
}
void curved()
{
    int cells=0;long double peak=0;
    for(bool rz:{false,true})for(int amr:{0,1}){
        auto mesh=makeMesh(4,amr);Case shape;shape.pureFirstJet=true;shape.radial=rz;shape.radius=rz?.34:.19;
        const auto phase=CchPhaseGeometry::build(*mesh,metrics(*mesh,rz),publication(*mesh,shape));
        CchRegularizedColorConfig cfg;cfg.epsilon=amr?.1875:.375;
        for(double epsilon:{cfg.epsilon,amr?.25:.5}){
        cfg.epsilon=epsilon;cfg.relativeTolerance=1e-9;
        const auto color=CchRegularizedColor::build(*mesh,phase,cfg);cfg.relativeTolerance=1e-10;
        const auto finer=CchRegularizedColor::build(*mesh,phase,cfg);
        if(rz){const Vector2D axis{0,shape.cy+.12};const int cell=phase.publication().state().geometry().snapshot().nativeCellCandidates(axis,axis).front();
            const int p=phase.layout().pressureDof(cell);const auto f=shape.profile(0);const auto v=f(axis);
            close(color.pointColor(p,axis),referenceH(-v.value/std::hypot(v.dx,v.dy),color.config().epsilon),3e-14,"RZ axis proxy limit");
            const double h=1e-5;const auto near=color.pointColor(p,{h,axis.y});
            check(std::abs(near-color.pointColor(p,axis))<1e-7,"RZ axis color lost even parity");}
        for(int p=0;p<phase.layout().pressureCount();++p){const auto b=phase.publication().state().geometry().snapshot().nativeCellBounds(phase.layout().nativeCell(p));
            const long double r0=referenceVolume(b[0],b[1],shape,cfg.epsilon,128),r1=referenceVolume(b[0],b[1],shape,cfg.epsilon,256);
            const long double exact=(4*r1-r0)/3/phase.metrics().cellVolumes()[p];
            peak=std::max(peak,close(color.colors()[p],exact,2e-7,"independent curved physical volume"));
            close(color.colors()[p],finer.colors()[p],4e-9,"quadrature refinement changed color");
            check(color.colors()[p]>=0&&color.colors()[p]<=1,"curved color outside bounds");++cells;}
        for(int p=0;p<phase.layout().pressureCount();++p)outputCell(rz?"sphere":"circle",rz,amr,shape.radius,color,p);
        std::cout<<"BF01_CURVE,"<<rz<<','<<amr<<','<<cfg.epsilon<<','<<color.report().sourceQueries<<','<<color.report().absoluteErrorEstimate<<",PASS\n";
        }
    }
    std::cout<<"BF01_CURVED_REFERENCE,"<<cells<<','<<peak<<",PASS\n";
}
void selected()
{
    auto mesh=makeMesh(4,1);Case shape;shape.pureFirstJet=true;shape.planar=true;shape.planeOffset=.43;
    const auto pub=publication(*mesh,shape);const auto full=CchPhaseGeometry::build(*mesh,metrics(*mesh,false),pub);
    const auto part=CchPhaseGeometry::build(*mesh,metrics(*mesh,false,true),pub);
    CchRegularizedColorConfig cfg;cfg.epsilon=.1875;
    const auto a=CchRegularizedColor::build(*mesh,full,cfg),b=CchRegularizedColor::build(*mesh,part,cfg);
    for(int p=0;p<part.layout().pressureCount();++p)check(b.colors()[p]==a.colors()[full.layout().pressureDof(part.layout().nativeCell(p))],"selected color borrowed or changed outside domain");
    rejected([&]{a.validateFor(*mesh,part);});
    cfg.epsilon=-1;rejected([&]{CchRegularizedColor::build(*mesh,full,cfg);});
}
void pure()
{
    for(bool rz:{false,true})for(int sign:{-1,1}){
        auto mesh=makeMesh(4);const Profile f=[sign](Vector2D){return Jet{double(sign),0,0,0};};
        const auto seed=stateFor(*mesh,f,rz,{},[](int,Vector2D,Vector2D){return 0.;},{},HermiteGalsTransport::FirstJet);
        HermiteGalsCandidate candidate(*mesh,{},seed.primary());const auto mixed=closeHermiteGalsMixedAmr(*mesh,candidate);
        const auto id=candidate.identity();const auto pub=publishHermiteGalsFirstJetGeometry(*mesh,std::move(candidate),mixed,{id,centers(*mesh,f)});
        const auto phase=CchPhaseGeometry::build(*mesh,metrics(*mesh,rz),pub);CchRegularizedColorConfig cfg;cfg.epsilon=.375;
        const auto color=CchRegularizedColor::build(*mesh,phase,cfg);
        for(double c:color.colors())check(c==(sign<0?1:0),"zero-gradient pure constant was not exact");
        check(color.report().sourceQueries==0&&color.report().pureCells==color.colors().size(),"pure cells sampled or unsupported");
    }
}
}
int main(){try{std::cout<<std::setprecision(17);
    std::cout<<"BF01_COLOR_CELL_HEADER,shape,rz,amr,parameter,epsilon,pressure,color,sharp_alpha,band_estimate,band_lower,band_upper,gradient_scale,min_gradient_ratio,min_certified_gradient_ratio,error_estimate,queries,enclosure_queries,accepted_panels,deepest_panel,certified_pure\n";
    kernel();planes();quartic();unresolvedFaceGradient();curved();selected();pure();std::cout<<"BF01_COLOR_CONTRACT,PASS\n";return 0;}
catch(const std::exception& e){std::cerr<<"BF01_COLOR_FAIL,"<<e.what()<<'\n';return 1;}}
