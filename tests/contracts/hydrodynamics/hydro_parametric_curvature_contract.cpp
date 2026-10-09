#include "Modules/Interface/Compatibility/InterfaceParametricMotion.h"
#include "Modules/Interface/Compatibility/InterfacePrescribedData.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>

namespace {
using namespace Quantum::Geometry;
using Quantum::Geometry::Parametric::Point;
using Quantum::Geometry::Parametric::SurfaceSample;
using Quantum::Geometry::Parametric::PolynomialArc;
using Quantum::Geometry::Parametric::PolynomialCurve;
using Quantum::Geometry::Parametric::gaussEight;
using Quantum::Geometry::Parametric::sampleSurface;
using Quantum::Geometry::Parametric::sampleSurfaceRZ;
using namespace InterfaceParametricTransport;
constexpr double dt=.05;
const double pi=std::acos(-1.);
void require(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
struct Reference {
    int kind;
    bool rz()const{return kind==2||kind==4;}
    bool ellipse()const{return kind==1||kind==3;}
    bool nonlinear()const{return kind>=3;}
    const char* name()const{return kind==0?"circleXY":kind==1?"shearEllipseXY":kind==2?"sphereRZ":kind==3?"quadraticEllipseXY":"quadraticSphereRZ";}
    double radius()const{return kind==0?.18:ellipse()?.22:.2;}
    double volume()const{return rz()?4*pi*.2*.2*.2/3:pi*radius()*(ellipse()?.14:.18);}
    MapJet map(Vector2D p,double t,double sign)const {
        const Vector2D v=rz()?Vector2D{0,sign*(kind==4?6*p.x*p.x:.5)}:
            Vector2D{sign*(kind==0?.5:kind==1?2*(p.y-.5):12*((p.y-.5)*(p.y-.5)-.01)),0};
        const auto move=v*t;
        return {p+move,{1,kind==4?12*sign*p.x*t:0},{kind==1?2*sign*t:kind==3?24*sign*(p.y-.5)*t:0,1},v,move};
    }
    SurfaceSample sample(double s,int step)const {
        const double q=dt*std::min(step,8-step),theta=(rz()?pi:2*pi)*s,w=rz()?pi:2*pi;
        Point x,d,dd;
        if(rz()){
            x={.2*std::sin(theta),.35-.2*std::cos(theta)+.5*q};
            d={.2*w*std::cos(theta),.2*w*std::sin(theta)};
            dd={-.2*w*w*std::sin(theta),.2*w*w*std::cos(theta)};
            if(kind==4){x.y-=.5*q;x.y+=6*q*x.x*x.x;dd.y+=12*q*(d.x*d.x+x.x*dd.x);d.y+=12*q*x.x*d.x;}
        } else {
            const double a=radius(),b=ellipse()?.14:.18;
            x={.4+a*std::cos(theta),.5+b*std::sin(theta)};
            d={-a*w*std::sin(theta),b*w*std::cos(theta)};
            dd={-a*w*w*std::cos(theta),-b*w*w*std::sin(theta)};
            if(kind==0)x.x+=.5*q;
            else if(kind==1){x.x+=2*q*(x.y-.5);d.x+=2*q*d.y;dd.x+=2*q*dd.y;}
            else{x.x+=12*q*((x.y-.5)*(x.y-.5)-.01);d.x+=24*q*(x.y-.5)*d.y;dd.x+=24*q*(d.y*d.y+(x.y-.5)*dd.y);}
        }
        const double speed=norm(d),meridian=cross(d,dd)/(speed*speed*speed);
        const double curvature=rz()?((s==0||s==1)?2*meridian:meridian+d.y/(speed*x.x)):meridian;
        return {x,d,{d.y/speed,-d.x/speed},speed,curvature};
    }
};
// The accepted curve uses the production curvature entry point. The transient
// degree-eight source has no persistent curvature field; differentiate it here
// solely to attribute the error introduced by constrained reconstruction.
template<std::size_t N> SurfaceSample surface(const PolynomialArc<N>& arc,double t,bool rz,bool pole,
        const Quantum::Geometry::Parametric::CurveBoundary& boundary){
    if constexpr(N==5)return rz?sampleSurfaceRZ(arc,t,boundary):sampleSurface(arc,t);
    else {
        const auto p=arc.sampleFraction(t),d=arc.sampleFraction(t,1),dd=arc.sampleFraction(t,2);
        const double speed=norm(d),meridian=cross(d,dd)/(speed*speed*speed);
        const Point normal{d.y/speed,-d.x/speed};
        return {p,d,normal,speed,rz?(pole?2*meridian:meridian+normal.x/p.x):meridian};
    }
}
struct Error {
    double position=0,normal=0,curvature=0,maximum=0,axis=0,join=0,positionJoin=0;
};
template<std::size_t N> Error measure(const PolynomialCurve<N>& curve,const Reference& ref,int step){
    Error out;long double area=0,p2=0,n2=0,k2=0,kr2=0;
    const auto& g=gaussEight();
    for(const auto& a:curve.arcs){
        for(int k=0;k<8;k++){
            const double t=.5*(1+g.x[k]),s=a.a+(a.b-a.a)*t;
            const auto v=surface(a,t,ref.rz(),false,curve.boundary),e=ref.sample(s,step);
            const double weight=.5*g.w[k]*(a.b-a.a)*e.speed*(ref.rz()?2*pi*e.point.x:1);
            area+=weight;p2+=weight*dot(v.point-e.point,v.point-e.point);
            n2+=weight*dot(v.normal-e.normal,v.normal-e.normal);
            k2+=weight*(v.curvature-e.curvature)*(v.curvature-e.curvature);kr2+=weight*e.curvature*e.curvature;
        }
        for(int end=0;end<2;end++){
            const double s=end?a.b:a.a;const bool pole=ref.rz()&&(s==0||s==1);
            const double v=surface(a,double(end),ref.rz(),pole,curve.boundary).curvature,e=ref.sample(s,step).curvature;
            const double error=std::abs(v-e)/std::abs(e);out.maximum=std::max(out.maximum,error);
            if(pole)out.axis=std::max(out.axis,error);
        }
    }
    for(std::size_t i=1;i<curve.arcs.size()+(ref.rz()?0:1);i++){
        const auto& a=curve.arcs[i-1];const auto& b=curve.arcs[i%curve.arcs.size()];
        const auto u=surface(a,1,ref.rz(),false,curve.boundary),v=surface(b,0,ref.rz(),false,curve.boundary);
        out.join=std::max(out.join,std::abs(u.curvature-v.curvature)*ref.radius());
        out.positionJoin=std::max(out.positionJoin,norm(u.point-v.point));
    }
    out.position=std::sqrt(double(p2/area))/ref.radius();out.normal=std::sqrt(double(n2/area));out.curvature=std::sqrt(double(k2/kr2));
    require(std::isfinite(out.position+out.normal+out.curvature+out.maximum+out.axis+out.join),"nonfinite analytic curvature comparison");return out;
}
InterfaceNormalFluxView exactFlux(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const Reference& ref,double sign){
    InterfaceNormalFluxView::Queries q;q.retainedMomentDegree=2;
    q.chart=[&mesh,&domain](int row){const auto f=domain.face(mesh,row);const auto& e=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
        const auto a=mesh.verts(e.edgeToVert[0]).coord(),b=mesh.verts(e.edgeToVert[1]).coord();return InterfaceNormalFluxView::FaceChart{a,(a+b)*.5,b};};
    q.polynomial=[chart=q.chart,&mesh,&domain,ref,sign](int row,int half){const auto f=domain.face(mesh,row);const auto c=chart(row);
        const auto a=half?c.midpoint:c.begin,b=half?c.end:c.midpoint;UnitPolynomialTerm term;
        const double width=f.horizontal?b.x-a.x:b.y-a.y;
        if(ref.kind==4){const double r=a.x,h=b.x-a.x;term.factor=f.horizontal?width*12*pi*sign:0;
            term.coefficients[0]=r*r*r;term.coefficients[1]=3*r*r*h;term.coefficients[2]=3*r*h*h;term.coefficients[3]=h*h*h;}
        else if(ref.rz()){term.factor=width*(f.horizontal?.5*sign:0);term.coefficients[0]=2*pi*a.x;term.coefficients[1]=2*pi*(b.x-a.x);}
        else if(ref.kind==3){const double y=a.y-.5,h=b.y-a.y;term.factor=f.horizontal?0:width*12*sign;
            term.coefficients[0]=y*y-.01;term.coefficients[1]=2*y*h;term.coefficients[2]=h*h;}
        else {term.factor=f.horizontal?0:width*sign;term.coefficients[0]=ref.kind==0?.5:2*(a.y-.5);term.coefficients[1]=ref.kind==0?0:2*(b.y-a.y);}
        return std::vector<UnitPolynomialTerm>{term};};
    q.integral=[p=q.polynomial](int r,int h,double a,double b,double s,double t){return integrateUnitPolynomial(p(r,h),a,b,s,t);};
    q.density=[p=q.polynomial](int r,int h,bool end,long double s){return evaluateUnitPolynomial(p(r,h),double(end?1-s:s));};
    q.phaseIntegral=[](int,int,const LevelSetMappedSegment&)->long double{throw std::logic_error("unused scalar phase query in analytic curvature study");};
    return InterfaceNormalFluxView(std::move(q));
}
struct Output {
    std::ofstream metrics,curves,projection;
    Output(){if(const char* prefix=std::getenv("QUANTUM_PARAMETRIC_CURVATURE_OUTPUT")){
        metrics.open(std::string(prefix)+"_metrics.csv");curves.open(std::string(prefix)+"_curves.csv");projection.open(std::string(prefix)+"_projection.csv");
        require(metrics&&curves&&projection,"cannot open analytic curvature output");
        metrics<<std::setprecision(17)<<"case,N,segments,step,time,stage,arcs,cells,position_l2_over_R,normal_l2,curvature_relative_l2,curvature_endpoint_relative_max,axis_relative_max,join_jump_times_R,position_join,inventory_drift_relative,initial_volume_error_relative,fraction_error\n";
        curves<<std::setprecision(17)<<"case,N,step,stage,s,x,y,kappa,reference_x,reference_y,reference_kappa\n";
        projection<<std::setprecision(17)<<"case,N,step,source_arcs,accepted_arcs,relative_h2,curvature_source_to_accepted_l2,position_source_to_accepted_l2_over_R\n";
    }}
    template<std::size_t N> Error record(const PolynomialCurve<N>& curve,const Reference& ref,int n,int step,const char* stage,
            std::size_t cells,double drift,double volumeError,double fraction){
        const auto e=measure(curve,ref,step);
        if(metrics)metrics<<ref.name()<<','<<n<<','<<2*n<<','<<step<<','<<step*dt<<','<<stage<<','<<curve.arcs.size()<<','<<cells<<','<<e.position<<','<<e.normal<<','<<e.curvature<<','<<e.maximum<<','<<e.axis<<','<<e.join<<','<<e.positionJoin<<','<<drift<<','<<volumeError<<','<<fraction<<std::endl;
        if(curves&&(step==0||step==4||step==8))for(int k=0;k<=512;k++){
            const double s=k/512.;const auto& a=curve.arcs[curve.index(s)];const auto v=surface(a,(s-a.a)/(a.b-a.a),ref.rz(),ref.rz()&&(k==0||k==512),curve.boundary),r=ref.sample(s,step);
            curves<<ref.name()<<','<<n<<','<<step<<','<<stage<<','<<s<<','<<v.point.x<<','<<v.point.y<<','<<v.curvature<<','<<r.point.x<<','<<r.point.y<<','<<r.curvature<<'\n';}
        return e;
    }
    void compare(const InterfaceParametricMotion::Result& next,const Reference& ref,int n,int step){
        const auto& source=next.transportSources().front();const auto& accepted=next.state.parametricGeometry();const auto& curve=accepted.geometry();
        std::vector<double> knots{1};for(const auto& a:source.arcs)knots.push_back(a.a);for(const auto& a:curve.arcs)knots.push_back(a.a);
        std::sort(knots.begin(),knots.end());knots.erase(std::unique(knots.begin(),knots.end()),knots.end());
        long double area=0,p2=0,k2=0,kr2=0;const auto& g=gaussEight();
        for(std::size_t i=1;i<knots.size();i++)for(int k=0;k<8;k++){
            const double h=knots[i]-knots[i-1],s=knots[i-1]+.5*h*(1+g.x[k]);const auto& a=source.arcs[source.index(s)];const auto& b=curve.arcs[curve.index(s)];
            const auto u=surface(a,(s-a.a)/(a.b-a.a),ref.rz(),false,source.boundary),v=surface(b,(s-b.a)/(b.b-b.a),ref.rz(),false,curve.boundary),r=ref.sample(s,step);
            const double weight=.5*h*g.w[k]*r.speed*(ref.rz()?2*pi*r.point.x:1);area+=weight;
            p2+=weight*dot(v.point-u.point,v.point-u.point);k2+=weight*(v.curvature-u.curvature)*(v.curvature-u.curvature);kr2+=weight*r.curvature*r.curvature;
        }
        if(projection)projection<<ref.name()<<','<<n<<','<<step<<','<<source.arcs.size()<<','<<curve.arcs.size()<<','<<accepted.report().relativeH2<<','<<std::sqrt(double(k2/kr2))<<','<<std::sqrt(double(p2/area))/ref.radius()<<std::endl;
    }
};
struct Summary {double curvature=0,source=0,drift=0,transfer=0;};
Summary study(const Reference& ref,int n,Output& out){
    DictionaryMesh input;input.minCoord={0,0};input.maxCoord={1,1};input.Nx=input.Ny=n;input.maxCellNum=16384;input.maxInitialIterations=0;
    QuadtreeMesh base(input);FaceSegmentGraph graph;graph.build(base);const auto baseDomain=NativeLeafDomain::build(base,graph,ref.rz());
    auto fine=base.clone();DictionaryMesh::RefineBox box;box.min=ref.rz()?Vector2D{0,.5}:Vector2D{.5,0};box.max={1,1};box.level=1;fine->applyRefineBoxes({box});
    FaceSegmentGraph fineGraph;fineGraph.build(*fine);const auto fineDomain=NativeLeafDomain::build(*fine,fineGraph,ref.rz());
    require(!fineDomain.uniformLeaves(),"curvature AMR study requires mixed leaves");
    const QuadtreeMesh* mesh=&base;const NativeLeafDomain* domain=&baseDomain;
    const auto initial=ref.rz()?InterfacePrescribedData::parametricSpheroid(.35,{.2,.2},2*n):InterfacePrescribedData::parametricEllipse({.4,.5},{ref.radius(),ref.ellipse()?.14:.18},2*n);
    auto state=InterfaceParametricMotion::initialize(*mesh,*domain,initial);
    const auto volume=[](const InterfaceState& s){const auto& v=s.parametricGeometry().inventory();return std::accumulate(v.begin(),v.end(),0.);};
    const double initialVolume=volume(state),volumeError=(initialVolume-ref.volume())/ref.volume();Summary summary;
    const auto record=[&](int step,const char* stage){const auto& g=state.parametricGeometry();const double drift=(volume(state)-initialVolume)/initialVolume;
        const auto e=out.record(g.geometry(),ref,n,step,stage,g.grid().size(),drift,volumeError,g.report().maximumInventoryFractionError);
        summary.curvature=std::max(summary.curvature,e.curvature);summary.drift=std::max(summary.drift,std::abs(drift));return e;};
    record(0,"initial");InterfaceParametricProjection::Options options;options.materialSegments=2*n;options.relativeH2=.001*std::pow(8./n,2);
    for(int step=1;step<=8;step++){
        const double sign=step<=4?1.:-1.;const auto flux=exactFlux(*mesh,*domain,ref,sign);
        const Flow flow=[=](Vector2D p,double t){return ref.map(p,t,sign);};
        const double begin=state.progress().endTime;const InterfaceStepTime time{step,begin,begin+dt,dt};
        auto next=InterfaceParametricMotion::prepareTimeLinear(*mesh,*domain,state,time,flow,flux,flux,options);
        require(next.transportSources().size()==1,"unexpected transient component count");
        if(ref.nonlinear())out.compare(next,ref,n,step);
        const auto e=out.record(next.transportSources().front(),ref,n,step,"source",domain->cellIndexByRow().size(),
            (volume(next.state)-initialVolume)/initialVolume,volumeError,next.state.parametricGeometry().report().maximumInventoryFractionError);
        summary.source=std::max(summary.source,e.curvature);state=std::move(next.state);record(step,"accepted");
        if(step==4||step==8){
            const auto* targetMesh=step==4?fine.get():&base;const auto* targetDomain=step==4?&fineDomain:&baseDomain;
            const auto handle=state.parametricGeometry().boundaryHandle();const auto oldTime=state.progress().endTime;
            auto transfer=InterfaceParametricMotion::prepareMeshTransfer(*mesh,*domain,state,*targetMesh,*targetDomain,options);
            require(transfer.state.parametricGeometry().boundaryHandle()==handle&&transfer.state.progress().endTime==oldTime,"mesh transfer changed curve ownership or physical clock");
            summary.transfer=std::max(summary.transfer,std::abs(transfer.report.inventoryBalance)/initialVolume);
            state=std::move(transfer.state);mesh=targetMesh;domain=targetDomain;record(step,step==4?"refine":"coarsen");
        }
    }
    std::cout<<std::setprecision(17)<<"PARAMETRIC_CURVATURE case="<<ref.name()<<" N="<<n<<" steps=8 time="<<state.progress().endTime
        <<" curvature_l2_max="<<summary.curvature<<" source_l2_max="<<summary.source<<" inventory_drift="<<summary.drift<<" transfer_balance="<<summary.transfer<<std::endl;
    require(summary.drift<1e-10&&summary.transfer<1e-10,"analytic curvature study lost conservative inventory");return summary;
}
}
void hydroParametricPoleContract(){
    namespace P=Quantum::Geometry::Parametric;
    // Exact binary64 coefficients saved at the RZ remesh rejection. Direct
    // evaluation gives zero radius at the still-interior sample below.
    P::Arc right{.9375,1,{{{.042922262749701692,.66684065552910965},
        {-.042365830794543112,.006514147699368162},{-.00082925096870428748,-.0032034699218040203},
        {.00026896706918988722,-.000047156884837961765},{.0000038519443558028123,.0000085656996884409549}}}};
    const auto boundary=P::CurveBoundary::axisMeridian();
    const auto recorded=P::sampleSurfaceRZ(right,.99999999999999967,boundary);
    require(recorded.point.x>0&&std::isfinite(recorded.curvature),"saved near-pole sample lost its positive radius or finite curvature");
    P::Arc left;left.a=0;left.b=.0625;
    for(int k=0;k<5;k++)for(int j=k;j<5;j++){
        int choose=1;for(int i=1;i<=k;i++)choose=choose*(j-i+1)/i;
        left.p[k]=left.p[k]+right.p[j]*(k%2?-choose:choose);
    }
    double error=0;int samples=0;
    for(bool end:{false,true}){const auto& arc=end?right:left;const double limit=2*P::sampleSurface(arc,end?1:0).curvature;
        for(int k=1;k<=52;k++){const double gap=std::ldexp(1.,-k),t=end?1-gap:gap;
            const auto q=P::sampleSurfaceRZ(arc,t,boundary);++samples;
            require(q.point.x>0&&q.speed>0&&std::isfinite(q.curvature),"factored pole neighborhood is not regular");
            if(k>=35)error=std::max(error,std::abs(q.curvature-limit)/std::max(1.,std::abs(limit)));
        }
    }
    require(error<1e-9,"factored curvature does not approach the analytic pole limit");
    P::Arc interior{.2,.3,{{{.3,.1},{.1,.2},{.02,-.01},{0,0},{0,0}}}};
    require(P::sampleSurfaceRZ(interior,.37,boundary).curvature==P::sampleSurfaceRZ(interior,.37).curvature,
        "non-pole arc changed its curvature evaluation");
    P::Arc negative{0,.5,{{{0,0},{1,0},{-2,.5},{0,0},{0,0}}}};
    bool rejected=false;try{P::sampleSurfaceRZ(negative,.75,boundary);}catch(const std::runtime_error&){rejected=true;}
    require(rejected,"declared pole hid a genuinely negative interior radius");
    std::cout<<std::setprecision(17)<<"PARAMETRIC_FACTORED_POLE samples="<<samples<<" recorded_radius="<<recorded.point.x
        <<" recorded_curvature="<<recorded.curvature<<" limit_error="<<error<<" negative_radius_rejected="<<rejected<<'\n';
}
void hydroParametricCurvatureContract(){
    hydroParametricPoleContract();
    Output output;bool accepted=true;
    for(int kind=0;kind<3;kind++){
        std::array<Summary,3> result;for(int level=0;level<3;level++)result[level]=study({kind},8<<level,output);
        const double order1=std::log2(result[0].curvature/result[1].curvature),order2=std::log2(result[1].curvature/result[2].curvature);
        std::cout<<"PARAMETRIC_CURVATURE_ORDER case="<<Reference{kind}.name()<<" coarse="<<order1<<" fine="<<order2
            <<" finest="<<result[2].curvature<<std::endl;
        accepted=accepted&&result[2].curvature<.01&&order1>=2&&order2>=2;
    }
    require(accepted,"smooth curvature has not reached the declared one-percent and second-order target");
}
void hydroParametricNonlinearCurvatureContract(){
    Output output;bool accepted=true;
    for(int kind:{3,4}){
        std::array<Summary,3> result;for(int level=0;level<3;level++)result[level]=study({kind},8<<level,output);
        const double order1=std::log2(result[0].curvature/result[1].curvature),order2=std::log2(result[1].curvature/result[2].curvature);
        std::cout<<"PARAMETRIC_NONLINEAR_CURVATURE_ORDER case="<<Reference{kind}.name()<<" coarse="<<order1<<" fine="<<order2<<" finest="<<result[2].curvature<<std::endl;
        accepted=accepted&&result[2].curvature<.01&&order1>=2&&order2>=2;
    }
    require(accepted,"nonlinear curvature has not reached the declared one-percent and second-order target");
}
