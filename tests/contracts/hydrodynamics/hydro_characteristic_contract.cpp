#include "Modules/Hydrodynamics/HydroCharacteristic.h" // Public API stands alone.
#include "Modules/Hydrodynamics/Detail/HydroRegularAdvectorGeometry.h"
#include <cstring>
#include "interface_first_jet_fixture.h"
#include <array>
#include <cfenv>
#include <fstream>
#include <sstream>
#if defined(_M_X64) || defined(__SSE__)
#include <xmmintrin.h>
#endif

namespace
{
using namespace LAMP::Hydrodynamics;
using Matrix=LevelSetBoundaryTrace::Jacobian;
double c(Vector2D p,int d){return d==0?p.x:p.y;}
double error(double a,double b){return std::abs(a-b)/(1+std::abs(b));}
std::array<double,4> entries(Matrix a){return {a.xx,a.xy,a.yx,a.yy};}
HydroPhysicalDivergence metrics(const QuadtreeMesh& mesh,bool rz,int selected=0)
{
    FaceSegmentGraph graph;graph.build(mesh,[selected](const QuadtreeMesh&,const Cell& cell){const auto p=cell.centerCoord();
        if(selected==1)return p.x<.5;
        if(selected==2)return !(p.x>.25&&p.x<.75&&p.y>.25&&p.y<.75);
        if(selected==3)return p.x<.25||p.x>.75;
        return true;});
    return HydroPhysicalDivergence::build(mesh,HydroIntegralDivergence::build(mesh,HydroDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,rz))));
}
using Field=std::function<Vector2D(Vector2D)>;
std::vector<double> nativeMeans(const QuadtreeMesh& mesh,const HydroDofLayout& layout,const Field& f)
{
    std::vector<double> v(layout.velocityCount());const bool rz=layout.domain().axisymmetric();
    for(int u=0;u<layout.velocityCount();++u) {
        const auto face=layout.domain().face(mesh,layout.nativeFace(u));const int d=face.horizontal?1:0;
        const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
        const auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
        long double integral=0,measure=0;
        for(double z:{-.57735026918962576451,.57735026918962576451}){const auto p=a+(b-a)*((1+z)*.5);const double weight=rz?p.x:1;
            integral+=c(f(p),d)*weight;measure+=weight;}
        v[u]=measure==0?0:static_cast<double>(integral/measure);
    }
    return v;
}
HydroCharacteristicSource sourceFor(const QuadtreeMesh& mesh,bool rz,Field f,bool varying=false,int selected=0,
    std::optional<HermiteGalsPublication>* retainedOriginal=nullptr)
{
    const auto g=HydroMomentumGeometry::build(mesh,metrics(mesh,rz,selected));const auto r=HydroAdvectorReconstruction::prepare(mesh,g);
    const auto search=HydroCharacteristicGeometry::prepare(mesh,r);Case shape;shape.pureFirstJet=true;shape.radial=rz;shape.planar=true;shape.planeOffset=2;
    const auto p=initial(mesh,shape),a=initial(mesh,shape);if(retainedOriginal)*retainedOriginal=p;
    const auto values=nativeMeans(mesh,g.layout(),f);
    auto before=values,after=values;for(auto& v:before)v*=varying?.7:1.;for(auto& v:after)v*=varying?1.3:1.;
    const std::vector<double> pressure(g.layout().pressureCount());
    return search.linear(mesh,r.evaluate(mesh,HydroNativeValues(g.layout(),before,pressure),p.state().identity(),p.state().timeLayer()),
        r.evaluate(mesh,HydroNativeValues(g.layout(),after,pressure),a.state().identity(),{1,0,.25,Position::End}));
}
Field smooth(bool rz)
{return [rz](Vector2D p){return rz?Vector2D{.2*p.x,.13}:Vector2D{-(p.y-.5),p.x-.5};};}
std::pair<Vector2D,Matrix> exact(Vector2D p,bool rz,double integral)
{
    if(rz){const double e=std::exp(-.2*integral);return {{p.x*e,p.y-.13*integral},{e,0,0,1}};}
    const double co=std::cos(integral),si=std::sin(integral);
    return {{.5+co*(p.x-.5)+si*(p.y-.5),.5-si*(p.x-.5)+co*(p.y-.5)},{co,si,-si,co}};
}
double finiteDifference(const HydroCharacteristicSource& source,Vector2D p,double at,double until,double step,std::optional<Matrix> supplied={})
{
    const auto result=source.traceBackward(p,at,until,{step});check(result.jacobian().has_value()&&!result.oneSided(),"FD reference not an interior map");
    const auto actual=entries(supplied?*supplied:*result.jacobian());double peak=0;const double h=2e-7;
    for(int d=0;d<2;++d) {
        auto left=p,right=p;(d==0?left.x:left.y)-=h;(d==0?right.x:right.y)+=h;
        const auto a=source.traceBackward(left,at,until,{step}),b=source.traceBackward(right,at,until,{step});
        check(a.jacobian().has_value()&&b.jacobian().has_value(),"FD crossed physical boundary");
        const auto measured=(b.point()-a.point())*(.5/h);
        peak=std::max({peak,error(actual[d],measured.x),error(actual[2+d],measured.y)});
    }
    check(peak<=5e-7,"actual event-map derivative differs from independent FD");return peak;
}
// The former GL9 point-jet path is an independent oracle for the new
// retained-piece powers. All sums include the original parameter span width.
Vector2D curvePoint(Vector2D a,Vector2D b,Vector2D d,double t)
{return a+(b*t+d*(t*t));}
Matrix curveMatrix(Matrix a,Matrix b,Matrix d,double t)
{return {a.xx+b.xx*t+d.xx*t*t,a.xy+b.xy*t+d.xy*t*t,
    a.yx+b.yx*t+d.yx*t*t,a.yy+b.yy*t+d.yy*t*t};}
std::array<double,4> curveProduct(Matrix a,Matrix b)
{return {a.xx*b.xx+a.xy*b.yx,a.xx*b.xy+a.xy*b.yy,
    a.yx*b.xx+a.yy*b.yx,a.yx*b.xy+a.yy*b.yy};}
template<class Query> HydroAdvectorCurveIntegral curveGL9(Query query,Vector2D r0,Vector2D r1,Vector2D r2,
    Matrix P0,Matrix P1,Matrix P2,double lo,double hi)
{
    constexpr double nodes[9]={-.96816023950762608984,-.83603110732663579430,-.61337143270059039731,-.32425342340380892904,0,.32425342340380892904,.61337143270059039731,.83603110732663579430,.96816023950762608984};
    constexpr double weights[9]={.08127438836157441197,.18064816069485740406,.26061069640293546232,.31234707704000284007,.33023935500125976316,.31234707704000284007,.26061069640293546232,.18064816069485740406,.08127438836157441197};
    HydroAdvectorCurveIntegral result;
    for(int q=0;q<9;++q){const double z=(1+nodes[q])*.5,t=(1-z)*lo+z*hi;
        const auto jet=query(curvePoint(r0,r1,r2,t),t);const auto product=curveProduct(jet.derivative,curveMatrix(P0,P1,P2,t));
        const long double weight=(static_cast<long double>(hi)-lo)*weights[q]*.5L;
        for(int k=0;k<2;++k){const long double v=c(jet.value,k);result.V0[k]+=weight*v;result.V1[k]+=weight*t*v;}
        for(int k=0;k<4;++k){result.K0[k]+=weight*product[k];result.K1[k]+=weight*t*product[k];}}
    return result;
}
double curveIntegralError(const HydroAdvectorCurveIntegral& a,const HydroAdvectorCurveIntegral& b)
{
    long double peak=0;const auto compare=[&](long double x,long double y){check(std::isfinite(x)&&std::isfinite(y),"nonfinite curve witness");peak=std::max(peak,std::abs(x-y)/(1+std::abs(y)));};
    for(int k=0;k<2;++k){compare(a.V0[k],b.V0[k]);compare(a.V1[k],b.V1[k]);}
    for(int k=0;k<4;++k){compare(a.K0[k],b.K0[k]);compare(a.K1[k],b.K1[k]);}
    // Same normalized 1e-12 affine-source/point gate as the original suite.
    check(peak<=1e-12L,"actual curve moments differ from independent reference");return static_cast<double>(peak);
}
HydroAdvectorCurveIntegral curveIntegralMean(HydroAdvectorCurveIntegral value,long double width)
{
    check(width>0,"curve witness needs positive actual span");
    for(auto* a:{&value.V0,&value.V1})for(auto& x:*a)x/=width;
    for(auto* a:{&value.K0,&value.K1})for(auto& x:*a)x/=width;
    return value;
}
bool sameCurveIntegral(const HydroAdvectorCurveIntegral& a,const HydroAdvectorCurveIntegral& b)
{return a.V0==b.V0&&a.V1==b.V1&&a.K0==b.K0&&a.K1==b.K1;}
void curveIntegralContract()
{
    double peak=0;int endpoints=0,pairedEndpoints=0,tinyEndpoints=0,constantCurveEndpoints=0;
    for(int scene=0;scene<3;++scene){const bool rz=scene==1;auto mesh=makeMesh(4,scene==2?1:0);
        const Field field=scene==0?Field([](Vector2D){return Vector2D{.23,-.17};}):scene==1?smooth(true):
            Field([](Vector2D p){return Vector2D{.2+.08*p.x*p.x+.03*p.x*p.y+.02*p.y*p.y*p.y,
                -.1+.05*p.y*p.y-.04*p.x*p.x*p.y+.07*p.x*p.x*p.x};});
        const auto source=sourceFor(*mesh,rz,field,true);const auto& pieces=source.geometry().reconstruction().pieces();
        int piece=-1;double minWidth=1,maxWidth=0;
        for(std::size_t i=0;i<pieces.size();++i){const auto& a=pieces[i];const double width=a.upper.x-a.lower.x;
            minWidth=std::min(minWidth,width);maxWidth=std::max(maxWidth,width);
            if(piece<0&&a.lower.x>0&&a.upper.x<1&&a.lower.y>0&&a.upper.y<1)piece=static_cast<int>(i);}
        check(piece>=0,"curve witness has no actual interior piece");if(scene==2)check(minWidth<maxWidth,"AMR curve witness lacks nonuniform pieces");
        const auto bounds=pieces[piece];const double hx=bounds.upper.x-bounds.lower.x,hy=bounds.upper.y-bounds.lower.y;
        const Vector2D r0=bounds.lower+Vector2D{.35*hx,.4*hy},r1{.12*hx,-.08*hy},r2{-.035*hx,.055*hy};
        const Matrix P0{1,.12,-.07,.9},P1{.15,-.11,.08,.05},P2{-.04,.03,.02,-.06};const double lo=.13,hi=.87;
        std::array<HydroAdvectorCurveIntegral,2> actual,references;
        for(int endpoint=0;endpoint<2;++endpoint){const auto& spatial=endpoint==0?source.previous():source.arrival();
            const auto id=spatial.identity();const auto geometry=spatial.geometryIdentity();const auto clock=spatial.timeLayer();
            const auto native=spatial.values().normalVelocity();const auto polynomials=spatial.regularAtomPolynomials();
            actual[endpoint]=spatial.pieceCurveIntegral(piece,r0,r1,r2,P0,P1,P2,lo,hi);
            const auto reference=curveGL9([&](Vector2D p,double){const auto a=spatial.pieceVelocityJet(piece,p);
                check(a.lower.x==bounds.lower.x&&a.lower.y==bounds.lower.y&&a.upper.x==bounds.upper.x&&a.upper.y==bounds.upper.y,"actual piece support changed");return a.jet;},r0,r1,r2,P0,P1,P2,lo,hi);
            references[endpoint]=reference;
            peak=std::max(peak,curveIntegralError(actual[endpoint],reference));
            const auto warm=spatial.pieceCurveIntegral(piece,r0,r1,r2,P0,P1,P2,lo,hi);
            check(actual[endpoint].V0==warm.V0&&actual[endpoint].V1==warm.V1&&actual[endpoint].K0==warm.K0&&actual[endpoint].K1==warm.K1,"warm curve moments changed");
            const auto now=spatial.timeLayer();check(spatial.identity()==id&&spatial.geometryIdentity()==geometry&&now.step==clock.step&&now.begin==clock.begin&&now.end==clock.end&&now.position==clock.position,"curve query changed actual owner/clock");
            check(native==spatial.values().normalVelocity()&&polynomials.size()==spatial.regularAtomPolynomials().size(),"curve query changed native source/face inventory");
            for(std::size_t i=0;i<polynomials.size();++i){const auto& a=polynomials[i];const auto& b=spatial.regularAtomPolynomials()[i];
                check(a.beginCoefficients==b.beginCoefficients&&a.endCoefficients==b.endCoefficients&&a.length==b.length&&a.radiusBegin==b.radiusBegin&&a.radiusEnd==b.radiusEnd&&a.radialMeasure==b.radialMeasure&&a.radialTangent==b.radialTangent,"curve query changed actual face polynomial");}
            if(scene<2){const long double factor=endpoint==0?.7L:1.3L;std::array<long double,4> m{};
                for(int k=0;k<4;++k)m[k]=(std::pow(static_cast<long double>(hi),k+1)-std::pow(static_cast<long double>(lo),k+1))/(k+1);
                HydroAdvectorCurveIntegral analytic;
                if(scene==0){for(int k=0;k<2;++k){analytic.V0[k]=factor*c(field(r0),k)*m[0];analytic.V1[k]=factor*c(field(r0),k)*m[1];}}
                else{analytic.V0={.2L*factor*(r0.x*m[0]+r1.x*m[1]+r2.x*m[2]),.13L*factor*m[0]};
                    analytic.V1={.2L*factor*(r0.x*m[1]+r1.x*m[2]+r2.x*m[3]),.13L*factor*m[1]};
                    const auto a=entries(P0),b=entries(P1),d=entries(P2);
                    for(int k=0;k<2;++k){analytic.K0[k]=.2L*factor*(a[k]*m[0]+b[k]*m[1]+d[k]*m[2]);analytic.K1[k]=.2L*factor*(a[k]*m[1]+b[k]*m[2]+d[k]*m[3]);}}
                peak=std::max(peak,curveIntegralError(actual[endpoint],analytic));}
            ++endpoints;}
        const auto previousId=source.previous().identity(),arrivalId=source.arrival().identity();
        const auto previousClock=source.previous().timeLayer(),arrivalClock=source.arrival().timeLayer();
        const auto paired=source.previous().pieceCurveIntegralPair(source.arrival(),piece,r0,r1,r2,P0,P1,P2,lo,hi);
        const std::array<HydroAdvectorCurveIntegral,2> pairValues{paired.previous,paired.arrival};
        for(int endpoint=0;endpoint<2;++endpoint){
            // GL9 uses the original actual endpoint jets, independently of both
            // retained single and shared-pair polynomial contractions.
            peak=std::max(peak,curveIntegralError(pairValues[endpoint],references[endpoint]));
            peak=std::max(peak,curveIntegralError(pairValues[endpoint],actual[endpoint]));++pairedEndpoints;}
        const auto pairWarm=source.previous().pieceCurveIntegralPair(source.arrival(),piece,r0,r1,r2,P0,P1,P2,lo,hi);
        check(sameCurveIntegral(paired.previous,pairWarm.previous)&&sameCurveIntegral(paired.arrival,pairWarm.arrival),"warm paired curve moments changed");
        const auto unchangedClock=[](const auto& a,const auto& b){return a.step==b.step&&a.begin==b.begin&&a.end==b.end&&a.position==b.position;};
        check(source.previous().identity()==previousId&&source.arrival().identity()==arrivalId&&
            unchangedClock(source.previous().timeLayer(),previousClock)&&unchangedClock(source.arrival().timeLayer(),arrivalClock),"paired query changed actual source/clock");
        // An actual time-dependent endpoint source, without a refit or owner.
        // alpha(theta)=.32+.44theta on the original [0,.25] clock.
        HydroAdvectorCurveIntegral combined,combinedReference;
        for(int k=0;k<2;++k)combined.V0[k]=.68L*pairValues[0].V0[k]-.44L*pairValues[0].V1[k]+.32L*pairValues[1].V0[k]+.44L*pairValues[1].V1[k];
        for(int k=0;k<4;++k)combined.K0[k]=.68L*pairValues[0].K0[k]-.44L*pairValues[0].K1[k]+.32L*pairValues[1].K0[k]+.44L*pairValues[1].K1[k];
        const auto timed=curveGL9([&](Vector2D p,double t){return source.sample(piece,p,.08+.11*t).spatial;},r0,r1,r2,P0,P1,P2,lo,hi);
        combinedReference.V0=timed.V0;combinedReference.K0=timed.K0;peak=std::max(peak,curveIntegralError(combined,combinedReference));
        if(scene==2){const auto a=entries(source.previous().pieceVelocityJet(piece,r0).jet.derivative);
            const auto b=entries(source.previous().pieceVelocityJet(piece,curvePoint(r0,r1,r2,hi)).jet.derivative);double variation=0;
            for(int k=0;k<4;++k)variation=std::max(variation,std::abs(a[k]-b[k]));check(variation>1e-8,"AMR curve source is accidentally affine");}
        if(scene==2){const auto& spatial=source.arrival();
            const double tinyLo=.37,tinyHi=std::nextafter(tinyLo,1.);
            const long double tinyWidth=static_cast<long double>(tinyHi)-tinyLo;
            const auto tiny=source.previous().pieceCurveIntegralPair(spatial,piece,r0,r1,r2,P0,P1,P2,tinyLo,tinyHi);
            const auto fixed=source.previous().pieceCurveIntegralPair(spatial,piece,r0,{0,0},{0,0},P0,P1,P2,lo,hi);
            const std::array<HydroAdvectorCurveIntegral,2> tinyValues{tiny.previous,tiny.arrival},fixedValues{fixed.previous,fixed.arrival};
            for(int endpoint=0;endpoint<2;++endpoint){const auto& endpointSource=endpoint==0?source.previous():spatial;
                const auto reference=curveGL9([&](Vector2D p,double){return endpointSource.pieceVelocityJet(piece,p).jet;},r0,r1,r2,P0,P1,P2,tinyLo,tinyHi);
                // Original 1e-12 gate on the mean, so a dropped tiny integral
                // cannot pass merely because its absolute span is tiny.
                peak=std::max(peak,curveIntegralError(curveIntegralMean(tinyValues[endpoint],tinyWidth),curveIntegralMean(reference,tinyWidth)));++tinyEndpoints;
                const auto fixedReference=curveGL9([&](Vector2D p,double){return endpointSource.pieceVelocityJet(piece,p).jet;},r0,{0,0},{0,0},P0,P1,P2,lo,hi);
                peak=std::max(peak,curveIntegralError(fixedValues[endpoint],fixedReference));++constantCurveEndpoints;}
            rejects([&]{source.previous().pieceCurveIntegralPair(spatial,-1,r0,r1,r2,P0,P1,P2,lo,hi);},"piece index out of range");
            rejects([&]{source.previous().pieceCurveIntegralPair(spatial,piece,r0,r1,r2,P0,P1,P2,lo,lo);},"invalid curve span");
            rejects([&]{source.previous().pieceCurveIntegralPair(spatial,piece,{bounds.upper.x+hx,r0.y},{0,0},{0,0},P0,P1,P2,lo,hi);},"curve outside actual piece");
            const auto foreign=sourceFor(*mesh,rz,field,true);
            rejects([&]{source.previous().pieceCurveIntegralPair(foreign.arrival(),piece,r0,r1,r2,P0,P1,P2,lo,hi);},"different reconstruction");
            rejects([&]{spatial.pieceCurveIntegral(-1,r0,r1,r2,P0,P1,P2,lo,hi);},"piece index out of range");
            rejects([&]{spatial.pieceCurveIntegral(piece,r0,r1,r2,P0,P1,P2,hi,lo);},"invalid curve span");
            rejects([&]{spatial.pieceCurveIntegral(piece,{bounds.upper.x+hx,r0.y},{0,0},{0,0},P0,P1,P2,lo,hi);},"curve outside actual piece");
            auto bad=r1;bad.x=std::numeric_limits<double>::infinity();rejects([&]{spatial.pieceCurveIntegral(piece,r0,bad,r2,P0,P1,P2,lo,hi);},"nonfinite curve data");
            rejects([&]{source.previous().pieceCurveIntegralPair(spatial,piece,r0,bad,r2,P0,P1,P2,lo,hi);},"nonfinite curve data");
            auto badP=P1;badP.yy=std::numeric_limits<double>::quiet_NaN();rejects([&]{spatial.pieceCurveIntegral(piece,r0,r1,r2,P0,badP,P2,lo,hi);},"nonfinite curve data");
            rejects([&]{source.previous().pieceCurveIntegralPair(spatial,piece,r0,r1,r2,P0,badP,P2,lo,hi);},"nonfinite curve data");
            const auto rounding=std::fegetround();std::fesetround(FE_DOWNWARD);
            try{rejects([&]{spatial.pieceCurveIntegral(piece,r0,r1,r2,P0,P1,P2,lo,hi);},"round-to-nearest");
                rejects([&]{source.previous().pieceCurveIntegralPair(spatial,piece,r0,r1,r2,P0,P1,P2,lo,hi);},"round-to-nearest");}catch(...){std::fesetround(rounding);throw;}std::fesetround(rounding);
#if defined(_M_X64) || defined(__SSE__)
            const auto csr=_mm_getcsr();_mm_setcsr(csr|0x8040u);
            try{rejects([&]{spatial.pieceCurveIntegral(piece,r0,r1,r2,P0,P1,P2,lo,hi);},"gradual underflow");
                rejects([&]{source.previous().pieceCurveIntegralPair(spatial,piece,r0,r1,r2,P0,P1,P2,lo,hi);},"gradual underflow");}catch(...){_mm_setcsr(csr);throw;}_mm_setcsr(csr);
#endif
            const auto restored=spatial.pieceCurveIntegral(piece,r0,r1,r2,P0,P1,P2,lo,hi);
            check(restored.V0==actual[1].V0&&restored.V1==actual[1].V1&&restored.K0==actual[1].K0&&restored.K1==actual[1].K1,"FP restoration changed warm curve source");
            const auto pairRestored=source.previous().pieceCurveIntegralPair(spatial,piece,r0,r1,r2,P0,P1,P2,lo,hi);
            check(std::fegetround()==rounding&&sameCurveIntegral(pairRestored.previous,paired.previous)&&sameCurveIntegral(pairRestored.arrival,paired.arrival),"FP restoration changed paired source/environment");}
        std::cout<<"CURVE_INTEGRAL_SCENE,"<<scene<<','<<rz<<','<<pieces.size()<<','<<piece<<','<<peak<<",PASS\n";}
    check(endpoints==6,"curve endpoint coverage incomplete");std::cout<<"CURVE_INTEGRAL,"<<endpoints<<','<<peak<<",PASS\n";
    check(pairedEndpoints==6&&tinyEndpoints==2&&constantCurveEndpoints==2,"paired curve endpoint/span coverage incomplete");
    std::cout<<"CURVE_INTEGRAL_PAIR,"<<pairedEndpoints<<','<<tinyEndpoints<<','<<constantCurveEndpoints<<','<<peak<<",PASS\n";
}
void certificateNegativeControls()
{
    int count=0;
    for(double scale:{1.,1e-20,1e-200})for(int component=0;component<6;++component)for(int degree:{0,3,4}) {
        LAMP::Hydrodynamics::Detail::HydroQuarticTrace a{},b{};
        for(int node=0;node<5;++node){const double t=node*.25;a[0][node]=scale;
            const double invisible=(t-.25)*(t-.5)*(t-.75);b[component][node]=scale*(degree==0?1:invisible*(degree==4?t:1));}
        for(int k=0;k<6;++k)for(int node=0;node<5;++node)b[k][node]+=a[k][node];
        check(LAMP::Hydrodynamics::Detail::hydroQuarticTraceAgreement(a,a).continuous,"identical trace rejected");
        check(!LAMP::Hydrodynamics::Detail::hydroQuarticTraceAgreement(a,b).continuous,"constant/high-degree trace mismatch hidden");
        if(degree>0)for(int node=1;node<4;++node)check(a[component][node]==b[component][node],"negative control visible to old three nodes");++count;
    }
    LAMP::Hydrodynamics::Detail::HydroQuarticTrace zero{};
    check(LAMP::Hydrodynamics::Detail::hydroQuarticTraceAgreement(zero,zero).continuous,"zero certificate failed");
    std::cout<<"COMPLETE_CERTIFICATE_NEGATIVE,"<<count<<",PASS\n";
}
void derivativeMatrix()
{
    double peak=0;int rows=0,events=0;
    for(bool rz:{false,true})for(int amr:{0,1})for(int n:{4,8,16})for(bool varying:{false,true}) {
        auto mesh=makeMesh(n,amr);const auto source=sourceFor(*mesh,rz,[rz](Vector2D p){return rz?Vector2D{.2*p.x*(1+.2*p.y*p.y),.22+.15*p.x*p.x-.04*p.y}:
            Vector2D{.17+.13*p.x*p.y,.21-.11*p.x*p.x+.07*p.y};},varying);
        for(int i=0;i<12;++i){const Vector2D p{.30137+.02671*(i%4),.32713+.03717*(i/4)};
            peak=std::max(peak,finiteDifference(source,p,.25,.07,.035));events+=source.traceBackward(p,.25,.07,{.035}).events();++rows;}
    }
    check(rows==288&&events>0,"general temporal source FD coverage missing");std::cout<<"EVENT_FD,"<<rows<<','<<events<<','<<peak<<",PASS\n";
}
bool originalFailure(bool rz,int amr,int n,bool varying)
{return rz&&((amr==0&&n==8)||(amr==1&&n==4)||(varying&&((amr==0&&n==16)||(amr==1&&n>=8))));}
void timeConvergence(bool originalsOnly=false)
{
    int groups=0,failures=0;double minPoint=100,minJet=100,maxFine=0;int events=0;
    for(bool rz:{false,true})for(int amr:{0,1})for(int n:{4,8,16})for(bool varying:{false,true}) {
        if(originalsOnly&&!originalFailure(rz,amr,n,varying))continue;
        auto mesh=makeMesh(n,amr);const auto source=sourceFor(*mesh,rz,smooth(rz),varying);double oldP=0,oldJ=0,orderP=0,orderJ=0;bool decreased=true;
        double peakFine=0;Vector2D worstPoint{0,0};Matrix worstMatrix;
        double sourceValue=0,sourceJet=0;const auto& pieces=source.geometry().reconstruction().pieces();
        for(std::size_t index=0;index<pieces.size();++index){const auto& piece=pieces[index];
            for(Vector2D p:std::vector<Vector2D>{piece.lower,piece.upper,{piece.lower.x,piece.upper.y},{piece.upper.x,piece.lower.y},piece.lower+(piece.upper-piece.lower)*.5})
                for(double time:{0.,.25}){const double factor=varying?(time==0?.7:1.3):1.;const auto actual=source.sample(static_cast<int>(index),p,time).spatial;
                    const auto reference=smooth(rz)(p)*factor;sourceValue=std::max({sourceValue,error(actual.value.x,reference.x),error(actual.value.y,reference.y)});
                    const auto got=entries(actual.derivative),expected=rz?std::array<double,4>{.2*factor,0,0,0}:std::array<double,4>{0,-factor,factor,0};
                    for(int k=0;k<4;++k)sourceJet=std::max(sourceJet,error(got[k],expected[k]));}}
        check(sourceValue<=1e-12&&sourceJet<=1e-12,"same affine source was not reproduced");
        const auto& continuity=source.continuityReport();
        check(continuity.connections>0&&continuity.connections==continuity.transparentConnections,"affine partitions are not roundoff C1");
        std::cout<<"CONTINUITY,"<<rz<<','<<amr<<','<<n<<','<<varying<<','<<continuity.connections<<','<<continuity.transparentConnections<<','<<continuity.maxTransparentRelativeJump<<'\n';
        std::cout<<"AFFINE_SOURCE,"<<rz<<','<<amr<<','<<n<<','<<varying<<','<<sourceValue<<','<<sourceJet<<'\n';
        for(int steps:{20,40,80}) {
            long double ep=0,ej=0;
            for(int i=0;i<32;++i) {
                const Vector2D p{.2631+.01413*(i%8),.3417+.01517*(i/8)};const auto reference=exact(p,rz,.25);
                const auto result=source.traceBackward(p,.25,0,{.25/steps});check(result.jacobian().has_value(),"smooth trajectory hit boundary");events+=result.events();
                check(result.splitEvents()==0,"artificial affine connection restarted the trajectory");
                ep+=std::pow(error(result.point().x,reference.first.x),2)+std::pow(error(result.point().y,reference.first.y),2);
                const auto actual=entries(*result.jacobian()),ref=entries(reference.second);for(int k=0;k<4;++k){const double e=error(actual[k],ref[k]);ej+=e*e;
                    if(steps==80&&e>peakFine){peakFine=e;worstPoint=p;worstMatrix=*result.jacobian();}}
                std::cout<<"TRACE_POINT,"<<rz<<','<<amr<<','<<n<<','<<varying<<','<<steps<<','<<i<<','<<result.events()<<','<<p.x<<','<<p.y<<','<<result.point().x<<','<<result.point().y;
                for(double v:actual)std::cout<<','<<v;std::cout<<'\n';
            }
            const double pe=std::sqrt(static_cast<double>(ep/64)),je=std::sqrt(static_cast<double>(ej/128));
            if(steps>20){decreased=decreased&&pe<oldP&&je<oldJ;orderP=std::log2(oldP/pe);orderJ=std::log2(oldJ/je);}oldP=pe;oldJ=je;
            std::cout<<"TRACE_TIME,"<<rz<<','<<amr<<','<<n<<','<<varying<<','<<steps<<','<<pe<<','<<je<<'\n';
        }
        const bool passed=decreased&&orderP>=1.8&&orderJ>=1.8&&oldP<=1e-4&&oldJ<=1e-4;if(!passed)++failures;
        minPoint=std::min(minPoint,orderP);minJet=std::min(minJet,orderJ);maxFine=std::max({maxFine,oldP,oldJ});++groups;
        std::cout<<"TRACE_ORDER,"<<rz<<','<<amr<<','<<n<<','<<varying<<','<<orderP<<','<<orderJ<<','<<decreased<<','<<(passed?"PASS":"FAIL")<<'\n';
        if(!passed){const double fd=finiteDifference(source,worstPoint,.25,0,.25/80);
            std::cout<<"TRACE_FAILURE_WITNESS,"<<rz<<','<<amr<<','<<n<<','<<varying<<','<<worstPoint.x<<','<<worstPoint.y<<','<<peakFine<<','<<fd;
            for(double v:entries(worstMatrix))std::cout<<','<<v;std::cout<<'\n';}
        if(originalFailure(rz,amr,n,varying)) {
            const Vector2D p{.36201,amr==1&&n==16?.35687:.38721};const auto result=source.traceBackward(p,.25,0,{.25/80});
            const double fd=finiteDifference(source,p,.25,0,.25/80);
            std::cout<<"ORIGINAL_WITNESS,"<<rz<<','<<amr<<','<<n<<','<<varying<<','<<p.x<<','<<p.y<<','<<result.jacobian()->xy<<','<<fd<<'\n';
        }
    }
    check(groups==(originalsOnly?7:24)&&events>0,"full temporal matrix incomplete");std::cout<<"TRACE_TIME_SUMMARY,"<<groups<<','<<events<<','<<minPoint<<','<<minJet<<','<<maxFine<<','<<failures<<'\n';
    check(failures==0,"smooth event-map temporal order gate failed");
}
void sourceAndBoundary()
{
    int rows=0;
    for(bool rz:{false,true})for(int amr:{0,1}) {
        auto mesh=makeMesh(4,amr);const auto source=sourceFor(*mesh,rz,smooth(rz),true);const auto copy=source;
        const auto p=Vector2D{.37123,.41137};const auto located=source.geometry().locatePieces(p);check(located.size()==1,"interior locate count");
        const auto begin=source.sample(located[0],p,0),end=source.sample(located[0],p,.25),mid=source.sample(located[0],p,.125);
        close(mid.spatial.value.x,.5*(begin.spatial.value.x+end.spatial.value.x),1e-15,"same endpoint interpolation");
        const auto old=source.previous().pieceVelocityJet(located[0],p).jet;close(begin.spatial.value.x,old.value.x,0,"endpoint value changed");
        source.validateFor(*mesh,source.previous().identity(),source.arrival().identity());
        rejects([&]{source.sample(located[0],p,-.01);},"outside actual slab");rejects([&]{source.sample(-1,p,0);},"piece index");
        rejects([&]{source.traceBackward(p,.25,0,{0});},"work budget");rejects([&]{source.traceBackward(p,.26,0,{.01});},"interval");
        rejects([&]{source.traceBackward({-1,0},.25,0,{.01});},"outside selected domain");
        rejects([&]{source.validateFor(*mesh,{},source.arrival().identity());},"endpoint identity");
        rejects([&]{source.traceBackward(p,.25,0,{.001,1,256});},"substep work budget");
        const auto mode=std::fegetround();std::fesetround(FE_UPWARD);try{rejects([&]{source.sample(located[0],p,.1);},"round-to-nearest");}catch(...){std::fesetround(mode);throw;}std::fesetround(mode);
        const auto result=source.traceBackward(p,.25,.24,{.01});
        LevelSetBoundaryTrace::Domain domain;domain.minCoord={0,0};domain.maxCoord={1,1};domain.cellSize={.125,.125};domain.radialAxisSymmetry=rz;
        const auto now=[&](const Vector2D& x){return source.sample(located[0],x,.25).spatial;};
        const auto then=[&](const Vector2D& x){return source.sample(located[0],x,.24).spatial;};
        const auto core=LevelSetBoundaryTrace::Detail::Tracer<true>(domain,.01,now,then).run(p);
        // Independent exact integral of the quadratic Picard predictor for
        // a(s)=a1-2.4*s. Static affine data give the original RK3 polynomial.
        const double h=.01,a1=1.3,a2=a1-2.4*h/3,a3=a1-2.4*2*h/3;
        const double slope=-2.4,mean=a1+slope*h/2;
        const double second=.5*mean*mean,third=.5*a1*a2*(a1/3+slope*h/4);
        const double radial=1-.2*h*mean+.04*h*h*second-.008*h*h*h*third;
        const double diagonal=1-h*h*second,off=h*mean-h*h*h*third;
        const Matrix reference=rz?Matrix{radial,0,0,1}:Matrix{diagonal,off,-off,diagonal};
        const Vector2D expected=rz?Vector2D{radial*p.x,p.y-.13*h*mean}:
            Vector2D{.5+diagonal*(p.x-.5)+off*(p.y-.5),.5-off*(p.x-.5)+diagonal*(p.y-.5)};
        close(result.point().x,expected.x,1e-12,"affine Picard polynomial x");close(result.point().y,expected.y,1e-12,"affine Picard polynomial y");
        const auto a=entries(*result.jacobian()),b=entries(reference);for(int k=0;k<4;++k)close(a[k],b[k],1e-12,"affine Picard stability polynomial");
        const auto truth=exact(p,rz,.7*h+1.2*(.25*.25-.24*.24));
        const double newError=std::hypot(result.point().x-truth.first.x,result.point().y-truth.first.y);
        const double coreError=std::hypot(core.trajectory.point.x-truth.first.x,core.trajectory.point.y-truth.first.y);
        std::cout<<"CORE_COMPARISON,"<<rz<<','<<amr<<','<<newError<<','<<coreError<<'\n';
        if(rz){const auto axis=source.traceBackward({0,.43137},.25,0,{.01});check(axis.jacobian().has_value()&&axis.oneSided()&&axis.point().x==0,"axis trace semantics");}
        auto equivalent=makeMesh(4,amr);mesh.reset();copy.validateFor(*equivalent,copy.previous().identity(),copy.arrival().identity());
        check(copy.sameIdentity(source)&&copy.traceBackward(p,.25,.24,{.01}).point().x==result.point().x,"source/copy lifetime changed result");
        DictionaryMesh::RefineBox box;box.level=1;box.min={0,0};box.max={1,1};equivalent->applyRefineBoxes({box});
        rejects([&]{copy.validateFor(*equivalent,copy.previous().identity(),copy.arrival().identity());},"NativeMeshIdentity");
        const auto fresh=sourceFor(*equivalent,rz,smooth(rz));fresh.validateFor(*equivalent,fresh.previous().identity(),fresh.arrival().identity());check(!fresh.sameIdentity(copy),"regrid reused temporal source");
        check(fresh.traceBackward(p,.25,.24,{.01}).jacobian().has_value(),"fresh regrid trace failed");++rows;
    }
    for(int selection:{0,1,2,3}) {
        auto mesh=makeMesh(8,1);const auto source=sourceFor(*mesh,false,[](Vector2D){return Vector2D{1,0};},false,selection);
        const Vector2D p=selection==2?Vector2D{.83,.43137}:Vector2D{.17,.18137};
        const auto result=source.traceBackward(p,.25,0,{.04});check(result.boundary().has_value()&&!result.jacobian().has_value(),"selected/outer boundary missing");
        close(result.point().x,selection==2?.75:0,1e-12,"actual boundary location");check(result.boundary()->normalAxis==0&&result.boundary()->outwardSign==-1,"wrong actual boundary normal");++rows;
        if(selection==3)rejects([&]{source.traceBackward({.5,.3},.25,0,{.01});},"outside selected domain");
    }
    std::cout<<"TRACE_SOURCE_BOUNDARY,"<<rows<<",PASS\n";
}
// Geometry/event assertions stay exact; a nonzero-time numerical Jacobian
// uses the existing affine matrix tolerance, with every entry checked.
void constantMapJacobian(const HydroCharacteristicResult& result)
{
    check(result.jacobian().has_value(),"constant map lacks interior Jacobian");
    const auto got=entries(*result.jacobian());const std::array<double,4> exact{1,0,0,1};
    for(int k=0;k<4;++k)close(got[k],exact[k],1e-12,"constant map affine Jacobian");
}
// Independent intersections with the actual support partition, including an
// exact terminal face and excluding the initial face. No fixed old-piece count.
int constantCrossings(const HydroCharacteristicSource& source,Vector2D start,double endX)
{
    std::set<double> crossed;
    for(const auto& piece:source.geometry().reconstruction().pieces())if(start.y>piece.lower.y&&start.y<piece.upper.y)
        for(double x:{piece.lower.x,piece.upper.x})if(x>=endX&&x<start.x)crossed.insert(x);
    return static_cast<int>(crossed.size());
}
void constantNumericalGuards()
{
    int rows=0,negative=0;double peak=0;
    for(bool rz:{false,true})for(int amr:{0,1})for(double magnitude:{0.,1e-20,1.}) {
        auto mesh=makeMesh(8,amr);const auto source=sourceFor(*mesh,rz,[&](Vector2D){return Vector2D{0,magnitude};});
        const Vector2D p{.313,.43137};const auto piece=source.geometry().locatePieces(p).front();const auto v=source.sample(piece,p,.25).spatial.value;
        check(magnitude==0?(v.x==0&&v.y==0):(std::abs(v.x)<=1e-12*magnitude&&std::abs(v.y-magnitude)<=1e-12*magnitude),"constant sample lost scale");
        const auto map=source.traceBackward(p,.25,.24,{.01});constantMapJacobian(map);
        peak=std::max(peak,finiteDifference(source,p,.25,.24,.01));
        for(int k=0;k<4;++k){auto wrong=*map.jacobian();(k==0?wrong.xx:k==1?wrong.xy:k==2?wrong.yx:wrong.yy)+=1e-3;
            rejects([&]{finiteDifference(source,p,.25,.24,.01,wrong);},"differs from independent FD");++negative;}++rows;
    }
    check(rows==12&&negative==48,"constant numerical guard coverage");std::cout<<"CONSTANT_NUMERICAL_GUARDS,"<<rows<<','<<negative<<','<<peak<<",PASS\n";
}
void eventAndSourceGuards()
{
    auto mesh=makeMesh(8);const auto constant=sourceFor(*mesh,false,[](Vector2D){return Vector2D{1,0};});
    const Vector2D p{.313,.43137};
    const auto zero=constant.traceBackward(p,.25,.25,{.1});
    check(zero.point().x==p.x&&zero.point().y==p.y&&zero.jacobian()->xx==1&&zero.jacobian()->yy==1&&zero.events()==0,"zero-duration identity");
    const auto terminal=constant.traceBackward({.3125,.43137},.25,.1875,{.0625});
    check(terminal.point().x==.25&&terminal.oneSided()&&terminal.events()==constantCrossings(constant,{.3125,.43137},.25)&&terminal.splitEvents()==0,"exact final connection semantics");constantMapJacobian(terminal);
    const auto node=constant.traceBackward({.3125,.43137},.25,.125,{.0625});
    check(node.point().x==.1875&&node.events()==constantCrossings(constant,{.3125,.43137},.1875)&&node.splitEvents()==0,"exact intermediate clock connection");constantMapJacobian(node);
    const auto tiny=constant.traceBackward({std::nextafter(.25,1.),.43137},.25,.24,{.01});
    check(tiny.events()==constantCrossings(constant,{std::nextafter(.25,1.),.43137},tiny.point().x)&&tiny.events()==1&&tiny.splitEvents()==0,"tiny transparent crossing lost");constantMapJacobian(tiny);
    close(tiny.point().x,std::nextafter(.25,1.)-(.25-.24),1e-15,"tiny crossing location");
    const auto jump=sourceFor(*mesh,false,[](Vector2D x){return Vector2D{.03*std::sin(15*x.y),1};});
    const Vector2D jumpStart{.43713,std::nextafter(.5,1.)};
    const auto tinyJump=jump.traceBackward(jumpStart,.25,.24,{.01});
    check(tinyJump.jacobian().has_value()&&tinyJump.events()>0&&tinyJump.splitEvents()==0,"tiny C1 crossing was restarted");
    finiteDifference(jump,jumpStart,.25,.24,.01);
    const auto outward=constant.traceBackward({0,.43137},.25,0,{.1});
    check(outward.boundary().has_value()&&!outward.jacobian()&&outward.time()==.25,"initial external hit");
    const auto ids=constant.geometry().locatePieces({.25,.43137});check(ids.size()==2,"one-sided fixture");
    int incoming=-1,outgoing=-1;for(int i:ids){if(constant.geometry().reconstruction().pieces()[i].lower.x==.25)outgoing=i;else incoming=i;}
    rejects([&]{constant.traceBackward({.25,.43137},.25,.24,{.01},outgoing);},"initial branch points outside");
    check(constant.traceBackward({.25,.43137},.25,.24,{.01},incoming).oneSided(),"declared one-sided start");
    rejects([&]{constant.traceBackward(p,.25,0,{.1},-2);},"initial piece selector");
    rejects([&]{constant.traceBackward(p,.25,0,{.1},incoming);},"does not contain");
    rejects([&]{constant.traceBackward({.73137,.43137},.25,0,{.04,4096,1});},"event work budget");
    const auto corner=sourceFor(*mesh,false,[](Vector2D){return Vector2D{1,1};});
    const auto crossing=corner.traceBackward({.313,.313},.25,.125,{.125});
    close(crossing.point().x,.188,1e-12,"C1 simultaneous x");close(crossing.point().y,.188,1e-12,"C1 simultaneous y");constantMapJacobian(crossing);
    const auto graze=sourceFor(*mesh,false,[](Vector2D x){return Vector2D{.3125-x.y,1};},false,1);
    rejects([&]{graze.traceBackward({.5-.5*.0625*.0625,.375},.25,.125,{.125});},"grazing");
    const auto stationary=sourceFor(*mesh,false,[](Vector2D){return Vector2D{0,0};});
    const auto resting=stationary.traceBackward({.25,.43137},.25,0,{.1});
    check(resting.point().x==.25&&resting.point().y==.43137,"stationary skeleton identity");constantMapJacobian(resting);
    check(stationary.traceBackward({.25,.43137},.25,0,{.1},ids.front()).jacobian().has_value(),"explicit stationary branch");
    std::cout<<"STATIONARY_BRANCH,PASS\n";
    for(bool rz:{false,true})for(int selection:{0,1,2,3}) {
        std::cout<<"SELECTED_AXIAL,"<<rz<<','<<selection<<'\n';
        const auto source=sourceFor(*mesh,rz,[](Vector2D){return Vector2D{0,1};},false,selection);
        const Vector2D start=selection==2?Vector2D{.43137,.83}:Vector2D{.17,.17};
        const auto result=source.traceBackward(start,.25,0,{.04});
        check(result.boundary().has_value()&&!result.jacobian()&&result.boundary()->normalAxis==1&&result.boundary()->outwardSign==-1,"selected axial boundary");
        close(result.point().y,selection==2?.75:0,1e-12,"selected axial location");
    }
    const auto& search=constant.geometry();const auto& r=search.reconstruction();const auto& g=r.momentumGeometry();
    const auto& pSource=constant.previous();const auto& aSource=constant.arrival();
    const auto endpoint=[&](LAMP::Interface::HermiteGalsTimeLayer clock){return r.evaluate(*mesh,aSource.values(),aSource.geometryIdentity(),clock);};
    for(auto bad:std::vector<LAMP::Interface::HermiteGalsTimeLayer>{{1,0,.25,Position::Midpoint},{2,0,.25,Position::End},{1,.01,.25,Position::End},{1,.25,.25,Position::Initial}})
        rejects([&]{search.linear(*mesh,pSource,endpoint(bad));},"endpoint clock");
    rejects([&]{search.linear(*mesh,endpoint({1,0,.1,Position::Midpoint}),aSource);},"endpoint clock");
    const auto r2=HydroAdvectorReconstruction::prepare(*mesh,g);
    rejects([&]{search.linear(*mesh,pSource,r2.evaluate(*mesh,aSource.values(),aSource.geometryIdentity(),aSource.timeLayer()));},"different actual reconstruction");
    const auto boundary=HydroAdvectorBoundary::impermeable(*mesh,g,{}),otherBoundary=HydroAdvectorBoundary::impermeable(*mesh,g,{});
    const auto bp=r.evaluate(*mesh,pSource.values(),pSource.geometryIdentity(),pSource.timeLayer(),boundary);
    const auto ba=r.evaluate(*mesh,aSource.values(),aSource.geometryIdentity(),aSource.timeLayer(),boundary);
    const auto other=r.evaluate(*mesh,aSource.values(),aSource.geometryIdentity(),aSource.timeLayer(),otherBoundary);
    check(search.linear(*mesh,bp,ba).traceBackward(p,.25,.24,{.01}).jacobian().has_value(),"same boundary identity refused");
    rejects([&]{search.linear(*mesh,bp,aSource);},"boundary law");rejects([&]{search.linear(*mesh,bp,other);},"boundary law");
    std::cout<<"EVENT_SOURCE_GUARDS,PASS\n";
}
void continuityAndPartition()
{
    int rows=0;double peak=0;
    for(bool rz:{false,true})for(bool varying:{false,true}) {
        auto base=makeMesh(4);const auto reference=sourceFor(*base,rz,smooth(rz),varying);
        for(int amr:{0,1})for(int n:{4,8,16}) {
            auto mesh=makeMesh(n,amr);const auto source=sourceFor(*mesh,rz,smooth(rz),varying);
            for(int steps:{20,40,80})for(int i=0;i<32;++i) {
                const Vector2D p{.2631+.01413*(i%8),.3417+.01517*(i/8)};
                const auto a=source.traceBackward(p,.25,0,{.25/steps}),b=reference.traceBackward(p,.25,0,{.25/steps});
                peak=std::max({peak,error(a.point().x,b.point().x),error(a.point().y,b.point().y)});
                const auto aj=entries(*a.jacobian()),bj=entries(*b.jacobian());for(int k=0;k<4;++k)peak=std::max(peak,error(aj[k],bj[k]));++rows;
            }
        }
    }
    check(rows==2304&&peak<=1e-12,"smooth partition consistency");std::cout<<"PARTITION_MAP,"<<rows<<','<<peak<<",PASS\n";
    for(bool rz:{false,true})for(int amr:{0,1})for(double magnitude:{0.,1e-20,1.}) {
        auto mesh=makeMesh(8,amr);const auto affine=sourceFor(*mesh,rz,[&](Vector2D p){return smooth(rz)(p)*magnitude;},true);
        const auto nonlinear=sourceFor(*mesh,rz,[&](Vector2D p){return Vector2D{(rz?p.x:1)*.03*std::sin(15*p.y),1}*magnitude;},true);
        const auto& a=affine.continuityReport();const auto& b=nonlinear.continuityReport();
        check(a.connections==a.transparentConnections,"scaled affine certificate not homogeneous");
        check(b.connections==b.transparentConnections,"regular nonlinear C1 certificate failed");
        // Zero previous and nonlinear arrival still share complete C1 traces.
        const auto& r=nonlinear.geometry().reconstruction();const auto& g=r.momentumGeometry();
        const HydroNativeValues zero(g.layout(),std::vector<double>(g.layout().velocityCount()),std::vector<double>(g.layout().pressureCount()));
        const auto old=r.evaluate(*mesh,zero,nonlinear.previous().geometryIdentity(),nonlinear.previous().timeLayer());
        const auto mixed=nonlinear.geometry().linear(*mesh,old,nonlinear.arrival());const auto& m=mixed.continuityReport();
        check(m.connections==m.transparentConnections,"mixed endpoint C1 certificate failed");
        std::cout<<"CERTIFICATE_SCALE,"<<rz<<','<<amr<<','<<magnitude<<','<<a.connections<<','<<a.transparentConnections<<','<<b.transparentConnections<<','<<m.transparentConnections<<",PASS\n";
    }
}

// Independent continuous-flow reference of these ORIGINAL native inputs:
// docs/task/H08a2r16_Resolved_Nonlinear_Motion_Reference_Report.md. Exact
// binary-coefficient piece polynomials, directed Taylor remainders and the
// bounded u_y perturbation give the point/J budgets below. Reference values
// never enter the source or accepted map. The 4e-16 allowance covers their
// conversion to binary64. No unresolved unpartitioned RK4 truth is retained.
void trustedRegularWitness(bool internal)
{
    auto mesh=makeMesh(8);
    if(internal){DictionaryMesh::RefineBox left;left.level=1;left.min={0,0};left.max={.375,1};auto right=left;right.min={.5,0};right.max={1,1};mesh->applyRefineBoxes({left,right});}
    const auto source=sourceFor(*mesh,false,[](Vector2D p){return Vector2D{.03*std::sin(15*p.y),1};});
    const double top=internal?.48137:.57137,bottom=internal?.36137:.41137,end=.25-(top-bottom);const Vector2D p{.43713,top};

    // Independent 100-digit Taylor/variational flow of actual binary Q4
    // coefficients. H14a3r2 updates only the changed uniform P3 source:
    // runs/checks/H14a3r2/20261010_192509_P3_source/characteristic_reference.
    // Internal AMR coefficients/reference remain byte-identical to H13
    // runs/checks/H13/20261009_003903; no production tracer supplies truth.
    const double refs[2][8]={{0.4353136501899857,0.41137,1.002583520296169,-0.011493854427436617,0,1,1.1688773415792806e-16,1.5689343994837832e-14},{0.4369693706869503,0.36137,0.9999910118882219,-0.03035256960438956,0,1,1.5078805650408715e-16,1.2608319572785833e-14}};const auto* r=refs[internal];
    const double referencePoint=r[6]+4e-16,referenceJac=r[7]+4e-16;
    check(std::max(referencePoint,referenceJac)<=1e-12,"trusted reference unresolved");
    double previous=0,order=0;int events=0;
    for(int steps:{40,80,160}){
        const auto result=source.traceBackward(p,.25,end,{(.25-end)/steps});
        check(result.jacobian().has_value()&&result.events()>0&&result.splitEvents()==0,"common C1 witness restarted");
        const double e=std::hypot(result.point().x-r[0],result.point().y-r[1]);
        if(previous>0)order=std::log2(previous/e);previous=e;events=result.events();
        std::cout<<"TRUSTED_SEQUENCE,"<<internal<<','<<steps<<','<<e<<','<<order<<'\n';
    }
    check(order>=1.8,"common-source witness convergence failed");
    const auto fine=source.traceBackward(p,.25,end,{(.25-end)/2560});check(fine.jacobian().has_value(),"fine witness lacks J");
    const double pointError=std::hypot(fine.point().x-r[0],fine.point().y-r[1]);
    double jacError=0;const auto fj=entries(*fine.jacobian());for(int k=0;k<4;++k)jacError=std::max(jacError,std::abs(fj[k]-r[k+2]));
    // Execute FD and the certificate even when the final accuracy gate fails.
    const double fd=finiteDifference(source,p,.25,end,(.25-end)/160);
    const auto& certificate=source.continuityReport();check(certificate.connections==certificate.transparentConnections,"common regular C1 missing");
    std::cout<<"TRUSTED_WITNESS,"<<internal<<','<<events<<','<<pointError<<','<<jacError<<','<<referencePoint<<','<<referenceJac<<','<<order<<','<<fd<<','<<certificate.connections<<','<<certificate.transparentConnections<<std::endl;
    check(pointError+referencePoint<=1e-11,"common-source point witness accuracy failed");
    check(jacError+referenceJac<=1e-11,"common-source Jacobian witness accuracy failed");
}
void firstJetSkeleton()
{
    int count=0;double peak=0;
    for(bool rz:{false,true})for(int n:{4,8})for(int amr:{0,1})for(int kind=0;kind<(rz?2:3);++kind){
        auto mesh=makeMesh(n,amr);const Vector2D velocity=kind==0?Vector2D{0,0}:kind==1?Vector2D{0,1}:Vector2D{1,1};
        const auto source=sourceFor(*mesh,rz,[velocity](Vector2D){return velocity;});
        Case shape;shape.pureFirstJet=true;shape.planar=true;shape.planeOffset=2;shape.radial=rz;
        const auto state=initial(*mesh,shape);std::vector<Vector2D> points;
        for(const auto& v:state.state().primary().vertices)points.push_back(point(*mesh,v.vertex));
        for(const auto& e:state.state().primary().edges){const auto a=point(*mesh,e.vertices[0]),b=point(*mesh,e.vertices[1]);points.push_back(a+(b-a)*.5);}
        for(int cell:state.state().primary().geometryCells)points.push_back(mesh->cells(mesh->cellIndexToLabel[cell]).centerCoord());
        for(auto p:points){const auto expected=p-velocity*.25;
            // Physical boundary hits deliberately have no interior Jacobian.
            // This matrix checks only fully internal straight trajectories.
            if(p.x<=0||p.x>=1||p.y<=0||p.y>=1||expected.x<=0||expected.x>=1||expected.y<=0||expected.y>=1)continue;
            const auto result=[&](){try{return source.traceBackward(p,.25,0,{.0125});}
                catch(const std::exception&){std::cout<<"SKELETON_FAILURE,"<<rz<<','<<n<<','<<amr<<','<<kind<<','<<p.x<<','<<p.y<<std::endl;throw;}}();
            check(result.jacobian().has_value(),"interior skeleton motion missing J");
            const auto j=entries(*result.jacobian());const auto id=entries(Matrix::identity());
            double e=std::max(std::abs(result.point().x-expected.x),std::abs(result.point().y-expected.y));
            for(int k=0;k<4;++k)e=std::max(e,std::abs(j[k]-id[k]));
            check(e<=1e-12,"analytic FirstJet skeleton point/J");peak=std::max(peak,e);++count;
        }
    }
    std::cout<<"H08A_SKELETON,"<<count<<','<<peak<<",PASS"<<std::endl;
}
// A real native/source path with an exactly integrable linear clock. Its
// tiny inward start precedes a near-start transverse C1 crossing; neither
// the endpoint nor its Jacobian may depend on the event classifier scale.
void nearStartCrossingRegression()
{
    auto mesh=makeMesh(8);const double tiny=std::ldexp(1.,-55),at=.25,h=1./64;
    const Vector2D p{.43137,.5};
    const auto endpoint=sourceFor(*mesh,false,[tiny](Vector2D){return Vector2D{0,tiny};});
    const auto& search=endpoint.geometry();const auto& r=search.reconstruction();const auto& layout=r.momentumGeometry().layout();
    const auto previousValues=nativeMeans(*mesh,layout,[](Vector2D){return Vector2D{0,-1};});
    const auto previous=r.evaluate(*mesh,HydroNativeValues(layout,previousValues,std::vector<double>(layout.pressureCount())),
        endpoint.previous().geometryIdentity(),endpoint.previous().timeLayer());
    const auto source=search.linear(*mesh,previous,endpoint.arrival());
    check(source.continuityReport().connections==source.continuityReport().transparentConnections,"near-start source lacks C1 certificate");
    const auto result=source.traceBackward(p,at,at-h,{h});
    check(result.jacobian().has_value()&&!result.boundary()&&result.events()==1&&result.splitEvents()==0&&result.substeps()==1,
        "near-start transverse connection lost or restarted");
    const double displacement=-tiny*h+(1+tiny)*h*h/(2*at);
    close(result.point().x,p.x,1e-12,"near-start exact translation x");
    close(result.point().y,p.y+displacement,1e-12,"near-start exact translation y");constantMapJacobian(result);
    // A skeleton start is marked one-sided by the public contract, but this
    // globally spatially constant source has the same smooth map on both sides.
    const auto actual=entries(*result.jacobian());double peak=0;const double delta=2e-7;
    for(int d=0;d<2;++d) {
        auto left=p,right=p;(d==0?left.x:left.y)-=delta;(d==0?right.x:right.y)+=delta;
        const auto a=source.traceBackward(left,at,at-h,{h}),b=source.traceBackward(right,at,at-h,{h});
        check(a.jacobian().has_value()&&b.jacobian().has_value(),"near-start FD hit physical boundary");
        const auto measured=(b.point()-a.point())*(.5/delta);
        peak=std::max({peak,error(actual[d],measured.x),error(actual[2+d],measured.y)});
    }
    check(peak<=5e-7,"near-start actual map Jacobian differs from independent FD");
    std::cout<<"NEAR_START_TRANSVERSE_REGRESSION,"<<result.events()<<','<<peak<<",PASS"<<std::endl;
}
// A real native endpoint field whose sub-ulp tangential excursion rounds
// to a grid corner. A transverse return can point into a diagonal C1 chart;
// rejecting its single face neighbour must not invent a physical boundary.
void roundedCornerCrossingRegression()
{
    auto mesh=makeMesh(8);const Vector2D point{.5,.75};const double at=.25;
    int count=0;
    for(int exponent:{-60,-90})for(bool swap:{false,true}) {
        const double tiny=std::ldexp(1.,exponent);
        const auto field=[&](Vector2D value){return Field([value,swap](Vector2D) {
            return swap?Vector2D{value.y,value.x}:value;});};
        const auto endpoint=sourceFor(*mesh,false,field({-tiny,-2*tiny}));
        const auto& search=endpoint.geometry();const auto& r=search.reconstruction();
        const auto& layout=r.momentumGeometry().layout();
        const auto before=nativeMeans(*mesh,layout,field({3*tiny,(10./3)*tiny}));
        const auto previous=r.evaluate(*mesh,HydroNativeValues(layout,before,
            std::vector<double>(layout.pressureCount())),endpoint.previous().geometryIdentity(),
            endpoint.previous().timeLayer());
        const auto source=search.linear(*mesh,previous,endpoint.arrival());
        check(source.continuityReport().connections==source.continuityReport().transparentConnections,
            "rounded corner source lacks C1 certificate");
        const auto result=source.traceBackward(point,at,0,{at});
        check(!result.boundary()&&result.jacobian().has_value()&&result.time()==0,
            "rounded internal corner fabricated a physical boundary");
        check(result.point().x==point.x&&result.point().y==point.y,
            "sub-ulp corner fixture changed its represented endpoint");
        constantMapJacobian(result);++count;
    }
    std::cout<<"ROUNDED_CORNER_CROSSING_REGRESSION,"<<count<<",PASS"<<std::endl;
}
// Original endpoint single queries remain independent witnesses for the
// shared geometric evaluation. Every field below is a real native source.
void pairedStageJetContract()
{
    int points=0,bubblePoints=0,nativeTraces=0,fallbackTraces=0;double peak=0;
    for(int scene=0;scene<3;++scene) {
        const bool rz=scene==1;auto mesh=makeMesh(4,scene==2?1:0);
        const Field field=scene==2?Field([](Vector2D p){return Vector2D{.03*std::sin(15*p.y),1};}):smooth(rz);
        std::optional<HermiteGalsPublication> nativeOriginal;
        const auto source=sourceFor(*mesh,rz,field,true,0,&nativeOriginal);
        const auto& previous=source.previous();const auto& arrival=source.arrival();
        check(previous.values().normalVelocity()!=arrival.values().normalVelocity(),
            "paired stage witness has identical endpoint fields");
        const auto previousId=previous.identity(),arrivalId=arrival.identity();
        const auto previousClock=previous.timeLayer(),arrivalClock=arrival.timeLayer();
        const auto nativeBefore=previous.values().normalVelocity(),nativeAfter=arrival.values().normalVelocity();
        const auto& regular=previous.reconstruction().regularGeometry();
        const auto coefficients=regular.apply(*mesh,previous.values(),previous.boundary());
        const auto& pieces=previous.reconstruction().pieces();int selected=-1;
        double minWidth=1,maxWidth=0;bool activeBubble=false;
        for(std::size_t i=0;i<pieces.size();++i) {
            const auto& piece=pieces[i];const auto span=piece.upper-piece.lower;
            minWidth=std::min(minWidth,span.x);maxWidth=std::max(maxWidth,span.x);
            if(!(piece.lower.x>0&&piece.upper.x<1&&piece.lower.y>0&&piece.upper.y<1))continue;
            const auto point=piece.lower+Vector2D{.371*span.x,.419*span.y};bool bubble=false;
            for(int u:regular.cells()[piece.pressure].faces) {
                const auto& face=regular.faces()[u];const int n=face.normalAxis,t=1-n;
                const double width=face.supportWidth[piece.pressure==face.positivePressure?1:0];
                const double length=c(face.upper-face.lower,t);
                const double local=c(point-face.lower,t)/length,normal=c(point-face.lower,n)/width;
                if(local>0&&local<1&&normal>-1&&normal<1&&std::abs(coefficients.faceCorrections[u])>1e-10)bubble=true;
            }
            if(selected<0||(scene==2&&bubble&&!activeBubble)){selected=static_cast<int>(i);activeBubble=bubble;}
        }
        check(selected>=0,"paired stage witness lacks real interior piece");
        if(scene==2)check(minWidth<maxWidth&&activeBubble,"paired stage AMR/bubble witness not active");
        if(scene==2){const auto& box=pieces[selected];const auto span=box.upper-box.lower;
            const auto a=entries(previous.pieceVelocityJet(selected,box.lower+span*.35).jet.derivative),
                z=entries(previous.pieceVelocityJet(selected,box.lower+span*.65).jet.derivative);
            double variation=0;for(int k=0;k<4;++k)variation=std::max(variation,std::abs(a[k]-z[k]));
            check(variation>1e-8,"native-start AMR source is accidentally affine");}
        const auto bounds=pieces[selected];const auto span=bounds.upper-bounds.lower;
        const Vector2D interior=bounds.lower+Vector2D{.371*span.x,.419*span.y};
        const std::array<Vector2D,5> queries{{interior,{bounds.lower.x,interior.y},
            {bounds.upper.x,interior.y},{interior.x,bounds.lower.y},{interior.x,bounds.upper.y}}};
        const auto compare=[&](const HydroAdvectorCellJet& actual,const HydroAdvectorCellJet& reference) {
            check(actual.lower.x==reference.lower.x&&actual.lower.y==reference.lower.y
                &&actual.upper.x==reference.upper.x&&actual.upper.y==reference.upper.y,
                "paired stage changed actual piece bounds");
            const std::array<double,6> a{{actual.jet.value.x,actual.jet.value.y,actual.jet.derivative.xx,
                actual.jet.derivative.xy,actual.jet.derivative.yx,actual.jet.derivative.yy}};
            const std::array<double,6> b{{reference.jet.value.x,reference.jet.value.y,reference.jet.derivative.xx,
                reference.jet.derivative.xy,reference.jet.derivative.yx,reference.jet.derivative.yy}};
            for(int k=0;k<6;++k){check(std::isfinite(a[k])&&std::isfinite(b[k]),"nonfinite paired stage witness");
                peak=std::max(peak,error(a[k],b[k]));
                check(std::memcmp(&a[k],&b[k],sizeof(double))==0,"paired stage changed original single component bits");}
        };
        for(const auto point:queries) {
            const auto originalPrevious=previous.pieceVelocityJet(selected,point),originalArrival=arrival.pieceVelocityJet(selected,point);
            const auto paired=previous.pieceVelocityJetPair(arrival,selected,point);
            const auto warm=previous.pieceVelocityJetPair(arrival,selected,point);
            compare(paired.previous,originalPrevious);compare(paired.arrival,originalArrival);
            compare(warm.previous,paired.previous);compare(warm.arrival,paired.arrival);++points;
            if(scene==2)++bubblePoints;
        }
        // Native arrival starts come from exactly the original FirstJet
        // inventory. The existing AMR/bubble source is genuinely non-affine;
        // the unchanged trace port remains the independent map/J oracle.
        check(nativeOriginal.has_value(),"native-start witness lost actual original publication");
        const auto& primary=nativeOriginal->state().primary();
        const auto starts=source.prepareNativeStarts(*mesh,*nativeOriginal),copiedStarts=starts;
        const auto copiedSource=source;
        std::set<std::pair<double,double>> nativeKeys;
        std::array<std::optional<Vector2D>,3> selectedPoints;
        std::array<double,3> nearest{{std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity()}};
        const auto addNative=[&](int kind,Vector2D q){
            nativeKeys.emplace(q.x,q.y);
            if(q.x<=.1||q.x>=.9||q.y<=.1||q.y>=.9)return;
            const double distance=std::hypot(q.x-interior.x,q.y-interior.y);
            if(distance<nearest[kind]){nearest[kind]=distance;selectedPoints[kind]=q;}
        };
        for(const auto& v:primary.vertices)addNative(0,point(*mesh,v.vertex));
        for(const auto& e:primary.edges){const auto a=point(*mesh,e.vertices[0]),b=point(*mesh,e.vertices[1]);
            addNative(1,a+(b-a)*.5);}
        for(int cell:primary.geometryCells)addNative(2,mesh->cells(mesh->cellIndexToLabel[cell]).centerCoord());
        check(starts.usesSource(source)&&starts.usesSource(copiedSource)&&copiedStarts.usesSource(source)
            &&starts.pointCount()==nativeKeys.size(),"native-start owner or actual point inventory changed");
        const auto sameBits=[](double actual,double reference){
            check(std::isfinite(actual)&&std::isfinite(reference)
                &&std::memcmp(&actual,&reference,sizeof(double))==0,"prepared native start changed map/J component bits");
        };
        const auto compareTrace=[&](const HydroCharacteristicResult& actual,const HydroCharacteristicResult& reference){
            sameBits(actual.point().x,reference.point().x);sameBits(actual.point().y,reference.point().y);
            sameBits(actual.time(),reference.time());
            check(actual.events()==reference.events()&&actual.splitEvents()==reference.splitEvents()
                &&actual.substeps()==reference.substeps()&&actual.initialPiece()==reference.initialPiece()
                &&actual.finalPiece()==reference.finalPiece()&&actual.oneSided()==reference.oneSided()
                &&actual.jacobian().has_value()==reference.jacobian().has_value()
                &&actual.boundary().has_value()==reference.boundary().has_value(),"prepared native start changed actual trace semantics");
            if(actual.jacobian()){const auto a=entries(*actual.jacobian()),b=entries(*reference.jacobian());
                for(int k=0;k<4;++k)sameBits(a[k],b[k]);}
            if(actual.boundary()){const auto a=*actual.boundary(),b=*reference.boundary();
                check(a.normalAxis==b.normalAxis&&a.outwardSign==b.outwardSign&&a.piece==b.piece,
                    "prepared native start changed physical boundary");}
        };
        const HydroCharacteristicConfig startConfig{.005};
        for(double at:{.2,.25}){
            const double until=at-.005;
            for(const auto& selectedPoint:selectedPoints){check(selectedPoint.has_value(),"native-start witness lacks point kind");
                const auto q=*selectedPoint;
                const auto originalTrace=source.traceBackward(q,at,until,startConfig);
                compareTrace(source.traceBackward(q,at,until,startConfig,starts),originalTrace);
                compareTrace(copiedSource.traceBackward(q,at,until,startConfig,copiedStarts),originalTrace);
                compareTrace(source.traceBackward(q,at,until,startConfig,starts,originalTrace.initialPiece()),
                    source.traceBackward(q,at,until,startConfig,originalTrace.initialPiece()));
                ++nativeTraces;
            }
            check(nativeKeys.count({interior.x,interior.y})==0,"fallback witness accidentally is a native point");
            compareTrace(source.traceBackward(interior,at,until,startConfig,starts),
                source.traceBackward(interior,at,until,startConfig));++fallbackTraces;
        }
        if(scene==2){
            const auto rebuilt=source.geometry().linear(*mesh,previous,arrival);
            check(!starts.usesSource(rebuilt)&&rebuilt.previous().identity()==previous.identity()
                &&rebuilt.arrival().identity()==arrival.identity(),"prepared starts accepted rebuilt equivalent source");
            const auto q=*selectedPoints[0];
            rejects([&]{rebuilt.traceBackward(q,.25,.245,startConfig,starts);},"source");
            rejects([&]{source.traceBackward(q,.26,.245,startConfig,starts);},"trace point or interval");
            rejects([&]{source.traceBackward(q,.25,-.001,startConfig,starts);},"trace point or interval");
            rejects([&]{source.traceBackward(q,.25,.245,{0},starts);},"trace work budget");
            {
                struct Restore {int saved=std::fegetround();~Restore(){std::fesetround(saved);}} restore;
                check(std::fesetround(FE_DOWNWARD)==0,"cannot exercise cached native-start FP guard");
                rejects([&]{source.traceBackward(q,.25,.245,startConfig,starts);},"round-to-nearest");
            }
#if defined(_M_X64) || defined(__SSE__)
            {
                struct Restore {unsigned saved=_mm_getcsr();~Restore(){_mm_setcsr(saved);}} restore;
                _mm_setcsr(restore.saved|0x8040u);
                rejects([&]{source.traceBackward(q,.25,.245,startConfig,starts);},"gradual underflow");
            }
#endif
            compareTrace(source.traceBackward(q,.25,.245,startConfig,starts),
                source.traceBackward(q,.25,.245,startConfig));
            check(starts.pointCount()==nativeKeys.size()&&starts.usesSource(source)
                &&nativeOriginal->state().identity()==previous.geometryIdentity(),"failed native trace changed actual prepared owner/source");
        }
        const auto sameClock=[](const auto& a,const auto& b){return a.step==b.step&&a.begin==b.begin&&a.end==b.end&&a.position==b.position;};
        check(previous.identity()==previousId&&arrival.identity()==arrivalId
            &&sameClock(previous.timeLayer(),previousClock)&&sameClock(arrival.timeLayer(),arrivalClock)
            &&previous.values().normalVelocity()==nativeBefore&&arrival.values().normalVelocity()==nativeAfter,
            "paired stage changed endpoint owner, full clock or native values");
        rejects([&]{source.validateFor(*mesh,previous.identity(),{});},"endpoint identity");
        if(scene==2) {
            const auto foreign=sourceFor(*mesh,rz,field,true);
            rejects([&]{previous.pieceVelocityJetPair(foreign.arrival(),selected,interior);},"different reconstruction");
            rejects([&]{previous.pieceVelocityJetPair(arrival,-1,interior);},"piece index");
            rejects([&]{previous.pieceVelocityJetPair(arrival,selected,{bounds.upper.x+span.x,interior.y});},"point outside");
            const auto rounding=std::fegetround();std::fesetround(FE_DOWNWARD);
            try{rejects([&]{previous.pieceVelocityJetPair(arrival,selected,interior);},"round-to-nearest");}
            catch(...){std::fesetround(rounding);throw;}std::fesetround(rounding);
#if defined(_M_X64) || defined(__SSE__)
            const auto csr=_mm_getcsr();_mm_setcsr(csr|0x8040u);
            try{rejects([&]{previous.pieceVelocityJetPair(arrival,selected,interior);},"gradual underflow");}
            catch(...){_mm_setcsr(csr);throw;}_mm_setcsr(csr);
#endif
            const auto restored=previous.pieceVelocityJetPair(arrival,selected,interior);
            compare(restored.previous,previous.pieceVelocityJet(selected,interior));
            compare(restored.arrival,arrival.pieceVelocityJet(selected,interior));
        }
    }
    check(points==15&&bubblePoints==5,"paired stage XY/RZ/AMR/bubble coverage");
    check(nativeTraces==18&&fallbackTraces==6,"native-start XY/RZ/AMR point/time coverage");
    std::cout<<"PREPARED_NATIVE_STARTS,"<<nativeTraces<<','<<fallbackTraces
        <<",REAL_NATIVE_VERTEX_EDGE_CENTER,NONAFFINE_AMR,BITWISE_MAP_J_EVENTS,OLD_PORT,OWNER_CLOCK_FP_FTZ,PASS"<<std::endl;
    std::cout<<"PAIRED_STAGE_JET,"<<points<<','<<bubblePoints<<','<<peak<<",SINGLE_BITS,CLOSED_PIECE,OWNER,FP,FTZ,PASS"<<std::endl;
}
// H31 caller time budget: independent affine flows assess position and J
// against their own spatial-truncation allowance. Existing fine-flow,
// event and sampling accuracy gates retain their original budgets.
void timeBudgetContract()
{
    const double at=.25,until=.245,dt=at-until,referenceSpacing=1./64;
    const double positionBudget=.01*referenceSpacing*referenceSpacing;
    const double jacobianBudget=.01*referenceSpacing;
    const std::array<Vector2D,5> points{{{.43713,.53137},{.5001,.43713},
        {.43713,.4999},{.56249,.43713},{.43713,.56251}}};
    int rows=0,totalEvents=0;
    double maxPosition=0,maxJacobian=0,maxFD=0;
    for(bool rz:{false,true})for(bool varying:{false,true}) {
        auto mesh=makeMesh(8);const auto source=sourceFor(*mesh,rz,smooth(rz),varying);
        const double integral=varying?1.3*dt-1.2*dt*dt:dt;
        double coarsePosition=0,coarseJacobian=0,finePosition=0,fineJacobian=0;
        double positionDifference=0,jacobianDifference=0,fd=0;
        int events=0,coarseSubsteps=0,fineSubsteps=0;
        for(const auto p:points) {
            const auto coarse=source.traceBackward(p,at,until,{dt});
            const auto fine=source.traceBackward(p,at,until,{dt/4});
            check(coarse.jacobian().has_value()&&fine.jacobian().has_value()
                &&!coarse.boundary()&&!fine.boundary(),"time-budget affine flow hit physical boundary");
            check(coarse.splitEvents()==0&&fine.splitEvents()==0,
                "time-budget affine connection restarted trajectory");
            if(p.x==points[0].x&&p.y==points[0].y) {
                const auto& piece=source.geometry().reconstruction().pieces()[coarse.initialPiece()];
                check(p.x>piece.lower.x&&p.x<piece.upper.x&&p.y>piece.lower.y&&p.y<piece.upper.y
                    &&coarse.point().x>piece.lower.x&&coarse.point().x<piece.upper.x
                    &&coarse.point().y>piece.lower.y&&coarse.point().y<piece.upper.y
                    &&coarse.events()==0&&coarse.initialPiece()==coarse.finalPiece(),
                    "strict-interior stage witness crossed its actual piece");
            }
            const auto ref=exact(p,rz,integral);
            coarsePosition=std::max(coarsePosition,std::hypot(coarse.point().x-ref.first.x,coarse.point().y-ref.first.y));
            finePosition=std::max(finePosition,std::hypot(fine.point().x-ref.first.x,fine.point().y-ref.first.y));
            positionDifference=std::max(positionDifference,std::hypot(coarse.point().x-fine.point().x,coarse.point().y-fine.point().y));
            const auto cj=entries(*coarse.jacobian()),fj=entries(*fine.jacobian()),rj=entries(ref.second);
            for(int k=0;k<4;++k) {
                coarseJacobian=std::max(coarseJacobian,std::abs(cj[k]-rj[k]));
                fineJacobian=std::max(fineJacobian,std::abs(fj[k]-rj[k]));
                jacobianDifference=std::max(jacobianDifference,std::abs(cj[k]-fj[k]));
            }
            fd=std::max(fd,finiteDifference(source,p,at,until,dt));
            events+=coarse.events();coarseSubsteps+=coarse.substeps();fineSubsteps+=fine.substeps();++rows;
        }
        check(events>0,"time-budget affine crossing coverage missing");
        std::cout<<"TIME_BUDGET_SCENE,"<<rz<<','<<varying<<','<<dt<<','<<integral
            <<','<<coarsePosition<<','<<coarseJacobian<<','<<finePosition<<','<<fineJacobian
            <<','<<positionDifference<<','<<jacobianDifference<<','<<positionBudget<<','<<jacobianBudget
            <<','<<fd<<','<<events<<','<<coarseSubsteps<<','<<fineSubsteps<<std::endl;
        check(coarsePosition<=positionBudget,"time-budget position allowance exceeded");
        check(coarseJacobian<=jacobianBudget,"time-budget Jacobian allowance exceeded");
        maxPosition=std::max(maxPosition,coarsePosition);maxJacobian=std::max(maxJacobian,coarseJacobian);
        maxFD=std::max(maxFD,fd);totalEvents+=events;
    }
    check(rows==20,"time-budget affine coverage incomplete");
    // Endpoint exactly on a real support face stays on the original walk;
    // its one-sided semantics and exact terminal event must survive preselection.
    { auto mesh=makeMesh(8);const auto source=sourceFor(*mesh,false,[](Vector2D){return Vector2D{1,0};});
      const auto terminal=source.traceBackward({.3125,.43137},.25,.1875,{.0625});
      check(terminal.point().x==.25&&terminal.oneSided()
        &&terminal.events()==constantCrossings(source,{.3125,.43137},.25)
        &&terminal.splitEvents()==0,"time-budget exact terminal tie changed");constantMapJacobian(terminal);
      const auto outgoing=source.traceBackward({0,.43137},.25,0,{.1});
      check(outgoing.boundary().has_value()&&!outgoing.jacobian()&&outgoing.time()==.25,
        "time-budget initial external tie changed"); }
    pairedStageJetContract();firstJetSkeleton();nearStartCrossingRegression();
    std::cout<<"TIME_BUDGET,"<<rows<<','<<totalEvents<<','<<maxPosition<<','<<maxJacobian
        <<','<<maxFD<<",PASS"<<std::endl;
}
void transparentPicardRegression()
{
    // Actual native input from the H13/H08a counterexample. Compare the same
    // point map with its returned J; keep the existing FD tolerance (5e-7).
    auto mesh=makeMesh(8);const auto m=metrics(*mesh,false);const auto g=HydroMomentumGeometry::build(*mesh,m);
    const auto& l=g.layout();std::vector<double> input(l.velocityCount());
    std::ifstream file(std::string(LAMP_CHARACTERISTIC_FIXTURE_DIR)+"near_tangent_native.csv");
    check(bool(file),"missing actual transparent-Picard regression input");std::string line;std::getline(file,line);int rows=0;
    while(std::getline(file,line)){std::replace(line.begin(),line.end(),',',' ');std::istringstream row(line);int u;double value;row>>u>>value;
        check(bool(row)&&u==rows&&u<int(input.size()),"invalid native regression row");input[u]=value;++rows;}
    check(rows==l.velocityCount(),"native regression layout changed");
    std::vector<int> walls;for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(*mesh,l.nativeFace(u));if(f.negativeCell<0||f.positiveCell<0)walls.push_back(u);}
    const auto wall=HydroAdvectorBoundary::impermeable(*mesh,g,walls);const auto r=HydroAdvectorReconstruction::prepare(*mesh,g);const auto search=HydroCharacteristicGeometry::prepare(*mesh,r);
    Case shape;shape.pureFirstJet=true;shape.planar=true;shape.planeOffset=2;const auto publication=initial(*mesh,shape);
    double peak=0;int count=0;
    for(double magnitude:{1e-4,1e-6,1e-8}){
        auto v=input;for(double& a:v)a*=magnitude;const HydroNativeValues values(l,v,std::vector<double>(l.pressureCount()));
        const auto begin=r.evaluate(*mesh,values,publication.state().identity(),publication.state().timeLayer(),wall);
        const auto end=r.evaluate(*mesh,values,publication.state().identity(),{1,0,.001,Position::End},wall);const auto source=search.linear(*mesh,begin,end);
        check(source.continuityReport().connections==source.continuityReport().transparentConnections,"regression lacks original C1 certificate");
        for(int i=1;i<8;++i){const Vector2D p{double(i)/8,.5};const auto t=source.traceBackward(p,.001,0,{.001});check(t.jacobian().has_value(),"interior regression lacks J");const auto actual=entries(*t.jacobian());
            for(double h:{1.25e-5,1.25e-6,1.25e-7})for(int d=0;d<2;++d){auto a=p,b=p;(d==0?a.x:a.y)-=h;(d==0?b.x:b.y)+=h;
                const auto ta=source.traceBackward(a,.001,0,{.001}),tb=source.traceBackward(b,.001,0,{.001});const auto measured=(tb.point()-ta.point())*(.5/h);
                peak=std::max({peak,error(actual[d],measured.x),error(actual[2+d],measured.y)});++count;}
        }
    }
    std::cout<<"TRANSPARENT_PICARD_REGRESSION,"<<rows<<','<<count<<','<<peak<<std::endl;
    check(peak<=5e-7,"actual near-tangent Picard J differs from original-map FD");
}
}
int main(int argc,char** argv){
 if(argc==2&&std::string(argv[1])=="--time-budget-only"){try{std::cout<<std::setprecision(17);timeBudgetContract();return 0;}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
 if(argc==2&&std::string(argv[1])=="--curve-integral-only"){try{std::cout<<std::setprecision(17);curveIntegralContract();return 0;}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
 if(argc==2&&std::string(argv[1])=="--original-failures"){try{timeConvergence(true);return 0;}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
 if(argc==2&&std::string(argv[1])=="--transparent-picard-regression"){try{transparentPicardRegression();return 0;}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
 if(argc==2&&std::string(argv[1])=="--near-start-crossing-regression"){try{nearStartCrossingRegression();return 0;}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
 if(argc==2&&std::string(argv[1])=="--rounded-corner-crossing-regression"){try{roundedCornerCrossingRegression();return 0;}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
 if(argc!=1){std::cerr<<"unknown characteristic contract argument"<<std::endl;return 1;}
 std::cout<<std::setprecision(17);int failures=0;
 const auto run=[&](const char* name,auto f){try{f();std::cout<<"MOTION_GROUP,"<<name<<",PASS\n";}catch(const std::exception& e){++failures;std::cout<<"MOTION_GROUP,"<<name<<",FAIL,"<<e.what()<<std::endl;}};
 run("curve_integral",curveIntegralContract);
 run("certificate_negative",[]{certificateNegativeControls();});run("derivative",[]{derivativeMatrix();});run("source_boundary",[]{sourceAndBoundary();});run("event_guards",[]{eventAndSourceGuards();});run("constant_numeric",[]{constantNumericalGuards();});run("partition",[]{continuityAndPartition();});run("original_time",[]{timeConvergence();});
 run("trusted_uniform",[]{trustedRegularWitness(false);});run("trusted_internal",[]{trustedRegularWitness(true);});
 run("first_jet_skeleton",firstJetSkeleton);
 run("near_start_crossing_regression",nearStartCrossingRegression);
 run("rounded_corner_crossing_regression",roundedCornerCrossingRegression);
 run("transparent_picard_regression",transparentPicardRegression);
 std::cout<<"MOTION_RESULT,"<<failures<<','<<(failures?"FAIL":"PASS")<<std::endl;return failures?1:0;
}
