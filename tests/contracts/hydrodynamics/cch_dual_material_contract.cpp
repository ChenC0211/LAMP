#include "Modules/Hydrodynamics/CchDualMaterial.h" // Public header stands alone.
#include "Modules/Hydrodynamics/CchSinglePhaseSaddle.h"
#include "interface_first_jet_fixture.h"
#include <chrono>
#include <functional>
#include <sstream>

namespace
{
using namespace Quantum::Hydrodynamics;
int scenarios=0,segmentChecks=0,dualChecks=0,basisChecks=0,linearChecks=0,crossings=0,saddles=0,rejected=0,timeRows=0,regrids=0;
int migrationFailures=0;
long double peakExchange=0,peakAnalytic=0,peakDefect=0,peakDecomposition=0,peakConstant=0,peakSaddleU=0,peakSaddleP=0;
long double omittedDrift=0,transferChange=0;
long double actualOmittedDrift=0,actualTransferChange=0;
static_assert(!std::is_default_constructible_v<CchDualMaterial>);
static_assert(std::is_same_v<decltype(std::declval<CchDualMaterial&>().faceRates()),const std::vector<long double>&>);
long double norm(long double a,long double b,long double scale,long double limit,const std::string& why)
{const auto e=std::abs(a-b)/(scale>0?scale:1);std::ostringstream detail;detail<<std::setprecision(17)<<why<<" error="<<e<<" actual="<<a<<" reference="<<b<<" scale="<<scale;check(std::isfinite(e) && e<=limit,detail.str());return e;}
long double cap(Vector2D a,Vector2D b,bool rz)
{return (static_cast<long double>(b.x)-a.x)*(static_cast<long double>(b.y)-a.y)*(rz?pi*(static_cast<long double>(a.x)+b.x):1.L);}
int axis(const CchDofLayout& layout,int u){return layout.normalAxis(u)==CchDofLayout::NormalAxis::First?0:1;}
CchPhysicalDivergence metricsFor(const QuadtreeMesh& mesh,bool rz,int selection=0)
{
    FaceSegmentGraph graph;graph.build(mesh,[selection](const QuadtreeMesh&,const Cell& c){const auto p=c.centerCoord();
        if(selection==1)return p.x<.5;if(selection==2)return !(p.x>.25 && p.x<.75 && p.y>.25 && p.y<.75);
        if(selection==3)return p.x<.25 || p.x>.75;return true;});
    return CchPhysicalDivergence::build(mesh,CchIntegralDivergence::build(mesh,CchDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,rz))));
}
HermiteGalsPublication publishAt(const QuadtreeMesh& mesh,const Case& s,HermiteGalsTimeLayer time)
{
    const auto f=s.profile(time.time());const auto seed=stateFor(mesh,f,s.radial,time,[](int,Vector2D,Vector2D){return 0.;},s.symmetry(mesh),HermiteGalsTransport::FirstJet);
    HermiteGalsCandidate candidate(mesh,time,seed.primary());const auto mixed=closeHermiteGalsMixedAmr(mesh,candidate,s.symmetry(mesh));
    HermiteGalsPublicationConfig config;config.volume.symmetryPlanes=s.symmetry(mesh);const auto id=candidate.identity();
    return publishHermiteGalsFirstJetGeometry(mesh,std::move(candidate),mixed,{id,centers(mesh,f)},config);
}
CchNativeValues velocity(const QuadtreeMesh& mesh,const CchDofLayout& layout,const std::function<Vector2D(Vector2D)>& field)
{
    std::vector<double> u;for(int f=0;f<layout.velocityCount();++f){const auto face=layout.domain().face(mesh,layout.nativeFace(f));const auto v=field(face.center);u.push_back(face.horizontal?v.y:v.x);}
    return CchNativeValues(layout,u,std::vector<double>(layout.pressureCount()));
}
struct Frames {CchPhaseGeometry a,m,b;CchMomentumMass ma,mm,mb;};
Frames frames(const QuadtreeMesh& mesh,const CchMomentumGeometry& geometry,HermiteGalsPublication a,HermiteGalsPublication m,HermiteGalsPublication b,double l,double g=1)
{
    const auto pa=CchPhaseGeometry::build(mesh,geometry.metrics(),std::move(a));const auto pm=CchPhaseGeometry::build(mesh,geometry.metrics(),std::move(m));const auto pb=CchPhaseGeometry::build(mesh,geometry.metrics(),std::move(b));
    return {pa,pm,pb,CchMomentumMass::build(mesh,geometry,pa,l,g),CchMomentumMass::build(mesh,geometry,pm,l,g),CchMomentumMass::build(mesh,geometry,pb,l,g)};
}
Frames staticFrames(const QuadtreeMesh& mesh,const CchMomentumGeometry& geometry,Case s,double rho=1)
{return frames(mesh,geometry,publishAt(mesh,s,{}),publishAt(mesh,s,{1,0,.04,Position::Midpoint}),publishAt(mesh,s,{1,0,.04,Position::End}),rho);}
CchDualMaterial prepare(const QuadtreeMesh& mesh,const Frames& f,const CchNativeValues& u)
{
    const auto actual=CchAdvectorReconstruction::prepare(mesh,f.mm.geometry()).evaluate(mesh,u,f.m.geometryIdentity(),f.m.timeLayer());
    const auto q=CchPhaseRate::prepareSpatial(mesh,f.m,actual);
    const auto transport=CchMaterialTransport::prepare(mesh,f.a,q,f.b,f.ma.liquidDensity(),f.ma.gasDensity());
    return CchDualMaterial::prepare(mesh,transport,f.ma,f.mm,f.mb);
}
// Independent intersections: native face -> actual adjacent cells -> half
// bounds. This never reads production overlap weights or cell diagnostics.
std::vector<CchMomentumOverlap> intersections(const QuadtreeMesh& mesh,const CchDofLayout& layout,int u)
{
    const auto face=layout.domain().face(mesh,layout.nativeFace(u));const int a=axis(layout,u);
    const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);auto e0=mesh.verts(edge.edgeToVert[0]).coord(),e1=mesh.verts(edge.edgeToVert[1]).coord();
    if(a==0?e0.y>e1.y:e0.x>e1.x)std::swap(e0,e1);std::vector<CchMomentumOverlap> result;
    for(int native:{face.negativeCell,face.positiveCell})if(native>=0){const auto& c=mesh.cells(mesh.cellIndexToLabel[native]);
        auto lo=mesh.verts(c.cellToVert[0]).coord(),hi=mesh.verts(c.cellToVert[2]).coord();const auto v=cap(lo,hi,layout.domain().axisymmetric());const auto mid=lo+(hi-lo)*.5;
        if(a==0){lo.y=e0.y;hi.y=e1.y;if(native==face.negativeCell)lo.x=mid.x;else hi.x=mid.x;}
        else{lo.x=e0.x;hi.x=e1.x;if(native==face.negativeCell)lo.y=mid.y;else hi.y=mid.y;}
        const auto volume=cap(lo,hi,layout.domain().axisymmetric());result.push_back({layout.pressureDof(native),lo,hi,volume,volume/v});}
    return result;
}
std::vector<long double> outward(const CchDualMaterial& d,bool integral)
{
    std::vector<long double> result(d.layout().velocityCount());
    const auto& values=integral?d.faceIntegrals():d.faceRates();
    for(std::size_t g=0;g<values.size();++g){const auto& s=d.geometry().segments()[g];
        if(s.negativeVelocity>=0)result[s.negativeVelocity]+=values[g];if(s.positiveVelocity>=0)result[s.positiveVelocity]-=values[g];}
    return result;
}
// Independent physical reference; phase/F outputs never enter its cuts or quadrature.
// Test-only polynomial reference, independent of production phase/root code.
// GL3 integrates degrees <= 5, including a quartic normal times RZ radius.
template<class Coefficients> long double referenceValue(const Coefficients& c,long double x)
{
    long double v=c.back();for(std::size_t k=c.size()-1;k>0;--k)v=c[k-1]+x*v;return v;
}
template<class Coefficients> long double referenceWork(const Coefficients& c,long double length)
{
    long double sum=0;for(std::size_t k=0;k<c.size();++k){long double term=std::abs(c[k]);for(std::size_t j=0;j<k;++j)term*=length;sum+=term;}return sum;
}
template<class Samples> long double referenceGauss(long double halfWidth,Samples sample)
{
    const std::array<long double,3> nodes{{-std::sqrt(3.L/5),0,std::sqrt(3.L/5)}},weights{{5.L/9,8.L/9,5.L/9}};
    long double value=0;for(int g=0;g<3;++g){const auto pair=sample(nodes[g]);value+=halfWidth*weights[g]*pair[0]*pair[1];}return value;
}
// Derivative roots partition higher-degree polynomials into monotone intervals.
// Known-factor fixtures below check this floating reference; it is not a
// certified phase-root authority and must never enter accepted geometry.
std::vector<long double> referenceRoots(std::vector<long double> c,long double lo,long double hi)
{
    while(c.size()>1 && c.back()==0)c.pop_back();std::vector<long double> out;
    const auto cut=[&](long double x){if(lo<x && x<hi)out.push_back(x);};
    if(c.size()==2){if(c[1]!=0)cut(-c[0]/c[1]);}
    else if(c.size()==3){const auto d=c[1]*c[1]-4*c[0]*c[2];if(d>0){cut((-c[1]-std::sqrt(d))/(2*c[2]));cut((-c[1]+std::sqrt(d))/(2*c[2]));}}
    else if(c.size()>3){std::vector<long double> derivative(c.size()-1);for(std::size_t k=1;k<c.size();++k)derivative[k-1]=k*c[k];
        auto critical=referenceRoots(derivative,lo,hi);critical.push_back(lo);critical.push_back(hi);std::sort(critical.begin(),critical.end());
        for(auto x:critical)if(std::abs(referenceValue(c,x))<=64*std::numeric_limits<long double>::epsilon()*referenceWork(c,std::abs(x)))cut(x);
        for(std::size_t k=1;k<critical.size();++k){auto a=critical[k-1],b=critical[k],fa=referenceValue(c,a),fb=referenceValue(c,b);
            if(fa==0 || fb==0 || std::signbit(fa)==std::signbit(fb))continue;
            for(int iteration=0;iteration<100;++iteration){const auto middle=(a+b)/2;if(middle==a || middle==b)break;const auto fm=referenceValue(c,middle);
                if(fm==0){a=b=middle;break;}if(std::signbit(fa)==std::signbit(fm)){a=middle;fa=fm;}else b=middle;}
            cut((a+b)/2);}}
    std::sort(out.begin(),out.end());out.erase(std::unique(out.begin(),out.end()),out.end());return out;
}
int referenceChoose(int n,int k)
{int out=1;for(int j=1;j<=k;++j)out=out*(n-j+1)/j;return out;}
template<std::size_t N> std::array<long double,N> referenceShift(const std::array<long double,N>& c,long double origin,int direction=1)
{
    static_assert(N>=1 && N<=5,"reference supports constant through quartic normals");std::array<long double,N> out{};
    for(std::size_t i=0;i<N;++i)for(std::size_t j=i;j<N;++j){long double term=c[j]*referenceChoose(static_cast<int>(j),static_cast<int>(i));for(std::size_t k=i;k<j;++k)term*=origin;out[i]+=(direction<0 && i%2?-term:term);}return out;
}
template<std::size_t N> long double referenceShiftedIntegral(const std::array<long double,N>& c,long double offset,long double width,
    bool rz,long double radiusOrigin,int radiusSlope,int normalSign)
{
    check(offset>=0 && width>0 && std::isfinite(offset) && std::isfinite(width) && (normalSign==1 || normalSign==-1),"invalid reference shifted range");
    const auto local=referenceShift(c,offset);const auto radius=radiusOrigin+radiusSlope*offset,half=width/2;
    // Keep width and the shifted polynomial independent of rounded endpoints.
    return normalSign*referenceGauss(half,[&](long double node){const auto t=half+half*node;
        return std::array<long double,2>{referenceValue(local,t),rz?2*pi*(radius+radiusSlope*t):1.L};});
}

template<std::size_t N>
CchPhaseMeasure referenceAtom(const CchAdvectorAtom& atom,const std::array<long double,N>& c,const Case& s,int normalSign=1)
{
    static_assert(N>=1 && N<=5,"reference supports constant through quartic normals");
    check(normalSign==1 || normalSign==-1,"invalid reference normal sign");
    // Independent analytic cuts and GL3 physical integration. Only basis tests
    // take c from the actual kinematic owner; MMS below constructs c independently.
    const int n=atom.normalAxis;const double lo=n==0?atom.begin.y:atom.begin.x,hi=n==0?atom.end.y:atom.end.x;
    const auto at=[&](long double t){return n==0?Vector2D{atom.begin.x,static_cast<double>(t)}:Vector2D{static_cast<double>(t),atom.begin.y};};
    std::vector<long double> cuts{lo,hi};const auto cut=[&](long double t){if(lo<t && t<hi)cuts.push_back(t);};
    const bool zero=s.planar && n==(s.axial()?1:0) && (n==0?atom.begin.x:atom.begin.y)==s.planeOffset;
    if(s.planar){if(n!=(s.axial()?1:0))cut(s.planeOffset);}
    else{const long double centerN=n==0?(s.radial?0:s.cx):s.cy,centerT=n==0?s.cy:(s.radial?0:s.cx);
        const long double delta=(n==0?atom.begin.x:atom.begin.y)-centerN,q=static_cast<long double>(s.radius)*s.radius-delta*delta;
        if(q>0){cut(centerT-std::sqrt(q));cut(centerT+std::sqrt(q));}}
    if(zero)for(auto root:referenceRoots(std::vector<long double>(c.begin(),c.end()),0,hi-lo))cut(lo+root);
    std::sort(cuts.begin(),cuts.end());CchPhaseMeasure out;
    for(std::size_t k=1;k<cuts.size();++k){const auto l=cuts[k-1],h=cuts[k],mid=(l+h)/2,width=(h-l)/2,d=mid-lo;const auto normal=referenceValue(c,d);
        const bool liquid=zero?s.planeSign*(-normal)<0:s.exact(at(mid),0).v<0;long double value=0;
        value=referenceGauss(width,[&](long double node){const auto t=mid+width*node,offset=t-lo;const auto p=at(t);
            return std::array<long double,2>{referenceValue(c,offset),s.radial?2*pi*p.x:1.L};});
        (liquid?out.liquid:out.gas)+=normalSign*value;
    }return out;
}
std::array<long double,3> mmsCoefficients(const CchAdvectorAtom& a,bool rz,int mode)
{
    const auto x=a.begin.x,y=a.begin.y;const int n=a.normalAxis;
    if(mode==0)return {};
    if(mode==1)return {n==0?(rz?-.13*x:-.3):-.35,0,0};
    if(mode==3)return rz?(n==0?std::array<long double,3>{}:std::array<long double,3>{.5*x*x-.015625,x,.5}):
        (n==0?std::array<long double,3>{y-.375,1,0}:std::array<long double,3>{});
    return rz?(n==0?std::array<long double,3>{.13*x,0,0}:std::array<long double,3>{.2+.155*x*x+.19*y,.31*x,.155}):
        (n==0?std::array<long double,3>{.3+.17*x-.23*y,-.23,0}:std::array<long double,3>{-.2+.31*x+.19*y,.31,0});
}
CchNativeValues mmsValues(const QuadtreeMesh& mesh,const CchDofLayout& layout,int mode)
{
    std::vector<double> values;const bool rz=layout.domain().axisymmetric();
    for(int u=0;u<layout.velocityCount();++u){const auto f=layout.domain().face(mesh,layout.nativeFace(u));const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
        auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();if(f.horizontal?a.x>b.x:a.y>b.y)std::swap(a,b);
        const auto c=mmsCoefficients({f.horizontal?1:0,-1,-1,f.center,f.center,f.area},rz,mode);long double v=c[0];
        if(rz && f.horizontal && (mode==2 || mode==3))v=mode==2?.2+.155*.5*(a.x*a.x+b.x*b.x)+.19*f.center.y:.25*(a.x*a.x+b.x*b.x)-.015625;
        values.push_back(static_cast<double>(v));}
    return CchNativeValues(layout,values,std::vector<double>(layout.pressureCount()));
}
void inspectTransfer(const QuadtreeMesh& mesh,const CchDualMaterial& d,const Case* shape=nullptr)
{
    check(d.material().isSpatial() && &d.faceRates()==&d.material().dualFaceRates(),"same-F dual route missing");
    if(shape){Case middle=*shape;const auto time=d.timeLayer().time();if(middle.motion==1){middle.planeOffset+=middle.speed*time;if(middle.radial)middle.cy+=middle.speed*time;else middle.cx+=middle.speed*time;}middle.motion=0;middle.speed=0;
        const auto& actual=d.material().midpointRate().spatialAdvector();const auto& reconstruction=actual.reconstruction();std::vector<long double> atomReference;
        for(std::size_t i=0;i<reconstruction.atoms().size();++i){const auto q=referenceAtom(reconstruction.atoms()[i],actual.atomPolynomials()[i].beginCoefficients,middle);const auto expected=d.material().liquidDensity()*q.liquid+d.material().gasDensity()*q.gas;atomReference.push_back(expected);
            peakAnalytic=std::max(peakAnalytic,norm(d.material().atomRates()[i],expected,(d.material().liquidDensity()+d.material().gasDensity())*reconstruction.atoms()[i].area,1e-12,"independent basis physical atom F"));}
        for(std::size_t g=0;g<d.faceRates().size();++g){long double expected=0;for(int id:reconstruction.dualAtoms()[g])expected+=atomReference[id];
            peakAnalytic=std::max(peakAnalytic,norm(d.faceRates()[g],expected,(d.material().liquidDensity()+d.material().gasDensity())*d.geometry().segments()[g].area,1e-12,"independent actual dual physical F"));}}
    const auto& layout=d.layout();const auto out=outward(d,false),integrated=outward(d,true);std::vector<long double> primal(layout.pressureCount()),gross(layout.pressureCount());std::array<long double,2> globalScale{};
    for(int u=0;u<layout.velocityCount();++u){const auto face=layout.domain().face(mesh,layout.nativeFace(u));const auto flux=d.material().faceRates()[u];
        for(auto owner:{std::pair<int,int>{face.negativeCell,1},{face.positiveCell,-1}})if(owner.first>=0){const int p=layout.pressureDof(owner.first);primal[p]+=owner.second*flux;gross[p]+=std::abs(flux);}}
    for(int u=0;u<layout.velocityCount();++u){long double expected=0,scale=0,mappedChange=0,mappedDefect=0;
        for(const auto& p:intersections(mesh,layout,u)){expected+=p.cellFraction*primal[p.pressure];scale+=p.cellFraction*gross[p.pressure];
            mappedChange+=p.cellFraction*d.material().cells()[p.pressure].massChange;mappedDefect+=p.cellFraction*d.material().cells()[p.pressure].integratedDefect;}
        for(int side=0;side<4;++side)for(int id:d.geometry().volumes()[u].sides[side])scale+=std::abs(d.faceRates()[id]);
        peakExchange=std::max(peakExchange,norm(d.cells()[u].spatialExchangeRate,out[u]-expected,scale,1e-12,"independent nonzero spatial exchange"));
        const auto& c=d.cells()[u];const auto integralScale=std::max(1.L,std::abs(c.massChange)+scale*d.deltaTime());
        globalScale[axis(layout,u)]+=std::abs(c.massChange)+scale*d.deltaTime();
        peakDefect=std::max(peakDefect,norm(c.outwardIntegral,integrated[u],integralScale,1e-12,"independent segment outward integral"));
        norm(c.mappedPrimalMassChange,mappedChange,integralScale,1e-12,"independent geometry change transfer");
        norm(c.mappedPrimalIntegratedDefect,mappedDefect,integralScale,1e-12,"independent cell defect transfer");
        norm(c.outwardRate,out[u],scale,1e-12,"actual dual outward rate");
        norm(c.mappedPrimalOutwardRate,expected,scale,1e-12,"actual mapped primal rate");
        peakDecomposition=std::max(peakDecomposition,norm(c.integratedDefect,mappedDefect+c.geometryTransferChange+d.deltaTime()*(out[u]-expected),std::max(1.L,std::abs(c.massChange)+scale*d.deltaTime()),1e-12,"independent defect decomposition with spatial exchange"));
        transferChange=std::max(transferChange,std::abs(c.geometryTransferChange));++dualChecks;
    }
    std::array<long double,2> boundary{},signedDefect{},change{};
    for(std::size_t g=0;g<d.faceRates().size();++g){const auto& s=d.geometry().segments()[g];norm(d.faceIntegrals()[g],d.deltaTime()*d.faceRates()[g],std::max(1.L,std::abs(d.faceIntegrals()[g])),1e-14,"dual midpoint integral");
        check(&d.faceRates()==&d.material().dualFaceRates(),"dual detached original shared F view");
        if(s.negativeVelocity<0 || s.positiveVelocity<0){const int u=s.negativeVelocity>=0?s.negativeVelocity:s.positiveVelocity;boundary[axis(d.layout(),u)]+=(s.negativeVelocity>=0?1:-1)*d.faceIntegrals()[g];}
        ++segmentChecks;
    }
    for(int u=0;u<layout.velocityCount();++u){const int a=axis(layout,u);signedDefect[a]+=d.cells()[u].integratedDefect;change[a]+=d.cells()[u].massChange;}
    for(int a=0;a<2;++a){norm(d.report().components[a].boundaryIntegral,boundary[a],std::max(1.L,globalScale[a]),1e-12,"independent dual boundary");norm(signedDefect[a],change[a]+boundary[a],std::max(1.L,globalScale[a]),1e-12,"global dual material balance");}
    d.validateFor(mesh,d.sources(),d.timeLayer());++scenarios;
}
void basisAndLinearity(bool rz,int amr)
{
    auto mesh=makeMesh(4,amr);Case s;s.pureFirstJet=true;s.radial=rz;s.planar=true;s.planeOffset=2;
    const auto geometry=CchMomentumGeometry::build(*mesh,metricsFor(*mesh,rz));const auto f=staticFrames(*mesh,geometry,s);const auto& l=geometry.layout();
    for(int chosen=0;chosen<l.velocityCount();++chosen){const auto a=geometry.metrics().faceMeasures()[chosen];if(a==0)continue;
        std::vector<double> values(l.velocityCount());values[chosen]=1/a;
        std::cout<<"BASIS_INPUT,"<<rz<<','<<amr<<','<<chosen<<','<<a<<','<<values[chosen]<<std::endl;
        const auto d=prepare(*mesh,f,CchNativeValues(l,values,std::vector<double>(l.pressureCount())));inspectTransfer(*mesh,d,&s);++basisChecks;
    }
    std::vector<double> a(l.velocityCount()),b=a,sum=a;
    for(int u=0;u<l.velocityCount();++u){const bool atAxis=geometry.metrics().faceMeasures()[u]==0;a[u]=atAxis?0:std::sin(1.3*u);b[u]=atAxis?0:std::cos(.7*u);sum[u]=a[u]+b[u];}
    const auto p=std::vector<double>(l.pressureCount());const auto da=prepare(*mesh,f,CchNativeValues(l,a,p)),db=prepare(*mesh,f,CchNativeValues(l,b,p)),ds=prepare(*mesh,f,CchNativeValues(l,sum,p));
    inspectTransfer(*mesh,ds,&s);for(std::size_t i=0;i<ds.faceRates().size();++i){norm(ds.faceRates()[i],da.faceRates()[i]+db.faceRates()[i],std::max(1.L,std::abs(da.faceRates()[i])+std::abs(db.faceRates()[i])),1e-12,"dual transfer linearity");++linearChecks;}
    std::cout<<"BASIS,"<<rz<<','<<amr<<','<<l.velocityCount()<<'\n';
}
void affine(bool rz,int amr,int selection)
{
    auto mesh=makeMesh(4,amr);Case s;s.pureFirstJet=true;s.radial=rz;s.planar=true;s.planeOffset=2;
    const auto geometry=CchMomentumGeometry::build(*mesh,metricsFor(*mesh,rz,selection));const auto f=staticFrames(*mesh,geometry,s);
    const auto field=[rz](Vector2D p){return Vector2D{rz?.3*p.x:.3+.2*p.x,-.1+.4*p.y};};
    const auto d=prepare(*mesh,f,velocity(*mesh,geometry.layout(),field));inspectTransfer(*mesh,d,&s);
    for(std::size_t i=0;i<d.faceRates().size();++i){const auto& g=d.geometry().segments()[i];const auto point=g.begin+(g.end-g.begin)*.5;const auto v=field(point);
        const auto expected=g.area*(g.normalAxis==0?v.x:v.y);peakAnalytic=std::max(peakAnalytic,norm(d.faceRates()[i],expected,std::max(1.L,std::abs(expected)),1e-12,"affine physical flux reproduction"));}
    if(selection!=0){
        const auto varying=velocity(*mesh,geometry.layout(),[rz](Vector2D p){return Vector2D{(rz?p.x:1)*(std::sin(3*p.x)+std::cos(2*p.y)),std::cos(5*p.x)-std::sin(p.y)};});
        inspectTransfer(*mesh,prepare(*mesh,f,varying),&s);std::cout<<"SELECTED_RANDOM,"<<rz<<','<<amr<<','<<selection<<'\n';
    }
    std::cout<<"AFFINE,"<<rz<<','<<amr<<','<<selection<<'\n';
}
long double liquid(Vector2D a,Vector2D b,Case s,double time)
{
    if(s.planar){const auto offset=s.planeOffset+s.speed*time;if(s.axial())b.y=std::min(b.y,offset);else b.x=std::min(b.x,offset);return a.x>=b.x || a.y>=b.y?0:cap(a,b,s.radial);}
    const auto center=Vector2D{s.radial?0:s.cx+(s.motion==1?s.speed*time:0),s.cy+(s.radial && s.motion==1?s.speed*time:0)};
    return circleCell(a,b,center,s.radius,s.radial);
}
std::vector<long double> constantRecovery(const QuadtreeMesh& mesh,const CchDualMaterial& d,const CchNativeValues& u,const Case* analytic=nullptr)
{
    const auto integrated=outward(d,true),impulse=d.consistencyImpulse(mesh,u,d.sources(),d.timeLayer());std::vector<long double> rhs;std::array<long double,2> kinetic{};
    for(int i=0;i<d.layout().velocityCount();++i){const long double v=u.normalVelocity()[i],ma=d.beginMass().masses()[i],mb=d.endMass().masses()[i];
        const auto corrected=ma*v-v*integrated[i]+impulse[i];rhs.push_back(corrected);
        peakConstant=std::max(peakConstant,norm(corrected/mb,v,std::max(1.L,std::abs(v)),1e-10,"no-force constant velocity recovery"));
        const auto omitted=std::abs((ma*v-v*integrated[i])/mb-v);omittedDrift=std::max(omittedDrift,omitted);
        if(analytic)actualOmittedDrift=std::max(actualOmittedDrift,omitted);
        kinetic[axis(d.layout(),i)]+=.5L*v*v*d.cells()[i].integratedDefect;
        if(analytic){long double va=0,vb=0,volume=0;for(const auto& box:intersections(mesh,d.layout(),i)){va+=liquid(box.lower,box.upper,*analytic,d.material().beginGeometry().timeLayer().time());vb+=liquid(box.lower,box.upper,*analytic,d.material().endGeometry().timeLayer().time());volume+=box.volume;}
            const auto change=(d.material().liquidDensity()-d.material().gasDensity())*(vb-va),r=change+integrated[i];
            const auto scale=(d.material().liquidDensity()+d.material().gasDensity())*volume;
            peakDefect=std::max({peakDefect,norm(d.cells()[i].massChange,change,scale,1e-10,"analytic dual geometry mass change"),norm(d.cells()[i].integratedDefect,r,scale+std::abs(integrated[i]),1e-10,"analytic nonconservative dual defect"),norm(impulse[i],v*r,scale*std::max(1.L,std::abs(v))+std::abs(integrated[i]),1e-10,"analytic consistency impulse")});}
    }
    const auto reported=d.materialKineticChange(mesh,u,d.sources(),d.timeLayer());for(int a=0;a<2;++a)norm(reported[a],kinetic[a],std::max(1.L,std::abs(kinetic[a])),1e-12,"material kinetic term");return rhs;
}
void saddle(const QuadtreeMesh& mesh,const CchDualMaterial& d,const CchNativeValues& u,const std::vector<long double>& corrected)
{
    const auto& l=d.layout();CchPressureLoad::Matrix a(l.velocityCount(),l.velocityCount());for(int i=0;i<l.velocityCount();++i)a.coeffRef(i,i)=d.endMass().masses()[i];
    const auto pressure=CchPressureLoad::build(d.geometry().metrics().divergence());const auto solver=CchSinglePhaseSaddle::buildClosed(pressure,CchSaddleVelocityBlock::build(l,a));
    std::vector<double> p(l.pressureCount());long double mean=0,total=0;for(int c=0;c<l.pressureCount();++c){p[c]=mesh.cells(mesh.cellIndexToLabel[l.nativeCell(c)]).centerCoord().y;const double volume=d.geometry().metrics().cellVolumes()[c];mean+=volume*p[c];total+=volume;}
    for(auto& value:p)value-=static_cast<double>(mean/total);
    const auto load=pressure.apply(mesh,CchNativeValues(l,std::vector<double>(l.velocityCount()),p));std::vector<double> rhs;for(int i=0;i<l.velocityCount();++i)rhs.push_back(static_cast<double>(corrected[i]+load[i]));
    auto prescribed=u.normalVelocity();for(int i:solver.freeVelocityDofs())prescribed[i]=0;
    const auto result=solver.solve(mesh,CchSaddleRhs(l,rhs,prescribed));
    for(int i=0;i<l.velocityCount();++i)peakSaddleU=std::max(peakSaddleU,norm(result.values().normalVelocity()[i],u.normalVelocity()[i],1,1e-10,"true material Mu saddle velocity"));
    for(int c=0;c<l.pressureCount();++c)peakSaddleP=std::max(peakSaddleP,norm(result.values().pressure()[c],p[c],1,1e-8,"true material Mu manufactured pressure"));++saddles;
}
void actualCrossing(bool rz,int amr,bool planar)
{
    auto mesh=makeMesh(4,amr);Case s;s.pureFirstJet=true;s.radial=rz;s.planar=planar;s.motion=1;s.speed=.4;s.planeOffset=.49;s.radius=rz?.34:.19;
    const auto geometry=CchMomentumGeometry::build(*mesh,metricsFor(*mesh,rz));const auto a=initial(*mesh,s);std::string stage;
    const auto m=advanceFirstJet(*mesh,a,s,{1,0,.1,Position::Midpoint},stage),b=advanceFirstJet(*mesh,a,s,{1,0,.1,Position::End},stage);
    const auto u=velocity(*mesh,geometry.layout(),[rz](Vector2D){return rz?Vector2D{0,.4}:Vector2D{.4,0};});
    for(double rho:{1.,1000.,1e6}){const auto d=prepare(*mesh,frames(*mesh,geometry,a,m,b,rho),u);inspectTransfer(*mesh,d,&s);const auto rhs=constantRecovery(*mesh,d,u,&s);
        for(const auto& c:d.cells())actualTransferChange=std::max(actualTransferChange,std::abs(c.geometryTransferChange));
        if(rho>1)check(d.report().components[0].l1IntegratedDefect+d.report().components[1].l1IntegratedDefect>0,"actual geometry defect filled zero");
        if(planar)saddle(*mesh,d,u,rhs);++crossings;std::cout<<"CROSSING,"<<rz<<','<<amr<<','<<planar<<','<<rho<<','<<d.report().components[0].l1IntegratedDefect+d.report().components[1].l1IntegratedDefect<<'\n';}
}
void syntheticConstant(bool rz,int amr)
{
    auto mesh=makeMesh(4,amr);Case s;s.pureFirstJet=true;s.radial=rz;s.planar=true;s.planeOffset=.43;
    const auto geometry=CchMomentumGeometry::build(*mesh,metricsFor(*mesh,rz));const auto f=staticFrames(*mesh,geometry,s,1e6);
    const auto u=velocity(*mesh,geometry.layout(),[rz](Vector2D){return rz?Vector2D{0,.27}:Vector2D{.31,-.27};});
    const auto d=prepare(*mesh,f,u);inspectTransfer(*mesh,d,&s);constantRecovery(*mesh,d,u);
    check(d.report().components[0].l1IntegratedDefect+d.report().components[1].l1IntegratedDefect>0,"synthetic continuity source disappeared");
    std::cout<<"SYNTHETIC,"<<rz<<','<<amr<<'\n';
}
void temporal(bool rz,int amr)
{
    auto mesh=makeMesh(4,amr);Case s;s.pureFirstJet=true;s.radial=rz;s.planar=true;s.planeOffset=2;
    const auto geometry=CchMomentumGeometry::build(*mesh,metricsFor(*mesh,rz));constexpr double total=.3;double previous=0,order=0;
    for(int n:{4,8,16}){std::vector<long double> accumulated(geometry.segments().size());auto begin=publishAt(*mesh,s,{});
        for(int k=0;k<n;++k){const double t0=total*k/n,t1=total*(k+1)/n,tm=(t0+t1)/2;
            const auto mid=publishAt(*mesh,s,{k+1,t0,t1,Position::Midpoint}),end=publishAt(*mesh,s,{k+1,t0,t1,Position::End});
            const auto u=velocity(*mesh,geometry.layout(),[tm](Vector2D p){return Vector2D{std::sin(tm)*p.x,std::cos(tm)*p.y};});
            const auto d=prepare(*mesh,frames(*mesh,geometry,begin,mid,end,1),u);for(std::size_t i=0;i<accumulated.size();++i)accumulated[i]+=d.faceIntegrals()[i];begin=end;
        }
        long double error=0,scale=0;for(std::size_t i=0;i<accumulated.size();++i){const auto& g=geometry.segments()[i];const auto mid=g.begin+(g.end-g.begin)*.5;
            const auto expected=g.area*(g.normalAxis==0?mid.x*(1-std::cos(total)):mid.y*std::sin(total));error+=std::abs(accumulated[i]-expected);scale+=std::abs(expected);}
        const auto e=static_cast<double>(error/scale);check(e>0 && std::isfinite(e),"invalid temporal reference error");if(n>4){check(e<previous,"dual midpoint error not decreasing");order=std::log2(previous/e);}previous=e;++timeRows;
        std::cout<<"TIME,"<<rz<<','<<amr<<','<<n<<','<<e<<'\n';
    }
    check(order>=1.8,"dual smooth midpoint order below gate");std::cout<<"TIME_ORDER,"<<rz<<','<<amr<<','<<order<<'\n';
}
template<class F>void reject(F f,const std::string& why){rejects(f,why);++rejected;}
void guards()
{
    auto mesh=makeMesh(4);Case s;s.pureFirstJet=true;s.planar=true;s.planeOffset=.43;
    const auto geometry=CchMomentumGeometry::build(*mesh,metricsFor(*mesh,false));const auto f=staticFrames(*mesh,geometry,s,1000);
    const auto u=velocity(*mesh,geometry.layout(),[](Vector2D p){return Vector2D{.1*p.x,.2*p.y};});const auto d=prepare(*mesh,f,u);const auto saved=d.faceRates();
    reject([&]{CchDualMaterial::prepare(*mesh,d.material(),f.ma,f.mm,CchMomentumMass::build(*mesh,geometry,f.b,1,1));},"different material densities");
    const auto other=CchMomentumGeometry::build(*mesh,geometry.metrics());
    reject([&]{CchDualMaterial::prepare(*mesh,d.material(),f.ma,CchMomentumMass::build(*mesh,other,f.m,1000,1),f.mb);},"different dual geometry identity");
    reject([&]{CchDualMaterial::prepare(*mesh,d.material(),CchMomentumMass::build(*mesh,other,f.a,1000,1),CchMomentumMass::build(*mesh,other,f.m,1000,1),CchMomentumMass::build(*mesh,other,f.b,1000,1));},"different momentum geometry instance");
    const auto foreign=staticFrames(*mesh,geometry,s,1000);
    reject([&]{CchDualMaterial::prepare(*mesh,d.material(),foreign.ma,f.mm,f.mb);},"different geometry identity");
    reject([&]{CchDualMaterial::prepare(*mesh,d.material(),f.ma,f.mb,f.mm);},"different geometry identity");
    const auto differentLayout=CchMomentumGeometry::build(*mesh,metricsFor(*mesh,false));
    const auto layoutA=CchPhaseGeometry::build(*mesh,differentLayout.metrics(),f.a.publication()),layoutM=CchPhaseGeometry::build(*mesh,differentLayout.metrics(),f.m.publication()),layoutB=CchPhaseGeometry::build(*mesh,differentLayout.metrics(),f.b.publication());
    reject([&]{CchDualMaterial::prepare(*mesh,d.material(),CchMomentumMass::build(*mesh,differentLayout,layoutA,1000,1),CchMomentumMass::build(*mesh,differentLayout,layoutM,1000,1),CchMomentumMass::build(*mesh,differentLayout,layoutB,1000,1));},"different native layout");
    const auto newDi=CchMomentumGeometry::build(*mesh,CchPhysicalDivergence::build(*mesh,CchIntegralDivergence::build(*mesh,geometry.layout())));
    const auto diA=CchPhaseGeometry::build(*mesh,newDi.metrics(),f.a.publication()),diM=CchPhaseGeometry::build(*mesh,newDi.metrics(),f.m.publication()),diB=CchPhaseGeometry::build(*mesh,newDi.metrics(),f.b.publication());
    reject([&]{CchDualMaterial::prepare(*mesh,d.material(),CchMomentumMass::build(*mesh,newDi,diA,1000,1),CchMomentumMass::build(*mesh,newDi,diM,1000,1),CchMomentumMass::build(*mesh,newDi,diB,1000,1));},"different retained D_I");
    for(int which=0;which<3;++which){auto source=d.sources();if(which==0)source.begin={};if(which==1)source.midpoint={};if(which==2)source.end={};reject([&]{d.validateFor(*mesh,source,d.timeLayer());},"different geometry identity");}
    auto wrongAdvector=d.sources();wrongAdvector.advector={};reject([&]{d.validateFor(*mesh,wrongAdvector,d.timeLayer());},"different actual advector identity");
    auto clock=d.timeLayer();clock.step+=1;reject([&]{d.consistencyImpulse(*mesh,u,d.sources(),clock);},"different full time layer");
    clock=d.timeLayer();clock.begin=-.01;reject([&]{d.materialKineticChange(*mesh,u,d.sources(),clock);},"different full time layer");
    const auto fresh=metricsFor(*mesh,false).layout();reject([&]{d.consistencyImpulse(*mesh,velocity(*mesh,fresh,[](Vector2D){return Vector2D{1,0};}),d.sources(),d.timeLayer());},"different layout identity");
    const auto huge=velocity(*mesh,geometry.layout(),[](Vector2D){return Vector2D{std::numeric_limits<double>::max(),std::numeric_limits<double>::max()};});
    const auto largeDefect=prepare(*mesh,f,velocity(*mesh,geometry.layout(),[](Vector2D p){return Vector2D{100*p.x,200*p.y};}));
    reject([&]{largeDefect.consistencyImpulse(*mesh,huge,largeDefect.sources(),largeDefect.timeLayer());},"nonfinite consistency impulse");
    const auto kineticHuge=velocity(*mesh,geometry.layout(),[](Vector2D){return Vector2D{1e155,1e155};});
    reject([&]{d.materialKineticChange(*mesh,kineticHuge,d.sources(),d.timeLayer());},"nonfinite material kinetic term");
    if(std::numeric_limits<long double>::digits==std::numeric_limits<double>::digits){
        int chosen=-1;for(int i=0;i<d.layout().velocityCount();++i)if(d.cells()[i].integratedDefect!=0 && std::abs(d.cells()[i].integratedDefect)<.25L){chosen=i;break;}
        check(chosen>=0,"missing underflow source fixture");std::vector<double> tiny(d.layout().velocityCount());tiny[chosen]=std::numeric_limits<double>::denorm_min();
        reject([&]{d.consistencyImpulse(*mesh,CchNativeValues(d.layout(),tiny,std::vector<double>(d.layout().pressureCount())),d.sources(),d.timeLayer());},"unrepresentable consistency impulse");
    }
    Case radial=s;radial.radial=true;const auto radialGeometry=CchMomentumGeometry::build(*mesh,metricsFor(*mesh,true));const auto radialFrames=staticFrames(*mesh,radialGeometry,radial,1000);
    const auto radialU=velocity(*mesh,radialGeometry.layout(),[](Vector2D){return Vector2D{0,.1};});const auto radialD=prepare(*mesh,radialFrames,radialU);
    auto wrongAxis=radialU.normalVelocity();wrongAxis[radialGeometry.metrics().axisVelocityDofs().front()]=1;
    reject([&]{radialD.consistencyImpulse(*mesh,CchNativeValues(radialGeometry.layout(),wrongAxis,std::vector<double>(radialGeometry.layout().pressureCount())),radialD.sources(),radialD.timeLayer());},"nonzero radial axis velocity");
    check(d.faceRates()==saved,"rejection changed retained unique flux");const auto copy=d;check(&copy.faceRates()==&d.faceRates(),"prepared copy detached data");
    auto equivalent=makeMesh(4);mesh.reset();copy.validateFor(*equivalent,copy.sources(),copy.timeLayer());inspectTransfer(*equivalent,copy,&s);
    auto& vertex=equivalent->verts(equivalent->vertIndexToLabel[0]);const auto old=vertex.coord();vertex.setCoord(old+Vector2D{.01,0});reject([&]{copy.validateFor(*equivalent,copy.sources(),copy.timeLayer());},"NativeMeshIdentity");vertex.setCoord(old);
    DictionaryMesh::RefineBox box;box.level=1;box.min={.25,.25};box.max={.75,.75};equivalent->applyRefineBoxes({box});reject([&]{copy.validateFor(*equivalent,copy.sources(),copy.timeLayer());},"NativeMeshIdentity");
}
void regrid(bool rz)
{
    auto oldMesh=makeMesh(4),fine=makeMesh(4,1);Case s;s.pureFirstJet=true;s.radial=rz;s.planar=true;s.planeOffset=.49;s.motion=1;s.speed=.4;
    const auto old=initial(*oldMesh,s);HermiteGalsPublicationConfig cfg;cfg.volume.symmetryPlanes=s.symmetry(*fine);
    const auto transferred=transferHermiteGalsFirstJetMesh(*oldMesh,old,*fine,HermiteGalsMeshTransferKind::ExactRefinement,cfg);
    check(transferred.sourceIdentity()==old.state().identity() && transferred.publication().state().timeLayer().time()==old.state().timeLayer().time(),"refinement lost physical origin/clock");
    const auto advance=[&](const QuadtreeMesh& mesh,const HermiteGalsPublication& start,HermiteGalsTimeLayer mid,HermiteGalsTimeLayer end){
        const auto geometry=CchMomentumGeometry::build(mesh,metricsFor(mesh,rz));std::string stage;
        const auto m=advanceFirstJet(mesh,start,s,mid,stage),b=advanceFirstJet(mesh,start,s,end,stage);
        const auto u=velocity(mesh,geometry.layout(),[rz](Vector2D){return rz?Vector2D{0,.4}:Vector2D{.4,0};});
        const auto d=prepare(mesh,frames(mesh,geometry,start,m,b,1e6),u);inspectTransfer(mesh,d,&s);constantRecovery(mesh,d,u,&s);return d;
    };
    const auto refined=advance(*fine,transferred.publication(),{1,0,.1,Position::Midpoint},{1,0,.1,Position::End});
    auto coarse=makeMesh(4);cfg.volume.symmetryPlanes=s.symmetry(*coarse);
    const auto projected=transferHermiteGalsFirstJetMesh(*fine,refined.endMass().phaseGeometry().publication(),*coarse,HermiteGalsMeshTransferKind::FirstJetProjection,cfg);
    check(projected.publication().state().timeLayer().time()==.1,"coarsening changed physical time");
    reject([&]{refined.validateFor(*coarse,refined.sources(),refined.timeLayer());},"NativeMeshIdentity");
    const auto continued=advance(*coarse,projected.publication(),{2,.1,.2,Position::Midpoint},{2,.1,.2,Position::End});
    check(continued.sources().begin==projected.publication().state().identity() && continued.timeLayer().step==2,"regrid silently retagged old interval");regrids+=2;
    std::cout<<"REGRID,"<<rz<<','<<transferred.report().relativeDomainVolumeChange<<','<<projected.report().relativeDomainVolumeChange<<'\n';
}
void spatialDual()
{
    int scenes=0,counterflow=0,tiny=0;std::size_t rows=0;long double peak=0,nonzeroExchange=0;
    for(bool rz:{false,true})for(int amr:{0,1}){auto mesh=makeMesh(4,amr);const auto geometry=CchMomentumGeometry::build(*mesh,metricsFor(*mesh,rz));const auto reconstruction=CchAdvectorReconstruction::prepare(*mesh,geometry);
        for(int kind=0;kind<10;++kind){Case s;s.pureFirstJet=true;s.radial=rz;s.planar=kind!=0;s.radius=rz?.34:.19;
            s.planeOffset=kind==1?.43:kind==2?2:kind==3?-1:kind==4?std::ldexp(1.,-40):kind==5?1-std::ldexp(1.,-40):kind<8?(amr?.125:.375):.5;s.planeSign=kind==7 || kind==9?-1:1;
            const auto f=staticFrames(*mesh,geometry,s);for(int mode=0;mode<4;++mode){if(kind<6 && mode!=2)continue;
                const auto actual=reconstruction.evaluate(*mesh,mmsValues(*mesh,geometry.layout(),mode),f.m.geometryIdentity(),f.m.timeLayer());const auto phase=CchPhaseRate::prepareSpatial(*mesh,f.m,actual);
                for(double rho:{1.,1000.,1e6}){const double l=kind==5?1:rho,g=kind==5?rho:1;
                    const auto flux=CchMaterialFlux::prepare(*mesh,phase,l,g);const auto transport=CchMaterialTransport::prepare(*mesh,f.a,flux,f.b);
                    const auto d=CchDualMaterial::prepare(*mesh,transport,CchMomentumMass::build(*mesh,geometry,f.a,l,g),CchMomentumMass::build(*mesh,geometry,f.m,l,g),CchMomentumMass::build(*mesh,geometry,f.b,l,g));
                    check(&d.faceRates()==&flux.dualFaceRates(),"H09 copied/recomputed single F");inspectTransfer(*mesh,d);
                    std::vector<long double> expected;for(std::size_t i=0;i<reconstruction.atoms().size();++i){const auto& a=reconstruction.atoms()[i];const auto q=referenceAtom(a,mmsCoefficients(a,rz,mode),s);expected.push_back(l*q.liquid+g*q.gas);++rows;
                        peak=std::max(peak,norm(flux.atomRates()[i],expected.back(),(l+g)*a.area,1e-12,"independent declared MMS material atom"));
                        if(rho>1 && (kind==4 || kind==5)){const auto small=kind==4?q.liquid:q.gas;if(small!=0 && std::abs(small)<1e-10L*a.area){const auto measured=kind==4?phase.atomRates()[i].liquid:phase.atomRates()[i].gas;
                            norm(measured,small,std::abs(small),1e-12,"dual tiny phase own scale");check(flux.atomRates()[i]!=(kind==4?g*phase.atomRates()[i].gas:l*phase.atomRates()[i].liquid),"dual tiny phase erased");++tiny;}}}
                    for(std::size_t id=0;id<d.faceRates().size();++id){long double sum=0;for(int a:reconstruction.dualAtoms()[id])sum+=expected[a];peak=std::max(peak,norm(d.faceRates()[id],sum,(l+g)*geometry.segments()[id].area,1e-12,"independent declared MMS dual F"));
                        const auto q=phase.dualFaceRates()[id];if(rho>1 && std::abs(q.liquid+q.gas)<1e-14L*geometry.segments()[id].area && q.liquid!=0 && q.gas!=0){check(d.faceRates()[id]!=0,"zero-net dual counterflow erased");++counterflow;}}
                    nonzeroExchange=std::max(nonzeroExchange,d.report().maxAbsSpatialExchange);++scenes;
                }
            }
        }
    }
    check(counterflow>0 && tiny>0 && nonzeroExchange>1,"spatial dual controls absent");
    std::cout<<"SPATIAL_DUAL_SUMMARY,"<<scenes<<','<<rows<<','<<counterflow<<','<<tiny<<','<<peak<<','<<nonzeroExchange<<",PASS\n";
}
template<class F>void migrationCheck(const std::string& name,F f)
{try{f();}catch(const std::exception& e){++migrationFailures;std::cout<<"MIGRATION_FAILURE,"<<name<<','<<e.what()<<std::endl;}}
// Closed power integrals are independent of the GL3 reference kernel.
template<std::size_t N> long double powerIntegral(const std::array<long double,N>& c,long double a,long double b,bool rz,long double r,int slope)
{
    long double out=0;for(std::size_t k=0;k<N;++k){const auto m=(std::pow(b,static_cast<int>(k+1))-std::pow(a,static_cast<int>(k+1)))/(k+1);
        const auto rm=(std::pow(b,static_cast<int>(k+2))-std::pow(a,static_cast<int>(k+2)))/(k+2);out+=c[k]*(rz?2*pi*(r*m+slope*rm):m);}return out;
}
template<std::size_t N> void polynomialRangeRows(const std::array<long double,N>& begin,int& rows,long double& peak)
{
    for(int metric=0;metric<3;++metric)for(bool fromEnd:{false,true})for(int sign:{-1,1})for(const auto& range:std::array<std::array<long double,2>,3>{{{{0,1}},{{.125L,.25L}},{{.875L,.125L}}}}){
        const bool rz=metric>0;const int physicalSlope=metric==2?1:0;const long double initialRadius=metric==1?2:.25L;
        const auto c=fromEnd?referenceShift(begin,1,-1):begin;const auto a=range[0],width=range[1];const auto offset=fromEnd?1-a-width:a;
        const auto radius=initialRadius+(fromEnd?physicalSlope:0);const int slope=fromEnd?-physicalSlope:physicalSlope;
        const auto result=referenceShiftedIntegral(c,offset,width,rz,radius,slope,sign);
        const auto expected=sign*powerIntegral(begin,a,a+width,rz,initialRadius,physicalSlope);
        const auto work=width*referenceWork(begin,1)*(rz?2*pi*(initialRadius+physicalSlope):1.L);
        const auto error=std::abs(result-expected)/work;check(error<=1e-12,"complete-degree shifted GL3 moment");peak=std::max(peak,error);++rows;
        std::cout<<"POLYNOMIAL_RANGE,"<<N<<','<<metric<<','<<fromEnd<<','<<sign<<','<<offset<<','<<width<<','<<radius<<','<<slope;
        for(int k=0;k<5;++k)std::cout<<','<<(k<N?c[k]:0);std::cout<<','<<result<<','<<expected<<','<<work<<'\n';
    }
}
void referencePolynomialContract()
{
    const auto started=std::chrono::steady_clock::now();int ranges=0,tiny=0,atoms=0,legacy=0;long double peak=0;
    polynomialRangeRows(std::array<long double,1>{1},ranges,peak);
    polynomialRangeRows(std::array<long double,2>{0,1},ranges,peak);
    polynomialRangeRows(std::array<long double,3>{0,0,1},ranges,peak);
    polynomialRangeRows(std::array<long double,4>{0,0,0,1},ranges,peak);
    polynomialRangeRows(std::array<long double,5>{0,0,0,0,1},ranges,peak);
    polynomialRangeRows(std::array<long double,5>{1,-2,3,-4,5},ranges,peak);
    // A fourth-order zero at offset=.5 leaves O(width^5/6) moments. Forming
    // global endpoints or sampling offset+t would erase every contribution.
    const std::array<long double,5> root{{.0625L,-.5L,1.5L,-2,1}};
    for(int metric=0;metric<3;++metric)for(bool fromEnd:{false,true})for(int sign:{-1,1})for(int exponent:{60,80}){
        const bool rz=metric>0;const int physicalSlope=metric==2?1:0,slope=fromEnd?-physicalSlope:physicalSlope;const long double width=std::ldexp(1.L,-exponent),offset=.5L;
        check(static_cast<double>(offset)+static_cast<double>(width)==static_cast<double>(offset),"tiny range must lose its double global endpoint");
        const long double origin=metric==1?2:.25L+(fromEnd?physicalSlope:0),radius=origin+slope*offset;
        const auto c=fromEnd?referenceShift(root,1,-1):root;const auto result=referenceShiftedIntegral(c,offset,width,rz,origin,slope,sign);
        const auto expected=sign*std::pow(width,5)*(rz?2*pi*(radius/5+slope*width/6):1.L/5);
        const auto e=std::abs((result-expected)/expected);check(result!=0 && e<=1e-12,"tiny quartic moment erased");peak=std::max(peak,e);++tiny;
        std::cout<<"POLYNOMIAL_TINY,5,"<<metric<<','<<fromEnd<<','<<sign<<','<<offset<<','<<width<<','<<origin<<','<<slope;
        for(auto x:c)std::cout<<','<<x;std::cout<<','<<result<<','<<expected<<','<<std::abs(expected)<<'\n';
    }
    // Known factors, not referenceRoots, define exact phase subintervals.
    const std::array<std::array<long double,4>,3> roots{{{{.125L,.375L,.625L,.875L}},{{.25L,.25L,.75L,.75L}},{{0,.25L,.75L,1}}}};
    for(int fixture=0;fixture<3;++fixture)for(bool rz:{false,true})for(int n:{0,1})for(int velocitySign:{-1,1})for(int phaseSign:{-1,1})for(int normalSign:{-1,1}){
        std::array<long double,5> c{};c[0]=velocitySign;int degree=0;
        for(auto r:roots[fixture]){for(int k=++degree;k>0;--k)c[k]=c[k-1]-r*c[k];c[0]*=-r;}
        CchAdvectorAtom atom;atom.normalAxis=n;atom.begin=n==0?Vector2D{.5,0}:Vector2D{0,.5};atom.end=n==0?Vector2D{.5,1}:Vector2D{1,.5};
        Case s;s.pureFirstJet=true;s.radial=rz;s.planar=true;s.planeOffset=.5;s.planeSign=phaseSign;
        // The fixture's planar direction follows its physical convention:
        // XY normal x, RZ normal z. Both normal axes also exercise true cuts.
        const bool coincident=n==(s.axial()?1:0);const auto actual=referenceAtom(atom,c,s,normalSign);CchPhaseMeasure expected;
        std::vector<long double> cuts{0,1};if(coincident)for(auto r:roots[fixture])cuts.push_back(r);else cuts.push_back(.5L);
        std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
        for(std::size_t k=1;k<cuts.size();++k){const auto a=cuts[k-1],b=cuts[k],mid=(a+b)/2;long double normal=velocitySign;for(auto r:roots[fixture])normal*=mid-r;
            const bool liquid=coincident?phaseSign*normal>0:phaseSign*(mid-.5L)<0;const auto integral=normalSign*powerIntegral(c,a,b,rz,n==0?.5L:0,n==0?0:1);
            (liquid?expected.liquid:expected.gas)+=integral;}
        const auto work=referenceWork(c,1)*(rz?2*pi:1.L);const auto e=std::max(std::abs(actual.liquid-expected.liquid),std::abs(actual.gas-expected.gas))/work;
        check(e<=1e-12,"quartic atom causal/physical phase split");peak=std::max(peak,e);++atoms;
        std::cout<<"POLYNOMIAL_ATOM,"<<fixture<<','<<rz<<','<<n<<','<<velocitySign<<','<<phaseSign<<','<<normalSign<<','<<coincident<<','<<actual.liquid<<','<<actual.gas<<','<<expected.liquid<<','<<expected.gas<<','<<work<<'\n';
    }
    for(int metric=0;metric<3;++metric)for(bool fromEnd:{false,true})for(int sign:{-1,1})for(const auto& range:std::array<std::array<long double,2>,2>{{{{0,1}},{{.125L,.25L}}}}){
        CchNormalPolynomial p;p.beginCoefficients={1,-2,3};p.endCoefficients={2,-4,3};p.length=1;p.radialMeasure=metric>0;p.radialTangent=metric==2;p.radiusBegin=metric==1?2:.25L;p.radiusEnd=p.radiusBegin+(p.radialTangent?1:0);
        const auto& c=fromEnd?p.endCoefficients:p.beginCoefficients;const auto actual=referenceShiftedIntegral(c,range[0],range[1],p.radialMeasure,fromEnd?p.radiusEnd:p.radiusBegin,p.radialTangent?(fromEnd?-1:1):0,sign);
        const auto expected=p.integral(fromEnd,range[0],range[1],sign);const auto error=std::abs((actual-expected)/expected);check(error<=1e-12,"quadratic reference/production regression");peak=std::max(peak,error);++legacy;
    }
    check(ranges==216 && tiny==24 && atoms==96 && legacy==24,"complete-degree reference inventory");
    std::cout<<"POLYNOMIAL_REFERENCE_SUMMARY,"<<ranges<<','<<tiny<<','<<atoms<<','<<legacy<<','<<peak<<",PASS\n";
    std::cout<<"POLYNOMIAL_REFERENCE_SECONDS,"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()<<'\n';
}

}
int main()
{
    try{
        std::cout<<std::setprecision(17);
        for(bool rz:{false,true})for(int amr:{0,1}){const auto name=std::to_string(rz)+","+std::to_string(amr);
            migrationCheck("basis,"+name,[&]{basisAndLinearity(rz,amr);});for(int selection:{0,1,2,3})migrationCheck("affine,"+name+","+std::to_string(selection),[&]{affine(rz,amr,selection);});
            for(bool planar:{false,true})migrationCheck("crossing,"+name+","+std::to_string(planar),[&]{actualCrossing(rz,amr,planar);});migrationCheck("synthetic,"+name,[&]{syntheticConstant(rz,amr);});migrationCheck("temporal,"+name,[&]{temporal(rz,amr);});}
        check(actualOmittedDrift>1e-5 && actualTransferChange>1e-5,"actual crossing controls did not expose missing consistency/geometry-transfer term");
        migrationCheck("guards",[&]{guards();});for(bool rz:{false,true})migrationCheck("regrid,"+std::to_string(rz),[&]{regrid(rz);});migrationCheck("spatialDual",[&]{spatialDual();});
        std::cout<<"CONTROL,"<<actualOmittedDrift<<','<<actualTransferChange<<','<<omittedDrift<<','<<transferChange<<'\n';
        std::cout<<"SUMMARY,"<<scenarios<<','<<segmentChecks<<','<<dualChecks<<','<<basisChecks<<','<<linearChecks<<','<<crossings<<','<<saddles<<','<<rejected<<','<<timeRows<<','<<regrids<<','<<peakExchange<<','<<peakAnalytic<<','<<peakDefect<<','<<peakDecomposition<<','<<peakConstant<<','<<peakSaddleU<<','<<peakSaddleP<<'\n';
        std::cout<<"MIGRATION_RESULT,"<<migrationFailures<<','<<(migrationFailures?"FAIL":"PASS")<<'\n';
        if(migrationFailures)return 1;
        std::cout<<"H09b CCH unique dual flux and nonconservative momentum consistency PASS\n";
        referencePolynomialContract();
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
