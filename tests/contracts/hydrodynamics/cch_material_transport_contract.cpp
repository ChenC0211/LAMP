#include "Modules/Hydrodynamics/CchMaterialTransport.h" // Public header stands alone.
#include "Modules/Hydrodynamics/CchKepConvection.h"
#include "interface_first_jet_fixture.h"
#include <chrono>
#include "Core/Mesh/CartesianMeasure.h"
#include <map>

namespace
{
using namespace Quantum::Hydrodynamics;
int scenarios=0,cellChecks=0,faceChecks=0,internalChecks=0,coarseChecks=0,rejected=0,spatialRows=0,temporalRows=0;
long double peakMass=0,peakFlux=0,peakDefect=0,peakBalance=0,peakCoarse=0;
static_assert(!std::is_default_constructible_v<CchMaterialTransport>);
static_assert(!std::is_default_constructible_v<CchMaterialFlux>);
static_assert(std::is_same_v<decltype(std::declval<CchMaterialFlux&>().atomRates()),const std::vector<long double>&>);
static_assert(std::is_same_v<decltype(std::declval<CchMaterialTransport&>().cells()),const std::vector<CchMaterialCell>&>);
CchPhysicalDivergence metricsFor(const QuadtreeMesh& mesh,bool rz,bool subset=false)
{
    FaceSegmentGraph graph;graph.build(mesh,[subset](const QuadtreeMesh&,const Cell& c){return !subset || c.centerCoord().x<.5;});
    return CchPhysicalDivergence::build(mesh,CchIntegralDivergence::build(mesh,CchDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,rz))));
}
HermiteGalsPublication publishAt(const QuadtreeMesh& mesh,const Case& s,HermiteGalsTimeLayer time)
{
    const auto f=s.profile(time.time());
    const auto seed=stateFor(mesh,f,s.radial,time,[](int,Vector2D,Vector2D){return 0.;},s.symmetry(mesh),HermiteGalsTransport::FirstJet);
    HermiteGalsCandidate candidate(mesh,time,seed.primary());const auto mixed=closeHermiteGalsMixedAmr(mesh,candidate,s.symmetry(mesh));
    HermiteGalsPublicationConfig config;config.volume.symmetryPlanes=s.symmetry(mesh);
    const auto predicted=centers(mesh,f);
    return publishHermiteGalsFirstJetGeometry(mesh,std::move(candidate),mixed,{candidate.identity(),predicted},config);
}
CchNativeValues velocity(const QuadtreeMesh& mesh,const CchDofLayout& layout,Vector2D v)
{
    std::vector<double> normal;
    for(int f=0;f<layout.velocityCount();++f){const auto face=layout.domain().face(mesh,layout.nativeFace(f));normal.push_back(face.horizontal?v.y:v.x*(layout.domain().axisymmetric()?face.center.x:1));}
    return CchNativeValues(layout,normal,std::vector<double>(layout.pressureCount()));
}
CchPhaseRate rate(const QuadtreeMesh& mesh,CchPhaseGeometry g,const CchNativeValues& u)
{
    const auto geometry=CchMomentumGeometry::build(mesh,g.metrics());
    const auto advector=CchAdvectorReconstruction::prepare(mesh,geometry).evaluate(mesh,u,g.geometryIdentity(),g.timeLayer());
    return CchPhaseRate::prepareSpatial(mesh,g,advector);
}
Case shapeAt(Case s,double time)
{
    if(s.motion==1){if(s.planar)s.planeOffset+=s.speed*time;else if(s.axial())s.cy+=s.speed*time;else s.cx+=s.speed*time;}
    return s;
}
long double capacity(Vector2D a,Vector2D b,bool rz)
{return (b.x-a.x)*(b.y-a.y)*(rz?CartesianMeasure::twoPi*(a.x+b.x)/2:1.L);}
long double liquidCell(Vector2D a,Vector2D b,const Case& s)
{
    if(!s.planar)return circleCell(a,b,{s.radial?0:s.cx,s.cy},s.radius,s.radial);
    if(s.axial()){if(s.planeSign>0)b.y=std::min(b.y,s.planeOffset);else a.y=std::max(a.y,s.planeOffset);}
    else {if(s.planeSign>0)b.x=std::min(b.x,s.planeOffset);else a.x=std::max(a.x,s.planeOffset);}
    if(b.x<=a.x || b.y<=a.y)return 0;return capacity(a,b,s.radial);
}
long double liquidInterval(double lo,double hi,double normal,bool horizontal,const Case& s)
{
    if(s.planar){
        if(horizontal!=s.axial()){if(s.planeSign>0)hi=std::min(hi,s.planeOffset);else lo=std::max(lo,s.planeOffset);}
        else if(s.planeSign*(normal-s.planeOffset)>=0)return 0;
    }else {
        const double nc=horizontal?s.cy:(s.radial?0:s.cx),tc=horizontal?(s.radial?0:s.cx):s.cy;
        const long double d=static_cast<long double>(normal)-nc,q=static_cast<long double>(s.radius)*s.radius-d*d;
        if(q<=0)return 0;const double width=static_cast<double>(std::sqrt(q));lo=std::max(lo,tc-width);hi=std::min(hi,tc+width);
    }
    if(hi<=lo)return 0;
    return (hi-lo)*(s.radial?CartesianMeasure::twoPi*(horizontal?(lo+hi)/2:normal):1.L);
}
std::array<Vector2D,2> faceBounds(const QuadtreeMesh& mesh,const NativeLeafFace& face)
{
    const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
    auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
    if(face.horizontal?a.x>b.x:a.y>b.y)std::swap(a,b);return {{a,b}};
}
long double referenceFace(const QuadtreeMesh& mesh,const NativeLeafFace& face,const Case& s)
{
    const auto b=faceBounds(mesh,face);return liquidInterval(face.horizontal?b[0].x:b[0].y,face.horizontal?b[1].x:b[1].y,face.horizontal?b[0].y:b[0].x,face.horizontal,s);
}
long double faceCapacity(const QuadtreeMesh& mesh,const NativeLeafFace& face,bool rz)
{
    const auto b=faceBounds(mesh,face);const double width=face.horizontal?b[1].x-b[0].x:b[1].y-b[0].y;
    return width*(rz?CartesianMeasure::twoPi*(b[0].x+b[1].x)/2:1.L);
}
long double normalized(long double a,long double b,long double scale,long double tolerance,const std::string& why)
{const auto e=std::abs(a-b)/(scale>0?scale:1);check(std::isfinite(e) && e<=tolerance,why+" error="+std::to_string(static_cast<double>(e)));return e;}
struct Error{long double mass=0,flux=0;};
Error inspect(const QuadtreeMesh& mesh,const Case& s,const CchMaterialTransport& material,const CchNativeValues& u,bool exact=true)
{
    const auto& layout=material.layout();const double l=material.liquidDensity(),g=material.gasDensity(),dt=material.deltaTime();
    check(material.isSpatial() && material.sources().advector==material.midpointRate().advectorIdentity(),"material fixture/source not migrated");
    const auto middle=shapeAt(s,material.timeLayer().time());
    std::vector<long double> outward(layout.pressureCount());
    long double boundary=0,sumAbs=0,totalVolume=0,totalArea=0;Error error;
    std::map<std::pair<int,int>,std::vector<int>> coarse;
    for(int f=0;f<layout.velocityCount();++f){
        const auto face=layout.domain().face(mesh,layout.nativeFace(f));
        const long double area=faceCapacity(mesh,face,s.radial),al=referenceFace(mesh,face,middle),v=u.normalVelocity()[f];
        const long double expected=(l*al+g*(area-al))*v;
        if(exact)peakFlux=std::max(peakFlux,normalized(material.faceRates()[f],expected,(l+g)*area*std::max(1.L,std::abs(v)),1e-10,"analytic mass face rate"));
        normalized(material.faceIntegrals()[f],dt*material.faceRates()[f],std::max(1.L,std::abs(dt*expected)),1e-14,"midpoint time integral");
        error.flux+=std::abs(material.faceRates()[f]-expected);totalArea+=area;sumAbs+=std::abs(dt*expected);
        for(auto side:{std::pair<int,int>{face.negativeCell,1},{face.positiveCell,-1}})if(side.first>=0){
            outward[layout.pressureDof(side.first)]+=side.second*dt*expected;
            const auto b=material.beginGeometry().publication().state().geometry().snapshot().nativeCellBounds(side.first);
            const double width=face.horizontal?b[1].x-b[0].x:b[1].y-b[0].y;
            if(face.tangentSpacing<width)coarse[{side.first,(face.horizontal?2:0)+(side.second>0?0:1)}].push_back(f);
        }
        boundary+=face.outwardSign*dt*expected;
        if(face.negativeCell>=0 && face.positiveCell>=0){
            check(face.outwardSign==0,"internal face marked boundary");
            ++internalChecks; // Incidence coverage; independent cell sums checked below.
        }
        ++faceChecks;
    }
    for(const auto& entry:coarse){
        long double actual=0;for(int f:entry.second)actual+=material.faceIntegrals()[f];
        const auto b=material.beginGeometry().publication().state().geometry().snapshot().nativeCellBounds(entry.first.first);
        const int side=entry.first.second;const bool horizontal=side>=2;
        const double normal=horizontal?(side==2?b[1].y:b[0].y):(side==0?b[1].x:b[0].x);
        const double lo=horizontal?b[0].x:b[0].y,hi=horizontal?b[1].x:b[1].y;
        const auto al=liquidInterval(lo,hi,normal,horizontal,middle);
        const auto area=(hi-lo)*(s.radial?CartesianMeasure::twoPi*(horizontal?(lo+hi)/2:normal):1.L);
        const auto expected=dt*(l*al+g*(area-al))*u.normalVelocity()[entry.second.front()];
        if(exact)peakCoarse=std::max(peakCoarse,normalized(actual,expected,std::max(1.L,std::abs(expected)),1e-12,"coarse side shared fine material flux"));++coarseChecks;
    }
    long double signedDefect=0,delta=0;
    for(int p=0;p<layout.pressureCount();++p){
        const auto bounds=material.beginGeometry().publication().state().geometry().snapshot().nativeCellBounds(layout.nativeCell(p));
        const auto cap=capacity(bounds[0],bounds[1],s.radial);const auto& cell=material.cells()[p];
        const auto va=liquidCell(bounds[0],bounds[1],shapeAt(s,material.beginGeometry().timeLayer().time()));
        const auto vm=liquidCell(bounds[0],bounds[1],middle);
        const auto vb=liquidCell(bounds[0],bounds[1],shapeAt(s,material.endGeometry().timeLayer().time()));
        const std::array<long double,3> expected{{l*va+g*(cap-va),l*vm+g*(cap-vm),l*vb+g*(cap-vb)}};
        const std::array<long double,3> actual{{cell.massBegin,cell.massMidpoint,cell.massEnd}};
        for(int i=0;i<3;++i){if(exact)peakMass=std::max(peakMass,normalized(actual[i],expected[i],(l+g)*cap,1e-10,"analytic geometry mass"));error.mass+=std::abs(actual[i]-expected[i]);}
        const auto change=(l-g)*(vb-va),defect=change+outward[p];
        if(exact){
            peakDefect=std::max(peakDefect,normalized(cell.massChange,change,(l+g)*cap,1e-10,"geometry mass change"));
            peakDefect=std::max(peakDefect,normalized(cell.outwardIntegral,outward[p],(l+g)*cap+sumAbs,1e-10,"independent native cell outward integral"));
            peakDefect=std::max(peakDefect,normalized(cell.integratedDefect,defect,(l+g)*cap+sumAbs,1e-10,"analytic nonconservative defect"));
        }
        normalized(cell.defectRate*dt,cell.integratedDefect,std::max(1.L,std::abs(cell.integratedDefect)),1e-14,"defect rate units");
        long double gross=std::abs(cell.massChange);for(int f=0;f<layout.velocityCount();++f){const auto face=layout.domain().face(mesh,layout.nativeFace(f));if(face.negativeCell==layout.nativeCell(p) || face.positiveCell==layout.nativeCell(p))gross+=std::abs(material.faceIntegrals()[f]);}
        normalized(cell.relativeDefect,gross>0?std::abs(cell.integratedDefect)/gross:0,1,1e-14,"relative defect scale");
        signedDefect+=cell.integratedDefect;delta+=cell.massChange;totalVolume+=cap;++cellChecks;
    }
    if(exact)peakBalance=std::max(peakBalance,normalized(material.report().boundaryIntegral,boundary,std::max(1.L,sumAbs),1e-12,"independent selected boundary integral"));
    peakBalance=std::max(peakBalance,normalized(signedDefect,delta+material.report().boundaryIntegral,std::max(1.L,sumAbs),1e-12,"global material defect identity"));
    check(material.report().balanceIdentity<=1e-12,"prepared global material balance defect");
    check(material.beginGeometry().publication().state().primary().phaseVolumes.empty() && material.endGeometry().publication().state().primary().phaseVolumes.empty(),"restored conservative inventory");
    material.validateFor(mesh,material.sources(),material.timeLayer());++scenarios;
    error.mass/=3*(l+g)*totalVolume;error.flux/=(l+g)*totalArea;return error;
}
CchMaterialTransport staticMaterial(const QuadtreeMesh& mesh,const CchPhysicalDivergence& metrics,const Case& s,const CchNativeValues& u,double l,double g)
{
    const auto a=CchPhaseGeometry::build(mesh,metrics,publishAt(mesh,s,{}));
    const auto m=CchPhaseGeometry::build(mesh,metrics,publishAt(mesh,s,{1,0,.04,Position::Midpoint}));
    const auto b=CchPhaseGeometry::build(mesh,metrics,publishAt(mesh,s,{1,0,.04,Position::End}));
    return CchMaterialTransport::prepare(mesh,a,rate(mesh,m,u),b,l,g);
}
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
void singleFace(const QuadtreeMesh& mesh,const Case& s,const CchPhaseGeometry& a,const CchPhaseGeometry& m,const CchPhaseGeometry& b,bool subset)
{
    const auto& layout=m.layout();int chosen=-1;
    for(int f=0;f<layout.velocityCount();++f){const auto face=layout.domain().face(mesh,layout.nativeFace(f));const auto phase=m.faceMeasures()[f];
        if(face.negativeCell>=0 && face.positiveCell>=0 && phase.liquid>0 && phase.gas>0){chosen=f;break;}}
    check(chosen>=0,"missing shared split phase basis face");
    const auto face=layout.domain().face(mesh,layout.nativeFace(chosen));
    std::vector<double> normal(layout.velocityCount());normal[chosen]=1;
    const auto material=CchMaterialTransport::prepare(mesh,a,rate(mesh,m,CchNativeValues(layout,normal,std::vector<double>(layout.pressureCount()))),b,1e6,1);
    const auto& actual=material.midpointRate().spatialAdvector();const auto& reconstruction=actual.reconstruction();std::vector<long double> atoms,outward(layout.pressureCount());long double boundary=0,gross=0;
    for(std::size_t i=0;i<reconstruction.atoms().size();++i){const auto phase=referenceAtom(reconstruction.atoms()[i],actual.atomPolynomials()[i].beginCoefficients,s);atoms.push_back(1e6*phase.liquid+phase.gas);}
    for(int f=0;f<layout.velocityCount();++f){long double expected=0;for(int i:reconstruction.primalAtoms()[f])expected+=atoms[i];expected*=material.deltaTime();const auto native=layout.domain().face(mesh,layout.nativeFace(f));
        normalized(material.faceIntegrals()[f],expected,1e6*faceCapacity(mesh,native,s.radial),1e-12,"basis independent varying material trace");gross+=std::abs(expected);boundary+=native.outwardSign*expected;
        if(native.negativeCell>=0)outward[layout.pressureDof(native.negativeCell)]+=expected;
        if(native.positiveCell>=0)outward[layout.pressureDof(native.positiveCell)]-=expected;}
    for(int p=0;p<layout.pressureCount();++p){
        normalized(material.cells()[p].outwardIntegral,outward[p],gross,1e-12,"basis actual neighbor contribution");
        check(material.cells()[p].massChange==0,"basis fixture geometry changed");
    }
    normalized(material.report().integratedDefect,boundary,gross,1e-12,"basis true shared cancellation and boundary");
    std::cout<<"BASIS,"<<s.radial<<','<<s.amr<<','<<subset<<','<<chosen<<'\n';
}
void crossing(bool rz,int amr)
{
    auto mesh=makeMesh(4,amr);Case s;s.pureFirstJet=true;s.planar=true;s.radial=rz;s.motion=1;s.speed=.4;s.planeOffset=.49;
    const auto metrics=metricsFor(*mesh,rz);const auto old=initial(*mesh,s);std::string stage;
    const auto mp=advanceFirstJet(*mesh,old,s,{1,0,.1,Position::Midpoint},stage);
    const auto ep=advanceFirstJet(*mesh,old,s,{1,0,.1,Position::End},stage);
    const auto a=CchPhaseGeometry::build(*mesh,metrics,old),m=CchPhaseGeometry::build(*mesh,metrics,mp),b=CchPhaseGeometry::build(*mesh,metrics,ep);
    const auto u=velocity(*mesh,metrics.layout(),rz?Vector2D{0,s.speed}:Vector2D{s.speed,0});
    const auto material=CchMaterialTransport::prepare(*mesh,a,rate(*mesh,m,u),b,1e6,1);
    inspect(*mesh,s,material,u);check(material.report().l1IntegratedDefect>100,"crossing defect silently removed");
    normalized(material.report().integratedDefect,0,1e6,1e-10,"crossing global mass closure");
    std::cout<<"CROSSING,"<<rz<<','<<amr<<','<<material.report().l1IntegratedDefect<<','<<material.report().maxAbsDefectRate<<'\n';
}
template<class F>void reject(F f,const std::string& why){rejects(f,why);++rejected;}
void guards()
{
    auto mesh=makeMesh(4);Case s;s.pureFirstJet=true;s.planar=true;s.planeOffset=.43;
    const auto metrics=metricsFor(*mesh,false);const auto u=velocity(*mesh,metrics.layout(),{.3,-.2});
    const auto a=CchPhaseGeometry::build(*mesh,metrics,publishAt(*mesh,s,{}));
    const auto m=CchPhaseGeometry::build(*mesh,metrics,publishAt(*mesh,s,{1,0,.04,Position::Midpoint}));
    const auto b=CchPhaseGeometry::build(*mesh,metrics,publishAt(*mesh,s,{1,0,.04,Position::End}));
    const auto q=rate(*mesh,m,u);const auto material=CchMaterialTransport::prepare(*mesh,a,q,b,1000,1);const auto saved=material.faceIntegrals();
    const auto& geometry=material.momentumGeometry();const auto& flux=material.midpointFlux();const auto copiedFlux=flux;
    const auto reused=CchMaterialTransport::prepare(*mesh,a,copiedFlux,b);
    check(&flux.atomRates()==&copiedFlux.atomRates() && &flux.faceRates()==&reused.faceRates() && &flux.dualFaceRates()==&reused.dualFaceRates(),"material interval recomputed unique flux");
    material.validateFor(*mesh,geometry,material.sources(),material.timeLayer());
    const auto freshGeometry=CchMomentumGeometry::build(*mesh,metrics);
    reject([&]{material.validateFor(*mesh,freshGeometry,material.sources(),material.timeLayer());},"different momentum geometry instance");
    const auto freshAdvector=CchAdvectorReconstruction::prepare(*mesh,geometry).evaluate(*mesh,u,m.geometryIdentity(),m.timeLayer());
    for(auto foreign:{CchAdvectorIdentity{},freshAdvector.identity()}){auto source=material.sources();source.advector=foreign;
        reject([&]{material.validateFor(*mesh,source,material.timeLayer());},"different actual advector identity");
        reject([&]{flux.validateFor(*mesh,geometry,m.geometryIdentity(),foreign,m.timeLayer());},"different actual advector identity");}
    auto wrongFluxTime=m.timeLayer();++wrongFluxTime.step;
    reject([&]{flux.validateFor(*mesh,geometry,m.geometryIdentity(),flux.advectorIdentity(),wrongFluxTime);},"different full time layer");
    reject([&]{flux.validateFor(*mesh,geometry,a.geometryIdentity(),flux.advectorIdentity(),m.timeLayer());},"different geometry identity");
    for(const auto* p:std::array<const CchPhaseGeometry*,3>{{&a,&m,&b}}){const auto instant=CchMaterialFlux::prepare(*mesh,rate(*mesh,*p,u),1000,1);
        check(instant.timeLayer().position==p->timeLayer().position && instant.geometryIdentity()==p->geometryIdentity(),"instantaneous flux retagged actual layer");}
    for(double bad:{0.,-1.,std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::infinity()}){
        reject([&]{CchMaterialTransport::prepare(*mesh,a,q,b,bad,1);},"nonpositive/nonfinite density");
        reject([&]{CchMaterialTransport::prepare(*mesh,a,q,b,1,bad);},"nonpositive/nonfinite density");
    }
    reject([&]{CchMaterialTransport::prepare(*mesh,a,rate(*mesh,b,u),b,1,1);},"requires Midpoint/End");
    reject([&]{CchMaterialTransport::prepare(*mesh,a,q,m,1,1);},"requires Midpoint/End");
    const auto wrong=CchPhaseGeometry::build(*mesh,metrics,publishAt(*mesh,s,{1,-.1,.14,Position::Midpoint}));
    reject([&]{CchMaterialTransport::prepare(*mesh,a,rate(*mesh,wrong,u),b,1,1);},"different full candidate interval");
    for(auto clock:{HermiteGalsTimeLayer{3,0,0,Position::Initial},HermiteGalsTimeLayer{0,.01,.01,Position::Initial},HermiteGalsTimeLayer{1,-.04,.04,Position::Midpoint}}){
        const auto wrongBegin=CchPhaseGeometry::build(*mesh,metrics,publishAt(*mesh,s,clock));
        reject([&]{CchMaterialTransport::prepare(*mesh,wrongBegin,q,b,1,1);},"previous Initial/End");
    }
    const auto fresh=CchPhaseGeometry::build(*mesh,metricsFor(*mesh,false),a.publication());
    reject([&]{CchMaterialTransport::prepare(*mesh,fresh,q,b,1,1);},"different native layout");
    const auto otherMetrics=CchPhysicalDivergence::build(*mesh,CchIntegralDivergence::build(*mesh,metrics.layout()));
    const auto otherDi=CchPhaseGeometry::build(*mesh,otherMetrics,a.publication());
    reject([&]{CchMaterialTransport::prepare(*mesh,otherDi,q,b,1,1);},"different retained D_I");
    for(int which=0;which<3;++which){auto id=material.sources();if(which==0)id.begin={};if(which==1)id.midpoint={};if(which==2)id.end={};reject([&]{material.validateFor(*mesh,id,material.timeLayer());},"different geometry identity");}
    auto time=material.timeLayer();time.begin=-.01;
    reject([&]{material.validateFor(*mesh,material.sources(),time);},"different full time layer");
    time=material.timeLayer();time.step=2;reject([&]{material.validateFor(*mesh,material.sources(),time);},"different full time layer");
    time=material.timeLayer();time.position=Position::End;reject([&]{material.validateFor(*mesh,material.sources(),time);},"different full time layer");
    reject([&]{CchMaterialTransport::prepare(*mesh,a,rate(*mesh,m,velocity(*mesh,metrics.layout(),{100,100})),b,std::numeric_limits<double>::max(),std::numeric_limits<double>::max());},"nonfinite atom material");
    reject([&]{CchMaterialTransport::prepare(*mesh,a,q,b,std::numeric_limits<double>::denorm_min(),1);},"unrepresentable atom material rate");
    const double tinyTime=4*std::numeric_limits<double>::denorm_min();
    const auto tinyMid=CchPhaseGeometry::build(*mesh,metrics,publishAt(*mesh,s,{1,0,tinyTime,Position::Midpoint}));
    const auto tinyEnd=CchPhaseGeometry::build(*mesh,metrics,publishAt(*mesh,s,{1,0,tinyTime,Position::End}));
    reject([&]{CchMaterialTransport::prepare(*mesh,a,rate(*mesh,tinyMid,u),tinyEnd,1,1);},"unrepresentable face material integral");
    const double adjacent=std::nextafter(1.,2.);
    const auto laterStart=CchPhaseGeometry::build(*mesh,metrics,publishAt(*mesh,s,{0,1,1,Position::Initial}));
    const auto collapsedMid=CchPhaseGeometry::build(*mesh,metrics,publishAt(*mesh,s,{1,1,adjacent,Position::Midpoint}));
    const auto laterEnd=CchPhaseGeometry::build(*mesh,metrics,publishAt(*mesh,s,{1,1,adjacent,Position::End}));
    reject([&]{CchMaterialTransport::prepare(*mesh,laterStart,rate(*mesh,collapsedMid,u),laterEnd,1,1);},"unrepresentable midpoint interval");
    check(material.faceIntegrals()==saved,"rejection polluted accepted material data");
    auto& vertex=mesh->verts(mesh->vertIndexToLabel[0]);const auto old=vertex.coord();vertex.setCoord(old+Vector2D{.01,0});
    reject([&]{material.validateFor(*mesh,material.sources(),material.timeLayer());},"NativeMeshIdentity");vertex.setCoord(old);
    DictionaryMesh::RefineBox box;box.level=1;box.min={.25,.25};box.max={.75,.75};mesh->applyRefineBoxes({box});
    reject([&]{material.validateFor(*mesh,material.sources(),material.timeLayer());},"NativeMeshIdentity");
    const auto rebuiltMetrics=metricsFor(*mesh,false);const auto copy=staticMaterial(*mesh,rebuiltMetrics,s,velocity(*mesh,rebuiltMetrics.layout(),{.3,-.2}),1000,1);
    inspect(*mesh,s,copy,velocity(*mesh,copy.layout(),{.3,-.2}));
    const auto retained=copy;check(&retained.faceRates()==&copy.faceRates(),"material copy rebuilt data");
    auto equivalent=makeMesh(4,1);retained.validateFor(*equivalent,retained.sources(),retained.timeLayer());
}
void spatial(bool rz,int amr)
{
    std::vector<Error> errors;
    for(int n:{4,8,16}){
        auto mesh=makeMesh(n,amr);Case s;s.pureFirstJet=true;s.radial=rz;s.amr=amr!=0;s.nonpolynomial=true;
        const auto metrics=metricsFor(*mesh,rz);const auto u=velocity(*mesh,metrics.layout(),rz?Vector2D{0,.2}:Vector2D{.3,-.2});
        const auto material=staticMaterial(*mesh,metrics,s,u,1000,1);const auto e=inspect(*mesh,s,material,u,false);errors.push_back(e);++spatialRows;
        std::cout<<"SPACE,"<<rz<<','<<amr<<','<<n<<','<<e.mass<<','<<e.flux<<'\n';
    }
    for(int i=1;i<3;++i)check(errors[i].mass<errors[i-1].mass && errors[i].flux<errors[i-1].flux,"material spatial errors did not decrease");
    const auto mo=std::log2(errors[1].mass/errors[2].mass),fo=std::log2(errors[1].flux/errors[2].flux);
    check(mo>=2 && fo>=2,"material spatial order below gate");std::cout<<"SPACE_ORDER,"<<rz<<','<<amr<<','<<mo<<','<<fo<<'\n';
}
void temporal(bool rz,int amr)
{
    auto mesh=makeMesh(4,amr);const auto metrics=metricsFor(*mesh,rz);Case s;s.pureFirstJet=true;s.radial=rz;s.amr=amr!=0;s.motion=1;s.speed=.2;
    s.cx=.43;s.cy=rz?.43:.375;s.radius=rz?.19:.08;
    // Aggregate real subfaces on a fixed physical side. The entire moving
    // disk section fits inside this side, away from temporal tangencies.
    std::vector<int> selected;long double area=0;
    for(int f=0;f<metrics.layout().velocityCount();++f){const auto face=metrics.layout().domain().face(*mesh,metrics.layout().nativeFace(f));const auto b=faceBounds(*mesh,face);
        if(rz ? face.horizontal && face.center.y==.5 && b[0].x>=0 && b[1].x<=.25 : !face.horizontal && face.center.x==.5 && b[0].y>=.25 && b[1].y<=.5){selected.push_back(f);area+=faceCapacity(*mesh,face,rz);}}
    check(!selected.empty(),"missing temporal physical side");
    constexpr double total=.1,l=1000,g=1;const long double d0=.5-(rz?s.cy:s.cx),d1=d0-s.speed*total,r=s.radius;
    const auto arc=[&](long double d){return (d*std::sqrt(r*r-d*d)+r*r*std::asin(d/r))/2;};
    const long double liquidIntegral=rz?pi*(r*r*(d0-d1)-(d0*d0*d0-d1*d1*d1)/3):2*(arc(d0)-arc(d1));
    const long double exact=l*liquidIntegral+g*(s.speed*total*area-liquidIntegral);double previous=0,order=0;
    for(int count:{4,8,16}){
        long double actual=0;auto previousGeometry=CchPhaseGeometry::build(*mesh,metrics,publishAt(*mesh,s,{}));
        for(int k=0;k<count;++k){
            const double begin=total*k/count,end=total*(k+1)/count;
            const auto mid=CchPhaseGeometry::build(*mesh,metrics,publishAt(*mesh,s,{k+1,begin,end,Position::Midpoint}));
            const auto arrival=CchPhaseGeometry::build(*mesh,metrics,publishAt(*mesh,s,{k+1,begin,end,Position::End}));
            const auto u=velocity(*mesh,metrics.layout(),rz?Vector2D{0,s.speed}:Vector2D{s.speed,0});
            const auto material=CchMaterialTransport::prepare(*mesh,previousGeometry,rate(*mesh,mid,u),arrival,l,g);
            for(int f:selected)actual+=material.faceIntegrals()[f];previousGeometry=arrival;
        }
        const double e=static_cast<double>(std::abs(actual-exact)/std::abs(exact));
        if(count>4){check(e<previous,"material temporal error did not decrease");order=std::log2(previous/e);}previous=e;++temporalRows;
        std::cout<<"TIME,"<<rz<<','<<amr<<','<<count<<','<<actual<<','<<exact<<','<<e<<'\n';
    }
    check(order>=1.8,"material temporal order below gate");std::cout<<"TIME_ORDER,"<<rz<<','<<amr<<','<<order<<'\n';
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
void spatialMaterial()
{
    int scenes=0,counterflow=0,tiny=0;std::size_t rows=0;long double peak=0;
    for(bool rz:{false,true})for(int amr:{0,1}){auto mesh=makeMesh(4,amr);const auto metrics=metricsFor(*mesh,rz);const auto geometry=CchMomentumGeometry::build(*mesh,metrics);const auto reconstruction=CchAdvectorReconstruction::prepare(*mesh,geometry);
        for(int kind=0;kind<10;++kind){Case s;s.pureFirstJet=true;s.radial=rz;s.planar=kind!=0;s.radius=rz?.34:.19;
            s.planeOffset=kind==1?.43:kind==2?2:kind==3?-1:kind==4?std::ldexp(1.,-40):kind==5?1-std::ldexp(1.,-40):kind<8?(amr?.125:.375):.5;s.planeSign=kind==7 || kind==9?-1:1;
            const auto a=CchPhaseGeometry::build(*mesh,metrics,publishAt(*mesh,s,{})),m=CchPhaseGeometry::build(*mesh,metrics,publishAt(*mesh,s,{1,0,.04,Position::Midpoint})),b=CchPhaseGeometry::build(*mesh,metrics,publishAt(*mesh,s,{1,0,.04,Position::End}));
            for(int mode=0;mode<4;++mode){if(kind<6 && mode!=2)continue;const auto advector=reconstruction.evaluate(*mesh,mmsValues(*mesh,metrics.layout(),mode),m.geometryIdentity(),m.timeLayer());const auto phase=CchPhaseRate::prepareSpatial(*mesh,m,advector);
                for(double rho:{1.,1000.,1e6}){const double l=kind==5?1:rho,g=kind==5?rho:1;const auto flux=CchMaterialFlux::prepare(*mesh,phase,l,g);const auto material=CchMaterialTransport::prepare(*mesh,a,flux,b);std::vector<long double> expected;
                    check(&flux.atomRates()==&material.atomRates() && &flux.faceRates()==&material.faceRates() && &flux.dualFaceRates()==&material.dualFaceRates(),"interval duplicated shared F storage");
                    material.validateFor(*mesh,geometry,material.sources(),material.timeLayer());
                    for(std::size_t i=0;i<reconstruction.atoms().size();++i){const auto& atom=reconstruction.atoms()[i];const auto ref=referenceAtom(atom,mmsCoefficients(atom,rz,mode),s);const auto q=l*ref.liquid+g*ref.gas;expected.push_back(q);
                        peak=std::max(peak,normalized(flux.atomRates()[i],q,(l+g)*atom.area,1e-12,"independent spatial atom F"));++rows;
                        if(rho>1 && (kind==4 || kind==5)){const auto small=kind==4?ref.liquid:ref.gas;if(std::abs(small)>0 && std::abs(small)<1e-10L*atom.area){const auto actual=kind==4?phase.atomRates()[i].liquid:phase.atomRates()[i].gas;
                            normalized(actual,small,std::abs(small),1e-12,"small phase actual rate scale");const auto majorDensity=kind==4?g:l;
                            check(flux.atomRates()[i]!=majorDensity*advector.atomPolynomials()[i].integral(false,0,advector.atomPolynomials()[i].length),"small weighted material contribution vanished");++tiny;}}
                    }
                    for(int u=0;u<metrics.layout().velocityCount();++u){long double ref=0;for(int i:reconstruction.primalAtoms()[u])ref+=expected[i];const auto area=metrics.faceMeasures()[u];
                        peak=std::max(peak,normalized(material.faceRates()[u],ref,(l+g)*area,1e-12,"independent spatial primal F"));
                        if(kind>=6 && mode==3 && phase.totalFaceRates()[u]==0 && phase.faceRates()[u].liquid!=0 && phase.faceRates()[u].gas!=0 && rho>1){check(material.faceRates()[u]!=0,"zero-net counterflow F vanished");++counterflow;}}
                    for(std::size_t i=0;i<geometry.segments().size();++i){long double ref=0;for(int j:reconstruction.dualAtoms()[i])ref+=expected[j];
                        peak=std::max(peak,normalized(material.dualFaceRates()[i],ref,(l+g)*geometry.segments()[i].area,1e-12,"independent spatial dual F"));}
                    const auto copied=material;check(&copied.atomRates()==&flux.atomRates(),"copied material recreated atom owner");check(material.report().balanceIdentity<=1e-12,"spatial material pressure balance");++scenes;
                }
            }
        }
    }
    check(scenes==264 && rows>0 && counterflow>0 && tiny>0,"spatial material MMS/counterflow/tiny coverage");
    std::cout<<"SPATIAL_MATERIAL_SUMMARY,"<<scenes<<','<<rows<<','<<counterflow<<','<<tiny<<','<<peak<<",PASS\n";
}
}
namespace {
void impermeableMaterial()
{
    int scenes=0,boundaryAtoms=0,boundarySegments=0,rawWitnesses=0,internalRows=0,selectiveCounterflow=0,phaseCounterflow=0;long double peak=0,rawPeak=0;
    for(bool rz:{false,true})for(int amr:{0,1})for(bool subset:{false,true}){auto mesh=makeMesh(4,amr);const auto metrics=metricsFor(*mesh,rz,subset);const auto g=CchMomentumGeometry::build(*mesh,metrics);const auto& l=g.layout();const auto r=CchAdvectorReconstruction::prepare(*mesh,g);
        std::vector<int> walls;std::vector<double> v(l.velocityCount());const double right=subset?.5:1;
        const auto psi=[&](Vector2D p){const double radial=rz?p.x*p.x*(right-p.x)*(right-p.x):p.x*(right-p.x);return radial*p.y*(1-p.y);};
        for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(*mesh,l.nativeFace(u));const auto ends=faceBounds(*mesh,f);const auto area=metrics.faceMeasures()[u];v[u]=area>0?(f.horizontal?-1:1)*(rz?static_cast<double>(2*pi):1)*(psi(ends[1])-psi(ends[0]))/area:0;if(f.outwardSign!=0){check(v[u]==0,"closed input native mean nonzero");walls.push_back(u);}}
        const auto law=CchAdvectorBoundary::impermeable(*mesh,g,walls);const CchNativeValues values(l,v,std::vector<double>(l.pressureCount()));
        for(bool cut:{false,true}){Case s;s.pureFirstJet=true;s.radial=rz;s.planar=true;s.planeOffset=cut?.43:2;
            const auto a=CchPhaseGeometry::build(*mesh,metrics,publishAt(*mesh,s,{})),m=CchPhaseGeometry::build(*mesh,metrics,publishAt(*mesh,s,{1,0,.04,Position::Midpoint})),b=CchPhaseGeometry::build(*mesh,metrics,publishAt(*mesh,s,{1,0,.04,Position::End}));
            const auto raw=r.evaluate(*mesh,values,m.geometryIdentity(),m.timeLayer()),actual=r.evaluate(*mesh,values,m.geometryIdentity(),m.timeLayer(),law);check(actual.boundary()->sameIdentity(law) && actual.identity()!=raw.identity(),"bound actual provenance");
            const auto phase=CchPhaseRate::prepareSpatial(*mesh,m,actual),rawPhase=CchPhaseRate::prepareSpatial(*mesh,m,raw);
            const auto original=metrics.apply(*mesh,values);for(int p=0;p<l.pressureCount();++p){long double scale=0;const int native=l.nativeCell(p);for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(*mesh,l.nativeFace(u));if(f.negativeCell==native || f.positiveCell==native)scale+=std::abs(static_cast<long double>(metrics.faceMeasures()[u])*v[u]);}
                peak=std::max(peak,normalized(actual.cellDivergence()[p]*metrics.cellVolumes()[p],original[p]*metrics.cellVolumes()[p],scale,1e-12,"closed actual B retains original DI"));}
            for(double rho:{1.,1000.,1e6}){const auto flux=CchMaterialFlux::prepare(*mesh,phase,rho,1),rawFlux=CchMaterialFlux::prepare(*mesh,rawPhase,rho,1);const auto transport=CchMaterialTransport::prepare(*mesh,a,flux,b);std::vector<long double> expected;
                check(&transport.atomRates()==&flux.atomRates() && transport.sources().advector==actual.identity(),"closed material unique source");
                for(std::size_t i=0;i<r.atoms().size();++i){const auto& atom=r.atoms()[i];const auto ref=referenceAtom(atom,actual.atomPolynomials()[i].beginCoefficients,s);const auto q=rho*ref.liquid+ref.gas;expected.push_back(q);const auto& c=actual.atomPolynomials()[i].beginCoefficients;const auto length=actual.atomPolynomials()[i].length;const auto work=(rho+1)*atom.area*referenceWork(c,length);
                    peak=std::max(peak,normalized(flux.atomRates()[i],q,work,1e-12,"independent bound phase/material integral"));
                    if(atom.primalVelocity>=0 && std::binary_search(walls.begin(),walls.end(),atom.primalVelocity)){for(auto k:c)check(k==0,"bound normal polynomial not exactly zero");for(auto k:actual.atomPolynomials()[i].endCoefficients)check(k==0,"reverse bound normal polynomial not exactly zero");check(flux.atomRates()[i]==0 && phase.atomRates()[i].liquid==0 && phase.atomRates()[i].gas==0,"impermeable atom leaked phase/material");++boundaryAtoms;
                        if(rawFlux.atomRates()[i]!=0){++rawWitnesses;rawPeak=std::max(rawPeak,std::abs(rawFlux.atomRates()[i]));}}
                }
                std::vector<long double> outward(l.pressureCount()),gross(l.pressureCount());for(int u=0;u<l.velocityCount();++u){long double sum=0,scale=0;for(int i:r.primalAtoms()[u]){sum+=expected[i];scale+=std::abs(expected[i]);}peak=std::max(peak,normalized(flux.faceRates()[u],sum,scale,1e-12,"bound native subface sum"));const auto f=l.domain().face(*mesh,l.nativeFace(u));
                    for(auto side:{std::pair<int,int>{f.negativeCell,1},{f.positiveCell,-1}})if(side.first>=0){const int p=l.pressureDof(side.first);outward[p]+=side.second*sum;gross[p]+=scale;}}
                for(std::size_t j=0;j<g.segments().size();++j){const auto& edge=g.segments()[j];long double sum=0,scale=0;for(int i:r.dualAtoms()[j]){sum+=expected[i];scale+=std::abs(expected[i]);}peak=std::max(peak,normalized(flux.dualFaceRates()[j],sum,scale,1e-12,"bound actual dual subface sum"));
                    if(edge.negativeVelocity<0 || edge.positiveVelocity<0){check(edge.primalVelocity>=0 && std::binary_search(walls.begin(),walls.end(),edge.primalVelocity),"dual boundary lacks actual prescribed native origin");check(flux.dualFaceRates()[j]==0,"bound dual boundary flux nonzero");++boundarySegments;}else ++internalRows;}
                for(int p=0;p<l.pressureCount();++p){check(transport.cells()[p].massChange==0,"stationary bound geometry changed");peak=std::max(peak,normalized(transport.cells()[p].integratedDefect,.04L*outward[p],.04L*gross[p],1e-12,"true bound R from independent incidence"));}
                check(transport.report().boundaryIntegral==0,"closed component material boundary nonzero");++scenes;
            }
        }
        if(!subset){Case s;s.pureFirstJet=true;s.radial=rz;s.planar=true;s.planeOffset=amr?.125:.375;auto input=mmsValues(*mesh,l,3);int chosen=-1;for(int u:walls)if(metrics.faceMeasures()[u]>0 && input.normalVelocity()[u]==0){chosen=u;break;}check(chosen>=0,"selective zero-mean boundary witness absent");
            // Constrain just one explicit boundary face. Other zero-mean faces
            // must retain their nonzero tangent trace. An interior coincident
            // plane also checks actual liquid/gas counterflow independently.
            const auto m=CchPhaseGeometry::build(*mesh,metrics,publishAt(*mesh,s,{}));
            const auto one=CchAdvectorBoundary::impermeable(*mesh,g,{chosen});const auto actual=r.evaluate(*mesh,input,m.geometryIdentity(),m.timeLayer(),one);
            for(int u=0;u<l.velocityCount();++u)if(u!=chosen && input.normalVelocity()[u]==0){bool varied=false;for(int i:r.primalAtoms()[u])for(auto c:actual.atomPolynomials()[i].beginCoefficients)varied=varied||c!=0;if(varied)++selectiveCounterflow;}
            for(int i:r.primalAtoms()[chosen])for(auto c:actual.atomPolynomials()[i].beginCoefficients)check(c==0,"selective wall not enforced");
            const auto phase=CchPhaseRate::prepareSpatial(*mesh,m,actual);const auto flux=CchMaterialFlux::prepare(*mesh,phase,1e6,1);
            for(std::size_t i=0;i<r.atoms().size();++i){const auto& atom=r.atoms()[i];const auto& polynomial=actual.atomPolynomials()[i];const auto& c=polynomial.beginCoefficients;const auto length=polynomial.length;const auto ref=referenceAtom(atom,c,s);
                // The signed phase integral can vanish by cancellation within
                // one phase. Bound GL3 rounding by actual polynomial work,
                // without a floor or a change to the tiny-phase own-scale gate.
                const auto scale=(1e6L+1)*atom.area*referenceWork(c,length);
                if(std::abs(flux.atomRates()[i]-(1e6*ref.liquid+ref.gas))>1e-12L*scale){const auto& a=r.atoms()[i];const auto& p=actual.atomPolynomials()[i];const auto q=phase.atomRates()[i];std::cout<<"SELECTIVE_WITNESS,"<<rz<<','<<amr<<','<<i<<','<<a.normalAxis<<','<<a.primalVelocity<<','<<a.ownerPressure<<','<<a.begin.x<<','<<a.begin.y<<','<<a.end.x<<','<<a.end.y;for(auto c:p.beginCoefficients)std::cout<<','<<c;for(auto c:p.endCoefficients)std::cout<<','<<c;std::cout<<','<<ref.liquid<<','<<ref.gas<<','<<q.liquid<<','<<q.gas<<','<<flux.atomRates()[i]<<std::endl;}
                peak=std::max(peak,normalized(flux.atomRates()[i],1e6*ref.liquid+ref.gas,scale,1e-12,"selective-law independent causal material counterflow"));}
            int local=0;for(std::size_t i=0;i<g.segments().size();++i){const auto q=phase.dualFaceRates()[i];if(q.liquid*q.gas<0 && std::abs(q.liquid+q.gas)<=1e-12L*(std::abs(q.liquid)+std::abs(q.gas))){check(flux.dualFaceRates()[i]!=0,"selective law erased phase counterflow");++local;}}
            // A true AMR common side can split the positive and negative
            // parts across separate shared subfaces. Check their real sum.
            for(std::size_t u=0;u<g.volumes().size();++u)for(int side=0;side<4;++side){const auto& ids=g.volumes()[u].sides[side];if(ids.size()<2)continue;CchPhaseMeasure q;long double material=0;for(int i:ids){q.liquid+=phase.dualFaceRates()[i].liquid;q.gas+=phase.dualFaceRates()[i].gas;material+=flux.dualFaceRates()[i];}
                if(q.liquid*q.gas<0 && std::abs(q.liquid+q.gas)<=1e-12L*(std::abs(q.liquid)+std::abs(q.gas))){check(material!=0,"selective law erased shared subface counterflow");normalized(material,1e6*q.liquid+q.gas,1e6*std::abs(q.liquid)+std::abs(q.gas),1e-12,"single F across actual split side");std::cout<<"BOUND_COUNTERFLOW_SIDE,"<<rz<<','<<amr<<','<<u<<','<<side<<','<<ids.size()<<','<<q.liquid<<','<<q.gas<<','<<material<<'\n';++local;}}
            check(local>0,"selective-law phase counterflow missing");phaseCounterflow+=local;
        }
    }
    check(scenes==48 && rawWitnesses>0 && selectiveCounterflow>0 && internalRows>0 && phaseCounterflow>0,"impermeable/counterflow coverage missing");std::cout<<"IMPERMEABLE_MATERIAL,"<<scenes<<','<<boundaryAtoms<<','<<boundarySegments<<','<<rawWitnesses<<','<<internalRows<<','<<selectiveCounterflow<<','<<peak<<','<<rawPeak<<",PASS\n";std::cout<<"IMPERMEABLE_PHASE_COUNTERFLOW,"<<phaseCounterflow<<",PASS\n";
}
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
void spaceTimeIdentities(const QuadtreeMesh& mesh,const CchMaterialTransport& m,const CchDualMaterial& dual)
{
    const auto& g=m.momentumGeometry();const auto& l=m.layout();const auto& r=m.midpointRate().spatialAdvector().reconstruction();
    long double peak=0;std::vector<long double> outward(l.pressureCount()),gross(l.pressureCount());
    for(std::size_t i=0;i<r.atoms().size();++i){const auto q=m.phaseAtomIntegrals()[i];const auto f=m.liquidDensity()*q.liquid+m.gasDensity()*q.gas;
        peak=std::max(peak,normalized(m.atomIntegrals()[i],f,m.liquidDensity()*std::abs(q.liquid)+m.gasDensity()*std::abs(q.gas),1e-12,"single spacetime density weighting"));}
    for(int u=0;u<l.velocityCount();++u){long double sum=0,scale=0;for(int i:r.primalAtoms()[u]){sum+=m.atomIntegrals()[i];scale+=std::abs(m.atomIntegrals()[i]);}
        peak=std::max(peak,normalized(m.faceIntegrals()[u],sum,scale,1e-12,"spacetime primal subface sum"));
        const auto f=l.domain().face(mesh,l.nativeFace(u));for(auto side:{std::pair<int,int>{f.negativeCell,1},{f.positiveCell,-1}})if(side.first>=0){const int p=l.pressureDof(side.first);outward[p]+=side.second*sum;gross[p]+=scale;}}
    for(int p=0;p<l.pressureCount();++p){const auto a=m.beginGeometry().cellVolumes()[p],b=m.endGeometry().cellVolumes()[p];
        const auto dm=m.liquidDensity()*(b.liquid-a.liquid)+m.gasDensity()*(b.gas-a.gas);
        peak=std::max(peak,normalized(m.cells()[p].integratedDefect,dm+outward[p],std::abs(dm)+gross[p],1e-12,"spacetime independent primal R"));}
    for(std::size_t k=0;k<g.segments().size();++k){long double sum=0,scale=0;for(int i:r.dualAtoms()[k]){sum+=m.atomIntegrals()[i];scale+=std::abs(m.atomIntegrals()[i]);}
        peak=std::max(peak,normalized(dual.faceIntegrals()[k],sum,scale,1e-12,"spacetime dual subface sum"));}
    const bool rz=l.domain().axisymmetric();const auto field=[rz](Vector2D p){return Vector2D{rz?p.x:1+p.x+p.y,.3+p.x-p.y};};
    std::vector<double> v(l.velocityCount()),traces(g.segments().size());
    for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(mesh,l.nativeFace(u));const auto q=field(f.center);v[u]=f.horizontal?q.y:q.x;}
    for(std::size_t k=0;k<g.segments().size();++k){const auto& e=g.segments()[k];if(e.negativeVelocity>=0&&e.positiveVelocity>=0)continue;
        const int u=e.negativeVelocity>=0?e.negativeVelocity:e.positiveVelocity;const auto q=field((e.begin+e.end)*.5);
        traces[k]=l.normalAxis(u)==CchDofLayout::NormalAxis::First?q.x:q.y;}
    const auto values=CchNativeValues(l,v,std::vector<double>(l.pressureCount()));
    const auto boundary=CchConvectionBoundary::bind(mesh,dual,traces,dual.sources(),dual.timeLayer());
    const auto result=CchKepConvection::build(mesh,dual).apply(mesh,values,boundary,dual.sources(),dual.timeLayer());
    std::array<long double,2> kinetic{},massPower{},materialPower{},boundaryPower{},work{};
    for(int u=0;u<l.velocityCount();++u){const int d=l.normalAxis(u)==CchDofLayout::NormalAxis::First?0:1;
        const auto a=dual.beginMass().phaseVolumes()[u],b=dual.endMass().phaseVolumes()[u];
        const auto dm=m.liquidDensity()*(b.liquid-a.liquid)+m.gasDensity()*(b.gas-a.gas);long double flow=0,scale=0;
        for(int side=0;side<4;++side)for(int k:g.volumes()[u].sides[side]){flow+=(side%2?1:-1)*m.dualFaceIntegrals()[k];scale+=std::abs(m.dualFaceIntegrals()[k]);}
        const auto R=(dm+flow)/m.deltaTime(),mdot=dm/m.deltaTime();const long double velocity=v[u];
        peak=std::max(peak,normalized(dual.cells()[u].integratedDefect,dm+flow,std::abs(dm)+scale,1e-12,"spacetime independent dual R"));
        const auto N=result.momentumRates()[u],load=result.momentumLoads()[u];
        peak=std::max(peak,normalized(load,-N+velocity*R,std::abs(N)+std::abs(velocity)*(std::abs(dm)+scale)/m.deltaTime(),1e-12,"spacetime independent uR load"));
        kinetic[d]+=velocity*load;massPower[d]+=.5L*velocity*velocity*mdot;materialPower[d]+=.5L*velocity*velocity*R;
        work[d]+=std::abs(velocity*load)+velocity*velocity*(std::abs(dm)+scale)/m.deltaTime();}
    for(std::size_t k=0;k<g.segments().size();++k){const auto& e=g.segments()[k];if(e.negativeVelocity>=0&&e.positiveVelocity>=0)continue;
        const int u=e.negativeVelocity>=0?e.negativeVelocity:e.positiveVelocity,d=l.normalAxis(u)==CchDofLayout::NormalAxis::First?0:1;
        const auto f=(e.negativeVelocity>=0?1:-1)*dual.faceRates()[k],b=f*(v[u]*traces[k]-.5L*v[u]*v[u]);boundaryPower[d]+=b;work[d]+=std::abs(b);}
    for(int d=0;d<2;++d)peak=std::max(peak,normalized(kinetic[d]-massPower[d],materialPower[d]-boundaryPower[d],work[d],1e-12,"spacetime independent material KE/skew"));
    std::cout<<"SPACETIME_IDENTITIES,"<<rz<<','<<m.liquidDensity()<<','<<peak<<",PASS"<<std::endl;
    // H11: the controlled operator consumes this actual integrated F and R.
    // The transported probe is independent of the prescribed material motion.
    std::vector<double> w(v.size());
    for(int u=0;u<l.velocityCount();++u){const auto p=l.domain().face(mesh,l.nativeFace(u)).center;
        w[u]=l.normalAxis(u)==CchDofLayout::NormalAxis::First?(rz?p.x:1)*std::sin(57*p.x+43*p.y):std::cos(61*p.x-47*p.y);}
    const CchNativeValues probe(l,w,std::vector<double>(l.pressureCount()));
    const auto controlled=CchKepConvection::buildControlled(mesh,dual,{}).apply(mesh,probe,boundary,dual.sources(),dual.timeLayer());
    std::vector<long double> L(w.size());std::array<long double,2> dissipation{},balance{},source{},boundaryWork{},scale{};
    long double controlledError=0;
    for(const auto& e:controlled.dissipationSamples()){
        const auto& s=g.segments()[e.firstSegment];check(s.negativeVelocity>=0 && s.positiveVelocity>=0,"spacetime stabilization must use actual internal segment");
        const int i=s.negativeVelocity,j=s.positiveVelocity,d=l.normalAxis(i)==CchDofLayout::NormalAxis::First?0:1;
        const auto flux=dual.faceIntegrals()[e.firstSegment]/m.deltaTime();const auto jump=static_cast<long double>(w[i])-w[j];
        const auto term=.5L*e.chi*std::abs(flux)*jump;L[i]+=term;L[j]-=term;dissipation[d]+=term*jump;scale[d]+=std::abs(term*w[i])+std::abs(term*w[j]);}
    for(int u=0;u<l.velocityCount();++u){const int d=l.normalAxis(u)==CchDofLayout::NormalAxis::First?0:1;
        const auto a=dual.beginMass().phaseVolumes()[u],b=dual.endMass().phaseVolumes()[u];
        const auto mdot=(m.liquidDensity()*(b.liquid-a.liquid)+m.gasDensity()*(b.gas-a.gas))/m.deltaTime();long double flow=0;
        for(int side=0;side<4;++side)for(int k:g.volumes()[u].sides[side])flow+=(side%2?1:-1)*dual.faceIntegrals()[k]/m.deltaTime();
        const long double velocity=w[u],load=controlled.momentumLoads()[u],N=controlled.momentumRates()[u];
        controlledError=std::max({controlledError,normalized(controlled.dissipationRates()[u],L[u],std::max(1.L,std::abs(L[u])),1e-12,"spacetime actual stabilization"),
            normalized(load,-N+velocity*(mdot+flow),std::max(1.L,std::abs(N)+std::abs(velocity)*(std::abs(mdot)+std::abs(flow))),1e-12,"spacetime controlled real uR")});
        balance[d]+=velocity*load-.5L*velocity*velocity*mdot;source[d]+=.5L*velocity*velocity*(mdot+flow);
        scale[d]+=std::abs(velocity*load)+velocity*velocity*(std::abs(mdot)+std::abs(flow));}
    for(std::size_t k=0;k<g.segments().size();++k){const auto& e=g.segments()[k];if(e.negativeVelocity>=0&&e.positiveVelocity>=0)continue;
        const int u=e.negativeVelocity>=0?e.negativeVelocity:e.positiveVelocity,d=l.normalAxis(u)==CchDofLayout::NormalAxis::First?0:1;
        const auto power=(e.negativeVelocity>=0?1:-1)*dual.faceRates()[k]*(w[u]*traces[k]-.5L*w[u]*w[u]);boundaryWork[d]+=power;scale[d]+=std::abs(power);}
    for(int d=0;d<2;++d){check(dissipation[d]>=0,"negative spacetime stabilization");controlledError=std::max(controlledError,
        normalized(balance[d],source[d]-boundaryWork[d]-dissipation[d],std::max(1.L,scale[d]),1e-12,"spacetime real material KE with controlled dissipation"));}
    std::cout<<"SPACETIME_CONTROLLED,"<<rz<<','<<m.liquidDensity()<<','<<controlled.dissipationSamples().size()<<','<<dissipation[0]<<','<<dissipation[1]<<','<<controlledError<<",PASS"<<std::endl;
}
void spaceTimeDenseField(const QuadtreeMesh& mesh,const CchCharacteristicSource& motion)
{
    const auto& r=motion.geometry().reconstruction();const auto& a=motion.previous();
    const auto other=mmsValues(mesh,r.momentumGeometry().layout(),3);
    const auto b=r.evaluate(mesh,other,motion.arrival().geometryIdentity(),motion.arrival().timeLayer());
    long double peak=0;int rows=0;
    for(double theta:{.125,.5,.875,1.}){auto clock=b.timeLayer();clock.end=clock.begin+theta*(clock.end-clock.begin);
        const auto dense=r.interpolate(mesh,a,b,b.geometryIdentity(),clock);
        for(std::size_t k=0;k<r.pieces().size();++k){const auto& box=r.pieces()[k];const auto h=box.upper-box.lower;
            const auto x=box.lower+Vector2D{.371*h.x,.613*h.y};const auto left=a.pieceVelocityJet(static_cast<int>(k),x).jet,right=b.pieceVelocityJet(static_cast<int>(k),x).jet,q=dense.pieceVelocityJet(static_cast<int>(k),x).jet;
            const std::array<double,6> lv{left.value.x,left.value.y,left.derivative.xx*h.x,left.derivative.xy*h.y,left.derivative.yx*h.x,left.derivative.yy*h.y};
            const std::array<double,6> rv{right.value.x,right.value.y,right.derivative.xx*h.x,right.derivative.xy*h.y,right.derivative.yx*h.x,right.derivative.yy*h.y};
            const std::array<double,6> actual{q.value.x,q.value.y,q.derivative.xx*h.x,q.derivative.xy*h.y,q.derivative.yx*h.x,q.derivative.yy*h.y};
            long double work=0;for(int d=0;d<6;++d)work+=std::abs(lv[d])+std::abs(rv[d]);
            for(int d=0;d<6;++d){peak=std::max(peak,normalized(actual[d],(1-theta)*lv[d]+theta*rv[d],work,1e-12,"same endpoint field value/J"));++rows;}}
    }
    auto invalid=b.timeLayer();invalid.begin=-.1;
    rejects([&]{r.interpolate(mesh,a,b,b.geometryIdentity(),invalid);},"dense geometry clock outside endpoint slab");
    std::cout<<"SPACETIME_DENSE_FIELD,"<<rows<<','<<peak<<",PASS"<<std::endl;
}
void spaceTimeGuards(const QuadtreeMesh& mesh,const HermiteGalsPublication& original,const CchCharacteristicSource& motion,
    CchSpaceTimeConfig config,const HermiteGalsBoundaryFirstJetQuery& boundary,const HermiteGalsPublicationConfig& publication)
{
    const auto saved=original.state().primary();const auto id=original.state().identity();int guards=0;
    const auto run=[&](auto f,const char* why){rejects(f,why);++guards;};
    const auto prepare=[&](CchSpaceTimeConfig c){return CchMaterialTransport::prepareSpaceTime(mesh,original,motion,1000,1,c,boundary,publication);};
    auto bad=config;bad.maxGeometrySamples=3;run([&]{prepare(bad);},"geometry sample budget exhausted");
    bad=config;bad.maxDepth=0;run([&]{prepare(bad);},"refinement depth exhausted");
    bad=config;bad.relativeTolerance=0;run([&]{prepare(bad);},"invalid spacetime quadrature budget");
    bad=config;bad.characteristic.maxStep=0;run([&]{prepare(bad);},"invalid trace work budget");
    Case shape;shape.pureFirstJet=true;shape.radial=original.state().primary().measure==HermiteGalsMeasure::AxisymmetricVolume;shape.planar=true;shape.planeOffset=2;
    const auto foreign=initial(mesh,shape);
    run([&]{CchMaterialTransport::prepareSpaceTime(mesh,foreign,motion,1000,1,config,boundary,publication);},"different velocity geometry identity");
    const auto& after=original.state().primary();check(original.state().identity()==id&&after.phaseVolumes.empty()&&after.vertices.size()==saved.vertices.size()&&after.edges.size()==saved.edges.size(),"failed spacetime changed source");
    for(std::size_t i=0;i<saved.vertices.size();++i){const auto a=saved.vertices[i],b=after.vertices[i];check(a.vertex==b.vertex&&a.phi==b.phi&&a.dx==b.dx&&a.dy==b.dy,"failed spacetime changed primary vertex");}
    for(std::size_t i=0;i<saved.edges.size();++i){const auto a=saved.edges[i],b=after.edges[i];check(a.vertices==b.vertices&&a.phi==b.phi&&a.dn==b.dn,"failed spacetime changed primary edge");}
    std::cout<<"SPACETIME_GUARDS,"<<guards<<",PASS"<<std::endl;
}
void spaceTimeMaterial(int onlyRz=-1,int onlyAmr=-1,bool coarseFineCrossing=false)
{
    for(bool rz:{false,true})for(int amr:{0,1}) {
        if((onlyRz>=0&&int(rz)!=onlyRz)||(onlyAmr>=0&&amr!=onlyAmr))continue;
        auto mesh=makeMesh(4,amr);Case s;s.pureFirstJet=true;s.radial=rz;s.planar=true;s.planeOffset=coarseFineCrossing?.23137:.48137;s.motion=1;s.speed=.2;
        const double dt=.13;const auto original=initial(*mesh,s),identity=original;
        const auto metric=metricsFor(*mesh,rz);const auto g=CchMomentumGeometry::build(*mesh,metric);
        const auto reconstruction=CchAdvectorReconstruction::prepare(*mesh,g);
        const auto u=velocity(*mesh,g.layout(),rz?Vector2D{0,.2}:Vector2D{.2,0});
        const auto a=reconstruction.evaluate(*mesh,u,original.state().identity(),original.state().timeLayer());
        const auto b=reconstruction.evaluate(*mesh,u,original.state().identity(),{1,0,dt,Position::End});
        const auto motion=CchCharacteristicGeometry::prepare(*mesh,reconstruction).linear(*mesh,a,b);
        CchSpaceTimeConfig config;config.characteristic.maxStep=.02;
        // The explicit outer scalar data declare the plane's transverse symmetry.
        HermiteGalsPublicationConfig publication;publication.volume.symmetryPlanes=s.symmetry(*mesh);
        const HermiteGalsBoundaryFirstJetQuery boundary=[s](Vector2D p,HermiteGalsTimeLayer t,Boundary){const auto j=s.profile(t.time())(p);return HermiteGalsBoundaryFirstJet{j.value,{j.dx,j.dy}};};
        for(double density:{1.,1000.,1e6}) {
            const auto material=CchMaterialTransport::prepareSpaceTime(*mesh,original,motion,density,1,config,boundary,publication);
            check(material.isSpaceTime()&&material.characteristicSource().sameIdentity(motion),"spacetime motion identity");
            check(original.state().identity()==identity.state().identity()&&original.state().primary().phaseVolumes.empty(),"spacetime modified source inventory");
            long double peak=0;
            for(std::size_t i=0;i<reconstruction.atoms().size();++i){const auto& atom=reconstruction.atoms()[i];
                const double speed=atom.normalAxis==(rz?1:0)?.2:0;
                const double normal=rz?atom.begin.y:atom.begin.x;
                const double liquidTime=std::clamp(dt-(normal-s.planeOffset)/.2,0.,dt);
                const auto expected=atom.area*speed*liquidTime;
                peak=std::max(peak,normalized(material.phaseAtomIntegrals()[i].liquid,expected,atom.area*.2*dt,1e-7,"analytic spacetime plane atom"));
                normalized(material.phaseAtomIntegrals()[i].liquid+material.phaseAtomIntegrals()[i].gas,atom.area*speed*dt,atom.area*.2*dt,1e-12,"spacetime total atom");
            }
            if(coarseFineCrossing){
                int crossed=0;long double facePeak=0;
                for(int f=0;f<g.layout().velocityCount();++f){const auto face=g.layout().domain().face(*mesh,g.layout().nativeFace(f));
                    if(face.negativeCell<0||face.positiveCell<0||face.horizontal!=rz)continue;
                    const auto& left=mesh->cells(mesh->cellIndexToLabel[face.negativeCell]);
                    const auto& right=mesh->cells(mesh->cellIndexToLabel[face.positiveCell]);
                    if(left.level()==right.level())continue;
                    const double normal=rz?face.center.y:face.center.x;
                    if(!(normal>s.planeOffset&&normal<s.planeOffset+s.speed*dt))continue;
                    check(material.beginGeometry().faceMeasures()[f].liquid==0&&material.endGeometry().faceMeasures()[f].gas==0,"coarse/fine face did not change phase");
                    const double liquidTime=dt-(normal-s.planeOffset)/s.speed;
                    const auto area=metric.faceMeasures()[f];
                    const auto expected=area*s.speed*(density*liquidTime+(dt-liquidTime));
                    facePeak=std::max(facePeak,normalized(material.faceIntegrals()[f],expected,area*s.speed*dt*density,1e-7,"moving material crossing actual coarse/fine face"));
                    ++crossed;
                }
                check(crossed>0,"no nonzero coarse/fine phase flux transition checked");
                std::cout<<"SPACETIME_COARSE_FINE,"<<rz<<','<<density<<','<<s.planeOffset<<','<<crossed<<','<<facePeak<<",PASS"<<std::endl;
            }
            const auto mass=[&](const CchPhaseGeometry& p){return CchMomentumMass::build(*mesh,g,p,density,1);};
            const auto dual=CchDualMaterial::prepare(*mesh,material,mass(material.beginGeometry()),mass(material.midpointRate().geometry()),mass(material.endGeometry()));
            check(dual.faceIntegrals()==material.dualFaceIntegrals(),"dual did not consume original time integral");
            spaceTimeIdentities(*mesh,material,dual);
            auto wrongClock=material.timeLayer();wrongClock.begin=-.01;
            rejects([&]{material.validateFor(*mesh,material.sources(),wrongClock);},"time layer");
            auto wrongSource=material.sources();wrongSource.begin=material.sources().end;
            rejects([&]{material.validateFor(*mesh,wrongSource,material.timeLayer());},"identity");
            const auto impulse=dual.consistencyImpulse(*mesh,u,dual.sources(),dual.timeLayer());
            for(std::size_t i=0;i<impulse.size();++i)normalized(impulse[i],u.normalVelocity()[i]*dual.cells()[i].integratedDefect,std::max(1.L,std::abs(impulse[i])),1e-12,"spacetime uR identity");
            std::cout<<"SPACETIME_PLANE,"<<rz<<','<<amr<<','<<density<<','<<material.spaceTimeReport().geometrySamples<<','<<material.spaceTimeReport().acceptedPanels<<','<<peak<<','<<material.report().l1IntegratedDefect/dt<<','<<dual.report().components[0].l1IntegratedDefect/dt<<','<<dual.report().components[1].l1IntegratedDefect/dt<<",PASS"<<std::endl;
        }
        if(coarseFineCrossing)continue; // Source/regrid guards run in the original four layouts.
        spaceTimeGuards(*mesh,original,motion,config,boundary,publication);
        spaceTimeDenseField(*mesh,motion);
        // Reuse the balanced nested layouts from the FirstJet transfer contract.
        auto target=makeMesh(4,amr?2:1);
        rejects([&]{CchMaterialTransport::prepareSpaceTime(*target,original,motion,1000,1,config,boundary,publication);},"mesh");
        const auto transferred=transferHermiteGalsFirstJetMesh(*mesh,original,*target,HermiteGalsMeshTransferKind::ExactRefinement,publication);
        const auto& source=transferred.publication();
        const auto nextGeometry=CchMomentumGeometry::build(*target,metricsFor(*target,rz));
        const auto nextReconstruction=CchAdvectorReconstruction::prepare(*target,nextGeometry);
        const auto nextVelocity=velocity(*target,nextGeometry.layout(),rz?Vector2D{0,.2}:Vector2D{.2,0});
        const auto nextBegin=nextReconstruction.evaluate(*target,nextVelocity,source.state().identity(),source.state().timeLayer());
        const auto nextEnd=nextReconstruction.evaluate(*target,nextVelocity,source.state().identity(),{1,0,.01,Position::End});
        const auto nextMotion=CchCharacteristicGeometry::prepare(*target,nextReconstruction).linear(*target,nextBegin,nextEnd);
        const auto resumed=CchMaterialTransport::prepareSpaceTime(*target,source,nextMotion,1000,1,config,boundary,publication);
        const auto nextMass=[&](const CchPhaseGeometry& p){return CchMomentumMass::build(*target,nextGeometry,p,1000,1);};
        const auto nextDual=CchDualMaterial::prepare(*target,resumed,nextMass(resumed.beginGeometry()),nextMass(resumed.midpointRate().geometry()),nextMass(resumed.endGeometry()));
        check(source.state().identity()!=original.state().identity()&&source.state().primary().phaseVolumes.empty(),"regrid reused old geometry or inventory");
        spaceTimeIdentities(*target,resumed,nextDual);
        std::cout<<"SPACETIME_REGRID,"<<rz<<','<<amr<<','<<transferred.report().refinedCells<<",PASS"<<std::endl;
    }
}
void spaceTimeEndpoints()
{
    int scenes=0;long double peak=0;
    for(bool rz:{false,true})for(int amr:{0,1})for(bool atEnd:{false,true}){
        auto mesh=makeMesh(4,amr);const double dt=.125,speed=.2;
        Case shape;shape.pureFirstJet=true;shape.planar=true;shape.radial=rz;shape.motion=1;shape.speed=speed;shape.planeOffset=atEnd?.5-speed*dt:.5;
        const auto original=initial(*mesh,shape);const auto g=CchMomentumGeometry::build(*mesh,metricsFor(*mesh,rz));
        const auto r=CchAdvectorReconstruction::prepare(*mesh,g);const auto u=velocity(*mesh,g.layout(),rz?Vector2D{0,speed}:Vector2D{speed,0});
        const auto a=r.evaluate(*mesh,u,original.state().identity(),original.state().timeLayer()),b=r.evaluate(*mesh,u,original.state().identity(),{1,0,dt,Position::End});
        const auto motion=CchCharacteristicGeometry::prepare(*mesh,r).linear(*mesh,a,b);
        CchSpaceTimeConfig c;c.characteristic.maxStep=.0125;HermiteGalsPublicationConfig p;p.volume.symmetryPlanes=shape.symmetry(*mesh);
        const HermiteGalsBoundaryFirstJetQuery boundary=[shape](Vector2D x,HermiteGalsTimeLayer t,Boundary){const auto j=shape.profile(t.time())(x);return HermiteGalsBoundaryFirstJet{j.value,{j.dx,j.dy}};};
        const auto material=CchMaterialTransport::prepareSpaceTime(*mesh,original,motion,1000,1,c,boundary,p);
        for(std::size_t i=0;i<r.atoms().size();++i){const auto& atom=r.atoms()[i];const auto normal=rz?atom.begin.y:atom.begin.x;
            const auto expected=atom.normalAxis==(rz?1:0)?atom.area*speed*std::clamp(dt-(normal-shape.planeOffset)/speed,0.,dt):0;
            peak=std::max(peak,normalized(material.phaseAtomIntegrals()[i].liquid,expected,atom.area*speed*dt,1e-7,"endpoint crossing spacetime integral"));}
        ++scenes;
    }
    std::cout<<"SPACETIME_ENDPOINTS,"<<scenes<<','<<peak<<",PASS"<<std::endl;
}
int main(int argc,char** argv)
{
    try{
        if(argc==3&&std::string(argv[1])=="--spacetime-amr-crossing"){
            check(std::string(argv[2])=="0"||std::string(argv[2])=="1","invalid spacetime coordinate system");
            std::cout<<std::setprecision(17);spaceTimeMaterial(std::stoi(argv[2]),1,true);return 0;}
        if(argc==2&&std::string(argv[1])=="--spacetime-endpoints"){std::cout<<std::setprecision(17);spaceTimeEndpoints();return 0;}
        if((argc==2||argc==4)&&std::string(argv[1])=="--spacetime"){
            if(argc==4)check((std::string(argv[2])=="0"||std::string(argv[2])=="1")&&(std::string(argv[3])=="0"||std::string(argv[3])=="1"),"invalid spacetime layout");
            std::cout<<std::setprecision(17);spaceTimeMaterial(argc==4?std::stoi(argv[2]):-1,argc==4?std::stoi(argv[3]):-1);return 0;}
        check(argc==1,"unknown material contract arguments");
        std::cout<<std::setprecision(17);impermeableMaterial();
        for(bool rz:{false,true})for(int amr:{0,1})for(bool subset:{false,true}){
            auto mesh=makeMesh(4,amr);const auto metrics=metricsFor(*mesh,rz,subset);
            for(int kind=0;kind<4;++kind){Case s;s.pureFirstJet=true;s.radial=rz;s.amr=amr!=0;s.planar=kind!=0;s.planeOffset=kind==1?.43:kind==2?2:-1;
                const auto a=CchPhaseGeometry::build(*mesh,metrics,publishAt(*mesh,s,{})),m=CchPhaseGeometry::build(*mesh,metrics,publishAt(*mesh,s,{1,0,.04,Position::Midpoint})),b=CchPhaseGeometry::build(*mesh,metrics,publishAt(*mesh,s,{1,0,.04,Position::End}));
                for(double rho:{1.,1000.,1e6})for(double sign:{0.,1.,-1.}){
                    const auto u=velocity(*mesh,metrics.layout(),rz?Vector2D{0,sign*.2}:Vector2D{sign*.3,-sign*.2});
                    const auto material=CchMaterialTransport::prepare(*mesh,a,rate(*mesh,m,u),b,rho,1);inspect(*mesh,s,material,u);
                    if(sign==0)check(material.report().l1IntegratedDefect==0,"static zero velocity defect not exact zero");
                }
                if(kind==1)singleFace(*mesh,s,a,m,b,subset);
            }
        }
        for(bool rz:{false,true})for(bool gas:{false,true}){
            auto mesh=makeMesh(4);Case s;s.pureFirstJet=true;s.planar=true;s.radial=rz;s.planeOffset=gas?1-std::ldexp(1.,-40):std::ldexp(1.,-40);
            const auto metrics=metricsFor(*mesh,rz);const auto u=velocity(*mesh,metrics.layout(),{.3,-.2});
            const auto material=staticMaterial(*mesh,metrics,s,u,gas?1:1e6,gas?1e6:1);inspect(*mesh,s,material,u);
            int tiny=0;for(int f=0;f<metrics.layout().velocityCount();++f){const auto phase=material.midpointRate().geometry().faceMeasures()[f];const auto small=gas?phase.gas:phase.liquid;
                if(small>0 && small<1e-10 && u.normalVelocity()[f]!=0){check(material.faceRates()[f]!=material.midpointRate().totalFaceRates()[f],"small phase material contribution vanished");++tiny;}}
            check(tiny>0,"no tiny phase mass flux checked");std::cout<<"TINY,"<<rz<<','<<gas<<','<<tiny<<'\n';
        }
        for(bool rz:{false,true})for(int amr:{0,1}){crossing(rz,amr);spatial(rz,amr);temporal(rz,amr);}
        guards();
        spatialMaterial();
        std::cout<<"SUMMARY,"<<scenarios<<','<<cellChecks<<','<<faceChecks<<','<<internalChecks<<','<<coarseChecks<<','<<rejected<<','<<spatialRows<<','<<temporalRows<<','<<peakMass<<','<<peakFlux<<','<<peakDefect<<','<<peakBalance<<','<<peakCoarse<<'\n';
        std::cout<<"H08 CCH midpoint material flux and nonconservative defect PASS\n";
        referencePolynomialContract();
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
