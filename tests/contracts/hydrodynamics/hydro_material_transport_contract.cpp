#include "Modules/Hydrodynamics/HydroMaterialTransport.h" // Public header stands alone.
#include "Modules/Hydrodynamics/HydroKepConvection.h"
#include "interface_first_jet_fixture.h"
#include <algorithm>
#include <chrono>
#include <cfenv>
#if defined(_MSC_VER) || defined(__SSE__)
#include <xmmintrin.h>
#endif
#include "Core/Mesh/CartesianMeasure.h"
#include <map>
#include <mutex>
#include <variant>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace
{
using namespace LAMP::Hydrodynamics;
int scenarios=0,cellChecks=0,faceChecks=0,internalChecks=0,coarseChecks=0,rejected=0,spatialRows=0,temporalRows=0;
long double peakMass=0,peakFlux=0,peakDefect=0,peakBalance=0,peakCoarse=0;
static_assert(!std::is_default_constructible_v<HydroMaterialTransport>);
static_assert(!std::is_default_constructible_v<HydroMaterialFlux>);
static_assert(std::is_same_v<decltype(std::declval<HydroMaterialFlux&>().atomRates()),const std::vector<long double>&>);
static_assert(std::is_same_v<decltype(std::declval<HydroMaterialTransport&>().cells()),const std::vector<HydroMaterialCell>&>);
HydroPhysicalDivergence metricsFor(const QuadtreeMesh& mesh,bool rz,bool subset=false)
{
    FaceSegmentGraph graph;graph.build(mesh,[subset](const QuadtreeMesh&,const Cell& c){return !subset || c.centerCoord().x<.5;});
    return HydroPhysicalDivergence::build(mesh,HydroIntegralDivergence::build(mesh,HydroDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,rz))));
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
HydroNativeValues velocity(const QuadtreeMesh& mesh,const HydroDofLayout& layout,Vector2D v)
{
    std::vector<double> normal;
    for(int f=0;f<layout.velocityCount();++f){const auto face=layout.domain().face(mesh,layout.nativeFace(f));normal.push_back(face.horizontal?v.y:v.x*(layout.domain().axisymmetric()?face.center.x:1));}
    return HydroNativeValues(layout,normal,std::vector<double>(layout.pressureCount()));
}
HydroPhaseRate rate(const QuadtreeMesh& mesh,HydroPhaseGeometry g,const HydroNativeValues& u)
{
    const auto geometry=HydroMomentumGeometry::build(mesh,g.metrics());
    const auto advector=HydroAdvectorReconstruction::prepare(mesh,geometry).evaluate(mesh,u,g.geometryIdentity(),g.timeLayer());
    return HydroPhaseRate::prepareSpatial(mesh,g,advector);
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
Error inspect(const QuadtreeMesh& mesh,const Case& s,const HydroMaterialTransport& material,const HydroNativeValues& u,bool exact=true)
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
HydroMaterialTransport staticMaterial(const QuadtreeMesh& mesh,const HydroPhysicalDivergence& metrics,const Case& s,const HydroNativeValues& u,double l,double g)
{
    const auto a=HydroPhaseGeometry::build(mesh,metrics,publishAt(mesh,s,{}));
    const auto m=HydroPhaseGeometry::build(mesh,metrics,publishAt(mesh,s,{1,0,.04,Position::Midpoint}));
    const auto b=HydroPhaseGeometry::build(mesh,metrics,publishAt(mesh,s,{1,0,.04,Position::End}));
    return HydroMaterialTransport::prepare(mesh,a,rate(mesh,m,u),b,l,g);
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
HydroPhaseMeasure referenceAtom(const HydroAdvectorAtom& atom,const std::array<long double,N>& c,const Case& s,int normalSign=1)
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
    std::sort(cuts.begin(),cuts.end());HydroPhaseMeasure out;
    for(std::size_t k=1;k<cuts.size();++k){const auto l=cuts[k-1],h=cuts[k],mid=(l+h)/2,width=(h-l)/2,d=mid-lo;const auto normal=referenceValue(c,d);
        const bool liquid=zero?s.planeSign*(-normal)<0:s.exact(at(mid),0).v<0;long double value=0;
        value=referenceGauss(width,[&](long double node){const auto t=mid+width*node,offset=t-lo;const auto p=at(t);
            return std::array<long double,2>{referenceValue(c,offset),s.radial?2*pi*p.x:1.L};});
        (liquid?out.liquid:out.gas)+=normalSign*value;
    }return out;
}
void singleFace(const QuadtreeMesh& mesh,const Case& s,const HydroPhaseGeometry& a,const HydroPhaseGeometry& m,const HydroPhaseGeometry& b,bool subset)
{
    const auto& layout=m.layout();int chosen=-1;
    for(int f=0;f<layout.velocityCount();++f){const auto face=layout.domain().face(mesh,layout.nativeFace(f));const auto phase=m.faceMeasures()[f];
        if(face.negativeCell>=0 && face.positiveCell>=0 && phase.liquid>0 && phase.gas>0){chosen=f;break;}}
    check(chosen>=0,"missing shared split phase basis face");
    const auto face=layout.domain().face(mesh,layout.nativeFace(chosen));
    std::vector<double> normal(layout.velocityCount());normal[chosen]=1;
    const auto material=HydroMaterialTransport::prepare(mesh,a,rate(mesh,m,HydroNativeValues(layout,normal,std::vector<double>(layout.pressureCount()))),b,1e6,1);
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
    const auto a=HydroPhaseGeometry::build(*mesh,metrics,old),m=HydroPhaseGeometry::build(*mesh,metrics,mp),b=HydroPhaseGeometry::build(*mesh,metrics,ep);
    const auto u=velocity(*mesh,metrics.layout(),rz?Vector2D{0,s.speed}:Vector2D{s.speed,0});
    const auto material=HydroMaterialTransport::prepare(*mesh,a,rate(*mesh,m,u),b,1e6,1);
    inspect(*mesh,s,material,u);check(material.report().l1IntegratedDefect>100,"crossing defect silently removed");
    normalized(material.report().integratedDefect,0,1e6,1e-10,"crossing global mass closure");
    std::cout<<"CROSSING,"<<rz<<','<<amr<<','<<material.report().l1IntegratedDefect<<','<<material.report().maxAbsDefectRate<<'\n';
}
template<class F>void reject(F f,const std::string& why){rejects(f,why);++rejected;}
void guards()
{
    auto mesh=makeMesh(4);Case s;s.pureFirstJet=true;s.planar=true;s.planeOffset=.43;
    const auto metrics=metricsFor(*mesh,false);const auto u=velocity(*mesh,metrics.layout(),{.3,-.2});
    const auto a=HydroPhaseGeometry::build(*mesh,metrics,publishAt(*mesh,s,{}));
    const auto m=HydroPhaseGeometry::build(*mesh,metrics,publishAt(*mesh,s,{1,0,.04,Position::Midpoint}));
    const auto b=HydroPhaseGeometry::build(*mesh,metrics,publishAt(*mesh,s,{1,0,.04,Position::End}));
    const auto q=rate(*mesh,m,u);const auto material=HydroMaterialTransport::prepare(*mesh,a,q,b,1000,1);const auto saved=material.faceIntegrals();
    const auto& geometry=material.momentumGeometry();const auto& flux=material.midpointFlux();const auto copiedFlux=flux;
    const auto reused=HydroMaterialTransport::prepare(*mesh,a,copiedFlux,b);
    check(&flux.atomRates()==&copiedFlux.atomRates() && &flux.faceRates()==&reused.faceRates() && &flux.dualFaceRates()==&reused.dualFaceRates(),"material interval recomputed unique flux");
    material.validateFor(*mesh,geometry,material.sources(),material.timeLayer());
    const auto freshGeometry=HydroMomentumGeometry::build(*mesh,metrics);
    reject([&]{material.validateFor(*mesh,freshGeometry,material.sources(),material.timeLayer());},"different momentum geometry instance");
    const auto freshAdvector=HydroAdvectorReconstruction::prepare(*mesh,geometry).evaluate(*mesh,u,m.geometryIdentity(),m.timeLayer());
    for(auto foreign:{HydroAdvectorIdentity{},freshAdvector.identity()}){auto source=material.sources();source.advector=foreign;
        reject([&]{material.validateFor(*mesh,source,material.timeLayer());},"different actual advector identity");
        reject([&]{flux.validateFor(*mesh,geometry,m.geometryIdentity(),foreign,m.timeLayer());},"different actual advector identity");}
    auto wrongFluxTime=m.timeLayer();++wrongFluxTime.step;
    reject([&]{flux.validateFor(*mesh,geometry,m.geometryIdentity(),flux.advectorIdentity(),wrongFluxTime);},"different full time layer");
    reject([&]{flux.validateFor(*mesh,geometry,a.geometryIdentity(),flux.advectorIdentity(),m.timeLayer());},"different geometry identity");
    for(const auto* p:std::array<const HydroPhaseGeometry*,3>{{&a,&m,&b}}){const auto instant=HydroMaterialFlux::prepare(*mesh,rate(*mesh,*p,u),1000,1);
        check(instant.timeLayer().position==p->timeLayer().position && instant.geometryIdentity()==p->geometryIdentity(),"instantaneous flux retagged actual layer");}
    for(double bad:{0.,-1.,std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::infinity()}){
        reject([&]{HydroMaterialTransport::prepare(*mesh,a,q,b,bad,1);},"nonpositive/nonfinite density");
        reject([&]{HydroMaterialTransport::prepare(*mesh,a,q,b,1,bad);},"nonpositive/nonfinite density");
    }
    reject([&]{HydroMaterialTransport::prepare(*mesh,a,rate(*mesh,b,u),b,1,1);},"requires Midpoint/End");
    reject([&]{HydroMaterialTransport::prepare(*mesh,a,q,m,1,1);},"requires Midpoint/End");
    const auto wrong=HydroPhaseGeometry::build(*mesh,metrics,publishAt(*mesh,s,{1,-.1,.14,Position::Midpoint}));
    reject([&]{HydroMaterialTransport::prepare(*mesh,a,rate(*mesh,wrong,u),b,1,1);},"different full candidate interval");
    for(auto clock:{HermiteGalsTimeLayer{3,0,0,Position::Initial},HermiteGalsTimeLayer{0,.01,.01,Position::Initial},HermiteGalsTimeLayer{1,-.04,.04,Position::Midpoint}}){
        const auto wrongBegin=HydroPhaseGeometry::build(*mesh,metrics,publishAt(*mesh,s,clock));
        reject([&]{HydroMaterialTransport::prepare(*mesh,wrongBegin,q,b,1,1);},"previous Initial/End");
    }
    const auto fresh=HydroPhaseGeometry::build(*mesh,metricsFor(*mesh,false),a.publication());
    reject([&]{HydroMaterialTransport::prepare(*mesh,fresh,q,b,1,1);},"different native layout");
    const auto otherMetrics=HydroPhysicalDivergence::build(*mesh,HydroIntegralDivergence::build(*mesh,metrics.layout()));
    const auto otherDi=HydroPhaseGeometry::build(*mesh,otherMetrics,a.publication());
    reject([&]{HydroMaterialTransport::prepare(*mesh,otherDi,q,b,1,1);},"different retained D_I");
    for(int which=0;which<3;++which){auto id=material.sources();if(which==0)id.begin={};if(which==1)id.midpoint={};if(which==2)id.end={};reject([&]{material.validateFor(*mesh,id,material.timeLayer());},"different geometry identity");}
    auto time=material.timeLayer();time.begin=-.01;
    reject([&]{material.validateFor(*mesh,material.sources(),time);},"different full time layer");
    time=material.timeLayer();time.step=2;reject([&]{material.validateFor(*mesh,material.sources(),time);},"different full time layer");
    time=material.timeLayer();time.position=Position::End;reject([&]{material.validateFor(*mesh,material.sources(),time);},"different full time layer");
    reject([&]{HydroMaterialTransport::prepare(*mesh,a,rate(*mesh,m,velocity(*mesh,metrics.layout(),{100,100})),b,std::numeric_limits<double>::max(),std::numeric_limits<double>::max());},"nonfinite atom material");
    reject([&]{HydroMaterialTransport::prepare(*mesh,a,q,b,std::numeric_limits<double>::denorm_min(),1);},"unrepresentable atom material rate");
    const double tinyTime=4*std::numeric_limits<double>::denorm_min();
    const auto tinyMid=HydroPhaseGeometry::build(*mesh,metrics,publishAt(*mesh,s,{1,0,tinyTime,Position::Midpoint}));
    const auto tinyEnd=HydroPhaseGeometry::build(*mesh,metrics,publishAt(*mesh,s,{1,0,tinyTime,Position::End}));
    reject([&]{HydroMaterialTransport::prepare(*mesh,a,rate(*mesh,tinyMid,u),tinyEnd,1,1);},"unrepresentable face material integral");
    const double adjacent=std::nextafter(1.,2.);
    const auto laterStart=HydroPhaseGeometry::build(*mesh,metrics,publishAt(*mesh,s,{0,1,1,Position::Initial}));
    const auto collapsedMid=HydroPhaseGeometry::build(*mesh,metrics,publishAt(*mesh,s,{1,1,adjacent,Position::Midpoint}));
    const auto laterEnd=HydroPhaseGeometry::build(*mesh,metrics,publishAt(*mesh,s,{1,1,adjacent,Position::End}));
    reject([&]{HydroMaterialTransport::prepare(*mesh,laterStart,rate(*mesh,collapsedMid,u),laterEnd,1,1);},"unrepresentable midpoint interval");
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
        long double actual=0;auto previousGeometry=HydroPhaseGeometry::build(*mesh,metrics,publishAt(*mesh,s,{}));
        for(int k=0;k<count;++k){
            const double begin=total*k/count,end=total*(k+1)/count;
            const auto mid=HydroPhaseGeometry::build(*mesh,metrics,publishAt(*mesh,s,{k+1,begin,end,Position::Midpoint}));
            const auto arrival=HydroPhaseGeometry::build(*mesh,metrics,publishAt(*mesh,s,{k+1,begin,end,Position::End}));
            const auto u=velocity(*mesh,metrics.layout(),rz?Vector2D{0,s.speed}:Vector2D{s.speed,0});
            const auto material=HydroMaterialTransport::prepare(*mesh,previousGeometry,rate(*mesh,mid,u),arrival,l,g);
            for(int f:selected)actual+=material.faceIntegrals()[f];previousGeometry=arrival;
        }
        const double e=static_cast<double>(std::abs(actual-exact)/std::abs(exact));
        if(count>4){check(e<previous,"material temporal error did not decrease");order=std::log2(previous/e);}previous=e;++temporalRows;
        std::cout<<"TIME,"<<rz<<','<<amr<<','<<count<<','<<actual<<','<<exact<<','<<e<<'\n';
    }
    check(order>=1.8,"material temporal order below gate");std::cout<<"TIME_ORDER,"<<rz<<','<<amr<<','<<order<<'\n';
}
std::array<long double,3> mmsCoefficients(const HydroAdvectorAtom& a,bool rz,int mode)
{
    const auto x=a.begin.x,y=a.begin.y;const int n=a.normalAxis;
    if(mode==0)return {};
    if(mode==1)return {n==0?(rz?-.13*x:-.3):-.35,0,0};
    if(mode==3)return rz?(n==0?std::array<long double,3>{}:std::array<long double,3>{.5*x*x-.015625,x,.5}):
        (n==0?std::array<long double,3>{y-.375,1,0}:std::array<long double,3>{});
    return rz?(n==0?std::array<long double,3>{.13*x,0,0}:std::array<long double,3>{.2+.155*x*x+.19*y,.31*x,.155}):
        (n==0?std::array<long double,3>{.3+.17*x-.23*y,-.23,0}:std::array<long double,3>{-.2+.31*x+.19*y,.31,0});
}
HydroNativeValues mmsValues(const QuadtreeMesh& mesh,const HydroDofLayout& layout,int mode)
{
    std::vector<double> values;const bool rz=layout.domain().axisymmetric();
    for(int u=0;u<layout.velocityCount();++u){const auto f=layout.domain().face(mesh,layout.nativeFace(u));const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
        auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();if(f.horizontal?a.x>b.x:a.y>b.y)std::swap(a,b);
        const auto c=mmsCoefficients({f.horizontal?1:0,-1,-1,f.center,f.center,f.area},rz,mode);long double v=c[0];
        if(rz && f.horizontal && (mode==2 || mode==3))v=mode==2?.2+.155*.5*(a.x*a.x+b.x*b.x)+.19*f.center.y:.25*(a.x*a.x+b.x*b.x)-.015625;
        values.push_back(static_cast<double>(v));}
    return HydroNativeValues(layout,values,std::vector<double>(layout.pressureCount()));
}
void spatialMaterial()
{
    int scenes=0,counterflow=0,tiny=0;std::size_t rows=0;long double peak=0;
    for(bool rz:{false,true})for(int amr:{0,1}){auto mesh=makeMesh(4,amr);const auto metrics=metricsFor(*mesh,rz);const auto geometry=HydroMomentumGeometry::build(*mesh,metrics);const auto reconstruction=HydroAdvectorReconstruction::prepare(*mesh,geometry);
        for(int kind=0;kind<10;++kind){Case s;s.pureFirstJet=true;s.radial=rz;s.planar=kind!=0;s.radius=rz?.34:.19;
            s.planeOffset=kind==1?.43:kind==2?2:kind==3?-1:kind==4?std::ldexp(1.,-40):kind==5?1-std::ldexp(1.,-40):kind<8?(amr?.125:.375):.5;s.planeSign=kind==7 || kind==9?-1:1;
            const auto a=HydroPhaseGeometry::build(*mesh,metrics,publishAt(*mesh,s,{})),m=HydroPhaseGeometry::build(*mesh,metrics,publishAt(*mesh,s,{1,0,.04,Position::Midpoint})),b=HydroPhaseGeometry::build(*mesh,metrics,publishAt(*mesh,s,{1,0,.04,Position::End}));
            for(int mode=0;mode<4;++mode){if(kind<6 && mode!=2)continue;const auto advector=reconstruction.evaluate(*mesh,mmsValues(*mesh,metrics.layout(),mode),m.geometryIdentity(),m.timeLayer());const auto phase=HydroPhaseRate::prepareSpatial(*mesh,m,advector);
                for(double rho:{1.,1000.,1e6}){const double l=kind==5?1:rho,g=kind==5?rho:1;const auto flux=HydroMaterialFlux::prepare(*mesh,phase,l,g);const auto material=HydroMaterialTransport::prepare(*mesh,a,flux,b);std::vector<long double> expected;
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
    for(bool rz:{false,true})for(int amr:{0,1})for(bool subset:{false,true}){auto mesh=makeMesh(4,amr);const auto metrics=metricsFor(*mesh,rz,subset);const auto g=HydroMomentumGeometry::build(*mesh,metrics);const auto& l=g.layout();const auto r=HydroAdvectorReconstruction::prepare(*mesh,g);
        std::vector<int> walls;std::vector<double> v(l.velocityCount());const double right=subset?.5:1;
        const auto psi=[&](Vector2D p){const double radial=rz?p.x*p.x*(right-p.x)*(right-p.x):p.x*(right-p.x);return radial*p.y*(1-p.y);};
        for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(*mesh,l.nativeFace(u));const auto ends=faceBounds(*mesh,f);const auto area=metrics.faceMeasures()[u];v[u]=area>0?(f.horizontal?-1:1)*(rz?static_cast<double>(2*pi):1)*(psi(ends[1])-psi(ends[0]))/area:0;if(f.outwardSign!=0){check(v[u]==0,"closed input native mean nonzero");walls.push_back(u);}}
        const auto law=HydroAdvectorBoundary::impermeable(*mesh,g,walls);const HydroNativeValues values(l,v,std::vector<double>(l.pressureCount()));
        for(bool cut:{false,true}){Case s;s.pureFirstJet=true;s.radial=rz;s.planar=true;s.planeOffset=cut?.43:2;
            const auto a=HydroPhaseGeometry::build(*mesh,metrics,publishAt(*mesh,s,{})),m=HydroPhaseGeometry::build(*mesh,metrics,publishAt(*mesh,s,{1,0,.04,Position::Midpoint})),b=HydroPhaseGeometry::build(*mesh,metrics,publishAt(*mesh,s,{1,0,.04,Position::End}));
            const auto raw=r.evaluate(*mesh,values,m.geometryIdentity(),m.timeLayer()),actual=r.evaluate(*mesh,values,m.geometryIdentity(),m.timeLayer(),law);check(actual.boundary()->sameIdentity(law) && actual.identity()!=raw.identity(),"bound actual provenance");
            const auto phase=HydroPhaseRate::prepareSpatial(*mesh,m,actual),rawPhase=HydroPhaseRate::prepareSpatial(*mesh,m,raw);
            const auto original=metrics.apply(*mesh,values);for(int p=0;p<l.pressureCount();++p){long double scale=0;const int native=l.nativeCell(p);for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(*mesh,l.nativeFace(u));if(f.negativeCell==native || f.positiveCell==native)scale+=std::abs(static_cast<long double>(metrics.faceMeasures()[u])*v[u]);}
                peak=std::max(peak,normalized(actual.cellDivergence()[p]*metrics.cellVolumes()[p],original[p]*metrics.cellVolumes()[p],scale,1e-12,"closed actual B retains original DI"));}
            for(double rho:{1.,1000.,1e6}){const auto flux=HydroMaterialFlux::prepare(*mesh,phase,rho,1),rawFlux=HydroMaterialFlux::prepare(*mesh,rawPhase,rho,1);const auto transport=HydroMaterialTransport::prepare(*mesh,a,flux,b);std::vector<long double> expected;
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
            const auto m=HydroPhaseGeometry::build(*mesh,metrics,publishAt(*mesh,s,{}));
            const auto one=HydroAdvectorBoundary::impermeable(*mesh,g,{chosen});
            const auto r=HydroAdvectorReconstruction::prepare(*mesh,g,one);const auto actual=r.evaluate(*mesh,input,m.geometryIdentity(),m.timeLayer(),one);
            for(int u=0;u<l.velocityCount();++u)if(u!=chosen && input.normalVelocity()[u]==0){bool varied=false;for(int i:r.primalAtoms()[u])for(auto c:actual.atomPolynomials()[i].beginCoefficients)varied=varied||c!=0;if(varied)++selectiveCounterflow;}
            for(int i:r.primalAtoms()[chosen])for(auto c:actual.atomPolynomials()[i].beginCoefficients)check(c==0,"selective wall not enforced");
            const auto phase=HydroPhaseRate::prepareSpatial(*mesh,m,actual);const auto flux=HydroMaterialFlux::prepare(*mesh,phase,1e6,1);
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
            for(std::size_t u=0;u<g.volumes().size();++u)for(int side=0;side<4;++side){const auto& ids=g.volumes()[u].sides[side];if(ids.size()<2)continue;HydroPhaseMeasure q;long double material=0;for(int i:ids){q.liquid+=phase.dualFaceRates()[i].liquid;q.gas+=phase.dualFaceRates()[i].gas;material+=flux.dualFaceRates()[i];}
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
        HydroAdvectorAtom atom;atom.normalAxis=n;atom.begin=n==0?Vector2D{.5,0}:Vector2D{0,.5};atom.end=n==0?Vector2D{.5,1}:Vector2D{1,.5};
        Case s;s.pureFirstJet=true;s.radial=rz;s.planar=true;s.planeOffset=.5;s.planeSign=phaseSign;
        // The fixture's planar direction follows its physical convention:
        // XY normal x, RZ normal z. Both normal axes also exercise true cuts.
        const bool coincident=n==(s.axial()?1:0);const auto actual=referenceAtom(atom,c,s,normalSign);HydroPhaseMeasure expected;
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
        HydroNormalPolynomial p;p.beginCoefficients={1,-2,3};p.endCoefficients={2,-4,3};p.length=1;p.radialMeasure=metric>0;p.radialTangent=metric==2;p.radiusBegin=metric==1?2:.25L;p.radiusEnd=p.radiusBegin+(p.radialTangent?1:0);
        const auto& c=fromEnd?p.endCoefficients:p.beginCoefficients;const auto actual=referenceShiftedIntegral(c,range[0],range[1],p.radialMeasure,fromEnd?p.radiusEnd:p.radiusBegin,p.radialTangent?(fromEnd?-1:1):0,sign);
        const auto expected=p.integral(fromEnd,range[0],range[1],sign);const auto error=std::abs((actual-expected)/expected);check(error<=1e-12,"quadratic reference/production regression");peak=std::max(peak,error);++legacy;
    }
    check(ranges==216 && tiny==24 && atoms==96 && legacy==24,"complete-degree reference inventory");
    std::cout<<"POLYNOMIAL_REFERENCE_SUMMARY,"<<ranges<<','<<tiny<<','<<atoms<<','<<legacy<<','<<peak<<",PASS\n";
    std::cout<<"POLYNOMIAL_REFERENCE_SECONDS,"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()<<'\n';
}

}
void spaceTimeIdentities(const QuadtreeMesh& mesh,const HydroMaterialTransport& m,const HydroDualMaterial& dual)
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
        traces[k]=l.normalAxis(u)==HydroDofLayout::NormalAxis::First?q.x:q.y;}
    const auto values=HydroNativeValues(l,v,std::vector<double>(l.pressureCount()));
    const auto boundary=HydroConvectionBoundary::bind(mesh,dual,traces,dual.sources(),dual.timeLayer());
    const auto result=HydroKepConvection::build(mesh,dual).apply(mesh,values,boundary,dual.sources(),dual.timeLayer());
    std::array<long double,2> kinetic{},massPower{},materialPower{},boundaryPower{},work{};
    for(int u=0;u<l.velocityCount();++u){const int d=l.normalAxis(u)==HydroDofLayout::NormalAxis::First?0:1;
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
        const int u=e.negativeVelocity>=0?e.negativeVelocity:e.positiveVelocity,d=l.normalAxis(u)==HydroDofLayout::NormalAxis::First?0:1;
        const auto f=(e.negativeVelocity>=0?1:-1)*dual.faceRates()[k],b=f*(v[u]*traces[k]-.5L*v[u]*v[u]);boundaryPower[d]+=b;work[d]+=std::abs(b);}
    for(int d=0;d<2;++d)peak=std::max(peak,normalized(kinetic[d]-massPower[d],materialPower[d]-boundaryPower[d],work[d],1e-12,"spacetime independent material KE/skew"));
    std::cout<<"SPACETIME_IDENTITIES,"<<rz<<','<<m.liquidDensity()<<','<<peak<<",PASS"<<std::endl;
    // H11: the controlled operator consumes this actual integrated F and R.
    // The transported probe is independent of the prescribed material motion.
    std::vector<double> w(v.size());
    for(int u=0;u<l.velocityCount();++u){const auto p=l.domain().face(mesh,l.nativeFace(u)).center;
        w[u]=l.normalAxis(u)==HydroDofLayout::NormalAxis::First?(rz?p.x:1)*std::sin(57*p.x+43*p.y):std::cos(61*p.x-47*p.y);}
    const HydroNativeValues probe(l,w,std::vector<double>(l.pressureCount()));
    const auto controlled=HydroKepConvection::buildControlled(mesh,dual,{}).apply(mesh,probe,boundary,dual.sources(),dual.timeLayer());
    std::vector<long double> L(w.size());std::array<long double,2> dissipation{},balance{},source{},boundaryWork{},scale{};
    long double controlledError=0;
    for(const auto& e:controlled.dissipationSamples()){
        const auto& s=g.segments()[e.firstSegment];check(s.negativeVelocity>=0 && s.positiveVelocity>=0,"spacetime stabilization must use actual internal segment");
        const int i=s.negativeVelocity,j=s.positiveVelocity,d=l.normalAxis(i)==HydroDofLayout::NormalAxis::First?0:1;
        const auto flux=dual.faceIntegrals()[e.firstSegment]/m.deltaTime();const auto jump=static_cast<long double>(w[i])-w[j];
        const auto term=.5L*e.chi*std::abs(flux)*jump;L[i]+=term;L[j]-=term;dissipation[d]+=term*jump;scale[d]+=std::abs(term*w[i])+std::abs(term*w[j]);}
    for(int u=0;u<l.velocityCount();++u){const int d=l.normalAxis(u)==HydroDofLayout::NormalAxis::First?0:1;
        const auto a=dual.beginMass().phaseVolumes()[u],b=dual.endMass().phaseVolumes()[u];
        const auto mdot=(m.liquidDensity()*(b.liquid-a.liquid)+m.gasDensity()*(b.gas-a.gas))/m.deltaTime();long double flow=0;
        for(int side=0;side<4;++side)for(int k:g.volumes()[u].sides[side])flow+=(side%2?1:-1)*dual.faceIntegrals()[k]/m.deltaTime();
        const long double velocity=w[u],load=controlled.momentumLoads()[u],N=controlled.momentumRates()[u];
        controlledError=std::max({controlledError,normalized(controlled.dissipationRates()[u],L[u],std::max(1.L,std::abs(L[u])),1e-12,"spacetime actual stabilization"),
            normalized(load,-N+velocity*(mdot+flow),std::max(1.L,std::abs(N)+std::abs(velocity)*(std::abs(mdot)+std::abs(flow))),1e-12,"spacetime controlled real uR")});
        balance[d]+=velocity*load-.5L*velocity*velocity*mdot;source[d]+=.5L*velocity*velocity*(mdot+flow);
        scale[d]+=std::abs(velocity*load)+velocity*velocity*(std::abs(mdot)+std::abs(flow));}
    for(std::size_t k=0;k<g.segments().size();++k){const auto& e=g.segments()[k];if(e.negativeVelocity>=0&&e.positiveVelocity>=0)continue;
        const int u=e.negativeVelocity>=0?e.negativeVelocity:e.positiveVelocity,d=l.normalAxis(u)==HydroDofLayout::NormalAxis::First?0:1;
        const auto power=(e.negativeVelocity>=0?1:-1)*dual.faceRates()[k]*(w[u]*traces[k]-.5L*w[u]*w[u]);boundaryWork[d]+=power;scale[d]+=std::abs(power);}
    for(int d=0;d<2;++d){check(dissipation[d]>=0,"negative spacetime stabilization");controlledError=std::max(controlledError,
        normalized(balance[d],source[d]-boundaryWork[d]-dissipation[d],std::max(1.L,scale[d]),1e-12,"spacetime real material KE with controlled dissipation"));}
    std::cout<<"SPACETIME_CONTROLLED,"<<rz<<','<<m.liquidDensity()<<','<<controlled.dissipationSamples().size()<<','<<dissipation[0]<<','<<dissipation[1]<<','<<controlledError<<",PASS"<<std::endl;
}
void spaceTimeDenseField(const QuadtreeMesh& mesh,const HydroCharacteristicSource& motion)
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
void spaceTimeGuards(const QuadtreeMesh& mesh,const HermiteGalsPublication& original,const HydroCharacteristicSource& motion,
    HydroSpaceTimeConfig config,const HermiteGalsBoundaryFirstJetQuery& boundary,const HermiteGalsPublicationConfig& publication)
{
    const auto saved=original.state().primary();const auto id=original.state().identity();int guards=0;
    const auto run=[&](auto f,const char* why){rejects(f,why);++guards;};
    const auto prepare=[&](HydroSpaceTimeConfig c){return HydroMaterialTransport::prepareSpaceTime(mesh,original,motion,1000,1,c,boundary,publication);};
    auto bad=config;bad.maxGeometrySamples=3;run([&]{prepare(bad);},"geometry sample budget exhausted");
    bad=config;bad.maxDepth=0;run([&]{prepare(bad);},"refinement depth exhausted");
    bad=config;bad.relativeTolerance=0;run([&]{prepare(bad);},"invalid spacetime quadrature budget");
    bad=config;bad.characteristic.maxStep=0;run([&]{prepare(bad);},"invalid trace work budget");
    Case shape;shape.pureFirstJet=true;shape.radial=original.state().primary().measure==HermiteGalsMeasure::AxisymmetricVolume;shape.planar=true;shape.planeOffset=2;
    const auto foreign=initial(mesh,shape);
    run([&]{HydroMaterialTransport::prepareSpaceTime(mesh,foreign,motion,1000,1,config,boundary,publication);},"different velocity geometry identity");
    const auto& after=original.state().primary();check(original.state().identity()==id&&after.phaseVolumes.empty()&&after.vertices.size()==saved.vertices.size()&&after.edges.size()==saved.edges.size(),"failed spacetime changed source");
    for(std::size_t i=0;i<saved.vertices.size();++i){const auto a=saved.vertices[i],b=after.vertices[i];check(a.vertex==b.vertex&&a.phi==b.phi&&a.dx==b.dx&&a.dy==b.dy,"failed spacetime changed primary vertex");}
    for(std::size_t i=0;i<saved.edges.size();++i){const auto a=saved.edges[i],b=after.edges[i];check(a.vertices==b.vertices&&a.phi==b.phi&&a.dn==b.dn,"failed spacetime changed primary edge");}
    std::cout<<"SPACETIME_GUARDS,"<<guards<<",PASS"<<std::endl;
}
// Compare the real restricted scalar/rate path with the unchanged full path
// at one time in the existing accelerating-plane fixture.
void spaceTimeLocalRateWitness(const QuadtreeMesh& mesh,const Case& shape,
    const HermiteGalsPublication& original,const HydroCharacteristicSource& motion,
    const HydroSpaceTimeConfig& config,const HermiteGalsBoundaryFirstJetQuery& boundary,
    const HermiteGalsPublicationConfig& publication)
{
    const auto& reconstruction=motion.geometry().reconstruction();const auto& atoms=reconstruction.atoms();
    const auto parent=motion.arrival().timeLayer();auto clock=parent;clock.end=.09;
    int liquid=-1,gas=-1,unsupported=-1;
    for(std::size_t i=0;i<atoms.size();++i){const auto& a=atoms[i];
        if(a.normalAxis==0&&a.begin.x==.5&&a.begin.y>=.25&&a.end.y<=.5&&liquid<0)liquid=static_cast<int>(i);
        if(a.normalAxis==1&&a.begin.y==0&&a.begin.x==0&&unsupported<0)unsupported=static_cast<int>(i);}
    check(liquid>=0&&unsupported>=0,"local rate witness lacks existing liquid/remote atoms");
    for(std::size_t i=0;i<atoms.size();++i){const auto& a=atoms[i];
        if(a.normalAxis==0&&a.begin.x>.5&&a.begin.y==atoms[liquid].begin.y&&a.end.y==atoms[liquid].end.y
            &&(gas<0||a.begin.x<atoms[gas].begin.x))gas=static_cast<int>(i);}
    check(gas>=0,"local rate witness lacks existing gas atom");
    const std::vector<int> requested{gas,liquid}; // Deliberately not native atom order.
    const auto cells=HydroPhaseRate::requiredNativeCells(mesh,reconstruction,requested);
    const HermiteGalsResolvedMotion resolved(original.state().identity(),original.state().timeLayer(),
        motion.arrival().geometryIdentity(),clock,[motion,clock,config,parent](Vector2D point){
            const auto t=motion.traceBackward(point,clock.time(),parent.begin,config.characteristic);
            std::optional<Boundary> side;
            if(t.boundary()){const auto b=*t.boundary();side=b.normalAxis==0?(b.outwardSign<0?Boundary::west:Boundary::east)
                :(b.outwardSign<0?Boundary::south:Boundary::north);}
            return LevelSetBoundaryTrace::JetResult{{t.point(),clock.time()-t.time(),side,t.boundary().has_value(),false},
                t.jacobian(),t.substeps()};});
    const auto full=prepareHermiteGalsFirstJetSample(mesh,original,resolved,boundary,publication,
        HermiteGalsBoundaryTreatment::PrescribedOuterFirstJet,HermiteGalsTimeProjection::CharacteristicPullback);
    const auto local=prepareHermiteGalsFirstJetLocalSample(mesh,original,resolved,cells,boundary,publication,
        HermiteGalsBoundaryTreatment::PrescribedOuterFirstJet,HermiteGalsTimeProjection::CharacteristicPullback);
    local.validateMesh(mesh);const auto view=local.geometry();const auto& scalar=view.snapshot();
    const auto sameClock=[&](HermiteGalsTimeLayer t){return t.step==clock.step&&t.begin==clock.begin
        &&t.end==clock.end&&t.position==clock.position;};
    check(local.identity()!=full.identity()&&view.identity()==local.identity()&&sameClock(local.timeLayer())&&sameClock(view.timeLayer())
        &&local.sourceIdentity()==original.state().identity()&&local.arrivalVelocityIdentity()==resolved.velocityIdentity(),
        "local rate sample detached actual scalar owner/source/full clock");
    check(!cells.empty()&&cells.size()<static_cast<std::size_t>(full.geometry().snapshot().nativeCellCount())
        &&local.nativeCells()==cells&&local.report().requestedNativeCells==cells.size()
        &&local.report().evaluatedNativeCells==cells.size(),"local sample did not actually restrict native Q4 evaluation");
    int active=0;for(int c=0;c<scalar.nativeCellCount();++c){const bool wanted=std::find(cells.begin(),cells.end(),c)!=cells.end();
        check(scalar.nativeCellActive(c)==wanted,"local scalar fabricated unsampled active support");if(wanted)++active;}
    check(active==static_cast<int>(cells.size()),"local scalar active-cell count changed");
    const auto fullAdvector=reconstruction.interpolate(mesh,motion.previous(),motion.arrival(),full.identity(),full.timeLayer());
    const auto localAdvector=reconstruction.interpolate(mesh,motion.previous(),motion.arrival(),local.identity(),local.timeLayer());
    const auto reference=HydroPhaseRate::sampleAtomRates(mesh,full,fullAdvector);
    const auto rates=HydroPhaseRate::sampleAtomRates(mesh,local,localAdvector,requested);
    check(rates.size()==requested.size(),"local rates did not return the requested ordered atom inventory");
    const double dt=parent.end-parent.begin,time=clock.time(),acceleration=(.24-.2)/dt;
    const double position=shape.planeOffset+.2*time+.5*acceleration*time*time,speed=.2+acceleration*time;
    long double peak=0;
    for(std::size_t k=0;k<requested.size();++k){const int i=requested[k];const auto a=reference[i],b=rates[k];
        const long double total=atoms[i].area*speed,expectedLiquid=atoms[i].begin.x<position?total:0,expectedGas=total-expectedLiquid;
        peak=std::max({peak,normalized(b.liquid,a.liquid,total,1e-12,"local/full ordered signed liquid rate"),
            normalized(b.gas,a.gas,total,1e-12,"local/full ordered signed gas rate"),
            normalized(b.liquid+b.gas,total,total,1e-12,"local original signed total rate")});
        for(const auto q:{a,b}){normalized(q.liquid,expectedLiquid,total,1e-12,"independent local plane signed liquid rate");
            normalized(q.gas,expectedGas,total,1e-12,"independent local plane signed gas rate");}
    }
    check(reference[liquid].liquid!=0&&reference[gas].gas!=0,"local subset did not exercise two distinct actual phases");
    check(reference[unsupported].liquid==0&&reference[unsupported].gas==0,"remote unsupported witness is not zero work");
    rejects([&]{HydroPhaseRate::sampleAtomRates(mesh,local,localAdvector,std::vector<int>{unsupported});},"missing local native support");
    rejects([&]{HydroPhaseRate::sampleAtomRates(mesh,local,localAdvector,std::vector<int>{liquid,liquid});},"duplicate atom");
    rejects([&]{HydroPhaseRate::sampleAtomRates(mesh,local,localAdvector,std::vector<int>{-1});},"invalid atom");
    std::cout<<"SPACETIME_LOCAL_RATE,"<<clock.time()<<','<<requested.size()<<','<<cells.size()<<','
        <<local.report().evaluatedNativeCells<<','<<local.report().stagedNativeCells<<','<<peak<<','<<reference.size()<<','
        <<full.geometry().snapshot().nativeCellCount()<<",PASS\n";
}
// Reuse only the accepted partition: the new endpoint shifts the plane's
// crossing time, so both fresh rates and further refinement are observable.
void spaceTimePreviousPlanWitness(const QuadtreeMesh& mesh,const Case& shape,
    const HermiteGalsPublication& original,const HydroCharacteristicSource& oldMotion,
    const HydroMaterialTransport& oldPlan,HydroSpaceTimeConfig config,
    const HermiteGalsPublicationConfig& publication)
{
    const auto& reconstruction=oldMotion.geometry().reconstruction();const auto& layout=oldPlan.layout();
    const auto parent=oldMotion.arrival().timeLayer();const double dt=parent.end-parent.begin,speed=.2,arrivalSpeed=.24;
    const double acceleration=(arrivalSpeed-speed)/dt;
    const auto travel=[=](double time){return speed*time+.5*acceleration*time*time;};
    const auto values=velocity(mesh,layout,{arrivalSpeed,0});
    const auto arrival=reconstruction.evaluate(mesh,values,original.state().identity(),parent);
    const auto motion=oldMotion.geometry().linear(mesh,oldMotion.previous(),arrival);
    const HermiteGalsBoundaryFirstJetQuery boundary=[=](Vector2D point,HermiteGalsTimeLayer clock,Boundary){
        return HermiteGalsBoundaryFirstJet{point.x-(shape.planeOffset+travel(clock.time())),{1,0}};};
    spaceTimeLocalRateWitness(mesh,shape,original,motion,config,boundary,publication);
    const auto cold=HydroMaterialTransport::prepareSpaceTime(mesh,original,motion,1000,1,config,boundary,publication);
    auto fullConfig=config;fullConfig.enableAtomLocalTimeIntegration=false;
    const auto full=HydroMaterialTransport::prepareSpaceTime(mesh,original,motion,1000,1,fullConfig,boundary,publication);
    std::size_t notifications=0,lastSamples=0,lastPanels=0;double lastElapsed=0;
    auto warmConfig=config;warmConfig.progress=[&](const HydroSpaceTimeProgress& p){
#ifdef _OPENMP
        check(!omp_in_parallel(),"warm plan progress entered worker team");
#endif
        check(p.timeLayer.step==parent.step&&p.timeLayer.begin==parent.begin&&p.timeLayer.end==parent.end
            &&p.timeLayer.position==Position::Midpoint,"warm plan progress detached full parent clock");
        check(p.geometrySamples>=lastSamples&&std::isfinite(p.elapsedSeconds)&&p.elapsedSeconds>=lastElapsed
            &&p.deepestPanel>=0&&std::isfinite(p.relativePhaseErrorEstimate)&&p.relativePhaseErrorEstimate>=0,
            "warm plan progress changed passive count/error/time semantics");
        ++notifications;lastSamples=p.geometrySamples;lastPanels=p.panels;lastElapsed=p.elapsedSeconds;
    };
    const auto warm=HydroMaterialTransport::prepareSpaceTime(mesh,original,motion,1000,1,warmConfig,boundary,publication,&oldPlan);
    const auto& localCost=cold.spaceTimeReport().samplingCost;const auto& warmCost=warm.spaceTimeReport().samplingCost;
    const auto& fullCost=full.spaceTimeReport().samplingCost;
    const auto actualSamplingOwner=[&](const auto& cost){return cost.sourceIdentity==original.state().identity()
        &&cost.previousAdvectorIdentity==motion.previous().identity()&&cost.arrivalAdvectorIdentity==motion.arrival().identity();};
    check(actualSamplingOwner(localCost)&&actualSamplingOwner(warmCost)&&actualSamplingOwner(fullCost),
        "atom-local sampling cost detached actual source/endpoint owners");
    check(localCost.atomLocal&&warmCost.atomLocal&&!fullCost.atomLocal
        &&localCost.localSamples>0&&localCost.localScalarPointQueries>0&&localCost.localPreparedNativeCells>0
        &&fullCost.localSamples==0&&fullCost.localAtomRateEvaluations==0,
        "declared atom-local/full routes did not perform their actual sampling work");
    check(localCost.localAtomRateEvaluations>0&&localCost.localAtomRateEvaluations<fullCost.fullAtomRateEvaluations
        &&localCost.localPreparedNativeCells+localCost.fullPreparedNativeCells<fullCost.fullPreparedNativeCells,
        "actual atom-local subset did not reduce atom/native-cell evaluation work");
    check(cold.spaceTimeReport().endpointSupportPanels>0&&full.spaceTimeReport().endpointSupportPanels>0,
        "accelerated plane lost the actual endpoint blindspot guard");
    check(!cold.spaceTimeReport().reusedTimePlan&&warm.spaceTimeReport().reusedTimePlan
        &&oldPlan.spaceTimeReport().acceptedPanels>1
        &&warm.spaceTimeReport().initialPlanPanels==oldPlan.spaceTimeReport().acceptedPanels,
        "warm plan did not reuse the actual accepted multi-panel partition");
    check(warm.spaceTimeReport().acceptedPanels>warm.spaceTimeReport().initialPlanPanels,
        "changed plane crossing did not refine the reused partition");
    check(notifications>0&&lastSamples==warm.spaceTimeReport().geometrySamples
        &&lastPanels==warm.spaceTimeReport().acceptedPanels,"warm plan progress omitted the completed global estimate");
    std::vector<long double> expected(reconstruction.atoms().size()),work(expected.size());
    long double peak=0,coldWarmDifference=0;int changed=0,crossed=0;
    for(std::size_t i=0;i<expected.size();++i){const auto& atom=reconstruction.atoms()[i];
        const long double total=atom.normalAxis==0?atom.area*travel(dt):0;
        const long double liquid=atom.normalAxis==0?atom.area*std::clamp(travel(dt)-(atom.begin.x-shape.planeOffset),0.,travel(dt)):0;
        const long double gas=total-liquid,phaseWork=atom.area*arrivalSpeed*dt;
        expected[i]=1000*liquid+gas;work[i]=1000*phaseWork;
        if(liquid>0&&gas>0)++crossed;
        for(const auto* material:{&cold,&warm,&full}){const auto q=material->phaseAtomIntegrals()[i];
            peak=std::max({peak,normalized(q.liquid,liquid,phaseWork,1e-7,"independent accelerated plane signed liquid atom"),
                normalized(q.gas,gas,phaseWork,1e-7,"independent accelerated plane signed gas atom"),
                normalized(material->atomIntegrals()[i],expected[i],work[i],1e-7,"independent accelerated plane signed material atom")});
            normalized(q.liquid+q.gas,total,phaseWork,1e-12,"accelerated plane original total phase closure");
            const auto reference=full.phaseAtomIntegrals()[i];
            normalized(q.liquid,reference.liquid,phaseWork,1e-7,"full/local signed liquid time integral");
            normalized(q.gas,reference.gas,phaseWork,1e-7,"full/local signed gas time integral");
            normalized(material->atomIntegrals()[i],full.atomIntegrals()[i],work[i],1e-7,"full/local signed material time integral");
        }
        coldWarmDifference=std::max(coldWarmDifference,std::abs(cold.atomIntegrals()[i]-warm.atomIntegrals()[i])/work[i]);
        if(std::abs(warm.atomIntegrals()[i]-oldPlan.atomIntegrals()[i])>1e-7L*work[i])++changed;
    }
    check(crossed>0&&changed>0,"warm plan witness did not change a genuine cut signed flux");
    for(const auto* material:{&cold,&warm,&full}){
        check(material->characteristicSource().sameIdentity(motion),"warm plan reused the old characteristic source");
        for(int u=0;u<layout.velocityCount();++u){long double sum=0,scale=0;
            for(int i:reconstruction.primalAtoms()[u]){sum+=expected[i];scale+=work[i];}
            peak=std::max(peak,normalized(material->faceIntegrals()[u],sum,scale,1e-7,"independent accelerated plane signed native face"));}
        for(std::size_t e=0;e<reconstruction.dualAtoms().size();++e){long double sum=0,scale=0;
            for(int i:reconstruction.dualAtoms()[e]){sum+=expected[i];scale+=work[i];}
            peak=std::max(peak,normalized(material->dualFaceIntegrals()[e],sum,scale,1e-7,"independent accelerated plane signed dual face"));}
    }
    rejects([&]{HydroMaterialTransport::prepareSpaceTime(mesh,original,motion,1000,1,fullConfig,boundary,publication,&oldPlan);},
        "previous time plan: different atom-local route");
    std::cout<<"SPACETIME_ATOM_LOCAL_TIME,"<<localCost.localSamples<<','<<fullCost.fullSamples<<','
        <<localCost.localAtomRateEvaluations<<','<<fullCost.fullAtomRateEvaluations<<','
        <<localCost.localPreparedNativeCells+localCost.fullPreparedNativeCells<<','<<fullCost.fullPreparedNativeCells<<','
        <<localCost.localScalarPointQueries<<','<<warmCost.localSamples<<",PASS\n";
    const auto foreign=initial(mesh,shape);
    const auto foreignPrevious=reconstruction.evaluate(mesh,oldMotion.previous().values(),foreign.state().identity(),foreign.state().timeLayer());
    const auto foreignArrival=reconstruction.evaluate(mesh,values,foreign.state().identity(),parent);
    const auto foreignMotion=oldMotion.geometry().linear(mesh,foreignPrevious,foreignArrival);
    rejects([&]{HydroMaterialTransport::prepareSpaceTime(mesh,foreign,foreignMotion,1000,1,config,boundary,publication,&oldPlan);},
        "previous time plan: different original publication");
    auto wrongClock=parent;wrongClock.end*=1.1;
    const auto otherArrival=reconstruction.evaluate(mesh,values,original.state().identity(),wrongClock);
    const auto otherClock=oldMotion.geometry().linear(mesh,oldMotion.previous(),otherArrival);
    rejects([&]{HydroMaterialTransport::prepareSpaceTime(mesh,original,otherClock,1000,1,config,boundary,publication,&oldPlan);},
        "previous time plan: different clock interval");
    const auto freshPrevious=reconstruction.evaluate(mesh,oldMotion.previous().values(),original.state().identity(),original.state().timeLayer());
    const auto otherSource=oldMotion.geometry().linear(mesh,freshPrevious,arrival);
    rejects([&]{HydroMaterialTransport::prepareSpaceTime(mesh,original,otherSource,1000,1,config,boundary,publication,&oldPlan);},
        "previous time plan: different previous source");
    std::cout<<"SPACETIME_PREVIOUS_PLAN,"<<oldPlan.spaceTimeReport().acceptedPanels<<','<<warm.spaceTimeReport().initialPlanPanels
        <<','<<warm.spaceTimeReport().acceptedPanels<<','<<crossed<<','<<changed<<','<<notifications<<','<<peak<<','<<coldWarmDifference<<",PASS\n";
}
void spaceTimeMaterial(int onlyRz=-1,int onlyAmr=-1,bool coarseFineCrossing=false)
{
    for(bool rz:{false,true})for(int amr:{0,1}) {
        if((onlyRz>=0&&int(rz)!=onlyRz)||(onlyAmr>=0&&amr!=onlyAmr))continue;
        auto mesh=makeMesh(4,amr);Case s;s.pureFirstJet=true;s.radial=rz;s.planar=true;s.planeOffset=coarseFineCrossing?.23137:.48137;s.motion=1;s.speed=.2;
        const double dt=.13;const auto original=initial(*mesh,s),identity=original;
        const auto metric=metricsFor(*mesh,rz);const auto g=HydroMomentumGeometry::build(*mesh,metric);
        const auto reconstruction=HydroAdvectorReconstruction::prepare(*mesh,g);
        const auto u=velocity(*mesh,g.layout(),rz?Vector2D{0,.2}:Vector2D{.2,0});
        const auto a=reconstruction.evaluate(*mesh,u,original.state().identity(),original.state().timeLayer());
        const auto b=reconstruction.evaluate(*mesh,u,original.state().identity(),{1,0,dt,Position::End});
        const auto motion=HydroCharacteristicGeometry::prepare(*mesh,reconstruction).linear(*mesh,a,b);
        HydroSpaceTimeConfig config;config.characteristic.maxStep=.02;
        // The explicit outer scalar data declare the plane's transverse symmetry.
        HermiteGalsPublicationConfig publication;publication.volume.symmetryPlanes=s.symmetry(*mesh);
        const HermiteGalsBoundaryFirstJetQuery boundary=[s](Vector2D p,HermiteGalsTimeLayer t,Boundary){const auto j=s.profile(t.time())(p);return HermiteGalsBoundaryFirstJet{j.value,{j.dx,j.dy}};};
        for(double density:{1.,1000.,1e6}) {
            const auto material=HydroMaterialTransport::prepareSpaceTime(*mesh,original,motion,density,1,config,boundary,publication);
            check(material.isSpaceTime()&&material.characteristicSource().sameIdentity(motion),"spacetime motion identity");
            if(!rz&&amr==0&&!coarseFineCrossing&&density==1000)
                spaceTimePreviousPlanWitness(*mesh,s,original,motion,material,config,publication);
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
            const auto mass=[&](const HydroPhaseGeometry& p){return HydroMomentumMass::build(*mesh,g,p,density,1);};
            const auto dual=HydroDualMaterial::prepare(*mesh,material,mass(material.beginGeometry()),mass(material.midpointRate().geometry()),mass(material.endGeometry()));
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
        rejects([&]{HydroMaterialTransport::prepareSpaceTime(*target,original,motion,1000,1,config,boundary,publication);},"mesh");
        const auto transferred=transferHermiteGalsFirstJetMesh(*mesh,original,*target,HermiteGalsMeshTransferKind::ExactRefinement,publication);
        const auto& source=transferred.publication();
        const auto nextGeometry=HydroMomentumGeometry::build(*target,metricsFor(*target,rz));
        const auto nextReconstruction=HydroAdvectorReconstruction::prepare(*target,nextGeometry);
        const auto nextVelocity=velocity(*target,nextGeometry.layout(),rz?Vector2D{0,.2}:Vector2D{.2,0});
        const auto nextBegin=nextReconstruction.evaluate(*target,nextVelocity,source.state().identity(),source.state().timeLayer());
        const auto nextEnd=nextReconstruction.evaluate(*target,nextVelocity,source.state().identity(),{1,0,.01,Position::End});
        const auto nextMotion=HydroCharacteristicGeometry::prepare(*target,nextReconstruction).linear(*target,nextBegin,nextEnd);
        const auto resumed=HydroMaterialTransport::prepareSpaceTime(*target,source,nextMotion,1000,1,config,boundary,publication);
        const auto nextMass=[&](const HydroPhaseGeometry& p){return HydroMomentumMass::build(*target,nextGeometry,p,1000,1);};
        const auto nextDual=HydroDualMaterial::prepare(*target,resumed,nextMass(resumed.beginGeometry()),nextMass(resumed.midpointRate().geometry()),nextMass(resumed.endGeometry()));
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
        const auto original=initial(*mesh,shape);const auto g=HydroMomentumGeometry::build(*mesh,metricsFor(*mesh,rz));
        const auto r=HydroAdvectorReconstruction::prepare(*mesh,g);const auto u=velocity(*mesh,g.layout(),rz?Vector2D{0,speed}:Vector2D{speed,0});
        const auto a=r.evaluate(*mesh,u,original.state().identity(),original.state().timeLayer()),b=r.evaluate(*mesh,u,original.state().identity(),{1,0,dt,Position::End});
        const auto motion=HydroCharacteristicGeometry::prepare(*mesh,r).linear(*mesh,a,b);
        HydroSpaceTimeConfig c;c.characteristic.maxStep=.0125;HermiteGalsPublicationConfig p;p.volume.symmetryPlanes=shape.symmetry(*mesh);
        const HermiteGalsBoundaryFirstJetQuery boundary=[shape](Vector2D x,HermiteGalsTimeLayer t,Boundary){const auto j=shape.profile(t.time())(x);return HermiteGalsBoundaryFirstJet{j.value,{j.dx,j.dy}};};
        const auto material=HydroMaterialTransport::prepareSpaceTime(*mesh,original,motion,1000,1,c,boundary,p);
        for(std::size_t i=0;i<r.atoms().size();++i){const auto& atom=r.atoms()[i];const auto normal=rz?atom.begin.y:atom.begin.x;
            const auto expected=atom.normalAxis==(rz?1:0)?atom.area*speed*std::clamp(dt-(normal-shape.planeOffset)/speed,0.,dt):0;
            peak=std::max(peak,normalized(material.phaseAtomIntegrals()[i].liquid,expected,atom.area*speed*dt,1e-7,"endpoint crossing spacetime integral"));}
        ++scenes;
    }
    std::cout<<"SPACETIME_ENDPOINTS,"<<scenes<<','<<peak<<",PASS"<<std::endl;
}
// Compact inputs are copied from the actual independent primary/closure,
// rather than reconstructed by querying or fitting the new subset snapshot.
Snapshot::NativeLeafSubset materialSubsetPayload(const QuadtreeMesh& mesh,
    const HermiteGalsState& state,const std::vector<int>& cells)
{
    Snapshot::NativeLeafSubset payload;payload.nativeCells=cells;
    const auto& primary=state.primary();const auto& closure=state.closure();
    std::set<int> vertices;std::set<Key> sides;
    for(int c:cells){
        const auto pos=std::lower_bound(primary.geometryCells.begin(),primary.geometryCells.end(),c);
        check(pos!=primary.geometryCells.end()&&*pos==c,"subset witness requested absent primary cell");
        payload.cellCenters.push_back(closure.cellBubbleCoefficients[static_cast<std::size_t>(pos-primary.geometryCells.begin())]);
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);
        for(int k=0;k<4;++k){
            const int a=mesh.verts(cell.cellToVert[k]).index(),b=mesh.verts(cell.cellToVert[(k+1)%4]).index();
            vertices.insert(a);Key key{{a,b}};const auto x=point(mesh,a),y=point(mesh,b);
            if(x.y==y.y?x.x>y.x:x.y>y.y)std::swap(key[0],key[1]);sides.insert(key);
        }
    }
    for(std::size_t k=0;k<primary.vertices.size();++k){const auto& v=primary.vertices[k];
        if(vertices.count(v.vertex)){payload.nativeVertices.push_back(v.vertex);
            payload.nodes.push_back({v.phi,v.dx,v.dy,closure.vertexMixedDerivatives[k]});}}
    for(const auto& edge:primary.edges)if(sides.count(edge.vertices))
        payload.edgeMidpoints.push_back({edge.vertices,edge.phi,edge.dn});
    return payload;
}
Snapshot materialSubsetGeometryWitness(const QuadtreeMesh& mesh,const Snapshot& meshSource,
    const HermiteGalsState& state,const std::vector<int>& cells)
{
    const auto payload=materialSubsetPayload(mesh,state,cells);
    auto compact=Snapshot::fromNativeLeafSubset(meshSource,payload);
    const auto referenceView=state.geometry();const auto& reference=referenceView.snapshot();
    check(compact.nativeCellCount()==reference.nativeCellCount(),"compact subset changed global native index space");
    int active=0;
    for(int c=0;c<mesh.cellIndexNum();++c){const bool requested=std::binary_search(cells.begin(),cells.end(),c);
        check(compact.nativeCellActive(c)==requested,"compact subset fabricated unrequested active geometry");
        if(requested){check(compact.identicalNativeQuarticCell(c,reference,c),
            "compact subset changed original full native Q4 exact components");++active;}}
    check(active==static_cast<int>(cells.size())&&!payload.nativeVertices.empty(),"compact witness lacks actual support");
    auto missing=payload;missing.nativeVertices.erase(missing.nativeVertices.begin());missing.nodes.erase(missing.nodes.begin());
    rejects([&]{Snapshot::fromNativeLeafSubset(meshSource,missing);},"subset support");
    auto unordered=payload;unordered.nativeCells.front()=-1;
    rejects([&]{Snapshot::fromNativeLeafSubset(meshSource,unordered);},"subset native cell indices");
    auto missingEdge=payload;check(!missingEdge.edgeMidpoints.empty(),"compact witness lacks independent side");
    missingEdge.edgeMidpoints.erase(missingEdge.edgeMidpoints.begin());
    rejects([&]{Snapshot::fromNativeLeafSubset(meshSource,missingEdge);},"midpoint");
    const int first=cells.front();const auto p=mesh.cells(mesh.cellIndexToLabel[first]).centerCoord();
    compact.scalarJetInNativeCell(first,p); // warm source-owned scalar data
    {
        struct Restore {int saved=std::fegetround();~Restore(){std::fesetround(saved);}} restore;
        check(restore.saved==FE_TONEAREST&&std::fesetround(FE_DOWNWARD)==0,"cannot exercise compact actual FP guard");
        rejects([&]{Snapshot::fromNativeLeafSubset(meshSource,payload);},"round-to-nearest");
        rejects([&]{compact.scalarJetInNativeCell(first,p);},"round-to-nearest");
    }
#if defined(_MSC_VER) || defined(__SSE__)
    {
        struct Restore {unsigned saved=_mm_getcsr();~Restore(){_mm_setcsr(saved);}} restore;
        _mm_setcsr(restore.saved|(1u<<15)|(1u<<6));
        rejects([&]{Snapshot::fromNativeLeafSubset(meshSource,payload);},"gradual underflow");
        rejects([&]{compact.scalarJetInNativeCell(first,p);},"gradual underflow");
    }
#endif
    check(compact.identicalNativeQuarticCell(first,reference,first),"compact FP rejection mutated retained exact Q4");
    for(int c=0;c<mesh.cellIndexNum();++c)if(!compact.nativeCellActive(c)){
        rejects([&]{compact.scalarJetInNativeCell(c,mesh.cells(mesh.cellIndexToLabel[c]).centerCoord());},"inactive leaf");break;}
    return compact;
}
void materialAmrSubsetWitness()
{
    auto mesh=makeMesh(4,1);Case shape;shape.pureFirstJet=true;shape.amr=true;
    const auto state=stateFor(*mesh,shape.profile(0),false,{0,0,0,Position::Initial},
        [](int,Vector2D,Vector2D){return 0.;},{},HermiteGalsTransport::FirstJet);
    const auto view=state.geometry();const auto& full=view.snapshot();
    materialSubsetGeometryWitness(*mesh,full,state,state.primary().geometryCells);
    auto missing=materialSubsetPayload(*mesh,state,state.primary().geometryCells);int coarse=-1;
    for(int v=0;v<mesh->vertIndexNum();++v){const auto label=mesh->verts(mesh->vertIndexToLabel[v]).supCell;
        if(label>=0){coarse=mesh->cells(label).index();break;}}
    check(coarse>=0,"compact AMR witness has no actual hanging coarse owner");
    const auto pos=std::lower_bound(missing.nativeCells.begin(),missing.nativeCells.end(),coarse);
    check(pos!=missing.nativeCells.end()&&*pos==coarse,"compact AMR coarse owner absent in full source");
    const auto offset=pos-missing.nativeCells.begin();missing.nativeCells.erase(pos);
    missing.cellCenters.erase(missing.cellCenters.begin()+offset);
    rejects([&]{Snapshot::fromNativeLeafSubset(full,missing);},"subset support");
    std::cout<<"COMPACT_NATIVE_Q4,AMR_COARSE_OWNER,EXACT_FULL_COMPONENTS,FP,SUPPORT,PASS\n";
}
// Compare one REAL internal GL3 time against the retained full-publication
// path. The two independently completed candidates keep their own IDs.
void firstJetScalarSampleWitness(const QuadtreeMesh& mesh,const HermiteGalsPublication& original,
    const HydroCharacteristicSource& motion,const HydroSpaceTimeConfig& config,const Case& scene)
{
    using Projection=HermiteGalsTimeProjection;using Jacobian=LevelSetBoundaryTrace::Jacobian;
    const auto& reconstruction=motion.geometry().reconstruction();
    const auto parent=motion.arrival().timeLayer();
    const double begin=parent.begin,width=parent.end-begin,center=begin+width*.5;
    auto clock=parent;clock.end=center+width*.5*(-.77459666924148337704);
    const HermiteGalsResolvedMotion resolved(original.state().identity(),original.state().timeLayer(),
        motion.arrival().geometryIdentity(),clock,[motion,clock,config,begin](Vector2D point){
            const auto t=motion.traceBackward(point,clock.time(),begin,config.characteristic);
            std::optional<Boundary> side;
            if(t.boundary()){const auto b=*t.boundary();
                side=b.normalAxis==0?(b.outwardSign<0?Boundary::west:Boundary::east)
                    :(b.outwardSign<0?Boundary::south:Boundary::north);}
            return LevelSetBoundaryTrace::JetResult{{t.point(),clock.time()-t.time(),side,t.boundary().has_value(),false},
                t.jacobian(),t.substeps()};
        });
    const auto full=prepareHermiteGalsFirstJetStep(mesh,original,resolved,{}, {},
        HermiteGalsBoundaryTreatment::InflowHistory,Projection::ClosedProjectionCompensated);
    const auto sample=prepareHermiteGalsFirstJetSample(mesh,original,resolved,{}, {},
        HermiteGalsBoundaryTreatment::InflowHistory,Projection::ClosedProjectionCompensated);
    sample.validateMesh(mesh);
    const auto view=sample.geometry();const auto& scalar=view.snapshot();
    const auto& referenceScalar=full.publication().state().geometry().snapshot();
    const auto sameClock=[](HermiteGalsTimeLayer a,HermiteGalsTimeLayer b){
        return a.step==b.step && a.begin==b.begin && a.end==b.end && a.position==b.position;
    };
    check(sample.sourceIdentity()==original.state().identity() && full.sourceIdentity()==sample.sourceIdentity()
        && sample.arrivalVelocityIdentity()==resolved.velocityIdentity()
        && full.arrivalVelocityIdentity()==sample.arrivalVelocityIdentity()
        && view.identity()==sample.identity() && sameClock(view.timeLayer(),clock)
        && sameClock(sample.timeLayer(),clock) && sameClock(full.publication().state().timeLayer(),clock),
        "scalar sample detached original/resolved owner or full clock");
    const auto copied=sample;
    check(copied.identity()==sample.identity() && &copied.geometry().snapshot()==&scalar,
        "copied scalar sample duplicated unique Q4 owner");
    check(scalar.usesNativeLeaves() && scalar.reconstruction()==LevelSetHermiteGeometry::Reconstruction::EnrichedQuartic
        && scalar.nativeCellCount()==referenceScalar.nativeCellCount(),"scalar sample is not the original native Q4");
    int cells=0;
    for(int c=0;c<scalar.nativeCellCount();++c){
        check(scalar.nativeCellActive(c)==referenceScalar.nativeCellActive(c),"scalar sample changed active support");
        if(!scalar.nativeCellActive(c))continue;
        // This exact native constrained-polynomial certificate covers all
        // primary/closure coefficients, not just displayed or sampled jets.
        check(scalar.identicalNativeQuarticCell(c,referenceScalar,c),
            "scalar sample differs from full-publication authoritative Q4");
        ++cells;
    }
    check(cells>0,"scalar sample has no actual active Q4 leaves");
    const auto geometry=HydroPhaseGeometry::build(mesh,reconstruction.momentumGeometry().metrics(),full.publication());
    const auto fullAdvector=reconstruction.interpolate(mesh,motion.previous(),motion.arrival(),
        geometry.geometryIdentity(),geometry.timeLayer());
    const auto sampledAdvector=reconstruction.interpolate(mesh,motion.previous(),motion.arrival(),
        sample.identity(),sample.timeLayer());
    const auto old=HydroPhaseRate::prepareSpatial(mesh,geometry,fullAdvector);
    const auto next=HydroPhaseRate::sampleAtomRates(mesh,sample,sampledAdvector);
    check(next.size()==old.atomRates().size() && next.size()==reconstruction.atoms().size(),
        "scalar sample changed actual phase-atom inventory");
    sampledAdvector.validateFor(mesh,reconstruction.momentumGeometry(),sample.identity(),sample.timeLayer());
    long double peak=0;int cut=0,tiny=0,nonzero=0;
    for(std::size_t i=0;i<next.size();++i){
        const auto& polynomial=sampledAdvector.atomPolynomials()[i];long double work=0,power=1;
        for(const auto c:polynomial.beginCoefficients){work+=std::abs(c)*power;power*=polynomial.length;}
        work*=reconstruction.atoms()[i].area;
        const auto a=old.atomRates()[i],b=next[i];
        check(std::isfinite(b.liquid) && std::isfinite(b.gas),"nonfinite sampled phase rate");
        if(work==0)check(a.liquid==0 && a.gas==0 && b.liquid==0 && b.gas==0,"scalar sample broke zero-work atom");
        else{
            peak=std::max({peak,
                normalized(b.liquid,a.liquid,work,1e-12,"scalar sample signed liquid atom"),
                normalized(b.gas,a.gas,work,1e-12,"scalar sample signed gas atom"),
                normalized(b.liquid+b.gas,polynomial.integral(false,0,polynomial.length),work,1e-12,
                    "scalar sample original atom closure")});
        }
        for(const auto phase:{std::array<long double,2>{a.liquid,b.liquid},std::array<long double,2>{a.gas,b.gas}})
            if(phase[0]!=0){
                check(phase[1]!=0 && (phase[1]>0)==(phase[0]>0),"scalar sample erased or reversed resolved phase");
                if(std::abs(phase[0])<1e-10L*work){
                    peak=std::max(peak,normalized(phase[1],phase[0],std::abs(phase[0]),1e-12,
                        "scalar sample small phase actual rate scale"));++tiny;
                }
                ++nonzero;
            }
        if(a.liquid!=0 && a.gas!=0)++cut;
    }
    check(cut>0 && nonzero>0,"scalar witness did not exercise actual cut/nonzero phase rates");
    // The same actual internal time exercises atom-only endpoint interpolation.
    // Exact sampling is explicit here; Numerical has its separate point budget.
    const auto slab=HydroAdvectorTimeSlab::prepare(mesh,original,reconstruction,motion.previous(),motion.arrival());
    std::vector<int> requested;
    for(std::size_t i=0;i<next.size();++i)if(next[i].liquid!=0&&next[i].gas!=0){requested.push_back(static_cast<int>(i));break;}
    check(!requested.empty()&&next.size()>1,"local slab witness has no cut atom");
    requested.push_back((requested.front()+1)%static_cast<int>(next.size()));
    const auto required=slab.nativeCells(requested);
    check(required==HydroPhaseRate::requiredNativeCells(mesh,reconstruction,requested),"time slab changed closed native support");
    const auto preparedSupport=slab.prepareAtomSupport(requested);
    check(preparedSupport.atomIds()==requested&&preparedSupport.nativeCells()==required,
        "prepared atom support changed caller order or original closed native support");
    for(std::size_t id=0;id<reconstruction.atoms().size();++id) {
        const std::vector<int> one{static_cast<int>(id)};
        check(slab.prepareAtomSupport(one).nativeCells()==HydroPhaseRate::requiredNativeCells(mesh,reconstruction,one),
            "prepared incidence changed actual mesh closed sides/corners");
    }
    const auto independentSlab=HydroAdvectorTimeSlab::prepare(mesh,original,reconstruction,motion.previous(),motion.arrival());
    rejects([&]{independentSlab.validateAtomSupport(preparedSupport);},"different prepared atom support owner");
    rejects([&]{slab.prepareAtomSupport({requested.front(),requested.front()});},"duplicate atom");
    rejects([&]{slab.prepareAtomSupport({-1});},"invalid atom");
    check(slab.sourceIdentity()==original.state().identity()
        &&slab.arrivalGeometryIdentity()==motion.arrival().geometryIdentity()
        &&slab.previousIdentity()==motion.previous().identity()&&slab.arrivalIdentity()==motion.arrival().identity()
        &&sameClock(slab.endpointTimeLayer(),parent),"time slab detached endpoint source/full clock");
    struct PointThreads {
#ifdef _OPENMP
        int count=omp_get_max_threads(),dynamic=omp_get_dynamic(),maxActiveLevels=omp_get_max_active_levels();
        explicit PointThreads(int n){omp_set_dynamic(0);omp_set_max_active_levels(std::min(maxActiveLevels,1));omp_set_num_threads(n);}
        ~PointThreads(){omp_set_num_threads(count);omp_set_dynamic(dynamic);omp_set_max_active_levels(maxActiveLevels);}
#else
        explicit PointThreads(int){}
#endif
    };
    const auto localAt=[&](int threads){
        PointThreads team(threads);
        return prepareHermiteGalsFirstJetLocalSample(mesh,original,resolved,required,{}, {},
            HermiteGalsBoundaryTreatment::InflowHistory,Projection::ClosedProjectionCompensated,
            HermiteGalsPointExecution::Parallel,HermiteGalsPointSamplingConfig{});
    };
    const auto localOne=localAt(1),localFour=localAt(4);
    // One immutable transaction is reused across real times/point teams;
    // the old free-local path remains the independent arithmetic oracle.
    const auto prepared=HermiteGalsFirstJetLocalPreparation::prepare(mesh,original);
    check(prepared.sourceIdentity()==original.state().identity(),"local preparation lost actual source owner");
    prepared.validateMesh(mesh);
    const auto preparedSerial=prepared.sample(resolved,required,{},Projection::ClosedProjectionCompensated);
    const auto preparedAt=[&](int threads){PointThreads team(threads);
        return prepared.sample(resolved,required,{},Projection::ClosedProjectionCompensated,
            HermiteGalsPointExecution::Parallel,HermiteGalsPointSamplingConfig{});};
    const auto preparedFour=preparedAt(4);
    const auto preparedMatches=[&](const HermiteGalsFirstJetLocalSample& candidate,
        const HermiteGalsFirstJetLocalSample& reference){
        check(candidate.sourceIdentity()==reference.sourceIdentity()
            &&candidate.arrivalVelocityIdentity()==reference.arrivalVelocityIdentity()
            &&sameClock(candidate.timeLayer(),reference.timeLayer())
            &&candidate.nativeCells()==reference.nativeCells()
            &&candidate.measure()==reference.measure(),"prepared local sample changed provenance/full clock/native support");
        const auto actualView=candidate.geometry(),referenceView=reference.geometry();
        const auto& actual=actualView.snapshot();const auto& expected=referenceView.snapshot();
        check(actual.nativeCellCount()==expected.nativeCellCount(),"prepared local sample changed native index space");
        for(int c=0;c<mesh.cellIndexNum();++c){check(actual.nativeCellActive(c)==expected.nativeCellActive(c),
            "prepared local sample fabricated or erased native support");
            if(actual.nativeCellActive(c))check(actual.identicalNativeQuarticCell(c,expected,c),
                "prepared local sample changed original exact primary/closure Q4 components");}
    };
    preparedMatches(preparedSerial,localOne);preparedMatches(preparedFour,localFour);
    check(preparedFour.report().pointThreads==localFour.report().pointThreads
        &&preparedFour.report().pointJobs==localFour.report().pointJobs,
        "prepared warm sample did not execute the actual parallel point team");
    const auto oldPullback=prepareHermiteGalsFirstJetLocalSample(mesh,original,resolved,required,{}, {},
        HermiteGalsBoundaryTreatment::InflowHistory,Projection::CharacteristicPullback);
    const auto preparedPullback=prepared.sample(resolved,required);
    preparedMatches(preparedPullback,oldPullback);
    // A closed identity map makes the two original-source departure points
    // equal mathematically. The staged snapshot remains a different owner.
    // Count motion calls independently, and compare to the original Q4 and
    // this existing circle's analytic first jet, not another memoized producer.
    const auto sourceGeometry=original.state().geometry();const auto& sourceScalar=sourceGeometry.snapshot();
    const auto profile=scene.profile(original.state().timeLayer().time());
    const double identityDuration=clock.time()-original.state().timeLayer().time();
    const auto memoSample=[&](bool numerical,bool parallel){
        PointThreads team(parallel?4:1);
        std::map<std::pair<double,double>,std::size_t> visits;std::mutex visitsMutex;
        const HermiteGalsResolvedMotion identityMap(original.state().identity(),original.state().timeLayer(),
            resolved.velocityIdentity(),clock,[&](Vector2D q){
                {std::lock_guard<std::mutex> lock(visitsMutex);++visits[{q.x,q.y}];}
                return LevelSetBoundaryTrace::JetResult{{q,identityDuration,{},false,false},Jacobian::identity(),1};
            });
        HermiteGalsPointSamplingConfig sampling;
        if(numerical){sampling.mode=HermiteGalsPointSampling::Numerical;
            sampling.valueAbsoluteTolerance=1e-13;sampling.gradientAbsoluteTolerance=1e-11;}
        const auto result=prepared.sample(identityMap,required,{},Projection::ClosedProjectionCompensated,
            parallel?HermiteGalsPointExecution::Parallel:HermiteGalsPointExecution::Serial,sampling);
        const auto& report=result.report();std::size_t traces=0,second=0;
        for(const auto& entry:visits){check(entry.second==1||entry.second==2,
            "identity memo witness changed the retained first/second trace sequence");
            traces+=entry.second;second+=entry.second==2;}
        check(second>0&&report.stagedNativeCells>0&&report.pointTraceQueries==traces,
            "source jet memo erased or falsely reported an actual second trace");
        // Each staged first-map point has one original scalar miss. Each
        // corrected point still queries its distinct projected snapshot.
        // The old uncached one/twice pair would add two scalar calls per point.
        const std::size_t originalUncached=visits.size()+3*second;
        check(report.evaluatedPointQueries==visits.size()+second
            &&report.evaluatedPointQueries<originalUncached,
            "source jet memo did not reduce real scalar API calls at equal departures");
        check(report.pointSampling==sampling.mode
            &&report.numericalPointQueries==(numerical?report.evaluatedPointQueries:0)
            &&report.numericalValueFallbacks<=report.numericalPointQueries
            &&report.numericalGradientFallbacks<=report.numericalPointQueries,
            "source jet memo mixed exact/numerical sampling or fabricated miss/fallback work");
        if(numerical)check(std::isfinite(report.maxAcceptedValueErrorBound)
            &&std::isfinite(report.maxAcceptedGradientErrorBound)
            &&report.maxAcceptedValueErrorBound<=sampling.valueAbsoluteTolerance
            &&report.maxAcceptedGradientErrorBound<=sampling.gradientAbsoluteTolerance,
            "memoized numerical source jet lost its actual physical error budget");
        check(result.sourceIdentity()==original.state().identity()
            &&result.arrivalVelocityIdentity()==resolved.velocityIdentity()
            &&sameClock(result.timeLayer(),clock)&&result.nativeCells()==required,
            "source jet memo detached original/arrival/full clock/native support");
        const auto resultGeometry=result.geometry();const auto& actual=resultGeometry.snapshot();
        for(int cell:required){const auto box=sourceScalar.nativeCellBounds(cell);
            for(double x:{0.,.5,1.})for(double y:{0.,.5,1.}){
                const Vector2D q{box[0].x+x*(box[1].x-box[0].x),box[0].y+y*(box[1].y-box[0].y)};
                const auto expected=profile(q);const auto originalJet=sourceScalar.scalarJetInNativeCell(cell,q);
                // These original 1e-12 direct FirstJet gates also establish
                // that the frozen native source is this polynomial oracle.
                close(originalJet.value,expected.value,1e-12,"memo witness original Q4 value oracle");
                close(originalJet.dx,expected.dx,1e-12,"memo witness original Q4 x-gradient oracle");
                close(originalJet.dy,expected.dy,1e-12,"memo witness original Q4 y-gradient oracle");
                const auto got=actual.scalarJetInNativeCell(cell,q);
                close(got.value,expected.value,1e-12,"memo identity compensated Q4 value");
                close(got.dx,expected.dx,1e-12,"memo identity compensated Q4 x gradient");
                close(got.dy,expected.dy,1e-12,"memo identity compensated Q4 y gradient");}}
        std::cout<<"SOURCE_JET_MEMO,"<<(numerical?"NUMERICAL":"EXACT")<<','
            <<(parallel?"POINT_4":"SERIAL")<<','<<visits.size()<<','<<second<<','
            <<report.evaluatedPointQueries<<','<<originalUncached<<','<<report.pointTraceQueries
            <<",INDEPENDENT_Q4_PROFILE_SECOND_TRACE,PASS\n";
        return result;
    };
    for(bool numerical:{false,true}){
        const auto memoSerial=memoSample(numerical,false),memoFour=memoSample(numerical,true);
        preparedMatches(memoFour,memoSerial);
        check(memoSerial.report().evaluatedPointQueries==memoFour.report().evaluatedPointQueries
            &&memoSerial.report().pointTraceQueries==memoFour.report().pointTraceQueries
            &&memoSerial.report().numericalPointQueries==memoFour.report().numericalPointQueries
            &&memoSerial.report().numericalValueFallbacks==memoFour.report().numericalValueFallbacks
            &&memoSerial.report().numericalGradientFallbacks==memoFour.report().numericalGradientFallbacks
            &&memoSerial.report().maxAcceptedValueErrorBound==memoFour.report().maxAcceptedValueErrorBound
            &&memoSerial.report().maxAcceptedGradientErrorBound==memoFour.report().maxAcceptedGradientErrorBound,
            "serial/parallel memo changed actual miss, fallback or certified bound data");}
    int rejectedQueries=0;
    const auto shouldNotTrace=[&](Vector2D q){++rejectedQueries;return resolved.trace(q);};
    const HermiteGalsResolvedMotion wrongSource(full.publication().state().identity(),resolved.sourceTime(),
        resolved.velocityIdentity(),resolved.targetTime(),shouldNotTrace);
    rejects([&]{prepared.sample(wrongSource,required);},"resolved motion source identity/clock differs");
    auto detachedClock=resolved.sourceTime();++detachedClock.step;
    const HermiteGalsResolvedMotion wrongSourceClock(resolved.sourceIdentity(),detachedClock,
        resolved.velocityIdentity(),resolved.targetTime(),shouldNotTrace);
    rejects([&]{prepared.sample(wrongSourceClock,required);},"resolved motion source identity/clock differs");
    check(rejectedQueries==0,"prepared source/full-clock refusal sampled a detached motion");
    auto duplicateCells=required;duplicateCells.push_back(required.front());
    rejects([&]{prepared.sample(resolved,duplicateCells);},"duplicate requested native cell");
    {
        struct Restore {int saved=std::fegetround();~Restore(){std::fesetround(saved);}} restore;
        check(restore.saved==FE_TONEAREST&&std::fesetround(FE_DOWNWARD)==0,"cannot exercise warm transaction FP guard");
        rejects([&]{prepared.sample(resolved,required);},"round-to-nearest");
    }
#if defined(_MSC_VER) || defined(__SSE__)
    {
        struct Restore {unsigned saved=_mm_getcsr();~Restore(){_mm_setcsr(saved);}} restore;
        _mm_setcsr(restore.saved|(1u<<15)|(1u<<6));
        rejects([&]{prepared.sample(resolved,required);},"gradual underflow");
    }
#endif
    // Core owns the immutable mesh index space. Live-mesh validation remains
    // explicit, while sample() must not read this subsequently changed mesh.
    auto changedMesh=makeMesh(4);
    const auto frozenPrepared=HermiteGalsFirstJetLocalPreparation::prepare(*changedMesh,original);
    auto& changedVertex=changedMesh->verts(changedMesh->vertIndexToLabel[0]);
    changedVertex.setCoord(changedVertex.coord()+Vector2D{.01,0});
    rejects([&]{frozenPrepared.validateMesh(*changedMesh);},"mesh");
    const auto frozenSample=frozenPrepared.sample(resolved,required,{},Projection::ClosedProjectionCompensated);
    preparedMatches(frozenSample,localOne);
    std::cout<<"LOCAL_PREPARED_TRANSACTION,EXACT_PULLBACK_COMPENSATED,WARM_POINT_TEAM_4,FROZEN_MESH,SOURCE_CLOCK_FP,PASS\n";
    const auto oneView=localOne.geometry(),fourView=localFour.geometry();
    const auto& oneScalar=oneView.snapshot();const auto& fourScalar=fourView.snapshot();
    const auto meshSourceView=original.state().geometry();
    const auto compact=materialSubsetGeometryWitness(mesh,meshSourceView.snapshot(),full.publication().state(),required);
    check(required.size()<static_cast<std::size_t>(mesh.cellIndexNum()),"compact local witness selected the full mesh");
    for(int c:required)check(oneScalar.identicalNativeQuarticCell(c,compact,c)
        &&fourScalar.identicalNativeQuarticCell(c,compact,c),"local point assembly differs from independent compact full oracle");
    auto regridded=makeMesh(4,1);
    rejects([&]{localOne.validateMesh(*regridded);},"mesh");
    rejects([&]{prepareHermiteGalsFirstJetLocalSample(*regridded,original,resolved,required,{}, {},
        HermiteGalsBoundaryTreatment::InflowHistory,Projection::ClosedProjectionCompensated);},"mesh");
    int actualFour=1;
#ifdef _OPENMP
    actualFour=4;
#endif
    check(localOne.report().pointExecution==HermiteGalsPointExecution::Parallel
        &&localFour.report().pointExecution==HermiteGalsPointExecution::Parallel
        &&localOne.report().pointThreads==1&&localFour.report().pointThreads==actualFour
        &&localOne.report().pointJobs>0&&localOne.report().pointJobs==localFour.report().pointJobs
        &&localOne.report().pointTraceQueries==localFour.report().pointTraceQueries
        &&localOne.report().evaluatedPointQueries==localFour.report().evaluatedPointQueries,
        "local point mode did not execute the actual 1/4 teams");
    check(localOne.report().pointSampling==HermiteGalsPointSampling::Exact
        &&localFour.report().pointSampling==HermiteGalsPointSampling::Exact
        &&localOne.report().numericalPointQueries==0&&localFour.report().numericalPointQueries==0,
        "exact point comparison used the separate numerical sampling route");
    for(int c:required)check(oneScalar.identicalNativeQuarticCell(c,fourScalar,c)
        &&oneScalar.identicalNativeQuarticCell(c,scalar,c),"point teams changed authoritative requested Q4");
    for(const auto* local:{&localOne,&localFour})
        check(local->sourceIdentity()==original.state().identity()&&local->arrivalVelocityIdentity()==resolved.velocityIdentity()
            &&sameClock(local->timeLayer(),clock),"point team lost source/arrival/full clock");
    const auto denseLocal=reconstruction.interpolate(mesh,motion.previous(),motion.arrival(),localOne.identity(),clock);
    const auto denseRates=HydroPhaseRate::sampleAtomRates(mesh,localOne,denseLocal,requested);
    const auto ratesAt=[&](const HermiteGalsFirstJetLocalSample& local,int threads){
        PointThreads team(threads);return HydroPhaseRate::sampleAtomRates(local,slab,requested);
    };
    const auto slabOne=ratesAt(localOne,1),slabFour=ratesAt(localFour,4);
    const auto preparedRates=ratesAt(preparedFour,4);
    const auto sameGroupRates=HydroPhaseRate::sampleAtomRates(preparedFour,slab,preparedSupport);
    check(sameGroupRates.size()==preparedRates.size(),"prepared rate group changed ordered inventory");
    for(std::size_t k=0;k<sameGroupRates.size();++k)check(sameGroupRates[k].liquid==preparedRates[k].liquid
        &&sameGroupRates[k].gas==preparedRates[k].gas,"prepared group changed original phase flux arithmetic");
    rejects([&]{HydroPhaseRate::sampleAtomRates(preparedFour,independentSlab,preparedSupport);},"different prepared atom support owner");
    for(std::size_t q=0;q<requested.size();++q){
        const int id=requested[q];const auto a=slab.atomPolynomial(id,clock),b=denseLocal.regularAtomPolynomials()[id];
        check(a.beginCoefficients==b.beginCoefficients&&a.endCoefficients==b.endCoefficients
            &&a.length==b.length&&a.radiusBegin==b.radiusBegin&&a.radiusEnd==b.radiusEnd
            &&a.radialMeasure==b.radialMeasure&&a.radialTangent==b.radialTangent,
            "atom-only time slab changed original polynomial arithmetic");
        check(denseRates[q].liquid==slabOne[q].liquid&&denseRates[q].gas==slabOne[q].gas
            &&slabOne[q].liquid==slabFour[q].liquid&&slabOne[q].gas==slabFour[q].gas,
            "atom-only time slab or point team changed signed phase rates");
        check(preparedRates[q].liquid==slabFour[q].liquid&&preparedRates[q].gas==slabFour[q].gas,
            "prepared transaction changed ordered signed phase atom rates");
    }
    rejects([&]{slab.validateSample(localOne.identity(),{},localOne.arrivalVelocityIdentity(),clock);},"different original publication");
    auto wrongClock=clock;wrongClock.end=parent.end+width;
    rejects([&]{slab.validateSample(localOne.identity(),localOne.sourceIdentity(),localOne.arrivalVelocityIdentity(),wrongClock);},
        "dense geometry clock outside endpoint slab");
    std::cout<<"LOCAL_PHASE_TIME_SLAB,"<<requested.size()<<','<<required.size()<<",EXACT_POINT_TEAMS,"
        <<localOne.report().pointThreads<<','<<localFour.report().pointThreads<<",ATOM_RATE_TEAMS_1_4,SIGNED_RATES_BITWISE,PASS\n";
    // The same real non-affine endpoint motion and compensated source now
    // supply every actual primal/dual/internal atom through a trace packet.
    // The complete unchanged sample above is the independent Q4 oracle.
    using TraceSample=HermiteGalsFirstJetTraceSample;
    std::vector<int> traceIds;std::vector<TraceSample::Segment> segments;
    for(std::size_t i=0;i<reconstruction.atoms().size();++i){
        traceIds.push_back(static_cast<int>(i));const auto& atom=reconstruction.atoms()[i];
        segments.push_back({atom.begin,atom.end,atom.normalAxis});
    }
    const auto traceCells=slab.nativeCells(traceIds);
    const auto tracesAt=[&](int threads,const std::vector<TraceSample::Segment>& actual){
        PointThreads team(threads);
        return prepared.sampleTraces(resolved,traceCells,actual,{},Projection::ClosedProjectionCompensated,
            HermiteGalsPointExecution::Parallel,HermiteGalsPointSamplingConfig{});
    };
    const auto traceOne=tracesAt(1,segments),traceFour=tracesAt(4,segments);
    check(traceOne.sourceIdentity()==original.state().identity()&&traceFour.sourceIdentity()==original.state().identity()
        &&traceOne.arrivalVelocityIdentity()==resolved.velocityIdentity()
        &&sameClock(traceOne.sourceTimeLayer(),original.state().timeLayer())
        &&sameClock(traceOne.timeLayer(),clock)&&sameClock(traceFour.timeLayer(),clock)
        &&traceOne.nativeCells()==traceCells&&traceOne.nativeCellCount()==scalar.nativeCellCount(),
        "trace sample detached real source/full clocks/native support");
    std::vector<HydroPhaseMeasure> traceRatesOne,traceRatesFour;
    {PointThreads team(1);traceRatesOne=HydroPhaseRate::sampleAtomRates(traceOne,slab,traceIds);}
    {PointThreads team(4);traceRatesFour=HydroPhaseRate::sampleAtomRates(traceFour,slab,traceIds);}
    const auto traceSupport=slab.prepareAtomSupport(traceIds);
    const auto preparedTraceRates=HydroPhaseRate::sampleAtomRates(traceFour,slab,traceSupport);
    check(preparedTraceRates.size()==traceRatesFour.size(),"prepared trace group changed full ordered inventory");
    for(std::size_t k=0;k<preparedTraceRates.size();++k)check(preparedTraceRates[k].liquid==traceRatesFour[k].liquid
        &&preparedTraceRates[k].gas==traceRatesFour[k].gas,"prepared trace group changed original phase flux arithmetic");
    rejects([&]{HydroPhaseRate::sampleAtomRates(traceFour,independentSlab,traceSupport);},"different prepared atom support owner");
    rejects([&]{HydroPhaseRate::sampleAtomRates(localOne,slab,traceSupport);},"missing local native support");
    check(traceRatesOne.size()==next.size()&&traceRatesFour.size()==next.size(),"trace packet changed actual atom inventory");
    int cardinalTraces=0;long double tracePeak=0;
    for(std::size_t q=0;q<traceIds.size();++q){
        const auto& atom=reconstruction.atoms()[q];const auto& polynomial=sampledAdvector.regularAtomPolynomials()[q];
        long double work=0,power=1;for(auto c:polynomial.beginCoefficients){work+=std::abs(c)*power;power*=polynomial.length;}
        work*=atom.area;const auto a=next[q],b=traceRatesOne[q],c=traceRatesFour[q];
        if(work==0)check(a.liquid==0&&a.gas==0&&b.liquid==0&&b.gas==0,"trace packet broke zero-work atom");
        else for(const auto values:{std::array<long double,2>{a.liquid,b.liquid},std::array<long double,2>{a.gas,b.gas}}){
            tracePeak=std::max(tracePeak,normalized(values[1],values[0],work,1e-12,"trace/full actual signed atom rate"));
            if(values[0]!=0&&std::abs(values[0])<1e-10L*work)
                normalized(values[1],values[0],std::abs(values[0]),1e-12,"trace/full resolved tiny phase rate");
        }
        check(b.liquid==c.liquid&&b.gas==c.gas,"trace point/rate teams changed ordered signed phase rates");
        auto actual=traceOne.trace(q);check(traceOne.ownsTrace(actual)&&!traceFour.ownsTrace(actual),"trace owner authentication changed");
        if(std::holds_alternative<LevelSetMappedSegment::CardinalQuarticTrace>(actual))++cardinalTraces;
        std::visit([&](auto& trace){if(!trace.events().empty())trace=trace.refined(0);},actual);
        check(traceOne.ownsTrace(actual),"trace refinement detached original packet owner");
    }
    check(cardinalTraces>0&&traceOne.report().pointThreads==1&&traceFour.report().pointThreads==actualFour,
        "real trace packet did not exercise cardinal data and actual point teams");
    const auto fullTraceOwner=std::make_shared<const HermiteGalsGeometryView>(view);
    const std::shared_ptr<const Snapshot> fullTraceScalar(fullTraceOwner,&fullTraceOwner->snapshot());
    const auto compareTrace=[&](const TraceSample& packet,std::size_t q,const TraceSample::Segment& segment){
        auto actual=packet.trace(q);
        auto reference=LevelSetMappedSegment::nativeQuarticTrace(fullTraceScalar,segment.begin,segment.end,segment.normalAxis);
        std::visit([&](const auto& trace){
            std::array<long double,2> a{},b{};
            for(const auto& r:trace.ordinaryRanges())a[r.negative?0:1]+=r.length;
            for(const auto& r:reference.ordinaryRanges())b[r.negative?0:1]+=r.length;
            for(int phase=0;phase<2;++phase)close(a[phase],b[phase],1e-12,"trace/full actual subsegment phase length");
            check(trace.coincidentRanges().size()==reference.coincidentRanges().size(),"trace packet changed zero-trace inventory");
        },actual);
    };
    const auto& cutAtom=reconstruction.atoms()[requested.front()];const auto delta=cutAtom.end-cutAtom.begin;
    const std::vector<TraceSample::Segment> subsegments{
        {cutAtom.begin+delta*.25,cutAtom.begin+delta*.75,cutAtom.normalAxis},
        {cutAtom.begin+delta*.75,cutAtom.begin+delta*.25,cutAtom.normalAxis}};
    const auto subpacket=tracesAt(1,subsegments);
    for(std::size_t q=0;q<subsegments.size();++q)compareTrace(subpacket,q,subsegments[q]);
    auto noncardinal=subsegments.front();
    if(noncardinal.normalAxis==0){noncardinal.begin.x+=.03125;noncardinal.end.x+=.03125;}
    else{noncardinal.begin.y+=.03125;noncardinal.end.y+=.03125;}
    const auto fallback=tracesAt(1,std::vector<TraceSample::Segment>{noncardinal});
    check(std::holds_alternative<LevelSetMappedSegment::NativeQuarticTrace>(fallback.trace(0)),
        "noncardinal segment did not retain real full-Q4 fallback");
    compareTrace(fallback,0,noncardinal);
    auto wrongSegments=segments;std::swap(wrongSegments[0],wrongSegments[1]);
    const auto misplaced=tracesAt(1,wrongSegments);
    rejects([&]{HydroPhaseRate::sampleAtomRates(misplaced,slab,traceIds);},"different requested physical atom trace");
    rejects([&]{prepared.sampleTraces(wrongSource,traceCells,segments);},"resolved motion source identity/clock differs");
    rejects([&]{prepared.sampleTraces(wrongSourceClock,traceCells,segments);},"resolved motion source identity/clock differs");
    {
        struct Restore {int saved=std::fegetround();~Restore(){std::fesetround(saved);}} restore;
        check(std::fesetround(FE_DOWNWARD)==0,"cannot exercise trace FP guard");
        rejects([&]{HydroPhaseRate::sampleAtomRates(traceOne,slab,traceIds);},"round-to-nearest");
        rejects([&]{prepared.sampleTraces(resolved,traceCells,segments);},"round-to-nearest");
    }
#if defined(_MSC_VER) || defined(__SSE__)
    {
        struct Restore {unsigned saved=_mm_getcsr();~Restore(){_mm_setcsr(saved);}} restore;
        _mm_setcsr(restore.saved|(1u<<15)|(1u<<6));
        rejects([&]{HydroPhaseRate::sampleAtomRates(traceOne,slab,traceIds);},"gradual underflow");
    }
#endif
    std::cout<<"LOCAL_PHASE_TRACE_PACKET,"<<traceIds.size()<<','<<cardinalTraces<<','<<tracePeak
        <<",REAL_COMPENSATED_FULL_ORACLE,PRIMAL_DUAL,POINT_RATE_TEAMS_1_4,SUBSEGMENT_REVERSE,FULL_FALLBACK,SOURCE_CLOCK_FP,PASS\n";
    rejects([&]{HydroPhaseRate::sampleAtomRates(mesh,sample,fullAdvector);},"geometry identity");
    // Reuse the existing n4 FirstJet second-query refusal on this SAME source.
    const double duration=clock.time()-original.state().timeLayer().time();
    const HermiteGalsResolvedMotion secondHit(original.state().identity(),original.state().timeLayer(),
        resolved.velocityIdentity(),clock,[duration](Vector2D p){
            if(p.x>0 && p.x<1 && std::abs(8*p.x-std::round(8*p.x))>1e-10)
                return LevelSetBoundaryTrace::JetResult{{{0,p.y},duration*.5,Boundary::west,true,false},{},1};
            return LevelSetBoundaryTrace::JetResult{{{p.x==0 || p.x==1?p.x:p.x-.01,p.y},duration,{},false,false},
                Jacobian::identity(),1};
        });
    const auto refusal=[&](bool temporary){
        if(temporary)prepareHermiteGalsFirstJetSample(mesh,original,secondHit,{}, {},
            HermiteGalsBoundaryTreatment::InflowHistory,Projection::ClosedProjectionCompensated);
        else prepareHermiteGalsFirstJetStep(mesh,original,secondHit,{}, {},
            HermiteGalsBoundaryTreatment::InflowHistory,Projection::ClosedProjectionCompensated);
    };
    for(bool temporary:{false,true})rejects([&]{refusal(temporary);},
        "closed projection second query requires interior Jacobian");
    {
        struct RestoreRound {
            int saved=std::fegetround();
            ~RestoreRound(){std::fesetround(saved);}
        } restore;
        check(restore.saved==FE_TONEAREST,"warm trajectory fixture did not start in round-to-nearest");
        const auto warmPoint=point(mesh,original.state().primary().vertices.front().vertex);
        int visits=0;bool poisoned=false;std::string failure;
        // First preparation warms the existing firstQuery(q) trajectory and
        // evaluates the original source Q4. The corrected firstQuery is cached;
        // its second trajectory is followed by the original scalar API FP check.
        const HermiteGalsResolvedMotion warm(original.state().identity(),original.state().timeLayer(),
            resolved.velocityIdentity(),clock,[&,warmPoint,duration](Vector2D p){
                if(p.x==warmPoint.x && p.y==warmPoint.y && ++visits==2){
                    check(std::fesetround(FE_DOWNWARD)==0,"cannot exercise corrected scalar FP guard");poisoned=true;
                }
                return LevelSetBoundaryTrace::JetResult{{p,duration,{},false,false},Jacobian::identity(),1};
            });
        try{prepareHermiteGalsFirstJetSample(mesh,original,warm,{}, {},
            HermiteGalsBoundaryTreatment::InflowHistory,Projection::ClosedProjectionCompensated);}
        catch(const std::exception& error){failure=error.what();}
        check(std::fesetround(FE_TONEAREST)==0,"cannot restore round-to-nearest after corrected scalar refusal");
        check(poisoned && visits==2,"warm trajectory fixture missed the corrected second trajectory");
        check(failure.find("stage=closed projection vertex correction")!=std::string::npos
            && failure.find("round-to-nearest")!=std::string::npos,"corrected scalar query ignored actual FP environment");
        // The same callback poisons the corrected second trajectory after
        // the local prepare has cached the original identity-departure jet.
        visits=0;poisoned=false;failure.clear();
        try{prepareHermiteGalsFirstJetLocalSample(mesh,original,warm,original.state().primary().geometryCells,{}, {},
            HermiteGalsBoundaryTreatment::InflowHistory,Projection::ClosedProjectionCompensated);}
        catch(const std::exception& error){failure=error.what();}
        check(std::fesetround(FE_TONEAREST)==0,"cannot restore round-to-nearest after local cached scalar refusal");
        check(poisoned&&visits==2,"local cached scalar witness missed the corrected second trajectory");
        check(failure.find("stage=local mixed/departure closure")!=std::string::npos
            &&failure.find("round-to-nearest")!=std::string::npos,"local cached scalar query ignored actual FP environment");
    }
    std::cout<<"FIRST_JET_SCALAR_SAMPLE,ACTUAL_GL3,Q4_EXACT,"<<cells<<','<<next.size()<<','<<cut<<','<<tiny<<','
        <<peak<<",OWNER,SECOND_QUERY,WARM_FP,PASS\n";
}
// Compare sorted support discovery with the unchanged closed-box predicate.
void nativeSupportWitness(const QuadtreeMesh& mesh,const HydroAdvectorReconstruction& reconstruction)
{
    // Independent original closed-box predicate, without sorted-grid lookup.
    const auto brute=[&](const std::vector<int>& ids){
        std::vector<unsigned char> hit(mesh.cellIndexNum(),0);
        for(int id:ids){const auto& a=reconstruction.atoms().at(id);
            const Vector2D lo{std::min(a.begin.x,a.end.x),std::min(a.begin.y,a.end.y)};
            const Vector2D hi{std::max(a.begin.x,a.end.x),std::max(a.begin.y,a.end.y)};
            for(int c=0;c<mesh.cellIndexNum();++c){const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);
                const auto p=mesh.verts(cell.cellToVert[0]).coord(),q=mesh.verts(cell.cellToVert[2]).coord();
                if(hi.x>=std::min(p.x,q.x)&&lo.x<=std::max(p.x,q.x)
                    &&hi.y>=std::min(p.y,q.y)&&lo.y<=std::max(p.y,q.y))hit[c]=1;}
        }
        std::vector<int> result;for(int c=0;c<mesh.cellIndexNum();++c)if(hit[c])result.push_back(c);
        return result;
    };
    std::vector<int> all;
    for(std::size_t k=0;k<reconstruction.atoms().size();++k){const int id=static_cast<int>(k);all.push_back(id);
        const std::vector<int> one{id};
        check(HydroPhaseRate::requiredNativeCells(mesh,reconstruction,one)==brute(one),
            "sorted native support changed original closed atom predicate");}
    check(HydroPhaseRate::requiredNativeCells(mesh,reconstruction,all)==brute(all),
        "sorted native support changed original ordered union");
    check(!all.empty(),"support witness has no actual atoms");
    std::cout<<"NATIVE_SUPPORT,"<<all.size()<<",BRUTE_CLOSED_SUPPORT,PASS\n";
}
// Independent GL3 of original manufactured normals, including true RZ
// radius, checks the consumed whole-phase total for both velocity signs.
void certifiedSinglePhaseRates()
{
    int scenes=0;std::size_t rows=0,reverseCharts=0;long double peak=0;
    for(bool rz:{false,true}) {
        auto mesh=makeMesh(4);const auto metrics=metricsFor(*mesh,rz);
        const auto geometry=HydroMomentumGeometry::build(*mesh,metrics);
        const auto reconstruction=HydroAdvectorReconstruction::prepare(*mesh,geometry);
        const auto& layout=metrics.layout();
        for(bool liquid:{false,true}) {
            Case shape;shape.pureFirstJet=true;shape.planar=true;shape.radial=rz;
            shape.planeOffset=liquid?2:-1;
            const auto phase=HydroPhaseGeometry::build(*mesh,metrics,publishAt(*mesh,shape,{1,0,.04,Position::Midpoint}));
            auto view=std::make_shared<const HermiteGalsGeometryView>(phase.publication().state().geometry());
            const std::shared_ptr<const Snapshot> scalar(view,&view->snapshot());
            for(double sign:{1.,-1.}) {
                auto normal=mmsValues(*mesh,layout,2).normalVelocity();
                for(auto& value:normal)value*=sign;
                const auto actual=reconstruction.evaluate(*mesh,HydroNativeValues(layout,normal,std::vector<double>(layout.pressureCount())),
                    phase.geometryIdentity(),phase.timeLayer());
                const auto rates=HydroPhaseRate::prepareSpatial(*mesh,phase,actual);
                const auto flux=HydroMaterialFlux::prepare(*mesh,rates,1000,1);
                for(std::size_t i=0;i<reconstruction.atoms().size();++i) {
                    const auto& atom=reconstruction.atoms()[i];
                    const auto trace=LevelSetMappedSegment::nativeQuarticTrace(scalar,atom.begin,atom.end,atom.normalAxis);
                    const auto reverse=LevelSetMappedSegment::nativeQuarticTrace(scalar,atom.end,atom.begin,atom.normalAxis);
                    for(const auto* chart:{&trace,&reverse}) {
                        const auto& r=chart->ordinaryRanges();
                        check(chart->usesSource(*scalar)&&r.size()==1&&r[0].offset==0&&r[0].length==1
                            &&r[0].negative==liquid&&chart->events().empty()&&chart->coincidentRanges().empty(),
                            "whole-phase atom lacks original full-range certificate");
                    }
                    reverseCharts+=reverse.ordinaryRanges()[0].fromEnd;
                    auto coefficients=mmsCoefficients(atom,rz,2);for(auto& value:coefficients)value*=sign;
                    const auto ref=referenceAtom(atom,coefficients,shape);
                    const auto& got=rates.atomRates()[i];
                    const auto total=actual.atomPolynomials()[i].integral(false,0,actual.atomPolynomials()[i].length);
                    check(liquid?(got.liquid==total&&got.gas==0):(got.gas==total&&got.liquid==0),
                        "whole-phase atom changed original signed physical total");
                    peak=std::max(peak,normalized(got.liquid,ref.liquid,atom.area,1e-12,"pure liquid independent GL3"));
                    peak=std::max(peak,normalized(got.gas,ref.gas,atom.area,1e-12,"pure gas independent GL3"));
                    peak=std::max(peak,normalized(flux.atomRates()[i],1000*ref.liquid+ref.gas,1001*atom.area,1e-12,
                        "pure true-density material independent GL3"));
                    const auto& polynomial=actual.atomPolynomials()[i];
                    const auto beginRef=referenceShiftedIntegral(polynomial.beginCoefficients,0,polynomial.length,rz,
                        polynomial.radiusBegin,polynomial.radialTangent?1:0,1);
                    const auto endRef=referenceShiftedIntegral(polynomial.endCoefficients,0,polynomial.length,rz,
                        polynomial.radiusEnd,polynomial.radialTangent?-1:0,1);
                    peak=std::max(peak,normalized(total,beginRef,atom.area,1e-12,"original polynomial begin-chart GL3"));
                    peak=std::max(peak,normalized(total,endRef,atom.area,1e-12,"original polynomial reverse-chart GL3"));++rows;
                }
                const auto& report=rates.spatialReport();
                check(report.sourceQueries==reconstruction.atoms().size()&&report.ordinaryRanges==report.sourceQueries
                    &&report.events==0&&report.coincidentRanges==0&&report.refinements==0
                    &&report.positiveIntegrationErrorEstimate==0,"pure atom report invents root/approximation work");
                ++scenes;
            }
        }
    }
    check(scenes==8&&rows>0&&reverseCharts>0,"pure signed XY/RZ atom coverage incomplete");
    std::cout<<"CERTIFIED_SINGLE_PHASE_RATES,"<<scenes<<','<<rows<<','<<reverseCharts<<','<<peak<<",PASS\n";
}
// Fresh publications for each thread count exercise cold same-scalar caches.
// Callback consumers remain serial; comparisons exclude only opaque IDs.
void parallelTimeSamples()
{
    materialAmrSubsetWitness();
    certifiedSinglePhaseRates();
    struct Threads {
#ifdef _OPENMP
        int count=omp_get_max_threads(),dynamic=omp_get_dynamic(),maxActiveLevels=omp_get_max_active_levels();
        explicit Threads(int requested){omp_set_dynamic(0);omp_set_max_active_levels(std::min(maxActiveLevels,1));omp_set_num_threads(requested);}
        ~Threads(){omp_set_num_threads(count);omp_set_dynamic(dynamic);omp_set_max_active_levels(maxActiveLevels);}
#else
        explicit Threads(int){}
#endif
    };
    struct Cold {
        std::unique_ptr<QuadtreeMesh> mesh;
        HermiteGalsPublication original;
        HydroMaterialTransport material;
        std::string budgetFailure;
        int callbacks=0;
    };
    using TimeRule=HydroSpaceTimeIntegrationRule;
    check(HydroSpaceTimeConfig{}.integrationRule==TimeRule::EmbeddedGaussLegendre35
        &&HydroSpaceTimeReport{}.integrationRule==TimeRule::EmbeddedGaussLegendre35
        &&HydroSpaceTimeProgress{}.integrationRule==TimeRule::EmbeddedGaussLegendre35,
        "generic time quadrature/report defaults changed");
    const auto prepare=[&](int threads,bool prescribed,TimeRule rule=TimeRule::EmbeddedGaussLegendre35){
        Threads team(threads);auto mesh=makeMesh(4);const auto metrics=metricsFor(*mesh,false);
        const auto geometry=HydroMomentumGeometry::build(*mesh,metrics);const auto& layout=geometry.layout();
        const auto reconstruction=HydroAdvectorReconstruction::prepare(*mesh,geometry);
        Case shape;shape.pureFirstJet=true;shape.cx=.43;shape.cy=.51;shape.radius=.19;
        if(prescribed){shape.planar=true;shape.planeOffset=2;}
        const auto original=initial(*mesh,shape);const auto saved=original.state().primary();const auto savedClosure=original.state().closure();
        const auto identity=original.state().identity();
        const auto psi=[](Vector2D p){const double x=p.x*(1-p.x),y=p.y*(1-p.y);return 16*x*x*y*y*(1+.3*p.x+.7*p.y);};
        std::vector<int> walls;std::vector<double> normal(layout.velocityCount());
        for(int u=0;u<layout.velocityCount();++u){const auto face=layout.domain().face(*mesh,layout.nativeFace(u));const auto ends=faceBounds(*mesh,face);
            normal[u]=(face.horizontal?-1:1)*(psi(ends[1])-psi(ends[0]))/metrics.faceMeasures()[u];
            if(face.outwardSign!=0){check(normal[u]==0,"parallel material closed input leaked");walls.push_back(u);}}
        const auto law=HydroAdvectorBoundary::impermeable(*mesh,geometry,walls);
        const HydroNativeValues before(layout,normal,std::vector<double>(layout.pressureCount()));
        for(auto& value:normal)value*=1.1;
        const HydroNativeValues after(layout,normal,std::vector<double>(layout.pressureCount()));
        const auto previous=reconstruction.evaluate(*mesh,before,identity,original.state().timeLayer(),law);
        const auto arrival=reconstruction.evaluate(*mesh,after,identity,{1,0,.002,Position::End},law);
        const auto motion=HydroCharacteristicGeometry::prepare(*mesh,reconstruction).linear(*mesh,previous,arrival);
        HydroSpaceTimeConfig config;config.initialPanels=1;config.characteristic.maxStep=.0005;config.integrationRule=rule;
        config.enablePointParallel=false;config.enableNumericalTransportSampling=false; // Keep original outer-time exact 1/4 contract.
        config.progress=[&](const HydroSpaceTimeProgress& q){
#ifdef _OPENMP
            check(!omp_in_parallel(),"quadrature progress callback entered parallel team");
#endif
            check(q.integrationRule==rule,"progress lost the selected actual time rule");
        };
        int callbacks=0;HermiteGalsBoundaryFirstJetQuery boundary;
        if(prescribed)boundary=[&](Vector2D point,HermiteGalsTimeLayer time,Boundary){
#ifdef _OPENMP
            check(!omp_in_parallel(),"external boundary callback entered parallel team");
#endif
            ++callbacks;const auto jet=shape.profile(time.time())(point);return HermiteGalsBoundaryFirstJet{jet.value,{jet.dx,jet.dy}};
        };
        const auto material=HydroMaterialTransport::prepareSpaceTime(*mesh,original,motion,1000,1,config,boundary);
        const auto sources=material.sources();const auto clock=material.timeLayer();
        check(material.isSpaceTime()&&material.characteristicSource().sameIdentity(motion)
            &&material.liquidDensity()==1000&&material.gasDensity()==1
            &&sources.begin==identity&&sources.midpoint==material.midpointRate().geometry().geometryIdentity()
            &&sources.end==material.endGeometry().geometryIdentity()
            &&sources.advector==material.midpointRate().spatialAdvector().identity()
            &&clock.step==1&&clock.begin==0&&clock.end==.002&&clock.position==Position::Midpoint,
            "cold time quadrature lost actual phase densities/source/full clock");
        if(threads==1 && !prescribed && rule==TimeRule::EmbeddedSelectiveBoole){
            auto fullConfig=config;fullConfig.enableAtomLocalTimeIntegration=false;
            const auto fullReference=HydroMaterialTransport::prepareSpaceTime(*mesh,original,motion,1000,1,fullConfig,boundary);
            const auto parent=motion.arrival().timeLayer();const long double actualDt=parent.end-parent.begin;
            const auto& atoms=reconstruction.atoms();
            const auto& previousPolynomials=motion.previous().atomPolynomials();
            const auto& arrivalPolynomials=motion.arrival().atomPolynomials();
            check(previousPolynomials.size()==atoms.size()&&arrivalPolynomials.size()==atoms.size(),
                "nonlinear full/local reference changed endpoint polynomial inventory");
            std::vector<long double> phaseWork(atoms.size());
            for(std::size_t i=0;i<atoms.size();++i){long double speed=0;
                for(const auto* polynomial:{&previousPolynomials[i],&arrivalPolynomials[i]}){
                    long double work=0,power=1;
                    for(const auto coefficient:polynomial->beginCoefficients){work+=std::abs(coefficient)*power;power*=polynomial->length;}
                    speed=std::max(speed,work);
                }
                phaseWork[i]=actualDt*atoms[i].area*speed;
                check(std::isfinite(phaseWork[i])&&phaseWork[i]>=0,"nonfinite actual nonlinear atom work");
            }
            check(config.relativeTolerance==1e-7 && actualDt==static_cast<long double>(.002),
                "nonlinear full/local reference changed original interval or phase budget");
            for(const auto* result:{&material,&fullReference}){
                const auto report=result->spaceTimeReport();const auto actual=result->timeLayer();const auto owner=result->sources();
                const auto& cost=report.samplingCost;
                check(std::isfinite(report.maxRelativePhaseErrorEstimate)&&report.maxRelativePhaseErrorEstimate>=0
                    &&report.maxRelativePhaseErrorEstimate<=config.relativeTolerance
                    &&std::isfinite(report.absolutePhaseErrorEstimate)&&report.absolutePhaseErrorEstimate>=0,
                    "nonlinear full/local route failed original cumulative phase budget");
                check(result->isSpaceTime()&&result->characteristicSource().sameIdentity(motion)
                    &&result->liquidDensity()==1000&&result->gasDensity()==1&&result->deltaTime()==parent.end-parent.begin
                    &&owner.begin==identity&&owner.midpoint==result->midpointRate().geometry().geometryIdentity()
                    &&owner.end==result->endGeometry().geometryIdentity()
                    &&owner.advector==result->midpointRate().spatialAdvector().identity()
                    &&actual.step==parent.step&&actual.begin==parent.begin&&actual.end==parent.end&&actual.position==Position::Midpoint
                    &&cost.sourceIdentity==identity&&cost.previousAdvectorIdentity==motion.previous().identity()
                    &&cost.arrivalAdvectorIdentity==motion.arrival().identity()
                    &&report.integrationRule==rule&&report.firstJetTimeProjection==HermiteGalsTimeProjection::ClosedProjectionCompensated,
                    "nonlinear full/local reference detached genuine source/full clock");
            }
            check(material.spaceTimeReport().samplingCost.atomLocal&&!fullReference.spaceTimeReport().samplingCost.atomLocal,
                "nonlinear full/local reference did not use independent integration routes");
            // Two independent routes retain their original cumulative gate.
            // Their difference uses the sum of those budgets and the original
            // 256-epsilon material rounding scale; zero work has no floor.
            const long double differenceTolerance=2*static_cast<long double>(config.relativeTolerance)
                +256*std::numeric_limits<double>::epsilon();
            long double peak=0;int cut=0;
            const auto compare=[&](long double candidate,long double reference,long double work,const std::string& description){
                check(std::isfinite(candidate)&&std::isfinite(reference)&&std::isfinite(work)&&work>=0,
                    "nonfinite nonlinear full/local "+description);
                const auto difference=std::abs(candidate-reference);
                if(work==0)check(difference==0,"nonlinear full/local zero-work "+description);
                else{const auto relative=difference/work;check(relative<=differenceTolerance,
                    "nonlinear full/local "+description+" exceeds original joint budgets");peak=std::max(peak,relative);}
            };
            check(material.phaseAtomIntegrals().size()==atoms.size()&&fullReference.phaseAtomIntegrals().size()==atoms.size()
                &&material.atomIntegrals().size()==atoms.size()&&fullReference.atomIntegrals().size()==atoms.size(),
                "nonlinear full/local atom integral inventory differs");
            const long double densityWork=material.liquidDensity()+material.gasDensity();
            for(std::size_t i=0;i<atoms.size();++i){const auto candidate=material.phaseAtomIntegrals()[i],reference=fullReference.phaseAtomIntegrals()[i];
                compare(candidate.liquid,reference.liquid,phaseWork[i],"liquid atom integral");
                compare(candidate.gas,reference.gas,phaseWork[i],"gas atom integral");
                compare(material.atomIntegrals()[i],fullReference.atomIntegrals()[i],densityWork*phaseWork[i],"material atom integral");
                if(reference.liquid!=0&&reference.gas!=0)++cut;
            }
            const auto compareFaces=[&](const auto& candidate,const auto& reference,const auto& groups,const std::string& description){
                check(candidate.size()==reference.size()&&candidate.size()==groups.size(),
                    "nonlinear full/local "+description+" inventory differs");
                for(std::size_t i=0;i<groups.size();++i){long double work=0;
                    for(const int atom:groups[i])work+=phaseWork.at(static_cast<std::size_t>(atom));
                    compare(candidate[i],reference[i],densityWork*work,description);
                }
            };
            compareFaces(material.faceIntegrals(),fullReference.faceIntegrals(),reconstruction.primalAtoms(),"primal material integral");
            compareFaces(material.dualFaceIntegrals(),fullReference.dualFaceIntegrals(),reconstruction.dualAtoms(),"dual material integral");
            check(cut>0,"nonlinear full/local reference did not exercise true two-phase atoms");
            const auto localDepth=material.spaceTimeReport().deepestPanel,fullDepth=fullReference.spaceTimeReport().deepestPanel;
            std::cout<<"NONLINEAR_LOCAL_FULL_TIME,"<<atoms.size()<<','<<cut<<','<<peak
                <<",LOCAL_DEPTH,"<<localDepth<<",FULL_DEPTH,"<<fullDepth<<','
                <<(localDepth>0?"LOCAL_REFINED":"LOCAL_NOT_REFINED")
                <<",ACTUAL_ENDPOINT_WORK,ORIGINAL_CUMULATIVE_BUDGET,SOURCE_FULL_CLOCK,PHASE_MATERIAL_PRIMAL_DUAL,PASS\n";
        }
        check(material.spaceTimeReport().integrationRule==rule,"report lost the selected actual time rule");
        check(material.spaceTimeReport().geometrySamples>=(rule==TimeRule::EmbeddedGaussLegendre35?8u:4u),
            "actual panel quadrature samples missing");
        std::size_t expectedBatch=0;
#ifdef _OPENMP
        if(threads>1&&!prescribed)expectedBatch=rule==TimeRule::EmbeddedGaussLegendre35?8:4;
#endif
        check(material.spaceTimeReport().timings.fullInternalBatchTasks==expectedBatch,"actual full/internal batch execution differs from closed-source eligibility");
        check(material.spaceTimeReport().firstJetTimeProjection==(prescribed?HermiteGalsTimeProjection::CharacteristicPullback:
            HermiteGalsTimeProjection::ClosedProjectionCompensated),"declared boundary time rule changed");
        if(threads==1 && !prescribed && rule==TimeRule::EmbeddedGaussLegendre35){firstJetScalarSampleWitness(*mesh,original,motion,config,shape);
            nativeSupportWitness(*mesh,reconstruction);}
        auto limited=config;limited.maxGeometrySamples=3;std::string failure;
        try{HydroMaterialTransport::prepareSpaceTime(*mesh,original,motion,1000,1,limited,boundary);}
        catch(const std::exception& e){failure=e.what();}
        check(failure.find("spacetime geometry sample budget exhausted")!=std::string::npos,"ordered geometry budget refusal changed");
        check(original.state().identity()==identity&&original.state().closure().vertexMixedDerivatives==savedClosure.vertexMixedDerivatives
            &&original.state().closure().cellBubbleCoefficients==savedClosure.cellBubbleCoefficients,"failed panel modified immutable source/closure");
        const auto& current=original.state().primary();
        for(std::size_t i=0;i<saved.vertices.size();++i){const auto& a=saved.vertices[i];const auto& b=current.vertices[i];
            check(a.vertex==b.vertex&&a.phi==b.phi&&a.dx==b.dx&&a.dy==b.dy,"failed panel modified original vertex");}
        for(std::size_t i=0;i<saved.edges.size();++i){const auto& a=saved.edges[i];const auto& b=current.edges[i];
            check(a.vertices==b.vertices&&a.phi==b.phi&&a.dn==b.dn,"failed panel modified original edge");}
        return Cold{std::move(mesh),original,material,failure,callbacks};
    };
    const auto serial=prepare(1,false),parallel=prepare(4,false);
    const auto sameGeometry=[](const HydroPhaseGeometry& a,const HydroPhaseGeometry& b){
        const auto& x=a.publication().state();const auto& y=b.publication().state();const auto& p=x.primary();const auto& q=y.primary();
        const auto tx=x.timeLayer(),ty=y.timeLayer();check(tx.step==ty.step&&tx.begin==ty.begin&&tx.end==ty.end&&tx.position==ty.position,"parallel geometry full clock changed");
        check(p.geometryCells==q.geometryCells&&p.transport==q.transport&&p.measure==q.measure&&p.phaseVolumes.empty()&&q.phaseVolumes.empty(),"parallel primary layout changed");
        check(p.vertices.size()==q.vertices.size()&&p.edges.size()==q.edges.size(),"parallel FirstJet size changed");
        for(std::size_t i=0;i<p.vertices.size();++i){const auto& v=p.vertices[i];const auto& w=q.vertices[i];check(v.vertex==w.vertex&&v.phi==w.phi&&v.dx==w.dx&&v.dy==w.dy,"parallel vertex changed");}
        for(std::size_t i=0;i<p.edges.size();++i){const auto& e=p.edges[i];const auto& f=q.edges[i];check(e.vertices==f.vertices&&e.phi==f.phi&&e.dn==f.dn,"parallel edge changed");}
        check(x.closure().vertexMixedDerivatives==y.closure().vertexMixedDerivatives&&x.closure().cellBubbleCoefficients==y.closure().cellBubbleCoefficients,"parallel full closure changed");
        check(a.cellVolumes().size()==b.cellVolumes().size(),"parallel phase volume inventory changed");
        for(std::size_t i=0;i<a.cellVolumes().size();++i)check(a.cellVolumes()[i].liquid==b.cellVolumes()[i].liquid&&a.cellVolumes()[i].gas==b.cellVolumes()[i].gas,"parallel true phase volume changed");
    };
    const auto sameCold=[&](const Cold& first,const Cold& second){
        const auto& a=first.material;const auto& b=second.material;
        sameGeometry(a.beginGeometry(),b.beginGeometry());sameGeometry(a.midpointRate().geometry(),b.midpointRate().geometry());sameGeometry(a.endGeometry(),b.endGeometry());
        check(a.atomIntegrals()==b.atomIntegrals()&&a.faceIntegrals()==b.faceIntegrals()&&a.dualFaceIntegrals()==b.dualFaceIntegrals(),"parallel unique material flux changed");
        check(a.phaseAtomIntegrals().size()==b.phaseAtomIntegrals().size(),"parallel phase integral inventory changed");
        for(std::size_t i=0;i<a.phaseAtomIntegrals().size();++i)check(a.phaseAtomIntegrals()[i].liquid==b.phaseAtomIntegrals()[i].liquid
            &&a.phaseAtomIntegrals()[i].gas==b.phaseAtomIntegrals()[i].gas,"parallel phase integrals changed");
        check(a.cells().size()==b.cells().size(),"parallel material cell inventory changed");
        for(std::size_t i=0;i<a.cells().size();++i){const auto& x=a.cells()[i];const auto& y=b.cells()[i];
            check(x.massBegin==y.massBegin&&x.massMidpoint==y.massMidpoint&&x.massEnd==y.massEnd&&x.massChange==y.massChange
                &&x.outwardIntegral==y.outwardIntegral&&x.integratedDefect==y.integratedDefect&&x.defectRate==y.defectRate,"parallel true mass/R changed");}
        check(a.spaceTimeReport().geometrySamples==b.spaceTimeReport().geometrySamples&&a.spaceTimeReport().acceptedPanels==b.spaceTimeReport().acceptedPanels
            &&a.spaceTimeReport().deepestPanel==b.spaceTimeReport().deepestPanel&&a.spaceTimeReport().absolutePhaseErrorEstimate==b.spaceTimeReport().absolutePhaseErrorEstimate
            &&a.spaceTimeReport().maxRelativePhaseErrorEstimate==b.spaceTimeReport().maxRelativePhaseErrorEstimate,"parallel adaptive decisions/counters changed");
        check(first.budgetFailure==second.budgetFailure,"thread count changed first geometry budget exception");
        check(a.spaceTimeReport().integrationRule==b.spaceTimeReport().integrationRule,"parallel actual time rule changed");
    };
    sameCold(serial,parallel);
    const auto booleSerial=prepare(1,false,TimeRule::EmbeddedSimpsonBoole),booleParallel=prepare(4,false,TimeRule::EmbeddedSimpsonBoole);
    sameCold(booleSerial,booleParallel);
    const auto selectiveSerial=prepare(1,false,TimeRule::EmbeddedSelectiveBoole);
    const auto selectiveParallel=prepare(4,false,TimeRule::EmbeddedSelectiveBoole);
    sameCold(selectiveSerial,selectiveParallel);
    const auto callback=prepare(4,true);check(callback.callbacks>0,"serial external boundary callback was not exercised");
    // Reuse the existing uniform-XY moving-plane analytic scene once at
    // rho_l/rho_g=1000/1. The front crosses x=.5 inside this actual interval;
    // independent phase residence times test temporal accuracy, while the
    // original identity helper tests unique native/dual I and R at 1e-12.
    {
        Threads team(4);auto mesh=makeMesh(4);const auto metrics=metricsFor(*mesh,false);
        const auto geometry=HydroMomentumGeometry::build(*mesh,metrics);
        const auto reconstruction=HydroAdvectorReconstruction::prepare(*mesh,geometry);
        Case plane;plane.pureFirstJet=true;plane.planar=true;plane.planeOffset=.48137;plane.motion=1;plane.speed=.2;
        const double dt=.13;const auto original=initial(*mesh,plane);const auto identity=original.state().identity();
        const auto u=velocity(*mesh,geometry.layout(),{.2,0});
        const auto before=reconstruction.evaluate(*mesh,u,identity,original.state().timeLayer());
        const auto after=reconstruction.evaluate(*mesh,u,identity,{1,0,dt,Position::End});
        const auto motion=HydroCharacteristicGeometry::prepare(*mesh,reconstruction).linear(*mesh,before,after);
        HydroSpaceTimeConfig config;config.initialPanels=1;config.characteristic.maxStep=.02;
        config.integrationRule=TimeRule::EmbeddedSimpsonBoole;
        HermiteGalsPublicationConfig publication;publication.volume.symmetryPlanes=plane.symmetry(*mesh);
        const HermiteGalsBoundaryFirstJetQuery boundary=[plane](Vector2D p,HermiteGalsTimeLayer t,Boundary){
#ifdef _OPENMP
            check(!omp_in_parallel(),"analytic plane outer callback entered parallel team");
#endif
            const auto j=plane.profile(t.time())(p);return HermiteGalsBoundaryFirstJet{j.value,{j.dx,j.dy}};
        };
        const auto material=HydroMaterialTransport::prepareSpaceTime(*mesh,original,motion,1000,1,config,boundary,publication);
        const auto sources=material.sources();const auto clock=material.timeLayer();
        check(material.spaceTimeReport().integrationRule==TimeRule::EmbeddedSimpsonBoole
            &&material.characteristicSource().sameIdentity(motion)&&sources.begin==identity
            &&sources.midpoint==material.midpointRate().geometry().geometryIdentity()
            &&sources.end==material.endGeometry().geometryIdentity()
            &&clock.step==1&&clock.begin==0&&clock.end==dt&&clock.position==Position::Midpoint,
            "analytic Boole plane lost actual rule/source/full clock");
        std::vector<long double> expectedMass(reconstruction.atoms().size());long double peak=0;int crossed=0;
        for(std::size_t i=0;i<reconstruction.atoms().size();++i){const auto& atom=reconstruction.atoms()[i];
            const double speed=atom.normalAxis==0?.2:0;
            const double liquidTime=std::clamp(dt-(atom.begin.x-plane.planeOffset)/.2,0.,dt);
            const long double total=atom.area*speed*dt,liquid=atom.area*speed*liquidTime,gas=total-liquid;
            const auto actual=material.phaseAtomIntegrals()[i];const long double phaseWork=atom.area*.2*dt;
            peak=std::max({peak,normalized(actual.liquid,liquid,phaseWork,1e-7,"Boole analytic moving-plane liquid integral"),
                normalized(actual.gas,gas,phaseWork,1e-7,"Boole analytic moving-plane gas integral")});
            normalized(actual.liquid+actual.gas,total,phaseWork,1e-12,"Boole analytic original signed total integral");
            expectedMass[i]=1000*liquid+gas;
            normalized(material.atomIntegrals()[i],expectedMass[i],1000*phaseWork,1e-7,"Boole analytic density-weighted atom integral");
            if(atom.normalAxis==0&&liquidTime>0&&liquidTime<dt)++crossed;
        }
        check(crossed>0,"Boole analytic witness did not cross an actual phase atom");
        for(int f=0;f<geometry.layout().velocityCount();++f){long double expected=0,work=0;
            for(int atom:reconstruction.primalAtoms()[f]){expected+=expectedMass[atom];work+=1000*reconstruction.atoms()[atom].area*.2*dt;}
            normalized(material.faceIntegrals()[f],expected,work,1e-7,"Boole analytic native mass integral");}
        const auto mass=[&](const HydroPhaseGeometry& phase){return HydroMomentumMass::build(*mesh,geometry,phase,1000,1);};
        const auto dual=HydroDualMaterial::prepare(*mesh,material,mass(material.beginGeometry()),mass(material.midpointRate().geometry()),mass(material.endGeometry()));
        for(std::size_t f=0;f<geometry.segments().size();++f){long double expected=0,work=0;
            for(int atom:reconstruction.dualAtoms()[f]){expected+=expectedMass[atom];work+=1000*reconstruction.atoms()[atom].area*.2*dt;}
            normalized(dual.faceIntegrals()[f],expected,work,1e-7,"Boole analytic dual mass integral");}
        spaceTimeIdentities(*mesh,material,dual);
        // The same ORIGINAL moving plane/source also checks the new rule's
        // full prescribed-boundary route. No manufactured time-rate callback
        // replaces a geometric query; phase residence times remain analytic.
        auto selectiveConfig=config;
        selectiveConfig.integrationRule=TimeRule::EmbeddedSelectiveBoole;
        selectiveConfig.enableAtomLocalTimeIntegration=false;
        const auto selective=HydroMaterialTransport::prepareSpaceTime(*mesh,original,motion,1000,1,
            selectiveConfig,boundary,publication);
        const auto selectiveSources=selective.sources();const auto selectiveClock=selective.timeLayer();
        check(selective.spaceTimeReport().integrationRule==TimeRule::EmbeddedSelectiveBoole
            &&!selective.spaceTimeReport().samplingCost.atomLocal
            &&selective.spaceTimeReport().firstJetTimeProjection==HermiteGalsTimeProjection::CharacteristicPullback
            &&selective.characteristicSource().sameIdentity(motion)&&selectiveSources.begin==identity
            &&selectiveSources.midpoint==selective.midpointRate().geometry().geometryIdentity()
            &&selectiveSources.end==selective.endGeometry().geometryIdentity()
            &&selectiveClock.step==1&&selectiveClock.begin==0&&selectiveClock.end==dt
            &&selectiveClock.position==Position::Midpoint,"selective Boole lost actual full boundary source/clock");
        long double selectivePeak=0,oldNewDifference=0;
        for(std::size_t i=0;i<reconstruction.atoms().size();++i) {
            const auto& atom=reconstruction.atoms()[i];const double speed=atom.normalAxis==0?.2:0;
            const double liquidTime=std::clamp(dt-(atom.begin.x-plane.planeOffset)/.2,0.,dt);
            const long double total=atom.area*speed*dt,liquid=atom.area*speed*liquidTime,gas=total-liquid;
            const long double phaseWork=atom.area*.2*dt;const auto q=selective.phaseAtomIntegrals()[i];
            selectivePeak=std::max({selectivePeak,
                normalized(q.liquid,liquid,phaseWork,1e-7,"selective Boole analytic moving-plane liquid integral"),
                normalized(q.gas,gas,phaseWork,1e-7,"selective Boole analytic moving-plane gas integral"),
                normalized(selective.atomIntegrals()[i],expectedMass[i],1000*phaseWork,1e-7,
                    "selective Boole analytic density-weighted atom integral")});
            normalized(q.liquid+q.gas,total,phaseWork,1e-12,"selective Boole original signed total integral");
            const auto old=material.phaseAtomIntegrals()[i];
            const auto scale=phaseWork>0?phaseWork:1;
            oldNewDifference=std::max({oldNewDifference,std::abs(q.liquid-old.liquid)/scale,
                std::abs(q.gas-old.gas)/scale});
        }
        // Each route has its independent unchanged 1e-7 analytic gate above;
        // their difference is reported separately, never used as a true flow.
        for(int f=0;f<geometry.layout().velocityCount();++f) {
            long double expected=0,work=0;
            for(int atom:reconstruction.primalAtoms()[f]) {
                expected+=expectedMass[atom];work+=1000*reconstruction.atoms()[atom].area*.2*dt;
            }
            normalized(selective.faceIntegrals()[f],expected,work,1e-7,"selective Boole analytic native mass integral");
        }
        for(std::size_t f=0;f<reconstruction.dualAtoms().size();++f) {
            long double expected=0,work=0;
            for(int atom:reconstruction.dualAtoms()[f]) {
                expected+=expectedMass[atom];work+=1000*reconstruction.atoms()[atom].area*.2*dt;
            }
            normalized(selective.dualFaceIntegrals()[f],expected,work,1e-7,"selective Boole analytic dual mass integral");
        }
        rejects([&]{HydroMaterialTransport::prepareSpaceTime(*mesh,original,motion,1000,1,
            selectiveConfig,boundary,publication,&material);},"previous time plan: different time integration rule");
        rejects([&]{HydroMaterialTransport::prepareSpaceTime(*mesh,original,motion,1000,1,
            config,boundary,publication,&selective);},"previous time plan: different time integration rule");
        std::cout<<"SELECTIVE_BOOLE_TIME_SAMPLES,COLD_1_4_EXACT,FULL_PRESCRIBED,RHO_1000,ANALYTIC_CROSSING,"
            <<crossed<<','<<selectivePeak<<','<<oldNewDifference<<','<<selective.spaceTimeReport().geometrySamples
            <<','<<selective.spaceTimeReport().acceptedPanels<<",PASS\n";
        std::cout<<"BOOLE_TIME_SAMPLES,COLD_1_4_EXACT,RHO_1000,ANALYTIC_CROSSING,"<<crossed<<','<<peak
            <<','<<material.spaceTimeReport().geometrySamples<<','<<material.spaceTimeReport().acceptedPanels<<",PASS\n";
    }
    std::cout<<"PARALLEL_TIME_SAMPLES,COLD_1_4,EXACT,ORDERED_BUDGET,SERIAL_CALLBACK,PASS\n";
}
// One actual moving gas bubble in each existing coordinate system. The callback
// supplies only inlet arrival data; the outlet must retain the original
// scalar's characteristic pullback, with no second material flux or R repair.
void inflowBoundaryOnly()
{
    using Treatment=HermiteGalsBoundaryTreatment;
    check(HydroSpaceTimeConfig{}.firstJetBoundaryTreatment==Treatment::PrescribedOuterFirstJet,
        "legacy callback boundary treatment default changed");
    const auto sameClock=[](HermiteGalsTimeLayer a,HermiteGalsTimeLayer b){
        return a.step==b.step&&a.begin==b.begin&&a.end==b.end&&a.position==b.position;};
    for(bool rz:{false,true}){
        auto mesh=makeMesh(4);Case shape;shape.pureFirstJet=true;shape.radial=rz;
        shape.radius=rz?.27:.19;shape.motion=1;shape.speed=.2;
        // A moving closed gas bubble in liquid: no non-axis wall contact or
        // exact scalar symmetry is asserted on naturally transported walls.
        const auto profile=[shape](double time)->Profile{return [shape,time](Vector2D q){
            const auto jet=shape.profile(time)(q);return Jet{-jet.value,-jet.dx,-jet.dy,-jet.dxy};};};
        const auto f=profile(0);
        const auto seed=stateFor(*mesh,f,rz,{},[](int,Vector2D,Vector2D){return 0.;},{},HermiteGalsTransport::FirstJet);
        HermiteGalsCandidate candidate(*mesh,{},seed.primary());const auto mixed=closeHermiteGalsMixedAmr(*mesh,candidate);
        const auto predicted=centers(*mesh,f);
        const auto original=publishHermiteGalsFirstJetGeometry(*mesh,std::move(candidate),mixed,{candidate.identity(),predicted});
        const auto retained=original;const auto sourceId=original.state().identity();
        const auto* sourceScalar=&original.state().geometry().snapshot();
        const auto metrics=metricsFor(*mesh,rz);const auto geometry=HydroMomentumGeometry::build(*mesh,metrics);
        const auto& layout=geometry.layout();const auto reconstruction=HydroAdvectorReconstruction::prepare(*mesh,geometry);
        const auto values=velocity(*mesh,layout,rz?Vector2D{0,shape.speed}:Vector2D{shape.speed,0});
        const double dt=.01;const HermiteGalsTimeLayer clock{1,0,dt,Position::End};
        check(dt*shape.speed<.5*std::min(layout.domain().spacing().x,layout.domain().spacing().y),
            "inflow witness permits an unsupported interior boundary event");
        const auto previous=reconstruction.evaluate(*mesh,values,sourceId,original.state().timeLayer());
        const auto arrival=reconstruction.evaluate(*mesh,values,sourceId,clock);
        const auto motion=HydroCharacteristicGeometry::prepare(*mesh,reconstruction).linear(*mesh,previous,arrival);
        HydroSpaceTimeConfig config;config.initialPanels=1;config.characteristic.maxStep=.0025;
        config.firstJetTimeScheme=HydroFirstJetTimeScheme::CharacteristicPullback;
        config.firstJetBoundaryTreatment=Treatment::InflowHistory;
        HermiteGalsPublicationConfig publication; // RZ axis parity comes from the original radial snapshot.
        std::size_t inletQueries=0;
        const HermiteGalsBoundaryFirstJetQuery inlet=[&](Vector2D q,HermiteGalsTimeLayer t,Boundary side){
            const bool coordinate=rz?q.y==mesh->minCoord.y:q.x==mesh->minCoord.x;
            const bool actualSide=rz?(side==Boundary::south||side==Boundary::southWest||side==Boundary::southEast)
                :(side==Boundary::west||side==Boundary::southWest||side==Boundary::northWest);
            check(coordinate&&actualSide,"inflow-only callback queried an outlet, wall, or interior point");
            check(t.step==clock.step&&t.begin==clock.begin&&t.end>t.begin&&t.end<=clock.end
                &&t.time()>clock.begin&&t.time()<=clock.end,"inflow-only callback lost its actual arrival full clock");
            ++inletQueries;const auto jet=profile(t.time())(q);
            return HermiteGalsBoundaryFirstJet{jet.value,{jet.dx,jet.dy}};
        };
        const auto material=HydroMaterialTransport::prepareSpaceTime(*mesh,original,motion,1000,1,config,inlet,publication);
        check(inletQueries>0&&material.spaceTimeReport().firstJetBoundaryTreatment==Treatment::InflowHistory
            &&material.spaceTimeReport().firstJetTimeProjection==HermiteGalsTimeProjection::CharacteristicPullback
            &&material.characteristicSource().sameIdentity(motion),"material did not retain the actual inlet-only boundary semantics/source");
        check(material.sources().begin==sourceId&&material.sources().midpoint==material.midpointRate().geometry().geometryIdentity()
            &&material.sources().end==material.endGeometry().geometryIdentity(),"inflow material lost its three actual scalar owners");
        long double peak=0;std::size_t outletJets=0;
        for(const auto* phase:std::array<const HydroPhaseGeometry*,2>{{&material.midpointRate().geometry(),&material.endGeometry()}}){
            const auto expectedClock=HermiteGalsTimeLayer{clock.step,clock.begin,clock.end,
                phase==&material.endGeometry()?Position::End:Position::Midpoint};
            check(sameClock(phase->timeLayer(),expectedClock),"inflow material lost its actual midpoint/end full clock");
            const auto actualShape=shapeAt(shape,phase->timeLayer().time());const auto& state=phase->publication().state();
            check(state.primary().transport==HermiteGalsTransport::FirstJet&&state.primary().phaseVolumes.empty(),
                "inflow transport introduced an inventory or another scalar method");
            for(const auto& v:state.primary().vertices){const auto q=point(*mesh,v.vertex);const auto jet=profile(phase->timeLayer().time())(q);
                peak=std::max({peak,normalized(v.phi,jet.value,1,1e-12,"inflow actual midpoint/end vertex phi"),
                    normalized(v.dx,jet.dx,1,1e-12,"inflow actual midpoint/end vertex dx"),
                    normalized(v.dy,jet.dy,1,1e-12,"inflow actual midpoint/end vertex dy")});
                if(rz?q.y==mesh->maxCoord.y:q.x==mesh->maxCoord.x){
                    const auto trace=motion.traceBackward(q,phase->timeLayer().time(),clock.begin,config.characteristic);
                    check(!trace.boundary()&&trace.jacobian(),"outlet did not use the original interior characteristic map");
                    const auto old=sourceScalar->scalarJetNative(trace.point());
                    const auto gradient=trace.jacobian()->transposeTimes({old.dx,old.dy});
                    normalized(v.phi,old.value,1,1e-12,"outlet natural original scalar pullback");
                    normalized(v.dx,gradient.x,1,1e-12,"outlet natural original dx pullback");
                    normalized(v.dy,gradient.y,1,1e-12,"outlet natural original dy pullback");++outletJets;}
            }
            for(const auto& e:state.primary().edges){const auto a=point(*mesh,e.vertices[0]),b=point(*mesh,e.vertices[1]);
                const auto q=a+(b-a)*.5;const auto jet=profile(phase->timeLayer().time())(q);
                normalized(e.phi,jet.value,1,1e-12,"inflow actual midpoint/end edge phi");
                normalized(e.dn,a.y==b.y?jet.dy:jet.dx,1,1e-12,"inflow actual midpoint/end edge normal derivative");}
            for(int p=0;p<layout.pressureCount();++p){const auto bounds=sourceScalar->nativeCellBounds(layout.nativeCell(p));
                const auto cap=capacity(bounds[0],bounds[1],rz),liquid=cap-liquidCell(bounds[0],bounds[1],actualShape);
                peak=std::max(peak,normalized(phase->cellVolumes()[p].liquid,liquid,cap,1e-7,"inflow analytic true moving phase volume"));
                normalized(phase->cellVolumes()[p].liquid+phase->cellVolumes()[p].gas,cap,cap,1e-12,"inflow actual true cell capacity");}
        }
        check(outletJets>0,"inflow witness did not inspect actual outlet FirstJets");
        long double volumeBoundary=0,volumeWork=0;int inletFaces=0,outletFaces=0;
        for(int u=0;u<layout.velocityCount();++u){const auto face=layout.domain().face(*mesh,layout.nativeFace(u));
            if(face.negativeCell>=0&&face.positiveCell>=0)continue;
            const auto area=faceCapacity(*mesh,face,rz);const auto normal=values.normalVelocity()[u];
            const int outwardSign=face.negativeCell>=0?1:-1;volumeBoundary+=outwardSign*area*normal;volumeWork+=std::abs(area*normal);
            const bool port=face.horizontal==rz;const bool entering=port&&outwardSign<0;
            const long double expected=port?1000.L*area*shape.speed*dt:0;
            peak=std::max(peak,normalized(material.faceIntegrals()[u],expected,std::abs(expected),1e-12,"independent actual inlet/outlet native material integral"));
            if(entering){check(expected>0&&material.faceIntegrals()[u]>0,"inlet has no actual liquid supply");++inletFaces;}
            else if(port){check(expected>0&&material.faceIntegrals()[u]>0,"outlet has no actual transported liquid");++outletFaces;}
        }
        normalized(volumeBoundary,0,volumeWork,1e-12,"independent true balanced volume boundary flow");
        long double phaseMovement=0;
        for(int p=0;p<layout.pressureCount();++p)phaseMovement+=std::abs(material.endGeometry().cellVolumes()[p].liquid-material.beginGeometry().cellVolumes()[p].liquid);
        check(inletFaces>0&&outletFaces>0&&phaseMovement>0,
            "inflow witness did not move a true two-phase gas bubble through balanced liquid ports");
        const auto mass=[&](const HydroPhaseGeometry& phase){return HydroMomentumMass::build(*mesh,geometry,phase,1000,1);};
        const auto dual=HydroDualMaterial::prepare(*mesh,material,mass(material.beginGeometry()),mass(material.midpointRate().geometry()),mass(material.endGeometry()));
        check(dual.faceIntegrals()==material.dualFaceIntegrals(),"inflow dual consumer replaced the original unique F");
        spaceTimeIdentities(*mesh,material,dual);
        auto outer=config;outer.firstJetBoundaryTreatment=Treatment::PrescribedOuterFirstJet;
        rejects([&]{HydroMaterialTransport::prepareSpaceTime(*mesh,original,motion,1000,1,outer,inlet,publication,&material);},
            "previous time plan: different FirstJet boundary treatment");
        auto invalid=config;invalid.firstJetBoundaryTreatment=static_cast<Treatment>(999);
        rejects([&]{HydroMaterialTransport::prepareSpaceTime(*mesh,original,motion,1000,1,invalid,inlet,publication);},"invalid FirstJet boundary treatment");
        if(!rz){
            std::size_t outerQueries=0;const HermiteGalsBoundaryFirstJetQuery allOuter=[&](Vector2D q,HermiteGalsTimeLayer t,Boundary){
                check(q.x==0||q.x==1||q.y==0||q.y==1,"legacy callback queried an interior point");
                if(q.x==1)++outerQueries;const auto j=profile(t.time())(q);return HermiteGalsBoundaryFirstJet{j.value,{j.dx,j.dy}};};
            HydroSpaceTimeConfig defaults;defaults.initialPanels=1;defaults.characteristic.maxStep=.0025;
            const auto legacy=HydroMaterialTransport::prepareSpaceTime(*mesh,original,motion,1000,1,defaults,allOuter,publication);
            check(outerQueries>0&&legacy.spaceTimeReport().firstJetBoundaryTreatment==Treatment::PrescribedOuterFirstJet,
                "default callback no longer prescribes all actual outer arrivals");
            for(int p=0;p<layout.pressureCount();++p){const auto bounds=sourceScalar->nativeCellBounds(layout.nativeCell(p));const auto cap=capacity(bounds[0],bounds[1],false);
                normalized(legacy.endGeometry().cellVolumes()[p].liquid,material.endGeometry().cellVolumes()[p].liquid,cap,1e-7,"legacy all-outer moving gas bubble compatibility");}
            const auto zero=velocity(*mesh,layout,{0,0});
            const auto zeroBefore=reconstruction.evaluate(*mesh,zero,sourceId,original.state().timeLayer());
            const auto zeroAfter=reconstruction.evaluate(*mesh,zero,sourceId,clock);
            const auto zeroMotion=motion.geometry().linear(*mesh,zeroBefore,zeroAfter);
            const auto noCallback=HydroMaterialTransport::prepareSpaceTime(*mesh,original,zeroMotion,1000,1,defaults,{},publication);
            check(noCallback.spaceTimeReport().firstJetBoundaryTreatment==Treatment::InflowHistory,
                "callback-free default did not retain effective InflowHistory");
        }
        check(original.state().identity()==retained.state().identity()&&sourceId==original.state().identity()
            &&sourceScalar==&original.state().geometry().snapshot()&&&original.state().primary()==&retained.state().primary(),
            "boundary strategy preparation/rejection changed the original accepted scalar authority");
        std::cout<<"INFLOW_BOUNDARY,"<<(rz?"RZ":"XY")<<",rho=1000:1,"<<inletQueries<<','<<outletJets<<','<<inletFaces<<','<<outletFaces<<','<<peak<<",PASS\n";
    }
    std::cout<<"INFLOW_BOUNDARY_ONLY,XY_RZ,ORIGINAL_F_R_KEP,NATURAL_OUTLET,PASS\n";
}

int main(int argc,char** argv)
{
    try{
        if(argc==2&&std::string(argv[1])=="--inflow-boundary-only"){
            std::cout<<std::setprecision(17);inflowBoundaryOnly();return 0;}
        if(argc==2&&std::string(argv[1])=="--parallel-time-only"){
            std::cout<<std::setprecision(17);parallelTimeSamples();return 0;}
        if(argc==3&&std::string(argv[1])=="--spacetime-amr-crossing"){
            check(std::string(argv[2])=="0"||std::string(argv[2])=="1","invalid spacetime coordinate system");
            std::cout<<std::setprecision(17);spaceTimeMaterial(std::stoi(argv[2]),1,true);return 0;}
        if(argc==2&&std::string(argv[1])=="--spacetime-endpoints"){std::cout<<std::setprecision(17);spaceTimeEndpoints();return 0;}
        if((argc==2||argc==4)&&std::string(argv[1])=="--spacetime"){
            if(argc==4)check((std::string(argv[2])=="0"||std::string(argv[2])=="1")&&(std::string(argv[3])=="0"||std::string(argv[3])=="1"),"invalid spacetime layout");
            std::cout<<std::setprecision(17);spaceTimeMaterial(argc==4?std::stoi(argv[2]):-1,argc==4?std::stoi(argv[3]):-1);return 0;}
        check(argc==1,"unknown material contract arguments");
        std::cout<<std::setprecision(17);inflowBoundaryOnly();parallelTimeSamples();impermeableMaterial();
        for(bool rz:{false,true})for(int amr:{0,1})for(bool subset:{false,true}){
            auto mesh=makeMesh(4,amr);const auto metrics=metricsFor(*mesh,rz,subset);
            for(int kind=0;kind<4;++kind){Case s;s.pureFirstJet=true;s.radial=rz;s.amr=amr!=0;s.planar=kind!=0;s.planeOffset=kind==1?.43:kind==2?2:-1;
                const auto a=HydroPhaseGeometry::build(*mesh,metrics,publishAt(*mesh,s,{})),m=HydroPhaseGeometry::build(*mesh,metrics,publishAt(*mesh,s,{1,0,.04,Position::Midpoint})),b=HydroPhaseGeometry::build(*mesh,metrics,publishAt(*mesh,s,{1,0,.04,Position::End}));
                for(double rho:{1.,1000.,1e6})for(double sign:{0.,1.,-1.}){
                    const auto u=velocity(*mesh,metrics.layout(),rz?Vector2D{0,sign*.2}:Vector2D{sign*.3,-sign*.2});
                    const auto material=HydroMaterialTransport::prepare(*mesh,a,rate(*mesh,m,u),b,rho,1);inspect(*mesh,s,material,u);
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
        std::cout<<"H08 HYDRO midpoint material flux and nonconservative defect PASS\n";
        referencePolynomialContract();
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
