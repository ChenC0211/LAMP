#include "Modules/Interface/HermiteGalsTemporalPhaseFlux.h"
#include "Modules/Interface/HermiteGalsInventoryRegion.h"
#include "Modules/Interface/HermiteGalsAmrClosure.h"
#include "Core/Geometry/ActiveQuarticHermitePatch.h"
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Discrete/BoundedCirculation.h"
#include <algorithm>
#include <array>
#include <cfenv>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <type_traits>
#if defined(_M_X64) || (defined(__x86_64__) && defined(__SSE2__))
#include <xmmintrin.h>
#endif

namespace
{
using namespace Quantum::Interface;
using Velocity=HermiteGalsTemporalVelocitySource;
using Characteristic=HermiteGalsTemporalCharacteristicSource;
using Interval=HermiteGalsPhaseFluxIntervalSource;
using Coverage=HermiteGalsCharacteristicPhaseCoverage;
using Phase=HermiteGalsCertifiedPhase;
using Position=HermiteGalsTimeLayer::Position;
using Patch=Quantum::Geometry::ActiveQuarticHermitePatch;
using Jet=Patch::NodalDof;
using Profile=std::function<Jet(Vector2D)>;
using Inventory=std::function<double(int,Vector2D,Vector2D,double)>;
using Key=std::array<int,2>;
static_assert(!std::is_default_constructible_v<Coverage>);
static_assert(std::is_same_v<decltype(std::declval<Coverage&>().cells()),const std::vector<HermiteGalsCharacteristicPhaseCell>&>);
int queries=0,samples=0,traces=0,guards=0,thinEvents=0,weightedCases=0;
std::size_t checkedLeaves=0,maximumLeaves=0;
void check(bool ok,const std::string& why){if(!ok)throw std::runtime_error(why);}
template<class F> void rejects(F action,const char* expected)
{try{action();}catch(const std::exception& e){check(std::string(e.what()).find(expected)!=std::string::npos,e.what());++guards;return;}
    throw std::runtime_error(std::string("missing rejection: ")+expected);}
std::unique_ptr<QuadtreeMesh> makeMesh(bool amr)
{DictionaryMesh d;d.Nx=d.Ny=4;d.maxCellNum=1024;d.maxInitialIterations=0;auto m=std::make_unique<QuadtreeMesh>(d);
    if(amr){DictionaryMesh::RefineBox b;b.level=1;b.min={.25,.25};b.max={.75,.75};m->applyRefineBoxes({b});}return m;}
Vector2D point(const QuadtreeMesh& m,int v){return m.verts(m.vertIndexToLabel[v]).coord();}
bool fineSide(const QuadtreeMesh& m,Key key)
{const auto a=point(m,key[0]),b=point(m,key[1]);const bool h=a.y==b.y;
    for(int v:key){const auto owner=m.verts(m.vertIndexToLabel[v]).supCell;if(owner<0)continue;
        const auto& c=m.cells(owner);const auto lo=m.verts(c.cellToVert[0]).coord(),hi=m.verts(c.cellToVert[2]).coord();
        if(h && (a.y==lo.y || a.y==hi.y) && a.x>=lo.x && b.x<=hi.x && b.x-a.x<hi.x-lo.x)return true;
        if(!h && (a.x==lo.x || a.x==hi.x) && a.y>=lo.y && b.y<=hi.y && b.y-a.y<hi.y-lo.y)return true;}return false;}
// Analytic primary/center samples provide a structural Q4 source. Only the
// explicitly weighted fixtures bind analytic M; no interface is advanced here.
HermiteGalsState stateFor(const QuadtreeMesh& m,const Profile& f,bool radial,const Inventory& inventory={},bool complete=true)
{
    FaceSegmentGraph graph;graph.build(m);const auto domain=NativeLeafDomain::build(m,graph,radial);
    HermiteGalsPrimaryState p;Patch patch;p.measure=radial?HermiteGalsMeasure::AxisymmetricVolume:HermiteGalsMeasure::CartesianArea;
    patch.activeCells.assign(m.cellIndexNum(),1);patch.nodes.resize(m.vertIndexNum());patch.cellCenters.resize(m.cellIndexNum());
    for(int v=0;v<m.vertIndexNum();++v)if(m.verts(m.vertIndexToLabel[v]).type()==TVertType::non){const auto j=f(point(m,v));
        patch.nodes[v]=j;p.vertices.push_back({v,j.value,j.dx,j.dy});}
    std::set<Key> edges;std::set<int> activeVertices;
    for(int c=0;c<m.cellIndexNum();++c){const auto& cell=m.cells(m.cellIndexToLabel[c]);const auto a=m.verts(cell.cellToVert[0]).coord(),b=m.verts(cell.cellToVert[2]).coord();
        if(complete || c+1<m.cellIndexNum()){p.geometryCells.push_back(c);p.phaseVolumes.push_back({c,inventory?inventory(c,a,b,domain.cellVolumes()[c]):0.});}
        else patch.activeCells[c]=0;
        patch.cellCenters[c]=f(a+(b-a)*.5).value;
        if(!patch.activeCells[c])continue;
        for(int k=0;k<4;++k){Key key{{m.verts(cell.cellToVert[k]).index(),m.verts(cell.cellToVert[(k+1)%4]).index()}};
            activeVertices.insert(key[0]);
            const auto x=point(m,key[0]),y=point(m,key[1]);if(x.y==y.y?x.x>y.x:x.y>y.y)std::swap(key[0],key[1]);edges.insert(key);}}
    p.vertices.erase(std::remove_if(p.vertices.begin(),p.vertices.end(),[&](const auto& v){return activeVertices.find(v.vertex)==activeVertices.end();}),p.vertices.end());
    for(const auto key:edges)if(!fineSide(m,key)){const auto a=point(m,key[0]),b=point(m,key[1]);const auto j=f(a+(b-a)*.5);const double dn=a.y==b.y?j.dy:j.dx;
        p.edges.push_back({key,j.value,dn});patch.edgeMidpoints.push_back({key,j.value,dn});}
    HermiteGalsCandidate candidate(m,{2,.125,.25,Position::End},p);const auto mixed=closeHermiteGalsMixedAmr(m,candidate);
    HermiteGalsClosure closure;closure.vertexMixedDerivatives=mixed.vertexMixedDerivatives();
    for(int c:p.geometryCells)closure.cellBubbleCoefficients.push_back(patch.cellCenters[c]);
    for(std::size_t k=0;k<p.vertices.size();++k)patch.nodes[p.vertices[k].vertex].dxy=closure.vertexMixedDerivatives[k];
    candidate.bindGeometry(m,candidate.identity(),closure,LevelSetHermiteGeometry::Snapshot::fromNativeLeaves(m,patch,radial));return std::move(candidate).finish();
}
using Data=Quantum::Geometry::PolynomialVelocityData;
using Slab=Quantum::Geometry::PolynomialVelocitySlab;
using Bounds=Quantum::Geometry::PolynomialFaceIntegralBounds;
using Result=HermiteGalsBoundedFaceIntegral;
static_assert(!std::is_default_constructible_v<Result>);
static_assert(std::is_same_v<decltype(std::declval<Result&>().bounds()),const Bounds&>);
constexpr double oldTime=.25,endTime=.375,duration=.125;
const double pi=std::acos(-1.);
std::size_t totalPanels=0,totalEvaluated=0,maxPanels=0;
int integrals=0,coreCases=0,matrixFailures=0,tightTraces=0;
template<class F> void matrixCase(const char* name,F action)
{try{action();}catch(const std::exception& e){++matrixFailures;std::cout<<"MATRIX_FAILURE,"<<name<<','<<e.what()<<'\n';}}
Profile pure(int sign){return [=](Vector2D){return Jet{double(sign),0,0,0};};}
Profile plane(bool horizontal,double intercept,int sign=1)
{return [=](Vector2D p){return Jet{sign*((horizontal?p.y:p.x)-intercept),horizontal?0.:double(sign),horizontal?double(sign):0.,0.};};}
Data translation(bool horizontal,double speed=.25)
{Data d;d.previous[0]=horizontal?Vector2D{0,speed}:Vector2D{speed,0};d.arrival=d.previous;return d;}
Interval intervalFor(const QuadtreeMesh& mesh,const HermiteGalsState& old,bool amr,const Data& data)
{
    HermiteGalsCandidate next(mesh,{3,oldTime,endTime,Position::End},old.primary());
    const auto velocity=Velocity::polynomial(old,next.identity(),next.timeLayer(),data);
    const double h=amr?.125:.25;const Characteristic characteristic(mesh,{mesh.minCoord,mesh.maxCoord,{h,h},old.primary().measure==HermiteGalsMeasure::AxisymmetricVolume});
    return Interval(mesh,old,velocity,characteristic,oldTime,endTime);
}
struct Face{int row;Vector2D a,b;};
Face findFace(const QuadtreeMesh& mesh,bool radial,bool horizontal,double normal=.5,double tangent=.25)
{
    FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,radial);
    for(const auto face:domain.faces(mesh)){if(face.horizontal!=horizontal)continue;
        const auto& edge=horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
        auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
        if(horizontal?a.x>b.x:a.y>b.y)std::swap(a,b);
        if((horizontal?a.y:a.x)==normal && (horizontal?a.x:a.y)==tangent)return {face.nativeRow,a,b};}
    throw std::runtime_error("native fixture face missing");
}
double measure(const Face& f,bool radial)
{return f.a.y==f.b.y?(radial?pi*(f.b.x*f.b.x-f.a.x*f.a.x):f.b.x-f.a.x):(f.b.y-f.a.y)*(radial?2*pi*f.a.x:1);}
void contains(const Bounds& bounds,const std::array<double,4>& expected,double referenceError=0)
{
    const auto values=bounds.representative();double error=0;
    for(int k=0;k<4;++k){const auto b=bounds.volumes[k];
        check(b.lower<=expected[k]+referenceError && b.upper>=expected[k]-referenceError,"independent reference outside bound, component "+std::to_string(k));
        check(b.lower<=values[k] && values[k]<=b.upper,"representative outside proven interval");error+=std::abs(values[k]-expected[k]);}
    check(error<=bounds.errorBound()+4*referenceError,"representative error exceeds certificate");
    check(values[2]<=values[0] && values[3]<=values[1],"liquid gross not subset of total");
}
void verify(const QuadtreeMesh& mesh,const Result& result,const std::array<double,4>& expected,const char* name,double referenceError=0)
{
    ++integrals;result.validateMesh(mesh);contains(result.bounds(),expected,referenceError);
    const auto& r=result.report();check(r.absoluteErrorBound==result.bounds().errorBound() && r.absoluteErrorBound<=result.config().absoluteVolumeTolerance,"reported budget mismatch");
    check(r.acceptedPanels>0 && r.evaluatedPanels<=1+4*(r.acceptedPanels-1),"candidate certificate added hidden trial work");
    const auto v=result.bounds().representative();const auto& transfer=result.transfer().integratedVolumes();
    check(v[0]==transfer.integratedTotalPositiveVolume && v[1]==transfer.integratedTotalNegativeVolume && v[2]==transfer.integratedLiquidPositiveVolume && v[3]==transfer.integratedLiquidNegativeVolume,"transfer not same bounded integral");
    totalPanels+=r.acceptedPanels;totalEvaluated+=r.evaluatedPanels;maxPanels=std::max(maxPanels,r.acceptedPanels);
    std::cout<<"INTEGRAL,"<<name<<",panels="<<r.acceptedPanels<<",evaluated="<<r.evaluatedPanels<<",E="<<r.absoluteErrorBound<<",reference="<<expected[2]<<','<<expected[3]<<'\n';
}
void coreMoments()
{
    // Independent antiderivatives in the physical coordinate, nonunit basis,
    // cubic tangential dependence and linearly changing endpoint amplitude.
    for(bool radial:{false,true})for(bool horizontal:{false,true})for(int sign:{-1,1}){
        Data d;d.origin={.125,.125};d.scale={.5,.5};
        const int axis=horizontal?1:0;auto set=[&](Vector2D& c,double v){if(axis)c.y=v;else c.x=v;};
        set(d.previous[0],sign*.25);set(d.previous[horizontal?12:3],sign*.125);
        d.arrival=d.previous;for(auto& c:d.arrival)c=c*2;
        const Slab slab(oldTime,endTime,d);const Vector2D a=horizontal?Vector2D{.25,.5}:Vector2D{.5,.25},b{.5,.5};
        const double q0=.25,q1=.75;
        const double ordinary=.5*(.25*(q1-q0)+.125*(std::pow(q1,4)-std::pow(q0,4))/4);
        const double radialMoment=.5*(.25*(.125*(q1-q0)+.25*(q1*q1-q0*q0))
            +.125*(.125*(std::pow(q1,4)-std::pow(q0,4))/4+.5*(std::pow(q1,5)-std::pow(q0,5))/5));
        const double expected=duration*1.5*(radial?(horizontal?2*pi*radialMoment:pi*ordinary):ordinary);
        const auto bounds=slab.encloseNormalFaceIntegral(a,b,oldTime,endTime,radial,.25L,.25L);
        contains(bounds,sign>0?std::array<double,4>{expected,0,expected*.25,0}:std::array<double,4>{0,expected,0,expected*.25},2e-18);
        const auto half1=slab.encloseNormalFaceIntegral(a,b,oldTime,oldTime+duration*.5,radial,.25L,.25L);
        const auto half2=slab.encloseNormalFaceIntegral(a,b,oldTime+duration*.5,endTime,radial,.25L,.25L);
        contains(bounds.intersected(half1.added(half2)),sign>0?std::array<double,4>{expected,0,expected*.25,0}:std::array<double,4>{0,expected,0,expected*.25},2e-18);++coreCases;
        // Transverse space/time partitions are not nested. Both must still
        // contain the independent physical antiderivative after combining.
        const auto middle=a+(b-a)*.5;
        const auto spatial=slab.encloseNormalFaceIntegral(a,middle,oldTime,endTime,radial,.25L,.25L).added(
            slab.encloseNormalFaceIntegral(middle,b,oldTime,endTime,radial,.25L,.25L));
        const auto temporal=half1.added(half2),common=bounds.intersected(spatial).intersected(temporal);
        contains(common,sign>0?std::array<double,4>{expected,0,expected*.25,0}:std::array<double,4>{0,expected,0,expected*.25},2e-18);
        check(common.errorBound()<=spatial.errorBound() && common.errorBound()<=temporal.errorBound(),"common proof lost a partition certificate");
    }
    auto d=translation(false);d.arrival[0].x=-.25;const Slab reverse(oldTime,endTime,d);
    const auto b=reverse.encloseNormalFaceIntegral({.5,.25},{.5,.5},oldTime,endTime,false,1,1);
    contains(b,{.25*.25*duration/4,.25*.25*duration/4,.25*.25*duration/4,.25*.25*duration/4});
    check(b.volumes[0].upper>0 && b.volumes[1].upper>0,"gross counterflow lost to mean velocity");
    check(!Quantum::Geometry::polynomialFaceMidpoint(1,std::nextafter(1.,2.)),"one-ulp interval has false midpoint");
    rejects([&]{reverse.encloseNormalFaceIntegral({0,0},{1,1},oldTime,endTime,false);},"axis-aligned");
    rejects([&]{reverse.encloseNormalFaceIntegral({-.5,0},{-.5,1},oldTime,endTime,true);},"negative radius");
    rejects([&]{b.withContentRange(.5L,.4L);},"content");
    rejects([&]{b.withContentRange(0,1.01L);},"content");
    rejects([&]{reverse.encloseNormalFaceIntegral({.5,.25},{.5,.5},endTime,oldTime,false);},"time");
    std::cout<<"CORE_FACE_MOMENTS,PASS cases="<<coreCases<<'\n';
}
void translations(bool radial,bool amr)
{
    auto mesh=makeMesh(amr);const bool horizontal=radial;const auto f=findFace(*mesh,radial,horizontal);const auto area=measure(f,radial);
    const HermiteGalsSpaceTimeIntegralConfig config{area*.25*duration*1e-5,4096,40};
    for(int sign:{-1,1})for(double q:{.25,.75}){
        const auto old=stateFor(*mesh,plane(horizontal,.5-sign*q*.25*duration,sign),radial);
        const auto interval=intervalFor(*mesh,old,amr,translation(horizontal,sign*.25));
        const double t=area*.25*duration,l=t*(1-q);const auto expected=sign>0?std::array<double,4>{t,0,l,0}:std::array<double,4>{0,t,0,l};
        const auto result=integrateHermiteGalsPolynomialFace(*mesh,interval,f.row,config);verify(*mesh,result,expected,"translation");
        check(result.transfer().interval().sameInterval(interval),"integral provenance changed");
        const auto repeated=integrateHermiteGalsPolynomialFace(*mesh,interval,f.row,config);
        for(int k=0;k<4;++k)check(result.bounds().volumes[k].lower==repeated.bounds().volumes[k].lower && result.bounds().volumes[k].upper==repeated.bounds().volumes[k].upper,"nondeterministic bounds");
        check(result.report().evaluatedPanels==repeated.report().evaluatedPanels,"nondeterministic work");
        const Interval first(*mesh,old,interval.velocitySource(),interval.characteristicSource(),oldTime,oldTime+duration*.5),
            second(*mesh,old,interval.velocitySource(),interval.characteristicSource(),oldTime+duration*.5,endTime);
        const auto a=integrateHermiteGalsPolynomialFace(*mesh,first,f.row,config),b=integrateHermiteGalsPolynomialFace(*mesh,second,f.row,config);
        contains(a.bounds().added(b.bounds()),expected);const auto combined=a.transfer().appended(*mesh,b.transfer());
        check(combined.interval().sameInterval(interval) && std::abs(combined.signedLiquidVolume()-sign*l)<=a.report().absoluteErrorBound+b.report().absoluteErrorBound,"time composition failed");
        const auto loose=integrateHermiteGalsPolynomialFace(*mesh,interval,f.row,{config.absoluteVolumeTolerance*16,4096,40});
        check(result.report().absoluteErrorBound<=loose.report().absoluteErrorBound && result.report().acceptedPanels>=loose.report().acceptedPanels,"tighter budget failed to refine");
    }
    auto d=translation(horizontal);d.arrival[0]=d.arrival[0]*(-1);
    const auto reverseOld=stateFor(*mesh,plane(horizontal,.5-3*.25*duration/16),radial);
    const auto reverse=integrateHermiteGalsPolynomialFace(*mesh,intervalFor(*mesh,reverseOld,amr,d),f.row,config);
    const double gross=area*.25*duration/4;verify(*mesh,reverse,{gross,gross,gross/4,gross/4},"reversing_velocity");
    check(reverse.transfer().integratedVolumes().integratedLiquidPositiveVolume>0 && reverse.transfer().integratedVolumes().integratedLiquidNegativeVolume>0,"hidden reversing liquid lost");
    constexpr double epsilon=.000244140625;const double center=.5-.25*duration*.75;
    const Profile band=[=](Vector2D p){const double z=(horizontal?p.y:p.x)-center;return Jet{z*z-epsilon*epsilon,horizontal?0.:2*z,horizontal?2*z:0.,0.};};
    const auto bandOld=stateFor(*mesh,band,radial);for(double q:{0.,.5,1.})check(band(horizontal?Vector2D{.3,.5-.25*duration*q}:Vector2D{.5-.25*duration*q,.3}).value>0,"band sparse samples not gas");
    const auto thin=integrateHermiteGalsPolynomialFace(*mesh,intervalFor(*mesh,bandOld,amr,translation(horizontal)),f.row,config);
    verify(*mesh,thin,{area*.25*duration,0,area*2*epsilon,0},"hidden_thin_band");++thinEvents;
    for(int sign:{-1,1}){
        const double weight=sign<0?1-4e-11:4e-11;
        const auto old=stateFor(*mesh,pure(sign),radial,[=](int,Vector2D,Vector2D,double cap){return cap*weight;});
        const auto original=old.primary().phaseVolumes;const HermiteGalsPhaseRegionSource regions(*mesh,old);const HermiteGalsInventoryRegionSource inventory(*mesh,regions);
        const auto interval=intervalFor(*mesh,old,amr,translation(horizontal));
        const auto result=integrateHermiteGalsPolynomialFace(*mesh,interval,f.row,inventory,{1e-12,4096,40});
        verify(*mesh,result,{area*.25*duration,0,area*.25*duration*weight,0},"weighted_pure",2e-17);++weightedCases;
        check(result.inventorySource()!=nullptr,"lost retained inventory");
        for(std::size_t k=0;k<original.size();++k)check(original[k].liquid==old.primary().phaseVolumes[k].liquid,"integration changed M");
        const auto pureResult=integrateHermiteGalsPolynomialFace(*mesh,interval,f.row,{1e-12,4096,40});
        verify(*mesh,pureResult,{area*.25*duration,0,sign<0?area*.25*duration:0.,0},"pure");
        const auto zero=integrateHermiteGalsPolynomialFace(*mesh,intervalFor(*mesh,old,amr,Data{}),f.row);
        verify(*mesh,zero,{0,0,0,0},"zero");
    }
    if(radial){const auto old=stateFor(*mesh,pure(-1),true);const auto axis=findFace(*mesh,true,false,0,0);
        verify(*mesh,integrateHermiteGalsPolynomialFace(*mesh,intervalFor(*mesh,old,amr,translation(true)),axis.row),{0,0,0,0},"axis_zero");}
    std::cout<<"TRANSLATION_MATRIX,"<<radial<<','<<amr<<",PASS\n";
}
void curvedVelocity(bool amr)
{
    auto mesh=makeMesh(amr);const auto f=findFace(*mesh,false,false);const double a=f.a.y,b=f.b.y;
    Data shear;shear.previous[1].x=.5;shear.arrival=shear.previous;
    const auto old=stateFor(*mesh,plane(false,.5-.015625),false);
    const double total=.5*duration*(b*b-a*a)/2,liquid=total-.015625*(b-a);
    matrixCase(amr?"shear_amr":"shear_uniform",[&]{verify(*mesh,integrateHermiteGalsPolynomialFace(*mesh,intervalFor(*mesh,old,amr,shear),f.row,{total*1e-3,4096,40}),{total,0,liquid,0},"shear");});
    Data swirl;swirl.origin={.5,.5};swirl.previous[9].x=-.5;swirl.previous[3].x=-.5;swirl.previous[12].y=.5;swirl.previous[6].y=.5;swirl.arrival=swirl.previous;
    const auto full=stateFor(*mesh,pure(-1),false);const double flux=.5*duration*(std::pow(a-.5,4)-std::pow(b-.5,4))/4;
    verify(*mesh,integrateHermiteGalsPolynomialFace(*mesh,intervalFor(*mesh,full,amr,swirl),f.row,{1e-12,4096,40}),{flux,0,flux,0},"nonlinear_pure_swirl");
    if(amr)return;
    constexpr double radius=.1875,omega=.5;
    const Profile circle=[=](Vector2D p){const double x=p.x-.5,y=p.y-.5;return Jet{x*x+y*y-radius*radius,2*x,2*y,0.};};
    const auto cut=stateFor(*mesh,circle,false);Data rotation;rotation.origin={.5,.5};rotation.previous[1].x=-omega;rotation.previous[4].y=omega;rotation.arrival=rotation.previous;
    double integral=0;const double lambda=std::pow(omega,4)/4;
    for(int k=0;k<6;++k)integral+=std::pow(-lambda,k)*std::pow(duration,4*k+1)/(4*k+1);
    const double ref=omega*radius*radius*.5*integral;
    matrixCase("rotating_cut_circle",[&]{verify(*mesh,integrateHermiteGalsPolynomialFace(*mesh,intervalFor(*mesh,cut,false,rotation),f.row,{2e-6,4096,40}),{omega*.25*.25*duration/2,0,ref,0},"rotating_cut_circle",1e-17);});
    // For this original discrete RK2 swirl, radius departure squared is
    // y^2*(1+a^2/2+3a^3/16+a^4/64), a=k^2*d^2*y^4.
    // Monotone root brackets + monotone time Riemann sums bound the reference.
    const auto root=[=](double t){double lo=0,hi=radius;
        for(int n=0;n<60;++n){const double y=(lo+hi)*.5,a=.25*t*t*std::pow(y,4);
            const double r2=y*y*(1+.5*a*a+.1875*a*a*a+.015625*a*a*a*a);
            if(r2<radius*radius)lo=y;else hi=y;}return std::array<double,2>{lo,hi};};
    double lower=0,upper=0;constexpr int n=1024;
    for(int j=0;j<n;++j){lower+=std::pow(root(duration*(j+1)/n)[0],4);upper+=std::pow(root(duration*j/n)[1],4);}
    lower*=.5/4*duration/n;upper*=.5/4*duration/n;
    matrixCase("nonlinear_cut_swirl",[&]{verify(*mesh,integrateHermiteGalsPolynomialFace(*mesh,intervalFor(*mesh,cut,false,swirl),f.row,{2e-7,4096,40}),{flux,0,(lower+upper)*.5,0},"nonlinear_cut_swirl",(upper-lower)*.5+1e-17);});
}
void guardsAndLifetime()
{
    auto mesh=makeMesh(false);const auto old=stateFor(*mesh,plane(false,.5-.25*duration*.25),false);const auto f=findFace(*mesh,false,false);
    auto interval=intervalFor(*mesh,old,false,translation(false));const auto before=old.primary().phaseVolumes;const auto* scalar=&old.geometry().snapshot();
    for(const auto config:std::array<HermiteGalsSpaceTimeIntegralConfig,4>{{{0,4096,40},{std::numeric_limits<double>::quiet_NaN(),4096,40},{1e-6,0,40},{1e-6,4096,0}}})
        rejects([&]{integrateHermiteGalsPolynomialFace(*mesh,interval,f.row,config);},"budget");
    // A translated plane now has an exact affine subset measure. Exercise
    // exhaustion with the original bilinear shear map, whose remainder is
    // nonzero, rather than requiring a resolved plane to fail gratuitously.
    Data shear;shear.previous[1].x=.5;shear.arrival=shear.previous;
    const auto unresolved=intervalFor(*mesh,old,false,shear);
    rejects([&]{integrateHermiteGalsPolynomialFace(*mesh,unresolved,f.row,{1e-10,1,40});},"panel budget exhausted");
    rejects([&]{integrateHermiteGalsPolynomialFace(*mesh,unresolved,f.row,{1e-10,4096,1});},"depth budget exhausted");
    const auto boundary=findFace(*mesh,false,false,0,.25);
    rejects([&]{integrateHermiteGalsPolynomialFace(*mesh,interval,boundary.row,{1,4,40});},"pending=");
    HermiteGalsCandidate next(*mesh,{3,oldTime,endTime,Position::End},old.primary());
    const auto constant=Velocity::constant(old,next.identity(),next.timeLayer(),{.25,0});
    const Interval constantInterval(*mesh,old,constant,interval.characteristicSource(),oldTime,endTime);
    rejects([&]{integrateHermiteGalsPolynomialFace(*mesh,constantInterval,f.row);},"owned polynomial");
    const auto opaque=Velocity(old,next.identity(),next.timeLayer(),[](Vector2D){return LevelSetBoundaryTrace::VelocityJet{{.25,0},0,0,0,0};},[](Vector2D){return LevelSetBoundaryTrace::VelocityJet{{.25,0},0,0,0,0};});
    const Interval opaqueInterval(*mesh,old,opaque,interval.characteristicSource(),oldTime,endTime);
    rejects([&]{integrateHermiteGalsPolynomialFace(*mesh,opaqueInterval,f.row);},"owned polynomial");
    auto domain=interval.characteristicSource().domain();domain.activeSegment=[](Vector2D,Vector2D){return true;};
    const Characteristic opaqueDomain(*mesh,domain);const Interval unsupported(*mesh,old,interval.velocitySource(),opaqueDomain,oldTime,endTime);
    rejects([&]{integrateHermiteGalsPolynomialFace(*mesh,unsupported,f.row);},"opaque");
    const auto foreignOld=stateFor(*mesh,pure(-1),false,[](int,Vector2D,Vector2D,double cap){return cap;});
    const HermiteGalsPhaseRegionSource regions(*mesh,foreignOld);auto inventory=HermiteGalsInventoryRegionSource(*mesh,regions);
    rejects([&]{integrateHermiteGalsPolynomialFace(*mesh,interval,f.row,inventory);},"foreign inventory");
    const auto keptInventory=std::move(inventory);rejects([&]{integrateHermiteGalsPolynomialFace(*mesh,interval,f.row,inventory);},"moved-from");
    {struct Restore{int previous=std::fegetround();~Restore(){std::fesetround(previous);}} restore;
        check(std::fesetround(FE_UPWARD)==0,"cannot change rounding");rejects([&]{integrateHermiteGalsPolynomialFace(*mesh,interval,f.row);},"round-to-nearest");}
#if defined(_M_X64) || (defined(__x86_64__) && defined(__SSE2__))
    {struct Restore{unsigned previous=_mm_getcsr();~Restore(){_mm_setcsr(previous);}} restore;
        _mm_setcsr(restore.previous|(1u<<15));rejects([&]{integrateHermiteGalsPolynomialFace(*mesh,interval,f.row);},"gradual underflow");}
#endif
    const auto radialOld=stateFor(*mesh,pure(-1),true);const auto radialInterval=intervalFor(*mesh,radialOld,false,translation(false));
    const auto axis=findFace(*mesh,true,false,0,0);rejects([&]{integrateHermiteGalsPolynomialFace(*mesh,radialInterval,axis.row);},"axis");
    std::optional<Result> lifetime;
    {const auto localOld=stateFor(*mesh,pure(-1),false);lifetime=integrateHermiteGalsPolynomialFace(*mesh,intervalFor(*mesh,localOld,false,translation(false)),f.row);}
    verify(*mesh,*lifetime,{.25*.25*duration,0,.25*.25*duration,0},"retained_owner");
    auto copied=*lifetime,moved=std::move(copied);rejects([&]{copied.bounds();},"moved-from bounded");check(&moved.bounds()==&lifetime->bounds(),"copy lost immutable owner");
    auto different=makeMesh(true);rejects([&]{integrateHermiteGalsPolynomialFace(*different,interval,f.row);},"mesh");
    const auto retainedInterval=std::move(interval);rejects([&]{integrateHermiteGalsPolynomialFace(*mesh,interval,f.row);},"moved-from");
    check(scalar==&old.geometry().snapshot(),"failure changed geometry");for(std::size_t k=0;k<before.size();++k)check(before[k].liquid==old.primary().phaseVolumes[k].liquid,"failure changed M");
    DictionaryMesh::RefineBox box;box.level=1;box.min={.25,.25};box.max={.75,.75};mesh->applyRefineBoxes({box});rejects([&]{lifetime->validateMesh(*mesh);},"mesh");
    std::cout<<"INTEGRAL_GUARDS_AND_LIFETIME,PASS checks="<<guards<<'\n';
}
void tightTimeTraces()
{
    const LevelSetBoundaryTrace::Domain domain{{0,0},{1,1},{.25,.25},false};
    for(int kind=0;kind<6;++kind){Data d;
        if(kind<3)d.previous[0]={kind==2?1e-310:.25,0};
        else {d.origin={.5,.5};d.previous[kind==5?3:1].x=-.5;d.previous[kind==5?12:4].y=.5;}
        d.arrival=d.previous;if(kind==1 || kind==4)for(auto& c:d.arrival)c=c*(-1);
        const Slab slab(oldTime,endTime,d);const Vector2D lo{.4,.4},hi{.45,.45};
        const auto bounds=slab.encloseInteriorRk2TightTime(domain,lo,hi,oldTime,endTime);
        for(int ti=1;ti<=16;++ti)for(int xi=0;xi<3;++xi)for(int yi=0;yi<3;++yi){const double t=oldTime+duration*ti/16;
            const Vector2D p{lo.x+(hi.x-lo.x)*xi/2,lo.y+(hi.y-lo.y)*yi/2};
            const auto trace=LevelSetBoundaryTrace::traceWithJacobian(domain,p,t-oldTime,
                [&](Vector2D x){return slab.evaluate(x,t);},[&](Vector2D x){return slab.evaluate(x,oldTime);});
            check(!trace.trajectory.boundaryHit && trace.acceptedSubsteps==1 && bool(trace.derivative),"tight range changed trace branch");
            const auto b=trace.trajectory.point;const auto j=*trace.derivative;const std::array<double,4> jac{j.xx,j.xy,j.yx,j.yy};
            check(b.x>=bounds.departureLower.x && b.x<=bounds.departureUpper.x && b.y>=bounds.departureLower.y && b.y<=bounds.departureUpper.y,"actual B escaped tight time range");
            for(int k=0;k<4;++k)check(bounds.jacobian[k].lower<=jac[k] && jac[k]<=bounds.jacobian[k].upper,"actual J escaped tight time range");++tightTraces;}
    }
    std::cout<<"TIGHT_TIME_ACTUAL_TRACES,PASS traces="<<tightTraces<<'\n';
}
void boundaryHistoryContract()
{
    using History=HermiteGalsBoundaryPhaseSource;using P=HermiteGalsInflowPhase;
    const std::array<Boundary,4> sides{Boundary::west,Boundary::south,Boundary::east,Boundary::north};
    const std::array<P,4> gas{P::Gas,P::Gas,P::Gas,P::Gas};
    int cases=0,actualTraces=0;const int oldGuards=guards;
    for(bool radial:{false,true})for(bool amr:{false,true})for(int side=0;side<4;++side){
        if(radial && side==0)continue; // the symmetry axis is not an inflow boundary
        const bool horizontal=side==1 || side==3;const double sign=side<2?1:-1;
        auto mesh=makeMesh(amr);const auto old=stateFor(*mesh,pure(1),radial,[](int,Vector2D,Vector2D,double c){return .25*c;});
        const auto before=old.primary().phaseVolumes;Data data;
        if(horizontal){data.previous[0].y=sign*.125;data.previous[4].y=sign*.0625;}
        else {data.previous[0].x=sign*.125;data.previous[1].x=sign*.0625;}
        data.arrival=data.previous;for(auto& c:data.arrival)c=c*1.5;
        const auto base=intervalFor(*mesh,old,amr,data);const auto face=findFace(*mesh,radial,horizontal,side<2?0:1,.25);
        const auto* field=base.velocitySource().polynomialVelocity();const auto& domain=base.characteristicSource().domain();
        check(field->encloseInitialBoundaryHit(domain,face.a,face.b,oldTime,endTime)==sides[side],"wrong certified boundary side");
        for(int ti=0;ti<=8;++ti)for(int si=0;si<=8;++si){const double t=oldTime+duration*ti/8;const auto p=face.a+(face.b-face.a)*(si/8.);
            const auto trace=LevelSetBoundaryTrace::traceWithJacobian(domain,p,t-oldTime,
                [&](Vector2D x){return field->evaluate(x,t);},[&](Vector2D x){return field->evaluate(x,oldTime);});
            check(trace.trajectory.point.x==p.x && trace.trajectory.point.y==p.y && trace.trajectory.elapsed==0
                && trace.acceptedSubsteps==0,"initial event changed point/time/substeps");
            if(ti)check(trace.trajectory.boundaryHit && trace.trajectory.side==sides[side] && !trace.derivative,"original initial event mismatch");
            else check(!trace.trajectory.boundaryHit && bool(trace.derivative),"zero duration identity changed");++actualTraces;
        }
        const HermiteGalsPhaseRegionSource regions(*mesh,old);const HermiteGalsInventoryRegionSource inventory(*mesh,regions,{.3});
        for(P phase:{P::Gas,P::Liquid}){
            auto phases=gas;phases[side]=phase;const History history(*mesh,base,phases);
            const Interval bound(*mesh,old,base.velocitySource(),base.characteristicSource(),oldTime,endTime,history);
            check(bound.boundaryPhaseSource()->sameSource(history) && !bound.inflowSource() && !bound.velocitySource().constantValue(),"boundary capability/source lost");
            const auto result=integrateHermiteGalsPolynomialFace(*mesh,bound,face.row,inventory,{1e-10,64,40});
            const double a=horizontal?face.a.x:face.a.y,b=horizontal?face.b.x:face.b.y;
            // Independent primitive of (.125+.0625*s)*(1+.5*theta), with true RZ measure.
            const double spatial=radial && horizontal?2*pi*(.125*(b*b-a*a)/2+.0625*(b*b*b-a*a*a)/3):
                (.125*(b-a)+.0625*(b*b-a*a)/2)*(radial?2*pi*face.a.x:1);
            const double total=spatial*duration*1.25;std::array<double,4> expected{};const int k=sign>0?0:1;
            expected[k]=total;expected[k+2]=phase==P::Liquid?total:0;contains(result.bounds(),expected,1e-17);
            check(result.report().initialBoundaryQueries==1 && result.report().phaseQueries==0
                && result.report().acceptedPanels==1 && result.transfer().interval().sameInterval(bound),"initial hit was not certified directly");
            const double split=oldTime+duration*.5;
            const Interval first(*mesh,old,base.velocitySource(),base.characteristicSource(),oldTime,split,history),
                second(*mesh,old,base.velocitySource(),base.characteristicSource(),split,endTime,history);
            const auto left=integrateHermiteGalsPolynomialFace(*mesh,first,face.row,inventory),right=integrateHermiteGalsPolynomialFace(*mesh,second,face.row,inventory);
            const auto combined=left.transfer().appended(*mesh,right.transfer());check(combined.interval().sameInterval(bound),"append dropped boundary history");
            check(std::abs(double(combined.signedLiquidVolume()-result.transfer().signedLiquidVolume()))<=1e-14,"split boundary integral changed");
            const History independent(*mesh,base,phases);const Interval foreign(*mesh,old,base.velocitySource(),base.characteristicSource(),split,endTime,independent);
            check(!second.sameInterval(foreign),"equal-valued histories merged");
            const HermiteGalsIntegratedFaceTransfer foreignTransfer(*mesh,foreign,face.row,right.transfer().orientation(),right.transfer().integratedVolumes());
            rejects([&]{left.transfer().appended(*mesh,foreignTransfer);},"different boundary history");++cases;
        }
        // A declared gas wall must not erase nonzero old gas inventory inside.
        const History history(*mesh,base,gas);const Interval bound(*mesh,old,base.velocitySource(),base.characteristicSource(),oldTime,endTime,history);
        const auto interior=findFace(*mesh,radial,horizontal,.5,.25);
        const auto result=integrateHermiteGalsPolynomialFace(*mesh,bound,interior.row,inventory);
        const auto v=result.bounds().representative();check(std::abs(v[2]+v[3]-.25*(v[0]+v[1]))<1e-14 && v[2]+v[3]>0
            && result.report().initialBoundaryQueries==0,"boundary gas replaced interior inventory content");
        for(std::size_t k=0;k<before.size();++k)check(before[k].liquid==old.primary().phaseVolumes[k].liquid,"history/integral changed old M");
    }
    auto mesh=makeMesh(false);const auto old=stateFor(*mesh,pure(1),false);const auto base=intervalFor(*mesh,old,false,translation(false));
    const History history(*mesh,base,gas);const Interval bound(*mesh,old,base.velocitySource(),base.characteristicSource(),oldTime,endTime,history);
    const auto face=findFace(*mesh,false,false,0,0);const auto& domain=base.characteristicSource().domain();
    rejects([&]{History(*mesh,bound,gas);},"unbound");
    rejects([&]{auto bad=gas;bad[0]=static_cast<P>(7);History(*mesh,base,bad);},"invalid boundary phase");
    rejects([&]{history.phaseFor(Boundary::southWest);},"invalid boundary side");
    const auto another=intervalFor(*mesh,old,false,translation(false));rejects([&]{history.validateInterval(*mesh,another);},"different velocity source");
    const auto otherOld=stateFor(*mesh,pure(1),false);const auto other=intervalFor(*mesh,otherOld,false,translation(false));
    rejects([&]{history.validateInterval(*mesh,other);},"different old source");
    const Characteristic otherMap(*mesh,domain);const Interval otherCharacteristic(*mesh,old,base.velocitySource(),otherMap,oldTime,endTime);
    rejects([&]{history.validateInterval(*mesh,otherCharacteristic);},"different characteristic source");
    const Interval shortExtent(*mesh,old,base.velocitySource(),base.characteristicSource(),oldTime,oldTime+duration*.5);
    const History shortHistory(*mesh,shortExtent,gas);rejects([&]{shortHistory.validateInterval(*mesh,base);},"does not cover interval");
    HermiteGalsCandidate next(*mesh,{3,oldTime,endTime,Position::End},old.primary());
    const auto constant=Velocity::constant(old,next.identity(),next.timeLayer(),{.25,0});
    const Interval constantInterval(*mesh,old,constant,base.characteristicSource(),oldTime,endTime);
    rejects([&]{History(*mesh,constantInterval,gas);},"owned polynomial");
    auto copy=history;auto moved=std::move(copy);check(moved.sameSource(history),"history copy changed owner");rejects([&]{copy.phaseFor(Boundary::west);},"moved-from");
    auto differentMesh=makeMesh(true);rejects([&]{history.validateInterval(*differentMesh,base);},"mesh");
    rejects([&]{integrateHermiteGalsPolynomialFace(*mesh,base,face.row,{1,4,40});},"pending=");
    const auto guard=[&](Data data,LevelSetBoundaryTrace::Domain d,Vector2D a,Vector2D b,const char* expected){
        const Slab slab(oldTime,endTime,data);rejects([&]{slab.encloseInitialBoundaryHit(d,a,b,oldTime,endTime);},expected);};
    guard(translation(true),domain,face.a,face.b,"unique immediate hit");
    Data corner=translation(false);corner.previous[0].y=.125;corner.arrival=corner.previous;
    guard(corner,domain,face.a,face.b,"unique immediate hit");
    guard(translation(false,2),domain,face.a,face.b,"half-cell");
    auto opaque=domain;opaque.activeSegment=[](Vector2D,Vector2D){return true;};guard(translation(false),opaque,face.a,face.b,"opaque");
    auto radial=domain;radial.radialAxisSymmetry=true;guard(translation(false),radial,face.a,face.b,"axis");
    Data reversal=translation(false);reversal.arrival[0].x=-.25;guard(reversal,domain,face.a,face.b,"unique immediate hit");
    guard(translation(false),domain,{.5,.25},{.5,.5},"unique immediate hit");
    const Slab slab(oldTime,endTime,translation(false));
    rejects([&]{slab.encloseInitialBoundaryHit(domain,face.a,face.b,oldTime,oldTime);},"positive duration");
    // Complete ordered algebraic transfers must retain this exact history.
    FaceSegmentGraph graph;graph.build(*mesh);const auto native=NativeLeafDomain::build(*mesh,graph,false);
    std::vector<HermiteGalsIntegratedFaceTransfer> zeros;
    for(const auto f:native.faces(*mesh))zeros.emplace_back(*mesh,bound,f.nativeRow,f.horizontal?HermiteGalsFaceTransferOrientation::PositiveY:HermiteGalsFaceTransferOrientation::PositiveX,HermiteGalsIntegratedFaceVolumes{});
    const auto assembly=assembleHermiteGalsIntegratedInventory(*mesh,bound,zeros);check(assembly.interval().sameInterval(bound),"assembly lost boundary owner");
    for(const auto& f:assembly.faces())check(f.interval().sameInterval(bound),"assembly transfer lost boundary owner");
    zeros[0]=HermiteGalsIntegratedFaceTransfer(*mesh,base,zeros[0].face().nativeRow,zeros[0].orientation(),{});
    rejects([&]{assembleHermiteGalsIntegratedInventory(*mesh,bound,zeros);},"different interval source");
    check(cases==28 && actualTraces==1134,"boundary matrix incomplete");
    std::cout<<"I12a16r1 boundary history PASS cases="<<cases<<" actual_traces="<<actualTraces<<" guards="<<guards-oldGuards<<'\n';
}
void monotoneWallContract()
{
    int boxes=0,actualTraces=0;const int oldGuards=guards;
    for(bool radial:{false,true})for(int side=0;side<4;++side){
        if(radial && side!=1 && side!=3)continue;
        const bool horizontal=side==1 || side==3;const double sign=side<2?-1:1;
        Vector2D lo=horizontal?Vector2D{radial?0:.25,side<2?0:.875}:Vector2D{side<2?0:.875,.25};
        Vector2D hi=lo+(horizontal?Vector2D{.25,.125}:Vector2D{.125,.25});
        const LevelSetBoundaryTrace::Domain domain{{0,0},{1,1},{.25,.25},radial};
        for(int kind=0;kind<3;++kind){
            Data data;const double speed=kind==2?1e-310:.125;
            if(horizontal){data.previous[0].y=sign*speed;if(kind==1)data.previous[1].y=sign*.0625;}
            else {data.previous[0].x=sign*speed;if(kind==1)data.previous[4].x=sign*.0625;}
            data.arrival=data.previous;if(kind==1)for(auto& c:data.arrival)c=c*1.25;
            const Slab slab(oldTime,endTime,data);const auto bounds=slab.encloseInteriorRk2TightTime(domain,lo,hi,oldTime,endTime);
            for(int ti=0;ti<=8;++ti)for(int xi=0;xi<=4;++xi)for(int yi=0;yi<=4;++yi){
                const double t=oldTime+duration*ti/8;const Vector2D p{lo.x+(hi.x-lo.x)*xi/4,lo.y+(hi.y-lo.y)*yi/4};
                const auto trace=LevelSetBoundaryTrace::traceWithJacobian(domain,p,t-oldTime,
                    [&](Vector2D x){return slab.evaluate(x,t);},[&](Vector2D x){return slab.evaluate(x,oldTime);});
                check(!trace.trajectory.boundaryHit && trace.derivative && trace.acceptedSubsteps==(ti?1:0),"monotone box missed original branch");
                const auto b=trace.trajectory.point;const auto j=*trace.derivative;const std::array<double,4> jac{j.xx,j.xy,j.yx,j.yy};
                check(b.x>=bounds.departureLower.x && b.x<=bounds.departureUpper.x && b.y>=bounds.departureLower.y
                    && b.y<=bounds.departureUpper.y,"monotone B escaped bound");
                for(int k=0;k<4;++k)check(bounds.jacobian[k].lower<=jac[k] && jac[k]<=bounds.jacobian[k].upper,"monotone J escaped bound");
                const double from=horizontal?p.y:p.x,to=horizontal?b.y:b.x;
                check(side<2?to>=from:to<=from,"actual floating target violated coordinate order");++actualTraces;
            }
            ++boxes;
        }
    }
    const LevelSetBoundaryTrace::Domain domain{{0,0},{1,1},{.25,.25},false};
    for(bool horizontal:{false,true})for(int sign:{-1,1}){
        const Slab slab(oldTime,endTime,translation(horizontal,sign*.25));
        const double normal=sign>0?.03125:.96875;const Vector2D p=horizontal?Vector2D{.5,normal}:Vector2D{normal,.5};
        const auto actual=LevelSetBoundaryTrace::traceWithJacobian(domain,p,duration,
            [&](Vector2D x){return slab.evaluate(x,endTime);},[&](Vector2D x){return slab.evaluate(x,oldTime);});
        check(actual.trajectory.boundaryHit && actual.trajectory.elapsed==duration,"later wall event fixture missing");
        rejects([&]{slab.encloseInteriorRk2TightTime(domain,p,p,oldTime,endTime);},"unresolved");
        const Vector2D wall=horizontal?Vector2D{.5,sign>0?0.:1.}:Vector2D{sign>0?0.:1.,.5};
        rejects([&]{slab.encloseInteriorRk2TightTime(domain,wall,wall,oldTime,endTime);},"initial boundary hit");
    }
    check(boxes==18 && actualTraces==4050 && guards-oldGuards==8,"monotone wall matrix incomplete");
    std::cout<<"I12a16r2 monotone wall PASS boxes="<<boxes<<" actual_traces="<<actualTraces<<" guards="<<guards-oldGuards<<'\n';
}
void tangentZeroContract()
{
    int cases=0,actualTraces=0;
    for(bool amr:{false,true})for(bool cubic:{false,true})for(bool horizontal:{false,true})for(bool upper:{false,true}){
        auto mesh=makeMesh(amr);const auto old=stateFor(*mesh,pure(1),false,[](int,Vector2D,Vector2D,double c){return .25*c;});
        Data data;data.origin={.5,.5};
        if(cubic){data.previous[3].x=-.5;data.previous[9].x=-.5;data.previous[12].y=.5;data.previous[6].y=.5;}
        else {data.previous[1].x=-.5;data.previous[4].y=.5;}data.arrival=data.previous;
        const auto interval=intervalFor(*mesh,old,amr,data);const auto face=findFace(*mesh,false,horizontal,.5,upper?.75:0);
        const auto* field=interval.velocitySource().polynomialVelocity();const auto& domain=interval.characteristicSource().domain();
        const auto bounds=field->encloseInteriorRk2TightTime(domain,face.a,face.b,oldTime,endTime);
        for(int ti=0;ti<=8;++ti)for(int si=0;si<=8;++si){const double t=oldTime+duration*ti/8;
            const auto p=face.a+(face.b-face.a)*(si/8.);const auto trace=LevelSetBoundaryTrace::traceWithJacobian(domain,p,t-oldTime,
                [&](Vector2D q){return field->evaluate(q,t);},[&](Vector2D q){return field->evaluate(q,oldTime);});
            const auto b=trace.trajectory.point;
            check(!trace.trajectory.boundaryHit && trace.derivative && trace.acceptedSubsteps==(ti?1:0)
                && b.x>=bounds.departureLower.x && b.x<=bounds.departureUpper.x
                && b.y>=bounds.departureLower.y && b.y<=bounds.departureUpper.y,"tangent-zero original trajectory escaped no-hit certificate");++actualTraces;
        }
        const HermiteGalsPhaseRegionSource regions(*mesh,old);const HermiteGalsInventoryRegionSource inventory(*mesh,regions,{.3});
        const auto result=integrateHermiteGalsPolynomialFace(*mesh,interval,face.row,inventory,{1e-10,4096,40});
        const double a=(horizontal?face.a.x:face.a.y)-.5,b=(horizontal?face.b.x:face.b.y)-.5;
        const int power=cubic?4:2;const double signedVolume=(horizontal?1:-1)*.5*duration*(std::pow(b,power)-std::pow(a,power))/power;
        std::array<double,4> expected{};const int k=signedVolume>0?0:1;expected[k]=std::abs(signedVolume);expected[k+2]=expected[k]*.25;
        contains(result.bounds(),expected,1e-17);check(result.report().branchRefusals==0 && result.report().acceptedPanels==1
            && result.report().initialBoundaryQueries==0,"tangent-zero face required unproved boundary history");++cases;
    }
    check(cases==16 && actualTraces==1296,"tangent-zero matrix incomplete");
    std::cout<<"I12a16r3 tangent zero PASS faces="<<cases<<" actual_traces="<<actualTraces<<'\n';
    // Diagnose the NEXT original driver face without claiming its integral.
    // This is the actual n=8 face x=.625, y=[0,.125], same full endpoint
    // coefficients and half arrival time as the rotation/vortex benchmark.
    int hits=0,interior=0,nearWallEvents=0;
    for(bool cubic:{false,true}){
        const double end=cubic?.001:.01,arrival=end*.5,k=cubic?40:2;
        Data data;data.origin={.5,.5};
        if(cubic){data.previous[3].x=-k;data.previous[9].x=-k;data.previous[12].y=k;data.previous[6].y=k;}
        else {data.previous[1].x=-k;data.previous[4].y=k;}
        data.arrival=data.previous;if(cubic)for(auto& c:data.arrival)c=c*std::cos(pi*end/.4);
        const Slab slab(0,end,data);const LevelSetBoundaryTrace::Domain domain{{0,0},{1,1},{.0625,.0625},false};
        for(double y:{0.,1e-8,.0001,.001,.0625,.125}){
            const auto trace=LevelSetBoundaryTrace::traceWithJacobian(domain,{.625,y},arrival,
                [&](Vector2D p){return slab.evaluate(p,arrival);},[&](Vector2D p){return slab.evaluate(p,0);});
            if(trace.trajectory.boundaryHit){check(trace.trajectory.side==Boundary::south,"driver face hit unexpected wall");++hits;}
            else {check(bool(trace.derivative),"interior driver trace lost derivative");++interior;}
        }
        const auto near=LevelSetBoundaryTrace::traceWithJacobian(domain,{.625,1e-100},arrival,
            [&](Vector2D p){check(p.x>=0 && p.x<=1 && p.y>=0 && p.y<=1,"near-wall exterior polynomial sample");return slab.evaluate(p,arrival);},
            [&](Vector2D p){check(p.x>=0 && p.x<=1 && p.y>=0 && p.y<=1,"near-wall exterior polynomial sample");return slab.evaluate(p,0);});
        check(near.trajectory.boundaryHit && near.trajectory.roundoffMidpointEvent && near.trajectory.side==Boundary::south
            && near.trajectory.point.y==0 && near.trajectory.elapsed>=0
            && near.trajectory.elapsed<=128*std::numeric_limits<double>::epsilon()*arrival && !near.derivative,
            "original polynomial near-wall shortening failure unresolved");++nearWallEvents;
    }
    check(hits==7 && interior==5 && nearWallEvents==2,"original mixed boundary branch diagnosis incomplete");
    std::cout<<"I12a16r4 next face DIAGNOSTIC boundary_hits="<<hits<<" interior="<<interior<<" near_wall_events="<<nearWallEvents<<'\n';
}
void uniformReachableContract()
{
    using History=HermiteGalsBoundaryPhaseSource;using P=HermiteGalsInflowPhase;
    using Reach=Quantum::Geometry::PolynomialRk2Reachability;
    static_assert(!std::is_default_constructible_v<Reach>);
    const std::array<P,4> gas{P::Gas,P::Gas,P::Gas,P::Gas};
    int cases=0,actual=0,hits=0,interior=0,mapRefusals=0,sourceBoxes=0;
    int initialMulti=0,bracketMulti=0,tangent=0;const int oldGuards=guards;
    for(bool radial:{false,true})for(bool amr:{false,true})for(int kind=0;kind<4;++kind){
        if(radial && kind<2)continue;
        auto mesh=makeMesh(amr);Data data;
        if(kind<2){data.origin={.5,.5};
            if(kind==0){data.previous[1].x=-.25;data.previous[4].y=.25;}
            else {data.previous[3].x=-.5;data.previous[9].x=-.5;data.previous[12].y=.5;data.previous[6].y=.5;}}
        else data.previous[0]=radial?Vector2D{0,kind==2?.125:-.125}:Vector2D{.125,.125};
        data.arrival=data.previous;if(kind==2)for(auto& c:data.arrival)c=c*(-1);
        const auto face=findFace(*mesh,radial,radial,radial || kind>=2?0:.75,0);
        for(bool liquid:{false,true}){
            const auto old=stateFor(*mesh,pure(liquid?-1:1),radial,[=](int,Vector2D,Vector2D,double c){return liquid?c:0.;});
            const auto before=old.primary().phaseVolumes;const auto base=intervalFor(*mesh,old,amr,data);
            auto phases=gas;if(liquid)phases.fill(P::Liquid);const History history(*mesh,base,phases);
            const Interval bound(*mesh,old,base.velocitySource(),base.characteristicSource(),oldTime,endTime,history);
            const HermiteGalsPhaseRegionSource regions(*mesh,old);const HermiteGalsInventoryRegionSource inventory(*mesh,regions,{1e-6});
            const auto result=integrateHermiteGalsPolynomialFace(*mesh,bound,face.row,inventory,{1e-10,256,40});
            std::array<double,4> expected{};
            if(kind<2){const double a=.5-face.a.y,b=.5-face.b.y;
                expected[0]=duration*(kind==0?.125*(a*a-b*b):.5*((std::pow(a,4)-std::pow(b,4))/4+.0625*(a*a-b*b)/2));}
            else if(kind==2)expected[0]=expected[1]=.125*duration*measure(face,radial)/4;
            else expected[radial?1:0]=.125*duration*measure(face,radial);
            if(liquid){expected[2]=expected[0];expected[3]=expected[1];}contains(result.bounds(),expected,2e-17);
            check(result.report().uniformReachableQueries>0 && result.report().branchRefusals==0
                && result.transfer().interval().sameInterval(bound),"uniform content not certified with original owners");
            const auto* field=base.velocitySource().polynomialVelocity();const auto& domain=base.characteristicSource().domain();
            const auto support=field->encloseRk2Reachability(domain,endTime).cover(face.a,face.b);
            if(radial)check(support.lower.x==face.a.x && support.upper.x==face.b.x
                && !support.boundary[0] && !support.boundary[2],"zero radial component acquired artificial wall source");
            if(!liquid)for(int ti=1;ti<=8;++ti)for(int si=0;si<=8;++si){
                const double t=oldTime+duration*ti/8;const auto point=face.a+(face.b-face.a)*(si/8.);
                std::optional<LevelSetBoundaryTrace::JetResult> trace;
                try{trace=LevelSetBoundaryTrace::traceWithJacobian(domain,point,t-oldTime,
                    [&](Vector2D p){check(p.x>=0 && p.x<=1 && p.y>=0 && p.y<=1,"exterior callback");return field->evaluate(p,t);},
                    [&](Vector2D p){check(p.x>=0 && p.x<=1 && p.y>=0 && p.y<=1,"exterior callback");return field->evaluate(p,oldTime);});}
                catch(const std::runtime_error& e){const std::string reason=e.what();
                    if(reason.find("simultaneous backward crossings")!=std::string::npos)++initialMulti;
                    else if(reason.find("multiple wall crossings cannot be separated")!=std::string::npos)++bracketMulti;
                    else if(reason.find("boundary event is tangential")!=std::string::npos)++tangent;
                    else throw std::runtime_error("unexpected actual map refusal: "+reason);
                    ++mapRefusals;continue;}
                const auto p=trace->trajectory.point;check(p.x>=support.lower.x && p.x<=support.upper.x
                    && p.y>=support.lower.y && p.y<=support.upper.y,"actual source escaped reachable box");
                if(trace->trajectory.boundaryHit){++hits;check(!trace->derivative,"boundary gained invented J");
                    const auto side=*trace->trajectory.side;const int k=side==Boundary::west?0:side==Boundary::south?1:side==Boundary::east?2:3;
                    check(support.boundary[k] && t-trace->trajectory.elapsed>=oldTime,"boundary support/history missed source");}
                else {++interior;check(bool(trace->derivative),"actual interior J lost");}++actual;
            }
            for(std::size_t k=0;k<before.size();++k)check(before[k].liquid==old.primary().phaseVolumes[k].liquid,"constant algebra changed M");
            ++cases;
        }
    }
    // Same mixed interior/wall geometry, but different physical contents MUST
    // remain unresolved even with an enormous integration-error tolerance.
    auto mesh=makeMesh(false);Data rotation;rotation.origin={.5,.5};rotation.previous[1].x=-.25;rotation.previous[4].y=.25;rotation.arrival=rotation.previous;
    const auto face=findFace(*mesh,false,false,.75,0);
    for(int variant=0;variant<6;++variant){
        const Profile hidden=[](Vector2D p){const double x=p.x-.75,y=p.y-.03125;
            return Jet{x*x+y*y-.015625*.015625,2*x,2*y,0};};
        const auto profile=variant==2?hidden:pure(1);
        const auto old=stateFor(*mesh,profile,false,[=](int,Vector2D,Vector2D,double c){return variant==0?.25*c:variant==5?1e-100*c:0.;});
        const auto base=intervalFor(*mesh,old,false,rotation);auto phases=gas;if(variant==1)phases[1]=P::Liquid;
        const double begin=variant==3?oldTime+duration*.5:oldTime;
        const Interval extent(*mesh,old,base.velocitySource(),base.characteristicSource(),begin,endTime);
        const History history(*mesh,extent,phases);
        const Interval bound(*mesh,old,base.velocitySource(),base.characteristicSource(),begin,endTime,history);
        const HermiteGalsPhaseRegionSource regions(*mesh,old);const HermiteGalsInventoryRegionSource inventory(*mesh,regions,{.3});
        if(variant==4){const auto other=stateFor(*mesh,pure(1),false);const HermiteGalsPhaseRegionSource foreignRegions(*mesh,other);
            const HermiteGalsInventoryRegionSource foreign(*mesh,foreignRegions,{.3});
            rejects([&]{integrateHermiteGalsPolynomialFace(*mesh,bound,face.row,foreign,{1,1,1});},"foreign inventory");}
        else if(variant==2)rejects([&]{integrateHermiteGalsPolynomialFace(*mesh,bound,face.row,{1,1,1});},"pending=");
        else rejects([&]{integrateHermiteGalsPolynomialFace(*mesh,bound,face.row,inventory,{1,1,1});},"pending=");
    }
    // Shifted walls, extreme distances/speeds, and zero coordinates. These
    // compare the support against actual successful native events/interiors.
    for(double wall:{0.,.3})for(double speed:{.125,1e-300,std::numeric_limits<double>::denorm_min()}){
        const Slab field(oldTime,endTime,translation(true,speed));
        const LevelSetBoundaryTrace::Domain domain{{wall,wall},{wall+1,wall+1},{.25,.25},false};
        const auto reach=field.encloseRk2Reachability(domain,endTime);
        for(double y:{wall,std::nextafter(wall,wall+1),wall+1e-100,wall+.001,wall+.5}){
            const Vector2D point{wall+.5,y};const auto support=reach.cover(point,point);
            const auto trace=LevelSetBoundaryTrace::traceWithJacobian(domain,point,duration,
                [&](Vector2D p){return field.evaluate(p,endTime);},[&](Vector2D p){return field.evaluate(p,oldTime);});
            const auto p=trace.trajectory.point;check(p.x==point.x && support.lower.x==point.x && support.upper.x==point.x
                && p.y>=support.lower.y && p.y<=support.upper.y,"extreme source escaped support");++sourceBoxes;
        }
    }
    const Slab slab(oldTime,endTime,translation(false,.125));
    const LevelSetBoundaryTrace::Domain domain{{0,0},{1,1},{.25,.25},false};
    auto opaque=domain;opaque.activeSegment=[](Vector2D,Vector2D){return true;};
    rejects([&]{slab.encloseRk2Reachability(opaque,endTime);},"opaque");
    auto radial=domain;radial.radialAxisSymmetry=true;rejects([&]{slab.encloseRk2Reachability(radial,endTime);},"radial");
    rejects([&]{Slab(oldTime,endTime,translation(false,2)).encloseRk2Reachability(domain,endTime);},"half-cell");
    rejects([&]{Slab(.1,.375,translation(false,.125)).encloseRk2Reachability(domain,.375);},"history time");
    rejects([&]{slab.encloseRk2Reachability(domain,oldTime);},"positive duration");
    rejects([&]{slab.encloseRk2Reachability(domain,endTime).cover({-.1,0},{0,0});},"outside closed");
    const int rounding=std::fegetround();std::fesetround(FE_DOWNWARD);
    try{rejects([&]{slab.encloseRk2Reachability(domain,endTime);},"round-to-nearest");}
    catch(...){std::fesetround(rounding);throw;}std::fesetround(rounding);
    check(cases==24 && hits>0 && interior>0 && mapRefusals>0 && sourceBoxes==30,"uniform source matrix incomplete");
    std::cout<<"I12a16r5 uniform reachable PASS cases="<<cases<<" actual="<<actual<<" hits="<<hits<<" interior="<<interior
        <<" map_refusals="<<mapRefusals<<" initial_multi="<<initialMulti<<" bracket_multi="<<bracketMulti
        <<" tangent="<<tangent<<" extreme_boxes="<<sourceBoxes<<" guards="<<guards-oldGuards<<'\n';
}
void certifiedCirculationContract()
{
    using Quantum::Discrete::BoundedCirculationArc;using Quantum::Discrete::solveBoundedCirculation;
    using Arc=BoundedCirculationArc;int networks=0,batches=0;const int oldGuards=guards;
    const double eta=std::numeric_limits<double>::denorm_min();
    const std::vector<std::vector<Arc>> examples{
        {{0,1,.5,.5,.5},{1,0,0,1,0}},
        {{0,1,.5,.5,.5},{1,2,0,1,0},{2,3,0,1,0},{3,0,0,1,0}},
        {{0,1,.5,.5,.5},{1,2,0,.25,0},{2,0,0,.25,0},{1,3,0,.25,0},{3,0,0,.25,0}},
        {{0,1,.125,.75,.5},{1,2,.125,.75,.25},{2,0,.125,.75,.125}},
        {{0,1,eta,eta,eta},{1,0,0,eta,0}},
        {{0,1,1,1,1},{0,1,eta,eta,eta},{1,0,1,std::nextafter(1.,2.),1}},
        {{0,1,.5,.5,.5},{1,0,.5,.5,.5}},{}
    };
    for(const auto& arcs:examples){std::size_t nodes=2;for(const auto a:arcs)nodes=std::max(nodes,1+std::max(a.from,a.to));
        const auto result=solveBoundedCirculation(nodes,arcs),again=solveBoundedCirculation(nodes,arcs);
        check(result.flows==again.flows && result.augmentations==again.augmentations
            && result.roundedBalance==again.roundedBalance,"circulation is not deterministic");
        for(std::size_t k=0;k<arcs.size();++k)check(result.flows[k]>=arcs[k].lower && result.flows[k]<=arcs[k].upper,"selected arc escaped bounds");
        if(networks==5){std::cout<<"EXACT_CIRCULATION_ROUNDING,"<<result.flows[2]<<','<<result.roundedBalance[0][0]<<','<<result.roundedBalance[0][1]
                <<','<<result.roundedBalance[1][0]<<','<<result.roundedBalance[1][1]<<','<<eta<<'\n';
            check(result.flows[2]==1 && result.roundedBalance[0][0]==-eta && result.roundedBalance[0][1]==-eta
            && result.roundedBalance[1][0]==eta && result.roundedBalance[1][1]==eta,"rounded exact residual was erased");}
        else {std::vector<long double> balance(nodes,0);for(std::size_t k=0;k<arcs.size();++k){balance[arcs[k].from]-=result.flows[k];balance[arcs[k].to]+=result.flows[k];}
            for(std::size_t k=0;k<nodes;++k)check(balance[k]==0 && result.roundedBalance[k][0]==0 && result.roundedBalance[k][1]==0,"dyadic circulation ledger mismatch");}
        if(networks==1)check(result.flows==std::vector<double>({.5,.5,.5,.5}),"multi-edge transfer reference mismatch");
        if(networks==2)check(result.flows==std::vector<double>({.5,.25,.25,.25,.25}),"split capacity reference mismatch");
        if(networks==4)check(result.flows[1]==eta,"subnormal residual disappeared");
        if(networks>=6)check(result.augmentations==0,"already feasible network changed");++networks;
    }
    rejects([&]{solveBoundedCirculation(2,{{0,1,.5,.5,.5},{1,0,0,.25,0}});},"infeasible");
    // Every individual node has a possible balancing interval, yet the cut
    // {0,1}|{2,3} has no return capacity for the fixed cross-cut flux.
    rejects([&]{solveBoundedCirculation(4,{{2,0,1,1,1},{0,1,0,1,0},{1,0,0,1,0},{2,3,0,1,0},{3,2,0,1,0}});},"infeasible");
    rejects([&]{solveBoundedCirculation(0,{});},"node count");
    rejects([&]{solveBoundedCirculation(2,{{0,0,0,1,0}});},"endpoints");
    rejects([&]{solveBoundedCirculation(2,{{0,2,0,1,0}});},"endpoints");
    rejects([&]{solveBoundedCirculation(2,{{0,1,-eta,1,0}});},"bounded arc");
    rejects([&]{solveBoundedCirculation(2,{{0,1,0,1,2}});},"bounded arc");
    rejects([&]{solveBoundedCirculation(2,{{0,1,0,std::numeric_limits<double>::infinity(),0}});},"bounded arc");
    const int mode=std::fegetround();std::fesetround(FE_UPWARD);
    try{rejects([&]{solveBoundedCirculation(2,{});},"round-to-nearest");}catch(...){std::fesetround(mode);throw;}std::fesetround(mode);
    static_assert(!std::is_default_constructible_v<HermiteGalsCertifiedInventoryResult>);
    static_assert(std::is_same_v<decltype(std::declval<HermiteGalsCertifiedInventoryResult&>().inventory()),const HermiteGalsIntegratedInventoryResult&>);
    using History=HermiteGalsBoundaryPhaseSource;using P=HermiteGalsInflowPhase;
    for(bool radial:{false,true})for(bool amr:{false,true})for(double fraction:{0.,.25,1.}){
        auto mesh=makeMesh(amr);const auto old=stateFor(*mesh,pure(fraction==1?-1:1),radial,
            [=](int,Vector2D,Vector2D,double c){return fraction*c;});
        const auto base=intervalFor(*mesh,old,amr,translation(radial,.125));const auto phase=fraction==1?P::Liquid:P::Gas;
        const History history(*mesh,base,{phase,phase,phase,phase});
        const Interval interval(*mesh,old,base.velocitySource(),base.characteristicSource(),oldTime,endTime,history);
        const HermiteGalsPhaseRegionSource regions(*mesh,old);const HermiteGalsInventoryRegionSource inventory(*mesh,regions,{.3});
        FaceSegmentGraph graph;graph.build(*mesh);const auto domain=NativeLeafDomain::build(*mesh,graph,radial);
        std::vector<Result> certificates;std::vector<HermiteGalsIntegratedFaceTransfer> transfers;
        for(const auto f:domain.faces(*mesh)){const auto result=integrateHermiteGalsPolynomialFace(*mesh,interval,f.nativeRow,inventory);
            certificates.push_back(result);transfers.push_back(result.transfer());}
        const auto original=assembleHermiteGalsIntegratedInventory(*mesh,interval,transfers);
        const auto result=assembleHermiteGalsCertifiedInventory(*mesh,interval,certificates);
        check(!result.report().selected && result.report().changedFaces==0 && result.report().augmentations==0
            && result.inventory().interval().sameInterval(interval) && result.certificates().size()==certificates.size(),"already legal certified batch changed");
        for(std::size_t k=0;k<original.phaseVolumes().size();++k)check(original.phaseVolumes()[k].liquid==result.inventory().phaseVolumes()[k].liquid,"legal inventory changed");
        for(std::size_t k=0;k<transfers.size();++k)check(original.faces()[k].signedLiquidVolume()==result.inventory().faces()[k].signedLiquidVolume(),"legal shared flux changed");
        if(!radial && !amr && fraction==0){
            auto bad=certificates;bad.pop_back();rejects([&]{assembleHermiteGalsCertifiedInventory(*mesh,interval,bad);},"missing/duplicate/unordered");
            bad=certificates;std::swap(bad[0],bad[1]);rejects([&]{assembleHermiteGalsCertifiedInventory(*mesh,interval,bad);},"missing/duplicate/unordered");
            bad=certificates;bad.push_back(certificates[0]);rejects([&]{assembleHermiteGalsCertifiedInventory(*mesh,interval,bad);},"foreign/inactive");
            const History otherHistory(*mesh,base,{phase,phase,phase,phase});
            const Interval other(*mesh,old,base.velocitySource(),base.characteristicSource(),oldTime,endTime,otherHistory);
            rejects([&]{assembleHermiteGalsCertifiedInventory(*mesh,other,certificates);},"different interval source");
            bad=certificates;bad[0]=integrateHermiteGalsPolynomialFace(*mesh,interval,transfers[0].face().nativeRow);
            rejects([&]{assembleHermiteGalsCertifiedInventory(*mesh,interval,bad);},"same old inventory source");
        }
        ++batches;
    }
    check(networks==8 && batches==12,"incomplete certified selection matrix");
    std::cout<<"I12a16r6 certified circulation PASS networks="<<networks<<" legal_batches="<<batches<<" guards="<<guards-oldGuards<<'\n';
}
}
int main() try
{
    std::cout<<std::scientific<<std::setprecision(16);coreMoments();tightTimeTraces();
    for(bool radial:{false,true})for(bool amr:{false,true})translations(radial,amr);
    curvedVelocity(false);curvedVelocity(true);guardsAndLifetime();
    check(integrals+matrixFailures==57 && thinEvents==4 && weightedCases==8 && guards>=24 && tightTraces==864,"incomplete integral matrix");
    std::cout<<"I12a15 bounded spacetime "<<(matrixFailures?"FAIL":"PASS")<<" integrals="<<integrals<<" core="<<coreCases<<" panels="<<totalPanels<<" evaluated="<<totalEvaluated<<" max_panels="<<maxPanels<<" hidden="<<thinEvents<<" weighted="<<weightedCases<<" guards="<<guards<<" matrix_failures="<<matrixFailures<<'\n';
    if(!matrixFailures)boundaryHistoryContract();
    if(!matrixFailures)monotoneWallContract();
    if(!matrixFailures)tangentZeroContract();
    if(!matrixFailures)uniformReachableContract();
    if(!matrixFailures)certifiedCirculationContract();
    return matrixFailures?1:0;
}catch(const std::exception& e){std::cerr<<"I12a15 bounded spacetime FAIL: "<<e.what()<<'\n';return 1;}
