#include "Modules/Interface/HermiteGalsTemporalPhaseFlux.h"
#include "Core/Geometry/ActiveQuarticHermitePatch.h"
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <array>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
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
using namespace LAMP::Interface;
using Data=LAMP::Geometry::PolynomialVelocityData;
using Slab=LAMP::Geometry::PolynomialVelocitySlab;
using Bounds=LAMP::Geometry::PolynomialVelocityBounds;
using Velocity=HermiteGalsTemporalVelocitySource;
using Position=HermiteGalsTimeLayer::Position;
using Jet=LevelSetBoundaryTrace::VelocityJet;
using Matrix=LevelSetBoundaryTrace::Jacobian;
static_assert(!std::is_default_constructible_v<Slab>);
static_assert(std::is_same_v<decltype(std::declval<Slab&>().data()),const Data&>);
static_assert(std::is_same_v<decltype(std::declval<Velocity&>().polynomialVelocity()),const Slab*>);
int rejections=0,samples=0,traces=0;
double maximumReferenceError=0,maximumTraceError=0,maximumJacobianDifference=0,maximumPointWidth=0;
void check(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
template<class F> void rejects(F action,const char* reason)
{
    try{action();}catch(const std::exception& e){check(std::string(e.what()).find(reason)!=std::string::npos,e.what());++rejections;return;}
    throw std::runtime_error(std::string("missing rejection: ")+reason);
}
std::unique_ptr<QuadtreeMesh> makeMesh()
{DictionaryMesh d;d.Nx=d.Ny=4;d.maxCellNum=64;d.maxInitialIterations=0;return std::make_unique<QuadtreeMesh>(d);}
// Structural source fixture only. No claim of advancing its interface or M.
HermiteGalsState stateFor(const QuadtreeMesh& mesh,bool radial,HermiteGalsTimeLayer time={})
{
    LAMP::Geometry::ActiveQuarticHermitePatch patch;HermiteGalsPrimaryState p;HermiteGalsClosure closure;
    p.measure=radial?HermiteGalsMeasure::AxisymmetricVolume:HermiteGalsMeasure::CartesianArea;
    patch.activeCells.assign(mesh.cellIndexNum(),1);patch.cellCenters.resize(mesh.cellIndexNum());
    for(int v=0;v<mesh.vertIndexNum();++v){const auto q=mesh.verts(mesh.vertIndexToLabel[v]).coord();
        p.vertices.push_back({v,q.y-.4,0,1});patch.nodes.push_back({q.y-.4,0,1,0});}
    std::set<std::array<int,2>> edges;
    for(int c=0;c<mesh.cellIndexNum();++c){const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);
        const auto lo=mesh.verts(cell.cellToVert[0]).coord(),hi=mesh.verts(cell.cellToVert[2]).coord();
        p.geometryCells.push_back(c);p.phaseVolumes.push_back({c,std::clamp(.4-lo.y,0.,hi.y-lo.y)
            *(radial?std::acos(-1.)*(hi.x*hi.x-lo.x*lo.x):hi.x-lo.x)});
        patch.cellCenters[c]=cell.centerCoord().y-.4;
        for(int k=0;k<4;++k){const auto& a=mesh.verts(cell.cellToVert[k]);const auto& b=mesh.verts(cell.cellToVert[(k+1)%4]);
            std::array<int,2> key{{a.index(),b.index()}};
            if(a.coord().y==b.coord().y?a.coord().x>b.coord().x:a.coord().y>b.coord().y)std::swap(key[0],key[1]);edges.insert(key);}}
    for(const auto key:edges){const auto a=mesh.verts(mesh.vertIndexToLabel[key[0]]).coord(),b=mesh.verts(mesh.vertIndexToLabel[key[1]]).coord();
        const double value=(a.y+b.y)*.5-.4,dn=a.y==b.y?1.:0.;p.edges.push_back({key,value,dn});patch.edgeMidpoints.push_back({key,value,dn});}
    closure.vertexMixedDerivatives.assign(p.vertices.size(),0);closure.cellBubbleCoefficients=patch.cellCenters;
    HermiteGalsCandidate candidate(mesh,time,p);candidate.bindGeometry(mesh,candidate.identity(),closure,
        LevelSetHermiteGeometry::Snapshot::fromNativeLeaves(mesh,patch,radial));return std::move(candidate).finish();
}
std::array<double,6> components(Jet j)
{return {j.value.x,j.value.y,j.derivative.xx,j.derivative.xy,j.derivative.yx,j.derivative.yy};}
double difference(Jet a,Jet b)
{const auto x=components(a),y=components(b);double error=0;for(int i=0;i<6;++i)error=std::max(error,std::abs(x[i]-y[i]));return error;}
std::array<LAMP::Geometry::PolynomialVelocityRange,6> components(const Bounds& b)
{return {b.value[0],b.value[1],b.derivative[0],b.derivative[1],b.derivative[2],b.derivative[3]};}
struct Flow
{
    int kind=0;double oldAmplitude=.8,newAmplitude=.8;bool radial=false;
    Data data() const
    {
        Data result;std::array<Vector2D,16> c{};
        if(kind==0 || kind==1){result.origin={.5,.5};
            if(kind==0){c[1].x=-1;c[4].y=1;}
            else{c[9].x=-1;c[3].x=-1;c[12].y=1;c[6].y=1;}}
        else if(kind==2){result.origin={.2,-.3};result.scale={2,.5};
            c[0]={.1,-.2};c[4].x=.2;c[1].x=-.3;c[6].x=.4;c[8].y=.15;c[3].y=.1;}
        else if(kind==3){c[4].x=1;c[0].y=1;c[1].y=-2;}
        else if(kind==4){c[0]={0,1};}
        for(int i=0;i<16;++i){result.previous[i]=c[i]*oldAmplitude;result.arrival[i]=c[i]*newAmplitude;}return result;
    }
    // Analytic expressions/differentiation independent of coefficient Horner.
    Jet reference(Vector2D p,double time,double begin=.2,double end=.3) const
    {
        const long double theta=(static_cast<long double>(time)-begin)/(static_cast<long double>(end)-begin),
            k=(1-theta)*oldAmplitude+theta*newAmplitude;
        long double u=0,v=0,ux=0,uy=0,vx=0,vy=0;
        if(kind==0 || kind==1){const long double x=static_cast<long double>(p.x)-.5L,y=static_cast<long double>(p.y)-.5L;
            if(kind==0){u=-y;v=x;uy=-1;vx=1;}
            else{u=-(x*x+y*y)*y;v=(x*x+y*y)*x;ux=-2*x*y;uy=-x*x-3*y*y;vx=3*x*x+y*y;vy=2*x*y;}}
        else if(kind==2){const long double x=(static_cast<long double>(p.x)-.2)/2,y=(static_cast<long double>(p.y)+.3)/.5;
            u=.1+.2*x-.3*y+.4*x*y*y;v=-.2+.15*x*x+.1*y*y*y;
            ux=(.2+.4*y*y)/2;uy=(-.3+.8*x*y)/.5;vx=.3*x/2;vy=.3*y*y/.5;}
        else if(kind==3){u=p.x;v=1-2*static_cast<long double>(p.y);ux=1;vy=-2;}
        else if(kind==4)v=1;
        return {{double(k*u),double(k*v)},{double(k*ux),double(k*uy),double(k*vx),double(k*vy)}};
    }
};
void contained(const Bounds& bounds,Jet value)
{const auto a=components(bounds);const auto b=components(value);for(int i=0;i<6;++i)check(b[i]>=a[i].lower && b[i]<=a[i].upper,"velocity/derivative escaped interval");}
void flowContract(const QuadtreeMesh& mesh,const Flow& f,int id)
{
    const auto old=stateFor(mesh,f.radial,{2,.1,.2,Position::End});const auto identity=old.identity();const auto inventory=old.primary().phaseVolumes;
    const auto* geometry=&old.geometry().snapshot();HermiteGalsCandidate candidate(mesh,{3,.2,.3,Position::End},old.primary());
    auto input=f.data();auto velocity=Velocity::polynomial(old,candidate.identity(),candidate.timeLayer(),input);
    check(velocity.polynomialVelocity() && !velocity.constantValue(),"polynomial capability missing/conflated with constant");
    const Slab& field=*velocity.polynomialVelocity();input.previous[0]={999,999};
    const auto bounds=field.enclose({.1,.1},{.9,.9},.2,.3);
    for(int it=0;it<=8;++it)for(int ix=0;ix<=8;++ix)for(int iy=0;iy<=8;++iy){
        const double time=it==0?.2:it==8?.3:.2+(.3-.2)*it/8;
        const Vector2D point{ix==8?.9:.1+.8*ix/8,iy==8?.9:.1+.8*iy/8};
        const auto actual=field.evaluate(point,time),exact=f.reference(point,time);++samples;
        const double error=difference(actual,exact);maximumReferenceError=std::max(maximumReferenceError,error);
        check(error<=2e-14,"analytic velocity/derivative error");contained(bounds,actual);contained(bounds,exact);
        const auto local=field.enclose(point,point,time,time);contained(local,actual);contained(local,exact);
        for(const auto r:components(local)){maximumPointWidth=std::max(maximumPointWidth,r.upper-r.lower);check(r.upper-r.lower<=2e-12,"point enclosure unnecessarily unresolved");}}
    const LevelSetBoundaryTrace::Domain domain{mesh.minCoord,mesh.maxCoord,{.25,.25},f.radial};
    const HermiteGalsTemporalCharacteristicSource characteristic(mesh,domain);
    const HermiteGalsPhaseFluxIntervalSource interval(mesh,old,velocity,characteristic,.2,.3);
    auto same=velocity;check(same.sameSource(velocity) && same.polynomialVelocity()==velocity.polynomialVelocity(),"copy lost polynomial owner");
    const auto other=Velocity::polynomial(old,candidate.identity(),candidate.timeLayer(),f.data());check(!other.sameSource(velocity),"independent polynomial sources conflated");
    const auto opaque=Velocity(old,candidate.identity(),candidate.timeLayer(),[f](Vector2D p){return f.reference(p,.2);},[f](Vector2D p){return f.reference(p,.3);});
    check(!opaque.polynomialVelocity(),"opaque samples acquired polynomial certificate");
    rejects([&]{opaque.polynomialArrivalView(candidate.identity(),candidate.timeLayer());},"owned polynomial");
    rejects([&]{integrateHermiteGalsPlanarFace(mesh,interval,0);},"owned constant");
    for(const HermiteGalsTimeLayer layer:std::array<HermiteGalsTimeLayer,4>{{{3,.2,.225,Position::End},
            {3,.2,.25,Position::End},{3,.2,.3,Position::End},{3,.2,.3,Position::Midpoint}}}){
        const double end=layer.time();HermiteGalsCandidate actualCandidate(mesh,layer,old.primary());
        const auto arrival=velocity.polynomialArrivalView(actualCandidate.identity(),actualCandidate.timeLayer());
        check(arrival.geometryIdentity()==actualCandidate.identity() && arrival.timeLayer().time()==end,"substage lost actual candidate/time");
        const HermiteGalsVelocityView oldReference(old.identity(),old.timeLayer(),[f](Vector2D p){return f.reference(p,.2);});
        const HermiteGalsVelocityView endReference(actualCandidate.identity(),actualCandidate.timeLayer(),[f,end](Vector2D p){return f.reference(p,end);});
        const Vector2D point{f.radial?0.:.62,.57};
        const auto a=traceHermiteGalsWithJacobian(domain,point,velocity.previous(),arrival),b=traceHermiteGalsWithJacobian(domain,point,oldReference,endReference);
        check(!a.boundaryHit() && !b.boundaryHit(),"unexpected interior trace boundary hit");++traces;
        const auto ja=a.interiorJacobian(),jb=b.interiorJacobian();
        const double error=difference({a.departurePoint(),ja},{b.departurePoint(),jb});maximumTraceError=std::max(maximumTraceError,error);
        check(error<=5e-14,"polynomial/callback discrete B/J mismatch");
        if(f.radial)check(a.departurePoint().x==0 && field.evaluate({0,.57},end).derivative.xy==0,"axis coherence lost");
        const Vector2D interior{.62,.57};const auto base=traceHermiteGalsWithJacobian(domain,interior,velocity.previous(),arrival).interiorJacobian();
        const double delta=1e-5;const std::array<double,4> expected{{base.xx,base.xy,base.yx,base.yy}};
        for(int axis=0;axis<2;++axis){auto p=interior,m=interior;if(axis==0){p.x+=delta;m.x-=delta;}else{p.y+=delta;m.y-=delta;}
            const auto plus=traceHermiteGalsWithJacobian(domain,p,velocity.previous(),arrival).departurePoint(),minus=traceHermiteGalsWithJacobian(domain,m,velocity.previous(),arrival).departurePoint();
            const double x=(plus.x-minus.x)/(2*delta),y=(plus.y-minus.y)/(2*delta);
            const double e=std::max(std::abs(x-expected[axis]),std::abs(y-expected[axis+2]));maximumJacobianDifference=std::max(maximumJacobianDifference,e);check(e<=2e-9,"finite-difference B/J mismatch");}
    }
    auto moved=std::move(velocity);rejects([&]{velocity.polynomialVelocity();},"moved-from velocity");
    check(moved.sameSource(same) && interval.velocitySource().sameSource(moved),"move/interval lost source owner");
    check(old.identity()==identity && &old.geometry().snapshot()==geometry,"source query changed geometry authority");
    for(std::size_t i=0;i<inventory.size();++i)check(old.primary().phaseVolumes[i].liquid==inventory[i].liquid,"source query changed M");
    std::cout<<"POLYNOMIAL_SOURCE,"<<id<<','<<f.kind<<','<<f.radial<<",PASS\n";
}
void analyticBounds()
{
    Data data;data.previous[4].x=1;data.previous[8].x=-1;data.previous[1].y=-1;data.previous[3].y=1;data.arrival=data.previous;
    const Slab field(0,1,data);const auto wide=field.enclose({0,-1},{1,1},0,1);
    check(wide.value[0].upper>=.25,"missed quadratic interior extremum");
    const double extremum=2/(3*std::sqrt(3.));check(wide.value[1].lower<=-extremum && wide.value[1].upper>=extremum,"missed cubic interior extrema");
    contained(wide,field.evaluate({.5,1/std::sqrt(3.)},.37));
    const auto small=field.enclose({.4,.4},{.6,.6},.4,.6);
    for(int i=0;i<2;++i)check(small.value[i].upper-small.value[i].lower<wide.value[i].upper-wide.value[i].lower,"box refinement did not narrow range");
    Data crossing;crossing.previous[0].x=-1;crossing.arrival[0].x=1;const Slab reverse(0,1,crossing);
    const auto bounds=reverse.enclose({0,0},{1,1},0,1);check(bounds.value[0].lower<=-1 && bounds.value[0].upper>=1,"time reversal excluded");
    check(reverse.evaluate({.4,.4},.5).value.x==0,"midpoint time reversal");
    const Slab zero(0,1,Data{});for(const auto r:components(zero.enclose({-1,-2},{3,4},0,1)))check(r.lower==0 && r.upper==0,"zero polynomial enclosure");
    std::cout<<"ANALYTIC_BOUNDS,quadratic_and_cubic_interior_extrema,time_reversal,box_refinement,zero,PASS\n";
}
void guardsAndLifetime(const QuadtreeMesh& mesh)
{
    const auto old=stateFor(mesh,false);const HermiteGalsTimeLayer time{1,0,.1,Position::End};HermiteGalsCandidate next(mesh,time,old.primary());
    Data data;data.previous[0]={.1,.2};data.arrival=data.previous;const Slab field(0,.1,data);
    for(int k=0;k<3;++k){auto bad=data;if(k==0)bad.scale.x=0;if(k==1)bad.origin.x=std::numeric_limits<double>::infinity();if(k==2)bad.scale.y=-1;
        rejects([&]{Slab invalid(0,.1,bad);},"coordinate basis");}
    auto bad=data;bad.arrival[15].x=std::numeric_limits<double>::quiet_NaN();rejects([&]{Slab invalid(0,.1,bad);},"coefficient");
    rejects([&]{Slab invalid(.1,.1,data);},"time slab");
    rejects([&]{field.evaluate({0,0},.2);},"outside slab");
    rejects([&]{field.evaluate({std::numeric_limits<double>::infinity(),0},.05);},"spatial query");
    rejects([&]{field.enclose({1,0},{0,1},0,.1);},"spatial box");
    rejects([&]{field.enclose({0,0},{1,1},.08,.03);},"outside slab");
    Data overflow;overflow.previous[12].x=std::numeric_limits<double>::max();overflow.arrival=overflow.previous;
    const Slab huge(0,1,overflow);rejects([&]{huge.evaluate({2,0},.5);},"nonfinite");rejects([&]{huge.enclose({1,0},{2,1},0,1);},"non-finite");
    {struct Restore{int old=std::fegetround();~Restore(){std::fesetround(old);}} restore;
        check(std::fesetround(FE_UPWARD)==0,"cannot set test rounding mode");rejects([&]{field.enclose({0,0},{1,1},0,.1);},"round-to-nearest");}
#if defined(_M_X64) || (defined(__x86_64__) && defined(__SSE2__))
    {struct Restore{unsigned old=_mm_getcsr();~Restore(){_mm_setcsr(old);}} restore;
        _mm_setcsr(restore.old|(1u<<15));rejects([&]{field.evaluate({.3,.4},.05);},"gradual underflow");}
#endif
    const auto source=Velocity::polynomial(old,next.identity(),time,data);
    rejects([&]{source.polynomialArrivalView({},time);},"unbound geometry");
    rejects([&]{source.polynomialArrivalView(next.identity(),{2,0,.05,Position::End});},"outside source");
    rejects([&]{source.polynomialArrivalView(next.identity(),{1,0,.2,Position::End});},"outside source");
    rejects([&]{source.polynomialArrivalView(next.identity(),{1,.01,.05,Position::End});},"outside source");
    rejects([&]{source.constantArrivalView(next.identity(),time);},"owned constant");
    const auto constant=Velocity::constant(old,next.identity(),time,{.1,.2});check(!constant.polynomialVelocity(),"constant source gained polynomial owner");
    rejects([&]{constant.polynomialArrivalView(next.identity(),time);},"owned polynomial");
    std::optional<HermiteGalsVelocityView> retainedOld,retainedEnd;
    {auto local=Velocity::polynomial(old,next.identity(),time,data);retainedOld=local.previous();retainedEnd=local.polynomialArrivalView(next.identity(),time);}
    const LevelSetBoundaryTrace::Domain domain{mesh.minCoord,mesh.maxCoord,{.25,.25},false};
    const auto trace=traceHermiteGalsWithJacobian(domain,{.5,.5},*retainedOld,*retainedEnd);
    check(std::abs(trace.departurePoint().x-.49)<=1e-15 && std::abs(trace.departurePoint().y-.48)<=1e-15,"retained substage lost lifetime");
    Data boundary;boundary.previous[0].x=1;boundary.arrival=boundary.previous;const auto inflow=Velocity::polynomial(old,next.identity(),time,boundary);
    const auto hit=traceHermiteGalsWithJacobian(domain,{0,.5},inflow.previous(),inflow.arrival());
    check(hit.boundaryHit() && !hit.jacobian(),"polynomial provider bypassed boundary event");
    boundary.previous[0].x=2;boundary.arrival=boundary.previous;const auto fast=Velocity::polynomial(old,next.identity(),time,boundary);
    rejects([&]{traceHermiteGalsWithJacobian(domain,{.5,.5},fast.previous(),fast.arrival());},"half a cell");
    const auto radialOld=stateFor(mesh,true);HermiteGalsCandidate radialNext(mesh,time,radialOld.primary());
    boundary.previous[0].x=.1;boundary.arrival=boundary.previous;const auto badAxis=Velocity::polynomial(radialOld,radialNext.identity(),time,boundary);
    auto radialDomain=domain;radialDomain.radialAxisSymmetry=true;
    rejects([&]{traceHermiteGalsWithJacobian(radialDomain,{0,.5},badAxis.previous(),badAxis.arrival());},"symmetry axis");
    std::cout<<"SOURCE_GUARDS_AND_LIFETIME,PASS rejections="<<rejections<<'\n';
}
void convergence(const QuadtreeMesh& mesh)
{
    double last=0;const Vector2D start{.7,.55};const double omega=.8,total=.2;
    const LevelSetBoundaryTrace::Domain domain{mesh.minCoord,mesh.maxCoord,{.25,.25},false};
    for(int steps:{4,8,16}){const double dt=total/steps;Vector2D point=start;
        for(int k=0;k<steps;++k){const double begin=k*dt,end=(k+1)*dt;
            const auto old=stateFor(mesh,false,k==0?HermiteGalsTimeLayer{}:HermiteGalsTimeLayer{std::int64_t(k),(k-1)*dt,begin,Position::End});
            HermiteGalsCandidate next(mesh,{std::int64_t(k+1),begin,end,Position::End},old.primary());
            Flow flow;flow.oldAmplitude=flow.newAmplitude=omega;
            const auto source=Velocity::polynomial(old,next.identity(),next.timeLayer(),flow.data());
            const auto trace=traceHermiteGalsWithJacobian(domain,point,source.previous(),source.arrival());
            check(!trace.boundaryHit() && trace.acceptedSubsteps()==1,"rotation convergence changed RK2 branch");point=trace.departurePoint();}
        const double c=std::cos(omega*total),s=std::sin(omega*total);
        const Vector2D exact{.5+c*(start.x-.5)+s*(start.y-.5),.5-s*(start.x-.5)+c*(start.y-.5)};
        const double error=std::hypot(point.x-exact.x,point.y-exact.y),order=last?std::log2(last/error):0;
        if(last)check(error<last && order>=1.95 && order<=2.1,"original RK2 rotation time order");last=error;
        std::cout<<"POLYNOMIAL_RK2_ROTATION,"<<steps<<','<<dt<<','<<total<<','<<error<<','<<order<<'\n';}
}
void signedBoundsContract()
{
    int boxes=0,values=0,zeros=0,signs=0;const int oldRejections=rejections;
    const std::array<std::array<Vector2D,2>,7> regions{{
        {{{.1,.1},{.9,.9}}},{{{.5,0},{.5,.125}}},{{{0,.5},{.125,.5}}},{{{.5,.5},{.5,.5}}},
        {{{0,0},{.5,.5}}},{{{.5,.5},{1,1}}},{{{.25,.5},{.5,.75}}}}};
    const std::array<Flow,8> flows{{{0,.8,.8,false},{0,.8,-.6,false},{1,.8,.8,false},{1,.8,-.6,false},
        {2,1.1,-.7,false},{3,.4,-.2,true},{4,.3,-.2,true},{5,0,0,false}}};
    for(const auto f:flows){const Slab slab(.2,.3,f.data());
        for(std::size_t region=0;region<regions.size();++region){const auto lo=regions[region][0],hi=regions[region][1];
            const auto bounds=slab.encloseSignPreserving(lo,hi,.2,.3),outer=slab.enclose(lo,hi,.2,.3);
            const auto a=components(bounds),b=components(outer);
            for(int k=0;k<6;++k)check(a[k].lower>=b[k].lower && a[k].upper<=b[k].upper,"signed bound not subset of old certificate");
            if(f.kind<=1 && (region==1 || region==2 || region==3)){
                const int k=region==2?0:1;check(bounds.value[k].lower==0 && bounds.value[k].upper==0,"center-coordinate exact zero lost");++zeros;
                check(bounds.derivative[0].lower==0 && bounds.derivative[0].upper==0
                    && bounds.derivative[3].lower==0 && bounds.derivative[3].upper==0,"center-coordinate derivative zero lost");}
            for(int ti=0;ti<=8;++ti)for(int xi=0;xi<=4;++xi)for(int yi=0;yi<=4;++yi){
                const double t=ti==8?.3:.2+(.3-.2)*ti/8;
                const Vector2D p{xi==4?hi.x:lo.x+(hi.x-lo.x)*xi/4,yi==4?hi.y:lo.y+(hi.y-lo.y)*yi/4};
                contained(bounds,slab.evaluate(p,t));contained(bounds,f.reference(p,t));++values;}
            ++boxes;
        }
    }
    // Full tensor powers with non-dyadic origin/scale on either side of each
    // chart axis. Signs follow the actual multiplication/Horner expression.
    for(int i=0;i<4;++i)for(int j=0;j<4;++j)for(int sx:{-1,1})for(int sy:{-1,1}){
        Data d;d.origin={.3,-.7};d.scale={.2,1.7};d.previous[4*i+j]={.25,-.125};d.arrival[4*i+j]={.5,-.25};
        const Slab slab(.25,.375,d);Vector2D lo{d.origin.x+(sx<0?-.1:0),d.origin.y+(sy<0?-.2:0)},
            hi{d.origin.x+(sx>0?.1:0),d.origin.y+(sy>0?.2:0)};
        const auto b=slab.encloseSignPreserving(lo,hi,.25,.375);const int sign=(i%2?sx:1)*(j%2?sy:1);
        check(sign>0?b.value[0].lower>=0 && b.value[1].upper<=0:b.value[0].upper<=0 && b.value[1].lower>=0,"Horner sign certificate lost");
        for(int ti=0;ti<=4;++ti)for(int xi=0;xi<=4;++xi)for(int yi=0;yi<=4;++yi){
            const Vector2D p{xi==4?hi.x:lo.x+(hi.x-lo.x)*xi/4,yi==4?hi.y:lo.y+(hi.y-lo.y)*yi/4};
            contained(b,slab.evaluate(p,.25+.125*ti/4));++values;}
        ++signs;++boxes;
    }
    for(bool horizontal:{false,true}){
        Data d;d.origin={.3,-.7};d.scale={.2,1.7};d.previous[horizontal?1:4]={1,1};d.arrival=d.previous;
        const Slab slab(.25,.375,d);const double origin=horizontal?d.origin.y:d.origin.x;
        for(int sign:{-1,1}){auto p=d.origin;
            (horizontal?p.y:p.x)=std::nextafter(origin,sign>0?std::numeric_limits<double>::infinity():-std::numeric_limits<double>::infinity());
            const auto b=slab.encloseSignPreserving(p,p,.25,.375);contained(b,slab.evaluate(p,.3125));
            check(sign>0?b.value[0].upper>0 && b.value[0].lower>=0:b.value[0].lower<0 && b.value[0].upper<=0,"adjacent representable nonzero collapsed");++signs;
        }
    }
    for(double tiny:{std::numeric_limits<double>::denorm_min(),1e-310}){
        Data d;d.previous[0]={tiny,-tiny};d.arrival=d.previous;const Slab slab(0,1,d);
        const auto b=slab.encloseSignPreserving({0,0},{1,1},0,1);
        check(b.value[0].upper>0 && b.value[0].lower>=0 && b.value[1].lower<0 && b.value[1].upper<=0,"subnormal sign/zero false certificate");
        for(int ti=0;ti<=16;++ti)contained(b,slab.evaluate({.5,.5},ti/16.));++signs;
    }
    Data d;d.previous[0]={1,-1};d.arrival[0]={-1,1};const Slab reversal(0,1,d);
    const auto reversed=reversal.encloseSignPreserving({0,0},{1,1},0,1);
    check(reversed.value[0].lower<=-1 && reversed.value[0].upper>=1 && reversed.value[1].lower<=-1 && reversed.value[1].upper>=1,"time reversal lost a sign");
    rejects([&]{reversal.encloseSignPreserving({1,0},{0,1},0,1);},"spatial box");
    rejects([&]{reversal.encloseSignPreserving({0,0},{1,1},-.1,1);},"outside slab");
    {struct Restore{int old=std::fegetround();~Restore(){std::fesetround(old);}} restore;
        check(std::fesetround(FE_UPWARD)==0,"cannot set rounding");rejects([&]{reversal.encloseSignPreserving({0,0},{1,1},0,1);},"round-to-nearest");}
    check(boxes==120 && values==20600 && zeros==12 && signs==70 && rejections-oldRejections==3,"signed velocity matrix incomplete");
    std::cout<<"I12a16r3 signed velocity PASS boxes="<<boxes<<" samples="<<values<<" exact_zeros="<<zeros<<" signs="<<signs<<" guards="<<rejections-oldRejections<<'\n';
}
}
int main() try
{
    std::cout<<std::scientific<<std::setprecision(16);auto mesh=makeMesh();
    const std::array<Flow,8> flows{{{0,.8,.8,false},{0,.8,-.6,false},{1,.8,.8,false},{1,.8,-.6,false},
        {2,1.1,-.7,false},{3,.4,-.2,true},{4,.3,-.2,true},{5,0,0,false}}};
    for(int i=0;i<int(flows.size());++i)flowContract(*mesh,flows[i],i);
    analyticBounds();guardsAndLifetime(*mesh);convergence(*mesh);
    check(samples==5832 && traces==32 && rejections>=44,"incomplete polynomial source matrix");
    std::cout<<"I12a12 polynomial velocity PASS flows="<<flows.size()<<" samples="<<samples<<" traces="<<traces<<" rejections="<<rejections
        <<" max_reference="<<maximumReferenceError<<" max_trace="<<maximumTraceError<<" max_fd_J="<<maximumJacobianDifference<<" max_point_width="<<maximumPointWidth<<'\n';
    signedBoundsContract();
}catch(const std::exception& e){std::cerr<<"I12a12 polynomial velocity FAIL: "<<e.what()<<'\n';return 1;}
