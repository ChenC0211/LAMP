#include "Core/Geometry/PolynomialVelocitySlab.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include <array>
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#if defined(_M_X64) || (defined(__x86_64__) && defined(__SSE2__))
#include <xmmintrin.h>
#endif
namespace
{
using namespace Quantum::Geometry;
using Data=PolynomialVelocityData;
using Slab=PolynomialVelocitySlab;
using Proof=PolynomialAffineRk2FaceBounds;
constexpr double oldTime=.25,endTime=.375;
int cases=0,traces=0,identity=0,guards=0;
double maximumRemainder=0,minimumShearRatio=1e100;
void check(bool ok,const std::string& why){if(!ok)throw std::runtime_error(why);}
template<class F> void rejects(F action,const char* message)
{try{action();}catch(const std::exception& e){check(std::string(e.what()).find(message)!=std::string::npos,e.what());++guards;return;}throw std::runtime_error(std::string("missing rejection: ")+message);}
double radius(const PolynomialAffineCoordinate& c){return std::max(std::abs(c.remainder.lower),std::abs(c.remainder.upper));}
void contains(const std::array<PolynomialVelocityRange,2>& box,Vector2D p)
{check(box[0].lower<=p.x && p.x<=box[0].upper && box[1].lower<=p.y && p.y<=box[1].upper,"departure escaped affine point range");}
Data flow(int kind)
{
    Data d;
    if(kind==1)d.previous[0]={.125,.0625};
    if(kind==2 || kind==3){d.origin={.5,.5};d.previous[1].x=-.5;d.previous[4].y=.5;}
    if(kind==4)d.previous[1].x=.5;
    if(kind==5 || kind==6){d.origin={.5,.5};d.previous[9].x=-.5;d.previous[3].x=-.5;d.previous[12].y=.5;d.previous[6].y=.5;}
    if(kind==7)d.previous[0].y=.125;
    if(kind==8 || kind==9){d.previous[4].x=.125;d.previous[0].y=.25;d.previous[1].y=-.25;}
    if(kind==10)d.previous[0]={1e-310,-2e-310};
    d.arrival=d.previous;if(kind==3 || kind==6 || kind==9)for(auto& c:d.arrival)c=c*(-.5);return d;
}
std::unique_ptr<QuadtreeMesh> meshFor(bool amr)
{DictionaryMesh d;d.Nx=d.Ny=4;d.maxCellNum=1024;d.maxInitialIterations=0;auto mesh=std::make_unique<QuadtreeMesh>(d);
    if(amr){DictionaryMesh::RefineBox b;b.level=1;b.min={.25,.25};b.max={.75,.75};mesh->applyRefineBoxes({b});}return mesh;}
std::array<Vector2D,2> faceFor(const QuadtreeMesh& mesh,bool radial,bool horizontal,bool axis)
{
    FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,radial);
    for(const auto f:domain.faces(mesh)){if(f.horizontal!=horizontal)continue;
        const auto& edge=horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
        auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
        if(horizontal?a.x>b.x:a.y>b.y)std::swap(a,b);
        if((horizontal?a.y:a.x)==(axis?0:.5) && (horizontal?a.x:a.y)==.25)return {a,b};}
    throw std::runtime_error("actual native face missing");
}
void samples(const Slab& slab,const LevelSetBoundaryTrace::Domain& domain,const Proof& proof,int kind,bool axis)
{
    ++cases;for(const auto& c:proof.departure)maximumRemainder=std::max(maximumRemainder,radius(c));
    for(int ti=0;ti<=16;++ti)for(int si=0;si<=8;++si){
        const double t=ti==16?proof.timeEnd:proof.timeBegin+(proof.timeEnd-proof.timeBegin)*ti/16;
        const auto p=si==8?proof.faceEnd:proof.faceBegin+(proof.faceEnd-proof.faceBegin)*(si/8.);
        const auto box=proof.encloseDeparture(p,t);const double dt=t-slab.timeBegin();
        if(dt==0){contains(box,p);++identity;continue;}
        const auto trace=LevelSetBoundaryTrace::traceWithJacobian(domain,p,dt,
            [&](Vector2D q){return slab.evaluate(q,t);},[&](Vector2D q){return slab.evaluate(q,slab.timeBegin());});
        check(!trace.trajectory.boundaryHit && trace.acceptedSubsteps==1 && bool(trace.derivative),"affine proof changed original branch");
        contains(box,trace.trajectory.point);++traces;
        // Independent original discrete algebra for constant/shear/rotation.
        if(kind==0)contains(box,p);
        if(kind==1)contains(box,{p.x-.125*dt,p.y-.0625*dt});
        if(kind==4)contains(box,{p.x-.5*p.y*dt,p.y});
        if(kind==7)contains(box,{p.x,p.y-.125*dt});
        if(kind==2){const double x=p.x-.5,y=p.y-.5,a=1-.125*dt*dt,s=.5*dt;contains(box,{.5+a*x+s*y,.5-s*x+a*y});}
        if(axis)check(trace.trajectory.point.x==0 && box[0].lower==0 && box[0].upper==0,"axis affine query displaced");
        if(kind==1)check(box[0].upper-box[0].lower<1e-10,"affine result disguised a full departure box");
    }
    if(axis){const auto& x=proof.departure[0];check(x.center==0 && x.tangent==0 && x.time==0 && x.remainder.lower==0 && x.remainder.upper==0,"raw affine axis lost exact zero");}
}
void matrix(bool radial,bool amr)
{
    auto mesh=meshFor(amr);const double h=amr?.125:.25;const LevelSetBoundaryTrace::Domain domain{mesh->minCoord,mesh->maxCoord,{h,h},radial};
    const auto kinds=radial?std::vector<int>{0,7,8,9}:std::vector<int>{0,1,2,3,4,5,6,10};
    for(int kind:kinds)for(bool horizontal:{false,true}){const auto face=faceFor(*mesh,radial,horizontal,false);const Slab slab(oldTime,endTime,flow(kind));
        const auto proof=slab.encloseAffineFaceRk2(domain,face[0],face[1],oldTime,endTime);samples(slab,domain,proof,kind,false);}
    if(radial)for(int kind:{0,7,8,9}){const Slab slab(oldTime,endTime,flow(kind));const auto face=faceFor(*mesh,true,false,true);
        samples(slab,domain,slab.encloseAffineFaceRk2(domain,face[0],face[1],oldTime,endTime),kind,true);}
    std::cout<<"AFFINE_MATRIX,"<<radial<<','<<amr<<",PASS\n";
}
void refinementAndNonDyadic()
{
    const LevelSetBoundaryTrace::Domain domain{{0,0},{1,1},{.25,.25},false};const Slab shear(oldTime,endTime,flow(4));
    double previous=0;
    for(int level=0;level<5;++level){const double space=std::ldexp(.125,-level),time=std::ldexp(.0625,-level);
        const auto proof=shear.encloseAffineFaceRk2(domain,{.5,.375-space},{.5,.375+space},.3125-time,.3125+time);
        const double error=radius(proof.departure[0]);check(error>=.5*space*time,"shear omitted its bilinear remainder");
        if(level){const double ratio=previous/error;minimumShearRatio=std::min(minimumShearRatio,ratio);check(ratio>3.5,"affine shear remainder did not shrink quadratically");}
        previous=error;samples(shear,domain,proof,4,false);std::cout<<"AFFINE_SHEAR_REMAINDER,"<<level<<','<<error<<'\n';}
    Data data;data.origin={.137,.221};data.scale={.73,.39};data.previous[0]={.015,-.01};data.previous[9].x=.007;data.previous[6].y=-.005;
    data.arrival=data.previous;data.arrival[0].x=-.02;data.arrival[4].y=.009;const Slab nonDyadic(.2,.3,data);
    for(bool horizontal:{false,true}){const Vector2D a=horizontal?Vector2D{.33,.43}:Vector2D{.43,.33},b=horizontal?Vector2D{.47,.43}:Vector2D{.43,.47};
        samples(nonDyadic,domain,nonDyadic.encloseAffineFaceRk2(domain,a,b,.217,.273),-1,false);}
}
void subnormalAndGuards()
{
    const LevelSetBoundaryTrace::Domain domain{{0,0},{1,1},{.25,.25},false};Data reversing;reversing.previous[4].x=.1;reversing.arrival[4].x=-.1;
    const Slab slab(0,.125,reversing);const double tiny=std::numeric_limits<double>::denorm_min();
    const auto tinyProof=slab.encloseAffineFaceRk2(domain,{.5,.4},{.5,.6},0,3*tiny);
    for(double t:{0.,tiny,2*tiny,3*tiny}){const Vector2D p{.5,.5};const auto box=tinyProof.encloseDeparture(p,t);
        if(t==0){contains(box,p);++identity;continue;}
        const auto trace=LevelSetBoundaryTrace::traceWithJacobian(domain,p,t,[&](Vector2D q){return slab.evaluate(q,t);},[&](Vector2D q){return slab.evaluate(q,0);});
        check(!trace.trajectory.boundaryHit && trace.acceptedSubsteps==1,"tiny original branch changed");contains(box,trace.trajectory.point);++traces;}
    const auto proof=slab.encloseAffineFaceRk2(domain,{.5,.4},{.5,.6},0,.125);
    rejects([&]{slab.encloseAffineFaceRk2(domain,{.5,.4},{.6,.5},0,.125);},"ascending");
    rejects([&]{slab.encloseAffineFaceRk2(domain,{.5,.6},{.5,.4},0,.125);},"ascending");
    rejects([&]{slab.encloseAffineFaceRk2(domain,{.5,.4},{.5,.6},0,0);},"positive");
    rejects([&]{slab.encloseAffineFaceRk2(domain,{.5,.4},{.5,.6},0,.25);},"outside slab");
    rejects([&]{proof.encloseDeparture({.6,.5},.1);},"retained face");
    rejects([&]{proof.encloseDeparture({.5,.7},.1);},"parameter extent");
    rejects([&]{proof.encloseDeparture({.5,.5},.2);},"parameter extent");
    const double nan=std::numeric_limits<double>::quiet_NaN(),inf=std::numeric_limits<double>::infinity();
    rejects([&]{slab.encloseAffineFaceRk2(domain,{nan,.4},{.5,.6},0,.125);},"ascending");
    rejects([&]{proof.encloseDeparture({.5,.5},nan);},"parameter extent");
    rejects([&]{proof.encloseDeparture({.5,inf},.1);},"retained face");
    rejects([&]{Proof{}.encloseDeparture({0,0},0);},"retained face");
    auto invalidDomain=domain;invalidDomain.cellSize.x=0;
    rejects([&]{slab.encloseAffineFaceRk2(invalidDomain,{.5,.4},{.5,.6},0,.125);},"invalid finite");
    auto bad=domain;bad.activeSegment=[](Vector2D,Vector2D){return true;};rejects([&]{slab.encloseAffineFaceRk2(bad,{.5,.4},{.5,.6},0,.125);},"opaque");
    const Slab moving(0,.125,flow(1));rejects([&]{moving.encloseAffineFaceRk2(domain,{0,.4},{0,.6},0,.125);},"boundary hit");
    bad=domain;bad.radialAxisSymmetry=true;rejects([&]{moving.encloseAffineFaceRk2(bad,{0,.4},{0,.6},0,.125);},"radial axis");
    Data fast;fast.previous[0].x=2;fast.arrival=fast.previous;const Slab speed(0,.125,fast);
    rejects([&]{speed.encloseAffineFaceRk2(domain,{.5,.4},{.5,.6},0,.125);},"half-cell speed guard");
    {struct Restore{int old=std::fegetround();~Restore(){std::fesetround(old);}} restore;
        check(std::fesetround(FE_UPWARD)==0,"cannot set rounding");rejects([&]{proof.encloseDeparture({.5,.5},.1);},"round-to-nearest");}
#if defined(_M_X64) || (defined(__x86_64__) && defined(__SSE2__))
    {struct Restore{unsigned old=_mm_getcsr();~Restore(){_mm_setcsr(old);}} restore;
        _mm_setcsr(restore.old|(1u<<15));rejects([&]{slab.encloseAffineFaceRk2(domain,{.5,.4},{.5,.6},0,.125);},"gradual underflow");}
#endif
    std::cout<<"AFFINE_SUBNORMAL_AND_GUARDS,PASS guards="<<guards<<'\n';
}
}
int main() try
{
    std::cout<<std::scientific<<std::setprecision(16);
    for(bool radial:{false,true})for(bool amr:{false,true})matrix(radial,amr);
    refinementAndNonDyadic();subnormalAndGuards();
    check(cases==63 && traces>9000 && guards==18,"incomplete affine matrix");
    std::cout<<"I12a15r2 affine characteristic PASS cases="<<cases<<" traces="<<traces<<" identity="<<identity<<" guards="<<guards
        <<" max_remainder="<<maximumRemainder<<" min_shear_ratio="<<minimumShearRatio<<'\n';
}catch(const std::exception& e){std::cerr<<"I12a15r2 affine characteristic FAIL: "<<e.what()<<'\n';return 1;}
