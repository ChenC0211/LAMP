#include "Core/Geometry/PolynomialVelocitySlab.h"
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Geometry/ActiveQuarticHermitePatch.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include <algorithm>
#include <array>
#include <cfenv>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#if defined(_M_X64) || (defined(__x86_64__) && defined(__SSE2__))
#include <xmmintrin.h>
#endif
namespace
{
using namespace LAMP::Geometry;
using Snapshot=LevelSetHermiteGeometry::Snapshot;
using Patch=ActiveQuarticHermitePatch;
using Jet=Patch::NodalDof;
using Profile=std::function<Jet(Vector2D)>;
using Form=PolynomialAffineCoordinate;
using Forms=std::array<Form,2>;
using Proof=PolynomialAffineRk2FaceBounds;
using Slab=PolynomialVelocitySlab;
using Data=PolynomialVelocityData;
using Key=std::array<int,2>;
constexpr double oldTime=.25,endTime=.375;
int cases=0,traces=0,identities=0,leafChecks=0,guards=0,coreChecks=0,hidden=0;
double minimumRatio=1e100;
const char* stage="setup";
void check(bool ok,const std::string& why){if(!ok)throw std::runtime_error(why);}
template<class F> void rejects(F action,const char* message)
{try{action();}catch(const std::exception& e){check(std::string(e.what()).find(message)!=std::string::npos,e.what());++guards;return;}
    throw std::runtime_error(std::string("missing rejection: ")+message);}
void contains(PolynomialVelocityRange r,double v)
{check(r.lower<=v && v<=r.upper,"original scalar escaped affine certificate");}
double radius(const Form& f){return std::max(std::abs(f.remainder.lower),std::abs(f.remainder.upper));}
std::unique_ptr<QuadtreeMesh> meshFor(bool amr,int n=4,Vector2D lo={0,0},Vector2D hi={1,1})
{DictionaryMesh d;d.Nx=d.Ny=n;d.minCoord=lo;d.maxCoord=hi;d.maxCellNum=4096;d.maxInitialIterations=0;
    auto mesh=std::make_unique<QuadtreeMesh>(d);
    if(amr){DictionaryMesh::RefineBox b;b.level=1;b.min=lo+(hi-lo)*.25;b.max=lo+(hi-lo)*.75;mesh->applyRefineBoxes({b});}return mesh;}
Patch patchFor(const QuadtreeMesh& mesh,const Profile& jet)
{
    Patch patch;patch.activeCells.assign(mesh.cellIndexNum(),1);patch.nodes.resize(mesh.vertIndexNum());patch.cellCenters.resize(mesh.cellIndexNum());
    for(int v=0;v<mesh.vertIndexNum();++v)patch.nodes[v]=jet(mesh.verts(mesh.vertIndexToLabel[v]).coord());
    std::map<Key,bool> sides;
    for(int c=0;c<mesh.cellIndexNum();++c){const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);patch.cellCenters[c]=jet(cell.centerCoord()).value;
        for(int k=0;k<4;++k){Key key{{mesh.verts(cell.cellToVert[k]).index(),mesh.verts(cell.cellToVert[(k+1)%4]).index()}};
            const auto a=mesh.verts(mesh.vertIndexToLabel[key[0]]).coord(),b=mesh.verts(mesh.vertIndexToLabel[key[1]]).coord();
            if(a.y==b.y?a.x>b.x:a.y>b.y)std::swap(key[0],key[1]);sides.emplace(key,true);}}
    for(const auto& side:sides){const auto a=mesh.verts(mesh.vertIndexToLabel[side.first[0]]).coord(),b=mesh.verts(mesh.vertIndexToLabel[side.first[1]]).coord();
        const auto j=jet(a+(b-a)*.5);patch.edgeMidpoints.push_back({side.first,j.value,a.y==b.y?j.dy:j.dx});}return patch;
}
Profile profile(int kind,bool radial=false)
{
    return [=](Vector2D p){const double x=radial?p.x:p.x-.5,y=p.y-.53125;
        if(kind==0)return radial?Jet{p.y-.484375,0,1,0}:Jet{p.x-.484375,1,0,0};
        if(kind==1)return Jet{x*x+y*y-.0625,2*x,2*y,0};
        if(kind==2)return Jet{x*x*x*x+y*y*y*y+.25*x*x*y*y+.125*x*x*x*x*y*y*y*y-.00390625,
            4*x*x*x+.5*x*y*y+.5*x*x*x*y*y*y*y,4*y*y*y+.5*x*x*y+.5*x*x*x*x*y*y*y,x*y+2*x*x*x*y*y*y};
        if(kind==3)return Jet{y*y-0x1p-20,0,2*y,0};
        if(kind==4){const double z=radial?p.x*p.x-.53125*.53125:p.x-.53125;
            return Jet{z*z+y*y-0x1p-14,radial?4*p.x*z:2*z,2*y,0};}
        return Jet{y*y,0,2*y,0};};
}
Data flow(int kind)
{
    Data d;
    if(kind==1)d.previous[0].y=.25;
    if(kind==2)d.previous[1].x=.5;
    if(kind==3 || kind==4){d.origin={.5,.5};d.previous[1].x=-.5;d.previous[4].y=.5;}
    if(kind==5){d.origin={.5,.5};d.previous[9].x=-.5;d.previous[3].x=-.5;d.previous[12].y=.5;d.previous[6].y=.5;}
    if(kind==6){d.previous[4].x=.125;d.previous[0].y=.25;d.previous[1].y=-.25;}
    d.arrival=d.previous;if(kind==4)for(auto& v:d.arrival)v=v*(-.5);return d;
}
std::array<Vector2D,2> faceFor(const QuadtreeMesh& mesh,bool radial,bool horizontal)
{
    FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,radial);
    for(const auto f:domain.faces(mesh)){if(f.horizontal!=horizontal)continue;
        const auto& edge=horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
        auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
        if(horizontal?a.x>b.x:a.y>b.y)std::swap(a,b);
        if((horizontal?a.y:a.x)==.5 && (horizontal?a.x:a.y)==.25)return {a,b};}
    throw std::runtime_error("actual native face missing");
}
bool inBox(Vector2D p,const std::array<Vector2D,2>& b)
{return b[0].x<=p.x && p.x<=b[1].x && b[0].y<=p.y && p.y<=b[1].y;}
void verify(const Snapshot& source,const Slab& slab,const LevelSetBoundaryTrace::Domain& domain,const Proof& proof,const Profile& reference={})
{
    struct Leaf {int index;std::array<Vector2D,2> box;Form scalar;std::array<double,25> before;};std::vector<Leaf> leaves;
    // Exhaustive closed intersection, independent of the candidate index.
    for(int c=0;c<static_cast<int>(source.nativeCellCount());++c){const auto b=source.nativeCellBounds(c);
        const auto lo=proof.branch.departureLower,hi=proof.branch.departureUpper;
        if(b[0].x>hi.x || b[1].x<lo.x || b[0].y>hi.y || b[1].y<lo.y)continue;
        leaves.push_back({c,b,source.nativeQuarticAffineEnclosure(c,proof.departure),source.nativeQuarticLocalDofs(c)});}
    check(!leaves.empty(),"empty closed leaf coverage");
    for(int j=0;j<=16;++j)for(int i=0;i<=8;++i){const double qs=-1+i/4.,qt=-1+j/8.;
        const auto p=proof.faceBegin+(proof.faceEnd-proof.faceBegin)*(i/8.);
        const double t=proof.timeBegin+(proof.timeEnd-proof.timeBegin)*(j/16.),dt=t-slab.timeBegin();Vector2D b=p;
        if(dt>0){const auto actual=LevelSetBoundaryTrace::traceWithJacobian(domain,p,dt,
                [&](Vector2D q){return slab.evaluate(q,t);},[&](Vector2D q){return slab.evaluate(q,slab.timeBegin());});
            check(!actual.trajectory.boundaryHit && actual.acceptedSubsteps==1 && bool(actual.derivative),"original RK2 branch changed");b=actual.trajectory.point;++traces;
        }else ++identities;
        int donors=0;
        for(const auto& leaf:leaves)if(inBox(b,leaf.box)){const auto range=leaf.scalar.enclose({qs,qs},{qt,qt});
            const double value=source.scalarJetInNativeCell(leaf.index,b).value;contains(range,value);
            if(reference)contains(range,reference(b).value);
            if(range.upper<0)check(source.negative(b),"negative certificate changed authoritative phase");
            if(range.lower>0)check(!source.negative(b),"positive certificate changed authoritative phase");++donors;++leafChecks;}
        check(donors>0,"actual departure missing from closed coverage");
    }
    for(const auto& leaf:leaves)check(source.nativeQuarticLocalDofs(leaf.index)==leaf.before,"query changed native source");++cases;
}
void matrix(bool radial,bool amr)
{
    auto mesh=meshFor(amr);const double h=amr?.125:.25;const LevelSetBoundaryTrace::Domain domain{mesh->minCoord,mesh->maxCoord,{h,h},radial};
    for(int shape=0;shape<6;++shape){const auto f=profile(shape,radial);const auto source=Snapshot::fromNativeLeaves(*mesh,patchFor(*mesh,f),radial);
        for(int kind:radial?std::vector<int>{0,1,6}:std::vector<int>{0,1,2,3,4,5})for(bool horizontal:{false,true}){
            const auto face=faceFor(*mesh,radial,horizontal);const Slab slab(oldTime,endTime,flow(kind));
            verify(source,slab,domain,slab.encloseAffineFaceRk2(domain,face[0],face[1],oldTime,endTime),f);}}
    std::cout<<"AFFINE_PHASE_MATRIX,"<<radial<<','<<amr<<",PASS\n";
}
void corePolynomial()
{
    std::array<std::array<double,2>,25> powers{};powers[0]={-.25,-.25};powers[4]={.5,.5};powers[20]={-.125,-.125};powers[12]={.25,.25};powers[24]={.125,.125};
    const Forms coordinates{{Form{.5,.0625,.03125,{-.0009765625,.0009765625}},Form{.375,-.03125,.0625,{0,.001953125}}}};
    const auto form=encloseTensorQuarticAffine(powers,{0,0},{1,1},coordinates);
    for(int j=0;j<=16;++j)for(int i=0;i<=16;++i)for(bool first:{false,true})for(bool second:{false,true}){
        const double s=-1+i/8.,t=-1+j/8.;const auto& a=coordinates[0];const auto& b=coordinates[1];
        const double x=a.center+a.tangent*s+a.time*t+(first?a.remainder.upper:a.remainder.lower),
            y=b.center+b.tangent*s+b.time*t+(second?b.remainder.upper:b.remainder.lower);
        const double exact=-.25+.5*x*x*x*x-.125*y*y*y*y+.25*x*x*y*y+.125*x*x*x*x*y*y*y*y;
        contains(form.enclose({s,s},{t,t}),exact);++coreChecks;}
    // No small coefficient deletion, including a subnormal scalar.
    powers={};const double tiny=std::numeric_limits<double>::denorm_min();powers[0]={tiny,tiny};
    contains(encloseTensorQuarticAffine(powers,{0,0},{1,1},coordinates).enclose({-1,1},{-1,1}),tiny);++coreChecks;
    // Plane cancellation: exact Q4 conversion must not invent curvature.
    auto mesh=meshFor(false);const auto source=Snapshot::fromNativeLeaves(*mesh,patchFor(*mesh,profile(0)));
    const auto plane=source.nativeQuarticAffineEnclosure(5,{{Form{.375,.0625,.03125,{0,0}},Form{.375,.03125,-.0625,{0,0}}}});
    check(radius(plane)<1e-12,"linear Q4 acquired nonlinear remainder");
    std::cout<<"AFFINE_PHASE_CORE,PASS samples="<<coreChecks<<" plane_R="<<radius(plane)<<'\n';
}
void refinement()
{
    auto mesh=meshFor(false);const LevelSetBoundaryTrace::Domain domain{{0,0},{1,1},{.25,.25},false};const Slab slab(oldTime,endTime,flow(2));
    for(int shape:{0,2}){const auto source=Snapshot::fromNativeLeaves(*mesh,patchFor(*mesh,profile(shape)));double previous=0;
        for(int level=0;level<5;++level){const double space=std::ldexp(.03125,-level),time=std::ldexp(.015625,-level);
            const auto proof=slab.encloseAffineFaceRk2(domain,{.5,.375-space},{.5,.375+space},.3125-time,.3125+time);
            const auto scalar=source.nativeQuarticAffineEnclosure(5,proof.departure);const double error=radius(scalar);
            if(level){const double ratio=previous/error;minimumRatio=std::min(minimumRatio,ratio);check(ratio>3.5,"Q4 remainder did not shrink quadratically");}
            previous=error;verify(source,slab,domain,proof,profile(shape));std::cout<<"AFFINE_PHASE_REMAINDER,"<<shape<<','<<level<<','<<error<<'\n';}}
    for(bool amr:{false,true}){auto originalMesh=meshFor(false,4,{.1,.2},{1.1,1.2}),targetMesh=meshFor(amr,amr?4:8,{.1,.2},{1.1,1.2});
        const auto original=Snapshot::fromNativeLeaves(*originalMesh,patchFor(*originalMesh,profile(2)));
        auto exact=Snapshot::refineNativeQuartic(*targetMesh,original);
        const auto again=Snapshot::refineNativeQuartic(*targetMesh,exact);
        // Same non-dyadic authoritative source, not an analytic reinitialization.
        const LevelSetBoundaryTrace::Domain shifted{{.1,.2},{1.1,1.2},{.125,.125},false};const Slab zero(oldTime,endTime,flow(0));
        const auto proof=zero.encloseAffineFaceRk2(shifted,{.6,.43},{.6,.57},oldTime,endTime);
        verify(original,zero,shifted,proof);verify(exact,zero,shifted,proof);verify(again,zero,shifted,proof);
        for(int c=0;c<static_cast<int>(again.nativeCellCount());++c){const auto box=again.nativeCellBounds(c);const auto p=box[0]+(box[1]-box[0])*.375;
            const auto v=again.nativeQuarticAffineEnclosure(c,{{Form{p.x,0,0,{0,0}},Form{p.y,0,0,{0,0}}}}).enclose({0,0},{0,0});
            contains(v,original.value(p));contains(v,again.value(p));++leafChecks;}}
    std::cout<<"AFFINE_PHASE_REFINEMENT,PASS min_ratio="<<minimumRatio<<'\n';
}
void hiddenEvents()
{
    auto mesh=meshFor(false);const LevelSetBoundaryTrace::Domain domain{{0,0},{1,1},{.25,.25},false};const Slab slab(oldTime,endTime,flow(1));
    for(int shape:{3,5}){const auto source=Snapshot::fromNativeLeaves(*mesh,patchFor(*mesh,profile(shape)));
        const auto proof=slab.encloseAffineFaceRk2(domain,{.5,.549},{.5,.551},oldTime,endTime);
        const auto full=source.nativeQuarticAffineEnclosure(9,proof.departure).enclose({-1,1},{-1,1});
        check(full.lower<=0 && full.upper>=0,"hidden band/tangency lost from full range");
        for(double t:{oldTime,.3125,endTime})check(!source.negative({.5,.55-.25*(t-oldTime)}),"hidden fixture has a sampled liquid node");
        if(shape==3){check(source.negative({.5,.53125}),"hidden band liquid missing");
            const auto inner=slab.encloseAffineFaceRk2(domain,{.5,.54999},{.5,.55001},.32499,.32501);
            check(source.nativeQuarticAffineEnclosure(9,inner.departure).enclose({-1,1},{-1,1}).upper<0,"resolved hidden band not certified");}
        ++hidden;}
    const auto bubble=Snapshot::fromNativeLeaves(*mesh,patchFor(*mesh,profile(4)));
    for(auto p:std::array<Vector2D,5>{{{.5,.5},{.625,.5},{.5,.625},{.625,.625},{.5625,.5625}}})check(!bubble.negative(p),"hidden bubble fixture node liquid");
    const Forms box{{Form{.5625,.0625,0,{0,0}},Form{.5625,0,.0625,{0,0}}}};
    const auto value=bubble.nativeQuarticAffineEnclosure(10,box).enclose({-1,1},{-1,1});
    check(value.lower<0 && value.upper>0 && bubble.negative({.53125,.53125}),"internal bubble missing");++hidden;
    std::cout<<"AFFINE_PHASE_HIDDEN,PASS cases="<<hidden<<'\n';
}
void rejectionCases()
{
    auto mesh=meshFor(false);const auto patch=patchFor(*mesh,profile(2));auto source=Snapshot::fromNativeLeaves(*mesh,patch);
    const Forms coordinates{{Form{.375,.03125,0,{0,0}},Form{.375,0,.03125,{0,0}}}};
    rejects([&]{source.nativeQuarticAffineEnclosure(-1,coordinates);},"index");
    rejects([&]{source.nativeQuarticAffineEnclosure(16,coordinates);},"index");
    const auto partial=Snapshot::fromNativeQuarticQueryCells(*mesh,patch,{5});
    rejects([&]{partial.nativeQuarticAffineEnclosure(0,coordinates);},"inactive");
    ActiveHermitePatch cubic;cubic.activeCells=patch.activeCells;cubic.nodes.assign(patch.nodes.size(),Jet{1,0,0,0});
    const auto cubicNative=Snapshot::fromNativeLeaves(*mesh,cubic),cubicGrid=Snapshot(*mesh,cubic);
    rejects([&]{cubicNative.nativeQuarticAffineEnclosure(5,coordinates);},"native quartic");
    rejects([&]{cubicGrid.nativeQuarticAffineEnclosure(5,coordinates);},"native quartic");
    auto moved=std::move(source);rejects([&]{source.nativeQuarticAffineEnclosure(5,coordinates);},"moved-from");
    const double inf=std::numeric_limits<double>::infinity(),nan=std::numeric_limits<double>::quiet_NaN();
    for(int k=0;k<4;++k){auto bad=coordinates;if(k==0)bad[0].center=nan;if(k==1)bad[1].tangent=inf;
        if(k==2)bad[0].remainder={1,-1};if(k==3)bad[1].time=nan;
        rejects([&]{moved.nativeQuarticAffineEnclosure(5,bad);},k==2?"interval":"coefficient");}
    std::array<std::array<double,2>,25> powers{};
    rejects([&]{encloseTensorQuarticAffine(powers,{0,0},{0,1},coordinates);},"bounds");
    rejects([&]{encloseTensorQuarticAffine(powers,{0,nan},{1,1},coordinates);},"bounds");
    powers[0]={2,1};rejects([&]{encloseTensorQuarticAffine(powers,{0,0},{1,1},coordinates);},"interval");
    powers[0]={0,inf};rejects([&]{encloseTensorQuarticAffine(powers,{0,0},{1,1},coordinates);},"interval");
    const auto value=moved.nativeQuarticAffineEnclosure(5,coordinates);
    rejects([&]{value.enclose({-2,1},{0,0});},"unit square");
    rejects([&]{value.enclose({0,0},{1,-1});},"interval");
    rejects([&]{value.enclose({nan,0},{0,0});},"interval");
    {struct Restore{int old=std::fegetround();~Restore(){std::fesetround(old);}} restore;
        check(std::fesetround(FE_UPWARD)==0,"cannot set rounding");
        rejects([&]{moved.nativeQuarticAffineEnclosure(5,coordinates);},"round-to-nearest");
        rejects([&]{value.enclose({0,0},{0,0});},"round-to-nearest");}
#if defined(_M_X64) || (defined(__x86_64__) && defined(__SSE2__))
    {struct Restore{unsigned old=_mm_getcsr();~Restore(){_mm_setcsr(old);}} restore;
        _mm_setcsr(restore.old|(1u<<15));rejects([&]{moved.nativeQuarticAffineEnclosure(5,coordinates);},"gradual underflow");}
#endif
    std::cout<<"AFFINE_PHASE_GUARDS,PASS guards="<<guards<<'\n';
}
void subsetIntegrals()
{
    int checks=0;const int oldGuards=guards;const double pi=std::acos(-1.);
    const auto verify=[](const PolynomialFaceIntegralBounds& result,double total,double liquid,int sign){
        for(int k=0;k<4;++k){const double expected=k==(sign>0?0:1)?total:k==(sign>0?2:3)?liquid:0;
            contains(result.volumes[k],expected);}
        const auto representative=result.representative();
        check(representative[2]<=representative[0] && representative[3]<=representative[1],"subset representative ordering");};
    for(bool radial:{false,true})for(bool horizontal:{false,true})for(int sign:{-1,1})for(double epsilon:{0.,.125})
        for(bool weighted:{false,true})for(bool hull:{false,true}){
        Data data;data.previous[0]=horizontal?Vector2D{0,sign*.25}:Vector2D{sign*.25,0};data.arrival=data.previous;
        const Slab slab(oldTime,endTime,data);const Vector2D a=horizontal?Vector2D{.25,.5}:Vector2D{.5,.25},b{.5,.5};
        const double wl=weighted?.75:1,wg=weighted?.25:0;
        const std::vector<Form> forms=hull?std::vector<Form>{{epsilon,1,1,{0,0}},{-epsilon,1,1,{0,0}}}:
            std::vector<Form>{{0,1,1,{-epsilon,epsilon}}};
        const auto result=slab.encloseNormalFaceIntegralWithAffinePhase(a,b,oldTime,endTime,radial,forms,wl,wl,wg,wg);
        const double scale=.25*.125*.25,fullWeight=radial?(horizontal?.375:.5):1,total=scale*(radial?2*pi:1)*fullWeight;
        for(double offset:{-epsilon,0.,epsilon}){
            const double k=1-offset/2,side=2-k;
            const double area=k<=1?k*k/2:1-side*side/2,moment=k<=1?k*k*k/6:.5-side*side/2+side*side*side/6;
            const double weight=radial?(horizontal?.25*area+.25*moment:.5*area):area;
            const double liquid=wg*total+(wl-wg)*scale*(radial?2*pi:1)*weight;
            verify(result,total,liquid,sign);
        }
        if(epsilon==0)check(result.errorBound()<1e-12*total,"exact affine half-plane remained unresolved");
        else check(result.errorBound()<.3*total,"certified subsets did not reduce gross uncertainty");++checks;
    }
    Data data;data.previous[0].x=.25;data.arrival=data.previous;const Slab slab(oldTime,endTime,data);
    for(int phase:{-1,0,1}){const auto result=slab.encloseNormalFaceIntegralWithAffinePhase({.5,.25},{.5,.5},oldTime,endTime,false,{{double(phase),0,0,{0,0}}});
        if(phase<0)verify(result,.0078125,.0078125,1);
        if(phase>0)verify(result,.0078125,0,1);
        if(phase==0)check(result.volumes[2].lower==0 && result.volumes[2].upper>=.0078125,"identically zero phase classified as pure");++checks;}
    const auto varyingWeights=slab.encloseNormalFaceIntegralWithAffinePhase({.5,.25},{.5,.5},oldTime,endTime,false,{{0,1,1,{0,0}}},.5,.75,0,.25);
    verify(varyingWeights,.0078125,.0078125*.25,1);verify(varyingWeights,.0078125,.0078125*.5,1);++checks;
    const auto differentSlopes=slab.encloseNormalFaceIntegralWithAffinePhase({.5,.25},{.5,.5},oldTime,endTime,false,{{0,1,1,{0,0}},{0,1.125,.875,{0,0}}});
    verify(differentSlopes,.0078125,.00390625,1);check(differentSlopes.errorBound()<.4*.0078125,"slope hull lost useful subset bounds");++checks;
    Data reverse=data;reverse.arrival[0].x=-.25;const Slab reversing(oldTime,endTime,reverse);
    const auto counter=reversing.encloseNormalFaceIntegralWithAffinePhase({.5,.25},{.5,.5},oldTime,endTime,false,{{0,1,1,{0,0}}});
    check(counter.volumes[0].upper>0 && counter.volumes[1].upper>0,"subset erased reversing gross flow");++checks;
    rejects([&]{slab.encloseNormalFaceIntegralWithAffinePhase({.5,.25},{.5,.5},oldTime,endTime,false,{});},"empty affine phase");
    rejects([&]{slab.encloseNormalFaceIntegralWithAffinePhase({.5,.25},{.5,.5},oldTime,endTime,false,{{0,1,1,{0,0}}},.25,.5,.5,.75);},"weights");
    rejects([&]{slab.encloseNormalFaceIntegralWithAffinePhase({.5,.25},{.5,.5},oldTime,endTime,false,{{0,1,1,{0,0}}},1.1,1.1);},"content");
    rejects([&]{slab.encloseNormalFaceIntegralWithAffinePhase({.5,.25},{.5,.5},oldTime,endTime,false,{{0,1,1,{1,-1}}});},"interval");
    rejects([&]{slab.encloseNormalFaceIntegralWithAffinePhase({.5,.25},{.5,.5},oldTime,endTime,false,{{0,std::numeric_limits<double>::infinity(),1,{0,0}}});},"coefficient");
    {struct Restore{int old=std::fegetround();~Restore(){std::fesetround(old);}} restore;check(std::fesetround(FE_UPWARD)==0,"cannot set rounding");
        rejects([&]{slab.encloseNormalFaceIntegralWithAffinePhase({.5,.25},{.5,.5},oldTime,endTime,false,{{0,1,1,{0,0}}});},"round-to-nearest");}
    check(checks==70 && guards-oldGuards==6,"incomplete subset matrix");
    std::cout<<"I12a15r4 exact phase subset PASS cases="<<checks<<" guards="<<guards-oldGuards<<'\n';
}
}
int main() try
{
    std::cout<<std::scientific<<std::setprecision(16);
    stage="core";corePolynomial();stage="matrix";
    for(bool radial:{false,true})for(bool amr:{false,true})matrix(radial,amr);
    stage="refinement";refinement();stage="hidden";hiddenEvents();stage="guards";rejectionCases();
    check(cases==232 && traces>30000 && leafChecks>traces && guards>=20 && hidden==3,"incomplete Q4 affine coverage");
    std::cout<<"I12a15r3 affine phase PASS cases="<<cases<<" traces="<<traces<<" identities="<<identities<<" leaf_checks="<<leafChecks
        <<" core_samples="<<coreChecks<<" guards="<<guards<<" hidden="<<hidden<<" min_ratio="<<minimumRatio<<'\n';
    stage="subset integrals";subsetIntegrals();
}catch(const std::exception& e){std::cerr<<"I12a15r3 affine phase FAIL stage="<<stage<<": "<<e.what()<<'\n';return 1;}
