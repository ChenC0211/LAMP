#include "Modules/Interface/HermiteGalsTemporalPhaseFlux.h"
#include "Modules/Interface/HermiteGalsAmrClosure.h"
#include "Core/Geometry/ActiveQuarticHermitePatch.h"
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <array>
#include <cfenv>
#include <cmath>
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
using namespace Quantum::Interface;
using Data=Quantum::Geometry::PolynomialVelocityData;
using Slab=Quantum::Geometry::PolynomialVelocitySlab;
using Bounds=Quantum::Geometry::PolynomialRk2Bounds;
using Range=Quantum::Geometry::PolynomialVelocityRange;
using Velocity=HermiteGalsTemporalVelocitySource;
using Characteristic=HermiteGalsTemporalCharacteristicSource;
using Interval=HermiteGalsPhaseFluxIntervalSource;
using Enclosure=HermiteGalsCharacteristicEnclosure;
using Position=HermiteGalsTimeLayer::Position;
using Key=std::array<int,2>;
static_assert(!std::is_default_constructible_v<Enclosure>);
static_assert(std::is_same_v<decltype(std::declval<Enclosure&>().bounds()),const Bounds&>);
int boxes=0,traces=0,rejections=0,identityChecks=0;
double maximumAnalyticError=0,maximumTravel=0;
void check(bool ok,const std::string& why){if(!ok)throw std::runtime_error(why);}
template<class F> void rejects(F action,const char* reason)
{
    try{action();}catch(const std::exception& e){check(std::string(e.what()).find(reason)!=std::string::npos,e.what());++rejections;return;}
    throw std::runtime_error(std::string("missing rejection: ")+reason);
}
std::unique_ptr<QuadtreeMesh> makeMesh(bool amr=false)
{
    DictionaryMesh d;d.Nx=d.Ny=4;d.maxCellNum=1024;d.maxInitialIterations=0;
    auto mesh=std::make_unique<QuadtreeMesh>(d);
    if(amr){DictionaryMesh::RefineBox b;b.level=1;b.min={.25,.25};b.max={.75,.75};mesh->applyRefineBoxes({b});}return mesh;
}
Vector2D point(const QuadtreeMesh& m,int v){return m.verts(m.vertIndexToLabel[v]).coord();}
bool fineSide(const QuadtreeMesh& m,Key key)
{
    const auto a=point(m,key[0]),b=point(m,key[1]);const bool horizontal=a.y==b.y;
    for(int v:key){const auto owner=m.verts(m.vertIndexToLabel[v]).supCell;if(owner<0)continue;
        const auto& c=m.cells(owner);const auto lo=m.verts(c.cellToVert[0]).coord(),hi=m.verts(c.cellToVert[2]).coord();
        if(horizontal && (a.y==lo.y || a.y==hi.y) && a.x>=lo.x && b.x<=hi.x && b.x-a.x<hi.x-lo.x)return true;
        if(!horizontal && (a.x==lo.x || a.x==hi.x) && a.y>=lo.y && b.y<=hi.y && b.y-a.y<hi.y-lo.y)return true;}return false;
}
// Complete native plane support is a structural source fixture, not an
// interface evolution or a reset of any transported inventory.
HermiteGalsState stateFor(const QuadtreeMesh& m,bool radial)
{
    FaceSegmentGraph graph;graph.build(m);const auto domain=NativeLeafDomain::build(m,graph,radial);
    Quantum::Geometry::ActiveQuarticHermitePatch patch;HermiteGalsPrimaryState p;
    p.measure=radial?HermiteGalsMeasure::AxisymmetricVolume:HermiteGalsMeasure::CartesianArea;
    patch.activeCells.assign(m.cellIndexNum(),1);patch.nodes.resize(m.vertIndexNum());patch.cellCenters.resize(m.cellIndexNum());
    for(int v=0;v<m.vertIndexNum();++v)if(m.verts(m.vertIndexToLabel[v]).type()==TVertType::non){const auto q=point(m,v);
        p.vertices.push_back({v,q.y-.4,0,1});patch.nodes[v]={q.y-.4,0,1,0};}
    std::set<Key> keys;
    for(int c=0;c<m.cellIndexNum();++c){const auto& cell=m.cells(m.cellIndexToLabel[c]);
        const auto lo=m.verts(cell.cellToVert[0]).coord(),hi=m.verts(cell.cellToVert[2]).coord();
        p.geometryCells.push_back(c);p.phaseVolumes.push_back({c,domain.cellVolumes()[c]*std::clamp((.4-lo.y)/(hi.y-lo.y),0.,1.)});
        patch.cellCenters[c]=cell.centerCoord().y-.4;
        for(int k=0;k<4;++k){Key key{{m.verts(cell.cellToVert[k]).index(),m.verts(cell.cellToVert[(k+1)%4]).index()}};
            const auto a=point(m,key[0]),b=point(m,key[1]);if(a.y==b.y?a.x>b.x:a.y>b.y)std::swap(key[0],key[1]);keys.insert(key);}}
    for(const auto key:keys)if(!fineSide(m,key)){const auto a=point(m,key[0]),b=point(m,key[1]);
        const double value=(a.y+b.y)*.5-.4,dn=a.y==b.y?1.:0.;p.edges.push_back({key,value,dn});patch.edgeMidpoints.push_back({key,value,dn});}
    HermiteGalsCandidate candidate(m,{2,.1,.2,Position::End},p);const auto mixed=closeHermiteGalsMixedAmr(m,candidate);
    HermiteGalsClosure closure;closure.vertexMixedDerivatives=mixed.vertexMixedDerivatives();closure.cellBubbleCoefficients=patch.cellCenters;
    for(std::size_t k=0;k<p.vertices.size();++k)patch.nodes[p.vertices[k].vertex].dxy=closure.vertexMixedDerivatives[k];
    candidate.bindGeometry(m,candidate.identity(),closure,LevelSetHermiteGeometry::Snapshot::fromNativeLeaves(m,patch,radial));return std::move(candidate).finish();
}
struct Flow
{
    int kind;double previous,arrival;bool radial;
    Data data() const
    {
        Data d;std::array<Vector2D,16> c{};
        if(kind<2){d.origin={.5,.5};if(kind==0){c[1].x=-1;c[4].y=1;}
            else{c[9].x=-1;c[3].x=-1;c[12].y=1;c[6].y=1;}}
        else if(kind==2)c[1].x=1;
        else if(kind==3){c[4].x=1;c[0].y=1;c[1].y=-2;}
        else if(kind==4)c[0].y=1;
        for(int i=0;i<16;++i){d.previous[i]=c[i]*previous;d.arrival[i]=c[i]*arrival;}return d;
    }
};
void contained(double value,Range range,const char* why)
{check(std::isfinite(value) && range.lower<=value && value<=range.upper,why);}
void contained(Vector2D value,Vector2D lo,Vector2D hi,const char* why)
{contained(value.x,{lo.x,hi.x},why);contained(value.y,{lo.y,hi.y},why);}
void contained(const Bounds& bounds,Vector2D departure,LevelSetBoundaryTrace::Jacobian j)
{
    contained(departure,bounds.departureLower,bounds.departureUpper,"departure escaped enclosure");
    const std::array<double,4> values{{j.xx,j.xy,j.yx,j.yy}};
    for(int k=0;k<4;++k)contained(values[k],bounds.jacobian[k],"Jacobian escaped enclosure");
}
// Closed-form *discrete* RK2 references for constant rotation, shear, axial
// translation and zero. Continuous rotation is intentionally not substituted.
void analytic(const Flow& flow,Vector2D p,double d,Vector2D b,LevelSetBoundaryTrace::Jacobian j)
{
    if(flow.previous!=flow.arrival || flow.kind==1 || flow.kind==3)return;
    const double k=flow.previous;Vector2D expected=p;LevelSetBoundaryTrace::Jacobian matrix{1,0,0,1};
    if(flow.kind==0){const double a=1-.5*k*k*d*d,s=k*d,x=p.x-.5,y=p.y-.5;
        expected={.5+a*x+s*y,.5-s*x+a*y};matrix={a,s,-s,a};}
    else if(flow.kind==2){expected.x=p.x-k*p.y*d;matrix.xy=-k*d;}
    else if(flow.kind==4)expected.y=p.y-k*d;
    const double error=std::max({std::abs(b.x-expected.x),std::abs(b.y-expected.y),std::abs(j.xx-matrix.xx),
        std::abs(j.xy-matrix.xy),std::abs(j.yx-matrix.yx),std::abs(j.yy-matrix.yy)});
    maximumAnalyticError=std::max(maximumAnalyticError,error);check(error<=8*std::numeric_limits<double>::epsilon(),"affine discrete reference mismatch");
}
void runBox(const QuadtreeMesh& mesh,bool amr,const Flow& flow,int id,bool axis=false)
{
    const auto old=stateFor(mesh,flow.radial);const auto inventory=old.primary().phaseVolumes;
    const auto identity=old.identity();const auto* geometry=&old.geometry().snapshot();
    HermiteGalsCandidate next(mesh,{3,.2,.3,Position::End},old.primary());
    const auto source=Velocity::polynomial(old,next.identity(),next.timeLayer(),flow.data());
    const double h=amr?.125:.25;const LevelSetBoundaryTrace::Domain domain{mesh.minCoord,mesh.maxCoord,{h,h},flow.radial};
    const Characteristic map(mesh,domain);const Interval interval(mesh,old,source,map,.2,.3);
    const Vector2D lo{axis?0.:.60,.55},hi{axis?0.:.62,.57};
    const auto result=encloseHermiteGalsCharacteristic(mesh,interval,lo,hi);const auto& b=result.bounds();++boxes;
    check(result.interval().sameInterval(interval) && result.interval().velocitySource().polynomialVelocity()==source.polynomialVelocity(),"lost source provenance");
    check(&result.interval().previousAccepted().geometry().snapshot()==geometry,"lost single geometry owner");
    maximumTravel=std::max(maximumTravel,b.maximumL1Travel);check(b.maximumL1Travel<h/2,"CFL witness not strict");
    for(int ti=0;ti<=8;++ti){const double t=ti==8?.3:.2+(.3-.2)*ti/8;
        for(int xi=0;xi<5;++xi)for(int yi=0;yi<5;++yi){const Vector2D p{lo.x+(hi.x-lo.x)*xi/4,lo.y+(hi.y-lo.y)*yi/4};
            if(ti==0){contained(b,p,{1,0,0,1});++identityChecks;continue;}
            HermiteGalsCandidate actual(mesh,{3,.2,t,Position::End},old.primary());
            const auto arrival=source.polynomialArrivalView(actual.identity(),actual.timeLayer());
            const auto trace=traceHermiteGalsWithJacobian(domain,p,source.previous(),arrival);++traces;
            check(!trace.boundaryHit() && trace.acceptedSubsteps()==1,"certificate accepted another trace branch");
            contained(b,trace.departurePoint(),trace.interiorJacobian());
            const auto half=p-source.polynomialVelocity()->evaluate(p,t).value*.5*(t-.2);
            contained(half,b.midpointLower,b.midpointUpper,"midpoint escaped enclosure");
            analytic(flow,p,t-.2,trace.departurePoint(),trace.interiorJacobian());
            if(axis)check(trace.departurePoint().x==0,"axis trace displaced");}}
    HermiteGalsCandidate mid(mesh,{3,.2,.3,Position::Midpoint},old.primary());
    const auto trace=traceHermiteGalsWithJacobian(domain,lo,source.previous(),source.polynomialArrivalView(mid.identity(),mid.timeLayer()));
    check(!trace.boundaryHit() && trace.acceptedSubsteps()==1,"midpoint view changed branch");contained(b,trace.departurePoint(),trace.interiorJacobian());++traces;
    // Closed spatial/time subset and point/face degeneracies keep provenance.
    const Interval narrow(mesh,old,source,map,.25,.26);
    const auto refined=encloseHermiteGalsCharacteristic(mesh,narrow,lo,{lo.x,lo.y+.001});const auto& r=refined.bounds();
    if(flow.kind!=5)check(r.departureUpper.y-r.departureLower.y<b.departureUpper.y-b.departureLower.y,"refinement failed to narrow departure");
    const auto pointBound=source.polynomialVelocity()->encloseInteriorRk2(domain,lo,lo,.3,.3);
    const auto last=traceHermiteGalsWithJacobian(domain,lo,source.previous(),source.arrival());contained(pointBound,last.departurePoint(),last.interiorJacobian());++traces;
    check(old.identity()==identity && &old.geometry().snapshot()==geometry,"enclosure modified geometry");
    for(std::size_t i=0;i<inventory.size();++i)check(old.primary().phaseVolumes[i].liquid==inventory[i].liquid,"enclosure modified M");
    std::cout<<"CHARACTERISTIC_BOX,"<<id<<','<<amr<<','<<flow.radial<<','<<axis<<','<<mesh.cellIndexNum()<<','<<b.maximumL1Travel<<",PASS\n";
}
void coreBranches()
{
    const LevelSetBoundaryTrace::Domain domain{{0,0},{1,1},{.25,.25},false};
    Data data;data.previous[0].x=.1;data.arrival=data.previous;const Slab field(0,.1,data);
    const auto zeroTime=field.encloseInteriorRk2(domain,{0,0},{1,1},0,0);
    check(zeroTime.departureLower.x==0 && zeroTime.departureUpper.y==1 && zeroTime.maximumL1Travel==0,"zero-time identity failed");
    const Slab zero(0,.1,Data{});const auto identity=zero.encloseInteriorRk2(domain,{0,0},{1,1},0,.1);
    check(identity.departureLower.x==0 && identity.departureUpper.x==1 && identity.departureLower.y==0 && identity.departureUpper.y==1,"zero-velocity box changed");
    for(int k=0;k<4;++k)check(identity.jacobian[k].lower==(k%3==0?1.:0.) && identity.jacobian[k].lower==identity.jacobian[k].upper,"zero Jacobian not identity");
    const auto outflow=field.encloseInteriorRk2(domain,{1,.4},{1,.6},.075,.1);
    for(double t:{.075,.0875,.1}){const auto q=LevelSetBoundaryTrace::traceWithJacobian(domain,{1,.5},t,
            [&](Vector2D p){return field.evaluate(p,t);},[&](Vector2D p){return field.evaluate(p,0);});
        check(!q.trajectory.boundaryHit && q.acceptedSubsteps==1,"outflow branch changed");contained(outflow,q.trajectory.point,*q.derivative);++traces;}
    rejects([&]{field.encloseInteriorRk2(domain,{0,.4},{0,.6},.075,.1);},"initial boundary hit");
    const auto hit=LevelSetBoundaryTrace::traceWithJacobian(domain,{0,.5},.1,
        [&](Vector2D p){return field.evaluate(p,.1);},[&](Vector2D p){return field.evaluate(p,0);});
    check(hit.trajectory.boundaryHit && !hit.derivative,"inflow refusal lost actual boundary hit");
    rejects([&]{field.encloseInteriorRk2(domain,{.001,.4},{.002,.6},.1,.1);},"midpoint domain unresolved");
    rejects([&]{field.encloseInteriorRk2(domain,{.007,.4},{.008,.6},.1,.1);},"target domain unresolved");
    // Endpoint contact with outward uncertainty still fails certification;
    // the exact native trace below independently confirms the hit.
    Data exact;exact.previous[0].x=.125;exact.arrival=exact.previous;const Slab contact(0,.125,exact);
    rejects([&]{contact.encloseInteriorRk2(domain,{.015625,.5},{.015625,.5},.125,.125);},"domain unresolved");
    const auto wall=LevelSetBoundaryTrace::traceWithJacobian(domain,{.015625,.5},.125,
        [&](Vector2D p){return contact.evaluate(p,.125);},[&](Vector2D p){return contact.evaluate(p,0);});
    check(wall.trajectory.boundaryHit,"exact contact did not hit native wall");
    auto active=domain;active.activeSegment=[](Vector2D,Vector2D){return true;};
    rejects([&]{field.encloseInteriorRk2(active,{.5,.5},{.6,.6},0,.1);},"opaque active coverage");
    auto radial=domain;radial.radialAxisSymmetry=true;
    rejects([&]{field.encloseInteriorRk2(radial,{0,.5},{0,.6},0,.1);},"radial axis");
    Data fast;fast.previous[0].x=2;fast.arrival=fast.previous;const Slab speed(0,.1,fast);
    rejects([&]{speed.encloseInteriorRk2(domain,{.5,.5},{.6,.6},0,.1);},"half-cell speed guard");
    rejects([&]{field.encloseInteriorRk2(domain,{-.1,.5},{.6,.6},0,.1);},"outside physical domain");
    rejects([&]{field.encloseInteriorRk2(domain,{.6,.5},{.5,.6},0,.1);},"spatial box");
    rejects([&]{field.encloseInteriorRk2(domain,{std::numeric_limits<double>::quiet_NaN(),.5},{.6,.6},0,.1);},"spatial box");
    rejects([&]{field.encloseInteriorRk2(domain,{.5,.5},{.6,.6},.05,.2);},"outside slab");
    rejects([&]{field.encloseInteriorRk2(domain,{.5,.5},{.6,.6},.08,.07);},"outside slab");
    auto invalid=domain;invalid.cellSize.x=0;
    rejects([&]{field.encloseInteriorRk2(invalid,{.5,.5},{.6,.6},0,.1);},"invalid finite");
    invalid=domain;invalid.cellSize={2,2};rejects([&]{field.encloseInteriorRk2(invalid,{.5,.5},{.6,.6},0,.1);},"cell size exceeds");
    const double tiny=std::numeric_limits<double>::denorm_min();
    Data reversing;reversing.previous[4].x=.1;reversing.arrival[4].x=-.1;const Slab subnormal(0,.1,reversing);
    const auto tinyBounds=subnormal.encloseInteriorRk2(domain,{.5,.5},{.6,.6},0,tiny);
    const auto tinyTrace=LevelSetBoundaryTrace::traceWithJacobian(domain,{.5,.5},tiny,
        [&](Vector2D p){return subnormal.evaluate(p,tiny);},[&](Vector2D p){return subnormal.evaluate(p,0);});
    check(!tinyTrace.trajectory.boundaryHit && tinyTrace.acceptedSubsteps==1,"subnormal trace branch changed");
    contained(tinyBounds,tinyTrace.trajectory.point,*tinyTrace.derivative);++traces;
    {struct Restore{int old=std::fegetround();~Restore(){std::fesetround(old);}} restore;
        check(std::fesetround(FE_UPWARD)==0,"cannot change test rounding");rejects([&]{field.encloseInteriorRk2(domain,{.5,.5},{.6,.6},0,.1);},"round-to-nearest");}
#if defined(_M_X64) || (defined(__x86_64__) && defined(__SSE2__))
    {struct Restore{unsigned old=_mm_getcsr();~Restore(){_mm_setcsr(old);}} restore;
        _mm_setcsr(restore.old|(1u<<15));rejects([&]{field.encloseInteriorRk2(domain,{.5,.5},{.6,.6},0,.1);},"gradual underflow");}
#endif
    std::cout<<"CORE_BRANCHES,zero_time,zero_field,outflow,inflow,midpoint,target,contact,coverage,axis,CFL,subnormal,environment,PASS\n";
}
void ownership(const QuadtreeMesh& mesh)
{
    const auto old=stateFor(mesh,false);const auto inventory=old.primary().phaseVolumes;const auto* scalar=&old.geometry().snapshot();
    HermiteGalsCandidate next(mesh,{3,.2,.3,Position::End},old.primary());
    const LevelSetBoundaryTrace::Domain domain{mesh.minCoord,mesh.maxCoord,{.25,.25},false};
    const Characteristic map(mesh,domain);Data data;data.previous[0].x=.1;data.arrival=data.previous;
    auto velocity=Velocity::polynomial(old,next.identity(),next.timeLayer(),data);
    auto interval=Interval(mesh,old,velocity,map,.2,.3);auto value=encloseHermiteGalsCharacteristic(mesh,interval,{.5,.5},{.6,.6});
    const auto same=value;check(&same.bounds()==&value.bounds() && same.interval().sameInterval(interval),"copy lost retained owner");
    const auto different=Velocity::polynomial(old,next.identity(),next.timeLayer(),data);
    const Interval other(mesh,old,different,map,.2,.3);const auto otherValue=encloseHermiteGalsCharacteristic(mesh,other,{.5,.5},{.6,.6});
    check(!otherValue.interval().sameInterval(value.interval()),"equal samples merged source identities");
    auto moved=std::move(value);rejects([&]{value.bounds();},"moved-from characteristic enclosure");
    auto movedInterval=std::move(interval);rejects([&]{encloseHermiteGalsCharacteristic(mesh,interval,{.5,.5},{.6,.6});},"moved-from interval");
    auto movedVelocity=std::move(velocity);rejects([&]{velocity.polynomialVelocity();},"moved-from velocity");
    check(moved.interval().sameInterval(movedInterval) && moved.interval().velocitySource().sameSource(movedVelocity),"move lost owner");
    // Existing mesh identity is exact geometry/connectivity, not an object
    // address. An independently allocated identical native mesh stays valid.
    auto identical=makeMesh();moved.validateMesh(*identical);
    auto foreign=makeMesh(true);rejects([&]{moved.validateMesh(*foreign);},"mesh");
    rejects([&]{encloseHermiteGalsCharacteristic(*foreign,movedInterval,{.5,.5},{.6,.6});},"mesh");
    const auto constant=Velocity::constant(old,next.identity(),next.timeLayer(),{.1,0});const Interval noPolynomial(mesh,old,constant,map,.2,.3);
    rejects([&]{encloseHermiteGalsCharacteristic(mesh,noPolynomial,{.5,.5},{.6,.6});},"owned polynomial");
    const auto callback=[](Vector2D){return LevelSetBoundaryTrace::VelocityJet{{.1,0},{0,0,0,0}};};
    const Velocity opaque(old,next.identity(),next.timeLayer(),callback,callback);const Interval noProof(mesh,old,opaque,map,.2,.3);
    rejects([&]{encloseHermiteGalsCharacteristic(mesh,noProof,{.5,.5},{.6,.6});},"owned polynomial");
    std::optional<Enclosure> retained;
    {const auto localOld=stateFor(mesh,false);HermiteGalsCandidate localNext(mesh,next.timeLayer(),localOld.primary());
        const auto localVelocity=Velocity::polynomial(localOld,localNext.identity(),localNext.timeLayer(),data);
        const Characteristic localMap(mesh,domain);const Interval local(mesh,localOld,localVelocity,localMap,.2,.3);
        retained=encloseHermiteGalsCharacteristic(mesh,local,{.5,.5},{.6,.6});}
    retained->validateMesh(mesh);contained(.49,{retained->bounds().departureLower.x,retained->bounds().departureUpper.x},"retained bounds lost lifetime");
    const auto cellBounds=retained->interval().previousAccepted().geometry().snapshot().nativeCellBounds(0);
    check(cellBounds[0].x<cellBounds[1].x && cellBounds[0].y<cellBounds[1].y,"retained geometry lost lifetime");
    check(&old.geometry().snapshot()==scalar,"failed query changed old geometry");
    for(std::size_t i=0;i<inventory.size();++i)check(old.primary().phaseVolumes[i].liquid==inventory[i].liquid,"failed query changed M");
    // Same mesh object after refinement must fail just like a foreign object.
    auto changing=makeMesh();const auto prior=stateFor(*changing,false);HermiteGalsCandidate after(*changing,next.timeLayer(),prior.primary());
    const auto source=Velocity::polynomial(prior,after.identity(),after.timeLayer(),data);const Characteristic c(*changing,domain);
    const Interval i(*changing,prior,source,c,.2,.3);const auto stale=encloseHermiteGalsCharacteristic(*changing,i,{.5,.5},{.6,.6});
    DictionaryMesh::RefineBox box;box.level=1;box.min={.25,.25};box.max={.75,.75};changing->applyRefineBoxes({box});
    rejects([&]{stale.validateMesh(*changing);},"mesh");
    std::cout<<"OWNER_AND_FAILURE_IMMUTABILITY,PASS rejections="<<rejections<<'\n';
}
}
int main() try
{
    std::cout<<std::scientific<<std::setprecision(16);
    const std::array<Flow,8> flows{{{0,.8,.8,false},{0,.8,-.6,false},{1,.8,.8,false},{1,.8,-.6,false},
        {2,.3,.3,false},{3,.4,-.2,true},{4,.1,.1,true},{5,0,0,false}}};
    for(bool amr:{false,true}){auto mesh=makeMesh(amr);for(int i=0;i<int(flows.size());++i){runBox(*mesh,amr,flows[i],i);
            if(flows[i].radial)runBox(*mesh,amr,flows[i],i,true);}}
    coreBranches();auto mesh=makeMesh();ownership(*mesh);
    check(boxes==20 && traces==4044 && identityChecks==500 && rejections>=24,"incomplete enclosure matrix");
    std::cout<<"I12a13 characteristic enclosure PASS boxes="<<boxes<<" traces="<<traces<<" identity_checks="<<identityChecks
        <<" rejections="<<rejections<<" max_affine_error="<<maximumAnalyticError<<" max_L1_travel="<<maximumTravel<<'\n';
}catch(const std::exception& e){std::cerr<<"I12a13 characteristic enclosure FAIL: "<<e.what()<<'\n';return 1;}
