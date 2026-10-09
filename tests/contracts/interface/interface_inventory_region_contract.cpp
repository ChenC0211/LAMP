#include "Modules/Interface/HermiteGalsInventoryRegion.h"
#include "Modules/Interface/HermiteGalsAmrClosure.h"

#include "Modules/Interface/HermiteGalsPhaseVolume.h"
#include "Modules/Interface/InterfaceCellMeasure.h"
#include "Core/Geometry/ActiveQuarticHermitePatch.h"
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace
{
using namespace Quantum::Interface;
using Patch=Quantum::Geometry::ActiveQuarticHermitePatch;
using Jet=Patch::NodalDof;
using Profile=std::function<Jet(Vector2D)>;
using Key=std::array<int,2>;
using Position=HermiteGalsTimeLayer::Position;
constexpr long double pi=3.1415926535897932384626433832795L;
void check(bool ok,const std::string& why){if(!ok)throw std::runtime_error(why);}
long double close(long double a,long double b,long double tolerance,const std::string& why)
{const auto error=std::abs(a-b);check(std::isfinite(a) && error<=tolerance,why+" error="+std::to_string(static_cast<double>(error)));return error;}
template<class F> void rejects(F f,const std::string& expected)
{try{f();}catch(const std::exception& e){check(std::string(e.what()).find(expected)!=std::string::npos,std::string("unexpected rejection: ")+e.what());std::cout<<"REJECT,"<<expected<<'\n';return;}throw std::runtime_error("missing rejection: "+expected);}
std::unique_ptr<QuadtreeMesh> makeMesh(int n,int amr=0)
{DictionaryMesh d;d.Nx=d.Ny=n;d.maxCellNum=16384;d.maxInitialIterations=0;auto mesh=std::make_unique<QuadtreeMesh>(d);
    if(amr>1){DictionaryMesh::RefineBox b;b.level=1;b.min={.125,.125};b.max={.875,.875};mesh->applyRefineBoxes({b});}
    if(amr){DictionaryMesh::RefineBox b;b.level=amr;b.min={.25,.25};b.max={.75,.75};mesh->applyRefineBoxes({b});}return mesh;}
Vector2D point(const QuadtreeMesh& mesh,int vertex){return mesh.verts(mesh.vertIndexToLabel[vertex]).coord();}
NativeLeafDomain domainFor(const QuadtreeMesh& mesh,bool radial)
{FaceSegmentGraph graph;graph.build(mesh);return NativeLeafDomain::build(mesh,graph,radial);}
bool fineSide(const QuadtreeMesh& m,Key key)
{const auto a=point(m,key[0]),b=point(m,key[1]);const bool h=a.y==b.y;
    for(int v:key){const auto owner=m.verts(m.vertIndexToLabel[v]).supCell;if(owner<0)continue;
        const auto& c=m.cells(owner);const auto lo=m.verts(c.cellToVert[0]).coord(),hi=m.verts(c.cellToVert[2]).coord();
        if(h && (a.y==lo.y || a.y==hi.y) && a.x>=lo.x && b.x<=hi.x && b.x-a.x<hi.x-lo.x)return true;
        if(!h && (a.x==lo.x || a.x==hi.x) && a.y>=lo.y && b.y<=hi.y && b.y-a.y<hi.y-lo.y)return true;}return false;}
using Inventory=std::function<double(int,Vector2D,Vector2D)>;
// Operator fixture with explicit primary M. Q4 geometry uses I06 mixed
// closure and analytic center samples; this does not implement I11.
HermiteGalsState stateFor(const QuadtreeMesh& mesh,const Profile& f,bool radial,HermiteGalsTimeLayer time,const Inventory& inventory,bool axisSymmetry=true)
{
    HermiteGalsPrimaryState p;p.measure=radial?HermiteGalsMeasure::AxisymmetricVolume:HermiteGalsMeasure::CartesianArea;
    Patch patch;patch.activeCells.assign(mesh.cellIndexNum(),1);patch.nodes.resize(mesh.vertIndexNum());patch.cellCenters.resize(mesh.cellIndexNum());
    for(int v=0;v<mesh.vertIndexNum();++v)if(mesh.verts(mesh.vertIndexToLabel[v]).type()==TVertType::non){const auto j=f(point(mesh,v));patch.nodes[v]=j;p.vertices.push_back({v,j.value,j.dx,j.dy});}
    std::set<Key> keys;
    for(int c=0;c<mesh.cellIndexNum();++c){const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto a=mesh.verts(cell.cellToVert[0]).coord(),b=mesh.verts(cell.cellToVert[2]).coord();
        p.geometryCells.push_back(c);p.phaseVolumes.push_back({c,inventory(c,a,b)});patch.cellCenters[c]=f(a+(b-a)*.5).value;
        for(int k=0;k<4;++k){Key key{{mesh.verts(cell.cellToVert[k]).index(),mesh.verts(cell.cellToVert[(k+1)%4]).index()}};
            const auto x=point(mesh,key[0]),y=point(mesh,key[1]);if(x.y==y.y?x.x>y.x:x.y>y.y)std::swap(key[0],key[1]);keys.insert(key);}}
    for(const auto& key:keys)if(!fineSide(mesh,key)){const auto a=point(mesh,key[0]),b=point(mesh,key[1]);const auto j=f(a+(b-a)*.5);const auto dn=a.y==b.y?j.dy:j.dx;
        p.edges.push_back({key,j.value,dn});patch.edgeMidpoints.push_back({key,j.value,dn});}
    HermiteGalsCandidate candidate(mesh,time,p);const auto mixed=closeHermiteGalsMixedAmr(mesh,candidate);
    HermiteGalsClosure closure;closure.vertexMixedDerivatives=mixed.vertexMixedDerivatives();closure.cellBubbleCoefficients=patch.cellCenters;
    for(std::size_t k=0;k<p.vertices.size();++k)patch.nodes[p.vertices[k].vertex].dxy=closure.vertexMixedDerivatives[k];
    candidate.bindGeometry(mesh,candidate.identity(),closure,LevelSetHermiteGeometry::Snapshot::fromNativeLeaves(mesh,patch,radial && axisSymmetry));return std::move(candidate).finish();
}
using Snapshot=LevelSetHermiteGeometry::Snapshot;
using Interval=Snapshot::QuarticCenterInterval;
Patch patchFor(const QuadtreeMesh& mesh,const HermiteGalsCandidate& candidate,const HermiteGalsAmrClosureResult& mixed,const std::vector<double>& c)
{
    Patch patch;const auto& p=candidate.primary();patch.activeCells.assign(mesh.cellIndexNum(),0);
    patch.nodes.resize(mesh.vertIndexNum());patch.cellCenters.resize(mesh.cellIndexNum());
    for(std::size_t k=0;k<p.vertices.size();++k){const auto& v=p.vertices[k];patch.nodes[v.vertex]={v.phi,v.dx,v.dy,mixed.vertexMixedDerivatives()[k]};}
    for(const auto& e:p.edges)patch.edgeMidpoints.push_back({e.vertices,e.phi,e.dn});
    for(std::size_t k=0;k<p.geometryCells.size();++k){patch.activeCells[p.geometryCells[k]]=1;patch.cellCenters[p.geometryCells[k]]=c[k];}
    return patch;
}
std::shared_ptr<const Snapshot> scalarFor(const QuadtreeMesh& mesh,const Patch& p,bool radial)
{return std::make_shared<const Snapshot>(Snapshot::fromNativeLeaves(mesh,p,radial));}
std::vector<double> centers(const QuadtreeMesh& mesh,const Profile& f)
{std::vector<double> c;for(int k=0;k<mesh.cellIndexNum();++k)c.push_back(f(mesh.cells(mesh.cellIndexToLabel[k]).centerCoord()).value);return c;}

constexpr long double eps=std::numeric_limits<double>::epsilon(),roundoff=64*eps;
int cases=0,rejections=0;std::size_t queries=0;
static_assert(!std::is_default_constructible_v<HermiteGalsInventoryRegionSource>);
static_assert(std::is_same_v<decltype(std::declval<HermiteGalsInventoryRegionSource&>().cells()),const std::vector<HermiteGalsInventoryRegionCell>&>);
struct CaseResult{long double volume=0,capacity=0;};
CaseResult curved(bool radial,bool amr,int sign,int perturb)
{
    auto mesh=makeMesh(4,amr?1:0);const auto domain=domainFor(*mesh,radial);
    const Profile profile=[=](Vector2D p){const double x=p.x-(radial?0:.4375),y=p.y-.53125,r=radial?.28125:.1875;
        return Jet{sign*(x*x+y*y-r*r),sign*2*x,sign*2*y,0};};
    const HermiteGalsTimeLayer time{7,.2,.21,Position::End};
    auto seed=stateFor(*mesh,profile,radial,time,[](int,Vector2D,Vector2D){return 0.;});
    const HermiteGalsPhaseRegionSource geometry(*mesh,seed);std::vector<double> inventory;
    for(int c=0;c<mesh->cellIndexNum();++c){const auto b=seed.geometry().snapshot().nativeCellBounds(c);
        const auto r=geometry.measure(*mesh,b[0],b[1]);const double cap=domain.cellVolumes()[c];
        double m=r.liquid==0?0:r.gas==0?cap:static_cast<double>(r.liquid);
        if(perturb){if(r.liquid==0)m=4e-11*cap;else if(r.gas==0)m=cap-4e-11*cap;else m+=perturb*4e-11*cap;}
        check(m>=0 && m<=cap,"fixture inventory out of bounds");inventory.push_back(m);
    }
    const auto source=stateFor(*mesh,profile,radial,time,[&](int c,Vector2D,Vector2D){return inventory[c];});
    const auto id=source.identity();const auto* scalar=&source.geometry().snapshot();
    const HermiteGalsPhaseRegionSource regions(*mesh,source);HermiteGalsInventoryRegionSource prepared(*mesh,regions);
    check(prepared.report().relativeGeometryInventoryTolerance==1e-10,"geometric budget not reported");
    check(&prepared.source().geometry().snapshot()==scalar && &prepared.regions().source().geometry().snapshot()==scalar
        && prepared.source().identity()==id && prepared.source().timeLayer().time()==time.time(),"prepared source changed");
    const auto full=prepared.measure(*mesh,{0,0},{1,1});++queries;
    long double expected=0,maxRecovery=0,maxPartition=0,partInventory=0;
    for(const auto& c:prepared.cells()){
        check(c.gasWeight>=0 && c.gasWeight<=c.liquidWeight && c.liquidWeight<=1,"weights not bounded");
        check(source.primary().phaseVolumes[c.cell].liquid==inventory[c.cell],"source M changed");expected+=inventory[c.cell];
        const auto b=scalar->nativeCellBounds(c.cell);const auto r=prepared.measure(*mesh,b[0],b[1]);++queries;
        check(r.liquid>=0 && r.gas>=0 && r.liquid<=r.geometry.capacity && r.gas<=r.geometry.capacity,"query out of strict bounds");
        maxRecovery=std::max(maxRecovery,std::abs(r.liquid-c.inventory)/c.capacity);
        const double mid=.5*b[0].x+.5*b[1].x;
        const auto left=prepared.measure(*mesh,b[0],{mid,b[1].y}),right=prepared.measure(*mesh,{mid,b[0].y},b[1]);queries+=2;
        maxPartition=std::max(maxPartition,std::abs(left.liquid+right.liquid-c.inventory)/c.capacity);
        close(left.l1Correction+right.l1Correction,std::abs(c.representedInventory-c.geometricLiquid),roundoff*c.capacity,"partition L1 correction");
        partInventory+=r.liquid;
    }
    check(maxRecovery<=roundoff && maxPartition<=roundoff,"whole/split inventory recovery exceeds roundoff");
    close(full.liquid,expected,roundoff*full.geometry.capacity,"global inventory recovery");
    close(full.liquid,partInventory,roundoff*full.geometry.capacity,"global leaf additivity");
    const Vector2D lo{.125,.1875},hi{.875,.8125};const auto query=prepared.measure(*mesh,lo,hi);++queries;
    const auto raw=regions.measure(*mesh,lo,hi);
    check(query.geometry.liquid==raw.liquid && query.geometry.gas==raw.gas,"old region geometry changed");
    check(query.liquid>=0 && query.gas>=0 && query.liquid<=query.geometry.capacity && query.gas<=query.geometry.capacity,"cross-leaf inventory not bounded");
    check(std::abs(query.signedCorrection)<=query.l1Correction+roundoff*query.geometry.capacity,"signed correction exceeds L1");
    check(query.l1Correction<=prepared.report().geometryInventoryL1+roundoff*full.geometry.capacity,"subregion invents correction mass");
    if(perturb)check(prepared.report().maxRelativeGeometryDefect>1000*eps,"manufactured non-roundoff gap disappeared");
    const auto repeated=prepared.measure(*mesh,lo,hi);check(repeated.liquid==query.liquid && repeated.gas==query.gas,"repeat query differs");++queries;
    auto copied=prepared;auto moved=std::move(copied);check(moved.measure(*mesh,lo,hi).liquid==query.liquid,"copy/move changed source");++queries;
    if(!radial && !amr && sign==1 && perturb==1){
        const auto reject=[&](const auto& f,const std::string& reason){rejects(f,reason);++rejections;};
        reject([&]{copied.report();},"moved-from");auto other=makeMesh(8);
        reject([&]{prepared.measure(*other,lo,hi);},"mesh");
        reject([&]{prepared.measure(*mesh,{-1,0},{.5,.5});},"leaves active scalar support");
        reject([&]{prepared.measure(*mesh,{.5,.5},{.5,.75});},"invalid rectangle bounds");
        reject([&]{HermiteGalsInventoryRegionSource x(*mesh,regions,{0});},"invalid geometric inventory tolerance");
        reject([&]{HermiteGalsInventoryRegionSource x(*mesh,regions,{std::numeric_limits<double>::quiet_NaN()});},"invalid geometric inventory tolerance");
        reject([&]{HermiteGalsInventoryRegionSource x(*mesh,regions,{1e-12});},"gap exceeds configured tolerance");
        const auto over=stateFor(*mesh,profile,radial,time,[&](int c,Vector2D,Vector2D){return c==0?2*domain.cellVolumes()[c]:inventory[c];});
        const HermiteGalsPhaseRegionSource overRegions(*mesh,over);
        reject([&]{HermiteGalsInventoryRegionSource x(*mesh,overRegions);},"outside strict capacity");
    }
    std::cout<<"INVENTORY_REGION,"<<radial<<','<<amr<<','<<sign<<','<<perturb<<','<<prepared.report().maxRelativeGeometryDefect
        <<','<<maxRecovery<<','<<maxPartition<<','<<prepared.report().maxRelativeCapacityDefect<<','<<query.signedCorrection<<'\n';++cases;
    return {query.liquid,query.geometry.capacity};
}
void constantReference(bool radial,bool amr,int phase)
{
    auto mesh=makeMesh(4,amr?1:0);const auto domain=domainFor(*mesh,radial);
    const Profile profile=[=](Vector2D){return Jet{double(phase),0,0,0};};
    const auto fraction=[](int c){return .125+.25*(c%4);};
    const auto source=stateFor(*mesh,profile,radial,{3,.1,.2,Position::End},[&](int c,Vector2D,Vector2D){return fraction(c)*domain.cellVolumes()[c];});
    const HermiteGalsPhaseRegionSource geometry(*mesh,source);
    const HermiteGalsInventoryRegionSource prepared(*mesh,geometry,{.95}); // explicit algebra stress, not a production geometry tolerance
    const Vector2D lo{.125,.1875},hi{.8125,.875};const auto query=prepared.measure(*mesh,lo,hi);++queries;
    long double reference=0;
    // Independent integration of the manufactured piecewise constant field.
    for(int c=0;c<mesh->cellIndexNum();++c){const auto b=source.geometry().snapshot().nativeCellBounds(c);
        const long double x0=std::max(lo.x,b[0].x),x1=std::min(hi.x,b[1].x),y0=std::max(lo.y,b[0].y),y1=std::min(hi.y,b[1].y);
        if(x0>=x1 || y0>=y1)continue;
        reference+=fraction(c)*(y1-y0)*(radial?pi*(x1*x1-x0*x0):x1-x0);
    }
    const auto error=std::abs(query.liquid-reference)/query.geometry.capacity;
    check(error<=16*eps && query.liquid>=0 && query.liquid<=query.geometry.capacity,"independent constant-field reference failed");
    std::cout<<"INVENTORY_CONSTANT,"<<radial<<','<<amr<<','<<phase<<','<<error<<'\n';++cases;
}
void ownershipAndIdentity()
{
    std::optional<HermiteGalsInventoryRegionSource> retained;
    {
        auto mesh=makeMesh(4);const Profile pure=[](Vector2D){return Jet{-1,0,0,0};};
        const auto source=stateFor(*mesh,pure,false,{},[](int,Vector2D lo,Vector2D hi){return (hi.x-lo.x)*(hi.y-lo.y);});
        const HermiteGalsPhaseRegionSource regions(*mesh,source);retained.emplace(*mesh,regions);
        const auto q=retained->measure(*mesh,{.125,.1875},{.875,.8125});++queries;
        check(q.liquid==q.geometry.liquid && q.gas==q.geometry.gas && q.signedCorrection==0,"identity changed geometry integral");
    }
    check(retained->cells().size()==16 && retained->report().inventory==1
        && retained->source().geometry().snapshot().scalarJetInNativeCell(0,{.125,.125}).value==-1,"prepared source lifetime failed");
    ++cases;std::cout<<"INVENTORY_OWNERSHIP,PASS\n";
}
void endpointWeights(bool radial,bool amr)
{
    auto mesh=makeMesh(4,amr?1:0);const auto domain=domainFor(*mesh,radial);
    const Profile profile=[=](Vector2D p){const double x=p.x-(radial?0:.4375),y=p.y-.53125,r=radial?.28125:.1875;
        return Jet{x*x+y*y-r*r,2*x,2*y,0};};
    const auto seed=stateFor(*mesh,profile,radial,{},[](int,Vector2D,Vector2D){return 0.;});
    const HermiteGalsPhaseRegionSource geometry(*mesh,seed);std::vector<double> inventory;int selected=-1;
    for(int c=0;c<mesh->cellIndexNum();++c){const auto b=seed.geometry().snapshot().nativeCellBounds(c);const auto q=geometry.measure(*mesh,b[0],b[1]);
        const double cap=domain.cellVolumes()[c],m=q.liquid==0?0:q.gas==0?cap:static_cast<double>(q.liquid);inventory.push_back(m);
        if(m>.15*cap && m<.85*cap)selected=c;
    }
    check(selected>=0,"missing endpoint cut");const auto b=seed.geometry().snapshot().nativeCellBounds(selected);
    for(bool full:{false,true}){inventory[selected]=full?domain.cellVolumes()[selected]:0;
        const auto source=stateFor(*mesh,profile,radial,{7,.2,.21,Position::End},[&](int c,Vector2D,Vector2D){return inventory[c];});
        const HermiteGalsPhaseRegionSource regions(*mesh,source);const HermiteGalsInventoryRegionSource prepared(*mesh,regions,{.95});
        const auto q=prepared.measure(*mesh,b[0],{.5*b[0].x+.5*b[1].x,b[1].y});++queries;
        check(q.liquid==(full?q.geometry.capacity:0) && q.gas==(full?0:q.geometry.capacity),"cut inventory endpoint not exact");
        check(source.primary().phaseVolumes[selected].liquid==inventory[selected],"endpoint M was clipped");
        std::cout<<"INVENTORY_ENDPOINT,"<<radial<<','<<amr<<','<<full<<','<<q.liquid<<','<<q.geometry.capacity<<'\n';++cases;
    }
}
void capacityNormalization(int phase)
{
    DictionaryMesh d;d.Nx=d.Ny=4;d.minCoord={.1,.2};d.maxCoord={1.3,1.7};d.maxCellNum=4096;d.maxInitialIterations=0;
    QuadtreeMesh mesh(d);const auto domain=domainFor(mesh,true);
    const Profile profile=[=](Vector2D){return Jet{double(phase),0,0,0};};
    const auto source=stateFor(mesh,profile,true,{7,.2,.21,Position::End},[&](int c,Vector2D,Vector2D){return .375*domain.cellVolumes()[c];},false);
    const HermiteGalsPhaseRegionSource geometry(mesh,source);const HermiteGalsInventoryRegionSource prepared(mesh,geometry,{.95});
    long double recovery=0;
    for(const auto& c:prepared.cells()){const auto b=source.geometry().snapshot().nativeCellBounds(c.cell);const auto q=prepared.measure(mesh,b[0],b[1]);++queries;
        recovery=std::max(recovery,std::abs(q.liquid-c.inventory)/c.capacity);}
    const Vector2D lo{.175,.33125},hi{1.175,1.45625};const auto q=prepared.measure(mesh,lo,hi);++queries;
    const long double reference=.375L*pi*((long double)hi.x*hi.x-(long double)lo.x*lo.x)*(hi.y-lo.y);
    const auto error=std::abs(q.liquid-reference)/q.geometry.capacity;
    check(prepared.report().maxRelativeCapacityDefect>0 && prepared.report().maxRelativeCapacityDefect<=roundoff
        && recovery<=roundoff && error<=16*eps,"nonbinary capacity normalization failed");
    std::cout<<"INVENTORY_CAPACITY,"<<phase<<','<<prepared.report().maxRelativeCapacityDefect<<','<<recovery<<','<<error<<'\n';++cases;
}
}
int main()try
{
    std::cout<<std::setprecision(17);
    for(bool radial:{false,true})for(bool amr:{false,true})for(int perturb:{-1,0,1}){
        const auto a=curved(radial,amr,1,perturb),b=curved(radial,amr,-1,-perturb);
        close(a.volume+b.volume,a.capacity,roundoff*a.capacity,"phase complement symmetry");
    }
    for(bool radial:{false,true})for(bool amr:{false,true})for(int phase:{-1,1})constantReference(radial,amr,phase);
    ownershipAndIdentity();for(bool radial:{false,true})for(bool amr:{false,true})endpointWeights(radial,amr);
    for(int phase:{-1,1})capacityNormalization(phase);
    check(cases==43 && rejections==8,"incomplete inventory region coverage");
    std::cout<<"I12a11r18 inventory region PASS cases="<<cases<<" rejections="<<rejections<<" queries="<<queries<<'\n';return 0;
}catch(const std::exception& e){std::cerr<<"inventory region FAIL: "<<e.what()<<'\n';return 1;}
