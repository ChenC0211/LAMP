#include "Modules/Interface/HermiteGalsPhaseRegion.h"
#include "Modules/Interface/HermiteGalsAmrClosure.h"
#include "Modules/Interface/HermiteGalsPhaseVolume.h"
#include "Modules/Interface/InterfaceCellMeasure.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Geometry/ActiveQuarticHermitePatch.h"
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Geometry/LevelSetMappedInterface.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <type_traits>
namespace {
using namespace Quantum::Interface;
using Patch=Quantum::Geometry::ActiveQuarticHermitePatch;
using Snapshot=LevelSetHermiteGeometry::Snapshot;
using Jet=Patch::NodalDof;
using Profile=std::function<Jet(Vector2D)>;
using Key=std::array<int,2>;
constexpr long double pi=3.1415926535897932384626433832795L;
int checks=0,rejections=0;
int sharedCutMeasures=0;
void check(bool ok,const std::string& why){if(!ok)throw std::runtime_error(why);}
template<class F> void rejects(F f,const std::string& reason){try{f();}catch(const std::exception& e){
    check(std::string(e.what()).find(reason)!=std::string::npos,std::string("unexpected rejection: ")+e.what());++rejections;return;}
    throw std::runtime_error("missing rejection: "+reason);}
std::unique_ptr<QuadtreeMesh> meshFor(int n,bool amr=false,double radialMin=0)
{
    DictionaryMesh d;d.Nx=d.Ny=n;d.minCoord={radialMin,0};d.maxCoord={1,1};d.maxCellNum=32768;d.maxInitialIterations=0;
    auto m=std::make_unique<QuadtreeMesh>(d);
    if(amr){DictionaryMesh::RefineBox b;b.level=1;b.min={radialMin,0.25};b.max={0.75,0.75};m->applyRefineBoxes({b});}
    return m;
}
Vector2D point(const QuadtreeMesh& m,int v){return m.verts(m.vertIndexToLabel[v]).coord();}
std::set<Key> sides(const QuadtreeMesh& m)
{
    std::set<Key> keys;
    for(int c=0;c<m.cellIndexNum();++c){const auto& cell=m.cells(m.cellIndexToLabel[c]);
        for(int k=0;k<4;++k){Key key{{m.verts(cell.cellToVert[k]).index(),m.verts(cell.cellToVert[(k+1)%4]).index()}};
            auto a=point(m,key[0]),b=point(m,key[1]);if(a.y==b.y?a.x>b.x:a.y>b.y)std::swap(key[0],key[1]);keys.insert(key);}}
    return keys;
}
bool fineSide(const QuadtreeMesh& m,Key key)
{
    const auto a=point(m,key[0]),b=point(m,key[1]);const bool h=a.y==b.y;
    for(int v:key){const int owner=m.verts(m.vertIndexToLabel[v]).supCell;if(owner<0)continue;
        const auto& c=m.cells(owner);const auto lo=m.verts(c.cellToVert[0]).coord(),hi=m.verts(c.cellToVert[2]).coord();
        if(h && (a.y==lo.y || a.y==hi.y) && a.x>=lo.x && b.x<=hi.x && b.x-a.x<hi.x-lo.x)return true;
        if(!h && (a.x==lo.x || a.x==hi.x) && a.y>=lo.y && b.y<=hi.y && b.y-a.y<hi.y-lo.y)return true;}
    return false;
}
// Initial GEOMETRY TEMPLATE only, never an accepted solution. Center samples
// supply the analytic Q4 interpolant. I06 supplies the actual compact mixed
// closure; I07 initializes the zero placeholder inventory. I11 remains absent.
HermiteGalsState initialTemplate(const QuadtreeMesh& mesh,const Profile& f,bool radial,
    HermiteGalsTimeLayer time={},bool mismatchedMeasure=false,bool extraInventory=false,
    HermiteGalsTransport transport=HermiteGalsTransport::ConservativeVolume)
{
    HermiteGalsPrimaryState p;p.measure=radial?HermiteGalsMeasure::AxisymmetricVolume:HermiteGalsMeasure::CartesianArea;p.transport=transport;
    Patch patch;patch.nodes.resize(mesh.vertIndexNum());patch.activeCells.assign(mesh.cellIndexNum(),1);
    patch.cellCenters.resize(mesh.cellIndexNum());
    for(int v=0;v<mesh.vertIndexNum();++v)if(mesh.verts(mesh.vertIndexToLabel[v]).type()==TVertType::non){
        const auto j=f(point(mesh,v));p.vertices.push_back({v,j.value,j.dx,j.dy});patch.nodes[v]=j;}
    for(const auto& key:sides(mesh))if(!fineSide(mesh,key)){
        const auto a=point(mesh,key[0]),b=point(mesh,key[1]);const auto j=f(a+(b-a)*.5);const auto dn=a.y==b.y?j.dy:j.dx;
        p.edges.push_back({key,j.value,dn});patch.edgeMidpoints.push_back({key,j.value,dn});}
    for(int c=0;c<mesh.cellIndexNum();++c){p.geometryCells.push_back(c);if(transport==HermiteGalsTransport::ConservativeVolume)p.phaseVolumes.push_back({c,0});
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto a=mesh.verts(cell.cellToVert[0]).coord(),b=mesh.verts(cell.cellToVert[2]).coord();
        patch.cellCenters[c]=f(a+(b-a)*.5).value;}
    HermiteGalsCandidate seed(mesh,time,p);const auto mixed=closeHermiteGalsMixedAmr(mesh,seed);
    HermiteGalsClosure closure;closure.vertexMixedDerivatives=mixed.vertexMixedDerivatives();closure.cellBubbleCoefficients=patch.cellCenters;
    for(std::size_t i=0;i<p.vertices.size();++i)patch.nodes[p.vertices[i].vertex].dxy=closure.vertexMixedDerivatives[i];
    if(extraInventory){p.geometryCells.pop_back();patch.activeCells.back()=0;closure.cellBubbleCoefficients.pop_back();
        std::set<Key> activeSides;
        for(int c:p.geometryCells){const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);
            for(int k=0;k<4;++k){Key key{{mesh.verts(cell.cellToVert[k]).index(),mesh.verts(cell.cellToVert[(k+1)%4]).index()}};
                const auto a=point(mesh,key[0]),b=point(mesh,key[1]);if(a.y==b.y?a.x>b.x:a.y>b.y)std::swap(key[0],key[1]);activeSides.insert(key);}}
        patch.edgeMidpoints.erase(std::remove_if(patch.edgeMidpoints.begin(),patch.edgeMidpoints.end(),
            [&](const auto& edge){return activeSides.count(edge.vertices)==0;}),patch.edgeMidpoints.end());}
    if(mismatchedMeasure)p.measure=radial?HermiteGalsMeasure::CartesianArea:HermiteGalsMeasure::AxisymmetricVolume;
    HermiteGalsCandidate candidate(mesh,time,p);
    candidate.bindGeometry(mesh,candidate.identity(),closure,Snapshot::fromNativeLeaves(mesh,patch,radial));
    return std::move(candidate).finish();
}
struct Circle{double x,y,r;};
Jet circleJet(Vector2D p,Circle c)
{const auto x=p.x-c.x,y=p.y-c.y;return {x*x+y*y-c.r*c.r,2*x,2*y,0};}
Jet product(Jet a,Jet b)
{return {a.value*b.value,a.dx*b.value+a.value*b.dx,a.dy*b.value+a.value*b.dy,
    a.dxy*b.value+a.dx*b.dy+a.dy*b.dx+a.value*b.dxy};}
long double diskPrimitive(long double u,long double r)
{u=std::clamp(u,-r,r);const auto s=std::sqrt(std::max(0.L,r*r-u*u));return .5L*(u*s+r*r*std::atan2(u,s));}
// Independent analytic rectangle/disk intersection: split only at exact
// circle crossings of the rectangle's horizontal bounds, then integrate
// constant + k*sqrt(R^2-(x-cx)^2). No production graph, roots or quadrature.
long double diskVolume(Vector2D lo,Vector2D hi,Circle c,bool radial)
{
    const long double cx=c.x,cy=c.y,r=c.r;
    const long double a=std::max<long double>(lo.x,cx-r),b=std::min<long double>(hi.x,cx+r);
    if(b<=a)return 0;
    std::vector<long double> points{a,b};
    for(long double y:{static_cast<long double>(lo.y),static_cast<long double>(hi.y)})
        if(std::abs(y-cy)<r){const auto u=std::sqrt(r*r-(y-cy)*(y-cy));for(auto x:{cx-u,cx+u})if(x>a && x<b)points.push_back(x);}
    std::sort(points.begin(),points.end());long double sum=0;
    const auto primitive=[&](long double x){const auto u=x-cx;
        return radial?cx*diskPrimitive(u,r)-std::pow(std::max(0.L,r*r-u*u),1.5L)/3:diskPrimitive(u,r);};
    for(std::size_t i=1;i<points.size();++i){const auto x0=points[i-1],x1=points[i],xm=(x0+x1)/2;
        const auto s=std::sqrt(std::max(0.L,r*r-(xm-cx)*(xm-cx)));
        if(std::min<long double>(hi.y,cy+s)<=std::max<long double>(lo.y,cy-s))continue;
        long double constant=0,k=0;
        if(cy+s<hi.y){constant+=cy;k+=1;}else constant+=hi.y;
        if(cy-s>lo.y){constant-=cy;k+=1;}else constant-=lo.y;
        sum+=constant*(radial?(x1*x1-x0*x0)/2:x1-x0)+k*(primitive(x1)-primitive(x0));}
    return radial?2*pi*sum:sum;
}

long double boxCapacity(Vector2D a,Vector2D b,bool radial)
{return (static_cast<long double>(b.x)-a.x)*(static_cast<long double>(b.y)-a.y)*(radial?pi*(static_cast<long double>(b.x)+a.x):1);}
void sharedWholeCells(const QuadtreeMesh& mesh,const HermiteGalsPhaseRegionSource& prepared){
    const auto& state=prepared.source();const auto id=state.identity();const auto before=state.primary().phaseVolumes;
    const auto volumes=integrateHermiteGalsPhaseVolumes(mesh,state);
    check(volumes.sourceIdentity()==id,"whole-cell measurement changed source");
    int cut=0;
    for(const auto& cell:volumes.cells())if(cell.certifiedPhase==0){
        const auto bounds=state.geometry().snapshot().nativeCellBounds(cell.cell);
        const auto region=prepared.measure(mesh,bounds[0],bounds[1]);
        check(region.intersectedLeaves==1 && region.cutLeaves==1,"whole-cell query lost cut leaf");
        check(cell.rawLiquid==region.liquid,"I07 and whole-cell region differ on same scalar");
        ++cut;++sharedCutMeasures;
    }
    check(cut>0 && state.identity()==id,"missing shared cut-cell comparison");
    for(std::size_t k=0;k<before.size();++k)check(state.primary().phaseVolumes[k].liquid==before[k].liquid,"whole-cell measurement changed M");
}
long double verify(const QuadtreeMesh& mesh,const HermiteGalsPhaseRegionSource& prepared,Vector2D a,Vector2D b,long double exact){
    const auto id=prepared.source().identity();const auto* scalar=&prepared.source().geometry().snapshot();
    const bool radial=prepared.source().primary().measure==HermiteGalsMeasure::AxisymmetricVolume;
    const auto cap=boxCapacity(a,b,radial);const auto m=prepared.measure(mesh,a,b);
    check(m.sourceIdentity==id && &prepared.source().geometry().snapshot()==scalar,"source identity changed");
    check(m.lower.x==a.x && m.lower.y==a.y && m.upper.x==b.x && m.upper.y==b.y && m.capacity==cap,"query bounds/capacity changed");
    check(m.liquid>=0 && m.liquid<=cap && m.gas>=0 && m.gas<=cap,"strict region capacity");
    check(std::abs(m.liquid-exact)<=1e-9L*cap,"analytic region error");
    check(std::abs(m.liquid+m.gas-cap)<=8*std::numeric_limits<double>::epsilon()*cap,"phase complement");
    check(std::abs(m.rawAxialLiquid-exact)<=1e-9L*cap && std::abs(m.rawRadialLiquid-exact)<=1e-9L*cap,"raw Green volume error");
    check(m.momentErrorEstimate<=1e-10L*cap && m.maxRelativeLedgerDefect<=1e-10L,"measure diagnostics");
    for(const auto& v:prepared.source().primary().phaseVolumes)check(v.liquid==0,"region measurement modified independent M");
    ++checks;return m.liquid;
}
void shapes(){
    for(bool radial:{false,true})for(bool amr:{false,true}){
        auto mesh=meshFor(4,amr);const Circle outer{radial?0:.53125,.515625,radial?.28125:.203125},inner{outer.x,outer.y,radial?.21875:.171875};
        for(bool ring:{false,true})for(int sign:{-1,1}){
            const auto profile=[&](Vector2D p){auto j=circleJet(p,outer);if(ring)j=product(j,circleJet(p,inner));
                j.value*=sign;j.dx*=sign;j.dy*=sign;j.dxy*=sign;return j;};
            auto state=initialTemplate(*mesh,profile,radial);HermiteGalsPhaseRegionSource prepared(*mesh,state);
            sharedWholeCells(*mesh,prepared);
            check(prepared.components()==(ring?2u:1u),"complete component inventory");
            const auto exact=[&](Vector2D a,Vector2D b){const auto liquid=diskVolume(a,b,outer,radial)-(ring?diskVolume(a,b,inner,radial):0);
                return sign>0?liquid:boxCapacity(a,b,radial)-liquid;};
            std::vector<std::array<Vector2D,2>> boxes{
                {{{0,0},{1,1}}},{{{.375,.375},{.6875,.65625}}},{{{.125,.125},{.875,.875}}},
                {{{.5,.56},{.5+1./1024,.625}}},{{{.4,.56},{.7,.56+1./1024}}},
                {{{.1,.1},{.2,.2}}},{{{.49,.49},{.55,.55}}}};
            if(radial)boxes={{{{0,0},{1,1}}},{{{0,.25},{.3,.75}}},{{{.125,.125},{.875,.875}}},
                {{{.125,.56},{.125+1./1024,.625}}},{{{.05,.56},{.27,.56+1./1024}}},
                {{{.1,.1},{.2,.2}}},{{{.05,.49},{.1,.55}}}};
            for(const auto& box:boxes){try{verify(*mesh,prepared,box[0],box[1],exact(box[0],box[1]));}
                catch(const std::exception& e){throw std::runtime_error("shape radial="+std::to_string(radial)+" amr="+std::to_string(amr)+" ring="+std::to_string(ring)+" sign="+std::to_string(sign)+" box="+std::to_string(box[0].x)+","+std::to_string(box[0].y)+": "+e.what());}}
            const Vector2D a{radial?.05:.3,.3},b{radial?.4:.8,.8};const auto whole=verify(*mesh,prepared,a,b,exact(a,b));
            const double splitX=radial?.2:.55;
            for(bool vertical:{false,true}){const Vector2D lowEnd=vertical?Vector2D{splitX,b.y}:Vector2D{b.x,.55};
                const Vector2D highBegin=vertical?Vector2D{splitX,a.y}:Vector2D{a.x,.55};
                const auto pieces=verify(*mesh,prepared,a,lowEnd,exact(a,lowEnd))+verify(*mesh,prepared,highBegin,b,exact(highBegin,b));
                check(std::abs(pieces-whole)<=1e-9L*boxCapacity(a,b,radial),"region partition additivity");}
        }
        for(int sign:{-1,1}){
            HermiteGalsPhaseRegionSource pure(*mesh,initialTemplate(*mesh,[&](Vector2D){return Jet{static_cast<double>(sign),0,0,0};},radial));
            for(const auto& bounds:std::vector<std::array<Vector2D,2>>{{{{.1,.2},{.9,.8}}},{{{.5,.3},{std::nextafter(.5,1.),.7}}}}){
                const auto cap=boxCapacity(bounds[0],bounds[1],radial);
                verify(*mesh,pure,bounds[0],bounds[1],sign<0?cap:0);
                const auto measured=pure.measure(*mesh,bounds[0],bounds[1]);
                check(measured.cutLeaves==0 && measured.liquid==(sign<0?cap:0),"pure phase exactness");
            }
        }
    }
}
void contracts(){
    auto mesh=meshFor(4);const Circle c{.53125,.515625,.203125};
    auto state=initialTemplate(*mesh,[&](Vector2D p){return circleJet(p,c);},false);
    HermiteGalsPhaseRegionSource source(*mesh,state),copy=source;
    check(&copy.source().geometry().snapshot()==&state.geometry().snapshot(),"prepared source copied scalar");
    static_assert(!std::is_default_constructible_v<HermiteGalsPhaseRegionSource>);
    auto clone=meshFor(4);source.validateMesh(*clone);verify(*clone,source,{.25,.25},{.75,.75},diskVolume({.25,.25},{.75,.75},c,false));
    auto foreign=meshFor(8);rejects([&]{source.measure(*foreign,{0,0},{1,1});},"mesh");
    for(const auto& b:std::vector<std::array<Vector2D,2>>{{{{0,0},{0,1}}},{{{.6,.5},{.5,.6}}},{{{0,0},{1,std::numeric_limits<double>::infinity()}}}}){
        rejects([&]{source.measure(*mesh,b[0],b[1]);},"invalid rectangle bounds");}
    rejects([&]{source.measure(*mesh,{-.01,0},{1,1});},"support");
    rejects([&]{source.measure(*mesh,{2,2},{3,3});},"support");
    HermiteGalsPhaseRegionConfig bad;bad.relativeLedgerTolerance=0;
    rejects([&]{HermiteGalsPhaseRegionSource invalid(*mesh,state,bad);},"invalid region relative tolerance");
    auto moved=std::move(copy);rejects([&]{copy.components();},"moved-from");verify(*mesh,moved,{.25,.25},{.75,.75},diskVolume({.25,.25},{.75,.75},c,false));
    const auto hole=initialTemplate(*mesh,[](Vector2D){return Jet{1,0,0,0};},false,{},false,true);
    // Existing complete discovery requires explicit contacts for a partial
    // active domain, even for a pure field. Preserve that earlier rejection.
    rejects([&]{HermiteGalsPhaseRegionSource sparse(*mesh,hole);},"active-subdomain");
    const auto plane=initialTemplate(*mesh,[](Vector2D p){return Jet{p.y-.53125,0,1,0};},false);
    HermiteGalsPhaseRegionConfig contact;contact.boundaryContacts={{0,.53125},{1,.53125}};
    HermiteGalsPhaseRegionSource planar(*mesh,plane,contact);
    verify(*mesh,planar,{.3,.4},{.7,.53125},.4L*(.53125L-.4L));
    verify(*mesh,planar,{.3,.4},{.7,.5},.4L*.1L);
    DictionaryMesh::RefineBox box;box.level=1;box.min={.25,.25};box.max={.75,.75};clone->applyRefineBoxes({box});
    rejects([&]{source.validateMesh(*clone);},"mesh");
    verify(*mesh,source,{.25,.25},{.75,.75},diskVolume({.25,.25},{.75,.75},c,false));
}
void coincidentClips(){
    using I=LevelSetMappedInterface;
    int scenarios=0,boxes=0,coreQueries=0;long double peak=0,ledgerPeak=0;
    for(bool radial:{false,true})for(bool amr:{false,true})for(bool factor:{false,true})for(int sign:{-1,1}){
        // Match the already-published H10i dual-position AMR fixture. Native
        // zero-edge publication itself remains the separate H10i1c task.
        auto mesh=meshFor(4);if(amr){DictionaryMesh::RefineBox refine;refine.level=1;
            refine.min={.25,.25};refine.max={.75,.75};mesh->applyRefineBoxes({refine});}
        const double plane=amr?.125:.375;
        const auto profile=[&](Vector2D p){
            const double normal=(radial?p.y:p.x)-plane,tangent=radial?p.x:p.y;
            const double f=factor?1+tangent*tangent:1,df=factor?2*tangent:0;
            return radial?Jet{sign*normal*f,sign*normal*df,sign*f,sign*df}:
                Jet{sign*normal*f,sign*f,sign*normal*df,sign*df};};
        auto state=initialTemplate(*mesh,profile,radial);HermiteGalsPhaseRegionConfig config;
        config.boundaryContacts=radial?std::vector<Vector2D>{{1,plane}}:std::vector<Vector2D>{{plane,0},{plane,1}};
        HermiteGalsPhaseRegionSource source(*mesh,state,config);
        const auto reference=[&](Vector2D a,Vector2D b){
            if(radial){if(sign>0)b.y=std::min(b.y,plane);else a.y=std::max(a.y,plane);}
            else{if(sign>0)b.x=std::min(b.x,plane);else a.x=std::max(a.x,plane);}
            return a.x>=b.x || a.y>=b.y?0.L:boxCapacity(a,b,radial);};
        const auto inspect=[&](Vector2D a,Vector2D b){const auto m=source.measure(*mesh,a,b);const auto expected=reference(a,b);
            check(m.liquid>=0 && m.gas>=0 && m.liquid<=m.capacity && m.gas<=m.capacity,"coincident clip strict capacity");
            for(auto v:{m.liquid,m.rawAxialLiquid,m.rawRadialLiquid})peak=std::max(peak,std::abs(v-expected)/m.capacity);
            for(auto v:{m.gas,m.rawAxialGas,m.rawRadialGas})peak=std::max(peak,std::abs(v-(m.capacity-expected))/m.capacity);
            if(expected>0)check(m.liquid>0,"lost nonzero short liquid phase");
            if(expected<m.capacity)check(m.gas>0,"lost nonzero short gas phase");
            check(peak<=1e-12L,"coincident clip closed analytic moments");++boxes;return m.liquid;};
        const auto box=[&](double nlo,double nhi,double tlo=0,double thi=1){return radial?
            std::array<Vector2D,2>{{{tlo,nlo},{thi,nhi}}}:std::array<Vector2D,2>{{{nlo,tlo},{nhi,thi}}};};
        for(const auto& b:{box(plane-.125,plane),box(plane,plane+.125),box(0,1),
            box(plane-.0625,plane+.0625,.5,.5+std::ldexp(1.,-40)),
            box(plane,plane+std::ldexp(1.,-40)),box(plane-std::ldexp(1.,-40),plane),
            box(plane-std::ldexp(1.,-40),plane+.0625),box(plane-.0625,plane+std::ldexp(1.,-40)),
            box(std::nextafter(plane,1.),plane+.0625),box(plane-.0625,std::nextafter(plane,0.))})inspect(b[0],b[1]);
        const auto whole=box(plane-.125,plane+.125),left=box(plane-.125,plane),right=box(plane,plane+.125);
        const auto sum=inspect(left[0],left[1])+inspect(right[0],right[1]);
        check(std::abs(sum-inspect(whole[0],whole[1]))<=1e-12L*boxCapacity(whole[0],whole[1],radial),"coincident split add-back");
        // Direct Core port checks reversed typed inputs and an explicit
        // fixed/exterior/coincident ledger with a nonzero common origin.
        const auto view=std::make_shared<const HermiteGalsGeometryView>(state.geometry());
        const std::shared_ptr<const Snapshot> scalar(view,&view->snapshot());
        const auto topology=LevelSetInterfacePanels::discoverNativeComponents(*mesh,scalar,radial,config.boundaryContacts,config.symmetryPlanes);
        for(int c=0;c<mesh->cellIndexNum();++c){const auto bounds=scalar->nativeCellBounds(c);const auto lo=bounds[0],hi=bounds[1];
            if(!(plane>(radial?lo.y:lo.x) && plane<(radial?hi.y:hi.x)))continue;
            const Vector2D cutLo=radial?Vector2D{lo.x,plane}:Vector2D{plane,lo.y},cutHi=radial?Vector2D{hi.x,plane}:Vector2D{plane,hi.y};
            const std::vector<I::Region> regions{{0,{{lo,{cutHi.x,lo.y},cutHi,{lo.x,cutHi.y}}}},
                {1,{{cutLo,{hi.x,cutLo.y},hi,{cutLo.x,hi.y}}}}};
            // XY first region is [lo.x,plane], RZ first is [lo.y,plane].
            std::map<std::array<double,4>,std::size_t> keys;std::vector<LevelSetMappedSegment::NativeQuarticTrace> cuts;
            for(const auto& r:regions)for(int e=0;e<4;++e){auto a=r.vertices[e],b=r.vertices[(e+1)%4];
                if(a.x>b.x || (a.x==b.x && a.y>b.y))std::swap(a,b);
                const std::array<double,4> key{{a.x,a.y,b.x,b.y}};
                if(keys.emplace(key,cuts.size()).second)cuts.push_back(LevelSetMappedSegment::nativeQuarticTrace(scalar,b,a,a.x==b.x?0:1));}
            std::vector<const LevelSetMappedSegment::NativeQuarticTrace*> refs;for(const auto& t:cuts)refs.push_back(&t);
            std::vector<I::Graph> graphs;for(const auto& chain:topology.curves)for(const auto& g:chain.graphs())if(g.nativeCellIndex==c)graphs.push_back(g);
            const Vector2D origin{-.25,-.125};
            const auto m=I::integrateNativeQuarticRegionMoments(*scalar,c,graphs,refs,regions,radial,origin);
            long double fixedA=0,fixedR=0;
            for(const auto& r:m.regions){const bool isLiquid=sign*(r.regionId?1:-1)<0;const auto expected=isLiquid?r.capacity:0;
                for(auto value:{r.axialLiquid,r.radialLiquid})check(std::abs(value-expected)<=1e-12L*r.capacity,"direct reversed Core closed phase");
                fixedA+=r.fixedAxialLiquid;fixedR+=r.fixedRadialLiquid;
                check(r.freeArcs==0 && r.freeAxialLiquid==0 && r.freeRadialLiquid==0,"coincident free arc was counted twice");}
            const auto cap=boxCapacity(lo,hi,radial);
            ledgerPeak=std::max({ledgerPeak,std::abs(fixedA-m.exteriorFixedAxialLiquid-m.coincidentAxialLiquid)/cap,
                std::abs(fixedR-m.exteriorFixedRadialLiquid-m.coincidentRadialLiquid)/cap});
            check(m.coincidentArcs>0 && ledgerPeak<=1e-12L,"explicit inside-trace coincident ledger");
            auto joinedRegions=regions;joinedRegions[1].regionId=0;
            const auto joined=I::integrateNativeQuarticRegionMoments(*scalar,c,graphs,refs,joinedRegions,radial,origin);
            check(joined.regions.size()==1 && joined.coincidentArcs==0 && joined.regions[0].freeArcs>0,
                "same-owner seam deleted its internal free interface");
            const auto expectedWhole=reference(lo,hi);
            for(auto value:{joined.regions[0].axialLiquid,joined.regions[0].radialLiquid})
                check(std::abs(value-expectedWhole)<=1e-12L*cap,"same-owner zero seam counted twice");
            check(joined.coincidentAxialLiquid==0 && joined.coincidentRadialLiquid==0,"same-owner artificial clip ledger");
            auto other=initialTemplate(*mesh,profile,radial);
            rejects([&]{I::integrateNativeQuarticRegionMoments(other.geometry().snapshot(),c,graphs,refs,regions,radial,origin);},"different snapshot");
            const auto saved=cuts.front().ordinaryRanges().size();auto moved=std::move(cuts.front());
            rejects([&]{I::integrateNativeQuarticRegionMoments(*scalar,c,graphs,refs,regions,radial,origin);},"moved-from");
            cuts.front()=std::move(moved);check(cuts.front().ordinaryRanges().size()==saved,"failed Core query mutated trace");
            auto missing=refs;missing.erase(missing.begin());
            rejects([&]{I::integrateNativeQuarticRegionMoments(*scalar,c,graphs,missing,regions,radial,origin);},"missing its actual complete trace");
            auto duplicate=refs;duplicate.push_back(refs.front());
            rejects([&]{I::integrateNativeQuarticRegionMoments(*scalar,c,graphs,duplicate,regions,radial,origin);},"duplicate actual");
            auto invalid=regions;invalid[0].vertices[2]=invalid[0].vertices[0];
            rejects([&]{I::integrateNativeQuarticRegionMoments(*scalar,c,graphs,refs,invalid,radial,origin);},"convex/cyclic");
            if(!radial && !amr && !factor && sign<0){
                const double cross=(lo.y+hi.y)*.5;
                const auto critical=initialTemplate(*mesh,[&](Vector2D p){return Jet{(p.x-plane)*(p.y-cross),p.y-cross,p.x-plane,1};},false);
                const auto criticalView=std::make_shared<const HermiteGalsGeometryView>(critical.geometry());
                const std::shared_ptr<const Snapshot> criticalScalar(criticalView,&criticalView->snapshot());
                std::vector<LevelSetMappedSegment::NativeQuarticTrace> criticalCuts;
                for(const auto& t:cuts)criticalCuts.push_back(LevelSetMappedSegment::nativeQuarticTrace(criticalScalar,t.begin(),t.end(),t.normalAxis()));
                std::vector<const LevelSetMappedSegment::NativeQuarticTrace*> criticalRefs;for(const auto& t:criticalCuts)criticalRefs.push_back(&t);
                rejects([&]{I::integrateNativeQuarticRegionMoments(*criticalScalar,c,graphs,criticalRefs,regions,false,origin);},"complete regular graph certificate");
            }
            const auto repeated=I::integrateNativeQuarticRegionMoments(*scalar,c,graphs,refs,regions,radial,origin);
            check(repeated.regions[0].axialLiquid==m.regions[0].axialLiquid && repeated.regions[1].radialLiquid==m.regions[1].radialLiquid,
                "failed native region call mutated its prepared source");
            ++coreQueries;break;
        }
        ++scenarios;
    }
    std::cout<<"COINCIDENT_CLIP,"<<scenarios<<','<<boxes<<','<<coreQueries<<','<<peak<<','<<ledgerPeak<<'\n';
}
void nativeTraceVolumes()
{
    int scenarios=0;std::size_t cells=0,halves=0,sides=0;long double peak=0;
    for(bool radial:{false,true})for(bool amr:{false,true})for(int axis:{0,1})for(int sign:{-1,1}) {
        auto mesh=meshFor(4,amr);const double line=.5;
        const Profile f=[=](Vector2D p){const double normal=axis==0?(radial?p.x*p.x-line*line:p.x-line):p.y-line;
            const double t=axis==0?p.y:p.x,factor=radial && axis==1?1+t*t:1+t+t*t,derivative=radial && axis==1?2*t:1+2*t;
            const double dn=axis==0 && radial?2*p.x:1;
            return Jet{sign*normal*factor,sign*(axis==0?dn*factor:normal*derivative),sign*(axis==1?factor:normal*derivative),sign*dn*derivative};};
        const auto state=initialTemplate(*mesh,f,radial,{},false,false,HermiteGalsTransport::FirstJet);
        check(state.primary().phaseVolumes.empty(),"trace fixture introduced inventory");
        auto owner=std::make_shared<const HermiteGalsGeometryView>(state.geometry());const std::shared_ptr<const Snapshot> scalar(owner,&owner->snapshot());
        FaceSegmentGraph graph;graph.build(*mesh);const auto domain=NativeLeafDomain::build(*mesh,graph,radial);
        const std::vector<Vector2D> contacts=axis==0?std::vector<Vector2D>{{line,0},{line,1}}:(radial?std::vector<Vector2D>{{1,line}}:std::vector<Vector2D>{{0,line},{1,line}});
        const auto trace=InterfaceCellMeasure::discoverNativeTraceVolumes(*mesh,domain,scalar,contacts);
        const auto topo=LevelSetInterfacePanels::discoverNativeTraceComponents(*mesh,scalar,radial,contacts);
        const auto cache=InterfacePhaseGeometry::MappedFacePartitions::prepareNativeTrace(*mesh,domain,scalar);
        const auto retained=InterfaceCellMeasure::measureNativeTraceVolumes(*mesh,domain,cache,topo,contacts);
        for(int cell:domain.cellIndexByRow()) {
            const auto bounds=scalar->nativeCellBounds(cell);const auto lo=bounds[0],hi=bounds[1];const long double capacity=boxCapacity(lo,hi,radial);
            auto clipped=hi;if(axis==0)clipped.x=std::max(lo.x,std::min(hi.x,line));else clipped.y=std::max(lo.y,std::min(hi.y,line));
            const long double prefix=boxCapacity(lo,clipped,radial),expected=sign>0?prefix:capacity-prefix;
            const auto axial=trace.measure.fixedBoundary.negativePhaseByCell[cell]+trace.measure.freeBoundaryMomentByCell[cell];
            const auto radialMoment=trace.ledger.fixedRadialByCell[cell]+trace.ledger.freeRadialByCell[cell];
            const auto trial=std::max({std::abs(trace.measure.negativeVolumeByCell[cell]-expected),std::abs(axial-expected),std::abs(radialMoment-expected)})/capacity;
            check(trial<=1e-12L && trace.measure.negativeVolumeByCell[cell]==retained.measure.negativeVolumeByCell[cell],"native trace actual CV two-Green reference");
            peak=std::max(peak,trial);++cells;
        }
        auto missing=topo;missing.curves.clear();rejects([&]{InterfaceCellMeasure::measureNativeTraceVolumes(*mesh,domain,cache,missing,contacts);},"complete same-source certificate");
        rejects([&]{InterfaceCellMeasure::discoverVolumes(*mesh,domain,scalar,contacts);},"coincident zero trace");
        std::weak_ptr<const Snapshot> lifetime=scalar;owner.reset();check(!lifetime.expired(),"retained trace source lifetime");
        halves+=trace.ledger.zeroHalfFaces;sides+=trace.ledger.sideGraphs;++scenarios;
    }
    check(scenarios==16 && halves>0 && sides>0,"native trace volume layout coverage");
    std::cout<<"NATIVE_TRACE_VOLUMES,"<<scenarios<<','<<cells<<','<<halves<<','<<sides<<','<<peak<<",PASS\n";
}
}
int main(int argc,char** argv){try{if(argc==2 && std::string(argv[1])=="--contracts"){contracts();coincidentClips();nativeTraceVolumes();}else{shapes();contracts();coincidentClips();nativeTraceVolumes();}
if(sharedCutMeasures)std::cout<<"SHARED_LOCAL_VOLUME,PASS,cut_cells="<<sharedCutMeasures<<'\n';
std::cout<<"I12a6 rectangular phase PASS cases="<<checks<<" rejections="<<rejections<<'\n';return 0;}
catch(const std::exception& e){std::cerr<<"I12a6 rectangular phase FAIL: "<<e.what()<<'\n';return 1;}}
