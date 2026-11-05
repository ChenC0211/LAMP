#include "Modules/Interface/HermiteGalsReflux.h"
#include "Modules/Interface/HermiteGalsTemporalPhaseFlux.h"
#include "Modules/Interface/HermiteGalsAmrClosure.h"
#include "Core/Geometry/ActiveQuarticHermitePatch.h"
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace
{
using namespace LAMP::Interface;
using Patch=LAMP::Geometry::ActiveQuarticHermitePatch;
using Jet=Patch::NodalDof;
using Profile=std::function<Jet(Vector2D)>;
using Key=std::array<int,2>;
using Position=HermiteGalsTimeLayer::Position;
using Velocity=HermiteGalsTemporalVelocitySource;
using Characteristic=HermiteGalsTemporalCharacteristicSource;
using Interval=HermiteGalsPhaseFluxIntervalSource;
using Transfer=HermiteGalsIntegratedFaceTransfer;
using Volumes=HermiteGalsIntegratedFaceVolumes;
using Orientation=HermiteGalsFaceTransferOrientation;
using Reflux=HermiteGalsIntegratedRefluxResult;
static_assert(!std::is_default_constructible_v<Reflux>);
static_assert(std::is_same_v<decltype(std::declval<Reflux&>().source()),const HermiteGalsIntegratedInventoryResult&>);
int cases=0,rejections=0;
void check(bool ok,const std::string& why){if(!ok)throw std::runtime_error(why);}
void near(long double a,long double b,long double scale,const char* why)
{check(std::abs(a-b)<=64*std::numeric_limits<double>::epsilon()*scale,why);}
template<class F>void rejects(F f,const char* why)
{
    try{f();}catch(const std::exception& e){check(std::string(e.what()).find(why)!=std::string::npos,e.what());++rejections;return;}
    throw std::runtime_error(std::string("missing rejection: ")+why);
}
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
HermiteGalsState stateFor(const QuadtreeMesh& mesh,const Profile& f,bool radial,HermiteGalsTimeLayer time,const Inventory& inventory)
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
    candidate.bindGeometry(mesh,candidate.identity(),closure,LevelSetHermiteGeometry::Snapshot::fromNativeLeaves(mesh,patch,radial));return std::move(candidate).finish();
}
LevelSetBoundaryTrace::Domain traceDomain(const QuadtreeMesh& mesh,bool radial)
{return {mesh.minCoord,mesh.maxCoord,{.125,.125},radial};}
double fraction(Vector2D lo,Vector2D hi,bool horizontal,int sign,double position)
{
    const double a=horizontal?lo.y:lo.x,b=horizontal?hi.y:hi.x;
    return sign>0?std::clamp((position-a)/(b-a),0.,1.):std::clamp((b-position)/(b-a),0.,1.);
}
HermiteGalsState plane(const QuadtreeMesh& mesh,bool radial,int sign,double position,bool empty=false)
{
    const auto domain=domainFor(mesh,radial);
    const Profile profile=[=](Vector2D x){return Jet{sign*((radial?x.y:x.x)-position),
        radial?0.:double(sign),radial?double(sign):0.,0};};
    return stateFor(mesh,profile,radial,{},[&](int cell,Vector2D lo,Vector2D hi){
        return empty?0.:domain.cellVolumes()[domain.cellRowByIndex()[cell]]*fraction(lo,hi,radial,sign,position);});
}
Transfer with(const QuadtreeMesh& mesh,const Transfer& t,const Volumes& v)
{return {mesh,t.interval(),t.face().nativeRow,t.orientation(),v};}
std::vector<Transfer> zeros(const QuadtreeMesh& mesh,const Interval& interval)
{
    const auto domain=domainFor(mesh,interval.measure()==HermiteGalsMeasure::AxisymmetricVolume);
    std::vector<Transfer> result;
    for(const auto f:domain.faces(mesh))result.emplace_back(mesh,interval,f.nativeRow,f.horizontal?Orientation::PositiveY:Orientation::PositiveX,Volumes{});
    return result;
}
std::vector<Transfer> crossingFaces(const QuadtreeMesh& mesh,const Interval& interval,bool radial,int sign,double position)
{
    auto result=zeros(mesh,interval);
    for(auto& t:result){
        const auto& f=t.face();if(f.horizontal!=radial)continue;
        if(f.negativeCell>=0 && f.positiveCell>=0){t=integrateHermiteGalsPlanarFace(mesh,interval,f.nativeRow);continue;}
        // Explicit analytic boundary HISTORY for this manufactured test.
        // I12a2 deliberately supplies no exterior extrapolation/history.
        const double coordinate=radial?f.center.y:f.center.x;
        const long double crossing=sign*(static_cast<long double>(coordinate)-position);
        const long double duration=interval.dt();
        const long double wet=std::max(0.L,static_cast<long double>(interval.timeEnd())-
            std::max(static_cast<long double>(interval.timeBegin()),crossing));
        const long double total=f.area*duration,liquid=f.area*wet;
        t=with(mesh,t,sign>0?Volumes{total,0,liquid,0}:Volumes{0,total,0,liquid});
    }
    return result;
}

void verify(const QuadtreeMesh& mesh,const HermiteGalsIntegratedInventoryResult& input,const Reflux& result,int level)
{
    result.validateMesh(mesh);
    check(&result.source().faces()==&input.faces() && &result.source().rawFaces()==&input.rawFaces() &&
        &result.source().phaseVolumes()==&input.phaseVolumes(),"duplicated source face/inventory owner");
    check(result.source().interval().sameInterval(input.interval()),"detached interval");
    check(&result.source().interval().previousAccepted().geometry().snapshot()==
        &input.interval().previousAccepted().geometry().snapshot(),"duplicated old geometry");
    const auto& faces=input.faces();check(result.sides().size()==4*input.phaseVolumes().size(),"missing sides");
    std::vector<int> incidence(faces.size(),0);std::vector<long double> exchange(mesh.cellIndexNum(),0);
    std::size_t coarse=0;int directions=0;
    for(std::size_t k=0;k<result.sides().size();++k){
        const auto& side=result.sides()[k];const auto orientation=static_cast<int>(side.side);
        check(side.cell==input.phaseVolumes()[k/4].cell && orientation==int(k%4),"side order");
        check(side.outwardSign==(orientation<2?-1:1),"side outward sign");
        // Independent oracle: actual owners/directions, no FaceSegmentGraph.
        std::vector<std::size_t> members;
        for(std::size_t j=0;j<faces.size();++j){const auto f=faces[j].face();
            if(f.horizontal!=(orientation==1 || orientation==3))continue;
            if((orientation<2?f.positiveCell:f.negativeCell)==side.cell)members.push_back(j);}
        std::sort(members.begin(),members.end(),[&](auto x,auto y){const auto a=faces[x].face(),b=faces[y].face();
            return a.horizontal?a.center.x<b.center.x:a.center.y<b.center.y;});
        check(members.size()==side.subfaceCount,"subface membership count");
        std::array<long double,6> expected{};
        for(std::size_t j=0;j<members.size();++j){
            check(side.subfaces[j]==members[j],"subface membership/order");++incidence[members[j]];
            const auto& f=faces[members[j]];const auto& v=f.integratedVolumes();
            expected[0]+=v.integratedTotalPositiveVolume;expected[1]+=v.integratedTotalNegativeVolume;
            expected[2]+=v.integratedLiquidPositiveVolume;expected[3]+=v.integratedLiquidNegativeVolume;
            expected[4]+=f.signedTotalVolume();expected[5]+=f.signedLiquidVolume();
        }
        const std::array<long double,6> actual{{side.totalPositiveVolume,side.totalNegativeVolume,
            side.liquidPositiveVolume,side.liquidNegativeVolume,side.totalVolume,side.liquidVolume}};
        check(actual==expected,"coarse integral differs from final actual subface sum");
        check(result.averageTotalFlux(k)==side.totalVolume/input.interval().dt() &&
            result.averageLiquidFlux(k)==side.liquidVolume/input.interval().dt(),"rate not derived from actual source interval");
        if(side.subfaceCount==2){++coarse;directions|=1<<orientation;}
        exchange[side.cell]-=side.outwardSign*side.liquidVolume;
    }
    for(std::size_t k=0;k<faces.size();++k){const auto f=faces[k].face();
        check(incidence[k]==int(f.negativeCell>=0)+int(f.positiveCell>=0),"orphan/duplicate actual face");}
    const auto domain=domainFor(mesh,input.interval().measure()==HermiteGalsMeasure::AxisymmetricVolume);
    const auto& old=input.interval().previousAccepted().primary().phaseVolumes;
    for(std::size_t k=0;k<input.phaseVolumes().size();++k){const auto& v=input.phaseVolumes()[k];
        near(v.liquid,old[k].liquid+exchange[v.cell],domain.cellVolumes()[domain.cellRowByIndex()[v.cell]],"side volume ledger");
        check(v.liquid==result.source().phaseVolumes()[k].liquid,"reflux changed M");}
    check(level?(coarse>0 && directions==15):coarse==0,"missing AMR four-direction coverage");
    check(result.report().cellSides==result.sides().size() && result.report().coarseSides==coarse,"incorrect count report");
    near(result.report().maxTotalVolumeResidual,0,1,"total coarse/fine volume residual");
    near(result.report().maxLiquidVolumeResidual,0,1,"liquid coarse/fine volume residual");
}
void casesWithMesh()
{
    for(bool radial:{false,true})for(int level:{0,1,2}){
        auto mesh=makeMesh(4,level);const auto domain=domainFor(*mesh,radial);
        // Manufactured algebraic integrated inputs exercise every side,
        // counterflow and metric. They are not claimed as quadrature results.
        for(int style:{0,1,2}){
            const auto old=stateFor(*mesh,[](Vector2D){return Jet{-1,0,0,0};},radial,{},
                [&](int c,Vector2D,Vector2D){return .5*domain.cellVolumes()[domain.cellRowByIndex()[c]];});
            const double dt=.015625;const auto v=Velocity::constant(old,old.identity(),{1,0,dt,Position::End},{0,0});
            const Characteristic characteristic(*mesh,traceDomain(*mesh,radial));const Interval interval(*mesh,old,v,characteristic,0,dt);
            auto input=zeros(*mesh,interval);
            for(auto& t:input)if(t.face().negativeCell>=0 && t.face().positiveCell>=0){
                const long double total=t.face().area*dt*.125,liquid=style==0?0:style==1?total:total*(.25+.25*(t.face().nativeRow%3));
                t=with(*mesh,t,{total,total,liquid,liquid});}
            const auto inventory=assembleHermiteGalsIntegratedInventory(*mesh,interval,input);
            verify(*mesh,inventory,assembleHermiteGalsReflux(*mesh,inventory),level);++cases;
        }
        // Actual I12a2 -> I12a3 -> I09 path, with explicit boundary history.
        for(int sign:{-1,1}){
            const double dt=1./64,position=.5-sign*.75*dt;const auto old=plane(*mesh,radial,sign,position);
            const auto v=Velocity::constant(old,old.identity(),{1,0,dt,Position::End},radial?Vector2D{0,double(sign)}:Vector2D{double(sign),0});
            const Characteristic characteristic(*mesh,traceDomain(*mesh,radial));const Interval interval(*mesh,old,v,characteristic,0,dt);
            const auto inventory=assembleHermiteGalsIntegratedInventory(*mesh,interval,crossingFaces(*mesh,interval,radial,sign,position));
            verify(*mesh,inventory,assembleHermiteGalsReflux(*mesh,inventory),level);++cases;
        }
        std::cout<<"INTEGRATED_REFLUX,"<<(radial?"RZ":"XY")<<','<<level<<",PASS\n";
    }
}
void repairedAndLifetime()
{
    auto mesh=makeMesh(4,2);const auto domain=domainFor(*mesh,false);std::optional<Reflux> retained;
    {
        const auto old=stateFor(*mesh,[](Vector2D){return Jet{1,0,0,0};},false,{},[](int,Vector2D,Vector2D){return 0.;});
        const auto v=Velocity::constant(old,old.identity(),{1,0,.125,Position::End},{0,0});
        const Characteristic characteristic(*mesh,traceDomain(*mesh,false));const Interval interval(*mesh,old,v,characteristic,0,.125);
        auto input=zeros(*mesh,interval);bool found=false;
        for(auto& t:input){const auto f=t.face();if(f.negativeCell<0 || f.positiveCell<0)continue;
            if(mesh->cells(mesh->cellIndexToLabel[f.negativeCell]).level()==mesh->cells(mesh->cellIndexToLabel[f.positiveCell]).level())continue;
            const long double cap=std::min(domain.cellVolumes()[domain.cellRowByIndex()[f.negativeCell]],domain.cellVolumes()[domain.cellRowByIndex()[f.positiveCell]]);
            const long double total=cap*.01,tiny=32*std::numeric_limits<double>::epsilon()*cap;
            t=with(*mesh,t,{total,total,tiny,0});found=true;break;}
        check(found,"no actual coarse/fine repair face");
        const auto inventory=assembleHermiteGalsIntegratedInventory(*mesh,interval,input);
        check(inventory.report().repairedFaces>0,"missing repaired input");
        retained=assembleHermiteGalsReflux(*mesh,inventory);verify(*mesh,inventory,*retained,2);++cases;
        auto moved=inventory;const auto owner=std::move(moved);
        rejects([&]{assembleHermiteGalsReflux(*mesh,moved);},"moved-from integrated inventory");
        check(&owner.faces()==&inventory.faces(),"inventory copy lost owner");
    }
    auto copy=*retained;auto moved=std::move(copy);
    rejects([&]{copy.sides();},"moved-from handle");
    check(&moved.source().faces()==&retained->source().faces(),"reflux copy lost owner");
    check(retained->source().interval().previousAccepted().geometry().snapshot().nativeCellCount()==std::size_t(mesh->cellIndexNum()),
        "old geometry source lifetime expired");
    auto clone=mesh->clone();retained->validateMesh(*clone);
    rejects([&]{retained->averageTotalFlux(retained->sides().size());},"invalid vector subscript");
    rejects([&]{retained->averageLiquidFlux(retained->sides().size());},"invalid vector subscript");
    mesh->verts(mesh->vertIndexToLabel[0]).setCoord({-.125,0});
    rejects([&]{retained->validateMesh(*mesh);},"mesh");
    rejects([&]{assembleHermiteGalsReflux(*mesh,retained->source());},"mesh");++cases;
    std::cout<<"REPAIRED_SUBFACE,LIFETIME,REJECTIONS,PASS\n";
}
void averageOverflow()
{
    auto mesh=makeMesh(4);const auto domain=domainFor(*mesh,false);
    const auto old=stateFor(*mesh,[](Vector2D){return Jet{-1,0,0,0};},false,{},
        [&](int c,Vector2D,Vector2D){return .5*domain.cellVolumes()[c];});
    const double dt=std::numeric_limits<double>::denorm_min();
    const auto v=Velocity::constant(old,old.identity(),{1,0,dt,Position::End},{0,0});
    const Characteristic characteristic(*mesh,traceDomain(*mesh,false));const Interval interval(*mesh,old,v,characteristic,0,dt);
    auto input=zeros(*mesh,interval);
    // Uniform prescribed directional volume: valid CFL/divergence despite
    // an unrepresentable derived rate. No quadrature is claimed here.
    for(auto& t:input)if(!t.face().horizontal)t=with(*mesh,t,{.001,0,.0005,0});
    const auto inventory=assembleHermiteGalsIntegratedInventory(*mesh,interval,input);const auto result=assembleHermiteGalsReflux(*mesh,inventory);
    rejects([&]{result.averageTotalFlux(0);},"unrepresentable derived average rate");
    rejects([&]{result.averageLiquidFlux(0);},"unrepresentable derived average rate");++cases;
}
}
int main()
{
    try{casesWithMesh();repairedAndLifetime();averageOverflow();
        std::cout<<"I12a4 integrated reflux PASS cases="<<cases<<" rejections="<<rejections<<'\n';return 0;}
    catch(const std::exception& e){std::cerr<<"I12a4 FAIL: "<<e.what()<<'\n';return 1;}
}
