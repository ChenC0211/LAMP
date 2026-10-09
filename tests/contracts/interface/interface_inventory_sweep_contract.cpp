#include "Modules/Interface/HermiteGalsInventoryRegion.h"
#include "Modules/Interface/HermiteGalsTemporalPhaseFlux.h"
#include "Modules/Interface/HermiteGalsReflux.h"
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

using Velocity=HermiteGalsTemporalVelocitySource;
using Characteristic=HermiteGalsTemporalCharacteristicSource;
using Interval=HermiteGalsPhaseFluxIntervalSource;
using Inflow=HermiteGalsConstantInflowSource;
using Transfer=HermiteGalsIntegratedFaceTransfer;
using Phase=HermiteGalsInflowPhase;
constexpr long double eps=std::numeric_limits<double>::epsilon(),roundoff=64*eps;
int cases=0,oldRejected=0,guardRejections=0;std::size_t testedFaces=0;
std::string currentCase;
long double physicalCapacity(Vector2D a,Vector2D b,bool radial)
{return (static_cast<long double>(b.y)-a.y)*(radial?pi*((long double)b.x*b.x-(long double)a.x*a.x):(long double)b.x-a.x);}
std::vector<HermiteGalsInflowFacePhase> inflowEntries(const QuadtreeMesh& mesh,const Interval& interval,bool liquid)
{
    std::vector<HermiteGalsInflowFacePhase> entries;const auto u=*interval.velocitySource().constantValue();
    const auto domain=domainFor(mesh,interval.measure()==HermiteGalsMeasure::AxisymmetricVolume);
    for(const auto f:domain.faces(mesh)){const double speed=f.horizontal?u.y:u.x;
        if((f.negativeCell<0 && speed>0) || (f.positiveCell<0 && speed<0))entries.push_back({f.nativeRow,liquid?Phase::Liquid:Phase::Gas});}
    return entries;
}
long double departureInventory(const QuadtreeMesh& mesh,const HermiteGalsInventoryRegionSource& inventory,int cell,
    Vector2D velocity,double dt,bool inflowLiquid,bool radial)
{
    const auto b=inventory.source().geometry().snapshot().nativeCellBounds(cell);const auto lo=b[0]-velocity*dt,hi=b[1]-velocity*dt;
    const Vector2D insideLo{std::max(lo.x,mesh.minCoord.x),std::max(lo.y,mesh.minCoord.y)},
        insideHi{std::min(hi.x,mesh.maxCoord.x),std::min(hi.y,mesh.maxCoord.y)};
    long double value=0,insideCapacity=0;
    if(insideLo.x<insideHi.x && insideLo.y<insideHi.y){const auto query=inventory.measure(mesh,insideLo,insideHi);
        value=query.liquid;insideCapacity=physicalCapacity(insideLo,insideHi,radial);}
    if(inflowLiquid)value+=physicalCapacity(lo,hi,radial)-insideCapacity;
    return value;
}
void run(bool radial,bool horizontal,bool amr,int sign,int direction,bool plane,double q)
{
    currentCase=std::to_string(radial)+","+std::to_string(horizontal)+","+std::to_string(amr)+","+
        std::to_string(sign)+","+std::to_string(direction)+","+std::to_string(plane)+","+std::to_string(q);
    auto mesh=makeMesh(4,amr?1:0);const auto domain=domainFor(*mesh,radial);constexpr double dt=.09375;
    const Vector2D velocity=horizontal?Vector2D{0,double(direction)}:Vector2D{double(direction),0};
    const double position=.5-direction*q*dt;
    const Profile profile=[=](Vector2D p){if(plane)return Jet{sign*((horizontal?p.y:p.x)-position),horizontal?0.:double(sign),horizontal?double(sign):0.,0};
        const double x=p.x-(radial?0:.4375),y=p.y-.53125,r=radial?.28125:.1875;
        return Jet{sign*(x*x+y*y-r*r),sign*2*x,sign*2*y,0};};
    const auto seed=stateFor(*mesh,profile,radial,{},[](int,Vector2D,Vector2D){return 0.;});
    HermiteGalsPhaseRegionConfig regionConfig;
    if(plane){if(horizontal){regionConfig.boundaryContacts={{1,position}};if(!radial)regionConfig.boundaryContacts.push_back({0,position});}
        else regionConfig.boundaryContacts={{position,0},{position,1}};}
    const HermiteGalsPhaseRegionSource initial(*mesh,seed,regionConfig);std::vector<double> sourceM;
    for(int c=0;c<mesh->cellIndexNum();++c){const auto b=seed.geometry().snapshot().nativeCellBounds(c);const auto v=initial.measure(*mesh,b[0],b[1]);
        const double cap=domain.cellVolumes()[c];double m=v.liquid==0?0:v.gas==0?cap:static_cast<double>(v.liquid);
        if(v.liquid>0 && v.gas>0)m+=sign*(direction?direction:1)*4e-11*cap;
        else if(!plane){m=v.liquid==0?4e-11*cap:cap-4e-11*cap;}
        check(m>=0 && m<=cap,"bad source fixture");sourceM.push_back(m);
    }
    const auto source=stateFor(*mesh,profile,radial,{},[&](int c,Vector2D,Vector2D){return sourceM[c];});
    const auto id=source.identity();const auto* scalar=&source.geometry().snapshot();
    const HermiteGalsPhaseRegionSource geometry(*mesh,source,regionConfig);const HermiteGalsInventoryRegionSource inventory(*mesh,geometry);
    const Characteristic map(*mesh,{mesh->minCoord,mesh->maxCoord,{.125,.125},radial});
    const auto v=Velocity::constant(source,id,{1,0,dt,Position::End},velocity);const Interval base(*mesh,source,v,map,0,dt);
    const bool incomingLiquid=plane?(sign*direction>0):(sign<0);
    const auto entries=inflowEntries(*mesh,base,incomingLiquid);const Inflow history(*mesh,base,entries);
    const Interval full(*mesh,source,v,map,0,dt,history),first(*mesh,source,v,map,0,dt/2,history),last(*mesh,source,v,map,dt/2,dt,history);
    std::vector<Transfer> transfers,oldTransfers;long double maxSplit=0,maxCorrection=0,boundary=0;int internalRow=-1;
    for(const auto face:domain.faces(*mesh)){
        const auto result=integrateHermiteGalsTranslationFace(*mesh,full,face.nativeRow,inventory);
        const auto old=integrateHermiteGalsTranslationFace(*mesh,full,face.nativeRow,geometry);
        check(result.inventoryWeighted && !old.inventoryWeighted && result.transfer.interval().sameInterval(full),"wrong producer/provenance");
        check(&result.transfer.interval().previousAccepted().geometry().snapshot()==scalar,"lost authoritative geometry");
        const auto& values=result.transfer.integratedVolumes();
        check(values.integratedLiquidPositiveVolume>=0 && values.integratedLiquidNegativeVolume>=0
            && values.integratedLiquidPositiveVolume<=values.integratedTotalPositiveVolume
            && values.integratedLiquidNegativeVolume<=values.integratedTotalNegativeVolume,"gross transfer bounds failed");
        check(result.transfer.signedTotalVolume()==old.transfer.signedTotalVolume(),"total flux changed");
        if(result.prescribedInflow)check(result.transfer.signedLiquidVolume()==old.transfer.signedLiquidVolume(),"boundary history changed");
        else{check(result.swept && result.swept->inventoryWeighted,"missing inventory sweep report");
            const auto& r=*result.swept;const double speed=face.horizontal?velocity.y:velocity.x;
            close(r.geometricLiquidVolume,std::abs(old.transfer.signedLiquidVolume()),0,"geometry diagnostic changed original integral");
            close(r.signedInventoryCorrection,result.transfer.signedLiquidVolume()-old.transfer.signedLiquidVolume(),0,"correction diagnostic mismatch");
            if(speed)check(std::abs(r.signedInventoryCorrection)<=r.l1InventoryCorrection+roundoff*std::abs(result.transfer.signedTotalVolume()),"face correction exceeds reported L1");
        }
        const auto left=integrateHermiteGalsTranslationFace(*mesh,first,face.nativeRow,inventory).transfer;
        const auto right=integrateHermiteGalsTranslationFace(*mesh,last,face.nativeRow,inventory).transfer;
        const auto combined=left.appended(*mesh,right);const auto total=std::abs(result.transfer.signedTotalVolume());
        const auto split=std::abs(combined.signedLiquidVolume()-result.transfer.signedLiquidVolume());
        if(total>0){maxSplit=std::max(maxSplit,split/total);check(split<=1e-10*total,"interval subdivision failed");}
        else check(split==0,"zero transfer split changed");
        maxCorrection=std::max(maxCorrection,std::abs(result.transfer.signedLiquidVolume()-old.transfer.signedLiquidVolume()));
        if(face.negativeCell<0 || face.positiveCell<0)boundary+=face.outwardSign*result.transfer.signedLiquidVolume();
        if(face.negativeCell>=0 && face.positiveCell>=0 && face.horizontal==horizontal)internalRow=face.nativeRow;
        transfers.push_back(result.transfer);oldTransfers.push_back(old.transfer);++testedFaces;
    }
    const auto advanced=assembleHermiteGalsIntegratedInventory(*mesh,full,transfers);const auto reflux=assembleHermiteGalsReflux(*mesh,advanced);
    check(&reflux.source().faces()==&advanced.faces() && &reflux.source().phaseVolumes()==&advanced.phaseVolumes(),"reflux duplicated shared owner");
    check(reflux.report().maxLiquidVolumeResidual==0 && reflux.report().maxTotalVolumeResidual==0,"reflux did not cancel exactly");
    check(bool(reflux.report().coarseSides)==amr,"coarse/fine fixture missing");
    long double maxDeparture=0,totalBefore=0,totalAfter=0,totalCapacity=0,maxLocal=0;std::vector<long double> exchange(mesh->cellIndexNum(),0);
    for(const auto& f:advanced.faces()){const auto volume=f.signedLiquidVolume();
        if(f.face().negativeCell>=0)exchange[f.face().negativeCell]-=volume;if(f.face().positiveCell>=0)exchange[f.face().positiveCell]+=volume;}
    for(std::size_t k=0;k<advanced.phaseVolumes().size();++k){const auto c=advanced.phaseVolumes()[k].cell;const auto m=advanced.phaseVolumes()[k].liquid,cap=domain.cellVolumes()[k];
        check(m>=0 && m<=cap && source.primary().phaseVolumes[k].liquid==sourceM[k],"strict M bounds/source immutability failed");
        if(direction==0)check(m==sourceM[k],"zero velocity changed M");
        const auto ref=departureInventory(*mesh,inventory,c,velocity,dt,incomingLiquid,radial);
        maxDeparture=std::max(maxDeparture,std::abs(m-ref)/cap);maxLocal=std::max(maxLocal,std::abs(m-sourceM[k]-exchange[c])/cap);
        totalBefore+=sourceM[k];totalAfter+=m;totalCapacity+=cap;
    }
    check(maxDeparture<=roundoff && maxLocal<=roundoff,"departure/shared local balance exceeds original roundoff budget");
    close(totalAfter-totalBefore+boundary,0,roundoff*totalCapacity,"global boundary ledger");
    check(source.identity()==id && &source.geometry().snapshot()==scalar,"old state changed");
    if(plane){rejects([&]{assembleHermiteGalsIntegratedInventory(*mesh,full,oldTransfers);},"inventory violation exceeds roundoff budget");++oldRejected;}
    if(!radial && !horizontal && !amr && sign==1 && direction==1 && plane && q==.25){
        const auto reject=[&](const auto& fn,const std::string& why){rejects(fn,why);++guardRejections;};
        auto copied=inventory;auto moved=std::move(copied);
        reject([&]{integrateHermiteGalsTranslationFace(*mesh,full,internalRow,copied);},"moved-from");
        auto wrongMesh=makeMesh(8);reject([&]{integrateHermiteGalsTranslationFace(*wrongMesh,full,internalRow,inventory);},"mesh");
        reject([&]{integrateHermiteGalsTranslationFace(*mesh,base,entries.front().nativeFace,inventory);},"missing explicit inflow history");
        const auto foreign=stateFor(*mesh,profile,radial,{},[&](int c,Vector2D,Vector2D){return sourceM[c];});
        const HermiteGalsPhaseRegionSource foreignGeometry(*mesh,foreign,regionConfig);const HermiteGalsInventoryRegionSource foreignInventory(*mesh,foreignGeometry);
        reject([&]{integrateHermiteGalsTranslationFace(*mesh,full,internalRow,foreignInventory);},"identical old region source");
        const auto diagonal=Velocity::constant(source,id,{1,0,dt,Position::End},{1,1});const Interval diagonalInterval(*mesh,source,diagonal,map,0,dt);
        reject([&]{integrateHermiteGalsTranslationFace(*mesh,diagonalInterval,internalRow,inventory);},"single-axis velocity");
        const auto query=[](const Vector2D&){return LevelSetBoundaryTrace::VelocityJet{{1,0},{0,0,0,0}};};
        const Velocity opaque(source,id,{1,0,dt,Position::End},query,query);const Interval opaqueInterval(*mesh,source,opaque,map,0,dt);
        reject([&]{integrateHermiteGalsTranslationFace(*mesh,opaqueInterval,internalRow,inventory);},"owned constant velocity");
        HermiteGalsSweptFaceConfig bad;bad.relativeCoordinateTolerance=0;
        reject([&]{integrateHermiteGalsTranslationFace(*mesh,full,internalRow,inventory,bad);},"invalid swept integral tolerance");
        check(moved.source().identity()==id,"moved source lost identity");
    }
    std::cout<<"INVENTORY_SWEEP,"<<currentCase<<','<<maxDeparture<<','<<maxLocal<<','<<maxSplit<<','<<maxCorrection<<','
        <<advanced.report().closure<<','<<reflux.report().coarseSides<<','<<advanced.report().repairedFaces<<'\n';++cases;
}
}
int main()try
{
    std::cout<<std::setprecision(17);
    for(bool radial:{false,true})for(bool horizontal:{false,true})if(!radial || horizontal)for(bool amr:{false,true})for(int sign:{-1,1}){
        for(int direction:{-1,1}){for(double q:{.25,.75})run(radial,horizontal,amr,sign,direction,true,q);
            run(radial,horizontal,amr,sign,direction,false,0);}
        run(radial,horizontal,amr,sign,0,false,0);
    }
    check(cases==84 && oldRejected==48 && guardRejections==7,"incomplete inventory sweep coverage");
    std::cout<<"I12a11r19 inventory sweep PASS cases="<<cases<<" geometric_rejections="<<oldRejected<<" guards="<<guardRejections<<" faces="<<testedFaces<<'\n';return 0;
}catch(const std::exception& e){std::cerr<<"inventory sweep FAIL case="<<currentCase<<": "<<e.what()<<'\n';return 1;}
