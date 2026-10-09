#include "Modules/Interface/HermiteGalsTemporalPhaseFlux.h"
#include "Modules/Interface/HermiteGalsPhaseRegion.h"
#include "Modules/Interface/HermiteGalsReflux.h"
#include "Modules/Interface/HermiteGalsAmrClosure.h"
#include "Core/Geometry/ActiveQuarticHermitePatch.h"
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Mesh/FaceSegmentGraph.h"
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
namespace {
using namespace Quantum::Interface;
using Patch=Quantum::Geometry::ActiveQuarticHermitePatch;
using Jet=Patch::NodalDof;
using Profile=std::function<Jet(Vector2D)>;
using Key=std::array<int,2>;
using Position=HermiteGalsTimeLayer::Position;
using Velocity=HermiteGalsTemporalVelocitySource;
using Characteristic=HermiteGalsTemporalCharacteristicSource;
using Interval=HermiteGalsPhaseFluxIntervalSource;
using Transfer=HermiteGalsIntegratedFaceTransfer;
using Inflow=HermiteGalsConstantInflowSource;
using Entry=HermiteGalsInflowFacePhase;
using Phase=HermiteGalsInflowPhase;
int cases=0,rejections=0;
void check(bool ok,const std::string& why){if(!ok)throw std::runtime_error(why);}
void near(long double a,long double b,long double scale,const char* why){check(std::abs(a-b)<=64*std::numeric_limits<double>::epsilon()*scale,why);}
template<class F>void rejects(F f,const char* why){try{f();}catch(const std::exception& e){check(std::string(e.what()).find(why)!=std::string::npos,e.what());++rejections;return;}
    throw std::runtime_error(std::string("missing rejection: ")+why);}
std::unique_ptr<QuadtreeMesh> makeMesh(int n,bool amr=false)
{DictionaryMesh d;d.Nx=d.Ny=n;d.maxCellNum=16384;d.maxInitialIterations=0;auto mesh=std::make_unique<QuadtreeMesh>(d);
    if(amr){DictionaryMesh::RefineBox b;b.level=1;b.min={0,.25};b.max={.75,.75};mesh->applyRefineBoxes({b});}return mesh;}
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

HermiteGalsState pure(const QuadtreeMesh& mesh,bool radial,bool liquid){
    const auto domain=domainFor(mesh,radial);
    return stateFor(mesh,[&](Vector2D){return Jet{liquid?-1.:1.,0,0,0};},radial,{},
        [&](int cell,Vector2D,Vector2D){return liquid?domain.cellVolumes().at(domain.cellRowByIndex().at(cell)):0.;});
}
std::vector<Entry> entries(const QuadtreeMesh& mesh,const Interval& base,int pattern){
    std::vector<Entry> result;const auto v=*base.velocitySource().constantValue();const auto domain=domainFor(mesh,base.measure()==HermiteGalsMeasure::AxisymmetricVolume);
    for(const auto face:domain.faces(mesh)){const double speed=face.horizontal?v.y:v.x;
        if((face.negativeCell<0 && speed>0) || (face.positiveCell<0 && speed<0))
            result.push_back({face.nativeRow,pattern==1 || (pattern==2 && result.size()%2==0)?Phase::Liquid:Phase::Gas});}
    return result;
}
void sameVolumes(const Transfer& a,const Transfer& b){
    const auto& x=a.integratedVolumes();const auto& y=b.integratedVolumes();
    check(x.integratedTotalPositiveVolume==y.integratedTotalPositiveVolume && x.integratedTotalNegativeVolume==y.integratedTotalNegativeVolume &&
        x.integratedLiquidPositiveVolume==y.integratedLiquidPositiveVolume && x.integratedLiquidNegativeVolume==y.integratedLiquidNegativeVolume,"changed swept volume");
}
void sameReport(const HermiteGalsSweptFaceReport& a,const HermiteGalsSweptFaceReport& b){
    check(a.coordinateVolumeBand==b.coordinateVolumeBand && a.momentErrorEstimate==b.momentErrorEstimate && a.innerLiquid==b.innerLiquid &&
        a.outerLiquid==b.outerLiquid && a.innerGas==b.innerGas && a.outerGas==b.outerGas && a.rectangleQueries==b.rectangleQueries &&
        a.integratedArcs==b.integratedArcs && a.rootRefinements==b.rootRefinements,"changed swept diagnostics");
}
void matrix(){
    for(bool radial:{false,true})for(bool horizontal:{false,true}){
        if(radial && !horizontal)continue;
        for(bool amr:{false,true}){auto mesh=makeMesh(4,amr);const auto domain=domainFor(*mesh,radial);
            const Characteristic map(*mesh,traceDomain(*mesh,radial));
            for(bool liquid:{false,true}){const auto state=pure(*mesh,radial,liquid);const HermiteGalsPhaseRegionSource regions(*mesh,state);
                for(int direction:{-1,1}){const double dt=.02;const Vector2D value=horizontal?Vector2D{0,direction*.25}:Vector2D{direction*.25,0};
                    const auto velocity=Velocity::constant(state,state.identity(),{1,0,dt,Position::End},value);
                    const Interval base(*mesh,state,velocity,map,0,dt);
                    for(int pattern:{0,1,2}){auto input=entries(*mesh,base,pattern);const Inflow history(*mesh,base,input);
                        check(!base.inflowSource() && history.baseExtent().sameInterval(base),"base interval lost identity");
                        const Inflow copy=history;check(copy.sameSource(history),"inflow copy lost owner");
                        const Interval full(*mesh,state,velocity,map,0,dt,history),first(*mesh,state,velocity,map,0,dt*.375,copy),
                            last(*mesh,state,velocity,map,dt*.375,dt,history);
                        check(full.inflowSource() && full.inflowSource()->sameSource(history) && !full.sameInterval(base),"BC missing from interval identity");
                        std::vector<Transfer> transfers;long double boundary=0,capacity=0;std::size_t inflows=0;
                        for(double c:domain.cellVolumes())capacity+=c;
                        for(const auto face:domain.faces(*mesh)){
                            const auto result=integrateHermiteGalsTranslationFace(*mesh,full,face.nativeRow,regions);
                            const double speed=face.horizontal?value.y:value.x;
                            const bool incoming=(face.negativeCell<0 && speed>0)||(face.positiveCell<0 && speed<0);
                            check(result.prescribedInflow==incoming && result.swept.has_value()!=incoming,"wrong translation source classification");
                            const long double expected=static_cast<long double>(face.area)*speed*dt;
                            near(result.transfer.signedTotalVolume(),expected,std::max(std::abs(expected),1e-30L),"native total measure");
                            if(incoming){++inflows;const auto phase=history.phaseFor(face.nativeRow);
                                check(result.transfer.signedLiquidVolume()==(phase==Phase::Liquid?result.transfer.signedTotalVolume():0),"prescribed pure phase");
                                const auto joined=integrateHermiteGalsTranslationFace(*mesh,first,face.nativeRow,regions).transfer.appended(*mesh,
                                    integrateHermiteGalsTranslationFace(*mesh,last,face.nativeRow,regions).transfer);
                                near(joined.signedLiquidVolume(),result.transfer.signedLiquidVolume(),std::abs(expected),"BC interval composition");
                                check(joined.interval().sameInterval(full) && joined.interval().inflowSource()->sameSource(history),"composition lost BC owner");
                            }else{
                                const auto swept=integrateHermiteGalsSweptFace(*mesh,full,face.nativeRow,regions);
                                sameVolumes(result.transfer,swept.transfer);sameReport(*result.swept,swept.report);
                                check(result.transfer.signedLiquidVolume()==(liquid?result.transfer.signedTotalVolume():0),"old pure phase identity");
                            }
                            if(face.negativeCell<0 || face.positiveCell<0)boundary+=face.outwardSign*result.transfer.signedLiquidVolume();
                            transfers.push_back(result.transfer);
                        }
                        check(inflows==input.size(),"inflow inventory size");
                        const auto assembled=assembleHermiteGalsIntegratedInventory(*mesh,full,transfers);
                        const auto reflux=assembleHermiteGalsReflux(*mesh,assembled);
                        check(&reflux.source().faces()==&assembled.faces() && reflux.source().interval().inflowSource()->sameSource(history),"reflux lost boundary owner/shared faces");
                        long double before=0,after=0;
                        for(std::size_t k=0;k<assembled.phaseVolumes().size();++k){
                            const auto c=domain.cellVolumes()[k],v=assembled.phaseVolumes()[k].liquid;before+=state.primary().phaseVolumes[k].liquid;after+=v;
                            check(v>=0 && v<=c,"strict assembled capacity");check(state.primary().phaseVolumes[k].liquid==(liquid?c:0),"source M mutated");
                            if(pattern==int(liquid))near(v,liquid?c:0,c,"uniform phase preservation");
                        }
                        near(after-before+boundary,0,capacity,"complete boundary/volume ledger");
                        near(assembled.report().boundaryOutward,boundary,capacity,"boundary report");
                        ++cases;
                    }
                }
            }
        }
    }
}
void guards(){
    auto mesh=makeMesh(4);const auto state=pure(*mesh,false,false);const HermiteGalsPhaseRegionSource regions(*mesh,state);
    const Characteristic map(*mesh,traceDomain(*mesh,false));const double dt=.02;
    const auto v=Velocity::constant(state,state.identity(),{1,0,dt,Position::End},{1,0});const Interval base(*mesh,state,v,map,0,dt);
    auto input=entries(*mesh,base,1);const Inflow history(*mesh,base,input);const int row=input.front().nativeFace;
    input.front().phase=Phase::Gas;check(history.phaseFor(row)==Phase::Liquid,"mutable input escaped inflow owner");
    input=history.entries();const Interval full(*mesh,state,v,map,0,dt,history);
    rejects([&]{integrateHermiteGalsTranslationFace(*mesh,base,row,regions);},"missing explicit inflow history");
    auto bad=input;bad.pop_back();rejects([&]{Inflow invalid(*mesh,base,bad);},"inflow face entry");
    bad=input;bad.insert(bad.begin(),bad.front());rejects([&]{Inflow invalid(*mesh,base,bad);},"inflow face entry");
    bad=input;std::reverse(bad.begin(),bad.end());rejects([&]{Inflow invalid(*mesh,base,bad);},"inflow face entry");
    bad=input;bad.push_back({999,Phase::Gas});rejects([&]{Inflow invalid(*mesh,base,bad);},"foreign/non-inflow");
    bad=input;bad.front().phase=static_cast<Phase>(7);rejects([&]{Inflow invalid(*mesh,base,bad);},"invalid inflow phase");
    rejects([&]{Inflow invalid(*mesh,full,input);},"unbound base interval");
    rejects([&]{history.phaseFor(-1);},"missing prescribed inflow face");
    const Inflow different(*mesh,base,input);check(!different.sameSource(history),"identical data forged owner");
    const Interval first(*mesh,state,v,map,0,dt/2,history),last(*mesh,state,v,map,dt/2,dt,different);
    const auto a=integrateHermiteGalsTranslationFace(*mesh,first,row,regions).transfer,b=integrateHermiteGalsTranslationFace(*mesh,last,row,regions).transfer;
    rejects([&]{a.appended(*mesh,b);},"different inflow history");
    const Interval missing(*mesh,state,v,map,dt/2,dt);
    const Transfer noBc(*mesh,missing,row,b.orientation(),b.integratedVolumes());
    rejects([&]{a.appended(*mesh,noBc);},"different inflow history");
    const Interval shortBase(*mesh,state,v,map,0,dt/2);const Inflow shorter(*mesh,shortBase,input);
    rejects([&]{Interval invalid(*mesh,state,v,map,0,dt,shorter);},"does not cover interval");
    const auto v2=Velocity::constant(state,state.identity(),{1,0,dt,Position::End},{1,0});
    rejects([&]{Interval invalid(*mesh,state,v2,map,0,dt,history);},"different velocity source");
    const Characteristic map2(*mesh,traceDomain(*mesh,false));
    rejects([&]{Interval invalid(*mesh,state,v,map2,0,dt,history);},"different characteristic source");
    const auto state2=pure(*mesh,false,false);const auto vv=Velocity::constant(state2,state2.identity(),{1,0,dt,Position::End},{1,0});
    rejects([&]{Interval invalid(*mesh,state2,vv,map,0,dt,history);},"different old source");
    auto copy=history;auto moved=std::move(copy);rejects([&]{copy.entries();},"moved-from inflow");
    check(moved.sameSource(history),"moved source lost owner");
    const auto query=[](Vector2D){return LevelSetBoundaryTrace::VelocityJet{{1,0},{0,0,0,0}};};
    const Velocity opaque(state,state.identity(),{1,0,dt,Position::End},query,query);const Interval oi(*mesh,state,opaque,map,0,dt);
    rejects([&]{Inflow invalid(*mesh,oi,input);},"owned constant velocity");
    const auto angled=Velocity::constant(state,state.identity(),{1,0,dt,Position::End},{1,1});const Interval ai(*mesh,state,angled,map,0,dt);
    rejects([&]{Inflow invalid(*mesh,ai,input);},"single-axis velocity");
    auto unknown=traceDomain(*mesh,false);unknown.activeSegment=[](Vector2D,Vector2D){return true;};const Characteristic um(*mesh,unknown);
    const Interval ui(*mesh,state,v,um,0,dt);rejects([&]{Inflow invalid(*mesh,ui,input);},"opaque coverage");
    const auto rz=pure(*mesh,true,false);const Characteristic rm(*mesh,traceDomain(*mesh,true));
    const auto rv=Velocity::constant(rz,rz.identity(),{1,0,dt,Position::End},{1,0});const Interval ri(*mesh,rz,rv,rm,0,dt);
    rejects([&]{Inflow invalid(*mesh,ri,input);},"axial velocity");
    auto differentMesh=makeMesh(8);rejects([&]{history.validateInterval(*differentMesh,full);},"mesh");
    std::vector<Transfer> all;const auto domain=domainFor(*mesh,false);
    for(const auto face:domain.faces(*mesh))all.push_back(integrateHermiteGalsTranslationFace(*mesh,full,face.nativeRow,regions).transfer);
    const auto original=all.front();all.front()=Transfer(*mesh,base,original.face().nativeRow,original.orientation(),original.integratedVolumes());
    rejects([&]{assembleHermiteGalsIntegratedInventory(*mesh,full,all);},"different interval source");
    const auto zero=Velocity::constant(state,state.identity(),{1,0,dt,Position::End},{0,0});const Interval zi(*mesh,state,zero,map,0,dt);
    const Inflow empty(*mesh,zi,{});const Interval boundZero(*mesh,state,zero,map,0,dt,empty);
    const auto noMotion=integrateHermiteGalsTranslationFace(*mesh,boundZero,row,regions);
    check(!noMotion.prescribedInflow && noMotion.transfer.signedLiquidVolume()==0 && noMotion.transfer.signedTotalVolume()==0,"empty zero-speed history");++cases;
    const auto retained=[&]{const Inflow temporary(*mesh,base,input);return Interval(*mesh,state,v,map,0,dt,temporary);}();
    check(integrateHermiteGalsTranslationFace(*mesh,retained,row,regions).transfer.signedLiquidVolume()>0,"history lifetime");++cases;
    mesh->verts(mesh->vertIndexToLabel[0]).setCoord({-.125,0});
    rejects([&]{integrateHermiteGalsTranslationFace(*mesh,full,row,regions);},"mesh");
}
}
int main(int argc,char** argv){try{if(argc==2 && std::string(argv[1])=="--guards")guards();else{matrix();guards();}
    std::cout<<"I12a8 constant inflow PASS cases="<<cases<<" rejections="<<rejections<<'\n';return 0;
}catch(const std::exception& e){std::cerr<<"I12a8 constant inflow FAIL: "<<e.what()<<'\n';return 1;}}
