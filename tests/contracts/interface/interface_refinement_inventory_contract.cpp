#include "Modules/Interface/HermiteGalsRefinementInventory.h"
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
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace
{
using namespace LAMP::Interface;
using Patch=LAMP::Geometry::ActiveQuarticHermitePatch;
using Snapshot=LevelSetHermiteGeometry::Snapshot;
using Jet=Patch::NodalDof;
using Key=std::array<int,2>;
using Profile=std::function<Jet(Vector2D)>;
constexpr double epsilon=std::numeric_limits<double>::epsilon();
int cases=0,rejections=0;
static_assert(!std::is_default_constructible_v<HermiteGalsRefinementInventory>);
static_assert(std::is_same_v<decltype(std::declval<HermiteGalsRefinementInventory&>().phaseVolumes()),const std::vector<HermiteGalsCellVolume>&>);
void check(bool ok,const std::string& why){if(!ok)throw std::runtime_error(why);}
template<class F> void rejects(F f,const std::string& expected){try{f();}catch(const std::exception& e){
    check(std::string(e.what()).find(expected)!=std::string::npos,std::string("wrong rejection: ")+e.what());
    if(expected=="no feasible child allocation"){
        const std::string message=e.what(),prefix="stage=parent allocation parent=";
        const auto start=message.find(prefix);
        check(start!=std::string::npos,"missing allocation context: "+message);
        const int parent=std::stoi(message.substr(start+prefix.size()));
        check(parent>=0 && parent<16 && message.find(" child=-1]")!=std::string::npos,
            "wrong infeasible-parent context: "+message);
    }
    ++rejections;std::cout<<"REJECT,"<<expected<<'\n';return;}throw std::runtime_error("missing rejection: "+expected);}
std::unique_ptr<QuadtreeMesh> meshFor(int n,int level=0,Vector2D origin={0,0}){
    DictionaryMesh d;d.Nx=d.Ny=n;d.minCoord=origin;d.maxCoord=origin+Vector2D{1,1};d.maxCellNum=4096;d.maxInitialIterations=0;
    auto mesh=std::make_unique<QuadtreeMesh>(d);
    if(level>1){DictionaryMesh::RefineBox b;b.level=1;b.min={.125,.125};b.max={.875,.875};mesh->applyRefineBoxes({b});}
    if(level){DictionaryMesh::RefineBox b;b.level=level;b.min={.25,.25};b.max={.75,.75};mesh->applyRefineBoxes({b});}return mesh;
}
Vector2D point(const QuadtreeMesh& mesh,int vertex){return mesh.verts(mesh.vertIndexToLabel[vertex]).coord();}
Profile profile(bool radial,int sign,bool pure=false){return [=](Vector2D p){
    if(pure)return Jet{double(sign),0,0,0};
    const double a=p.x-(radial?0:.4375),b=p.y-.53125,r=radial?.28125:.1875;
    return Jet{sign*(a*a+b*b-r*r+.0078125*std::pow(p.x,4)*std::pow(p.y,4)),
        sign*(2*a+.03125*std::pow(p.x,3)*std::pow(p.y,4)),
        sign*(2*b+.03125*std::pow(p.x,4)*std::pow(p.y,3)),sign*.125*std::pow(p.x,3)*std::pow(p.y,3)};};}
// Explicit old-state operator fixture. Initial M is measured once, then the
// requested independent inventory perturbation is retained with the same Q4.
HermiteGalsState sourceFor(const QuadtreeMesh& mesh,const Profile& f,bool radial,double epsShift=0,
    bool partial=false,const HermiteGalsPhaseVolumeConfig& config={}){
    Patch patch;patch.activeCells.assign(mesh.cellIndexNum(),0);patch.nodes.resize(mesh.vertIndexNum());patch.cellCenters.resize(mesh.cellIndexNum());
    HermiteGalsPrimaryState p;p.measure=radial?HermiteGalsMeasure::AxisymmetricVolume:HermiteGalsMeasure::CartesianArea;
    std::vector<bool> used(mesh.vertIndexNum(),false);std::map<Key,bool> sides;
    for(int v=0;v<mesh.vertIndexNum();++v)patch.nodes[v]=f(point(mesh,v));
    for(int c=0;c<mesh.cellIndexNum();++c){const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);
        if(partial && cell.centerCoord().x>.5)continue;
        patch.activeCells[c]=1;p.geometryCells.push_back(c);p.phaseVolumes.push_back({c,0});patch.cellCenters[c]=f(cell.centerCoord()).value;
        for(int k=0;k<4;++k){Key key{{mesh.verts(cell.cellToVert[k]).index(),mesh.verts(cell.cellToVert[(k+1)%4]).index()}};
            used[key[0]]=used[key[1]]=true;const auto a=point(mesh,key[0]),b=point(mesh,key[1]);
            if(a.y==b.y?a.x>b.x:a.y>b.y)std::swap(key[0],key[1]);sides.emplace(key,true);}}
    HermiteGalsClosure closure;
    for(int v=0;v<mesh.vertIndexNum();++v)if(used[v] && mesh.verts(mesh.vertIndexToLabel[v]).type()==TVertType::non){
        const auto j=patch.nodes[v];p.vertices.push_back({v,j.value,j.dx,j.dy});closure.vertexMixedDerivatives.push_back(j.dxy);}
    for(const auto& side:sides){const auto a=point(mesh,side.first[0]),b=point(mesh,side.first[1]);const bool horizontal=a.y==b.y;
        const auto j=f(a+(b-a)*.5);const double dn=horizontal?j.dy:j.dx;patch.edgeMidpoints.push_back({side.first,j.value,dn});
        bool contained=false;for(const auto& other:sides)if(side.first!=other.first){const auto lo=point(mesh,other.first[0]),hi=point(mesh,other.first[1]);
            if(horizontal && lo.y==hi.y && a.y==lo.y && a.x>=lo.x && b.x<=hi.x)contained=true;
            if(!horizontal && lo.x==hi.x && a.x==lo.x && a.y>=lo.y && b.y<=hi.y)contained=true;}
        if(!contained)p.edges.push_back({side.first,j.value,dn});}
    for(int c:p.geometryCells)closure.cellBubbleCoefficients.push_back(patch.cellCenters[c]);
    const auto make=[&]{HermiteGalsCandidate candidate(mesh,{7,.2,.21,HermiteGalsTimeLayer::Position::End},p);
        candidate.bindGeometry(mesh,candidate.identity(),closure,Snapshot::fromNativeLeaves(mesh,patch,radial));return std::move(candidate).finish();};
    const auto seed=make();
    // This rejection fixture is strictly constant gas on partial support;
    // M=0 is exact. It does not invoke I07's full/boundary-contact discovery.
    if(partial){for(int c:p.geometryCells)check(patch.cellCenters[c]==1,"partial fixture must be constant gas");return seed;}
    const auto measured=integrateHermiteGalsPhaseVolumes(mesh,seed,config);p.phaseVolumes=measured.phaseVolumes();
    for(std::size_t k=0;k<p.phaseVolumes.size();++k)if(measured.cells()[k].certifiedPhase==0)
        p.phaseVolumes[k].liquid+=epsShift*epsilon*measured.cells()[k].capacity;
    return make();
}
void scenario(bool radial,int sign,int perturb,int mode,bool pure=false){
    std::cout<<"SCENARIO,"<<radial<<','<<sign<<','<<perturb<<','<<mode<<','<<pure<<'\n';
    auto oldMesh=meshFor(4,mode==3?1:0),newMesh=meshFor(mode==1?8:mode==4?16:4,mode==2?1:mode==3?2:0);
    auto source=std::make_shared<const HermiteGalsState>(sourceFor(*oldMesh,profile(radial,sign,pure),radial,4*perturb));
    const auto before=source->primary();const auto id=source->identity();const auto time=source->timeLayer();
    const auto result=refineHermiteGalsInventory(*oldMesh,*source,*newMesh);result.validateMeshes(*oldMesh,*newMesh);
    check(result.sourceIdentity()==id && result.timeLayer().step==time.step && result.timeLayer().time()==time.time(),"lost source/time provenance");
    check(result.phaseVolumes().size()==std::size_t(newMesh->cellIndexNum()),"incomplete child inventory");
    auto owner=std::make_shared<const HermiteGalsRefinementInventory>(result);
    std::shared_ptr<const Snapshot> scalar(owner,&owner->geometry());
    FaceSegmentGraph graph;graph.build(*newMesh);const auto domain=NativeLeafDomain::build(*newMesh,graph,radial);
    const auto fresh=InterfaceCellMeasure::discoverVolumes(*newMesh,domain,scalar);
    std::vector<long double> totals(oldMesh->cellIndexNum(),0);double parentError=0;bool retainsDifference=false;
    for(std::size_t k=0;k<result.cells().size();++k){const auto& cell=result.cells()[k];const auto m=result.phaseVolumes()[k].liquid;
        check(cell.cell==int(k) && m>=0 && m<=cell.capacity,"child inventory owner/bounds");
        check(cell.rawLiquid==fresh.measure.negativeVolumeByCell[k] && cell.certifiedPhase==fresh.cellPhase[k],"measurement changed scalar source");
        totals.at(result.parentCells()[k])+=m;
        const auto childBox=scalar->nativeCellBounds(int(k)),parentBox=source->geometry().snapshot().nativeCellBounds(result.parentCells()[k]);
        check(childBox[0].x>=parentBox[0].x && childBox[0].y>=parentBox[0].y && childBox[1].x<=parentBox[1].x && childBox[1].y<=parentBox[1].y,"false parent provenance");}
    for(const auto& balance:result.parentBalances()){
        const auto original=before.phaseVolumes.at(balance.parent).liquid;
        check(balance.originalInventory==original,"parent target was reset");
        const double error=double(std::abs(totals[balance.parent]-original)/balance.capacity);parentError=std::max(parentError,error);
        check(error<=4*epsilon,"parent inventory loss beyond summation roundoff");
        if(balance.children==1){const auto found=std::find(result.parentCells().begin(),result.parentCells().end(),balance.parent);
            check(result.phaseVolumes()[found-result.parentCells().begin()].liquid==original,"unchanged leaf inventory changed");}
        retainsDifference=retainsDifference || std::abs(balance.childInventory-balance.childGeometry)>epsilon*balance.capacity;
    }
    if(perturb)check(retainsDifference,"perturbed independent inventory was reset to geometry");
    check(source->identity()==id && source->primary().phaseVolumes.size()==before.phaseVolumes.size(),"source authority mutated");
    for(std::size_t k=0;k<before.phaseVolumes.size();++k)check(source->primary().phaseVolumes[k].liquid==before.phaseVolumes[k].liquid,"source inventory mutated");
    if(mode==1 && !perturb){const auto repeated=refineHermiteGalsInventory(*oldMesh,*source,*newMesh);
        for(std::size_t k=0;k<result.phaseVolumes().size();++k)check(result.phaseVolumes()[k].liquid==repeated.phaseVolumes()[k].liquid,"repeated transfer differs");}
    auto movedFrom=result;auto moved=std::move(movedFrom);check(&moved.geometry()==&result.geometry(),"copy/move changed scalar owner");
    if(!radial && sign==1 && mode==0 && !perturb && !pure)rejects([&]{movedFrom.report();},"moved-from");
    const auto value=result.geometry().scalarJetNative({.37,.61}).value;oldMesh.reset();newMesh.reset();source.reset();
    check(result.geometry().scalarJetNative({.37,.61}).value==value,"result borrowed external mesh");
    const auto& r=result.report();std::cout<<"REFINE_INVENTORY,"<<radial<<','<<sign<<','<<perturb<<','<<mode<<','<<pure<<','
        <<r.refinedParents<<','<<r.unchangedParents<<','<<r.childCells<<','<<parentError<<','<<r.maxRelativeSourceGeometryDefect<<','
        <<r.maxRelativeChildGeometryDefect<<','<<r.maxRelativeParentReferenceDifference<<'\n';++cases;
}
void rejectionCases(){
    auto oldMesh=meshFor(4),fine=meshFor(8),coarse=meshFor(2),foreign=meshFor(4,0,{1,0});
    const auto source=sourceFor(*oldMesh,profile(false,1),false);const auto id=source.identity();
    const auto prepared=refineHermiteGalsInventory(*oldMesh,source,*fine);
    rejects([&]{prepared.validateMeshes(*foreign,*fine);},"mesh");
    rejects([&]{prepared.validateMeshes(*oldMesh,*foreign);},"mesh");
    rejects([&]{refineHermiteGalsInventory(*foreign,source,*fine);},"mesh");
    rejects([&]{refineHermiteGalsInventory(*oldMesh,source,*foreign);},"changes domain");
    rejects([&]{refineHermiteGalsInventory(*oldMesh,source,*coarse);},"cannot coarsen or span");
    const auto over=sourceFor(*oldMesh,profile(false,1),false,1e17);
    rejects([&]{refineHermiteGalsInventory(*oldMesh,over,*fine);},"outside strict capacity bounds");
    const auto inconsistent=sourceFor(*oldMesh,profile(false,1),false,1e8);
    rejects([&]{refineHermiteGalsInventory(*oldMesh,inconsistent,*fine);},"original geometry roundoff budget");
    const auto partial=sourceFor(*oldMesh,profile(false,1,true),false,0,true);
    rejects([&]{refineHermiteGalsInventory(*oldMesh,partial,*fine);},"full geometry and inventory coverage");
    HermiteGalsPhaseVolumeConfig config;constexpr double intercept=.5+.25/256;
    config.boundaryContacts={{intercept,0},{0,intercept}};
    const auto plane=[=](Vector2D p){return Jet{intercept-p.x-p.y,-1,-1,0};};
    const auto weak=sourceFor(*oldMesh,plane,false,32,false,config);auto deep=meshFor(16);
    rejects([&]{refineHermiteGalsInventory(*oldMesh,weak,*deep,config);},"no feasible child allocation");
    check(source.identity()==id,"rejected transfer changed source");
}
}
int main()try{
    std::cout<<std::scientific<<std::setprecision(16);
    for(bool radial:{false,true})for(int sign:{1,-1}){
        for(int mode=0;mode<4;++mode)for(int perturb:{0,1,-1})scenario(radial,sign,perturb,mode);
        scenario(radial,sign,0,4);for(int mode:{1,2})scenario(radial,sign,0,mode,true);}
    rejectionCases();std::cout<<"I12a11r15 refinement inventory PASS cases="<<cases<<" rejections="<<rejections<<'\n';return 0;
}catch(const std::exception& e){std::cerr<<"refinement inventory FAIL: "<<e.what()<<'\n';return 1;}
