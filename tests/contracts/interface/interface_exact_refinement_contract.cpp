#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Geometry/LevelSetInterfacePanels.h"
#include "Core/Geometry/ActiveQuarticHermitePatch.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Modules/Interface/InterfaceCellMeasure.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>

namespace {
using Snapshot=LevelSetHermiteGeometry::Snapshot;
using Patch=Quantum::Geometry::ActiveQuarticHermitePatch;
using Key=std::array<int,2>;
int cases=0,rejections=0;
const char* stage="setup";
void check(bool ok,const std::string& why){if(!ok)throw std::runtime_error(why);}
template<class F> void rejects(F f,const std::string& expected){try{f();}catch(const std::exception& e){
    check(std::string(e.what()).find(expected)!=std::string::npos,std::string("wrong rejection: ")+e.what());
    ++rejections;std::cout<<"REJECT,"<<expected<<'\n';return;}throw std::runtime_error("missing rejection: "+expected);}
std::unique_ptr<QuadtreeMesh> meshFor(int n,int level=0,Vector2D origin={0,0}){
    DictionaryMesh d;d.Nx=d.Ny=n;d.minCoord=origin;d.maxCoord=origin+Vector2D{1,1};d.maxCellNum=4096;d.maxInitialIterations=0;
    auto mesh=std::make_unique<QuadtreeMesh>(d);
    if(level>1){DictionaryMesh::RefineBox b;b.level=1;b.min=origin+Vector2D{.125,.125};b.max=origin+Vector2D{.875,.875};mesh->applyRefineBoxes({b});}
    if(level){DictionaryMesh::RefineBox b;b.level=level;b.min=origin+Vector2D{.25,.25};b.max=origin+Vector2D{.75,.75};mesh->applyRefineBoxes({b});}
    return mesh;
}
Patch patchFor(const QuadtreeMesh& mesh,bool radial,int sign,bool bubbles){
    const auto jet=[&](Vector2D p){const double x=p.x-mesh.minCoord.x,y=p.y-mesh.minCoord.y,
        a=x-(radial?0:.4375),b=y-.53125,r=radial?.28125:.1875;
        return Patch::NodalDof{sign*(a*a+b*b-r*r+.0078125*x*x*x*x*y*y*y*y),
            sign*(2*a+.03125*x*x*x*y*y*y*y),sign*(2*b+.03125*x*x*x*x*y*y*y),sign*.125*x*x*x*y*y*y};};
    Patch patch;patch.activeCells.assign(mesh.cellIndexNum(),1);patch.nodes.resize(mesh.vertIndexNum());patch.cellCenters.resize(mesh.cellIndexNum());
    for(int v=0;v<mesh.vertIndexNum();++v)patch.nodes[v]=jet(mesh.verts(mesh.vertIndexToLabel[v]).coord());
    std::map<Key,bool> sides;
    for(int c=0;c<mesh.cellIndexNum();++c){const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);
        patch.cellCenters[c]=jet(cell.centerCoord()).value+(bubbles?sign*std::ldexp(std::sin(double(c+1)),-18):0);
        for(int k=0;k<4;++k){Key key{{mesh.verts(cell.cellToVert[k]).index(),mesh.verts(cell.cellToVert[(k+1)%4]).index()}};
            const auto a=mesh.verts(mesh.vertIndexToLabel[key[0]]).coord(),b=mesh.verts(mesh.vertIndexToLabel[key[1]]).coord();
            if(a.y==b.y?a.x>b.x:a.y>b.y)std::swap(key[0],key[1]);sides.emplace(key,true);}}
    for(const auto& side:sides){const auto a=mesh.verts(mesh.vertIndexToLabel[side.first[0]]).coord(),b=mesh.verts(mesh.vertIndexToLabel[side.first[1]]).coord();
        const auto j=jet(a+(b-a)*.5);patch.edgeMidpoints.push_back({side.first,j.value,a.y==b.y?j.dy:j.dx});}
    return patch;
}
double curvature(const Snapshot::ScalarJet2& j){const double n=std::hypot(j.dx,j.dy);check(n>0,"zero reference gradient");
    return (j.dy*j.dy*j.dxx-2*j.dx*j.dy*j.dxy+j.dx*j.dx*j.dyy)/(n*n*n);}
void scenario(bool radial,int sign,bool bubbles,int mode){
    std::cout<<"SCENARIO,"<<radial<<','<<sign<<','<<bubbles<<','<<mode<<'\n';
    const auto origin=mode==4?Vector2D{radial?0:0x1p20,0x1p20}:mode==5?Vector2D{radial?0:.1,.2}:Vector2D{0,0};
    auto oldMesh=meshFor(mode==6?8:4,mode==3?1:0,origin);
    auto newMesh=meshFor(mode==6?16:mode==1 || mode>=4?8:4,mode==2?1:mode==3?2:0,origin);
    stage="source";auto old=std::make_shared<const Snapshot>(Snapshot::fromNativeLeaves(*oldMesh,patchFor(*oldMesh,radial,sign,bubbles),radial));
    if(mode==6){auto seedMesh=meshFor(4);const auto seed=Snapshot::fromNativeLeaves(*seedMesh,patchFor(*seedMesh,radial,sign,bubbles),radial);
        old=std::make_shared<const Snapshot>(Snapshot::refineNativeQuartic(*oldMesh,seed));}
    stage="first refinement";
    const auto refined=std::make_shared<const Snapshot>(Snapshot::refineNativeQuartic(*newMesh,*old));
    check(refined->nativeCellCount()==std::size_t(newMesh->cellIndexNum()) && refined->radialAxisSymmetry()==radial,"wrong refined metadata");
    stage="jets";double jetError=0,traceError=0,curvatureError=0;std::size_t samples=0,quadrature=0;
    for(int c=0;c<newMesh->cellIndexNum();++c){const auto b=refined->nativeCellBounds(c);const auto h=b[1]-b[0];
        for(double u:{.125,.5,.875})for(double v:{.125,.5,.875}){const auto p=b[0]+Vector2D{h.x*u,h.y*v};
            const auto a=old->scalarJetNative(p),r=refined->scalarJetInNativeCell(c,p);
            const double av[]={a.value,a.dx,a.dy,a.dxx,a.dxy,a.dyy},rv[]={r.value,r.dx,r.dy,r.dxx,r.dxy,r.dyy};
            for(int k=0;k<6;++k)jetError=std::max(jetError,std::abs(av[k]-rv[k])/(1+std::abs(av[k])));++samples;}
        for(const auto p:std::vector<Vector2D>{b[0]+Vector2D{0,h.y/2},b[0]+Vector2D{h.x,h.y/2},b[0]+Vector2D{h.x/2,0},b[0]+Vector2D{h.x/2,h.y}}){
            const auto a=old->scalarJetNative(p),r=refined->scalarJetInNativeCell(c,p);
            traceError=std::max({traceError,std::abs(a.value-r.value),std::abs(a.dx-r.dx),std::abs(a.dy-r.dy)});}}
    check(jetError<=1e-10 && traceError<=1e-12,"refinement jet/trace mismatch");
    FaceSegmentGraph oldGraph,newGraph;oldGraph.build(*oldMesh);newGraph.build(*newMesh);
    const auto oldDomain=NativeLeafDomain::build(*oldMesh,oldGraph,radial),newDomain=NativeLeafDomain::build(*newMesh,newGraph,radial);
    stage="volumes";const auto oldVolumes=InterfaceCellMeasure::discoverVolumes(*oldMesh,oldDomain,old),newVolumes=InterfaceCellMeasure::discoverVolumes(*newMesh,newDomain,refined);
    std::vector<long double> totals(oldMesh->cellIndexNum(),0);std::vector<int> parents(newMesh->cellIndexNum(),-1);double volumeError=0;
    const auto volume=[](const auto& m,const auto& d,int c){const int phase=m.cellPhase[c];return phase<0?static_cast<long double>(d.cellVolumes()[c]):phase>0?0.L:m.measure.negativeVolumeByCell[c];};
    for(int c=0;c<newMesh->cellIndexNum();++c){const auto b=refined->nativeCellBounds(c);const auto center=b[0]+(b[1]-b[0])*.5;int parent=-1;
        for(int p:old->nativeCellCandidates(center,center)){const auto q=old->nativeCellBounds(p);
            if(center.x>q[0].x && center.x<q[1].x && center.y>q[0].y && center.y<q[1].y){check(parent<0,"multiple refinement parents");parent=p;}}
        check(parent>=0,"missing refinement parent");parents[c]=parent;totals[parent]+=volume(newVolumes,newDomain,c);}
    for(int c=0;c<oldMesh->cellIndexNum();++c)volumeError=std::max(volumeError,double(std::abs(totals[c]-volume(oldVolumes,oldDomain,c))/oldDomain.cellVolumes()[c]));
    check(volumeError<=1e-10,"refinement phase-volume mismatch");
    stage="curvature";const auto topology=LevelSetInterfacePanels::discoverNativeComponents(*newMesh,refined,radial);
    for(const auto& chain:topology.curves)for(std::size_t g=0;g<chain.graphs().size();++g)for(const auto& q:chain.quadrature(static_cast<int>(g))){
        const int cell=q.frame.nativeCellIndex,parent=parents.at(cell);
        const auto childBox=refined->nativeCellBounds(cell),parentBox=old->nativeCellBounds(parent);
        const auto chart=LevelSetHermiteGeometry::cellChartPoint(q.frame);
        // Keep the retained curve chart. Its rounded global display position
        // loses sub-ULP offsets on translated domains and is not a reference.
        const auto pullback=[](Quantum::Geometry::CellChartCoordinate t,double lo,double hi,double pLo,double pHi){
            return Quantum::Geometry::CellChartCoordinate{std::fma(hi-lo,t.offset,t.fromUpper?pHi-hi:lo-pLo)/(pHi-pLo),t.fromUpper};};
        const Quantum::Geometry::CellChartPoint parentChart{
            pullback(chart.x,childBox[0].x,childBox[1].x,parentBox[0].x,parentBox[1].x),
            pullback(chart.y,childBox[0].y,childBox[1].y,parentBox[0].y,parentBox[1].y)};
        const auto a=old->scalarJetAtNativeCellChart(parent,parentChart);
        curvatureError=std::max(curvatureError,std::abs(curvature(a)-q.frame.curvature));++quadrature;}
    check(quadrature>0 && curvatureError<=1e-9,"refinement curvature consumer mismatch");
    // A second restriction must accept an already exact (non-double) owner.
    stage="repeated refinement";const auto again=Snapshot::refineNativeQuartic(*newMesh,*refined);
    const auto before=refined->scalarJetNative(origin+Vector2D{.37,.61});old.reset();oldMesh.reset();newMesh.reset();
    const auto after=refined->scalarJetNative(origin+Vector2D{.37,.61});check(before.value==after.value && before.dxx==after.dxx,"refined scalar borrowed source lifetime");
    const auto twice=again.scalarJetNative(origin+Vector2D{.37,.61});check(before.value==twice.value && before.dxx==twice.dxx,"repeated restriction changed scalar");
    std::cout<<"REFINEMENT,"<<radial<<','<<sign<<','<<bubbles<<','<<mode<<','<<refined->nativeCellCount()<<','<<samples<<','<<quadrature
        <<','<<jetError<<','<<traceError<<','<<volumeError<<','<<curvatureError<<'\n';++cases;
}
Snapshot sectorSource(const QuadtreeMesh& mesh){
    const auto jet=[](Vector2D p){return Patch::NodalDof{p.x+p.y-1,1,1,0};};
    Patch patch;patch.activeCells.assign(mesh.cellIndexNum(),0);patch.nodes.resize(mesh.vertIndexNum());patch.cellCenters.resize(mesh.cellIndexNum());
    int sectorCell=-1,pin=-1;std::map<Key,bool> sides;
    for(int v=0;v<mesh.vertIndexNum();++v){const auto p=mesh.verts(mesh.vertIndexToLabel[v]).coord();patch.nodes[v]=jet(p);if(p.x==.5 && p.y==.5)pin=v;}
    for(int c=0;c<mesh.cellIndexNum();++c){const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto center=cell.centerCoord();
        if(center.x>.5 && center.y<.5)continue;
        patch.activeCells[c]=1;patch.cellCenters[c]=jet(center).value;
        const auto hi=mesh.verts(cell.cellToVert[2]).coord();if(hi.x==.5 && hi.y==.5)sectorCell=c;
        for(int k=0;k<4;++k){Key key{{mesh.verts(cell.cellToVert[k]).index(),mesh.verts(cell.cellToVert[(k+1)%4]).index()}};
            const auto a=mesh.verts(mesh.vertIndexToLabel[key[0]]).coord(),b=mesh.verts(mesh.vertIndexToLabel[key[1]]).coord();
            if(a.y==b.y?a.x>b.x:a.y>b.y)std::swap(key[0],key[1]);sides.emplace(key,true);}}
    for(const auto& side:sides){const auto a=mesh.verts(mesh.vertIndexToLabel[side.first[0]]).coord(),b=mesh.verts(mesh.vertIndexToLabel[side.first[1]]).coord();
        const auto j=jet(a+(b-a)*.5);patch.edgeMidpoints.push_back({side.first,j.value,a.y==b.y?j.dy:j.dx});}
    check(sectorCell>=0 && pin>=0,"missing corner sector fixture");auto corner=patch.nodes[pin];corner.dy=0;
    patch.cornerSectors.push_back({sectorCell,pin,corner,-1});return Snapshot::fromNativeLeaves(mesh,patch,false);
}
void rejectionCases(){
    auto mesh=meshFor(4),fine=meshFor(8),coarse=meshFor(2),foreign=meshFor(4,0,{1,0});const auto patch=patchFor(*mesh,false,1,false);
    auto source=Snapshot::fromNativeLeaves(*mesh,patch,false);
    rejects([&]{Snapshot::refineNativeQuartic(*coarse,source);},"cannot coarsen or span");
    rejects([&]{Snapshot::refineNativeQuartic(*foreign,source);},"changes domain");
    const auto partialSource=Snapshot::fromNativeQuarticQueryCells(*mesh,patch,{5},false);
    rejects([&]{Snapshot::refineNativeQuartic(*fine,partialSource);},"requires full active source");
    Quantum::Geometry::ActiveHermitePatch cubic;cubic.activeCells=patch.activeCells;cubic.nodes=patch.nodes;
    const auto cubicSource=Snapshot::fromNativeLeaves(*mesh,cubic,false);
    rejects([&]{Snapshot::refineNativeQuartic(*fine,cubicSource);},"requires a native Q4 source");
    const auto sector=sectorSource(*mesh);rejects([&]{Snapshot::refineNativeQuartic(*fine,sector);},"corner sectors");
    auto moved=std::move(source);rejects([&]{Snapshot::refineNativeQuartic(*fine,source);},"moved-from");
    check(std::isfinite(moved.scalarJetNative({.37,.61}).value),"rejection mutated source");
}
}
int main()try{
    std::cout<<std::scientific<<std::setprecision(16);
    for(bool radial:{false,true})for(int sign:{1,-1})for(bool bubbles:{false,true})for(int mode=0;mode<7;++mode)scenario(radial,sign,bubbles,mode);
    rejectionCases();std::cout<<"I12a11r14 exact refinement PASS cases="<<cases<<" rejections="<<rejections<<'\n';return 0;
}catch(const std::exception& e){std::cerr<<"exact refinement FAIL stage="<<stage<<": "<<e.what()<<'\n';return 1;}
