#include "Modules/Interface/GeometricVofPublication.h"
#include "Modules/Interface/GeometricVofTransport.h"
#include "Modules/Interface/GeometricVofMeshTransfer.h"
#include "Core/Geometry/Detail/FiniteVolumeClosure.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>

namespace {
using namespace LAMP::Interface;
using Publication=GeometricVofPublication;
using Clock=GeometricVofTimeLayer;
using Slab=GeometricVofAdvectorSlab;
using Transport=GeometricVofTransport;
using Transfer=GeometricVofMeshTransfer;
namespace Exact=LAMP::Geometry::Detail::FiniteVolume;
std::size_t scenes=0,guards=0,continued=0,microphases=0,overcapacity=0,nonzeroCertificates=0;
void check(bool value,const char* message){if(!value)throw std::runtime_error(std::string("VOF_CHECKPOINT: ")+message);}
template<class F> void reject(F&& f){bool rejected=false;try{f();}catch(const std::exception&){rejected=true;}check(rejected,"negative control accepted");++guards;}
NativeLeafDomain domain(const QuadtreeMesh& mesh,bool rz){FaceSegmentGraph graph;graph.build(mesh);return NativeLeafDomain::build(mesh,graph,rz);}
DictionaryMesh input(){DictionaryMesh cfg;cfg.Nx=cfg.Ny=4;cfg.minCoord={0,0};cfg.maxCoord={1,1};cfg.maxInitialIterations=0;cfg.levelSetCount=0;cfg.maxCellNum=2048;cfg.maxLevel_local={2};return cfg;}
std::unique_ptr<QuadtreeMesh> meshFor(bool left=false){auto mesh=std::make_unique<QuadtreeMesh>(input());DictionaryMesh::RefineBox box;
    box.min={left?0.:.5,0};box.max={left?.5:1.,1};box.level=1;mesh->applyRefineBoxes({box});return mesh;}
bool bits(long double a,long double b){if constexpr(sizeof(long double)==sizeof(double))return std::memcmp(&a,&b,sizeof a)==0;
    return a==b&&std::signbit(a)==std::signbit(b);}
void same(const Publication& a,const Publication& b){
    check(a.identity()!=b.identity(),"restore retained publication identity");check(a.checkpoint()==b.checkpoint(),"complete sealed state changed across restore");
    const auto x=a.timeLayer(),y=b.timeLayer();check(x.step==y.step&&bits(x.begin,y.begin)&&bits(x.end,y.end)&&x.position==y.position,"full clock changed");
    check(a.plicCells().size()==b.plicCells().size(),"native PLIC size changed");
    for(std::size_t ci=0;ci<a.plicCells().size();++ci){const auto& p=a.plicCells()[ci];const auto& q=b.plicCells()[ci];
        check(bits(p.offset,q.offset)&&bits(p.normal.x,q.normal.x)&&bits(p.normal.y,q.normal.y)&&p.resolved==q.resolved,"actual PLIC plane/resolution changed");
        for(unsigned phase=0;phase<2;++phase)check(bits(p.inventory[phase],q.inventory[phase])&&bits(p.measured[phase],q.measured[phase]),"raw independent phase or measure changed");
        check(a.capacityDefects()[ci].compare(b.capacityDefects()[ci])==0,"exact capacity certificate changed");
        microphases+=p.inventory[0]>0&&p.inventory[1]>0&&std::min(p.inventory[0],p.inventory[1])<1e-16L*a.cells()[ci].capacity;
        overcapacity+=p.inventory[0]>a.cells()[ci].capacity||p.inventory[1]>a.cells()[ci].capacity;
        nonzeroCertificates+=a.capacityDefects()[ci].compare(Exact::Dyadic{})!=0;
    }
    ++scenes;
}
std::vector<double> circulation(const QuadtreeMesh& mesh,const Publication& source,double dt){
    const auto psi=[](Vector2D p){const long double x=p.x,y=p.y;return x*x*(1-x)*(1-x)*y*y*(1-y)*(1-y);};
    std::vector<double> normal(source.domain().nativeFaceCount());
    for(const auto f:source.domain().faces(mesh)){const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
        const auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();const auto volume=(f.horizontal?-1:1)*(psi(b)-psi(a));
        if(f.area==0||f.outwardSign!=0){check(volume==0,"closed streamfunction leaks");normal[f.nativeRow]=0;}
        else normal[f.nativeRow]=static_cast<double>(volume/f.area);}
    const auto rates=Transport::localPlicCourantRates(mesh,source,normal);const auto peak=std::max({rates.maximumLocalColumn,rates.maximumIncoming,rates.maximumOutflow});
    check(peak>0,"circulation fixture is empty");for(auto& value:normal)value=static_cast<double>(value*(.125L/(dt*peak)));return normal;
}
void continuation(const QuadtreeMesh& mesh,const Publication& source,bool moving){
    auto restored=Publication::restoreCheckpoint(mesh,domain(mesh,source.domain().axisymmetric()),source.checkpoint());same(source,restored);
    const auto clock=source.timeLayer();const Clock end{clock.step+1,clock.time(),clock.time()+1e-4,Clock::Position::End};
    const auto normal=moving?circulation(mesh,source,end.end-end.begin):std::vector<double>(source.domain().nativeFaceCount());
    const auto a=Slab::prepare(mesh,source,end,normal,normal),b=Slab::prepare(mesh,restored,end,normal,normal);
    const auto actual=Transport::prepareLocalPlicSweep(mesh,source,a,2,1),resumed=Transport::prepareLocalPlicSweep(mesh,restored,b,2,1);
    same(actual.endpoint(),resumed.endpoint());++continued;
    const auto restoredEnd=Publication::restoreCheckpoint(mesh,domain(mesh,source.domain().axisymmetric()),actual.endpoint().checkpoint());same(actual.endpoint(),restoredEnd);
    if(moving){const Clock next{end.step+1,end.end,end.end+1e-4,Clock::Position::End};const auto velocity=circulation(mesh,actual.endpoint(),next.end-next.begin);
        const auto nextA=Slab::prepare(mesh,actual.endpoint(),next,velocity,velocity),nextB=Slab::prepare(mesh,restoredEnd,next,velocity,velocity);
        const auto uninterrupted=Transport::prepareLocalPlicSweep(mesh,actual.endpoint(),nextA,2,1),continuedEnd=Transport::prepareLocalPlicSweep(mesh,restoredEnd,nextB,2,1);
        same(uninterrupted.endpoint(),continuedEnd.endpoint());++continued;}
}
void ordinary(bool rz){auto mesh=meshFor();const auto d=domain(*mesh,rz);const Clock clock{0,.2,.2,Clock::Position::Initial};
    std::vector<std::array<long double,2>> stocks(mesh->cellIndexNum());
    for(int ci=0;ci<mesh->cellIndexNum();++ci){const auto center=mesh->cells(mesh->cellIndexToLabel[ci]).centerCoord();const auto capacity=static_cast<long double>(d.cellVolumes()[ci]);
        const long double fraction=center.x+center.y<1?.25L:.75L;stocks[ci]={capacity*fraction,capacity*(1-fraction)};}
    stocks.front()={-0.L,d.cellVolumes().front()};stocks[1]={d.cellVolumes()[1],7*std::numeric_limits<long double>::denorm_min()};
    const auto source=Publication::prepareLocalPlicInventory(*mesh,d,clock,stocks,{});const auto image=source.checkpoint();
    const auto restored=Publication::restoreCheckpoint(*mesh,domain(*mesh,rz),image);same(source,restored);
    check(std::signbit(restored.plicCells().front().inventory[0]),"negative-zero raw phase lost its bits");
    continuation(*mesh,source,true);
    auto foreign=meshFor(true);reject([&]{Publication::restoreCheckpoint(*foreign,domain(*foreign,rz),image);});
    reject([&]{Publication::restoreCheckpoint(*mesh,domain(*mesh,!rz),image);});
    auto corrupt=image;corrupt[corrupt.size()/2]^=1;reject([&]{Publication::restoreCheckpoint(*mesh,d,corrupt);});
    reject([&]{Publication::restoreCheckpoint(*mesh,d,image.substr(0,image.size()-1));});
    reject([&]{Publication::restoreCheckpoint(*mesh,d,image+"x");});
    FaceSegmentGraph graph;graph.build(*mesh,[](const auto&,const Cell& cell){return cell.centerCoord().x<.5;});
    const auto selected=NativeLeafDomain::build(*mesh,graph,rz);reject([&]{Publication::restoreCheckpoint(*mesh,selected,image);});
    std::vector<std::array<long double,2>> selectedStocks(stocks.size());for(int ci:selected.cellIndexByRow())selectedStocks[ci]={selected.cellVolumes()[selected.cellRowByIndex()[ci]],0};
    const auto excluded=Publication::prepareLocalPlicInventory(*mesh,selected,clock,selectedStocks,{});reject([&]{excluded.checkpoint();});
    const Clock midpoint{1,.2,.3,Clock::Position::Midpoint};const auto middle=Publication::prepareLocalPlicInventory(*mesh,d,midpoint,stocks,{});reject([&]{middle.checkpoint();});
    check(source.checkpoint()==image,"failed restore changed accepted source");
}
void finiteCertificate(){if constexpr(std::numeric_limits<long double>::digits!=53)return;
    DictionaryMesh cfg;cfg.Nx=cfg.Ny=1;cfg.minCoord={0,9};cfg.maxCoord={1,10};cfg.maxInitialIterations=0;cfg.levelSetCount=0;cfg.maxCellNum=2048;cfg.maxLevel_local={4};
    QuadtreeMesh coarse(cfg);DictionaryMesh::RefineBox box;box.min=cfg.minCoord;box.max=cfg.maxCoord;box.level=1;
    {auto once=coarse.clone();once->applyRefineBoxes({box});check(once->cellIndexNum()==4,"single-cell refinement did not publish four actual leaves");
        for(int ci=0;ci<once->cellIndexNum();++ci){const auto& child=once->cells(once->cellIndexToLabel[ci]);
            check(child.relatedPreIndex.size()==1&&child.relatedPreIndex[0]==0&&child.relatedPreCoord.size()==1,
                "single-cell refinement invented a missing opposite donor");}}
    auto fine=coarse.clone();box.level=4;fine->applyRefineBoxes({box});check(fine->cellIndexNum()==256,"actual L4 checkpoint witness changed tree");
    const auto d=domain(coarse,true),fd=domain(*fine,true);const Clock clock{7,.125,.25,Clock::Position::End};
    for(unsigned phase:{0u,1u}){std::vector<std::array<long double,2>> initial(1);initial[0][phase]=d.cellVolumes()[0];
        const auto root=Publication::prepareLocalPlicInventory(coarse,d,clock,initial,{});const auto remap=Transfer::prepare(coarse,root,*fine,fd,clock);const auto& source=remap.candidate();
        std::vector<std::array<long double,2>> raw(source.plicCells().size());bool excess=false;
        for(std::size_t ci=0;ci<raw.size();++ci){raw[ci]=source.plicCells()[ci].inventory;excess=excess||raw[ci][phase]>source.cells()[ci].capacity;}
        check(excess,"actual native CDF did not produce certified overcapacity witness");reject([&]{Publication::prepareLocalPlicInventory(*fine,fd,clock,raw,{});});
        const auto restored=Publication::restoreCheckpoint(*fine,domain(*fine,true),source.checkpoint());same(source,restored);continuation(*fine,source,false);
    }
}
}
int main()try{ordinary(false);ordinary(true);finiteCertificate();
    check(scenes>=6&&continued>=2&&microphases>0&&guards>=16,"required checkpoint witnesses incomplete");
    if constexpr(std::numeric_limits<long double>::digits==53)check(overcapacity>0&&nonzeroCertificates>0,"actual typed capacity certificate coverage missing");
    std::cout<<"VOF_CHECKPOINT_PASS,scenes,"<<scenes<<",continuations,"<<continued<<",guards,"<<guards<<",microphases,"<<microphases
        <<",certified_overcapacity_cells,"<<overcapacity<<",nonzero_certificates,"<<nonzeroCertificates
        <<",XY_RZ,AMR,raw_planes_stocks_and_support,Initial_End,fresh_identity,no_reconstruction\n";return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
