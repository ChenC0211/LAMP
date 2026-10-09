#include "Modules/Interface/HermiteGalsTemporalPhaseFlux.h"
#include "Modules/Interface/HermiteGalsInventoryRegion.h"
#include "Modules/Interface/HermiteGalsAmrClosure.h"
#include "Core/Geometry/ActiveQuarticHermitePatch.h"
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <array>
#include <cfenv>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
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
using Velocity=HermiteGalsTemporalVelocitySource;
using Characteristic=HermiteGalsTemporalCharacteristicSource;
using Interval=HermiteGalsPhaseFluxIntervalSource;
using Coverage=HermiteGalsCharacteristicPhaseCoverage;
using Phase=HermiteGalsCertifiedPhase;
using Position=HermiteGalsTimeLayer::Position;
using Patch=Quantum::Geometry::ActiveQuarticHermitePatch;
using Jet=Patch::NodalDof;
using Profile=std::function<Jet(Vector2D)>;
using Inventory=std::function<double(int,Vector2D,Vector2D,double)>;
using Key=std::array<int,2>;
static_assert(!std::is_default_constructible_v<Coverage>);
static_assert(std::is_same_v<decltype(std::declval<Coverage&>().cells()),const std::vector<HermiteGalsCharacteristicPhaseCell>&>);
int queries=0,samples=0,traces=0,guards=0,thinEvents=0,weightedCases=0;
std::size_t checkedLeaves=0,maximumLeaves=0;
void check(bool ok,const std::string& why){if(!ok)throw std::runtime_error(why);}
template<class F> void rejects(F action,const char* expected)
{try{action();}catch(const std::exception& e){check(std::string(e.what()).find(expected)!=std::string::npos,e.what());++guards;return;}
    throw std::runtime_error(std::string("missing rejection: ")+expected);}
std::unique_ptr<QuadtreeMesh> makeMesh(bool amr)
{DictionaryMesh d;d.Nx=d.Ny=4;d.maxCellNum=1024;d.maxInitialIterations=0;auto m=std::make_unique<QuadtreeMesh>(d);
    if(amr){DictionaryMesh::RefineBox b;b.level=1;b.min={.25,.25};b.max={.75,.75};m->applyRefineBoxes({b});}return m;}
Vector2D point(const QuadtreeMesh& m,int v){return m.verts(m.vertIndexToLabel[v]).coord();}
bool fineSide(const QuadtreeMesh& m,Key key)
{const auto a=point(m,key[0]),b=point(m,key[1]);const bool h=a.y==b.y;
    for(int v:key){const auto owner=m.verts(m.vertIndexToLabel[v]).supCell;if(owner<0)continue;
        const auto& c=m.cells(owner);const auto lo=m.verts(c.cellToVert[0]).coord(),hi=m.verts(c.cellToVert[2]).coord();
        if(h && (a.y==lo.y || a.y==hi.y) && a.x>=lo.x && b.x<=hi.x && b.x-a.x<hi.x-lo.x)return true;
        if(!h && (a.x==lo.x || a.x==hi.x) && a.y>=lo.y && b.y<=hi.y && b.y-a.y<hi.y-lo.y)return true;}return false;}
// Analytic primary/center samples provide a structural Q4 source. Only the
// explicitly weighted fixtures bind analytic M; no interface is advanced here.
HermiteGalsState stateFor(const QuadtreeMesh& m,const Profile& f,bool radial,const Inventory& inventory={},bool complete=true)
{
    FaceSegmentGraph graph;graph.build(m);const auto domain=NativeLeafDomain::build(m,graph,radial);
    HermiteGalsPrimaryState p;Patch patch;p.measure=radial?HermiteGalsMeasure::AxisymmetricVolume:HermiteGalsMeasure::CartesianArea;
    patch.activeCells.assign(m.cellIndexNum(),1);patch.nodes.resize(m.vertIndexNum());patch.cellCenters.resize(m.cellIndexNum());
    for(int v=0;v<m.vertIndexNum();++v)if(m.verts(m.vertIndexToLabel[v]).type()==TVertType::non){const auto j=f(point(m,v));
        patch.nodes[v]=j;p.vertices.push_back({v,j.value,j.dx,j.dy});}
    std::set<Key> edges;std::set<int> activeVertices;
    for(int c=0;c<m.cellIndexNum();++c){const auto& cell=m.cells(m.cellIndexToLabel[c]);const auto a=m.verts(cell.cellToVert[0]).coord(),b=m.verts(cell.cellToVert[2]).coord();
        if(complete || c+1<m.cellIndexNum()){p.geometryCells.push_back(c);p.phaseVolumes.push_back({c,inventory?inventory(c,a,b,domain.cellVolumes()[c]):0.});}
        else patch.activeCells[c]=0;
        patch.cellCenters[c]=f(a+(b-a)*.5).value;
        if(!patch.activeCells[c])continue;
        for(int k=0;k<4;++k){Key key{{m.verts(cell.cellToVert[k]).index(),m.verts(cell.cellToVert[(k+1)%4]).index()}};
            activeVertices.insert(key[0]);
            const auto x=point(m,key[0]),y=point(m,key[1]);if(x.y==y.y?x.x>y.x:x.y>y.y)std::swap(key[0],key[1]);edges.insert(key);}}
    p.vertices.erase(std::remove_if(p.vertices.begin(),p.vertices.end(),[&](const auto& v){return activeVertices.find(v.vertex)==activeVertices.end();}),p.vertices.end());
    for(const auto key:edges)if(!fineSide(m,key)){const auto a=point(m,key[0]),b=point(m,key[1]);const auto j=f(a+(b-a)*.5);const double dn=a.y==b.y?j.dy:j.dx;
        p.edges.push_back({key,j.value,dn});patch.edgeMidpoints.push_back({key,j.value,dn});}
    HermiteGalsCandidate candidate(m,{2,.1,.2,Position::End},p);const auto mixed=closeHermiteGalsMixedAmr(m,candidate);
    HermiteGalsClosure closure;closure.vertexMixedDerivatives=mixed.vertexMixedDerivatives();
    for(int c:p.geometryCells)closure.cellBubbleCoefficients.push_back(patch.cellCenters[c]);
    for(std::size_t k=0;k<p.vertices.size();++k)patch.nodes[p.vertices[k].vertex].dxy=closure.vertexMixedDerivatives[k];
    candidate.bindGeometry(m,candidate.identity(),closure,LevelSetHermiteGeometry::Snapshot::fromNativeLeaves(m,patch,radial));return std::move(candidate).finish();
}
Profile shape(int kind,bool radial=false,int sign=1)
{
    return [=](Vector2D p){
        if(kind==0)return Jet{sign*(p.y-.4),0.,double(sign),0.};
        if(kind==5)return Jet{double(sign),0,0,0};
        const double x=p.x-(radial?0.:kind==4?.25:kind==6?.53125:.5),y=p.y-((kind==2 || kind==3 || kind==6)?.53125:.5);
        const double r=kind==2?.0009765625:kind==3?0:kind==4?.125:kind==6?.0078125:.25;
        if(kind==2 || kind==3)return Jet{sign*(y*y-r*r),0,sign*2*y,0};
        return Jet{sign*(x*x+y*y-r*r),sign*2*x,sign*2*y,0};
    };
}
Quantum::Geometry::PolynomialVelocityData flowData(int kind)
{
    Quantum::Geometry::PolynomialVelocityData d;
    if(kind==1 || kind==2 || kind==3){d.origin={.5,.5};
        if(kind<3){d.previous[1].x=-.8;d.previous[4].y=.8;}
        else{d.previous[9].x=-.8;d.previous[3].x=-.8;d.previous[12].y=.8;d.previous[6].y=.8;}}
    if(kind==4)d.previous[1].x=.25;
    if(kind==5)d.previous[0].y=.25;
    if(kind==6){d.previous[4].x=.2;d.previous[0].y=.2;d.previous[1].y=-.4;}
    d.arrival=d.previous;
    if(kind==2 || kind==3 || kind==6)for(auto& c:d.arrival)c=c*(-.5);
    return d;
}
Interval intervalFor(const QuadtreeMesh& mesh,const HermiteGalsState& old,bool amr,int flow,double begin=.2,double end=.3)
{
    HermiteGalsCandidate next(mesh,{3,.2,.3,Position::End},old.primary());
    const auto velocity=Velocity::polynomial(old,next.identity(),next.timeLayer(),flowData(flow));
    const double h=amr?.125:.25;const Characteristic characteristic(mesh,{mesh.minCoord,mesh.maxCoord,{h,h},old.primary().measure==HermiteGalsMeasure::AxisymmetricVolume});
    return Interval(mesh,old,velocity,characteristic,begin,end);
}
bool inside(Vector2D p,Vector2D lo,Vector2D hi)
{return p.x>=lo.x && p.x<=hi.x && p.y>=lo.y && p.y<=hi.y;}
void verify(const QuadtreeMesh& mesh,const Coverage& result)
{
    ++queries;result.validateMesh(mesh);const auto& characteristic=result.characteristic();const auto& box=characteristic.bounds();
    const auto& interval=characteristic.interval();const auto& old=interval.previousAccepted();const auto& scalar=old.geometry().snapshot();
    const auto& cells=result.cells();std::vector<int> expected;
    // Independent exhaustive native scan checks the indexed closed incidence,
    // including leaves touching Q only along an edge or at one point.
    for(int c=0;c<mesh.cellIndexNum();++c){const auto b=scalar.nativeCellBounds(c);
        if(std::max(b[0].x,box.departureLower.x)<=std::min(b[1].x,box.departureUpper.x)
            && std::max(b[0].y,box.departureLower.y)<=std::min(b[1].y,box.departureUpper.y))expected.push_back(c);}
    check(cells.size()==expected.size() && !cells.empty(),"closed leaf incidence count");
    for(std::size_t i=0;i<cells.size();++i){check(cells[i].cell==expected[i],"closed leaf incidence order");
        check(0<=cells[i].minimumLiquidContent && cells[i].minimumLiquidContent<=cells[i].maximumLiquidContent && cells[i].maximumLiquidContent<=1,"leaf content bounds");}
    checkedLeaves+=cells.size();maximumLeaves=std::max(maximumLeaves,cells.size());
    for(int ti=0;ti<5;++ti){const double t=ti==4?interval.timeEnd():interval.timeBegin()+(interval.timeEnd()-interval.timeBegin())*ti/4;
        for(int xi=0;xi<3;++xi)for(int yi=0;yi<3;++yi){const Vector2D p{box.arrivalLower.x+(box.arrivalUpper.x-box.arrivalLower.x)*xi/2,
                box.arrivalLower.y+(box.arrivalUpper.y-box.arrivalLower.y)*yi/2};Vector2D departure=p;
            if(t>old.timeLayer().time()){HermiteGalsCandidate actual(mesh,{3,.2,t,Position::End},old.primary());
                const auto& source=interval.velocitySource();const auto trace=traceHermiteGalsWithJacobian(interval.characteristicSource().domain(),p,
                    source.previous(),source.polynomialArrivalView(actual.identity(),actual.timeLayer()));
                check(!trace.boundaryHit() && trace.acceptedSubsteps()==1,"phase coverage lost certified RK branch");departure=trace.departurePoint();++traces;}
            const bool negative=scalar.negative(departure);++samples;
            if(result.phase()!=Phase::Unresolved)check(negative==(result.phase()==Phase::Liquid),"global phase excluded actual scalar sign");
            bool found=false;
            for(const auto& c:cells)if(inside(departure,c.lower,c.upper)){found=true;
                if(c.phase!=Phase::Unresolved)check(negative==(c.phase==Phase::Liquid),"leaf phase excluded actual scalar sign");
                long double content=negative?1:0;
                if(result.inventorySource()){const auto& weights=result.inventorySource()->cells().at(c.cell);content=negative?weights.liquidWeight:weights.gasWeight;}
                check(content>=c.minimumLiquidContent && content<=c.maximumLiquidContent,"actual content escaped leaf range");
                check(content>=result.minimumLiquidContent() && content<=result.maximumLiquidContent(),"actual content escaped global range");}
            check(found,"actual departure has no covered native leaf");}}
}
Coverage cover(const QuadtreeMesh& mesh,const Interval& interval,Vector2D lo,Vector2D hi)
{return coverHermiteGalsCharacteristicPhase(mesh,encloseHermiteGalsCharacteristic(mesh,interval,lo,hi));}
void ordinary(bool radial,bool amr)
{
    auto mesh=makeMesh(amr);const std::vector<int> flows=radial?std::vector<int>{0,5,6}:std::vector<int>{0,1,2,3,4,5};
    for(int kind:{0,1}){const auto old=stateFor(*mesh,shape(kind,radial),radial);const auto inventory=old.primary().phaseVolumes;
        const auto id=old.identity();const auto* scalar=&old.geometry().snapshot();
        for(int f:flows){const auto interval=intervalFor(*mesh,old,amr,f);
            const Vector2D liquidLo{radial?.10:.49,kind==0?.23:.49},liquidHi{radial?.12:.51,kind==0?.25:.51};
            // The original wide-time far-corner rotation box remains an
            // explicit conservative CFL refusal on the fine mesh. Use a
            // resolved pure-gas region for the positive classification test.
            if(amr && !radial && kind==1 && f==1)
                rejects([&]{cover(*mesh,interval,{.78,.72},{.80,.74});},"half-cell speed guard");
            const Vector2D gasLo{kind==0?.49:radial?.49:.78,kind==0?.60:radial?.72:.49},
                gasHi{kind==0?.51:radial?.51:.80,kind==0?.62:radial?.74:.51};
            for(bool liquid:{true,false}){const auto result=cover(*mesh,interval,liquid?liquidLo:gasLo,liquid?liquidHi:gasHi);
                check(result.phase()==(liquid?Phase::Liquid:Phase::Gas),"resolved ordinary phase missing");verify(*mesh,result);
                check(!result.inventorySource() && result.characteristic().interval().sameInterval(interval),"geometry provenance changed");}
            if(radial && kind==1){const auto axis=cover(*mesh,interval,{0,.49},{0,.51});check(axis.phase()==Phase::Liquid,"axis pure phase unresolved");verify(*mesh,axis);}}
        check(old.identity()==id && &old.geometry().snapshot()==scalar,"phase query changed old geometry");
        for(std::size_t i=0;i<inventory.size();++i)check(old.primary().phaseVolumes[i].liquid==inventory[i].liquid,"phase query changed M");}
    // Exact point/line incidence at uniform and coarse/fine seams.
    const auto seamOld=stateFor(*mesh,shape(radial?0:4,radial),radial);const auto interval=intervalFor(*mesh,seamOld,amr,0);
    for(const auto bounds:std::array<std::array<Vector2D,2>,3>{{{{{.24,.49},{.26,.51}}},{{{.25,.49},{.25,.51}}},{{{.25,.5},{.25,.5}}}}}){
        const auto result=cover(*mesh,interval,bounds[0],bounds[1]);check(result.phase()==(radial?Phase::Gas:Phase::Liquid) && result.cells().size()>=3,"seam incidence missing");verify(*mesh,result);}
    std::cout<<"PHASE_MATRIX,"<<radial<<','<<amr<<','<<flows.size()<<",PASS\n";
}
void hiddenEvents(bool amr)
{
    auto mesh=makeMesh(amr);
    for(int kind:{2,3}){const auto old=stateFor(*mesh,shape(kind),false);const auto interval=intervalFor(*mesh,old,amr,5);
        const Vector2D arrival{.5,.55};const auto result=cover(*mesh,interval,arrival,arrival);
        check(result.phase()==Phase::Unresolved,"hidden temporal zero incorrectly certified absent");verify(*mesh,result);
        for(double t:{.2,.25,.3}){const Vector2D p{arrival.x,arrival.y-.25*(t-.2)};check(old.geometry().snapshot().value(p)>0,"sparse time nodes must all be gas");}
        const auto refined=cover(*mesh,intervalFor(*mesh,old,amr,5,.274,.276),arrival,arrival);
        if(kind==2){check(refined.phase()==Phase::Liquid,"resolved hidden thin liquid interval missing");++thinEvents;}
        else {check(refined.phase()==Phase::Unresolved,"same-sign tangency incorrectly excluded");
            check(old.geometry().snapshot().value({arrival.x,.53125})==0,"exact tangency fixture lost its zero");}verify(*mesh,refined);
        for(const auto times:std::array<std::array<double,2>,2>{{{{.2,.22}},{{.29,.3}}}}){
            const auto pure=cover(*mesh,intervalFor(*mesh,old,amr,5,times[0],times[1]),arrival,arrival);check(pure.phase()==Phase::Gas,"refined pure time range missing");verify(*mesh,pure);}}
    const auto bubble=stateFor(*mesh,shape(6),false);const auto interval=intervalFor(*mesh,bubble,amr,0);
    const auto result=cover(*mesh,interval,{.5,.5},{.625,.625});check(result.phase()==Phase::Unresolved,"hidden spatial liquid falsely excluded");verify(*mesh,result);
    for(const Vector2D p:std::array<Vector2D,5>{{{.5,.5},{.625,.5},{.625,.625},{.5,.625},{.5625,.5625}}})
        check(!bubble.geometry().snapshot().negative(p),"spatial corner/center fixture not all gas");
    check(bubble.geometry().snapshot().negative({.53125,.53125}),"hidden spatial bubble lost");
    const auto liquid=cover(*mesh,interval,{.530,.530},{.532,.532});check(liquid.phase()==Phase::Liquid,"refined hidden bubble not certified");verify(*mesh,liquid);++thinEvents;
    std::cout<<"HIDDEN_EVENTS,"<<amr<<",temporal_thin_layer,tangency,spatial_bubble,PASS\n";
}
void weighted(bool radial,bool amr)
{
    auto mesh=makeMesh(amr);
    for(int sign:{-1,1}){const auto old=stateFor(*mesh,shape(5,radial,sign),radial,[=](int c,Vector2D,Vector2D,double cap){
            const double defect=2e-11*(1+c%3);return cap*(sign<0?1-defect:defect);});
        const auto inventoryBefore=old.primary().phaseVolumes;const HermiteGalsPhaseRegionSource regions(*mesh,old);
        const HermiteGalsInventoryRegionSource inventory(*mesh,regions);const auto interval=intervalFor(*mesh,old,amr,0);
        const auto enclosed=encloseHermiteGalsCharacteristic(*mesh,interval,{0,0},{1,1});
        const auto result=coverHermiteGalsCharacteristicPhase(*mesh,enclosed,inventory);verify(*mesh,result);++weightedCases;
        check(result.phase()==(sign<0?Phase::Liquid:Phase::Gas) && result.cells().size()==static_cast<std::size_t>(mesh->cellIndexNum()),"weighted pure coverage");
        check(result.minimumLiquidContent()>0 && result.maximumLiquidContent()<1,"geometry purity erased independent inventory tail");
        check(result.minimumLiquidContent()<result.maximumLiquidContent(),"cross-leaf weight variation erased");
        check(&result.inventorySource()->source().geometry().snapshot()==&old.geometry().snapshot(),"weighted geometry authority changed");
        for(std::size_t k=0;k<inventoryBefore.size();++k)check(old.primary().phaseVolumes[k].liquid==inventoryBefore[k].liquid,"weighted coverage modified M");}
    const auto old=stateFor(*mesh,shape(0),radial,[](int,Vector2D a,Vector2D b,double cap){
        const double f=std::clamp((.4-a.y)/(b.y-a.y),0.,1.);return cap*(f==1?1-4e-11:f+4e-11);});
    HermiteGalsPhaseRegionConfig config;config.boundaryContacts={{1,.4}};if(!radial)config.boundaryContacts.push_back({0,.4});
    const HermiteGalsPhaseRegionSource regions(*mesh,old,config);const HermiteGalsInventoryRegionSource inventory(*mesh,regions);
    const auto interval=intervalFor(*mesh,old,amr,0);
    int region=0;
    for(const auto bounds:std::array<std::array<Vector2D,2>,3>{{{{{.24,.27},{.51,.30}}},{{{.24,.43},{.51,.56}}},{{{.24,.39},{.51,.41}}}}}){
        const auto result=coverHermiteGalsCharacteristicPhase(*mesh,encloseHermiteGalsCharacteristic(*mesh,interval,bounds[0],bounds[1]),inventory);
        check(result.phase()==(region==0?Phase::Liquid:region==1?Phase::Gas:Phase::Unresolved),"weighted phase classification changed");
        if(region==1)check(result.minimumLiquidContent()>0 && result.minimumLiquidContent()<result.maximumLiquidContent(),"gas/cut inventory weight hull lost");
        verify(*mesh,result);++region;}
    ++weightedCases;std::cout<<"WEIGHTED_COVERAGE,"<<radial<<','<<amr<<",pure_tails,cut_leaf_weights,PASS\n";
}
void guardsAndLifetime()
{
    auto mesh=makeMesh(false);const auto old=stateFor(*mesh,shape(5),false,[](int,Vector2D,Vector2D,double cap){return 4e-11*cap;});
    const auto originalM=old.primary().phaseVolumes;const auto* scalar=&old.geometry().snapshot();
    const HermiteGalsPhaseRegionSource regions(*mesh,old);auto inventory=HermiteGalsInventoryRegionSource(*mesh,regions);
    auto enclosure=encloseHermiteGalsCharacteristic(*mesh,intervalFor(*mesh,old,false,0),{.4,.4},{.6,.6});
    auto result=coverHermiteGalsCharacteristicPhase(*mesh,enclosure,inventory);const auto same=result;auto moved=std::move(result);
    rejects([&]{result.cells();},"moved-from characteristic phase coverage");check(&same.cells()==&moved.cells(),"copy lost coverage owner");
    const auto retained=std::move(enclosure);rejects([&]{coverHermiteGalsCharacteristicPhase(*mesh,enclosure);},"moved-from characteristic enclosure");
    const auto retainedInventory=std::move(inventory);rejects([&]{coverHermiteGalsCharacteristicPhase(*mesh,retained,inventory);},"moved-from inventory region");
    const auto other=stateFor(*mesh,shape(5),false,[](int,Vector2D,Vector2D,double cap){return 4e-11*cap;});
    const HermiteGalsPhaseRegionSource otherRegions(*mesh,other);const HermiteGalsInventoryRegionSource foreign(*mesh,otherRegions);
    rejects([&]{coverHermiteGalsCharacteristicPhase(*mesh,retained,foreign);},"foreign inventory/geometry");
    auto different=makeMesh(true);rejects([&]{moved.validateMesh(*different);},"mesh");
    rejects([&]{coverHermiteGalsCharacteristicPhase(*different,retained);},"mesh");
    const auto partial=stateFor(*mesh,shape(5),false,{},false);const auto partialInterval=intervalFor(*mesh,partial,false,0);
    rejects([&]{encloseHermiteGalsCharacteristic(*mesh,partialInterval,{.4,.4},{.6,.6});},"full old scalar support");
    {struct Restore{int previous=std::fegetround();~Restore(){std::fesetround(previous);}} restore;
        check(std::fesetround(FE_UPWARD)==0,"cannot change test rounding");rejects([&]{coverHermiteGalsCharacteristicPhase(*mesh,retained);},"round-to-nearest");}
#if defined(_M_X64) || (defined(__x86_64__) && defined(__SSE2__))
    {struct Restore{unsigned previous=_mm_getcsr();~Restore(){_mm_setcsr(previous);}} restore;
        _mm_setcsr(restore.previous|(1u<<15));rejects([&]{coverHermiteGalsCharacteristicPhase(*mesh,retained);},"gradual underflow");}
#endif
    std::optional<Coverage> lifetime;
    {const auto localOld=stateFor(*mesh,shape(5),false,[](int,Vector2D,Vector2D,double cap){return 4e-11*cap;});
        const HermiteGalsPhaseRegionSource localRegions(*mesh,localOld);const HermiteGalsInventoryRegionSource localInventory(*mesh,localRegions);
        const auto local=encloseHermiteGalsCharacteristic(*mesh,intervalFor(*mesh,localOld,false,0),{.4,.4},{.6,.6});
        lifetime=coverHermiteGalsCharacteristicPhase(*mesh,local,localInventory);}
    verify(*mesh,*lifetime);check(lifetime->inventorySource() && lifetime->minimumLiquidContent()>0,"retained weighted owner lost");
    check(&old.geometry().snapshot()==scalar,"failure changed old geometry");
    for(std::size_t k=0;k<originalM.size();++k)check(old.primary().phaseVolumes[k].liquid==originalM[k].liquid,"failure changed old M");
    DictionaryMesh::RefineBox box;box.level=1;box.min={.25,.25};box.max={.75,.75};mesh->applyRefineBoxes({box});
    rejects([&]{lifetime->validateMesh(*mesh);},"mesh");
    std::cout<<"PHASE_SOURCE_GUARDS_AND_LIFETIME,PASS checks="<<guards<<'\n';
}
}
int main() try
{
    std::cout<<std::scientific<<std::setprecision(16);
    for(bool radial:{false,true})for(bool amr:{false,true})ordinary(radial,amr);
    for(bool amr:{false,true})hiddenEvents(amr);
    for(bool radial:{false,true})for(bool amr:{false,true})weighted(radial,amr);
    guardsAndLifetime();check(queries>=125 && samples==45*queries && traces>4000 && thinEvents==4 && weightedCases==12 && guards==11,"incomplete phase matrix");
    std::cout<<"I12a14 characteristic phase PASS queries="<<queries<<" samples="<<samples<<" traces="<<traces<<" leaf_intersections="<<checkedLeaves
        <<" max_leaves="<<maximumLeaves<<" hidden_liquid_cases="<<thinEvents<<" weighted_cases="<<weightedCases<<" guards="<<guards<<'\n';
}catch(const std::exception& e){std::cerr<<"I12a14 characteristic phase FAIL: "<<e.what()<<'\n';return 1;}
