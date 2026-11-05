#include "Modules/Interface/HermiteGalsPublication.h"
#include "Modules/Interface/HermiteGalsFirstJetStep.h"
#include "Modules/Interface/HermiteGalsPhaseFlux.h"
#include "Core/Geometry/LevelSetMappedSegment.h"

#include "Modules/Interface/HermiteGalsPhaseVolume.h"
#include "Modules/Interface/InterfaceCellMeasure.h"
#include "Core/Geometry/ActiveQuarticHermitePatch.h"
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
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
using namespace LAMP::Interface;
using Patch=LAMP::Geometry::ActiveQuarticHermitePatch;
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
HermiteGalsState stateFor(const QuadtreeMesh& mesh,const Profile& f,bool radial,HermiteGalsTimeLayer time,const Inventory& inventory,
    HermiteGalsTransport transport=HermiteGalsTransport::ConservativeVolume)
{
    HermiteGalsPrimaryState p;p.measure=radial?HermiteGalsMeasure::AxisymmetricVolume:HermiteGalsMeasure::CartesianArea;p.transport=transport;
    Patch patch;patch.activeCells.assign(mesh.cellIndexNum(),1);patch.nodes.resize(mesh.vertIndexNum());patch.cellCenters.resize(mesh.cellIndexNum());
    for(int v=0;v<mesh.vertIndexNum();++v)if(mesh.verts(mesh.vertIndexToLabel[v]).type()==TVertType::non){const auto j=f(point(mesh,v));patch.nodes[v]=j;p.vertices.push_back({v,j.value,j.dx,j.dy});}
    std::set<Key> keys;
    for(int c=0;c<mesh.cellIndexNum();++c){const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto a=mesh.verts(cell.cellToVert[0]).coord(),b=mesh.verts(cell.cellToVert[2]).coord();
        p.geometryCells.push_back(c);if(transport==HermiteGalsTransport::ConservativeVolume)p.phaseVolumes.push_back({c,inventory(c,a,b)});patch.cellCenters[c]=f(a+(b-a)*.5).value;
        for(int k=0;k<4;++k){Key key{{mesh.verts(cell.cellToVert[k]).index(),mesh.verts(cell.cellToVert[(k+1)%4]).index()}};
            const auto x=point(mesh,key[0]),y=point(mesh,key[1]);if(x.y==y.y?x.x>y.x:x.y>y.y)std::swap(key[0],key[1]);keys.insert(key);}}
    for(const auto& key:keys)if(!fineSide(mesh,key)){const auto a=point(mesh,key[0]),b=point(mesh,key[1]);const auto j=f(a+(b-a)*.5);const auto dn=a.y==b.y?j.dy:j.dx;
        p.edges.push_back({key,j.value,dn});patch.edgeMidpoints.push_back({key,j.value,dn});}
    HermiteGalsCandidate candidate(mesh,time,p);const auto mixed=closeHermiteGalsMixedAmr(mesh,candidate);
    HermiteGalsClosure closure;closure.vertexMixedDerivatives=mixed.vertexMixedDerivatives();closure.cellBubbleCoefficients=patch.cellCenters;
    for(std::size_t k=0;k<p.vertices.size();++k)patch.nodes[p.vertices[k].vertex].dxy=closure.vertexMixedDerivatives[k];
    candidate.bindGeometry(mesh,candidate.identity(),closure,LevelSetHermiteGeometry::Snapshot::fromNativeLeaves(mesh,patch,radial));return std::move(candidate).finish();
}
using Snapshot=LevelSetHermiteGeometry::Snapshot;
using Interval=Snapshot::QuarticCenterInterval;
static_assert(!std::is_default_constructible_v<HermiteGalsBubbleResult>);
static_assert(std::is_same_v<decltype(std::declval<HermiteGalsBubbleResult&>().closure()),const HermiteGalsClosure&>);
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
static_assert(!std::is_default_constructible_v<HermiteGalsPublication>);
static_assert(std::is_same_v<decltype(std::declval<HermiteGalsPublication&>().state()),const HermiteGalsState&>);
HermiteGalsPrimaryState measuredPrimary(const QuadtreeMesh& mesh,const Profile& f,bool radial)
{
    const auto seed=stateFor(mesh,f,radial,{},[](int,Vector2D,Vector2D){return 0.;});
    const auto volume=integrateHermiteGalsPhaseVolumes(mesh,seed);
    auto p=seed.primary();p.phaseVolumes=volume.phaseVolumes();return p;
}
HermiteGalsPublication publish(const QuadtreeMesh& mesh,const Profile& f,const HermiteGalsPrimaryState& p,HermiteGalsTimeLayer time={})
{
    HermiteGalsCandidate candidate(mesh,time,p);const auto id=candidate.identity();
    const auto mixed=closeHermiteGalsMixedAmr(mesh,candidate);
    const auto bubble=closeHermiteGalsBubble(mesh,candidate,mixed,{id,centers(mesh,f)});
    auto result=publishHermiteGalsGeometry(mesh,std::move(candidate),bubble);
    check(result.state().identity()==id,"publication minted another identity");
    rejects([&]{candidate.identity();},"consumed");
    rejects([&]{publishHermiteGalsGeometry(mesh,std::move(candidate),bubble);},"consumed");
    return result;
}
void macroScalarBinding()
{
    auto mesh=makeMesh(2);const Profile f=[](Vector2D p){return Jet{(p.x-.5)*(p.x-.5)+(p.y-.5)*(p.y-.5)-.13*.13,2*(p.x-.5),2*(p.y-.5),0};};
    const auto seed=stateFor(*mesh,f,false,{},[](int,Vector2D,Vector2D){return 0.;},HermiteGalsTransport::FirstJet);
    const HermiteGalsTimeLayer clock{1,0,.01,Position::Midpoint};
    HermiteGalsCandidate candidate(*mesh,clock,seed.primary());
    const auto mixed=closeHermiteGalsMixedAmr(*mesh,candidate);
    HermiteGalsPublicationConfig cfg;cfg.scalarRepresentation=HermiteGalsScalarRepresentation::C2MacroQuartic;
    HermiteGalsCenterPrediction predicted{candidate.identity(),centers(*mesh,f)};
    predicted.values[0]+=.001;const auto expected=predicted.values;
    HermiteGalsCandidate foreign(*mesh,clock,seed.primary());
    rejects([&]{prepareHermiteGalsFirstJetScalar(*mesh,foreign,mixed,predicted,cfg);},"revision");
    auto bad=predicted;bad.values.pop_back();
    rejects([&]{prepareHermiteGalsFirstJetScalar(*mesh,candidate,mixed,bad,cfg);},"layout");
    auto scalar=prepareHermiteGalsFirstJetScalar(*mesh,candidate,mixed,predicted,cfg);
    check(!candidate.hasGeometry(),"scalar preparation consumed/bound candidate prematurely");
    predicted.values[0]+=1;
    for(int c=0;c<mesh->cellIndexNum();++c)
        close(scalar.scalarJetAtNativeCellChart(c,{{.5,false},{.5,false}}).value,expected[c],2e-13,"macro actual immutable center");
    candidate.bindGeometry(*mesh,candidate.identity(),{mixed.vertexMixedDerivatives(),expected},std::move(scalar));
    const auto state=std::move(candidate).finish();
    check(state.geometry().snapshot().usesC2Macro()&&state.timeLayer().step==clock.step
        &&state.timeLayer().position==clock.position,"macro binding lost representation/full clock");
    const Profile radialProfile=[](Vector2D p){return Jet{p.x*p.x+(p.y-.5)*(p.y-.5)-.13*.13,2*p.x,2*(p.y-.5),0};};
    const auto radialSeed=stateFor(*mesh,radialProfile,true,{},[](int,Vector2D,Vector2D){return 0.;},HermiteGalsTransport::FirstJet);
    HermiteGalsCandidate rz(*mesh,{},radialSeed.primary());const auto rzMixed=closeHermiteGalsMixedAmr(*mesh,rz);
    rejects([&]{prepareHermiteGalsFirstJetScalar(*mesh,rz,rzMixed,{rz.identity(),centers(*mesh,radialProfile)},cfg);},"XY");
    std::cout<<"C2_MACRO_SCALAR_BINDING,PASS\n";
}
void macroPublicationTransport()
{
    auto mesh=makeMesh(4);const Profile f=[](Vector2D p){return Jet{(p.x-.5)*(p.x-.5)+(p.y-.5)*(p.y-.5)-.19*.19,2*(p.x-.5),2*(p.y-.5),0};};
    HermiteGalsPublicationConfig cfg;cfg.scalarRepresentation=HermiteGalsScalarRepresentation::C2MacroQuartic;
    const auto seed=stateFor(*mesh,f,false,{},[](int,Vector2D,Vector2D){return 0.;},HermiteGalsTransport::FirstJet);
    HermiteGalsCandidate candidate(*mesh,{},seed.primary());const auto mixed=closeHermiteGalsMixedAmr(*mesh,candidate);
    const auto initialId=candidate.identity();
    const auto initial=publishHermiteGalsFirstJetGeometry(*mesh,std::move(candidate),mixed,{initialId,centers(*mesh,f)},cfg);
    const HermiteGalsTimeLayer clock{1,0,.01,Position::End};
    const double factor=.99;
    const HermiteGalsResolvedMotion motion(initial.state().identity(),initial.state().timeLayer(),initial.state().identity(),clock,
        [=](Vector2D p){
            const Vector2D source{.5+factor*(p.x-.5),.5+factor*(p.y-.5)};
            const LevelSetBoundaryTrace::Jacobian jacobian{factor,0,0,factor};
            return LevelSetBoundaryTrace::JetResult{{source,.01,{},false,false},jacobian,1};});
    for(auto projection:{HermiteGalsTimeProjection::CharacteristicPullback,HermiteGalsTimeProjection::ClosedProjectionCompensated}) {
        const auto step=prepareHermiteGalsFirstJetStep(*mesh,initial,motion,{},cfg,HermiteGalsBoundaryTreatment::InflowHistory,projection);
        const auto sample=prepareHermiteGalsFirstJetSample(*mesh,initial,motion,{},cfg,HermiteGalsBoundaryTreatment::InflowHistory,projection);
        check(step.publication().state().geometry().snapshot().usesC2Macro()&&sample.geometry().snapshot().usesC2Macro(),
            "full macro transport rebuilt a different scalar family");
        check(step.publication().state().primary().phaseVolumes.empty(),"macro added transported inventory");
        for(int c=0;c<mesh->cellIndexNum();++c) {
            const auto bounds=step.publication().state().geometry().snapshot().nativeCellBounds(c);
            const Vector2D p=bounds[0]+(bounds[1]-bounds[0])*.37,departure{.5+factor*(p.x-.5),.5+factor*(p.y-.5)};
            const auto actual=step.publication().state().geometry().snapshot().scalarJetInNativeCell(c,p);
            const auto ref=f(departure);
            close(actual.value,ref.value,2e-12,"macro analytic actual value transport");
            close(actual.dx,factor*ref.dx,2e-11,"macro analytic actual dx transport");
            close(actual.dy,factor*ref.dy,2e-11,"macro analytic actual dy transport");
            const auto sampled=sample.geometry().snapshot().scalarJetInNativeCell(c,p);
            close(sampled.value,actual.value,2e-13,"publication/rate macro scalar consistency");
        }
        check(step.publication().state().timeLayer().step==clock.step&&step.publication().state().timeLayer().position==clock.position,
            "macro full transport lost target clock");
    }
    const auto old=initial.state().identity();auto wrong=cfg;wrong.scalarRepresentation=HermiteGalsScalarRepresentation::SingleQuartic;
    rejects([&]{prepareHermiteGalsFirstJetStep(*mesh,initial,motion,{},wrong);},"representations differ");
    check(initial.state().identity()==old&&initial.state().geometry().snapshot().usesC2Macro(),"failed macro transaction changed source");
    std::cout<<"C2_MACRO_PUBLICATION_TRANSPORT,PASS\n";
}
void sameFactoryTraceVolumes(const QuadtreeMesh& mesh,const HermiteGalsPublication& pub,
    const HermiteGalsPhaseVolumeConfig& config)
{
    // Public retained input remains independently rediscovered/certified.
    const auto reference=integrateHermiteGalsNativeTracePhaseVolumes(mesh,pub.state(),pub.domain(),
        pub.facePartitions(),pub.topology(),config);
    const auto& actual=pub.volumes();const auto aTime=actual.timeLayer(),bTime=reference.timeLayer();
    check(actual.sourceIdentity()==reference.sourceIdentity() && aTime.step==bTime.step
        && aTime.begin==bTime.begin && aTime.end==bTime.end && aTime.position==bTime.position
        && &actual.geometry().snapshot()==&reference.geometry().snapshot(),"factory/public trace source/full clock changed");
    const auto same=[](auto a,auto b){using T=decltype(a);
        if constexpr(sizeof(T)==sizeof(double))return std::memcmp(&a,&b,sizeof(T))==0;
        else return a==b || (std::isnan(a)&&std::isnan(b));};
    const auto vectors=[&](const auto& a,const auto& b){check(a.size()==b.size(),"factory/public trace vector size");
        for(std::size_t k=0;k<a.size();++k)check(same(a[k],b[k]),"factory/public trace ledger component changed");};
    const auto& av=actual.phaseVolumes();const auto& bv=reference.phaseVolumes();
    check(av.size()==bv.size() && actual.cells().size()==reference.cells().size(),"factory/public trace volume inventory changed");
    for(std::size_t k=0;k<av.size();++k){check(av[k].cell==bv[k].cell && same(av[k].liquid,bv[k].liquid),"factory/public phase volume changed");
        const auto a=actual.cells()[k],b=reference.cells()[k];
        check(a.cell==b.cell && a.certifiedPhase==b.certifiedPhase && same(a.capacity,b.capacity)
            && same(a.rawLiquid,b.rawLiquid),"factory/public raw native cell changed");}
    const auto& ar=actual.report();const auto& br=reference.report();
    check(ar.liquidCells==br.liquidCells && ar.gasCells==br.gasCells && ar.cutCells==br.cutCells
        && ar.components==br.components && ar.graphs==br.graphs && ar.faceRootIncidences==br.faceRootIncidences,
        "factory/public graph/phase inventory changed");
    const std::array<long double,10> ra{{ar.totalCapacity,ar.totalRawLiquid,ar.totalInitialLiquid,ar.relativeBoundaryClosure,
        ar.relativeMomentEstimate,ar.maxRelativeCapacityDefect,ar.maxRelativePureDefect,ar.totalSourceInventory,
        ar.geometricInventoryDifference,ar.maxCellRelativeInventoryDifference}},
        rb{{br.totalCapacity,br.totalRawLiquid,br.totalInitialLiquid,br.relativeBoundaryClosure,
        br.relativeMomentEstimate,br.maxRelativeCapacityDefect,br.maxRelativePureDefect,br.totalSourceInventory,
        br.geometricInventoryDifference,br.maxCellRelativeInventoryDifference}};
    vectors(ra,rb);
    const auto& al=actual.nativeTraceLedger();const auto& bl=reference.nativeTraceLedger();
    vectors(al.fixedRadialByCell,bl.fixedRadialByCell);vectors(al.freeRadialByCell,bl.freeRadialByCell);
    vectors(al.axialVolumeByCell,bl.axialVolumeByCell);vectors(al.radialVolumeByCell,bl.radialVolumeByCell);
    vectors(al.radialCapacityByCell,bl.radialCapacityByCell);
    check(same(al.exteriorFixedRadial,bl.exteriorFixedRadial) && same(al.fixedRadialClosure,bl.fixedRadialClosure)
        && same(al.maxRelativeGreenDefect,bl.maxRelativeGreenDefect) && same(al.maxRelativeCapacityDefect,bl.maxRelativeCapacityDefect)
        && al.sideGraphs==bl.sideGraphs && al.sideLiquidOwners==bl.sideLiquidOwners && al.zeroHalfFaces==bl.zeroHalfFaces,
        "factory/public complete trace ledger changed");
    std::cout<<"FACTORY_TRACE_REUSE,SAME_SOURCE_FULL_CLOCK,EXACT_VOLUMES_LEDGER_GRAPHS,PASS\n";
}
void nativeSideTransactions(bool quickReuse=false)
{
    for(bool radial:{false,true})for(bool amr:{false,true})for(int sign:{-1,1}){
        if(quickReuse && (amr || sign<0))continue;
        auto mesh=makeMesh(4,amr?1:0);const double line=amr?.5:.25;
        const Profile f=[=](Vector2D p){const double normal=(radial?p.y:p.x)-line,t=radial?p.x:p.y;
            const double factor=radial?1+t*t:1+t+t*t,derivative=radial?2*t:1+2*t;
            return Jet{sign*normal*factor,sign*(radial?normal*derivative:factor),sign*(radial?factor:normal*derivative),sign*derivative};};
        const HermiteGalsTimeLayer clock{7,.2,.3,Position::Midpoint};
        const auto seed=stateFor(*mesh,f,radial,clock,[](int,Vector2D,Vector2D){return 0.;},HermiteGalsTransport::FirstJet);
        const auto primary=seed.primary();const auto id=seed.identity();
        auto view=std::make_shared<const HermiteGalsGeometryView>(seed.geometry());
        std::shared_ptr<const Snapshot> source(view,&view->snapshot());
        const std::vector<Vector2D> contacts=radial?std::vector<Vector2D>{{1,line}}:std::vector<Vector2D>{{line,0},{line,1}};
        const auto topo=LevelSetInterfacePanels::discoverNativeTraceComponents(*mesh,source,radial,contacts);
        check(topo.curves.size()==1 && primary.phaseVolumes.empty(),"FirstJet side source inventory/component");
        for(const auto& g:topo.curves[0].graphs())check(g.nativeSide && g.nativeSide->usesSource(seed.geometry().snapshot()),"FirstJet side source changed");
        HermiteGalsCandidate candidate(*mesh,clock,primary);const auto nextId=candidate.identity();
        const auto mixed=closeHermiteGalsMixedAmr(*mesh,candidate);HermiteGalsPublicationConfig config;config.volume.boundaryContacts=contacts;
        auto invalid=config;invalid.volume.relativeLedgerTolerance=0;
        rejects([&]{publishHermiteGalsFirstJetGeometry(*mesh,std::move(candidate),mixed,{nextId,centers(*mesh,f)},invalid);},"relative tolerances");
        check(candidate.identity()==nextId && !candidate.hasGeometry() && candidate.primary().phaseVolumes.empty()
            && candidate.timeLayer().step==clock.step && candidate.timeLayer().time()==clock.time(),"failed publication consumed or advanced FirstJet candidate");
        for(std::size_t k=0;k<primary.vertices.size();++k){const auto a=primary.vertices[k],b=candidate.primary().vertices[k];
            check(a.vertex==b.vertex && a.phi==b.phi && a.dx==b.dx && a.dy==b.dy,"failed side publication changed primary vertex");}
        for(std::size_t k=0;k<primary.edges.size();++k){const auto a=primary.edges[k],b=candidate.primary().edges[k];
            check(a.vertices==b.vertices && a.phi==b.phi && a.dn==b.dn,"failed side publication changed primary edge");}
        check(seed.identity()==id && seed.timeLayer().time()==clock.time() && seed.primary().phaseVolumes.empty(),"Core topology changed accepted clock/state");
        // The default historical discovery/cache still rejects this source.
        rejects([&]{LevelSetInterfacePanels::discoverNativeComponents(*mesh,source,radial,contacts);},"coincident zero trace");
        const auto predicted=centers(*mesh,f);
        const auto pub=publishHermiteGalsFirstJetGeometry(*mesh,std::move(candidate),mixed,{nextId,predicted},config);
        const auto& scalar=pub.state().geometry().snapshot();const auto& ledger=pub.volumes().nativeTraceLedger();
        check(pub.state().identity()==nextId && pub.state().timeLayer().step==clock.step && pub.state().timeLayer().time()==clock.time()
            && pub.state().primary().phaseVolumes.empty() && pub.state().closure().cellBubbleCoefficients==predicted,"trace publication changed primary/clock/centers");
        check(&pub.facePartitions().source()==&scalar && &pub.volumes().geometry().snapshot()==&scalar
            && pub.volumes().sourceIdentity()==nextId && pub.volumes().timeLayer().time()==clock.time(),"trace publication source/revision/time mismatch");
        check(ledger.sideGraphs==pub.report().graphs && ledger.sideLiquidOwners==ledger.sideGraphs && ledger.zeroHalfFaces>0,"trace publication side ledger incomplete");
        long double peak=0,whole=0;std::size_t samples=0,gasOwners=0,hanging=0;
        for(const auto& cell:pub.volumes().cells()) {
            const auto bounds=scalar.nativeCellBounds(cell.cell);const auto lo=bounds[0],hi=bounds[1];
            const auto low=radial?lo.y:lo.x,high=radial?hi.y:hi.x;
            const long double length=std::max(0.,std::min(line,high)-low);
            const long double cross=radial?pi*(static_cast<long double>(hi.x)-lo.x)*(static_cast<long double>(hi.x)+lo.x):hi.y-lo.y;
            const auto expected=sign>0?cross*length:cell.capacity-cross*length;
            peak=std::max({peak,std::abs(cell.rawLiquid-expected)/cell.capacity,std::abs(ledger.axialVolumeByCell[cell.cell]-expected)/cell.capacity,
                std::abs(ledger.radialVolumeByCell[cell.cell]-expected)/cell.capacity});whole+=cell.rawLiquid;
        }
        for(std::size_t c=0;c<pub.topology().curves.size();++c) {
            const auto& chain=pub.topology().curves[c];check(chain.usesSource(scalar),"published side chain changed owner");
            for(std::size_t g=0;g<chain.graphs().size();++g) {
                const auto& graph=chain.graphs()[g];const auto& support=*graph.nativeSide;
                gasOwners+=pub.topology().cellPhase[support.ownerCell()]>0;
                const auto& edge=support.horizontal()?mesh->hEdges(support.edgeLabel()):mesh->vEdges(support.edgeLabel());
                hanging+=mesh->verts(edge.edgeToVert[0]).type()==TVertType::non || mesh->verts(edge.edgeToVert[1]).type()==TVertType::non;
                for(double t:{0.,.5,1.}) {
                    const auto sample=pub.sample(c,g,t);const auto raw=chain.frameAt(static_cast<int>(g),t);
                    check(sample.geometryIdentity==nextId && sample.timeLayer.time()==clock.time() && raw.point.x==sample.frame.point.x
                        && raw.point.y==sample.frame.point.y && raw.nativeCellIndex==sample.frame.nativeCellIndex,"published side sample lost provenance");++samples;
                }
                for(const auto& node:pub.quadrature(c,g))check(node.weight>0 && node.physicalWeight>0 && node.sample.geometryIdentity==nextId,"side quadrature source/weight");
            }
        }
        check(peak<=1e-12L && ledger.maxRelativeGreenDefect<=1e-12L && ledger.maxRelativeCapacityDefect<=1e-12L,"actual trace publication closed reference");
        if(radial)check(pub.report().axisEndpoints==1,"trace publication lost retained axis endpoint");
        check(samples==3*pub.report().graphs,"trace sample coverage");
        sameFactoryTraceVolumes(*mesh,pub,config.volume);
        const auto standalone=integrateHermiteGalsPhaseVolumes(*mesh,pub.state(),config.volume);
        close(standalone.report().totalRawLiquid,whole,1e-12,"standalone trace volume differs from retained publication cache");
        rejects([&]{integrateHermiteGalsNativeTracePhaseVolumes(*mesh,seed,pub.domain(),pub.facePartitions(),pub.topology(),config.volume);},"changed scalar owner");
        auto tampered=pub.topology();tampered.cellPhase[primary.geometryCells.front()]=0;
        rejects([&]{InterfaceCellMeasure::measureNativeTraceVolumes(*mesh,pub.domain(),pub.facePartitions(),tampered,contacts);},"complete same-source certificate");
        auto missing=pub.topology();missing.curves.clear();
        rejects([&]{InterfaceCellMeasure::measureNativeTraceVolumes(*mesh,pub.domain(),pub.facePartitions(),missing,contacts);},"complete same-source certificate");
        rejects([&]{pub.volumes().initializeCandidate(*mesh);},"step-0 time-0");
        rejects([&]{candidate.identity();},"consumed");
        std::cout<<"FIRST_JET_NATIVE_TRACE,"<<radial<<','<<amr<<','<<sign<<','<<pub.report().graphs<<','<<ledger.zeroHalfFaces<<','<<samples<<','<<gasOwners<<','<<hanging<<','<<peak<<",PASS\n";
        std::weak_ptr<const Snapshot> lifetime=source;source.reset();view.reset();check(!lifetime.expired(),"Core side source view lifetime");
        std::cout<<"NATIVE_SIDE_TRANSACTION,"<<radial<<','<<amr<<','<<sign<<",PASS\n";
    }
}
void nativeTraceExtended()
{
    int cases=0,cornerRejections=0;long double peak=0;std::size_t missingLiquid=0,joined=0,multiple=0,selected=0;
    for(bool radial:{false,true})for(int scenario=0;scenario<6;++scenario)for(int sign:{-1,1}) {
        if(radial && (scenario==0 || scenario==1))continue;
        auto mesh=makeMesh(4);Profile f;std::vector<Vector2D> contacts;long double liquid=0;
        if(scenario==0) {
            f=[=](Vector2D p){const double t=std::max(0.,p.y-.5);return Jet{sign*(p.x-.5+t*t),static_cast<double>(sign),sign*2*t,0};};
            contacts={{.5,0},{.25,1}};liquid=sign>0?11.L/24:13.L/24;
        } else if(scenario==1) {
            f=[=](Vector2D p){return Jet{sign*(p.x-.25)*(p.x-.75),sign*(2*p.x-1),0,0};};
            contacts={{.25,0},{.25,1},{.75,0},{.75,1}};liquid=.5L;
        } else {
            const double line=scenario==5?.5:(scenario==2?(radial?1:0):1);
            const bool axial=radial?(scenario==3 || scenario==4):scenario==4;
            f=[=](Vector2D p){if(radial && !axial)return Jet{sign*(p.x*p.x-line*line),sign*2*p.x,0,0};
                return Jet{sign*((axial?p.y:p.x)-line),static_cast<double>(axial?0:sign),static_cast<double>(axial?sign:0),0};};
            contacts=axial?(radial?std::vector<Vector2D>{{1,line}}:std::vector<Vector2D>{{0,line},{1,line}}):std::vector<Vector2D>{{line,0},{line,1}};
            const bool interiorNegative=sign*(.5-line)<0;
            liquid=interiorNegative?(radial?pi:1.L):0;
            if(scenario==5)liquid=sign>0?(radial?.25L*pi:.5L):0;
        }
        const HermiteGalsTimeLayer clock{2,.1,.2,Position::End};
        const auto seed=stateFor(*mesh,f,radial,clock,[](int,Vector2D,Vector2D){return 0.;},HermiteGalsTransport::FirstJet);
        auto p=seed.primary();
        if(scenario==5) {
            p.geometryCells.erase(std::remove_if(p.geometryCells.begin(),p.geometryCells.end(),[&](int c){return seed.geometry().snapshot().nativeCellBounds(c)[0].x>=.5;}),p.geometryCells.end());
            std::set<int> vertices;std::set<Key> edges;
            for(int c:p.geometryCells){const auto& cell=mesh->cells(mesh->cellIndexToLabel[c]);for(int k=0;k<4;++k){
                const int v=mesh->verts(cell.cellToVert[k]).index();vertices.insert(v);Key key{v,mesh->verts(cell.cellToVert[(k+1)%4]).index()};
                const auto a=point(*mesh,key[0]),b=point(*mesh,key[1]);if(a.y==b.y?a.x>b.x:a.y>b.y)std::swap(key[0],key[1]);edges.insert(key);}}
            p.vertices.erase(std::remove_if(p.vertices.begin(),p.vertices.end(),[&](const auto& v){return !vertices.count(v.vertex);}),p.vertices.end());
            p.edges.erase(std::remove_if(p.edges.begin(),p.edges.end(),[&](const auto& e){return !edges.count(e.vertices);}),p.edges.end());++selected;
        }
        HermiteGalsCandidate candidate(*mesh,clock,p);const auto id=candidate.identity();const auto mixed=closeHermiteGalsMixedAmr(*mesh,candidate);
        std::vector<double> center;for(int c:p.geometryCells)center.push_back(f(mesh->cells(mesh->cellIndexToLabel[c]).centerCoord()).value);
        HermiteGalsPublicationConfig config;config.volume.boundaryContacts=contacts;
        if(radial && (scenario==3 || scenario==4)) {
            // Existing Core explicitly excludes an axis endpoint at a physical
            // domain corner. Preserve its boundary-identity rejection rather
            // than expanding topology scope in this publication task.
            rejects([&]{publishHermiteGalsFirstJetGeometry(*mesh,std::move(candidate),mixed,{id,center},config);},"undeclared active boundary");
            check(candidate.identity()==id && !candidate.hasGeometry() && candidate.timeLayer().time()==clock.time(),"uncertified axis corner consumed candidate");
            ++cornerRejections;std::cout<<"FIRST_JET_NATIVE_TRACE_AXIS_CORNER_REJECT,"<<scenario<<','<<sign<<",PASS\n";continue;
        }
        const auto pub=publishHermiteGalsFirstJetGeometry(*mesh,std::move(candidate),mixed,{id,center},config);
        const auto& ledger=pub.volumes().nativeTraceLedger();const auto cap=pub.volumes().report().totalCapacity;
        const auto trial=std::abs(pub.volumes().report().totalRawLiquid-liquid)/cap;peak=std::max(peak,trial);
        check(trial<=1e-12L && ledger.maxRelativeGreenDefect<=1e-12L,"extended trace publication closed reference");
        check(pub.state().identity()==id && pub.state().primary().phaseVolumes.empty() && pub.state().closure().cellBubbleCoefficients==center,"extended trace modified primary");
        if(scenario==0){bool side=false,ordinary=false;for(const auto& chain:pub.topology().curves)for(const auto& g:chain.graphs())(g.nativeSide?side:ordinary)=true;
            check(side && ordinary && pub.topology().curves.size()==1,"actual publication ordinary/side join missing");++joined;}
        if(scenario==1){check(pub.topology().curves.size()==2,"actual side multiple components missing");++multiple;}
        missingLiquid+=ledger.sideGraphs-ledger.sideLiquidOwners;
        std::cout<<"FIRST_JET_NATIVE_TRACE_EXTRA,"<<radial<<','<<scenario<<','<<sign<<','<<ledger.sideGraphs<<','<<ledger.sideLiquidOwners<<','<<trial<<",PASS\n";++cases;
    }
    check(missingLiquid>0 && joined==2 && multiple==2 && selected==4 && cornerRejections==4,"extended zero trace coverage incomplete");
    std::cout<<"FIRST_JET_NATIVE_TRACE_EXTRA_SUMMARY,"<<cases<<','<<joined<<','<<multiple<<','<<selected<<','<<missingLiquid<<','<<peak<<','<<cornerRejections<<",PASS\n";
}
void sourceAndCuts(const QuadtreeMesh& mesh,const HermiteGalsPublication& pub)
{
    pub.validateMesh(mesh);const auto& scalar=pub.state().geometry().snapshot();
    check(&pub.volumes().geometry().snapshot()==&scalar && &pub.facePartitions().source()==&scalar,"different scalar allocation");
    check(pub.volumes().sourceIdentity()==pub.state().identity(),"volume revision mismatch");
    for(const auto& c:pub.topology().curves)check(c.usesSource(scalar),"curve source mismatch");
    std::size_t count=0;
    for(const auto face:pub.domain().faces(mesh))if(face.negativeCell>=0 || face.positiveCell>=0)for(int half=0;half<2;++half){
        const auto& cut=pub.facePartitions().halfFace(face.nativeRow,half);check(cut.usesSource(scalar),"cut source mismatch");
        long double length=0;for(const auto& range:cut.ranges()){check(range.length>0 && range.offset>=0,"invalid retained range");length+=range.length;}
        close(length,1,2e-13,"retained range coverage");
        for(const auto& root:cut.roots())check(root.nativeCellIndex>=0 && root.lower<=root.parameter && root.parameter<=root.upper,"root provenance/enclosure");++count;}
    check(count==pub.report().halfFaces,"half-face count mismatch");
}
void analytic(bool radial,bool amr,int sign,bool multiple)
{
    auto mesh=makeMesh(8,amr?1:0);const double r=multiple?.13:.27;
    const std::vector<Vector2D> centersXY=multiple?std::vector<Vector2D>{{.31,.49},{.73,.49}}:std::vector<Vector2D>{{.47,.51}};
    const std::vector<Vector2D> centersRZ=multiple?std::vector<Vector2D>{{0,.31},{0,.73}}:std::vector<Vector2D>{{0,.51}};
    const auto origin=radial?centersRZ:centersXY;
    const Profile f=[=](Vector2D p){const auto a=p-origin[0];const double A=a.x*a.x+a.y*a.y-r*r;
        if(!multiple)return Jet{sign*A,sign*2*a.x,sign*2*a.y,0};
        const auto b=p-origin[1];const double B=b.x*b.x+b.y*b.y-r*r;
        return Jet{sign*A*B,sign*(2*a.x*B+2*b.x*A),sign*(2*a.y*B+2*b.y*A),sign*4*(a.x*b.y+b.x*a.y)};};
    const auto p=measuredPrimary(*mesh,f,radial);const auto pub=publish(*mesh,f,p,{1,0,.01,Position::Midpoint});sourceAndCuts(*mesh,pub);
    check(pub.topology().curves.size()==origin.size(),"analytic component count");
    long double normalError=0,curvatureError=0;std::size_t nodes=0,axes=0;
    const auto sample=[&](const HermiteGalsSurfaceSample& s){
        check(s.geometryIdentity==pub.state().identity() && s.timeLayer.position==Position::Midpoint && s.timeLayer.time()==.005,"sample identity/time");
        const auto& frame=s.frame;auto center=origin[0];double best=1e100;
        for(const auto c:origin){const auto d=frame.point-c;const double error=std::abs(std::hypot(d.x,d.y)-r);if(error<best){best=error;center=c;}}
        const auto d=frame.point-center;const auto n=d*(sign/r);
        normalError=std::max(normalError,close(std::hypot(frame.normal.x-n.x,frame.normal.y-n.y),0,2e-10,"analytic normal"));
        curvatureError=std::max(curvatureError,close(frame.curvature,sign/r,2e-8,"analytic meridional curvature"));
        curvatureError=std::max(curvatureError,close(s.curvature,sign*(radial?2:1)/r,2e-8,"analytic physical curvature"));
        const auto jet=pub.state().geometry().snapshot().scalarJetAtNativeCellChart(frame.nativeCellIndex,frame.cellChartPoint);
        close(jet.dx/std::hypot(jet.dx,jet.dy),frame.normal.x,2e-13,"normal source gradient x");
        close(jet.dy/std::hypot(jet.dx,jet.dy),frame.normal.y,2e-13,"normal source gradient y");
        if(radial && frame.point.x==0){check(frame.normal.x==0 && s.azimuthalCurvature==frame.curvature && s.curvature==2*frame.curvature,"axis even limit");++axes;}
    };
    for(std::size_t c=0;c<pub.topology().curves.size();++c)for(std::size_t g=0;g<pub.topology().curves[c].graphs().size();++g){
        for(double t:{0.,.37,1.})sample(pub.sample(c,g,t));
        for(const auto& q:pub.quadrature(c,g)){sample(q.sample);check(q.weight>0 && q.physicalWeight>0,"nonpositive quadrature");++nodes;}}
    check(axes==(radial?2*origin.size():0) && axes==pub.report().axisEndpoints,"axis endpoint count");
    check(nodes==pub.report().quadratureNodes && nodes>0,"quadrature count");
    const long double cap=radial?pi:1,inner=origin.size()*(radial?4*pi*r*r*r/3:pi*r*r);
    const auto volumeError=close(pub.volumes().report().totalRawLiquid/cap,(sign>0?inner:cap-inner)/cap,1e-10,"analytic total volume");
    for(std::size_t k=0;k<p.phaseVolumes.size();++k)check(pub.state().primary().phaseVolumes[k].liquid==p.phaseVolumes[k].liquid,"publication changed M");
    check(pub.report().maxRelativeInventoryResidual<=1e-10,"inventory gate");
    auto copy=pub;check(&copy.state().geometry().snapshot()==&pub.state().geometry().snapshot(),"copy cloned scalar");
    auto moved=std::move(copy);rejects([&]{copy.state();},"moved-from");sourceAndCuts(*mesh,moved);
    rejects([&]{pub.sample(999,0,.5);},"sample component");rejects([&]{pub.sample(0,0,-.1);},"sample component");
    if(sign>0 && !multiple){const auto previous=publish(*mesh,f,p);std::vector<HermiteGalsFaceVelocity> velocity;
        for(const auto face:pub.domain().faces(*mesh))if(face.negativeCell>=0 || face.positiveCell>=0){HermiteGalsFaceVelocity v;v.nativeFace=face.nativeRow;
            v.boundary=face.negativeCell>=0 && face.positiveCell>=0?HermiteGalsFaceVelocity::Boundary::Interior:HermiteGalsFaceVelocity::Boundary::Impermeable;velocity.push_back(v);}
        const HermiteGalsFluxVelocity zero(*mesh,pub.state(),velocity);
        const auto flux=advanceHermiteGalsPhaseVolumes(*mesh,previous.state(),pub.state(),zero);
        check(&flux.geometry().snapshot()==&pub.state().geometry().snapshot() && flux.geometryIdentity()==pub.state().identity(),"phase flux source mismatch");
        for(std::size_t k=0;k<p.phaseVolumes.size();++k)check(flux.phaseVolumes()[k].liquid==p.phaseVolumes[k].liquid,"zero velocity changed M");
        std::cout<<"ZERO_FLOW,"<<(radial?"RZ":"XY")<<','<<amr<<",PASS\n";}
    std::cout<<"PUBLICATION,"<<(radial?"RZ":"XY")<<','<<amr<<','<<sign<<','<<origin.size()<<','<<nodes<<','<<axes<<','<<normalError<<','<<curvatureError<<','<<volumeError<<','<<pub.report().maxRelativeInventoryResidual<<'\n';
}
// Independent analytic disk/rectangle integral. Breakpoints come from the
// circle's intersections with horizontal bounds, without production geometry.
long double tangentDiskVolume(Vector2D lo,Vector2D hi,Vector2D center,double radius,bool radial)
{
    const long double cx=center.x,cy=center.y,r=radius;
    const long double a=std::max<long double>(lo.x,cx-r),b=std::min<long double>(hi.x,cx+r);
    if(b<=a)return 0;
    std::vector<long double> points{a,b};
    for(long double y:{static_cast<long double>(lo.y),static_cast<long double>(hi.y)})
        if(std::abs(y-cy)<r){const auto u=std::sqrt(r*r-(y-cy)*(y-cy));for(auto x:{cx-u,cx+u})if(x>a && x<b)points.push_back(x);}
    std::sort(points.begin(),points.end());long double sum=0;
    const auto primitive=[&](long double x){const auto u=std::clamp(x-cx,-r,r),s=std::sqrt(std::max(0.L,r*r-u*u));
        const auto disk=.5L*(u*s+r*r*std::atan2(u,s));return radial?cx*disk-s*s*s/3:disk;};
    for(std::size_t i=1;i<points.size();++i){const auto x0=points[i-1],x1=points[i],xm=(x0+x1)/2;
        const auto s=std::sqrt(std::max(0.L,r*r-(xm-cx)*(xm-cx)));
        if(std::min<long double>(hi.y,cy+s)<=std::max<long double>(lo.y,cy-s))continue;
        long double constant=0,k=0;
        if(cy+s<hi.y){constant+=cy;k+=1;}else constant+=hi.y;
        if(cy-s>lo.y){constant-=cy;k+=1;}else constant-=lo.y;
        sum+=constant*(radial?(x1*x1-x0*x0)/2:x1-x0)+k*(primitive(x1)-primitive(x0));}
    return radial?2*pi*sum:sum;
}
// One full native-cell local gauge and the same complete graphs at global
// origin. Reference calls keep their original independent quadrature paths.
void nativeCellDualMoments()
{
    using I=LevelSetMappedInterface;
    using Trace=LevelSetMappedSegment::NativeQuarticTrace;
    for(int scene=0;scene<3;++scene) {
        const bool radial=scene!=0,amr=scene==2;
        auto mesh=makeMesh(4,amr?1:0);
        const Vector2D center=radial?Vector2D{0,.51}:Vector2D{.47,.51};
        constexpr double radius=.27;
        const Profile f=[=](Vector2D p){const auto d=p-center;
            return Jet{d.x*d.x+d.y*d.y-radius*radius,2*d.x,2*d.y,0};};
        const HermiteGalsTimeLayer clock{7,.2,.3,Position::Midpoint};
        const auto seed=stateFor(*mesh,f,radial,clock,[](int,Vector2D,Vector2D){return 0.;},HermiteGalsTransport::FirstJet);
        HermiteGalsCandidate candidate(*mesh,clock,seed.primary());const auto id=candidate.identity();
        const auto mixed=closeHermiteGalsMixedAmr(*mesh,candidate);
        const auto predictedCenters=centers(*mesh,f);
        // Same actual master/mixed/center data, independently reconstructed
        // before the candidate is consumed. Old reference roots never share
        // the publication source's initial-root cache or guarded policy.
        const auto referenceSource=scalarFor(*mesh,patchFor(*mesh,candidate,mixed,predictedCenters),radial);
        const auto pub=publishHermiteGalsFirstJetGeometry(*mesh,std::move(candidate),mixed,{id,predictedCenters});
        sourceAndCuts(*mesh,pub);
        // Reuse only this immutable publication's derived binding. Its
        // actual live-mesh identity still rejects physical coordinate edits
        // for both the original handle and an ordinary handle copy.
        const auto copied=pub;const auto* boundScalar=&pub.state().geometry().snapshot();
        pub.validateMesh(*mesh);copied.validateMesh(*mesh);
        auto& actualCell=mesh->cells(mesh->cellIndexToLabel[pub.state().primary().geometryCells.front()]);
        const auto savedCenter=actualCell.centerCoord();
        actualCell.setCenterCoord(savedCenter+Vector2D{.001,0});
        rejects([&]{pub.validateMesh(*mesh);},"NativeMeshIdentity");
        rejects([&]{copied.validateMesh(*mesh);},"NativeMeshIdentity");
        actualCell.setCenterCoord(savedCenter);
        pub.validateMesh(*mesh);copied.validateMesh(*mesh);
        const int actualVertexIndex=pub.state().primary().vertices.front().vertex;
        auto& actualVertex=mesh->verts(mesh->vertIndexToLabel[actualVertexIndex]);
        const auto savedCoordinate=actualVertex.coord();
        actualVertex.setCoord(savedCoordinate+Vector2D{.001,0});
        rejects([&]{pub.validateMesh(*mesh);},"NativeMeshIdentity");
        rejects([&]{copied.validateMesh(*mesh);},"NativeMeshIdentity");
        actualVertex.setCoord(savedCoordinate);
        pub.validateMesh(*mesh);copied.validateMesh(*mesh);
        check(&pub.state().geometry().snapshot()==boundScalar
            && &copied.state().geometry().snapshot()==boundScalar,
            "restored publication validation changed the actual immutable scalar owner");
        auto view=std::make_shared<const HermiteGalsGeometryView>(pub.state().geometry());
        const std::shared_ptr<const Snapshot> source(view,&view->snapshot());
        const auto* originalScalar=&pub.state().geometry().snapshot();
        std::vector<std::vector<I::Graph>> byCell(mesh->cellIndexNum());
        std::array<bool,2> foundOrientation{{false,false}};
        std::set<int> selected;
        for(const auto& chain:pub.topology().curves) {
            check(chain.usesSource(*source),"dual moments changed actual graph source");
            for(const auto& graph:chain.graphs())if(!graph.nativeSide) {
                byCell.at(graph.nativeCellIndex).push_back(graph);
                const int orientation=graph.coordinateIsX?0:1;
                if(!foundOrientation[orientation]) {
                    foundOrientation[orientation]=true;selected.insert(graph.nativeCellIndex);
                }
            }
        }
        check(foundOrientation[0] && foundOrientation[1],"dual moment fixture missed a coordinate orientation");
        bool hangingCut=false;
        if(amr)for(int cell:pub.state().primary().geometryCells)if(!byCell[cell].empty()) {
            const auto& native=mesh->cells(mesh->cellIndexToLabel[cell]);bool hanging=false;
            for(int vertex:native.cellToVert)hanging|=mesh->verts(vertex).supCell>=0;
            for(int direction=0;direction<2;++direction)
                hanging|=native.cellToVEdge[direction][1]>=0 || native.cellToHEdge[direction][1]>=0;
            if(hanging){selected.insert(cell);hangingCut=true;break;}
        }
        check(!amr || hangingCut,"dual moment AMR fixture missed actual hanging cut support");
        std::size_t graphsChecked=0;std::array<std::size_t,2> orientationChecks{{0,0}};
        long double peakGlobal=0,peakLocal=0,peakGreen=0,peakAnalytic=0;
        for(int cell:selected) {
            check(referenceSource.get()!=source.get() && referenceSource->identicalNativeQuarticCell(cell,*source,cell),
                "dual moment independent reference changed the actual exact native Q4 cell");
            const auto bounds=source->nativeCellBounds(cell);const auto lo=bounds[0],hi=bounds[1];
            const Vector2D greenOrigin{hi.x,lo.y};
            const long double capacity=(static_cast<long double>(hi.x)-lo.x)*(static_cast<long double>(hi.y)-lo.y)
                *(radial?pi*(static_cast<long double>(lo.x)+hi.x):1);
            const I::Region region{0,{{lo,{hi.x,lo.y},hi,{lo.x,hi.y}}}};
            const std::vector<I::NativeMomentTarget> targets{{{region},greenOrigin}};
            const auto makeCuts=[&](const std::shared_ptr<const Snapshot>& owner) {
                std::vector<Trace> cuts;cuts.reserve(4);
                for(int edge=0;edge<4;++edge) {
                    auto a=region.vertices[edge],b=region.vertices[(edge+1)%4];
                    if(a.x>b.x || (a.x==b.x && a.y>b.y))std::swap(a,b);
                    cuts.push_back(LevelSetMappedSegment::nativeQuarticTrace(owner,a,b,a.x==b.x?0:1));
                    check(cuts.back().usesSource(*owner),"dual moment face cut changed actual source");
                }
                return cuts;
            };
            // Rebuild independent move-only traces from identical actual inputs.
            // Each consumer performs its own mandatory face-root refinements.
            const auto retry=[&](const std::shared_ptr<const Snapshot>& owner,auto invoke) {
                auto working=makeCuts(owner);
                for(int attempt=0;attempt<=256;++attempt) {
                    std::vector<const Trace*> refs;for(const auto& cut:working)refs.push_back(&cut);
                    try{return invoke(refs);}
                    catch(const I::RootResolutionError& failure) {
                        check(attempt<256 && !failure.requests.empty(),"dual moment face-root refinement limit");
                        std::set<std::size_t> replaced;
                        for(const auto& request:failure.requests) {
                            check(request.segmentIndex<working.size(),"dual moment foreign refinement segment");
                            if(replaced.insert(request.segmentIndex).second) {
                                auto& cut=working[request.segmentIndex];
                                check(request.rootIndex<cut.events().size(),"dual moment foreign refinement root");
                                cut=cut.refined(request.rootIndex);
                            }
                        }
                    }
                }
                throw std::runtime_error("unreachable dual moment refinement exhaustion");
            };
            std::vector<I::GraphMoments> global;
            const auto actual=retry(source,[&](const auto& refs){return I::integrateNativeQuarticRegionMomentsBatch(
                *source,cell,byCell[cell],refs,targets,radial,&global);});
            const auto reference=retry(referenceSource,[&](const auto& refs){return I::integrateNativeQuarticRegionMoments(
                *referenceSource,cell,byCell[cell],refs,{region},radial,greenOrigin);});
            check(actual.size()==1 && actual[0].regions.size()==1 && reference.regions.size()==1
                && global.size()==byCell[cell].size(),"dual moment full-cell/global inventory changed");
            for(std::size_t g=0;g<global.size();++g) {
                const auto old=I::nativeMoments(*referenceSource,byCell[cell][g],radial);const auto fresh=global[g];
                peakGlobal=std::max({peakGlobal,
                    close(fresh.axial,old.axial,2e-12L*capacity,"dual global axial moment")/capacity,
                    close(fresh.radial,old.radial,2e-12L*capacity,"dual global radial moment")/capacity,
                    close(fresh.errorEstimate,old.errorEstimate,2e-12L*capacity,"dual global moment estimate")/capacity});
                check(std::isfinite(fresh.errorEstimate) && fresh.errorEstimate>=0 && fresh.errorEstimate<=2e-12L*capacity
                    && old.errorEstimate>=0 && old.errorEstimate<=2e-12L*capacity,"dual global original controlled-moment budget");
                ++orientationChecks[byCell[cell][g].coordinateIsX?0:1];++graphsChecked;
            }
            const auto& a=actual[0].regions[0];const auto& b=reference.regions[0];
            check(a.regionId==0 && b.regionId==0 && a.freeArcs>0 && b.freeArcs>0,"dual local free-arc owner lost");
            const std::array<long double,10> av{{a.capacity,a.axialLiquid,a.radialLiquid,a.axialGas,a.radialGas,
                a.fixedAxialLiquid,a.fixedRadialLiquid,a.freeAxialLiquid,a.freeRadialLiquid,a.errorEstimate}},
                bv{{b.capacity,b.axialLiquid,b.radialLiquid,b.axialGas,b.radialGas,
                b.fixedAxialLiquid,b.fixedRadialLiquid,b.freeAxialLiquid,b.freeRadialLiquid,b.errorEstimate}};
            for(std::size_t k=0;k<av.size();++k)
                peakLocal=std::max(peakLocal,close(av[k],bv[k],1e-10L*capacity,"dual local original moment")/capacity);
            const auto disk=tangentDiskVolume(lo,hi,center,radius,radial);
            for(const auto* r:{&a,&b}) {
                check(std::isfinite(r->capacity) && r->capacity>0 && r->axialLiquid>=0 && r->radialLiquid>=0
                    && r->axialGas>=0 && r->radialGas>=0 && r->axialLiquid<=r->capacity && r->radialLiquid<=r->capacity
                    && r->axialGas<=r->capacity && r->radialGas<=r->capacity,"dual local strict capacity; no clipping");
                check(std::isfinite(r->errorEstimate) && r->errorEstimate>=0 && r->errorEstimate<=1e-10L*capacity,
                    "dual local original moment estimate budget");
                peakGreen=std::max({peakGreen,
                    close(r->capacity,capacity,1e-10L*capacity,"dual local physical capacity")/capacity,
                    close(r->axialLiquid,r->radialLiquid,1e-10L*capacity,"dual local liquid Green ledger")/capacity,
                    close(r->axialGas,r->radialGas,1e-10L*capacity,"dual local gas Green ledger")/capacity,
                    close(r->axialLiquid+r->axialGas,capacity,1e-10L*capacity,"dual local axial closure")/capacity,
                    close(r->radialLiquid+r->radialGas,capacity,1e-10L*capacity,"dual local radial closure")/capacity});
                peakAnalytic=std::max({peakAnalytic,
                    close(r->axialLiquid,disk,1e-10L*capacity,"dual local analytic liquid")/capacity,
                    close(r->axialGas,capacity-disk,1e-10L*capacity,"dual local analytic gas")/capacity});
            }
            close(a.fixedAxialLiquid,actual[0].exteriorFixedAxialLiquid+actual[0].coincidentAxialLiquid,
                1e-10L*capacity,"dual local fixed axial assembly");
            close(a.fixedRadialLiquid,actual[0].exteriorFixedRadialLiquid+actual[0].coincidentRadialLiquid,
                1e-10L*capacity,"dual local fixed radial assembly");
        }
        const auto retainedClock=pub.state().timeLayer();
        check(pub.state().identity()==id && &pub.state().geometry().snapshot()==originalScalar && pub.state().primary().phaseVolumes.empty()
            && retainedClock.step==clock.step && retainedClock.begin==clock.begin && retainedClock.end==clock.end
            && retainedClock.position==clock.position && pub.volumes().sourceIdentity()==id
            && pub.volumes().timeLayer().step==clock.step && pub.volumes().timeLayer().begin==clock.begin
            && pub.volumes().timeLayer().end==clock.end && pub.volumes().timeLayer().position==clock.position,
            "dual moment test changed actual source/full clock");
        check(orientationChecks[0]>0 && orientationChecks[1]>0,"dual moments missed retained endpoint recovery orientation");
        std::cout<<"NATIVE_CELL_DUAL_MOMENTS,"<<(radial?"RZ":"XY")<<','<<amr<<','<<selected.size()<<','<<graphsChecked
            <<','<<orientationChecks[0]<<','<<orientationChecks[1]<<','<<hangingCut<<','<<peakGlobal<<','<<peakLocal
            <<','<<peakGreen<<','<<peakAnalytic<<",PASS\n";
    }
}
void isolatedContacts()
{
    for(bool radial:{false,true})for(bool amr:{false,true})for(int sign:{-1,1}) {
        const std::vector<Vector2D> origins=radial?std::vector<Vector2D>{{0,.53125},{0,.5}}:
            std::vector<Vector2D>{{.5,.53125},{.53125,.5},{.5,.5},{.5-std::ldexp(1.,-12),.53125},{.5+std::ldexp(1.,-12),.53125}};
        for(std::size_t fixture=0;fixture<origins.size();++fixture)try {
            auto mesh=makeMesh(4,amr?1:0);const auto center=origins[fixture];constexpr double r=.25;
            check(amr?mesh->cellIndexNum()>16:mesh->cellIndexNum()==16,"contact AMR fixture was not refined");
            const Profile f=[=](Vector2D p){const auto d=p-center;return Jet{sign*(d.x*d.x+d.y*d.y-r*r),sign*2*d.x,sign*2*d.y,0};};
            const auto primary=measuredPrimary(*mesh,f,radial);
            const auto pub=publish(*mesh,f,primary,{1,0,.01,Position::Midpoint});sourceAndCuts(*mesh,pub);
            check(pub.topology().curves.size()==1 && pub.topology().curves[0].closed()==!radial,"contact curve count/closure");
            check(pub.topology().negativeRegions.size()==1,"contact spuriously split/joined global liquid regions");
            check(pub.topology().negativeRegions[0].touchesNonAxisBoundary==(sign<0),"contact changed liquid boundary ownership");
            long double maxCellError=0,maxNormalError=0,maxCurvatureError=0,length=0,area=0;int contactSamples=0;
            for(const auto& cell:pub.volumes().cells()) {
                const auto& c=mesh->cells(mesh->cellIndexToLabel[cell.cell]);
                const auto lo=mesh->verts(c.cellToVert[0]).coord(),hi=mesh->verts(c.cellToVert[2]).coord();
                const auto disk=tangentDiskVolume(lo,hi,center,r,radial);
                const auto reference=sign>0?disk:cell.capacity-disk;
                maxCellError=std::max(maxCellError,close(cell.rawLiquid,reference,1e-10L*cell.capacity,"contact cell volume"));
            }
            const auto sample=[&](const HermiteGalsSurfaceSample& s,bool endpoint){
                const auto d=s.frame.point-center;const auto n=d*(sign/r);
                maxNormalError=std::max(maxNormalError,close(std::hypot(s.frame.normal.x-n.x,s.frame.normal.y-n.y),0,2e-10,"contact analytic normal"));
                maxCurvatureError=std::max(maxCurvatureError,close(s.curvature,sign*(radial?2:1)/r,2e-8,"contact analytic curvature"));
                check(s.geometryIdentity==pub.state().identity() && s.timeLayer.position==Position::Midpoint,"contact sample source/time");
                if(endpoint && ((std::abs(d.y)<1e-13 && std::abs(std::abs(d.x)-r)<1e-13)
                    || (!radial && std::abs(d.x)<1e-13 && std::abs(std::abs(d.y)-r)<1e-13)))++contactSamples;
            };
            const auto& curve=pub.topology().curves[0];
            for(std::size_t g=0;g<curve.graphs().size();++g) {
                for(double t:{0.,.37,1.})sample(pub.sample(0,g,t),t==0 || t==1);
                for(const auto& q:pub.quadrature(0,g)){sample(q.sample,false);check(q.weight>0 && q.physicalWeight>0,"contact positive quadrature");}
            }
            if(fixture<3)check(contactSamples>=2,"complete curve did not consume contact endpoints");
            // Publication retains one GL8 rule per graph; its contract has no
            // adaptive arclength bound. Use the existing controlled material
            // panel consumer of the SAME chain for length/area verification.
            const LevelSetInterfacePanels material(curve);
            for(const auto& panel:material.panels()){length+=panel.measure;area+=panel.physicalMeasure;}
            close(length,radial?pi*r:2*pi*r,2e-9,"contact curve arclength");
            close(area,radial?4*pi*r*r:2*pi*r,2e-9,"contact physical area");
            const long double cap=radial?pi:1,inner=radial?4*pi*r*r*r/3:pi*r*r;
            const auto volumeError=close(pub.volumes().report().totalRawLiquid,sign>0?inner:cap-inner,1e-10L*cap,"contact total volume");
            for(std::size_t k=0;k<primary.phaseVolumes.size();++k)check(pub.state().primary().phaseVolumes[k].liquid==primary.phaseVolumes[k].liquid,"contact publication changed M");
            std::cout<<"ISOLATED_CONTACT,"<<(radial?"RZ":"XY")<<','<<amr<<','<<sign<<','<<fixture<<','<<curve.graphs().size()<<','<<contactSamples
                <<','<<maxCellError<<','<<volumeError<<','<<maxNormalError<<','<<maxCurvatureError<<",PASS\n";
        }catch(const std::exception& e){throw std::runtime_error("isolated contact radial="+std::to_string(radial)+", AMR="+std::to_string(amr)
            +", sign="+std::to_string(sign)+", fixture="+std::to_string(fixture)+": "+e.what());}
    }
    auto mesh=makeMesh(4);const Profile singular=[](Vector2D p){
        const auto x=p.x-.5,y=p.y-.5;return Jet{x*x+y*y,2*x,2*y,0};};
    rejects([&]{measuredPrimary(*mesh,singular,false);},"contact lacks a nonzero original normal derivative certificate");
}
void nonbinaryContacts()
{
    for(int power:{2,4})for(bool radial:{false,true})for(bool amr:{false,true})for(int sign:{-1,1})try {
        auto mesh=makeMesh(8,amr?1:0);const double cx=radial?0:.5;constexpr double r=.25;
        const Profile f=[=](Vector2D p){const auto x=p.x-cx,u=3*p.y-1;
            return power==2?Jet{sign*(x*x+u*u-r*r),sign*2*x,sign*6*u,0}:
                Jet{sign*(x*x+u*u*u*u-r*r),sign*2*x,sign*12*u*u*u,0};};
        const auto primary=measuredPrimary(*mesh,f,radial);const auto pub=publish(*mesh,f,primary);sourceAndCuts(*mesh,pub);
        check(pub.topology().curves.size()==1 && pub.topology().curves[0].closed()==!radial
            && pub.topology().negativeRegions.size()==1,"nonbinary contact topology");
        const auto& curve=pub.topology().curves[0];long double normalError=0,curvatureError=0;int nearSamples=0,contacts=0;
        for(std::size_t g=0;g<curve.graphs().size();++g) {
            for(double t:{0.,std::ldexp(1.,-48),.37,1-std::ldexp(1.,-48),1.}) {
                const auto s=pub.sample(0,g,t);const auto x=s.frame.point.x-cx,u=3*s.frame.point.y-1;
                const double px=2*x,py=power==2?6*u:12*u*u*u,pyy=power==2?18:108*u*u;
                const auto norm=std::hypot(px,py);
                normalError=std::max(normalError,close(std::hypot(s.frame.normal.x-sign*px/norm,s.frame.normal.y-sign*py/norm),0,2e-10,"nonbinary original normal"));
                const auto k=sign*(2*py*py+pyy*px*px)/(norm*norm*norm)+(radial?sign*2/norm:0);
                curvatureError=std::max(curvatureError,close(s.curvature,k,2e-8,"nonbinary original curvature"));
                if(t>0 && t<1 && t!=.37)++nearSamples;
                if((t==0 || t==1) && std::abs(u)<1e-13 && std::abs(std::abs(x)-r)<1e-13) {
                    ++contacts;
                    if(power==2 && !radial && !amr && sign>0 && t==0) {
                        // Force sampling inside the original nonbinary root
                        // enclosure. Refinement must retain the multiple root;
                        // only the final binary64 chart representation fails.
                        rejects([&]{pub.sample(0,g,std::ldexp(1.,-60));},"not representably interior");
                        check(pub.sample(0,g,0).geometryIdentity==pub.state().identity(),"failed contact refinement mutated source");
                    }
                }
            }
        }
        check(contacts>=2 && nearSamples>0,"nonbinary contact/refinement samples missing");
        // Independent affine ellipse volume or beta integral of a quartic oval.
        const long double inner=power==2?(radial?4*pi*r*r*r/9:pi*r*r/3):
            radial?8*pi*std::pow(static_cast<long double>(r),2.5L)/15:
                std::pow(static_cast<long double>(r),1.5L)*std::tgamma(.25L)*std::tgamma(1.5L)/(3*std::tgamma(1.75L));
        const long double cap=radial?pi:1;
        const auto error=close(pub.volumes().report().totalRawLiquid,sign>0?inner:cap-inner,1e-10L*cap,"nonbinary analytic volume");
        std::cout<<"NONBINARY_CONTACT,"<<power<<','<<(radial?"RZ":"XY")<<','<<amr<<','<<sign<<','<<contacts<<','<<nearSamples
            <<','<<error<<','<<normalError<<','<<curvatureError<<",PASS\n";
    }catch(const std::exception& e){throw std::runtime_error("nonbinary contact power="+std::to_string(power)+", radial="+std::to_string(radial)
        +", AMR="+std::to_string(amr)+", sign="+std::to_string(sign)+": "+e.what());}
}
void pureAndLifetime()
{
    const auto retained=[] {
        auto mesh=makeMesh(8);const Profile f=[](Vector2D p){const auto x=p.x-.47,y=p.y-.51;return Jet{x*x+y*y-.27*.27,2*x,2*y,0};};
        return publish(*mesh,f,measuredPrimary(*mesh,f,false));
    }();
    close(retained.sample(0,0,.4).curvature,1/.27,2e-8,"surface after mesh/input destruction");
    check(!retained.quadrature(0,0).empty(),"lost retained surface quadrature");
    for(bool radial:{false,true})for(int sign:{-1,1}){
        const auto pub=[&]{auto mesh=makeMesh(4,1);const auto domain=domainFor(*mesh,radial);const Profile f=[=](Vector2D){return Jet{double(sign),0,0,0};};
            const auto p=measuredPrimary(*mesh,f,radial);const auto pub=publish(*mesh,f,p);sourceAndCuts(*mesh,pub);return pub;}();
        check(pub.topology().curves.empty() && pub.report().graphs==0 && pub.report().quadratureNodes==0 && pub.report().maxAbsCurvature==0,"pure phase fabricated interface");
        close(pub.volumes().report().totalRawLiquid,sign<0?(radial?pi:1):0,1e-13,"pure volume");
        std::cout<<"PURE_LIFETIME,"<<(radial?"RZ":"XY")<<','<<sign<<",PASS\n";}
}
void transactions()
{
    auto mesh=makeMesh(8);const auto domain=domainFor(*mesh,false);
    const Profile f=[](Vector2D p){const auto x=p.x-.47,y=p.y-.51;return Jet{x*x+y*y-.27*.27,2*x,2*y,0};};
    auto primary=measuredPrimary(*mesh,f,false);const auto accepted=publish(*mesh,f,primary);
    const auto* oldScalar=&accepted.state().geometry().snapshot();const auto oldId=accepted.state().identity();
    HermiteGalsCandidate seed(*mesh,{},primary);const auto mixed0=closeHermiteGalsMixedAmr(*mesh,seed);const auto pred=centers(*mesh,f);
    auto patch=patchFor(*mesh,seed,mixed0,pred);const auto scalar=scalarFor(*mesh,patch,false);
    const auto initial=InterfaceCellMeasure::discoverBubbleVolumes(*mesh,domain,scalar);
    for(int c=0;c<mesh->cellIndexNum();++c)if(!initial.cellPhase[c]){const auto range=scalar->nativeQuarticCenterInterval(c,.125,.25);patch.cellCenters[c]+=.3*((c%2?range.lower:range.upper)-pred[c]);}
    const auto target=InterfaceCellMeasure::discoverVolumes(*mesh,domain,scalarFor(*mesh,patch,false));
    for(auto& m:primary.phaseVolumes)m.liquid=target.cellPhase[m.cell]?(target.cellPhase[m.cell]<0?domain.cellVolumes()[m.cell]:0):double(target.measure.negativeVolumeByCell[m.cell]);
    HermiteGalsCandidate candidate(*mesh,{1,0,.01,Position::End},primary);const auto id=candidate.identity();
    const auto mixed=closeHermiteGalsMixedAmr(*mesh,candidate);const auto bubble=closeHermiteGalsBubble(*mesh,candidate,mixed,{id,pred});
    check(bubble.report().correctedCells>0,"nonzero bubble not exercised");
    // Candidate must have copied the caller arrays before this mutation.
    const auto retained=candidate.primary();primary.vertices[0].phi+=1;primary.phaseVolumes[0].liquid+=1;
    HermiteGalsPublicationConfig cfg;cfg.relativeInventoryTolerance=1e-18;
    rejects([&]{publishHermiteGalsGeometry(*mesh,std::move(candidate),bubble,cfg);},"geometric inventory residual");
    check(candidate.identity()==id && !candidate.hasGeometry(),"late failure consumed/bound candidate");
    for(std::size_t k=0;k<retained.vertices.size();++k)check(candidate.primary().vertices[k].phi==retained.vertices[k].phi,"input alias/failure changed primary");
    for(std::size_t k=0;k<retained.phaseVolumes.size();++k)check(candidate.primary().phaseVolumes[k].liquid==retained.phaseVolumes[k].liquid,"late failure changed M");
    cfg={};cfg.unitNormalTolerance=0;rejects([&]{publishHermiteGalsGeometry(*mesh,std::move(candidate),bubble,cfg);},"invalid publication tolerance");
    const auto label=mesh->cellIndexToLabel[0];const auto center=mesh->cells(label).centerCoord();mesh->cells(label).setCenterCoord(center+Vector2D(.001,0));
    rejects([&]{publishHermiteGalsGeometry(*mesh,std::move(candidate),bubble);},"NativeMeshIdentity");
    rejects([&]{accepted.validateMesh(*mesh);},"NativeMeshIdentity");mesh->cells(label).setCenterCoord(center);
    HermiteGalsCandidate foreign(*mesh,candidate.timeLayer(),retained);
    rejects([&]{publishHermiteGalsGeometry(*mesh,std::move(foreign),bubble);},"revision mismatch");
    check(!foreign.hasGeometry(),"stale closure consumed foreign candidate");
    const auto pub=publishHermiteGalsGeometry(*mesh,std::move(candidate),bubble);sourceAndCuts(*mesh,pub);
    check(pub.state().identity()==id && pub.state().identity()!=oldId && &pub.state().geometry().snapshot()!=oldScalar,"new revision lacks new source");
    const auto obsolete=foreign.identity();const auto fm=closeHermiteGalsMixedAmr(*mesh,foreign);
    const auto fb=closeHermiteGalsBubble(*mesh,foreign,fm,{obsolete,pred});foreign.replacePrimary(*mesh,retained);
    check(foreign.identity()!=obsolete,"replacePrimary reused revision");
    rejects([&]{publishHermiteGalsGeometry(*mesh,std::move(foreign),fb);},"revision mismatch");
    const auto currentId=foreign.identity();const auto newMixed=closeHermiteGalsMixedAmr(*mesh,foreign);
    const auto newBubble=closeHermiteGalsBubble(*mesh,foreign,newMixed,{currentId,pred});
    const auto replacement=publishHermiteGalsGeometry(*mesh,std::move(foreign),newBubble);
    check(replacement.state().identity()==currentId && &replacement.state().geometry().snapshot()!=&pub.state().geometry().snapshot(),"replacement publication source");
    check(accepted.state().identity()==oldId && &accepted.state().geometry().snapshot()==oldScalar,"accepted source changed");
    sourceAndCuts(*mesh,accepted);
    std::cout<<"TRANSACTION,"<<bubble.report().correctedCells<<','<<pub.report().maxRelativeInventoryResidual<<",PASS\n";
}
}
int main(int argc,char** argv) try
{
    std::cout<<std::scientific<<std::setprecision(16);
    if(argc==2 && std::string(argv[1])=="--macro-scalar-only"){macroScalarBinding();return 0;}
    if(argc==2 && std::string(argv[1])=="--macro-only"){macroScalarBinding();macroPublicationTransport();return 0;}
    if(argc==2 && std::string(argv[1])=="--native-trace-reuse"){
        nativeSideTransactions(true);nativeCellDualMoments();std::cout<<"I11 factory trace reuse contract: PASS\n";return 0;}
    for(bool radial:{false,true})for(bool amr:{false,true})for(int sign:{-1,1})analytic(radial,amr,sign,false);
    for(bool radial:{false,true})analytic(radial,true,1,true);
    macroScalarBinding();macroPublicationTransport();pureAndLifetime();transactions();isolatedContacts();nonbinaryContacts();nativeSideTransactions();nativeCellDualMoments();nativeTraceExtended();std::cout<<"I11 authoritative geometry publication contract: PASS\n";return 0;
}
catch(const std::exception& e){std::cerr<<"I11 authoritative geometry publication contract: FAIL: "<<e.what()<<'\n';return 1;}
