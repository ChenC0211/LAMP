#include "Modules/Interface/HermiteGalsPhaseVolume.h"

#include "Modules/Interface/HermiteGalsAmrClosure.h"
#include "Core/Geometry/ActiveQuarticHermitePatch.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace
{
using namespace Quantum::Interface;
using Patch=Quantum::Geometry::ActiveQuarticHermitePatch;
using Snapshot=LevelSetHermiteGeometry::Snapshot;
using Jet=Patch::NodalDof;
using Profile=std::function<Jet(Vector2D)>;
using Key=std::array<int,2>;
constexpr long double pi=3.1415926535897932384626433832795L;
static_assert(!std::is_default_constructible_v<HermiteGalsPhaseVolumeResult>);
static_assert(std::is_same_v<decltype(std::declval<HermiteGalsPhaseVolumeResult&>().phaseVolumes()),const std::vector<HermiteGalsCellVolume>&>);
static_assert(std::is_same_v<decltype(std::declval<HermiteGalsPhaseVolumeResult&>().cells()),const std::vector<HermiteGalsPhaseVolumeCell>&>);
void check(bool ok,const std::string& why){if(!ok)throw std::runtime_error(why);}
template<class F> void rejects(F f,const std::string& reason)
{
    try{f();}catch(const std::exception& e){check(std::string(e.what()).find(reason)!=std::string::npos,
        std::string("unexpected rejection: ")+e.what());std::cout<<"REJECT,"<<reason<<'\n';return;}
    throw std::runtime_error("missing rejection: "+reason);
}
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
    HermiteGalsTimeLayer time={},bool mismatchedMeasure=false,bool extraInventory=false)
{
    HermiteGalsPrimaryState p;p.measure=radial?HermiteGalsMeasure::AxisymmetricVolume:HermiteGalsMeasure::CartesianArea;
    Patch patch;patch.nodes.resize(mesh.vertIndexNum());patch.activeCells.assign(mesh.cellIndexNum(),1);
    patch.cellCenters.resize(mesh.cellIndexNum());
    for(int v=0;v<mesh.vertIndexNum();++v)if(mesh.verts(mesh.vertIndexToLabel[v]).type()==TVertType::non){
        const auto j=f(point(mesh,v));p.vertices.push_back({v,j.value,j.dx,j.dy});patch.nodes[v]=j;}
    for(const auto& key:sides(mesh))if(!fineSide(mesh,key)){
        const auto a=point(mesh,key[0]),b=point(mesh,key[1]);const auto j=f(a+(b-a)*.5);const auto dn=a.y==b.y?j.dy:j.dx;
        p.edges.push_back({key,j.value,dn});patch.edgeMidpoints.push_back({key,j.value,dn});}
    for(int c=0;c<mesh.cellIndexNum();++c){p.geometryCells.push_back(c);p.phaseVolumes.push_back({c,0});
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
using ExactVolume=std::function<long double(Vector2D,Vector2D)>;
long double cellCapacity(Vector2D a,Vector2D b,bool radial)
{return radial?pi*(static_cast<long double>(b.x)*b.x-static_cast<long double>(a.x)*a.x)*(b.y-a.y):(static_cast<long double>(b.x)-a.x)*(b.y-a.y);}
long double verify(const std::string& name,const QuadtreeMesh& mesh,const HermiteGalsState& state,
    const ExactVolume& exact,const HermiteGalsPhaseVolumeConfig& config={},double tolerance=1e-9,std::size_t expectedComponents=0)
{
    const auto id=state.identity();const auto* scalar=&state.geometry().snapshot();const auto result=integrateHermiteGalsPhaseVolumes(mesh,state,config);
    check(result.sourceIdentity()==id && &result.geometry().snapshot()==scalar,"result changed source identity/geometry");
    check(result.timeLayer().time()==state.timeLayer().time(),"result changed source time");
    const auto& report=result.report();check(report.components==expectedComponents,"missing/extra interface component: "+name);
    check(report.totalSourceInventory==0 && report.geometricInventoryDifference==report.totalRawLiquid,"inventory discrepancy diagnostic is not raw");
    auto initialized=result.initializeCandidate(mesh);
    check(initialized.identity()!=id && !initialized.hasGeometry(),"initialization did not invalidate derived revision");
    check(initialized.primary().vertices.size()==state.primary().vertices.size() && initialized.primary().edges.size()==state.primary().edges.size(),"initialization changed master layout");
    for(std::size_t i=0;i<state.primary().vertices.size();++i){const auto a=state.primary().vertices[i],b=initialized.primary().vertices[i];
        check(a.vertex==b.vertex && a.phi==b.phi && a.dx==b.dx && a.dy==b.dy,"initialization changed first jet");}
    for(std::size_t i=0;i<state.primary().edges.size();++i){const auto a=state.primary().edges[i],b=initialized.primary().edges[i];check(a.vertices==b.vertices && a.phi==b.phi && a.dn==b.dn,"initialization changed edge");}
    long double exactSum=0,maxError=0,rawError=0;
    for(std::size_t i=0;i<result.cells().size();++i){const auto& c=result.cells()[i];const auto& cell=mesh.cells(mesh.cellIndexToLabel[c.cell]);
        const auto lo=mesh.verts(cell.cellToVert[0]).coord(),hi=mesh.verts(cell.cellToVert[2]).coord();
        const auto reference=exact(lo,hi);exactSum+=reference;
        const auto error=std::abs(result.phaseVolumes()[i].liquid-reference);maxError=std::max(maxError,error);rawError=std::max(rawError,std::abs(c.rawLiquid-reference));
        check(error<=tolerance*report.totalCapacity && std::abs(c.rawLiquid-reference)<=tolerance*report.totalCapacity,"cell volume error "+name+" cell="+std::to_string(c.cell));
        check(initialized.primary().phaseVolumes[i].liquid==result.phaseVolumes()[i].liquid,"inventory differs from measured initialization");
        check(state.primary().phaseVolumes[i].liquid==0,"source inventory was overwritten");
        if(c.certifiedPhase!=0)check(result.phaseVolumes()[i].liquid==(c.certifiedPhase<0?c.capacity:0),"pure inventory not exact");
        std::cout<<"CELL,"<<name<<','<<c.cell<<','<<c.certifiedPhase<<','<<c.capacity<<','<<c.rawLiquid<<','<<result.phaseVolumes()[i].liquid<<','<<reference<<','<<error<<'\n';
    }
    const auto totalError=std::abs(report.totalInitialLiquid-exactSum);
    check(totalError<=tolerance*report.totalCapacity,"global volume error "+name);
    check(state.identity()==id && &state.geometry().snapshot()==scalar,"measurement mutated source");
    std::cout<<"VOLUME,"<<name<<','<<result.cells().size()<<','<<report.components<<','<<report.cutCells<<','<<report.totalInitialLiquid<<','<<exactSum
        <<','<<totalError<<','<<maxError<<','<<rawError<<','<<report.relativeMomentEstimate<<','<<report.relativeBoundaryClosure<<'\n';
    return totalError;
}
void analyticCases()
{
    for(bool radial:{false,true})for(bool amr:{false,true}){
        auto mesh=meshFor(8,amr);const std::string tag=std::string(radial?"RZ":"XY")+(amr?"_AMR":"_uniform");
        for(int sign:{-1,1}){const auto state=initialTemplate(*mesh,[&](Vector2D){return Jet{static_cast<double>(sign),0,0,0};},radial);
            verify(tag+(sign<0?"_liquid":"_gas"),*mesh,state,[&](Vector2D a,Vector2D b){return sign<0?cellCapacity(a,b,radial):0.L;},{},3e-12);}
        const double z=.421875;const auto planar=initialTemplate(*mesh,[&](Vector2D p){return Jet{p.y-z,0,1,0};},radial);
        HermiteGalsPhaseVolumeConfig boundary;boundary.symmetryPlanes={{0,1}};if(!radial)boundary.symmetryPlanes.push_back({0,0});
        verify(tag+"_plane",*mesh,planar,[&](Vector2D a,Vector2D b){b.y=std::clamp(z,a.y,b.y);return cellCapacity(a,b,radial);},boundary,3e-12,1);
        const Circle circle{radial?0:.47,.53,.231};const auto state=initialTemplate(*mesh,[&](Vector2D p){return circleJet(p,circle);},radial);
        verify(tag+"_disk",*mesh,state,[&](Vector2D a,Vector2D b){return diskVolume(a,b,circle,radial);},{},1e-9,1);
        const auto exterior=initialTemplate(*mesh,[&](Vector2D p){auto j=circleJet(p,circle);return Jet{-j.value,-j.dx,-j.dy,-j.dxy};},radial);
        verify(tag+"_exterior",*mesh,exterior,[&](Vector2D a,Vector2D b){return cellCapacity(a,b,radial)-diskVolume(a,b,circle,radial);},{},1e-9,1);
        const Circle first{radial?0:.32,.31,.145},second{radial?0:.68,.72,.12};
        const auto multiple=initialTemplate(*mesh,[&](Vector2D p){return product(circleJet(p,first),circleJet(p,second));},radial);
        verify(tag+"_multiple",*mesh,multiple,[&](Vector2D a,Vector2D b){return diskVolume(a,b,first,radial)+diskVolume(a,b,second,radial);},{},1e-9,2);
    }
    auto torusMesh=meshFor(8,true);Circle torus{.55,.5,.19},mirror{-.55,.5,.19};
    // Even quartic torus scalar on the full r>=0 domain. The second factor
    // is strictly positive there, so the negative region is the one disk.
    const auto state=initialTemplate(*torusMesh,[&](Vector2D p){return product(circleJet(p,torus),circleJet(p,mirror));},true);
    verify("RZ_torus",*torusMesh,state,[&](Vector2D a,Vector2D b){return diskVolume(a,b,torus,true);},{},1e-9,1);
    for(bool radial:{false,true}){
        DictionaryMesh d;d.Nx=d.Ny=8;d.maxCellNum=2048;d.maxInitialIterations=0;
        d.minCoord={radial?0:.125,-.25};d.maxCoord={1.625,1.25};QuadtreeMesh shifted(d);
        DictionaryMesh::RefineBox b;b.level=1;b.min={d.minCoord.x,.125};b.max={.94,.688};shifted.applyRefineBoxes({b});
        const Circle c{radial?0:.71,.37,.19};
        const auto shape=initialTemplate(shifted,[&](Vector2D p){return circleJet(p,c);},radial);
        verify(radial?"RZ_shifted_anisotropic":"XY_shifted_anisotropic",shifted,shape,
            [&](Vector2D a,Vector2D z){return diskVolume(a,z,c,radial);},{},1e-9,1);
    }
    // A one-graph corner cap is below the existing complete-topology contract:
    // retain its explicit resolution rejection, never discard the component.
    auto mesh=meshFor(8);const double cut=1.0/65536;
    const auto small=initialTemplate(*mesh,[&](Vector2D p){return Jet{p.x+p.y-cut,1,1,0};},false);
    HermiteGalsPhaseVolumeConfig boundary;boundary.boundaryContacts={{cut,0},{0,cut}};
    rejects([&]{integrateHermiteGalsPhaseVolumes(*mesh,small,boundary);},"fewer than three retained graphs");
    // A resolved eight-graph film still has arbitrarily small cell fractions.
    // Its area must be retained with relative accuracy, not just the global
    // 1e-9 absolute gate that would otherwise allow deletion.
    const double filmHeight=std::ldexp(1.,-30);
    const auto film=initialTemplate(*mesh,[&](Vector2D p){return Jet{p.y-filmHeight,0,1,0};},false);
    boundary={};boundary.symmetryPlanes={{0,0},{0,1}};
    verify("XY_small_film",*mesh,film,[&](Vector2D a,Vector2D b){return a.y==0?(b.x-a.x)*filmHeight:0.L;},boundary,3e-12,1);
    const auto result=integrateHermiteGalsPhaseVolumes(*mesh,film,boundary);
    check(result.report().totalInitialLiquid>0 && std::abs(result.report().totalInitialLiquid/filmHeight-1)<1e-10,"small film was lost");
}
void refinementAndConvergence()
{
    for(bool radial:{false,true}){
        const Circle circle{radial?0:.47,.53,.231};const Profile f=[&](Vector2D p){return circleJet(p,circle);};
        auto coarse=meshFor(8),fine=meshFor(8,true);
        auto c=integrateHermiteGalsPhaseVolumes(*coarse,initialTemplate(*coarse,f,radial));
        auto r=integrateHermiteGalsPhaseVolumes(*fine,initialTemplate(*fine,f,radial));
        long double maximum=0;
        for(const auto& parent:c.cells()){const auto& pc=coarse->cells(coarse->cellIndexToLabel[parent.cell]);
            const auto lo=coarse->verts(pc.cellToVert[0]).coord(),hi=coarse->verts(pc.cellToVert[2]).coord();long double sum=0;
            for(const auto& child:r.phaseVolumes()){const auto& cc=fine->cells(fine->cellIndexToLabel[child.cell]);const auto x=cc.centerCoord();
                if(x.x>lo.x && x.x<hi.x && x.y>lo.y && x.y<hi.y)sum+=child.liquid;}
            maximum=std::max(maximum,std::abs(sum-c.phaseVolumes()[parent.cell].liquid));}
        check(maximum<=1e-9*c.report().totalCapacity,"coarse/fine phase volume sum mismatch");
        std::cout<<"REFINEMENT,"<<(radial?"RZ":"XY")<<','<<c.cells().size()<<','<<r.cells().size()<<','<<maximum<<'\n';
        long double previous=0,lastRate=0;
        for(int n:{8,16,32}){auto mesh=meshFor(n);const auto state=initialTemplate(*mesh,[&](Vector2D p){
                const auto q=circleJet(p,circle);const auto e=std::exp(q.value);return Jet{e-1,e*q.dx,e*q.dy,e*q.dx*q.dy};},radial);
            const auto result=integrateHermiteGalsPhaseVolumes(*mesh,state);const auto reference=diskVolume(mesh->minCoord,mesh->maxCoord,circle,radial);
            const auto error=std::abs(result.report().totalInitialLiquid-reference);const auto rate=previous>0?std::log2(previous/error):0;
            if(previous>0)check(error<previous,"smooth geometric volume error not decreasing");
            std::cout<<"CONVERGENCE,"<<(radial?"RZ":"XY")<<','<<n<<','<<error<<','<<rate<<','<<result.report().relativeMomentEstimate<<'\n';
            previous=error;lastRate=rate;}
        check(lastRate>=3,"smooth geometric volume last order below 3");
    }
}
void lifecycleAndFailures()
{
    auto mesh=meshFor(8);const Circle circle{.47,.53,.231};const Profile f=[&](Vector2D p){return circleJet(p,circle);};
    const auto source=initialTemplate(*mesh,f,false);const auto id=source.identity();
    HermiteGalsPhaseVolumeConfig config;config.relativeLedgerTolerance=0;
    rejects([&]{integrateHermiteGalsPhaseVolumes(*mesh,source,config);},"relative tolerances");
    config={};config.relativeMomentEstimateTolerance=std::numeric_limits<double>::quiet_NaN();
    rejects([&]{integrateHermiteGalsPhaseVolumes(*mesh,source,config);},"relative tolerances");
    auto changed=meshFor(8,true);rejects([&]{integrateHermiteGalsPhaseVolumes(*changed,source);},"mesh");
    const auto mismatch=initialTemplate(*mesh,f,false,{},true);
    rejects([&]{integrateHermiteGalsPhaseVolumes(*mesh,mismatch);},"measure mismatch");
    const auto incomplete=initialTemplate(*mesh,[](Vector2D){return Jet{1,0,0,0};},false,{},false,true);
    rejects([&]{integrateHermiteGalsPhaseVolumes(*mesh,incomplete);},"scalar support");
    for(const auto time:{HermiteGalsTimeLayer{1,0,.1,HermiteGalsTimeLayer::Position::Midpoint},
        HermiteGalsTimeLayer{1,0,.1,HermiteGalsTimeLayer::Position::End},
        HermiteGalsTimeLayer{4,1,1,HermiteGalsTimeLayer::Position::Initial}}){
        const auto evolving=initialTemplate(*mesh,f,false,time);const auto result=integrateHermiteGalsPhaseVolumes(*mesh,evolving);
        check(result.timeLayer().position==time.position && result.timeLayer().time()==time.time(),"measurement lost time layer");
        rejects([&]{result.initializeCandidate(*mesh);},"restart/advancement forbidden");
        for(const auto& volume:evolving.primary().phaseVolumes)check(volume.liquid==0,"measurement replaced evolving M");}
    const auto ambiguous=initialTemplate(*mesh,[](Vector2D){return Jet{0,0,0,0};},false);
    rejects([&]{integrateHermiteGalsPhaseVolumes(*mesh,ambiguous);},"HermiteGalsPhaseVolume source=");
    auto kept=[&]{return integrateHermiteGalsPhaseVolumes(*mesh,initialTemplate(*mesh,f,false));}();
    auto copy=kept;const auto oldId=copy.sourceIdentity();const auto* oldGeometry=&copy.geometry().snapshot();
    auto candidate=copy.initializeCandidate(*mesh);auto transported=candidate.primary();transported.phaseVolumes.front().liquid=.01;
    candidate.replacePrimary(*mesh,transported);
    check(copy.sourceIdentity()==oldId && &copy.geometry().snapshot()==oldGeometry,"volume result lifetime lost authority");
    check(candidate.primary().phaseVolumes.front().liquid==.01,"M failed independent persistence");
    rejects([&]{copy.initializeCandidate(*changed);},"mesh");
    check(source.identity()==id,"failure changed source identity");
    std::cout<<"LIFECYCLE,retained_source,new_revision,independent_inventory,PASS\n";
}
}
int main() try
{
    std::cout<<std::scientific<<std::setprecision(16);
    analyticCases();refinementAndConvergence();lifecycleAndFailures();
    std::cout<<"I07 phase-volume contract: PASS\n";return 0;
}catch(const std::exception& e){std::cerr<<"I07 phase-volume contract: FAIL: "<<e.what()<<'\n';return 1;}
