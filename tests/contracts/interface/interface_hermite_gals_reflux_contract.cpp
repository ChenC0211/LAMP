#include "Modules/Interface/HermiteGalsReflux.h"

#include "Modules/Interface/HermiteGalsAmrClosure.h"
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
using Coeff=std::array<double,11>;
using Boundary=HermiteGalsFaceVelocity::Boundary;
using Position=HermiteGalsTimeLayer::Position;
constexpr long double pi=3.1415926535897932384626433832795L;
static_assert(!std::is_default_constructible_v<HermiteGalsPhaseFluxResult>);
static_assert(std::is_same_v<decltype(std::declval<HermiteGalsPhaseFluxResult&>().faces()),const std::vector<HermiteGalsPhaseFaceFlux>&>);
static_assert(std::is_same_v<decltype(std::declval<HermiteGalsFluxVelocity&>().faces()),const std::vector<HermiteGalsFaceVelocity>&>);
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
using Velocity=std::function<Coeff(const NativeLeafFace&,int,Vector2D,Vector2D)>;
std::vector<HermiteGalsFaceVelocity> velocityFor(const QuadtreeMesh& mesh,bool radial,const Velocity& field,bool closed=false)
{
    const auto domain=domainFor(mesh,radial);std::vector<HermiteGalsFaceVelocity> result;
    for(const auto face:domain.faces(mesh)){HermiteGalsFaceVelocity f;f.nativeFace=face.nativeRow;
        const bool internal=face.negativeCell>=0 && face.positiveCell>=0;
        f.boundary=internal?Boundary::Interior:((closed || (radial && !face.horizontal && face.center.x==0))?Boundary::Impermeable:Boundary::Open);
        const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
        auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();if(face.horizontal?a.x>b.x:a.y>b.y)std::swap(a,b);
        for(int half=0;half<2;++half)f.halfPolynomial[half]=field(face,half,half==0?a:face.center,half==0?face.center:b);
        result.push_back(f);}return result;
}
Velocity constantVelocity(double x,double y)
{return [=](const NativeLeafFace& face,int,Vector2D,Vector2D){Coeff p{};p[0]=face.horizontal?y:x;return p;};}
static_assert(!std::is_default_constructible_v<HermiteGalsRefluxResult>);
static_assert(std::is_same_v<decltype(std::declval<HermiteGalsRefluxResult&>().source()),const HermiteGalsPhaseFluxResult&>);
static_assert(std::is_same_v<decltype(std::declval<HermiteGalsRefluxResult&>().sides()),const std::vector<HermiteGalsSideFlux>&>);
void verify(const QuadtreeMesh& mesh,const HermiteGalsState& old,
    const HermiteGalsPhaseFluxResult& flux,const HermiteGalsRefluxResult& reflux,int expectedAmr)
{
    reflux.validateMesh(mesh);
    check(&reflux.source().faces()==&flux.faces() && &reflux.source().phaseVolumes()==&flux.phaseVolumes(),"duplicate flux or inventory owner");
    check(&reflux.source().geometry().snapshot()==&flux.geometry().snapshot()
        && reflux.source().geometryIdentity()==flux.geometryIdentity()
        && reflux.source().previousIdentity()==old.identity(),"identity/geometry changed");
    check(reflux.source().timeLayer().step==flux.timeLayer().step && reflux.source().timeLayer().time()==flux.timeLayer().time(),"time detached");
    check(reflux.sides().size()==4*flux.phaseVolumes().size(),"missing complete sides");
    std::vector<int> incidence(flux.faces().size(),0);
    std::vector<long double> change(mesh.cellIndexNum(),0);
    int directions=0;std::size_t coarse=0;
    for(std::size_t k=0;k<reflux.sides().size();++k){const auto& side=reflux.sides()[k];const int s=static_cast<int>(side.side);
        check(side.cell==flux.phaseVolumes()[k/4].cell && s==int(k%4),"side order");
        check(side.outwardSign==(s<2?-1:1),"outward sign");
        // Independent oracle: scan actual face owners and direction, without
        // using the production FaceSegmentGraph grouping or side metadata.
        std::vector<std::size_t> expected;
        for(std::size_t j=0;j<flux.faces().size();++j){const auto& f=flux.faces()[j];const bool horizontal=f.nativeFace>=mesh.vEdgeIndexNum();
            if(horizontal!=(s==1 || s==3))continue;
            if((s<2?f.positiveCell:f.negativeCell)==side.cell)expected.push_back(j);}
        const auto coord=[&](std::size_t j){const int face=flux.faces()[j].nativeFace;return face>=mesh.vEdgeIndexNum()?
            mesh.hEdges(mesh.hEdgeIndexToLabel[face-mesh.vEdgeIndexNum()]).midpointCoord().x:
            mesh.vEdges(mesh.vEdgeIndexToLabel[face]).midpointCoord().y;};
        std::sort(expected.begin(),expected.end(),[&](auto a,auto b){return coord(a)<coord(b);});
        check(expected.size()==side.subfaceCount && (expected.size()==1 || expected.size()==2),"subface membership");
        std::array<long double,8> sums{};
        for(std::size_t j=0;j<expected.size();++j){check(expected[j]==side.subfaces[j],"subface identity/order");
            ++incidence[expected[j]];const auto& f=flux.faces()[expected[j]];
            sums[0]+=f.totalPositiveRate;sums[1]+=f.totalNegativeRate;sums[2]+=f.liquidPositiveRate;sums[3]+=f.liquidNegativeRate;
            sums[4]+=f.totalRate;sums[5]+=f.liquidRate;sums[6]+=f.liquidVolume;sums[7]+=f.integrationErrorEstimate;}
        const std::array<long double,8> actual{{side.totalPositiveRate,side.totalNegativeRate,side.liquidPositiveRate,side.liquidNegativeRate,
            side.totalRate,side.liquidRate,side.liquidVolume,side.integrationErrorEstimate}};
        check(actual==sums,"coarse flux differs from final fine subflux sum");
        if(expected.size()==2){++coarse;directions|=1<<s;
            const auto& c=mesh.cells(mesh.cellIndexToLabel[side.cell]);
            for(auto j:expected){const auto& f=flux.faces()[j];const auto n=s<2?f.negativeCell:f.positiveCell;
                if(n>=0)check(mesh.cells(mesh.cellIndexToLabel[n]).level()==c.level()+1,"coarse/fine levels");}}
        change[side.cell]-=side.outwardSign*side.liquidVolume;
    }
    for(std::size_t j=0;j<flux.faces().size();++j){const auto& f=flux.faces()[j];
        check(incidence[j]==int(f.negativeCell>=0)+int(f.positiveCell>=0),"orphan/duplicate shared subface");}
    const auto domain=domainFor(mesh,old.primary().measure==HermiteGalsMeasure::AxisymmetricVolume);
    for(std::size_t k=0;k<flux.phaseVolumes().size();++k){const auto m=flux.phaseVolumes()[k];
        close(m.liquid,old.primary().phaseVolumes[k].liquid+change[m.cell],2e-15*domain.cellVolumes()[m.cell],"side-based M ledger");
        check(m.liquid>=0 && m.liquid<=domain.cellVolumes()[m.cell],"M unbounded");}
    if(expectedAmr>=0)check(expectedAmr?(coarse>0 && directions==15):coarse==0,"AMR coverage/all four orientations");
    check(reflux.report().coarseSides==coarse,"coarse diagnostic count");
    close(reflux.report().maxTotalRateResidual,0,1e-13,"total reflux residual");
    close(reflux.report().maxLiquidRateResidual,0,1e-13,"liquid reflux residual");
    close(reflux.report().maxLiquidVolumeResidual,0,1e-13,"volume reflux residual");
    close(flux.report().closure,0,1e-12,"global volume closure");
}

void analyticSides()
{
    // Avoid the separately documented unsupported positive-length zero trace.
    const double lo=.277,hi=.349,turn=.3125;
    for(bool radial:{false,true})for(int level:{0,1,2})for(int orientation:{0,1})for(int sign:{-1,1}){
        auto mesh=makeMesh(8,level);const auto domain=domainFor(*mesh,radial);
        const Profile band=[=](Vector2D p){if(orientation)return Jet{(p.y-lo)*(p.y-hi),0,2*p.y-lo-hi,0};
            if(radial)return Jet{(p.x*p.x-lo*lo)*(p.x*p.x-hi*hi),2*p.x*(2*p.x*p.x-lo*lo-hi*hi),0,0};
            return Jet{(p.x-lo)*(p.x-hi),2*p.x-lo-hi,0,0};};
        const auto half=[&](int c,Vector2D,Vector2D){return .5*domain.cellVolumes()[c];};
        const auto old=stateFor(*mesh,band,radial,{},half),mid=stateFor(*mesh,band,radial,{1,0,.001,Position::Midpoint},half);
        const Velocity field=[=](const NativeLeafFace& f,int,Vector2D a,Vector2D b){Coeff c{};
            if(orientation){c[0]=sign*(radial?(f.horizontal?-2*f.center.y:f.center.x):(f.horizontal?0:1));}
            else if(f.horizontal){const auto h=b.x-a.x;c[0]=sign*(radial?a.x*a.x-turn*turn:a.x-turn);
                c[1]=sign*(radial?2*a.x*h:h);if(radial)c[2]=sign*h*h;}return c;};
        const HermiteGalsFluxVelocity velocity(*mesh,mid,velocityFor(*mesh,radial,field));
        const auto flux=advanceHermiteGalsPhaseVolumes(*mesh,old,mid,velocity);const auto result=assembleHermiteGalsReflux(*mesh,flux);
        verify(*mesh,old,flux,result,level);
        long double error=0;int counter=0,cutCoarse=0;
        for(const auto& side:result.sides()){
            const auto& cell=mesh->cells(mesh->cellIndexToLabel[side.cell]);const auto a=mesh->verts(cell.cellToVert[0]).coord(),b=mesh->verts(cell.cellToVert[2]).coord();
            const int s=static_cast<int>(side.side);const bool horizontal=s==1 || s==3;
            const double fixed=horizontal?(s==1?a.y:b.y):(s==0?a.x:b.x);
            const long double x0=orientation?a.x:std::max(a.x,lo),x1=orientation?b.x:std::min(b.x,hi);
            const long double y0=orientation?std::max(a.y,lo):a.y,y1=orientation?std::min(b.y,hi):b.y;
            long double expected=0;
            if(!orientation && horizontal && x1>x0)expected=sign*(radial?
                pi*((std::pow(x1,4)-std::pow(x0,4))/2-turn*turn*(x1*x1-x0*x0)):
                (x1*x1-x0*x0)/2-turn*(x1-x0));
            if(orientation){if(horizontal && fixed>lo && fixed<hi && radial)expected=-2*sign*fixed*pi*(x1*x1-x0*x0);
                if(!horizontal && y1>y0)expected=sign*(radial?2*pi*fixed*fixed:1)*(y1-y0);}
            error=std::max(error,close(side.liquidRate,expected,3e-11,"analytic complete-side flux"));
            counter+=side.liquidPositiveRate>0 && side.liquidNegativeRate>0;
            const auto liquidGross=side.liquidPositiveRate+side.liquidNegativeRate;
            const auto totalGross=side.totalPositiveRate+side.totalNegativeRate;
            cutCoarse+=side.subfaceCount==2 && liquidGross>0 && liquidGross<totalGross;
        }
        if(!orientation)check(counter>0,"missing liquid counterflow");
        if(level)check(cutCoarse>0,"missing cut coarse face");
        std::cout<<"ANALYTIC,"<<(radial?"RZ":"XY")<<','<<level<<','<<orientation<<','<<sign<<','<<mesh->cellIndexNum()<<','
            <<result.report().coarseSides<<','<<counter<<','<<cutCoarse<<','<<error<<','<<result.report().maxLiquidRateResidual<<','<<flux.report().closure<<'\n';
    }
}
void pureAndStationary()
{
    for(bool radial:{false,true})for(int phase:{-1,1}){
        auto mesh=makeMesh(8,2);const auto domain=domainFor(*mesh,radial);const Profile scalar=[=](Vector2D){return Jet{double(phase),0,0,0};};
        const auto inventory=[&](int c,Vector2D,Vector2D){return phase<0?domain.cellVolumes()[c]:0.;};
        const auto old=stateFor(*mesh,scalar,radial,{},inventory),mid=stateFor(*mesh,scalar,radial,{1,0,.01,Position::Midpoint},inventory);
        for(double speed:{-.25,0.,.25}){const HermiteGalsFluxVelocity velocity(*mesh,mid,velocityFor(*mesh,radial,constantVelocity(0,speed),speed==0));
            const auto flux=advanceHermiteGalsPhaseVolumes(*mesh,old,mid,velocity);const auto result=assembleHermiteGalsReflux(*mesh,flux);verify(*mesh,old,flux,result,2);
            for(const auto& s:result.sides())check(s.liquidRate==(phase<0?s.totalRate:0),"pure-phase side identity");
            for(std::size_t k=0;k<flux.phaseVolumes().size();++k)check(flux.phaseVolumes()[k].liquid==old.primary().phaseVolumes[k].liquid,"pure/zero inventory changed");
            std::cout<<"PURE,"<<(radial?"RZ":"XY")<<','<<phase<<','<<speed<<','<<result.report().coarseSides<<'\n';}
    }
}
void translation()
{
    for(bool radial:{false,true})for(int sign:{-1,1}){
        // Cross face planes at step endpoints: the prescribed midpoint rule
        // then integrates each constant-in-time face phase exactly. Off-step
        // crossings require a temporal convergence study in I12, not an
        // assertion that a Heaviside discontinuity has exact midpoint integral.
        auto mesh=makeMesh(8,2);const double dt=.0625,u=.25*sign,start=sign>0?.21875:.34375;
        const auto scalar=[=](double height){return Profile([=](Vector2D p){return radial?Jet{p.y-height,0,1,0}:Jet{p.x-height,1,0,0};});};
        const auto amount=[=](double height,Vector2D a,Vector2D b){return radial?
            pi*(b.x*b.x-a.x*a.x)*std::clamp(height-a.y,0.,b.y-a.y):std::clamp(height-a.x,0.,b.x-a.x)*(b.y-a.y);};
        auto old=stateFor(*mesh,scalar(start),radial,{},[&](int,Vector2D a,Vector2D b){return double(amount(start,a,b));});
        long double maxError=0,maxClosure=0;
        for(int step=1;step<=8;++step){const double begin=(step-1)*dt,end=step*dt,height=start+u*end;
            const auto mid=stateFor(*mesh,scalar(start+u*(begin+end)/2),radial,{step,begin,end,Position::Midpoint},[](int,Vector2D,Vector2D){return 0.;});
            const HermiteGalsFluxVelocity velocity(*mesh,mid,velocityFor(*mesh,radial,constantVelocity(radial?0:u,radial?u:0)));
            const auto flux=advanceHermiteGalsPhaseVolumes(*mesh,old,mid,velocity);const auto result=assembleHermiteGalsReflux(*mesh,flux);verify(*mesh,old,flux,result,2);
            for(const auto m:flux.phaseVolumes()){const auto& c=mesh->cells(mesh->cellIndexToLabel[m.cell]);const auto a=mesh->verts(c.cellToVert[0]).coord(),b=mesh->verts(c.cellToVert[2]).coord();
                maxError=std::max(maxError,close(m.liquid,amount(height,a,b),3e-11,"multistep AMR translation"));}
            maxClosure=std::max(maxClosure,std::abs(flux.report().closure));
            // Carry the conserved result into the next operator step. Analytic
            // endpoint/midpoint geometry is supplied, not an I11 publisher.
            old=stateFor(*mesh,scalar(height),radial,{step,begin,end,Position::End},[&](int c,Vector2D,Vector2D){return flux.phaseVolumes()[c].liquid;});
        }
        std::cout<<"AMR_TRANSLATION,"<<(radial?"RZ":"XY")<<','<<sign<<",8,"<<maxError<<','<<maxClosure<<'\n';
    }
}
void repairedCoarseFace()
{
    auto mesh=makeMesh(8,1);const auto domain=domainFor(*mesh,false);const double dt=.125,u=.125,split=.25;
    const Profile scalar=[=](Vector2D p){return Jet{split-u*dt/2-p.x,-1,0,0};};
    const auto inventory=[&](int c,Vector2D a,Vector2D b){return a.x>=split?domain.cellVolumes()[c]:
        (b.x==split?std::nextafter(dt*u*(b.y-a.y),0.):0.);};
    const auto old=stateFor(*mesh,scalar,false,{},inventory),mid=stateFor(*mesh,scalar,false,{1,0,dt,Position::Midpoint},inventory);
    const HermiteGalsFluxVelocity velocity(*mesh,mid,velocityFor(*mesh,false,constantVelocity(u,0)));
    const auto flux=advanceHermiteGalsPhaseVolumes(*mesh,old,mid,velocity);const auto result=assembleHermiteGalsReflux(*mesh,flux);verify(*mesh,old,flux,result,1);
    std::size_t repairedCoarse=0;long double maxDifference=0;
    for(const auto& s:result.sides())if(s.subfaceCount==2){long double raw=0;
        for(auto row:s.subfaces){const auto& f=flux.faces()[row];raw+=f.rawLiquidPositiveRate-f.rawLiquidNegativeRate;}
        if(raw!=s.liquidRate){++repairedCoarse;maxDifference=std::max(maxDifference,std::abs(raw-s.liquidRate));}}
    check(repairedCoarse>0 && flux.report().repairedFaces>0,"no repaired coarse/fine face exercised");
    std::cout<<"REPAIRED_COARSE,"<<repairedCoarse<<','<<flux.report().repairedFaces<<','<<maxDifference<<','<<flux.report().closure<<'\n';
}
void lifetimeAndFailures()
{
    auto mesh=makeMesh(8,1);std::optional<HermiteGalsRefluxResult> retained;
    const LevelSetHermiteGeometry::Snapshot* geometry=nullptr;const void* faces=nullptr;HermiteGalsIdentity id;
    {
        const auto domain=domainFor(*mesh,false);const Profile gas=[](Vector2D){return Jet{1,0,0,0};};
        const auto inventory=[](int,Vector2D,Vector2D){return 0.;};
        const auto old=stateFor(*mesh,gas,false,{},inventory),mid=stateFor(*mesh,gas,false,{1,0,.01,Position::Midpoint},inventory);
        const HermiteGalsFluxVelocity zero(*mesh,mid,velocityFor(*mesh,false,constantVelocity(0,0),true));
        auto flux=advanceHermiteGalsPhaseVolumes(*mesh,old,mid,zero);faces=&flux.faces();geometry=&mid.geometry().snapshot();id=mid.identity();
        auto moved=std::move(flux);rejects([&]{assembleHermiteGalsReflux(*mesh,flux);},"moved-from");
        retained=assembleHermiteGalsReflux(*mesh,moved);auto copy=*retained;
        auto movedResult=std::move(copy);rejects([&]{copy.sides();},"moved-from");
        check(&movedResult.source().faces()==faces,"result move changed owner");
    }
    check(retained->source().geometryIdentity()==id && &retained->source().geometry().snapshot()==geometry
        && &retained->source().faces()==faces,"retained result lifetime");
    const auto label=mesh->cellIndexToLabel[0];const auto center=mesh->cells(label).centerCoord();
    mesh->cells(label).setCenterCoord(center+Vector2D(.001,0));
    rejects([&]{retained->validateMesh(*mesh);},"NativeMeshIdentity");
    rejects([&]{assembleHermiteGalsReflux(*mesh,retained->source());},"NativeMeshIdentity");
    mesh->cells(label).setCenterCoord(center);retained->validateMesh(*mesh);
    auto foreign=makeMesh(8);rejects([&]{assembleHermiteGalsReflux(*foreign,retained->source());},"NativeMeshIdentity");
    std::cout<<"OWNERSHIP,PASS\n";
}
}
int main()
{
    try{std::cout<<std::setprecision(17);analyticSides();pureAndStationary();translation();repairedCoarseFace();lifetimeAndFailures();
        std::cout<<"Hermite GALS 2:1 reflux contract PASS\n";return 0;}
    catch(const std::exception& e){std::cerr<<"Hermite GALS reflux FAIL: "<<e.what()<<'\n';return 1;}
}
