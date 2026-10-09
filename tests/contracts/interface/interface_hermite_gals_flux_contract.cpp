#include "Modules/Interface/HermiteGalsPhaseFlux.h"

#include "Modules/Interface/HermiteGalsAmrClosure.h"
#include "Modules/Interface/HermiteGalsPhaseVolume.h"
#include "Modules/Interface/InterfacePhaseGeometry.h"
#include "Core/Geometry/ActiveQuarticHermitePatch.h"
#include "Core/Geometry/UnitPolynomialRoots.h"
#include "Core/Geometry/LevelSetMappedSegment.h"
#include "Core/Geometry/LevelSetInterfacePanels.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <cmath>
#include <cfenv>
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
using Term=Quantum::Geometry::UnitPolynomialTerm;
constexpr long double pi=3.1415926535897932384626433832795L;
static_assert(!std::is_default_constructible_v<HermiteGalsPhaseFluxResult>);
static_assert(std::is_same_v<decltype(std::declval<HermiteGalsPhaseFluxResult&>().faces()),const std::vector<HermiteGalsPhaseFaceFlux>&>);
static_assert(std::is_same_v<decltype(std::declval<HermiteGalsFluxVelocity&>().faces()),const std::vector<HermiteGalsFaceVelocity>&>);
void check(bool ok,const std::string& why){if(!ok)throw std::runtime_error(why);}
long double close(long double a,long double b,long double tolerance,const std::string& why)
{const auto error=std::abs(a-b);check(std::isfinite(a) && error<=tolerance,why+" error="+std::to_string(static_cast<double>(error)));return error;}
template<class F> void rejects(F f,const std::string& expected)
{try{f();}catch(const std::exception& e){check(std::string(e.what()).find(expected)!=std::string::npos,std::string("unexpected rejection: ")+e.what());std::cout<<"REJECT,"<<expected<<'\n';return;}throw std::runtime_error("missing rejection: "+expected);}
std::unique_ptr<QuadtreeMesh> makeMesh(int n,bool amr=false)
{DictionaryMesh d;d.Nx=d.Ny=n;d.maxCellNum=16384;d.maxInitialIterations=0;auto mesh=std::make_unique<QuadtreeMesh>(d);
    if(amr){DictionaryMesh::RefineBox b;b.level=1;b.min={.25,.25};b.max={.75,.75};mesh->applyRefineBoxes({b});}return mesh;}
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
long double capacity(Vector2D a,Vector2D b,bool radial)
{return radial?pi*(b.x*b.x-a.x*a.x)*(b.y-a.y):(b.x-a.x)*(b.y-a.y);}
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
void coreRanges()
{
    namespace G=Quantum::Geometry;
    Term constant;constant.factor=1;constant.coefficients[0]=3;
    const long double tiny=std::ldexp(1.L,-60);
    for(bool reverse:{false,true}){
        const auto value=G::integratePositiveUnitPolynomialRange({constant},0,1,reverse,.5L,tiny,1e-14L*tiny);
        close(value.value,3*tiny,1e-13L*tiny,"sub-ulp constant range");
        Term slope;slope.factor=1;slope.coefficients[0]=reverse?.5:-.5;slope.coefficients[1]=reverse?-1:1;
        const auto root=G::integratePositiveUnitPolynomialRange({slope},0,1,reverse,.5L,tiny,1e-13L*tiny*tiny);
        close(root.value,tiny*tiny/2,1e-12L*tiny*tiny,"sub-ulp linear root");
        close(G::integrateUnitPolynomialRange({slope},0,1,reverse,.5L,tiny),tiny*tiny/2,1e-12L*tiny*tiny,"old signed range regression");
        std::cout<<"CORE_TINY,"<<reverse<<','<<value.value<<','<<root.value<<','<<root.absoluteError<<'\n';}
    Term degree;degree.factor=1;degree.coefficients[10]=1;
    const auto p=G::integratePositiveUnitPolynomialRange({degree},0,1,false,.25,.5,1e-13);
    close(p.value,(std::pow(.75L,11)-std::pow(.25L,11))/11,3e-12,"degree-10 exact integral");
    Term crossing;crossing.factor=1;crossing.coefficients[0]=-.375;crossing.coefficients[1]=1;
    const auto newValue=G::integratePositiveUnitPolynomialRange({crossing},0,1,false,.125,.75,1e-13);
    const auto oldValue=G::integratePositiveUnitPolynomial({crossing},.125,.875,1e-13);
    close(newValue.value,.125L,3e-12,"positive crossing analytic");close(newValue.value,oldValue.value,3e-12,"old positive primitive");
    rejects([&]{G::integratePositiveUnitPolynomialRange({constant},0,1,false,-.1,.2,1e-13);},"invalid mapped polynomial range");
    std::cout<<"CORE_RANGE,PASS\n";
}
void verifyLedger(const QuadtreeMesh& mesh,const HermiteGalsState& previous,const HermiteGalsState& mid,const HermiteGalsPhaseFluxResult& result,bool expectRepair=false)
{
    const auto domain=domainFor(mesh,mid.primary().measure==HermiteGalsMeasure::AxisymmetricVolume);
    check(result.previousIdentity()==previous.identity() && result.geometryIdentity()==mid.identity() && &result.geometry().snapshot()==&mid.geometry().snapshot(),"flux result lost authority");
    std::vector<long double> change(mesh.cellIndexNum(),0);long double boundary=0,cap=0;int last=-1;
    for(const auto& f:result.faces()){
        check(f.nativeFace>last,"duplicate face result");last=f.nativeFace;
        check(f.liquidPositiveRate>=0 && f.liquidPositiveRate<=f.totalPositiveRate && f.liquidNegativeRate>=0 && f.liquidNegativeRate<=f.totalNegativeRate,"gross phase bounds");
        check(f.liquidVolume==(mid.timeLayer().end-mid.timeLayer().begin)*f.liquidRate,"cached flux volume differs");
        if(f.negativeCell>=0)change[f.negativeCell]-=f.liquidVolume;else boundary-=f.liquidVolume;
        if(f.positiveCell>=0)change[f.positiveCell]+=f.liquidVolume;else boundary+=f.liquidVolume;
        if(f.negativeCell>=0 && f.positiveCell>=0)check((-f.liquidVolume)+f.liquidVolume==0,"internal face pair does not cancel exactly");
        if(!expectRepair)check(f.liquidPositiveRate==f.rawLiquidPositiveRate && f.liquidNegativeRate==f.rawLiquidNegativeRate,"admissible flux modified");
        if(f.negativeCell<0 || f.positiveCell<0)check(f.liquidPositiveRate==f.rawLiquidPositiveRate && f.liquidNegativeRate==f.rawLiquidNegativeRate,"boundary flux modified");}
    for(std::size_t i=0;i<result.phaseVolumes().size();++i){const auto cell=result.phaseVolumes()[i].cell;cap+=domain.cellVolumes()[i];
        check(result.phaseVolumes()[i].liquid>=0 && result.phaseVolumes()[i].liquid<=domain.cellVolumes()[i],"updated inventory unbounded");
        close(result.phaseVolumes()[i].liquid,previous.primary().phaseVolumes[i].liquid+change[cell],2e-16L*domain.cellVolumes()[i],"M does not use shared flux");}
    close(result.report().boundaryOutward,boundary,1e-14L*cap,"boundary ledger");close(result.report().closure,0,1e-12L*cap,"global closure");
    check(expectRepair?result.report().repairedFaces>0:result.report().repairedFaces==0,"repair count mismatch");
}
long double underLine(Vector2D a,Vector2D b,double height,double slope,bool radial)
{
    std::vector<long double> points{a.x,b.x};if(slope!=0)for(double y:{a.y,b.y}){const long double x=(y-height)/slope;if(x>a.x && x<b.x)points.push_back(x);}
    std::sort(points.begin(),points.end());long double volume=0;
    for(std::size_t k=1;k<points.size();++k){const auto x0=points[k-1],x1=points[k],y=height+slope*(x0+x1)/2;
        if(y<=a.y)continue;const auto c=y>=b.y?b.y-a.y:height-a.y;const auto q=y>=b.y?0:slope;
        volume+=radial?2*pi*(c*(x1*x1-x0*x0)/2+q*(x1*x1*x1-x0*x0*x0)/3):c*(x1-x0)+q*(x1*x1-x0*x0)/2;}return volume;
}
void translationAndPhases()
{
    for(bool radial:{false,true})for(bool amr:{false,true})for(int sign:{-1,1}){
        auto mesh=makeMesh(8,amr);const double dt=.0125,v=.25*sign,height=.421875,slope=radial?0:.25;
        const auto time=HermiteGalsTimeLayer{1,0,dt,Position::Midpoint};
        const double midHeight=height+(radial?v:-slope*v)*dt/2,newHeight=height+(radial?v:-slope*v)*dt;
        const auto profile=[&](double h){return [=](Vector2D p){return Jet{p.y-h-slope*p.x,-slope,1,0};};};
        const auto old=stateFor(*mesh,profile(height),radial,{},[&](int,Vector2D a,Vector2D b){return static_cast<double>(underLine(a,b,height,slope,radial));});
        const auto mid=stateFor(*mesh,profile(midHeight),radial,time,[](int,Vector2D,Vector2D){return 0.;});
        auto input=velocityFor(*mesh,radial,constantVelocity(radial?0:v,radial?v:0));const HermiteGalsFluxVelocity velocity(*mesh,mid,input);
        for(auto& face:input)face.halfPolynomial[0][0]=99; // no alias into owned input
        const auto result=advanceHermiteGalsPhaseVolumes(*mesh,old,mid,velocity);verifyLedger(*mesh,old,mid,result);
        long double error=0,faceError=0;
        const auto domain=domainFor(*mesh,radial);
        for(const auto& m:result.phaseVolumes()){const auto& c=mesh->cells(mesh->cellIndexToLabel[m.cell]);const auto a=mesh->verts(c.cellToVert[0]).coord(),b=mesh->verts(c.cellToVert[2]).coord();
            error=std::max(error,close(m.liquid,underLine(a,b,newHeight,slope,radial),3e-11,"analytic midpoint translation"));}
        for(const auto& flux:result.faces()){const auto face=domain.faces(*mesh)[flux.nativeFace];long double expected=0;
            const auto& edge=face.horizontal?mesh->hEdges(face.edgeLabel):mesh->vEdges(face.edgeLabel);const auto a=mesh->verts(edge.edgeToVert[0]).coord(),b=mesh->verts(edge.edgeToVert[1]).coord();
            if(radial && face.horizontal)expected=face.center.y<midHeight?v*face.area:0;
            if(!radial && !face.horizontal)expected=v*std::clamp(midHeight+slope*face.center.x-std::min(a.y,b.y),0.,std::abs(b.y-a.y));
            faceError=std::max(faceError,close(flux.liquidRate,expected,3e-11,"analytic face rate"));}
        std::cout<<"TRANSLATION,"<<(radial?"RZ":"XY")<<','<<amr<<','<<sign<<','<<result.faces().size()<<','<<error<<','<<faceError<<','<<result.report().closure<<'\n';
        for(int phase:{-1,1}){const Profile constant=[=](Vector2D){return Jet{static_cast<double>(phase),0,0,0};};
            const auto pureOld=stateFor(*mesh,constant,radial,{},[&](int c,Vector2D,Vector2D){return phase<0?domain.cellVolumes()[c]:0.;});
            const auto pureMid=stateFor(*mesh,constant,radial,time,[](int,Vector2D,Vector2D){return 0.;});
            const HermiteGalsFluxVelocity flow(*mesh,pureMid,velocityFor(*mesh,radial,constantVelocity(radial?0:v,radial?v:0)));
            const auto pure=advanceHermiteGalsPhaseVolumes(*mesh,pureOld,pureMid,flow);verifyLedger(*mesh,pureOld,pureMid,pure);
            for(const auto& face:pure.faces())check(face.liquidRate==(phase<0?face.totalRate:0),"pure phase flux is not exact");}
        const HermiteGalsFluxVelocity zero(*mesh,mid,velocityFor(*mesh,radial,constantVelocity(0,0),true));
        const auto stationary=advanceHermiteGalsPhaseVolumes(*mesh,old,mid,zero);verifyLedger(*mesh,old,mid,stationary);
        for(std::size_t i=0;i<stationary.phaseVolumes().size();++i)check(stationary.phaseVolumes()[i].liquid==old.primary().phaseVolumes[i].liquid,"zero velocity changed inventory");
    }
}
void counterflowAndRepair()
{
    auto mesh=makeMesh(4);const auto domain=domainFor(*mesh,false);const Profile liquid=[](Vector2D){return Jet{-1,0,0,0};};
    const auto halfInventory=[&](int c,Vector2D,Vector2D){return .5*domain.cellVolumes()[c];};
    const auto old=stateFor(*mesh,liquid,false,{},halfInventory);
    const auto makeCounter=[&](double dt){return stateFor(*mesh,liquid,false,{1,0,dt,Position::Midpoint},halfInventory);};
    const Velocity counter=[](const NativeLeafFace& f,int half,Vector2D,Vector2D){Coeff p{};if(!f.horizontal)p[0]=half==0?1:-1;return p;};
    const auto mid=makeCounter(.03125);const HermiteGalsFluxVelocity velocity(*mesh,mid,velocityFor(*mesh,false,counter));
    const auto result=advanceHermiteGalsPhaseVolumes(*mesh,old,mid,velocity);verifyLedger(*mesh,old,mid,result);
    check(result.report().counterflowFaces>0 && result.report().maxOutflowCfl>0,"counterflow hidden by signed net");
    for(const auto& f:result.faces())check(f.totalRate==0 && f.liquidRate==0,"symmetric counterflow net not zero");
    const auto large=makeCounter(1.);const HermiteGalsFluxVelocity tooLarge(*mesh,large,velocityFor(*mesh,false,counter));
    rejects([&]{advanceHermiteGalsPhaseVolumes(*mesh,old,large,tooLarge);},"gross outflow CFL");
    std::cout<<"COUNTERFLOW,"<<result.report().counterflowFaces<<','<<result.report().maxOutflowCfl<<','<<result.report().closure<<'\n';
    const double dt=.125,u=.125,cut=.5-u*dt/2;
    const Profile interface=[&](Vector2D p){return Jet{cut-p.x,-1,0,0};};
    const auto previous=stateFor(*mesh,interface,false,{},[&](int c,Vector2D a,Vector2D b){
        if(a.x==.25)return std::nextafter(dt*u*(b.y-a.y),0.);
        return a.x>=.5?domain.cellVolumes()[c]:0.;});
    const auto midpoint=stateFor(*mesh,interface,false,{1,0,dt,Position::Midpoint},halfInventory);
    const HermiteGalsFluxVelocity flow(*mesh,midpoint,velocityFor(*mesh,false,constantVelocity(u,0)));
    const auto repaired=advanceHermiteGalsPhaseVolumes(*mesh,previous,midpoint,flow);verifyLedger(*mesh,previous,midpoint,repaired,true);
    check(repaired.report().maxRelativeFluxCorrection<=64*std::numeric_limits<double>::epsilon(),"roundoff repair exceeds budget");
    std::cout<<"REPAIR,"<<repaired.report().roundoffRepairSweeps<<','<<repaired.report().repairedFaces<<','<<repaired.report().maxRelativeFluxCorrection<<','<<repaired.report().closure<<'\n';
    const auto empty=stateFor(*mesh,interface,false,{},[&](int c,Vector2D a,Vector2D){return a.x>=.5?domain.cellVolumes()[c]:0.;});
    rejects([&]{advanceHermiteGalsPhaseVolumes(*mesh,empty,midpoint,flow);},"violation exceeds roundoff budget");
    const Profile boundaryCut=[&](Vector2D p){return Jet{1-u*dt/2-p.x,-1,0,0};};
    const auto lastOld=stateFor(*mesh,boundaryCut,false,{},[&](int,Vector2D a,Vector2D b){return a.x==.75?std::nextafter(dt*u*(b.y-a.y),0.):0.;});
    const auto lastMid=stateFor(*mesh,boundaryCut,false,{1,0,dt,Position::Midpoint},halfInventory);
    const HermiteGalsFluxVelocity prescribed(*mesh,lastMid,velocityFor(*mesh,false,constantVelocity(u,0)));
    rejects([&]{advanceHermiteGalsPhaseVolumes(*mesh,lastOld,lastMid,prescribed);},"no interior flux capacity");
}
void multiRootAndConvergence()
{
    for(bool radial:{false,true}){
        auto mesh=makeMesh(4,true);const auto domain=domainFor(*mesh,radial);
        // Two roots in one actual horizontal face, with liquid between them.
        const double a=.27,b=.31;const Profile band=[=](Vector2D p){
            if(!radial)return Jet{(p.x-a)*(p.x-b),2*p.x-a-b,0,0};
            return Jet{(p.x*p.x-a*a)*(p.x*p.x-b*b),2*p.x*(2*p.x*p.x-a*a-b*b),0,0};};
        const auto inventory=[&](int,Vector2D lo,Vector2D hi){const auto x0=std::max(lo.x,a),x1=std::min(hi.x,b);
            return x1>x0?static_cast<double>((radial?pi*(x1*x1-x0*x0):x1-x0)*(hi.y-lo.y)):0.;};
        const auto old=stateFor(*mesh,band,radial,{},inventory),mid=stateFor(*mesh,band,radial,{1,0,.01,Position::Midpoint},inventory);
        const HermiteGalsFluxVelocity flow(*mesh,mid,velocityFor(*mesh,radial,constantVelocity(0,1)));
        const auto result=advanceHermiteGalsPhaseVolumes(*mesh,old,mid,flow);verifyLedger(*mesh,old,mid,result);
        long double maximum=0;
        for(const auto& f:result.faces()){const auto face=domain.faces(*mesh)[f.nativeFace];long double ref=0;
            if(face.horizontal){const auto& e=mesh->hEdges(face.edgeLabel);const auto x0=std::max(a,mesh->verts(e.edgeToVert[0]).coord().x),x1=std::min(b,mesh->verts(e.edgeToVert[1]).coord().x);
                if(x1>x0)ref=radial?pi*(x1*x1-x0*x0):x1-x0;}maximum=std::max(maximum,close(f.liquidRate,ref,3e-11,"multi-root phase flux"));}
        std::cout<<"MULTI_ROOT,"<<(radial?"RZ":"XY")<<','<<maximum<<','<<result.report().phaseRanges<<'\n';
        const double reversal=.2875;
        const Velocity reversed=[=](const NativeLeafFace& face,int,Vector2D lo,Vector2D hi){Coeff p{};if(face.horizontal){const auto h=hi.x-lo.x;
            if(radial){p[0]=lo.x*lo.x-reversal*reversal;p[1]=2*lo.x*h;p[2]=h*h;}
            else{p[0]=lo.x-reversal;p[1]=h;}}return p;};
        const HermiteGalsFluxVelocity counter(*mesh,mid,velocityFor(*mesh,radial,reversed));
        const auto bidirectional=advanceHermiteGalsPhaseVolumes(*mesh,old,mid,counter);verifyLedger(*mesh,old,mid,bidirectional);
        long double counterError=0;std::size_t mixedDirections=0;
        for(const auto& f:bidirectional.faces()){const auto face=domain.faces(*mesh)[f.nativeFace];long double ref=0;
            if(face.horizontal){const auto& e=mesh->hEdges(face.edgeLabel);const auto x0=std::max(a,mesh->verts(e.edgeToVert[0]).coord().x),x1=std::min(b,mesh->verts(e.edgeToVert[1]).coord().x);
                if(x1>x0)ref=radial?pi*((std::pow(x1,4)-std::pow(x0,4))/2-reversal*reversal*(x1*x1-x0*x0)):
                    (x1*x1-x0*x0)/2-reversal*(x1-x0);}
            counterError=std::max(counterError,close(f.liquidRate,ref,3e-11,"within-half counterflow phase flux"));
            mixedDirections+=f.liquidPositiveRate>0 && f.liquidNegativeRate>0;}
        check(mixedDirections>0,"phase counterflow roots not exercised");
        std::cout<<"POLYNOMIAL_COUNTERFLOW,"<<(radial?"RZ":"XY")<<','<<mixedDirections<<','<<counterError<<','<<bidirectional.report().maxRelativeDivergence<<'\n';
        long double previousError=0,lastOrder=0;
        for(int n:{8,16,32}){auto grid=makeMesh(n);const double radius=1./3;
            const Profile scalar=[=](Vector2D p){const double q=radial?p.x*p.x:.7*p.x,target=radial?radius*radius:.7*radius,e=std::exp(q);
                return Jet{e-std::exp(target),e*(radial?2*p.x:.7),0,0};};
            const auto gd=domainFor(*grid,radial);const auto inventoryM=[&](int c,Vector2D,Vector2D){return .5*gd.cellVolumes()[c];};
            const auto source=stateFor(*grid,scalar,radial,{},inventoryM),middle=stateFor(*grid,scalar,radial,{1,0,.001,Position::Midpoint},inventoryM);
            const HermiteGalsFluxVelocity v(*grid,middle,velocityFor(*grid,radial,constantVelocity(0,1)));
            const auto measured=advanceHermiteGalsPhaseVolumes(*grid,source,middle,v);verifyLedger(*grid,source,middle,measured);
            long double topFlux=0;for(const auto& f:measured.faces()){const auto face=gd.faces(*grid)[f.nativeFace];if(face.horizontal && face.center.y==1)topFlux+=f.liquidRate;}
            const auto error=std::abs(topFlux-(radial?pi*radius*radius:radius));const auto order=previousError>0?std::log2(previousError/error):0;
            if(previousError>0)check(error<previousError,"face-cut error not decreasing");
            std::cout<<"CONVERGENCE,"<<(radial?"RZ":"XY")<<','<<n<<','<<error<<','<<order<<'\n';previousError=error;lastOrder=order;
        }check(lastOrder>=3,"face-cut last order below three");
    }
}
void failuresAndLifetime()
{
    auto mesh=makeMesh(4);const auto domain=domainFor(*mesh,false);const Profile liquid=[](Vector2D){return Jet{-1,0,0,0};};
    const auto half=[&](int c,Vector2D,Vector2D){return .5*domain.cellVolumes()[c];};
    const auto old=stateFor(*mesh,liquid,false,{},half),mid=stateFor(*mesh,liquid,false,{1,0,.01,Position::Midpoint},half);
    const auto input=velocityFor(*mesh,false,constantVelocity(0,0),true);const auto oldId=old.identity(),midId=mid.identity();
    auto bad=input;bad.pop_back();rejects([&]{HermiteGalsFluxVelocity v(*mesh,mid,bad);},"missing/duplicate/unordered");
    bad=input;bad[1].nativeFace=bad[0].nativeFace;rejects([&]{HermiteGalsFluxVelocity v(*mesh,mid,bad);},"missing/duplicate/unordered");
    bad=input;bad[0].halfPolynomial[0][0]=std::numeric_limits<double>::quiet_NaN();rejects([&]{HermiteGalsFluxVelocity v(*mesh,mid,bad);},"nonfinite");
    bad=input;const auto interior=std::find_if(bad.begin(),bad.end(),[](const auto& f){return f.boundary==Boundary::Interior;});
    check(interior!=bad.end(),"missing interior boundary fixture");interior->boundary=Boundary::Open;
    rejects([&]{HermiteGalsFluxVelocity v(*mesh,mid,bad);},"boundary kind");
    // A zero NET boundary flow is not impermeability.
    bad=input;auto outer=std::find_if(bad.begin(),bad.end(),[](const auto& f){return f.boundary==Boundary::Impermeable;});
    outer->halfPolynomial[0][0]=1;outer->halfPolynomial[1][0]=-1;
    rejects([&]{HermiteGalsFluxVelocity v(*mesh,mid,bad);},"identically zero");
    rejects([&]{HermiteGalsFluxVelocity v(*mesh,old,input);},"requires midpoint");
    const HermiteGalsFluxVelocity velocity(*mesh,mid,input);
    const auto foreign=stateFor(*mesh,liquid,false,{1,0,.01,Position::Midpoint},half);
    rejects([&]{advanceHermiteGalsPhaseVolumes(*mesh,old,foreign,velocity);},"foreign velocity");
    auto changed=makeMesh(4,true);rejects([&]{advanceHermiteGalsPhaseVolumes(*changed,old,mid,velocity);},"mesh");
    const auto wrongTime=stateFor(*mesh,liquid,false,{3,0,.01,Position::Midpoint},half);
    const HermiteGalsFluxVelocity wrongVelocity(*mesh,wrongTime,input);
    rejects([&]{advanceHermiteGalsPhaseVolumes(*mesh,old,wrongTime,wrongVelocity);},"time layer");
    HermiteGalsPhaseFluxConfig config;config.relativeIntegrationTolerance=0;
    rejects([&]{advanceHermiteGalsPhaseVolumes(*mesh,old,mid,velocity,config);},"relative tolerance");
    auto divergent=velocityFor(*mesh,false,[](const NativeLeafFace& f,int,Vector2D,Vector2D){Coeff p{};if(!f.horizontal)p[0]=f.center.x;return p;});
    const HermiteGalsFluxVelocity divergence(*mesh,mid,divergent);
    rejects([&]{advanceHermiteGalsPhaseVolumes(*mesh,old,mid,divergence);},"not divergence free");
    const auto overfilled=stateFor(*mesh,liquid,false,{},[&](int c,Vector2D,Vector2D){return 1.1*domain.cellVolumes()[c];});
    rejects([&]{advanceHermiteGalsPhaseVolumes(*mesh,overfilled,mid,velocity);},"old inventory outside capacity");
    auto kept=[&]{const auto local=stateFor(*mesh,liquid,false,{1,0,.01,Position::Midpoint},half);return advanceHermiteGalsPhaseVolumes(*mesh,old,local,HermiteGalsFluxVelocity(*mesh,local,input));}();
    auto copy=kept;check(&copy.geometry().snapshot()==&kept.geometry().snapshot(),"retained result copied geometry");
    check(old.identity()==oldId && mid.identity()==midId,"failure mutated source identity");
    for(std::size_t k=0;k<old.primary().phaseVolumes.size();++k)check(old.primary().phaseVolumes[k].liquid==.5*domain.cellVolumes()[k],"failure changed old M");
    std::cout<<"LIFECYCLE,PASS\n";
}
void tangentFullPhase()
{
    for(bool radial:{false,true})for(bool amr:{false,true})for(bool interiorRoot:{false,true}){auto mesh=makeMesh(4,amr);const auto domain=domainFor(*mesh,radial);const double r=interiorRoot?.28125:.25;
        const Profile shape=[=](Vector2D p){const auto q=radial?p.x*p.x-r*r:p.x-r;return Jet{p.y-q*q,radial?-4*p.x*q:-2*q,1,0};};
        const auto inventory=[&](int c,Vector2D,Vector2D){return .5*domain.cellVolumes()[c];};
        const auto old=stateFor(*mesh,shape,radial,{},inventory),mid=stateFor(*mesh,shape,radial,{1,0,.001,Position::Midpoint},inventory);
        const HermiteGalsFluxVelocity flow(*mesh,mid,velocityFor(*mesh,radial,constantVelocity(0,1)));
        const auto result=advanceHermiteGalsPhaseVolumes(*mesh,old,mid,flow);verifyLedger(*mesh,old,mid,result);
        for(const auto& f:result.faces()){const auto face=domain.faces(*mesh)[f.nativeFace];
            if(face.horizontal && face.center.y==0)check(f.liquidRate==f.totalRate,"isolated tangent root lost full-liquid identity");}
        std::cout<<"TANGENT_FULL_PHASE,"<<(radial?"RZ":"XY")<<','<<amr<<','<<interiorRoot<<",PASS\n";
    }
}
void tangentMappedRoots()
{
    std::size_t configurations=0;
    for(bool amr:{false,true})for(int multiplicity:{2,3,4})for(int sign:{-1,1})for(bool reverse:{false,true}){
        auto mesh=makeMesh(4,amr);const Profile f=[=](Vector2D p){const auto q=3*p.x-1;double value=1,derivative=1;
            for(int k=0;k<multiplicity;++k)value*=q;for(int k=1;k<multiplicity;++k)derivative*=q;
            return Jet{sign*(p.y+value),sign*3*multiplicity*derivative,double(sign),0};};
        const auto state=stateFor(*mesh,f,false,{},[](int,Vector2D a,Vector2D b){return .5*capacity(a,b,false);});
        const auto owner=std::make_shared<const HermiteGalsGeometryView>(state.geometry());
        const std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source(owner,&owner->snapshot());
        const Vector2D a{.25,0},b{.5,0};const LevelSetMappedSegment segment(source,reverse?b:a,reverse?a:b);
        check(segment.roots().size()==1 && !segment.rootAtInputEndpoint(0,false) && !segment.rootAtInputEndpoint(0,true),"rational multiple-root identity");
        long double negative=0,total=0;for(const auto& r:segment.ranges()){total+=r.length;if(r.negative)negative+=r.length;}
        const long double expected=multiplicity%2?(sign>0?1.L/3:2.L/3):(sign<0?1.L:0.L);
        close(total,1,1e-15,"tangent phase coverage");close(negative,expected,1e-14,"multiple-root phase parity");
        const auto root=segment.roots()[0];check(root.lower<root.upper,"nonbinary root must retain a bracket");
        const auto refined=segment.refined(0);check(refined.roots().size()==1 && segment.roots()[0].lower==root.lower && segment.roots()[0].upper==root.upper,"multiple-root refinement mutated input");
        negative=0;for(const auto& r:refined.ranges())if(r.negative)negative+=r.length;close(negative,expected,1e-14,"refinement changed tangent phase measure");
        rejects([&]{segment.nativeQuarticNegativePhaseVariation();},"transverse");
        ++configurations;
    }
    for(int sign:{-1,1}){auto mesh=makeMesh(4);const double gap=sign*std::ldexp(1.,-20);
        const Profile f=[=](Vector2D p){const auto q=3*p.x-1;return Jet{p.y+q*q+gap,6*q,1,0};};
        const auto state=stateFor(*mesh,f,false,{},[](int,Vector2D a,Vector2D b){return .5*capacity(a,b,false);});
        const auto owner=std::make_shared<const HermiteGalsGeometryView>(state.geometry());const std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source(owner,&owner->snapshot());
        const LevelSetMappedSegment segment(source,{.25,0},{.5,0});check(segment.roots().size()==(sign<0?2u:0u),"near tangent root inventory changed");
        long double negative=0;for(const auto& r:segment.ranges())if(r.negative)negative+=r.length;
        close(negative,sign<0?8*std::sqrt(-gap)/3:0,1e-13,"near tangent phase width");
        check(segment.nativeQuarticNegativePhaseVariation().size()==segment.roots().size(),"simple near-tangent variation unavailable");}
    for(bool reverse:{false,true}){auto mesh=makeMesh(4);const Profile f=[](Vector2D p){const auto q=p.x*p.x-p.x+.125;
        return Jet{p.y+q*q,2*q*(2*p.x-1),1,0};};
        const auto state=stateFor(*mesh,f,false,{},[](int,Vector2D a,Vector2D b){return .5*capacity(a,b,false);});
        const auto owner=std::make_shared<const HermiteGalsGeometryView>(state.geometry());const std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source(owner,&owner->snapshot());
        const LevelSetMappedSegment segment(source,reverse?Vector2D{1,0}:Vector2D{0,0},reverse?Vector2D{0,0}:Vector2D{1,0});
        check(segment.roots().size()==2,"irrational double roots lost");long double total=0;
        for(const auto& r:segment.ranges()){check(!r.negative,"even irrational roots flipped phase");total+=r.length;}close(total,1,1e-15,"irrational phase coverage");
        const auto refined=segment.refined(0);check(refined.roots().size()==2,"irrational root refinement lost inventory");
        rejects([&]{refined.nativeQuarticNegativePhaseVariation();},"transverse");}
    {auto mesh=makeMesh(4);const auto volume=[](int,Vector2D a,Vector2D b){return .5*capacity(a,b,false);};
        const auto zero=stateFor(*mesh,[](Vector2D p){return Jet{p.y,0,1,0};},false,{},volume);
        const auto zeroOwner=std::make_shared<const HermiteGalsGeometryView>(zero.geometry());const std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> z(zeroOwner,&zeroOwner->snapshot());
        rejects([&]{LevelSetMappedSegment line(z,{.25,0},{.5,0});},"positive-length coincident zero trace");
        // Isolated contact belongs to the same complete closed zero set.
        const auto circle=stateFor(*mesh,[](Vector2D p){const auto x=p.x-.5,y=p.y-.53125;return Jet{x*x+y*y-.25*.25,2*x,2*y,0};},false,{},volume);
        const auto circleOwner=std::make_shared<const HermiteGalsGeometryView>(circle.geometry());const std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> c(circleOwner,&circleOwner->snapshot());
        const auto topology=LevelSetInterfacePanels::discoverNativeComponents(*mesh,c,false);
        check(topology.curves.size()==1 && topology.curves.front().closed(),"tangent circle lost complete topology");
        const auto& curve=topology.curves.front();
        for(std::size_t g=0;g<curve.graphs().size();++g)close(curve.frameAt(static_cast<int>(g),.5).curvature,4,1e-10,"tangent circle curvature");}
    std::cout<<"TANGENT_MAPPED_ROOTS,"<<configurations<<",irrational=2,zero_trace_guard=1,closed_tangent=1,PASS\n";
}
void directionalSubsetIdentity()
{
    std::size_t identities=0,reconciliations=0,nearRootPositive=0;long double maxAdjustment=0,maxComplementError=0;
    const double k=40*std::cos(double(pi)*.00025/.4);
    const std::vector<double> cuts{.4375+.0625*.1,.4375+.0625*.1528854156359154,.4375+.0625*.7,
        std::nextafter(.5,0.),std::nextafter(.5,1.),.57};
    for(bool radial:{false,true})for(bool amr:{false,true})for(int sign:{-1,1})for(double cut:cuts){
        auto mesh=makeMesh(8,amr);const auto domain=domainFor(*mesh,radial);
        const Profile shape=[=](Vector2D p){return radial?Jet{cut*cut-p.x*p.x,-2*p.x,0,0}:Jet{cut-p.y,0,-1,0};};
        const auto inventory=[&](int c,Vector2D,Vector2D){return .5*domain.cellVolumes()[c];};
        const HermiteGalsTimeLayer time{1,0,1e-5,Position::Midpoint};
        const auto old=stateFor(*mesh,shape,radial,{},inventory),mid=stateFor(*mesh,shape,radial,time,inventory);
        const auto input=velocityFor(*mesh,radial,[=](const NativeLeafFace& f,int,Vector2D a,Vector2D b){Coeff c{};
            if(f.horizontal!=radial)return c;
            const double h=radial?b.x-a.x:b.y-a.y,t=radial?a.x:a.y-.5,q=sign*k;
            if(radial){c[0]=-q*(t*t-.25)*(t*t+.0625);c[1]=q*(.375*t-4*t*t*t)*h;
                c[2]=q*(.1875-6*t*t)*h*h;c[3]=-4*q*t*h*h*h;c[4]=-q*h*h*h*h;}
            else {c[0]=-q*(.0625+t*t)*t;c[1]=-q*(.0625+3*t*t)*h;c[2]=-3*q*t*h*h;c[3]=-q*h*h*h;}
            return c;});
        const auto result=advanceHermiteGalsPhaseVolumes(*mesh,old,mid,HermiteGalsFluxVelocity(*mesh,mid,input));verifyLedger(*mesh,old,mid,result);
        identities+=result.report().directionalIdentityCount;reconciliations+=result.report().subsetReconciliationCount;
        maxAdjustment=std::max(maxAdjustment,result.report().maxRelativeSubsetAdjustment);
        const Profile complement=[=](Vector2D p){auto j=shape(p);return Jet{-j.value,-j.dx,-j.dy,-j.dxy};};
        const auto opposite=stateFor(*mesh,complement,radial,time,inventory);
        const auto gas=advanceHermiteGalsPhaseVolumes(*mesh,old,opposite,HermiteGalsFluxVelocity(*mesh,opposite,input));
        verifyLedger(*mesh,old,opposite,gas);
        identities+=gas.report().directionalIdentityCount;reconciliations+=gas.report().subsetReconciliationCount;
        maxAdjustment=std::max(maxAdjustment,gas.report().maxRelativeSubsetAdjustment);
        const Profile empty=[](Vector2D){return Jet{1,0,0,0};};const auto pure=stateFor(*mesh,empty,radial,time,inventory);
        const auto reference=advanceHermiteGalsPhaseVolumes(*mesh,old,pure,HermiteGalsFluxVelocity(*mesh,pure,input));
        for(std::size_t i=0;i<result.faces().size();++i){const auto& f=result.faces()[i];const auto& g=gas.faces()[i];const auto& ref=reference.faces()[i];
            check(f.totalPositiveRate==ref.totalPositiveRate && f.totalNegativeRate==ref.totalNegativeRate && f.totalRate==ref.totalRate,"phase changed prescribed total flux");
            check(g.totalPositiveRate==ref.totalPositiveRate && g.totalNegativeRate==ref.totalNegativeRate && g.totalRate==ref.totalRate,"complement changed prescribed total flux");
            const auto gross=f.totalPositiveRate+f.totalNegativeRate;
            const auto error=std::max(std::abs(f.liquidPositiveRate+g.liquidPositiveRate-f.totalPositiveRate),std::abs(f.liquidNegativeRate+g.liquidNegativeRate-f.totalNegativeRate));
            check(error<=4e-13L*gross,"complementary phase integral discrepancy");if(gross>0)maxComplementError=std::max(maxComplementError,error/gross);
            const auto cap=f.negativeCell<0?domain.cellVolumes()[f.positiveCell]:f.positiveCell<0?domain.cellVolumes()[f.negativeCell]:
                std::min(domain.cellVolumes()[f.negativeCell],domain.cellVolumes()[f.positiveCell]);
            for(const auto* phase:{&f,&g}){
                const auto adjusted=std::max(phase->integratedLiquidPositiveRate-phase->rawLiquidPositiveRate,phase->integratedLiquidNegativeRate-phase->rawLiquidNegativeRate);
                check(adjusted>=0 && adjusted<=64*std::numeric_limits<double>::epsilon()*gross,"unrecorded/excessive subset adjustment");
                check((time.end-time.begin)*adjusted<=64*std::numeric_limits<double>::epsilon()*cap,"subset adjustment exceeds inventory roundoff budget");
            }
            if(std::abs(cut-.5)<1e-15 && (f.liquidPositiveRate>0 || f.liquidNegativeRate>0))++nearRootPositive;
        }
        std::cout<<"DIRECTIONAL_CASE,"<<(radial?"RZ":"XY")<<','<<amr<<','<<sign<<','<<cut<<','<<result.report().directionalIdentityCount<<','<<result.report().subsetReconciliationCount<<','<<gas.report().directionalIdentityCount<<','<<gas.report().subsetReconciliationCount<<','<<std::max(result.report().maxRelativeSubsetAdjustment,gas.report().maxRelativeSubsetAdjustment)<<'\n';
    }
    check(identities>0 && reconciliations>identities && nearRootPositive>0,"directional/complementary subset cases not exercised");
    check(maxAdjustment<=64*std::numeric_limits<double>::epsilon(),"subset adjustment budget");
    std::cout<<"DIRECTIONAL_SUMMARY,"<<identities<<','<<reconciliations<<','<<nearRootPositive<<','<<maxAdjustment<<','<<maxComplementError<<'\n';
}
using Cache=InterfacePhaseGeometry::MappedFacePartitions;
using Trace=LevelSetMappedSegment::NativeQuarticTrace;
using Snapshot=LevelSetHermiteGeometry::Snapshot;
static_assert(!std::is_copy_constructible_v<Cache> && !std::is_copy_constructible_v<Trace>);
long double traceMeasure(const std::vector<LevelSetMappedSegment::Range>& ranges,Vector2D a,Vector2D b,bool rz,bool negativeOnly=true)
{
    const long double dx=static_cast<long double>(b.x)-a.x,dy=static_cast<long double>(b.y)-a.y,ds=std::hypot(dx,dy);
    long double total=0;for(const auto& r:ranges)if(!negativeOnly || r.negative){const auto t=r.fromEnd?1-r.offset-r.length/2:r.offset+r.length/2;
        total+=ds*r.length*(rz?2*pi*(a.x+dx*t):1);}return total;
}
std::shared_ptr<const Snapshot> selectedTraceSource(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const Profile& f)
{
    Patch patch;patch.activeCells.assign(mesh.cellIndexNum(),0);patch.nodes.resize(mesh.vertIndexNum());patch.cellCenters.resize(mesh.cellIndexNum());
    for(int v=0;v<mesh.vertIndexNum();++v)patch.nodes[v]=f(point(mesh,v));std::set<Key> keys;
    for(int c:domain.cellIndexByRow()){const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);patch.activeCells[c]=1;patch.cellCenters[c]=f(cell.centerCoord()).value;
        for(int k=0;k<4;++k){Key e{{mesh.verts(cell.cellToVert[k]).index(),mesh.verts(cell.cellToVert[(k+1)%4]).index()}};
            const auto a=point(mesh,e[0]),b=point(mesh,e[1]);if(a.y==b.y?a.x>b.x:a.y>b.y)std::swap(e[0],e[1]);keys.insert(e);}}
    for(const auto& e:keys){const auto a=point(mesh,e[0]),b=point(mesh,e[1]),mid=a+(b-a)*.5;const auto j=f(mid);
        patch.edgeMidpoints.push_back({e,j.value,a.y==b.y?j.dy:j.dx});}
    return std::make_shared<const Snapshot>(Snapshot::fromNativeLeaves(mesh,patch,domain.axisymmetric()));
}
void traceCacheContract()
{
    int layouts=0,halves=0,zeros=0,hanging=0,refinements=0;long double peak=0;
    for(bool rz:{false,true})for(bool amr:{false,true})for(int axis:{0,1})for(int sign:{-1,1}){
        auto mesh=makeMesh(4,amr);auto domain=domainFor(*mesh,rz);const double line=.25;
        const Profile f=[=](Vector2D p){const double t=axis==0?p.y:p.x,n=rz && axis==0?p.x*p.x-line*line:(axis==0?p.x:p.y)-line;
            const double g=rz && axis==1?1+t*t:1+t+t*t,dg=rz && axis==1?2*t:1+2*t,dn=rz && axis==0?2*p.x:1;
            return Jet{sign*n*g,sign*(axis==0?dn*g:n*dg),sign*(axis==0?n*dg:dn*g),sign*dn*dg};};
        const auto seed=stateFor(*mesh,f,rz,{},[](int,Vector2D,Vector2D){return 0.;},HermiteGalsTransport::FirstJet);
        auto owner=std::make_shared<const HermiteGalsGeometryView>(seed.geometry());std::shared_ptr<const Snapshot> source(owner,&owner->snapshot());
        rejects([&]{Cache old(*mesh,domain,source);},"positive-length coincident zero trace");
        auto cache=Cache::prepareNativeTrace(*mesh,domain,source);cache.validate(*mesh,domain);
        check(&cache.source()==source.get() && seed.primary().phaseVolumes.empty(),"typed cache lost pure source");
        long double measure=0,minus=0,plus=0;std::size_t events=0,values=0;
        for(const auto face:domain.faces(*mesh))for(int half=0;half<2;++half){const auto& q=cache.traceHalfFace(face.nativeRow,half);++halves;
            const auto& edge=face.horizontal?mesh->hEdges(face.edgeLabel):mesh->vEdges(face.edgeLabel);
            auto a=mesh->verts(edge.edgeToVert[0]).coord(),b=mesh->verts(edge.edgeToVert[1]).coord();if(face.horizontal?a.x>b.x:a.y>b.y)std::swap(a,b);
            const auto begin=half==0?a:face.center,end=half==0?face.center:b;
            const long double halfArea=rz?(face.horizontal?pi*(static_cast<long double>(end.x)-begin.x)*(static_cast<long double>(end.x)+begin.x):2*pi*begin.x*(static_cast<long double>(end.y)-begin.y)):
                std::hypot(static_cast<long double>(end.x)-begin.x,static_cast<long double>(end.y)-begin.y);
            check(q.usesSource(*source) && q.normalAxis()==(face.horizontal?1:0) && q.begin().x==begin.x && q.begin().y==begin.y && q.end().x==end.x && q.end().y==end.y,"actual typed half-face addressing");
            events+=q.events().size();for(const auto& e:q.events())if(e.normalOrder==0)++values;
            auto reverse=LevelSetMappedSegment::nativeQuarticTrace(source,end,begin,q.normalAxis());
            if(q.coincidentRanges().empty()){const LevelSetMappedSegment old(source,begin,end);const auto& legacy=cache.halfFace(face.nativeRow,half);
                check(&legacy==&q.ordinary() && legacy.ranges().size()==old.ranges().size() && legacy.roots().size()==old.roots().size(),"ordinary cache forked query/inventory");
                for(std::size_t k=0;k<old.ranges().size();++k){const auto& a=old.ranges()[k];const auto& b=legacy.ranges()[k];check(a.fromEnd==b.fromEnd && a.offset==b.offset && a.length==b.length && a.negative==b.negative,"ordinary cached range changed");}
                close(traceMeasure(reverse.ordinaryRanges(),end,begin,rz),traceMeasure(q.ordinaryRanges(),begin,end,rz),1e-12L*face.area,"ordinary reversed physical measure");
            }else{++zeros;rejects([&]{cache.halfFace(face.nativeRow,half);},"no legacy phase view");
                long double ml=0,pl=0,rml=0,rpl=0;
                for(const auto& z:q.coincidentRanges()){check(z.minus.available && z.plus.available && z.minus.normalOrder==1 && z.plus.normalOrder==1,"regular actual two-sided half-face phase");
                    ml+=traceMeasure(z.minus.ranges,begin,end,rz);pl+=traceMeasure(z.plus.ranges,begin,end,rz);
                    measure+=traceMeasure(z.minus.ranges,begin,end,rz,false);
                    const auto m=source->nativeCellBounds(z.minus.nativeCellIndex),p=source->nativeCellBounds(z.plus.nativeCellIndex);
                    if(m[1].x-m[0].x!=p[1].x-p[0].x)++hanging;}
                for(const auto& z:reverse.coincidentRanges()){rml+=traceMeasure(z.minus.ranges,end,begin,rz);rpl+=traceMeasure(z.plus.ranges,end,begin,rz);}
                peak=std::max(peak,close(ml,sign>0?halfArea:0,1e-12L*halfArea,"minus actual half-face reference")/halfArea);
                peak=std::max(peak,close(pl,sign<0?halfArea:0,1e-12L*halfArea,"plus actual half-face reference")/halfArea);
                close(rml,ml,1e-12L*face.area,"reverse minus phase");close(rpl,pl,1e-12L*face.area,"reverse plus phase");minus+=ml;plus+=pl;
            }
        }
        const long double expected=rz?(axis==0?2*pi*line:pi):1;
        close(measure,expected,1e-12L*expected,"real zero half/subface partition additivity");close(minus,sign>0?expected:0,1e-12L*expected,"complete minus normal phase");close(plus,sign<0?expected:0,1e-12L*expected,"complete plus normal phase");
        check(cache.rootIncidenceCount()==values && cache.traceEventIncidenceCount()==events,"typed value/normal count confused");
        auto duplicate=domain;rejects([&]{cache.validate(*mesh,duplicate);},"different fitted domain");
        rejects([&]{cache.traceHalfFace(-1,0);},"invalid mapped");rejects([&]{cache.refineTraceEvent(0,0,events+1);},"event index");
        std::weak_ptr<const Snapshot> life=source;source.reset();owner.reset();check(!life.expired(),"typed cache lost owning view");
        auto moved=std::move(cache);rejects([&]{cache.source();},"moved-from");moved.validate(*mesh,domain);++layouts;
    }
    check(hanging>0,"typed cache omitted actual hanging subfaces");
    auto mesh=makeMesh(4);auto domain=domainFor(*mesh,false);
    // A nonbinary normal coefficient root, without a regular-interface claim.
    const Profile f=[](Vector2D p){return Jet{(p.x-.5)*(3*p.y-1),3*p.y-1,3*(p.x-.5),3};};
    auto source=selectedTraceSource(*mesh,domain,f);auto cache=Cache::prepareNativeTrace(*mesh,domain,source);
    bool normal=false,value=false;
    for(const auto face:domain.faces(*mesh))for(int half=0;half<2;++half){const auto& q=cache.traceHalfFace(face.nativeRow,half);
        for(std::size_t k=0;k<q.events().size();++k){const auto e=q.events()[k];if(e.root.pointOnly || !(e.root.lower<e.root.upper))continue;
            if((e.normalOrder>0 && normal) || (e.normalOrder==0 && value))continue;
            const auto count=cache.traceEventIncidenceCount(),rootCount=cache.rootIncidenceCount();const auto before=&q;
            const auto ml=q.coincidentRanges().empty()?0:traceMeasure(q.coincidentRanges()[0].plus.ranges,q.begin(),q.end(),false);
            rejects([&]{cache.refineTraceEvent(face.nativeRow,half,q.events().size());},"event index");check(&cache.traceHalfFace(face.nativeRow,half)==before,"failed typed refinement replaced original");
            if(e.normalOrder==0)cache.refineRoot(face.nativeRow,half,k);else cache.refineTraceEvent(face.nativeRow,half,k);
            const auto& newQuery=cache.traceHalfFace(face.nativeRow,half);check(newQuery.usesSource(*source) && &newQuery!=before && cache.traceEventIncidenceCount()==count && cache.rootIncidenceCount()==rootCount,"typed refinement provenance/count/commit");
            check(newQuery.events()[k].root.chart>e.root.chart,"typed refinement did not retain its new same-root child chart");
            if(e.normalOrder>0){normal=true;const auto actual=traceMeasure(newQuery.coincidentRanges()[0].plus.ranges,newQuery.begin(),newQuery.end(),false);
                const long double expected=std::max(0.L,std::min(static_cast<long double>(newQuery.end().y),1.L/3)-newQuery.begin().y);
                close(actual,ml,1e-12,"recomputed normal phase after refinement");close(actual,expected,1e-12L*(newQuery.end().y-newQuery.begin().y),"independent nonbinary normal phase reference");}else{
                value=true;const auto child=newQuery.events()[k].root.chart;cache.refineRoot(face.nativeRow,half,k);++refinements;
                check(cache.traceHalfFace(face.nativeRow,half).events()[k].root.chart>child && cache.traceEventIncidenceCount()==count,"repeated ordinary refinement forked the typed root");}
            ++refinements;break;
        }}
    check(normal && value,"typed normal and ordinary same-source refinements not exercised");
    for(int order:{2,3,4}){const Profile power=[=](Vector2D p){const double n=p.x-.5;return Jet{std::pow(n,order),order*std::pow(n,order-1),0,0};};
        const auto scalar=selectedTraceSource(*mesh,domain,power);auto typed=Cache::prepareNativeTrace(*mesh,domain,scalar);bool found=false;
        for(const auto face:domain.faces(*mesh))for(int half=0;half<2;++half)for(const auto& z:typed.traceHalfFace(face.nativeRow,half).coincidentRanges()){
            check(z.minus.normalOrder==order && z.plus.normalOrder==order,"typed higher normal order lost");found=true;}check(found,"higher normal support absent");}
    const Profile mixed=[](Vector2D p){const double t=std::max(0.,p.y-.5);return Jet{p.x-.5+t*t,1,2*t,0};};
    DictionaryMesh md;md.Nx=md.Ny=2;md.maxCellNum=64;md.maxInitialIterations=0;QuadtreeMesh mixedMesh(md);auto mixedDomain=domainFor(mixedMesh,false);
    const auto mixedSource=selectedTraceSource(mixedMesh,mixedDomain,mixed);auto mixedCache=Cache::prepareNativeTrace(mixedMesh,mixedDomain,mixedSource);
    bool below=false,above=false;
    for(const auto face:mixedDomain.faces(mixedMesh))for(int half=0;half<2;++half){const auto& q=mixedCache.traceHalfFace(face.nativeRow,half);
        if(q.normalAxis()!=0 || q.begin().x!=.5)continue;
        below=below || (!q.coincidentRanges().empty() && q.end().y==.5);
        above=above || (!q.ordinaryRanges().empty() && q.begin().y==.5);}
    check(below && above,"ordinary-to-zero actual half-face incidence not exercised");
    FaceSegmentGraph selectedGraph;selectedGraph.build(*mesh,[](const QuadtreeMesh&,const Cell& c){return c.centerCoord().x<.5;});
    auto selected=NativeLeafDomain::build(*mesh,selectedGraph,false);const Profile plane=[](Vector2D p){return Jet{p.x-.5,1,0,0};};
    const auto partial=selectedTraceSource(*mesh,selected,plane);auto selectedCache=Cache::prepareNativeTrace(*mesh,selected,partial);int unavailable=0,solid=0;
    for(const auto face:selected.faces(*mesh))for(int half=0;half<2;++half){if(face.negativeCell<0 && face.positiveCell<0){++solid;rejects([&]{selectedCache.traceHalfFace(face.nativeRow,half);},"solid-only");continue;}
        for(const auto& z:selectedCache.traceHalfFace(face.nativeRow,half).coincidentRanges()){check(z.minus.available && !z.plus.available,"selected outside side extrapolated");++unavailable;}}
    check(unavailable>0 && solid>0,"selected availability/solid slots not exercised");
    rejects([&]{Cache::prepareNativeTrace(*mesh,selected,source);},"activity differ");
    const auto zero=selectedTraceSource(*mesh,domain,[](Vector2D){return Jet{};});rejects([&]{Cache::prepareNativeTrace(*mesh,domain,zero);},"normal Taylor coefficients vanish");
    rejects([&]{Cache::prepareNativeTrace(*mesh,domain,nullptr);},"null mapped");
    const auto scalar=selectedTraceSource(*mesh,domain,plane);const auto full=Cache::prepareNativeTrace(*mesh,domain,scalar);check(!full.traceHalfFace(0,0).usesSource(*source),"foreign scalar accepted");
    const auto label=mesh->cellIndexToLabel[0];const auto center=mesh->cells(label).centerCoord();mesh->cells(label).setCenterCoord(center+Vector2D(.001,0));
    rejects([&]{full.validate(*mesh,domain);},"NativeMeshIdentity");rejects([&]{Cache::prepareNativeTrace(*mesh,domain,scalar);},"NativeMeshIdentity");
    mesh->cells(label).setCenterCoord(center);full.validate(*mesh,domain);
    const int rounding=std::fegetround();std::fesetround(FE_DOWNWARD);bool fp=false;try{full.traceHalfFace(0,0).normalAxis();}catch(const std::exception&){fp=true;}std::fesetround(rounding);check(fp,"typed cache bypassed strict FP source gate");
    const auto ordinarySource=selectedTraceSource(*mesh,domain,[](Vector2D p){return Jet{p.x-.375,1,0,0};});
    Cache legacy(*mesh,domain,ordinarySource);rejects([&]{legacy.traceHalfFace(0,0);},"requires prepareNativeTrace");rejects([&]{legacy.traceEventIncidenceCount();},"requires prepareNativeTrace");
    auto ordinaryCache=Cache::prepareNativeTrace(*mesh,domain,ordinarySource);bool lateFailure=false;
    for(const auto face:domain.faces(*mesh))for(int half=0;half<2 && !lateFailure;++half){const auto& q=ordinaryCache.traceHalfFace(face.nativeRow,half);
        for(std::size_t k=0;k<q.events().size();++k)if(q.events()[k].root.lower==q.events()[k].root.upper){const auto oldAddress=&q;const auto count=ordinaryCache.traceEventIncidenceCount();
            rejects([&]{ordinaryCache.refineTraceEvent(face.nativeRow,half,k);},"exact point root");
            check(&ordinaryCache.traceHalfFace(face.nativeRow,half)==oldAddress && ordinaryCache.traceEventIncidenceCount()==count,"late refinement failure changed the cache");lateFailure=true;break;}}
    check(lateFailure,"late typed refinement failure not exercised");
    Quantum::Geometry::ActiveHermitePatch cubic;cubic.activeCells.assign(mesh->cellIndexNum(),1);cubic.nodes.resize(mesh->vertIndexNum());
    for(int v=0;v<mesh->vertIndexNum();++v)cubic.nodes[v]={point(*mesh,v).x-.375,1,0,0};
    const auto wrongFamily=std::make_shared<const Snapshot>(Snapshot::fromNativeLeaves(*mesh,cubic,false));
    rejects([&]{Cache::prepareNativeTrace(*mesh,domain,wrongFamily);},"native Q4 source");
    int exterior=0;
    for(bool rz:{false,true})for(int axis:{0,1})for(double line:{0.,1.}){
        auto m=makeMesh(4);auto d=domainFor(*m,rz);const Profile boundary=[=](Vector2D p){return rz && axis==0?Jet{p.x*p.x-line*line,2*p.x,0,0}:
            axis==0?Jet{p.x-line,1,0,0}:Jet{p.y-line,0,1,0};};
        const auto s=selectedTraceSource(*m,d,boundary);const auto c=Cache::prepareNativeTrace(*m,d,s);int found=0;
        for(const auto face:d.faces(*m))for(int half=0;half<2;++half)for(const auto& z:c.traceHalfFace(face.nativeRow,half).coincidentRanges()){
            check(z.minus.available==(line==1) && z.plus.available==(line==0),"physical exterior side extrapolated");
            const auto& active=line==1?z.minus:z.plus;check(active.normalOrder==(rz && axis==0 && line==0?2:1),"axis/physical normal order changed");++found;}
        check(found==8,"physical zero edge half-face coverage");exterior+=found;
    }
    auto parallelMesh=makeMesh(8,true);auto parallelDomain=domainFor(*parallelMesh,true);
    const auto parallelSource=selectedTraceSource(*parallelMesh,parallelDomain,[](Vector2D p){return Jet{p.y-.5,0,1,0};});
    const auto parallel=Cache::prepareNativeTrace(*parallelMesh,parallelDomain,parallelSource);int parallelHalves=0;
    check(parallelDomain.nativeFaceCount()>128,"parallel cache construction threshold not exercised");
    for(const auto face:parallelDomain.faces(*parallelMesh))for(int half=0;half<2;++half){check(parallel.traceHalfFace(face.nativeRow,half).usesSource(*parallelSource),"parallel typed cache source");++parallelHalves;}
    std::cout<<"NATIVE_HALF_TRACE_CACHE_EXTRA,"<<exterior<<','<<unavailable<<','<<solid<<','<<parallelHalves<<",PASS\n";
    std::cout<<"NATIVE_HALF_TRACE_CACHE,"<<layouts<<','<<halves<<','<<zeros<<','<<hanging<<','<<refinements<<','<<peak<<",PASS\n";
}
}
int main(int argc,char** argv) try
{std::cout<<std::scientific<<std::setprecision(16);if(argc==2 && std::string(argv[1])=="--trace-cache"){traceCacheContract();return 0;}
    coreRanges();translationAndPhases();counterflowAndRepair();multiRootAndConvergence();tangentFullPhase();tangentMappedRoots();failuresAndLifetime();directionalSubsetIdentity();traceCacheContract();std::cout<<"I08 unique phase-flux contract: PASS\n";return 0;}
catch(const std::exception& e){std::cerr<<"I08 unique phase-flux contract: FAIL: "<<e.what()<<'\n';return 1;}
