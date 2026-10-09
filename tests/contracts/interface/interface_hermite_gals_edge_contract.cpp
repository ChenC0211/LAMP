#include "Modules/Interface/HermiteGalsEdgeTransport.h"

#include "Core/Geometry/ActiveQuarticHermitePatch.h"
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace
{
using namespace Quantum::Interface;
using Key=std::array<int,2>;
using Jet=Quantum::Geometry::ActiveHermitePatch::NodalDof;
using Patch=Quantum::Geometry::ActiveQuarticHermitePatch;
using Snapshot=LevelSetHermiteGeometry::Snapshot;
using Query=HermiteGalsVelocityView::Query;
using Position=HermiteGalsTimeLayer::Position;
static_assert(std::is_same_v<decltype(std::declval<HermiteGalsEdgeTransportResult&>().edges()),
    const std::vector<HermiteGalsEdgeJet>&>,"edge results must be read-only");
static_assert(!std::is_default_constructible_v<HermiteGalsEdgeTransportResult>,"no unbound stage result");
void check(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
void close(double a,double b,double tolerance,const char* why)
{check(std::isfinite(a) && std::abs(a-b)<=tolerance*(1+std::abs(b)),why);}
template<class F> void rejects(F action,const char* message)
{
    try{action();}catch(const std::exception& error){
        const std::string s=error.what();
        check(s.find("HermiteGalsEdgeTransport source=")!=std::string::npos && s.find("arrival=")!=std::string::npos &&
            s.find("edge=(")!=std::string::npos,"edge failure lacks identities/endpoints");
        if(s.find(message)==std::string::npos)throw std::runtime_error("unexpected rejection: "+s);
        return;
    }
    throw std::runtime_error(std::string("missing edge rejection: ")+message);
}
DictionaryMesh input(int n=8)
{DictionaryMesh d;d.Nx=d.Ny=n;d.maxCellNum=4096;d.maxInitialIterations=0;return d;}
Vector2D point(const QuadtreeMesh& mesh,int v){return mesh.verts(mesh.vertIndexToLabel[v]).coord();}
Vector2D midpoint(const QuadtreeMesh& mesh,Key key)
{const auto a=point(mesh,key[0]),b=point(mesh,key[1]);return a+(b-a)*.5;}
bool horizontal(const QuadtreeMesh& mesh,Key key){return point(mesh,key[0]).y==point(mesh,key[1]).y;}
LevelSetBoundaryTrace::Domain domainFor(const QuadtreeMesh& mesh,bool radial=false)
{
    Vector2D size=mesh.maxCoord-mesh.minCoord;
    for(int c=0;c<mesh.cellIndexNum();++c){const auto h=mesh.cells(mesh.cellIndexToLabel[c]).size();
        size.x=std::min(size.x,h.x);size.y=std::min(size.y,h.y);}
    return {mesh.minCoord,mesh.maxCoord,size,radial};
}
Jet polynomial(Vector2D p,bool radial=false)
{
    const double x=p.x,y=p.y;
    const double f=radial?1+.2*x*x+.1*x*x*x*x:1+.3*x+.2*x*x+.1*x*x*x+.05*x*x*x*x;
    const double fx=radial?.4*x+.4*x*x*x:.3+.4*x+.3*x*x+.2*x*x*x;
    const double g=.7-.2*y+.3*y*y-.1*y*y*y+.05*y*y*y*y;
    const double gy=-.2+.6*y-.3*y*y+.2*y*y*y;
    return {.1+f*g,fx*g,f*gy,fx*gy}; // positive Q4, exact negative-phase M=0
}
struct Options {
    bool radial=false,sector=false,includeSubedges=false;
    int support=0;
    std::vector<Key> extra;
    Key omit{{-1,-1}};
};
// Analytic test fixture only. Full Q4 factory input and independently filtered
// maximal sides; this does not implement V2 initialization/AMR closure.
HermiteGalsState sourceFor(const QuadtreeMesh& mesh,const Options& options={})
{
    const auto profile=[&](Vector2D p){return options.sector?Jet{p.x+p.y-1,1,1,0}:polynomial(p,options.radial);};
    Patch patch;patch.activeCells.assign(mesh.cellIndexNum(),0);patch.nodes.resize(mesh.vertIndexNum());
    patch.cellCenters.resize(mesh.cellIndexNum());HermiteGalsPrimaryState primary;HermiteGalsClosure closure;
    primary.measure=options.radial?HermiteGalsMeasure::AxisymmetricVolume:HermiteGalsMeasure::CartesianArea;
    std::vector<bool> used(mesh.vertIndexNum(),false);std::set<Key> seen;
    int pin=-1,sectorCell=-1;
    for(int c=0;c<mesh.cellIndexNum();++c){
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto center=cell.centerCoord();
        const auto lo=mesh.verts(cell.cellToVert[0]).coord(),hi=mesh.verts(cell.cellToVert[2]).coord();
        double liquid=0;
        if(options.sector){const auto a=[](double s){return .5*std::pow(std::max(s,0.0),2);};
            liquid=a(1-lo.x-lo.y)-a(1-hi.x-lo.y)-a(1-lo.x-hi.y)+a(1-hi.x-hi.y);}
        primary.phaseVolumes.push_back({c,liquid});
        if((options.support==1 && center.x>.5) || (options.sector && center.x>.5 && center.y<.5))continue;
        patch.activeCells[c]=1;primary.geometryCells.push_back(c);patch.cellCenters[c]=profile(center).value;
        closure.cellBubbleCoefficients.push_back(patch.cellCenters[c]); // zero-center boundary extension
        if(hi.x==.5 && hi.y==.5)sectorCell=c;
        for(int k=0;k<4;++k){
            const auto& a=mesh.verts(cell.cellToVert[k]);const auto& b=mesh.verts(cell.cellToVert[(k+1)%4]);
            used[a.index()]=true;Key key{{a.index(),b.index()}};const bool h=a.coord().y==b.coord().y;
            if(h?a.coord().x>b.coord().x:a.coord().y>b.coord().y)std::swap(key[0],key[1]);
            if(seen.insert(key).second){const auto j=profile(midpoint(mesh,key));
                patch.edgeMidpoints.push_back({key,j.value,h?j.dy:j.dx});}
        }
    }
    for(int v=0;v<mesh.vertIndexNum();++v){const auto p=point(mesh,v);patch.nodes[v]=profile(p);
        if(p.x==.5 && p.y==.5)pin=v;
        if(used[v] && mesh.verts(mesh.vertIndexToLabel[v]).type()==TVertType::non){const auto j=patch.nodes[v];
            primary.vertices.push_back({v,j.value,j.dx,j.dy});closure.vertexMixedDerivatives.push_back(j.dxy);}}
    for(const auto& e:patch.edgeMidpoints){bool contained=false;const auto a=point(mesh,e.vertices[0]),b=point(mesh,e.vertices[1]);
        for(const auto& other:patch.edgeMidpoints)if(e.vertices!=other.vertices){
            const auto lo=point(mesh,other.vertices[0]),hi=point(mesh,other.vertices[1]);
            if(a.y==b.y && lo.y==hi.y && a.y==lo.y && a.x>=lo.x && b.x<=hi.x)contained=true;
            if(a.x==b.x && lo.x==hi.x && a.x==lo.x && a.y>=lo.y && b.y<=hi.y)contained=true;
        }
        if((!contained || options.includeSubedges) && e.vertices!=options.omit)
            primary.edges.push_back({e.vertices,e.value,e.normalDerivative});
    }
    for(const auto& key:options.extra){const auto j=profile(midpoint(mesh,key));
        primary.edges.push_back({key,j.value,horizontal(mesh,key)?j.dy:j.dx});}
    if(options.sector){check(pin>=0 && sectorCell>=0,"missing one-sided fixture topology");
        auto j=patch.nodes[pin];j.dy=0;patch.cornerSectors.push_back({sectorCell,pin,j,-1});}
    HermiteGalsCandidate candidate(mesh,{},primary);
    candidate.bindGeometry(mesh,candidate.identity(),closure,Snapshot::fromNativeLeaves(mesh,patch,options.radial));
    return std::move(candidate).finish();
}
std::vector<Key> selection(const QuadtreeMesh& mesh,const HermiteGalsState& source,bool interior=false)
{
    std::vector<Key> keys;for(const auto& e:source.primary().edges){const auto p=midpoint(mesh,e.vertices);
        if(!interior || (p.x>0 && p.x<1 && p.y>0 && p.y<1))keys.push_back(e.vertices);}
    std::sort(keys.begin(),keys.end());return keys;
}
HermiteGalsVelocityView oldView(const HermiteGalsState& source,Query q)
{return {source.identity(),source.timeLayer(),std::move(q)};}
HermiteGalsVelocityView newView(const HermiteGalsState& source,double dt,Query q)
{return {source.identity(),{1,0,dt,Position::End},std::move(q)};}
const Query zero=[](const Vector2D&){return LevelSetBoundaryTrace::VelocityJet{{0,0},{}};};
void unchanged(const HermiteGalsPrimaryState& a,const HermiteGalsPrimaryState& b)
{
    check(a.measure==b.measure && a.geometryCells==b.geometryCells && a.vertices.size()==b.vertices.size() &&
        a.edges.size()==b.edges.size() && a.phaseVolumes.size()==b.phaseVolumes.size(),"source layout changed");
    for(std::size_t k=0;k<a.vertices.size();++k){const auto& x=a.vertices[k];const auto& y=b.vertices[k];
        check(x.vertex==y.vertex && x.phi==y.phi && x.dx==y.dx && x.dy==y.dy,"source vertices changed");}
    for(std::size_t k=0;k<a.edges.size();++k){const auto& x=a.edges[k];const auto& y=b.edges[k];
        check(x.vertices==y.vertices && x.phi==y.phi && x.dn==y.dn,"source edges changed");}
    for(std::size_t k=0;k<a.phaseVolumes.size();++k)
        check(a.phaseVolumes[k].cell==b.phaseVolumes[k].cell && a.phaseVolumes[k].liquid==b.phaseVolumes[k].liquid,"source M changed");
}

void affineAndDifference()
{
    QuadtreeMesh mesh(input());const auto source=sourceFor(mesh);const auto domain=domainFor(mesh);
    const auto keys=selection(mesh,source,true),all=selection(mesh,source);const auto before=source.primary();
    const auto beforeClosure=source.closure();const auto view=source.geometry();const auto* scalar=&view.snapshot();
    int nh=0,nv=0;for(const auto& key:keys)(horizontal(mesh,key)?nh:nv)++;
    check(nh==56 && nv==56,"affine test must cover every interior edge in both directions");
    double valueError=0,normalError=0;
    for(int kind:{0,1,2}){
        const double a=kind==0?0:.12,b=-.03,v=.04,dt=.1;
        const Query velocity=[=](const Vector2D& p){return kind==2?LevelSetBoundaryTrace::VelocityJet{{b,a*p.x+v},{0,0,a,0}}:
            LevelSetBoundaryTrace::VelocityJet{{a*p.y+b,v},{0,a,0,0}};};
        const auto old=oldView(source,velocity),now=newView(source,dt,velocity);
        const auto result=transportHermiteGalsEdges(mesh,source,keys,domain,old,now);
        check(result.edges().size()==keys.size() && result.sourceIdentity()==source.identity() &&
            result.arrivalVelocityGeometryIdentity()==now.geometryIdentity() && result.sourceTimeLayer().time()==0 &&
            result.arrivalTimeLayer().time()==dt,"edge selection/provenance changed");
        for(const auto& e:result.edges()){
            const auto p=midpoint(mesh,e.vertices);
            const Vector2D departure=kind==2?Vector2D{p.x-dt*b,p.y-dt*(a*p.x+v)+.5*dt*dt*a*b}:
                Vector2D{p.x-dt*(a*p.y+b)+.5*dt*dt*a*v,p.y-dt*v};
            const auto exact=polynomial(departure);
            const double dx=kind==2?exact.dx-dt*a*exact.dy:exact.dx;
            const double dy=kind==2?exact.dy:exact.dy-dt*a*exact.dx;
            const double dn=horizontal(mesh,e.vertices)?dy:dx;
            close(e.phi,exact.value,2e-12,"affine edge value");close(e.dn,dn,2e-12,"affine edge normal/transpose");
            valueError=std::max(valueError,std::abs(e.phi-exact.value));normalError=std::max(normalError,std::abs(e.dn-dn));
        }
    }
    const Query nonlinear=[](const Vector2D& p){return LevelSetBoundaryTrace::VelocityJet{
        {.12*p.x*p.y+.02,.04*p.x*p.x-.03*p.y},{.12*p.y,.12*p.x,.08*p.x,-.03}};};
    const auto old=oldView(source,nonlinear),now=newView(source,.08,nonlinear);
    const auto result=transportHermiteGalsEdges(mesh,source,keys,domain,old,now);
    const auto composed=[&](Vector2D p){return scalar->value(traceHermiteGalsWithJacobian(domain,p,old,now).departurePoint());};
    const double h=1e-5;double fdError=0;
    for(const auto& e:result.edges()){const auto p=midpoint(mesh,e.vertices);const Vector2D n=horizontal(mesh,e.vertices)?Vector2D{0,1}:Vector2D{1,0};
        fdError=std::max(fdError,std::abs(e.dn-(composed(p+n*h)-composed(p-n*h))/(2*h)));}
    check(fdError<2e-8,"normal derivative differs from actual composed scalar");
    double stationaryValue=0,stationaryNormal=0;
    for(bool zeroStep:{false,true}){
        const auto stopped=oldView(source,zero),next=zeroStep?stopped:newView(source,.1,zero);
        const auto fixed=transportHermiteGalsEdges(mesh,source,all,domain,stopped,next);
        for(const auto& e:fixed.edges()){
            const auto j=scalar->scalarJetNative(midpoint(mesh,e.vertices));const double dn=horizontal(mesh,e.vertices)?j.dy:j.dx;
            check(e.phi==j.value && e.dn==dn,"stationary edge pullback changed the authoritative sample");
            const auto it=std::find_if(before.edges.begin(),before.edges.end(),[&](const auto& b){return b.vertices==e.vertices;});
            close(e.phi,it->phi,8*std::numeric_limits<double>::epsilon(),"stationary edge value beyond source roundoff");
            close(e.dn,it->dn,8*std::numeric_limits<double>::epsilon(),"stationary edge normal beyond source roundoff");
            stationaryValue=std::max(stationaryValue,std::abs(e.phi-it->phi));stationaryNormal=std::max(stationaryNormal,std::abs(e.dn-it->dn));
        }
    }
    for(bool horizontalTarget:{false,true}){
        const auto key=*std::find_if(keys.begin(),keys.end(),[&](Key k){return horizontal(mesh,k)==horizontalTarget;});
        const auto center=midpoint(mesh,key);const double a=.7,b=-.4,dt=.1;
        const Query stationaryPoint=[=](const Vector2D& p){return LevelSetBoundaryTrace::VelocityJet{
            {a*(p.x-center.x),b*(p.y-center.y)},{a,0,0,b}};};
        const auto previous=oldView(source,stationaryPoint),next=newView(source,dt,stationaryPoint);
        const auto fixed=transportHermiteGalsEdges(mesh,source,{key},domain,previous,next).edges()[0];
        const auto j=scalar->scalarJetNative(center);const double rate=horizontalTarget?b:a,original=horizontalTarget?j.dy:j.dx;
        check(fixed.phi==j.value && fixed.dn!=original,"fixed midpoint incorrectly implies identity J");
        close(fixed.dn,(1-dt*rate+.5*dt*dt*rate*rate)*original,2e-14,"fixed midpoint tangent");
        const HermiteGalsVelocityView mid(source.identity(),{1,0,2*dt,Position::Midpoint},stationaryPoint);
        const auto half=transportHermiteGalsEdges(mesh,source,{key},domain,previous,mid);
        check(half.edges()[0].phi==fixed.phi && half.edges()[0].dn==fixed.dn &&
            half.arrivalTimeLayer().position==Position::Midpoint,"midpoint physical time differs from end stage");
    }
    unchanged(before,source.primary());check(source.identity()==view.identity() && &source.geometry().snapshot()==scalar &&
        source.closure().vertexMixedDerivatives==beforeClosure.vertexMixedDerivatives &&
        source.closure().cellBubbleCoefficients==beforeClosure.cellBubbleCoefficients,"source closure/snapshot changed");
    std::cout<<"EDGE_AFFINE,"<<nh<<','<<nv<<','<<valueError<<','<<normalError<<'\n';
    std::cout<<"EDGE_SCALAR_FD,"<<keys.size()<<','<<h<<','<<fdError<<'\n';
    std::cout<<"EDGE_STATIONARY,"<<all.size()<<','<<stationaryValue<<','<<stationaryNormal<<'\n';
}

void nonlinearOrder(bool radial)
{
    QuadtreeMesh mesh(input());Options options;options.radial=radial;const auto source=sourceFor(mesh,options);
    const auto keys=selection(mesh,source);const auto domain=domainFor(mesh,radial);const double a=.7,b=-.4;
    const Query velocity=[=](const Vector2D& p){return LevelSetBoundaryTrace::VelocityJet{
        {radial?a*p.x*(1-p.x*p.x):a*p.x*(1-p.x),b*p.y*(1-p.y)},
        {radial?a*(1-3*p.x*p.x):a*(1-2*p.x),0,0,b*(1-2*p.y)}};};
    double previousValue=0,previousNormal=0;
    for(double dt:{.2,.1,.05,.025}){
        const auto old=oldView(source,velocity),now=newView(source,dt,velocity);
        const auto result=transportHermiteGalsEdges(mesh,source,keys,domain,old,now);
        double valueError=0,normalError=0;int axis=0;
        for(const auto& e:result.edges()){
            const auto p=midpoint(mesh,e.vertices);const double ex=std::exp((radial?2:1)*a*dt),ey=std::exp(b*dt);
            const double sx=radial?p.x*p.x+(1-p.x*p.x)*ex:p.x+(1-p.x)*ex,sy=p.y+(1-p.y)*ey;
            const auto exact=polynomial({p.x/(radial?std::sqrt(sx):sx),p.y/sy},radial);
            const double jx=radial?ex/(sx*std::sqrt(sx)):ex/(sx*sx),jy=ey/(sy*sy);
            const double dn=horizontal(mesh,e.vertices)?jy*exact.dy:jx*exact.dx;
            valueError=std::max(valueError,std::abs(e.phi-exact.value));normalError=std::max(normalError,std::abs(e.dn-dn));
            if(radial && p.x==0){check(e.dn==0,"axis edge lost scalar parity");++axis;}
        }
        check(keys.size()==144 && (!radial || axis==8),"nonlinear test missed boundary/axis edges");
        const double po=previousValue>0?std::log2(previousValue/valueError):0,no=previousNormal>0?std::log2(previousNormal/normalError):0;
        std::cout<<"EDGE_LOCAL_ORDER,"<<(radial?"RZ":"XY")<<','<<keys.size()<<','<<dt<<','<<valueError<<','<<normalError<<','<<po<<','<<no<<'\n';
        if(dt<=.05)check(po>=2.8 && po<=3.2 && no>=2.8 && no<=3.2,"edge RK2 local order");
        previousValue=valueError;previousNormal=normalError;
    }
}

void failuresAndNativeSides()
{
    QuadtreeMesh mesh(input());const auto source=sourceFor(mesh),foreign=sourceFor(mesh);
    const auto keys=selection(mesh,source,true),all=selection(mesh,source);const auto domain=domainFor(mesh);
    const auto before=source.primary();const auto closure=source.closure();const auto view=source.geometry();
    const auto old=oldView(source,zero),now=newView(source,.1,zero);
    const auto run=[&](const std::vector<Key>& k){return transportHermiteGalsEdges(mesh,source,k,domain,old,now);};
    rejects([&]{run({});},"nonempty");rejects([&]{run({Key{{-1,0}}});},"native vertex indices");
    rejects([&]{run({Key{{0,mesh.vertIndexNum()}}});},"native vertex indices");
    rejects([&]{run({keys[0],keys[0]});},"ascending and unique");rejects([&]{run({keys[1],keys[0]});},"ascending and unique");
    int validationQueries=0;
    const Query counted=[&](const Vector2D&){++validationQueries;return LevelSetBoundaryTrace::VelocityJet{{0,0},{}};};
    const auto countedOld=oldView(source,counted),countedNew=newView(source,.1,counted);
    rejects([&]{transportHermiteGalsEdges(mesh,source,{keys[0],keys[0]},domain,countedOld,countedNew);},"ascending and unique");
    check(validationQueries==0,"selection validation must finish before querying velocity");
    rejects([&]{run({Key{{keys[0][1],keys[0][0]}}});},"positive Cartesian tangent");
    const auto wrongOld=oldView(foreign,zero);
    rejects([&]{transportHermiteGalsEdges(mesh,source,keys,domain,wrongOld,now);},"geometry differs");
    const HermiteGalsVelocityView wrongTime(source.identity(),{1,0,0,Position::Initial},zero);
    rejects([&]{transportHermiteGalsEdges(mesh,source,keys,domain,wrongTime,now);},"time layer differs");
    const HermiteGalsVelocityView wrongBegin(source.identity(),{1,.01,.1,Position::End},zero);
    rejects([&]{transportHermiteGalsEdges(mesh,source,keys,domain,old,wrongBegin);},"begin at source time");
    auto wrongDomain=domain;wrongDomain.maxCoord.x=2;
    rejects([&]{transportHermiteGalsEdges(mesh,source,keys,wrongDomain,old,now);},"bounds");
    wrongDomain=domain;wrongDomain.radialAxisSymmetry=true;
    rejects([&]{transportHermiteGalsEdges(mesh,source,keys,wrongDomain,old,now);},"axis semantics");
    wrongDomain=domain;wrongDomain.cellSize={1,1};
    rejects([&]{transportHermiteGalsEdges(mesh,source,keys,wrongDomain,old,now);},"minimum native leaf size");
    auto changed=mesh.clone();changed->verts(changed->vertIndexToLabel[0]).setCoord({-.1,0});
    rejects([&]{transportHermiteGalsEdges(*changed,source,keys,domain,old,now);},"mesh geometry/topology changed");
    Options missingOptions;missingOptions.omit=keys[0];const auto missing=sourceFor(mesh,missingOptions);
    const auto missingOld=oldView(missing,zero),missingNew=newView(missing,.1,zero);
    rejects([&]{transportHermiteGalsEdges(mesh,missing,{keys[0]},domain,missingOld,missingNew);},"not a source master edge");
    int left=-1,right=-1;
    for(int v=0;v<mesh.vertIndexNum();++v){const auto p=point(mesh,v);if(p.x==0 && p.y==.5)left=v;if(p.x==.25 && p.y==.5)right=v;}
    Options fakeOptions;fakeOptions.extra.push_back({left,right});const auto fake=sourceFor(mesh,fakeOptions);
    const auto fakeOld=oldView(fake,zero),fakeNew=newView(fake,.1,zero);
    rejects([&]{transportHermiteGalsEdges(mesh,fake,fakeOptions.extra,domain,fakeOld,fakeNew);},"not an independent native master side");
    const Query inflow=[](const Vector2D&){return LevelSetBoundaryTrace::VelocityJet{{.2,0},{}};};
    const auto inOld=oldView(source,inflow),inNow=newView(source,.1,inflow);
    rejects([&]{transportHermiteGalsEdges(mesh,source,all,domain,inOld,inNow);},"boundary-history tangent");
    int queries=0;const Query fails=[&](const Vector2D&){if(++queries==5)throw std::runtime_error("injected mid-batch failure");
        return LevelSetBoundaryTrace::VelocityJet{{0,0},{}};};
    const auto failOld=oldView(source,fails),failNew=newView(source,.1,fails);
    rejects([&]{transportHermiteGalsEdges(mesh,source,keys,domain,failOld,failNew);},"injected mid-batch failure");
    check(queries==5,"failure must follow a successful edge");
    unchanged(before,source.primary());check(source.identity()==view.identity() && &source.geometry().snapshot()==&view.snapshot() &&
        closure.vertexMixedDerivatives==source.closure().vertexMixedDerivatives && closure.cellBubbleCoefficients==source.closure().cellBubbleCoefficients,
        "failed edge batch polluted source");
    Options bandOptions;bandOptions.support=1;const auto band=sourceFor(mesh,bandOptions);
    const auto bandKeys=selection(mesh,band);const auto boundary=*std::find_if(bandKeys.begin(),bandKeys.end(),[&](Key k){return midpoint(mesh,k).x==.5;});
    const Query escape=[](const Vector2D&){return LevelSetBoundaryTrace::VelocityJet{{-.2,0},{}};};
    const auto escapeOld=oldView(band,escape),escapeNew=newView(band,.1,escape);
    rejects([&]{transportHermiteGalsEdges(mesh,band,{boundary},domain,escapeOld,escapeNew);},"outside the active domain");
    Options sectorOptions;sectorOptions.sector=true;const auto sector=sourceFor(mesh,sectorOptions);
    const auto sectorOld=oldView(sector,zero),sectorNew=newView(sector,.1,zero);
    rejects([&]{transportHermiteGalsEdges(mesh,sector,selection(mesh,sector),domain,sectorOld,sectorNew);},"one-sided edge transport");

    QuadtreeMesh amr(input(4));DictionaryMesh::RefineBox box;box.min={.25,.25};box.max={.5,.5};box.level=1;amr.applyRefineBoxes({box});
    const auto adaptive=sourceFor(amr);const auto masterKeys=selection(amr,adaptive);const auto amrDomain=domainFor(amr);
    const auto amrOld=oldView(adaptive,zero),amrNew=newView(adaptive,.1,zero);
    const auto advanced=transportHermiteGalsEdges(amr,adaptive,masterKeys,amrDomain,amrOld,amrNew);
    int hangingEnds=0;for(const auto& e:advanced.edges()){
        if(amr.verts(amr.vertIndexToLabel[e.vertices[0]]).type()!=TVertType::non ||
            amr.verts(amr.vertIndexToLabel[e.vertices[1]]).type()!=TVertType::non)++hangingEnds;
        const auto j=adaptive.geometry().snapshot().scalarJetNative(midpoint(amr,e.vertices));
        check(e.phi==j.value && e.dn==(horizontal(amr,e.vertices)?j.dy:j.dx),"AMR master stationary sample");
    }
    check(hangingEnds>0,"missing perpendicular master edges ending at hanging vertices");
    Options subedgeOptions;subedgeOptions.includeSubedges=true;const auto subedges=sourceFor(amr,subedgeOptions);
    const auto rawKeys=selection(amr,subedges);int rejected=0;
    const auto subOld=oldView(subedges,zero),subNew=newView(subedges,.1,zero);
    for(const auto& key:rawKeys)if(!std::binary_search(masterKeys.begin(),masterKeys.end(),key)){
        rejects([&]{transportHermiteGalsEdges(amr,subedges,{key},amrDomain,subOld,subNew);},"constrained subedge");++rejected;
    }
    check(rejected>0,"AMR fixture lacks constrained subedges");
    std::cout<<"EDGE_AMR,"<<masterKeys.size()<<','<<hangingEnds<<','<<rejected<<'\n';

    auto tinyInput=input(1);tinyInput.minCoord={1,0};tinyInput.maxCoord={std::nextafter(1.0,2.0),1};
    // Core rejects this grid before a state can exist; do not fabricate an
    // invalid snapshot merely to reach the edge kernel's defensive check.
    bool rejectedMesh=false;
    try{QuadtreeMesh tiny(tinyInput);}catch(const std::exception& error){
        if(std::string(error.what()).find("adjacent logical coordinates collapse in binary64")==std::string::npos)throw;
        rejectedMesh=true;
    }
    check(rejectedMesh,"unrepresentable mesh must be rejected before source construction");
    std::cout<<"EDGE_UNRESOLVED_REJECT,Core mesh construction\n";
    std::cout<<"edge rejection/source isolation/native master classification PASS\n";
}
}
int main()
{
    try{std::cout<<std::setprecision(17);affineAndDifference();nonlinearOrder(false);nonlinearOrder(true);failuresAndNativeSides();}
    catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
    std::cout<<"I04 edge first-jet contract PASS (no closure, mass update or geometry publication)\n";return 0;
}
