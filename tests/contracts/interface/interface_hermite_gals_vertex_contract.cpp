#include "Modules/Interface/HermiteGalsVertexTransport.h"

#include "Core/Geometry/ActiveQuarticHermitePatch.h"
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace
{
using namespace LAMP::Interface;
using Position=HermiteGalsTimeLayer::Position;
using Query=HermiteGalsVelocityView::Query;
using Jet=LAMP::Geometry::ActiveHermitePatch::NodalDof;
using Patch=LAMP::Geometry::ActiveQuarticHermitePatch;
using Snapshot=LevelSetHermiteGeometry::Snapshot;
using Domain=LevelSetBoundaryTrace::Domain;
using Profile=std::function<Jet(Vector2D)>;
static_assert(std::is_same_v<decltype(std::declval<HermiteGalsVertexTransportResult&>().vertices()),
    const std::vector<HermiteGalsVertexJet>&>,"vertex results must be read-only");
static_assert(!std::is_default_constructible_v<HermiteGalsVertexTransportResult>,"no fabricated empty stage result");
void check(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
void close(double a,double b,double tolerance,const char* why)
{check(std::isfinite(a) && std::abs(a-b)<=tolerance*(1+std::abs(b)),why);}
template<class F> void rejects(F action,const char* message)
{
    try {action();}
    catch(const std::exception& error) {
        const std::string text=error.what();
        check(text.find("HermiteGalsVertexTransport source=")!=std::string::npos && text.find("vertex=")!=std::string::npos,
            "vertex failure lacks stage/source/vertex context");
        if(text.find(message)==std::string::npos) throw std::runtime_error("unexpected rejection: "+text);
        return;
    }
    throw std::runtime_error(std::string("missing vertex rejection: ")+message);
}
DictionaryMesh input(int n=8)
{DictionaryMesh d;d.Nx=d.Ny=n;d.maxCellNum=4096;d.maxInitialIterations=0;return d;}
Domain domainFor(const QuadtreeMesh& mesh,bool radial=false)
{
    Vector2D h{1,1};
    for(int c=0;c<mesh.cellIndexNum();++c){const auto s=mesh.cells(mesh.cellIndexToLabel[c]).size();h.x=std::min(h.x,s.x);h.y=std::min(h.y,s.y);}
    return {mesh.minCoord,mesh.maxCoord,h,radial};
}
Jet polynomial(Vector2D p,bool radial=false)
{
    const double x=p.x,y=p.y;
    const double f=radial?1+.2*x*x+.1*x*x*x*x:1+.3*x+.2*x*x+.1*x*x*x+.05*x*x*x*x;
    const double fx=radial?.4*x+.4*x*x*x:.3+.4*x+.3*x*x+.2*x*x*x;
    const double g=.7-.2*y+.3*y*y-.1*y*y*y+.05*y*y*y*y;
    const double gy=-.2+.6*y-.3*y*y+.2*y*y*y;
    return {.1+f*g,fx*g,f*gy,fx*gy}; // strictly positive: exact M_K=0
}

// Analytic fixtures feed the existing Core native Q4 factory. They are not a
// new V2 initializer/closure or an AMR transport implementation.
HermiteGalsState sourceFor(const QuadtreeMesh& mesh,const Profile& profile,bool radial=false,
    int support=0,bool sector=false)
{
    Patch patch;patch.activeCells.assign(mesh.cellIndexNum(),0);patch.nodes.resize(mesh.vertIndexNum());
    patch.cellCenters.resize(mesh.cellIndexNum());
    HermiteGalsPrimaryState primary;
    primary.measure=radial?HermiteGalsMeasure::AxisymmetricVolume:HermiteGalsMeasure::CartesianArea;
    std::vector<bool> used(mesh.vertIndexNum(),false);
    std::map<std::array<int,2>,bool> sides;
    int sectorCell=-1,pin=-1;
    for(int c=0;c<mesh.cellIndexNum();++c){
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto center=cell.centerCoord();
        double liquid=0;
        if(sector){ // This fixture is x+y-1; its modified corner cell stays negative.
            const auto lo=mesh.verts(cell.cellToVert[0]).coord(),hi=mesh.verts(cell.cellToVert[2]).coord();
            const auto area=[](double t){return .5*std::pow(std::max(t,0.0),2);};
            liquid=area(1-lo.x-lo.y)-area(1-hi.x-lo.y)-area(1-lo.x-hi.y)+area(1-hi.x-hi.y);
        }
        primary.phaseVolumes.push_back({c,liquid});
        if((support==1 && center.x>.5) || (support==2 && center.x>.5 && center.y<.5))continue;
        patch.activeCells[c]=1;primary.geometryCells.push_back(c);patch.cellCenters[c]=profile(center).value;
        const auto hi=mesh.verts(cell.cellToVert[2]).coord();
        if(hi.x==.5 && hi.y==.5)sectorCell=c;
        for(int k=0;k<4;++k){
            const auto& a=mesh.verts(cell.cellToVert[k]);const auto& b=mesh.verts(cell.cellToVert[(k+1)%4]);
            used[a.index()]=true;
            const bool horizontal=a.coord().y==b.coord().y;
            std::array<int,2> ids{{a.index(),b.index()}};
            if(horizontal?a.coord().x>b.coord().x:a.coord().y>b.coord().y)std::swap(ids[0],ids[1]);
            if(!sides.emplace(ids,horizontal).second)continue;
            const auto j=profile(a.coord()+(b.coord()-a.coord())*.5);
            patch.edgeMidpoints.push_back({ids,j.value,horizontal?j.dy:j.dx});
        }
    }
    HermiteGalsClosure closure;
    for(int v=0;v<mesh.vertIndexNum();++v){
        const auto& vertex=mesh.verts(mesh.vertIndexToLabel[v]);const auto p=vertex.coord();
        patch.nodes[v]=profile(p);
        if(p.x==.5 && p.y==.5)pin=v;
        if(used[v] && vertex.type()==TVertType::non){const auto j=patch.nodes[v];
            primary.vertices.push_back({v,j.value,j.dx,j.dy});closure.vertexMixedDerivatives.push_back(j.dxy);}
    }
    // Only maximal sides are primary master midpoint DOFs. Fine-side payloads
    // in patch are fixture inputs to Core's existing algebraic constraints.
    for(const auto& e:patch.edgeMidpoints){
        const auto a=mesh.verts(mesh.vertIndexToLabel[e.vertices[0]]).coord();
        const auto b=mesh.verts(mesh.vertIndexToLabel[e.vertices[1]]).coord();bool contained=false;
        for(const auto& other:patch.edgeMidpoints)if(other.vertices!=e.vertices){
            const auto lo=mesh.verts(mesh.vertIndexToLabel[other.vertices[0]]).coord();
            const auto hi=mesh.verts(mesh.vertIndexToLabel[other.vertices[1]]).coord();
            if(a.y==b.y && lo.y==hi.y && a.y==lo.y && a.x>=lo.x && b.x<=hi.x)contained=true;
            if(a.x==b.x && lo.x==hi.x && a.x==lo.x && a.y>=lo.y && b.y<=hi.y)contained=true;
        }
        if(!contained)primary.edges.push_back({e.vertices,e.value,e.normalDerivative});
    }
    // Fixture boundary extension uses zero center datum, so its remaining Q4
    // cardinal center coefficient equals this analytic center sample.
    for(int c:primary.geometryCells)closure.cellBubbleCoefficients.push_back(patch.cellCenters[c]);
    if(sector){
        check(pin>=0 && sectorCell>=0,"missing sector fixture topology");
        auto jet=patch.nodes[pin];jet.dy=0;patch.cornerSectors.push_back({sectorCell,pin,jet,-1});
    }
    HermiteGalsCandidate candidate(mesh,{},primary);
    candidate.bindGeometry(mesh,candidate.identity(),closure,Snapshot::fromNativeLeaves(mesh,patch,radial));
    return std::move(candidate).finish();
}
std::vector<int> selection(const QuadtreeMesh& mesh,const HermiteGalsState& source,bool interior)
{
    std::vector<int> ids;
    for(const auto& jet:source.primary().vertices){const auto p=mesh.verts(mesh.vertIndexToLabel[jet.vertex]).coord();
        if(!interior || (p.x>0 && p.x<1 && p.y>0 && p.y<1))ids.push_back(jet.vertex);}
    return ids;
}
void unchanged(const HermiteGalsPrimaryState& a,const HermiteGalsPrimaryState& b)
{
    check(a.measure==b.measure && a.geometryCells==b.geometryCells && a.vertices.size()==b.vertices.size() &&
        a.edges.size()==b.edges.size() && a.phaseVolumes.size()==b.phaseVolumes.size(),"source layout changed");
    for(std::size_t i=0;i<a.vertices.size();++i){const auto& x=a.vertices[i];const auto& y=b.vertices[i];
        check(x.vertex==y.vertex && x.phi==y.phi && x.dx==y.dx && x.dy==y.dy,"source first jet changed");}
    for(std::size_t i=0;i<a.edges.size();++i){const auto& x=a.edges[i];const auto& y=b.edges[i];
        check(x.vertices==y.vertices && x.phi==y.phi && x.dn==y.dn,"edge changed in vertex-only stage");}
    for(std::size_t i=0;i<a.phaseVolumes.size();++i)
        check(a.phaseVolumes[i].cell==b.phaseVolumes[i].cell && a.phaseVolumes[i].liquid==b.phaseVolumes[i].liquid,"inventory changed");
}
HermiteGalsVelocityView oldView(const HermiteGalsState& source,Query query)
{return {source.identity(),source.timeLayer(),std::move(query)};}
HermiteGalsVelocityView newView(const HermiteGalsState& source,double dt,Query query)
{return {source.identity(),{1,0,dt,Position::End},std::move(query)};}

void affineAndScalarDifference()
{
    QuadtreeMesh mesh(input());const auto source=sourceFor(mesh,[](Vector2D p){return polynomial(p);});
    const auto vertices=selection(mesh,source,true),all=selection(mesh,source,false);
    const auto before=source.primary();const auto beforeClosure=source.closure();const auto id=source.identity();
    const auto view=source.geometry();const auto* pointer=&view.snapshot();const auto domain=domainFor(mesh);
    double maxPhi=0,maxGradient=0;
    for(bool shear:{false,true}){
        const double a=shear?.12:0,b=-.03,v=.04,dt=.1;
        const Query velocity=[=](const Vector2D& p){return LevelSetBoundaryTrace::VelocityJet{{a*p.y+b,v},{0,a,0,0}};};
        const auto old=oldView(source,velocity),now=newView(source,dt,velocity);
        const auto result=transportHermiteGalsVertices(mesh,source,vertices,domain,old,now);
        check(result.vertices().size()==49 && result.sourceIdentity()==id && result.arrivalVelocityGeometryIdentity()==now.geometryIdentity()
            && result.sourceTimeLayer().time()==0 && result.arrivalTimeLayer().time()==dt,"vertex result lost selection or provenance");
        for(const auto& jet:result.vertices()){
            const auto p=mesh.verts(mesh.vertIndexToLabel[jet.vertex]).coord();
            const auto exact=polynomial({p.x-dt*(a*p.y+b)+.5*dt*dt*a*v,p.y-dt*v});
            const double dx=exact.dx,dy=exact.dy-dt*a*exact.dx;
            close(jet.phi,exact.value,2e-12,"affine transported value");close(jet.dx,dx,2e-12,"affine transported dx");
            close(jet.dy,dy,2e-12,"affine transported dy / transpose convention");
            maxPhi=std::max(maxPhi,std::abs(jet.phi-exact.value));
            maxGradient=std::max({maxGradient,std::abs(jet.dx-dx),std::abs(jet.dy-dy)});
        }
    }
    const Query nonlinear=[](const Vector2D& p){return LevelSetBoundaryTrace::VelocityJet{
        {.12*p.x*p.y+.02,.04*p.x*p.x-.03*p.y},{.12*p.y,.12*p.x,.08*p.x,-.03}};};
    const auto old=oldView(source,nonlinear),now=newView(source,.08,nonlinear);
    const auto result=transportHermiteGalsVertices(mesh,source,vertices,domain,old,now);
    const auto scalar=[&](Vector2D p){const auto trace=traceHermiteGalsWithJacobian(domain,p,old,now);
        return pointer->value(trace.departurePoint());};
    double fdError=0;const double h=1e-5;
    for(const auto& jet:result.vertices()){
        const auto p=mesh.verts(mesh.vertIndexToLabel[jet.vertex]).coord();
        const double dx=(scalar({p.x+h,p.y})-scalar({p.x-h,p.y}))/(2*h);
        const double dy=(scalar({p.x,p.y+h})-scalar({p.x,p.y-h}))/(2*h);
        fdError=std::max({fdError,std::abs(jet.dx-dx),std::abs(jet.dy-dy)});
    }
    check(fdError<2e-8,"vertex gradient differs from the derivative of the actual composed scalar");
    const Query zero=[](const Vector2D&){return LevelSetBoundaryTrace::VelocityJet{{0,0},{}};};
    double stationaryValueError=0,stationaryGradientError=0;
    for(bool zeroStep:{false,true}){
        const auto stopped=oldView(source,zero);
        const auto next=zeroStep?stopped:newView(source,.1,zero);
        const auto fixed=transportHermiteGalsVertices(mesh,source,all,domain,stopped,next);
        for(std::size_t i=0;i<all.size();++i){const auto& a=fixed.vertices()[i];const auto& b=before.vertices[i];
            const auto sampled=pointer->scalarJetNative(mesh.verts(mesh.vertIndexToLabel[a.vertex]).coord());
            // Core's normalized metric cache may round a nodal derivative.
            // An identity pullback must exactly retain the authoritative sample,
            // with no extra transport error or bypass through primary data.
            check(a.phi==sampled.value && a.dx==sampled.dx && a.dy==sampled.dy,
                "stationary/zero-step pullback changed the source sample");
            const double roundoff=8*std::numeric_limits<double>::epsilon();
            close(a.phi,b.phi,roundoff,"stationary value exceeds source roundoff");
            close(a.dx,b.dx,roundoff,"stationary dx exceeds source roundoff");
            close(a.dy,b.dy,roundoff,"stationary dy exceeds source roundoff");
            stationaryValueError=std::max(stationaryValueError,std::abs(a.phi-b.phi));
            stationaryGradientError=std::max({stationaryGradientError,std::abs(a.dx-b.dx),std::abs(a.dy-b.dy)});
        }
    }
    unchanged(before,source.primary());
    check(source.identity()==id && &source.geometry().snapshot()==pointer && source.timeLayer().time()==0 &&
        beforeClosure.vertexMixedDerivatives==source.closure().vertexMixedDerivatives &&
        beforeClosure.cellBubbleCoefficients==source.closure().cellBubbleCoefficients,"source closure or authority changed");
    std::cout<<"AFFINE,"<<maxPhi<<','<<maxGradient<<'\n';
    std::cout<<"SCALAR_FD,"<<vertices.size()<<','<<h<<','<<fdError<<'\n';
    std::cout<<"STATIONARY,"<<all.size()<<','<<stationaryValueError<<','<<stationaryGradientError<<'\n';
    std::cout<<"stationary/isolation: source samples bitwise fixed; edges, inventory, closure and snapshot retained PASS\n";
}

void analyticNonlinear(bool radial)
{
    QuadtreeMesh mesh(input());const auto source=sourceFor(mesh,[=](Vector2D p){return polynomial(p,radial);},radial);
    const auto vertices=selection(mesh,source,false);const auto domain=domainFor(mesh,radial);
    const double a=.7,b=-.4;
    const Query velocity=[=](const Vector2D& p){
        return LevelSetBoundaryTrace::VelocityJet{{radial?a*p.x*(1-p.x*p.x):a*p.x*(1-p.x),b*p.y*(1-p.y)},
            {radial?a*(1-3*p.x*p.x):a*(1-2*p.x),0,0,b*(1-2*p.y)}};
    };
    double previousValue=0,previousGradient=0;
    for(double dt:{.2,.1,.05,.025}){
        const auto old=oldView(source,velocity),now=newView(source,dt,velocity);
        const auto result=transportHermiteGalsVertices(mesh,source,vertices,domain,old,now);
        double valueError=0,gradientError=0;int fixedNontrivial=0;
        for(const auto& jet:result.vertices()){
            const auto p=mesh.verts(mesh.vertIndexToLabel[jet.vertex]).coord();
            const double ex=std::exp((radial?2:1)*a*dt),ey=std::exp(b*dt);
            const double denominator=radial?p.x*p.x+(1-p.x*p.x)*ex:p.x+(1-p.x)*ex;
            const double scale=radial?std::sqrt(denominator):denominator;
            const Vector2D departure{p.x/scale,p.y/(p.y+(1-p.y)*ey)};
            const double jx=radial?ex/(denominator*scale):ex/(denominator*denominator);
            const double sy=p.y+(1-p.y)*ey,jy=ey/(sy*sy);
            const auto exact=polynomial(departure,radial);
            valueError=std::max(valueError,std::abs(jet.phi-exact.value));
            gradientError=std::max({gradientError,std::abs(jet.dx-jx*exact.dx),std::abs(jet.dy-jy*exact.dy)});
            if((p.x==0 || p.x==1) && (p.y==0 || p.y==1)){
                const auto original=polynomial(p,radial);
                check(jet.phi==original.value,"fixed corner value changed");
                if(jet.dx!=original.dx || jet.dy!=original.dy)++fixedNontrivial;
            }
            if(radial && p.x==0)check(jet.dx==0,"axis scalar lost even first-jet parity");
        }
        check(vertices.size()==81 && fixedNontrivial>0,"nonlinear test missed full skeleton or fixed-point tangent");
        const double valueOrder=previousValue>0?std::log2(previousValue/valueError):0;
        const double gradientOrder=previousGradient>0?std::log2(previousGradient/gradientError):0;
        std::cout<<"LOCAL_ORDER,"<<(radial?"RZ":"XY")<<','<<vertices.size()<<','<<dt<<','<<valueError<<','
            <<gradientError<<','<<valueOrder<<','<<gradientOrder<<'\n';
        if(dt<=.05)check(valueOrder>=2.8 && valueOrder<=3.2 && gradientOrder>=2.8 && gradientOrder<=3.2,
            "vertex first jet failed expected RK2 local third-order error");
        previousValue=valueError;previousGradient=gradientError;
    }
}

void rejectionAndAtomicity()
{
    QuadtreeMesh mesh(input());const auto profile=[](Vector2D p){return polynomial(p);};
    const auto source=sourceFor(mesh,profile);const auto vertices=selection(mesh,source,true);
    const auto all=selection(mesh,source,false);const auto domain=domainFor(mesh);const auto before=source.primary();
    const auto beforeClosure=source.closure();const auto beforeView=source.geometry();
    const Query zero=[](const Vector2D&){return LevelSetBoundaryTrace::VelocityJet{{0,0},{}};};
    const auto old=oldView(source,zero),now=newView(source,.1,zero);
    const auto run=[&](const std::vector<int>& ids){return transportHermiteGalsVertices(mesh,source,ids,domain,old,now);};
    rejects([&]{run({});},"nonempty");rejects([&]{run({-1});},"ascending, unique and native");
    rejects([&]{run({vertices[0],vertices[0]});},"ascending, unique and native");
    rejects([&]{run({mesh.vertIndexNum()});},"ascending, unique and native");
    const auto foreign=sourceFor(mesh,profile);const auto wrongOld=oldView(foreign,zero);
    rejects([&]{transportHermiteGalsVertices(mesh,source,vertices,domain,wrongOld,now);},"geometry differs");
    const HermiteGalsVelocityView wrongLayer(source.identity(),{1,0,0,Position::Initial},zero);
    rejects([&]{transportHermiteGalsVertices(mesh,source,vertices,domain,wrongLayer,now);},"time layer differs");
    const HermiteGalsVelocityView wrongBegin(source.identity(),{1,.01,.1,Position::End},zero);
    rejects([&]{transportHermiteGalsVertices(mesh,source,vertices,domain,old,wrongBegin);},"begin at source time");
    auto wrongDomain=domain;wrongDomain.maxCoord.x=2;
    rejects([&]{transportHermiteGalsVertices(mesh,source,vertices,wrongDomain,old,now);},"bounds");
    wrongDomain=domain;wrongDomain.radialAxisSymmetry=true;
    rejects([&]{transportHermiteGalsVertices(mesh,source,vertices,wrongDomain,old,now);},"axis semantics");
    wrongDomain=domain;wrongDomain.cellSize={1,1};
    rejects([&]{transportHermiteGalsVertices(mesh,source,vertices,wrongDomain,old,now);},"minimum native leaf size");
    auto changed=mesh.clone();changed->verts(changed->vertIndexToLabel[0]).setCoord({-.1,0});
    rejects([&]{transportHermiteGalsVertices(*changed,source,vertices,domain,old,now);},"mesh geometry/topology changed");
    const Query inflow=[](const Vector2D&){return LevelSetBoundaryTrace::VelocityJet{{.2,0},{}};};
    const auto inOld=oldView(source,inflow),inNow=newView(source,.1,inflow);
    rejects([&]{transportHermiteGalsVertices(mesh,source,all,domain,inOld,inNow);},"boundary-history tangent");
    int queries=0;
    const Query fails=[&](const Vector2D&){if(++queries==5)throw std::runtime_error("injected provider failure");
        return LevelSetBoundaryTrace::VelocityJet{{0,0},{}};};
    const auto failOld=oldView(source,fails),failNow=newView(source,.1,fails);
    rejects([&]{transportHermiteGalsVertices(mesh,source,vertices,domain,failOld,failNow);},"injected provider failure");
    check(queries==5,"failure did not occur after a successful first vertex");
    unchanged(before,source.primary());
    check(source.identity()==beforeView.identity() && &source.geometry().snapshot()==&beforeView.snapshot() &&
        source.closure().vertexMixedDerivatives==beforeClosure.vertexMixedDerivatives &&
        source.closure().cellBubbleCoefficients==beforeClosure.cellBubbleCoefficients,"failed batch polluted source");

    const auto band=sourceFor(mesh,profile,false,1);const auto bandOld=oldView(band,inflow);
    int boundary=-1,missing=-1;
    for(const auto& jet:source.primary().vertices){const auto p=mesh.verts(mesh.vertIndexToLabel[jet.vertex]).coord();
        if(p.x==.5 && p.y==.5)boundary=jet.vertex;if(p.x==.75 && p.y==.5)missing=jet.vertex;}
    const auto bandNow=newView(band,.1,inflow);
    rejects([&]{transportHermiteGalsVertices(mesh,band,{missing},domain,bandOld,bandNow);},"not a source master");
    const Query escape=[](const Vector2D&){return LevelSetBoundaryTrace::VelocityJet{{-.2,0},{}};};
    const auto escapeOld=oldView(band,escape),escapeNew=newView(band,.1,escape);
    rejects([&]{transportHermiteGalsVertices(mesh,band,{boundary},domain,escapeOld,escapeNew);},"outside the active domain");

    const auto giant=sourceFor(mesh,[](Vector2D p){return Jet{1e200*(1+p.x),1e200,0,0};});
    const Query steep=[](const Vector2D& p){return LevelSetBoundaryTrace::VelocityJet{{1e200*(p.x-.5),0},{1e200,0,0,0}};};
    const auto steepOld=oldView(giant,steep),steepNew=newView(giant,.1,zero);
    rejects([&]{transportHermiteGalsVertices(mesh,giant,{boundary},domain,steepOld,steepNew);},"nonfinite transported gradient");

    QuadtreeMesh amr(input(4));DictionaryMesh::RefineBox box;box.min={.25,.25};box.max={.5,.5};box.level=1;
    amr.applyRefineBoxes({box});const auto adaptive=sourceFor(amr,profile);int hanging=-1;
    for(int v=0;v<amr.vertIndexNum();++v)if(amr.verts(amr.vertIndexToLabel[v]).type()!=TVertType::non){hanging=v;break;}
    check(hanging>=0,"AMR fixture lacks a hanging vertex");
    const auto adaptiveOld=oldView(adaptive,zero),adaptiveNew=newView(adaptive,.1,zero);
    rejects([&]{transportHermiteGalsVertices(amr,adaptive,{hanging},domainFor(amr),adaptiveOld,adaptiveNew);},"hanging vertex");

    // A real Core one-sided sector is not a unique shared vertex first jet.
    const auto sector=sourceFor(mesh,[](Vector2D p){return Jet{p.x+p.y-1,1,1,0};},false,2,true);
    const auto sectorOld=oldView(sector,zero),sectorNew=newView(sector,.1,zero);
    rejects([&]{transportHermiteGalsVertices(mesh,sector,{boundary},domain,sectorOld,sectorNew);},"one-sided vertex transport");
    std::cout<<"failure isolation: provenance, topology, geometry coverage, boundary, overflow and mid-batch rejection PASS\n";
}
}

int main()
{
    try{std::cout<<std::setprecision(17);affineAndScalarDifference();analyticNonlinear(false);analyticNonlinear(true);rejectionAndAtomicity();}
    catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
    std::cout<<"I03 vertex first-jet contract PASS (no edge/closure/conservative geometry advancement)\n";
    return 0;
}
