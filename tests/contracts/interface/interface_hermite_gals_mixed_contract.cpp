#include "Modules/Interface/HermiteGalsMixedClosure.h"

#include "Modules/Interface/HermiteGalsVertexTransport.h"
#include "Modules/Interface/HermiteGalsEdgeTransport.h"
#include "Core/Geometry/CompactQuarticDerivative.h"
#include "Core/Geometry/ActiveQuarticHermitePatch.h"
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/UniformLeafGrid.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <type_traits>

void interfaceGalsContract();
namespace
{
using namespace Quantum::Interface;
using Key=std::array<int,2>;
using Jet=Quantum::Geometry::ActiveHermitePatch::NodalDof;
using Patch=Quantum::Geometry::ActiveQuarticHermitePatch;
using Profile=std::function<Jet(Vector2D)>;
using Position=HermiteGalsTimeLayer::Position;
static_assert(std::is_same_v<decltype(std::declval<HermiteGalsMixedClosureResult&>().vertexMixedDerivatives()),
    const std::vector<double>&>,"mixed derivative result must be read-only");
static_assert(!std::is_default_constructible_v<HermiteGalsMixedClosureResult>,"no fabricated unbound closure");
void check(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
void close(double a,double b,double tolerance,const char* why)
{check(std::isfinite(a) && std::abs(a-b)<=tolerance*(1+std::abs(b)),why);}
template<class F> void rejects(F action,const char* fragment,bool module=true)
{
    try{action();}catch(const std::exception& error){const std::string s=error.what();
        if(module)check(s.find("HermiteGalsMixedClosure primary=")!=std::string::npos && s.find("vertex=")!=std::string::npos,
            "mixed closure rejection lacks provenance");
        if(s.find(fragment)==std::string::npos)throw std::runtime_error("unexpected rejection: "+s);
        return;
    }
    throw std::runtime_error(std::string("missing closure rejection: ")+fragment);
}
DictionaryMesh input(int nx=8,int ny=8)
{DictionaryMesh d;d.Nx=nx;d.Ny=ny;d.maxCellNum=32768;d.maxInitialIterations=0;return d;}
Vector2D point(const QuadtreeMesh& mesh,int id){return mesh.verts(mesh.vertIndexToLabel[id]).coord();}
// Unpublished analytic first-jet samples; no geometric volume is asserted.
HermiteGalsPrimaryState primaryFor(const QuadtreeMesh& mesh,const Profile& f,bool radial=false)
{
    HermiteGalsPrimaryState p;p.measure=radial?HermiteGalsMeasure::AxisymmetricVolume:HermiteGalsMeasure::CartesianArea;
    for(int v=0;v<mesh.vertIndexNum();++v)if(mesh.verts(mesh.vertIndexToLabel[v]).type()==TVertType::non){
        const auto j=f(point(mesh,v));p.vertices.push_back({v,j.value,j.dx,j.dy});}
    std::set<Key> sides;
    for(int c=0;c<mesh.cellIndexNum();++c){const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);
        for(int k=0;k<4;++k){const auto& a=mesh.verts(cell.cellToVert[k]);const auto& b=mesh.verts(cell.cellToVert[(k+1)%4]);
            const bool horizontal=a.coord().y==b.coord().y;Key key{{a.index(),b.index()}};
            if(horizontal?a.coord().x>b.coord().x:a.coord().y>b.coord().y)std::swap(key[0],key[1]);
            if(sides.insert(key).second){const auto lo=point(mesh,key[0]),hi=point(mesh,key[1]);const auto j=f(lo+(hi-lo)*.5);
                p.edges.push_back({key,j.value,horizontal?j.dy:j.dx});}
        }
    }
    return p;
}
Jet monomial(Vector2D p,int a,int b)
{
    const auto power=[](double x,int n){return std::pow(x,n);};
    return {power(p.x,a)*power(p.y,b),a?a*power(p.x,a-1)*power(p.y,b):0,
        b?b*power(p.x,a)*power(p.y,b-1):0,a && b?a*b*power(p.x,a-1)*power(p.y,b-1):0};
}
Jet combination(Vector2D p)
{
    Jet result{};
    for(int a=0;a<=4;++a)for(int b=0;b<=4;++b){const auto j=monomial(p,a,b);
        const double w=((a+b)%2?-1:1)*(a+2*b+1)/32.0;
        result.value+=w*j.value;result.dx+=w*j.dx;result.dy+=w*j.dy;result.dxy+=w*j.dxy;}
    return result;
}
Jet positive(Vector2D p)
{
    const double f=1+.2*p.x+.1*std::pow(p.x,4),fx=.2+.4*std::pow(p.x,3);
    const double g=1+.3*p.y+.05*std::pow(p.y,4),gy=.3+.2*std::pow(p.y,3);
    return {f*g,fx*g,f*gy,fx*gy};
}
HermiteGalsCandidate boundSource(const QuadtreeMesh& mesh)
{
    auto primary=primaryFor(mesh,positive);Patch patch;HermiteGalsClosure closure;
    patch.activeCells.assign(mesh.cellIndexNum(),1);patch.nodes.resize(mesh.vertIndexNum());
    for(int v=0;v<mesh.vertIndexNum();++v){patch.nodes[v]=positive(point(mesh,v));closure.vertexMixedDerivatives.push_back(patch.nodes[v].dxy);}
    for(const auto& e:primary.edges)patch.edgeMidpoints.push_back({e.vertices,e.phi,e.dn});
    for(int c=0;c<mesh.cellIndexNum();++c){primary.geometryCells.push_back(c);primary.phaseVolumes.push_back({c,0});
        const auto j=positive(mesh.cells(mesh.cellIndexToLabel[c]).centerCoord());patch.cellCenters.push_back(j.value);
        closure.cellBubbleCoefficients.push_back(j.value);}
    HermiteGalsCandidate candidate(mesh,{},primary);
    candidate.bindGeometry(mesh,candidate.identity(),closure,LevelSetHermiteGeometry::Snapshot::fromNativeLeaves(mesh,patch,false));
    return candidate;
}
void unchanged(const HermiteGalsPrimaryState& a,const HermiteGalsPrimaryState& b)
{
    check(a.measure==b.measure && a.geometryCells==b.geometryCells && a.vertices.size()==b.vertices.size() &&
        a.edges.size()==b.edges.size() && a.phaseVolumes.size()==b.phaseVolumes.size(),"closure changed primary layout");
    for(std::size_t k=0;k<a.vertices.size();++k){const auto& x=a.vertices[k];const auto& y=b.vertices[k];
        check(x.vertex==y.vertex && x.phi==y.phi && x.dx==y.dx && x.dy==y.dy,"closure modified first jet");}
    for(std::size_t k=0;k<a.edges.size();++k){const auto& x=a.edges[k];const auto& y=b.edges[k];
        check(x.vertices==y.vertices && x.phi==y.phi && x.dn==y.dn,"closure modified edge data");}
    for(std::size_t k=0;k<a.phaseVolumes.size();++k)
        check(a.phaseVolumes[k].cell==b.phaseVolumes[k].cell && a.phaseVolumes[k].liquid==b.phaseVolumes[k].liquid,"closure modified inventory");
}

void scalarStencil()
{
    std::array<double,9> x{};for(int i=0;i<9;++i)x[i]=2.13+.17*i;
    double maximum=0;int samples=0;
    for(int degree=0;degree<=4;++degree)for(int target:{0,2,4,6,8}){
        const auto f=[&](int i){return std::pow(x[i]-2.13,degree);};
        const double got=Quantum::Geometry::compactQuarticDerivativeAtNode(target,9,[&](int i){return x[i];},f,check);
        const double exact=degree?degree*std::pow(x[target]-2.13,degree-1):0;
        close(got,exact,2e-12,"quartic one-dimensional derivative reproduction");
        if(!degree)check(got==0,"constant derivative must be exactly zero");
        maximum=std::max(maximum,std::abs(got-exact));++samples;
    }
    rejects([&]{Quantum::Geometry::compactQuarticDerivativeAtNode(0,4,[&](int i){return x[i];},[](int){return 1.;},check);},"insufficient",false);
    x[4]=x[3];rejects([&]{Quantum::Geometry::compactQuarticDerivativeAtNode(4,9,[&](int i){return x[i];},[](int){return 1.;},check);},"non-increasing",false);
    std::cout<<"MIXED_STENCIL,"<<samples<<','<<maximum<<'\n';
}

void reproduction(bool refined)
{
    auto d=input(refined?3:5,refined?2:3);d.minCoord={.13,-.31};d.maxCoord={1.27,.82};
    if(refined){d.minCoord={.125,-.25};d.maxCoord={1.125,.75};}
    QuadtreeMesh mesh(d);
    if(refined){DictionaryMesh::RefineBox box;box.min=d.minCoord;box.max=d.maxCoord;box.level=1;mesh.applyRefineBoxes({box});}
    const UniformLeafGrid grid(mesh,CartesianGridScope::UniformLeaves);
    double maximum=0;int fields=0;
    const auto test=[&](const Profile& f){
        const auto p=primaryFor(mesh,f);HermiteGalsCandidate candidate(mesh,{3,.25,.5,Position::Midpoint},p);
        const auto id=candidate.identity();const auto result=closeHermiteGalsMixedUniform(mesh,candidate);
        check(result.primaryIdentity()==id && result.timeLayer().time()==.375 && result.report().vertices==p.vertices.size(),
            "closure lost primary identity/time/layout");
        for(std::size_t k=0;k<p.vertices.size();++k){const double exact=f(point(mesh,p.vertices[k].vertex)).dxy;
            close(result.vertexMixedDerivatives()[k],exact,3e-11,"full Q4 mixed reproduction");
            maximum=std::max(maximum,std::abs(result.vertexMixedDerivatives()[k]-exact));}
        unchanged(p,candidate.primary());check(!candidate.hasGeometry() && candidate.identity()==id,"closure published a state");
        ++fields;
    };
    for(int a=0;a<=4;++a)for(int b=0;b<=4;++b)test([=](Vector2D p){return monomial(p,a,b);});
    test(combination);
    check(!refined || grid.level()==1,"refinement fixture did not create uniform actual leaves");
    std::cout<<"MIXED_Q4,"<<grid.nx()<<','<<grid.ny()<<','<<grid.level()<<','<<fields<<','<<maximum<<'\n';
}

void convergence(bool radial)
{
    double previous=0;
    for(int n:{8,16,32,64}){
        QuadtreeMesh mesh(input(n,n));
        const Profile f=[=](Vector2D p){const double a=.7,b=.4,e=std::exp(b*p.y);
            if(radial)return Jet{std::cos(a*p.x)*e,-a*std::sin(a*p.x)*e,b*std::cos(a*p.x)*e,-a*b*std::sin(a*p.x)*e};
            const double v=std::exp(a*p.x+b*p.y);return Jet{v,a*v,b*v,a*b*v};};
        const auto primary=primaryFor(mesh,f,radial);HermiteGalsCandidate candidate(mesh,{},primary);
        const auto result=closeHermiteGalsMixedUniform(mesh,candidate);
        double interior=0,boundary=0,axis=0;
        for(const auto& v:primary.vertices){const auto p=point(mesh,v.vertex);
            const double error=std::abs(result.vertexMixedDerivatives()[v.vertex]-f(p).dxy);
            if(p.x==0 || p.x==1 || p.y==0 || p.y==1)boundary=std::max(boundary,error);else interior=std::max(interior,error);
            if(radial && p.x==0)axis=std::max(axis,std::abs(result.vertexMixedDerivatives()[v.vertex]));}
        const double error=std::max(interior,boundary),order=previous?std::log2(previous/error):0;
        check(!radial || axis==0,"RZ mixed derivative axis parity must be exact");
        std::cout<<"MIXED_ORDER,"<<(radial?"RZ":"XY")<<','<<n<<','<<interior<<','<<boundary<<','<<error<<','<<order<<','<<axis<<'\n';
        if(n>=32)check(order>=3.7 && order<=4.3,"mixed closure fourth-order convergence");
        previous=error;
    }
}
void projectionAndFailures()
{
    QuadtreeMesh mesh(input(4,4));
    const Profile inconsistent=[](Vector2D p){return Jet{p.x*p.y,4*p.y,2*p.x,0};};
    auto primary=primaryFor(mesh,inconsistent);
    for(int c=0;c<mesh.cellIndexNum();++c){primary.geometryCells.push_back(c);primary.phaseVolumes.push_back({c,.01});}
    HermiteGalsCandidate candidate(mesh,{},primary);const auto originalId=candidate.identity();
    const auto result=closeHermiteGalsMixedUniform(mesh,candidate);
    for(double value:result.vertexMixedDerivatives())close(value,3,2e-12,"directional estimates must have fixed equal weights");
    close(result.report().maxRelativeDirectionalDisagreement,.5,2e-12,"directional disagreement diagnostic");
    auto changed=primary;for(auto& v:changed.vertices)v.phi=1000+v.vertex;for(auto& e:changed.edges)e.phi=-1000;
    candidate.replacePrimary(mesh,changed);const auto independent=closeHermiteGalsMixedUniform(mesh,candidate);
    check(independent.vertexMixedDerivatives()==result.vertexMixedDerivatives(),"scalar values entered mixed closure");
    check(result.primaryIdentity()==originalId && candidate.identity()!=originalId,"old result rebound to new primary revision");
    unchanged(changed,candidate.primary());check(!candidate.hasGeometry(),"mixed closure published geometry");
    const Profile affine=[](Vector2D p){return Jet{1+.7*p.x-.3*p.y,.7,-.3,0};};
    HermiteGalsCandidate constant(mesh,{},primaryFor(mesh,affine));
    const auto constantResult=closeHermiteGalsMixedUniform(mesh,constant);
    for(double value:constantResult.vertexMixedDerivatives())check(value==0,"constant gradients did not close to exact zero");
    auto wrongMesh=mesh.clone();wrongMesh->verts(wrongMesh->vertIndexToLabel[0]).setCoord({-.1,0});
    rejects([&]{closeHermiteGalsMixedUniform(*wrongMesh,candidate);},"mesh geometry/topology changed");
    const auto rejectPrimary=[&](const HermiteGalsPrimaryState& p,const char* reason){
        HermiteGalsCandidate bad(mesh,{},p);const auto id=bad.identity();
        rejects([&]{closeHermiteGalsMixedUniform(mesh,bad);},reason);
        unchanged(p,bad.primary());check(bad.identity()==id && !bad.hasGeometry(),"failure polluted primary revision");
    };
    auto bad=primary;bad.vertices.erase(bad.vertices.begin());rejectPrimary(bad,"complete uniform vertex skeleton");
    bad=primary;bad.edges.pop_back();rejectPrimary(bad,"complete uniform edge skeleton");
    const UniformLeafGrid grid(mesh,CartesianGridScope::UniformLeaves);
    const auto at=[&](int i,int j){return mesh.verts(grid.vertices()[j*(grid.nx()+1)+i]).index();};
    bad=primary;const Key first{{at(0,0),at(1,0)}};
    auto edge=std::find_if(bad.edges.begin(),bad.edges.end(),[&](const auto& e){return e.vertices==first;});
    edge->vertices={{at(0,0),at(2,0)}};rejectPrimary(bad,"missing native uniform master edge");
    bad=primary;bad.vertices[0].dy=std::numeric_limits<double>::max();
    for(auto& e:bad.edges)e.dn=-std::numeric_limits<double>::max();rejectPrimary(bad,"nonfinite compact mixed derivative");
    auto published=boundSource(mesh);const auto boundId=published.identity();
    rejects([&]{closeHermiteGalsMixedUniform(mesh,published);},"unpublished primary revision");
    check(published.hasGeometry() && published.identity()==boundId,"rejected published state was altered");
    QuadtreeMesh small(input(1,4));HermiteGalsCandidate tooSmall(small,{},primaryFor(small,affine));
    rejects([&]{closeHermiteGalsMixedUniform(small,tooSmall);},"at least two cells");
    QuadtreeMesh amr(input(4,4));DictionaryMesh::RefineBox box;box.min={.25,.25};box.max={.5,.5};box.level=1;
    amr.applyRefineBoxes({box});HermiteGalsCandidate adaptive(amr,{},primaryFor(amr,affine));
    rejects([&]{closeHermiteGalsMixedUniform(amr,adaptive);},"UniformLeafGrid");

    const Profile even=[](Vector2D p){const double f=p.x*p.x+.1*std::pow(p.x,4),fx=2*p.x+.4*std::pow(p.x,3);
        const double g=1+p.y+p.y*p.y+std::pow(p.y,4),gy=1+2*p.y+4*std::pow(p.y,3);
        return Jet{1+f*g,fx*g,f*gy,fx*gy};};
    const auto radial=primaryFor(mesh,even,true);HermiteGalsCandidate rz(mesh,{},radial);
    const auto rzResult=closeHermiteGalsMixedUniform(mesh,rz);double maxError=0;
    for(const auto& v:radial.vertices){const auto p=point(mesh,v.vertex);const double got=rzResult.vertexMixedDerivatives()[v.vertex];
        close(got,even(p).dxy,3e-11,"RZ even Q4 reproduction");maxError=std::max(maxError,std::abs(got-even(p).dxy));
        if(p.x==0)check(got==0,"RZ axis mixed derivative not exactly zero");}
    bad=radial;bad.vertices[at(0,1)].dx=1e-20;rejectPrimary(bad,"nonzero radial first jet on axis");
    bad=radial;const Key axial{{at(0,0),at(0,1)}};
    std::find_if(bad.edges.begin(),bad.edges.end(),[&](const auto& e){return e.vertices==axial;})->dn=1e-20;
    rejectPrimary(bad,"nonzero radial edge derivative on axis");
    auto negativeInput=input(4,4);negativeInput.minCoord={-.5,0};negativeInput.maxCoord={.5,1};
    QuadtreeMesh negative(negativeInput);
    // I01 rejects a negative-radius candidate before this closure can run.
    rejects([&]{HermiteGalsCandidate negativeR(negative,{},primaryFor(negative,even,true));},
        "RZ requires nonnegative radius",false);
    std::cout<<"MIXED_PROJECTION,3,"<<result.report().maxRelativeDirectionalDisagreement<<'\n';
    std::cout<<"MIXED_RZ_Q4,"<<radial.vertices.size()<<','<<maxError<<'\n';
    std::cout<<"mixed closure provenance, coverage, topology, overflow and failure isolation PASS\n";
}

// Independent closed-form uniform stencils, not the production Lagrange
// routine. h is the full leaf spacing; samples are half a leaf apart.
double referenceDerivative(const std::array<double,5>& f,int center,double h)
{
    if(center==0)return (-25*f[0]+48*f[1]-36*f[2]+16*f[3]-3*f[4])/(6*h);
    if(center==4)return (3*f[0]-16*f[1]+36*f[2]-48*f[3]+25*f[4])/(6*h);
    return (f[0]-8*f[1]+8*f[3]-f[4])/(6*h);
}
void transportedSkeleton()
{
    QuadtreeMesh mesh(input(8,8));auto initial=boundSource(mesh);const auto source=std::move(initial).finish();
    const auto primaryBefore=source.primary();const auto sourceView=source.geometry();
    const UniformLeafGrid grid(mesh,CartesianGridScope::UniformLeaves);
    const LevelSetBoundaryTrace::Domain domain{mesh.minCoord,mesh.maxCoord,grid.spacing(),false};
    const HermiteGalsVelocityView::Query velocity=[](const Vector2D& p){return LevelSetBoundaryTrace::VelocityJet{
        {.12*p.x*(1-p.x),-.09*p.y*(1-p.y)},{.12*(1-2*p.x),0,0,-.09*(1-2*p.y)}};};
    const HermiteGalsVelocityView old(source.identity(),source.timeLayer(),velocity);
    const HermiteGalsVelocityView now(source.identity(),{1,0,.1,Position::End},velocity);
    std::vector<int> vertexIds;for(const auto& v:source.primary().vertices)vertexIds.push_back(v.vertex);
    std::vector<Key> edgeIds;for(const auto& e:source.primary().edges)edgeIds.push_back(e.vertices);std::sort(edgeIds.begin(),edgeIds.end());
    const auto vertices=transportHermiteGalsVertices(mesh,source,vertexIds,domain,old,now);
    const auto edges=transportHermiteGalsEdges(mesh,source,edgeIds,domain,old,now);
    auto advanced=source.primary();advanced.vertices=vertices.vertices();advanced.edges=edges.edges();
    HermiteGalsCandidate candidate(mesh,now.timeLayer(),advanced);const auto closed=closeHermiteGalsMixedUniform(mesh,candidate);
    check(closed.primaryIdentity()==candidate.identity() && closed.primaryIdentity()!=source.identity() &&
        closed.timeLayer().position==Position::End && closed.timeLayer().time()==.1,"transport closure revision/time is ambiguous");
    std::map<Key,double> normals;for(const auto& e:advanced.edges)normals[e.vertices]=e.dn;
    const int nx=grid.nx(),ny=grid.ny();
    const auto id=[&](int i,int j){return mesh.verts(grid.vertices()[j*(nx+1)+i]).index();};
    double maximum=0;
    for(int j=0;j<=ny;++j)for(int i=0;i<=nx;++i){
        std::array<double,5> row{},column{};const int startX=std::clamp(2*i-2,0,2*nx-4),startY=std::clamp(2*j-2,0,2*ny-4);
        for(int k=0;k<5;++k){const int x=startX+k,y=startY+k;
            row[k]=x%2?normals.at({{id(x/2,j),id(x/2+1,j)}}):advanced.vertices[id(x/2,j)].dy;
            column[k]=y%2?normals.at({{id(i,y/2),id(i,y/2+1)}}):advanced.vertices[id(i,y/2)].dx;}
        const double exact=.5*referenceDerivative(row,2*i-startX,grid.spacing().x)+
            .5*referenceDerivative(column,2*j-startY,grid.spacing().y);
        const double got=closed.vertexMixedDerivatives()[id(i,j)];close(got,exact,3e-11,"transported skeleton closure/formula wiring");
        maximum=std::max(maximum,std::abs(got-exact));
    }
    unchanged(advanced,candidate.primary());unchanged(primaryBefore,source.primary());
    check(!candidate.hasGeometry() && source.identity()==sourceView.identity() && &source.geometry().snapshot()==&sourceView.snapshot(),
        "closure changed old geometry or prematurely published candidate");
    std::cout<<"MIXED_TRANSPORT,"<<advanced.vertices.size()<<','<<advanced.edges.size()<<','<<maximum<<'\n';
}
} // namespace

int main(int argc,char** argv)
{
    try{
        std::cout<<std::setprecision(17);
        if(argc==2 && std::string(argv[1])=="--legacy-extraction-control"){interfaceGalsContract();return 0;}
        scalarStencil();reproduction(false);reproduction(true);projectionAndFailures();transportedSkeleton();
        convergence(false);convergence(true);interfaceGalsContract();
    }catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
    std::cout<<"I05 uniform mixed closure contract PASS (no AMR closure or geometry publication)\n";
    return 0;
}
