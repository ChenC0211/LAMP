#include "Modules/Interface/HermiteGalsAmrClosure.h"

#include "Modules/Interface/HermiteGalsMixedClosure.h"
#include "Modules/Interface/HermiteGalsVertexTransport.h"
#include "Modules/Interface/HermiteGalsEdgeTransport.h"
#include "Core/Geometry/ActiveQuarticHermitePatch.h"
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <cstdint>
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
#include <tuple>

namespace
{
using namespace Quantum::Interface;
using Skeleton=Quantum::Geometry::ActiveQuarticHermiteSkeleton;
using Input=Quantum::Geometry::ActiveQuarticHermiteSkeletonInput;
using Sensitivity=Quantum::Geometry::ActiveQuarticHermiteSensitivity;
using Patch=Quantum::Geometry::ActiveQuarticHermitePatch;
using Jet=Patch::NodalDof;
using Key=std::array<int,2>;
using Profile=std::function<Jet(Vector2D)>;
constexpr Key none{{-1,-1}};
static_assert(!std::is_default_constructible_v<HermiteGalsAmrClosureResult>);
static_assert(std::is_same_v<decltype(std::declval<HermiteGalsAmrClosureResult&>().vertexMixedDerivatives()),const std::vector<double>&>);
static_assert(std::is_same_v<decltype(std::declval<HermiteGalsAmrClosureResult&>().skeleton()),const Skeleton&>);
void check(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
double close(double a,double b,const char* why,double tolerance=3e-11)
{if(!std::isfinite(a) || std::abs(a-b)>tolerance*(1+std::abs(b)))throw std::runtime_error(std::string(why)+" got="+std::to_string(a)+" reference="+std::to_string(b));return std::abs(a-b);}
template<class F> void rejects(F action,const char* why,const char* owner="HermiteGalsAmrClosure primary=")
{
    try{action();}catch(const std::exception& e){const std::string s=e.what();
        if(s.find(why)==std::string::npos || s.find(owner)==std::string::npos)throw std::runtime_error("unexpected rejection: "+s);return;}
    throw std::runtime_error(std::string("missing rejection: ")+why);
}
Vector2D point(const QuadtreeMesh& m,int v){return m.verts(m.vertIndexToLabel[v]).coord();}
bool horizontal(const QuadtreeMesh& m,Key k){return point(m,k[0]).y==point(m,k[1]).y;}
std::unique_ptr<QuadtreeMesh> makeMesh(int n,int pattern,bool radial=false,bool shifted=false)
{
    DictionaryMesh d;d.Nx=d.Ny=n;d.maxCellNum=32768;d.maxInitialIterations=0;
    if(shifted){d.minCoord={radial?0:.13,-.27};d.maxCoord={1.87,1.36};}
    auto m=std::make_unique<QuadtreeMesh>(d);
    if(pattern){const auto lo=m->minCoord,span=m->maxCoord-lo;std::vector<DictionaryMesh::RefineBox> boxes;
        const int levels=pattern==3?3:1;
        for(int level=1;level<=levels;++level){DictionaryMesh::RefineBox b;b.level=level;
            const double a=level==1?.26:(level==2?.38:.45),z=level==1?(pattern==1?.49:.74):(level==2?.62:.55);
            b.min={radial?lo.x:lo.x+a*span.x,lo.y+a*span.y};
            b.max={lo.x+(radial?(level==1?(pattern==3?.49:.24):(level==2?.24:.12)):z)*span.x,lo.y+z*span.y};
            if(pattern==4){b.min={lo.x+.25*span.x,lo.y+.25*span.y};b.max={lo.x+.75*span.x,lo.y+.75*span.y};}
            boxes.push_back(b);
        }
        m->applyRefineBoxes(boxes);int maxLevel=0;
        for(int c=0;c<m->cellIndexNum();++c)maxLevel=std::max(maxLevel,m->cells(m->cellIndexToLabel[c]).level());
        check(maxLevel==levels,"AMR test refinement missed cells");
    }
    return m;
}
std::set<Key> actualSides(const QuadtreeMesh& m)
{
    std::set<Key> result;
    for(int c=0;c<m.cellIndexNum();++c){const auto& cell=m.cells(m.cellIndexToLabel[c]);
        for(int k=0;k<4;++k){Key key{{m.verts(cell.cellToVert[k]).index(),m.verts(cell.cellToVert[(k+1)%4]).index()}};
            const auto a=point(m,key[0]),b=point(m,key[1]);if(a.y==b.y?a.x>b.x:a.y>b.y)std::swap(key[0],key[1]);result.insert(key);}}
    return result;
}
bool constrained(const QuadtreeMesh& m,Key key)
{
    const auto a=point(m,key[0]),b=point(m,key[1]);const bool h=horizontal(m,key);
    for(int v:key){const int owner=m.verts(m.vertIndexToLabel[v]).supCell;if(owner<0)continue;
        const auto& cell=m.cells(owner);const auto lo=m.verts(cell.cellToVert[0]).coord(),hi=m.verts(cell.cellToVert[2]).coord();
        if(h && (a.y==lo.y || a.y==hi.y) && a.x>=lo.x && b.x<=hi.x && b.x-a.x<hi.x-lo.x)return true;
        if(!h && (a.x==lo.x || a.x==hi.x) && a.y>=lo.y && b.y<=hi.y && b.y-a.y<hi.y-lo.y)return true;}
    return false;
}
HermiteGalsPrimaryState primaryFor(const QuadtreeMesh& m,const Profile& f,bool radial=false)
{
    HermiteGalsPrimaryState p;p.measure=radial?HermiteGalsMeasure::AxisymmetricVolume:HermiteGalsMeasure::CartesianArea;
    for(int v=0;v<m.vertIndexNum();++v)if(m.verts(m.vertIndexToLabel[v]).type()==TVertType::non){const auto j=f(point(m,v));p.vertices.push_back({v,j.value,j.dx,j.dy});}
    for(const auto& key:actualSides(m))if(!constrained(m,key)){const auto a=point(m,key[0]),b=point(m,key[1]);const auto j=f(a+(b-a)*.5);
        p.edges.push_back({key,j.value,horizontal(m,key)?j.dy:j.dx});}
    for(int c=0;c<m.cellIndexNum();++c){p.geometryCells.push_back(c);p.phaseVolumes.push_back({c,0});}
    return p; // unpublished analytic samples; no volume initialization claim
}
Jet q4(Vector2D p,int term=-1,bool radial=false)
{
    Jet r{};
    for(int a=0;a<5;++a)for(int b=0;b<5;++b){if(radial && a%2)continue;if(term>=0 && term!=5*a+b)continue;
        const double w=term>=0?1:((a+b)%2?-1:1)*(a+2*b+1)/50.;
        const double x=std::pow(p.x,a),y=std::pow(p.y,b),dx=a?a*std::pow(p.x,a-1):0,dy=b?b*std::pow(p.y,b-1):0;
        r.value+=w*x*y;r.dx+=w*dx*y;r.dy+=w*x*dy;r.dxy+=w*dx*dy;}
    return r;
}
void unchanged(const HermiteGalsPrimaryState& a,const HermiteGalsPrimaryState& b)
{
    check(a.vertices.size()==b.vertices.size() && a.edges.size()==b.edges.size() && a.geometryCells==b.geometryCells && a.phaseVolumes.size()==b.phaseVolumes.size(),"primary topology changed");
    for(std::size_t i=0;i<a.vertices.size();++i){const auto& x=a.vertices[i];const auto& y=b.vertices[i];check(x.vertex==y.vertex && x.phi==y.phi && x.dx==y.dx && x.dy==y.dy,"master jet changed");}
    for(std::size_t i=0;i<a.edges.size();++i){const auto& x=a.edges[i];const auto& y=b.edges[i];check(x.vertices==y.vertices && x.phi==y.phi && x.dn==y.dn,"master edge changed");}
    for(std::size_t i=0;i<a.phaseVolumes.size();++i)check(a.phaseVolumes[i].cell==b.phaseVolumes[i].cell && a.phaseVolumes[i].liquid==b.phaseVolumes[i].liquid,"conserved inventory changed");
}
double c1(const QuadtreeMesh& m,const Skeleton& skeleton)
{
    double error=0;std::size_t fine=0;
    for(const auto& e:skeleton.edges())if(e.parentSide!=none){++fine;const auto a=point(m,e.vertices[0]),b=point(m,e.vertices[1]);const bool h=horizontal(m,e.vertices);
        for(int i=0;i<=8;++i){const auto p=i==0?a:(i==8?b:a+(b-a)*(i/8.));const double s=h?p.x:p.y;
            const auto x=skeleton.edgeJet(e.vertices,s),y=skeleton.edgeJet(e.parentSide,s);
            for(int k=0;k<4;++k)error=std::max(error,close(x[k],y[k],"AMR coarse/fine C1"));}}
    check(fine>0,"C1 fixture has no fine sides");return error;
}
// Independent Vandermonde solve, using coordinates shifted to target/scaled by
// support span; coefficient of t is the derivative. No Core stencil call.
double independentDerivative(const std::vector<std::pair<double,double>>& samples,int target)
{
    const int begin=std::clamp(target-2,0,static_cast<int>(samples.size())-5);long double matrix[5][6]{};
    const long double x0=samples[target].first,scale=static_cast<long double>(samples[begin+4].first)-samples[begin].first;
    for(int i=0;i<5;++i){const long double t=(samples[begin+i].first-x0)/scale;matrix[i][0]=1;
        for(int j=1;j<5;++j)matrix[i][j]=matrix[i][j-1]*t;matrix[i][5]=samples[begin+i].second;}
    for(int k=0;k<5;++k){int pivot=k;for(int i=k+1;i<5;++i)if(std::abs(matrix[i][k])>std::abs(matrix[pivot][k]))pivot=i;
        check(matrix[pivot][k]!=0,"singular reference support");for(int j=k;j<=5;++j)std::swap(matrix[pivot][j],matrix[k][j]);
        const auto d=matrix[k][k];for(int j=k;j<=5;++j)matrix[k][j]/=d;
        for(int i=0;i<5;++i)if(i!=k){const auto f=matrix[i][k];for(int j=k;j<=5;++j)matrix[i][j]-=f*matrix[k][j];}}
    return static_cast<double>(matrix[1][5]/scale);
}
double directCompatibility(const QuadtreeMesh& m,const HermiteGalsPrimaryState& p,const HermiteGalsAmrClosureResult& result)
{
    // Build connected components independently as point adjacency on each
    // primary side. A component may terminate at a hanging perpendicular trace.
    std::map<std::pair<int,int>,double> estimates;
    for(int direction=0;direction<2;++direction){const bool h=direction==0;
        std::map<double,std::map<double,Key>> rows;
        for(const auto& e:p.edges)if(horizontal(m,e.vertices)==h){const auto a=point(m,e.vertices[0]);rows[h?a.y:a.x][h?a.x:a.y]=e.vertices;}
        for(const auto& row:rows){std::vector<std::pair<double,double>> values;std::vector<std::pair<int,int>> vertices;
            double end=0;bool started=false;
            const auto flush=[&]{if(values.empty())return;check(values.size()>=5,"reference component too short");
                for(const auto& v:vertices)if(m.verts(m.vertIndexToLabel[v.first]).type()==TVertType::non)
                    estimates[{v.first,direction}]=(p.measure==HermiteGalsMeasure::AxisymmetricVolume && h && point(m,v.first).x==0)?0:independentDerivative(values,v.second);};
            for(const auto& side:row.second){const auto key=side.second;const auto a=point(m,key[0]),b=point(m,key[1]);const double low=h?a.x:a.y,high=h?b.x:b.y;
                if(!started || low!=end){flush();values.clear();vertices.clear();values.push_back({low,result.skeleton().nodeJet(key[0])[h?2:1]});vertices.push_back({key[0],0});started=true;}
                values.push_back({low+(high-low)*.5,result.skeleton().edgeMidpoint(key)[1]});
                vertices.push_back({key[1],static_cast<int>(values.size())});values.push_back({high,result.skeleton().nodeJet(key[1])[h?2:1]});end=high;
            }flush();
        }
    }
    double error=0;for(std::size_t i=0;i<p.vertices.size();++i){const int v=p.vertices[i].vertex;
        error=std::max(error,close(result.vertexMixedDerivatives()[i],.5*estimates.at({v,0})+.5*estimates.at({v,1}),"independent compact compatibility"));}
    return error;
}
void reproduction()
{
    for(bool radial:{false,true})for(int pattern=0;pattern<4;++pattern){auto m=makeMesh(4,pattern,radial,true);double master=0,hanging=0,trace=0,uniform=0;std::size_t nonzeros=0;int fields=0;
        for(int term=-1;term<(pattern==1?25:0);++term){if(radial && term>=0 && (term/5)%2)continue;++fields;
            const Profile f=[=](Vector2D p){return q4(p,term,radial);};const auto p=primaryFor(*m,f,radial);HermiteGalsCandidate candidate(*m,{},p);
            const auto result=closeHermiteGalsMixedAmr(*m,candidate);nonzeros=std::max(nonzeros,result.report().compatibilityNonzeros);
            for(std::size_t i=0;i<p.vertices.size();++i)master=std::max(master,close(result.vertexMixedDerivatives()[i],f(point(*m,p.vertices[i].vertex)).dxy,"AMR master Q4"));
            for(const auto& v:result.skeleton().vertices()){const auto j=result.skeleton().nodeJet(v.index);const auto exact=f(point(*m,v.index));
                hanging=std::max(hanging,close(j[3],exact.dxy,"AMR constrained Q4"));
                if(radial && point(*m,v.index).x==0)check(j[1]==0 && j[3]==0,"AMR axis parity");}
            if(pattern)trace=std::max(trace,c1(*m,result.skeleton()));
            else {const auto old=closeHermiteGalsMixedUniform(*m,candidate);for(std::size_t i=0;i<p.vertices.size();++i)
                uniform=std::max(uniform,close(result.vertexMixedDerivatives()[i],old.vertexMixedDerivatives()[i],"uniform I05 equivalence"));}
            unchanged(p,candidate.primary());check(!candidate.hasGeometry() && result.primaryIdentity()==candidate.identity(),"closure publication/identity");
        }
        std::cout<<"AMR_Q4,"<<(radial?"RZ":"XY")<<','<<pattern<<','<<m->cellIndexNum()<<','<<fields<<','<<master<<','<<hanging<<','<<trace<<','<<uniform<<','<<nonzeros<<'\n';
    }
    auto m=makeMesh(4,3);const Profile affine=[](Vector2D p){return Jet{2+p.x+2*p.y,1,2,0};};HermiteGalsCandidate candidate(*m,{},primaryFor(*m,affine));
    const auto closed=closeHermiteGalsMixedAmr(*m,candidate);for(double q:closed.vertexMixedDerivatives())check(q==0,"constant gradients must close to exactly zero");
    std::cout<<"AMR_CONSTANT,"<<closed.vertexMixedDerivatives().size()<<",0\n";
}
void arbitraryAndSensitivity()
{
    auto m=makeMesh(4,3,false,true);auto p=primaryFor(*m,[](Vector2D q){return q4(q);});
    std::uint32_t seed=89123;const auto random=[&]{seed=1664525u*seed+1013904223u;return static_cast<double>(seed)/4294967296.0-.5;};
    for(auto& v:p.vertices){v.phi=random();v.dx=random();v.dy=random();}for(auto& e:p.edges){e.phi=random();e.dn=random();}
    HermiteGalsCandidate candidate(*m,{},p);const auto result=closeHermiteGalsMixedAmr(*m,candidate);
    const double defect=c1(*m,result.skeleton()),reference=directCompatibility(*m,p,result);
    Input input;input.activeCells.assign(m->cellIndexNum(),1);input.nodes.resize(m->vertIndexNum());
    for(std::size_t i=0;i<p.vertices.size();++i){const auto& v=p.vertices[i];input.nodes[v.vertex]={v.phi,v.dx,v.dy,result.vertexMixedDerivatives()[i]};}
    for(const auto& e:p.edges)input.edgeMidpoints.push_back({e.vertices,e.phi,e.dn});
    auto perturbed=input;std::map<std::tuple<int,int,int>,double> delta;
    for(const auto& v:p.vertices){const std::array<double,4> change{{random(),random(),random(),1e7*random()}};
        auto& j=perturbed.nodes[v.vertex];j.value+=change[0];j.dx+=change[1];j.dy+=change[2];j.dxy+=change[3];
        const auto& old=input.nodes[v.vertex];const std::array<double,4> actual{{j.value-old.value,j.dx-old.dx,j.dy-old.dy,j.dxy-old.dxy}};
        for(int k=0;k<4;++k)delta[{0,v.vertex,k}]=actual[k];}
    for(std::size_t i=0;i<perturbed.edgeMidpoints.size();++i){auto& e=perturbed.edgeMidpoints[i];e.value+=random();e.normalDerivative+=random();
        delta[{1,static_cast<int>(i),0}]=e.value-input.edgeMidpoints[i].value;delta[{1,static_cast<int>(i),1}]=e.normalDerivative-input.edgeMidpoints[i].normalDerivative;}
    const Skeleton changed(*m,perturbed);double sensitivityError=0,scaledSensitivity=0;std::size_t queries=0,mixedDependencies=0;
    for(const auto& v:result.skeleton().vertices())for(int k=0;k<4;++k){long double prediction=0;
        for(const auto& s:result.skeleton().nodeSensitivity(v.index,k)){
            check(s.kind!=Sensitivity::Kind::CellCenter,"node sensitivity depends on interior center");
            prediction+=s.coefficient*delta.at({static_cast<int>(s.kind),s.index,s.component});
            if(v.masterSide!=none && k<3 && s.kind==Sensitivity::Kind::Vertex && s.component==3)++mixedDependencies;
        }
        const double expected=result.skeleton().nodeJet(v.index)[k]+static_cast<double>(prediction);
        const double error=close(changed.nodeJet(v.index)[k],expected,"constraint sensitivity superposition");
        sensitivityError=std::max(sensitivityError,error);scaledSensitivity=std::max(scaledSensitivity,error/(1+std::abs(expected)));++queries;
    }
    std::cout<<"AMR_DEPENDENCIES,"<<mixedDependencies<<','<<result.report().compatibilityNonzeros<<','<<p.vertices.size()<<'\n';
    // For these balanced traces, tangential hanging samples can be independent
    // of mixed inputs even when the full node jet is not. Record actual support;
    // do not claim an off-diagonal compatibility test when it reduces to I.
    rejects([&]{result.skeleton().nodeSensitivity(-1,0);},"component 0..3","ActiveQuarticHermiteSkeleton");
    rejects([&]{result.skeleton().nodeSensitivity(p.vertices[0].vertex,4);},"component 0..3","ActiveQuarticHermiteSkeleton");
    unchanged(p,candidate.primary());
    std::cout<<"AMR_RANDOM,"<<p.vertices.size()<<','<<result.report().hangingVertices<<','<<result.report().constrainedEdges<<','<<result.report().compatibilityNonzeros<<','<<defect<<','<<reference<<','<<result.report().algebraicResidual<<','<<result.report().compatibilityResidual<<'\n';
    std::cout<<"AMR_SENSITIVITY,"<<queries<<','<<mixedDependencies<<','<<sensitivityError<<','<<scaledSensitivity<<'\n';
}
HermiteGalsPrimaryState subset(const QuadtreeMesh& m,HermiteGalsPrimaryState p,const std::function<bool(Vector2D)>& active)
{
    p.geometryCells.clear();std::set<int> vertices;std::set<Key> sides;
    for(int c=0;c<m.cellIndexNum();++c){const auto& cell=m.cells(m.cellIndexToLabel[c]);if(!active(cell.centerCoord()))continue;
        p.geometryCells.push_back(c);
        for(int k=0;k<4;++k){Key key{{m.verts(cell.cellToVert[k]).index(),m.verts(cell.cellToVert[(k+1)%4]).index()}};
            vertices.insert(key[0]);const auto a=point(m,key[0]),b=point(m,key[1]);if(a.y==b.y?a.x>b.x:a.y>b.y)std::swap(key[0],key[1]);sides.insert(key);}}
    p.vertices.erase(std::remove_if(p.vertices.begin(),p.vertices.end(),[&](const auto& v){return !vertices.count(v.vertex);}),p.vertices.end());
    p.edges.erase(std::remove_if(p.edges.begin(),p.edges.end(),[&](const auto& e){return !sides.count(e.vertices);}),p.edges.end());return p;
}
void activeSupport()
{
    const Profile field=[](Vector2D p){return q4(p);};auto m=makeMesh(8,1);const auto primary=primaryFor(*m,field);
    const auto p=subset(*m,primary,[](Vector2D q){return q.x<.75 && q.y<.75;});HermiteGalsCandidate candidate(*m,{},p);
    const auto result=closeHermiteGalsMixedAmr(*m,candidate);double error=0;
    for(std::size_t i=0;i<p.vertices.size();++i)error=std::max(error,close(result.vertexMixedDerivatives()[i],field(point(*m,p.vertices[i].vertex)).dxy,"supported AMR subset Q4"));
    const double trace=c1(*m,result.skeleton());
    auto uniform=makeMesh(8,0);const auto full=primaryFor(*uniform,field);
    const auto separate=subset(*uniform,full,[](Vector2D q){return q.x<.25 || q.x>.75;});HermiteGalsCandidate disconnected(*uniform,{},separate);
    const auto split=closeHermiteGalsMixedAmr(*uniform,disconnected);const double reference=directCompatibility(*uniform,separate,split);
    const auto shortLine=subset(*uniform,full,[](Vector2D q){return q.x<.125 || q.x>.75;});HermiteGalsCandidate broken(*uniform,{},shortLine);
    rejects([&]{closeHermiteGalsMixedAmr(*uniform,broken);},"insufficient five-point support");
    unchanged(shortLine,broken.primary());
    std::cout<<"AMR_ACTIVE_SUPPORT,"<<p.geometryCells.size()<<','<<p.vertices.size()<<','<<error<<','<<trace<<','<<split.report().lineComponents<<','<<reference<<",gap_rejected\n";
}
void convergence(bool radial)
{
    const Profile f=[=](Vector2D p){if(radial){const double x=std::cos(.7*p.x),dx=-.7*std::sin(.7*p.x),y=std::exp(.4*p.y);
        return Jet{x*y,dx*y,.4*x*y,.4*dx*y};}const double v=std::exp(.7*p.x+.4*p.y);return Jet{v,.7*v,.4*v,.28*v};};
    double previous=0;
    for(int n:{8,16,32,64}){auto m=makeMesh(n,4,radial);const auto p=primaryFor(*m,f,radial);HermiteGalsCandidate candidate(*m,{},p);const auto result=closeHermiteGalsMixedAmr(*m,candidate);
        double master=0,hanging=0;std::size_t axis=0;
        for(std::size_t i=0;i<p.vertices.size();++i)master=std::max(master,std::abs(result.vertexMixedDerivatives()[i]-f(point(*m,p.vertices[i].vertex)).dxy));
        for(const auto& v:result.skeleton().vertices()){const auto q=point(*m,v.index);const double dxy=result.skeleton().nodeJet(v.index)[3];
            if(v.masterSide!=none)hanging=std::max(hanging,std::abs(dxy-f(q).dxy));if(radial && q.x==0){check(dxy==0,"MMS axis mixed not zero");++axis;}}
        const double order=previous?std::log2(previous/master):0;
        std::cout<<"AMR_ORDER,"<<(radial?"RZ":"XY")<<','<<n<<','<<p.vertices.size()<<','<<result.report().hangingVertices<<','<<master<<','<<hanging<<','<<order<<','<<axis<<','<<result.report().compatibilityResidual<<'\n';
        if(n>=32)check(order>=2.7 && order<=4.3,"AMR master mixed convergence order");previous=master;
    }
}
Jet positive(Vector2D p,bool radial=false)
{
    const double x=radial?1+.2*p.x*p.x+.1*std::pow(p.x,4):1+.2*p.x+.1*std::pow(p.x,4);
    const double dx=radial?.4*p.x+.4*std::pow(p.x,3):.2+.4*std::pow(p.x,3);
    const double y=1+.3*p.y+.05*std::pow(p.y,4),dy=.3+.2*std::pow(p.y,3);return {x*y,dx*y,x*dy,dx*dy};
}
HermiteGalsCandidate bound(const QuadtreeMesh& m,bool radial=false)
{
    const auto f=[=](Vector2D q){return positive(q,radial);};const auto p=primaryFor(m,f,radial);HermiteGalsClosure closure;Patch patch;
    patch.activeCells.assign(m.cellIndexNum(),1);patch.nodes.resize(m.vertIndexNum());patch.cellCenters.resize(m.cellIndexNum());
    for(int v=0;v<m.vertIndexNum();++v)patch.nodes[v]=f(point(m,v));for(const auto& v:p.vertices)closure.vertexMixedDerivatives.push_back(f(point(m,v.vertex)).dxy);
    for(const auto& e:p.edges)patch.edgeMidpoints.push_back({e.vertices,e.phi,e.dn});
    for(int c=0;c<m.cellIndexNum();++c){const double value=f(m.cells(m.cellIndexToLabel[c]).centerCoord()).value;patch.cellCenters[c]=value;closure.cellBubbleCoefficients.push_back(value);}
    HermiteGalsCandidate candidate(m,{},p);auto snapshot=LevelSetHermiteGeometry::Snapshot::fromNativeLeaves(m,patch,radial);
    candidate.bindGeometry(m,candidate.identity(),closure,std::move(snapshot));return candidate;
}
void transport(bool radial)
{
    auto m=makeMesh(4,3,radial);auto initial=bound(*m,radial);const auto source=std::move(initial).finish();const auto before=source.primary();const auto view=source.geometry();
    Vector2D size=m->maxCoord-m->minCoord;for(int c=0;c<m->cellIndexNum();++c){const auto h=m->cells(m->cellIndexToLabel[c]).size();size.x=std::min(size.x,h.x);size.y=std::min(size.y,h.y);}
    LevelSetBoundaryTrace::Domain domain{m->minCoord,m->maxCoord,size,radial};
    const HermiteGalsVelocityView::Query velocity=[=](const Vector2D& p){return LevelSetBoundaryTrace::VelocityJet{
        {radial?.12*p.x*(1-p.x*p.x):.12*p.x*(1-p.x),-.09*p.y*(1-p.y)},
        {radial?.12*(1-3*p.x*p.x):.12*(1-2*p.x),0,0,-.09*(1-2*p.y)}};};
    const HermiteGalsVelocityView old(source.identity(),source.timeLayer(),velocity),now(source.identity(),{1,0,.05,HermiteGalsTimeLayer::Position::End},velocity);
    std::vector<int> vertices;for(const auto& v:before.vertices)vertices.push_back(v.vertex);std::vector<Key> edges;for(const auto& e:before.edges)edges.push_back(e.vertices);
    const auto movedVertices=transportHermiteGalsVertices(*m,source,vertices,domain,old,now);const auto movedEdges=transportHermiteGalsEdges(*m,source,edges,domain,old,now);
    auto p=before;p.vertices=movedVertices.vertices();p.edges=movedEdges.edges();HermiteGalsCandidate candidate(*m,now.timeLayer(),p);const auto result=closeHermiteGalsMixedAmr(*m,candidate);
    const double reference=directCompatibility(*m,p,result),defect=c1(*m,result.skeleton());
    check(result.primaryIdentity()==candidate.identity() && result.primaryIdentity()!=source.identity() && result.timeLayer().time()==.05,"AMR result provenance");
    check(vertices.size()+result.report().hangingVertices==static_cast<std::size_t>(m->vertIndexNum()),"hanging target was independently transported");
    unchanged(before,source.primary());unchanged(p,candidate.primary());check(!candidate.hasGeometry() && &source.geometry().snapshot()==&view.snapshot(),"AMR closure altered published geometry");
    const auto id=result.primaryIdentity();candidate.replacePrimary(*m,p);check(result.primaryIdentity()==id && candidate.identity()!=id,"result silently rebound after replacement");
    std::cout<<"AMR_TRANSPORT,"<<(radial?"RZ":"XY")<<','<<vertices.size()<<','<<edges.size()<<','<<result.report().hangingVertices<<','<<result.report().constrainedEdges<<','<<reference<<','<<defect<<'\n';
}
void failures()
{
    auto m=makeMesh(4,1);const auto p=primaryFor(*m,[](Vector2D q){return q4(q);});
    const auto fail=[&](const HermiteGalsPrimaryState& bad,const char* message){HermiteGalsCandidate candidate(*m,{},bad);const auto id=candidate.identity();rejects([&]{closeHermiteGalsMixedAmr(*m,candidate);},message);unchanged(bad,candidate.primary());check(id==candidate.identity() && !candidate.hasGeometry(),"failed closure polluted candidate");};
    auto bad=p;bad.vertices.pop_back();fail(bad,"non-finite");bad=p;bad.edges.pop_back();fail(bad,"missing independent");
    for(const auto& key:actualSides(*m))if(constrained(*m,key)){bad=p;bad.edges.push_back({key,0,0});fail(bad,"constrained fine edge");break;}
    bad=p;bad.geometryCells.clear();fail(bad,"active cells");
    for(int v=0;v<m->vertIndexNum();++v){const int coarse=m->verts(m->vertIndexToLabel[v]).supCell;if(coarse>=0){bad=p;const int c=m->cells(coarse).index();bad.geometryCells.erase(std::find(bad.geometryCells.begin(),bad.geometryCells.end(),c));fail(bad,"inactive leaf");break;}}
    auto small=makeMesh(1,0);HermiteGalsCandidate unsupported(*small,{},primaryFor(*small,[](Vector2D q){return q4(q);}));rejects([&]{closeHermiteGalsMixedAmr(*small,unsupported);},"five-point support");
    HermiteGalsCandidate stale(*m,{},p);auto changed=m->clone();changed->verts(changed->vertIndexToLabel[0]).setCoord({-.1,0});rejects([&]{closeHermiteGalsMixedAmr(*changed,stale);},"mesh");
    auto published=bound(*m);rejects([&]{closeHermiteGalsMixedAmr(*m,published);},"unpublished");
    auto rz=makeMesh(4,1,true);auto radial=primaryFor(*rz,[](Vector2D q){return positive(q,true);},true);radial.vertices[0].dx=1e-20;
    HermiteGalsCandidate brokenAxis(*rz,{},radial);rejects([&]{closeHermiteGalsMixedAmr(*rz,brokenAxis);},"axis dx");
    bad=p;bad.vertices[0].dy=std::numeric_limits<double>::max();for(auto& e:bad.edges)e.dn=-std::numeric_limits<double>::max();fail(bad,"overflow");
    std::cout<<"AMR_REJECTIONS,coverage_support_identity_axis_overflow_isolation,PASS\n";
}
void declaredSymmetry()
{
    using Plane=LevelSetInterfacePanels::SymmetryPlane;
    int cases=0,guards=0;double maxMixed=0,maxTrace=0;
    const auto verify=[&](const QuadtreeMesh& m,const Profile& field,bool radial,const std::vector<Plane>& planes){
        const auto p=primaryFor(m,field,radial);HermiteGalsCandidate candidate(m,{},p);const auto id=candidate.identity();
        const auto result=closeHermiteGalsMixedAmr(m,candidate,planes);check(result.report().symmetryConstrainedVertices>0,"missing declared symmetry constraints");
        Patch patch;patch.activeCells.assign(m.cellIndexNum(),1);patch.nodes.resize(m.vertIndexNum());patch.cellCenters.resize(m.cellIndexNum());
        for(std::size_t k=0;k<p.vertices.size();++k){const auto& v=p.vertices[k];const auto q=point(m,v.vertex);
            maxMixed=std::max(maxMixed,close(result.vertexMixedDerivatives()[k],field(q).dxy,"symmetry Q4 reproduction"));
            patch.nodes[v.vertex]={v.phi,v.dx,v.dy,result.vertexMixedDerivatives()[k]};}
        for(const auto& e:p.edges)patch.edgeMidpoints.push_back({e.vertices,e.phi,e.dn});
        for(int c:p.geometryCells)patch.cellCenters[c]=field(m.cells(m.cellIndexToLabel[c]).centerCoord()).value;
        const auto scalar=LevelSetHermiteGeometry::Snapshot::fromNativeLeaves(m,patch,radial);
        for(const auto& plane:planes){const auto coordinate=[&](Vector2D q){return plane.normalAxis==0?q.x:q.y;};
            for(const auto& v:result.skeleton().vertices())if(coordinate(point(m,v.index))==plane.coordinate){const auto j=result.skeleton().nodeJet(v.index);
                check(j[plane.normalAxis+1]==0 && j[3]==0,"boundary constrained normal/mixed not exact zero");}
            for(int c:p.geometryCells){const auto b=scalar.nativeCellBounds(c);if(coordinate(b[0])!=plane.coordinate && coordinate(b[1])!=plane.coordinate)continue;
                const auto dofs=scalar.nativeQuarticLocalDofs(c);const int side=coordinate(b[0])==plane.coordinate?2:3;
                for(int k=0;k<5;++k)check(dofs[plane.normalAxis==0?5*k+side:5*side+k]==0,"whole Q4 boundary normal not identically zero");}}
        if(result.report().constrainedEdges)maxTrace=std::max(maxTrace,c1(m,result.skeleton()));unchanged(p,candidate.primary());
        check(candidate.identity()==id && !candidate.hasGeometry() && result.primaryIdentity()==id,"symmetry closure changed primary revision");++cases;
    };
    // Independent nonplanar Q4 fields, even about each tested outer plane.
    for(bool shifted:{false,true})for(int pattern:{0,1,3})for(int axis=0;axis<2;++axis)for(int side=0;side<2;++side){
        auto m=makeMesh(4,pattern,false,shifted);const double wall=axis==0?(side?m->maxCoord.x:m->minCoord.x):(side?m->maxCoord.y:m->minCoord.y);
        const Profile field=[=](Vector2D p){const double n=(axis==0?p.x:p.y)-wall,t=axis==0?p.y:p.x;
            const double a=1+n*n+.125*n*n*n*n,an=2*n+.5*n*n*n,b=1+.3*t+.1*t*t*t*t,bt=.3+.4*t*t*t;
            return axis==0?Jet{a*b,an*b,a*bt,an*bt}:Jet{a*b,a*bt,an*b,an*bt};};
        verify(*m,field,false,{{axis,wall}});
    }
    for(bool radial:{false,true})for(int pattern:{0,1,3}){
        auto m=makeMesh(4,pattern,radial);const Profile field=[](Vector2D p){
            const auto g=[](double t){return 1+t*t*(t-1)*(t-1);};const auto dg=[](double t){return 2*t*(t-1)*(2*t-1);};
            return Jet{g(p.x)*g(p.y),dg(p.x)*g(p.y),g(p.x)*dg(p.y),dg(p.x)*dg(p.y)};};
        std::vector<Plane> planes{{0,1},{1,0},{1,1}};if(!radial)planes.push_back({0,0});verify(*m,field,radial,planes);
    }
    auto m=makeMesh(8,0);const Profile one=[](Vector2D){return Jet{1,0,0,0};};const auto p=primaryFor(*m,one);
    const auto reject=[&](const QuadtreeMesh& mesh,const HermiteGalsPrimaryState& primary,const std::vector<Plane>& planes,const char* why){
        HermiteGalsCandidate candidate(mesh,{},primary);const auto id=candidate.identity();
        rejects([&]{closeHermiteGalsMixedAmr(mesh,candidate,planes);},why);unchanged(primary,candidate.primary());
        check(candidate.identity()==id && !candidate.hasGeometry(),"rejected symmetry polluted primary");++guards;};
    reject(*m,p,{{2,0}},"invalid symmetry plane");reject(*m,p,{{0,std::numeric_limits<double>::quiet_NaN()}},"invalid symmetry plane");
    reject(*m,p,{{0,.5}},"actual outer domain section");reject(*m,p,{{0,0},{0,0}},"duplicate symmetry plane");
    auto bad=p;for(auto& v:bad.vertices)if(point(*m,v.vertex).x==0){v.dx=1e-30;break;}
    reject(*m,bad,{{0,0}},"exactly zero primary vertex normal derivative");
    bad=p;for(auto& e:bad.edges)if(point(*m,e.vertices[0]).x==0 && point(*m,e.vertices[1]).x==0){e.dn=1e-30;break;}
    reject(*m,bad,{{0,0}},"exactly zero primary edge normal derivative");
    const auto active=subset(*m,p,[](Vector2D q){return q.x>.25;});reject(*m,active,{{0,0}},"no active boundary trace");
    auto rz=makeMesh(4,1,true);reject(*rz,primaryFor(*rz,one,true),{{0,0}},"radial axis already has its own symmetry identity");
    check(cases==30 && guards==8,"incomplete declared symmetry coverage");
    std::cout<<"AMR_DECLARED_SYMMETRY,PASS,cases="<<cases<<",rejects="<<guards<<",max_mixed="<<maxMixed<<",max_c1="<<maxTrace<<'\n';
}
}
int main()
{
    try{std::cout<<std::setprecision(17);reproduction();arbitraryAndSensitivity();activeSupport();transport(false);transport(true);failures();convergence(false);convergence(true);declaredSymmetry();}
    catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
    std::cout<<"I06b AMR compact closure/master-only transport contract PASS; no phase/bubble/geometry publication\n";return 0;
}
