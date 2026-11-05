#include "Modules/Interface/HermiteGalsAmrClosure.h"

#include "Modules/Interface/HermiteGalsMixedClosure.h"
#include "Modules/Interface/HermiteGalsVertexTransport.h"
#include "Modules/Interface/HermiteGalsEdgeTransport.h"
#include "Core/Geometry/ActiveQuarticHermitePatch.h"
#include "Core/Geometry/ActiveQuarticHermiteCompletion.h"
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <cstdint>
#include <cfenv>
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
using namespace LAMP::Interface;
using Skeleton=LAMP::Geometry::ActiveQuarticHermiteSkeleton;
using Input=LAMP::Geometry::ActiveQuarticHermiteSkeletonInput;
using Sensitivity=LAMP::Geometry::ActiveQuarticHermiteSensitivity;
using Patch=LAMP::Geometry::ActiveQuarticHermitePatch;
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
void declaredSymmetry(bool quick=false)
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
    if(quick){
        auto m=makeMesh(4,0);const Profile field=[](Vector2D p){
            const auto g=[](double t){return 1+t*t*(t-1)*(t-1);};const auto dg=[](double t){return 2*t*(t-1)*(2*t-1);};
            return Jet{g(p.x)*g(p.y),dg(p.x)*g(p.y),g(p.x)*dg(p.y),dg(p.x)*dg(p.y)};};
        verify(*m,field,false,{{0,0}});
        check(cases==1 && guards==0,"incomplete short declared symmetry coverage");
        std::cout<<"AMR_DECLARED_SYMMETRY_SHORT,PASS,cases="<<cases<<",max_mixed="<<maxMixed<<",max_c1="<<maxTrace<<'\n';return;
    }
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
HermiteGalsState sensitivityState(const QuadtreeMesh& mesh,const HermiteGalsPrimaryState& primary,
    HermiteGalsTimeLayer time={},const std::vector<LevelSetInterfacePanels::SymmetryPlane>& planes={})
{
    HermiteGalsCandidate candidate(mesh,time,primary);
    const auto mixed=closeHermiteGalsMixedAmr(mesh,candidate,planes);
    HermiteGalsClosure closure;closure.vertexMixedDerivatives=mixed.vertexMixedDerivatives();
    Patch patch;patch.activeCells.assign(mesh.cellIndexNum(),0);patch.nodes.resize(mesh.vertIndexNum());
    patch.cellCenters.resize(mesh.cellIndexNum());
    for(const auto& v:mixed.skeleton().vertices()){const auto jet=mixed.skeleton().nodeJet(v.index);
        patch.nodes[v.index]={jet[0],jet[1],jet[2],jet[3]};}
    for(const auto& edge:primary.edges)patch.edgeMidpoints.push_back({edge.vertices,edge.phi,edge.dn});
    for(int c:primary.geometryCells){patch.activeCells[c]=1;patch.cellCenters[c]=1;closure.cellBubbleCoefficients.push_back(1);}
    auto snapshot=LevelSetHermiteGeometry::Snapshot::fromNativeLeaves(mesh,patch,primary.measure==HermiteGalsMeasure::AxisymmetricVolume);
    candidate.bindGeometry(mesh,candidate.identity(),closure,std::move(snapshot));return std::move(candidate).finish();
}
void mixedSensitivityContract()
{
    using Component=HermiteGalsPrimaryDerivative;using Plane=LevelSetInterfacePanels::SymmetryPlane;
    constexpr const char* owner="HermiteGalsAmrClosureSensitivity primary=";
    static_assert(!std::is_default_constructible_v<HermiteGalsAmrClosureSensitivity>);
    const HermiteGalsTimeLayer clock{9,.2,.3,HermiteGalsTimeLayer::Position::Midpoint};
    std::size_t polynomialRows=0,perturbedRows=0,zeroRows=0,edgeTerms=0,guards=0;double peak=0;
    const auto firstJet=[](HermiteGalsPrimaryState primary){primary.transport=HermiteGalsTransport::FirstJet;primary.phaseVolumes.clear();return primary;};
    const auto value=[](const HermiteGalsPrimaryState& primary,const HermiteGalsPrimaryDerivativeWeight& term){
        if(term.component==Component::EdgeNormal)return primary.edges.at(term.primaryOrdinal).dn;
        const auto& vertex=primary.vertices.at(term.primaryOrdinal);return term.component==Component::VertexDx?vertex.dx:vertex.dy;
    };
    for(bool shifted:{false,true}){
        auto mesh=makeMesh(4,0,false,shifted);
        for(int term:{0,6,9,17,24}){const Profile field=[=](Vector2D p){return q4(p,term);};
            const auto primary=firstJet(primaryFor(*mesh,field));const auto source=sensitivityState(*mesh,primary,clock);
            const auto weights=HermiteGalsAmrClosureSensitivity::buildUniformFirstJet(*mesh,source);
            weights.validateFor(*mesh,source);weights.validateFor(*mesh,source.identity(),clock);
            check(weights.primaryIdentity()==source.identity()&&weights.rows().size()==primary.vertices.size()
                &&weights.symmetryPlanes().empty(),"sensitivity source/row binding");
            for(std::size_t row=0;row<weights.rows().size();++row){long double actual=0;
                check(weights.rows()[row].size()<=10,"uniform five-point sensitivity lost sparse support");
                for(const auto& weight:weights.rows()[row]){check(std::isfinite(weight.weight),"nonfinite sensitivity coefficient");
                    actual+=static_cast<long double>(weight.weight)*value(primary,weight);edgeTerms+=weight.component==Component::EdgeNormal;}
                peak=std::max({peak,close(static_cast<double>(actual),source.closure().vertexMixedDerivatives[row],"sparse sensitivity actual compact closure",1e-12),
                    close(static_cast<double>(actual),field(point(*mesh,primary.vertices[row].vertex)).dxy,"sparse sensitivity Q4 reproduction")});++polynomialRows;}
            unchanged(primary,source.primary());
        }
    }
    auto mesh=makeMesh(4,0,false,true);const auto primary=firstJet(primaryFor(*mesh,[](Vector2D p){return q4(p);}));
    const auto source=sensitivityState(*mesh,primary,clock);
    const auto weights=HermiteGalsAmrClosureSensitivity::buildUniformFirstJet(*mesh,source);
    for(Component component:{Component::VertexDx,Component::VertexDy,Component::EdgeNormal}){
        std::vector<std::size_t> ordinals{0,primary.vertices.size()/2,primary.vertices.size()-1};
        if(component==Component::EdgeNormal){ordinals.clear();
            for(bool h:{false,true})for(std::size_t k=0;k<primary.edges.size();++k)if(horizontal(*mesh,primary.edges[k].vertices)==h){ordinals.push_back(k);break;}}
        for(std::size_t ordinal:ordinals){auto changed=primary;constexpr double delta=.125;
            if(component==Component::EdgeNormal)changed.edges[ordinal].dn+=delta;
            else if(component==Component::VertexDx)changed.vertices[ordinal].dx+=delta;else changed.vertices[ordinal].dy+=delta;
            HermiteGalsCandidate perturbed(*mesh,clock,changed);const auto actual=closeHermiteGalsMixedAmr(*mesh,perturbed);
            for(std::size_t row=0;row<weights.rows().size();++row){long double predicted=0;
                for(const auto& weight:weights.rows()[row])if(weight.component==component&&weight.primaryOrdinal==ordinal)predicted+=delta*weight.weight;
                peak=std::max(peak,close(actual.vertexMixedDerivatives()[row]-source.closure().vertexMixedDerivatives[row],
                    static_cast<double>(predicted),"sparse sensitivity independent primary perturbation",1e-12));++perturbedRows;}
        }
    }
    const auto foreign=sensitivityState(*mesh,primary,clock);
    rejects([&]{weights.validateFor(*mesh,foreign);},"different primary revision",owner);++guards;
    for(int change=0;change<4;++change){auto wrong=clock;
        if(change==0)++wrong.step;if(change==1)wrong.begin+=.01;if(change==2)wrong.end+=.01;
        if(change==3)wrong.position=HermiteGalsTimeLayer::Position::End;
        rejects([&]{weights.validateFor(*mesh,source.identity(),wrong);},"different full time layer",owner);++guards;}
    auto changedMesh=mesh->clone();changedMesh->verts(changedMesh->vertIndexToLabel[0]).setCoord({-.1,0});
    rejects([&]{weights.validateFor(*changedMesh,source.identity(),clock);},"mesh",owner);++guards;
    auto symmetryMesh=makeMesh(4,0);const Profile even=[](Vector2D p){
        const auto g=[](double t){return 1+t*t*(t-1)*(t-1);};const auto dg=[](double t){return 2*t*(t-1)*(2*t-1);};
        return Jet{g(p.x)*g(p.y),dg(p.x)*g(p.y),g(p.x)*dg(p.y),dg(p.x)*dg(p.y)};};
    const auto symmetricPrimary=firstJet(primaryFor(*symmetryMesh,even));const std::vector<Plane> planes{{0,0},{0,1},{1,0},{1,1}};
    const auto symmetric=sensitivityState(*symmetryMesh,symmetricPrimary,clock,planes);
    const auto symmetricWeights=HermiteGalsAmrClosureSensitivity::buildUniformFirstJet(*symmetryMesh,symmetric,planes);
    check(symmetricWeights.symmetryPlanes()==planes,"sensitivity lost declared normal-rate constraints");
    for(std::size_t row=0;row<symmetricPrimary.vertices.size();++row){const auto p=point(*symmetryMesh,symmetricPrimary.vertices[row].vertex);
        long double predicted=0;for(const auto& weight:symmetricWeights.rows()[row])predicted+=static_cast<long double>(weight.weight)*value(symmetricPrimary,weight);
        peak=std::max(peak,close(static_cast<double>(predicted),symmetric.closure().vertexMixedDerivatives[row],"symmetry sensitivity actual closure",1e-12));
        if(p.x==0||p.x==1||p.y==0||p.y==1){check(symmetricWeights.rows()[row].empty()
            &&symmetric.closure().vertexMixedDerivatives[row]==0,"symmetry sensitivity row not exact zero");++zeroRows;}}
    for(bool edge:{false,true}){auto bad=symmetricPrimary;
        if(edge){for(auto& e:bad.edges)if(point(*symmetryMesh,e.vertices[0]).x==0&&point(*symmetryMesh,e.vertices[1]).x==0){e.dn=1e-30;break;}}
        else bad.vertices[0].dx=1e-30;
        const auto badSource=sensitivityState(*symmetryMesh,bad,clock);
        rejects([&]{HermiteGalsAmrClosureSensitivity::buildUniformFirstJet(*symmetryMesh,badSource,planes);},
            edge?"exactly zero primary edge normal derivative":"exactly zero primary vertex normal derivative",owner);++guards;
    }
    auto rz=makeMesh(4,0,true);const auto radial=sensitivityState(*rz,firstJet(primaryFor(*rz,[](Vector2D p){return positive(p,true);},true)),clock);
    rejects([&]{HermiteGalsAmrClosureSensitivity::buildUniformFirstJet(*rz,radial);},"RZ unsupported",owner);++guards;
    auto adaptive=makeMesh(4,1);const auto amr=sensitivityState(*adaptive,firstJet(primaryFor(*adaptive,[](Vector2D p){return positive(p);})),clock);
    rejects([&]{HermiteGalsAmrClosureSensitivity::buildUniformFirstJet(*adaptive,amr);},"same-level",owner);++guards;
    const auto partial=sensitivityState(*mesh,subset(*mesh,primary,[](Vector2D p){return p.x<1.4;}),clock);
    rejects([&]{HermiteGalsAmrClosureSensitivity::buildUniformFirstJet(*mesh,partial);},"complete uniform geometry cells",owner);++guards;
    const auto historical=std::move(bound(*symmetryMesh)).finish();
    rejects([&]{HermiteGalsAmrClosureSensitivity::buildUniformFirstJet(*symmetryMesh,historical);},"current FirstJet",owner);++guards;
    unchanged(primary,source.primary());weights.validateFor(*mesh,source);
    check(polynomialRows>0&&perturbedRows>0&&zeroRows>0&&edgeTerms>0&&guards==12,"incomplete sparse closure sensitivity coverage");
    std::cout<<"AMR_MIXED_SENSITIVITY,"<<polynomialRows<<','<<perturbedRows<<','<<zeroRows<<','<<edgeTerms<<','<<guards<<','<<peak<<",PASS\n";
}


// Completion A exercises the neutral E operator only. It neither creates a
// macro snapshot nor substitutes a typed Interface revision/full-clock owner.
using Completion=LAMP::Geometry::ActiveQuarticHermiteCompletion;
using CompletionSource=std::tuple<int,int,int>;
static_assert(!std::is_default_constructible_v<Completion>);
static_assert(std::is_same_v<decltype(std::declval<Completion&>().sourcePatch()),const Patch&>);
static_assert(std::is_same_v<decltype(std::declval<Completion&>().rows()),const std::vector<Completion::Row>&>);
CompletionSource sourceKey(Sensitivity::Kind kind,int index,int component)
{return {static_cast<int>(kind),index,component};}
CompletionSource sourceKey(const Sensitivity& source)
{return sourceKey(source.kind,source.index,source.component);}
long double formalSource(const Patch& patch,const CompletionSource& source)
{
    const auto [kind,index,component]=source;
    if(kind==static_cast<int>(Sensitivity::Kind::Vertex)){
        const auto& j=patch.nodes.at(index);
        const std::array<double,4> v{{j.value,j.dx,j.dy,j.dxy}};
        return v.at(component);
    }
    if(kind==static_cast<int>(Sensitivity::Kind::EdgeMidpoint)){
        const auto& edge=patch.edgeMidpoints.at(index);
        return component==0?edge.value:edge.normalDerivative;
    }
    return patch.cellCenters.at(index);
}
void addFormal(Patch& patch,const CompletionSource& source,double value)
{
    const auto [kind,index,component]=source;
    if(kind==static_cast<int>(Sensitivity::Kind::Vertex)){
        auto& j=patch.nodes.at(index);
        if(component==0)j.value+=value;if(component==1)j.dx+=value;
        if(component==2)j.dy+=value;if(component==3)j.dxy+=value;
    }else if(kind==static_cast<int>(Sensitivity::Kind::EdgeMidpoint)){
        auto& edge=patch.edgeMidpoints.at(index);
        if(component==0)edge.value+=value;else edge.normalDerivative+=value;
    }else patch.cellCenters.at(index)+=value;
}
std::vector<CompletionSource> formalSources(const Patch& patch)
{
    std::vector<CompletionSource> result;
    for(int v=0;v<static_cast<int>(patch.nodes.size());++v)
        for(int k=0;k<4;++k)result.push_back(sourceKey(Sensitivity::Kind::Vertex,v,k));
    for(int e=0;e<static_cast<int>(patch.edgeMidpoints.size());++e)
        for(int k=0;k<2;++k)result.push_back(sourceKey(Sensitivity::Kind::EdgeMidpoint,e,k));
    for(int c=0;c<static_cast<int>(patch.cellCenters.size());++c)
        result.push_back(sourceKey(Sensitivity::Kind::CellCenter,c,0));
    return result;
}
std::unique_ptr<QuadtreeMesh> completionMesh(int nx,int ny,bool shifted=true)
{
    DictionaryMesh d;d.Nx=nx;d.Ny=ny;d.maxCellNum=1024;d.maxInitialIterations=0;
    d.minCoord=shifted?Vector2D{.13,-.27}:Vector2D{0,0};
    d.maxCoord=shifted?Vector2D{1.87,1.36}:Vector2D{1.7,1.5};
    return std::make_unique<QuadtreeMesh>(d);
}
Patch completionPatch(const QuadtreeMesh& mesh,bool arbitrary)
{
    Patch patch;patch.activeCells.assign(mesh.cellIndexNum(),1);
    patch.nodes.resize(mesh.vertIndexNum());patch.cellCenters.resize(mesh.cellIndexNum());
    for(int v=0;v<mesh.vertIndexNum();++v)patch.nodes[v]=q4(point(mesh,v));
    for(const auto& key:actualSides(mesh)){
        const auto j=q4(point(mesh,key[0])+(point(mesh,key[1])-point(mesh,key[0]))*.5);
        patch.edgeMidpoints.push_back({key,j.value,horizontal(mesh,key)?j.dy:j.dx});
    }
    // Source ordinal is deliberately unrelated to canonical side/face order.
    std::reverse(patch.edgeMidpoints.begin(),patch.edgeMidpoints.end());
    for(int c=0;c<mesh.cellIndexNum();++c)
        patch.cellCenters[c]=q4(mesh.cells(mesh.cellIndexToLabel[c]).centerCoord()).value;
    if(arbitrary){
        std::uint32_t seed=287919;
        const auto random=[&]{seed=1664525u*seed+1013904223u;
            return static_cast<double>(seed)/4294967296.0-.5;};
        for(const auto& source:formalSources(patch)){
            addFormal(patch,source,static_cast<double>(-formalSource(patch,source)));
            addFormal(patch,source,random());
        }
    }
    return patch;
}
struct CompletionReference
{
    const QuadtreeMesh& mesh;
    const Patch& patch;
    int nx,ny;
    double dx,dy;
    using Lattice=std::pair<int,int>;
    std::map<Lattice,CompletionSource> values,gx,gy,mixed;
    CompletionReference(const QuadtreeMesh& m,const Patch& p,int x,int y)
        :mesh(m),patch(p),nx(x),ny(y),
         dx((m.maxCoord.x-m.minCoord.x)/(2*x)),
         dy((m.maxCoord.y-m.minCoord.y)/(2*y))
    {
        const auto insert=[&](auto& table,Vector2D position,CompletionSource source){
            const int i=static_cast<int>(std::llround((position.x-m.minCoord.x)/dx));
            const int j=static_cast<int>(std::llround((position.y-m.minCoord.y)/dy));
            close(position.x,m.minCoord.x+i*dx,"reference half-grid x position",1e-13);
            close(position.y,m.minCoord.y+j*dy,"reference half-grid y position",1e-13);
            check(table.emplace(Lattice{i,j},source).second,"duplicate independent half-grid source");
        };
        for(int v=0;v<m.vertIndexNum();++v){
            const auto x0=point(m,v);insert(values,x0,sourceKey(Sensitivity::Kind::Vertex,v,0));
            insert(gx,x0,sourceKey(Sensitivity::Kind::Vertex,v,1));
            insert(gy,x0,sourceKey(Sensitivity::Kind::Vertex,v,2));
            insert(mixed,x0,sourceKey(Sensitivity::Kind::Vertex,v,3));
        }
        for(int e=0;e<static_cast<int>(p.edgeMidpoints.size());++e){
            const auto key=p.edgeMidpoints[e].vertices;
            const auto x0=point(m,key[0])+(point(m,key[1])-point(m,key[0]))*.5;
            insert(values,x0,sourceKey(Sensitivity::Kind::EdgeMidpoint,e,0));
            insert(horizontal(m,key)?gy:gx,x0,sourceKey(Sensitivity::Kind::EdgeMidpoint,e,1));
        }
        for(int c=0;c<m.cellIndexNum();++c)
            insert(values,m.cells(m.cellIndexToLabel[c]).centerCoord(),sourceKey(Sensitivity::Kind::CellCenter,c,0));
        check(values.size()==static_cast<std::size_t>((2*nx+1)*(2*ny+1)),
            "incomplete independent value half-grid");
    }
    // Independent moment-system solve in physical coordinates. Choose nearest
    // support by distance, then order its coordinates; no production FD table.
    std::vector<std::pair<int,long double>> weights(bool x,int target) const
    {
        const int n=2*(x?nx:ny)+1;const long double h=x?dx:dy;
        std::vector<int> candidates;for(int i=0;i<n;++i)candidates.push_back(i);
        std::sort(candidates.begin(),candidates.end(),[=](int a,int b){
            return std::make_pair(std::abs(a-target),a)<std::make_pair(std::abs(b-target),b);});
        candidates.resize(5);std::sort(candidates.begin(),candidates.end());
        const long double scale=h*(candidates.back()-candidates.front());
        long double a[5][6]{};
        for(int column=0;column<5;++column){
            const long double t=h*(candidates[column]-target)/scale;a[0][column]=1;
            for(int row=1;row<5;++row)a[row][column]=a[row-1][column]*t;
        }
        a[2][5]=2/(scale*scale);
        for(int k=0;k<5;++k){
            int pivot=k;for(int i=k+1;i<5;++i)if(std::abs(a[i][k])>std::abs(a[pivot][k]))pivot=i;
            check(a[pivot][k]!=0,"singular independent completion moment system");
            for(int j=k;j<6;++j)std::swap(a[pivot][j],a[k][j]);
            const auto d=a[k][k];for(int j=k;j<6;++j)a[k][j]/=d;
            for(int i=0;i<5;++i)if(i!=k){const auto f=a[i][k];
                for(int j=k;j<6;++j)a[i][j]-=f*a[k][j];}
        }
        std::vector<std::pair<int,long double>> result;
        for(int i=0;i<5;++i)result.push_back({candidates[i],a[i][5]});
        return result;
    }
    long double at(const std::map<Lattice,CompletionSource>& table,int i,int j) const
    {return formalSource(patch,table.at({i,j}));}
    long double derivative(const std::map<Lattice,CompletionSource>& table,bool x,int i,int j) const
    {
        long double result=0;for(const auto& w:weights(x,x?i:j))
            result+=w.second*at(table,x?w.first:i,x?j:w.first);
        return result;
    }
    std::array<long double,49> cell(int c) const
    {
        const auto lo=point(mesh,mesh.verts(mesh.cells(mesh.cellIndexToLabel[c]).cellToVert[0]).index());
        const int baseX=static_cast<int>(std::llround((lo.x-mesh.minCoord.x)/dx));
        const int baseY=static_cast<int>(std::llround((lo.y-mesh.minCoord.y)/dy));
        constexpr int location[7]{0,1,2,0,2,0,2},order[7]{0,0,0,1,1,2,2};
        std::array<long double,49> result{};
        for(int y=0;y<7;++y)for(int x=0;x<7;++x){
            const int i=baseX+location[x],j=baseY+location[y],a=order[x],b=order[y];
            auto& value=result[x+7*y];
            if(a==0 && b==0)value=at(values,i,j);
            else if(a==1 && b==0)value=at(gx,i,j);
            else if(a==0 && b==1)value=at(gy,i,j);
            else if(a==1 && b==1)value=at(mixed,i,j);
            else if(a==2 && b==0)value=derivative(values,true,i,j);
            else if(a==0 && b==2)value=derivative(values,false,i,j);
            else if(a==2 && b==1)value=derivative(gy,true,i,j);
            else if(a==1 && b==2)value=derivative(gx,false,i,j);
            else {
                // One tensor contraction of original V, with no rounded
                // intermediate D2 values or regenerated first derivative.
                for(const auto& wx:weights(true,i))for(const auto& wy:weights(false,j))
                    value+=wx.second*wy.second*at(values,wx.first,wy.first);
            }
        }
        return result;
    }
};
long double completionPolynomialDerivative(Vector2D p,int dx,int dy)
{
    long double result=0;
    for(int a=dx;a<5;++a)for(int b=dy;b<5;++b){
        long double coefficient=((a+b)%2?-1:1)*(a+2*b+1)/50.L;
        for(int k=0;k<dx;++k)coefficient*=a-k;
        for(int k=0;k<dy;++k)coefficient*=b-k;
        result+=coefficient*std::pow(static_cast<long double>(p.x),a-dx)
            *std::pow(static_cast<long double>(p.y),b-dy);
    }
    return result;
}
template<class F> void completionRejects(F action)
{
    try{action();}catch(const std::exception& e){
        check(std::string(e.what()).find("ActiveQuarticHermiteCompletion")!=std::string::npos,
            "completion rejection lacks responsible owner");return;}
    throw std::runtime_error("missing completion preparation/owner rejection");
}
std::map<CompletionSource,long double> completionWeights(const std::vector<Sensitivity>& weights)
{
    std::map<CompletionSource,long double> result;
    for(const auto& w:weights){
        check(std::isfinite(w.coefficient)&&w.coefficient!=0,"invalid completion sparse weight");
        check(result.emplace(sourceKey(w),w.coefficient).second,"duplicate unmerged completion source");
    }
    return result;
}
void c2CompletionContract()
{
    std::size_t cells=0,newRows=0,perturbations=0,transposeCases=0,nonlocalCenters=0,sharedReferences=0,guards=0,zeroRows=0;
    double maximum=0;
    for(const auto size:{std::array<int,2>{{2,2}},std::array<int,2>{{4,3}}}){
        const int nx=size[0],ny=size[1];auto mesh=completionMesh(nx,ny);
        const auto patch=completionPatch(*mesh,true);
        const auto completion=Completion::fromUniformPatch(*mesh,patch);
        completion.validateFor(*mesh,patch);
        check(completion.symmetryPlanes().empty(),"completion invented scalar symmetry");
        check(completion.rows().size()==9*patch.nodes.size()+3*patch.edgeMidpoints.size()+patch.cellCenters.size(),
            "completion rows lack unique shared ownership");
        std::map<CompletionSource,std::size_t> outputRows;std::set<std::size_t> usedRows;
        for(std::size_t r=0;r<completion.rows().size();++r){
            const auto& row=completion.rows()[r];check(outputRows.emplace(sourceKey(row.kind,row.index,row.component),r).second,
                "completion repeated shared output identity");
            completionWeights(row.sourceWeights);
            long double predicted=0;for(const auto& weight:row.sourceWeights)
                predicted+=weight.coefficient*formalSource(patch,sourceKey(weight));
            maximum=std::max(maximum,close(row.value,static_cast<double>(predicted),"completion retained sparse value",1e-12));
        }
        const CompletionReference reference(*mesh,patch,nx,ny);
        std::vector<long double> rowWeights(completion.rows().size(),0);
        for(int c=0;c<mesh->cellIndexNum();++c){
            const auto expected=reference.cell(c);const auto actual=completion.physicalCellDofs(c);
            const auto& indices=completion.cellRows(c);std::array<long double,49> lambda{};
            const auto& cell=mesh->cells(mesh->cellIndexToLabel[c]);const auto center=cell.centerCoord();
            for(int k=0;k<49;++k){
                check(indices[k]<completion.rows().size(),"cell references foreign completion row");
                const auto& row=completion.rows()[indices[k]];
                check(actual[k]==row.value,"cell copied a different shared completion value");
                if(!usedRows.insert(indices[k]).second)++sharedReferences;
                maximum=std::max(maximum,close(actual[k],static_cast<double>(expected[k]),
                    "independent arbitrary completion datum",3e-11));
                lambda[k]=((17*k+13*c)%23-11)/13.L;rowWeights[indices[k]]+=lambda[k];
                if(k%7>=5 || k/7>=5){
                    ++newRows;for(const auto& source:row.sourceWeights)
                        if(source.kind==Sensitivity::Kind::CellCenter){
                            const auto q=mesh->cells(mesh->cellIndexToLabel[source.index]).centerCoord();
                            if(q.x!=center.x || q.y!=center.y)++nonlocalCenters;
                        }
                }else {
                    check(row.sourceWeights.size()==1 && row.sourceWeights[0].coefficient==1,
                        "completion changed original 25 selectors");
                }
            }
            const auto transpose=completionWeights(completion.pullbackCell(c,lambda));
            Patch rate=patch;for(const auto& source:formalSources(rate))
                addFormal(rate,source,static_cast<double>(-formalSource(rate,source)));
            for(const auto& source:formalSources(rate)){
                const auto [kind,index,component]=source;
                addFormal(rate,source,((kind*19+index*7+component*11)%29-14)/17.);
            }
            const auto direct=CompletionReference(*mesh,rate,nx,ny).cell(c);
            long double left=0,right=0;for(int k=0;k<49;++k)left+=lambda[k]*direct[k];
            for(const auto& weight:transpose)right+=weight.second*formalSource(rate,weight.first);
            maximum=std::max(maximum,close(static_cast<double>(right),static_cast<double>(left),
                "independent completion cell transpose",3e-11));++transposeCases;++cells;
        }
        check(usedRows.size()==completion.rows().size(),"completion retains unused shared output rows");
        Patch rate=patch;for(const auto& source:formalSources(rate))
            addFormal(rate,source,static_cast<double>(-formalSource(rate,source)));
        for(const auto& source:formalSources(rate)){
            const auto [kind,index,component]=source;
            addFormal(rate,source,((kind*11+index*13+component*3)%31-15)/19.);
        }
        const auto rates=completion.applyRates(rate);
        check(rates.size()==completion.rows().size(),"completion rate layout differs from value layout");
        const auto transpose=completionWeights(completion.pullbackRows(rowWeights));
        long double left=0,right=0;
        for(int c=0;c<mesh->cellIndexNum();++c){
            const auto direct=CompletionReference(*mesh,rate,nx,ny).cell(c);
            const auto& indices=completion.cellRows(c);
            for(int k=0;k<49;++k)maximum=std::max(maximum,close(rates[indices[k]],
                static_cast<double>(direct[k]),"independent completion applied source rate",3e-11));
        }
        for(std::size_t r=0;r<rates.size();++r)left+=rowWeights[r]*rates[r];
        for(const auto& weight:transpose)right+=weight.second*formalSource(rate,weight.first);
        maximum=std::max(maximum,close(static_cast<double>(right),static_cast<double>(left),
            "completion shared-row transpose",3e-11));++transposeCases;
        const auto sources=formalSources(patch);
        for(std::size_t ordinal=0;ordinal<sources.size();++ordinal){
            // Exhaust all 64 minimal-mesh sources; larger support only needs a
            // representative per-kind/component set in addition to full E/rate.
            if(nx!=2 && ordinal%11!=0 && ordinal+1!=sources.size())continue;
            auto changed=patch;addFormal(changed,sources[ordinal],.125);
            const auto perturbed=Completion::fromUniformPatch(*mesh,changed);
            Patch delta=patch;for(const auto& source:sources)
                addFormal(delta,source,static_cast<double>(-formalSource(delta,source)));
            addFormal(delta,sources[ordinal],static_cast<double>(formalSource(changed,sources[ordinal])-formalSource(patch,sources[ordinal])));
            for(int c=0;c<mesh->cellIndexNum();++c){
                const auto direct=CompletionReference(*mesh,delta,nx,ny).cell(c);
                const auto before=completion.physicalCellDofs(c),after=perturbed.physicalCellDofs(c);
                const auto& indices=completion.cellRows(c);
                for(int k=0;k<49;++k){
                    long double sparse=0;for(const auto& weight:completion.rows()[indices[k]].sourceWeights)
                        if(sourceKey(weight)==sources[ordinal])sparse+=weight.coefficient*formalSource(delta,sources[ordinal]);
                    maximum=std::max({maximum,close(after[k]-before[k],static_cast<double>(direct[k]),
                        "independent completion unit source perturbation",3e-11),
                        close(static_cast<double>(sparse),static_cast<double>(direct[k]),
                        "completion sparse source perturbation",3e-11)});++perturbations;
                }
            }
        }
        auto ownedInput=patch;const auto frozen=Completion::fromUniformPatch(*mesh,ownedInput);
        const auto frozenCell=frozen.physicalCellDofs(0);
        ownedInput.nodes[0].value+=.125;ownedInput.edgeMidpoints[0].value+=.125;
        ownedInput.cellCenters.back()+=.125;
        for(const auto& source:formalSources(patch))check(formalSource(frozen.sourcePatch(),source)==formalSource(patch,source),
            "immutable completion borrows mutable formal source");
        check(frozen.physicalCellDofs(0)==frozenCell,"mutated caller data changes immutable completion");
        auto changed=patch;changed.cellCenters.back()+=.125;
        completionRejects([&]{completion.validateFor(*mesh,changed);});++guards;
        changed=patch;std::swap(changed.edgeMidpoints[0],changed.edgeMidpoints[1]);
        completionRejects([&]{completion.validateFor(*mesh,changed);});++guards;
        auto foreign=mesh->clone();foreign->verts(foreign->vertIndexToLabel[0]).setCoord({-.19,-.27});
        completionRejects([&]{completion.validateFor(*foreign,patch);});++guards;
        completionRejects([&]{completion.validateFor(*mesh,patch,{{0,mesh->minCoord.x}});});++guards;
        auto shortRates=rate;shortRates.cellCenters.pop_back();
        completionRejects([&]{completion.applyRates(shortRates);});++guards;
    }
    // Tensor Q4 has independent analytic high derivative data at exact sites.
    auto mesh=completionMesh(4,3);const auto polynomial=completionPatch(*mesh,false);
    const auto completed=Completion::fromUniformPatch(*mesh,polynomial);
    constexpr int location[7]{0,1,2,0,2,0,2},order[7]{0,0,0,1,1,2,2};
    for(int c=0;c<mesh->cellIndexNum();++c){
        const auto& cell=mesh->cells(mesh->cellIndexToLabel[c]);
        const auto lo=point(*mesh,mesh->verts(cell.cellToVert[0]).index());
        const auto hi=point(*mesh,mesh->verts(cell.cellToVert[2]).index());const auto values=completed.physicalCellDofs(c);
        for(int y=0;y<7;++y)for(int x=0;x<7;++x){
            const Vector2D q{lo.x+(hi.x-lo.x)*location[x]*.5,lo.y+(hi.y-lo.y)*location[y]*.5};
            maximum=std::max(maximum,close(values[x+7*y],static_cast<double>(
                completionPolynomialDerivative(q,order[x],order[y])),"completion analytic tensor Q4 reproduction",3e-11));
        }
    }
    auto symmetricMesh=completionMesh(2,2,false);auto symmetric=completionPatch(*symmetricMesh,true);
    using Plane=Completion::SymmetryPlane;
    const std::vector<Plane> planes{{0,symmetricMesh->minCoord.x},{0,symmetricMesh->maxCoord.x},
        {1,symmetricMesh->minCoord.y},{1,symmetricMesh->maxCoord.y}};
    for(int v=0;v<symmetricMesh->vertIndexNum();++v){const auto q=point(*symmetricMesh,v);
        const bool x=q.x==symmetricMesh->minCoord.x || q.x==symmetricMesh->maxCoord.x;
        const bool y=q.y==symmetricMesh->minCoord.y || q.y==symmetricMesh->maxCoord.y;
        if(x)symmetric.nodes[v].dx=0;if(y)symmetric.nodes[v].dy=0;if(x||y)symmetric.nodes[v].dxy=0;
    }
    for(auto& edge:symmetric.edgeMidpoints){const auto q=point(*symmetricMesh,edge.vertices[0]);
        if(horizontal(*symmetricMesh,edge.vertices)
            ?q.y==symmetricMesh->minCoord.y || q.y==symmetricMesh->maxCoord.y
            :q.x==symmetricMesh->minCoord.x || q.x==symmetricMesh->maxCoord.x)edge.normalDerivative=0;
    }
    const auto symmetry=Completion::fromUniformPatch(*symmetricMesh,symmetric,planes);
    symmetry.validateFor(*symmetricMesh,symmetric,planes);
    check(symmetry.symmetryPlanes()==planes,"completion lost explicit scalar symmetry");
    for(const auto& row:symmetry.rows()){
        if(row.kind!=Sensitivity::Kind::Vertex)continue;
        const auto q=point(*symmetricMesh,row.index);
        const bool oddX=row.component==7 && (q.x==symmetricMesh->minCoord.x || q.x==symmetricMesh->maxCoord.x);
        const bool oddY=row.component==6 && (q.y==symmetricMesh->minCoord.y || q.y==symmetricMesh->maxCoord.y);
        if(oddX||oddY){check(row.symmetryConstrained && row.value==0 && row.sourceWeights.empty(),
            "completion odd high derivative is not an exact zero row");++zeroRows;}
    }
    const auto symmetricRates=symmetry.applyRates(symmetric);
    for(std::size_t r=0;r<symmetry.rows().size();++r)if(symmetry.rows()[r].symmetryConstrained
        &&symmetry.rows()[r].sourceWeights.empty())check(symmetricRates[r]==0,"constrained completion rate not exact zero");
    auto bad=symmetric;bad.nodes[0].dx=1e-30;
    completionRejects([&]{Completion::fromUniformPatch(*symmetricMesh,bad,planes);});++guards;
    completionRejects([&]{symmetry.applyRates(bad);});++guards;
    bad=symmetric;bad.nodes[0].dxy=1e-30;
    completionRejects([&]{Completion::fromUniformPatch(*symmetricMesh,bad,planes);});++guards;
    bad=symmetric;
    for(auto& edge:bad.edgeMidpoints)if(!horizontal(*symmetricMesh,edge.vertices)
        &&point(*symmetricMesh,edge.vertices[0]).x==symmetricMesh->minCoord.x){edge.normalDerivative=1e-30;break;}
    completionRejects([&]{Completion::fromUniformPatch(*symmetricMesh,bad,planes);});++guards;
    auto invalid=polynomial;invalid.activeCells.back()=0;
    completionRejects([&]{Completion::fromUniformPatch(*mesh,invalid);});++guards;
    invalid=polynomial;invalid.cellCenters[0]=std::numeric_limits<double>::quiet_NaN();
    completionRejects([&]{Completion::fromUniformPatch(*mesh,invalid);});++guards;
    invalid=polynomial;invalid.edgeMidpoints.push_back(invalid.edgeMidpoints.front());
    completionRejects([&]{Completion::fromUniformPatch(*mesh,invalid);});++guards;
    invalid=polynomial;invalid.cornerSectors.push_back({});
    completionRejects([&]{Completion::fromUniformPatch(*mesh,invalid);});++guards;
    auto insufficient=completionMesh(1,2);
    const auto shortPatch=completionPatch(*insufficient,true);
    completionRejects([&]{Completion::fromUniformPatch(*insufficient,shortPatch);});++guards;
    auto adaptive=makeMesh(4,1);
    const auto adaptivePatch=completionPatch(*adaptive,true);
    completionRejects([&]{Completion::fromUniformPatch(*adaptive,adaptivePatch);});++guards;
    check(cells==16&&newRows==24*cells&&perturbations>12544&&transposeCases==18
        &&nonlocalCenters>0&&sharedReferences>0&&guards==20&&zeroRows==12,
        "incomplete completion/sensitivity/constraint coverage");
    std::cout<<"C2_COMPLETION_A,PASS,cells="<<cells<<",new24="<<newRows
        <<",perturbed_rows="<<perturbations<<",transpose="<<transposeCases
        <<",neighbor_centers="<<nonlocalCenters<<",shared_refs="<<sharedReferences
        <<",zero_rows="<<zeroRows<<",rejects="<<guards<<",max_absolute="<<maximum<<'\n';
}

void uncoupledClosureOnly()
{
    for(bool radial:{false,true})for(int pattern:{0,1}){
        auto m=makeMesh(4,pattern,radial,true);const Profile field=[=](Vector2D p){return q4(p,-1,radial);};
        const auto p=primaryFor(*m,field,radial);const HermiteGalsTimeLayer clock{9,.2,.3,HermiteGalsTimeLayer::Position::Midpoint};
        HermiteGalsCandidate candidate(*m,clock,p);const auto id=candidate.identity();
        const auto result=closeHermiteGalsMixedAmr(*m,candidate);double mixed=0,skeleton=0,native=0;std::size_t axis=0;
        check(result.report().masterVertices==p.vertices.size() && result.report().masterEdges==p.edges.size(),"short closure master inventory");
        check(pattern?result.report().hangingVertices>0 && result.report().constrainedEdges>0:
            result.report().hangingVertices==0 && result.report().constrainedEdges==0,"short closure topology route witness");
        Patch patch;patch.activeCells.assign(m->cellIndexNum(),1);patch.nodes.resize(m->vertIndexNum());patch.cellCenters.resize(m->cellIndexNum());
        for(std::size_t k=0;k<p.vertices.size();++k){const auto& v=p.vertices[k];const auto jet=result.skeleton().nodeJet(v.vertex);
            check(jet[0]==v.phi && jet[1]==v.dx && jet[2]==v.dy,"short closure changed original master first jet");
            mixed=std::max(mixed,close(result.vertexMixedDerivatives()[k],field(point(*m,v.vertex)).dxy,"short closure analytic mixed"));
            patch.nodes[v.vertex]={v.phi,v.dx,v.dy,result.vertexMixedDerivatives()[k]};}
        for(const auto& e:p.edges){const auto midpoint=result.skeleton().edgeMidpoint(e.vertices);
            check(midpoint[0]==e.phi && midpoint[1]==e.dn,"short closure changed original master midpoint");
            patch.edgeMidpoints.push_back({e.vertices,e.phi,e.dn});}
        for(const auto& v:result.skeleton().vertices()){const auto actual=result.skeleton().nodeJet(v.index);const auto exact=field(point(*m,v.index));
            skeleton=std::max({skeleton,close(actual[0],exact.value,"short skeleton analytic phi"),close(actual[1],exact.dx,"short skeleton analytic dx"),
                close(actual[2],exact.dy,"short skeleton analytic dy"),close(actual[3],exact.dxy,"short skeleton analytic mixed")});
            if(radial && point(*m,v.index).x==0){check(actual[1]==0 && actual[3]==0,"short closure lost exact axis parity");++axis;}}
        if(radial)check(axis>0,"short closure missed actual radial axis");
        for(int c:p.geometryCells)patch.cellCenters[c]=field(m->cells(m->cellIndexToLabel[c]).centerCoord()).value;
        const auto scalar=LevelSetHermiteGeometry::Snapshot::fromNativeLeaves(*m,patch,radial);
        for(int c:p.geometryCells){const auto bounds=scalar.nativeCellBounds(c);const auto width=bounds[1]-bounds[0];
            const auto center=m->cells(m->cellIndexToLabel[c]).centerCoord();
            close(scalar.scalarJetInNativeCell(c,center).value,patch.cellCenters[c],"short native Q4 original center");
            for(const Vector2D x:{center,bounds[0]+Vector2D{.3125*width.x,.6875*width.y}}){const auto actual=scalar.scalarJetInNativeCell(c,x);const auto exact=field(x);
                native=std::max({native,close(actual.value,exact.value,"short final native Q4 phi"),close(actual.dx,exact.dx,"short final native Q4 dx"),
                    close(actual.dy,exact.dy,"short final native Q4 dy"),close(actual.dxy,exact.dxy,"short final native Q4 mixed")});}}
        const double compatibility=directCompatibility(*m,p,result),trace=pattern?c1(*m,result.skeleton()):0;
        check(result.report().algebraicResidual<=1e-12 && result.report().compatibilityResidual<=1e-12,"short closure original residual gates");
        unchanged(p,candidate.primary());check(candidate.identity()==id && result.primaryIdentity()==id && !candidate.hasGeometry()
            && result.timeLayer().step==clock.step && result.timeLayer().time()==clock.time(),"short closure source/clock isolation");
        std::cout<<"AMR_UNCOUPLED_Q4,"<<(radial?"RZ":"XY")<<','<<pattern<<','<<m->cellIndexNum()<<','<<mixed<<','<<skeleton<<','<<native
            <<','<<compatibility<<','<<trace<<','<<result.report().compatibilityNonzeros<<",PASS\n";
    }
    declaredSymmetry(true);failures();
    auto m=makeMesh(4,0);const auto p=primaryFor(*m,[](Vector2D x){return q4(x);});HermiteGalsCandidate candidate(*m,{},p);const auto id=candidate.identity();
    const int rounding=std::fegetround();check(rounding==FE_TONEAREST,"short closure requires initial nearest rounding");
    check(std::fesetround(FE_DOWNWARD)==0 && std::fegetround()==FE_DOWNWARD,"short closure could not set actual downward rounding");
    try{rejects([&]{closeHermiteGalsMixedAmr(*m,candidate);},"round-to-nearest");}
    catch(...){std::fesetround(rounding);throw;}
    check(std::fesetround(rounding)==0 && std::fegetround()==rounding,"short closure could not restore rounding");
    unchanged(p,candidate.primary());check(candidate.identity()==id && !candidate.hasGeometry(),"FP rejected closure polluted source");
    std::cout<<"AMR_UNCOUPLED_FP,ACTUAL_DOWNWARD_REJECTED,SOURCE_UNCHANGED,PASS\n";
}
}
int main(int argc,char** argv)
{
    try{std::cout<<std::setprecision(17);
        if(argc==2 && std::string(argv[1])=="--c2-completion-only"){c2CompletionContract();return 0;}
        if(argc==2 && std::string(argv[1])=="--mixed-sensitivity-only"){mixedSensitivityContract();return 0;}
        if(argc==2 && std::string(argv[1])=="--uncoupled-closure-only"){uncoupledClosureOnly();std::cout<<"I06b uncoupled/general closure short contract PASS\n";return 0;}
        reproduction();arbitraryAndSensitivity();activeSupport();transport(false);transport(true);failures();convergence(false);convergence(true);declaredSymmetry();mixedSensitivityContract();c2CompletionContract();}
    catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
    std::cout<<"I06b AMR compact closure/master-only transport contract PASS; no phase/bubble/geometry publication\n";return 0;
}
