#include "Core/Geometry/ActiveQuarticHermiteSkeleton.h"

#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <cmath>
#include <future>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

namespace
{
using Skeleton=LAMP::Geometry::ActiveQuarticHermiteSkeleton;
using Input=LAMP::Geometry::ActiveQuarticHermiteSkeletonInput;
using Patch=LAMP::Geometry::ActiveQuarticHermitePatch;
using Key=Skeleton::EdgeKey;
using Jet=Skeleton::Jet;
constexpr Key none{{-1,-1}};
constexpr double tolerance=3e-11;
static_assert(!std::is_default_constructible_v<Skeleton>);
static_assert(std::is_same_v<decltype(std::declval<Skeleton&>().vertices()),const std::vector<Skeleton::Vertex>&>);
static_assert(std::is_same_v<decltype(std::declval<Skeleton&>().edges()),const std::vector<Skeleton::Edge>&>);
void check(bool ok,const char* reason){if(!ok)throw std::runtime_error(reason);}
double difference(double a,double b)
{
    if(!std::isfinite(a) || !std::isfinite(b) || std::abs(a-b)>tolerance*(1+std::abs(b)))
        throw std::runtime_error("skeleton mismatch: actual="+std::to_string(a)+" expected="+std::to_string(b));
    return std::abs(a-b);
}
double difference(const Jet& a,const Jet& b)
{double e=0;for(int k=0;k<4;++k)e=std::max(e,difference(a[k],b[k]));return e;}
template<class F> void rejects(F action,const char* message)
{
    try{action();}catch(const std::exception& error){
        const std::string s=error.what();
        check(s.find("ActiveQuarticHermiteSkeleton")!=std::string::npos,"missing skeleton failure context");
        if(s.find(message)==std::string::npos)throw std::runtime_error("unexpected rejection: "+s);
        return;
    }
    throw std::runtime_error(std::string("missing rejection: ")+message);
}
Vector2D point(const QuadtreeMesh& mesh,int v){return mesh.verts(mesh.vertIndexToLabel[v]).coord();}
bool horizontal(const QuadtreeMesh& mesh,Key k){return point(mesh,k[0]).y==point(mesh,k[1]).y;}
Vector2D middle(const QuadtreeMesh& mesh,Key k){const auto a=point(mesh,k[0]);return a+(point(mesh,k[1])-a)*.5;}
std::unique_ptr<QuadtreeMesh> makeMesh(int pattern,bool radial=false,bool shifted=false)
{
    DictionaryMesh d;d.Nx=d.Ny=4;d.maxCellNum=2048;d.maxInitialIterations=0;
    if(shifted){d.minCoord={radial?0:.13,-.27};d.maxCoord={1.87,1.36};}
    auto mesh=std::make_unique<QuadtreeMesh>(d);
    if(pattern){
        const auto lo=d.minCoord,span=d.maxCoord-lo;
        DictionaryMesh::RefineBox b;
        b.min={lo.x+.26*span.x,lo.y+.26*span.y};
        const double upper=pattern==1?.49:.74;
        b.max={lo.x+upper*span.x,lo.y+upper*span.y};b.level=1;
        if(radial){b.min.x=lo.x;b.max.x=lo.x+.24*span.x;}
        std::vector<DictionaryMesh::RefineBox> boxes{b};
        if(pattern==3){
            if(radial)boxes[0].max.x=lo.x+.49*span.x;
            for(int level=2;level<=3;++level){
                const double low=level==2?.38:.45,high=level==2?.62:.55;
                b.min={radial?lo.x:lo.x+low*span.x,lo.y+low*span.y};
                b.max={lo.x+(radial?(level==2?.24:.12):high)*span.x,lo.y+high*span.y};
                b.level=level;boxes.push_back(b);
            }
        }
        mesh->applyRefineBoxes(boxes);
        int maxLevel=0;
        for(int c=0;c<mesh->cellIndexNum();++c)maxLevel=std::max(maxLevel,mesh->cells(mesh->cellIndexToLabel[c]).level());
        check(maxLevel==(pattern==3?3:1),"refinement fixture did not reach requested level");
    }
    return mesh;
}
// Independent fixture topology: enumerate all sides then test geometric strict
// containment, instead of using the production skeleton's parent classification.
struct Topology
{
    std::map<Key,Key> parent;
    std::map<int,Key> vertexParent;
    std::set<int> used;
    explicit Topology(const QuadtreeMesh& mesh,const std::vector<std::uint8_t>& active)
    {
        for(int c=0;c<mesh.cellIndexNum();++c)if(active[c]){
            const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);
            for(int k=0;k<4;++k){
                const auto& a=mesh.verts(cell.cellToVert[k]);const auto& b=mesh.verts(cell.cellToVert[(k+1)%4]);
                Key key{{a.index(),b.index()}};
                if(a.coord().y==b.coord().y?a.coord().x>b.coord().x:a.coord().y>b.coord().y)std::swap(key[0],key[1]);
                parent.emplace(key,none);used.insert(a.index());
            }
        }
        for(auto& e:parent){const auto a=point(mesh,e.first[0]),b=point(mesh,e.first[1]);
            for(const auto& other:parent)if(other.first!=e.first){const auto lo=point(mesh,other.first[0]),hi=point(mesh,other.first[1]);
                if((a.y==b.y && lo.y==hi.y && a.y==lo.y && a.x>=lo.x && b.x<=hi.x) ||
                   (a.x==b.x && lo.x==hi.x && a.x==lo.x && a.y>=lo.y && b.y<=hi.y))e.second=other.first;
            }
        }
        for(int v:used){Key parentSide=none;const auto p=point(mesh,v);
            if(mesh.verts(mesh.vertIndexToLabel[v]).type()!=TVertType::non)
                for(const auto& e:parent){const auto a=point(mesh,e.first[0]),b=point(mesh,e.first[1]);
                    if((a.y==b.y && p.y==a.y && a.x<p.x && p.x<b.x) ||
                       (a.x==b.x && p.x==a.x && a.y<p.y && p.y<b.y))parentSide=e.first;}
            vertexParent.emplace(v,parentSide);
        }
    }
};
Jet polynomial(Vector2D p,int term,bool radial)
{
    Jet result{};
    for(int i=0;i<5;++i)for(int j=0;j<5;++j){
        if(radial && i%2)continue;
        if(term>=0 && term!=5*i+j)continue;
        const double a=term>=0?1:((i+j)%2?-1:1)*(i+2*j+1)/50.0;
        const double x=std::pow(p.x,i),y=std::pow(p.y,j);
        const double dx=i?i*std::pow(p.x,i-1):0,dy=j?j*std::pow(p.y,j-1):0;
        result[0]+=a*x*y;result[1]+=a*dx*y;result[2]+=a*x*dy;result[3]+=a*dx*dy;
    }
    return result;
}
template<class Field> Input analyticInput(const QuadtreeMesh& mesh,const Field& field,bool allSides=true)
{
    Input input;input.activeCells.assign(mesh.cellIndexNum(),1);input.nodes.resize(mesh.vertIndexNum());
    const Topology topology(mesh,input.activeCells);
    for(int v:topology.used){const auto j=field(point(mesh,v));input.nodes[v]={j[0],j[1],j[2],j[3]};}
    for(const auto& e:topology.parent)if(allSides || e.second==none){const auto j=field(middle(mesh,e.first));
        input.edgeMidpoints.push_back({e.first,j[0],j[horizontal(mesh,e.first)?2:1]});}
    return input;
}
void checkTopology(const Skeleton& skeleton,const Topology& topology)
{
    check(skeleton.vertices().size()==topology.used.size() && skeleton.edges().size()==topology.parent.size(),"skeleton topology size");
    std::size_t masters=0,edges=0;
    for(const auto& v:skeleton.vertices()){check(v.masterSide==topology.vertexParent.at(v.index),"vertex master identity");if(v.masterSide==none)++masters;}
    for(const auto& e:skeleton.edges()){check(e.parentSide==topology.parent.at(e.vertices),"edge parent identity");if(e.parentSide==none)++edges;}
    check(skeleton.independentDofCount()==4*masters+2*edges,"hanging or center entered independent DOF inventory");
}
void q4Reproduction()
{
    for(bool radial:{false,true})for(int pattern=0;pattern<4;++pattern){
        auto mesh=makeMesh(pattern,radial,true);double error=0,traceError=0;std::size_t queries=0,axis=0;int fields=0;
        // All 25 monomials on the isolated refinement, plus a mixed Q4 on
        // every topology (including multilevel and boundary/axis refinement).
        for(int term=-1;term<(pattern==1?25:0);++term){if(radial && term>=0 && (term/5)%2)continue;++fields;
            const auto field=[&](Vector2D p){return polynomial(p,term,radial);};
            const auto input=analyticInput(*mesh,field,false);const Topology topology(*mesh,input.activeCells);
            const Skeleton skeleton(*mesh,input,radial);checkTopology(skeleton,topology);
            for(const auto& v:skeleton.vertices()){const auto p=point(*mesh,v.index);const auto j=skeleton.nodeJet(v.index);
                error=std::max(error,difference(j,field(p)));++queries;
                if(radial && p.x==0){check(j[1]==0 && j[3]==0,"axis nodal parity");++axis;}}
            for(const auto& e:skeleton.edges()){
                const auto a=point(*mesh,e.vertices[0]),b=point(*mesh,e.vertices[1]);const bool h=horizontal(*mesh,e.vertices);
                const auto midpoint=skeleton.edgeMidpoint(e.vertices);const auto exact=field(middle(*mesh,e.vertices));
                error=std::max({error,difference(midpoint[0],exact[0]),difference(midpoint[1],exact[h?2:1])});
                for(double fraction:{0.,.125,.375,.5,.625,.875,1.}){
                    const auto p=fraction==0?a:(fraction==1?b:a+(b-a)*fraction);const double s=h?p.x:p.y;const auto j=skeleton.edgeJet(e.vertices,s);
                    error=std::max(error,difference(j,field(p)));++queries;
                    if(e.parentSide!=none)traceError=std::max(traceError,difference(j,skeleton.edgeJet(e.parentSide,s)));
                    if(radial && p.x==0){check(j[1]==0 && j[3]==0,"axis trace parity");++axis;}
                }
            }
        }
        std::cout<<"SKELETON_Q4,"<<(radial?"RZ":"XY")<<','<<pattern<<','<<mesh->cellIndexNum()<<','<<fields<<','<<queries<<','<<axis<<','<<error<<','<<traceError<<'\n';
    }
}
// Independent power-basis solve for the five Hermite interpolation conditions.
// This uses no production basis, rational restriction or internal engine API.
std::array<double,2> quartic(double f0,double f1,double g0,double g1,double fm,double width,double s)
{
    const long double a=f0,b=static_cast<long double>(width)*g0;
    const long double A=f1-a-b,B=static_cast<long double>(width)*g1-b,C=fm-a-.5L*b;
    const long double c=-5*A+B+16*C,d=14*A-3*B-32*C,e=-8*A+2*B+16*C,t=s;
    return {{static_cast<double>(a+t*(b+t*(c+t*(d+t*e)))),static_cast<double>((b+t*(2*c+t*(3*d+4*t*e)))/width)}};
}
struct Reference
{
    const QuadtreeMesh& mesh;const Input& input;const Topology& topology;
    std::map<Key,Patch::EdgeMidpoint> midpoints;
    Reference(const QuadtreeMesh& m,const Input& i,const Topology& t):mesh(m),input(i),topology(t)
    {for(const auto& e:input.edgeMidpoints)midpoints.emplace(e.vertices,e);}
    Jet node(int v) const
    {
        const auto parent=topology.vertexParent.at(v);
        if(parent!=none){const auto p=point(mesh,v);return edge(parent,horizontal(mesh,parent)?p.x:p.y);}
        const auto& j=input.nodes[v];return {{j.value,j.dx,j.dy,j.dxy}};
    }
    Jet edge(Key key,double coordinate) const
    {
        const auto parent=topology.parent.at(key);if(parent!=none)return edge(parent,coordinate);
        const auto a=point(mesh,key[0]),b=point(mesh,key[1]);const bool h=horizontal(mesh,key);
        const auto lo=node(key[0]),hi=node(key[1]);const auto middleValue=midpoints.at(key);
        const int tangent=h?1:2,normal=h?2:1;const double width=h?b.x-a.x:b.y-a.y;
        const double t=(coordinate-(h?a.x:a.y))/width;
        const auto value=quartic(lo[0],hi[0],lo[tangent],hi[tangent],middleValue.value,width,t);
        const auto dn=quartic(lo[normal],hi[normal],lo[3],hi[3],middleValue.normalDerivative,width,t);
        Jet result{};result[0]=value[0];result[tangent]=value[1];result[normal]=dn[0];result[3]=dn[1];return result;
    }
};
void randomC1AndOwnership(int pattern,bool radial)
{
    auto mesh=makeMesh(pattern,radial);auto input=analyticInput(*mesh,[](Vector2D){return Jet{};});
    const Topology topology(*mesh,input.activeCells);std::uint32_t seed=1234567;
    const auto random=[&](){seed=1664525u*seed+1013904223u;return (static_cast<double>(seed)/4294967296.0-.5)*2;};
    for(int v:topology.used){input.nodes[v]={random(),random(),random(),random()};if(radial && point(*mesh,v).x==0){input.nodes[v].dx=0;input.nodes[v].dxy=0;}}
    for(auto& e:input.edgeMidpoints){e.value=random();e.normalDerivative=random();if(radial && !horizontal(*mesh,e.vertices) && point(*mesh,e.vertices[0]).x==0)e.normalDerivative=0;}
    const Skeleton skeleton(*mesh,input,radial);checkTopology(skeleton,topology);const Reference reference(*mesh,input,topology);
    std::size_t constrained=0,hanging=0,horizontalFine=0,verticalFine=0,hangingEnds=0;double defect=0,referenceError=0,snapshotError=0;
    Patch patch;patch.activeCells=input.activeCells;patch.nodes=input.nodes;patch.edgeMidpoints=input.edgeMidpoints;
    patch.cellCenters.assign(mesh->cellIndexNum(),.125); // test geometry only, arbitrary independent interiors
    const auto snapshot=LevelSetHermiteGeometry::Snapshot::fromNativeLeaves(*mesh,patch,radial);
    for(const auto& v:skeleton.vertices()){referenceError=std::max(referenceError,difference(skeleton.nodeJet(v.index),reference.node(v.index)));if(v.masterSide!=none)++hanging;}
    for(const auto& e:skeleton.edges()){
        const auto a=point(*mesh,e.vertices[0]),b=point(*mesh,e.vertices[1]);const bool h=horizontal(*mesh,e.vertices);
        if(e.parentSide==none && (topology.vertexParent.at(e.vertices[0])!=none || topology.vertexParent.at(e.vertices[1])!=none))++hangingEnds;
        if(e.parentSide!=none){++constrained;(h?horizontalFine:verticalFine)++;}
        for(int i=0;i<=8;++i){const auto p=a+(b-a)*(i/8.);const double s=h?p.x:p.y;const auto jet=skeleton.edgeJet(e.vertices,s);
            referenceError=std::max(referenceError,difference(jet,reference.edge(e.vertices,s)));
            if(e.parentSide!=none)defect=std::max(defect,difference(jet,skeleton.edgeJet(e.parentSide,s)));
        }
    }
    // Old cell evaluator on EVERY side of EVERY cell: the skeleton is the same
    // boundary of the existing Q4 geometry despite arbitrary interior centers.
    for(int c=0;c<mesh->cellIndexNum();++c){const auto& cell=mesh->cells(mesh->cellIndexToLabel[c]);
        for(int k=0;k<4;++k){Key key{{mesh->verts(cell.cellToVert[k]).index(),mesh->verts(cell.cellToVert[(k+1)%4]).index()}};
            const bool h=horizontal(*mesh,key);auto a=point(*mesh,key[0]),b=point(*mesh,key[1]);
            if(h?a.x>b.x:a.y>b.y){std::swap(key[0],key[1]);std::swap(a,b);}
            for(double t:{.125,.5,.875}){const auto p=a+(b-a)*t;const auto actual=snapshot.scalarJetInNativeCell(c,p);
                snapshotError=std::max(snapshotError,difference({{actual.value,actual.dx,actual.dy,actual.dxy}},skeleton.edgeJet(key,h?p.x:p.y)));}
        }
    }
    check(hanging && horizontalFine && verticalFine && hangingEnds,"random fixture misses AMR constraint cases");
    auto ignored=input;const double nan=std::numeric_limits<double>::quiet_NaN();
    for(const auto& v:topology.vertexParent)if(v.second!=none)ignored.nodes[v.first]={nan,nan,nan,nan};
    for(auto& e:ignored.edgeMidpoints)if(topology.parent.at(e.vertices)!=none){e.value=nan;e.normalDerivative=nan;}
    const Skeleton poisoned(*mesh,ignored,radial);
    for(const auto& v:skeleton.vertices())check(poisoned.nodeJet(v.index)==skeleton.nodeJet(v.index),"hanging payload entered result");
    for(const auto& e:skeleton.edges())check(poisoned.edgeMidpoint(e.vertices)==skeleton.edgeMidpoint(e.vertices),"fine payload entered result");
    auto mastersOnly=ignored;
    mastersOnly.edgeMidpoints.erase(std::remove_if(mastersOnly.edgeMidpoints.begin(),mastersOnly.edgeMidpoints.end(),[&](const auto& e){return topology.parent.at(e.vertices)!=none;}),mastersOnly.edgeMidpoints.end());
    const Skeleton minimal(*mesh,mastersOnly,radial);
    for(const auto& v:skeleton.vertices())check(minimal.nodeJet(v.index)==skeleton.nodeJet(v.index),"omitting dependent inputs changed result");
    const auto key=skeleton.edges().front().vertices;const auto p=middle(*mesh,key);const double s=horizontal(*mesh,key)?p.x:p.y;
    const auto expected=skeleton.edgeJet(key,s);const auto copy=skeleton;std::vector<std::future<Jet>> queries;
    for(int i=0;i<8;++i)queries.push_back(std::async(std::launch::async,[copy,key,s]{Jet last{};for(int j=0;j<20;++j)last=copy.edgeJet(key,s);return last;}));
    for(auto& q:queries)check(q.get()==expected,"concurrent query changed immutable trace");
    input.nodes.clear();input.edgeMidpoints.clear();skeleton.validateMesh(*mesh);
    auto changed=mesh->clone();changed->verts(changed->vertIndexToLabel[0]).setCoord({-.5,0});
    rejects([&]{skeleton.validateMesh(*changed);},"mesh");
    const int cells=mesh->cellIndexNum();mesh.reset();check(copy.edgeJet(key,s)==expected,"owner borrowed mutable input/mesh");
    std::cout<<"SKELETON_C1,"<<(radial?"RZ":"XY")<<','<<pattern<<','<<cells<<','<<hanging<<','<<constrained<<','<<horizontalFine<<','<<verticalFine<<','<<hangingEnds<<','<<skeleton.independentDofCount()<<','<<defect<<','<<referenceError<<','<<snapshotError<<'\n';
}
void rejectionAndSubset()
{
    auto mesh=makeMesh(1);const auto field=[](Vector2D p){return polynomial(p,-1,false);};
    const auto input=analyticInput(*mesh,field,false);const Skeleton owner(*mesh,input);
    auto bad=input;bad.nodes.pop_back();rejects([&]{Skeleton x(*mesh,bad);},"full vertex/cell arrays");
    bad=input;bad.activeCells[0]=2;rejects([&]{Skeleton x(*mesh,bad);},"binary");
    bad=input;bad.edgeMidpoints.pop_back();rejects([&]{Skeleton x(*mesh,bad);},"missing independent");
    bad=input;bad.edgeMidpoints.push_back(bad.edgeMidpoints.front());rejects([&]{Skeleton x(*mesh,bad);},"duplicate");
    bad=input;std::swap(bad.edgeMidpoints[0].vertices[0],bad.edgeMidpoints[0].vertices[1]);rejects([&]{Skeleton x(*mesh,bad);},"ordered active leaf side");
    bad=input;bad.edgeMidpoints[0].vertices={{0,0}};rejects([&]{Skeleton x(*mesh,bad);},"ordered active leaf side");
    const double nan=std::numeric_limits<double>::quiet_NaN();
    bad=input;bad.nodes[0].dxy=nan;rejects([&]{Skeleton x(*mesh,bad);},"non-finite");
    bad=input;bad.edgeMidpoints[0].normalDerivative=nan;rejects([&]{Skeleton x(*mesh,bad);},"non-finite");
    bad=input;std::fill(bad.activeCells.begin(),bad.activeCells.end(),0);rejects([&]{Skeleton x(*mesh,bad);},"no active leaves");
    for(const auto& v:owner.vertices())if(v.masterSide!=none){
        const int coarse=mesh->verts(mesh->vertIndexToLabel[v.index]).supCell;bad=input;bad.activeCells[mesh->cells(coarse).index()]=0;
        bad.edgeMidpoints.clear();rejects([&]{Skeleton x(*mesh,bad);},"inactive leaf");break;
    }
    rejects([&]{owner.nodeJet(-1);},"used native vertex");rejects([&]{owner.edgeJet({{0,0}},0);},"ordered active native side");
    const auto edge=owner.edges().front().vertices;
    rejects([&]{owner.edgeJet(edge,-1);},"leaves its actual cell");rejects([&]{owner.edgeJet(edge,nan);},"leaves its actual cell");
    auto rz=makeMesh(2,true);auto regular=analyticInput(*rz,[](Vector2D p){return polynomial(p,-1,true);},false);
    const Skeleton axis(*rz,regular,true);regular.nodes[0].dx=1;rejects([&]{Skeleton x(*rz,regular,true);},"axis dx and dxy");
    regular.nodes[0].dx=0;regular.nodes[0].dxy=1;rejects([&]{Skeleton x(*rz,regular,true);},"axis dx and dxy");regular.nodes[0].dxy=0;
    for(auto& e:regular.edgeMidpoints)if(!horizontal(*rz,e.vertices) && point(*rz,e.vertices[0]).x==0){e.normalDerivative=1;break;}
    rejects([&]{Skeleton x(*rz,regular,true);},"axis midpoint");
    auto shifted=makeMesh(0,false,true);const auto shiftedInput=analyticInput(*shifted,field,false);
    rejects([&]{Skeleton x(*shifted,shiftedInput,true);},"requires r=0");
    auto uniform=makeMesh(0);auto subset=analyticInput(*uniform,field);subset.activeCells.assign(uniform->cellIndexNum(),0);subset.activeCells[0]=1;
    const Topology subTopology(*uniform,subset.activeCells);
    subset.edgeMidpoints.erase(std::remove_if(subset.edgeMidpoints.begin(),subset.edgeMidpoints.end(),[&](const auto& e){return !subTopology.parent.count(e.vertices);}),subset.edgeMidpoints.end());
    for(int v=0;v<uniform->vertIndexNum();++v)if(!subTopology.used.count(v))subset.nodes[v]={nan,nan,nan,nan};
    const Skeleton active(*uniform,subset);checkTopology(active,subTopology);check(active.independentDofCount()==24,"subset unexpectedly includes cell center");
    rejects([&]{active.nodeJet(uniform->vertIndexNum()-1);},"used native vertex");
    for(const auto& v:active.vertices())difference(active.nodeJet(v.index),field(point(*uniform,v.index)));
    // An AMR subset must retain the real coarse neighbors of its hanging
    // vertices, while unrelated leaves and their payloads remain absent.
    auto amrSubset=input;amrSubset.activeCells.assign(mesh->cellIndexNum(),0);
    for(int c=0;c<mesh->cellIndexNum();++c){const auto& cell=mesh->cells(mesh->cellIndexToLabel[c]);
        if(cell.level()==0)continue;amrSubset.activeCells[c]=1;
        for(int label:cell.cellToVert){const int coarse=mesh->verts(label).supCell;
            if(coarse>=0)amrSubset.activeCells[mesh->cells(coarse).index()]=1;}
    }
    const Topology amrTopology(*mesh,amrSubset.activeCells);
    amrSubset.edgeMidpoints.erase(std::remove_if(amrSubset.edgeMidpoints.begin(),amrSubset.edgeMidpoints.end(),
        [&](const auto& e){return !amrTopology.parent.count(e.vertices);}),amrSubset.edgeMidpoints.end());
    const Skeleton amrActive(*mesh,amrSubset);checkTopology(amrActive,amrTopology);
    check(amrActive.vertices().size()<owner.vertices().size(),"AMR subset accidentally covered the full mesh");
    for(const auto& v:amrActive.vertices())check(amrActive.nodeJet(v.index)==owner.nodeJet(v.index),"supported AMR subset changed constraints");
    for(const auto& e:amrActive.edges())check(amrActive.edgeMidpoint(e.vertices)==owner.edgeMidpoint(e.vertices),"supported AMR subset changed fine midpoint");
    auto moved=active;const auto retained=std::move(moved);rejects([&]{moved.nodeJet(0);},"moved-from");check(retained.vertices().size()==4,"move lost owner");
    std::cout<<"SKELETON_REJECTION,full_and_subset_input_identity_axis_and_query_guards,PASS\n";
}
}
int main()
{
    try{std::cout<<std::setprecision(17);q4Reproduction();for(bool radial:{false,true})for(int p:{1,2,3})randomC1AndOwnership(p,radial);rejectionAndSubset();}
    catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
    std::cout<<"I06a exact skeleton constraint API PASS; I06b AMR mixed closure remains\n";return 0;
}
