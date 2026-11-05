#include "Modules/Interface/HermiteGalsAmrClosure.h"

#include "Core/Geometry/CompactQuarticDerivative.h"
#include "Core/Geometry/Detail/ExactExpansion.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/UniformLeafGrid.h"
#include <Eigen/SparseLU>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>

namespace LAMP::Interface
{
namespace
{
using Skeleton=LAMP::Geometry::ActiveQuarticHermiteSkeleton;
using Key=Skeleton::EdgeKey;
using Sensitivity=LAMP::Geometry::ActiveQuarticHermiteSensitivity;
constexpr Key none{{-1,-1}};
constexpr double compatibilityTolerance=1e-12;
void require(bool ok,const char* why){if(!ok)throw std::invalid_argument(why);}
struct Sample
{
    double coordinate=0,constant=0;
    int vertex=-1;Key edge=none;
    std::map<int,double> mixed;
};
struct Direction
{
    std::array<Sample,5> samples;
    int target=-1;
    bool zero=false;
    template<class Value> double derivative(const Value& value) const
    {
        if(zero)return 0;
        return LAMP::Geometry::compactQuarticDerivativeAtNode(target,5,
            [&](int i){return samples[i].coordinate;},value,require);
    }
    double constant() const {return derivative([&](int i){return samples[i].constant;});}
    std::map<int,double> coefficients() const
    {
        std::set<int> columns;for(const auto& sample:samples)for(const auto& c:sample.mixed)columns.insert(c.first);
        std::map<int,double> result;
        for(int column:columns){const double c=derivative([&](int i){const auto found=samples[i].mixed.find(column);
            return found==samples[i].mixed.end()?0:found->second;});if(c!=0)result.emplace(column,c);}
        return result;
    }
    double evaluate(const Skeleton& skeleton,int component) const
    {
        return derivative([&](int i){const auto& s=samples[i];
            return s.vertex>=0?skeleton.nodeJet(s.vertex)[component]:skeleton.edgeMidpoint(s.edge)[1];});
    }
};
// Uncoupled actual topology needs no affine mixed-DOF maps. Retain the
// identical five physical samples and derivative accumulation arithmetic.
struct DirectSample
{
    double coordinate=0,constant=0;
    int vertex=-1;Key edge=none;
};
struct DirectDirection
{
    std::array<DirectSample,5> samples;
    int target=-1;
    bool zero=false;
    template<class Value> double derivative(const Value& value) const
    {
        if(zero)return 0;
        return LAMP::Geometry::compactQuarticDerivativeAtNode(target,5,
            [&](int i){return samples[i].coordinate;},value,require);
    }
    double constant() const {return derivative([&](int i){return samples[i].constant;});}
    double evaluate(const Skeleton& skeleton,int component) const
    {
        return derivative([&](int i){const auto& s=samples[i];
            return s.vertex>=0?skeleton.nodeJet(s.vertex)[component]:skeleton.edgeMidpoint(s.edge)[1];});
    }
};
struct Side {Key key;bool horizontal=false;double line=0,lo=0,hi=0;};
struct DirectClosurePreparation
{
    std::vector<std::array<DirectDirection,2>> directions;
    std::vector<std::array<bool,2>> reflection;
    HermiteGalsAmrClosureReport report;
};
DirectClosurePreparation prepareDirectClosure(const QuadtreeMesh& mesh,const HermiteGalsPrimaryState& primary,
    const std::vector<LevelSetInterfacePanels::SymmetryPlane>& symmetryPlanes,bool axis,
    const std::map<int,int>& index,const std::map<Key,double>& normals,const char*& stage,int& vertex)
{
    const auto point=[&](int v){return mesh.verts(mesh.vertIndexToLabel[v]).coord();};
    std::set<int> usedVertices;std::set<Key> usedSides;std::vector<Side> sides;
    constexpr int cornerOrder[4]={0,1,3,2};
    constexpr int sideCorners[4][2]={{0,1},{2,3},{0,2},{1,3}};
    for(int c:primary.geometryCells){
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);std::array<int,4> vertices;
        for(int k=0;k<4;++k){vertices[k]=mesh.verts(cell.cellToVert[cornerOrder[k]]).index();usedVertices.insert(vertices[k]);}
        for(int side=0;side<4;++side){
            const Key key{{vertices[sideCorners[side][0]],vertices[sideCorners[side][1]]}};
            if(!usedSides.insert(key).second)continue;
            require(normals.count(key)!=0,"missing active master midpoint");
            const auto a=point(key[0]),b=point(key[1]);const bool h=a.y==b.y;
            sides.push_back({key,h,h?a.y:a.x,h?a.x:a.y,h?b.x:b.y});
        }
    }
    for(int v:usedVertices)require(index.count(v)!=0,"missing active master vertex");
    require(usedVertices.size()==primary.vertices.size(),"primary includes foreign or unused master vertices");
    require(usedSides.size()==primary.edges.size(),"primary includes foreign or unused master edges");
    HermiteGalsAmrClosureReport report;
    report.masterVertices=usedVertices.size();report.masterEdges=usedSides.size();
    stage="declared symmetry constraints";
    std::vector<std::array<bool,2>> reflection(primary.vertices.size(),{{false,false}});
    for(std::size_t i=0;i<primary.vertices.size();++i)reflection[i][0]=axis && point(primary.vertices[i].vertex).x==0;
    for(std::size_t k=0;k<symmetryPlanes.size();++k){const auto& plane=symmetryPlanes[k];
        require((plane.normalAxis==0 || plane.normalAxis==1) && std::isfinite(plane.coordinate),"invalid symmetry plane");
        const auto coordinate=[&](Vector2D q){return plane.normalAxis==0?q.x:q.y;};
        require(plane.coordinate==coordinate(mesh.minCoord) || plane.coordinate==coordinate(mesh.maxCoord),
            "symmetry plane must be an actual outer domain section");
        require(!(axis && plane.normalAxis==0 && plane.coordinate==0),"radial axis already has its own symmetry identity");
        for(std::size_t j=0;j<k;++j)require(!(plane==symmetryPlanes[j]),"duplicate symmetry plane");
        bool found=false;
        for(const auto& e:primary.edges)if(coordinate(point(e.vertices[0]))==plane.coordinate && coordinate(point(e.vertices[1]))==plane.coordinate){
            found=true;require(e.dn==0,"symmetry requires exactly zero primary edge normal derivative");}
        require(found,"symmetry plane has no active boundary trace");
        for(std::size_t i=0;i<primary.vertices.size();++i){const auto& v=primary.vertices[i];
            if(coordinate(point(v.vertex))!=plane.coordinate)continue;vertex=v.vertex;
            require((plane.normalAxis==0?v.dx:v.dy)==0,"symmetry requires exactly zero primary vertex normal derivative");
            if(!reflection[i][0] && !reflection[i][1])++report.symmetryConstrainedVertices;
            reflection[i][plane.normalAxis]=true;
        }
    }
    std::sort(sides.begin(),sides.end(),[](const Side& a,const Side& b){
        if(a.horizontal!=b.horizontal)return a.horizontal<b.horizontal;
        if(a.line!=b.line)return a.line<b.line;return a.lo<b.lo;});
    std::vector<std::array<DirectDirection,2>> directions(primary.vertices.size());
    const auto node=[&](int v,bool h){
        DirectSample s;const auto p=point(v);s.coordinate=h?p.x:p.y;s.vertex=v;
        const auto& datum=primary.vertices[index.at(v)];s.constant=h?datum.dy:datum.dx;
        return s;
    };
    stage="connected five-point support";
    for(std::size_t first=0;first<sides.size();){
        std::size_t last=first+1;
        while(last<sides.size() && sides[last].horizontal==sides[first].horizontal && sides[last].line==sides[first].line
            && sides[last].lo==sides[last-1].hi){
            require(sides[last].key[0]==sides[last-1].key[1],"ambiguous connected skeleton endpoint");++last;
        }
        const bool h=sides[first].horizontal;const int direction=h?0:1;
        std::vector<DirectSample> line;line.push_back(node(sides[first].key[0],h));
        for(std::size_t e=first;e<last;++e){const auto& side=sides[e];const double mid=side.lo+(side.hi-side.lo)*.5;
            require(mid>side.lo && mid<side.hi,"unresolved master midpoint coordinate");
            DirectSample s;s.coordinate=mid;s.constant=normals.at(side.key);s.edge=side.key;line.push_back(std::move(s));
            line.push_back(node(side.key[1],h));
        }
        ++report.lineComponents;
        for(int j=0;j<static_cast<int>(line.size());++j){const auto found=index.find(line[j].vertex);if(found==index.end())continue;
            vertex=line[j].vertex;require(line.size()>=5,"connected master skeleton has insufficient five-point support");
            auto& d=directions[found->second][direction];require(d.target<0,"duplicate target on disconnected skeleton");
            const int begin=std::clamp(j-2,0,static_cast<int>(line.size())-5);d.target=j-begin;
            std::copy_n(line.begin()+begin,5,d.samples.begin());
            d.zero=reflection[found->second][direction];
        }
        first=last;
    }
    return {std::move(directions),std::move(reflection),report};
}
}
HermiteGalsAmrClosureResult::HermiteGalsAmrClosureResult(HermiteGalsIdentity primary,HermiteGalsTimeLayer time,
    std::vector<double> mixed,Skeleton skeleton,HermiteGalsAmrClosureReport report)
    :primary_(primary),time_(time),mixed_(std::move(mixed)),skeleton_(std::move(skeleton)),report_(report) {}


HermiteGalsAmrClosureSensitivity::HermiteGalsAmrClosureSensitivity(HermiteGalsState source,
    std::vector<std::vector<HermiteGalsPrimaryDerivativeWeight>> rows,
    std::vector<LevelSetInterfacePanels::SymmetryPlane> symmetry)
    :source_(std::move(source)),primary_(source_.identity()),time_(source_.timeLayer()),
        rows_(std::move(rows)),symmetry_(std::move(symmetry)) {}

void HermiteGalsAmrClosureSensitivity::validateFor(const QuadtreeMesh& mesh,const HermiteGalsState& source) const
{
    validateFor(mesh,source.identity(),source.timeLayer());
}
void HermiteGalsAmrClosureSensitivity::validateFor(const QuadtreeMesh& mesh,HermiteGalsIdentity identity,
    HermiteGalsTimeLayer time) const
{
    try{
        source_.validateMesh(mesh);time.time();
        require(primary_==identity,"different primary revision");
        require(time_.step==time.step&&time_.begin==time.begin&&time_.end==time.end&&time_.position==time.position,
            "different full time layer");
    }catch(const std::exception& error){throw std::runtime_error("[Interface/HermiteGalsAmrClosureSensitivity primary="
        +std::to_string(primary_.value())+"] binding: "+error.what());}
}

HermiteGalsAmrClosureSensitivity HermiteGalsAmrClosureSensitivity::buildUniformFirstJet(const QuadtreeMesh& mesh,
    const HermiteGalsState& source,const std::vector<LevelSetInterfacePanels::SymmetryPlane>& symmetryPlanes)
{
    HermiteGalsIdentity identity;int vertex=-1;const char* stage="input";
    try{
        identity=source.identity();source.validateMesh(mesh);
        const auto& primary=source.primary();
        require(primary.transport==HermiteGalsTransport::FirstJet,"sensitivity requires current FirstJet");
        require(primary.measure==HermiteGalsMeasure::CartesianArea,"sensitivity requires Cartesian XY; RZ unsupported");
        const UniformLeafGrid grid(mesh,CartesianGridScope::UniformLeaves);
        require(grid.nx()>=2&&grid.ny()>=2,"sensitivity requires five-point uniform support");
        require(primary.geometryCells.size()==static_cast<std::size_t>(mesh.cellIndexNum()),
            "sensitivity requires complete uniform geometry cells");
        for(int c=0;c<mesh.cellIndexNum();++c)require(primary.geometryCells[c]==c,
            "sensitivity requires complete uniform geometry cells");
        require(primary.vertices.size()==static_cast<std::size_t>(mesh.vertIndexNum())&&
            primary.edges.size()==static_cast<std::size_t>(mesh.hEdgeIndexNum()+mesh.vEdgeIndexNum()),
            "sensitivity requires complete nonhanging primary skeleton");
        std::map<int,int> index;std::map<Key,double> normals;std::map<Key,std::size_t> edges;
        for(std::size_t i=0;i<primary.vertices.size();++i){const int v=primary.vertices[i].vertex;
            require(mesh.verts(mesh.vertIndexToLabel[v]).supCell<0,"sensitivity does not support hanging vertices");
            index.emplace(v,static_cast<int>(i));}
        for(std::size_t i=0;i<primary.edges.size();++i){const auto& edge=primary.edges[i];
            normals.emplace(edge.vertices,edge.dn);edges.emplace(edge.vertices,i);}
        LAMP::Geometry::Detail::ExactExpansion::validateEnvironment();stage="connected five-point sensitivity";
        const auto direct=prepareDirectClosure(mesh,primary,symmetryPlanes,false,index,normals,stage,vertex);
        std::vector<std::vector<HermiteGalsPrimaryDerivativeWeight>> rows(primary.vertices.size());
        require(source.closure().vertexMixedDerivatives.size()==rows.size(),"different retained mixed closure inventory");
        for(std::size_t i=0;i<rows.size();++i){vertex=primary.vertices[i].vertex;
            const auto& directions=direct.directions[i];
            require(directions[0].target>=0&&directions[1].target>=0,"master target lacks both directional supports");
            const bool constrained=direct.reflection[i][0]||direct.reflection[i][1];
            if(!constrained){
                std::map<std::pair<HermiteGalsPrimaryDerivative,std::size_t>,double> weights;
                for(int direction=0;direction<2;++direction){const auto& d=directions[direction];
                    for(int slot=0;slot<5;++slot){const auto& sample=d.samples[slot];
                        const double weight=.5*d.derivative([&](int j){return j==slot?1.:0.;});
                        if(weight==0)continue;
                        const auto component=sample.vertex>=0?(direction==0?HermiteGalsPrimaryDerivative::VertexDy:
                            HermiteGalsPrimaryDerivative::VertexDx):HermiteGalsPrimaryDerivative::EdgeNormal;
                        const std::size_t ordinal=sample.vertex>=0?static_cast<std::size_t>(index.at(sample.vertex)):edges.at(sample.edge);
                        weights[{component,ordinal}]+=weight;
                    }}
                for(const auto& term:weights)if(term.second!=0)rows[i].push_back({term.first.first,term.first.second,term.second});
            }
            const double expected=constrained?0:.5*directions[0].constant()+.5*directions[1].constant();
            const double actual=source.closure().vertexMixedDerivatives[i];
            const double scale=std::max({1.,std::abs(expected),std::abs(actual)});
            require(std::isfinite(expected)&&std::isfinite(actual)&&std::abs(expected/scale-actual/scale)<=compatibilityTolerance,
                "sensitivity differs from retained same-revision compact closure");
        }
        source.validateMesh(mesh);LAMP::Geometry::Detail::ExactExpansion::validateEnvironment();
        return HermiteGalsAmrClosureSensitivity(source,std::move(rows),symmetryPlanes);
    }catch(const std::exception& error){throw std::runtime_error("[Interface/HermiteGalsAmrClosureSensitivity primary="
        +std::to_string(identity.value())+" vertex="+std::to_string(vertex)+"] "+stage+": "+error.what());}
}

HermiteGalsAmrClosureResult closeHermiteGalsMixedAmr(const QuadtreeMesh& mesh,const HermiteGalsCandidate& candidate,
    const std::vector<LevelSetInterfacePanels::SymmetryPlane>& symmetryPlanes)
{
    HermiteGalsIdentity identity;int vertex=-1;const char* stage="input";
    try{
        identity=candidate.identity();candidate.validateMesh(mesh);
        require(!candidate.hasGeometry(),"AMR closure requires an unpublished primary revision");
        const auto& primary=candidate.primary();
        require(!primary.geometryCells.empty() && !primary.vertices.empty(),"AMR closure requires active cells and master first jets");
        const auto point=[&](int v){return mesh.verts(mesh.vertIndexToLabel[v]).coord();};
        const bool axis=primary.measure==HermiteGalsMeasure::AxisymmetricVolume && mesh.minCoord.x==0;
        LAMP::Geometry::ActiveQuarticHermiteSkeletonInput input;
        input.activeCells.assign(mesh.cellIndexNum(),0);for(int c:primary.geometryCells)input.activeCells[c]=1;
        const double nan=std::numeric_limits<double>::quiet_NaN();input.nodes.assign(mesh.vertIndexNum(),{nan,nan,nan,nan});
        std::map<int,int> index;std::map<Key,double> normals;
        for(std::size_t i=0;i<primary.vertices.size();++i){const auto& v=primary.vertices[i];
            input.nodes[v.vertex]={v.phi,v.dx,v.dy,0};index.emplace(v.vertex,static_cast<int>(i));}
        for(const auto& e:primary.edges){input.edgeMidpoints.push_back({e.vertices,e.phi,e.dn});normals.emplace(e.vertices,e.dn);}
        // Core assigns a hanging constraint only from an actually used
        // vertex's supCell, and an edge parent only from such a vertex's
        // coarseEdge. With none of these vertices, the mixed operator is I.
        bool uncoupled=true;
        for(int c:primary.geometryCells){
            for(int label:mesh.cells(mesh.cellIndexToLabel[c]).cellToVert)
                if(mesh.verts(label).supCell>=0){uncoupled=false;break;}
            if(!uncoupled)break;
        }
        if(uncoupled){
            // The removed origin constructor checked this environment before
            // any compact arithmetic; retain the same entry check here.
            LAMP::Geometry::Detail::ExactExpansion::validateEnvironment();
            stage="affine constraint assembly";
            auto direct=prepareDirectClosure(mesh,primary,symmetryPlanes,axis,index,normals,stage,vertex);
            auto& report=direct.report;const auto& directions=direct.directions;const auto& reflection=direct.reflection;
            const int count=static_cast<int>(primary.vertices.size());std::vector<double> mixed(count);
            stage="compact compatibility system";
            for(int i=0;i<count;++i){vertex=primary.vertices[i].vertex;const auto& d=directions[i];
                require(d[0].target>=0 && d[1].target>=0,"master target lacks both directional supports");
                mixed[i]=reflection[i][0] || reflection[i][1]?0:.5*d[0].constant()+.5*d[1].constant();
                require(std::isfinite(mixed[i]),"nonfinite compatibility right-hand side");
                input.nodes[vertex].dxy=mixed[i];
            }
            // I has exactly one nonzero per row. solution=rhs gives the exact
            // algebraic residual I*rhs-rhs=0; no matrix or solver was executed.
            report.compatibilityNonzeros=static_cast<std::size_t>(count);report.algebraicResidual=0;
            stage="final exact constraints";const Skeleton skeleton(mesh,input,axis);
            for(const auto& plane:symmetryPlanes)for(const auto& v:skeleton.vertices()){
                const auto q=point(v.index);if((plane.normalAxis==0?q.x:q.y)!=plane.coordinate)continue;
                const auto jet=skeleton.nodeJet(v.index);
                require(jet[plane.normalAxis+1]==0 && jet[3]==0,"final skeleton violates declared symmetry");
            }
            double compatibility=0;
            for(int i=0;i<count;++i){vertex=primary.vertices[i].vertex;
                const double dx=directions[i][0].evaluate(skeleton,2),dy=directions[i][1].evaluate(skeleton,1);
                const double expected=.5*dx+.5*dy;require(std::isfinite(expected),"nonfinite direct compact closure");
                const double scale=std::max({1.0,std::abs(dx),std::abs(dy),std::abs(mixed[i])});
                compatibility=std::max(compatibility,std::abs(mixed[i]/scale-expected/scale));
                const double directionalScale=std::max({1.0,std::abs(dx),std::abs(dy)});
                report.maxRelativeDirectionalDisagreement=std::max(report.maxRelativeDirectionalDisagreement,
                    2*std::abs(.5*dx/directionalScale-.5*dy/directionalScale));
            }
            require(compatibility<=compatibilityTolerance,"direct compact/constraint residual exceeds tolerance");report.compatibilityResidual=compatibility;
            stage="final mesh binding";vertex=-1;candidate.validateMesh(mesh);
            LAMP::Geometry::Detail::ExactExpansion::validateEnvironment();
            return HermiteGalsAmrClosureResult(identity,candidate.timeLayer(),std::move(mixed),skeleton,report);
        }
        stage="affine constraint assembly";
        // q=0 is ONLY the origin of the affine constraint map a+Cq. This
        // private operator is never returned as the completed closure.
        const Skeleton origin(mesh,input,axis);HermiteGalsAmrClosureReport report;
        for(const auto& v:origin.vertices()){
            if(v.masterSide==none){require(index.count(v.index)!=0,"missing active master vertex");++report.masterVertices;}
            else ++report.hangingVertices;
        }
        require(report.masterVertices==primary.vertices.size(),"primary includes foreign or unused master vertices");
        std::vector<Side> sides;
        for(const auto& e:origin.edges()){
            if(e.parentSide!=none){++report.constrainedEdges;require(!normals.count(e.vertices),"constrained fine edge cannot be a primary master");continue;}
            require(normals.count(e.vertices)!=0,"missing active master midpoint");++report.masterEdges;
            const auto a=point(e.vertices[0]),b=point(e.vertices[1]);const bool h=a.y==b.y;
            sides.push_back({e.vertices,h,h?a.y:a.x,h?a.x:a.y,h?b.x:b.y});
        }
        require(report.masterEdges==primary.edges.size(),"primary includes foreign or unused master edges");
        stage="declared symmetry constraints";
        std::vector<std::array<bool,2>> reflection(primary.vertices.size(),{{false,false}});
        for(std::size_t i=0;i<primary.vertices.size();++i)reflection[i][0]=axis && point(primary.vertices[i].vertex).x==0;
        for(std::size_t k=0;k<symmetryPlanes.size();++k){const auto& plane=symmetryPlanes[k];
            require((plane.normalAxis==0 || plane.normalAxis==1) && std::isfinite(plane.coordinate),"invalid symmetry plane");
            const auto coordinate=[&](Vector2D q){return plane.normalAxis==0?q.x:q.y;};
            require(plane.coordinate==coordinate(mesh.minCoord) || plane.coordinate==coordinate(mesh.maxCoord),
                "symmetry plane must be an actual outer domain section");
            require(!(axis && plane.normalAxis==0 && plane.coordinate==0),"radial axis already has its own symmetry identity");
            for(std::size_t j=0;j<k;++j)require(!(plane==symmetryPlanes[j]),"duplicate symmetry plane");
            bool found=false;
            for(const auto& e:primary.edges)if(coordinate(point(e.vertices[0]))==plane.coordinate && coordinate(point(e.vertices[1]))==plane.coordinate){
                found=true;require(e.dn==0,"symmetry requires exactly zero primary edge normal derivative");}
            require(found,"symmetry plane has no active boundary trace");
            for(std::size_t i=0;i<primary.vertices.size();++i){const auto& v=primary.vertices[i];
                if(coordinate(point(v.vertex))!=plane.coordinate)continue;vertex=v.vertex;
                require((plane.normalAxis==0?v.dx:v.dy)==0,"symmetry requires exactly zero primary vertex normal derivative");
                if(!reflection[i][0] && !reflection[i][1])++report.symmetryConstrainedVertices;
                reflection[i][plane.normalAxis]=true;
            }
        }
        std::sort(sides.begin(),sides.end(),[](const Side& a,const Side& b){
            if(a.horizontal!=b.horizontal)return a.horizontal<b.horizontal;
            if(a.line!=b.line)return a.line<b.line;return a.lo<b.lo;});
        std::vector<std::array<Direction,2>> directions(primary.vertices.size());
        const auto node=[&](int v,bool h){
            Sample s;const auto p=point(v);s.coordinate=h?p.x:p.y;s.vertex=v;const int component=h?2:1;
            s.constant=origin.nodeJet(v)[component];
            for(const auto& c:origin.nodeSensitivity(v,component))if(c.kind==Sensitivity::Kind::Vertex && c.component==3){
                const auto found=index.find(c.index);require(found!=index.end(),"constraint depends on absent primary master");
                const double coefficient=static_cast<double>(c.coefficient);
                require(std::isfinite(coefficient) && (c.coefficient==0 || coefficient!=0),"unrepresentable mixed constraint coefficient");
                if(coefficient!=0)s.mixed.emplace(found->second,coefficient);
            }
            return s;
        };
        stage="connected five-point support";
        for(std::size_t first=0;first<sides.size();){
            std::size_t last=first+1;
            while(last<sides.size() && sides[last].horizontal==sides[first].horizontal && sides[last].line==sides[first].line
                && sides[last].lo==sides[last-1].hi){
                require(sides[last].key[0]==sides[last-1].key[1],"ambiguous connected skeleton endpoint");++last;
            }
            const bool h=sides[first].horizontal;const int direction=h?0:1;
            std::vector<Sample> line;line.push_back(node(sides[first].key[0],h));
            for(std::size_t e=first;e<last;++e){const auto& side=sides[e];const double mid=side.lo+(side.hi-side.lo)*.5;
                require(mid>side.lo && mid<side.hi,"unresolved master midpoint coordinate");
                Sample s;s.coordinate=mid;s.constant=normals.at(side.key);s.edge=side.key;line.push_back(std::move(s));
                line.push_back(node(side.key[1],h));
            }
            ++report.lineComponents;
            for(int j=0;j<static_cast<int>(line.size());++j){const auto found=index.find(line[j].vertex);if(found==index.end())continue;
                vertex=line[j].vertex;require(line.size()>=5,"connected master skeleton has insufficient five-point support");
                auto& d=directions[found->second][direction];require(d.target<0,"duplicate target on disconnected skeleton");
                const int begin=std::clamp(j-2,0,static_cast<int>(line.size())-5);d.target=j-begin;
                std::copy_n(line.begin()+begin,5,d.samples.begin());
                d.zero=reflection[found->second][direction];
            }
            first=last;
        }
        const int count=static_cast<int>(primary.vertices.size());
        std::vector<Eigen::Triplet<double>> entries;Eigen::VectorXd rhs(count);
        stage="compact compatibility system";
        for(int i=0;i<count;++i){vertex=primary.vertices[i].vertex;const auto& d=directions[i];
            require(d[0].target>=0 && d[1].target>=0,"master target lacks both directional supports");
            if(reflection[i][0] || reflection[i][1]){rhs[i]=0;entries.emplace_back(i,i,1);continue;}
            rhs[i]=.5*d[0].constant()+.5*d[1].constant();require(std::isfinite(rhs[i]),"nonfinite compatibility right-hand side");
            std::map<int,double> row{{i,1}};
            for(int k=0;k<2;++k)for(const auto& c:d[k].coefficients())row[c.first]-=.5*c.second;
            for(const auto& c:row){require(std::isfinite(c.second),"nonfinite compatibility matrix");if(c.second!=0)entries.emplace_back(i,c.first,c.second);}
        }
        Eigen::SparseMatrix<double> matrix(count,count);matrix.setFromTriplets(entries.begin(),entries.end());matrix.makeCompressed();
        report.compatibilityNonzeros=static_cast<std::size_t>(matrix.nonZeros());
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;solver.compute(matrix);
        require(solver.info()==Eigen::Success,"singular compact compatibility system");
        const Eigen::VectorXd solution=solver.solve(rhs);require(solver.info()==Eigen::Success && solution.allFinite(),"compact compatibility solve failed");
        const Eigen::VectorXd residual=matrix*solution-rhs;
        require(residual.allFinite(),"nonfinite compact algebraic residual");
        double algebraic=0;
        for(int i=0;i<count;++i)algebraic=std::max(algebraic,std::abs(residual[i])/std::max({1.0,std::abs(rhs[i]),std::abs(solution[i])}));
        require(algebraic<=compatibilityTolerance,"compact algebraic residual exceeds tolerance");report.algebraicResidual=algebraic;
        std::vector<double> mixed(count);
        for(int i=0;i<count;++i){const int v=primary.vertices[i].vertex;
            mixed[i]=reflection[i][0] || reflection[i][1]?0:solution[i];input.nodes[v].dxy=mixed[i];}
        stage="final exact constraints";const Skeleton skeleton(mesh,input,axis);
        for(const auto& plane:symmetryPlanes)for(const auto& v:skeleton.vertices()){
            const auto q=point(v.index);if((plane.normalAxis==0?q.x:q.y)!=plane.coordinate)continue;
            const auto jet=skeleton.nodeJet(v.index);
            require(jet[plane.normalAxis+1]==0 && jet[3]==0,"final skeleton violates declared symmetry");
        }
        double compatibility=0;
        for(int i=0;i<count;++i){vertex=primary.vertices[i].vertex;
            const double dx=directions[i][0].evaluate(skeleton,2),dy=directions[i][1].evaluate(skeleton,1);
            const double expected=.5*dx+.5*dy;require(std::isfinite(expected),"nonfinite direct compact closure");
            const double scale=std::max({1.0,std::abs(dx),std::abs(dy),std::abs(mixed[i])});
            compatibility=std::max(compatibility,std::abs(mixed[i]/scale-expected/scale));
            const double directionalScale=std::max({1.0,std::abs(dx),std::abs(dy)});
            report.maxRelativeDirectionalDisagreement=std::max(report.maxRelativeDirectionalDisagreement,
                2*std::abs(.5*dx/directionalScale-.5*dy/directionalScale));
        }
        require(compatibility<=compatibilityTolerance,"direct compact/constraint residual exceeds tolerance");report.compatibilityResidual=compatibility;
        stage="final mesh binding";vertex=-1;candidate.validateMesh(mesh);
        return HermiteGalsAmrClosureResult(identity,candidate.timeLayer(),std::move(mixed),skeleton,report);
    }catch(const std::exception& error){throw std::runtime_error("[Interface/HermiteGalsAmrClosure primary="+std::to_string(identity.value())
        +" vertex="+std::to_string(vertex)+"] "+stage+": "+error.what());}
}
} // namespace LAMP::Interface
