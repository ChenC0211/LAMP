#include "Modules/Interface/HermiteGalsEdgeTransport.h"

#include "Modules/Interface/detail/HermiteGalsFirstJet.h"
#include <algorithm>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>

namespace LAMP::Interface
{
namespace
{
using Key=std::array<int,2>;
using Detail::requireFirstJet;
struct Side { Key key; bool horizontal; double line,lo,hi; };

// Maximal collinear leaf sides are the independent edge midpoint DOFs.
// Classify by native topology, not by endpoint vertex type: a perpendicular
// master edge can legitimately end at a hanging vertex.
std::set<Key> nativeMasterSides(const QuadtreeMesh& mesh)
{
    std::map<Key,Side> unique;
    for(int c=0;c<mesh.cellIndexNum();++c){
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);
        for(int k=0;k<4;++k){
            const auto& va=mesh.verts(cell.cellToVert[k]);
            const auto& vb=mesh.verts(cell.cellToVert[(k+1)%4]);
            auto a=va.coord(),b=vb.coord();Key key{{va.index(),vb.index()}};
            const bool horizontal=a.y==b.y;
            requireFirstJet((horizontal && a.x!=b.x) || (a.x==b.x && a.y!=b.y),"invalid native Cartesian side");
            if(horizontal?a.x>b.x:a.y>b.y){std::swap(a,b);std::swap(key[0],key[1]);}
            unique.emplace(key,Side{key,horizontal,horizontal?a.y:a.x,horizontal?a.x:a.y,horizontal?b.x:b.y});
        }
    }
    std::vector<Side> sides;sides.reserve(unique.size());
    for(const auto& entry:unique)sides.push_back(entry.second);
    std::sort(sides.begin(),sides.end(),[](const Side& a,const Side& b){
        if(a.horizontal!=b.horizontal)return a.horizontal<b.horizontal;
        if(a.line!=b.line)return a.line<b.line;
        if(a.lo!=b.lo)return a.lo<b.lo;
        return a.hi>b.hi;
    });
    std::set<Key> masters;
    bool first=true,horizontal=false;double line=0,upper=0;
    for(const auto& side:sides){
        if(first || side.horizontal!=horizontal || side.line!=line){
            first=false;horizontal=side.horizontal;line=side.line;upper=side.lo;
        }
        if(side.hi>upper)masters.insert(side.key);
        upper=std::max(upper,side.hi);
    }
    return masters;
}
struct Target { Key key; Vector2D point; bool horizontal; };
}

HermiteGalsEdgeTransportResult::HermiteGalsEdgeTransportResult(
    HermiteGalsIdentity source,HermiteGalsTimeLayer sourceTime,
    const HermiteGalsVelocityView& arrival,std::vector<HermiteGalsEdgeJet> edges)
    : source_(source),arrival_(arrival.geometryIdentity()),sourceTime_(sourceTime),
      arrivalTime_(arrival.timeLayer()),edges_(std::move(edges)) {}

HermiteGalsEdgeTransportResult transportHermiteGalsEdges(const QuadtreeMesh& mesh,const HermiteGalsState& source,
    const std::vector<std::array<int,2>>& masterEdges,const LevelSetBoundaryTrace::Domain& domain,
    const HermiteGalsVelocityView& previousVelocity,const HermiteGalsVelocityView& arrivalVelocity)
{
    HermiteGalsIdentity sourceId;Key edge{{-1,-1}};const char* stage="input";
    try{
        sourceId=source.identity();source.validateMesh(mesh);
        requireFirstJet(!masterEdges.empty(),"a nonempty master edge selection is required");
        const auto geometry=source.geometry();const auto& scalar=geometry.snapshot();
        Detail::validateFirstJetSource(mesh,source,scalar,domain,previousVelocity,arrivalVelocity);
        requireFirstJet(!scalar.hasCornerSectors(),"corner-sector geometry requires explicit one-sided edge transport");
        const auto masters=nativeMasterSides(mesh);
        std::set<Key> sourceEdges;
        for(const auto& jet:source.primary().edges)sourceEdges.insert(jet.vertices);
        std::vector<Target> targets;targets.reserve(masterEdges.size());
        Key previous{{-1,-1}};
        for(const auto& key:masterEdges){
            edge=key;
            requireFirstJet(key[0]>=0 && key[1]>=0 && key[0]<mesh.vertIndexNum() && key[1]<mesh.vertIndexNum(),
                "edge endpoints must be native vertex indices");
            requireFirstJet(previous<key,"target edges must be lexicographically ascending and unique");
            const auto a=mesh.verts(mesh.vertIndexToLabel[key[0]]).coord();
            const auto b=mesh.verts(mesh.vertIndexToLabel[key[1]]).coord();
            const bool horizontal=a.y==b.y;
            requireFirstJet((horizontal && a.x<b.x) || (a.x==b.x && a.y<b.y),
                "edge endpoints must follow the positive Cartesian tangent");
            requireFirstJet(masters.count(key)!=0,"target is not an independent native master side (non-side or constrained subedge)");
            requireFirstJet(sourceEdges.count(key)!=0,"target is not a source master edge");
            const Vector2D point=horizontal?Vector2D{a.x+(b.x-a.x)*.5,a.y}:Vector2D{a.x,a.y+(b.y-a.y)*.5};
            requireFirstJet(horizontal?(point.x>a.x && point.x<b.x):(point.y>a.y && point.y<b.y),
                "edge midpoint is unresolved in floating point");
            targets.push_back({key,point,horizontal});previous=key;
        }
        std::vector<HermiteGalsEdgeJet> values;values.reserve(targets.size());
        for(const auto& target:targets){
            edge=target.key;
            const auto jet=Detail::pullbackFirstJet(scalar,domain,target.point,previousVelocity,arrivalVelocity,stage);
            values.push_back({target.key,jet.phi,target.horizontal?jet.gradient.y:jet.gradient.x});
        }
        stage="final mesh binding";edge={{-1,-1}};source.validateMesh(mesh);
        return HermiteGalsEdgeTransportResult(sourceId,source.timeLayer(),arrivalVelocity,std::move(values));
    }catch(const std::exception& error){
        throw std::runtime_error("[Interface/HermiteGalsEdgeTransport source="+std::to_string(sourceId.value())
            +" arrival="+std::to_string(arrivalVelocity.geometryIdentity().value())+" edge=("+std::to_string(edge[0])
            +","+std::to_string(edge[1])+")] "+stage+": "+error.what());
    }
}
} // namespace LAMP::Interface
