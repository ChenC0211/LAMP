#include "Core/Geometry/ParametricBoundary.h"
#include <set>

namespace Quantum::Geometry::Parametric {
namespace {
struct Box {std::array<double,2> x,y;};
Box box(const Arc& a){return {range(a.coordinate(0)),range(a.coordinate(1))};}
bool separated(const Box& a,const Box& b){
    const double reserve=64*std::numeric_limits<double>::epsilon()*std::max({1.,std::abs(a.x[0]),std::abs(a.x[1]),std::abs(a.y[0]),std::abs(a.y[1]),std::abs(b.x[0]),std::abs(b.x[1]),std::abs(b.y[0]),std::abs(b.y[1])});
    return a.x[1]+reserve<b.x[0]||b.x[1]+reserve<a.x[0]||a.y[1]+reserve<b.y[0]||b.y[1]+reserve<a.y[0];
}
double diameter(const Box& a){return std::max(a.x[1]-a.x[0],a.y[1]-a.y[0]);}
void disjoint(const Arc& a,const Arc& b,int depth=0){
    const auto x=box(a),y=box(b);if(separated(x,y))return;
    check(depth<48,"material components touch/intersect or their separation is unresolved");
    if(diameter(x)>=diameter(y)){
        const double m=.5*(a.a+a.b);check(m>a.a&&m<a.b,"component separation exhausted material resolution");
        disjoint(restrictedArc(a,a.a,m),b,depth+1);disjoint(restrictedArc(a,m,a.b),b,depth+1);
    }else{
        const double m=.5*(b.a+b.b);check(m>b.a&&m<b.b,"component separation exhausted material resolution");
        disjoint(a,restrictedArc(b,b.a,m),depth+1);disjoint(a,restrictedArc(b,m,b.b),depth+1);
    }
}
}
Boundary::Boundary(std::vector<Component> components):components_(std::move(components)){
    std::set<std::uint64_t> ids;
    for(const auto& c:components_)check(c.id&&ids.insert(c.id).second&&c.curve&&!c.curve->arcs.empty(),"invalid/duplicate material component identity");
    for(std::size_t i=0;i<components_.size();i++)for(std::size_t j=0;j<i;j++){
        const auto& a=*components_[i].curve;const auto& b=*components_[j].curve;
        auto aa=a.arcs,bb=b.arcs;const auto ac=a.closure(),bc=b.closure();aa.insert(aa.end(),ac.begin(),ac.end());bb.insert(bb.end(),bc.begin(),bc.end());
        for(const auto& x:aa)for(const auto& y:bb)disjoint(x,y);
        // After boundary separation, one representative suffices to reject
        // strict containment. Nested liquid regions require a hole/phase
        // ownership model, not two positive overlapping liquid inventories.
        check(signedDistance(a,b.arcs.front().p[0]).value>0&&signedDistance(b,a.arcs.front().p[0]).value>0,
            "material components overlap by containment");
    }
}
}
