#include "Core/Mesh/NativeLeafCoverage.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Geometry/Detail/LevelSetMappedCoordinates.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

NativeLeafCoverage::NativeLeafCoverage(const QuadtreeMesh& mesh,const NativeLeafDomain& domain)
    : mesh_(mesh),domain_(domain),locator_(mesh),
      complete_(domain.cellIndexByRow().size()==static_cast<std::size_t>(mesh.cellIndexNum()))
{
    domain.validateMeshIdentity(mesh);
}

bool NativeLeafCoverage::covers(Vector2D a,Vector2D b) const
{
    namespace R=Quantum::Geometry::Detail::MappedHermiteSegment;
    R::Exact::validateEnvironment();
    const auto lower=mesh_.minCoord,upper=mesh_.maxCoord;
    const auto inside=[&](Vector2D p) {
        return std::isfinite(p.x) && std::isfinite(p.y) && p.x>=lower.x && p.x<=upper.x
            && p.y>=lower.y && p.y<=upper.y;
    };
    if (!inside(a) || !inside(b)) return false;
    if (complete_) return true;
    // Native closed-box search includes both incident coarse/fine owners.
    // Classification and merging below retain exact rational intervals;
    // no midpoint sampling, virtual cells or gap tolerance is introduced.
    const auto candidates=locator_.intersectingBox({std::min(a.x,b.x),std::min(a.y,b.y)},
        {std::max(a.x,b.x),std::max(a.y,b.y)});
    const auto zero=R::constant(0),one=R::constant(1);
    std::vector<R::Chart> intervals;
    for (int c:candidates) {
        const auto& cell=mesh_.cells(mesh_.cellIndexToLabel[c]);
        if (domain_.cellRowByIndex()[cell.index()]<0) continue;
        const auto lo=mesh_.verts(cell.cellToVert[0]).coord(),hi=mesh_.verts(cell.cellToVert[2]).coord();
        auto begin=zero,end=one;
        const auto clip=[&](double x,double y,double l,double u) {
            if (x==y) return x>=l && x<=u;
            auto p=R::ratio(R::difference(l,x),R::difference(y,x));
            auto q=R::ratio(R::difference(u,x),R::difference(y,x));
            if (R::compare(p,q)>0) std::swap(p,q);
            if (R::compare(p,begin)>0) begin=std::move(p);
            if (R::compare(q,end)<0) end=std::move(q);
            return R::compare(begin,end)<=0;
        };
        if (clip(a.x,b.x,lo.x,hi.x) && clip(a.y,b.y,lo.y,hi.y))
            intervals.push_back({std::move(begin),std::move(end)});
    }
    std::sort(intervals.begin(),intervals.end(),[](const auto& x,const auto& y){return R::compare(x.begin,y.begin)<0;});
    auto end=zero;
    for (const auto& interval:intervals) {
        if (R::compare(interval.begin,end)>0) return false;
        if (R::compare(interval.end,end)>0) end=interval.end;
        if (R::compare(end,one)==0) return true;
    }
    return false;
}
