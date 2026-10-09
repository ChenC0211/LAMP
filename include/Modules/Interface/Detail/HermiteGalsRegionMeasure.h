#pragma once

#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Geometry/LevelSetMappedInterface.h"
#include "Core/Geometry/LevelSetMappedSegment.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Shared read-only integration of a native Q4 scalar over a leaf rectangle.
// The caller owns complete topology and supplies graphs from this same source.
// No inventory, reconstruction or accepted state is created here.
namespace Quantum::Interface::Detail::RegionMeasure {
using I=LevelSetMappedInterface;
using Snapshot=LevelSetHermiteGeometry::Snapshot;
constexpr long double pi=3.1415926535897932384626433832795L;
inline void require(bool ok,const std::string& message){if(!ok)throw std::invalid_argument("[HermiteGalsPhaseRegion] "+message);}
struct Sum {long double value=0,correction=0;void add(long double x){const auto a=x-correction,b=value+a;correction=(b-value)-a;value=b;}};
struct Box {Vector2D lower,upper;int cell=-1;};
inline long double capacity(Vector2D a,Vector2D b,bool radial){return (static_cast<long double>(b.x)-a.x)*(static_cast<long double>(b.y)-a.y)*
    (radial?pi*(static_cast<long double>(b.x)+a.x):1);}
struct CellMeasure {long double axialLiquid=0,radialLiquid=0,axialGas=0,radialGas=0,error=0;
    std::size_t arcs=0,refinements=0,coincidentArcs=0;
    long double fixedAxialAssemblyDefect=0,fixedRadialAssemblyDefect=0,assemblyCapacity=0;};
using CutKey=std::array<double,4>;
inline CutKey key(Vector2D a,Vector2D b){if(a.x>b.x || (a.x==b.x && a.y>b.y))std::swap(a,b);return {{a.x,a.y,b.x,b.y}};}
// Integrals use each retained phase interval's translated offset and length;
// never subtract two rounded physical phase endpoints.
inline std::array<long double,2> phaseWeights(const LevelSetMappedSegment& cut,bool radial){
    const auto a=cut.begin(),b=cut.end();const bool horizontal=a.y==b.y;
    const long double length=horizontal?static_cast<long double>(b.x)-a.x:static_cast<long double>(b.y)-a.y;
    require(length>0,"unordered region boundary cut");Sum negative,positive;
    for(const auto& r:cut.ranges()){
        long double value=length*r.length;
        if(horizontal && radial){const long double anchor=r.fromEnd?b.x:a.x,delta=r.fromEnd?-length:length;
            const long double mean=anchor+delta*r.offset+delta*(r.length*.5L);value*=2*pi*mean;}
        require(std::isfinite(value) && value>=0,"invalid phase boundary measure");
        (r.negative?negative:positive).add(value);
    }
    return {{negative.value,positive.value}};
}
inline CellMeasure integrateCell(const std::shared_ptr<const Snapshot>& scalar,const Box& box,
    const std::vector<I::Graph>& graphs,bool radial){
    const auto bounds=scalar->nativeCellBounds(box.cell);const auto lo=bounds[0],hi=bounds[1];
    std::vector<double> xs{lo.x,box.lower.x,box.upper.x,hi.x},ys{lo.y,box.lower.y,box.upper.y,hi.y};
    for(auto* v:{&xs,&ys}){std::sort(v->begin(),v->end());v->erase(std::unique(v->begin(),v->end()),v->end());}
    std::vector<I::Region> regions;std::map<CutKey,std::size_t> ids;
    std::vector<LevelSetMappedSegment::NativeQuarticTrace> cuts;
    for(std::size_t j=1;j<ys.size();++j)for(std::size_t k=1;k<xs.size();++k){
        const Vector2D a{xs[k-1],ys[j-1]},b{xs[k],ys[j]};
        const int owner=a.x>=box.lower.x && b.x<=box.upper.x && a.y>=box.lower.y && b.y<=box.upper.y?0:1;
        const I::Region region{owner,{{a,{b.x,a.y},b,{a.x,b.y}}}};regions.push_back(region);
        for(int e=0;e<4;++e){const auto edge=key(region.vertices[e],region.vertices[(e+1)%4]);
            if(ids.emplace(edge,cuts.size()).second)cuts.push_back(LevelSetMappedSegment::nativeQuarticTrace(
                scalar,Vector2D{edge[0],edge[1]},Vector2D{edge[2],edge[3]},edge[0]==edge[2]?0:1));}
    }
    // Green potentials admit any common origin. Keep it near the interface
    // to avoid subtracting full-capacity moments for a tiny phase at a face.
    // This is only an arithmetic gauge: no scalar, root, graph or region is
    // moved. The estimate need not be a zero of the retained scalar.
    const auto middle=(box.lower+box.upper)*.5;
    Vector2D estimate=middle;
    const auto gaugeJet=scalar->scalarJetNative(estimate);
    const double norm=gaugeJet.dx*gaugeJet.dx+gaugeJet.dy*gaugeJet.dy;
    if(norm>0 && std::isfinite(norm)) {
        estimate.x=std::clamp(estimate.x-gaugeJet.value*gaugeJet.dx/norm,box.lower.x,box.upper.x);
        estimate.y=std::clamp(estimate.y-gaugeJet.value*gaugeJet.dy/norm,box.lower.y,box.upper.y);
    }
    // Use an existing box corner: an interior gauge would introduce another
    // rounded offset in the fixed-edge full-capacity sum. Only the potential's
    // origin changes; the integration domain and all retained roots stay put.
    const Vector2D greenOrigin{estimate.x<=middle.x?box.lower.x:box.upper.x,
        estimate.y<=middle.y?box.lower.y:box.upper.y};
    std::size_t refinements=0;
    for(int attempt=0;attempt<=256;++attempt){
        try{
            std::vector<const LevelSetMappedSegment::NativeQuarticTrace*> refs;for(const auto& cut:cuts)refs.push_back(&cut);
            const auto integrated=I::integrateNativeQuarticRegionMoments(*scalar,box.cell,graphs,refs,regions,radial,greenOrigin);
            CellMeasure result;Sum fixedAxial,fixedRadial,totalCapacity;
            for(const auto& r:integrated.regions){fixedAxial.add(r.fixedAxialLiquid);fixedRadial.add(r.fixedRadialLiquid);
                totalCapacity.add(r.capacity);
                if(r.regionId==0){result.axialLiquid=r.axialLiquid;result.radialLiquid=r.radialLiquid;
                    result.axialGas=r.axialGas;result.radialGas=r.radialGas;result.error=r.errorEstimate;result.arcs=r.freeArcs;}}
            result.refinements=refinements+integrated.localRootRefinements;result.coincidentArcs=integrated.coincidentArcs;
            result.fixedAxialAssemblyDefect=fixedAxial.value-integrated.exteriorFixedAxialLiquid-integrated.coincidentAxialLiquid;
            result.fixedRadialAssemblyDefect=fixedRadial.value-integrated.exteriorFixedRadialLiquid-integrated.coincidentRadialLiquid;
            result.assemblyCapacity=totalCapacity.value;return result;
        }catch(const I::RootResolutionError& failure){
            require(attempt<256,"region root refinement work limit");
            require(!failure.requests.empty(),"empty region root refinement request");
            // One root per query is replaced in this retry. Every subsequent
            // request is refreshed after all dependent moments restart.
            std::set<std::size_t> replaced;
            for(const auto& request:failure.requests){require(request.segmentIndex<cuts.size(),"foreign refinement cut");
                if(replaced.insert(request.segmentIndex).second){auto& cut=cuts[request.segmentIndex];
                    require(request.rootIndex<cut.events().size(),"foreign refinement event");cut=cut.refined(request.rootIndex);++refinements;}}
        }
    }
    throw std::invalid_argument("[HermiteGalsPhaseRegion] unreachable refinement exhaustion");
}
// Select the directly integrated smaller phase, then form its complement.
// Inputs must first pass the caller's strict capacity/ledger checks.
inline std::array<long double,2> complement(long double cap,long double liquid,long double gas){
    if(liquid<=gas)return {{liquid,cap-liquid}};
    return {{cap-gas,gas}};
}
}
