#pragma once

#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Geometry/LevelSetMappedInterface.h"
#include "Core/Geometry/LevelSetMappedSegment.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <chrono>
#include <cstdint>
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
namespace LAMP::Interface::Detail::RegionMeasure {
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
// Local passive work profile for ONE actual pair call. No query-global writes
// occur inside geometry preparation, Core integration or a root retry.
struct HalfPairWorkProfile {
    std::uint64_t regionGaugeNanoseconds=0,cutPreparationNanoseconds=0,coreNanoseconds=0,
        retryNanoseconds=0,assemblyNanoseconds=0;
    std::size_t cutConstructions=0,coreAttempts=0,faceCutRefinements=0;
    I::NativeMomentWorkProfile coreWork;
};
// Diagnostic-only ordered accumulation. Each worker owns its local profile;
// callers merge once after the synchronous actual cell/group work completes.
inline void accumulateCoreWorkProfile(I::NativeMomentWorkProfile& target,
    const I::NativeMomentWorkProfile& source)noexcept{
    target.restrictionConstruction+=source.restrictionConstruction;
    target.initialRootCacheLookup+=source.initialRootCacheLookup;
    target.initialRootCacheRemember+=source.initialRootCacheRemember;
    target.initialRootIsolation+=source.initialRootIsolation;
    target.rootMappingRefinement+=source.rootMappingRefinement;
    target.fixedAndSetup+=source.fixedAndSetup;
    target.remainingIntegration+=source.remainingIntegration;
    target.total+=source.total;
    target.restrictionQueries+=source.restrictionQueries;
    target.initialRootCacheQueries+=source.initialRootCacheQueries;
    target.initialRootCacheHits+=source.initialRootCacheHits;
    target.initialRootCacheMisses+=source.initialRootCacheMisses;
    target.initialRootCacheRememberQueries+=source.initialRootCacheRememberQueries;
    target.initialRootIsolations+=source.initialRootIsolations;
    target.rootMappings+=source.rootMappings;
    target.rootRefinements+=source.rootRefinements;
}
class HalfPairWorkTimer {
public:
    explicit HalfPairWorkTimer(std::uint64_t* elapsed):elapsed_(elapsed),begin_(elapsed?Clock::now():Clock::time_point{}){}
    ~HalfPairWorkTimer()noexcept{if(elapsed_)record(Clock::now());}
    void next(std::uint64_t* elapsed)noexcept{
        if(elapsed_){const auto now=Clock::now();record(now);begin_=now;}
        else if(elapsed)begin_=Clock::now();
        elapsed_=elapsed;
    }
private:
    using Clock=std::chrono::steady_clock;
    void record(Clock::time_point now)noexcept{
        *elapsed_+=static_cast<std::uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(now-begin_).count());
    }
    std::uint64_t* elapsed_;
    Clock::time_point begin_;
};
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
    const std::vector<I::Graph>& graphs,bool radial,std::vector<I::GraphMoments>* graphMoments=nullptr,
    I::NativeMomentWorkProfile* coreWork=nullptr){
    const auto bounds=scalar->nativeCellBounds(box.cell);const auto lo=bounds[0],hi=bounds[1];
    require(!graphMoments || (box.lower.x==lo.x && box.lower.y==lo.y && box.upper.x==hi.x && box.upper.y==hi.y),
        "global graph moments require the complete actual native cell");
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
            std::vector<I::GraphMoments> preparedGraphMoments;
            const auto integrated=graphMoments
                ? I::integrateNativeQuarticRegionMomentsBatch(*scalar,box.cell,graphs,refs,
                    {{regions,greenOrigin}},radial,&preparedGraphMoments,coreWork).at(0)
                : I::integrateNativeQuarticRegionMoments(*scalar,box.cell,graphs,refs,regions,radial,greenOrigin);
            CellMeasure result;Sum fixedAxial,fixedRadial,totalCapacity;
            for(const auto& r:integrated.regions){fixedAxial.add(r.fixedAxialLiquid);fixedRadial.add(r.fixedRadialLiquid);
                totalCapacity.add(r.capacity);
                if(r.regionId==0){result.axialLiquid=r.axialLiquid;result.radialLiquid=r.radialLiquid;
                    result.axialGas=r.axialGas;result.radialGas=r.radialGas;result.error=r.errorEstimate;result.arcs=r.freeArcs;}}
            result.refinements=refinements+integrated.localRootRefinements;result.coincidentArcs=integrated.coincidentArcs;
            result.fixedAxialAssemblyDefect=fixedAxial.value-integrated.exteriorFixedAxialLiquid-integrated.coincidentAxialLiquid;
            result.fixedRadialAssemblyDefect=fixedRadial.value-integrated.exteriorFixedRadialLiquid-integrated.coincidentRadialLiquid;
            result.assemblyCapacity=totalCapacity.value;
            if(graphMoments)*graphMoments=std::move(preparedGraphMoments);
            return result;
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
// One common partition and graph traversal, with each original box's own
// directly evaluated Green potential. The atom labels locate arcs only;
// no quarter phase volume is subtracted to recover a small half-box phase.
inline std::array<Box,4> halfBoxes(const Snapshot& scalar,int cell,Vector2D middle){
    const auto b=scalar.nativeCellBounds(cell);const auto lo=b[0],hi=b[1];
    require(lo.x<middle.x && middle.x<hi.x && lo.y<middle.y && middle.y<hi.y,"invalid actual leaf center");
    return {{{lo,{middle.x,hi.y},cell},{{middle.x,lo.y},hi,cell},
        {lo,{hi.x,middle.y},cell},{{lo.x,middle.y},hi,cell}}};
}
inline std::vector<CellMeasure> integrateBoxesBatch(const std::shared_ptr<const Snapshot>& scalar,
    const std::vector<Box>& boxes,const std::vector<I::Graph>& graphs,bool radial){
    require(!boxes.empty(),"empty region target batch");
    const int cell=boxes[0].cell;const auto b=scalar->nativeCellBounds(cell);
    std::vector<double> xs{b[0].x,b[1].x},ys{b[0].y,b[1].y};
    for(const auto& box:boxes){require(box.cell==cell,"batch targets use different native cells");
        xs.push_back(box.lower.x);xs.push_back(box.upper.x);ys.push_back(box.lower.y);ys.push_back(box.upper.y);}
    for(auto* values:{&xs,&ys}){std::sort(values->begin(),values->end());values->erase(std::unique(values->begin(),values->end()),values->end());}
    std::vector<I::NativeMomentTarget> targets;
    for(const auto& box:boxes){
        I::NativeMomentTarget target;
        for(std::size_t j=1;j<ys.size();++j)for(std::size_t k=1;k<xs.size();++k){
            const Vector2D a{xs[k-1],ys[j-1]},z{xs[k],ys[j]};
            const int owner=a.x>=box.lower.x && z.x<=box.upper.x && a.y>=box.lower.y && z.y<=box.upper.y?0:1;
            target.regions.push_back({owner,{{a,{z.x,a.y},z,{a.x,z.y}}}});
        }
        // Identical arithmetic and original half-box corner gauge as single.
        const auto middle=(box.lower+box.upper)*.5;Vector2D estimate=middle;
        const auto jet=scalar->scalarJetNative(estimate);const double norm=jet.dx*jet.dx+jet.dy*jet.dy;
        if(norm>0 && std::isfinite(norm)){
            estimate.x=std::clamp(estimate.x-jet.value*jet.dx/norm,box.lower.x,box.upper.x);
            estimate.y=std::clamp(estimate.y-jet.value*jet.dy/norm,box.lower.y,box.upper.y);
        }
        target.greenOrigin={estimate.x<=middle.x?box.lower.x:box.upper.x,
            estimate.y<=middle.y?box.lower.y:box.upper.y};
        targets.push_back(std::move(target));
    }
    std::map<CutKey,std::size_t> ids;std::vector<LevelSetMappedSegment::NativeQuarticTrace> cuts;
    for(const auto& region:targets[0].regions)for(int e=0;e<4;++e){
        const auto edge=key(region.vertices[e],region.vertices[(e+1)%4]);
        if(ids.emplace(edge,cuts.size()).second)cuts.push_back(LevelSetMappedSegment::nativeQuarticTrace(
            scalar,{edge[0],edge[1]},{edge[2],edge[3]},edge[0]==edge[2]?0:1));
    }
    std::size_t refinements=0;
    for(int attempt=0;attempt<=256;++attempt){
        try{
            std::vector<const LevelSetMappedSegment::NativeQuarticTrace*> refs;for(const auto& cut:cuts)refs.push_back(&cut);
            const auto integrated=I::integrateNativeQuarticRegionMomentsBatch(*scalar,cell,graphs,refs,targets,radial);
            require(integrated.size()==boxes.size(),"batch lost an actual target");
            std::vector<CellMeasure> result(boxes.size());
            for(std::size_t t=0;t<result.size();++t){const auto& moments=integrated[t];auto& r=result[t];
                Sum fixedAxial,fixedRadial,totalCapacity;
                for(const auto& region:moments.regions){fixedAxial.add(region.fixedAxialLiquid);
                    fixedRadial.add(region.fixedRadialLiquid);totalCapacity.add(region.capacity);
                    if(region.regionId==0){r.axialLiquid=region.axialLiquid;r.radialLiquid=region.radialLiquid;
                        r.axialGas=region.axialGas;r.radialGas=region.radialGas;r.error=region.errorEstimate;r.arcs=region.freeArcs;}}
                r.refinements=refinements+moments.localRootRefinements;r.coincidentArcs=moments.coincidentArcs;
                r.fixedAxialAssemblyDefect=fixedAxial.value-moments.exteriorFixedAxialLiquid-moments.coincidentAxialLiquid;
                r.fixedRadialAssemblyDefect=fixedRadial.value-moments.exteriorFixedRadialLiquid-moments.coincidentRadialLiquid;
                r.assemblyCapacity=totalCapacity.value;
            }
            return result;
        }catch(const I::RootResolutionError& failure){
            require(attempt<256,"batch region root refinement work limit");require(!failure.requests.empty(),"empty batch root request");
            std::set<std::size_t> replaced;
            for(const auto& request:failure.requests){require(request.segmentIndex<cuts.size(),"foreign batch cut");
                if(replaced.insert(request.segmentIndex).second){auto& cut=cuts[request.segmentIndex];
                    require(request.rootIndex<cut.events().size(),"foreign batch event");cut=cut.refined(request.rootIndex);++refinements;}}
        }
    }
    throw std::invalid_argument("[HermiteGalsPhaseRegion] unreachable batch refinement exhaustion");
}
// Opposite half targets already have an identical original partition.
// Pairing by direction avoids introducing the other direction's extra cuts.
inline std::array<CellMeasure,2> integrateHalfPair(const std::shared_ptr<const Snapshot>& scalar,
    const std::array<Box,2>& boxes,const std::vector<I::Graph>& graphs,bool radial,
    HalfPairWorkProfile* work=nullptr){
    HalfPairWorkTimer timer(work?&work->regionGaugeNanoseconds:nullptr);
    require(!boxes.empty(),"empty region target batch");
    const int cell=boxes[0].cell;const auto b=scalar->nativeCellBounds(cell);
    std::vector<double> xs{b[0].x,b[1].x},ys{b[0].y,b[1].y};
    for(const auto& box:boxes){require(box.cell==cell,"batch targets use different native cells");
        xs.push_back(box.lower.x);xs.push_back(box.upper.x);ys.push_back(box.lower.y);ys.push_back(box.upper.y);}
    for(auto* values:{&xs,&ys}){std::sort(values->begin(),values->end());values->erase(std::unique(values->begin(),values->end()),values->end());}
    std::vector<I::NativeSelectedMomentTarget> targets;
    for(const auto& box:boxes){
        I::NativeSelectedMomentTarget target;
        target.selectedRegionIds={0};
        for(std::size_t j=1;j<ys.size();++j)for(std::size_t k=1;k<xs.size();++k){
            const Vector2D a{xs[k-1],ys[j-1]},z{xs[k],ys[j]};
            const int owner=a.x>=box.lower.x && z.x<=box.upper.x && a.y>=box.lower.y && z.y<=box.upper.y?0:1;
            target.regions.push_back({owner,{{a,{z.x,a.y},z,{a.x,z.y}}}});
        }
        // Identical arithmetic and original half-box corner gauge as single.
        const auto middle=(box.lower+box.upper)*.5;Vector2D estimate=middle;
        const auto jet=scalar->scalarJetNative(estimate);const double norm=jet.dx*jet.dx+jet.dy*jet.dy;
        if(norm>0 && std::isfinite(norm)){
            estimate.x=std::clamp(estimate.x-jet.value*jet.dx/norm,box.lower.x,box.upper.x);
            estimate.y=std::clamp(estimate.y-jet.value*jet.dy/norm,box.lower.y,box.upper.y);
        }
        target.greenOrigin={estimate.x<=middle.x?box.lower.x:box.upper.x,
            estimate.y<=middle.y?box.lower.y:box.upper.y};
        targets.push_back(std::move(target));
    }
    timer.next(work?&work->cutPreparationNanoseconds:nullptr);
    std::map<CutKey,std::size_t> ids;std::vector<LevelSetMappedSegment::NativeQuarticTrace> cuts;
    for(const auto& region:targets[0].regions)for(int e=0;e<4;++e){
        const auto edge=key(region.vertices[e],region.vertices[(e+1)%4]);
        if(ids.emplace(edge,cuts.size()).second){
            if(work)++work->cutConstructions;
            cuts.push_back(LevelSetMappedSegment::nativeQuarticTrace(
                scalar,{edge[0],edge[1]},{edge[2],edge[3]},edge[0]==edge[2]?0:1));
        }
    }
    timer.next(work?&work->coreNanoseconds:nullptr);
    std::size_t refinements=0;
    for(int attempt=0;attempt<=256;++attempt){
        try{
            std::vector<const LevelSetMappedSegment::NativeQuarticTrace*> refs;for(const auto& cut:cuts)refs.push_back(&cut);
            if(work)++work->coreAttempts;
            const auto integrated=I::integrateNativeQuarticSelectedRegionMomentsBatch(*scalar,cell,graphs,refs,targets,radial,
                work?&work->coreWork:nullptr);
            timer.next(work?&work->assemblyNanoseconds:nullptr);
            require(integrated.size()==boxes.size(),"batch lost an actual target");
            std::array<CellMeasure,2> result{};
            for(std::size_t t=0;t<result.size();++t){const auto& moments=integrated[t];auto& r=result[t];
                Sum fixedAxial,fixedRadial,totalCapacity;
                // Complete fixed/coincident assembly uses BOTH original halves.
                // Only the requested half needs free phase-volume channels.
                for(const auto& region:moments.fixedRegions){fixedAxial.add(region.fixedAxialLiquid);
                    fixedRadial.add(region.fixedRadialLiquid);totalCapacity.add(region.capacity);}
                require(moments.selectedRegions.size()==1 && moments.selectedRegions[0].regionId==0,
                    "selected batch lost the actual half target");
                const auto& region=moments.selectedRegions[0];
                r.axialLiquid=region.axialLiquid;r.radialLiquid=region.radialLiquid;
                r.axialGas=region.axialGas;r.radialGas=region.radialGas;r.error=region.errorEstimate;r.arcs=region.freeArcs;
                r.refinements=refinements+moments.localRootRefinements;r.coincidentArcs=moments.coincidentArcs;
                r.fixedAxialAssemblyDefect=fixedAxial.value-moments.exteriorFixedAxialLiquid-moments.coincidentAxialLiquid;
                r.fixedRadialAssemblyDefect=fixedRadial.value-moments.exteriorFixedRadialLiquid-moments.coincidentRadialLiquid;
                r.assemblyCapacity=totalCapacity.value;
            }
            return result;
        }catch(const I::RootResolutionError& failure){
            timer.next(work?&work->retryNanoseconds:nullptr);
            require(attempt<256,"batch region root refinement work limit");require(!failure.requests.empty(),"empty batch root request");
            std::set<std::size_t> replaced;
            for(const auto& request:failure.requests){require(request.segmentIndex<cuts.size(),"foreign batch cut");
                if(replaced.insert(request.segmentIndex).second){auto& cut=cuts[request.segmentIndex];
                    require(request.rootIndex<cut.events().size(),"foreign batch event");
                    if(work)++work->faceCutRefinements;
                    cut=cut.refined(request.rootIndex);++refinements;}}
            timer.next(work?&work->coreNanoseconds:nullptr);
        }
    }
    throw std::invalid_argument("[HermiteGalsPhaseRegion] unreachable batch refinement exhaustion");
}
inline std::array<CellMeasure,4> integrateHalfCells(const std::shared_ptr<const Snapshot>& scalar,
    const std::array<Box,4>& boxes,const std::vector<I::Graph>& graphs,bool radial){
    const auto result=integrateBoxesBatch(scalar,{boxes.begin(),boxes.end()},graphs,radial);
    return {{result[0],result[1],result[2],result[3]}};
}
// Select the directly integrated smaller phase, then form its complement.
// Inputs must first pass the caller's strict capacity/ledger checks.
inline std::array<long double,2> complement(long double cap,long double liquid,long double gas){
    if(liquid<=gas)return {{liquid,cap-liquid}};
    return {{cap-gas,gas}};
}
}
