#include "Hydrodynamics/HydrodynamicsVofRegrid.h"
#include "Modules/Hydrodynamics/HydroVofRefinement.h"
#include "Modules/Hydrodynamics/HydroTwoPhaseStep.h"
#include "Core/Field/CellField.h"
#include "Core/Field/NativeTargetLevelField.h"
#include "Core/Mesh/DyadicGridCoordinates.h"
#include "Core/Field/HEdgeField.h"
#include "Core/Field/VEdgeField.h"
#include "Core/Field/VertField.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/NativeLeafOverlap.h"
#include <algorithm>
#include <chrono>
#include <cfenv>
#include <cerrno>
#include <exception>
#if defined(_MSC_VER) && defined(_M_X64)
#include <float.h>
#include <xmmintrin.h>
#endif
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <tuple>
#include <vector>
namespace
{
void require(bool ok,const char* why)
{if(!ok)throw std::invalid_argument(std::string("[Hydrodynamics VOF regrid] ")+why);}
using WallClock=std::chrono::steady_clock;
double seconds(WallClock::time_point start)
{return std::chrono::duration<double>(WallClock::now()-start).count();}
// New observations restore only their own clock/conversion arithmetic.
// No scope below ever encloses a mesh/geometry/adaptation business call.
struct TreeTimingEnvironment {
    std::fenv_t native{};int error=errno;bool captured=false,eligible=false;
#if defined(_MSC_VER) && defined(_M_X64)
    unsigned csr=_mm_getcsr();
#endif
    TreeTimingEnvironment() noexcept {
        captured=std::fegetenv(&native)==0;
#if defined(_MSC_VER) && defined(_M_X64)
        unsigned control=0;
        eligible=captured&&_controlfp_s(&control,0,0)==0
            &&(control&_MCW_EM)==_MCW_EM&&(csr&0x1f80u)==0x1f80u
            &&std::fegetround()==FE_TONEAREST;
#endif
    }
    ~TreeTimingEnvironment() noexcept {
        if(captured&&std::fesetenv(&native)!=0)std::terminate();
#if defined(_MSC_VER) && defined(_M_X64)
        _mm_setcsr(csr);
#endif
        errno=error;
    }
};
struct TreeTimingStamp {WallClock::time_point value{};bool valid=false;};
TreeTimingStamp treeTimingStart() noexcept {
    const TreeTimingEnvironment environment;
    if(!environment.eligible)return {};
    return {WallClock::now(),true};
}
void recordTreeTiming(TreeTimingStamp start,double& sum,std::size_t& calls,bool& available) noexcept {
    const TreeTimingEnvironment environment;++calls;
    if(!start.valid||!environment.eligible){available=false;return;}
    const auto elapsed=std::chrono::duration<double>(WallClock::now()-start.value).count();
    if(!std::isfinite(elapsed)||elapsed<0){available=false;return;}
    sum+=elapsed;if(!std::isfinite(sum))available=false;
}
void finishTreeTiming(HydrodynamicsVofRegridReport& report) noexcept {
    const TreeTimingEnvironment environment;auto& timing=report.treePreparation;
    if(!timing.available||!environment.eligible){timing.available=false;return;}
    timing.unattributedSeconds=report.adaptationSeconds-timing.cloneSeconds-timing.graphSeconds
        -timing.domainSeconds-timing.overlapSeconds-timing.nativeAdaptSeconds;
    if(!std::isfinite(timing.unattributedSeconds))timing.available=false;
}
using LeafSignature=std::vector<std::tuple<unsigned,std::uint64_t,std::uint64_t>>;
LeafSignature leafSignature(const QuadtreeMesh& mesh)
{
    LeafSignature result;result.reserve(mesh.cellIndexNum());
    for(int i=0;i<mesh.cellIndexNum();++i){const auto q=LAMP::Mesh::cellAddress(mesh,mesh.cells(mesh.cellIndexToLabel[i]));
        result.emplace_back(q.level,q.i,q.j);}
    std::sort(result.begin(),result.end());
    require(std::adjacent_find(result.begin(),result.end())==result.end(),"duplicate actual dyadic leaf signature");
    return result;
}
std::string targetCounts(const QuadtreeMesh& mesh,const CellField<double>& targets)
{
    std::size_t refine=0,coarsen=0,pending=0;
    for(int i=0;i<mesh.cellIndexNum();++i){const auto& c=mesh.cells(mesh.cellIndexToLabel[i]);
        refine+=c.level()<targets.value[i];coarsen+=c.level()>targets.value[i];pending+=c.refineFlag!=0||c.coarsenFlag!=0;}
    return "; unmet_refine="+std::to_string(refine)+" unmet_coarsen="+std::to_string(coarsen)+" pending_votes="+std::to_string(pending);
}
void validateRequirementsSource(const QuadtreeMesh& mesh,
    const LAMP::Hydrodynamics::HydroTwoPhaseState& accepted)
{
    using LAMP::Hydrodynamics::HydroPhaseGeometry;
    require(accepted.model()==HydroPhaseGeometry::Model::GeometricVof&&!accepted.prescribedFlow(),
        "cell requirements need a genuine VOF source");
    const auto clock=accepted.timeLayer();(void)clock.time();
    require(clock.position==LAMP::TimeLayer::Position::End,"cell requirements need the complete accepted End clock");
    const auto& publication=accepted.phaseGeometry().vofPublication();
    accepted.validateVofFor(mesh,publication.identity(),clock);
    const auto& domain=publication.domain();const auto n=static_cast<std::size_t>(mesh.cellIndexNum());
    const auto& layout=accepted.values().layout();
    require(domain.cellRowByIndex().size()==n
        &&domain.cellIndexByRow().size()==static_cast<std::size_t>(layout.pressureCount())
        &&layout.domain().cellRowByIndex()==domain.cellRowByIndex()
        &&layout.domain().cellIndexByRow()==domain.cellIndexByRow(),
        "cell requirements need the exact selected native maps");
    for(std::size_t ci=0;ci<n;++ci){const int row=domain.cellRowByIndex()[ci];
        require(layout.pressureDof(static_cast<int>(ci))==row,"cell requirements have a different native pressure map");
        if(row<0)continue;
        require(row<layout.pressureCount()&&domain.cellIndexByRow()[row]==static_cast<int>(ci)
            &&layout.nativeCell(row)==static_cast<int>(ci),"cell requirements selected native maps are not invertible");
    }
}
}
bool hydrodynamicsVofRegridDue(const HydrodynamicsVofRegridPolicy& policy,const LAMP::TimeLayer& target)
{
    require(policy.interval>0,"regrid interval must be a positive integer");
    (void)target.time();
    require(target.position==LAMP::TimeLayer::Position::End,"regrid interval requires the complete End clock");
    return target.step%static_cast<std::int64_t>(policy.interval)==0;
}
struct HydrodynamicsVofRegridRequirements::Data
{
    LAMP::Hydrodynamics::HydroTwoPhaseState source;
    std::vector<int> levels;
    Data(const LAMP::Hydrodynamics::HydroTwoPhaseState& accepted,std::vector<int> minimumLevels)
        :source(accepted),levels(std::move(minimumLevels)){}
};
HydrodynamicsVofRegridRequirements::HydrodynamicsVofRegridRequirements(std::shared_ptr<const Data> data)
    :data_(std::move(data)){}
HydrodynamicsVofRegridRequirements HydrodynamicsVofRegridRequirements::fromCellMinimumLevels(
    const QuadtreeMesh& mesh,const LAMP::Hydrodynamics::HydroTwoPhaseState& accepted,std::vector<int> levels)
{
    validateRequirementsSource(mesh,accepted);
    require(levels.size()==static_cast<std::size_t>(mesh.cellIndexNum()),"cell requirements need one level for every native leaf");
    for(const int level:levels)require(level>=0&&level<=mesh.maxLevel(),"cell requirement outside allocated native levels");
    return HydrodynamicsVofRegridRequirements(std::make_shared<Data>(accepted,std::move(levels)));
}
const std::vector<int>& HydrodynamicsVofRegridRequirements::minimumLevels() const
{require(data_!=nullptr,"empty cell requirements");return data_->levels;}
LAMP::TimeLayer HydrodynamicsVofRegridRequirements::sourceClock() const
{require(data_!=nullptr,"empty cell requirements");return data_->source.timeLayer();}
void HydrodynamicsVofRegridRequirements::validateFor(const QuadtreeMesh& mesh,
    const LAMP::Hydrodynamics::HydroTwoPhaseState& accepted,const HydrodynamicsVofRegridPolicy& policy) const
{
    require(data_!=nullptr,"empty cell requirements");
    const auto actual=accepted.timeLayer(),expected=data_->source.timeLayer();
    require(actual.step==expected.step&&actual.begin==expected.begin&&actual.end==expected.end
        &&actual.position==expected.position,"cell requirements have a different full source clock");
    require(data_->source.sameIdentity(accepted),"cell requirements have a different actual source owner");
    validateRequirementsSource(mesh,data_->source);
    require(data_->levels.size()==static_cast<std::size_t>(mesh.cellIndexNum()),"cell requirements have a different native leaf domain");
    const auto& active=data_->source.phaseGeometry().vofPublication().domain().cellRowByIndex();
    for(std::size_t ci=0;ci<data_->levels.size();++ci){const int level=data_->levels[ci];
        require(level>=policy.minLevel&&level<=policy.maxLevel,"cell requirement outside active regrid min/max levels");
        if(active[ci]<0)require(level==policy.minLevel,"excluded cell requirement must equal policy minLevel");
    }
}
HydrodynamicsVofRegridDecision prepareHydrodynamicsVofRegrid(
    const QuadtreeMesh& source,const LAMP::Hydrodynamics::HydroTwoPhaseState& accepted,
    const HydrodynamicsVofRegridPolicy& policy,const HydrodynamicsVofRegridRequirements* requirements,
    const HydrodynamicsVofTargetDomainProvider& targetDomainProvider)
{
    const auto started=WallClock::now();
    const int liquidMinimum=policy.liquidMinLevel.value_or(policy.minLevel),gasMinimum=policy.gasMinLevel.value_or(policy.minLevel);
    require(liquidMinimum>=policy.minLevel&&liquidMinimum<=policy.maxLevel
        &&gasMinimum>=policy.minLevel&&gasMinimum<=policy.maxLevel
        &&policy.minLevel>=0&&policy.minLevel<=policy.maxLevel&&policy.maxLevel<=source.maxLevel()
        &&policy.interfaceBandLayers>=0&&std::isfinite(policy.vorticityTolerance)&&policy.vorticityTolerance>0
        &&std::isfinite(policy.coarsenRatio)&&policy.coarsenRatio>0&&policy.coarsenRatio<.5&&policy.interval>0,
        "invalid min/max/allocation, interface halo, vorticity hysteresis or interval policy");
    if(requirements)requirements->validateFor(source,accepted,policy);
    const auto* minimumLevels=requirements?&requirements->minimumLevels():nullptr;
    HydrodynamicsVofRegridReport report;report.sourceClock=accepted.timeLayer();
    report.sourceLeaves=report.targetLeaves=source.cellIndexNum();
    const auto indicatorStarted=WallClock::now();
    const auto indicators=LAMP::Hydrodynamics::HydroVofRefinement::sample(source,accepted);
    report.indicatorSeconds=seconds(indicatorStarted);report.maximumHAbsVorticity=indicators.maximumHAbsVorticity;
    const auto planningStarted=WallClock::now();
    const auto& publication=accepted.phaseGeometry().vofPublication();const auto& domain=publication.domain();
    const auto& active=domain.cellRowByIndex();
    const bool selected=domain.cellIndexByRow().size()!=static_cast<std::size_t>(source.cellIndexNum());
    require(!selected||bool(targetDomainProvider),"selected source requires an explicit target domain provider");
    if(targetDomainProvider){const auto certified=targetDomainProvider(source);certified.validateMeshIdentity(source);
        require(certified.axisymmetric()==domain.axisymmetric()
            &&certified.cellRowByIndex()==active&&certified.cellIndexByRow()==domain.cellIndexByRow(),
            "target domain provider has different actual source selected maps or measure");
    }
    const auto& stocks=publication.plicCells();
    std::vector<std::uint8_t> band(stocks.size());std::vector<int> pure(stocks.size(),-1);
    for(std::size_t ci=0;ci<stocks.size();++ci){
        if(active[ci]<0)continue; // Excluded solid has no phase inventory or interface seed.
        const auto& cell=stocks[ci];const auto& q=cell.inventory;
        if(q[0]>0&&q[1]>0){
            ++report.mixedCells;if(!cell.resolved)++report.unresolvedMixedCells;
        }else require((q[0]>0&&q[1]==0)||(q[1]>0&&q[0]==0),
            "empty or invalid independent phase inventory");
        // A retained unresolved stock is not a resolved geometric interface.
        // No threshold, clipping, transfer-law or inventory change occurs here.
        if(cell.resolved&&cell.measured[0]>0&&cell.measured[1]>0){
            band[ci]=1;++report.geometricMixedCells;
        }
        // Actual majority-pure planes remain valid pure-jump evidence even
        // when their independent minor inventory is nonzero and unresolved.
        if(cell.measured[0]>0&&cell.measured[1]==0)pure[ci]=0;
        else if(cell.measured[1]>0&&cell.measured[0]==0)pure[ci]=1;
    }
    // Protect true geometric pure L/G jumps even without a resolved mixed cell.
    // Native coarse/fine subfaces identify both real neighboring leaves.
    for(const auto face:domain.faces(source))
        if(face.negativeCell>=0&&face.positiveCell>=0&&pure[face.negativeCell]>=0
            &&pure[face.positiveCell]>=0&&pure[face.negativeCell]!=pure[face.positiveCell]){
            band[face.negativeCell]=band[face.positiveCell]=1;++report.pureJumpFaces;
        }
    report.interfaceSeedCells=std::count(band.begin(),band.end(),std::uint8_t{1});
    if(!selected)source.expandCellMask(band,policy.interfaceBandLayers);
    else for(int layer=0;layer<policy.interfaceBandLayers;++layer){
        source.expandCellMask(band,1);
        for(std::size_t ci=0;ci<band.size();++ci)if(active[ci]<0)band[ci]=0;
    }
    report.interfaceBandCells=std::count(band.begin(),band.end(),std::uint8_t{1});
    std::vector<int> desired(stocks.size());CellField<double> oldTargets(source,0.);
    const long double low=static_cast<long double>(policy.coarsenRatio)*policy.vorticityTolerance;
    for(int ci=0;ci<source.cellIndexNum();++ci){
        if(active[ci]<0){desired[ci]=policy.minLevel;oldTargets.value[ci]=policy.minLevel;continue;}
        const int floor=pure[ci]==0?liquidMinimum:pure[ci]==1?gasMinimum:std::max(liquidMinimum,gasMinimum);
        const int current=source.cells(source.cellIndexToLabel[ci]).level();int next=current;
        if(band[ci])next=policy.maxLevel;
        else if(indicators.hAbsVorticity[ci]>policy.vorticityTolerance){
            next=std::clamp(current+1,floor,policy.maxLevel);
            if(next>current)++report.flowRefineCells;
        }else if(indicators.hAbsVorticity[ci]<low){
            next=std::clamp(current-1,floor,policy.maxLevel);
            if(next<current)++report.flowCoarsenCells;
        }
        desired[ci]=std::clamp(next,floor,policy.maxLevel);
        if(minimumLevels&&(*minimumLevels)[ci]>desired[ci]){
            desired[ci]=(*minimumLevels)[ci];
            if(desired[ci]>current)++report.additionalRefineCells;else ++report.additionalHoldCells;
        }
        oldTargets.value[ci]=desired[ci];
    }
    // Complete only this temporary target wave; accepted phase/flow and
    // the source topology remain immutable during planning.
    source.completeRefinementPrerequisites(oldTargets);
    report.requested=source.mayAdaptToTargets(oldTargets);report.planningSeconds=seconds(planningStarted);
    if(!report.requested)for(int ci=0;ci<source.cellIndexNum();++ci)
        require(source.cells(source.cellIndexToLabel[ci]).level()>=desired[ci],
            "native preflight stopped before an original requested refinement level");
    const auto done=[&](std::optional<HydrodynamicsVofGridCandidate> candidate){
        report.totalSeconds=seconds(started);return HydrodynamicsVofRegridDecision{std::move(candidate),report};
    };
    // False preflight creates zero clone/overlap/physical transfer/solve.
    if(!report.requested)return done(std::nullopt);
    const auto adaptationStarted=WallClock::now();
    auto& preparation=report.treePreparation;preparation.available=true;
    const auto sourceStorage=source.entityPoolStorageReport();
    const auto cloneStarted=treeTimingStart();auto candidate=source.clone();++report.clones;
    recordTreeTiming(cloneStarted,preparation.cloneSeconds,preparation.cloneCalls,preparation.available);
    const auto cloneStorage=candidate->entityPoolStorageReport();
    for(std::size_t pool=0;pool<4;++pool){
        report.cloneStorage.logicalCapacity[pool]=sourceStorage[pool].logicalCapacity;
        report.cloneStorage.sourceConstructed[pool]=sourceStorage[pool].constructedSlots;
        report.cloneStorage.cloneConstructed[pool]=cloneStorage[pool].constructedSlots;
    }
    report.cloneStorage.available=true;
    const int maximumSweeps=source.maxLevel()+2;
    std::vector<LeafSignature> visited;
    for(int sweep=0;sweep<=maximumSweeps;++sweep){
        NativeLeafDomain targetDomain;
        if(targetDomainProvider){const auto domainStarted=treeTimingStart();
            targetDomain=targetDomainProvider(*candidate);targetDomain.validateMeshIdentity(*candidate);
            require(targetDomain.axisymmetric()==domain.axisymmetric(),"target domain provider changed physical measure");
            // Includes the provider's actual CSG certificate and selected graph.
            recordTreeTiming(domainStarted,preparation.domainSeconds,preparation.domainCalls,preparation.available);
        }else{
            FaceSegmentGraph graph;const auto graphStarted=treeTimingStart();graph.build(*candidate);
            recordTreeTiming(graphStarted,preparation.graphSeconds,preparation.graphCalls,preparation.available);
            const auto domainStarted=treeTimingStart();
            targetDomain=NativeLeafDomain::build(*candidate,graph,domain.axisymmetric());
            recordTreeTiming(domainStarted,preparation.domainSeconds,preparation.domainCalls,preparation.available);
        }
        const auto overlapStarted=treeTimingStart();
        const auto overlap=LAMP::Mesh::NativeLeafOverlap::build(source,domain,*candidate,targetDomain);
        recordTreeTiming(overlapStarted,preparation.overlapSeconds,preparation.overlapCalls,preparation.available);
        LAMP::Mesh::NativeTargetLevelField targets(*candidate,0.),gradient(*candidate,0.);
        for(int ci=0;ci<candidate->cellIndexNum();++ci){
            int level=policy.minLevel;
            // Max over ALL old donors protects fine siblings that a parent
            // center query could miss when deciding a coarse target.
            for(const auto index:overlap.byTarget()[ci])
                level=std::max(level,desired[overlap.pieces()[index].sourceCell]);
            targets.value[ci]=level;
        }
        candidate->completeRefinementPrerequisites(targets);
        if(!candidate->mayAdaptToTargets(targets)){
            // Balancing can retain fine leaves instead of coarsening. It must
            // never silently leave an overlapping original refinement unmet.
            for(const auto& piece:overlap.pieces())
                require(candidate->cells(candidate->cellIndexToLabel[piece.targetCell]).level()>=desired[piece.sourceCell],
                    "settled native mesh did not reach an original requested refinement level");
            report.changed=std::any_of(overlap.pieces().begin(),overlap.pieces().end(),
                [](const auto& q){return q.kind!=LAMP::Mesh::NativeLeafOverlapPiece::Kind::Unchanged;});
            // Excluded leaves are still real background topology. A solid-only
            // change must not be discarded and retried at every later step.
            if(selected&&!report.changed)report.changed=leafSignature(source)!=leafSignature(*candidate);
            report.targetLeaves=candidate->cellIndexNum();report.adaptationSeconds=seconds(adaptationStarted);
            finishTreeTiming(report);
            if(!report.changed)return done(std::nullopt);
            return done(HydrodynamicsVofGridCandidate{std::move(candidate),std::move(targetDomain)});
        }
        // Preflight is conservative. Only actual dyadic topology change is
        // progress; a rebuilt frozen target plus false preflight above is the
        // independent certificate for a stable mesh (including consumed votes).
        auto signature=leafSignature(*candidate);
        if(!visited.empty()&&signature==visited.back())
            throw std::invalid_argument("[Hydrodynamics VOF regrid] native target adaptation made no topology progress at sweep="+std::to_string(sweep)+targetCounts(*candidate,targets));
        if(std::find(visited.begin(),visited.end(),signature)!=visited.end())
            throw std::invalid_argument("[Hydrodynamics VOF regrid] native target adaptation repeated actual topology at sweep="+std::to_string(sweep)+targetCounts(*candidate,targets));
        visited.push_back(std::move(signature));
        require(sweep<maximumSweeps,"native target adaptation did not settle within the original allocated-level bound");
        std::vector<VertField<double>> levelSets;std::vector<VertField<double>*> vertices;
        std::vector<HEdgeField<double>*> horizontal;std::vector<VEdgeField<double>*> vertical;
        std::vector<LAMP::Mesh::NativeTargetLevelField*> cells{&targets};
        const auto nativeAdaptStarted=treeTimingStart();
        candidate->adapt(levelSets,gradient,vertices,horizontal,vertical,cells,&targets);++report.sweeps;
        recordTreeTiming(nativeAdaptStarted,preparation.nativeAdaptSeconds,preparation.nativeAdaptCalls,preparation.available);
    }
    throw std::logic_error("[Hydrodynamics VOF regrid] unreachable adaptation state");
}
