#include "Modules/Interface/HermiteGalsPhaseRegion.h"
#include "Modules/Interface/HermiteGalsPublication.h"
#include "Modules/Interface/Detail/HermiteGalsRegionMeasure.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Geometry/Detail/ExactExpansion.h"
#include <mutex>
#include <optional>
#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace LAMP::Interface {
namespace {
using I=LevelSetMappedInterface;
using Snapshot=LevelSetHermiteGeometry::Snapshot;
constexpr long double pi=3.1415926535897932384626433832795L;
void require(bool ok,const std::string& message){if(!ok)throw std::invalid_argument("[HermiteGalsPhaseRegion] "+message);}
struct Sum {long double value=0,correction=0;void add(long double x){const auto a=x-correction,b=value+a;correction=(b-value)-a;value=b;}};
long double capacity(Vector2D a,Vector2D b,bool radial){return (static_cast<long double>(b.x)-a.x)*(static_cast<long double>(b.y)-a.y)*
    (radial?pi*(static_cast<long double>(b.x)+a.x):1);}
using Box=Detail::RegionMeasure::Box;
// Exact ordering/coverage of represented rectangle endpoints. No area-sum
// tolerance or rounded interior sample can conceal a small inactive gap.
void coverage(Vector2D lo,Vector2D hi,const std::vector<Box>& boxes){
    std::vector<double> xs{lo.x,hi.x};for(const auto& b:boxes){xs.push_back(b.lower.x);xs.push_back(b.upper.x);}
    std::sort(xs.begin(),xs.end());xs.erase(std::unique(xs.begin(),xs.end()),xs.end());
    for(std::size_t k=1;k<xs.size();++k){std::vector<std::pair<double,double>> ys;
        for(const auto& b:boxes)if(b.lower.x<=xs[k-1] && b.upper.x>=xs[k])ys.emplace_back(b.lower.y,b.upper.y);
        std::sort(ys.begin(),ys.end());double reach=lo.y;
        for(const auto& y:ys){require(y.first==reach,"rectangle has a gap or overlapping active leaves");reach=y.second;}
        require(reach==hi.y,"rectangle leaves active scalar support");}
}
// Per-native-cell execution slots. Only real cut integrations write these;
// relaxed atomics also cover concurrent independent fallback queries.
struct CutWork {
    std::atomic<std::uint64_t> nanoseconds{0};
    std::atomic<std::size_t> calls{0};
};
class CutWorkTimer {
public:
    explicit CutWorkTimer(CutWork& work):work_(work),begin_(Clock::now()){}
    ~CutWorkTimer()noexcept{
        const auto elapsed=std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now()-begin_).count();
        work_.nanoseconds.fetch_add(static_cast<std::uint64_t>(elapsed),std::memory_order_relaxed);
        work_.calls.fetch_add(1,std::memory_order_relaxed);
    }
private:
    using Clock=std::chrono::steady_clock;
    CutWork& work_;
    Clock::time_point begin_;
};
struct PairDetailWork {
    std::atomic<std::uint64_t> regionGaugeNanoseconds{0},cutPreparationNanoseconds{0},coreNanoseconds{0},
        retryNanoseconds{0},assemblyNanoseconds{0};
    std::atomic<std::size_t> cutConstructions{0},coreAttempts{0},faceCutRefinements{0};
    // One writer under this actual group's call_once; read after queries join.
    I::NativeMomentWorkProfile coreWork;
};
// Commit once when an actual pair finishes or unwinds. The measured work stays
// local throughout every cut/Core/retry; only the native-cell slot is updated.
class PairDetailCommit {
public:
    PairDetailCommit(PairDetailWork& target,const Detail::RegionMeasure::HalfPairWorkProfile& local)
        :target_(target),local_(local){}
    ~PairDetailCommit()noexcept{
        target_.regionGaugeNanoseconds.fetch_add(local_.regionGaugeNanoseconds,std::memory_order_relaxed);
        target_.cutPreparationNanoseconds.fetch_add(local_.cutPreparationNanoseconds,std::memory_order_relaxed);
        target_.coreNanoseconds.fetch_add(local_.coreNanoseconds,std::memory_order_relaxed);
        target_.retryNanoseconds.fetch_add(local_.retryNanoseconds,std::memory_order_relaxed);
        target_.assemblyNanoseconds.fetch_add(local_.assemblyNanoseconds,std::memory_order_relaxed);
        target_.cutConstructions.fetch_add(local_.cutConstructions,std::memory_order_relaxed);
        target_.coreAttempts.fetch_add(local_.coreAttempts,std::memory_order_relaxed);
        target_.faceCutRefinements.fetch_add(local_.faceCutRefinements,std::memory_order_relaxed);
        Detail::RegionMeasure::accumulateCoreWorkProfile(target_.coreWork,local_.coreWork);
    }
private:
    PairDetailWork& target_;
    const Detail::RegionMeasure::HalfPairWorkProfile& local_;
};
using CellMeasure=Detail::RegionMeasure::CellMeasure;
using Detail::RegionMeasure::integrateCell;
}
struct HermiteGalsPhaseRegionSource::Data {
    explicit Data(const HermiteGalsState& s,const HermiteGalsPhaseRegionConfig& c):source(s),config(c){}
    HermiteGalsState source;
    HermiteGalsPhaseRegionConfig config;
    std::shared_ptr<const Snapshot> scalar;
    LevelSetInterfacePanels::NativeComponents ownedTopology;
    std::optional<HermiteGalsPublication> publication;
    const LevelSetInterfacePanels::NativeComponents& topology() const
    {return publication?publication->topology():ownedTopology;}
    std::vector<std::vector<I::Graph>> graphs;
    std::size_t nativeSideGraphs=0;
    struct HalfBatch {
        std::array<Box,4> boxes;
        mutable std::array<std::once_flag,2> ready;
        mutable std::array<std::optional<std::array<CellMeasure,2>>,2> moments;
        mutable std::array<CutWork,2> pairWork;
        mutable std::array<PairDetailWork,2> pairDetailWork;
        mutable CutWork fallbackWork;
    };
    std::vector<std::unique_ptr<HalfBatch>> halfBatches;
    CellMeasure cutMeasure(const Box& box,bool radial) const {
        LAMP::Geometry::Detail::ExactExpansion::validateEnvironment();
        const auto& batch=*halfBatches.at(box.cell);int slot=-1;
        for(int k=0;k<4;++k){const auto& b=batch.boxes[k];
            if(b.lower.x==box.lower.x && b.lower.y==box.lower.y && b.upper.x==box.upper.x && b.upper.y==box.upper.y)slot=k;}
        if(slot>=0){
            const int group=slot/2;
            std::call_once(batch.ready[group],[&] {
                try{batch.moments[group]=[&]{
                    const CutWorkTimer timer(batch.pairWork[group]);
                    Detail::RegionMeasure::HalfPairWorkProfile work;
                    const PairDetailCommit commit(batch.pairDetailWork[group],work);
                    return Detail::RegionMeasure::integrateHalfPair(scalar,
                        {{batch.boxes[2*group],batch.boxes[2*group+1]}},graphs[box.cell],radial,&work);
                }();}
                catch(...){ // Additional batch work cannot replace a requested single-query failure.
                    batch.moments[group].reset();
                }
            });
            if(batch.moments[group]){const auto& m=(*batch.moments[group])[slot%2];
                const auto cap=capacity(box.lower,box.upper,radial);
                const auto defect=std::max({std::abs(m.axialLiquid-m.radialLiquid),std::abs(m.axialGas-m.radialGas),
                    std::abs(m.axialLiquid+m.axialGas-cap),std::abs(m.radialLiquid+m.radialGas-cap)});
                bool valid=std::isfinite(defect) && defect<=config.relativeLedgerTolerance*cap
                    && std::isfinite(m.error) && m.error>=0 && m.error<=config.relativeMomentEstimateTolerance*cap;
                for(long double v:{m.axialLiquid,m.radialLiquid,m.axialGas,m.radialGas})valid=valid && std::isfinite(v) && v>=0 && v<=cap;
                valid=valid && std::isfinite(m.fixedAxialAssemblyDefect) && std::isfinite(m.fixedRadialAssemblyDefect)
                    && m.assemblyCapacity>0 && std::max(std::abs(m.fixedAxialAssemblyDefect),std::abs(m.fixedRadialAssemblyDefect))
                        <=config.relativeLedgerTolerance*m.assemblyCapacity;
                if(valid)return m;
            }
        }
        // Generic cuts, unsupported batches or a target outside the original
        // gates use that target's original independent integration and budget.
        const CutWorkTimer timer(batch.fallbackWork);
        return integrateCell(scalar,box,graphs[box.cell],radial);
    }
};
HermiteGalsPhaseRegionSource::HermiteGalsPhaseRegionSource(const QuadtreeMesh& mesh,const HermiteGalsState& source,
    const HermiteGalsPhaseRegionConfig& config)
    :HermiteGalsPhaseRegionSource(mesh,source,config,nullptr){}
HermiteGalsPhaseRegionSource::HermiteGalsPhaseRegionSource(const QuadtreeMesh& mesh,const HermiteGalsPublication& publication,
    const HermiteGalsPhaseRegionConfig& config)
    :HermiteGalsPhaseRegionSource(mesh,publication.state(),config,&publication){}
HermiteGalsPhaseRegionSource::HermiteGalsPhaseRegionSource(const QuadtreeMesh& mesh,const HermiteGalsState& source,
    const HermiteGalsPhaseRegionConfig& config,const HermiteGalsPublication* publication){
    source.validateMesh(mesh);
    for(double v:{config.relativeLedgerTolerance,config.relativeMomentEstimateTolerance})
        require(std::isfinite(v) && v>0 && v<1,"invalid region relative tolerance");
    auto data=std::make_shared<Data>(source,config);auto view=std::make_shared<const HermiteGalsGeometryView>(source.geometry());
    data->scalar=std::shared_ptr<const Snapshot>(view,&view->snapshot());const auto& scalar=*data->scalar;
    require(scalar.usesNativeLeaves() &&
        (scalar.reconstruction()==LevelSetHermiteGeometry::Reconstruction::EnrichedQuartic || scalar.usesC2Macro()),
        "region requires a native quartic scalar");
    if(scalar.usesC2Macro())require(source.primary().transport==HermiteGalsTransport::FirstJet
        && source.primary().phaseVolumes.empty()
        && source.primary().measure==HermiteGalsMeasure::CartesianArea
        && source.primary().geometryCells.size()==static_cast<std::size_t>(mesh.cellIndexNum()),
        "C2 region requires complete XY FirstJet without transported inventory");
    const bool radial=source.primary().measure==HermiteGalsMeasure::AxisymmetricVolume;
    require(!scalar.radialAxisSymmetry() || radial,"radial scalar with Cartesian region measure");
    const bool nativeTrace=source.primary().transport==HermiteGalsTransport::FirstJet;
    if(publication){
        publication->validateMesh(mesh);
        const auto& published=publication->state();
        const auto clock=source.timeLayer(),publishedClock=published.timeLayer();
        require(published.identity()==source.identity() && &published.geometry().snapshot()==data->scalar.get()
            && clock.step==publishedClock.step && clock.begin==publishedClock.begin && clock.end==publishedClock.end
            && clock.position==publishedClock.position,"different retained publication scalar/clock");
        std::vector<Vector2D> contacts;
        std::vector<LevelSetInterfacePanels::SymmetryPlane> planes;
        const auto samePoint=[](Vector2D a,Vector2D b){return a.x==b.x && a.y==b.y;};
        for(const auto& chain:publication->topology().curves){
            for(Vector2D point:chain.boundaryContacts())
                if(std::none_of(contacts.begin(),contacts.end(),[&](Vector2D p){return samePoint(p,point);}))contacts.push_back(point);
            for(const auto& plane:chain.symmetryPlanes())
                if(std::find(planes.begin(),planes.end(),plane)==planes.end())planes.push_back(plane);
        }
        require(config.boundaryContacts.size()==contacts.size() && std::all_of(contacts.begin(),contacts.end(),[&](Vector2D point){
            return std::any_of(config.boundaryContacts.begin(),config.boundaryContacts.end(),[&](Vector2D p){return samePoint(p,point);});}),
            "region contacts differ from retained publication declarations");
        require(config.symmetryPlanes.size()==planes.size() && std::all_of(planes.begin(),planes.end(),[&](const auto& plane){
            return std::find(config.symmetryPlanes.begin(),config.symmetryPlanes.end(),plane)!=config.symmetryPlanes.end();}),
            "region symmetry planes differ from retained publication declarations");
        data->publication=*publication;
    }else{
        data->ownedTopology=nativeTrace?LevelSetInterfacePanels::discoverNativeTraceComponents(mesh,data->scalar,radial,config.boundaryContacts,config.symmetryPlanes)
            :LevelSetInterfacePanels::discoverNativeComponents(mesh,data->scalar,radial,config.boundaryContacts,config.symmetryPlanes);
    }
    data->graphs.resize(mesh.cellIndexNum());
    for(const auto& chain:data->topology().curves){require(chain.usesSource(scalar) && chain.axisymmetric()==radial,"foreign region component");
        for(const auto& graph:chain.graphs()){require(graph.nativeCellIndex>=0 && graph.nativeCellIndex<mesh.cellIndexNum(),"foreign native region graph");
            if(graph.nativeSide) {
                require(nativeTrace && graph.nativeSide->usesSource(scalar),"foreign/legacy region side support");
                for(const auto* side:{&graph.nativeSide->minus(),&graph.nativeSide->plus()})if(side->available())
                    require((side->phase==-1 || side->phase==1) && data->topology().cellPhase.at(side->nativeCellIndex)==side->phase,"region side lacks complete pure interior certificate");
                ++data->nativeSideGraphs;continue; // Boundary geometry adds no interior box volume.
            }
            data->graphs[graph.nativeCellIndex].push_back(graph);}}
    data->halfBatches.resize(mesh.cellIndexNum());
    for(int cell=0;cell<mesh.cellIndexNum();++cell)if(scalar.nativeCellActive(cell)){
        auto batch=std::make_unique<Data::HalfBatch>();
        const auto center=mesh.cells(mesh.cellIndexToLabel[cell]).centerCoord();
        batch->boxes=Detail::RegionMeasure::halfBoxes(scalar,cell,center);
        data->halfBatches[cell]=std::move(batch);
    }
    data_=std::move(data);
}
const HermiteGalsPhaseRegionSource::Data& HermiteGalsPhaseRegionSource::checked()const{require(bool(data_),"moved-from region source");return *data_;}
const HermiteGalsState& HermiteGalsPhaseRegionSource::source()const{return checked().source;}
std::size_t HermiteGalsPhaseRegionSource::components()const{return checked().topology().curves.size();}
std::size_t HermiteGalsPhaseRegionSource::nativeSideGraphs()const{return checked().nativeSideGraphs;}
HermiteGalsPhaseRegionWorkProfile HermiteGalsPhaseRegionSource::workProfile()const{
    HermiteGalsPhaseRegionWorkProfile profile;
    const auto add=[&](const CutWork& work,bool paired){
        const auto seconds=static_cast<double>(work.nanoseconds.load(std::memory_order_relaxed))*1e-9;
        const auto calls=work.calls.load(std::memory_order_relaxed);
        if(paired){profile.pairWorkerSeconds+=seconds;profile.pairCalls+=calls;}
        else{profile.fallbackWorkerSeconds+=seconds;profile.fallbackCalls+=calls;}
    };
    for(const auto& batch:checked().halfBatches)if(batch){
        for(std::size_t group=0;group<batch->pairWork.size();++group){
            add(batch->pairWork[group],true);const auto& work=batch->pairDetailWork[group];
            profile.pairRegionGaugeWorkerSeconds+=static_cast<double>(work.regionGaugeNanoseconds.load(std::memory_order_relaxed))*1e-9;
            profile.pairCutPreparationWorkerSeconds+=static_cast<double>(work.cutPreparationNanoseconds.load(std::memory_order_relaxed))*1e-9;
            profile.pairCoreWorkerSeconds+=static_cast<double>(work.coreNanoseconds.load(std::memory_order_relaxed))*1e-9;
            profile.pairRetryWorkerSeconds+=static_cast<double>(work.retryNanoseconds.load(std::memory_order_relaxed))*1e-9;
            profile.pairAssemblyWorkerSeconds+=static_cast<double>(work.assemblyNanoseconds.load(std::memory_order_relaxed))*1e-9;
            profile.pairCutConstructions+=work.cutConstructions.load(std::memory_order_relaxed);
            profile.pairCoreAttempts+=work.coreAttempts.load(std::memory_order_relaxed);
            profile.pairFaceCutRefinements+=work.faceCutRefinements.load(std::memory_order_relaxed);
            Detail::RegionMeasure::accumulateCoreWorkProfile(profile.pairCoreWork,work.coreWork);
        }
        add(batch->fallbackWork,false);
    }
    return profile;
}
void HermiteGalsPhaseRegionSource::validateMesh(const QuadtreeMesh& mesh)const{checked().source.validateMesh(mesh);}
HermiteGalsPhaseRegionMeasure HermiteGalsPhaseRegionSource::measure(const QuadtreeMesh& mesh,Vector2D lo,Vector2D hi)const{
    validateMesh(mesh);return measureValidatedMesh(mesh,lo,hi);
}
std::vector<HermiteGalsPhaseRegionBatchMeasure> HermiteGalsPhaseRegionSource::measureBatch(const QuadtreeMesh& mesh,
    const std::vector<std::array<Vector2D,2>>& orderedBoxes)const{
    validateMesh(mesh);
    std::vector<HermiteGalsPhaseRegionBatchMeasure> results(orderedBoxes.size());
    const auto measureOne=[&](std::size_t q){
        try{results[q].measure=measureValidatedMesh(mesh,orderedBoxes[q][0],orderedBoxes[q][1]);}
        catch(...){results[q].failure=std::current_exception();}
    };
    if(!orderedBoxes.empty()){
        measureOne(0);
        if(results[0].failure)std::rethrow_exception(results[0].failure);
    }
    // Later query failures remain beside their results: the material consumer
    // interleaves their rethrow with its original ordered assembly checks.
    // No final validation can overtake an earlier query/assembly failure.
    const auto queryCount=static_cast<std::ptrdiff_t>(orderedBoxes.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1) if(orderedBoxes.size()>1 && !omp_in_parallel() && omp_get_max_threads()>1)
#endif
    for(std::ptrdiff_t q=1;q<queryCount;++q)measureOne(static_cast<std::size_t>(q));
    return results;
}
HermiteGalsPhaseRegionMeasure HermiteGalsPhaseRegionSource::measureValidatedMesh(const QuadtreeMesh& mesh,Vector2D lo,Vector2D hi)const{
    LAMP::Geometry::Detail::ExactExpansion::validateEnvironment();
    const auto& data=checked();const auto& scalar=*data.scalar;
    const bool radial=data.source.primary().measure==HermiteGalsMeasure::AxisymmetricVolume;
    require(std::isfinite(lo.x) && std::isfinite(lo.y) && std::isfinite(hi.x) && std::isfinite(hi.y) && lo.x<hi.x && lo.y<hi.y && (!radial || lo.x>=0),"invalid rectangle bounds");
    std::vector<Box> boxes;
    for(int cell:scalar.nativeCellCandidates(lo,hi)){const auto b=scalar.nativeCellBounds(cell);
        const Vector2D a{std::max(lo.x,b[0].x),std::max(lo.y,b[0].y)},z{std::min(hi.x,b[1].x),std::min(hi.y,b[1].y)};
        if(a.x>=z.x || a.y>=z.y)continue;require(scalar.nativeCellActive(cell),"rectangle intersects inactive scalar support");
        boxes.push_back({a,z,cell});}
    coverage(lo,hi,boxes);
    HermiteGalsPhaseRegionMeasure result;result.sourceIdentity=data.source.identity();result.lower=lo;result.upper=hi;
    result.capacity=capacity(lo,hi,radial);require(std::isfinite(result.capacity) && result.capacity>0,"unrepresentable rectangle capacity");
    Sum axialLiquid,radialLiquid,axialGas,radialGas,negative,positive,error;bool allNegative=true,allPositive=true;
    for(const auto& box:boxes){++result.intersectedLeaves;const auto cap=capacity(box.lower,box.upper,radial);
        require(std::isfinite(cap) && cap>0,"unrepresentable leaf intersection capacity");
        int phase=data.topology().cellPhase.at(box.cell);if(phase==0)phase=scalar.nativeQuarticBoxSign(box.cell,box.lower,box.upper);
        CellMeasure m;
        if(phase!=0){require(phase==-1 || phase==1,"invalid phase certificate");
            m.axialLiquid=m.radialLiquid=phase<0?cap:0;m.axialGas=m.radialGas=phase>0?cap:0;}
        else{require(!data.graphs[box.cell].empty(),"uncertified graph-free leaf");++result.cutLeaves;
            m=data.cutMeasure(box,radial);}
        allNegative=allNegative && phase<0;allPositive=allPositive && phase>0;
        for(long double value:{m.axialLiquid,m.radialLiquid,m.axialGas,m.radialGas,m.error,
            m.fixedAxialAssemblyDefect,m.fixedRadialAssemblyDefect})require(std::isfinite(value),"nonfinite region moment");
        // The new inside-trace clip ledger includes the actual coincident
        // interface term. It is independent of the legacy shared-face ledger
        // and imposes no material-flux or per-cell inventory constraint.
        if(m.assemblyCapacity>0)require(std::max(std::abs(m.fixedAxialAssemblyDefect),std::abs(m.fixedRadialAssemblyDefect))
            <=data.config.relativeLedgerTolerance*m.assemblyCapacity,"region inside-trace fixed/coincident ledger defect");
        const auto defect=std::max({std::abs(m.axialLiquid-m.radialLiquid),std::abs(m.axialGas-m.radialGas),
            std::abs(m.axialLiquid+m.axialGas-cap),std::abs(m.radialLiquid+m.radialGas-cap)})/cap;
        result.maxRelativeLedgerDefect=std::max(result.maxRelativeLedgerDefect,defect);
        require(defect<=data.config.relativeLedgerTolerance,"region Green/capacity ledger defect");
        const bool axial=box.upper.y-box.lower.y<=box.upper.x-box.lower.x;
        const auto liquid=axial?m.axialLiquid:m.radialLiquid,gas=axial?m.axialGas:m.radialGas;
        require(liquid>=0 && liquid<=cap && gas>=0 && gas<=cap,"region phase volume outside strict capacity; no clipping");
        const auto phases=Detail::RegionMeasure::complement(cap,liquid,gas);
        result.leaves.push_back({box.cell,box.lower,box.upper,cap,phases[0],phases[1]});
        negative.add(liquid);positive.add(gas);axialLiquid.add(m.axialLiquid);radialLiquid.add(m.radialLiquid);
        axialGas.add(m.axialGas);radialGas.add(m.radialGas);error.add(m.error);result.integratedArcs+=m.arcs;result.rootRefinements+=m.refinements;
    }
    result.rawAxialLiquid=axialLiquid.value;result.rawRadialLiquid=radialLiquid.value;result.rawAxialGas=axialGas.value;result.rawRadialGas=radialGas.value;
    result.momentErrorEstimate=error.value;
    require(error.value>=0 && error.value<=data.config.relativeMomentEstimateTolerance*result.capacity,"region moment error estimate exceeds tolerance");
    require(std::abs(negative.value+positive.value-result.capacity)<=data.config.relativeLedgerTolerance*result.capacity,"region assembled capacity defect");
    if(allNegative){result.liquid=result.capacity;result.gas=0;}
    else if(allPositive){result.liquid=0;result.gas=result.capacity;}
    else{const auto phases=Detail::RegionMeasure::complement(result.capacity,negative.value,positive.value);result.liquid=phases[0];result.gas=phases[1];}
    require(std::isfinite(result.liquid) && std::isfinite(result.gas) && result.liquid>=0 && result.gas>=0 && result.liquid<=result.capacity && result.gas<=result.capacity,"assembled region outside strict capacity");
    return result;
}
}
