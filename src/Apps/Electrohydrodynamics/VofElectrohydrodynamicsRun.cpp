#include "Electrohydrodynamics/VofElectrohydrodynamicsRun.h"
#include "Electrohydrodynamics/IO/VofElectrohydrodynamicsOutput.h"
#include "Hydrodynamics/HydrodynamicsVofRegrid.h"
#include "Modules/Hydrodynamics/HydroVofMeshTransfer.h"
#include "Modules/Hydrodynamics/HydroMaxwellLoad.h"
#include "Modules/Electric/ElectricVofChargeTransport.h"
#include "Modules/Electric/ElectricVofRefinement.h"
#include "Modules/Electric/ElectricVofMaxwellStress.h"
#include "Modules/Electric/ElectricSurfaceChargeTransport.h"
#include "Modules/Electric/ElectricVofSharpMaxwellStress.h"
#include "Core/Geometry/Detail/FiniteVolumeClosure.h"
#include "Core/IO/BinaryArchive.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/NativeMeshIdentity.h"
#include "Core/Geometry/FittedDomainGeometry.h"
#include <algorithm>
#include <array>
#include <cfenv>
#include <cerrno>
#include <exception>
#include <cstring>
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
#include <float.h>
#include <xmmintrin.h>
#endif
#include <cmath>
#include <vector>
#include <optional>
#include <stdexcept>
#include <utility>

namespace {
namespace H=LAMP::Hydrodynamics;
namespace E=LAMP::Electric;
namespace FV=LAMP::Geometry::Detail::FiniteVolume;
using Vof=LAMP::Interface::GeometricVofPublication;
using EhdClock=LAMP::TimeLayer;
// This is a business execution cache, not an observer. Only an unchanged
// COMPLETE floating state can certify skipping the geometry construction.
struct MaxwellGeometryFpState {
    std::array<unsigned char,sizeof(std::fenv_t)> bytes{};
    unsigned csr=0;
    static bool read(MaxwellGeometryFpState& result) noexcept {
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
        unsigned int control=0;
        if(_controlfp_s(&control,0,0)!=0||(control&_MCW_EM)!=_MCW_EM)return false;
        std::fenv_t state;std::memset(&state,0,sizeof(state));
        if(std::fegetenv(&state)!=0)return false;
        result.csr=_mm_getcsr();
        if((result.csr&0x1f80u)!=0x1f80u)return false;
        std::memcpy(result.bytes.data(),&state,sizeof(state));return true;
#else
        (void)result;return false;
#endif
    }
    bool same(const MaxwellGeometryFpState& other) const noexcept
    {return bytes==other.bytes&&csr==other.csr;}
};
// Each timestamp/conversion is passive: never span a business call with an
// environment guard, and never turn an unsupported timer into a physical error.
struct PassiveAmrTimingEnvironment {
    std::fenv_t native{};int error=errno;bool captured=false,eligible=false;
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
    unsigned csr=_mm_getcsr();
#endif
    PassiveAmrTimingEnvironment() noexcept {
        captured=std::fegetenv(&native)==0;
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
        unsigned control=0;
        eligible=captured&&_controlfp_s(&control,0,0)==0
            &&(control&_MCW_EM)==_MCW_EM&&(csr&0x1f80u)==0x1f80u;
#endif
    }
    ~PassiveAmrTimingEnvironment() noexcept {
        if(captured&&std::fesetenv(&native)!=0)std::terminate();
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
        _mm_setcsr(csr);
#endif
        errno=error;
    }
};
#if defined(_MSC_VER)
#pragma float_control(precise,on,push)
#pragma fenv_access(on)
#pragma float_control(except,on)
#endif
using AmrTimingClock=std::chrono::steady_clock;
std::optional<AmrTimingClock::time_point> amrTimingStart() noexcept {
    PassiveAmrTimingEnvironment environment;
    if(!environment.eligible)return std::nullopt;
    return AmrTimingClock::now();
}
std::optional<double> amrTimingFinish(const std::optional<AmrTimingClock::time_point>& started) noexcept {
    if(!started)return std::nullopt;
    PassiveAmrTimingEnvironment environment;
    if(!environment.eligible)return std::nullopt;
    const double seconds=std::chrono::duration<double>(AmrTimingClock::now()-*started).count();
    if(!std::isfinite(seconds)||seconds<0)return std::nullopt;
    return seconds;
}
#if defined(_MSC_VER)
#if !defined(_M_FP_STRICT) || !_M_FP_STRICT
#pragma fenv_access(off)
#endif
#pragma float_control(pop)
#endif
struct RetainedMaxwellGeometry {
    H::HydroStrainOperator geometry;
    MaxwellGeometryFpState environment;
};
struct AttemptMaxwellGeometry {
    explicit AttemptMaxwellGeometry(std::optional<RetainedMaxwellGeometry> source):offered(std::move(source)){}
    std::optional<RetainedMaxwellGeometry> offered,registered;
};
bool maxwellGeometryMode(const H::HydroStrainOperator& geometry,bool greenOnly) noexcept {
    return greenOnly?geometry.isRt0GreenCompatible()&&!geometry.isStressCompatible():geometry.isStressCompatible();
}
void require(bool ok,const char* why){if(!ok)throw std::logic_error(std::string("[VofEHD transaction] ")+why);}
bool sameClock(EhdClock a,EhdClock b){return a.step==b.step&&a.begin==b.begin&&a.end==b.end&&a.position==b.position;}
constexpr std::size_t csgBoundaryByteLimit=32u*1024u*1024u;
struct RetainedCsgBoundaryIds {
    std::weak_ptr<const NativeMeshIdentity> domain;
    std::weak_ptr<const DictionaryGeometry> definition;
    std::weak_ptr<const FittedDomainGeometry::InitialPhaseSampler> sampler;
    std::vector<int> excludedPhases,ids;
    MaxwellGeometryFpState environment;
    int error=0;
};
struct CsgBoundaryCache { std::optional<RetainedCsgBoundaryIds> retained; };
struct CsgBoundaryWork {
    std::shared_ptr<CsgBoundaryCache> cache;
    VofEhdCsgBoundaryStageReport* report=nullptr;
    ~CsgBoundaryWork() noexcept {
        PassiveAmrTimingEnvironment current;cache.reset();
    }
};
// Probes and optional stock management preserve their CURRENT environment;
// the original selector/graph/classification is never inside this guard.
bool csgBoundaryEnvironment(MaxwellGeometryFpState& state,int& error) noexcept {
    PassiveAmrTimingEnvironment current;
    error=current.error;
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
    if(!current.eligible)return false;
    std::memcpy(state.bytes.data(),&current.native,sizeof(current.native));
    state.csr=current.csr;return true;
#else
    (void)state;return false;
#endif
}
template<class T> bool sameCsgOwner(const std::weak_ptr<const T>& retained,
    const std::shared_ptr<const T>& actual) noexcept {
    const auto owner=retained.lock();
    return owner&&owner.get()==actual.get()&&!retained.owner_before(actual)&&!actual.owner_before(retained);
}
std::size_t csgBoundaryRetainedBytes(const RetainedCsgBoundaryIds& retained) noexcept {
    return sizeof(CsgBoundaryCache)+256u+sizeof(int)*(retained.ids.capacity()+retained.excludedPhases.capacity());
}
CsgBoundaryWork csgBoundaryWork(std::shared_ptr<CsgBoundaryCache>& cache,
    VofEhdCsgBoundaryStageReport& report,bool enabled) noexcept {
    if(enabled&&!cache){
        PassiveAmrTimingEnvironment current;
        if(current.eligible)try{cache=std::make_shared<CsgBoundaryCache>();}
        catch(const std::bad_alloc&){}catch(const std::length_error&){}
    }
    return {cache,&report};
}
void dropCsgBoundaryCache(std::shared_ptr<CsgBoundaryCache>& cache) noexcept {
    PassiveAmrTimingEnvironment current;cache.reset();
}
std::optional<std::vector<int>> cachedCsgBoundaryIds(const CsgBoundaryWork& work,
    const NativeLeafDomain& domain,const HydrodynamicsInput& input,
    const MaxwellGeometryFpState& environment,int error) noexcept {
    PassiveAmrTimingEnvironment current;
    if(!current.eligible||!work.cache||!work.cache->retained)return std::nullopt;
    auto& retained=*work.cache->retained;
    if(!sameCsgOwner(retained.domain,domain.meshIdentity())
        ||!sameCsgOwner(retained.definition,input.csgGeometry)
        ||!sameCsgOwner(retained.sampler,input.csgSampler)
        ||retained.excludedPhases!=input.csgExcludedPhases
        ||!retained.environment.same(environment)||retained.error!=error){
        work.cache->retained.reset();return std::nullopt;
    }
    try{
        std::vector<int> ids=retained.ids;
        if(work.report)work.report->peakRetainedBytes=std::max(work.report->peakRetainedBytes,csgBoundaryRetainedBytes(retained));
        return ids;
    }catch(const std::bad_alloc&){}catch(const std::length_error&){}
    return std::nullopt;
}
void retainCsgBoundaryIds(const CsgBoundaryWork& work,const NativeLeafDomain& domain,
    const HydrodynamicsInput& input,const std::vector<int>& ids,
    const MaxwellGeometryFpState& environment,int error) noexcept {
    PassiveAmrTimingEnvironment current;
    if(!current.eligible||!work.cache)return;
    work.cache->retained.reset();
    const auto limit=(csgBoundaryByteLimit-sizeof(CsgBoundaryCache)-256u)/sizeof(int);
    if(ids.size()>limit||input.csgExcludedPhases.size()>limit-ids.size())return;
    try{
        RetainedCsgBoundaryIds retained;
        retained.domain=domain.meshIdentity();retained.definition=input.csgGeometry;retained.sampler=input.csgSampler;
        retained.excludedPhases=input.csgExcludedPhases;retained.ids=ids;
        retained.environment=environment;retained.error=error;
        const auto bytes=csgBoundaryRetainedBytes(retained);
        if(bytes>csgBoundaryByteLimit)return;
        work.cache->retained.emplace(std::move(retained));
        if(work.report){++work.report->retainedBuilds;
            work.report->peakRetainedBytes=std::max(work.report->peakRetainedBytes,bytes);}
    }catch(const std::bad_alloc&){}catch(const std::length_error&){}
}
void registerAcceptedCsgBoundary(std::shared_ptr<CsgBoundaryCache>& accepted,
    std::shared_ptr<CsgBoundaryCache>& attempt,const NativeLeafDomain& actual) noexcept {
    PassiveAmrTimingEnvironment current;
    accepted.reset();
    if(attempt&&attempt->retained&&sameCsgOwner(attempt->retained->domain,actual.meshIdentity()))accepted=std::move(attempt);
    attempt.reset();
}
#if defined(_MSC_VER)
#pragma float_control(precise,on,push)
#pragma fenv_access(on)
#pragma float_control(except,on)
#endif
struct CsgBoundaryTiming {
    VofEhdCsgBoundaryStageReport* report;
    std::optional<AmrTimingClock::time_point> total,rebuild;
    bool rebuildRequested=false;
    explicit CsgBoundaryTiming(VofEhdCsgBoundaryStageReport* value) noexcept
        :report(value),total(value?amrTimingStart():std::nullopt){}
    void rebuilding() noexcept {if(report){rebuildRequested=true;rebuild=amrTimingStart();}}
    ~CsgBoundaryTiming() noexcept {
        if(!report)return;
        const auto rebuilt=rebuild?amrTimingFinish(rebuild):std::optional<double>{};
        const auto elapsed=amrTimingFinish(total);
        PassiveAmrTimingEnvironment current;
        if(!elapsed||(rebuildRequested&&!rebuilt))report->timingAvailable=false;
        if(elapsed)report->totalSeconds+=*elapsed;
        if(rebuilt)report->rebuildSeconds+=*rebuilt;
    }
};
#if defined(_MSC_VER)
#if !defined(_M_FP_STRICT) || !_M_FP_STRICT
#pragma fenv_access(off)
#endif
#pragma float_control(pop)
#endif
// One immutable Core CSG selector certifies every actual native physical boundary.
// Shared by electric conditions, pressure ports and stationary no-slip walls.
std::vector<int> csgBoundaryIds(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const HydrodynamicsInput& input,CsgBoundaryWork work)
{
    if(work.report)++work.report->requests;
    CsgBoundaryTiming timing(work.report);
    require(input.csgGeometry&&input.csgSampler,"actual CSG binding requires immutable geometry and sampler");
    domain.validateMeshIdentity(mesh);
    MaxwellGeometryFpState before;int error=0;
    const bool qualified=csgBoundaryEnvironment(before,error);
    if(qualified)if(auto cached=cachedCsgBoundaryIds(work,domain,input,before,error)){
        if(work.report)++work.report->hits;return std::move(*cached);
    }
    if(work.report){++work.report->rebuilds;++work.report->selectorCalls;}
    timing.rebuilding();
    auto ids=[&]{
    const auto solids=FittedDomainGeometry::selectFixedSolids(mesh,*input.csgGeometry,
        *input.csgSampler,input.csgExcludedPhases);
    require(solids.excludedPhaseByCell.size()==domain.cellRowByIndex().size(),"different CSG/native cell map size");
    for(std::size_t ci=0;ci<solids.excludedPhaseByCell.size();++ci)
        require((solids.excludedPhaseByCell[ci]==0)==(domain.cellRowByIndex()[ci]>=0),"actual VOF selected map differs from fixed CSG solid mask");
    FaceSegmentGraph complete;if(work.report)++work.report->completeGraphBuilds;complete.build(mesh);
    require(solids.boundaryIdsBySegment.size()==complete.segments().size(),"different actual Core CSG boundary segment support");
    std::vector<int> ids(domain.nativeFaceCount());
    for(std::size_t i=0;i<complete.segments().size();++i){const int id=solids.boundaryIdsBySegment[i];if(id==0)continue;
        const auto& segment=complete.segments()[i];const auto& edge=segment.isHorizontal?mesh.hEdges(segment.edgeLabel):mesh.vEdges(segment.edgeLabel);
        const int row=edge.index()+(segment.isHorizontal?mesh.vEdgeIndexNum():0);
        require(row>=0&&row<domain.nativeFaceCount()&&(ids[row]==0||ids[row]==id),"ambiguous actual native Core boundary binding");ids[row]=id;
    }
    return ids;
    }(); // All original classification temporaries are gone before the quiet receipt.
    MaxwellGeometryFpState after;int finalError=0;
    if(qualified&&csgBoundaryEnvironment(after,finalError)&&before.same(after)&&error==finalError)
        retainCsgBoundaryIds(work,domain,input,ids,before,error);
    return ids;
}
std::vector<E::ElectricVofBoundary> boundary(const QuadtreeMesh& mesh,const Vof& geometry,
    const VofElectrohydrodynamicsInput& input,CsgBoundaryWork work={})
{
    const auto& domain=geometry.domain();std::vector<E::ElectricVofBoundary> result(domain.nativeFaceCount());
    if(input.fluid.csgGeometry){
        require(input.fluid.csgSampler&&!input.csgElectricBoundaries.empty(),"CSG electric domain requires explicit boundary-ID conditions and immutable source sampler");
        const auto ids=csgBoundaryIds(mesh,domain,input.fluid,work);
        std::vector<unsigned> used(input.csgElectricBoundaries.size());
        for(int row=0;row<domain.nativeFaceCount();++row){const auto face=domain.face(mesh,row);
            if((face.negativeCell<0&&face.positiveCell<0)||(face.negativeCell>=0&&face.positiveCell>=0))continue;
            require(ids[row]>0,"actual exposed CSG fluid face has no declared Core boundary ID");
            const auto found=std::find_if(input.csgElectricBoundaries.begin(),input.csgElectricBoundaries.end(),
                [&](const auto& patch){return patch.boundaryId==ids[row];});
            require(found!=input.csgElectricBoundaries.end(),"actual exposed CSG boundary ID has no electric condition");
            result[row]=found->condition;++used[static_cast<std::size_t>(found-input.csgElectricBoundaries.begin())];
            if(face.area==0)require(domain.axisymmetric()&&!face.horizontal&&face.center.x==0&&
                result[row].type==E::ElectricVofBoundary::Type::HomogeneousNeumann&&result[row].value==0,
                "true zero-area RZ axis requires homogeneous Neumann0");
        }
        for(unsigned count:used)require(count>0,"electric boundary-ID condition has no actual selected-fluid coverage");
        return result;
    }
    require(input.csgElectricBoundaries.empty(),"boundary-ID conditions require an actual selected CSG domain");
    for(int row=0;row<domain.nativeFaceCount();++row){const auto face=domain.face(mesh,row);
        if(face.negativeCell<0&&face.positiveCell<0)continue;
        if(face.negativeCell>=0&&face.positiveCell>=0)continue;
        const int side=face.horizontal?(face.negativeCell<0?1:3):(face.negativeCell<0?0:2);
        result[row]=input.boundaries[side];}
    return result;
}
E::ElectricVofSharpProblem sharpProblem(const VofElectrohydrodynamicsInput& input)
{
    E::ElectricVofSharpProblem result;
    result.epsilon={{input.materials.epsilonLiquid,input.materials.epsilonGas}};
    result.conductivity={{input.materials.conductivityLiquid,input.materials.conductivityGas}};
    for(std::size_t side=0;side<result.boundaries.size();++side){const auto source=input.boundaries[side];auto& target=result.boundaries[side];
        require(source.type!=E::ElectricVofBoundary::Type::Interior,"sharp physical side has no explicit boundary condition");
        target.type=source.type==E::ElectricVofBoundary::Type::Dirichlet?E::ElectricVofSharpBoundary::Type::Dirichlet:E::ElectricVofSharpBoundary::Type::HomogeneousNeumann;
        target.value=[value=source.value](Vector2D){return value;};}
    return result;
}
E::ElectricVofSharpControls sharpControls(const VofElectrohydrodynamicsInput& input)
{return {input.controls.relativeLinearTolerance,input.controls.relativeConservationTolerance};}
void qualifySharpScope(const VofElectrohydrodynamicsInput& input,const LAMP::Interface::VofSmoothGeometry& smooth)
{
    require(!input.fluid.radial&&!input.fluid.vofRegrid&&!input.fluid.initialInterfaceBand
        &&!input.fluid.vofSouthEven90Interval&&!input.fluid.vofPressureOutletBackflow&&!input.smoothGeometry.allowBoundaryEndpoints
        &&input.fluid.vofBoundaryPatches.empty()&&input.fluid.vofSymmetryPlanes.empty()&&!input.fluid.vof.southEven90Boundary,
        "sharp coupled entry requires closed XY geometry without regrid or moving boundary endpoints");
    require(smooth.source().domain().uniformLeaves()&&!smooth.source().domain().axisymmetric(),"sharp coupled entry requires actual uniform XY leaves");
    require(!smooth.arcs().empty()&&!smooth.nodes().empty(),"sharp coupled entry requires one resolved closed component");
    std::vector<std::vector<std::size_t>> adjacent(smooth.nodes().size());
    for(const auto& arc:smooth.arcs()){adjacent.at(arc.nodes[0]).push_back(arc.nodes[1]);adjacent.at(arc.nodes[1]).push_back(arc.nodes[0]);}
    for(const auto& neighbors:adjacent)require(neighbors.size()==2,"sharp coupled entry rejects boundary endpoints and branches");
    std::vector<bool> visited(adjacent.size());std::vector<std::size_t> pending{0};visited[0]=true;std::size_t count=0;
    while(!pending.empty()){const auto node=pending.back();pending.pop_back();++count;
        for(auto other:adjacent[node])if(!visited[other]){visited[other]=true;pending.push_back(other);}}
    require(count==adjacent.size(),"sharp coupled entry supports one closed surface component only");
}
E::ElectricSurfaceChargeTransport moveSurface(const QuadtreeMesh& mesh,const E::ElectricVofSurfaceCharge& accepted,
    const LAMP::Interface::GeometricVofTransport& transport,const LAMP::Interface::VofSmoothGeometry& target,
    const E::ElectricSurfaceChargeTransferConfig& controls)
{
    try{return E::ElectricSurfaceChargeTransport::prepareStage(mesh,accepted,transport,target,controls);}
    catch(const E::ElectricSurfaceChargeTransferError& error){const auto reason=error.reason();
        if(reason==E::ElectricSurfaceChargeTransferFailure::ExcessiveDisplacement||reason==E::ElectricSurfaceChargeTransferFailure::ProjectionDistance)
            throw H::HydroTwoPhaseStepRetry(H::HydroTwoPhaseRetryReason::Transport,error.what());
        throw;}
}
struct JointSharpAccepted {
    H::HydroTwoPhaseState fluid;E::ElectricVofSurfaceCharge surface;E::ElectricVofSharpField field;
    JointSharpAccepted(H::HydroTwoPhaseState h,E::ElectricVofSurfaceCharge q,E::ElectricVofSharpField e)
        :fluid(std::move(h)),surface(std::move(q)),field(std::move(e)){}
    void validate(const QuadtreeMesh& mesh) const {
        fluid.validateFor(mesh);require(!fluid.vofPrescribedFlow(),"sharp accepted state cannot carry actual prescribedFlow capability");
        const auto clock=fluid.timeLayer();const auto& publication=fluid.phaseGeometry().vofPublication();
        const auto& smooth=surface.geometry();smooth.validateFor(mesh,publication);surface.validateFor(mesh,smooth,clock);
        require((clock.position==EhdClock::Position::Initial||clock.position==EhdClock::Position::End)
            &&sameClock(smooth.sourceTime(),clock)&&smooth.sourceIdentity()==publication.identity()
            &&field.geometry().sameInstance(smooth)&&sameClock(field.candidateTime(),clock)
            &&field.chargeCandidate().sameInstance(surface),"accepted sharp fluid/surface/field lost actual owner or full Initial/End clock");
        if(clock.position==EhdClock::Position::Initial)
            require(field.effectiveDt()==0&&field.sourceCharge().sameInstance(surface)&&field.advection()==nullptr,"sharp Initial field requires the real dt0 accepted snapshot");
        else {
            require(field.advection()!=nullptr&&field.effectiveDt()==clock.time()-clock.begin
                &&field.sourceCharge().timeLayer().time()==clock.begin,"sharp End field lacks actual full-stage Ohmic provenance");
            field.implicitInputCharge().validateFor(mesh,smooth,clock);
        }
    }
};
using Patch=HydrodynamicsInput::VofBoundaryPatch;
using Mode=LAMP::Interface::GeometricVofAdvectorSlab::Boundary::Mode;
std::vector<int> physicalPatches(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const HydrodynamicsInput& input,std::vector<bool>* stationaryRows=nullptr,CsgBoundaryWork work={})
{
    domain.validateMeshIdentity(mesh);std::vector<int> result(domain.nativeFaceCount(),-1);
    if(stationaryRows)stationaryRows->assign(domain.nativeFaceCount(),false);
    std::vector<unsigned> covered(input.vofBoundaryPatches.size());
    if(input.csgGeometry){
        require(input.vofPressureDrivenReservoirs,"CSG physical ports require actual pressure-driven reservoirs");
        const auto ids=csgBoundaryIds(mesh,domain,input,work);
        std::vector<unsigned> walls(input.vofStationaryNoSlipBoundaryIds.size());
        for(int row=0;row<domain.nativeFaceCount();++row){const auto face=domain.face(mesh,row);
            if((face.negativeCell<0&&face.positiveCell<0)||(face.negativeCell>=0&&face.positiveCell>=0))continue;
            require(ids[row]>0,"actual exposed CSG fluid face has no Core boundary ID");
            unsigned matches=0;
            for(std::size_t pi=0;pi<input.vofBoundaryPatches.size();++pi){const auto& patch=input.vofBoundaryPatches[pi];
                require(patch.boundaryId.has_value(),"actual CSG port cannot use exterior side/interval inference");
                if(*patch.boundaryId!=ids[row])continue;
                require(face.area>0,"CSG pressure reservoir selects the zero-area RZ axis");
                result[row]=static_cast<int>(pi);++covered[pi];++matches;
            }
            for(std::size_t wi=0;wi<input.vofStationaryNoSlipBoundaryIds.size();++wi){
                if(input.vofStationaryNoSlipBoundaryIds[wi]!=ids[row])continue;
                require(face.area>0,"stationary no-slip selects the zero-area RZ axis");
                if(stationaryRows)(*stationaryRows)[row]=true;++walls[wi];++matches;
            }
            if(face.area==0)require(matches==0&&domain.axisymmetric()&&!face.horizontal&&face.center.x==0,
                "only the genuine RZ axis may remain zero-area Closed");
            else require(matches==1,"every actual positive-area exposed CSG face needs one disjoint pressure-port or stationary-no-slip declaration");
        }
        for(unsigned count:covered)require(count>0,"CSG physical pressure patch has no actual selected-fluid coverage");
        for(unsigned count:walls)require(count>0,"CSG stationary no-slip ID has no actual positive-area selected-fluid coverage");
        return result;
    }
    require(input.vofStationaryNoSlipBoundaryIds.empty(),"stationary no-slip physical IDs require actual CSG geometry");
    for(int row=0;row<domain.nativeFaceCount();++row){const auto face=domain.face(mesh,row);
        if((face.negativeCell<0&&face.positiveCell<0)||(face.negativeCell>=0&&face.positiveCell>=0))continue;
        const int side=face.horizontal?(face.negativeCell<0?1:3):(face.negativeCell<0?0:2);
        const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
        const auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
        const double lower=face.horizontal?std::min(a.x,b.x):std::min(a.y,b.y);
        const double upper=face.horizontal?std::max(a.x,b.x):std::max(a.y,b.y);
        for(std::size_t pi=0;pi<input.vofBoundaryPatches.size();++pi){const auto& patch=input.vofBoundaryPatches[pi];
            require(!patch.boundaryId,"physical boundary ID requires actual CSG geometry");
            if(static_cast<int>(patch.side)!=side||!(std::max(lower,patch.interval.x)<std::min(upper,patch.interval.y)))continue;
            require(lower>=patch.interval.x&&upper<=patch.interval.y,"physical patch cuts an actual native exterior face");
            require(result[row]<0,"overlapping physical boundary patches on an actual native face");
            require(face.area>0,"physical patch selects the zero-area RZ axis");
            result[row]=static_cast<int>(pi);++covered[pi];
        }
    }
    for(const auto count:covered)require(count>0,"physical patch has no actual native exterior coverage");
    return result;
}
H::HydroVofPrescribedFlow prescribedFlow(const QuadtreeMesh& mesh,const H::HydroDofLayout& layout,
    const VofElectrohydrodynamicsInput& input,const std::vector<double>* fixedMeans=nullptr,CsgBoundaryWork work={})
{
    const auto& domain=layout.domain();std::vector<bool> stationaryRows;
    const auto patches=physicalPatches(mesh,domain,input.fluid,&stationaryRows,work);
    require(!input.fluid.vofBoundaryPatches.empty(),"prescribed flow has no physical patches");
    if(fixedMeans)require(fixedMeans->size()==patches.size(),"different source normal-mean native rows");
    H::HydroVofPrescribedFlow flow;flow.rulesByNativeRow.resize(patches.size(),Mode::Interior);
    flow.pressureDrivenReservoirs=input.fluid.vofPressureDrivenReservoirs;
    flow.normalVelocity.resize(layout.velocityCount());
    std::vector<Vector2D> vectors(flow.pressureDrivenReservoirs?0:patches.size());
    for(int row=0;row<domain.nativeFaceCount();++row){const auto face=domain.face(mesh,row);
        if(face.negativeCell>=0&&face.positiveCell>=0)continue;
        flow.rulesByNativeRow[row]=Mode::Closed;
        if(face.negativeCell<0&&face.positiveCell<0)continue;
        if(stationaryRows[row]){const int u=layout.velocityDof(row);
            require(face.area>0&&u>=0&&layout.nativeFace(u)==row,"actual stationary no-slip wall lacks its positive-area native Hydro DOF");
            flow.stationaryNoSlipVelocityDofs.push_back(u);
        }
        const int pi=patches[row];if(pi<0)continue;
        const auto& patch=input.fluid.vofBoundaryPatches[pi];
        require(patch.reservoir==Mode::PureLiquidReservoir||patch.reservoir==Mode::PureGasReservoir,"physical port requires a pure reservoir rule");
        const int u=layout.velocityDof(row);require(u>=0&&layout.nativeFace(u)==row,"physical native row has no inverse Hydro velocity DOF");
        flow.rulesByNativeRow[row]=patch.reservoir;
        if(patch.kind==Patch::Kind::PressureOutlet){flow.pressureOutlets.push_back({u,patch.pressure});continue;}
        require(!flow.pressureDrivenReservoirs,"pressure-driven reservoir cannot prescribe a fixed velocity feed");
        const double normal=fixedMeans?fixedMeans->at(row):(face.horizontal?patch.velocity.y:patch.velocity.x);
        require(std::isfinite(normal)&&normal!=0,"nonfinite/zero declared physical normal");
        flow.prescribedVelocityDofs.push_back(u);flow.normalVelocity[u]=normal;vectors[row]=patch.velocity;
        (face.horizontal?vectors[row].y:vectors[row].x)=normal;
    }
    std::sort(flow.stationaryNoSlipVelocityDofs.begin(),flow.stationaryNoSlipVelocityDofs.end());
    require(std::adjacent_find(flow.stationaryNoSlipVelocityDofs.begin(),flow.stationaryNoSlipVelocityDofs.end())==flow.stationaryNoSlipVelocityDofs.end(),
        "duplicate actual stationary no-slip native Hydro DOF");
    if(input.fluid.vofPressureOutletBackflow){const auto& declaration=*input.fluid.vofPressureOutletBackflow;
        using Model=H::HydroVofPressureOutletMomentumModel;
        require(!flow.pressureOutlets.empty(),"explicit backflow has no actual native pressure outlets");
        require(declaration.model==Model::Ambient||declaration.model==Model::InteriorContinuationZeroShear,
            "unknown declared pressure-outlet momentum model");
        flow.pressureOutletBackflowModel=declaration.model;
        if(declaration.model==Model::Ambient){
            require(declaration.velocity.has_value(),"ambient backflow requires its actual complete vector");
            const auto velocity=*declaration.velocity;
            require(std::isfinite(velocity.x)&&std::isfinite(velocity.y),"nonfinite declared ambient momentum vector");
            std::vector<bool> allowed(domain.nativeFaceCount());
            for(const auto& port:flow.pressureOutlets)allowed[layout.nativeFace(port.velocityDof)]=true;
            std::vector<Vector2D> ambient(domain.nativeFaceCount());
            for(const auto& port:flow.pressureOutlets)ambient[layout.nativeFace(port.velocityDof)]=velocity;
            flow.pressureOutletBackflowMomentum=H::HydroVofSteadyVelocity::prepare(mesh,layout,
                std::move(ambient),std::move(allowed),H::HydroVofSteadyVelocity::ClockRole::FullMidpoint);
        }else require(!declaration.velocity.has_value(),"interior continuation forbids an unused ambient vector");
    }
    // An all-pressure reservoir has no prescribed complete endpoint vector.
    // The module derives its free normal flow from the actual physical pressures.
    if(flow.pressureDrivenReservoirs){
        require(input.fluid.vofPressureOutletBackflow.has_value(),"pressure-driven reservoirs require an explicit backflow model");
        require(flow.prescribedVelocityDofs.empty()&&!flow.endpointVelocity,"pressure-driven reservoirs cannot carry prescribed endpoint velocity");
        return flow;
    }
    std::vector<bool> covered(domain.nativeFaceCount());
    for(const int u:flow.prescribedVelocityDofs)covered[layout.nativeFace(u)]=true;
    flow.endpointVelocity=H::HydroVofSteadyVelocity::prepare(mesh,layout,std::move(vectors),
        std::move(covered),H::HydroVofSteadyVelocity::ClockRole::Endpoint);
    return flow;
}
H::HydroVofSouthEven90Boundary southEven90Geometry(const QuadtreeMesh& mesh,const Vof& publication,Vector2D interval)
{
    NativeMeshIdentity::ValidationBatch validationBatch(mesh);
    publication.validateFor(mesh,publication.identity(),publication.timeLayer());
    require(publication.timeLayer().position==EhdClock::Position::Midpoint,"S Even90 provider requires the actual stage Midpoint");
    const auto& domain=publication.domain();std::vector<int> rows;
    for(int row=0;row<domain.nativeFaceCount();++row){const auto face=domain.face(mesh,row);
        if(!face.horizontal||face.negativeCell>=0||face.positiveCell<0||face.center.y!=mesh.minCoord.y)continue;
        const auto& edge=mesh.hEdges(face.edgeLabel);const auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
        const double lower=std::min(a.x,b.x),upper=std::max(a.x,b.x);
        if(!(std::max(lower,interval.x)<std::min(upper,interval.y)))continue;
        require(lower>=interval.x&&upper<=interval.y,"S Even90 physical declaration cuts an actual native face");
        require(face.area>0,"S Even90 physical face has zero measure");rows.push_back(row);
    }
    const auto boundary=H::HydroVofSouthEven90Boundary::prepare(mesh,publication,rows);
    const auto& covered=boundary.intervals();
    require(covered.size()==1&&covered.front()[0]==interval.x&&covered.front()[1]==interval.y,
        "S Even90 actual native faces do not completely cover the declared physical interval");
    return boundary;
}
E::ElectricVofAdvectedCharge advectMidpointCharge(const QuadtreeMesh& mesh,const E::ElectricVofCharge& source,
    const LAMP::Interface::GeometricVofTransport& transport,const VofElectrohydrodynamicsInput& input,CsgBoundaryWork work)
{
    // Fixed physical declarations and the owned immutable Core CSG sampler
    // contain no endpoint/provider callback in this charge preparation.
    NativeMeshIdentity::ValidationBatch validationBatch(mesh);
    if(input.fluid.vofBoundaryPatches.empty())return E::ElectricVofChargeTransport::prepareMidpointClosed(mesh,source,transport,input.controls);
    const auto& domain=transport.source().domain();const auto patches=physicalPatches(mesh,domain,input.fluid,nullptr,work);
    require(input.inletChargeDensityByPatch.size()==input.fluid.vofBoundaryPatches.size(),"different physical inlet-charge patch table");
    const auto* fluidBoundary=transport.advector().boundary();
    require(fluidBoundary!=nullptr,"physical charge patches have no actual fluid boundary owner");
    fluidBoundary->validateFor(mesh,transport.source());
    std::vector<double> rho(domain.nativeFaceCount());
    for(const auto& flux:transport.faces()){const auto face=domain.face(mesh,flux.nativeRow);
        if(face.outwardSign==0)continue;
        const auto incomingGross=face.outwardSign>0?flux.full.negativeTotal:flux.full.positiveTotal;
        if(!(incomingGross>0))continue;
        const int pi=patches.at(flux.nativeRow);require(pi>=0,"actual boundary inflow has no declared physical patch");
        const auto& patch=input.fluid.vofBoundaryPatches[pi];
        require((patch.reservoir==Mode::PureLiquidReservoir||patch.reservoir==Mode::PureGasReservoir)&&
            fluidBoundary->rules().at(flux.nativeRow)==patch.reservoir,"actual incoming phase reservoir differs from its physical patch");
        const bool steadyInlet=patch.kind==Patch::Kind::PrescribedVelocity&&
            face.outwardSign*(face.horizontal?patch.velocity.y:patch.velocity.x)<0;
        const bool pressureBackflow=patch.kind==Patch::Kind::PressureOutlet&&input.fluid.vofPressureOutletBackflow.has_value();
        require(steadyInlet||pressureBackflow,
            "actual boundary inflow is not a declared steady inlet or explicit pressure-backflow charge reservoir");
        rho[flux.nativeRow]=input.inletChargeDensityByPatch[pi];
    }
    const auto chargeBoundary=E::ElectricVofChargeBoundary::prepare(mesh,source,transport,rho);
    return E::ElectricVofChargeTransport::prepareMidpointBoundary(mesh,source,transport,chargeBoundary,input.controls);
}
struct JointAccepted {
    H::HydroTwoPhaseState fluid;
    E::ElectricVofCharge charge;
    E::ElectricVofField field;
    JointAccepted(H::HydroTwoPhaseState h,E::ElectricVofCharge q,E::ElectricVofField e)
        :fluid(std::move(h)),charge(std::move(q)),field(std::move(e)){}
    void validate(const QuadtreeMesh& mesh) const {
        // State validation can invoke a retained endpoint law. Only the
        // following immutable charge/field checks share a local receipt.
        {NativeMeshIdentity::ValidationBatch::Pause pause;fluid.validateFor(mesh);}
        NativeMeshIdentity::ValidationBatch validationBatch(mesh);
        const auto clock=fluid.timeLayer();
        const auto& geometry=fluid.phaseGeometry().vofPublication();charge.validateFor(mesh,clock);
        field.validateFor(mesh,geometry,field.sourceCharge(),clock);
        require(sameClock(field.geometryTime(),clock)&&sameClock(field.candidateTime(),clock)&&
            field.geometry().identity()==geometry.identity()&&field.chargeCandidate().sameInstance(charge),
            "accepted fluid/charge/field must share actual geometry and complete Initial/End clock");
    }
};
}
struct VofElectrohydrodynamicsTransaction::Data {
    const VofElectrohydrodynamicsInput input;
    VofElectrohydrodynamicsOutput& writer;
    std::shared_ptr<JointAccepted> accepted,prepared;
    std::shared_ptr<JointSharpAccepted> acceptedSharp,preparedSharp;
    std::optional<H::HydroVofStep> step;
    std::optional<E::ElectricVofChargeTransport> physicalTransport;
    std::optional<E::ElectricVofField> physicalEnd;
    std::optional<HydrodynamicsVofRegridDecision> regrid;
    std::optional<H::HydroVofMeshTransfer> fluidRemap;
    std::optional<E::ElectricVofChargeMeshTransfer> chargeRemap;
    HydrodynamicsVofRemapReceipt receipt;
    std::optional<RetainedMaxwellGeometry> acceptedMaxwellGeometry,pendingMaxwellGeometry;
    std::optional<RetainedMaxwellGeometry> targetMaxwellGeometry; // fresh target pure geometry only
    std::shared_ptr<CsgBoundaryCache> acceptedCsgBoundary,attemptCsgBoundary;
    VofEhdCsgBoundaryReport csgBoundaryReport;
    bool committed=false;
    Data(const VofElectrohydrodynamicsInput& in,VofElectrohydrodynamicsOutput& out):input(in),writer(out){}
};
VofElectrohydrodynamicsTransaction::VofElectrohydrodynamicsTransaction(
    const VofElectrohydrodynamicsInput& input,VofElectrohydrodynamicsOutput& writer)
    :data_(std::make_unique<Data>(input,writer)){}
VofElectrohydrodynamicsTransaction::~VofElectrohydrodynamicsTransaction(){clear();}
const char* VofElectrohydrodynamicsTransaction::name() const{return "VofElectrohydrodynamics";}
const char* VofElectrohydrodynamicsTransaction::bannerSuffix() const{return data_->input.method==VofElectrohydrodynamicsInput::Method::VofSharpSurface?" / B1 smooth surface charge / sharp coeval Maxwell":" / conservative Eulerian charge / coeval Maxwell";}
const char* VofElectrohydrodynamicsTransaction::inputOkMessage() const{return "VOF_EHD_INPUT_OK";}
const char* VofElectrohydrodynamicsTransaction::completionSuffix() const{return " joint fluid/charge/field accepted";}
const char* VofElectrohydrodynamicsTransaction::prepareStage() const{return "VOF-EHD AC preparation";}
H::HydroTwoPhaseState VofElectrohydrodynamicsTransaction::makeInitialVofState(const QuadtreeMesh& mesh,
    H::HydroPhaseGeometry phase,const H::HydroNativeValues& seed) const
{
    if(data_->input.method==VofElectrohydrodynamicsInput::Method::VofSharpSurface){
        require(!data_->input.fluid.radial&&!data_->input.fluid.vofRegrid&&!data_->input.fluid.initialInterfaceBand
            &&data_->input.fluid.vofBoundaryPatches.empty()&&data_->input.fluid.vofSymmetryPlanes.empty()
            &&!data_->input.fluid.vofSouthEven90Interval&&!data_->input.fluid.vof.southEven90Boundary
            &&!data_->input.fluid.vofPressureOutletBackflow,"sharp Initial hook requires its actual closed uniform XY capability");
        return HydrodynamicsTransaction::makeInitialVofState(mesh,std::move(phase),seed);
    }
    if(data_->input.fluid.vofBoundaryPatches.empty()){
        require(!data_->input.fluid.vofPressureOutletBackflow,"closed initialization cannot consume pressure-outlet backflow data");
        return HydrodynamicsTransaction::makeInitialVofState(mesh,std::move(phase),seed);
    }
    require(phase.timeLayer().position==EhdClock::Position::Initial,"open initialization requires the actual Initial layer");
    NativeMeshIdentity::ValidationBatch validationBatch(mesh);
    const auto& layout=phase.layout();seed.validateFor(layout,mesh);
    auto work=csgBoundaryWork(data_->acceptedCsgBoundary,data_->csgBoundaryReport.stages[0],bool(data_->input.fluid.csgGeometry));
    auto flow=prescribedFlow(mesh,layout,data_->input,nullptr,work);
    const auto geometry=H::HydroMomentumGeometry::build(mesh,phase.metrics());const auto& config=data_->input.fluid.vof;
    const auto mass=config.momentumDensityRepresentation==H::HydroMomentumMass::DensityRepresentation::CellAverageProjection
        ?H::HydroMomentumMass::buildCellAverageProjection(mesh,geometry,phase,config.liquidDensity,config.gasDensity)
        :H::HydroMomentumMass::build(mesh,geometry,phase,config.liquidDensity,config.gasDensity);
    H::HydroSaddleVelocityBlock::Matrix diagonal(layout.velocityCount(),layout.velocityCount());
    std::vector<double> rhs(layout.velocityCount()),fixed(layout.velocityCount());
    for(int u=0;u<layout.velocityCount();++u){diagonal.insert(u,u)=mass.masses()[u];
        const auto value=static_cast<long double>(mass.masses()[u])*seed.normalVelocity()[u];rhs[u]=static_cast<double>(value);
        require(std::isfinite(rhs[u])&&(value==0||rhs[u]!=0),"unrepresentable original Initial M-projection RHS");}
    diagonal.makeCompressed();for(const int u:flow.prescribedVelocityDofs)fixed[u]=flow.normalVelocity[u];
    std::vector<H::HydroSaddlePressureOutlet> projectionOutlets;
    for(const auto& outlet:flow.pressureOutlets)projectionOutlets.push_back({outlet.velocityDof,0});
    const auto pressure=H::HydroPressureLoad::build(phase.metrics().divergence());
    const auto block=H::HydroSaddleVelocityBlock::build(layout,diagonal);
    H::HydroSaddleConfig projectionConfig;projectionConfig.gmresOnlyParallel=true;projectionConfig.gmresGlobalPressure=true;
    const auto saddle=projectionOutlets.empty()?H::HydroSinglePhaseSaddle::buildClosed(pressure,block,projectionConfig)
        :H::HydroSinglePhaseSaddle::buildWithPressureOutlets(pressure,block,projectionOutlets,projectionConfig);
    const auto projected=saddle.solve(mesh,H::HydroSaddleRhs(layout,rhs,fixed));
    // Constraint lambda is not physical initial p or a fabricated time step.
    // The state factory may query its newly prepared endpoint law.
    NativeMeshIdentity::ValidationBatch::Pause pause;
    return H::HydroTwoPhaseState::createVofPrescribedFlow(mesh,std::move(phase),
        H::HydroNativeValues(layout,projected.values().normalVelocity(),seed.pressure()),std::move(flow));
}
const H::HydroTwoPhaseState* VofElectrohydrodynamicsTransaction::accepted() const
{return data_->input.method==VofElectrohydrodynamicsInput::Method::VofSharpSurface?(data_->acceptedSharp?&data_->acceptedSharp->fluid:nullptr):(data_->accepted?&data_->accepted->fluid:nullptr);}
bool VofElectrohydrodynamicsTransaction::supportsCheckpoint() const noexcept
{return data_->input.method==VofElectrohydrodynamicsInput::Method::VofFiniteVolume;}
std::string VofElectrohydrodynamicsTransaction::checkpoint(const QuadtreeMesh& mesh) const
{
    require(supportsCheckpoint()&&data_->accepted,"checkpoint requires an accepted Eulerian-charge state");
    const auto& a=*data_->accepted;a.validate(mesh);
    LAMP::IO::Checkpoint::Writer out("LAMP/VofElectrohydrodynamics/transaction/1");
    out.text(a.fluid.checkpoint(mesh,data_->step?&data_->step->beginState():nullptr));
    out.values(a.charge.integratedCharge());
    const auto* flow=a.fluid.vofPrescribedFlow();out.boolean(flow!=nullptr);
    if(flow){
        const auto& layout=a.fluid.values().layout();
        std::vector<double> means(layout.domain().nativeFaceCount());
        for(const int u:flow->prescribedVelocityDofs)means.at(layout.nativeFace(u))=flow->normalVelocity.at(u);
        out.values(means);
    }
    out.text(data_->writer.checkpoint());return out.finish();
}
void VofElectrohydrodynamicsTransaction::restoreCheckpoint(const QuadtreeMesh& mesh,const std::string& bytes)
{
    require(supportsCheckpoint()&&!data_->accepted&&!data_->acceptedSharp,"checkpoint restore requires a fresh Eulerian-charge transaction");
    LAMP::IO::Checkpoint::Reader in(bytes,"LAMP/VofElectrohydrodynamics/transaction/1");
    const auto fluidBytes=in.text();const auto charge=in.values<double>();const bool hasFlow=in.boolean();
    const auto means=hasFlow?in.values<double>():std::vector<double>{};const auto observer=in.text();in.finish();
    require(hasFlow==!data_->input.fluid.vofBoundaryPatches.empty(),"checkpoint flow capability differs from physical declarations");
    std::optional<H::HydroVofPrescribedFlow> flow;
    auto work=csgBoundaryWork(data_->acceptedCsgBoundary,data_->csgBoundaryReport.stages[0],bool(data_->input.fluid.csgGeometry));
    if(!data_->input.fluid.vofBoundaryPatches.empty()){
        FaceSegmentGraph graph;graph.build(mesh);
        const auto domain=NativeLeafDomain::build(mesh,graph,data_->input.fluid.radial);
        flow=prescribedFlow(mesh,H::HydroDofLayout::build(mesh,domain),data_->input,&means,work);
    }
    auto fluid=H::HydroTwoPhaseState::restoreCheckpoint(mesh,fluidBytes,std::move(flow));
    const auto& geometry=fluid.phaseGeometry().vofPublication();
    auto q=E::ElectricVofCharge::fromIntegrated(mesh,geometry,charge);
    // The potential is an algebraic snapshot of the accepted Q. It is never
    // used as a fabricated Ohmic step or as predecessor charge history.
    auto field=E::ElectricVofField::prepareElectrostaticSnapshot(mesh,geometry,q,data_->input.materials,
        boundary(mesh,geometry,data_->input,work),data_->input.controls,E::ElectricVofField::CellGradientRecovery::Omit);
    auto restored=std::make_shared<JointAccepted>(std::move(fluid),std::move(q),std::move(field));restored->validate(mesh);
    data_->writer.restoreCheckpoint(mesh,observer,restored->fluid,restored->charge);data_->accepted.swap(restored);
}
const E::ElectricVofCharge* VofElectrohydrodynamicsTransaction::acceptedCharge() const noexcept
{return data_->accepted?&data_->accepted->charge:nullptr;}
const E::ElectricVofField* VofElectrohydrodynamicsTransaction::acceptedField() const noexcept
{return data_->accepted?&data_->accepted->field:nullptr;}
const E::ElectricVofSurfaceCharge* VofElectrohydrodynamicsTransaction::acceptedSurfaceCharge() const noexcept
{return data_->acceptedSharp?&data_->acceptedSharp->surface:nullptr;}
const E::ElectricVofSharpField* VofElectrohydrodynamicsTransaction::acceptedSharpField() const noexcept
{return data_->acceptedSharp?&data_->acceptedSharp->field:nullptr;}
const H::HydroVofStep* VofElectrohydrodynamicsTransaction::preparedStep() const noexcept
{return data_->step&&(data_->prepared||data_->preparedSharp)?&*data_->step:nullptr;}
const E::ElectricVofChargeTransport* VofElectrohydrodynamicsTransaction::physicalChargeTransport() const noexcept
{return data_->physicalTransport?&*data_->physicalTransport:nullptr;}
const E::ElectricVofField* VofElectrohydrodynamicsTransaction::physicalEndField() const noexcept
{return data_->physicalEnd?&*data_->physicalEnd:nullptr;}
const VofEhdCsgBoundaryReport* VofElectrohydrodynamicsTransaction::csgBoundaryPreparation() const noexcept
{
    for(const auto& stage:data_->csgBoundaryReport.stages)if(stage.requests)return &data_->csgBoundaryReport;
    return nullptr;
}
void VofElectrohydrodynamicsTransaction::initializeVof(const QuadtreeMesh& mesh,H::HydroTwoPhaseState fluid,
    Clock::time_point start,const char*& stage)
{
    require(!data_->accepted&&!data_->acceptedSharp,"initialization attempted over an accepted joint state");
    if(data_->input.method==VofElectrohydrodynamicsInput::Method::VofSharpSurface){
        fluid.validateFor(mesh);require(!fluid.vofPrescribedFlow(),"sharp Initial cannot accept actual prescribedFlow state");
        const auto& publication=fluid.phaseGeometry().vofPublication();
        stage="B1 actual Initial smooth geometry";
        const auto smooth=LAMP::Interface::VofSmoothGeometry::prepare(mesh,publication,data_->input.smoothGeometry);qualifySharpScope(data_->input,smooth);
        const auto q=[&]{NativeMeshIdentity::ValidationBatch::Pause pause;
            return E::ElectricVofSurfaceCharge::fromDensity(mesh,smooth,[&](std::size_t arc,double parameter){const auto point=smooth.arcs()[arc].sample(parameter);
                return data_->input.initialSurfaceChargeDensity({static_cast<double>(point.x),static_cast<double>(point.y)});});}();
        stage="B1 initial electrostatic snapshot";
        auto field=E::ElectricVofSharpField::prepareElectrostatic(mesh,smooth,q,sharpProblem(data_->input),sharpControls(data_->input));
        auto initial=std::make_shared<JointSharpAccepted>(std::move(fluid),q,std::move(field));initial->validate(mesh);
        data_->acceptedSharp.swap(initial);stage="accepted B1 output";
        data_->writer.accepted(mesh,data_->acceptedSharp->fluid,data_->acceptedSharp->surface,data_->acceptedSharp->field,nullptr,start,
            std::chrono::duration<double>(Clock::now()-start).count(),false);return;
    }
    const auto& geometry=fluid.phaseGeometry().vofPublication();const auto& domain=geometry.domain();
    auto work=csgBoundaryWork(data_->acceptedCsgBoundary,data_->csgBoundaryReport.stages[0],bool(data_->input.fluid.csgGeometry));
    std::vector<double> charge(mesh.cellIndexNum());
    if(data_->input.initialChargeType==VofElectrohydrodynamicsInput::InitialChargeType::SurfaceChargeProjection){
        fluid.validateFor(mesh);
        require(geometry.timeLayer().position==EhdClock::Position::Initial&&!fluid.vofPrescribedFlow()
            &&!data_->input.fluid.radial&&domain.uniformLeaves()&&!domain.axisymmetric()
            &&data_->input.fluid.vofBoundaryPatches.empty()&&data_->input.fluid.vofSymmetryPlanes.empty()
            &&!data_->input.fluid.vofSouthEven90Interval&&!data_->input.fluid.vof.southEven90Boundary
            &&!data_->input.fluid.vofPressureOutletBackflow&&!data_->input.smoothGeometry.allowBoundaryEndpoints,
            "surfaceChargeProjection requires actual closed uniform XY Initial geometry without contact boundary data");
        stage="A actual Initial surface-charge projection geometry";
        const auto smooth=LAMP::Interface::VofSmoothGeometry::prepare(mesh,geometry,data_->input.smoothGeometry);
        require(!smooth.arcs().empty()&&!smooth.nodes().empty()&&smooth.report().boundaryEndpoints==0,
            "surfaceChargeProjection requires a resolved closed Initial interface");
        const auto surface=[&]{NativeMeshIdentity::ValidationBatch::Pause pause;
            return E::ElectricVofSurfaceCharge::fromDensity(mesh,smooth,[&](std::size_t arc,double parameter){
                const auto point=smooth.arcs()[arc].sample(parameter);
                return data_->input.initialSurfaceChargeDensity({static_cast<double>(point.x),static_cast<double>(point.y)});});}();
        surface.validateFor(mesh,smooth,geometry.timeLayer());
        std::vector<FV::Dyadic> cellCharge(mesh.cellIndexNum());
        for(std::size_t arc=0;arc<smooth.arcs().size();++arc){const int native=smooth.arcs()[arc].nativeCell;
            require(native>=0&&native<mesh.cellIndexNum(),"surfaceChargeProjection arc outside actual native cells");
            cellCharge[native].add(surface.moments()[arc][0]);}
        for(int ci=0;ci<mesh.cellIndexNum();++ci){charge[ci]=FV::roundDouble(cellCharge[ci]).value;
            require(std::isfinite(charge[ci]),"unrepresentable projected Initial integrated surface charge");}
    }else{
        for(int ci=0;ci<mesh.cellIndexNum();++ci){
            const auto center=mesh.cells(mesh.cellIndexToLabel[ci]).centerCoord();
            const auto density=[&]{NativeMeshIdentity::ValidationBatch::Pause pause;
                return data_->input.initialChargeDensity(center);}();
            const auto volume=domain.cellVolumes()[domain.cellRowByIndex()[ci]];
            charge[ci]=FV::roundDouble(FV::Dyadic::product({density,volume})).value;
            require(std::isfinite(charge[ci]),"unrepresentable actual initial integrated charge");
        }
    }
    stage="VOF-EHD initial electrostatic snapshot";
    {
    // All initial density callbacks have finished; output remains outside.
    NativeMeshIdentity::ValidationBatch validationBatch(mesh);
    auto q=E::ElectricVofCharge::fromIntegrated(mesh,geometry,charge);
    auto field=E::ElectricVofField::prepareElectrostaticSnapshot(mesh,geometry,q,data_->input.materials,
        boundary(mesh,geometry,data_->input,work),data_->input.controls,E::ElectricVofField::CellGradientRecovery::Omit);
    auto initial=std::make_shared<JointAccepted>(std::move(fluid),std::move(q),std::move(field));initial->validate(mesh);
    data_->accepted.swap(initial);stage="accepted VOF-EHD output";
    }
    data_->writer.accepted(mesh,data_->accepted->fluid,data_->accepted->charge,data_->accepted->field,nullptr,start,
        std::chrono::duration<double>(Clock::now()-start).count(),false,nullptr,nullptr,nullptr,nullptr,csgBoundaryPreparation());
}
void VofElectrohydrodynamicsTransaction::prepare(const QuadtreeMesh& mesh,const EhdClock& target,
    H::HydroTwoPhaseStepConfig,bool,const char*& stage)
{
    releaseStep();
    if(data_->input.method==VofElectrohydrodynamicsInput::Method::VofSharpSurface){
        require(bool(data_->acceptedSharp),"preparation has no accepted sharp joint state");data_->acceptedSharp->validate(mesh);
        auto config=data_->input.fluid.vof;config.requireMaxwellCorrection=true;const auto source=data_->acceptedSharp;const auto* input=&data_->input;
        config.transportProvider=[source,input](const QuadtreeMesh& grid,const Vof& old,
            const LAMP::Interface::GeometricVofAdvectorSlab& slab,long double rhoL,long double rhoG,
            const LAMP::Interface::GeometricVofTransportConfig& controls){
            return LAMP::Interface::GeometricVofTransport::prepareLocalSmoothSweep(grid,old,
                source->surface.geometry(),slab,rhoL,rhoG,controls,input->smoothGeometry);
        };
        config.maxwellLoad=[source,input](const QuadtreeMesh& grid,const LAMP::Interface::GeometricVofTransport& transport,
            const H::HydroPhaseGeometry& mid,const H::HydroPressureLoad& pressure){
            require(transport.smoothMidpoint()!=nullptr,"B1 requires its actual smooth phase-flux Midpoint owner");
            const auto& smooth=*transport.smoothMidpoint();smooth.validateFor(grid,mid.vofPublication());qualifySharpScope(*input,smooth);
            const auto moved=moveSurface(grid,source->surface,transport,smooth,input->surfaceTransport);
            const auto field=E::ElectricVofSharpField::prepareOhmic(grid,moved,sharpProblem(*input),sharpControls(*input));
            const auto stress=E::ElectricVofSharpMaxwellStress::prepare(grid,field);
            return H::HydroMaxwellLoad::prepareSharp(grid,pressure,mid,stress);
        };
        stage="B1 actual AC sharp Maxwell preparation";data_->step=H::HydroVofStep::prepare(mesh,source->fluid,target,config);
        const auto& physical=data_->step->candidateState();stage="B1 actual End smooth geometry";
        require(data_->step->transport().smoothEndpoint()!=nullptr,"B1 requires its actual smooth phase-flux End owner");
        const auto& smooth=*data_->step->transport().smoothEndpoint();smooth.validateFor(mesh,physical.phaseGeometry().vofPublication());qualifySharpScope(*input,smooth);
        stage="B1 final transport charge and End Ohmic";
        const auto moved=moveSurface(mesh,source->surface,data_->step->transport(),smooth,input->surfaceTransport);
        auto field=E::ElectricVofSharpField::prepareOhmic(mesh,moved,sharpProblem(*input),sharpControls(*input));
        auto q=field.chargeCandidate();data_->preparedSharp=std::make_shared<JointSharpAccepted>(physical,std::move(q),std::move(field));
        data_->preparedSharp->validate(mesh);return;
    }
    // One source-mesh preparation only. Destroy this receipt before AMR;
    // callback-capable helpers suspend every enclosing receipt separately.
    std::optional<NativeMeshIdentity::ValidationBatch> sourceMeshBatch;
    sourceMeshBatch.emplace(mesh);
    // Consume before any A preparation can fail. A retry/rejected candidate
    // cannot reuse the prior accepted offer; only a successful commit returns it.
    auto offeredMaxwellGeometry=std::exchange(data_->acceptedMaxwellGeometry,std::nullopt);
    data_->attemptCsgBoundary=std::exchange(data_->acceptedCsgBoundary,{});
    auto physicalCsg=csgBoundaryWork(data_->attemptCsgBoundary,data_->csgBoundaryReport.stages[1],bool(data_->input.fluid.csgGeometry));
    data_->targetMaxwellGeometry.reset();
    require(bool(data_->accepted),"preparation has no accepted joint state");data_->accepted->validate(mesh);
    auto config=data_->input.fluid.vof;
    config.saddle.gmresOnlyParallel=true;config.saddle.gmresGlobalPressure=true;
    const auto source=data_->accepted;const auto* input=&data_->input;
    const auto* sourceFlow=source->fluid.vofPrescribedFlow();
    const bool greenOnly=(config.liquidViscosity==0&&config.gasViscosity==0)
        ||(sourceFlow&&!sourceFlow->stationaryNoSlipVelocityDofs.empty());
    if(data_->input.fluid.vofSouthEven90Interval){
        // Geometry is rebuilt at each actual stage on its actual tree. The
        // capture owns only the physical interval, never a mesh or old pub.
        const auto interval=*data_->input.fluid.vofSouthEven90Interval;
        config.southEven90Boundary=[interval](const QuadtreeMesh& grid,const Vof& actual){return southEven90Geometry(grid,actual,interval);};
    }

    // The original within-attempt reuse is unchanged. A previous accepted
    // pure-geometry offer is separate, consumed once, and never contains Q,
    // phase, stress, load, solver or a previous accepted-state handle.
    const auto maxwellGeometry=std::make_shared<std::optional<H::HydroStrainOperator>>();
    const auto retainedGeometry=std::make_shared<AttemptMaxwellGeometry>(std::move(offeredMaxwellGeometry));
    config.maxwellLoad=[source,input,maxwellGeometry,retainedGeometry,physicalCsg,greenOnly](const QuadtreeMesh& grid,const LAMP::Interface::GeometricVofTransport& transport,
        const H::HydroPhaseGeometry& mid,const H::HydroPressureLoad& pressure){
        // This App-owned provider uses fixed physical declarations and typed
        // kernels. Hydro suspends its opaque caller receipt before entry.
        NativeMeshIdentity::ValidationBatch validationBatch(grid);
        const auto moved=advectMidpointCharge(grid,source->charge,transport,*input,physicalCsg);
        const auto field=E::ElectricVofField::prepareOhmic(grid,mid.vofPublication(),moved,input->materials,
            boundary(grid,mid.vofPublication(),*input,physicalCsg),input->controls,E::ElectricVofField::CellGradientRecovery::Omit,&source->field);
        const auto stress=E::ElectricVofMaxwellStress::prepare(grid,field);
        require(source->fluid.values().layout().sameIdentity(mid.layout())
            && &source->fluid.phaseGeometry().metrics().divergence().matrix()==&mid.metrics().divergence().matrix(),
            "Maxwell stage lost the actual accepted native layout/divergence");
        if(maxwellGeometry->has_value()){
            require(maxwellGeometryMode(maxwellGeometry->value(),greenOnly),"Maxwell stage geometry has a different requested capability");
            return H::HydroMaxwellLoad::prepare(grid,pressure,mid,stress,maxwellGeometry->value());
        }
        auto offered=std::exchange(retainedGeometry->offered,std::nullopt);
        MaxwellGeometryFpState before;const bool qualifiedFp=MaxwellGeometryFpState::read(before);
        const bool reuse=qualifiedFp&&offered&&before.same(offered->environment)
            &&maxwellGeometryMode(offered->geometry,greenOnly)&&offered->geometry.layout().sameIdentity(mid.layout())
            && &offered->geometry.metrics().divergence().matrix()==&mid.metrics().divergence().matrix();
        if(!reuse)offered.reset(); // Drop an obsolete geometry BEFORE a fresh build.
        auto load=reuse?H::HydroMaxwellLoad::prepare(grid,pressure,mid,stress,offered->geometry)
            :greenOnly?H::HydroMaxwellLoad::prepareRt0Green(grid,pressure,mid,stress)
                :H::HydroMaxwellLoad::prepare(grid,pressure,mid,stress);
        MaxwellGeometryFpState after;
        if(qualifiedFp&&MaxwellGeometryFpState::read(after)&&before.same(after))
            retainedGeometry->registered.emplace(RetainedMaxwellGeometry{load.geometry(),before});
        maxwellGeometry->emplace(load.geometry());
        return load;
    };
    stage="VOF-EHD AC preparation";data_->step=source->fluid.vofPrescribedFlow()
        ?H::HydroVofStep::preparePrescribedFlow(mesh,source->fluid,target,config)
        :H::HydroVofStep::prepare(mesh,source->fluid,target,config);
    data_->pendingMaxwellGeometry=std::move(retainedGeometry->registered);
    const auto& physical=data_->step->candidateState();const auto& end=physical.phaseGeometry().vofPublication();
    stage="VOF-EHD final transport charge and End Ohmic";
    const auto& finalLoad=data_->step->frozenMidpoint().forces().maxwellLoad;
    require(finalLoad&&!finalLoad->isSharp(),"A End completion requires the actual final Frozen Maxwell load");
    const auto& finalField=finalLoad->stress().field();const auto* half=finalField.advection();const auto& finalTransport=data_->step->transport();
    const auto midClock=finalTransport.midpoint().timeLayer();
    require(half&&half->sourceCharge().sameInstance(source->charge)&&half->transport().sameInstance(finalTransport)
        &&sameClock(half->timeLayer(),midClock)&&sameClock(finalField.geometryTime(),midClock)
        &&sameClock(finalField.candidateTime(),midClock)&&sameClock(finalTransport.endTime(),target)
        &&half->geometry().identity()==finalTransport.midpoint().identity()
        &&finalField.geometry().identity()==finalTransport.midpoint().identity(),
        "A End completion requires the same accepted charge/final transport owner and complete Midpoint clock");
    data_->physicalTransport=E::ElectricVofChargeTransport::completeEndpoint(mesh,*half);
    auto field=E::ElectricVofField::prepareOhmic(mesh,end,data_->physicalTransport->endpoint(),input->materials,
        boundary(mesh,end,*input,physicalCsg),input->controls,E::ElectricVofField::CellGradientRecovery::Omit,&finalField);
    data_->physicalEnd=field;auto q=field.chargeCandidate();
    if(!input->fluid.vofRegrid||!hydrodynamicsVofRegridDue(*input->fluid.vofRegrid,target)){
        data_->prepared=std::make_shared<JointAccepted>(physical,q,std::move(field));data_->prepared->validate(mesh);return;
    }
    stage="VOF-EHD regrid tree preparation";const auto regridStart=Clock::now();
    std::optional<HydrodynamicsVofRegridRequirements> electricRequirements;
    if(input->chargeDensityJumpTolerance){
        const auto indicators=E::ElectricVofRefinement::sample(mesh,end,q);
        const auto& policy=*input->fluid.vofRegrid;
        const long double high=*input->chargeDensityJumpTolerance;
        const long double low=high*policy.coarsenRatio;
        std::vector<int> minimum(static_cast<std::size_t>(mesh.cellIndexNum()),policy.minLevel);
        for(const int ci:end.domain().cellIndexByRow()){
            const int level=mesh.cells(mesh.cellIndexToLabel[ci]).level();
            const auto jump=indicators.densityJump()[static_cast<std::size_t>(ci)];
            if(jump>high)minimum[static_cast<std::size_t>(ci)]=std::clamp(level+1,policy.minLevel,policy.maxLevel);
            else if(jump>0&&jump>=low)minimum[static_cast<std::size_t>(ci)]=std::clamp(level,policy.minLevel,policy.maxLevel);
        }
        electricRequirements=[&]{NativeMeshIdentity::ValidationBatch::Pause pause;
            return HydrodynamicsVofRegridRequirements::fromCellMinimumLevels(mesh,physical,std::move(minimum));}();
    }
    HydrodynamicsVofTargetDomainProvider targetDomain;
    if(input->fluid.csgGeometry){
        const auto definition=input->fluid.csgGeometry;const auto sampler=input->fluid.csgSampler;
        const auto excluded=input->fluid.csgExcludedPhases;const bool radial=input->fluid.radial;
        require(bool(sampler),"CSG target AMR requires its immutable Core source sampler");
        targetDomain=[definition,sampler,excluded,radial](const QuadtreeMesh& current){
            const auto solids=FittedDomainGeometry::selectFixedSolids(current,*definition,*sampler,excluded);
            FaceSegmentGraph graph;graph.build(current,[&](const QuadtreeMesh&,const Cell& cell){
                return solids.excludedPhaseByCell.at(cell.index())==0;});
            return NativeLeafDomain::build(current,graph,radial);
        };
    }
    // Tree preparation owns mutable candidates and an opaque target-domain
    // provider. No source proof survives it, including an enclosing caller.
    sourceMeshBatch.reset();
    data_->regrid=[&]{NativeMeshIdentity::ValidationBatch::Pause pause;
        return prepareHydrodynamicsVofRegrid(mesh,physical,*input->fluid.vofRegrid,
            electricRequirements?&*electricRequirements:nullptr,targetDomain);}();
    if(!data_->regrid->candidate){
        NativeMeshIdentity::ValidationBatch unchangedMeshBatch(mesh);
        data_->prepared=std::make_shared<JointAccepted>(physical,q,std::move(field));data_->prepared->validate(mesh);
        data_->receipt={&mesh,&physical,&*data_->step,nullptr,nullptr,&data_->regrid->report};
        data_->receipt.regridSeconds=std::chrono::duration<double>(Clock::now()-regridStart).count();return;
    }
    auto& grid=*data_->regrid->candidate;
    NativeMeshIdentity::ValidationBatch targetMeshBatch(*grid.mesh);
    stage="VOF-EHD same-End inventory transfer";const auto inventoryStart=Clock::now();
    auto targetCsg=csgBoundaryWork(data_->attemptCsgBoundary,data_->csgBoundaryReport.stages[2],bool(input->fluid.csgGeometry));
    auto transfer=LAMP::Interface::GeometricVofMeshTransfer::prepare(mesh,end,*grid.mesh,grid.domain,target,config.transport.geometry);
    const auto inventorySeconds=std::chrono::duration<double>(Clock::now()-inventoryStart).count();
    stage="VOF-EHD same-End fluid transfer";H::HydroVofMeshTransferConfig remapConfig;
    remapConfig.momentumDensityRepresentation=config.momentumDensityRepresentation;
    remapConfig.saddle=config.saddle;
    if(physical.vofPrescribedFlow()){
        const auto means=[&]{NativeMeshIdentity::ValidationBatch::Pause pause;
            return H::HydroVofMeshTransfer::boundaryNormalMeans(mesh,physical,*grid.mesh,grid.domain);}();
        const auto targetLayout=H::HydroDofLayout::build(*grid.mesh,grid.domain);
        auto targetFlow=prescribedFlow(*grid.mesh,targetLayout,*input,&means,targetCsg);
        // Remap authenticates source/target endpoint and ambient providers
        // and constructs a State; keep its before/after live checks literal.
        NativeMeshIdentity::ValidationBatch::Pause pause;
        data_->fluidRemap=H::HydroVofMeshTransfer::preparePrescribedFlow(mesh,physical,*grid.mesh,std::move(transfer),
            config.liquidDensity,config.gasDensity,std::move(targetFlow),remapConfig);
    }else{NativeMeshIdentity::ValidationBatch::Pause pause;
        data_->fluidRemap=H::HydroVofMeshTransfer::prepare(mesh,physical,*grid.mesh,std::move(transfer),config.liquidDensity,config.gasDensity,remapConfig);}
    stage="VOF-EHD same-End signed charge transfer";
    data_->chargeRemap=E::ElectricVofChargeMeshTransfer::prepare(mesh,end,q,*grid.mesh,data_->fluidRemap->interfaceTransfer(),input->controls);
    const auto& targetFluid=data_->fluidRemap->candidateState();const auto& targetVof=targetFluid.phaseGeometry().vofPublication();
    // Both actual meshes still exist here. Construct fresh target CSR now;
    // only scalar local blocks pass from the source, never an old owner chain.
    // The original full FP offer gate and next actual callback's typed checks
    // remain; a retry or failed joint commit never registers this candidate.
    std::optional<double> targetStressGeometrySeconds;
    const auto* targetFlow=targetFluid.vofPrescribedFlow();
    const bool targetGreenOnly=(config.liquidViscosity==0&&config.gasViscosity==0)
        ||(targetFlow&&!targetFlow->stationaryNoSlipVelocityDofs.empty());
    if(data_->pendingMaxwellGeometry){
        MaxwellGeometryFpState before;
        const auto& oldGeometry=data_->pendingMaxwellGeometry->geometry;
        const bool qualified=MaxwellGeometryFpState::read(before)&&before.same(data_->pendingMaxwellGeometry->environment)
            &&maxwellGeometryMode(oldGeometry,greenOnly)&&oldGeometry.layout().sameIdentity(physical.phaseGeometry().layout())
            && &oldGeometry.metrics().divergence().matrix()==&physical.phaseGeometry().metrics().divergence().matrix();
        if(qualified){
            stage="VOF-EHD target local stress geometry";
            const auto targetStressStarted=amrTimingStart();
            const bool sameMode=maxwellGeometryMode(oldGeometry,targetGreenOnly);
            auto targetGeometry=targetGreenOnly
                ?sameMode?H::HydroStrainOperator::buildRt0GreenRemapped(mesh,oldGeometry,*grid.mesh,
                    targetFluid.phaseGeometry().metrics(),data_->fluidRemap->interfaceTransfer().nativeOverlaps())
                    :H::HydroStrainOperator::buildRt0Green(*grid.mesh,targetFluid.phaseGeometry().metrics())
                :sameMode?H::HydroStrainOperator::buildStressCompatibleRemapped(mesh,oldGeometry,*grid.mesh,
                    targetFluid.phaseGeometry().metrics(),data_->fluidRemap->interfaceTransfer().nativeOverlaps())
                    :H::HydroStrainOperator::buildStressCompatible(*grid.mesh,targetFluid.phaseGeometry().metrics());
            const auto targetStressElapsed=amrTimingFinish(targetStressStarted);
            MaxwellGeometryFpState after;
            if(MaxwellGeometryFpState::read(after)&&before.same(after)){
                data_->targetMaxwellGeometry.emplace(RetainedMaxwellGeometry{std::move(targetGeometry),before});
                targetStressGeometrySeconds=targetStressElapsed;
            }
        }
    }
    stage="VOF-EHD target electrostatic snapshot";
    auto targetQ=data_->chargeRemap->chargeCandidate();
    auto snapshot=E::ElectricVofField::prepareElectrostaticSnapshot(*grid.mesh,targetVof,targetQ,input->materials,
        boundary(*grid.mesh,targetVof,*input,targetCsg),input->controls,E::ElectricVofField::CellGradientRecovery::Omit,&*data_->physicalEnd);
    data_->prepared=std::make_shared<JointAccepted>(targetFluid,targetQ,std::move(snapshot));
    const auto finalJointValidationStarted=amrTimingStart();
    data_->prepared->validate(*grid.mesh);
    const auto finalJointValidationSeconds=amrTimingFinish(finalJointValidationStarted);
    data_->receipt={grid.mesh.get(),&data_->fluidRemap->sourceState(),&*data_->step,&data_->fluidRemap->interfaceTransfer(),
        &*data_->fluidRemap,&data_->regrid->report};
    data_->receipt.inventoryTransferSeconds=inventorySeconds;
    data_->receipt.targetStressGeometrySeconds=targetStressGeometrySeconds;
    data_->receipt.finalJointValidationSeconds=finalJointValidationSeconds;
    data_->receipt.regridSeconds=std::chrono::duration<double>(Clock::now()-regridStart).count();
}
void VofElectrohydrodynamicsTransaction::commit(QuadtreeMesh& mesh)
{
    // No receipt from any caller can cross callback validation or swapStorage.
    NativeMeshIdentity::ValidationBatch::Pause pause;
    if(data_->input.method==VofElectrohydrodynamicsInput::Method::VofSharpSurface){
        require(data_->step.has_value()&&bool(data_->preparedSharp)&&!data_->committed,"no uncommitted complete sharp joint candidate");
        H::HydroTwoPhaseState physical=data_->acceptedSharp->fluid;data_->step->commit(mesh,physical);
        require(physical.sameIdentity(data_->step->candidateState())&&data_->preparedSharp->fluid.sameIdentity(physical),"sharp joint commit has stale physical End");
        data_->preparedSharp->validate(mesh);data_->acceptedSharp.swap(data_->preparedSharp);data_->committed=true;return;
    }
    require(data_->step.has_value()&&bool(data_->prepared)&&!data_->committed,"no uncommitted complete joint candidate");
    H::HydroTwoPhaseState physical=data_->accepted->fluid;data_->step->commit(mesh,physical);
    require(physical.sameIdentity(data_->step->candidateState()),"joint commit has stale physical End");
    if(data_->fluidRemap){
        require(data_->regrid&&data_->regrid->candidate&&data_->chargeRemap&&
            physical.sameIdentity(data_->fluidRemap->sourceState())&&
            data_->prepared->fluid.sameIdentity(data_->fluidRemap->candidateState())&&
            data_->prepared->charge.sameInstance(data_->chargeRemap->chargeCandidate()),"joint remap has unmatched actual candidate");
        auto& grid=*data_->regrid->candidate;data_->prepared->validate(*grid.mesh);
        require(sameClock(physical.timeLayer(),data_->prepared->fluid.timeLayer()),"joint regrid changed complete End clock");
        // All allocations, solves and throwing owner checks precede this swap.
        // swapStorage validates compatibility before its nonallocating swaps.
        mesh.swapStorage(*grid.mesh);
    }else{
        require(data_->prepared->fluid.sameIdentity(physical),"joint candidate differs from actual physical End");data_->prepared->validate(mesh);
    }
    data_->accepted.swap(data_->prepared);
    registerAcceptedCsgBoundary(data_->acceptedCsgBoundary,data_->attemptCsgBoundary,
        data_->accepted->fluid.phaseGeometry().vofPublication().domain());
    if(data_->fluidRemap)data_->acceptedMaxwellGeometry.swap(data_->targetMaxwellGeometry);
    else data_->acceptedMaxwellGeometry.swap(data_->pendingMaxwellGeometry);
    data_->targetMaxwellGeometry.reset();data_->pendingMaxwellGeometry.reset();data_->committed=true;
}
const HydrodynamicsVofRemapReceipt* VofElectrohydrodynamicsTransaction::remapReceipt() const noexcept
{return data_->committed&&data_->regrid?&data_->receipt:nullptr;}
void VofElectrohydrodynamicsTransaction::output(const QuadtreeMesh& mesh,Clock::time_point start,double compute,bool final)
{
    require(data_->committed,"output requires a committed joint state");
    if(data_->input.method==VofElectrohydrodynamicsInput::Method::VofSharpSurface){
        data_->writer.accepted(mesh,data_->acceptedSharp->fluid,data_->acceptedSharp->surface,data_->acceptedSharp->field,&*data_->step,start,compute,final);return;
    }
    data_->writer.accepted(mesh,data_->accepted->fluid,data_->accepted->charge,data_->accepted->field,&*data_->step,
        start,compute,final,remapReceipt(),data_->chargeRemap?&*data_->chargeRemap:nullptr,
        physicalChargeTransport(),physicalEndField(),csgBoundaryPreparation());
}
void VofElectrohydrodynamicsTransaction::releaseStep()
{
    data_->committed=false;data_->receipt={};data_->targetMaxwellGeometry.reset();data_->pendingMaxwellGeometry.reset();
    dropCsgBoundaryCache(data_->attemptCsgBoundary);data_->csgBoundaryReport={};
    // After a changed commit prepared holds the OLD accepted bundle. Destroy
    // every borrowed old field/receipt before regrid's old mesh storage.
    data_->preparedSharp.reset();data_->prepared.reset();data_->physicalEnd.reset();data_->physicalTransport.reset();data_->chargeRemap.reset();data_->fluidRemap.reset();data_->step.reset();data_->regrid.reset();
}
void VofElectrohydrodynamicsTransaction::clear(){if(data_){releaseStep();data_->acceptedSharp.reset();data_->accepted.reset();data_->acceptedMaxwellGeometry.reset();dropCsgBoundaryCache(data_->acceptedCsgBoundary);}}


