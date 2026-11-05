#include "Modules/Hydrodynamics/HydroVofStep.h"
#include "Modules/Hydrodynamics/Detail/HydroFiniteContraction.h"
#include "Core/Geometry/Detail/FiniteVolumeClosure.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/NativeMeshIdentity.h"
#include <algorithm>
#include <chrono>
#include <cfenv>
#include <cerrno>
#include <cstring>
#include <memory>
#if defined(_MSC_VER)
#include <float.h>
#include <intrin.h>
#endif
#if defined(_OPENMP)
#include <omp.h>
#endif
#if defined(_M_X64) || defined(__SSE2__)
#include <xmmintrin.h>
#endif
#include <exception>
#include <future>
#include <system_error>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>
#include <sstream>
#include <utility>

namespace LAMP::Hydrodynamics {
namespace {
namespace I=LAMP::Interface;
namespace FV=LAMP::Geometry::Detail::FiniteVolume;
using Clock=I::GeometricVofTimeLayer;
void require(bool ok,const char* why){if(!ok)throw std::invalid_argument(std::string("[HydroVofStep] ")+why);}
// Optional rate-row preparation only. The original ordered diagnostics below
// remains the literal fallback and owns every norm, work sum and gate.
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE) \
    && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST) \
    && !defined(__AVX__) && !defined(__AVX2__) && !defined(__AVX512F__) && !defined(__FMA__)
#pragma float_control(precise,on,push)
#pragma float_control(except,on)
#pragma fenv_access(on)
struct DiagnosticRateState {
    std::fenv_t native{};unsigned csr=0;int error=0;bool valid=false;
    bool capture() noexcept {error=errno;csr=_mm_getcsr();valid=std::fegetenv(&native)==0;return valid;}
    void restore() const noexcept {
        if(valid&&std::fesetenv(&native)!=0)std::terminate();
        _mm_setcsr(csr);errno=error;
    }
    bool enterCaller(const DiagnosticRateState& caller) const noexcept {
        if(std::fesetenv(&caller.native)!=0)return false;
        _mm_setcsr(caller.csr);errno=caller.error;return true;
    }
};
struct DiagnosticRateRestore {
    DiagnosticRateState state;
    DiagnosticRateRestore() noexcept {state.capture();}
    ~DiagnosticRateRestore() noexcept {state.restore();}
};
// One coarse HF task consumes the App-configured thread budget; no outer
// OpenMP region may disable the existing viscosity kernels.
struct CapillaryTaskThreads {
    int previous=0,dynamic=0;bool active=false,fixed=false;
    explicit CapillaryTaskThreads(int count,bool fixDynamic=false):fixed(fixDynamic){
        if(count<=0)return;previous=omp_get_max_threads();dynamic=omp_get_dynamic();
        omp_set_num_threads(count);if(fixed)omp_set_dynamic(0);active=true;
    }
    ~CapillaryTaskThreads() noexcept {if(active){omp_set_num_threads(previous);if(fixed)omp_set_dynamic(dynamic);}}
};
bool capillaryOverlapQualified(DiagnosticRateState& entry,int& budget) noexcept {
    DiagnosticRateRestore restore;entry=restore.state;budget=omp_get_max_threads();unsigned control=0;
    if(!entry.valid||omp_in_parallel()||budget<2||sizeof(long double)!=8
        ||std::numeric_limits<long double>::radix!=2||std::numeric_limits<long double>::digits!=53
        ||std::fegetround()!=FE_TONEAREST||_controlfp_s(&control,0,0)!=0
        ||(control&_MCW_EM)!=_MCW_EM||(control&_MCW_RC)!=_RC_NEAR
        ||(entry.csr&0x1f80u)!=0x1f80u||(entry.csr&0xe040u)!=0)return false;
#ifdef _EM_AMBIGUOUS
    if(control&_EM_AMBIGUOUS)return false;
#endif
    return true;
}
bool capillarySameControls(const DiagnosticRateState& a,const DiagnosticRateState& b) noexcept {
    auto x=a.native,y=b.native;x._Fe_stat=0;y._Fe_stat=0;
    return a.valid&&b.valid&&std::memcmp(&x,&y,sizeof x)==0&&(a.csr&~0x3fu)==(b.csr&~0x3fu);
}
void capillaryMergeFlags(const DiagnosticRateState& main,const DiagnosticRateState& worker) noexcept {
    auto merged=main.native;if(worker.valid)merged._Fe_stat|=worker.native._Fe_stat;
    if(main.valid&&std::fesetenv(&merged)!=0)std::terminate();
    _mm_setcsr(main.csr|(worker.valid?(worker.csr&0x3fu):0));errno=main.error;
}
class DiagnosticMomentumRates {
    std::unique_ptr<long double[]> values_;
public:
    DiagnosticMomentumRates()=default;
    DiagnosticMomentumRates(DiagnosticMomentumRates&&) noexcept=default;
    DiagnosticMomentumRates(const DiagnosticMomentumRates&)=delete;
    DiagnosticMomentumRates& operator=(const DiagnosticMomentumRates&)=delete;
    explicit operator bool() const noexcept {return bool(values_);}
    long double operator[](int u) const noexcept {return values_[static_cast<std::size_t>(u)];}
    void allocate(std::size_t n){values_=std::make_unique<long double[]>(n);}
    void store(int u,long double value) noexcept {values_[static_cast<std::size_t>(u)]=value;}
    __declspec(noinline) static void releasePreservingState(long double* values) noexcept {
        // Baseline x64 FXSAVE/FXRSTOR has no capture-failure return. Keep the
        // complete x87/MXCSR state across valid POD deallocation without a new
        // late fegetenv/terminate path. This isolated void helper returns no
        // XMM value whose register could conflict with the restored snapshot.
        const int error=errno;alignas(16) unsigned char state[512];
        _fxsave64(state);delete[] values;_fxrstor64(state);errno=error;
    }
    ~DiagnosticMomentumRates() noexcept {
        if(values_)releasePreservingState(values_.release());
    }
};
bool diagnosticRateSameState(const DiagnosticRateState& a,const DiagnosticRateState& b) noexcept {
    return std::memcmp(&a.native,&b.native,sizeof a.native)==0&&a.csr==b.csr&&a.error==b.error;
}
DiagnosticMomentumRates prepareDiagnosticMomentumRates(int count,const double* beginMass,const double* endMass,
    const double* beginVelocity,const double* u1,double dt,HydroVofStageReport::MomentumRateRows& receipt) noexcept
{
    DiagnosticRateRestore restore;
    if(!restore.state.valid)return {}; // No clock/FP timing without a restorable entry.
    const auto start=std::chrono::steady_clock::now();
    const auto finish=[&] {receipt.seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();};
    DiagnosticMomentumRates result;
    receipt.requested=count>1&&omp_get_max_threads()>1;
    if(!receipt.requested){finish();return result;}
    receipt.fallback=true;
    const auto& caller=restore.state;
    unsigned control=0;
    if(!caller.valid||caller.error==EDOM||caller.error==ERANGE||caller.error==EINVAL
        ||omp_in_parallel()||sizeof(double)!=8||sizeof(long double)!=8
        ||std::numeric_limits<double>::radix!=2||std::numeric_limits<double>::digits!=53
        ||std::numeric_limits<long double>::digits!=53||std::numeric_limits<long double>::max_exponent!=1024
        ||std::numeric_limits<long double>::min_exponent!=-1021||!std::numeric_limits<double>::is_iec559
        ||std::fegetround()!=FE_TONEAREST||_controlfp_s(&control,0,0)!=0
        ||(control&_MCW_EM)!=_MCW_EM||(control&_MCW_RC)!=_RC_NEAR
        ||(caller.csr&0x1f80u)!=0x1f80u||(caller.csr&0xe040u)!=0){finish();return result;}
#ifdef _EM_AMBIGUOUS
    if(control&_EM_AMBIGUOUS){finish();return result;}
#endif
    receipt.eligible=true;
    struct Worker {DiagnosticRateState original;bool failed=false;};
    const int team=std::min(count,omp_get_max_threads());
    constexpr std::size_t byteLimit=64u*1024u*1024u,integerReserve=128u*1024u;
    std::size_t bytes=1024;
    const auto account=[&](std::size_t n,std::size_t width) {
        if(bytes>byteLimit||n>(byteLimit-bytes)/width)return false;bytes+=n*width;return true;
    };
    // No whole-domain arbitrary-integer inventory. Two native products over
    // one dt need less integer space than the existing four-factor KEP bound;
    // retain that conservative per-active-worker reserve in the peak cap.
    if(team<2||!account(static_cast<std::size_t>(count),sizeof(long double))
        ||!account(static_cast<std::size_t>(team),sizeof(Worker)+integerReserve)){finish();return result;}
    try {
        auto workers=std::make_unique<Worker[]>(static_cast<std::size_t>(team));
        result.allocate(static_cast<std::size_t>(count));
        int actualTeam=0;
#pragma omp parallel num_threads(team)
        {
            const int tid=omp_get_thread_num();
            if(tid==0)actualTeam=omp_get_num_threads();
            auto& local=workers[static_cast<std::size_t>(tid)];
            const bool captured=local.original.capture();
            const bool ready=captured&&local.original.enterCaller(caller);
            if(!ready)local.failed=true;
#pragma omp for schedule(static)
            for(int u=0;u<count;++u){
                if(!ready||local.failed)continue;
                try {
                    const long double m0=beginMass[u],m1=endMass[u],u0=beginVelocity[u];
                    Detail::HydroFiniteContraction momentumRate;
                    momentumRate.addQuotientProduct({m1,u1[u]},dt);
                    momentumRate.addQuotientProduct({-1.L,m0,u0},dt);
                    result.store(u,momentumRate.rounded());
                }catch(...){local.failed=true;}
            }
            if(ready){DiagnosticRateState after;
                if(!after.capture()||!diagnosticRateSameState(after,caller))local.failed=true;}
            local.original.restore();
        }
        bool failed=actualTeam<2||actualTeam>team;
        for(int t=0;t<actualTeam&&t<team;++t)failed=failed||workers[static_cast<std::size_t>(t)].failed;
        if(!failed){receipt.used=true;receipt.fallback=false;receipt.rows=static_cast<std::size_t>(count);
            receipt.workers=actualTeam;receipt.scratchBytes=bytes;finish();return result;}
    }catch(...){/* Optional preparation never owns the original first error. */}
    // All workers have joined. Scratch destruction precedes restoration of the
    // exact entry state; the caller then executes the complete original u loop.
    finish();return DiagnosticMomentumRates{};
}
#if !defined(_M_FP_STRICT) || !_M_FP_STRICT
#pragma fenv_access(off)
#endif
#pragma float_control(pop)
#else
class DiagnosticMomentumRates {
public:
    explicit operator bool() const noexcept {return false;}
    long double operator[](int) const noexcept {return 0;}
};
DiagnosticMomentumRates prepareDiagnosticMomentumRates(int,const double*,const double*,const double*,const double*,double,
    HydroVofStageReport::MomentumRateRows&) noexcept {return {};}
#endif
struct CapillaryScale {bool available=false;double minimumSpacing=0,logRate=0;};
CapillaryScale capillaryScale(const I::GeometricVofPublication& source,const HydroVofAcConfig& config)
{
    CapillaryScale result;
    if(config.surfaceTension>0){double h=std::numeric_limits<double>::infinity();
        for(const auto& c:source.plicCells())if(c.inventory[0]>0&&c.inventory[1]>0)h=std::min(h,static_cast<double>(std::min(c.rectangle.upper.x-c.rectangle.lower.x,c.rectangle.upper.y-c.rectangle.lower.y)));
        if(std::isfinite(h)){const double rho=std::max(config.liquidDensity,config.gasDensity),logMean=std::log(rho)+std::log(.5*(config.liquidDensity/rho+config.gasDensity/rho));
            result.available=true;result.minimumSpacing=h;
            result.logRate=.5*(std::log(std::acos(-1.))+std::log(config.surfaceTension)-logMean-3*std::log(h));
        }
    }
    return result;
}

double relative(long double error,long double scale)
{
    require(std::isfinite(error)&&std::isfinite(scale)&&error>=0&&scale>=0,"nonfinite residual scale");
    if(scale==0){require(error==0,"nonzero residual with zero scale");return 0;}
    const auto value=static_cast<double>(error/scale);require(std::isfinite(value),"unrepresentable relative residual");return value;
}
HydroNativeValues midpoint(const HydroNativeValues& before,const HydroNativeValues& endpoint)
{
    std::vector<double> u(before.normalVelocity().size());
    for(std::size_t i=0;i<u.size();++i)u[i]=.5*before.normalVelocity()[i]+.5*endpoint.normalVelocity()[i];
    // The original saddle pressure is already the stage midpoint pressure.
    return {before.layout(),std::move(u),endpoint.pressure()};
}
HydroNativeValues endpointGuess(const HydroNativeValues& before,const HydroNativeValues& advector)
{
    std::vector<double> u(before.normalVelocity().size());
    for(std::size_t i=0;i<u.size();++i){u[i]=2*advector.normalVelocity()[i]-before.normalVelocity()[i];require(std::isfinite(u[i]),"nonfinite endpoint guess");}
    // The advector carries the preceding actual AC midpoint pressure.
    // Transport consumes only U; the native solve uses this pressure as a hint.
    return {before.layout(),std::move(u),advector.pressure()};
}
double correctionRelative(const HydroNativeValues& previous,const HydroNativeValues& solved,const HydroMomentumMass& mass,double reference)
{
    long double error=0,norm=0,total=0;
    for(std::size_t u=0;u<mass.masses().size();++u){const long double m=mass.masses()[u],v=solved.normalVelocity()[u],d=v-previous.normalVelocity()[u];
        require(m>0&&std::isfinite(m),"nonpositive correction metric mass");error+=m*d*d;norm+=m*v*v;total+=m;}
    const auto scale=std::max(std::sqrt(norm),static_cast<long double>(reference)*std::sqrt(total));
    require(std::isfinite(scale)&&scale>0,"nonpositive correction norm scale");return relative(std::sqrt(error),scale);
}

double mismatch(const QuadtreeMesh& mesh,const HydroPhaseGeometry& phase,
    const HydroNativeValues& advector,const HydroNativeValues& solved,double dt)
{
    const auto& layout=phase.layout();
    std::vector<Vector2D> first(layout.pressureCount()),last(layout.pressureCount()),weights(layout.pressureCount());
    for(int u=0;u<layout.velocityCount();++u){const auto f=layout.domain().face(mesh,layout.nativeFace(u));
        const double area=phase.metrics().faceMeasures()[u];const bool x=layout.normalAxis(u)==HydroDofLayout::NormalAxis::First;
        for(int native:{f.negativeCell,f.positiveCell})if(native>=0){const int p=layout.pressureDof(native);if(p<0)continue;
            if(x){first[p].x+=area*advector.normalVelocity()[u];last[p].x+=area*solved.normalVelocity()[u];weights[p].x+=area;}
            else{first[p].y+=area*advector.normalVelocity()[u];last[p].y+=area*solved.normalVelocity()[u];weights[p].y+=area;}}
    }
    double result=0;for(int p=0;p<layout.pressureCount();++p){const auto& c=phase.vofPublication().plicCells()[layout.nativeCell(p)];
        if(!(c.inventory[0]>0&&c.inventory[1]>0))continue;
        const double h=std::min(c.rectangle.upper.x-c.rectangle.lower.x,c.rectangle.upper.y-c.rectangle.lower.y);
        require(h>0&&weights[p].x>0&&weights[p].y>0,"incomplete native interface velocity support");
        const auto dx=(last[p].x-first[p].x)/weights[p].x,dy=(last[p].y-first[p].y)/weights[p].y;
        const auto value=dt/h*std::abs(dx*c.normal.x+dy*c.normal.y);require(std::isfinite(value),"nonfinite interface mismatch");result=std::max(result,value);
    }return result;
}
// Explicit open native mapping; the old closed helper remains untouched.
I::GeometricVofAdvectorSlab prescribedSlab(const QuadtreeMesh& mesh,const I::GeometricVofPublication& source,
    Clock end,const HydroNativeValues& before,const HydroNativeValues& guess,const HydroVofPrescribedFlow& flow,
    const I::GeometricVofAdvectorSlab::Boundary* retainedBoundary)
{
    const auto& layout=before.layout();before.validateFor(layout,mesh);guess.validateFor(layout,mesh);
    require(layout.domain().meshIdentity()==source.domain().meshIdentity()&&
        layout.domain().axisymmetric()==source.domain().axisymmetric(),"different open VOF native source binding");
    require(layout.domain().cellRowByIndex()==source.domain().cellRowByIndex()&&
        layout.domain().cellIndexByRow()==source.domain().cellIndexByRow(),"different open VOF active native map");
    require(layout.domain().nativeFaceCount()==source.domain().nativeFaceCount(),"different open VOF native face layout");
    std::vector<double> first(source.domain().nativeFaceCount()),last(first.size());
    for(int row=0;row<source.domain().nativeFaceCount();++row){const int u=layout.velocityDof(row);
        const auto face=source.domain().face(mesh,row);
        if(face.negativeCell<0&&face.positiveCell<0){require(u<0,"excluded open VOF face has a Hydro velocity DOF");continue;}
        require(u>=0&&layout.nativeFace(u)==row,"invalid open VOF native normal mapping");first[row]=before.normalVelocity()[u];last[row]=guess.normalVelocity()[u];}
    if(retainedBoundary)return I::GeometricVofAdvectorSlab::prepareBoundary(mesh,source,*retainedBoundary,end,first,last);
    const auto boundary=I::GeometricVofAdvectorSlab::Boundary::prepare(mesh,source,flow.rulesByNativeRow);
    return I::GeometricVofAdvectorSlab::prepareBoundary(mesh,source,boundary,end,first,last);
}
// Only actual immutable value owners can bypass an opaque callback boundary.
// The State's retained Boundary separately binds its source and full clock.
bool steadyPrescribedFlow(const QuadtreeMesh& mesh,const HydroDofLayout& layout,const HydroVofPrescribedFlow& flow)
{
    using Role=HydroVofSteadyVelocity::ClockRole;
    const auto typed=[&](const auto& query,Role role){
        if(!query)return true;
        const auto* owner=query.template target<HydroVofSteadyVelocity>();
        if(!owner)return false;
        require(owner->matches(layout,role),"steady VOF vector has another layout or clock role");
        owner->validateFor(mesh,layout,role);return true;
    };
    return typed(flow.endpointVelocity,Role::Endpoint)&&typed(flow.pressureOutletBackflowMomentum,Role::FullMidpoint);
}
// Explicit F16 opt-in: weak incoming momentum is not prescribed normal data.
HydroConvectionBoundary prescribedOpenMomentum(const QuadtreeMesh& mesh,const HydroDualMaterial& dual,
    const HydroVofPrescribedFlow& flow)
{
    require(flow.pressureOutletBackflowModel.has_value(),"missing explicit VOF open momentum model");
    const auto model=*flow.pressureOutletBackflowModel;const auto& layout=dual.layout();
    std::vector<Vector2D> ambient(layout.domain().nativeFaceCount());
    std::vector<int> ports;
    for(const auto& port:flow.pressureOutlets){require(port.velocityDof>=0&&port.velocityDof<layout.velocityCount(),"invalid VOF open momentum port");
        ports.push_back(port.velocityDof);
        if(model==HydroVofPressureOutletMomentumModel::Ambient){
            require(bool(flow.pressureOutletBackflowMomentum),"missing actual VOF ambient momentum query");
            ambient[layout.nativeFace(port.velocityDof)]=flow.pressureOutletBackflowMomentum(layout.nativeFace(port.velocityDof),dual.timeLayer());}}
    // The production producer queries the actual full gross once and invokes
    // these same endpoint callbacks only for genuinely incoming nonfree rows.
    return HydroConvectionBoundary::bindVofOpenProductionMomentum(mesh,dual,ports,ambient,
        flow.endpointVelocity,dual.sources(),dual.timeLayer(),model);
}
HydroConvectionBoundary prescribedInlet(const QuadtreeMesh& mesh,const HydroDualMaterial& dual,
    const HydroVofPrescribedFlow& flow)
{
    if(flow.pressureOutletBackflowModel)return prescribedOpenMomentum(mesh,dual,flow);
    const auto& material=dual.material();const auto& volume=material.vofTransport();const auto& layout=dual.layout();
    const auto& segments=dual.geometry().segments();std::vector<Vector2D> first(segments.size()),last(segments.size());
    std::vector<bool> queried(layout.domain().nativeFaceCount()),free(layout.velocityCount());
    std::vector<Vector2D> sourceVector(queried.size()),endVector(queried.size());
    for(const auto& port:flow.pressureOutlets)free[port.velocityDof]=true;
    for(std::size_t g=0;g<segments.size();++g){const auto& segment=segments[g];
        if(segment.negativeVelocity>=0&&segment.positiveVelocity>=0)continue;
        const auto q=(segment.negativeVelocity>=0?1.L:-1.L)*dual.faceRates()[g];if(!(q<0))continue;
        require(segment.primalVelocity>=0&&segment.primalVelocity<layout.velocityCount(),"inflow segment has no actual native face source");
        require(!free[segment.primalVelocity],"free VOF pressure outlet cannot consume reservoir backflow");
        const int row=layout.nativeFace(segment.primalVelocity);
        if(!queried[row]){sourceVector[row]=flow.endpointVelocity(row,volume.sourceTime());endVector[row]=flow.endpointVelocity(row,volume.endTime());queried[row]=true;}
        first[g]=sourceVector[row];last[g]=endVector[row];
    }
    return HydroConvectionBoundary::bindVofEndpointVelocity(mesh,dual,first,last,dual.sources(),dual.timeLayer());
}
HydroConvectionBoundary prescribedUpwind(const QuadtreeMesh& mesh,const HydroDualMaterial& dual,
    const HydroConvectionBoundary& inlet,const HydroNativeValues& values)
{
    const auto& segments=dual.geometry().segments();std::vector<double> traces(segments.size());
    for(std::size_t g=0;g<segments.size();++g){const auto& segment=segments[g];
        if(segment.negativeVelocity>=0&&segment.positiveVelocity>=0)continue;
        const auto q=(segment.negativeVelocity>=0?1.L:-1.L)*dual.faceRates()[g];const int owner=segment.negativeVelocity>=0?segment.negativeVelocity:segment.positiveVelocity;
        require(owner>=0,"exterior VOF segment has no actual native owner");
        traces[g]=q<0?inlet.traces()[g]:q>0?values.normalVelocity()[owner]:0.;}
    return HydroConvectionBoundary::bind(mesh,dual,traces,dual.sources(),dual.timeLayer());
}
struct Stage {
    I::GeometricVofTransport volume;
    HydroFrozenMidpoint frozen;
    HydroNativeValues solvedMidpoint;
    std::optional<HydroBalancedCsfLoad> capillary;
    HydroVofStageReport report;
};
Stage stage(const QuadtreeMesh& mesh,const HydroTwoPhaseState& accepted,Clock end,
    const HydroVofAcConfig& config,const HydroMomentumMass& beginMass,const HydroNativeValues& advector,
    const HydroPressureLoad& pressure,bool prescribed,const HydroPhaseViscosity* previousViscosity=nullptr)
{
    const auto* flow=prescribed?accepted.vofPrescribedFlow():nullptr;
    const bool openMomentum=flow&&flow->pressureOutletBackflowModel.has_value();
    using Wall=std::chrono::steady_clock;const auto start=Wall::now();auto mark=start;
    const auto seconds=[](Wall::duration duration){
        const int error=errno;std::fenv_t environment{};const bool saved=std::fegetenv(&environment)==0;
#if defined(_M_X64) || defined(__SSE2__)
        const unsigned csr=_mm_getcsr();
#endif
        const auto value=std::chrono::duration<double>(duration).count();
        if(saved&&std::fesetenv(&environment)!=0)std::terminate();
#if defined(_M_X64) || defined(__SSE2__)
        _mm_setcsr(csr);
#endif
        errno=error;return value;};
    const auto part=[&]{const auto now=Wall::now();const auto value=seconds(now-mark);mark=now;return value;};
    HydroVofStageReport report;const auto& geometry=accepted.momentumGeometry();const auto& before=accepted.values();const auto& layout=before.layout();
    const auto& old=accepted.phaseGeometry().vofPublication();
    const auto* retainedBoundary=accepted.vofTransportBoundary();
    const bool steadyFlow=flow&&retainedBoundary&&steadyPrescribedFlow(mesh,layout,*flow);
    const auto predictedEndpoint=endpointGuess(before,advector);
    const auto slab=flow?prescribedSlab(mesh,old,end,before,predictedEndpoint,*flow,retainedBoundary)
        :prepareHydroVofAdvectorSlab(mesh,old,end,before,predictedEndpoint);
    const auto volume=[&]{try{
        if(config.transportProvider){NativeMeshIdentity::ValidationBatch::Pause pause;
            auto prepared=config.transportProvider(mesh,old,slab,config.liquidDensity,config.gasDensity,config.transport);
            prepared.validateFor(mesh,old,slab,config.liquidDensity,config.gasDensity);return prepared;}
        return flow
        ?I::GeometricVofTransport::prepareLocalPlicBoundarySweep(mesh,old,slab,config.liquidDensity,config.gasDensity,config.transport)
        :I::GeometricVofTransport::prepareLocalPlicSweep(mesh,old,slab,config.liquidDensity,config.gasDensity,config.transport);}
        catch(const I::GeometricVofCourantExceeded& error){throw HydroTwoPhaseStepRetry(HydroTwoPhaseRetryReason::Transport,error.what());}}();
    report.timing.transport=part();
    const auto materialStart=mark;
    const auto middle=HydroPhaseGeometry::build(mesh,geometry.metrics(),volume.midpoint());
    const auto middleEnd=Wall::now();
    std::optional<HydroVofSouthEven90Boundary> preparedCapillaryBoundary;
    std::optional<HydroBalancedCsfLoad> capillary;
    double earlyBoundarySeconds=0,backgroundCapillarySeconds=0;
    const auto capillaryKernel=[&]{
        return preparedCapillaryBoundary
            ?HydroBalancedCsfLoad::localPlicHeightFunctionCapillary(mesh,pressure,middle,config.surfaceTension,*preparedCapillaryBoundary,config.curvature)
            :HydroBalancedCsfLoad::localPlicHeightFunctionCapillary(mesh,pressure,middle,config.surfaceTension,config.curvature);
    };
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE) && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST) && !defined(__AVX__) && !defined(__AVX2__) && !defined(__AVX512F__) && !defined(__FMA__)
    DiagnosticRateState overlapEntry,capillaryExit;int overlapBudget=0;bool workerEntered=false;
    std::future<std::optional<HydroBalancedCsfLoad>> pendingCapillary;
    // Drain a task on any main-path failure before its borrowed owners die.
    // Retain that main failure's FP/errno rather than a speculative HF failure.
    struct CapillaryFailureJoin {
        std::future<std::optional<HydroBalancedCsfLoad>>& pending;
        int failures=std::uncaught_exceptions();
        ~CapillaryFailureJoin() noexcept {
            if(std::uncaught_exceptions()>failures&&pending.valid()){
                DiagnosticRateRestore restore;try{(void)pending.get();}catch(...){}
            }
        }
    } failureJoin{pendingCapillary};
    std::optional<HydroBalancedCsfLoad> backgroundCapillary;
    const auto startCapillary=[&]{
        if(config.surfaceTension>0&&(config.liquidViscosity!=0||config.gasViscosity!=0)
            &&capillaryOverlapQualified(overlapEntry,overlapBudget)){
            // Opaque boundary laws run on the main thread before the task.
            // In open flow this deliberately precedes material preparation;
            // independent failure priority changes, not the physical kernels.
            if(config.southEven90Boundary){const auto boundaryStarted=Wall::now();
                {NativeMeshIdentity::ValidationBatch::Pause pause;
                    preparedCapillaryBoundary=config.southEven90Boundary(mesh,volume.midpoint());}
                earlyBoundarySeconds=seconds(Wall::now()-boundaryStarted);
            }
            if(capillaryOverlapQualified(overlapEntry,overlapBudget))try{
                pendingCapillary=std::async(std::launch::async,[&,entry=overlapEntry,workers=overlapBudget/2]() ->std::optional<HydroBalancedCsfLoad>{
                    DiagnosticRateRestore restore;CapillaryTaskThreads threads(workers,true);
                    if(!restore.state.valid||!restore.state.enterCaller(entry)||omp_in_parallel())return std::nullopt;
                    workerEntered=true;NativeMeshIdentity::ValidationBatch validationBatch(mesh);
                    const auto started=Wall::now();
                    try{auto load=capillaryKernel();backgroundCapillarySeconds=seconds(Wall::now()-started);
                        capillaryExit.capture();return load;}
                    catch(...){backgroundCapillarySeconds=seconds(Wall::now()-started);capillaryExit.capture();throw;}
                });
            }catch(const std::system_error& error){
                if(error.code()!=std::make_error_code(std::errc::resource_unavailable_try_again))throw;
            }catch(const std::bad_alloc&){/* No task started: retain the serial calculation. */}
        }
    };
    const auto joinCapillary=[&]{
        if(!pendingCapillary.valid())return;
        DiagnosticRateState afterMain;afterMain.capture();
        try{backgroundCapillary=pendingCapillary.get();afterMain.restore();}
        catch(...){
            capillaryMergeFlags(afterMain,capillaryExit);
            if(capillaryExit.valid&&capillaryExit.error!=overlapEntry.error)errno=capillaryExit.error;
            throw;
        }
    };
    // Generic open providers still join before their first opaque query.
    // Actual steady value owners also permit inlet/convection consumption.
    // Closed flow retains the same task's later viscosity overlap.
    if(flow)startCapillary();
#endif
    const auto preparationStarted=Wall::now();
    Wall::time_point phaseGeometryEnd,materialTransportEnd;
    const auto dual=[&]{
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE) && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST) && !defined(__AVX__) && !defined(__AVX2__) && !defined(__AVX512F__) && !defined(__FMA__)
        CapillaryTaskThreads threads(pendingCapillary.valid()?overlapBudget-overlapBudget/2:0);
        try{
#endif
        const auto endpoint=HydroPhaseGeometry::build(mesh,geometry.metrics(),volume.endpoint());
        phaseGeometryEnd=Wall::now();
        const auto material=HydroMaterialTransport::prepareGeometricVof(mesh,geometry,accepted.phaseGeometry(),middle,endpoint,volume,accepted.boundary());
        materialTransportEnd=Wall::now();
        const auto mass=[&](const HydroPhaseGeometry& phase){
            return HydroMomentumMass::buildCellAverageProjection(mesh,geometry,phase,config.liquidDensity,config.gasDensity);
        };
        return HydroDualMaterial::prepare(mesh,material,beginMass,mass(middle),mass(endpoint));
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE) && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST) && !defined(__AVX__) && !defined(__AVX2__) && !defined(__AVX512F__) && !defined(__FMA__)
        }catch(...){
            DiagnosticRateRestore failure;
            if(pendingCapillary.valid())try{(void)pendingCapillary.get();}catch(...){}
            throw;
        }
#endif
    }();
    const auto dualMassEnd=Wall::now();
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE) && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST) && !defined(__AVX__) && !defined(__AVX2__) && !defined(__AVX512F__) && !defined(__FMA__)
    if(flow&&!steadyFlow)joinCapillary();
    std::optional<CapillaryTaskThreads> convectionThreads;
    if(pendingCapillary.valid())convectionThreads.emplace(overlapBudget-overlapBudget/2);
#endif
    const auto convectionStarted=Wall::now();
    const auto sources=dual.sources();const auto clock=dual.timeLayer();
    std::optional<HydroConvectionBoundary> inlet;std::vector<int> freeOutlets;
    const auto bindInlet=[&]{inlet=prescribedInlet(mesh,dual,*flow);
        for(const auto& port:flow->pressureOutlets)freeOutlets.push_back(port.velocityDof);};
    if(flow){
        if(steadyFlow)bindInlet();
        else {NativeMeshIdentity::ValidationBatch::Pause pause;bindInlet();}
    }
    const auto convectionBindingEnd=Wall::now();
    const auto convection=openMomentum?HydroKepConvection::buildVofOpenProductionFlow(mesh,dual,accepted.boundary(),*inlet)
        :flow?HydroKepConvection::buildVofPrescribedFlow(mesh,dual,accepted.boundary(),*inlet,freeOutlets)
        :HydroKepConvection::buildImpermeable(mesh,dual,accepted.boundary());
    const auto convectionOperatorEnd=Wall::now();
    const auto boundary=openMomentum?*inlet:flow?prescribedUpwind(mesh,dual,*inlet,advector)
        :HydroConvectionBoundary::bind(mesh,dual,std::vector<double>(dual.faceRates().size()),sources,clock);
    const auto convectionBoundaryEnd=Wall::now();
    // Arbitrary incoming endpoint providers still run without a receipt.
    // A role/layout-validated steady owner has no callback or mutable data.
    const auto apply=[&]{return openMomentum?convection.applyVofOpenProductionLoads(mesh,advector,sources,clock)
        :convection.apply(mesh,advector,boundary,sources,clock);};
    const auto applied=[&]{if(steadyFlow)return apply();
        NativeMeshIdentity::ValidationBatch::Pause pause;return apply();}();
    const auto convectionEnd=Wall::now();
    if(openMomentum){auto& r=report.openConvection;r.available=true;
        r.boundaryBindingSeconds=seconds(convectionBindingEnd-convectionStarted);
        r.operatorSeconds=seconds(convectionOperatorEnd-convectionBindingEnd);
        r.advectorBoundarySeconds=seconds(convectionBoundaryEnd-convectionOperatorEnd);
        r.advectorApplySeconds=seconds(convectionEnd-convectionBoundaryEnd);r.advector=applied.openApplyReport();}
    // Deliberate order change: actual Maxwell preparation precedes viscosity so
    // canonical A geometry can lend its already computed raw CSR. Error priority
    // changes; neither physical operator or any acceptance gate changes.
    std::optional<HydroMaxwellLoad> maxwell;
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE) && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST) && !defined(__AVX__) && !defined(__AVX2__) && !defined(__AVX512F__) && !defined(__FMA__)
    convectionThreads.reset();
    if(config.maxwellLoad)joinCapillary(); // No task can span an opaque Maxwell provider.
#endif
    if(config.maxwellLoad){NativeMeshIdentity::ValidationBatch::Pause pause;
        maxwell=config.maxwellLoad(mesh,volume,middle,pressure);}
    const auto maxwellEnd=config.maxwellLoad?Wall::now():convectionEnd;
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE) && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST) && !defined(__AVX__) && !defined(__AVX2__) && !defined(__AVX512F__) && !defined(__FMA__)
    if(!flow)startCapillary();
#endif
    const auto viscosityStarted=Wall::now();
    const auto viscosity=[&]{
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE) && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST) && !defined(__AVX__) && !defined(__AVX2__) && !defined(__AVX512F__) && !defined(__FMA__)
        CapillaryTaskThreads threads(pendingCapillary.valid()?overlapBudget-overlapBudget/2:0);
        try{
#endif
        return config.liquidViscosity==0&&config.gasViscosity==0
        ?HydroPhaseViscosity::zero(mesh,middle)
        :flow&&!flow->stationaryNoSlipVelocityDofs.empty()
            ?previousViscosity?HydroPhaseViscosity::buildVofCellAverage(mesh,middle,config.liquidViscosity,config.gasViscosity,
                flow->stationaryNoSlipVelocityDofs,*previousViscosity)
                :HydroPhaseViscosity::buildVofCellAverage(mesh,middle,config.liquidViscosity,config.gasViscosity,
                    flow->stationaryNoSlipVelocityDofs,geometry)
        :previousViscosity?HydroPhaseViscosity::buildVofCellAverage(mesh,middle,config.liquidViscosity,config.gasViscosity,*previousViscosity)
            :maxwell&&!maxwell->isSharp()&&maxwell->geometry().hasCanonicalRawStrain()
                ?HydroPhaseViscosity::buildVofCellAverage(mesh,middle,config.liquidViscosity,config.gasViscosity,maxwell->geometry())
                :HydroPhaseViscosity::buildVofCellAverage(mesh,middle,config.liquidViscosity,config.gasViscosity);
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE) && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST) && !defined(__AVX__) && !defined(__AVX2__) && !defined(__AVX512F__) && !defined(__FMA__)
        }catch(...){
            // Future destruction also joins, but explicitly drain it while
            // retaining the original viscosity failure and its FP/errno state.
            DiagnosticRateRestore failure;
            if(pendingCapillary.valid())try{(void)pendingCapillary.get();}catch(...){}
            throw;
        }
#endif
    }();
    report.viscosityProduct=viscosity.viscosity().productReport();
    report.viscosityProductReportAvailable=true;
    report.viscosityProductPhaseIdentity=viscosity.geometry().vofIdentity().value();
    report.viscosityProductTime=viscosity.geometry().timeLayer();
    const auto viscosityEnd=Wall::now();mark=viscosityEnd;
    {
        // Material spans both sides of Maxwell. Count each logical material
        // interval once and put Maxwell in its own separate bucket. Diagnostic
        // conversions/sums preserve the current numerical FP flags.
        std::fenv_t environment{};const bool saved=std::fegetenv(&environment)==0;
        report.timing.phaseGeometry=std::chrono::duration<double>(middleEnd-materialStart).count()
            +std::chrono::duration<double>(phaseGeometryEnd-preparationStarted).count();
        report.timing.materialTransport=std::chrono::duration<double>(materialTransportEnd-phaseGeometryEnd).count();
        report.timing.dualMass=std::chrono::duration<double>(dualMassEnd-materialTransportEnd).count();
        report.timing.convection=std::chrono::duration<double>(convectionEnd-convectionStarted).count();
        report.timing.viscosity=std::chrono::duration<double>(viscosityEnd-viscosityStarted).count();
        report.timing.material=report.timing.phaseGeometry+report.timing.materialTransport+report.timing.dualMass
            +report.timing.convection+report.timing.viscosity;
        report.timing.maxwell=std::chrono::duration<double>(maxwellEnd-convectionEnd).count();
        if(saved&&std::fesetenv(&environment)!=0)std::terminate();
    }
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE) && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST) && !defined(__AVX__) && !defined(__AVX2__) && !defined(__AVX512F__) && !defined(__FMA__)
    joinCapillary();
    if(backgroundCapillary){
        DiagnosticRateState afterViscosity;afterViscosity.capture();
        if(workerEntered&&capillaryExit.error==overlapEntry.error
            &&capillarySameControls(overlapEntry,capillaryExit)&&capillarySameControls(overlapEntry,afterViscosity)
            &&(!flow||afterViscosity.error==overlapEntry.error)){
            // An early successful task keeps its flags private until the old
            // consumption point. Opaque providers must retain its controls and
            // errno; otherwise discard only this success and run the real HF.
            capillaryMergeFlags(afterViscosity,capillaryExit);capillary=std::move(backgroundCapillary);
            report.timing.concurrentViscosityCapillary=!flow||(steadyFlow&&!config.maxwellLoad);
        }
    }
#endif
    HydroFrozenMidpointForces forces;forces.bodyAcceleration=config.bodyAcceleration;
    if(maxwell){report.maxwellLoadReportAvailable=true;report.maxwellLoadReport=maxwell->report();}
    forces.maxwellLoad=std::move(maxwell);
    const auto serialCapillaryStarted=Wall::now();
    if(config.surfaceTension>0){
        if(!capillary){
            if(config.southEven90Boundary&&!preparedCapillaryBoundary){NativeMeshIdentity::ValidationBatch::Pause pause;
                preparedCapillaryBoundary=config.southEven90Boundary(mesh,volume.midpoint());}
            capillary=capillaryKernel();
        }
        forces.interfaceLoad=*capillary;
    }
    // These are complete kernel durations, including discarded successful
    // work. Concurrent child times overlap; the actual total remains wall time.
    report.timing.capillary=earlyBoundarySeconds+backgroundCapillarySeconds+seconds(Wall::now()-serialCapillaryStarted);
    mark=Wall::now();
    const auto frozen=[&]{try{return HydroFrozenMidpoint::prepare(mesh,dual,convection,applied,viscosity,pressure,before,std::move(forces),
        flow?HydroFrozenMidpointBoundary::prescribed(accepted.boundary(),flow->normalVelocity,flow->pressureOutlets)
            :HydroFrozenMidpointBoundary::closed(accepted.boundary()),HydroFrozenMidpointScheme::FrozenConvection,sources,clock,config.saddle,
            config.saddle.gmresOnlyParallel?&predictedEndpoint:nullptr);}
        catch(const std::exception& e){if(std::string(e.what()).find("velocity block is not positive definite")!=std::string::npos)
            throw HydroTwoPhaseStepRetry(HydroTwoPhaseRetryReason::VelocityBlock,e.what());throw;}}();
    report.timing.frozenNs=part();
    const auto& solved=frozen.solveResult().values();const auto um=midpoint(before,solved);
    const auto reassessmentStart=Wall::now();
    const auto reassess=[&]{return openMomentum?convection.applyVofOpenBoundaryBudget(mesh,um,sources,clock)
        :convection.apply(mesh,um,flow?prescribedUpwind(mesh,dual,*inlet,um):boundary,sources,clock);};
    const auto reassessed=[&]{if(steadyFlow)return reassess();
        NativeMeshIdentity::ValidationBatch::Pause pause;return reassess();}();
    const auto reassessmentEnd=Wall::now();
    if(openMomentum){report.openConvection.solvedMidpointApplySeconds=seconds(reassessmentEnd-reassessmentStart);
        report.openConvection.solvedMidpoint=reassessed.openApplyReport();}
    report.interfaceMismatch=mismatch(mesh,middle,advector,um,dual.deltaTime());report.saddle=frozen.solveResult().report();
    const Eigen::Map<const Eigen::VectorXd> u1(solved.normalVelocity().data(),layout.velocityCount()),vm(um.normalVelocity().data(),layout.velocityCount());
    const Eigen::VectorXd physicalAu=frozen.physicalVelocityMatrix()*u1,kum=frozen.viscosityMatrix()*vm;
    const auto cp=pressure.apply(mesh,solved);const auto dt=dual.deltaTime();long double physical=0,physicalScale=0,power=0,powerScale=0;
    std::vector<bool> free(layout.velocityCount());for(int u:frozen.saddle().freeVelocityDofs())free[u]=true;
    auto& energy=report.energy;long double energyScale=0;
    const bool maxwellActive=frozen.hasMaxwellLoad();
    std::optional<FV::Dyadic> maxwellPower;
    if(maxwellActive){
        maxwellPower.emplace();
        for(int axis:frozen.forces().maxwellLoad->geometry().metrics().axisVelocityDofs())
            require(vm[axis]==0,"nonzero actual radial axis velocity");
    }
    const auto* traction=flow?&frozen.saddle().pressureOutletLoads():nullptr;
    long double prescribedReactionPower=0;
    std::vector<long double> exteriorAbsolute;
    if(flow&&!openMomentum){exteriorAbsolute.resize(layout.velocityCount());const auto& segments=dual.geometry().segments();
        for(std::size_t g=0;g<segments.size();++g){const auto& segment=segments[g];
            if(segment.negativeVelocity>=0&&segment.positiveVelocity>=0)continue;
            const int owner=segment.negativeVelocity>=0?segment.negativeVelocity:segment.positiveVelocity;
            require(owner>=0,"open energy segment has no actual owner");exteriorAbsolute[owner]+=std::abs(dual.faceRates()[g]);}}
    const auto& beginVolumes=dual.beginMass().representedPhaseVolumes();
    const auto& midpointVolumes=dual.midpointMass().representedPhaseVolumes();
    const auto& endVolumes=dual.endMass().representedPhaseVolumes();
    const auto& h=frozen.frozenConvectionMatrix();
    const auto momentumRates=prepareDiagnosticMomentumRates(layout.velocityCount(),dual.beginMass().masses().data(),
        dual.endMass().masses().data(),before.normalVelocity().data(),u1.data(),dt,report.momentumRateRows);
    for(int u=0;u<layout.velocityCount();++u){const long double v=vm[u],u0=before.normalVelocity()[u],m0=dual.beginMass().masses()[u],m1=dual.endMass().masses()[u];
        // Contract the complete endpoint momentum difference before rounding.
        // Near steady flow, rounding M1*U1 and M0*U0 separately can swamp the
        // physical rate; this changes no matrix, RHS or accepted state.
        std::optional<Detail::HydroFiniteContraction> literalRate;
        if(!momentumRates){auto& momentumRate=literalRate.emplace();
            momentumRate.addQuotientProduct({m1,u1[u]},dt);
            momentumRate.addQuotientProduct({-1.L,m0,u0},dt);}
        const long double force=capillary?capillary->momentumRhsLoad()[u]:0,body=frozen.bodyLoads()[u],rate=momentumRates?momentumRates[u]:literalRate->rounded(),load=reassessed.momentumLoads()[u];
        const long double electric=maxwellActive?frozen.maxwellLoads()[u]:0;
        // Same represented operands and exact ascending-row sum as the old
        // separate power pass; round only the complete dot, then apply dt.
        if(maxwellActive)maxwellPower->add(FV::Dyadic::product({v,electric}));
        report.maximumAbsDualDefect=std::max(report.maximumAbsDualDefect,std::abs(dual.cells()[u].integratedDefect));
        report.kineticChange+=(m1*u1[u]*u1[u]-m0*u0*u0)/2;
        report.capillaryWork+=dt*force*v;
        if(free[u]){
            long double actualRhs=m0*u0/dt-.5L*frozen.beginViscousAction()[u];actualRhs+=load;actualRhs+=force;if(maxwellActive)actualRhs+=electric;actualRhs+=body;if(flow)actualRhs+=(*traction)[u];
            physical=std::hypot(physical,static_cast<long double>(physicalAu[u])+cp[u]-actualRhs);physicalScale=std::hypot(physicalScale,std::abs(static_cast<long double>(physicalAu[u]))+std::abs(static_cast<long double>(cp[u]))+std::abs(actualRhs));}
        if(flow){const long double outlet=(*traction)[u];
            const long double residual=v*(rate+kum[u]+cp[u]-load-force-body-electric-outlet);
            power+=residual;powerScale+=std::abs(v*rate)+std::abs(v*kum[u])+std::abs(v*cp[u])+std::abs(v*load)+std::abs(v*force)+std::abs(v*body)+std::abs(v*electric)+std::abs(v*outlet);
            report.pressureOutletWork+=dt*v*outlet;if(!free[u])prescribedReactionPower+=residual;
        }else if(maxwellActive){
            power+=v*(rate+kum[u]+cp[u]-load-force-body-electric);
            powerScale+=std::abs(v*rate)+std::abs(v*kum[u])+std::abs(v*cp[u])+std::abs(v*load)+std::abs(v*force)+std::abs(v*body)+std::abs(v*electric);
        }else{
        power+=v*(rate+kum[u]+cp[u]-load-force-body);
        powerScale+=std::abs(v*rate)+std::abs(v*kum[u])+std::abs(v*cp[u])+std::abs(v*load)+std::abs(v*force)+std::abs(v*body);
        }
        // Conservative VOF transport has no uR source. The retained endpoint
        // KE residual is -1/2*u_mid^2*r. H excludes M/dt and viscosity.
        const auto& cell=dual.cells()[u];const long double deltaMass=m1-m0,deltaU=static_cast<long double>(u1[u])-u0;
        constexpr long double halfSign=-.5L;
        const long double halfR=halfSign*cell.integratedDefect*v*v;
        const long double mappedHalfR=halfSign*cell.mappedPrimalIntegratedDefect*v*v;
        const long double geometryHalfR=halfSign*cell.geometryTransferChange*v*v;
        const long double spatialHalfR=halfSign*dt*cell.spatialExchangeRate*v*v;
        const long double timeMass=deltaMass*deltaU*deltaU/8;
        const long double roundedPhaseMass=.5L*(cell.massChange-deltaMass)*v*v;
        long double operatorRoundoff=(.5L*cell.outwardIntegral-static_cast<long double>(dt)*h.coeff(u,u))*v*v;
        // Open Hii has +.5*sum(exterior |q|). Remove that declared boundary
        // diagonal here; actual boundary kinetic work is accounted separately.
        if(flow&&!openMomentum)operatorRoundoff+=.5L*dt*exteriorAbsolute[u]*v*v;
        const long double bodyWork=dt*body*v,viscousWork=dt*static_cast<long double>(kum[u])*v,pressureWork=dt*static_cast<long double>(cp[u])*v;
        energy.halfRWork+=halfR;energy.mappedPrimalHalfRWork+=mappedHalfR;
        energy.geometryTransferHalfRWork+=geometryHalfR;energy.spatialExchangeHalfRWork+=spatialHalfR;
        energy.timeMassWork+=timeMass;energy.roundedPhaseMassWork+=roundedPhaseMass;energy.operatorRateRoundoffWork+=operatorRoundoff;
        // Raw rate sums are a separate rounding diagnostic, not extra work.
        energy.rawMinusCanonicalOutflowHalfWork+=.5L*(static_cast<long double>(dt)*cell.rawOutwardRate-cell.outwardIntegral)*v*v;
        energy.bodyWork+=bodyWork;energy.viscousWork+=viscousWork;energy.pressureWork+=pressureWork;
        energyScale+=std::abs((m1*u1[u]*u1[u]-m0*u0*u0)/2)+std::abs(dt*force*v)+std::abs(bodyWork)+std::abs(viscousWork)
            +std::abs(pressureWork)+std::abs(halfR)+std::abs(timeMass)+std::abs(roundedPhaseMass)+std::abs(operatorRoundoff);
        if(maxwellActive)energyScale+=std::abs(dt*electric*v);
        if(beginVolumes[u].liquid==0&&midpointVolumes[u].liquid==0&&endVolumes[u].liquid==0){
            ++energy.stableGasDofs;energy.stableGasHalfRWork+=halfR;energy.stableGasSpatialExchangeHalfRWork+=spatialHalfR;
            const long double midMass=dual.midpointMass().masses()[u],ratio=std::abs(cell.integratedDefect)/midMass;
            energy.stableGasMaxAbsIntegratedDefectOverMidMass=std::max(energy.stableGasMaxAbsIntegratedDefectOverMidMass,ratio);
            const long double speed=std::abs(static_cast<long double>(u1[u]));
            if(energy.peakGasNativeFace<0||speed>energy.peakGasAbsEndNormal){
                energy.peakGasNativeFace=layout.nativeFace(u);energy.peakGasAbsEndNormal=speed;
                energy.peakGasIntegratedDefectOverMidMass=cell.integratedDefect/midMass;
                energy.peakGasMassChangeOverMidMass=deltaMass/midMass;
                energy.peakGasOutwardIntegralOverMidMass=cell.outwardIntegral/midMass;
                energy.peakGasGeometryTransferOverMidMass=cell.geometryTransferChange/midMass;
                energy.peakGasSpatialExchangeOverMidMass=static_cast<long double>(dt)*cell.spatialExchangeRate/midMass;
            }
        }
    }
    if(maxwellActive)report.maxwellWork=static_cast<long double>(dt)*maxwellPower->rounded();
    require(!openMomentum||reassessed.openBoundaryBudget()!=nullptr,"open production convection has no boundary budget");
    if(flow){power-=prescribedReactionPower;powerScale+=std::abs(prescribedReactionPower);report.prescribedReactionWork=dt*prescribedReactionPower;
        for(int axis=0;axis<2;++axis){
            // Open production requests publish only the six consumed contractions.
            const auto* budget=reassessed.openBoundaryBudget();
            const auto& legacy=reassessed.report().components[axis];
            const auto& component=budget?budget->components[axis]:HydroKepOpenBoundaryBudget::Component{
                legacy.boundaryMomentumRate,legacy.boundaryKineticPower,legacy.boundaryPhysicalKineticPower,
                legacy.inletDissipationPower,legacy.boundaryFluxRoundingPower,legacy.boundaryDiagonalPower};
            report.boundaryMomentumRate[axis]=component.boundaryMomentumRate;
            report.boundaryKineticWork+=dt*component.boundaryKineticPower;
            report.boundaryPhysicalKineticWork+=dt*component.boundaryPhysicalKineticPower;
            report.inletDissipationWork+=dt*component.inletDissipationPower;
            report.boundaryFluxRoundingWork+=dt*component.boundaryFluxRoundingPower;
            if(openMomentum){
                const auto diagonalWork=static_cast<long double>(dt)*component.boundaryDiagonalPower;
                energy.operatorRateRoundoffWork+=diagonalWork;energyScale+=std::abs(diagonalWork);
            }
        }}
    const auto physicalDefect=relative(physical,physicalScale),powerDefect=relative(std::abs(power),powerScale);
    report.momentumPowerRelativeDefect=powerDefect;
    if(physicalDefect>1e-11||powerDefect>config.momentumPowerTolerance){
        std::ostringstream why;why.precision(std::numeric_limits<long double>::max_digits10);
        why<<"[HydroVofStep] same-owner physical momentum power gate failed physical="<<physicalDefect
            <<" power="<<powerDefect<<" power_residual="<<power<<" power_scale="<<powerScale
            <<" power_tolerance="<<config.momentumPowerTolerance;
        throw std::invalid_argument(why.str());
    }
    energy.halfRDecompositionResidual=energy.halfRWork-energy.mappedPrimalHalfRWork-energy.geometryTransferHalfRWork-energy.spatialExchangeHalfRWork;
    energy.momentumEquationWork=static_cast<long double>(dt)*power;
    long double budget=report.capillaryWork+energy.bodyWork-energy.viscousWork-energy.pressureWork
        +energy.halfRWork+energy.timeMassWork+energy.roundedPhaseMassWork+energy.operatorRateRoundoffWork;
    if(maxwellActive)budget+=report.maxwellWork;
    if(flow){budget+=report.pressureOutletWork+report.prescribedReactionWork-report.boundaryKineticWork;
        energyScale+=std::abs(report.pressureOutletWork)+std::abs(report.prescribedReactionWork)+std::abs(report.boundaryKineticWork);}
    energy.kineticBudgetResidual=report.kineticChange-budget;
    energy.closureResidual=energy.kineticBudgetResidual-energy.momentumEquationWork;
    const auto unavailable=std::numeric_limits<long double>::quiet_NaN();
    energy.closureRelativeResidual=energyScale>0?std::abs(energy.closureResidual)/energyScale:(energy.closureResidual==0?0:unavailable);
    if(energy.stableGasDofs==0){
        energy.stableGasMaxAbsIntegratedDefectOverMidMass=unavailable;energy.peakGasAbsEndNormal=unavailable;
        energy.peakGasIntegratedDefectOverMidMass=unavailable;energy.peakGasMassChangeOverMidMass=unavailable;
        energy.peakGasOutwardIntegralOverMidMass=unavailable;energy.peakGasGeometryTransferOverMidMass=unavailable;
        energy.peakGasSpatialExchangeOverMidMass=unavailable;
    }
    report.timing.diagnostics=part();report.timing.total=seconds(Wall::now()-start);
    return {volume,frozen,um,std::move(capillary),report};
}
}
struct HydroVofStep::Data {
    Data(HydroTwoPhaseState a,HydroTwoPhaseState c,std::vector<Stage> s,HydroVofStepReport r)
        :accepted(std::move(a)),candidate(std::move(c)),stages(std::move(s)),report(std::move(r)){}
    HydroTwoPhaseState accepted,candidate;std::vector<Stage> stages;HydroVofStepReport report;
};
HydroVofStep::HydroVofStep(std::shared_ptr<const Data> d):data_(std::move(d)){}
HydroVofStabilityRates HydroVofStep::stabilityRates(const QuadtreeMesh& mesh,
    const HydroTwoPhaseState& accepted,const HydroVofAcConfig& config)
{return stabilityRatesImpl(mesh,accepted,config,false);}
HydroVofStabilityRates HydroVofStep::stabilityRatesPrescribedFlow(const QuadtreeMesh& mesh,
    const HydroTwoPhaseState& accepted,const HydroVofAcConfig& config)
{return stabilityRatesImpl(mesh,accepted,config,true);}
HydroVofStabilityRates HydroVofStep::stabilityRatesImpl(const QuadtreeMesh& mesh,
    const HydroTwoPhaseState& accepted,const HydroVofAcConfig& config,bool prescribed)
{
    // Only an actual retained typed-flow Boundary excludes opaque laws.
    if(accepted.vofTransportBoundary())accepted.validateFor(mesh);
    else {NativeMeshIdentity::ValidationBatch::Pause pause;accepted.validateFor(mesh);}
    NativeMeshIdentity::ValidationBatch validationBatch(mesh);
    require(accepted.model()==HydroPhaseGeometry::Model::GeometricVof&&!accepted.prescribedFlow(),"stability rates require a genuine alpha-only VOF state");
    require(bool(accepted.vofPrescribedFlow())==prescribed,"stability rates physical VOF boundary capability differs from explicit entry");
    require(!prescribed||!config.transportProvider,"explicit transport provider supports closed VOF flow only");
    const auto& source=accepted.phaseGeometry().vofPublication();
    require(source.geometryKind()==I::GeometricVofGeometryKind::LocalPlic,"stability rates require actual local PLIC publication");
    require(config.liquidDensity>0&&std::isfinite(config.liquidDensity)&&config.gasDensity>0&&std::isfinite(config.gasDensity),"invalid phase density");
    require(config.surfaceTension>=0&&std::isfinite(config.surfaceTension),"invalid surface tension");
    const auto& values=accepted.values();const auto& layout=values.layout();const auto& native=source.domain();
    require(layout.domain().axisymmetric()==native.axisymmetric()&&layout.domain().meshIdentity()==native.meshIdentity(),"stability rates have a different native mesh/measure");
    require(static_cast<std::size_t>(layout.pressureCount())==native.cellIndexByRow().size()
        &&layout.domain().cellIndexByRow()==native.cellIndexByRow()
        &&layout.domain().cellRowByIndex()==native.cellRowByIndex(),"stability rates require actual selected VOF leaf coverage");
    require(layout.domain().nativeFaceCount()==native.nativeFaceCount(),"stability rates require the same actual native face layout");
    std::vector<double> normal(static_cast<std::size_t>(native.nativeFaceCount()));
    for(int row=0;row<native.nativeFaceCount();++row){const auto f=native.face(mesh,row);const auto u=layout.velocityDof(row);
        if(f.negativeCell<0&&f.positiveCell<0){require(u<0,"excluded stability face has a Hydro velocity DOF");continue;}
        require(f.nativeRow==row&&u>=0&&layout.nativeFace(u)==row,"invalid stability native normal mapping");
        const auto axis=f.horizontal?HydroDofLayout::NormalAxis::Second:HydroDofLayout::NormalAxis::First;
        require(layout.normalAxis(u)==axis,"stability native normal axis differs");
        normal[row]=values.normalVelocity()[u];
    }
    HydroVofStabilityRates result;
    if(prescribed){
        if(const auto* boundary=accepted.vofTransportBoundary())
            result.transport=I::GeometricVofTransport::localPlicBoundaryCourantRates(mesh,source,*boundary,normal);
        else {const auto preparedBoundary=I::GeometricVofAdvectorSlab::Boundary::prepare(mesh,source,accepted.vofPrescribedFlow()->rulesByNativeRow);
            result.transport=I::GeometricVofTransport::localPlicBoundaryCourantRates(mesh,source,preparedBoundary,normal);}
    }
    else result.transport=I::GeometricVofTransport::localPlicCourantRates(mesh,source,normal);
    result.flowRate=std::max({result.transport.maximumLocalColumn,result.transport.maximumIncoming,result.transport.maximumOutflow});
    const auto scale=capillaryScale(source,config);result.capillaryAvailable=scale.available;
    if(scale.available){result.interfaceMinimumSpacing=scale.minimumSpacing;
        result.capillaryRate=std::exp(static_cast<long double>(scale.logRate));
        require(std::isfinite(result.capillaryRate)&&result.capillaryRate>0,"unrepresentable positive interface capillary rate");}
    return result;
}
HydroVofStep HydroVofStep::prepare(const QuadtreeMesh& mesh,const HydroTwoPhaseState& accepted,Clock end,const HydroVofAcConfig& config)
{return prepareImpl(mesh,accepted,end,config,false);}
HydroVofStep HydroVofStep::preparePrescribedFlow(const QuadtreeMesh& mesh,const HydroTwoPhaseState& accepted,Clock end,const HydroVofAcConfig& config)
{return prepareImpl(mesh,accepted,end,config,true);}
HydroVofStep HydroVofStep::prepareImpl(const QuadtreeMesh& mesh,const HydroTwoPhaseState& accepted,Clock end,const HydroVofAcConfig& config,bool prescribed)
{
    // Only repeated exact mesh scans share this callback-free preparation
    // receipt. Every opaque callback below suspends and clears the whole chain.
    NativeMeshIdentity::ValidationBatch validationBatch(mesh);
    if(accepted.vofTransportBoundary())accepted.validateFor(mesh);
    else {NativeMeshIdentity::ValidationBatch::Pause pause;accepted.validateFor(mesh);}
    require(accepted.model()==HydroPhaseGeometry::Model::GeometricVof&&!accepted.prescribedFlow(),"requires a genuine alpha-only VOF state");
    require(bool(accepted.vofPrescribedFlow())==prescribed,"VOF physical boundary capability differs from explicit step entry");
    require(!prescribed||!config.transportProvider,"explicit transport provider supports closed VOF flow only");
    const auto& source=accepted.phaseGeometry().vofPublication();require(source.geometryKind()==I::GeometricVofGeometryKind::LocalPlic,"requires actual local PLIC publication");
    require(config.momentumDensityRepresentation==HydroMomentumMass::DensityRepresentation::CellAverageProjection,
        "conservative VOF momentum requires inventory cellAverageProjection mass; sharpGeometry is a reference-only metric");
    require(config.liquidDensity>0&&std::isfinite(config.liquidDensity)&&config.gasDensity>0&&std::isfinite(config.gasDensity),"invalid phase density");
    require(config.liquidViscosity>=0&&std::isfinite(config.liquidViscosity)
        &&config.gasViscosity>=0&&std::isfinite(config.gasViscosity),"invalid phase viscosity");
    require(config.surfaceTension>=0&&std::isfinite(config.surfaceTension),"invalid surface tension");
    require(std::isfinite(config.bodyAcceleration.x)&&std::isfinite(config.bodyAcceleration.y),"invalid body acceleration");
    require(!source.domain().axisymmetric()||config.bodyAcceleration.x==0,"RZ body acceleration must be axial");
    require(config.scheme==HydroTwoPhaseTimeScheme::AdaptiveLinearMidpoint||config.scheme==HydroTwoPhaseTimeScheme::ExtrapolatedMidpoint,"unsupported pure VOF time scheme");
    require(config.interfaceMismatchTolerance>0&&std::isfinite(config.interfaceMismatchTolerance)&&config.maximumCorrectionRelative>0&&std::isfinite(config.maximumCorrectionRelative)
        &&config.correctionReferenceVelocity>0&&std::isfinite(config.correctionReferenceVelocity),"invalid AC correction controls");
    require(config.capillaryCorrectionCourant>0&&config.maximumCapillaryCourant>=config.capillaryCorrectionCourant&&std::isfinite(config.maximumCapillaryCourant)
        &&config.momentumPowerTolerance>0&&config.momentumPowerTolerance<1,"invalid capillary or power controls");
    if(config.requireMaxwellCorrection)require(bool(config.maxwellLoad)&&std::isfinite(config.maximumMaxwellLoadChange)
        &&config.maximumMaxwellLoadChange>0&&config.maximumMaxwellLoadChange<=.25&&std::isfinite(config.maximumMaxwellWorkChange)
        &&config.maximumMaxwellWorkChange>0&&config.maximumMaxwellWorkChange<=.25,"invalid explicit electric correction controls");
    const auto old=accepted.timeLayer();const double dt=end.end-end.begin;
    require(end.position==Clock::Position::End&&old.step<std::numeric_limits<std::int64_t>::max()&&end.step==old.step+1&&end.begin==old.time()&&std::isfinite(dt)&&dt>0,"requires contiguous complete End clock");
    end.time();const auto& before=accepted.values();const auto& geometry=accepted.momentumGeometry();
    HydroVofStepReport report;report.deltaTime=dt;
    const auto* history=accepted.velocityHistory();const bool useHistory=history&&!config.resetVelocityHistory;
    report.usedVelocityHistory=useHistory;report.startupCorrection=!useHistory;
    std::vector<double> predicted=before.normalVelocity();
    if(useHistory){require(history->model==HydroPhaseGeometry::Model::GeometricVof&&history->vofSource.value()!=0,"different VOF velocity history");
        const double previousDt=old.end-old.begin;require(previousDt>0&&std::isfinite(previousDt),"invalid previous time interval");
        for(std::size_t u=0;u<predicted.size();++u){predicted[u]+=dt/(2*previousDt)*(before.normalVelocity()[u]-history->values.normalVelocity()[u]);require(std::isfinite(predicted[u]),"nonfinite midpoint extrapolation");}}
    bool forced=false;const auto scale=capillaryScale(source,config);
    if(scale.available){
        report.capillaryCourant=std::exp(std::log(dt)+scale.logRate);
        if(!std::isfinite(report.capillaryCourant)||report.capillaryCourant>config.maximumCapillaryCourant)throw HydroTwoPhaseStepRetry(HydroTwoPhaseRetryReason::Capillary,"[HydroVofStep] interface-leaf capillary Courant exceeds configured bound");
        forced=report.capillaryCourant>=config.capillaryCorrectionCourant;}
    report.capillaryForcedCorrection=forced;
    const auto beginMass=[&]{const auto* retained=accepted.retainedMass();if(retained&&retained->liquidDensity()==config.liquidDensity&&retained->gasDensity()==config.gasDensity&&retained->representation()==config.momentumDensityRepresentation){retained->validateVofFor(mesh,source.identity(),old);return *retained;}
        return HydroMomentumMass::buildCellAverageProjection(mesh,geometry,accepted.phaseGeometry(),config.liquidDensity,config.gasDensity);}();
    const auto pressure=HydroPressureLoad::build(geometry.metrics().divergence());std::vector<Stage> stages;stages.reserve(2);
    // Retain only the first AC values actually consumed by the correction.
    // A PhaseViscosity copy shares its original raw strain/geometry owner;
    // none of these receipts retains the first Frozen/Saddle/LU.
    std::optional<HydroNativeValues> firstMidpoint;
    std::optional<HydroPhaseViscosity> firstViscosity;
    bool firstHasMaxwellLoad=false;std::vector<double> firstMaxwellLoads;
    HydroNativeValues advector(before.layout(),std::move(predicted),before.pressure());
    for(int iteration=0;iteration<2;++iteration){
        // The retained first viscosity owns canonical raw P2 geometry for this attempt.
        // No viscosity coefficients, K or phase clock are reused or retagged.
        const auto* previousViscosity=firstViscosity?&*firstViscosity:nullptr;
        stages.push_back(stage(mesh,accepted,end,config,beginMass,advector,pressure,prescribed,previousViscosity));const auto& current=stages.back();
        report.stages.push_back(current.report);++report.linearSolves;report.finalInterfaceMismatch=current.report.interfaceMismatch;
        if(iteration==0)report.interfaceMismatch=current.report.interfaceMismatch;
        else{report.corrected=true;report.correctionRelative=correctionRelative(*firstMidpoint,current.solvedMidpoint,current.frozen.material().midpointMass(),config.correctionReferenceVelocity);
            if(report.correctionRelative>config.maximumCorrectionRelative)throw HydroTwoPhaseStepRetry(HydroTwoPhaseRetryReason::Correction,"[HydroVofStep] fixed second-stage correction exceeds bound");
            if(config.requireMaxwellCorrection){
                require(firstHasMaxwellLoad&&current.frozen.hasMaxwellLoad(),"electric correction lost an actual stage load");
                const auto& f0=firstMaxwellLoads;const auto& f1=current.frozen.maxwellLoads();const auto& u=current.solvedMidpoint.normalVelocity();
                require(f0.size()==f1.size()&&u.size()==f1.size(),"electric correction native dimensions differ");
                long double forceError=0,forceScale=0,workError=0,workScale=0;
                for(std::size_t k=0;k<f0.size();++k){const long double a=f0[k],b=f1[k],v=u[k];
                    forceError=std::hypot(forceError,b-a);forceScale=std::hypot(forceScale,std::abs(a)+std::abs(b));
                    workError+=(b-a)*v;workScale+=(std::abs(a)+std::abs(b))*std::abs(v);}
                report.maxwellLoadChange=relative(forceError,forceScale);report.maxwellWorkChange=relative(std::abs(workError),workScale);
                report.maxwellCorrectionChecked=true;
                if(report.maxwellLoadChange>config.maximumMaxwellLoadChange||report.maxwellWorkChange>config.maximumMaxwellWorkChange)
                    throw HydroTwoPhaseStepRetry(HydroTwoPhaseRetryReason::Correction,"[HydroVofStep] same-velocity electric load/work stage change exceeds bound");
            }}
        const bool correct=iteration==0&&(config.requireMaxwellCorrection||!useHistory||forced||(config.scheme==HydroTwoPhaseTimeScheme::AdaptiveLinearMidpoint&&current.report.interfaceMismatch>config.interfaceMismatchTolerance));
        if(!correct)break;advector=current.solvedMidpoint;
        firstMidpoint.emplace(current.solvedMidpoint);
        firstViscosity.emplace(current.frozen.viscosity());
        firstHasMaxwellLoad=current.frozen.hasMaxwellLoad();
        if(firstHasMaxwellLoad)firstMaxwellLoads=current.frozen.maxwellLoads();
        // All current references end here. Release the first complete Frozen
        // (including its LU) before preparing the actual second AC stage.
        stages.clear();
    }
    const auto& final=stages.back();auto candidate=[&]{
        // Accepted-state construction can invoke the physical endpoint law.
        NativeMeshIdentity::ValidationBatch::Pause pause;
        return HydroTwoPhaseState::createVofAccepted(mesh,final.frozen.material().material().endGeometry(),
            final.frozen.solveResult().values(),accepted,final.frozen.material().endMass());}();
    return HydroVofStep(std::make_shared<Data>(accepted,std::move(candidate),std::move(stages),std::move(report)));
}
const HydroTwoPhaseState& HydroVofStep::beginState()const noexcept{return data_->accepted;}
const HydroTwoPhaseState& HydroVofStep::candidateState()const noexcept{return data_->candidate;}
const HydroFrozenMidpoint& HydroVofStep::frozenMidpoint()const noexcept{return data_->stages.back().frozen;}
const HydroNativeValues& HydroVofStep::solvedMidpointValues()const noexcept{return data_->stages.back().solvedMidpoint;}
const I::GeometricVofTransport& HydroVofStep::transport()const noexcept{return data_->stages.back().volume;}
bool HydroVofStep::hasCapillaryLoad()const noexcept{return data_->stages.back().capillary.has_value();}
const HydroBalancedCsfLoad& HydroVofStep::capillaryLoad()const{require(hasCapillaryLoad(),"zero sigma has no capillary load");return *data_->stages.back().capillary;}
const HydroVofStepReport& HydroVofStep::report()const noexcept{return data_->report;}
void HydroVofStep::commit(const QuadtreeMesh& mesh,HydroTwoPhaseState& accepted)const
{require(accepted.sameIdentity(beginState()),"stale VOF step cannot replace accepted state");
    NativeMeshIdentity::ValidationBatch::Pause pause;
    accepted.validateFor(mesh);candidateState().validateFor(mesh);accepted=data_->candidate;}
}
