#include "Modules/Interface/GeometricVofTransport.h"
#include "Modules/Interface/VofSmoothGeometry.h"
#include "Modules/Interface/GeometricVofPhaseBudget.h"
#include "Core/Geometry/PlanarFaceTimeIntegral.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/CartesianMeasure.h"
#include "Core/Geometry/Detail/FiniteVolumeClosure.h"
#include <Eigen/QR>
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <atomic>
#include <cfenv>
#include <cerrno>
#include <cstring>
#include <cstdint>
#include <exception>
#include <type_traits>
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
#include <omp.h>
#include <float.h>
#include <xmmintrin.h>
#endif

namespace LAMP::Interface
{
namespace
{
namespace Cut=Geometry::PlanarCellCut;
namespace FaceTime=Geometry::PlanarFaceTime;
namespace FiniteVolume=Geometry::Detail::FiniteVolume;
using LocalDyadic=FiniteVolume::Dyadic;
using Point=Cut::Point;
using AnchoredPoint=Cut::AnchoredPoint;
using Clock=GeometricVofTimeLayer;
using Poly=std::vector<AnchoredPoint>;
constexpr long double arithmetic=512*std::numeric_limits<long double>::epsilon();
void require(bool ok,const char* reason)
{if(!ok)throw std::invalid_argument(std::string("[Interface/GeometricVofTransport] ")+reason);}
void require(bool ok,const std::string& reason)
{if(!ok)throw std::invalid_argument("[Interface/GeometricVofTransport] "+reason);}
Point nativePoint(const GeometricVofPublication& source,std::size_t vertex)
{
    const auto& sites=source.anchoredVertices();require(vertex<sites.size(),"native face site outside publication");
    const auto& site=sites[vertex];
    require(site.offset.x==0&&site.offset.y==0&&std::isfinite(site.anchor.x)&&std::isfinite(site.anchor.y),
        "native face endpoint must retain its original zero-offset Point");
    return site.anchor;
}
struct Sum{long double value=0,error=0;void add(long double x){const auto y=x-error,t=value+y;error=(t-value)-y;value=t;}};
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
#pragma float_control(precise,on,push)
#pragma float_control(except,on)
#pragma fenv_access(on)
// Optional execution only. The literal face loop remains the failure replay.
struct LocalFaceFloatingState
{
    std::fenv_t native{};unsigned csr=0;bool valid=false;
    bool capture() noexcept {csr=_mm_getcsr();valid=std::fegetenv(&native)==0;return valid;}
    void restore() const noexcept
    {if(valid&&std::fesetenv(&native)!=0)std::terminate();_mm_setcsr(csr);}
    bool enterFresh(const LocalFaceFloatingState& from) const noexcept
    {auto clean=from.native;clean._Fe_stat=0;if(std::fesetenv(&clean)!=0)return false;
        _mm_setcsr(from.csr&~0x3fu);return true;}
};
bool localFaceSameControl(const LocalFaceFloatingState& a,const LocalFaceFloatingState& b) noexcept
{
    auto x=a.native,y=b.native;x._Fe_stat=0;y._Fe_stat=0;
    return std::memcmp(&x,&y,sizeof x)==0&&(a.csr&~0x3fu)==(b.csr&~0x3fu);
}
bool localFaceQualified(LocalFaceFloatingState& caller) noexcept
{
#if !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST)
    if(!caller.capture()||omp_in_parallel()||omp_get_max_threads()<2)return false;
    unsigned control=0;
    if(_controlfp_s(&control,0,0)!=0||(control&_MCW_EM)!=_MCW_EM
        ||(control&_MCW_RC)!=_RC_NEAR||std::fegetround()!=FE_TONEAREST
        ||(caller.csr&0x1f80u)!=0x1f80u||(caller.csr&0xe040u)!=0)return false;
#ifdef _EM_AMBIGUOUS
    if(control&_EM_AMBIGUOUS)return false;
#endif
    return sizeof(double)==8&&sizeof(long double)==8&&std::numeric_limits<double>::is_iec559
        &&std::numeric_limits<double>::digits==53&&std::numeric_limits<long double>::digits==53
        &&std::is_nothrow_swappable_v<FiniteVolume::Dyadic>&&std::is_nothrow_move_assignable_v<FiniteVolume::Dyadic>;
#else
    (void)caller;return false;
#endif
}
struct LocalFaceCleanup
{
    LocalFaceFloatingState finalState;int error=0;bool armed=false;
    void arm() noexcept
    {error=errno;armed=finalState.capture();errno=error;if(!armed)std::terminate();}
    ~LocalFaceCleanup() noexcept {if(armed){finalState.restore();errno=error;}}
};
struct LocalFaceWorkerFrame
{
    LocalFaceFloatingState before,after;int error=0;bool captured=false,ready=false;
};
struct LocalFaceIncidence {std::uint32_t row=0,side=0;};
struct LocalFaceObservation {long double correction=0;bool limited=false;};
struct LocalFaceRoute
{
    // First member dies last, after optional-array destruction.
    LocalFaceCleanup cleanup;
    std::unique_ptr<std::uint32_t[]> offsets;
    std::unique_ptr<LocalFaceIncidence[]> incidence;
    std::unique_ptr<LocalFaceObservation[]> observations;
    std::size_t cells=0,faces=0,incidences=0,bytes=0;
    std::array<std::size_t,2> orientedFaces{};int team=1;
    ~LocalFaceRoute() noexcept {cleanup.arm();}
};
struct LocalFaceOptionalFailure {};
struct LocalFaceBudget
{
    static constexpr std::size_t limit=64u*1024u*1024u;std::size_t used=0;
    void add(std::size_t count,std::size_t width)
    {if(width==0||used>limit||count>(limit-used)/width)throw std::length_error("optional local face scratch cap");used+=count*width;}
};
std::unique_ptr<LocalFaceRoute> prepareLocalFaceRoute(const std::vector<NativeLeafFace>& native,
    std::size_t n,bool enabled,GeometricVofLocalFaceExecutionReport& report)
{
    if(!enabled)return {};
    const int error=errno;LocalFaceFloatingState entry;const bool eligible=localFaceQualified(entry);
    report.requested=!omp_in_parallel()&&omp_get_max_threads()>1;
    if(entry.valid)entry.restore();errno=error;
    if(!eligible)return {};
    try{
        if(n>=std::numeric_limits<std::uint32_t>::max()||native.size()>std::numeric_limits<std::uint32_t>::max())
            throw std::length_error("optional local face index width");
        std::size_t incidences=0;
        for(std::size_t row=0;row<native.size();++row){const auto& f=native[row];
            if(f.nativeRow<0||static_cast<std::size_t>(f.nativeRow)!=row)throw std::length_error("optional local face row");
            for(int cell:{f.negativeCell,f.positiveCell})if(cell>=0){
                if(static_cast<std::size_t>(cell)>=n||incidences==std::numeric_limits<std::uint32_t>::max())
                    throw std::length_error("optional local face incidence");++incidences;}}
        const auto workers=std::min<std::size_t>(static_cast<std::size_t>(omp_get_max_threads()),
            std::max<std::size_t>(1,std::max(n,native.size())));
        LocalFaceBudget budget;budget.add(1,sizeof(LocalFaceRoute));
        budget.add(n+1,sizeof(std::uint32_t));budget.add(n,sizeof(std::uint32_t));
        budget.add(incidences,sizeof(LocalFaceIncidence));budget.add(native.size(),sizeof(LocalFaceObservation));
        // At most uint32_t native terms per cell: binary64 exponent span +
        // 32 carry bits is < 2132 bits. 8KiB/worker conservatively covers the
        // transient UInt/add objects, worker frame and fixed scalar scratch.
        // LocalDyadic is the original Core storage; no extra filter payload.
        // Keep the doubled worker reserve for transient UInt work unchanged.
        budget.add(workers,sizeof(LocalFaceWorkerFrame)+16384u);budget.add(8,64);
        budget.add(2,sizeof(FiniteVolume::Dyadic)+2048u);budget.add(1,sizeof(GeometricVofTransportReport));
        auto out=std::make_unique<LocalFaceRoute>();out->cells=n;out->faces=native.size();out->incidences=incidences;
        out->bytes=budget.used;out->team=static_cast<int>(workers);
        out->offsets=std::make_unique<std::uint32_t[]>(n+1);
        out->incidence=std::make_unique<LocalFaceIncidence[]>(incidences);
        out->observations=std::make_unique<LocalFaceObservation[]>(native.size());
        auto cursor=std::make_unique<std::uint32_t[]>(n);
        for(const auto& f:native){++out->orientedFaces[f.horizontal?1:0];
            for(int cell:{f.negativeCell,f.positiveCell})if(cell>=0)++out->offsets[static_cast<std::size_t>(cell)+1];}
        for(std::size_t i=0;i<n;++i){out->offsets[i+1]+=out->offsets[i];cursor[i]=out->offsets[i];}
        for(const auto& f:native){
            if(f.negativeCell>=0)out->incidence[cursor[f.negativeCell]++]={static_cast<std::uint32_t>(f.nativeRow),0};
            if(f.positiveCell>=0)out->incidence[cursor[f.positiveCell]++]={static_cast<std::uint32_t>(f.nativeRow),1};}
        cursor.reset();entry.restore();errno=error;
        report.eligible=true;report.nativeFaces=out->faces;report.incidences=out->incidences;report.scratchBytes=out->bytes;
        return out;
    }catch(...){entry.restore();errno=error;report.fallback=true;return {};}
}
template<class Kernel>
bool localFaceTeam(int team,Kernel&& kernel,unsigned& actualWorkers)
{
    LocalFaceFloatingState prepared;const int parentErrno=errno;if(!prepared.capture())return false;
    std::atomic<int> failed{0};unsigned long nativeFlags=0;unsigned csrFlags=0;
#pragma omp parallel num_threads(team)
    {
        LocalFaceWorkerFrame worker;worker.error=errno;worker.captured=worker.before.capture();
        worker.ready=worker.captured&&worker.before.enterFresh(prepared);
        if(!worker.ready)failed.store(1,std::memory_order_relaxed);else errno=parentErrno;
#pragma omp single
        {actualWorkers=std::max(actualWorkers,static_cast<unsigned>(omp_get_num_threads()));}
        kernel(worker.ready,failed);
        if(worker.ready){const int numericalErrno=errno;const bool got=worker.after.capture();
            if(!got||numericalErrno!=parentErrno||!localFaceSameControl(worker.after,prepared))failed.store(1,std::memory_order_relaxed);
#pragma omp critical(lamp_local_face_floating_flags)
            {if(got){nativeFlags|=worker.after.native._Fe_stat;csrFlags|=worker.after.csr&0x3fu;}}}
        if(worker.captured)worker.before.restore();errno=worker.error;
    }
    if(failed.load(std::memory_order_relaxed))return false;
    auto merged=prepared.native;merged._Fe_stat|=nativeFlags;if(std::fesetenv(&merged)!=0)std::terminate();
    _mm_setcsr(prepared.csr|csrFlags);errno=parentErrno;return true;
}
template<class Kernel>
bool localFaceRows(std::size_t count,int team,Kernel&& kernel,unsigned& actualWorkers)
{
    return localFaceTeam(team,[&](bool ready,std::atomic<int>& failed){
#pragma omp for schedule(static) nowait
        for(std::ptrdiff_t row=0;row<static_cast<std::ptrdiff_t>(count);++row){
            if(!ready)continue;
            try{kernel(static_cast<std::size_t>(row));}catch(...){failed.store(1,std::memory_order_relaxed);}
        }
    },actualWorkers);
}
struct LocalFaceInputs
{
    const std::vector<NativeLeafFace>& native;
    const std::vector<long double>& total;
    const std::vector<bool>& counterflow;
    const std::vector<std::array<long double,2>>& counterGross;
    const std::vector<std::array<long double,2>>& gross;
    const std::map<int,std::array<long double,2>>& wholePhaseFlux;
    const std::vector<int>& counterPhase;
    const std::vector<long double>& limiter;
    const std::vector<unsigned>& limitedPhase;
    const std::vector<int>& exactLimitedFace;
    const std::vector<std::array<long double,2>>& stocks;
    const GeometricVofAdvectorSlab::Boundary* boundary;
    const GeometricVofPublication& source;
    std::vector<GeometricVofTransportFace>& faces;
    bool horizontal=false,half=false;long double rhoL=0,rhoG=0;
};
struct LocalFaceAccumulators
{
    std::vector<std::array<Sum,2>>& net;
    std::vector<Sum>& directionalTotal;
    std::vector<std::array<LocalDyadic,2>>& exactNet;
    std::vector<std::array<LocalDyadic,2>>& actualOut;
    std::vector<std::array<LocalDyadic,2>>& actualIn;
    std::vector<LocalDyadic>& exactTotal;
};
void localFaceCompute(const LocalFaceInputs& in,std::size_t row,LocalFaceObservation& observation)
{
    const auto& f=in.native[row];if(f.horizontal!=in.horizontal)return;
    observation={};const auto q=in.total[row];auto& flux=in.half?in.faces[row].half:in.faces[row].full;
    flux.total=q;flux.positiveTotal=std::max(q,0.L);flux.negativeTotal=std::max(-q,0.L);
    if(in.counterflow[row]){
        const int inside=f.negativeCell>=0?f.negativeCell:f.positiveCell;
        if(in.counterPhase[inside]==2){
            const bool insideNegative=f.negativeCell>=0;
            const auto positive=in.counterGross[row][0],negative=in.counterGross[row][1];
            const auto outTotal=insideNegative?positive:negative,inTotal=insideNegative?negative:positive;
            auto outgoingPhases=in.gross[row];
            const auto mode=in.boundary->rules()[row];
            const unsigned incomingPhase=mode==GeometricVofAdvectorSlab::Boundary::Mode::PureLiquidReservoir?0:1;
            std::array<long double,2> incomingPhases{};incomingPhases[incomingPhase]=inTotal;
            flux.rawLiquid=insideNegative?outgoingPhases[0]-incomingPhases[0]:incomingPhases[0]-outgoingPhases[0];
            if(in.limiter[inside]<1){const auto t=in.limitedPhase[inside];const auto original=outgoingPhases[t];
                outgoingPhases[t]*=in.limiter[inside];outgoingPhases[1-t]=outTotal-outgoingPhases[t];
                observation.correction=std::abs(outgoingPhases[t]-original)/in.source.cells()[inside].capacity;observation.limited=true;}
            flux.limiter=in.limiter[inside];flux.positiveTotal=positive;flux.negativeTotal=negative;
            const auto& positivePhases=insideNegative?outgoingPhases:incomingPhases;
            const auto& negativePhases=insideNegative?incomingPhases:outgoingPhases;
            flux.positiveLiquid=positivePhases[0];flux.positiveGas=positivePhases[1];
            flux.negativeLiquid=negativePhases[0];flux.negativeGas=negativePhases[1];
            flux.liquid=positivePhases[0]-negativePhases[0];flux.gas=positivePhases[1]-negativePhases[1];
            flux.mass=in.rhoL*flux.liquid+in.rhoG*flux.gas;
            require(std::isfinite(flux.mass),"unrepresentable counterflow P0 material net");
            for(unsigned p=0;p<2;++p)
                require(std::isfinite(outgoingPhases[p])&&outgoingPhases[p]>=0,"invalid counterflow P0 outgoing phase");
            return;
        }
        const unsigned phase=static_cast<unsigned>(in.counterPhase[inside]);
        require(phase<2&&in.limiter[inside]==1,"counterflow pure donor requires unmodified gross availability");
        const auto positive=in.counterGross[row][0],negative=in.counterGross[row][1];
        flux.positiveTotal=positive;flux.negativeTotal=negative;
        if(phase==0){flux.positiveLiquid=positive;flux.negativeLiquid=negative;flux.liquid=q;flux.rawLiquid=q;}
        else{flux.positiveGas=positive;flux.negativeGas=negative;flux.gas=q;}
        flux.mass=in.rhoL*flux.liquid+in.rhoG*flux.gas;
        require(std::isfinite(flux.mass),"unrepresentable counterflow material net");return;
    }
    if(q==0)return;const auto donor=q>0?f.negativeCell:f.positiveCell;
    const auto t=donor>=0?in.limitedPhase[donor]:0;auto phases=in.gross[row];flux.rawLiquid=(q>0?1:-1)*phases[0];
    if(const auto receipt=in.wholePhaseFlux.find(f.nativeRow);receipt!=in.wholePhaseFlux.end()){
        const auto sign=q>0?1.L:-1.L;
        flux.exhaustionLiquidCorrection=sign*(receipt->second[0]-phases[0]);
        flux.exhaustionGasCorrection=sign*(receipt->second[1]-phases[1]);phases=receipt->second;
    }
    auto faceLimiter=donor>=0?in.limiter[donor]:1;
    if(donor>=0&&in.limiter[donor]<1){const auto original=phases[t];
        if(in.exactLimitedFace[donor]==f.nativeRow)phases[t]=in.stocks[donor][t];
        else {require(in.exactLimitedFace[donor]<0||original==0,"unique phase outflow has another nonzero face");phases[t]*=in.limiter[donor];}
        phases[1-t]=std::abs(q)-phases[t];
        if(in.exactLimitedFace[donor]>=0)faceLimiter=original==0?1:phases[t]/original;
        observation.correction=std::abs(phases[t]-original)/in.source.cells()[donor].capacity;observation.limited=true;}
    flux.limiter=faceLimiter;flux.liquid=(q>0?1:-1)*phases[0];flux.gas=(q>0?1:-1)*phases[1];flux.mass=in.rhoL*flux.liquid+in.rhoG*flux.gas;
    require(std::isfinite(flux.mass),"unrepresentable local material flux");
    if(q>0){flux.positiveLiquid=phases[0];flux.positiveGas=phases[1];}
    else {flux.negativeLiquid=phases[0];flux.negativeGas=phases[1];}
    for(unsigned phase=0;phase<2;++phase)
        require(std::isfinite(phases[phase])&&phases[phase]>=0,"invalid shared local phase flux");
}
void localFaceGather(const LocalFaceInputs& in,const LocalFaceRoute& route,LocalFaceAccumulators& out,std::size_t cell)
{
    for(auto at=route.offsets[cell];at<route.offsets[cell+1];++at){const auto entry=route.incidence[at];
        const auto& f=in.native[entry.row];if(f.horizontal!=in.horizontal)continue;
        const auto q=in.total[entry.row];const auto& flux=in.half?in.faces[entry.row].half:in.faces[entry.row].full;
        if(in.counterflow[entry.row]){
            const int inside=f.negativeCell>=0?f.negativeCell:f.positiveCell;
            if(in.counterPhase[inside]==2){
                for(unsigned phase=0;phase<2;++phase){
                    const auto positive=phase==0?flux.positiveLiquid:flux.positiveGas;
                    const auto negative=phase==0?flux.negativeLiquid:flux.negativeGas;
                    out.actualOut[cell][phase].add(f.negativeCell>=0?positive:negative);
                    out.actualIn[cell][phase].add(f.negativeCell>=0?negative:positive);
                    const auto signedPhase=phase==0?flux.liquid:flux.gas;
                    if(entry.side==0){out.net[cell][phase].add(signedPhase);out.exactNet[cell][phase].add(signedPhase);}
                    else {out.net[cell][phase].add(-signedPhase);out.exactNet[cell][phase].add(signedPhase,-1);}
                }
            }else{
                const unsigned phase=static_cast<unsigned>(in.counterPhase[inside]);
                out.actualOut[cell][phase].add(f.negativeCell>=0?flux.positiveTotal:flux.negativeTotal);
                out.actualIn[cell][phase].add(f.negativeCell>=0?flux.negativeTotal:flux.positiveTotal);
                if(entry.side==0){out.net[cell][phase].add(q);out.exactNet[cell][phase].add(q);}
                else {out.net[cell][phase].add(-q);out.exactNet[cell][phase].add(q,-1);}
            }
        }else{
            if(q==0)continue;
            const unsigned donorSide=q>0?0:1,receiverSide=q>0?1:0;
            for(unsigned phase=0;phase<2;++phase){
                const auto amount=q>0?(phase==0?flux.positiveLiquid:flux.positiveGas):
                    (phase==0?flux.negativeLiquid:flux.negativeGas);
                const auto signedPhase=phase==0?flux.liquid:flux.gas;
                if(entry.side==donorSide)out.actualOut[cell][phase].add(amount);
                if(entry.side==receiverSide&&in.counterPhase[cell]>=0)out.actualIn[cell][phase].add(amount);
                if(entry.side==0){out.net[cell][phase].add(signedPhase);out.exactNet[cell][phase].add(signedPhase);}
                else {out.net[cell][phase].add(-signedPhase);out.exactNet[cell][phase].add(signedPhase,-1);}
            }
        }
        if(entry.side==0){out.directionalTotal[cell].add(q);out.exactTotal[cell].add(q);}
        else {out.directionalTotal[cell].add(-q);out.exactTotal[cell].add(q,-1);}
    }
}
#pragma fenv_access(off)
#pragma float_control(pop)
#endif
using Gross=FaceTime::Gross;
using Box=FaceTime::Bounds;
void faceTimeDiagnostics(const Gross& measure,GeometricVofTransportReport& report)
{
    report.temporalEvents+=measure.temporalEvents;
    report.temporalIntervals+=measure.temporalIntervals;
    report.maximumTemporalRefinements=std::max(report.maximumTemporalRefinements,measure.maximumRootRefinements);
    report.maximumFaceTimeRelativeUncertainty=std::max(report.maximumFaceTimeRelativeUncertainty,measure.maximumRelativeUncertainty);
}
bool overlaps(Box a,Box b){return a.x0<b.x1&&b.x0<a.x1&&a.y0<b.y1&&b.y0<a.y1;}
struct Image{int cell;int sx,sy;long double tx,ty;};
class SourceGeometry
{
public:
    SourceGeometry(const QuadtreeMesh& mesh,const GeometricVofPublication& source):mesh_(mesh),source_(source),
        boxes_(mesh.cellLabelCapacity()),valid_(boxes_.size(),0),polygons_(source.cells().size())
    {
        const auto o=source.frameOrigin();
        lowerX_=static_cast<long double>(mesh.minCoord.x)-o.x;upperX_=static_cast<long double>(mesh.maxCoord.x)-o.x;
        lowerY_=static_cast<long double>(mesh.minCoord.y)-o.y;upperY_=static_cast<long double>(mesh.maxCoord.y)-o.y;
        width_=upperX_-lowerX_;height_=upperY_-lowerY_;
        for(const auto& c:source.cells()){
            const auto& native=mesh.cells(mesh.cellIndexToLabel[c.nativeCell]);const auto a=mesh.verts(native.cellToVert[0]).coord(),b=mesh.verts(native.cellToVert[2]).coord();
            boxes_[native.label()]={static_cast<long double>(a.x)-o.x,static_cast<long double>(b.x)-o.x,
                static_cast<long double>(a.y)-o.y,static_cast<long double>(b.y)-o.y};valid_[native.label()]=1;
            if(source.alpha()[c.nativeCell]>0&&source.alpha()[c.nativeCell]<1){
                polygons_[c.nativeCell][0]=source.anchoredPhasePolygons(c.nativeCell,true);polygons_[c.nativeCell][1]=source.anchoredPhasePolygons(c.nativeCell,false);}
        }
        for(int label=0;label<mesh.Nx*mesh.Ny;++label)nodeBox(label);
    }
    std::vector<Image> query(Box box,GeometricVofTransportReport& report) const
    {
        std::vector<Image> result;
        for(int ix=-1;ix<=1;++ix)for(int iy=-1;iy<=1;++iy){const int sx=ix==0?1:-1,sy=iy==0?1:-1;
            const long double tx=ix==0?0:(ix<0?2*lowerX_:2*upperX_),ty=iy==0?0:(iy<0?2*lowerY_:2*upperY_);
            const auto down=[](long double x){return std::nextafter(x,-std::numeric_limits<long double>::infinity());};
            const auto up=[](long double x){return std::nextafter(x,std::numeric_limits<long double>::infinity());};
            Box base{sx>0?box.x0:down(tx-box.x1),sx>0?box.x1:up(tx-box.x0),
                sy>0?box.y0:down(ty-box.y1),sy>0?box.y1:up(ty-box.y0)};
            const Box tile{down(ix<0?lowerX_-width_:(ix>0?upperX_:lowerX_)),up(ix<0?lowerX_:(ix>0?upperX_+width_:upperX_)),
                down(iy<0?lowerY_-height_:(iy>0?upperY_:lowerY_)),up(iy<0?lowerY_:(iy>0?upperY_+height_:upperY_))};
            if(!overlaps(box,tile))continue;
            const auto x0=std::max(lowerX_,base.x0),x1=std::min(upperX_,base.x1),y0=std::max(lowerY_,base.y0),y1=std::min(upperY_,base.y1);
            if(x1<=x0||y1<=y0)continue;
            const int imin=std::max(0,static_cast<int>(std::floor((x0-lowerX_)/width_*mesh_.Nx))-1);
            const int imax=std::min(mesh_.Nx-1,static_cast<int>(std::floor((x1-lowerX_)/width_*mesh_.Nx))+1);
            const int jmin=std::max(0,static_cast<int>(std::floor((y0-lowerY_)/height_*mesh_.Ny))-1);
            const int jmax=std::min(mesh_.Ny-1,static_cast<int>(std::floor((y1-lowerY_)/height_*mesh_.Ny))+1);
            for(int j=jmin;j<=jmax;++j)for(int i=imin;i<=imax;++i)visit(j*mesh_.Nx+i,base,sx,sy,tx,ty,result);
        }
        report.donorQueries+=result.size();require(!result.empty(),"actual departure support is empty");return result;
    }
    Box originalDomainBounds() const{return {lowerX_,upperX_,lowerY_,upperY_};}
    // Return the original owner points. The known reflection is applied by
    // the Core integral before coordinate differences, preserving its tails.
    std::vector<Poly> phase(const Image& image,bool liquid) const
    {
        std::vector<Poly> result;const auto a=source_.alpha()[image.cell];
        if(a==0||a==1){if((a==1)!=liquid)return result;
            const auto& c=source_.cells()[image.cell];Poly p;for(auto v:c.boundaryVertices)p.push_back(source_.anchoredVertices()[v]);result.push_back(std::move(p));}
        else result=polygons_[image.cell][liquid?0:1];
        return result;
    }
private:
    const QuadtreeMesh& mesh_;const GeometricVofPublication& source_;
    std::vector<Box> boxes_;std::vector<unsigned char> valid_;
    std::vector<std::array<std::vector<Poly>,2>> polygons_;
    long double lowerX_,upperX_,lowerY_,upperY_,width_,height_;
    Box nodeBox(int label)
    {
        require(label>=0&&static_cast<std::size_t>(label)<boxes_.size(),"invalid actual tree label");if(valid_[label])return boxes_[label];
        const auto& c=mesh_.cells(label);require(c.son[0]>=0,"unselected actual source leaf");auto box=nodeBox(c.son[0]);
        for(unsigned j=1;j<4;++j){const auto b=nodeBox(c.son[j]);box={std::min(box.x0,b.x0),std::max(box.x1,b.x1),std::min(box.y0,b.y0),std::max(box.y1,b.y1)};}
        valid_[label]=1;return boxes_[label]=box;
    }
    void visit(int label,Box box,int sx,int sy,long double tx,long double ty,std::vector<Image>& out) const
    {
        if(!overlaps(boxes_[label],box))return;const auto& c=mesh_.cells(label);
        if(c.son[0]<0){require(c.index()>=0&&static_cast<std::size_t>(c.index())<source_.alpha().size(),"source tree/alpha indexing mismatch");out.push_back({c.index(),sx,sy,tx,ty});}
        else for(int child:c.son)visit(child,box,sx,sy,tx,ty,out);
    }
};
Point samplingPoint(const GeometricVofPublication& source,const NativeLeafFace& face)
{
    Point p{face.center.x,face.center.y};if(source.domain().axisymmetric()&&face.horizontal){
        const auto& f=source.faces()[face.nativeRow];const auto origin=source.frameOrigin();
        const auto a=origin.x+nativePoint(source,f.vertices[0]).x,b=origin.x+nativePoint(source,f.vertices[1]).x;
        p.x=2*(a*a+a*b+b*b)/(3*(a+b));}
    return p;
}
}
struct GeometricVofAdvectorSlab::Boundary::Data
{
    GeometricVofIdentity source;Clock sourceTime;NativeLeafDomain domain;
    std::vector<Mode> rules;bool allClosed=true;
};
GeometricVofAdvectorSlab::Boundary::Boundary(std::shared_ptr<const Data> data):data_(std::move(data)){}
GeometricVofAdvectorSlab::Boundary GeometricVofAdvectorSlab::Boundary::prepare(
    const QuadtreeMesh& mesh,const GeometricVofPublication& source,const std::vector<Mode>& rules)
{
    source.validateFor(mesh,source.identity(),source.timeLayer());const auto clock=source.timeLayer();
    require(source.geometryKind()==GeometricVofGeometryKind::LocalPlic,"boundary requires LocalPlic source");
    require(clock.position==Clock::Position::Initial||clock.position==Clock::Position::End,
        "boundary requires accepted Initial/End source");
    const auto& domain=source.domain();
    require(domain.cellRowByIndex().size()==static_cast<std::size_t>(mesh.cellIndexNum()),
        "boundary requires actual native selected-domain layout");
    require(rules.size()==static_cast<std::size_t>(domain.nativeFaceCount()),"boundary rule dimensions differ");
    auto data=std::make_shared<Data>();data->source=source.identity();data->sourceTime=clock;data->domain=domain;data->rules=rules;
    for(int row=0;row<domain.nativeFaceCount();++row){const auto f=domain.face(mesh,row);const auto rule=rules[row];
        const bool interior=f.negativeCell>=0&&f.positiveCell>=0;
        if(interior)require(rule==Mode::Interior,"interior face must have Interior boundary rule");
        else if(f.negativeCell<0&&f.positiveCell<0)
            require(rule==Mode::Closed,"excluded face must have Closed zero-flux rule");
        else {
            require((f.negativeCell>=0)!=(f.positiveCell>=0),"boundary face must have one actual adjacent cell");
            require(rule==Mode::Closed||rule==Mode::PureLiquidReservoir||rule==Mode::PureGasReservoir,
                "exterior face requires explicit Closed or pure reservoir rule");
            if(f.area==0)require(domain.axisymmetric()&&!f.horizontal&&f.center.x==0&&rule==Mode::Closed,
                "true zero-area axis must be Closed");
            if(rule!=Mode::Closed)data->allClosed=false;
        }
    }
    return Boundary(std::move(data));
}
GeometricVofIdentity GeometricVofAdvectorSlab::Boundary::sourceIdentity() const
{require(bool(data_),"consumed boundary");return data_->source;}
Clock GeometricVofAdvectorSlab::Boundary::sourceTime() const
{require(bool(data_),"consumed boundary");return data_->sourceTime;}
const NativeLeafDomain& GeometricVofAdvectorSlab::Boundary::domain() const
{require(bool(data_),"consumed boundary");return data_->domain;}
const std::vector<GeometricVofAdvectorSlab::Boundary::Mode>& GeometricVofAdvectorSlab::Boundary::rules() const
{require(bool(data_),"consumed boundary");return data_->rules;}
bool GeometricVofAdvectorSlab::Boundary::allClosed() const
{require(bool(data_),"consumed boundary");return data_->allClosed;}
void GeometricVofAdvectorSlab::Boundary::validateFor(const QuadtreeMesh& mesh,const GeometricVofPublication& source) const
{
    require(bool(data_),"consumed boundary");source.validateFor(mesh,data_->source,data_->sourceTime);
    const GeometricVofPublication::MeshValidationBatch batch(mesh,source.domain().meshIdentity());
    source.validateDomainInBatch(mesh,data_->domain,batch);
}
void GeometricVofAdvectorSlab::Boundary::validateForInBatch(const QuadtreeMesh& mesh,
    const GeometricVofPublication& source,const GeometricVofPublication::MeshValidationBatch& batch) const
{
    require(bool(data_),"consumed boundary");
    source.validateForInBatch(mesh,data_->source,data_->sourceTime,batch);
    source.validateDomainInBatch(mesh,data_->domain,batch);
}
bool GeometricVofAdvectorSlab::Boundary::sameInstance(const Boundary& other) const noexcept
{return data_&&data_==other.data_;}
struct GeometricVofAdvectorSlab::Data
{
    GeometricVofIdentity source;Clock sourceTime,endTime;NativeLeafDomain domain;
    std::vector<double> beginNormal,endNormal;std::vector<Point> beginVelocity,endVelocity;
    std::optional<Boundary> boundary;
};
GeometricVofAdvectorSlab::GeometricVofAdvectorSlab(std::shared_ptr<const Data> data):data_(std::move(data)){}
GeometricVofAdvectorSlab GeometricVofAdvectorSlab::prepare(const QuadtreeMesh& mesh,
    const GeometricVofPublication& source,Clock nextEnd,const std::vector<double>& begin,const std::vector<double>& end)
{return prepareImpl(mesh,source,nextEnd,begin,end,nullptr);}
GeometricVofAdvectorSlab GeometricVofAdvectorSlab::prepareBoundary(const QuadtreeMesh& mesh,
    const GeometricVofPublication& source,const Boundary& boundary,Clock nextEnd,
    const std::vector<double>& begin,const std::vector<double>& end)
{
    boundary.validateFor(mesh,source);
    if(boundary.allClosed())return prepare(mesh,source,nextEnd,begin,end);
    return prepareImpl(mesh,source,nextEnd,begin,end,&boundary);
}
const GeometricVofAdvectorSlab::Boundary* GeometricVofAdvectorSlab::boundary() const
{require(bool(data_),"consumed slab");return data_->boundary?&*data_->boundary:nullptr;}
GeometricVofAdvectorSlab GeometricVofAdvectorSlab::prepareImpl(const QuadtreeMesh& mesh,const GeometricVofPublication& source,
    Clock nextEnd,const std::vector<double>& begin,const std::vector<double>& end,const Boundary* boundary)
{
    source.validateMesh(mesh);nextEnd.time();const auto old=source.timeLayer();old.time();
    require(old.position!=Clock::Position::Midpoint&&nextEnd.position==Clock::Position::End&&
        old.step<std::numeric_limits<std::int64_t>::max()&&nextEnd.step==old.step+1&&nextEnd.begin==old.time()&&nextEnd.end>nextEnd.begin,
        "slab requires accepted source and complete next End clock");
    auto midpoint=nextEnd;midpoint.position=Clock::Position::Midpoint;
    require(midpoint.time()>nextEnd.begin&&midpoint.time()<nextEnd.end,"unrepresentable complete interval midpoint");
    const auto& domain=source.domain();const auto count=static_cast<std::size_t>(domain.nativeFaceCount());
    require(begin.size()==count&&end.size()==count,"native normal input does not cover actual shared subfaces");
    auto data=std::make_shared<Data>();data->source=source.identity();data->sourceTime=old;data->endTime=nextEnd;data->domain=domain;
    data->beginNormal=begin;data->endNormal=end;data->beginVelocity.resize(count);data->endVelocity.resize(count);
    if(boundary){boundary->validateFor(mesh,source);data->boundary=*boundary;}
    const bool local=source.geometryKind()==GeometricVofGeometryKind::LocalPlic;
    std::vector<std::vector<int>> incident(local?0:source.cells().size());
    for(int row=0;row<domain.nativeFaceCount();++row){const auto f=domain.face(mesh,row);
        require(std::isfinite(begin[row])&&std::isfinite(end[row]),"nonfinite native endpoint mean");
        if(!boundary){
        if(f.negativeCell<0||f.positiveCell<0||f.area==0)require(begin[row]==0&&end[row]==0,"VOF stage currently requires impermeable exterior/axis normals");
        }else if(f.negativeCell<0||f.positiveCell<0||f.area==0){
            const auto mode=boundary->rules()[row];
            if(mode==Boundary::Mode::Closed)require(begin[row]==0&&end[row]==0,"Closed boundary/axis normals must be exactly zero");
            else if(!local&&((begin[row]<0&&end[row]>0)||(begin[row]>0&&end[row]<0))){
                // Failure-only: actual slab snapshots, not an accepted future
                // state or a cell-centred velocity estimate.
                const char* side=f.horizontal?(f.outwardSign<0?"S":"N"):(f.outwardSign<0?"W":"E");
                std::ostringstream detail;detail.precision(std::numeric_limits<long double>::max_digits10);
                detail<<"same-slab boundary velocity sign reversal is unsupported"
                    <<" native_row="<<row<<" side="<<side<<" normal_axis="<<(f.horizontal?1:0)
                    <<" outward_sign="<<f.outwardSign<<" negative_cell="<<f.negativeCell<<" positive_cell="<<f.positiveCell
                    <<" reservoir_mode="<<static_cast<int>(mode)<<" source_identity="<<source.identity().value()
                    <<" source_clock=("<<old.step<<','<<old.begin<<','<<old.end<<','<<static_cast<int>(old.position)<<')'
                    <<" slab_end_clock=("<<nextEnd.step<<','<<nextEnd.begin<<','<<nextEnd.end<<','<<static_cast<int>(nextEnd.position)<<')'
                    <<" face_center=("<<f.center.x<<','<<f.center.y<<") face_area="<<f.area
                    <<" begin_normal="<<begin[row]<<" end_normal="<<end[row]
                    <<" begin_outward="<<f.outwardSign*begin[row]<<" end_outward="<<f.outwardSign*end[row]
                    <<" begin_normal_hex="<<std::hexfloat<<begin[row]<<" end_normal_hex="<<end[row];
                require(false,detail.str());
            }
        }
        if(!local){if(f.negativeCell>=0)incident[f.negativeCell].push_back(row);if(f.positiveCell>=0)incident[f.positiveCell].push_back(row);}
    }
    // Local swept strips need only the original canonical normal snapshot.
    // Do not build a tangential P1 fit or expose normal-only data as vectors.
    if(source.geometryKind()==GeometricVofGeometryKind::LocalPlic){
        data->beginVelocity.clear();data->endVelocity.clear();
        return GeometricVofAdvectorSlab(std::move(data));
    }
    for(int row=0;row<domain.nativeFaceCount();++row){const auto f=domain.face(mesh,row);
        Point a{f.horizontal?0:begin[row],f.horizontal?begin[row]:0},b{f.horizontal?0:end[row],f.horizontal?end[row]:0};
        if(f.negativeCell>=0&&f.positiveCell>=0){std::set<int> samples;
            for(int cell:{f.negativeCell,f.positiveCell})for(int other:incident[cell])if(domain.face(mesh,other).horizontal!=f.horizontal)samples.insert(other);
            require(samples.size()>=3,"insufficient native tangential support row="+std::to_string(row));
            const auto point=samplingPoint(source,f);const auto scale=std::min(f.normalSpacing,f.tangentSpacing);require(scale>0,"invalid native reconstruction scale");
            Eigen::MatrixXd design(samples.size(),3),values(samples.size(),2);int i=0;
            for(int other:samples){const auto p=samplingPoint(source,domain.face(mesh,other));const auto x=(p.x-point.x)/scale,y=(p.y-point.y)/scale;
                const auto weight=1/std::sqrt(1+x*x+y*y);design(i,0)=static_cast<double>(weight);design(i,1)=static_cast<double>(weight*x);design(i,2)=static_cast<double>(weight*y);
                values(i,0)=static_cast<double>(weight)*begin[other];values(i,1)=static_cast<double>(weight)*end[other];++i;}
            Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(design);qr.setThreshold(64*std::numeric_limits<double>::epsilon());
            require(qr.rank()==3,"rank-deficient native tangential support row="+std::to_string(row));const Eigen::MatrixXd coefficients=qr.solve(values);
            require(coefficients.allFinite(),"nonfinite local tangential reconstruction");
            if(f.horizontal){a.x=coefficients(0,0);b.x=coefficients(0,1);}else {a.y=coefficients(0,0);b.y=coefficients(0,1);}}
        data->beginVelocity[row]=a;data->endVelocity[row]=b;
    }
    return GeometricVofAdvectorSlab(std::move(data));
}
GeometricVofIdentity GeometricVofAdvectorSlab::sourceIdentity() const{require(bool(data_),"consumed slab");return data_->source;}
Clock GeometricVofAdvectorSlab::sourceTime() const{require(bool(data_),"consumed slab");return data_->sourceTime;}
Clock GeometricVofAdvectorSlab::endTime() const{require(bool(data_),"consumed slab");return data_->endTime;}
const NativeLeafDomain& GeometricVofAdvectorSlab::domain() const{require(bool(data_),"consumed slab");return data_->domain;}
const std::vector<double>& GeometricVofAdvectorSlab::beginNormal() const{require(bool(data_),"consumed slab");return data_->beginNormal;}
const std::vector<double>& GeometricVofAdvectorSlab::endNormal() const{require(bool(data_),"consumed slab");return data_->endNormal;}
const std::vector<Point>& GeometricVofAdvectorSlab::beginVelocity() const{require(bool(data_),"consumed slab");require(!data_->beginVelocity.empty(),"LocalPlic slab retains normal means only");return data_->beginVelocity;}
const std::vector<Point>& GeometricVofAdvectorSlab::endVelocity() const{require(bool(data_),"consumed slab");require(!data_->endVelocity.empty(),"LocalPlic slab retains normal means only");return data_->endVelocity;}
void GeometricVofAdvectorSlab::validateFor(const QuadtreeMesh& mesh,const GeometricVofPublication& source) const
{
    require(bool(data_),"consumed slab");source.validateFor(mesh,data_->source,data_->sourceTime);
    const GeometricVofPublication::MeshValidationBatch batch(mesh,source.domain().meshIdentity());
    source.validateDomainInBatch(mesh,data_->domain,batch);
    if(data_->boundary)data_->boundary->validateForInBatch(mesh,source,batch);
}
void GeometricVofAdvectorSlab::validateForInBatch(const QuadtreeMesh& mesh,
    const GeometricVofPublication& source,const GeometricVofPublication::MeshValidationBatch& batch) const
{
    require(bool(data_),"consumed slab");
    source.validateForInBatch(mesh,data_->source,data_->sourceTime,batch);
    source.validateDomainInBatch(mesh,data_->domain,batch);
    if(data_->boundary)data_->boundary->validateForInBatch(mesh,source,batch);
}
bool GeometricVofAdvectorSlab::sameInstance(const GeometricVofAdvectorSlab& other) const noexcept
{return data_&&data_==other.data_;}
namespace
{
struct StageResult
{
    std::vector<double> alpha;
    std::shared_ptr<const Detail::GeometricVofPhaseBudgetRestriction> restriction;
};
StageResult boundedStage(const GeometricVofPublication& source,std::vector<GeometricVofTransportFace>& faces,bool half,
    long double rhoL,long double rhoG,const GeometricVofTransportConfig& config,GeometricVofTransportReport& report,
    const std::vector<GeometricVofTransportFace>& geometry)
{
    const auto size=source.cells().size();std::vector<Detail::GeometricVofPhaseBudgetCell> cells(size);
    for(std::size_t c=0;c<size;++c)cells[c]={source.cells()[c].capacity,source.alpha()[c]};
    const auto budget=Detail::prepareGeometricVofPhaseBudget(cells,faces,half,config.maximumOutflowCfl,&geometry,config.relativeConservationTolerance);
    report.maximumOutflowCfl=std::max(report.maximumOutflowCfl,budget.maximumOutflowCfl);
    report.maximumTotalDivergence=std::max(report.maximumTotalDivergence,budget.maximumTotalDivergence);
    report.maximumCapacityReconciliation=std::max(report.maximumCapacityReconciliation,budget.maximumCapacityReconciliation);
    report.limitedFaces+=budget.limitedFaceExtents;std::vector<Sum> phaseDivergence(size),gasDivergence(size);
    for(std::size_t row=0;row<faces.size();++row){auto& face=faces[row];auto& f=half?face.half:face.full;f=budget.flux[row];
        if(!(f.positiveLiquid>=0&&f.negativeLiquid>=0&&f.positiveGas>=0&&f.negativeGas>=0&&
            f.positiveLiquid<=f.positiveTotal&&f.negativeLiquid<=f.negativeTotal&&
            f.positiveGas<=f.positiveTotal&&f.negativeGas<=f.negativeTotal)){
            std::ostringstream detail;detail.precision(17);detail<<"accepted phase gross exceeds native total row="<<row
                <<" extent="<<(half?"half":"full")<<" positiveTotal="<<f.positiveTotal<<" positiveLiquid="<<f.positiveLiquid
                <<" positiveGas="<<f.positiveGas<<" negativeTotal="<<f.negativeTotal<<" negativeLiquid="<<f.negativeLiquid
                <<" negativeGas="<<f.negativeGas<<" limiter="<<f.limiter;
            require(false,detail.str());}
        f.mass=rhoL*f.liquid+rhoG*f.gas;
        require(std::isfinite(f.mass),"unrepresentable shared material integral row="+std::to_string(row));
        const auto gross=f.positiveTotal+f.negativeTotal;
        require((f.positiveTotal==0||(std::abs(f.positiveLiquid+f.positiveGas-f.positiveTotal)/f.positiveTotal<=arithmetic))&&
            (f.negativeTotal==0||(std::abs(f.negativeLiquid+f.negativeGas-f.negativeTotal)/f.negativeTotal<=arithmetic)),"accepted gross two-phase flux lost total closure");
        require(gross==0||(std::abs(f.total-f.liquid-f.gas)/gross<=arithmetic),"accepted shared two-phase flux lost total closure");
        report.maximumFluxCorrection=std::max(report.maximumFluxCorrection,std::abs(f.liquid-f.rawLiquid));
        if(face.negativeCell>=0){phaseDivergence[face.negativeCell].add(f.liquid);gasDivergence[face.negativeCell].add(f.gas);}
        if(face.positiveCell>=0){phaseDivergence[face.positiveCell].add(-f.liquid);gasDivergence[face.positiveCell].add(-f.gas);}}
    StageResult result;result.restriction=budget.restriction;result.alpha.resize(size);Sum oldL,oldG,newL,newG,totalCapacity;
    for(std::size_t c=0;c<size;++c){const auto capacity=source.cells()[c].capacity;
        const auto l=budget.phaseInventory[c][0],g=budget.phaseInventory[c][1];
        require(std::isfinite(l)&&std::isfinite(g)&&l>=0&&g>=0&&l+g>0,"synchronized phase budgets lost bounds cell="+std::to_string(c));
        const auto defect=std::abs(l+g-capacity)/capacity;report.maximumCapacityReconciliation=std::max(report.maximumCapacityReconciliation,defect);
        require(defect<=arithmetic,"phase capacity reconciliation exceeds arithmetic budget");
        if(l==0)result.alpha[c]=0;else if(g==0)result.alpha[c]=1;
        else {result.alpha[c]=static_cast<double>(l/(l+g));if(!(result.alpha[c]>0&&result.alpha[c]<1)){std::ostringstream detail;detail.precision(17);detail<<"double alpha cannot retain two positive transported phases cell="<<c<<" extent="<<(half?"half":"full")<<" oldAlpha="<<source.alpha()[c]<<" capacity="<<capacity<<" liquid="<<l<<" gas="<<g;require(false,detail.str());}}
        const auto actual=static_cast<long double>(result.alpha[c])*capacity;
        const auto residual=std::max(std::abs(actual-source.cells()[c].inventory+phaseDivergence[c].value),
            std::abs(source.cells()[c].inventory-actual+gasDivergence[c].value))/capacity;
        report.maximumInventoryResidual=std::max(report.maximumInventoryResidual,residual);
        require(residual<=config.relativeConservationTolerance,"inventory/accepted shared flux ledger mismatch");
        oldL.add(source.cells()[c].inventory);oldG.add(capacity-source.cells()[c].inventory);newL.add(actual);newG.add(capacity-actual);totalCapacity.add(capacity);}
    const auto relative=[&](long double a,long double b){return std::abs(a-b)/(b>0?b:totalCapacity.value);};
    const auto error=std::max(relative(newL.value,oldL.value),relative(newG.value,oldG.value));
    require(error<=config.relativeConservationTolerance,"closed two-phase global inventory failed");
    (half?report.halfRelativeConservation:report.fullRelativeConservation)=error;return result;
}
}
struct GeometricVofTransport::Data
{
    GeometricVofIdentity source;Clock sourceTime,endTime;
    GeometricVofPublication midpoint,endpoint;
    std::vector<GeometricVofTransportFace> faces;GeometricVofTransportReport report;
    GeometricVofPublication original;
    GeometricVofAdvectorSlab slab;
    long double rhoL,rhoG;
    GeometricVofTransportConfig config;
    std::optional<VofSmoothGeometry> smoothSource,smoothMidpoint,smoothEndpoint;
    std::shared_ptr<const Detail::GeometricVofPhaseBudgetRestriction> halfRestriction,fullRestriction;
    Data(const GeometricVofPublication& old,const GeometricVofAdvectorSlab& velocity,long double l,long double g,
        GeometricVofTransportConfig cfg,GeometricVofPublication mid,GeometricVofPublication last,
        std::vector<GeometricVofTransportFace> f,GeometricVofTransportReport r,
        std::shared_ptr<const Detail::GeometricVofPhaseBudgetRestriction> h,
        std::shared_ptr<const Detail::GeometricVofPhaseBudgetRestriction> e):
        source(old.identity()),sourceTime(old.timeLayer()),endTime(velocity.endTime()),
        midpoint(std::move(mid)),endpoint(std::move(last)),faces(std::move(f)),report(r),
        original(old),slab(velocity),rhoL(l),rhoG(g),config(cfg),halfRestriction(std::move(h)),fullRestriction(std::move(e)){}
};
GeometricVofTransport::GeometricVofTransport(std::shared_ptr<const Data> data):data_(std::move(data)){}
GeometricVofTransport GeometricVofTransport::prepare(const QuadtreeMesh& mesh,const GeometricVofPublication& source,
    const GeometricVofAdvectorSlab& slab,long double rhoL,long double rhoG,const GeometricVofTransportConfig& config)
{
    require(source.geometryKind()==GeometricVofGeometryKind::CompatibleFan,"fan transport requires CompatibleFan geometry; use explicit local sweep");
    require(config.geometry.localResolution==GeometricVofLocalResolution::VolumeOnly,
        "native-interface resolution requires the explicit local sweep");
    slab.validateFor(mesh,source);require(std::isfinite(rhoL)&&rhoL>0&&std::isfinite(rhoG)&&rhoG>0,"densities must be finite and positive");
    require(std::isfinite(config.maximumOutflowCfl)&&config.maximumOutflowCfl>0&&config.maximumOutflowCfl<=1,"invalid gross outflow CFL");
    require(std::isfinite(config.relativeConservationTolerance)&&config.relativeConservationTolerance>0&&config.relativeConservationTolerance<=1e-12L,"cannot loosen conservation gate");
    require(config.maximumTemporalEvents>=3&&config.maximumTemporalEvents<=4096,"invalid finite temporal budget");
    const auto end=slab.endTime();GeometricVofTransportReport report;
    SourceGeometry geometry(mesh,source);std::vector<GeometricVofTransportFace> faces,rawGeometry;faces.reserve(source.faces().size());rawGeometry.reserve(source.faces().size());const bool rz=source.domain().axisymmetric();
    for(const auto& actual:source.faces()){
        const auto native=source.domain().face(mesh,actual.nativeRow);GeometricVofTransportFace face;face.nativeRow=actual.nativeRow;face.negativeCell=native.negativeCell;face.positiveCell=native.positiveCell;
        const auto v0=slab.beginVelocity()[actual.nativeRow],v1=slab.endVelocity()[actual.nativeRow];
        const auto total=FaceTime::integrateNormal(slab.beginNormal()[actual.nativeRow],slab.endNormal()[actual.nativeRow],end.begin,end.end,native.area);
        report.temporalEvents+=total.temporalEvents;
        report.maximumTemporalRefinements=std::max(report.maximumTemporalRefinements,total.maximumRootRefinements);
        report.maximumFaceTimeRelativeUncertainty=std::max(report.maximumFaceTimeRelativeUncertainty,total.maximumRelativeUncertainty);
        face.half.positiveTotal=total.halfPositive;face.half.negativeTotal=total.halfNegative;face.full.positiveTotal=total.fullPositive;face.full.negativeTotal=total.fullNegative;
        face.half.total=total.halfPositive-total.halfNegative;face.full.total=total.fullPositive-total.fullNegative;
        auto raw=face;
        if(total.fullPositive+total.fullNegative>0){
            const auto a=nativePoint(source,actual.vertices[0]),b=nativePoint(source,actual.vertices[1]);
            FaceTime::validateDepartureSupport(a,b,v0,v1,end.begin,end.end,geometry.originalDomainBounds());
            const auto images=geometry.query(FaceTime::departureBounds(a,b,v0,v1,end.begin,end.end),report);
            bool pureGas=true,pureLiquid=true;for(auto image:images){pureGas&=source.alpha()[image.cell]==0;pureLiquid&=source.alpha()[image.cell]==1;}
            Gross liquid,gas;
            if(pureGas)gas=total;else if(pureLiquid)liquid=total;
            else for(auto image:images)for(bool phase:{true,false})for(const auto& p:geometry.phase(image,phase)){
                Gross measure;
                try {measure=FaceTime::integrateAnchoredPolygon(p,a,b,v0,v1,end.begin,end.end,source.frameOrigin(),rz,
                    config.maximumTemporalEvents,{image.sx,image.sy,{image.tx,image.ty}});}
                catch(const std::exception& error){std::ostringstream detail;detail.precision(17);
                    detail<<"[Interface/GeometricVofTransport] actual face-time refused row="<<actual.nativeRow<<" sourceCell="<<image.cell
                        <<" liquid="<<phase<<" clock="<<end.begin<<","<<end.end<<" reflection="<<image.sx<<","<<image.sy<<","<<image.tx<<","<<image.ty
                        <<": "<<error.what();throw std::invalid_argument(detail.str());}
                ++report.polygonIntegrals;faceTimeDiagnostics(measure,report);auto& aggregate=phase?liquid:gas;
                aggregate.halfPositive+=measure.halfPositive;aggregate.halfNegative+=measure.halfNegative;aggregate.fullPositive+=measure.fullPositive;aggregate.fullNegative+=measure.fullNegative;}
            raw.half.positiveLiquid=liquid.halfPositive;raw.half.negativeLiquid=liquid.halfNegative;
            raw.half.positiveGas=gas.halfPositive;raw.half.negativeGas=gas.halfNegative;
            raw.full.positiveLiquid=liquid.fullPositive;raw.full.negativeLiquid=liquid.fullNegative;
            raw.full.positiveGas=gas.fullPositive;raw.full.negativeGas=gas.fullNegative;
            const auto normalize=[&](long double l,long double g,long double t){
                using Pair=std::pair<long double,long double>;
                if(t==0){require(l==0&&g==0,"zero total flux has geometric phase transport");return Pair{0,0};}
                require(l>=0&&g>=0&&std::isfinite(l+g)&&l+g>0,"invalid geometric two-phase face integral");
                const auto residual=std::abs(l+g-t)/t;report.maximumCoverageResidual=std::max(report.maximumCoverageResidual,residual);
                require(residual<=config.relativeConservationTolerance,"incomplete/overlapping actual fan face-time support row="+std::to_string(actual.nativeRow));
                if(l==0)return Pair{0,t};if(g==0)return Pair{t,0};
                const auto minor=t*((l<=g?l:g)/(l+g));
                require(std::isfinite(minor)&&minor>0&&t-minor>0,"unrepresentable positive geometric phase integral");
                return l<=g?Pair{minor,t-minor}:Pair{t-minor,minor};};
            const auto hp=normalize(liquid.halfPositive,gas.halfPositive,total.halfPositive),hn=normalize(liquid.halfNegative,gas.halfNegative,total.halfNegative);
            const auto fp=normalize(liquid.fullPositive,gas.fullPositive,total.fullPositive),fn=normalize(liquid.fullNegative,gas.fullNegative,total.fullNegative);
            face.half.positiveLiquid=hp.first;face.half.negativeLiquid=hn.first;face.half.positiveGas=hp.second;face.half.negativeGas=hn.second;
            face.full.positiveLiquid=fp.first;face.full.negativeLiquid=fn.first;face.full.positiveGas=fp.second;face.full.negativeGas=fn.second;
            face.half.rawLiquid=face.half.positiveLiquid-face.half.negativeLiquid;face.full.rawLiquid=face.full.positiveLiquid-face.full.negativeLiquid;
        }
        faces.push_back(face);rawGeometry.push_back(raw);
    }
    const auto half=boundedStage(source,faces,true,rhoL,rhoG,config,report,rawGeometry),full=boundedStage(source,faces,false,rhoL,rhoG,config,report,rawGeometry);
    auto mid=end;mid.position=Clock::Position::Midpoint;
    const auto publication=[&](Clock clock,const std::vector<double>& alpha,const char* extent){
        try{return GeometricVofPublication::prepare(mesh,source.domain(),clock,alpha,config.geometry);}
        catch(const std::exception& error){std::ostringstream detail;detail.precision(17);
            detail<<"[Interface/GeometricVofTransport] "<<extent<<" reconstruction refused: "<<error.what();
            const std::string reason=error.what();const auto pos=reason.find("cell=");
            if(pos!=std::string::npos){std::istringstream index(reason.substr(pos+5));std::size_t cell;
                if(index>>cell&&cell<alpha.size())detail<<" transportedAlpha="<<alpha[cell]<<" sourceAlpha="<<source.alpha()[cell]<<" capacity="<<source.cells()[cell].capacity;}
            throw std::invalid_argument(detail.str());}};
    auto midpoint=publication(mid,half.alpha,"half"),endpoint=publication(end,full.alpha,"full");
    return GeometricVofTransport(std::make_shared<Data>(source,slab,rhoL,rhoG,config,std::move(midpoint),std::move(endpoint),std::move(faces),report,half.restriction,full.restriction));
}
namespace
{
std::array<Point,2> localFaceEndpoints(const QuadtreeMesh& mesh,const NativeLeafFace& f)
{
    const auto& e=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
    const auto a=mesh.verts(e.edgeToVert[0]).coord(),b=mesh.verts(e.edgeToVert[1]).coord();
    return {{{a.x,a.y},{b.x,b.y}}};
}
Cut::Rectangle localColumn(const QuadtreeMesh& mesh,const NativeLeafFace& face,
    const Cut::Rectangle& donor)
{
    auto column=donor;const auto e=localFaceEndpoints(mesh,face);
    if(face.horizontal){column.lower.x=static_cast<double>(std::min(e[0].x,e[1].x));column.upper.x=static_cast<double>(std::max(e[0].x,e[1].x));}
    else {column.lower.y=static_cast<double>(std::min(e[0].y,e[1].y));column.upper.y=static_cast<double>(std::max(e[0].y,e[1].y));}
    require(column.lower.x>=donor.lower.x&&column.upper.x<=donor.upper.x&&
        column.lower.y>=donor.lower.y&&column.upper.y<=donor.upper.y,"actual subface column outside donor");
    return column;
}
long double localColumnRatio(const QuadtreeMesh& mesh,const NativeLeafFace& face,
    const Cut::Rectangle& donor,long double absoluteVolume,bool rz)
{
    const auto column=localColumn(mesh,face,donor);
    return absoluteVolume/Cut::capacity(column,rz);
}
using LocalPin=GeometricVofPublication::PinCertificate;
std::optional<Point> localPinFarEndpoint(const QuadtreeMesh& mesh,const NativeLeafFace& face,const LocalPin* pin)
{
    if(!pin||face.negativeCell<0||face.positiveCell<0||face.area<=0)return std::nullopt;
    const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
    const int a=mesh.verts(edge.edgeToVert[0]).index(),b=mesh.verts(edge.edgeToVert[1]).index();
    if(a!=pin->nativeVertex()&&b!=pin->nativeVertex())return std::nullopt;
    const auto e=localFaceEndpoints(mesh,face);const auto point=pin->physicalPoint();
    const auto at=a==pin->nativeVertex()?e[0]:e[1];
    require(at.x==point.x&&at.y==point.y,"sealed pin native endpoint and physical point differ");
    return a==pin->nativeVertex()?e[1]:e[0];
}
std::vector<Point> localPinPolygon(const NativeLeafFace& face,Point origin,Point far,long double depth,bool positive)
{
    require(std::isfinite(depth)&&depth>0,"unrepresentable pinned sweep depth");
    const Point tangent{far.x-origin.x,far.y-origin.y};auto inner=tangent;
    if(face.horizontal)inner.y+=positive?-depth:depth;
    else inner.x+=positive?-depth:depth;
    return {{0,0},tangent,inner};
}
long double localPinWidth(const QuadtreeMesh& mesh,const NativeLeafFace& face,const Cut::Rectangle& donor,
    const LocalPin& pin,bool positive)
{
    const auto column=localColumn(mesh,face,donor);const auto point=pin.physicalPoint();
    require(face.horizontal?point.y==(positive?column.upper.y:column.lower.y):
        point.x==(positive?column.upper.x:column.lower.x),"pin sweep is not on the actual donor side");
    return face.horizontal?static_cast<long double>(column.upper.y)-column.lower.y:
        static_cast<long double>(column.upper.x)-column.lower.x;
}
// For every permitted c in (0,.5], V(c*w)>=c*capacityBound. Radial
// inward flow has a positive quadratic: its infinitesimal bound is needed,
// since its peak/mean factor can exceed three near the actual axis.
long double localPinRateCapacity(const QuadtreeMesh& mesh,const NativeLeafFace& face,
    const Cut::Rectangle& donor,const LocalPin& pin,Point far,bool positive,bool rz)
{
    const auto width=localPinWidth(mesh,face,donor,pin,positive);
    const auto site=pin.physicalPoint();const Point origin{site.x,site.y};
    const auto bound=rz&&!face.horizontal&&!positive?static_cast<long double>(face.area)*width/2:
        2*Cut::polygonMeasureInFrame(localPinPolygon(face,origin,far,width/2,positive),origin,rz);
    require(std::isfinite(bound)&&bound>0,"unrepresentable pinned peak-depth rate capacity");return bound;
}
struct LocalPinnedSweep {std::vector<Point> polygon;Point origin;long double actual=0,peakCourant=0;};
LocalPinnedSweep localPinnedSweep(const QuadtreeMesh& mesh,const NativeLeafFace& face,const Cut::Rectangle& donor,
    const LocalPin& pin,Point far,long double total,bool positive,bool rz,long double maximumCourant,long double tolerance)
{
    const auto width=localPinWidth(mesh,face,donor,pin,positive),limit=maximumCourant*width;
    const auto site=pin.physicalPoint();const Point origin{site.x,site.y};
    const auto measure=[&](long double depth){return Cut::polygonMeasureInFrame(localPinPolygon(face,origin,far,depth,positive),origin,rz);};
    const auto maximum=measure(limit);
    require(std::isfinite(total)&&total>0&&std::isfinite(maximum)&&maximum>0,"invalid pinned swept volume");
    if(total>maximum)throw GeometricVofCourantExceeded("pinned actual peak donor-depth Courant exceeded");
    long double low=0,high=limit,depth=(total/maximum)*limit;
    require(std::isfinite(depth)&&depth>0,"positive pinned sweep depth is unrepresentable");
    // Tighten in the depth's own scale before bisection; tiny nonzero FT is
    // never rounded to an empty strip by a capacity-scaled stopping rule.
    auto value=measure(depth);
    if(value<total){low=depth;
        while(value<total&&depth<limit){const auto next=std::min(limit,2*depth);
            require(next>depth,"pinned sweep cannot bracket positive volume");depth=next;value=measure(depth);if(value<total)low=depth;}
        high=depth;
    }else high=depth;
    long double best=depth,bestValue=value,bestError=std::abs(value/total-1);
    for(unsigned iteration=0;iteration<128&&low<high&&bestError!=0;++iteration){
        depth=low+(high-low)/2;if(depth==low||depth==high)break;value=measure(depth);
        const auto error=std::abs(value/total-1);if(error<bestError){best=depth;bestValue=value;bestError=error;}
        if(value==total)break;if(value<total)low=depth;else high=depth;
    }
    auto polygon=localPinPolygon(face,origin,far,best,positive);
    const auto actual=Cut::polygonMeasureInFrame(polygon,origin,rz);
    const auto column=localColumn(mesh,face,donor);const auto capacity=Cut::capacity(donor,rz);
    require(actual==bestValue&&std::isfinite(actual)&&actual>0&&best/width<=maximumCourant,
        "invalid final pinned swept polygon");
    require(std::abs(actual-total)/capacity<=tolerance&&bestError<=128*std::numeric_limits<long double>::epsilon(),
        "pinned polygon cannot represent canonical requested volume");
    const auto clipped=Cut::clipPolygonInFrame(column,polygon,origin);
    require(std::abs(Cut::polygonMeasureInFrame(clipped,origin,rz)-actual)<=tolerance*capacity,
        "pinned polygon leaves actual native donor column");
    return {std::move(polygon),origin,actual,best/width};
}
std::array<long double,2> localPinnedPhases(const QuadtreeMesh& mesh,const NativeLeafFace& face,
    const GeometricVofPlicCell& cell,const LocalPin& pin,Point far,long double total,bool positive,bool rz,
    long double maximumCourant,long double tolerance,GeometricVofTransportReport& report)
{
    const auto sweep=localPinnedSweep(mesh,face,cell.rectangle,pin,far,total,positive,rz,maximumCourant,tolerance);
    ++report.localPinnedSweeps;report.maximumPinnedPeakCfl=std::max(report.maximumPinnedPeakCfl,sweep.peakCourant);
    const auto capacity=Cut::capacity(cell.rectangle,rz);
    report.maximumLocalStripResidual=std::max(report.maximumLocalStripResidual,std::abs(sweep.actual-total)/capacity);
    if(cell.inventory[0]==0)return {0,total};if(cell.inventory[1]==0)return {total,0};
    if(!cell.resolved){
        ++report.localUnresolvedDonors;const bool liquidMinor=cell.inventory[0]<=cell.inventory[1];
        const auto minor=total*((liquidMinor?cell.inventory[0]:cell.inventory[1])/capacity);
        return liquidMinor?std::array<long double,2>{minor,total-minor}:std::array<long double,2>{total-minor,minor};
    }
    long double l,g;
    const auto donor=positive?face.negativeCell:face.positiveCell;
    if(donor==pin.nativeCell()){
        // The sealed corner plane's same normalization and same corner
        // expression; this is not a small-value snap in a translated frame.
        const auto nx=static_cast<long double>(cell.normal.x),ny=static_cast<long double>(cell.normal.y),norm=std::hypot(nx,ny);
        const auto hx=(static_cast<long double>(cell.rectangle.upper.x)-cell.rectangle.lower.x)/2;
        const auto hy=(static_cast<long double>(cell.rectangle.upper.y)-cell.rectangle.lower.y)/2;
        const auto c=pin.corner();const auto sx=c==Cut::RectangleCorner::LowerLeft||c==Cut::RectangleCorner::UpperLeft?-1.L:1.L;
        const auto sy=c==Cut::RectangleCorner::LowerLeft||c==Cut::RectangleCorner::LowerRight?-1.L:1.L;
        require(sweep.origin.x==(sx<0?cell.rectangle.lower.x:cell.rectangle.upper.x)&&
            sweep.origin.y==(sy<0?cell.rectangle.lower.y:cell.rectangle.upper.y),"sealed pin rectangle corner differs");
        const Point normal{nx/norm,ny/norm};const auto threshold=cell.offset-(normal.x*(sx*hx)+normal.y*(sy*hy));
        require(threshold==0,"sealed pin final normal/offset lost exact corner expression");
        const auto polygon=Cut::clipPolygonInFrame(cell.rectangle,sweep.polygon,sweep.origin);
        l=Cut::polygonMeasureInFrame(Cut::clipConvexPolygon(polygon,normal,threshold),sweep.origin,rz);
        g=Cut::polygonMeasureInFrame(Cut::clipConvexPolygon(polygon,{-normal.x,-normal.y},-threshold),sweep.origin,rz);
    }else{
        const auto liquid=Cut::overlapInFrame(cell.rectangle,sweep.polygon,sweep.origin,cell.normal,cell.offset,rz);
        const auto gas=Cut::overlapInFrame(cell.rectangle,sweep.polygon,sweep.origin,{-cell.normal.x,-cell.normal.y},-cell.offset,rz);
        require(std::abs(liquid.total-sweep.actual)<=tolerance*capacity&&std::abs(gas.total-sweep.actual)<=tolerance*capacity,
            "pinned independent phase polygons lost donor coverage");l=liquid.liquid;g=gas.liquid;
    }
    ++report.localStripClips;
    require(l>=0&&g>=0&&std::isfinite(l+g)&&l+g>0&&std::abs(l+g-sweep.actual)<=tolerance*capacity,
        "invalid pinned independent two-phase polygon coverage");
    const bool liquidMinor=l<=g;const auto minor=total*((liquidMinor?l:g)/(l+g));
    return liquidMinor?std::array<long double,2>{minor,total-minor}:std::array<long double,2>{total-minor,minor};
}
Cut::Rectangle localStripRectangle(const NativeLeafFace& face,const Cut::Rectangle& column,
    long double total,bool positive,bool rz)
{
    auto strip=column;
    if(face.horizontal){
        const auto section=(static_cast<long double>(column.upper.x)-column.lower.x)*
            (rz?CartesianMeasure::twoPi*(static_cast<long double>(column.lower.x)+column.upper.x)/2:1.L);
        const auto depth=total/section;
        if(positive)strip.lower.y=static_cast<double>(static_cast<long double>(column.upper.y)-depth);
        else strip.upper.y=static_cast<double>(static_cast<long double>(column.lower.y)+depth);
    }else{
        const auto length=static_cast<long double>(column.upper.y)-column.lower.y;
        long double depth;
        if(!rz)depth=total/length;
        else {
            const auto radius=positive?static_cast<long double>(column.upper.x):static_cast<long double>(column.lower.x);
            const auto delta=total/((CartesianMeasure::twoPi/2)*length);
            const auto square=radius*radius+(positive?-delta:delta);
            require(square>=0&&std::isfinite(square),"radial strip exceeds actual donor");
            depth=delta/(radius+std::sqrt(square)); // Stable true-volume inversion.
        }
        if(positive)strip.lower.x=static_cast<double>(static_cast<long double>(column.upper.x)-depth);
        else strip.upper.x=static_cast<double>(static_cast<long double>(column.lower.x)+depth);
    }
    require(strip.lower.x>=column.lower.x&&strip.upper.x<=column.upper.x&&
        strip.lower.y>=column.lower.y&&strip.upper.y<=column.upper.y,"swept strip outside actual column");
    return strip;
}
bool localWholePhaseCovered(const QuadtreeMesh& mesh,const NativeLeafFace& face,const GeometricVofPlicCell& cell,
    long double total,bool positive,bool rz,unsigned phase,long double tolerance)
{
    // The complete phase polygon and consumed strip use the actual final
    // source plane and stored rectangle. No stock-size threshold proves this.
    const auto strip=localStripRectangle(face,localColumn(mesh,face,cell.rectangle),total,positive,rz);
    const auto hx=(static_cast<long double>(cell.rectangle.upper.x)-cell.rectangle.lower.x)/2;
    const auto hy=(static_cast<long double>(cell.rectangle.upper.y)-cell.rectangle.lower.y)/2;
    const Point origin{static_cast<long double>(cell.rectangle.lower.x)+hx,static_cast<long double>(cell.rectangle.lower.y)+hy};
    const auto nx=static_cast<long double>(cell.normal.x),ny=static_cast<long double>(cell.normal.y),norm=std::hypot(nx,ny);
    const Point normal{nx/norm,ny/norm};
    const std::vector<Point> rectangle{{-hx,-hy},{hx,-hy},{hx,hy},{-hx,hy}};
    const auto liquid=Cut::clipConvexPolygon(rectangle,normal,cell.offset);
    const auto gas=Cut::clipConvexPolygon(rectangle,{-normal.x,-normal.y},-cell.offset);
    const std::array<long double,2> measured{{Cut::polygonMeasureInFrame(liquid,origin,rz),Cut::polygonMeasureInFrame(gas,origin,rz)}};
    for(unsigned p=0;p<2;++p)if(!(measured[p]>0&&cell.inventory[p]>0&&
        std::abs(measured[p]-cell.inventory[p])/cell.inventory[p]<=tolerance))return false;
    const auto& polygon=phase==0?liquid:gas;
    const auto lx=static_cast<long double>(strip.lower.x)-origin.x,ux=static_cast<long double>(strip.upper.x)-origin.x;
    const auto ly=static_cast<long double>(strip.lower.y)-origin.y,uy=static_cast<long double>(strip.upper.y)-origin.y;
    if(polygon.size()<3)return false;
    for(const auto p:polygon)if(p.x<lx||p.x>ux||p.y<ly||p.y>uy)return false;
    // Four actual complement halfplanes independently confirm zero unswept
    // phase measure; the strict vertex containment above is tolerance-free.
    const std::array<Point,4> outsideNormals{{{1,0},{-1,0},{0,1},{0,-1}}};
    const std::array<long double,4> outsideOffsets{{lx,-ux,ly,-uy}};
    for(unsigned side=0;side<4;++side)
        require(Cut::polygonMeasureInFrame(Cut::clipConvexPolygon(polygon,outsideNormals[side],outsideOffsets[side]),origin,rz)==0,
            "whole-phase sweep certificate has positive actual unswept geometry");
    return true;
}
std::array<long double,2> localStripPhases(const QuadtreeMesh& mesh,const NativeLeafFace& face,
    const GeometricVofPlicCell& cell,long double total,bool positive,bool rz,
    long double tolerance,GeometricVofTransportReport& report,const LocalPin* pin=nullptr,long double maximumCourant=.5L)
{
    if(const auto far=localPinFarEndpoint(mesh,face,pin))return localPinnedPhases(mesh,face,cell,*pin,*far,total,positive,rz,maximumCourant,tolerance,report);
    const auto column=localColumn(mesh,face,cell.rectangle);
    const auto columnCapacity=Cut::capacity(column,rz),capacity=Cut::capacity(cell.rectangle,rz);
    require(total>0&&columnCapacity>0,"invalid local swept column");
    const auto strip=localStripRectangle(face,column,total,positive,rz);
    // Actual resolved planes use their geometric strip regardless of the
    // capillary/normal resolution mask. Unresolved or unrepresentable strips
    // retain the declared donor-fraction approximation without phase erasure.
    const auto fractionPhases=[&]{
        ++report.localUnresolvedDonors;
        const bool liquidMinor=cell.inventory[0]<=cell.inventory[1];
        const auto minor=total*((liquidMinor?cell.inventory[0]:cell.inventory[1])/capacity);
        return liquidMinor?std::array<long double,2>{minor,total-minor}:std::array<long double,2>{total-minor,minor};
    };
    if(cell.inventory[0]==0)return {0,total};
    if(cell.inventory[1]==0)return {total,0};
    if(!cell.resolved||strip.lower.x==strip.upper.x||strip.lower.y==strip.upper.y)return fractionPhases();
    const auto actual=Cut::capacity(strip,rz);
    const auto residual=std::abs(actual-total)/capacity;
    report.maximumLocalStripResidual=std::max(report.maximumLocalStripResidual,residual);
    require(residual<=tolerance,"local actual strip cannot represent requested volume");
    const auto nx=static_cast<long double>(cell.normal.x),ny=static_cast<long double>(cell.normal.y),norm=std::hypot(nx,ny);
    const auto cx=static_cast<long double>(cell.rectangle.lower.x)+(static_cast<long double>(cell.rectangle.upper.x)-cell.rectangle.lower.x)/2;
    const auto cy=static_cast<long double>(cell.rectangle.lower.y)+(static_cast<long double>(cell.rectangle.upper.y)-cell.rectangle.lower.y)/2;
    const auto sx=static_cast<long double>(strip.lower.x)+(static_cast<long double>(strip.upper.x)-strip.lower.x)/2;
    const auto sy=static_cast<long double>(strip.lower.y)+(static_cast<long double>(strip.upper.y)-strip.lower.y)/2;
    const auto offset=cell.offset+(nx*(cx-sx)+ny*(cy-sy))/norm;
    const auto l=Cut::volume(strip,cell.normal,offset,rz);
    const auto g=Cut::volume(strip,Vector2D(-cell.normal.x,-cell.normal.y),-offset,rz);
    ++report.localStripClips;
    require(l>=0&&g>=0&&std::isfinite(l+g)&&l+g>0&&std::abs(l+g-actual)<=tolerance*capacity,"invalid local two-phase strip coverage");
    // Normalize one shared strip measure to FT, preserving the measured minor.
    const bool liquidMinor=l<=g;const auto minor=total*((liquidMinor?l:g)/(l+g));
    return liquidMinor?std::array<long double,2>{minor,total-minor}:std::array<long double,2>{total-minor,minor};
}
}
GeometricVofLocalCourantRates GeometricVofTransport::localPlicCourantRates(
    const QuadtreeMesh& mesh,const GeometricVofPublication& source,const std::vector<double>& nativeNormal)
{return localPlicCourantRatesImpl(mesh,source,nativeNormal,nullptr);}
GeometricVofLocalCourantRates GeometricVofTransport::localPlicBoundaryCourantRates(
    const QuadtreeMesh& mesh,const GeometricVofPublication& source,
    const GeometricVofAdvectorSlab::Boundary& boundary,const std::vector<double>& nativeNormal)
{
    boundary.validateFor(mesh,source);
    if(boundary.allClosed())return localPlicCourantRates(mesh,source,nativeNormal);
    return localPlicCourantRatesImpl(mesh,source,nativeNormal,&boundary);
}
GeometricVofLocalCourantRates GeometricVofTransport::localPlicCourantRatesImpl(
    const QuadtreeMesh& mesh,const GeometricVofPublication& source,
    const std::vector<double>& nativeNormal,const GeometricVofAdvectorSlab::Boundary* boundary)
{
    source.validateFor(mesh,source.identity(),source.timeLayer());
    require(source.geometryKind()==GeometricVofGeometryKind::LocalPlic,"local rates require LocalPlic source");
    if(boundary)boundary->validateFor(mesh,source);
    const auto& domain=source.domain();const auto n=source.cells().size();const bool rz=domain.axisymmetric();
    require(nativeNormal.size()==static_cast<std::size_t>(domain.nativeFaceCount()),"local rate native normal dimensions differ");
    std::vector<Sum> incoming(n),outgoing(n);GeometricVofLocalCourantRates result;
    const auto pin=source.pinnedCorner();if(pin)pin->validateFor(mesh,source);
    for(int row=0;row<domain.nativeFaceCount();++row){
        require(std::isfinite(nativeNormal[row]),"nonfinite local rate native normal");
        if(boundary&&boundary->rules()[row]==GeometricVofAdvectorSlab::Boundary::Mode::Closed)
            require(nativeNormal[row]==0,"Closed boundary/axis rate normal must be exactly zero");
        const auto f=domain.face(mesh,row);
        if(f.negativeCell<0&&f.positiveCell<0){
            require(nativeNormal[row]==0,"excluded face rate normal must be exactly zero");continue;}
        const auto q=static_cast<long double>(f.area)*nativeNormal[row];
        require(std::isfinite(q),"unrepresentable local native volume rate");
        require(nativeNormal[row]==0||f.area==0||q!=0,"positive local native volume rate underflow");
        if(q==0)continue;const auto donor=q>0?f.negativeCell:f.positiveCell,receiver=q>0?f.positiveCell:f.negativeCell;
        if(!boundary)require(donor>=0&&receiver>=0,"local exterior transport is unsupported");
        if(donor>=0)outgoing[donor].add(std::abs(q));if(receiver>=0)incoming[receiver].add(std::abs(q));
        const auto adjacent=donor>=0?donor:receiver;
        auto rate=localColumnRatio(mesh,f,source.plicCells()[adjacent].rectangle,std::abs(q),rz);
        if(const auto far=localPinFarEndpoint(mesh,f,pin?&*pin:nullptr))
            rate=std::max(rate,std::abs(q)/localPinRateCapacity(mesh,f,source.plicCells()[adjacent].rectangle,*pin,*far,q>0,rz));
        require(std::isfinite(rate)&&rate>0,"unrepresentable local donor-column rate");
        result.maximumLocalColumn=std::max(result.maximumLocalColumn,rate);
    }
    for(std::size_t i=0;i<n;++i){if(domain.cellRowByIndex()[i]<0)continue;
        const auto capacity=source.cells()[i].capacity;
        require(std::isfinite(capacity)&&capacity>0,"invalid local rate cell capacity");
        const auto in=incoming[i].value/capacity,out=outgoing[i].value/capacity;
        require(std::isfinite(in)&&in>=0&&std::isfinite(out)&&out>=0,"unrepresentable local complete-cell rate");
        require((incoming[i].value==0||in>0)&&(outgoing[i].value==0||out>0),"positive local complete-cell rate underflow");
        result.maximumIncoming=std::max(result.maximumIncoming,in);result.maximumOutflow=std::max(result.maximumOutflow,out);
    }
    return result;
}
GeometricVofTransport GeometricVofTransport::prepareLocalPlicSweep(const QuadtreeMesh& mesh,
    const GeometricVofPublication& source,const GeometricVofAdvectorSlab& slab,
    long double rhoL,long double rhoG,const GeometricVofTransportConfig& config)
{
    require(!slab.boundary()||slab.boundary()->allClosed(),"closed local sweep rejects open boundary slab");
    return prepareLocalPlicSweepImpl(mesh,source,slab,rhoL,rhoG,config,false);
}
GeometricVofTransport GeometricVofTransport::prepareLocalPlicBoundarySweep(const QuadtreeMesh& mesh,
    const GeometricVofPublication& source,const GeometricVofAdvectorSlab& slab,
    long double rhoL,long double rhoG,const GeometricVofTransportConfig& config)
{
    if(!slab.boundary()||slab.boundary()->allClosed())return prepareLocalPlicSweep(mesh,source,slab,rhoL,rhoG,config);
    slab.boundary()->validateFor(mesh,source);
    const GeometricVofPublication::MeshValidationBatch batch(mesh,source.domain().meshIdentity());
    return prepareLocalPlicSweepImpl(mesh,source,slab,rhoL,rhoG,config,true,nullptr,nullptr,&batch);
}
GeometricVofTransport GeometricVofTransport::prepareLocalSmoothSweep(const QuadtreeMesh& mesh,
    const GeometricVofPublication& source,const VofSmoothGeometry& smooth,
    const GeometricVofAdvectorSlab& slab,long double rhoL,long double rhoG,
    const GeometricVofTransportConfig& config,const VofSmoothGeometryConfig& smoothConfig)
{
    require(!slab.boundary()||slab.boundary()->allClosed(),"smooth sweep requires a closed actual slab");
    return prepareLocalPlicSweepImpl(mesh,source,slab,rhoL,rhoG,config,false,&smooth,&smoothConfig);
}
GeometricVofTransport GeometricVofTransport::prepareLocalPlicSweepImpl(const QuadtreeMesh& mesh,
    const GeometricVofPublication& source,const GeometricVofAdvectorSlab& slab,
    long double rhoL,long double rhoG,const GeometricVofTransportConfig& config,bool open,
    const VofSmoothGeometry* smooth,const VofSmoothGeometryConfig* smoothConfig,
    const GeometricVofPublication::MeshValidationBatch* entryBatch)
{
    require(source.geometryKind()==GeometricVofGeometryKind::LocalPlic,"local sweep requires LocalPlic source");
    if(entryBatch)slab.validateForInBatch(mesh,source,*entryBatch);else slab.validateFor(mesh,source);
    // No callback or mesh mutation occurs before this sweep returns. The
    // four derived publications may reuse only this exact live binding/parent.
    const GeometricVofPublication::MeshValidationBatch batch(mesh,source.domain().meshIdentity());
    const GeometricVofPublication::ConstructionMeshReceipt constructionReceipt(mesh,source);
    const auto* boundary=open?slab.boundary():nullptr;
    require(!open||boundary,"open local sweep requires retained boundary owner");
    require(std::isfinite(rhoL)&&rhoL>0&&std::isfinite(rhoG)&&rhoG>0,"densities must be finite and positive");
    const auto tolerance=config.relativeConservationTolerance;
    require(std::isfinite(tolerance)&&tolerance>0&&tolerance<=1e-12L,"cannot loosen local conservation gate");
    require(std::isfinite(config.maximumLocalCourant)&&config.maximumLocalCourant>0&&config.maximumLocalCourant<=.5L,"local Courant must be in (0,.5]");
    const auto end=slab.endTime();const long double dt=end.end-end.begin;
    const auto& domain=source.domain();const auto n=source.cells().size();const bool rz=domain.axisymmetric();
    require((smooth!=nullptr)==(smoothConfig!=nullptr),"smooth sweep configuration/owner pair incomplete");
    if(smooth){
        require(!open&&!rz&&domain.uniformLeaves()&&domain.cellIndexByRow().size()==n,
            "smooth sweep requires actual closed complete uniform XY leaves");
        smooth->validateFor(mesh,source);
    }
    const auto pin=source.pinnedCorner();if(pin)pin->validateForInBatch(mesh,source,batch);
    if(config.geometry.localResolution==GeometricVofLocalResolution::NativeInterfaceSupport)
        require(!rz&&!open&&!pin&&!smooth,"native-interface resolution requires closed unpinned XY LocalPlic stocks");
    GeometricVofTransportReport report;
    std::vector<GeometricVofTransportFace> faces(domain.nativeFaceCount());
    std::vector<NativeLeafFace> native;native.reserve(faces.size());
    for(int row=0;row<domain.nativeFaceCount();++row){native.push_back(domain.face(mesh,row));auto& f=faces[row];
        f.nativeRow=row;f.negativeCell=native.back().negativeCell;f.positiveCell=native.back().positiveCell;}
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
    auto localFaceRoute=prepareLocalFaceRoute(native,n,smooth==nullptr,report.localFaceExecution);
#endif
    std::vector<Gross> boundaryMeasure(open?native.size():0);
    if(open)for(const auto& f:native)if((f.negativeCell>=0)!=(f.positiveCell>=0)){
        boundaryMeasure[f.nativeRow]=FaceTime::integrateNormal(slab.beginNormal()[f.nativeRow],slab.endNormal()[f.nativeRow],end.begin,end.end,f.area);
        faceTimeDiagnostics(boundaryMeasure[f.nativeRow],report);
    }
    // Require the actual endpoint snapshots to be closed and divergence-free,
    // not only a cancelling time average. Retain their numerical residues.
    for(bool last:{false,true}){
        std::vector<Sum> div(n);const auto& velocity=last?slab.endNormal():slab.beginNormal();
        for(const auto& f:native){const auto q=static_cast<long double>(f.area)*velocity[f.nativeRow]*dt;
            if(f.negativeCell>=0)div[f.negativeCell].add(q);if(f.positiveCell>=0)div[f.positiveCell].add(-q);}
        for(std::size_t i=0;i<n;++i){if(domain.cellRowByIndex()[i]<0)continue;
            const auto relative=std::abs(div[i].value)/source.cells()[i].capacity;
            report.maximumTotalDivergence=std::max(report.maximumTotalDivergence,relative);
            require(relative<=tolerance,"local sweep requires discrete incompressible endpoint snapshots");}
    }
    const auto extent=[&](bool half){
        const auto duration=half?dt/2:dt;std::vector<long double> total(native.size());
        std::vector<bool> counterflow(native.size());
        std::vector<std::array<long double,2>> counterGross(native.size());bool hasCounterflow=false;
        std::vector<Sum> incoming(n),outgoing(n),div(n);
        for(const auto& f:native){const auto a=static_cast<long double>(slab.beginNormal()[f.nativeRow]),b=static_cast<long double>(slab.endNormal()[f.nativeRow]);
            const auto mean=half?(3*a+b)/4:(a+b)/2;
            auto q=static_cast<long double>(f.area)*mean*duration;
            if(open&&(f.negativeCell<0||f.positiveCell<0)){const auto& measure=boundaryMeasure[f.nativeRow];
                q=half?measure.halfPositive-measure.halfNegative:measure.fullPositive-measure.fullNegative;}
            total[f.nativeRow]=q;
            if(open&&(f.negativeCell<0||f.positiveCell<0)){
                const auto& measure=boundaryMeasure[f.nativeRow];
                const auto positive=half?measure.halfPositive:measure.fullPositive;
                const auto negative=half?measure.halfNegative:measure.fullNegative;
                if(positive>0&&negative>0){
                    counterflow[f.nativeRow]=true;hasCounterflow=true;counterGross[f.nativeRow]={positive,negative};
                    const int inside=f.negativeCell>=0?f.negativeCell:f.positiveCell;
                    require(inside>=0,"counterflow boundary has no actual inside cell");
                    const auto out=f.negativeCell>=0?positive:negative,in=f.negativeCell>=0?negative:positive;
                    outgoing[inside].add(out);incoming[inside].add(in);
                    for(const auto amount:{positive,negative}){
                        const auto cfl=localColumnRatio(mesh,f,source.plicCells()[inside].rectangle,amount,rz);
                        report.maximumLocalColumnCfl=std::max(report.maximumLocalColumnCfl,cfl);
                        if(!(cfl<=config.maximumLocalCourant))throw GeometricVofCourantExceeded("counterflow boundary subface-column Courant exceeded");
                    }
                    if(f.negativeCell>=0)div[f.negativeCell].add(q);if(f.positiveCell>=0)div[f.positiveCell].add(-q);
                    continue;
                }
            }
            if(q==0)continue;const auto donor=q>0?f.negativeCell:f.positiveCell,receiver=q>0?f.positiveCell:f.negativeCell;
            if(!open)require(donor>=0&&receiver>=0,"local exterior transport is unsupported");
            if(donor>=0)outgoing[donor].add(std::abs(q));if(receiver>=0)incoming[receiver].add(std::abs(q));
            if(f.negativeCell>=0)div[f.negativeCell].add(q);if(f.positiveCell>=0)div[f.positiveCell].add(-q);
            const auto adjacent=donor>=0?donor:receiver;
            const auto cfl=localColumnRatio(mesh,f,source.plicCells()[adjacent].rectangle,std::abs(q),rz);
            report.maximumLocalColumnCfl=std::max(report.maximumLocalColumnCfl,cfl);
            if(!(cfl<=config.maximumLocalCourant))throw GeometricVofCourantExceeded("local donor subface-column Courant exceeded");
            if(const auto far=localPinFarEndpoint(mesh,f,pin?&*pin:nullptr)){
                const auto sweep=localPinnedSweep(mesh,f,source.plicCells()[adjacent].rectangle,*pin,*far,std::abs(q),q>0,rz,config.maximumLocalCourant,tolerance);
                report.maximumPinnedPeakCfl=std::max(report.maximumPinnedPeakCfl,sweep.peakCourant);}
        }
        for(std::size_t i=0;i<n;++i){if(domain.cellRowByIndex()[i]<0)continue;
            const auto c=source.cells()[i].capacity;
            const auto in=incoming[i].value/c,out=outgoing[i].value/c;
            report.maximumIncomingCfl=std::max(report.maximumIncomingCfl,in);report.maximumOutflowCfl=std::max(report.maximumOutflowCfl,out);
            if(!(in<=config.maximumLocalCourant&&out<=config.maximumLocalCourant))throw GeometricVofCourantExceeded("local complete-cell incoming/outgoing Courant exceeded");
            const auto r=std::abs(div[i].value)/c;report.maximumTotalDivergence=std::max(report.maximumTotalDivergence,r);
            report.maximumCompressionResidual=std::max(report.maximumCompressionResidual,r);
            require(r<=tolerance,"local total compression residue exceeds conservation budget");
        }
        std::vector<std::array<long double,2>> stocks;stocks.reserve(n);
        std::vector<unsigned> tracked;tracked.reserve(n);
        for(const auto& c:source.plicCells()){stocks.push_back(c.inventory);tracked.push_back(c.inventory[0]>c.inventory[1]?1:0);}
        std::array<FiniteVolume::Dyadic,2> compressionSum,roundingSum,netFluxRoundingSum;
        std::array<LocalDyadic,2> finalBefore,finalAfter;bool finalTotalsReady=false;
        auto defects=source.capacityDefects();
        require(defects.size()==n,"local closure witness dimensions differ");
        std::string pinTransferTrace;
        const auto publish=[&](Clock clock,GeometricVofPublication::CertifiedOperation operation){
            try{return GeometricVofPublication::prepareCertifiedLocalPlicInventory(mesh,source.domain(),clock,source,
                operation,stocks,defects,config.geometry,&constructionReceipt);}
            catch(const std::exception& error){
                if(pinTransferTrace.empty())throw;
                throw std::invalid_argument(std::string(error.what())+pinTransferTrace);
            }};
        auto scratch=source;const bool firstHorizontal=(end.step%2)==0;
        for(unsigned direction=0;direction<2;++direction){
            const bool horizontal=direction==0?firstHorizontal:!firstHorizontal;
            if(direction!=0)scratch=publish(source.timeLayer(),GeometricVofPublication::CertifiedOperation::DirectionalScratch);
            const auto directionalPin=scratch.pinnedCorner();if(directionalPin)directionalPin->validateForInBatch(mesh,scratch,batch);
            std::optional<VofSmoothGeometry> directionalSmooth;
            if(smooth)directionalSmooth=direction==0?*smooth:VofSmoothGeometry::prepare(mesh,scratch,*smoothConfig);
            std::vector<std::array<long double,2>> gross(native.size());
            std::vector<int> counterPhase(n,-1);
            std::vector<std::array<LocalDyadic,2>> demand(n);
            std::vector<LocalDyadic> demandedTotal(n);
            // -1: none; -2: multiple or a boundary gross outflow. Exact
            // saturation is eligible only on a genuine internal shared face.
            std::vector<std::array<int,2>> soleOutgoing(n,std::array<int,2>{{-1,-1}});
            const auto grossFace=[&](const NativeLeafFace& f,GeometricVofTransportReport& localReport){
                if(f.horizontal!=horizontal||(total[f.nativeRow]==0&&!counterflow[f.nativeRow]))return;
                const auto q=total[f.nativeRow];const auto donor=q>0?f.negativeCell:f.positiveCell;
                if(counterflow[f.nativeRow]){
                    const int inside=f.negativeCell>=0?f.negativeCell:f.positiveCell;
                    const auto mode=boundary->rules()[f.nativeRow];
                    require(mode==GeometricVofAdvectorSlab::Boundary::Mode::PureLiquidReservoir||
                        mode==GeometricVofAdvectorSlab::Boundary::Mode::PureGasReservoir,"counterflow requires an explicit pure reservoir");
                    const unsigned phase=mode==GeometricVofAdvectorSlab::Boundary::Mode::PureLiquidReservoir?0:1;
                    const auto& cell=scratch.plicCells()[inside];
                    // This is the producer's existing P0 discretization, selected
                    // before any cut/exception; it never relaxes resolved geometry.
                    if(!cell.resolved&&cell.inventory[0]>0&&cell.inventory[1]>0){
                        require(counterPhase[inside]<0||counterPhase[inside]==2,"counterflow donor discretizations disagree");
                        counterPhase[inside]=2;
                        const auto out=counterGross[f.nativeRow][f.negativeCell>=0?0:1];
                        gross[f.nativeRow]=localStripPhases(mesh,f,cell,out,f.negativeCell>=0,rz,tolerance,localReport);
                        for(unsigned p=0;p<2;++p){demand[inside][p].add(gross[f.nativeRow][p]);
                            if(gross[f.nativeRow][p]>0)soleOutgoing[inside][p]=-2;}
                        demandedTotal[inside].add(out);return;
                    }
                    require(cell.resolved&&cell.inventory[1-phase]==0&&cell.measured[1-phase]==0&&cell.measured[phase]>0,
                        "counterflow boundary requires a strict scratch-pure inside cell matching its reservoir");
                    require(counterPhase[inside]<0||counterPhase[inside]==static_cast<int>(phase),"counterflow reservoirs disagree on actual inside phase");
                    counterPhase[inside]=static_cast<int>(phase);
                    const auto out=counterGross[f.nativeRow][f.negativeCell>=0?0:1];
                    demand[inside][phase].add(out);demandedTotal[inside].add(out);
                    if(out>0)soleOutgoing[inside][phase]=-2;
                    return;
                }
                if(donor>=0){
                if(directionalSmooth){
                    const auto strip=directionalSmooth->measureNativeStrip(mesh,scratch,donor,
                        f.horizontal?1U:0U,q>0,std::abs(q));
                    gross[f.nativeRow]=strip.phase;++localReport.localStripClips;
                    localReport.maximumLocalStripResidual=std::max(localReport.maximumLocalStripResidual,strip.relativeCoverageResidual);
                    require(strip.relativeCoverageResidual<=tolerance,"smooth actual strip exceeds original coverage budget");
                }else gross[f.nativeRow]=localStripPhases(mesh,f,scratch.plicCells()[donor],std::abs(q),q>0,rz,tolerance,localReport,directionalPin?&*directionalPin:nullptr,config.maximumLocalCourant);
                for(unsigned phase=0;phase<2;++phase){
                    const auto value=gross[f.nativeRow][phase];demand[donor][phase].add(value);
                    if(value>0){auto& sole=soleOutgoing[donor][phase];
                        sole=f.negativeCell>=0&&f.positiveCell>=0&&sole==-1?f.nativeRow:-2;}
                }
                demandedTotal[donor].add(std::abs(q));
                }else {
                    require(open&&boundary,"external donor requires explicit reservoir");
                    const auto mode=boundary->rules()[f.nativeRow];
                    require(mode==GeometricVofAdvectorSlab::Boundary::Mode::PureLiquidReservoir||
                        mode==GeometricVofAdvectorSlab::Boundary::Mode::PureGasReservoir,"inflow requires explicit pure reservoir");
                    gross[f.nativeRow]=mode==GeometricVofAdvectorSlab::Boundary::Mode::PureLiquidReservoir?
                        std::array<long double,2>{std::abs(q),0}:std::array<long double,2>{0,std::abs(q)};
                }
            };
            std::map<int,std::array<long double,2>> wholePhaseFlux;
            const auto wholePhase=[&]{
            if(directionalPin&&!directionalSmooth)for(std::size_t i=0;i<n;++i){
                const auto phase=tracked[i],other=1U-phase;const auto& cell=scratch.plicCells()[i];
                const auto stock=stocks[i][phase],requested=demand[i][phase].rounded();const int row=soleOutgoing[i][phase];
                if(!cell.resolved||!(stock>0&&stock<=stocks[i][other]&&requested>0&&requested<stock)||
                    (stock-requested)/stock>tolerance||row<0)continue;
                const auto& f=native[row];const auto q=total[row];
                // The triangular pin profile has its own geometric producer;
                // this certificate concerns the original actual rectangular strip.
                if(localPinFarEndpoint(mesh,f,&*directionalPin)||stock>std::abs(q))continue;
                if(!localWholePhaseCovered(mesh,f,cell,std::abs(q),q>0,rz,phase,tolerance))continue;
                auto replacement=gross[row];const auto before=replacement;
                replacement[phase]=stock;replacement[other]=std::abs(q)-stock;
                require(std::isfinite(replacement[other])&&replacement[other]>=0&&
                    (before[other]==0||replacement[other]>0),"whole-phase shared flux cannot retain its other-phase measure");
                demand[i][phase].add(before[phase],-1);demand[i][phase].add(replacement[phase]);
                demand[i][other].add(before[other],-1);demand[i][other].add(replacement[other]);
                require(demand[i][phase].compare(LocalDyadic(stock))==0,
                    "whole-phase shared demand does not consume the exact independent stock");
                wholePhaseFlux.emplace(row,replacement);++report.localWholePhaseExhaustions;
                report.maximumWholePhaseExhaustionRelativeCorrection=std::max(
                    report.maximumWholePhaseExhaustionRelativeCorrection,(stock-requested)/stock);
                report.maximumFluxCorrection=std::max(report.maximumFluxCorrection,(stock-requested)/source.cells()[i].capacity);
            }
            };
            std::vector<long double> limiter(n,1);
            std::vector<unsigned> limitedPhase=tracked;
            std::vector<int> exactLimitedFace(n,-1);
            const auto limitCell=[&](std::size_t i){
                LocalDyadic available(stocks[i][0]);available.add(stocks[i][1]);
                require(demandedTotal[i].compare(available)<=0,"local total outflow exceeds independent donor stocks");
                int deficitPhase=-1;
                for(unsigned phase=0;phase<2;++phase)if(demand[i][phase].compare(LocalDyadic(stocks[i][phase]))>0){
                    require(deficitPhase<0,"both local phase outflows exceed donor stocks");deficitPhase=static_cast<int>(phase);}
                if(deficitPhase<0)return;
                const auto phase=static_cast<unsigned>(deficitPhase);const auto stock=stocks[i][phase];
                auto deficit=demand[i][phase];deficit.add(stock,-1);
                require(deficit.compare(FiniteVolume::Dyadic::product({tolerance,source.cells()[i].capacity}))<=0,
                    "local geometric outflow exceeds independent donor phase stock");
                const auto requested=demand[i][phase].rounded();
                require(std::isfinite(requested)&&requested>0,"invalid local phase demand");
                if(soleOutgoing[i][phase]>=0){
                    // A unique internal phase outflow has the exact feasible
                    // saturated solution F=stock before the shared ledger.
                    exactLimitedFace[i]=soleOutgoing[i][phase];limiter[i]=stock/requested;
                }else limiter[i]=stock==0?0:std::nextafter(std::min(1.L,stock/requested),0.L);
                require(std::isfinite(limiter[i])&&limiter[i]>=0&&limiter[i]<1,"invalid shared local donor limiter");
                limitedPhase[i]=phase;
            };
            std::vector<std::array<Sum,2>> net(n);std::vector<Sum> directionalTotal(n);
            std::vector<std::array<LocalDyadic,2>> exactNet(n),actualOut(n),actualIn(hasCounterflow?n:0);
            std::vector<LocalDyadic> exactTotal(n);
            const auto literalFaces=[&]{
            for(const auto& f:native)if(f.horizontal==horizontal){
                const auto q=total[f.nativeRow];auto& flux=half?faces[f.nativeRow].half:faces[f.nativeRow].full;
                flux.total=q;flux.positiveTotal=std::max(q,0.L);flux.negativeTotal=std::max(-q,0.L);
                if(counterflow[f.nativeRow]){
                    const int inside=f.negativeCell>=0?f.negativeCell:f.positiveCell;
                    if(counterPhase[inside]==2){
                        const bool insideNegative=f.negativeCell>=0;
                        const auto positive=counterGross[f.nativeRow][0],negative=counterGross[f.nativeRow][1];
                        const auto outTotal=insideNegative?positive:negative,inTotal=insideNegative?negative:positive;
                        auto outgoingPhases=gross[f.nativeRow];
                        const auto mode=boundary->rules()[f.nativeRow];
                        const unsigned incomingPhase=mode==GeometricVofAdvectorSlab::Boundary::Mode::PureLiquidReservoir?0:1;
                        std::array<long double,2> incomingPhases{};incomingPhases[incomingPhase]=inTotal;
                        flux.rawLiquid=insideNegative?outgoingPhases[0]-incomingPhases[0]:incomingPhases[0]-outgoingPhases[0];
                        if(limiter[inside]<1){const auto t=limitedPhase[inside];const auto original=outgoingPhases[t];
                            outgoingPhases[t]*=limiter[inside];outgoingPhases[1-t]=outTotal-outgoingPhases[t];
                            report.maximumFluxCorrection=std::max(report.maximumFluxCorrection,std::abs(outgoingPhases[t]-original)/source.cells()[inside].capacity);++report.limitedFaces;}
                        flux.limiter=limiter[inside];flux.positiveTotal=positive;flux.negativeTotal=negative;
                        const auto& positivePhases=insideNegative?outgoingPhases:incomingPhases;
                        const auto& negativePhases=insideNegative?incomingPhases:outgoingPhases;
                        flux.positiveLiquid=positivePhases[0];flux.positiveGas=positivePhases[1];
                        flux.negativeLiquid=negativePhases[0];flux.negativeGas=negativePhases[1];
                        flux.liquid=positivePhases[0]-negativePhases[0];flux.gas=positivePhases[1]-negativePhases[1];
                        flux.mass=rhoL*flux.liquid+rhoG*flux.gas;
                        require(std::isfinite(flux.mass),"unrepresentable counterflow P0 material net");
                        for(unsigned p=0;p<2;++p){
                            require(std::isfinite(outgoingPhases[p])&&outgoingPhases[p]>=0,"invalid counterflow P0 outgoing phase");
                            actualOut[inside][p].add(outgoingPhases[p]);actualIn[inside][p].add(incomingPhases[p]);
                            const auto signedPhase=p==0?flux.liquid:flux.gas;
                            if(f.negativeCell>=0){net[f.negativeCell][p].add(signedPhase);exactNet[f.negativeCell][p].add(signedPhase);}
                            if(f.positiveCell>=0){net[f.positiveCell][p].add(-signedPhase);exactNet[f.positiveCell][p].add(signedPhase,-1);}
                            FiniteVolume::Dyadic error(signedPhase);error.add(positivePhases[p],-1);error.add(negativePhases[p]);
                            netFluxRoundingSum[p].add(error,insideNegative?-1:1);
                        }
                        if(f.negativeCell>=0){directionalTotal[f.negativeCell].add(q);exactTotal[f.negativeCell].add(q);}
                        if(f.positiveCell>=0){directionalTotal[f.positiveCell].add(-q);exactTotal[f.positiveCell].add(q,-1);}
                        continue;
                    }
                    const unsigned phase=static_cast<unsigned>(counterPhase[inside]);
                    require(phase<2&&limiter[inside]==1,"counterflow pure donor requires unmodified gross availability");
                    const auto positive=counterGross[f.nativeRow][0],negative=counterGross[f.nativeRow][1];
                    flux.positiveTotal=positive;flux.negativeTotal=negative;
                    if(phase==0){flux.positiveLiquid=positive;flux.negativeLiquid=negative;flux.liquid=q;flux.rawLiquid=q;}
                    else{flux.positiveGas=positive;flux.negativeGas=negative;flux.gas=q;}
                    flux.mass=rhoL*flux.liquid+rhoG*flux.gas;
                    require(std::isfinite(flux.mass),"unrepresentable counterflow material net");
                    actualOut[inside][phase].add(f.negativeCell>=0?positive:negative);
                    actualIn[inside][phase].add(f.negativeCell>=0?negative:positive);
                    if(f.negativeCell>=0){net[f.negativeCell][phase].add(q);exactNet[f.negativeCell][phase].add(q);}
                    if(f.positiveCell>=0){net[f.positiveCell][phase].add(-q);exactNet[f.positiveCell][phase].add(q,-1);}
                    FiniteVolume::Dyadic error(q);error.add(positive,-1);error.add(negative);
                    netFluxRoundingSum[phase].add(error,f.negativeCell>=0?-1:1);
                    if(f.negativeCell>=0){directionalTotal[f.negativeCell].add(q);exactTotal[f.negativeCell].add(q);}
                    if(f.positiveCell>=0){directionalTotal[f.positiveCell].add(-q);exactTotal[f.positiveCell].add(q,-1);}
                    continue;
                }
                if(q==0)continue;const auto donor=q>0?f.negativeCell:f.positiveCell;
                const auto t=donor>=0?limitedPhase[donor]:0;auto phases=gross[f.nativeRow];flux.rawLiquid=(q>0?1:-1)*phases[0];
                if(const auto receipt=wholePhaseFlux.find(f.nativeRow);receipt!=wholePhaseFlux.end()){
                    const auto sign=q>0?1.L:-1.L;
                    flux.exhaustionLiquidCorrection=sign*(receipt->second[0]-phases[0]);
                    flux.exhaustionGasCorrection=sign*(receipt->second[1]-phases[1]);phases=receipt->second;
                }
                auto faceLimiter=donor>=0?limiter[donor]:1;
                if(donor>=0&&limiter[donor]<1){const auto original=phases[t];
                    if(exactLimitedFace[donor]==f.nativeRow)phases[t]=stocks[donor][t];
                    else {require(exactLimitedFace[donor]<0||original==0,"unique phase outflow has another nonzero face");phases[t]*=limiter[donor];}
                    phases[1-t]=std::abs(q)-phases[t];
                    if(exactLimitedFace[donor]>=0)faceLimiter=original==0?1:phases[t]/original;
                    report.maximumFluxCorrection=std::max(report.maximumFluxCorrection,std::abs(phases[t]-original)/source.cells()[donor].capacity);++report.limitedFaces;}
                flux.limiter=faceLimiter;flux.liquid=(q>0?1:-1)*phases[0];flux.gas=(q>0?1:-1)*phases[1];flux.mass=rhoL*flux.liquid+rhoG*flux.gas;
                require(std::isfinite(flux.mass),"unrepresentable local material flux");
                if(q>0){flux.positiveLiquid=phases[0];flux.positiveGas=phases[1];}
                else {flux.negativeLiquid=phases[0];flux.negativeGas=phases[1];}
                for(unsigned phase=0;phase<2;++phase){
                    require(std::isfinite(phases[phase])&&phases[phase]>=0,"invalid shared local phase flux");
                    const auto signedPhase=phase==0?flux.liquid:flux.gas;
                    if(donor>=0)actualOut[donor][phase].add(phases[phase]);
                    const int receiver=q>0?f.positiveCell:f.negativeCell;
                    if(receiver>=0&&counterPhase[receiver]>=0)actualIn[receiver][phase].add(phases[phase]);
                    if(f.negativeCell>=0)net[f.negativeCell][phase].add(signedPhase);if(f.positiveCell>=0)net[f.positiveCell][phase].add(-signedPhase);
                    if(f.negativeCell>=0)exactNet[f.negativeCell][phase].add(signedPhase);if(f.positiveCell>=0)exactNet[f.positiveCell][phase].add(signedPhase,-1);
                }
                if(f.negativeCell>=0)directionalTotal[f.negativeCell].add(q);if(f.positiveCell>=0)directionalTotal[f.positiveCell].add(-q);
                if(f.negativeCell>=0)exactTotal[f.negativeCell].add(q);if(f.positiveCell>=0)exactTotal[f.positiveCell].add(q,-1);
            }
            };
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
            const auto foldObservations=[&]{
                        for(const auto& f:native)if(f.horizontal==horizontal){
                            const auto& observation=localFaceRoute->observations[f.nativeRow];
                            if(observation.limited){report.maximumFluxCorrection=std::max(report.maximumFluxCorrection,observation.correction);++report.limitedFaces;}
                            if(!counterflow[f.nativeRow])continue;
                            const auto& flux=half?faces[f.nativeRow].half:faces[f.nativeRow].full;
                            const int inside=f.negativeCell>=0?f.negativeCell:f.positiveCell;
                            if(counterPhase[inside]==2){
                                for(unsigned p=0;p<2;++p){const auto signedPhase=p==0?flux.liquid:flux.gas;
                                    const auto positive=p==0?flux.positiveLiquid:flux.positiveGas;
                                    const auto negative=p==0?flux.negativeLiquid:flux.negativeGas;
                                    FiniteVolume::Dyadic error(signedPhase);error.add(positive,-1);error.add(negative);
                                    netFluxRoundingSum[p].add(error,f.negativeCell>=0?-1:1);}
                            }else{
                                const auto phase=static_cast<unsigned>(counterPhase[inside]);
                                FiniteVolume::Dyadic error(total[f.nativeRow]);error.add(counterGross[f.nativeRow][0],-1);error.add(counterGross[f.nativeRow][1]);
                                netFluxRoundingSum[phase].add(error,f.negativeCell>=0?-1:1);
                            }
                        }
            };
#endif
            const auto tracePin=[&]{
            if(directionalPin){
                const auto wallCell=[&](int row){const auto f=domain.face(mesh,row);return f.negativeCell>=0?f.negativeCell:f.positiveCell;};
                const int wet=wallCell(directionalPin->wetNativeRow()),dry=wallCell(directionalPin->dryNativeRow());
                for(const auto& f:native)if(f.horizontal==horizontal&&total[f.nativeRow]!=0){
                    const auto& flux=half?faces[f.nativeRow].half:faces[f.nativeRow].full;
                    for(unsigned forbidden=0;forbidden<2;++forbidden){
                        const int receiver=forbidden==0?dry:wet;
                        if(receiver!=f.negativeCell&&receiver!=f.positiveCell)continue;
                        const auto phase=forbidden==0?flux.liquid:flux.gas;
                        const auto incoming=receiver==f.positiveCell?phase:-phase;if(!(incoming>0))continue;
                        const int donor=total[f.nativeRow]>0?f.negativeCell:f.positiveCell;
                        if(donor<0)continue;const auto& cell=scratch.plicCells()[donor];const auto e=localFaceEndpoints(mesh,f);
                        std::ostringstream detail;detail.precision(17);
                        detail<<" [pin_phase_entry extent="<<(half?"half":"full")<<" direction="<<direction
                            <<" horizontal="<<horizontal<<" face="<<f.nativeRow<<" receiver="<<receiver<<" forbidden_phase="<<forbidden
                            <<" donor="<<donor<<" pin_profile="<<localPinFarEndpoint(mesh,f,&*directionalPin).has_value()
                            <<" FT="<<total[f.nativeRow]<<" raw_L="<<gross[f.nativeRow][0]<<" raw_G="<<gross[f.nativeRow][1]
                            <<" accepted_L="<<flux.liquid<<" accepted_G="<<flux.gas<<" limiter="<<flux.limiter
                            <<" donor_L="<<cell.inventory[0]<<" donor_G="<<cell.inventory[1]<<" donor_resolved="<<cell.resolved
                            <<" normal_x="<<cell.normal.x<<" normal_y="<<cell.normal.y<<" offset="<<cell.offset
                            <<" lower_x="<<cell.rectangle.lower.x<<" lower_y="<<cell.rectangle.lower.y
                            <<" upper_x="<<cell.rectangle.upper.x<<" upper_y="<<cell.rectangle.upper.y
                            <<" face_x0="<<e[0].x<<" face_y0="<<e[0].y<<" face_x1="<<e[1].x<<" face_y1="<<e[1].y<<"]";
                        pinTransferTrace+=detail.str();
                    }
                }
            }
            };
            const auto updateCell=[&](std::size_t i,GeometricVofTransportReport& localReport,
                std::array<LocalDyadic,2>& localCompression,std::array<LocalDyadic,2>& localRounding){
                if(domain.cellRowByIndex()[i]<0)return;
                const auto t=tracked[i];const auto c=source.cells()[i].capacity;
                if(counterPhase[i]>=0&&counterPhase[i]<2)require(actualIn[i][1-static_cast<unsigned>(counterPhase[i])].compare(LocalDyadic{})==0,
                    "counterflow scratch-pure cell receives another phase in the actual direction");
                const auto prior=stocks[i];auto predicted=defects[i];
                // The original WY directional source cancels the geometric
                // total divergence in the joint ledger. Retain the exact
                // exported phase-pair defect and actual final RN errors.
                predicted.add(exactNet[i][0],-1);predicted.add(exactNet[i][1],-1);predicted.add(exactTotal[i]);
                for(unsigned phase=0;phase<2;++phase){
                    require(actualOut[i][phase].compare(LocalDyadic(prior[phase]))<=0,
                        "accepted shared local phase outflow exceeds actual donor stock");
                    const auto compression=phase==t?LocalDyadic{}:exactTotal[i];
                    const auto updated=FiniteVolume::phaseUpdate(prior[phase],exactNet[i][phase],compression);
                    ++localReport.localPhaseUpdateCalls;++localReport.localPhaseUpdateFallback;
                    require(std::isfinite(updated.value)&&updated.value>=0,"negative or nonfinite independent local stock; no clipping");
                    stocks[i][phase]=updated.value;predicted.add(updated.error);
                    if(open){localCompression[phase].add(compression);localRounding[phase].add(updated.error);}
                    const auto cc=phase==t?0.L:1.L;
                    const auto expected=prior[phase]-net[i][phase].value+cc*directionalTotal[i].value;
                    const auto residual=std::abs(stocks[i][phase]-expected)/c;
                    localReport.maximumCapacityReconciliation=std::max(localReport.maximumCapacityReconciliation,residual);
                    require(residual<=tolerance,"local directional phase/source ledger mismatch");
                }
                defects[i]=std::move(predicted);
            };
            bool indexedFacesUsed=false;
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
            if(localFaceRoute){
                const int entryErrno=errno;LocalFaceFloatingState entry;const bool qualified=localFaceQualified(entry);
                if(entry.valid)entry.restore();errno=entryErrno;
                if(qualified){
                    const auto reportBefore=report;const auto traceSize=pinTransferTrace.size();
                    std::optional<std::array<LocalDyadic,6>> globalBefore;
                    try{
                        globalBefore.emplace(std::array<LocalDyadic,6>{{compressionSum[0],compressionSum[1],
                            roundingSum[0],roundingSum[1],netFluxRoundingSum[0],netFluxRoundingSum[1]}});
                        LocalFaceCleanup completed;
                        struct Totals {GeometricVofTransportReport report;std::array<LocalDyadic,2> compression,rounding,before,after;};
                        LocalFaceBudget budget;budget.used=localFaceRoute->bytes;
                        budget.add(static_cast<std::size_t>(localFaceRoute->team),sizeof(Totals));
                        std::vector<Totals> totals(static_cast<std::size_t>(localFaceRoute->team));
                        const LocalFaceInputs input{native,total,counterflow,counterGross,gross,wholePhaseFlux,counterPhase,
                            limiter,limitedPhase,exactLimitedFace,stocks,boundary,source,faces,horizontal,half,rhoL,rhoG};
                        LocalFaceAccumulators output{net,directionalTotal,exactNet,actualOut,actualIn,exactTotal};
                        entry.restore();errno=entryErrno;unsigned workers=0;
                        if(!localFaceTeam(localFaceRoute->team,[&](bool ready,std::atomic<int>& failed){
                            auto& local=totals[static_cast<std::size_t>(omp_get_thread_num())];
                            const auto rows=[&](std::size_t count,const auto& kernel){
#pragma omp for schedule(static)
                                for(std::ptrdiff_t row=0;row<static_cast<std::ptrdiff_t>(count);++row){
                                    if(!ready||failed.load(std::memory_order_relaxed))continue;
                                    try{kernel(static_cast<std::size_t>(row));}catch(...){failed.store(1,std::memory_order_relaxed);}
                                }
                            };
                            const auto owners=[&](const auto& kernel){
#pragma omp for schedule(dynamic,16)
                                for(std::ptrdiff_t row=0;row<static_cast<std::ptrdiff_t>(n);++row){
                                    if(!ready||failed.load(std::memory_order_relaxed))continue;
                                    try{kernel(static_cast<std::size_t>(row));}catch(...){failed.store(1,std::memory_order_relaxed);}
                                }
                            };
                            const auto ordered=[&](const auto& kernel){
#pragma omp master
                                {if(!failed.load(std::memory_order_relaxed))try{kernel();}catch(...){failed.store(1,std::memory_order_relaxed);}}
#pragma omp barrier
                            };
                            owners([&](std::size_t cell){
                                // Preserve native-face order within the unique
                                // donor/inside owner, including zero-net gross.
                                for(auto k=localFaceRoute->offsets[cell];k<localFaceRoute->offsets[cell+1];++k){
                                    const auto& f=native[localFaceRoute->incidence[k].row];const auto q=total[f.nativeRow];
                                    const int donor=q>0?f.negativeCell:f.positiveCell;
                                    const int owner=counterflow[f.nativeRow]?(f.negativeCell>=0?f.negativeCell:f.positiveCell)
                                        :donor>=0?donor:(q>0?f.positiveCell:f.negativeCell);
                                    if(owner==static_cast<int>(cell))grossFace(f,local.report);
                                }
                            });
                            ordered(wholePhase);rows(n,limitCell);
                            rows(native.size(),[&](std::size_t row){localFaceCompute(input,row,localFaceRoute->observations[row]);});
                            rows(n,[&](std::size_t cell){localFaceGather(input,*localFaceRoute,output,cell);});
                            ordered([&]{foldObservations();tracePin();});
                            rows(n,[&](std::size_t cell){
                                updateCell(cell,local.report,local.compression,local.rounding);
                                if(open&&direction==1)for(unsigned phase=0;phase<2;++phase){
                                    local.before[phase].add(source.plicCells()[cell].inventory[phase]);local.after[phase].add(stocks[cell][phase]);}
                            });
                        },workers))throw LocalFaceOptionalFailure{};
                        for(const auto& local:totals){
                            report.localStripClips+=local.report.localStripClips;report.localUnresolvedDonors+=local.report.localUnresolvedDonors;
                            report.localPinnedSweeps+=local.report.localPinnedSweeps;report.localPhaseUpdateCalls+=local.report.localPhaseUpdateCalls;
                            report.localPhaseUpdateFallback+=local.report.localPhaseUpdateFallback;
                            report.maximumLocalStripResidual=std::max(report.maximumLocalStripResidual,local.report.maximumLocalStripResidual);
                            report.maximumPinnedPeakCfl=std::max(report.maximumPinnedPeakCfl,local.report.maximumPinnedPeakCfl);
                            report.maximumCapacityReconciliation=std::max(report.maximumCapacityReconciliation,local.report.maximumCapacityReconciliation);
                            if(open)for(unsigned phase=0;phase<2;++phase){
                                compressionSum[phase].add(local.compression[phase]);roundingSum[phase].add(local.rounding[phase]);
                                if(direction==1){finalBefore[phase].add(local.before[phase]);finalAfter[phase].add(local.after[phase]);}}
                        }
                        if(open&&direction==1)finalTotalsReady=true;
                        const int finalErrno=errno;LocalFaceFloatingState finalState;
                        if(!finalState.capture())throw LocalFaceOptionalFailure{};
                        completed.finalState=finalState;completed.error=finalErrno;completed.armed=true;
                        globalBefore.reset();finalState.restore();errno=finalErrno;
                        auto& execution=report.localFaceExecution;
                        execution.faceValuesParallelUsed|=workers>1;execution.cellGatherParallelUsed|=workers>1;
                        execution.maximumWorkers=std::max(execution.maximumWorkers,workers);execution.scratchBytes=std::max(execution.scratchBytes,budget.used);
                        ++execution.directionsUsed;execution.faceEvaluations+=localFaceRoute->orientedFaces[horizontal?1:0];execution.cellGathers+=n;
                        indexedFacesUsed=true;
                    }catch(...){
                        LocalFaceCleanup replayRestore;replayRestore.finalState=entry;replayRestore.error=entryErrno;replayRestore.armed=true;
                        // Scratch is the immutable direction-entry stocks and
                        // certified defects; only failure copies them back.
                        for(std::size_t i=0;i<n;++i){stocks[i]=scratch.plicCells()[i].inventory;defects[i]=scratch.capacityDefects()[i];}
                        if(globalBefore){using std::swap;for(unsigned phase=0;phase<2;++phase){
                            swap(compressionSum[phase],(*globalBefore)[phase]);swap(roundingSum[phase],(*globalBefore)[2+phase]);
                            swap(netFluxRoundingSum[phase],(*globalBefore)[4+phase]);}}
                        finalBefore={};finalAfter={};finalTotalsReady=false;report=reportBefore;pinTransferTrace.resize(traceSize);wholePhaseFlux.clear();
                        for(auto& value:gross)value={};std::fill(counterPhase.begin(),counterPhase.end(),-1);
                        for(auto& cell:demand)for(auto& phase:cell)phase=LocalDyadic{};for(auto& cell:demandedTotal)cell=LocalDyadic{};
                        std::fill(soleOutgoing.begin(),soleOutgoing.end(),std::array<int,2>{{-1,-1}});
                        std::fill(limiter.begin(),limiter.end(),1.L);limitedPhase=tracked;std::fill(exactLimitedFace.begin(),exactLimitedFace.end(),-1);
                        for(const auto& f:native)if(f.horizontal==horizontal)(half?faces[f.nativeRow].half:faces[f.nativeRow].full)={};
                        for(auto& cell:net)cell={};for(auto& cell:directionalTotal)cell={};
                        for(auto& cell:exactNet)for(auto& phase:cell)phase=LocalDyadic{};
                        for(auto& cell:actualOut)for(auto& phase:cell)phase=LocalDyadic{};
                        for(auto& cell:actualIn)for(auto& phase:cell)phase=LocalDyadic{};for(auto& cell:exactTotal)cell=LocalDyadic{};
                        globalBefore.reset();entry.restore();errno=entryErrno;report.localFaceExecution.fallback=true;localFaceRoute.reset();
                    }
                }else {report.localFaceExecution.fallback=true;localFaceRoute.reset();}
            }
#endif
            if(!indexedFacesUsed){
                for(const auto& f:native)grossFace(f,report);
                wholePhase();for(std::size_t i=0;i<n;++i)limitCell(i);literalFaces();tracePin();
                for(std::size_t i=0;i<n;++i)updateCell(i,report,compressionSum,roundingSum);
            }
            ++report.localDirectionalSweeps;
        }
        Sum oldL,newL,oldG,newG,capacitySum;
        std::vector<std::array<Sum,2>> net(n);std::vector<Sum> totalNet(n);
        for(const auto& f:native){const auto& flux=half?faces[f.nativeRow].half:faces[f.nativeRow].full;
            for(int cell:{f.negativeCell,f.positiveCell})if(cell>=0){const auto sign=cell==f.negativeCell?1.L:-1.L;
                net[cell][0].add(sign*flux.liquid);net[cell][1].add(sign*flux.gas);totalNet[cell].add(sign*flux.total);}}
        for(std::size_t i=0;i<n;++i){if(domain.cellRowByIndex()[i]<0)continue;
            const auto& old=source.plicCells()[i].inventory;const auto c=source.cells()[i].capacity;
            oldL.add(old[0]);newL.add(stocks[i][0]);oldG.add(old[1]);newG.add(stocks[i][1]);capacitySum.add(c);
            for(unsigned phase=0;phase<2;++phase){const auto cc=phase==tracked[i]?0.L:1.L;
                const auto expected=old[phase]-net[i][phase].value+cc*totalNet[i].value;
                const auto residual=std::abs(stocks[i][phase]-expected)/c;
                report.maximumInventoryResidual=std::max(report.maximumInventoryResidual,residual);
                require(residual<=tolerance,"local complete phase/flux/compression ledger mismatch");}
        }
        long double relative;
        if(!open){
            relative=std::max(std::abs(oldL.value-newL.value),std::abs(oldG.value-newG.value))/capacitySum.value;
            require(relative<=tolerance,"local closed-domain phase conservation failed");
        }else {
            auto& ledger=half?report.boundaryHalf:report.boundaryFull;ledger.enabled=true;
            ledger.sourceIdentity=source.identity();ledger.sourceTime=source.timeLayer();ledger.targetTime=end;
            if(half)ledger.targetTime.position=Clock::Position::Midpoint;
            std::array<FiniteVolume::Dyadic,2> before,after,in,out;
            if(finalTotalsReady){before=std::move(finalBefore);after=std::move(finalAfter);}
            else for(std::size_t i=0;i<n;++i)for(unsigned phase=0;phase<2;++phase){
                before[phase].add(source.plicCells()[i].inventory[phase]);after[phase].add(stocks[i][phase]);}
            for(const auto& f:native)if((f.negativeCell>=0)!=(f.positiveCell>=0)){
                const auto& flux=half?faces[f.nativeRow].half:faces[f.nativeRow].full;
                const auto sign=f.negativeCell>=0?1.L:-1.L;
                for(unsigned phase=0;phase<2;++phase){
                    const auto positive=phase==0?flux.positiveLiquid:flux.positiveGas;
                    const auto negative=phase==0?flux.negativeLiquid:flux.negativeGas;
                    out[phase].add(sign>0?positive:negative);in[phase].add(sign>0?negative:positive);
                }
            }
            FiniteVolume::Dyadic mass;long double massScale=0;
            for(unsigned phase=0;phase<2;++phase){
                ledger.source[phase]=before[phase].rounded();ledger.target[phase]=after[phase].rounded();
                ledger.inflow[phase]=in[phase].rounded();ledger.outflow[phase]=out[phase].rounded();
                ledger.compression[phase]=compressionSum[phase].rounded();ledger.roundingError[phase]=roundingSum[phase].rounded();
                ledger.netFluxRounding[phase]=netFluxRoundingSum[phase].rounded();
                auto physical=after[phase];physical.add(before[phase],-1);physical.add(out[phase]);physical.add(in[phase],-1);
                ledger.physicalResidual[phase]=physical.rounded();
                const auto phaseRelative=std::abs(ledger.physicalResidual[phase])/capacitySum.value;
                ledger.physicalRelativeConservation=std::max(ledger.physicalRelativeConservation,phaseRelative);
                require(physical.absolute().compare(FiniteVolume::Dyadic::product({tolerance,capacitySum.value}))<=0,
                    "local open-domain physical phase boundary conservation failed");
                if(before[phase].compare(FiniteVolume::Dyadic{})==0&&in[phase].compare(FiniteVolume::Dyadic{})==0)
                    require(after[phase].compare(FiniteVolume::Dyadic{})==0&&out[phase].compare(FiniteVolume::Dyadic{})==0,
                        "zero source/inflow phase must remain exactly zero");
                auto algebraic=physical;algebraic.add(compressionSum[phase],-1);algebraic.add(roundingSum[phase],-1);
                algebraic.add(netFluxRoundingSum[phase],-1);
                ledger.algebraicResidual[phase]=algebraic.rounded();
                ledger.algebraicRelativeConservation=std::max(ledger.algebraicRelativeConservation,
                    std::abs(ledger.algebraicResidual[phase])/capacitySum.value);
                require(algebraic.absolute().compare(FiniteVolume::Dyadic::product({tolerance,capacitySum.value}))<=0,
                    "local open-domain WY/RN boundary algebraic ledger mismatch");
                const auto density=phase==0?rhoL:rhoG;
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE) && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST)
                if constexpr(sizeof(long double)==8&&sizeof(double)==8
                    &&std::numeric_limits<long double>::radix==2&&std::numeric_limits<long double>::digits==53
                    &&std::numeric_limits<long double>::min_exponent==-1021&&std::numeric_limits<long double>::max_exponent==1024
                    &&std::numeric_limits<long double>::has_denorm==std::denorm_present){
                    // The retained exact phase residual already contains these
                    // same stocks and oriented gross boundary integrals.
                    mass.add(physical.multiplied(density));
                }else
#endif
                {
                for(std::size_t i=0;i<n;++i){
                    mass.add(FiniteVolume::Dyadic::product({density,stocks[i][phase]}));
                    mass.add(FiniteVolume::Dyadic::product({density,source.plicCells()[i].inventory[phase]}),-1);
                }
                for(const auto& f:native)if((f.negativeCell>=0)!=(f.positiveCell>=0)){
                    const auto& flux=half?faces[f.nativeRow].half:faces[f.nativeRow].full;
                    const auto positive=phase==0?flux.positiveLiquid:flux.positiveGas;
                    const auto negative=phase==0?flux.negativeLiquid:flux.negativeGas;
                    const bool forwardOut=f.negativeCell>=0;
                    mass.add(FiniteVolume::Dyadic::product({density,forwardOut?positive:negative}));
                    mass.add(FiniteVolume::Dyadic::product({density,forwardOut?negative:positive}),-1);
                }
                }
                massScale+=density*(ledger.source[phase]+ledger.inflow[phase]);
            }
            require(std::isfinite(massScale)&&massScale>0,"invalid open boundary physical mass scale");
            ledger.massResidual=mass.rounded();ledger.massRelativeConservation=std::abs(ledger.massResidual)/massScale;
            require(ledger.massRelativeConservation<=tolerance,"local open-domain physical mass boundary conservation failed");
            relative=ledger.physicalRelativeConservation;
        }
        if(half)report.halfRelativeConservation=relative;else report.fullRelativeConservation=relative;
        auto clock=end;if(half)clock.position=Clock::Position::Midpoint;
        return publish(clock,half?GeometricVofPublication::CertifiedOperation::HalfTransport:
            GeometricVofPublication::CertifiedOperation::FullTransport);
    };
    auto midpoint=extent(true),endpoint=extent(false);
    auto data=std::make_shared<Data>(source,slab,rhoL,rhoG,config,std::move(midpoint),std::move(endpoint),std::move(faces),report,nullptr,nullptr);
    if(smooth){
        data->smoothSource=*smooth;
        data->smoothMidpoint=VofSmoothGeometry::prepare(mesh,data->midpoint,*smoothConfig);
        data->smoothEndpoint=VofSmoothGeometry::prepare(mesh,data->endpoint,*smoothConfig);
    }
    return GeometricVofTransport(std::move(data));
}
GeometricVofIdentity GeometricVofTransport::sourceIdentity() const{require(bool(data_),"consumed transport");return data_->source;}
Clock GeometricVofTransport::sourceTime() const{require(bool(data_),"consumed transport");return data_->sourceTime;}
Clock GeometricVofTransport::endTime() const{require(bool(data_),"consumed transport");return data_->endTime;}
const GeometricVofPublication& GeometricVofTransport::midpoint() const{require(bool(data_),"consumed transport");return data_->midpoint;}
const GeometricVofPublication& GeometricVofTransport::endpoint() const{require(bool(data_),"consumed transport");return data_->endpoint;}
const std::vector<GeometricVofTransportFace>& GeometricVofTransport::faces() const{require(bool(data_),"consumed transport");return data_->faces;}
const GeometricVofTransportReport& GeometricVofTransport::report() const{require(bool(data_),"consumed transport");return data_->report;}
const GeometricVofPublication& GeometricVofTransport::source() const
{require(bool(data_),"consumed transport");return data_->original;}
const GeometricVofAdvectorSlab& GeometricVofTransport::advector() const
{require(bool(data_),"consumed transport");return data_->slab;}
const VofSmoothGeometry* GeometricVofTransport::smoothSource() const{require(bool(data_),"consumed transport");return data_->smoothSource?&*data_->smoothSource:nullptr;}
const VofSmoothGeometry* GeometricVofTransport::smoothMidpoint() const{require(bool(data_),"consumed transport");return data_->smoothMidpoint?&*data_->smoothMidpoint:nullptr;}
const VofSmoothGeometry* GeometricVofTransport::smoothEndpoint() const{require(bool(data_),"consumed transport");return data_->smoothEndpoint?&*data_->smoothEndpoint:nullptr;}
long double GeometricVofTransport::liquidDensity() const
{require(bool(data_),"consumed transport");return data_->rhoL;}
long double GeometricVofTransport::gasDensity() const
{require(bool(data_),"consumed transport");return data_->rhoG;}
bool GeometricVofTransport::sameInstance(const GeometricVofTransport& other) const noexcept
{return data_&&data_==other.data_;}
void GeometricVofTransport::validateFor(const QuadtreeMesh& mesh,const GeometricVofPublication& old,
    const GeometricVofAdvectorSlab& slab,long double rhoL,long double rhoG) const
{
    require(bool(data_),"consumed transport");old.validateFor(mesh,data_->source,data_->sourceTime);
    // All following operations are private, read-only and callback-free. No
    // witness escapes this invocation or suppresses a later public scan.
    const GeometricVofPublication::MeshValidationBatch batch(mesh,old.domain().meshIdentity());
    slab.validateForInBatch(mesh,old,batch);require(data_->slab.sameInstance(slab),"different advector slab instance");
    require(rhoL==data_->rhoL&&rhoG==data_->rhoG,"different material densities");
    if(old.geometryKind()==GeometricVofGeometryKind::CompatibleFan)
        require(bool(data_->halfRestriction)&&bool(data_->fullRestriction),"missing exact accepted native restriction");
}
namespace
{
// LocalPlic material consumers use the ALREADY ACCEPTED split flux as a P0
// numerical face trace. They do not perform a second PLIC transport.
long double localTraceMeasure(bool horizontal,long double lower,long double upper,bool rz)
{
    require(std::isfinite(lower)&&std::isfinite(upper)&&upper>=lower,"invalid local accepted-trace interval");
    if(rz&&horizontal){require(lower>=0,"negative radial accepted-trace interval");return (upper-lower)*(upper+lower)/2;}
    return upper-lower;
}
bool localSubnormal(long double value)
{return std::fpclassify(value)==FP_SUBNORMAL;}
bool localTraceQuantizedPhase(const GeometricVofIntegratedFlux& parent,bool liquid)
{
    return liquid?(localSubnormal(parent.positiveLiquid)||localSubnormal(parent.negativeLiquid)||localSubnormal(parent.liquid)):
        (localSubnormal(parent.positiveGas)||localSubnormal(parent.negativeGas)||localSubnormal(parent.gas));
}
bool localTraceNeedsQuantization(const GeometricVofIntegratedFlux& parent)
{return localTraceQuantizedPhase(parent,true)||localTraceQuantizedPhase(parent,false)||localSubnormal(parent.rawLiquid);}
long double localTraceCumulativeBranch(long double parent,long double prefixLower,long double prefixUpper)
{
    require(std::isfinite(parent)&&std::isfinite(prefixLower)&&std::isfinite(prefixUpper)
        &&prefixLower>=0&&prefixLower<=prefixUpper&&prefixUpper<=1,"invalid cumulative exceptional-phase P0 interval");
    // A common endpoint is evaluated by ONE identical primitive. Endpoints
    // carry the exact retained branch, not a repaired final query/atom.
    const auto primitive=[&](long double fraction){
        if(fraction==0)return 0.L;if(fraction==1)return parent;
        const auto value=parent*fraction;
        require(std::isfinite(value)&&(parent>0?value>=0&&value<=parent:value<=0&&value>=parent),
            "nonmonotone cumulative exceptional-phase P0 primitive");return value;
    };
    const auto first=primitive(prefixLower),last=primitive(prefixUpper),value=last-first;
    // A subinterval can quantize to zero. This declared finite P0 measure
    // redistributes the SAME native branch across shared endpoints; no VOF
    // stock or accepted parent is erased, clipped or complemented. In an
    // exceptional phase BOTH gross branches use this primitive, even when
    // a branch is zero/normal: normal-gross cancellation can have a subnormal
    // net, so separately rounded gross pieces would otherwise erase that net.
    require(std::isfinite(value)&&(parent>0?last>=first&&value>=0&&value<=parent:
        last<=first&&value<=0&&value>=parent),"invalid cumulative exceptional-phase P0 branch difference");
    return value;
}
GeometricVofIntegratedFlux localTraceQuantizedPart(const GeometricVofIntegratedFlux& parent,long double ratio,
    long double prefixLower,long double prefixUpper,long double rhoL,long double rhoG,bool native)
{
    require(std::isfinite(ratio)&&ratio>=0&&ratio<=1,"invalid local accepted-trace restriction ratio");
    if(ratio==1)return parent;
    const bool liquid=localTraceQuantizedPhase(parent,true),gas=localTraceQuantizedPhase(parent,false);
    if(liquid)require(parent.liquid==parent.positiveLiquid-parent.negativeLiquid,
        "subnormal native P0 liquid net differs from retained gross authority");
    if(gas)require(parent.gas==parent.positiveGas-parent.negativeGas,
        "subnormal native P0 gas net differs from retained gross authority");
    GeometricVofIntegratedFlux result;result.limiter=parent.limiter;
    const std::array<long double GeometricVofIntegratedFlux::*,10> fields{{
        &GeometricVofIntegratedFlux::total,&GeometricVofIntegratedFlux::positiveTotal,&GeometricVofIntegratedFlux::negativeTotal,
        &GeometricVofIntegratedFlux::rawLiquid,&GeometricVofIntegratedFlux::liquid,&GeometricVofIntegratedFlux::gas,
        &GeometricVofIntegratedFlux::positiveLiquid,&GeometricVofIntegratedFlux::negativeLiquid,
        &GeometricVofIntegratedFlux::positiveGas,&GeometricVofIntegratedFlux::negativeGas}};
    for(std::size_t k=0;k<fields.size();++k){
        if((k==4&&liquid)||(k==5&&gas))continue;const auto field=fields[k];
        const bool phaseGross=(liquid&&(k==6||k==7))||(gas&&(k==8||k==9));
        if(phaseGross||(k==3&&localSubnormal(parent.*field)))
            result.*field=localTraceCumulativeBranch(parent.*field,prefixLower,prefixUpper);
        else{result.*field=(parent.*field)*ratio;
            // Original physical-total representability remains mandatory for
            // native atoms. The normal Coons path had no nonzero gate.
            require(std::isfinite(result.*field)&&(!native||ratio==0||parent.*field==0||result.*field!=0),
                "unrepresentable local accepted-trace restriction");}
    }
    // In the exceptional phase its positive/negative gross currents are the
    // authority. Independently rounding its signed net would disagree on a
    // piece (e.g. parent pos=3 denorm, neg=2 denorm).
    if(liquid)result.liquid=result.positiveLiquid-result.negativeLiquid;
    if(gas)result.gas=result.positiveGas-result.negativeGas;
    result.mass=rhoL*result.liquid+rhoG*result.gas;
    require(std::isfinite(result.mass),"unrepresentable local accepted-trace mass");return result;
}
GeometricVofIntegratedFlux localTraceScale(const GeometricVofIntegratedFlux& parent,long double ratio,
    long double rhoL,long double rhoG)
{
    require(std::isfinite(ratio)&&ratio>=0&&ratio<=1,"invalid local accepted-trace restriction ratio");
    if(ratio==1)return parent;
    GeometricVofIntegratedFlux result;result.limiter=parent.limiter;
    const std::array<long double GeometricVofIntegratedFlux::*,10> fields{{
        &GeometricVofIntegratedFlux::total,&GeometricVofIntegratedFlux::positiveTotal,&GeometricVofIntegratedFlux::negativeTotal,
        &GeometricVofIntegratedFlux::rawLiquid,&GeometricVofIntegratedFlux::liquid,&GeometricVofIntegratedFlux::gas,
        &GeometricVofIntegratedFlux::positiveLiquid,&GeometricVofIntegratedFlux::negativeLiquid,
        &GeometricVofIntegratedFlux::positiveGas,&GeometricVofIntegratedFlux::negativeGas}};
    for(auto field:fields){result.*field=(parent.*field)*ratio;require(std::isfinite(result.*field)&&
        (ratio==0||parent.*field==0||result.*field!=0),"unrepresentable local accepted-trace restriction");}
    result.mass=rhoL*result.liquid+rhoG*result.gas;require(std::isfinite(result.mass),"unrepresentable local accepted-trace mass");return result;
}
template<bool FullOnly> using LocalNativeTraceResult=std::conditional_t<FullOnly,GeometricVofNativeFullFlux,GeometricVofTransportFace>;
template<bool FullOnly> using LocalInteriorTraceResult=std::conditional_t<FullOnly,GeometricVofInteriorFullFlux,GeometricVofInteriorFlux>;
template<bool FullOnly>
std::vector<LocalNativeTraceResult<FullOnly>> localNativeTrace(const QuadtreeMesh& mesh,const GeometricVofPublication& old,
    const std::vector<GeometricVofTransportFace>& retained,const std::vector<GeometricVofNativeSubface>& queries,
    long double rhoL,long double rhoG,GeometricVofTransportReport& report)
{
    std::vector<LocalNativeTraceResult<FullOnly>> result;result.reserve(queries.size());const bool rz=old.domain().axisymmetric();
    for(const auto& q:queries){
        require(q.nativeRow>=0&&static_cast<std::size_t>(q.nativeRow)<retained.size(),"local restriction row outside source");
        require(std::isfinite(q.begin.x)&&std::isfinite(q.begin.y)&&std::isfinite(q.end.x)&&std::isfinite(q.end.y),"nonfinite local native subface");
        const auto native=old.domain().face(mesh,q.nativeRow);
        require(native.negativeCell>=0||native.positiveCell>=0,"local restriction requires an active fluid face");
        const auto endpoints=localFaceEndpoints(mesh,native);
        const auto tangent=[&](Point p){return native.horizontal?p.x:p.y;};
        const auto normal=[&](Point p){return native.horizontal?p.y:p.x;};
        const auto lo=std::min(tangent(endpoints[0]),tangent(endpoints[1])),hi=std::max(tangent(endpoints[0]),tangent(endpoints[1]));
        require(normal(q.begin)==normal(endpoints[0])&&normal(q.end)==normal(endpoints[0])&&
            tangent(q.begin)>=lo&&tangent(q.begin)<=hi&&tangent(q.end)>=lo&&tangent(q.end)<=hi,"local subface outside actual native face");
        const auto a=std::min(tangent(q.begin),tangent(q.end)),b=std::max(tangent(q.begin),tangent(q.end));
        if(a==lo&&b==hi){
            if constexpr(FullOnly){const auto& parent=retained[q.nativeRow];
                result.push_back({parent.nativeRow,parent.negativeCell,parent.positiveCell,parent.full});}
            else result.push_back(retained[q.nativeRow]);continue;}
        const auto whole=localTraceMeasure(native.horizontal,lo,hi,rz),part=localTraceMeasure(native.horizontal,a,b,rz);
        require(whole>0,"degenerate local native source face");const auto ratio=part/whole;
        LocalNativeTraceResult<FullOnly> face;face.nativeRow=q.nativeRow;face.negativeCell=native.negativeCell;face.positiveCell=native.positiveCell;
        const auto restrict=[&](const GeometricVofIntegratedFlux& parent){
            if(!localTraceNeedsQuantization(parent))return localTraceScale(parent,ratio,rhoL,rhoG);
            const auto prefixLower=a==lo?0.L:localTraceMeasure(native.horizontal,lo,a,rz)/whole;
            const auto prefixUpper=b==hi?1.L:localTraceMeasure(native.horizontal,lo,b,rz)/whole;
            return localTraceQuantizedPart(parent,ratio,prefixLower,prefixUpper,rhoL,rhoG,true);
        };
        if constexpr(!FullOnly)face.half=restrict(retained[q.nativeRow].half);
        face.full=restrict(retained[q.nativeRow].full);result.push_back(face);++report.localNativeP0Restrictions;
    }return result;
}
template<bool FullOnly>
std::vector<LocalInteriorTraceResult<FullOnly>> localInteriorTrace(const QuadtreeMesh& mesh,const GeometricVofPublication& old,
    const GeometricVofAdvectorSlab& slab,const std::vector<GeometricVofTransportFace>& retained,
    const std::vector<GeometricVofInteriorSubface>& queries,long double rhoL,long double rhoG,long double tolerance,
    GeometricVofTransportReport& report)
{
    struct SidePiece {int row;long double lower,upper,measure,area;};
    const auto cells=old.plicCells().size();const bool rz=old.domain().axisymmetric();
    if(cells>(std::numeric_limits<std::size_t>::max()-1)/4)throw std::length_error("local Coons side index extent");
    const auto none=std::numeric_limits<std::size_t>::max();
    std::vector<SidePiece> pieces(static_cast<std::size_t>(old.domain().nativeFaceCount()));
    std::vector<std::array<std::size_t,2>> routes(pieces.size(),{{none,none}});
    std::vector<std::size_t> sideOffsets(4*cells+1);
    // Each native payload is stored once. Cell-side incidence retains the
    // native-row order without four allocating side vectors per source cell.
    for(int row=0;row<old.domain().nativeFaceCount();++row){const auto f=old.domain().face(mesh,row);const int axis=f.horizontal?1:0;
        const auto endpoints=localFaceEndpoints(mesh,f);
        const auto a=std::min(f.horizontal?endpoints[0].x:endpoints[0].y,f.horizontal?endpoints[1].x:endpoints[1].y);
        const auto b=std::max(f.horizontal?endpoints[0].x:endpoints[0].y,f.horizontal?endpoints[1].x:endpoints[1].y);
        const auto& piece=pieces[row]=SidePiece{row,a,b,localTraceMeasure(f.horizontal,a,b,rz),f.area};
        require(piece.measure>0,"zero local Coons native tangent measure");
        int incidence=0;
        for(int cell:{f.negativeCell,f.positiveCell}){
            if(cell>=0){require(static_cast<std::size_t>(cell)<cells,"local native incidence outside source");
                const auto& box=old.plicCells()[cell].rectangle;const auto c=axis?f.center.y:f.center.x;
                const auto lower=axis?box.lower.y:box.lower.x,upper=axis?box.upper.y:box.upper.x;
                require(c==lower||c==upper,"local native incidence is not an actual cell side");
                const auto key=4*static_cast<std::size_t>(cell)+2*axis+(c==upper?1:0);
                routes[row][incidence]=key;++sideOffsets[key+1];}
            ++incidence;
        }
    }
    for(std::size_t key=0;key<4*cells;++key)sideOffsets[key+1]+=sideOffsets[key];
    std::vector<std::size_t> sideRows(sideOffsets.back());
    {auto cursor=sideOffsets;
        for(std::size_t row=0;row<routes.size();++row)for(auto key:routes[row])if(key!=none)sideRows[cursor[key]++]=row;}
    std::vector<std::array<std::size_t,2>>().swap(routes);
    const std::array<long double GeometricVofIntegratedFlux::*,10> fields{{
        &GeometricVofIntegratedFlux::total,&GeometricVofIntegratedFlux::positiveTotal,&GeometricVofIntegratedFlux::negativeTotal,
        &GeometricVofIntegratedFlux::rawLiquid,&GeometricVofIntegratedFlux::liquid,&GeometricVofIntegratedFlux::gas,
        &GeometricVofIntegratedFlux::positiveLiquid,&GeometricVofIntegratedFlux::negativeLiquid,
        &GeometricVofIntegratedFlux::positiveGas,&GeometricVofIntegratedFlux::negativeGas}};
    struct SideIntegral {std::array<std::array<Sum,10>,2> flux;std::array<Sum,2> normal;long double limiter=1;};
    const auto sideIntegral=[&](int cell,int axis,int side,long double lo,long double hi){
        SideIntegral value;Sum coverage;const bool horizontal=axis==1;
        const auto queryMeasure=localTraceMeasure(horizontal,lo,hi,rz);
        const auto key=4*static_cast<std::size_t>(cell)+2*axis+side;
        for(auto at=sideOffsets[key];at<sideOffsets[key+1];++at){const auto& piece=pieces[sideRows[at]];
            const auto first=std::max(lo,piece.lower),last=std::min(hi,piece.upper);if(last<=first)continue;
            const auto part=localTraceMeasure(horizontal,first,last,rz),ratio=part/piece.measure;coverage.add(part);
            // The paired result's limiter includes half before full even
            // when only the full phase integrals are requested.
            if constexpr(FullOnly)value.limiter=std::min(value.limiter,retained[piece.row].half.limiter);
            for(int extent=FullOnly?1:0;extent<2;++extent){const auto& flux=extent?retained[piece.row].full:retained[piece.row].half;
                if(localTraceNeedsQuantization(flux)){
                    const auto prefixLower=first==piece.lower?0.L:localTraceMeasure(horizontal,piece.lower,first,rz)/piece.measure;
                    const auto prefixUpper=last==piece.upper?1.L:localTraceMeasure(horizontal,piece.lower,last,rz)/piece.measure;
                    const auto restricted=localTraceQuantizedPart(flux,ratio,prefixLower,prefixUpper,rhoL,rhoG,false);
                    for(std::size_t k=0;k<fields.size();++k)value.flux[extent][k].add(restricted.*fields[k]);
                }else{
                    for(std::size_t k=0;k<fields.size();++k)value.flux[extent][k].add(ratio*(flux.*fields[k]));
                }
                value.limiter=std::min(value.limiter,flux.limiter);}
            if constexpr(!FullOnly){
                value.normal[0].add(ratio*piece.area*slab.beginNormal()[piece.row]);
                value.normal[1].add(ratio*piece.area*slab.endNormal()[piece.row]);
            }
        }
        require(std::abs(coverage.value-queryMeasure)<=tolerance*queryMeasure,"local Coons side subfaces do not cover the query");
        return value;
    };
    using WholeSides=std::array<std::array<std::optional<SideIntegral>,2>,2>;
    const auto evaluate=[&](const GeometricVofInteriorSubface& q,WholeSides* wholeCache,
        std::vector<WholeSides>* domainCache,std::vector<long double>& cuts){
        const auto wholeSideIntegral=[&](int cell,int axis,int side,long double lo,long double hi)->const SideIntegral&{
            // Serial replay allocates at the same original first side use.
            if(!wholeCache){if(domainCache->empty())domainCache->resize(cells);wholeCache=&(*domainCache)[cell];}
            auto& retainedWhole=(*wholeCache)[axis][side];
            if(!retainedWhole)retainedWhole=sideIntegral(cell,axis,side,lo,hi);
            return *retainedWhole;
        };
        require(q.nativeCell>=0&&static_cast<std::size_t>(q.nativeCell)<cells&&old.cells()[q.nativeCell].nativeCell==q.nativeCell,"local interior native cell outside source");
        require(std::isfinite(q.begin.x)&&std::isfinite(q.begin.y)&&std::isfinite(q.end.x)&&std::isfinite(q.end.y),"nonfinite local interior endpoints");
        const bool horizontal=q.begin.y==q.end.y,vertical=q.begin.x==q.end.x;require(horizontal!=vertical,"local interior face must be nondegenerate and axis aligned");
        const int axis=horizontal?1:0;const auto& box=old.plicCells()[q.nativeCell].rectangle;
        const long double lower=axis?box.lower.y:box.lower.x,upper=axis?box.upper.y:box.upper.x;
        const auto coordinate=axis?q.begin.y:q.begin.x;
        const auto ta=horizontal?q.begin.x:q.begin.y,tb=horizontal?q.end.x:q.end.y;
        const auto lo=std::min(ta,tb),hi=std::max(ta,tb);
        const long double tlo=horizontal?box.lower.x:box.lower.y,thi=horizontal?box.upper.x:box.upper.y;
        require(coordinate>lower&&coordinate<upper&&lo>=tlo&&hi<=thi,"local interior face outside actual owner cell");
        const auto tangentMeasure=localTraceMeasure(horizontal,tlo,thi,rz);
        const auto queryMeasure=localTraceMeasure(horizontal,lo,hi,rz);require(queryMeasure>0,"zero local interior tangent measure");
        // Both coordinates are normalized physical VOLUME coordinates. For
        // RZ the radial coordinate uses stable (r-r0)*(r+r0), never r*r-r0*r0.
        const auto normalMeasure=localTraceMeasure(!horizontal,lower,upper,rz);
        const auto prefixMeasure=localTraceMeasure(!horizontal,lower,coordinate,rz);
        const auto xi=prefixMeasure/normalMeasure;require(std::isfinite(xi)&&xi>0&&xi<1,"invalid local Coons normal coordinate");
        std::array<const SideIntegral*,2> whole;
        std::array<SideIntegral,2> prefix;
        for(int side=0;side<2;++side){whole[side]=&wholeSideIntegral(q.nativeCell,1-axis,side,lower,upper);
            prefix[side]=sideIntegral(q.nativeCell,1-axis,side,lower,coordinate);}
        std::array<std::array<long double,10>,2> adjustment{};std::array<long double,2> normalAdjustment{};
        const auto departure=[&](long double p0,long double p1,long double w0,long double w1){
            Sum value;value.add(p1);value.add(-p0);value.add(-xi*w1);value.add(xi*w0);
            require(std::isfinite(value.value),"nonfinite local Coons prefix correction");return value.value;};
        long double limiter=std::min(whole[0]->limiter,whole[1]->limiter);
        for(int extent=FullOnly?1:0;extent<2;++extent){
            for(std::size_t k=0;k<fields.size();++k)adjustment[extent][k]=departure(prefix[0].flux[extent][k].value,
                prefix[1].flux[extent][k].value,whole[0]->flux[extent][k].value,whole[1]->flux[extent][k].value);
            if constexpr(!FullOnly)normalAdjustment[extent]=departure(prefix[0].normal[extent].value,prefix[1].normal[extent].value,
                whole[0]->normal[extent].value,whole[1]->normal[extent].value);
        }
        // Within each original side breakpoint the Coons normal current is
        // constant in the tangent volume coordinate. Integrate its sign there,
        // preserving BOTH accepted positive/negative phase source currents.
        cuts.assign({lo,hi});
        for(int side=0;side<2;++side){const auto key=4*static_cast<std::size_t>(q.nativeCell)+2*axis+side;
            for(auto at=sideOffsets[key];at<sideOffsets[key+1];++at){const auto& piece=pieces[sideRows[at]];
                for(auto endpoint:{piece.lower,piece.upper})if(endpoint>lo&&endpoint<hi)cuts.push_back(endpoint);}}
        std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
        std::array<std::array<Sum,10>,2> accum;std::array<Sum,2> normals;
        for(std::size_t part=1;part<cuts.size();++part){
            const auto first=cuts[part-1],last=cuts[part];
            const auto low=sideIntegral(q.nativeCell,axis,0,first,last),high=sideIntegral(q.nativeCell,axis,1,first,last);
            const auto fraction=localTraceMeasure(horizontal,first,last,rz)/tangentMeasure;
            limiter=std::min({limiter,low.limiter,high.limiter});
            for(int extent=FullOnly?1:0;extent<2;++extent){std::array<long double,10> current{};
                for(std::size_t k=0;k<fields.size();++k){Sum value;value.add((1-xi)*low.flux[extent][k].value);
                    value.add(xi*high.flux[extent][k].value);value.add(-fraction*adjustment[extent][k]);current[k]=value.value;
                    require(std::isfinite(current[k]),"nonfinite local Coons phase current");}
                // Authoritative signed nets are extended directly, rather than
                // subtracting rounded exported gross values after counterflow.
                for(std::size_t k:{0U,3U,4U,5U})accum[extent][k].add(current[k]);
                for(std::size_t positive:{6U,8U}){const auto negative=positive+1;
                    if(current[positive]>=0)accum[extent][positive].add(current[positive]);
                    else accum[extent][negative].add(-current[positive]);
                    if(current[negative]>=0)accum[extent][negative].add(current[negative]);
                    else accum[extent][positive].add(-current[negative]);}
                if constexpr(!FullOnly){Sum value;value.add((1-xi)*low.normal[extent].value);value.add(xi*high.normal[extent].value);
                    value.add(-fraction*normalAdjustment[extent]);normals[extent].add(value.value);}
            }
        }
        LocalInteriorTraceResult<FullOnly> face;face.nativeCell=q.nativeCell;
        for(int extent=FullOnly?1:0;extent<2;++extent){auto& flux=[&]() -> GeometricVofIntegratedFlux& {
                if constexpr(FullOnly)return face.full;else return extent?face.full:face.half;}();
            for(std::size_t k=0;k<fields.size();++k){flux.*fields[k]=accum[extent][k].value;require(std::isfinite(flux.*fields[k]),"nonfinite local Coons material flux");}
            // The finite subnormal phase allocation has gross authority on
            // every atom, including this Coons extension. Keep the original
            // normal signed-current arithmetic outside that exceptional case.
            if(localTraceQuantizedPhase(flux,true))flux.liquid=flux.positiveLiquid-flux.negativeLiquid;
            if(localTraceQuantizedPhase(flux,false))flux.gas=flux.positiveGas-flux.negativeGas;
            flux.positiveTotal=flux.positiveLiquid+flux.positiveGas;flux.negativeTotal=flux.negativeLiquid+flux.negativeGas;
            flux.mass=rhoL*flux.liquid+rhoG*flux.gas;flux.limiter=limiter;
            require(std::isfinite(flux.mass)&&std::isfinite(flux.positiveTotal)&&std::isfinite(flux.negativeTotal),"nonfinite local Coons mass or gross");}
        const auto area=(hi-lo)*(rz?CartesianMeasure::twoPi*(horizontal?(lo+hi)/2:coordinate):1.L);require(std::isfinite(area)&&area>0,"invalid local Coons physical atom area");
        if constexpr(!FullOnly){
            const auto u0=normals[0].value/area,u1=normals[1].value/area;require(std::isfinite(u0)&&std::isfinite(u1),"nonfinite local Coons normal sample");
            if(axis){face.beginVelocity.y=u0;face.endVelocity.y=u1;}else{face.beginVelocity.x=u0;face.endVelocity.x=u1;}
        }
        return face;
    };
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
    // Stable owner-cell routing keeps each cell's original query/side order.
    // A worker retains only four whole-side integrals and its current cuts;
    // no full-domain table of large optional SideIntegral objects is needed.
    const int entryErrno=errno;LocalFaceFloatingState entry;
    const bool eligible=queries.size()>=128&&cells<=static_cast<std::size_t>(std::numeric_limits<std::ptrdiff_t>::max())
        &&localFaceQualified(entry);
    if(entry.valid)entry.restore();errno=entryErrno;
    if(eligible){
        try{
            LocalFaceCleanup cleanup;
            const auto workers=std::min<std::size_t>(cells,static_cast<std::size_t>(omp_get_max_threads()));
            LocalFaceBudget budget;budget.add(cells+1,2*sizeof(std::size_t));budget.add(queries.size(),sizeof(std::size_t));
            budget.add(workers,sizeof(WholeSides)+sizeof(std::vector<long double>)+2*sizeof(long double));
            budget.add(workers,sizeof(LocalFaceWorkerFrame)+16384u);
            budget.add(sideRows.size(),4*sizeof(long double));
            std::vector<std::size_t> offsets(cells+1),order(queries.size());
            for(const auto& q:queries){
                require(q.nativeCell>=0&&static_cast<std::size_t>(q.nativeCell)<cells,"local interior native cell outside source");
                ++offsets[static_cast<std::size_t>(q.nativeCell)+1];
            }
            for(std::size_t cell=0;cell<cells;++cell)offsets[cell+1]+=offsets[cell];
            auto cursor=offsets;
            for(std::size_t i=0;i<queries.size();++i)order[cursor[queries[i].nativeCell]++]=i;
            std::vector<LocalInteriorTraceResult<FullOnly>> result(queries.size());unsigned actualWorkers=0;
            const int setupErrno=errno;LocalFaceFloatingState setup;
            const bool setupCompatible=setup.capture()&&setupErrno==entryErrno&&localFaceSameControl(setup,entry);
            entry.restore();errno=entryErrno;
            const int team=static_cast<int>(workers);
            if(setupCompatible&&team>1&&localFaceRows(cells,team,[&](std::size_t cell){
                WholeSides wholeCache;std::vector<long double> cuts;
                for(auto at=offsets[cell];at<offsets[cell+1];++at){const auto i=order[at];result[i]=evaluate(queries[i],&wholeCache,nullptr,cuts);}
            },actualWorkers)){
                report.localInteriorP0Extensions+=queries.size();cleanup.arm();return result;
            }
        }catch(...){/* Replay the original query order and first failure. */}
        entry.restore();errno=entryErrno;
    }
#endif
    std::vector<LocalInteriorTraceResult<FullOnly>> result;result.reserve(queries.size());
    std::vector<WholeSides> wholeCache;std::vector<long double> cuts;
    for(const auto& q:queries){
        result.push_back(evaluate(q,nullptr,&wholeCache,cuts));++report.localInteriorP0Extensions;
    }return result;
}
}
std::vector<GeometricVofNativeFullFlux> GeometricVofTransport::restrictNativeFacesFull(const QuadtreeMesh& mesh,
    const std::vector<GeometricVofNativeSubface>& queries,GeometricVofTransportReport* diagnostics) const
{
    require(bool(data_),"consumed transport");validateFor(mesh,data_->original,data_->slab,data_->rhoL,data_->rhoG);
    require(data_->original.geometryKind()==GeometricVofGeometryKind::LocalPlic,"full-only restriction requires LocalPlic");
    GeometricVofTransportReport report;
    auto result=localNativeTrace<true>(mesh,data_->original,data_->faces,queries,data_->rhoL,data_->rhoG,report);
    if(diagnostics)*diagnostics=report;return result;
}
std::vector<GeometricVofInteriorFullFlux> GeometricVofTransport::measureInteriorFacesFull(const QuadtreeMesh& mesh,
    const std::vector<GeometricVofInteriorSubface>& queries,GeometricVofTransportReport* diagnostics) const
{
    require(bool(data_),"consumed transport");validateFor(mesh,data_->original,data_->slab,data_->rhoL,data_->rhoG);
    require(data_->original.geometryKind()==GeometricVofGeometryKind::LocalPlic,"full-only extension requires LocalPlic");
    GeometricVofTransportReport report;
    if(queries.empty()){if(diagnostics)*diagnostics=report;return {};}
    auto result=localInteriorTrace<true>(mesh,data_->original,data_->slab,data_->faces,queries,data_->rhoL,data_->rhoG,data_->config.relativeConservationTolerance,report);
    if(diagnostics)*diagnostics=report;return result;
}
std::vector<GeometricVofTransportFace> GeometricVofTransport::restrictNativeFaces(const QuadtreeMesh& mesh,
    const std::vector<GeometricVofNativeSubface>& queries,GeometricVofTransportReport* diagnostics) const
{
    require(bool(data_),"consumed transport");validateFor(mesh,data_->original,data_->slab,data_->rhoL,data_->rhoG);
    if(data_->original.geometryKind()==GeometricVofGeometryKind::LocalPlic){GeometricVofTransportReport report;
        auto result=localNativeTrace<false>(mesh,data_->original,data_->faces,queries,data_->rhoL,data_->rhoG,report);
        if(diagnostics)*diagnostics=report;return result;}
    const auto& old=data_->original;const auto& slab=data_->slab;const auto end=data_->endTime;
    const bool rz=old.domain().axisymmetric();GeometricVofTransportReport report;
    if(queries.empty()){if(diagnostics)*diagnostics=report;return {};}
    SourceGeometry geometry(mesh,old);std::vector<GeometricVofTransportFace> result;result.reserve(queries.size());
    const auto same=[](Point a,Point b){return a.x==b.x&&a.y==b.y;};
    const auto finite=[](Point p){return std::isfinite(p.x)&&std::isfinite(p.y);};
    for(const auto& q:queries){
        require(q.nativeRow>=0&&static_cast<std::size_t>(q.nativeRow)<data_->faces.size(),"native restriction row outside source");
        const auto& actual=old.faces()[q.nativeRow];const auto native=old.domain().face(mesh,q.nativeRow);
        const auto a=nativePoint(old,actual.vertices[0]),b=nativePoint(old,actual.vertices[1]);
        require(finite(q.begin)&&finite(q.end),"nonfinite native subface endpoint");
        const auto origin=old.frameOrigin();
        const auto local=[&](Point world){
            const auto value=Cut::anchoredDifference({world,{}},{origin,{}});
            require(Cut::sameAnchoredPoint({value,origin},{world,{}}),
                "native subface world-to-frame conversion lost an exact coordinate tail");return value;};
        const auto first=local(q.begin),last=local(q.end);
        const auto inside=[&](Point p){return native.horizontal?
            p.y==a.y&&p.x>=std::min(a.x,b.x)&&p.x<=std::max(a.x,b.x):
            p.x==a.x&&p.y>=std::min(a.y,b.y)&&p.y<=std::max(a.y,b.y);};
        require(inside(first)&&inside(last),"subface is not an exact interval of its native source");
        if((same(first,a)&&same(last,b))||(same(first,b)&&same(last,a))){result.push_back(data_->faces[q.nativeRow]);continue;}
        GeometricVofTransportFace face;face.nativeRow=q.nativeRow;face.negativeCell=native.negativeCell;face.positiveCell=native.positiveCell;
        if(same(first,last)){result.push_back(face);continue;}
        const auto delta=Cut::anchoredDifference({q.end,{}},{q.begin,{}});
        const auto length=native.horizontal?std::abs(delta.x):std::abs(delta.y);
        const auto area=length*(rz?CartesianMeasure::twoPi*(native.horizontal?(q.begin.x/2+q.end.x/2):q.begin.x):1.L);
        require(length>0&&std::isfinite(area)&&area>=0&&(area>0||rz&&!native.horizontal&&q.begin.x==0),
            "unrepresentable physical native subface measure");
        const auto v0=slab.beginVelocity()[q.nativeRow],v1=slab.endVelocity()[q.nativeRow];
        const auto total=FaceTime::integrateNormal(slab.beginNormal()[q.nativeRow],slab.endNormal()[q.nativeRow],end.begin,end.end,area);
        faceTimeDiagnostics(total,report);Gross liquid,gas;
        if(total.fullPositive+total.fullNegative>0){
            FaceTime::validateDepartureSupport(first,last,v0,v1,end.begin,end.end,geometry.originalDomainBounds());
            const auto images=geometry.query(FaceTime::departureBounds(first,last,v0,v1,end.begin,end.end),report);
            bool pureGas=true,pureLiquid=true;for(auto image:images){pureGas&=old.alpha()[image.cell]==0;pureLiquid&=old.alpha()[image.cell]==1;}
            if(pureGas)gas=total;else if(pureLiquid)liquid=total;
            else for(auto image:images)for(bool phase:{true,false})for(const auto& p:geometry.phase(image,phase)){
                const auto measure=FaceTime::integrateAnchoredPolygon(p,first,last,v0,v1,end.begin,end.end,
                    old.frameOrigin(),rz,data_->config.maximumTemporalEvents,{image.sx,image.sy,{image.tx,image.ty}});
                ++report.polygonIntegrals;faceTimeDiagnostics(measure,report);auto& aggregate=phase?liquid:gas;
                aggregate.halfPositive+=measure.halfPositive;aggregate.halfNegative+=measure.halfNegative;
                aggregate.fullPositive+=measure.fullPositive;aggregate.fullNegative+=measure.fullNegative;
            }
        }
        const auto raw=[&](bool half){GeometricVofIntegratedFlux f;
            f.positiveTotal=half?total.halfPositive:total.fullPositive;f.negativeTotal=half?total.halfNegative:total.fullNegative;
            f.positiveLiquid=half?liquid.halfPositive:liquid.fullPositive;f.negativeLiquid=half?liquid.halfNegative:liquid.fullNegative;
            f.positiveGas=half?gas.halfPositive:gas.fullPositive;f.negativeGas=half?gas.halfNegative:gas.fullNegative;
            const auto coverage=[&](long double t,long double l,long double g){
                require(t==0?l==0&&g==0:l>=0&&g>=0&&std::isfinite(l+g)&&l+g>0,"invalid subface raw geometry");
                if(t>0){const auto error=std::abs(l+g-t)/t;report.maximumCoverageResidual=std::max(report.maximumCoverageResidual,error);
                    require(error<=data_->config.relativeConservationTolerance,"subface actual fan coverage exceeds original budget");}};
            coverage(f.positiveTotal,f.positiveLiquid,f.positiveGas);coverage(f.negativeTotal,f.negativeLiquid,f.negativeGas);return f;};
        const auto measure=[&](bool half){const auto r=raw(half);auto f=(half?data_->halfRestriction:data_->fullRestriction)->measure(q.nativeRow,r);
            f.mass=data_->rhoL*f.liquid+data_->rhoG*f.gas;require(std::isfinite(f.mass),"unrepresentable restricted material integral");
            for(const auto& t:{std::pair<long double,long double>{f.positiveTotal,r.positiveTotal},{f.negativeTotal,r.negativeTotal}})
                if(t.second>0)report.maximumCoverageResidual=std::max(report.maximumCoverageResidual,std::abs(t.first-t.second)/t.second);
            return f;};
        face.half=measure(true);face.full=measure(false);result.push_back(face);
    }
    if(diagnostics)*diagnostics=report;return result;
}
namespace
{
struct InteriorVelocityFit
{
    Point center{};
    long double scale=0;
    Eigen::Matrix<double,3,4> coefficients; // begin x/y, end x/y.
};
InteriorVelocityFit prepareInteriorVelocity(const QuadtreeMesh& mesh,
    const GeometricVofPublication& source,const GeometricVofAdvectorSlab& slab,
    int nativeCell,const std::vector<int>& incident)
{
    require(incident.size()>=3,"insufficient fixed-cell interior velocity support");
    const auto& cell=mesh.cells(mesh.cellIndexToLabel[nativeCell]);
    const auto lo=mesh.verts(cell.cellToVert[0]).coord(),hi=mesh.verts(cell.cellToVert[2]).coord();
    const auto center=cell.centerCoord();InteriorVelocityFit result;result.center={center.x,center.y};
    result.scale=std::min(static_cast<long double>(hi.x)-lo.x,static_cast<long double>(hi.y)-lo.y);
    require(std::isfinite(result.scale)&&result.scale>0,"invalid fixed-cell interior velocity scale");
    Eigen::MatrixXd design(incident.size(),3),values(incident.size(),4);Eigen::Index index=0;
    for(int row:incident){const auto point=samplingPoint(source,source.domain().face(mesh,row));
        const auto x=(point.x-result.center.x)/result.scale,y=(point.y-result.center.y)/result.scale;
        const auto weight=1/std::sqrt(1+x*x+y*y);
        require(std::isfinite(x)&&std::isfinite(y)&&std::isfinite(weight)&&weight>0,"unrepresentable fixed-cell interior velocity support");
        design(index,0)=static_cast<double>(weight);design(index,1)=static_cast<double>(weight*x);design(index,2)=static_cast<double>(weight*y);
        const auto a=slab.beginVelocity()[row],b=slab.endVelocity()[row];
        values(index,0)=static_cast<double>(weight)*static_cast<double>(a.x);values(index,1)=static_cast<double>(weight)*static_cast<double>(a.y);
        values(index,2)=static_cast<double>(weight)*static_cast<double>(b.x);values(index,3)=static_cast<double>(weight)*static_cast<double>(b.y);++index;
    }
    require(design.allFinite()&&values.allFinite(),"nonfinite fixed-cell interior velocity observations");
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(design);qr.setThreshold(64*std::numeric_limits<double>::epsilon());
    require(qr.rank()==3,"rank-deficient fixed-cell interior velocity support");
    result.coefficients=qr.solve(values);require(result.coefficients.allFinite(),"nonfinite fixed-cell interior velocity fit");return result;
}
std::array<Point,2> evaluateInteriorVelocity(const InteriorVelocityFit& fit,Point point)
{
    const auto x=(point.x-fit.center.x)/fit.scale,y=(point.y-fit.center.y)/fit.scale;
    const auto component=[&](int column){return static_cast<long double>(fit.coefficients(0,column))+
        static_cast<long double>(fit.coefficients(1,column))*x+static_cast<long double>(fit.coefficients(2,column))*y;};
    const std::array<Point,2> result{{{component(0),component(1)},{component(2),component(3)}}};
    for(auto p:result)require(std::isfinite(p.x)&&std::isfinite(p.y),"unrepresentable interior endpoint velocity");return result;
}
// Two terms of a phase INTEGRAL, in volume units. The existing strict Core
// primitive provides exact (head+tail)-(head+tail) with one final RN. These
// arithmetic carriers are never clipped, measured or published as geometry.
struct InteriorPhaseIntegral {long double head=0,tail=0;};
long double accuratePhaseDifference(InteriorPhaseIntegral a,InteriorPhaseIntegral b)
{return Cut::anchoredDifference({{a.head,0},{a.tail,0}},{{b.head,0},{b.tail,0}}).x;}
struct InteriorPhasePair {InteriorPhaseIntegral liquid,gas;};
InteriorPhasePair normalizeInteriorPhase(long double liquid,long double gas,long double total,
    long double tolerance,GeometricVofTransportReport& report)
{
    if(total==0){require(liquid==0&&gas==0,"zero interior normal branch has geometric phase transport");return {};}
    require(total>0&&std::isfinite(total)&&liquid>=0&&gas>=0&&std::isfinite(liquid+gas)&&liquid+gas>0,
        "invalid actual interior two-phase integral");
    const auto residual=std::abs(liquid+gas-total)/total;
    report.maximumCoverageResidual=std::max(report.maximumCoverageResidual,residual);
    require(residual<=tolerance,"actual interior fan coverage exceeds original conservation certificate");
    if(liquid==0)return {{},{total,0}};if(gas==0)return {{total,0},{}};
    // Exactly the original native minor-normalization scalar arithmetic.
    // Preserve the dominant complement as two terms instead of materializing
    // its rounded copy before a positive-minus-negative micro-net difference.
    const auto minor=total*((liquid<=gas?liquid:gas)/(liquid+gas));
    require(std::isfinite(minor)&&minor>0&&total-minor>0,"unrepresentable positive interior normalized phase");
    return liquid<=gas?InteriorPhasePair{{minor,0},{total,-minor}}:InteriorPhasePair{{total,-minor},{minor,0}};
}
struct InteriorGrossSum
{
    std::array<Sum,4> sums;
    void add(const Gross& measure)
    {
        const std::array<long double,4> values{{measure.halfPositive,measure.halfNegative,measure.fullPositive,measure.fullNegative}};
        for(std::size_t i=0;i<values.size();++i){require(std::isfinite(values[i])&&values[i]>=0,"invalid positive interior geometric contribution");
            sums[i].add(values[i]);require(std::isfinite(sums[i].value),"overflowing positive interior geometric sum");}
    }
    Gross values() const
    {Gross result;result.halfPositive=sums[0].value;result.halfNegative=sums[1].value;result.fullPositive=sums[2].value;result.fullNegative=sums[3].value;return result;}
};
}
std::vector<GeometricVofInteriorFlux> GeometricVofTransport::measureInteriorFaces(const QuadtreeMesh& mesh,
    const std::vector<GeometricVofInteriorSubface>& queries,GeometricVofTransportReport* diagnostics) const
{
    require(bool(data_),"consumed transport");validateFor(mesh,data_->original,data_->slab,data_->rhoL,data_->rhoG);
    if(data_->original.geometryKind()==GeometricVofGeometryKind::LocalPlic){GeometricVofTransportReport report;
        if(queries.empty()){if(diagnostics)*diagnostics=report;return {};}
        auto result=localInteriorTrace<false>(mesh,data_->original,data_->slab,data_->faces,queries,data_->rhoL,data_->rhoG,data_->config.relativeConservationTolerance,report);
        if(diagnostics)*diagnostics=report;return result;}
    GeometricVofTransportReport report;if(queries.empty()){if(diagnostics)*diagnostics=report;return {};}
    const auto& old=data_->original;const auto& slab=data_->slab;const auto end=data_->endTime;const auto origin=old.frameOrigin();
    const bool rz=old.domain().axisymmetric();SourceGeometry geometry(mesh,old);
    std::vector<std::vector<int>> incident(old.cells().size());
    for(const auto& face:old.faces()){const auto native=old.domain().face(mesh,face.nativeRow);
        if(native.negativeCell>=0)incident[native.negativeCell].push_back(face.nativeRow);
        if(native.positiveCell>=0)incident[native.positiveCell].push_back(face.nativeRow);}
    std::map<int,InteriorVelocityFit> fits;std::vector<GeometricVofInteriorFlux> result;result.reserve(queries.size());
    for(const auto& q:queries){
        require(q.nativeCell>=0&&static_cast<std::size_t>(q.nativeCell)<old.cells().size()&&old.cells()[q.nativeCell].nativeCell==q.nativeCell,
            "interior query source native cell outside original owner");
        require(std::isfinite(q.begin.x)&&std::isfinite(q.begin.y)&&std::isfinite(q.end.x)&&std::isfinite(q.end.y),"nonfinite interior face endpoints");
        const bool horizontal=q.begin.y==q.end.y;const bool vertical=q.begin.x==q.end.x;
        require(horizontal!=vertical,"interior face must be nondegenerate and axis aligned");
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[q.nativeCell]);
        const auto lo=mesh.verts(cell.cellToVert[0]).coord(),hi=mesh.verts(cell.cellToVert[2]).coord();
        const auto inside=[&](Point p){return horizontal?p.y>lo.y&&p.y<hi.y&&p.x>=lo.x&&p.x<=hi.x:
            p.x>lo.x&&p.x<hi.x&&p.y>=lo.y&&p.y<=hi.y;};
        require(inside(q.begin)&&inside(q.end),"interior face is not strictly inside its source native cell");
        const auto local=[&](Point world){const auto value=Cut::anchoredDifference({world,{}},{origin,{}});
            require(Cut::sameAnchoredPoint({value,origin},{world,{}}),"interior world-to-frame conversion lost an exact coordinate tail");return value;};
        const auto a=local(q.begin),b=local(q.end),delta=Cut::anchoredDifference({q.end,{}},{q.begin,{}});
        const auto length=horizontal?std::abs(delta.x):std::abs(delta.y);
        const auto radius=horizontal?q.begin.x/2+q.end.x/2:q.begin.x;
        const auto area=length*(rz?CartesianMeasure::twoPi*radius:1.L);
        require(length>0&&std::isfinite(area)&&area>0&&(!rz||(q.begin.x>=0&&q.end.x>=0)),"unrepresentable positive physical interior face measure");
        Point centroid{q.begin.x/2+q.end.x/2,q.begin.y/2+q.end.y/2};
        if(rz&&horizontal){const auto r0=q.begin.x,r1=q.end.x;
            centroid.x=2*(r0*r0+r0*r1+r1*r1)/(3*(r0+r1));}
        require(std::isfinite(centroid.x)&&std::isfinite(centroid.y)&&(!rz||centroid.x>0),"unrepresentable physical interior face centroid");
        auto fit=fits.find(q.nativeCell);if(fit==fits.end())fit=fits.emplace(q.nativeCell,prepareInteriorVelocity(mesh,old,slab,q.nativeCell,incident[q.nativeCell])).first;
        const auto velocity=evaluateInteriorVelocity(fit->second,centroid);const auto v0=velocity[0],v1=velocity[1];
        const auto total=FaceTime::integrateNormal(horizontal?v0.y:v0.x,horizontal?v1.y:v1.x,end.begin,end.end,area);
        faceTimeDiagnostics(total,report);Gross liquid,gas;
        if(total.fullPositive+total.fullNegative>0){
            FaceTime::validateDepartureSupport(a,b,v0,v1,end.begin,end.end,geometry.originalDomainBounds());
            const auto images=geometry.query(FaceTime::departureBounds(a,b,v0,v1,end.begin,end.end),report);
            bool pureGas=true,pureLiquid=true;for(auto image:images){pureGas&=old.alpha()[image.cell]==0;pureLiquid&=old.alpha()[image.cell]==1;}
            if(pureGas)gas=total;else if(pureLiquid)liquid=total;
            else {InteriorGrossSum liquidSum,gasSum;
                for(auto image:images)for(bool phase:{true,false})for(const auto& polygon:geometry.phase(image,phase)){
                    Gross measure;
                    try {measure=FaceTime::integrateAnchoredPolygon(polygon,a,b,v0,v1,end.begin,end.end,origin,rz,
                        data_->config.maximumTemporalEvents,{image.sx,image.sy,{image.tx,image.ty}});}
                    catch(const std::exception& error){std::ostringstream detail;detail.precision(std::numeric_limits<long double>::max_digits10);
                        detail<<"actual interior face-time refused ownerCell="<<q.nativeCell<<" sourceCell="<<image.cell<<" liquid="<<phase
                            <<" clock="<<end.begin<<","<<end.end<<": "<<error.what();require(false,detail.str());}
                    ++report.polygonIntegrals;faceTimeDiagnostics(measure,report);(phase?liquidSum:gasSum).add(measure);
                }
                liquid=liquidSum.values();gas=gasSum.values();
            }
        }
        const auto accepted=[&](bool half){GeometricVofIntegratedFlux flux;
            flux.positiveTotal=half?total.halfPositive:total.fullPositive;flux.negativeTotal=half?total.halfNegative:total.fullNegative;
            const auto positive=normalizeInteriorPhase(half?liquid.halfPositive:liquid.fullPositive,half?gas.halfPositive:gas.fullPositive,
                flux.positiveTotal,data_->config.relativeConservationTolerance,report);
            const auto negative=normalizeInteriorPhase(half?liquid.halfNegative:liquid.fullNegative,half?gas.halfNegative:gas.fullNegative,
                flux.negativeTotal,data_->config.relativeConservationTolerance,report);
            flux.positiveLiquid=accuratePhaseDifference(positive.liquid,{});flux.negativeLiquid=accuratePhaseDifference(negative.liquid,{});
            flux.positiveGas=accuratePhaseDifference(positive.gas,{});flux.negativeGas=accuratePhaseDifference(negative.gas,{});
            flux.total=accuratePhaseDifference({flux.positiveTotal,0},{flux.negativeTotal,0});
            flux.rawLiquid=flux.liquid=accuratePhaseDifference(positive.liquid,negative.liquid);flux.gas=accuratePhaseDifference(positive.gas,negative.gas);
            flux.mass=data_->rhoL*flux.liquid+data_->rhoG*flux.gas;
            require(std::isfinite(flux.mass)&&(!flux.liquid||data_->rhoL*flux.liquid!=0)&&(!flux.gas||data_->rhoG*flux.gas!=0),"unrepresentable interior material integral");
            flux.limiter=1;return flux;};
        GeometricVofInteriorFlux face;face.nativeCell=q.nativeCell;face.beginVelocity=v0;face.endVelocity=v1;
        face.half=accepted(true);face.full=accepted(false);result.push_back(face);
    }
    if(diagnostics)*diagnostics=report;return result;
}
}
