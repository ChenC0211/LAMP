#include "Modules/Hydrodynamics/HydroFrozenMidpoint.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <chrono>
#include <cfenv>
#include <cerrno>
#include <exception>
#include <cstdint>
#include <cstring>
#include <limits>
#include <mutex>
#if defined(_OPENMP)
#include <omp.h>
#endif
#if defined(_MSC_VER)
#include <float.h>
#endif
#if defined(_M_X64) || defined(__SSE2__)
#include <xmmintrin.h>
#endif

namespace LAMP::Hydrodynamics
{
namespace Detail
{
using FrozenAssemblyMatrix=HydroFrozenMidpoint::Matrix;
// Optional allocation/lookup only. Never encloses fresh numerical operations.
class FrozenAssemblyFloatingRestore
{
    std::fenv_t native_{};int error_=errno;
#if defined(_M_X64) || defined(__SSE2__)
    unsigned csr_=_mm_getcsr();
#endif
public:
    bool captured=std::fegetenv(&native_)==0;
    ~FrozenAssemblyFloatingRestore() noexcept
    {
        if(captured&&std::fesetenv(&native_)!=0)std::terminate();
#if defined(_M_X64) || defined(__SSE2__)
        _mm_setcsr(csr_);
#endif
        errno=error_;
    }
};
bool frozenAssemblyQualified() noexcept
{
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE) && !defined(__FAST_MATH__) && !defined(_M_FP_FAST) && defined(FE_TONEAREST) && !defined(EIGEN_SCALAR_BINARY_OP_PLUGIN) && EIGEN_WORLD_VERSION==3 && EIGEN_MAJOR_VERSION==4 && EIGEN_MINOR_VERSION==1
    if(!std::numeric_limits<double>::is_iec559||sizeof(double)!=sizeof(std::uint64_t)
        ||std::numeric_limits<double>::radix!=2||std::numeric_limits<double>::digits!=53
        ||std::fegetround()!=FE_TONEAREST)return false;
    unsigned control=0;if(_controlfp_s(&control,0,0)!=0)return false;
    const unsigned csr=_mm_getcsr();
    if((control&_MCW_EM)!=_MCW_EM||(control&_MCW_RC)!=_RC_NEAR
        ||(csr&0x1f80u)!=0x1f80u||(csr&0x6000u)!=0||(csr&0x8040u)!=0)return false;
#ifdef _EM_AMBIGUOUS
    if((control&_EM_AMBIGUOUS)!=0)return false;
#endif
    return true;
#else
    return false;
#endif
}
constexpr std::size_t frozenAssemblyCap=128u*1024u*1024u;
constexpr std::uint32_t frozenNoSlot=std::numeric_limits<std::uint32_t>::max();
struct FrozenAssemblyBudget
{
    std::size_t used=0;
    bool take(std::size_t count,std::size_t width,std::size_t overhead=0) noexcept
    {
        if(used>frozenAssemblyCap||overhead>frozenAssemblyCap-used)return false;
        const auto left=frozenAssemblyCap-used-overhead;
        if(width&&count>left/width)return false;
        used+=overhead+count*width;return true;
    }
};
// Exact requested length: no hidden vector spare-capacity or escaping lease.
struct FrozenIndexArray
{
    std::unique_ptr<std::uint32_t[]> values;std::size_t count=0;
    bool allocate(std::size_t size,FrozenAssemblyBudget& budget)
    {
        if(!budget.take(size,sizeof(std::uint32_t),size?64u:0u))return false;
        if(size)values=std::make_unique<std::uint32_t[]>(size);
        count=size;return true;
    }
    std::uint32_t& operator[](std::size_t i) noexcept{return values[i];}
    std::uint32_t operator[](std::size_t i) const noexcept{return values[i];}
};
bool frozenFinite(double value) noexcept
{
    std::uint64_t bits=0;std::memcpy(&bits,&value,sizeof bits);
    return (bits&0x7ff0000000000000ULL)!=0x7ff0000000000000ULL;
}
bool frozenShape(const FrozenAssemblyMatrix& m,Eigen::Index n) noexcept
{
    return n>=0&&static_cast<std::uint64_t>(n)<frozenNoSlot&&m.rows()==n&&m.cols()==n
        &&m.isCompressed()&&m.nonZeros()>=0&&static_cast<std::uint64_t>(m.nonZeros())<frozenNoSlot
        &&m.outerIndexPtr()&&m.outerIndexPtr()[0]==0&&m.outerIndexPtr()[n]==m.nonZeros();
}
bool frozenStoredMatrix(const FrozenAssemblyMatrix& m,Eigen::Index n) noexcept
{
    if(!frozenShape(m,n))return false;
    for(Eigen::Index row=0;row<n;++row){
        const auto begin=m.outerIndexPtr()[row],end=m.outerIndexPtr()[row+1];
        if(begin<0||end<begin||end>m.nonZeros())return false;
        Eigen::Index previous=-1;
        for(auto q=begin;q<end;++q){const auto column=m.innerIndexPtr()[q];
            if(column<=previous||column>=n||!frozenFinite(m.valuePtr()[q]))return false;previous=column;}
    }
    return true;
}
struct FrozenStoredPattern
{
    std::uint32_t rows=0;FrozenIndexArray outer,inner;
    bool capture(const FrozenAssemblyMatrix& m,FrozenAssemblyBudget& budget,bool validated=false)
    {
        if(!validated&&!frozenStoredMatrix(m,m.rows()))return false;
        rows=static_cast<std::uint32_t>(m.rows());
        if(!outer.allocate(static_cast<std::size_t>(rows)+1,budget)
            ||!inner.allocate(static_cast<std::size_t>(m.nonZeros()),budget))return false;
        for(std::size_t q=0;q<outer.count;++q)outer[q]=static_cast<std::uint32_t>(m.outerIndexPtr()[q]);
        for(std::size_t q=0;q<inner.count;++q)inner[q]=static_cast<std::uint32_t>(m.innerIndexPtr()[q]);
        return true;
    }
    bool matches(const FrozenAssemblyMatrix& m) const noexcept
    {
        if(!frozenShape(m,rows)||static_cast<std::size_t>(m.nonZeros())!=inner.count)return false;
        for(std::size_t q=0;q<outer.count;++q)
            if(m.outerIndexPtr()[q]!=static_cast<Eigen::Index>(outer[q]))return false;
        for(std::size_t q=0;q<inner.count;++q)
            if(m.innerIndexPtr()[q]!=static_cast<Eigen::Index>(inner[q])||!frozenFinite(m.valuePtr()[q]))return false;
        return true;
    }
    void storage(FrozenAssemblyMatrix& result) const
    {
        result.resize(rows,rows);result.resizeNonZeros(static_cast<Eigen::Index>(inner.count));
        // The current row kernel installs indices together with fresh values.
        // No second full CSR copy precedes that same mandatory write.
    }
};
struct FrozenConvectionIndexPlan
{
    FrozenStoredPattern h;
    FrozenIndexArray exceptionOuter,exceptionColumns,exceptionSlots,hSlots;
    std::size_t retainedBytes=0;
};
struct FrozenMatrixIndexPlan
{
    FrozenStoredPattern k,physical;
    FrozenIndexArray kSlots,diagonalSlots,diagonalSource;
    std::unique_ptr<const FrozenConvectionIndexPlan> convection;
    std::size_t fixedBytes=0;
    std::size_t retainedBytes() const noexcept
    {return fixedBytes+(convection?convection->retainedBytes:0u);}
};
struct HydroFrozenMatrixPlanCache
{
    std::mutex mutex;std::unique_ptr<FrozenMatrixIndexPlan> plan;
};
std::unique_ptr<FrozenMatrixIndexPlan> compileFrozenPhysicalPlan(const FrozenAssemblyMatrix& k)
{
    // The caller has checked the complete current K. Its exact stored graph,
    // including explicit zero, unions with the actually required diagonal.
    const auto n=k.rows();std::size_t missing=0;
    for(Eigen::Index row=0;row<n;++row){
        bool diagonal=false;
        for(auto q=k.outerIndexPtr()[row];q<k.outerIndexPtr()[row+1];++q)
            if(k.innerIndexPtr()[q]==row){diagonal=true;break;}
        if(!diagonal)++missing;
    }
    const auto entries=static_cast<std::size_t>(k.nonZeros());
    if(missing>=frozenNoSlot||entries>=frozenNoSlot-missing)return {};
    FrozenAssemblyBudget budget;
    if(!budget.take(1,sizeof(HydroFrozenMatrixPlanCache)+sizeof(FrozenMatrixIndexPlan)+256u))return {};
    auto plan=std::make_unique<FrozenMatrixIndexPlan>();plan->physical.rows=static_cast<std::uint32_t>(n);
    if(!plan->k.capture(k,budget,true)
        ||!plan->physical.outer.allocate(static_cast<std::size_t>(n)+1,budget)
        ||!plan->physical.inner.allocate(entries+missing,budget)
        ||!plan->kSlots.allocate(entries,budget)
        ||!plan->diagonalSlots.allocate(static_cast<std::size_t>(n),budget)
        ||!plan->diagonalSource.allocate(static_cast<std::size_t>(n),budget))return {};
    std::size_t slot=0;
    for(Eigen::Index row=0;row<n;++row){
        plan->physical.outer[static_cast<std::size_t>(row)]=static_cast<std::uint32_t>(slot);bool diagonal=false;
        const auto insertDiagonal=[&](bool fromK){
            plan->diagonalSlots[static_cast<std::size_t>(row)]=static_cast<std::uint32_t>(slot);
            plan->diagonalSource[static_cast<std::size_t>(row)]=fromK?1u:0u;diagonal=true;
            plan->physical.inner[slot++]=static_cast<std::uint32_t>(row);
        };
        for(auto q=k.outerIndexPtr()[row];q<k.outerIndexPtr()[row+1];++q){
            const auto col=k.innerIndexPtr()[q];
            if(!diagonal&&col>row)insertDiagonal(false);
            plan->kSlots[static_cast<std::size_t>(q)]=static_cast<std::uint32_t>(slot);
            if(col==row)insertDiagonal(true);
            else plan->physical.inner[slot++]=static_cast<std::uint32_t>(col);
        }
        if(!diagonal)insertDiagonal(false);
    }
    plan->physical.outer[static_cast<std::size_t>(n)]=static_cast<std::uint32_t>(slot);
    if(slot!=entries+missing)return {};
    plan->fixedBytes=budget.used;return plan;
}
// Search actual current H coordinates only; never insert a potential zero.
std::unique_ptr<const FrozenConvectionIndexPlan> compileFrozenConvectionPlan(
    const FrozenStoredPattern& physical,const FrozenAssemblyMatrix& h,std::size_t liveBytes)
{
    if(!frozenStoredMatrix(h,physical.rows))return {};
    std::size_t exceptions=0;
    for(std::uint32_t row=0;row<physical.rows;++row){
        const auto first=physical.inner.values.get()+physical.outer[row];
        const auto last=physical.inner.values.get()+physical.outer[row+1];
        for(auto hi=h.outerIndexPtr()[row];hi<h.outerIndexPtr()[row+1];++hi){
            const auto col=static_cast<std::uint32_t>(h.innerIndexPtr()[hi]);
            const auto found=std::lower_bound(first,last,col);if(found==last||*found!=col)++exceptions;
        }
    }
    if(exceptions>=frozenNoSlot||physical.inner.count>=frozenNoSlot-exceptions)return {};
    FrozenAssemblyBudget budget;budget.used=liveBytes;
    if(!budget.take(1,sizeof(FrozenConvectionIndexPlan)+256u))return {};
    auto plan=std::make_unique<FrozenConvectionIndexPlan>();
    if(!plan->h.capture(h,budget,true)
        ||!plan->exceptionOuter.allocate(static_cast<std::size_t>(physical.rows)+1,budget)
        ||!plan->exceptionColumns.allocate(exceptions,budget)
        ||!plan->exceptionSlots.allocate(exceptions,budget)
        ||!plan->hSlots.allocate(static_cast<std::size_t>(h.nonZeros()),budget))return {};
    std::size_t extra=0;
    for(std::uint32_t row=0;row<physical.rows;++row){
        plan->exceptionOuter[row]=static_cast<std::uint32_t>(extra);
        const auto first=physical.inner.values.get()+physical.outer[row];
        const auto last=physical.inner.values.get()+physical.outer[row+1];
        for(auto hi=h.outerIndexPtr()[row];hi<h.outerIndexPtr()[row+1];++hi){
            const auto col=static_cast<std::uint32_t>(h.innerIndexPtr()[hi]);const auto found=std::lower_bound(first,last,col);
            const auto pi=static_cast<std::size_t>(found-physical.inner.values.get());
            const auto slot=static_cast<std::uint32_t>(pi+extra);
            if(hi>0&&plan->hSlots[static_cast<std::size_t>(hi)-1]>=slot)return {};
            plan->hSlots[static_cast<std::size_t>(hi)]=slot;
            if(found==last||*found!=col){plan->exceptionColumns[extra]=col;plan->exceptionSlots[extra]=slot;++extra;}
        }
    }
    plan->exceptionOuter[physical.rows]=static_cast<std::uint32_t>(extra);if(extra!=exceptions)return {};
    plan->retainedBytes=budget.used-liveBytes;return plan;
}
bool frozenSameExceptions(const FrozenConvectionIndexPlan* a,const FrozenConvectionIndexPlan& b) noexcept
{
    if(!a)return b.exceptionColumns.count==0;
    if(a->exceptionOuter.count!=b.exceptionOuter.count||a->exceptionColumns.count!=b.exceptionColumns.count)return false;
    for(std::size_t q=0;q<a->exceptionOuter.count;++q)if(a->exceptionOuter[q]!=b.exceptionOuter[q])return false;
    for(std::size_t q=0;q<a->exceptionColumns.count;++q)if(a->exceptionColumns[q]!=b.exceptionColumns[q])return false;
    return true;
}
void frozenConvectionStorage(const FrozenStoredPattern& physical,const FrozenConvectionIndexPlan& route,
    FrozenAssemblyMatrix& result)
{
    result.resize(physical.rows,physical.rows);
    result.resizeNonZeros(static_cast<Eigen::Index>(physical.inner.count+route.exceptionColumns.count));
    // Current union rows install both the exact stored graph and fresh values.
}
bool frozenConvectionMatches(const FrozenStoredPattern& physical,const FrozenConvectionIndexPlan& route,const FrozenAssemblyMatrix& actual) noexcept
{
    if(!frozenStoredMatrix(actual,physical.rows)||static_cast<std::size_t>(actual.nonZeros())!=physical.inner.count+route.exceptionColumns.count)return false;
    for(std::uint32_t row=0;row<physical.rows;++row){
        if(actual.outerIndexPtr()[row]!=static_cast<Eigen::Index>(physical.outer[row])+route.exceptionOuter[row])return false;
        std::size_t pi=physical.outer[row],ei=route.exceptionOuter[row];const auto pe=physical.outer[row+1],ee=route.exceptionOuter[row+1];auto ai=actual.outerIndexPtr()[row];
        while(pi<pe||ei<ee){const auto col=ei==ee||(pi<pe&&physical.inner[pi]<route.exceptionColumns[ei])?physical.inner[pi++]:route.exceptionColumns[ei++];
            if(actual.innerIndexPtr()[ai++]!=static_cast<Eigen::Index>(col))return false;}
    }
    return true;
}
#if defined(_MSC_VER)
#pragma float_control(precise,on,push)
#pragma float_control(except,on)
#pragma fenv_access(on)
#endif
void frozenPhysicalRow(const FrozenMatrixIndexPlan& plan,const FrozenAssemblyMatrix& k,
    FrozenAssemblyMatrix& result,std::size_t row) noexcept
{
    result.outerIndexPtr()[row]=static_cast<Eigen::Index>(plan.physical.outer[row]);
    for(auto q=plan.physical.outer[row];q<plan.physical.outer[row+1];++q)
        result.innerIndexPtr()[q]=static_cast<Eigen::Index>(plan.physical.inner[q]);
    for(auto q=k.outerIndexPtr()[row];q<k.outerIndexPtr()[row+1];++q)
        result.valuePtr()[plan.kSlots[static_cast<std::size_t>(q)]]=.5*k.valuePtr()[q];
}
void frozenDiagonalRow(const FrozenMatrixIndexPlan& plan,const std::vector<double>& masses,double dt,
    FrozenAssemblyMatrix& result,std::size_t row) noexcept
{
    auto& value=result.valuePtr()[plan.diagonalSlots[row]];
    if(!plan.diagonalSource[row])value=0.;
    value+=masses[row]/dt;
}
bool frozenUnionRow(const FrozenStoredPattern& pattern,const FrozenConvectionIndexPlan& route,
    const FrozenAssemblyMatrix& physical,const FrozenAssemblyMatrix& h,
    FrozenAssemblyMatrix& result,std::size_t row) noexcept
{
    auto pi=static_cast<std::size_t>(pattern.outer[row]),ei=static_cast<std::size_t>(route.exceptionOuter[row]);
    const auto pe=pattern.outer[row+1],ee=route.exceptionOuter[row+1];
    auto hi=static_cast<std::size_t>(h.outerIndexPtr()[row]);const auto he=static_cast<std::size_t>(h.outerIndexPtr()[row+1]);
    auto out=pi+ei;result.outerIndexPtr()[row]=static_cast<Eigen::Index>(out);
    const Eigen::internal::scalar_sum_op<double,double> sum;bool finite=true;
    while(pi<pe||ei<ee){
        const bool exception=ei<ee&&(pi==pe||route.exceptionColumns[ei]<pattern.inner[pi]);
        const auto col=exception?route.exceptionColumns[ei++]:pattern.inner[pi];
        const double left=exception?0.:physical.valuePtr()[pi++];double right=0.;
        if(hi<he&&route.hSlots[hi]==out){volatile double half=.5*h.valuePtr()[hi++];right=half;}
        result.innerIndexPtr()[out]=static_cast<Eigen::Index>(col);auto& value=result.valuePtr()[out++];value=sum(left,right);
        finite=frozenFinite(value)&&finite;
    }
    return finite;
}
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE) \
    && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST) \
    && !defined(__AVX__) && !defined(__AVX2__) && !defined(__AVX512F__) && !defined(__FMA__)
struct FrozenValueFloatingState
{
    std::fenv_t native{};unsigned csr=0;bool valid=false;
    bool capture() noexcept {csr=_mm_getcsr();valid=std::fegetenv(&native)==0;return valid;}
    void restore() const noexcept {if(valid&&std::fesetenv(&native)!=0)std::terminate();_mm_setcsr(csr);}
    bool enterFresh(const FrozenValueFloatingState& caller) const noexcept {
        auto fresh=caller.native;fresh._Fe_stat=0;if(std::fesetenv(&fresh)!=0)return false;
        _mm_setcsr(caller.csr&~0x3fu);return true;
    }
    bool sameControl(const FrozenValueFloatingState& other) const noexcept {
        auto a=native,b=other.native;a._Fe_stat=0;b._Fe_stat=0;
        return std::memcmp(&a,&b,sizeof a)==0&&(csr&~0x3fu)==(other.csr&~0x3fu);
    }
};
struct FrozenValueWorker
{
    FrozenValueFloatingState original;unsigned long nativeFlags=0;unsigned csrFlags=0;
};
template<class First,class Second>
bool frozenParallelValues(std::size_t firstCount,std::size_t secondCount,std::size_t liveBytes,
    HydroFrozenMatrixPlanReport& report,First first,Second second) noexcept
{
    const int callerError=errno;const bool previousFallback=report.valuesParallelFallback;
    FrozenValueFloatingState caller;std::size_t workers=0,bytes=0;
    std::unique_ptr<FrozenValueWorker[]> slots;
    {
        FrozenAssemblyFloatingRestore optional;
        if(!optional.captured||!caller.capture()||!frozenAssemblyQualified()||omp_in_parallel()
            ||omp_get_max_threads()<=1)return false;
        workers=std::min(std::max(firstCount,secondCount),static_cast<std::size_t>(omp_get_max_threads()));
        if(workers<2||workers>static_cast<std::size_t>(std::numeric_limits<int>::max()))return false;
        report.valuesParallelFallback=true;
        constexpr std::size_t overhead=sizeof(slots)+64u;
        if(liveBytes>frozenAssemblyCap||overhead>frozenAssemblyCap-liveBytes
            ||workers>(frozenAssemblyCap-liveBytes-overhead)/sizeof(FrozenValueWorker))return false;
        bytes=overhead+workers*sizeof(FrozenValueWorker);
        try{slots=std::make_unique<FrozenValueWorker[]>(workers);}catch(...){return false;}
        report.valuesScratchBytes=std::max(report.valuesScratchBytes,bytes);
    }
    errno=callerError;int failure=0,actual=0;const int team=static_cast<int>(workers);
#pragma omp parallel num_threads(team)
    {
        const int tid=omp_get_thread_num(),workerError=errno;
        if(tid==0)actual=omp_get_num_threads();
        auto& slot=slots[static_cast<std::size_t>(tid)];
        const bool captured=slot.original.capture(),ready=captured&&slot.original.enterFresh(caller);
        if(!ready){
#pragma omp atomic write
            failure=1;
        }
        if(ready)errno=callerError;
#pragma omp for schedule(static)
        for(std::size_t row=0;row<firstCount;++row)if(ready)first(row);
        // Keep the complete original .5*K pass before any M/dt diagonal add.
#pragma omp for schedule(static)
        for(std::size_t row=0;row<secondCount;++row)if(ready)second(row);
        if(ready){FrozenValueFloatingState after;const int arithmeticError=errno;
            if(!after.capture()||arithmeticError!=callerError||!after.sameControl(caller)){
#pragma omp atomic write
                failure=1;
            }else{slot.nativeFlags=after.native._Fe_stat;slot.csrFlags=after.csr&0x3fu;}
        }
        slot.original.restore();errno=workerError;
    }
    if(failure||actual<2){
        {FrozenAssemblyFloatingRestore current;slots.reset();}
        caller.restore();errno=callerError;return false;
    }
    unsigned long nativeFlags=0;unsigned csrFlags=0;
    for(int t=0;t<actual;++t){nativeFlags|=slots[static_cast<std::size_t>(t)].nativeFlags;csrFlags|=slots[static_cast<std::size_t>(t)].csrFlags;}
    auto merged=caller.native;merged._Fe_stat|=nativeFlags;
    if(std::fesetenv(&merged)!=0)std::terminate();_mm_setcsr(caller.csr|csrFlags);errno=callerError;
    {FrozenAssemblyFloatingRestore current;slots.reset();}
    report.valuesParallelUsed=true;report.valuesParallelFallback=previousFallback;
    report.valuesMaximumWorkers=std::max(report.valuesMaximumWorkers,actual);return true;
}
#else
template<class First,class Second>
bool frozenParallelValues(std::size_t,std::size_t,std::size_t,HydroFrozenMatrixPlanReport&,First,Second) noexcept{return false;}
#endif
#if defined(_MSC_VER)
#if !defined(_M_FP_STRICT) || !_M_FP_STRICT
#pragma fenv_access(off)
#endif
#pragma float_control(pop)
#endif
// One exclusive lease. H changes never discard a certified K/P backbone.
class FrozenMatrixAssemblyContext
{
    std::shared_ptr<HydroFrozenMatrixPlanCache> cache_;
    std::unique_lock<std::mutex> lock_;
    HydroFrozenMatrixPlanReport& report_;
    const FrozenMatrixIndexPlan* source_=nullptr;
    const FrozenConvectionIndexPlan* activeRoute_=nullptr;
    std::unique_ptr<const FrozenConvectionIndexPlan> pending_;
    std::unique_ptr<FrozenMatrixIndexPlan> cold_;
    const FrozenAssemblyMatrix* physicalResult_=nullptr;
    bool unionStorageReady_=false;
public:
    FrozenMatrixAssemblyContext(std::shared_ptr<HydroFrozenMatrixPlanCache> cache,
        bool requested,const FrozenAssemblyMatrix& k,Eigen::Index n,HydroFrozenMatrixPlanReport& report) noexcept
        :cache_(std::move(cache)),report_(report)
    {
        FrozenAssemblyFloatingRestore restore;report_.requested=requested;
        report_.viscosityNonZeros=static_cast<std::size_t>(k.nonZeros());
        if(!requested||!restore.captured||!frozenAssemblyQualified()||!frozenShape(k,n))return;
        try{
            if(cache_){lock_=std::unique_lock<std::mutex>(cache_->mutex,std::try_to_lock);
                if(lock_.owns_lock()){source_=cache_->plan.get();if(source_)report_.retainedBytes=source_->retainedBytes();}}
            if(source_&&source_->k.matches(k)){report_.physicalPlanHit=true;report_.eligible=true;return;}
            if(!frozenStoredMatrix(k,n))return;
            report_.eligible=true;
            // A busy foreign lease keeps the original literal path. Otherwise
            // discard the obsolete index slot before constructing its replacement.
            if(!cache_||!lock_.owns_lock())return;
            source_=nullptr;
            if(cache_){cache_->plan.reset();report_.retainedBytes=0;}
            cold_=compileFrozenPhysicalPlan(k);
            if(cold_){source_=cold_.get();report_.physicalColdCompiled=true;
                report_.retainedBytes=source_->retainedBytes();}
        }catch(...){source_=nullptr;}
    }
    ~FrozenMatrixAssemblyContext() noexcept
    {
        FrozenAssemblyFloatingRestore restore;pending_.reset();cold_.reset();if(lock_.owns_lock())lock_.unlock();
    }
    bool physicalStorage(FrozenAssemblyMatrix& result) noexcept
    {
        physicalResult_=nullptr;
        FrozenAssemblyFloatingRestore restore;
        if(!restore.captured||!report_.eligible||!source_||!frozenAssemblyQualified())return false;
        try{source_->physical.storage(result);return true;}catch(...){return false;}
    }
    void physicalValues(const FrozenAssemblyMatrix& k,const std::vector<double>& masses,double dt,FrozenAssemblyMatrix& result)
    {
        const auto count=static_cast<std::size_t>(source_->physical.rows);
        const auto scale=[&](std::size_t row) noexcept {frozenPhysicalRow(*source_,k,result,row);};
        const auto diagonal=[&](std::size_t row) noexcept {frozenDiagonalRow(*source_,masses,dt,result,row);};
        if(!frozenParallelValues(count,count,source_->retainedBytes(),report_,scale,diagonal)){
            for(std::size_t row=0;row<count;++row)scale(row);
            for(std::size_t row=0;row<count;++row)diagonal(row);
        }
        result.outerIndexPtr()[count]=static_cast<Eigen::Index>(source_->physical.inner.count);
        // This local result has not escaped the callback-free assembly chain.
        // Its exact row/column graph was installed above by this same context.
        physicalResult_=&result;
        report_.physicalValuesUsed=true;report_.physicalColdUsed=bool(cold_);
    }
    bool unionStorage(const FrozenAssemblyMatrix& physical,FrozenAssemblyMatrix& result) noexcept
    {
        FrozenAssemblyFloatingRestore restore;
        if(!restore.captured||!report_.eligible||!source_||!frozenAssemblyQualified())return false;
        if(physicalResult_!=&physical&&!source_->physical.matches(physical))return false;
        try{
            if(source_->convection)frozenConvectionStorage(source_->physical,*source_->convection,result);
            else source_->physical.storage(result);
            unionStorageReady_=true;return true;
        }catch(...){return false;}
    }
    bool unionMatches(const FrozenAssemblyMatrix& h,FrozenAssemblyMatrix& storage) noexcept
    {
        FrozenAssemblyFloatingRestore restore;report_.convectionNonZeros=static_cast<std::size_t>(h.nonZeros());
        if(!restore.captured||!report_.eligible||!source_||!frozenAssemblyQualified())return false;
        report_.unionPlanHit=source_->convection&&source_->convection->h.matches(h);
        if(!unionStorageReady_)return false;
        if(report_.unionPlanHit){activeRoute_=source_->convection.get();return true;}
        try{
            pending_=compileFrozenConvectionPlan(source_->physical,h,source_->retainedBytes());
            if(!pending_)return false;
            report_.unionColdCompiled=true;
            if(!frozenSameExceptions(source_->convection.get(),*pending_)){
                frozenConvectionStorage(source_->physical,*pending_,storage);
            }
            activeRoute_=pending_.get();return true;
        }catch(...){pending_.reset();activeRoute_=nullptr;return false;}
    }
    bool unionValues(const FrozenAssemblyMatrix& physical,const FrozenAssemblyMatrix& h,FrozenAssemblyMatrix& result)
    {
        const auto& route=*activeRoute_;const auto count=static_cast<std::size_t>(source_->physical.rows);
        std::atomic<bool> finite{true};
        const auto row=[&](std::size_t i) noexcept {
            if(!frozenUnionRow(source_->physical,route,physical,h,result,i))finite.store(false,std::memory_order_relaxed);
        };
        const auto unused=[](std::size_t) noexcept {};
        const auto live=source_->retainedBytes()+(pending_?pending_->retainedBytes:0u);
        if(!frozenParallelValues(count,0,live,report_,row,unused)){
            finite.store(true,std::memory_order_relaxed);for(std::size_t i=0;i<count;++i)row(i);
        }
        result.outerIndexPtr()[count]=static_cast<Eigen::Index>(source_->physical.inner.count+route.exceptionColumns.count);
        report_.unionValuesUsed=true;report_.unionColdUsed=bool(pending_);
        return finite.load(std::memory_order_relaxed);
    }
    void finish(const FrozenAssemblyMatrix& k,const FrozenAssemblyMatrix& h,
        const FrozenAssemblyMatrix& physical,const FrozenAssemblyMatrix& combined) noexcept
    {
        FrozenAssemblyFloatingRestore restore;
        report_.convectionNonZeros=static_cast<std::size_t>(h.nonZeros());
        report_.physicalNonZeros=static_cast<std::size_t>(physical.nonZeros());
        report_.velocityNonZeros=static_cast<std::size_t>(combined.nonZeros());
        if(!restore.captured||!report_.eligible||!cache_||!lock_.owns_lock()||!frozenAssemblyQualified())return;
        if(cold_){
            if(physicalResult_!=&physical&&!source_->physical.matches(physical))return;
            if(pending_)cold_->convection=std::move(pending_);
            else if(!cold_->convection){
                try{auto route=compileFrozenConvectionPlan(cold_->physical,h,cold_->fixedBytes);
                    if(route&&frozenConvectionMatches(cold_->physical,*route,combined)){
                        report_.unionColdCompiled=true;cold_->convection=std::move(route);}
                }catch(...){ }
            }
            source_=nullptr;activeRoute_=nullptr;cache_->plan=std::move(cold_);
            report_.retainedBytes=cache_->plan->retainedBytes();return;
        }
        if(report_.physicalPlanHit&&source_){
            // Sparse H rebind cannot invalidate this fixed backbone.
            try{
                if(!report_.unionPlanHit){
                    if(!pending_)pending_=compileFrozenConvectionPlan(source_->physical,h,source_->retainedBytes());
                    if(pending_&&(report_.unionValuesUsed||frozenConvectionMatches(source_->physical,*pending_,combined)))
                        cache_->plan->convection=std::move(pending_);
                }
                report_.retainedBytes=cache_->plan->retainedBytes();
            }catch(...){ } // Optional integer route only; keep the fixed plan.
            return;
        }
        // Only a REAL K/P change replaces the one bounded slot.
        source_=nullptr;activeRoute_=nullptr;pending_.reset();cache_->plan.reset();report_.retainedBytes=0;
        try{
            auto candidate=frozenStoredMatrix(k,k.rows())?compileFrozenPhysicalPlan(k):nullptr;
            if(candidate&&candidate->physical.matches(physical)){
                report_.physicalColdCompiled=true;
                auto route=compileFrozenConvectionPlan(candidate->physical,h,candidate->fixedBytes);
                if(route&&frozenConvectionMatches(candidate->physical,*route,combined))candidate->convection=std::move(route);
                report_.retainedBytes=candidate->retainedBytes();cache_->plan=std::move(candidate);
            }
        }catch(...){ } // No original arithmetic, checks, or exceptions enclosed.
    }
};
}
std::shared_ptr<Detail::HydroFrozenMatrixPlanCache> HydroMomentumGeometry::makeFrozenMatrixAssemblyPlans() noexcept
{
    Detail::FrozenAssemblyFloatingRestore restore;if(!restore.captured)return {};
    try{return std::make_shared<Detail::HydroFrozenMatrixPlanCache>();}catch(...){return {};}
}
namespace
{
// New report-only clock reads/conversions preserve the numerical caller's
// flags/controls/errno. Physical operations are never inside this scope.
class PassiveAssemblySampler
{
    using Clock=std::chrono::steady_clock;
    struct FloatingRestore {
        std::fenv_t native{};int error=errno;
#if defined(_M_X64) || defined(__SSE2__)
        unsigned csr=_mm_getcsr();
#endif
        bool captured=std::fegetenv(&native)==0;
        ~FloatingRestore() noexcept {
            if(captured&&std::fesetenv(&native)!=0)std::terminate();
#if defined(_M_X64) || defined(__SSE2__)
            _mm_setcsr(csr);
#endif
            errno=error;
        }
    };
public:
    bool available=true;
    Clock::time_point mark(){
        if(!available)return {};
        FloatingRestore environment;if(!environment.captured){available=false;return {};}
        return Clock::now();
    }
    double seconds(Clock::time_point start){
        if(!available)return 0;
        FloatingRestore environment;if(!environment.captured){available=false;return 0;}
        const double value=std::chrono::duration<double>(Clock::now()-start).count();
        if(!std::isfinite(value)||value<0){available=false;return 0;}
        return value;
    }
};
void require(bool ok,const char* why){if(!ok)throw std::invalid_argument(std::string("[HydroFrozenMidpoint] ")+why);}
bool sameTime(LAMP::TimeLayer a,LAMP::TimeLayer b)
{a.time();b.time();return a.step==b.step&&a.begin==b.begin&&a.end==b.end&&a.position==b.position;}
void validatePhysicalInputs(const QuadtreeMesh& mesh,const HydroDualMaterial& material,const HydroKepConvection& convection,
    const HydroKepConvectionResult& applied,const HydroPhaseViscosity& viscosity,const HydroPressureLoad& pressure,
    const HydroNativeValues& before,const HydroFrozenMidpointBoundary& boundary,HydroMaterialSources sources,
    LAMP::TimeLayer clock)
{
    material.validateFor(mesh,sources,clock);const auto& geometry=material.geometry();const auto& layout=material.layout();
    before.validateFor(layout,mesh);boundary.walls.validateFor(mesh,geometry);
    require(convection.isImpermeable()||convection.isPrescribedFlow(),"requires genuine typed physical convection");
    require(&convection.material().faceRates()==&material.faceRates()&&&applied.material().faceRates()==&material.faceRates(),"different actual material/convection/result owners");
    require(convection.material().geometry().sameIdentity(geometry),"different actual convection geometry");
    applied.boundary().validateFor(mesh,material,sources,clock);
    require(viscosity.geometry().model()==material.material().model(),"different viscosity model");
    if(material.material().model()==HydroPhaseGeometry::Model::GeometricVof){
        require(viscosity.materialModel()==HydroPhaseViscosity::MaterialModel::VofCellAverage||
            (viscosity.materialModel()==HydroPhaseViscosity::MaterialModel::ExactZero&&
             viscosity.liquidViscosity()==0&&viscosity.gasViscosity()==0),
            "VOF frozen midpoint requires genuine inventory P0 or exact-zero phase viscosity");
        viscosity.geometry().validateVofFor(mesh,sources.vofIds[1],clock);
    }else viscosity.validateFor(mesh,sources.midpoint,clock);
    require(&viscosity.geometry().metrics().divergence().matrix()==&geometry.metrics().divergence().matrix()
        &&&pressure.divergence().matrix()==&geometry.metrics().divergence().matrix(),"different authoritative viscosity/pressure D_I");
    require(boundary.normalVelocity.size()==static_cast<std::size_t>(layout.velocityCount()),"different prescribed normal dimensions");
    std::vector<bool> outlet(layout.velocityCount());
    for(const auto& port:boundary.pressureOutlets){require(port.velocityDof>=0&&port.velocityDof<layout.velocityCount()&&!outlet[port.velocityDof],"duplicate or invalid physical pressure outlet");
        const auto face=layout.domain().face(mesh,layout.nativeFace(port.velocityDof));
        require(face.outwardSign!=0&&geometry.metrics().faceMeasures()[port.velocityDof]>0&&std::isfinite(port.pressure),"invalid physical pressure outlet");
        require(!std::binary_search(boundary.walls.velocityDofs().begin(),boundary.walls.velocityDofs().end(),port.velocityDof),"pressure outlet overlaps a physical wall");
        require(boundary.normalVelocity[port.velocityDof]==0,"pressure outlet has a prescribed normal payload");outlet[port.velocityDof]=true;}
    if(material.material().model()==HydroPhaseGeometry::Model::FirstJet){const auto& actual=material.material().midpointRate().spatialAdvector();
        require(actual.boundary()&&actual.boundary()->sameIdentity(boundary.walls),"different actual FirstJet wall owner");}
    std::vector<int> fixed;const auto& normal=before.normalVelocity();
    for(int u=0;u<layout.velocityCount();++u){const auto face=layout.domain().face(mesh,layout.nativeFace(u));require(std::isfinite(boundary.normalVelocity[u]),"nonfinite prescribed normal");
        const bool exterior=face.negativeCell<0||face.positiveCell<0;
        const bool wall=std::binary_search(boundary.walls.velocityDofs().begin(),boundary.walls.velocityDofs().end(),u);
        if(!exterior){require(!outlet[u]&&boundary.normalVelocity[u]==0,"physical boundary payload on an interior normal");continue;}
        if(wall)require(normal[u]==0&&boundary.normalVelocity[u]==0,"nonzero actual impermeable normal");
        if(!outlet[u]){fixed.push_back(u);require(normal[u]==boundary.normalVelocity[u],"begin normal differs from retained prescribed law");
            if(material.material().model()==HydroPhaseGeometry::Model::GeometricVof){const auto& slab=material.material().vofTransport().advector();const auto row=layout.nativeFace(u);
                require(slab.beginNormal()[row]==boundary.normalVelocity[u]&&slab.endNormal()[row]==boundary.normalVelocity[u],
                    "steady prescribed normal differs from original VOF endpoint slab");}
            else require(material.material().midpointRate().spatialAdvector().values().normalVelocity()[u]==boundary.normalVelocity[u],"prescribed normal differs from actual FirstJet source");}
    }
    require(fixed==convection.fixedVelocityDofs(),"physical saddle and KEP boundary partitions differ");
    if(convection.isImpermeable())require(boundary.pressureOutlets.empty()&&fixed==boundary.walls.velocityDofs(),"closed convection has open physical boundary data");
    const auto divergence=geometry.metrics().divergence().apply(mesh,before);const auto& di=geometry.metrics().divergence().matrix();
    for(int p=0;p<layout.pressureCount();++p){long double scale=0;for(HydroIntegralDivergence::Matrix::InnerIterator it(di,p);it;++it)scale+=std::abs(static_cast<long double>(it.value())*normal[it.col()]);
        require(scale>0?std::abs(divergence[p])/scale<=1e-11L:divergence[p]==0,"begin native velocity is not divergence free");}
    if(boundary.pressureOutlets.empty()){long double mean=0,gross=0;for(int p=0;p<layout.pressureCount();++p){const auto value=static_cast<long double>(geometry.metrics().cellVolumes()[p])*before.pressure()[p];mean+=value;gross+=std::abs(value);}
        require(gross>0?std::abs(mean)/gross<=1e-11L:mean==0,"begin pressure gauge is not mean zero");}
}
}
HydroFrozenMidpointBoundary HydroFrozenMidpointBoundary::closed(HydroAdvectorBoundary walls)
{return {walls,std::vector<double>(walls.momentumGeometry().layout().velocityCount()),{}};}
HydroFrozenMidpointBoundary HydroFrozenMidpointBoundary::prescribed(HydroAdvectorBoundary walls,std::vector<double> normal,std::vector<HydroSaddlePressureOutlet> outlets)
{return {std::move(walls),std::move(normal),std::move(outlets)};}
struct HydroFrozenMidpoint::Data
{
    Data(HydroDualMaterial d,HydroKepConvection c,HydroKepConvectionResult r,HydroPhaseViscosity v,HydroPressureLoad p,
        HydroNativeValues b,HydroFrozenMidpointForces f,HydroFrozenMidpointBoundary w,HydroFrozenMidpointScheme s)
        :material(std::move(d)),convection(std::move(c)),applied(std::move(r)),viscosity(std::move(v)),pressure(std::move(p)),before(std::move(b)),forces(std::move(f)),boundary(std::move(w)),scheme(s){}
    HydroDualMaterial material;HydroKepConvection convection;HydroKepConvectionResult applied;HydroPhaseViscosity viscosity;
    HydroPressureLoad pressure;HydroNativeValues before;HydroFrozenMidpointForces forces;HydroFrozenMidpointBoundary boundary;HydroFrozenMidpointScheme scheme;
    Matrix physicalA,h;Eigen::VectorXd ku0,hu0;
    HydroFrozenMatrixPlanReport matrixPlan;
    std::vector<long double> bodyLoads;std::vector<double> interfaceRhs,maxwellRhs,rhs;
    std::optional<HydroSinglePhaseSaddle> saddle;std::shared_ptr<const HydroSaddleResult> solved;
    std::size_t solveCount=0;
};
HydroFrozenMidpoint::HydroFrozenMidpoint(std::shared_ptr<const Data> data):data_(std::move(data)){}
std::shared_ptr<HydroFrozenMidpoint::Data> HydroFrozenMidpoint::assemble(const QuadtreeMesh& mesh,HydroDualMaterial material,HydroKepConvection convection,
    HydroKepConvectionResult applied,HydroPhaseViscosity viscosity,HydroPressureLoad pressure,HydroNativeValues before,
    HydroFrozenMidpointForces forces,HydroFrozenMidpointBoundary boundary,HydroFrozenMidpointScheme scheme,
    HydroMaterialSources sources,LAMP::TimeLayer clock,HydroSaddleConfig saddleConfig)
{
    PassiveAssemblySampler sampler;HydroSaddlePreparationTimings preparation;const auto inputsStarted=sampler.mark();
    require(scheme==HydroFrozenMidpointScheme::LaggedConvection||scheme==HydroFrozenMidpointScheme::FrozenConvection,"invalid midpoint convection scheme");
    validatePhysicalInputs(mesh,material,convection,applied,viscosity,pressure,before,boundary,sources,clock);
    require(std::isfinite(forces.bodyAcceleration.x)&&std::isfinite(forces.bodyAcceleration.y),"nonfinite physical body acceleration");
    require(!material.layout().domain().axisymmetric()||forces.bodyAcceleration.x==0,"RZ physical body acceleration must be axial");
    auto data=std::make_shared<Data>(std::move(material),std::move(convection),std::move(applied),std::move(viscosity),std::move(pressure),std::move(before),std::move(forces),std::move(boundary),scheme);
    const auto& dual=data->material;const auto& l=dual.layout();const auto dt=dual.deltaTime();const bool linear=scheme==HydroFrozenMidpointScheme::FrozenConvection;
    preparation.inputValidationSeconds=sampler.seconds(inputsStarted);const auto forcesStarted=sampler.mark();
    std::visit([&](const auto& load){using T=std::decay_t<decltype(load)>;
        if constexpr(std::is_same_v<T,std::monostate>){}
        else if constexpr(std::is_same_v<T,HydroBalancedCsfLoad>){
            if(dual.material().model()==HydroPhaseGeometry::Model::GeometricVof)
                load.validateVofFor(mesh,data->pressure,sources.vofIds[1],clock);
            else{require(dual.material().model()==HydroPhaseGeometry::Model::FirstJet,"CSF owner requires the FirstJet material model");
                load.validateFor(mesh,data->pressure,sources.midpoint,clock);}
            require(&load.phaseGeometry().cellVolumes()==&dual.material().midpointGeometry().cellVolumes(),"different actual CSF midpoint owner");data->interfaceRhs=load.momentumRhsLoad();
        }else if constexpr(std::is_same_v<T,HydroInterfaceNormalLoad>){
            require(dual.material().model()==HydroPhaseGeometry::Model::FirstJet,"sharp FirstJet owner requires the FirstJet material model");
            load.validateFor(mesh,data->pressure,sources.midpoint,clock);require(&load.geometry().cellVolumes()==&dual.material().midpointGeometry().cellVolumes(),"different actual sharp midpoint owner");
            data->interfaceRhs=load.balanceLoads();for(double& value:data->interfaceRhs)value=-value;
        }
    },data->forces.interfaceLoad);
    const bool capillaryActive=!std::holds_alternative<std::monostate>(data->forces.interfaceLoad);
    if(capillaryActive)require(data->interfaceRhs.size()==static_cast<std::size_t>(l.velocityCount()),"different actual interface force dimensions");
    const auto* capillaryRhs=capillaryActive?&data->interfaceRhs:nullptr;
    if(data->forces.maxwellLoad){
        require(dual.material().model()==HydroPhaseGeometry::Model::GeometricVof,"Maxwell owner requires actual VOF material");
        const auto& load=*data->forces.maxwellLoad;
        if(load.isSharp())load.validateSharpFor(mesh,data->pressure,dual.material().midpointGeometry(),load.sharpStress(),clock);
        else load.validateFor(mesh,data->pressure,dual.material().midpointGeometry(),load.stress(),clock);
        data->maxwellRhs=load.momentumRhsLoad();
        require(data->maxwellRhs.size()==static_cast<std::size_t>(l.velocityCount()),"different actual Maxwell force dimensions");
    }
    const auto* maxwellRhs=data->forces.maxwellLoad?&data->maxwellRhs:nullptr;
    preparation.forceValidationSeconds=sampler.seconds(forcesStarted);const auto matrixStarted=sampler.mark();
    // Original TwoPhaseStep matrix/RHS arithmetic and operation order.
    const Matrix& k=data->viscosity.viscosity().matrix();require(data->viscosity.liquidViscosity()!=0||data->viscosity.gasViscosity()!=0||k.nonZeros()==0,"mu zero viscosity is not exactly zero");
    Detail::FrozenMatrixAssemblyContext matrixPlan(linear?dual.endMass().geometry().frozenMatrixAssemblyPlans():std::shared_ptr<Detail::HydroFrozenMatrixPlanCache>{},
        linear,k,l.velocityCount(),data->matrixPlan);
    Matrix physicalA;
    if(matrixPlan.physicalStorage(physicalA))matrixPlan.physicalValues(k,dual.endMass().masses(),dt,physicalA);
    else{
        physicalA=.5*k;for(int u=0;u<l.velocityCount();++u)physicalA.coeffRef(u,u)+=dual.endMass().masses()[u]/dt;physicalA.makeCompressed();
    }
    Matrix a;bool currentValuesFinite=false;const bool preparedUnion=matrixPlan.unionStorage(physicalA,a);
    if(!preparedUnion)a=physicalA;
    Matrix h(l.velocityCount(),l.velocityCount());
    const Eigen::Map<const Eigen::VectorXd> u0(data->before.normalVelocity().data(),l.velocityCount());const Eigen::VectorXd ku0=k*u0;
    Eigen::VectorXd hu0=Eigen::VectorXd::Zero(l.velocityCount());
    std::vector<long double> boundaryLoad;
    if(linear){auto frozen=data->convection.prepareFrozenMomentum(mesh,dual,dual.sources(),dual.timeLayer());
        h=std::move(frozen.first);boundaryLoad=std::move(frozen.second);
        const bool unionMatches=matrixPlan.unionMatches(h,a);
        if(preparedUnion&&unionMatches)currentValuesFinite=matrixPlan.unionValues(physicalA,h,a);
        else{if(preparedUnion)a=physicalA;a+=.5*h;}
        a.makeCompressed();hu0=h*u0;}
    preparation.matrixActionSeconds=sampler.seconds(matrixStarted);const auto rhsStarted=sampler.mark();
    if(!linear)boundaryLoad.resize(l.velocityCount());
    const auto acceleration=data->forces.bodyAcceleration;const bool bodyActive=acceleration.x!=0||acceleration.y!=0;
    std::vector<long double> bodyLoads(l.velocityCount());
    if(bodyActive)for(int u=0;u<l.velocityCount();++u){
        const double g=l.normalAxis(u)==HydroDofLayout::NormalAxis::First?acceleration.x:acceleration.y;
        if(g==0)continue;const auto mass=dual.midpointMass().masses()[u];
        const long double load=static_cast<long double>(mass)*g;const double represented=static_cast<double>(load);
        require(std::isfinite(load)&&std::isfinite(represented)&&(mass==0||load!=0)&&(represented!=0||load==0),"unrepresentable body load");bodyLoads[u]=load;
    }
    std::vector<double> rhs(l.velocityCount());
    for(int u=0;u<l.velocityCount();++u){long double value=static_cast<long double>(dual.beginMass().masses()[u])*u0[u]/dt-.5L*ku0[u];
        if(linear){value-=.5L*hu0[u];value+=boundaryLoad[u];}else value+=data->applied.momentumLoads()[u];
        if(capillaryRhs)value+=(*capillaryRhs)[u];if(maxwellRhs)value+=(*maxwellRhs)[u];if(bodyActive)value+=bodyLoads[u];rhs[u]=static_cast<double>(value);
        require(std::isfinite(rhs[u])&&(rhs[u]!=0||value==0),"unrepresentable two-phase momentum RHS");}
    preparation.rhsSeconds=sampler.seconds(rhsStarted);const auto velocityStarted=sampler.mark();
    // The actual momentum matrix already represents the physical equation.
    // Its sole immutable owner also serves the public Frozen matrix view.
    const auto velocityBlock=HydroSaddleVelocityBlock::buildFrozenActual(l,std::move(a),currentValuesFinite);
    preparation.velocityPrepareSeconds=sampler.seconds(velocityStarted);preparation.available=sampler.available;
    if(dual.usesConservativeVofMomentum()){
        // Match the existing momentum refinement target to physical loads,
        // rather than accepting early against a large transported background.
        // Only the inner target is tightened; all original physical gates and
        // the full matrix/RHS remain unchanged. Avoid demanding epsilon when
        // the endpoint background is zero or the physical load dominates it.
        long double inertiaScale=0,loadScale=0;
        for(int u=0;u<l.velocityCount();++u){
            inertiaScale=std::hypot(inertiaScale,static_cast<long double>(dual.beginMass().masses()[u])*u0[u]/dt);
            loadScale=std::hypot(loadScale,data->applied.momentumLoads()[u]);
            loadScale=std::hypot(loadScale,static_cast<long double>(ku0[u]));
            loadScale=std::hypot(loadScale,bodyLoads[u]);
            if(capillaryRhs)loadScale=std::hypot(loadScale,static_cast<long double>((*capillaryRhs)[u]));
            if(maxwellRhs)loadScale=std::hypot(loadScale,static_cast<long double>((*maxwellRhs)[u]));
        }
        if(inertiaScale>0&&loadScale<inertiaScale){
            const auto target=saddleConfig.momentumRelativeTolerance*(loadScale/inertiaScale);
            saddleConfig.momentumRelativeTolerance=static_cast<double>(std::max(
                static_cast<long double>(std::numeric_limits<double>::epsilon()),target));
        }
    }
    const auto saddle=data->boundary.pressureOutlets.empty()?HydroSinglePhaseSaddle::buildClosed(data->pressure,velocityBlock,saddleConfig,preparation,dual.endMass().geometry())
        :HydroSinglePhaseSaddle::buildWithPressureOutlets(data->pressure,velocityBlock,data->boundary.pressureOutlets,saddleConfig,preparation,dual.endMass().geometry());

    matrixPlan.finish(k,h,physicalA,velocityBlock.matrix());
    data->physicalA.swap(physicalA);data->h.swap(h);data->ku0=ku0;data->hu0=hu0;
    data->bodyLoads=std::move(bodyLoads);data->rhs=std::move(rhs);data->saddle=saddle;
    return data;
}
HydroFrozenMidpoint HydroFrozenMidpoint::prepare(const QuadtreeMesh& mesh,HydroDualMaterial material,HydroKepConvection convection,
    HydroKepConvectionResult applied,HydroPhaseViscosity viscosity,HydroPressureLoad pressure,HydroNativeValues before,
    HydroFrozenMidpointForces forces,HydroFrozenMidpointBoundary boundary,HydroFrozenMidpointScheme scheme,
    HydroMaterialSources sources,LAMP::TimeLayer clock,HydroSaddleConfig saddleConfig,const HydroNativeValues* initialGuess)
{
    require(!initialGuess||saddleConfig.gmresOnlyParallel,"native initial guess requires explicit GMRES-only solve");
    auto data=assemble(mesh,std::move(material),std::move(convection),std::move(applied),std::move(viscosity),
        std::move(pressure),std::move(before),std::move(forces),std::move(boundary),scheme,sources,clock,saddleConfig);
    const HydroSaddleRhs rhs(data->material.layout(),data->rhs,data->boundary.normalVelocity);
    data->solved=std::make_shared<const HydroSaddleResult>(saddleConfig.gmresOnlyParallel
        ?data->saddle->solve(mesh,rhs,initialGuess?*initialGuess:data->before)
        :data->saddle->solve(mesh,rhs));
    data->solveCount=1;
    return HydroFrozenMidpoint(std::move(data));
}
std::size_t HydroFrozenMidpoint::solveCount() const noexcept{return data_->solveCount;}
const HydroDualMaterial& HydroFrozenMidpoint::material() const noexcept{return data_->material;}
const HydroKepConvection& HydroFrozenMidpoint::convection() const noexcept{return data_->convection;}
const HydroKepConvectionResult& HydroFrozenMidpoint::appliedConvection() const noexcept{return data_->applied;}
const HydroPhaseViscosity& HydroFrozenMidpoint::viscosity() const noexcept{return data_->viscosity;}
const HydroPressureLoad& HydroFrozenMidpoint::pressureLoad() const noexcept{return data_->pressure;}
const HydroNativeValues& HydroFrozenMidpoint::beginValues() const noexcept{return data_->before;}
const HydroFrozenMidpointForces& HydroFrozenMidpoint::forces() const noexcept{return data_->forces;}
const HydroFrozenMidpointBoundary& HydroFrozenMidpoint::boundary() const noexcept{return data_->boundary;}
HydroFrozenMidpointScheme HydroFrozenMidpoint::scheme() const noexcept{return data_->scheme;}
const HydroFrozenMatrixPlanReport& HydroFrozenMidpoint::matrixAssemblyPlanReport() const noexcept{return data_->matrixPlan;}
const HydroFrozenMidpoint::Matrix& HydroFrozenMidpoint::viscosityMatrix() const noexcept{return data_->viscosity.viscosity().matrix();}
const HydroFrozenMidpoint::Matrix& HydroFrozenMidpoint::physicalVelocityMatrix() const noexcept{return data_->physicalA;}
const HydroFrozenMidpoint::Matrix& HydroFrozenMidpoint::velocityMatrix() const noexcept{return data_->saddle->velocityBlock().matrix();}
const HydroFrozenMidpoint::Matrix& HydroFrozenMidpoint::frozenConvectionMatrix() const noexcept{return data_->h;}
const Eigen::VectorXd& HydroFrozenMidpoint::beginViscousAction() const noexcept{return data_->ku0;}
const Eigen::VectorXd& HydroFrozenMidpoint::beginFrozenAction() const noexcept{return data_->hu0;}
const std::vector<long double>& HydroFrozenMidpoint::bodyLoads() const noexcept{return data_->bodyLoads;}
const std::vector<double>& HydroFrozenMidpoint::interfaceLoads() const noexcept{return data_->interfaceRhs;}
bool HydroFrozenMidpoint::hasInterfaceLoad() const noexcept{return !std::holds_alternative<std::monostate>(data_->forces.interfaceLoad);}
const std::vector<double>& HydroFrozenMidpoint::maxwellLoads() const noexcept{return data_->maxwellRhs;}
bool HydroFrozenMidpoint::hasMaxwellLoad() const noexcept{return data_->forces.maxwellLoad.has_value();}
const std::vector<double>& HydroFrozenMidpoint::integralMomentumRhs() const noexcept{return data_->rhs;}
const HydroSinglePhaseSaddle& HydroFrozenMidpoint::saddle() const noexcept{return *data_->saddle;}
const HydroSaddleResult& HydroFrozenMidpoint::solveResult() const noexcept{return *data_->solved;}
void HydroFrozenMidpoint::validateFor(const QuadtreeMesh& mesh,HydroMaterialSources sources,LAMP::TimeLayer clock) const
{material().validateFor(mesh,sources,clock);require(sameTime(clock,material().timeLayer()),"different actual frozen midpoint clock");}
}


