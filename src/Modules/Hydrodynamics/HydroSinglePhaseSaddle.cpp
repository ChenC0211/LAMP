#include "Modules/Hydrodynamics/HydroSinglePhaseSaddle.h"
#include "Modules/Hydrodynamics/HydroEnergyViscosity.h"
#include "Modules/Hydrodynamics/HydroMomentumGeometry.h"
#include "Modules/Hydrodynamics/Detail/HydroPressureCoarseCorrection.h"
#include "Solvers/Linear/ParallelBlockJacobi.h"

#include <Eigen/SparseCholesky>
#include <Eigen/SuperLUSupport>
#include <chrono>
#include <atomic>
#include <cfenv>
#include <cstdint>
#include <cstring>
#include <limits>
#include <mutex>
#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <string>
#include <sstream>
#include <utility>
#include <type_traits>
#include <cerrno>
#include <exception>
#include <optional>
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
#include <omp.h>
#endif
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
#include <float.h>
#include <xmmintrin.h>
#endif
#if defined(_M_X64) || defined(__SSE2__)
#include <xmmintrin.h>
#endif

namespace LAMP::Hydrodynamics
{
namespace
{
using Sparse = HydroSinglePhaseSaddle::Matrix;
using DirectSparse = Eigen::SparseMatrix<double,Eigen::ColMajor,int>;
using WallClock = std::chrono::steady_clock;
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
    // New solve sub-timers must not introduce the first unmasked FP trap.
    void requireMaskedExceptions() noexcept {
        const int error=errno;
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
        unsigned control=0;
        if(_controlfp_s(&control,0,0)!=0||(control&_MCW_EM)!=_MCW_EM
            ||(_mm_getcsr()&0x1f80u)!=0x1f80u)available=false;
#else
        available=false;
#endif
        errno=error;
    }
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
    void addSeconds(Clock::time_point start,double& total){
        if(!available)return;
        FloatingRestore environment;if(!environment.captured){available=false;return;}
        const double value=std::chrono::duration<double>(Clock::now()-start).count(),next=total+value;
        if(!std::isfinite(value)||value<0||!std::isfinite(next)||next<total){available=false;return;}
        total=next;
    }
};
double elapsedSeconds(WallClock::time_point start)
{return std::chrono::duration<double>(WallClock::now()-start).count();}
// Compensated residual of the SAME binary64 augmented system. At nearly
// stationary states, ordinary dot products can lose tiny local continuity
// residuals even though the global force backward error is already small.
Eigen::VectorXd correctionResidual(const Sparse& matrix,const Eigen::VectorXd& x,
    const Eigen::VectorXd& rhs)
{
    Eigen::VectorXd high=rhs,low=Eigen::VectorXd::Zero(rhs.size());
    const auto add=[&](Eigen::Index row,double value){
        const double previous=high[row],next=previous+value;
        low[row]+=std::abs(previous)>=std::abs(value)?(previous-next)+value:(value-next)+previous;
        high[row]=next;
    };
    for(Eigen::Index col=0;col<matrix.outerSize();++col)
        for(Sparse::InnerIterator it(matrix,col);it;++it){
            const double product=-it.value()*x[it.col()];
            add(it.row(),product);
            add(it.row(),std::fma(-it.value(),x[it.col()],-product));
        }
    return high+low;
}
void require(bool ok, const char* reason)
{
    if (!ok) throw std::invalid_argument(std::string("[HydroSinglePhaseSaddle] ")+reason);
}
// The public operator retains its original Eigen::Index storage. SuperLU
// requires int indices: check before conversion, retain every stored entry,
// and verify the binary64 value array including signed and explicit zeros.
DirectSparse directQuotient(const Sparse& matrix)
{
    const auto limit=static_cast<Eigen::Index>(std::numeric_limits<int>::max());
    require(matrix.isCompressed(),"SuperLU quotient must be compressed");
    require(matrix.rows()>=0&&matrix.cols()>=0&&matrix.nonZeros()>=0
        &&matrix.rows()<=limit&&matrix.cols()<=limit&&matrix.nonZeros()<=limit,
        "SuperLU quotient exceeds checked int index range");
    for(Eigen::Index col=0;col<=matrix.outerSize();++col)
        require(matrix.outerIndexPtr()[col]>=0&&matrix.outerIndexPtr()[col]<=limit,
            "SuperLU quotient outer index exceeds checked int range");
    for(Eigen::Index entry=0;entry<matrix.nonZeros();++entry)
        require(matrix.innerIndexPtr()[entry]>=0&&matrix.innerIndexPtr()[entry]<matrix.rows(),
            "SuperLU quotient inner index is invalid");
    DirectSparse result=matrix;result.makeCompressed();
    require(result.rows()==matrix.rows()&&result.cols()==matrix.cols()
        &&result.nonZeros()==matrix.nonZeros(),"SuperLU quotient conversion changed dimensions or stored entries");
    for(Eigen::Index col=0;col<=matrix.outerSize();++col)
        require(result.outerIndexPtr()[col]==matrix.outerIndexPtr()[col],
            "SuperLU quotient conversion changed outer indices");
    for(Eigen::Index entry=0;entry<matrix.nonZeros();++entry)
        require(result.innerIndexPtr()[entry]==matrix.innerIndexPtr()[entry],
            "SuperLU quotient conversion changed inner indices");
    require(matrix.nonZeros()==0||std::memcmp(result.valuePtr(),matrix.valuePtr(),
        static_cast<std::size_t>(matrix.nonZeros())*sizeof(double))==0,
        "SuperLU quotient conversion changed binary64 coefficients");
    // All ABI conversion checks above precede this private factor-only
    // compaction. Retain every nonzero bit pattern; remove only IEEE +/-0.
    static_assert(sizeof(double)==sizeof(std::uint64_t) && std::numeric_limits<double>::is_iec559,
        "SuperLU exact-zero compaction requires binary64");
    const auto magnitudeBits=[](const double& value) {
        std::uint64_t word=0;std::memcpy(&word,&value,sizeof(word));
        return word&UINT64_C(0x7fffffffffffffff);
    };
    result.prune([&](Eigen::Index,Eigen::Index,const double& value) {return magnitudeBits(value)!=0;});
    require(result.isCompressed() && result.rows()==matrix.rows() && result.cols()==matrix.cols(),
        "SuperLU exact-zero compaction changed dimensions or compression");
    Eigen::Index retained=0;
    for(Eigen::Index col=0;col<matrix.outerSize();++col) {
        int cursor=result.outerIndexPtr()[col];
        const int end=result.outerIndexPtr()[col+1];
        for(Eigen::Index entry=matrix.outerIndexPtr()[col];entry<matrix.outerIndexPtr()[col+1];++entry) {
            if(magnitudeBits(matrix.valuePtr()[entry])==0)continue;
            require(cursor<end && result.innerIndexPtr()[cursor]==matrix.innerIndexPtr()[entry],
                "SuperLU exact-zero compaction changed nonzero row or column");
            require(std::memcmp(result.valuePtr()+cursor,matrix.valuePtr()+entry,sizeof(double))==0,
                "SuperLU exact-zero compaction changed nonzero binary64 coefficient");
            ++cursor;++retained;
        }
        require(cursor==end,"SuperLU exact-zero compaction inserted nonzero entries");
    }
    require(retained==result.nonZeros(),"SuperLU exact-zero compaction changed retained entry count");
    return result;
}
// Complete per-factory input certificate, not a same-pattern or old-value cache.
// On unsupported Eigen/plugin/compiler configurations the original eager path
// below remains literal. Source guards retain their original order and text.
bool deferredDirectInputCertified(const Sparse& matrix)
{
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE) \
    && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST) \
    && EIGEN_WORLD_VERSION==3 && EIGEN_MAJOR_VERSION==4 && EIGEN_MINOR_VERSION==1 \
    && !defined(EIGEN_SPARSEMATRIX_PLUGIN) && !defined(EIGEN_SPARSEMATRIXBASE_PLUGIN) \
    && !defined(EIGEN_SPARSE_CREATE_TEMPORARY_PLUGIN) \
    && !defined(EIGEN_SPARSE_TRANSPOSED_COPY_PLUGIN) \
    && !defined(EIGEN_SPARSE_COMPRESSED_STORAGE_REALLOCATE_PLUGIN) \
    && !defined(EIGEN_SCALAR_BINARY_OP_PLUGIN)
    static_assert(sizeof(double)==sizeof(std::uint64_t) && std::numeric_limits<double>::is_iec559,
        "deferred SuperLU input requires binary64");
    const auto limit=static_cast<Eigen::Index>(std::numeric_limits<int>::max());
    require(matrix.isCompressed(),"SuperLU quotient must be compressed");
    require(matrix.rows()>=0&&matrix.cols()>=0&&matrix.nonZeros()>=0
        &&matrix.rows()<=limit&&matrix.cols()<=limit&&matrix.nonZeros()<=limit,
        "SuperLU quotient exceeds checked int index range");
    for(Eigen::Index col=0;col<=matrix.outerSize();++col)
        require(matrix.outerIndexPtr()[col]>=0&&matrix.outerIndexPtr()[col]<=limit,
            "SuperLU quotient outer index exceeds checked int range");
    for(Eigen::Index entry=0;entry<matrix.nonZeros();++entry)
        require(matrix.innerIndexPtr()[entry]>=0&&matrix.innerIndexPtr()[entry]<matrix.rows(),
            "SuperLU quotient inner index is invalid");
    // Extra canonical qualification is optional: failure executes the original
    // complete conversion, rather than introducing another input rejection.
    const auto count=matrix.nonZeros();
    if(matrix.outerIndexPtr()[0]!=0||matrix.outerIndexPtr()[matrix.outerSize()]!=count)return false;
    for(Eigen::Index col=0;col<=matrix.outerSize();++col){
        const auto value=matrix.outerIndexPtr()[col];
        if(static_cast<Eigen::Index>(static_cast<int>(value))!=value)return false;
    }
    int retained=0;
    for(Eigen::Index col=0;col<matrix.outerSize();++col){
        const auto first=matrix.outerIndexPtr()[col],last=matrix.outerIndexPtr()[col+1];
        if(first>last||last>count)return false;
        Eigen::Index previous=-1;
        for(Eigen::Index entry=first;entry<last;++entry){
            const auto row=matrix.innerIndexPtr()[entry];
            if(row<=previous||static_cast<Eigen::Index>(static_cast<int>(row))!=row)return false;
            previous=row;
            // Audited same-double assignment preserves these exact bytes;
            // no finite classification or floating zero comparison is used.
            std::uint64_t word=0;double shadow=0.;
            std::memcpy(&word,matrix.valuePtr()+entry,sizeof(word));
            std::memcpy(&shadow,&word,sizeof(shadow));
            if(std::memcmp(&shadow,matrix.valuePtr()+entry,sizeof(shadow))!=0)return false;
            if((word&UINT64_C(0x7fffffffffffffff))!=0)++retained;
        }
    }
    // Stable pruning can only remove +/-0. Counts remain int-representable;
    // actual Direct fallback reruns the unchanged conversion and verification.
    return static_cast<Eigen::Index>(retained)<=count;
#else
    (void)matrix;return false;
#endif
}
struct QuotientSolveTrace
{
    int gmresSolves=0,directSolves=0,iterations=0,restarts=0,fallbackReason=0;
    int lastSolveIterations=0,lastSolveBudget=0,lastSolveRestartLength=0;
    std::size_t pressureCoarseApplications=0,pressureGlobalSolveCalls=0;
    double pressureGlobalSolveSeconds=0;bool pressureGlobalTimingsAvailable=true;
    const char* lastSolveExit="not_started";
    bool directOnly=false,stagnated=false;
    std::size_t orderedSpmvCalls=0,orderedSpmvParallelCalls=0;
    int orderedSpmvMaximumWorkers=0;bool orderedSpmvFallback=false;
    double lastRelativeResidual=0,lastScaledRelativeResidual=0,lazyFactorSeconds=0;
    bool solveBreakdownTimingsAvailable=true;
    double krylovPreconditionApplySeconds=0,krylovSpmvSeconds=0,krylovOrthogonalizationSeconds=0;
    double krylovTrueResidualSeconds=0,physicalAssessmentSeconds=0;
    int orthogonalizationParallelCalls=0,orthogonalizationMaximumWorkers=0;
};
// One numerical CGS2 program. Each pass freezes all complete column dots,
// then subtracts one ascending-column sum per row. Team size changes only
// ownership of independent columns/rows, never their arithmetic order.
constexpr Eigen::Index maximumCgsColumns=80,cgsTileRows=128;
enum class Cgs2Execution {NotStarted,Completed,FailedAfterStart};
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
#pragma float_control(precise,on,push)
#pragma float_control(except,on)
#pragma fenv_access(on)
#endif
bool cgs2Projection(const Eigen::MatrixXd& basis,const Eigen::VectorXd& work,
    Eigen::MatrixXd& h,Eigen::Index column,Eigen::Index k,double* projections) noexcept
{
    const double projection=basis.col(k).dot(work);
    if(!std::isfinite(projection))return false;
    projections[k]=projection;h(k,column)+=projection;return true;
}
void cgs2UpdateTile(const Eigen::MatrixXd& basis,Eigen::VectorXd& work,
    Eigen::Index column,Eigen::Index first,const double* projections,double* sums) noexcept
{
    const auto count=std::min(cgsTileRows,work.size()-first);
    for(Eigen::Index i=0;i<count;++i)sums[i]=0.;
    for(Eigen::Index k=0;k<=column;++k)
        for(Eigen::Index i=0;i<count;++i){volatile double product=basis(first+i,k)*projections[k];sums[i]+=product;}
    for(Eigen::Index i=0;i<count;++i)work[first+i]-=sums[i];
}
bool serialCgs2(const Eigen::MatrixXd& basis,Eigen::VectorXd& work,
    Eigen::MatrixXd& h,Eigen::Index column) noexcept
{
    std::array<double,static_cast<std::size_t>(maximumCgsColumns)> projections;
    std::array<double,static_cast<std::size_t>(cgsTileRows)> sums;
    const auto tiles=work.size()>0?(work.size()-1)/cgsTileRows+1:0;
    for(int pass=0;pass<2;++pass){
        for(Eigen::Index k=0;k<=column;++k)
            if(!cgs2Projection(basis,work,h,column,k,projections.data()))return false;
        for(Eigen::Index tile=0;tile<tiles;++tile)cgs2UpdateTile(basis,work,column,tile*cgsTileRows,projections.data(),sums.data());
    }
    return work.allFinite();
}
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
#if !defined(_M_FP_STRICT) || !_M_FP_STRICT
#pragma fenv_access(off)
#endif
#pragma float_control(pop)
#endif
// Optional CURRENT-factor execution layout. No old physical coefficient owner,
// cross-factor registry, factor reuse or changed Krylov acceptance criterion.
class OrderedSpmvOptionalRestore
{
    std::fenv_t native_{};int error_=errno;
#if defined(_M_X64) || defined(__SSE2__)
    unsigned csr_=_mm_getcsr();
#endif
    bool captured_=std::fegetenv(&native_)==0;
public:
    ~OrderedSpmvOptionalRestore() noexcept {
        if(captured_&&std::fesetenv(&native_)!=0)std::terminate();
#if defined(_M_X64) || defined(__SSE2__)
        _mm_setcsr(csr_);
#endif
        errno=error_;
    }
};
class OrderedQuotientRows;
struct OrderedQuotientRowsDeleter {void operator()(OrderedQuotientRows*) const noexcept;};
using OrderedQuotientRowsHandle=std::unique_ptr<OrderedQuotientRows,OrderedQuotientRowsDeleter>;
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE) \
    && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST) \
    && !defined(__AVX__) && !defined(__AVX2__) && !defined(__AVX512F__) && !defined(__FMA__) \
    && !defined(EIGEN_SCALAR_BINARY_OP_PLUGIN) \
    && EIGEN_WORLD_VERSION==3 && EIGEN_MAJOR_VERSION==4 && EIGEN_MINOR_VERSION==1
// The original TU stays Precise/SSE2. Only the new complete worker helper
// observes status and requires separately rounded product/add operations.
#pragma float_control(precise,on,push)
#pragma float_control(except,on)
#pragma fenv_access(on)
class OrderedQuotientRows
{
    struct FloatingState {
        std::fenv_t native{};unsigned csr=0;bool valid=false;
        bool capture() noexcept {csr=_mm_getcsr();valid=std::fegetenv(&native)==0;return valid;}
        void restore() const noexcept {
            if(!valid)return;
            if(std::fesetenv(&native)!=0)std::terminate();_mm_setcsr(csr);
        }
        bool enterFresh(const FloatingState& caller) const noexcept {
            auto clean=caller.native;clean._Fe_stat=0;
            if(std::fesetenv(&clean)!=0)return false;
            _mm_setcsr(caller.csr&~0x3fu);return true;
        }
    };
    struct Worker {
        FloatingState original;
        unsigned long raisedNative=0;unsigned raisedCsr=0;
        std::array<double,static_cast<std::size_t>(cgsTileRows)> rowSums;
    };
    Eigen::Index size_=0;
    std::size_t entries_=0,workerCapacity_=0,bytes_=0;
    std::unique_ptr<Eigen::Index[]> outer_;
    std::unique_ptr<std::uint32_t[]> columns_;
    std::unique_ptr<double[]> values_;
    // Caller holds the existing factorMutex for the entire GMRES solve.
    mutable std::unique_ptr<Worker[]> workers_;
    mutable std::array<double,static_cast<std::size_t>(maximumCgsColumns)> projections_{};
    static bool sameControl(const FloatingState& a,const FloatingState& b) noexcept {
        auto x=a.native,y=b.native;x._Fe_stat=0;y._Fe_stat=0;
        return std::memcmp(&x,&y,sizeof x)==0&&(a.csr&~0x3fu)==(b.csr&~0x3fu);
    }
    static bool normalOrZero(const double& value) noexcept {
        std::uint64_t bits=0;std::memcpy(&bits,&value,sizeof bits);
        const auto exponent=bits&UINT64_C(0x7ff0000000000000);
        return exponent!=UINT64_C(0x7ff0000000000000)
            &&(exponent!=0||(bits&UINT64_C(0x000fffffffffffff))==0);
    }
    static bool qualified(FloatingState& caller) noexcept {
        if(!caller.capture()||omp_in_parallel()||omp_get_max_threads()<=1)return false;
        unsigned control=0;
        if(_controlfp_s(&control,0,0)!=0||(control&_MCW_EM)!=_MCW_EM
            ||(control&_MCW_RC)!=_RC_NEAR||(caller.csr&0x1f80u)!=0x1f80u
            ||(caller.csr&0x6000u)!=0||(caller.csr&((1u<<15)|(1u<<6)))!=0)return false;
#if defined(_EM_AMBIGUOUS)
        if(control&_EM_AMBIGUOUS)return false;
#endif
        return sizeof(double)==8&&std::numeric_limits<double>::digits==53;
    }
public:
    static OrderedQuotientRowsHandle prepare(const Sparse& matrix) noexcept {
        OrderedSpmvOptionalRestore restore;FloatingState caller;
        if(!qualified(caller))return {};
        try {
            static_assert(!Sparse::IsRowMajor,"ordered quotient source must be CSC");
            if(!matrix.isCompressed()||matrix.rows()<=0||matrix.rows()!=matrix.cols()
                ||static_cast<std::uint64_t>(matrix.rows())>std::numeric_limits<std::uint32_t>::max()
                ||matrix.nonZeros()<0)return {};
            const auto n=static_cast<std::size_t>(matrix.rows()),e=static_cast<std::size_t>(matrix.nonZeros());
            const auto team=std::min(n,static_cast<std::size_t>(omp_get_max_threads()));
            if(team<2||team>static_cast<std::size_t>(std::numeric_limits<int>::max()))return {};
            constexpr std::size_t limit=64u*1024u*1024u;
            std::size_t stock=sizeof(OrderedQuotientRows)+512u;
            const auto add=[&](std::size_t& sum,std::size_t count,std::size_t width) {
                if(sum>limit||count>(limit-sum)/width)return false;sum+=count*width;return true;
            };
            if(!add(stock,n+1,sizeof(Eigen::Index))||!add(stock,e,sizeof(std::uint32_t))
                ||!add(stock,e,sizeof(double))||!add(stock,team,sizeof(Worker)))return {};
            std::size_t peak=stock;
            if(!add(peak,n,sizeof(Eigen::Index)))return {}; // temporary construction cursor also counts
            OrderedQuotientRowsHandle result(new OrderedQuotientRows);
            result->size_=matrix.rows();result->entries_=e;result->workerCapacity_=team;result->bytes_=stock;
            result->outer_=std::make_unique<Eigen::Index[]>(n+1);
            result->columns_=std::unique_ptr<std::uint32_t[]>(new std::uint32_t[e]);
            result->values_=std::unique_ptr<double[]>(new double[e]);
            result->workers_=std::make_unique<Worker[]>(team);
            auto cursor=std::unique_ptr<Eigen::Index[]>(new Eigen::Index[n]);
            // Integer-only counting; no coefficient arithmetic or zero pruning.
            for(Eigen::Index col=0;col<matrix.outerSize();++col)
                for(Sparse::InnerIterator entry(matrix,col);entry;++entry) {
                    if(entry.row()<0||entry.row()>=matrix.rows()||!normalOrZero(entry.value()))return {};
                    ++result->outer_[static_cast<std::size_t>(entry.row())+1];
                }
            for(std::size_t row=0;row<n;++row){
                result->outer_[row+1]+=result->outer_[row];cursor[row]=result->outer_[row];
            }
            if(result->outer_[n]!=matrix.nonZeros())return {};
            // Stable append in EXACT original CSC col/inner order. Thus every
            // row sees its original sequence, including explicit signed zeros.
            for(Eigen::Index col=0;col<matrix.outerSize();++col)
                for(Sparse::InnerIterator entry(matrix,col);entry;++entry) {
                    const auto slot=static_cast<std::size_t>(cursor[static_cast<std::size_t>(entry.row())]++);
                    result->columns_[slot]=static_cast<std::uint32_t>(col);
                    result->values_[slot]=entry.value();
                }
            for(std::size_t row=0;row<n;++row)if(cursor[row]!=result->outer_[row+1])return {};
            return result;
        }catch(...){return {};}// Optional layout only; original factor/numerical work is outside this scope.
    }
    std::size_t bytes() const noexcept {return bytes_;}
    bool multiply(const double* rhs,Eigen::Index count,double* output,int& actualWorkers) const noexcept {
        actualWorkers=0;const int callerErrno=errno;FloatingState caller;
        const bool capable=qualified(caller);errno=callerErrno;
        if(!capable){caller.restore();errno=callerErrno;return false;}
        if(count!=size_||!rhs||!output){caller.restore();errno=callerErrno;return false;}
        for(Eigen::Index i=0;i<count;++i)if(!normalOrZero(rhs[i])){
            caller.restore();errno=callerErrno;return false;
        }
        const auto requested=std::min(workerCapacity_,static_cast<std::size_t>(omp_get_max_threads()));
        if(requested<2){caller.restore();errno=callerErrno;return false;}
        const int team=static_cast<int>(requested);int setupFailure=0,actualTeam=0;
        for(std::size_t t=0;t<requested;++t){workers_[t].raisedNative=0;workers_[t].raisedCsr=0;}
#pragma omp parallel num_threads(team)
        {
            const int tid=omp_get_thread_num(),seen=omp_get_num_threads();
            if(tid==0)actualTeam=seen;
            auto& local=workers_[static_cast<std::size_t>(tid)];const int workerErrno=errno;
            const bool captured=local.original.capture();
            const bool ready=captured&&local.original.enterFresh(caller);
            if(!ready){
#pragma omp atomic write
                setupFailure=1;
            }
            if(ready)errno=callerErrno;
#pragma omp for schedule(static)
            for(Eigen::Index row=0;row<count;++row){
                if(!ready)continue;
                double sum=0.; // Original Eigen assignment evalTo first setZero(), not first-term assignment.
                for(Eigen::Index slot=outer_[static_cast<std::size_t>(row)];slot<outer_[static_cast<std::size_t>(row)+1];++slot){
                    // Actual Eigen3.4.1: rhs_j(alpha*rhs[col]); then res[row]+=value*rhs_j.
                    // normal/zero RHS makes repeating alpha==1 exact and flag-neutral.
                    volatile double y=1.0*rhs[columns_[static_cast<std::size_t>(slot)]];
                    volatile double product=values_[static_cast<std::size_t>(slot)]*y;
                    sum+=product;
                }
                output[row]=sum;
            }
            if(ready){
                const int arithmeticErrno=errno;FloatingState after;
                if(!after.capture()||arithmeticErrno!=callerErrno||!sameControl(after,caller)){
#pragma omp atomic write
                    setupFailure=1;
                }else{local.raisedNative=after.native._Fe_stat;local.raisedCsr=after.csr&0x3fu;}
            }
            local.original.restore();errno=workerErrno;
        }
        if(setupFailure||actualTeam<2){caller.restore();errno=callerErrno;return false;}
        unsigned long raisedNative=0;unsigned raisedCsr=0;
        for(int t=0;t<actualTeam;++t){raisedNative|=workers_[static_cast<std::size_t>(t)].raisedNative;raisedCsr|=workers_[static_cast<std::size_t>(t)].raisedCsr;}
        auto merged=caller.native;merged._Fe_stat|=raisedNative;
        if(std::fesetenv(&merged)!=0)std::terminate();_mm_setcsr(caller.csr|raisedCsr);errno=callerErrno;
        actualWorkers=actualTeam;return true;
    }
    Cgs2Execution orthogonalize(const Eigen::MatrixXd& basis,Eigen::VectorXd& work,
        Eigen::MatrixXd& h,Eigen::Index column,int& actualWorkers) const noexcept
    {
        actualWorkers=0;const int callerErrno=errno;FloatingState caller;
        const bool capable=qualified(caller);errno=callerErrno;
        if(!capable||work.size()!=size_||basis.rows()!=size_||column<0
            ||column>=maximumCgsColumns||column>=basis.cols()||column>=h.cols()
            ||column>=h.rows()){
            caller.restore();errno=callerErrno;return Cgs2Execution::NotStarted;
        }
        const auto requested=std::min(workerCapacity_,static_cast<std::size_t>(omp_get_max_threads()));
        if(requested<2){caller.restore();errno=callerErrno;return Cgs2Execution::NotStarted;}
        const int team=static_cast<int>(requested);int setupFailure=0,controlFailure=0,actualTeam=0;
        std::atomic<int> arithmeticFailure{0};bool started=false;
        const auto tiles=work.size()>0?(work.size()-1)/cgsTileRows+1:0;
        for(std::size_t t=0;t<requested;++t){workers_[t].raisedNative=0;workers_[t].raisedCsr=0;}
#pragma omp parallel num_threads(team)
        {
            const int tid=omp_get_thread_num();if(tid==0)actualTeam=omp_get_num_threads();
            auto& local=workers_[static_cast<std::size_t>(tid)];const int workerErrno=errno;
            const bool captured=local.original.capture(),ready=captured&&local.original.enterFresh(caller);
            if(!ready){
#pragma omp atomic write
                setupFailure=1;
            }
            if(ready)errno=callerErrno;
#pragma omp barrier
            const bool execute=setupFailure==0;
            if(tid==0)started=execute;
            for(int pass=0;pass<2;++pass){
#pragma omp for schedule(static)
                for(Eigen::Index k=0;k<=column;++k){
                    if(execute&&arithmeticFailure.load(std::memory_order_relaxed)==0
                        &&!cgs2Projection(basis,work,h,column,k,projections_.data()))
                        arithmeticFailure.store(1,std::memory_order_relaxed);
                }
#pragma omp for schedule(static)
                for(Eigen::Index tile=0;tile<tiles;++tile)
                    if(execute&&arithmeticFailure.load(std::memory_order_relaxed)==0)
                        cgs2UpdateTile(basis,work,column,tile*cgsTileRows,projections_.data(),local.rowSums.data());
            }
            if(execute){
                const int arithmeticErrno=errno;FloatingState after;
                const bool observed=after.capture();
                if(observed){local.raisedNative=after.native._Fe_stat;local.raisedCsr=after.csr&0x3fu;}
                if(!observed||arithmeticErrno!=callerErrno||!sameControl(after,caller)){
#pragma omp atomic write
                    controlFailure=1;
                }
            }
            local.original.restore();errno=workerErrno;
        }
        if(!started){caller.restore();errno=callerErrno;return Cgs2Execution::NotStarted;}
        // Once either pass has written work/H there is NEVER serial replay.
        unsigned long raisedNative=0;unsigned raisedCsr=0;
        for(int t=0;t<actualTeam;++t){raisedNative|=workers_[static_cast<std::size_t>(t)].raisedNative;raisedCsr|=workers_[static_cast<std::size_t>(t)].raisedCsr;}
        auto merged=caller.native;merged._Fe_stat|=raisedNative;
        if(std::fesetenv(&merged)!=0)std::terminate();_mm_setcsr(caller.csr|raisedCsr);errno=callerErrno;
        if(controlFailure||arithmeticFailure.load(std::memory_order_relaxed)!=0||!work.allFinite())
            return Cgs2Execution::FailedAfterStart;
        actualWorkers=actualTeam;return Cgs2Execution::Completed;
    }
};
#if !defined(_M_FP_STRICT) || !_M_FP_STRICT
#pragma fenv_access(off)
#endif
#pragma float_control(pop)
#else
class OrderedQuotientRows
{
public:
    static OrderedQuotientRowsHandle prepare(const Sparse&) noexcept {return {};}
    std::size_t bytes() const noexcept {return 0;}
    bool multiply(const double*,Eigen::Index,double*,int&) const noexcept {return false;}
    Cgs2Execution orthogonalize(const Eigen::MatrixXd&,Eigen::VectorXd&,Eigen::MatrixXd&,
        Eigen::Index,int&) const noexcept {return Cgs2Execution::NotStarted;}
};
#endif
void OrderedQuotientRowsDeleter::operator()(OrderedQuotientRows* value) const noexcept
{
    if(!value)return;OrderedSpmvOptionalRestore restore;delete value;
}

// A right preconditioner for the ORIGINAL quotient. Its diagonal/Schur
// approximation never replaces a coefficient in the Krylov matrix-vector.
// Private, version-audited adapter for Eigen3.4.1 only. Public analyzePattern
// allocates L before returning; the count-only dry run below precedes that
// allocation. Numeric LDLT and solve remain the original Eigen implementation.
constexpr std::size_t globalPressureByteLimit=128u*1024u*1024u;
class BoundedPressureLdlt : public Eigen::SimplicialLDLT<Sparse>
{
    using Ldlt=Eigen::SimplicialLDLT<Sparse>;
    using I=typename Ldlt::StorageIndex;
public:
    bool analyzeBounded(const Sparse& s,std::size_t& admission,std::size_t& fill){
#if EIGEN_WORLD_VERSION==3 && EIGEN_MAJOR_VERSION==4 && EIGEN_MINOR_VERSION==1
        if(s.rows()<=0||s.rows()!=s.cols()||s.nonZeros()<0||this->m_matrix.rows()!=0
            ||this->m_matrix.data().allocatedSize()!=0)return false;
        const auto n=static_cast<std::uint64_t>(s.rows()),nz=static_cast<std::uint64_t>(s.nonZeros());
        const auto indexLimit=static_cast<std::uint64_t>(std::numeric_limits<I>::max());
        if(n>indexLimit||nz>indexLimit/2||n+1>indexLimit/8)return false;
        // AMD uses c+c/5+2*n storage indices, with c<=2*nnz(S).
        const auto c=2*nz;
        if(c>indexLimit-c/5||2*n>indexLimit-c-c/5)return false;
        std::size_t bytes=128u*1024u;
        const auto claim=[&](std::uint64_t count,std::size_t width){
            if(count>(globalPressureByteLimit-bytes)/width)return false;
            bytes+=static_cast<std::size_t>(count)*width;return true;};
        // Declared conservative input/order/scratch payload model, not a
        // vendor allocator/RSS bound. Strict-L payload below is exact.
        if(!claim(nz,8*(sizeof(double)+sizeof(I)))||!claim(n+1,64*sizeof(I))
            ||!claim(n,16*sizeof(double)))return false;
        const auto limit=(globalPressureByteLimit-bytes)/(sizeof(double)+sizeof(I));
        typename Ldlt::CholMatrixType ap;
        typename Ldlt::ConstCholMatrixPtr pmat=nullptr;
        this->ordering(s,pmat,ap);
        if(!pmat||pmat->rows()!=s.rows()||pmat->cols()!=s.cols())return false;
        const auto size=static_cast<I>(n);
        std::vector<I> parent(static_cast<std::size_t>(n)),counts(parent.size()),tags(parent.size());
        std::size_t nnz=0;
        // Count-only adaptation of Eigen SimplicialCholesky_impl.h:35-57
        // (MPL2.0; Gael Guennebaud/Eigen, adapted from Timothy A.Davis LDL).
        // It reads the SAME preordered upper triangle; no numeric work or
        // private flags are copied. Every visit is one stored strict-L entry.
        for(I k=0;k<size;++k){parent[k]=-1;tags[k]=k;counts[k]=0;
            for(typename Ldlt::CholMatrixType::InnerIterator it(*pmat,k);it;++it){I i=it.index();
                if(i<0||i>=size)return false;
                if(i<k)for(;tags[i]!=k;i=parent[i]){
                    if(parent[i]==-1)parent[i]=k;
                    if(nnz>=limit||nnz>=indexLimit||counts[i]>=std::numeric_limits<I>::max())return false;
                    ++nnz;++counts[i];tags[i]=k;
                }}
        }
        // Original Eigen analysis repeats the same integer count, then does
        // its first strict-L fill allocation, now bounded by this certificate.
        this->analyzePattern_preordered(*pmat,true);
        if(this->info()!=Eigen::Success||this->m_matrix.nonZeros()<0
            ||static_cast<std::size_t>(this->m_matrix.nonZeros())!=nnz)return false;
        admission=bytes+nnz*(sizeof(double)+sizeof(I));fill=nnz;return true;
#else
        (void)s;(void)admission;(void)fill;return false;
#endif
    }
    std::size_t storageBytes() const{
        const auto entries=static_cast<std::size_t>(this->m_matrix.data().allocatedSize());
        const auto indices=static_cast<std::size_t>(this->m_matrix.outerSize()+1)+this->m_parent.size()
            +this->m_nonZerosPerCol.size()+this->m_P.indices().size()+this->m_Pinv.indices().size();
        return entries*(sizeof(double)+sizeof(I))+indices*sizeof(I)+static_cast<std::size_t>(this->m_diag.size())*sizeof(double);
    }
};
struct GlobalPressureFactorReceipt
{
    bool available=false,prepared=false,timingsAvailable=false;
    std::size_t rows=0,nonZeros=0,factorNonZeros=0,storageBytes=0,admissionBytes=0,byteLimit=globalPressureByteLimit;
    double symbolicSeconds=0,numericSeconds=0;
};
class BlockQuotientFactor
{
    const Sparse* matrix_=nullptr;
    OrderedQuotientRowsHandle orderedRows_;
    bool orderedSpmvRequested_=false;
    mutable DirectSparse directInput_;
    bool directInputDeferred_=false;
    mutable std::atomic<bool> directInputMaterialized_{false};
    Sparse coupling_;
    Eigen::VectorXd inverseDiagonal_,rowScale_,inverseRowScale_;
    BoundedPressureLdlt schur_;
    GlobalPressureFactorReceipt pressureGlobal_;
    LampLinear::ParallelBlockJacobi parallelSchur_;
    Detail::HydroPressureCoarseCorrection pressureCoarse_;
    Eigen::Index nf_=0,nq_=0;
    bool eligible_=false,gmresOnlyParallel_=false,gmresGlobalPressure_=false;
    double preconditionerSeconds_=0,preparedDirectSeconds_=0;
    mutable std::unique_ptr<Eigen::SuperLU<DirectSparse>> direct_;
    bool precondition(const Eigen::Ref<const Eigen::VectorXd>& rhs,Eigen::Ref<Eigen::VectorXd> out,QuotientSolveTrace& trace) const
    {
        const Eigen::VectorXd t=inverseDiagonal_.array()*rhs.head(nf_).array();
        if(!t.allFinite())return false;
        if(nq_==0){out=t;return out.allFinite();}
        const Eigen::VectorXd pressureRhs=coupling_.transpose()*t-rhs.tail(nq_);
        if(!pressureRhs.allFinite())return false;
        Eigen::VectorXd p(nq_);
        if(gmresOnlyParallel_){
            if(gmresGlobalPressure_){
                PassiveAssemblySampler sampler;const auto started=sampler.mark();
                p=schur_.solve(pressureRhs);const bool solved=schur_.info()==Eigen::Success;
                sampler.addSeconds(started,trace.pressureGlobalSolveSeconds);
                trace.pressureGlobalTimingsAvailable=trace.pressureGlobalTimingsAvailable&&sampler.available;
                ++trace.pressureGlobalSolveCalls;if(!solved)return false;
            }else if(pressureCoarse_.ready()){
                const auto local=[this](const Eigen::Ref<const Eigen::VectorXd>& in,Eigen::Ref<Eigen::VectorXd> value){
                    return parallelSchur_.solve(in,value);};
                if(!pressureCoarse_.apply(pressureRhs,p,local))return false;
                ++trace.pressureCoarseApplications;
            }else if(!parallelSchur_.solve(pressureRhs,p))return false;
        }
        else{p=schur_.solve(pressureRhs);if(schur_.info()!=Eigen::Success)return false;}
        if(!p.allFinite())return false;
        const Eigen::VectorXd cp=coupling_*p;
        out.head(nf_)=t-(inverseDiagonal_.array()*cp.array()).matrix();
        out.tail(nq_)=p;
        return out.allFinite();
    }
    bool gmres(const Eigen::VectorXd& rhs,Eigen::VectorXd& x,QuotientSolveTrace& trace,const Eigen::VectorXd* initialGuess=nullptr) const
    {
        ++trace.gmresSolves;PassiveAssemblySampler timings;timings.requireMaskedExceptions();
        trace.solveBreakdownTimingsAvailable=trace.solveBreakdownTimingsAvailable&&timings.available;
        const Eigen::Index n=rhs.size(),restart=std::min<Eigen::Index>(gmresOnlyParallel_?maximumCgsColumns:50,n);
        // More bounded numerical work for the explicit GMRES-only route;
        // the original true residual target and physical gates remain unchanged.
        const int originalBudget=300,budget=gmresOnlyParallel_?2*originalBudget:originalBudget;
        const double target=1e-12; // Inner only; original local gates and compensated outer refinements decide acceptance.
        trace.lastSolveIterations=0;trace.lastSolveBudget=budget;trace.lastSolveRestartLength=static_cast<int>(restart);
        trace.lastSolveExit="arithmetic_or_preconditioner_failure";
        const double rhsNorm=rhs.stableNorm();
        if(!std::isfinite(rhsNorm))return false;
        x=Eigen::VectorXd::Zero(n);
        if(rhsNorm==0){trace.lastRelativeResidual=0;trace.lastScaledRelativeResidual=0;const bool zero=rhs.isZero(0);
            if(zero)trace.lastSolveExit="converged";return zero;}
        const Eigen::VectorXd scaledRhs=(rowScale_.array()*rhs.array()).matrix();
        const double scaledRhsNorm=scaledRhs.stableNorm();
        if(!scaledRhs.allFinite()||!(scaledRhsNorm>0)||!std::isfinite(scaledRhsNorm))return false;
        Eigen::VectorXd residual=rhs,scaledResidual=scaledRhs;
        if(initialGuess){
            if(initialGuess->size()!=n||!initialGuess->allFinite())return false;
            x=*initialGuess;
            const auto trueResidualStarted=timings.mark();residual=rhs-(*matrix_)*x;
            scaledResidual=(rowScale_.array()*residual.array()).matrix();
            if(!residual.allFinite()||!scaledResidual.allFinite())return false;
            timings.addSeconds(trueResidualStarted,trace.krylovTrueResidualSeconds);
            trace.solveBreakdownTimingsAvailable=trace.solveBreakdownTimingsAvailable&&timings.available;
            trace.lastRelativeResidual=residual.stableNorm()/rhsNorm;
            trace.lastScaledRelativeResidual=scaledResidual.stableNorm()/scaledRhsNorm;
            if(!std::isfinite(trace.lastRelativeResidual)||!std::isfinite(trace.lastScaledRelativeResidual))return false;
            if(trace.lastScaledRelativeResidual<=target){trace.lastSolveExit="converged";return true;}
        }
        Eigen::MatrixXd v(n,restart+1),z(n,restart),h(restart+1,restart);
        Eigen::VectorXd cs(restart),sn(restart),g(restart+1),w(n),physicalBasis(n);
        int total=0,cycle=0,noProgress=0;
        while(total<budget){
            const double beta=scaledResidual.stableNorm();
            trace.lastRelativeResidual=residual.stableNorm()/rhsNorm;
            trace.lastScaledRelativeResidual=beta/scaledRhsNorm;
            if(!std::isfinite(beta)||!std::isfinite(trace.lastRelativeResidual))return false;
            if(trace.lastScaledRelativeResidual<=target){const bool valid=x.allFinite();
                if(valid)trace.lastSolveExit="converged";return valid;}
            if(cycle++>0)++trace.restarts;
            v.col(0)=scaledResidual/beta;h.setZero();g.setZero();g[0]=beta;
            Eigen::Index used=0;
            // Retain the original partial cycle and true-residual assessment
            // at300 before spending any additional GMRES-only work.
            const int cycleBudget=total<originalBudget?originalBudget:budget;
            for(Eigen::Index j=0;j<restart&&total<cycleBudget;++j){
                // Arnoldi acts on R Q P^-1 R^-1. z is ALWAYS a physical
                // correction, so x and every original outer gate are unscaled.
                physicalBasis=(inverseRowScale_.array()*v.col(j).array()).matrix();
                if(!physicalBasis.allFinite())return false;
                const auto preconditionStarted=timings.mark();const bool conditioned=precondition(physicalBasis,z.col(j),trace);
                timings.addSeconds(preconditionStarted,trace.krylovPreconditionApplySeconds);
                trace.solveBreakdownTimingsAvailable=trace.solveBreakdownTimingsAvailable&&timings.available;
                if(!conditioned)return false;
                const auto spmvStarted=timings.mark();++trace.orderedSpmvCalls;int actualWorkers=0;
                const auto column=z.col(j);
                const bool ordered=orderedRows_&&column.innerStride()==1
                    &&orderedRows_->multiply(column.data(),n,w.data(),actualWorkers);
                if(ordered){++trace.orderedSpmvParallelCalls;
                    trace.orderedSpmvMaximumWorkers=std::max(trace.orderedSpmvMaximumWorkers,actualWorkers);}
                else{trace.orderedSpmvFallback=true;w.noalias()=(*matrix_)*z.col(j);}
                w.array()*=rowScale_.array();const bool finiteProduct=w.allFinite();
                timings.addSeconds(spmvStarted,trace.krylovSpmvSeconds);
                trace.solveBreakdownTimingsAvailable=trace.solveBreakdownTimingsAvailable&&timings.available;
                if(!finiteProduct)return false;
                const auto orthogonalizationStarted=timings.mark();int orthogonalizationWorkers=0;
                const auto execution=orderedRows_?orderedRows_->orthogonalize(v,w,h,j,orthogonalizationWorkers)
                    :Cgs2Execution::NotStarted;
                const bool orthogonalized=execution==Cgs2Execution::Completed
                    ||(execution==Cgs2Execution::NotStarted&&serialCgs2(v,w,h,j));
                timings.addSeconds(orthogonalizationStarted,trace.krylovOrthogonalizationSeconds);
                trace.solveBreakdownTimingsAvailable=trace.solveBreakdownTimingsAvailable&&timings.available;
                if(!orthogonalized)return false;
                if(execution==Cgs2Execution::NotStarted)orthogonalizationWorkers=1;
                if(orthogonalizationWorkers>1)++trace.orthogonalizationParallelCalls;
                trace.orthogonalizationMaximumWorkers=std::max(trace.orthogonalizationMaximumWorkers,orthogonalizationWorkers);
                const double next=w.stableNorm();
                if(!std::isfinite(next))return false;
                h(j+1,j)=next;
                if(next>0)v.col(j+1)=w/next;
                for(Eigen::Index i=0;i<j;++i){
                    const double upper=cs[i]*h(i,j)+sn[i]*h(i+1,j);
                    h(i+1,j)=-sn[i]*h(i,j)+cs[i]*h(i+1,j);h(i,j)=upper;
                }
                const double norm=std::hypot(h(j,j),h(j+1,j));
                if(!(norm>0)||!std::isfinite(norm))return false;
                cs[j]=h(j,j)/norm;sn[j]=h(j+1,j)/norm;
                h(j,j)=norm;h(j+1,j)=0;
                g[j+1]=-sn[j]*g[j];g[j]=cs[j]*g[j];
                ++total;++trace.iterations;trace.lastSolveIterations=total;used=j+1;
                if(!std::isfinite(g[j+1]))return false;
                if(next==0||std::abs(g[j+1])/scaledRhsNorm<=target)break;
            }
            if(used==0)return false;
            const Eigen::VectorXd y=h.topLeftCorner(used,used).template triangularView<Eigen::Upper>().solve(g.head(used));
            if(!y.allFinite())return false;
            x.noalias()+=z.leftCols(used)*y;
            if(!x.allFinite())return false;
            const auto trueResidualStarted=timings.mark();residual=rhs-(*matrix_)*x;
            scaledResidual=(rowScale_.array()*residual.array()).matrix();
            if(!residual.allFinite()||!scaledResidual.allFinite())return false;
            const double nextBeta=scaledResidual.stableNorm();
            trace.lastRelativeResidual=residual.stableNorm()/rhsNorm;
            trace.lastScaledRelativeResidual=nextBeta/scaledRhsNorm;
            timings.addSeconds(trueResidualStarted,trace.krylovTrueResidualSeconds);
            trace.solveBreakdownTimingsAvailable=trace.solveBreakdownTimingsAvailable&&timings.available;
            if(!std::isfinite(trace.lastRelativeResidual)||!std::isfinite(trace.lastScaledRelativeResidual))return false;
            if(trace.lastScaledRelativeResidual<=target){trace.lastSolveExit="converged";return true;}
            // A numerical-work budget only; stagnation NEVER means success.
            noProgress=nextBeta>=beta?noProgress+1:0;
            if(noProgress>=3){trace.stagnated=true;trace.lastSolveExit="stagnation";return false;}
        }
        trace.lastSolveExit="iteration_budget";return false;
    }
public:
    void compute(const Sparse& matrix,Eigen::Index nf,bool gmresOnlyParallel,bool gmresGlobalPressure)
    {
        orderedRows_.reset();orderedSpmvRequested_=false;
        matrix_=&matrix;nf_=nf;nq_=matrix.rows()-nf;gmresOnlyParallel_=gmresOnlyParallel;
        gmresGlobalPressure_=gmresGlobalPressure;
        // All actual input/ABI gates still precede preconditioner arithmetic.
        // Only the certified unused private allocation is deferred.
        directInputDeferred_=!gmresOnlyParallel_&&deferredDirectInputCertified(matrix);
        if(!gmresOnlyParallel_&&!directInputDeferred_){
            directInput_=directQuotient(matrix);
            directInputMaterialized_.store(true,std::memory_order_release);
        }
        const auto started=WallClock::now();
        eligible_=matrix.rows()==matrix.cols()&&nf_>0&&nq_>=0;
        inverseDiagonal_.resize(nf_);
        rowScale_.resize(matrix.rows());inverseRowScale_.resize(matrix.rows());
        if(eligible_)for(Eigen::Index i=0;i<nf_;++i){
            const double d=matrix.coeff(i,i);
            if(!(d>0)||!std::isfinite(d)){eligible_=false;break;}
            inverseDiagonal_[i]=1/d;
            if(!(inverseDiagonal_[i]>0)||!std::isfinite(inverseDiagonal_[i])){eligible_=false;break;}
            rowScale_[i]=std::sqrt(inverseDiagonal_[i]);inverseRowScale_[i]=1/rowScale_[i];
            if(!(rowScale_[i]>0)||!std::isfinite(rowScale_[i])||!(inverseRowScale_[i]>0)||!std::isfinite(inverseRowScale_[i])){eligible_=false;break;}
        }
        if(eligible_&&nq_>0){
            coupling_=matrix.block(0,nf_,nf_,nq_);coupling_.makeCompressed();
            Sparse weighted=coupling_;
            for(Eigen::Index col=0;col<weighted.outerSize();++col)
                for(Sparse::InnerIterator e(weighted,col);e;++e)e.valueRef()*=inverseDiagonal_[e.row()];
            Sparse s=coupling_.transpose()*weighted;s.makeCompressed();
            for(Eigen::Index i=0;i<s.nonZeros();++i)if(!std::isfinite(s.valuePtr()[i]))eligible_=false;
            if(eligible_)for(Eigen::Index i=0;i<nq_;++i){
                const double diagonal=s.coeff(i,i);
                if(!(diagonal>0)||!std::isfinite(diagonal)){eligible_=false;break;}
                rowScale_[nf_+i]=1/std::sqrt(diagonal);inverseRowScale_[nf_+i]=1/rowScale_[nf_+i];
                if(!(rowScale_[nf_+i]>0)||!std::isfinite(rowScale_[nf_+i])||!(inverseRowScale_[nf_+i]>0)||!std::isfinite(inverseRowScale_[nf_+i])){eligible_=false;break;}
            }
            if(eligible_){
                if(gmresOnlyParallel_&&gmresGlobalPressure_){
                    pressureGlobal_.rows=static_cast<std::size_t>(nq_);pressureGlobal_.nonZeros=static_cast<std::size_t>(s.nonZeros());
                    PassiveAssemblySampler sampler;const auto symbolic=sampler.mark();
                    eligible_=schur_.analyzeBounded(s,pressureGlobal_.admissionBytes,pressureGlobal_.factorNonZeros);
                    pressureGlobal_.symbolicSeconds=sampler.seconds(symbolic);
                    if(eligible_){const auto numeric=sampler.mark();schur_.factorize(s);
                        pressureGlobal_.numericSeconds=sampler.seconds(numeric);eligible_=schur_.info()==Eigen::Success;
                        if(eligible_){const Eigen::VectorXd pivots=schur_.vectorD();
                            eligible_=pivots.size()==nq_&&pivots.allFinite()&&(pivots.array()>0).all();}
                        if(eligible_){pressureGlobal_.storageBytes=schur_.storageBytes();
                            eligible_=pressureGlobal_.storageBytes<=pressureGlobal_.admissionBytes
                                &&pressureGlobal_.admissionBytes<=pressureGlobal_.byteLimit;}}
                    pressureGlobal_.prepared=eligible_;pressureGlobal_.available=eligible_;
                    pressureGlobal_.timingsAvailable=eligible_&&sampler.available;
                }else if(gmresOnlyParallel_){
                    eligible_=s.rows()<=std::numeric_limits<int>::max()&&s.nonZeros()<=std::numeric_limits<int>::max();
                    if(eligible_){const Eigen::SparseMatrix<double> localInput=s;
                        eligible_=parallelSchur_.compute(localInput,256,true);
                        // Optional policy is FIXED before the first Krylov call.
                        // Unsupported coarse preparation retains original RAS;
                        // an active coarse apply failure never selects RAS.
                        if(eligible_)(void)pressureCoarse_.compute(localInput);}
                }else{
                    schur_.compute(s);
                    eligible_=schur_.info()==Eigen::Success;
                    if(eligible_){const Eigen::VectorXd pivots=schur_.vectorD();
                        eligible_=pivots.size()==nq_&&pivots.allFinite()&&(pivots.array()>0).all();}
                }
            }
        }
        preconditionerSeconds_=elapsedSeconds(started);
        if(!eligible_&&!gmresOnlyParallel_){QuotientSolveTrace trace;forceDirect(trace,1);preparedDirectSeconds_=trace.lazyFactorSeconds;}
        orderedSpmvRequested_=eligible_&&matrix.rows()>0;
        if(orderedSpmvRequested_)orderedRows_=OrderedQuotientRows::prepare(matrix);
    }
    Eigen::ComputationInfo info() const
    {return direct_?direct_->info():(eligible_?Eigen::Success:Eigen::InvalidInput);}
    bool usingDirect() const noexcept {return bool(direct_);}
    bool directInputDeferred() const noexcept {return directInputDeferred_;}
    bool directInputMaterialized() const noexcept {return directInputMaterialized_.load(std::memory_order_acquire);}
    double preconditionerSeconds() const noexcept {return preconditionerSeconds_;}
    double preparedDirectSeconds() const noexcept {return preparedDirectSeconds_;}
    bool orderedSpmvRequested() const noexcept {return orderedSpmvRequested_;}
    const GlobalPressureFactorReceipt& pressureGlobalReport() const noexcept {return pressureGlobal_;}
    const Detail::HydroPressureCoarseCorrection::Report& pressureCoarseReport() const noexcept {return pressureCoarse_.report();}
    bool pressureCoarsePrepared() const noexcept {return pressureCoarse_.ready();}
    int parallelPressureBlockCount() const noexcept {return gmresOnlyParallel_&&!gmresGlobalPressure_&&nq_>0?parallelSchur_.blockCount():0;}
    int parallelPressureMaximumWorkers() const noexcept {return gmresOnlyParallel_&&!gmresGlobalPressure_&&nq_>0?parallelSchur_.maximumWorkers():0;}
    bool orderedSpmvEligible() const noexcept {return bool(orderedRows_);}
    std::size_t orderedSpmvViewBytes() const noexcept {return orderedRows_?orderedRows_->bytes():0;}
    // Caller serializes all mutable native workspace and lazy construction.
    void forceDirect(QuotientSolveTrace& trace,int reason) const
    {
        if(gmresOnlyParallel_)throw std::runtime_error("[HydroSinglePhaseSaddle] GMRES-only candidate failed; global direct fallback is disabled");
        trace.directOnly=true;trace.fallbackReason=reason;
        if(direct_)return;
        const auto started=WallClock::now();
        if(!directInputMaterialized_.load(std::memory_order_acquire)){
            // Exact same current private quotient; full original conversion and
            // checks are repeated immediately before actual native factor use.
            directInput_=directQuotient(*matrix_);
            directInputMaterialized_.store(true,std::memory_order_release);
        }
        direct_=std::make_unique<Eigen::SuperLU<DirectSparse>>();
        direct_->compute(directInput_);
        if(direct_->info()!=Eigen::Success)
            throw std::runtime_error("[HydroSinglePhaseSaddle] quotient saddle factorization failed");
        trace.lazyFactorSeconds+=elapsedSeconds(started);
    }
    Eigen::VectorXd solve(const Eigen::VectorXd& rhs,QuotientSolveTrace& trace,const Eigen::VectorXd* initialGuess=nullptr) const
    {
        if(gmresOnlyParallel_){
            if(!eligible_)throw std::runtime_error("[HydroSinglePhaseSaddle] parallel pressure preconditioner preparation failed");
            Eigen::VectorXd x;
            if(gmres(rhs,x,trace,initialGuess))return x;
            std::ostringstream message;message.precision(17);message.setf(std::ios::scientific,std::ios::floatfield);
            message<<"[HydroSinglePhaseSaddle] bounded block GMRES failed: iterations="<<trace.iterations
                <<", restarts="<<trace.restarts<<", physical_residual="<<trace.lastRelativeResidual
                <<", scaled_residual="<<trace.lastScaledRelativeResidual<<", stagnated="<<trace.stagnated
                <<", pressure_blocks="<<parallelPressureBlockCount()
                <<", inner_solves="<<trace.gmresSolves<<", last_solve_iterations="<<trace.lastSolveIterations
                <<", last_solve_budget="<<trace.lastSolveBudget<<", restart_length="<<trace.lastSolveRestartLength
                <<", termination="<<trace.lastSolveExit
                <<", coarse_failure_stage="<<pressureCoarse_.failureStage()<<", coarse_failure_row="<<pressureCoarse_.failureRow()
                <<", coarse_failure_defect="<<pressureCoarse_.failureDefect()<<", coarse_failure_allowed="<<pressureCoarse_.failureAllowed();
            throw std::runtime_error(message.str());
        }
        if(!trace.directOnly){
            if(eligible_){Eigen::VectorXd x;if(gmres(rhs,x,trace))return x;forceDirect(trace,2);}
            else forceDirect(trace,1);
        }
        ++trace.directSolves;
        const Eigen::VectorXd x=direct_->solve(rhs);
        require(direct_->info()==Eigen::Success,"SuperLU quotient saddle backsolve failed");
        return x;
    }
};
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
#pragma float_control(precise,on,push)
#pragma float_control(except,on)
#pragma fenv_access(on)
#endif
// Optional ordered-row SPD scans. They use the existing team/state owner;
// unsupported arithmetic or storage leaves the original literal path intact.
int spdRowTeam(Eigen::Index count,bool stickyInexact) noexcept
{
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE) && defined(FE_TONEAREST) && defined(FE_INEXACT) && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST)
    if(count<512||!std::numeric_limits<double>::is_iec559||sizeof(double)!=sizeof(std::uint64_t)
        ||std::numeric_limits<double>::radix!=2||std::numeric_limits<double>::digits!=53
        ||std::fegetround()!=FE_TONEAREST)return 1;
    const unsigned csr=_mm_getcsr();unsigned int control=0;std::fenv_t native{};
    if(std::fegetenv(&native)!=0||_controlfp_s(&control,0,0)!=0
        ||(control&_MCW_EM)!=_MCW_EM||(control&_MCW_RC)!=_RC_NEAR
        ||(csr&0x1f80u)!=0x1f80u||(csr&0xe040u)!=0)return 1;
#ifdef _EM_AMBIGUOUS
    if((control&_EM_AMBIGUOUS)!=0)return 1;
#endif
    // The literal dominance scan can set inexact. Suppressing worker flags is
    // invisible only when BOTH original status registers already contain it.
    if(stickyInexact&&((native._Fe_stat&FE_INEXACT)==0||(csr&0x20u)==0))return 1;
    return LampLinear::BlockJacobiDetail::workersFor(count);
#else
    (void)count;(void)stickyInexact;return 1;
#endif
}
bool spdNormalOrZero(double value) noexcept
{
    std::uint64_t word=0;std::memcpy(&word,&value,sizeof(word));word&=UINT64_C(0x7fffffffffffffff);
    return word==0||(word>=UINT64_C(0x0010000000000000)&&word<UINT64_C(0x7ff0000000000000));
}
std::optional<bool> gatheredLowerPositiveDominance(const HydroSaddleVelocityBlock::Matrix& matrix)
{
    static_assert(HydroSaddleVelocityBlock::Matrix::IsRowMajor,"ordered lower slots require RowMajor velocity storage");
    const LampLinear::BlockJacobiDetail::FloatingState restore;
    if(!restore.enter())return std::nullopt;
    const int team=spdRowTeam(matrix.rows(),true);if(team<=1)return std::nullopt;
    try {
        const auto count=static_cast<std::size_t>(matrix.rows());
        const auto nonzeros=static_cast<std::size_t>(matrix.nonZeros());
        constexpr std::size_t payloadLimit=64u*1024u*1024u,overhead=512u;
        constexpr std::size_t rowBytes=2*sizeof(double)+2*sizeof(std::uint32_t);
        if(count>UINT32_MAX||nonzeros>UINT32_MAX||count>(payloadLimit-overhead)/rowBytes)return std::nullopt;
        const auto rowPayload=count*rowBytes;
        if(nonzeros>(payloadLimit-overhead-rowPayload)/sizeof(std::uint32_t))return std::nullopt;
        // Integer positions into immutable RowMajor storage: no matrix/value
        // copy, no coefficient change, and no global Q/E-like dense object.
        std::vector<std::uint32_t> offsets(count+1,0),cursor(count,0),slots;
        const auto* outer=matrix.outerIndexPtr();const auto* inner=matrix.innerIndexPtr();
        const auto* values=matrix.valuePtr();
        if(outer[0]!=0||outer[count]!=matrix.nonZeros())return std::nullopt;
        for(std::size_t row=0;row<count;++row){
            if(outer[row]<0||outer[row+1]<outer[row]||outer[row+1]>matrix.nonZeros())return std::nullopt;
            for(Eigen::Index p=outer[row];p<outer[row+1];++p){
                const auto col=inner[p];if(col<0||col>=matrix.cols())return std::nullopt;
                if(col<static_cast<Eigen::Index>(row))++offsets[static_cast<std::size_t>(col)+1];
            }
        }
        for(std::size_t row=0;row<count;++row){offsets[row+1]+=offsets[row];cursor[row]=offsets[row];}
        slots.resize(offsets[count]);
        for(std::size_t row=0;row<count;++row)
            for(Eigen::Index p=outer[row];p<outer[row+1];++p)
                if(inner[p]<static_cast<Eigen::Index>(row))
                    slots[cursor[static_cast<std::size_t>(inner[p])]++]=static_cast<std::uint32_t>(p);
        std::vector<double> diagonal(count,0.),radii(count,0.);int observed=1;
        const bool complete=LampLinear::BlockJacobiDetail::parallelJobs(matrix.rows(),team,[&](Eigen::Index row){
            double diagonalValue=0,radius=0;bool hasDiagonal=false;
            const auto add=[&](double value){
                if(!spdNormalOrZero(value))return false;
                const double magnitude=std::abs(value);
                if(magnitude==0)return true;
                volatile double sum=radius+magnitude;
                if(!spdNormalOrZero(sum)||!(sum>0))return false;
                const double upper=std::nextafter(sum,std::numeric_limits<double>::infinity());
                if(!spdNormalOrZero(upper)||!(upper>0))return false;
                radius=upper;return true;
            };
            for(Eigen::Index p=outer[row];p<outer[row+1];++p){
                const auto col=inner[p];const double value=values[p];
                // Upper entries were also validated by the literal scan.
                if(!spdNormalOrZero(value))return false;
                if(row==col){if(hasDiagonal||!(value>0))return false;hasDiagonal=true;diagonalValue=value;}
                else if(row>col&&!add(value))return false;
            }
            // Original outer rows i>row visit their col==row contributions in
            // ascending i, after this row's own lower entries. Continue the
            // SAME accumulator; separately summing either segment is invalid.
            for(std::uint32_t p=offsets[static_cast<std::size_t>(row)];
                p<offsets[static_cast<std::size_t>(row)+1];++p)
                if(!add(values[slots[p]]))return false;
            if(!hasDiagonal)return false;
            diagonal[static_cast<std::size_t>(row)]=diagonalValue;radii[static_cast<std::size_t>(row)]=radius;
            return true;
        },observed);
        if(!complete)return std::nullopt;
        // All rows/operands completed normally, even when an early row is
        // not dominant. A mathematical false does not require literal replay.
        for(std::size_t row=0;row<count;++row)if(!(diagonal[row]>radii[row]))return false;
        return true;
    }catch(...){return std::nullopt;}
}
template<class Positive,class MultiplyLower>
bool parallelSpdMassScales(const std::vector<double>& masses,double dt,
    std::vector<double>& diagonal,std::vector<double>& scales,double& minimum,
    const Positive& positive,const MultiplyLower& multiplyLower)
{
    const LampLinear::BlockJacobiDetail::FloatingState restore;
    if(!restore.enter())return false;
    const auto count=static_cast<Eigen::Index>(masses.size());const int team=spdRowTeam(count,false);
    if(team<=1||masses.size()>(64u*1024u*1024u-512u)/(4*sizeof(double)))return false;
    try {
        std::vector<double> lowerValues(masses.size(),0.);int observed=1;
        if(!LampLinear::BlockJacobiDetail::parallelJobs(count,team,[&](Eigen::Index i){
            const auto row=static_cast<std::size_t>(i);
            if(!positive(masses[row]))return false;
            volatile double d=masses[row]/dt;if(!positive(d))return false;diagonal[row]=d;
            volatile double root=std::sqrt(d);if(!positive(root))return false;
            volatile double s=1./root;if(!positive(s))return false;scales[row]=s;
            double first=0,lower=0;if(!multiplyLower(s,d,first)||!multiplyLower(first,s,lower))return false;
            lowerValues[row]=lower;return true;
        },observed))return false;
        double result=std::numeric_limits<double>::infinity();
        for(const double lower:lowerValues)result=std::min(result,lower);
        minimum=result;return true;
    }catch(...){return false;}
}
bool parallelSpdEnergyExtrema(const HydroSaddleVelocityBlock::Matrix& matrix,
    const HydroEnergyViscosity::Matrix& viscosity,const std::vector<double>& diagonal,
    std::size_t& entries,double& maximumS,double& maximumK,double& maximumT,double& maximumR)
{
    static_assert(HydroSaddleVelocityBlock::Matrix::IsRowMajor,"energy row scan requires RowMajor velocity storage");
    static_assert(HydroEnergyViscosity::Matrix::IsRowMajor,"energy row scan requires RowMajor viscosity storage");
    const LampLinear::BlockJacobiDetail::FloatingState restore;
    if(!restore.enter())return false;
    const auto count=matrix.rows();const int team=spdRowTeam(count,false);if(team<=1)return false;
    struct Row {std::size_t entries=0;double s=0,k=0,t=0,r=0;};
    constexpr std::size_t exactLimit=4503599627370494ULL;
    if(static_cast<std::size_t>(count)>(64u*1024u*1024u-512u)/(sizeof(Row)+3*sizeof(double)))return false;
    try {
        std::vector<Row> rows(static_cast<std::size_t>(count));int observed=1;
        if(!LampLinear::BlockJacobiDetail::parallelJobs(count,team,[&](Eigen::Index i){
            auto& row=rows[static_cast<std::size_t>(i)];
            const auto visit=[&](Eigen::Index j,double sv,double kv){
                if(!spdNormalOrZero(sv)||!spdNormalOrZero(kv))return false;
                volatile double half=.5*kv;
                if(!spdNormalOrZero(half)||(kv!=0&&half==0))return false;
                volatile double t=sv-half;
                volatile double r=t-(i==j?diagonal[static_cast<std::size_t>(i)]:0.);
                if(!spdNormalOrZero(t)||!spdNormalOrZero(r))return false;
                const std::size_t multiplicity=i==j?1u:2u;
                if(row.entries>exactLimit-multiplicity)return false;row.entries+=multiplicity;
                row.s=std::max(row.s,std::abs(sv));row.k=std::max(row.k,std::abs(kv));
                row.t=std::max(row.t,std::abs(static_cast<double>(t)));
                row.r=std::max(row.r,i==j?std::max(0.,-static_cast<double>(r)):std::abs(static_cast<double>(r)));return true;
            };
            HydroSaddleVelocityBlock::Matrix::InnerIterator s(matrix,i);
            HydroEnergyViscosity::Matrix::InnerIterator v(viscosity,i);bool seenDiagonal=false;
            while((s&&s.col()<=i)||(v&&v.col()<=i)){
                const Eigen::Index sc=s&&s.col()<=i?s.col():count;
                const Eigen::Index kc=v&&v.col()<=i?v.col():count,j=std::min(sc,kc);
                if(!visit(j,sc==j?s.value():0.,kc==j?v.value():0.))return false;
                if(i==j)seenDiagonal=true;if(sc==j)++s;if(kc==j)++v;
            }
            if(!seenDiagonal&&!visit(i,0.,0.))return false;
            return true;
        },observed))return false;
        std::size_t combinedEntries=0;double s=0,k=0,t=0,r=0;
        for(const auto& row:rows){
            if(combinedEntries>exactLimit-row.entries)return false;combinedEntries+=row.entries;
            s=std::max(s,row.s);k=std::max(k,row.k);t=std::max(t,row.t);r=std::max(r,row.r);
        }
        entries=combinedEntries;maximumS=s;maximumK=k;maximumT=t;maximumR=r;return true;
    }catch(...){return false;}
}
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
#if !defined(_M_FP_STRICT) || !_M_FP_STRICT
#pragma fenv_access(off)
#endif
#pragma float_control(pop)
#endif

// A sufficient SPD certificate for exactly the Lower selfadjoint completion
// consumed by the original SimplicialLLT, without modifying any stored entry.
// Each strict-lower |a_ij| contributes to BOTH symmetric row radii. With
// certified upper bounds r_i, d_i>r_i implies
// x^T S x >= sum_i (d_i-r_i)*x_i^2 > 0 for every nonzero real x.
// This does not assert exact symmetry of the caller's approximately symmetric A.
bool strictLowerPositiveDominance(const HydroSaddleVelocityBlock::Matrix& matrix)
{
#if !defined(FE_TONEAREST) || defined(__FAST_MATH__) || defined(_M_FP_FAST)
    return false;
#else
    if(!std::numeric_limits<double>::is_iec559||sizeof(double)!=sizeof(std::uint64_t)
        ||std::numeric_limits<double>::radix!=2||std::numeric_limits<double>::digits!=53
        ||std::fegetround()!=FE_TONEAREST)return false;
    if(!matrix.isCompressed()||matrix.rows()!=matrix.cols()||matrix.rows()<=0)return false;
    if(const auto gathered=gatheredLowerPositiveDominance(matrix))return *gathered;
    const auto size=static_cast<std::size_t>(matrix.rows());
    std::vector<double> diagonal(size,0.),radius(size,0.);
    std::vector<unsigned char> hasDiagonal(size,0);
    const double infinity=std::numeric_limits<double>::infinity();
    const auto magnitudeBits=[](double value){
        std::uint64_t word=0;std::memcpy(&word,&value,sizeof(word));
        return word&UINT64_C(0x7fffffffffffffff);
    };
    const auto subnormal=[&](double value){const auto word=magnitudeBits(value);
        return word!=0&&word<UINT64_C(0x0010000000000000);};
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE) && defined(FE_INEXACT)
    // Only pre-existing sticky inexact permits an invisible successful probe:
    // both this bound and the old normal-only certificate can then add no new
    // status bit. Restore the complete native/raw state on EVERY probe result.
    const int savedErrno=errno;std::fenv_t caller{};const unsigned callerCsr=_mm_getcsr();
    const bool captured=std::fegetenv(&caller)==0;unsigned int control=0;
    bool fastQualified=captured&&_controlfp_s(&control,0,0)==0
        &&(control&_MCW_EM)==_MCW_EM&&(control&_MCW_RC)==_RC_NEAR
        &&(callerCsr&0x1f80u)==0x1f80u&&(callerCsr&0x6000u)==0
        &&(caller._Fe_stat&FE_INEXACT)!=0&&(callerCsr&0x20u)!=0;
#ifdef _EM_AMBIGUOUS
    fastQualified=fastQualified&&(control&_EM_AMBIGUOUS)==0;
#endif
    if(fastQualified){
        bool certified=false;
        try{
            const auto fastCertificate=[&](){
                std::vector<std::uint64_t> terms(size,0);
                constexpr std::uint64_t maximumTerms=UINT64_C(1)<<49;
                const auto addRounded=[&](Eigen::Index row,double magnitude){
                    if(magnitudeBits(magnitude)==0)return true;
                    const auto index=static_cast<std::size_t>(row);
                    if(terms[index]+1>=maximumTerms)return false;
                    volatile double sum=radius[index]+magnitude;
                    if(!std::isfinite(sum)||subnormal(sum))return false;
                    radius[index]=sum;++terms[index];return true;
                };
                for(Eigen::Index outer=0;outer<matrix.outerSize();++outer)
                    for(HydroSaddleVelocityBlock::Matrix::InnerIterator entry(matrix,outer);entry;++entry){
                        const auto row=entry.row(),col=entry.col();const double value=entry.value();
                        if(!std::isfinite(value)||subnormal(value))return false;
                        if(row==col){
                            if(hasDiagonal[static_cast<std::size_t>(row)]||!(value>0))return false;
                            hasDiagonal[static_cast<std::size_t>(row)]=1;diagonal[static_cast<std::size_t>(row)]=value;
                        }else if(row>col){const double magnitude=std::abs(value);
                            if(!addRounded(row,magnitude)||!addRounded(col,magnitude))return false;
                        }
                    }
                for(std::size_t row=0;row<size;++row){
                    if(!hasDiagonal[row])return false;
                    if(terms[row]==0){if(!(diagonal[row]>0))return false;continue;}
                    // n<2^49: 3*n, (3*n)*2^-52 and 1-(3*n)*2^-52
                    // are exact binary64. This bounds the ORIGINAL repeated
                    // RN-plus-successor radius, not merely the real sum.
                    volatile double penalty=static_cast<double>(3*terms[row])*std::numeric_limits<double>::epsilon();
                    volatile double denominator=1.-penalty;
                    if(!(denominator>0.625))return false;
                    volatile double quotient=radius[row]/denominator;
                    if(!std::isfinite(quotient)||subnormal(quotient))return false;
                    const double upper=std::nextafter(quotient,infinity);
                    if(!std::isfinite(upper)||subnormal(upper)||!(diagonal[row]>upper))return false;
                }
                return true;
            };
            certified=fastCertificate();
        }catch(...){certified=false;} // Uncertain work takes the original path.
        if(std::fesetenv(&caller)!=0)std::terminate();_mm_setcsr(callerCsr);errno=savedErrno;
        if(certified)return true;
        std::fill(diagonal.begin(),diagonal.end(),0.);std::fill(radius.begin(),radius.end(),0.);
        std::fill(hasDiagonal.begin(),hasDiagonal.end(),0);
    }else errno=savedErrno;
#endif
    // Reject all nonzero subnormal operands/results by bits, never through
    // a DAZ-sensitive floating zero comparison. All accepted positive sums
    // then contain normal operands only, so FTZ/DAZ cannot shrink a radius.
    const auto addUpper=[&](Eigen::Index row,double magnitude){
        // For nonnegative binary64 operands in qualified RNE arithmetic,
        // one successor of the rounded sum encloses their exact real sum.
        // A zero term is exact and must not inflate an empty row radius.
        if(magnitudeBits(magnitude)==0)return true;
        volatile double sum=radius[static_cast<std::size_t>(row)]+magnitude;
        if(!std::isfinite(sum)||subnormal(sum))return false;
        const double upper=std::nextafter(sum,infinity);
        if(!std::isfinite(upper)||subnormal(upper))return false;
        radius[static_cast<std::size_t>(row)]=upper;return true;
    };
    for(Eigen::Index outer=0;outer<matrix.outerSize();++outer)
        for(HydroSaddleVelocityBlock::Matrix::InnerIterator entry(matrix,outer);entry;++entry){
            const auto row=entry.row(),col=entry.col();const double value=entry.value();
            if(!std::isfinite(value)||subnormal(value))return false;
            if(row==col){
                if(hasDiagonal[static_cast<std::size_t>(row)]||!(value>0))return false;
                hasDiagonal[static_cast<std::size_t>(row)]=1;diagonal[static_cast<std::size_t>(row)]=value;
            }else if(row>col){
                const double magnitude=std::abs(value);
                if(!addUpper(row,magnitude)||!addUpper(col,magnitude))return false;
            }
        }
    for(std::size_t row=0;row<size;++row)
        if(!hasDiagonal[row]||!(diagonal[row]>radius[row]))return false;
    return true;
#endif
}

struct EnergySpdInputs
{
    const HydroMomentumMass& mass;
    const HydroEnergyViscosity& viscosity;
    double dt;
};
// A sufficient certificate for exactly the Lower selfadjoint completion used
// by the original LLT. This optional probe cannot change its caller's FP state.
bool energyLowerPositiveCertificate(const HydroSaddleVelocityBlock::Matrix& matrix,
    const EnergySpdInputs& input,double& boundUpper,double& massLower,bool& boundsAvailable)
{
#if defined(__FAST_MATH__) || (defined(_M_FP_FAST) && _M_FP_FAST)
    return false;
#elif defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE) && defined(FE_TONEAREST)
    struct FloatingRestore {
        std::fenv_t native{};unsigned csr=_mm_getcsr();int error=errno;
        bool captured=std::fegetenv(&native)==0;
        ~FloatingRestore() noexcept {
            if(captured&&std::fesetenv(&native)!=0)std::terminate();
            _mm_setcsr(csr);errno=error;
        }
    } restore;
    if(!restore.captured||!std::numeric_limits<double>::is_iec559
        ||sizeof(double)!=sizeof(std::uint64_t)||std::numeric_limits<double>::digits!=53
        ||std::numeric_limits<double>::radix!=2||std::fegetround()!=FE_TONEAREST)return false;
    unsigned int control=0;
    if(_controlfp_s(&control,0,0)!=0||(control&_MCW_EM)!=_MCW_EM||(control&_MCW_RC)!=_RC_NEAR
        ||(restore.csr&0x1f80u)!=0x1f80u||(restore.csr&0xe040u)!=0)return false;
#ifdef _EM_AMBIGUOUS
    if((control&_EM_AMBIGUOUS)!=0)return false;
#endif
    try {
        const auto normalOrZero=[](double value){
            std::uint64_t bits=0;std::memcpy(&bits,&value,sizeof(bits));bits&=UINT64_C(0x7fffffffffffffff);
            return bits==0||(bits>=UINT64_C(0x0010000000000000)&&bits<UINT64_C(0x7ff0000000000000));
        };
        const auto positive=[&](double value){return normalOrZero(value)&&value>0;};
        const double infinity=std::numeric_limits<double>::infinity();
        // Every nontrivial operation is a separate volatile binary64 RN step,
        // followed by its outward successor; no multiply/add is contracted.
        const auto addUpper=[&](double a,double b,double& out){
            if(a==0){out=b;return normalOrZero(out)&&out>=0;}
            if(b==0){out=a;return normalOrZero(out)&&out>=0;}
            volatile double value=a+b;if(!positive(value))return false;
            out=std::nextafter(value,infinity);return positive(out);
        };
        const auto multiplyUpper=[&](double a,double b,double& out){
            if(a==0||b==0){out=0;return true;}
            if(!positive(a)||!positive(b))return false;
            volatile double value=a*b;if(!positive(value))return false;
            out=std::nextafter(value,infinity);return positive(out);
        };
        const auto multiplyLower=[&](double a,double b,double& out){
            if(!positive(a)||!positive(b))return false;
            volatile double value=a*b;if(!positive(value))return false;
            out=std::nextafter(value,0.);return positive(out);
        };
        const auto subtract=[&](double a,double b,bool lower,double& out){
            if(!normalOrZero(a)||!normalOrZero(b))return false;
            if(a==b){out=0;return true;} // Exact binary64 endpoint subtraction.
            if(b==0){out=a;return true;}if(a==0){out=-b;return true;}
            volatile double value=a-b;if(!normalOrZero(value)||value==0)return false;
            out=std::nextafter(value,lower?-infinity:infinity);return normalOrZero(out);
        };
        const auto& k=input.viscosity.matrix();const auto& masses=input.mass.masses();
        const Eigen::Index count=matrix.rows();
        if(count<=0||matrix.cols()!=count||k.rows()!=count||k.cols()!=count
            ||!matrix.isCompressed()||!k.isCompressed()||masses.size()!=static_cast<std::size_t>(count)
            ||!positive(input.dt))return false;
        std::vector<double> diagonal(count),scales(count),radii(count,0.);
        double minimum=std::numeric_limits<double>::infinity();
        if(!parallelSpdMassScales(masses,input.dt,diagonal,scales,minimum,positive,multiplyLower)){
        for(Eigen::Index i=0;i<count;++i){
            if(!positive(masses[static_cast<std::size_t>(i)]))return false;
            volatile double d=masses[static_cast<std::size_t>(i)]/input.dt;if(!positive(d))return false;
            diagonal[static_cast<std::size_t>(i)]=d;
            volatile double root=std::sqrt(d);if(!positive(root))return false;
            volatile double s=1./root;if(!positive(s))return false;scales[static_cast<std::size_t>(i)]=s;
            double first=0,lower=0;if(!multiplyLower(s,d,first)||!multiplyLower(first,s,lower))return false;
            minimum=std::min(minimum,lower);
        }
        }
        // Producer checks BOTH original K-build FP controls and this query's
        // controls, and bounds symmetric Lower(K)-exact(E^T weights E).
        // Coarse sufficient bound: cheap exact-support maxima first. The old
        // tighter interval path remains the fallback for every uncertainty.
        const auto coarseCertificate=[&](){
            const double scale=*std::max_element(scales.begin(),scales.end());
            const auto gramCoarse=input.viscosity.coarseScaledGramLowerErrorBound(scale);
            if(!gramCoarse||!normalOrZero(*gramCoarse)||*gramCoarse<0)return false;
            constexpr std::size_t exactLimit=4503599627370494ULL;
            std::size_t entries=0;double maximumS=0,maximumK=0,maximumT=0,maximumR=0;
            const double maximumD=*std::max_element(diagonal.begin(),diagonal.end());
            const auto visit=[&](Eigen::Index i,Eigen::Index j,double sv,double kv){
                if(!normalOrZero(sv)||!normalOrZero(kv))return false;
                volatile double half=.5*kv;
                // Normal binary scaling by one half is exact. Subnormal
                // uncertainty is deliberately left to the original path.
                if(!normalOrZero(half)||(kv!=0&&half==0))return false;
                volatile double t=sv-half;
                volatile double r=t-(i==j?diagonal[static_cast<std::size_t>(i)]:0.);
                if(!normalOrZero(t)||!normalOrZero(r))return false;
                const std::size_t multiplicity=i==j?1u:2u;
                if(entries>exactLimit-multiplicity)return false;entries+=multiplicity;
                maximumS=std::max(maximumS,std::abs(sv));maximumK=std::max(maximumK,std::abs(kv));
                maximumT=std::max(maximumT,std::abs(static_cast<double>(t)));
                maximumR=std::max(maximumR,i==j?std::max(0.,-static_cast<double>(r)):std::abs(static_cast<double>(r)));return true;
            };
            if(!parallelSpdEnergyExtrema(matrix,k,diagonal,entries,maximumS,maximumK,maximumT,maximumR)){
            for(Eigen::Index i=0;i<count;++i){
                HydroSaddleVelocityBlock::Matrix::InnerIterator s(matrix,i);
                HydroEnergyViscosity::Matrix::InnerIterator v(k,i);bool seenDiagonal=false;
                while((s&&s.col()<=i)||(v&&v.col()<=i)){
                    const Eigen::Index sc=s&&s.col()<=i?s.col():count;
                    const Eigen::Index kc=v&&v.col()<=i?v.col():count,j=std::min(sc,kc);
                    if(!visit(i,j,sc==j?s.value():0.,kc==j?v.value():0.))return false;
                    if(i==j)seenDiagonal=true;if(sc==j)++s;if(kc==j)++v;
                }
                if(!seenDiagonal&&!visit(i,i,0.,0.))return false;
            }
            }
            double halfK=0,terms=0,next=0,error=0,residual=0,scaleSquared=0;
            if(!multiplyUpper(.5,maximumK,halfK)||!addUpper(maximumS,halfK,terms)
                ||!addUpper(terms,maximumT,next)||!addUpper(next,maximumD,terms)
                ||!multiplyUpper(std::numeric_limits<double>::epsilon(),terms,error)
                ||!addUpper(error,2.*std::numeric_limits<double>::denorm_min(),next)
                ||!addUpper(maximumR,next,residual)||!multiplyUpper(scale,scale,scaleSquared))return false;
            // R differs from Rbad only by a nonnegative diagonal. The same
            // entrywise error still encloses min(R_ii,0) (a 1-Lipschitz map).
            // Count BOTH symmetric offdiagonal entries: ||D Rbad D||2 is at
            // most its Frobenius bound smax^2*sqrt(entries)*entryBound.
            // Existing exactLimit makes this count conversion exact. On this
            // qualified MSVC x64 RNE binary64 path sqrt is correctly rounded;
            // store before taking its outward successor, never a fast rsqrt.
            volatile double supportRoot=std::sqrt(static_cast<double>(entries));
            if(!positive(supportRoot))return false;
            const double supportUpper=std::nextafter(supportRoot,infinity);
            if(!positive(supportUpper))return false;
            double allResidual=0,scaledResidual=0,halfGram=0,total=0;
            if(!multiplyUpper(supportUpper,residual,allResidual)
                ||!multiplyUpper(scaleSquared,allResidual,scaledResidual)
                ||!multiplyUpper(.5,*gramCoarse,halfGram)||!addUpper(scaledResidual,halfGram,total))return false;
            if(!(total<minimum))return false;
            boundUpper=total;massLower=minimum;boundsAvailable=true;return true;
        };
        if(coarseCertificate())return true;
        const auto gram=input.viscosity.scaledGramLowerErrorBounds(scales);
        if(!gram||gram->size()!=static_cast<std::size_t>(count))return false;
        const auto accumulate=[&](Eigen::Index i,Eigen::Index j,double sv,double kv){
            if(!normalOrZero(sv)||!normalOrZero(kv))return false;
            double halfLow=0,halfHigh=0;
            if(kv!=0){
                volatile double half=.5*kv;if(!normalOrZero(half)||half==0)return false;
                halfLow=std::nextafter(half,-infinity);halfHigh=std::nextafter(half,infinity);
                if(!normalOrZero(halfLow)||!normalOrZero(halfHigh))return false;
            }
            double low=0,high=0;
            if(!subtract(sv,halfHigh,true,low)||!subtract(sv,halfLow,false,high))return false;
            if(i==j){double nextLow=0,nextHigh=0;
                if(!subtract(low,diagonal[static_cast<std::size_t>(i)],true,nextLow)
                    ||!subtract(high,diagonal[static_cast<std::size_t>(i)],false,nextHigh))return false;
                low=nextLow;high=nextHigh;
            }
            double first=0,magnitude=0;const auto ui=static_cast<std::size_t>(i),uj=static_cast<std::size_t>(j);
            if(!multiplyUpper(std::max(std::abs(low),std::abs(high)),scales[ui],first)
                ||!multiplyUpper(first,scales[uj],magnitude))return false;
            double sum=0;if(!addUpper(radii[ui],magnitude,sum))return false;radii[ui]=sum;
            if(i!=j){if(!addUpper(radii[uj],magnitude,sum))return false;radii[uj]=sum;}
            return true;
        };
        // Merge BOTH complete Lower supports; missing S/K entries are exact
        // zero. A mass diagonal is processed even if absent from both matrices.
        for(Eigen::Index i=0;i<count;++i){
            HydroSaddleVelocityBlock::Matrix::InnerIterator s(matrix,i);
            HydroEnergyViscosity::Matrix::InnerIterator v(k,i);bool seenDiagonal=false;
            while((s&&s.col()<=i)||(v&&v.col()<=i)){
                const Eigen::Index sc=s&&s.col()<=i?s.col():count;
                const Eigen::Index kc=v&&v.col()<=i?v.col():count;const Eigen::Index j=std::min(sc,kc);
                const double sv=sc==j?s.value():0.,kv=kc==j?v.value():0.;
                if(!accumulate(i,j,sv,kv))return false;if(i==j)seenDiagonal=true;
                if(sc==j)++s;if(kc==j)++v;
            }
            if(!seenDiagonal&&!accumulate(i,i,0.,0.))return false;
        }
        double maximum=0;
        for(Eigen::Index i=0;i<count;++i){
            const double error=(*gram)[static_cast<std::size_t>(i)];if(!normalOrZero(error)||error<0)return false;
            double half=0,total=0;if(!multiplyUpper(.5,error,half)||!addUpper(radii[static_cast<std::size_t>(i)],half,total))return false;
            maximum=std::max(maximum,total);
        }
        boundUpper=maximum;massLower=minimum;boundsAvailable=true;
        return maximum<minimum;
    }catch(...){return false;} // Optional probe uncertainty keeps original LLT.
#else
    (void)matrix;(void)input;(void)boundUpper;(void)massLower;(void)boundsAvailable;return false;
#endif
}

double ratio(double error, double scale)
{
    if (!std::isfinite(error) || !std::isfinite(scale) || scale<0)
        throw std::runtime_error("[HydroSinglePhaseSaddle] nonfinite residual/scale");
    if (scale==0) {
        if (error!=0) throw std::runtime_error("[HydroSinglePhaseSaddle] nonzero residual at zero scale");
        return 0;
    }
    const double result=std::abs(error)/scale;
    if (!std::isfinite(result)) throw std::runtime_error("[HydroSinglePhaseSaddle] nonfinite normalized residual");
    return result;
}
void accept(double defect, double tolerance, const char* stage)
{
    if (defect>tolerance) throw std::runtime_error(std::string("[HydroSinglePhaseSaddle] ")+stage+" residual exceeds tolerance");
}
struct CompensatedSum
{
    double high=0,low=0;
    void add(double value)
    {
        if(!std::isfinite(value))throw std::runtime_error("[HydroSinglePhaseSaddle] nonfinite quotient component term");
        const double next=high+value;
        low+=std::abs(high)>=std::abs(value)?(high-next)+value:(value-next)+high;
        high=next;
        if(!std::isfinite(high)||!std::isfinite(low))
            throw std::runtime_error("[HydroSinglePhaseSaddle] nonfinite quotient component sum");
    }
    double value() const
    {
        const double result=high+low;
        if(!std::isfinite(result))throw std::runtime_error("[HydroSinglePhaseSaddle] nonfinite quotient component sum");
        return result;
    }
};
int root(std::vector<int>& parents, int p)
{
    while (parents[p]!=p) {parents[p]=parents[parents[p]];p=parents[p];}
    return p;
}
}
namespace Detail
{
using ValidationMatrix=HydroSaddleVelocityBlock::Matrix;
// Only private execution probes are restored. Original matrix/norm/gate calls
// remain outside this scope and retain their own flags and exception order.
class SparseValidationFloatingRestore
{
    std::fenv_t native_{};int error_=errno;
#if defined(_M_X64) || defined(__SSE2__)
    unsigned csr_=_mm_getcsr();
#endif
public:
    bool captured=std::fegetenv(&native_)==0;
    ~SparseValidationFloatingRestore() noexcept
    {
        if(captured&&std::fesetenv(&native_)!=0)std::terminate();
#if defined(_M_X64) || defined(__SSE2__)
        _mm_setcsr(csr_);
#endif
        errno=error_;
    }
};
bool sparseValidationFloatingQualified() noexcept
{
    SparseValidationFloatingRestore restore;
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE) && !defined(__FAST_MATH__) && !defined(_M_FP_FAST) && defined(FE_TONEAREST)
    if(!restore.captured||!std::numeric_limits<double>::is_iec559||sizeof(double)!=sizeof(std::uint64_t)
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

// One current numerical stage, not a retained coefficient or execution cache.
// Bodies are allocation-free/noexcept. Setup failure restores the entry state
// and the caller executes the unchanged whole stage in its original order.
struct MatrixStageExecution
{
    bool parallelUsed=false,parallelFallback=false;
    int maximumWorkers=0;
};
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE) \
    && !defined(__FAST_MATH__) && !defined(_M_FP_FAST) \
    && !defined(__AVX__) && !defined(__AVX2__) && !defined(__AVX512F__) && !defined(__FMA__)
#pragma float_control(precise,on,push)
#pragma float_control(except,on)
#pragma fenv_access(on)
template<class Body> bool tryMatrixStage(Eigen::Index count,Body& body,
    std::size_t availableBytes,int& actualWorkers,bool& attempted) noexcept
{
    static_assert(std::is_nothrow_invocable_v<Body&,Eigen::Index>,"current matrix body must not throw");
    actualWorkers=0;attempted=false;
    const int callerErrno=errno;
    struct State {
        std::fenv_t native{};unsigned csr=0;bool valid=false;
        bool capture() noexcept {csr=_mm_getcsr();valid=std::fegetenv(&native)==0;return valid;}
        void restore() const noexcept {if(valid&&std::fesetenv(&native)!=0)std::terminate();_mm_setcsr(csr);}
        bool enter(const State& caller) const noexcept {
            auto clean=caller.native;clean._Fe_stat=0;
            if(std::fesetenv(&clean)!=0)return false;_mm_setcsr(caller.csr&~0x3fu);return true;
        }
    } caller;
    if(count<=1||omp_in_parallel()||omp_get_max_threads()<=1||!sparseValidationFloatingQualified())return false;
    if(!caller.capture()){errno=callerErrno;return false;}
    attempted=true;
    struct Worker {State entry;unsigned long nativeFlags=0;unsigned csrFlags=0;};
    const auto team=static_cast<int>(std::min<Eigen::Index>(count,omp_get_max_threads()));
    if(availableBytes<sizeof(std::vector<Worker>)
        ||static_cast<std::size_t>(team)>(availableBytes-sizeof(std::vector<Worker>))/sizeof(Worker)){
        caller.restore();errno=callerErrno;return false;
    }
    try {
        std::vector<Worker> workers(static_cast<std::size_t>(team));
        if(workers.capacity()>(availableBytes-sizeof(workers))/sizeof(Worker)){
            std::vector<Worker>().swap(workers);
            caller.restore();errno=callerErrno;return false;
        }
        int failed=0,actualTeam=0;
#pragma omp parallel num_threads(team)
        {
            const int tid=omp_get_thread_num();auto& local=workers[static_cast<std::size_t>(tid)];
            if(tid==0)actualTeam=omp_get_num_threads();
            const int entryErrno=errno;const bool captured=local.entry.capture();
            const bool ready=captured&&local.entry.enter(caller);
            if(!ready){
#pragma omp atomic write
                failed=1;
            }
            if(ready)errno=callerErrno;
#pragma omp for schedule(static)
            for(Eigen::Index k=0;k<count;++k)if(ready)body(k);
            if(ready){
                State after;const int afterErrno=errno;
                bool same=after.capture();
                if(same){auto x=after.native,y=caller.native;x._Fe_stat=0;y._Fe_stat=0;
                    same=std::memcmp(&x,&y,sizeof x)==0&&(after.csr&~0x3fu)==(caller.csr&~0x3fu);}
                if(!same||afterErrno!=callerErrno){
#pragma omp atomic write
                    failed=1;
                }else{local.nativeFlags=after.native._Fe_stat;local.csrFlags=after.csr&0x3fu;}
            }
            local.entry.restore();errno=entryErrno;
        }
        if(failed||actualTeam<2){std::vector<Worker>().swap(workers);caller.restore();errno=callerErrno;return false;}
        unsigned long nativeFlags=0;unsigned csrFlags=0;
        for(int t=0;t<actualTeam;++t){nativeFlags|=workers[static_cast<std::size_t>(t)].nativeFlags;csrFlags|=workers[static_cast<std::size_t>(t)].csrFlags;}
        // Release optional storage before publishing the physical stage flags.
        std::vector<Worker>().swap(workers);
        auto merged=caller.native;merged._Fe_stat|=nativeFlags;
        if(std::fesetenv(&merged)!=0)std::terminate();_mm_setcsr(caller.csr|csrFlags);errno=callerErrno;
        actualWorkers=actualTeam;return true;
    }catch(const std::bad_alloc&){caller.restore();errno=callerErrno;return false;}
    catch(const std::length_error&){caller.restore();errno=callerErrno;return false;}
}
#if !defined(_M_FP_STRICT) || !_M_FP_STRICT
#pragma fenv_access(off)
#endif
#pragma float_control(pop)
#else
template<class Body> bool tryMatrixStage(Eigen::Index,Body&,std::size_t,int&,bool& attempted) noexcept
{attempted=false;return false;}
#endif
template<class Body> void matrixStage(Eigen::Index count,Body body,bool allow,
    std::size_t availableBytes,MatrixStageExecution& execution)
{
    int workers=0;bool attempted=false;
    if(allow&&tryMatrixStage(count,body,availableBytes,workers,attempted)){
        execution.parallelUsed=true;execution.maximumWorkers=std::max(execution.maximumWorkers,workers);return;
    }
    execution.parallelFallback=execution.parallelFallback||attempted;
    for(Eigen::Index k=0;k<count;++k)body(k);
    if(count>0)execution.maximumWorkers=std::max(execution.maximumWorkers,1);
}

struct HydroSparseTransposePlan
{
    using Index=ValidationMatrix::StorageIndex;
    Eigen::Index rows=0;
    std::vector<Index> outer,inner;
    std::vector<std::uint32_t> transpose;
    std::size_t retainedBytes=0;
    bool matches(const ValidationMatrix& matrix) const noexcept
    {
        if(!matrix.isCompressed()||matrix.rows()!=rows||matrix.cols()!=rows
            ||matrix.outerSize()!=rows||matrix.nonZeros()!=static_cast<Eigen::Index>(inner.size()))return false;
        return std::memcmp(matrix.outerIndexPtr(),outer.data(),outer.size()*sizeof(Index))==0
            &&std::memcmp(matrix.innerIndexPtr(),inner.data(),inner.size()*sizeof(Index))==0;
    }
    struct Proof {bool normal=true,exact=true;};
    Proof zeroDefect(const ValidationMatrix& matrix,bool skew) const noexcept
    {
        Proof result;
        for(std::size_t k=0;k<inner.size();++k){
            std::uint64_t a=0,b=0;std::memcpy(&a,matrix.valuePtr()+k,sizeof(a));
            std::memcpy(&b,matrix.valuePtr()+transpose[k],sizeof(b));
            const auto am=a&UINT64_C(0x7fffffffffffffff),bm=b&UINT64_C(0x7fffffffffffffff);
            const auto normalOrZero=[](std::uint64_t value){return value==0
                ||(value>=UINT64_C(0x0010000000000000)&&value<UINT64_C(0x7ff0000000000000));};
            if(!normalOrZero(am)||!normalOrZero(bm))return {false,false};
            // Either sign of exact zero gives an exact zero defect. Nonzero
            // S partners must match bits; W partners must differ only in sign.
            if(am==0&&bm==0)continue;
            if(a!=(skew?(b^UINT64_C(0x8000000000000000)):b))result.exact=false;
        }
        return result;
    }
};
using SaddleMatrix=HydroSinglePhaseSaddle::Matrix;
struct SaddleScatterBinding
{
    const ValidationMatrix &a,&c;
    const std::vector<int> &freeByVelocity,&velocityByFree,&component,&gaugeByComponent;
    const std::vector<HydroSaddlePressureOutlet>& outlets;
    int normalVelocities=0,pressures=0,tangential=0,gauges=0;
};
struct HydroSaddleScatterPlan
{
    using Index=SaddleMatrix::StorageIndex;
    static constexpr std::size_t byteLimit=128u*1024u*1024u;
    Eigen::Index velocityRows=0;
    std::vector<std::uint32_t> aOuter,aInner;
    std::vector<Index> cOuter,cInner,systemOuter,systemInner,quotientOuter,quotientInner;
    std::vector<std::uint32_t> systemSources,quotientSources;
    std::vector<unsigned char> systemKinds;
    bool identityQuotient=false;
    std::vector<int> freeByVelocity,velocityByFree,component,gaugeByComponent,outletRoles;
    std::vector<int> quotientByPressure,pressureByQuotient;
    int normalVelocities=0,pressures=0,tangential=0,gauges=0;
    Eigen::Index systemRows=0,quotientRows=0;
    std::size_t retainedBytes=0;
    template<class MatrixT> static bool words(const MatrixT& matrix,Eigen::Index rows,Eigen::Index cols,
        const std::vector<Index>& outer,const std::vector<Index>& inner) noexcept
    {
        return matrix.isCompressed()&&matrix.outerIndexPtr()!=nullptr&&matrix.rows()==rows&&matrix.cols()==cols
            &&outer.size()==static_cast<std::size_t>(matrix.outerSize()+1)
            &&inner.size()==static_cast<std::size_t>(matrix.nonZeros())
            &&std::memcmp(matrix.outerIndexPtr(),outer.data(),outer.size()*sizeof(Index))==0
            &&(inner.empty()||std::memcmp(matrix.innerIndexPtr(),inner.data(),inner.size()*sizeof(Index))==0);
    }
    template<class MatrixT> static bool ordered(const MatrixT& matrix) noexcept
    {
        if(!matrix.isCompressed()||matrix.rows()<0||matrix.cols()<0||matrix.nonZeros()<0)return false;
        const auto* outer=matrix.outerIndexPtr();const auto* inner=matrix.innerIndexPtr();
        if(!outer||(matrix.nonZeros()!=0&&!inner)||outer[0]!=0||outer[matrix.outerSize()]!=matrix.nonZeros())return false;
        for(Eigen::Index p=0;p<matrix.outerSize();++p){
            if(outer[p]<0||outer[p]>outer[p+1]||outer[p+1]>matrix.nonZeros())return false;
            Index previous=-1;for(Index k=outer[p];k<outer[p+1];++k){
                if(inner[k]<=previous||inner[k]<0||inner[k]>=matrix.innerSize())return false;
                previous=inner[k];
            }
        }
        return true;
    }
    bool actualPattern(const ValidationMatrix& a) const noexcept
    {
        if(!a.isCompressed()||a.rows()!=velocityRows||a.cols()!=velocityRows||!a.outerIndexPtr()
            ||aOuter.size()!=static_cast<std::size_t>(a.outerSize()+1)||aInner.size()!=static_cast<std::size_t>(a.nonZeros()))return false;
        for(std::size_t k=0;k<aOuter.size();++k)if(a.outerIndexPtr()[k]!=aOuter[k])return false;
        for(std::size_t k=0;k<aInner.size();++k)if(a.innerIndexPtr()[k]!=aInner[k])return false;
        return true;
    }
    bool matches(const SaddleScatterBinding& b) const noexcept
    {
        if(!actualPattern(b.a)
            ||!words(b.c,normalVelocities,pressures,cOuter,cInner)
            ||normalVelocities!=b.normalVelocities||pressures!=b.pressures||tangential!=b.tangential||gauges!=b.gauges
            ||freeByVelocity!=b.freeByVelocity||velocityByFree!=b.velocityByFree
            ||component!=b.component||gaugeByComponent!=b.gaugeByComponent||outletRoles.size()!=b.outlets.size())return false;
        for(std::size_t k=0;k<outletRoles.size();++k)if(outletRoles[k]!=b.outlets[k].velocityDof)return false;
        return true;
    }
    bool matchesQuotient(const std::vector<int>& qp,const std::vector<int>& pq) const noexcept
    {return quotientByPressure==qp&&pressureByQuotient==pq;}
    static std::shared_ptr<const HydroSaddleScatterPlan> compile(const SaddleScatterBinding& b,
        int& symbolicWorkers)
    {
        symbolicWorkers=0;
        if(!ordered(b.a)||b.a.rows()!=b.a.cols()||!ordered(b.c)
            ||b.c.rows()!=b.normalVelocities||b.c.cols()!=b.pressures)return {};
        constexpr auto missing=std::numeric_limits<std::uint32_t>::max();
        const auto nf=b.velocityByFree.size(),np=static_cast<std::size_t>(b.pressures);
        const auto n=static_cast<std::size_t>(b.a.rows()),nz=static_cast<std::size_t>(b.a.nonZeros());
        const auto cnz=static_cast<std::size_t>(b.c.nonZeros()),groups=static_cast<std::size_t>(b.gauges);
        if(n>=missing||nz>=missing||cnz>=missing||np>=missing||nf>missing-np||groups>missing-nf-np)return {};
        const auto rows=nf+np+groups;
        std::size_t stock=sizeof(HydroSaddleScatterPlan)+256u;
        const auto add=[&](std::size_t& bytes,std::size_t count,std::size_t width){
            if(bytes>byteLimit||count>(byteLimit-bytes)/width)return false;bytes+=count*width;return true;};
        // Current A/C column routes only contain integers. They are discarded
        // before numerical matrices are allocated; no double snapshot exists.
        std::size_t temporary=8*sizeof(std::vector<std::uint32_t>)+sizeof(std::vector<int>);
        if(!add(temporary,n+1,sizeof(std::uint32_t))||!add(temporary,nz,2*sizeof(std::uint32_t))
            ||!add(temporary,np+1,sizeof(std::uint32_t))||!add(temporary,cnz,2*sizeof(std::uint32_t))
            ||!add(temporary,groups+1,sizeof(std::uint32_t))||!add(temporary,np,sizeof(std::uint32_t)))return {};
        if(!add(stock,n+1,sizeof(std::uint32_t))||!add(stock,nz,sizeof(std::uint32_t))
            ||!add(stock,static_cast<std::size_t>(b.c.outerSize()+1),sizeof(Index))||!add(stock,cnz,sizeof(Index))
            ||!add(stock,rows+1,sizeof(Index)))return {};
        for(const auto* v:{&b.freeByVelocity,&b.velocityByFree,&b.component,&b.gaugeByComponent})
            if(!add(stock,v->size(),sizeof(int)))return {};
        if(!add(stock,b.outlets.size()+2*np,sizeof(int))||!add(temporary,groups,sizeof(int))||temporary>byteLimit-stock)return {};
        std::vector<std::uint32_t> ao(n+1,0),ar(nz),as(nz),co(np+1,0),cr(cnz),cs(cnz),go(groups+1,0),gp(np);
        std::vector<int> representative(groups,-1);
        temporary=8*sizeof(std::vector<std::uint32_t>)+sizeof(representative);
        for(const auto* v:{&ao,&ar,&as,&co,&cr,&cs,&go,&gp})if(!add(temporary,v->capacity(),sizeof(std::uint32_t)))return {};
        if(!add(temporary,representative.capacity(),sizeof(int))||temporary>byteLimit-stock)return {};
        // Reuse a single current cursor; it never overlaps with the per-column
        // worker allocation beyond the explicitly checked cold peak.
        std::size_t cursorBytes=sizeof(std::vector<std::uint32_t>);
        if(!add(cursorBytes,std::max(n,std::max(np,groups)),sizeof(std::uint32_t))
            ||cursorBytes>byteLimit-stock-temporary)return {};
        {
            std::vector<std::uint32_t> next(std::max(n,std::max(np,groups)),0);
            if(next.capacity()>(byteLimit-stock-temporary-sizeof(next))/sizeof(std::uint32_t))return {};
            for(std::size_t k=0;k<nz;++k)++ao[static_cast<std::size_t>(b.a.innerIndexPtr()[k])+1];
            for(std::size_t i=0;i<n;++i){ao[i+1]+=ao[i];next[i]=ao[i];}
            for(std::size_t row=0;row<n;++row)for(Index k=b.a.outerIndexPtr()[row];k<b.a.outerIndexPtr()[row+1];++k){
                const auto slot=next[static_cast<std::size_t>(b.a.innerIndexPtr()[k])]++;
                ar[slot]=static_cast<std::uint32_t>(row);as[slot]=static_cast<std::uint32_t>(k);}
            for(std::size_t k=0;k<cnz;++k)++co[static_cast<std::size_t>(b.c.innerIndexPtr()[k])+1];
            for(std::size_t i=0;i<np;++i){co[i+1]+=co[i];next[i]=co[i];}
            for(Eigen::Index row=0;row<b.c.outerSize();++row)for(Index k=b.c.outerIndexPtr()[row];k<b.c.outerIndexPtr()[row+1];++k){
                const auto slot=next[static_cast<std::size_t>(b.c.innerIndexPtr()[k])]++;
                cr[slot]=static_cast<std::uint32_t>(row);cs[slot]=static_cast<std::uint32_t>(k);}
            for(std::size_t p=0;p<np;++p){const int g=b.gaugeByComponent[b.component[p]];if(g>=0)++go[static_cast<std::size_t>(g)+1];}
            for(std::size_t g=0;g<groups;++g){go[g+1]+=go[g];next[g]=go[g];}
            for(std::size_t p=0;p<np;++p){const int g=b.gaugeByComponent[b.component[p]];if(g>=0)gp[next[static_cast<std::size_t>(g)]++]=static_cast<std::uint32_t>(p);}
        }
        auto p=std::make_shared<HydroSaddleScatterPlan>();p->velocityRows=b.a.rows();
        p->normalVelocities=b.normalVelocities;p->pressures=b.pressures;p->tangential=b.tangential;p->gauges=b.gauges;
        p->freeByVelocity=b.freeByVelocity;p->velocityByFree=b.velocityByFree;p->component=b.component;p->gaugeByComponent=b.gaugeByComponent;
        p->outletRoles.resize(b.outlets.size());
        for(std::size_t k=0;k<b.outlets.size();++k)p->outletRoles[k]=b.outlets[k].velocityDof;
        p->aOuter.resize(n+1);p->aInner.resize(nz);
        for(std::size_t k=0;k<=n;++k)p->aOuter[k]=static_cast<std::uint32_t>(b.a.outerIndexPtr()[k]);
        for(std::size_t k=0;k<nz;++k)p->aInner[k]=static_cast<std::uint32_t>(b.a.innerIndexPtr()[k]);
        p->cOuter.assign(b.c.outerIndexPtr(),b.c.outerIndexPtr()+b.c.outerSize()+1);
        if(cnz)p->cInner.assign(b.c.innerIndexPtr(),b.c.innerIndexPtr()+cnz);
        p->systemRows=static_cast<Eigen::Index>(rows);p->systemOuter.resize(rows+1);
        p->quotientByPressure.assign(np,-1);p->pressureByQuotient.reserve(np);
        for(std::size_t q=0;q<np;++q){const int group=b.gaugeByComponent[b.component[q]];
            if(group>=0&&representative[static_cast<std::size_t>(group)]<0){representative[static_cast<std::size_t>(group)]=static_cast<int>(q);continue;}
            p->quotientByPressure[q]=static_cast<int>(p->pressureByQuotient.size());p->pressureByQuotient.push_back(static_cast<int>(q));}
        // Check actual base capacities before allocating the O(nnz) final routes.
        stock=sizeof(HydroSaddleScatterPlan)+256u;
        for(const auto* v:{&p->cOuter,&p->cInner,&p->systemOuter})if(!add(stock,v->capacity(),sizeof(Index)))return {};
        for(const auto* v:{&p->aOuter,&p->aInner})if(!add(stock,v->capacity(),sizeof(std::uint32_t)))return {};
        for(const auto* v:{&p->freeByVelocity,&p->velocityByFree,&p->component,&p->gaugeByComponent,&p->outletRoles,&p->quotientByPressure,&p->pressureByQuotient})if(!add(stock,v->capacity(),sizeof(int)))return {};
        if(temporary>byteLimit-stock)return {};
        // Columns are emitted in exactly increasing row order: free A, then C;
        // free C-transpose, then gauge; and finally the ordered gauge column.
        const auto visit=[&](Eigen::Index column,auto&& emit) noexcept {
            const auto col=static_cast<std::size_t>(column);
            if(col<nf){const int u=b.velocityByFree[col];
                for(auto k=ao[static_cast<std::size_t>(u)];k<ao[static_cast<std::size_t>(u)+1];++k){const int row=b.freeByVelocity[ar[k]];
                    if(row>=0)emit(static_cast<Index>(row),static_cast<unsigned char>(0),as[k]);}
                if(u<b.normalVelocities)for(Index k=b.c.outerIndexPtr()[u];k<b.c.outerIndexPtr()[u+1];++k)
                    emit(static_cast<Index>(nf)+b.c.innerIndexPtr()[k],static_cast<unsigned char>(1),static_cast<std::uint32_t>(k));
            }else if(col<nf+np){const auto q=col-nf;
                for(auto k=co[q];k<co[q+1];++k){const int row=b.freeByVelocity[cr[k]];
                    if(row>=0)emit(static_cast<Index>(row),static_cast<unsigned char>(1),cs[k]);}
                const int group=b.gaugeByComponent[b.component[q]];
                if(group>=0)emit(static_cast<Index>(nf+np)+group,static_cast<unsigned char>(2),static_cast<std::uint32_t>(q));
            }else{const auto group=col-nf-np;
                for(auto k=go[group];k<go[group+1];++k)emit(static_cast<Index>(nf)+gp[k],static_cast<unsigned char>(2),gp[k]);}
        };
        MatrixStageExecution countExecution;
        matrixStage(static_cast<Eigen::Index>(rows),[&](Eigen::Index col) noexcept {
            Index count=0;visit(col,[&](Index,unsigned char,std::uint32_t) noexcept {++count;});
            p->systemOuter[static_cast<std::size_t>(col)+1]=count;
        },true,byteLimit-stock-temporary,countExecution);
        p->systemOuter[0]=0;
        for(std::size_t col=0;col<rows;++col){
            if(p->systemOuter[col+1]>static_cast<Index>(missing)-p->systemOuter[col])return {};
            p->systemOuter[col+1]+=p->systemOuter[col];}
        const auto ns=static_cast<std::size_t>(p->systemOuter[rows]);
        if(!add(stock,ns,sizeof(Index)+sizeof(std::uint32_t)+sizeof(unsigned char))||temporary>byteLimit-stock)return {};
        p->systemInner.resize(ns);p->systemSources.resize(ns);p->systemKinds.resize(ns);
        // Reconcile any allocator capacity rounding before the worker scratch.
        const auto systemRequested=ns*(sizeof(Index)+sizeof(std::uint32_t)+sizeof(unsigned char));
        stock-=systemRequested;
        if(!add(stock,p->systemInner.capacity(),sizeof(Index))||!add(stock,p->systemSources.capacity(),sizeof(std::uint32_t))
            ||!add(stock,p->systemKinds.capacity(),sizeof(unsigned char))||temporary>byteLimit-stock)return {};
        std::atomic<bool> valid{true};MatrixStageExecution fillExecution;
        matrixStage(static_cast<Eigen::Index>(rows),[&](Eigen::Index column) noexcept {
            const auto col=static_cast<std::size_t>(column);Index slot=p->systemOuter[col],previous=-1;
            visit(column,[&](Index row,unsigned char kind,std::uint32_t source) noexcept {
                if(row<0||row>=p->systemRows||row<=previous||slot>=p->systemOuter[col+1]
                    ||kind>2||(kind==0&&source>=nz)||(kind==1&&source>=cnz)||(kind==2&&source>=np)){
                    valid.store(false,std::memory_order_relaxed);return;}
                previous=row;p->systemInner[static_cast<std::size_t>(slot)]=row;
                p->systemSources[static_cast<std::size_t>(slot)]=source;p->systemKinds[static_cast<std::size_t>(slot)]=kind;++slot;
            });
            if(slot!=p->systemOuter[col+1])valid.store(false,std::memory_order_relaxed);
        },true,byteLimit-stock-temporary,fillExecution);
        if(!valid.load(std::memory_order_relaxed))return {};
        symbolicWorkers=std::min(countExecution.maximumWorkers,fillExecution.maximumWorkers);
        p->quotientRows=static_cast<Eigen::Index>(nf+p->pressureByQuotient.size());
        p->identityQuotient=groups==0&&p->quotientRows==p->systemRows;
        if(p->identityQuotient){for(std::size_t q=0;q<np;++q)
            if(p->quotientByPressure[q]!=static_cast<int>(q)||p->pressureByQuotient[q]!=static_cast<int>(q))return {};}
        else{
            const auto reduced=[&](Index row)->Index {
                if(row<static_cast<Index>(nf))return row;
                if(row<static_cast<Index>(nf+np)){const int q=p->quotientByPressure[static_cast<std::size_t>(row)-nf];return q<0?-1:static_cast<Index>(nf)+q;}return -1;};
            std::size_t nq=0;
            for(std::size_t col=0;col<rows;++col)if(reduced(static_cast<Index>(col))>=0)
                for(Index k=p->systemOuter[col];k<p->systemOuter[col+1];++k)if(reduced(p->systemInner[static_cast<std::size_t>(k)])>=0)++nq;
            if(!add(stock,static_cast<std::size_t>(p->quotientRows)+1,sizeof(Index))
                ||!add(stock,nq,sizeof(Index)+sizeof(std::uint32_t))||temporary>byteLimit-stock)return {};
            p->quotientOuter.resize(static_cast<std::size_t>(p->quotientRows)+1);p->quotientInner.reserve(nq);p->quotientSources.reserve(nq);
            stock-=(static_cast<std::size_t>(p->quotientRows)+1)*sizeof(Index)+nq*(sizeof(Index)+sizeof(std::uint32_t));
            if(!add(stock,p->quotientOuter.capacity(),sizeof(Index))||!add(stock,p->quotientInner.capacity(),sizeof(Index))
                ||!add(stock,p->quotientSources.capacity(),sizeof(std::uint32_t))||temporary>byteLimit-stock)return {};
            for(std::size_t col=0;col<rows;++col){const auto target=reduced(static_cast<Index>(col));if(target<0)continue;
                p->quotientOuter[static_cast<std::size_t>(target)]=static_cast<Index>(p->quotientInner.size());
                for(Index k=p->systemOuter[col];k<p->systemOuter[col+1];++k){const auto row=reduced(p->systemInner[static_cast<std::size_t>(k)]);if(row<0)continue;
                    p->quotientInner.push_back(row);p->quotientSources.push_back(static_cast<std::uint32_t>(k));}}
            p->quotientOuter.back()=static_cast<Index>(p->quotientInner.size());
        }
        std::size_t bytes=sizeof(HydroSaddleScatterPlan)+256u;
        for(const auto* v:{&p->cOuter,&p->cInner,&p->systemOuter,&p->systemInner,&p->quotientOuter,&p->quotientInner})if(!add(bytes,v->capacity(),sizeof(Index)))return {};
        for(const auto* v:{&p->aOuter,&p->aInner,&p->systemSources,&p->quotientSources})if(!add(bytes,v->capacity(),sizeof(std::uint32_t)))return {};
        if(!add(bytes,p->systemKinds.capacity(),sizeof(unsigned char)))return {};
        for(const auto* v:{&p->freeByVelocity,&p->velocityByFree,&p->component,&p->gaugeByComponent,&p->outletRoles,&p->quotientByPressure,&p->pressureByQuotient})if(!add(bytes,v->capacity(),sizeof(int)))return {};
        std::size_t actualTemporary=8*sizeof(std::vector<std::uint32_t>)+sizeof(representative);
        for(const auto* v:{&ao,&ar,&as,&co,&cr,&cs,&go,&gp})if(!add(actualTemporary,v->capacity(),sizeof(std::uint32_t)))return {};
        if(!add(actualTemporary,representative.capacity(),sizeof(int))||actualTemporary>byteLimit-bytes)return {};
        p->retainedBytes=bytes;return p;
    }
};
// Pure pattern storage: no matrix values, state, right-hand side or factor.
struct HydroVelocitySplitIndexPlan
{
    using Index=ValidationMatrix::StorageIndex;
    static constexpr std::size_t byteLimit=128u*1024u*1024u;
    static constexpr std::uint32_t missing=std::numeric_limits<std::uint32_t>::max();
    Eigen::Index rows=0;
    std::vector<std::uint32_t> inputOuter,inputInner,source,transposedSource;
    std::vector<Index> unionOuter,unionInner;
    std::size_t retainedBytes=0;
    bool matches(const ValidationMatrix& a) const noexcept
    {
        if(!a.isCompressed()||a.rows()!=rows||a.cols()!=rows||!a.outerIndexPtr()
            ||a.nonZeros()!=static_cast<Eigen::Index>(inputInner.size()))return false;
        for(std::size_t k=0;k<inputOuter.size();++k)
            if(a.outerIndexPtr()[k]!=static_cast<Index>(inputOuter[k]))return false;
        for(std::size_t k=0;k<inputInner.size();++k)
            if(a.innerIndexPtr()[k]!=static_cast<Index>(inputInner[k]))return false;
        return true;
    }
    static bool normalValues(const ValidationMatrix& a) noexcept
    {
        for(Eigen::Index k=0;k<a.nonZeros();++k){
            std::uint64_t bits=0;std::memcpy(&bits,a.valuePtr()+k,sizeof(bits));
            const auto magnitude=bits&UINT64_C(0x7fffffffffffffff);
            if(magnitude!=0&&(magnitude<UINT64_C(0x0010000000000000)
                ||magnitude>=UINT64_C(0x7ff0000000000000)))return false;
        }
        return true;
    }
    static std::shared_ptr<const HydroVelocitySplitIndexPlan> compile(
        const ValidationMatrix& a,const ValidationMatrix& unionPattern)
    {
        if(!HydroSaddleScatterPlan::ordered(a)||a.rows()<=0||a.cols()!=a.rows()
            ||static_cast<std::uint64_t>(a.rows())>missing||static_cast<std::uint64_t>(a.nonZeros())>missing)return {};
        if(!HydroSaddleScatterPlan::ordered(unionPattern)
            ||unionPattern.rows()!=a.rows()||unionPattern.cols()!=a.rows())return {};
        const auto n=static_cast<std::size_t>(a.rows()),nz=static_cast<std::size_t>(a.nonZeros());
        constexpr std::size_t overhead=sizeof(HydroVelocitySplitIndexPlan)+256u;
        // Integer-only cold preparation uses the SAME transpose/union plan.
        // No numerical transpose, coefficient snapshot or second graph owner.
        // All source indices fit uint32_t; the output CSR retains Eigen Index.
        std::size_t estimate=overhead;
        const auto estimatePart=[&](std::size_t count,std::size_t width){
            if(count>(byteLimit-estimate)/width)return false;estimate+=count*width;return true;};
        if(!estimatePart(n+1,2*sizeof(std::uint32_t))||!estimatePart(nz,2*sizeof(std::uint32_t)))return {};
        std::vector<std::uint32_t> transposeOuter(n+1,0),transposeInner(nz),transposeSlots(nz);
        for(std::size_t k=0;k<nz;++k)++transposeOuter[static_cast<std::size_t>(a.innerIndexPtr()[k])+1];
        for(std::size_t row=0;row<n;++row)transposeOuter[row+1]+=transposeOuter[row];
        auto next=transposeOuter;
        for(std::size_t row=0;row<n;++row)
            for(Index k=a.outerIndexPtr()[row];k<a.outerIndexPtr()[row+1];++k){
                const auto column=static_cast<std::size_t>(a.innerIndexPtr()[k]);const auto slot=next[column]++;
                transposeInner[slot]=static_cast<std::uint32_t>(row);transposeSlots[slot]=static_cast<std::uint32_t>(k);
            }
        std::size_t unz=0;
        for(std::size_t row=0;row<n;++row){
            Index left=a.outerIndexPtr()[row];std::uint32_t right=transposeOuter[row];
            const Index leftEnd=a.outerIndexPtr()[row+1];const auto rightEnd=transposeOuter[row+1];
            while(left<leftEnd||right<rightEnd){
                const Index column=std::min(left<leftEnd?a.innerIndexPtr()[left]:a.rows(),
                    right<rightEnd?static_cast<Index>(transposeInner[right]):a.rows());
                if(unz==missing)return {};++unz;
                if(left<leftEnd&&a.innerIndexPtr()[left]==column)++left;
                if(right<rightEnd&&transposeInner[right]==static_cast<std::uint32_t>(column))++right;
            }
        }
        if(unionPattern.nonZeros()!=static_cast<Eigen::Index>(unz))return {};
        estimate=overhead;
        if(!estimatePart(n+1,sizeof(std::uint32_t)+sizeof(Index))||!estimatePart(nz,sizeof(std::uint32_t))
            ||!estimatePart(unz,sizeof(Index)+2*sizeof(std::uint32_t)))return {};
        const auto temporaryBytes=[&](std::size_t& bytes){
            constexpr auto metadata=4*sizeof(std::vector<std::uint32_t>);
            if(metadata>byteLimit-bytes)return false;bytes+=metadata;
            for(const auto* v:{&transposeOuter,&transposeInner,&transposeSlots,&next}){
                if(v->capacity()>(byteLimit-bytes)/sizeof(std::uint32_t))return false;
                bytes+=v->capacity()*sizeof(std::uint32_t);
            }return true;};
        auto peak=estimate;if(!temporaryBytes(peak))return {};
        auto result=std::make_shared<HydroVelocitySplitIndexPlan>();result->rows=a.rows();
        result->inputOuter.resize(n+1);result->inputInner.resize(nz);
        result->unionOuter.resize(n+1);result->unionInner.resize(unz);
        result->source.resize(unz);result->transposedSource.resize(unz);
        std::size_t bytes=overhead;
        const auto account=[&](std::size_t count,std::size_t width){
            if(count>(byteLimit-bytes)/width)return false;bytes+=count*width;return true;};
        if(!account(result->inputOuter.capacity(),sizeof(std::uint32_t))
            ||!account(result->inputInner.capacity(),sizeof(std::uint32_t))
            ||!account(result->unionOuter.capacity(),sizeof(Index))||!account(result->unionInner.capacity(),sizeof(Index))
            ||!account(result->source.capacity(),sizeof(std::uint32_t))
            ||!account(result->transposedSource.capacity(),sizeof(std::uint32_t)))return {};
        peak=bytes;if(!temporaryBytes(peak))return {};
        for(std::size_t k=0;k<n+1;++k)result->inputOuter[k]=static_cast<std::uint32_t>(a.outerIndexPtr()[k]);
        for(std::size_t k=0;k<nz;++k)result->inputInner[k]=static_cast<std::uint32_t>(a.innerIndexPtr()[k]);
        std::size_t slot=0;
        for(std::size_t row=0;row<n;++row){
            result->unionOuter[row]=static_cast<Index>(slot);
            if(unionPattern.outerIndexPtr()[row]!=static_cast<Index>(slot))return {};
            Index left=a.outerIndexPtr()[row];std::uint32_t right=transposeOuter[row];
            const Index leftEnd=a.outerIndexPtr()[row+1];const auto rightEnd=transposeOuter[row+1];
            while(left<leftEnd||right<rightEnd){
                const Index column=std::min(left<leftEnd?a.innerIndexPtr()[left]:a.rows(),
                    right<rightEnd?static_cast<Index>(transposeInner[right]):a.rows());
                const bool hasLeft=left<leftEnd&&a.innerIndexPtr()[left]==column;
                const bool hasRight=right<rightEnd&&transposeInner[right]==static_cast<std::uint32_t>(column);
                if(slot>=unz||unionPattern.innerIndexPtr()[slot]!=column)return {};
                result->unionInner[slot]=column;
                result->source[slot]=hasLeft?static_cast<std::uint32_t>(left++):missing;
                result->transposedSource[slot]=hasRight?transposeSlots[right++]:missing;++slot;
            }
        }
        result->unionOuter[n]=static_cast<Index>(slot);
        if(slot!=unz||unionPattern.outerIndexPtr()[n]!=static_cast<Index>(slot))return {};
        result->retainedBytes=bytes;return result;
    }
    using Matrices=std::pair<ValidationMatrix,ValidationMatrix>;
    std::unique_ptr<Matrices> freshMatrices() const noexcept
    {
        SparseValidationFloatingRestore restore;if(!restore.captured)return {};
        try{
            auto result=std::make_unique<Matrices>(ValidationMatrix(rows,rows),ValidationMatrix(rows,rows));
            for(auto* matrix:{&result->first,&result->second}){
                matrix->resizeNonZeros(static_cast<Eigen::Index>(unionInner.size()));
                std::memcpy(matrix->outerIndexPtr(),unionOuter.data(),unionOuter.size()*sizeof(Index));
                if(!unionInner.empty())std::memcpy(matrix->innerIndexPtr(),unionInner.data(),unionInner.size()*sizeof(Index));
            }
            return result;
        }catch(...){return {};}
    }
    void freshValues(const ValidationMatrix& a,Matrices& result) const
    {
        const Eigen::internal::scalar_sum_op<double,double> sum;
        const Eigen::internal::scalar_difference_op<double,double> difference;
        const Eigen::internal::scalar_product_op<double,double> multiply;
        // Explicit public certificate: complete original S, then complete W.
        for(std::size_t k=0;k<unionInner.size();++k){
            const double left=source[k]==missing?0.:a.valuePtr()[source[k]];
            const double right=transposedSource[k]==missing?0.:a.valuePtr()[transposedSource[k]];
            volatile double rounded=sum(left,right);result.first.valuePtr()[k]=multiply(.5,static_cast<double>(rounded));
        }
        for(std::size_t k=0;k<unionInner.size();++k){
            const double left=source[k]==missing?0.:a.valuePtr()[source[k]];
            const double right=transposedSource[k]==missing?0.:a.valuePtr()[transposedSource[k]];
            volatile double rounded=difference(left,right);result.second.valuePtr()[k]=multiply(.5,static_cast<double>(rounded));
        }
    }
};
struct HydroSparseValidationPlanCache
{
    static constexpr std::size_t byteLimit=64u*1024u*1024u;
    std::mutex mutex;
    std::shared_ptr<const HydroSparseTransposePlan> plan;
    std::shared_ptr<const HydroSaddleScatterPlan> scatter;
    std::shared_ptr<const HydroVelocitySplitIndexPlan> velocitySplit;
    std::shared_ptr<const HydroSparseTransposePlan> findOrPrepare(const ValidationMatrix& matrix,bool& hit)
    {
        std::lock_guard<std::mutex> lock(mutex);hit=false;
        if(plan){if(plan->matches(matrix)){hit=true;return plan;}return {};}
        using Index=HydroSparseTransposePlan::Index;
        const auto n=matrix.rows(),nz=matrix.nonZeros();
        if(!matrix.isCompressed()||n<=0||matrix.cols()!=n||matrix.outerSize()!=n||nz<=0
            ||static_cast<std::uint64_t>(nz)>std::numeric_limits<std::uint32_t>::max())return {};
        const auto count=static_cast<std::size_t>(nz),size=static_cast<std::size_t>(n);
        constexpr auto overhead=sizeof(HydroSparseTransposePlan)+sizeof(HydroSparseValidationPlanCache)+256u;
        if(size>=byteLimit/sizeof(Index)||count>byteLimit/(sizeof(Index)+sizeof(std::uint32_t)))return {};
        const auto estimate=overhead+(size+1)*sizeof(Index)+count*(sizeof(Index)+sizeof(std::uint32_t));
        if(estimate>byteLimit)return {};
        const Index* outer=matrix.outerIndexPtr();const Index* inner=matrix.innerIndexPtr();
        if(outer[0]!=0||outer[n]!=nz)return {};
        for(Eigen::Index row=0;row<n;++row){
            if(outer[row]<0||outer[row]>outer[row+1]||outer[row+1]>nz)return {};
            Index previous=-1;
            for(Index k=outer[row];k<outer[row+1];++k){
                if(inner[k]<0||inner[k]>=n||inner[k]<=previous)return {};
                previous=inner[k];
            }
        }
        auto candidate=std::make_shared<HydroSparseTransposePlan>();candidate->rows=n;
        candidate->outer.assign(outer,outer+n+1);candidate->inner.assign(inner,inner+nz);
        candidate->transpose.resize(count);
        std::size_t bytes=overhead;
        const auto account=[&](std::size_t count,std::size_t width){
            if(count>(byteLimit-bytes)/width)return false;
            bytes+=count*width;return true;
        };
        if(!account(candidate->outer.capacity(),sizeof(Index))||!account(candidate->inner.capacity(),sizeof(Index))
            ||!account(candidate->transpose.capacity(),sizeof(std::uint32_t)))return {};
        // CSR rows are sorted. Visiting source rows in order also visits the
        // transposed inner labels in order: no numeric values or sorting keys.
        std::vector<Index> next(outer,outer+n);
        for(Eigen::Index row=0;row<n;++row)
            for(Index k=outer[row];k<outer[row+1];++k){
                const Index column=inner[k],partner=next[column];
                if(partner>=outer[column+1]||inner[partner]!=row)return {};
                candidate->transpose[static_cast<std::size_t>(k)]=static_cast<std::uint32_t>(partner);
                next[column]=partner+1;
            }
        for(Eigen::Index row=0;row<n;++row)if(next[row]!=outer[row+1])return {};
        for(std::size_t k=0;k<count;++k)
            if(candidate->transpose[candidate->transpose[k]]!=k)return {};
        candidate->retainedBytes=bytes;plan=std::move(candidate);return plan;
    }
};
bool velocitySplitQualified() noexcept
{
#if EIGEN_WORLD_VERSION==3 && EIGEN_MAJOR_VERSION==4 && EIGEN_MINOR_VERSION==1 && !defined(EIGEN_SCALAR_BINARY_OP_PLUGIN)
    return sparseValidationFloatingQualified();
#else
    return false;
#endif
}
std::shared_ptr<const HydroVelocitySplitIndexPlan> findVelocitySplit(
    const std::shared_ptr<HydroSparseValidationPlanCache>& cache,const ValidationMatrix& a,
    bool sourceQualified,HydroVelocitySplitPlanReport& report) noexcept
{
    SparseValidationFloatingRestore restore;
    if(!restore.captured||!sourceQualified||!cache||!velocitySplitQualified())return {};
    try{
        std::shared_ptr<const HydroVelocitySplitIndexPlan> plan;
        {std::lock_guard<std::mutex> lock(cache->mutex);plan=cache->velocitySplit;}
        if(!plan||!plan->matches(a))return {};
        report.planHit=true;report.retainedBytes=plan->retainedBytes;
        report.eligible=HydroVelocitySplitIndexPlan::normalValues(a);
        return report.eligible?plan:std::shared_ptr<const HydroVelocitySplitIndexPlan>{};
    }catch(...){return {};}
}
void publishVelocitySplit(const std::shared_ptr<HydroSparseValidationPlanCache>& cache,
    const ValidationMatrix& a,const ValidationMatrix& unionPattern,bool sourceQualified,
    HydroVelocitySplitPlanReport& report) noexcept
{
    SparseValidationFloatingRestore restore;
    if(!restore.captured||!sourceQualified||!cache||!velocitySplitQualified()
        ||!HydroVelocitySplitIndexPlan::normalValues(a))return;
    try{
        {std::lock_guard<std::mutex> lock(cache->mutex);
            if(cache->velocitySplit&&cache->velocitySplit->matches(a)){
                report.eligible=true;report.retainedBytes=cache->velocitySplit->retainedBytes;return;}
            // One current input structure, not an accumulating pattern history.
            cache->velocitySplit.reset();}
        auto candidate=HydroVelocitySplitIndexPlan::compile(a,unionPattern);if(!candidate)return;
        std::lock_guard<std::mutex> lock(cache->mutex);
        if(!cache->velocitySplit||!cache->velocitySplit->matches(a))cache->velocitySplit=std::move(candidate);
        report.eligible=true;report.retainedBytes=cache->velocitySplit->retainedBytes;
    }catch(...){return;}
}
bool saddleScatterQualified() noexcept
{
#if EIGEN_WORLD_VERSION==3 && EIGEN_MAJOR_VERSION==4 && EIGEN_MINOR_VERSION==1
    return sparseValidationFloatingQualified();
#else
    return false;
#endif
}
std::shared_ptr<const HydroSaddleScatterPlan> prepareSaddleScatter(
    const std::shared_ptr<HydroSparseValidationPlanCache>& cache,const SaddleScatterBinding& binding,
    bool sourceQualified,HydroSaddleAssemblyPlanReport& report) noexcept
{
    SparseValidationFloatingRestore restore;
    if(!restore.captured||!sourceQualified||!cache||!saddleScatterQualified())return {};
    try{
        {std::lock_guard<std::mutex> lock(cache->mutex);
            if(cache->scatter&&cache->scatter->matches(binding)){
                report.eligible=true;report.planHit=true;report.retainedBytes=cache->scatter->retainedBytes;return cache->scatter;}
            // One current pattern, as before: release stale stock before cold scratch.
            cache->scatter.reset();}
        int workers=0;auto candidate=HydroSaddleScatterPlan::compile(binding,workers);
        if(candidate){report.eligible=true;report.symbolicBuilt=true;report.symbolicWorkers=workers;report.retainedBytes=candidate->retainedBytes;}
        return candidate;
    }catch(...){return {};}
}
void publishSaddleScatter(const std::shared_ptr<HydroSparseValidationPlanCache>& cache,
    const std::shared_ptr<const HydroSaddleScatterPlan>& candidate) noexcept
{
    SparseValidationFloatingRestore restore;if(!restore.captured||!cache||!candidate)return;
    try{std::lock_guard<std::mutex> lock(cache->mutex);cache->scatter=candidate;}catch(...){return;}
}
// The authenticated plan owns integer structure only. Allocate the fresh final
// CSC before the original scatter, so each current coefficient is written once,
// without an intermediate double array followed by a Map-to-Matrix value copy.
std::unique_ptr<SaddleMatrix> saddleFreshMatrix(Eigen::Index rows,
    const std::vector<SaddleMatrix::StorageIndex>& outer,
    const std::vector<SaddleMatrix::StorageIndex>& inner) noexcept
{
    SparseValidationFloatingRestore restore;if(!restore.captured)return {};
    if(rows<0||outer.size()!=static_cast<std::size_t>(rows)+1
        ||inner.size()>static_cast<std::size_t>(std::numeric_limits<Eigen::Index>::max()))return {};
    try{
        auto matrix=std::make_unique<SaddleMatrix>(rows,rows);
        matrix->resizeNonZeros(static_cast<Eigen::Index>(inner.size()));
        std::copy(outer.begin(),outer.end(),matrix->outerIndexPtr());
        std::copy(inner.begin(),inner.end(),matrix->innerIndexPtr());
        return matrix;
    }catch(const std::bad_alloc&){return {};}
    catch(const std::length_error&){return {};}
}
class SparseValidationContext
{
    std::shared_ptr<HydroSparseValidationPlanCache> cache_;
    const ValidationMatrix *symmetric_,*skew_;
    std::shared_ptr<const HydroSparseTransposePlan> appendPlan_;
    std::unique_ptr<double[]> appendDifference_;
    bool sourceQualified_=false,prepared_=false,symmetricZero_=false,skewZero_=false;
    HydroSparseValidationReport& report_;
    void prepare() noexcept
    {
        if(prepared_)return;prepared_=true;
        SparseValidationFloatingRestore restore;
        if(!restore.captured||!sourceQualified_||!cache_||!sparseValidationFloatingQualified())return;
        try{
            bool hit=false;const auto plan=cache_->findOrPrepare(*symmetric_,hit);
            report_.planHit=hit;if(!plan)return;
            report_.retainedBytes=plan->retainedBytes;
            if(!plan->matches(*skew_))return;
            const auto s=plan->zeroDefect(*symmetric_,false),w=plan->zeroDefect(*skew_,true);
            report_.eligible=s.normal&&w.normal;if(!report_.eligible)return;
            symmetricZero_=s.exact;skewZero_=w.exact;appendPlan_=plan;
        }catch(...){report_.eligible=false;} // Optional execution/index work only.
    }
public:
    SparseValidationContext(std::shared_ptr<HydroSparseValidationPlanCache> cache,
        const ValidationMatrix& symmetric,const ValidationMatrix& skew,bool sourceQualified,bool requested,
        HydroSparseValidationReport& report)
        :cache_(std::move(cache)),symmetric_(&symmetric),skew_(&skew),sourceQualified_(sourceQualified),report_(report)
    {report_.requested=requested;if(requested)report_.nonZeros=static_cast<std::size_t>(symmetric.nonZeros());}
    bool symmetryZero(const ValidationMatrix& matrix) noexcept
    {
        if(&matrix!=symmetric_)return false;prepare();
        if(!report_.eligible||!symmetricZero_||!sparseValidationFloatingQualified())return false;
        report_.exactSymmetryUsed=true;return true;
    }
    bool skewZero(const ValidationMatrix& matrix) noexcept
    {
        if(&matrix!=skew_||!prepared_||!report_.eligible||!skewZero_||!sparseValidationFloatingQualified())return false;
        report_.exactSkewUsed=true;return true;
    }
    ~SparseValidationContext() noexcept
    {
        // Capture the CURRENT numerical environment only for optional cleanup.
        // No sum/difference/norm/gate is enclosed by this restoring scope.
        SparseValidationFloatingRestore restore;
        appendDifference_.reset();appendPlan_.reset();
    }
    bool appendSkew(ValidationMatrix& integral,const ValidationMatrix& skew,double& symmetry)
    {
        if(&integral!=symmetric_||&skew!=skew_||!prepared_||!report_.eligible
            ||!report_.exactSkewUsed||!appendPlan_)return false;
        const auto count=static_cast<std::size_t>(integral.nonZeros());
        {
            SparseValidationFloatingRestore restore;
            constexpr std::size_t scratchLimit=64u*1024u*1024u;
            if(!restore.captured||!velocitySplitQualified()
                ||!appendPlan_->matches(integral)||!appendPlan_->matches(skew)
                ||count==0||count>scratchLimit/sizeof(double))return false;
            try{appendDifference_=std::unique_ptr<double[]>(new double[count]);}
            catch(...){return false;} // No numerical operation or mutation yet.
        }
        const Eigen::internal::scalar_sum_op<double,double> sum;
        const Eigen::internal::scalar_difference_op<double,double> difference;
        // The validated S/W share the complete transpose-closed CSR union.
        // Eigen 3.4.1 evaluates all S+W before the original full finite scan.
        for(std::size_t k=0;k<count;++k)
            integral.valuePtr()[k]=sum(integral.valuePtr()[k],skew.valuePtr()[k]);
        for(Eigen::Index k=0;k<integral.nonZeros();++k)
            require(std::isfinite(integral.valuePtr()[k]),"nonfinite combined velocity block");
        // Preserve rounded S+W, explicit +/-0, original transpose subtraction
        // operand order and the complete row/inner difference expression.
        for(std::size_t k=0;k<count;++k)
            appendDifference_[k]=difference(integral.valuePtr()[k],
                integral.valuePtr()[appendPlan_->transpose[k]]);
        const Eigen::Map<const ValidationMatrix> defect(integral.rows(),integral.cols(),integral.nonZeros(),
            integral.outerIndexPtr(),integral.innerIndexPtr(),appendDifference_.get());
        // Both norms use SparseDot -> cwiseAbs2().sum(), the original sparse
        // scalar row/inner reduction, followed by the same sqrt and ratio.
        symmetry=ratio(defect.norm(),integral.norm());
        report_.appendPlanUsed=true;report_.appendScratchBytes=count*sizeof(double);
        return true;
    }
};
}
std::shared_ptr<Detail::HydroSparseValidationPlanCache> HydroMomentumGeometry::makeSparseValidationPlans() noexcept
{
    Detail::SparseValidationFloatingRestore restore;
    if(!restore.captured)return {};
    try{return std::make_shared<Detail::HydroSparseValidationPlanCache>();}
    catch(...){return {};}
}
struct HydroSaddleVelocityBlock::Data
{
    Data(HydroDofLayout input, Matrix matrix) : layout(std::move(input)),integral(std::move(matrix)) {}
    // A unique prepared matrix installs its fresh storage without a copy.
    struct SwapOwnedSplit {};
    Data(HydroDofLayout input,Matrix& fresh,SwapOwnedSplit) : layout(std::move(input))
    {integral.swap(fresh);}
    HydroDofLayout layout;
    Matrix integral;
    double symmetry=0,validationSeconds=0;
    int tangential=0;
    HydroVelocitySpdReport spd;
    HydroSparseValidationReport sparseValidation;
    HydroVelocitySplitPlanReport splitPlan;
    void validateIntegral(const EnergySpdInputs* energy=nullptr,Detail::SparseValidationContext* sparse=nullptr)
    {
        bool diagonal=true;Eigen::Index positiveDiagonals=0;
        for (Eigen::Index outer=0;outer<integral.outerSize();++outer)
            for (Matrix::InnerIterator it(integral,outer);it;++it) {
                require(std::isfinite(it.value()),"nonfinite velocity block");
                if(it.row()==it.col()){if(it.value()>0)++positiveDiagonals;}
                else if(it.value()!=0)diagonal=false;
            }
        if(sparse && sparse->symmetryZero(integral)) {
            const double scale=integral.norm();
            require(std::isfinite(scale) && scale>0,"unresolved velocity block norm");
            symmetry=ratio(0.,scale);
            require(symmetry<=1e-12,"velocity block is not symmetric");
        } else {
            const Matrix difference=integral-Matrix(integral.transpose());
            const double scale=integral.norm();
            require(std::isfinite(scale) && scale>0,"unresolved velocity block norm");
            symmetry=ratio(difference.norm(),scale);
            require(symmetry<=1e-12,"velocity block is not symmetric");
        }
        // An exactly diagonal matrix is SPD iff every diagonal is positive. The
        // compressed matrix has at most one stored diagonal per row; keep all
        // caller entries, including explicit zeros, and the original gates above.
        if(diagonal) {
            require(positiveDiagonals==integral.rows(),"velocity block is not positive definite");
            spd.tier=HydroVelocitySpdTier::Diagonal;
            return;
        }
        const auto certificateStarted=WallClock::now();
        if(strictLowerPositiveDominance(integral)){
            spd.certificateSeconds=elapsedSeconds(certificateStarted);
            spd.tier=HydroVelocitySpdTier::StrictDiagonalDominance;
            return;
        }
        if(energy){
            spd.gramAttempted=true;
            if(energyLowerPositiveCertificate(integral,*energy,spd.gramBoundUpper,spd.scaledMassLower,spd.gramBoundsAvailable)){
                spd.certificateSeconds=elapsedSeconds(certificateStarted);
                spd.tier=HydroVelocitySpdTier::EnergyGram;
                return;
            }
        }
        spd.certificateSeconds=elapsedSeconds(certificateStarted);
        const auto lltStarted=WallClock::now();
        Eigen::SimplicialLLT<Sparse> check;
        check.compute(Sparse(integral));
        spd.lltSeconds=elapsedSeconds(lltStarted);spd.tier=HydroVelocitySpdTier::LLT;
        require(check.info()==Eigen::Success,"velocity block is not positive definite");
        const Sparse factor=check.matrixL();
        for (Eigen::Index i=0;i<factor.nonZeros();++i)
            require(std::isfinite(factor.valuePtr()[i]),"nonfinite velocity block factorization");
    }
    static void validateSkew(const Matrix& integral,const Matrix& skew,Detail::SparseValidationContext* sparse=nullptr)
    {
        require(skew.rows()==integral.rows()&&skew.cols()==integral.cols(),"different skew block dimensions");
        for(Eigen::Index i=0;i<skew.nonZeros();++i)require(std::isfinite(skew.valuePtr()[i]),"nonfinite skew transport block");
        const double scale=skew.norm();require(std::isfinite(scale),"nonfinite skew block norm");
        if(sparse && sparse->skewZero(skew)) {
            require(ratio(0.,scale)<=1e-12,"transport block is not skew; material/boundary pairing required");
            return;
        }
        const Matrix defect=skew+Matrix(skew.transpose());
        require(ratio(defect.norm(),scale)<=1e-12,"transport block is not skew; material/boundary pairing required");
    }
    void appendSkew(const Matrix& skew,Detail::SparseValidationContext* sparse=nullptr)
    {
        if(sparse&&sparse->appendSkew(integral,skew,symmetry))return;
        integral+=skew;integral.makeCompressed();
        for(Eigen::Index i=0;i<integral.nonZeros();++i)require(std::isfinite(integral.valuePtr()[i]),"nonfinite combined velocity block");
        const Matrix difference=integral-Matrix(integral.transpose());
        symmetry=ratio(difference.norm(),integral.norm());
    }
};
HydroSaddleVelocityBlock::HydroSaddleVelocityBlock(std::shared_ptr<const Data> data) : data_(std::move(data)) {}
HydroSaddleVelocityBlock HydroSaddleVelocityBlock::build(HydroDofLayout layout, const Matrix& integralVelocity,int tangentialCount)
{
    const auto started=WallClock::now();
    require(tangentialCount>=0,"negative tangential block size");
    require(integralVelocity.rows()==layout.velocityCount()+tangentialCount && integralVelocity.cols()==integralVelocity.rows(),
        "velocity block dimensions do not match layout");
    auto data=std::make_shared<Data>(std::move(layout),integralVelocity);
    data->tangential=tangentialCount;
    data->integral.makeCompressed();
    data->validateIntegral();
    data->validationSeconds=elapsedSeconds(started);
    return HydroSaddleVelocityBlock(std::move(data));
}
HydroSaddleVelocityBlock HydroSaddleVelocityBlock::buildWithTransportSplit(HydroDofLayout layout,
    const Matrix& integralVelocity,int tangentialCount)
{
    require(tangentialCount>=0,"negative tangential block size");
    require(integralVelocity.rows()==layout.velocityCount()+tangentialCount && integralVelocity.cols()==integralVelocity.rows(),
        "velocity block dimensions do not match layout");
    // Both original expressions consume the same unmodified A transpose.
    // The temporary transpose is released before the original validation path.
    auto split=[&] {
        const Matrix transpose(integralVelocity.transpose());
        Matrix symmetric=.5*(integralVelocity+transpose),skew=.5*(integralVelocity-transpose);
        return std::make_pair(std::move(symmetric),std::move(skew));
    }();
    const auto started=WallClock::now();
    auto data=std::make_shared<Data>(std::move(layout),std::move(split.first));
    data->tangential=tangentialCount;
    data->integral.makeCompressed();
    data->validateIntegral();
    data->validationSeconds=elapsedSeconds(started);
    const auto skewStarted=WallClock::now();
    Data::validateSkew(data->integral,split.second);
    data->appendSkew(split.second);
    data->validationSeconds+=elapsedSeconds(skewStarted);
    return HydroSaddleVelocityBlock(std::move(data));
}
HydroSaddleVelocityBlock HydroSaddleVelocityBlock::buildWithEnergyTransportSplit(const QuadtreeMesh& mesh,
    HydroDofLayout layout,const Matrix& actualVelocity,const HydroMomentumMass& mass,
    const HydroEnergyViscosity& viscosity,double dt)
{
    return buildEnergyTransportSplitImpl(mesh,std::move(layout),actualVelocity,mass,viscosity,dt,nullptr);
}
HydroSaddleVelocityBlock HydroSaddleVelocityBlock::buildWithEnergyTransportSplit(const QuadtreeMesh& mesh,
    HydroDofLayout layout,const Matrix& actualVelocity,const HydroMomentumMass& mass,
    const HydroEnergyViscosity& viscosity,double dt,const HydroMomentumGeometry& sparsePlanOwner)
{
    return buildEnergyTransportSplitImpl(mesh,std::move(layout),actualVelocity,mass,viscosity,dt,&sparsePlanOwner);
}
HydroSaddleVelocityBlock HydroSaddleVelocityBlock::buildFrozenActual(HydroDofLayout layout,Matrix&& actualVelocity,bool currentValuesFinite)
{
    PassiveAssemblySampler sampler;const auto started=sampler.mark();
    require(actualVelocity.rows()==layout.velocityCount()&&actualVelocity.cols()==actualVelocity.rows(),
        "velocity block dimensions do not match layout");
    auto data=std::make_shared<Data>(std::move(layout),actualVelocity,Data::SwapOwnedSplit{});
    data->integral.makeCompressed();
    if(!currentValuesFinite)for(Eigen::Index k=0;k<data->integral.nonZeros();++k)
        require(std::isfinite(data->integral.valuePtr()[k]),"nonfinite velocity block");
    // Actual nonsymmetric momentum A: no requested symmetric/SPD/split audit.
    // The Frozen caller validated this actual layout/source/full clock/live
    // mesh in the same callback-free invocation before constructing A.
    data->symmetry=std::numeric_limits<double>::quiet_NaN();
    data->validationSeconds=sampler.seconds(started);
    return HydroSaddleVelocityBlock(std::move(data));
}
HydroSaddleVelocityBlock HydroSaddleVelocityBlock::buildEnergyTransportSplitImpl(const QuadtreeMesh& mesh,
    HydroDofLayout layout,const Matrix& actualVelocity,const HydroMomentumMass& mass,
    const HydroEnergyViscosity& viscosity,double dt,const HydroMomentumGeometry* sparsePlanOwner)
{
    const bool sourceQualified=sparsePlanOwner&&Detail::sparseValidationFloatingQualified();
    require(std::isfinite(dt)&&dt>0,"invalid energy-certificate dt");
    layout.validateMesh(mesh);
    require(layout.sameIdentity(mass.layout())&&layout.sameIdentity(viscosity.strain().layout()),"different energy-certificate actual layout");
    require(&mass.geometry().metrics().divergence().matrix()==&viscosity.strain().metrics().divergence().matrix(),
        "different energy-certificate authoritative D_I");
    if(mass.model()==HydroPhaseGeometry::Model::GeometricVof)mass.validateVofFor(mesh,mass.vofIdentity(),mass.timeLayer());
    else mass.validateFor(mesh,mass.geometryIdentity(),mass.timeLayer());
    require(actualVelocity.rows()==layout.velocityCount()&&actualVelocity.cols()==actualVelocity.rows(),
        "velocity block dimensions do not match layout");
    if(sparsePlanOwner)require(sparsePlanOwner->sameIdentity(mass.geometry()),
        "different sparse-validation actual geometry owner");
    HydroVelocitySplitPlanReport splitReport;splitReport.requested=sparsePlanOwner!=nullptr;
    if(splitReport.requested)splitReport.inputNonZeros=static_cast<std::size_t>(actualVelocity.nonZeros());
    const auto splitCache=sparsePlanOwner?sparsePlanOwner->sparseValidationPlans():nullptr;
    auto splitPlan=Detail::findVelocitySplit(splitCache,actualVelocity,sourceQualified,splitReport);
    auto preparedSplit=splitPlan?splitPlan->freshMatrices():nullptr;
    const auto finishSplit=[&](Matrix& symmetric,const Matrix& skew,auto ownership) {
        if(splitReport.requested)splitReport.unionNonZeros=static_cast<std::size_t>(skew.nonZeros());
        const auto started=WallClock::now();
        auto data=[&] {
            if constexpr(std::is_same_v<decltype(ownership),Data::SwapOwnedSplit>)
                return std::make_shared<Data>(std::move(layout),symmetric,ownership);
            else return std::make_shared<Data>(std::move(layout),std::move(symmetric));
        }();
        data->integral.makeCompressed();data->splitPlan=splitReport;
        Detail::SparseValidationContext sparse(sparsePlanOwner?sparsePlanOwner->sparseValidationPlans():nullptr,
            data->integral,skew,sourceQualified,sparsePlanOwner!=nullptr,data->sparseValidation);
        auto* validation=sparsePlanOwner?&sparse:nullptr;
        const EnergySpdInputs energy{mass,viscosity,dt};data->validateIntegral(&energy,validation);
        data->validationSeconds=elapsedSeconds(started);const auto skewStarted=WallClock::now();
        Data::validateSkew(data->integral,skew,validation);data->appendSkew(skew,validation);
        data->validationSeconds+=elapsedSeconds(skewStarted);
        if(sparsePlanOwner)Detail::publishVelocitySplit(splitCache,actualVelocity,skew,sourceQualified,data->splitPlan);
        return HydroSaddleVelocityBlock(std::move(data));
    };
#if EIGEN_WORLD_VERSION==3 && EIGEN_MAJOR_VERSION==4 && EIGEN_MINOR_VERSION==1 && !defined(EIGEN_SCALAR_BINARY_OP_PLUGIN) && !defined(EIGEN_SPARSE_CREATE_TEMPORARY_PLUGIN) && !defined(EIGEN_SPARSE_TRANSPOSED_COPY_PLUGIN) && !defined(EIGEN_SPARSE_COMPRESSED_STORAGE_REALLOCATE_PLUGIN) && !defined(EIGEN_SPARSEMATRIX_PLUGIN)
    if(preparedSplit){
        // Preserve the original whole S expression, then the whole W expression.
        splitPlan->freshValues(actualVelocity,*preparedSplit);
        splitReport.symmetricValuesUsed=true;splitReport.skewValuesUsed=true;
        auto result=finishSplit(preparedSplit->first,preparedSplit->second,Data::SwapOwnedSplit{});
        return result;
    }
#endif
    auto split=[&] {
        if(preparedSplit){
            splitPlan->freshValues(actualVelocity,*preparedSplit);
            splitReport.symmetricValuesUsed=true;splitReport.skewValuesUsed=true;
            return std::move(*preparedSplit);
        }
        const Matrix transpose(actualVelocity.transpose());
        Matrix symmetric=.5*(actualVelocity+transpose),skew=.5*(actualVelocity-transpose);
        return std::make_pair(std::move(symmetric),std::move(skew));
    }();
    return finishSplit(split.first,split.second,std::false_type{});
}
const HydroVelocitySpdReport& HydroSaddleVelocityBlock::spdReport() const noexcept {return data_->spd;}
const HydroSparseValidationReport& HydroSaddleVelocityBlock::sparseValidationReport() const noexcept {return data_->sparseValidation;}
const HydroVelocitySplitPlanReport& HydroSaddleVelocityBlock::splitPlanReport() const noexcept {return data_->splitPlan;}
double HydroSaddleVelocityBlock::validationSeconds() const noexcept {return data_->validationSeconds;}
const HydroDofLayout& HydroSaddleVelocityBlock::layout() const noexcept {return data_->layout;}
const HydroSaddleVelocityBlock::Matrix& HydroSaddleVelocityBlock::matrix() const noexcept {return data_->integral;}
HydroSaddleVelocityBlock HydroSaddleVelocityBlock::withSkew(const Matrix& skew) const
{
    const auto started=WallClock::now();
    Data::validateSkew(matrix(),skew);
    auto data=std::make_shared<Data>(*data_);
    data->sparseValidation=HydroSparseValidationReport{};
    data->splitPlan=HydroVelocitySplitPlanReport{};
    data->appendSkew(skew);
    data->validationSeconds+=elapsedSeconds(started);
    return HydroSaddleVelocityBlock(std::move(data));
}
double HydroSaddleVelocityBlock::symmetryDefect() const noexcept {return data_->symmetry;}
int HydroSaddleVelocityBlock::tangentialCount() const noexcept{return data_->tangential;}

HydroSaddleRhs::HydroSaddleRhs(HydroDofLayout layout,const std::vector<double>& integralMomentum,
    const std::vector<double>& prescribedNormalVelocity,const std::vector<double>& tangentialMomentum)
    : layout_(std::move(layout)),momentum_(integralMomentum),prescribed_(prescribedNormalVelocity),tangential_(tangentialMomentum)
{
    require(momentum_.size()==static_cast<std::size_t>(layout_.velocityCount()) && prescribed_.size()==momentum_.size(),
        "RHS dimensions do not match layout");
    for (const auto* values : {&momentum_,&prescribed_,&tangential_})
        for (double value : *values) require(std::isfinite(value),"nonfinite RHS/prescribed velocity");
}
void HydroSaddleRhs::validateFor(const HydroDofLayout& expected,const QuadtreeMesh& mesh) const
{
    require(layout_.sameIdentity(expected),"different RHS layout identity");
    expected.validateMesh(mesh);
}
HydroSaddleResult::HydroSaddleResult(HydroNativeValues values,HydroSaddleReport report,std::vector<double> tangential)
    : values_(std::move(values)),report_(report),tangential_(std::move(tangential)) {}

struct HydroSinglePhaseSaddle::Data
{
    Data(HydroPressureLoad p,HydroSaddleVelocityBlock a,HydroSaddleConfig options)
        : pressure(std::move(p)),velocity(std::move(a)),config(options) {}
    HydroPressureLoad pressure;
    HydroSaddleVelocityBlock velocity;
    HydroSaddleConfig config;
    std::vector<int> freeByVelocity,velocityByFree,component,gaugeByComponent;
    std::vector<HydroSaddlePressureOutlet> outlets;
    std::vector<double> outletLoads;
    int gauges=0;
    Matrix system,quotient;
    bool quotientAliasesSystem=false;
    const Matrix& factorInput() const noexcept {return quotientAliasesSystem?system:quotient;}
    HydroSaddlePreparationTimings preparation;
    HydroSaddleAssemblyPlanReport assemblyPlan;
    bool assemblyTimingsAvailable=false;
    double partitionSeconds=0,systemAssemblySeconds=0,quotientAssemblySeconds=0;
    std::vector<int> quotientByPressure,pressureByQuotient;
    std::vector<double> volumes,componentVolumes;
    BlockQuotientFactor factor;
    // SuperLU const backsolve mutates its private B/X/status workspace.
    // Copies share this owner; only its native backsolve is serialized.
    mutable std::mutex factorMutex;
    double factorSeconds=0;

    // Invert the ORIGINAL mixed equations for arbitrary RHS/corrections.
    // Only closed components use a quotient representative and volume lift;
    // every outlet-anchored pressure coordinate stays in the same factor.
    Eigen::VectorXd solveAugmented(const Eigen::VectorXd& rhs,QuotientSolveTrace& trace,const HydroNativeValues* initialGuess=nullptr) const
    {
        require(rhs.size()==system.rows()&&rhs.allFinite(),"invalid augmented quotient RHS");
        const Eigen::Index nf=static_cast<Eigen::Index>(velocityByFree.size());
        const int np=static_cast<int>(component.size());
        Eigen::VectorXd result=Eigen::VectorXd::Zero(system.rows());
        std::vector<CompensatedSum> continuity(gauges);
        for(int p=0;p<np;++p) {
            const int group=gaugeByComponent[component[p]];
            if(group>=0)continuity[group].add(rhs[nf+p]);
        }
        for(int group=0;group<gauges;++group) {
            const double lambda=continuity[group].value()/componentVolumes[group];
            require(std::isfinite(lambda),"nonfinite quotient gauge multiplier");
            result[nf+np+group]=lambda;
        }
        Eigen::VectorXd reduced(factorInput().rows());
        reduced.head(nf)=rhs.head(nf);
        for(std::size_t q=0;q<pressureByQuotient.size();++q) {
            const int p=pressureByQuotient[q],group=gaugeByComponent[component[p]];
            reduced[nf+static_cast<Eigen::Index>(q)]=group<0?rhs[nf+p]:
                std::fma(-volumes[p],result[nf+np+group],rhs[nf+p]);
        }
        require(reduced.allFinite(),"nonfinite reduced quotient RHS");
        Eigen::VectorXd reducedGuess;
        const bool useInitialGuess=initialGuess&&!reduced.isZero(0);
        if(useInitialGuess){
            reducedGuess.resize(factorInput().rows());
            for(Eigen::Index i=0;i<nf;++i)reducedGuess[i]=initialGuess->normalVelocity()[velocityByFree[static_cast<std::size_t>(i)]];
            std::vector<double> representative(gauges,0);
            for(int p=0;p<np;++p){const int group=gaugeByComponent[component[p]];
                if(group>=0&&quotientByPressure[p]<0)representative[group]=initialGuess->pressure()[p];}
            for(std::size_t q=0;q<pressureByQuotient.size();++q){
                const int p=pressureByQuotient[q],group=gaugeByComponent[component[p]];
                reducedGuess[nf+static_cast<Eigen::Index>(q)]=group<0?initialGuess->pressure()[p]:
                    initialGuess->pressure()[p]-representative[group];
            }
            require(reducedGuess.allFinite(),"nonfinite reduced native initial guess");
        }
        // Isolated singleton pressure components can leave an empty factor.
        if(factorInput().rows()!=0) {
            const Eigen::VectorXd solved=[&]{
                std::lock_guard<std::mutex> guard(factorMutex);
                const Eigen::VectorXd candidate=factor.solve(reduced,trace,useInitialGuess?&reducedGuess:nullptr);
                require(factor.info()==Eigen::Success,"quotient saddle solve failed");
                return candidate;
            }();
            require(solved.allFinite(),"nonfinite quotient saddle solution");
            result.head(nf)=solved.head(nf);
            for(std::size_t q=0;q<pressureByQuotient.size();++q)
                result[nf+pressureByQuotient[q]]=solved[nf+static_cast<Eigen::Index>(q)];
        }
        std::vector<CompensatedSum> meanRhs(gauges);
        for(int group=0;group<gauges;++group)meanRhs[group].add(rhs[nf+np+group]);
        for(int p=0;p<np;++p) {
            const int group=gaugeByComponent[component[p]];if(group<0)continue;
            const double product=-volumes[p]*result[nf+p];
            meanRhs[group].add(product);
            meanRhs[group].add(std::fma(-volumes[p],result[nf+p],-product));
        }
        std::vector<double> shift(gauges);
        for(int group=0;group<gauges;++group) {
            shift[group]=meanRhs[group].value()/componentVolumes[group];
            require(std::isfinite(shift[group]),"nonfinite quotient pressure mean lift");
        }
        for(int p=0;p<np;++p) {
            const int group=gaugeByComponent[component[p]];if(group>=0)result[nf+p]+=shift[group];
        }
        require(result.allFinite(),"nonfinite lifted augmented solution");
        return result;
    }
};
HydroSinglePhaseSaddle::HydroSinglePhaseSaddle(std::shared_ptr<const Data> data) : data_(std::move(data)) {}
HydroSinglePhaseSaddle HydroSinglePhaseSaddle::buildClosed(HydroPressureLoad pressure,
    HydroSaddleVelocityBlock velocity,HydroSaddleConfig config,const HydroSaddlePreparationTimings& preparation)
{return buildImpl(std::move(pressure),std::move(velocity),nullptr,config,preparation);}
HydroSinglePhaseSaddle HydroSinglePhaseSaddle::buildWithPressureOutlets(HydroPressureLoad pressure,
    HydroSaddleVelocityBlock velocity,const std::vector<HydroSaddlePressureOutlet>& outlets,HydroSaddleConfig config,
    const HydroSaddlePreparationTimings& preparation)
{
    require(!outlets.empty(),"pressure outlet list is empty");
    return buildImpl(std::move(pressure),std::move(velocity),&outlets,config,preparation);
}
HydroSinglePhaseSaddle HydroSinglePhaseSaddle::buildClosed(HydroPressureLoad pressure,
    HydroSaddleVelocityBlock velocity,HydroSaddleConfig config,const HydroSaddlePreparationTimings& preparation,
    const HydroMomentumGeometry& sparsePlanOwner)
{return buildImpl(std::move(pressure),std::move(velocity),nullptr,config,preparation,&sparsePlanOwner);}
HydroSinglePhaseSaddle HydroSinglePhaseSaddle::buildWithPressureOutlets(HydroPressureLoad pressure,
    HydroSaddleVelocityBlock velocity,const std::vector<HydroSaddlePressureOutlet>& outlets,HydroSaddleConfig config,
    const HydroSaddlePreparationTimings& preparation,const HydroMomentumGeometry& sparsePlanOwner)
{
    require(!outlets.empty(),"pressure outlet list is empty");
    return buildImpl(std::move(pressure),std::move(velocity),&outlets,config,preparation,&sparsePlanOwner);
}
HydroSinglePhaseSaddle HydroSinglePhaseSaddle::buildImpl(HydroPressureLoad pressure,
    HydroSaddleVelocityBlock velocity,const std::vector<HydroSaddlePressureOutlet>* outlets,HydroSaddleConfig config,
    const HydroSaddlePreparationTimings& preparation,const HydroMomentumGeometry* sparsePlanOwner)
{
    const bool sourceQualified=sparsePlanOwner&&Detail::saddleScatterQualified();
    PassiveAssemblySampler sampler;const auto partitionStarted=sampler.mark();
    require(pressure.layout().sameIdentity(velocity.layout()),"different velocity block layout identity");
    require(!config.gmresGlobalPressure||config.gmresOnlyParallel,
        "global pressure preconditioner requires explicit GMRES-only mode");
    for (double tolerance : {config.algebraicRelativeTolerance,config.momentumRelativeTolerance,
        config.continuityRelativeTolerance,config.pressureMeanRelativeTolerance,config.boundaryFluxRelativeTolerance})
        require(std::isfinite(tolerance) && tolerance>0 && tolerance<1,"invalid residual tolerance");
    if(sparsePlanOwner){
        require(sparsePlanOwner->layout().sameIdentity(pressure.layout()),"different sparse-assembly actual layout");
        require(&sparsePlanOwner->metrics().divergence().matrix()==&pressure.divergence().matrix(),
            "different sparse-assembly authoritative D_I");
    }
    auto data=std::make_shared<Data>(std::move(pressure),std::move(velocity),config);
    data->preparation=preparation;
    const auto& c=data->pressure.matrix(); const int nu=static_cast<int>(c.rows()),np=static_cast<int>(c.cols());
    const int nm=data->velocity.tangentialCount();data->freeByVelocity.assign(nu+nm,-1);
    data->outletLoads.assign(nu,0);
    std::vector<int> outletByVelocity(nu,-1);std::vector<bool> anchoredPressure(np,false);
    if(outlets) {
        data->outlets=*outlets;
        for(std::size_t k=0;k<outlets->size();++k) {
            const auto& outlet=(*outlets)[k];const int u=outlet.velocityDof;
            require(u>=0&&u<nu,"pressure outlet velocity index out of range");
            require(std::isfinite(outlet.pressure),"nonfinite outlet pressure");
            require(outletByVelocity[u]<0,"duplicate pressure outlet");outletByVelocity[u]=static_cast<int>(k);
        }
    }
    std::vector<int> parents(np);std::iota(parents.begin(),parents.end(),0);
    for (int u=0;u<nu;++u)
    {
        int first=-1,last=-1,count=0;double firstValue=0,lastValue=0;
        for (HydroPressureLoad::Matrix::InnerIterator it(c,u);it;++it) {
            if (!count) {first=static_cast<int>(it.col());firstValue=it.value();}
            last=static_cast<int>(it.col());lastValue=it.value();++count;
        }
        require(count<=2,"invalid native pressure incidence");
        if(outletByVelocity[u]>=0) {
            require(count==1&&std::isfinite(firstValue)&&firstValue!=0,
                "pressure outlet requires a nonzero physical boundary incidence");
            const auto& outlet=data->outlets[outletByVelocity[u]];
            const double load=firstValue*outlet.pressure;
            require(std::isfinite(load)&&(outlet.pressure==0||load!=0),"unrepresentable outlet traction load");
            data->outletLoads[u]=load;anchoredPressure[first]=true;
            data->freeByVelocity[u]=static_cast<int>(data->velocityByFree.size());data->velocityByFree.push_back(u);
        }
        else if (count==2) {
            require(std::isfinite(firstValue)&&firstValue!=0&&lastValue==-firstValue,
                "free pressure incidence does not annihilate component constants");
            data->freeByVelocity[u]=static_cast<int>(data->velocityByFree.size());data->velocityByFree.push_back(u);
            parents[root(parents,last)]=root(parents,first);
        }
    }
    for(int m=0;m<nm;++m){const int u=nu+m;
        data->freeByVelocity[u]=static_cast<int>(data->velocityByFree.size());data->velocityByFree.push_back(u);}
    data->component.resize(np);std::vector<int> numbered(np,-1);int components=0;
    for (int p=0;p<np;++p) {
        const int r=root(parents,p);if (numbered[r]<0) numbered[r]=components++;
        data->component[p]=numbered[r];
    }
    std::vector<bool> anchoredComponent(components,false);
    for(int p=0;p<np;++p)if(anchoredPressure[p])anchoredComponent[data->component[p]]=true;
    data->gaugeByComponent.assign(components,-1);
    for(int group=0;group<components;++group)if(!anchoredComponent[group])data->gaugeByComponent[group]=data->gauges++;
    const Eigen::Index nf=static_cast<Eigen::Index>(data->velocityByFree.size()),size=nf+np+data->gauges;
    data->partitionSeconds=sampler.seconds(partitionStarted);const auto systemStarted=sampler.mark();
    const auto& actualVelocity=data->velocity.matrix();
    const Detail::SaddleScatterBinding binding{actualVelocity,c,data->freeByVelocity,data->velocityByFree,
        data->component,data->gaugeByComponent,data->outlets,nu,np,nm,data->gauges};
    auto scatterCache=sparsePlanOwner?sparsePlanOwner->sparseValidationPlans():nullptr;
    data->assemblyPlan.requested=sparsePlanOwner!=nullptr;
    const auto scatterPlan=Detail::prepareSaddleScatter(scatterCache,binding,sourceQualified,data->assemblyPlan);
    auto freshSystem=scatterPlan?Detail::saddleFreshMatrix(scatterPlan->systemRows,
        scatterPlan->systemOuter,scatterPlan->systemInner):std::unique_ptr<Matrix>{};
    const bool systemPrepared=scatterPlan&&freshSystem;
    if(systemPrepared){
        const auto& currentVolumes=data->pressure.layout().domain().cellVolumes();
        for(int p=0;p<np;++p)require(std::isfinite(currentVolumes[p])&&currentVolumes[p]>0,"invalid pressure mean volume");
        Detail::MatrixStageExecution execution;
        Detail::matrixStage(scatterPlan->systemRows,[&](Eigen::Index col) noexcept {
            for(auto slot=scatterPlan->systemOuter[static_cast<std::size_t>(col)];slot<scatterPlan->systemOuter[static_cast<std::size_t>(col)+1];++slot){
                const auto k=static_cast<std::size_t>(slot);const auto source=scatterPlan->systemSources[k];
                const auto kind=scatterPlan->systemKinds[k];
                const double* value=kind==0?actualVelocity.valuePtr()+source:kind==1?c.valuePtr()+source:currentVolumes.data()+source;
                std::memcpy(freshSystem->valuePtr()+slot,value,sizeof(double));
            }
        },true,Detail::HydroSaddleScatterPlan::byteLimit-scatterPlan->retainedBytes,execution);
        data->assemblyPlan.scatterParallelUsed=execution.parallelUsed;
        data->assemblyPlan.scatterParallelFallback=execution.parallelFallback;
        data->assemblyPlan.scatterMaximumWorkers=execution.maximumWorkers;
        data->system.swap(*freshSystem);data->system.makeCompressed();data->assemblyPlan.systemScatterUsed=true;
    }else{
        std::vector<Eigen::Triplet<double,Eigen::Index>> entries;
        const auto& a=data->velocity.matrix();
        entries.reserve(static_cast<std::size_t>(a.nonZeros()+2*c.nonZeros()+2*np));
        for (int u : data->velocityByFree)
            for (HydroSaddleVelocityBlock::Matrix::InnerIterator it(a,u);it;++it)
                if (data->freeByVelocity[it.col()]>=0) entries.emplace_back(data->freeByVelocity[u],data->freeByVelocity[it.col()],it.value());
        for (int u : data->velocityByFree)if(u<nu)
            for (HydroPressureLoad::Matrix::InnerIterator it(c,u);it;++it) {
                entries.emplace_back(data->freeByVelocity[u],nf+it.col(),it.value());
                entries.emplace_back(nf+it.col(),data->freeByVelocity[u],it.value());
            }
        const auto& volumes=data->pressure.layout().domain().cellVolumes();
        for (int p=0;p<np;++p) {
            require(std::isfinite(volumes[p]) && volumes[p]>0,"invalid pressure mean volume");
            const int group=data->gaugeByComponent[data->component[p]];if(group<0)continue;
            entries.emplace_back(nf+p,nf+np+group,volumes[p]);
            entries.emplace_back(nf+np+group,nf+p,volumes[p]);
        }
        data->system.resize(size,size);data->system.setFromTriplets(entries.begin(),entries.end());data->system.makeCompressed();
    }
    freshSystem.reset();
    const auto& volumes=data->pressure.layout().domain().cellVolumes();
    data->systemAssemblySeconds=sampler.seconds(systemStarted);const auto quotientStarted=sampler.mark();
    // Keep the original mixed K above. Remove one redundant pressure/continuity
    // coordinate only per closed component from the INTERNAL factor input;
    // anchored components keep all coordinates, without a volume-gauge border.
    data->volumes=volumes;
    data->componentVolumes.resize(data->gauges);
    std::vector<CompensatedSum> volumeSums(data->gauges);
    std::vector<int> representative(data->gauges,-1);
    data->quotientByPressure.assign(np,-1);
    for(int p=0;p<np;++p) {
        const int group=data->gaugeByComponent[data->component[p]];
        if(group>=0) {
            volumeSums[group].add(volumes[p]);
            if(representative[group]<0){representative[group]=p;continue;}
        }
        data->quotientByPressure[p]=static_cast<int>(data->pressureByQuotient.size());
        data->pressureByQuotient.push_back(p);
    }
    for(int group=0;group<data->gauges;++group) {
        data->componentVolumes[group]=volumeSums[group].value();
        require(data->componentVolumes[group]>0,"invalid quotient component volume");
    }
    const auto reducedIndex=[&](Eigen::Index index)->Eigen::Index {
        if(index<nf)return index;
        if(index<nf+np) {
            const int q=data->quotientByPressure[static_cast<std::size_t>(index-nf)];
            return q<0?-1:nf+q;
        }
        return -1;
    };
    bool identityQuotient=systemPrepared&&scatterPlan->identityQuotient
        &&data->gauges==0&&scatterPlan->matchesQuotient(data->quotientByPressure,data->pressureByQuotient)
        &&data->system.rows()==nf+np&&data->pressureByQuotient.size()==static_cast<std::size_t>(np);
    if(identityQuotient)for(int p=0;p<np;++p)
        identityQuotient=identityQuotient&&data->quotientByPressure[p]==p&&data->pressureByQuotient[p]==p;
    data->quotientAliasesSystem=identityQuotient;data->assemblyPlan.quotientAliasesSystem=identityQuotient;
    auto freshQuotient=!identityQuotient&&scatterPlan&&!scatterPlan->identityQuotient&&scatterPlan->matchesQuotient(data->quotientByPressure,data->pressureByQuotient)
        ?Detail::saddleFreshMatrix(scatterPlan->quotientRows,scatterPlan->quotientOuter,scatterPlan->quotientInner)
        :std::unique_ptr<Matrix>{};
    const bool quotientPrepared=!identityQuotient&&scatterPlan&&freshQuotient;
    if(identityQuotient){
        // No pressure coordinate or gauge was removed: the current factor input
        // is the same immutable system, including every explicit signed zero.
    }else if(quotientPrepared){
        Detail::MatrixStageExecution execution;
        Detail::matrixStage(scatterPlan->quotientRows,[&](Eigen::Index col) noexcept {
            for(auto k=scatterPlan->quotientOuter[static_cast<std::size_t>(col)];k<scatterPlan->quotientOuter[static_cast<std::size_t>(col)+1];++k)
                std::memcpy(freshQuotient->valuePtr()+k,data->system.valuePtr()+scatterPlan->quotientSources[static_cast<std::size_t>(k)],sizeof(double));
        },true,Detail::HydroSaddleScatterPlan::byteLimit-scatterPlan->retainedBytes,execution);
        data->assemblyPlan.scatterParallelUsed=data->assemblyPlan.scatterParallelUsed||execution.parallelUsed;
        data->assemblyPlan.scatterParallelFallback=data->assemblyPlan.scatterParallelFallback||execution.parallelFallback;
        data->assemblyPlan.scatterMaximumWorkers=std::max(data->assemblyPlan.scatterMaximumWorkers,execution.maximumWorkers);
        data->quotient.swap(*freshQuotient);data->quotient.makeCompressed();data->assemblyPlan.quotientScatterUsed=true;
    }else{
        std::vector<Eigen::Triplet<double,Eigen::Index>> reducedEntries;
        reducedEntries.reserve(static_cast<std::size_t>(data->system.nonZeros()));
        for(Eigen::Index col=0;col<data->system.outerSize();++col) {
            const Eigen::Index reducedCol=reducedIndex(col);if(reducedCol<0)continue;
            for(Matrix::InnerIterator it(data->system,col);it;++it) {
                const Eigen::Index reducedRow=reducedIndex(it.row());
                if(reducedRow>=0)reducedEntries.emplace_back(reducedRow,reducedCol,it.value());
            }
        }
        const Eigen::Index factorSize=nf+static_cast<Eigen::Index>(data->pressureByQuotient.size());
        data->quotient.resize(factorSize,factorSize);
        data->quotient.setFromTriplets(reducedEntries.begin(),reducedEntries.end());data->quotient.makeCompressed();
    }
    freshQuotient.reset();
    const Eigen::Index factorSize=nf+static_cast<Eigen::Index>(data->pressureByQuotient.size());
    if(sparsePlanOwner){
        data->assemblyPlan.systemNonZeros=static_cast<std::size_t>(data->system.nonZeros());
        data->assemblyPlan.quotientNonZeros=static_cast<std::size_t>(data->factorInput().nonZeros());
        if(systemPrepared&&(quotientPrepared||data->quotientAliasesSystem))
            Detail::publishSaddleScatter(scatterCache,scatterPlan);
    }
    data->quotientAssemblySeconds=sampler.seconds(quotientStarted);data->assemblyTimingsAvailable=sampler.available;
    if(factorSize!=0) {
        const auto started=WallClock::now();
        data->factor.compute(data->factorInput(),nf,config.gmresOnlyParallel,config.gmresGlobalPressure);
        if(data->factor.info()!=Eigen::Success)
            throw std::runtime_error("[HydroSinglePhaseSaddle] quotient saddle preconditioner/factor preparation failed");
        data->factorSeconds=elapsedSeconds(started);
    }
    return HydroSinglePhaseSaddle(std::move(data));
}
const HydroPressureLoad& HydroSinglePhaseSaddle::pressureLoad() const noexcept {return data_->pressure;}
const HydroSaddleVelocityBlock& HydroSinglePhaseSaddle::velocityBlock() const noexcept {return data_->velocity;}
const HydroDofLayout& HydroSinglePhaseSaddle::layout() const noexcept {return data_->pressure.layout();}
const HydroSinglePhaseSaddle::Matrix& HydroSinglePhaseSaddle::matrix() const noexcept {return data_->system;}
const std::vector<int>& HydroSinglePhaseSaddle::freeVelocityDofs() const noexcept {return data_->velocityByFree;}
const std::vector<int>& HydroSinglePhaseSaddle::pressureComponents() const noexcept {return data_->component;}
const std::vector<int>& HydroSinglePhaseSaddle::pressureGaugeByComponent() const noexcept {return data_->gaugeByComponent;}
const std::vector<HydroSaddlePressureOutlet>& HydroSinglePhaseSaddle::pressureOutlets() const noexcept {return data_->outlets;}
const std::vector<double>& HydroSinglePhaseSaddle::pressureOutletLoads() const noexcept {return data_->outletLoads;}
int HydroSinglePhaseSaddle::pressureGaugeCount() const noexcept {return data_->gauges;}
HydroSaddleResult HydroSinglePhaseSaddle::solve(const QuadtreeMesh& mesh,const HydroSaddleRhs& rhs) const
{return solveImpl(mesh,rhs,nullptr);}
HydroSaddleResult HydroSinglePhaseSaddle::solve(const QuadtreeMesh& mesh,const HydroSaddleRhs& rhs,
    const HydroNativeValues& initialGuess) const
{return solveImpl(mesh,rhs,&initialGuess);}
HydroSaddleResult HydroSinglePhaseSaddle::solveImpl(const QuadtreeMesh& mesh,const HydroSaddleRhs& rhs,
    const HydroNativeValues* initialGuess) const
{
    const auto started=WallClock::now();
    PassiveAssemblySampler rhsSampler;const auto rhsStarted=rhsSampler.mark();
    rhs.validateFor(layout(),mesh);
    const int nu=layout().velocityCount(),nm=data_->velocity.tangentialCount(),nv=nu+nm,
        np=layout().pressureCount(),nf=static_cast<int>(data_->velocityByFree.size());
    const auto& c=data_->pressure.matrix(); const auto& a=data_->velocity.matrix();
    if(initialGuess){
        require(data_->config.gmresOnlyParallel,"native initial guess requires explicit GMRES-only mode");
        require(nm==0,"native initial guess does not support tangential unknowns");
        initialGuess->validateFor(layout(),mesh);
        for(const auto* values:{&initialGuess->normalVelocity(),&initialGuess->pressure()})
            for(double value:*values)require(std::isfinite(value),"nonfinite native initial guess");
    }
    require(rhs.tangentialMomentum().size()==static_cast<std::size_t>(nm),"tangential RHS dimensions do not match velocity block");
    Eigen::VectorXd known=Eigen::VectorXd::Zero(nv),momentum(nv);
    for(int u=0;u<nu;++u){known[u]=rhs.prescribedNormalVelocity()[u];momentum[u]=rhs.integralMomentum()[u];}
    if(!data_->outlets.empty())for(int u=0;u<nu;++u) {
        momentum[u]+=data_->outletLoads[u];require(std::isfinite(momentum[u]),"nonfinite traction-inclusive momentum RHS");
    }
    for(int m=0;m<nm;++m)momentum[nu+m]=rhs.tangentialMomentum()[m];
    std::vector<double> net(data_->gauges,0),gross(data_->gauges,0);HydroSaddleReport report;
    report.freeVelocityCount=nf-nm;report.freeTangentialCount=nm;report.pressureGaugeCount=data_->gauges;
    report.factorDimension=static_cast<int>(data_->factorInput().rows());report.factorNonZeros=data_->factorInput().nonZeros();
    // Actual backend flags are filled only after the original gates pass.
    report.velocityBlockSeconds=data_->velocity.validationSeconds();report.factorSeconds=data_->factorSeconds;
    report.velocitySparseValidation=data_->velocity.sparseValidationReport();
    report.velocitySplitPlan=data_->velocity.splitPlanReport();
    report.sparseAssemblyPlan=data_->assemblyPlan;
    const auto& spd=data_->velocity.spdReport();report.velocitySpdTier=static_cast<int>(spd.tier);
    report.velocitySpdGramAttempted=spd.gramAttempted;report.velocitySpdGramBoundsAvailable=spd.gramBoundsAvailable;
    report.velocitySpdGramBoundUpper=spd.gramBoundUpper;report.velocitySpdScaledMassLower=spd.scaledMassLower;
    report.velocitySpdCertificateSeconds=spd.certificateSeconds;report.velocitySpdLltSeconds=spd.lltSeconds;
    const auto& preparation=data_->preparation;report.frozenPreparationTimingsAvailable=preparation.available;
    report.frozenInputValidationSeconds=preparation.inputValidationSeconds;report.frozenForceValidationSeconds=preparation.forceValidationSeconds;
    report.frozenMatrixActionSeconds=preparation.matrixActionSeconds;report.frozenRhsSeconds=preparation.rhsSeconds;
    report.frozenVelocityPrepareSeconds=preparation.velocityPrepareSeconds;
    report.saddlePartitionSeconds=data_->partitionSeconds;report.saddleSystemAssemblySeconds=data_->systemAssemblySeconds;
    report.saddleQuotientAssemblySeconds=data_->quotientAssemblySeconds;
    for (int u=0;u<nu;++u) {
        if (data_->freeByVelocity[u]>=0) require(known[u]==0,"prescribed velocity payload on a free face");
        else {
            int count=0;
            for (HydroPressureLoad::Matrix::InnerIterator it(c,u);it;++it) {
                const int group=data_->gaugeByComponent[data_->component[it.col()]];
                const double flow=-it.value()*known[u];
                if(group>=0){net[group]+=flow;gross[group]+=std::abs(flow);}++count;
            }
            if (!count) require(known[u]==0,"nonzero RZ axis velocity");
        }
    }
    for (int component=0;component<data_->gauges;++component) {
        const double defect=ratio(net[component],gross[component]);
        report.maxBoundaryFluxDefect=std::max(report.maxBoundaryFluxDefect,defect);
        accept(defect,data_->config.boundaryFluxRelativeTolerance,"component boundary flux compatibility");
    }
    const Eigen::VectorXd knownMomentum=a*known,knownContinuity=c.transpose()*known.head(nu);
    Eigen::VectorXd b=Eigen::VectorXd::Zero(data_->system.rows());
    for (int i=0;i<nf;++i) b[i]=momentum[data_->velocityByFree[i]]-knownMomentum[data_->velocityByFree[i]];
    b.segment(nf,np)=-knownContinuity;
    if (!b.allFinite()) throw std::runtime_error("[HydroSinglePhaseSaddle] nonfinite eliminated RHS");
    report.saddleRhsEliminationSeconds=rhsSampler.seconds(rhsStarted);
    report.saddleAssemblyTimingsAvailable=data_->assemblyTimingsAvailable&&rhsSampler.available;
    // Assess the original equations and tolerances before publishing a report.
    // No geometry, material or manufactured equilibrium enters this solve.
    struct Assessment {
        std::vector<double> u,p;
        std::array<double,4> residual{};
        double gaugeMultiplier=0,boundaryFlux=0;
    };
    const std::array<double,4> tolerance{{data_->config.algebraicRelativeTolerance,
        data_->config.momentumRelativeTolerance,data_->config.continuityRelativeTolerance,
        data_->config.pressureMeanRelativeTolerance}};
    const auto assess=[&](const Eigen::VectorXd& x){
        if(!x.allFinite())throw std::runtime_error("[HydroSinglePhaseSaddle] nonfinite reference saddle solution");
        Assessment value;value.u=rhs.prescribedNormalVelocity();value.u.resize(nv);value.p.resize(np);
        value.residual[0]=ratio((data_->system*x-b).stableNorm(),b.stableNorm());
        for(int i=0;i<nf;++i)value.u[data_->velocityByFree[i]]=x[i];
        for(int i=0;i<np;++i)value.p[i]=x[nf+i];
        Eigen::Map<const Eigen::VectorXd> velocity(value.u.data(),nv),pressure(value.p.data(),np);
        const Eigen::VectorXd inertia=a*velocity;
        Eigen::VectorXd load=Eigen::VectorXd::Zero(nv);load.head(nu)=c*pressure;
        double error=0,scale=0;
        for(int i:data_->velocityByFree){
            error=std::hypot(error,inertia[i]+load[i]-momentum[i]);
            scale=std::hypot(scale,std::abs(inertia[i])+std::abs(load[i])+std::abs(momentum[i]));
        }
        value.residual[1]=ratio(error,scale);
        const auto& d=data_->pressure.divergence().matrix();const Eigen::VectorXd continuity=d*velocity.head(nu);
        for(int i=0;i<np;++i){
            double flowScale=0;
            for(HydroIntegralDivergence::Matrix::InnerIterator it(d,i);it;++it)flowScale+=std::abs(it.value()*velocity[it.col()]);
            value.residual[2]=std::max(value.residual[2],ratio(continuity[i],flowScale));
        }
        std::vector<double> mean(data_->gauges,0),meanScale(data_->gauges,0);
        const auto& volumes=layout().domain().cellVolumes();
        for(int i=0;i<np;++i){const int group=data_->gaugeByComponent[data_->component[i]];if(group<0)continue;
            const double term=volumes[i]*value.p[i];mean[group]+=term;meanScale[group]+=std::abs(term);}
        for(int i=0;i<data_->gauges;++i){
            value.residual[3]=std::max(value.residual[3],ratio(mean[i],meanScale[i]));
            value.gaugeMultiplier=std::max(value.gaugeMultiplier,std::abs(x[nf+np+i]));
        }
        if(!data_->outlets.empty()) {
            // Assess actual fixed inlet plus solved outlet flux, with original
            // signed native incidence and the original boundary-flux budget.
            const auto components=data_->gaugeByComponent.size();
            std::vector<double> boundaryNet(components,0),boundaryGross(components,0);
            for(int u=0;u<nu;++u) {
                int count=0,p=-1;double coefficient=0;
                for(HydroPressureLoad::Matrix::InnerIterator it(c,u);it;++it){++count;p=static_cast<int>(it.col());coefficient=it.value();}
                if(count!=1)continue;const int group=data_->component[p];
                if(data_->gaugeByComponent[group]>=0)continue;
                const double flow=-coefficient*value.u[u];
                boundaryNet[group]+=flow;boundaryGross[group]+=std::abs(flow);
            }
            for(std::size_t group=0;group<components;++group)if(data_->gaugeByComponent[group]<0)
                value.boundaryFlux=std::max(value.boundaryFlux,ratio(boundaryNet[group],boundaryGross[group]));
        }
        return value;
    };
    QuotientSolveTrace trace;
    for(;;){
      bool directBaseline=false;
      try {
    const int directBefore=trace.directSolves;
    Eigen::VectorXd x=data_->solveAugmented(b,trace,initialGuess);
    directBaseline=trace.directSolves>directBefore;
    for(int iteration=0;;++iteration){
        PassiveAssemblySampler assessmentTimings;assessmentTimings.requireMaskedExceptions();
        const auto assessmentStarted=assessmentTimings.mark();const auto value=assess(x);
        assessmentTimings.addSeconds(assessmentStarted,trace.physicalAssessmentSeconds);
        trace.solveBreakdownTimingsAvailable=trace.solveBreakdownTimingsAvailable&&assessmentTimings.available;
        bool accepted=true;
        for(int k=0;k<4;++k)accepted=accepted && value.residual[k]<=tolerance[k];
        accepted=accepted && value.boundaryFlux<=data_->config.boundaryFluxRelativeTolerance;
        if(accepted){
            report.algebraicRelativeResidual=value.residual[0];report.freeMomentumDefect=value.residual[1];
            report.maxCellContinuityDefect=value.residual[2];report.maxPressureMeanDefect=value.residual[3];
            report.maxAbsGaugeMultiplier=value.gaugeMultiplier;report.linearRefinements=iteration;
            if(!data_->outlets.empty())report.maxBoundaryFluxDefect=std::max(report.maxBoundaryFluxDefect,value.boundaryFlux);
            report.usedSuperLU=trace.directSolves!=0;report.usedBlockGmres=trace.gmresSolves!=0;
            report.usedDirectFallback=trace.fallbackReason!=0;report.directFallbackReason=trace.fallbackReason;
            report.krylovIterations=trace.iterations;report.krylovRestarts=trace.restarts;report.krylovSolves=trace.gmresSolves;
            report.krylovRelativeResidual=trace.lastRelativeResidual;
            report.krylovScaledRelativeResidual=trace.lastScaledRelativeResidual;report.krylovStagnated=trace.stagnated;
            report.solveBreakdownTimingsAvailable=trace.solveBreakdownTimingsAvailable&&trace.gmresSolves>0&&trace.fallbackReason==0;
            report.krylovPreconditionApplySeconds=trace.krylovPreconditionApplySeconds;report.krylovSpmvSeconds=trace.krylovSpmvSeconds;
            report.krylovOrthogonalizationSeconds=trace.krylovOrthogonalizationSeconds;
            report.krylovTrueResidualSeconds=trace.krylovTrueResidualSeconds;report.physicalAssessmentSeconds=trace.physicalAssessmentSeconds;
            report.orthogonalizationParallelCalls=trace.orthogonalizationParallelCalls;
            report.orthogonalizationMaximumWorkers=trace.orthogonalizationMaximumWorkers;
            const auto& global=data_->factor.pressureGlobalReport();
            report.pressureGlobalRequested=data_->config.gmresOnlyParallel&&data_->config.gmresGlobalPressure;
            report.pressureGlobalAvailable=global.available;
            report.pressureGlobalPrepared=global.prepared;
            report.pressureGlobalTimingsAvailable=global.timingsAvailable&&trace.pressureGlobalTimingsAvailable;
            report.pressureGlobalRows=global.rows;report.pressureGlobalNonZeros=global.nonZeros;
            report.pressureGlobalFactorNonZeros=global.factorNonZeros;report.pressureGlobalStorageBytes=global.storageBytes;
            report.pressureGlobalAdmissionBytes=global.admissionBytes;report.pressureGlobalByteLimit=global.byteLimit;
            report.pressureGlobalSymbolicSeconds=global.symbolicSeconds;report.pressureGlobalNumericSeconds=global.numericSeconds;
            report.pressureGlobalSolveCalls=trace.pressureGlobalSolveCalls;report.pressureGlobalSolveSeconds=trace.pressureGlobalSolveSeconds;
            report.parallelPressureBlockCount=data_->factor.parallelPressureBlockCount();
            report.parallelPressureMaximumWorkers=data_->factor.parallelPressureMaximumWorkers();
            const auto& coarse=data_->factor.pressureCoarseReport();
            report.pressureCoarsePrepared=data_->factor.pressureCoarsePrepared();
            report.pressureCoarseDimension=coarse.dimension;report.pressureCoarseReason=static_cast<int>(coarse.reason);
            report.pressureCoarseApplications=trace.pressureCoarseApplications;
            report.pressureCoarseConditionUpper=coarse.conditionUpper;
            report.pressureCoarseFormationNormUpper=coarse.formationNormUpper;
            report.orderedSpmvRequested=data_->factor.orderedSpmvRequested();
            report.orderedSpmvEligible=data_->factor.orderedSpmvEligible();
            report.orderedSpmvCalls=trace.orderedSpmvCalls;
            report.orderedSpmvParallelCalls=trace.orderedSpmvParallelCalls;
            report.orderedSpmvUsed=trace.orderedSpmvParallelCalls!=0;
            report.orderedSpmvFallback=trace.orderedSpmvFallback;
            report.orderedSpmvMaximumWorkers=trace.orderedSpmvMaximumWorkers;
            report.orderedSpmvViewBytes=data_->factor.orderedSpmvViewBytes();
            report.directInputDeferred=data_->factor.directInputDeferred();
            report.directInputMaterialized=data_->factor.directInputMaterialized();
            report.preconditionerSeconds=data_->factor.preconditionerSeconds();
            report.directFactorSeconds=data_->factor.preparedDirectSeconds()+trace.lazyFactorSeconds;
            report.factorSeconds=data_->factorSeconds+trace.lazyFactorSeconds;
            report.solveSeconds=elapsedSeconds(started);
            const std::vector<double> normal(value.u.begin(),value.u.begin()+nu),tangential(value.u.begin()+nu,value.u.end());
            return HydroSaddleResult(HydroNativeValues(layout(),normal,value.p),report,tangential);
        }
        if(iteration==8){
            const std::array<const char*,4> stages{{"complete algebraic","free momentum","native cell continuity","component pressure mean"}};
            for(int k=0;k<4;++k)accept(value.residual[k],tolerance[k],stages[k]);
            accept(value.boundaryFlux,data_->config.boundaryFluxRelativeTolerance,"actual outlet component boundary flux");
        }
        Eigen::VectorXd residual=correctionResidual(data_->system,x,b);
        // An already accurate force balance is preserved by A*du+C*dp=0.
        // Reinjecting its roundoff can obscure relative continuity at u~0.
        // Every correction still solves the SAME full augmented matrix.
        if(value.residual[0]<=tolerance[0] && value.residual[1]<=tolerance[1])residual.head(nf).setZero();
        if(value.residual[3]<=tolerance[3])residual.tail(data_->gauges).setZero();
        x+=data_->solveAugmented(residual,trace);
    }
      } catch(const std::runtime_error&) {
        if(data_->config.gmresOnlyParallel)throw;
        // If Krylov passed its own norm but the original local gates or
        // eight compensated corrections fail, retry the ENTIRE original
        // assessment/refinement from b with the original lazy direct factor.
        std::lock_guard<std::mutex> guard(data_->factorMutex);
        if(data_->factorInput().rows()==0||directBaseline||data_->factor.info()!=Eigen::Success)throw;
        data_->factor.forceDirect(trace,3);
      }
    }
}
}
