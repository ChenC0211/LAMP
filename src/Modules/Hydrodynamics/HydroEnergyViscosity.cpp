#include "Modules/Hydrodynamics/HydroEnergyViscosity.h"
#include "Modules/Hydrodynamics/HydroMomentumGeometry.h"

#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/NativeLeafOverlap.h"
#include <tuple>
#include <type_traits>
#include <Eigen/QR>
#include <Eigen/SVD>
#include <algorithm>
#include <atomic>
#include <array>
#include <cfenv>
#include <chrono>
#include <cstdint>
#include <system_error>
#include <cerrno>
#include <exception>
#include <memory>
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
#include <omp.h>
#endif
#include <cmath>
#include <cstring>
#include <limits>
#include <map>
#include <mutex>
#include <stdexcept>
#include <string>
#include <utility>
#include <iterator>
#if defined(_MSC_VER) && defined(_M_X64)
#include <xmmintrin.h>
#include <float.h>
#endif

namespace LAMP::Hydrodynamics
{
namespace
{
using Matrix=HydroStrainOperator::Matrix;
using Axis=HydroDofLayout::NormalAxis;
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
// Same qualified native/raw worker FP protocol as Maxwell. UCRT's
// FE_ALL_EXCEPT omits denormal; never synthesize status with feraiseexcept.
struct StrainFloatingState
{
    std::fenv_t native{};unsigned csr=0;bool valid=false;
    bool capture() noexcept
    {csr=_mm_getcsr();valid=std::fegetenv(&native)==0;return valid;}
    void restore() const noexcept
    {if(valid&&std::fesetenv(&native)!=0)std::terminate();_mm_setcsr(csr);}
    bool enterFresh(const StrainFloatingState& source) const noexcept
    {
        auto clean=source.native;clean._Fe_stat=0;
        if(std::fesetenv(&clean)!=0)return false;
        _mm_setcsr(source.csr&~0x3fu);return true;
    }
};
bool strainParallelQualified(StrainFloatingState& caller) noexcept
{
    if(omp_in_parallel()||omp_get_max_threads()<=1||!caller.capture())return false;
    unsigned int control=0;
    if(_controlfp_s(&control,0,0)!=0||(control&_MCW_EM)!=_MCW_EM
        ||(caller.csr&0x1f80u)!=0x1f80u)return false;
#if defined(_EM_AMBIGUOUS)
    if(control&_EM_AMBIGUOUS)return false;
#endif
    return true;
}
#endif
void require(bool ok,const std::string& reason)
{
    if (!ok) throw std::invalid_argument("[HydroEnergyViscosity] "+reason);
}
// Literal success checks do not allocate a temporary diagnostic string.
// Failure still executes the original string helper and exception exactly.
void require(bool ok,const char* reason)
{
    if(!ok)require(false,std::string(reason));
}
void validate(const QuadtreeMesh& mesh,const HydroNativeValues& values,const HydroPhysicalDivergence& metrics)
{
    values.validateFor(metrics.layout(),mesh);
    for (int u : metrics.axisVelocityDofs())
        require(values.normalVelocity()[u]==0,"nonzero radial axis velocity");
}
std::vector<double> multiply(const Matrix& matrix,const std::vector<double>& input)
{
    if(matrix.rows()==0)return {};
    const Eigen::Map<const Eigen::VectorXd> u(input.data(),static_cast<Eigen::Index>(input.size()));
    const Eigen::VectorXd result=matrix*u;
    require(result.allFinite(),"nonfinite operator result");
    return {result.data(),result.data()+result.size()};
}
struct LocalSample {int dof;double sign;Vector2D point;};
struct Fit {std::vector<LocalSample> samples;Eigen::MatrixXd coefficients;};
// Rows retain map key order and explicit zero entries. Only storage changes:
// every caller still performs its original long-double update in that order.
class StressGeometryRow
{
    using Entry=std::pair<int,long double>;
    std::vector<Entry> entries_;
public:
    long double& operator[](int key)
    {
        auto found=std::lower_bound(entries_.begin(),entries_.end(),key,
            [](const Entry& entry,int value) {return entry.first<value;});
        if (found==entries_.end() || found->first!=key)
            found=entries_.insert(found,Entry{key,0.L});
        return found->second;
    }
    auto begin() const noexcept {return entries_.begin();}
    auto end() const noexcept {return entries_.end();}
};
#if defined(_MSC_VER) && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST)
#pragma float_control(precise,on,push)
#pragma float_control(except,on)
#pragma fenv_access(on)
#endif
// Each clock access restores its CURRENT entry, including flags from any
// intervening physical assembly. A timer never restores its start environment.
struct ViscosityPreparationClock
{
    using Clock=std::chrono::steady_clock;
    Clock::time_point start{};bool available=false;
    ViscosityPreparationClock() noexcept
    {
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
        const int error=errno;StrainFloatingState current;
        if(!current.capture()){errno=error;return;}
        try{start=Clock::now();available=true;}catch(...){}
        current.restore();errno=error;
#endif
    }
    double seconds() noexcept
    {
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
        if(!available)return 0;
        const int error=errno;StrainFloatingState current;
        if(!current.capture()){available=false;errno=error;return 0;}
        double value=0;
        try{value=std::chrono::duration<double>(Clock::now()-start).count();
            if(!std::isfinite(value)||value<0){available=false;value=0;}}
        catch(...){available=false;value=0;}
        current.restore();errno=error;return value;
#else
        return 0;
#endif
    }
};
#if defined(_MSC_VER) && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST)
#if !defined(_M_FP_STRICT) || !_M_FP_STRICT
#pragma fenv_access(off)
#endif
#pragma float_control(pop)
#endif
// Eigen 3.3.9's ColMajor/RowMajor -> RowMajor conservative product
// converts the RHS to ColMajor, then visits each coefficient's k in ascending
// order. Only independent output coefficients are reordered by this plan.
// The symbolic indices are phase-free; an optional bounded numeric stock below
// retains exact current inputs/values/row flags, never a physical owner or history link.
bool viscosityProductPlanEnvironment() noexcept
{
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_FP_FAST) && !defined(__FAST_MATH__) && EIGEN_WORLD_VERSION==3 && ((EIGEN_MAJOR_VERSION==3 && EIGEN_MINOR_VERSION==9) || (EIGEN_MAJOR_VERSION==4 && EIGEN_MINOR_VERSION==1))
    unsigned control=0;
    // mask=0 is a read-only query. Do not modify the caller's FP environment.
    if(_controlfp_s(&control,0,0)!=0||(control&_MCW_EM)!=_MCW_EM)return false;
#if defined(_EM_AMBIGUOUS)
    if(control&_EM_AMBIGUOUS)return false;
#endif
    // SSE exception masks (invalid, denormal, zero-divide, overflow,
    // underflow, precision) must ALL be set before independent reordering.
    return (_mm_getcsr()&0x1f80u)==0x1f80u;
#else
    return false; // Unknown Eigen/order/trap contract keeps the original path.
#endif
}
struct ViscosityProductPlan
{
    using Matrix=HydroStrainOperator::Matrix;
    using Index=Matrix::StorageIndex;
    using Route=std::uint32_t;
    std::vector<Route> columnOffsets,strainRows,entryOffsets,outer,inner;
    std::size_t retainedBytes=0,preparationPeakBytes=0;
    bool timingAvailable=false;int symbolicWorkers=0;
    double columnRoutingSeconds=0,symbolicGraphSeconds=0;
    explicit ViscosityProductPlan(const Matrix& strain)
    {
        static_assert(Matrix::IsRowMajor,"viscosity product plan requires row-major raw strain");
        ViscosityPreparationClock columnClock;
        constexpr auto routeLimit=std::numeric_limits<Route>::max();
        if(!strain.isCompressed()||strain.cols()<0||strain.rows()<0||strain.nonZeros()<0
            ||static_cast<std::uint64_t>(strain.cols())>routeLimit
            ||static_cast<std::uint64_t>(strain.rows())>routeLimit
            ||static_cast<std::uint64_t>(strain.nonZeros())>routeLimit)
            throw std::length_error("optional viscosity symbolic dimensions");
        constexpr std::size_t limit=128u*1024u*1024u;
        const auto add=[&](std::size_t& bytes,std::size_t count,std::size_t width){
            if(bytes>limit||count>(limit-bytes)/width)
                throw std::length_error("optional viscosity symbolic memory scope");
            bytes+=count*width;
        };
        const auto count=static_cast<std::size_t>(strain.cols());
        const auto nonzeros=static_cast<std::size_t>(strain.nonZeros());
        std::size_t requestedBytes=sizeof(ViscosityProductPlan)+512u;
        add(requestedBytes,count+1,sizeof(Route));
        add(requestedBytes,nonzeros,2*sizeof(Route));
        add(requestedBytes,count,sizeof(Route)); // Temporary column cursors.
        columnOffsets.assign(count+1,0);strainRows.resize(nonzeros);entryOffsets.resize(nonzeros);
        retainedBytes=sizeof(ViscosityProductPlan)+512u;
        add(retainedBytes,columnOffsets.capacity(),sizeof(Route));
        add(retainedBytes,strainRows.capacity(),sizeof(Route));add(retainedBytes,entryOffsets.capacity(),sizeof(Route));
        const auto* sourceOuter=strain.outerIndexPtr();const auto* sourceInner=strain.innerIndexPtr();
        for(Eigen::Index k=0;k<strain.rows();++k)
            for(Index p=sourceOuter[k];p<sourceOuter[k+1];++p)
                ++columnOffsets[static_cast<std::size_t>(sourceInner[p])+1];
        for(std::size_t i=0;i<count;++i)columnOffsets[i+1]+=columnOffsets[i];
        {
            std::vector<Route> cursor(columnOffsets.begin(),columnOffsets.end()-1);
            std::size_t bytes=retainedBytes;add(bytes,cursor.capacity(),sizeof(Route));
            preparationPeakBytes=bytes;
            // Ascending original k gives each column the original lhs order.
            // Values remain in this immutable raw E; only row/entry IDs persist.
            for(Eigen::Index k=0;k<strain.rows();++k)
                for(Index p=sourceOuter[k];p<sourceOuter[k+1];++p){
                    const auto position=cursor[static_cast<std::size_t>(sourceInner[p])]++;
                    strainRows[position]=static_cast<Route>(k);entryOffsets[position]=static_cast<Route>(p);
                }
        } // Release cursors before allocating K routing and its markers.
        columnRoutingSeconds=columnClock.seconds();
        ViscosityPreparationClock symbolicClock;
        requestedBytes=retainedBytes;add(requestedBytes,count+1,sizeof(Route));
        outer.resize(count+1);add(retainedBytes,outer.capacity(),sizeof(Route));
        struct SymbolicWorker {
            std::vector<Route> marks;
            std::exception_ptr failure;
        };
        const auto runRows=[&](bool fill){
            if(count==0)return;
            std::size_t perWorker=sizeof(SymbolicWorker)+256u;
            add(perWorker,count,sizeof(Route));
            int requestedTeam=1;
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
            requestedTeam=std::max(1,omp_get_max_threads());
#endif
            const auto available=(limit-retainedBytes)/perWorker;
            const auto teamSize=std::min<std::size_t>(static_cast<std::size_t>(requestedTeam),std::min(count,available));
            if(teamSize==0||teamSize>static_cast<std::size_t>(std::numeric_limits<int>::max()))
                throw std::length_error("optional viscosity symbolic worker scope");
            std::vector<SymbolicWorker> workers(teamSize);
            std::size_t bytes=retainedBytes;add(bytes,workers.capacity(),sizeof(SymbolicWorker));
            for(auto& worker:workers){worker.marks.assign(count,routeLimit);
                add(bytes,worker.marks.capacity(),sizeof(Route));add(bytes,256u,1);}
            preparationPeakBytes=std::max(preparationPeakBytes,bytes);
            const auto row=[&](std::size_t i,SymbolicWorker& worker) noexcept {
                try{
                    Route position=fill?outer[i]:0;
                    for(Route q=columnOffsets[i];q<columnOffsets[i+1];++q)
                        for(Matrix::InnerIterator rhs(strain,strainRows[q]);rhs;++rhs){
                            auto& mark=worker.marks[static_cast<std::size_t>(rhs.col())];
                            if(mark==static_cast<Route>(i))continue;
                            mark=static_cast<Route>(i);
                            if(fill)inner[position++]=static_cast<Route>(rhs.col());else ++position;
                        }
                    if(fill)std::sort(inner.begin()+outer[i],inner.begin()+outer[i+1]);
                    else outer[i+1]=position;
                }catch(...){if(!worker.failure)worker.failure=std::current_exception();}
            };
            int actualWorkers=1,setupFailure=0;
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
            const int team=static_cast<int>(teamSize);
#pragma omp parallel num_threads(team)
            {
                auto& local=workers[static_cast<std::size_t>(omp_get_thread_num())];
                const int error=errno;StrainFloatingState current;
                const bool captured=current.capture();
                if(!captured){
#pragma omp atomic write
                    setupFailure=1;
                }
#pragma omp master
                {actualWorkers=omp_get_num_threads();}
#pragma omp for schedule(static) nowait
                for(Eigen::Index i=0;i<strain.cols();++i)if(captured)row(static_cast<std::size_t>(i),local);
                if(captured)current.restore();errno=error;
            }
#else
            for(std::size_t i=0;i<count;++i)row(i,workers[0]);
#endif
            if(setupFailure)throw std::runtime_error("optional viscosity symbolic worker environment");
            for(const auto& worker:workers)if(worker.failure)std::rethrow_exception(worker.failure);
            symbolicWorkers=symbolicWorkers==0?actualWorkers:std::min(symbolicWorkers,actualWorkers);
        };
        // Independent rows count all actual incidences, including explicit
        // zeros. Only this checked serial prefix combines row counts.
        runRows(false);outer[0]=0;
        std::uint64_t values=0;
        for(std::size_t i=0;i<count;++i){
            if(values>static_cast<std::uint64_t>(routeLimit)-outer[i+1])
                throw std::length_error("optional viscosity symbolic index scope");
            values+=outer[i+1];outer[i+1]=static_cast<Route>(values);
        }
        requestedBytes=retainedBytes;add(requestedBytes,static_cast<std::size_t>(values),sizeof(Route));
        inner.resize(static_cast<std::size_t>(values));add(retainedBytes,inner.capacity(),sizeof(Route));
        // Count markers are gone; fill is admitted against the final K routes.
        runRows(true);
        symbolicGraphSeconds=symbolicClock.seconds();
        timingAvailable=columnClock.available&&symbolicClock.available;
    }
    // One optional numeric slot belongs to this exact immutable E/index plan.
    // No Energy/State/phase/old Data pointer or lease escapes the exclusive lock.
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
    struct RowFlags {unsigned long native=0;unsigned csr=0;};
    struct NumericMemo {
        std::unique_ptr<double[]> weights;
        std::unique_ptr<RowFlags[]> flags;
        std::shared_ptr<const Matrix> matrix;
        RowFlags canonicalFlags;
        std::size_t weightCount=0,valueCount=0,rowCount=0,bytes=0;
        StrainFloatingState environment;unsigned control=0;
    };
    mutable std::mutex numericMutex;
    mutable std::unique_ptr<NumericMemo> numericMemo;
    static bool sameNumericControl(const StrainFloatingState& a,const StrainFloatingState& b) noexcept
    {
        auto x=a.native,y=b.native;x._Fe_stat=0;y._Fe_stat=0;
        return std::memcmp(&x,&y,sizeof x)==0&&(a.csr&~0x3fu)==(b.csr&~0x3fu);
    }
    static bool numericEnvironment(StrainFloatingState& state,unsigned& control) noexcept
    {
#if !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST) && !defined(__AVX__) && !defined(__AVX2__) && !defined(__AVX512F__) && !defined(__FMA__)
        if(omp_in_parallel()||!state.capture()||!viscosityProductPlanEnvironment()
            ||_controlfp_s(&control,0,0)!=0)return false;
        return sizeof(double)==8&&std::numeric_limits<double>::digits==53
            &&(control&_MCW_RC)==_RC_NEAR&&(state.csr&0x6000u)==0
            &&(state.csr&((1u<<15)|(1u<<6)))==0;
#else
        (void)state;(void)control;return false;
#endif
    }
    struct NumericPending {
        const ViscosityProductPlan& owner;
        std::unique_lock<std::mutex> lock;
        std::unique_ptr<NumericMemo> replacement;
        explicit NumericPending(const ViscosityProductPlan& plan)
            :owner(plan),lock(plan.numericMutex,std::try_to_lock){}
        // Publish only AFTER original prune, finite and source checks. This
        // optional stock update must not erase any flags raised by those gates.
        void publish(bool originalSourceQualified,HydroViscosityProductReport& report) noexcept
        {
            const int error=errno;StrainFloatingState current;unsigned control=0;
            const bool qualified=numericEnvironment(current,control);
            if(originalSourceQualified&&qualified&&replacement&&control==replacement->control
                &&sameNumericControl(current,replacement->environment)){
                owner.numericMemo=std::move(replacement);report.memoPublished=true;
                report.retainedBytes=owner.numericMemo->bytes;
            }
            if(current.valid)current.restore();errno=error;
        }
    };
    struct NumericPendingDeleter {
        void operator()(NumericPending* pending) const noexcept {
            const int error=errno;StrainFloatingState current;current.capture();
            delete pending;current.restore();errno=error;
        }
    };
    using NumericPendingHandle=std::unique_ptr<NumericPending,NumericPendingDeleter>;
#else
    struct NumericPending {void publish(bool,HydroViscosityProductReport&) noexcept {}};
    using NumericPendingHandle=std::unique_ptr<NumericPending>;
#endif
    NumericPendingHandle assembleRemembered(const Matrix& strain,
        const std::vector<double>& weights,std::shared_ptr<const Matrix>& result,HydroViscosityProductReport& report) const;
};
#if defined(_MSC_VER) && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST)
    // The surrounding legacy TU is /fp:precise. Only this new row-receipt
    // function observes/clears exception state between rows; make that local
    // ordering part of its compiler contract, then restore the original mode.
#pragma float_control(precise,on,push)
#pragma float_control(except,on)
#pragma fenv_access(on)
#endif
ViscosityProductPlan::NumericPendingHandle ViscosityProductPlan::assembleRemembered(const Matrix& strain,
    const std::vector<double>& weights,std::shared_ptr<const Matrix>& result,HydroViscosityProductReport& report) const
    {
        report.requested=true;report.fullAssembly=true;report.weightRows=weights.size();
        report.symbolicStoredValues=inner.size();
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
        const int callerErrno=errno;StrainFloatingState caller;unsigned control=0;
        const bool qualified=numericEnvironment(caller,control);errno=callerErrno;
        report.currentFpQualified=qualified;
        if(!qualified){if(caller.valid)caller.restore();errno=callerErrno;report.fallbackReason=2;return {};}
        try {
            NumericPendingHandle pending(new NumericPending(*this));
            if(!pending->lock.owns_lock()){report.fallbackReason=4;throw std::runtime_error("optional numeric slot busy");}
            const auto n=strain.cols();
            if(n<=0||strain.rows()<0||weights.size()!=static_cast<std::size_t>(strain.rows()))
                throw std::length_error("optional numeric product dimensions");
            const auto count=static_cast<std::size_t>(n),valueCount=inner.size();
            const auto* source=numericMemo.get();report.sourceAvailable=source!=nullptr;
            const bool matching=source&&source->matrix&&source->weightCount==weights.size()&&source->valueCount==valueCount
                &&source->rowCount==count&&source->matrix->rows()==n&&source->matrix->cols()==n&&source->matrix->isCompressed()
                &&source->control==control&&sameNumericControl(caller,source->environment);
            report.sourceFpReceipt=matching;
            constexpr std::size_t limit=128u*1024u*1024u;
            const auto addBytes=[&](std::size_t& total,std::size_t items,std::size_t width){
                if(total>limit||items>(limit-total)/width){report.fallbackReason=3;throw std::length_error("optional numeric product memory scope");}
                total+=items*width;
            };
            std::size_t memoBytes=sizeof(NumericMemo)+512u;
            addBytes(memoBytes,weights.size(),sizeof(double));addBytes(memoBytes,count,sizeof(RowFlags));
            const auto matrixBytes=[&](const Matrix& matrix){
                if(!matrix.isCompressed()||matrix.rows()!=n||matrix.cols()!=n||matrix.data().allocatedSize()<0)
                    throw std::length_error("optional numeric actual matrix stock");
                std::size_t bytes=sizeof(Matrix)+512u;
                addBytes(bytes,static_cast<std::size_t>(matrix.data().allocatedSize()),sizeof(double)+sizeof(Index));
                addBytes(bytes,static_cast<std::size_t>(matrix.outerSize())+1,sizeof(Index));return bytes;
            };
            const auto oldMatrixBytes=source&&source->matrix?matrixBytes(*source->matrix):0;
            std::size_t newMatrixBytes=sizeof(Matrix)+512u;
            addBytes(newMatrixBytes,valueCount,sizeof(double)+sizeof(Index));addBytes(newMatrixBytes,count+1,sizeof(Index));
            struct Scratch {
                std::unique_ptr<Route[]> slot;std::unique_ptr<unsigned char[]> touched;
                std::exception_ptr failure;
            };
            // Admit both metadata stocks before allocating optional buffers.
            std::size_t bytes=sizeof(NumericPending)+512u;addBytes(bytes,retainedBytes,1);addBytes(bytes,2,memoBytes);
            addBytes(bytes,oldMatrixBytes,1);addBytes(bytes,count,sizeof(unsigned char));
            auto next=std::make_unique<NumericMemo>();next->weightCount=weights.size();next->valueCount=valueCount;
            next->rowCount=count;next->environment=caller;next->control=control;
            next->weights=std::make_unique<double[]>(weights.size());next->flags=std::make_unique<RowFlags[]>(count);
            std::copy(weights.begin(),weights.end(),next->weights.get());
            auto dirty=std::make_unique<unsigned char[]>(count);
            if(!matching)std::fill_n(dirty.get(),count,1);
            std::size_t changed=0;
            for(std::size_t k=0;k<weights.size();++k){
                std::uint64_t oldBits=0,newBits=0;std::memcpy(&newBits,&weights[k],sizeof newBits);
                if(matching)std::memcpy(&oldBits,source->weights.get()+k,sizeof oldBits);
                if(!matching||oldBits!=newBits){++changed;
                    for(Matrix::InnerIterator e(strain,static_cast<Eigen::Index>(k));e;++e)
                        dirty[static_cast<std::size_t>(e.col())]=1;
                }
            }
            std::size_t recomputed=0;for(std::size_t i=0;i<count;++i)recomputed+=dirty[i]!=0;
            // The fresh K is mandatory physical output. Only one distinct old
            // matrix is optional stock; at publication that role switches to K.
            // Retain the larger actual footprint, plus both simultaneous small
            // metadata stocks. Shared indices/values/outer are never charged zero.
            if(recomputed){bytes-=oldMatrixBytes;addBytes(bytes,std::max(oldMatrixBytes,newMatrixBytes),1);}
            report.peakOptionalBytes=bytes;
            if(!recomputed){
                next->matrix=source->matrix;next->canonicalFlags=source->canonicalFlags;
                std::copy_n(source->flags.get(),count,next->flags.get());
                unsigned long raisedNative=next->canonicalFlags.native;unsigned raisedCsr=next->canonicalFlags.csr;
                for(std::size_t i=0;i<count;++i){raisedNative|=next->flags[i].native;raisedCsr|=next->flags[i].csr;}
                next->bytes=memoBytes;addBytes(next->bytes,oldMatrixBytes,1);
                dirty.reset();pending->replacement=std::move(next);result=pending->replacement->matrix;
                auto merged=caller.native;merged._Fe_stat|=raisedNative;
                if(std::fesetenv(&merged)!=0)std::terminate();_mm_setcsr(caller.csr|raisedCsr);errno=callerErrno;
                report.eligible=true;report.localUsed=true;report.fullAssembly=false;
                report.changedWeightRows=changed;report.outputRows=count;report.recomputedRows=0;report.reusedRows=count;
                report.retainedBytes=source->bytes;report.fallbackReason=0;return pending;
            }
            std::size_t perWorker=sizeof(Scratch)+256u;addBytes(perWorker,count,sizeof(Route)+sizeof(unsigned char));
            const auto available=(limit-bytes)/perWorker;
            const auto teamSize=std::min<std::size_t>(static_cast<std::size_t>(std::max(1,omp_get_max_threads())),std::min(count,available));
            if(teamSize==0||teamSize>static_cast<std::size_t>(std::numeric_limits<int>::max())){
                report.fallbackReason=3;throw std::length_error("optional numeric product scratch scope");}
            addBytes(bytes,teamSize,perWorker);report.peakOptionalBytes=bytes;
            auto scratch=std::make_unique<Scratch[]>(teamSize);
            for(std::size_t t=0;t<teamSize;++t){scratch[t].slot=std::make_unique<Route[]>(count);scratch[t].touched=std::make_unique<unsigned char[]>(count);}
            auto actual=std::make_shared<Matrix>();actual->resize(n,n);actual->makeCompressed();
            actual->resizeNonZeros(static_cast<Eigen::Index>(valueCount));
            std::copy(outer.begin(),outer.end(),actual->outerIndexPtr());std::copy(inner.begin(),inner.end(),actual->innerIndexPtr());
            const auto actualMatrixBytes=matrixBytes(*actual);
            bytes-=std::max(oldMatrixBytes,newMatrixBytes);
            addBytes(bytes,std::max(oldMatrixBytes,actualMatrixBytes),1);
            report.peakOptionalBytes=bytes;
            const int team=static_cast<int>(teamSize);int setupFailure=0;
#pragma omp parallel num_threads(team)
            {
                auto& local=scratch[static_cast<std::size_t>(omp_get_thread_num())];
                const int workerErrno=errno;StrainFloatingState worker;const bool captured=worker.capture();
                const bool ready=captured&&worker.enterFresh(caller);
                if(!ready){
#pragma omp atomic write
                    setupFailure=1;
                }
                auto& slot=local.slot;auto& touched=local.touched;
#pragma omp for schedule(static) nowait
                for(Eigen::Index i=0;i<n;++i){
                    if(!ready)continue;
                    try {
                        const auto row=static_cast<std::size_t>(i);
                        if(!dirty[row]){
                            const auto& old=*source->matrix;Index q=old.outerIndexPtr()[row],end=old.outerIndexPtr()[row+1];
                            for(Index p=outer[row];p<outer[row+1];++p){
                                const auto column=inner[static_cast<std::size_t>(p)];
                                if(q<end&&old.innerIndexPtr()[q]==column)actual->valuePtr()[p]=old.valuePtr()[q++];
                                else actual->valuePtr()[p]=0.;
                            }
                            if(q!=end)throw std::runtime_error("optional numeric pruned source route");
                            next->flags[row]=source->flags[row];continue;
                        }
                        if(!worker.enterFresh(caller))throw std::runtime_error("optional numeric row fresh environment");
                        errno=callerErrno;
            for(Index p=outer[static_cast<std::size_t>(i)];p<outer[static_cast<std::size_t>(i+1)];++p)
                slot[static_cast<std::size_t>(inner[static_cast<std::size_t>(p)])]=static_cast<Route>(p);
            // Integer column routes visit k in the same ascending order used
            // by Eigen. Keep the original separate binary64 row scaling before
            // x*y. Its flags join this row receipt too.
            for(Route q=columnOffsets[static_cast<std::size_t>(i)];q<columnOffsets[static_cast<std::size_t>(i)+1];++q){
                const Eigen::Index k=strainRows[q];const double x=strain.valuePtr()[entryOffsets[q]];
                for(Matrix::InnerIterator rhs(strain,k);rhs;++rhs){
                    const Eigen::Index j=rhs.col();const double y=rhs.value()*weights[static_cast<std::size_t>(k)];
                    auto& value=actual->valuePtr()[slot[static_cast<std::size_t>(j)]];
                    if(!touched[static_cast<std::size_t>(j)]){
                        touched[static_cast<std::size_t>(j)]=1;value=x*y;
                    }else value+=x*y;
                }
            }
            for(Index p=outer[static_cast<std::size_t>(i)];p<outer[static_cast<std::size_t>(i+1)];++p){
                touched[static_cast<std::size_t>(inner[static_cast<std::size_t>(p)])]=0;
                std::uint64_t bits=0;std::memcpy(&bits,actual->valuePtr()+p,sizeof bits);
                if((bits&UINT64_C(0x7ff0000000000000))==UINT64_C(0x7ff0000000000000))
                    throw std::runtime_error("optional numeric unpruned nonfinite");
            }
                        const int arithmeticErrno=errno;StrainFloatingState after;
                        if(!after.capture()||arithmeticErrno!=callerErrno||!sameNumericControl(after,caller))
                            throw std::runtime_error("optional numeric row status receipt");
                        next->flags[row]={after.native._Fe_stat,after.csr&0x3fu};
                    }catch(...){if(!local.failure)local.failure=std::current_exception();}
                }
                if(captured)worker.restore();errno=workerErrno;
            }
            if(setupFailure)throw std::runtime_error("optional numeric worker capture");
            for(std::size_t t=0;t<teamSize;++t)if(scratch[t].failure)std::rethrow_exception(scratch[t].failure);
            unsigned long raisedNative=0;unsigned raisedCsr=0;
            for(std::size_t i=0;i<count;++i){raisedNative|=next->flags[i].native;raisedCsr|=next->flags[i].csr;}
            scratch.reset();dirty.reset();
            auto merged=caller.native;merged._Fe_stat|=raisedNative;
            if(std::fesetenv(&merged)!=0)std::terminate();_mm_setcsr(caller.csr|raisedCsr);errno=callerErrno;
            // One actual canonicalization, before any const stock escapes.
            // Capture its genuine raises once so an all-clean call can share
            // the immutable K without a dummy prune/finite scan.
            StrainFloatingState beforeCanonical;
            if(!beforeCanonical.capture()||!beforeCanonical.enterFresh(caller))
                throw std::runtime_error("optional numeric canonical environment");
            actual->prune(0.0);actual->makeCompressed();
            StrainFloatingState afterCanonical;
            if(!afterCanonical.capture()||!sameNumericControl(afterCanonical,caller))
                throw std::runtime_error("optional numeric canonical receipt");
            const int canonicalErrno=errno;
            next->canonicalFlags={afterCanonical.native._Fe_stat,afterCanonical.csr&0x3fu};
            auto completed=beforeCanonical.native;completed._Fe_stat|=next->canonicalFlags.native;
            if(std::fesetenv(&completed)!=0)std::terminate();_mm_setcsr(beforeCanonical.csr|next->canonicalFlags.csr);errno=canonicalErrno;
            next->bytes=memoBytes;addBytes(next->bytes,matrixBytes(*actual),1);
            next->matrix=std::move(actual);result=next->matrix;pending->replacement=std::move(next);
            report.eligible=true;report.localUsed=recomputed<count;report.fullAssembly=recomputed==count;
            report.changedWeightRows=changed;report.outputRows=count;report.recomputedRows=recomputed;report.reusedRows=count-recomputed;
            report.retainedBytes=source?source->bytes:0;report.fallbackReason=0;
            return pending;
        }catch(...){
            // Optional state has fully unwound before the caller executes the
            // entire original Task32b path. This also restores original error order.
            caller.restore();errno=callerErrno;if(!report.fallbackReason)report.fallbackReason=5;
            report.eligible=false;report.localUsed=false;report.fullAssembly=true;return {};
        }
#else
        (void)strain;(void)weights;(void)result;report.fallbackReason=2;return {};
#endif
    }
#if defined(_MSC_VER) && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST)
#if !defined(_M_FP_STRICT)
#pragma fenv_access(off)
#endif
#pragma float_control(pop)
#endif
// Complete producer numerical-environment receipt for local AMR block reuse.
// A source computed under UPWARD/DAZ/FTZ/traps cannot be relabelled as RNE.
struct StressReuseEnvironment
{
    std::fenv_t native{};unsigned csr=0,control=0;
    bool read() noexcept
    {
#if defined(_MSC_VER) && defined(_M_X64) && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST)
        const int saved=errno;std::memset(&native,0,sizeof(native));csr=_mm_getcsr();
        const bool captured=std::fegetenv(&native)==0;
        const bool queried=_controlfp_s(&control,0,0)==0;errno=saved;
        if(!captured||!queried||sizeof(double)!=8||std::numeric_limits<double>::digits!=53)return false;
#if defined(_EM_AMBIGUOUS)
        if(control&_EM_AMBIGUOUS)return false;
#endif
        return std::fegetround()==FE_TONEAREST&&(control&_MCW_RC)==_RC_NEAR
            &&(control&_MCW_EM)==_MCW_EM&&(csr&0x1f80u)==0x1f80u
            &&(csr&0x6000u)==0&&(csr&((1u<<15)|(1u<<6)))==0;
#else
        return false;
#endif
    }
    void restore() const noexcept
    {
        if(std::fesetenv(&native)!=0)std::terminate();
#if defined(_MSC_VER) && defined(_M_X64)
        _mm_setcsr(csr);
#endif
    }
    bool same(const StressReuseEnvironment& other) const noexcept
    {return csr==other.csr&&control==other.control&&std::memcmp(&native,&other.native,sizeof(native))==0;}
    bool sameControl(const StressReuseEnvironment& other) const noexcept
    {
#if defined(_MSC_VER) && defined(_M_X64)
        return control==other.control&&(csr&~0x3fu)==(other.csr&~0x3fu)&&native._Fe_ctl==other.native._Fe_ctl;
#else
        return false;
#endif
    }
};
// This cache belongs to one buildImpl only. It never retains actual sample
// owners: an identical normalized system may reuse only its successful fit.
struct FitTemplates
{
#if defined(_MSC_VER) && defined(_M_X64)
    static_assert(sizeof(double)==8 && std::numeric_limits<double>::is_iec559,
        "P2 templates require binary64 input bits");
    using Key=std::vector<unsigned char>;
    struct Environment {
        std::array<unsigned char,sizeof(std::fenv_t)> bytes{};
        unsigned csr=0;
    };
    std::map<Key,Eigen::MatrixXd> coefficients;
    static bool environment(Environment& result)
    {
        std::fenv_t state;
        std::memset(&state,0,sizeof(state));
        if (std::fegetenv(&state)!=0) return false;
        std::memcpy(result.bytes.data(),&state,sizeof(state));
        result.csr=_mm_getcsr();
        return true;
    }
    static bool same(const Environment& a,const Environment& b)
    {
        return a.bytes==b.bytes && a.csr==b.csr;
    }
    static void append(Key& key,const void* data,std::size_t size)
    {
        if (size==0) return;
        const auto* first=static_cast<const unsigned char*>(data);
        key.insert(key.end(),first,first+size);
    }
    template<class T> static void append(Key& key,const T& value)
    {
        append(key,&value,sizeof(value));
    }
#endif
};
Eigen::Matrix<double,1,6> basis(double x,double y)
{
    Eigen::Matrix<double,1,6> result;result<<1,x,y,x*x,x*y,y*y;return result;
}
Fit fit(const HydroPhysicalDivergence& metrics,const std::vector<NativeLeafFace>& faces,
    const std::vector<int>& support,const std::vector<int>& nativeFaces,Axis axis,
    Vector2D center,Vector2D width,bool reflected,int cell,FitTemplates& templates)
{
    Fit result;
    for (int u : support) if (metrics.layout().normalAxis(u)==axis) {
        result.samples.push_back({u,1,faces[u].center});
        if (reflected && faces[u].center.x>0)
            result.samples.push_back({u,metrics.reflectAxisNormalVelocity(axis,1),{-faces[u].center.x,faces[u].center.y}});
    }
    Eigen::MatrixXd a(result.samples.size(),6);
    for (Eigen::Index i=0;i<a.rows();++i) {
        const auto p=result.samples[static_cast<std::size_t>(i)].point;
        a.row(i)=basis((p.x-center.x)/width.x,(p.y-center.y)/width.y);
    }
#if defined(_MSC_VER) && defined(_M_X64)
    // Preserve a's original storage order and the native constraint order.
    // The original miss path below computes its own constraints unchanged.
    std::vector<Eigen::Index> templateConstraints;
    if (!reflected)
        for (int u : nativeFaces) if (metrics.layout().normalAxis(u)==axis)
            for (Eigen::Index j=0;j<a.rows();++j)
                if (result.samples[static_cast<std::size_t>(j)].dof==u) {templateConstraints.push_back(j);break;}
    FitTemplates::Key key;
    const auto rows=a.rows(),columns=a.cols();
    const auto count=templateConstraints.size();
    const int axisCode=static_cast<int>(axis);
    const unsigned char reflection=reflected ? 1 : 0;
    key.reserve(sizeof(rows)+sizeof(columns)+sizeof(count)+sizeof(axisCode)+sizeof(reflection)
        +count*sizeof(Eigen::Index)+static_cast<std::size_t>(a.size())*sizeof(double)
        +sizeof(std::fenv_t)+sizeof(unsigned));
    FitTemplates::append(key,rows);FitTemplates::append(key,columns);
    FitTemplates::append(key,count);FitTemplates::append(key,axisCode);FitTemplates::append(key,reflection);
    FitTemplates::append(key,templateConstraints.data(),count*sizeof(Eigen::Index));
    FitTemplates::append(key,a.data(),static_cast<std::size_t>(a.size())*sizeof(double));
    FitTemplates::Environment before;
    const bool cacheable=FitTemplates::environment(before);
    const Eigen::MatrixXd* cached=nullptr;
    if (cacheable) {
        FitTemplates::append(key,before.bytes.data(),before.bytes.size());
        FitTemplates::append(key,before.csr);
        const auto found=templates.coefficients.find(key);
        if (found!=templates.coefficients.end()) cached=&found->second;
    }
    if (cached) result.coefficients=*cached;
    else
#else
    (void)templates;
#endif
    {
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(a);
    qr.setThreshold(1e-12);
    require(qr.rank()==6,"rank-deficient local P2 support at cell "+std::to_string(cell));
    const Eigen::MatrixXd identity=Eigen::MatrixXd::Identity(a.rows(),a.rows());
    if (reflected) {
        // Enforce regular parity algebraically, not by cancellation of large
        // 1/r coefficients. Odd radial P2: r,r*z. Even axial P2: 1,z,r*r,z*z.
        const std::vector<int> columns=axis==Axis::First ? std::vector<int>{1,4} : std::vector<int>{0,2,3,5};
        Eigen::MatrixXd regular(a.rows(),columns.size());
        for (Eigen::Index j=0;j<regular.cols();++j) regular.col(j)=a.col(columns[static_cast<std::size_t>(j)]);
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> parity(regular);parity.setThreshold(1e-12);
        require(parity.rank()==regular.cols(),"rank-deficient axis-regular strain support");
        const Eigen::MatrixXd coefficients=parity.solve(identity);
        result.coefficients=Eigen::MatrixXd::Zero(6,a.rows());
        for (Eigen::Index j=0;j<regular.cols();++j) result.coefficients.row(columns[static_cast<std::size_t>(j)])=coefficients.row(j);
    }
    else {
        // Interpolate this cell's actual normal faces. On a regular cell this
        // preserves the central normal difference, rather than filtering it.
        std::vector<Eigen::Index> constraints;
        for (int u : nativeFaces) if (metrics.layout().normalAxis(u)==axis)
            for (Eigen::Index j=0;j<a.rows();++j)
                if (result.samples[static_cast<std::size_t>(j)].dof==u) {constraints.push_back(j);break;}
        Eigen::MatrixXd c(constraints.size(),6),selector=Eigen::MatrixXd::Zero(constraints.size(),a.rows());
        for (Eigen::Index i=0;i<c.rows();++i) {c.row(i)=a.row(constraints[static_cast<std::size_t>(i)]);selector(i,constraints[static_cast<std::size_t>(i)])=1;}
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(c,Eigen::ComputeThinU|Eigen::ComputeFullV);
        svd.setThreshold(1e-12);
        require(svd.rank()==c.rows(),"dependent native face constraints at cell "+std::to_string(cell));
        const Eigen::MatrixXd particular=svd.solve(selector);
        const Eigen::MatrixXd nullspace=svd.matrixV().rightCols(6-c.rows());
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> reduced(a*nullspace);reduced.setThreshold(1e-12);
        require(reduced.rank()==nullspace.cols(),"rank-deficient constrained strain support");
        result.coefficients=particular+nullspace*reduced.solve(identity-a*particular);
    }
    }
    require(result.coefficients.allFinite(),"nonfinite local reconstruction");
#if defined(_MSC_VER) && defined(_M_X64)
    // Do not erase or restore caller flags. A miss that raises any flag is not
    // reusable; only an already stable complete FP state can certify a hit.
    FitTemplates::Environment after;
    if (!cached && cacheable && FitTemplates::environment(after) && FitTemplates::same(before,after))
        templates.coefficients.emplace(std::move(key),result.coefficients);
#endif
    return result;
}
}
namespace
{
// Integer incidence/width and actual stored-coefficient magnitude belong only
// to one immutable strain CSR. No W, mu, K, phase, clock or previous owner.
// Called only after the original strict source/current Gram FP qualification;
// the enclosing query restores every caller FP/status/errno bit on all paths.
struct GramGeometryStatistics
{
    static constexpr std::size_t exactLimit=4503599627370494ULL;
    bool valid=false;
    std::size_t maximumIncidence=0,maximumWidth=0,total=0;
    double maximumEntry=0;
    static GramGeometryStatistics prepare(const Matrix& e)
    {
        GramGeometryStatistics result;
        if(!e.isCompressed()||e.rows()<0||e.cols()<0)return result;
        std::vector<std::size_t> incidence(static_cast<std::size_t>(e.cols()),0);
        for(Eigen::Index row=0;row<e.outerSize();++row){
            std::size_t width=0;
            for(Matrix::InnerIterator entry(e,row);entry;++entry){
                if(entry.col()<0||static_cast<std::size_t>(entry.col())>=incidence.size()
                    ||!std::isfinite(entry.value())||result.total>=exactLimit)return result;
                auto& count=incidence[static_cast<std::size_t>(entry.col())];
                if(count>=exactLimit)return result;
                ++count;++width;++result.total;result.maximumIncidence=std::max(result.maximumIncidence,count);
                result.maximumEntry=std::max(result.maximumEntry,std::abs(entry.value()));
            }
            result.maximumWidth=std::max(result.maximumWidth,width);
        }
        result.valid=true;return result;
    }
};
}
namespace
{
// Immutable actual geometry only; no Data/State/mesh/phase/clock/history link.
struct StressCanonicalPreparation
{
    static constexpr std::size_t byteLimit=64u*1024u*1024u;
    std::vector<Vector2D> centers,widths;
    std::vector<NativeLeafFace> faces;
    std::vector<std::vector<int>> cellFaces,neighbors,patches,supports;
    std::vector<std::vector<std::size_t>> pointRows;
    std::vector<unsigned char> complete;
    std::size_t pointCount=0,retainedBytes=0;
    bool matches(int cells,int velocities,std::size_t points) const noexcept
    {
        if(cells<0||velocities<0)return false;const auto n=static_cast<std::size_t>(cells);
        return centers.size()==n&&widths.size()==n&&cellFaces.size()==n&&neighbors.size()==n
            &&pointRows.size()==n&&patches.size()==n&&supports.size()==n&&complete.size()==n
            &&faces.size()==static_cast<std::size_t>(velocities)&&pointCount==points;
    }
};
struct StressPreparationScope
{
    int savedErrno=errno;std::fenv_t native{};bool captured=false;
#if defined(_MSC_VER) && defined(_M_X64)
    unsigned csr=_mm_getcsr();
#endif
    StressPreparationScope() noexcept {captured=std::fegetenv(&native)==0;errno=savedErrno;}
    ~StressPreparationScope() noexcept {
        if(captured&&std::fesetenv(&native)!=0)std::terminate();
#if defined(_MSC_VER) && defined(_M_X64)
        _mm_setcsr(csr);
#endif
        errno=savedErrno;
    }
};
// Each cell owns distinct vector/byte storage. Only budget/failure are shared.
struct StressPreparationRecorder
{
    using Preparation=StressCanonicalPreparation;
    std::vector<std::vector<int>> patches,supports;
    std::vector<unsigned char> complete;
    std::atomic<std::size_t> bytes{0};
    std::atomic<bool> disabled{false};
    std::size_t pointCount=0;
    ~StressPreparationRecorder() noexcept {
        StressPreparationScope scope;
        std::vector<std::vector<int>>().swap(patches);
        std::vector<std::vector<int>>().swap(supports);
        std::vector<unsigned char>().swap(complete);
    }
    static bool account(std::size_t& value,std::size_t count,std::size_t width) noexcept
    {
        if(value>Preparation::byteLimit||count>(Preparation::byteLimit-value)/width)return false;
        value+=count*width;return true;
    }
    template<class T>static bool accountVector(std::size_t& value,const std::vector<T>& v) noexcept
    {return account(value,v.capacity(),sizeof(T));}
    template<class T>static bool accountNested(std::size_t& value,const std::vector<std::vector<T>>& v) noexcept
    {
        if(!accountVector(value,v))return false;
        for(const auto& row:v)if(!accountVector(value,row))return false;
        return true;
    }
    bool reserve(std::size_t extra) noexcept
    {
        auto value=bytes.load(std::memory_order_relaxed);
        do{if(value>Preparation::byteLimit||extra>Preparation::byteLimit-value)return false;}
        while(!bytes.compare_exchange_weak(value,value+extra,std::memory_order_relaxed));
        return true;
    }
    static std::unique_ptr<StressPreparationRecorder> make(
        const std::vector<Vector2D>& centers,const std::vector<Vector2D>& widths,
        const std::vector<NativeLeafFace>& faces,const std::vector<std::vector<int>>& cellFaces,
        const std::vector<std::vector<int>>& neighbors,const std::vector<std::vector<std::size_t>>& pointRows,
        std::size_t points) noexcept
    {
        StressPreparationScope scope;if(!scope.captured)return {};
        try{
            std::size_t value=sizeof(Preparation)+256u;
            if(!accountVector(value,centers)||!accountVector(value,widths)||!accountVector(value,faces)
                ||!accountNested(value,cellFaces)||!accountNested(value,neighbors)||!accountNested(value,pointRows)
                ||centers.size()>(Preparation::byteLimit-value)/(2*sizeof(std::vector<int>)+sizeof(unsigned char)))return {};
            auto result=std::make_unique<StressPreparationRecorder>();
            result->patches.resize(centers.size());result->supports.resize(centers.size());
            result->complete.assign(centers.size(),0);result->pointCount=points;
            if(!accountVector(value,result->patches)||!accountVector(value,result->supports)||!accountVector(value,result->complete))return {};
            result->bytes.store(value,std::memory_order_relaxed);return result;
        }catch(const std::bad_alloc&){return {};}catch(const std::length_error&){return {};}
    }
    void save(int cell,const std::vector<int>& patch,const std::vector<int>& support) noexcept
    {
        if(disabled.load(std::memory_order_relaxed)||cell<0||static_cast<std::size_t>(cell)>=complete.size()||complete[cell])return;
        StressPreparationScope scope;if(!scope.captured){disabled.store(true,std::memory_order_relaxed);return;}
        try{
            const auto maximum=Preparation::byteLimit/sizeof(int);
            if(patch.size()>maximum||support.size()>maximum-patch.size()
                ||!reserve((patch.size()+support.size())*sizeof(int))){disabled.store(true,std::memory_order_relaxed);return;}
            std::vector<int> p(patch),s(support);
            const auto extra=(p.capacity()-patch.size())+(s.capacity()-support.size());
            if(extra>maximum||!reserve(extra*sizeof(int))){disabled.store(true,std::memory_order_relaxed);return;}
            patches[cell]=std::move(p);supports[cell]=std::move(s);complete[cell]=1;
        }catch(const std::bad_alloc&){disabled.store(true,std::memory_order_relaxed);}
        catch(const std::length_error&){disabled.store(true,std::memory_order_relaxed);}
    }
    std::shared_ptr<const Preparation> finish(
        std::vector<Vector2D>& centers,std::vector<Vector2D>& widths,std::vector<NativeLeafFace>& faces,
        std::vector<std::vector<int>>& cellFaces,std::vector<std::vector<int>>& neighbors,
        std::vector<std::vector<std::size_t>>& pointRows) noexcept
    {
        if(disabled.load(std::memory_order_relaxed))return {};
        try{
            // Allocate before moving; resource failure leaves business arrays
            // unchanged. std::allocator vector moves below are nonthrowing.
            auto result=std::make_shared<Preparation>();
            result->centers=std::move(centers);result->widths=std::move(widths);result->faces=std::move(faces);
            result->cellFaces=std::move(cellFaces);result->neighbors=std::move(neighbors);result->pointRows=std::move(pointRows);
            result->patches=std::move(patches);result->supports=std::move(supports);result->complete=std::move(complete);
            result->pointCount=pointCount;result->retainedBytes=bytes.load(std::memory_order_relaxed);return result;
        }catch(const std::bad_alloc&){return {};}catch(const std::length_error&){return {};}
    }
};
}
struct HydroStrainOperator::Data
{
    explicit Data(HydroPhysicalDivergence input) : metrics(std::move(input)) {}
    HydroPhysicalDivergence metrics;
    std::vector<HydroStrainPoint> points;
    std::vector<int> stationaryNoSlipWalls;
    Matrix matrix;
    // This actual raw geometry owns its symbolic plan and optional single numeric stock.
    // No EnergyViscosity/phase/previous-K link can form a history chain.
    mutable std::mutex viscosityProductMutex;
    mutable std::shared_ptr<const ViscosityProductPlan> viscosityProductPlan;
    mutable std::optional<GramGeometryStatistics> gramGeometryStatistics;
    int components=3;
    bool stressCompatible=false,rt0GreenOnly=false,localReuseEligible=false;
    StressReuseEnvironment localReuseBuildEntry,localReuseEnvironment;
    HydroStressGeometryReuseReport reuse;
    std::vector<HydroStressFacePoint> facePoints;
    Matrix faceVelocity,volumeVelocity,componentIntegrals;
    HydroStressCompatibilityReport compatibility;
    // Immutable phase-free raw CSR: no Fit, stress, K, clock or back-reference.
    std::shared_ptr<const Data> canonicalRaw;
    std::shared_ptr<const StressCanonicalPreparation> preparation;
};
namespace Detail
{
// One optional phase-free stock belongs to ONE actual Momentum Data. It never
// owns a State, material/clock, EnergyViscosity, or a Momentum back-reference.
struct HydroStationaryNoSlipStrainCache
{
    static constexpr std::size_t byteLimit=512u*1024u*1024u;
    static constexpr std::size_t numericReservation=128u*1024u*1024u;
    std::mutex mutex;
    std::optional<HydroStrainOperator> raw;
    StressReuseEnvironment environment;
    int error=0;
    std::size_t retainedBytes=0;
};
}
#if defined(_MSC_VER) && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST)
#pragma float_control(precise,on,push)
#pragma float_control(except,on)
#pragma fenv_access(on)
#endif
std::shared_ptr<Detail::HydroStationaryNoSlipStrainCache>
HydroMomentumGeometry::makeStationaryNoSlipStrainCache() noexcept
{
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
    const int error=errno;StrainFloatingState caller;
    if(!caller.capture()){errno=error;return {};}
    std::shared_ptr<Detail::HydroStationaryNoSlipStrainCache> result;
    try{result=std::make_shared<Detail::HydroStationaryNoSlipStrainCache>();}catch(...){}
    caller.restore();errno=error;return result;
#else
    return {};
#endif
}
#if defined(_MSC_VER) && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST)
#if !defined(_M_FP_STRICT) || !_M_FP_STRICT
#pragma fenv_access(off)
#endif
#pragma float_control(pop)
#endif

namespace
{
// Transient source views and exact native rebindings. This plan is destroyed
// before the fresh target operator returns; no old Data/State chain is stored.
struct StressLocalReuseMismatch {}; // Optional certificate mismatch, never a physical guard.
struct StressLocalReuse
{
    const Matrix *raw=nullptr,*stress=nullptr,*volume=nullptr,*faceVelocity=nullptr,*integrals=nullptr;
    const std::vector<HydroStressFacePoint>* sourceFacePoints=nullptr;
    std::vector<int> oldCell,oldFace,toTarget;
    std::vector<Eigen::Index> oldPoint;
    std::vector<std::vector<std::size_t>> oldFacePoints;
    std::vector<unsigned char> reuseRaw,reuseStress,reuseFace,needFit;
    bool enabled=false;
    bool columns(const Matrix& m,Eigen::Index row) const
    {
        for(Matrix::InnerIterator e(m,row);e;++e)
            if(e.col()<0||static_cast<std::size_t>(e.col())>=toTarget.size()||toTarget[e.col()]<0||!std::isfinite(e.value()))return false;
        return true;
    }
    template<class Emit>void copy(const Matrix& m,Eigen::Index row,Emit emit) const
    {
        for(Matrix::InnerIterator e(m,row);e;++e){
            require(e.col()>=0&&static_cast<std::size_t>(e.col())<toTarget.size()&&toTarget[e.col()]>=0,
                "certified stress block has an unbound actual native column");
            require(std::isfinite(e.value()),"nonfinite copied stress geometry coefficient");
            emit(toTarget[e.col()],e.value());}
    }
    StressGeometryRow row(const Matrix& m,Eigen::Index index) const
    {
        // The certified global face map is injective and the source CSR has
        // unique columns. Assign, never +=: no new rounded sum is introduced.
        StressGeometryRow result;copy(m,index,[&](int u,double value){result[u]=value;});return result;
    }
};

// Task41: transient representation routing only. The actual local-reuse
// certificate above remains responsible for geometry, source FP and columns.
// Row storage is already a sorted vector; no new fit/moment is cached here.
struct StressCsrRebindBudget
{
    static constexpr std::size_t limit=64u*1024u*1024u;
    std::size_t live=0,peak=0;
    bool claim(std::size_t bytes) noexcept {
        if(bytes>limit-live)return false;live+=bytes;return true;
    }
    void used() noexcept {peak=std::max(peak,live);}
    void release(std::size_t bytes) noexcept {live-=bytes;}
};
bool stressCsrRebindQualified() noexcept
{
#if defined(_MSC_VER) && defined(_M_X64) && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST) && EIGEN_WORLD_VERSION==3 && EIGEN_MAJOR_VERSION==4 && EIGEN_MINOR_VERSION==1 && !defined(EIGEN_SPARSEMATRIX_PLUGIN) && !defined(EIGEN_SPARSE_COMPRESSEDBASE_PLUGIN)
    return sizeof(long double)==sizeof(double)
        &&std::numeric_limits<long double>::digits==53;
#else
    return false;
#endif
}
bool stressCsrOrderedRow(const StressLocalReuse& reuse,const Matrix& source,Eigen::Index row) noexcept
{
    if(!source.isCompressed()||row<0||row>=source.rows())return false;
    int previous=-1;
    for(Matrix::InnerIterator e(source,row);e;++e){
        if(e.col()<0||static_cast<std::size_t>(e.col())>=reuse.toTarget.size())return false;
        const int target=reuse.toTarget[static_cast<std::size_t>(e.col())];
        if(target<=previous||static_cast<std::size_t>(target)>=reuse.oldFace.size())return false;previous=target;
        std::uint64_t bits=0;const double value=e.value();std::memcpy(&bits,&value,sizeof(bits));
        const auto magnitude=bits&UINT64_C(0x7fffffffffffffff);
        const auto exponent=magnitude&UINT64_C(0x7ff0000000000000);
        if((exponent==0&&magnitude!=0)||exponent==UINT64_C(0x7ff0000000000000))return false;
    }
    return true;
}
struct StressCsrCopyRows
{
    const StressLocalReuse* reuse=nullptr;const Matrix* source=nullptr;
    std::unique_ptr<Eigen::Index[]> sourceRows;
    std::size_t size=0,count=0,bytes=0;
    StressCsrRebindBudget* budget=nullptr;
    void reset() noexcept {
        if(!sourceRows&&!budget)return;
        StressPreparationScope scope;sourceRows.reset();
        if(budget)budget->release(bytes);
        reuse=nullptr;source=nullptr;size=0;count=0;bytes=0;budget=nullptr;
    }
    ~StressCsrCopyRows() noexcept {reset();}
    void prepare(std::size_t rows,const StressLocalReuse& input,const Matrix& matrix,
        StressCsrRebindBudget& memory,bool enabled) noexcept {
        if(!enabled||rows>(StressCsrRebindBudget::limit-256u)/sizeof(Eigen::Index))return;
        const auto requested=rows*sizeof(Eigen::Index)+256u;
        if(!memory.claim(requested))return;
        StressPreparationScope scope;
        if(!scope.captured){memory.release(requested);return;}
        try{
            auto storage=std::make_unique<Eigen::Index[]>(rows);
            std::fill_n(storage.get(),rows,Eigen::Index(-1));
            sourceRows=std::move(storage);size=rows;bytes=requested;
            reuse=&input;source=&matrix;budget=&memory;memory.used();
        }catch(const std::bad_alloc&){memory.release(requested);}
        catch(const std::length_error&){memory.release(requested);}
    }
    bool record(std::size_t row,Eigen::Index oldRow) noexcept {
        if(!sourceRows||row>=size||sourceRows[row]>=0||!stressCsrOrderedRow(*reuse,*source,oldRow))return false;
        sourceRows[row]=oldRow;++count;return true;
    }
    bool has(std::size_t row) const noexcept {return sourceRows&&row<size&&sourceRows[row]>=0;}
    template<class Emit>void copy(std::size_t row,Emit emit) const {
        reuse->copy(*source,sourceRows[row],emit); // Keep original column/finite guards.
    }
};
// Allocation is the only optional part. It completes before any original
// conversion/finite arithmetic. Failed allocation restores its entire entry.
// Destruction restores the CURRENT final state, not an earlier numeric state.
struct StressCsrPreparedMatrix
{
    std::unique_ptr<Matrix> matrix;StressCsrRebindBudget* budget=nullptr;std::size_t bytes=0;
    ~StressCsrPreparedMatrix() noexcept {
        if(!matrix&&!budget)return;
        StressPreparationScope scope;matrix.reset();if(budget)budget->release(bytes);
    }
    bool prepare(Eigen::Index rows,int cols,std::size_t nonzeros,StressCsrRebindBudget& memory,
        bool requiredPhysicalResult=false) noexcept {
        if(rows<0||cols<0||nonzeros>static_cast<std::size_t>(std::numeric_limits<Eigen::Index>::max()))return false;
        const auto rowCount=static_cast<std::size_t>(rows);
        constexpr auto slot=sizeof(double)+sizeof(Matrix::StorageIndex);
        const auto allocationLimit=requiredPhysicalResult?std::numeric_limits<std::size_t>::max():StressCsrRebindBudget::limit;
        if(rowCount>(allocationLimit-256u)/sizeof(Matrix::StorageIndex)-1)return false;
        const auto outer=(rowCount+1)*sizeof(Matrix::StorageIndex)+256u;
        if(nonzeros>(allocationLimit-outer)/slot)return false;
        // For Task49 this is the required fresh result, replacing the original
        // setFromTriplets output. Its coefficient payload is not extra scratch.
        // Charge the additional outer/object plus all separate literal/routing
        // workspace; Task41 callers retain the original full-storage accounting.
        const auto requested=requiredPhysicalResult?outer:outer+nonzeros*slot;
        if(!memory.claim(requested))return false;
        StressPreparationScope scope;
        if(!scope.captured){memory.release(requested);return false;}
        try{
            auto prepared=std::make_unique<Matrix>(rows,cols);
            prepared->reserve(static_cast<Eigen::Index>(nonzeros));
            if(prepared->data().allocatedSize()<0
                ||static_cast<std::size_t>(prepared->data().allocatedSize())>nonzeros){
                memory.release(requested);return false;
            }
            matrix=std::move(prepared);budget=&memory;bytes=requested;memory.used();return true;
        }catch(const std::bad_alloc&){memory.release(requested);return false;}
        catch(const std::length_error&){memory.release(requested);return false;}
    }
};

// Task49: merge original dirty triplets and disjoint certified source rows.
// This is a representation path only; no fit, moment or coefficient is cached.
bool stressCopiedTripletMatrix(Matrix& result,Eigen::Index rows,int columns,
    const std::vector<Eigen::Triplet<double,Eigen::Index>>& entries,
    const StressCsrCopyRows& copies,StressCsrRebindBudget& budget,std::size_t& directRows)
{
    if(rows<0||columns<0||!copies.count||!copies.sourceRows||copies.size!=static_cast<std::size_t>(rows))return false;
#if defined(EIGEN_SPARSE_COMPRESSED_STORAGE_REALLOCATE_PLUGIN) || defined(EIGEN_SCALAR_BINARY_OP_PLUGIN)
    return false; // The literal-workspace and duplicate arithmetic proof excludes custom hooks.
#endif
    StressReuseEnvironment entry;bool qualified=false;
    {StressPreparationScope observation;if(observation.captured)qualified=entry.read();}
    if(!qualified)return false;
    const int entryErrno=errno;
    try{
        // Omitted source rows and literal rows must be disjoint. Preserve every
        // dirty triplet and its original encounter order for Eigen's reduction.
        for(const auto& e:entries)
            if(e.row()<0||e.row()>=rows||e.col()<0||e.col()>=columns||copies.has(static_cast<std::size_t>(e.row())))
                throw std::length_error("optional raw/value row partition");
        std::size_t upper=entries.size();
        for(std::size_t row=0;row<copies.size;++row)if(copies.has(row)){
            const auto sourceRow=copies.sourceRows[row];
            const auto length=static_cast<std::size_t>(copies.source->outerIndexPtr()[sourceRow+1]
                -copies.source->outerIndexPtr()[sourceRow]);
            if(length>std::numeric_limits<std::size_t>::max()-upper)
                throw std::length_error("optional raw/value nonzero range");
            upper+=length;
        }
        // Bound the literal matrix plus Eigen3.4.1's opposite-major temporary,
        // duplicate index vector and assignment storage BEFORE entering Eigen.
        // Four coefficient slots/term and eight index arrays conservatively
        // cover its exact-size reserve/copy implementation (no custom plugins).
        StressCsrPreparedMatrix literal;
        {
            StressPreparationScope allocation;
            if(!allocation.captured)throw std::length_error("optional raw/value environment");
            constexpr std::size_t slot=sizeof(double)+sizeof(Matrix::StorageIndex);
            const auto n=static_cast<std::size_t>(rows),m=static_cast<std::size_t>(columns);
            std::size_t bytes=1024u;
            const auto take=[&](std::size_t count,std::size_t width){
                if(bytes>StressCsrRebindBudget::limit||count>(StressCsrRebindBudget::limit-bytes)/width)
                    throw std::length_error("optional raw/value literal workspace");
                bytes+=count*width;
            };
            take(n+1,8u*sizeof(Matrix::StorageIndex));take(m+1,8u*sizeof(Matrix::StorageIndex));
            take(entries.size(),4u*slot);
            if(!budget.claim(bytes))throw std::length_error("optional raw/value literal budget");
            literal.budget=&budget;literal.bytes=bytes;
            literal.matrix=std::make_unique<Matrix>(rows,columns);budget.used();
        }
        StressCsrPreparedMatrix prepared;
        if(!prepared.prepare(rows,columns,upper,budget,true))
            throw std::length_error("optional raw/value final workspace");
        // Same first-assignment/+= order for each literal coordinate. Do not
        // prune or finite-check here: those stay at the original full boundary.
        literal.matrix->setFromTriplets(entries.begin(),entries.end());
        Matrix& output=*prepared.matrix;
        for(Eigen::Index row=0;row<rows;++row){
            output.startVec(row);
            if(copies.has(static_cast<std::size_t>(row)))
                copies.copy(static_cast<std::size_t>(row),[&](int column,double value){output.insertBack(row,column)=value;});
            else for(Matrix::InnerIterator e(*literal.matrix,row);e;++e)
                output.insertBack(row,e.col())=e.value();
        }
        output.finalize();output.makeCompressed();
        if(errno!=entryErrno)throw std::length_error("optional raw/value errno changed");
        result.swap(output);directRows=copies.count;return true;
    }catch(...){
        // Private matrices unwind preserving their CURRENT FP state. Replay
        // only this representation, never dirty fits or a physical guard.
        entry.restore();errno=entryErrno;return false;
    }
}

bool stressSameDouble(double a,double b) noexcept
{return std::memcmp(&a,&b,sizeof(double))==0;}
bool stressSamePoint(Vector2D a,Vector2D b) noexcept
{return stressSameDouble(a.x,b.x)&&stressSameDouble(a.y,b.y);}
struct StressReuseFaceKey
{
    bool horizontal=false;Vector2D center{0,0};double tangent=0;long double area=0;double outward=0;
    bool negative=false,positive=false;
    auto words() const {return std::tie(horizontal,center.x,center.y,tangent,area,outward,negative,positive);}
    bool operator<(const StressReuseFaceKey& other) const {return words()<other.words();}
};
StressReuseFaceKey stressFaceKey(const NativeLeafFace& f)
{return {f.horizontal,f.center,f.tangentSpacing,f.area,f.outwardSign,f.negativeCell>=0,f.positiveCell>=0};}
bool stressSameFace(const NativeLeafFace& a,const NativeLeafFace& b)
{
    return a.horizontal==b.horizontal&&stressSamePoint(a.center,b.center)&&stressSameDouble(a.tangentSpacing,b.tangentSpacing)
        &&a.area==b.area&&a.outwardSign==b.outwardSign&&(a.negativeCell>=0)==(b.negativeCell>=0)&&(a.positiveCell>=0)==(b.positiveCell>=0);
}
std::vector<int> stressPatch(int p,const std::vector<std::vector<int>>& neighbors)
{
    std::vector<int> patch{p},frontier{p};
    for(int ring=0;ring<2;++ring){std::vector<int> next;
        for(int c:frontier)for(int n:neighbors[c])if(std::find(patch.begin(),patch.end(),n)==patch.end()){patch.push_back(n);next.push_back(n);}
        frontier=std::move(next);}
    return patch;
}
std::vector<int> stressSupport(const std::vector<int>& patch,const std::vector<std::vector<int>>& cellFaces)
{
    std::vector<int> support;for(int p:patch)support.insert(support.end(),cellFaces[p].begin(),cellFaces[p].end());
    std::sort(support.begin(),support.end());support.erase(std::unique(support.begin(),support.end()),support.end());return support;
}
// An overlap owns selected native groups, while validateFor authenticates only
// the two live meshes. Authenticate BOTH complete consumer selections before
// copying any block; equal selected counts do not establish equal membership.
// These integer/vector reads allocate nothing and perform no FP/errno work.
bool stressSelectedOverlapMatches(const std::vector<std::vector<std::size_t>>& groups,
    const std::vector<int>& selectedRows,int nativeCount) noexcept
{
    if(nativeCount<0||groups.size()!=static_cast<std::size_t>(nativeCount)
        ||selectedRows.size()!=static_cast<std::size_t>(nativeCount))return false;
    for(std::size_t native=0;native<groups.size();++native)
        if(groups[native].empty()!=(selectedRows[native]<0))return false;
    return true;
}
// The certificate checks the ORIGINAL ordered inputs, not a halo label/hash.
// Reflection, actual native constraints and each sampled functional all survive.
template<class SourceData>
StressLocalReuse prepareStressLocalReuse(const QuadtreeMesh& sourceMesh,const SourceData& source,
    const QuadtreeMesh& targetMesh,const HydroPhysicalDivergence& metrics,const LAMP::Mesh::NativeLeafOverlap& overlap,
    const std::vector<HydroStrainPoint>& points,const std::vector<Vector2D>& centers,const std::vector<Vector2D>& widths,
    const std::vector<NativeLeafFace>& faces,const std::vector<std::vector<int>>& cellFaces,
    const std::vector<std::vector<int>>& neighbors,const std::vector<std::vector<std::size_t>>& pointRows,
    HydroStressGeometryReuseReport& report,StressPreparationRecorder* recording,bool rt0GreenOnly)
{
    StressLocalReuse plan;
    const auto& old=source.metrics.layout();const auto& next=metrics.layout();
    const int count=next.pressureCount(),oldCount=old.pressureCount();
    StressReuseEnvironment current;
    if(!source.localReuseEligible
        ||(!rt0GreenOnly&&(!source.stressCompatible||!source.canonicalRaw))
        ||(rt0GreenOnly&&!source.stressCompatible&&!source.rt0GreenOnly)
        ||!current.read()||!current.same(source.localReuseEnvironment)
        ||!stressSelectedOverlapMatches(overlap.bySource(),old.domain().cellRowByIndex(),sourceMesh.cellIndexNum())
        ||!stressSelectedOverlapMatches(overlap.byTarget(),next.domain().cellRowByIndex(),targetMesh.cellIndexNum())
        ||source.components!=(next.domain().axisymmetric()?4:3))return plan;
    if(!stressSamePoint(sourceMesh.minCoord,targetMesh.minCoord)||!stressSamePoint(sourceMesh.maxCoord,targetMesh.maxCoord)
        ||old.domain().axisymmetric()!=next.domain().axisymmetric())return plan;
    plan.oldCell.assign(count,-1);plan.oldFace.assign(next.velocityCount(),-1);plan.toTarget.assign(old.velocityCount(),-1);
    plan.oldPoint.assign(points.size(),-1);plan.reuseRaw.assign(count,0);plan.reuseStress.assign(count,0);
    plan.reuseFace.assign(next.velocityCount(),0);plan.needFit.assign(count,0);
    plan.raw=source.canonicalRaw?&source.canonicalRaw->matrix:nullptr;plan.stress=&source.matrix;plan.volume=&source.volumeVelocity;
    plan.faceVelocity=&source.faceVelocity;plan.integrals=&source.componentIntegrals;plan.sourceFacePoints=&source.facePoints;
    // Original source mesh/overlap/full-FP gates above remain mandatory.
    // Private source preparation belongs to these same immutable Data owners.
    const auto* previous=source.preparation&&(rt0GreenOnly||source.preparation==source.canonicalRaw->preparation)
        &&source.preparation->matches(oldCount,old.velocityCount(),source.points.size())?source.preparation.get():nullptr;
    std::vector<Vector2D> fallbackCenters,fallbackWidths;
    std::vector<std::vector<std::size_t>> fallbackPoints;
    if(!previous){
        fallbackCenters.resize(oldCount);fallbackWidths.resize(oldCount);fallbackPoints.resize(oldCount);
        for(int p=0;p<oldCount;++p){const auto& c=sourceMesh.cells(sourceMesh.cellIndexToLabel[old.nativeCell(p)]);
            fallbackCenters[p]=c.centerCoord();fallbackWidths[p]=c.size();}
        for(std::size_t q=0;q<source.points.size();++q){const int p=source.points[q].cellDof;if(p<0||p>=oldCount)return StressLocalReuse{};fallbackPoints[p].push_back(q);}
    }
    const auto& oldCenters=previous?previous->centers:fallbackCenters;
    const auto& oldWidths=previous?previous->widths:fallbackWidths;
    const auto& oldPoints=previous?previous->pointRows:fallbackPoints;
    for(int p=0;p<count;++p){const int native=next.nativeCell(p);if(native<0||static_cast<std::size_t>(native)>=overlap.byTarget().size())continue;
        const auto& ids=overlap.byTarget()[native];if(ids.size()!=1)continue;const auto& piece=overlap.pieces()[ids[0]];
        if(piece.kind!=LAMP::Mesh::NativeLeafOverlapPiece::Kind::Unchanged)continue;
        const int prior=old.pressureDof(piece.sourceCell);if(prior<0)continue;
        if(!stressSamePoint(centers[p],oldCenters[prior])||!stressSamePoint(widths[p],oldWidths[prior])
            ||metrics.cellVolumes()[p]!=source.metrics.cellVolumes()[prior])continue;
        plan.oldCell[p]=prior;}
    std::vector<NativeLeafFace> fallbackFaces;
    std::vector<std::vector<int>> fallbackCellFaces,fallbackNeighbors;
    if(!previous){fallbackFaces.reserve(old.velocityCount());fallbackCellFaces.resize(oldCount);fallbackNeighbors.resize(oldCount);}
    const auto& oldFaces=previous?previous->faces:fallbackFaces;
    const auto& oldCellFaces=previous?previous->cellFaces:fallbackCellFaces;
    const auto& oldNeighbors=previous?previous->neighbors:fallbackNeighbors;
    std::map<StressReuseFaceKey,int> faceIndex;
    for(int u=0;u<old.velocityCount();++u){
        const auto f=previous?oldFaces[u]:old.domain().face(sourceMesh,old.nativeFace(u));
        if(!previous)fallbackFaces.push_back(f);
        // Keep the original floating-key insertion/comparison sequence.
        if(!faceIndex.emplace(stressFaceKey(f),u).second)return StressLocalReuse{};
        if(!previous){
            const int a=f.negativeCell<0?-1:old.pressureDof(f.negativeCell),b=f.positiveCell<0?-1:old.pressureDof(f.positiveCell);
            if(a>=0)fallbackCellFaces[a].push_back(u);if(b>=0)fallbackCellFaces[b].push_back(u);
            if(a>=0&&b>=0){fallbackNeighbors[a].push_back(b);fallbackNeighbors[b].push_back(a);}
        }
    }
    for(int u=0;u<next.velocityCount();++u){const auto found=faceIndex.find(stressFaceKey(faces[u]));if(found==faceIndex.end())continue;
        const int prior=found->second;if(!stressSameFace(faces[u],oldFaces[prior])
            ||metrics.faceMeasures()[u]!=source.metrics.faceMeasures()[prior])continue;
        if(plan.toTarget[prior]>=0)return StressLocalReuse{};plan.oldFace[u]=prior;plan.toTarget[prior]=u;}
    for(int p=0;p<count;++p){const int prior=plan.oldCell[p];if(prior<0)continue;
        bool same=oldPoints[prior].size()==4&&pointRows[p].size()==4;
        if(same)for(std::size_t k=0;k<4;++k){const auto a=oldPoints[prior][k],b=pointRows[p][k];
            if(!stressSamePoint(source.points[a].point,points[b].point)||!stressSameDouble(source.points[a].physicalWeight,points[b].physicalWeight)){same=false;break;}
            plan.oldPoint[b]=static_cast<Eigen::Index>(a);}
        if(!same){++report.quadratureRejectedCells;continue;}
        const auto patch=stressPatch(p,neighbors);
        const bool prepared=previous&&previous->complete[prior];
        const auto fallbackPatch=prepared?std::vector<int>{}:stressPatch(prior,oldNeighbors);
        const auto& oldPatch=prepared?previous->patches[prior]:fallbackPatch;
        if(patch.size()!=oldPatch.size()){++report.stencilRejectedCells;continue;}
        for(int c:patch)if(plan.oldCell[c]<0||std::find(oldPatch.begin(),oldPatch.end(),plan.oldCell[c])==oldPatch.end()){same=false;break;}
        const auto support=stressSupport(patch,cellFaces);
        const auto fallbackSupport=prepared?std::vector<int>{}:stressSupport(oldPatch,oldCellFaces);
        const auto& oldSupport=prepared?previous->supports[prior]:fallbackSupport;
        if(recording)recording->save(p,patch,support);
        if(support.size()!=oldSupport.size())same=false;
        if(same)for(std::size_t i=0;i<support.size();++i)if(plan.oldFace[support[i]]!=oldSupport[i]){same=false;break;}
        if(cellFaces[p].size()!=oldCellFaces[prior].size())same=false;
        if(same)for(std::size_t i=0;i<cellFaces[p].size();++i)if(plan.oldFace[cellFaces[p][i]]!=oldCellFaces[prior][i]){same=false;break;}
        // Same selected-patch incidence; an outside cell's index is not a fit datum.
        if(same)for(int u:support){const auto& a=faces[u];const auto& b=oldFaces[plan.oldFace[u]];
            for(int side=0;side<2;++side){const int c=side==0?a.negativeCell:a.positiveCell,oc=side==0?b.negativeCell:b.positiveCell;
                const int cp=c<0?-1:next.pressureDof(c),op=oc<0?-1:old.pressureDof(oc);
                const int selected=cp<0?-1:std::find(patch.begin(),patch.end(),cp)==patch.end()?-2:plan.oldCell[cp];
                const int oldSelected=op<0?-1:std::find(oldPatch.begin(),oldPatch.end(),op)==oldPatch.end()?-2:op;
                if(selected!=oldSelected){same=false;break;}}
            if(!same)break;}
        if(!same){++report.stencilRejectedCells;continue;}
        for(std::size_t q:pointRows[p]){const auto oq=plan.oldPoint[q];
            if(!rt0GreenOnly)for(int k=0;k<source.components;++k)if(!plan.columns(*plan.raw,static_cast<Eigen::Index>(oq)*source.components+k))same=false;
            for(int k=0;k<2;++k)if(!plan.columns(*plan.volume,static_cast<Eigen::Index>(oq)*2+k))same=false;}
        if(!same){++report.columnRejectedCells;continue;}
        plan.reuseRaw[p]=1;if(!rt0GreenOnly)++report.reusedRawCells;
    }
    plan.oldFacePoints.resize(old.velocityCount());
    for(std::size_t q=0;q<source.facePoints.size();++q){const int u=source.facePoints[q].velocityDof;
        if(u<0||u>=old.velocityCount())return StressLocalReuse{};plan.oldFacePoints[u].push_back(q);}
    const bool rz=next.domain().axisymmetric();const double gauss=1/std::sqrt(3.0);
    for(int u=0;u<next.velocityCount();++u){const int prior=plan.oldFace[u];if(prior<0||plan.oldFacePoints[prior].size()!=2)continue;
        const auto& f=faces[u];bool same=true;
        for(int side=0;side<2;++side){const int c=side==0?f.negativeCell:f.positiveCell;
            const int oc=side==0?oldFaces[prior].negativeCell:oldFaces[prior].positiveCell;
            const int p=c<0?-1:next.pressureDof(c),op=oc<0?-1:old.pressureDof(oc);
            if(p>=0&&(!plan.reuseRaw[p]||plan.oldCell[p]!=op))same=false;}
        for(int k=0;k<2&&same;++k){const double g=k==0?-gauss:gauss;auto point=f.center;
            if(f.horizontal)point.x+=g*f.tangentSpacing/2;else point.y+=g*f.tangentSpacing/2;
            const double weight=f.area/2*(rz&&f.horizontal?point.x/f.center.x:1);
            const auto oldQ=plan.oldFacePoints[prior][k];const auto& oldSample=source.facePoints[oldQ];
            if(!stressSamePoint(point,oldSample.point)||!stressSameDouble(weight,oldSample.physicalWeight)
                ||!plan.columns(*plan.faceVelocity,static_cast<Eigen::Index>(oldQ)*2)
                ||!plan.columns(*plan.faceVelocity,static_cast<Eigen::Index>(oldQ)*2+1))same=false;}
        if(same){plan.reuseFace[u]=1;++report.reusedFaces;}
    }
    for(int p=0;p<count;++p){bool same=plan.reuseRaw[p]!=0;const int prior=plan.oldCell[p];
        if(!rt0GreenOnly){
        if(same)for(int u:cellFaces[p])if(!plan.reuseFace[u])same=false;
        if(same){for(std::size_t q:pointRows[p])for(int k=0;k<source.components;++k)
                if(!plan.columns(*plan.stress,static_cast<Eigen::Index>(plan.oldPoint[q])*source.components+k))same=false;
            for(int k=0;k<source.components;++k)if(!plan.columns(*plan.integrals,static_cast<Eigen::Index>(prior)*source.components+k))same=false;
            if(!same)++report.columnRejectedCells;}
        if(same){plan.reuseStress[p]=1;++report.reusedStressCells;}
        }
        if(!plan.reuseRaw[p])plan.needFit[p]=1;
    }
    // Actual nonreused face endpoints, not a radius guess, define the trace halo.
    for(int u=0;u<next.velocityCount();++u)if(!plan.reuseFace[u])for(int cell:{faces[u].negativeCell,faces[u].positiveCell})
        if(cell>=0)plan.needFit[next.pressureDof(cell)]=1;
    report.recomputedFitCells=static_cast<std::size_t>(std::count(plan.needFit.begin(),plan.needFit.end(),static_cast<unsigned char>(1)));
    report.rebuiltFaces=static_cast<std::size_t>(next.velocityCount())-report.reusedFaces;
    report.enabled=true;plan.enabled=true;return plan;
}
}
HydroStrainOperator::HydroStrainOperator(std::shared_ptr<const Data> data) : data_(std::move(data)) {}
const HydroPhysicalDivergence& HydroStrainOperator::metrics() const noexcept {return data_->metrics;}
const HydroDofLayout& HydroStrainOperator::layout() const noexcept {return metrics().layout();}
const std::vector<HydroStrainPoint>& HydroStrainOperator::points() const noexcept {return data_->points;}
int HydroStrainOperator::components() const noexcept {return data_->components;}
const Matrix& HydroStrainOperator::matrix() const
{require(!data_->rt0GreenOnly,"RT0 Green geometry has no requested strain matrix");return data_->matrix;}
bool HydroStrainOperator::isStressCompatible() const noexcept {return data_->stressCompatible;}
bool HydroStrainOperator::isRt0GreenCompatible() const noexcept {return data_->stressCompatible||data_->rt0GreenOnly;}
const HydroStressGeometryReuseReport& HydroStrainOperator::geometryReuseReport() const noexcept {return data_->reuse;}
bool HydroStrainOperator::hasCanonicalRawStrain() const noexcept
{return data_->stressCompatible&&bool(data_->canonicalRaw);}
HydroStrainOperator HydroStrainOperator::canonicalRawStrain() const
{
    require(hasCanonicalRawStrain(),"stress geometry has no builtin canonical raw strain");
    return HydroStrainOperator(data_->canonicalRaw);
}
const std::vector<HydroStressFacePoint>& HydroStrainOperator::stressFacePoints() const
{require(isRt0GreenCompatible(),"raw strain has no stress-compatible face trace");return data_->facePoints;}
const Matrix& HydroStrainOperator::stressFaceVelocityMatrix() const
{require(isRt0GreenCompatible(),"raw strain has no stress-compatible face trace");return data_->faceVelocity;}
const Matrix& HydroStrainOperator::stressVolumeVelocityMatrix() const
{require(isRt0GreenCompatible(),"raw strain has no stress-compatible volume value map");return data_->volumeVelocity;}
const Matrix& HydroStrainOperator::stressComponentIntegralMatrix() const
{require(data_->stressCompatible,"raw strain has no stress-compatible component integrals");return data_->componentIntegrals;}
const HydroStressCompatibilityReport& HydroStrainOperator::stressCompatibilityReport() const
{require(data_->stressCompatible,"raw strain has no stress compatibility report");return data_->compatibility;}
std::vector<HydroStrainPoint> HydroStrainOperator::volumeQuadrature(const QuadtreeMesh& mesh,
    const HydroPhysicalDivergence& metrics)
{
    const auto& layout=metrics.layout();layout.validateMesh(mesh);
    std::vector<HydroStrainPoint> result;result.reserve(static_cast<std::size_t>(layout.pressureCount())*4);
    const double g=1/std::sqrt(3.0);
    for (int p=0;p<layout.pressureCount();++p) {
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[layout.nativeCell(p)]);
        const auto c=cell.centerCoord(),w=cell.size();
        for (double x : {-g,g}) for (double y : {-g,g}) {
            const Vector2D point{c.x+x*w.x/2,c.y+y*w.y/2};
            const double weight=metrics.cellVolumes()[p]/4*(layout.domain().axisymmetric() ? point.x/c.x : 1);
            require(std::isfinite(weight) && weight>0,"invalid physical volume quadrature");
            result.push_back({p,point,weight});
        }
    }
    return result;
}
HydroStrainOperator HydroStrainOperator::build(const QuadtreeMesh& mesh,HydroPhysicalDivergence metrics,
    const std::vector<HydroStrainPoint>& points)
{
    return buildImpl(mesh,std::move(metrics),points,false);
}
const std::vector<int>& HydroStrainOperator::stationaryNoSlipVelocityDofs() const noexcept
{return data_->stationaryNoSlipWalls;}
namespace {
using NoSlipTriplet=Eigen::Triplet<double,Eigen::Index>;
struct NoSlipFitCounts {
    std::size_t requests=0,hits=0,qr=0,svd=0,bankPeak=0;
    void add(const NoSlipFitCounts& x) noexcept
    {requests+=x.requests;hits+=x.hits;qr+=x.qr;svd+=x.svd;bankPeak+=x.bankPeak;}
};
struct NoSlipCellEntries {std::vector<NoSlipTriplet> values;};
// Eigen's two input passes see precisely the old p/axis/q/component/j order.
// No second complete flat triplet payload is allocated.
class NoSlipTripletIterator {
    const std::vector<NoSlipCellEntries>* cells_=nullptr;
    std::size_t cell_=0,entry_=0;
    void advance() noexcept {
        while(cell_<cells_->size()&&entry_==(*cells_)[cell_].values.size()){++cell_;entry_=0;}
    }
public:
    using iterator_category=std::forward_iterator_tag;
    using value_type=NoSlipTriplet;using difference_type=std::ptrdiff_t;
    using pointer=const NoSlipTriplet*;using reference=const NoSlipTriplet&;
    NoSlipTripletIterator() noexcept=default;
    NoSlipTripletIterator(const std::vector<NoSlipCellEntries>& cells,bool end) noexcept
        :cells_(&cells),cell_(end?cells.size():0){advance();}
    reference operator*() const noexcept {return (*cells_)[cell_].values[entry_];}
    pointer operator->() const noexcept {return &**this;}
    NoSlipTripletIterator& operator++() noexcept {++entry_;advance();return *this;}
    NoSlipTripletIterator operator++(int) noexcept {auto previous=*this;++*this;return previous;}
    bool operator==(const NoSlipTripletIterator& b) const noexcept
    {return cells_==b.cells_&&cell_==b.cell_&&entry_==b.entry_;}
    bool operator!=(const NoSlipTripletIterator& b) const noexcept {return !(*this==b);}
};
struct NoSlipAssemblyBudget {
    static constexpr std::size_t limit=256u*1024u*1024u;
    std::atomic<std::size_t> bytes{0};
    void add(std::size_t count,std::size_t width=1) {
        if(width&&count>limit/width)throw std::length_error("optional no-slip assembly capacity");
        const auto amount=count*width;auto previous=bytes.load(std::memory_order_relaxed);
        do{if(previous>limit||amount>limit-previous)throw std::length_error("optional no-slip assembly capacity");}
        while(!bytes.compare_exchange_weak(previous,previous+amount,std::memory_order_relaxed));
    }
};
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE) \
    && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST)
#pragma float_control(precise,on,push)
#pragma float_control(except,on)
#pragma fenv_access(on)
struct NoSlipOptionalState {
    StrainFloatingState state;int error=errno;
    NoSlipOptionalState() noexcept {state.capture();errno=error;}
    ~NoSlipOptionalState() noexcept {state.restore();errno=error;}
};
struct NoSlipFitRequest {
    FitTemplates::Key key;
    FitTemplates::Environment environment;
    int error=0;bool qualified=false;
    ~NoSlipFitRequest() noexcept {if(key.capacity()){NoSlipOptionalState current;FitTemplates::Key().swap(key);}}
};
struct NoSlipFitBank {
    const std::size_t limit;
    std::map<FitTemplates::Key,Eigen::MatrixXd> fits;
    std::size_t bytes=sizeof(NoSlipFitBank);
    explicit NoSlipFitBank(std::size_t capacity) noexcept:limit(capacity){}
    ~NoSlipFitBank() noexcept {NoSlipOptionalState current;fits.clear();}
};
bool noSlipQualified(StrainFloatingState& caller) noexcept {
    if(omp_in_parallel()||!caller.capture())return false;
    unsigned control=0;
    if(_controlfp_s(&control,0,0)!=0||(control&_MCW_EM)!=_MCW_EM
        ||(caller.csr&0x1f80u)!=0x1f80u)return false;
#if defined(_EM_AMBIGUOUS)
    if(control&_EM_AMBIGUOUS)return false;
#endif
    return true;
}
template<class Trace> __declspec(noinline) bool noSlipFitLookup(NoSlipFitBank* bank,const Eigen::MatrixXd& a,
    const std::vector<int>& columns,bool trueAxis,int component,bool wall,Trace&& trace,
    NoSlipFitRequest& request,Eigen::MatrixXd& result) noexcept {
    if(!bank)return false;
    NoSlipOptionalState current;
    if(!current.state.valid)return false;
    try{
        request.error=errno;
        if(!FitTemplates::environment(request.environment))return false;
        const auto rows=a.rows(),cols=a.cols();const auto count=columns.size();
        FitTemplates::append(request.key,rows);FitTemplates::append(request.key,cols);
        FitTemplates::append(request.key,trueAxis);FitTemplates::append(request.key,component);
        FitTemplates::append(request.key,wall);FitTemplates::append(request.key,count);
        FitTemplates::append(request.key,columns.data(),count*sizeof(int));
        FitTemplates::append(request.key,a.data(),static_cast<std::size_t>(a.size())*sizeof(double));
        if(wall){const Eigen::MatrixXd actualTrace=trace();const auto n=actualTrace.rows(),m=actualTrace.cols();
            FitTemplates::append(request.key,n);FitTemplates::append(request.key,m);
            FitTemplates::append(request.key,actualTrace.data(),static_cast<std::size_t>(actualTrace.size())*sizeof(double));}
        FitTemplates::append(request.key,request.environment.bytes.data(),request.environment.bytes.size());
        FitTemplates::append(request.key,request.environment.csr);FitTemplates::append(request.key,request.error);
        request.qualified=true;
        const auto found=bank->fits.find(request.key);
        if(found==bank->fits.end())return false;
        result=found->second;return true;
    }catch(...){request.qualified=false;return false;}
}
__declspec(noinline) void noSlipFitRemember(NoSlipFitBank* bank,NoSlipFitRequest& request,const Eigen::MatrixXd& coefficients) noexcept {
    if(!bank||!request.qualified)return;
    NoSlipOptionalState current;
    if(!current.state.valid)return;
    try{
        FitTemplates::Environment after;
        if(errno!=request.error||!FitTemplates::environment(after)
            ||!FitTemplates::same(request.environment,after)||!coefficients.allFinite())return;
        // Account exact key capacity, owned coefficient elements, object/node
        // metadata and a conservative per-allocation/control allowance.
        const auto keyBytes=request.key.capacity();const auto values=static_cast<std::size_t>(coefficients.size());
        constexpr auto metadata=sizeof(std::pair<const FitTemplates::Key,Eigen::MatrixXd>)+256u;
        if(bank->bytes>bank->limit||metadata>bank->limit-bank->bytes)return;
        auto remaining=bank->limit-bank->bytes-metadata;
        if(keyBytes>remaining)return;remaining-=keyBytes;
        if(values>remaining/sizeof(double))return;
        const auto extra=metadata+keyBytes+values*sizeof(double);
        const auto inserted=bank->fits.emplace(std::move(request.key),coefficients);
        if(inserted.second)bank->bytes+=extra;
    }catch(...){} // Optional math reuse never creates a physical failure.
}
#if !defined(_M_FP_STRICT) || !_M_FP_STRICT
#pragma fenv_access(off)
#endif
#pragma float_control(pop)
#else
struct NoSlipFitBank {};
struct NoSlipFitRequest {};
template<class Trace> bool noSlipFitLookup(NoSlipFitBank*,const Eigen::MatrixXd&,
    const std::vector<int>&,bool,int,bool,Trace&&,NoSlipFitRequest&,Eigen::MatrixXd&) noexcept {return false;}
void noSlipFitRemember(NoSlipFitBank*,NoSlipFitRequest&,const Eigen::MatrixXd&) noexcept {}
#endif
}
HydroStrainOperator HydroStrainOperator::buildStationaryNoSlip(const QuadtreeMesh& mesh,
    HydroPhysicalDivergence metrics,const std::vector<HydroStrainPoint>& points,
    const std::vector<int>& walls)
{
    return buildStationaryNoSlipImpl(mesh,std::move(metrics),points,walls,nullptr);
}
HydroStrainOperator HydroStrainOperator::buildStationaryNoSlipImpl(const QuadtreeMesh& mesh,
    HydroPhysicalDivergence metrics,const std::vector<HydroStrainPoint>& points,
    const std::vector<int>& wallVelocityDofs,HydroNoSlipPreparationReport* report)
{
    metrics.layout().validateMesh(mesh);
    require(!points.empty()&&!wallVelocityDofs.empty(),"stationary no-slip requires actual points and wall rows");
    const auto coordinate=[](Vector2D x,int axis){return axis==0?x.x:x.y;};
    auto data=std::make_shared<Data>(std::move(metrics));data->points=points;
    const auto& layout=data->metrics.layout();const bool rz=layout.domain().axisymmetric();data->components=rz?4:3;
    auto& walls=data->stationaryNoSlipWalls;walls=wallVelocityDofs;std::sort(walls.begin(),walls.end());
    require(std::adjacent_find(walls.begin(),walls.end())==walls.end(),"duplicate stationary no-slip wall");
    std::vector<NativeLeafFace> faces;faces.reserve(layout.velocityCount());
    std::vector<std::vector<int>> cellFaces(layout.pressureCount()),neighbors(layout.pressureCount()),cellWalls(layout.pressureCount());
    for(int u=0;u<layout.velocityCount();++u){
        const auto f=layout.domain().face(mesh,layout.nativeFace(u));faces.push_back(f);
        const int a=f.negativeCell<0?-1:layout.pressureDof(f.negativeCell),b=f.positiveCell<0?-1:layout.pressureDof(f.positiveCell);
        if(a>=0)cellFaces[a].push_back(u);if(b>=0)cellFaces[b].push_back(u);
        if(a>=0&&b>=0){neighbors[a].push_back(b);neighbors[b].push_back(a);}
    }
    for(int u:walls){
        require(u>=0&&u<layout.velocityCount(),"stationary no-slip wall DOF out of range");const auto& f=faces[u];
        require(f.outwardSign!=0&&f.area>0,"stationary no-slip requires positive-area physical wall, not interior/axis");
        const int native=f.negativeCell>=0?f.negativeCell:f.positiveCell,p=layout.pressureDof(native);
        require(p>=0,"stationary no-slip has no active owner");cellWalls[p].push_back(u);
    }
    std::vector<std::vector<std::size_t>> byCell(layout.pressureCount());
    for(std::size_t q=0;q<points.size();++q){const auto& x=points[q];const int p=x.cellDof;
        require(p>=0&&p<layout.pressureCount()&&std::isfinite(x.point.x)&&std::isfinite(x.point.y)
            &&std::isfinite(x.physicalWeight)&&x.physicalWeight>=0,"invalid no-slip volume sample");
        const auto& c=mesh.cells(mesh.cellIndexToLabel[layout.nativeCell(p)]);const auto center=c.centerCoord(),h=c.size();
        require(x.point.x>=center.x-h.x/2&&x.point.x<=center.x+h.x/2&&x.point.y>=center.y-h.y/2&&x.point.y<=center.y+h.y/2,
            "no-slip volume sample outside active leaf");byCell[p].push_back(q);
    }
    const auto cellKernel=[&](int p,std::vector<NoSlipTriplet>& entries,NoSlipFitBank* bank,
        NoSlipFitCounts& counts,NoSlipAssemblyBudget* budget,std::size_t* workspace){
        if(byCell[p].empty())return;
        std::vector<int> patch{p},frontier{p};
        for(int ring=0;ring<2;++ring){std::vector<int> next;for(int c:frontier)for(int n:neighbors[c])
            if(std::find(patch.begin(),patch.end(),n)==patch.end()){patch.push_back(n);next.push_back(n);}frontier=std::move(next);}
        std::vector<int> support;for(int c:patch)support.insert(support.end(),cellFaces[c].begin(),cellFaces[c].end());
        std::sort(support.begin(),support.end());support.erase(std::unique(support.begin(),support.end()),support.end());
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[layout.nativeCell(p)]);const auto h=cell.size();
        const bool axis=rz&&std::any_of(support.begin(),support.end(),[&](int u){return faces[u].area==0;});
        const auto center=axis?Vector2D{0,cell.centerCoord().y}:cell.centerCoord();
        // A no-slip line is a full trace, not an extra velocity unknown. Split
        // native wall subfaces on this cell must lie on the same actual line.
        if(!cellWalls[p].empty()){const auto& first=faces[cellWalls[p][0]];
            for(int u:cellWalls[p])require(faces[u].horizontal==first.horizontal
                &&coordinate(faces[u].center,first.horizontal?1:0)==coordinate(first.center,first.horizontal?1:0),
                "stationary no-slip P2 rejects a genuine multiwall corner");}
        if(budget){
            const auto n0=static_cast<std::size_t>(std::count_if(support.begin(),support.end(),[&](int u){return !faces[u].horizontal;}));
            const auto n1=support.size()-n0;
            const auto perPoint=(rz?3u:2u)*n0+2u*n1;
            if(perPoint&&byCell[p].size()>entries.max_size()/perPoint)throw std::length_error("optional no-slip triplet range");
            const auto required=byCell[p].size()*perPoint;
            budget->add(required,sizeof(NoSlipTriplet));
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE) && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST)
            NoSlipOptionalState current;
#endif
            entries.reserve(required);
            if(entries.capacity()>required)budget->add(entries.capacity()-required,sizeof(NoSlipTriplet));
            // Conservative Eigen QR/SVD/identity/solve workspace from this
            // actual support dimension, plus actual local integer capacities.
            const auto side=support.size()+8u;
            if(side>NoSlipAssemblyBudget::limit/side/512u)throw std::length_error("optional no-slip fit workspace");
            std::size_t work=0;
            const auto addWork=[&](std::size_t count,std::size_t width){
                if(count>(NoSlipAssemblyBudget::limit-work)/width)throw std::length_error("optional no-slip fit workspace");
                work+=count*width;
            };
            addWork(side*side,512u);addWork(patch.capacity(),sizeof(int));
            addWork(frontier.capacity(),sizeof(int));addWork(support.capacity(),sizeof(int));addWork(4096u,1u);
            if(work>*workspace){budget->add(work-*workspace);*workspace=work;}
        }
        for(int componentAxis=0;componentAxis<2;++componentAxis){
            std::vector<int> columns=axis?(componentAxis==0?std::vector<int>{1,4}:std::vector<int>{0,2,3,5}):std::vector<int>{0,1,2,3,4,5};
            std::vector<int> samples;for(int u:support)if((faces[u].horizontal?1:0)==componentAxis)samples.push_back(u);
            Eigen::MatrixXd a(samples.size(),columns.size());
            for(std::size_t j=0;j<samples.size();++j){const auto& f=faces[samples[j]];
                const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
                const auto v0=mesh.verts(edge.edgeToVert[0]).coord(),v1=mesh.verts(edge.edgeToVert[1]).coord();
                const double length=f.horizontal?std::abs(v1.x-v0.x):std::abs(v1.y-v0.y);
                const double x=(f.center.x-center.x)/h.x,y=(f.center.y-center.y)/h.y;
                Eigen::Matrix<double,1,6> mean=basis(x,y);
                if(f.horizontal){const double variance=length*length/(12*h.x*h.x);
                    if(rz){require(f.center.x>0,"horizontal no-slip mean at nonpositive radius");
                        const double shift=length*length/(12*f.center.x*h.x);mean[1]+=shift;mean[3]+=2*x*shift+variance;mean[4]+=y*shift;}
                    else mean[3]+=variance;
                }else mean[5]+=length*length/(12*h.y*h.y);
                for(std::size_t k=0;k<columns.size();++k)a(j,k)=mean[columns[k]];
            }
            const auto makeTrace=[&]{
                const auto& f=faces[cellWalls[p][0]];Eigen::MatrixXd trace(3,columns.size());int j=0;
                for(double offset:{-.5,0.,.5}){auto point=f.center;
                    if(f.horizontal){point.x=cell.centerCoord().x+offset*h.x;point.y=f.center.y;}
                    else{point.x=f.center.x;point.y=cell.centerCoord().y+offset*h.y;}
                    const auto v=basis((point.x-center.x)/h.x,(point.y-center.y)/h.y);
                    for(std::size_t k=0;k<columns.size();++k)trace(j,k)=v[columns[k]];++j;
                }
                return trace;
            };
            ++counts.requests;NoSlipFitRequest request;Eigen::MatrixXd coefficients;
            if(noSlipFitLookup(bank,a,columns,axis,componentAxis,!cellWalls[p].empty(),makeTrace,request,coefficients))++counts.hits;
            else{
                const auto literalFit=[&]{
            ++counts.qr;Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(a);qr.setThreshold(1e-12);
            require(qr.rank()==a.cols(),"rank-deficient actual face-mean no-slip P2 support");
            Eigen::MatrixXd coefficients;
            if(cellWalls[p].empty())coefficients=qr.solve(Eigen::MatrixXd::Identity(a.rows(),a.rows()));
            else{
                // Strain is a least-squares differential of actual face means.
                // Its physical Dirichlet trace is hard; native face means are
                // samples, not redundant exact constraints on split AMR faces.
                const Eigen::MatrixXd trace=makeTrace();
                ++counts.svd;Eigen::JacobiSVD<Eigen::MatrixXd> svd(trace,Eigen::ComputeThinU|Eigen::ComputeFullV);svd.setThreshold(1e-12);
                const Eigen::MatrixXd nullspace=svd.matrixV().rightCols(trace.cols()-svd.rank());
                require(nullspace.cols()>0,"no-slip P2/true-axis constraints leave no velocity space; refine actual support");
                const Eigen::MatrixXd reduced=a*nullspace;++counts.qr;Eigen::ColPivHouseholderQR<Eigen::MatrixXd> freeFit(reduced);freeFit.setThreshold(1e-12);
                require(freeFit.rank()==nullspace.cols(),"rank-deficient constrained no-slip fit");
                coefficients=nullspace*freeFit.solve(Eigen::MatrixXd::Identity(a.rows(),a.rows()));
                require(coefficients.allFinite()&&(trace*coefficients).cwiseAbs().maxCoeff()<=1e-12,
                    "inconsistent stationary no-slip homogeneous trace");
            }
            return coefficients;
                };
                coefficients=literalFit(); // QR/SVD temporaries are gone before quiet admission.
                noSlipFitRemember(bank,request,coefficients);
            }
            for(std::size_t q:byCell[p]){const auto point=points[q].point;const double x=(point.x-center.x)/h.x,y=(point.y-center.y)/h.y;
                Eigen::Matrix<double,1,6> dx,dy;dx<<0,1/h.x,0,2*x/h.x,y/h.x,0;dy<<0,0,1/h.y,0,x/h.y,2*y/h.y;
                for(int k=0;k<data->components;++k){Eigen::Matrix<double,1,6> functional=Eigen::Matrix<double,1,6>::Zero();
                    if(k==0&&componentAxis==0)functional=dx;if(k==1&&componentAxis==1)functional=dy;
                    if(rz&&k==2&&componentAxis==0){require(point.x>0||axis,"missing true-axis no-slip parity");functional=point.x==0?dx:basis(x,y)/point.x;}
                    if(k==data->components-1)functional=(componentAxis==0?dy:dx)/std::sqrt(2.);
                    Eigen::RowVectorXd restricted(columns.size());for(std::size_t j=0;j<columns.size();++j)restricted[j]=functional[columns[j]];
                    const Eigen::RowVectorXd row=restricted*coefficients;require(row.allFinite(),"nonfinite stationary no-slip strain");
                    for(Eigen::Index j=0;j<row.size();++j)if(row[j]!=0)entries.emplace_back(static_cast<Eigen::Index>(q)*data->components+k,samples[j],row[j]);
                }
            }
        }
    };
    std::vector<NoSlipTriplet> entries;NoSlipFitCounts totals;
    bool completed=false;
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE) \
    && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST)
    const int callerError=errno;StrainFloatingState caller;
    const bool qualified=noSlipQualified(caller);errno=callerError;
    if(qualified){
        std::vector<NoSlipFitCounts> attempts;int actualWorkers=0;
        try{
            const auto workCells=static_cast<int>(std::count_if(byCell.begin(),byCell.end(),[](const auto& cell){return !cell.empty();}));
            const int team=std::min(std::max(1,omp_get_max_threads()),workCells);
            const auto bankAllowance=std::min<std::size_t>(4u*1024u*1024u,32u*1024u*1024u/static_cast<std::size_t>(team));
            if(bankAllowance<sizeof(NoSlipFitBank)+256u)throw std::length_error("optional no-slip bank metadata capacity");
            NoSlipAssemblyBudget budget;
            budget.add(team,bankAllowance);
            budget.add(1,sizeof(budget)+3u*sizeof(std::vector<std::size_t>)+512u);
            budget.add(team,sizeof(NoSlipFitCounts)+sizeof(std::unique_ptr<NoSlipFitBank>)+sizeof(std::size_t)+4096u);
            attempts.resize(static_cast<std::size_t>(team));
            if(attempts.capacity()>attempts.size())budget.add(attempts.capacity()-attempts.size(),sizeof(NoSlipFitCounts));
            budget.add(static_cast<std::size_t>(layout.pressureCount()),sizeof(NoSlipCellEntries));
            std::vector<NoSlipCellEntries> assembled(static_cast<std::size_t>(layout.pressureCount()));
            if(assembled.capacity()>assembled.size())budget.add(assembled.capacity()-assembled.size(),sizeof(NoSlipCellEntries));
            std::vector<std::unique_ptr<NoSlipFitBank>> banks(static_cast<std::size_t>(team));
            if(banks.capacity()>banks.size())budget.add(banks.capacity()-banks.size(),sizeof(std::unique_ptr<NoSlipFitBank>));
            for(int worker=0;worker<team;++worker){banks[worker]=std::make_unique<NoSlipFitBank>(bankAllowance);attempts[worker].bankPeak=banks[worker]->bytes;}
            std::vector<std::size_t> workspace(static_cast<std::size_t>(team));
            if(workspace.capacity()>workspace.size())budget.add(workspace.capacity()-workspace.size(),sizeof(std::size_t));
            std::atomic<bool> failed{false};unsigned long raisedNative=0;unsigned raisedCsr=0;
            caller.restore();errno=callerError;
#pragma omp parallel num_threads(team)
            {
                const int workerId=omp_get_thread_num(),workerError=errno;
                StrainFloatingState worker;const bool captured=worker.capture();
                const bool ready=captured&&worker.enterFresh(caller);errno=callerError;
                if(!ready)failed.store(true,std::memory_order_relaxed);
#pragma omp single
                actualWorkers=omp_get_num_threads();
#pragma omp for schedule(static) nowait
                for(int p=0;p<layout.pressureCount();++p){
                    if(!ready||failed.load(std::memory_order_relaxed))continue;
                    try{cellKernel(p,assembled[static_cast<std::size_t>(p)].values,banks[workerId].get(),
                        attempts[workerId],&budget,&workspace[workerId]);}
                    catch(...){failed.store(true,std::memory_order_relaxed);}
                }
                if(ready){
                    std::fenv_t after{};const bool got=std::fegetenv(&after)==0;
                    if(!got||errno!=callerError)failed.store(true,std::memory_order_relaxed);
                    const unsigned flags=_mm_getcsr()&0x3fu;
#pragma omp critical(lamp_noslip_floating_flags)
                    {if(got)raisedNative|=after._Fe_stat;raisedCsr|=flags;}
                }
                if(captured)worker.restore();errno=workerError;
            }
            for(int worker=0;worker<team;++worker)attempts[worker].bankPeak=banks[worker]->bytes;
            if(report){report->maximumWorkers=std::max(report->maximumWorkers,actualWorkers);
                report->parallelUsed=report->parallelUsed||actualWorkers>1;}
            if(failed.load(std::memory_order_relaxed))throw std::length_error("optional no-slip worker replay");
            auto merged=caller.native;merged._Fe_stat|=raisedNative;
            if(std::fesetenv(&merged)!=0)std::terminate();_mm_setcsr(caller.csr|raisedCsr);errno=callerError;
            data->matrix.resize(static_cast<Eigen::Index>(points.size())*data->components,layout.velocityCount());
            data->matrix.setFromTriplets(NoSlipTripletIterator(assembled,false),NoSlipTripletIterator(assembled,true));
            data->matrix.makeCompressed();
            {NoSlipOptionalState current;
                std::vector<NoSlipCellEntries>().swap(assembled);
                for(auto& bank:banks)bank.reset();
                std::vector<std::unique_ptr<NoSlipFitBank>>().swap(banks);
                std::vector<std::size_t>().swap(workspace);}
            completed=true;
        }catch(...){
            // All optional banks/cell buffers have unwound. Restore the exact
            // original cell-loop entry before reproducing its first error.
            {NoSlipOptionalState current;Matrix().swap(data->matrix);}
            caller.restore();errno=callerError;
            if(report){report->parallelFallback=true;report->maximumWorkers=std::max(report->maximumWorkers,actualWorkers);
                report->parallelUsed=report->parallelUsed||actualWorkers>1;}
        }
        for(const auto& count:attempts)totals.add(count);
        {NoSlipOptionalState current;std::vector<NoSlipFitCounts>().swap(attempts);}
    }else {if(caller.valid)caller.restore();errno=callerError;}
#endif
    if(!completed){
        if(report)report->maximumWorkers=std::max(report->maximumWorkers,1);
        // The literal reference/replay never consults or publishes a bank.
        for(int p=0;p<layout.pressureCount();++p)cellKernel(p,entries,nullptr,totals,nullptr,nullptr);
        data->matrix.resize(static_cast<Eigen::Index>(points.size())*data->components,layout.velocityCount());
        data->matrix.setFromTriplets(entries.begin(),entries.end());data->matrix.makeCompressed();
    }
    if(report){report->fitRequests+=totals.requests;report->fitHits+=totals.hits;
        report->qrFactorizations+=totals.qr;report->svdFactorizations+=totals.svd;
        report->bankPeakBytes=std::max(report->bankPeakBytes,totals.bankPeak);}
    return HydroStrainOperator(std::move(data));
}

HydroStrainOperator HydroStrainOperator::buildRt0Green(const QuadtreeMesh& mesh,HydroPhysicalDivergence metrics)
{
    const auto points=volumeQuadrature(mesh,metrics);
    return buildImpl(mesh,std::move(metrics),points,true,nullptr,false,nullptr,nullptr,nullptr,true);
}
HydroStrainOperator HydroStrainOperator::buildRt0GreenRemapped(const QuadtreeMesh& sourceMesh,
    const HydroStrainOperator& source,const QuadtreeMesh& targetMesh,HydroPhysicalDivergence metrics,
    const LAMP::Mesh::NativeLeafOverlap& overlap)
{
    source.layout().validateMesh(sourceMesh);overlap.validateFor(sourceMesh,targetMesh);
    require(source.layout().domain().axisymmetric()==metrics.layout().domain().axisymmetric(),
        "remapped stress source/target physical measure differs");
    const auto points=volumeQuadrature(targetMesh,metrics);
    StressReuseEnvironment entry;const bool captured=entry.read();const int savedErrno=errno;
    try{return buildImpl(targetMesh,metrics,points,true,nullptr,false,&sourceMesh,&source,&overlap,true);}
    catch(const StressLocalReuseMismatch&){
        if(!captured)throw;entry.restore();errno=savedErrno;
        return buildImpl(targetMesh,std::move(metrics),points,true,nullptr,false,nullptr,&source,nullptr,true);
    }
}
HydroStrainOperator HydroStrainOperator::buildStressCompatible(const QuadtreeMesh& mesh,HydroPhysicalDivergence metrics)
{
    const auto points=volumeQuadrature(mesh,metrics);
    return buildImpl(mesh,std::move(metrics),points,true,nullptr,true);
}
HydroStrainOperator HydroStrainOperator::buildStressCompatibleRemapped(const QuadtreeMesh& sourceMesh,
    const HydroStrainOperator& source,const QuadtreeMesh& targetMesh,HydroPhysicalDivergence metrics,
    const LAMP::Mesh::NativeLeafOverlap& overlap)
{
    source.layout().validateMesh(sourceMesh);overlap.validateFor(sourceMesh,targetMesh);
    require(source.layout().domain().axisymmetric()==metrics.layout().domain().axisymmetric(),
        "remapped stress source/target physical measure differs");
    const auto points=volumeQuadrature(targetMesh,metrics);
    StressReuseEnvironment entry;const bool captured=entry.read();const int savedErrno=errno;
    try{return buildImpl(targetMesh,metrics,points,true,nullptr,true,&sourceMesh,&source,&overlap);}
    catch(const StressLocalReuseMismatch&){
        // Only this optional certificate may select a literal cold replay.
        // Original numerical/owner exceptions are never caught here.
        if(!captured)throw;entry.restore();errno=savedErrno;
        return buildImpl(targetMesh,std::move(metrics),points,true,nullptr,true,nullptr,&source,nullptr);
    }
}
HydroStrainOperator HydroStrainOperator::buildStressCompatible(const QuadtreeMesh& mesh,
    HydroPhysicalDivergence metrics,const std::vector<HydroStrainPoint>& points,
    const std::vector<HydroStressFacePoint>& faces)
{
    require(!faces.empty(),"empty partitioned stress face quadrature");
    auto result=buildImpl(mesh,std::move(metrics),points,true,&faces);
    require(result.stressCompatibilityReport().maximumRelativeQuadratureWeightDefect<=1e-12,
        "partitioned stress cell quadrature does not have the actual full native measure");
    return result;
}
HydroStrainOperator HydroStrainOperator::buildImpl(const QuadtreeMesh& mesh,HydroPhysicalDivergence metrics,
    const std::vector<HydroStrainPoint>& points,bool stressCompatible,
    const std::vector<HydroStressFacePoint>* partitionedFaces,bool retainCanonicalRaw,
    const QuadtreeMesh* remapSourceMesh,const HydroStrainOperator* remapSource,
    const LAMP::Mesh::NativeLeafOverlap* remapOverlaps,bool rt0GreenOnly)
{
    metrics.layout().validateMesh(mesh);
    require(!points.empty(),"empty strain quadrature");
    auto data=std::make_shared<Data>(std::move(metrics));data->points=points;
    const auto& layout=data->metrics.layout();const bool rz=layout.domain().axisymmetric();
    data->components=rz ? 4 : 3;
    StressReuseEnvironment localReuseEntry;
    const bool localReuseEligible=stressCompatible&&(retainCanonicalRaw||rt0GreenOnly)&&!partitionedFaces&&localReuseEntry.read();
    data->localReuseBuildEntry=localReuseEntry;
    data->reuse.requested=remapSource!=nullptr;data->reuse.targetCells=static_cast<std::size_t>(layout.pressureCount());
    data->reuse.recomputedFitCells=data->reuse.targetCells;data->reuse.rebuiltFaces=static_cast<std::size_t>(layout.velocityCount());
    require(points.size()<=static_cast<std::size_t>(std::numeric_limits<Eigen::Index>::max()/data->components),"too many strain points");
    std::vector<Vector2D> centers(layout.pressureCount()),widths(layout.pressureCount());
    std::vector<std::vector<std::size_t>> pointRows(layout.pressureCount());
    for (int p=0;p<layout.pressureCount();++p) {
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[layout.nativeCell(p)]);
        centers[p]=cell.centerCoord();widths[p]=cell.size();
    }
    for (std::size_t q=0;q<points.size();++q) {
        const auto& sample=points[q];const int p=sample.cellDof;
        require(p>=0 && p<layout.pressureCount(),"invalid quadrature cell");
        require(std::isfinite(sample.point.x) && std::isfinite(sample.point.y)
            && std::isfinite(sample.physicalWeight) && sample.physicalWeight>=0,"invalid point/physical weight");
        const auto c=centers[p],w=widths[p];
        require(sample.point.x>=c.x-w.x/2 && sample.point.x<=c.x+w.x/2
            && sample.point.y>=c.y-w.y/2 && sample.point.y<=c.y+w.y/2,"quadrature point outside native cell");
        pointRows[p].push_back(q);
    }
    std::vector<NativeLeafFace> faces;faces.reserve(layout.velocityCount());
    std::vector<std::vector<int>> cellFaces(layout.pressureCount()),neighbors(layout.pressureCount());
    for (int u=0;u<layout.velocityCount();++u) {
        const auto f=layout.domain().face(mesh,layout.nativeFace(u));faces.push_back(f);
        const int a=f.negativeCell<0 ? -1 : layout.pressureDof(f.negativeCell);
        const int b=f.positiveCell<0 ? -1 : layout.pressureDof(f.positiveCell);
        if (a>=0) cellFaces[a].push_back(u);
        if (b>=0) cellFaces[b].push_back(u);
        if (a>=0 && b>=0) {neighbors[a].push_back(b);neighbors[b].push_back(a);}
    }
    auto recording=localReuseEligible?StressPreparationRecorder::make(centers,widths,faces,cellFaces,neighbors,pointRows,points.size())
        :std::unique_ptr<StressPreparationRecorder>{};
    StressLocalReuse reuse;
    if(remapSource&&remapSourceMesh&&remapOverlaps&&localReuseEligible)
        reuse=prepareStressLocalReuse(*remapSourceMesh,*remapSource->data_,mesh,data->metrics,*remapOverlaps,
            points,centers,widths,faces,cellFaces,neighbors,pointRows,data->reuse,recording.get(),rt0GreenOnly);
    if(!reuse.enabled){const bool requested=data->reuse.requested;data->reuse={};data->reuse.requested=requested;
        data->reuse.targetCells=static_cast<std::size_t>(layout.pressureCount());data->reuse.recomputedFitCells=data->reuse.targetCells;
        data->reuse.rebuiltFaces=static_cast<std::size_t>(layout.velocityCount());}
    StressCsrRebindBudget csrRebindBudget;
    const bool csrRebindEligible=reuse.enabled&&stressCsrRebindQualified();
    std::vector<std::array<Fit,2>> fitByCell;std::vector<Vector2D> fitCenters;
    if(stressCompatible){fitByCell.resize(layout.pressureCount());fitCenters.resize(layout.pressureCount());}
    std::vector<Eigen::Triplet<double,Eigen::Index>> volumeEntries;
    std::vector<Eigen::Triplet<double,Eigen::Index>> entries;
    const auto serialCells=[&] {
    FitTemplates templates;
    for (int p=0;p<layout.pressureCount();++p) if (!pointRows[p].empty()&&(!reuse.enabled||reuse.needFit[p])) {
        std::vector<int> patch{p},frontier{p};
        for (int ring=0;ring<2;++ring) {
            std::vector<int> next;
            for (int c : frontier) for (int n : neighbors[c])
                if (std::find(patch.begin(),patch.end(),n)==patch.end()) {patch.push_back(n);next.push_back(n);}
            frontier=std::move(next);
        }
        std::vector<int> support;
        for (int c : patch) support.insert(support.end(),cellFaces[c].begin(),cellFaces[c].end());
        std::sort(support.begin(),support.end());support.erase(std::unique(support.begin(),support.end()),support.end());
        if(recording)recording->save(p,patch,support);
        const bool reflected=rz && std::any_of(support.begin(),support.end(),[&](int u) {return faces[u].area==0;});
        const auto c=reflected ? Vector2D{0,centers[p].y} : centers[p];const auto w=widths[p];
        const std::array<Fit,2> fits{fit(data->metrics,faces,support,cellFaces[p],Axis::First,c,w,reflected,p,templates),
            fit(data->metrics,faces,support,cellFaces[p],Axis::Second,c,w,reflected,p,templates)};
        if(stressCompatible){fitByCell[p]=fits;fitCenters[p]=c;}
        if(reuse.enabled&&reuse.reuseRaw[p])continue; // fit was needed only by a rebuilt shared trace.
        for (std::size_t q : pointRows[p]) {
            const auto point=points[q].point;const double x=(point.x-c.x)/w.x,y=(point.y-c.y)/w.y;
            Eigen::Matrix<double,1,6> dx,dy;
            if(!rt0GreenOnly){
            dx<<0,1/w.x,0,2*x/w.x,y/w.x,0;
            dy<<0,0,1/w.y,0,x/w.y,2*y/w.y;
            }
            const auto value=basis(x,y);
            if (stressCompatible) for (int axis=0;axis<2;++axis) {
                const Eigen::RowVectorXd coefficients=value*fits[axis].coefficients;
                require(coefficients.allFinite(),"nonfinite stress volume value coefficients");
                for (Eigen::Index j=0;j<coefficients.size();++j) if (coefficients[j]!=0)
                    volumeEntries.emplace_back(static_cast<Eigen::Index>(q)*2+axis,
                        fits[axis].samples[static_cast<std::size_t>(j)].dof,
                        coefficients[j]*fits[axis].samples[static_cast<std::size_t>(j)].sign);
            }
            if(!rt0GreenOnly)for (int axis=0;axis<2;++axis) for (int component=0;component<data->components;++component) {
                Eigen::Matrix<double,1,6> functional=Eigen::Matrix<double,1,6>::Zero();
                if (component==0 && axis==0) functional=dx;
                if (component==1 && axis==1) functional=dy;
                if (rz && component==2 && axis==0) {
                    require(point.x>0 || reflected,"hoop axis limit requires retained axis");
                    functional=point.x==0 ? dx : value/point.x;
                }
                if (component==data->components-1) functional=(axis==0 ? dy : dx)/std::sqrt(2.0);
                const Eigen::RowVectorXd coefficients=functional*fits[axis].coefficients;
                require(coefficients.allFinite(),"nonfinite strain coefficients");
                for (Eigen::Index j=0;j<coefficients.size();++j) if (coefficients[j]!=0)
                    entries.emplace_back(static_cast<Eigen::Index>(q)*data->components+component,
                        fits[axis].samples[static_cast<std::size_t>(j)].dof,
                        coefficients[j]*fits[axis].samples[static_cast<std::size_t>(j)].sign);
            }
        }
    }
    };
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
    StrainFloatingState caller;
    if(strainParallelQualified(caller)) {
        try {
            using Triplet=Eigen::Triplet<double,Eigen::Index>;
            struct CellAssembly {std::vector<Triplet> entries,volumeEntries;};
            const int count=layout.pressureCount();
            std::vector<CellAssembly> assembled(static_cast<std::size_t>(count));
            std::vector<std::exception_ptr> failures(static_cast<std::size_t>(count));
            StrainFloatingState prepared;
            if(!prepared.capture())throw std::runtime_error("strain parallel environment capture failed");
            unsigned long raisedNative=0;unsigned raisedCsr=0;int setupFailure=0;
#pragma omp parallel
            {
                StrainFloatingState worker;
                const bool captured=worker.capture();
                bool ready=captured&&worker.enterFresh(prepared);
                std::unique_ptr<FitTemplates> workerTemplates;
                if(ready)try {workerTemplates=std::make_unique<FitTemplates>();}
                catch(...) {ready=false;}
                if(!ready) {
#pragma omp atomic write
                    setupFailure=1;
                }
#pragma omp for schedule(static) nowait
                for(int p=0;p<count;++p) {
                    if(!ready||pointRows[p].empty()||(reuse.enabled&&!reuse.needFit[p]))continue;
                    try {
                        auto& templates=*workerTemplates;
                        auto& entries=assembled[static_cast<std::size_t>(p)].entries;
                        auto& volumeEntries=assembled[static_cast<std::size_t>(p)].volumeEntries;
        std::vector<int> patch{p},frontier{p};
        for (int ring=0;ring<2;++ring) {
            std::vector<int> next;
            for (int c : frontier) for (int n : neighbors[c])
                if (std::find(patch.begin(),patch.end(),n)==patch.end()) {patch.push_back(n);next.push_back(n);}
            frontier=std::move(next);
        }
        std::vector<int> support;
        for (int c : patch) support.insert(support.end(),cellFaces[c].begin(),cellFaces[c].end());
        std::sort(support.begin(),support.end());support.erase(std::unique(support.begin(),support.end()),support.end());
        if(recording)recording->save(p,patch,support);
        const bool reflected=rz && std::any_of(support.begin(),support.end(),[&](int u) {return faces[u].area==0;});
        const auto c=reflected ? Vector2D{0,centers[p].y} : centers[p];const auto w=widths[p];
        const std::array<Fit,2> fits{fit(data->metrics,faces,support,cellFaces[p],Axis::First,c,w,reflected,p,templates),
            fit(data->metrics,faces,support,cellFaces[p],Axis::Second,c,w,reflected,p,templates)};
        if(stressCompatible){fitByCell[p]=fits;fitCenters[p]=c;}
        if(reuse.enabled&&reuse.reuseRaw[p])continue; // fit was needed only by a rebuilt shared trace.
        for (std::size_t q : pointRows[p]) {
            const auto point=points[q].point;const double x=(point.x-c.x)/w.x,y=(point.y-c.y)/w.y;
            Eigen::Matrix<double,1,6> dx,dy;
            if(!rt0GreenOnly){
            dx<<0,1/w.x,0,2*x/w.x,y/w.x,0;
            dy<<0,0,1/w.y,0,x/w.y,2*y/w.y;
            }
            const auto value=basis(x,y);
            if (stressCompatible) for (int axis=0;axis<2;++axis) {
                const Eigen::RowVectorXd coefficients=value*fits[axis].coefficients;
                require(coefficients.allFinite(),"nonfinite stress volume value coefficients");
                for (Eigen::Index j=0;j<coefficients.size();++j) if (coefficients[j]!=0)
                    volumeEntries.emplace_back(static_cast<Eigen::Index>(q)*2+axis,
                        fits[axis].samples[static_cast<std::size_t>(j)].dof,
                        coefficients[j]*fits[axis].samples[static_cast<std::size_t>(j)].sign);
            }
            if(!rt0GreenOnly)for (int axis=0;axis<2;++axis) for (int component=0;component<data->components;++component) {
                Eigen::Matrix<double,1,6> functional=Eigen::Matrix<double,1,6>::Zero();
                if (component==0 && axis==0) functional=dx;
                if (component==1 && axis==1) functional=dy;
                if (rz && component==2 && axis==0) {
                    require(point.x>0 || reflected,"hoop axis limit requires retained axis");
                    functional=point.x==0 ? dx : value/point.x;
                }
                if (component==data->components-1) functional=(axis==0 ? dy : dx)/std::sqrt(2.0);
                const Eigen::RowVectorXd coefficients=functional*fits[axis].coefficients;
                require(coefficients.allFinite(),"nonfinite strain coefficients");
                for (Eigen::Index j=0;j<coefficients.size();++j) if (coefficients[j]!=0)
                    entries.emplace_back(static_cast<Eigen::Index>(q)*data->components+component,
                        fits[axis].samples[static_cast<std::size_t>(j)].dof,
                        coefficients[j]*fits[axis].samples[static_cast<std::size_t>(j)].sign);
            }
        }
                    }catch(...) {failures[static_cast<std::size_t>(p)]=std::current_exception();}
                }
                if(ready) {
                    std::fenv_t after{};const bool got=std::fegetenv(&after)==0;const unsigned flags=_mm_getcsr()&0x3fu;
                    if(!got) {
#pragma omp atomic write
                        setupFailure=1;
                    }
#pragma omp critical(lamp_strain_floating_flags)
                    {if(got)raisedNative|=after._Fe_stat;raisedCsr|=flags;}
                }
                if(captured)worker.restore();
            }
            auto merged=prepared.native;merged._Fe_stat|=raisedNative;
            if(std::fesetenv(&merged)!=0)std::terminate();
            _mm_setcsr(prepared.csr|raisedCsr);
            if(setupFailure)throw std::runtime_error("strain parallel worker environment setup failed");
            for(const auto& failure:failures)if(failure)std::rethrow_exception(failure);
            std::size_t strainSize=0,volumeSize=0;
            for(const auto& cell:assembled) {
                if(cell.entries.size()>entries.max_size()-strainSize
                    ||cell.volumeEntries.size()>volumeEntries.max_size()-volumeSize)
                    throw std::length_error("strain parallel triplet allocation range");
                strainSize+=cell.entries.size();volumeSize+=cell.volumeEntries.size();
            }
            entries.reserve(strainSize);volumeEntries.reserve(volumeSize);
            // Restore original p/q/axis/component/j triplet order, including
            // duplicate coefficients; no arithmetic is performed by copying.
            for(auto& cell:assembled) {
                entries.insert(entries.end(),cell.entries.begin(),cell.entries.end());
                volumeEntries.insert(volumeEntries.end(),cell.volumeEntries.begin(),cell.volumeEntries.end());
                std::vector<Triplet>().swap(cell.entries);
                std::vector<Triplet>().swap(cell.volumeEntries);
            }
        }catch(...) {
            // Execution/resource fallback only: worker buffers have unwound.
            // Discard partial output, restore caller entry FP state, and
            // reproduce numerical guards in their original global order.
            entries.clear();volumeEntries.clear();
            if(stressCompatible)for(int p=0;p<layout.pressureCount();++p) {
                fitByCell[p]={};fitCenters[p]={};
            }
            caller.restore();
            serialCells();
        }
    }else serialCells();
#else
    serialCells();
#endif
    StressCsrCopyRows rawCopies,volumeCopies;
    if(csrRebindEligible){
        if(!rt0GreenOnly)rawCopies.prepare(points.size()*data->components,reuse,*reuse.raw,csrRebindBudget,true);
        volumeCopies.prepare(points.size()*2,reuse,*reuse.volume,csrRebindBudget,true);
    }
    if(reuse.enabled)for(int p=0;p<layout.pressureCount();++p)if(reuse.reuseRaw[p])for(std::size_t q:pointRows[p]){
        const auto oldQ=reuse.oldPoint[q];
        if(!rt0GreenOnly)for(int k=0;k<data->components;++k){
            const auto row=static_cast<Eigen::Index>(q)*data->components+k,oldRow=static_cast<Eigen::Index>(oldQ)*data->components+k;
            if(rawCopies.record(static_cast<std::size_t>(row),oldRow))reuse.copy(*reuse.raw,oldRow,[](int,double){});
            else reuse.copy(*reuse.raw,oldRow,[&](int u,double value){entries.emplace_back(row,u,value);});
        }
        for(int k=0;k<2;++k){
            const auto row=static_cast<Eigen::Index>(q)*2+k,oldRow=static_cast<Eigen::Index>(oldQ)*2+k;
            if(volumeCopies.record(static_cast<std::size_t>(row),oldRow))reuse.copy(*reuse.volume,oldRow,[](int,double){});
            else reuse.copy(*reuse.volume,oldRow,[&](int u,double value){volumeEntries.emplace_back(row,u,value);});
        }
    }
    const auto materializeCopies=[&](std::vector<Eigen::Triplet<double,Eigen::Index>>& triplets,
        const StressCsrCopyRows& copies,int components){
        if(!copies.count)return;
        // Same original p/q/k encounter order; no mapped-column sorting.
        for(int p=0;p<layout.pressureCount();++p)if(reuse.reuseRaw[p])for(std::size_t q:pointRows[p])for(int k=0;k<components;++k){
            const auto row=static_cast<Eigen::Index>(q)*components+k;
            if(copies.has(static_cast<std::size_t>(row)))
                copies.copy(static_cast<std::size_t>(row),[&](int u,double value){triplets.emplace_back(row,u,value);});
        }
    };
    if(!rt0GreenOnly){
    data->matrix.resize(static_cast<Eigen::Index>(points.size())*data->components,layout.velocityCount());
    if(!stressCopiedTripletMatrix(data->matrix,data->matrix.rows(),layout.velocityCount(),entries,
        rawCopies,csrRebindBudget,data->reuse.rawDirectReboundRows)){
        materializeCopies(entries,rawCopies,data->components);
        data->matrix.setFromTriplets(entries.begin(),entries.end());
    }
    data->matrix.prune(0.0);data->matrix.makeCompressed();
    for (Eigen::Index i=0;i<data->matrix.nonZeros();++i) require(std::isfinite(data->matrix.valuePtr()[i]),"nonfinite assembled strain");
    }
    rawCopies.reset(); // Consumed routing is not retained into later face/integral work.
    std::shared_ptr<Data> canonicalRawOwner;
    if(retainCanonicalRaw){
        // Only the builtin factory sets this private flag after generating the
        // canonical p/x/y sequence. Preserve the ACTUAL matrix before Ehat's
        // original mean-replacement paths mutate it. No fits survive the build.
        auto raw=std::make_shared<Data>(data->metrics);
        raw->points=data->points;raw->components=data->components;raw->matrix=data->matrix;
        data->canonicalRaw=raw;canonicalRawOwner=std::move(raw);
    }
    if (stressCompatible) {
        using Row=StressGeometryRow;
        const int count=layout.pressureCount(),components=data->components;
        const long double root2=rt0GreenOnly?0.L:std::sqrt(2.0L);
        const auto addRow=[](Row& target,const Row& source,long double factor) {
            for(const auto& e:source)target[e.first]+=factor*e.second;
        };
        const auto matrixRow=[](const Matrix& matrix,Eigen::Index row) {
            Row result;for(Matrix::InnerIterator e(matrix,row);e;++e)result[static_cast<int>(e.col())]+=e.value();return result;
        };
        const auto literalSetMatrix=[](Matrix& matrix,Eigen::Index rows,int cols,const std::vector<Row>& input) {
            std::vector<Eigen::Triplet<double,Eigen::Index>> triplets;
            for(std::size_t row=0;row<input.size();++row)for(const auto& e:input[row])if(e.second!=0) {
                const double value=static_cast<double>(e.second);require(std::isfinite(value),"nonfinite stress geometry coefficient");
                if(value!=0)triplets.emplace_back(static_cast<Eigen::Index>(row),e.first,value);
            }
            matrix.resize(rows,cols);matrix.setFromTriplets(triplets.begin(),triplets.end());matrix.prune(0.0);matrix.makeCompressed();
        };
        const auto setMatrix=[&](Matrix& matrix,Eigen::Index rows,int cols,std::vector<Row>& input,
            const StressCsrCopyRows& copies) {
            if(!copies.count){literalSetMatrix(matrix,rows,cols,input);return;}
            std::size_t upper=0;bool representable=true;
            if(copies.count)for(std::size_t row=0;row<input.size();++row){
                const auto size=copies.has(row)
                    ?static_cast<std::size_t>(copies.source->outerIndexPtr()[copies.sourceRows[row]+1]
                        -copies.source->outerIndexPtr()[copies.sourceRows[row]])
                    :static_cast<std::size_t>(input[row].end()-input[row].begin());
                if(size>std::numeric_limits<std::size_t>::max()-upper){representable=false;break;}
                upper+=size;
            }
            StressCsrPreparedMatrix prepared;
            if(copies.count&&representable&&prepared.prepare(rows,cols,upper,csrRebindBudget)){
                Matrix& output=*prepared.matrix;
                for(std::size_t row=0;row<input.size();++row){
                    output.startVec(static_cast<Eigen::Index>(row));
                    if(copies.has(row)){
                        copies.copy(row,[&](int u,double value){output.insertBack(static_cast<Eigen::Index>(row),u)=value;});
                    }else for(const auto& e:input[row])if(e.second!=0){
                        const double value=static_cast<double>(e.second);require(std::isfinite(value),"nonfinite stress geometry coefficient");
                        if(value!=0)output.insertBack(static_cast<Eigen::Index>(row),e.first)=value;
                    }
                }
                output.finalize();output.prune(0.0);output.makeCompressed();
                matrix.swap(output);data->reuse.directReboundRows+=copies.count;return;
            }
            // Resource/unsupported path materializes the same original ordered
            // Row values. Dirty Row arithmetic above is never recomputed.
            for(std::size_t row=0;row<input.size();++row)if(copies.has(row)){
                input[row]=copies.reuse->row(*copies.source,copies.sourceRows[row]);
                ++data->reuse.literalReboundRows;
            }
            literalSetMatrix(matrix,rows,cols,input);
        };
        StressCsrCopyRows faceCopies,integralCopies;
        // A source view is installed only after the original reuse certificate.
        if(csrRebindEligible){
            faceCopies.prepare(faces.size()*4,reuse,*reuse.faceVelocity,csrRebindBudget,true);
            if(!rt0GreenOnly)integralCopies.prepare(static_cast<std::size_t>(count)*components,reuse,*reuse.integrals,csrRebindBudget,true);
        }
        data->volumeVelocity.resize(static_cast<Eigen::Index>(points.size())*2,layout.velocityCount());
        if(!stressCopiedTripletMatrix(data->volumeVelocity,data->volumeVelocity.rows(),layout.velocityCount(),volumeEntries,
            volumeCopies,csrRebindBudget,data->reuse.volumeDirectReboundRows)){
            materializeCopies(volumeEntries,volumeCopies,2);
            data->volumeVelocity.setFromTriplets(volumeEntries.begin(),volumeEntries.end());
        }
        data->volumeVelocity.prune(0.0);data->volumeVelocity.makeCompressed();
        volumeCopies.reset();
        // ONE physical face partition/trace, shared by both incident cells.
        std::vector<Row> faceRows;faceRows.reserve(faces.size()*4);
        const double gauss=1/std::sqrt(3.0);
        const auto addValue=[&](Row& row,int p,int axis,Vector2D point,long double factor) {
            const auto center=fitCenters[p],width=widths[p];const auto& f=fitByCell[p][axis];
            const auto value=basis((point.x-center.x)/width.x,(point.y-center.y)/width.y);
            const Eigen::RowVectorXd coefficients=value*f.coefficients;
            require(coefficients.allFinite(),"nonfinite shared stress face value coefficients");
            for(Eigen::Index j=0;j<coefficients.size();++j)if(coefficients[j]!=0)
                row[f.samples[static_cast<std::size_t>(j)].dof]+=factor*coefficients[j]*f.samples[static_cast<std::size_t>(j)].sign;
        };
        if(partitionedFaces){
            std::vector<long double> faceWeights(static_cast<std::size_t>(layout.velocityCount()));
            for(const auto& sample:*partitionedFaces){
                const int u=sample.velocityDof;
                require(u>=0&&u<layout.velocityCount(),"partitioned stress face velocity DOF out of range");
                const auto& f=faces[u];const int a=f.negativeCell<0?-1:layout.pressureDof(f.negativeCell),b=f.positiveCell<0?-1:layout.pressureDof(f.positiveCell);
                require(f.tangentSpacing>0&&std::isfinite(f.tangentSpacing),"invalid actual stress face length");
                const int normal=f.horizontal?1:0,tangent=1-normal,sides=(a>=0?1:0)+(b>=0?1:0);
                require(sides>0,"stress face has no actual adjacent cell");
                require(sample.nativeFace==layout.nativeFace(u)&&sample.negativeCellDof==a&&sample.positiveCellDof==b,
                    "partitioned stress face has foreign native incidence");
                const auto point=sample.point;const auto weight=sample.physicalWeight;
                require(std::isfinite(point.x)&&std::isfinite(point.y)&&std::isfinite(weight)&&weight>=0,
                    "invalid actual stress face quadrature point/weight");
                require((f.horizontal?point.y==f.center.y:point.x==f.center.x)
                    &&std::abs(f.horizontal?point.x-f.center.x:point.y-f.center.y)<=f.tangentSpacing/2,
                    "partitioned stress face point outside its actual native face");
                require(!rz||point.x>=0,"stress face quadrature outside physical RZ radius");
                require(f.area==0?weight==0:weight>0,"partitioned face has a nonphysical zero/axis weight");
                faceWeights[static_cast<std::size_t>(u)]+=weight;
                const std::size_t q=data->facePoints.size();data->facePoints.push_back(sample);faceRows.resize((q+1)*2);
                faceRows[2*q+normal][u]=1;
                if(a>=0)addValue(faceRows[2*q+tangent],a,tangent,point,1.L/sides);
                if(b>=0)addValue(faceRows[2*q+tangent],b,tangent,point,1.L/sides);
            }
            for(int u=0;u<layout.velocityCount();++u){
                const auto actual=faces[u].area;const auto sum=faceWeights[static_cast<std::size_t>(u)];
                require(actual==0?sum==0:std::abs(sum-actual)<=1e-12L*actual,
                    "stress face quadrature does not have the actual full native measure");
                if(actual>0)data->compatibility.maximumRelativeFaceQuadratureWeightDefect=std::max(
                    data->compatibility.maximumRelativeFaceQuadratureWeightDefect,static_cast<double>(std::abs(sum-actual)/actual));
            }
        }
        else for(int u=0;u<layout.velocityCount();++u) {
            const auto& f=faces[u];const int a=f.negativeCell<0?-1:layout.pressureDof(f.negativeCell),b=f.positiveCell<0?-1:layout.pressureDof(f.positiveCell);
            require(f.tangentSpacing>0&&std::isfinite(f.tangentSpacing),"invalid actual stress face length");
            const int normal=f.horizontal?1:0,tangent=1-normal,sides=(a>=0?1:0)+(b>=0?1:0);
            require(sides>0,"stress face has no actual adjacent cell");
            for(double g:{-gauss,gauss}) {
                auto point=f.center;if(f.horizontal)point.x+=g*f.tangentSpacing/2;else point.y+=g*f.tangentSpacing/2;
                require(!rz||point.x>=0,"stress face quadrature outside physical RZ radius");
                const double weight=f.area/2*(rz&&f.horizontal?point.x/f.center.x:1);
                require(std::isfinite(weight)&&weight>=0,"invalid actual stress face quadrature weight");
                const std::size_t q=data->facePoints.size();data->facePoints.push_back({u,layout.nativeFace(u),a,b,point,weight});faceRows.resize((q+1)*2);
                if(reuse.enabled&&reuse.reuseFace[u]){
                    const auto oldQ=reuse.oldFacePoints[reuse.oldFace[u]][g<0?0:1];
                    const auto& actualSource=(*reuse.sourceFacePoints)[oldQ];
                    if(!stressSamePoint(point,actualSource.point)||!stressSameDouble(weight,actualSource.physicalWeight))
                        throw StressLocalReuseMismatch{};
                    for(int k=0;k<2;++k){
                        const auto oldRow=static_cast<Eigen::Index>(oldQ)*2+k;
                        if(faceCopies.record(2*q+k,oldRow))
                            reuse.copy(*reuse.faceVelocity,oldRow,[](int,double){});
                        else{faceRows[2*q+k]=reuse.row(*reuse.faceVelocity,oldRow);++data->reuse.literalReboundRows;}
                    }
                }else{
                faceRows[2*q+normal][u]=1; // Native normal mean, not a fitted normal value.
                if(a>=0)addValue(faceRows[2*q+tangent],a,tangent,point,1.L/sides);
                if(b>=0)addValue(faceRows[2*q+tangent],b,tangent,point,1.L/sides);
                }
            }
        }
        setMatrix(data->faceVelocity,static_cast<Eigen::Index>(faceRows.size()),layout.velocityCount(),faceRows,faceCopies);
        if(rt0GreenOnly){
            // Preserve builtin point/RZ and positive finite cell-measure gates;
            // no unrequested component moment or compatibility report is built.
            std::vector<long double> weights(count);
            for(const auto& point:points){
                weights[point.cellDof]+=point.physicalWeight;
                if(rz)require(point.point.x>0,"stress volume Gauss point must have positive physical radius");
            }
            for(int p=0;p<count;++p)
                require(std::isfinite(weights[p])&&weights[p]>0,"nonpositive actual stress quadrature weight sum");
        }else{
        std::vector<Row> target(static_cast<std::size_t>(count)*components),rawSums(target.size()),dx(count),dy(count),hr(count),hz(count);
        std::vector<long double> weights(count);
        const auto originalMomentAssembly=[&] {
        for(int u=0;u<layout.velocityCount();++u) {
            const auto& f=faces[u];for(int side=0;side<2;++side) {
                const int cell=side==0?f.negativeCell:f.positiveCell;if(cell<0)continue;
                const int p=layout.pressureDof(cell);if(reuse.enabled&&reuse.reuseStress[p])continue;
                (f.horizontal?dy[p]:dx[p])[u]+=(side==0?1.L:-1.L)*f.area;
            }
        }
        for(std::size_t q=0;q<points.size();++q) {
            const int p=points[q].cellDof;const long double w=points[q].physicalWeight;weights[p]+=w;
            if(reuse.enabled&&reuse.reuseStress[p]){if(rz)require(points[q].point.x>0,"stress volume Gauss point must have positive physical radius");continue;}
            for(int k=0;k<components;++k)addRow(rawSums[static_cast<std::size_t>(p)*components+k],matrixRow(data->matrix,static_cast<Eigen::Index>(q)*components+k),w);
            if(rz) {
                require(points[q].point.x>0,"stress volume Gauss point must have positive physical radius");
                addRow(hr[p],matrixRow(data->volumeVelocity,static_cast<Eigen::Index>(q)*2),w/points[q].point.x);
                addRow(hz[p],matrixRow(data->volumeVelocity,static_cast<Eigen::Index>(q)*2+1),w/points[q].point.x);
            }
        }
        for(std::size_t q=0;q<data->facePoints.size();++q) {
            const auto& f=data->facePoints[q];const int tangent=layout.normalAxis(f.velocityDof)==Axis::First?1:0;
            if(reuse.enabled&&(f.negativeCellDof<0||reuse.reuseStress[f.negativeCellDof])
                &&(f.positiveCellDof<0||reuse.reuseStress[f.positiveCellDof]))continue;
            // div(2*pi*r*v) needs the additional -integral(v_z/r) below.
            const auto trace=matrixRow(data->faceVelocity,static_cast<Eigen::Index>(q)*2+tangent);
            for(int side=0;side<2;++side) {const int p=side==0?f.negativeCellDof:f.positiveCellDof;if(p<0||(reuse.enabled&&reuse.reuseStress[p]))continue;
                addRow(target[static_cast<std::size_t>(p)*components+components-1],trace,(side==0?1.L:-1.L)*f.physicalWeight/root2);}
        }
        };
        bool parallelMomentsComplete=false;
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
        StrainFloatingState momentCaller;
        if(strainParallelQualified(momentCaller)) {
            try {
                // Preserve every incident contribution in the old u/side and
                // q/side order. No worker searches the whole native domain.
                std::vector<std::vector<std::pair<int,int>>> normalIncidence(count);
                std::vector<std::vector<std::pair<std::size_t,int>>> traceIncidence(count);
                for(int u=0;u<layout.velocityCount();++u) {
                    const auto& f=faces[u];for(int side=0;side<2;++side) {
                        const int cell=side==0?f.negativeCell:f.positiveCell;if(cell<0)continue;
                        const int p=layout.pressureDof(cell);normalIncidence[p].emplace_back(u,side);
                    }
                }
                for(std::size_t q=0;q<data->facePoints.size();++q) {
                    const auto& f=data->facePoints[q];
                    for(int side=0;side<2;++side) {
                        const int p=side==0?f.negativeCellDof:f.positiveCellDof;if(p<0)continue;
                        traceIncidence[p].emplace_back(q,side);
                    }
                }
                std::vector<std::exception_ptr> momentFailures(count);
                StrainFloatingState momentPrepared;
                if(!momentPrepared.capture())throw std::runtime_error("stress moment parallel environment capture failed");
                unsigned long raisedNative=0;unsigned raisedCsr=0;int setupFailure=0;
#pragma omp parallel
                {
                    StrainFloatingState worker;
                    const bool captured=worker.capture();
                    const bool ready=captured&&worker.enterFresh(momentPrepared);
                    if(!ready){
#pragma omp atomic write
                        setupFailure=1;
                    }
#pragma omp for schedule(static) nowait
                    for(int p=0;p<count;++p) {
                        if(!ready)continue;
                        try {
                            if(reuse.enabled&&reuse.reuseStress[p]){
                                for(std::size_t q:pointRows[p]){weights[p]+=points[q].physicalWeight;
                                    if(rz)require(points[q].point.x>0,"stress volume Gauss point must have positive physical radius");}
                                continue;
                            }
                            for(const auto& incident:normalIncidence[p]) {
                                const int u=incident.first,side=incident.second;const auto& f=faces[u];
                                (f.horizontal?dy[p]:dx[p])[u]+=(side==0?1.L:-1.L)*f.area;
                            }
                            for(std::size_t q:pointRows[p]) {
                                const long double w=points[q].physicalWeight;weights[p]+=w;
                                for(int k=0;k<components;++k)addRow(rawSums[static_cast<std::size_t>(p)*components+k],matrixRow(data->matrix,static_cast<Eigen::Index>(q)*components+k),w);
                                if(rz) {
                                    require(points[q].point.x>0,"stress volume Gauss point must have positive physical radius");
                                    addRow(hr[p],matrixRow(data->volumeVelocity,static_cast<Eigen::Index>(q)*2),w/points[q].point.x);
                                    addRow(hz[p],matrixRow(data->volumeVelocity,static_cast<Eigen::Index>(q)*2+1),w/points[q].point.x);
                                }
                            }
                            for(const auto& incident:traceIncidence[p]) {
                                const auto q=incident.first;const int side=incident.second;const auto& f=data->facePoints[q];
                                const int tangent=layout.normalAxis(f.velocityDof)==Axis::First?1:0;
                                const auto trace=matrixRow(data->faceVelocity,static_cast<Eigen::Index>(q)*2+tangent);
                                addRow(target[static_cast<std::size_t>(p)*components+components-1],trace,(side==0?1.L:-1.L)*f.physicalWeight/root2);
                            }
                        }catch(...){momentFailures[p]=std::current_exception();}
                    }
                    if(ready) {
                        std::fenv_t after{};const bool got=std::fegetenv(&after)==0;const unsigned flags=_mm_getcsr()&0x3fu;
                        if(!got){
#pragma omp atomic write
                            setupFailure=1;
                        }
#pragma omp critical(lamp_stress_moment_floating_flags)
                        {if(got)raisedNative|=after._Fe_stat;raisedCsr|=flags;}
                    }
                    if(captured)worker.restore();
                }
                auto merged=momentPrepared.native;merged._Fe_stat|=raisedNative;
                if(std::fesetenv(&merged)!=0)std::terminate();
                _mm_setcsr(momentPrepared.csr|raisedCsr);
                if(setupFailure)throw std::runtime_error("stress moment parallel worker environment setup failed");
                for(const auto& failure:momentFailures)if(failure)std::rethrow_exception(failure);
                parallelMomentsComplete=true;
            }catch(...) {
                // Discard partial independent rows before original global-order
                // replay. No physical inventory/matrix owner has been changed.
                momentCaller.restore();
                for(auto& row:target)row=Row{};for(auto& row:rawSums)row=Row{};
                for(auto& row:dx)row=Row{};for(auto& row:dy)row=Row{};
                for(auto& row:hr)row=Row{};for(auto& row:hz)row=Row{};
                std::fill(weights.begin(),weights.end(),0.L);
                originalMomentAssembly();parallelMomentsComplete=true;
            }
        }
#endif
        if(!parallelMomentsComplete)originalMomentAssembly();
        data->compatibility.cellQuadratureWeightSums.resize(count);data->compatibility.cellQuadratureWeightDefects.resize(count);
        for(int p=0;p<count;++p) {
            require(std::isfinite(weights[p])&&weights[p]>0,"nonpositive actual stress quadrature weight sum");
            const auto base=static_cast<std::size_t>(p)*components;
            if(reuse.enabled&&reuse.reuseStress[p])for(int k=0;k<components;++k){
                const auto oldRow=static_cast<Eigen::Index>(reuse.oldCell[p])*components+k;
                if(integralCopies.record(base+k,oldRow))
                    reuse.copy(*reuse.integrals,oldRow,[](int,double){});
                else{target[base+k]=reuse.row(*reuse.integrals,oldRow);++data->reuse.literalReboundRows;}
            }
            else{target[base]=dx[p];target[base+1]=dy[p];
                if(rz) {addRow(target[base],hr[p],-1);target[base+2]=hr[p];addRow(target[base+3],hz[p],-1/root2);}}
            const long double defect=weights[p]-data->metrics.cellVolumes()[p];
            data->compatibility.cellQuadratureWeightSums[p]=static_cast<double>(weights[p]);
            data->compatibility.cellQuadratureWeightDefects[p]=static_cast<double>(defect);
            data->compatibility.maximumRelativeQuadratureWeightDefect=std::max(data->compatibility.maximumRelativeQuadratureWeightDefect,
                static_cast<double>(std::abs(defect)/data->metrics.cellVolumes()[p]));
        }
        setMatrix(data->componentIntegrals,static_cast<Eigen::Index>(target.size()),layout.velocityCount(),target,integralCopies);
        // Pure geometry mean replacement; no T, physical force or residual is
        // consulted. The original raw zero-mean strain stays at every point.
        // Keep the literal triplet implementation for unsupported execution
        // contracts and replay it after any optional direct-assembly failure.
        const auto originalCompatibleAssembly=[&] {
        std::vector<Eigen::Triplet<double,Eigen::Index>> compatibleEntries;
        for(Eigen::Index row=0;row<data->matrix.rows();++row){
            const auto q=static_cast<std::size_t>(row/components),p=static_cast<std::size_t>(points[q].cellDof);
            if(reuse.enabled&&reuse.reuseStress[p]){
                reuse.copy(*reuse.stress,static_cast<Eigen::Index>(reuse.oldPoint[q])*components+row%components,
                    [&](int u,double value){compatibleEntries.emplace_back(row,u,value);});
                ++data->reuse.literalReboundRows;
            }else for(Matrix::InnerIterator e(data->matrix,row);e;++e)compatibleEntries.emplace_back(row,e.col(),e.value());
        }
        // This build's current-cell mean correction is identical at each of
        // its volume points. Keep first-use arithmetic/column order and retain
        // rounded zeros for the original per-point finite/zero checks below.
        // Unknown or unmasked native environments use the literal old loop.
        bool reuseMeanRows=viscosityProductPlanEnvironment();
        int meanRowCell=-1;
        std::array<bool,4> meanRowReady{};
        std::array<std::vector<std::pair<int,double>>,4> meanRows;
        for(std::size_t q=0;q<points.size();++q)for(int k=0;k<components;++k) {
            const int p=points[q].cellDof;const auto component=static_cast<std::size_t>(p)*components+k;
            if(reuse.enabled&&reuse.reuseStress[p])continue;
            if(!reuseMeanRows) {
                Row delta=target[component];addRow(delta,rawSums[component],-1);
                for(const auto& e:delta)if(e.second!=0) {
                    const double value=static_cast<double>(e.second/weights[p]);require(std::isfinite(value),"nonfinite stress-compatible mean coefficient");
                    if(value!=0)compatibleEntries.emplace_back(static_cast<Eigen::Index>(q)*components+k,e.first,value);
                }
            }else {
                if(p!=meanRowCell) {
                    meanRowCell=p;meanRowReady.fill(false);
                    for(auto& row:meanRows)row.clear();
                }
                auto& row=meanRows[static_cast<std::size_t>(k)];
                if(!meanRowReady[static_cast<std::size_t>(k)]) {
                    Row delta=target[component];addRow(delta,rawSums[component],-1);
                    for(const auto& e:delta)if(e.second!=0) {
                        const double value=static_cast<double>(e.second/weights[p]);require(std::isfinite(value),"nonfinite stress-compatible mean coefficient");
                        if(value!=0)compatibleEntries.emplace_back(static_cast<Eigen::Index>(q)*components+k,e.first,value);
                        // Only the optional record allocation may fall back;
                        // never catch the original arithmetic/guard/triplet.
                        if(reuseMeanRows)try {row.emplace_back(e.first,value);}
                        catch(const std::bad_alloc&) {
                            reuseMeanRows=false;for(auto& cached:meanRows)cached.clear();
                        }catch(const std::length_error&) {
                            reuseMeanRows=false;for(auto& cached:meanRows)cached.clear();
                        }
                    }
                    meanRowReady[static_cast<std::size_t>(k)]=reuseMeanRows;
                }else for(const auto& e:row) {
                    const double value=e.second;require(std::isfinite(value),"nonfinite stress-compatible mean coefficient");
                    if(value!=0)compatibleEntries.emplace_back(static_cast<Eigen::Index>(q)*components+k,e.first,value);
                }
            }
        }
        data->matrix.setZero();data->matrix.setFromTriplets(compatibleEntries.begin(),compatibleEntries.end());data->matrix.prune(0.0);data->matrix.makeCompressed();
        for(Eigen::Index i=0;i<data->matrix.nonZeros();++i)require(std::isfinite(data->matrix.valuePtr()[i]),"nonfinite compatible stress strain");
        };
        bool directCompatibleComplete=false;
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE) && EIGEN_WORLD_VERSION==3 && EIGEN_MAJOR_VERSION==4 && EIGEN_MINOR_VERSION==1
        StrainFloatingState compatibleEntry;
        if(viscosityProductPlanEnvironment()&&compatibleEntry.capture()) {
            const auto reboundBeforeCompatible=data->reuse.directReboundRows;
            const auto literalBeforeCompatible=data->reuse.literalReboundRows;
            try {
                const Matrix& raw=data->matrix;
                Matrix corrected(raw.rows(),raw.cols());corrected.reserve(raw.nonZeros());
                int meanCell=-1;
                std::array<bool,4> ready{};
                std::array<std::vector<std::pair<int,double>>,4> corrections;
                for(std::size_t q=0;q<points.size();++q)for(int k=0;k<components;++k) {
                    const int p=points[q].cellDof;
                    const auto component=static_cast<std::size_t>(p)*components+k;
                    if(reuse.enabled&&reuse.reuseStress[p]){
                        const auto output=static_cast<Eigen::Index>(q)*components+k;
                        const auto oldRow=static_cast<Eigen::Index>(reuse.oldPoint[q])*components+k;
                        if(csrRebindEligible&&stressCsrOrderedRow(reuse,*reuse.stress,oldRow)){
                            corrected.startVec(output);
                            reuse.copy(*reuse.stress,oldRow,[&](int u,double value){corrected.insertBack(output,u)=value;});
                            ++data->reuse.directReboundRows;
                        }else{
                            const auto copied=reuse.row(*reuse.stress,oldRow);
                            corrected.startVec(output);for(const auto& entry:copied)corrected.insertBack(output,entry.first)=static_cast<double>(entry.second);
                            ++data->reuse.literalReboundRows;
                        }
                        continue;
                    }
                    if(p!=meanCell) {
                        meanCell=p;ready.fill(false);
                        for(auto& row:corrections)row.clear();
                    }
                    auto& row=corrections[static_cast<std::size_t>(k)];
                    if(!ready[static_cast<std::size_t>(k)]) {
                        Row delta=target[component];addRow(delta,rawSums[component],-1);
                        for(const auto& e:delta)if(e.second!=0) {
                            const double value=static_cast<double>(e.second/weights[p]);require(std::isfinite(value),"nonfinite stress-compatible mean coefficient");
                            row.emplace_back(e.first,value); // Retain RN zeros for the original guard/filter.
                        }
                        ready[static_cast<std::size_t>(k)]=true;
                    }else for(const auto& e:row)
                        require(std::isfinite(e.second),"nonfinite stress-compatible mean coefficient");
                    const auto output=static_cast<Eigen::Index>(q)*components+k;
                    corrected.startVec(output);
                    Matrix::InnerIterator entry(raw,output);
                    std::size_t at=0;
                    // Original correction triplets omit rounded zero; do not
                    // add it to a raw value (including a signed raw zero).
                    const auto skipZeros=[&] {while(at<row.size()&&row[at].second==0)++at;};
                    skipZeros();
                    while(entry||at<row.size()) {
                        if(at==row.size()||(entry&&entry.col()<row[at].first)) {
                            corrected.insertBack(output,entry.col())=entry.value();++entry;
                        }else if(!entry||row[at].first<entry.col()) {
                            corrected.insertBack(output,row[at].first)=row[at].second;++at;skipZeros();
                        }else {
                            // Eigen 3.4.1 collapseDuplicates receives the raw
                            // triplet first, then exactly one correction.
                            corrected.insertBack(output,entry.col())=entry.value()+row[at].second;
                            ++entry;++at;skipZeros();
                        }
                    }
                }
                corrected.finalize();corrected.prune(0.0);corrected.makeCompressed();
                for(Eigen::Index i=0;i<corrected.nonZeros();++i)require(std::isfinite(corrected.valuePtr()[i]),"nonfinite compatible stress strain");
                data->matrix.swap(corrected);directCompatibleComplete=true;
            }catch(...) {
                // The raw matrix is untouched until the final nonallocating
                // swap. Restore controls/status before original-order replay.
                compatibleEntry.restore();
                data->reuse.directReboundRows=reboundBeforeCompatible;
                data->reuse.literalReboundRows=literalBeforeCompatible;
                originalCompatibleAssembly();directCompatibleComplete=true;
            }
        }
#endif
        if(!directCompatibleComplete)originalCompatibleAssembly();
        }
        data->stressCompatible=!rt0GreenOnly;data->rt0GreenOnly=rt0GreenOnly;
    }
    data->reuse.csrRebindScratchBytes=csrRebindBudget.peak;
    StressReuseEnvironment localReuseExit;
    data->localReuseEligible=localReuseEligible&&localReuseExit.read()&&localReuseEntry.sameControl(localReuseExit);
    data->localReuseEnvironment=localReuseExit;
    if(recording){
        // Only optional storage is restored; original numerical work and its
        // complete producer exit receipt have already finished.
        StressPreparationScope optional;
        if(optional.captured&&data->localReuseEligible){
            auto prepared=recording->finish(centers,widths,faces,cellFaces,neighbors,pointRows);
            if(prepared){data->preparation=prepared;if(canonicalRawOwner)canonicalRawOwner->preparation=std::move(prepared);}
        }
        recording.reset(); // Dispose failed copies before restoring FP/errno.
    }
    return HydroStrainOperator(std::move(data));
}
std::vector<double> HydroStrainOperator::apply(const QuadtreeMesh& mesh,const HydroNativeValues& values) const
{
    require(!data_->rt0GreenOnly,"RT0 Green geometry has no requested strain application");
    validate(mesh,values,metrics());return multiply(matrix(),values.normalVelocity());
}
namespace
{
// Conservative certificate of the already assembled K, never a replacement
// product. Known conservative Eigen/plan products use one weighted multiply,
// one product, and at most m-1 additions for m common stored strain rows.
bool gramProductImplementationKnown() noexcept
{
#if defined(_MSC_VER) && defined(_M_X64) && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST) && EIGEN_WORLD_VERSION==3 && ((EIGEN_MAJOR_VERSION==3 && EIGEN_MINOR_VERSION==9) || (EIGEN_MAJOR_VERSION==4 && EIGEN_MINOR_VERSION==1))
    return std::numeric_limits<double>::is_iec559
        &&std::numeric_limits<double>::digits==53
        &&std::numeric_limits<double>::has_denorm==std::denorm_present;
#else
    return false;
#endif
}
struct GramBoundUnavailable {};
struct GramBoundScope
{
    std::fenv_t native{};unsigned csr=0;int savedErrno=0;bool captured=false;
    GramBoundScope() noexcept:savedErrno(errno)
    {
#if defined(_MSC_VER) && defined(_M_X64)
        csr=_mm_getcsr();
#endif
        captured=std::fegetenv(&native)==0;errno=savedErrno;
    }
    ~GramBoundScope() noexcept
    {
        if(captured&&std::fesetenv(&native)!=0)std::terminate();
#if defined(_MSC_VER) && defined(_M_X64)
        _mm_setcsr(csr);
#endif
        errno=savedErrno;
    }
};
// Each positive arithmetic result is enclosed independently. In particular,
// a positive exact product rounded to zero advances to denorm_min, not zero.
double gramUpperSuccessor(double value)
{
    if(!std::isfinite(value)||value<0)throw GramBoundUnavailable{};
    const double result=std::nextafter(value,std::numeric_limits<double>::infinity());
    if(!std::isfinite(result))throw GramBoundUnavailable{};
    return result;
}
double gramUpperAdd(double a,double b)
{
    if(a==0)return b;if(b==0)return a;
    return gramUpperSuccessor(a+b);
}
double gramUpperMultiply(double a,double b)
{
    if(a==0||b==0)return 0;
    return gramUpperSuccessor(a*b);
}

}
struct HydroEnergyViscosity::Data
{
    explicit Data(HydroStrainOperator input) : strain(std::move(input)) {}
    HydroStrainOperator strain;
    std::vector<double> viscosity,weights;
    HydroViscosityProductReport product;
    std::shared_ptr<const Matrix> matrix;
    StressReuseEnvironment gramBuildEntry,gramBuildExit;
    bool gramSourceQualified=false;
    // Fresh actual rowWeights only. Integer magnitude comparison introduces no
    // additional floating operation/denormal exception into the original build.
    double gramMaximumWeight=0;
    bool gramWeightsValid=true;
    void rememberGramWeight(double value) noexcept
    {
        if constexpr(sizeof(double)==sizeof(std::uint64_t)) {
            std::uint64_t bits=0,maximum=0;
            std::memcpy(&bits,&value,sizeof bits);
            std::memcpy(&maximum,&gramMaximumWeight,sizeof maximum);
            const auto magnitude=bits&UINT64_C(0x7fffffffffffffff);
            if(magnitude>=UINT64_C(0x7ff0000000000000)
                ||((bits&UINT64_C(0x8000000000000000))!=0&&magnitude!=0))gramWeightsValid=false;
            if(magnitude>maximum)std::memcpy(&gramMaximumWeight,&magnitude,sizeof magnitude);
        }else gramWeightsValid=false;
    }
};
HydroEnergyViscosity::HydroEnergyViscosity(std::shared_ptr<const Data> data) : data_(std::move(data)) {}
const HydroStrainOperator& HydroEnergyViscosity::strain() const noexcept {return data_->strain;}
const std::vector<double>& HydroEnergyViscosity::viscosity() const noexcept {return data_->viscosity;}
const std::vector<double>& HydroEnergyViscosity::rowWeights() const noexcept {return data_->weights;}
const Matrix& HydroEnergyViscosity::matrix() const noexcept {return *data_->matrix;}
const HydroViscosityProductReport& HydroEnergyViscosity::productReport() const noexcept {return data_->product;}
std::optional<double> HydroEnergyViscosity::coarseScaledGramLowerErrorBound(double maximumScale) const
{
    const GramBoundScope saved;
    if(!saved.captured||!data_->gramSourceQualified)return std::nullopt;
    StressReuseEnvironment caller;
    if(!caller.read()||!caller.sameControl(data_->gramBuildExit))return std::nullopt;
    try {
        const auto& e=data_->strain.matrix();
        if(!std::isfinite(maximumScale)||maximumScale<=0||!e.isCompressed()||e.rows()<0||e.cols()<0
            ||static_cast<std::size_t>(e.rows())!=data_->weights.size()
            ||data_->matrix->rows()!=e.cols()||data_->matrix->cols()!=e.cols())return std::nullopt;
        constexpr std::size_t exactLimit=GramGeometryStatistics::exactLimit;
        if(!data_->gramWeightsValid)return std::nullopt;
        GramGeometryStatistics geometry;
        {
            const std::lock_guard<std::mutex> guard(data_->strain.data_->viscosityProductMutex);
            auto& savedGeometry=data_->strain.data_->gramGeometryStatistics;
            if(!savedGeometry)savedGeometry=GramGeometryStatistics::prepare(e);
            geometry=*savedGeometry;
        }
        if(!geometry.valid)return std::nullopt;
        const auto maximumIncidence=geometry.maximumIncidence,maximumWidth=geometry.maximumWidth,total=geometry.total;
        const double maximumEntry=geometry.maximumEntry,maximumWeight=data_->gramMaximumWeight;
        if(total==0||maximumEntry==0||maximumWeight==0)return 0.;
        if(maximumWidth==0||total>exactLimit/maximumWidth)return std::nullopt;
        const std::size_t pairs=total*maximumWidth;
        if(maximumIncidence>exactLimit-2)return std::nullopt;
        const double qe=gramUpperMultiply(static_cast<double>(maximumIncidence+2),std::numeric_limits<double>::epsilon());
        if(!(qe<1))return std::nullopt;
        const double denominator=std::nextafter(1-qe,-std::numeric_limits<double>::infinity());
        if(!std::isfinite(denominator)||denominator<=0)return std::nullopt;
        const double beta=gramUpperSuccessor(1/denominator),gamma=gramUpperMultiply(qe,beta);
        const double relative=gramUpperMultiply(gramUpperMultiply(gamma,maximumWeight),gramUpperMultiply(maximumEntry,maximumEntry));
        const double absolute=gramUpperMultiply(gramUpperMultiply(beta,std::numeric_limits<double>::denorm_min()),
            gramUpperAdd(gramUpperMultiply(2.,maximumEntry),1.));
        return gramUpperMultiply(gramUpperMultiply(maximumScale,maximumScale),
            gramUpperMultiply(static_cast<double>(pairs),gramUpperAdd(relative,absolute)));
    }catch(const GramBoundUnavailable&){return std::nullopt;}
    catch(const std::bad_alloc&){return std::nullopt;}
    catch(const std::length_error&){return std::nullopt;}
    catch(const std::system_error&){return std::nullopt;}
}
std::optional<std::vector<double>> HydroEnergyViscosity::scaledGramLowerErrorBounds(
    const std::vector<double>& positiveScale) const
{
    const GramBoundScope saved;
    // Capturing first also preserves invalid/NaN-input query side effects.
    if(!saved.captured||!data_->gramSourceQualified)return std::nullopt;
    StressReuseEnvironment caller;
    if(!caller.read()||!caller.sameControl(data_->gramBuildExit))return std::nullopt;
    try {
        const auto& e=data_->strain.matrix();
        if(!e.isCompressed()||e.rows()<0||e.cols()<0
            ||static_cast<std::size_t>(e.cols())!=positiveScale.size()
            ||static_cast<std::size_t>(e.rows())!=data_->weights.size()
            ||data_->matrix->rows()!=e.cols()||data_->matrix->cols()!=e.cols())return std::nullopt;
        const auto n=positiveScale.size();
        for(double s:positiveScale)if(!std::isfinite(s)||s<=0)return std::nullopt;
        std::vector<std::size_t> counts(n,0);
        std::vector<double> relative(n,0),absolute(n,0),bound(n,0);
        // n_i bounds m_ij for BOTH actual orientations. Exact-zero-weight
        // rows generate exact zero contributions and can be omitted.
        for(Eigen::Index row=0;row<e.outerSize();++row) {
            const double w=data_->weights[static_cast<std::size_t>(row)];
            if(!std::isfinite(w)||w<0)return std::nullopt;
            double a=0,s=0;
            for(Matrix::InnerIterator entry(e,row);entry;++entry) {
                if(entry.col()<0||static_cast<std::size_t>(entry.col())>=n
                    ||!std::isfinite(entry.value()))return std::nullopt;
                if(w==0||entry.value()==0)continue;
                const auto i=static_cast<std::size_t>(entry.col());
                // q=n_i+2 must remain an exactly representable integer.
                if(counts[i]>=4503599627370494ULL)return std::nullopt;
                ++counts[i];
                const double magnitude=std::abs(entry.value());
                a=gramUpperAdd(a,gramUpperMultiply(magnitude,positiveScale[i]));
                s=gramUpperAdd(s,positiveScale[i]);
            }
            if(w==0||s==0)continue;
            const double weightedA=gramUpperMultiply(w,a);
            const double sharedAbsolute=gramUpperAdd(s,a);
            for(Matrix::InnerIterator entry(e,row);entry;++entry) {
                if(entry.value()==0)continue;
                const auto i=static_cast<std::size_t>(entry.col());
                const double magnitude=std::abs(entry.value());
                relative[i]=gramUpperAdd(relative[i],gramUpperMultiply(magnitude,weightedA));
                // (|e_i|+1)*S+A = |e_i|*S+(S+A), with no cancellation.
                absolute[i]=gramUpperAdd(absolute[i],
                    gramUpperAdd(gramUpperMultiply(magnitude,s),sharedAbsolute));
            }
        }
        const double eps=std::numeric_limits<double>::epsilon();
        const double eta=std::numeric_limits<double>::denorm_min();
        for(std::size_t i=0;i<n;++i) {
            if(counts[i]==0)continue; // exact zero row; no nextafter(0) inflation.
            const double q=static_cast<double>(counts[i]+2);
            const double qe=gramUpperMultiply(q,eps);
            if(!(qe<1))return std::nullopt;
            const double denominator=std::nextafter(1-qe,-std::numeric_limits<double>::infinity());
            if(!std::isfinite(denominator)||!(denominator>0))return std::nullopt;
            const double beta=gramUpperSuccessor(1/denominator);
            const double gamma=gramUpperMultiply(qe,beta);
            // eta/2 is NEVER evaluated: that value is not binary64.
            const double underflow=gramUpperMultiply(gramUpperMultiply(beta,eta),absolute[i]);
            bound[i]=gramUpperMultiply(positiveScale[i],
                gramUpperAdd(gramUpperMultiply(gamma,relative[i]),underflow));
        }
        return bound;
    }catch(const GramBoundUnavailable&) {return std::nullopt;}
    catch(const std::bad_alloc&) {return std::nullopt;}
    catch(const std::length_error&) {return std::nullopt;}
}
HydroEnergyViscosity HydroEnergyViscosity::zero(const QuadtreeMesh& mesh,HydroPhysicalDivergence metrics)
{
    metrics.layout().validateMesh(mesh);
    auto strain=std::make_shared<HydroStrainOperator::Data>(std::move(metrics));
    strain->components=strain->metrics.layout().domain().axisymmetric()?4:3;
    strain->matrix.resize(0,strain->metrics.layout().velocityCount());
    strain->matrix.makeCompressed();
    return build(HydroStrainOperator(std::move(strain)),{});
}
HydroEnergyViscosity HydroEnergyViscosity::build(HydroStrainOperator strain,const std::vector<double>& viscosity)
{
    return buildImpl(std::move(strain),viscosity,nullptr);
}
HydroEnergyViscosity HydroEnergyViscosity::buildImpl(HydroStrainOperator strain,const std::vector<double>& viscosity,
    const HydroNoSlipPreparationReport* noSlip)
{
    require(!strain.data_->rt0GreenOnly,"RT0 Green geometry has no requested viscosity strain");
    require(viscosity.size()==strain.points().size(),"viscosity dimensions do not match quadrature");
    StressReuseEnvironment gramEntry;
    const bool gramEntryQualified=gramProductImplementationKnown()&&gramEntry.read();
    auto data=std::make_shared<Data>(std::move(strain));data->viscosity=viscosity;
    if(noSlip)data->product.noSlip=*noSlip;
    for (std::size_t q=0;q<viscosity.size();++q) {
        require(std::isfinite(viscosity[q]) && viscosity[q]>=0,"negative/nonfinite viscosity");
        const long double product=2*static_cast<long double>(viscosity[q])*data->strain.points()[q].physicalWeight;
        const double weight=static_cast<double>(product);
        require(std::isfinite(weight),"nonfinite viscosity quadrature weight");
        for (int k=0;k<data->strain.components();++k) data->weights.push_back(weight);
        data->rememberGramWeight(weight); // The same validated value just appended at every component.
    }
    ViscosityProductPlan::NumericPendingHandle numericPending;
    const auto& e=data->strain.matrix();
    if (std::any_of(data->weights.begin(),data->weights.end(),[](double x) {return x!=0;})) {
        bool assembled=false;std::size_t preparationPeakBytes=0;
        if(!data->strain.isStressCompatible()&&viscosityProductPlanEnvironment()){
            data->product.requested=true;data->product.fullAssembly=true;data->product.weightRows=data->weights.size();
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
            std::shared_ptr<const ViscosityProductPlan> plan;
            const int setupErrno=errno;StrainFloatingState setup;unsigned control=0;
            const bool qualified=ViscosityProductPlan::numericEnvironment(setup,control);
            if(qualified){
                try{
                    const std::lock_guard<std::mutex> guard(data->strain.data_->viscosityProductMutex);
                    plan=data->strain.data_->viscosityProductPlan;
                    if(!plan){
                        plan=std::make_shared<const ViscosityProductPlan>(e);
                        preparationPeakBytes=plan->preparationPeakBytes;
                        data->product.symbolicBuilt=true;data->product.symbolicWorkers=plan->symbolicWorkers;
                        data->product.columnRoutingSeconds=plan->columnRoutingSeconds;
                        data->product.symbolicGraphSeconds=plan->symbolicGraphSeconds;
                        data->strain.data_->viscosityProductPlan=plan;
                    }
                }catch(const std::length_error&){data->product.fallbackReason=3;}
                catch(const std::bad_alloc&){data->product.fallbackReason=3;}
                catch(...){data->product.fallbackReason=5;}
            }else data->product.fallbackReason=2;
            if(setup.valid)setup.restore();errno=setupErrno;
            if(plan){
                ViscosityPreparationClock numericClock;
                numericPending=plan->assembleRemembered(e,data->weights,data->matrix,data->product);
                data->product.numericAssemblySeconds=numericClock.seconds();
                data->product.timingAvailable=numericClock.available&&(!data->product.symbolicBuilt||plan->timingAvailable);
                assembled=static_cast<bool>(numericPending);
            }
#else
            data->product.fallbackReason=2;
#endif
        }
        data->product.peakOptionalBytes=std::max(data->product.peakOptionalBytes,preparationPeakBytes);
        if(!assembled){
            // The complete original numerical path is the only fallback.
            auto actual=std::make_shared<Matrix>();actual->resize(e.cols(),e.cols());
            Matrix weighted=e;
            for (Eigen::Index row=0;row<weighted.outerSize();++row)
                for (Matrix::InnerIterator entry(weighted,row);entry;++entry) entry.valueRef()*=data->weights[static_cast<std::size_t>(row)];
            *actual=e.transpose()*weighted;
            actual->prune(0.0);actual->makeCompressed();
            for(Eigen::Index i=0;i<actual->nonZeros();++i)require(std::isfinite(actual->valuePtr()[i]),"nonfinite assembled viscosity");
            data->matrix=std::move(actual);
        }
    }else{
        auto actual=std::make_shared<Matrix>();actual->resize(e.cols(),e.cols());actual->makeCompressed();data->matrix=std::move(actual);
    }
    data->gramBuildEntry=gramEntry;
    data->gramSourceQualified=gramEntryQualified&&data->gramBuildExit.read()
        &&gramEntry.sameControl(data->gramBuildExit);
    if(numericPending)numericPending->publish(data->gramSourceQualified,data->product);
    return HydroEnergyViscosity(std::move(data));
}
#if defined(_MSC_VER) && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST)
#pragma float_control(precise,on,push)
#pragma float_control(except,on)
#pragma fenv_access(on)
#endif
HydroEnergyViscosity HydroEnergyViscosity::buildStationaryNoSlipOwned(const QuadtreeMesh& mesh,
    const HydroMomentumGeometry& owner,HydroPhysicalDivergence metrics,
    const std::vector<HydroStrainPoint>& points,const std::vector<int>& walls,
    const std::vector<double>& viscosity)
{
    // The caller has already executed the original live phase/wall/current-mu
    // checks and volumeQuadrature. The physical cold factory stays at exactly
    // its old position; optional work never catches its exceptions.
    HydroNoSlipPreparationReport preparation;preparation.available=true;
    const auto cold=[&]{
        ViscosityPreparationClock clock;
        auto raw=HydroStrainOperator::buildStationaryNoSlipImpl(mesh,metrics,points,walls,&preparation);
        preparation.coldBuilt=true;preparation.coldSeconds=clock.seconds();preparation.timingAvailable=clock.available;
        return raw;
    };
    const auto literal=[&]{return buildImpl(cold(),viscosity,&preparation);};
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE) \
    && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST)
    if(!owner.layout().sameIdentity(metrics.layout())
        || &owner.metrics().divergence().matrix()!=&metrics.divergence().matrix())return literal();
    const auto cache=owner.stationaryNoSlipStrainCache();if(!cache)return literal();
    // Reuse the original raw worker state representation. Each small optional
    // lookup/admission/cleanup transaction restores its CURRENT entry, including
    // errno; it cannot erase flags from the intervening physical raw/K builds.
    const auto passive=[](const auto& function) noexcept {
        const int error=errno;StrainFloatingState caller;
        if(!caller.capture()){errno=error;return;}
        try{function();}catch(...){}
        caller.restore();errno=error;
    };
    StressReuseEnvironment sourceEntry;bool sourceQualified=false;int sourceErrno=errno;
    std::optional<HydroStrainOperator> retained;
    passive([&]{
        sourceErrno=errno;sourceQualified=sourceEntry.read();
        std::unique_lock<std::mutex> lock(cache->mutex,std::try_to_lock);
        if(!lock.owns_lock())return;
        bool same=sourceQualified&&cache->raw&&sourceErrno==cache->error
            &&sourceEntry.same(cache->environment);
        if(same){
            const auto& raw=*cache->raw;
            raw.layout().validateMesh(mesh);
            same=!raw.isStressCompatible()&&raw.layout().sameIdentity(metrics.layout())
                && &raw.metrics().divergence().matrix()==&metrics.divergence().matrix()
                &&raw.stationaryNoSlipVelocityDofs()==walls&&raw.points().size()==points.size();
            // The current builtin quadrature remains literal every call.
            // Compare exact value bits rather than struct padding/tolerances.
            for(std::size_t i=0;same&&i<points.size();++i){const auto& a=raw.points()[i];const auto& b=points[i];
                same=a.cellDof==b.cellDof&&std::memcmp(&a.point.x,&b.point.x,sizeof(double))==0
                    &&std::memcmp(&a.point.y,&b.point.y,sizeof(double))==0
                    &&std::memcmp(&a.physicalWeight,&b.physicalWeight,sizeof(double))==0;}
            if(same){retained.emplace(raw);return;}
        }
        // Drop an obsolete optional strong stock before a cold physical build.
        // Other actual users may still own their immutable raw; this is only
        // the cache's logical payload bound, not a process-RSS guarantee.
        cache->raw.reset();cache->retainedBytes=0;
    });
    if(retained){preparation.rawCacheHit=true;preparation.timingAvailable=true;
        return buildImpl(*retained,viscosity,&preparation);} // ALL current weights/K refresh.
    auto raw=cold();
    bool quiet=false;
    passive([&]{
        StressReuseEnvironment exit;const int error=errno;
        quiet=sourceQualified&&exit.read()&&error==sourceErrno&&sourceEntry.same(exit);
    });
    auto result=buildImpl(raw,viscosity,&preparation); // Original prune/finite/source gates first.
    if(quiet)passive([&]{
        using Cache=Detail::HydroStationaryNoSlipStrainCache;
        // Quiet immutable raw geometry may be retained before any nonzero K.
        // The full future reservation below covers a later plan and numeric work.
        std::unique_lock<std::mutex> planLock(raw.data_->viscosityProductMutex,std::try_to_lock);
        if(!planLock.owns_lock())return;
        const auto plan=raw.data_->viscosityProductPlan;
        std::size_t bytes=sizeof(Cache)+1536u;
        const auto add=[&](std::size_t count,std::size_t width){
            if(bytes>Cache::byteLimit||count>(Cache::byteLimit-bytes)/width)return false;
            bytes+=count*width;return true;
        };
        const auto sparse=[&](const auto& matrix){
            if(!matrix.isCompressed()||matrix.rows()<0||matrix.cols()<0||matrix.data().allocatedSize()<0)return false;
            const auto allocated=static_cast<std::size_t>(matrix.data().allocatedSize());
            if(allocated<static_cast<std::size_t>(matrix.nonZeros()))return false;
            using Index=typename std::decay_t<decltype(matrix)>::StorageIndex;
            return add(allocated,sizeof(double)+sizeof(Index))
                &&add(static_cast<std::size_t>(matrix.outerSize())+1,sizeof(Index));
        };
        // Raw Data owns only its metrics/points/walls/CSR in this builtin
        // factory. Metrics is the SAME already-owned layout/D; no charge or
        // retention of another Momentum/phase owner is introduced.
        if(raw.data_->stressCompatible||raw.data_->canonicalRaw||raw.data_->preparation
            ||!raw.data_->facePoints.empty()||raw.data_->faceVelocity.nonZeros()!=0
            ||raw.data_->volumeVelocity.nonZeros()!=0||raw.data_->componentIntegrals.nonZeros()!=0)return;
        if(!add(1,sizeof(HydroStrainOperator::Data))||!sparse(raw.data_->matrix)
            ||!add(raw.data_->points.capacity(),sizeof(HydroStrainPoint))
            ||!add(raw.data_->stationaryNoSlipWalls.capacity(),sizeof(int)))return;
        if(plan&&(!add(1,sizeof(ViscosityProductPlan))
            ||!add(plan->columnOffsets.capacity(),sizeof(ViscosityProductPlan::Route))
            ||!add(plan->strainRows.capacity(),sizeof(ViscosityProductPlan::Route))
            ||!add(plan->entryOffsets.capacity(),sizeof(ViscosityProductPlan::Route))
            ||!add(plan->outer.capacity(),sizeof(ViscosityProductPlan::Route))
            ||!add(plan->inner.capacity(),sizeof(ViscosityProductPlan::Route))))return;
        // No-plan retention reserves the complete future route+numeric scope.
        // With a plan already present, its actual storage plus this full reserve
        // is deliberately conservative; this is not a minimal process peak.
        if(!add(Cache::numericReservation,1))return;
        std::unique_lock<std::mutex> lock(cache->mutex,std::try_to_lock);
        if(!lock.owns_lock()||cache->raw)return;
        cache->environment=sourceEntry;cache->error=sourceErrno;
        cache->raw.emplace(raw);cache->retainedBytes=bytes;
    });
    return result;
#else
    (void)owner;return literal();
#endif
}
#if defined(_MSC_VER) && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST)
#if !defined(_M_FP_STRICT) || !_M_FP_STRICT
#pragma fenv_access(off)
#endif
#pragma float_control(pop)
#endif

std::vector<double> HydroEnergyViscosity::apply(const QuadtreeMesh& mesh,const HydroNativeValues& values) const
{
    validate(mesh,values,strain().metrics());return multiply(matrix(),values.normalVelocity());
}
double HydroEnergyViscosity::dissipation(const QuadtreeMesh& mesh,const HydroNativeValues& values) const
{
    if (std::none_of(rowWeights().begin(),rowWeights().end(),[](double x) {return x!=0;})) {
        validate(mesh,values,strain().metrics());return 0;
    }
    const auto e=strain().apply(mesh,values);long double total=0;
    for (std::size_t i=0;i<e.size();++i) total+=static_cast<long double>(rowWeights()[i])*e[i]*e[i];
    const double result=static_cast<double>(total);require(std::isfinite(result),"nonfinite viscous dissipation");return result;
}
}

