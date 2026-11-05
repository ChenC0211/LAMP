#include "Solvers/Poisson/PolarizationPoissonFVSolverCellCentered.h"

#include <algorithm>
#include <atomic>
#include <cfenv>
#include <cerrno>
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
#include <omp.h>
#include <float.h>
#include <xmmintrin.h>
#endif
#include <chrono>
#include <cmath>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <new>
#include <queue>
#include <sstream>

#include <Eigen/SparseLU>
#include <Eigen/IterativeLinearSolvers>

#ifdef EIGEN_SUPERLU_SUPPORT
#include <Eigen/SuperLUSupport>
#endif

#include "Core/IO/ConsoleControl.h"
#include "Core/Geometry/Detail/FiniteVolumeClosure.h"
#include "Core/Mesh/PeriodicFaceUtils.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/NativeMeshIdentity.h"
#include "Core/Mesh/CartesianMeasure.h"
#include "Solvers/Poisson/CutCellMetrics.h"
#include "Solvers/Linear/ParallelBlockJacobi.h"
#include "Solvers/Poisson/detail/PolarizationFluxCorrectionValidation.h"
#include "Solvers/Poisson/detail/PolarizationFittedBoundaryValidation.h"
#include "Solvers/Poisson/detail/PolarizationFittedRhsPlan.h"

// Like the fitted Hydro cache, this compares actual compressed coefficients,
// not a pattern hash. It never reuses the legacy matrix-freezing path.
struct PolarizationPoissonFVSolverCellCentered::FittedFactorCache
{
    Eigen::SparseMatrix<double> matrix;
    const QuadtreeMesh* mesh = nullptr; // identity only; never dereferenced
    std::vector<int> cellIndices, cellLabels;
    std::unique_ptr<const PolarizationFittedBoundaryDetail::RhsPlan> rhsPlan;
#ifdef EIGEN_SUPERLU_SUPPORT
    Eigen::SuperLU<Eigen::SparseMatrix<double>> solver;
#endif

    bool matches(const Eigen::SparseMatrix<double>& candidate,
        const QuadtreeMesh& candidateMesh, const std::vector<int>& rows) const
    {
        if (mesh != &candidateMesh || cellIndices != rows ||
            matrix.rows() != candidate.rows() || matrix.cols() != candidate.cols() ||
            !matrix.isCompressed() || !candidate.isCompressed() ||
            matrix.nonZeros() != candidate.nonZeros()) return false;
        for (std::size_t i = 0; i < rows.size(); ++i)
            if (cellLabels[i] != candidateMesh.cellIndexToLabel[rows[i]]) return false;
        using Index = Eigen::SparseMatrix<double>::StorageIndex;
        return std::memcmp(matrix.outerIndexPtr(), candidate.outerIndexPtr(),
                   (candidate.outerSize() + 1) * sizeof(Index)) == 0 &&
               std::memcmp(matrix.innerIndexPtr(), candidate.innerIndexPtr(),
                   candidate.nonZeros() * sizeof(Index)) == 0 &&
               std::memcmp(matrix.valuePtr(), candidate.valuePtr(),
                   candidate.nonZeros() * sizeof(double)) == 0;
    }
};

PolarizationPoissonFVSolverCellCentered::PolarizationPoissonFVSolverCellCentered() = default;
PolarizationPoissonFVSolverCellCentered::~PolarizationPoissonFVSolverCellCentered() = default;

struct PolarizationPoissonFVSolverCellCentered::PreparedIntegratedStructure::Data
{
    explicit Data(const QuadtreeMesh& actual,std::unique_ptr<FaceSegmentGraph> owned)
        :identity(actual),graph(std::move(owned)){}
    NativeMeshIdentity identity;
    std::unique_ptr<const FaceSegmentGraph> graph;
    bool axisymmetric=false;
    std::vector<FittedBoundaryType> types;
    std::vector<std::vector<int>> correctionColumns,affineColumns;
    std::vector<int> rowOfCell,cellOfRow;
    std::vector<std::size_t> reverse,termOffsets,slots;
    std::vector<unsigned char> first;
    Eigen::SparseMatrix<double> pattern;
    std::vector<std::uint32_t> residualOffsets;
    std::vector<std::array<std::uint32_t,2>> residualEntries;
};

namespace
{
using Clock = std::chrono::steady_clock;

auto secondsBetween(const Clock::time_point& a, const Clock::time_point& b) -> double
{
    return std::chrono::duration_cast<std::chrono::duration<double>>(b - a).count();
}

// Current-call exact integrated rows only; no new solver/backend or cache.
struct IntegratedRowReceipt
{
    bool available=false;
    std::size_t defectPhases=0,projectionPhases=0,parallelDefectRows=0,parallelProjectionRows=0,fallbackPhases=0;
    unsigned maximumWorkers=0;
};
enum class IntegratedRowStatus {Literal,Parallel,Fallback};
struct ExpandedResidualReceipt
{
    bool available=false,parallelUsed=false,parallelFallback=false,timingAvailable=false;
    unsigned workers=0;
    std::size_t rows=0,scratchBytes=0;
    double routeSeconds=0,seconds=0;
};
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
// These private functions alone have strict local FP semantics. The old TU
// and literal algorithms retain their original compiler options. Explicit
// products are stored before exact conversion; no fused add/product exists.
#pragma float_control(precise,on,push)
#pragma float_control(except,on)
#pragma fenv_access(on)
struct IntegratedRowFloatingState
{
    std::fenv_t native{};unsigned csr=0;bool valid=false;
    bool capture() noexcept {csr=_mm_getcsr();valid=std::fegetenv(&native)==0;return valid;}
    void restore() const noexcept
    {if(valid&&std::fesetenv(&native)!=0)std::terminate();_mm_setcsr(csr);}
    bool enterFresh(const IntegratedRowFloatingState& from) const noexcept
    {auto clean=from.native;clean._Fe_stat=0;if(std::fesetenv(&clean)!=0)return false;
        _mm_setcsr(from.csr&~0x3fu);return true;}
};
bool integratedRowSameControl(const IntegratedRowFloatingState& a,const IntegratedRowFloatingState& b) noexcept
{
    auto x=a.native,y=b.native;x._Fe_stat=0;y._Fe_stat=0;
    return std::memcmp(&x,&y,sizeof x)==0&&(a.csr&~0x3fu)==(b.csr&~0x3fu);
}
bool integratedRowQualified(IntegratedRowFloatingState& caller) noexcept
{
#if !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST)
    if(omp_in_parallel()||omp_get_max_threads()<2||!caller.capture())return false;
    unsigned control=0;
    if(_controlfp_s(&control,0,0)!=0||(control&_MCW_EM)!=_MCW_EM
        ||(control&_MCW_RC)!=_RC_NEAR||std::fegetround()!=FE_TONEAREST
        ||(caller.csr&0x1f80u)!=0x1f80u||(caller.csr&0xe040u)!=0)return false;
#ifdef _EM_AMBIGUOUS
    if(control&_EM_AMBIGUOUS)return false;
#endif
    return sizeof(double)==8&&sizeof(long double)==8
        &&std::numeric_limits<double>::is_iec559&&std::numeric_limits<double>::digits==53
        &&std::numeric_limits<long double>::digits==53;
#else
    (void)caller;return false;
#endif
}
// Declare before optional phase storage and arm before no-fail commits;
// scratch destruction cannot leak allocator FP/errno into the caller.
struct IntegratedRowCleanup
{
    IntegratedRowFloatingState finalState;int error=0;bool armed=false;
    bool arm() noexcept
    {error=errno;armed=finalState.capture();errno=error;return armed;}
    ~IntegratedRowCleanup() noexcept
    {if(armed){finalState.restore();errno=error;}}
};
template<class Kernel>
bool integratedRowRows(std::size_t count,int team,Kernel&& kernel,unsigned& actualWorkers)
{
    IntegratedRowFloatingState prepared;const int parentErrno=errno;
    if(!prepared.capture())return false;
    std::atomic<int> failed{0};unsigned long raisedNative=0;unsigned raisedCsr=0;
#pragma omp parallel num_threads(team)
    {
        IntegratedRowFloatingState worker;const int workerErrno=errno;
        const bool captured=worker.capture(),ready=captured&&worker.enterFresh(prepared);
        if(!ready)failed.store(1,std::memory_order_relaxed);
        if(ready)errno=parentErrno;
#pragma omp single
        {actualWorkers=std::max(actualWorkers,static_cast<unsigned>(omp_get_num_threads()));}
#pragma omp for schedule(static) nowait
        for(std::ptrdiff_t row=0;row<static_cast<std::ptrdiff_t>(count);++row){
            if(!ready)continue;
            try{kernel(static_cast<std::size_t>(row));}
            catch(...){failed.store(1,std::memory_order_relaxed);}
        }
        if(ready){
            const int numericalErrno=errno;IntegratedRowFloatingState after;
            const bool got=after.capture();
            if(!got||numericalErrno!=parentErrno||!integratedRowSameControl(after,prepared))
                failed.store(1,std::memory_order_relaxed);
#pragma omp critical(lamp_integrated_poisson_row_flags)
            {if(got){raisedNative|=after.native._Fe_stat;raisedCsr|=after.csr&0x3fu;}}
        }
        if(captured)worker.restore();errno=workerErrno;
    }
    if(failed.load(std::memory_order_relaxed)||actualWorkers<2)return false;
    auto merged=prepared.native;merged._Fe_stat|=raisedNative;
    if(std::fesetenv(&merged)!=0)std::terminate();
    _mm_setcsr(prepared.csr|raisedCsr);errno=parentErrno;return true;
}
// These observations preserve the caller's current flags/control and errno.
// They never surround or erase a numerical operation from the original gate.
struct IntegratedCurrentEnvironment
{
    IntegratedRowFloatingState current;int error=errno;
    IntegratedCurrentEnvironment() noexcept {current.capture();errno=error;}
    ~IntegratedCurrentEnvironment() noexcept {if(current.valid)current.restore();errno=error;}
};
__declspec(noinline) bool integratedClockSample(Clock::time_point& result) noexcept
{
    IntegratedCurrentEnvironment restore;if(!restore.current.valid)return false;
    unsigned control=0;
    if(_controlfp_s(&control,0,0)!=0||(control&_MCW_EM)!=_MCW_EM||(restore.current.csr&0x1f80u)!=0x1f80u)return false;
    result=Clock::now();return true;
}
__declspec(noinline) bool integratedClockElapsed(const Clock::time_point& start,double& result) noexcept
{
    IntegratedCurrentEnvironment restore;if(!restore.current.valid)return false;
    unsigned control=0;
    if(_controlfp_s(&control,0,0)!=0||(control&_MCW_EM)!=_MCW_EM||(restore.current.csr&0x1f80u)!=0x1f80u)return false;
    const auto stop=Clock::now();const double value=std::chrono::duration_cast<std::chrono::duration<double>>(stop-start).count();
    if(!std::isfinite(value)||value<0)return false;result=value;return true;
}
__declspec(noinline) double integratedExpandedResidualRow(std::size_t row,
    const Eigen::SparseMatrix<double>& matrix,const std::vector<std::uint32_t>& offsets,
    const std::vector<std::array<std::uint32_t,2>>& routes,
    const std::vector<std::array<double,3>>& potential,const Eigen::VectorXd& rhs)
{
    namespace Exact=LAMP::Geometry::Detail::FiniteVolume;
    Exact::Dyadic exact;exact.add(rhs[static_cast<Eigen::Index>(row)],-1);
    for(auto at=offsets[row];at<offsets[row+1];++at){const auto& item=routes[at];
        for(double component:potential[item[0]])
            exact.add(Exact::Dyadic::product({matrix.valuePtr()[item[1]],component}));}
    return Exact::roundDouble(exact).value;
}
__declspec(noinline) double integratedDefectRow(std::size_t compact,const QuadtreeMesh& mesh,
    const FaceSegmentGraph& graph,const std::vector<PolarizationPoissonFVSolverCellCentered::IntegratedFluxFactor>& factors,
    const std::vector<std::array<double,3>>& expandedPotential,const std::vector<int>& rowOfCellIndex,
    const std::vector<int>& cellIndexOfRow,const CellField<double>& rhs)
{
    namespace Exact=LAMP::Geometry::Detail::FiniteVolume;
    using Transfer=PolarizationPoissonFVSolverCellCentered::IntegratedFluxFactor::SecondaryTransfer;
    const auto& segs=graph.segments();
                    const auto row=static_cast<int>(compact);const int cell=cellIndexOfRow[row];
                    const auto& offsets=graph.cellOffsets();
                    if(cell<0||static_cast<std::size_t>(cell)+1>=offsets.size())
                        throw std::runtime_error("optional integrated row native range");
                    const auto first=offsets[cell],last=offsets[cell+1];
                    if(first>last||last>segs.size())throw std::runtime_error("optional integrated row SID range");
                    Exact::Dyadic exact;exact.add(rhs.value[cell]);
                    for(auto sid=first;sid<last;++sid){
                        if(mesh.cells(segs[sid].K).index()!=cell)throw std::runtime_error("optional integrated row owner mismatch");
                        const auto& factor=factors[sid];Exact::Dyadic h(factor.affine.constant);
                        for(const auto& coefficient:factor.affine.coefficients)
                            for(double component:expandedPotential[rowOfCellIndex[coefficient.cellIndex]])
                                h.add(Exact::Dyadic::product({coefficient.value,component}));
                        const double value=static_cast<double>(h.rounded());
                        const double primary=static_cast<double>(static_cast<long double>(factor.primary)*value);
                        const double secondary=factor.secondary==0?0:static_cast<double>(static_cast<long double>(factor.secondary)*value);
                        if(!(std::isfinite(value)&&std::isfinite(primary)&&std::isfinite(secondary)))
                            throw std::runtime_error("nonfinite uncompressed integrated face action");
                        exact.add(primary,-1);auto transfer=Exact::Dyadic::product({factor.weight,secondary});
                        if(factor.secondaryTransfer==Transfer::MinimumBinary64QuantumNearestEven)
                            transfer=transfer.roundedBinary64Quantum();
                        exact.add(transfer,-1);
                    }
                    return Exact::roundDouble(exact).value;
}
__declspec(noinline) double integratedProjectionRow(std::size_t row,int pass,
    const std::vector<std::array<double,3>>& expandedPotential,const Eigen::VectorXd& delta)
{
    namespace Exact=LAMP::Geometry::Detail::FiniteVolume;
                    auto terms=expandedPotential[row];terms[pass+1]=delta[static_cast<int>(row)];
                    Exact::Dyadic sum;for(double component:terms)sum.add(component);
                    return Exact::roundDouble(sum).value;
}
#if !defined(_M_FP_STRICT) || !_M_FP_STRICT
#pragma fenv_access(off)
#endif
#pragma float_control(pop)

#endif
template<class Kernel,class Commit>
IntegratedRowStatus tryIntegratedRows(std::size_t count,Kernel&& kernel,Commit&& commit,unsigned& workers,
    std::size_t retainedBytes=0,std::size_t* scratchBytes=nullptr)
{
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
    const int error=errno;IntegratedRowFloatingState entry;const bool qualified=count>1&&integratedRowQualified(entry);
    if(entry.valid)entry.restore();errno=error;
    if(!qualified)return IntegratedRowStatus::Literal;
    try{
        IntegratedRowCleanup cleanup;
        constexpr std::size_t limit=64u*1024u*1024u,workerReserve=1024u*1024u;
        const auto teamCount=std::min<std::size_t>(count,static_cast<std::size_t>(omp_get_max_threads()));
        if(teamCount+1>limit/workerReserve||count>static_cast<std::size_t>(std::numeric_limits<std::ptrdiff_t>::max()))
            throw std::length_error("optional integrated row worker range");
        const auto reserve=(teamCount+1)*workerReserve;
        // Count all three coexisting double buffers (phase values, retained
        // parallelDefect and original defect) at projection peak. The common
        // original exactDefect stock/owners remain outside this extra cap.
        if(retainedBytes>limit-reserve||count>(limit-reserve-retainedBytes)/(3*sizeof(double)))
            throw std::length_error("optional integrated row scratch cap");
        std::vector<double> values(count);
        if(scratchBytes){
            const auto extra=sizeof(values)+2*sizeof(Eigen::VectorXd);
            if(extra>limit-reserve-retainedBytes||
                values.capacity()>(limit-reserve-retainedBytes-extra)/(sizeof(double))||
                2*count>(limit-reserve-retainedBytes-extra-values.capacity()*sizeof(double))/sizeof(double))
                throw std::length_error("optional integrated row actual capacity cap");
            *scratchBytes=std::max(*scratchBytes,retainedBytes+reserve+extra+
                (values.capacity()+2*count)*sizeof(double));
        }
        if(errno!=error)throw std::runtime_error("optional integrated row allocation errno");
        entry.restore();errno=error;
        unsigned actual=0;
        if(!integratedRowRows(count,static_cast<int>(teamCount),[&](std::size_t row){values[row]=kernel(row);},actual))
            throw std::runtime_error("optional integrated row worker unavailable");
        if(!cleanup.arm())throw std::runtime_error("optional integrated row final FP receipt");
        commit(values);workers=actual;return IntegratedRowStatus::Parallel;
    }catch(...){entry.restore();errno=error;return IntegratedRowStatus::Fallback;}
#else
    (void)count;(void)kernel;(void)commit;(void)workers;(void)retainedBytes;(void)scratchBytes;return IntegratedRowStatus::Literal;
#endif
}

struct ExpandedResidualRoute
{
    static constexpr std::size_t limit=64u*1024u*1024u,workerReserve=1024u*1024u;
    std::vector<std::uint32_t> offsets;
    std::vector<std::array<std::uint32_t,2>> entries; // column, actual CSC value offset
    const std::vector<std::uint32_t>* borrowedOffsets=nullptr;
    const std::vector<std::array<std::uint32_t,2>>* borrowedEntries=nullptr;
    const auto& rowOffsets() const noexcept {return borrowedOffsets?*borrowedOffsets:offsets;}
    const auto& rowEntries() const noexcept {return borrowedEntries?*borrowedEntries:entries;}
    void borrow(const std::vector<std::uint32_t>& rows,const std::vector<std::array<std::uint32_t,2>>& routes) noexcept
    {borrowedOffsets=&rows;borrowedEntries=&routes;ready=true;peakBytes=retainedBytes();}
    bool attempted=false,counting=false,ready=false,failed=false,timingAvailable=true;
    std::size_t peakBytes=0;
    double seconds=0;
    Clock::time_point started{};
    std::size_t retainedBytes() const noexcept
    {return sizeof(*this)+rowOffsets().capacity()*sizeof(std::uint32_t)+rowEntries().capacity()*sizeof(entries[0]);}
    void discard() noexcept
    {std::vector<std::uint32_t>().swap(offsets);std::vector<std::array<std::uint32_t,2>>().swap(entries);counting=ready=false;}
    ~ExpandedResidualRoute() noexcept
    {
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
        if(!offsets.empty()||!entries.empty()){IntegratedCurrentEnvironment restore;discard();}
#endif
    }
    void begin(const Eigen::SparseMatrix<double>& matrix)
    {
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
        IntegratedCurrentEnvironment restore;if(!restore.current.valid||matrix.rows()<2)return;
        IntegratedRowFloatingState entry;if(!integratedRowQualified(entry))return;
        attempted=true;timingAvailable=integratedClockSample(started);
        try{
            const auto n=static_cast<std::size_t>(matrix.rows()),nnz=static_cast<std::size_t>(matrix.nonZeros());
            if(!matrix.isCompressed()||matrix.rows()!=matrix.cols()||n<2||
                n>=std::numeric_limits<std::uint32_t>::max()||nnz>std::numeric_limits<std::uint32_t>::max())
                throw std::length_error("optional expanded residual index range");
            const auto team=std::min(n,static_cast<std::size_t>(omp_get_max_threads()));
            if(team+1>limit/workerReserve||n>limit/(3*sizeof(double)+2*sizeof(std::uint32_t))||
                nnz>limit/sizeof(entries[0]))throw std::length_error("optional expanded residual route cap");
            const auto expected=sizeof(*this)+(2*n+1)*sizeof(std::uint32_t)+nnz*sizeof(entries[0])+sizeof(std::vector<std::uint32_t>);
            const auto future=(team+1)*workerReserve+3*n*sizeof(double)+sizeof(std::vector<double>)+2*sizeof(Eigen::VectorXd);
            if(expected>limit||future>limit-expected)throw std::length_error("optional expanded residual combined cap");
            offsets.resize(n+1,0);peakBytes=retainedBytes();
            if(errno!=restore.error)throw std::runtime_error("optional expanded residual allocation errno");
            counting=true;
        }catch(...){failed=true;discard();}
#else
        (void)matrix;
#endif
    }
    void count(Eigen::Index row) noexcept
    {
        if(!counting)return;
        if(row<0||static_cast<std::size_t>(row)+1>=offsets.size()||offsets[static_cast<std::size_t>(row)+1]==std::numeric_limits<std::uint32_t>::max())
        {failed=true;counting=false;return;}
        ++offsets[static_cast<std::size_t>(row)+1];
    }
    void finish(const Eigen::SparseMatrix<double>& matrix)
    {
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
        if(!attempted)return;
        IntegratedCurrentEnvironment restore;
        if(!restore.current.valid){failed=true;counting=false;return;}
        if(counting)try{
            const auto n=static_cast<std::size_t>(matrix.rows()),nnz=static_cast<std::size_t>(matrix.nonZeros());
            std::uint64_t total=0;
            for(std::size_t row=0;row<n;++row){total+=offsets[row+1];
                if(total>std::numeric_limits<std::uint32_t>::max())throw std::length_error("optional expanded residual prefix range");
                offsets[row+1]=static_cast<std::uint32_t>(total);}
            if(total!=nnz)throw std::runtime_error("optional expanded residual incomplete scan");
            std::vector<std::uint32_t> cursor(offsets.begin(),offsets.end()-1);
            entries.resize(nnz);
            const auto bytes=retainedBytes()+sizeof(cursor)+cursor.capacity()*sizeof(cursor[0]);
            peakBytes=std::max(peakBytes,bytes);
            const auto team=std::min(n,static_cast<std::size_t>(omp_get_max_threads()));
            const auto future=(team+1)*workerReserve+3*n*sizeof(double)+sizeof(std::vector<double>)+2*sizeof(Eigen::VectorXd);
            if(bytes>limit||future>limit-bytes||errno!=restore.error)
                throw std::length_error("optional expanded residual actual capacity cap");
            for(int column=0;column<matrix.outerSize();++column)
                for(auto at=matrix.outerIndexPtr()[column];at<matrix.outerIndexPtr()[column+1];++at){
                    const auto row=static_cast<std::size_t>(matrix.innerIndexPtr()[at]);
                    entries[cursor[row]++]={static_cast<std::uint32_t>(column),static_cast<std::uint32_t>(at)};}
            ready=true;counting=false;
        }catch(...){failed=true;discard();}
        else if(failed)discard();
        if(timingAvailable)timingAvailable=integratedClockElapsed(started,seconds);
#else
        (void)matrix;
#endif
    }
};

inline bool isDirichletLike(const PicBCType t)
{
    return (t == PicBCType::Dirichlet) || (t == PicBCType::Robin);
}

inline const PicBoundarySide& bcSide(const PicDomainBC& bc, const Boundary b)
{
    switch (b)
    {
    case Boundary::west:  return bc.left;
    case Boundary::east:  return bc.right;
    case Boundary::south: return bc.bottom;
    case Boundary::north: return bc.top;
    default:              return bc.left;
    }
}

inline const char* boundaryToString(const Boundary b)
{
    switch (b)
    {
    case Boundary::west:  return "west";
    case Boundary::south: return "south";
    case Boundary::east:  return "east";
    case Boundary::north: return "north";
    default:              return "non";
    }
}

inline double safeCellArea(const QuadtreeMesh& mesh, int cellIndex)
{
    if (cellIndex < 0 || cellIndex >= mesh.cellIndexNum())
        return 0.0;
    const int cLbl = mesh.cellIndexToLabel[cellIndex];
    const Cell& c = mesh.cells(cLbl);
    const double a = c.area();
    if (!std::isfinite(a) || !(a > 0.0))
        return 0.0;
    return a;
}

inline double stabilizedEmbeddedDistance(const double dGamma, const double cellArea)
{
    const double h = (std::isfinite(cellArea) && cellArea > 0.0) ? std::sqrt(cellArea) : 1.0;
    const double floor = std::max(1e-12, 2.0e-2 * h);
    if (!std::isfinite(dGamma))
        return floor;
    return std::max(dGamma, floor);
}

inline double clampKappa(double v)
{
    if (!std::isfinite(v))
        return 1.0;
    if (v < 0.0)
        return 0.0;
    return v;
}

inline std::uint64_t sparsePatternSignature(const Eigen::SparseMatrix<double>& A)
{
    constexpr std::uint64_t kOffset = 1469598103934665603ull;
    constexpr std::uint64_t kPrime = 1099511628211ull;
    auto mix = [kPrime](std::uint64_t h, std::uint64_t v) -> std::uint64_t {
        h ^= v;
        h *= kPrime;
        return h;
    };

    std::uint64_t h = kOffset;
    h = mix(h, static_cast<std::uint64_t>(A.rows()));
    h = mix(h, static_cast<std::uint64_t>(A.cols()));
    h = mix(h, static_cast<std::uint64_t>(A.nonZeros()));

    const Eigen::Index outerSize = A.outerSize();
    if (outerSize > 0)
    {
        const int* outer = A.outerIndexPtr();
        const Eigen::Index stride = std::max<Eigen::Index>(1, outerSize / 64);
        for (Eigen::Index i = 0; i < outerSize; i += stride)
            h = mix(h, static_cast<std::uint64_t>(outer[i]));
        h = mix(h, static_cast<std::uint64_t>(outer[outerSize]));
    }

    const Eigen::Index nnz = A.nonZeros();
    if (nnz > 0)
    {
        const int* inner = A.innerIndexPtr();
        const Eigen::Index stride = std::max<Eigen::Index>(1, nnz / 256);
        for (Eigen::Index k = 0; k < nnz; k += stride)
            h = mix(h, static_cast<std::uint64_t>(inner[k]));
        h = mix(h, static_cast<std::uint64_t>(inner[nnz - 1]));
    }

    return h;
}

inline bool sameSparsePattern(const Eigen::SparseMatrix<double>& a,
                              const Eigen::SparseMatrix<double>& b)
{
    if (a.rows() != b.rows() || a.cols() != b.cols() ||
        !a.isCompressed() || !b.isCompressed() || a.nonZeros() != b.nonZeros())
        return false;
    using Index = Eigen::SparseMatrix<double>::StorageIndex;
    return std::memcmp(a.outerIndexPtr(), b.outerIndexPtr(),
               (a.outerSize() + 1) * sizeof(Index)) == 0 &&
           std::memcmp(a.innerIndexPtr(), b.innerIndexPtr(),
               a.nonZeros() * sizeof(Index)) == 0;
}

inline bool sameSparseValues(const Eigen::SparseMatrix<double>& a,
                             const Eigen::SparseMatrix<double>& b)
{
    // Called only after exact pattern equality has been established.
    return std::memcmp(a.valuePtr(), b.valuePtr(), a.nonZeros() * sizeof(double)) == 0;
}

inline bool envBool(const char* name, bool fallback)
{
    if (!name)
        return fallback;
    std::string v;
#if defined(_MSC_VER)
    char* s = nullptr;
    size_t len = 0;
    if (_dupenv_s(&s, &len, name) != 0 || !s)
        return fallback;
    v.assign(s);
    std::free(s);
#else
    const char* s = std::getenv(name);
    if (!s)
        return fallback;
    v.assign(s);
#endif
    std::transform(v.begin(), v.end(), v.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (v == "1" || v == "true" || v == "yes" || v == "on")
        return true;
    if (v == "0" || v == "false" || v == "no" || v == "off")
        return false;
    return fallback;
}

// Shared by the full fitted/corrected solve and the owned-operator RHS solve.
// Keep the original complete residual expression and tolerance semantics.
bool fullResidualAccepted(const Eigen::SparseMatrix<double>& matrix, const Eigen::VectorXd& x,
    const Eigen::VectorXd& b, bool requiredBackendUsed, double relTol,
    double& relativeResidual, std::string& reason)
{
    const Eigen::VectorXd residual = matrix * x - b;
    const double rhsNorm = b.norm();
    const double residualNorm = residual.norm();
    if (!requiredBackendUsed || !x.allFinite() || !residual.allFinite() ||
        !std::isfinite(rhsNorm) || !std::isfinite(residualNorm))
    {
        reason = "fluxCorrections: non-finite solution/full-operator residual norm or required backend was not used";
        return false;
    }
    relativeResidual = rhsNorm > 0.0 ? residualNorm / rhsNorm : residualNorm;
    if (!std::isfinite(relativeResidual) || relativeResidual > relTol)
    {
        std::ostringstream message;
        message << "fluxCorrections: full-operator residual " << relativeResidual
                << " exceeds relTol=" << relTol;
        reason = message.str();
        return false;
    }
    return true;
}

// Called ONLY after the fitted validator has checked the complete actual
// active-leaf graph, fields, coefficients, boundary types and reverse geometry.
// The legacy flux validator deliberately remains uniform/full-domain-only;
// this supplement adds its algebraic contract without substituting its mesh.
void validateFittedFluxCorrections(const QuadtreeMesh& mesh, const FaceSegmentGraph& graph,
    const std::vector<PolarizationPoissonFVSolverCellCentered::FittedBoundary>& boundaries,
    const std::vector<PolarizationPoissonFVSolverCellCentered::SegmentFluxCorrection>& corrections)
{
    using Type = PolarizationPoissonFVSolverCellCentered::FittedBoundaryType;
    const auto require = [](bool condition, const char* message)
    {
        if (!condition) throw std::runtime_error(std::string("fitted fluxCorrections: ") + message);
    };
    const auto& segments = graph.segments();
    const auto& offsets = graph.cellOffsets();
    const int count = mesh.cellIndexNum();
    require(corrections.size() == segments.size(), "correction count must equal oriented segment count");
    for (std::size_t sid = 0; sid < segments.size(); ++sid)
    {
        const auto& correction = corrections[sid];
        require(std::isfinite(correction.constant), "non-finite correction constant");
        int previous = -1;
        for (const auto& coefficient : correction.coefficients)
        {
            const int index = coefficient.cellIndex;
            require(index > previous && index < count && std::isfinite(coefficient.value),
                "coefficients must have sorted unique active cell indices and finite values");
            require(offsets[index] < offsets[index + 1], "correction coefficient references an inactive cell index");
            previous = index;
        }
        if (boundaries[sid].type == Type::HomogeneousNeumann)
        {
            // The fitted geometry validator requires this same type on the
            // zero-area axis. Neither cancellation nor a small-value tolerance
            // can turn a prescribed zero-flux boundary into another operator.
            require(correction.constant == 0.0,
                "homogeneous Neumann boundary correction constant must be exactly zero");
            for (const auto& coefficient : correction.coefficients)
                require(coefficient.value == 0.0,
                    "homogeneous Neumann boundary correction coefficients must each be exactly zero");
        }
    }
    for (std::size_t sid = 0; sid < segments.size(); ++sid)
    {
        if (boundaries[sid].type != Type::Interior) continue;
        const auto& face = segments[sid];
        const int neighbor = mesh.cells(face.L).index(); // Already validated as active.
        std::size_t reverseId = segments.size();
        for (std::size_t candidate = offsets[neighbor]; candidate < offsets[neighbor + 1]; ++candidate)
        {
            const auto& reverse = segments[candidate];
            if (reverse.isHorizontal != face.isHorizontal || reverse.edgeLabel != face.edgeLabel) continue;
            require(reverseId == segments.size() && reverse.K == face.L && reverse.L == face.K,
                "interior correction requires one reverse incidence of the same actual face");
            reverseId = candidate;
        }
        require(reverseId < segments.size(), "interior correction has no reverse actual face");
        if (sid > reverseId) continue; // Check every physical pair once, without coalescing split faces.
        const auto& forward = corrections[sid];
        const auto& reverse = corrections[reverseId];
        require(reverse.constant == -forward.constant && reverse.coefficients.size() == forward.coefficients.size(),
            "interior correction pair must be exactly antisymmetric");
        for (std::size_t j = 0; j < forward.coefficients.size(); ++j)
            require(reverse.coefficients[j].cellIndex == forward.coefficients[j].cellIndex
                && reverse.coefficients[j].value == -forward.coefficients[j].value,
                "interior correction coefficients must be exactly antisymmetric");
    }
}

template<class Get>
void validatePreparedCorrections(
    const std::vector<PolarizationPoissonFVSolverCellCentered::FittedBoundary>& boundaries,
    const std::vector<std::vector<int>>& columns,const std::vector<std::size_t>& reverse,
    std::size_t count,const Get& get)
{
    using Type=PolarizationPoissonFVSolverCellCentered::FittedBoundaryType;
    const auto require=[](bool ok,const char* what){if(!ok)throw std::runtime_error(std::string("fitted fluxCorrections: ")+what);};
    require(count==columns.size(),"coefficient schema must cover every oriented segment");
    for(std::size_t sid=0;sid<count;++sid){const auto& correction=get(sid);
        require(std::isfinite(correction.constant),"non-finite correction constant");
        require(correction.coefficients.size()==columns[sid].size(),"coefficient count differs from owned integrated schema");
        for(std::size_t k=0;k<columns[sid].size();++k){const auto& c=correction.coefficients[k];
            require(c.cellIndex==columns[sid][k]&&std::isfinite(c.value),"coefficient index/value differs from owned active schema");}
        if(boundaries[sid].type==Type::HomogeneousNeumann){require(correction.constant==0,"homogeneous Neumann boundary correction constant must be exactly zero");
            for(const auto& c:correction.coefficients)require(c.value==0,"homogeneous Neumann boundary correction coefficients must each be exactly zero");}
    }
    for(std::size_t sid=0;sid<count;++sid)if(boundaries[sid].type==Type::Interior&&sid<reverse[sid]){
        const auto& a=get(sid);const auto& b=get(reverse[sid]);
        require(b.constant==-a.constant,"interior correction pair must be exactly antisymmetric");
        for(std::size_t k=0;k<a.coefficients.size();++k)
            require(b.coefficients[k].value==-a.coefficients[k].value,"interior correction coefficients must be exactly antisymmetric");
    }
}

} // namespace

bool PolarizationPoissonFVSolverCellCentered::solve(const QuadtreeMesh& mesh,
                                                   const FaceSegmentGraph& graph,
                                                   const PicDomainBC& bc,
                                                   const CellField<double>& rhs,
                                                   CellField<double>& phi,
                                                   const std::vector<double>* kappaBySegment,
                                                   const std::vector<double>* kappaByCellIndex,
                                                   SolveReport* report,
                                                   const Options& opt,
                                                   const Solvers::Poisson::CutCellMetrics* cutMetrics,
                                                   const PicEmbeddedBC* embeddedBC,
                                                   const PicEmbeddedBCMap* embeddedBCById,
                                                   const std::vector<SegmentFluxCorrection>* fluxCorrections) const
{
    // This option has a consumer only on the explicit fitted-boundary port.
    fittedFactorCache_.reset();
    if (opt.reuseFittedFactorization)
    {
        if (report)
        {
            *report = {};
            report->message = "reuseFittedFactorization requires solveFittedBoundary";
        }
        return false;
    }
    if (opt.axisymmetric && !fluxCorrections)
    {
        if (report)
        {
            *report = {};
            report->message = "axisymmetric mode requires the strict integrated fluxCorrections port";
        }
        return false;
    }
    if (!fluxCorrections)
        return solveImpl(mesh, graph, bc, rhs, phi, kappaBySegment, kappaByCellIndex,
                         report, opt, cutMetrics, embeddedBC, embeddedBCById, nullptr);

    const auto begin = Clock::now();
    // The corrected operator has no cache identity contract. Never consult old
    // pattern/value hashes (including environment-enabled freeze modes).
    clearDirectCache();
    try
    {
        PolarizationFluxCorrectionValidation::validate(mesh, graph, bc, rhs, phi,
            kappaBySegment, kappaByCellIndex, opt, cutMetrics, embeddedBC,
            embeddedBCById, *fluxCorrections);
        Options strictOptions = opt;
        strictOptions.reuseDirectPattern = false;
        strictOptions.freezeDirectMatrixUntilPatternChange = false;
        strictOptions.reuseIterativePattern = false;
        strictOptions.useTemporalExtrapolatedGuess = false;
        const bool success = solveImpl(mesh, graph, bc, rhs, phi, kappaBySegment,
            kappaByCellIndex, report, strictOptions, nullptr, nullptr, nullptr, fluxCorrections);
        if (!success) clearDirectCache();
        return success;
    }
    catch (const std::exception& error)
    {
        clearDirectCache();
        if (report)
        {
            *report = {};
            report->message = std::string("Integrated flux-correction solve failed: ") + error.what();
            report->timeTotal_s = secondsBetween(begin, Clock::now());
        }
        return false;
    }
}

bool PolarizationPoissonFVSolverCellCentered::solveFittedBoundary(
    const QuadtreeMesh& mesh, const FaceSegmentGraph& graph,
    const CellField<double>& rhs, CellField<double>& phi,
    const std::vector<double>& kappaBySegment,
    const std::vector<FittedBoundary>& boundaryBySegment,
    const Options& opt, SolveReport* report) const
{
    const auto begin = Clock::now();
    if (!opt.reuseFittedFactorization) clearDirectCache();
    try
    {
        PolarizationFittedBoundaryValidation::validate(mesh, graph, rhs, phi,
            kappaBySegment, boundaryBySegment, opt);
        // The shared solver publishes into this owned stage. Inactive values
        // are retained, and final publication needs no allocation or relabeling.
        CellField<double> nextPhi(phi);
        SolveReport nextReport;
        Options strictOptions = opt;
        strictOptions.reuseDirectPattern = false;
        strictOptions.freezeDirectMatrixUntilPatternChange = false;
        strictOptions.reuseIterativePattern = false;
        strictOptions.useTemporalExtrapolatedGuess = false;
        const PicDomainBC unusedDomainBC{};
        const bool success = solveImpl(mesh, graph, unusedDomainBC, rhs, nextPhi,
            &kappaBySegment, nullptr, &nextReport, strictOptions, nullptr, nullptr,
            nullptr, nullptr, &boundaryBySegment);
        if (!success) clearDirectCache();
        // Report construction/copy also precedes the no-throw value commit.
        if (report) *report = std::move(nextReport);
        if (!success) return false;
        std::copy_n(nextPhi.value.get(), phi.size(), phi.value.get());
        return true;
    }
    catch (const std::exception& error)
    {
        clearDirectCache();
        if (report)
        {
            *report = {};
            report->message = std::string("Fitted-boundary solve failed: ") + error.what();
            report->timeTotal_s = secondsBetween(begin, Clock::now());
        }
        return false;
    }
}

bool PolarizationPoissonFVSolverCellCentered::solveFittedBoundary(
    const QuadtreeMesh& mesh, const FaceSegmentGraph& graph,
    const CellField<double>& rhs, CellField<double>& phi,
    const std::vector<double>& kappaBySegment,
    const std::vector<FittedBoundary>& boundaryBySegment,
    const std::vector<SegmentFluxCorrection>& fluxCorrections,
    const Options& opt, SolveReport* report) const
{
    const auto begin = Clock::now();
    if (!opt.reuseFittedFactorization) clearDirectCache();
    try
    {
        PolarizationFittedBoundaryValidation::validate(mesh, graph, rhs, phi,
            kappaBySegment, boundaryBySegment, opt);
        validateFittedFluxCorrections(mesh, graph, boundaryBySegment, fluxCorrections);
        // Keep the old fitted entry untouched. Both explicit inputs use the
        // existing combined assembly and publish only this owned staged field.
        CellField<double> nextPhi(phi);
        SolveReport nextReport;
        Options strictOptions = opt;
        strictOptions.reuseDirectPattern = false;
        strictOptions.freezeDirectMatrixUntilPatternChange = false;
        strictOptions.reuseIterativePattern = false;
        strictOptions.useTemporalExtrapolatedGuess = false;
        const PicDomainBC unusedDomainBC{};
        const bool success = solveImpl(mesh, graph, unusedDomainBC, rhs, nextPhi,
            &kappaBySegment, nullptr, &nextReport, strictOptions, nullptr, nullptr,
            nullptr, &fluxCorrections, &boundaryBySegment);
        if (!success) clearDirectCache();
        // No observable phi write precedes report construction, complete
        // corrected-operator residual verification, or a possible exception.
        if (report) *report = std::move(nextReport);
        if (!success) return false;
        std::copy_n(nextPhi.value.get(), phi.size(), phi.value.get());
        return true;
    }
    catch (const std::exception& error)
    {
        clearDirectCache();
        if (report)
        {
            *report = {};
            report->message = std::string("Fitted-boundary integrated flux-correction solve failed: ") + error.what();
            report->timeTotal_s = secondsBetween(begin, Clock::now());
        }
        return false;
    }
}

PolarizationPoissonFVSolverCellCentered::PreparedIntegratedStructure::PreparedIntegratedStructure(
    std::shared_ptr<const Data> data):data_(std::move(data)){}
const FaceSegmentGraph& PolarizationPoissonFVSolverCellCentered::PreparedIntegratedStructure::graph() const noexcept
{return *data_->graph;}
std::shared_ptr<const PolarizationPoissonFVSolverCellCentered::PreparedIntegratedStructure>
PolarizationPoissonFVSolverCellCentered::prepareIntegratedStructure(const QuadtreeMesh& mesh,
    std::unique_ptr<FaceSegmentGraph> graph,const std::vector<FittedBoundaryType>& types,
    const std::vector<std::vector<int>>& correctionColumns,const std::vector<std::vector<int>>& affineColumns,
    bool axisymmetric)
{
    const auto require=[](bool ok,const char* what){if(!ok)throw std::runtime_error(what);};
    require(bool(graph),"prepared integrated structure requires an actual owned graph");
    std::vector<FittedBoundary> boundaries(types.size());
    for(std::size_t sid=0;sid<types.size();++sid)boundaries[sid].type=types[sid];
    Options geometryOptions;geometryOptions.axisymmetric=axisymmetric;CellField<double> unused;
    PolarizationFittedBoundaryValidation::validate<true>(mesh,*graph,unused,unused,{},boundaries,geometryOptions);
    auto data=std::make_shared<PreparedIntegratedStructure::Data>(mesh,std::move(graph));
    data->axisymmetric=axisymmetric;data->types=types;
    data->correctionColumns=correctionColumns;data->affineColumns=affineColumns;
    const auto& segments=data->graph->segments();const auto& offsets=data->graph->cellOffsets();
    const auto cells=static_cast<std::size_t>(mesh.cellIndexNum());
    data->rowOfCell.assign(cells,-1);
    for(std::size_t ci=0;ci<cells;++ci)if(offsets[ci]<offsets[ci+1]){
        data->rowOfCell[ci]=static_cast<int>(data->cellOfRow.size());data->cellOfRow.push_back(static_cast<int>(ci));}
    const auto n=data->cellOfRow.size();require(n>0,"No active cells in FaceSegmentGraph (N=0).");
    const auto validateColumns=[&](const auto& columns){
        require(columns.size()==segments.size(),"prepared integrated column schema dimensions differ");
        for(const auto& row:columns){int previous=-1;for(int ci:row){
            require(ci>previous&&static_cast<std::size_t>(ci)<cells&&data->rowOfCell[ci]>=0,
                "prepared integrated columns must be sorted unique active native cells");previous=ci;}}
    };
    validateColumns(correctionColumns);if(!affineColumns.empty())validateColumns(affineColumns);
    data->reverse.assign(segments.size(),segments.size());
    std::map<std::pair<bool,int>,std::size_t> firstFace;
    for(std::size_t sid=0;sid<segments.size();++sid){const auto& f=segments[sid];
        const auto inserted=firstFace.emplace(std::make_pair(f.isHorizontal,f.edgeLabel),sid);
        if(!inserted.second){const auto other=inserted.first->second;data->reverse[sid]=other;data->reverse[other]=sid;}}
    for(std::size_t sid=0;sid<segments.size();++sid)if(types[sid]==FittedBoundaryType::Interior){
        const auto other=data->reverse[sid];require(other<segments.size(),"prepared integrated interior has no reverse");
        require(correctionColumns[sid]==correctionColumns[other],"prepared integrated reverse correction columns differ");
        if(!affineColumns.empty())require(affineColumns[sid]==affineColumns[other],"prepared integrated reverse affine columns differ");}
    std::vector<Eigen::Triplet<double>> terms;data->termOffsets.resize(segments.size()+1);
    std::vector<std::vector<int>> neighbors(n);std::vector<unsigned char> anchored(n,0);
    for(std::size_t sid=0;sid<segments.size();++sid){data->termOffsets[sid]=terms.size();const auto& f=segments[sid];
        const auto i=data->rowOfCell[mesh.cells(f.K).index()];
        for(int ci:correctionColumns[sid])terms.emplace_back(i,data->rowOfCell[ci],0.);
        const double area=axisymmetric?CartesianMeasure::face(f.Af,f.faceCenter.x,true):f.Af;
        if(area==0||types[sid]==FittedBoundaryType::HomogeneousNeumann)continue;
        terms.emplace_back(i,i,0.);
        if(types[sid]==FittedBoundaryType::Dirichlet)anchored[i]=1;
        else {const auto j=data->rowOfCell[mesh.cells(f.L).index()];terms.emplace_back(i,j,0.);
            neighbors[i].push_back(j);neighbors[j].push_back(i);}
    }
    data->termOffsets.back()=terms.size();
    std::vector<unsigned char> visited(n,0);std::vector<int> queue;queue.reserve(n);
    for(std::size_t seed=0;seed<n;++seed)if(!visited[seed]){queue.clear();queue.push_back(static_cast<int>(seed));visited[seed]=1;bool hasAnchor=false;
        for(std::size_t at=0;at<queue.size();++at){const auto row=queue[at];hasAnchor=hasAnchor||anchored[row]!=0;
            for(int other:neighbors[row])if(!visited[other]){visited[other]=1;queue.push_back(other);}}
        require(hasAnchor,"fluxCorrections: unanchored component; compatibility correction/gauge pin is forbidden");}
    data->pattern.resize(static_cast<int>(n),static_cast<int>(n));
    data->pattern.setFromTriplets(terms.begin(),terms.end());data->pattern.makeCompressed();
    data->slots.reserve(terms.size());data->first.reserve(terms.size());
    std::vector<unsigned char> seen(static_cast<std::size_t>(data->pattern.nonZeros()),0);
    for(const auto& term:terms){const auto first=data->pattern.outerIndexPtr()[term.col()],last=data->pattern.outerIndexPtr()[term.col()+1];
        const auto at=std::lower_bound(data->pattern.innerIndexPtr()+first,data->pattern.innerIndexPtr()+last,term.row());
        require(at!=data->pattern.innerIndexPtr()+last&&*at==term.row(),"prepared integrated CSC slot is missing");
        const auto slot=static_cast<std::size_t>(at-data->pattern.innerIndexPtr());data->slots.push_back(slot);
        data->first.push_back(seen[slot]?0u:1u);seen[slot]=1;
    }
    // Stable CSC-to-row routing can be borrowed by each same-pattern residual.
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
    IntegratedCurrentEnvironment optionalRouteEnvironment;
#endif
    const auto routeNnz=static_cast<std::size_t>(data->pattern.nonZeros());
    const auto routeLimit=ExpandedResidualRoute::limit;
    if(n<std::numeric_limits<std::uint32_t>::max()&&routeNnz<=std::numeric_limits<std::uint32_t>::max()
        &&n+1<=routeLimit/(2*sizeof(std::uint32_t))
        &&routeNnz<=(routeLimit-2*(n+1)*sizeof(std::uint32_t))/sizeof(std::array<std::uint32_t,2>))try{
        data->residualOffsets.assign(n+1,0);data->residualEntries.resize(static_cast<std::size_t>(data->pattern.nonZeros()));
        for(Eigen::Index at=0;at<data->pattern.nonZeros();++at)++data->residualOffsets[static_cast<std::size_t>(data->pattern.innerIndexPtr()[at])+1];
        std::partial_sum(data->residualOffsets.begin(),data->residualOffsets.end(),data->residualOffsets.begin());auto cursor=data->residualOffsets;
        for(int column=0;column<data->pattern.outerSize();++column)
            for(auto at=data->pattern.outerIndexPtr()[column];at<data->pattern.outerIndexPtr()[column+1];++at)
                data->residualEntries[cursor[static_cast<std::size_t>(data->pattern.innerIndexPtr()[at])]++]={static_cast<std::uint32_t>(column),static_cast<std::uint32_t>(at)};
        if((data->residualOffsets.capacity()+cursor.capacity())*sizeof(std::uint32_t)
            +data->residualEntries.capacity()*sizeof(data->residualEntries[0])>routeLimit)
            throw std::length_error("optional prepared residual actual capacity cap");
    }catch(...){std::vector<std::uint32_t>().swap(data->residualOffsets);
        std::vector<std::array<std::uint32_t,2>>().swap(data->residualEntries);}
    return std::shared_ptr<const PreparedIntegratedStructure>(new PreparedIntegratedStructure(std::move(data)));
}

bool PolarizationPoissonFVSolverCellCentered::solveFittedIntegratedBoundary(
    const QuadtreeMesh& mesh, const FaceSegmentGraph& graph,
    const CellField<double>& integratedRhs, CellField<double>& phi,
    const std::vector<double>& kappaBySegment,
    const std::vector<FittedBoundary>& boundaryBySegment,
    const std::vector<SegmentFluxCorrection>& fluxCorrections,
    const Options& opt, SolveReport* report,
    const std::vector<IntegratedFluxFactor>* integratedFlux,
    std::vector<std::array<double,3>>* potentialTerms) const
{
    return solveFittedIntegratedBoundaryImpl(mesh,graph,integratedRhs,phi,kappaBySegment,boundaryBySegment,
        fluxCorrections,opt,report,integratedFlux,potentialTerms,nullptr);
}
bool PolarizationPoissonFVSolverCellCentered::solvePreparedIntegratedBoundary(
    const QuadtreeMesh& mesh,const PreparedIntegratedStructure& structure,
    const CellField<double>& integratedRhs,CellField<double>& phi,
    const std::vector<double>& kappaBySegment,const std::vector<FittedBoundary>& boundaryBySegment,
    const std::vector<SegmentFluxCorrection>& fluxCorrections,const Options& opt,SolveReport* report,
    const std::vector<IntegratedFluxFactor>* integratedFlux,std::vector<std::array<double,3>>* potentialTerms) const
{
    return solveFittedIntegratedBoundaryImpl(mesh,structure.graph(),integratedRhs,phi,kappaBySegment,boundaryBySegment,
        fluxCorrections,opt,report,integratedFlux,potentialTerms,&structure);
}
bool PolarizationPoissonFVSolverCellCentered::solveFittedIntegratedBoundaryImpl(
    const QuadtreeMesh& mesh, const FaceSegmentGraph& graph,
    const CellField<double>& integratedRhs, CellField<double>& phi,
    const std::vector<double>& kappaBySegment,
    const std::vector<FittedBoundary>& boundaryBySegment,
    const std::vector<SegmentFluxCorrection>& fluxCorrections,
    const Options& opt, SolveReport* report,
    const std::vector<IntegratedFluxFactor>* integratedFlux,
    std::vector<std::array<double,3>>* potentialTerms, const PreparedIntegratedStructure* prepared) const
{
    const auto begin = Clock::now();
    // This port owns no volumetric RHS plan or retained factorization.
    clearDirectCache();
    try
    {
        if (opt.reuseFittedFactorization)
            throw std::runtime_error("integrated RHS fitted port requires reuseFittedFactorization=false");
        if ((integratedFlux != nullptr) != (potentialTerms != nullptr))
            throw std::runtime_error("integrated face factors require the paired native potential expansion output");
        if (integratedFlux && !opt.requireSuperLU && !opt.useGMRES)
            throw std::runtime_error("integrated face factors require an explicit strict SuperLU or GMRES contract");
        if(prepared){const auto& data=*prepared->data_;
            // AMR commits move validated storage to another mesh wrapper.
            // The exact geometry/connectivity identity owns this value-only plan.
            data.identity.validate(mesh);
            if(data.axisymmetric!=opt.axisymmetric)throw std::runtime_error("prepared integrated structure has a different measure");
            PolarizationFittedBoundaryValidation::validatePreparedValues(mesh,graph,integratedRhs,phi,
                kappaBySegment,boundaryBySegment,opt,data.types,data.reverse);
            validatePreparedCorrections(boundaryBySegment,data.correctionColumns,data.reverse,fluxCorrections.size(),
                [&](std::size_t sid)->const SegmentFluxCorrection&{return fluxCorrections[sid];});
        }else{
            PolarizationFittedBoundaryValidation::validate(mesh, graph, integratedRhs, phi,
                kappaBySegment, boundaryBySegment, opt);
            validateFittedFluxCorrections(mesh, graph, boundaryBySegment, fluxCorrections);
        }
        if (integratedFlux)
        {
            const auto& segments = graph.segments();
            if (integratedFlux->size() != segments.size())
                throw std::runtime_error("integrated face factor count differs from actual oriented graph");
            std::vector<SegmentFluxCorrection> affine;if(!prepared)affine.reserve(segments.size());
            for (std::size_t sid = 0; sid < segments.size(); ++sid)
            {
                const auto& factor = (*integratedFlux)[sid];
                if (!(std::isfinite(factor.primary) && factor.primary > 0 &&
                    std::isfinite(factor.secondary) && factor.secondary >= 0 &&
                    std::isfinite(factor.weight) && factor.weight >= 0 &&
                    (factor.secondaryTransfer == IntegratedFluxFactor::SecondaryTransfer::Exact ||
                     factor.secondaryTransfer == IntegratedFluxFactor::SecondaryTransfer::MinimumBinary64QuantumNearestEven)))
                    throw std::runtime_error("invalid integrated face action factors");
                const double combined = static_cast<double>(static_cast<long double>(factor.primary) +
                    static_cast<long double>(factor.weight) * factor.secondary);
                if (combined != kappaBySegment[sid])
                    throw std::runtime_error("integrated face factors differ from actual implicit coefficient");
                const auto& segment = segments[sid];
                const auto& boundary = boundaryBySegment[sid];
                if (boundary.type == FittedBoundaryType::Interior && factor.affine.constant != 0)
                    throw std::runtime_error("interior integrated face action has a foreign constant");
                if (boundary.type == FittedBoundaryType::Dirichlet)
                {
                    const auto center = mesh.cells(segment.K).centerCoord();
                    const long double distance = segment.isHorizontal
                        ? static_cast<long double>(segment.n.y) * (static_cast<long double>(segment.faceCenter.y) - center.y)
                        : static_cast<long double>(segment.n.x) * (static_cast<long double>(segment.faceCenter.x) - center.x);
                    const double area = opt.axisymmetric
                        ? CartesianMeasure::face(segment.Af, segment.faceCenter.x, true) : segment.Af;
                    const double expected = static_cast<double>(-(static_cast<long double>(area) / distance) * boundary.value);
                    if (factor.affine.constant != expected)
                        throw std::runtime_error("integrated face action uses a different fitted boundary constant");
                }
                if(!prepared)affine.push_back(factor.affine);
            }
            // Includes actual active-map, zero-Neumann and exact reverse
            // coefficient/constant checks; no second graph is constructed.
            if(prepared)validatePreparedCorrections(boundaryBySegment,prepared->data_->affineColumns,prepared->data_->reverse,integratedFlux->size(),
                [&](std::size_t sid)->const SegmentFluxCorrection&{return (*integratedFlux)[sid].affine;});
            else validateFittedFluxCorrections(mesh, graph, boundaryBySegment, affine);
            const auto& offsets = graph.cellOffsets();
            for (std::size_t sid = 0; sid < segments.size(); ++sid)
            {
                if (boundaryBySegment[sid].type != FittedBoundaryType::Interior) continue;
                const auto& segment = segments[sid]; const auto& forward = (*integratedFlux)[sid];
                const int neighbor = mesh.cells(segment.L).index();
                const auto first=prepared?prepared->data_->reverse[sid]:offsets[neighbor];
                const auto last=prepared?first+1:offsets[neighbor+1];
                for (std::size_t reverse=first;reverse<last;++reverse)
                    if (segments[reverse].isHorizontal == segment.isHorizontal && segments[reverse].edgeLabel == segment.edgeLabel)
                    {
                        const auto& other = (*integratedFlux)[reverse];
                        if (other.primary != forward.primary || other.secondary != forward.secondary || other.weight != forward.weight ||
                            other.secondaryTransfer != forward.secondaryTransfer)
                            throw std::runtime_error("reverse integrated face factors differ on the same physical face");
                    }
            }
        }
        CellField<double> nextPhi(phi);
        std::vector<std::array<double,3>> nextPotentialTerms;
        SolveReport nextReport;
        Options strictOptions = opt;
        strictOptions.reuseDirectPattern = false;
        strictOptions.freezeDirectMatrixUntilPatternChange = false;
        strictOptions.reuseIterativePattern = false;
        strictOptions.useTemporalExtrapolatedGuess = false;
        const PicDomainBC unusedDomainBC{};
        const bool success = solveImpl(mesh, graph, unusedDomainBC, integratedRhs, nextPhi,
            &kappaBySegment, nullptr, &nextReport, strictOptions, nullptr, nullptr,
            nullptr, &fluxCorrections, &boundaryBySegment, true, integratedFlux,
            potentialTerms ? &nextPotentialTerms : nullptr, prepared);
        if (!success) clearDirectCache();
        if (success)
            nextReport.message = "Solved face-fitted active-leaf Poisson with integrated cell RHS and flux corrections; full-operator residual verified.";
        // All validation, residual gates and report construction precede the
        // same no-allocation publication used by the corrected fitted entry.
        if (report) *report = std::move(nextReport);
        if (!success) return false;
        if (potentialTerms) potentialTerms->swap(nextPotentialTerms);
        std::copy_n(nextPhi.value.get(), phi.size(), phi.value.get());
        return true;
    }
    catch (const std::exception& error)
    {
        clearDirectCache();
        if (report)
        {
            *report = {};
            report->message = std::string("Fitted-boundary integrated RHS solve failed: ") + error.what();
            report->timeTotal_s = secondsBetween(begin, Clock::now());
        }
        return false;
    }
}

bool PolarizationPoissonFVSolverCellCentered::solvePreparedFittedRhs(
    const CellField<double>& rhs, CellField<double>& phi,
    const std::vector<double>* constants, double relTol, SolveReport* report) const
{
    const auto begin = Clock::now();
    try
    {
        if (!fittedFactorCache_ || !fittedFactorCache_->rhsPlan)
            throw std::runtime_error("no successfully prepared fitted operator");
        if (!std::isfinite(relTol) || relTol <= 0.0 || relTol >= 1.0)
            throw std::runtime_error("relTol must be finite and in (0,1)");
        auto& prepared = *fittedFactorCache_;
        const auto b = prepared.rhsPlan->assemble(rhs, phi, constants);
        SolveReport next;
        next.N = static_cast<int>(prepared.matrix.rows());
        next.nnz = static_cast<long long>(prepared.matrix.nonZeros());
        next.fittedFactorizationReused = next.fittedRhsOnlySolve = true;
        next.timeAssemble_s = secondsBetween(begin, Clock::now());
        Eigen::VectorXd x;
#ifdef EIGEN_SUPERLU_SUPPORT
        const auto solveBegin = Clock::now();
        x = prepared.solver.solve(b);
        next.timeSolve_s = secondsBetween(solveBegin, Clock::now());
        next.usedSuperLU = prepared.solver.info() == Eigen::Success;
#endif
        std::string reason = "prepared fitted SuperLU backsolve failed; no fallback is allowed";
        next.success = next.usedSuperLU && fullResidualAccepted(prepared.matrix, x, b,
            next.usedSuperLU, relTol, next.finalResidual, reason);
        next.timeTotal_s = secondsBetween(begin, Clock::now());
        next.message = next.success
            ? "Solved owned face-fitted operator with updated RHS; full-operator residual verified."
            : reason;
        const bool accepted = next.success;
        // Report construction and all fallible work precede phi publication.
        if (report) *report = std::move(next);
        if (!accepted)
        {
            clearDirectCache();
            return false;
        }
        for (std::size_t row = 0; row < prepared.cellIndices.size(); ++row)
            phi.value[prepared.cellIndices[row]] = x[row];
        return true;
    }
    catch (const std::exception& error)
    {
        clearDirectCache();
        if (report)
        {
            *report = {};
            report->message = std::string("Prepared fitted RHS solve failed: ") + error.what();
            report->timeTotal_s = secondsBetween(begin, Clock::now());
        }
        return false;
    }
    catch (...)
    {
        clearDirectCache();
        throw;
    }
}

bool PolarizationPoissonFVSolverCellCentered::solveImpl(const QuadtreeMesh& mesh,
                                                   const FaceSegmentGraph& graph,
                                                   const PicDomainBC& bc,
                                                   const CellField<double>& rhs,
                                                   CellField<double>& phi,
                                                   const std::vector<double>* kappaBySegment,
                                                   const std::vector<double>* kappaByCellIndex,
                                                   SolveReport* report,
                                                   const Options& opt,
                                                   const Solvers::Poisson::CutCellMetrics* cutMetrics,
                                                   const PicEmbeddedBC* embeddedBC,
                                                   const PicEmbeddedBCMap* embeddedBCById,
                                                   const std::vector<SegmentFluxCorrection>* fluxCorrections,
                                                   const std::vector<FittedBoundary>* fittedBoundaries,
                                                   bool integratedRhs,
                                                   const std::vector<IntegratedFluxFactor>* integratedFlux,
                                                   std::vector<std::array<double,3>>* potentialTerms, const PreparedIntegratedStructure* prepared) const
{
    const auto tBegin = Clock::now();
    const auto* structure=prepared?prepared->data_.get():nullptr;
    const bool correctedOperator = fluxCorrections != nullptr;
    const bool fittedOperator = fittedBoundaries != nullptr;
    const bool checkedOperator = correctedOperator || fittedOperator;

    // Validate the strict backend contract before changing phi or solver state.
    auto requiredBackendFailure = [&](const std::string& message) {
        if (report)
        {
            *report = {};
            report->message = message;
        }
        return false;
    };
    if (opt.useGMRES && (!integratedRhs || !fittedOperator || opt.requireSuperLU || opt.useSuperLU ||
        opt.useCG || opt.useIC || opt.reuseFittedFactorization || opt.maxIter <= 0 ||
        opt.gmresRestart <= 0 || opt.gmresBlockRows <= 0))
        return requiredBackendFailure("explicit GMRES requires the integrated fitted port, positive iteration/block limits and no alternative backend or fallback");
    if (opt.requireSuperLU && (!opt.useSuperLU || opt.useGMRES || opt.useCG || opt.useIC))
        return requiredBackendFailure("requireSuperLU requires useSuperLU=true, useCG=false and useIC=false; no alternative backend is allowed");
#ifndef EIGEN_SUPERLU_SUPPORT
    if (opt.requireSuperLU)
        return requiredBackendFailure("requireSuperLU requested but SuperLU support is unavailable in this build; SparseLU fallback is disabled");
#endif

    auto secondsSince = [&](const Clock::time_point& t0) -> double {
        return secondsBetween(t0, Clock::now());
    };

    if (!checkedOperator && (phi.size() != mesh.cellIndexNum() || !phi.value))
        phi.initialize(mesh, 0.0);

    // Active unknown set: cell indices with at least one face segment.
    const int nCellIndex = mesh.cellIndexNum();
    std::vector<int> freshRowOfCellIndex(structure?0:static_cast<size_t>(nCellIndex), -1);
    std::vector<int> freshCellIndexOfRow;
    if(!structure)freshCellIndexOfRow.reserve(static_cast<size_t>(nCellIndex));

    if(!structure)
    for (int ci = 0; ci < nCellIndex; ++ci)
    {
        const auto range = graph.segmentRangeForCellIndex(ci);
        if (range.first < range.second)
        {
            freshRowOfCellIndex[static_cast<size_t>(ci)] = static_cast<int>(freshCellIndexOfRow.size());
            freshCellIndexOfRow.push_back(ci);
        }
    }

    const auto& rowOfCellIndex=structure?structure->rowOfCell:freshRowOfCellIndex;
    const auto& cellIndexOfRow=structure?structure->cellOfRow:freshCellIndexOfRow;
    const int N = static_cast<int>(cellIndexOfRow.size());
    if (N == 0)
    {
        if (report)
        {
            *report = {};
            report->success = false;
            report->message = "No active cells in FaceSegmentGraph (N=0).";
            report->timeTotal_s = secondsSince(tBegin);
            report->timeAssemble_s = report->timeTotal_s;
        }
        return false;
    }

    const auto& segs = graph.segments();
    const bool useKappaBySegment =
        (kappaBySegment != nullptr) && (kappaBySegment->size() == segs.size());
    const bool useKappaByCell =
        (kappaByCellIndex != nullptr) && (kappaByCellIndex->size() == static_cast<size_t>(mesh.cellIndexNum()));

    const bool useCutMetrics =
        (cutMetrics != nullptr) &&
        (cutMetrics->betaBySegment().size() == segs.size()) &&
        (cutMetrics->cellMetricsByIndex().size() == static_cast<size_t>(mesh.cellIndexNum()));

    const bool useEmbeddedBCMap =
        useCutMetrics && (embeddedBCById != nullptr) && (!embeddedBCById->empty());
    const bool useEmbeddedBCFallback =
        useCutMetrics && (embeddedBC != nullptr) && embeddedBC->enabled;
    const bool useEmbeddedBC = useEmbeddedBCMap || useEmbeddedBCFallback;

    // Per-row volume weight: legacy alpha*area, or revolved full-cell volume
    // for the validated strict axisymmetric port (which forbids cut metrics).
    std::vector<double> volOfRow(static_cast<size_t>(N), 0.0);
    for (int i = 0; i < N; ++i)
    {
        const int cellIndex = cellIndexOfRow[i];
        const double area = opt.axisymmetric
            ? CartesianMeasure::cell(mesh.cells(mesh.cellIndexToLabel[cellIndex]).area(),
                mesh.cells(mesh.cellIndexToLabel[cellIndex]).centerCoord().x, true)
            : safeCellArea(mesh, cellIndex);

        double alpha = 1.0;
        if (useCutMetrics)
        {
            alpha = cutMetrics->cellMetricsByIndex()[static_cast<size_t>(cellIndex)].alpha;
            if (!std::isfinite(alpha))
                alpha = 1.0;
            if (alpha < 0.0)
                alpha = 0.0;
            else if (alpha > 1.0)
                alpha = 1.0;
        }
        volOfRow[static_cast<size_t>(i)] = area * alpha;
    }

    std::unique_ptr<PolarizationFittedBoundaryDetail::RhsPlan> preparedRhs;
    if (fittedOperator && opt.reuseFittedFactorization)
        preparedRhs = std::make_unique<PolarizationFittedBoundaryDetail::RhsPlan>(mesh, graph,
            *fittedBoundaries, cellIndexOfRow, rowOfCellIndex, volOfRow, correctedOperator);

    // b_i = volumetric rhs * physical cell measure.
    Eigen::VectorXd bVec(N);
    bVec.setZero();
    const bool rhsOk = (rhs.size() == mesh.cellIndexNum()) && (rhs.value != nullptr);
    for (int i = 0; i < N; ++i)
    {
        const int cellIndex = cellIndexOfRow[i];
        const double f = rhsOk ? rhs.value[cellIndex] : 0.0;
        if (integratedRhs)
            bVec[i] = f;
        else
            bVec[i] = f * volOfRow[static_cast<size_t>(i)];
    }

    const auto tAsm0 = Clock::now();

    Eigen::SparseMatrix<double> A;
    if(!structure){
    std::vector<Eigen::Triplet<double>> trips;
    trips.reserve(segs.size() * 2u + static_cast<std::size_t>(N) + 32u);

    std::vector<double> diag(static_cast<size_t>(N), 0.0);
    std::vector<std::vector<int>> neighbors(static_cast<size_t>(N));
    std::vector<unsigned char> anchoredRow(static_cast<size_t>(N), 0u);

    auto kappaForSegment = [&](size_t sid) -> double {
        if (!useKappaBySegment)
            return 1.0;
        return clampKappa((*kappaBySegment)[sid]);
    };

    auto kappaForCellIndex = [&](int cellIndex) -> double {
        if (!useKappaByCell || cellIndex < 0 || cellIndex >= mesh.cellIndexNum())
            return 1.0;
        return clampKappa((*kappaByCellIndex)[static_cast<size_t>(cellIndex)]);
    };

    bool warnedRobin = false;
    bool warnedPeriodicMismatch = false;

    for (size_t sid = 0; sid < segs.size(); ++sid)
    {
        const auto& seg = segs[sid];

        const int Kindex = mesh.cells(seg.K).index();
        if (Kindex < 0)
            continue;
        const int i = rowOfCellIndex[static_cast<size_t>(Kindex)];
        if (i < 0)
            continue;

        if (correctedOperator)
        {
            // Corrections already represent integrated physical flux. Their
            // RZ metric belongs to the caller and must not be applied twice.
            const auto& correction = (*fluxCorrections)[sid];
            for (const auto& coefficient : correction.coefficients)
            {
                const int j = rowOfCellIndex[static_cast<std::size_t>(coefficient.cellIndex)];
                trips.emplace_back(i, j, coefficient.value);
                if (j == i) diag[static_cast<std::size_t>(i)] += coefficient.value;
            }
            bVec[i] -= correction.constant;
        }

        double beta = 1.0;
        if (useCutMetrics)
        {
            beta = cutMetrics->betaBySegment()[sid];
            if (!std::isfinite(beta))
                beta = 1.0;
            if (beta < 0.0)
                beta = 0.0;
            else if (beta > 1.0)
                beta = 1.0;
        }

        const double AfOpen = opt.axisymmetric
            ? CartesianMeasure::face(seg.Af, seg.faceCenter.x, true)
            : beta * seg.Af;
        if (!(seg.df > 0.0) || !(AfOpen >= 0.0))
            continue;
        if (AfOpen == 0.0)
            continue;

        const double kappa_f = kappaForSegment(sid);

        if (fittedOperator)
        {
            const auto& boundary = (*fittedBoundaries)[sid];
            if (boundary.type == FittedBoundaryType::HomogeneousNeumann)
                continue; // Prescribed integrated flux is exactly zero.
            if (boundary.type == FittedBoundaryType::Dirichlet)
            {
                // An inactive neighbor is not an unknown at its center. The
                // supplied potential is on this actual physical face instead.
                const double distance = (seg.faceCenter - mesh.cells(seg.K).centerCoord()).dotProduct(seg.n);
                const double c = kappa_f * AfOpen / distance;
                if (preparedRhs) preparedRhs->recordDirichlet(sid, c, boundary.value);
                trips.emplace_back(i, i, c);
                diag[static_cast<std::size_t>(i)] += c;
                bVec[i] += c * boundary.value;
                anchoredRow[static_cast<std::size_t>(i)] = 1u;
                continue;
            }
        }

        if (seg.boundary == Boundary::non)
        {
            if (seg.L == -1)
                continue;

            const int Lindex = mesh.cells(seg.L).index();
            if (Lindex < 0)
                continue;

            const int j = rowOfCellIndex[static_cast<size_t>(Lindex)];
            if (j < 0)
                continue;

            const double c = kappa_f * AfOpen / seg.df;
            trips.emplace_back(i, i, c);
            trips.emplace_back(i, j, -c);
            diag[static_cast<size_t>(i)] += c;
            neighbors[static_cast<size_t>(i)].push_back(j);
            neighbors[static_cast<size_t>(j)].push_back(i);
        }
        else
        {
            const PicBoundarySide& s = bcSide(bc, seg.boundary);
            PicBCType t = s.type;

            if (t == PicBCType::Robin)
            {
                if (!warnedRobin)
                {
                    warnedRobin = true;
                    if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
                    {
                        std::cout << "[PolarizationPoissonFVSolverCellCentered] Warning: Robin BC encountered on domain boundary; "
                                  << "treating as Dirichlet(value).\n";
                    }
                }
                t = PicBCType::Dirichlet;
            }

            if (t == PicBCType::Dirichlet)
            {
                const double c = kappa_f * AfOpen / seg.df;
                trips.emplace_back(i, i, c);
                diag[static_cast<size_t>(i)] += c;
                bVec[i] += c * s.value;
                anchoredRow[static_cast<size_t>(i)] = 1u;
            }
            else if (t == PicBCType::Neumann)
            {
                // Interpret Neumann value as dphi/dn = g.
                // For -div(kappa grad phi)=rhs, boundary flux contribution scales by kappa.
                bVec[i] += (s.value * kappa_f) * AfOpen;
            }
            else if (t == PicBCType::Periodic)
            {
                const Boundary mateB = PeriodicFaceUtils::mateBoundary(seg.boundary);
                if (!warnedPeriodicMismatch)
                {
                    if (mateB == Boundary::non || bcSide(bc, mateB).type != PicBCType::Periodic)
                    {
                        warnedPeriodicMismatch = true;
                        if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
                        {
                            std::cout << "[PolarizationPoissonFVSolverCellCentered] Warning: Periodic BC on '"
                                      << boundaryToString(seg.boundary)
                                      << "' but mate side is not Periodic. Assembly will still wrap.\n";
                        }
                    }
                }

                const Vector2D shift = PeriodicFaceUtils::shiftToMate(mesh, seg.boundary);
                int Llabel = -1;
                if (!PeriodicFaceUtils::findMateCellLabel(mesh, seg, shift, Llabel))
                    continue;

                const int Lindex = mesh.cells(Llabel).index();
                if (Lindex < 0)
                    continue;

                const int j = rowOfCellIndex[static_cast<size_t>(Lindex)];
                if (j < 0)
                    continue;

                const double df = PeriodicFaceUtils::centerDistance(mesh, seg, Llabel, shift);
                if (!(df > 0.0))
                    continue;

                const double c = kappa_f * AfOpen / df;
                trips.emplace_back(i, i, c);
                trips.emplace_back(i, j, -c);
                diag[static_cast<size_t>(i)] += c;
                neighbors[static_cast<size_t>(i)].push_back(j);
                neighbors[static_cast<size_t>(j)].push_back(i);
            }
        }
    }

    // Embedded-boundary closure on cut cells.
    if (useEmbeddedBC)
    {
        for (int i = 0; i < N; ++i)
        {
            const int cellIndex = cellIndexOfRow[i];
            const auto& cm = cutMetrics->cellMetricsByIndex()[static_cast<size_t>(cellIndex)];
            if (!cm.isCut || !(cm.A_gamma > 0.0))
                continue;
            if (!(volOfRow[static_cast<size_t>(i)] > 0.0))
                continue;

            const PicEmbeddedBC* eb = nullptr;
            if (useEmbeddedBCMap && cm.gammaBoundaryId > 0)
            {
                auto it = embeddedBCById->find(cm.gammaBoundaryId);
                if (it != embeddedBCById->end() && it->second.enabled)
                    eb = &it->second;
            }
            if (!eb && useEmbeddedBCFallback)
                eb = embeddedBC;

            if (!eb)
            {
                if (opt.embeddedBCMissingPolicy == Options::EmbeddedBCMissingPolicy::Error)
                {
                    if (report)
                    {
                        *report = {};
                        report->success = false;
                        report->message = "Missing embedded BC for boundaryId " + std::to_string(cm.gammaBoundaryId) +
                                          " (no map entry and no fallback).";
                        report->timeTotal_s = secondsSince(tBegin);
                        report->timeAssemble_s = report->timeTotal_s;
                        report->N = N;
                    }
                    return false;
                }
                continue;
            }

            const double kappa_g = kappaForCellIndex(cellIndex);
            const double a = eb->a;
            const double b = eb->b;
            const double g = eb->g;
            const PicBCType t = eb->type;
            const double A_gamma = cm.A_gamma;

            if (t == PicBCType::Robin)
            {
                if (std::abs(b) > 0.0)
                {
                    const double addDiag = kappa_g * (a / b) * A_gamma;
                    const double addRhs = kappa_g * (g / b) * A_gamma;
                    if (!std::isfinite(addDiag) || !std::isfinite(addRhs))
                        continue;
                    if (addDiag != 0.0)
                    {
                        trips.emplace_back(i, i, addDiag);
                        diag[static_cast<size_t>(i)] += addDiag;
                        if (a != 0.0)
                        {
                            anchoredRow[static_cast<size_t>(i)] = 1u;
                        }
                    }
                    bVec[i] += addRhs;
                }
                else
                {
                    if (std::abs(a) <= 0.0)
                        continue;
                    const double phiD = g / a;
                    const double d = stabilizedEmbeddedDistance(cm.d_gamma,
                                                                volOfRow[static_cast<size_t>(i)]);

                    const double c = kappa_g * (A_gamma / d);
                    if (!std::isfinite(c) || c <= 0.0)
                        continue;

                    trips.emplace_back(i, i, c);
                    diag[static_cast<size_t>(i)] += c;
                    bVec[i] += c * phiD;
                    anchoredRow[static_cast<size_t>(i)] = 1u;
                }
            }
            else if (t == PicBCType::Dirichlet)
            {
                const double phiD = g;
                const double d = stabilizedEmbeddedDistance(cm.d_gamma,
                                                            volOfRow[static_cast<size_t>(i)]);

                const double c = kappa_g * (A_gamma / d);
                if (!std::isfinite(c) || c <= 0.0)
                    continue;

                trips.emplace_back(i, i, c);
                diag[static_cast<size_t>(i)] += c;
                bVec[i] += c * phiD;
                anchoredRow[static_cast<size_t>(i)] = 1u;
            }
            else if (t == PicBCType::Neumann)
            {
                const double addRhs = (g * kappa_g) * A_gamma;
                if (!std::isfinite(addRhs))
                    continue;
                bVec[i] += addRhs;
            }
        }
    }

    // Stabilize any floating connected component. This covers fully-Neumann
    // problems and mixed-BC cases where moving EB temporarily disconnects a
    // plasma island from all Dirichlet/Robin anchors.
    std::vector<unsigned char> visited(static_cast<size_t>(N), 0u);
    std::vector<int> componentRows;
    componentRows.reserve(static_cast<size_t>(N));
    size_t floatingComponentCount = 0u;

    for (int seed = 0; seed < N; ++seed)
    {
        if (visited[static_cast<size_t>(seed)] != 0u)
            continue;

        componentRows.clear();
        bool hasAnchor = false;
        std::queue<int> q;
        q.push(seed);
        visited[static_cast<size_t>(seed)] = 1u;

        while (!q.empty())
        {
            const int row = q.front();
            q.pop();
            componentRows.push_back(row);
            hasAnchor = hasAnchor || (anchoredRow[static_cast<size_t>(row)] != 0u);

            for (const int nbr : neighbors[static_cast<size_t>(row)])
            {
                if (nbr < 0 || nbr >= N)
                    continue;
                if (visited[static_cast<size_t>(nbr)] != 0u)
                    continue;
                visited[static_cast<size_t>(nbr)] = 1u;
                q.push(nbr);
            }
        }

        if (hasAnchor)
            continue;

        if (checkedOperator)
            throw std::runtime_error("fluxCorrections: unanchored component; compatibility correction/gauge pin is forbidden");

        ++floatingComponentCount;

        double totalB = 0.0;
        double totalV = 0.0;
        for (const int row : componentRows)
        {
            totalB += bVec[row];
            totalV += volOfRow[static_cast<size_t>(row)];
        }
        if (totalV > 0.0)
        {
            const double meanBperV = totalB / totalV;
            for (const int row : componentRows)
                bVec[row] -= meanBperV * volOfRow[static_cast<size_t>(row)];
        }

        int pinRow = componentRows.front();
        for (const int row : componentRows)
        {
            if (volOfRow[static_cast<size_t>(row)] > 0.0)
            {
                pinRow = row;
                break;
            }
        }

        trips.emplace_back(pinRow, pinRow, 1.0);
        diag[static_cast<size_t>(pinRow)] += 1.0;
        bVec[pinRow] = 0.0;
        anchoredRow[static_cast<size_t>(pinRow)] = 1u;
    }

    if ((floatingComponentCount > 0u) && (opt.verbose || ConsoleControl::enabled(ConsoleControl::Level::kDebug)))
    {
        std::cout << "[PolarizationPoissonFVSolverCellCentered] stabilized "
                  << floatingComponentCount
                  << " floating component(s) with compatibility correction + gauge pin.\n";
    }

    // Safety pin for zero/invalid diagonal.
    for (int i = 0; i < N; ++i)
    {
        const double d = diag[static_cast<size_t>(i)];
        if (!std::isfinite(d) || !(d > 0.0))
        {
            if (checkedOperator)
                throw std::runtime_error("fluxCorrections: non-positive/non-finite assembled diagonal; safety pin is forbidden");
            trips.emplace_back(i, i, 1.0);
            diag[static_cast<size_t>(i)] += 1.0;
            bVec[i] = 0.0;
        }
    }

    A.resize(N,N);
    A.setFromTriplets(trips.begin(), trips.end());
    A.makeCompressed();
    }else{
        A=structure->pattern;
        // Eigen 3.4.1 first assigns a duplicate key, then adds subsequent
        // triplets in input order. Each row owns all of its CSC value slots.
        const auto fillRow=[&](std::size_t row){const int ci=cellIndexOfRow[row];double right=rhs.value[ci],diagonal=0;
            const auto range=graph.segmentRangeForCellIndex(ci);
            for(auto sid=range.first;sid<range.second;++sid){const auto& seg=segs[sid];const auto& correction=(*fluxCorrections)[sid];
                auto term=structure->termOffsets[sid];
                const auto put=[&](double value){const auto slot=structure->slots[term];
                    if(structure->first[term])A.valuePtr()[slot]=value;else A.valuePtr()[slot]+=value;++term;};
                for(const auto& coefficient:correction.coefficients){put(coefficient.value);
                    if(coefficient.cellIndex==ci)diagonal+=coefficient.value;}
                right-=correction.constant;
                const double area=opt.axisymmetric?CartesianMeasure::face(seg.Af,seg.faceCenter.x,true):1.0*seg.Af;
                if(!(seg.df>0)||!(area>=0)||area==0)continue;
                const double kappa=clampKappa((*kappaBySegment)[sid]);const auto& boundary=(*fittedBoundaries)[sid];
                if(boundary.type==FittedBoundaryType::HomogeneousNeumann)continue;
                if(boundary.type==FittedBoundaryType::Dirichlet){
                    const double distance=(seg.faceCenter-mesh.cells(seg.K).centerCoord()).dotProduct(seg.n);
                    const double c=kappa*area/distance;put(c);diagonal+=c;right+=c*boundary.value;
                }else{const double c=kappa*area/seg.df;put(c);put(-c);diagonal+=c;}
            }
            if(!std::isfinite(diagonal)||!(diagonal>0))
                throw std::runtime_error("fluxCorrections: non-positive/non-finite assembled diagonal; safety pin is forbidden");
            return right;
        };
        auto status=IntegratedRowStatus::Literal;
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
        unsigned workers=0;
        status=tryIntegratedRows(static_cast<std::size_t>(N),fillRow,[&](const std::vector<double>& values){
            for(int row=0;row<N;++row)bVec[row]=values[static_cast<std::size_t>(row)];
        },workers);
#endif
        if(status!=IntegratedRowStatus::Parallel)
            for(int row=0;row<N;++row)bVec[row]=fillRow(static_cast<std::size_t>(row));
    }
    ExpandedResidualRoute expandedResidualRoute;
    if(checkedOperator&&integratedFlux){
        if(structure&&!structure->residualOffsets.empty())expandedResidualRoute.borrow(structure->residualOffsets,structure->residualEntries);
        else expandedResidualRoute.begin(A);
    }
    if (checkedOperator)
    {
        for (Eigen::Index entry = 0; entry < A.nonZeros(); ++entry){
            if (!std::isfinite(A.valuePtr()[entry]))
                throw std::runtime_error("fluxCorrections: assembled matrix contains non-finite entries");
            expandedResidualRoute.count(A.innerIndexPtr()[entry]);
        }
        if (!bVec.allFinite())
            throw std::runtime_error("fluxCorrections: assembled RHS contains non-finite entries");
        for (int row = 0; row < N; ++row)
            if (!(A.coeff(row, row) > 0.0))
                throw std::runtime_error("fluxCorrections: complete assembled diagonal must be positive; no safety pin is allowed");
    }
    expandedResidualRoute.finish(A);
    const bool useIterativeCG = opt.useCG || opt.useIC;
    const bool needPatternSig = !checkedOperator && (
        (useIterativeCG && (opt.reuseIterativePattern || opt.useTemporalExtrapolatedGuess)) ||
        ((!useIterativeCG) && (opt.reuseDirectPattern || opt.freezeDirectMatrixUntilPatternChange || opt.useSuperLU)));
    const std::uint64_t patternSig = needPatternSig ? sparsePatternSignature(A) : 0ull;

    const auto tAsm1 = Clock::now();
    const double tAssemble_s = secondsBetween(tAsm0, tAsm1);

    double tFactor_s = 0.0;
    double tSolve_s = 0.0;
    bool usedSuperLU = false;
    bool usedGMRES = false;
    int gmresIterations = 0, gmresRestarts = 0, gmresMaximumWorkers = 0;
    bool usedCG = false;
    const bool cgAttempted = useIterativeCG;
    bool cgFallbackToDirect = false;
    int cgIterations = 0;
    double cgResidual = 0.0;
    bool cgPatternReused = false;
    bool cgTemporalGuessUsed = false;
    bool directPatternReused = false;
    bool directMatrixFrozen = false;
    bool fittedFactorizationReused = false;
    std::string directFailReason;
    bool ok = false;
    Eigen::VectorXd x;
    std::vector<std::array<double,3>> expandedPotential;
    IntegratedRowReceipt orderedRows;
    ExpandedResidualReceipt expandedResidual;

    // Both explicit backends retain the original two exact-face defect passes.
    auto correctIntegratedPotential = [&](auto&& solveDelta) {
        if (!(integratedRhs && checkedOperator)) return true;
        // Integrated charge can have a much smaller local
        // scale than the anchored boundary RHS.  The global
        // residual alone need not resolve those interior rows.
        // Correct the SAME matrix with its SAME fixed operator; no
        // geometry, coefficients, charge or boundary is rebuilt.
        constexpr int maximumDefectCorrections = 2;
        if (integratedFlux)
        {
            expandedPotential.resize(static_cast<std::size_t>(N));
            for (int row = 0; row < N; ++row) expandedPotential[row] = {x[row], 0, 0};
        }
        for (int pass = 0; pass < maximumDefectCorrections; ++pass)
        {
            // Preserve the supplied operator's residual before
            // sparse row compression or dot-product rounding.
            // Each defect rounds once for the SAME fixed operator
            // solve; the original two passes and full-A
            // residual gate remain unchanged.
            namespace Exact = LAMP::Geometry::Detail::FiniteVolume;
            std::vector<Exact::Dyadic> exactDefect(static_cast<std::size_t>(N));
            Eigen::VectorXd parallelDefect;
            if (integratedFlux)
            {
                // Preserve the actual uncompressed shared-face
                // action. Rounded sparse diagonal sums need not
                // retain its exact constant-potential kernel.
                // The original integrated RHS excludes fitted
                // boundary constants, already present in H.
                Eigen::VectorXd orderedDefect;
                unsigned actualWorkers=0;
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
                const auto status=tryIntegratedRows(static_cast<std::size_t>(N),[&](std::size_t compact){
                    return integratedDefectRow(compact,mesh,graph,*integratedFlux,expandedPotential,rowOfCellIndex,cellIndexOfRow,rhs);
                },[&](const std::vector<double>& values){
                    orderedDefect.resize(N);for(int row=0;row<N;++row)orderedDefect[row]=values[static_cast<std::size_t>(row)];
                },actualWorkers,expandedResidualRoute.retainedBytes());
#else
                const auto status=IntegratedRowStatus::Literal;
#endif
                orderedRows.available=true;
                if(status==IntegratedRowStatus::Parallel){
                    parallelDefect=std::move(orderedDefect);orderedRows.parallelDefectRows+=static_cast<std::size_t>(N);
                    orderedRows.maximumWorkers=std::max(orderedRows.maximumWorkers,actualWorkers);
                }else{
                    if(status==IntegratedRowStatus::Fallback)++orderedRows.fallbackPhases;
                for (int row = 0; row < N; ++row) exactDefect[row].add(rhs.value[cellIndexOfRow[row]]);
                for (std::size_t sid = 0; sid < segs.size(); ++sid)
                {
                    const auto& factor = (*integratedFlux)[sid];
                    Exact::Dyadic h(factor.affine.constant);
                    for (const auto& coefficient : factor.affine.coefficients)
                        for (double component : expandedPotential[rowOfCellIndex[coefficient.cellIndex]])
                            h.add(Exact::Dyadic::product({coefficient.value, component}));
                    const double value = static_cast<double>(h.rounded());
                    const double primary = static_cast<double>(static_cast<long double>(factor.primary) * value);
                    const double secondary = factor.secondary == 0 ? 0 :
                        static_cast<double>(static_cast<long double>(factor.secondary) * value);
                    if (!(std::isfinite(value) && std::isfinite(primary) && std::isfinite(secondary)))
                        throw std::runtime_error("nonfinite uncompressed integrated face action");
                    const int row = rowOfCellIndex[mesh.cells(segs[sid].K).index()];
                    exactDefect[row].add(primary, -1);
                    auto transfer = Exact::Dyadic::product({factor.weight, secondary});
                    if (factor.secondaryTransfer == IntegratedFluxFactor::SecondaryTransfer::MinimumBinary64QuantumNearestEven)
                        transfer = transfer.roundedBinary64Quantum();
                    exactDefect[row].add(transfer, -1);
                }
                }
                ++orderedRows.defectPhases;
            }
            else
            {
                for (int row = 0; row < N; ++row) exactDefect[row].add(bVec[row]);
                for (int column = 0; column < A.outerSize(); ++column)
                    for (Eigen::SparseMatrix<double>::InnerIterator entry(A, column); entry; ++entry)
                        exactDefect[entry.row()].add(Exact::Dyadic::product({entry.value(), x[entry.col()]}), -1);
            }
            Eigen::VectorXd defect(N);
            if(parallelDefect.size()==N){for(int row=0;row<N;++row)defect[row]=parallelDefect[row];}
            else{
            for (int row = 0; row < N; ++row)
                defect[row] = Exact::roundDouble(exactDefect[row]).value;
            }
            if (!defect.allFinite())
            {
                ok = false;
                directFailReason = "integrated Poisson defect is non-finite";
                break;
            }
            if (defect.isZero(0.0)) break;
            Eigen::VectorXd delta;
            if (!solveDelta(defect, delta) || !delta.allFinite())
            {
                ok = false;
                directFailReason = "same-operator integrated Poisson defect solve failed" +
                    (directFailReason.empty() ? std::string() : ": " + directFailReason);
                break;
            }
            if (integratedFlux)
            {
                // Retain each of the already computed two
                // binary64 deltas, including sub-ULP values.
                // The fixed three-term sum is authoritative;
                // x remains its nearest scalar projection for
                // the original assembled-matrix/output gates.
                unsigned actualWorkers=0;
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
                const auto status=tryIntegratedRows(static_cast<std::size_t>(N),[&](std::size_t row){
                    return integratedProjectionRow(row,pass,expandedPotential,delta);
                },[&](const std::vector<double>& values){
                    // Plain scalar assignments cannot fail after this join.
                    for(int row=0;row<N;++row){expandedPotential[row][pass+1]=delta[row];x[row]=values[static_cast<std::size_t>(row)];}
                },actualWorkers,expandedResidualRoute.retainedBytes());
#else
                const auto status=IntegratedRowStatus::Literal;
#endif
                if(status==IntegratedRowStatus::Parallel){
                    orderedRows.parallelProjectionRows+=static_cast<std::size_t>(N);
                    orderedRows.maximumWorkers=std::max(orderedRows.maximumWorkers,actualWorkers);
                }else{
                    if(status==IntegratedRowStatus::Fallback)++orderedRows.fallbackPhases;
                for (int row = 0; row < N; ++row)
                {
                    auto& terms = expandedPotential[row]; terms[pass + 1] = delta[row];
                    Exact::Dyadic sum; for (double component : terms) sum.add(component);
                    x[row] = Exact::roundDouble(sum).value;
                }
                }
                ++orderedRows.projectionPhases;
            }
            else x += delta;
            if (!x.allFinite())
            {
                ok = false;
                directFailReason = "integrated Poisson defect update is non-finite";
                break;
            }
        }
        return ok;
    };

    if (opt.useGMRES)
    {
        LampLinear::ParallelBlockJacobi preconditioner;
        const auto factorBegin = Clock::now();
        ok = preconditioner.compute(A, opt.gmresBlockRows);
        tFactor_s += secondsBetween(factorBegin, Clock::now());
        gmresMaximumWorkers = preconditioner.maximumWorkers();
        if (!ok) directFailReason = "GMRES parallel block-Jacobi factorization failed; no direct fallback is allowed";
        auto solveSameOperator = [&](const Eigen::VectorXd& systemRhs, Eigen::VectorXd& solution) {
            solution = Eigen::VectorXd::Zero(N);
            LampLinear::GmresReport gmres;
            usedGMRES = true;
            const bool success = LampLinear::solveGmres(A, preconditioner, systemRhs, solution,
                opt.relTol, opt.maxIter, opt.gmresRestart, gmres);
            gmresIterations += gmres.iterations;
            gmresRestarts += gmres.restarts;
            gmresMaximumWorkers = std::max(gmresMaximumWorkers, gmres.maximumWorkers);
            if (!success)
            {
                const double rhsNorm = systemRhs.stableNorm();
                double trueRelative = std::numeric_limits<double>::quiet_NaN();
                if (solution.size() == N && solution.allFinite() && std::isfinite(rhsNorm))
                {
                    const double norm = (A * solution - systemRhs).stableNorm();
                    trueRelative = rhsNorm > 0 ? norm / rhsNorm : norm;
                }
                std::ostringstream reason;
                reason << gmres.message << "; iterations=" << gmres.iterations
                       << " restarts=" << gmres.restarts
                       << " true_relative_residual=" << trueRelative
                       << " reported_relative_residual=" << gmres.relativeResidual;
                directFailReason = reason.str();
            }
            return success;
        };
        if (ok)
        {
            const auto solveBegin = Clock::now();
            ok = solveSameOperator(bVec, x);
            if (ok) ok = correctIntegratedPotential(solveSameOperator);
            tSolve_s += secondsBetween(solveBegin, Clock::now());
        }
    }

    if (useIterativeCG)
    {
        Eigen::VectorXd x0(N);
        const bool canUseTemporalGuess =
            opt.useTemporalExtrapolatedGuess &&
            iterCache_.hasLastX &&
            (iterCache_.xN == N) &&
            (iterCache_.xPatternSig == patternSig) &&
            (iterCache_.lastX.size() == N);

        if (canUseTemporalGuess)
        {
            x0 = iterCache_.lastX;
            if (iterCache_.hasPrevX && (iterCache_.prevX.size() == N))
            {
                const double w = std::max(0.0, opt.temporalExtrapolationOmega);
                if (w > 0.0)
                {
                    x0.noalias() += w * (iterCache_.lastX - iterCache_.prevX);
                    cgTemporalGuessUsed = true;
                }
            }
        }
        else
        {
            for (int i = 0; i < N; ++i)
            {
                const int cellIndex = cellIndexOfRow[i];
                x0[i] = (cellIndex >= 0 && cellIndex < mesh.cellIndexNum()) ? phi.value[cellIndex] : 0.0;
            }
        }

        if (opt.useIC)
        {
            auto& cg = iterCache_.cgIC;
            cg.setMaxIterations(std::max(1, opt.maxIter));
            cg.setTolerance(std::max(1e-16, opt.relTol));

            const bool canReuse =
                opt.reuseIterativePattern &&
                iterCache_.icAnalyzed &&
                iterCache_.icN == N &&
                iterCache_.icNnz == static_cast<long long>(A.nonZeros()) &&
                iterCache_.icPatternSig == patternSig;

            const auto tCgFact0 = Clock::now();
            bool haveFactorization = false;
            if (canReuse)
            {
                cg.factorize(A);
                if (cg.info() == Eigen::Success)
                {
                    haveFactorization = true;
                    cgPatternReused = true;
                }
            }
            if (!haveFactorization)
            {
                cg.compute(A);
                haveFactorization = (cg.info() == Eigen::Success);
                iterCache_.icAnalyzed = haveFactorization;
                if (haveFactorization)
                {
                    iterCache_.icN = N;
                    iterCache_.icNnz = static_cast<long long>(A.nonZeros());
                    iterCache_.icPatternSig = patternSig;
                }
            }
            const auto tCgFact1 = Clock::now();
            const auto tCgSol0 = tCgFact1;
            if (haveFactorization)
            {
                x = cg.solveWithGuess(bVec, x0);
                ok = (cg.info() == Eigen::Success);
                cgIterations = static_cast<int>(cg.iterations());
                cgResidual = cg.error();
            }
            const auto tCgSol1 = Clock::now();
            tFactor_s += secondsBetween(tCgFact0, tCgFact1);
            tSolve_s += secondsBetween(tCgSol0, tCgSol1);
        }
        else
        {
            auto& cg = iterCache_.cg;
            cg.setMaxIterations(std::max(1, opt.maxIter));
            cg.setTolerance(std::max(1e-16, opt.relTol));

            const bool canReuse =
                opt.reuseIterativePattern &&
                iterCache_.cgAnalyzed &&
                iterCache_.cgN == N &&
                iterCache_.cgNnz == static_cast<long long>(A.nonZeros()) &&
                iterCache_.cgPatternSig == patternSig;

            const auto tCgFact0 = Clock::now();
            bool haveFactorization = false;
            if (canReuse)
            {
                cg.factorize(A);
                if (cg.info() == Eigen::Success)
                {
                    haveFactorization = true;
                    cgPatternReused = true;
                }
            }
            if (!haveFactorization)
            {
                cg.compute(A);
                haveFactorization = (cg.info() == Eigen::Success);
                iterCache_.cgAnalyzed = haveFactorization;
                if (haveFactorization)
                {
                    iterCache_.cgN = N;
                    iterCache_.cgNnz = static_cast<long long>(A.nonZeros());
                    iterCache_.cgPatternSig = patternSig;
                }
            }
            const auto tCgFact1 = Clock::now();
            const auto tCgSol0 = tCgFact1;
            if (haveFactorization)
            {
                x = cg.solveWithGuess(bVec, x0);
                ok = (cg.info() == Eigen::Success);
                cgIterations = static_cast<int>(cg.iterations());
                cgResidual = cg.error();
            }
            const auto tCgSol1 = Clock::now();
            tFactor_s += secondsBetween(tCgFact0, tCgFact1);
            tSolve_s += secondsBetween(tCgSol0, tCgSol1);
        }

        if (ok)
        {
            for (int i = 0; i < N; ++i)
            {
                if (!std::isfinite(x[i]))
                {
                    ok = false;
                    break;
                }
            }
        }
        if (ok)
            usedCG = true;
    }

    if (!ok && !opt.useGMRES)
    {
        if (cgAttempted)
            cgFallbackToDirect = true;
        const auto tDirect0 = Clock::now();
        double directSolveWall_s = 0.0;
        const bool canReuseDirectPattern = opt.reuseDirectPattern && needPatternSig;
        const bool canFreezeDirect = opt.freezeDirectMatrixUntilPatternChange && needPatternSig;
        static int s_directDebugCached = -1;
        if (s_directDebugCached < 0)
            s_directDebugCached = envBool("LAMP_POL_DIRECT_FREEZE_DEBUG", false) ? 1 : 0;
        const bool directDebug = (s_directDebugCached == 1);
        static std::uint64_t s_directCallId = 0;
        const std::uint64_t directCallId = directDebug ? (++s_directCallId) : 0;
#ifdef EIGEN_SUPERLU_SUPPORT
        if (opt.useSuperLU)
        {
            // Default to a per-call SuperLU instance unless direct matrix freezing is enabled.
            // This avoids long-lived internal factorization workspaces accumulating across
            // many macro-steps in the common "refactorize every step" path.
            const bool superluLocal = checkedOperator || envBool("LAMP_POL_SUPERLU_LOCAL", !canFreezeDirect);

            auto configureSuperlu = [](auto& solver) {
                solver.options().ColPerm = COLAMD;
                solver.options().DiagPivotThresh = 0.1;
            };

            if (fittedOperator && opt.reuseFittedFactorization)
            {
                fittedFactorizationReused = fittedFactorCache_ &&
                    fittedFactorCache_->matches(A, mesh, cellIndexOfRow);
                if (!fittedFactorizationReused)
                {
                    // Destroy old factors before allocating/factorizing a new
                    // operator. A failed factorization never becomes reusable.
                    fittedFactorCache_.reset();
                    auto candidate = std::make_unique<FittedFactorCache>();
                    candidate->matrix = A;
                    candidate->mesh = &mesh;
                    candidate->cellIndices = cellIndexOfRow;
                    candidate->cellLabels.reserve(cellIndexOfRow.size());
                    for (int cell : cellIndexOfRow)
                        candidate->cellLabels.push_back(mesh.cellIndexToLabel[cell]);
                    configureSuperlu(candidate->solver);
                    candidate->solver.analyzePattern(candidate->matrix);
                    candidate->solver.factorize(candidate->matrix);
                    if (candidate->solver.info() == Eigen::Success)
                        fittedFactorCache_ = std::move(candidate);
                    else
                        directFailReason = "fitted SuperLU factorize failed, info=" +
                            std::to_string(static_cast<int>(candidate->solver.info()));
                }
                if (fittedFactorCache_)
                {
                    const auto tSolve0 = Clock::now();
                    x = fittedFactorCache_->solver.solve(bVec);
                    ok = fittedFactorCache_->solver.info() == Eigen::Success;
                    if (!ok)
                        directFailReason = "fitted SuperLU solve failed, info=" +
                            std::to_string(static_cast<int>(fittedFactorCache_->solver.info()));
                    directSolveWall_s += secondsBetween(tSolve0, Clock::now());
                    usedSuperLU = ok;
                }
            }
            else if (superluLocal)
            {
                Eigen::SuperLU<Eigen::SparseMatrix<double>> superlu;
                configureSuperlu(superlu);
                if (directDebug)
                    std::cout << "[PolPoisson][direct][SuperLU] call=" << directCallId << " LOCAL solver mode" << std::endl;

                if (directDebug)
                    std::cout << "[PolPoisson][direct][SuperLU] call=" << directCallId << " analyzePattern..." << std::endl;
                superlu.analyzePattern(A);
                if (directDebug)
                    std::cout << "[PolPoisson][direct][SuperLU] call=" << directCallId << " analyzePattern done (info="
                              << static_cast<int>(superlu.info()) << ")" << std::endl;

                if (directDebug)
                    std::cout << "[PolPoisson][direct][SuperLU] call=" << directCallId << " factorize..." << std::endl;
                superlu.factorize(A);
                if (directDebug)
                    std::cout << "[PolPoisson][direct][SuperLU] call=" << directCallId << " factorize done (info="
                              << static_cast<int>(superlu.info()) << ")" << std::endl;

                const bool factored = (superlu.info() == Eigen::Success);
                if (!factored)
                    directFailReason = "SuperLU factorize failed, info=" + std::to_string(static_cast<int>(superlu.info()));
                if (factored)
                {
                    const auto tSolve0 = Clock::now();
                    if (directDebug)
                        std::cout << "[PolPoisson][direct][SuperLU] call=" << directCallId << " solve..." << std::endl;
                    x = superlu.solve(bVec);
                    if (directDebug)
                        std::cout << "[PolPoisson][direct][SuperLU] call=" << directCallId << " solve done (info="
                                  << static_cast<int>(superlu.info()) << ")" << std::endl;
                    ok = (superlu.info() == Eigen::Success);
                    if (!ok)
                        directFailReason = "SuperLU solve failed, info=" + std::to_string(static_cast<int>(superlu.info()));
                    if (ok) ok = correctIntegratedPotential([&](const Eigen::VectorXd& residual, Eigen::VectorXd& delta) {
                        delta = superlu.solve(residual);
                        return superlu.info() == Eigen::Success;
                    });
                    directSolveWall_s += secondsBetween(tSolve0, Clock::now());
                }
                if (ok)
                {
                    usedSuperLU = true;
                }
            }
            else
            {
                auto& superlu = directCache_.superLU;
                configureSuperlu(superlu);
            const bool samePattern =
                directCache_.superAnalyzed &&
                sameSparsePattern(directCache_.superMatrix, A);

            const bool reusePatternNow = samePattern && canReuseDirectPattern;
            const bool freezeNow =
                samePattern && canFreezeDirect && directCache_.superFactorized &&
                sameSparseValues(directCache_.superMatrix, A);
            bool factored = directCache_.superFactorized;

            if (directDebug)
            {
                std::cout << "[PolPoisson][direct][SuperLU]"
                          << " call=" << directCallId
                          << " N=" << N
                          << " nnz=" << static_cast<long long>(A.nonZeros())
                          << " samePattern=" << (samePattern ? 1 : 0)
                          << " reusePatternNow=" << (reusePatternNow ? 1 : 0)
                          << " hadFactor=" << (directCache_.superFactorized ? 1 : 0)
                          << " canFreeze=" << (canFreezeDirect ? 1 : 0)
                          << " freezeNow=" << (freezeNow ? 1 : 0)
                          << std::endl;
            }

            if (!freezeNow)
            {
                bool analyzed = reusePatternNow;
                directCache_.superAnalyzed = false;
                directCache_.superFactorized = false;
                if (!analyzed)
                {
                    if (directDebug)
                        std::cout << "[PolPoisson][direct][SuperLU] call=" << directCallId << " analyzePattern..." << std::endl;
                    superlu.analyzePattern(A);
                    if (directDebug)
                        std::cout << "[PolPoisson][direct][SuperLU] call=" << directCallId << " analyzePattern done (info="
                                  << static_cast<int>(superlu.info()) << ")" << std::endl;
                    analyzed = true;
                }
                if (analyzed)
                {
                    if (directDebug)
                    {
                        std::size_t nonFinite = 0u;
                        const double* vals = A.valuePtr();
                        const Eigen::Index nnz = A.nonZeros();
                        for (Eigen::Index k = 0; k < nnz; ++k)
                        {
                            if (!std::isfinite(vals[k]))
                                ++nonFinite;
                        }
                        if (nonFinite > 0u)
                        {
                            std::cout << "[PolPoisson][direct][SuperLU] call=" << directCallId
                                      << " WARNING: matrix has non-finite values, skipping SuperLU" << std::endl;
                            analyzed = false;
                        }
                    }
                    if (analyzed && directDebug)
                        std::cout << "[PolPoisson][direct][SuperLU] call=" << directCallId << " factorize..." << std::endl;
                    superlu.factorize(A);
                    if (directDebug)
                        std::cout << "[PolPoisson][direct][SuperLU] call=" << directCallId << " factorize done (info="
                                  << static_cast<int>(superlu.info()) << ")" << std::endl;
                    factored = (superlu.info() == Eigen::Success);
                    if (!factored && directFailReason.empty())
                    {
                        std::ostringstream oss;
                        oss << "SuperLU factorize failed, info=" << static_cast<int>(superlu.info());
                        directFailReason = oss.str();
                    }
                    if (factored)
                        directCache_.superMatrix = A;
                    directCache_.superAnalyzed = factored;
                    directCache_.superFactorized = factored;
                    if (factored && reusePatternNow)
                        directPatternReused = true;
                }
                else
                {
                    factored = false;
                    directCache_.superFactorized = false;
                }
            }
            else
            {
                directPatternReused = true;
                directMatrixFrozen = true;
            }

            if (factored)
            {
                const auto tSolve0 = Clock::now();
                if (directDebug)
                    std::cout << "[PolPoisson][direct][SuperLU] call=" << directCallId << " solve..." << std::endl;
                x = superlu.solve(bVec);
                if (directDebug)
                    std::cout << "[PolPoisson][direct][SuperLU] call=" << directCallId << " solve done (info="
                              << static_cast<int>(superlu.info()) << ")" << std::endl;
                ok = (superlu.info() == Eigen::Success);
                directSolveWall_s += secondsBetween(tSolve0, Clock::now());
                if (!ok && directFailReason.empty())
                {
                    std::ostringstream oss;
                    oss << "SuperLU solve failed, info=" << static_cast<int>(superlu.info());
                    directFailReason = oss.str();
                }
            }
            if (ok)
            {
                usedSuperLU = true;
            }
            else
            {
                usedSuperLU = false;
                directCache_.superFactorized = false;
            }
            }
        }
#endif

        if (opt.requireSuperLU && ok && (!usedSuperLU || !x.allFinite()))
        {
            ok = false;
            usedSuperLU = false;
            directFailReason = "required SuperLU solve returned non-finite values or did not use SuperLU";
        }

        if (!ok && !opt.requireSuperLU)
        {
            auto& sparselu = directCache_.sparseLU;
            const bool samePattern =
                directCache_.sparseAnalyzed &&
                sameSparsePattern(directCache_.sparseMatrix, A);

            const bool reusePatternNow = samePattern && canReuseDirectPattern;
            const bool freezeNow =
                samePattern && canFreezeDirect && directCache_.sparseFactorized &&
                sameSparseValues(directCache_.sparseMatrix, A);
            bool factored = directCache_.sparseFactorized;

            if (directDebug)
            {
                std::cout << "[PolPoisson][direct][SparseLU]"
                          << " call=" << directCallId
                          << " N=" << N
                          << " nnz=" << static_cast<long long>(A.nonZeros())
                          << " samePattern=" << (samePattern ? 1 : 0)
                          << " reusePatternNow=" << (reusePatternNow ? 1 : 0)
                          << " hadFactor=" << (directCache_.sparseFactorized ? 1 : 0)
                          << " canFreeze=" << (canFreezeDirect ? 1 : 0)
                          << " freezeNow=" << (freezeNow ? 1 : 0)
                          << std::endl;
            }

            if (!freezeNow)
            {
                bool analyzed = reusePatternNow;
                directCache_.sparseAnalyzed = false;
                directCache_.sparseFactorized = false;
                if (!analyzed)
                {
                    sparselu.analyzePattern(A);
                    analyzed = true;
                }
                if (analyzed)
                {
                    sparselu.factorize(A);
                    factored = (sparselu.info() == Eigen::Success);
                    if (!factored && directFailReason.empty())
                    {
                        std::ostringstream oss;
                        oss << "SparseLU factorize failed, info=" << static_cast<int>(sparselu.info());
                        directFailReason = oss.str();
                    }
                    if (factored)
                        directCache_.sparseMatrix = A;
                    directCache_.sparseAnalyzed = factored;
                    directCache_.sparseFactorized = factored;
                    if (factored && reusePatternNow)
                        directPatternReused = true;
                }
                else
                {
                    factored = false;
                    directCache_.sparseFactorized = false;
                }
            }
            else
            {
                directPatternReused = true;
                directMatrixFrozen = true;
            }

            if (factored)
            {
                const auto tSolve0 = Clock::now();
                x = sparselu.solve(bVec);
                ok = (sparselu.info() == Eigen::Success);
                directSolveWall_s += secondsBetween(tSolve0, Clock::now());
                if (!ok && directFailReason.empty())
                {
                    std::ostringstream oss;
                    oss << "SparseLU solve failed, info=" << static_cast<int>(sparselu.info());
                    directFailReason = oss.str();
                }
            }
            if (!ok)
            {
                directCache_.sparseFactorized = false;
            }
        }
        const double directTotal_s = secondsBetween(tDirect0, Clock::now());
        tSolve_s += directSolveWall_s;
        tFactor_s += std::max(0.0, directTotal_s - directSolveWall_s);
        usedCG = false;
        if (!cgAttempted)
        {
            cgIterations = 0;
            cgResidual = 0.0;
        }
    }

    double correctedResidual = 0.0;
    if (checkedOperator && ok)
    {
        // Check the operator actually factorized, including every integrated
        // flux correction and boundary constant. Do not use CG's estimate or a
        // residual of only the underlying two-point diffusion stencil.
        ok = fullResidualAccepted(A, x, bVec, usedSuperLU || usedGMRES, opt.relTol,
            correctedResidual, directFailReason);
        if (ok && integratedFlux)
        {
            // Bind the same relTol to the actual physical potential sum as
            // well as retaining the old projected-x check above.
            namespace Exact = LAMP::Geometry::Detail::FiniteVolume;
            expandedResidual.available=true;expandedResidual.rows=static_cast<std::size_t>(N);
            expandedResidual.routeSeconds=expandedResidualRoute.seconds;
            expandedResidual.scratchBytes=expandedResidualRoute.peakBytes;
            Eigen::VectorXd residual;
            auto status=IntegratedRowStatus::Literal;
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
            Clock::time_point residualStarted{};
            const bool residualClock=integratedClockSample(residualStarted);
            if(expandedResidualRoute.ready){
                status=tryIntegratedRows(static_cast<std::size_t>(N),[&](std::size_t row){
                    return integratedExpandedResidualRow(row,A,expandedResidualRoute.rowOffsets(),expandedResidualRoute.rowEntries(),expandedPotential,bVec);
                },[&](const std::vector<double>& values){
                    residual.resize(N);for(int row=0;row<N;++row)residual[row]=values[static_cast<std::size_t>(row)];
                },expandedResidual.workers,expandedResidualRoute.retainedBytes(),&expandedResidual.scratchBytes);
            }
#endif
            expandedResidual.parallelUsed=status==IntegratedRowStatus::Parallel;
            expandedResidual.parallelFallback=expandedResidualRoute.failed||status==IntegratedRowStatus::Fallback;
            if(!expandedResidual.parallelUsed){
            std::vector<Exact::Dyadic> exactResidual(static_cast<std::size_t>(N));
            for (int row = 0; row < N; ++row) exactResidual[row].add(bVec[row], -1);
            for (int column = 0; column < A.outerSize(); ++column)
                for (Eigen::SparseMatrix<double>::InnerIterator entry(A, column); entry; ++entry)
                    for (double component : expandedPotential[entry.col()])
                        exactResidual[entry.row()].add(Exact::Dyadic::product({entry.value(), component}));
            Eigen::VectorXd literalResidual(N);
            for (int row = 0; row < N; ++row) literalResidual[row] = Exact::roundDouble(exactResidual[row]).value;
            residual=std::move(literalResidual);
            }
#if defined(_OPENMP) && defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
            expandedResidual.timingAvailable=residualClock&&expandedResidualRoute.timingAvailable&&
                integratedClockElapsed(residualStarted,expandedResidual.seconds);
#endif
            const double norm = residual.norm(), rhsNorm = bVec.norm();
            const double relative = rhsNorm > 0 ? norm / rhsNorm : norm;
            if (!residual.allFinite() || !std::isfinite(relative) || relative > opt.relTol)
            {
                ok = false;
                directFailReason = "integrated expanded potential full-A residual exceeds original relTol";
            }
            correctedResidual = std::max(correctedResidual, relative);
        }
    }

    if (!ok)
    {
        if (report)
        {
            *report = {};
            report->success = false;
            report->finalResidual = correctedResidual;
            if (opt.useGMRES)
                report->message = "GMRES solve failed; direct fallback is disabled: " +
                    (directFailReason.empty() ? std::string("block-Jacobi/iteration failure") : directFailReason);
            else if (opt.requireSuperLU)
                report->message = "requireSuperLU solve failed; SparseLU fallback is disabled: " +
                    (directFailReason.empty() ? std::string("SuperLU factorization/solve failed") : directFailReason);
            else if (!directFailReason.empty())
                report->message = "Polarization Poisson solve failed (CG and direct fallback): " + directFailReason;
            else
                report->message = "Polarization Poisson solve failed (CG and direct fallback).";
            report->timeAssemble_s = tAssemble_s;
            report->timeFactor_s = tFactor_s;
            report->timeSolve_s = tSolve_s;
            report->timeTotal_s = secondsSince(tBegin);
            report->N = N;
            report->nnz = static_cast<long long>(A.nonZeros());
            report->usedSuperLU = usedSuperLU;
            report->usedGMRES = usedGMRES;
            report->iterations = report->gmresIterations = gmresIterations;
            report->gmresRestarts = gmresRestarts;
            report->gmresMaximumWorkers = gmresMaximumWorkers;
            report->usedCG = usedCG;
            report->cgAttempted = cgAttempted;
            report->cgFallbackToDirect = cgFallbackToDirect;
            report->cgIterations = cgIterations;
            report->cgResidual = cgResidual;
            report->cgPatternReused = cgPatternReused;
            report->cgTemporalGuessUsed = cgTemporalGuessUsed;
            report->directPatternReused = directPatternReused;
            report->directMatrixFrozen = directMatrixFrozen;
            report->fittedFactorizationReused = fittedFactorizationReused;
            report->integratedExpandedResidualAvailable=expandedResidual.available;
            report->integratedExpandedResidualParallelUsed=expandedResidual.parallelUsed;
            report->integratedExpandedResidualParallelFallback=expandedResidual.parallelFallback;
            report->integratedExpandedResidualWorkers=expandedResidual.workers;
            report->integratedExpandedResidualRows=expandedResidual.rows;
            report->integratedExpandedResidualScratchBytes=expandedResidual.scratchBytes;
            report->integratedExpandedResidualTimingAvailable=expandedResidual.timingAvailable;
            report->integratedExpandedResidualRouteSeconds=expandedResidual.routeSeconds;
            report->integratedExpandedResidualSeconds=expandedResidual.seconds;
        }
        return false;
    }

    if (useIterativeCG && opt.useTemporalExtrapolatedGuess)
    {
        const bool haveCompatibleLast =
            iterCache_.hasLastX &&
            (iterCache_.xN == N) &&
            (iterCache_.xPatternSig == patternSig) &&
            (iterCache_.lastX.size() == N);
        if (haveCompatibleLast)
        {
            iterCache_.prevX = iterCache_.lastX;
            iterCache_.hasPrevX = true;
        }
        else
        {
            iterCache_.prevX.resize(0);
            iterCache_.hasPrevX = false;
        }
        iterCache_.lastX = x;
        iterCache_.xN = N;
        iterCache_.xPatternSig = patternSig;
        iterCache_.hasLastX = true;
    }

    if (preparedRhs)
    {
        if (!fittedFactorCache_)
            throw std::runtime_error("fitted RHS plan has no successful SuperLU factorization");
        fittedFactorCache_->rhsPlan = std::move(preparedRhs);
    }

    if (integratedFlux)
    {
        std::vector<std::array<double,3>> nativeTerms(static_cast<std::size_t>(mesh.cellIndexNum()));
        for (int row = 0; row < N; ++row) nativeTerms[cellIndexOfRow[row]] = expandedPotential[row];
        potentialTerms->swap(nativeTerms); // The public wrapper still owns this staged output.
    }

    // New-port validation, factorization, solve and residual failures have not
    // touched phi. Allocate a missing output only at the successful commit.
    if (correctedOperator && (phi.size() != mesh.cellIndexNum() || !phi.value))
        phi.initialize(mesh, 0.0);

    for (int i = 0; i < N; ++i)
    {
        const int cellIndex = cellIndexOfRow[i];
        if (cellIndex >= 0 && cellIndex < mesh.cellIndexNum())
            phi.value[cellIndex] = x[i];
    }

    if (report)
    {
        *report = {};
        report->success = true;
        report->finalResidual = correctedResidual;
        report->message = fittedOperator && correctedOperator ?
            "Solved face-fitted active-leaf Poisson with integrated flux corrections; full-operator residual verified." :
            correctedOperator ?
            "Solved variable-coefficient Poisson with integrated flux corrections; full-operator residual verified." :
            fittedOperator ? "Solved face-fitted active-leaf Poisson; full-operator residual verified." :
            "Solved variable-coefficient Poisson (one-shot).";
        report->timeAssemble_s = tAssemble_s;
        report->timeFactor_s = tFactor_s;
        report->timeSolve_s = tSolve_s;
        report->timeTotal_s = secondsSince(tBegin);
        report->N = N;
        report->nnz = static_cast<long long>(A.nonZeros());
        report->usedSuperLU = usedSuperLU;
        report->usedGMRES = usedGMRES;
        report->iterations = report->gmresIterations = gmresIterations;
        report->gmresRestarts = gmresRestarts;
        report->gmresMaximumWorkers = gmresMaximumWorkers;
        report->integratedOrderedRowsAvailable=orderedRows.available;
        report->integratedDefectPhases=orderedRows.defectPhases;
        report->integratedProjectionPhases=orderedRows.projectionPhases;
        report->integratedParallelDefectRows=orderedRows.parallelDefectRows;
        report->integratedParallelProjectionRows=orderedRows.parallelProjectionRows;
        report->integratedRowsMaximumWorkers=orderedRows.maximumWorkers;
        report->integratedRowsFallbackPhases=orderedRows.fallbackPhases;
        report->integratedExpandedResidualAvailable=expandedResidual.available;
        report->integratedExpandedResidualParallelUsed=expandedResidual.parallelUsed;
        report->integratedExpandedResidualParallelFallback=expandedResidual.parallelFallback;
        report->integratedExpandedResidualWorkers=expandedResidual.workers;
        report->integratedExpandedResidualRows=expandedResidual.rows;
        report->integratedExpandedResidualScratchBytes=expandedResidual.scratchBytes;
        report->integratedExpandedResidualTimingAvailable=expandedResidual.timingAvailable;
        report->integratedExpandedResidualRouteSeconds=expandedResidual.routeSeconds;
        report->integratedExpandedResidualSeconds=expandedResidual.seconds;

        report->usedCG = usedCG;
        report->cgAttempted = cgAttempted;
        report->cgFallbackToDirect = cgFallbackToDirect;
        report->cgIterations = cgIterations;
        report->cgResidual = cgResidual;
        report->cgPatternReused = cgPatternReused;
        report->cgTemporalGuessUsed = cgTemporalGuessUsed;
        report->directPatternReused = directPatternReused;
        report->directMatrixFrozen = directMatrixFrozen;
        report->fittedFactorizationReused = fittedFactorizationReused;
    }

    return true;
}

void PolarizationPoissonFVSolverCellCentered::clearDirectCache() const
{
    fittedFactorCache_.reset();
    iterCache_.icN = -1;
    iterCache_.icNnz = -1;
    iterCache_.icPatternSig = 0;
    iterCache_.icAnalyzed = false;
    iterCache_.cgN = -1;
    iterCache_.cgNnz = -1;
    iterCache_.cgPatternSig = 0;
    iterCache_.cgAnalyzed = false;
    iterCache_.xN = -1;
    iterCache_.xPatternSig = 0;
    iterCache_.hasLastX = false;
    iterCache_.hasPrevX = false;
    iterCache_.lastX.resize(0);
    iterCache_.prevX.resize(0);

#ifdef EIGEN_SUPERLU_SUPPORT
    directCache_.superLU.~SuperLU();
    new (&directCache_.superLU) Eigen::SuperLU<Eigen::SparseMatrix<double>>();
    directCache_.superMatrix.resize(0, 0);
    directCache_.superAnalyzed = false;
    directCache_.superFactorized = false;
#endif
    directCache_.sparseLU.~SparseLU();
    new (&directCache_.sparseLU) Eigen::SparseLU<Eigen::SparseMatrix<double>>();
    directCache_.sparseMatrix.resize(0, 0);
    directCache_.sparseAnalyzed = false;
    directCache_.sparseFactorized = false;
}
