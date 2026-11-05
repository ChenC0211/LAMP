#pragma once

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <algorithm>
#include <atomic>
#include <cerrno>
#include <cfenv>
#include <cmath>
#include <memory>
#include <string>
#include <vector>
#if defined(_OPENMP)
#include <omp.h>
#endif
#if defined(_M_X64) || defined(__SSE2__)
#include <xmmintrin.h>
#endif

namespace LampLinear {

using Sparse = Eigen::SparseMatrix<double>;

struct GmresReport {
    int iterations=0,restarts=0,maximumWorkers=0;
    double relativeResidual=0;
    std::string message;
};

namespace BlockJacobiDetail {

// Numerical work may change arithmetic flags, but must not leave a worker's
// rounding controls, floating-point state or errno changed for its next owner.
class FloatingState {
    std::fenv_t environment_{};
    int error_=errno;
    bool valid_=std::fegetenv(&environment_)==0;
#if defined(_M_X64) || defined(__SSE2__)
    unsigned csr_=_mm_getcsr();
#endif
public:
    bool enter() const noexcept {
        if(!valid_||std::fesetenv(&environment_)!=0)return false;
#if defined(_M_X64) || defined(__SSE2__)
        _mm_setcsr(csr_);
#endif
        errno=error_;
        return true;
    }
    ~FloatingState() noexcept {enter();}
};

inline int workersFor(Eigen::Index jobs) noexcept {
#if defined(_OPENMP)
    if(omp_in_parallel())return 1;
    return static_cast<int>(std::max<Eigen::Index>(1,std::min<Eigen::Index>(jobs,omp_get_max_threads())));
#else
    (void)jobs;return 1;
#endif
}

// Each job owns its complete output slot. Small local LU/triangular operations
// do not call a threaded BLAS; the application owns the OpenMP thread budget.
template<class Function>
inline bool parallelJobs(Eigen::Index jobs,int team,Function&& function,int& observed) {
    const FloatingState caller;
    std::atomic<bool> failed{false};
    observed=1;
#if defined(_OPENMP)
#pragma omp parallel num_threads(team) if(team>1)
#else
    (void)team;
#endif
    {
        const FloatingState restore;
        const bool ready=caller.enter();
        if(!ready)failed.store(true,std::memory_order_relaxed);
#if defined(_OPENMP)
#pragma omp single
        {observed=omp_get_num_threads();}
#pragma omp for schedule(static)
#endif
        for(Eigen::Index job=0;job<jobs;++job){
            if(!ready)continue;
            try {
                if(!function(job))failed.store(true,std::memory_order_relaxed);
            }catch(...){failed.store(true,std::memory_order_relaxed);}
        }
    }
    return !failed.load(std::memory_order_relaxed);
}

inline void recordWorkers(std::atomic<int>& maximum,int observed) noexcept {
    int prior=maximum.load(std::memory_order_relaxed);
    while(prior<observed&&!maximum.compare_exchange_weak(prior,observed,std::memory_order_relaxed)){}
}

using RowSparse = Eigen::SparseMatrix<double,Eigen::RowMajor>;

inline bool multiply(const RowSparse& matrix,const Eigen::Ref<const Eigen::VectorXd>& rhs,
                     Eigen::Ref<Eigen::VectorXd> out,int& observed) {
    // About 512 rows per worker avoids launching every configured thread for
    // the small adaptive systems. Conversion to row storage is once per solve.
    const int team=workersFor((matrix.rows()+511)/512);
    return parallelJobs(matrix.rows(),team,[&](Eigen::Index row){
        double sum=0;
        for(RowSparse::InnerIterator e(matrix,row);e;++e)sum+=e.value()*rhs[e.col()];
        out[row]=sum;
        return std::isfinite(sum);
    },observed);
}

} // namespace BlockJacobiDetail

// Fixed linear right preconditioner: contiguous core blocks, optionally with
// one fixed graph-neighbour layer (restricted additive Schwarz). Each local
// sparse LU writes only its core, without accumulating or weighting overlaps.
// This is local block LU, not a global direct fallback and not ILU.
class ParallelBlockJacobi {
    using Factor=Eigen::SparseLU<Sparse,Eigen::COLAMDOrdering<int>>;
    struct Block {
        Eigen::Index begin=0,size=0;
        std::unique_ptr<Factor> factor;
        std::vector<Eigen::Index> neighbours;
        std::vector<int> corePositions;
        mutable Eigen::VectorXd rhsScratch,solutionScratch;
    };
    Eigen::Index size_=0;
    std::vector<Block> blocks_;
    mutable std::atomic<int> maximumWorkers_{0};
    bool ready_=false,overlap_=false;
public:
    bool compute(const Sparse& matrix,int blockRows=256,bool oneLayerOverlap=false) {
        const BlockJacobiDetail::FloatingState restore;
        ready_=false;blocks_.clear();size_=matrix.rows();overlap_=oneLayerOverlap;
        maximumWorkers_.store(0,std::memory_order_relaxed);
        if(blockRows<=0||matrix.rows()!=matrix.cols())return false;
        try {
            for(Eigen::Index column=0;column<matrix.outerSize();++column)
                for(Sparse::InnerIterator e(matrix,column);e;++e)
                    if(!std::isfinite(e.value()))return false;
            if(size_==0){ready_=true;return true;}
            const BlockJacobiDetail::RowSparse rows(matrix);
            const Eigen::Index count=(size_-1)/blockRows+1;
            blocks_.resize(static_cast<std::size_t>(count));
            int observed=0;
            const bool success=BlockJacobiDetail::parallelJobs(count,BlockJacobiDetail::workersFor(count),
                [&](Eigen::Index id){
                    auto& block=blocks_[static_cast<std::size_t>(id)];
                    block.begin=id*blockRows;block.size=std::min<Eigen::Index>(blockRows,size_-block.begin);
                    const Eigen::Index end=block.begin+block.size;
                    std::vector<Eigen::Triplet<double>> entries;
                    Eigen::Index localSize=block.size;
                    if(overlap_){
                        auto& neighbours=block.neighbours;
                        neighbours.reserve(static_cast<std::size_t>(block.size));
                        // Expand from the original core exactly once, using
                        // real nonzero row/column adjacency, never new halos.
                        for(Eigen::Index row=block.begin;row<end;++row){
                            neighbours.push_back(row);
                            for(BlockJacobiDetail::RowSparse::InnerIterator e(rows,row);e;++e)
                                if(e.value()!=0)neighbours.push_back(e.col());
                            for(Sparse::InnerIterator e(matrix,row);e;++e)
                                if(e.value()!=0)neighbours.push_back(e.row());
                        }
                        std::sort(neighbours.begin(),neighbours.end());
                        neighbours.erase(std::unique(neighbours.begin(),neighbours.end()),neighbours.end());
                        localSize=static_cast<Eigen::Index>(neighbours.size());
                        block.corePositions.resize(static_cast<std::size_t>(block.size));
                        for(Eigen::Index i=0;i<block.size;++i)
                            block.corePositions[static_cast<std::size_t>(i)]=static_cast<int>(
                                std::lower_bound(neighbours.begin(),neighbours.end(),block.begin+i)-neighbours.begin());
                        for(Eigen::Index localRow=0;localRow<localSize;++localRow)
                            for(BlockJacobiDetail::RowSparse::InnerIterator e(rows,neighbours[static_cast<std::size_t>(localRow)]);e;++e){
                                const auto found=std::lower_bound(neighbours.begin(),neighbours.end(),e.col());
                                if(found!=neighbours.end()&&*found==e.col())
                                    entries.emplace_back(static_cast<int>(localRow),
                                        static_cast<int>(found-neighbours.begin()),e.value());
                            }
                        block.rhsScratch.resize(localSize);block.solutionScratch.resize(localSize);
                    }else{
                        for(Eigen::Index row=block.begin;row<end;++row)
                            for(BlockJacobiDetail::RowSparse::InnerIterator e(rows,row);e;++e)
                                if(e.col()>=block.begin&&e.col()<end)
                                    entries.emplace_back(static_cast<int>(row-block.begin),
                                        static_cast<int>(e.col()-block.begin),e.value());
                    }
                    Sparse local(localSize,localSize);
                    local.setFromTriplets(entries.begin(),entries.end());local.makeCompressed();
                    block.factor=std::make_unique<Factor>();
                    block.factor->compute(local);
                    return block.factor->info()==Eigen::Success;
                },observed);
            maximumWorkers_.store(observed,std::memory_order_relaxed);
            ready_=success;
            if(!success)blocks_.clear();
            return success;
        }catch(...){blocks_.clear();return false;}
    }

    // rhs/out are caller-owned and distinct. The factor owner serializes solve
    // calls; every parallel block reuses only its own private RAS scratch.
    // On failure out is not an accepted solution.
    bool solve(const Eigen::Ref<const Eigen::VectorXd>& rhs,Eigen::Ref<Eigen::VectorXd> out) const {
        const BlockJacobiDetail::FloatingState restore;
        if(!ready_||rhs.size()!=size_||out.size()!=size_||!rhs.allFinite()
           ||(overlap_&&rhs.data()==out.data()))return false;
        if(size_==0)return true;
        int observed=0;
        const auto count=static_cast<Eigen::Index>(blocks_.size());
        const bool success=BlockJacobiDetail::parallelJobs(count,BlockJacobiDetail::workersFor(count),
            [&](Eigen::Index id){
                const auto& block=blocks_[static_cast<std::size_t>(id)];
                if(overlap_){
                    for(Eigen::Index i=0;i<block.rhsScratch.size();++i)
                        block.rhsScratch[i]=rhs[block.neighbours[static_cast<std::size_t>(i)]];
                    block.solutionScratch=block.factor->solve(block.rhsScratch);
                    if(block.factor->info()!=Eigen::Success||!block.solutionScratch.allFinite())return false;
                    for(Eigen::Index i=0;i<block.size;++i)
                        out[block.begin+i]=block.solutionScratch[block.corePositions[static_cast<std::size_t>(i)]];
                    return true;
                }
                out.segment(block.begin,block.size)=block.factor->solve(rhs.segment(block.begin,block.size));
                return block.factor->info()==Eigen::Success&&out.segment(block.begin,block.size).allFinite();
            },observed);
        BlockJacobiDetail::recordWorkers(maximumWorkers_,observed);
        return success;
    }

    int maximumWorkers() const noexcept {return maximumWorkers_.load(std::memory_order_relaxed);}
    int blockCount() const noexcept {return static_cast<int>(blocks_.size());}
};

// Restarted right-preconditioned GMRES with two MGS passes. Diagonal row
// equilibration changes only the Krylov coordinates; x and the final residual
// always refer to the original matrix. No failure selects a direct solver.
inline bool solveGmres(const Sparse& matrix,const ParallelBlockJacobi& preconditioner,
                       const Eigen::VectorXd& rhs,Eigen::VectorXd& x,double tolerance,
                       int maxIterations,int restart,GmresReport& report) {
    const BlockJacobiDetail::FloatingState restore;
    report=GmresReport{};
    const auto fail=[&](const char* message){report.message=message;return false;};
    const Eigen::Index n=matrix.rows();
    if(matrix.cols()!=n||rhs.size()!=n||!rhs.allFinite()||!(tolerance>0)||!std::isfinite(tolerance)
       ||maxIterations<=0||restart<=0)return fail("invalid GMRES input");
    try {
        for(Eigen::Index column=0;column<matrix.outerSize();++column)
            for(Sparse::InnerIterator e(matrix,column);e;++e)
                if(!std::isfinite(e.value()))return fail("nonfinite GMRES matrix");
        x=Eigen::VectorXd::Zero(n);
        const double rhsNorm=rhs.stableNorm();
        if(!std::isfinite(rhsNorm))return fail("nonfinite GMRES RHS norm");
        if(rhsNorm==0){report.message="zero RHS";return rhs.isZero(0);}
        const BlockJacobiDetail::RowSparse rows(matrix);
        const Eigen::Index width=std::min<Eigen::Index>(n,std::min(restart,maxIterations));
        Eigen::VectorXd scale(n),inverseScale(n);
        for(Eigen::Index i=0;i<n;++i){
            const double diagonal=std::abs(matrix.coeff(i,i));
            inverseScale[i]=diagonal>0?std::sqrt(diagonal):1;
            scale[i]=1/inverseScale[i];
        }
        if(!scale.allFinite()||!inverseScale.allFinite())return fail("nonfinite GMRES row scale");
        const Eigen::VectorXd scaledRhs=(scale.array()*rhs.array()).matrix();
        const double scaledRhsNorm=scaledRhs.stableNorm();
        if(!scaledRhs.allFinite()||!(scaledRhsNorm>0)||!std::isfinite(scaledRhsNorm))
            return fail("nonfinite GMRES scaled RHS");
        // ||r|| <= max(R^-1) ||R r||. Use the original RHS scale for
        // Arnoldi's early exit, with internal margin for its residual estimate.
        // A scaled-relative test can otherwise restart after every single step
        // while the actual unscaled residual still fails the original gate.
        const double estimatedAbsoluteTarget=(0.5*tolerance)*(rhsNorm/inverseScale.maxCoeff());
        Eigen::MatrixXd v(n,width+1),z(n,width),h(width+1,width);
        Eigen::VectorXd cosine(width),sine(width),g(width+1),work(n),physicalBasis(n);
        Eigen::VectorXd residual=rhs,scaledResidual=scaledRhs;
        int cycle=0,noProgress=0;
        while(report.iterations<maxIterations){
            const double beta=scaledResidual.stableNorm();
            report.relativeResidual=residual.stableNorm()/rhsNorm;
            if(!std::isfinite(beta)||!std::isfinite(report.relativeResidual))
                return fail("nonfinite GMRES residual");
            if(report.relativeResidual<=tolerance){report.message="converged";return true;}
            if(!(beta>0))return fail("GMRES scaled residual underflow");
            if(cycle++>0)++report.restarts;
            v.col(0)=scaledResidual/beta;h.setZero();g.setZero();g[0]=beta;
            Eigen::Index used=0;
            for(Eigen::Index j=0;j<width&&report.iterations<maxIterations;++j){
                physicalBasis=(inverseScale.array()*v.col(j).array()).matrix();
                if(!physicalBasis.allFinite()||!preconditioner.solve(physicalBasis,z.col(j)))
                    return fail("GMRES block Jacobi application failed");
                report.maximumWorkers=std::max(report.maximumWorkers,preconditioner.maximumWorkers());
                int observed=0;
                if(!BlockJacobiDetail::multiply(rows,z.col(j),work,observed))
                    return fail("nonfinite GMRES matrix product");
                report.maximumWorkers=std::max(report.maximumWorkers,observed);
                work.array()*=scale.array();
                if(!work.allFinite())return fail("nonfinite GMRES scaled matrix product");
                for(int pass=0;pass<2;++pass)for(Eigen::Index i=0;i<=j;++i){
                    const double projection=v.col(i).dot(work);
                    if(!std::isfinite(projection))return fail("nonfinite GMRES orthogonalization");
                    h(i,j)+=projection;work.noalias()-=projection*v.col(i);
                }
                const double next=work.stableNorm();
                if(!std::isfinite(next))return fail("nonfinite GMRES Arnoldi norm");
                h(j+1,j)=next;
                if(next>0)v.col(j+1)=work/next;
                for(Eigen::Index i=0;i<j;++i){
                    const double upper=cosine[i]*h(i,j)+sine[i]*h(i+1,j);
                    h(i+1,j)=-sine[i]*h(i,j)+cosine[i]*h(i+1,j);h(i,j)=upper;
                }
                const double norm=std::hypot(h(j,j),h(j+1,j));
                if(!(norm>0)||!std::isfinite(norm))return fail("GMRES Arnoldi breakdown");
                cosine[j]=h(j,j)/norm;sine[j]=h(j+1,j)/norm;
                h(j,j)=norm;h(j+1,j)=0;
                g[j+1]=-sine[j]*g[j];g[j]=cosine[j]*g[j];
                ++report.iterations;used=j+1;
                if(!std::isfinite(g[j+1]))return fail("nonfinite GMRES residual estimate");
                if(next==0||std::abs(g[j+1])<=estimatedAbsoluteTarget)break;
            }
            if(used==0)return fail("GMRES made no Krylov step");
            const Eigen::VectorXd y=h.topLeftCorner(used,used).template triangularView<Eigen::Upper>().solve(g.head(used));
            if(!y.allFinite())return fail("nonfinite GMRES Hessenberg solution");
            x.noalias()+=z.leftCols(used)*y;
            if(!x.allFinite())return fail("nonfinite GMRES solution");
            // Final acceptance uses the actual original matrix, not an
            // Arnoldi estimate or the block approximation.
            residual=rhs-matrix*x;
            scaledResidual=(scale.array()*residual.array()).matrix();
            if(!residual.allFinite()||!scaledResidual.allFinite())return fail("nonfinite GMRES true residual");
            report.relativeResidual=residual.stableNorm()/rhsNorm;
            if(std::isfinite(report.relativeResidual)&&report.relativeResidual<=tolerance){
                report.message="converged";return true;
            }
            noProgress=scaledResidual.stableNorm()>=beta?noProgress+1:0;
            if(noProgress>=3)return fail("GMRES stagnated");
        }
        return fail("GMRES iteration limit");
    }catch(...){return fail("GMRES numerical workspace failure");}
}

} // namespace LampLinear
