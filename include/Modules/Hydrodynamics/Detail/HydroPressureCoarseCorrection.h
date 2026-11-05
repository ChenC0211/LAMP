#pragma once

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/SparseCore>
#include <algorithm>
#include <cerrno>
#include <cfenv>
#include <cmath>
#include <cstddef>
#include <exception>
#include <limits>
#include <vector>
#if defined(_MSC_VER) && defined(_M_X64)
#include <float.h>
#endif
#if defined(_M_X64) || defined(__SSE2__)
#include <xmmintrin.h>
#endif

namespace LAMP::Hydrodynamics::Detail
{
// Current-factor pressure-only action. The caller serializes mutable scratch.
// No physical operator, current RHS, quotient/gauge or boundary is modified.
class HydroPressureCoarseCorrection
{
public:
    using Sparse=Eigen::SparseMatrix<double>;
    enum class Reason : int { NotPrepared=0, TooFewCores=1, UnsupportedEnvironment=2,
        InvalidMatrix=3, StorageLimit=4, Rank=5, InverseCertificate=6,
        ConditionOrFormation=7, Workspace=8, Prepared=9 };
    struct Report {
        Reason reason=Reason::NotPrepared;
        int dimension=0;
        std::size_t ownedPayloadBudget=0;
        double inverseResidualUpper=0,inverseNormUpper=0,conditionUpper=0;
        double formationNormUpper=0,coarseQualityUpper=0;
    };
private:
    static constexpr int coreRows_=256;
    static constexpr std::size_t storageLimit_=64u*1024u*1024u;
    struct FloatingState {
        std::fenv_t state{};int error=errno;
#if defined(_M_X64) || defined(__SSE2__)
        unsigned csr=_mm_getcsr();
#endif
        bool valid=std::fegetenv(&state)==0;
        ~FloatingState() noexcept {
            if(valid&&std::fesetenv(&state)!=0)std::terminate();
#if defined(_M_X64) || defined(__SSE2__)
            _mm_setcsr(csr);
#endif
            errno=error;
        }
    };
    Sparse matrix_;
    Eigen::MatrixXd coarse_;
    Eigen::FullPivLU<Eigen::MatrixXd> factor_;
    mutable Eigen::VectorXd coarseRhs_,coarseSolution_,permuted_;
    mutable Eigen::VectorXd initial_,remainder_,local_,action_;
    Report report_;
    bool ready_=false;
    mutable int failureStage_=0,failureRow_=-1;
    mutable double failureDefect_=0,failureAllowed_=0;
    static double up(double x) noexcept {
        return std::nextafter(x,std::numeric_limits<double>::infinity());
    }
    static double addUp(double a,double b) noexcept {return up(a+b);}
    static double multiplyUp(double a,double b) noexcept {return up(a*b);}
    static double divideUp(double a,double b) noexcept {return up(a/b);}
    static double denominatorDown(double x) noexcept {
        return std::nextafter(1-x,-std::numeric_limits<double>::infinity());
    }
    static double gamma(std::size_t operations) noexcept {
        if(operations==0)return 0;
        const double e=multiplyUp(static_cast<double>(operations),std::numeric_limits<double>::epsilon());
        const double d=denominatorDown(e);
        return d>0?divideUp(e,d):std::numeric_limits<double>::infinity();
    }
    // The absolute gradual-underflow term also propagates through later
    // rounded operations: m*q/(1-m*eps), all bounds rounded outwards.
    static double underflowBound(std::size_t operations) noexcept {
        if(operations==0)return 0;
        const double e=multiplyUp(static_cast<double>(operations),std::numeric_limits<double>::epsilon());
        const double d=denominatorDown(e);
        return d>0?divideUp(multiplyUp(static_cast<double>(operations),
            std::numeric_limits<double>::denorm_min()),d):std::numeric_limits<double>::infinity();
    }
    static bool environment() noexcept {
        if(std::numeric_limits<double>::radix!=2||std::numeric_limits<double>::digits!=53
            ||std::numeric_limits<double>::has_denorm!=std::denorm_present
            ||std::fegetround()!=FE_TONEAREST)return false;
#if defined(__FAST_MATH__) || (defined(_M_FP_FAST) && _M_FP_FAST)
        return false;
#endif
#if defined(_M_X64) || defined(__SSE2__)
        const unsigned csr=_mm_getcsr();
        if((csr&0x8040u)!=0||(csr&0x1f80u)!=0x1f80u)return false;
#endif
#if defined(_MSC_VER) && defined(_M_X64)
        unsigned control=0;
        if(_controlfp_s(&control,0,0)!=0||(control&_MCW_RC)!=_RC_NEAR
            ||(control&_MCW_EM)!=_MCW_EM)return false;
#endif
        volatile double normal=std::numeric_limits<double>::min(),two=2;
        const double half=normal/two;
        return half>0&&std::fpclassify(half)==FP_SUBNORMAL;
    }
    bool reject(Reason reason) noexcept {ready_=false;report_.reason=reason;return false;}
    static double normUpper(const Eigen::MatrixXd& matrix) noexcept {
        double norm=0;
        for(Eigen::Index i=0;i<matrix.rows();++i){
            double sum=0;
            for(Eigen::Index j=0;j<matrix.cols();++j)sum=addUp(sum,std::abs(matrix(i,j)));
            norm=std::max(norm,sum);
        }
        return norm;
    }
    // Bounds the exact-real residual of a stored-double dot product, including
    // the actual rounded product/add/subtract and gradual-underflow error.
    static double residualUpper(double rhs,double dot,double magnitude,std::size_t length) noexcept {
        const auto operations=2*length+2;
        const double scale=addUp(std::abs(rhs),magnitude);
        return addUp(std::abs(rhs-dot),addUp(multiplyUp(gamma(operations),scale),
            underflowBound(operations)));
    }
    static double inverseDefectUpper(const Eigen::MatrixXd& a,const Eigen::MatrixXd& b) noexcept {
        double norm=0;
        for(Eigen::Index i=0;i<a.rows();++i){
            double row=0;
            for(Eigen::Index j=0;j<b.cols();++j){
                double dot=0,magnitude=0;
                for(Eigen::Index k=0;k<a.cols();++k){
                    const double product=a(i,k)*b(k,j);dot+=product;
                    magnitude=addUp(magnitude,multiplyUp(std::abs(a(i,k)),std::abs(b(k,j))));
                }
                const double bound=residualUpper(i==j?1.:0.,dot,magnitude,static_cast<std::size_t>(a.cols()));
                if(!std::isfinite(bound))return std::numeric_limits<double>::infinity();
                row=addUp(row,bound);
            }
            norm=std::max(norm,row);
        }
        return norm;
    }
    // Same full-rank pivoted LU as Eigen's solve, with caller-factor scratch.
    // PAQ=LU: P*c -> unit-L solve -> U solve -> Q scatter.
    bool rawSolve(const Eigen::Ref<const Eigen::VectorXd>& rhs,Eigen::Ref<Eigen::VectorXd> out) const {
        permuted_=factor_.permutationP()*rhs;
        factor_.matrixLU().template triangularView<Eigen::UnitLower>().solveInPlace(permuted_);
        factor_.matrixLU().template triangularView<Eigen::Upper>().solveInPlace(permuted_);
        for(Eigen::Index i=0;i<coarse_.rows();++i)out[factor_.permutationQ().indices()[i]]=permuted_[i];
        return out.allFinite();
    }
    bool solveCoarse(const Eigen::Ref<const Eigen::VectorXd>& rhs,Eigen::Ref<Eigen::VectorXd> out) const {
        if(!rhs.allFinite()||!rawSolve(rhs,out))return false;
        // The exact zero equation has the exact zero solution; its relative
        // forward error is undefined, so handle this identity separately.
        if(rhs.isZero(0)&&out.isZero(0))return true;
        const auto failed=[&](double errorUpper,double allowedLower){
            failureRow_=-2;failureDefect_=errorUpper;failureAllowed_=allowedLower;return false;
        };
        double solutionNorm=0;
        for(Eigen::Index i=0;i<out.size();++i)solutionNorm=std::max(solutionNorm,std::abs(out[i]));
        const double target=std::sqrt(std::numeric_limits<double>::epsilon());
        // E/(X-E)<=target follows from E<=target*X/(1+target).
        // The threshold is rounded down; all error estimates are rounded up.
        const double fractionLower=std::nextafter(target/addUp(1.,target),0.);
        const double allowedLower=std::nextafter(fractionLower*solutionNorm,0.);
        if(!(solutionNorm>0)||!std::isfinite(allowedLower))
            return failed(std::numeric_limits<double>::infinity(),allowedLower);
        const auto length=static_cast<std::size_t>(coarse_.rows());
        double residualNormUpper=0;
        for(Eigen::Index i=0;i<coarse_.rows();++i){
            double dot=0,magnitude=0;
            for(Eigen::Index j=0;j<coarse_.cols();++j){
                const double product=coarse_(i,j)*out[j];dot+=product;
                magnitude=addUp(magnitude,multiplyUp(std::abs(coarse_(i,j)),std::abs(out[j])));
            }
            const double defect=residualUpper(rhs[i],dot,magnitude,length);
            if(!std::isfinite(defect))
                return failed(std::numeric_limits<double>::infinity(),allowedLower);
            residualNormUpper=std::max(residualNormUpper,defect);
        }
        // K bounds ||G_stored^-1||_inf and F bounds ||G_stored-G_ideal||_inf.
        // With zeta=K*F<1, K/(1-zeta) bounds ||G_ideal^-1||_inf.
        // This certifies the ideal coarse equation for the CURRENT stored rhs;
        // no exact-real restriction identity is inferred from rounded sums.
        const double formation=report_.formationNormUpper,inverse=report_.inverseNormUpper;
        const double zeta=multiplyUp(inverse,formation),denominator=denominatorDown(zeta);
        if(!std::isfinite(zeta)||!(denominator>0))
            return failed(std::numeric_limits<double>::infinity(),allowedLower);
        const double idealInverseUpper=divideUp(inverse,denominator);
        const double idealResidualUpper=addUp(residualNormUpper,multiplyUp(formation,solutionNorm));
        const double errorUpper=multiplyUp(idealInverseUpper,idealResidualUpper);
        if(!std::isfinite(errorUpper)||errorUpper>allowedLower)return failed(errorUpper,allowedLower);
        return true;
    }
    void restrictCore(const Eigen::Ref<const Eigen::VectorXd>& values) const {
        coarseRhs_.setZero();
        for(Eigen::Index i=0;i<values.size();++i)coarseRhs_[i/coreRows_]+=values[i];
    }
public:
    // False is a preparation result only: caller fixes the original RAS policy
    // BEFORE any GMRES iteration. False from apply never permits such a switch.
    bool compute(const Sparse& matrix) {
        FloatingState restore;ready_=false;report_=Report{};
        if(!restore.valid||!environment())return reject(Reason::UnsupportedEnvironment);
        if(matrix.rows()!=matrix.cols()||matrix.rows()<0
            ||matrix.rows()>std::numeric_limits<int>::max())return reject(Reason::InvalidMatrix);
        const Eigen::Index n=matrix.rows(),k=n?((n-1)/coreRows_+1):0;
        report_.dimension=static_cast<int>(k);
        if(k<2)return reject(Reason::TooFewCores);
        // Bound explicit retained/qualification payload before dense allocation.
        // This is an owner payload cap, not a vendor allocator/RSS measurement.
        const auto nn=static_cast<std::size_t>(n),kk=static_cast<std::size_t>(k);
        const auto nz=static_cast<std::size_t>(matrix.nonZeros());
        const auto bytes=nz*(sizeof(double)+sizeof(int))+(nn+1)*sizeof(int)
            +4*nn*sizeof(double)+8*kk*kk*sizeof(double)+6*kk*sizeof(double);
        report_.ownedPayloadBudget=bytes;
        if(bytes>storageLimit_)return reject(Reason::StorageLimit);
        try {
            for(Eigen::Index col=0;col<matrix.outerSize();++col)
                for(Sparse::InnerIterator e(matrix,col);e;++e)
                    if(!std::isfinite(e.value()))return reject(Reason::InvalidMatrix);
            matrix_=matrix;matrix_.makeCompressed();
            coarse_=Eigen::MatrixXd::Zero(k,k);
            Eigen::MatrixXd absolute=Eigen::MatrixXd::Zero(k,k);
            std::vector<std::size_t> counts(kk*kk,0);
            for(Eigen::Index col=0;col<matrix_.outerSize();++col)
                for(Sparse::InnerIterator e(matrix_,col);e;++e){
                    const Eigen::Index i=e.row()/coreRows_,j=col/coreRows_;
                    coarse_(i,j)+=e.value();
                    absolute(i,j)=addUp(absolute(i,j),std::abs(e.value()));
                    ++counts[static_cast<std::size_t>(i)*kk+static_cast<std::size_t>(j)];
                }
            if(!coarse_.allFinite()||!absolute.allFinite())return reject(Reason::InvalidMatrix);
            double formation=0;
            for(Eigen::Index i=0;i<k;++i){
                double row=0;
                for(Eigen::Index j=0;j<k;++j){
                    const auto count=counts[static_cast<std::size_t>(i)*kk+static_cast<std::size_t>(j)];
                    row=addUp(row,addUp(multiplyUp(gamma(count),absolute(i,j)),underflowBound(count)));
                }
                formation=std::max(formation,row);
            }
            report_.formationNormUpper=formation;
            factor_.setThreshold(gamma(kk));
            factor_.compute(coarse_);
            if(!factor_.matrixLU().allFinite())return reject(Reason::InverseCertificate);
            if(factor_.rank()!=k)return reject(Reason::Rank);
            coarseRhs_.resize(k);coarseSolution_.resize(k);permuted_.resize(k);
            Eigen::MatrixXd inverse(k,k);
            for(Eigen::Index j=0;j<k;++j){
                coarseRhs_.setZero();coarseRhs_[j]=1;
                if(!rawSolve(coarseRhs_,inverse.col(j)))return reject(Reason::InverseCertificate);
            }
            const double eta=std::max(inverseDefectUpper(coarse_,inverse),inverseDefectUpper(inverse,coarse_));
            const double denominator=denominatorDown(eta);
            if(!std::isfinite(eta)||!(denominator>0))return reject(Reason::InverseCertificate);
            const double inverseNorm=divideUp(normUpper(inverse),denominator);
            const double norm=normUpper(coarse_);
            const double condition=multiplyUp(norm,inverseNorm);
            const double quality=multiplyUp(inverseNorm,
                addUp(multiplyUp(gamma(8*kk+8),norm),formation));
            report_.inverseResidualUpper=eta;report_.inverseNormUpper=inverseNorm;
            report_.conditionUpper=condition;report_.coarseQualityUpper=quality;
            // For normal-range RHS, theta <= sqrt(eps)/(2+sqrt(eps)) is the
            // admission counterpart of 2*theta/(1-theta)<=sqrt(eps).
            // Tiny RHS retain the separate absolute-quantum residual bound.
            const double target=std::sqrt(std::numeric_limits<double>::epsilon());
            const double admission=std::nextafter(target/(2+target),0.);
            if(!std::isfinite(formation)||!std::isfinite(inverseNorm)||!std::isfinite(condition)
                ||!std::isfinite(quality)||!(norm>0)||quality>admission)
                return reject(Reason::ConditionOrFormation);
            initial_.resize(n);remainder_.resize(n);local_.resize(n);action_.resize(n);
            ready_=true;report_.reason=Reason::Prepared;return true;
        }catch(...){return reject(Reason::Workspace);}
    }
    template<class LocalSolve>
    bool apply(const Eigen::Ref<const Eigen::VectorXd>& rhs,Eigen::Ref<Eigen::VectorXd> out,
        const LocalSolve& localSolve) const {
        FloatingState restore;failureStage_=1;failureRow_=-1;failureDefect_=failureAllowed_=0;
        if(!ready_||!restore.valid||!environment()||rhs.size()!=matrix_.rows()
            ||out.size()!=matrix_.rows()||rhs.data()==out.data()||!rhs.allFinite())return false;
        try {
            restrictCore(rhs);failureStage_=2;
            if(!solveCoarse(coarseRhs_,coarseSolution_))return false;
            for(Eigen::Index i=0;i<rhs.size();++i)initial_[i]=coarseSolution_[i/coreRows_];
            failureStage_=3;action_.noalias()=matrix_*initial_;
            remainder_=rhs-action_;
            if(!remainder_.allFinite()||!localSolve(remainder_,local_)||!local_.allFinite())return false;
            action_.noalias()=matrix_*local_;
            if(!action_.allFinite())return false;
            restrictCore(action_);failureStage_=4;
            if(!solveCoarse(coarseRhs_,coarseSolution_))return false;
            for(Eigen::Index i=0;i<rhs.size();++i)
                out[i]=(initial_[i]+local_[i])-coarseSolution_[i/coreRows_];
            failureStage_=5;const bool finite=out.allFinite();if(finite)failureStage_=0;return finite;
        }catch(...){return false;}
    }
    int failureStage() const noexcept {return failureStage_;}
    int failureRow() const noexcept {return failureRow_;}
    double failureDefect() const noexcept {return failureDefect_;}
    double failureAllowed() const noexcept {return failureAllowed_;}
    bool ready() const noexcept {return ready_;}
    const Report& report() const noexcept {return report_;}
    const Sparse& matrix() const noexcept {return matrix_;}
    const Eigen::MatrixXd& coarseMatrix() const noexcept {return coarse_;}
    static constexpr int coreRows() noexcept {return coreRows_;}
};
}
