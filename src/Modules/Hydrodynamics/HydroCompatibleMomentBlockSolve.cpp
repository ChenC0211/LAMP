#include "Modules/Hydrodynamics/detail/HydroCompatibleMomentBlockSolve.h"
#include <Eigen/SuperLUSupport>
#include <Eigen/LU>
#include <unsupported/Eigen/IterativeSolvers>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace Quantum::Hydrodynamics {class MomentBlockOperator;}
namespace Eigen::internal {
template<> struct traits<Quantum::Hydrodynamics::MomentBlockOperator>:traits<Quantum::Hydrodynamics::HydroCompatibleMomentBlockSolve::Matrix>{};
}
namespace Quantum::Hydrodynamics {
// Exact cell-interior elimination in the complete velocity/boundary block.
// Local Hermite bubbles have support on one native leaf. Their response is
// retained in the Schur complement and recovered on EVERY solve, including
// applications of a retained factor as a right preconditioner.
class MomentVelocityFactor {
    using Matrix=HydroCompatibleMomentBlockSolve::Matrix;
    using Vector9=Eigen::Matrix<double,9,1>;
    struct Local {
        std::array<int,9> interior;
        std::vector<int> rows,columns;
        Eigen::PartialPivLU<Eigen::MatrixXd> diagonal;
        Eigen::MatrixXd left,inverseRight;
    };
    static void require(bool ok,const char* message){if(!ok)throw std::invalid_argument(message);}
    int size_=0;
    std::vector<int> retained_;
    std::vector<Local> locals_;
    Eigen::SuperLU<Matrix> factor_;
public:
    int retainedSize() const{return int(retained_.size());}
    void compute(const Matrix& matrix,const std::vector<std::array<int,9>>& interiors){
        const auto start=std::chrono::steady_clock::now();
        size_=int(matrix.rows());require(matrix.cols()==size_,"nonsquare complete velocity block");
        std::vector<int> owner(size_,-1),position(size_,-1),retainedIndex(size_,-1);
        locals_.clear();locals_.resize(interiors.size());retained_.clear();
        for(int k=0;k<int(interiors.size());k++)for(int j=0;j<9;j++){
            const int i=interiors[k][j];require(i>=0&&i<size_&&owner[i]<0,"invalid or repeated cell-interior coordinate");
            owner[i]=k;position[i]=j;locals_[k].interior[j]=i;
        }
        for(int i=0;i<size_;i++)if(owner[i]<0){retainedIndex[i]=int(retained_.size());retained_.push_back(i);}
        Matrix schur(retainedSize(),retainedSize());
        {
            struct Work {
                Eigen::Matrix<double,9,9> diagonal=Eigen::Matrix<double,9,9>::Zero();
                std::vector<Eigen::Triplet<double>> left,right;
            };
            std::vector<Work> work(interiors.size());
            std::vector<Eigen::Triplet<double>> entries;entries.reserve(matrix.nonZeros());
            for(int col=0;col<matrix.outerSize();col++)for(Matrix::InnerIterator it(matrix,col);it;++it){
                if(it.value()==0)continue;
                const int row=int(it.row()),a=owner[row],b=owner[col];
                if(a<0&&b<0)entries.emplace_back(retainedIndex[row],retainedIndex[col],it.value());
                else if(a>=0&&b>=0){
                    // Do not discard small cross-cell entries: that would be
                    // a different operator, rather than exact condensation.
                    require(a==b,"cell-interior coordinates couple across native leaves");
                    work[a].diagonal(position[row],position[col])=it.value();
                }else if(a<0)work[b].left.emplace_back(retainedIndex[row],position[col],it.value());
                else work[a].right.emplace_back(position[row],retainedIndex[col],it.value());
            }
            for(int k=0;k<int(locals_.size());k++){
                auto& local=locals_[k];const auto& w=work[k];
                for(const auto& v:w.left)local.rows.push_back(v.row());
                for(const auto& v:w.right)local.columns.push_back(v.col());
                const auto unique=[](std::vector<int>& v){std::sort(v.begin(),v.end());v.erase(std::unique(v.begin(),v.end()),v.end());};
                unique(local.rows);unique(local.columns);
                local.left=Eigen::MatrixXd::Zero(local.rows.size(),9);
                Eigen::MatrixXd right=Eigen::MatrixXd::Zero(9,local.columns.size());
                for(const auto& v:w.left)local.left(std::lower_bound(local.rows.begin(),local.rows.end(),v.row())-local.rows.begin(),v.col())=v.value();
                for(const auto& v:w.right)right(v.row(),std::lower_bound(local.columns.begin(),local.columns.end(),v.col())-local.columns.begin())=v.value();
                local.diagonal.compute(w.diagonal);
                require(local.diagonal.matrixLU().allFinite()&&(local.diagonal.matrixLU().diagonal().array()!=0).all(),"singular cell-interior velocity block");
                local.inverseRight=local.diagonal.solve(right);
                require(local.inverseRight.allFinite(),"nonfinite cell-interior elimination");
                const Eigen::MatrixXd correction=local.left*local.inverseRight;
                for(int j=0;j<correction.cols();j++)for(int i=0;i<correction.rows();i++)
                    if(correction(i,j)!=0)entries.emplace_back(local.rows[i],local.columns[j],-correction(i,j));
            }
            schur.setFromTriplets(entries.begin(),entries.end());
        } // Release assembly storage before the global sparse factorization.
        const auto assembled=std::chrono::steady_clock::now();
        if(retainedSize()){factor_.compute(schur);require(factor_.info()==Eigen::Success,"condensed velocity/boundary factorization failed");}
        if(std::getenv("QUANTUM_MOMENT_STEP_PROFILE"))std::cout<<"[native_velocity_condensation] full="<<size_
            <<" retained="<<retainedSize()<<" cells="<<locals_.size()<<" input_nnz="<<matrix.nonZeros()<<" schur_nnz="<<schur.nonZeros()
            <<" assembly_seconds="<<std::chrono::duration<double>(assembled-start).count()
            <<" factor_seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-assembled).count()<<std::endl;
    }
    Eigen::VectorXd solve(const Eigen::VectorXd& rhs) const{
        require(rhs.size()==size_&&rhs.allFinite(),"invalid condensed velocity right hand side");
        Eigen::VectorXd condensed(retainedSize()),interior(9*locals_.size());
        for(int i=0;i<retainedSize();i++)condensed[i]=rhs[retained_[i]];
        for(int k=0;k<int(locals_.size());k++){
            const auto& local=locals_[k];Vector9 load;for(int i=0;i<9;i++)load[i]=rhs[local.interior[i]];
            const Vector9 response=local.diagonal.solve(load);interior.segment<9>(9*k)=response;
            for(int i=0;i<int(local.rows.size());i++)condensed[local.rows[i]]-=local.left.row(i).dot(response);
        }
        Eigen::VectorXd trace(retainedSize());
        if(retainedSize()){trace=factor_.solve(condensed);require(factor_.info()==Eigen::Success&&trace.allFinite(),"condensed velocity solve failed");}
        Eigen::VectorXd result(size_);
        for(int i=0;i<retainedSize();i++)result[retained_[i]]=trace[i];
        for(int k=0;k<int(locals_.size());k++){
            const auto& local=locals_[k];Vector9 value=interior.segment<9>(9*k);
            for(int j=0;j<int(local.columns.size());j++)value-=local.inverseRight.col(j)*trace[local.columns[j]];
            for(int i=0;i<9;i++)result[local.interior[i]]=value[i];
        }
        require(result.allFinite(),"nonfinite recovered complete velocity");return result;
    }
};
// Matrix-free A P^{-1}; right preconditioning keeps the Krylov residual in
// the original momentum/boundary equations, with their original scaling.
class MomentBlockOperator:public Eigen::EigenBase<MomentBlockOperator> {
public:
    using Scalar=double;using RealScalar=double;using StorageIndex=int;
    using Matrix=HydroCompatibleMomentBlockSolve::Matrix;
    enum {ColsAtCompileTime=Eigen::Dynamic,MaxColsAtCompileTime=Eigen::Dynamic,IsRowMajor=false};
    MomentBlockOperator(const Matrix& matrix,MomentVelocityFactor& factor):matrix_(matrix),factor_(factor){}
    Eigen::Index rows() const{return matrix_.rows();}
    Eigen::Index cols() const{return matrix_.cols();}
    template<class Rhs> Eigen::Product<MomentBlockOperator,Rhs,Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs>& rhs) const{
        return Eigen::Product<MomentBlockOperator,Rhs,Eigen::AliasFreeProduct>(*this,rhs.derived());
    }
    template<class Rhs> Eigen::VectorXd apply(const Eigen::MatrixBase<Rhs>& rhs) const{
        const Eigen::VectorXd velocity=factor_.solve(rhs);
        return matrix_*velocity;
    }
private:
    const Matrix& matrix_;
    MomentVelocityFactor& factor_;
};
}
namespace Eigen::internal {
template<class Rhs> struct generic_product_impl<Quantum::Hydrodynamics::MomentBlockOperator,Rhs,SparseShape,DenseShape,GemvProduct>
    :generic_product_impl_base<Quantum::Hydrodynamics::MomentBlockOperator,Rhs,generic_product_impl<Quantum::Hydrodynamics::MomentBlockOperator,Rhs>> {
    using Operator=Quantum::Hydrodynamics::MomentBlockOperator;
    template<class Dest> static void scaleAndAddTo(Dest& result,const Operator& matrix,const Rhs& rhs,const double& alpha){
        result+=alpha*matrix.apply(rhs);
    }
};
}
namespace Quantum::Hydrodynamics {
namespace {
void require(bool ok,const char* message){if(!ok)throw std::invalid_argument(message);}
using Matrix=HydroCompatibleMomentBlockSolve::Matrix;
bool identical(const Matrix& a,const Matrix& b){
    if(a.rows()!=b.rows()||a.cols()!=b.cols()||a.nonZeros()!=b.nonZeros())return false;
    for(int col=0;col<a.outerSize();col++){
        Matrix::InnerIterator i(a,col),j(b,col);
        for(;i&&j;++i,++j)if(i.row()!=j.row()||i.value()!=j.value())return false;
        if(i||j)return false;
    }
    return true;
}
}
struct HydroCompatibleMomentBlockSolve::Impl {
    MomentVelocityFactor velocity;
    Eigen::SuperLU<Matrix> pressure;
    std::vector<std::array<int,9>> cellInteriors;
    Matrix coupling,velocityMatrix,pressureMatrix;
    int nq=0,ne=0,np=0,nb=0,ng=0,total=0;
    bool ready=false,velocityCurrent=false;
    Statistics statistics;
    void factorVelocity(){
        velocity.compute(velocityMatrix,cellInteriors);statistics.retainedVelocityUnknowns=velocity.retainedSize();
        velocityCurrent=true;++statistics.velocityFactorizations;
    }
};
HydroCompatibleMomentBlockSolve::HydroCompatibleMomentBlockSolve():impl_(std::make_unique<Impl>()){}
HydroCompatibleMomentBlockSolve::~HydroCompatibleMomentBlockSolve()=default;
void HydroCompatibleMomentBlockSolve::compute(const Matrix& system,int nq,int ne,int np,int nb,int ng,const std::vector<std::array<int,9>>& cellInteriors){
    auto& p=*impl_;const bool compatible=p.ready&&p.nq==nq&&p.ne==ne&&p.np==np&&p.nb==nb&&p.ng==ng&&p.cellInteriors==cellInteriors;p.ready=false;
    require(nq>0&&np>0&&nb>=0&&(ng==0||ng==1)&&ne+ng==np,"invalid complete saddle block dimensions");
    const int nv=nq+ne,total=nv+np+nb+ng,boundaryBegin=nv+np;
    require(system.rows()==total&&system.cols()==total,"complete saddle matrix has inconsistent dimensions");
    for(const auto& cell:cellInteriors)for(int i:cell)require(i>=0&&i<nq,"cell-interior coordinate outside solenoidal velocity");
    p.cellInteriors=cellInteriors;
    p.nq=nq;p.ne=ne;p.np=np;p.nb=nb;p.ng=ng;p.total=total;
    const auto reduced=[&](int i){return i<nq?i:(i>=boundaryBegin&&i<boundaryBegin+nb?nq+i-boundaryBegin:-1);};
    Matrix v(nq+nb,nq+nb),c(np,np);
    {
    std::vector<Eigen::Triplet<double>> velocity,pressure,coupling;
    for(int col=0;col<system.outerSize();col++){
        const int j=reduced(col);
        for(Matrix::InnerIterator it(system,col);it;++it){
            const int row=int(it.row()),i=reduced(row);
            if(i>=0&&j>=0)velocity.emplace_back(i,j,it.value());
            if(row>=nq&&row<nv&&col<nq)coupling.emplace_back(row-nq,col,it.value());
            if(col>=nv&&col<nv+np){
                if(row>=nq&&row<nv)pressure.emplace_back(row-nq,col-nv,it.value());
                else if(ng&&row==total-1)pressure.emplace_back(ne,col-nv,it.value());
            }
        }
    }
    v.setFromTriplets(velocity.begin(),velocity.end());c.setFromTriplets(pressure.begin(),pressure.end());
    p.coupling.resize(ne,nq);p.coupling.setFromTriplets(coupling.begin(),coupling.end());
    }
    // The pressure system consists of -G^T and, for a closed domain, the
    // ORIGINAL physical mean-pressure gauge. No cell pressure is pinned.
    if(!compatible||!identical(c,p.pressureMatrix)){
        p.pressure.compute(c);require(p.pressure.info()==Eigen::Success,"complete pressure block factorization failed");
        ++p.statistics.pressureFactorizations;p.pressureMatrix=std::move(c);
    }
    if(!compatible||!identical(v,p.velocityMatrix))p.velocityCurrent=false;
    p.velocityMatrix=std::move(v);
    if(!compatible)p.factorVelocity();
    p.ready=true;
}
const HydroCompatibleMomentBlockSolve::Statistics& HydroCompatibleMomentBlockSolve::statistics() const noexcept{return impl_->statistics;}
Eigen::VectorXd HydroCompatibleMomentBlockSolve::solve(const Eigen::VectorXd& rhs,double budget){
    auto& p=*impl_;require(p.ready&&rhs.size()==p.total&&rhs.allFinite(),"invalid complete block-solve right hand side");
    require(std::isfinite(budget)&&budget>=0,"invalid block-solve residual budget");p.statistics.krylovIterations=0;
    const int nv=p.nq+p.ne,boundaryBegin=nv+p.np;
    require(rhs.segment(nv,p.np).isZero(0),"complete solenoidal block solve requires zero divergence source");
    Eigen::VectorXd reduced(p.nq+p.nb);reduced.head(p.nq)=rhs.head(p.nq);reduced.tail(p.nb)=rhs.segment(boundaryBegin,p.nb);
    if(!p.velocityCurrent&&budget==0)p.factorVelocity();
    Eigen::VectorXd velocity=p.velocity.solve(reduced);
    if(!p.velocityCurrent){
        const Eigen::VectorXd residual=reduced-p.velocityMatrix*velocity;
        if(residual.stableNorm()>budget){
            const MomentBlockOperator current(p.velocityMatrix,p.velocity);
            Eigen::GMRES<MomentBlockOperator,Eigen::IdentityPreconditioner> krylov;
            krylov.setMaxIterations(24);krylov.set_restart(std::min<Eigen::Index>(24,reduced.size()));
            krylov.setTolerance(.5*budget/residual.stableNorm());
            krylov.compute(current);
            const Eigen::VectorXd correction=krylov.solveWithGuess(residual,Eigen::VectorXd::Zero(reduced.size()));
            velocity+=p.velocity.solve(correction);
            p.statistics.krylovIterations=int(krylov.iterations());
        }
        // A preconditioned stopping estimate is not the acceptance criterion.
        // If the old factor is no longer useful, refresh it for this SAME
        // equation; never accept a preceding Picard matrix or change dt.
        if(!velocity.allFinite()||(p.velocityMatrix*velocity-reduced).stableNorm()>budget){
            p.factorVelocity();velocity=p.velocity.solve(reduced);
        }
    }
    Eigen::VectorXd result=Eigen::VectorXd::Zero(p.total);result.head(p.nq)=velocity.head(p.nq);result.segment(boundaryBegin,p.nb)=velocity.tail(p.nb);
    Eigen::VectorXd pressureRhs(p.np);pressureRhs.head(p.ne)=rhs.segment(p.nq,p.ne)-p.coupling*result.head(p.nq);
    if(p.ng)pressureRhs[p.ne]=rhs[p.total-1];
    result.segment(nv,p.np)=p.pressure.solve(pressureRhs);
    require(p.pressure.info()==Eigen::Success&&result.allFinite(),"complete pressure recovery failed");
    return result;
}
}
