#include "Modules/Hydrodynamics/CchSinglePhaseSaddle.h"

#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>

namespace Quantum::Hydrodynamics
{
namespace
{
using Sparse = CchSinglePhaseSaddle::Matrix;
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
    if (!ok) throw std::invalid_argument(std::string("[CchSinglePhaseSaddle] ")+reason);
}
double ratio(double error, double scale)
{
    if (!std::isfinite(error) || !std::isfinite(scale) || scale<0)
        throw std::runtime_error("[CchSinglePhaseSaddle] nonfinite residual/scale");
    if (scale==0) {
        if (error!=0) throw std::runtime_error("[CchSinglePhaseSaddle] nonzero residual at zero scale");
        return 0;
    }
    const double result=std::abs(error)/scale;
    if (!std::isfinite(result)) throw std::runtime_error("[CchSinglePhaseSaddle] nonfinite normalized residual");
    return result;
}
void accept(double defect, double tolerance, const char* stage)
{
    if (defect>tolerance) throw std::runtime_error(std::string("[CchSinglePhaseSaddle] ")+stage+" residual exceeds tolerance");
}
int root(std::vector<int>& parents, int p)
{
    while (parents[p]!=p) {parents[p]=parents[parents[p]];p=parents[p];}
    return p;
}
}
struct CchSaddleVelocityBlock::Data
{
    Data(CchDofLayout input, Matrix matrix) : layout(std::move(input)),integral(std::move(matrix)) {}
    CchDofLayout layout;
    Matrix integral;
    double symmetry=0;
    int tangential=0;
};
CchSaddleVelocityBlock::CchSaddleVelocityBlock(std::shared_ptr<const Data> data) : data_(std::move(data)) {}
CchSaddleVelocityBlock CchSaddleVelocityBlock::build(CchDofLayout layout, const Matrix& integralVelocity,int tangentialCount)
{
    require(tangentialCount>=0,"negative tangential block size");
    require(integralVelocity.rows()==layout.velocityCount()+tangentialCount && integralVelocity.cols()==integralVelocity.rows(),
        "velocity block dimensions do not match layout");
    auto data=std::make_shared<Data>(std::move(layout),integralVelocity);
    data->tangential=tangentialCount;
    data->integral.makeCompressed();
    for (Eigen::Index i=0;i<data->integral.nonZeros();++i)
        require(std::isfinite(data->integral.valuePtr()[i]),"nonfinite velocity block");
    const Matrix difference=data->integral-Matrix(data->integral.transpose());
    const double scale=data->integral.norm();
    require(std::isfinite(scale) && scale>0,"unresolved velocity block norm");
    data->symmetry=ratio(difference.norm(),scale);
    require(data->symmetry<=1e-12,"velocity block is not symmetric");
    Eigen::SimplicialLLT<Sparse> check;
    check.compute(Sparse(data->integral));
    require(check.info()==Eigen::Success,"velocity block is not positive definite");
    const Sparse factor=check.matrixL();
    for (Eigen::Index i=0;i<factor.nonZeros();++i)
        require(std::isfinite(factor.valuePtr()[i]),"nonfinite velocity block factorization");
    return CchSaddleVelocityBlock(std::move(data));
}
const CchDofLayout& CchSaddleVelocityBlock::layout() const noexcept {return data_->layout;}
const CchSaddleVelocityBlock::Matrix& CchSaddleVelocityBlock::matrix() const noexcept {return data_->integral;}
CchSaddleVelocityBlock CchSaddleVelocityBlock::withSkew(const Matrix& skew) const
{
    require(skew.rows()==matrix().rows()&&skew.cols()==matrix().cols(),"different skew block dimensions");
    for(Eigen::Index i=0;i<skew.nonZeros();++i)require(std::isfinite(skew.valuePtr()[i]),"nonfinite skew transport block");
    const double scale=skew.norm();require(std::isfinite(scale),"nonfinite skew block norm");
    const Matrix defect=skew+Matrix(skew.transpose());
    require(ratio(defect.norm(),scale)<=1e-12,"transport block is not skew; material/boundary pairing required");
    auto data=std::make_shared<Data>(*data_);data->integral+=skew;data->integral.makeCompressed();
    for(Eigen::Index i=0;i<data->integral.nonZeros();++i)require(std::isfinite(data->integral.valuePtr()[i]),"nonfinite combined velocity block");
    const Matrix difference=data->integral-Matrix(data->integral.transpose());
    data->symmetry=ratio(difference.norm(),data->integral.norm());
    return CchSaddleVelocityBlock(std::move(data));
}
double CchSaddleVelocityBlock::symmetryDefect() const noexcept {return data_->symmetry;}
int CchSaddleVelocityBlock::tangentialCount() const noexcept{return data_->tangential;}

CchSaddleRhs::CchSaddleRhs(CchDofLayout layout,const std::vector<double>& integralMomentum,
    const std::vector<double>& prescribedNormalVelocity,const std::vector<double>& tangentialMomentum)
    : layout_(std::move(layout)),momentum_(integralMomentum),prescribed_(prescribedNormalVelocity),tangential_(tangentialMomentum)
{
    require(momentum_.size()==static_cast<std::size_t>(layout_.velocityCount()) && prescribed_.size()==momentum_.size(),
        "RHS dimensions do not match layout");
    for (const auto* values : {&momentum_,&prescribed_,&tangential_})
        for (double value : *values) require(std::isfinite(value),"nonfinite RHS/prescribed velocity");
}
void CchSaddleRhs::validateFor(const CchDofLayout& expected,const QuadtreeMesh& mesh) const
{
    require(layout_.sameIdentity(expected),"different RHS layout identity");
    expected.validateMesh(mesh);
}
CchSaddleResult::CchSaddleResult(CchNativeValues values,CchSaddleReport report,std::vector<double> tangential)
    : values_(std::move(values)),report_(report),tangential_(std::move(tangential)) {}

struct CchSinglePhaseSaddle::Data
{
    Data(CchPressureLoad p,CchSaddleVelocityBlock a,CchSaddleConfig options)
        : pressure(std::move(p)),velocity(std::move(a)),config(options) {}
    CchPressureLoad pressure;
    CchSaddleVelocityBlock velocity;
    CchSaddleConfig config;
    std::vector<int> freeByVelocity,velocityByFree,component;
    int gauges=0;
    Matrix system;
    Eigen::SparseLU<Matrix> factor;
};
CchSinglePhaseSaddle::CchSinglePhaseSaddle(std::shared_ptr<const Data> data) : data_(std::move(data)) {}
CchSinglePhaseSaddle CchSinglePhaseSaddle::buildClosed(CchPressureLoad pressure,
    CchSaddleVelocityBlock velocity,CchSaddleConfig config)
{
    require(pressure.layout().sameIdentity(velocity.layout()),"different velocity block layout identity");
    for (double tolerance : {config.algebraicRelativeTolerance,config.momentumRelativeTolerance,
        config.continuityRelativeTolerance,config.pressureMeanRelativeTolerance,config.boundaryFluxRelativeTolerance})
        require(std::isfinite(tolerance) && tolerance>0 && tolerance<1,"invalid residual tolerance");
    auto data=std::make_shared<Data>(std::move(pressure),std::move(velocity),config);
    const auto& c=data->pressure.matrix(); const int nu=static_cast<int>(c.rows()),np=static_cast<int>(c.cols());
    const int nm=data->velocity.tangentialCount();data->freeByVelocity.assign(nu+nm,-1);
    std::vector<int> parents(np);std::iota(parents.begin(),parents.end(),0);
    for (int u=0;u<nu;++u)
    {
        int first=-1,last=-1,count=0;
        for (CchPressureLoad::Matrix::InnerIterator it(c,u);it;++it) {if (!count) first=static_cast<int>(it.col());last=static_cast<int>(it.col());++count;}
        require(count<=2,"invalid native pressure incidence");
        if (count==2) {
            data->freeByVelocity[u]=static_cast<int>(data->velocityByFree.size());data->velocityByFree.push_back(u);
            parents[root(parents,last)]=root(parents,first);
        }
    }
    for(int m=0;m<nm;++m){const int u=nu+m;
        data->freeByVelocity[u]=static_cast<int>(data->velocityByFree.size());data->velocityByFree.push_back(u);}
    data->component.resize(np);std::vector<int> numbered(np,-1);
    for (int p=0;p<np;++p) {
        const int r=root(parents,p);if (numbered[r]<0) numbered[r]=data->gauges++;
        data->component[p]=numbered[r];
    }
    const Eigen::Index nf=static_cast<Eigen::Index>(data->velocityByFree.size()),size=nf+np+data->gauges;
    std::vector<Eigen::Triplet<double,Eigen::Index>> entries;
    const auto& a=data->velocity.matrix();
    entries.reserve(static_cast<std::size_t>(a.nonZeros()+2*c.nonZeros()+2*np));
    for (int u : data->velocityByFree)
        for (CchSaddleVelocityBlock::Matrix::InnerIterator it(a,u);it;++it)
            if (data->freeByVelocity[it.col()]>=0) entries.emplace_back(data->freeByVelocity[u],data->freeByVelocity[it.col()],it.value());
    for (int u : data->velocityByFree)if(u<nu)
        for (CchPressureLoad::Matrix::InnerIterator it(c,u);it;++it) {
            entries.emplace_back(data->freeByVelocity[u],nf+it.col(),it.value());
            entries.emplace_back(nf+it.col(),data->freeByVelocity[u],it.value());
        }
    const auto& volumes=data->pressure.layout().domain().cellVolumes();
    for (int p=0;p<np;++p) {
        require(std::isfinite(volumes[p]) && volumes[p]>0,"invalid pressure mean volume");
        entries.emplace_back(nf+p,nf+np+data->component[p],volumes[p]);
        entries.emplace_back(nf+np+data->component[p],nf+p,volumes[p]);
    }
    data->system.resize(size,size);data->system.setFromTriplets(entries.begin(),entries.end());data->system.makeCompressed();
    data->factor.compute(data->system);
    if (data->factor.info()!=Eigen::Success) throw std::runtime_error("[CchSinglePhaseSaddle] reference saddle factorization failed");
    return CchSinglePhaseSaddle(std::move(data));
}
const CchPressureLoad& CchSinglePhaseSaddle::pressureLoad() const noexcept {return data_->pressure;}
const CchSaddleVelocityBlock& CchSinglePhaseSaddle::velocityBlock() const noexcept {return data_->velocity;}
const CchDofLayout& CchSinglePhaseSaddle::layout() const noexcept {return data_->pressure.layout();}
const CchSinglePhaseSaddle::Matrix& CchSinglePhaseSaddle::matrix() const noexcept {return data_->system;}
const std::vector<int>& CchSinglePhaseSaddle::freeVelocityDofs() const noexcept {return data_->velocityByFree;}
const std::vector<int>& CchSinglePhaseSaddle::pressureComponents() const noexcept {return data_->component;}
int CchSinglePhaseSaddle::pressureGaugeCount() const noexcept {return data_->gauges;}
CchSaddleResult CchSinglePhaseSaddle::solve(const QuadtreeMesh& mesh,const CchSaddleRhs& rhs) const
{
    rhs.validateFor(layout(),mesh);
    const int nu=layout().velocityCount(),nm=data_->velocity.tangentialCount(),nv=nu+nm,
        np=layout().pressureCount(),nf=static_cast<int>(data_->velocityByFree.size());
    const auto& c=data_->pressure.matrix(); const auto& a=data_->velocity.matrix();
    require(rhs.tangentialMomentum().size()==static_cast<std::size_t>(nm),"tangential RHS dimensions do not match velocity block");
    Eigen::VectorXd known=Eigen::VectorXd::Zero(nv),momentum(nv);
    for(int u=0;u<nu;++u){known[u]=rhs.prescribedNormalVelocity()[u];momentum[u]=rhs.integralMomentum()[u];}
    for(int m=0;m<nm;++m)momentum[nu+m]=rhs.tangentialMomentum()[m];
    std::vector<double> net(data_->gauges,0),gross(data_->gauges,0);CchSaddleReport report;
    report.freeVelocityCount=nf-nm;report.freeTangentialCount=nm;report.pressureGaugeCount=data_->gauges;
    for (int u=0;u<nu;++u) {
        if (data_->freeByVelocity[u]>=0) require(known[u]==0,"prescribed velocity payload on a free face");
        else {
            int count=0;
            for (CchPressureLoad::Matrix::InnerIterator it(c,u);it;++it) {
                const int component=data_->component[it.col()];const double flow=-it.value()*known[u];
                net[component]+=flow;gross[component]+=std::abs(flow);++count;
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
    if (!b.allFinite()) throw std::runtime_error("[CchSinglePhaseSaddle] nonfinite eliminated RHS");
    // Assess the original equations and tolerances before publishing a report.
    // No geometry, material or manufactured equilibrium enters this solve.
    struct Assessment {
        std::vector<double> u,p;
        std::array<double,4> residual{};
        double gaugeMultiplier=0;
    };
    const std::array<double,4> tolerance{{data_->config.algebraicRelativeTolerance,
        data_->config.momentumRelativeTolerance,data_->config.continuityRelativeTolerance,
        data_->config.pressureMeanRelativeTolerance}};
    const auto assess=[&](const Eigen::VectorXd& x){
        if(!x.allFinite())throw std::runtime_error("[CchSinglePhaseSaddle] nonfinite reference saddle solution");
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
            for(CchIntegralDivergence::Matrix::InnerIterator it(d,i);it;++it)flowScale+=std::abs(it.value()*velocity[it.col()]);
            value.residual[2]=std::max(value.residual[2],ratio(continuity[i],flowScale));
        }
        std::vector<double> mean(data_->gauges,0),meanScale(data_->gauges,0);
        const auto& volumes=layout().domain().cellVolumes();
        for(int i=0;i<np;++i){const double term=volumes[i]*value.p[i];mean[data_->component[i]]+=term;meanScale[data_->component[i]]+=std::abs(term);}
        for(int i=0;i<data_->gauges;++i){
            value.residual[3]=std::max(value.residual[3],ratio(mean[i],meanScale[i]));
            value.gaugeMultiplier=std::max(value.gaugeMultiplier,std::abs(x[nf+np+i]));
        }
        return value;
    };
    Eigen::VectorXd x=data_->factor.solve(b);
    for(int iteration=0;;++iteration){
        const auto value=assess(x);bool accepted=true;
        for(int k=0;k<4;++k)accepted=accepted && value.residual[k]<=tolerance[k];
        if(accepted){
            report.algebraicRelativeResidual=value.residual[0];report.freeMomentumDefect=value.residual[1];
            report.maxCellContinuityDefect=value.residual[2];report.maxPressureMeanDefect=value.residual[3];
            report.maxAbsGaugeMultiplier=value.gaugeMultiplier;report.linearRefinements=iteration;
            const std::vector<double> normal(value.u.begin(),value.u.begin()+nu),tangential(value.u.begin()+nu,value.u.end());
            return CchSaddleResult(CchNativeValues(layout(),normal,value.p),report,tangential);
        }
        if(iteration==8){
            const std::array<const char*,4> stages{{"complete algebraic","free momentum","native cell continuity","component pressure mean"}};
            for(int k=0;k<4;++k)accept(value.residual[k],tolerance[k],stages[k]);
        }
        Eigen::VectorXd residual=correctionResidual(data_->system,x,b);
        // An already accurate force balance is preserved by A*du+C*dp=0.
        // Reinjecting its roundoff can obscure relative continuity at u~0.
        // Every correction still solves the SAME full augmented matrix.
        if(value.residual[0]<=tolerance[0] && value.residual[1]<=tolerance[1])residual.head(nf).setZero();
        if(value.residual[3]<=tolerance[3])residual.tail(data_->gauges).setZero();
        x+=data_->factor.solve(residual);
    }
}
}
