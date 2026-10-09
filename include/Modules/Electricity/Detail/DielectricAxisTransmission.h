#pragma once

#include "Modules/Electricity/Detail/DielectricCubicTrace.h"
#include <Eigen/SVD>

namespace Quantum::Electricity::Detail
{
namespace
{
namespace DielectricAxisTransmission
{
using Phase=DielectricCubicTrace::Phase;
using SamplePoint=DielectricCubicTrace::SamplePoint;
using Jump=DielectricCubicTrace::Jump;
using State=DielectricCubicTrace::State;
using SegmentPiece=DielectricCubicTrace::SegmentPiece;
using LinkGeometry=DielectricCubicTrace::LinkGeometry;
inline void require(bool ok,const char* reason)
{if(!ok)throw std::runtime_error(std::string("[DielectricAxisTransmission] ")+reason);}

// Local thin-column model, not a generic replacement for P3. The caller must
// certify negative phase on the complete axis interval of its fixed window,
// and positive-radius positive-phase queries. A detached cap/gas axis is NOT
// in this model's geometric scope. Geometry, support and global FV assembly
// are deliberately not inferred by this algebraic reconstruction.
//
// Interior: four regular axisymmetric harmonics through axial degree three.
// Exterior: the same four plus their four logarithmic partners. Seven local
// curved-interface transmission constraints leave five homogeneous modes.
// Both phase polynomials solve the cylindrical PDE throughout their domain;
// continuity through d_s^3 and Dn jumps through d_s^2 hold at the anchor only.
// This finite local expansion is a model, not an accuracy claim for a neck.
class Reconstruction
{
    using Row=Eigen::Matrix<double,1,12>;
    struct Jet {double v=0,r=0,z=0,rr=0,rz=0,zz=0,rrr=0,rrz=0,rzz=0,zzz=0;};
    Vector2D anchor_,normal_;
    double h_,k_,ks_,em_,ep_,condition_=0,reproduction_=0,constraintCondition_=0;
    Eigen::MatrixXd homogeneous_,particular_,sampleParticular_,sampling_;
    static double derivative(int power,int order,double x,bool logarithm)
    {
        if(!logarithm) {
            if(order>power)return 0;
            double value=std::pow(x,power-order);
            for(int j=0;j<order;++j)value*=power-j;
            return value;
        }
        require(x>0,"exterior logarithmic harmonic touched the axis");
        const double p=power,l=std::log(x);
        const double factor=order==0?l:order==1?p*l+1:order==2?p*(p-1)*l+2*p-1:
            p*(p-1)*(p-2)*l+3*p*p-6*p+2;
        return std::pow(x,power-order)*factor;
    }
    template<class F>static void terms(int mode,F term)
    {
        const bool logarithm=mode>=4;const int degree=mode%4;
        term(1,0,degree,logarithm);
        if(degree>=2)term(degree==2?-.5:-1.5,2,degree-2,logarithm);
        if(logarithm && degree>=2)term(degree==2?.5:1.5,2,degree-2,false);
    }
    static Jet harmonic(int mode,double r,double z)
    {
        Jet result;
        const auto term=[&](double coefficient,int rp,int zp,bool logarithm=false) {
            const auto d=[&](int a,int b){return coefficient*derivative(rp,a,r,logarithm)*derivative(zp,b,z,false);};
            result.v+=d(0,0);result.r+=d(1,0);result.z+=d(0,1);
            result.rr+=d(2,0);result.rz+=d(1,1);result.zz+=d(0,2);
            result.rrr+=d(3,0);result.rrz+=d(2,1);result.rzz+=d(1,2);result.zzz+=d(0,3);
        };
        terms(mode,term);
        return result;
    }
    static long double meanPower(int power,long double midpoint,long double halfspan)
    {
        require(power>=0 && power<=3,"unsupported harmonic face monomial");
        if(power==0)return 1;
        if(power==1)return midpoint;
        if(power==2)return midpoint*midpoint+halfspan*halfspan/3;
        return midpoint*(midpoint*midpoint+halfspan*halfspan);
    }
    // Analytic logarithmic moments in a bounded log1p series. Subdivision
    // ensures |u/m|<=1/2; 64 terms leave an absolute series remainder below
    // 2e-20*m^p for p<=3. Rounding still applies; this is not an exact-arithmetic
    // certificate. No polynomial quadrature rule is applied to log(r).
    static long double meanLogPower(int power,long double midpoint,long double halfspan,int depth=0)
    {
        require(midpoint>std::abs(halfspan) && depth<=64,"unresolved positive-radius logarithmic face interval");
        const long double delta=std::abs(halfspan)/midpoint;
        if(delta>.5L)return .5L*(meanLogPower(power,midpoint-halfspan/2,halfspan/2,depth+1)
            +meanLogPower(power,midpoint+halfspan/2,halfspan/2,depth+1));
        long double value=std::log(midpoint)*meanPower(power,midpoint,halfspan),series=0;
        for(int k=0;k<=power;++k) {
            const int binomial=k==0 || k==power?1:power;
            for(int j=1;j<=64;++j)if((k+j)%2==0)
                series+=binomial*(j%2?1.L:-1.L)*std::pow(delta,k+j)/(j*(k+j+1.L));
        }
        return value+std::pow(midpoint,power)*series;
    }
    Row integratedBasis(const std::array<double,2>& midpoint,const std::array<long double,2>& span,
        double length,Phase phase,Vector2D direction) const
    {
        require((direction.x==0 && std::abs(direction.y)==1 && span[1]==0)
            || (direction.y==0 && std::abs(direction.x)==1 && span[0]==0),"harmonic flux needs an actual Cartesian face");
        require(length>0 && midpoint[0]>=0,"invalid harmonic physical face");
        Row result=Row::Zero();const double epsilon=phase==Phase::Negative?em_:ep_;
        const int first=phase==Phase::Negative?0:4,count=phase==Phase::Negative?4:8;
        for(int j=0;j<count;++j) {
            long double mean=0;
            terms(j,[&](double coefficient,int rp,int zp,bool logarithm) {
                if(direction.x)mean+=coefficient*midpoint[0]*direction.x*derivative(rp,1,midpoint[0],logarithm)
                    *meanPower(zp,midpoint[1],span[1]/2);
                else if(zp)mean+=coefficient*zp*std::pow(midpoint[1],zp-1)*direction.y
                    *(logarithm?meanLogPower(rp+1,midpoint[0],span[0]/2):meanPower(rp+1,midpoint[0],span[0]/2));
            });
            result[first+j]=static_cast<double>(-CartesianMeasure::twoPi*epsilon*h_*length*mean);
        }
        require(result.allFinite(),"axis transmission integrated flux overflow");return result;
    }
    Row pieceBasis(const SegmentPiece& piece,Vector2D direction) const
    {
        std::array<double,2> midpoint{};std::array<long double,2> span{};
        const std::array<double,2> a{{piece.begin.x/h_,(piece.begin.y-anchor_.y)/h_}},
            b{{piece.end.x/h_,(piece.end.y-anchor_.y)/h_}};
        const double fraction=piece.retainedRange?piece.fraction:1;
        require(fraction>0 && fraction<=1,"invalid retained harmonic face fraction");
        for(int axis=0;axis<2;++axis) {
            const long double whole=static_cast<long double>(b[axis])-a[axis];span[axis]=whole*fraction;
            midpoint[axis]=static_cast<double>(piece.retainedRange && piece.fromEnd?
                static_cast<long double>(b[axis])-whole*piece.offset-span[axis]/2:
                static_cast<long double>(a[axis])+whole*(piece.retainedRange?piece.offset:0)+span[axis]/2);
        }
        return integratedBasis(midpoint,span,std::hypot(static_cast<double>(span[0]),static_cast<double>(span[1])),piece.phase,direction);
    }
    Row correctionBasis(Vector2D owner,Phase ownerPhase,Vector2D neighbor,Phase neighborPhase,Vector2D normal,
        double length,double distance,double epsilon,const std::vector<SegmentPiece>& pieces,LinkGeometry geometry) const
    {
        Row integral=Row::Zero();long double measure=0;
        DielectricTraceFaceGeometry::visitValidatedFace(owner,ownerPhase,neighbor,neighborPhase,normal,length,distance,
            epsilon,pieces,h_,[&](Vector2D point){return std::array<double,2>{{point.x/h_,(point.y-anchor_.y)/h_}};},
            [&](const auto&,const auto&,const std::array<double,2>& midpoint,double partLength,Phase phase,const std::array<long double,2>& span) {
                integral+=integratedBasis(midpoint,span,partLength,phase,normal);
                measure+=CartesianMeasure::twoPi*h_*h_*partLength*midpoint[0];
            },"[DielectricAxisTransmission] ",geometry);
        const double conductance=static_cast<double>(epsilon*measure/distance);
        require(std::isfinite(conductance) && conductance>=0,"invalid harmonic two-point conductance");
        return integral-conductance*(potentialBasis(owner,ownerPhase)-potentialBasis(neighbor,neighborPhase));
    }
    std::array<Jet,12> jets(Vector2D point,Phase phase) const
    {
        require(std::isfinite(point.x) && std::isfinite(point.y) && point.x>=0,"invalid axis transmission query");
        DielectricCubicTrace::Detail::phase(phase);
        const double r=point.x/h_,z=(point.y-anchor_.y)/h_;
        std::array<Jet,12> result{};
        const int first=phase==Phase::Negative?0:4,count=phase==Phase::Negative?4:8;
        for(int j=0;j<count;++j)result[first+j]=harmonic(j,r,z);
        return result;
    }
    std::array<double,7> traces(const Jet& a,double epsilon) const
    {
        const auto n=normal_;const Vector2D t{-n.y,n.x};
        const auto H=[&](Vector2D u,Vector2D v){return u.x*v.x*a.rr+(u.x*v.y+u.y*v.x)*a.rz+u.y*v.y*a.zz;};
        const auto T=[&](Vector2D u,Vector2D v,Vector2D w){return u.x*v.x*w.x*a.rrr
            +(u.y*v.x*w.x+u.x*v.y*w.x+u.x*v.x*w.y)*a.rrz
            +(u.x*v.y*w.y+u.y*v.x*w.y+u.y*v.y*w.x)*a.rzz+u.y*v.y*w.y*a.zzz;};
        const double gn=n.x*a.r+n.y*a.z,gt=t.x*a.r+t.y*a.z;
        return {{a.v,gt,H(t,t)-k_*gn,T(t,t,t)-3*k_*H(n,t)-ks_*gn-k_*k_*gt,
            -epsilon*gn,-epsilon*(H(n,t)+k_*gt),
            -epsilon*(T(n,t,t)-k_*H(n,n)+2*k_*H(t,t)+ks_*gt-k_*k_*gn)}};
    }
    Row potentialBasis(Vector2D point,Phase phase) const
    {Row row;const auto a=jets(point,phase);for(int j=0;j<12;++j)row[j]=a[j].v;
        require(row.allFinite(),"axis transmission potential basis overflow");return row;}
    Row displacementBasis(Vector2D point,Phase phase,Vector2D direction) const
    {
        DielectricCubicTrace::Detail::unit(direction,"axis transmission displacement direction is not unit");
        Row row;const auto a=jets(point,phase);const double epsilon=phase==Phase::Negative?em_:ep_;
        for(int j=0;j<12;++j)row[j]=-epsilon*(direction.x*a[j].r+direction.y*a[j].z)/h_;
        require(row.allFinite(),"axis transmission displacement basis overflow");
        return row;
    }
    Eigen::Vector3d scaledJump(const Jump& q) const
    {
        const Eigen::Vector3d result(q.value*h_,q.firstDerivative*h_*h_,q.secondDerivative*h_*h_*h_);
        require(result.allFinite(),"nonfinite scaled transmission jump");return result;
    }
    Eigen::VectorXd coefficients(const std::vector<double>& values,const Jump& q) const
    {
        require(values.size()==sampleCount(),"axis transmission potential/sample count mismatch");
        Eigen::VectorXd differences(values.size());
        for(std::size_t i=0;i<values.size();++i) {
            require(std::isfinite(values[i]),"nonfinite axis transmission potential");differences[i]=values[i]-values[0];
        }
        const auto jump=scaledJump(q);
        Eigen::VectorXd result=homogeneous_*(differences-sampleParticular_*jump)+particular_*jump;
        require(result.allFinite(),"nonfinite axis transmission coefficients");return result;
    }
    std::vector<double> project(const Row& basis,bool potential) const
    {
        const Eigen::RowVectorXd raw=basis*homogeneous_;require(raw.allFinite(),"nonfinite axis transmission weights");
        std::vector<double> result(raw.data(),raw.data()+raw.size());long double sum=0;
        for(std::size_t i=1;i<result.size();++i)sum+=result[i];
        result[0]=static_cast<double>((potential?1.L:0.L)-sum);return result;
    }
public:
    Reconstruction(Vector2D point,Vector2D normal,double curvature,double curvatureDerivative,
        double epsilonNegative,double epsilonPositive,const std::vector<SamplePoint>& samples,double lengthScale)
        :anchor_(point),normal_(normal),h_(lengthScale),k_(curvature*lengthScale),
        ks_(curvatureDerivative*lengthScale*lengthScale),em_(epsilonNegative),ep_(epsilonPositive)
    {
        require(h_>0 && std::isfinite(h_) && point.x>0 && std::isfinite(point.x) && std::isfinite(point.y)
            && std::isfinite(k_) && std::isfinite(ks_) && em_>0 && ep_>0 && std::isfinite(em_) && std::isfinite(ep_),
            "invalid axis transmission geometry/material");
        DielectricCubicTrace::Detail::unit(normal,"axis transmission normal is not unit");
        require(samples.size()>=5,"axis transmission needs at least five actual samples");
        Eigen::MatrixXd constraint(7,12);Eigen::VectorXd scaling(7);
        for(Phase phase:{Phase::Negative,Phase::Positive}) {
            const auto a=jets(point,phase);const double sign=phase==Phase::Negative?-1:1;
            const int first=phase==Phase::Negative?0:4,last=phase==Phase::Negative?4:12;
            for(int j=first;j<last;++j) {
                const auto trace=traces(a[j],phase==Phase::Negative?em_:ep_);
                for(int row=0;row<7;++row)constraint(row,j)=sign*trace[row];
            }
        }
        for(int row=0;row<7;++row) {
            scaling[row]=1/constraint.row(row).norm();constraint.row(row)*=scaling[row];
        }
        require(constraint.allFinite(),"invalid axis transmission constraint scaling");
        const double tolerance=128*std::numeric_limits<double>::epsilon();
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(constraint,Eigen::ComputeFullU|Eigen::ComputeFullV);
        svd.setThreshold(tolerance);require(svd.rank()==7,"rank-deficient local transmission constraints");
        constraintCondition_=svd.singularValues()[0]/svd.singularValues()[6];
        require(std::isfinite(constraintCondition_) && constraintCondition_<=1e8,"ill-conditioned local transmission constraints");
        Eigen::MatrixXd rhs=Eigen::MatrixXd::Zero(7,3);
        for(int j=0;j<3;++j)rhs(j+4,j)=scaling[j+4];
        particular_=svd.solve(rhs);
        const Eigen::MatrixXd nullspace=svd.matrixV().rightCols(5);
        Eigen::MatrixXd sampling(samples.size(),12);
        for(std::size_t i=0;i<samples.size();++i) {
            for(std::size_t j=0;j<i;++j)require(std::hypot(samples[i].point.x-samples[j].point.x,
                samples[i].point.y-samples[j].point.y)>tolerance*h_,"duplicate/unresolved axis transmission samples");
            sampling.row(i)=potentialBasis(samples[i].point,samples[i].phase);
        }
        const Eigen::MatrixXd design=sampling*nullspace;
        require(design.allFinite(),"nonfinite axis transmission sampling design");
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(design);qr.setThreshold(tolerance);
        require(qr.rank()==5,"rank-deficient axis transmission sampling design");
        const Eigen::MatrixXd inverse=qr.solve(Eigen::MatrixXd::Identity(samples.size(),samples.size()));
        condition_=design.cwiseAbs().rowwise().sum().maxCoeff()*inverse.cwiseAbs().rowwise().sum().maxCoeff();
        reproduction_=(inverse*design-Eigen::MatrixXd::Identity(5,5)).cwiseAbs().maxCoeff();
        require(inverse.allFinite() && std::isfinite(condition_) && condition_>=1 && condition_<=1e8
            && reproduction_<=256*std::numeric_limits<double>::epsilon()*condition_,
            "axis transmission sampling condition/reproduction gate failed");
        const double constraintError=std::max((constraint*nullspace).cwiseAbs().maxCoeff(),
            (constraint*particular_-rhs).cwiseAbs().maxCoeff());
        require(constraintError<=256*std::numeric_limits<double>::epsilon()*constraintCondition_*std::max(1.,rhs.cwiseAbs().maxCoeff()),
            "axis transmission constraints lost reproduction");
        homogeneous_=nullspace*inverse;sampleParticular_=sampling*particular_;sampling_=std::move(sampling);
    }
    std::size_t sampleCount() const{return homogeneous_.cols();}
    double conditionEstimate() const{return condition_;}
    double constraintConditionEstimate() const{return constraintCondition_;}
    double normalizedReproductionError() const{return reproduction_;}
    double relativeResidual(const std::vector<double>& values,const Jump& q={}) const
    {
        const auto c=coefficients(values,q);Eigen::VectorXd actual(values.size());
        for(std::size_t i=0;i<values.size();++i)actual[i]=values[i]-values[0];
        const double residual=(sampling_*c-actual).norm(),scale=actual.norm();
        return scale>0?residual/scale:residual;
    }
    State evaluate(const std::vector<double>& values,const Jump& q={}) const
    {
        const auto c=coefficients(values,q);const auto a=jets(anchor_,Phase::Negative);
        std::array<double,7> trace{};
        for(int j=0;j<4;++j) {const auto t=traces(a[j],em_);for(int k=0;k<7;++k)trace[k]+=c[j]*t[k];}
        return {values[0]+trace[0],-trace[1]/h_,trace[4]/h_,-trace[2]/(h_*h_),trace[5]/(h_*h_),
            -trace[3]/(h_*h_*h_),trace[6]/(h_*h_*h_)};
    }
    double potential(Vector2D point,Phase phase,const std::vector<double>& values,const Jump& q={}) const
    {const auto c=coefficients(values,q);return values[0]+(potentialBasis(point,phase)*c)[0];}
    double displacementNormal(Vector2D point,Phase phase,Vector2D direction,const std::vector<double>& values,const Jump& q={}) const
    {return (displacementBasis(point,phase,direction)*coefficients(values,q))[0];}
    std::vector<double> potentialWeights(Vector2D point,Phase phase) const{return project(potentialBasis(point,phase),true);}
    std::vector<double> displacementNormalWeights(Vector2D point,Phase phase,Vector2D n) const{return project(displacementBasis(point,phase,n),false);}
    double particularPotential(Vector2D point,Phase phase,const Jump& q) const
    {return (potentialBasis(point,phase)*particular_*scaledJump(q))[0];}
    double particularDisplacementNormal(Vector2D point,Phase phase,Vector2D n,const Jump& q) const
    {return (displacementBasis(point,phase,n)*particular_*scaledJump(q))[0];}
    std::vector<double> integratedDisplacementWeights(const SegmentPiece& piece,Vector2D normal) const
    {return project(pieceBasis(piece,normal),false);}
    double particularIntegratedDisplacement(const SegmentPiece& piece,Vector2D normal,const Jump& q) const
    {return (pieceBasis(piece,normal)*particular_*scaledJump(q))[0];}
    std::vector<double> correctionWeights(Vector2D owner,Phase op,Vector2D neighbor,Phase np,Vector2D normal,
        double length,double distance,double epsilon,const std::vector<SegmentPiece>& pieces,
        LinkGeometry geometry=LinkGeometry::Orthogonal) const
    {return project(correctionBasis(owner,op,neighbor,np,normal,length,distance,epsilon,pieces,geometry),false);}
    DielectricCubicTrace::AffineCorrection affineCorrection(Vector2D owner,Phase op,Vector2D neighbor,Phase np,
        Vector2D normal,double length,double distance,double epsilon,const std::vector<SegmentPiece>& pieces,
        const Jump& q,LinkGeometry geometry=LinkGeometry::Orthogonal) const
    {
        const auto basis=correctionBasis(owner,op,neighbor,np,normal,length,distance,epsilon,pieces,geometry);
        DielectricCubicTrace::AffineCorrection result{project(basis,false),0};
        const auto jump=scaledJump(q);const Eigen::VectorXd sampled=sampleParticular_*jump;
        long double constant=(basis*particular_*jump)[0];
        for(std::size_t i=1;i<result.weights.size();++i)constant-=result.weights[i]*(static_cast<long double>(sampled[i])-sampled[0]);
        result.constant=static_cast<double>(constant);require(std::isfinite(result.constant),"harmonic affine flux overflow");return result;
    }
};
}
}
}
