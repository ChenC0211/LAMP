#pragma once

#include "Modules/Electricity/Detail/DielectricCubicTrace.h"
#include <Eigen/SVD>

namespace Quantum::Electricity::Detail
{
namespace
{
namespace DielectricCapTransmission
{
using Jump=DielectricCubicTrace::Jump;
using State=DielectricCubicTrace::State;
inline void require(bool value,const char* message)
{if(!value)throw std::runtime_error(std::string("[DielectricCapTransmission] ")+message);}

// Algebraic transmission operator for one or two axisymmetric caps. It does
// NOT discover geometry or authorize a subgrid model. Its caller must certify
// each expansion centre lies strictly inside its labelled liquid, every
// sample/query label, and all gas integration segments exclude those centres.
// No physical/global entry currently enables this operator.
//
// One regular degree-three harmonic potential per liquid. The COMMON gas
// potential has four regular harmonics and four decaying solid harmonics per
// liquid centre. Seven local curved transmission-jet conditions per interface
// leave 4+N homogeneous modes. Constraints hold at the supplied anchors, not
// pointwise along arbitrary cap shapes. Independent charge jets are evaluated
// anew; the retained geometry operator never owns an evolving charge callback.
struct Anchor {
    Vector2D point,normal;
    double curvature=0,curvatureDerivative=0,centerZ=0;
};
struct Sample {Vector2D point;int liquid=-1;}; // -1 is the common gas.
struct Piece {DielectricCubicTrace::SegmentPiece geometry;int liquid=-1;};

class Reconstruction
{
    using Row=Eigen::RowVectorXd;
    // Taylor coefficients, i.e. derivatives divided by i!j!, through order 3.
    struct Jet {
        std::array<double,16> a{};
        explicit Jet(double value=0){a[0]=value;}
        static Jet variable(double value,int direction){Jet q(value);q.a[direction?1:4]=1;return q;}
        friend Jet operator+(Jet u,const Jet& v){for(int i=0;i<16;++i)u.a[i]+=v.a[i];return u;}
        friend Jet operator*(const Jet& u,const Jet& v){
            Jet q;
            for(int i=0;i<=3;++i)for(int j=0;j<=3-i;++j)
                for(int k=0;k<=i;++k)for(int l=0;l<=j;++l)q.a[4*i+j]+=u.a[4*k+l]*v.a[4*(i-k)+j-l];
            return q;
        }
        friend Jet operator*(double s,Jet q){for(auto& v:q.a)v*=s;return q;}
        double d(int i,int j) const {constexpr int factorial[4]={1,1,2,6};return a[4*i+j]*factorial[i]*factorial[j];}
        static Jet power(const Jet& q,double exponent){
            require(q.a[0]>0 && std::isfinite(q.a[0]),"gas harmonic touched an interior expansion centre");
            Jet u=(1/q.a[0])*q;u.a[0]=0;
            return std::pow(q.a[0],exponent)*(Jet(1)+exponent*u
                +(.5*exponent*(exponent-1))*(u*u)
                +(exponent*(exponent-1)*(exponent-2)/6)*(u*u*u));
        }
    };
    std::vector<Anchor> anchors_;
    double h_,em_,ep_,origin_,condition_=0,constraintCondition_=0,reproduction_=0;
    int columns_=0;
    Eigen::MatrixXd homogeneous_,particular_,sampleParticular_,sampling_;

    static std::array<Jet,4> regular(const Jet& r,const Jet& z)
    {return {Jet(1),z,z*z+(-.5)*(r*r),z*z*z+(-1.5)*(z*r*r)};}
    std::vector<Jet> basis(Vector2D point,int liquid) const {
        require(std::isfinite(point.x)&&std::isfinite(point.y)&&point.x>=0
            &&liquid>=-1&&liquid<static_cast<int>(anchors_.size()),"invalid cap sample or liquid label");
        std::vector<Jet> result(columns_);
        const auto r=Jet::variable(point.x/h_,0);
        if(liquid>=0){
            const auto p=regular(r,Jet::variable((point.y-anchors_[liquid].centerZ)/h_,1));
            for(int j=0;j<4;++j)result[4*liquid+j]=p[j];
        }else{
            const int start=4*static_cast<int>(anchors_.size());
            const auto p=regular(r,Jet::variable((point.y-origin_)/h_,1));
            for(int j=0;j<4;++j)result[start+j]=p[j];
            for(std::size_t k=0;k<anchors_.size();++k){
                const auto z=Jet::variable((point.y-anchors_[k].centerZ)/h_,1),square=r*r+z*z;
                const auto p=regular(r,z);
                for(int j=0;j<4;++j)result[start+4+4*k+j]=p[j]*Jet::power(square,-j-.5);
            }
        }
        for(const auto& q:result)for(double value:q.a)require(std::isfinite(value),"nonfinite cap harmonic jet");
        return result;
    }
    std::array<double,7> trace(const Jet& q,const Anchor& a,double epsilon) const {
        const auto n=a.normal;const Vector2D t{-n.y,n.x};const double k=a.curvature*h_,ks=a.curvatureDerivative*h_*h_;
        const auto H=[&](Vector2D u,Vector2D v){return u.x*v.x*q.d(2,0)+(u.x*v.y+u.y*v.x)*q.d(1,1)+u.y*v.y*q.d(0,2);};
        const auto T=[&](Vector2D u,Vector2D v,Vector2D w){return u.x*v.x*w.x*q.d(3,0)
            +(u.y*v.x*w.x+u.x*v.y*w.x+u.x*v.x*w.y)*q.d(2,1)
            +(u.x*v.y*w.y+u.y*v.x*w.y+u.y*v.y*w.x)*q.d(1,2)+u.y*v.y*w.y*q.d(0,3);};
        const double gn=n.x*q.d(1,0)+n.y*q.d(0,1),gt=t.x*q.d(1,0)+t.y*q.d(0,1);
        return {q.a[0],gt,H(t,t)-k*gn,T(t,t,t)-3*k*H(n,t)-ks*gn-k*k*gt,
            -epsilon*gn,-epsilon*(H(n,t)+k*gt),
            -epsilon*(T(n,t,t)-k*H(n,n)+2*k*H(t,t)+ks*gt-k*k*gn)};
    }
    Row potentialBasis(Vector2D point,int liquid) const {
        const auto q=basis(point,liquid);Row row(columns_);for(int i=0;i<columns_;++i)row[i]=q[i].a[0];return row;
    }
    Row displacementBasis(Vector2D point,int liquid,Vector2D normal) const {
        DielectricCubicTrace::Detail::unit(normal,"cap flux direction is not unit");
        const auto q=basis(point,liquid);Row row(columns_);const double epsilon=liquid<0?ep_:em_;
        for(int i=0;i<columns_;++i)row[i]=-epsilon*(normal.x*q[i].d(1,0)+normal.y*q[i].d(0,1))/h_;
        return row;
    }
    Eigen::VectorXd scaledJumps(const std::vector<Jump>& jumps) const {
        require(jumps.size()==anchors_.size(),"cap charge/anchor count mismatch");Eigen::VectorXd q(3*jumps.size());
        for(std::size_t i=0;i<jumps.size();++i){q[3*i]=jumps[i].value*h_;q[3*i+1]=jumps[i].firstDerivative*h_*h_;
            q[3*i+2]=jumps[i].secondDerivative*h_*h_*h_;}
        require(q.allFinite(),"nonfinite cap charge jet");return q;
    }
    Eigen::VectorXd coefficients(const std::vector<double>& values,const std::vector<Jump>& jumps) const {
        require(values.size()==sampleCount(),"cap potential/sample count mismatch");Eigen::VectorXd v(values.size());
        for(std::size_t i=0;i<values.size();++i){require(std::isfinite(values[i]),"nonfinite cap potential");v[i]=values[i]-values[0];}
        const auto q=scaledJumps(jumps);Eigen::VectorXd c=homogeneous_*(v-sampleParticular_*q)+particular_*q;
        require(c.allFinite(),"nonfinite cap coefficients");return c;
    }
    std::vector<double> weights(const Row& row,bool potential) const {
        const Row raw=row*homogeneous_;require(raw.allFinite(),"nonfinite cap functional");
        std::vector<double> w(raw.data(),raw.data()+raw.size());long double sum=0;
        for(std::size_t i=1;i<w.size();++i)sum+=w[i];w[0]=static_cast<double>((potential?1.L:0.L)-sum);return w;
    }
    // Stream function satisfies Psi_r=r*V_z, Psi_z=-r*V_r. Differences
    // integrate an entire Cartesian physical face, including the 2*pi*r
    // measure. Midpoint and span remain separate for sub-ULP retained pieces.
    static long double streamDifference(int mode,bool decaying,long double r,long double z,
        long double dr,long double dz) {
        if(!decaying){
            if(mode==0)return 0;
            if(mode==1)return r*dr;
            if(mode==2)return dz*r*r+2*z*r*dr;
            return 3*z*dz*r*r+dr*r*(3*z*z-1.5L*r*r-.375L*dr*dr);
        }
        const long double ra=r-dr/2,rb=r+dr/2,za=z-dz/2,zb=z+dz/2;
        const long double rhoA=std::hypot(ra,za),rhoB=std::hypot(rb,zb);
        require(rhoA>0&&rhoB>0 && !(r==0&&dr==0&&std::abs(z)<=std::abs(dz)/2),
            "gas face intersects a liquid expansion centre");
        if(dr==0&&r==0)return 0;
        if(mode==0&&dr==0){
            // Rationalization avoids cancellation of z/rho on a thin axial
            // face far from the source. Across z=0 the direct difference is
            // well conditioned; no small width is recovered by subtraction.
            if((za>0&&zb>0)||(za<0&&zb<0))
                return r*r*dz*(2*z)/(rhoA*rhoB*(zb*rhoA+za*rhoB));
            return zb/rhoB-za/rhoA;
        }
        const int power=2*mode+1;
        const long double deltaRho=(2*r*dr+2*z*dz)/(rhoA+rhoB);
        long double dividedPower=0;
        for(int j=0;j<power;++j)dividedPower+=std::pow(rhoB,power-1-j)*std::pow(rhoA,j);
        const long double pa=std::pow(rhoA,power),pb=std::pow(rhoB,power);
        const long double inverseMean=.5L*(1/pa+1/pb),inverseDifference=-deltaRho*dividedPower/(pa*pb);
        long double mean=0,difference=0;
        if(mode==0){mean=z;difference=dz;}
        else if(mode==1){mean=-(r*r+dr*dr/4);difference=-2*r*dr;}
        else if(mode==2){mean=-1.5L*z*(r*r+dr*dr/4);difference=-1.5L*(dz*r*r+2*z*r*dr);}
        else {
            mean=.5L*(r*r*r*r+1.5L*r*r*dr*dr+dr*dr*dr*dr/16)
                -2*(z*z+dz*dz/4)*(r*r+dr*dr/4);
            difference=dr*r*(2*r*r+.5L*dr*dr-4*z*z)-4*z*dz*r*r;
        }
        return difference*inverseMean+mean*inverseDifference;
    }
    Row integratedBasis(const std::array<double,2>& midpoint,const std::array<long double,2>& span,
        int liquid,Vector2D normal) const {
        require(liquid>=-1&&liquid<static_cast<int>(anchors_.size()),"invalid cap face liquid label");
        require((normal.x==0&&std::abs(normal.y)==1&&span[1]==0)
            ||(normal.y==0&&std::abs(normal.x)==1&&span[0]==0),"cap flux requires an actual Cartesian face");
        const long double r=midpoint[0],dr=std::abs(span[0]),dz=std::abs(span[1]);
        require(std::isfinite(r)&&r>=0&&(dr>0||dz>0)&&r>=dr/2,"invalid cap face span/radius");
        Row row=Row::Zero(columns_);const double epsilon=liquid<0?ep_:em_;
        const long double factor=CartesianMeasure::twoPi*epsilon*h_*(normal.x-normal.y);
        const auto fill=[&](int start,bool decaying,double center){
            const long double z=static_cast<long double>(midpoint[1])-(static_cast<long double>(center)-origin_)/h_;
            for(int j=0;j<4;++j)row[start+j]=static_cast<double>(factor*streamDifference(j,decaying,r,z,dr,dz));
        };
        if(liquid>=0)fill(4*liquid,false,anchors_[liquid].centerZ);
        else {
            const int start=4*static_cast<int>(anchors_.size());fill(start,false,origin_);
            for(std::size_t i=0;i<anchors_.size();++i)fill(start+4+4*static_cast<int>(i),true,anchors_[i].centerZ);
        }
        require(row.allFinite(),"nonfinite integrated cap flux");return row;
    }
    Row pieceBasis(const Piece& piece,Vector2D normal) const {
        const auto& p=piece.geometry;
        require((piece.liquid<0)==(p.phase==DielectricCubicTrace::Phase::Positive),"cap face phase/region mismatch");
        const std::array<double,2> a{{p.begin.x/h_,(p.begin.y-origin_)/h_}},b{{p.end.x/h_,(p.end.y-origin_)/h_}};
        const long double fraction=p.retainedRange?p.fraction:1;
        require(std::isfinite(fraction)&&fraction>0&&fraction<=1
            &&(!p.retainedRange||(std::isfinite(p.offset)&&p.offset>=0
                &&p.offset+fraction<=1+256*std::numeric_limits<double>::epsilon())),"invalid retained cap face piece");
        std::array<double,2> midpoint{};std::array<long double,2> span{};
        for(int i=0;i<2;++i){const long double whole=static_cast<long double>(b[i])-a[i];span[i]=whole*fraction;
            midpoint[i]=static_cast<double>(p.retainedRange&&p.fromEnd?
                static_cast<long double>(b[i])-whole*p.offset-span[i]/2:
                static_cast<long double>(a[i])+whole*(p.retainedRange?p.offset:0)+span[i]/2);}
        return integratedBasis(midpoint,span,piece.liquid,normal);
    }
    Row correctionBasis(Vector2D owner,int ownerLiquid,Vector2D neighbor,int neighborLiquid,Vector2D normal,
        double length,double distance,double epsilon,const std::vector<Piece>& pieces,DielectricCubicTrace::LinkGeometry geometry) const {
        using Phase=DielectricCubicTrace::Phase;std::vector<DielectricCubicTrace::SegmentPiece> geometryPieces;
        for(const auto& piece:pieces){require((piece.liquid<0)==(piece.geometry.phase==Phase::Positive),"cap correction phase/region mismatch");
            geometryPieces.push_back(piece.geometry);}
        Row row=Row::Zero(columns_);long double area=0;std::size_t index=0;
        DielectricTraceFaceGeometry::visitValidatedFace(owner,ownerLiquid<0?Phase::Positive:Phase::Negative,
            neighbor,neighborLiquid<0?Phase::Positive:Phase::Negative,normal,length,distance,epsilon,geometryPieces,h_,
            [&](Vector2D p){return std::array<double,2>{{p.x/h_,(p.y-origin_)/h_}};},
            [&](const auto&,const auto&,const std::array<double,2>& midpoint,double partLength,Phase,
                const std::array<long double,2>& span){
                row+=integratedBasis(midpoint,span,pieces.at(index++).liquid,normal);
                area+=CartesianMeasure::twoPi*h_*h_*partLength*midpoint[0];
            },"[DielectricCapTransmission] ",geometry);
        const double conductance=static_cast<double>(epsilon*area/distance);
        require(std::isfinite(conductance)&&conductance>=0,"invalid cap two-point conductance");
        return row-conductance*(potentialBasis(owner,ownerLiquid)-potentialBasis(neighbor,neighborLiquid));
    }
public:
    Reconstruction(std::vector<Anchor> anchors,const std::vector<Sample>& samples,double h,double em,double ep)
        :anchors_(std::move(anchors)),h_(h),em_(em),ep_(ep),origin_(anchors_.empty()?0:anchors_[0].point.y)
    {
        const int count=static_cast<int>(anchors_.size()),constraints=7*count,modes=4+count;columns_=8*count+4;
        require(count>=1&&count<=2&&std::isfinite(h_)&&h_>0&&std::isfinite(em_)&&em_>0
            &&std::isfinite(ep_)&&ep_>0&&samples.size()>=static_cast<std::size_t>(modes),"invalid cap transmission support/material");
        Eigen::MatrixXd c=Eigen::MatrixXd::Zero(constraints,columns_),rhs=Eigen::MatrixXd::Zero(constraints,3*count);
        for(int i=0;i<count;++i){
            const auto& a=anchors_[i];
            require(a.point.x>0&&std::isfinite(a.point.x)&&std::isfinite(a.point.y)&&std::isfinite(a.centerZ)
                &&std::isfinite(a.curvature)&&std::isfinite(a.curvatureDerivative),"invalid nonpolar cap anchor");
            DielectricCubicTrace::Detail::unit(a.normal,"cap normal is not unit");
            for(int j=0;j<i;++j)require(std::abs(a.centerZ-anchors_[j].centerZ)>128*std::numeric_limits<double>::epsilon()*h_,
                "duplicate cap expansion centres");
            for(int side:{i,-1}){
                const auto q=basis(a.point,side);
                for(int j=0;j<columns_;++j){const auto t=trace(q[j],a,side<0?ep_:em_);
                    for(int k=0;k<7;++k)c(7*i+k,j)+=(side<0?1:-1)*t[k];}
            }
            for(int k=0;k<7;++k){const double scale=c.row(7*i+k).norm();
                require(std::isfinite(scale)&&scale>0,"degenerate cap transmission constraint");c.row(7*i+k)/=scale;
                if(k>=4)rhs(7*i+k,3*i+k-4)=1/scale;}
        }
        const double tolerance=128*std::numeric_limits<double>::epsilon();
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(c,Eigen::ComputeFullU|Eigen::ComputeFullV);svd.setThreshold(tolerance);
        require(svd.rank()==constraints,"rank-deficient coupled cap constraints");
        constraintCondition_=svd.singularValues()[0]/svd.singularValues()[constraints-1];
        require(std::isfinite(constraintCondition_)&&constraintCondition_<=1e8,"ill-conditioned coupled cap constraints");
        const Eigen::MatrixXd nullspace=svd.matrixV().rightCols(modes);particular_=svd.solve(rhs);
        sampling_.resize(samples.size(),columns_);
        for(std::size_t i=0;i<samples.size();++i){
            for(std::size_t j=0;j<i;++j)require(std::hypot(samples[i].point.x-samples[j].point.x,
                samples[i].point.y-samples[j].point.y)>tolerance*h_,"duplicate cap samples");
            sampling_.row(i)=potentialBasis(samples[i].point,samples[i].liquid);
        }
        const Eigen::MatrixXd design=sampling_*nullspace;
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(design);qr.setThreshold(tolerance);
        require(qr.rank()==modes,"rank-deficient coupled cap sampling");
        const Eigen::MatrixXd inverse=qr.solve(Eigen::MatrixXd::Identity(samples.size(),samples.size()));
        condition_=design.cwiseAbs().rowwise().sum().maxCoeff()*inverse.cwiseAbs().rowwise().sum().maxCoeff();
        reproduction_=(inverse*design-Eigen::MatrixXd::Identity(modes,modes)).cwiseAbs().maxCoeff();
        require(inverse.allFinite()&&std::isfinite(condition_)&&condition_<=1e8
            &&reproduction_<=256*std::numeric_limits<double>::epsilon()*std::max(1.,condition_),"cap fit condition/reproduction failed");
        require(std::max((c*nullspace).cwiseAbs().maxCoeff(),(c*particular_-rhs).cwiseAbs().maxCoeff())
            <=256*std::numeric_limits<double>::epsilon()*constraintCondition_*std::max(1.,rhs.cwiseAbs().maxCoeff()),
            "cap constraints lost reproduction");
        homogeneous_=nullspace*inverse;sampleParticular_=sampling_*particular_;
    }
    std::size_t sampleCount() const{return homogeneous_.cols();}
    double conditionEstimate() const{return condition_;}
    double constraintConditionEstimate() const{return constraintCondition_;}
    double normalizedReproductionError() const{return reproduction_;}
    double potential(Vector2D point,int liquid,const std::vector<double>& values,const std::vector<Jump>& jumps) const
    {return values.at(0)+(potentialBasis(point,liquid)*coefficients(values,jumps))[0];}
    double displacementNormal(Vector2D point,int liquid,Vector2D normal,const std::vector<double>& values,const std::vector<Jump>& jumps) const
    {return (displacementBasis(point,liquid,normal)*coefficients(values,jumps))[0];}
    State evaluate(std::size_t anchor,const std::vector<double>& values,const std::vector<Jump>& jumps) const {
        const auto& a=anchors_.at(anchor);const auto q=basis(a.point,static_cast<int>(anchor));const auto c=coefficients(values,jumps);
        std::array<double,7> t{};for(int j=0;j<columns_;++j){const auto v=trace(q[j],a,em_);for(int k=0;k<7;++k)t[k]+=c[j]*v[k];}
        return {values[0]+t[0],-t[1]/h_,t[4]/h_,-t[2]/(h_*h_),t[5]/(h_*h_),-t[3]/(h_*h_*h_),t[6]/(h_*h_*h_)};
    }
    std::vector<double> potentialWeights(Vector2D point,int liquid) const{return weights(potentialBasis(point,liquid),true);}
    std::vector<double> displacementNormalWeights(Vector2D point,int liquid,Vector2D normal) const{return weights(displacementBasis(point,liquid,normal),false);}
    double particularPotential(Vector2D point,int liquid,const std::vector<Jump>& jumps) const
    {return (potentialBasis(point,liquid)*particular_*scaledJumps(jumps))[0];}
    double particularDisplacementNormal(Vector2D point,int liquid,Vector2D normal,const std::vector<Jump>& jumps) const
    {return (displacementBasis(point,liquid,normal)*particular_*scaledJumps(jumps))[0];}
    double integratedDisplacement(const Piece& piece,Vector2D normal,const std::vector<double>& values,const std::vector<Jump>& jumps) const
    {return (pieceBasis(piece,normal)*coefficients(values,jumps))[0];}
    std::vector<double> integratedDisplacementWeights(const Piece& piece,Vector2D normal) const
    {return weights(pieceBasis(piece,normal),false);}
    double particularIntegratedDisplacement(const Piece& piece,Vector2D normal,const std::vector<Jump>& jumps) const
    {return (pieceBasis(piece,normal)*particular_*scaledJumps(jumps))[0];}
    DielectricCubicTrace::AffineCorrection affineCorrection(Vector2D owner,int ownerLiquid,Vector2D neighbor,int neighborLiquid,
        Vector2D normal,double length,double distance,double epsilon,const std::vector<Piece>& pieces,
        const std::vector<Jump>& jumps,DielectricCubicTrace::LinkGeometry geometry=DielectricCubicTrace::LinkGeometry::Orthogonal) const {
        const auto row=correctionBasis(owner,ownerLiquid,neighbor,neighborLiquid,normal,length,distance,epsilon,pieces,geometry);
        DielectricCubicTrace::AffineCorrection result{weights(row,false),0};const auto q=scaledJumps(jumps);
        const Eigen::VectorXd sampled=sampleParticular_*q;long double constant=(row*particular_*q)[0];
        for(std::size_t i=1;i<result.weights.size();++i)constant-=result.weights[i]*(static_cast<long double>(sampled[i])-sampled[0]);
        result.constant=static_cast<double>(constant);require(std::isfinite(result.constant),"nonfinite affine cap flux");return result;
    }
    double relativeResidual(const std::vector<double>& values,const std::vector<Jump>& jumps) const {
        const auto c=coefficients(values,jumps);Eigen::VectorXd v(values.size());for(std::size_t i=0;i<values.size();++i)v[i]=values[i]-values[0];
        const double error=(sampling_*c-v).norm(),scale=v.norm();return scale>0?error/scale:error;
    }
};
}
}
}
