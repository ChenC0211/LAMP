#pragma once
// Implementation-only volume-moment fit shared by geometry synchronization
// and height curvature. Physical columns/branch admissibility belong to each
// caller; this algebra never selects a surface or modifies the inventory.
#include "Core/Mesh/CartesianMeasure.h"
#include "Core/Geometry/PlanarCellCut.h"
#include <Eigen/QR>
#include <array>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>

namespace InterfaceHeightMoments
{
struct Column {long double tlo,thi,nlo,nhi,volume;};
// Connected normal support containing the anchor row. Truncate at an actual
// inactive cell; never sum through solid or replace it with a pure phase.
// Physical column bounds remain different across a stepped wall.
template<class Active> bool activeSpan(int centre,int lower,int upper,Active active,int& lo,int& hi)
{
    if (!active(centre)) return false;
    lo=hi=centre;
    while (lo>lower && active(lo-1)) --lo;
    while (hi<upper && active(hi+1)) ++hi;
    return true;
}
// Same absolute volume resolution as the PLIC inverse. This classifies
// geometric support only; every column sum still uses its full inventory.
inline long double volumeResolution(long double capacity)
{return Quantum::Geometry::PlanarCellCut::relativeVolumeResolution*capacity;}
inline int columnPhase(long double volume,long double capacity)
{
    const auto tolerance=volumeResolution(capacity);
    if (volume<=tolerance) return 1; // Geometrically gas.
    if (capacity-volume<=tolerance) return -1; // Geometrically liquid.
    return 0;
}
inline long double power(long double x,int degree)
{long double value=1;while (degree-->0) value*=x;return value;}
inline long double primitive(long double lo,long double hi,int degree)
{return (power(hi,degree+1)-power(lo,degree+1))/(degree+1);}
struct Jet {long double height=0,first=0,second=0;bool valid=false;};
struct Fit
{
    std::array<double,5> coefficients{};
    int count=0;bool even=false,squaredRadius=false,valid=false,implicitEven=false;
    long double origin=0,reference=0,spacing=0,scale=0;
    double rankRatio=0,scaledResidual=0;
    Jet evaluate(long double tangent) const
    {
        if (!valid) return {};
        const auto x=(tangent-origin)/spacing;long double f=0,df=0,ddf=0;
        if (implicitEven) {
            // y=a+b*x^2+c*y^2, with the graph root continuous at c=0.
            const long double a=coefficients[0],b=coefficients[1],c=coefficients[2];
            const auto q=a+b*x*x,D2=1-4*c*q;if (!(D2>0)) return {};
            const auto D=std::sqrt(D2);f=2*q/(1+D);df=2*b*x/D;
            ddf=2*b/D+2*c*df*df/D;
            f=reference+scale*f;df*=scale/spacing;ddf*=scale/(spacing*spacing);
            return {f,df,ddf,std::isfinite(f) && std::isfinite(df) && std::isfinite(ddf)};
        }
        for (int k=0;k<count;++k) {
            const int degree=even?2*k:k;const auto coefficient=coefficients[k];
            f+=coefficient*power(x,degree);
            if (degree>=1) df+=degree*coefficient*power(x,degree-1);
            if (degree>=2) ddf+=degree*(degree-1)*coefficient*power(x,degree-2);
        }
        f*=scale;df*=scale/spacing;ddf*=scale/(spacing*spacing);
        if (squaredRadius) {
            const auto q=reference*reference+f;if (!(q>0)) return {};
            f=std::sqrt(q);ddf=ddf/(2*f)-df*df/(4*f*f*f);df/=2*f;
        } else f+=reference;
        return {f,df,ddf,std::isfinite(f) && std::isfinite(df) && std::isfinite(ddf)};
    }
};
// Read-only axis curvature refinement. Three annular means determine a local
// even conic, including paraboloids and ellipsoids, without a sphere radius or
// centre input. Geometry synchronization retains its polynomial moment chart.
// An inadmissible nonlinear graph returns the existing valid polynomial fit.
inline Fit axisCurvatureFit(const std::vector<Column>& columns,const Fit& polynomial,int sign)
{
    if (!polynomial.valid || !polynomial.even || polynomial.squaredRadius || columns.size()!=3) return polynomial;
    Eigen::Vector3d target;std::array<double,3> lo{},hi{};
    // Use the original physical volume sums, not sampled graph heights.
    for (int i=0;i<3;++i) {
        lo[i]=std::pow(static_cast<double>(columns[i].tlo/polynomial.spacing),2);
        hi[i]=std::pow(static_cast<double>(columns[i].thi/polynomial.spacing),2);
        const auto& column=columns[i];
        const auto area=(CartesianMeasure::twoPi/2)*(column.thi-column.tlo)*(column.thi+column.tlo);
        target[i]=static_cast<double>(((sign>0?column.nlo+column.volume/area:column.nhi-column.volume/area)
            -polynomial.reference)/polynomial.scale);
    }
    // Parabolic start keeps the normal derivative positive everywhere; a
    // bounded Newton line search then matches the original annular means.
    Eigen::Vector3d coefficients(polynomial.coefficients[0],polynomial.coefficients[1],0);
    const auto moments=[&](const Eigen::Vector3d& c,Eigen::Vector3d& residual,Eigen::Matrix3d& jacobian) {
        for (int i=0;i<3;++i) {
            const double q0=c[0]+c[1]*lo[i],q1=c[0]+c[1]*hi[i];
            const double d0=1-4*c[2]*q0,d1=1-4*c[2]*q1;
            if (!(d0>0) || !(d1>0)) return false;
            const double D0=std::sqrt(d0),D1=std::sqrt(d1),sum=D0+D1;
            const double y0=2*q0/(1+D0),y1=2*q1/(1+D1);
            // Exact annular integral: r dr is proportional to d(x^2).
            // This divided-difference form also remains regular at c=0/b=0.
            residual[i]=(y0+y1)/2-c[2]*(y1-y0)*(y1-y0)/(3*sum)-target[i];
            jacobian(i,0)=2/sum;
            jacobian(i,1)=2/sum*(lo[i]+(hi[i]-lo[i])*(2*D0+D1)/(3*sum));
            jacobian(i,2)=2/sum*(y0*y0+y0*y1+y1*y1)/3;
        }
        return residual.allFinite() && jacobian.allFinite();
    };
    for (int iteration=0;iteration<20;++iteration) {
        Eigen::Vector3d residual;Eigen::Matrix3d jacobian;
        if (!moments(coefficients,residual,jacobian)) return polynomial;
        const double error=residual.cwiseAbs().maxCoeff();
        Eigen::ColPivHouseholderQR<Eigen::Matrix3d> qr(jacobian);qr.setThreshold(1e-10);
        if (qr.rank()!=3) return polynomial;
        const auto diagonal=qr.matrixR().diagonal().cwiseAbs().eval();
        const double ratio=diagonal.minCoeff()/diagonal.maxCoeff();
        if (ratio<1e-8) return polynomial;
        if (error<=128*std::numeric_limits<double>::epsilon()*(1+target.cwiseAbs().maxCoeff())) {
            Fit result=polynomial;result.implicitEven=true;result.rankRatio=ratio;result.scaledResidual=error;
            for(int k=0;k<3;++k) result.coefficients[k]=coefficients[k];
            return result;
        }
        const Eigen::Vector3d increment=qr.solve(-residual);if(!increment.allFinite()) return polynomial;
        bool accepted=false;
        for (double step=1;step>=1./1024;step*=.5) {
            const Eigen::Vector3d next=coefficients+step*increment;
            Eigen::Vector3d nextResidual;Eigen::Matrix3d nextJacobian;
            if (moments(next,nextResidual,nextJacobian) && nextResidual.cwiseAbs().maxCoeff()<error) {
                coefficients=next;accepted=true;break;
            }
        }
        if (!accepted) return polynomial;
    }
    return polynomial;
}
inline Fit fit(const std::vector<Column>& columns,long double origin,long double reference,
    double ht,double hn,int axis,int sign,bool rz,bool even,int coefficientCount=0,bool fixedOrigin=false)
{
    Fit result;const int rows=static_cast<int>(columns.size());
    result.count=coefficientCount?coefficientCount:rows;result.even=even;
    result.squaredRadius=rz && axis==0;result.origin=origin;result.reference=reference;result.spacing=ht;
    result.scale=result.squaredRadius?hn*std::max(static_cast<double>(reference),hn):hn;
    const bool supported=even?(rows==3 && result.count==3):
        (((rows==3 || rows==4) && result.count==rows) || ((rows==5 || rows==7) && result.count==5));
    if (!supported || (fixedOrigin && even) || !(ht>0) || !(hn>0) || !(result.scale>0)) return result;
    // A declared wall endpoint supplies H(origin)=reference. Eliminate the
    // constant coefficient instead of adding a large penalty or ghost volume.
    const int first=fixedOrigin?1:0,unknowns=result.count-first;
    Eigen::MatrixXd matrix(rows,unknowns);Eigen::VectorXd rhs(rows);
    for (int i=0;i<rows;++i) {
        const auto& c=columns[i];long double value=0;
        if (result.squaredRadius) {
            const auto q=c.volume/((CartesianMeasure::twoPi/2)*(c.thi-c.tlo));
            value=sign>0?(c.nlo-reference)*(c.nlo+reference)+q:(c.nhi-reference)*(c.nhi+reference)-q;
        } else {
            const long double width=rz?(CartesianMeasure::twoPi/2)*(c.thi-c.tlo)*(c.thi+c.tlo):c.thi-c.tlo;
            value=sign>0?(c.nlo-reference)+c.volume/width:(c.nhi-reference)-c.volume/width;
        }
        rhs[i]=static_cast<double>(value/result.scale);
        const auto xlo=(c.tlo-origin)/ht,xhi=(c.thi-origin)/ht;
        const bool weighted=rz && axis==1;
        const auto denominator=weighted?origin*(xhi-xlo)+ht*primitive(xlo,xhi,1):xhi-xlo;
        for (int k=first;k<result.count;++k) {
            const int degree=even?2*k:k;
            const auto moment=weighted?origin*primitive(xlo,xhi,degree)+ht*primitive(xlo,xhi,degree+1):primitive(xlo,xhi,degree);
            matrix(i,k-first)=static_cast<double>(moment/denominator);
        }
    }
    if (!matrix.allFinite() || !rhs.allFinite()) return result;
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(matrix);qr.setThreshold(1e-10);
    if (qr.rank()!=unknowns) return result;
    const auto diagonal=qr.matrixR().topLeftCorner(unknowns,unknowns).diagonal().cwiseAbs().eval();
    result.rankRatio=diagonal.minCoeff()/diagonal.maxCoeff();
    const Eigen::VectorXd coefficients=qr.solve(rhs);
    if (result.rankRatio<1e-8 || !coefficients.allFinite()) return result;
    const auto residual=(matrix*coefficients-rhs).eval();
    result.scaledResidual=residual.cwiseAbs().maxCoeff();
    // Exact column interpolation remains the synchronization contract. An
    // explicitly requested overdetermined fit retains its residual for the
    // caller; algebraic fit quality alone does not prove a restoring symbol.
    if (rows==unknowns && residual.norm()>1e-10*(1+rhs.norm())) return result;
    for (int k=first;k<result.count;++k) result.coefficients[k]=coefficients[k-first];
    result.valid=true;return result;
}
}
