#include "Core/Geometry/ArcScalarRegularization.h"
#include <Eigen/SparseQR>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace Quantum::Geometry::ArcScalarRegularization
{
namespace {
Result solve(const std::vector<double>& lengths,const std::vector<double>& values,double h,
    const std::vector<int>& columns,int unknowns,bool periodic=true)
{
    const auto require=[](bool ok,const char* s){if(!ok)throw std::invalid_argument(s);};
    require(lengths.size()+(periodic?0:1)==values.size() && values.size()>=4 && std::isfinite(h) && h>0,
        "[ArcScalarRegularization] invalid support or scale");
    const int n=static_cast<int>(values.size());std::vector<double> widths(n),mass(n);
    for(double value:values)require(std::isfinite(value),"[ArcScalarRegularization] nonfinite scalar");
    for(int i=0;i<int(lengths.size());++i){require(std::isfinite(lengths[i]) && lengths[i]>0,
        "[ArcScalarRegularization] nonpositive edge");widths[i]=lengths[i]/h;
        require(std::isfinite(widths[i]) && widths[i]>0,"[ArcScalarRegularization] unrepresentable normalized edge");}
    long double volume=0,sum=0;for(int i=0;i<n;++i){mass[i]=(widths[(i+n-1)%n]+widths[i])/2;volume+=mass[i];sum+=mass[i]*static_cast<long double>(values[i]);}
    require(std::isfinite(volume) && volume>0 && std::isfinite(sum),"[ArcScalarRegularization] unrepresentable mass");
    const double mean=static_cast<double>(sum/volume);
    const int derivativeRows=periodic?n:n-3;
    Eigen::VectorXd rhs=Eigen::VectorXd::Zero(n+derivativeRows);std::vector<Eigen::Triplet<double>> entries;entries.reserve(5*n);
    const auto entry=[&](int row,int knot,double value){
        if(columns[knot]>=0)entries.emplace_back(row,columns[knot],value);
        else rhs[row]-=value*(values[knot]-mean);
    };
    for(int i=0;i<n;++i){rhs[i]+=std::sqrt(mass[i])*(values[i]-mean);entry(i,i,std::sqrt(mass[i]));}
    for(int row=0;row<derivativeRows;++row){const int i=periodic?row:row+1;
        // The third divided difference annihilates every quadratic on an
        // arbitrary positive knot sequence. A gradient of a lumped FV
        // Laplacian does not: rapidly varying cut-cell arc lengths turn its
        // second-derivative truncation error into spurious roughness.
        // Integrate the cubic's constant third derivative over the central
        // edge. Reflection and orientation reversal retain the same rows.
        const int knots[4]={(i+n-1)%n,i,(i+1)%n,(i+2)%n};
        const double x[4]={-widths[knots[0]],0,widths[i],widths[i]+widths[knots[2]]};
        double measure=widths[i];
        if(!periodic){if(i==1)measure+=widths[0];if(i==n-3)measure+=widths[n-2];}
        for(int k=0;k<4;++k){double coefficient=6*std::sqrt(measure);
            for(int j=0;j<4;++j)if(j!=k)coefficient/=x[k]-x[j];
            require(std::isfinite(coefficient),"[ArcScalarRegularization] unrepresentable derivative row");
            entry(n+row,knots[k],coefficient);}}
    Eigen::SparseMatrix<double> matrix(n+derivativeRows,unknowns);matrix.setFromTriplets(entries.begin(),entries.end());
    Eigen::SparseQR<Eigen::SparseMatrix<double>,Eigen::COLAMDOrdering<int>> qr;
    qr.compute(matrix);require(qr.info()==Eigen::Success && qr.rank()==unknowns,"[ArcScalarRegularization] unresolved numerical rank");
    const Eigen::VectorXd solved=qr.solve(rhs);require(qr.info()==Eigen::Success && solved.allFinite(),"[ArcScalarRegularization] invalid solve");
    const Eigen::VectorXd residual=matrix*solved-rhs,normal=matrix.transpose()*residual,normalScale=matrix.transpose()*rhs;
    Result result;result.relativeNormalResidual=normal.lpNorm<Eigen::Infinity>()/std::max(1.,normalScale.lpNorm<Eigen::Infinity>());
    result.dataEnergy=residual.head(n).squaredNorm();result.derivativeEnergy=residual.tail(derivativeRows).squaredNorm();result.values.resize(n);
    long double change=0;for(int i=0;i<n;++i){result.values[i]=columns[i]<0?values[i]:mean+solved[columns[i]];
        change+=mass[i]*(static_cast<long double>(result.values[i])-values[i]);}
    result.weightedMeanChange=static_cast<double>(change/volume);
    require(std::isfinite(result.relativeNormalResidual) && std::isfinite(result.weightedMeanChange),
        "[ArcScalarRegularization] nonfinite solve diagnostics");
    return result;
}
}
Result closedThirdDerivative(const std::vector<double>& lengths,const std::vector<double>& values,double h)
{
    std::vector<int> columns(values.size());for(int i=0;i<int(values.size());++i)columns[i]=i;
    return solve(lengths,values,h,columns,static_cast<int>(values.size()));
}
Result evenOpenThirdDerivative(const std::vector<double>& lengths,const std::vector<double>& values,double h)
{
    if(values.size()<3 || lengths.size()+1!=values.size())
        throw std::invalid_argument("[ArcScalarRegularization] invalid even interval");
    const int count=static_cast<int>(values.size()),n=2*(count-1);
    std::vector<double> reflected(n),edges(n);std::vector<int> columns(n);
    for(int i=0;i<n;++i){const int index=i<count?i:n-i;reflected[i]=values[index];
        columns[i]=(index==0 || index==count-1)?-1:index-1;edges[i]=lengths[i<count-1?i:n-1-i];}
    auto result=solve(edges,reflected,h,columns,count-2);result.values.resize(count);
    result.dataEnergy/=2;result.derivativeEnergy/=2;return result;
}
Result openThirdDerivative(const std::vector<double>& lengths,const std::vector<double>& values,double h,
    Endpoint begin,Endpoint end)
{
    if(values.size()<4 || lengths.size()+1!=values.size()
        || (begin!=Endpoint::FixedValue && begin!=Endpoint::FixedEven)
        || (end!=Endpoint::FixedValue && end!=Endpoint::FixedEven))
        throw std::invalid_argument("[ArcScalarRegularization] invalid open interval or endpoint");
    if(begin==Endpoint::FixedEven && end==Endpoint::FixedEven)return evenOpenThirdDerivative(lengths,values,h);
    if(end==Endpoint::FixedEven) {
        auto result=openThirdDerivative(std::vector<double>(lengths.rbegin(),lengths.rend()),
            std::vector<double>(values.rbegin(),values.rend()),h,end,begin);
        std::reverse(result.values.begin(),result.values.end());return result;
    }
    const int count=static_cast<int>(values.size());const bool even=begin==Endpoint::FixedEven;
    const int n=even?2*count-1:count;
    std::vector<double> field(n),edges(n-1);std::vector<int> columns(n);
    for(int i=0;i<n;++i){const int index=even?std::abs(i-(count-1)):i;field[i]=values[index];
        columns[i]=(index==0 || index==count-1)?-1:index-1;
        if(i+1<n)edges[i]=lengths[even?(i<count-1?count-2-i:i-(count-1)):i];}
    auto result=solve(edges,field,h,columns,count-2,false);
    if(even){result.values.erase(result.values.begin(),result.values.begin()+count-1);
        result.dataEnergy/=2;result.derivativeEnergy/=2;}
    return result;
}
}
