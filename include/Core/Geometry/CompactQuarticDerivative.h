#pragma once

#include <algorithm>
#include <array>
#include <cmath>

namespace Quantum::Geometry
{
// Differentiate the degree-four interpolant on five consecutive skeleton
// samples, centered where possible and one-sided at an endpoint. Actual
// physical abscissae are shifted/scaled; subtracting the target value gives
// an exactly zero derivative for constant data. Check must throw on false
// and supplies the owning caller's error context.
// Extracted without changing the arithmetic/order of the existing GALS kernel.
template<class Coordinate,class Value,class Check>
double compactQuarticDerivativeAtNode(int target,int count,const Coordinate& coordinate,
    const Value& value,const Check& require)
{
    require(count>=5 && target>=0 && target<count,"insufficient compact derivative support");
    const int begin=std::clamp(target-2,0,count-5), center=target-begin;
    const long double origin=coordinate(target);
    const long double scale=static_cast<long double>(coordinate(begin+4))-coordinate(begin);
    require(std::isfinite(scale) && scale>0,"unresolved skeleton derivative scale");
    std::array<long double,5> t{};
    for(int j=0;j<5;++j) {
        t[j]=(static_cast<long double>(coordinate(begin+j))-origin)/scale;
        require(std::isfinite(t[j]) && (j==0 || t[j]>t[j-1]),"non-increasing skeleton abscissae");
    }
    require(t[center]==0,"compact derivative lost its target node");
    const long double reference=value(target);long double sum=0;
    for(int j=0;j<5;++j) if(j!=center) {
        long double numerator=1,denominator=1;
        for(int k=0;k<5;++k) if(k!=j) {
            denominator*=t[j]-t[k];
            if(k!=center) numerator*=-t[k];
        }
        require(denominator!=0 && std::isfinite(denominator),"singular compact derivative support");
        sum+=(numerator/denominator)*(static_cast<long double>(value(begin+j))-reference);
    }
    const double result=static_cast<double>(sum/scale);
    require(std::isfinite(result),"nonfinite compact mixed derivative");
    return result;
}
} // namespace Quantum::Geometry
