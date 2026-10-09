#pragma once

#include "Core/Geometry/Detail/ExactExpansion.h"
#include <array>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>

namespace Quantum::Geometry::Detail
{
namespace
{
// Algebraic certificate on the SAME exact coefficients. No approximate GCD,
// division by a rounded coefficient, coefficient deletion or root merging.
namespace PolynomialSquareFree
{
namespace E=ExactExpansion;
template<std::size_t N> using Polynomial=std::array<E::Expansion,N>;
inline void require(bool ok,const char* reason)
{if(!ok)throw std::runtime_error(std::string("[PolynomialSquareFree] ")+reason);}
template<std::size_t N> int degree(const Polynomial<N>& p)
{int d=static_cast<int>(N)-1;while(d>=0 && E::sign(p[d])==0)--d;return d;}
template<std::size_t N> void normalize(Polynomial<N>& p)
{E::ExponentRange range;for(const auto& c:p)range.include(c);const int shift=range.centeredShift();
    for(auto& c:p)c=E::scalePow2(c,shift);}
template<std::size_t N> Polynomial<N> derivative(const Polynomial<N>& p)
{Polynomial<N> d{};for(std::size_t k=1;k<N;++k)d[k-1]=E::multiply(p[k],E::fromDouble(static_cast<double>(k)));return d;}
template<std::size_t N> Polynomial<N> remainder(Polynomial<N> a,Polynomial<N> b)
{
    normalize(a);normalize(b);const int db=degree(b);require(db>=0,"zero divisor");
    for(int da=degree(a);da>=db;da=degree(a)){
        const auto leading=a[da];const int shift=da-db;
        for(auto& c:a)c=E::multiply(c,b[db]);
        for(int j=0;j<=db;++j)a[j+shift]=E::subtract(a[j+shift],E::multiply(leading,b[j]));
        require(degree(a)<da,"pseudo-remainder failed to reduce degree");normalize(a);
    }return a;
}
template<std::size_t N> Polynomial<N> gcd(Polynomial<N> a,Polynomial<N> b)
{
    normalize(a);normalize(b);
    while(degree(b)>=0){if(degree(b)==0)return b;auto r=remainder(a,b);a=std::move(b);b=std::move(r);}
    return a;
}
template<std::size_t N> Polynomial<N> exactQuotient(const Polynomial<N>& input,Polynomial<N> divisor)
{
    auto r=input;Polynomial<N> q{};normalize(r);normalize(divisor);const int d=degree(divisor);require(d>=0,"zero quotient divisor");
    for(int dr=degree(r);dr>=d;dr=degree(r)){
        const auto leading=r[dr];const int shift=dr-d;
        for(std::size_t k=0;k<N;++k){r[k]=E::multiply(r[k],divisor[d]);q[k]=E::multiply(q[k],divisor[d]);}
        for(int j=0;j<=d;++j)r[j+shift]=E::subtract(r[j+shift],E::multiply(leading,divisor[j]));
        q[shift]=E::add(q[shift],leading);require(degree(r)<dr,"pseudo-quotient failed to reduce degree");
        // Q and remainder must share ONE scaling, preserving S*P=Q*G+R.
        E::ExponentRange range;for(std::size_t k=0;k<N;++k){range.include(q[k]);range.include(r[k]);}
        const int exponent=range.centeredShift();for(std::size_t k=0;k<N;++k){q[k]=E::scalePow2(q[k],exponent);r[k]=E::scalePow2(r[k],exponent);}
    }
    require(degree(r)<0,"claimed factor leaves a nonzero exact remainder");
    require(degree(q)+d==degree(input),"quotient degree mismatch");
    Polynomial<N> product{};for(int i=0;i<=degree(q);++i)for(int j=0;j<=d;++j)product[i+j]=E::add(product[i+j],E::multiply(q[i],divisor[j]));
    auto p=input;normalize(p);normalize(product);int pivot=0;while(pivot<static_cast<int>(N) && E::sign(p[pivot])==0)++pivot;
    require(pivot<static_cast<int>(N) && E::sign(product[pivot])!=0,"zero recomposition pivot");
    for(std::size_t k=0;k<N;++k)require(E::sign(E::subtract(E::multiply(p[k],product[pivot]),E::multiply(product[k],p[pivot])))==0,
        "quotient recomposition is not exactly proportional to original polynomial");
    normalize(q);return q;
}
template<std::size_t N> struct Result {Polynomial<N> polynomial;bool deflated=false;};
template<std::size_t N> Result<N> squareFree(const Polynomial<N>& p)
{
    static_assert(N<=9,"bounded native Q4 restriction only");E::validateEnvironment();
    if(degree(p)<2)return {p,false};auto normalized=p;normalize(normalized);
    const auto common=gcd(normalized,derivative(normalized));
    if(degree(common)<=0)return {p,false};
    return {exactQuotient(normalized,common),true};
}
} // namespace PolynomialSquareFree
} // TU-local exact arithmetic
}
