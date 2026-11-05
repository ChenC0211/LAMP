#pragma once
// Private derived evaluator: same rational linear polynomial and native rounding.
#include "Core/Geometry/Detail/QuadraticGraphExact.h"
#include <algorithm>
#include <cmath>
#include <cfenv>
#include <limits>
#include <optional>
#if defined(_M_X64) || defined(_M_IX86) || defined(__SSE__)
#include <xmmintrin.h>
#endif

namespace LAMP::Geometry::Detail::ExactLinearNative
{
// A filter of the ORIGINAL rational polynomial, never a replacement curve.
// All uncertain rounding cells return to QuadraticExact::rounded below.
#ifdef _MSC_VER
#pragma float_control(precise,on,push)
#endif
struct NativePair {long double high=0,low=0;};
inline NativePair twoSum(long double a,long double b)
{
    const auto high=a+b,virtualB=high-a;
    return {high,(a-(high-virtualB))+(b-virtualB)};
}
inline bool slopeFilterEnvironment()
{
#if defined(__FAST_MATH__) || defined(_M_FP_FAST)
    return false;
#else
    if(!std::numeric_limits<long double>::is_iec559||std::numeric_limits<long double>::radix!=2||
        std::numeric_limits<long double>::has_denorm!=std::denorm_present||std::fegetround()!=FE_TONEAREST)return false;
#if defined(_M_X64) || defined(_M_IX86) || defined(__SSE__)
    if((_mm_getcsr()&0x8040u)!=0)return false; // FTZ/DAZ break gradual-underflow bounds.
#endif
    return true;
#endif
}
struct SlopeCoefficient {long double high=0,low=0,error=0;};
inline std::optional<SlopeCoefficient> packSlopeCoefficient(const QuadraticExact::Scalar& value)
{
    namespace E=QuadraticExact;
    static const E::Scalar minimum=E::scalar(std::numeric_limits<long double>::min()),
        maximum=E::scalar(std::numeric_limits<long double>::max());
    const auto inNormalRange=[](const E::Scalar& x){const auto magnitude=x<0?-x:x;
        return x==0||(magnitude>=minimum&&magnitude<=maximum);};
    if(!inNormalRange(value))return {};
    SlopeCoefficient result;result.high=E::rounded(value);
    const E::Scalar remainder=value-E::scalar(result.high);
    if(!inNormalRange(remainder))return {};
    result.low=E::rounded(remainder);
    E::Scalar remaining=remainder-E::scalar(result.low);if(remaining<0)remaining=-remaining;
    if(!inNormalRange(remaining))return {};
    result.error=E::rounded(remaining);
    if(E::scalar(result.error)<remaining)result.error=std::nextafter(result.error,std::numeric_limits<long double>::infinity());
    if(!std::isfinite(result.error))return {};
    return result;
}
inline long double upwardSum(long double a,long double b)
{
    if(a==0)return b;if(b==0)return a;
    return std::nextafter(a+b,std::numeric_limits<long double>::infinity());
}
inline long double upwardProduct(long double a,long double b)
{
    if(a==0||b==0)return 0;
    return std::nextafter(a*b,std::numeric_limits<long double>::infinity());
}
inline long double nativeRoundBound(long double value)
{
    const auto magnitude=std::abs(value);
    return std::nextafter(magnitude,std::numeric_limits<long double>::infinity())-magnitude;
}
struct Filter
{
    SlopeCoefficient constant,linear;
    bool available=false;
    static Filter prepare(const QuadraticExact::Scalar& constantValue,const QuadraticExact::Scalar& linearValue)
    {
        if(!slopeFilterEnvironment())return {};
        try{
            const auto a=packSlopeCoefficient(constantValue),b=packSlopeCoefficient(linearValue);
            if(a&&b)return {*a,*b,true};
        }catch(const std::invalid_argument&){
            // Packing must not reject a graph whose actual slope may be valid
            // after cancellation. Keep the original pointwise exact query.
        }
        return {};
    }
    bool evaluate(long double t,long double& result) const
    {
        const auto normalOrZero=[](long double x){return x==0||std::isnormal(x);};
        if(!available||!slopeFilterEnvironment()||!normalOrZero(t))return false;
        const auto product=linear.high*t;
        if(!normalOrZero(product)||(product==0&&linear.high!=0&&t!=0))return false;
        const auto leading=twoSum(constant.high,product);
        if(!normalOrZero(leading.high)||!normalOrZero(leading.low))return false;
        if(constant.high!=0&&product!=0&&((constant.high<0)!=(product<0))&&
            std::abs(leading.high)<.25L*std::max(std::abs(constant.high),std::abs(product)))return false;
        const auto productLow=std::fma(linear.high,t,-product),low=std::fma(linear.low,t,constant.low);
        const auto tail=leading.low+productLow,tailWithLow=tail+low;
        if(!normalOrZero(productLow)||!normalOrZero(low)||!normalOrZero(tail)||!normalOrZero(tailWithLow))return false;
        const auto rounded=twoSum(leading.high,tailWithLow);
        if(!std::isnormal(rounded.high)||!normalOrZero(rounded.low))return false;
        // TwoSum terms are exact. Bound the four rounded low-order operations
        // and the once-packed coefficient remainder, all outward in native.
        auto error=upwardSum(constant.error,upwardProduct(linear.error,std::abs(t)));
        for(const auto value:{productLow,low,tail,tailWithLow})error=upwardSum(error,nativeRoundBound(value));
        if(!std::isfinite(error))return false;
        const auto previous=std::nextafter(rounded.high,-std::numeric_limits<long double>::infinity()),
            next=std::nextafter(rounded.high,std::numeric_limits<long double>::infinity());
        if(!std::isfinite(previous)||!std::isfinite(next))return false;
        const auto lowerHalf=(rounded.high-previous)/2,upperHalf=(next-rounded.high)/2;
        if(!std::isnormal(lowerHalf)||!std::isnormal(upperHalf))return false;
        const auto lower=std::nextafter(rounded.low-error,-std::numeric_limits<long double>::infinity()),
            upper=std::nextafter(rounded.low+error,std::numeric_limits<long double>::infinity());
        if(!(lower>-lowerHalf&&upper<upperHalf))return false;
        result=rounded.high;return true;
    }
};
#ifdef _MSC_VER
#pragma float_control(pop)
#endif
}
