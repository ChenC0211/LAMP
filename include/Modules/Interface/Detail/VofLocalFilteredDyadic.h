#pragma once
// VOF-only implementation detail. No public ABI and no generic arithmetic API.
#include "Core/Geometry/Detail/FiniteVolumeClosure.h"
#include <cstdint>
#include <cstring>
#include <optional>
#include <type_traits>
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
#include <xmmintrin.h>
#endif
namespace LAMP::Interface::Detail {
namespace {
namespace VofFilteredFV=LAMP::Geometry::Detail::FiniteVolume;
constexpr bool vofFilteredNativeFormat=sizeof(double)==8&&sizeof(long double)==8
    &&std::numeric_limits<double>::is_iec559&&std::numeric_limits<long double>::radix==2
    &&std::numeric_limits<long double>::digits==53&&std::numeric_limits<long double>::min_exponent==-1021
    &&std::numeric_limits<long double>::max_exponent==1024;
enum class VofLocalPromotion {None,ThirdTail,NonNormal,Unsupported,ExactOperand};
enum class VofLocalPhaseUpdatePath {Noop,Fast,Promoted,Fallback};
inline std::uint64_t vofFilteredBits(double value) noexcept
{std::uint64_t result=0;std::memcpy(&result,&value,sizeof result);return result;}
inline bool vofFilteredNormalOrZero(double value) noexcept
{const auto bits=vofFilteredBits(value)&UINT64_C(0x7fffffffffffffff);
 return bits==0||(bits>=UINT64_C(0x0010000000000000)&&bits<UINT64_C(0x7ff0000000000000));}
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE) \
    && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST)
#pragma float_control(precise,on,push)
#pragma float_control(except,on)
#pragma fenv_access(on)
struct VofFilteredPair {double head=0,tail=0;};
// Only volatile SSE binary64 add/sub and integer classification occur here.
// The caller restores raw MXCSR even when an optional intermediate overflows.
__declspec(noinline) inline bool vofFilteredTwoSum(double a,double b,VofFilteredPair& result) noexcept
{
    volatile double head=a+b;if(!vofFilteredNormalOrZero(head))return false;
    volatile double back=head-a;if(!vofFilteredNormalOrZero(back))return false;
    volatile double inner=head-back;if(!vofFilteredNormalOrZero(inner))return false;
    volatile double left=a-inner;if(!vofFilteredNormalOrZero(left))return false;
    volatile double right=b-back;if(!vofFilteredNormalOrZero(right))return false;
    volatile double tail=left+right;if(!vofFilteredNormalOrZero(tail))return false;
    result.head=head;result.tail=tail;return true;
}
__declspec(noinline) inline bool vofFilteredTryAdd(double h,double l,double x,
    double& nextH,double& nextL,VofLocalPromotion& reason) noexcept
{
    const unsigned csr=_mm_getcsr();
    // RNE, masked exceptions, gradual underflow. Existing sticky flags remain.
    if((csr&0x1f80u)!=0x1f80u||(csr&0xe040u)!=0){reason=VofLocalPromotion::Unsupported;return false;}
    if(!vofFilteredNormalOrZero(h)||!vofFilteredNormalOrZero(l)||!vofFilteredNormalOrZero(x)){
        reason=VofLocalPromotion::NonNormal;return false;}
    VofFilteredPair first,second,last;
    const bool firstOk=vofFilteredTwoSum(h,x,first);
    const bool secondOk=firstOk&&vofFilteredTwoSum(l,first.tail,second);
    const bool third=secondOk&&second.tail!=0;
    const bool finalOk=secondOk&&!third&&vofFilteredTwoSum(first.head,second.head,last);
    _mm_setcsr(csr); // No CRT, errno, x87, callback, or old scalar work is shielded.
    if(!finalOk){reason=third?VofLocalPromotion::ThirdTail:VofLocalPromotion::NonNormal;return false;}
    nextH=last.head==0?0.:last.head;nextL=last.tail==0?0.:last.tail;return true;
}
#if !defined(_M_FP_STRICT) || !_M_FP_STRICT
#pragma fenv_access(off)
#endif
#pragma float_control(pop)
#endif
class VofLocalFilteredDyadic
{
    double high_=0,low_=0;
    std::optional<VofFilteredFV::Dyadic> exact_;
    VofLocalPromotion promotion_=VofLocalPromotion::None;
    void promote(VofLocalPromotion why)
    {if(exact_)return;exact_.emplace(high_);if(low_!=0)exact_->add(low_);promotion_=why;}
public:
    VofLocalFilteredDyadic()=default;
    explicit VofLocalFilteredDyadic(long double value){add(value);}
    bool isFast() const noexcept{return !exact_;}
    VofLocalPromotion promotion() const noexcept{return promotion_;}
    double high() const noexcept{return high_;}
    double low() const noexcept{return low_;}
    void add(long double value,int sign=1)
    {
        // Preserve native-overload validation order for invalid values/signs.
        if(sign!=1&&sign!=-1){promote(VofLocalPromotion::ExactOperand);exact_->add(value,sign);return;}
        if(exact_){exact_->add(value,sign);return;}
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE) \
    && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST)
        if constexpr(sizeof(long double)==8&&std::numeric_limits<long double>::digits==53
            &&std::numeric_limits<long double>::radix==2&&std::numeric_limits<long double>::min_exponent==-1021
            &&std::numeric_limits<long double>::max_exponent==1024&&std::numeric_limits<double>::is_iec559){
            double x=static_cast<double>(value);if(sign<0){auto bits=vofFilteredBits(x)^UINT64_C(0x8000000000000000);std::memcpy(&x,&bits,sizeof x);}
            double h=0,l=0;VofLocalPromotion why=VofLocalPromotion::None;
            if(vofFilteredTryAdd(high_,low_,x,h,l,why)){high_=h;low_=l;return;}
            promote(why);exact_->add(value,sign);return;
        }
#endif
        promote(VofLocalPromotion::Unsupported);exact_->add(value,sign);
    }
    void add(const VofLocalFilteredDyadic& value,int sign=1)
    {
        VofFilteredFV::Internal::require(sign==1||sign==-1,"[FiniteVolume] invalid dyadic addition sign");
        if(this==&value){const auto copy=value;add(copy,sign);return;}
        if(value.exact_){promote(value.promotion_);exact_->add(*value.exact_,sign);return;}
        const double h=value.high_,l=value.low_;add(h,sign);if(l!=0)add(l,sign);
    }
    VofFilteredFV::Dyadic dyadic() const
    {if(exact_)return *exact_;VofFilteredFV::Dyadic result(high_);if(low_!=0)result.add(low_);return result;}
    int compare(const VofLocalFilteredDyadic& value) const
    {
        if(exact_||value.exact_)return dyadic().compare(value.dyadic());
        if(high_!=value.high_)return high_<value.high_?-1:1;
        return low_==value.low_?0:(low_<value.low_?-1:1);
    }
    int compare(const VofFilteredFV::Dyadic& value) const{return dyadic().compare(value);}
    long double rounded() const
    {
        // Preserve the original environment and native scalbn suffix. Avoid a
        // changed exact-top/quantum at binade carry, and nonnoop exact zero.
        if(exact_)return exact_->rounded();
        if constexpr(!vofFilteredNativeFormat)return dyadic().rounded();
        const auto magnitude=vofFilteredBits(high_)&UINT64_C(0x7fffffffffffffff);
        const bool opposite=low_!=0&&((vofFilteredBits(high_)^vofFilteredBits(low_))>>63)!=0;
        if(magnitude==0||((magnitude&UINT64_C(0x000fffffffffffff))==0&&opposite))return dyadic().rounded();
        return VofFilteredFV::Dyadic(high_).rounded();
    }
};
inline VofFilteredFV::Rounded vofLocalPhaseUpdate(long double oldStock,const VofLocalFilteredDyadic& net,
    const VofLocalFilteredDyadic& compression,VofLocalPhaseUpdatePath& path)
{
    VofFilteredFV::Internal::require(std::isfinite(oldStock)&&oldStock>=0,"[FiniteVolume] invalid source phase stock");
    const VofLocalFilteredDyadic zero;
    if(net.compare(zero)==0&&compression.compare(zero)==0){path=VofLocalPhaseUpdatePath::Noop;return {oldStock,{}};}
    if constexpr(!vofFilteredNativeFormat){path=VofLocalPhaseUpdatePath::Fallback;
        return VofFilteredFV::phaseUpdate(oldStock,net.dyadic(),compression.dyadic());}
    VofLocalFilteredDyadic result(oldStock);result.add(net,-1);result.add(compression);
    const auto exactFallback=[&]{return VofFilteredFV::phaseUpdate(oldStock,net.dyadic(),compression.dyadic());};
    if(!result.isFast()){
        path=result.promotion()==VofLocalPromotion::ThirdTail||result.promotion()==VofLocalPromotion::ExactOperand
            ?VofLocalPhaseUpdatePath::Promoted:VofLocalPhaseUpdatePath::Fallback;return exactFallback();}
    const auto magnitude=vofFilteredBits(result.high())&UINT64_C(0x7fffffffffffffff);
    const bool opposite=result.low()!=0&&((vofFilteredBits(result.high())^vofFilteredBits(result.low()))>>63)!=0;
    if(magnitude==0||((magnitude&UINT64_C(0x000fffffffffffff))==0&&opposite)){
        path=VofLocalPhaseUpdatePath::Fallback;return exactFallback();}
    path=VofLocalPhaseUpdatePath::Fast;
    VofFilteredFV::Rounded answer;answer.value=VofFilteredFV::Dyadic(result.high()).rounded();
    if(result.low()!=0)answer.error=VofFilteredFV::Dyadic(-result.low());return answer;
}
} // TU-local implementation: production and the independent oracle use this same header.
} // namespace LAMP::Interface::Detail
