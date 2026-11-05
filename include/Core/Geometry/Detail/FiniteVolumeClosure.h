#pragma once
#include "Core/Geometry/PlanarCellCut.h"
#include "Core/Mesh/CartesianMeasure.h"
#include <boost/multiprecision/cpp_int.hpp>
#include <algorithm>
#include <array>
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
#include <intrin.h>
#endif
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <initializer_list>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace LAMP::Geometry::Detail::FiniteVolume
{
// Exact binary arithmetic for finite-volume closure receipts. No physical
// inventory is stored here. A receipt is publishable only by typed owners.
namespace Internal
{
inline void require(bool condition,const char* message)
{ if(!condition)throw std::invalid_argument(message); }
// Unsigned facade over the installed multiprecision integer. Dyadic sums
// never use a rational denominator or gcd; integer division is CDF-only.
struct UInt
{
    // MinBits changes only inline storage. MaxBits=0 retains the allocator-backed
    // unbounded signed-magnitude integer and the original exact arithmetic.
    using Integer=boost::multiprecision::number<boost::multiprecision::cpp_int_backend<
        256,0,boost::multiprecision::signed_magnitude,boost::multiprecision::unchecked>>;
    Integer n=0;
    UInt()=default;
    explicit UInt(std::uint64_t value):n(value){}
    explicit UInt(Integer value):n(std::move(value)){}
    bool zero() const{return n==0;}
    int bits() const{return zero()?0:static_cast<int>(boost::multiprecision::msb(n))+1;}
    bool odd() const{return (n&1)!=0;}
    int compare(const UInt& b) const{return n<b.n?-1:n>b.n?1:0;}
    UInt shifted(int shift) const
    {require(shift>=0,"[FiniteVolume] negative integer left shift");return UInt(Integer(n<<shift));}
    UInt right(int shift) const
    {require(shift>=0,"[FiniteVolume] negative integer right shift");return UInt(Integer(n>>shift));}
    void add(const UInt& b){n+=b.n;}
    void subtract(const UInt& b)
    {require(n>=b.n,"[FiniteVolume] negative unsigned difference");n-=b.n;}
    UInt multiplied(const UInt& b) const{return UInt(Integer(n*b.n));}
    int trailing() const{return zero()?0:static_cast<int>(boost::multiprecision::lsb(n));}
};
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
// Read-only terminal RNE for a nonnegative exact magnitude. The caller's
// original top/quantum selection proves the rounded integer is <=2^53,
// and a nonpositive shift has magnitude <=52. No UInt is shifted or copied.
template<class WordAt>
inline std::uint64_t nativeTerminalRneWords(const WordAt& word,std::size_t size,int shift) noexcept
{
    if(shift<=0)return word(0)<<static_cast<unsigned>(-shift);
    const auto first=static_cast<std::size_t>(shift/64);
    const auto bits=static_cast<unsigned>(shift%64);
    auto integer=word(first)>>bits;
    if(bits)integer|=word(first+1)<<(64-bits);
    const auto halfWord=static_cast<std::size_t>((shift-1)/64);
    const auto halfBit=static_cast<unsigned>((shift-1)%64);
    const auto halfMask=UINT64_C(1)<<halfBit;
    if((word(halfWord)&halfMask)==0)return integer;
    bool sticky=(word(halfWord)&(halfMask-1))!=0;
    for(std::size_t k=0;k<halfWord&&k<size&&!sticky;++k)sticky=word(k)!=0;
    if(sticky||(integer&1u))++integer;
    return integer;
}
inline std::uint64_t nativeTerminalRneInteger(const UInt& magnitude,int shift) noexcept
{
    const auto* words=magnitude.n.backend().limbs();
    static constexpr auto limbBytes=sizeof(*words);
    static_assert(limbBytes==4||limbBytes==8,
        "native terminal RNE requires complete 32-bit or 64-bit unsigned limbs");
    const auto limbCount=static_cast<std::size_t>(magnitude.n.backend().size());
    const auto size=limbBytes==8?limbCount:limbCount/2+limbCount%2;
    const auto word=[&](std::size_t at) noexcept ->std::uint64_t {
        if constexpr(limbBytes==8)
            return at<limbCount?static_cast<std::uint64_t>(words[at]):UINT64_C(0);
        else {
            const auto first=at*2;
            const auto lower=first<limbCount?static_cast<std::uint64_t>(words[first]):UINT64_C(0);
            const auto upper=first+1<limbCount?static_cast<std::uint64_t>(words[first+1]):UINT64_C(0);
            return lower|(upper<<32);
        }
    };
    return nativeTerminalRneWords(word,size,shift);
}
#endif
inline UInt nearestInteger(const UInt& n,const UInt& d)
{
    require(!d.zero(),"[FiniteVolume] zero integer divisor");
    UInt::Integer q=n.n/d.n,remainder=n.n%d.n;
    const UInt::Integer twice=remainder<<1;
    if(twice>d.n||(twice==d.n&&(q&1)!=0))++q;
    return UInt(std::move(q));
}
inline void environment()
{
    require(std::numeric_limits<long double>::radix==2&&std::numeric_limits<long double>::has_denorm==std::denorm_present,
        "[FiniteVolume] native binary gradual-underflow representation required");
    require(std::fegetround()==FE_TONEAREST,"[FiniteVolume] native round-to-nearest-even required");
    volatile long double tiny=std::numeric_limits<long double>::denorm_min();volatile long double one=1;
    require(tiny*one>0,"[FiniteVolume] gradual underflow required");
}
// This port stores binary64, independently of the native-long-double format.
inline void binary64Environment()
{
    require(std::numeric_limits<double>::radix==2&&std::numeric_limits<double>::digits==53&&
        std::numeric_limits<double>::min_exponent==-1021&&std::numeric_limits<double>::max_exponent==1024&&
        std::numeric_limits<double>::has_denorm==std::denorm_present,
        "[FiniteVolume] binary64 gradual-underflow representation required");
    require(std::fegetround()==FE_TONEAREST,"[FiniteVolume] binary64 round-to-nearest-even required");
    volatile double tiny=std::numeric_limits<double>::denorm_min();volatile double one=1.;
    require(tiny*one>0,"[FiniteVolume] binary64 gradual underflow required");
}
// The original scalar guards run before this optional decomposition. Only
// normal native binary64 values use bits; every other input stays literal.
inline bool decodeNormalBinary64(long double value,std::uint64_t& significand,int& exponent)
{
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE) && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST)
    if constexpr(sizeof(long double)==sizeof(std::uint64_t)&&std::numeric_limits<long double>::is_iec559
        &&std::numeric_limits<long double>::radix==2&&std::numeric_limits<long double>::digits==53
        &&std::numeric_limits<long double>::min_exponent==-1021&&std::numeric_limits<long double>::max_exponent==1024
        &&std::numeric_limits<long double>::has_denorm==std::denorm_present){
        std::uint64_t bits=0;std::memcpy(&bits,&value,sizeof(bits));
        const auto field=static_cast<unsigned>((bits>>52)&0x7ffu);
        if(field!=0&&field!=0x7ffu){
            significand=(bits&UINT64_C(0x000fffffffffffff))|UINT64_C(0x0010000000000000);
            exponent=static_cast<int>(field)-1075;return true;
        }
    }
#endif
    (void)value;(void)significand;(void)exponent;return false;
}
}
class Dyadic;
class Binary64ExactProducts;
struct Partition;
Partition partition(long double,const std::vector<Dyadic>&);
struct PartitionDouble;
PartitionDouble partitionDouble(double,const std::vector<Dyadic>&);
class Dyadic
{
public:
    Dyadic()=default;
    explicit Dyadic(long double value)
    {
        Internal::require(std::isfinite(value),"[FiniteVolume] nonfinite dyadic input");
        Internal::require(std::numeric_limits<long double>::radix==2,"[FiniteVolume] binary native arithmetic required");
        if(value==0)return;negative_=std::signbit(value);
        std::uint64_t significand=0;int exponent=0;
        if(Internal::decodeNormalBinary64(value,significand,exponent)){
            magnitude_=Internal::UInt(significand);exponent_=exponent;normalize();return;
        }
        int e=0;auto f=std::frexp(std::abs(value),&e);
        constexpr int p=std::numeric_limits<long double>::digits;exponent_=e-p;
        if constexpr(p<=64){const auto integer=std::scalbn(f,p);
            Internal::require(integer>=1&&integer<std::scalbn(1.L,64),"[FiniteVolume] native significand outside uint64");
            magnitude_=Internal::UInt(static_cast<std::uint64_t>(integer));
        }else{int remaining=p;while(remaining){const int take=std::min(remaining,32);const auto expanded=std::scalbn(f,take);
            const auto chunk=static_cast<std::uint32_t>(expanded);magnitude_=magnitude_.shifted(take);magnitude_.add(Internal::UInt(chunk));
            f=expanded-chunk;remaining-=take;}}
        normalize();
    }
    void add(long double value,int sign=1){add(Dyadic(value),sign);}
    void add(const Dyadic& value,int sign=1)
    {
        Internal::require(sign==1||sign==-1,"[FiniteVolume] invalid dyadic addition sign");if(value.magnitude_.zero())return;
        const bool neg=value.negative_!=(sign<0);
        if(magnitude_.zero()){Dyadic result=value;result.negative_=neg;*this=std::move(result);return;}
        const bool leftHigher=exponent_>=value.exponent_;
        const auto& higher=leftHigher?*this:value;const auto& lower=leftHigher?value:*this;
        const bool higherNegative=leftHigher?negative_:neg,lowerNegative=leftHigher?neg:negative_;
        Dyadic result;result.magnitude_=higher.magnitude_.shifted(higher.exponent_-lower.exponent_);
        result.exponent_=lower.exponent_;result.negative_=higherNegative;
        if(higherNegative==lowerNegative)result.magnitude_.add(lower.magnitude_);
        else{
            // Only this unpublished local integer may temporarily be signed.
            result.magnitude_.n-=lower.magnitude_.n;
            if(result.magnitude_.n<0){result.magnitude_.n=-result.magnitude_.n;result.negative_=lowerNegative;}
        }
        result.normalize();*this=std::move(result);
    }
    static Dyadic product(std::initializer_list<long double> factors)
    {
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE) && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST)
        if constexpr(sizeof(long double)==8&&sizeof(double)==8
            &&std::numeric_limits<long double>::radix==2&&std::numeric_limits<long double>::digits==53
            &&std::numeric_limits<long double>::min_exponent==-1021&&std::numeric_limits<long double>::max_exponent==1024
            &&std::numeric_limits<long double>::has_denorm==std::denorm_present){
            if(factors.size()>=2&&factors.size()<=4)return smallNativeProduct(factors);
        }
#endif
        Dyadic result(1);for(const auto f:factors){const Dyadic x(f);if(x.magnitude_.zero())return {};
            result.negative_=result.negative_!=x.negative_;result.exponent_+=x.exponent_;result.magnitude_=result.magnitude_.multiplied(x.magnitude_);}
        result.normalize();return result;
    }
    // Exact native-factor scaling of a retained contraction, without
    // rounding either input or modifying this value.
    Dyadic multiplied(long double factor) const
    {
        const Dyadic scale(factor);
        if(magnitude_.zero()||scale.magnitude_.zero())return {};
        Dyadic result=*this;
        result.negative_=result.negative_!=scale.negative_;
        result.exponent_+=scale.exponent_;
        result.magnitude_=result.magnitude_.multiplied(scale.magnitude_);
        result.normalize();return result;
    }
    // Fixed minimum binary64 charge quantum, WITHOUT a 53-bit significand
    // restriction. An eta-integral value retains every bit and any exact range.
    // Only a fractional eta is rounded, with symmetric nearest-even ties.
    Dyadic roundedBinary64Quantum() const
    {
        Internal::binary64Environment();
        constexpr int quantum=std::numeric_limits<double>::min_exponent-std::numeric_limits<double>::digits;
        if(magnitude_.zero()||exponent_>=quantum)return *this;
        const int shift=quantum-exponent_;auto integer=magnitude_.right(shift);
        auto remainder=magnitude_;remainder.subtract(integer.shifted(shift));
        const int half=remainder.compare(Internal::UInt(1).shifted(shift-1));
        if(half>0||(half==0&&integer.odd()))integer.add(Internal::UInt(1));
        Dyadic result;result.magnitude_=std::move(integer);result.exponent_=quantum;
        result.negative_=negative_;result.normalize();return result;
    }
    int compare(const Dyadic& b) const
    {
        if(magnitude_.zero())return b.magnitude_.zero()?0:b.negative_?1:-1;
        if(b.magnitude_.zero())return negative_?-1:1;
        if(negative_!=b.negative_)return negative_?-1:1;
        int result=0;
        if(exponent_==b.exponent_)result=magnitude_.compare(b.magnitude_);
        else{
            const int bits=magnitude_.bits(),otherBits=b.magnitude_.bits();
            const auto top=static_cast<std::int64_t>(exponent_)+bits;
            const auto otherTop=static_cast<std::int64_t>(b.exponent_)+otherBits;
            if(top!=otherTop)result=top<otherTop?-1:1;
            else if(bits<otherBits)result=magnitude_.shifted(otherBits-bits).compare(b.magnitude_);
            else result=-b.magnitude_.shifted(bits-otherBits).compare(magnitude_);
        }
        return negative_?-result:result;
    }
    Dyadic absolute() const{auto r=*this;r.negative_=false;return r;}
    long double rounded() const
    {
        Internal::environment();if(magnitude_.zero())return 0;
        constexpr int p=std::numeric_limits<long double>::digits;
        constexpr int least=std::numeric_limits<long double>::min_exponent-p;
        const int top=magnitude_.bits()-1+exponent_;
        Internal::require(top<std::numeric_limits<long double>::max_exponent,"[FiniteVolume] final native dyadic overflow");
        const int quantum=std::max(top-(p-1),least);const int shift=quantum-exponent_;
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
        if constexpr(std::numeric_limits<long double>::digits==53){
            const auto integer=Internal::nativeTerminalRneInteger(magnitude_,shift);
            if(integer==0)return std::copysign(0.L,negative_?-1.L:1.L);
            const long double significand=static_cast<long double>(integer);
            const auto value=std::scalbn(significand,quantum);
            Internal::require(std::isfinite(value),"[FiniteVolume] final native dyadic overflow");
            Internal::require(value>0,"[FiniteVolume] native dyadic conversion lost positive quantum");return negative_?-value:value;
        }
#endif
        Internal::UInt integer;
        if(shift>0){integer=magnitude_.right(shift);auto remainder=magnitude_;remainder.subtract(integer.shifted(shift));
            const int half=remainder.compare(Internal::UInt(1).shifted(shift-1));if(half>0||(half==0&&integer.odd()))integer.add(Internal::UInt(1));
        }else integer=magnitude_.shifted(-shift);
        if(integer.zero())return std::copysign(0.L,negative_?-1.L:1.L);
        // The rounded integer has <=p bits, or is exactly the carry 2^p;
        // conversion is exact in the native significand before power scaling.
        const long double significand=integer.n.convert_to<long double>();
        const auto value=std::scalbn(significand,quantum);
        Internal::require(std::isfinite(value),"[FiniteVolume] final native dyadic overflow");
        Internal::require(value>0,"[FiniteVolume] native dyadic conversion lost positive quantum");return negative_?-value:value;
    }
    // Direct binary64 RNE, never native-LD rounding followed by a cast.
    double roundedDouble() const
    {
        Internal::binary64Environment();if(magnitude_.zero())return 0;
        constexpr int p=std::numeric_limits<double>::digits;
        constexpr int least=std::numeric_limits<double>::min_exponent-p;
        const int top=magnitude_.bits()-1+exponent_;
        Internal::require(top<std::numeric_limits<double>::max_exponent,"[FiniteVolume] final binary64 dyadic overflow");
        const int quantum=std::max(top-(p-1),least);const int shift=quantum-exponent_;
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
        if constexpr(std::numeric_limits<double>::digits==53){
            const auto integer=Internal::nativeTerminalRneInteger(magnitude_,shift);
            if(integer==0)return std::copysign(0.,negative_?-1.:1.);
            const double significand=static_cast<double>(integer);
            const auto value=std::scalbn(significand,quantum);
            Internal::require(std::isfinite(value),"[FiniteVolume] final binary64 dyadic overflow");
            Internal::require(value>0,"[FiniteVolume] binary64 dyadic conversion lost positive quantum");return negative_?-value:value;
        }
#endif
        Internal::UInt integer;
        if(shift>0){integer=magnitude_.right(shift);auto remainder=magnitude_;remainder.subtract(integer.shifted(shift));
            const int half=remainder.compare(Internal::UInt(1).shifted(shift-1));if(half>0||(half==0&&integer.odd()))integer.add(Internal::UInt(1));
        }else integer=magnitude_.shifted(-shift);
        if(integer.zero())return std::copysign(0.,negative_?-1.:1.);
        // The rounded integer has <=p bits, or is exactly the carry 2^p;
        // conversion is exact in the binary64 significand before power scaling.
        const double significand=integer.n.convert_to<double>();
        const auto value=std::scalbn(significand,quantum);
        Internal::require(std::isfinite(value),"[FiniteVolume] final binary64 dyadic overflow");
        Internal::require(value>0,"[FiniteVolume] binary64 dyadic conversion lost positive quantum");return negative_?-value:value;
    }
private:
    Internal::UInt magnitude_;
    int exponent_=0;
    bool negative_=false;
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE) && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST)
    // A bounded INTEGER alternative to constructing/multiplying one UInt per
    // native binary64 factor. The original scalar FP prefix, first-zero short
    // circuit and final Dyadic normalization/terminal RNE are unchanged.
    static Dyadic smallNativeProduct(std::initializer_list<long double> factors)
    {
        Dyadic result(1);std::array<std::uint64_t,4> words{{1,0,0,0}};std::size_t used=1;
        for(const auto value:factors){
            Internal::require(std::isfinite(value),"[FiniteVolume] nonfinite dyadic input");
            Internal::require(std::numeric_limits<long double>::radix==2,"[FiniteVolume] binary native arithmetic required");
            if(value==0)return {};const bool negative=std::signbit(value);
            std::uint64_t significand=0;int exponent=0;
            if(!Internal::decodeNormalBinary64(value,significand,exponent)){
                int e=0;auto f=std::frexp(std::abs(value),&e);
                constexpr int p=std::numeric_limits<long double>::digits;
                const auto integer=std::scalbn(f,p);
                Internal::require(integer>=1&&integer<std::scalbn(1.L,64),"[FiniteVolume] native significand outside uint64");
                significand=static_cast<std::uint64_t>(integer);exponent=e-p;
            }
            unsigned long trailing=0;_BitScanForward64(&trailing,significand);
            significand>>=trailing;
            result.negative_=result.negative_!=negative;
            result.exponent_+=exponent+static_cast<int>(trailing);
            std::uint64_t carry=0;
            for(std::size_t k=0;k<used;++k){
                std::uint64_t high=0;const auto low=_umul128(words[k],significand,&high);
                const auto next=low+carry;
                carry=high+static_cast<std::uint64_t>(next<low);words[k]=next;
            }
            if(carry){
                // Each normalized factor has <=53 bits. At most four factors
                // have <=212 bits: this fifth-word guard is unreachable for
                // every finite input admitted by the complete format gate.
                Internal::require(used<words.size(),"[FiniteVolume] proved native product limb bound exceeded");
                words[used++]=carry;
            }
        }
        // Import once into the EXISTING unbounded/MinBits256 UInt. There is no
        // product alignment, floating multiply, truncation or early rounding.
        Internal::UInt::Integer magnitude=0;
        for(std::size_t k=used;k>0;--k){magnitude<<=64;magnitude+=words[k-1];}
        result.magnitude_=Internal::UInt(std::move(magnitude));result.normalize();return result;
    }
#endif
    void normalize(){if(magnitude_.zero()){exponent_=0;negative_=false;return;}const int trailing=magnitude_.trailing();if(trailing==0)return;magnitude_.n>>=trailing;exponent_+=trailing;}
    static Dyadic allocation(Internal::UInt count,int exponent){Dyadic r;r.magnitude_=std::move(count);r.exponent_=exponent;r.normalize();return r;}
    friend class Binary64ExactProducts;
    friend Partition partition(long double,const std::vector<Dyadic>&);
    friend PartitionDouble partitionDouble(double,const std::vector<Dyadic>&);
};
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE) && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST)
// Current-scalar local arithmetic only. No geometry, field, stock or cache.
// Literal Dyadic remains the consumer fallback on unsupported native formats.
class NativeSmallExact
{
    struct Scalar {std::uint64_t significand=0;int exponent=0;bool negative=false;};
    static Scalar scalar(long double value)
    {
        // Keep the original scalar guards and first-zero short circuit.
        Internal::require(std::isfinite(value),"[FiniteVolume] nonfinite dyadic input");
        Internal::require(std::numeric_limits<long double>::radix==2,"[FiniteVolume] binary native arithmetic required");
        Scalar result;if(value==0)return result;result.negative=std::signbit(value);
        if(!Internal::decodeNormalBinary64(value,result.significand,result.exponent)){
            int e=0;auto f=std::frexp(std::abs(value),&e);
            constexpr int p=std::numeric_limits<long double>::digits;
            const auto integer=std::scalbn(f,p);
            Internal::require(integer>=1&&integer<std::scalbn(1.L,64),"[FiniteVolume] native significand outside uint64");
            result.significand=static_cast<std::uint64_t>(integer);result.exponent=e-p;
        }
        unsigned long trailing=0;_BitScanForward64(&trailing,result.significand);
        result.significand>>=trailing;result.exponent+=static_cast<int>(trailing);return result;
    }
    static int bitCount(const std::uint64_t* words,std::size_t used) noexcept
    {
        if(!used)return 0;
        unsigned long top=0;_BitScanReverse64(&top,words[used-1]);
        return static_cast<int>((used-1)*64)+static_cast<int>(top)+1;
    }
    static double roundWords(const std::uint64_t* words,std::size_t used,int exponent,bool negative)
    {
        // Same environment, overflow/quantum and floating suffix as Dyadic.
        Internal::binary64Environment();if(!used)return 0;
        constexpr int p=std::numeric_limits<double>::digits;
        constexpr int least=std::numeric_limits<double>::min_exponent-p;
        const int top=bitCount(words,used)-1+exponent;
        Internal::require(top<std::numeric_limits<double>::max_exponent,"[FiniteVolume] final binary64 dyadic overflow");
        const int quantum=std::max(top-(p-1),least);const int shift=quantum-exponent;
        const auto word=[&](std::size_t at) noexcept ->std::uint64_t {return at<used?words[at]:UINT64_C(0);};
        const auto integer=Internal::nativeTerminalRneWords(word,used,shift);
        if(integer==0)return std::copysign(0.,negative?-1.:1.);
        const double significand=static_cast<double>(integer);
        const auto value=std::scalbn(significand,quantum);
        Internal::require(std::isfinite(value),"[FiniteVolume] final binary64 dyadic overflow");
        Internal::require(value>0,"[FiniteVolume] binary64 dyadic conversion lost positive quantum");return negative?-value:value;
    }
public:
    static constexpr bool supported=sizeof(long double)==8&&sizeof(double)==8
        &&std::numeric_limits<long double>::radix==2&&std::numeric_limits<long double>::digits==53
        &&std::numeric_limits<long double>::min_exponent==-1021&&std::numeric_limits<long double>::max_exponent==1024
        &&std::numeric_limits<long double>::has_denorm==std::denorm_present;
    class Product
    {
        std::array<std::uint64_t,4> words_{};std::size_t used_=0;int exponent_=0;bool negative_=false;
        friend class NativeSmallExact;
    public:
        Product()=default;
        double roundedDouble() const {return roundWords(words_.data(),used_,exponent_,negative_);}
    };
    class Sum
    {
        // Unused words are never read, including on a value copy/move.
        std::array<std::uint64_t,132> words_;std::size_t used_=0;int exponent_=0;bool negative_=false;
        friend class NativeSmallExact;
    public:
        Sum() noexcept=default;
        Sum(const Sum& other) noexcept:used_(other.used_),exponent_(other.exponent_),negative_(other.negative_)
        {std::copy_n(other.words_.data(),used_,words_.data());}
        Sum& operator=(const Sum& other) noexcept
        {if(this!=&other){std::copy_n(other.words_.data(),other.used_,words_.data());used_=other.used_;
            exponent_=other.exponent_;negative_=other.negative_;}return *this;}
        Sum(Sum&& other) noexcept:Sum(static_cast<const Sum&>(other)){}
        Sum& operator=(Sum&& other) noexcept{return *this=static_cast<const Sum&>(other);}
        double roundedDouble() const {return roundWords(words_.data(),used_,exponent_,negative_);}
    };
    static Product product(std::initializer_list<long double> factors)
    {
        Internal::require(supported&&factors.size()>=2&&factors.size()<=4,
            "[FiniteVolume] native small exact product requires supported two-to-four-factor format");
        // smallNativeProduct begins with Dyadic(1); preserve its scalar FP prefix.
        const auto one=scalar(1.L);Product result;result.words_[0]=one.significand;
        result.used_=1;result.exponent_=one.exponent;result.negative_=one.negative;
        for(const auto value:factors){
            const auto part=scalar(value);if(part.significand==0)return {};
            result.negative_=result.negative_!=part.negative;result.exponent_+=part.exponent;
            std::uint64_t carry=0;
            for(std::size_t k=0;k<result.used_;++k){
                std::uint64_t high=0;const auto low=_umul128(result.words_[k],part.significand,&high);
                const auto next=low+carry;carry=high+static_cast<std::uint64_t>(next<low);result.words_[k]=next;
            }
            if(carry){Internal::require(result.used_<result.words_.size(),
                "[FiniteVolume] proved native product limb bound exceeded");result.words_[result.used_++]=carry;}
        }
        return result; // All nonzero factors are odd: no further normalization.
    }
    static Sum combine(const Product& a,int signA,const Product& b,int signB)
    {
        Internal::require(signA==1||signA==-1,"[FiniteVolume] invalid dyadic addition sign");
        Internal::require(signB==1||signB==-1,"[FiniteVolume] invalid dyadic addition sign");
        Sum result;
        if(!a.used_&&!b.used_)return result;
        const auto copy=[&](const Product& p,int sign){
            result.used_=p.used_;result.exponent_=p.exponent_;result.negative_=p.negative_!=(sign<0);
            std::copy_n(p.words_.data(),p.used_,result.words_.data());};
        if(!a.used_){copy(b,signB);return result;}if(!b.used_){copy(a,signA);return result;}
        result.exponent_=std::min(a.exponent_,b.exponent_);
        const int shiftA=a.exponent_-result.exponent_,shiftB=b.exponent_-result.exponent_;
        const int high=std::max(bitCount(a.words_.data(),a.used_)-1+a.exponent_,
            bitCount(b.words_.data(),b.used_)-1+b.exponent_);
        const auto count=static_cast<std::size_t>((high-result.exponent_+2+63)/64);
        Internal::require(count<=result.words_.size(),"[FiniteVolume] proved native small sum limb bound exceeded");
        const auto shifted=[](const Product& p,int shift,std::size_t at) noexcept {
            const auto base=static_cast<std::size_t>(shift/64);if(at<base)return UINT64_C(0);
            const auto i=at-base;const auto bits=static_cast<unsigned>(shift%64);std::uint64_t value=0;
            if(i<p.used_)value=p.words_[i]<<bits;
            if(bits&&i>0&&i-1<p.used_)value|=p.words_[i-1]>>(64-bits);return value;};
        const bool negativeA=a.negative_!=(signA<0),negativeB=b.negative_!=(signB<0);
        if(negativeA==negativeB){
            result.negative_=negativeA;std::uint64_t carry=0;
            for(std::size_t k=0;k<count;++k){const auto x=shifted(a,shiftA,k),y=shifted(b,shiftB,k);
                const auto added=x+y,next=added+carry;carry=static_cast<std::uint64_t>(added<x||next<added);result.words_[k]=next;}
            Internal::require(carry==0,"[FiniteVolume] proved native small sum carry exceeded");
        }else{
            int comparison=0;
            for(std::size_t k=count;k>0&&!comparison;--k){const auto x=shifted(a,shiftA,k-1),y=shifted(b,shiftB,k-1);
                if(x!=y)comparison=x<y?-1:1;}
            if(!comparison){result.exponent_=0;return result;}
            result.negative_=comparison>0?negativeA:negativeB;std::uint64_t borrow=0;
            for(std::size_t k=0;k<count;++k){const auto x=comparison>0?shifted(a,shiftA,k):shifted(b,shiftB,k);
                const auto y=comparison>0?shifted(b,shiftB,k):shifted(a,shiftA,k);
                const auto subtracted=x-y,next=subtracted-borrow;borrow=static_cast<std::uint64_t>(x<y||subtracted<borrow);result.words_[k]=next;}
            Internal::require(borrow==0,"[FiniteVolume] proved native small sum borrow exceeded");
        }
        result.used_=count;while(result.used_&&result.words_[result.used_-1]==0)--result.used_;
        if(!result.used_){result.exponent_=0;result.negative_=false;}return result;
    }
};
#endif
// Compiled two-binary64 factors and fresh binary64 values. This neutral
// exact-arithmetic port owns no physical inventory/field/geometry receipt.
// It deliberately reuses Dyadic's original terminal rounding implementation.
class Binary64ExactProducts
{
public:
    class Factor
    {
    public:
        Factor()=default;
        // Exact sign reversal; magnitude/exponent and zero factors stay exact.
        Factor negated() const noexcept
        {auto result=*this;result.packedHigh_^=UINT64_C(1)<<55;return result;}
    private:
        std::uint64_t low_=0,packedHigh_=0;
        friend class Binary64ExactProducts;
    };
    class Value
    {
    public: Value()=default;
    private:
        std::uint64_t significand_=0;int exponent_=0;bool negative_=false;
        friend class Binary64ExactProducts;
    };
    struct Result {Dyadic sum,absolute;};
    static Value decode(double value)
    {
        Internal::require(sizeof(double)==sizeof(std::uint64_t)&&std::numeric_limits<double>::is_iec559
            &&std::numeric_limits<double>::digits==53,"[FiniteVolume] exact product requires binary64");
        std::uint64_t bits=0;
        if constexpr(sizeof(double)==sizeof(bits))std::memcpy(&bits,&value,sizeof(bits));
        const auto field=static_cast<unsigned>((bits>>52)&0x7ffu);
        Internal::require(field!=0x7ffu,"[FiniteVolume] nonfinite compiled product input");
        Value result;result.significand_=bits&UINT64_C(0x000fffffffffffff);
        if(field==0){if(result.significand_==0)return result;result.exponent_=-1074;}
        else {result.significand_|=UINT64_C(0x0010000000000000);result.exponent_=static_cast<int>(field)-1023-52;}
        result.negative_=(bits>>63)!=0;
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
        unsigned long trailing=0;_BitScanForward64(&trailing,result.significand_);
        result.significand_>>=trailing;result.exponent_+=static_cast<int>(trailing);
#else
        while((result.significand_&1u)==0){result.significand_>>=1;++result.exponent_;}
#endif
        return result;
    }
    static Factor compile(double first,double second)
    {
        const auto a=decode(first),b=decode(second);Factor result;
        if(a.significand_==0||b.significand_==0)return result;
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
        // The decoded factors each contain at most 53 bits. Their complete
        // unsigned product is at most 106 bits; no floating multiply occurs.
        std::uint64_t high=0;
        const auto low=_umul128(a.significand_,b.significand_,&high);
#else
        using Integer=Internal::UInt::Integer;
        const Integer product=Integer(a.significand_)*b.significand_;
        const Integer mask=(Integer(1)<<64)-1;
        const Integer highInteger=product>>64,lowInteger=product&mask;
        const auto high=highInteger.convert_to<std::uint64_t>();
#endif
        const int exponent=a.exponent_+b.exponent_;
        Internal::require(high<(UINT64_C(1)<<42)&&exponent>=-2148&&exponent<=2046,
            "[FiniteVolume] compiled factor outside proved binary64 range");
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
        result.low_=low;
#else
        result.low_=lowInteger.convert_to<std::uint64_t>();
#endif
        result.packedHigh_=high|(static_cast<std::uint64_t>(exponent+2148)<<42)
            |(static_cast<std::uint64_t>(a.negative_!=b.negative_)<<55);
        return result;
    }
    static Result sum(const Factor* factors,const std::uint32_t* rows,std::size_t count,
        const Value* values,std::size_t valueCount)
    {
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
        return sumLimbs(factors,rows,count,values,valueCount);
#else
        Internal::require(count==0||(factors&&rows),"[FiniteVolume] missing compiled product terms");
        Internal::require(valueCount==0||values,"[FiniteVolume] missing compiled product values");
        constexpr auto highMask=(UINT64_C(1)<<42)-1;
        const auto magnitudeHigh=[&](const Factor& f){return f.packedHigh_&highMask;};
        const auto exponent=[](const Factor& f){return static_cast<int>((f.packedHigh_>>42)&0x1fffu)-2148;};
        bool present=false;int quantum=0;
        for(std::size_t t=0;t<count;++t){
            Internal::require(rows[t]<valueCount,"[FiniteVolume] compiled product row out of range");
            const auto& factor=factors[t];const auto& value=values[rows[t]];
            if(value.significand_==0||(factor.low_==0&&magnitudeHigh(factor)==0))continue;
            const int termExponent=exponent(factor)+value.exponent_;
            Internal::require(termExponent>=-3222&&termExponent<=3069,"[FiniteVolume] compiled product exponent out of range");
            if(!present||termExponent<quantum)quantum=termExponent;present=true;
        }
        Result result;if(!present)return result;
        using Integer=Internal::UInt::Integer;Integer signedSum=0,absoluteSum=0;
        for(std::size_t t=0;t<count;++t){
            const auto& factor=factors[t];const auto& value=values[rows[t]];
            if(value.significand_==0||(factor.low_==0&&magnitudeHigh(factor)==0))continue;
            Integer term=magnitudeHigh(factor);term<<=64;term+=factor.low_;term*=value.significand_;
            const int shift=exponent(factor)+value.exponent_-quantum;
            Internal::require(shift>=0&&shift<=6291,"[FiniteVolume] compiled product alignment out of range");
            if(shift)term<<=shift;
            absoluteSum+=term;
            if((((factor.packedHigh_>>55)&1u)!=0)!=value.negative_)signedSum-=term;else signedSum+=term;
        }
        result.sum.negative_=signedSum<0;
        if(result.sum.negative_)signedSum=-signedSum;
        result.sum.magnitude_=Internal::UInt(std::move(signedSum));result.sum.exponent_=quantum;result.sum.normalize();
        result.absolute.magnitude_=Internal::UInt(std::move(absoluteSum));result.absolute.exponent_=quantum;result.absolute.normalize();
        return result;
#endif
    }
private:
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE)
    static Result sumLimbs(const Factor* factors,const std::uint32_t* rows,std::size_t count,
        const Value* values,std::size_t valueCount)
    {
        static_assert(std::numeric_limits<std::size_t>::digits<=64,"fixed accumulator size_t proof");
        Internal::require(count==0||(factors&&rows),"[FiniteVolume] missing compiled product terms");
        Internal::require(valueCount==0||values,"[FiniteVolume] missing compiled product values");
        constexpr auto highMask=(UINT64_C(1)<<42)-1;
        const auto magnitudeHigh=[&](const Factor& f){return f.packedHigh_&highMask;};
        const auto exponent=[](const Factor& f){return static_cast<int>((f.packedHigh_>>42)&0x1fffu)-2148;};
        bool present=false;int quantum=0,maximum=0;
        // Keep the original complete row/exponent validation pass and its
        // first-error order, before any product accumulation.
        for(std::size_t t=0;t<count;++t){
            Internal::require(rows[t]<valueCount,"[FiniteVolume] compiled product row out of range");
            const auto& factor=factors[t];const auto& value=values[rows[t]];
            if(value.significand_==0||(factor.low_==0&&magnitudeHigh(factor)==0))continue;
            const int termExponent=exponent(factor)+value.exponent_;
            Internal::require(termExponent>=-3222&&termExponent<=3069,"[FiniteVolume] compiled product exponent out of range");
            if(!present||termExponent<quantum)quantum=termExponent;
            if(!present||termExponent>maximum)maximum=termExponent;present=true;
        }
        Result result;if(!present)return result;
        // 106-bit compiled factor * 53-bit fresh value is <2^159.
        // Alignment is <=6291 bits. Even count=SIZE_MAX (<2^64) gives
        // either nonnegative subtotal <2^(159+6291+64)=2^6514.
        // 102*64=6528 bits therefore includes every possible final carry.
        constexpr std::size_t capacity=102;
        unsigned carryBits=0;
        for(auto remaining=count-1;remaining;remaining>>=1)++carryBits;
        const auto limbCount=static_cast<std::size_t>((159+maximum-quantum+carryBits+63)/64);
        Internal::require(limbCount<=capacity,"[FiniteVolume] proved product accumulator span exceeded");
        std::array<std::uint64_t,capacity> positive,negative;
        std::fill_n(positive.begin(),limbCount,UINT64_C(0));
        std::fill_n(negative.begin(),limbCount,UINT64_C(0));
        for(std::size_t t=0;t<count;++t){
            const auto& factor=factors[t];const auto& value=values[rows[t]];
            if(value.significand_==0||(factor.low_==0&&magnitudeHigh(factor)==0))continue;
            std::uint64_t lowerHigh=0,upperHigh=0;
            const auto lower=_umul128(factor.low_,value.significand_,&lowerHigh);
            const auto upperLow=_umul128(magnitudeHigh(factor),value.significand_,&upperHigh);
            const std::uint64_t middle=lowerHigh+upperLow;
            // upperHigh <2^31; this carry cannot overflow its uint64 word.
            const std::uint64_t high=upperHigh+static_cast<std::uint64_t>(middle<lowerHigh);
            const int shift=exponent(factor)+value.exponent_-quantum;
            Internal::require(shift>=0&&shift<=6291,"[FiniteVolume] compiled product alignment out of range");
            const auto word=static_cast<std::size_t>(shift/64);
            const unsigned bits=static_cast<unsigned>(shift%64);
            std::array<std::uint64_t,4> part{{lower,middle,high,0}};
            if(bits){
                part[3]=high>>(64-bits);
                part[2]=(high<<bits)|(middle>>(64-bits));
                part[1]=(middle<<bits)|(lower>>(64-bits));
                part[0]=lower<<bits;
            }
            std::size_t width=4;while(width&&part[width-1]==0)--width;
            auto& accumulator=(((((factor.packedHigh_>>55)&1u)!=0)!=value.negative_)?negative:positive);
            Internal::require(word+width<=limbCount,"[FiniteVolume] proved product word span exceeded");
            std::uint64_t carry=0;std::size_t at=word;
            for(std::size_t k=0;k<width;++k,++at){
                const auto before=accumulator[at],added=before+part[k];
                const auto next=added+carry;
                carry=static_cast<std::uint64_t>(added<before||next<added);
                accumulator[at]=next;
            }
            while(carry){
                Internal::require(at<limbCount,"[FiniteVolume] proved product carry span exceeded");
                const auto before=accumulator[at];accumulator[at]=before+1;
                carry=static_cast<std::uint64_t>(accumulator[at]==0);++at;
            }
        }
        using Integer=Internal::UInt::Integer;
        const auto integer=[&](const std::array<std::uint64_t,capacity>& words){
            Integer value=0;
            for(std::size_t k=limbCount;k>0;--k){value<<=64;value+=words[k-1];}
            return value;
        };
        // Only terminal per-column conversion uses multiprecision; no
        // per-term cpp_int construction, multiplication or alignment shift.
        Integer signedSum=integer(positive),absoluteSum=signedSum;
        const Integer negativeSum=integer(negative);
        signedSum-=negativeSum;absoluteSum+=negativeSum;
        result.sum.negative_=signedSum<0;
        if(result.sum.negative_)signedSum=-signedSum;
        result.sum.magnitude_=Internal::UInt(std::move(signedSum));result.sum.exponent_=quantum;result.sum.normalize();
        result.absolute.magnitude_=Internal::UInt(std::move(absoluteSum));result.absolute.exponent_=quantum;result.absolute.normalize();
        return result;
    }
#endif
};
static_assert(sizeof(Binary64ExactProducts::Factor)==16,"compiled binary64 factor storage contract");

struct Rounded{long double value=0;Dyadic error;}; // error = represented - exact.
inline Rounded round(const Dyadic& exact)
{Rounded r;r.value=exact.rounded();r.error=Dyadic(r.value);r.error.add(exact,-1);return r;}
inline Rounded phaseUpdate(long double oldStock,const Dyadic& exactPhaseNet,const Dyadic& exactCompression)
{
    Internal::require(std::isfinite(oldStock)&&oldStock>=0,"[FiniteVolume] invalid source phase stock");
    if(exactPhaseNet.compare(Dyadic{})==0&&exactCompression.compare(Dyadic{})==0)return {oldStock,{}};
    Dyadic exact(oldStock);exact.add(exactPhaseNet,-1);exact.add(exactCompression);return round(exact);
}
inline Rounded phaseUpdate(long double oldStock,long double phaseNet,long double compression)
{return phaseUpdate(oldStock,Dyadic(phaseNet),Dyadic(compression));}
inline Dyadic rectangleCapacityExact(const PlanarCellCut::Rectangle& box,bool rz)
{
    Internal::require(std::isfinite(box.lower.x)&&std::isfinite(box.lower.y)&&std::isfinite(box.upper.x)&&std::isfinite(box.upper.y)&&
        box.upper.x>box.lower.x&&box.upper.y>box.lower.y&&(!rz||box.lower.x>=0),"[FiniteVolume] invalid exact capacity rectangle");
    const long double x0=box.lower.x,x1=box.upper.x,y0=box.lower.y,y1=box.upper.y;Dyadic result;
    if(rz){const long double halfTwoPi=static_cast<long double>(CartesianMeasure::twoPi)*.5L;
        result.add(Dyadic::product({halfTwoPi,x1,x1,y1}));result.add(Dyadic::product({halfTwoPi,x1,x1,y0}),-1);
        result.add(Dyadic::product({halfTwoPi,x0,x0,y1}),-1);result.add(Dyadic::product({halfTwoPi,x0,x0,y0}));
    }else{result.add(Dyadic::product({x1,y1}));result.add(Dyadic::product({x1,y0}),-1);
        result.add(Dyadic::product({x0,y1}),-1);result.add(Dyadic::product({x0,y0}));}
    Internal::require(result.compare(Dyadic{})>0,"[FiniteVolume] nonpositive exact rectangle capacity");return result;
}
struct Partition
{
    std::vector<long double> values;
    std::vector<Dyadic> exactAllocations;
    Dyadic sumError;
};
inline Partition partition(long double sourceStock,const std::vector<Dyadic>& geometry)
{
    Internal::environment();Internal::require(std::isfinite(sourceStock)&&sourceStock>=0,"[FiniteVolume] invalid partition source stock");
    Dyadic whole;for(const auto& g:geometry){Internal::require(g.compare(Dyadic{})>=0,"[FiniteVolume] negative partition geometry");whole.add(g);}
    Partition r;r.values.reserve(geometry.size());r.exactAllocations.reserve(geometry.size());
    if(sourceStock==0){r.values.assign(geometry.size(),0);r.exactAllocations.resize(geometry.size());return r;}
    Internal::require(whole.compare(Dyadic{})>0,"[FiniteVolume] positive stock has no partition geometry");
    int sourceExponent=0;std::frexp(sourceStock,&sourceExponent);
    constexpr int p=std::numeric_limits<long double>::digits;const int quantum=std::max(sourceExponent-p,std::numeric_limits<long double>::min_exponent-p);
    const Dyadic source(sourceStock);Internal::require(source.exponent_>=quantum,"[FiniteVolume] source stock below native quantum");
    const auto count=source.magnitude_.shifted(source.exponent_-quantum);Internal::UInt prior;Dyadic prefix,sum;
    for(const auto& g:geometry){prefix.add(g);Internal::UInt cumulative;
        const int endpoint=prefix.compare(whole);Internal::require(endpoint<=0,"[FiniteVolume] partition prefix exceeds whole");
        if(endpoint==0)cumulative=count;
        else if(prefix.compare(Dyadic{})>0){auto numerator=count.multiplied(prefix.magnitude_);auto denominator=whole.magnitude_;
            const int shift=prefix.exponent_-whole.exponent_;if(shift>=0)numerator=numerator.shifted(shift);else denominator=denominator.shifted(-shift);
            cumulative=Internal::nearestInteger(numerator,denominator);}
        Internal::require(cumulative.compare(prior)>=0&&cumulative.compare(count)<=0,"[FiniteVolume] partition CDF is not monotone/bounded");
        auto pieceCount=cumulative;pieceCount.subtract(prior);auto allocation=Dyadic::allocation(std::move(pieceCount),quantum);
        const auto value=allocation.rounded();Internal::require(Dyadic(value).compare(allocation)==0,"[FiniteVolume] partition allocation is not native representable");
        r.values.push_back(value);r.exactAllocations.push_back(allocation);sum.add(allocation);prior=std::move(cumulative);
    }
    Internal::require(prefix.compare(whole)==0&&prior.compare(count)==0,"[FiniteVolume] incomplete partition endpoints");
    sum.add(source,-1);Internal::require(sum.compare(Dyadic{})==0,"[FiniteVolume] partition lost source phase");r.sumError=std::move(sum);return r;
}

// Error is represented binary64 value minus complete exact input. No physical
// charge is owned here; signed physical consumers keep their own provenance.
struct RoundedDouble{double value=0;Dyadic error;};
inline RoundedDouble roundDouble(const Dyadic& exact)
{RoundedDouble r;r.value=exact.roundedDouble();r.error=Dyadic(r.value);r.error.add(exact,-1);return r;}
struct PartitionDouble
{
    std::vector<double> values;
    std::vector<Dyadic> exactAllocations;
    Dyadic sumError;
};
// Nonnegative magnitude only. Source binary64 normal ULP or denorm_min defines
// every integer allocation, regardless of long-double precision/range.
inline PartitionDouble partitionDouble(double sourceStock,const std::vector<Dyadic>& geometry)
{
    Internal::binary64Environment();Internal::require(std::isfinite(sourceStock)&&sourceStock>=0,"[FiniteVolume] invalid partition source stock");
    Dyadic whole;for(const auto& g:geometry){Internal::require(g.compare(Dyadic{})>=0,"[FiniteVolume] negative partition geometry");whole.add(g);}
    PartitionDouble r;r.values.reserve(geometry.size());r.exactAllocations.reserve(geometry.size());
    if(sourceStock==0){r.values.assign(geometry.size(),0);r.exactAllocations.resize(geometry.size());return r;}
    Internal::require(whole.compare(Dyadic{})>0,"[FiniteVolume] positive stock has no partition geometry");
    int sourceExponent=0;std::frexp(sourceStock,&sourceExponent);
    constexpr int p=std::numeric_limits<double>::digits;const int quantum=std::max(sourceExponent-p,std::numeric_limits<double>::min_exponent-p);
    const Dyadic source(sourceStock);Internal::require(source.exponent_>=quantum,"[FiniteVolume] source stock below binary64 quantum");
    const auto count=source.magnitude_.shifted(source.exponent_-quantum);Internal::UInt prior;Dyadic prefix,sum;
    for(const auto& g:geometry){prefix.add(g);Internal::UInt cumulative;
        const int endpoint=prefix.compare(whole);Internal::require(endpoint<=0,"[FiniteVolume] partition prefix exceeds whole");
        if(endpoint==0)cumulative=count;
        else if(prefix.compare(Dyadic{})>0){auto numerator=count.multiplied(prefix.magnitude_);auto denominator=whole.magnitude_;
            const int shift=prefix.exponent_-whole.exponent_;if(shift>=0)numerator=numerator.shifted(shift);else denominator=denominator.shifted(-shift);
            cumulative=Internal::nearestInteger(numerator,denominator);}
        Internal::require(cumulative.compare(prior)>=0&&cumulative.compare(count)<=0,"[FiniteVolume] partition CDF is not monotone/bounded");
        auto pieceCount=cumulative;pieceCount.subtract(prior);auto allocation=Dyadic::allocation(std::move(pieceCount),quantum);
        const auto value=allocation.roundedDouble();Internal::require(Dyadic(value).compare(allocation)==0,"[FiniteVolume] partition allocation is not binary64 representable");
        r.values.push_back(value);r.exactAllocations.push_back(allocation);sum.add(allocation);prior=std::move(cumulative);
    }
    Internal::require(prefix.compare(whole)==0&&prior.compare(count)==0,"[FiniteVolume] incomplete partition endpoints");
    sum.add(source,-1);Internal::require(sum.compare(Dyadic{})==0,"[FiniteVolume] partition lost binary64 source stock");r.sumError=std::move(sum);return r;
}
}