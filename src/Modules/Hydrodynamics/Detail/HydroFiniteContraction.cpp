#include "Modules/Hydrodynamics/Detail/HydroFiniteContraction.h"
#include <boost/multiprecision/cpp_int.hpp>
#include <algorithm>
#include <string>
#include <cfenv>
#include <cstdint>
#include <cstddef>
#include <cstring>
#include <cmath>
#include <limits>
#include <utility>
#include <tuple>

namespace LAMP::Hydrodynamics::Detail
{
namespace
{
using I=boost::multiprecision::cpp_int;
thread_local bool rangeEnabled=false;
[[noreturn]] void fail(const char* message)
{throw std::invalid_argument(std::string("[HydroFiniteContraction] ")+message);}
I absolute(I v){return v<0?-v:v;}
I gcd(I a,I b)
{
    a=absolute(a);b=absolute(b);if(a==0)return b;if(b==0)return a;
    if(a==b)return a;
    // gcd(2^i*aOdd,2^j*bOdd)=2^min(i,j)*gcd(aOdd,bOdd).
    // In particular, a denominator odd(dt)*2^k no longer sends all k
    // binary scale bits through general integer remainder operations.
    const auto aTwos=boost::multiprecision::lsb(a),bTwos=boost::multiprecision::lsb(b);
    const auto commonTwos=std::min(aTwos,bTwos);a>>=aTwos;b>>=bTwos;
    if(a==1||b==1)return I(1)<<commonTwos;
    while(b!=0){I r=a%b;a=std::move(b);b=std::move(r);}
    return a<<commonTwos;
}
struct Dyadic{I significand;int exponent=0;};
Dyadic exact(long double value)
{
    if(!std::isfinite(value))fail("nonfinite input factor");
    if(value==0)return {};
    int exponent=0;auto fraction=std::frexp(std::abs(value),&exponent);I bits=0;
    constexpr int precision=std::numeric_limits<long double>::digits;
    if constexpr(precision<=64){
        // frexp supplies fraction in [1/2,1), with at most p bits. Scaling
        // by 2^p is exact and gives an integer in [2^(p-1),2^p-1]. Even
        // p=64 is strictly below 2^64, so the uint64_t cast cannot overflow.
        const auto integer=std::scalbn(fraction,precision);
        if(!(integer>=1&&integer<std::scalbn(1.L,64)))fail("unsupported floating significand");
        const auto word=static_cast<std::uint64_t>(integer);
        if(static_cast<long double>(word)!=integer)fail("unsupported floating significand");
        bits=word;
    }else{
        // Wider binary formats remain exact: remove at most32 leading bits
        // at once. Each chunk is an integer below2^32; its subtraction is
        // exact, and the final residual must be exactly zero as before.
        for(int remaining=precision;remaining>0;){const int count=std::min(remaining,32);
            const auto shifted=std::scalbn(fraction,count);const auto word=static_cast<std::uint32_t>(shifted);
            bits<<=count;bits+=word;fraction=shifted-static_cast<long double>(word);remaining-=count;}
        if(fraction!=0)fail("unsupported floating significand");
    }
    // Remove only exact binary factors. This changes the representation,
    // not its value, and avoids building p-bit powers of two for 1 and .5.
    const auto trailing=boost::multiprecision::lsb(bits);bits>>=trailing;
    if(value<0)bits=-bits;
    return {std::move(bits),exponent-std::numeric_limits<long double>::digits+static_cast<int>(trailing)};
}
void environment()
{
    if(std::numeric_limits<long double>::radix!=2 || std::fegetround()!=FE_TONEAREST
        || std::numeric_limits<long double>::has_denorm!=std::denorm_present
        || std::scalbn(1.L,std::numeric_limits<long double>::min_exponent-std::numeric_limits<long double>::digits)
            !=std::numeric_limits<long double>::denorm_min())fail("requires binary gradual-underflow round-to-nearest environment");
}
}
struct HydroFiniteContraction::Data
{
    I numerator=0,denominator=1;
    // The fast tag is native binary64 only. General remains canonical and
    // never switches back except through the original reset/fresh Data.
    static constexpr bool binary64=std::numeric_limits<long double>::radix==2
        &&std::numeric_limits<long double>::digits==53
        &&std::numeric_limits<long double>::max_exponent==1024
        &&std::numeric_limits<long double>::min_exponent==-1021;
    bool fast=binary64;
    int binaryExponent=0;
    std::uint64_t oddDenominator=1;
    void add(I n,I d)
    {
        if(d==0)fail("zero contraction denominator");if(d<0){n=-n;d=-d;}if(n==0)return;
        // These are exact rational identities. Retain the original reduced
        // normal form after every nonzero sum, including equal denominators.
        if(numerator==0){numerator=std::move(n);denominator=std::move(d);}
        else if(denominator==d)numerator+=n;
        else if(denominator==1){numerator=numerator*d+n;denominator=std::move(d);}
        else if(d==1)numerator+=n*denominator;
        else{
            const I common=gcd(denominator,d),a=d/common,b=denominator/common;
            numerator=numerator*a+n*b;denominator*=a;
        }
        if(numerator==0){denominator=1;return;}
        if(denominator==1)return;
        const I reduce=gcd(numerator,denominator);numerator/=reduce;denominator/=reduce;
    }
    void expanded(I& n,I& d) const
    {
        n=numerator;d=fast?I(oddDenominator):denominator;
        if(fast){
            if(binaryExponent>=0)n<<=static_cast<unsigned long long>(binaryExponent);
            else d<<=static_cast<unsigned long long>(-static_cast<std::int64_t>(binaryExponent));
        }
    }
    void makeGeneral()
    {
        if(!fast)return;
        I n,d;expanded(n,d);
        numerator=0;denominator=1;fast=false;binaryExponent=0;oddDenominator=1;
        add(std::move(n),std::move(d));
    }
    void normalizeFast()
    {
        if(numerator==0){binaryExponent=0;oddDenominator=1;return;}
        std::size_t trailing=0;
        if(numerator<0){const I magnitude=-numerator;trailing=boost::multiprecision::lsb(magnitude);}
        else trailing=boost::multiprecision::lsb(numerator);
        const auto room=static_cast<std::int64_t>(std::numeric_limits<int>::max())-binaryExponent;
        if(trailing>static_cast<std::uint64_t>(room)){makeGeneral();return;}
        const auto next=static_cast<std::int64_t>(binaryExponent)+static_cast<std::int64_t>(trailing);
        if(trailing)numerator>>=trailing;
        binaryExponent=static_cast<int>(next);
    }
    void addFast(I n,int exponent,std::uint64_t odd)
    {
        // All native factors/divisor were already validated in original order.
        // n is a value snapshot, including when the source is *this.
        if(n==0)return;
        if(!fast){
            I d=odd;
            if(exponent>=0)n<<=static_cast<unsigned long long>(exponent);
            else d<<=static_cast<unsigned long long>(-static_cast<std::int64_t>(exponent));
            add(std::move(n),std::move(d));return;
        }
        if(numerator==0){
            numerator=std::move(n);binaryExponent=exponent;oddDenominator=odd;
            normalizeFast();return;
        }
        // Canonicalize BEFORE two distinct nonunit odd denominators multiply.
        // General then stays canonical; no lazy LCM history can accumulate.
        if(oddDenominator!=odd&&oddDenominator!=1&&odd!=1){
            makeGeneral();addFast(std::move(n),exponent,odd);return;
        }
        const auto common=std::min(binaryExponent,exponent);
        const auto targetShift=static_cast<unsigned long long>(
            static_cast<std::int64_t>(binaryExponent)-common);
        const auto sourceShift=static_cast<unsigned long long>(
            static_cast<std::int64_t>(exponent)-common);
        if(oddDenominator==1&&odd!=1)numerator*=odd;
        else if(odd==1&&oddDenominator!=1)n*=oddDenominator;
        numerator<<=targetShift;n<<=sourceShift;numerator+=n;
        binaryExponent=common;
        if(oddDenominator==1)oddDenominator=odd;
        normalizeFast();
    }
};
HydroFiniteContraction::HydroFiniteContraction():data_(std::make_unique<Data>()){environment();}
HydroFiniteContraction::~HydroFiniteContraction()=default;
void HydroFiniteContraction::reset()
{environment();*data_=Data{};}
void HydroFiniteContraction::resetQualifiedStorage()
{
    // Arithmetic assignment retains the cpp_int limb capacities; no fresh Data move.
    data_->numerator=0;data_->denominator=1;data_->fast=Data::binary64;
    data_->binaryExponent=0;data_->oddDenominator=1;
}
void HydroFiniteContraction::addBinary64Product(long double first,long double second,int binaryScale)
{
    if(!Data::binary64||sizeof(long double)!=8||(binaryScale!=0&&binaryScale!=-1))fail("invalid qualified binary64 product");
    const auto decode=[](long double value){
        std::uint64_t bits=0;std::memcpy(&bits,&value,sizeof bits);
        const auto field=static_cast<unsigned>((bits>>52)&0x7ffu);
        if(field==0x7ffu)fail("nonfinite input factor");
        std::uint64_t significand=bits&UINT64_C(0x000fffffffffffff);
        int exponent=-1074;if(field){significand|=UINT64_C(0x0010000000000000);exponent=static_cast<int>(field)-1023-52;}
        if(significand){while((significand&1u)==0){significand>>=1;++exponent;}}
        return std::tuple<std::uint64_t,int,bool>{significand,exponent,(bits>>63)!=0};
    };
    const auto [a,ae,an]=decode(first);const auto [b,be,bn]=decode(second);
    // At most106 exact significand bits; the exponent spans [-2149,2046].
    I numerator=I(a)*b;if(an!=bn)numerator=-numerator;
    data_->addFast(std::move(numerator),ae+be+binaryScale,1);
}
HydroFiniteContraction::HydroFiniteContraction(const HydroFiniteContraction& v):data_(std::make_unique<Data>(*v.data_)){}
HydroFiniteContraction& HydroFiniteContraction::operator=(const HydroFiniteContraction& v)
{if(this!=&v)*data_=*v.data_;return *this;}
HydroFiniteContraction::HydroFiniteContraction(HydroFiniteContraction&&) noexcept=default;
HydroFiniteContraction& HydroFiniteContraction::operator=(HydroFiniteContraction&&) noexcept=default;
void HydroFiniteContraction::add(long double v){addProduct({v});}
void HydroFiniteContraction::addProduct(std::initializer_list<long double> factors){addQuotientProduct(factors,1.L);}
void HydroFiniteContraction::addQuotientProduct(std::initializer_list<long double> factors,long double divisor)
{
    if(factors.size()==0||factors.size()>8)fail("invalid bounded contraction arity");
    auto den=exact(divisor);if(den.significand==0)fail("zero contraction denominator");
    I n=1;int exponent=-den.exponent;
    for(auto v:factors){auto x=exact(v);n*=x.significand;exponent+=x.exponent;}
    I d=std::move(den.significand);
    if constexpr(Data::binary64){
        if(d<0){n=-n;d=-d;}
        if(d>0&&d<=std::numeric_limits<std::uint64_t>::max()&&boost::multiprecision::bit_test(d,0)){
            data_->addFast(std::move(n),exponent,d.convert_to<std::uint64_t>());return;
        }
    }
    if(exponent>=0)n<<=exponent;else d<<=-exponent;
    data_->makeGeneral();data_->add(std::move(n),std::move(d));
}
void HydroFiniteContraction::add(const HydroFiniteContraction& v,int sign)
{if(sign!=1&&sign!=-1)fail("invalid contraction sum sign");
    if(v.data_->fast)data_->addFast(sign*v.data_->numerator,v.data_->binaryExponent,v.data_->oddDenominator);
    else{data_->makeGeneral();data_->add(sign*v.data_->numerator,v.data_->denominator);}}
void HydroFiniteContraction::addAbsolute(const HydroFiniteContraction& v)
{if(v.data_->fast)data_->addFast(absolute(v.data_->numerator),v.data_->binaryExponent,v.data_->oddDenominator);
    else{data_->makeGeneral();data_->add(absolute(v.data_->numerator),v.data_->denominator);}}
int HydroFiniteContraction::compare(const HydroFiniteContraction& other) const
{
    if(data_->fast&&other.data_->fast&&data_->oddDenominator==other.data_->oddDenominator){
        I left=data_->numerator,right=other.data_->numerator;
        const auto exponent=std::min(data_->binaryExponent,other.data_->binaryExponent);
        left<<=static_cast<unsigned long long>(static_cast<std::int64_t>(data_->binaryExponent)-exponent);
        right<<=static_cast<unsigned long long>(static_cast<std::int64_t>(other.data_->binaryExponent)-exponent);
        const I difference=left-right;return difference<0?-1:difference>0?1:0;
    }
    I n,d,otherN,otherD;data_->expanded(n,d);other.data_->expanded(otherN,otherD);
    const I difference=n*otherD-otherN*d;
    return difference<0?-1:difference>0?1:0;
}
long double HydroFiniteContraction::scaledQuotient(long double scale,const HydroFiniteContraction& denominator) const
{
    if(denominator.data_->numerator==0)fail("zero exact contraction denominator");
    const auto factor=exact(scale);
    I sourceN,sourceD,denominatorN,denominatorD;
    data_->expanded(sourceN,sourceD);denominator.data_->expanded(denominatorN,denominatorD);
    I n=sourceN*factor.significand*denominatorD;
    I d=sourceD*denominatorN;
    if(factor.exponent>=0)n<<=factor.exponent;else d<<=-factor.exponent;
    HydroFiniteContraction result;result.data_->makeGeneral();result.data_->add(std::move(n),std::move(d));return result.rounded();
}
long double HydroFiniteContraction::rounded() const
{return roundedImpl(true);}
long double HydroFiniteContraction::roundedQualified() const
{return roundedImpl(false);}
long double HydroFiniteContraction::roundedImpl(bool checkEnvironment) const
{
    if(checkEnvironment)environment();if(data_->numerator==0)return 0.L;
    const bool negative=data_->numerator<0;I n=absolute(data_->numerator),d=data_->denominator;
    if(data_->fast){
        d=data_->oddDenominator;
        if(data_->binaryExponent>=0)n<<=static_cast<unsigned long long>(data_->binaryExponent);
        else d<<=static_cast<unsigned long long>(-static_cast<std::int64_t>(data_->binaryExponent));
    }
    int exponent=static_cast<int>(boost::multiprecision::msb(n))-static_cast<int>(boost::multiprecision::msb(d));
    if(exponent>=0){if(n<(d<<exponent))--exponent;}else if((n<<-exponent)<d)--exponent;
    if(exponent>=std::numeric_limits<long double>::max_exponent)fail("final contraction overflow");
    const int quantum=std::max(exponent-(std::numeric_limits<long double>::digits-1),
        std::numeric_limits<long double>::min_exponent-std::numeric_limits<long double>::digits);
    if(quantum>=0)d<<=quantum;else n<<=-quantum;
    I rounded=n/d;const I remainder=n%d,twice=remainder<<1;
    if(twice>d || (twice==d && static_cast<bool>(rounded&1)))++rounded;
    if(rounded==0)return std::copysign(0.L,negative?-1.L:1.L);
    // At most p+1 bits; every construction operation below is exact.
    long double significand=0;bool nativeSignificand=false;
#if defined(_MSC_VER) && defined(_M_X64) && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST)
    if constexpr(Data::binary64&&sizeof(long double)==8){
        // Integer RNE is already complete. Every integer through 2^53 is
        // exactly representable, including the one-bit carry at 2^53. Use
        // signed conversion to avoid an unsigned-to-double lowering path.
        // Keep the original scaler and every nonnormal/overflow path below.
        if(rounded<=(std::uint64_t(1)<<53)){
            const int nativeExponent=static_cast<int>(boost::multiprecision::msb(rounded))+quantum;
            if(nativeExponent>=std::numeric_limits<long double>::min_exponent-1
                &&nativeExponent<std::numeric_limits<long double>::max_exponent){
                const auto word=rounded.convert_to<std::uint64_t>();
                significand=static_cast<long double>(static_cast<std::int64_t>(word));
                nativeSignificand=true;
            }
        }
    }
#endif
    if(!nativeSignificand)
    for(int bit=static_cast<int>(boost::multiprecision::msb(rounded));bit>=0;--bit)
        significand=significand*2+(boost::multiprecision::bit_test(rounded,static_cast<unsigned>(bit))?1:0);
    const auto result=std::scalbn(significand,quantum);
    if(!std::isfinite(result))fail("final contraction overflow");
    return negative?-result:result;
}
HydroFiniteRangeScope::HydroFiniteRangeScope(bool enabled) noexcept:previous_(rangeEnabled){rangeEnabled=enabled;}
HydroFiniteRangeScope::~HydroFiniteRangeScope(){rangeEnabled=previous_;}
bool hydroFiniteRangeEnabled() noexcept{return rangeEnabled;}
void hydroFiniteProductRange(long double a,long double b,long double result)
{if(rangeEnabled && (!std::isfinite(result)||(a!=0&&b!=0&&(result==0||std::fpclassify(result)==FP_SUBNORMAL))))throw HydroFiniteRange{};}
void hydroFiniteQuotientRange(long double a,long double b,long double result)
{if(rangeEnabled && (!std::isfinite(result)||(a!=0&&b!=0&&(result==0||std::fpclassify(result)==FP_SUBNORMAL))))throw HydroFiniteRange{};}
}
