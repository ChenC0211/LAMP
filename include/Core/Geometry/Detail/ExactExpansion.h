#pragma once

#include <algorithm>
#include <array>
#include <cfenv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#if defined(_M_X64) || (defined(__x86_64__) && defined(__SSE2__))
#include <xmmintrin.h>
#endif
#if defined(_MSC_VER) && defined(_M_X64)
#include <intrin.h>
#endif

// Implementation-only engine. Exact sums/products of finite binary64 inputs, stored
// as nonoverlapping low-to-high floating-point expansions. This is not an
// arbitrary-exponent number system: unrepresentable product tails and scale
// changes are rejected, never silently rounded away. Unsafe fast-math or
// reassociation invalidates error-free transforms and is unsupported.
// Implementation-only: the anonymous ancestor deliberately gives this engine
// TU-local linkage. No type or inline definition here is a cross-TU public ABI.
// Keep every include outside these namespaces; public callers use the opaque
// out-of-line snapshot, while legacy tests may inspect this local engine.
namespace LAMP::Geometry::Detail
{
namespace
{
namespace ExactExpansion
{
namespace Detail
{
// Keep exception construction on the failing path. Inlining its string/stack
// machinery into every scalar check prevents the hot arithmetic from inlining.
// This does not disable stack protection or change the diagnostic/exception.
[[noreturn]]
#if defined(_MSC_VER)
__declspec(noinline)
#elif defined(__GNUC__)
__attribute__((noinline,cold))
#endif
inline void fail(const char* message)
{ throw std::runtime_error(std::string("[ExactExpansionCandidate] ")+message); }
inline void require(bool condition,const char* message)
{
    if(!condition) fail(message);
}
inline double checked(double value)
{
    // The engine already requires IEEE binary64. Classify its exponent bits,
    // retaining infinities and every NaN payload as failures, without a CRT
    // classification call for every exact sum/product component. No arithmetic
    // or FP-mode change is involved; return the original finite value bits.
    static_assert(sizeof(double)==sizeof(std::uint64_t)&&std::numeric_limits<double>::is_iec559
        &&std::numeric_limits<double>::digits==53,"exact arithmetic requires IEEE binary64");
    std::uint64_t bits=0;std::memcpy(&bits,&value,sizeof bits);
    require((bits&UINT64_C(0x7ff0000000000000))!=UINT64_C(0x7ff0000000000000),
        "non-finite or overflowing expansion arithmetic");
    return value;
}
struct Access;
}

class Expansion
{
public:
    // Read-only component range, invalidated by mutation/move of its owner.
    // Most error-free transforms have only one or two components. Store a
    // short expansion in the value object, while retaining every component
    // of larger expansions in ordinary owned storage (no capacity cutoff).
    class Components
    {
    public:
        const double* begin() const {return data_;}
        const double* end() const {return data_+size_;}
        std::size_t size() const {return size_;}
        bool empty() const {return size_==0;}
        double back() const {return data_[size_-1];}
        bool operator==(Components other) const
        {return size_==other.size_ && std::equal(begin(),end(),other.begin());}
        bool operator!=(Components other) const {return !(*this==other);}
    private:
        friend class Expansion;
        Components(const double* data,std::size_t size):data_(data),size_(size){}
        const double* data_;
        std::size_t size_;
    };
    Expansion()=default; // Exact zero; no signed-zero distinction.
    Components components() const
    {return spilled_.empty()?Components(inline_.data(),inlineSize_):Components(spilled_.data(),spilled_.size());}
private:
    static constexpr std::size_t inlineCapacity=4;
    std::array<double,inlineCapacity> inline_{};
    std::size_t inlineSize_=0;
    std::vector<double> spilled_;
    friend struct Detail::Access;
};

struct Interval {double lo=0.,hi=0.;};

namespace Detail
{
inline bool& nativeQuarticQueryActive()
{
    // TU-local and actual-thread-local, with stack lifetime only. This flag
    // never grants permission to a later query or another worker thread.
    static thread_local bool active=false;
    return active;
}
inline bool& preparedUnitPolynomialQueryActive()
{
    // Separate permission for one synchronous immutable prepared-polynomial
    // call. No callback, later query, worker or TU inherits this stack scope.
    static thread_local bool active=false;
    return active;
}
}

// Always inspect the actual environment, including inside a nested scope.
inline void validateActualEnvironment()
{
    static_assert(sizeof(double)==sizeof(std::uint64_t)
        &&std::numeric_limits<double>::is_iec559
        &&std::numeric_limits<double>::radix==2&&std::numeric_limits<double>::digits==53,
        "ExactExpansionCandidate requires IEEE binary64 double");
    Detail::require(std::fegetround()==FE_TONEAREST,
        "exact arithmetic requires round-to-nearest; the floating-point mode was not changed");
#if defined(_M_X64) || (defined(__x86_64__) && defined(__SSE2__))
    // Binary64 arithmetic on these targets uses SSE2. Read the CURRENT
    // thread's control register, rather than manufacturing subnormal operands
    // at every exact sum/product/sign check. This is not a cached permission:
    // callbacks may change the mode between queries. Keep the portable
    // rounding check above (including the x87 environment where present).
    // MXCSR: bit 15 flushes tiny results; bit 6 treats tiny inputs as zero.
    constexpr unsigned flushToZero = 1u << 15, denormalsAreZero = 1u << 6;
    Detail::require((_mm_getcsr() & (flushToZero | denormalsAreZero)) == 0,
        "exact arithmetic requires gradual underflow with FTZ/DAZ disabled; the floating-point mode was not changed");
#else
    volatile double minimumNormal=std::numeric_limits<double>::min();
    volatile double half=0.5;
    volatile double subnormal=minimumNormal*half;
    volatile double restored=subnormal+subnormal;
    Detail::require(subnormal>0.&&subnormal<minimumNormal&&restored==minimumNormal,
        "exact arithmetic requires gradual underflow with FTZ/DAZ disabled; the floating-point mode was not changed");
#endif
}
inline void validateEnvironment()
{
    if(!Detail::nativeQuarticQueryActive()&&!Detail::preparedUnitPolynomialQueryActive())validateActualEnvironment();
}
// Only the synchronous, immutable native-Q4 point-query entries use this
// internal scope. They invoke no motion/boundary/user callbacks. Arithmetic,
// component ordering and every finite/range/tail check remain unchanged;
// internal FP checks delegate to the actual entry check until validateReturn.
class NativeQuarticQueryScope
{
public:
    NativeQuarticQueryScope():previous_(Detail::nativeQuarticQueryActive())
    {
        validateActualEnvironment();
        Detail::nativeQuarticQueryActive()=true;
    }
    NativeQuarticQueryScope(const NativeQuarticQueryScope&)=delete;
    NativeQuarticQueryScope& operator=(const NativeQuarticQueryScope&)=delete;
    ~NativeQuarticQueryScope() noexcept
    {Detail::nativeQuarticQueryActive()=previous_;}
    void validateReturn() const {validateActualEnvironment();}
private:
    bool previous_;
};
// Only PreparedUnitPolynomial's constructor/isolate/enclose use this scope.
// Their complete bodies contain pure arithmetic/container work, no callbacks
// or FP-mode mutation. Entry and return inspect the CURRENT actual environment;
// every finite/range/tail guard still executes inside the original kernels.
class PreparedUnitPolynomialQueryScope
{
public:
    PreparedUnitPolynomialQueryScope():previous_(Detail::preparedUnitPolynomialQueryActive())
    {
        validateActualEnvironment();
        Detail::preparedUnitPolynomialQueryActive()=true;
    }
    PreparedUnitPolynomialQueryScope(const PreparedUnitPolynomialQueryScope&)=delete;
    PreparedUnitPolynomialQueryScope& operator=(const PreparedUnitPolynomialQueryScope&)=delete;
    ~PreparedUnitPolynomialQueryScope() noexcept
    {Detail::preparedUnitPolynomialQueryActive()=previous_;}
    void validateReturn() const {validateActualEnvironment();}
private:
    bool previous_;
};

namespace Detail
{
struct Access
{
    static Expansion make(double value)
    {Expansion result;if(value!=0.)push(result,value);return result;}
    // Only internal scratch compaction borrows mutable storage. Its address
    // remains stable until shrink/push; do not retain it past either operation.
    static double* data(Expansion& value)
    {return value.spilled_.empty()?value.inline_.data():value.spilled_.data();}
    // Only append's consumed prefix is resized; it never grows here.
    static void shrink(Expansion& value,std::size_t size)
    {
        if(value.spilled_.empty()) value.inlineSize_=size;
        else value.spilled_.resize(size);
    }
    static void push(Expansion& value,double component)
    {
        if(!value.spilled_.empty()) {value.spilled_.push_back(component);return;}
        if(value.inlineSize_<Expansion::inlineCapacity)
            value.inline_[value.inlineSize_++]=component;
        else
        {
            // Prepare before publishing so allocation failure leaves a valid
            // short expansion; no input operand shares this scratch storage.
            std::vector<double> spilled;
            spilled.reserve(2*Expansion::inlineCapacity);
            spilled.insert(spilled.end(),value.inline_.begin(),value.inline_.end());
            spilled.push_back(component);
            value.spilled_=std::move(spilled);
            value.inlineSize_=0;
        }
    }
};
struct Pair {double head=0.,tail=0.;};
inline Pair twoSum(double a,double b)
{
    const double head=checked(a+b),back=head-a;
    const double tail=checked((a-(head-back))+(b-back));
    return {head,tail};
}
// Exponent of the least significant NONZERO bit, including subnormals.
// For a*b this exponent adds exactly (the integer significands are odd).
inline int leastBitExponent(double value)
{
    std::uint64_t bits=0; std::memcpy(&bits,&value,sizeof bits);
    const int encoded=static_cast<int>((bits>>52)&0x7ffu);
    std::uint64_t significand=bits&UINT64_C(0x000fffffffffffff);
    int exponent=-1074;
    if(encoded!=0) {significand|=UINT64_C(0x0010000000000000);exponent=encoded-1023-52;}
    require(significand!=0,"least-bit exponent requested for zero");
#if defined(_MSC_VER) && defined(_M_X64)
    unsigned long trailingZeros=0;
    _BitScanForward64(&trailingZeros,significand);
    exponent+=static_cast<int>(trailingZeros);
#elif defined(__GNUC__) || defined(__clang__)
    exponent+=__builtin_ctzll(significand);
#else
    while((significand&1u)==0u) {significand>>=1;++exponent;}
#endif
    return exponent;
}
// Exponent of the most significant nonzero bit. A normal binary64's
// encoded exponent already provides ilogb exactly; retain the CRT path for
// subnormals and every value outside that finite-normal case.
inline int mostBitExponent(double value)
{
    std::uint64_t bits=0;std::memcpy(&bits,&value,sizeof bits);
    const int encoded=static_cast<int>((bits>>52)&0x7ffu);
    if(encoded>0&&encoded<0x7ff)return encoded-1023;
    return std::ilogb(value);
}
inline Pair twoProduct(double a,double b)
{
    if(a==0.||b==0.) return {};
    require(leastBitExponent(a)+leastBitExponent(b)>=-1074,
        "exact product tail would underflow binary64; unsupported arithmetic scale");
    const double head=checked(a*b);
    const double tail=checked(std::fma(a,b,-head));
    return {head,tail};
}
inline void append(Expansion& expansion,double value)
{
    if(value==0.) return;
    // Compact in place: at most one tail is written per consumed component,
    // so writes never overtake unread input. Preserve the old TwoSum order
    // and zero elimination exactly; only a possible final head grows storage.
    // This scratch value belongs to sum/product, never to an input Expansion.
    const auto count=expansion.components().size();
    double* const components=Access::data(expansion);
    std::size_t written=0;
    double q=value;
    for(std::size_t i=0;i<count;++i)
    {
        const auto sum=twoSum(q,components[i]);
        if(sum.tail!=0.) components[written++]=sum.tail;
        q=sum.head;
    }
    Access::shrink(expansion,written);
    if(q!=0.) Access::push(expansion,q);
}
inline Expansion sum(const Expansion& a,const Expansion& b,bool subtract)
{
    auto result=a;
    for(double value:b.components()) append(result,subtract?-value:value);
    return result;
}
inline Expansion product(const Expansion& a,const Expansion& b)
{
    Expansion result;
    for(double left:a.components())
        for(double right:b.components())
        {
            const auto p=twoProduct(left,right);
            append(result,p.tail);append(result,p.head);
        }
    return result;
}
}

inline Expansion fromDouble(double value)
{
    validateEnvironment(); Detail::checked(value);
    return Detail::Access::make(value);
}
inline Expansion add(const Expansion& a,const Expansion& b)
{validateEnvironment();return Detail::sum(a,b,false);}
// Accumulate into query-local scratch with the original add component order.
// No published input is mutated: callers discard partial scratch on failure.
// Self-aliasing retains the original term before append can invalidate it.
inline void accumulate(Expansion& target,const Expansion& term)
{
    validateEnvironment();
    if(&target==&term) {
        const auto retained=term;
        for(double value:retained.components())Detail::append(target,value);
    } else {
        for(double value:term.components())Detail::append(target,value);
    }
}
inline Expansion subtract(const Expansion& a,const Expansion& b)
{validateEnvironment();return Detail::sum(a,b,true);}
inline Expansion multiply(const Expansion& a,const Expansion& b)
{validateEnvironment();return Detail::product(a,b);}
inline Expansion scalePow2(const Expansion& value,int exponent)
{
    validateEnvironment();
    Expansion result;
    for(double component:value.components())
    {
        const auto lowest=static_cast<long long>(Detail::leastBitExponent(component))+exponent;
        const auto highest=static_cast<long long>(Detail::mostBitExponent(component))+exponent;
        Detail::require(lowest>=-1074&&highest<=1023,
            "power-of-two scaling would discard or overflow expansion bits; unsupported arithmetic scale");
        std::uint64_t bits=0;std::memcpy(&bits,&component,sizeof bits);
        const int encoded=static_cast<int>((bits>>52)&0x7ffu);
        const auto scaledEncoded=static_cast<long long>(encoded)+static_cast<long long>(exponent);
        double scaled=0.;
        if(encoded>0&&encoded<0x7ff&&scaledEncoded>0&&scaledEncoded<0x7ff)
        {
            // Normal-to-normal power-of-two scaling changes only the encoded
            // exponent. Preserve every sign/significand bit, with no rounding
            // or arithmetic reassociation. Verify the inverse on those bits
            // instead of invoking the CRT for an exact scale and round trip.
            constexpr auto exponentMask=UINT64_C(0x7ff0000000000000);
            const auto scaledBits=(bits&~exponentMask)
                |(static_cast<std::uint64_t>(scaledEncoded)<<52);
            std::memcpy(&scaled,&scaledBits,sizeof scaled);
            scaled=Detail::checked(scaled);
            const auto restoredBits=(scaledBits&~exponentMask)
                |(static_cast<std::uint64_t>(scaledEncoded-static_cast<long long>(exponent))<<52);
            Detail::require(scaled!=0.&&restoredBits==bits,
                "power-of-two scaling would discard expansion bits; unsupported arithmetic scale");
        }
        else
        {
            // The representability bounds constrain exponent to [-2097,
            // 2097], so negating it for this original round trip is safe.
            // Keep the original arithmetic for subnormal sources or results.
            scaled=Detail::checked(std::scalbn(component,exponent));
            Detail::require(scaled!=0.&&std::scalbn(scaled,-exponent)==component,
                "power-of-two scaling would discard expansion bits; unsupported arithmetic scale");
        }
        Detail::Access::push(result,scaled);
    }
    return result;
}
// A homogeneous object has ONE arbitrary positive scale. Center the bit range
// of all its numerators AND denominator together, including every low tail.
// Scaling each coefficient separately would change the represented object.
// This is exponent management, not arbitrary-exponent arithmetic: an object
// or later product wider than binary64's window still fails explicitly.
struct ExponentRange
{
    int lowest=std::numeric_limits<int>::max(),highest=std::numeric_limits<int>::min();
    void include(const Expansion& value)
    {
        for(double component:value.components()) {
            lowest=std::min(lowest,Detail::leastBitExponent(component));
            highest=std::max(highest,Detail::mostBitExponent(component));
        }
    }
    int centeredShift() const
    {
        if(lowest==std::numeric_limits<int>::max())return 0;
        const int lower=-1074-lowest,upper=1023-highest;
        Detail::require(lower<=upper,"homogeneous bit range exceeds binary64 exponent window");
        return std::clamp(-(lowest+highest)/2,lower,upper);
    }
};
template<std::size_t N>
inline void centerHomogeneous(std::array<Expansion,N>& numerator,Expansion& denominator)
{
    ExponentRange range;range.include(denominator);
    for(const auto& value:numerator)range.include(value);
    const int shift=range.centeredShift();
    denominator=scalePow2(denominator,shift);
    for(auto& value:numerator)value=scalePow2(value,shift);
}
inline int sign(const Expansion& value)
{
    validateEnvironment(); // DAZ can otherwise even change a subnormal comparison.
    if(value.components().empty()) return 0;
    return value.components().back()>0.?1:-1;
}
// A finite approximation, NOT a correctly-rounded conversion contract.
// Consumers needing a sign or bound must use sign/enclosure instead.
inline double toDouble(const Expansion& value)
{
    validateEnvironment();
    double result=0.;
    for(double component:value.components()) result=Detail::checked(result+component);
    return result;
}
inline double estimate(const Expansion& value)
{
    const double result=toDouble(value);
    const int expected=sign(value),actual=result>0.?1:result<0.?-1:0;
    Detail::require(actual==expected,"finite estimate would erase or reverse an exact nonzero sign");
    return result;
}
inline Interval enclosure(const Expansion& value)
{
    validateEnvironment();
    Interval result;
    for(double component:value.components())
    {
        const auto lower=Detail::twoSum(result.lo,component);
        const auto upper=Detail::twoSum(result.hi,component);
        result.lo=lower.tail<0.?Detail::checked(std::nextafter(lower.head,-std::numeric_limits<double>::infinity())):lower.head;
        result.hi=upper.tail>0.?Detail::checked(std::nextafter(upper.head,std::numeric_limits<double>::infinity())):upper.head;
    }
    return result;
}
template<std::size_t Count>
inline Expansion evaluate(const std::array<Expansion,Count>& coefficients,double t)
{
    static_assert(Count == 7 || Count == 9 || Count == 12, "unsupported retained polynomial coefficient count");
    validateEnvironment(); Detail::checked(t);
    const auto argument=Detail::Access::make(t);
    Expansion result=coefficients[Count-1];
    for(int k=static_cast<int>(Count)-2;k>=0;--k) result=Detail::sum(Detail::product(result,argument),coefficients[k],false);
    return result;
}
}
} // TU-local engine
} // namespace LAMP::Geometry::Detail
