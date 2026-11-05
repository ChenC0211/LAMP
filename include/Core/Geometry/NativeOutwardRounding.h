#pragma once

#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <type_traits>

#if defined(_MSC_VER)
#define LAMP_NATIVE_OUTWARD_INLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__)
#define LAMP_NATIVE_OUTWARD_INLINE inline __attribute__((always_inline))
#else
#define LAMP_NATIVE_OUTWARD_INLINE inline
#endif

namespace LAMP::Geometry
{
// Same-type adjacent value toward the requested infinity. Interior binary64
// normal values need only one encoding step and raise no arithmetic status.
// Keep the library's boundary/nonfinite behavior, including errno, exceptions,
// signed zero and NaN payloads. Environment state is not inspected, cleared or
// repaired; the same library fallback may raise its original arithmetic status.
template<class Real> LAMP_NATIVE_OUTWARD_INLINE Real nativeOutwardNeighbor(Real value,bool upward)
{
    static_assert(std::is_same_v<Real,double>||std::is_same_v<Real,long double>);
    if constexpr(sizeof(Real)==sizeof(std::uint64_t)&&std::numeric_limits<Real>::is_iec559
        &&std::numeric_limits<Real>::radix==2&&std::numeric_limits<Real>::digits==53)
    {
        std::uint64_t bits=0;
        std::memcpy(&bits,&value,sizeof bits);
        const unsigned exponent=static_cast<unsigned>((bits>>52)&0x7ffu);
        if(exponent>1&&exponent<0x7feu)
        {
            const bool negative=(bits>>63)!=0;
            if(upward!=negative)++bits;else --bits;
            Real result;
            std::memcpy(&result,&bits,sizeof result);
            return result;
        }
    }
    return std::nextafter(value,upward?std::numeric_limits<Real>::infinity()
        :-std::numeric_limits<Real>::infinity());
}
}

#undef LAMP_NATIVE_OUTWARD_INLINE
