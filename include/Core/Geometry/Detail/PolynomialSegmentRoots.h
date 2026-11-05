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
#include <vector>

#include "Core/Geometry/Detail/ExactExpansion.h"
#include "Core/Geometry/NativeOutwardRounding.h"

// Implementation-only engine. Bounded floating-point root isolation on [0,1], not
// exact-rational certification. Uses outward-rounded arithmetic enclosures,
// Bernstein sign exclusion and recursively isolated derivative roots. Normal
// IEEE arithmetic/gradual underflow is assumed; unsafe fast-math is unsupported.
// Interior multiple/near-critical roots whose sign cannot be resolved THROW;
// they are never silently called tangencies, deleted, or reported as no roots.
// Implementation-only: the anonymous ancestor deliberately gives this engine
// TU-local linkage. No type or inline definition here is a cross-TU public ABI.
// Keep every include outside these namespaces; public callers use the opaque
// out-of-line snapshot, while legacy tests may inspect this local engine.
namespace LAMP::Geometry::Detail
{
namespace
{
namespace PolynomialSegmentRoots
{
using ExpandedPolynomial=std::array<ExactExpansion::Expansion,7>;
using ExpandedQuarticPolynomial=std::array<ExactExpansion::Expansion,9>;
// Query-private evidence for the explicitly selected engine port. Production
// callers leave this null; no process-global counters or observer are used.
struct GuardedNewtonWork
{
    std::size_t qualifiedRoots=0,suggestions=0,signedSteps=0,enclosedRoots=0,midpointFallbacks=0;
};
struct RootEnclosure {double lo=0.,hi=0.;};
struct Result
{
    std::vector<double> roots;
    bool identicallyZero=false;
    // Same order/count as roots; each representative belongs to its complete
    // enclosure. Phase samples belong between prev.hi and next.lo, NOT merely
    // between representatives, which can leave the sample inside an enclosure.
    std::vector<RootEnclosure> enclosures;
};
namespace Detail
{
inline void require(bool condition,const char* message)
{ if(!condition) throw std::runtime_error(std::string("[PolynomialSegmentRootsCandidate] ")+message); }
inline void checkFloatingPointEnvironment()
{
    require(std::fegetround()==FE_TONEAREST,
        "root isolation requires round-to-nearest; the floating-point mode was not changed");
    // Volatile operands/results force runtime operations, not constant-folded
    // subnormal literals. The first operation detects flushed results (FTZ);
    // reusing that result as both operands detects flushed inputs (DAZ).
    volatile double minimumNormal=std::numeric_limits<double>::min();
    volatile double half=0.5;
    volatile double subnormal=minimumNormal*half;
    volatile double restored=subnormal+subnormal;
    require(subnormal>0.&&subnormal<minimumNormal&&restored==minimumNormal,
        "root isolation requires gradual underflow with FTZ/DAZ disabled; the floating-point mode was not changed");
}
inline double checked(double value)
{ require(std::isfinite(value),"non-finite polynomial arithmetic"); return value; }
inline double outwardNeighbor(double value,bool upward)
{
    checked(value);
    return LAMP::Geometry::nativeOutwardNeighbor(value,upward);
}
inline double down(double value) {return outwardNeighbor(value,false);}
inline double up(double value) {return outwardNeighbor(value,true);}
struct Interval {double lo=0.,hi=0.;};
inline bool zero(Interval value) {return value.lo==0.&&value.hi==0.;}
inline Interval add(Interval a,Interval b)
{
    if(zero(a)) return b; if(zero(b)) return a;
    return {down(a.lo+b.lo),up(a.hi+b.hi)};
}
inline Interval multiply(Interval a,Interval b)
{
    if(zero(a)||zero(b)) return {};
    if(a.lo==1.&&a.hi==1.) return b; if(b.lo==1.&&b.hi==1.) return a;
    // Point/nonnegative factors repeat two of the four endpoint products.
    // Monotonicity selects the SAME extrema; retain checked arithmetic and
    // the original outward rounding, including zero/subnormal boundaries.
    if(b.lo==b.hi&&b.lo>=0.&&a.lo<=a.hi)
    {
        const double lower=checked(a.lo*b.lo),upper=checked(a.hi*b.hi);
        return {down(lower),up(upper)};
    }
    if(a.lo==a.hi&&a.lo>=0.&&b.lo<=b.hi)
    {
        const double lower=checked(a.lo*b.lo),upper=checked(a.lo*b.hi);
        return {down(lower),up(upper)};
    }
    const std::array<double,4> products{{checked(a.lo*b.lo),checked(a.lo*b.hi),checked(a.hi*b.lo),checked(a.hi*b.hi)}};
    return {down(*std::min_element(products.begin(),products.end())),up(*std::max_element(products.begin(),products.end()))};
}
inline Interval divide(Interval a,double positive)
{
    require(positive>0.&&std::isfinite(positive),"invalid interval divisor");
    if(zero(a)||positive==1.) return a;
    return {down(a.lo/positive),up(a.hi/positive)};
}
inline int sign(Interval value) {return value.lo>0.?1:value.hi<0.?-1:0;}
inline double magnitude(Interval value) {return std::max(std::abs(value.lo),std::abs(value.hi));}
inline int choose(int n,int k)
{
    int result=1; for(int i=1;i<=k;++i) result=result*(n-i+1)/i; return result;
}
struct Bracket {double lo=0.,hi=0.;};

// The same isolator handles the legacy degree-six Hermite restriction and
// degree-eight scalar restrictions and degree-eleven native face moments.
// Storage size is a compile-time property; no coefficient cutoff, changed
// enclosure rule or higher iteration budget is selected for the new size.
template<std::size_t CoefficientCount>
class Isolator
{
    static_assert(CoefficientCount==7||CoefficientCount==9||CoefficientCount==12,
        "unsupported retained polynomial coefficient count");
    using ExpandedInput=std::array<ExactExpansion::Expansion,CoefficientCount>;
public:
    explicit Isolator(const std::array<double,CoefficientCount>& input)
    {
        double largest=0.;
        for(double value:input) {require(std::isfinite(value),"polynomial coefficients must be finite");largest=std::max(largest,std::abs(value));}
        if(largest==0.) {degree_=-1;return;}
        int exponent=0; std::frexp(largest,&exponent);
        for(std::size_t i=0;i<input.size();++i)
        {
            coefficients_[i]=std::scalbn(input[i],-exponent);
            require(std::isfinite(coefficients_[i])&&(input[i]==0.||coefficients_[i]!=0.)
                &&std::scalbn(coefficients_[i],exponent)==input[i],
                "coefficient dynamic range cannot be preserved by power-of-two normalization");
        }
        degree_=static_cast<int>(CoefficientCount)-1; while(degree_>0&&coefficients_[degree_]==0.) --degree_;
    }
    explicit Isolator(const ExpandedInput& input,bool endpointRelative=false,
        bool guardedNewton=false,GuardedNewtonWork* guardedWork=nullptr,
        double absoluteResolution=64.*std::numeric_limits<double>::epsilon())
        :usesExpansion_(true),endpointRelative_(endpointRelative),guardedNewton_(guardedNewton),guardedWork_(guardedWork),
        absoluteResolution_(absoluteResolution)
    {
        namespace Exact=ExactExpansion;
        if(absoluteResolution!=64.*std::numeric_limits<double>::epsilon())
            require(std::isfinite(absoluteResolution)&&absoluteResolution>0.
                &&absoluteResolution<=64.*std::numeric_limits<double>::epsilon(),"invalid explicit absolute root resolution");
        Exact::ExponentRange range;
        for(const auto& coefficient:input)range.include(coefficient);
        if(range.lowest==std::numeric_limits<int>::max()) {degree_=-1;return;}
        const int shift=range.centeredShift();
        ExpandedInput normalized;
        for(std::size_t i=0;i<input.size();++i) normalized[i]=Exact::scalePow2(input[i],shift);
        degree_=static_cast<int>(CoefficientCount)-1;while(degree_>0&&Exact::sign(normalized[degree_])==0) --degree_;
        // Cache coefficient enclosures/signs for the ordinary interval path;
        // retain normalized expansions for rare exact point-sign fallback.
        // Endpoint sums retain every expansion component, including exact
        // cancellation that a rounded double conversion would destroy.
        for(int order=0;order<=degree_;++order)
        {
            Exact::Expansion endpoint;
            for(int power=0;power<=degree_-order;++power)
            {
                const auto coefficient=Exact::multiply(normalized[power+order],
                    Exact::fromDouble(falling(power+order,order)));
                const auto interval=Exact::enclosure(coefficient);
                expandedDerivativeCoefficients_[order][power]={interval.lo,interval.hi};
                endpoint=Exact::add(endpoint,coefficient);
            }
            expandedEndpointSigns_[order][0]=Exact::sign(normalized[order]);
            expandedEndpointSigns_[order][1]=Exact::sign(endpoint);
        }
        exactCoefficients_=std::move(normalized);
    }
    Result run() const
    {
        Result result;
        if(degree_<0) {result.identicallyZero=true;return result;}
        const auto brackets=roots(0);
        for(const auto& bracket:brackets)
        {
            require(std::isfinite(bracket.lo)&&std::isfinite(bracket.hi)
                &&bracket.lo>=0.&&bracket.hi<=1.&&bracket.lo<=bracket.hi
                &&(result.enclosures.empty()||result.enclosures.back().hi<bracket.lo),
                "invalid or overlapping public root enclosures");
            const double root=checked(.5*bracket.lo+.5*bracket.hi);
            require(root>=0.&&root<=1.&&(result.roots.empty()||root>result.roots.back()),
                "cannot resolve distinct root locations in the segment");
            require(root>=bracket.lo&&root<=bracket.hi,"root representative lies outside its enclosure");
            result.roots.push_back(root);
            result.enclosures.push_back({bracket.lo,bracket.hi});
        }
        require(result.enclosures.size()==result.roots.size(),"root and enclosure counts disagree");
        return result;
    }
private:
    std::array<double,CoefficientCount> coefficients_{};
    ExpandedInput exactCoefficients_{};
    bool usesExpansion_=false;
    bool endpointRelative_=false;
    bool guardedNewton_=false;
    GuardedNewtonWork* guardedWork_=nullptr;
    double absoluteResolution_=64.*std::numeric_limits<double>::epsilon();
    std::array<std::array<Interval,CoefficientCount>,CoefficientCount> expandedDerivativeCoefficients_{};
    std::array<std::array<int,2>,CoefficientCount> expandedEndpointSigns_{};
    int degree_=0;
    static double falling(int n,int order)
    {double result=1.;for(int k=0;k<order;++k) result*=n-k;return result;}
    Interval coefficient(int power,int order) const
    {
        if(usesExpansion_) return expandedDerivativeCoefficients_[order][power];
        const double value=coefficients_[power+order],factor=falling(power+order,order);
        if(value==0.||factor==1.) return {value,value};
        const double product=checked(value*factor),error=checked(std::fma(value,factor,-product));
        if(error==0.) return {product,product};
        return {error<0.?down(product):product,error>0.?up(product):product};
    }
    Interval evaluate(int order,Interval x) const
    {
        Interval result=coefficient(degree_-order,order);
        for(int k=degree_-order-1;k>=0;--k) result=add(multiply(result,x),coefficient(k,order));
        return result;
    }
    int exactPointSign(int order,double x) const
    {
        namespace E=ExactExpansion;
        const auto coefficientAt=[&](int power) {
            return E::multiply(usesExpansion_?exactCoefficients_[power]:E::fromDouble(coefficients_[power]),
                E::fromDouble(falling(power,order)));
        };
        auto value=coefficientAt(degree_);const auto point=E::fromDouble(x);
        for(int power=degree_-1;power>=order;--power)value=E::add(E::multiply(value,point),coefficientAt(power));
        return E::sign(value);
    }
    // Convex-hull enclosure of the Bernstein coefficients on [a,b]. This is
    // a sufficient no-root test, unlike equal signs at the two endpoints.
    Interval range(int order,double a,double b) const
    {
        if(a==b) return evaluate(order,{a,a});
        const int degree=degree_-order;
        std::array<Interval,CoefficientCount> ap{},wp{},power{};
        ap[0]=wp[0]={1.,1.};
        const Interval width{down(b-a),up(b-a)};
        for(int k=1;k<=degree;++k) {ap[k]=multiply(ap[k-1],{a,a});wp[k]=multiply(wp[k-1],width);}
        for(int j=0;j<=degree;++j)
            for(int k=j;k<=degree;++k)
                power[j]=add(power[j],multiply(multiply(multiply(coefficient(k,order),
                    {static_cast<double>(choose(k,j)),static_cast<double>(choose(k,j))}),ap[k-j]),wp[j]));
        Interval result{std::numeric_limits<double>::infinity(),-std::numeric_limits<double>::infinity()};
        for(int i=0;i<=degree;++i)
        {
            Interval bernstein{};
            for(int j=0;j<=i;++j) bernstein=add(bernstein,divide(multiply(power[j],
                {static_cast<double>(choose(i,j)),static_cast<double>(choose(i,j))}),choose(degree,j)));
            result.lo=std::min(result.lo,bernstein.lo); result.hi=std::max(result.hi,bernstein.hi);
        }
        checked(result.lo);checked(result.hi);return result;
    }
    // Endpoint 1 is a sum of binary coefficients times small exact integers.
    // Error-free sums/FMA product tails recognize its EXACT zero without
    // declaring a merely small Horner residual to be an endpoint root.
    static void appendExpansion(std::vector<double>& expansion,double value)
    {
        std::vector<double> next; next.reserve(expansion.size()+1);
        double q=value;
        for(double item:expansion)
        {
            const double sum=checked(q+item),back=sum-q;
            const double error=(q-(sum-back))+(item-back);
            if(error!=0.) next.push_back(error); q=sum;
        }
        if(q!=0.) next.push_back(q); expansion.swap(next);
    }
    int endpointSign(int order,bool upper) const
    {
        if(usesExpansion_) return expandedEndpointSigns_[order][upper?1:0];
        if(!upper) return coefficients_[order]>0.?1:coefficients_[order]<0.?-1:0;
        std::vector<double> expansion;
        for(int k=order;k<=degree_;++k)
        {
            const double factor=falling(k,order),product=checked(coefficients_[k]*factor);
            const double error=checked(std::fma(coefficients_[k],factor,-product));
            appendExpansion(expansion,error);appendExpansion(expansion,product);
        }
        return expansion.empty()?0:expansion.back()>0.?1:-1;
    }
    Bracket checkedBracket(int order,double lo,double hi) const
    {
        require(lo>=0.&&hi<=1.&&lo<=hi,"invalid root enclosure");
        const auto slope=range(order+1,lo,hi);
        require(sign(slope)!=0,"cannot resolve root near a critical point");
        const double middle=.5*lo+.5*hi;
        const auto residual=evaluate(order,{middle,middle});
        const double bound=up(magnitude(slope)*up(hi-lo));
        require(residual.lo<=bound&&residual.hi>=-bound,"root enclosure fails its polynomial residual check");
        constexpr double nativeResolution=64.*std::numeric_limits<double>::epsilon();
        if(!endpointRelative_&&absoluteResolution_!=nativeResolution&&order==0&&up(hi-lo)>nativeResolution)
        {
            // The original MVT certificate can return before the requested
            // absolute width. Preserve its actual interval and polynomial;
            // do not inherit signs from the earlier monotonic interval.
            const int left=exactPointSign(order,lo),right=exactPointSign(order,hi);
            if(left==0)return checkedBracket(order,lo,lo);
            if(right==0)return checkedBracket(order,hi,hi);
            require(left!=right,"absolute root refinement lacks opposite exact endpoint signs");
            for(int iteration=0;iteration<256;++iteration)
            {
                if(up(hi-lo)<=nativeResolution)return checkedBracket(order,lo,hi);
                const double midpoint=.5*lo+.5*hi;
                require(midpoint>lo&&midpoint<hi,"absolute root refinement stalled above native resolution");
                const int current=exactPointSign(order,midpoint);
                if(current==0)return checkedBracket(order,midpoint,midpoint);
                if(current==left)lo=midpoint;else hi=midpoint;
            }
            require(up(hi-lo)<=nativeResolution,"absolute root refinement exceeded its iteration bound");
            return checkedBracket(order,lo,hi);
        }
        return {lo,hi};
    }
    Bracket simpleRoot(int order,double lo,double hi,int leftSign) const
    {
        const double initialLo=lo,initialHi=hi;
        // The actual retained polynomial's strict Bernstein derivative hull
        // proves this original signed interval monotone. Higher derivative
        // inventories and all unqualified roots keep the original engine.
        const bool guarded=guardedNewton_ && order==0 && endpointRelative_
            && sign(range(order+1,initialLo,initialHi))!=0;
        if(guarded && guardedWork_)++guardedWork_->qualifiedRoots;
        double lastSignedPoint=guarded?.5*lo+.5*hi:0.;
        constexpr double resolution=64.*std::numeric_limits<double>::epsilon();
        for(int iteration=0;iteration<256;++iteration)
        {
            const double middle=.5*lo+.5*hi;
            if(endpointRelative_)
            {
                // Local charts can contain resolved roots far below one ULP
                // of a global coordinate. Refine relative to the nearest
                // endpoint without inventing an absolute resolution floor.
                // The same policy applies to recursively isolated critical
                // roots. Exact endpoint roots are handled separately in roots().
                const double endpointScale=std::min(middle,1.-middle);
                const double localResolution=resolution*endpointScale;
                require(endpointScale>0.&&localResolution>0.&&std::isfinite(localResolution),
                    "endpoint-relative root resolution is unrepresentable");
                if(hi-lo<=localResolution||middle==lo||middle==hi) return checkedBracket(order,lo,hi);
            }
            else if(hi-lo<=absoluteResolution_||middle==lo||middle==hi) return checkedBracket(order,lo,hi);
            if(guarded)
            {
                // Newton supplies a query point only. The current bracket is
                // changed exclusively by the SAME interval sign, or the
                // original residual/slope containment and checkedBracket.
                const auto previousValue=evaluate(order,{lastSignedPoint,lastSignedPoint});
                const auto previousSlope=evaluate(order+1,{lastSignedPoint,lastSignedPoint});
                const double residual=.5*previousValue.lo+.5*previousValue.hi;
                const double derivative=.5*previousSlope.lo+.5*previousSlope.hi;
                const double suggestion=sign(previousSlope)!=0 && derivative!=0
                    ?lastSignedPoint-residual/derivative:middle;
                if(std::isfinite(suggestion) && suggestion>lo && suggestion<hi && suggestion!=middle)
                {
                    if(guardedWork_)++guardedWork_->suggestions;
                    const auto proposedValue=evaluate(order,{suggestion,suggestion});
                    const int proposedSign=sign(proposedValue);
                    if(proposedSign!=0)
                    {
                        if(proposedSign==leftSign)lo=suggestion;else hi=suggestion;
                        lastSignedPoint=suggestion;
                        if(guardedWork_)++guardedWork_->signedSteps;
                        continue;
                    }
                    // This is the original MVT certificate's first local
                    // interval, not a Newton residual stopping criterion.
                    const auto slope=range(order+1,initialLo,initialHi);
                    const double lowerSlope=std::min(std::abs(slope.lo),std::abs(slope.hi));
                    const double rawRadius=lowerSlope>0?magnitude(proposedValue)/lowerSlope
                        :std::numeric_limits<double>::infinity();
                    if(sign(slope)!=0 && std::isfinite(rawRadius))
                    {
                        const double radius=up(rawRadius);
                        const double rawLo=suggestion-radius,rawHi=suggestion+radius;
                        if(std::isfinite(rawLo) && std::isfinite(rawHi))
                        {
                            const double nextLo=down(rawLo),nextHi=up(rawHi);
                            if(nextLo>initialLo && nextHi<initialHi)
                            {
                                const double enclosedLo=std::max(lo,nextLo),enclosedHi=std::min(hi,nextHi);
                                if(enclosedLo<enclosedHi)
                                {
                                    // A rejected candidate certificate falls
                                    // straight through to the old midpoint.
                                    // This checker has no source/FP operations;
                                    // their validation is outside this scope.
                                    try {
                                        const auto bracket=checkedBracket(order,enclosedLo,enclosedHi);
                                        if(guardedWork_)++guardedWork_->enclosedRoots;
                                        return bracket;
                                    } catch(const std::runtime_error&) {}
                                }
                            }
                        }
                    }
                }
                if(guardedWork_)++guardedWork_->midpointFallbacks;
            }
            const auto value=evaluate(order,{middle,middle}); const int current=sign(value);
            if(current!=0)
            {
                if(current==leftSign) lo=middle;else hi=middle;
                if(guarded)lastSignedPoint=middle;
                continue;
            }
            // An entire monotonic interval can touch a derivative-root
            // enclosure at either end. Its Bernstein derivative bound may
            // therefore contain zero even at a perfectly regular middle root.
            // Seek a LOCAL nonzero-slope certificate instead of mistaking that
            // endpoint uncertainty for an unresolved multiple root.
            // The signed bisection bracket may already be narrower than the
            // point-evaluation roundoff enclosure. Start this certificate in
            // the original monotonic interval, not that narrowed bracket;
            // otherwise further shrinking could never contain the bound.
            double localLo=initialLo,localHi=initialHi;
            for(int localIteration=0;localIteration<256;++localIteration)
            {
                const auto slope=range(order+1,localLo,localHi);
                if(sign(slope)!=0)
                {
                    const double lowerSlope=std::min(std::abs(slope.lo),std::abs(slope.hi));
                    const double radius=up(magnitude(value)/lowerSlope);
                    const double nextLo=down(middle-radius),nextHi=up(middle+radius);
                    if(nextLo>localLo&&nextHi<localHi&&nextLo>initialLo&&nextHi<initialHi)
                    {
                        // Strict containment proves opposite endpoint signs
                        // on the local interval by the mean-value theorem.
                        // Hence a unique root exists there and lies in this
                        // whole residual/slope enclosure; no zero is guessed.
                        // Signed bisection has already enclosed this same
                        // unique root. Intersect both certificates: replacing
                        // the bracket by the MVT bound can widen it severely
                        // when the local slope lower bound is small, making
                        // a recursively isolated critical point unusable.
                        const double enclosedLo=std::max(lo,nextLo);
                        const double enclosedHi=std::min(hi,nextHi);
                        require(enclosedLo<enclosedHi,"inconsistent signed and residual root enclosures");
                        return checkedBracket(order,enclosedLo,enclosedHi);
                    }
                }
                const double nextLocalLo=.5*localLo+.5*middle;
                const double nextLocalHi=.5*middle+.5*localHi;
                if(nextLocalLo<=localLo||nextLocalHi>=localHi
                    ||nextLocalLo>=middle||nextLocalHi<=middle) break;
                localLo=nextLocalLo;localHi=nextLocalHi;
            }
            // A regular root only one or two representable points from an
            // endpoint can defeat the residual/slope enclosure: its roundoff
            // radius crosses the endpoint even though the polynomial changes
            // sign strictly inside. Retain the original cheap certificate
            // first; only on its failure use the exact stored coefficients to
            // continue signed bisection. No coefficient or root is discarded.
            const int exact=exactPointSign(order,middle);
            if(exact!=0) {
                if(exact==leftSign)lo=middle;else hi=middle;
                if(guarded)lastSignedPoint=middle;
                continue;
            }
            if(sign(range(order+1,middle,middle))!=0)return checkedBracket(order,middle,middle);
            throw std::runtime_error("[PolynomialSegmentRootsCandidate] cannot resolve root sign near a multiple or critical root");
        }
        throw std::runtime_error("[PolynomialSegmentRootsCandidate] bounded root refinement did not converge");
    }
    std::vector<Bracket> roots(int order) const
    {
        const int degree=degree_-order;
        if(degree==0||sign(range(order,0.,1.))!=0) return {};
        const int firstSign=endpointSign(order,false),lastSign=endpointSign(order,true);
        std::vector<Bracket> result;
        if(firstSign==0) result.push_back({0.,0.});
        const auto critical=degree>1?roots(order+1):std::vector<Bracket>{};
        auto monotonic=[&](double a,double b,int sa,int sb)
        {
            if(a==b) return;
            require(!(sa==0&&sb==0),"unresolved polynomial with zero endpoints and no isolated critical point");
            if(sa!=0&&sb!=0&&sa!=sb) result.push_back(simpleRoot(order,a,b,sa));
        };
        double previous=0.;int previousSign=firstSign;
        for(const auto& bracket:critical)
        {
            if(bracket.lo==bracket.hi&&(bracket.lo==0.||bracket.lo==1.)) continue;
            require(bracket.lo>=previous&&bracket.lo<bracket.hi,"cannot resolve overlapping critical-root enclosures");
            const int criticalSign=sign(range(order,bracket.lo,bracket.hi));
            require(criticalSign!=0,"cannot resolve an interior multiple, tangential, or near-critical root");
            monotonic(previous,bracket.lo,previousSign,criticalSign);
            previous=bracket.hi;previousSign=criticalSign;
        }
        monotonic(previous,1.,previousSign,lastSign);
        if(lastSign==0) result.push_back({1.,1.});
        require(result.size()<=static_cast<std::size_t>(degree),"isolated root count exceeds polynomial degree");
        for(std::size_t i=1;i<result.size();++i)
            require(result[i-1].hi<result[i].lo,"cannot resolve distinct root enclosures without merging phase intervals");
        return result;
    }
};
}

// Coefficients are ascending powers c[k]*t^k. Only EXACT trailing zeros reduce
// degree; no small nonzero high-order term is dropped. Roots are ordered
// representatives of disjoint internally checked enclosures. Exact endpoint
// roots are retained, including endpoint multiplicity; interior unresolvable
// multiplicity is an explicit unsupported case. No shape/observed-error tuning.
inline Result isolate(const std::array<double,7>& coefficients)
{
    Detail::checkFloatingPointEnvironment();
    return Detail::Isolator<7>(coefficients).run();
}
inline Result isolate(const ExpandedPolynomial& coefficients)
{
    ExactExpansion::validateEnvironment();
    return Detail::Isolator<7>(coefficients).run();
}
// Explicit local-chart policy; the default isolate contract is unchanged.
// Keeps the original signed brackets, residual/slope enclosure certificates
// and bounded iteration count. The MVT roundoff enclosure can terminate first;
// this does not promise arbitrary relative accuracy for ill-conditioned roots.
inline Result isolateEndpointRelative(const ExpandedPolynomial& coefficients)
{
    ExactExpansion::validateEnvironment();
    return Detail::Isolator<7>(coefficients,true).run();
}
// Explicit degree-eight ports share every normalization, enclosure and failure
// rule above. The original non-template degree-six overloads remain available.
inline Result isolate(const std::array<double,9>& coefficients)
{
    Detail::checkFloatingPointEnvironment();
    return Detail::Isolator<9>(coefficients).run();
}
inline Result isolate(const ExpandedQuarticPolynomial& coefficients)
{
    ExactExpansion::validateEnvironment();
    return Detail::Isolator<9>(coefficients).run();
}
inline Result isolateEndpointRelative(const ExpandedQuarticPolynomial& coefficients)
{
    ExactExpansion::validateEnvironment();
    return Detail::Isolator<9>(coefficients,true).run();
}
// Explicit native-transverse candidate. Complete inventories, coefficient
// normalization, interval/exact signs, MVT checks and all work/FP gates remain
// those of the original engine; default public ports never select this policy.
inline Result isolateEndpointRelativeGuardedNewton(const ExpandedQuarticPolynomial& coefficients,
    GuardedNewtonWork* work=nullptr)
{
    ExactExpansion::validateEnvironment();
    return Detail::Isolator<9>(coefficients,true,true,work).run();
}
inline Result isolateEndpointRelative(const std::array<ExactExpansion::Expansion,12>& coefficients)
{
    ExactExpansion::validateEnvironment();
    return Detail::Isolator<12>(coefficients,true).run();
}
}
} // TU-local engine
} // namespace LAMP::Geometry::Detail
