#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Geometry/Detail/ExactFaceTimeRational.h"
#include "Core/Geometry/QuarticBernsteinRoots.h"
#include "Core/Geometry/UnitPolynomialRoots.h"
#include "Core/Geometry/Detail/LevelSetHermiteGeometry.h"
#include "Core/Geometry/Detail/LevelSetHermiteQuarticGridData.h"
#include "Core/Geometry/LevelSetMappedSegment.h"
#include "Core/Geometry/Detail/LevelSetMappedSegment.h"
#include "Core/Geometry/Detail/LevelSetMappedSegmentSensitivity.h"
#include "Core/Geometry/LevelSetMappedInterface.h"
#include "Core/Geometry/Detail/LevelSetMappedInterface.h"
#include "Core/Geometry/Detail/LevelSetNativeInterfaceChain.h"
#include "Core/Geometry/Detail/LevelSetInterfaceQuadrature.h"

#include <cmath>
#include <array>
#include <cstdint>
#include <cstring>
#include <exception>
#include <functional>
#include <map>
#include <mutex>
#include <stdexcept>
#include <sstream>
#include <iomanip>
#include <string>
#include <utility>
#include <variant>
#include <type_traits>

namespace Quantum::Geometry
{
static std::array<Detail::ExactExpansion::Expansion,12> expandedUnitPolynomial(const std::vector<UnitPolynomialTerm>& terms)
{
    namespace E=Detail::ExactExpansion;E::validateEnvironment();
    std::array<Detail::ExactExpansion::Expansion,12> polynomial{};
    for (const auto& term:terms) {
        E::Expansion factor;
        long double remainder=term.factor;
        if (!std::isfinite(remainder)) throw std::invalid_argument("nonfinite unit polynomial coefficient");
        for (int part=0;remainder!=0 && part<4;++part) {
            const double value=static_cast<double>(remainder);
            if (!std::isfinite(value) || value==0)
                throw std::invalid_argument("unit polynomial coefficient is outside binary64 expansion range");
            factor=E::add(factor,E::fromDouble(value));
            remainder-=static_cast<long double>(value);
        }
        if (remainder!=0) throw std::invalid_argument("unit polynomial coefficient expansion is unresolved");
        if (!std::isfinite(term.weight)) throw std::invalid_argument("nonfinite unit polynomial weight");
        factor=E::multiply(factor,E::fromDouble(term.weight));
        for (std::size_t i=0;i<term.coefficients.size();++i) {
            if (!std::isfinite(term.coefficients[i])) throw std::invalid_argument("nonfinite unit polynomial basis coefficient");
            polynomial[i]=E::add(polynomial[i],E::multiply(factor,E::fromDouble(term.coefficients[i])));
        }
    }
    return polynomial;
}
UnitPolynomialRoots isolateUnitPolynomial(const std::vector<UnitPolynomialTerm>& terms)
{
    const auto polynomial=expandedUnitPolynomial(terms);
    const auto isolated=Detail::PolynomialSegmentRoots::isolateEndpointRelative(polynomial);
    UnitPolynomialRoots result;result.identicallyZero=isolated.identicallyZero;
    for (std::size_t i=0;i<isolated.roots.size();++i)
        result.roots.push_back({isolated.roots[i],isolated.enclosures[i].lo,isolated.enclosures[i].hi});
    return result;
}
long double evaluateUnitPolynomial(const std::vector<UnitPolynomialTerm>& terms,double parameter)
{
    namespace E=Detail::ExactExpansion;
    if(!std::isfinite(parameter) || parameter<0 || parameter>1)
        throw std::invalid_argument("unit polynomial value parameter outside [0,1]");
    const auto exact=E::evaluate(expandedUnitPolynomial(terms),parameter);long double value=0;
    for(double part:exact.components())value+=part;
    if(!std::isfinite(value) || (E::sign(exact)!=0 && value==0))
        throw std::invalid_argument("unit polynomial value is unrepresentable");
    return value;
}
std::array<double,2> encloseUnitPolynomial(const std::vector<UnitPolynomialTerm>& terms,double a,double b)
{
    namespace E=Detail::ExactExpansion;namespace B=Detail::PolynomialSegmentRoots::Detail;
    if(!std::isfinite(a) || !std::isfinite(b) || a<0 || b>1 || a>b)
        throw std::invalid_argument("unit polynomial range outside [0,1]");
    const auto polynomial=expandedUnitPolynomial(terms);
    if(a==b){const auto value=E::enclosure(E::evaluate(polynomial,a));return {value.lo,value.hi};}
    int degree=11;while(degree>0 && E::sign(polynomial[degree])==0)--degree;
    std::array<E::Expansion,12> origin{},width{},power{};
    origin[0]=width[0]=E::fromDouble(1);
    const auto span=E::subtract(E::fromDouble(b),E::fromDouble(a));
    for(int k=1;k<=degree;++k){origin[k]=E::multiply(origin[k-1],E::fromDouble(a));width[k]=E::multiply(width[k-1],span);}
    for(int j=0;j<=degree;++j)for(int k=j;k<=degree;++k)
        power[j]=E::add(power[j],E::multiply(E::multiply(E::multiply(polynomial[k],E::fromDouble(B::choose(k,j))),origin[k-j]),width[j]));
    // Common denominators preserve cancellation during the power-to-Bernstein
    // conversion. Only the final positive division uses outward rounding.
    constexpr int denominator[]{1,1,2,3,12,10,60,105,280,252,2520,2310};const int scale=denominator[degree];
    double lo=std::numeric_limits<double>::infinity(),hi=-lo;
    for(int i=0;i<=degree;++i){E::Expansion numerator;
        for(int j=0;j<=i;++j)numerator=E::add(numerator,E::multiply(power[j],E::fromDouble(scale*B::choose(i,j)/B::choose(degree,j))));
        const auto exact=E::enclosure(numerator);const auto bound=B::divide({exact.lo,exact.hi},scale);
        lo=std::min(lo,bound.lo);hi=std::max(hi,bound.hi);
    }
    return {lo,hi};
}
namespace {
// An exact polynomial antiderivative has no quadrature error. Measure the
// rounding of its represented quotient against the retained expansion,
// instead of assigning a relative error even to an exact constant integral.
Detail::ExactExpansion::Expansion expandedIntegralValue(long double value) {
    namespace E=Detail::ExactExpansion;E::Expansion result;
    while(value!=0) {
        const double part=static_cast<double>(value);
        if(!std::isfinite(part)||part==0)throw std::invalid_argument("unrepresentable integral expansion");
        result=E::add(result,E::fromDouble(part));value-=part;
    }
    return result;
}
long double integralRoundingError(const Detail::ExactExpansion::Expansion& numerator,
    double denominator,long double value) {
    namespace E=Detail::ExactExpansion;
    const auto residual=E::subtract(numerator,E::multiply(expandedIntegralValue(value),E::fromDouble(denominator)));
    if(E::sign(residual)==0)return 0;
    const auto bound=E::enclosure(residual);
    return std::nextafter(static_cast<long double>(std::max(std::abs(bound.lo),std::abs(bound.hi)))/denominator,
        std::numeric_limits<long double>::infinity());
}
UnitPolynomialIntegral signedUnitPolynomialIntegral(const std::vector<UnitPolynomialTerm>& terms,double a0,double a1,double s0,double s1)
{
    namespace E=Detail::ExactExpansion;
    for(double v:{a0,a1,s0,s1}) if(!std::isfinite(v) || v<0 || v>1)
        throw std::invalid_argument("unit polynomial integral parameter outside [0,1]");
    const auto polynomial=expandedUnitPolynomial(terms);
    const auto origin=E::fromDouble(a0),width=E::subtract(E::fromDouble(a1),origin);
    const auto t0=E::add(origin,E::multiply(width,E::fromDouble(s0)));
    const auto t1=E::add(origin,E::multiply(width,E::fromDouble(s1)));
    auto power0=E::fromDouble(1),power1=power0;E::Expansion numerator;
    // Retain original low-degree denominators; degree eleven uses lcm(1,...,12).
    int degree=11;while(degree>0 && E::sign(polynomial[degree])==0)--degree;
    // Do not form unused high powers of a tiny endpoint. Exact trailing
    // zeros reduce the polynomial degree; no small coefficient is discarded.
    for(int k=0;k<=degree;++k) {
        power0=E::multiply(power0,t0);power1=E::multiply(power1,t1);
        numerator=E::add(numerator,E::multiply(E::multiply(polynomial[k],E::fromDouble((degree<=6?420:(degree<=9?2520:27720))/(k+1))),E::subtract(power1,power0)));
    }
    long double value=0;for(double part:numerator.components()) value+=part;
    value/=(degree<=6?420:(degree<=9?2520:27720));
    if(!std::isfinite(value) || (E::sign(numerator)!=0 && value==0))
        throw std::invalid_argument("unit polynomial integral is unrepresentable");
    return {value,integralRoundingError(numerator,degree<=6?420:(degree<=9?2520:27720),value)};
}
}
long double integrateUnitPolynomial(const std::vector<UnitPolynomialTerm>& terms,double a0,double a1,double s0,double s1)
{
    return signedUnitPolynomialIntegral(terms,a0,a1,s0,s1).value;
}
namespace {
void integrateBernsteinRange(const std::vector<UnitPolynomialTerm>& terms,int degree,int first,int last,int half,double lo,double hi,UnitPolynomialIntegral* result){
    namespace E=Detail::ExactExpansion;
    if(degree<0||degree>20||first<0||last<first||last>degree||(half!=0&&half!=1)||!std::isfinite(lo)||!std::isfinite(hi)||lo<0||hi>1||lo>hi)
        throw std::invalid_argument("invalid weighted unit polynomial chart");
    // 16 nodes integrate the complete degree <=31 product. Keep the weight
    // factored: expansion into degree-31 powers creates huge cancellation
    // and unsupported tiny exact product tails without adding information.
    static const auto rule=[](){std::array<std::array<double,2>,16> out{};
        for(int k=0;k<16;k++){double x=std::cos(std::acos(-1.)*(k+.75)/16.5),d=0;
            for(int iteration=0;iteration<32;iteration++){double a=1,b=x;for(int n=2;n<=16;n++){const double c=((2*n-1)*x*b-(n-1)*a)/n;a=b;b=c;}
                d=16*(x*b-a)/(x*x-1);const double next=x-b/d;if(next==x)break;x=next;}
            out[k]={x,2/((1-x*x)*d*d)};}return out;}();
    const auto polynomial=expandedUnitPolynomial(terms);std::array<double,21> binomial{};
    for(int index=first;index<=last;index++){binomial[index]=1;for(int k=1;k<=index;k++)binomial[index]=binomial[index]*(degree-k+1)/k;}
    for(const auto& node:rule){const double t=lo+.5*(hi-lo)*(1+node[0]);const double xi=.5*(half+t);
        const auto evaluated=E::evaluate(polynomial,t);long double q=0;for(double part:evaluated.components())q+=part;
        for(int index=first;index<=last;index++){
            const long double value=.5L*(hi-lo)*node[1]*q*binomial[index]*std::pow(xi,index)*std::pow(1-xi,degree-index);
            auto& out=result[index-first];out.value+=value;out.absoluteError+=std::abs(value);}}
    for(int index=first;index<=last;index++){
        auto& out=result[index-first];out.absoluteError*=64*std::numeric_limits<double>::epsilon();
        if(!std::isfinite(out.value)||!std::isfinite(out.absoluteError))throw std::invalid_argument("weighted unit polynomial integral is unrepresentable");}
}
}
UnitPolynomialIntegral integrateUnitPolynomialBernstein(const std::vector<UnitPolynomialTerm>& terms,int degree,int index,int half,double lo,double hi){
    UnitPolynomialIntegral result;integrateBernsteinRange(terms,degree,index,index,half,lo,hi,&result);return result;
}
std::vector<UnitPolynomialIntegral> integrateUnitPolynomialBernsteinBasis(const std::vector<UnitPolynomialTerm>& terms,int degree,int half,double lo,double hi){
    if(degree<0||degree>20)throw std::invalid_argument("invalid weighted unit polynomial degree");
    std::vector<UnitPolynomialIntegral> result(degree+1);integrateBernsteinRange(terms,degree,0,degree,half,lo,hi,result.data());return result;
}
namespace {
std::array<Detail::ExactExpansion::Expansion,2> mappedUnitPolynomialEndpoints(
    double a0,double a1,bool fromEnd,long double offset,long double length)
{
    namespace E=Detail::ExactExpansion;
    if(!std::isfinite(a0) || !std::isfinite(a1) || a0<0 || a1>1 || a0>a1 ||
        !std::isfinite(offset) || !std::isfinite(length) || offset<0 || length<0 || offset>1 ||
        offset+length>1+8*std::numeric_limits<long double>::epsilon())
        throw std::invalid_argument("invalid mapped polynomial range");
    const auto expanded=[](long double value) {
        E::Expansion result;
        for(int k=0;value!=0 && k<4;++k) {
            const double part=static_cast<double>(value);
            if(!std::isfinite(part) || part==0)throw std::invalid_argument("unrepresentable polynomial range");
            result=E::add(result,E::fromDouble(part));value-=part;
        }
        if(value!=0)throw std::invalid_argument("unresolved polynomial range expansion");
        return result;
    };
    const auto width=E::subtract(E::fromDouble(a1),E::fromDouble(a0));
    auto s0=expanded(offset),s1=E::add(s0,expanded(length));
    if(fromEnd) {const auto one=E::fromDouble(1);const auto low=E::subtract(one,s1);s1=E::subtract(one,s0);s0=low;}
    const auto t0=E::add(E::fromDouble(a0),E::multiply(width,s0));
    const auto t1=E::add(E::fromDouble(a0),E::multiply(width,s1));
    return {{t0,t1}};
}
}
long double integrateUnitPolynomialRange(const std::vector<UnitPolynomialTerm>& terms,
    double a0,double a1,bool fromEnd,long double offset,long double length)
{
    namespace E=Detail::ExactExpansion;
    const auto endpoints=mappedUnitPolynomialEndpoints(a0,a1,fromEnd,offset,length);
    const auto& t0=endpoints[0];const auto& t1=endpoints[1];
    const auto polynomial=expandedUnitPolynomial(terms);
    auto p0=E::fromDouble(1),p1=p0;E::Expansion numerator;
    int degree=11;while(degree>0 && E::sign(polynomial[degree])==0)--degree;
    for(int k=0;k<=degree;++k) {
        p0=E::multiply(p0,t0);p1=E::multiply(p1,t1);
        numerator=E::add(numerator,E::multiply(E::multiply(polynomial[k],E::fromDouble((degree<=6?420:(degree<=9?2520:27720))/(k+1))),E::subtract(p1,p0)));
    }
    long double result=0;for(double part:numerator.components())result+=part;result/=(degree<=6?420:(degree<=9?2520:27720));
    if(!std::isfinite(result) || (E::sign(numerator)!=0 && result==0))
        throw std::invalid_argument("mapped polynomial integral is unrepresentable");
    return result;
}
namespace {
long double positiveLinearAverage(long double a,long double b) {
    if(a>=0 && b>=0)return .5L*a+.5L*b;
    if(a<=0 && b<=0)return 0;
    const auto positive=std::max(a,b);return .5L*positive*(positive/std::abs(b-a));
}
UnitPolynomialIntegral positivePairIntegral(const std::vector<UnitPolynomialTerm>& first,
    const std::vector<UnitPolynomialTerm>& last,double lo,double hi,long double tolerance,int depth) {
    constexpr long double eps=std::numeric_limits<double>::epsilon();
    if(lo==hi)return {};
    const auto a=encloseUnitPolynomial(first,lo,hi),b=&first==&last?a:encloseUnitPolynomial(last,lo,hi);
    const long double width=static_cast<long double>(hi)-lo;
    if(a[1]<=0 && b[1]<=0)return {};
    const auto ip=signedUnitPolynomialIntegral(first,0,1,lo,hi);
    const auto iq=&first==&last?ip:signedUnitPolynomialIntegral(last,0,1,lo,hi);
    const auto signedAverage=&first==&last?ip.value:.5L*ip.value+.5L*iq.value;
    namespace E=Detail::ExactExpansion;
    const auto averageError=&first==&last?ip.absoluteError:
        .5L*ip.absoluteError+.5L*iq.absoluteError+integralRoundingError(
            E::add(expandedIntegralValue(ip.value),expandedIntegralValue(iq.value)),2,signedAverage);
    if(a[0]>=0 && b[0]>=0)return {signedAverage,averageError};
    const double middle=.5*lo+.5*hi;
    if((a[0]>=0 && b[1]<=0)||(a[1]<=0 && b[0]>=0)) {
        // Signs are certified throughout the interval. The time-integrated
        // reversal is a smooth rational function; keep monitored Gauss
        // integration, without first asking for global coordinates of roots.
        constexpr double x[]{-.9061798459386639928,-.5384693101056830910,0,.5384693101056830910,.9061798459386639928};
        constexpr double w[]{.2369268850561890875,.4786286704993664680,.5688888888888888889,.4786286704993664680,.2369268850561890875};
        const auto rule=[&](double lower,double upper){long double sum=0;
            for(int k=0;k<5;k++){const double s=.5*(lower+upper+(upper-lower)*x[k]);
                sum+=w[k]*positiveLinearAverage(evaluateUnitPolynomial(first,s),evaluateUnitPolynomial(last,s));}
            return .5L*(upper-lower)*sum;};
        const auto coarse=rule(lo,hi),fine=rule(lo,middle)+rule(middle,hi);
        const auto error=std::abs(fine-coarse)+8*eps*std::abs(fine);
        if(error<=tolerance)return {fine,error};
    }else{
        // max(q,0) is monotone and convex. Bound its integral, including an
        // arbitrarily thin positive/negative interval missed by quadrature.
        // When the time mean stays nonnegative, only negative endpoint
        // triangles differ from the exact signed integral.
        long double lower=std::max(0.L,signedAverage),upper;
        if(static_cast<long double>(a[0])+b[0]>=0){
            const auto correction=[](long double negative,long double positive){
                return negative==0?0.L:.5L*negative*(negative/(negative+positive));};
            upper=signedAverage+width*(correction(std::max(0.L,-static_cast<long double>(a[0])),std::max(0.L,static_cast<long double>(b[0])))
                +correction(std::max(0.L,-static_cast<long double>(b[0])),std::max(0.L,static_cast<long double>(a[0]))));
        }else upper=width*positiveLinearAverage(a[1],b[1]);
        const auto rounding=averageError+8*eps*(std::abs(lower)+std::abs(upper));
        lower=std::max(0.L,lower-rounding);upper=std::max(lower,upper+rounding);
        const auto error=.5L*(upper-lower);
        if(error<=tolerance)return {.5L*lower+.5L*upper,error};
    }
    if(depth>=64 || middle<=lo || middle>=hi)
        throw std::runtime_error("positive polynomial integral cannot meet its requested flux budget");
    const auto left=positivePairIntegral(first,last,lo,middle,.5L*tolerance,depth+1);
    const auto right=positivePairIntegral(first,last,middle,hi,.5L*tolerance,depth+1);
    const auto value=left.value+right.value;
    const auto rounding=integralRoundingError(E::add(expandedIntegralValue(left.value),expandedIntegralValue(right.value)),1,value);
    const auto error=left.absoluteError+right.absoluteError+rounding;
    return {value,error==0?0:std::nextafter(error,std::numeric_limits<long double>::infinity())};
}
}
UnitPolynomialIntegral integratePositiveTimeLinearUnitPolynomial(const std::vector<UnitPolynomialTerm>& first,
    const std::vector<UnitPolynomialTerm>& last,double lo,double hi,long double tolerance) {
    if(!std::isfinite(lo)||!std::isfinite(hi)||lo<0||hi>1||lo>hi||!std::isfinite(tolerance)||tolerance<0)
        throw std::invalid_argument("invalid positive polynomial integral interval/budget");
    const auto result=positivePairIntegral(first,last,lo,hi,tolerance,0);
    if(!std::isfinite(result.value)||!std::isfinite(result.absoluteError)||result.absoluteError>tolerance) {
        std::ostringstream message;message<<std::setprecision(17)
            <<"positive polynomial integral exceeds its requested flux budget: interval="<<lo<<":"<<hi
            <<" value="<<result.value<<" error="<<result.absoluteError<<" budget="<<tolerance;
        throw std::runtime_error(message.str());
    }
    return result;
}
UnitPolynomialIntegral integratePositiveUnitPolynomial(const std::vector<UnitPolynomialTerm>& terms,
    double lo,double hi,long double tolerance) {
    return integratePositiveTimeLinearUnitPolynomial(terms,terms,lo,hi,tolerance);
}
UnitPolynomialIntegral integratePositiveUnitPolynomialRange(const std::vector<UnitPolynomialTerm>& terms,
    double a0,double a1,bool fromEnd,long double offset,long double length,long double tolerance)
{
    namespace E=Detail::ExactExpansion;
    namespace B=Detail::PolynomialSegmentRoots::Detail;
    const auto endpoints=mappedUnitPolynomialEndpoints(a0,a1,fromEnd,offset,length);
    const auto polynomial=expandedUnitPolynomial(terms);
    if(!std::isfinite(tolerance) || tolerance<0)throw std::invalid_argument("invalid mapped positive polynomial budget");
    const auto width=E::subtract(endpoints[1],endpoints[0]);
    if(E::sign(width)==0)return {};
    int degree=11;while(degree>0 && E::sign(polynomial[degree])==0)--degree;
    std::array<E::Expansion,12> originPower{},widthPower{},power{};
    originPower[0]=widthPower[0]=E::fromDouble(1);
    for(int k=1;k<=degree;++k){originPower[k]=E::multiply(originPower[k-1],endpoints[0]);widthPower[k]=E::multiply(widthPower[k-1],width);}
    // P(t0+width*x)*width, integrated over x in [0,1]. Exact coefficient
    // expansions are re-expressed as retained terms, never rounded/fit.
    for(int j=0;j<=degree;++j)for(int k=j;k<=degree;++k)
        power[j]=E::add(power[j],E::multiply(E::multiply(E::multiply(polynomial[k],
            E::fromDouble(B::choose(k,j))),originPower[k-j]),widthPower[j]));
    std::vector<UnitPolynomialTerm> transformed;
    for(int j=0;j<=degree;++j){const auto coefficient=E::multiply(power[j],width);
        for(double part:coefficient.components()) {
            UnitPolynomialTerm term;term.coefficients[j]=part;term.factor=1;transformed.push_back(term);
        }}
    return integratePositiveUnitPolynomial(transformed,0,1,tolerance);
}
QuarticBernsteinRoots isolateQuarticBernstein(const std::array<double,5>& controls)
{
    namespace E=Detail::ExactExpansion;
    E::validateEnvironment();
    Detail::PolynomialSegmentRoots::ExpandedPolynomial polynomial{};
    constexpr int coefficients[5][5]={{1,-4,6,-4,1},{0,4,-12,12,-4},
        {0,0,6,-12,6},{0,0,0,4,-4},{0,0,0,0,1}};
    for (int i=0;i<5;++i) {
        if (!std::isfinite(controls[i])) throw std::invalid_argument("nonfinite Bernstein control");
        for (int k=i;k<5;++k) polynomial[k]=E::add(polynomial[k],
            E::multiply(E::fromDouble(controls[i]),E::fromDouble(coefficients[i][k])));
    }
    const auto result=Detail::PolynomialSegmentRoots::isolate(polynomial);
    return {result.roots,result.identicallyZero};
}
namespace
{
namespace ChartExact = Detail::ExactExpansion;
void requireChart(bool valid, const char* reason)
{
    if (!valid) throw std::runtime_error(std::string("[CellChartPoint] ") + reason);
}
ChartExact::Expansion coordinateExpansion(CellChartCoordinate coordinate)
{
    requireChart(std::isfinite(coordinate.offset) && coordinate.offset >= 0.0 && coordinate.offset <= 1.0,
        "offset must be finite and in [0,1]");
    const auto offset = ChartExact::fromDouble(coordinate.offset);
    return coordinate.fromUpper ? ChartExact::subtract(ChartExact::fromDouble(1.0), offset) : offset;
}
}

CellChartPoint cellChartFromLocal(Vector2D local)
{
    const CellChartPoint result{{local.x, false}, {local.y, false}};
    validateCellChart(result);
    return result;
}

void validateCellChart(const CellChartPoint& point)
{
    ChartExact::validateEnvironment();
    coordinateExpansion(point.x); coordinateExpansion(point.y);
}

bool cellChartInterior(const CellChartPoint& point)
{
    validateCellChart(point);
    return point.x.offset > 0.0 && point.x.offset < 1.0
        && point.y.offset > 0.0 && point.y.offset < 1.0;
}

int compareCellChartCoordinate(CellChartCoordinate a, CellChartCoordinate b)
{
    ChartExact::validateEnvironment();
    return ChartExact::sign(ChartExact::subtract(coordinateExpansion(a), coordinateExpansion(b)));
}

bool cellChartEqual(const CellChartPoint& a, const CellChartPoint& b)
{
    validateCellChart(a); validateCellChart(b);
    return compareCellChartCoordinate(a.x, b.x) == 0 && compareCellChartCoordinate(a.y, b.y) == 0;
}

bool cellChartEqual(const CellChartPoint& a, Vector2D local)
{
    return cellChartEqual(a, cellChartFromLocal(local));
}

double cellChartCoordinateValue(CellChartCoordinate coordinate)
{
    ChartExact::validateEnvironment();
    const auto exact = coordinateExpansion(coordinate);
    const double value = ChartExact::estimate(exact);
    requireChart(std::isfinite(value)
        && ChartExact::sign(ChartExact::subtract(exact, ChartExact::fromDouble(value))) == 0,
        "coordinate is not exactly representable as one double; retain its endpoint-relative chart");
    return value;
}

double cellChartIntervalFraction(CellChartCoordinate coordinate, double begin, double end)
{
    ChartExact::validateEnvironment();
    const auto value = coordinateExpansion(coordinate);
    const auto first = coordinateExpansion({begin, false}), last = coordinateExpansion({end, false});
    const auto span = ChartExact::subtract(last, first), numerator = ChartExact::subtract(value, first);
    const int direction = ChartExact::sign(span);
    requireChart(direction != 0, "graph parameter interval must have distinct endpoints");
    const int firstSign = ChartExact::sign(numerator);
    const int lastSign = ChartExact::sign(ChartExact::subtract(last, value));
    requireChart(firstSign * direction >= 0 && lastSign * direction >= 0,
        "chart coordinate is outside the original graph interval");
    if (firstSign == 0) return 0.0;
    if (lastSign == 0) return 1.0;
    const auto fraction = Detail::MappedHermiteSegment::ratio(numerator, span);
    const double result = static_cast<double>(Detail::MappedHermiteSegment::approximate(fraction));
    requireChart(std::isfinite(result) && result > 0.0 && result < 1.0,
        "positive interior graph parameter is not representable without endpoint collapse");
    return result;
}
}

// This TU has an isolated strict FP build contract. The Detail engine has
// TU-local linkage, so precise legacy-test instantiations cannot replace its
// exact-arithmetic functions through external-inline/COMDAT selection.
// Existing Q4/Eigen construction remains shared approximate DOF generation;
// it is not advertised as part of that isolated exact-arithmetic closure.
namespace LevelSetHermiteGeometry
{
namespace
{
namespace Engine = Quantum::Geometry::Detail::LevelSetHermiteGeometry;
namespace Exact = Quantum::Geometry::Detail::ExactExpansion;

Frame resultFrame(const Engine::Frame& frame)
{
    return {frame.point, frame.normal, frame.curvature, frame.levelSetResidual,
        frame.stationarityResidual, frame.iterations, frame.attemptedSeeds,
        frame.convergedSeeds, frame.cellLocalPoint, frame.cellI, frame.cellJ, frame.hasCellLocalPoint,
        frame.cellChartPoint, frame.hasCellChartPoint};
}
}

Quantum::Geometry::CellChartPoint cellChartPoint(const Frame& frame)
{
    if (frame.hasCellLocalPoint == frame.hasCellChartPoint)
        throw std::runtime_error("[LevelSetHermiteGeometry] frame requires exactly one explicit local provenance representation");
    if (frame.hasCellLocalPoint) return Quantum::Geometry::cellChartFromLocal(frame.cellLocalPoint);
    Quantum::Geometry::validateCellChart(frame.cellChartPoint);
    return frame.cellChartPoint;
}

struct Snapshot::Impl
{
    using NativeScalar = Quantum::Geometry::Detail::LevelSetHermiteScalar::LeafSnapshot;
    using NativeQuarticScalar = Quantum::Geometry::Detail::LevelSetHermiteScalar::QuarticLeafData;
    using QuarticScalar = Quantum::Geometry::Detail::LevelSetHermiteScalar::QuarticGridData;
    using QuarticGeometry = Engine::BasicSnapshot<QuarticScalar>;
    std::unique_ptr<Engine::Snapshot> uniform;
    std::unique_ptr<QuarticGeometry> quartic;
    std::unique_ptr<NativeScalar> native;
    std::unique_ptr<NativeQuarticScalar> nativeQuartic;
    Impl(const QuadtreeMesh& mesh, const VertField<double>& phi, bool radialAxisSymmetry,
        CartesianGridScope scope, Reconstruction reconstruction)
    {
        switch (reconstruction)
        {
        case Reconstruction::Cubic:
            uniform = std::make_unique<Engine::Snapshot>(mesh,phi,radialAxisSymmetry,scope);
            break;
        case Reconstruction::EnrichedQuartic:
            quartic = std::make_unique<QuarticGeometry>(QuarticScalar::fromEnrichedQ4(mesh,phi,radialAxisSymmetry,scope));
            break;
        default:
            throw std::invalid_argument("[LevelSetHermiteGeometry] invalid scalar reconstruction");
        }
    }
    Impl(const QuadtreeMesh& mesh, const Quantum::Geometry::ActiveHermitePatch& patch,
        bool radialAxisSymmetry, CartesianGridScope scope)
        : uniform(std::make_unique<Engine::Snapshot>(mesh,patch,radialAxisSymmetry,scope)) {}
    Impl(const QuadtreeMesh& mesh, const Quantum::Geometry::ActiveHermitePatch& patch, bool radialAxisSymmetry)
        : native(std::make_unique<NativeScalar>(mesh,patch,radialAxisSymmetry)) {}
    Impl(const QuadtreeMesh& mesh,const Quantum::Geometry::ActiveQuarticHermitePatch& patch,bool radialAxisSymmetry)
        : nativeQuartic(std::make_unique<NativeQuarticScalar>(mesh,patch,radialAxisSymmetry)) {}
    Impl(const QuadtreeMesh& mesh,const NativeQuarticScalar& source)
        : nativeQuartic(NativeQuarticScalar::refined(mesh,source)) {}
    Impl(const QuadtreeMesh& mesh,const Quantum::Geometry::ActiveQuarticHermitePatch& patch,
        const std::vector<int>& cells,bool radialAxisSymmetry)
        : nativeQuartic(std::make_unique<NativeQuarticScalar>(mesh,patch,radialAxisSymmetry,&cells)) {}
    bool usesNative() const { return native!=nullptr || nativeQuartic!=nullptr; }
    // The legacy mapped entry owns 16 corner-jet sensitivities. Quartic
    // curves must use the explicitly typed 25-column integration entry.
    const Engine::Snapshot& cubicGeometry() const
    {
        if (quartic || nativeQuartic) throw std::runtime_error(
            "[LevelSetHermiteGeometry] EnrichedQuartic requires the typed 25-column mapped integration entry");
        if (!uniform) throw std::runtime_error(
            "[LevelSetHermiteGeometry] legacy logical-grid/global-jet/projection operation is not connected for native leaves");
        return *uniform;
    }
    template<class Query> decltype(auto) visitGeometry(Query&& query) const
    {
        if (quartic) return query(*quartic);
        return query(cubicGeometry());
    }
    const NativeScalar& leaf() const
    {
        if (nativeQuartic) throw std::runtime_error(
            "[LevelSetHermiteGeometry] legacy 16-column native operation cannot consume quartic data");
        if (!native) throw std::runtime_error("[LevelSetHermiteGeometry] native-cell query requires a native-leaf snapshot");
        return *native;
    }
    template<class Query> decltype(auto) visitLeaf(Query&& query) const
    {
        if (nativeQuartic) return query(*nativeQuartic);
        return query(leaf());
    }
};

Snapshot::Snapshot(const QuadtreeMesh& mesh, const VertField<double>& phi, bool radialAxisSymmetry,
    CartesianGridScope scope, Reconstruction reconstruction)
{
    // Check before the engine constructor invokes the original Q4 fit.
    Exact::validateEnvironment();
    impl_ = std::make_unique<Impl>(mesh, phi, radialAxisSymmetry, scope, reconstruction);
}

Reconstruction Snapshot::reconstruction() const
{ const auto& backend=checked(); return backend.quartic || backend.nativeQuartic
    ? Reconstruction::EnrichedQuartic : Reconstruction::Cubic; }

Snapshot::Snapshot(const QuadtreeMesh& mesh, const Quantum::Geometry::ActiveHermitePatch& patch,
    bool radialAxisSymmetry, CartesianGridScope scope)
{
    Exact::validateEnvironment();
    impl_ = std::make_unique<Impl>(mesh, patch, radialAxisSymmetry, scope);
}

Snapshot::Snapshot(std::unique_ptr<const Impl> impl) : impl_(std::move(impl)) {}
Snapshot Snapshot::fromNativeLeaves(const QuadtreeMesh& mesh,
    const Quantum::Geometry::ActiveHermitePatch& patch, bool radialAxisSymmetry)
{
    Exact::validateEnvironment();
    return Snapshot(std::make_unique<Impl>(mesh,patch,radialAxisSymmetry));
}
Snapshot Snapshot::fromNativeLeaves(const QuadtreeMesh& mesh,
    const Quantum::Geometry::ActiveQuarticHermitePatch& patch,bool radialAxisSymmetry)
{
    Exact::validateEnvironment();
    return Snapshot(std::make_unique<Impl>(mesh,patch,radialAxisSymmetry));
}
Snapshot Snapshot::refineNativeQuartic(const QuadtreeMesh& targetMesh,const Snapshot& source)
{
    Exact::validateEnvironment();const auto& data=source.checked();
    if(!data.nativeQuartic)throw std::invalid_argument("[LevelSetHermiteGeometry] exact refinement requires a native Q4 source");
    return Snapshot(std::make_unique<Impl>(targetMesh,*data.nativeQuartic));
}
bool Snapshot::usesNativeLeaves() const { return checked().usesNative(); }
Snapshot Snapshot::fromNativeQuarticQueryCells(const QuadtreeMesh& mesh,
    const Quantum::Geometry::ActiveQuarticHermitePatch& patch,const std::vector<int>& cells,bool radialAxisSymmetry)
{
    Exact::validateEnvironment();
    return Snapshot(std::make_unique<Impl>(mesh,patch,cells,radialAxisSymmetry));
}
std::size_t Snapshot::nativeCellCount() const
{ return checked().visitLeaf([](const auto& leaf){return leaf.cellCount();}); }
bool Snapshot::nativeCellActive(int index) const
{ return checked().visitLeaf([&](const auto& leaf){return leaf.cellActive(index);}); }
std::array<Vector2D,2> Snapshot::nativeCellBounds(int index) const
{ return checked().visitLeaf([&](const auto& leaf){return leaf.leafBounds(index);}); }
bool Snapshot::nativeZeroSetPartialResolved(int index,int axis,Vector2D lower,Vector2D upper) const
{
    if(axis<0 || axis>1)throw std::invalid_argument("native zero-set partial requires axis 0 or 1");
    return checked().visitLeaf([&](const auto& leaf){return leaf.resolvedZeroSetPartial(index,axis,lower,upper);});
}
std::array<double,25> Snapshot::nativeQuarticLocalDofs(int index) const
{
    const auto& backend=checked();
    if(!backend.nativeQuartic) throw std::invalid_argument(
        "[LevelSetHermiteGeometry] local Q4 DOFs require a native quartic source");
    return backend.nativeQuartic->mappedCellLocalDofs(index);
}
Quantum::Geometry::PolynomialAffineCoordinate Snapshot::nativeQuarticAffineEnclosure(int index,
    const std::array<Quantum::Geometry::PolynomialAffineCoordinate,2>& physical) const
{
    const auto& backend=checked();
    if(!backend.nativeQuartic)throw std::invalid_argument("[LevelSetHermiteGeometry] affine Q4 enclosure requires a native quartic source");
    return backend.nativeQuartic->affineEnclosure(index,physical);
}
bool Snapshot::identicalNativeQuarticCell(int index,const Snapshot& other,int otherIndex) const
{
    const auto& a=checked();const auto& b=other.checked();
    if(!a.nativeQuartic || !b.nativeQuartic) throw std::invalid_argument(
        "[LevelSetHermiteGeometry] exact native cell identity requires two native quartic sources");
    return a.nativeQuartic->identicalCell(index,*b.nativeQuartic,otherIndex);
}
namespace
{
using Quantum::Geometry::Detail::FaceTime::FaceTimeRational;
using Quantum::Geometry::Detail::FaceTime::faceTimeValue;
}
long double constantNormalFaceTransferVolume(Vector2D a,Vector2D b,double normalSpeed,
    double begin,double end,bool radial)
{
    Exact::validateEnvironment();
    const auto require=[](bool ok,const char* message){if(!ok)throw std::invalid_argument(
        std::string("[LevelSetHermiteGeometry] ")+message);};
    for(double value:{a.x,a.y,b.x,b.y,normalSpeed,begin,end})require(std::isfinite(value),"nonfinite constant face transfer input");
    require(end>begin,"invalid constant face transfer interval");
    const bool horizontal=a.y==b.y && a.x<b.x;
    require(horizontal || (a.x==b.x && a.y<b.y),"constant transfer requires ascending axis-aligned face");
    require(!radial || (a.x>=0 && b.x>=0),"negative constant face radius");
    const auto number=[](double v){return Exact::fromDouble(v);};
    const auto difference=[&](double x,double y){return Exact::subtract(number(x),number(y));};
    using R=FaceTimeRational;
    const R one(1),length(horizontal?difference(b.x,a.x):difference(b.y,a.y),number(1));
    const auto measure=length*(radial?R(a.x)+R(b.x):one);
    const R width(Exact::multiply(number(std::abs(normalSpeed)),difference(end,begin)),number(1));
    const long double volume=faceTimeValue(width*measure)*(radial?std::acos(-1.L):1);
    require(std::isfinite(volume) && volume>=0,"unrepresentable constant face transfer volume");
    return volume;
}
Snapshot::AffineFaceTimeMeasure Snapshot::integrateNativeAffineTranslatedFace(Vector2D a,Vector2D b,
    Vector2D velocity,double referenceTime,double begin,double end) const
{
    Exact::validateEnvironment();
    const auto require=[](bool ok,const char* message){if(!ok)throw std::invalid_argument(
        std::string("[LevelSetHermiteGeometry] ")+message);};
    const auto& backend=checked();
    require(bool(backend.nativeQuartic),"affine face-time integral requires native Q4");
    for(double v:{a.x,a.y,b.x,b.y,velocity.x,velocity.y,referenceTime,begin,end})
        require(std::isfinite(v),"nonfinite affine face-time input");
    require(begin>=referenceTime && end>begin,"invalid affine face-time interval");
    const bool horizontal=a.y==b.y && a.x<b.x;
    require(horizontal || (a.x==b.x && a.y<b.y),"affine face-time requires an ascending axis-aligned segment");
    const bool radial=backend.nativeQuartic->radialAxisSymmetry();
    require(!radial || velocity.x==0,"RZ affine face-time requires axial translation");
    const auto number=[](double v){return Exact::fromDouble(v);};
    const auto difference=[&](double x,double y){return Exact::subtract(number(x),number(y));};
    const auto last=difference(end,referenceTime),first=difference(begin,referenceTime),duration=difference(end,begin);
    Vector2D lower{std::min(a.x,b.x),std::min(a.y,b.y)},upper{std::max(a.x,b.x),std::max(a.y,b.y)};
    // Convex hull of every trajectory, from the old authority to interval end.
    // Outward enclosures ensure a rounded departure never hides a touched leaf.
    for(const auto p:{a,b}) {
        const auto x=Exact::enclosure(Exact::subtract(number(p.x),Exact::multiply(number(velocity.x),last)));
        const auto y=Exact::enclosure(Exact::subtract(number(p.y),Exact::multiply(number(velocity.y),last)));
        lower.x=std::min(lower.x,x.lo);lower.y=std::min(lower.y,y.lo);
        upper.x=std::max(upper.x,x.hi);upper.y=std::max(upper.y,y.hi);
    }
    const auto affine=backend.nativeQuartic->affineOnBox(lower,upper,a);
    const auto& n=affine.numerator;
    const auto motion=Exact::add(Exact::multiply(n[1],number(velocity.x)),Exact::multiply(n[2],number(velocity.y)));
    const auto f0=Exact::subtract(n[0],Exact::multiply(motion,first));
    const auto fs=Exact::add(Exact::multiply(n[1],difference(b.x,a.x)),Exact::multiply(n[2],difference(b.y,a.y)));
    const auto ft=Exact::subtract(Exact::Expansion{},Exact::multiply(motion,duration));
    require(Exact::sign(f0)!=0 || Exact::sign(fs)!=0 || Exact::sign(ft)!=0,
        "identically zero face-time trace has positive measure");
    using R=FaceTimeRational;
    const auto moments=Quantum::Geometry::Detail::FaceTime::negativeHalfPlane(f0,fs,ft);
    const auto area=moments.area,moment=moments.first;
    const auto fullWeight=radial?(R(a.x)+R(b.x))/R(2):R(1);
    const auto liquidWeight=radial?R(a.x)*area+(R(b.x)-R(a.x))*moment:area;
    require(Exact::sign(liquidWeight.n)>=0 && Exact::sign((fullWeight-liquidWeight).n)>=0,
        "invalid exact face-time subset measure");
    if(Exact::sign(fullWeight.n)==0)return {};
    const auto fraction=liquidWeight/fullWeight;
    // Evaluate near-one subsets through their exact complement. This is the
    // same rational measure, not clipping or a roundoff inventory correction.
    const auto complement=R(1)-fraction;
    const long double fractionValue=Exact::sign((fraction-R(.5)).n)>0
        ?1-faceTimeValue(complement):faceTimeValue(fraction);
    const R length(horizontal?difference(b.x,a.x):difference(b.y,a.y),number(1));
    const auto total=faceTimeValue(length*R(duration,number(1))*fullWeight)*(radial?2*std::acos(-1.L):1);
    require(std::isfinite(total) && total>0 && fractionValue>=0 && fractionValue<=1,
        "unrepresentable affine face-time result");
    return {total,fractionValue};
}
Snapshot::NormalTranslationSweep Snapshot::nativeNormalTranslationSweep(Vector2D a,Vector2D b,
    double normalSpeed,double referenceTime,double begin,double end,bool radial) const
{
    Exact::validateEnvironment();
    const auto require=[](bool ok,const char* message){if(!ok)throw std::invalid_argument(
        std::string("[LevelSetHermiteGeometry] ")+message);};
    const auto& backend=checked();require(bool(backend.nativeQuartic),"normal sweep requires native Q4");
    for(double value:{a.x,a.y,b.x,b.y,normalSpeed,referenceTime,begin,end})
        require(std::isfinite(value),"nonfinite normal sweep input");
    require(normalSpeed!=0 && begin>=referenceTime && end>begin,"invalid normal sweep speed/interval");
    const bool horizontal=a.y==b.y && a.x<b.x;
    require(horizontal || (a.x==b.x && a.y<b.y),"normal sweep requires ascending axis-aligned face");
    require(!radial || (horizontal && a.x>=0),"RZ normal sweep requires axial translation");
    for(std::size_t cell=0;cell<backend.nativeQuartic->cellCount();++cell)
        require(backend.nativeQuartic->cellActive(static_cast<int>(cell)),"normal sweep requires full active scalar support");
    const auto number=[](double v){return Exact::fromDouble(v);};
    const auto difference=[&](double x,double y){return Exact::subtract(number(x),number(y));};
    const auto normal=number(horizontal?a.y:a.x),speed=number(normalSpeed);
    const auto first=Exact::subtract(normal,Exact::multiply(speed,difference(begin,referenceTime)));
    const auto last=Exact::subtract(normal,Exact::multiply(speed,difference(end,referenceTime)));
    const auto low=Exact::enclosure(normalSpeed>0?last:first),high=Exact::enclosure(normalSpeed>0?first:last);
    NormalTranslationSweep result;
    if(horizontal){result.innerLower={a.x,low.hi};result.innerUpper={b.x,high.lo};
        result.outerLower={a.x,low.lo};result.outerUpper={b.x,high.hi};}
    else{result.innerLower={low.hi,a.y};result.innerUpper={high.lo,b.y};
        result.outerLower={low.lo,a.y};result.outerUpper={high.hi,b.y};}
    require(low.hi<high.lo,"normal sweep inner rectangle is not distinguishable");
    const auto lower=backend.nativeQuartic->minCoord(),upper=backend.nativeQuartic->maxCoord();
    for(const auto p:{a,b,result.outerLower,result.outerUpper})
        require(p.x>=lower.x && p.x<=upper.x && p.y>=lower.y && p.y<=upper.y,
            "normal sweep leaves original scalar domain; exterior history required");
    using R=FaceTimeRational;
    const R one(1),length(horizontal?difference(b.x,a.x):difference(b.y,a.y),number(1));
    const auto measure=length*(radial?R(a.x)+R(b.x):one);
    const R width(Exact::multiply(number(std::abs(normalSpeed)),difference(end,begin)),number(1));
    const R outerWidth(difference(high.hi,low.lo),number(1)),innerWidth(difference(high.lo,low.hi),number(1));
    const auto band=outerWidth-innerWidth;
    require(Exact::sign((width-innerWidth).n)>=0 && Exact::sign((outerWidth-width).n)>=0,
        "invalid exact normal sweep containment");
    const long double metric=radial?std::acos(-1.L):1;
    result.totalVolume=constantNormalFaceTransferVolume(a,b,normalSpeed,begin,end,radial);
    result.coordinateVolumeBand=faceTimeValue(band*measure)*metric;
    require(std::isfinite(result.totalVolume) && result.totalVolume>0 &&
        std::isfinite(result.coordinateVolumeBand) && result.coordinateVolumeBand>=0,
        "unrepresentable normal sweep volume");
    return result;
}
bool Snapshot::identicalNativeQuarticRepresentation(int index,const Snapshot& other,int otherIndex) const
{
    const auto& a=checked(); const auto& b=other.checked();
    if(!a.nativeQuartic || !b.nativeQuartic) return false;
    return a.nativeQuartic->identicalCellRepresentation(index,*b.nativeQuartic,otherIndex);
}
std::vector<int> Snapshot::nativeCellCandidates(Vector2D lower, Vector2D upper) const
{ return checked().visitLeaf([&](const auto& leaf){return leaf.leafCandidates(lower,upper);}); }
std::vector<Quantum::Geometry::ActiveHermiteNodalSensitivity> Snapshot::nativeCellSensitivity(
    int index, const std::array<long double,16>& coefficients,
    Quantum::Geometry::CornerSectorVariation sectors) const
{ return checked().leaf().pullbackCellSensitivity(index, coefficients, sectors); }
std::vector<Quantum::Geometry::ActiveQuarticHermiteSensitivity> Snapshot::nativeQuarticCellSensitivity(
    int index,const std::array<long double,25>& coefficients,Quantum::Geometry::CornerSectorVariation sectors) const
{
    const auto& backend=checked();
    if(!backend.nativeQuartic) throw std::invalid_argument(
        "[LevelSetHermiteGeometry] 25-column native sensitivity requires a native quartic source");
    return backend.nativeQuartic->pullbackCellSensitivity(index,coefficients,sectors);
}
Snapshot::ScalarJet2 Snapshot::scalarJetInNativeCell(int index, Vector2D point) const
{
    return checked().visitLeaf([&](const auto& native)->ScalarJet2 {
        Exact::validateEnvironment();
        auto sample = native.evaluateCell(index,point);
        sample.value = native.valueInCell(index,point);
        return {sample.value,sample.dx,sample.dy,sample.dxx,sample.dxy,sample.dyy};
    });
}
Snapshot::ScalarJet2 Snapshot::scalarJetAtNativeCellChart(int index,
    const Quantum::Geometry::CellChartPoint& point) const
{
    const auto sample=checked().visitLeaf([&](const auto& native){return native.evaluateCellChart(index,point);});
    return {sample.value,sample.dx,sample.dy,sample.dxx,sample.dxy,sample.dyy};
}
Frame Snapshot::jetAtNativeCellChart(int index,const Quantum::Geometry::CellChartPoint& point) const
{
    return checked().visitLeaf([&](const auto& native) {
        const auto sample=native.evaluateCellChart(index,point);
        const long double gradient=std::hypot(static_cast<long double>(sample.dx),static_cast<long double>(sample.dy));
        const auto bounds=native.leafBounds(index);
        Engine::Detail::require(Quantum::Geometry::Detail::HermiteDifferential::usableGradient(
            sample,gradient,native.derivativeResolutionScale(index),bounds[1]-bounds[0]),
            "unresolved native cell-chart interface gradient");
        auto result=Quantum::Geometry::Detail::HermiteDifferential::frame<Frame>(
            native.chartDisplayPoint(index,point),sample,gradient,Engine::Detail::checked);
        result.nativeCellIndex=index; result.cellChartPoint=point; result.hasCellChartPoint=true;
        return result;
    });
}
double Snapshot::curvatureTangentialDerivativeAtNativeCellChart(int index,
    const Quantum::Geometry::CellChartPoint& point) const
{
    return checked().visitLeaf([&](const auto& native) {
        const auto sample=native.evaluateCellChart(index,point);
        const long double gradient=std::hypot(static_cast<long double>(sample.dx),static_cast<long double>(sample.dy));
        const auto bounds=native.leafBounds(index);
        Engine::Detail::require(Quantum::Geometry::Detail::HermiteDifferential::usableGradient(
            sample,gradient,native.derivativeResolutionScale(index),bounds[1]-bounds[0]),
            "unresolved native cell-chart curvature derivative");
        return Quantum::Geometry::Detail::HermiteDifferential::curvatureDerivative(
            sample,native.thirdDerivativesCellChart(index,point),gradient,Engine::Detail::checked);
    });
}

Snapshot::~Snapshot() = default;
Snapshot::Snapshot(Snapshot&&) noexcept = default;
Snapshot& Snapshot::operator=(Snapshot&&) noexcept = default;

const Snapshot::Impl& Snapshot::checked() const
{
    if (!impl_) throw std::runtime_error("[LevelSetHermiteGeometry] query on a moved-from snapshot");
    return *impl_;
}

int Snapshot::nx() const { return checked().visitGeometry([&](const auto& geometry) { return geometry.nx(); }); }
int Snapshot::ny() const { return checked().visitGeometry([&](const auto& geometry) { return geometry.ny(); }); }
bool Snapshot::radialAxisSymmetry() const
{ const auto& impl = checked(); return impl.usesNative() ? impl.visitLeaf([](const auto& leaf){return leaf.radialAxisSymmetry();}) : impl.visitGeometry([&](const auto& geometry) { return geometry.radialAxisSymmetry(); }); }
bool Snapshot::restrictedDomain() const
{ const auto& impl = checked(); return impl.usesNative() ? true : impl.visitGeometry([&](const auto& geometry) { return geometry.restrictedDomain(); }); }
bool Snapshot::hasCornerSectors() const
{ const auto& impl = checked(); return impl.usesNative() ? impl.visitLeaf([](const auto& leaf){return leaf.hasCornerSectors();}) : impl.visitGeometry([&](const auto& geometry) { return geometry.hasCornerSectors(); }); }
bool Snapshot::cellActive(int i, int j) const { return checked().visitGeometry([&](const auto& geometry) { return geometry.cellActive(i, j); }); }
Vector2D Snapshot::minCoord() const
{ const auto& impl = checked(); return impl.usesNative() ? impl.visitLeaf([](const auto& leaf){return leaf.minCoord();}) : impl.visitGeometry([&](const auto& geometry) { return geometry.minCoord(); }); }
Vector2D Snapshot::maxCoord() const
{ const auto& impl = checked(); return impl.usesNative() ? impl.visitLeaf([](const auto& leaf){return leaf.maxCoord();}) : impl.visitGeometry([&](const auto& geometry) { return geometry.maxCoord(); }); }
Vector2D Snapshot::cellSize() const { return checked().visitGeometry([&](const auto& geometry) { return geometry.cellSize(); }); }
std::size_t Snapshot::candidateCellCount() const { return checked().visitGeometry([&](const auto& geometry) { return geometry.candidateCellCount(); }); }
std::size_t Snapshot::seedCount() const { return checked().visitGeometry([&](const auto& geometry) { return geometry.seedCount(); }); }
std::size_t Snapshot::vertexCount() const
{ const auto& impl = checked(); return impl.usesNative() ? impl.visitLeaf([](const auto& leaf){return leaf.vertexCount();}) : impl.visitGeometry([&](const auto& geometry) { return geometry.vertexCount(); }); }

std::array<int, 2> Snapshot::zeroVertexGradientSignsInCell(int i, int j, Vector2D point) const
{
    return checked().visitGeometry([&](const auto& geometry) -> std::array<int, 2>
    {
        Exact::validateEnvironment();
        return geometry.zeroVertexGradientSignsInCell(i, j, point);
    });
}

std::vector<Vector2D> Snapshot::intersectCellSegment(int i, int j, Vector2D a, Vector2D b) const
{
    return checked().visitGeometry([&](const auto& geometry) -> std::vector<Vector2D>
    {
        Exact::validateEnvironment();
        return geometry.intersectCellSegment(i, j, a, b);
    });
}

Frame Snapshot::jetAtCellLocal(int i, int j, Vector2D local) const
{
    return checked().visitGeometry([&](const auto& geometry) -> Frame
    {
        Exact::validateEnvironment();
        return resultFrame(geometry.jetAtCellLocal(i, j, local));
    });
}

Frame Snapshot::jetAtCellChart(int i, int j, const Quantum::Geometry::CellChartPoint& point) const
{
    return checked().visitGeometry([&](const auto& geometry) -> Frame
    {
        Exact::validateEnvironment();
        return resultFrame(geometry.jetAtCellChart(i, j, point));
    });
}

std::vector<Quantum::Geometry::CellChartPoint> Snapshot::intersectCellSection(int i, int j,
    bool coordinateIsX, Quantum::Geometry::CellChartCoordinate fixedCoordinate) const
{
    return checked().visitGeometry([&](const auto& geometry) -> std::vector<Quantum::Geometry::CellChartPoint>
    {
        Exact::validateEnvironment();
        return geometry.intersectCellSection(i, j, coordinateIsX, fixedCoordinate);
    });
}

Snapshot::ScalarJet2 Snapshot::scalarJetAtCellLocal(int i, int j, Vector2D local) const
{
    return checked().visitGeometry([&](const auto& geometry) -> Snapshot::ScalarJet2
    {
        Exact::validateEnvironment();
        const auto sample = geometry.scalarJetAtCellLocal(i, j, local);
        return {sample.value, sample.dx, sample.dy, sample.dxx, sample.dxy, sample.dyy};
    });
}

Snapshot::ScalarJet2 Snapshot::scalarJetAtCellChart(int i, int j,
    const Quantum::Geometry::CellChartPoint& point) const
{
    return checked().visitGeometry([&](const auto& geometry) -> Snapshot::ScalarJet2
    {
        Exact::validateEnvironment();
        const auto sample = geometry.scalarJetAtCellChart(i, j, point);
        return {sample.value, sample.dx, sample.dy, sample.dxx, sample.dxy, sample.dyy};
    });
}

double Snapshot::curvatureTangentialDerivativeAtCellChart(int i, int j,
    const Quantum::Geometry::CellChartPoint& point) const
{
    return checked().visitGeometry([&](const auto& geometry) -> double
    {
        Exact::validateEnvironment();
        return geometry.curvatureTangentialDerivativeAtCellChart(i, j, point);
    });
}

std::vector<Quantum::Geometry::ActiveQuarticHermiteSensitivity> Snapshot::nativeQuarticDifferentialSensitivity(
    int index,const Quantum::Geometry::CellChartPoint& point,const std::array<long double,6>& weights) const
{
    Exact::validateEnvironment();
    Quantum::Geometry::validateCellChart(point);
    const auto& data=checked();
    if (!data.nativeQuartic) throw std::invalid_argument("Q4 differential sensitivity requires native quartic data");
    const auto& scalar=*data.nativeQuartic;
    const auto box=scalar.leafBounds(index);
    const long double hx=static_cast<long double>(box[1].x)-box[0].x,hy=static_cast<long double>(box[1].y)-box[0].y;
    const auto basis=[](Quantum::Geometry::CellChartCoordinate chart) {
        const auto one=Exact::fromDouble(1),offset=Exact::fromDouble(chart.offset);
        const auto values=Quantum::Geometry::Detail::LevelSetHermiteScalar::QuarticCell::rationalBasisJet(
            chart.fromUpper?Exact::subtract(one,offset):offset,one,true);
        std::array<std::array<long double,5>,3> result{};
        for (int k=0;k<5;++k) {
            result[0][k]=Exact::estimate(values.value[k])/Exact::estimate(values.denominator);
            result[1][k]=Exact::estimate(values.derivative[k])/Exact::estimate(values.denominator);
            result[2][k]=Exact::estimate(values.second[k])/Exact::estimate(values.denominator);
        }
        return result;
    };
    for (auto weight:weights) if (!std::isfinite(weight)) throw std::invalid_argument("nonfinite differential sensitivity weight");
    const auto bx=basis(point.x),by=basis(point.y);
    std::array<long double,25> coefficients{};
    for (int y=0;y<5;++y) for (int x=0;x<5;++x)
        coefficients[5*y+x]=(x==2 || x==3?hx:1.L)*(y==2 || y==3?hy:1.L)*
            (weights[0]*bx[0][x]*by[0][y]+weights[1]*bx[1][x]*by[0][y]/hx
                +weights[2]*bx[0][x]*by[1][y]/hy+weights[3]*bx[2][x]*by[0][y]/(hx*hx)
                +weights[4]*bx[1][x]*by[1][y]/(hx*hy)+weights[5]*bx[0][x]*by[2][y]/(hy*hy));
    return scalar.pullbackCellSensitivity(index,coefficients,Quantum::Geometry::CornerSectorVariation::RequireAbsent);
}
Snapshot::ScalarJet2 Snapshot::scalarJetNative(Vector2D point) const
{
    return checked().visitLeaf([&](const auto& native)->ScalarJet2 {
        Exact::validateEnvironment();
        const int cell=native.selectedCell(point);
        const auto sample=native.evaluateCell(cell,point);
        return {native.valueInCell(cell,point),sample.dx,sample.dy,sample.dxx,sample.dxy,sample.dyy};
    });
}

double Snapshot::value(Vector2D point) const
{
    const auto& impl = checked();
    Exact::validateEnvironment();
    if (impl.usesNative())
        return impl.visitLeaf([&](const auto& leaf){return leaf.valueInCell(leaf.selectedCell(point),point);});
    return impl.visitGeometry([&](const auto& geometry) { return geometry.value(point); });
}

bool Snapshot::negative(Vector2D point) const
{
    const auto& impl = checked();
    Exact::validateEnvironment();
    if (impl.usesNative()) return impl.visitLeaf([&](const auto& leaf)
    {
        const auto p = leaf.restrictPhysical(leaf.selectedCell(point),point,point);
        return Exact::sign(p.numerator[0]) <= 0;
    });
    return impl.visitGeometry([&](const auto& geometry) { return geometry.negative(point); });
}

int Snapshot::closedCellSignExceptVertices(Vector2D point, unsigned allowedCornerMask) const
{
    const auto& impl = checked();
    Exact::validateEnvironment();
    if (impl.usesNative()) return impl.visitLeaf([&](const auto& leaf)
        {return leaf.cellSignExceptVertices(leaf.selectedCell(point),allowedCornerMask);});
    return impl.visitGeometry([&](const auto& geometry) { return geometry.closedCellSignExceptVertices(point, allowedCornerMask); });
}

int Snapshot::closedCellSign(Vector2D point) const
{
    const auto& impl = checked();
    Exact::validateEnvironment();
    if (impl.usesNative()) return impl.visitLeaf([&](const auto& leaf)
        {return leaf.cellSignExceptVertices(leaf.selectedCell(point),0u);});
    return impl.visitGeometry([&](const auto& geometry) { return geometry.closedCellSign(point); });
}

std::array<int, 2> Snapshot::zeroVertexGradientSigns(Vector2D point) const
{
    return checked().visitGeometry([&](const auto& geometry) -> std::array<int, 2>
    {
        Exact::validateEnvironment();
        return geometry.zeroVertexGradientSigns(point);
    });
}

Frame Snapshot::jet(Vector2D point) const
{
    return checked().visitGeometry([&](const auto& geometry) -> Frame
    {
        Exact::validateEnvironment();
        return resultFrame(geometry.jet(point));
    });
}

Frame Snapshot::jetInCell(int i, int j, Vector2D point) const
{
    return checked().visitGeometry([&](const auto& geometry) -> Frame
    {
        Exact::validateEnvironment();
        return resultFrame(geometry.jetInCell(i, j, point));
    });
}

bool Snapshot::closedCellOutsideBand(Vector2D point, double halfWidth) const
{
    const auto& impl = checked();
    Exact::validateEnvironment();
    if (impl.usesNative()) return impl.visitLeaf([&](const auto& leaf)
        {return leaf.closedCellOutsideBand(leaf.selectedCell(point),halfWidth);});
    return impl.visitGeometry([&](const auto& geometry) { return geometry.closedCellOutsideBand(point, halfWidth); });
}

bool Snapshot::closedCellInsideBand(Vector2D point, double halfWidth) const
{
    const auto& impl = checked();
    Exact::validateEnvironment();
    if (impl.usesNative()) return impl.visitLeaf([&](const auto& leaf)
        {return leaf.closedCellInsideBand(leaf.selectedCell(point),halfWidth);});
    return impl.visitGeometry([&](const auto& geometry) { return geometry.closedCellInsideBand(point, halfWidth); });
}

std::size_t Snapshot::firstNativeBandViolation(const std::vector<int>& cells,
    double halfWidth, BandRequirement requirement) const
{
    const auto& impl=checked();
    Exact::validateEnvironment();
    if (!impl.usesNative() || !std::isfinite(halfWidth) || halfWidth<0
        || (requirement!=BandRequirement::EntireCellInside && requirement!=BandRequirement::EntireCellOutside))
        throw std::invalid_argument("native band batch requires native leaves, finite nonnegative width and a valid requirement");
    return impl.visitLeaf([&](const auto& leaf) {
        std::vector<unsigned char> satisfied(cells.size());
        std::vector<std::exception_ptr> errors(cells.size());
        const auto evaluate=[&](std::size_t i) {
            satisfied[i]=requirement==BandRequirement::EntireCellInside
                ? leaf.closedCellInsideBand(cells[i],halfWidth) : leaf.closedCellOutsideBand(cells[i],halfWidth);
        };
        if (!cells.empty())
        {
            evaluate(0);
            if (!satisfied[0]) return std::size_t{0};
        }
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 16) if(cells.size()>32)
#endif
        for (std::ptrdiff_t i=1;i<static_cast<std::ptrdiff_t>(cells.size());++i)
        {
            try {evaluate(static_cast<std::size_t>(i));}
            catch (...) {errors[i]=std::current_exception();}
        }
        // Match the original serial predicate/exception order, even when a
        // speculative query beyond the first violation rejected its input.
        for (std::size_t i=1;i<cells.size();++i)
        {
            if (errors[i]) std::rethrow_exception(errors[i]);
            if (!satisfied[i]) return i;
        }
        return cells.size();
    });
}

Frame Snapshot::projectedFrame(Vector2D point) const
{
    return checked().visitGeometry([&](const auto& geometry) -> Frame
    {
        Exact::validateEnvironment();
        return resultFrame(geometry.projectedFrame(point));
    });
}

Frame Snapshot::projectedCellLocalFrame(Vector2D point) const
{
    return checked().visitGeometry([&](const auto& geometry) -> Frame
    {
        Exact::validateEnvironment();
        return resultFrame(geometry.projectedCellLocalFrame(point));
    });
}

Snapshot::LocalStationaryCandidates Snapshot::stationaryCellLocalFrames(Vector2D point) const
{
    return checked().visitGeometry([&](const auto& geometry) -> Snapshot::LocalStationaryCandidates
    {
        Exact::validateEnvironment();
        Engine::Frame diagnostics;
        const auto frames = geometry.stationaryCellLocalFrames(point, diagnostics);
        LocalStationaryCandidates result;
        result.attemptedSeeds = diagnostics.attemptedSeeds;
        result.convergedSeeds = diagnostics.convergedSeeds;
        result.frames.reserve(frames.size());
        for (const auto& frame : frames) result.frames.push_back(resultFrame(frame));
        return result;
    });
}

BoundaryProjection Snapshot::boundaryProjectedFrame(Vector2D query,
    const std::vector<Vector2D>& boundaryEndpoints) const
{
    return checked().visitGeometry([&](const auto& geometry) -> BoundaryProjection
    {
        Exact::validateEnvironment();
        const auto result = geometry.boundaryProjectedFrame(query, boundaryEndpoints);
        return {resultFrame(result.frame), result.endpoint, result.inwardDerivative, result.kktResidual};
    });
}

BoundaryProjection Snapshot::patchBoundaryProjectedFrame(Vector2D query,
    const std::vector<PatchEndpoint>& boundaryEndpoints) const
{
    return checked().visitGeometry([&](const auto& geometry) -> BoundaryProjection
    {
        Exact::validateEnvironment();
        std::vector<Engine::PatchEndpoint> endpoints;
        endpoints.reserve(boundaryEndpoints.size());
        for (const auto& endpoint : boundaryEndpoints)
            endpoints.push_back({endpoint.point, endpoint.cellI, endpoint.cellJ, endpoint.inwardTangentSign});
        const auto result = geometry.patchBoundaryProjectedFrame(query, endpoints);
        return {resultFrame(result.frame), result.endpoint, result.inwardDerivative, result.kktResidual};
    });
}

double Snapshot::curvatureTangentialDerivative(Vector2D point) const
{
    return checked().visitGeometry([&](const auto& geometry) -> double
    {
        Exact::validateEnvironment();
        return geometry.curvatureTangentialDerivative(point);
    });
}

double Snapshot::curvatureTangentialDerivativeAtCellLocal(int i, int j, Vector2D local) const
{
    return checked().visitGeometry([&](const auto& geometry) -> double
    {
        Exact::validateEnvironment();
        return geometry.curvatureTangentialDerivativeAtCellLocal(i, j, local);
    });
}

double Snapshot::reconstructedCurvatureTangentialDerivative(Vector2D point) const
{
    return checked().visitGeometry([&](const auto& geometry) -> double
    {
        Exact::validateEnvironment();
        return geometry.reconstructedCurvatureTangentialDerivative(point);
    });
}

bool Snapshot::segmentInsideActiveDomain(Vector2D a, Vector2D b) const
{
    Exact::validateEnvironment();
    if (!checked().usesNative()) return checked().visitGeometry([&](const auto& geometry) { return geometry.segmentInsideActiveDomain(a, b); });
    namespace Mapped = Quantum::Geometry::Detail::MappedHermiteSegment;
    return checked().visitLeaf([&](const auto& source) {
        const auto lo = source.minCoord(), hi = source.maxCoord();
        for (auto p : {a, b})
            Mapped::require(std::isfinite(p.x) && std::isfinite(p.y) && p.x >= lo.x && p.x <= hi.x
                && p.y >= lo.y && p.y <= hi.y, "activity segment must lie in the finite native rectangle");
        Mapped::BasicNativeCellLayout<std::decay_t<decltype(source)>> layout(source);
        auto cuts = layout.cuts(a, b);
        if (layout.incidentCellsAtInputPoint(a).empty() || layout.incidentCellsAtInputPoint(b).empty()) return false;
        std::sort(cuts.begin(), cuts.end(), [](const auto& x, const auto& y) { return Mapped::compare(x,y) < 0; });
        // Exact clip intervals, not phase samples or a rounded global midpoint.
        // A solid gap of any resolved width rejects; touching a fluid wall is legal.
        for (std::size_t i = 1; i < cuts.size(); ++i)
            if (Mapped::compare(cuts[i-1], cuts[i]) < 0
                && layout.incidentCells(Mapped::between(cuts[i-1], cuts[i], .5)).empty()) return false;
        return true;
    });
}

int Snapshot::certifiedSegmentSign(Vector2D a, Vector2D b) const
{
    return checked().visitGeometry([&](const auto& geometry) -> int
    {
        Exact::validateEnvironment();
        return geometry.certifiedSegmentSign(a, b);
    });
}

std::vector<Piece> Snapshot::partitionSegment(Vector2D a, Vector2D b) const
{
    return checked().visitGeometry([&](const auto& geometry) -> std::vector<Piece>
    {
        Exact::validateEnvironment();
        const auto pieces = geometry.partitionSegment(a, b);
        std::vector<Piece> result;
        result.reserve(pieces.size());
        for (const auto& piece : pieces) result.push_back({piece.begin, piece.end, piece.negative});
        return result;
    });
}

std::vector<ParametricPiece> Snapshot::partitionSegmentParametric(Vector2D a, Vector2D b, SegmentPartitionPolicy policy) const
{
    return checked().visitGeometry([&](const auto& geometry) -> std::vector<ParametricPiece>
    {
        Exact::validateEnvironment();
        if (policy != SegmentPartitionPolicy::RequireRepresentableGridPieces
            && policy != SegmentPartitionPolicy::CertifyUnresolvedGridPieces)
            throw std::invalid_argument("[LevelSetHermiteGeometry] invalid segment partition policy");
        const auto pieces = geometry.partitionSegmentParametric(a, b,
            policy == SegmentPartitionPolicy::CertifyUnresolvedGridPieces);
        std::vector<ParametricPiece> result;
        result.reserve(pieces.size());
        for (const auto& piece : pieces)
            result.push_back({piece.anchorBegin, piece.anchorEnd, piece.tBegin, piece.tEnd, piece.negative});
        return result;
    });
}

std::vector<SegmentRootIncidence> Snapshot::intersectSegmentWithProvenance(
    Vector2D a, Vector2D b, SegmentPartitionPolicy policy) const
{
    return checked().visitGeometry([&](const auto& geometry) -> std::vector<SegmentRootIncidence>
    {
        Exact::validateEnvironment();
        if (policy != SegmentPartitionPolicy::RequireRepresentableGridPieces
            && policy != SegmentPartitionPolicy::CertifyUnresolvedGridPieces)
            throw std::invalid_argument("[LevelSetHermiteGeometry] invalid segment intersection policy");
        const auto incidences = geometry.intersectSegmentWithProvenance(a, b,
            policy == SegmentPartitionPolicy::CertifyUnresolvedGridPieces);
        std::vector<SegmentRootIncidence> result;
        result.reserve(incidences.size());
        for (const auto& root : incidences)
            result.push_back({root.cellI, root.cellJ, root.anchorBegin, root.anchorEnd,
                root.localBegin, root.localEnd, root.localPoint, root.parameter, root.pointOnly});
        return result;
    });
}

GridEdgeZeros Snapshot::intersectGridEdge(Vector2D a, Vector2D b) const
{
    return checked().visitGeometry([&](const auto& geometry) -> GridEdgeZeros
    {
        Exact::validateEnvironment();
        auto result = geometry.intersectGridEdge(a, b);
        return {std::move(result.points), result.coincident, std::move(result.parameters)};
    });
}

double Snapshot::negativeFraction(Vector2D a, Vector2D b) const
{
    return checked().visitGeometry([&](const auto& geometry) -> double
    {
        Exact::validateEnvironment();
        return geometry.negativeFraction(a, b);
    });
}
}

// Keep the mapped query beside the opaque snapshot in this strict-FP TU:
// its implementation-only expansion/engine types must not cross a TU ABI.
struct LevelSetMappedSegment::Impl
{
    static const NativeQuarticTrace::Impl& retainedTrace(const NativeQuarticTrace& trace);
    using UniformQuery = Quantum::Geometry::Detail::MappedHermiteSegment::Query;
    using NativeQuery = Quantum::Geometry::Detail::MappedHermiteSegment::NativeQuery;
    using QuarticQuery = Quantum::Geometry::Detail::MappedHermiteSegment::QuarticQuery;
    using NativeQuarticQuery = Quantum::Geometry::Detail::MappedHermiteSegment::NativeQuarticQuery;
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source;
    std::variant<UniformQuery,NativeQuery,QuarticQuery,NativeQuarticQuery> query;
    Impl(std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> retained,
        const Quantum::Geometry::Detail::LevelSetHermiteGeometry::Snapshot& geometry,
        Vector2D begin, Vector2D end)
        : source(std::move(retained)), query(std::in_place_type<UniformQuery>,geometry,begin,end) {}
    Impl(std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> retained,
        const Quantum::Geometry::Detail::LevelSetHermiteScalar::LeafSnapshot& scalar,
        Vector2D begin, Vector2D end)
        : source(std::move(retained)), query(std::in_place_type<NativeQuery>,scalar,begin,end) {}
    Impl(std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> retained,
        const Quantum::Geometry::Detail::MappedHermiteSegment::QuarticGeometry& geometry,
        Vector2D begin, Vector2D end)
        : source(std::move(retained)), query(std::in_place_type<QuarticQuery>,geometry,begin,end) {}
    Impl(std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> retained,
        const Quantum::Geometry::Detail::LevelSetHermiteScalar::QuarticLeafData& scalar,
        Vector2D begin,Vector2D end)
        : source(std::move(retained)),query(std::in_place_type<NativeQuarticQuery>,scalar,begin,end) {}
    Impl(std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> retained,NativeQuarticQuery nativeQuery)
        : source(std::move(retained)),query(std::in_place_type<NativeQuarticQuery>,std::move(nativeQuery)) {}
    const UniformQuery& uniformQuery() const
    {
        const auto* result = std::get_if<UniformQuery>(&query);
        if (!result) throw std::runtime_error("[LevelSetMappedInterface] uniform cut access requires a uniform query");
        return *result;
    }
};

LevelSetMappedSegment::LevelSetMappedSegment(
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source, Vector2D begin, Vector2D end)
{
    if (!source) throw std::invalid_argument("[LevelSetMappedSegment] null retained source");
    Quantum::Geometry::Detail::ExactExpansion::validateEnvironment();
    const auto& backend = source->checked();
    if (backend.native) impl_ = std::make_unique<Impl>(std::move(source),*backend.native,begin,end);
    else if (backend.nativeQuartic) impl_ = std::make_unique<Impl>(std::move(source),*backend.nativeQuartic,begin,end);
    else if (backend.quartic) impl_ = std::make_unique<Impl>(std::move(source),*backend.quartic,begin,end);
    else impl_ = std::make_unique<Impl>(std::move(source),backend.cubicGeometry(),begin,end);
}
LevelSetMappedSegment::~LevelSetMappedSegment() = default;
LevelSetMappedSegment::LevelSetMappedSegment(std::unique_ptr<const Impl> impl) : impl_(std::move(impl)) {}
LevelSetMappedSegment::LevelSetMappedSegment(LevelSetMappedSegment&&) noexcept = default;
LevelSetMappedSegment& LevelSetMappedSegment::operator=(LevelSetMappedSegment&&) noexcept = default;
const LevelSetMappedSegment::Impl& LevelSetMappedSegment::checked() const
{
    if (!impl_) throw std::runtime_error("[LevelSetMappedSegment] query on a moved-from segment");
    Quantum::Geometry::Detail::ExactExpansion::validateEnvironment();
    return *impl_;
}
bool LevelSetMappedSegment::usesSource(const LevelSetHermiteGeometry::Snapshot& source) const
{ return checked().source.get() == &source; }
Vector2D LevelSetMappedSegment::begin() const { return std::visit([](const auto& query) { return query.inputBegin; },checked().query); }
Vector2D LevelSetMappedSegment::end() const { return std::visit([](const auto& query) { return query.inputEnd; },checked().query); }
const std::vector<LevelSetMappedSegment::Range>& LevelSetMappedSegment::ranges() const
{ return std::visit([](const auto& query) -> const auto& { return query.ranges; },checked().query); }
const std::vector<LevelSetMappedSegment::Root>& LevelSetMappedSegment::roots() const
{ return std::visit([](const auto& query) -> const auto& { return query.rootInfo; },checked().query); }
bool LevelSetMappedSegment::rootAtInputEndpoint(std::size_t rootIndex, bool atEnd) const
{
    namespace Mapped = Quantum::Geometry::Detail::MappedHermiteSegment;
    return std::visit([&](const auto& query) {
        Mapped::require(rootIndex < query.rootInfo.size(), "root endpoint query index is out of range");
        const bool reverse = std::make_pair(query.inputEnd.x,query.inputEnd.y)
            < std::make_pair(query.inputBegin.x,query.inputBegin.y);
        const auto& root = query.roots[reverse ? query.roots.size()-1-rootIndex : rootIndex];
        return Mapped::compare(root.parameter,Mapped::constant(atEnd != reverse ? 1.0 : 0.0)) == 0;
    },checked().query);
}
std::size_t LevelSetMappedSegment::gridPieceCount() const { return std::visit([](const auto& query) { return query.gridPieces; },checked().query); }
bool LevelSetMappedSegment::sameNativeQuarticTraceZeros(const LevelSetMappedSegment& other) const
{
    const auto* a=std::get_if<Impl::NativeQuarticQuery>(&checked().query);
    const auto* b=std::get_if<Impl::NativeQuarticQuery>(&other.checked().query);
    if (!a || !b || begin().x!=other.begin().x || begin().y!=other.begin().y
        || end().x!=other.end().x || end().y!=other.end().y)
        throw std::invalid_argument("native quartic trace comparison requires the same oriented segment");
    Quantum::Geometry::Detail::ExactExpansion::validateEnvironment();
    return a->proportionalTrace(*b);
}
LevelSetMappedSegment::QuarticValueVariation LevelSetMappedSegment::nativeQuarticValueVariation(
    bool fromEnd,long double offset,Quantum::Geometry::CornerSectorVariation sectors) const
{
    const auto& retained=checked();
    const auto* query=std::get_if<Impl::NativeQuarticQuery>(&retained.query);
    if (!query) throw std::invalid_argument("trace value variation requires native quartic data");
    Quantum::Geometry::Detail::ExactExpansion::validateEnvironment();
    return Quantum::Geometry::Detail::MappedHermiteSegment::valueVariation(
        *query,*retained.source->checked().nativeQuartic,fromEnd,offset,sectors);
}
bool LevelSetMappedSegment::sameNativeCubicTraceZeros(const LevelSetMappedSegment& other) const
{
    const auto* a=std::get_if<Impl::NativeQuery>(&checked().query);
    const auto* b=std::get_if<Impl::NativeQuery>(&other.checked().query);
    if (!a || !b || begin().x!=other.begin().x || begin().y!=other.begin().y
        || end().x!=other.end().x || end().y!=other.end().y)
        throw std::invalid_argument("native cubic trace comparison requires the same oriented segment");
    Quantum::Geometry::Detail::ExactExpansion::validateEnvironment();
    return a->proportionalTrace(*b);
}
LevelSetMappedSegment::CubicValueVariation LevelSetMappedSegment::nativeCubicValueVariation(
    bool fromEnd,long double offset,Quantum::Geometry::CornerSectorVariation sectors) const
{
    const auto& retained=checked();
    const auto* query=std::get_if<Impl::NativeQuery>(&retained.query);
    if (!query) throw std::invalid_argument("trace value variation requires native cubic data");
    Quantum::Geometry::Detail::ExactExpansion::validateEnvironment();
    return Quantum::Geometry::Detail::MappedHermiteSegment::valueVariation(
        *query,retained.source->checked().leaf(),fromEnd,offset,sectors);
}
LevelSetMappedSegment LevelSetMappedSegment::refined(std::size_t rootIndex) const
{
    auto candidate=std::make_unique<Impl>(checked());
    std::visit([&](auto& query) { query.refineRoot(rootIndex); },candidate->query);
    return LevelSetMappedSegment(std::move(candidate));
}

struct LevelSetMappedSegment::NativeQuarticTrace::Impl
{
    using Query=LevelSetMappedSegment::Impl::NativeQuarticQuery;
    using Scalar=Quantum::Geometry::Detail::LevelSetHermiteScalar::QuarticLeafData;
    using Ratio=Quantum::Geometry::Detail::MappedHermiteSegment::Ratio;
    struct SideWork {std::size_t zero;int side,order,native;Query query;};
    struct Binding {int work;std::size_t root;Ratio parameter;};
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source;
    Query value;
    int axis;
    std::vector<SideWork> sides;
    std::vector<CoincidentRange> coincident;
    std::vector<Event> events;
    std::vector<Binding> bindings;
    std::unique_ptr<LevelSetMappedSegment> legacy;

    Impl(std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> retained,const Scalar& scalar,
        Vector2D a,Vector2D b,int normalAxis)
        : source(std::move(retained)),value(scalar,a,b,
            Quantum::Geometry::Detail::MappedHermiteSegment::NativeQuarticCellLayout(scalar),true),axis(normalAxis)
    {
        namespace M=Quantum::Geometry::Detail::MappedHermiteSegment;
        const bool reverse=reversed();
        const Vector2D first=reverse?b:a,last=reverse?a:b;
        const double line=axis==0?a.x:a.y;
        for(std::size_t z=0;z<value.zeroTraces.size();++z){
            const auto& zero=value.zeroTraces[z];
            for(int side:{-1,1}){
                M::CellAddress owner;int count=0;
                for(const auto& cell:zero.cells){
                    const auto box=scalar.leafBounds(cell.native);
                    const double lo=axis==0?box[0].x:box[0].y,hi=axis==0?box[1].x:box[1].y;
                    if(side<0?(line>lo && line<=hi):(line>=lo && line<hi)){owner=cell;++count;}
                }
                M::require(count<=1,"normal trace side has overlapping active owners");
                if(count==0)continue;
                int order=1;
                for(;order<=4;++order){
                    M::NativeQuarticCellLayout layout(scalar,owner,zero.chart,axis,order);
                    layout.cuts(first,last);
                    const auto p=layout.polynomial(owner,M::constant(0),M::constant(1));
                    if(std::any_of(p.begin(),p.end(),[](const auto& c){return M::Exact::sign(c)!=0;})){
                        sides.push_back({z,side,order,owner.native,Query(scalar,a,b,std::move(layout))});break;
                    }
                }
                M::require(order<=4,"all normal Taylor coefficients vanish in an active native leaf");
            }
        }
        finish();
    }
    Impl(const Impl& other)
        : source(other.source),value(other.value),axis(other.axis),sides(other.sides) {}
    bool reversed() const{return std::make_pair(value.inputEnd.x,value.inputEnd.y)
        <std::make_pair(value.inputBegin.x,value.inputBegin.y);}
    static Ratio compose(const Quantum::Geometry::Detail::MappedHermiteSegment::Chart& c,const Ratio& t)
    {
        namespace M=Quantum::Geometry::Detail::MappedHermiteSegment;
        const auto w=M::subtract(c.end,c.begin);
        return M::add(c.begin,M::ratio(M::Exact::multiply(w.numerator,t.numerator),
            M::Exact::multiply(w.denominator,t.denominator)));
    }
    Range interval(const Ratio& a,const Ratio& b,bool negative) const
    {
        namespace M=Quantum::Geometry::Detail::MappedHermiteSegment;
        M::require(M::compare(a,b)<0,"native trace range has collapsed");
        const bool fromEnd=M::compare(M::add(a,b),M::constant(1))>0;
        const auto o=M::approximate(fromEnd?M::subtract(M::constant(1),b):a),w=M::approximate(M::subtract(b,a));
        M::require(o>=0 && o<=1 && w>0 && w<=1,"native trace integration chart is unrepresentable");
        return {fromEnd!=reversed(),o,w,negative};
    }
    void appendEvents(const Query& query,int work,std::size_t zero,int side,int order)
    {
        const bool reverse=reversed();
        for(std::size_t k=0;k<query.rootInfo.size();++k){
            const auto& r=query.roots[reverse?query.roots.size()-1-k:k];
            const auto p=work<0?r.parameter:compose(value.zeroTraces[zero].chart,r.parameter);
            events.push_back({query.rootInfo[k],order,side,zero});bindings.push_back({work,k,p});
        }
    }
    void finish()
    {
        coincident.clear();events.clear();bindings.clear();legacy.reset();
        for(std::size_t z=0;z<value.zeroTraces.size();++z){const auto& c=value.zeroTraces[z].chart;
            const auto r=interval(c.begin,c.end,false);
            coincident.push_back({z,{r.fromEnd,r.offset,r.length},{}, {}});}
        appendEvents(value,-1,0,0,0);
        for(std::size_t w=0;w<sides.size();++w){const auto& work=sides[w];
            auto& out=work.side<0?coincident[work.zero].minus:coincident[work.zero].plus;
            out.available=true;out.nativeCellIndex=work.native;out.normalOrder=work.order;
            const auto& chart=value.zeroTraces[work.zero].chart;
            for(const auto& phase:work.query.phases){const bool negative=phase.negative!=(work.side<0 && work.order%2==1);
                out.ranges.push_back(interval(compose(chart,phase.begin),compose(chart,phase.end),negative));}
            if(reversed())std::reverse(out.ranges.begin(),out.ranges.end());
            appendEvents(work.query,static_cast<int>(w),work.zero,work.side,work.order);
        }
        if(coincident.empty())legacy=std::unique_ptr<LevelSetMappedSegment>(new LevelSetMappedSegment(
            std::make_unique<LevelSetMappedSegment::Impl>(source,value)));
    }
};
const LevelSetMappedSegment::NativeQuarticTrace::Impl& LevelSetMappedSegment::Impl::retainedTrace(
    const NativeQuarticTrace& trace) { return trace.checked(); }
LevelSetMappedSegment::NativeQuarticTrace LevelSetMappedSegment::nativeQuarticTrace(
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source,Vector2D a,Vector2D b,int axis)
{
    if(!source)throw std::invalid_argument("[NativeQuarticTrace] null retained source");
    Quantum::Geometry::Detail::ExactExpansion::validateEnvironment();
    if((axis!=0 && axis!=1) || (axis==0?(a.x!=b.x || a.y==b.y):(a.y!=b.y || a.x==b.x)))
        throw std::invalid_argument("[NativeQuarticTrace] normal axis requires a perpendicular nonzero axis-aligned segment");
    const auto& scalar=source->checked().nativeQuartic;
    if(!scalar)throw std::invalid_argument("[NativeQuarticTrace] native Q4 source required");
    return NativeQuarticTrace(std::make_unique<NativeQuarticTrace::Impl>(std::move(source),*scalar,a,b,axis));
}
LevelSetMappedSegment::NativeQuarticTrace::NativeQuarticTrace(std::unique_ptr<const Impl> impl):impl_(std::move(impl)){}
LevelSetMappedSegment::NativeQuarticTrace::~NativeQuarticTrace()=default;
LevelSetMappedSegment::NativeQuarticTrace::NativeQuarticTrace(NativeQuarticTrace&&) noexcept=default;
LevelSetMappedSegment::NativeQuarticTrace& LevelSetMappedSegment::NativeQuarticTrace::operator=(NativeQuarticTrace&&) noexcept=default;
const LevelSetMappedSegment::NativeQuarticTrace::Impl& LevelSetMappedSegment::NativeQuarticTrace::checked() const
{
    if(!impl_)throw std::runtime_error("[NativeQuarticTrace] query on a moved-from trace");
    Quantum::Geometry::Detail::ExactExpansion::validateEnvironment();return *impl_;
}
bool LevelSetMappedSegment::NativeQuarticTrace::usesSource(const LevelSetHermiteGeometry::Snapshot& source) const{return checked().source.get()==&source;}
Vector2D LevelSetMappedSegment::NativeQuarticTrace::begin() const{return checked().value.inputBegin;}
Vector2D LevelSetMappedSegment::NativeQuarticTrace::end() const{return checked().value.inputEnd;}
int LevelSetMappedSegment::NativeQuarticTrace::normalAxis() const{return checked().axis;}
std::size_t LevelSetMappedSegment::NativeQuarticTrace::gridPieceCount() const{return checked().value.gridPieces;}
const std::vector<LevelSetMappedSegment::Range>& LevelSetMappedSegment::NativeQuarticTrace::ordinaryRanges() const{return checked().value.ranges;}
const std::vector<LevelSetMappedSegment::NativeQuarticTrace::CoincidentRange>& LevelSetMappedSegment::NativeQuarticTrace::coincidentRanges() const{return checked().coincident;}
const std::vector<LevelSetMappedSegment::NativeQuarticTrace::Event>& LevelSetMappedSegment::NativeQuarticTrace::events() const{return checked().events;}
bool LevelSetMappedSegment::NativeQuarticTrace::rootAtInputEndpoint(std::size_t index,bool atEnd) const
{
    namespace M=Quantum::Geometry::Detail::MappedHermiteSegment;const auto& q=checked();
    M::require(index<q.bindings.size(),"native trace event index is out of range");
    return M::compare(q.bindings[index].parameter,M::constant(atEnd!=q.reversed()?1:0))==0;
}
LevelSetMappedSegment::NativeQuarticTrace LevelSetMappedSegment::NativeQuarticTrace::refined(std::size_t index) const
{
    namespace M=Quantum::Geometry::Detail::MappedHermiteSegment;const auto& q=checked();
    M::require(index<q.bindings.size(),"native trace event refinement index is out of range");
    auto candidate=std::make_unique<Impl>(q);const auto binding=q.bindings[index];
    if(binding.work<0)candidate->value.refineRoot(binding.root);
    else candidate->sides[binding.work].query.refineRoot(binding.root);
    candidate->finish();return NativeQuarticTrace(std::move(candidate));
}
const LevelSetMappedSegment& LevelSetMappedSegment::NativeQuarticTrace::ordinary() const
{
    const auto& q=checked();if(!q.legacy)throw std::runtime_error("[NativeQuarticTrace] coincident trace has no legacy phase view");
    return *q.legacy;
}

std::vector<LevelSetMappedSegment::NegativePhaseVariation> LevelSetMappedSegment::
    nativeNegativePhaseVariation(Quantum::Geometry::CornerSectorVariation sectors) const
{
    const auto &retained = checked();
    const auto *query = std::get_if<Impl::NativeQuery>(&retained.query);
    if (!query)
        throw std::runtime_error("[LevelSetMappedSegment] legacy phase variation requires cubic native data; no 25-to-16 conversion");
    return Quantum::Geometry::Detail::MappedHermiteSegment::negativePhaseVariation(
        *query, retained.source->checked().leaf(), sectors);
}

std::vector<LevelSetMappedSegment::QuarticNegativePhaseVariation> LevelSetMappedSegment::
    nativeQuarticNegativePhaseVariation(Quantum::Geometry::CornerSectorVariation sectors) const
{
    const auto &retained = checked();
    const auto *query = std::get_if<Impl::NativeQuarticQuery>(&retained.query);
    if (!query)
        throw std::runtime_error("[LevelSetMappedSegment] quartic phase variation requires native quartic data");
    return Quantum::Geometry::Detail::MappedHermiteSegment::negativePhaseVariation(
        *query, *retained.source->checked().nativeQuartic, sectors);
}

struct LevelSetInterfacePanels::NativeGraphChain::SamplingData
{
    using Sample=Quantum::Geometry::Detail::MappedHermiteInterface::GraphSample;
    std::vector<std::function<Sample(double,double,double)>> graphs;
    // Query-only results on this immutable source/chain. No field, velocity,
    // charge or mesh pointer, and no error/partially computed result is kept.
    mutable std::mutex projectionMutex;
    mutable std::map<std::array<std::uint64_t,2>,LocalProjection> projections;
    mutable std::map<std::array<std::uint64_t,2>,std::vector<LocalProjection>> loadCandidates;
    mutable std::size_t retainedLoadCandidates=0;
};

LevelSetInterfacePanels::LocalProjection LevelSetInterfacePanels::NativeGraphChain::projectedFrame(
    Vector2D query) const
{
    namespace Q=Quantum::Geometry::Detail::InterfaceQuadrature;
    Q::require(source_ && !graphs_.empty(),"native load projection requires an owned retained curve");
    Quantum::Geometry::Detail::ExactExpansion::validateEnvironment();
    Q::require(static_cast<bool>(sampling_),"native graph lost its prepared sampling geometry");
    std::array<std::uint64_t,2> key;
    std::memcpy(&key[0],&query.x,sizeof(double));
    std::memcpy(&key[1],&query.y,sizeof(double));
    {
        const std::lock_guard<std::mutex> lock(sampling_->projectionMutex);
        const auto found=sampling_->projections.find(key);
        if(found!=sampling_->projections.end()) return found->second;
    }
    // Original search and every original acceptance gate remain in one body.
    // Compute outside the lock: independent faces must not serialize. Two
    // simultaneous misses may repeat work but cannot publish incomplete data.
    auto result=projectedFrameUncached(query);
    {
        const std::lock_guard<std::mutex> lock(sampling_->projectionMutex);
        constexpr std::size_t maximumRetainedQueries=4096;
        if(sampling_->projections.size()<maximumRetainedQueries)
            sampling_->projections.emplace(key,result);
    }
    return result;
}

std::vector<LevelSetInterfacePanels::LocalProjection>
LevelSetInterfacePanels::NativeGraphChain::loadProjectionCandidates(Vector2D query) const
{
    namespace Q=Quantum::Geometry::Detail::InterfaceQuadrature;
    Q::require(source_ && !graphs_.empty(),"native load projection requires an owned retained curve");
    Quantum::Geometry::Detail::ExactExpansion::validateEnvironment();
    Q::require(static_cast<bool>(sampling_),"native graph lost its prepared sampling geometry");
    std::array<std::uint64_t,2> key;
    std::memcpy(&key[0],&query.x,sizeof(double));
    std::memcpy(&key[1],&query.y,sizeof(double));
    {
        const std::lock_guard<std::mutex> lock(sampling_->projectionMutex);
        const auto found=sampling_->loadCandidates.find(key);
        if(found!=sampling_->loadCandidates.end()) return found->second;
    }
    // Preserve the entire original stationary/endpoint inventory, including
    // competing branches. Neither a chosen nearest point nor a force-provider
    // decision can replace it. Failed searches publish nothing; independent
    // faces compute outside the lock, as for projectedFrame above.
    auto result=loadProjectionCandidatesUncached(query);
    {
        const std::lock_guard<std::mutex> lock(sampling_->projectionMutex);
        constexpr std::size_t maximumQueries=4096,maximumCandidates=16384;
        if(sampling_->loadCandidates.size()<maximumQueries
            && result.size()<=maximumCandidates-sampling_->retainedLoadCandidates) {
            const auto inserted=sampling_->loadCandidates.emplace(key,result);
            if(inserted.second)sampling_->retainedLoadCandidates+=result.size();
        }
    }
    return result;
}

LevelSetInterfacePanels::NativeGraphChain LevelSetInterfacePanels::discoverNativeGraphs(
    const QuadtreeMesh& mesh,std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source,
    bool axisymmetric,const std::vector<Vector2D>& boundaryContacts)
{
    auto components = discoverNativeComponentsImpl(mesh, std::move(source), axisymmetric, boundaryContacts, true);
    return std::move(components.curves.front());
}

LevelSetInterfacePanels::NativeComponents LevelSetInterfacePanels::discoverNativeComponents(
    const QuadtreeMesh& mesh,std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source,
    bool axisymmetric,const std::vector<Vector2D>& boundaryContacts,const std::vector<SymmetryPlane>& symmetryPlanes)
{
    return discoverNativeComponentsImpl(mesh, std::move(source), axisymmetric, boundaryContacts, false,symmetryPlanes);
}
LevelSetInterfacePanels::NativeComponents LevelSetInterfacePanels::discoverNativeTraceComponents(
    const QuadtreeMesh& mesh,std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source,
    bool axisymmetric,const std::vector<Vector2D>& boundaryContacts,const std::vector<SymmetryPlane>& symmetryPlanes)
{
    if(!source || !source->checked().nativeQuartic)
        throw std::invalid_argument("[LevelSetInterfacePanels] native trace components require an owned native Q4 snapshot");
    return discoverNativeComponentsImpl(mesh,std::move(source),axisymmetric,boundaryContacts,false,symmetryPlanes,true);
}

LevelSetInterfacePanels::NativeComponents LevelSetInterfacePanels::discoverNativeComponentsImpl(
    const QuadtreeMesh& mesh,std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source,
    bool axisymmetric,const std::vector<Vector2D>& boundaryContacts, bool legacy,const std::vector<SymmetryPlane>& symmetryPlanes,bool nativeTrace)
{
    if(!source) throw std::invalid_argument("[LevelSetInterfacePanels] null native interface source");
    auto all=source->checked().visitLeaf([&](const auto& scalar) {
        return Quantum::Geometry::Detail::NativeInterfaceChain::discoverAll(mesh,scalar,axisymmetric,boundaryContacts,legacy,symmetryPlanes,nativeTrace);
    });
    NativeComponents collection;
    collection.negativeRegions = std::move(all.negativeRegions);
    collection.cellPhase = std::move(all.cellPhase);
    auto& results = collection.curves;
    results.reserve(all.curves.size());
    for (const auto& discovered : all.curves)
    {
        NativeGraphChain result;
        result.source_=source; result.graphs_=discovered.graphs; result.mapping_=discovered.mapping;
        for(std::size_t g=0;g<result.graphs_.size();++g){const auto& s=discovered.sideSupports.at(g);if(s.owner<0)continue;
            result.graphs_[g].nativeSide=std::shared_ptr<const LevelSetMappedInterface::NativeSideSupport>(
                new LevelSetMappedInterface::NativeSideSupport(source,source->checked().nativeQuartic.get(),s.owner,s.side,
                    s.horizontal,s.edgeLabel,s.begin,s.end,s.beginVertex,s.endVertex,s.minus,s.plus));
        }
        result.cellGraphs_=discovered.cellGraphs;
        result.negativeRegion_=discovered.negativeRegion;
        result.rootNodes_=discovered.rootNodes; result.minimumSpacing_=discovered.minimumSpacing;
        result.boundaryContacts_=discovered.boundaryContacts;
        result.symmetryPlanes_=symmetryPlanes;result.endpointSymmetryPlanes_=discovered.endpointSymmetryPlanes;
        result.closed_=discovered.closed; result.axisymmetric_=axisymmetric;
        result.beginKind_=discovered.begin; result.endKind_=discovered.end;
        auto sampling=std::make_shared<NativeGraphChain::SamplingData>();
        sampling->graphs.reserve(result.graphs_.size());
        result.source_->checked().visitLeaf([&](const auto& scalar) {
            using Scalar=std::decay_t<decltype(scalar)>;
            using Integrator=Quantum::Geometry::Detail::MappedHermiteInterface::NativeIntegratorFor<Scalar>;
            // The alias owns the opaque snapshot, not a borrowed visitor argument.
            // The prepared object owns its Graph/Cuts and has no chain-vector borrow.
            const std::shared_ptr<const Scalar> retained(result.source_,&scalar);
            for (const auto& graph:result.graphs_)
            {
                auto prepared=Integrator::prepareGraphSampling(scalar,graph);
                sampling->graphs.emplace_back([retained,prepared=std::move(prepared)](double begin,double end,double fraction) {
                    return Integrator::graphSample(*retained,prepared,begin,end,fraction);
                });
            }
        });
        result.sampling_=std::move(sampling);
        results.push_back(std::move(result));
    }
    return collection;
}

LevelSetHermiteGeometry::Frame LevelSetInterfacePanels::NativeGraphChain::frameAt(
    int graphIndex,double fraction) const
{
    // Reuse only exact fixed fractions on this retained immutable source.
    // Discovery remains lazy; all other fractions retain the original query.
    Quantum::Geometry::Detail::ExactExpansion::validateEnvironment();
    if(source_ && graphIndex>=0 && static_cast<std::size_t>(graphIndex)<graphs_.size() &&
        fixedProjectionFrames_.size()==graphs_.size())
    {
        const auto& frames=fixedProjectionFrames_[graphIndex];
        if(fraction==0) return frames[0];
        if(fraction==1) return frames[1];
        const auto& nodes=Quantum::Geometry::Detail::InterfaceQuadrature::gaussX;
        for(std::size_t k=0;k<nodes.size();++k)
            if(fraction==(1.0+nodes[k])/2.0) return frames[k+2];
    }
    return frameAtInterval(graphIndex,0,1,fraction);
}

LevelSetHermiteGeometry::Frame LevelSetInterfacePanels::NativeGraphChain::frameAtInterval(
    int graphIndex,double begin,double end,double fraction) const
{
    namespace Q=Quantum::Geometry::Detail::InterfaceQuadrature;
    Q::require(source_ && graphIndex>=0 && static_cast<std::size_t>(graphIndex)<graphs_.size(),
        "native frame requires an owned source and valid graph index");
    const auto& graph=graphs_[graphIndex];
    Q::require(sampling_ && sampling_->graphs.size()==graphs_.size(),"native graph lost its prepared sampling geometry");
    const auto sampled=sampling_->graphs[graphIndex](begin,end,fraction);
    auto frame=source_->jetAtNativeCellChart(graph.nativeCellIndex,sampled.chart);
    const double component=graph.coordinateIsX?-frame.normal.y:frame.normal.x;
    Q::require(std::isfinite(component) && component*sampled.coordinateDirection>
        256.0*std::numeric_limits<double>::epsilon(),"chosen native graph has an unresolved or reversed tangent");
    const bool endpoint=(begin==0 && fraction==0) || (end==1 && fraction==1);
    if(!endpoint) {
        if(!graph.nativeSide)Q::require(Quantum::Geometry::cellChartInterior(sampled.chart),
            "native graph sample is not strictly inside its local cell");
        else {
            const auto& side=*graph.nativeSide;const auto tangent=graph.coordinateIsX?sampled.chart.x:sampled.chart.y;
            const auto normal=graph.coordinateIsX?sampled.chart.y:sampled.chart.x;
            Q::require(side.usesSource(*source_) && Quantum::Geometry::compareCellChartCoordinate(
                normal,{side.side()<2?0.:1.,false})==0 && Quantum::Geometry::compareCellChartCoordinate(tangent,{0,false})>0
                && Quantum::Geometry::compareCellChartCoordinate(tangent,{1,false})<0,
                "native side sample left its certified support/tangent interior");
        }
    }
    if(axisymmetric_)
    {
        const auto bounds=source_->nativeCellBounds(graph.nativeCellIndex);
        const bool atAxis=bounds[0].x==0 && Quantum::Geometry::compareCellChartCoordinate(
            sampled.chart.x,{0,false})==0;
        if(atAxis && endpoint)
            Q::require(source_->radialAxisSymmetry() && frame.normal.x==0,"native axis endpoint is not symmetric");
        else Q::offAxisWeight(frame);
    }
    return frame;
}

std::vector<LevelSetInterfacePanels::ArcSample> LevelSetInterfacePanels::NativeGraphChain::quadrature(
    int graphIndex) const
{ return quadrature(graphIndex,0,1); }

std::vector<LevelSetInterfacePanels::ArcSample> LevelSetInterfacePanels::NativeGraphChain::quadrature(
    int graphIndex,double begin,double end) const
{
    namespace Q=Quantum::Geometry::Detail::InterfaceQuadrature;
    Q::require(source_ && graphIndex>=0 && static_cast<std::size_t>(graphIndex)<graphs_.size(),
        "native quadrature requires an owned source and valid graph index");
    Q::require(std::isfinite(begin) && std::isfinite(end) && begin>=0 && end<=1 && begin<=end,
        "invalid native subarc interval");
    if(begin==end) return {};
    const auto& graph=graphs_[graphIndex];
    Q::require(sampling_ && sampling_->graphs.size()==graphs_.size(),"native graph lost its prepared sampling geometry");
    const auto middle=sampling_->graphs[graphIndex](begin,end,.5);
    return Q::sampleGraph(middle.physicalSpan,axisymmetric_,
        [&](double fraction) { return frameAtInterval(graphIndex,begin,end,fraction); },
        [&](ArcSample& node) {
            node.localGraphCoordinate=std::numeric_limits<double>::quiet_NaN();
            node.localGraphChart=graph.coordinateIsX?node.frame.cellChartPoint.x:node.frame.cellChartPoint.y;
            node.hasLocalGraphChart=true;
        },
        [&](const LevelSetHermiteGeometry::Frame& frame) {
            const double component=graph.coordinateIsX?-frame.normal.y:frame.normal.x;
            Q::require(std::isfinite(component) && component*middle.coordinateDirection>
                256.0*std::numeric_limits<double>::epsilon(),"chosen native graph has an unresolved or reversed tangent");
            return Q::checked(1.0L/std::abs(component),"graph arc metric overflow");
        },
        [&](const LevelSetHermiteGeometry::Frame& frame) {
            return source_->scalarJetAtNativeCellChart(frame.nativeCellIndex,frame.cellChartPoint);
        });
}

double LevelSetInterfacePanels::NativeGraphChain::physicalCoordinateSpan(int graphIndex) const
{
    namespace Q=Quantum::Geometry::Detail::InterfaceQuadrature;
    Q::require(source_ && graphIndex>=0 && static_cast<std::size_t>(graphIndex)<graphs_.size(),
        "native span requires an owned source and valid graph index");
    Q::require(sampling_ && sampling_->graphs.size()==graphs_.size(),"native graph lost its prepared sampling geometry");
    const auto sample=sampling_->graphs[graphIndex](0,1,.5);
    return Q::checked(sample.physicalSpan,"native material coordinate span is unrepresentable");
}

double LevelSetInterfacePanels::NativeGraphChain::graphFraction(
    int graphIndex,const Quantum::Geometry::CellChartPoint& point) const
{
    namespace Q=Quantum::Geometry::Detail::InterfaceQuadrature;
    Q::require(source_ && graphIndex>=0 && static_cast<std::size_t>(graphIndex)<graphs_.size(),
        "native coordinate requires an owned source and valid graph index");
    return source_->checked().visitLeaf([&](const auto& scalar) {
        return Quantum::Geometry::Detail::MappedHermiteInterface::NativeIntegratorFor<
            std::decay_t<decltype(scalar)>>::graphFraction(scalar,graphs_[graphIndex],point);
    });
}

LevelSetMappedInterface::Result LevelSetMappedInterface::integrate(
    const LevelSetHermiteGeometry::Snapshot& source, const Graph& graph,
    const std::vector<const LevelSetMappedSegment*>& faceCuts,
    const std::vector<Region>& regions, bool axisymmetric)
{
    if(graph.nativeSide)throw std::invalid_argument("legacy sensitivity does not admit native side-attached support");
    Quantum::Geometry::Detail::ExactExpansion::validateEnvironment();
    const auto& backend = source.checked();
    if (backend.quartic || backend.nativeQuartic)
        throw std::invalid_argument("[LevelSetMappedInterface] EnrichedQuartic requires integrateQuartic and 25 data columns");
    std::vector<const Quantum::Geometry::Detail::MappedHermiteSegment::Query*> retainedCuts;
    std::vector<const Quantum::Geometry::Detail::MappedHermiteSegment::NativeQuery*> nativeCuts;
    retainedCuts.reserve(faceCuts.size());
    if (backend.native) nativeCuts.reserve(faceCuts.size());
    for (const auto* cut : faceCuts)
    {
        if (!cut) throw std::invalid_argument("[LevelSetMappedInterface] null retained face cut");
        const auto& retained = cut->checked();
        if (retained.source.get() != &source)
            throw std::invalid_argument("[LevelSetMappedInterface] face cut uses a different snapshot");
        if (backend.native)
        {
            const auto* query=std::get_if<LevelSetMappedSegment::Impl::NativeQuery>(&retained.query);
            if (!query) throw std::invalid_argument("[LevelSetMappedInterface] native source requires native face cuts");
            nativeCuts.push_back(query);
        }
        else retainedCuts.push_back(&retained.uniformQuery());
    }
    if (backend.native)
    {
        auto result = Quantum::Geometry::Detail::MappedHermiteInterface::NativeIntegrator(
            *backend.native, graph, nativeCuts, regions, axisymmetric).run();
        for (auto& arc : result.arcs)
        {
            arc.independentSensitivity = backend.native->pullbackCellSensitivity(
                arc.nativeCellIndex, arc.sensitivity);
            arc.requiresNativeDofConstraints = false;
        }
        return result;
    }
    return Quantum::Geometry::Detail::MappedHermiteInterface::Integrator(
        backend.cubicGeometry(), graph, retainedCuts, regions, axisymmetric).run();
}

LevelSetHermiteGeometry::Snapshot::QuarticCenterInterval
LevelSetHermiteGeometry::Snapshot::nativeQuarticCenterInterval(int index,double displacementCells,
    double retainedFraction,QuarticCenterDirection direction) const
{
    const auto& data=checked();
    if(!data.nativeQuartic) throw std::invalid_argument("Q4 center interval requires a native quartic scalar");
    return data.nativeQuartic->centerInterval(index,displacementCells,retainedFraction,direction);
}

int LevelSetHermiteGeometry::Snapshot::nativeQuarticBoxSign(int index,Vector2D lower,Vector2D upper) const
{
    const auto& data=checked();
    if(!data.nativeQuartic)throw std::invalid_argument("Q4 box sign requires a native quartic scalar");
    return data.nativeQuartic->boxSign(index,lower,upper);
}

LevelSetHermiteGeometry::Snapshot::QuarticChangeCertificate
LevelSetHermiteGeometry::Snapshot::nativeQuarticChangeCertificate(int index,const Snapshot& other,double budget) const
{
    const auto& a=checked();const auto& b=other.checked();
    if(!a.nativeQuartic || !b.nativeQuartic)throw std::invalid_argument("Q4 change certificate requires two native quartic scalars");
    return a.nativeQuartic->changeCertificate(index,*b.nativeQuartic,budget);
}

LevelSetMappedInterface::GraphMoments LevelSetMappedInterface::nativeMoments(
    const LevelSetHermiteGeometry::Snapshot& source,const Graph& graph,bool axisymmetric)
{
    if(graph.nativeSide && !graph.nativeSide->usesSource(source))
        throw std::invalid_argument("native side moments use a different snapshot");
    Quantum::Geometry::Detail::ExactExpansion::validateEnvironment();
    const auto& backend=source.checked();
    if (graph.nativeCellIndex<0 || (!backend.native && !backend.nativeQuartic))
        throw std::invalid_argument("nativeMoments requires an original native cubic/Q4 graph");
    if (backend.native) return Quantum::Geometry::Detail::MappedHermiteInterface::NativeIntegrator::graphMoments(
        *backend.native,graph,axisymmetric);
    return Quantum::Geometry::Detail::MappedHermiteInterface::NativeQuarticIntegrator::graphMoments(
        *backend.nativeQuartic,graph,axisymmetric);
}

int LevelSetInterfacePanels::NativeGraphChain::graphAt(
    int nativeCellIndex,const Quantum::Geometry::CellChartPoint& point) const
{
    namespace Q=Quantum::Geometry::Detail::InterfaceQuadrature;
    namespace C=Quantum::Geometry;
    C::validateCellChart(point);
    Q::require(source_ && nativeCellIndex>=0 && static_cast<std::size_t>(nativeCellIndex)<cellGraphs_.size(),
        "native graph lookup requires an owned source and valid leaf index");
    int found=-1;bool foundEndpoint=false;
    for(int g:cellGraphs_[nativeCellIndex]) {
        const auto& graph=graphs_.at(g);
        const auto inside=[](C::CellChartCoordinate value,double lo,double hi) {
            return C::compareCellChartCoordinate(value,{lo,false})>=0
                && C::compareCellChartCoordinate(value,{hi,false})<=0;
        };
        if(!inside(point.x,graph.localLower.x,graph.localUpper.x)
            || !inside(point.y,graph.localLower.y,graph.localUpper.y))continue;
        const auto first=frameAt(g,0).cellChartPoint,last=frameAt(g,1).cellChartPoint;
        const bool endpoint=C::cellChartEqual(point,first) || C::cellChartEqual(point,last);
        if(!endpoint) {
            const auto q=graph.coordinateIsX?point.x:point.y;
            const auto a=graph.coordinateIsX?first.x:first.y,b=graph.coordinateIsX?last.x:last.y;
            const int from=C::compareCellChartCoordinate(q,a),to=C::compareCellChartCoordinate(q,b);
            if(!((from>0 && to<0)||(from<0 && to>0)))continue;
        }
        if(found>=0) {
            Q::require(endpoint && foundEndpoint && (std::abs(g-found)==1
                || (closed_ && std::abs(g-found)==static_cast<int>(graphs_.size())-1)),
                "native chart point has competing retained graph incidences");
            found=std::min(found,g);
        } else {found=g;foundEndpoint=endpoint;}
    }
    Q::require(found>=0,"native chart point has no retained graph incidence");
    graphFraction(found,point); // Preserve the represented-interior and root identity gates.
    return found;
}

LevelSetMappedInterface::QuarticResult LevelSetMappedInterface::integrateQuartic(
    const LevelSetHermiteGeometry::Snapshot& source, const Graph& graph,
    const std::vector<const LevelSetMappedSegment*>& faceCuts,
    const std::vector<Region>& regions, bool axisymmetric)
{
    return integrateQuartic(source,graph,faceCuts,regions,axisymmetric,{0,0});
}

LevelSetMappedInterface::QuarticResult LevelSetMappedInterface::integrateQuartic(
    const LevelSetHermiteGeometry::Snapshot& source, const Graph& graph,
    const std::vector<const LevelSetMappedSegment*>& faceCuts,
    const std::vector<Region>& regions, bool axisymmetric, Vector2D greenOrigin)
{
    if(graph.nativeSide)throw std::invalid_argument("legacy quartic sensitivity does not admit native side-attached support");
    Quantum::Geometry::Detail::ExactExpansion::validateEnvironment();
    const auto& backend=source.checked();
    if (!backend.quartic && !backend.nativeQuartic)
        throw std::invalid_argument("[LevelSetMappedInterface] integrateQuartic requires an EnrichedQuartic snapshot");
    const auto integrate=[&](const auto& scalar)
    {
        namespace M=Quantum::Geometry::Detail::MappedHermiteSegment;
        namespace I=Quantum::Geometry::Detail::MappedHermiteInterface;
        constexpr bool native=std::is_same_v<std::decay_t<decltype(scalar)>,
            Quantum::Geometry::Detail::LevelSetHermiteScalar::QuarticLeafData>;
        using Query=std::conditional_t<native,M::NativeQuarticQuery,M::QuarticQuery>;
        using Integrator=std::conditional_t<native,I::NativeQuarticIntegrator,I::QuarticIntegrator>;
        std::vector<const Query*> cuts;
        cuts.reserve(faceCuts.size());
        for (const auto* cut:faceCuts)
        {
            if (!cut) throw std::invalid_argument("[LevelSetMappedInterface] null retained face cut");
            const auto& retained=cut->checked();
            if (retained.source.get()!=&source)
                throw std::invalid_argument("[LevelSetMappedInterface] face cut uses a different snapshot");
            const auto* query=std::get_if<Query>(&retained.query);
            if (!query) throw std::invalid_argument("[LevelSetMappedInterface] quartic source requires matching quartic face cuts");
            cuts.push_back(query);
        }
        return Integrator(scalar,graph,cuts,regions,axisymmetric,greenOrigin).run();
    };
    if(backend.nativeQuartic) return integrate(*backend.nativeQuartic);
    return integrate(*backend.quartic);
}

LevelSetMappedInterface::NativeRegionMoments LevelSetMappedInterface::integrateNativeQuarticRegionMoments(
    const LevelSetHermiteGeometry::Snapshot& source,int cell,const std::vector<Graph>& graphs,
    const std::vector<const LevelSetMappedSegment::NativeQuarticTrace*>& faceCuts,
    const std::vector<Region>& regions,bool radial,Vector2D greenOrigin)
{
    namespace M=Quantum::Geometry::Detail::MappedHermiteSegment;
    namespace I=Quantum::Geometry::Detail::MappedHermiteInterface;
    constexpr long double pi=3.1415926535897932384626433832795L;
    M::Exact::validateEnvironment();
    const auto& scalar=source.checked().nativeQuartic;
    M::require(scalar!=nullptr,"native region moments require a native Q4 snapshot");
    M::require(std::isfinite(greenOrigin.x) && std::isfinite(greenOrigin.y),"nonfinite region Green origin");
    const auto bounds=scalar->leafBounds(cell); // active source cell only
    M::require(!regions.empty(),"native region moments require actual rectangles");
    std::vector<const M::NativeQuarticQuery*> cuts;
    for(const auto* cut:faceCuts) {
        M::require(cut!=nullptr,"null native region trace");
        const auto& q=LevelSetMappedSegment::Impl::retainedTrace(*cut);
        M::require(q.source.get()==&source,"native region trace uses a different snapshot");
        cuts.push_back(&q.value);
    }
    const auto equal=[](Vector2D a,Vector2D b){return a.x==b.x && a.y==b.y;};
    struct Rectangle {Vector2D lo,hi;};std::vector<Rectangle> rectangles;
    std::vector<double> xs{bounds[0].x,bounds[1].x},ys{bounds[0].y,bounds[1].y};
    NativeRegionMoments result;std::map<int,std::size_t> owners;
    for(const auto& r:regions) {
        M::require(r.regionId>=0,"negative native region id");
        Vector2D lo=r.vertices[0],hi=lo;int orientation=0;
        for(int e=0;e<4;++e) {
            const auto a=r.vertices[e],b=r.vertices[(e+1)%4];
            M::require(std::isfinite(a.x) && std::isfinite(a.y) && (!radial || a.x>=0),"invalid native region vertex");
            M::require((a.x==b.x)!=(a.y==b.y),"native region moments require nonzero axis-aligned rectangle sides");
            const int turn=I::orientation(a,b,I::point(r.vertices[(e+2)%4]));
            M::require(turn!=0 && (orientation==0 || orientation==turn),"native rectangle is not strictly convex/cyclic");orientation=turn;
            lo.x=std::min(lo.x,a.x);lo.y=std::min(lo.y,a.y);hi.x=std::max(hi.x,a.x);hi.y=std::max(hi.y,a.y);
        }
        M::require(lo.x>=bounds[0].x && lo.y>=bounds[0].y && hi.x<=bounds[1].x && hi.y<=bounds[1].y,
            "native region rectangle left its actual cell");
        for(const auto& v:r.vertices)M::require((v.x==lo.x || v.x==hi.x) && (v.y==lo.y || v.y==hi.y),"invalid rectangle corners");
        rectangles.push_back({lo,hi});xs.push_back(lo.x);xs.push_back(hi.x);ys.push_back(lo.y);ys.push_back(hi.y);
        const auto inserted=owners.emplace(r.regionId,result.regions.size());
        if(inserted.second){RegionMoments moments;moments.regionId=r.regionId;result.regions.push_back(moments);}
        result.regions[inserted.first->second].capacity+=(static_cast<long double>(hi.x)-lo.x)
            *(static_cast<long double>(hi.y)-lo.y)*(radial?pi*(static_cast<long double>(hi.x)+lo.x):1.L);
    }
    for(auto* v:{&xs,&ys}){std::sort(v->begin(),v->end());v->erase(std::unique(v->begin(),v->end()),v->end());}
    for(std::size_t y=1;y<ys.size();++y)for(std::size_t x=1;x<xs.size();++x) {
        int covers=0;for(const auto& box:rectangles)if(box.lo.x<=xs[x-1] && box.hi.x>=xs[x]
            && box.lo.y<=ys[y-1] && box.hi.y>=ys[y])++covers;
        M::require(covers==1,"native region partition has a gap or positive-area overlap");
    }
    // Fixed sides use retained ordinary/normal phases in canonical input
    // coordinates. Positive length is converted only after exact clipping.
    for(std::size_t k=0;k<regions.size();++k)for(int e=0;e<4;++e) {
        auto a=regions[k].vertices[e],b=regions[k].vertices[(e+1)%4];
        if(std::make_pair(b.x,b.y)<std::make_pair(a.x,a.y))std::swap(a,b);
        std::size_t id=faceCuts.size();
        for(std::size_t s=0;s<cuts.size();++s)if((equal(a,cuts[s]->inputBegin) && equal(b,cuts[s]->inputEnd))
            || (equal(b,cuts[s]->inputBegin) && equal(a,cuts[s]->inputEnd))) {
            M::require(id==faceCuts.size(),"duplicate actual native region trace");id=s;
        }
        M::require(id<faceCuts.size(),"native rectangle side is missing its actual complete trace");
        const auto& q=LevelSetMappedSegment::Impl::retainedTrace(*faceCuts[id]);
        const bool horizontal=a.y==b.y;const int axis=horizontal?1:0;
        M::require(q.axis==axis,"native region trace has the wrong normal axis");
        const auto& box=rectangles[k];
        const int outward=horizontal?(a.y==box.hi.y?1:-1):(a.x==box.hi.x?1:-1);
        const int inside=-outward;
        const auto delta=M::ratio(M::difference(horizontal?b.x:b.y,horizontal?a.x:a.y),M::number(1));
        const auto anchor=M::constant(horizontal?a.x:a.y);
        std::array<long double,2> weights{{0,0}},zeroWeights{{0,0}};
        const auto add=[&](const M::Ratio& begin,const M::Ratio& end,bool negative,bool zero) {
            const auto width=I::product(delta,I::difference(end,begin));
            long double measure=M::approximate(width);
            if(horizontal && radial) {
                const auto mean=I::sum(anchor,I::product(delta,I::product(I::sum(begin,end),M::constant(.5))));
                measure*=2*pi*M::approximate(mean);
            }
            M::require(std::isfinite(measure) && measure>=0,"invalid native region phase side measure");
            weights[negative?0:1]+=measure;if(zero)zeroWeights[negative?0:1]+=measure;
        };
        for(const auto& phase:q.value.phases)add(phase.begin,phase.end,phase.negative,false);
        for(std::size_t z=0;z<q.value.zeroTraces.size();++z) {
            const auto work=std::find_if(q.sides.begin(),q.sides.end(),[&](const auto& s){return s.zero==z && s.side==inside;});
            M::require(work!=q.sides.end(),"native region clip lacks its actual inside normal trace");
            M::require(work->native==cell,"native region inside zero trace belongs to a different cell");
            const auto& chart=q.value.zeroTraces[z].chart;
            for(const auto& phase:work->query.phases)add(q.compose(chart,phase.begin),q.compose(chart,phase.end),
                phase.negative!=(inside<0 && work->order%2==1),true);
        }
        const auto coordinate=M::constant(horizontal?a.y:a.x),origin=M::constant(horizontal?greenOrigin.y:greenOrigin.x);
        const auto shifted=I::difference(coordinate,origin);
        const long double potential=outward*M::approximate(horizontal || !radial?shifted:
            I::product(shifted,I::sum(coordinate,origin)))*(horizontal || !radial?1.L:pi);
        int incidences=0;bool sameOwner=true;
        for(const auto& other:regions)for(int f=0;f<4;++f) {
            const auto c=other.vertices[f],d=other.vertices[(f+1)%4];
            if((equal(a,c)&&equal(b,d)) || (equal(a,d)&&equal(b,c))) {
                ++incidences;sameOwner=sameOwner && other.regionId==regions[k].regionId;
            }
        }
        M::require(incidences==1 || incidences==2,"native fixed side has invalid incidence");
        // A seam inside the SAME union is not a fixed boundary of that
        // region. Its zero curve remains a single free interface instead.
        if(incidences==2 && sameOwner)continue;
        auto& r=result.regions[owners.at(regions[k].regionId)];
        if(horizontal){r.fixedAxialLiquid+=potential*weights[0];r.axialGas+=potential*weights[1];}
        else {r.fixedRadialLiquid+=potential*weights[0];r.radialGas+=potential*weights[1];}
        if(incidences==1) {
            if(horizontal)result.exteriorFixedAxialLiquid+=potential*weights[0];
            else result.exteriorFixedRadialLiquid+=potential*weights[0];
        } else {
            if(horizontal)result.coincidentAxialLiquid+=potential*zeroWeights[0];
            else result.coincidentRadialLiquid+=potential*zeroWeights[0];
        }
    }
    for(const auto& graph:graphs) {
        M::require(graph.nativeCellIndex==cell && graph.cellI==-1 && graph.cellJ==-1,"native region graph belongs to a different cell");
        try {
            const auto integrated=I::NativeQuarticIntegrator::regionGeometry(*scalar,graph,cuts,regions,radial,greenOrigin);
            result.coincidentArcs+=integrated.coincidentArcs;result.localRootRefinements+=integrated.localRootRefinements;
            result.sourceRootIncidences+=integrated.integral.sourceRootIncidences;
            for(const auto& arc:integrated.integral.arcs) {
                auto& r=result.regions[owners.at(arc.regionId)];
                r.freeAxialLiquid+=arc.axialMoment;r.freeRadialLiquid+=arc.radialMoment;
                r.errorEstimate+=arc.momentErrorEstimate;++r.freeArcs;
            }
        } catch(const RootResolutionError& failure) {
            auto requests=failure.requests;
            for(auto& request:requests) {
                const auto& q=LevelSetMappedSegment::Impl::retainedTrace(*faceCuts.at(request.segmentIndex));
                // Value events precede side events and use public input order.
                M::require(request.rootIndex<q.value.roots.size(),"foreign region value-root resolution request");
            }
            throw RootResolutionError(failure.what(),std::move(requests));
        }
    }
    for(auto& r:result.regions) {
        r.axialLiquid=r.fixedAxialLiquid+r.freeAxialLiquid;r.radialLiquid=r.fixedRadialLiquid+r.freeRadialLiquid;
        r.axialGas-=r.freeAxialLiquid;r.radialGas-=r.freeRadialLiquid;
    }
    return result;
}
