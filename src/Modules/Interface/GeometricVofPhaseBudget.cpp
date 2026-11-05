#include "Modules/Interface/GeometricVofPhaseBudget.h"
#include "Core/Geometry/Detail/ExactFaceTimeRational.h"
#include "Core/Geometry/Detail/ExactRationalRound.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace LAMP::Interface::Detail
{
namespace
{
namespace Exact=Geometry::Detail::ExactExpansion;
using E=Exact::Expansion;
using Rational=Geometry::Detail::FaceTime::FaceTimeRational;
constexpr long double arithmetic=512*std::numeric_limits<long double>::epsilon();
void require(bool ok,const std::string& reason)
{if(!ok)throw std::invalid_argument("[Interface/GeometricVofPhaseBudget] "+reason);}
// Literal success checks need no temporary error string. The original string
// overload constructs the identical failure only after the condition fails.
void require(bool ok,const char* reason)
{if(!ok)require(false,std::string(reason));}
E negative(const E& x){return Exact::subtract(E{},x);}
E scalar(long double x)
{
    require(std::isfinite(x),"nonfinite arithmetic input");E result;
    // Preserve extended long-double inputs as binary64 components as well.
    constexpr unsigned count=(std::numeric_limits<long double>::digits+52)/53+1;
    for(unsigned i=0;x!=0&&i<count;++i){const double part=static_cast<double>(x);
        require(std::isfinite(part)&&part!=0,"input exceeds exact binary64 exponent window");
        Exact::accumulate(result,Exact::fromDouble(part));const auto rest=x-static_cast<long double>(part);
        require(rest!=x,"input decomposition did not advance");x=rest;}
    require(x==0,"input decomposition did not retain every bit");return result;
}
long double value(const E& x)
{
    long double result=0;for(double component:x.components())result+=component;
    require(std::isfinite(result)&&(Exact::sign(x)==0||result!=0),"unrepresentable exact inventory");return result;
}
long double quotient(const E& numerator,const E& denominator)
{
    // The retained blend/ledger is rounded once with an exact midpoint
    // certificate. In particular, D*T / D exports the original binary T.
    return Geometry::Detail::ExactRationalRound::quotient(numerator,denominator);
}
struct Ratio{E n,d;};
Ratio unit(){return {Exact::fromDouble(1),Exact::fromDouble(1)};}
Ratio limited(const E& budget,const E& demand)
{
    require(Exact::sign(budget)>=0&&Exact::sign(demand)>=0,"negative phase budget/demand");
    if(Exact::sign(demand)==0||Exact::sign(Exact::subtract(budget,demand))>=0)return unit();
    const Rational normalized(budget,demand);return {normalized.n,normalized.d};
}
bool less(const Ratio& a,const Ratio& b)
{return Exact::sign(Exact::subtract(Exact::multiply(a.n,b.d),Exact::multiply(b.n,a.d)))<0;}
E slackNumerator(const Ratio& cell,const Ratio& face)
{
    auto result=Exact::subtract(Exact::multiply(cell.n,face.d),Exact::multiply(face.n,cell.d));
    require(Exact::sign(result)>=0,"negative exact synchronized limiter slack");return result;
}
struct Pair{E liquid,gas;};
Pair highPair(long double total,long double liquid,long double gas)
{
    require(total>=0&&liquid>=0&&gas>=0&&std::isfinite(total)&&std::isfinite(liquid)&&std::isfinite(gas),"invalid gross phase input");
    require(total==0?(liquid==0&&gas==0):(std::abs(liquid+gas-total)/total<=arithmetic),"raw gross phase input lost total closure");
    const auto t=scalar(total);Pair result;
    // The independent measured minor is authoritative. Reconstruct the
    // dominant complement exactly, rather than adding its rounded copy.
    if(liquid<=gas){result.liquid=scalar(liquid);result.gas=Exact::subtract(t,result.liquid);}
    else {result.gas=scalar(gas);result.liquid=Exact::subtract(t,result.gas);}
    require(Exact::sign(result.liquid)>=0&&Exact::sign(result.gas)>=0,"minor phase exceeds native total");return result;
}
struct NativeGross {E positive,negative;};
NativeGross nativeGross(const GeometricVofIntegratedFlux& raw)
{
    require(std::isfinite(raw.total)&&std::isfinite(raw.positiveTotal)&&std::isfinite(raw.negativeTotal)&&
        raw.positiveTotal>=0&&raw.negativeTotal>=0,"invalid native signed total input");
    NativeGross result{scalar(raw.positiveTotal),scalar(raw.negativeTotal)};
    const auto net=scalar(raw.total),gross=Exact::add(result.positive,result.negative);
    const auto defect=Exact::subtract(Exact::subtract(result.positive,result.negative),net);
    const auto magnitude=Exact::sign(defect)<0?negative(defect):defect;
    require(Exact::sign(Exact::subtract(magnitude,Exact::multiply(scalar(arithmetic),gross)))<=0,
        "native signed total exceeds gross arithmetic certificate");
    // Anchor the opposite branch. Reconstructing the larger branch by adding
    // a signed mean avoids manufacturing a negative roundoff counterflow.
    if(Exact::sign(net)>=0)result.positive=Exact::add(result.negative,net);
    else result.negative=Exact::subtract(result.positive,net);
    require(Exact::sign(result.positive)>=0&&Exact::sign(result.negative)>=0,"negative reconciled native gross");
    return result;
}
Pair nativeHighPair(const E& total,long double originalTotal,long double liquid,long double gas)
{
    auto result=highPair(originalTotal,liquid,gas); // Keep the original raw-coverage gate.
    if(liquid<=gas)result.gas=Exact::subtract(total,result.liquid);
    else result.liquid=Exact::subtract(total,result.gas);
    require(Exact::sign(result.liquid)>=0&&Exact::sign(result.gas)>=0,"native signed total cannot retain both measured phases");
    require((liquid==0||Exact::sign(result.liquid)>0)&&(gas==0||Exact::sign(result.gas)>0),
        "native signed total lost a positive measured phase");
    require((liquid!=0||Exact::sign(result.liquid)==0)&&(gas!=0||Exact::sign(result.gas)==0),
        "native signed total created an unsupported phase");
    return result;
}
Pair lowPair(const E& total,const E& alpha)
{auto l=Exact::multiply(total,alpha);return {l,Exact::subtract(total,l)};}
struct Face{Pair highP,highN,lowP,lowN;E correction;};
struct PositiveSum
{
    long double sum=0,error=0;
    void add(long double x){require(x>=0&&std::isfinite(x),"invalid positive inventory contribution");
        const auto y=x-error,t=sum+y;error=(t-sum)-y;sum=t;}
};
}
namespace
{
using RestrictionRational=Geometry::Detail::FaceTime::FaceTimeRational;
enum class RestrictionOperation {Add,Subtract,Multiply,Divide};
// Both complete homogeneous operands receive powers of two BEFORE products.
// Include all least bits, the leading product carry and input scaling bounds;
// leave two further leading bits for sums and signed differences. This does
// not extend the exact engine's exponent window or discard unsupported tails.
RestrictionRational restrictionBinary(const RestrictionRational& a,
    const RestrictionRational& b,RestrictionOperation operation)
{
    const int signA=Exact::sign(a.n),signB=Exact::sign(b.n);
    if(operation==RestrictionOperation::Divide)require(signB!=0,"zero restriction divisor");
    if(operation==RestrictionOperation::Multiply || operation==RestrictionOperation::Divide){
        if(signA==0 || (operation==RestrictionOperation::Multiply && signB==0))return RestrictionRational{};}
    else {
        if(signB==0)return a;
        if(signA==0)return operation==RestrictionOperation::Add?b:RestrictionRational(negative(b.n),b.d);}
    std::vector<std::pair<const E*,const E*>> products;
    if(operation==RestrictionOperation::Multiply)products={{&a.n,&b.n},{&a.d,&b.d}};
    else if(operation==RestrictionOperation::Divide)products={{&a.n,&b.d},{&a.d,&b.n}};
    else products={{&a.n,&b.d},{&b.n,&a.d},{&a.d,&b.d}};
    Exact::ExponentRange wholeA,wholeB;wholeA.include(a.n);wholeA.include(a.d);wholeB.include(b.n);wholeB.include(b.d);
    int lowest=std::numeric_limits<int>::max(),highest=std::numeric_limits<int>::min();
    for(const auto& p:products){Exact::ExponentRange x,y;x.include(*p.first);y.include(*p.second);
        if(x.lowest==std::numeric_limits<int>::max() || y.lowest==std::numeric_limits<int>::max())continue;
        lowest=std::min(lowest,x.lowest+y.lowest);highest=std::max(highest,x.highest+y.highest+1);}
    require(lowest!=std::numeric_limits<int>::max(),"empty restriction product plan");
    const int aLow=-1074-wholeA.lowest,aHigh=1023-wholeA.highest;
    const int bLow=-1074-wholeB.lowest,bHigh=1023-wholeB.highest;
    const int commonLow=std::max(-1074-lowest,aLow+bLow),commonHigh=std::min(1020-highest,aHigh+bHigh);
    require(commonLow<=commonHigh,"restriction exact products exceed supported exponent window");
    const int common=std::clamp(-(lowest+highest)/2,commonLow,commonHigh);
    const int leftLow=std::max(aLow,common-bHigh),leftHigh=std::min(aHigh,common-bLow);
    require(leftLow<=leftHigh,"restriction products have no lossless homogeneous scale");
    const int left=std::clamp(common/2,leftLow,leftHigh),right=common-left;
    const auto an=Exact::scalePow2(a.n,left),ad=Exact::scalePow2(a.d,left);
    const auto bn=Exact::scalePow2(b.n,right),bd=Exact::scalePow2(b.d,right);
    if(operation==RestrictionOperation::Multiply)return {Exact::multiply(an,bn),Exact::multiply(ad,bd)};
    if(operation==RestrictionOperation::Divide)return {Exact::multiply(an,bd),Exact::multiply(ad,bn)};
    const auto lhs=Exact::multiply(an,bd),rhs=Exact::multiply(bn,ad);
    return {operation==RestrictionOperation::Add?Exact::add(lhs,rhs):Exact::subtract(lhs,rhs),Exact::multiply(ad,bd)};
}
RestrictionRational restrictionAdd(const RestrictionRational& a,const RestrictionRational& b)
{return restrictionBinary(a,b,RestrictionOperation::Add);}
RestrictionRational restrictionSubtract(const RestrictionRational& a,const RestrictionRational& b)
{return restrictionBinary(a,b,RestrictionOperation::Subtract);}
RestrictionRational restrictionMultiply(const RestrictionRational& a,const RestrictionRational& b)
{return restrictionBinary(a,b,RestrictionOperation::Multiply);}
RestrictionRational restrictionDivide(const RestrictionRational& a,const RestrictionRational& b)
{return restrictionBinary(a,b,RestrictionOperation::Divide);}
RestrictionRational restrictionNumber(const E& value)
{return RestrictionRational(value,Exact::fromDouble(1));}
int restrictionCompare(const RestrictionRational& a,const RestrictionRational& b)
{
    const int x=Exact::sign(a.n),y=Exact::sign(b.n);
    if(x!=y)return x<y?-1:1;if(x==0)return 0;
    const auto positiveA=x<0?RestrictionRational(negative(a.n),a.d):a;
    const auto positiveB=y<0?RestrictionRational(negative(b.n),b.d):b;
    const auto comparison=Geometry::Detail::ExactRationalRound::compare(positiveA,positiveB);
    return x<0?-comparison:comparison;
}
bool restrictionRelative(const RestrictionRational& actual,
    const RestrictionRational& expected,const RestrictionRational& tolerance)
{
    require(Exact::sign(actual.n)>=0&&Exact::sign(expected.n)>=0,"negative restriction positive measure");
    if(Exact::sign(expected.n)==0)return Exact::sign(actual.n)==0;
    const auto difference=restrictionSubtract(actual,expected);
    const auto magnitude=Exact::sign(difference.n)<0?RestrictionRational(negative(difference.n),difference.d):difference;
    return restrictionCompare(restrictionDivide(magnitude,expected),tolerance)<=0;
}
Pair restrictionRawPair(long double liquid,long double gas)
{
    require(std::isfinite(liquid)&&std::isfinite(gas)&&liquid>=0&&gas>=0,"invalid raw restriction phase input");
    // These are two independent actual geometrical integrals. In particular,
    // do not apply highPair's dominant-complement reconstruction to raw A/a.
    return {scalar(liquid),scalar(gas)};
}
void restrictionRawClosure(const Pair& pair,const E& total,const RestrictionRational& tolerance)
{
    require(Exact::sign(total)>=0,"negative raw restriction gross total");
    const auto phases=restrictionAdd(restrictionNumber(pair.liquid),restrictionNumber(pair.gas));
    require(restrictionRelative(phases,restrictionNumber(total),tolerance),"raw restriction phases lost analytic gross coverage");
}
void restrictionParentCoefficient(const E& high,const E& raw,const RestrictionRational& tolerance)
{
    require(Exact::sign(high)>=0&&Exact::sign(raw)>=0,"negative parent restriction phase");
    if(Exact::sign(raw)==0){require(Exact::sign(high)==0,"parent high phase has absent raw support");return;}
    require(Exact::sign(high)>0,"positive raw parent phase lost numerical high phase");
    require(restrictionRelative(restrictionNumber(high),restrictionNumber(raw),tolerance),
        "parent phase calibration exceeds own-phase conservation certificate");
}
void restrictionSubset(const E& part,const E& parent,const RestrictionRational& tolerance)
{
    require(Exact::sign(part)>=0&&Exact::sign(parent)>=0,"negative raw restriction subset");
    if(Exact::sign(parent)==0){require(Exact::sign(part)==0,"raw restriction subset has absent parent support");return;}
    const auto a=restrictionNumber(part),b=restrictionNumber(parent);
    require(restrictionCompare(a,b)<=0 || restrictionRelative(a,b,tolerance),
        "raw restriction subset exceeds retained parent measure");
}
RestrictionRational restrictionHigh(const E& high,const E& parent,const E& part)
{
    if(Exact::sign(part)==0)return RestrictionRational{};
    require(Exact::sign(parent)>0&&Exact::sign(high)>0,"positive subset has no calibrated parent phase");
    return restrictionMultiply(restrictionNumber(high),restrictionDivide(restrictionNumber(part),restrictionNumber(parent)));
}
struct RestrictionBranch
{
    RestrictionRational highLiquid,highGas,total,liquid,gas;
};
RestrictionBranch restrictionBranch(const Pair& high,const Pair& parent,const Pair& part,
    const E& donorAlpha,const RestrictionRational& lambda,const RestrictionRational& complement)
{
    RestrictionBranch result;result.highLiquid=restrictionHigh(high.liquid,parent.liquid,part.liquid);
    result.highGas=restrictionHigh(high.gas,parent.gas,part.gas);
    result.total=restrictionAdd(result.highLiquid,result.highGas);
    const auto lowLiquid=restrictionMultiply(restrictionNumber(donorAlpha),result.total);
    const auto lowGas=restrictionMultiply(restrictionNumber(Exact::subtract(Exact::fromDouble(1),donorAlpha)),result.total);
    result.liquid=restrictionAdd(restrictionMultiply(complement,lowLiquid),restrictionMultiply(lambda,result.highLiquid));
    result.gas=restrictionAdd(restrictionMultiply(complement,lowGas),restrictionMultiply(lambda,result.highGas));
    return result;
}
RestrictionBranch nativeRestrictionBranch(const Pair& high,const Pair& parent,const Pair& part,
    const E& total,const E& donorAlpha,const RestrictionRational& lambda,
    const RestrictionRational& complement,const RestrictionRational& tolerance)
{
    RestrictionBranch result;result.highLiquid=restrictionHigh(high.liquid,parent.liquid,part.liquid);
    result.highGas=restrictionHigh(high.gas,parent.gas,part.gas);
    result.total=restrictionNumber(total);
    // Keep the independently measured restricted minor (after the retained
    // parent calibration). Restore its dominant complement to the actual
    // query gross, rather than letting two parent ratios change the mean.
    const bool liquidMinor=Exact::sign(Exact::subtract(part.liquid,part.gas))<=0;
    auto& dominant=liquidMinor?result.highGas:result.highLiquid;
    const auto& minor=liquidMinor?result.highLiquid:result.highGas;const auto original=dominant;
    dominant=restrictionSubtract(result.total,minor);
    require(Exact::sign(dominant.n)>=0,"native restriction cannot retain its measured minor phase");
    require((Exact::sign(original.n)==0)==(Exact::sign(dominant.n)==0),
        "native restriction changed phase support");
    require(restrictionRelative(dominant,original,tolerance),
        "native restriction phase reconciliation exceeds own-phase certificate");
    const auto lowLiquid=restrictionMultiply(restrictionNumber(donorAlpha),result.total);
    const auto lowGas=restrictionMultiply(restrictionNumber(Exact::subtract(Exact::fromDouble(1),donorAlpha)),result.total);
    result.liquid=restrictionAdd(restrictionMultiply(complement,lowLiquid),restrictionMultiply(lambda,result.highLiquid));
    result.gas=restrictionAdd(restrictionMultiply(complement,lowGas),restrictionMultiply(lambda,result.highGas));
    return result;
}
long double restrictionValue(const RestrictionRational& value)
{return quotient(value.n,value.d);}
}
struct GeometricVofPhaseBudgetRestriction::Data
{
    struct Row
    {
        Pair highP,highN,rawP,rawN;
        E rawTotalP,rawTotalN,donorAlphaP,donorAlphaN;
        Ratio lambda;
        long double limiterDiagnostic=1;
    };
    std::vector<Row> rows;
    RestrictionRational tolerance;
    bool preserveNativeSignedTotal=false;
};
GeometricVofIntegratedFlux reconcileGeometricVofNativeSignedTotal(const GeometricVofIntegratedFlux& raw)
{
    Exact::validateActualEnvironment();const auto gross=nativeGross(raw);
    const auto p=nativeHighPair(gross.positive,raw.positiveTotal,raw.positiveLiquid,raw.positiveGas);
    const auto n=nativeHighPair(gross.negative,raw.negativeTotal,raw.negativeLiquid,raw.negativeGas);
    auto result=raw;result.positiveTotal=value(gross.positive);result.negativeTotal=value(gross.negative);
    result.positiveLiquid=value(p.liquid);result.positiveGas=value(p.gas);
    result.negativeLiquid=value(n.liquid);result.negativeGas=value(n.gas);
    result.liquid=value(Exact::subtract(p.liquid,n.liquid));result.gas=value(Exact::subtract(p.gas,n.gas));
    result.mass=0;Exact::validateActualEnvironment();return result;
}
GeometricVofPhaseBudgetRestriction::GeometricVofPhaseBudgetRestriction(std::shared_ptr<const Data> data)
    :data_(std::move(data))
{require(static_cast<bool>(data_),"empty retained phase restriction");}
GeometricVofIntegratedFlux GeometricVofPhaseBudgetRestriction::measure(
    std::size_t row,const GeometricVofIntegratedFlux& raw) const
{
    Exact::validateActualEnvironment();require(static_cast<bool>(data_)&&row<data_->rows.size(),"invalid retained restriction row");
    require(std::isfinite(raw.total)&&std::isfinite(raw.positiveTotal)&&std::isfinite(raw.negativeTotal)&&
        raw.positiveTotal>=0&&raw.negativeTotal>=0,"invalid raw restriction gross total");
    const auto& source=data_->rows[row];const auto& tolerance=data_->tolerance;
    const auto rp=restrictionRawPair(raw.positiveLiquid,raw.positiveGas),rn=restrictionRawPair(raw.negativeLiquid,raw.negativeGas);
    const auto tp=scalar(raw.positiveTotal),tn=scalar(raw.negativeTotal);
    restrictionRawClosure(rp,tp,tolerance);restrictionRawClosure(rn,tn,tolerance);
    restrictionSubset(rp.liquid,source.rawP.liquid,tolerance);restrictionSubset(rp.gas,source.rawP.gas,tolerance);
    restrictionSubset(rn.liquid,source.rawN.liquid,tolerance);restrictionSubset(rn.gas,source.rawN.gas,tolerance);
    restrictionSubset(tp,source.rawTotalP,tolerance);restrictionSubset(tn,source.rawTotalN,tolerance);
    const RestrictionRational lambda(source.lambda.n,source.lambda.d);
    const RestrictionRational complement(Exact::subtract(source.lambda.d,source.lambda.n),source.lambda.d);
    const auto gross=data_->preserveNativeSignedTotal?nativeGross(raw):NativeGross{tp,tn};
    const auto positive=data_->preserveNativeSignedTotal?
        nativeRestrictionBranch(source.highP,source.rawP,rp,gross.positive,source.donorAlphaP,lambda,complement,tolerance):
        restrictionBranch(source.highP,source.rawP,rp,source.donorAlphaP,lambda,complement);
    const auto negativeBranch=data_->preserveNativeSignedTotal?
        nativeRestrictionBranch(source.highN,source.rawN,rn,gross.negative,source.donorAlphaN,lambda,complement,tolerance):
        restrictionBranch(source.highN,source.rawN,rn,source.donorAlphaN,lambda,complement);
    require(restrictionRelative(positive.total,restrictionNumber(tp),tolerance)&&
        restrictionRelative(negativeBranch.total,restrictionNumber(tn),tolerance),
        "calibrated restriction total exceeds raw analytic normal certificate");
    GeometricVofIntegratedFlux result;
    result.positiveTotal=restrictionValue(positive.total);result.negativeTotal=restrictionValue(negativeBranch.total);
    require(restrictionRelative(restrictionNumber(scalar(result.positiveTotal)),restrictionNumber(tp),tolerance)&&
        restrictionRelative(restrictionNumber(scalar(result.negativeTotal)),restrictionNumber(tn),tolerance),
        "rounded calibrated restriction total exceeds raw analytic normal certificate");
    result.total=restrictionValue(restrictionSubtract(positive.total,negativeBranch.total));
    result.rawLiquid=restrictionValue(restrictionSubtract(positive.highLiquid,negativeBranch.highLiquid));
    result.positiveLiquid=restrictionValue(positive.liquid);result.negativeLiquid=restrictionValue(negativeBranch.liquid);
    result.positiveGas=restrictionValue(positive.gas);result.negativeGas=restrictionValue(negativeBranch.gas);
    result.liquid=restrictionValue(restrictionSubtract(positive.liquid,negativeBranch.liquid));
    result.gas=restrictionValue(restrictionSubtract(positive.gas,negativeBranch.gas));
    result.mass=0;result.limiter=source.limiterDiagnostic;
    Exact::validateActualEnvironment();return result;
}
GeometricVofAcceptedWorkMoments GeometricVofPhaseBudgetRestriction::liquidWorkMoments(
    std::size_t row,const GeometricVofDirectionalWorkMoments& input,
    std::optional<std::array<std::array<long double,2>,2>> acceptedP0Phases) const
{
    Exact::validateActualEnvironment();
    require(static_cast<bool>(data_)&&row<data_->rows.size(),"invalid retained work-moment row");
    require(std::isfinite(input.mean)&&std::isfinite(input.halfVariation)&&
        std::isfinite(input.duration)&&input.duration>0,"invalid directional work trace/duration");
    const auto validMoments=[&](const std::array<long double,2>& moments,const char* reason){
        require(std::isfinite(moments[0])&&std::isfinite(moments[1])&&moments[0]>=0&&
            std::abs(moments[1])<=moments[0],reason);
    };
    validMoments(input.zeroTraceLiquid,"invalid zero-trace liquid support moments");
    for(unsigned branch=0;branch<2;++branch){
        validMoments(input.gross[branch],"invalid directional gross support moments");
        validMoments(input.liquid[branch],"invalid directional liquid support moments");
        require(input.liquid[branch][0]<=input.gross[branch][0],"raw liquid work support exceeds gross support");
        if(acceptedP0Phases)for(const auto phase:(*acceptedP0Phases)[branch])
            require(std::isfinite(phase)&&phase>=0,"invalid declared accepted P0 branch phase");
    }
    const auto& source=data_->rows[row];
    const std::array<const Pair*,2> high{{&source.highP,&source.highN}},rawPhase{{&source.rawP,&source.rawN}};
    const std::array<const E*,2> donorAlpha{{&source.donorAlphaP,&source.donorAlphaN}};
    const std::array<E,2> gross{{Exact::add(source.highP.liquid,source.highP.gas),
        Exact::add(source.highN.liquid,source.highN.gas)}};
    GeometricVofAcceptedWorkMoments result;
    if(input.mean==0&&input.halfVariation==0){
        for(unsigned branch=0;branch<2;++branch){
            require(Exact::sign(gross[branch])==0,"zero trace has nonzero retained gross extent");
            require(input.gross[branch][0]==0&&input.gross[branch][1]==0&&
                input.liquid[branch][0]==0&&input.liquid[branch][1]==0,"zero trace has active directional support");
            if(acceptedP0Phases)require((*acceptedP0Phases)[branch][0]==0&&(*acceptedP0Phases)[branch][1]==0,
                "zero trace has a nonzero accepted P0 extent");
        }
        // A zero normal trace still has a real sharp liquid face variation.
        // Keep it for the adjoint; its work and signed liquid extent are zero.
        result.value=input.zeroTraceLiquid;result.zeroTrace=true;
        Exact::validateActualEnvironment();return result;
    }
    const auto number=[](long double x){return restrictionNumber(scalar(x));};
    const auto magnitude=[](const RestrictionRational& x){
        return Exact::sign(x.n)<0?RestrictionRational(negative(x.n),x.d):x;
    };
    const auto mean=number(input.mean),variation=number(input.halfVariation),duration=number(input.duration);
    // Normalize the actual trace BEFORE products. These ratios retain its
    // original bits; neither tiny components nor a large velocity are squared
    // in native arithmetic. P1 physical trace L2 metric is diag(1,1/3).
    const auto scale=number(std::max(std::abs(input.mean),std::abs(input.halfVariation)));
    const std::array<RestrictionRational,2> normalized{{restrictionDivide(mean,scale),restrictionDivide(variation,scale)}};
    const auto third=restrictionDivide(number(1),number(3));
    const std::array<RestrictionRational,2> metricDirection{{normalized[0],restrictionMultiply(normalized[1],third)}};
    const auto normSquare=restrictionAdd(restrictionMultiply(normalized[0],metricDirection[0]),
        restrictionMultiply(normalized[1],metricDirection[1]));
    require(Exact::sign(normSquare.n)>0,"nonzero work trace has no positive physical L2 metric");
    const RestrictionRational lambda(source.lambda.n,source.lambda.d);
    const RestrictionRational complement(Exact::subtract(source.lambda.d,source.lambda.n),source.lambda.d);
    const auto dot=[&](const std::array<RestrictionRational,2>& moments){
        return restrictionAdd(restrictionMultiply(mean,moments[0]),restrictionMultiply(variation,moments[1]));
    };
    const auto recordRawResidual=[&](const RestrictionRational& measured,const RestrictionRational& expected,
        const RestrictionRational& zeroSupportScale){
        require(Exact::sign(expected.n)>=0&&Exact::sign(zeroSupportScale.n)>=0,"negative retained work diagnostic scale");
        const auto normalization=Exact::sign(expected.n)>0?expected:zeroSupportScale;
        if(Exact::sign(normalization.n)==0){require(Exact::sign(measured.n)==0,"inactive gross has directional work");return;}
        const auto relative=restrictionDivide(magnitude(restrictionSubtract(measured,expected)),normalization);
        // This is the UNCORRECTED midpoint work discrepancy, not a quadrature
        // certificate. For absent raw liquid its branch gross supplies a
        // finite diagnostic scale; no nonzero high contribution is invented.
        result.rawWorkRelativeResidual=std::max(result.rawWorkRelativeResidual,restrictionValue(relative));
    };
    std::array<std::array<RestrictionRational,2>,2> accepted{};
    std::array<RestrictionRational,2> acceptedLiquid{};
    for(unsigned branch=0;branch<2;++branch){
        const auto actualGross=restrictionNumber(gross[branch]);
        const std::array<RestrictionRational,2> grossMoments{{number(input.gross[branch][0]),number(input.gross[branch][1])}},
            liquidMoments{{number(input.liquid[branch][0]),number(input.liquid[branch][1])}};
        const auto signedDuration=branch==0?duration:restrictionSubtract(RestrictionRational{},duration);
        const auto grossWork=restrictionMultiply(signedDuration,dot(grossMoments));
        const auto rawLiquidWork=restrictionMultiply(signedDuration,dot(liquidMoments));
        recordRawResidual(grossWork,actualGross,RestrictionRational{});
        recordRawResidual(rawLiquidWork,restrictionNumber(rawPhase[branch]->liquid),actualGross);
        if(Exact::sign(gross[branch])==0){
            require(input.gross[branch][0]==0&&input.gross[branch][1]==0&&
                input.liquid[branch][0]==0&&input.liquid[branch][1]==0,"inactive branch has directional support");
            if(acceptedP0Phases)require((*acceptedP0Phases)[branch][0]==0&&(*acceptedP0Phases)[branch][1]==0,
                "inactive branch has declared accepted P0 phase");
            continue;
        }
        require(Exact::sign(grossWork.n)>0,"positive retained gross has no directional work support");
        const auto grossScale=restrictionDivide(actualGross,grossWork);
        std::array<RestrictionRational,2> calibratedGross{},calibratedHigh{};
        for(unsigned k=0;k<2;++k)calibratedGross[k]=restrictionMultiply(grossMoments[k],grossScale);
        if(Exact::sign(high[branch]->liquid)>0){
            const auto rawLiquid=restrictionNumber(rawPhase[branch]->liquid);
            require(Exact::sign(rawLiquid.n)>0,"positive retained high liquid has absent original raw phase");
            const auto signedRaw=branch==0?rawLiquid:restrictionSubtract(RestrictionRational{},rawLiquid);
            const auto target=restrictionDivide(signedRaw,restrictionMultiply(duration,scale));
            const auto midpointDot=restrictionAdd(restrictionMultiply(liquidMoments[0],normalized[0]),
                restrictionMultiply(liquidMoments[1],normalized[1]));
            const auto correction=restrictionDivide(restrictionSubtract(target,midpointDot),normSquare);
            // Gonzalez discrete gradient of the ORIGINAL sharp finite extent.
            // Retain midpoint shape information orthogonal to the trace; do
            // not replace this by an accepted/raw uniform face coefficient.
            std::array<RestrictionRational,2> rawGradient{};
            for(unsigned k=0;k<2;++k)rawGradient[k]=restrictionAdd(liquidMoments[k],
                restrictionMultiply(correction,metricDirection[k]));
            const auto highScale=restrictionDivide(restrictionNumber(high[branch]->liquid),rawLiquid);
            for(unsigned k=0;k<2;++k)calibratedHigh[k]=restrictionMultiply(rawGradient[k],highScale);
        }
        // A zero raw/high phase contributes no sharp high-order accepted
        // channel. Its midpoint discrepancy stays visible above; the original
        // donor-alpha low-order contribution may still be strictly positive.
        if(acceptedP0Phases){
            const auto& declared=(*acceptedP0Phases)[branch];
            const auto phases=nativeHighPair(gross[branch],restrictionValue(actualGross),declared[0],declared[1]);
            restrictionParentCoefficient(phases.liquid,scalar(declared[0]),data_->tolerance);
            restrictionParentCoefficient(phases.gas,scalar(declared[1]),data_->tolerance);
            acceptedLiquid[branch]=restrictionNumber(phases.liquid);
            const auto fraction=restrictionDivide(acceptedLiquid[branch],actualGross);
            for(unsigned k=0;k<2;++k)accepted[branch][k]=restrictionMultiply(fraction,calibratedGross[k]);
        }else{
            const auto lowFactor=restrictionMultiply(complement,restrictionNumber(*donorAlpha[branch]));
            for(unsigned k=0;k<2;++k)accepted[branch][k]=restrictionAdd(
                restrictionMultiply(lowFactor,calibratedGross[k]),restrictionMultiply(lambda,calibratedHigh[k]));
            acceptedLiquid[branch]=restrictionAdd(restrictionMultiply(lowFactor,actualGross),
                restrictionMultiply(lambda,restrictionNumber(high[branch]->liquid)));
        }
        for(unsigned k=0;k<2;++k)result.branch[branch][k]=restrictionValue(accepted[branch][k]);
    }
    for(unsigned k=0;k<2;++k)result.value[k]=restrictionValue(restrictionAdd(accepted[0][k],accepted[1][k]));
    // Sum the retained branch extents in exact arithmetic, independently of
    // rounded branch views and of the final rounded work coefficients.
    const auto exactSignedLiquid=restrictionSubtract(acceptedLiquid[0],acceptedLiquid[1]);
    const std::array<RestrictionRational,2> roundedMoments{{number(result.value[0]),number(result.value[1])}};
    const auto roundedWork=restrictionMultiply(duration,dot(roundedMoments));
    result.roundingWorkResidual=restrictionValue(restrictionSubtract(roundedWork,exactSignedLiquid));
    Exact::validateActualEnvironment();return result;
}
GeometricVofPhaseBudgetResult prepareGeometricVofPhaseBudget(
    const std::vector<GeometricVofPhaseBudgetCell>& cells,
    const std::vector<GeometricVofTransportFace>& faces,bool half,long double maximumOutflowCfl,
    const std::vector<GeometricVofTransportFace>* geometry,long double relativeConservationTolerance,
    const std::vector<double>* exteriorDonorAlpha,bool preserveNativeSignedTotal)
{
    Exact::validateActualEnvironment();
    require(!cells.empty(),"empty phase inventory domain");
    require(std::isfinite(maximumOutflowCfl)&&maximumOutflowCfl>0&&maximumOutflowCfl<=1,"invalid gross outflow CFL");
    require(std::isfinite(relativeConservationTolerance)&&relativeConservationTolerance>0&&
        relativeConservationTolerance<=1e-12L,"invalid restriction conservation certificate");
    require(geometry==nullptr || geometry->size()==faces.size(),"raw geometry restriction row coverage differs");
    require(exteriorDonorAlpha==nullptr || exteriorDonorAlpha->size()==faces.size(),"exterior donor row coverage differs");
    auto restrictionData=std::make_shared<GeometricVofPhaseBudgetRestriction::Data>();
    restrictionData->rows.resize(faces.size());restrictionData->tolerance=restrictionNumber(scalar(relativeConservationTolerance));
    restrictionData->preserveNativeSignedTotal=preserveNativeSignedTotal;
    const auto count=cells.size();const auto one=Exact::fromDouble(1);
    std::vector<E> capacity(count),alpha(count),outgoing(count),divergence(count),lowL(count),lowG(count),plus(count),minus(count);
    for(std::size_t c=0;c<count;++c){require(cells[c].capacity>0&&std::isfinite(cells[c].capacity)&&
        std::isfinite(cells[c].alpha)&&cells[c].alpha>=0&&cells[c].alpha<=1,"invalid source phase cell");
        capacity[c]=scalar(cells[c].capacity);
        if(cells[c].fraction){const auto& primary=*cells[c].fraction;const auto minor=primary.minorFraction();
            const bool liquidMinor=primary.minorPhase()==VofPhaseFraction::Phase::Liquid;
            require(cells[c].alpha==(liquidMinor?minor:1-minor),"typed fraction and rounded donor view disagree");
            // Same exact alpha as alphaComponents(), directly in this original
            // kernel's arithmetic. No new link dependency for legacy kernels.
            alpha[c]=liquidMinor?Exact::fromDouble(minor):Exact::subtract(one,Exact::fromDouble(minor));}
        else alpha[c]=Exact::fromDouble(cells[c].alpha);}
    std::vector<Face> work(faces.size());GeometricVofPhaseBudgetResult result;result.flux.reserve(faces.size());
    for(std::size_t row=0;row<faces.size();++row){const auto& face=faces[row];const auto& f=half?face.half:face.full;
        require(face.nativeRow==static_cast<int>(row)&&face.negativeCell!=face.positiveCell&&
            face.negativeCell>=-1&&face.positiveCell>=-1&&
            (face.negativeCell<0||static_cast<std::size_t>(face.negativeCell)<count)&&
            (face.positiveCell<0||static_cast<std::size_t>(face.positiveCell)<count),"invalid native face incidence");
        const bool exterior=face.negativeCell<0 || face.positiveCell<0;
        require(!exterior || exteriorDonorAlpha || (f.positiveTotal==0&&f.negativeTotal==0),"budget requires an explicit exterior donor input for an open boundary");
        E exteriorAlpha;
        if(exterior && exteriorDonorAlpha){const auto a=(*exteriorDonorAlpha)[row];
            require(std::isfinite(a)&&a>=0&&a<=1,"invalid explicit exterior donor alpha");exteriorAlpha=Exact::fromDouble(a);}
        const auto gross=preserveNativeSignedTotal?nativeGross(f):NativeGross{scalar(f.positiveTotal),scalar(f.negativeTotal)};
        const auto& tp=gross.positive;const auto& tn=gross.negative;const auto total=Exact::subtract(tp,tn);
        auto& w=work[row];
        w.highP=preserveNativeSignedTotal?nativeHighPair(tp,f.positiveTotal,f.positiveLiquid,f.positiveGas):highPair(f.positiveTotal,f.positiveLiquid,f.positiveGas);
        w.highN=preserveNativeSignedTotal?nativeHighPair(tn,f.negativeTotal,f.negativeLiquid,f.negativeGas):highPair(f.negativeTotal,f.negativeLiquid,f.negativeGas);
        w.lowP=lowPair(tp,face.negativeCell>=0?alpha[face.negativeCell]:exteriorAlpha);
        w.lowN=lowPair(tn,face.positiveCell>=0?alpha[face.positiveCell]:exteriorAlpha);
        if(face.negativeCell>=0){Exact::accumulate(outgoing[face.negativeCell],tp);Exact::accumulate(divergence[face.negativeCell],total);}
        if(face.positiveCell>=0){Exact::accumulate(outgoing[face.positiveCell],tn);Exact::accumulate(divergence[face.positiveCell],negative(total));}
        w.correction=Exact::subtract(Exact::subtract(w.highP.liquid,w.lowP.liquid),Exact::subtract(w.highN.liquid,w.lowN.liquid));
        result.flux.push_back(f);
        auto& retained=restrictionData->rows[row];retained.highP=w.highP;retained.highN=w.highN;
        retained.donorAlphaP=face.negativeCell>=0?alpha[face.negativeCell]:exteriorAlpha;
        retained.donorAlphaN=face.positiveCell>=0?alpha[face.positiveCell]:exteriorAlpha;
        const auto& rawFace=geometry?(*geometry)[row]:face;
        require(rawFace.nativeRow==face.nativeRow&&rawFace.negativeCell==face.negativeCell&&rawFace.positiveCell==face.positiveCell,
            "raw geometry restriction native incidence differs");
        const auto& raw=half?rawFace.half:rawFace.full;
        require(std::isfinite(raw.total)&&std::isfinite(raw.positiveTotal)&&std::isfinite(raw.negativeTotal)&&
            raw.total==f.total&&raw.positiveTotal==f.positiveTotal&&raw.negativeTotal==f.negativeTotal,
            "raw geometry restriction native total differs");
        retained.rawP=restrictionRawPair(raw.positiveLiquid,raw.positiveGas);retained.rawN=restrictionRawPair(raw.negativeLiquid,raw.negativeGas);
        retained.rawTotalP=scalar(raw.positiveTotal);retained.rawTotalN=scalar(raw.negativeTotal);
        restrictionRawClosure(retained.rawP,retained.rawTotalP,restrictionData->tolerance);
        restrictionRawClosure(retained.rawN,retained.rawTotalN,restrictionData->tolerance);
        restrictionParentCoefficient(retained.highP.liquid,retained.rawP.liquid,restrictionData->tolerance);
        restrictionParentCoefficient(retained.highP.gas,retained.rawP.gas,restrictionData->tolerance);
        restrictionParentCoefficient(retained.highN.liquid,retained.rawN.liquid,restrictionData->tolerance);
        restrictionParentCoefficient(retained.highN.gas,retained.rawN.gas,restrictionData->tolerance);
    }
    for(std::size_t c=0;c<count;++c){const auto cfl=value(outgoing[c])/cells[c].capacity;
        result.maximumOutflowCfl=std::max(result.maximumOutflowCfl,cfl);
        require(Exact::sign(Exact::subtract(outgoing[c],Exact::multiply(capacity[c],scalar(maximumOutflowCfl))))<=0,"gross outflow CFL exceeds configured bound cell="+std::to_string(c));
        const auto defect=std::abs(value(divergence[c]))/cells[c].capacity;
        result.maximumTotalDivergence=std::max(result.maximumTotalDivergence,defect);
        require(defect<=arithmetic,"native total flux is not discretely incompressible cell="+std::to_string(c));
        const auto remaining=Exact::subtract(capacity[c],outgoing[c]);require(Exact::sign(remaining)>=0,"negative remaining native capacity");
        lowL[c]=Exact::multiply(alpha[c],remaining);lowG[c]=Exact::multiply(Exact::subtract(one,alpha[c]),remaining);
    }
    for(std::size_t row=0;row<faces.size();++row){const auto& face=faces[row];const auto& w=work[row];
        if(face.negativeCell>=0){Exact::accumulate(lowL[face.negativeCell],w.lowN.liquid);Exact::accumulate(lowG[face.negativeCell],w.lowN.gas);}
        if(face.positiveCell>=0){Exact::accumulate(lowL[face.positiveCell],w.lowP.liquid);Exact::accumulate(lowG[face.positiveCell],w.lowP.gas);}
        const auto sign=Exact::sign(w.correction);const auto magnitude=sign<0?negative(w.correction):w.correction;
        if(face.negativeCell>=0)Exact::accumulate(sign>0?minus[face.negativeCell]:plus[face.negativeCell],magnitude);
        if(face.positiveCell>=0)Exact::accumulate(sign>0?plus[face.positiveCell]:minus[face.positiveCell],magnitude);
    }
    std::vector<Ratio> rPlus(count),rMinus(count);std::vector<PositiveSum> liquid(count),gas(count);
    for(std::size_t c=0;c<count;++c){rPlus[c]=limited(lowG[c],plus[c]);rMinus[c]=limited(lowL[c],minus[c]);
        const auto l=Exact::subtract(lowL[c],minus[c]),g=Exact::subtract(lowG[c],plus[c]);
        // Only saturated mathematical bases are zero. A sub-ulp unsaturated
        // remainder is still a nonzero expansion and survives conversion.
        if(Exact::sign(l)>=0)liquid[c].add(value(l));if(Exact::sign(g)>=0)gas[c].add(value(g));}
    for(std::size_t row=0;row<faces.size();++row){const auto& face=faces[row];const auto& w=work[row];auto& f=result.flux[row];
        Ratio lambda=unit();const auto sign=Exact::sign(w.correction);
        if(sign!=0){if(face.negativeCell>=0){const auto& r=sign>0?rMinus[face.negativeCell]:rPlus[face.negativeCell];if(less(r,lambda))lambda=r;}
            if(face.positiveCell>=0){const auto& r=sign>0?rPlus[face.positiveCell]:rMinus[face.positiveCell];if(less(r,lambda))lambda=r;}}
        f.limiter=quotient(lambda.n,lambda.d);require(std::isfinite(f.limiter)&&f.limiter>=0&&f.limiter<=1,"unrepresentable bounded limiter diagnostic");
        if(less(lambda,unit()))++result.limitedFaceExtents;
        const auto complement=Exact::subtract(lambda.d,lambda.n);
        const auto blend=[&](const E& low,const E& high){return Exact::add(Exact::multiply(complement,low),Exact::multiply(lambda.n,high));};
        const auto lp=blend(w.lowP.liquid,w.highP.liquid),ln=blend(w.lowN.liquid,w.highN.liquid);
        const auto gp=blend(w.lowP.gas,w.highP.gas),gn=blend(w.lowN.gas,w.highN.gas);
        f.positiveLiquid=quotient(lp,lambda.d);f.negativeLiquid=quotient(ln,lambda.d);
        f.positiveGas=quotient(gp,lambda.d);f.negativeGas=quotient(gn,lambda.d);
        f.liquid=quotient(Exact::subtract(lp,ln),lambda.d);f.gas=quotient(Exact::subtract(gp,gn),lambda.d);
        // The interval totals and raw candidate remain the original native
        // records. Densities are applied by Transport to these accepted signs.
        f.mass=0;
        restrictionData->rows[row].lambda=lambda;restrictionData->rows[row].limiterDiagnostic=f.limiter;
        if(sign!=0){const auto magnitude=sign<0?negative(w.correction):w.correction;
            for(auto incidence:{std::pair<int,int>{face.negativeCell,-sign},std::pair<int,int>{face.positiveCell,sign}})
                if(incidence.first>=0){const auto c=incidence.first;const auto& r=incidence.second>0?rPlus[c]:rMinus[c];
                    const auto retained=quotient(Exact::multiply(slackNumerator(r,lambda),magnitude),Exact::multiply(r.d,lambda.d));
                    const auto accepted=quotient(Exact::multiply(lambda.n,magnitude),lambda.d);
                    if(incidence.second>0){liquid[c].add(accepted);gas[c].add(retained);}
                    else {liquid[c].add(retained);gas[c].add(accepted);}}}
    }
    result.phaseInventory.resize(count);
    for(std::size_t c=0;c<count;++c){const auto l=liquid[c].sum,g=gas[c].sum;
        require(std::isfinite(l)&&std::isfinite(g)&&l>=0&&g>=0&&l+g>0,"exact phase budget lost bounds cell="+std::to_string(c));
        const auto defect=std::abs(l+g-cells[c].capacity)/cells[c].capacity;
        result.maximumCapacityReconciliation=std::max(result.maximumCapacityReconciliation,defect);
        require(defect<=arithmetic,"phase capacity reconciliation exceeds arithmetic budget");result.phaseInventory[c]={{l,g}};}
    result.restriction=std::shared_ptr<const GeometricVofPhaseBudgetRestriction>(
        new GeometricVofPhaseBudgetRestriction(std::move(restrictionData)));
    Exact::validateActualEnvironment();return result;
}
}
