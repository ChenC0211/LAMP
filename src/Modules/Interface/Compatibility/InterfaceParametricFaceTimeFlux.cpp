#include "Modules/Interface/Compatibility/InterfaceParametricTransport.h"
#include "Modules/Interface/Compatibility/Detail/InterfaceParametricWorldSheet.h"
#include <numeric>
#include <sstream>
#include <tuple>
#include <exception>
#include <optional>
#include <memory>

namespace InterfaceParametricTransport {
namespace {
using namespace Quantum::Geometry::Parametric;
void require(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
using namespace WorldSheet;
struct Piece {
    std::function<Crossing(double)> curve;
    std::function<GraphEnvelope(double,double)> envelope;
    double a=0,b=1;int strip=0,depth=0;double value=0,error=0,vertexError=0;
};
// One geometric subdivision tree per grid plane, shared by all its scalar
// face functionals. Children have exactly the original restriction path;
// equal parameter bounds from different paths are never merged by a key.
struct SheetGeometry {
    Sheet sheet;SheetBounds bounds;
    std::array<std::unique_ptr<SheetGeometry>,2> children;
    std::optional<std::vector<Chart>> charts;
    SheetGeometry(const Sheet& s,const SheetBounds& b):sheet(s),bounds(b){}
};
using PlaneGeometry=std::vector<std::unique_ptr<SheetGeometry>>;
// Positive flux is an integral functional, not a demand to represent every
// zero in global coordinates. Shared Bernstein bounds retain sign changes
// and charge unresolved sign regions to the existing face-volume budget.
class PositivePrimitive {
    using Term=Quantum::Geometry::UnitPolynomialTerm;
    std::array<std::vector<Term>,2> first_,last_;
    std::array<bool,2> frozen_{};
    double lower_,halfWidth_;
    mutable long double integralError_=0;
    std::function<double(double)> weight_;
    bool polynomialWeight_=false;
    int productIndex_=-1;
    bool positive_=true;
    std::vector<Term> terms(int half,double theta) const {
        if(frozen_[half])return first_[half];
        auto out=first_[half];for(auto& term:out)term.weight*=1-theta;
        for(auto term:last_[half]){term.weight*=theta;out.push_back(term);}return out;
    }
    static long double magnitude(const std::vector<Term>& terms) {
        long double bound=0;for(const auto& term:terms)for(double c:term.coefficients)bound+=std::abs(term.weight*term.factor*c);return bound;
    }
    Quantum::Geometry::UnitPolynomialIntegral weighted(int half,const std::vector<Term>& p,const std::vector<Term>& q,double lo,double hi,
        long double maximumError=std::numeric_limits<long double>::infinity()) const{
        using Integral=Quantum::Geometry::UnitPolynomialIntegral;const auto& rule=gaussEight();
        const bool frozen=p.size()==q.size()&&std::equal(p.begin(),p.end(),q.begin(),[](const Term& a,const Term& b){return a.coefficients==b.coefficients&&a.factor==b.factor&&a.weight==b.weight;});
        auto averageTerms=p;for(auto& term:averageTerms)term.weight*=.5;
        for(auto term:q){term.weight*=.5;averageTerms.push_back(term);}
        const auto exact=[&](const std::vector<Term>& terms,double a,double b){return Quantum::Geometry::integrateUnitPolynomialBernstein(terms,FaceVolumes::productDegree,productIndex_,half,a,b);};
        if(polynomialWeight_&&!positive_)return exact(averageTerms,lo,hi);
        const auto constant=[](const std::vector<Term>& terms){for(const auto& term:terms)if(term.factor!=0&&term.weight!=0)for(int k=1;k<int(term.coefficients.size());k++)if(term.coefficients[k]!=0)return false;return true;};
        const auto positiveAverage=[](long double a,long double b){if(a<=0&&b<=0)return 0.L;if(a>=0&&b>=0)return .5L*(a+b);const auto v=std::max(a,b);return .5L*v*(v/std::abs(b-a));};
        if(polynomialWeight_&&positive_&&constant(p)&&constant(q)){Term term;term.coefficients[0]=1;
            term.factor=positiveAverage(Quantum::Geometry::evaluateUnitPolynomial(p,0),Quantum::Geometry::evaluateUnitPolynomial(q,0));return exact({term},lo,hi);}
        const auto value=[&](double x){
            const long double first=Quantum::Geometry::evaluateUnitPolynomial(p,x),last=frozen?first:Quantum::Geometry::evaluateUnitPolynomial(q,x);
            long double average=.5L*(first+last);
            if(positive_){
                if(first<=0&&last<=0)average=0;
                else if(first>0&&last<0)average=.5L*first*(first/(first-last));
                else if(last>0&&first<0)average=.5L*last*(last/(last-first));
            }
            const double w=weight_(lower_+halfWidth_*(half+x));require(w>=0&&w<=1&&std::isfinite(w),"face functional weight outside its declared unit bound");
            return average*w;
        };
        const auto integrate=[&](double a,double b){Integral out;
            for(int i=0;i<8;i++){const long double v=.5L*(b-a)*rule.w[i]*value(a+.5*(b-a)*(1+rule.x[i]));out.value+=v;out.absoluteError+=std::abs(v);}
            out.absoluteError*=4*std::numeric_limits<double>::epsilon();return out;};
        std::function<Integral(double,double,long double,int)> adaptive;
        adaptive=[&](double a,double b,long double budget,int depth){const double mid=.5*(a+b);bool quadrature=true;
            if(positive_){const auto first=Quantum::Geometry::encloseUnitPolynomial(p,a,b),last=frozen?first:Quantum::Geometry::encloseUnitPolynomial(q,a,b);
                if(first[1]<=0&&last[1]<=0)return Integral{};
                if(first[0]>=0&&last[0]>=0&&polynomialWeight_)return exact(averageTerms,a,b);
                const bool samePositive=first[0]>=0&&last[0]>=0;
                const bool opposite=(first[0]>=0&&last[1]<=0)||(first[1]<=0&&last[0]>=0);
                if(!samePositive&&!opposite){
                    // The multiplier is nonnegative and <=1. Bound unresolved
                    // sign regions in integral units, so a narrow counterflow
                    // interval cannot be missed by point quadrature.
                    const long double upper=(b-a)*positiveAverage(first[1],last[1]);
                    if(.5L*upper<=budget)return Integral{.5L*upper,.5L*upper};quadrature=false;
                }
            }
            if(quadrature){const auto coarse=integrate(a,b);auto fine=integrate(a,mid);const auto right=integrate(mid,b);
                fine.value+=right.value;fine.absoluteError+=right.absoluteError+std::abs(fine.value-coarse.value);
                if(fine.absoluteError<=budget)return fine;}
            require(depth<64&&mid>a&&mid<b,"weighted face primitive quadrature unresolved");
            auto left=adaptive(a,mid,.5L*budget,depth+1);const auto other=adaptive(mid,b,.5L*budget,depth+1);left.value+=other.value;left.absoluteError+=other.absoluteError;return left;};
        return adaptive(lo,hi,std::min(maximumError,32*std::numeric_limits<double>::epsilon()*(hi-lo)*std::max(magnitude(p),magnitude(q))),0);
    }
public:
    PositivePrimitive(const TimeLinearFlux& native,bool horizontal,int line,int strip,bool positive=true)
        :lower_(native.begin.face(horizontal,line,strip).lower),halfWidth_(.5*(native.begin.face(horizontal,line,strip).upper-lower_)),positive_(positive) {
        require(bool(native.begin.functionalWeight)==bool(native.end.functionalWeight),"face functional changed across time");
        productIndex_=native.begin.productFunctional;require(productIndex_==native.end.productFunctional,"face product basis changed across time");
        polynomialWeight_=productIndex_>=0&&!(horizontal&&native.begin.grid.axisymmetric&&lower_>0);
        if(native.begin.functionalWeight)weight_=[w=native.begin.functionalWeight,horizontal,line,strip](double x){return w(horizontal,line,strip,x);};
        for(int half=0;half<2;++half){first_[half]=native.begin.polynomial(horizontal,line,strip,half);
            last_[half]=native.end.polynomial(horizontal,line,strip,half);
            frozen_[half]=first_[half].size()==last_[half].size() && std::equal(first_[half].begin(),first_[half].end(),last_[half].begin(),
                [](const Term& a,const Term& b){return a.coefficients==b.coefficients && a.factor==b.factor && a.weight==b.weight;});}
    }
    double spatial(double a,double b,double theta,double errorWeight=1) const {
        if(a==b)return 0;if(a>b)return -spatial(b,a,theta,errorWeight);long double result=0;
        for(int half=0;half<2;++half){const double origin=lower_+half*halfWidth_,lo=std::clamp((a-origin)/halfWidth_,0.,1.),hi=std::clamp((b-origin)/halfWidth_,0.,1.);
            if(lo==hi)continue;const auto p=terms(half,theta);
            const auto integral=weight_?weighted(half,p,p,lo,hi):Quantum::Geometry::integratePositiveUnitPolynomial(p,lo,hi,32*std::numeric_limits<double>::epsilon()*(hi-lo)*magnitude(p));
            result+=integral.value;integralError_+=std::abs(errorWeight)*integral.absoluteError;}
        return static_cast<double>(result);
    }
    double time(double a,double b,double begin,double end,long double maximumError=std::numeric_limits<long double>::infinity()) const {
        if(a==b || begin==end)return 0;if(a>b)return -time(b,a,begin,end,maximumError);if(begin>end)return -time(a,b,end,begin,maximumError);
        long double result=0;
        for(int half=0;half<2;++half){const double origin=lower_+half*halfWidth_,lo=std::clamp((a-origin)/halfWidth_,0.,1.),hi=std::clamp((b-origin)/halfWidth_,0.,1.);
            if(lo==hi)continue;const auto p=terms(half,begin),q=terms(half,end);
            const long double limit=.5L*maximumError/(end-begin);
            const auto tolerance=std::min(limit,32*std::numeric_limits<double>::epsilon()*(hi-lo)*std::max(magnitude(p),magnitude(q)));
            const auto integral=weight_?weighted(half,p,q,lo,hi,limit):(frozen_[half]?Quantum::Geometry::integratePositiveUnitPolynomial(p,lo,hi,tolerance):
                Quantum::Geometry::integratePositiveTimeLinearUnitPolynomial(p,q,lo,hi,tolerance));
            result+=(end-begin)*integral.value;integralError_+=(end-begin)*integral.absoluteError;
        }
        return static_cast<double>(result);
    }
    long double integralError() const noexcept {return integralError_;}
    std::array<Quantum::Geometry::UnitPolynomialIntegral,FaceVolumes::productDegree+1> polynomialSignedBasis() const {
        // One frozen or time-linear normal trace supplies ALL work weights.
        // The caller only uses this polynomial primitive where the physical
        // weight has no off-axis RZ denominator. Sign-split flux is separate.
        require(!positive_,"signed basis integral requested from a positive primitive");
        std::array<Quantum::Geometry::UnitPolynomialIntegral,FaceVolumes::productDegree+1> result{};
        for(int half=0;half<2;half++){
            const auto p=terms(half,0),q=terms(half,1);auto average=p;for(auto& term:average)term.weight*=.5;
            for(auto term:q){term.weight*=.5;average.push_back(term);}
            const auto values=Quantum::Geometry::integrateUnitPolynomialBernsteinBasis(average,FaceVolumes::productDegree,half,0,1);
            for(int k=0;k<=FaceVolumes::productDegree;k++){result[k].value+=values[k].value;result[k].absoluteError+=values[k].absoluteError;}
        }
        return result;
    }
};
std::vector<double> planeFlux(const std::vector<Sheet>& sheets,const std::vector<SheetBounds>& sheetBounds,Sampler& sampler,PlaneGeometry& geometry,const TimeLinearFlux& native,
    double dt,bool horizontal,int line,FluxReport& report,bool positivePart=false,double budgetScale=1) {
    const int n=native.begin.stripCount(horizontal,line),axis=horizontal?1:0;const double plane=native.begin.planeCoordinate(horizontal,line);
    const auto loFace=[&](int strip){return native.begin.face(horizontal,line,strip).lower;};
    const auto hiFace=[&](int strip){return native.begin.face(horizontal,line,strip).upper;};
    const double minimumCapacity=native.begin.planeCapacity(horizontal,line),budget=budgetScale*1e-13*minimumCapacity;
    std::vector<double> errorWeight(n);double maximumCapacity=0;
    for(int strip=0;strip<n;strip++){const double capacity=native.begin.faceCapacity(horizontal,line,strip);
        errorWeight[strip]=minimumCapacity/capacity;maximumCapacity=std::max(maximumCapacity,capacity);}
    std::vector<PositivePrimitive> positive;const bool usePrimitive=positivePart||bool(native.begin.functionalWeight);
    if(usePrimitive)for(int strip=0;strip<n;++strip)positive.emplace_back(native,horizontal,line,strip,positivePart);
    const auto spatialIntegral=[&](int strip,double a,double b,double theta,double errorWeight=1){return usePrimitive?positive[strip].spatial(a,b,theta,errorWeight):native.spatialIntegral(horizontal,line,strip,a,b,theta);};
    const auto timeIntegral=[&](int strip,double a,double b,double begin,double end,long double allowance=std::numeric_limits<long double>::infinity()){
        return usePrimitive?positive[strip].time(a,b,begin,end,allowance):native.timeIntegral(horizontal,line,strip,a,b,begin,end);};
    // Direct wet-boundary and uncut-face rectangles have disjoint transverse
    // support. Reserve half the common plane budget for their primitives,
    // distributed by face count and occupied width, in the same V_min/V_f
    // norm used for final acceptance. The remaining budget covers crossings.
    const auto rectangleAllowance=[&](int strip,double lo,double hi){
        return .5L*budget*(hi-lo)/(hiFace(strip)-loFace(strip))/(n*std::abs(dt)*errorWeight[strip]);};
    const auto normal=[axis](Point p){return axis?p.y:p.x;};const auto transverse=[axis](Point p){return axis?p.x:p.y;};
    std::vector<double> flux(n);
    const auto primitiveError=[&](int strip=-1){long double sum=0;if(usePrimitive){
        if(strip>=0)return positive[strip].integralError()*errorWeight[strip];for(int k=0;k<n;k++)sum+=positive[k].integralError()*errorWeight[k];}return sum;};
    long double acceptedPrimitiveError=0;
    const auto finish=[&](){
        report.quadratureError+=static_cast<double>(acceptedPrimitiveError);
        if(!(report.quadratureError<=budget)){std::ostringstream reason;reason<<std::setprecision(17)
            <<"positive primitives exceed their shared plane integral budget: horizontal="<<horizontal<<" plane="<<plane<<" functional="<<native.begin.productFunctional
            <<" positive="<<positivePart<<" primitive_error="<<double(acceptedPrimitiveError)<<" total_error="<<report.quadratureError<<" budget="<<budget;
            throw std::runtime_error(reason.str());}
        report.quadratureFractionError=report.quadratureError/minimumCapacity;
        // During integration use min(V)/V_f weights, i.e. the sum of local
        // fraction errors. Convert diagnostics back to a conservative bound
        // in physical volume; acceptance uses the fractional norm itself.
        const double physicalBound=maximumCapacity/minimumCapacity;
        report.quadratureError*=physicalBound;report.contactIntegralError*=physicalBound;report.vertexIntervalError*=physicalBound;
        return flux;};
    const bool exterior=native.begin.exteriorPlane(horizontal,line);
    if(exterior&&!horizontal&&plane==0&&native.begin.grid.axisymmetric){
        for(int strip=0;strip<n;strip++)require(native.begin.densityMagnitudeBound(horizontal,line,strip)==0&&native.end.densityMagnitudeBound(horizontal,line,strip)==0,"RZ axis has nonzero supplied phase flux");
        return finish();
    }
    const auto closure=sampler.old.closure();
    const bool fixedWetTrace=exterior&&std::any_of(closure.begin(),closure.end(),[&](const auto& edge){
        return normal(edge.p[0])==plane&&normal(edge.p[0]+edge.p[1])==plane;});
    if(fixedWetTrace){
        // Only a coincident fixed domain closure owns a stationary wet
        // trace (e.g. a pinned inlet disk). An unrelated axis closure must
        // not suppress free-surface world-sheet crossings at an outlet.
        for(const auto& edge:closure){
            const Point a=edge.p[0],b=a+edge.p[1];if(normal(a)!=plane||normal(b)!=plane)continue;
            const double first=transverse(a),last=transverse(b);
            const double orientation=horizontal?(line==0?1.:-1.):(line==0?-1.:1.);
            require(orientation*(last-first)>0,"domain phase closure has reversed wet boundary orientation");
            for(int strip=0;strip<n;strip++){
                const double lo=std::max(loFace(strip),std::min(first,last)),hi=std::min(hiFace(strip),std::max(first,last));if(!(lo<hi))continue;
                const auto before=primitiveError(strip);flux[strip]+=dt*timeIntegral(strip,lo,hi,0,1,rectangleAllowance(strip,lo,hi));
                acceptedPrimitiveError+=std::abs(dt)*(primitiveError(strip)-before);
            }
        }
        return finish();
    }
    std::vector<bool> cut(n,false);
    for(std::size_t i=0;i<sheets.size();++i){const auto& sheet=sheets[i];const auto& bounds=sheetBounds[i];const auto f=axis?bounds.y:bounds.x,z=axis?bounds.x:bounds.y;const double pad=std::max(1e-12,sheet.positionEnclosure);
        if(f.lo>plane+pad || f.hi<plane-pad)continue;
        for(int strip=0;strip<n;++strip)if(z.lo<=hiFace(strip)+pad && z.hi>=loFace(strip)-pad)cut[strip]=true;}
    // A face-time rectangle with no interface crossing has constant phase.
    // Evaluate that integral directly before quadrature; this is neither an
    // inventory clamp nor a correction after cancellation of large primitives.
    PhaseMeasure oldPhase(sampler.old);
    // A structurally zero trace (in particular the RZ axis) carries no
    // volume for either phase, even if a physical pole lies on the face.
    for(int strip=0;strip<n;strip++)if(native.begin.densityMagnitudeBound(horizontal,line,strip)==0&&native.end.densityMagnitudeBound(horizontal,line,strip)==0)cut[strip]=false;
    if(positivePart)for(int strip=0;strip<n;++strip)
        if(timeIntegral(strip,loFace(strip),hiFace(strip),0,1)==0)cut[strip]=false;
    const bool fromLower=exterior&&line+1==native.begin.planeCount(horizontal);
    for(int strip=0;strip<n;++strip)if(!cut[strip]){
        const double coordinate=.5*(loFace(strip)+hiFace(strip));bool liquid=false;
        for(const auto& interval:oldPhase.intervals(axis,plane,fromLower))if(coordinate>interval[0] && coordinate<interval[1]){liquid=true;break;}
        const auto before=primitiveError(strip);
        flux[strip]=liquid?dt*timeIntegral(strip,loFace(strip),hiFace(strip),0,1,rectangleAllowance(strip,loFace(strip),hiFace(strip))):0;
        acceptedPrimitiveError+=std::abs(dt)*(primitiveError(strip)-before);
    }
    if(std::none_of(cut.begin(),cut.end(),[](bool value){return value;}))return finish();
    std::vector<Piece> pieces;
    std::vector<double> densityBounds(n);
    const bool hasPolynomial=native.begin.polynomial && native.end.polynomial;
    using HalfDensities=std::array<std::vector<Quantum::Geometry::UnitPolynomialTerm>,2>;
    std::vector<std::array<HalfDensities,2>> densityPolynomials(hasPolynomial?n:0);
    for(int strip=0;strip<n;++strip){densityBounds[strip]=std::max(native.begin.densityMagnitudeBound(horizontal,line,strip),native.end.densityMagnitudeBound(horizontal,line,strip));
        require(densityBounds[strip]>=0 && std::isfinite(densityBounds[strip]),"invalid native density bound");
        if(hasPolynomial)for(int half=0;half<2;++half){
            densityPolynomials[strip][0][half]=native.begin.polynomial(horizontal,line,strip,half);
            densityPolynomials[strip][1][half]=native.end.polynomial(horizontal,line,strip,half);}}
    const auto densityBound=[&](int strip,double a,double b){
        if(!(a<b))return 0.;if(!hasPolynomial)return densityBounds[strip];double bound=0;
        for(int half=0;half<2;++half){const double scale=2/(hiFace(strip)-loFace(strip)),origin=loFace(strip)+.5*half*(hiFace(strip)-loFace(strip));
            const double lo=std::clamp((a-origin)*scale,0.,1.),hi=std::clamp((b-origin)*scale,0.,1.);
            if(lo==hi)continue;
            const double lower=std::max(0.,std::nextafter(lo,-std::numeric_limits<double>::infinity()));
            const double upper=std::min(1.,std::nextafter(hi,std::numeric_limits<double>::infinity()));
            for(int endpoint=0;endpoint<2;++endpoint){
                const auto range=Quantum::Geometry::encloseUnitPolynomial(densityPolynomials[strip][endpoint][half],lower,upper);
                bound=std::max(bound,scale*std::max(std::abs(range[0]),std::abs(range[1])));}}
        // Convexity covers every time-linear trace, including its positive
        // part. Bound the retained density only on the uncertain wet interval;
        // preserve cancellation before the outward Bernstein conversion.
        require(std::isfinite(bound),"retained density bound overflow");
        if(native.begin.productFunctional>=0){const double a0=std::clamp((a-loFace(strip))/(hiFace(strip)-loFace(strip)),0.,1.),b0=std::clamp((b-loFace(strip))/(hiFace(strip)-loFace(strip)),0.,1.);
            bound*=FaceVolumes::productWeightBound(native.begin.productFunctional,a0,b0,horizontal,native.begin.grid.axisymmetric,loFace(strip),hiFace(strip));}
        return bound==0?0.:std::nextafter(bound,std::numeric_limits<double>::infinity());
    };
    struct Occupation {long double zeroth=0,first=0;double uncertainty=0;std::vector<std::array<double,3>> intervals;};
    const auto occupation=[&](const Sheet& sheet,double side,double tolerance) {
        // Moments of 1_{F(side,theta)<0} dtheta. The actual face-relative
        // map determines the sign; the retained derivative enclosure bounds
        // a whole interval, including even-multiplicity endpoint roots.
        // Unknown intervals carry half their measure as an error, rather
        // than being classified by a single sample or silently discarded.
        Occupation result;std::vector<Sheet> unknown;
        const auto add=[&](const Sheet& part,double weight) {
            const long double width=static_cast<long double>(part.d)-part.c;
            result.zeroth+=weight*width;result.first+=weight*width*(static_cast<long double>(part.c)+part.d)*.5L;
            if(positivePart)result.intervals.push_back({part.c,part.d,weight});
        };
        const auto classify=[&](const Sheet& part) {
            // Rebuild the scalar edge jet on THIS time interval. Carrying
            // the full 2D sheet's fixed vector interpolation error down to
            // a stationary endpoint would never resolve its O(time^2)
            // normal motion. Retain face-relative values and the derivatives
            // of the same map, without subtracting absolute coordinates.
            const auto start=sampler.at(part,side,0,axis,plane),finish=sampler.at(part,side,1,axis,plane);
            const double left=start.f,right=finish.f,d0=normal(start.dt),d1=normal(finish.dt),delta=right-left;
            const double p2=3*delta-2*d0-d1,p3=d0+d1-2*delta;
            Bounds derivative{std::min({d0,3*delta-d0-d1,d1}),std::max({d0,3*delta-d0-d1,d1})};
            double derivativeError=0,positionError=0,middle=0;
            for(double t:{.2113248654051871,.5,.7886751345948129}) {
                const auto exact=sampler.at(part,side,t,axis,plane);
                positionError=std::max(positionError,std::abs(left+t*(d0+t*(p2+t*p3))-exact.f));
                derivativeError=std::max(derivativeError,std::abs(d0+t*(2*p2+t*3*p3)-normal(exact.dt)));
                if(t==.5)middle=exact.f;
            }
            // As for the retained world sheet, these are interpolation-error
            // estimates, not an interval certificate for arbitrary callbacks.
            derivative.lo-=4*derivativeError;derivative.hi+=4*derivativeError;
            const double magnitude=std::max(std::abs(derivative.lo),std::abs(derivative.hi));
            const double roundoff=128*std::numeric_limits<double>::epsilon()*std::max({std::abs(left),std::abs(middle),std::abs(right),magnitude});
            const double lo=std::max({middle-.5*magnitude,left+std::min(0.,derivative.lo),right-std::max(0.,derivative.hi),
                std::min({left,left+d0/3,right-d1/3,right})-4*positionError})-roundoff;
            const double hi=std::min({middle+.5*magnitude,left+std::max(0.,derivative.hi),right-std::min(0.,derivative.lo),
                std::max({left,left+d0/3,right-d1/3,right})+4*positionError})+roundoff;
            if(lo>=0)return 1;if(hi<0)return -1;return 0;
        };
        const auto insert=[&](const Sheet& part) {const int phase=classify(part);if(phase<0)add(part,1);else if(phase==0)unknown.push_back(part);};
        insert(sheet);
        for(int iteration=0;!unknown.empty();++iteration) {
            double width=0;std::size_t worst=0;
            for(std::size_t k=0;k<unknown.size();++k){width+=unknown[k].d-unknown[k].c;
                if(unknown[k].d-unknown[k].c>unknown[worst].d-unknown[worst].c)worst=k;}
            if(.5*width<=tolerance){result.uncertainty=.5*width;for(const auto& part:unknown)add(part,.5);break;}
            if(iteration>=4096){const auto& unresolved=unknown[worst];
                std::ostringstream message;message<<std::setprecision(17)<<"nonregular contact edge occupation unresolved: axis="<<axis<<" line="<<line
                    <<" material="<<sheet.a<<':'<<sheet.b<<" side="<<side<<" time="<<sheet.c<<':'<<sheet.d
                    <<" unresolved_measure="<<width<<" tolerance="<<tolerance<<" intervals="<<unknown.size()
                    <<" worst="<<unresolved.c<<':'<<unresolved.d<<" f="<<sampler.at(unresolved,side,0,axis,plane).f<<':'
                    <<sampler.at(unresolved,side,.5,axis,plane).f<<':'<<sampler.at(unresolved,side,1,axis,plane).f;throw std::runtime_error(message.str());}
            const auto part=unknown[worst];unknown[worst]=unknown.back();unknown.pop_back();
            const double middle=.5*(part.c+part.d);
            require(middle>part.c && middle<part.d,"nonregular contact edge exhausted time resolution");
            insert(part.restrict(0,1,0,.5));insert(part.restrict(0,1,.5,1));++report.contactEdgeRefinements;
        }
        return result;
    };
    const auto contact=[&](const Sheet& sheet,const SheetBounds& bounds) {
        const auto errorBefore=primitiveError();
        // On Omega={F<0}, Green's theorem gives the SAME oriented face
        // contour integral as the regular charts:
        //   sigma * [ integral_edges 1_Omega P dt - integral_Omega q z_s ds dt ].
        // Write P=P(anchor,t)+(P-P(anchor,t)). The first term uses only the
        // two material-edge occupation moments. Bound both residual edge
        // terms and the area term, with no regularity assumption on F=0.
        // In particular an isolated touch, a birth, or a returning branch
        // has its constant-primitive contribution; none is simply dropped.
        const double allowance=(.25*budget-report.contactIntegralError)*.125;
        if(!(allowance>0))return false;
        auto z=axis?bounds.x:bounds.y;const auto zs=axis?bounds.xs:bounds.ys;
        const double pad=sheet.positionEnclosure+64*std::numeric_limits<double>::epsilon();z.lo-=pad;z.hi+=pad;
        const double variation=std::max(std::abs(zs.lo),std::abs(zs.hi)),timeWidth=sheet.d-sheet.c;
        std::vector<double> first(n),last(n),anchors(n);double residualBound=0,primitiveBound=0;
        const double center=transverse(sampler.at(sheet,.5,.5).p);
        for(int strip=0;strip<n;++strip)if(cut[strip]) {
            const double lower=loFace(strip),upper=hiFace(strip),lo=std::clamp(z.lo,lower,upper),hi=std::clamp(z.hi,lower,upper);
            const double anchor=std::clamp(center,lo,hi);
            anchors[strip]=anchor;
            const auto boundErrorBefore=primitiveError(strip);
            first[strip]=spatialIntegral(strip,lower,anchor,0);
            last[strip]=spatialIntegral(strip,lower,anchor,1);
            primitiveBound+=std::max(std::abs(first[strip]),std::abs(last[strip]))*errorWeight[strip]+static_cast<double>(primitiveError(strip)-boundErrorBefore);
            if(lo<hi)residualBound+=std::abs(dt)*densityBound(strip,lo,hi)*timeWidth*(2*std::max(anchor-lo,hi-anchor)+variation)*errorWeight[strip];
        }
        if(residualBound>.5*allowance)return false;
        const double coefficient=std::abs(dt)*primitiveBound;
        const double tolerance=coefficient>0?(allowance-residualBound)/(4*coefficient):timeWidth;
        const auto left=occupation(sheet,0,tolerance),right=occupation(sheet,1,tolerance);
        const double error=residualBound+coefficient*(left.uncertainty+right.uncertainty);
        require(error<=allowance,"nonregular contact exceeded its assigned flux budget");
        const long double m0=right.zeroth-left.zeroth,m1=right.first-left.first;
        for(int strip=0;strip<n;++strip)if(cut[strip]) {
            long double contribution=first[strip]*m0+(last[strip]-first[strip])*m1;
            if(positivePart){contribution=0;const double lower=loFace(strip);
                for(const auto& part:right.intervals)contribution+=part[2]*timeIntegral(strip,lower,anchors[strip],part[0],part[1]);
                for(const auto& part:left.intervals)contribution-=part[2]*timeIntegral(strip,lower,anchors[strip],part[0],part[1]);}
            flux[strip]+=static_cast<double>((horizontal?-1.L:1.L)*dt*contribution);
        }
        acceptedPrimitiveError+=std::abs(dt)*(primitiveError()-errorBefore);
        report.contactIntegralError+=error;++report.boundedContactPatches;return true;
    };
    const auto rule=[&](const Piece& piece,double a,double b,double anchor) {
        const auto errorBefore=primitiveError(piece.strip);
        const auto first=piece.curve(a),last=piece.curve(b);const double lower=loFace(piece.strip),upper=hiFace(piece.strip);
        const double base=dt*first.orientation*timeIntegral(piece.strip,lower,anchor,first.time,last.time);
        double value=0;for(int k=0;k<5;++k){const auto q=piece.curve(.5*(a+b+(b-a)*gx[k]));
            const double primitive=spatialIntegral(piece.strip,anchor,std::clamp(q.z,lower,upper),q.time,.5*(b-a)*gw[k]*q.dtime);
            value+=gw[k]*dt*q.orientation*primitive*q.dtime;}
        return Quantum::Geometry::UnitPolynomialIntegral{base+.5*(b-a)*value,std::abs(dt)*(primitiveError(piece.strip)-errorBefore)};
    };
    const auto evaluate=[&](Piece piece){
        const double lower=loFace(piece.strip),upper=hiFace(piece.strip),middle=.5*(piece.a+piece.b);
        const auto range=piece.envelope(piece.a,piece.b);const double lo=std::clamp(range.z.lo,lower,upper),hi=std::clamp(range.z.hi,lower,upper);
        require(lo<=hi && range.timeVariation>=0,"invalid regular-graph enclosure");
        if(lo==hi){const auto first=piece.curve(piece.a),last=piece.curve(piece.b);
            const auto before=primitiveError(piece.strip);
            piece.value=dt*first.orientation*timeIntegral(piece.strip,lower,lo,first.time,last.time);
            piece.error=static_cast<double>(std::abs(dt)*(primitiveError(piece.strip)-before));piece.vertexError=0;return piece;}
        const double anchor=std::clamp(piece.curve(middle).z,lo,hi);const auto coarse=rule(piece,piece.a,piece.b,anchor);
        const auto left=rule(piece,piece.a,middle,anchor),right=rule(piece,middle,piece.b,anchor);
        const double ruleErrors=static_cast<double>(coarse.absoluteError+2*(left.absoluteError+right.absoluteError));
        piece.value=static_cast<double>(left.value+right.value);piece.error=std::abs(piece.value-static_cast<double>(coarse.value))*errorWeight[piece.strip]+ruleErrors;piece.vertexError=0;
        // A root need not be simple or even change sign. Until the graph's
        // whole range lies in one half face, bound the omitted-resolution
        // uncertainty in flux units. Anchor subtraction integrates the large
        // constant primitive exactly. Both the exact residual integral and
        // its positive-weight quadrature are <= Q*(hi-lo)*totalVariation(t).
        if(std::floor((range.z.lo-lower)*2/(upper-lower))!=std::floor((range.z.hi-lower)*2/(upper-lower))) {
            piece.vertexError=2*std::abs(dt)*densityBound(piece.strip,lo,hi)*(hi-lo)*range.timeVariation*errorWeight[piece.strip];
            piece.error=std::max(piece.error,piece.vertexError+ruleErrors);
        }
        return piece;
    };
    std::function<void(SheetGeometry&,int)> append;
    append=[&](SheetGeometry& node,int depth) {
        const auto& sheet=node.sheet;const auto& bounds=node.bounds;
        const auto f=axis?bounds.y:bounds.x,fs=axis?bounds.ys:bounds.xs,ft=axis?bounds.yt:bounds.xt;
        const double enclosure=std::max(1e-12,sheet.positionEnclosure);
        if(f.lo>plane+enclosure || f.hi<plane-enclosure)return;
        const bool sGraph=strict(fs),tGraph=strict(ft);
        if(!sGraph && !tGraph) {
            if(contact(sheet,bounds))return;
            if(depth>=64){std::ostringstream message;message<<std::setprecision(17)<<"space-time chart unresolved: axis="<<axis<<" plane="<<plane<<" material="<<sheet.a<<':'<<sheet.b<<" time="<<sheet.c<<':'<<sheet.d;
                throw std::runtime_error(message.str());}
            ++report.chartSubdivisions;
            const bool splitS=std::max(std::abs(fs.lo),std::abs(fs.hi))>=std::max(std::abs(ft.lo),std::abs(ft.hi));
            const double middle=splitS?.5*(sheet.a+sheet.b):.5*(sheet.c+sheet.d);
            require(splitS?(middle>sheet.a && middle<sheet.b):(middle>sheet.c && middle<sheet.d),"space-time chart exhausted parameter resolution");
            if(!node.children[0]){
                const auto first=sheet.restrict(0,splitS?.5:1,0,splitS?1:.5);
                node.children[0]=std::make_unique<SheetGeometry>(first,enclose(first));
                const auto last=sheet.restrict(splitS?.5:0,1,splitS?0:.5,1);
                node.children[1]=std::make_unique<SheetGeometry>(last,enclose(last));
            }
            append(*node.children[0],depth+1);append(*node.children[1],depth+1);return;
        }
        // The chart's callbacks outlive this functional's temporary report.
        // Bind them to the plane-owned sampler report, shared for the whole
        // batch. Integral errors and budgets remain local to each query.
        if(!node.charts)node.charts=regularCharts(sheet,bounds,sampler,axis,plane,sampler.report);
        for(const auto& chart:*node.charts) {
            const auto& crossing=chart.curve;const auto& envelope=chart.envelope;
            const double a=chart.a,b=chart.b;const auto range=envelope(a,b);const auto first=crossing(a),last=crossing(b);
            // One clipped primitive per strip works for crossings, returns
            // and endpoint tangencies alike. No transverse monotonicity is
            // assumed, and an unresolved interval is never silently dropped.
            for(int strip=0;strip<n;++strip)if(cut[strip] && range.z.hi>loFace(strip)) {
                if(range.z.lo>=hiFace(strip)){const auto before=primitiveError(strip);
                    flux[strip]+=dt*first.orientation*timeIntegral(strip,loFace(strip),hiFace(strip),first.time,last.time);
                    acceptedPrimitiveError+=std::abs(dt)*(primitiveError(strip)-before);}
                else pieces.push_back(evaluate({crossing,envelope,a,b,strip}));}
        }
    };
    for(std::size_t i=0;i<sheets.size();++i){
        const auto f=axis?sheetBounds[i].y:sheetBounds[i].x;
        const double enclosure=std::max(1e-12,sheets[i].positionEnclosure);
        if(f.lo>plane+enclosure || f.hi<plane-enclosure)continue;
        if(!geometry[i])geometry[i]=std::make_unique<SheetGeometry>(sheets[i],sheetBounds[i]);
        append(*geometry[i],0);
    }
    for(int iteration=0;!pieces.empty();++iteration) {
        double error=report.contactIntegralError;std::size_t worst=0;for(std::size_t i=0;i<pieces.size();++i){error+=pieces[i].error;if(pieces[i].error>pieces[worst].error)worst=i;}
        if(error+acceptedPrimitiveError<=budget){report.quadratureError=std::max(report.quadratureError,error);break;}
        if(iteration>=100000){const auto& part=pieces[worst];std::ostringstream message;message<<std::setprecision(17)
            <<"space-time flux exhausted the plane quadrature budget: axis="<<axis<<" line="<<line<<" dt="<<dt
            <<" positive="<<positivePart<<" pieces="<<pieces.size()<<" error="<<error<<" budget="<<budget
            <<" contact_error="<<report.contactIntegralError<<" worst="<<part.a<<':'<<part.b
            <<" primitive_error="<<double(acceptedPrimitiveError)<<" functional="<<native.begin.productFunctional
            <<" piece_error="<<part.error<<" vertex_error="<<part.vertexError;
            const auto range=part.envelope(part.a,part.b);
            message<<" strip="<<part.strip<<" transverse_enclosure="<<range.z.lo<<':'<<range.z.hi
                <<" time_variation="<<range.timeVariation<<" density_bound="
                <<densityBound(part.strip,std::clamp(range.z.lo,loFace(part.strip),hiFace(part.strip)),
                    std::clamp(range.z.hi,loFace(part.strip),hiFace(part.strip)));
            message<<" chart=";for(double value:range.chart)message<<value<<',';
            for(double q:{part.a,.5*(part.a+part.b),part.b}){const auto point=part.curve(q);
                message<<" graph_sample="<<q<<','<<point.z<<','<<point.dz<<','<<point.time<<','<<point.dependent;}
            throw std::runtime_error(message.str());}
        const auto old=pieces[worst];const double middle=.5*(old.a+old.b);
        if(!(old.depth<52 && middle>old.a && middle<old.b)) {
            std::ostringstream message;message<<std::setprecision(17)<<"space-time flux quadrature unresolved: axis="<<axis<<" line="<<line
                <<" strip="<<old.strip<<" interval="<<old.a<<':'<<old.b<<" error="<<error<<" piece_error="<<old.error<<" budget="<<budget;
            throw std::runtime_error(message.str());
        }
        auto first=old,last=old;first.b=last.a=middle;++first.depth;++last.depth;
        pieces[worst]=evaluate(first);pieces.push_back(evaluate(last));++report.quadratureRefinements;
    }
    for(const auto& piece:pieces){flux[piece.strip]+=piece.value;
        if(piece.vertexError>0){++report.boundedVertexIntervals;report.vertexIntervalError+=piece.vertexError;}}
    report.quadratureError=std::max(report.quadratureError,report.contactIntegralError);
    return finish();
}
NativeFlux weightedTrace(const NativeFlux& native,int degree) {
    require(native.polynomial && degree>=1 && degree<=2,"face moments require the retained polynomial trace");
    using Terms=std::vector<Quantum::Geometry::UnitPolynomialTerm>;
    NativeFlux layout;layout.grid=native.grid;
    auto saved=std::make_shared<std::vector<std::array<Terms,2>>>(layout.faceCount(false)+layout.faceCount(true));
    const auto index=[layout](bool h,int line,int strip){return (h?layout.faceCount(false):0)+layout.faceOffset(h,line)+strip;};
    for(bool horizontal:{false,true})for(int line=0;line<layout.planeCount(horizontal);line++)for(int strip=0;strip<layout.stripCount(horizontal,line);strip++)for(int half=0;half<2;half++){
        const std::array<double,3> weight=degree==1?std::array<double,3>{.5*half,.5,0}:std::array<double,3>{.25*half,.5*half,.25};
        auto& out=(*saved)[index(horizontal,line,strip)][half];
        for(const auto& original:native.polynomial(horizontal,line,strip,half))for(int power=0;power<=degree;power++)if(weight[power]!=0){
            auto term=original;term.coefficients.fill(0);term.factor*=weight[power];
            require(original.factor==0 || term.factor!=0,"weighted face-moment scale is unrepresentable");
            for(int k=0;k<int(original.coefficients.size());k++)if(original.coefficients[k]!=0){
                require(k+power<int(term.coefficients.size()),"weighted face trace exceeds retained polynomial degree");term.coefficients[k+power]=original.coefficients[k];}
            out.push_back(term);
        }
    }
    NativeFlux out=layout;out.densityMagnitudeBound=native.densityMagnitudeBound; // 0 <= xi^k <= 1.
    out.polynomial=[saved,index](bool h,int line,int strip,int half){return saved->at(index(h,line,strip)).at(half);};
    out.integral=[saved,index,layout](bool h,int line,int strip,double a,double b){long double sum=0;const auto face=layout.face(h,line,strip);
        for(int half=0;half<2;half++){const double origin=face.lower+.5*half*(face.upper-face.lower),scale=2/(face.upper-face.lower);
            sum+=Quantum::Geometry::integrateUnitPolynomial(saved->at(index(h,line,strip))[half],0,1,
                std::clamp((a-origin)*scale,0.,1.),std::clamp((b-origin)*scale,0.,1.));}
        return static_cast<double>(sum);};
    out.density=[saved,index,layout](bool h,int line,int strip,double s){const auto face=layout.face(h,line,strip);const int half=s<=.5*(face.lower+face.upper)?0:1;const double scale=2/(face.upper-face.lower);
        return static_cast<double>(scale*Quantum::Geometry::evaluateUnitPolynomial(saved->at(index(h,line,strip))[half],std::clamp(scale*(s-face.lower)-half,0.,1.)));};
    return out;
}
NativeFlux productTrace(const NativeFlux& native,int basis){
    NativeFlux result=native;result.retainedVelocityProducts=false;result.productFunctional=basis;
    result.functionalWeight=[layout=native,basis](bool h,int line,int strip,double coordinate){const auto f=layout.face(h,line,strip);
        const double xi=std::clamp((coordinate-f.lower)/(f.upper-f.lower),0.,1.);
        return FaceVolumes::productWeight(basis,xi,h,layout.grid.axisymmetric,f.lower,f.upper);};
    // Keep the exact base q polynomial for sign decisions and its rigorous
    // envelope. PositivePrimitive integrates q times the bounded multiplier.
    return result;
}
}

std::vector<double> advanceInventory(const Curve& old,const Source& end,const Flow& flow,const NativeFlux& native,
    double dt,const std::vector<double>& previous,FluxReport& report) {
    return advanceInventory(old,end,flow,TimeLinearFlux{native,native},dt,previous,report);
}
namespace {
std::vector<double> integratePhaseFlux(const Curve* old,const Source* end,const Flow& flow,const TimeLinearFlux& native,
    double dt,const std::vector<double>& previous,FluxReport& report,FaceVolumes* retained) {
    const auto& grid=native.begin.grid;
    require(grid.sameLayout(native.end.grid),"face flux physical grid or native identity changes across time");
    require(grid.nx>0 && grid.ny>0 && grid.width>0 && grid.height>0 && previous.size()==native.begin.cellCount() && bool(native.begin.integral) && bool(native.end.integral)
        && bool(native.begin.density) && bool(native.end.density) && bool(native.begin.densityMagnitudeBound) && bool(native.end.densityMagnitudeBound)
        && (!old||(end&&bool(flow))) && std::isfinite(dt)&&dt!=0,"invalid parametric space-time flux input");
    const int momentDegree=retained?native.begin.retainedMomentDegree:0;
    require(native.begin.retainedVelocityProducts==native.end.retainedVelocityProducts,"complete face functional requirements changed across time");
    const bool products=retained&&native.begin.retainedVelocityProducts;
    require(!products||(native.begin.polynomial&&native.end.polynomial),"complete face products require the retained normal trace");
    require(native.begin.retainedMomentDegree==native.end.retainedMomentDegree && (momentDegree==0 || momentDegree==2),
        "face moment requirements changed across the physical interval");
    std::vector<TimeLinearFlux> momentFlux;
    for(int k=1;k<=momentDegree;k++)momentFlux.push_back({weightedTrace(native.begin,k),weightedTrace(native.end,k)});
    std::optional<Sampler> sampler;std::vector<Sheet> sheets;
    if(old){sampler.emplace(Sampler{*old,flow,dt,report});sheets=worldSheet(*old,*end,grid,*sampler,report);}
    // These polynomial bounds depend only on the retained sheet, not on the
    // queried grid plane. Reuse exactly the same bounds on every face; newly
    // subdivided charts still compute their own bounds before being used.
    std::vector<SheetBounds> bounds;bounds.reserve(sheets.size());for(const auto& sheet:sheets)bounds.push_back(enclose(sheet));
    const bool directions=retained && bool(native.begin.polynomial) && bool(native.end.polynomial);
    struct Plane {std::vector<double> flux;std::vector<FaceVolumes::Face> volumes;FluxReport report;std::exception_ptr failure;};
    const int verticalCount=native.begin.planeCount(false);std::vector<Plane> planes(verticalCount+native.begin.planeCount(true));
#pragma omp parallel for schedule(dynamic,1) if(planes.size()>16)
    for(int index=0;index<int(planes.size());++index){auto& plane=planes[index];try{
        std::optional<Sampler> local;if(old)local.emplace(Sampler{*old,flow,dt,plane.report,{},&sampler->cache});
        const bool horizontal=index>=verticalCount;const int line=index-(horizontal?verticalCount:0),strips=native.begin.stripCount(horizontal,line);
        PlaneGeometry geometry(sheets.size());
        const auto wetFlux=[&](const TimeLinearFlux& trace,FluxReport& measured,bool positive=false,double share=1.){
            return local?planeFlux(sheets,bounds,*local,geometry,trace,dt,horizontal,line,measured,positive,share):std::vector<double>(strips,0.);
        };
        plane.flux=wetFlux(native,plane.report,false,directions?.25:1);
        if(retained){plane.volumes.resize(strips);FluxReport positiveReport;std::vector<double> positive;
            if(directions)positive=wetFlux(native,positiveReport,true,.25);
            double netAdjustment=0,positiveAdjustment=0,netFraction=0,positiveFraction=0;
            for(int strip=0;strip<strips;++strip){auto& face=plane.volumes[strip];const auto chart=native.begin.face(horizontal,line,strip);const double lo=chart.lower,hi=chart.upper;
                face.total=dt*native.timeIntegral(horizontal,line,strip,lo,hi,0,1);face.liquid=plane.flux[strip];
                if(directions){const PositivePrimitive primitive(native,horizontal,line,strip);
                    const double forward=std::abs(dt)*primitive.time(lo,hi,0,1);
                    // This is a whole-face TOTAL-flow integral, not a cut
                    // liquid contribution to the common plane quadrature.
                    // Charge its error to the neighboring cell capacity of
                    // this face. Distant large gas faces must not exhaust the
                    // smallest axis cell's liquid-interface error allowance.
                    const double fractionError=std::abs(dt)*static_cast<double>(primitive.integralError())/native.begin.faceCapacity(horizontal,line,strip);
                    require(fractionError<=1e-13,"total directional face integral exceeds its local volume budget");
                    plane.report.totalFluxFractionError=std::max(plane.report.totalFluxFractionError,fractionError);
                    double tp=dt>=0?forward:face.total+forward;
                    const double floor=std::max(face.total,0.);
                    require(tp>=floor-128*std::numeric_limits<double>::epsilon()*std::max(std::abs(face.total),forward),
                        "directional total primitive disagrees with the native signed integral");
                    tp=std::max(tp,floor);face.totalPositive=tp;face.totalNegative=tp-face.total;
                    face.total=face.totalPositive-face.totalNegative;
                    // The known phase-subvolume bounds constrain the numerical
                    // quadrature itself, BEFORE this one shared flux updates
                    // either cell. Spend only the original integral budget;
                    // never clip an inventory or repair a geometric volume.
                    const double net=std::clamp(face.liquid,-face.totalNegative,face.totalPositive);
                    const double rawPositive=dt>=0?positive[strip]:net-positive[strip];
                    face.liquidPositive=std::clamp(rawPositive,std::max(net,0.),std::min(face.totalPositive,face.totalNegative+net));
                    face.liquidNegative=face.liquidPositive-net;face.liquid=face.liquidPositive-face.liquidNegative;
                    netAdjustment+=std::abs(face.liquid-plane.flux[strip]);
                    positiveAdjustment+=std::abs((dt>=0?face.liquidPositive:-face.liquidNegative)-positive[strip]);
                    const double capacity=native.begin.faceCapacity(horizontal,line,strip);
                    netFraction+=std::abs(face.liquid-plane.flux[strip])/capacity;
                    positiveFraction+=std::abs((dt>=0?face.liquidPositive:-face.liquidNegative)-positive[strip])/capacity;
                    plane.flux[strip]=face.liquid;
                }
            }
            const double budget=1e-13;
            if(directions && (plane.report.quadratureFractionError+netFraction>budget || positiveReport.quadratureFractionError+positiveFraction>budget)){
                std::ostringstream message;message<<std::setprecision(17)<<"directional phase subvolume bounds exceed the common face integral budget: axis="<<horizontal<<" plane="<<native.begin.planeCoordinate(horizontal,line)
                    <<" net_fraction="<<netFraction<<" positive_fraction="<<positiveFraction<<" signed_fraction="<<plane.report.quadratureFractionError<<" positive_fraction_error="<<positiveReport.quadratureFractionError<<" budget="<<budget;
                throw std::runtime_error(message.str());}
            plane.report.quadratureError+=netAdjustment;
            plane.report.quadratureFractionError+=netFraction;
            plane.report.directionalQuadratureError=positiveReport.quadratureError+positiveAdjustment;
            plane.report.directionalFractionError=positiveReport.quadratureFractionError+positiveFraction;
            plane.report.directionalBoundAdjustment=std::max(netAdjustment,positiveAdjustment);
            for(int power=1;power<=momentDegree;power++){
                FluxReport momentReport;const auto& weighted=momentFlux[power-1];
                const auto liquid=wetFlux(weighted,momentReport);
                plane.report.momentQuadratureError=std::max(plane.report.momentQuadratureError,momentReport.quadratureError);
                const bool exterior=directions&&native.begin.exteriorPlane(horizontal,line);
                FluxReport directionReport;std::vector<double> liquidPositive;
                if(exterior)liquidPositive=wetFlux(weighted,directionReport,true);
                plane.report.momentQuadratureError=std::max(plane.report.momentQuadratureError,directionReport.quadratureError);
                for(int strip=0;strip<strips;strip++){
                    const auto face=native.begin.face(horizontal,line,strip);
                    auto& saved=plane.volumes[strip];saved.totalHigher[power-1]=dt*weighted.timeIntegral(horizontal,line,strip,face.lower,face.upper,0,1);
                    saved.liquidHigher[power-1]=liquid[strip];
                    if(exterior){const PositivePrimitive primitive(weighted,horizontal,line,strip);
                        const double forward=std::abs(dt)*primitive.time(face.lower,face.upper,0,1);
                        const double error=std::abs(dt)*double(primitive.integralError())/native.begin.faceCapacity(horizontal,line,strip);
                        require(error<=1e-13,"directional face moment exceeds its local volume budget");
                        saved.totalPositiveHigher[power-1]=dt>=0?forward:saved.totalHigher[power-1]+forward;
                        saved.liquidPositiveHigher[power-1]=dt>=0?liquidPositive[strip]:liquid[strip]-liquidPositive[strip];
                        plane.report.totalFluxFractionError=std::max(plane.report.totalFluxFractionError,error);
                    }
                }
            }
            using ProductIntegrals=std::array<Quantum::Geometry::UnitPolynomialIntegral,FaceVolumes::productDegree+1>;
            std::vector<std::optional<ProductIntegrals>> wholeProducts(products?strips:0);
            if(products)for(int strip=0;strip<strips;strip++){
                const auto f=native.begin.face(horizontal,line,strip);
                if(!(horizontal&&grid.axisymmetric&&f.lower>0))
                    wholeProducts[strip]=PositivePrimitive(native,horizontal,line,strip,false).polynomialSignedBasis();
            }
            if(products)for(int k=0;k<=FaceVolumes::productDegree;k++){
                const TimeLinearFlux weighted{productTrace(native.begin,k),productTrace(native.end,k)};
                constexpr double share=1./(FaceVolumes::productDegree+1);FluxReport signedReport,positiveWeightedReport;
                const auto liquid=wetFlux(weighted,signedReport,false,share);
                const bool exterior=native.begin.exteriorPlane(horizontal,line);
                std::vector<double> liquidPositive;if(exterior)liquidPositive=wetFlux(weighted,positiveWeightedReport,true,share);
                plane.report.momentQuadratureError=std::max({plane.report.momentQuadratureError,signedReport.quadratureError,positiveWeightedReport.quadratureError});
                for(int strip=0;strip<strips;strip++){
                    const auto f=native.begin.face(horizontal,line,strip);auto& saved=plane.volumes[strip];
                    // Give the primitive the SAME physical error allowance
                    // checked below. A power-coefficient magnitude can be a
                    // very loose scale for cancelling high-order traces; it
                    // must not decide the accuracy of the final face flux.
                    // Exterior signed/positive queries split one allowance,
                    // and each query splits its allowance over two halves.
                    const long double allowance=share*1e-13L*native.begin.faceCapacity(horizontal,line,strip)/std::abs(dt);
                    const long double partAllowance=allowance/(exterior?2:1);
                    double integral,integralError;
                    if(wholeProducts[strip]){const auto& value=(*wholeProducts[strip])[k];integral=static_cast<double>(value.value);integralError=static_cast<double>(value.absoluteError);}
                    else{const PositivePrimitive whole(weighted,horizontal,line,strip,false);integral=whole.time(f.lower,f.upper,0,1,partAllowance);integralError=static_cast<double>(whole.integralError());}
                    saved.totalProducts[k]=dt*integral;saved.liquidProducts[k]=liquid[strip];
                    double error=std::abs(dt)*integralError/native.begin.faceCapacity(horizontal,line,strip);
                    if(exterior){const PositivePrimitive forward(weighted,horizontal,line,strip);
                        const double positive=std::abs(dt)*forward.time(f.lower,f.upper,0,1,partAllowance);
                        saved.totalPositiveProducts[k]=dt>=0?positive:saved.totalProducts[k]+positive;
                        saved.liquidPositiveProducts[k]=dt>=0?liquidPositive[strip]:liquid[strip]-liquidPositive[strip];
                        error+=std::abs(dt)*double(forward.integralError())/native.begin.faceCapacity(horizontal,line,strip);
                    }
                    if(!(error<=share*1e-13)){std::ostringstream reason;reason<<std::setprecision(17)
                        <<"complete face functional exceeds its local primitive budget: horizontal="<<horizontal<<" plane="<<native.begin.planeCoordinate(horizontal,line)
                        <<" strip="<<strip<<" basis="<<k<<" error="<<error<<" budget="<<share*1e-13;throw std::runtime_error(reason.str());}
                    plane.report.totalFluxFractionError=std::max(plane.report.totalFluxFractionError,error);
                }
            }
        }
    }catch(...){plane.failure=std::current_exception();}}
    for(const auto& plane:planes)if(plane.failure)std::rethrow_exception(plane.failure);
    // Only independent plane calculations run concurrently. Flux additions
    // to each inventory and the exterior balance retain the original order.
    auto inventory=previous;
    FaceVolumes published;published.hasDirections=directions;published.momentDegree=momentDegree;published.hasExteriorDirectionMoments=directions&&momentDegree==2;
    published.hasVelocityProducts=products;
    if(retained&&sampler)published.surface=measureSurfaceFlux(sheets,bounds,*sampler,grid);
    if(retained)for(bool h:{false,true})published.faces[h].resize(native.begin.faceCount(h));
    for(bool horizontal:{false,true})for(int line=0;line<native.begin.planeCount(horizontal);++line) {
        const auto& plane=planes[(horizontal?verticalCount:0)+line];const auto& flux=plane.flux;const auto& r=plane.report;
        report.crossingIntervals+=r.crossingIntervals;report.vertexEvents+=r.vertexEvents;report.tangentTurns+=r.tangentTurns;
        report.quadratureRefinements+=r.quadratureRefinements;report.timeCharts+=r.timeCharts;report.materialCharts+=r.materialCharts;report.chartSubdivisions+=r.chartSubdivisions;
        report.boundedVertexIntervals+=r.boundedVertexIntervals;report.vertexIntervalError=std::max(report.vertexIntervalError,r.vertexIntervalError);
        report.boundedContactPatches+=r.boundedContactPatches;report.contactEdgeRefinements+=r.contactEdgeRefinements;
        report.contactIntegralError=std::max(report.contactIntegralError,r.contactIntegralError);
        report.mapQueries+=r.mapQueries;report.mapCacheHits+=r.mapCacheHits;
        report.crossingResidual=std::max(report.crossingResidual,r.crossingResidual);report.quadratureError=std::max(report.quadratureError,r.quadratureError);
        report.directionalQuadratureError=std::max(report.directionalQuadratureError,r.directionalQuadratureError);
        report.directionalBoundAdjustment=std::max(report.directionalBoundAdjustment,r.directionalBoundAdjustment);
        report.momentQuadratureError=std::max(report.momentQuadratureError,r.momentQuadratureError);
        report.totalFluxFractionError=std::max(report.totalFluxFractionError,r.totalFluxFractionError);
        report.quadratureFractionError=std::max(report.quadratureFractionError,r.quadratureFractionError);
        report.directionalFractionError=std::max(report.directionalFractionError,r.directionalFractionError);
        for(int strip=0;strip<native.begin.stripCount(horizontal,line);++strip){const double q=flux[strip];
            if(retained)published.faces[horizontal][native.begin.faceOffset(horizontal,line)+strip]=plane.volumes[strip];
            const auto face=native.begin.face(horizontal,line,strip);
            if(face.negative>=0)inventory[face.negative]-=q;if(face.positive>=0)inventory[face.positive]+=q;
            if(face.negative<0)report.exteriorFlux-=q;if(face.positive<0)report.exteriorFlux+=q;}
    }
    if(retained)*retained=std::move(published);
    return inventory;
}
}
std::vector<double> advanceInventory(const Curve& old,const Source& end,const Flow& flow,const TimeLinearFlux& native,
    double dt,const std::vector<double>& previous,FluxReport& report,FaceVolumes* retained){
    return integratePhaseFlux(&old,&end,flow,native,dt,previous,report,retained);
}
void backgroundFlux(const TimeLinearFlux& native,double dt,FluxReport& report,FaceVolumes& retained){
    integratePhaseFlux(nullptr,nullptr,{},native,dt,std::vector<double>(native.begin.cellCount(),0.),report,&retained);
}
} // namespace InterfaceParametricTransport
