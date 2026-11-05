#include "Core/Geometry/AxisymmetricNeckCap.h"
#include "Core/Geometry/UnitPolynomialRoots.h"
#include "Core/Mesh/CartesianMeasure.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace LAMP::Geometry
{
namespace
{
using Jet=AxisymmetricNeckCap::AxialJet;
using Terms=std::vector<UnitPolynomialTerm>;
void require(bool ok,const char* why)
{if(!ok)throw std::runtime_error(std::string("[AxisymmetricNeckCap] ")+why);}
// Retain endpoint Hermite data as separate exact-expansion terms. Rounded
// power coefficients would destroy the exact shared and reflection jets.
Terms terms(Jet a,Jet b,int derivative=0)
{
    Terms result(4);
    result[0].coefficients={{1,0,-3,2,0,0,0}};result[0].factor=a.value;
    result[1].coefficients={{0,1,-2,1,0,0,0}};result[1].factor=a.derivative;
    result[2].coefficients={{0,0,3,-2,0,0,0}};result[2].factor=b.value;
    result[3].coefficients={{0,0,-1,1,0,0,0}};result[3].factor=b.derivative;
    for(int order=0;order<derivative;++order)for(auto& term:result) {
        for(int k=0;k<6;++k)term.coefficients[k]=(k+1)*term.coefficients[k+1];
        term.coefficients[6]=0;
    }
    return result;
}
struct Interval {double lo=0,hi=0;};
double down(double x){return std::nextafter(x,-std::numeric_limits<double>::infinity());}
double up(double x){return std::nextafter(x,std::numeric_limits<double>::infinity());}
Interval add(Interval a,Interval b){return {down(a.lo+b.lo),up(a.hi+b.hi)};}
Interval multiply(Interval a,Interval b)
{
    const double v[]{a.lo*b.lo,a.lo*b.hi,a.hi*b.lo,a.hi*b.hi};
    return {down(*std::min_element(v,v+4)),up(*std::max_element(v,v+4))};
}
// Only bounds the represented Hermite polynomial. The source fitting and
// subsequent native-patch sampling have separate, explicitly checked errors.
Interval range(const Terms& p,Interval t)
{
    Interval sum;
    for(const auto& term:p) {
        Interval v;
        for(int k=3;k>=0;--k)v=add(multiply(v,t),{term.coefficients[k],term.coefficients[k]});
        sum=add(sum,multiply(v,{static_cast<double>(term.factor),static_cast<double>(term.factor)}));
    }
    require(std::isfinite(sum.lo)&&std::isfinite(sum.hi),"unrepresentable polynomial range");
    return sum;
}
double value(const Terms& p,double t)
{
    long double v=0;
    for(const auto& term:p) {
        long double part=term.coefficients[3];
        for(int k=2;k>=0;--k)part=part*t+term.coefficients[k];
        v+=part*term.factor;
    }
    return static_cast<double>(v);
}
Interval extrema(Jet a,Jet b)
{
    const auto p=terms(a,b);auto bounds=range(p,{0,0});
    const auto include=[&](Interval t){const auto q=range(p,t);bounds.lo=std::min(bounds.lo,q.lo);bounds.hi=std::max(bounds.hi,q.hi);};
    include({1,1});
    const auto roots=isolateUnitPolynomial(terms(a,b,1));
    for(const auto& root:roots.roots)include({root.lower,root.upper});
    return bounds;
}
double hermite(Jet a,Jet b,double t,int order=0)
{
    if(order==0&&(t==0||t==1))return t==0?a.value:b.value;
    if(order==1&&(t==0||t==1))return t==0?a.derivative:b.derivative;
    return value(terms(a,b,order),t);
}
struct Bridge
{
    Jet a,b,c,middle;
    Bridge(const AxisymmetricNeckCap::LeftTrace& left,Jet knot):a(left.axis),middle(knot)
    {
        const double dv=left.outer.value-left.axis.value,dd=left.outer.derivative-left.axis.derivative;
        b={3*dv-left.outerRadial,3*dd-left.outerRadialAxial};
        c={left.outerRadial-2*dv,left.outerRadialAxial-2*dd};
    }
    std::array<double,3> coefficients(double t,int order=0)const
    {return {{hermite(a,{-middle.value,-middle.derivative},t,order),hermite(b,{1,0},t,order),hermite(c,{0,0},t,order)}};}
    void requireMonotone()const
    {
        require(extrema(b,{1,0}).lo>0
            &&extrema({2*b.value+3*c.value,2*b.derivative+3*c.derivative},{2,0}).lo>0,
            "source scalar bridge is not strictly increasing in radius");
    }
    double radialFraction(double t)const
    {
        const auto p=coefficients(t);
        if(p[0]>=0)return 0;
        if(p[0]+p[1]+p[2]<=0)return 1;
        UnitPolynomialTerm term;term.factor=1;term.coefficients={{p[0],0,p[1],p[2],0,0,0}};
        const auto roots=isolateUnitPolynomial({term});
        require(!roots.identicallyZero&&roots.roots.size()==1,"bridge radial slice lacks a single simple root");
        const double r=roots.roots.front().point;
        require(r>0&&r<1,"bridge radial root is unresolved at a cell boundary");
        return r*r;
    }
};
struct Quadrature {long double value=0,error=0;};
long double gauss(const Bridge& bridge,double lo,double hi)
{
    constexpr double x[]{-.96028985649753623168,-.79666647741362673959,-.52553240991632898582,-.18343464249564980494,
        .18343464249564980494,.52553240991632898582,.79666647741362673959,.96028985649753623168};
    constexpr double w[]{.10122853629037625915,.22238103445337447054,.31370664587788728734,.36268378337836198297,
        .36268378337836198297,.31370664587788728734,.22238103445337447054,.10122853629037625915};
    long double sum=0;
    // A radial cubic gives r^2 = A + O(A^(3/2)) near an axis crossing.
    // Smoothstep approaches BOTH endpoints quadratically, so these fractional
    // powers become regular in the quadrature parameter. The Jacobian is
    // retained; phase roots and volume tolerances are not altered.
    for(int k=0;k<8;++k) {
        const double s=(x[k]+1)/2,t=lo+(hi-lo)*s*s*(3-2*s);
        sum+=w[k]*bridge.radialFraction(t)*6*s*(1-s);
    }
    return sum*(hi-lo)/2;
}
Quadrature integrateBridge(const Bridge& bridge,double lo,double hi,long double coarse,long double tolerance,int depth)
{
    const double mid=lo+(hi-lo)/2;
    require(mid>lo&&mid<hi,"bridge quadrature cannot resolve an axial interval");
    const auto a=gauss(bridge,lo,mid),b=gauss(bridge,mid,hi),fine=a+b,error=std::abs(fine-coarse);
    if(error<=tolerance)return {fine,error};
    if(depth>=14) {
        std::ostringstream message;message.precision(17);
        message<<"[AxisymmetricNeckCap] bridge volume quadrature did not converge; interval="<<lo<<','<<hi
            <<" value="<<fine<<" error="<<error<<" tolerance="<<tolerance<<" middle_radius_squared="<<bridge.middle.value;
        throw std::runtime_error(message.str());
    }
    const auto left=integrateBridge(bridge,lo,mid,a,tolerance/2,depth+1);
    const auto right=integrateBridge(bridge,mid,hi,b,tolerance/2,depth+1);
    return {left.value+right.value,left.error+right.error};
}
Quadrature bridgeVolume(const Bridge& bridge,long double tolerance)
{
    // Bracketing iterates may be empty/full over part of the first cell.
    // Split at both axis and radial-wall crossings before smooth quadrature;
    // the final cap separately requires strict interior radial containment.
    std::vector<double> cuts{0,1};
    for(double r:{0.,1.}) {
        const Jet left{bridge.a.value+r*(bridge.b.value+bridge.c.value),
            bridge.a.derivative+r*(bridge.b.derivative+bridge.c.derivative)};
        const Jet right{r-bridge.middle.value,-bridge.middle.derivative};
        const auto roots=isolateUnitPolynomial(terms(left,right));
        require(!roots.identicallyZero,"bridge lies identically on a cell boundary");
        for(const auto& root:roots.roots)if(root.point>0&&root.point<1)cuts.push_back(root.point);
    }
    std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
    Quadrature result;
    for(std::size_t i=1;i<cuts.size();++i) {
        const auto q=integrateBridge(bridge,cuts[i-1],cuts[i],gauss(bridge,cuts[i-1],cuts[i]),
            tolerance*(cuts[i]-cuts[i-1]),0);
        result.value+=q.value;result.error+=q.error;
    }
    return result;
}
}
AxisymmetricNeckCap AxisymmetricNeckCap::reconstruct(const Input& input)
{
    const double hr=input.radialWidth,h=input.axialWidth;
    require(std::isfinite(hr)&&hr>0&&std::isfinite(h)&&h>0&&std::isfinite(2*h),"invalid cell widths");
    require(std::isfinite(input.leftRadiusSquared)&&input.leftRadiusSquared>0
        &&std::isfinite(input.leftRadiusSquaredDerivative),"invalid left radius jet");
    const long double radialSquared=static_cast<long double>(hr)*hr;
    const long double capacity=(CartesianMeasure::twoPi/2)*radialSquared*h;
    require(std::isfinite(capacity)&&capacity>0&&std::isfinite(1/hr)&&std::isfinite(1/h),"unrepresentable cell scale");
    for(auto v:input.liquidVolumes)require(std::isfinite(v)&&v>0&&v<capacity,"each cell requires a strictly mixed positive inventory");
    const long double a=input.leftRadiusSquared/radialSquared;
    const long double d=static_cast<long double>(input.leftRadiusSquaredDerivative)*h/radialSquared;
    const long double m0=input.liquidVolumes[0]/capacity,m1=input.liquidVolumes[1]/capacity;
    const long double q=2*m0-a-d/6;
    AxisymmetricNeckCap cap;cap.input_=input;
    cap.jets_={{{static_cast<double>(a),static_cast<double>(d)},
        {static_cast<double>((64*m1+13*q)/33),static_cast<double>((384*m1-120*q)/33)},
        {static_cast<double>((17*q-160*m1)/33),0}}};
    cap.completeReport();return cap;
}
AxisymmetricNeckCap AxisymmetricNeckCap::reconstructMatched(const Input& input,const LeftTrace& left)
{
    const double hr=input.radialWidth,h=input.axialWidth,eps=std::numeric_limits<double>::epsilon();
    require(std::isfinite(hr)&&hr>0&&std::isfinite(h)&&h>0&&std::isfinite(2*h),"invalid matched cell widths");
    for(double v:{left.axis.value,left.axis.derivative,left.outer.value,left.outer.derivative,left.outerRadial,left.outerRadialAxial})
        require(std::isfinite(v),"nonfinite matched scalar trace");
    require(left.axis.value<0&&left.outer.value>0,"matched source trace must cross the first radial cell");
    const long double capacity=(CartesianMeasure::twoPi/2)*static_cast<long double>(hr)*hr*h;
    require(std::isfinite(capacity)&&capacity>0,"unrepresentable matched cell scale");
    require(std::isfinite(input.leftRadiusSquared)&&input.leftRadiusSquared>0&&input.leftRadiusSquared<hr*hr
        &&std::isfinite(input.leftRadiusSquaredDerivative),"invalid matched radius jet");
    for(auto v:input.liquidVolumes)require(std::isfinite(v)&&v>0&&v<capacity,"each matched cell requires strictly mixed inventory");
    const long double m0=input.liquidVolumes[0]/capacity,m1=input.liquidVolumes[1]/capacity;
    const auto knot=[&](double b)->Jet{return {b,static_cast<double>((384*m1-120*b)/13)};};
    const Bridge incoming(left,knot(0));incoming.requireMonotone();
    const double radiusSquared=incoming.radialFraction(0),radius=std::sqrt(radiusSquared);
    const auto incomingValue=incoming.coefficients(0),incomingDerivative=incoming.coefficients(0,1);
    const double radiusSlope=-2*(incomingDerivative[0]+radiusSquared*(incomingDerivative[1]+radius*incomingDerivative[2]))
        /(2*incomingValue[1]+3*radius*incomingValue[2]);
    const double suppliedRadius=input.leftRadiusSquared/(hr*hr),suppliedSlope=input.leftRadiusSquaredDerivative*h/(hr*hr);
    require(std::abs(radiusSquared-suppliedRadius)<=4096*eps*std::max(radiusSquared,suppliedRadius)
        &&std::abs(radiusSlope-suppliedSlope)<=4096*eps*std::max({1.,std::abs(radiusSlope),std::abs(suppliedSlope)}),
        "supplied radius jet disagrees with the actual matched scalar trace");
    // The last-cell root, volume and even right derivative leave one free B.
    // dB=(384*M1-120*B)/13, E=(17*B-96*M1)/13. For 0<B<96*M1/17,
    // the last cell has exactly one interior cap and a reflected gas gap.
    // In the first cell dphi/dB=-t^2*(159-146*t)/13 < 0; radial monotonicity
    // therefore makes its liquid volume a monotone function of B.
    double lo=0,hi=static_cast<double>(96*m1/17);
    const long double tolerance=256*eps*m0;
    const auto v0=bridgeVolume(Bridge(left,knot(lo)),tolerance);
    const auto v1=bridgeVolume(Bridge(left,knot(hi)),tolerance);
    require(v0.value+v0.error<m0&&v1.value-v1.error>m0,"source trace and two inventories have no interior cap bracket");
    AxisymmetricNeckCap cap;cap.input_=input;cap.leftTrace_=left;
    cap.jets_[0]={input.leftRadiusSquared/(hr*hr),input.leftRadiusSquaredDerivative*h/(hr*hr)};
    for(int iteration=1;iteration<=60;++iteration) {
        const double b=lo+(hi-lo)/2;require(b>lo&&b<hi,"matched cap radius bracket is exhausted");
        const auto measure=bridgeVolume(Bridge(left,knot(b)),tolerance);
        if(std::abs(measure.value-m0)+measure.error<=1024*eps*m0) {
            cap.jets_[1]=knot(b);cap.jets_[2]={static_cast<double>((17*b-96*m1)/13),0};
            cap.report_.matchingIterations=iteration;cap.report_.matchingQuadratureError=measure.error*capacity;
            cap.completeReport();return cap;
        }
        if(measure.value<m0)lo=b;else hi=b;
    }
    throw std::runtime_error("[AxisymmetricNeckCap] source-matched volume solve did not converge");
}
AxisymmetricNeckCap AxisymmetricNeckCap::reconstructMatchedGap(const Input& input,const LeftTrace& left,double gasValue)
{
    const double hr=input.radialWidth,h=input.axialWidth,eps=std::numeric_limits<double>::epsilon();
    require(std::isfinite(hr)&&hr>0&&std::isfinite(h)&&h>0&&std::isfinite(2*h),"invalid matched-gap cell widths");
    require(std::isfinite(gasValue)&&gasValue>0,"matched gap requires a positive common gas value");
    for(double v:{left.axis.value,left.axis.derivative,left.outer.value,left.outer.derivative,left.outerRadial,left.outerRadialAxial})
        require(std::isfinite(v),"nonfinite matched-gap source trace");
    require(left.axis.value<0&&left.outer.value>0,"matched-gap source trace must cross the first radial cell");
    const long double capacity=(CartesianMeasure::twoPi/2)*static_cast<long double>(hr)*hr*h;
    require(std::isfinite(capacity)&&capacity>0,"unrepresentable matched-gap cell scale");
    for(auto v:input.liquidVolumes)require(std::isfinite(v)&&v>0&&v<capacity,"each matched-gap cell requires strictly mixed inventory");
    require(std::isfinite(input.leftRadiusSquared)&&input.leftRadiusSquared>0&&input.leftRadiusSquared<hr*hr
        &&std::isfinite(input.leftRadiusSquaredDerivative),"invalid matched-gap radius jet");
    const Bridge incoming(left,{0,0});incoming.requireMonotone();
    const double radiusSquared=incoming.radialFraction(0),radius=std::sqrt(radiusSquared);
    const auto p=incoming.coefficients(0),dp=incoming.coefficients(0,1);
    const double slope=-2*(dp[0]+radiusSquared*(dp[1]+radius*dp[2]))/(2*p[1]+3*radius*p[2]);
    const double supplied=input.leftRadiusSquared/(hr*hr),suppliedSlope=input.leftRadiusSquaredDerivative*h/(hr*hr);
    require(std::abs(radiusSquared-supplied)<=4096*eps*std::max(radiusSquared,supplied)
        &&std::abs(slope-suppliedSlope)<=4096*eps*std::max({1.,std::abs(slope),std::abs(suppliedSlope)}),
        "matched-gap radius jet disagrees with its scalar trace");
    const long double m0=input.liquidVolumes[0]/capacity,m1=input.liquidVolumes[1]/capacity;
    const long double tolerance=256*eps*m0;
    // For fixed B, increasing D raises phi in the first-cell interior and
    // decreases its liquid volume. Radial monotonicity is independent of
    // B,D. The clipped bracketing iterates need not be admissible caps; the
    // final candidate passes the complete geometry and volume checks below.
    const auto firstVolume=[&](double b,double d){
        return bridgeVolume(Bridge(left,{b,d}),tolerance);};
    const auto middle=[&](double b) {
        double lo=-1,hi=1;auto a=firstVolume(b,lo),z=firstVolume(b,hi);
        for(int k=0;a.value-a.error<=m0 && k<40;++k){lo*=2;a=firstVolume(b,lo);}
        for(int k=0;z.value+z.error>=m0 && k<40;++k){hi*=2;z=firstVolume(b,hi);}
        require(a.value-a.error>m0&&z.value+z.error<m0,"matched-gap slope has no resolved volume bracket");
        for(int k=0;k<64;++k) {
            const double d=lo+(hi-lo)/2;require(d>lo&&d<hi,"matched-gap slope bracket exhausted");
            const auto v=firstVolume(b,d);
            if(std::abs(v.value-m0)+v.error<=1024*eps*m0)return Jet{b,d};
            if(v.value>m0)lo=d;else hi=d;
        }
        throw std::runtime_error("[AxisymmetricNeckCap] matched-gap slope solve did not converge");
    };
    // q=B H00 + D H10 - gasValue H01. H00,H10 are nonnegative.
    // The first-volume constraint gives increasing D(B), hence the positive
    // last-cell volume is increasing in B as well. Its zero and primitive
    // are evaluated independently of the first-cell radial-root quadrature.
    const auto lastVolume=[&](Jet knot) {
        const auto polynomial=terms(knot,{-gasValue,0});
        const auto roots=isolateUnitPolynomial(polynomial);
        require(!roots.identicallyZero,"matched-gap terminal polynomial vanished");
        std::vector<double> cuts{0,1};
        for(const auto& root:roots.roots)if(root.point>0&&root.point<1)cuts.push_back(root.point);
        std::sort(cuts.begin(),cuts.end());long double volume=0;
        for(std::size_t i=1;i<cuts.size();++i)if(value(polynomial,(cuts[i-1]+cuts[i])/2)>0)
            volume+=integrateUnitPolynomial(polynomial,0,1,cuts[i-1],cuts[i]);
        return volume;
    };
    double lo=0,hi=1;
    require(lastVolume(middle(lo))<m1&&lastVolume(middle(hi))>m1,"matched-gap radius has no interior volume bracket");
    for(int k=0;k<64;++k) {
        const double b=lo+(hi-lo)/2;require(b>lo&&b<hi,"matched-gap radius bracket exhausted");
        const auto knot=middle(b);const auto volume=lastVolume(knot);
        if(std::abs(volume-m1)<=1024*eps*m1) {
            AxisymmetricNeckCap cap;cap.input_=input;cap.leftTrace_=left;cap.prescribedGap_=true;
            cap.jets_={{{supplied,suppliedSlope},knot,{-gasValue,0}}};
            cap.report_.matchingIterations=k+1;cap.completeReport();return cap;
        }
        if(volume<m1)lo=b;else hi=b;
    }
    throw std::runtime_error("[AxisymmetricNeckCap] matched-gap radius solve did not converge");
}
void AxisymmetricNeckCap::completeReport()
{
    const auto& cap=*this;const auto& input=input_;
    const double hr=input.radialWidth,h=input.axialWidth;
    const long double radialSquared=static_cast<long double>(hr)*hr,capacity=(CartesianMeasure::twoPi/2)*radialSquared*h;
    for(const auto& j:cap.jets_)require(std::isfinite(j.value)&&std::isfinite(j.derivative),"unrepresentable cap jet");
    const auto first=leftTrace_?Interval{}:extrema(cap.jets_[0],cap.jets_[1]);
    const auto last=extrema(cap.jets_[1],cap.jets_[2]);
    if(leftTrace_) {
        const auto& left=*leftTrace_;const auto middle=jets_[1];
        Bridge(left,middle).requireMonotone();
        require(extrema(left.axis,{-middle.value,-middle.derivative}).hi<0,"matched bridge develops an extra axis cap");
        require(extrema(left.outer,{1-middle.value,-middle.derivative}).lo>0,"matched bridge crosses its radial wall");
    } else {
        require(first.lo>0,"first cell would develop an extra cap or negative radius squared");
        require(first.hi<1,"cap reaches the radial wall of its inventory cells");
    }
    require(last.hi<1,"cap reaches the radial wall of its inventory cells");
    require(cap.jets_[2].value<0,"last cell has no reflected gas gap");
    const auto p=terms(cap.jets_[1],cap.jets_[2]);
    const auto roots=isolateUnitPolynomial(p);
    require(!roots.identicallyZero&&roots.roots.size()==1,"cap requires exactly one simple axial root");
    const auto root=roots.roots.front();
    const double eps=std::numeric_limits<double>::epsilon();
    require(root.lower>0&&root.upper<1&&(prescribedGap_||std::abs(root.point-closureFraction)<512*eps)
        &&root.upper-root.lower<1024*eps,prescribedGap_?"closure root is not resolved inside the terminal cell":
            "closure root is not resolved at the declared half-cell cutoff");
    const auto slope=range(terms(cap.jets_[1],cap.jets_[2],1),{root.lower,root.upper});
    require(slope.hi<0,"cap has a degenerate or wrongly oriented axis tip");
    auto& report=report_;
    const auto tipBounds=multiply({h,h},add({1,1},{root.lower,root.upper}));
    report.tipLocalZ=h*(1+root.point);report.tipLowerZ=tipBounds.lo;report.tipUpperZ=tipBounds.hi;
    report.maxRadiusUpperBound=leftTrace_?hr:up(hr*up(std::sqrt(std::max(first.hi,last.hi))));
    const double tipSlope=value(terms(cap.jets_[1],cap.jets_[2],1),root.point);
    report.axisCurvature=static_cast<double>(4*h/(radialSquared*(-tipSlope)));
    require(std::isfinite(report.axisCurvature)&&report.axisCurvature>0,"unrepresentable axis curvature");
    report.scalarSecondDerivativeJump=leftTrace_?
        (-Bridge(*leftTrace_,jets_[1]).coefficients(1,2)[0]-value(terms(jets_[1],jets_[2],2),0))/(h*h):
        (value(terms(cap.jets_[0],cap.jets_[1],2),1)-value(terms(cap.jets_[1],cap.jets_[2],2),0))/(h*h);
    require(std::isfinite(report.scalarSecondDerivativeJump),"unrepresentable seam curvature jump");
    const auto rootValues=range(p,{root.lower,root.upper});
    const long double rootVolumeBound=capacity*(root.upper-root.lower)*std::max(std::abs(rootValues.lo),std::abs(rootValues.hi));
    for(int i=0;i<2;++i) {
        report.representedVolumes[i]=capacity*integrateUnitPolynomial(terms(cap.jets_[i],cap.jets_[i+1]),0,1,0,i?root.point:1);
        if(i==0&&leftTrace_) {
            const auto q=bridgeVolume(Bridge(*leftTrace_,jets_[1]),256*eps*input.liquidVolumes[0]/capacity);
            report.representedVolumes[0]=q.value*capacity;report.matchingQuadratureError=q.error*capacity;
        }
        report.volumeResiduals[i]=report.representedVolumes[i]-input.liquidVolumes[i];
        report.volumeTolerances[i]=4096*eps*input.liquidVolumes[i];
        require(std::isfinite(report.representedVolumes[i])
            &&std::abs(report.volumeResiduals[i])+(i?rootVolumeBound:report.matchingQuadratureError)<=report.volumeTolerances[i],
            "represented cap cannot retain each cell inventory within relative arithmetic tolerance");
    }
}
AxisymmetricNeckCap::ScalarJet AxisymmetricNeckCap::scalar(double r,double z) const
{
    const double hr=input_.radialWidth,h=input_.axialWidth;
    require(std::isfinite(r)&&r>=0&&r<=hr&&std::isfinite(z)&&z>=0&&z<=2*h,"query outside cap band");
    const int cell=z<h?0:1;const double t=(z-cell*h)/h;
    const auto a=jets_[cell],b=jets_[cell+1];
    const double q=t==0?a.value:t==1?b.value:value(terms(a,b),t);
    const double dq=t==0?a.derivative:t==1?b.derivative:value(terms(a,b,1),t);
    ScalarJet result{(r/hr)*(r/hr)-q,2*(r/hr)/hr,-dq/h,2/hr/hr,0,-value(terms(a,b,2),t)/h/h};
    if(cell==0&&leftTrace_) {
        const Bridge bridge(*leftTrace_,jets_[1]);const double u=r/hr;
        const auto p=bridge.coefficients(t),dp=bridge.coefficients(t,1),ddp=bridge.coefficients(t,2);
        result={p[0]+u*u*(p[1]+u*p[2]),u*(2*p[1]+3*u*p[2])/hr,
            (dp[0]+u*u*(dp[1]+u*dp[2]))/h,(2*p[1]+6*u*p[2])/hr/hr,
            u*(2*dp[1]+3*u*dp[2])/hr/h,(ddp[0]+u*u*(ddp[1]+u*ddp[2]))/h/h};
    }
    require(std::isfinite(result.value)&&std::isfinite(result.dr)&&std::isfinite(result.dz)
        &&std::isfinite(result.drr)&&std::isfinite(result.drz)&&std::isfinite(result.dzz),"unrepresentable scalar jet");
    return result;
}
}
