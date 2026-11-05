#include "Core/Geometry/QuadraticCellCut.h"
#include "Core/Geometry/QuadraticSurfaceMeasure.h"
#include "Core/Geometry/QuadraticFixedCornerSensitivity.h"
#include "Core/Geometry/Detail/QuadraticGraphExact.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cfenv>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>
#include <sstream>
namespace Q=LAMP::Geometry::QuadraticCellCut;
namespace S=LAMP::Geometry::QuadraticSurfaceMeasure;
namespace E=LAMP::Geometry::Detail::QuadraticExact;
namespace {
unsigned checks=0,guards=0;
void check(bool ok,const char* message){++checks;if(!ok)throw std::runtime_error(message);}
void near(long double a,long double b,long double tolerance,const char* message){check(std::isfinite(a)&&std::isfinite(b)&&std::abs(a-b)<=tolerance,message);}
template<class F>void rejects(F f,const char* message){bool rejected=false;try{f();}catch(const std::invalid_argument&){rejected=true;}check(rejected,message);++guards;}
using Values=std::vector<long double>;using Jets=std::vector<Q::Jet2>;
long double width(const Q::Rectangle& b){return std::max(static_cast<long double>(b.upper.x)-b.lower.x,static_cast<long double>(b.upper.y)-b.lower.y);}
Q::Point moved(Q::Point p,Q::Point v,long double lambda){return {p.x+lambda*v.x,p.y+lambda*v.y};}
bool same(Q::Point p,Q::Point q){return p.x==q.x&&p.y==q.y;}
Jets entries(const Q::FixedPortVolumeBoundaryJet& j){Jets out(11);out[0]=j.eta;out[1]=j.slope;out[2]=j.span;
    for(int i=0;i<4;++i)out[3+i]=j.cornerImplicit[i];for(int i=0;i<2;++i)for(int a=0;a<2;++a)out[7+2*i+a]=j.portTangent[i][a];
    for(const auto& stationary:j.stationaryFace)if(stationary)for(const auto& side:stationary->slack)out.push_back(side);return out;}
// Independent VALUE probes of the actual trial polynomial. No jet rule,
// numerical derivative or approximate candidate enters production here.
Values values(const Q::Rectangle& box,const Q::Reconstruction& graph,Q::Point p0,Q::Point p1,bool stationaryEvents=true){
    const auto& g=*graph.exact;const auto h=width(box);const auto parameter=[&](Q::Point p){return E::dot(E::subtract(E::vector(p),g.origin),g.tangent)/g.determinant;};
    const auto t0=parameter(p0),t1=parameter(p1);const E::Scalar span=t1>t0?t1-t0:t0-t1;
    Values out(11);out[0]=E::rounded(g.q[2]*E::scalar(h));out[1]=E::rounded(g.q[1]);out[2]=E::rounded(span/E::scalar(h));
    const std::array<Q::Point,4> corners{{{box.lower.x,box.lower.y},{box.upper.x,box.lower.y},{box.upper.x,box.upper.y},{box.lower.x,box.upper.y}}};
    for(int i=0;i<4;++i)out[3+i]=Q::implicitValue(box,graph,corners[i])/h;
    for(int i=0;i<2;++i){const auto t=i?t1:t0;const E::Scalar slope=g.q[1]+2*g.q[2]*t;
        out[7+2*i]=E::rounded(g.tangent.x+g.normal.x*slope);out[8+2*i]=E::rounded(g.tangent.y+g.normal.y*slope);}
    if(stationaryEvents)for(const auto& stationary:Q::stationaryFaceBoundaries(box,graph))if(stationary)for(const auto value:stationary->slack)out.push_back(value);return out;
}
void events(const Q::Rectangle& box,const Q::Reconstruction& graph,Q::Point first,Q::Point last,bool rz){
    const auto s=S::evaluate(box,graph,box,rz);check(s.patches.size()==1,"actual trial changed the original one-arc topology");const auto& p=s.patches.front();
    check((same(p.begin,first)&&same(p.end,last))||(same(p.begin,last)&&same(p.end,first)),"actual trial lost the original input port events");
}
void sameBoundaryJet(const Q::FixedPortVolumeBoundaryJet& retained,const Q::FixedPortVolumeBoundaryJet& strict){
    const auto a=entries(retained),b=entries(strict);check(a.size()==b.size(),"retained certificate changed boundary jet event coverage");
    for(std::size_t i=0;i<a.size();++i)check(a[i].value==b[i].value&&a[i].first==b[i].first&&a[i].second==b[i].second,
        "retained certificate changed a strict jet value or directional derivative");
    for(int axis=0;axis<2;++axis){check(bool(retained.stationaryFace[axis])==bool(strict.stationaryFace[axis]),
            "retained certificate changed actual stationary eligibility");
        if(!retained.stationaryFace[axis])continue;
        const auto& x=retained.stationaryFace[axis]->current;const auto& y=strict.stationaryFace[axis]->current;
        check(x.parameter==y.parameter&&same(x.position,y.position)&&x.inOriginalArc==y.inOriginalArc&&
            x.side==y.side&&x.branchSign==y.branchSign&&x.slack==y.slack&&x.normalizedCoeffGradient==y.normalizedCoeffGradient,
            "retained certificate changed actual stationary boundary metadata");}
}
void certificateEquivalent(const Q::Rectangle& box,const Q::Reconstruction& current,Q::Point p0,Q::Point p1,
    Q::Point v0,Q::Point v1,long double target,bool rz,long double eta){
    const auto solved=Q::fixedPortVolumeGraph(box,current.shape,p0,p1,target,rz,eta);
    check(bool(solved.graph.exact->fixedPortCertificate),"fully certified fixed-port solve did not retain its original certificate");
    auto strict=solved.graph;
    strict.exact=std::make_shared<LAMP::Geometry::Detail::QuadraticGraphExact>(*solved.graph.exact);
    check(strict.exact.get()!=solved.graph.exact.get()&&strict.exact->fixedPortCertificate==solved.graph.exact->fixedPortCertificate,
        "strict reference must copy the polynomial payload while retaining the opaque token value");
    // The certificate binds the original payload address. A deep copy with
    // the same opaque token therefore exercises the original strict path.
    sameBoundaryJet(Q::fixedPortVolumeBoundaryJet(box,solved.graph,p0,p1,v0,v1,target,rz),
        Q::fixedPortVolumeBoundaryJet(box,strict,p0,p1,v0,v1,target,rz));
}
void stationaryGradient(const Q::Rectangle&,const Q::Reconstruction&);
void taylor(const Q::Rectangle& box,Q::Shape shape,Q::Point p0,Q::Point p1,Q::Point v0,Q::Point v1,long double eta,bool rz,bool rotated){
    const auto current=Q::fromFixedPorts(box,shape,p0,p1,eta);if(rotated)check(current.exact->determinant!=1,"rotated native frame must exercise true D != 1");
    const auto target=Q::graphMeasure(box,current,rz).volume;const auto completeJet=Q::fixedPortVolumeBoundaryJet(box,current,p0,p1,v0,v1,target,rz);const auto jet=entries(completeJet);
    certificateEquivalent(box,current,p0,p1,v0,v1,target,rz,eta);
    const auto center=values(box,current,p0,p1);events(box,current,p0,p1,rz);
    for(std::size_t i=0;i<jet.size();++i)near(jet[i].value,center[i],128*std::numeric_limits<long double>::epsilon()*std::max(1.L,std::abs(center[i])),"jet value differs from actual current polynomial");
    // Two fixed small steps distinguish true second derivatives from lambda^2
    // coefficients: the actual Taylor remainder must have cubic reduction.
    constexpr long double lambda=.01L;std::vector<long double> firstError(jet.size()),secondError(jet.size());
    for(int level=0;level<2;++level){const auto step=std::ldexp(lambda,-level);auto& error=level?secondError:firstError;
        for(int direction:{-1,1}){const auto l=direction*step;const auto a=moved(p0,v0,l),b=moved(p1,v1,l);
            const auto trial=Q::fixedPortVolumeGraph(box,current.shape,a,b,target,rz,eta);events(box,trial.graph,a,b,rz);const auto actual=values(box,trial.graph,a,b);check(actual.size()==jet.size(),"small trial changed stationary-event eligibility");
            for(std::size_t i=0;i<jet.size();++i){const auto prediction=jet[i].value+l*jet[i].first+l*l*jet[i].second/2;
                error[i]=std::max(error[i],std::abs(actual[i]-prediction));}}
    }
    for(std::size_t i=0;i<jet.size();++i){const auto scale=std::max({1.L,std::abs(jet[i].value),std::abs(jet[i].first),std::abs(jet[i].second)});
        const auto verify=[&](bool ok,const char*gate){
            if(ok){check(true,gate);return;}
            std::ostringstream why;why.precision(std::numeric_limits<long double>::max_digits10);
            why<<gate<<" fixture(rz="<<rz<<",sign="<<shape.sign<<",rotated="<<rotated<<",eta="<<eta<<") component="<<i
                <<" jet=("<<jet[i].value<<','<<jet[i].first<<','<<jet[i].second<<") error=("<<firstError[i]<<','<<secondError[i]<<") scale="<<scale
                <<" p0=("<<p0.x<<','<<p0.y<<") p1=("<<p1.x<<','<<p1.y<<") v0=("<<v0.x<<','<<v0.y<<") v1=("<<v1.x<<','<<v1.y<<')'
                <<" q=("<<current.offset<<','<<current.shape.slope<<','<<current.shape.second<<")";
            if(i>=11){std::size_t component=11;for(int axis=0;axis<2;++axis)if(completeJet.stationaryFace[axis]){
                const auto&metadata=completeJet.stationaryFace[axis]->current;
                for(int side=0;side<2;++side,++component)if(component==i)why<<" stationary(axis="<<axis<<",side="<<side<<",id="<<metadata.side[side]
                    <<",branchSign="<<metadata.branchSign[side]<<",t="<<metadata.parameter<<",position="<<metadata.position.x<<','<<metadata.position.y
                    <<",inOriginalArc="<<metadata.inOriginalArc<<')';}}
            check(false,why.str().c_str());
        };
        // The original 11 fields retain their original absolute Taylor gate.
        // A stationary coordinate contains -b^2/(4a); its third derivative is
        // not bounded by max(value,first,second,1). No such C3 bound was part
        // of the new event specification. All fields retain the unchanged
        // two-scale cubic reduction below, plus center/zero-motion checks.
        if(i<11)verify(firstError[i]<=2e-6L*scale,"second-order model disagrees with actual fixed-volume trial");
        verify(secondError[i]<=.2L*firstError[i]+1024*std::numeric_limits<long double>::epsilon()*scale,"fixed-volume boundary Taylor remainder lacks cubic reduction");}
    const auto stationary=entries(Q::fixedPortVolumeBoundaryJet(box,current,p0,p1,{0,0},{0,0},target,rz));
    for(const auto& j:stationary)check(j.first==0&&j.second==0,"zero port motion creates nonzero boundary derivatives");
    stationaryGradient(box,current);
}
Q::Reconstruction changedCoefficient(const Q::Rectangle& box,const Q::Reconstruction& graph,int coefficient,long double delta){
    // Public legacy lift supplies an independent value probe of the chosen
    // physical polynomial coefficient. It is never a production jet path.
    const auto& e=*graph.exact;const auto h=E::scalar(width(box));auto q=e.q;
    q[coefficient]+=E::scalar(delta)*(coefficient==0?h:coefficient==1?E::Scalar(1):E::Scalar(1)/h);
    Q::Reconstruction native;native.shape=graph.shape;native.offset=E::rounded(q[0]*e.determinant);
    native.shape.slope=E::rounded(q[1]);native.shape.second=E::rounded(2*q[2]/e.determinant);
    return Q::recenter(box,box,native);
}
void stationaryGradient(const Q::Rectangle& box,const Q::Reconstruction& graph){
    const auto current=Q::stationaryFaceBoundaries(box,graph);auto phaseFlip=graph;phaseFlip.shape.sign=-phaseFlip.shape.sign;const auto opposite=Q::stationaryFaceBoundaries(box,phaseFlip);
    for(int axis=0;axis<2;++axis){check(bool(current[axis])==bool(opposite[axis]),"stationary geometry depends on phase label");if(!current[axis])continue;
        const auto& event=*current[axis];check(event.side[0]==2*axis&&event.side[1]==2*axis+1,"physical stationary side ids");
        check(event.position.x==opposite[axis]->position.x&&event.position.y==opposite[axis]->position.y&&event.branchSign==opposite[axis]->branchSign&&event.slack==opposite[axis]->slack,"phase flip changed stationary branch identity");
        constexpr long double step=1e-5L;
        for(int c=0;c<3;++c){const auto plus=Q::stationaryFaceBoundaries(box,changedCoefficient(box,graph,c,step));const auto minus=Q::stationaryFaceBoundaries(box,changedCoefficient(box,graph,c,-step));
            check(plus[axis]&&minus[axis],"coefficient probe changed actual stationary eligibility");for(int side=0;side<2;++side){
                check(plus[axis]->branchSign[side]==event.branchSign[side]&&minus[axis]->branchSign[side]==event.branchSign[side],"coefficient probe crossed the original stationary branch");
                // X*_j/h = O_j/h + N_j*theta0 - (T_j+N_j*theta1)^2/(4*N_j*eta).
                // Only eta has a reciprocal term: its finite central difference
                // equals the derivative divided by 1-(step/eta)^2 exactly.
                const auto centralDifference=(plus[axis]->slack[side]-minus[axis]->slack[side])/(2*step);
                const auto eta=E::rounded(graph.exact->q[2]*E::scalar(width(box)));
                const auto referenceFactor=c==2?1-(step/eta)*(step/eta):1.L;
                const auto actual=centralDifference*referenceFactor,expected=event.normalizedCoeffGradient[side][c];
                if(!(std::isfinite(actual)&&std::isfinite(expected)&&std::abs(actual-expected)<=2e-7L)){
                    std::ostringstream why;why.precision(std::numeric_limits<long double>::max_digits10);const auto&e=*graph.exact;
                    why<<"stationary normalized coefficient gradient axis="<<axis<<" side="<<side<<" coefficient="<<c
                        <<" actual="<<actual<<" expected="<<expected<<" residual="<<actual-expected<<" tolerance="<<2e-7L<<" step="<<step
                        <<" centralDifference="<<centralDifference<<" referenceFactor="<<referenceFactor
                        <<" plusSlack="<<plus[axis]->slack[side]<<" minusSlack="<<minus[axis]->slack[side]<<" currentSlack="<<event.slack[side]
                        <<" current=(t="<<event.parameter<<",p="<<event.position.x<<','<<event.position.y<<",inArc="<<event.inOriginalArc<<",branchSign="<<event.branchSign[side]<<')'
                        <<" q=("<<E::rounded(e.q[0])<<','<<E::rounded(e.q[1])<<','<<E::rounded(e.q[2])<<") eta="<<E::rounded(e.q[2]*E::scalar(width(box)))
                        <<" D="<<E::rounded(e.determinant)<<" T=("<<E::rounded(e.tangent.x)<<','<<E::rounded(e.tangent.y)<<") N=("<<E::rounded(e.normal.x)<<','<<E::rounded(e.normal.y)<<')';
                    for(int direction:{-1,1}){const auto probe=changedCoefficient(box,graph,c,direction*step);const auto&pe=*probe.exact;const auto&metadata=direction>0?*plus[axis]:*minus[axis];
                        why<<" probe["<<direction<<"]=(q="<<E::rounded(pe.q[0])<<','<<E::rounded(pe.q[1])<<','<<E::rounded(pe.q[2])
                            <<";t="<<metadata.parameter<<";p="<<metadata.position.x<<','<<metadata.position.y<<";branchSign="<<metadata.branchSign[side]<<')';}
                    check(false,why.str().c_str());
                }else check(true,"stationary normalized coefficient gradient");}}
    }
}
void faceInteriorBirth(){
    // repair18 actual cell24 frame/ports. This single fixed shallower source
    // is legal; no search over eta, lambda or floating neighbours is used.
    const Q::Rectangle box{{0,.375},{.125,.5}};Q::Shape shape;shape.frameNormal=Vector2D{.056153609131858642,.99842214127164997};
    const Q::Point right{.125,.49930440509411561L},left{0,.499984740652413L},v0{0,0},v1{0,1e-5L};constexpr long double eta=-.0070L;
    const auto source=Q::fromFixedPorts(box,shape,right,left,eta);const auto target=Q::graphMeasure(box,source,false).volume;events(box,source,right,left,false);
    const auto jet=Q::fixedPortVolumeBoundaryJet(box,source,right,left,v0,v1,target,false);check(jet.stationaryFace[1].has_value(),"actual near-top stationary event missing");
    const auto& top=jet.stationaryFace[1]->slack[1];const auto& metadata=jet.stationaryFace[1]->current;
    check(metadata.inOriginalArc&&metadata.side[1]==3&&metadata.branchSign[1]==1&&top.value>0,"legal source lacks its actual inward Top stationary branch");
    check(top.value+top.first+top.second/2<0,"new stationary jet does not preannounce actual face-interior birth at lambda=1");
    // This continuation is a MODEL DIAGNOSTIC only, not a certified inventory
    // candidate. The strict original fixed-volume API below must reject it.
    const auto a=moved(right,v0,1),b=moved(left,v1,1);const auto diagnostic=Q::fromFixedPorts(box,shape,a,b,jet.eta.value+jet.eta.first+jet.eta.second/2);
    const auto old=values(box,source,right,left,false),continued=values(box,diagnostic,a,b,false);
    for(int corner=0;corner<4;++corner)check(old[3+corner]*continued[3+corner]>0,"old corner bound already detects this face-interior event");
    for(int port=0;port<2;++port){const auto index=7+2*port;check(old[index]*continued[index]>0,"old native-normal port tangent already detects this event");}
    check(S::evaluate(box,diagnostic,box,false).patches.size()==2,"diagnostic source does not expose the actual two Top crossings");
    rejects([&]{(void)Q::fixedPortVolumeGraph(box,source.shape,a,b,target,false,eta);},"original actual lambda=1 fixed-volume gate accepted the two-arc branch");
    stationaryGradient(box,source);
}
void zeroGapConvention(){
    const Q::Rectangle box{{.25,0},{.375,.125}};Q::Shape shape;const Q::Point a{.25,.1171875},b{.375,.1171875};
    const auto graph=Q::fromFixedPorts(box,shape,a,b,-.25L);events(box,graph,a,b,false);const auto stationary=Q::stationaryFaceBoundaries(box,graph);
    check(stationary[1]&&stationary[1]->inOriginalArc&&stationary[1]->position.y==.125L&&stationary[1]->slack[1]==0&&stationary[1]->branchSign[1]==1,"exact in-arc zero-gap convention lost inward branch");
    auto flipped=graph;flipped.shape.sign=-1;const auto other=Q::stationaryFaceBoundaries(box,flipped);check(other[1]&&other[1]->branchSign==stationary[1]->branchSign,"zero-gap convention changes under phase flip");
    auto legacy=graph;legacy.exact.reset();rejects([&]{(void)Q::stationaryFaceBoundaries(box,legacy);},"stationary event query accepted unencoded legacy graph");
}
void cap(bool rz,int sign){
    const Q::Rectangle box{{.25,0},{.375,.125}};const Q::Point first{.28125,.125},last{.34375,.125},v0{-.0078125,0},v1{.0078125,0};
    Q::Shape shape;shape.axis=1;shape.sign=sign;constexpr long double eta=.125L;
    const auto current=Q::fromFixedPorts(box,shape,first,last,eta);const auto target=Q::graphMeasure(box,current,rz).volume;
    const auto j=Q::fixedPortVolumeBoundaryJet(box,current,first,last,v0,v1,target,rz);
    const auto rate=(v1.x-v0.x)/(last.x-first.x);
    // Same top face: chord selected volume is exactly 0 or capacity.
    // Antisymmetric velocities keep r0+r1 fixed, so XY and true RZ cap
    // stock both imply eta proportional to span^-3 on this branch.
    near(j.eta.first,-3*rate*j.eta.value,1e-12L*std::abs(j.eta.value),"same-face cap eta first derivative violates fixed inventory");
    near(j.eta.second,12*rate*rate*j.eta.value,1e-12L*std::abs(j.eta.value),"same-face cap eta second derivative is not a true second derivative");
    near(j.span.first,(v1.x-v0.x)/width(box),1e-14L,"cap span first derivative");check(j.span.second==0,"affine same-face span has nonzero second derivative");
    taylor(box,shape,first,last,v0,v1,eta,rz,false);
}
void fixedCorner(){
    const Q::Rectangle box{{0,0},{1,1}};Q::Shape shape;
    const Q::Point first{0,0},last{1,.5L},zero{0,0},velocity{0,1.L/8};
    const auto graph=Q::fromFixedPorts(box,shape,first,last,0);
    constexpr long double target=1.L/4;
    near(Q::graphMeasure(box,graph,false).volume,target,1e-12L*target,"fixed LL source has the independent true quarter-cell phase volume");
    events(box,graph,first,last,false);
    const auto jet=Q::fixedPortVolumeBoundaryJet(box,graph,first,last,zero,velocity,target,false);
    certificateEquivalent(box,graph,first,last,zero,velocity,target,false,0);
    check(jet.eta.value==0&&jet.eta.first==3.L/8&&jet.eta.second==0,"fixed-corner constant-volume curvature jet has its exact dyadic derivative");
    check(jet.slope.value==.5L&&jet.slope.first==1.L/8&&jet.slope.second==0,"fixed-corner centered slope jet has its exact dyadic derivative");
    const auto& corner=jet.cornerImplicit[0];
    check(corner.value==0&&corner.first==0&&corner.second==0,"fixed LL implicit value and both directional derivatives remain exact zero");
    const auto sensitivity=S::differentiateWithFixedCornerEndpoints(box,graph,box,false,S::All);
    check(sensitivity.surface.patches.size()==1,"fixed-corner sensitivity measures the original one actual arc");
    const auto referenceArea=std::sqrt(5.L/4),referenceDerivative=(.5L/referenceArea)/8;
    near(sensitivity.surface.integral.physicalMeasure,referenceArea,1e-12L*referenceArea,"fixed-corner sensitivity retains original actual physical area");
    // Centered physical (offset,slope,second) chain follows the SAME fixed
    // LL and exact volume: (-1/32,1/8,3/4). At t=-1/2 its value is zero.
    constexpr std::array<long double,3> chain{{-1.L/32,1.L/8,3.L/4}};
    check(chain[0]-chain[1]/2+chain[2]/8==0,"analytic coefficient chain exactly preserves the fixed corner");
    long double pulled=0;for(int k=0;k<3;++k)pulled+=sensitivity.derivative[k]*chain[k];
    near(pulled,referenceDerivative,1e-12L*referenceDerivative,"fixed-corner physical covector pulls back to the independent area derivative");
    // One actual dyadic trial, not a derivative probe: y' = 1/2 +
    // lambda*(3x/4-1/4), and integral(3x/4-1/4)^2 dx=1/16. Since
    // sqrt(1+u^2)''<=1, the first-order area remainder is <=lambda^2/32.
    constexpr long double lambda=1.L/16;
    const auto movedLast=moved(last,velocity,lambda);
    const auto trial=Q::fromFixedPorts(box,shape,first,movedLast,3*lambda/8);
    events(box,trial,first,movedLast,false);
    near(Q::graphMeasure(box,trial,false).volume,target,1e-12L*target,"actual fixed-corner trial retains the independent quarter-cell stock");
    check(trial.offset-graph.offset==lambda*chain[0]&&trial.shape.slope-graph.shape.slope==lambda*chain[1]&&
        trial.shape.second-graph.shape.second==lambda*chain[2],"actual dyadic trial follows the independent centered physical coefficient chain");
    const auto actual=S::evaluate(box,trial,box,false).integral.physicalMeasure;
    check(std::abs(actual-referenceArea-lambda*referenceDerivative)<=lambda*lambda/32+1e-12L*referenceArea,
        "actual fixed-corner area obeys the independent analytic first-order remainder bound");
    rejects([&]{(void)S::differentiate(box,graph,box,false);},"original unconstrained surface sensitivity must still reject a corner endpoint");
    rejects([&]{(void)Q::fixedPortVolumeBoundaryJet(box,graph,first,last,{1.L/8,1.L/16},velocity,target,false);},
        "nonzero corner motion is rejected even when tangent to the original physical line");
}
void invalid(){
    const Q::Rectangle box{{.25,0},{.375,.125}};Q::Shape shape;const Q::Point p0{.25,.0625},p1{.375,.0625};
    const auto graph=Q::fromFixedPorts(box,shape,p0,p1,.03L);const auto target=Q::graphMeasure(box,graph,false).volume;
    const auto solved=Q::fixedPortVolumeGraph(box,shape,p0,p1,target,false,.03L);
    check(bool(solved.graph.exact->fixedPortCertificate),"invalid-binding probes require the actual certified payload");
    auto stale=solved.graph;stale.offset=std::nextafter(stale.offset,std::numeric_limits<long double>::infinity());
    rejects([&]{(void)Q::fixedPortVolumeBoundaryJet(box,stale,p0,p1,{0,0},{0,0},target,false);},"retained token bypassed stale native-view rejection");
    rejects([&]{(void)Q::fixedPortVolumeBoundaryJet(box,solved.graph,p0,p1,{0,0},{0,0},target/2,false);},"retained token bypassed a different target-volume certificate");
    auto differentBox=box;differentBox.upper.x+=.125L;
    rejects([&]{(void)Q::fixedPortVolumeBoundaryJet(differentBox,solved.graph,p0,p1,{0,0},{0,0},target,false);},"retained token bypassed changed actual box/port incidence");
    rejects([&]{(void)Q::fixedPortVolumeBoundaryJet(box,graph,p0,p1,{.01,0},{0,0},target,false);},"normal port motion accepted as actual face sliding");
    auto legacy=graph;legacy.exact.reset();rejects([&]{(void)Q::fixedPortVolumeBoundaryJet(box,legacy,p0,p1,{0,0},{0,0},target,false);},"legacy graph accepted by exact boundary jet");
    rejects([&]{(void)Q::fixedPortVolumeBoundaryJet(box,graph,p0,p0,{0,0},{0,0},target,false);},"duplicate port events accepted");
    rejects([&]{(void)Q::fixedPortVolumeBoundaryJet(box,graph,p0,p1,{0,std::numeric_limits<long double>::quiet_NaN()},{0,0},target,false);},"nonfinite port direction accepted");
    rejects([&]{(void)Q::fixedPortVolumeBoundaryJet(box,graph,p0,p1,{0,0},{0,0},0,false);},"pure target accepted on current mixed branch");
    const Q::Rectangle axisBox{{0,0},{.125,.25}};const Q::Point a{0,.1},b{.125,.15};const auto axis=Q::fromAxisPorts(axisBox,1,a,b);const auto volume=Q::graphMeasure(axisBox,axis,true).volume;
    rejects([&]{(void)Q::fixedPortVolumeBoundaryJet(axisBox,axis,a,b,{0,0},{0,0},volume,true);},"physical-axis graph accepted by ordinary free-port jet");
}
}
void filteredSlopeSamples()
{
    const Q::Rectangle box{{0,0},{1,1}};Q::Shape shape;shape.axis=1;
    const auto ordinary=Q::fromFixedPorts(box,shape,{0,.25},{1,.75},.125L);
    Q::Shape rotated;rotated.frameNormal=Vector2D{.27488060480878,.9614783685034};
    const auto framed=Q::fromFixedPorts(box,rotated,{0,.25},{1,.75},.125L);
    const auto reference=[&](const Q::Reconstruction& graph,long double t){
        const auto& g=*graph.exact;const auto slope=E::rounded(LAMP::Geometry::Detail::graphSlope(g,E::scalar(t)));
        const auto tangent=E::native(g.tangent),normal=E::native(g.normal);
        const auto metric=std::sqrt(E::rounded(g.determinant)),length=std::hypot(1.L,slope);
        const auto q=S::sample(box,graph,t,false);
        const auto identical=[](long double a,long double b){return a==b&&std::signbit(a)==std::signbit(b);};
        check(identical(q.unitNormal.x,graph.shape.sign*(normal.x-slope*tangent.x)/(metric*length)),"filtered slope changed correctly rounded normal x");
        check(identical(q.unitNormal.y,graph.shape.sign*(normal.y-slope*tangent.y)/(metric*length)),"filtered slope changed correctly rounded normal y");
        check(identical(q.meridionalCurvature,-graph.shape.sign*graph.shape.second/(metric*length*length*length)),"filtered slope changed correctly rounded curvature");
    };
    // Ordinary and nonunit-frame queries, an exact midpoint between adjacent
    // native slopes, and strong cancellation use one independent rational rule.
    for(const auto* graph:{&ordinary,&framed})for(const auto t:{-.375L,-.125L,0.L,.125L,.375L,
        std::nextafter(1.L,2.L),-2.L,std::nextafter(-2.L,0.L)})reference(*graph,t);
    const auto small=Q::fromFixedPorts(box,shape,{0,.25},{1,.25},.5L);
    reference(small,std::numeric_limits<long double>::denorm_min());
    const auto underflow=Q::fromFixedPorts(box,shape,{0,.25},{1,.25},.125L);
    rejects([&]{(void)E::rounded(LAMP::Geometry::Detail::graphSlope(*underflow.exact,E::scalar(std::numeric_limits<long double>::denorm_min())));},"independent original exact query rejects a strictly sub-half-denorm slope");
    rejects([&]{(void)S::sample(box,underflow,std::numeric_limits<long double>::denorm_min(),false);},"exact nonzero slope underflow must retain its original rejection");
    const auto oldRound=std::fegetround();
    if(std::fesetround(FE_DOWNWARD)==0){
        struct RestoreRound{int value;~RestoreRound(){std::fesetround(value);}}restore{oldRound};
        reference(ordinary,.125L);
    }
}
void nativeRoundingCases()
{
    namespace E=LAMP::Geometry::Detail::QuadraticExact;
    // Independent binary rounding cells, including both tie parities.
    const int precision=std::numeric_limits<long double>::digits;
    for(int exponent:{-1000,0,1000}){
        const auto base=std::scalbn(1.L,exponent),ulp=std::scalbn(1.L,exponent-precision+1);
        for(int offset:{0,1,2}){
            const auto lower=base+offset*ulp,upper=std::nextafter(lower,std::numeric_limits<long double>::infinity());
            const E::Scalar middle=(E::scalar(lower)+E::scalar(upper))/2;
            const E::Scalar delta=(E::scalar(upper)-E::scalar(lower))/8;
            const auto tied=offset%2?upper:lower;
            for(int sign:{-1,1}){
                check(E::rounded(sign*(middle-delta))==sign*lower,"normal exact ratio below midpoint");
                check(E::rounded(sign*middle)==sign*tied,"normal exact ratio ties-to-even");
                check(E::rounded(sign*(middle+delta))==sign*upper,"normal exact ratio above midpoint");
            }
        }
    }
    check(E::rounded(E::Scalar(1)/3)==1.L/3,"nonbinary normal rational rounds once");
    check(E::rounded(E::scalar(std::numeric_limits<long double>::min()))==std::numeric_limits<long double>::min(),"minimum normal view");
    check(E::rounded(E::scalar(std::numeric_limits<long double>::denorm_min()))==std::numeric_limits<long double>::denorm_min(),"original subnormal fallback");
    check(E::rounded(E::scalar(std::numeric_limits<long double>::max()))==std::numeric_limits<long double>::max(),"original overflow-binade fallback");
}
int main(){try{
    nativeRoundingCases();
    filteredSlopeSamples();
    const Q::Rectangle box{{.25,0},{.375,.125}};const Q::Point p0{.25,.0625},p1{.375,.0625},v0{0,.015625},v1{0,-.0078125};
    Q::Shape shape;shape.frameNormal=Vector2D{.27488060480878,.9614783685034};
    for(bool rz:{false,true})for(int sign:{1,-1}){shape.sign=sign;taylor(box,shape,p0,p1,v0,v1,.03L,rz,true);cap(rz,sign);}invalid();faceInteriorBirth();zeroGapConvention();fixedCorner();
    std::cout<<"quadratic fixed-port boundary jet PASS checks="<<checks<<" guards="<<guards<<" XY/RZ rotated and same-face cap\n";return 0;
}catch(const std::exception& error){std::cerr<<"quadratic fixed-port boundary jet FAIL "<<error.what()<<'\n';return 1;}}
