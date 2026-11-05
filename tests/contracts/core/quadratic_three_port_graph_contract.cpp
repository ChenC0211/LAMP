#include "Core/Geometry/QuadraticCellCut.h"
#include "Core/Geometry/QuadraticCornerBoundary.h"
#include "Core/Geometry/QuadraticSurfaceMeasure.h"
#include "Core/Geometry/QuadraticCurveMeasure.h"
#include "Core/Geometry/Detail/QuadraticGraphExact.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
namespace Q=LAMP::Geometry::QuadraticCellCut;
namespace S=LAMP::Geometry::QuadraticSurfaceMeasure;
namespace C=LAMP::Geometry::QuadraticCurveMeasure;
namespace E=LAMP::Geometry::Detail::QuadraticExact;
namespace {
unsigned checks=0,guards=0;
void check(bool ok,const char* message){++checks;if(!ok)throw std::runtime_error(message);}
void near(long double a,long double b,long double tolerance,const char* message){check(std::isfinite(a)&&std::isfinite(b)&&std::abs(a-b)<=tolerance,message);}
template<class F>void rejects(F f,const char* message){bool rejected=false;try{f();}catch(const std::invalid_argument&){rejected=true;}check(rejected,message);++guards;}
bool same(Q::Point a,Q::Point b){return a.x==b.x&&a.y==b.y;}
long double width(const Q::Rectangle& b){return std::max(static_cast<long double>(b.upper.x)-b.lower.x,static_cast<long double>(b.upper.y)-b.lower.y);}
E::Scalar parameter(const Q::Reconstruction& g,Q::Point p){const auto& e=*g.exact;return E::dot(E::subtract(E::vector(p),e.origin),e.tangent)/e.determinant;}
std::array<E::Scalar,3> normalized(const Q::Rectangle& b,const Q::Reconstruction& graph){const auto& g=*graph.exact;const auto h=E::scalar(width(b));return {{g.q[0]/h,g.q[1],g.q[2]*h}};}
void jacobian(const Q::Rectangle& box,const Q::ThreePortGraphResult& result,const std::array<Q::Point,3>& ports){
    const auto& g=*result.graph.exact;const auto h=E::scalar(width(box));std::array<E::Scalar,3> t;
    for(int i=0;i<3;++i)t[i]=parameter(result.graph,ports[i]);
    // Differentiate q(t_l)=n_l. These nine scalar constraints test EACH
    // physical-coordinate column without moving a port off its real face.
    for(int moved=0;moved<3;++moved)for(int coordinate=0;coordinate<2;++coordinate)for(int equation=0;equation<3;++equation){
        const auto& j=result.normalizedCoefficientJacobian[moved][coordinate];
        const E::Scalar lhs=E::scalar(j[0])*h+E::scalar(j[1])*t[equation]+E::scalar(j[2])*t[equation]*t[equation]/h;
        E::Scalar rhs=0;if(equation==moved){const auto tangent=coordinate?g.tangent.y:g.tangent.x,normal=coordinate?g.normal.y:g.normal.x;
            rhs=(normal-(g.q[1]+2*g.q[2]*t[equation])*tangent)/g.determinant;}
        near(E::rounded(lhs),E::rounded(rhs),128*std::numeric_limits<long double>::epsilon()*std::max(1.L,std::abs(E::rounded(rhs))),"three-port full Jacobian violates interpolation differential identity");
    }
    // Only one REAL tangent port is displaced at a time. For two fixed
    // ports, q_new=q+c*L_i exactly, with c=[dn-q'(ti)dt-q2*dt^2]/L_i(ti+dt).
    // This finite-step identity includes all truncation terms; no arbitrary
    // central-difference error threshold or derivative step extrapolation.
    constexpr long double lambda=1.L/64;
    for(int moved=0;moved<3;++moved){Q::Point velocity{};
        if(ports[moved].x==box.lower.x||ports[moved].x==box.upper.x)velocity.y=width(box)/64;
        else velocity.x=width(box)/64;
        const int j=(moved+1)%3,k=(moved+2)%3;const E::Scalar denominator=(t[moved]-t[j])*(t[moved]-t[k]);
        const std::array<E::Scalar,3> l{{t[j]*t[k]/denominator,-(t[j]+t[k])/denominator,E::Scalar(1)/denominator}};
        for(int direction:{-1,1}){auto displaced=ports;const auto step=direction*lambda;
            displaced[moved]={ports[moved].x+step*velocity.x,ports[moved].y+step*velocity.y};
            const auto actual=Q::fromThreePorts(box,result.graph.shape,displaced);const auto delta=E::subtract(E::vector(displaced[moved]),E::vector(ports[moved]));
            const E::Scalar dt=E::dot(delta,g.tangent)/g.determinant,dn=E::dot(delta,g.normal)/g.determinant,u=t[moved]+dt;
            const E::Scalar value=l[0]+u*(l[1]+u*l[2]),correction=(dn-(g.q[1]+2*g.q[2]*t[moved])*dt-g.q[2]*dt*dt)/value;
            const auto q=normalized(box,actual.graph);
            const std::array<E::Scalar,3> scale{{E::Scalar(1)/h,1,h}};
            for(int c=0;c<3;++c)near(E::rounded(q[c]),E::rounded((g.q[c]+correction*l[c])*scale[c]),
                128*std::numeric_limits<long double>::epsilon()*std::max(1.L,std::abs(E::rounded(q[c]))),"actual same-face tangent displacement violates exact finite-step interpolation identity");
        }
    }
}
void events(const Q::Rectangle& box,const Q::Reconstruction& graph,const std::array<Q::Point,3>& ports,bool rz){
    const auto surface=S::evaluate(box,graph,box,rz);check(!surface.patches.empty(),"three-port graph lost all actual arcs");
    std::array<unsigned,3> uses{};
    for(const auto& arc:surface.patches){const auto curve=C::Curve::fromGraph(Q::materialCurveData(box,graph,arc.lower,arc.upper));
        for(int end=0;end<2;++end){const auto physical=end?arc.end:arc.begin;const auto material=C::position(curve,end);
            check(material.x==static_cast<double>(physical.x)&&material.y==static_cast<double>(physical.y),"material export changed original three-port side-event bits");
            for(int i=0;i<3;++i)if(same(physical,ports[i]))++uses[i];}}
    for(int i=0;i<3;++i){check(Q::implicitValue(box,graph,ports[i])==0,"exact three-port interpolation residual");check(uses[i]==1,"one original input port did not retain its own actual side event");}
}
constexpr Q::Point left{1.25L,.49853515625L},top0{41.L/32,.5L},top1{43.L/32,.5L},right{1.5L,.48291015625L};
const Q::Rectangle coarse{{1.25,.25},{1.5,.5}},fine{{1.25,.5},{1.375,.625}};
void dividedIncidence(const Q::Rectangle& box,const Q::ThreePortGraphResult& source,const std::array<Q::Point,3>& ports,bool rotated){
    const auto h=width(box);const Q::Point known=ports[1],companion=ports[2],fourth{(known.x+companion.x)/2,box.upper.y};
    const auto exact=Q::sameFaceDividedIncidence(box,source.graph,known,companion);
    check(exact.face==3&&exact.residual==0&&same(exact.companion,companion),"same-face Vieta bridge lost the actual distinct Top companion event");
    const auto value=Q::sameFaceDividedIncidence(box,source.graph,known,fourth);
    near(value.residual,(fourth.x-companion.x)/h,64*std::numeric_limits<long double>::epsilon(),"non-root fourth node divided incidence");
    const auto coincident=Q::sameFaceDividedIncidence(box,source.graph,known,known);
    near(coincident.residual,(known.x-companion.x)/h,64*std::numeric_limits<long double>::epsilon(),"coincident nodes must not create a fake zero incidence");
    check(coincident.residual!=0,"distinct quadratic roots were silently merged by coincidence");
    const auto& g=*source.graph.exact;const E::Scalar eh=E::scalar(h),normal=g.normal.y,tangent=g.tangent.y,projection=g.tangent.x;
    const E::Scalar b=tangent+normal*g.q[1];
    const std::array<E::Scalar,3> expected{{0,g.determinant/(projection*eh*g.q[2]),-g.determinant*b/(projection*eh*eh*normal*g.q[2]*g.q[2])}};
    for(int k=0;k<3;++k)near(value.normalizedCoefficientJacobian[k],E::rounded(expected[k]),128*std::numeric_limits<long double>::epsilon()*std::max(1.L,std::abs(E::rounded(expected[k]))),"normalized divided-incidence coefficient Jacobian units");
    for(int point=0;point<2;++point)for(int coordinate=0;coordinate<2;++coordinate){const E::Scalar derivative=(coordinate?g.tangent.y:g.tangent.x)/(projection*eh);
        near(value.pointJacobian[point][coordinate],E::rounded(derivative),128*std::numeric_limits<long double>::epsilon()*std::max(1.L,std::abs(E::rounded(derivative))),"same-face point Jacobian units");}
    check(value.pointJacobian[1][0]==1/h,"fourth-node actual Top tangent derivative must be exactly 1/h");
    // Only the constrained chart is differentiated: known remains an actual
    // factory root. Moving root2 changes the real companion by the same amount;
    // moving known/root1 or the Left port leaves the other Top root fixed.
    for(int moved=0;moved<3;++moved){const int coordinate=moved==0?1:0;long double chain=0,magnitude=0;
        for(int k=0;k<3;++k){const auto term=value.normalizedCoefficientJacobian[k]*source.normalizedCoefficientJacobian[moved][coordinate][k];chain+=term;magnitude+=std::abs(term);}
        if(moved==1){chain+=value.pointJacobian[0][coordinate];magnitude+=std::abs(value.pointJacobian[0][coordinate]);}
        const auto reference=moved==2?-1/h:0.L;
        near(chain,reference,256*std::numeric_limits<long double>::epsilon()*(1+magnitude+std::abs(reference)),"known-root-constrained three-port incidence chain rule");
        for(int direction:{-1,1}){auto changed=ports;const auto step=direction*h/4096;
            if(coordinate)changed[moved].y+=step;else changed[moved].x+=step;
            const auto next=Q::fromThreePorts(box,source.graph.shape,changed);
            const auto actual=Q::sameFaceDividedIncidence(box,next.graph,changed[1],fourth);
            // Exact finite-step relation from the TWO supplied true Top roots,
            // not a central-difference approximation to a derivative.
            near(actual.residual,(fourth.x-changed[2].x)/h,128*std::numeric_limits<long double>::epsilon(),"actual finite three-port change violates constrained companion incidence");}
    }
    const auto shifted=Q::sameFaceDividedIncidence(box,source.graph,known,{fourth.x+h/4096,fourth.y});
    near(shifted.residual,(fourth.x+h/4096-companion.x)/h,128*std::numeric_limits<long double>::epsilon(),"actual fourth-node tangent finite displacement");
    auto flip=source.graph;flip.shape.sign=-flip.shape.sign;const auto opposite=Q::sameFaceDividedIncidence(box,flip,known,fourth);
    check(opposite.residual==value.residual&&same(opposite.companion,value.companion)&&opposite.normalizedCoefficientJacobian==value.normalizedCoefficientJacobian&&opposite.pointJacobian==value.pointJacobian,"phase label changed same-curve divided incidence");
    if(rotated){const Q::Point onLeft{box.lower.x,(static_cast<long double>(box.lower.y)+box.upper.y)/2};
        const auto vertical=Q::sameFaceDividedIncidence(box,source.graph,ports[0],onLeft);
        check(vertical.face==0&&vertical.pointJacobian[1][1]==1/h,"actual Left tangent divided incidence derivative");
        const auto moved=Q::sameFaceDividedIncidence(box,source.graph,ports[0],{onLeft.x,onLeft.y+h/4096});
        near(moved.residual-vertical.residual,1.L/4096,256*std::numeric_limits<long double>::epsilon()*std::max(1.L,std::abs(vertical.residual)),"supporting-line companion view must retain exact fourth-node tangent change");}
    rejects([&]{(void)Q::sameFaceDividedIncidence(box,source.graph,fourth,companion);},"non-root known point accepted by same-face incidence");
    rejects([&]{(void)Q::sameFaceDividedIncidence(box,source.graph,known,ports[0]);},"different original faces accepted by divided incidence");
    rejects([&]{(void)Q::sameFaceDividedIncidence(box,source.graph,known,{box.lower.x,box.upper.y});},"fourth corner accepted as an open-face node");
    rejects([&]{(void)Q::sameFaceDividedIncidence(box,source.graph,known,{std::numeric_limits<long double>::infinity(),box.upper.y});},"nonfinite fourth node accepted");
    auto stale=source.graph;stale.offset=std::nextafter(stale.offset,std::numeric_limits<long double>::infinity());
    rejects([&]{(void)Q::sameFaceDividedIncidence(box,stale,known,fourth);},"stale encoded curve accepted by divided incidence");
    auto legacy=source.graph;legacy.exact.reset();rejects([&]{(void)Q::sameFaceDividedIncidence(box,legacy,known,fourth);},"unencoded curve accepted by divided incidence");
}
void dividedLimits(){
    const Q::Rectangle box{{.25,0},{.375,.125}};Q::Shape shape;const Q::Point a{.25,.1171875},b{.375,.1171875},touch{.3125,.125};
    const auto doubleRoot=Q::fromFixedPorts(box,shape,a,b,-.25L);const auto limit=Q::sameFaceDividedIncidence(box,doubleRoot,touch,touch);
    check(limit.residual==0&&same(limit.companion,touch),"true double-root limit was confused with coincident non-root factor");
    const auto plane=Q::fromFixedPorts(box,shape,{.28125,.125},{.34375,.125},0);
    rejects([&]{(void)Q::sameFaceDividedIncidence(box,plane,{.28125,.125},{.34375,.125});},"zero quadratic same-face polynomial accepted");
    Q::Shape transverse;transverse.axis=0;
    const auto noProjection=Q::fromFixedPorts(box,transverse,{.3125,.125},{.375,.0625},.03L);
    // In an orthogonal physical frame T_other=0 also implies N_face=0:
    // this is the same legitimate projection/a degeneracy, not a fake chart.
    rejects([&]{(void)Q::sameFaceDividedIncidence(box,noProjection,{.3125,.125},{.34375,.125});},"zero actual face-tangent projection accepted");
}
void fixture(bool rotated){
    Q::Shape shape;shape.axis=1;if(rotated)shape.frameNormal=Vector2D{.27488060480878,.9614783685034};
    const std::array<Q::Point,3> ports{{left,top0,top1}};const auto result=Q::fromThreePorts(coarse,shape,ports);
    auto oppositeShape=shape;oppositeShape.sign=-shape.sign;const auto oppositeFactory=Q::fromThreePorts(coarse,oppositeShape,ports);
    check(oppositeFactory.graph.shape.sign==-shape.sign&&oppositeFactory.normalizedCoefficientJacobian==result.normalizedCoefficientJacobian,"three-port factory phase flip changed coefficient geometry or Jacobian");
    check(result.graph.exact!=nullptr,"three-port factory missing same exact payload");
    if(rotated)check(result.graph.exact->determinant!=1,"rotated test does not exercise actual frame D != 1");
    const auto reference=rotated?Q::fromFixedPorts(coarse,result.graph.shape,left,top1,E::rounded(result.graph.exact->q[2]*E::scalar(width(coarse)))):
        Q::fromFixedPorts(coarse,shape,left,right,-.125L);
    if(!rotated)for(int c=0;c<3;++c)check(result.graph.exact->q[c]==reference.exact->q[c],"dyadic three ports do not recover the exact original fixed-port parabola");
    for(bool rz:{false,true}){
        const auto actual=Q::graphMeasure(coarse,result.graph,rz).volume,expected=Q::graphMeasure(coarse,reference,rz).volume;
        near(actual,expected,1e-12L*std::min(expected,LAMP::Geometry::PlanarCellCut::capacity(coarse,rz)-expected),"three-port same polynomial original phase-volume certificate");
        auto referenceOther=reference;referenceOther.shape.sign=-referenceOther.shape.sign;
        const auto other=Q::graphMeasure(coarse,oppositeFactory.graph,rz).volume,expectedOther=Q::graphMeasure(coarse,referenceOther,rz).volume;
        near(other,expectedOther,1e-12L*expectedOther,"three-port independent opposite-phase own volume");
        const auto a=S::evaluate(coarse,result.graph,coarse,rz),b=S::evaluate(coarse,reference,coarse,rz);
        check(a.patches.size()==2&&b.patches.size()==2,"coarse original parabola requires two real arcs");
        near(a.integral.physicalMeasure,b.integral.physicalMeasure,1e-12L*b.integral.physicalMeasure,"three-port true XY/RZ coarse area");events(coarse,result.graph,ports,rz);
        if(!rotated){const auto child=Q::recenter(coarse,fine,result.graph);const auto exactChild=Q::fromFixedPorts(fine,shape,top0,top1,-.0625L);
            const auto cap=S::evaluate(fine,child,fine,rz);check(cap.patches.size()==1&&(cap.patches[0].beginBoundary&S::Bottom)&&(cap.patches[0].endBoundary&S::Bottom),"same original coarse polynomial does not make the true adjacent fine cap");
            const auto volume=Q::graphMeasure(fine,exactChild,rz).volume;near(Q::graphMeasure(fine,child,rz).volume,volume,1e-12L*volume,"fine same-curve cap physical volume");
            near(cap.integral.physicalMeasure,S::evaluate(fine,exactChild,fine,rz).integral.physicalMeasure,1e-12L*cap.integral.physicalMeasure,"fine same-curve cap physical area");}
    }
    jacobian(coarse,result,ports);
    dividedIncidence(coarse,result,ports,rotated);
    const auto stationary=Q::stationaryFaceBoundaries(coarse,result.graph);check(stationary[1]&&!stationary[1]->inOriginalArc&&stationary[1]->position.y>coarse.upper.y,"two-arc peak must be outside every actual retained arc");
    auto flip=result.graph;flip.shape.sign=-flip.shape.sign;const auto flipped=Q::stationaryFaceBoundaries(coarse,flip);
    check(flipped[1]&&flipped[1]->branchSign==stationary[1]->branchSign&&flipped[1]->slack==stationary[1]->slack,"pure stationary arc-union metadata changed with phase label");
    rejects([&]{(void)Q::fixedPortVolumeBoundaryJet(coarse,result.graph,ports[0],ports[2],{0,0},{0,0},Q::graphMeasure(coarse,result.graph,false).volume,false);},"ordinary fixed-volume boundary jet silently accepted two arcs");
    auto stale=result.graph;stale.offset=std::nextafter(stale.offset,std::numeric_limits<long double>::infinity());
    rejects([&]{(void)Q::graphMeasure(coarse,stale,false);},"three-port stale coefficient view accepted");
}
void cornerKinds(){
    // Same dyadic box, frame and curvature as the existing exact parabola.
    // These factories are pure geometry charts, not accepted phase owners.
    Q::Shape shape;shape.axis=1;const auto h=width(coarse);
    const Q::Point ul{coarse.lower.x,coarse.upper.y},ur{coarse.upper.x,coarse.upper.y};
    const auto contact=Q::fromFixedPorts(coarse,shape,ul,ur,-.125L);
    const auto zero=Q::cornerBoundaryKind(coarse,contact,3);
    check(zero.value==0&&zero.inwardDerivative[0]==-.125L&&zero.inwardDerivative[1]==-1&&zero.isolatedClosed,
        "exact upper-left contact with two strictly negative inward derivatives is isolated closed");
    const auto oppositeCorner=Q::cornerBoundaryKind(coarse,contact,2);
    check(oppositeCorner.value==0&&oppositeCorner.inwardDerivative==zero.inwardDerivative&&oppositeCorner.isolatedClosed,
        "physical upper-right inward x orientation is retained in exact corner classification");
    const auto shifted=[&](long double dy){return Q::fromFixedPorts(coarse,shape,{ul.x,ul.y+dy},{ur.x,ur.y+dy},-.125L);};
    const auto sameSide=Q::cornerBoundaryKind(coarse,shifted(h/8),3);
    check(sameSide.value==-.125L&&sameSide.inwardDerivative==zero.inwardDerivative&&sameSide.isolatedClosed,
        "same-sign exact corner value and inward derivatives remain closed after q0 shift");
    const auto otherSide=Q::cornerBoundaryKind(coarse,shifted(-h/8),3);
    check(otherSide.value==.125L&&otherSide.inwardDerivative==zero.inwardDerivative&&!otherSide.isolatedClosed,
        "opposite-sign corner value is open even when both inward derivatives agree");
    const auto crossing=Q::fromFixedPorts(coarse,shape,ul,{ur.x,ur.y-h/4},-.125L);
    const auto cross=Q::cornerBoundaryKind(coarse,crossing,3);
    check(cross.value==0&&cross.inwardDerivative[0]==.125L&&cross.inwardDerivative[1]==-1&&!cross.isolatedClosed,
        "opposite inward derivative signs identify a nonisolated crossing corner");
    const auto tangent=Q::fromFixedPorts(coarse,shape,ul,{ur.x,ur.y-h/8},-.125L);
    const auto degenerate=Q::cornerBoundaryKind(coarse,tangent,3);
    check(degenerate.value==0&&degenerate.inwardDerivative[0]==0&&degenerate.inwardDerivative[1]==-1&&!degenerate.isolatedClosed,
        "zero exact inward derivative cannot be promoted to strictly isolated closed");
    const std::array<Q::Point,3> ports{{left,top0,top1}};
    const auto original=Q::fromThreePorts(coarse,shape,ports);
    const auto actual=Q::cornerBoundaryKind(coarse,original.graph,3);
    check(actual.value==3.L/512&&actual.inwardDerivative[0]==-1.L/16&&actual.inwardDerivative[1]==-1&&!actual.isolatedClosed,
        "existing actual three-port curve has an open opposite-sign upper-left corner");
    auto flip=contact;flip.shape.sign=-flip.shape.sign;
    const auto flipped=Q::cornerBoundaryKind(coarse,flip,3);
    check(flipped.value==-zero.value&&flipped.inwardDerivative[0]==-zero.inwardDerivative[0]&&
        flipped.inwardDerivative[1]==-zero.inwardDerivative[1]&&flipped.isolatedClosed==zero.isolatedClosed,
        "phase flip reverses signed corner data without changing geometric closed classification");
    auto legacy=contact;legacy.exact.reset();
    rejects([&]{(void)Q::cornerBoundaryKind(coarse,legacy,3);},"unencoded legacy graph accepted by exact corner classification");
    auto stale=contact;stale.offset=std::nextafter(stale.offset,std::numeric_limits<long double>::infinity());
    rejects([&]{(void)Q::cornerBoundaryKind(coarse,stale,3);},"stale graph view accepted by exact corner classification");
    rejects([&]{(void)Q::cornerBoundaryKind(coarse,contact,4);},"invalid physical corner index accepted");
}
void invalid(){
    Q::Shape shape;const std::array<Q::Point,3> ports{{left,top0,top1}};
    auto bad=ports;bad[1]=bad[0];rejects([&]{(void)Q::fromThreePorts(coarse,shape,bad);},"duplicate actual port accepted");
    bad=ports;bad[1]={left.x,.4375L};rejects([&]{(void)Q::fromThreePorts(coarse,shape,bad);},"different points with equal exact tangent accepted");
    bad=ports;bad[0]={1.25,.25};rejects([&]{(void)Q::fromThreePorts(coarse,shape,bad);},"corner accepted as an open-face port");
    bad=ports;bad[0]={1.3,.4};rejects([&]{(void)Q::fromThreePorts(coarse,shape,bad);},"interior point accepted as an actual port");
    bad=ports;bad[0].x=std::numeric_limits<long double>::quiet_NaN();rejects([&]{(void)Q::fromThreePorts(coarse,shape,bad);},"nonfinite port accepted");
    auto frame=shape;frame.frameNormal=Vector2D{0,0};rejects([&]{(void)Q::fromThreePorts(coarse,frame,ports);},"zero frame accepted");
    frame=shape;frame.axis=2;rejects([&]{(void)Q::fromThreePorts(coarse,frame,ports);},"invalid graph axis accepted");
    const Q::Rectangle zeroBox{{.25,0},{.375,.125}};const auto touching=Q::fromFixedPorts(zeroBox,shape,{.25,.1171875},{.375,.1171875},-.25L);
    const auto touch=Q::stationaryFaceBoundaries(zeroBox,touching);check(touch[1]&&touch[1]->inOriginalArc&&touch[1]->branchSign[1]==1&&touch[1]->slack[1]==0,"pure stationary nonempty arc union lost exact in-arc zero-gap convention");
}
}
int main(){try{fixture(false);fixture(true);invalid();dividedLimits();cornerKinds();std::cout<<"quadratic three-port graph PASS checks="<<checks<<" guards="<<guards<<" exact events, full Jacobian, XY/RZ coarse-two-arc/fine-cap\n";return 0;}
catch(const std::exception& error){std::cerr<<"quadratic three-port graph FAIL "<<error.what()<<'\n';return 1;}}
