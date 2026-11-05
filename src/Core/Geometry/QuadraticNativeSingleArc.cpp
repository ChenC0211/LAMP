#include "Core/Geometry/QuadraticNativeSingleArc.h"
#include <algorithm>
#include <cfenv>
#include <cfloat>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace LAMP::Geometry::QuadraticNativeSingleArc
{
namespace
{
using Real=long double;
void require(bool ok,const char* why){if(!ok)throw std::invalid_argument(std::string("[QuadraticNativeSingleArc] ")+why);}
Real down(Real x){return std::nextafter(x,-std::numeric_limits<Real>::infinity());}
Real up(Real x){return std::nextafter(x,std::numeric_limits<Real>::infinity());}
struct R
{
    Real v=0,lo=0,hi=0;
    R()=default;R(Real a):v(a),lo(a),hi(a){}
    R(Real a,Real l,Real h):v(a),lo(l),hi(h){}
};
R bad(){const auto inf=std::numeric_limits<Real>::infinity();return {std::numeric_limits<Real>::quiet_NaN(),-inf,inf};}
bool valid(const R&a){return std::isfinite(a.v)&&std::isfinite(a.lo)&&std::isfinite(a.hi)&&a.lo<=a.v&&a.v<=a.hi;}
bool zero(const R&a){return a.lo==0&&a.hi==0;}
R operator+(const R&a,const R&b){if(zero(a))return b;if(zero(b))return a;return {a.v+b.v,down(a.lo+b.lo),up(a.hi+b.hi)};}
R operator-(const R&a){return {-a.v,-a.hi,-a.lo};}
R operator-(const R&a,const R&b){return a+(-b);}
R operator*(const R&a,const R&b){if(zero(a)||zero(b))return R(0);const std::array<Real,4> p{{a.lo*b.lo,a.lo*b.hi,a.hi*b.lo,a.hi*b.hi}};
    return {a.v*b.v,down(*std::min_element(p.begin(),p.end())),up(*std::max_element(p.begin(),p.end()))};}
R operator/(const R&a,const R&b){if(b.lo<=0&&b.hi>=0)return bad();if(zero(a))return R(0);
    auto value=a*R(1/b.v,down(1/b.hi),up(1/b.lo));value.v=a.v/b.v;return value;}
R square(const R&a){return a*a;}
R root(const R&a){if(a.lo<=0)return bad();return {std::sqrt(a.v),down(std::sqrt(a.lo)),up(std::sqrt(a.hi))};}
int sign(const R&a){return a.lo>0?1:a.hi<0?-1:0;}
Real error(const R&a){return std::max(std::abs(a.v-a.lo),std::abs(a.hi-a.v));}
// Error-free sum is used only to recognize exactly representable AMR centres
// and local coordinate differences. Uncertain inputs still use intervals.
R sum(Real a,Real b){const Real s=a+b;if(!std::isfinite(s))return bad();const Real z=s-a;
    const Real tail=(a-(s-z))+(b-z);if(tail==0)return R(s);return R(a)+R(b);}
struct P {R x,y;};
P add(P a,P b){return {a.x+b.x,a.y+b.y};}
P scale(P a,R b){return {a.x*b,a.y*b};}
R dot(P a,P b){return a.x*b.x+a.y*b.y;}
P input(Point p){return {R(p.x),R(p.y)};}
R coord(P p,int axis){return axis?p.y:p.x;}
Real coord(Point p,int axis){return axis?p.y:p.x;}
P local(Point p,Point origin){return {sum(p.x,-origin.x),sum(p.y,-origin.y)};}
Point value(P p){return {p.x.v,p.y.v};}
struct J
{
    Real v=0,d=0,dd=0;
    J()=default;J(Real a):v(a){}J(Real a,Real b,Real c=0):v(a),d(b),dd(c){}
};
J operator+(J a,J b){return {a.v+b.v,a.d+b.d,a.dd+b.dd};}
J operator-(J a){return {-a.v,-a.d,-a.dd};}
J operator-(J a,J b){return a+(-b);}
J operator*(J a,J b){return {a.v*b.v,a.d*b.v+a.v*b.d,a.dd*b.v+2*a.d*b.d+a.v*b.dd};}
J operator/(J a,J b){require(b.v!=0,"singular native direction quotient");return a*J(1/b.v,-b.d/(b.v*b.v),2*b.d*b.d/(b.v*b.v*b.v)-b.dd/(b.v*b.v));}
struct JP {J x,y;};
Q::Jet2 jet(J a){require(std::isfinite(a.v)&&std::isfinite(a.d)&&std::isfinite(a.dd),"nonfinite native direction jet");return {a.v,a.d,a.dd};}
unsigned mask(const Q::Rectangle&b,Point p){unsigned out=0;if(p.x==b.lower.x)out|=1;if(p.x==b.upper.x)out|=2;if(p.y==b.lower.y)out|=4;if(p.y==b.upper.y)out|=8;return out;}
int faceAxis(unsigned face){return (face&3)?0:1;}
int faceIndex(unsigned face){return face==1?0:face==2?1:face==4?2:3;}
bool same(Point a,Point b){return a.x==b.x&&a.y==b.y;}
bool sameShape(const Q::Shape&a,const Q::Shape&b){return a.axis==b.axis&&a.sign==b.sign&&a.slope==b.slope&&a.second==b.second&&
    a.frameNormal.has_value()==b.frameNormal.has_value()&&(!a.frameNormal||(a.frameNormal->x==b.frameNormal->x&&a.frameNormal->y==b.frameNormal->y));}
bool sameRequest(const Request&a,const Request&b){return a.box.lower.x==b.box.lower.x&&a.box.lower.y==b.box.lower.y&&a.box.upper.x==b.box.upper.x&&a.box.upper.y==b.box.upper.y&&
    sameShape(a.shape,b.shape)&&same(a.ports[0],b.ports[0])&&same(a.ports[1],b.ports[1])&&a.targetPhaseVolume==b.targetPhaseVolume&&a.volumeLiquid==b.volumeLiquid&&
    a.branchEta==b.branchEta&&a.relativeVolumeTolerance==b.relativeVolumeTolerance&&a.minimumSlope==b.minimumSlope&&a.maximumSlope==b.maximumSlope;}
}
namespace Detail
{
struct Data
{
    Request input;Point centre{},tangent{},normal{};
    Real D=0,h=0,q0=0,q1=0,q2=0,eta=0,chordArea=0,capacity=0;
    std::array<unsigned,2> face{};std::array<Real,2> t{},n{};
    std::array<bool,4> chordInside{};bool forward=true;int selectedSign=1;
    S::MomentPatch moment;std::array<std::optional<Q::StationaryFaceBoundary>,2> stationary{};
};
}
namespace
{
void check(const Trial&t,const Detail::Data&d)
{
    require(std::fegetround()==FE_TONEAREST,"native packet queried in a different rounding environment");
    require(sameRequest(t.input,d.input)&&sameShape(t.shape,Q::Shape{d.input.shape.axis,d.input.shape.sign,d.q1,2*d.q2,d.input.shape.frameNormal})&&
        t.offset==d.q0&&t.eta==d.eta&&t.h==d.h&&t.frameD==d.D&&t.capacity==d.capacity,
        "stale or mutated native numerical packet");
    require(t.moment.lower==d.moment.lower&&t.moment.upper==d.moment.upper&&t.moment.m0==d.moment.m0&&t.moment.m1==d.moment.m1&&t.moment.scale==d.moment.scale&&
        same(t.moment.begin,d.moment.begin)&&same(t.moment.end,d.moment.end)&&t.moment.beginBoundary==d.moment.beginBoundary&&t.moment.endBoundary==d.moment.endBoundary,
        "stale native moment packet");
}
void checkMoments(const S::MomentValues&m){for(int i=0;i<3;++i)require(std::isfinite(m.moments[i])&&m.moments[i]>0&&std::isfinite(m.firstDerivative[i])&&std::isfinite(m.lastDerivative[i]),"invalid same-model moments/partials");}
Preparation no(Reason reason,const char*why){return {std::nullopt,reason,why};}
}
Preparation prepare(const Request&r)
{
    require(std::isfinite(r.box.lower.x)&&std::isfinite(r.box.lower.y)&&std::isfinite(r.box.upper.x)&&std::isfinite(r.box.upper.y)&&
        r.box.upper.x>r.box.lower.x&&r.box.upper.y>r.box.lower.y,"invalid physical rectangle");
    require(std::isfinite(r.targetPhaseVolume)&&std::isfinite(r.branchEta)&&std::isfinite(r.relativeVolumeTolerance)&&r.relativeVolumeTolerance>0&&r.relativeVolumeTolerance<1&&
        std::isfinite(r.minimumSlope)&&std::isfinite(r.maximumSlope)&&r.minimumSlope<r.maximumSlope,"invalid native model request");
    // basis() performs the SAME native frame normalization as the exact factory.
    const auto basis=Q::basis(r.shape);
    if(!std::numeric_limits<Real>::is_iec559||FLT_EVAL_METHOD!=0||std::fegetround()!=FE_TONEAREST)
        return no(Reason::Unsupported,"native interval environment is not nearest IEEE arithmetic");
    const auto safeHalf=[](Real a){return a==0||std::fpclassify(a/2)==FP_NORMAL;};
    if(!safeHalf(r.box.lower.x)||!safeHalf(r.box.upper.x)||!safeHalf(r.box.lower.y)||!safeHalf(r.box.upper.y))
        return no(Reason::UncertainRoundoff,"cell midpoint exponent window");
    const auto cx=sum(r.box.lower.x/2,r.box.upper.x/2),cy=sum(r.box.lower.y/2,r.box.upper.y/2);
    if(!valid(cx)||!valid(cy)||cx.lo!=cx.hi||cy.lo!=cy.hi)return no(Reason::UncertainRoundoff,"cell midpoint not exactly native representable");
    const Point centre{cx.v,cy.v};const P T=input(basis.tangent),N=input(basis.normal);
    const R D=dot(T,T),hx=sum(r.box.upper.x,-r.box.lower.x),hy=sum(r.box.upper.y,-r.box.lower.y);
    if(!valid(D)||D.lo<=0||!valid(hx)||!valid(hy))return no(Reason::UncertainRoundoff,"frame/cell metric");
    const R h(std::max(hx.v,hy.v),std::max(hx.lo,hy.lo),std::max(hx.hi,hy.hi)),capacity=hx*hy;
    require(r.targetPhaseVolume>0&&r.targetPhaseVolume<capacity.v,"native fixed-volume chart requires mixed target");
    require(!same(r.ports[0],r.ports[1]),"duplicate physical ports");
    std::array<unsigned,2> faces{};std::array<R,2> t,n;
    for(int i=0;i<2;++i){const auto p=r.ports[i];require(std::isfinite(p.x)&&std::isfinite(p.y)&&p.x>=r.box.lower.x&&p.x<=r.box.upper.x&&p.y>=r.box.lower.y&&p.y<=r.box.upper.y,"port outside physical cell");
        faces[i]=mask(r.box,p);require(faces[i]!=0,"port is not on physical boundary");
        if((faces[i]&(faces[i]-1))!=0)return no(Reason::Unsupported,"corner port belongs to strict chart");
        const auto pLocal=local(p,centre);t[i]=dot(pLocal,T)/D;n[i]=dot(pLocal,N)/D;
        if(!valid(t[i])||!valid(n[i]))return no(Reason::UncertainRoundoff,"port frame coordinates");}
    if(faces[0]==faces[1])return no(Reason::NearEvent,"same-face cap requires strict nonmonotone chart");
    const R difference=t[1]-t[0];const int order=sign(difference);if(!order)return no(Reason::NearEvent,"unresolved tangent order");
    const int first=order>0?0:1,last=1-first;const R lo=t[first],hi=t[last],width=hi-lo;
    const R chordSlope=(n[1]-n[0])/difference,chordOffset=n[0]-chordSlope*t[0];
    if(!valid(chordSlope)||!valid(chordOffset)||width.lo<=0)return no(Reason::UncertainRoundoff,"fixed-port chord condition");
    const int selectedSign=r.shape.sign*(r.volumeLiquid?1:-1);
    const std::array<Point,4> corners{{{r.box.lower.x,r.box.lower.y},{r.box.upper.x,r.box.lower.y},{r.box.upper.x,r.box.upper.y},{r.box.lower.x,r.box.upper.y}}};
    std::array<R,4> ct,cn;std::array<bool,4> inside{};R boxLo,boxHi;
    for(int k=0;k<4;++k){const auto pLocal=local(corners[k],centre);ct[k]=dot(pLocal,T)/D;cn[k]=dot(pLocal,N)/D;
        const auto f=R(selectedSign)*(cn[k]-chordOffset-chordSlope*ct[k]);const int s=sign(f);
        if(!valid(ct[k])||!valid(cn[k])||!valid(f)||!s)return no(Reason::NearEvent,"chord corner predicate unresolved");inside[k]=s<0;
        if(k==0){boxLo=boxHi=ct[k];}else{boxLo={std::min(boxLo.v,ct[k].v),std::min(boxLo.lo,ct[k].lo),std::min(boxLo.hi,ct[k].hi)};
            boxHi={std::max(boxHi.v,ct[k].v),std::max(boxHi.lo,ct[k].lo),std::max(boxHi.hi,ct[k].hi)};}}
    // Actual known crossing identities, not newly rounded line intersections.
    const std::array<unsigned,4> edges{{4,2,8,1}};std::vector<P> polygon;
    for(int k=0;k<4;++k){if(inside[k])polygon.push_back(local(corners[k],{r.box.lower.x,r.box.lower.y}));
        if(inside[k]!=inside[(k+1)%4]){int port=-1;for(int j=0;j<2;++j)if(faces[j]==edges[k])port=j;
            if(port<0)return no(Reason::NearEvent,"chord crossing identities differ from ports");polygon.push_back(local(r.ports[port],{r.box.lower.x,r.box.lower.y}));}}
    R twiceArea;for(std::size_t k=0;k<polygon.size();++k){const auto&a=polygon[k];const auto&b=polygon[(k+1)%polygon.size()];twiceArea=twiceArea+a.x*b.y-b.x*a.y;}
    const R chordArea=twiceArea/R(2),b=-R(selectedSign)*D*width*width*width/(R(6)*h);
    const R eta=(R(r.targetPhaseVolume)-chordArea)/b;
    const R q2=eta/h,q1=chordSlope-q2*(t[0]+t[1]),q0=chordOffset+q2*t[0]*t[1];
    if(!valid(chordArea)||chordArea.lo<=0||!valid(b)||!sign(b)||!valid(eta)||!valid(q0)||!valid(q1)||!valid(q2))
        return no(Reason::UncertainRoundoff,"analytic volume coefficients unresolved");
    // Factored independent bubble measure on the certified original arc.
    const R measured=chordArea-R(selectedSign)*D*eta*width*width*width/(R(6)*h);
    const R defect=measured-R(r.targetPhaseVolume);const Real own=std::min(r.targetPhaseVolume,capacity.lo-r.targetPhaseVolume);
    const Real bound=std::max(std::abs(defect.lo),std::abs(defect.hi));
    if(!valid(measured)||own<=0||!std::isfinite(bound)||bound>r.relativeVolumeTolerance*own||bound>PlanarCellCut::relativeVolumeResolution*capacity.lo)
        return no(Reason::UncertainRoundoff,"native volume enclosure exceeds original own-phase budget");
    // A sufficient interval one-arc domain: both physical coordinates are
    // monotone over EVERY parameter whose projection can occur in the cell.
    // Thus their rectangle preimages are intervals. Known entering/exiting
    // ports give their sole intersection; no missed second side-root/contact.
    std::array<int,2> monotone{};const R tBox(0,boxLo.lo,boxHi.hi);
    for(int axis=0;axis<2;++axis){const R derivative=coord(T,axis)+coord(N,axis)*(q1+R(2)*q2*tBox);monotone[axis]=sign(derivative);
        if(!valid(derivative)||!monotone[axis])return no(Reason::NearEvent,"physical-coordinate monotonicity not separated from zero");}
    const auto inward=[](unsigned face){return face==1||face==4?1:-1;};
    if(monotone[faceAxis(faces[first])]*inward(faces[first])<=0||monotone[faceAxis(faces[last])]*inward(faces[last])>=0)
        return no(Reason::NearEvent,"known ports do not enter/exit original cell");
    const R metric=root(D),m0=q1+R(2)*q2*lo,m1=q1+R(2)*q2*hi;
    if(!valid(metric)||!valid(m0)||!valid(m1)||m0.lo<r.minimumSlope||m0.hi>r.maximumSlope||m1.lo<r.minimumSlope||m1.hi>r.maximumSlope)
        return no(Reason::NearEvent,"slope interval outside supplied moment model");
    if(!std::isfinite(2*q2.v))return no(Reason::UncertainRoundoff,"physical second derivative exponent window");
    Trial result;result.input=r;result.shape=r.shape;result.offset=q0.v;result.shape.slope=q1.v;result.shape.second=2*q2.v;
    result.eta=eta.v;result.h=h.v;result.frameD=D.v;result.capacity=capacity.v;result.selectedVolume=measured.v;result.relativeVolumeErrorBound=bound/own;
    const R dv0=R(selectedSign)*D*h*width,dv1=R(selectedSign)*D*(hi*hi-lo*lo)/R(2),dv2=R(selectedSign)*D*(hi*hi*hi-lo*lo*lo)/(R(3)*h);
    if(!valid(dv0)||!valid(dv1)||!valid(dv2))return no(Reason::UncertainRoundoff,"volume sensitivity enclosure");
    result.selectedDerivative={{dv0.v,dv1.v,dv2.v}};
    auto& moment=result.moment;moment.lower=lo.v;moment.upper=hi.v;moment.begin=r.ports[first];moment.end=r.ports[last];
    moment.beginBoundary=faces[first];moment.endBoundary=faces[last];moment.m0=m0.v;moment.m1=m1.v;moment.scale=(metric*width).v;
    const R qa=q0+q1*lo+q2*lo*lo;
    const std::array<R,3> radius{{R(centre.x)+T.x*lo+N.x*qa,width*(T.x+N.x*m0),N.x*q2*width*width}};
    for(int j=0;j<3;++j){if(!valid(radius[j]))return no(Reason::UncertainRoundoff,"moment radius metadata");moment.radius[j]=radius[j].v;}
    std::array<std::array<R,3>,2> dt{};
    for(int e=0;e<2;++e){const int inputPort=e==0?first:last;const R a=e==0?lo:hi,m=e==0?m0:m1;
        const int axis=faceAxis(faces[inputPort]);const R denominator=coord(T,axis)+coord(N,axis)*m;
        if(!valid(denominator)||!sign(denominator))return no(Reason::NearEvent,"endpoint sensitivity is not simple");
        const std::array<R,3> height{{R(1),a,a*a/R(2)}};for(int k=0;k<3;++k)dt[e][k]=-coord(N,axis)*height[k]/denominator;}
    for(int k=0;k<3;++k){const R dw=dt[1][k]-dt[0][k],dma=R(k==1?1:0)+(k==2?lo:R(0))+R(2)*q2*dt[0][k],
            dmb=R(k==1?1:0)+(k==2?hi:R(0))+R(2)*q2*dt[1][k];
        const R height=k==0?R(1):k==1?lo:lo*lo/R(2),dqa=height+m0*dt[0][k];
        const std::array<R,6> derivatives{{dma,dmb,metric*dw,T.x*dt[0][k]+N.x*dqa,
            dw*(T.x+N.x*m0)+width*N.x*dma,N.x*(R(k==2?.5L:0)*width*width+R(2)*q2*width*dw)}};
        for(const auto&v:derivatives)if(!valid(v))return no(Reason::UncertainRoundoff,"moment derivative enclosure");
        moment.m0Derivative[k]=derivatives[0].v;moment.m1Derivative[k]=derivatives[1].v;moment.scaleDerivative[k]=derivatives[2].v;
        for(int j=0;j<3;++j)moment.radiusDerivative[j][k]=derivatives[3+j].v;}
    const Real determinant=basis.tangent.x*basis.normal.y-basis.tangent.y*basis.normal.x;
    const Real orientation=r.shape.sign*(determinant>0?-1.L:1.L);
    for(int e=0;e<2;++e){const R a=e==0?lo:hi,m=e==0?m0:m1,g=root(R(1)+m*m);const P tangent=add(T,scale(N,m));
        const R length=metric*g,tx=R(orientation)*tangent.x/length,ty=R(orientation)*tangent.y/length,curvature=-R(r.shape.sign)*R(2)*eta/(metric*g*g*g);
        if(!valid(tx)||!valid(ty)||!valid(curvature))return no(Reason::UncertainRoundoff,"seed shape interval");
        result.seed[e]={{{tx.v,ty.v,curvature.v}},m.v,(a/h).v,e==0?moment.begin.x:moment.end.x};}
    for(int k=0;k<4;++k){const R m=q1+R(2)*q2*ct[k],f=R(r.shape.sign)*(cn[k]-q0-q1*ct[k]-q2*ct[k]*ct[k])/h;
        const R gx=R(r.shape.sign*(k==0||k==3?1:-1))*(N.x-m*T.x)/D,gy=R(r.shape.sign*(k<2?1:-1))*(N.y-m*T.y)/D;
        const int sf=sign(f),sx=sign(gx),sy=sign(gy);if(!valid(f)||!valid(gx)||!valid(gy)||!sf||!sx||!sy)
            return no(Reason::NearEvent,"corner gap/classification interval unresolved");
        result.corners[k]={f.v,{{gx.v,gy.v}},sx==sy&&sf==sx};}
    for(int axis=0;axis<2;++axis){if(zero(coord(N,axis)))continue;
        const R a=coord(N,axis)*q2;if(!valid(a)||!sign(a))return no(Reason::UncertainRoundoff,"stationary coefficient interval");
        const R parameter=-(coord(T,axis)+coord(N,axis)*q1)/(R(2)*a),height=q0+q1*parameter+q2*parameter*parameter;
        const P position=add(input(centre),add(scale(T,parameter),scale(N,height)));const R other=coord(position,1-axis);
        const Real otherLo=axis?r.box.lower.x:r.box.lower.y,otherHi=axis?r.box.upper.x:r.box.upper.y;
        if(!valid(parameter)||!valid(position.x)||!valid(position.y))return no(Reason::UncertainRoundoff,"stationary position interval");
        if(other.hi<otherLo||other.lo>otherHi)continue;
        if(!(other.lo>otherLo&&other.hi<otherHi))return no(Reason::NearEvent,"stationary face eligibility interval");
        if(!(parameter.hi<lo.lo||parameter.lo>hi.hi))return no(Reason::NearEvent,"stationary parameter overlaps original arc");
        Q::StationaryFaceBoundary event;event.parameter=parameter.v;event.position=value(position);event.inOriginalArc=false;
        const R coordinate=coord(position,axis),lower(axis?r.box.lower.y:r.box.lower.x),upper(axis?r.box.upper.y:r.box.upper.x),scaled=parameter/h;
        const std::array<R,3> coefficient{{coord(N,axis),coord(N,axis)*scaled,coord(N,axis)*scaled*scaled}};
        for(int side=0;side<2;++side){const R gap=side?upper-coordinate:coordinate-lower;const int branch=sign(gap);
            if(!valid(gap)||!branch)return no(Reason::NearEvent,"stationary side gap interval");
            event.side[side]=2*axis+side;event.branchSign[side]=branch;event.slack[side]=(R(branch)*gap/h).v;
            for(int k=0;k<3;++k)event.normalizedCoeffGradient[side][k]=(R(branch*(side?-1:1))*coefficient[k]).v;}
        result.stationary[axis]=event;}
    auto data=std::make_shared<Detail::Data>();data->input=r;data->centre=centre;data->tangent=basis.tangent;data->normal=basis.normal;
    data->D=D.v;data->h=h.v;data->q0=q0.v;data->q1=q1.v;data->q2=q2.v;data->eta=eta.v;data->chordArea=chordArea.v;data->capacity=capacity.v;
    data->face=faces;data->t={{t[0].v,t[1].v}};data->n={{n[0].v,n[1].v}};data->chordInside=inside;data->forward=order>0;data->selectedSign=selectedSign;data->moment=moment;data->stationary=result.stationary;
    result.data_=std::move(data);return {std::move(result),Reason::Ready,{}};
}
S::Sensitivity area(const Trial&t,const S::MomentValues&m)
{
    require(t.data_!=nullptr,"unprepared native packet");check(t,*t.data_);checkMoments(m);const auto&p=t.data_->moment;
    const Real length=p.scale*m.moments[0];require(std::isfinite(length)&&length>0,"invalid native moment-model area");S::Sensitivity out;
    out.surface.patches.push_back({p.lower,p.upper,p.begin,p.end,p.beginBoundary,p.endBoundary,{length,length,0}});out.surface.integral={length,length,0};
    for(int k=0;k<3;++k){out.derivative[k]=p.scaleDerivative[k]*m.moments[0]+p.scale*(m.firstDerivative[0]*p.m0Derivative[k]+m.lastDerivative[0]*p.m1Derivative[k]);
        require(std::isfinite(out.derivative[k]),"nonfinite native area sensitivity");}return out;
}
Q::FixedPortVolumeBoundaryJet direction(const Trial&t,Point v0,Point v1)
{
    require(t.data_!=nullptr,"unprepared native packet");const auto&d=*t.data_;check(t,d);const std::array<Point,2> velocity{{v0,v1}};
    std::array<JP,2> p;std::array<J,2> a,n;
    for(int i=0;i<2;++i){const auto v=velocity[i];require(std::isfinite(v.x)&&std::isfinite(v.y)&&(faceAxis(d.face[i])==0?v.x==0:v.y==0),"direction leaves actual open face");
        p[i]={{d.input.ports[i].x,v.x},{d.input.ports[i].y,v.y}};
        a[i]={d.t[i],(v.x*d.tangent.x+v.y*d.tangent.y)/d.D};n[i]={d.n[i],(v.x*d.normal.x+v.y*d.normal.y)/d.D};}
    const int first=d.forward?0:1,last=1-first;const J width=a[last]-a[first],chordSlope=(n[1]-n[0])/(a[1]-a[0]),chordOffset=n[0]-chordSlope*a[0];
    const std::array<Point,4> corners{{{d.input.box.lower.x,d.input.box.lower.y},{d.input.box.upper.x,d.input.box.lower.y},{d.input.box.upper.x,d.input.box.upper.y},{d.input.box.lower.x,d.input.box.upper.y}}};
    const std::array<unsigned,4> edges{{4,2,8,1}};std::vector<JP> polygon;
    for(int k=0;k<4;++k){if(d.chordInside[k])polygon.push_back({J(corners[k].x-d.input.box.lower.x),J(corners[k].y-d.input.box.lower.y)});
        if(d.chordInside[k]!=d.chordInside[(k+1)%4]){const int i=d.face[0]==edges[k]?0:1;polygon.push_back({p[i].x-J(d.input.box.lower.x),p[i].y-J(d.input.box.lower.y)});}}
    J chordArea;for(std::size_t k=0;k<polygon.size();++k){const auto&x=polygon[k];const auto&y=polygon[(k+1)%polygon.size()];chordArea=chordArea+(x.x*y.y-y.x*x.y)/J(2);}
    const J b=-J(d.selectedSign*d.D/(6*d.h))*width*width*width,c=chordArea-J(d.input.targetPhaseVolume);
    require(b.v!=0,"singular native fixed-volume branch");const Real etaFirst=-(b.d*d.eta+c.d)/b.v,etaSecond=-(b.dd*d.eta+2*b.d*etaFirst+c.dd)/b.v;
    const J eta{d.eta,etaFirst,etaSecond},q2=eta/J(d.h),q1=chordSlope-q2*(a[0]+a[1]),q0=chordOffset+q2*a[0]*a[1];
    Q::FixedPortVolumeBoundaryJet out;out.eta=jet(eta);out.slope=jet(q1);out.span=jet(width/J(d.h));
    for(int k=0;k<4;++k){const Real x=corners[k].x-d.centre.x,y=corners[k].y-d.centre.y,ct=(x*d.tangent.x+y*d.tangent.y)/d.D,cn=(x*d.normal.x+y*d.normal.y)/d.D;
        out.cornerImplicit[k]=jet(J(d.selectedSign/d.h)*(J(cn)-q0-q1*J(ct)-q2*J(ct*ct)));}
    for(int i=0;i<2;++i){const J m=q1+J(2)*q2*a[i];out.portTangent[i][0]=jet(J(d.tangent.x)+J(d.normal.x)*m);out.portTangent[i][1]=jet(J(d.tangent.y)+J(d.normal.y)*m);}
    for(int axis=0;axis<2;++axis)if(d.stationary[axis]){Q::StationaryFaceBoundaryJet event;event.current=*d.stationary[axis];
        const Real N=axis?d.normal.y:d.normal.x,T=axis?d.tangent.y:d.tangent.x,O=axis?d.centre.y:d.centre.x;
        const J parameter=-(J(T)+J(N)*q1)/(J(2*N)*q2),position=J(O)+J(T)*parameter+J(N)*(q0+q1*parameter+q2*parameter*parameter);
        const J lower(axis?d.input.box.lower.y:d.input.box.lower.x),upper(axis?d.input.box.upper.y:d.input.box.upper.x);
        for(int side=0;side<2;++side)event.slack[side]=jet(J(event.current.branchSign[side]/d.h)*(side?upper-position:position-lower));out.stationaryFace[axis]=event;}
    return out;
}
long double pairedAreaDifference(const Trial&before,const S::MomentValues&m0,const Trial&after,const S::MomentValues&m1)
{
    require(before.data_&&after.data_,"unprepared paired native packet");check(before,*before.data_);check(after,*after.data_);checkMoments(m0);checkMoments(m1);
    require(before.moment.beginBoundary==after.moment.beginBoundary&&before.moment.endBoundary==after.moment.endBoundary,"paired native trial changed actual patch sides");
    const Real a=before.moment.scale,b=after.moment.scale,j0=m0.moments[0],j1=m1.moments[0];
    const Real delta=(b-a)*(j1+j0)/2+(b+a)*(j1-j0)/2;require(std::isfinite(delta),"nonfinite paired native area difference");return delta;
}
}
