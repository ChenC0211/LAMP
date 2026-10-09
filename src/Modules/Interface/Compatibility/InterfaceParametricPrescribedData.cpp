#include "Modules/Interface/Compatibility/InterfacePrescribedData.h"
#include "Modules/Interface/Compatibility/Detail/InterfaceParametricInitialCurve.h"
#include "Core/Geometry/LevelSetGaussTrace.h"
#include <cmath>

namespace InterfacePrescribedData {
InterfaceParametricProjection::Curve parametricSupplyCap(double radius,double height,int segments,double bottom){
    using namespace Quantum::Geometry::Parametric;
    check(std::isfinite(radius)&&std::isfinite(height)&&std::isfinite(bottom)&&std::isfinite(bottom+height)&&radius>0&&height>0&&segments>0,"invalid supplied cap");
    Curve curve;curve.boundary=CurveBoundary::suppliedMeridian({radius,bottom},{0,bottom});
    Arc shape;shape.a=0;shape.b=1;shape.p[0]={radius,bottom};shape.p[1]={-radius,2*height};shape.p[2]={0,-height};
    for(int k=0;k<segments;k++)curve.arcs.push_back(restrictedArc(shape,double(k)/segments,double(k+1)/segments));return curve;
}
InterfaceParametricTransport::Flow ReversibleVortex::parametricSpatialFlow() const {
    const auto field=*this;
    return [field](Vector2D point,double duration) {
        const auto query=[&](Vector2D p){const auto j=field.sampleJet(p,field.profile_==Profile::SlipReversal?.25*field.period_:0);return LevelSetGaussTrace::VelocityJet{j.value,{j.derivative.xx,j.derivative.yx},{j.derivative.xy,j.derivative.yy}};};
        const double lower=field.profile_==Profile::SlipReversal?-.5:0;
        const auto map=LevelSetGaussTrace::advance(query,[lower](Vector2D,Vector2D){return LevelSetGaussTrace::SmoothPatch{{lower,lower},{lower+1,lower+1}};},point,duration,8);
        return InterfaceParametricTransport::MapJet{map.point,{map.derivative(0,0),map.derivative(1,0)},
            {map.derivative(0,1),map.derivative(1,1)},map.timeDerivative,map.displacement,map.coordinateResolution};
    };
}
double ReversibleVortex::parameterIncrement(double begin,double end) const {
    if(!(std::isfinite(begin) && std::isfinite(end) && begin>=0 && end>begin && end<=period_) || (begin<.5*period_ && end>.5*period_))
        throw std::invalid_argument("reversible parametric interval must stay within one temporal direction");
    const double pi=std::acos(-1.);
    if(profile_==Profile::SlipReversal)return period_/pi*std::sin(pi*(begin+end)/period_)*std::sin(pi*(end-begin)/period_);
    return 2*period_/pi*std::cos(pi*(begin+end)/(2*period_))*std::sin(pi*(end-begin)/(2*period_));
}
InterfaceParametricProjection::Curve parametricEllipse(Vector2D center,Vector2D radii,int segments) {
    using namespace Quantum::Geometry::Parametric;
    if(!(radii.x>0 && radii.y>0 && std::isfinite(radii.x) && std::isfinite(radii.y) && std::isfinite(center.x) && std::isfinite(center.y) && segments>=8 && segments<=65536))
        throw std::invalid_argument("invalid parametric initial ellipse");
    Curve curve;const double pi=std::acos(-1.),angle=2*pi/segments;
    for(int i=0;i<segments;++i) {
        const double a=double(i)/segments,b=double(i+1)/segments,t0=2*pi*a,t1=i+1==segments?0:2*pi*b,mid=pi*(a+b);
        const Point p0{center.x+radii.x*std::cos(t0),center.y+radii.y*std::sin(t0)};
        const Point delta{-2*radii.x*std::sin(mid)*std::sin(.5*angle),2*radii.y*std::cos(mid)*std::sin(.5*angle)};
        const Point d0{-angle*radii.x*std::sin(t0),angle*radii.y*std::cos(t0)},d1{-angle*radii.x*std::sin(t1),angle*radii.y*std::cos(t1)};
        Arc arc{a,b,{{p0,d0,delta*3-d0*2-d1,delta*(-2)+d0+d1,{}}}};
        double area=0;for(int q=0;q<5;++q){const double t=.5*(1+gx[q]);area+=.25*(b-a)*gw[q]*cross(arc.sampleFraction(t)-Point{center.x,center.y},arc.sampleFraction(t,1));}
        const Point tangent=arc.p[1]/30+arc.p[2]/30+arc.p[3]/35,normal{tangent.y,-tangent.x};
        const Point bubble=normal*((.5*radii.x*radii.y*angle-area)/dot(normal,normal));
        arc.p[2]+=bubble;arc.p[3]+=bubble*(-2);arc.p[4]=bubble;curve.arcs.push_back(arc);
    }
    return curve;
}
InterfaceParametricProjection::Curve parametricAnalytic(const Quantum::ScalarExpression& first,
    const Quantum::ScalarExpression& second,int segments,bool rz){
    using namespace Quantum::Geometry::Parametric;
    return detail::fromInitialJet([&](double s){
        const auto x=first.sample(s),y=second.sample(s);
        check(std::hypot(x.derivative,y.derivative)>0,"analytic curve has a zero tangent");
        return std::array<Point,2>{{{x.value,y.value},{x.derivative,y.derivative}}};
    },segments,rz);
}
InterfaceParametricProjection::Curve detail::fromInitialJet(const InitialMap& jet,int segments,bool rz){
    using namespace Quantum::Geometry::Parametric;
    check(segments>=8&&segments<=65536,"analytic curve requires 8..65536 segments");
    auto begin=jet(0),end=jet(1);
    const double scale=std::max({norm(begin[0]),norm(end[0]),norm(begin[1]),norm(end[1])});
    const double tolerance=256*std::numeric_limits<double>::epsilon()*scale;
    Curve curve;
    if(rz){
        check(std::abs(begin[0].x)<=tolerance&&std::abs(end[0].x)<=tolerance&&
            std::abs(begin[1].y)<=tolerance&&std::abs(end[1].y)<=tolerance&&
            begin[1].x>0&&end[1].x<0&&begin[0].y<end[0].y,
            "analytic RZ curve requires regular lower/upper axis poles: r=0, z'=0, r'(0)>0, r'(1)<0");
        begin[0].x=end[0].x=0;begin[1].y=end[1].y=0;curve.boundary=CurveBoundary::axisMeridian();
    }else {
        check(norm(begin[0]-end[0])<=tolerance&&norm(begin[1]-end[1])<=tolerance,
            "analytic XY curve must be periodic in position and first derivative");
        end=begin;
    }
    for(int i=0;i<segments;i++){
        const double a=double(i)/segments,b=double(i+1)/segments,h=b-a;
        const auto left=i==0?begin:jet(a),right=i+1==segments?end:jet(b);
        const auto delta=right[0]-left[0],d0=left[1]*h,d1=right[1]*h;
        Arc arc{a,b,{{left[0],d0,delta*3-d0*2-d1,delta*(-2)+d0+d1,{}}}};
        const auto bubble=(jet(.5*(a+b))[0]-arc.sampleFraction(.5))*16;
        arc.p[2]+=bubble;arc.p[3]+=bubble*(-2);arc.p[4]=bubble;curve.arcs.push_back(arc);
    }
    double signedArea=0;
    for(const auto& arc:curve.arcs)for(int q=0;q<5;q++){
        const double t=.5*(1+gx[q]);
        signedArea+=.25*(arc.b-arc.a)*gw[q]*cross(arc.sampleFraction(t)-begin[0],arc.sampleFraction(t,1));
    }
    check(std::isfinite(signedArea)&&signedArea>0,"analytic curve must enclose a positive area with counterclockwise orientation");
    return curve;
}
InterfaceParametricProjection::Curve parametricNeckedSpheroid(double centerZ,Vector2D radii,double ratio,int segments){
    using namespace Quantum::Geometry::Parametric;
    check(std::isfinite(centerZ)&&std::isfinite(radii.x)&&std::isfinite(radii.y)&&radii.x>0&&radii.y>0&&std::isfinite(ratio)&&ratio>0&&ratio<2./3&&segments>=8&&segments<=65536,"invalid necked spheroid");
    Curve curve;curve.boundary=CurveBoundary::axisMeridian();const double pi=std::acos(-1.),amplitude=1-ratio;
    const auto position=[&](double s){const double q=std::sin(pi*s);return Point{s==0||s==1?0:radii.x*q*(1-amplitude*q*q),centerZ-radii.y*std::cos(pi*s)};};
    const auto derivative=[&](double s){const double q=std::sin(pi*s);return Point{pi*radii.x*std::cos(pi*s)*(1-3*amplitude*q*q),s==0||s==1?0:radii.y*pi*q};};
    for(int k=0;k<segments;k++){const double a=double(k)/segments,b=double(k+1)/segments,h=b-a;const auto p=position(a),delta=position(b)-p,d0=derivative(a)*h,d1=derivative(b)*h;
        Arc arc{a,b,{{p,d0,delta*3-d0*2-d1,delta*(-2)+d0+d1,{}}}};const auto bubble=(position(.5*(a+b))-arc.sampleFraction(.5))*16;
        arc.p[2]+=bubble;arc.p[3]+=bubble*(-2);arc.p[4]=bubble;curve.arcs.push_back(arc);}return curve;
}
InterfaceParametricProjection::Curve parametricSpheroid(double centerZ,Vector2D radii,int segments) {
    using namespace Quantum::Geometry::Parametric;
    if(!(std::isfinite(centerZ)&&radii.x>0&&radii.y>0&&std::isfinite(radii.x)&&std::isfinite(radii.y)&&segments>=8&&segments<=65536))
        throw std::invalid_argument("invalid parametric RZ spheroid");
    Curve curve;curve.boundary=CurveBoundary::axisMeridian();const double pi=std::acos(-1.),angle=pi/segments;
    const auto position=[&](int i){return Point{i==0||i==segments?0:radii.x*std::sin(pi*i/segments),centerZ-radii.y*std::cos(pi*i/segments)};};
    const auto tangent=[&](int i){return Point{angle*radii.x*std::cos(pi*i/segments),i==0||i==segments?0:angle*radii.y*std::sin(pi*i/segments)};};
    for(int i=0;i<segments;i++){
        const double a=double(i)/segments,b=double(i+1)/segments;const Point p=position(i),delta=position(i+1)-p,d0=tangent(i),d1=tangent(i+1);
        Arc arc{a,b,{{p,d0,delta*3-d0*2-d1,delta*(-2)+d0+d1,{}}}};
        const double mid=pi*(i+.5)/segments;const Point center{radii.x*std::sin(mid),centerZ-radii.y*std::cos(mid)};
        const Point bubble=(center-arc.sampleFraction(.5))*16;arc.p[2]+=bubble;arc.p[3]+=bubble*(-2);arc.p[4]=bubble;curve.arcs.push_back(arc);
    }
    return curve;
}
} // namespace InterfacePrescribedData
