#include "Modules/Interface/Compatibility/InterfacePrescribedData.h"
#include "Modules/Interface/Compatibility/Detail/InterfaceParametricInitialCurve.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>

namespace InterfacePrescribedData {
namespace {
namespace P=Quantum::Geometry::Parametric;
using Point=P::Point;
struct FieldPoint {Point point,normal,tangent;double distance;};
struct TracePoint {FieldPoint field;double length;};
class ImplicitTrace {
public:
    ImplicitTrace(const Quantum::ScalarExpression& field,const P::DomainRectangle& box,bool rz,double tolerance):
        field_(field),box_(box),rz_(rz),tolerance_(tolerance),maximumStep_(std::min(box.right-box.left,box.top-box.bottom)/32){
        box.validate();
        const double scale=std::max({std::abs(box.left),std::abs(box.right),std::abs(box.bottom),std::abs(box.top),maximumStep_});
        P::check(std::isfinite(tolerance)&&tolerance>=256*std::numeric_limits<double>::epsilon()*scale&&tolerance<maximumStep_/64,
            "implicit traceTolerance must resolve the domain and exceed coordinate roundoff");
        P::check(!rz||box.left==0,"implicit RZ trace requires an r=0 domain boundary");
    }
    void trace(Point seed){
        P::check(box_.contains(seed),"implicit seed must be inside Domain");
        P::check(!rz_||seed.x==0,"implicit RZ seed must lie on the lower axis pole (r=0)");
        auto first=project(seed,maximumStep_);
        if(rz_)first=axisPoint(first.point,true);
        else P::check(interior(first.point),"implicit XY seed must lie strictly inside Domain");
        points_.push_back({first,0});double step=maximumStep_;
        for(int count=0;count<20000;count++){
            const auto old=points_.back();
            auto accepted=advance(old.field.point,step);const auto next=accepted.first;
            const double used=accepted.second;
            P::check(old.length+used>old.length,"implicit continuation made no arclength progress");
            if(rz_&&next.point.x<=0){
                P::check(old.field.point.x>0,"implicit meridian returned to its starting pole");
                auto event=intersect(old.field.point,used,[](Point p){return p.x;});
                const auto pole=axisPoint(event.first.point,false);
                P::check(pole.point.y>first.point.y,"implicit meridian must end at an upper pole");
                points_.push_back({pole,old.length+event.second});return;
            }
            if(!rz_&&points_.size()>5&&P::dot(old.field.point-first.point,first.tangent)<0&&
                P::dot(next.point-first.point,first.tangent)>=0&&P::dot(next.tangent,first.tangent)>.5){
                auto event=intersect(old.field.point,used,[&](Point p){return P::dot(p-first.point,first.tangent);});
                if(P::norm(event.first.point-first.point)<=64*tolerance_&&P::dot(event.first.tangent,first.tangent)>.5){
                    points_.push_back({first,old.length+event.second});return;
                }
            }
            P::check(interior(next.point),"implicit component reaches an undeclared physical boundary");
            points_.push_back({next,old.length+used});step=std::min(maximumStep_,1.5*used);
        }
        throw std::runtime_error("implicit component did not close within 20000 continuation steps");
    }
    detail::InitialJet sample(double s) const {
        const double total=points_.back().length;
        if(s==0)return {points_.front().field.point,points_.front().field.tangent*total};
        if(s==1)return {points_.back().field.point,points_.back().field.tangent*total};
        const double target=s*total;
        auto right=std::upper_bound(points_.begin(),points_.end(),target,[](double length,const TracePoint& p){return length<p.length;});
        P::check(right!=points_.begin()&&right!=points_.end(),"invalid implicit resampling interval");
        const auto& left=*(right-1);const double h=right->length-left.length,u=(target-left.length)/h;
        const auto p=left.field.point,delta=right->field.point-p,d0=left.field.tangent*h,d1=right->field.tangent*h;
        const auto guess=p+d0*u+(delta*3-d0*2-d1)*(u*u)+(delta*(-2)+d0+d1)*(u*u*u);
        const auto value=project(guess,maximumStep_);
        return {value.point,value.tangent*total};
    }
private:
    const Quantum::ScalarExpression& field_;
    P::DomainRectangle box_;
    bool rz_;
    double tolerance_,maximumStep_;
    std::vector<TracePoint> points_;
    bool interior(Point p)const{return p.x>box_.left&&p.x<box_.right&&p.y>box_.bottom&&p.y<box_.top;}
    FieldPoint evaluate(Point p)const {
        const auto j=field_.sample(p.x,p.y);const double norm=std::hypot(j.derivative,j.partialSecond);
        if(!(norm>0)||!std::isfinite(norm))throw std::domain_error("implicit zero-set continuation requires a finite nonzero gradient");
        const Point normal{j.derivative/norm,j.partialSecond/norm};const double distance=j.value/norm;
        if(!std::isfinite(distance))throw std::domain_error("implicit zero-set residual is not finite");
        return {p,normal,{-normal.y,normal.x},distance};
    }
    FieldPoint project(Point seed,double limit)const {
        auto value=evaluate(seed);
        for(int i=0;i<24;i++){
            if(std::abs(value.distance)<=tolerance_)return value;
            const auto next=value.point-value.normal*value.distance;
            if(P::norm(next-seed)>limit)throw std::domain_error("implicit seed/corrector left its local zero-set neighborhood");
            value=evaluate(next);
        }
        throw std::domain_error("implicit zero-set corrector did not converge");
    }
    FieldPoint axisPoint(Point p,bool lower)const {
        p.x=0;
        for(int i=0;i<24;i++){
            auto f=evaluate(p);
            if(std::abs(f.distance)<=tolerance_){
                P::check(std::abs(f.normal.x)<=1e-10&&(lower?f.normal.y<0:f.normal.y>0)&&p.y>box_.bottom&&p.y<box_.top,
                    "implicit RZ pole must be regular with F_r=0 and the declared lower/upper orientation");
                f.normal={0,lower?-1.:1.};f.tangent={lower?1.:-1.,0};return f;
            }
            P::check(std::abs(f.normal.y)>.5,"implicit RZ pole cannot be corrected along the axis");
            p.y-=f.distance/f.normal.y;
        }
        throw std::runtime_error("implicit RZ axis root did not converge");
    }
    Point rk4(Point p,double h)const {
        const auto a=evaluate(p).tangent,b=evaluate(p+a*(.5*h)).tangent;
        const auto c=evaluate(p+b*(.5*h)).tangent,d=evaluate(p+c*h).tangent;
        if(std::min({P::dot(a,b),P::dot(a,c),P::dot(a,d)})<.5)
            throw std::domain_error("implicit continuation step turns too far");
        return p+(a+b*2+c*2+d)*(h/6);
    }
    // RK4 step doubling and normal correction control the initial trace;
    // this is not a second physical-time interface transport scheme.
    std::optional<FieldPoint> trial(Point p,double h)const {
        try{
            const auto full=rk4(p,h),half=rk4(rk4(p,.5*h),.5*h);
            const auto projected=project(half,.5*h);
            const double error=P::norm(full-half)/15+P::norm(projected.point-half);
            if(error<=tolerance_)return projected;
        }catch(const std::domain_error&){} // Shorten a geometric continuation step only.
        return {};
    }
    std::pair<FieldPoint,double> advance(Point p,double h)const {
        for(int i=0;i<60;i++){
            if(auto value=trial(p,h))return {*value,h};
            h*=.5;if(h<8*tolerance_)break;
        }
        throw std::runtime_error("implicit continuation cannot resolve a regular local branch; check the formula, seed and traceTolerance");
    }
    std::pair<FieldPoint,double> intersect(Point p,double h,const std::function<double(Point)>& section)const {
        double lo=0,hi=h;const double initial=section(p);
        auto endpoint=trial(p,h);P::check(bool(endpoint),"implicit event lost its accepted step");
        P::check(initial*section(endpoint->point)<=0,"implicit closure event is not bracketed");
        for(int i=0;i<60;i++){
            const double mid=.5*(lo+hi);auto value=trial(p,mid);
            P::check(bool(value),"implicit event substep cannot be resolved");
            const double residual=section(value->point);
            if(std::abs(residual)<=tolerance_||hi-lo<=tolerance_)return {*value,mid};
            if((residual<0)==(initial<0))lo=mid;else hi=mid;
        }
        throw std::runtime_error("implicit closure event did not converge");
    }
};
}
InterfaceParametricProjection::Curve parametricImplicit(const Quantum::ScalarExpression& field,
    const P::DomainRectangle& domain,Vector2D seed,int segments,bool rz,double tolerance){
    ImplicitTrace trace(field,domain,rz,tolerance);trace.trace({seed.x,seed.y});
    auto curve=detail::fromInitialJet([&](double s){return trace.sample(s);},segments,rz);
    P::validateContainedSurface(curve,domain,rz);return curve;
}
}
