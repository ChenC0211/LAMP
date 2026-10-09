#pragma once
#include "Core/Geometry/LevelSetBoundaryTrace.h"

namespace LevelSetBoundaryTrace
{
// One classical RK4 characteristic of the SAME linearly interpolated two-time
// velocity, with its exact discrete tangent. Complete closed rectangular
// domain only: every stage is checked; an exterior stage is rejected. This
// does not invent an inflow history or reuse an interior tangent at an event.
// Fourth order describes the resolved smooth interpolated ODE, not the time
// accuracy of the supplied Hydro fields or arbitrary nonsmooth callbacks.
inline JetResult traceClosedRk4WithJacobian(const Domain& domain,Vector2D point,double dt,
    const DifferentiableVelocity& current,const DifferentiableVelocity& previous)
{
    using namespace Detail;
    require(finite(domain.minCoord) && finite(domain.maxCoord) && finite(domain.cellSize) &&
        domain.minCoord.x<domain.maxCoord.x && domain.minCoord.y<domain.maxCoord.y &&
        domain.cellSize.x>0 && domain.cellSize.y>0 && std::isfinite(dt) && dt>=0 && current && previous,
        "invalid closed RK4 domain or time");
    require(!domain.radialAxisSymmetry || domain.minCoord.x==0,"closed RK4 axis requires r=0");
    const auto inside=[&](Vector2D p) {return finite(p) && p.x>=domain.minCoord.x && p.x<=domain.maxCoord.x &&
        p.y>=domain.minCoord.y && p.y<=domain.maxCoord.y;};
    require(inside(point),"closed RK4 starts outside domain");
    require(!domain.activeSegment || domain.activeSegment(point,point),"closed RK4 starts outside active coverage");
    const double h=std::min(domain.cellSize.x,domain.cellSize.y);
    const auto valid=[&](Vector2D p,const VelocityJet& v) {
        require(finite(v.value) && finite(v.derivative) && std::hypot(v.value.x,v.value.y)*dt<=.5*h,
            "invalid closed RK4 velocity/Jacobian or half-cell CFL exceeded");
        require(!domain.radialAxisSymmetry || p.x!=0 || (v.value.x==0 && v.derivative.xy==0),
            "closed RK4 radial velocity parity violation");
    };
    const auto sample=[&](Vector2D p,double fraction) {
        require(inside(p),"closed RK4 stage crosses physical boundary; boundary history not enabled");
        require(!domain.activeSegment || domain.activeSegment(point,p),"closed RK4 stage crosses inactive coverage");
        if(fraction==1){auto v=current(p);valid(p,v);return v;}
        if(fraction==0){auto v=previous(p);valid(p,v);return v;}
        const auto a=current(p),b=previous(p);valid(p,a);valid(p,b);
        VelocityJet v{(a.value+b.value)*.5,{.5*(a.derivative.xx+b.derivative.xx),.5*(a.derivative.xy+b.derivative.xy),
            .5*(a.derivative.yx+b.derivative.yx),.5*(a.derivative.yy+b.derivative.yy)}};
        valid(p,v);return v;
    };
    if(dt==0)return {{point,0,std::nullopt,false},Jacobian::identity(),0};
    const auto a=sample(point,1);
    const auto b=sample(point-a.value*(dt*.5),.5);
    const auto db=multiply(b.derivative,identityMinus(dt*.5,a.derivative));
    const auto c=sample(point-b.value*(dt*.5),.5);
    const auto dc=multiply(c.derivative,identityMinus(dt*.5,db));
    const auto d=sample(point-c.value*dt,0);
    const auto dd=multiply(d.derivative,identityMinus(dt,dc));
    const Vector2D result=point-(a.value+b.value*2+c.value*2+d.value)*(dt/6);
    const Jacobian sum{a.derivative.xx+2*db.xx+2*dc.xx+dd.xx,a.derivative.xy+2*db.xy+2*dc.xy+dd.xy,
        a.derivative.yx+2*db.yx+2*dc.yx+dd.yx,a.derivative.yy+2*db.yy+2*dc.yy+dd.yy};
    const auto tangent=identityMinus(dt/6,sum);
    require(inside(result) && finite(tangent),"closed RK4 endpoint leaves domain or tangent is nonfinite");
    require(!domain.activeSegment || domain.activeSegment(point,result),"closed RK4 endpoint crosses inactive coverage");
    return {{result,dt,std::nullopt,false},tangent,1};
}
}
