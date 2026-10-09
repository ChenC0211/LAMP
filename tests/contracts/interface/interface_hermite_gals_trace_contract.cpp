#include "Modules/Interface/HermiteGalsTrace.h"

#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace
{
using namespace Quantum::Interface;
using Matrix = LevelSetBoundaryTrace::Jacobian;
using Jet = LevelSetBoundaryTrace::VelocityJet;
using Query = HermiteGalsVelocityView::Query;
using Domain = LevelSetBoundaryTrace::Domain;
using Position = HermiteGalsTimeLayer::Position;
static_assert(!std::is_default_constructible_v<HermiteGalsTraceResult>, "no fabricated empty trace result");
static_assert(std::is_same_v<decltype(std::declval<HermiteGalsTraceResult&>().jacobian()),
    const std::optional<Matrix>&>, "published tangent is read-only");

void check(bool ok, const char* reason) { if (!ok) throw std::runtime_error(reason); }
void close(double actual, double expected, double tolerance, const char* reason)
{ check(std::isfinite(actual) && std::abs(actual - expected) <= tolerance * (1 + std::abs(expected)), reason); }
template<class F> void rejects(F action, const char* fragment)
{
    try { action(); }
    catch (const std::exception& error) {
        check(std::string(error.what()).find("Interface/HermiteGalsTrace") != std::string::npos,
            "trace failure lacks module/stage context");
        check(std::string(error.what()).find(fragment) != std::string::npos, "unexpected rejection reason");
        return;
    }
    throw std::runtime_error(std::string("missing rejection: ") + fragment);
}
double error(Matrix a, Matrix b)
{ return std::max({std::abs(a.xx-b.xx),std::abs(a.xy-b.xy),std::abs(a.yx-b.yx),std::abs(a.yy-b.yy)}); }
Matrix product(Matrix a, Matrix b)
{ return {a.xx*b.xx+a.xy*b.yx,a.xx*b.xy+a.xy*b.yy,a.yx*b.xx+a.yy*b.yx,a.yx*b.xy+a.yy*b.yy}; }
double norm(Vector2D p) { return std::hypot(p.x, p.y); }
template<class Map> Matrix numericalDerivative(Map map, Vector2D p, double h)
{
    const auto xp = map(Vector2D{p.x+h,p.y}), xm = map(Vector2D{p.x-h,p.y});
    const auto yp = map(Vector2D{p.x,p.y+h}), ym = map(Vector2D{p.x,p.y-h});
    return {(xp.x-xm.x)/(2*h),(yp.x-ym.x)/(2*h),(xp.y-xm.y)/(2*h),(yp.y-ym.y)/(2*h)};
}

// Standalone analytic velocity tests need legitimate I01 revision tokens,
// not a production geometry publisher. No Interface transport is performed.
struct Context
{
    std::unique_ptr<QuadtreeMesh> mesh;
    HermiteGalsIdentity oldId, newId;
    Domain domain{{0,0},{1,1},{.25,.25},false};
    Context()
    {
        DictionaryMesh input; input.Nx=input.Ny=1; input.maxCellNum=16; input.maxInitialIterations=0;
        mesh = std::make_unique<QuadtreeMesh>(input);
        HermiteGalsPrimaryState data; data.phaseVolumes={{0,.5}};
        HermiteGalsCandidate a(*mesh, {}, data), b(*mesh, {}, data);
        oldId=a.identity(); newId=b.identity();
    }
    HermiteGalsVelocityView old(double t, Query query) const
    { return {oldId,{0,t,t,Position::Initial},std::move(query)}; }
    HermiteGalsVelocityView next(double a, double b, Query query) const
    { return {newId,{1,a,b,Position::End},std::move(query)}; }
};

void affineAndTime(const Context& c)
{
    const double begin=1.25, end=1.45, dt=end-begin, shear=.3;
    const auto field=[=](double t) -> Query {
        return [=](const Vector2D& p) {return Jet{{shear*p.y+.02*t,0},{0,shear,0,0}};};
    };
    const Vector2D p{.61,.37};
    const auto a=c.old(begin,field(begin)), b=c.next(begin,end,field(end));
    const auto result=traceHermiteGalsWithJacobian(c.domain,p,a,b);
    close(result.departurePoint().x,p.x-dt*(shear*p.y+.01*(begin+end)),2e-14,"affine shear map");
    close(result.departurePoint().y,p.y,2e-14,"affine shear y");
    check(error(result.interiorJacobian(),{1,-dt*shear,0,1})<2e-14,"affine shear Jacobian orientation");
    check(result.previousGeometryIdentity()==c.oldId && result.arrivalGeometryIdentity()==c.newId &&
        result.previousTimeLayer().time()==begin && result.arrivalTimeLayer().time()==end &&
        result.departureTime()==begin && result.acceptedSubsteps()==1,"trace lost endpoint provenance");

    const HermiteGalsTimeLayer midpoint{2,1.3,1.7,Position::Midpoint};
    const double mid=midpoint.time();
    const auto translation=[](double t)->Query {return [=](const Vector2D&){return Jet{{.02*t,0},{}};};};
    const auto old=c.old(midpoint.begin,translation(midpoint.begin));
    const HermiteGalsVelocityView middle(c.newId,midpoint,translation(mid));
    const auto half=traceHermiteGalsWithJacobian(c.domain,p,old,middle);
    close(half.departurePoint().x,p.x-.01*(mid*mid-midpoint.begin*midpoint.begin),2e-14,"midpoint physical interval");
    check(half.elapsed()==mid-midpoint.begin && half.arrivalTimeLayer().position==Position::Midpoint,
        "midpoint incorrectly traced the full step");

    int calls=0;
    const Query unused=[&](const Vector2D&){++calls;return Jet{{99,99},{99,99,99,99}};};
    auto z=c.old(5,unused);
    const auto zero=traceHermiteGalsWithJacobian(c.domain,p,z,z);
    check(calls==0 && zero.departurePoint().x==p.x && zero.departurePoint().y==p.y &&
        zero.departureTime()==5 && zero.acceptedSubsteps()==0 && error(zero.interiorJacobian(),Matrix::identity())==0,
        "zero interval is not exact unsampled identity");
    const Query stopped=[](const Vector2D&){return Jet{{0,0},{}};};
    auto stoppedOld=c.old(0,stopped), stoppedNew=c.next(0,.1,stopped);
    const auto stationary=traceHermiteGalsWithJacobian(c.domain,p,stoppedOld,stoppedNew);
    check(stationary.departurePoint().x==p.x && stationary.departurePoint().y==p.y &&
        error(stationary.interiorJacobian(),Matrix::identity())==0,"stationary map drift");
    std::cout<<"affine/time/zero: analytic B/J and endpoint provenance PASS\n";
}

Query nonlinear(double t)
{
    return [=](const Vector2D& p) {
        const double a=.08+.02*t,b=.06-.01*t,s=1.3*p.x+.4*p.y,r=.7*p.x-1.1*p.y;
        return Jet{{a*std::sin(s)+.03*p.y*p.y,b*std::cos(r)+.02*p.x*p.y},
            {1.3*a*std::cos(s),.4*a*std::cos(s)+.06*p.y,
             -.7*b*std::sin(r)+.02*p.y,1.1*b*std::sin(r)+.02*p.x}};
    };
}
void differenceSweep(const Context& c)
{
    for(double duration:{.037,.113,.21}) {
        const double begin=.7,end=begin+duration;
        const auto a=c.old(begin,nonlinear(begin)),b=c.next(begin,end,nonlinear(end));
        double previous=0,lastOrder=0,finest=0;
        for(double h:{.02,.01,.005,.0025,.00125,.000625}) {
            double maximum=0;
            for(double x:{.23,.41,.68,.79}) for(double y:{.24,.46,.71}) {
                const Vector2D p{x,y};
                const auto reference=traceHermiteGalsWithJacobian(c.domain,p,a,b);
                const auto numeric=numericalDerivative([&](Vector2D q) {
                    const auto result=traceHermiteGalsWithJacobian(c.domain,q,a,b);
                    check(!result.boundaryHit() && result.acceptedSubsteps()==reference.acceptedSubsteps(),
                        "finite differences crossed a trace branch");
                    return result.departurePoint();
                },p,h);
                maximum=std::max(maximum,error(reference.interiorJacobian(),numeric));
            }
            lastOrder=previous>0?std::log2(previous/maximum):0;
            std::cout<<"FD,"<<duration<<','<<h<<','<<maximum<<','<<lastOrder<<'\n';
            previous=maximum;finest=maximum;
        }
        check(finest<2e-8 && lastOrder>=1.85 && lastOrder<=2.15,
            "nonlinear four-component finite differences did not converge quadratically to the discrete tangent");
    }
}

void substepWitness(const Context& c)
{
    std::vector<std::array<double,3>> valueQueries,jetQueries;
    const auto value=[](Vector2D p,bool now) {
        return now?Vector2D{.04*p.x+.03*p.y+.3,.02*p.x-.01*p.y}:
            Vector2D{-.02*p.x+.04*p.y-.3,-.03*p.x+.01*p.y};
    };
    const auto query=[&](bool now)->Query {return [&,now](const Vector2D& p) {
        jetQueries.push_back({p.x,p.y,now?1.0:0.0});
        return Jet{value(p,now),now?Matrix{.04,.03,.02,-.01}:Matrix{-.02,.04,-.03,.01}};
    };};
    const LevelSetBoundaryTrace::Velocity oldValue=[&](const Vector2D& p) {
        valueQueries.push_back({p.x,p.y,0});return value(p,false);};
    const LevelSetBoundaryTrace::Velocity newValue=[&](const Vector2D& p) {
        valueQueries.push_back({p.x,p.y,1});return value(p,true);};
    const auto a=c.old(0,query(false)),b=c.next(0,.1,query(true));
    const Vector2D p{.01,.45};
    const auto primal=LevelSetBoundaryTrace::trace(c.domain,p,.1,newValue,oldValue);
    const auto tangent=traceHermiteGalsWithJacobian(c.domain,p,a,b);
    check(valueQueries==jetQueries && primal.point.x==tangent.departurePoint().x &&
        primal.point.y==tangent.departurePoint().y && primal.elapsed==tangent.elapsed() &&
        primal.side==tangent.boundarySide() && tangent.acceptedSubsteps()==2,
        "tangent path changed the actual primal map, sampling order or selected substeps");
    const auto numeric=numericalDerivative([&](Vector2D q) {
        const auto jet=traceHermiteGalsWithJacobian(c.domain,q,a,b);
        check(jet.acceptedSubsteps()==2,"substep witness changed branch");
        return LevelSetBoundaryTrace::trace(c.domain,q,.1,newValue,oldValue).point;
    },p,1e-5);
    const double difference=error(tangent.interiorJacobian(),numeric);
    check(difference<2e-9,"composite tangent differs from independent primal differences");
    std::cout<<"SUBSTEP,2,"<<difference<<'\n';
}

void rotationConvergence(const Context& c)
{
    // omega(t)=.6+.4*t: supplied linear-in-time interpolant is EXACT here.
    // Total angle on [0,1] is .8. Matrix comparison includes all four entries.
    const auto velocity=[](double t)->Query { return [=](const Vector2D& p) {
        const double w=.6+.4*t;
        return Jet{{-w*(p.y-.5),w*(p.x-.5)},{0,-w,w,0}};
    };};
    const Vector2D start{.68,.43};
    const double cosine=std::cos(.8),sine=std::sin(.8);
    const Matrix exact{cosine,sine,-sine,cosine};
    const Vector2D expected{.5+cosine*(start.x-.5)+sine*(start.y-.5),
        .5-sine*(start.x-.5)+cosine*(start.y-.5)};
    for(bool alternating:{false,true}) {
        double previousPoint=0,previousMatrix=0;
        for(int count:{8,16,32,64}) {
            std::vector<double> times(static_cast<std::size_t>(count)+1,0.0);
            for(int i=0;i<count;++i)
                times[i+1]=times[i]+(alternating?(i%2==0?1.0:2.0)/(1.5*count):1.0/count);
            times.back()=1;
            auto point=start;auto jacobian=Matrix::identity();
            for(int i=count-1;i>=0;--i) {
                const HermiteGalsVelocityView old(c.oldId,{i,times[i],times[i],Position::Initial},velocity(times[i]));
                const HermiteGalsVelocityView now(c.newId,{i+1,times[i],times[i+1],Position::End},velocity(times[i+1]));
                const auto traced=traceHermiteGalsWithJacobian(c.domain,point,old,now);
                jacobian=product(traced.interiorJacobian(),jacobian);
                point=traced.departurePoint();
            }
            const double pointError=norm(point-expected),matrixError=error(jacobian,exact);
            const double pointOrder=previousPoint>0?std::log2(previousPoint/pointError):0;
            const double matrixOrder=previousMatrix>0?std::log2(previousMatrix/matrixError):0;
            std::cout<<"CONVERGENCE,"<<(alternating?"alternating":"uniform")<<','<<count<<','<<pointError<<','
                <<matrixError<<','<<pointOrder<<','<<matrixOrder<<'\n';
            if(count>=32) check(pointOrder>=1.85 && pointOrder<=2.15 && matrixOrder>=1.85 && matrixOrder<=2.15,
                "B/J failed global second-order time convergence");
            previousPoint=pointError;previousMatrix=matrixError;
        }
    }
}

void boundaryAxisAndFailures(const Context& c)
{
    const auto constant=[](Vector2D v)->Query {return [=](const Vector2D& p) {
        check(p.x>=0 && p.x<=1 && p.y>=0 && p.y<=1,"callback sampled outside physical domain");
        return Jet{v,{}};
    };};
    const std::array<Vector2D,4> points{{{.005,.5},{.5,.005},{.995,.5},{.5,.995}}};
    const std::array<Vector2D,4> speeds{{{.2,0},{0,.2},{-.2,0},{0,-.2}}};
    const std::array<Boundary,4> sides{{Boundary::west,Boundary::south,Boundary::east,Boundary::north}};
    for(std::size_t i=0;i<points.size();++i) {
        const auto a=c.old(0,constant(speeds[i])),b=c.next(0,.1,constant(speeds[i]));
        const auto hit=traceHermiteGalsWithJacobian(c.domain,points[i],a,b);
        check(hit.boundaryHit() && hit.boundarySide()==sides[i] && !hit.jacobian(),"boundary hit invented interior tangent");
        close(hit.elapsed(),.025,1e-12,"boundary hit time");
        close(hit.departureTime(),.075,1e-12,"physical boundary time");
        rejects([&]{hit.interiorJacobian();},"boundary-history tangent");
    }
    auto axis=c.domain;axis.radialAxisSymmetry=true;
    const Query radial=[](const Vector2D& p){return Jet{{.3*p.x,-.2*(p.y-.5)},{.3,0,0,-.2}};};
    const auto a=c.old(0,radial),b=c.next(0,.1,radial);
    const auto pole=traceHermiteGalsWithJacobian(axis,{0,.5},a,b);
    check(!pole.boundaryHit() && pole.departurePoint().x==0 && pole.departurePoint().y==.5,
        "regular axis point moved or was treated as inflow");
    close(pole.interiorJacobian().xx,1-.03+.00045,2e-14,"axis transverse tangent");
    close(pole.interiorJacobian().yy,1+.02+.0002,2e-14,"axis longitudinal tangent");
    check(error(pole.interiorJacobian(),Matrix::identity())>.02,"stationary point lost nontrivial derivative");
    const Query parity=[](const Vector2D&){return Jet{{0,0},{0,.1,0,0}};};
    const auto badAxis=c.next(0,.1,parity);
    rejects([&]{traceHermiteGalsWithJacobian(axis,{0,.5},a,badAxis);},"tangential derivative");
    const auto flowOld=c.old(0,constant({.2,0})),flowNew=c.next(0,.1,constant({.2,0}));
    rejects([&]{traceHermiteGalsWithJacobian(axis,{.005,.5},flowOld,flowNew);},"negative r");
    auto covered=c.domain;covered.activeSegment=[](Vector2D p,Vector2D q){return p.x>=.4 && q.x>=.4;};
    rejects([&]{traceHermiteGalsWithJacobian(covered,{.405,.5},flowOld,flowNew);},"selected fluid domain");
    auto tiny=c.domain;tiny.cellSize={.01,.01};
    rejects([&]{traceHermiteGalsWithJacobian(tiny,{.5,.5},flowOld,flowNew);},"half a cell");
    auto invalid=c.domain;invalid.cellSize={2,2};
    rejects([&]{traceHermiteGalsWithJacobian(invalid,{.5,.5},flowOld,flowNew);},"cell size exceeds");
    rejects([&]{traceHermiteGalsWithJacobian(c.domain,{-1,.5},flowOld,flowNew);},"initial point");
    const Query badValue=[](const Vector2D&){return Jet{{std::numeric_limits<double>::quiet_NaN(),0},{}};};
    const Query badDerivative=[](const Vector2D&){return Jet{{0,0},{0,0,std::numeric_limits<double>::infinity(),0}};};
    const auto v=c.next(0,.1,badValue),g=c.next(0,.1,badDerivative);
    rejects([&]{traceHermiteGalsWithJacobian(c.domain,{.5,.5},flowOld,v);},"nonfinite sampled velocity");
    rejects([&]{traceHermiteGalsWithJacobian(c.domain,{.5,.5},flowOld,g);},"velocity Jacobian");
    rejects([&]{HermiteGalsVelocityView view({}, {}, constant({0,0}));},"unbound geometry");
    rejects([&]{HermiteGalsVelocityView view(c.oldId, {}, {});},"missing paired");
    rejects([&]{HermiteGalsVelocityView view(c.oldId, {0,0,1,Position::Initial}, radial);},"initial/restart layer");
    const auto wrongTime=c.old(.01,radial);
    rejects([&]{traceHermiteGalsWithJacobian(c.domain,{.5,.5},wrongTime,b);},"previous velocity time");
    const double large=1e16, adjacent=std::nextafter(large,std::numeric_limits<double>::infinity());
    const auto largeOld=c.old(large,radial);
    const HermiteGalsVelocityView unresolved(c.newId,{1,large,adjacent,Position::Midpoint},radial);
    rejects([&]{traceHermiteGalsWithJacobian(c.domain,{.5,.5},largeOld,unresolved);},"unresolved in floating point");
    std::cout<<"events/axis/failures: physical times, no invented tangent, domain/CFL and invalid input guards PASS\n";
}
void roundoffMidpointEvents()
{
    const double eps=std::numeric_limits<double>::epsilon(),eta=std::numeric_limits<double>::denorm_min();
    const std::array<Boundary,4> sides{{Boundary::west,Boundary::south,Boundary::east,Boundary::north}};
    int cases=0,zeroTimes=0,guards=0;std::size_t maximumQueries=0;
    const auto inside=[](const Domain& d,Vector2D p){return p.x>=d.minCoord.x && p.x<=d.maxCoord.x && p.y>=d.minCoord.y && p.y<=d.maxCoord.y;};
    const auto run=[&](const Domain& domain,Vector2D point,Vector2D speed,double dt,int side){
        std::vector<std::array<double,3>> primalSamples,jetSamples;
        const auto record=[&](auto& samples,Vector2D p,int layer){check(inside(domain,p),"near-wall callback sampled exterior midpoint");samples.push_back({double(layer),p.x,p.y});};
        const auto primal=LevelSetBoundaryTrace::trace(domain,point,dt,
            [&](Vector2D p){record(primalSamples,p,1);return speed;},[&](Vector2D p){record(primalSamples,p,0);return speed;});
        const auto tangent=LevelSetBoundaryTrace::traceWithJacobian(domain,point,dt,
            [&](Vector2D p){record(jetSamples,p,1);return Jet{speed,{}};},[&](Vector2D p){record(jetSamples,p,0);return Jet{speed,{}};});
        const auto& result=tangent.trajectory;const bool horizontal=side==1 || side==3;
        const double wall=horizontal?(side==1?domain.minCoord.y:domain.maxCoord.y):(side==0?domain.minCoord.x:domain.maxCoord.x);
        const double distance=std::abs((horizontal?point.y:point.x)-wall),normal=std::abs(horizontal?speed.y:speed.x),exactTime=distance/normal;
        const double scale=std::max({domain.cellSize.x,domain.cellSize.y,std::abs(domain.minCoord.x),std::abs(domain.minCoord.y),std::abs(domain.maxCoord.x),std::abs(domain.maxCoord.y)});
        const double timeTolerance=128*eps*dt,coordinateTolerance=128*eps*scale;
        check(result.boundaryHit && result.roundoffMidpointEvent && result.side==sides[side]
            && !tangent.derivative && tangent.acceptedSubsteps==0,"near-wall event missing or fabricated interior tangent");
        check(primal.boundaryHit && primal.roundoffMidpointEvent && primal.side==result.side
            && primal.point.x==result.point.x && primal.point.y==result.point.y && primal.elapsed==result.elapsed
            && primalSamples==jetSamples,"primal/jet near-wall map or sampling order differs");
        check((horizontal?result.point.y:result.point.x)==wall && result.elapsed>=0 && result.elapsed<=dt
            && std::abs(result.elapsed-exactTime)<=timeTolerance,"near-wall physical hit time/side reference");
        const auto reference=point-speed*exactTime;
        check(std::abs(result.point.x-reference.x)<=coordinateTolerance && std::abs(result.point.y-reference.y)<=coordinateTolerance,"near-wall position outside original event resolution");
        check(primalSamples.size()<=3,"unsafe midpoint shortening sampled velocity");maximumQueries=std::max(maximumQueries,primalSamples.size());
        if(result.elapsed==0)++zeroTimes;++cases;
    };
    for(bool radial:{false,true})for(int side=0;side<4;++side){
        if(radial && side!=1 && side!=3)continue;
        const bool horizontal=side==1 || side==3,upper=side>=2;
        Domain domain{{0,0},{1,1},{.25,.25},radial};
        if(upper){if(horizontal){domain.minCoord.y=-1;domain.maxCoord.y=0;}else{domain.minCoord.x=-1;domain.maxCoord.x=0;}}
        for(double normal:{.5,8.}){
            const double dt=normal==.5?.01:.001;Vector2D speed=horizontal?Vector2D{radial?0:.125,upper?-normal:normal}:Vector2D{upper?-normal:normal,.125};
            for(double distance:{eta,1e-300,1e-100,1e-16}){
                const double p=upper?-distance:distance;run(domain,horizontal?Vector2D{.5,p}:Vector2D{p,.5},speed,dt,side);}
            auto shifted=domain;
            if(horizontal){shifted.minCoord.y+=.3;shifted.maxCoord.y+=.3;}else{shifted.minCoord.x+=.3;shifted.maxCoord.x+=.3;}
            const double wall=.3,p=std::nextafter(wall,upper?-std::numeric_limits<double>::infinity():std::numeric_limits<double>::infinity());
            run(shifted,horizontal?Vector2D{.5,p}:Vector2D{p,.5},speed,dt,side);
        }
        const double p=upper?-eta:eta;const Vector2D speed=horizontal?Vector2D{0,upper?-8.:8.}:Vector2D{upper?-8.:8.,0};
        run(domain,horizontal?Vector2D{.5,p}:Vector2D{p,.5},speed,eta,side);
    }
    const Domain domain{{0,0},{1,1},{.25,.25},false};
    const auto rawReject=[&](auto action,const char* fragment){
        try{action();}catch(const std::exception& e){check(std::string(e.what()).find(fragment)!=std::string::npos,"new event guard returned wrong reason");++guards;return;}
        throw std::runtime_error(std::string("missing raw event rejection: ")+fragment);
    };
    const auto constant=[](Vector2D speed){return [speed](Vector2D){return speed;};};
    rawReject([&]{LevelSetBoundaryTrace::trace(domain,{eta,eta},.01,constant({1,1}),constant({1,1}));},"multiple wall crossings");
    const auto tangent=[](Vector2D p){return Vector2D{0,std::sqrt(p.y)};};
    rawReject([&]{LevelSetBoundaryTrace::trace(domain,{.5,1e-100},.01,tangent,tangent);},"tangential, non-inflow");
    const auto reversed=[](Vector2D p){return Vector2D{0,p.y==0?-.5:.5};};
    rawReject([&]{LevelSetBoundaryTrace::trace(domain,{.5,eta},.01,reversed,reversed);},"tangential, non-inflow");
    auto covered=domain;covered.activeSegment=[](Vector2D,Vector2D q){return q.y>0;};
    rawReject([&]{LevelSetBoundaryTrace::trace(covered,{.5,eta},.01,constant({0,.5}),constant({0,.5}));},"selected fluid domain");
    auto radial=domain;radial.radialAxisSymmetry=true;
    rawReject([&]{LevelSetBoundaryTrace::trace(radial,{eta,.5},.01,constant({.5,0}),constant({.5,0}));},"negative r");
    rawReject([&]{LevelSetBoundaryTrace::trace(domain,{.5,eta},.01,constant({0,100}),constant({0,100}));},"half a cell");
    const auto invalid=[](Vector2D p){return Vector2D{0,p.y==0?std::numeric_limits<double>::quiet_NaN():.5};};
    rawReject([&]{LevelSetBoundaryTrace::trace(domain,{.5,eta},.01,invalid,invalid);},"nonfinite sampled velocity");
    const auto invalidJet=[](Vector2D p){return Jet{{0,.5},{0,0,0,p.y==0?std::numeric_limits<double>::quiet_NaN():0}};};
    rawReject([&]{LevelSetBoundaryTrace::traceWithJacobian(domain,{.5,eta},.01,invalidJet,invalidJet);},"velocity Jacobian");
    // A resolved first substep can also end very near a wall. Continuous
    // piecewise-linear callback; queried plateaus have their exact zero slope.
    const double nearHalf=std::nextafter(std::nextafter(.5,0.),0.);
    const auto field=[&](Vector2D p){check(inside(domain,p),"post-substep exterior query");
        const double blend=std::clamp((p.y-.001)/.001,0.,1.),slope=p.y>.001 && p.y<.002?(1-nearHalf)/.001:0.;
        return Jet{{0,nearHalf+(1-nearHalf)*blend},{0,0,0,slope}};};
    const auto after=LevelSetBoundaryTrace::traceWithJacobian(domain,{.5,.03125},.125,field,field);
    check(after.trajectory.roundoffMidpointEvent && after.trajectory.side==Boundary::south && !after.derivative
        && after.acceptedSubsteps==1 && after.trajectory.elapsed>=.0625 && after.trajectory.elapsed<=.125,"near-wall event after an accepted substep failed");
    const auto ordinary=LevelSetBoundaryTrace::trace(domain,{.5,.005},.01,constant({0,.5}),constant({0,.5}));
    const auto initial=LevelSetBoundaryTrace::trace(domain,{.5,0},.01,constant({0,.5}),constant({0,.5}));
    const auto identity=LevelSetBoundaryTrace::trace(domain,{.5,eta},0,constant({0,.5}),constant({0,.5}));
    check(ordinary.boundaryHit && !ordinary.roundoffMidpointEvent && initial.boundaryHit && !initial.roundoffMidpointEvent
        && !identity.boundaryHit && !identity.roundoffMidpointEvent,"new branch changed ordinary/initial/identity paths");
    check(cases==66 && zeroTimes>=6 && guards==8,"near-wall event matrix incomplete");
    std::cout<<"I12a16r4 roundoff midpoint PASS cases="<<cases<<" post_substep=1 zero_times="<<zeroTimes
        <<" max_queries="<<maximumQueries<<" guards="<<guards<<'\n';
}
}

int main()
{
    try {
        std::cout<<std::setprecision(17);
        const Context context;
        affineAndTime(context);differenceSweep(context);substepWitness(context);
        rotationConvergence(context);boundaryAxisAndFailures(context);
        std::cout<<"I02 characteristic B/J contract PASS (tracer only; Gate I remains incomplete)\n";
        roundoffMidpointEvents();
    } catch(const std::exception& error) {std::cerr<<error.what()<<'\n';return 1;}
    return 0;
}
