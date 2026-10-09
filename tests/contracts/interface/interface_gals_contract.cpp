#include "Core/Geometry/ActiveQuarticHermiteGals.h"
#include "Core/Geometry/LevelSetClosedRk4Trace.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <vector>

namespace
{
namespace Trace=LevelSetBoundaryTrace;
using Patch=Quantum::Geometry::ActiveQuarticHermitePatch;
using Jet=LevelSetHermiteGeometry::Snapshot::ScalarJet2;
void require(bool valid,const char* message) {if(!valid) throw std::runtime_error(message);}
void close(double actual,double expected,const char* message,double tolerance=2e-11)
{ require(std::isfinite(actual) && std::abs(actual-expected)<=tolerance*(1+std::abs(expected)),message); }
template<class Profile> Patch prescribed(const QuadtreeMesh& mesh,const Profile& profile)
{
    Patch patch;patch.activeCells.assign(mesh.cellIndexNum(),1);patch.nodes.resize(mesh.vertIndexNum());
    patch.cellCenters.resize(mesh.cellIndexNum());std::map<std::array<int,2>,bool> edges;
    for(int v=0;v<mesh.vertIndexNum();++v) {
        const auto j=profile(mesh.verts(mesh.vertIndexToLabel[v]).coord());patch.nodes[v]={j.value,j.dx,j.dy,j.dxy};
    }
    for(int c=0;c<mesh.cellIndexNum();++c) {
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);patch.cellCenters[c]=profile(cell.centerCoord()).value;
        for(int k=0;k<4;++k) {
            const auto& a=mesh.verts(cell.cellToVert[k]);const auto& b=mesh.verts(cell.cellToVert[(k+1)%4]);
            const bool horizontal=a.coord().y==b.coord().y;std::array<int,2> key{{a.index(),b.index()}};
            if(horizontal?a.coord().x>b.coord().x:a.coord().y>b.coord().y) std::swap(key[0],key[1]);
            if(!edges.emplace(key,true).second) continue;
            const auto j=profile(a.coord()+(b.coord()-a.coord())*.5);
            patch.edgeMidpoints.push_back({key,j.value,horizontal?j.dy:j.dx});
        }
    }
    return patch;
}
void samePatch(const Patch& a,const Patch& b)
{
    require(a.nodes.size()==b.nodes.size() && a.edgeMidpoints.size()==b.edgeMidpoints.size(),"GALS patch shape changed");
    for(std::size_t k=0;k<a.nodes.size();++k) {
        const auto& p=a.nodes[k];const auto& q=b.nodes[k];
        require(p.value==q.value && p.dx==q.dx && p.dy==q.dy && p.dxy==q.dxy,"GALS stationary/parallel vertex mismatch");
    }
    for(std::size_t k=0;k<a.edgeMidpoints.size();++k) {
        const auto& p=a.edgeMidpoints[k];const auto& q=b.edgeMidpoints[k];
        require(p.vertices==q.vertices && p.value==q.value && p.normalDerivative==q.normalDerivative,
            "GALS stationary/parallel edge mismatch");
    }
    require(a.cellCenters==b.cellCenters,"GALS stationary/parallel center mismatch");
}
void traceContract()
{
    const Trace::Domain domain{{0,0},{1,1},{1,1},false};
    // The midpoint crosses the wall for a full step, but actual half steps
    // stay inside as the velocity reverses. Test tangent composition, not
    // merely a one-stage formula. The primal call sequence must be identical.
    constexpr double a=.04,dt=.1;
    std::vector<std::array<double,3>> primalQueries,jetQueries;
    const auto value=[=](Vector2D p,bool current){return Vector2D{a*p.x+(current?.3:-.3),.02*p.y};};
    const Trace::Velocity current=[&](const Vector2D& p){primalQueries.push_back({p.x,p.y,1});return value(p,true);};
    const Trace::Velocity old=[&](const Vector2D& p){primalQueries.push_back({p.x,p.y,0});return value(p,false);};
    const Trace::DifferentiableVelocity currentJet=[&](const Vector2D& p){jetQueries.push_back({p.x,p.y,1});return Trace::VelocityJet{value(p,true),{a,0,0,.02}};};
    const Trace::DifferentiableVelocity oldJet=[&](const Vector2D& p){jetQueries.push_back({p.x,p.y,0});return Trace::VelocityJet{value(p,false),{a,0,0,.02}};};
    const Vector2D point{.01,.45};
    const auto primal=Trace::trace(domain,point,dt,current,old);
    const auto jet=Trace::traceWithJacobian(domain,point,dt,currentJet,oldJet);
    require(primalQueries==jetQueries,"differentiated trace changed primal sample order");
    require(primal.point.x==jet.trajectory.point.x && primal.point.y==jet.trajectory.point.y &&
        primal.elapsed==jet.trajectory.elapsed && primal.side==jet.trajectory.side,"differentiated trace changed primal map");
    require(jet.derivative.has_value() && jet.acceptedSubsteps==2,"two-substep GALS control did not exercise composition");
    const double bx=1-.05*a+.5*.05*.05*a*a,by=1-.05*.02+.5*.05*.05*.02*.02;
    close(jet.derivative->xx,bx*bx,"GALS substep tangent xx",2e-14);
    close(jet.derivative->yy,by*by,"GALS substep tangent yy",2e-14);
    // Finite differences are an independent CHECK, not production transport.
    constexpr double epsilon=1e-5;
    const auto plus=Trace::trace(domain,{point.x+epsilon,point.y},dt,current,old);
    const auto minus=Trace::trace(domain,{point.x-epsilon,point.y},dt,current,old);
    close((plus.point.x-minus.point.x)/(2*epsilon),jet.derivative->xx,"GALS derivative differs from actual map",2e-10);
    const Trace::DifferentiableVelocity inflow=[](const Vector2D&){return Trace::VelocityJet{{.2,0},{}};};
    const auto boundary=Trace::traceWithJacobian(domain,{.005,.5},dt,inflow,inflow);
    require(boundary.trajectory.boundaryHit && boundary.trajectory.side==Boundary::west && !boundary.derivative,
        "GALS boundary hit invented an interior Jacobian");
    const auto zeroTime=Trace::traceWithJacobian(domain,point,0,currentJet,oldJet);
    require(zeroTime.derivative && zeroTime.derivative->xx==1 && zeroTime.derivative->yy==1 && zeroTime.acceptedSubsteps==0,
        "zero-time GALS map is not identity");
    bool rejected=false;
    try {
        const Trace::DifferentiableVelocity bad=[](const Vector2D&){return Trace::VelocityJet{{0,0},
            {std::numeric_limits<double>::quiet_NaN(),0,0,0}};};
        Trace::traceWithJacobian(domain,point,dt,bad,bad);
    } catch(const std::exception&) {rejected=true;}
    require(rejected,"GALS accepted a nonfinite velocity derivative");
}
Jet polynomial(Vector2D p,bool radial)
{
    const double x=p.x,y=p.y;
    const double f=radial?1+x*x+.2*x*x*x*x:1+.3*x+.7*x*x-.4*x*x*x+.2*x*x*x*x;
    const double fx=radial?2*x+.8*x*x*x:.3+1.4*x-1.2*x*x+.8*x*x*x;
    const double fxx=radial?2+2.4*x*x:1.4-2.4*x+2.4*x*x;
    const double g=.2+.3*y-.6*y*y+.2*y*y*y+.1*y*y*y*y;
    const double gy=.3-1.2*y+.6*y*y+.4*y*y*y,gyy=-1.2+1.2*y+1.2*y*y;
    return {f*g-.15,fx*g,f*gy,fxx*g,fx*gy,f*gyy};
}
void polynomialPullback(bool radial,int nx,int ny,bool fourth=false)
{
    DictionaryMesh input;input.minCoord={0,0};input.maxCoord={1,1};input.Nx=nx;input.Ny=ny;
    input.maxCellNum=512;input.maxInitialIterations=0;QuadtreeMesh mesh(input);
    constexpr double dt=.0125;
    const auto velocity=[=](bool current) -> Trace::DifferentiableVelocity {
        const double ax=current?.13:.07,ay=current?.19:.12;
        return [=](const Vector2D& p){return Trace::VelocityJet{{ax*(p.x-(radial?0:.5)),ay*(p.y-.5)},{ax,0,0,ay}};};
    };
    const double bx=fourth?std::exp(-dt*.10):1-dt*.10*(1-.5*dt*.13);
    const double by=fourth?std::exp(-dt*.155):1-dt*.155*(1-.5*dt*.19);
    const auto expected=[=](Vector2D p) {
        const double ox=radial?0:.5;
        auto j=polynomial({ox+bx*(p.x-ox),.5+by*(p.y-.5)},radial);
        j.dx*=bx;j.dy*=by;j.dxx*=bx*bx;j.dxy*=bx*by;j.dyy*=by*by;return j;
    };
    const auto previous=prescribed(mesh,[=](Vector2D p){return polynomial(p,radial);});
    ActiveQuarticHermiteGals::Options options;options.radialAxisSymmetry=radial;
    options.integrator=fourth?ActiveQuarticHermiteGals::Integrator::ClosedRk4:ActiveQuarticHermiteGals::Integrator::Midpoint;
    const auto candidate=ActiveQuarticHermiteGals::predictUniform(mesh,previous,dt,velocity(true),velocity(false),options);
    for(double theta:{0.,.3,1.}) {
        const Vector2D p{.37,.42};const auto traced=candidate.backwardAt(theta,p);
        const double ax=.07+theta*.06,ay=.12+theta*.07,s=theta*dt;
        const double xx=fourth?std::exp(-s*(.07+.5*theta*.06)):1-s*(.07+.5*theta*.06)*(1-.5*s*ax);
        const double yy=fourth?std::exp(-s*(.12+.5*theta*.07)):1-s*(.12+.5*theta*.07)*(1-.5*s*ay);
        close(traced.trajectory.point.x,(radial?0:.5)+xx*(p.x-(radial?0:.5)),"GALS fractional trace x");
        close(traced.trajectory.point.y,.5+yy*(p.y-.5),"GALS fractional trace y");
        close(traced.derivative->xx,xx,"GALS fractional trace Jxx");close(traced.derivative->yy,yy,"GALS fractional trace Jyy");
    }
    require(candidate.report.tracedPoints==static_cast<std::size_t>((2*nx+1)*(2*ny+1)),"GALS traced extra stencil points");
    for(int v=0;v<mesh.vertIndexNum();++v) {
        const auto e=expected(mesh.verts(mesh.vertIndexToLabel[v]).coord());const auto& j=candidate.patch.nodes[v];
        close(j.value,e.value,"Q4 GALS vertex value");close(j.dx,e.dx,"Q4 GALS vertex dx");
        close(j.dy,e.dy,"Q4 GALS vertex dy");close(j.dxy,e.dxy,"Q4 GALS compact mixed derivative");
    }
    for(const auto& edge:candidate.patch.edgeMidpoints) {
        const auto a=mesh.verts(mesh.vertIndexToLabel[edge.vertices[0]]).coord();
        const auto b=mesh.verts(mesh.vertIndexToLabel[edge.vertices[1]]).coord();const auto e=expected(a+(b-a)*.5);
        close(edge.value,e.value,"Q4 GALS edge value");
        close(edge.normalDerivative,a.y==b.y?e.dy:e.dx,"Q4 GALS edge normal derivative");
    }
    for(int c=0;c<mesh.cellIndexNum();++c) {
        const auto p=mesh.cells(mesh.cellIndexToLabel[c]).centerCoord();const auto e=expected(p);
        close(candidate.patch.cellCenters[c],e.value,"Q4 GALS center prediction");
        const auto j=candidate.geometry->scalarJetInNativeCell(c,p);
        close(j.dxx,e.dxx,"Q4 GALS Hessian xx",1e-10);close(j.dxy,e.dxy,"Q4 GALS Hessian xy",1e-10);
        close(j.dyy,e.dyy,"Q4 GALS Hessian yy",1e-10);
    }
    const Trace::DifferentiableVelocity zero=[](const Vector2D&){return Trace::VelocityJet{{0,0},{}};};
    const auto stationary=ActiveQuarticHermiteGals::predictUniform(mesh,candidate.patch,dt,zero,zero,options);
    samePatch(stationary.patch,candidate.patch);
    options.independentVelocityQueries=true;
    const auto parallel=ActiveQuarticHermiteGals::predictUniform(mesh,previous,dt,velocity(true),velocity(false),options);
    samePatch(candidate.patch,parallel.patch);
    require(candidate.report.maxScaledMixedDisagreement<1e-12,"Q4 GALS mixed closures lost polynomial compatibility");
    std::cout<<"GALS Q4 "<<(radial?"RZ":"XY")<<' '<<nx<<'x'<<ny<<" integrator="<<(fourth?4:2)<<" traces="<<candidate.report.tracedPoints
        <<" mixed_disagreement="<<candidate.report.maxScaledMixedDisagreement<<" stationary_exact=1 parallel_exact=1\n";
}
void nonlinearPullback()
{
    DictionaryMesh input;input.minCoord={0,0};input.maxCoord={1,1};input.Nx=input.Ny=4;
    input.maxCellNum=64;input.maxInitialIterations=0;QuadtreeMesh mesh(input);
    const auto profile=[](Vector2D p)->Jet{return {p.x+2*p.y-.9,1,2,0,0,0};};
    constexpr double dt=.01,a=.1,c=.06,d=.2,q=1-.5*dt*d;
    const Trace::DifferentiableVelocity velocity=[=](const Vector2D& p) {
        return Trace::VelocityJet{{(a+c*(p.y-.5))*(p.x-.5),d*(p.y-.5)},
            {a+c*(p.y-.5),c*(p.x-.5),0,d}};
    };
    const auto source=prescribed(mesh,profile);
    ActiveQuarticHermiteGals::Options options;options.integrator=ActiveQuarticHermiteGals::Integrator::Midpoint;
    const auto moved=ActiveQuarticHermiteGals::predictUniform(mesh,source,dt,velocity,velocity,options);
    for(int v=0;v<mesh.vertIndexNum();++v) {
        const auto p=mesh.verts(mesh.vertIndexToLabel[v]).coord();const double x=p.x-.5,y=p.y-.5;
        const double f=1-dt*(a+c*q*y)*(1-.5*dt*(a+c*y));
        const double df=-dt*(c*q*(1-.5*dt*(a+c*y))-.5*dt*c*(a+c*q*y));
        const double by=1-d*dt+.5*d*d*dt*dt;const auto& j=moved.patch.nodes[v];
        close(j.value,x*f+.5+2*(by*y+.5)-.9,"nonlinear GALS value");
        close(j.dx,f,"nonlinear GALS dx");close(j.dy,x*df+2*by,"nonlinear GALS dy");
        close(j.dxy,df,"nonlinear GALS compact mixed derivative");
        if(p.x==.5 && p.y==.5) require(j.value==source.nodes[v].value && j.dx!=source.nodes[v].dx,
            "GALS froze the gradient at a stationary point");
    }
    const Trace::DifferentiableVelocity inflow=[](const Vector2D&){return Trace::VelocityJet{{.02,0},{}};};
    bool rejected=false;
    try {ActiveQuarticHermiteGals::predictUniform(mesh,source,dt,inflow,inflow);}
    catch(const std::exception&) {rejected=true;}
    require(rejected,"GALS predictor silently supplied boundary-history derivatives");
    samePatch(source,prescribed(mesh,profile));
}
void closedRk4Contract()
{
    const Trace::Domain domain{{0,0},{1,1},{1,1},false};const Vector2D p{.65,.43};
    const Trace::DifferentiableVelocity rotation=[](const Vector2D& p){return Trace::VelocityJet{{-(p.y-.5),p.x-.5},{0,-1,1,0}};};
    double coarse=0,fine=0;
    for(double dt:{.1,.05}) {
        const auto q=Trace::traceClosedRk4WithJacobian(domain,p,dt,rotation,rotation);
        const Vector2D exact{.5+(p.x-.5)*std::cos(dt)+(p.y-.5)*std::sin(dt),.5-(p.x-.5)*std::sin(dt)+(p.y-.5)*std::cos(dt)};
        const double error=std::hypot(q.trajectory.point.x-exact.x,q.trajectory.point.y-exact.y);
        (dt==.1?coarse:fine)=error;
        close(q.derivative->xx,std::cos(dt),"RK4 rotation tangent xx",2e-9);
        close(q.derivative->xy,std::sin(dt),"RK4 rotation tangent xy",1e-7);
    }
    require(coarse/fine>30 && coarse/fine<34,"RK4 map did not have fifth-order local rotation error");
    const Trace::DifferentiableVelocity nonlinear=[](const Vector2D& p){return Trace::VelocityJet{{.1*p.x*p.y,.05*p.x*p.x},{.1*p.y,.1*p.x,.1*p.x,0}};};
    const auto q=Trace::traceClosedRk4WithJacobian(domain,p,.04,nonlinear,nonlinear);const double epsilon=1e-5;
    for(int axis=0;axis<2;++axis) {
        auto plus=p,minus=p;(axis?plus.y:plus.x)+=epsilon;(axis?minus.y:minus.x)-=epsilon;
        const auto a=Trace::traceClosedRk4WithJacobian(domain,plus,.04,nonlinear,nonlinear).trajectory.point;
        const auto b=Trace::traceClosedRk4WithJacobian(domain,minus,.04,nonlinear,nonlinear).trajectory.point;
        close((a.x-b.x)/(2*epsilon),axis?q.derivative->xy:q.derivative->xx,"RK4 nonlinear discrete tangent x",2e-10);
        close((a.y-b.y)/(2*epsilon),axis?q.derivative->yy:q.derivative->yx,"RK4 nonlinear discrete tangent y",2e-10);
    }
    std::cout<<"Closed RK4 local_rotation_error_ratio="<<coarse/fine<<" nonlinear_discrete_tangent=1\n";
}
}

void interfaceGalsContract()
{
    traceContract();polynomialPullback(false,2,2);polynomialPullback(false,6,3);
    polynomialPullback(true,4,4);nonlinearPullback();closedRk4Contract();
    polynomialPullback(false,4,4,true);polynomialPullback(true,4,4,true);
    std::cout<<"Coherent trace and compact Q4 first-jet predictor passed; predictor_controls_only=1\n";
}
