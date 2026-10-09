#include "Core/Geometry/ActiveQuarticHermiteGals.h"
#include "Core/Geometry/CompactQuarticDerivative.h"
#include "Core/Geometry/LevelSetClosedRk4Trace.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/UniformLeafGrid.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace ActiveQuarticHermiteGals
{
namespace
{
using Patch=Quantum::Geometry::ActiveQuarticHermitePatch;
using Trace=LevelSetBoundaryTrace::JetResult;
void require(bool valid,const char* reason)
{ if(!valid) throw std::runtime_error(std::string("[ActiveQuarticHermiteGals] ")+reason); }
bool same(Vector2D a,Vector2D b) { return a.x==b.x && a.y==b.y; }
bool identity(const LevelSetBoundaryTrace::Jacobian& a)
{ return a.xx==1 && a.xy==0 && a.yx==0 && a.yy==1; }
struct Point
{
    Vector2D position;
    double value=0.0;
    Vector2D gradient;
};
struct Advanced
{
    double value=0.0;
    Vector2D gradient;
    double traceRatio=0.0;
    int substeps=0;
    bool identity=false;
};

// Retain this predictor's error context while sharing the same pure stencil.
template<class Coordinate,class Value>
double derivativeAtNode(int target,int count,const Coordinate& coordinate,const Value& value)
{
    return Quantum::Geometry::compactQuarticDerivativeAtNode(target,count,coordinate,value,require);
}
}

Candidate predictUniform(const QuadtreeMesh& mesh,const Patch& previous,double dt,
    const LevelSetBoundaryTrace::DifferentiableVelocity& current,
    const LevelSetBoundaryTrace::DifferentiableVelocity& old,const Options& options)
{
    require(std::isfinite(dt) && dt>=0 && current && old,"finite nonnegative dt and both velocity jets required");
    require(previous.cornerSectors.empty(),"contact sectors require the boundary GALS entry");
    const UniformLeafGrid grid(mesh,CartesianGridScope::UniformLeaves);
    require(grid.nx()>=2 && grid.ny()>=2,"compact Q4 closure needs at least two cells per direction");
    require(previous.activeCells.size()==static_cast<std::size_t>(mesh.cellIndexNum()) &&
        std::all_of(previous.activeCells.begin(),previous.activeCells.end(),[](auto active){return active==1;}),
        "uniform GALS currently requires a complete active domain");
    // The common scalar owner validates all independent data and native edge
    // identities before any callback. No second geometric representation.
    const auto source=std::make_shared<const LevelSetHermiteGeometry::Snapshot>(
        LevelSetHermiteGeometry::Snapshot::fromNativeLeaves(mesh,previous,options.radialAxisSymmetry));
    const LevelSetBoundaryTrace::Domain domain{grid.minCoord(),grid.maxCoord(),grid.spacing(),options.radialAxisSymmetry};
    Candidate result;result.patch=previous;
    result.maximumDisplacement=dt==0?0:.5*std::min(grid.spacing().x,grid.spacing().y)*
        (1+64*std::numeric_limits<double>::epsilon());
    result.independentVelocityQueries=options.independentVelocityQueries;
    require(options.integrator==Integrator::Midpoint || options.integrator==Integrator::ClosedRk4,"unknown characteristic integrator");
    const auto trace=[method=options.integrator](const auto& domain,Vector2D p,double dt,const auto& current,const auto& old) {
        return method==Integrator::Midpoint?LevelSetBoundaryTrace::traceWithJacobian(domain,p,dt,current,old):
            LevelSetBoundaryTrace::traceClosedRk4WithJacobian(domain,p,dt,current,old);
    };
    result.backward=[domain,dt,current,old,trace](Vector2D point) {
        return trace(domain,point,dt,current,old);
    };
    result.backwardAt=[domain,dt,current,old,trace,full=result.backward](double theta,Vector2D point) {
        require(std::isfinite(theta) && theta>=0 && theta<=1,"invalid characteristic time fraction");
        if(theta==1)return full(point);
        const auto stage=[current,old,theta](const Vector2D& p) {
            const auto a=old(p),b=current(p);LevelSetBoundaryTrace::VelocityJet value;
            value.value=a.value*(1-theta)+b.value*theta;
            value.derivative={a.derivative.xx*(1-theta)+b.derivative.xx*theta,
                a.derivative.xy*(1-theta)+b.derivative.xy*theta,
                a.derivative.yx*(1-theta)+b.derivative.yx*theta,
                a.derivative.yy*(1-theta)+b.derivative.yy*theta};
            return value;
        };
        return trace(domain,point,dt*theta,stage,old);
    };
    const int nv=mesh.vertIndexNum(), ne=static_cast<int>(previous.edgeMidpoints.size()), nc=mesh.cellIndexNum();
    std::vector<Point> points(static_cast<std::size_t>(nv)+ne+nc);
    for(int v=0;v<nv;++v) {
        const auto& node=previous.nodes[v];
        points[v]={mesh.verts(mesh.vertIndexToLabel[v]).coord(),node.value,{node.dx,node.dy}};
    }
    std::map<std::array<int,2>,int> edges;
    std::vector<bool> horizontal(ne);
    for(int k=0;k<ne;++k) {
        const auto& edge=previous.edgeMidpoints[k];
        require(edges.emplace(edge.vertices,k).second,"duplicate skeleton edge");
        const auto a=mesh.verts(mesh.vertIndexToLabel[edge.vertices[0]]).coord();
        const auto b=mesh.verts(mesh.vertIndexToLabel[edge.vertices[1]]).coord();
        horizontal[k]=a.y==b.y;
        points[nv+k].position=a+(b-a)*0.5;
        points[nv+k].value=edge.value;
        // Only dn is a stored midpoint datum; the other component is not an
        // independent constraint. An identity trace retains this datum below.
        if(horizontal[k]) points[nv+k].gradient.y=edge.normalDerivative;
        else points[nv+k].gradient.x=edge.normalDerivative;
    }
    for(int c=0;c<nc;++c) points[nv+ne+c]={mesh.cells(mesh.cellIndexToLabel[c]).centerCoord(),previous.cellCenters[c],{}};
    std::vector<Advanced> advanced(points.size());
    const auto evaluate=[&](std::size_t k) {
        const auto& point=points[k];const Trace trace=result.backward(point.position);
        require(!trace.trajectory.boundaryHit && trace.derivative.has_value(),
            "first-jet predictor encountered physical inflow; boundary-history derivative is not enabled");
        Advanced a;a.identity=same(trace.trajectory.point,point.position) && identity(*trace.derivative);
        a.substeps=trace.acceptedSubsteps;
        a.traceRatio=std::hypot(trace.trajectory.point.x-point.position.x,trace.trajectory.point.y-point.position.y)
            /std::min(grid.spacing().x,grid.spacing().y);
        if(a.identity) {a.value=point.value;a.gradient=point.gradient;}
        else {
            const auto jet=source->scalarJetNative(trace.trajectory.point);
            // A stationary point can still have a nonidentity map Jacobian.
            // Its stored scalar value remains exact while its gradient evolves.
            a.value=same(trace.trajectory.point,point.position)?point.value:jet.value;
            a.gradient=trace.derivative->transposeTimes({jet.dx,jet.dy});
        }
        require(std::isfinite(a.value) && std::isfinite(a.gradient.x) && std::isfinite(a.gradient.y),
            "nonfinite transported first jet");
        advanced[k]=a;
    };
    if(options.independentVelocityQueries) {
        // Queries write only their own slots; failures publish no candidate.
        std::vector<std::exception_ptr> failures(points.size());
#pragma omp parallel for schedule(static)
        for(int k=0;k<static_cast<int>(points.size());++k) {
            try {evaluate(static_cast<std::size_t>(k));}
            catch(...) {failures[k]=std::current_exception();}
        }
        for(const auto& failure:failures) if(failure) std::rethrow_exception(failure);
    } else for(std::size_t k=0;k<points.size();++k) evaluate(k);
    bool skeletonChanged=false;
    for(std::size_t k=0;k<points.size();++k) {
        const auto& a=advanced[k];result.report.acceptedSubsteps+=a.substeps;
        result.report.maxTraceRatio=std::max(result.report.maxTraceRatio,a.traceRatio);
        bool changed=false;
        if(k<static_cast<std::size_t>(nv)) {
            auto& node=result.patch.nodes[k];
            changed=node.value!=a.value || node.dx!=a.gradient.x || node.dy!=a.gradient.y;
            node.value=a.value;node.dx=a.gradient.x;node.dy=a.gradient.y;
        } else if(k<static_cast<std::size_t>(nv+ne)) {
            const int e=static_cast<int>(k)-nv;auto& edge=result.patch.edgeMidpoints[e];
            const double dn=horizontal[e]?a.gradient.y:a.gradient.x;
            changed=edge.value!=a.value || edge.normalDerivative!=dn;
            edge.value=a.value;edge.normalDerivative=dn;
        } else result.patch.cellCenters[k-nv-ne]=a.value;
        skeletonChanged=skeletonChanged || changed;
        result.report.changedFirstJetPoints+=changed;
    }
    result.report.tracedPoints=points.size();
    if(skeletonChanged) {
        // Uniform layout contains only independent regular vertices/edges.
        // Each row samples phi_y at vertices and horizontal edge midpoints;
        // each column samples phi_x at vertices and vertical edge midpoints.
        const int nx=grid.nx(),ny=grid.ny();
        const auto vertex=[&](int i,int j){return mesh.verts(grid.vertices()[j*(nx+1)+i]).index();};
        const auto edgePoint=[&](int a,int b) {
            const auto found=edges.find({{a,b}});
            require(found!=edges.end(),"missing maximal skeleton edge");
            return nv+found->second;
        };
        for(int j=0;j<=ny;++j) for(int i=0;i<=nx;++i) {
            const auto rowPoint=[&](int s) {
                return s%2?edgePoint(vertex(s/2,j),vertex(s/2+1,j)):vertex(s/2,j);
            };
            const auto columnPoint=[&](int s) {
                return s%2?edgePoint(vertex(i,s/2),vertex(i,s/2+1)):vertex(i,s/2);
            };
            const auto x=derivativeAtNode(2*i,2*nx+1,[&](int s){return points[rowPoint(s)].position.x;},
                [&](int s){return advanced[rowPoint(s)].gradient.y;});
            const auto y=derivativeAtNode(2*j,2*ny+1,[&](int s){return points[columnPoint(s)].position.y;},
                [&](int s){return advanced[columnPoint(s)].gradient.x;});
            result.report.maxScaledMixedDisagreement=std::max(result.report.maxScaledMixedDisagreement,
                std::abs(x-y)*grid.spacing().x*grid.spacing().y);
            auto& node=result.patch.nodes[vertex(i,j)];node.dxy=0.5*x+0.5*y;
            // The source and velocity enforce even scalar/odd radial velocity.
            // At r=0 the mixed derivative is analytically zero. The one-sided
            // off-axis closure is not an independent boundary datum.
            if(options.radialAxisSymmetry && points[vertex(i,j)].position.x==0) {
                require(node.dx==0,"GALS pullback violated the radial-axis first jet");
                node.dxy=0;
            }
        }
    }
    result.geometry=std::make_unique<LevelSetHermiteGeometry::Snapshot>(
        LevelSetHermiteGeometry::Snapshot::fromNativeLeaves(mesh,result.patch,options.radialAxisSymmetry));
    return result;
}
}
