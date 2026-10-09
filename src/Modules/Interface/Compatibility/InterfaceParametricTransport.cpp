#include "Modules/Interface/Compatibility/InterfaceParametricTransport.h"
#include <numeric>
#include <sstream>
#include <chrono>
#include <atomic>
#include <cstdlib>

namespace InterfaceParametricTransport {
using namespace Quantum::Geometry::Parametric;
using V=Point;
using TransportCurve=Source;
namespace {
void require(bool ok,const char* reason){if(!ok)throw std::runtime_error(reason);}
}
V applyJet(const MapJet& jet,V tangent) {return {jet.dx.x*tangent.x+jet.dy.x*tangent.y,jet.dx.y*tangent.x+jet.dy.y*tangent.y};}
TransportCurve transportSource(const Curve& source,const InterfaceParametricProjection::Grid& grid,double dt,const Flow& flow,SourceReport& diagnostics) {
    require(grid.nx>0 && grid.ny>0 && grid.width>0 && grid.height>0 && bool(flow) && std::isfinite(dt) && !source.arcs.empty(),"invalid parametric first-jet source input");
    TransportCurve out;out.boundary=source.boundary;std::array<std::vector<double>,2> levels;
    for(int axis=0;axis<2;axis++){
        if(grid.leaves){for(const auto& p:grid.leaves->planes[axis])levels[axis].push_back(p.coordinate);
            for(const auto& c:grid.leaves->cells)levels[axis].push_back(axis?.5*(c.bottom+c.top):.5*(c.left+c.right));
            std::sort(levels[axis].begin(),levels[axis].end());levels[axis].erase(std::unique(levels[axis].begin(),levels[axis].end()),levels[axis].end());
        }else {const int n=axis?grid.ny:grid.nx;const double origin=axis?grid.bottom:grid.left,step=(axis?grid.height:grid.width)/(2*n);
            for(int k=0;k<=2*n;k++)levels[axis].push_back(origin+k*step);}
    }
    // Native plane/half-chart coordinates suffice on mixed leaves. Never
    // construct a globally finest shadow mesh just to split this curve.
    std::vector<double> cuts;
    for(const auto& arc:source.arcs)for(int axis=0;axis<2;axis++){
        const auto coordinate=arc.coordinate(axis);const auto bound=range(coordinate);
        for(auto level=std::lower_bound(levels[axis].begin(),levels[axis].end(),bound[0]);level!=levels[axis].end()&&*level<=bound[1];++level){
            auto p=coordinate;p[0]-=*level;for(double t:roots(p,0,1))cuts.push_back(arc.a+(arc.b-arc.a)*t);}
    }
    const double begin=source.arcs.front().a,finish=source.arcs.back().b;
    cuts.erase(std::remove_if(cuts.begin(),cuts.end(),[&](double s){return s<begin || s>finish;}),cuts.end());
    for(const auto& arc:source.arcs){
        // The characteristic map resolves the same geometric features as
        // curvature planning and sharp traction. These are exact source
        // restrictions, not new material unknowns or a geometry refit.
        for(double t:surfacePartition(arc,source.boundary,grid.axisymmetric))cuts.push_back(arc.a+(arc.b-arc.a)*t);
        cuts.push_back(arc.a);cuts.push_back(arc.b);
    }
    std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
    const auto mapped=[&](const Arc& arc,double t){V p=arc.sampleFraction(t);
        const bool pole=(source.boundary.axis(0)&&arc.a==0&&t==0)||(source.boundary.axis(1)&&arc.b==1&&t==1);if(pole)p.x=0;
        const bool pinned=(source.boundary.fixed(0)&&arc.a==0&&t==0)||(source.boundary.fixed(1)&&arc.b==1&&t==1);
        if(pinned)p=source.boundary.position[arc.b==1&&t==1];
        MapJet out;try{out=flow({p.x,p.y},dt);}catch(const std::exception& e){std::ostringstream message;message<<std::setprecision(17)<<e.what()<<"; material arc ["<<arc.a<<", "<<arc.b<<"], fraction="<<t<<", pole="<<pole;throw std::runtime_error(message.str());}
        require(!pinned||(out.point.x==p.x&&out.point.y==p.y&&out.timeDerivative.x==0&&out.timeDerivative.y==0),"flow map moves a fixed supply endpoint");
        require(!pole||(out.point.x==0&&out.dx.y==0&&out.dy.x==0),"flow map does not preserve a regular RZ axis pole");return out;};
    const auto make=[&](const Arc& base) {
        // Interpolate the map Jacobian along the source chart and integrate
        // J(R)*R'. This retains the original polynomial first jet exactly
        // for a constant map matrix, including arbitrarily short restrictions.
        // No near-equal mapped endpoint positions are subtracted to infer a
        // derivative. The independent mapped positions check integration error.
        std::array<V,5> columnX{},columnY{};const auto origin=mapped(base,0);
        columnX[0]={origin.dx.x,origin.dx.y};columnY[0]={origin.dy.x,origin.dy.y};
        for(int node=1;node<5;++node) {
            const double t=.25*node;const auto jet=mapped(base,t);const V dx{jet.dx.x-origin.dx.x,jet.dx.y-origin.dx.y},dy{jet.dy.x-origin.dy.x,jet.dy.y-origin.dy.y};
            std::array<double,5> basis{};basis[0]=1;int degree=0;
            for(int other=0;other<5;++other)if(other!=node){std::array<double,5> next{};const double u=.25*other;
                for(int k=0;k<=degree;++k){next[k]-=basis[k]*u/(t-u);next[k+1]+=basis[k]/(t-u);}basis=next;++degree;}
            for(int k=0;k<5;++k){columnX[k]+=dx*basis[k];columnY[k]+=dy*basis[k];}
        }
        PolynomialArc<9> moved;moved.a=base.a;moved.b=base.b;moved.p[0]={origin.point.x,origin.point.y};
        for(int j=0;j<5;++j)for(int k=1;k<5;++k)moved.p[j+k]+=(columnX[j]*base.p[k].x+columnY[j]*base.p[k].y)*(double(k)/(j+k));
        // The source approximation has shared endpoint values AND endpoint
        // tangents as constraints. Interpret the endpoint value in its fixed
        // floating coordinate enclosure, so an unresolvable absolute-position
        // difference does not destroy a retained short chart's first jet.
        const auto end=mapped(base,1);const V represented=moved.sampleFraction(1);
        require(origin.coordinateResolution>=0 && end.coordinateResolution>=0 && std::isfinite(origin.coordinateResolution) && std::isfinite(end.coordinateResolution),"invalid characteristic coordinate resolution");
        const double enclosure=std::max({origin.coordinateResolution,end.coordinateResolution,
            4*std::numeric_limits<double>::epsilon()*std::max({1.,std::abs(end.point.x),std::abs(end.point.y)})});
        const auto outside=[&](double residual){return std::abs(residual)<=enclosure?0.:residual-std::copysign(enclosure,residual);};
        const V constraint{outside(end.point.x-represented.x),outside(end.point.y-represented.y)};
        moved.p[2]+=constraint*3;moved.p[3]+=constraint*(-2);
        return moved;
    };
    std::function<void(const Arc&,int)> append;
    append=[&](const Arc& arc,int depth) {
        const auto moved=make(arc);double positionError=0,tangentError=0;
        for(double t:{.1127016653792583,.276393202250021,.723606797749979,.8872983346207417}) {
            const auto jet=mapped(arc,t);const V exact{jet.point.x,jet.point.y},tangent=applyJet(jet,arc.sampleFraction(t,1));
            positionError=std::max(positionError,norm(moved.sampleFraction(t)-exact));
            tangentError=std::max(tangentError,norm(moved.sampleFraction(t,1)-tangent)/std::max(1.,norm(tangent)));
        }
        if(positionError<=5e-13 && tangentError<=2e-9) {
            diagnostics.positionApproximation=std::max(diagnostics.positionApproximation,positionError);diagnostics.tangentApproximation=std::max(diagnostics.tangentApproximation,tangentError);out.arcs.push_back(moved);return;
        }
        if(depth>=20){std::cerr<<std::setprecision(17)<<"NATIVE_CURVE_APPROXIMATION a="<<arc.a<<" b="<<arc.b<<" position="<<positionError<<" tangent="<<tangentError<<std::endl;throw std::runtime_error("first-jet transport approximation unresolved");}
        const double middle=.5*(arc.a+arc.b);
        if(!(middle>arc.a && middle<arc.b))std::cerr<<std::setprecision(17)<<"NATIVE_CURVE_PARAMETER a="<<arc.a<<" b="<<arc.b<<" position="<<positionError<<" tangent="<<tangentError<<std::endl;
        require(middle>arc.a && middle<arc.b,"transport refinement exhausted material parameter");++diagnostics.subdivisions;
        append(restrictedArc(arc,arc.a,middle),depth+1);append(restrictedArc(arc,middle,arc.b),depth+1);
    };
    for(std::size_t i=1;i<cuts.size();++i) {
        const Arc base=restrictedArc(source.arcs[source.index(cuts[i-1])],cuts[i-1],cuts[i]);
        // Split also where the transported endpoint chart crosses a native
        // velocity branch. This is a numerical simple-crossing control; it
        // is not a certification of arbitrary tangencies or critical events.
        std::vector<double> pieces{0,1};std::array<MapJet,9> samples;
        for(int k=0;k<=8;++k)samples[k]=mapped(base,k/8.);
        for(int axis=0;axis<2;++axis)for(int panel=0;panel<8;++panel) {
            const double v0=axis?samples[panel].point.y:samples[panel].point.x,v1=axis?samples[panel+1].point.y:samples[panel+1].point.x;
            const auto first=std::lower_bound(levels[axis].begin(),levels[axis].end(),std::min(v0,v1));
            for(auto edge=first;edge!=levels[axis].end()&&*edge<=std::max(v0,v1);++edge) {
                const double level=*edge;double lo=panel/8.,hi=(panel+1)/8.,fl=v0-level,fh=v1-level;
                if(fl==0){pieces.push_back(lo);continue;}if(fh==0){pieces.push_back(hi);continue;}if((fl<0)==(fh<0))continue;
                for(int n=0;n<52;++n){const double middle=.5*(lo+hi);if(middle==lo || middle==hi)break;const auto value=mapped(base,middle);const double fm=(axis?value.point.y:value.point.x)-level;
                    if((fm<0)==(fl<0)){lo=middle;fl=fm;}else hi=middle;}
                pieces.push_back(.5*(lo+hi));
            }
        }
        std::sort(pieces.begin(),pieces.end());pieces.erase(std::unique(pieces.begin(),pieces.end()),pieces.end());
        std::vector<double> parameters{base.a};for(double p:pieces)if(p>0 && p<1){const double s=base.a+(base.b-base.a)*p;if(s>base.a && s<base.b)parameters.push_back(s);}parameters.push_back(base.b);
        std::sort(parameters.begin(),parameters.end());parameters.erase(std::unique(parameters.begin(),parameters.end()),parameters.end());
        for(std::size_t j=1;j<parameters.size();++j)append(restrictedArc(base,parameters[j-1],parameters[j]),0);
    }
    diagnostics.maximumSourcePieces=std::max(diagnostics.maximumSourcePieces,out.arcs.size());
    for(std::size_t i=0;i<out.arcs.size();++i){if(i+1==out.arcs.size() && (source.boundary.open()||begin!=0 || finish!=1))break;const auto& a=out.arcs[i];const auto& b=out.arcs[(i+1)%out.arcs.size()];
        diagnostics.sourceJoin=std::max(diagnostics.sourceJoin,norm(a.sampleFraction(1)-b.sampleFraction(0)));
        diagnostics.sourceTangentJoin=std::max(diagnostics.sourceTangentJoin,norm(a.sampleFraction(1,1)-b.sampleFraction(0,1)));}
    return out;
}
Step prepare(const Candidate& old,double dt,const Flow& flow,const NativeFlux& flux,
    const InterfaceParametricProjection::Options& options,const Observer& observer) {
    return prepare(old,dt,flow,TimeLinearFlux{flux,flux},options,observer);
}
Step prepare(const Candidate& old,double dt,const Flow& flow,const TimeLinearFlux& flux,
    const InterfaceParametricProjection::Options& options,const Observer& observer) {
    require(old.grid().sameLayout(flux.begin.grid)&&old.grid().sameLayout(flux.end.grid),"parametric flux/inventory physical grid or native leaf ownership differs");
    require(dt!=0&&std::isfinite(dt),"invalid parametric step duration");
    if(!old.componentCount()){
        InterfaceParametricProjection::activeOptions(old,options);FluxReport report;FaceVolumes volumes;
        backgroundFlux(flux,dt,report,volumes);
        // Frozen signed remainders remain in old.inventory(): no curve, no
        // liquid flux, no artificial redistribution and no new source chart.
        return {old,{},report,0,std::move(volumes)};
    }
    if(old.componentCount()>1||!old.retiredComponents().empty()){
        const auto local=InterfaceParametricProjection::activeOptions(old,options);
        std::vector<Candidate> components;std::vector<std::uint64_t> ids;std::vector<Step> steps;
        for(std::size_t i=0;i<old.componentCount();i++){
            steps.push_back(prepare(old.component(i),dt,flow,flux,local,observer));components.push_back(steps.back().candidate);ids.push_back(old.boundaryHandle()->components()[i].id);
        }
        auto candidate=InterfaceParametricProjection::combine(old.grid(),std::move(components),ids,old.retiredComponents());
        require(candidate.report().maximumInventoryFractionError<=options.inventoryFractionUncertainty+options.measureFractionReserve,
            "active transport and retired inventory exceed the common fraction budget");
        require(std::abs(candidate.report().totalGeometryInventoryDifference)<=options.measureFractionReserve*old.grid().totalCapacity(),
            "active transport and retired inventory exceed the common total-volume reserve");
        auto out=std::move(steps.front());out.candidate=std::move(candidate);
        for(std::size_t i=1;i<steps.size();i++){
            const auto& b=steps[i];out.inventoryBalance+=b.inventoryBalance;
            auto& surface=out.faceVolumes.surface;const auto& other=b.faceVolumes.surface;
            surface.advection+=other.advection;surface.endpointWork+=other.endpointWork;surface.quadratureError+=other.quadratureError;
            surface.criticalError+=other.criticalError;surface.budget+=other.budget;surface.timeCharts+=other.timeCharts;
            surface.materialCharts+=other.materialCharts;surface.refinements+=other.refinements;surface.criticalPatches+=other.criticalPatches;
            auto& s=out.transport;const auto& t=b.transport;
            s.positionApproximation=std::max(s.positionApproximation,t.positionApproximation);s.tangentApproximation=std::max(s.tangentApproximation,t.tangentApproximation);
            s.sourceJoin=std::max(s.sourceJoin,t.sourceJoin);s.sourceTangentJoin=std::max(s.sourceTangentJoin,t.sourceTangentJoin);s.maximumSourcePieces+=t.maximumSourcePieces;s.subdivisions+=t.subdivisions;
            auto& r=out.flux;const auto& q=b.flux;
            r.crossingIntervals+=q.crossingIntervals;r.vertexEvents+=q.vertexEvents;r.tangentTurns+=q.tangentTurns;r.quadratureRefinements+=q.quadratureRefinements;
            r.spaceTimePatches+=q.spaceTimePatches;r.timeCharts+=q.timeCharts;r.materialCharts+=q.materialCharts;r.chartSubdivisions+=q.chartSubdivisions;r.boundedVertexIntervals+=q.boundedVertexIntervals;r.boundedContactPatches+=q.boundedContactPatches;r.contactEdgeRefinements+=q.contactEdgeRefinements;r.mapQueries+=q.mapQueries;r.mapCacheHits+=q.mapCacheHits;
            r.crossingResidual=std::max(r.crossingResidual,q.crossingResidual);r.sheetPositionError=std::max(r.sheetPositionError,q.sheetPositionError);
            r.quadratureError+=q.quadratureError;r.exteriorFlux+=q.exteriorFlux;r.vertexIntervalError+=q.vertexIntervalError;r.contactIntegralError+=q.contactIntegralError;r.directionalQuadratureError+=q.directionalQuadratureError;r.directionalBoundAdjustment+=q.directionalBoundAdjustment;r.momentQuadratureError+=q.momentQuadratureError;r.quadratureFractionError+=q.quadratureFractionError;r.directionalFractionError+=q.directionalFractionError;r.totalFluxFractionError+=q.totalFluxFractionError;
            for(int axis=0;axis<2;axis++)for(std::size_t k=0;k<out.faceVolumes.faces[axis].size();k++){
                auto& a=out.faceVolumes.faces[axis][k];const auto& f=b.faceVolumes.faces[axis][k];
                require(a.total==f.total&&a.totalPositive==f.totalPositive&&a.totalNegative==f.totalNegative&&a.totalHigher==f.totalHigher&&a.totalPositiveHigher==f.totalPositiveHigher,"components do not share one total face flux");
                require(a.totalProducts==f.totalProducts&&a.totalPositiveProducts==f.totalPositiveProducts,"components do not share one weighted total flux");
                a.liquid+=f.liquid;a.liquidPositive+=f.liquidPositive;a.liquidNegative+=f.liquidNegative;
                for(int order=0;order<2;order++){a.liquidHigher[order]+=f.liquidHigher[order];a.liquidPositiveHigher[order]+=f.liquidPositiveHigher[order];}
                for(int k=0;k<=FaceVolumes::productDegree;k++){a.liquidProducts[k]+=f.liquidProducts[k];a.liquidPositiveProducts[k]+=f.liquidPositiveProducts[k];}
            }
        }
        return out;
    }
    const auto& grid=old.grid();
    // Optional, synchronous observations of the unchanged characteristic map.
    // Sum of map time may exceed wall time when flux charts run in parallel.
    using Clock=std::chrono::steady_clock;
    const bool profile=std::getenv("QUANTUM_MOMENT_STEP_PROFILE")!=nullptr;
    auto last=Clock::now();std::atomic<long long> mapNanoseconds{0},mapCalls{0};
    Flow measured=flow;
    if(profile)measured=[&](Vector2D p,double duration){const auto start=Clock::now();const auto result=flow(p,duration);
        mapNanoseconds.fetch_add(std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now()-start).count(),std::memory_order_relaxed);
        mapCalls.fetch_add(1,std::memory_order_relaxed);return result;};
    const auto mark=[&](const char* stage){if(profile){const auto now=Clock::now();std::ostringstream out;
        out<<"[moment_motion_profile] dt="<<dt<<" work="<<stage<<" seconds="<<std::chrono::duration<double>(now-last).count()
            <<" map_seconds="<<1e-9*mapNanoseconds.exchange(0)<<" map_calls="<<mapCalls.exchange(0)<<'\n';
        std::cout<<out.str()<<std::flush;last=now;}};
    SourceReport transport;FluxReport report;
    const auto source=transportSource(old.geometry(),grid,dt,measured,transport);
    mark("source");
    auto projection=options;
    projection.minimumMaterialKnots=InterfaceParametricProjection::resolveMaterialSpace(source,grid,options);
    mark("free_geometry");
    FaceVolumes volumes;
    const auto inventory=advanceInventory(old.geometry(),source,measured,flux,dt,old.inventory(),report,&volumes);
    mark("inventory_flux");
    auto candidate=InterfaceParametricProjection::prepare(source,grid,inventory,projection);
    mark("projection");
    if(old.boundaryHandle()->components()[0].id!=1)candidate=InterfaceParametricProjection::combine({std::move(candidate)},{old.boundaryHandle()->components()[0].id});
    const double before=std::accumulate(old.inventory().begin(),old.inventory().end(),0.);
    const double after=std::accumulate(inventory.begin(),inventory.end(),0.);
    const double balance=after-before+report.exteriorFlux;
    require(std::abs(balance)<=64*std::numeric_limits<double>::epsilon()*std::max({1.,std::abs(before),std::abs(after),std::abs(report.exteriorFlux)}),
        "shared native phase fluxes do not close the global inventory");
    if(observer)observer(source,candidate);
    return {std::move(candidate),transport,report,balance,std::move(volumes)};
}
} // namespace InterfaceParametricTransport
