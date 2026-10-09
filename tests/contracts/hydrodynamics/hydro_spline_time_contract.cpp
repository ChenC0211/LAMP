#include "Modules/Hydrodynamics/HydroSplineVelocity.h"
#include "Modules/Hydrodynamics/HydroInterfaceTransport.h"
#include "Modules/Interface/Compatibility/InterfaceParametricMotion.h"
#include "Modules/Interface/Compatibility/InterfacePrescribedData.h"
#include "Core/Geometry/LevelSetGaussTrace.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Field/VEdgeField.h"
#include "Core/Field/HEdgeField.h"
#include <cmath>
#include <iostream>
#include <iomanip>

namespace {
void require(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
void exactTimeMap() {
    // Noncommuting affine fields: u=( (1+t)y, 1 ). The exact flow is cubic
    // in elapsed time, so its derivative tests the stage-time terms directly.
    const auto query=[](Vector2D p,double t){return LevelSetGaussTrace::TimeVelocityJet{{(1+t)*p.y,1},{0,0},{1+t,0},{p.y,0}};};
    double error=0,fdError=0;
    for(int blocks:{1,8})for(double t:{0.,.17,.4}) {
        const Vector2D p{.23,.31};const auto map=LevelSetGaussTrace::advanceTimeDependent(query,{{-2,-2},{2,2}},p,t,blocks);
        const Vector2D exact{p.x+p.y*(t+.5*t*t)+.5*t*t+t*t*t/3,p.y+t};
        error=std::max({error,(map.point-exact).mod(),(map.timeDerivative-Vector2D{(1+t)*(p.y+t),1}).mod(),
            std::abs(map.derivative(0,0)-1),std::abs(map.derivative(1,0)),std::abs(map.derivative(0,1)-t-.5*t*t),std::abs(map.derivative(1,1)-1)});
        if(t>0){constexpr double h=1e-6;const auto a=LevelSetGaussTrace::advanceTimeDependent(query,{{-2,-2},{2,2}},p,t-h,blocks);
            const auto b=LevelSetGaussTrace::advanceTimeDependent(query,{{-2,-2},{2,2}},p,t+h,blocks);
            fdError=std::max(fdError,((b.point-a.point)*(.5/h)-map.timeDerivative).mod());}
    }
    require(error<3e-13 && fdError<2e-9,"physical-time collocation lost the exact noncommuting flow or map derivative");
    std::cout<<std::setprecision(17)<<"HYDRO_SPLINE_TIME_MAP exact_jet_error="<<error<<" duration_fd_error="<<fdError<<'\n';
}

void returningFlux() {
    using namespace InterfaceParametricTransport;
    using namespace Quantum::Geometry::Parametric;
    constexpr int n=8;constexpr double duration=.2;
    const auto curve=InterfacePrescribedData::parametricEllipse({.45,.5},{.15,.15},16);
    const auto old=InterfaceParametricProjection::initialize(curve,{n,n});
    const Vector2D start{.7,.2},end{-.7,.3};
    const Flow flow=[=](Vector2D p,double t){const auto move=start*t+(end-start)*(.5*t*t/duration);
        return MapJet{p+move,{1,0},{0,1},start+(end-start)*(t/duration),move};};
    const auto fluxFor=[&](Vector2D v){NativeFlux f;f.grid.nx=f.grid.ny=n;
        f.integral=[v](bool horizontal,int,int,double a,double b){return (horizontal?v.y:v.x)*(b-a);};
        f.density=[v](bool horizontal,int,int,double){return horizontal?v.y:v.x;};
        f.densityMagnitudeBound=[v](bool horizontal,int,int){return std::abs(horizontal?v.y:v.x);};return f;};
    SourceReport sr;FluxReport fr;const auto source=transportSource(curve,{n,n},duration,flow,sr);
    const auto ledger=advanceInventory(curve,source,flow,TimeLinearFlux{fluxFor(start),fluxFor(end)},duration,old.inventory(),fr);
    auto exact=curve;const auto shift=(start+end)*(.5*duration);
    for(auto& arc:exact.arcs)arc.p[0]+=Point{shift.x,shift.y};
    PhaseMeasure measure(exact);double fraction=0,balance=0;
    for(int i=0;i<n;++i)for(int j=0;j<n;++j){const int k=i*n+j;
        const double volume=measure.simpleRectangle(double(i)/n,double(i+1)/n,double(j)/n,double(j+1)/n);
        fraction=std::max(fraction,std::abs(volume-ledger[k])*n*n);balance+=ledger[k]-old.inventory()[k];}
    InterfaceParametricProjection::Options options;const auto candidate=InterfaceParametricProjection::prepare(source,{n,n},ledger,options);
    require(fraction<3e-12 && std::abs(balance)<2e-15 && candidate.report().maximumInventoryFractionError<=3e-12,
        "time-linear returning flux disagrees with the independent translated geometry");
    std::cout<<"HYDRO_SPLINE_TIME_RETURN native_fraction_error="<<fraction<<" projected_fraction_error="<<candidate.report().maximumInventoryFractionError
        <<" balance="<<balance<<" within_step_normal_reversal=1 time_charts="<<fr.timeCharts<<" material_charts="<<fr.materialCharts<<'\n';
}

void stationaryContactFlux() {
    using namespace InterfaceParametricTransport;
    using namespace Quantum::Geometry::Parametric;
    constexpr int n=8;constexpr double duration=.1;
    // The circle touches x=.5 and y=.5 at rest. Acceleration either removes
    // the contact or creates two crossings; exchanging x/y checks orientation.
    // The independently translated original curve supplies every cell volume.
    const auto curve=InterfacePrescribedData::parametricEllipse({.65,.65},{.15,.15},16);
    const auto old=InterfaceParametricProjection::initialize(curve,{n,n});
    for(const Vector2D acceleration:{Vector2D{0,2},Vector2D{0,-2},Vector2D{2,0},Vector2D{-2,0}}) {
        const Flow flow=[=](Vector2D p,double t){const auto move=acceleration*(.5*t*t);
            return MapJet{p+move,{1,0},{0,1},acceleration*t,move};};
        const auto fluxFor=[&](double time){NativeFlux f;f.grid.nx=f.grid.ny=n;const auto v=acceleration*time;
            f.integral=[v](bool horizontal,int,int,double a,double b){return (horizontal?v.y:v.x)*(b-a);};
            f.density=[v](bool horizontal,int,int,double){return horizontal?v.y:v.x;};
            f.densityMagnitudeBound=[v](bool horizontal,int,int){return std::abs(horizontal?v.y:v.x);};return f;};
        SourceReport sr;FluxReport fr;const auto source=transportSource(curve,{n,n},duration,flow,sr);
        const auto ledger=advanceInventory(curve,source,flow,TimeLinearFlux{fluxFor(0),fluxFor(duration)},duration,old.inventory(),fr);
        auto exact=curve;const auto shift=acceleration*(.5*duration*duration);
        for(auto& arc:exact.arcs)arc.p[0]+=Point{shift.x,shift.y};
        PhaseMeasure measure(exact);double fraction=0,balance=fr.exteriorFlux;
        for(int i=0;i<n;++i)for(int j=0;j<n;++j){const int k=i*n+j;
            const double volume=measure.simpleRectangle(double(i)/n,double(i+1)/n,double(j)/n,double(j+1)/n);
            fraction=std::max(fraction,std::abs(volume-ledger[k])*n*n);balance+=ledger[k]-old.inventory()[k];}
        require(fraction<3e-12 && std::abs(balance)<2e-15 && fr.boundedContactPatches>0 && fr.contactIntegralError<=1e-13/(n*n),
            "stationary contact flux disagrees with the independent translated geometry or its declared budget");
        std::cout<<"HYDRO_SPLINE_CONTACT acceleration="<<acceleration.x<<','<<acceleration.y<<" native_fraction_error="<<fraction
            <<" balance="<<balance<<" contact_patches="<<fr.boundedContactPatches<<" contact_integral_error="<<fr.contactIntegralError
            <<" edge_refinements="<<fr.contactEdgeRefinements<<'\n';
    }
}

void splineTimeRefinement(bool fullInterface) {
    DictionaryMesh input;input.minCoord={0,0};input.maxCoord={1,1};input.Nx=input.Ny=8;input.maxCellNum=1024;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,false);
    const double pi=std::acos(-1.);constexpr double finish=.4;
    const auto snapshot=[&](double a,double b){VEdgeField<double> ux(mesh,0);HEdgeField<double> uy(mesh,0);
        const auto stream=[&](Vector2D p){return .3*(a*std::sin(pi*p.x)*std::sin(pi*p.y)/pi+b*std::sin(2*pi*p.x)*std::sin(pi*p.y)/(2*pi));};
        for(const auto face:domain.faces(mesh)){const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
            const auto lo=mesh.verts(edge.edgeToVert[0]).coord(),hi=mesh.verts(edge.edgeToVert[1]).coord();
            (face.horizontal?uy.value:ux.value)[edge.index()]=face.outwardSign?0:(face.horizontal?-1:1)*(stream(hi)-stream(lo))/face.area;}
        return HydroSplineVelocity(mesh,domain,ux,uy);};
    const auto a=snapshot(1,0),b=snapshot(0,1);
    const auto amplitude=[](double t){return Vector2D{1+t*t,.4*std::sin(3*t)};};
    const auto velocity=[&](Vector2D p,double t){const auto weights=amplitude(t);return a.firstJet(p).value*weights.x+b.firstJet(p).value*weights.y;};
    const std::vector<Vector2D> seeds{{.23,.31},{.63,.27},{.37,.67},{.72,.71},{.49,.57},{.81,.42}};
    const auto reference=[&](int steps,std::vector<Vector2D> points){const double h=finish/steps;
        for(int k=0;k<steps;++k)for(auto& p:points){const double t=k*h;const auto k1=velocity(p,t),k2=velocity(p+k1*(.5*h),t+.5*h),
            k3=velocity(p+k2*(.5*h),t+.5*h),k4=velocity(p+k3*h,t+h);p+=(k1+k2*2+k3*2+k4)*(h/6);}
        return points;};
    const auto ref=reference(4096,seeds),coarseRef=reference(2048,seeds);double referenceDifference=0;
    for(std::size_t i=0;i<ref.size();++i)referenceDifference=std::max(referenceDifference,(ref[i]-coarseRef[i]).mod());
    require(referenceDifference<2e-11,"independent time reference is unresolved");
    double previousError=0;
    for(int steps:{4,8,16}) {auto points=seeds;const double h=finish/steps;double jetError=0;
        const auto initial=InterfacePrescribedData::parametricEllipse({.5,.65},{.15,.15},16);
        auto state=InterfaceParametricMotion::initialize(mesh,domain,initial);
        InterfaceParametricProjection::Options options;double maximumFraction=0,maximumBalance=0,areaDrift=0;
        for(int k=0;k<steps;++k){const auto w0=amplitude(k*h),w1=amplitude((k+1)*h);const auto old=snapshot(w0.x,w0.y),next=snapshot(w1.x,w1.y);
            const auto flow=HydroInterfaceTransport::timeLinearParametricFlow(old,next,h);
            if(k==0){const auto p=seeds[0];constexpr double e=1e-6;const auto m=flow(p,h*.7);
                jetError=std::max({(m.dx-(flow(p+Vector2D{e,0},h*.7).point-flow(p-Vector2D{e,0},h*.7).point)*(.5/e)).mod(),
                    (m.dy-(flow(p+Vector2D{0,e},h*.7).point-flow(p-Vector2D{0,e},h*.7).point)*(.5/e)).mod(),
                    (m.timeDerivative-(flow(p,h*.7+e).point-flow(p,h*.7-e).point)*(.5/e)).mod()});}
            for(auto& p:points)p=flow(p,h).point;
            if(steps==16 && fullInterface){const auto firstFlux=HydroInterfaceTransport::normalFlux(old),lastFlux=HydroInterfaceTransport::normalFlux(next);
                const auto inventory=state.parametricGeometry().inventory();
                if(k==0){const auto p=initial.sample(.75);const auto m0=flow({p.x,p.y},0),m1=flow({p.x,p.y},h/8);const auto tangent=initial.sample(.75,1);
                    std::cout<<"HYDRO_SPLINE_TIME_VERTEX point="<<p.x<<','<<p.y<<" material_tangent="<<tangent.x<<','<<tangent.y
                        <<" initial_velocity="<<m0.timeDerivative.x<<','<<m0.timeDerivative.y<<" elapsed="<<h/8
                        <<" moved="<<m1.displacement.x<<','<<m1.displacement.y<<" current_velocity="<<m1.timeDerivative.x<<','<<m1.timeDerivative.y<<std::endl;}
                std::cout<<"HYDRO_SPLINE_TIME_INTERFACE_BEGIN step="<<k+1<<" begin_time="<<k*h<<" end_time="<<(k+1)*h<<std::endl;
                auto candidate=InterfaceParametricMotion::prepareTimeLinear(mesh,domain,state,{k+1,k*h,(k+1)*h,h},flow,firstFlux,lastFlux,options);
                require(state.progress().completedSteps==k && state.parametricGeometry().inventory()==inventory,"time-varying multi-step update modified the old state");
                maximumFraction=std::max(maximumFraction,candidate.state.parametricGeometry().report().maximumInventoryFractionError);
                maximumBalance=std::max(maximumBalance,std::abs(candidate.inventoryBalance));
                areaDrift=std::max(areaDrift,std::abs(candidate.state.parametricGeometry().geometry().area()-initial.area()));
                state=std::move(candidate.state);}
        }
        double error=0;for(std::size_t i=0;i<points.size();++i){const auto d=points[i]-ref[i];error+=d.x*d.x+d.y*d.y;}error=std::sqrt(error/points.size());
        require(jetError<2e-8 && (!previousError || previousError/error>3.7 && previousError/error<4.3),"two-layer physical-time transport did not show second-order refinement");
        std::cout<<"HYDRO_SPLINE_TIME_ACCURACY steps="<<steps<<" end_time="<<finish<<" position_rms="<<error<<" ratio="<<(previousError?previousError/error:0)
            <<" reference_difference="<<referenceDifference<<" map_jet_fd_error="<<jetError<<'\n';previousError=error;
        if(steps==16 && fullInterface){std::vector<Vector2D> boundary;for(int i=0;i<64;++i){const auto p=initial.sample(double(i)/64);boundary.push_back({p.x,p.y});}
            const auto exact=reference(4096,boundary),check=reference(2048,boundary);double shape=0,boundaryReference=0;
            for(int i=0;i<64;++i){const auto p=state.parametricGeometry().geometry().sample(double(i)/64);const auto d=Vector2D{p.x,p.y}-exact[i];
                shape+=d.x*d.x+d.y*d.y;boundaryReference=std::max(boundaryReference,(exact[i]-check[i]).mod());}
            shape=std::sqrt(shape/64);
            require(state.progress().completedSteps==16 && state.progress().endTime==finish && maximumFraction<=3e-12 && maximumBalance<2e-15
                && areaDrift<2e-10 && boundaryReference<2e-11 && shape<1e-4,"complete time-varying interface trajectory exceeded its original geometry/inventory budgets");
            std::cout<<"HYDRO_SPLINE_TIME_INTERFACE steps=16 end_time="<<finish<<" maximum_fraction_error="<<maximumFraction<<" inventory_balance="<<maximumBalance
                <<" maximum_area_drift="<<areaDrift<<" material_position_rms="<<shape<<" reference_difference="<<boundaryReference<<" curvature_order_claim=0\n";}
    }
}
}
void hydroSplineTimeContract(bool fullInterface){exactTimeMap();returningFlux();stationaryContactFlux();splineTimeRefinement(fullInterface);}
