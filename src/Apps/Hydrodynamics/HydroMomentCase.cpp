#include "Hydrodynamics/HydroMomentCase.h"
#include "Modules/Hydrodynamics/HydrodynamicsConfig.h"
#include "Modules/Hydrodynamics/HydroCompatibleMomentTraction.h"
#include "Modules/Interface/Compatibility/InterfaceParametricMotion.h"
#include "Modules/Interface/Compatibility/InterfacePrescribedData.h"
#include "Modules/Interface/Compatibility/IO/DictionaryInterface.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/IO/OutputLayout.h"
#include "Core/IO/TecplotWriter.h"
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <map>
#include <limits>

namespace {
using Flow=Quantum::Hydrodynamics::HydroCompatibleMomentFlow;
using Space=Quantum::Hydrodynamics::HydroCompatibleMomentSpace;
using Stress=Quantum::Hydrodynamics::HydroCompatibleMomentViscosity;
std::string number(double v){std::ostringstream s;s<<std::setprecision(17)<<v;return s.str();}
std::string stem(const std::string& prefix,int step){std::ostringstream s;s<<prefix<<"_step"<<std::setw(4)<<std::setfill('0')<<step;return s.str();}
const char* neckStatus(Flow::NeckReadiness status){
    switch(status){case Flow::NeckReadiness::ScalePending:return "scale_pending";case Flow::NeckReadiness::CapPending:return "cap_pending";
    case Flow::NeckReadiness::MomentumPending:return "momentum_pending";case Flow::NeckReadiness::EnergyPending:return "energy_pending";case Flow::NeckReadiness::Ready:return "ready";}
    throw std::logic_error("unknown coupled neck readiness");
}
const char* transferStatus(Flow::MeshTransferReadiness status){
    switch(status){case Flow::MeshTransferReadiness::MomentumPending:return "momentum_pending";
    case Flow::MeshTransferReadiness::EnergyPending:return "energy_pending";
    case Flow::MeshTransferReadiness::AccuracyPending:return "accuracy_pending";
    case Flow::MeshTransferReadiness::Ready:return "ready";}
    throw std::logic_error("unknown mesh transfer readiness");
}
}
struct HydroMomentCase::Impl {
    std::unique_ptr<NativeLeafDomain> domain;
    std::unique_ptr<Flow> flow;
    Flow::Options options;
    Vector2D acceleration;
    double viscosity=0,mass=0,initialEnergy=0,transport=0,viscous=0,external=0,pressure=0,central=0;
    bool twoPhase=false;
    double sigma=0,initialSurface=0,initialLiquid=0,capillary=0,temporal=0,metric=0;
    double surfaceAdvection=0,surfaceEndpoint=0,surfaceFluxError=0,surfaceFluxBudget=0,surfaceCriticalError=0;
    double meshEnergyChange=0,geometricWork=0;
    double initialPotential=std::numeric_limits<double>::quiet_NaN();
    std::array<double,2> initialPhaseInventory{};
    double neckRadiusCells=0,eventKineticChange=0,eventSurfaceChange=0,eventDissipation=0;
    double boundaryTraction=0,boundaryReaction=0,boundaryKinetic=0,boundaryDissipation=0;
    std::ofstream history,curves,adaptationHistory,coarseningHistory,componentHistory,phaseHistory,eventHistory,neckHistory,resolutionHistory,terminalHistory;
    std::vector<Vector2D> protectedBoundaryPoints;
    std::optional<InterfaceAdaptationControls> adaptation;
    InterfaceResolutionPlan::Result resolution;
    int resolutionReview=0;
    InterfaceResolutionPlan::Result review(const QuadtreeMesh&,const InterfaceState&,const char*);
    std::map<std::uint64_t,double> initialComponentVolumes;
    int writtenStep=-1;
};
HydroMomentCase::HydroMomentCase(std::unique_ptr<QuadtreeMesh>& meshOwner,const HydrodynamicsConfig& config,
    const std::optional<HydroGeometryInput::ParametricLoop>& shape,const InterfaceAdaptationControls* adaptation):impl_(std::make_unique<Impl>()){
    auto& p=*impl_;const bool rz=config.coordinateSystem=="axisymmetricRZ";
    p.twoPhase=bool(shape);p.sigma=config.surfaceTension;p.neckRadiusCells=config.momentNeckRadiusCells;
    if(adaptation)p.adaptation=*adaptation;
    for(const auto& patch:config.boundaryPatches)for(double s:{patch.start,patch.end}){
        switch(patch.side){
        case HydroBoundarySide::West:p.protectedBoundaryPoints.push_back({meshOwner->minCoord.x,s});break;
        case HydroBoundarySide::East:p.protectedBoundaryPoints.push_back({meshOwner->maxCoord.x,s});break;
        case HydroBoundarySide::South:p.protectedBoundaryPoints.push_back({s,meshOwner->minCoord.y});break;
        case HydroBoundarySide::North:p.protectedBoundaryPoints.push_back({s,meshOwner->maxCoord.y});break;
        }
    }
    FaceSegmentGraph graph;graph.build(*meshOwner);p.domain=std::make_unique<NativeLeafDomain>(NativeLeafDomain::build(*meshOwner,graph,rz));
    std::optional<InterfaceState> initialInterface;
    if(shape){
        const auto make=[&](const HydroGeometryInput::ParametricShape& s){
            if(s.analytic)return InterfacePrescribedData::parametricAnalytic(Quantum::ScalarExpression(s.firstExpression,s.constants),
                Quantum::ScalarExpression(s.secondExpression,s.constants),s.segments,rz);
            if(s.implicit){
                const auto lo=meshOwner->minCoord,hi=meshOwner->maxCoord;
                const double tolerance=s.traceTolerance>0?s.traceTolerance:1e-9*std::min(hi.x-lo.x,hi.y-lo.y);
                return InterfacePrescribedData::parametricImplicit(Quantum::ScalarExpression(s.implicitExpression,s.constants,
                    rz?Quantum::ScalarExpression::Coordinates::Meridional:Quantum::ScalarExpression::Coordinates::Planar),
                    {lo.x,hi.x,lo.y,hi.y},s.seed,s.segments,rz,tolerance);
            }
            if(s.necked)return InterfacePrescribedData::parametricNeckedSpheroid(s.center.y,s.radii,s.neckRatio,s.segments);
            if(s.suppliedCap)return InterfacePrescribedData::parametricSupplyCap(s.radii.x,s.radii.y,s.segments,meshOwner->minCoord.y);
            if(rz)return InterfacePrescribedData::parametricSpheroid(s.center.y,s.radii,s.segments);
            return InterfacePrescribedData::parametricEllipse(s.center,s.radii,s.segments);
        };
        std::vector<InterfaceParametricProjection::Curve> curves;
        if(shape->components.empty())curves.push_back(make(*shape));else for(const auto& component:shape->components)curves.push_back(make(component));
        initialInterface=InterfaceParametricMotion::initialize(*meshOwner,*p.domain,std::move(curves));
        p.options.interfaceProjection.materialSegments=shape->segments;p.options.interfaceProjection.relativeH2=.001*std::pow(16./shape->segments,2);
        if(p.adaptation){
            for(int sweep=0;;sweep++){
                p.resolution=p.review(*meshOwner,*initialInterface,"initial");
                if(p.resolution.saturated)throw std::length_error("initial moment resolution request exceeds maximumLevel outside terminal scope");
                if(!p.resolution.changesMesh())break;
                if(sweep>p.adaptation->policy.maximumLevel)throw std::runtime_error("initial local resolution did not reach a fixed mesh");
                auto target=InterfaceResolutionPlan::refine(*meshOwner,p.resolution,p.adaptation->cellBudget);
                FaceSegmentGraph nextGraph;nextGraph.build(*target);auto domain=std::make_unique<NativeLeafDomain>(NativeLeafDomain::build(*target,nextGraph,rz));
                auto next=InterfaceParametricMotion::prepareMeshTransfer(*meshOwner,*p.domain,*initialInterface,*target,*domain,p.options.interfaceProjection);
                initialInterface=std::move(next.state);meshOwner.swap(target);p.domain.swap(domain);
            }
        }
    }else if(p.adaptation)throw std::invalid_argument("initial moment adaptation requires a parametric interface");
    // Interpolate the prescribed physical initial velocity only AFTER the
    // initial geometry/mesh is resolved; no coarse fluid predictor is taken.
    const auto& mesh=*meshOwner;auto space=Space::build(mesh,*p.domain);
    const double width=mesh.maxCoord.x-mesh.minCoord.x,height=mesh.maxCoord.y-mesh.minCoord.y,length=std::min(width,height);
    const auto initial=space.interpolateVelocity(mesh,*p.domain,[&](Vector2D v)->Eigen::Vector2d{
        if(config.initialVelocity=="zero")return Eigen::Vector2d::Zero();
        const double x=(v.x-mesh.minCoord.x)/width,y=(v.y-mesh.minCoord.y)/height;
        const double f=x*x*(1-x)*(1-x),g=y*y*(1-y)*(1-y),fx=2*x*(1-x)*(1-2*x),gy=2*y*(1-y)*(1-2*y);
        if(rz){const double x2=x*x;
            return {8*config.velocityAmplitude*length/height*x*(1-x2)*(1-x2)*gy,
                8*config.velocityAmplitude*length/width*(-2+8*x2-6*x2*x2)*g};}
        return {32*config.velocityAmplitude*length/height*f*gy,-32*config.velocityAmplitude*length/width*fx*g};
    });
    p.acceleration={config.accelerationX,config.accelerationY};
    if(shape){
        auto interface=std::move(*initialInterface);
        p.initialSurface=p.sigma*Quantum::Hydrodynamics::HydroCompatibleMomentTraction::measure(*interface.parametricGeometry().boundaryHandle(),mesh,*p.domain).length;
        for(double v:interface.parametricGeometry().inventory())p.initialLiquid+=v;
        for(std::size_t i=0;i<interface.parametricGeometry().componentCount();i++){double volume=0;for(double v:interface.parametricGeometry().componentInventory(i))volume+=v;p.initialComponentVolumes.emplace(interface.parametricGeometry().boundaryHandle()->components()[i].id,volume);}
        p.componentHistory.open(OutputLayout::series()/"hydro_moment_components.csv");if(!p.componentHistory)throw std::runtime_error("cannot open component history");
        p.componentHistory<<std::setprecision(17)<<"step,time,component,liquid_volume,volume_change,geometry_volume,arcs,status,retired_at\n";
        p.flow=std::make_unique<Flow>(Flow::createTwoPhase(mesh,*p.domain,std::move(space),std::move(interface),
            Quantum::Hydrodynamics::HydroCompatibleMomentMaterial::Properties{config.densityNegative,config.densityPositive,config.viscosityNegative,config.viscosityPositive},initial,p.sigma,config.boundaryPatches));
        p.curves.open(OutputLayout::series()/"hydro_moment_interface.csv");if(!p.curves)throw std::runtime_error("cannot open moment interface history");
        p.curves<<std::setprecision(17)<<"step,time,parameter,x,y,component,curvature_model\n";
    } else {
        auto stress=Stress::buildCellwise(mesh,*p.domain,space,std::vector<double>(mesh.cellIndexNum(),config.viscosity));
        for(double v:p.domain->cellVolumes())p.mass+=config.density*v;
        p.flow=std::make_unique<Flow>(Flow::create(mesh,*p.domain,std::move(space),std::move(stress),config.density,initial,0,config.boundaryPatches));
    }
    if(!config.boundaryPatches.empty())std::cout<<"[Hydrodynamics] Initial boundary correction in moment mass norm="<<p.flow->initialBoundaryCorrection()<<"; physical time remains zero.\n";
    p.initialEnergy=p.flow->kineticEnergy();
    p.viscosity=config.viscosity;p.options.maximumIterations=config.momentMaximumIterations;p.options.relativeResidualTolerance=config.linearResidualTolerance;
    p.options.convectionFlux=config.momentConvection=="central"?Quantum::Hydrodynamics::HydroCompatibleMomentTransport::Flux::Central:Quantum::Hydrodynamics::HydroCompatibleMomentTransport::Flux::TraceUpwind;
    p.history.open(OutputLayout::series()/"hydro_moment_energy.csv");if(!p.history)throw std::runtime_error("cannot open moment energy history");
    p.history<<std::setprecision(17)<<"step,time,dt,pressure_time,pressure_available,mass,kinetic_energy,initial_energy,cumulative_transport_dissipation,cumulative_viscous_dissipation,cumulative_external_work,cumulative_pressure_work,cumulative_central_work,energy_balance,nonlinear_iterations,equation_residual,divergence_residual,trace_residual,local_momentum_defect,initial_surface_energy,surface_energy,cumulative_capillary_work,cumulative_temporal_mass_work,cumulative_material_metric_work,total_energy_balance,liquid_volume,liquid_volume_drift,fraction_error,shape_moment,cumulative_mesh_energy_change,leaf_cells,cumulative_boundary_traction_work,cumulative_boundary_reaction_work,cumulative_boundary_kinetic_flux,cumulative_inflow_dissipation,inflow_volume_rate,outflow_volume_rate,net_boundary_volume_rate,cumulative_event_kinetic_change,cumulative_event_surface_change,cumulative_event_dissipation,event_count,cumulative_geometric_mass_work,potential_energy,initial_potential_energy,mechanical_energy,mechanical_energy_change,cumulative_surface_advection,cumulative_surface_endpoint_work,cumulative_surface_flux_error,cumulative_surface_flux_budget,cumulative_surface_critical_estimate,surface_work_defect\n";
    if(p.twoPhase){
        std::ofstream meta(OutputLayout::meta()/"hydro_moment_case.json");
        meta<<std::setprecision(17)<<"{\n\"format_version\":1,\n\"coordinate\":"<<std::quoted(rz?"RZ":"XY")
            <<",\n\"domain\":["<<mesh.minCoord.x<<','<<mesh.minCoord.y<<','<<mesh.maxCoord.x<<','<<mesh.maxCoord.y
            <<"],\n\"density_negative\":"<<config.densityNegative<<",\n\"density_positive\":"<<config.densityPositive
            <<",\n\"viscosity_negative\":"<<config.viscosityNegative<<",\n\"viscosity_positive\":"<<config.viscosityPositive
            <<",\n\"surface_tension\":"<<p.sigma<<",\n\"acceleration\":["<<p.acceleration.x<<','<<p.acceleration.y
            <<"],\n\"initial_velocity\":"<<std::quoted(config.initialVelocity)
            <<",\n\"convection\":"<<std::quoted(config.momentConvection)
            <<",\n\"pressure_time_level\":\"step_midpoint\",\n\"geometry_time_level\":\"accepted_endpoint\",\n\"arc_basis\":\"power_in_fraction_0_1\",\n\"boundary_patches\":[";
        const std::array<const char*,4> sides{"west","south","east","north"};
        const std::array<const char*,3> types{"normalVelocity","pressureOutlet","freeSlip"};
        for(std::size_t i=0;i<config.boundaryPatches.size();i++){const auto& patch=config.boundaryPatches[i];
            if(i)meta<<',';meta<<"{\"side\":"<<std::quoted(sides.at(int(patch.side)))<<",\"type\":"<<std::quoted(types.at(int(patch.type)))
                <<",\"start\":"<<patch.start<<",\"end\":"<<patch.end<<",\"value\":"<<patch.value<<'}';}
        meta<<"]\n}\n";
        meta.flush();if(!meta)throw std::runtime_error("cannot write resolved moment case metadata");
    }
}
HydroMomentCase::~HydroMomentCase()=default;
const Flow::State& HydroMomentCase::state() const{return impl_->flow->state();}
InterfaceResolutionPlan::Result HydroMomentCase::Impl::review(const QuadtreeMesh& mesh,const InterfaceState& state,const char* stage){
    if(!adaptation)throw std::logic_error("resolution review requires declared adaptation controls");
    auto policy=adaptation->policy;policy.parametricNeckRadiusCells=neckRadiusCells;
    auto plan=InterfaceResolutionPlan::evaluate(mesh,*domain,state,policy,adaptation->cellBudget);
    if(!resolutionHistory.is_open()){
        resolutionHistory.open(OutputLayout::series()/"hydro_moment_resolution.csv");if(!resolutionHistory)throw std::runtime_error("cannot open moment resolution history");
        resolutionHistory<<std::setprecision(17)<<"review,step,time,stage,cells,required_leaves,planned_leaves,interface_max_level,saturated,terminal_radius,finest_event_scale,terminal_regions,terminal_geometric_volume,terminal_curvature_requests,terminal_gap_requests,max_terminal_curvature_cells\n";
    }
    const auto clock=state.progress();const int index=resolutionReview++;
    resolutionHistory<<index<<','<<clock.completedSteps<<','<<clock.endTime<<','<<stage<<','<<mesh.cellIndexNum()<<','<<plan.requiredLeafCount<<','<<plan.plannedLeafCount<<','<<plan.interfaceLevel<<','<<plan.saturated<<','<<plan.terminalRadius<<','<<plan.finestNeckScale<<','<<plan.terminalRanges.size()<<','<<plan.terminalGeometricVolume<<','<<plan.terminalCurvatureRequests<<','<<plan.terminalGapRequests<<','<<plan.maxTerminalCurvatureCells<<'\n';resolutionHistory.flush();
    if(!resolutionHistory)throw std::runtime_error("cannot append moment resolution history");
    if(!plan.terminalRanges.empty()){
        if(!terminalHistory.is_open()){terminalHistory.open(OutputLayout::series()/"hydro_moment_terminal_regions.csv");if(!terminalHistory)throw std::runtime_error("cannot open terminal-region history");
            terminalHistory<<std::setprecision(17)<<"review,step,time,stage,component,birth_parent,begin,end,lower_z,upper_z,pole,terminal_radius\n";}
        for(const auto& r:plan.terminalRanges)terminalHistory<<index<<','<<clock.completedSteps<<','<<clock.endTime<<','<<stage<<','<<r.component<<','<<r.birthParent<<','<<r.begin<<','<<r.end<<','<<r.lowerZ<<','<<r.upperZ<<','<<r.pole<<','<<plan.terminalRadius<<'\n';
        terminalHistory.flush();if(!terminalHistory)throw std::runtime_error("cannot append terminal-region history");
    }
    if(plan.saturated){
        // Preserve the exact accepted coefficients/inventory at the failed
        // zero-time review, independent of the normal field-write interval.
        const auto directory=OutputLayout::meta();const auto& geometry=state.parametricGeometry();
        const auto& components=geometry.boundaryHandle()->components();const auto& grid=geometry.grid();
        const auto& indices=InterfaceParametricMotion::cellIndices(state);
        std::ofstream snapshot(directory/"moment_resolution_rejected.dat");snapshot<<std::setprecision(17)
            <<"QUANTUM_MOMENT_RESOLUTION_DIAGNOSTIC_V1\n"<<clock.completedSteps<<' '<<clock.endTime<<' '<<domain->axisymmetric()<<' '<<components.size()<<' '<<grid.size()<<' '<<state.parametricNeckEvents().size()<<'\n'
            <<mesh.minCoord.x<<' '<<mesh.minCoord.y<<' '<<mesh.maxCoord.x<<' '<<mesh.maxCoord.y<<' '<<mesh.cellSize.x<<' '<<mesh.cellSize.y<<'\n'
            <<policy.maximumLevel<<' '<<policy.maxCurvatureCells<<' '<<policy.minimumOpposedDistanceCells<<' '<<policy.maximumGeometryFractionDifference<<' '<<policy.parametricNeckRadiusCells<<'\n';
        for(std::size_t i=0;i<components.size();i++){
            const auto& c=components[i];const auto& curve=*c.curve;const auto& b=curve.boundary;
            snapshot<<c.id<<' '<<curve.arcs.size()<<' '<<int(b.kind[0])<<' '<<int(b.kind[1])<<' '<<b.position[0].x<<' '<<b.position[0].y<<' '<<b.position[1].x<<' '<<b.position[1].y<<' '<<b.corners.size()<<'\n';
            for(const auto& p:b.corners)snapshot<<p.x<<' '<<p.y<<'\n';
            for(const auto& arc:curve.arcs){snapshot<<arc.a<<' '<<arc.b;for(const auto& p:arc.p)snapshot<<' '<<p.x<<' '<<p.y;snapshot<<'\n';}
            for(double v:geometry.componentInventory(i))snapshot<<v<<'\n';
        }
        for(std::size_t k=0;k<grid.size();k++){const int c=indices[k];const auto box=grid.rectangle(k);
            snapshot<<c<<' '<<mesh.cells(mesh.cellIndexToLabel[c]).level()<<' '<<box.left<<' '<<box.right<<' '<<box.bottom<<' '<<box.top<<' '<<geometry.inventory()[k]<<'\n';}
        std::ofstream witnesses(directory/"moment_resolution_unresolved.csv");witnesses<<std::setprecision(17)<<"kind,cell,component,other_component,begin,end,time,x,y,other_x,other_y,observed,allowed\n";
        for(const auto& w:plan.unresolved){using K=InterfaceResolutionPlan::Result::UnresolvedFeature::Kind;
            witnesses<<(w.kind==K::Curvature?"curvature":w.kind==K::OpposedGap?"opposed_gap":"inventory")<<','<<w.cell<<','<<w.component<<','<<w.otherComponent<<','<<w.begin<<','<<w.end<<','<<w.time<<','<<w.point.x<<','<<w.point.y<<','<<w.otherPoint.x<<','<<w.otherPoint.y<<','<<w.observed<<','<<w.allowed<<'\n';}
        snapshot.flush();witnesses.flush();if(!snapshot||!witnesses)throw std::runtime_error("cannot preserve rejected resolution observations");
    }
    return plan;
}
void HydroMomentCase::adapt(std::unique_ptr<QuadtreeMesh>& mesh,const char* stage){
    auto& p=*impl_;
    if(!p.adaptation)return;const auto& controls=*p.adaptation;
    if(!p.twoPhase)throw std::invalid_argument("moment AMR requires the paired parametric interface");
    const auto accept=[&](std::unique_ptr<QuadtreeMesh> target,std::unique_ptr<NativeLeafDomain> domain,
        Flow::MeshTransferCandidate candidate,std::size_t requiredLeaves,int interfaceLevel,const char* operation){
    const auto report=candidate.report();const int sourceCells=mesh->cellIndexNum();
    if(!p.adaptationHistory.is_open()){
        p.adaptationHistory.open(OutputLayout::series()/"hydro_moment_adaptation.csv");
        if(!p.adaptationHistory)throw std::runtime_error("cannot open moment mesh transfer history");
        p.adaptationHistory<<std::setprecision(17)<<"step,time,source_cells,target_cells,required_leaf_count,interface_level,inventory_balance,fraction_error,kinetic_before,kinetic_after,mesh_energy_change,minimum_kinetic,unconstrained_kinetic,retained_nullspace,momentum_residual,divergence_residual,operation,maximum_relative_velocity_error\n";
    }
    p.flow->acceptMeshTransfer(*mesh,*p.domain,*target,*domain,std::move(candidate));
    // noexcept exchanges retain both bound objects at their prepared addresses.
    mesh.swap(target);p.domain.swap(domain);
    const auto& level=report.levels.front();const double change=level.kineticAfter-level.kineticBefore;
    double maximumError=0;for(const auto& r:report.levels)maximumError=std::max(maximumError,r.relativeProjectionError);
    p.meshEnergyChange+=change;
    p.adaptationHistory<<p.flow->state().acceptedSteps<<','<<p.flow->state().time<<','<<sourceCells<<','<<mesh->cellIndexNum()<<','
        <<requiredLeaves<<','<<interfaceLevel<<','<<report.interface.inventoryBalance<<','<<report.interface.maximumFractionError<<','
        <<level.kineticBefore<<','<<level.kineticAfter<<','<<change<<','<<level.minimumKinetic<<','<<level.unconstrainedKinetic<<','<<level.retainedNullspace<<','
        <<level.momentumResidual<<','<<level.divergenceResidual<<','<<operation<<','<<maximumError<<'\n';p.adaptationHistory.flush();
    if(!p.adaptationHistory)throw std::runtime_error("cannot append accepted moment mesh transfer");
    std::cout<<std::setprecision(17)<<"[Hydrodynamics] Moment AMR accepted at step="<<p.flow->state().acceptedSteps<<" time="<<p.flow->state().time
        <<" operation="<<operation<<" cells="<<sourceCells<<"->"<<mesh->cellIndexNum()<<" mesh_energy_change="<<change<<" relative_velocity_error="<<maximumError<<" momentum_residual="<<level.momentumResidual<<'\n';
    };
    // Mandatory resolution is satisfied first. Optional coarsening is one
    // generation per review; a declined candidate cannot undo refinement.
    for(int sweep=0;;sweep++){
        const auto plan=p.review(*mesh,p.flow->interfaceState(),stage);p.resolution=plan;
        if(plan.saturated)throw std::length_error("moment AMR resolution request exceeds maximumLevel outside terminal scope");
        if(!plan.changesMesh())break;
        if(sweep>controls.policy.maximumLevel)throw std::runtime_error("local moment resolution did not reach a fixed mesh");
        auto target=InterfaceResolutionPlan::refine(*mesh,plan,controls.cellBudget);
        FaceSegmentGraph graph;graph.build(*target);
        auto domain=std::make_unique<NativeLeafDomain>(NativeLeafDomain::build(*target,graph,p.domain->axisymmetric()));
        auto candidate=p.flow->prepareMeshTransfer(*mesh,*p.domain,*target,*domain,p.options);
        accept(std::move(target),std::move(domain),std::move(candidate),plan.requiredLeafCount,plan.interfaceLevel,"refine");
    }
    if(!controls.coarsening)return;
    auto policy=controls.policy;policy.parametricNeckRadiusCells=p.neckRadiusCells;
    const auto plan=InterfaceResolutionPlan::planCoarsening(*mesh,*p.domain,p.flow->interfaceState(),policy,
        controls.coarsening->hysteresis,p.protectedBoundaryPoints);
    if(!plan.families)return;
    auto target=InterfaceResolutionPlan::coarsen(*mesh,plan);FaceSegmentGraph graph;graph.build(*target);
    auto domain=std::make_unique<NativeLeafDomain>(NativeLeafDomain::build(*target,graph,p.domain->axisymmetric()));
    auto assessment=p.flow->assessMeshTransfer(*mesh,*p.domain,*target,*domain,p.options,controls.coarsening->relativeVelocityError);
    if(!p.coarseningHistory.is_open()){
        p.coarseningHistory.open(OutputLayout::series()/"hydro_moment_coarsening.csv");
        if(!p.coarseningHistory)throw std::runtime_error("cannot open moment coarsening history");
        p.coarseningHistory<<std::setprecision(17)<<"step,time,source_cells,target_cells,families,status,level,level_time,relative_velocity_error,allowed_error,kinetic_before,kinetic_after,minimum_kinetic,momentum_residual,divergence_residual\n";
    }
    const int assessed=assessment.rejectedLevel>=0?assessment.rejectedLevel+1:int(assessment.report.levels.size());
    for(int i=0;i<assessed;i++){const auto& r=assessment.report.levels[i];
        p.coarseningHistory<<p.flow->state().acceptedSteps<<','<<p.flow->state().time<<','<<mesh->cellIndexNum()<<','<<target->cellIndexNum()<<','<<plan.families<<','
            <<transferStatus(assessment.readiness)<<','<<i<<','<<r.time<<','<<r.relativeProjectionError<<','<<controls.coarsening->relativeVelocityError<<','
            <<r.kineticBefore<<','<<r.kineticAfter<<','<<r.minimumKinetic<<','<<r.momentumResidual<<','<<r.divergenceResidual<<'\n';
    }
    p.coarseningHistory.flush();if(!p.coarseningHistory)throw std::runtime_error("cannot append moment coarsening history");
    if(!assessment.candidate)return;
    accept(std::move(target),std::move(domain),std::move(*assessment.candidate),plan.plannedLeafCount,p.resolution.interfaceLevel,"coarsen");
    p.resolution=p.review(*mesh,p.flow->interfaceState(),"coarsened");
    if(p.resolution.saturated||p.resolution.changesMesh())throw std::logic_error("coarsening violated the mandatory geometric retention band");
}
Flow::Report HydroMomentCase::advance(std::unique_ptr<QuadtreeMesh>& meshOwner,int step,double begin,double dt){
    auto& p=*impl_;
    if(step>1)adapt(meshOwner,"physical");
    if(p.neckRadiusCells>0){
        int events=0;
        for(;;){const auto sites=InterfaceParametricProjection::inspectNecks(p.flow->interfaceState().parametricGeometry(),p.neckRadiusCells);bool accepted=false;
          for(const auto& site:sites){
            const bool admitted=!p.adaptation||InterfaceResolutionPlan::admitsNeck(*meshOwner,p.flow->interfaceState(),p.resolution,p.adaptation->policy.maximumLevel,site);
            Flow::NeckAssessment assessment;const bool inspectGeometry=admitted&&site.radius<=p.neckRadiusCells*site.localScale;
            if(inspectGeometry)assessment=p.flow->assessNeckSplit(*meshOwner,*p.domain,site.options);
            else {auto& g=assessment.report.geometry;g.parent=site.options.parent;g.section=site.options.section;g.neckRadius=site.radius;g.neckPosition=site.position;g.localScale=site.localScale;}
            if(!p.neckHistory.is_open()){
                p.neckHistory.open(OutputLayout::series()/"hydro_moment_necks.csv");if(!p.neckHistory)throw std::runtime_error("cannot open neck assessment history");
                p.neckHistory<<std::setprecision(17)<<"step,time,parent,section,neck_radius,neck_position,local_scale,status,geometry_assessed,caps_assessed,fluid_assessed,section_radius,lower_pole,upper_pole,model_volume,minimum_kinetic,kinetic_ceiling,constraint_residual,lower_tangent_r,upper_tangent_r,lower_curvature,upper_curvature,lower_cap_length,upper_cap_length,lower_cap_discriminant,upper_cap_discriminant\n";
            }
            const auto& q=assessment.report;const auto& observed=q.geometry;
            const bool caps=assessment.readiness!=Flow::NeckReadiness::ScalePending;
            const bool fluid=caps&&assessment.readiness!=Flow::NeckReadiness::CapPending;
            p.neckHistory<<p.flow->state().acceptedSteps<<','<<p.flow->state().time<<','<<observed.parent<<','<<observed.section<<','<<observed.neckRadius<<','<<observed.neckPosition<<','<<observed.localScale<<','<<neckStatus(assessment.readiness)<<','<<inspectGeometry<<','<<caps<<','<<fluid<<','<<observed.sectionRadius<<','<<observed.lowerPole<<','<<observed.upperPole<<','<<observed.modelRegionVolume<<','<<q.momentum.minimumKinetic<<','<<q.momentum.kineticCeiling<<','<<q.momentum.constraintResidual<<','<<observed.seamRadialTangent[0]<<','<<observed.seamRadialTangent[1]<<','<<observed.seamSignedCurvature[0]<<','<<observed.seamSignedCurvature[1]<<','<<observed.capLength[0]<<','<<observed.capLength[1]<<','<<observed.capDiscriminant[0]<<','<<observed.capDiscriminant[1]<<'\n';p.neckHistory.flush();
            if(!p.neckHistory)throw std::runtime_error("cannot append neck assessment history");
            if(!assessment.candidate)continue;
            if(++events>meshOwner->cellIndexNum())throw std::runtime_error("neck event sequence exceeds finite native support");
            auto candidate=std::move(*assessment.candidate);const auto r=candidate.report();
            if(!p.eventHistory.is_open()){
                p.eventHistory.open(OutputLayout::series()/"hydro_moment_events.csv");if(!p.eventHistory)throw std::runtime_error("cannot open neck event history");
                p.eventHistory<<std::setprecision(17)<<"step,time,parent,lower_child,upper_child,section,neck_radius,local_scale,lower_pole,upper_pole,model_volume,parent_volume,lower_volume,upper_volume,cell_balance,fraction_error,kinetic_before,kinetic_after,surface_before,surface_after,model_dissipation,momentum_residual,divergence_residual,mass_matrix_change\n";
            }
            p.flow->acceptNeckSplit(*meshOwner,*p.domain,std::move(candidate));
            p.eventKineticChange+=r.kineticChange;p.eventSurfaceChange+=r.surfaceAfter-r.surfaceBefore;p.eventDissipation+=r.modelDissipation;
            const auto& g=r.geometry;for(int child=0;child<2;child++)p.initialComponentVolumes.emplace(g.children[child],g.childVolumes[child]);
            p.eventHistory<<p.flow->state().acceptedSteps<<','<<p.flow->state().time<<','<<g.parent<<','<<g.children[0]<<','<<g.children[1]<<','<<g.section<<','<<g.neckRadius<<','<<g.localScale<<','<<g.lowerPole<<','<<g.upperPole<<','<<g.modelRegionVolume<<','<<g.parentVolume<<','<<g.childVolumes[0]<<','<<g.childVolumes[1]<<','<<g.maximumCellBalance<<','<<g.maximumInventoryFractionError<<','<<r.momentum.kineticBefore<<','<<r.momentum.kineticAfter<<','<<r.surfaceBefore<<','<<r.surfaceAfter<<','<<r.modelDissipation<<','<<r.momentum.momentumResidual<<','<<r.momentum.divergenceResidual<<','<<r.relativeMassMatrixChange<<'\n';p.eventHistory.flush();
            if(!p.eventHistory)throw std::runtime_error("cannot append accepted neck event");
            std::cout<<std::setprecision(17)<<"[Hydrodynamics] Moment neck event at time="<<p.flow->state().time<<" parent="<<g.parent<<" children="<<g.children[0]<<','<<g.children[1]<<" model_dissipation="<<r.modelDissipation<<" momentum_residual="<<r.momentum.momentumResidual<<'\n';
            adapt(meshOwner,"event");accepted=true;break;
          }
          if(!accepted)break;
        }
    }
    const auto& mesh=*meshOwner;
    const auto load=p.flow->makeUniformAccelerationLoad(mesh,*p.domain,p.acceleration,begin+.5*dt);
    auto candidate=p.flow->prepareNavierStokesStep(mesh,*p.domain,{step,begin,dt,HydroTimePolicy::Fixed},load,p.options);
    const auto report=candidate.report();p.flow->accept(mesh,*p.domain,std::move(candidate));return report;
}
void HydroMomentCase::write(const QuadtreeMesh& mesh,const std::string& prefix,bool fields,const Flow::Report* report){
    auto& p=*impl_;const auto& state=p.flow->state();const auto& space=p.flow->space();space.validate(mesh,*p.domain);
    if(state.acceptedSteps!=p.writtenStep+1 || (state.acceptedSteps>0)!=bool(report))throw std::logic_error("moment output must follow each accepted state exactly once");
    const double energy=p.flow->kineticEnergy();
    if(report){p.transport+=report->transportDissipation;p.viscous+=report->viscousDissipation;p.external+=report->externalWork;p.pressure+=report->pressureWork;p.central+=report->centralConvectionWork;
        p.capillary+=report->capillaryWork;p.temporal+=report->temporalMassWork;p.metric+=report->materialMetricWork;p.geometricWork+=report->geometricMassWork;
        p.surfaceAdvection+=report->surfaceAdvection;p.surfaceEndpoint+=report->surfaceEndpointWork;
        p.surfaceFluxError+=report->surfaceFluxError;p.surfaceFluxBudget+=report->surfaceFluxBudget;p.surfaceCriticalError+=report->surfaceCriticalError;
        p.boundaryTraction+=report->boundaryTractionWork;p.boundaryReaction+=report->boundaryReactionWork;
        p.boundaryKinetic+=report->boundaryKineticFlux;p.boundaryDissipation+=report->boundaryInflowDissipation;}
    const double undefined=std::numeric_limits<double>::quiet_NaN();
    double liquid=0,fraction=0,surface=0,shape=0,potential=undefined;std::vector<double> alpha(mesh.cellIndexNum(),0);
    if(p.twoPhase){const auto& interface=p.flow->interfaceState();const auto& geometry=interface.parametricGeometry();const auto& cells=InterfaceParametricMotion::cellIndices(interface);
        const auto material=p.flow->material().properties();p.mass=0;double positiveInventory=0;
        for(int slot=0;slot<int(cells.size());slot++){const int c=cells[slot];const double volume=p.domain->cellVolumes()[p.domain->cellRowByIndex()[c]],v=geometry.inventory()[slot];
            liquid+=v;alpha[c]=v/volume;fraction=std::max(fraction,std::abs(v-p.flow->material().cells()[c].liquid(0,0))/volume);
            positiveInventory+=volume-v;p.mass+=material.liquidDensity*v+material.gasDensity*(volume-v);}
        const bool rz=p.domain->axisymmetric();const auto surfaceReport=Quantum::Hydrodynamics::HydroCompatibleMomentTraction::measure(*geometry.boundaryHandle(),mesh,*p.domain);const double surfaceMeasure=surfaceReport.length;
        surface=p.sigma*surfaceMeasure;const auto integrals=p.flow->phaseIntegrals(mesh,*p.domain);
        const std::array<double,2> inventory{liquid,positiveInventory};
        if(state.acceptedSteps==0)p.initialPhaseInventory=inventory;
        bool closedBody=!surfaceReport.clipped;
        for(std::size_t i=0;i<geometry.componentCount();i++){
            const auto& component=geometry.boundaryHandle()->components()[i];const auto& boundary=component.curve->boundary;
            closedBody=closedBody&&(rz?(!boundary.fixed(0)&&!boundary.fixed(1)):!boundary.open());
            double geometric=0;for(const auto& cell:p.flow->material().cells())geometric+=cell.componentLiquid.at(i)(0,0);
            double v=0;for(double x:geometry.componentInventory(i))v+=x;
            p.componentHistory<<state.acceptedSteps<<','<<state.time<<','<<component.id<<','<<v<<','<<v-p.initialComponentVolumes.at(component.id)<<','<<geometric<<','<<component.curve->arcs.size()<<",active,nan\n";
        }
        for(const auto& component:geometry.retiredComponents()){
            double remainder=0;for(double v:*component.inventory)remainder+=v;
            p.componentHistory<<state.acceptedSteps<<','<<state.time<<','<<component.id<<','<<remainder<<','<<remainder-p.initialComponentVolumes.at(component.id)<<",0,0,retired,"<<component.time<<'\n';
        }
        p.componentHistory.flush();if(!p.componentHistory)throw std::runtime_error("cannot append component inventory");
        const auto& negative=integrals[0];shape=undefined;
        if(negative.present){const auto& m=negative.centralSecondMoment;
            const double transverse=rz?.5*(m[0]+negative.volume*negative.coordinateMean.x*negative.coordinateMean.x):m[0],axial=m[2];
            if(transverse+axial>0)shape=(transverse-axial)/(transverse+axial);
        }
        if(!p.phaseHistory.is_open()){
            p.phaseHistory.open(OutputLayout::series()/"hydro_moment_phase_integrals.csv");
            if(!p.phaseHistory)throw std::runtime_error("cannot open moment phase integral history");
            p.phaseHistory<<std::setprecision(17)<<"step,time,phase,present,volume_fv,volume_geometry,geometry_minus_fv,volume_fv_drift,"
                <<(rz?"mean_radius,centroid_z,integral_u_r,integral_u_z,mean_u_r,mean_u_z,central_rr,central_rz,central_zz,":"centroid_x,centroid_y,integral_u_x,integral_u_y,mean_u_x,mean_u_y,central_xx,central_xy,central_yy,")
                <<"interface_measure,circularity,sphericity,mean_curvature,curvature_squared_integral,surface_quadrature_relative_estimate,surface_clipped\n";
        }
        const double pi=std::acos(-1.);potential=0;
        for(int phase=0;phase<2;phase++){const auto& r=integrals[phase];
            if(r.present)potential-=(phase==0?material.liquidDensity:material.gasDensity)*r.volume*(p.acceleration.x*r.coordinateMean.x+p.acceleration.y*r.coordinateMean.y);
            const bool closed=phase==0&&closedBody&&r.present&&surfaceMeasure>0;
            const double circularity=closed&&!rz?2*std::sqrt(pi*r.volume)/surfaceMeasure:undefined;
            const double sphericity=closed&&rz?std::cbrt(pi)*std::pow(6*r.volume,2./3)/surfaceMeasure:undefined;
            p.phaseHistory<<state.acceptedSteps<<','<<r.time<<','<<(phase==0?"negative":"positive")<<','<<r.present<<','<<inventory[phase]<<','<<r.volume<<','<<r.volume-inventory[phase]<<','<<inventory[phase]-p.initialPhaseInventory[phase]<<','
                <<r.coordinateMean.x<<','<<r.coordinateMean.y<<','<<r.velocityIntegral.x<<','<<r.velocityIntegral.y<<','<<r.meanVelocity.x<<','<<r.meanVelocity.y<<','
                <<r.centralSecondMoment[0]<<','<<r.centralSecondMoment[1]<<','<<r.centralSecondMoment[2]<<','<<surfaceMeasure<<','<<circularity<<','<<sphericity<<','
                <<surfaceReport.meanCurvature<<','<<surfaceReport.curvatureSquaredIntegral<<','<<surfaceReport.quadratureError<<','<<surfaceReport.clipped<<'\n';
        }
        if(state.acceptedSteps==0)p.initialPotential=potential;
        p.phaseHistory.flush();if(!p.phaseHistory)throw std::runtime_error("cannot append moment phase integrals");
    }
    const double mechanical=energy+surface+potential,mechanicalChange=mechanical-p.initialEnergy-p.initialSurface-p.initialPotential;
    const double balance=energy-p.initialEnergy-p.meshEnergyChange-p.eventKineticChange+p.transport+p.viscous-p.external-p.pressure
        -p.boundaryTraction-p.boundaryReaction+(p.twoPhase?p.metric-p.temporal-p.capillary:p.central);
    p.history<<state.acceptedSteps<<','<<state.time<<','<<state.lastDt<<','<<state.stageTime.value_or(0)<<','<<bool(state.stageTime)<<','<<p.mass<<','<<energy<<','<<p.initialEnergy<<','
        <<p.transport<<','<<p.viscous<<','<<p.external<<','<<p.pressure<<','<<p.central<<','<<balance<<','
        <<(report?report->nonlinearIterations:0)<<','<<(report?report->equationResidual:0)<<','<<(report?report->divergenceResidual:0)<<','<<(report?report->viscousTraceResidual:0)<<','<<(report?report->localMomentumDefect:0)<<','
        <<p.initialSurface<<','<<surface<<','<<p.capillary<<','<<p.temporal<<','<<p.metric<<','<<balance+surface-p.initialSurface+p.capillary-p.eventSurfaceChange+p.surfaceAdvection-p.surfaceEndpoint<<','<<liquid<<','<<(p.twoPhase?liquid-p.initialLiquid:0)<<','<<fraction<<','<<shape<<','<<p.meshEnergyChange<<','<<mesh.cellIndexNum()<<','<<p.boundaryTraction<<','<<p.boundaryReaction<<','<<p.boundaryKinetic<<','<<p.boundaryDissipation<<','<<(report?report->boundaryInflowRate:0)<<','<<(report?report->boundaryOutflowRate:0)<<','<<(report?report->boundaryVolumeRate:0)<<','<<p.eventKineticChange<<','<<p.eventSurfaceChange<<','<<p.eventDissipation<<','<<(p.twoPhase?p.flow->interfaceState().parametricNeckEvents().size():0)<<','<<p.geometricWork<<','<<potential<<','<<p.initialPotential<<','<<mechanical<<','<<mechanicalChange<<','<<p.surfaceAdvection<<','<<p.surfaceEndpoint<<','<<p.surfaceFluxError<<','<<p.surfaceFluxBudget<<','<<p.surfaceCriticalError<<','<<(report?report->surfaceWorkDefect:0)<<'\n';
    p.history.flush();if(!p.history)throw std::runtime_error("cannot append moment energy history");
    if(fields){
        std::vector<std::array<double,2>> sampledSpeed(mesh.cellIndexNum());
        if(state.physicalVelocity)for(const auto& row:state.physicalVelocity->cellSpeedSamples())sampledSpeed.at(row.cell)={row.coarseMaximum,row.fineMaximum};
        TecplotWriter::Zone zone;zone.title="Compatible moment flow";zone.solutionTime=state.time;
        zone.auxiliary={{"TimeScheme","cchMomentMidpoint"},{"VelocityTime",number(state.time)},{"PressureTime",number(state.stageTime.value_or(0))},
            {"PressureTimeLevel",state.stageTime?"step_midpoint":"unavailable"},{"PressureAvailable",state.stageTime?"1":"0"},
            {"FieldRepresentation",state.physicalVelocity?"complete native velocity and physical first jets; P1 pressure":"cell center P1 velocity/pressure with physical gradients; native leaves"},
            {"Coordinates",p.domain->axisymmetric()?"RZ: x=r, y=z":"XY"},{"Units","code"},
            {"VelocitySampling",state.physicalVelocity?"nested 5x5 and 9x9 complete-field samples per leaf; not continuous bounds":"affine field vertex maxima"}};
        for(const char* name:{"u_x","u_y","speed","p","du_x_dx","du_x_dy","du_y_dx","du_y_dy","dp_dx","dp_dy","rho","mu","liquid_fraction_FV","liquid_fraction_geometry","speed_sampled_max_5","speed_sampled_max_9","p_volume_mean"})zone.fields.push_back({name,TecplotWriter::Location::CellCentered,{}});
        for(int c=0;c<mesh.cellIndexNum();c++){
            const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto center=cell.centerCoord(),h=cell.size();const int offset=int(zone.points.size());
            zone.points.insert(zone.points.end(),{{center.x-.5*h.x,center.y-.5*h.y},{center.x+.5*h.x,center.y-.5*h.y},{center.x+.5*h.x,center.y+.5*h.y},{center.x-.5*h.x,center.y+.5*h.y}});
            zone.connectivity.push_back({offset,offset+1,offset+2,offset+3});const auto& b=space.cells()[c];Eigen::VectorXd local(b.rows.size());for(int j=0;j<local.size();j++)local[j]=state.velocity[b.rows[j]];
            const Eigen::VectorXd coefficient=b.affineProjection*local;Eigen::Vector2d u=Space::affineVelocityBasis({0,0},h)*coefficient;
            Eigen::Vector2d dx=(Space::affineVelocityBasis({.5,0},h)-Space::affineVelocityBasis({-.5,0},h))*coefficient/h.x;
            Eigen::Vector2d dy=(Space::affineVelocityBasis({0,.5},h)-Space::affineVelocityBasis({0,-.5},h))*coefficient/h.y;
            if(state.physicalVelocity){const auto jet=state.physicalVelocity->firstJet(center);u={jet.value.x,jet.value.y};dx={jet.dx.x,jet.dx.y};dy={jet.dy.x,jet.dy.y};}
            else for(double x:{-.5,.5})for(double y:{-.5,.5}){
                const double speed=(Space::affineVelocityBasis({x,y},h)*coefficient).norm();
                sampledSpeed[c][0]=sampledSpeed[c][1]=std::max(sampledSpeed[c][0],speed);
            }
            const Eigen::Vector3d pressure=state.stageTime?Eigen::Vector3d(state.pressure.segment<3>(3*c)):Eigen::Vector3d::Zero();
            double density=0,viscosity=p.viscosity,geometric=0;
            if(p.twoPhase){const auto& material=p.flow->material().cells()[c];const double volume=material.volume(0,0);density=material.density(0,0)/volume;viscosity=material.viscosity(0,0)/volume;geometric=material.liquid(0,0)/volume;}
            else density=p.flow->density();
            // RZ P1 pressure is centered geometrically, not in physical volume.
            // Its exact full-cell mean includes integral(X*2*pi*r)/volume.
            const double meanPressure=pressure[0]+(p.domain->axisymmetric()?pressure[1]*h.x/(12*center.x):0);
            const double data[]={u[0],u[1],u.norm(),pressure[0],dx[0],dy[0],dx[1],dy[1],pressure[1]/h.x,pressure[2]/h.y,density,viscosity,alpha[c],geometric,sampledSpeed[c][0],sampledSpeed[c][1],meanPressure};
            for(int j=0;j<17;j++)zone.fields[j].values.push_back(data[j]);
        }
        TecplotWriter::write(OutputLayout::fields()/(stem(prefix,state.acceptedSteps)+".dat"),zone);
        TecplotWriter::Zone faces;faces.title="Owned face normal velocity moments";faces.solutionTime=state.time;faces.element=TecplotWriter::Element::Line;faces.auxiliary={{"TimeScheme","cchMomentMidpoint"},{"NormalOrientation","+x on vertical faces, +y on horizontal faces"},{"Units","code"}};
        for(const char* name:{"normal_mean","normal_linear","tangent_center","tangent_scale","area","volume_rate"})faces.fields.push_back({name,TecplotWriter::Location::CellCentered,{}});
        for(const auto face:p.domain->faces(mesh)){
            const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);const int offset=int(faces.points.size());
            faces.points.push_back(mesh.verts(edge.edgeToVert[0]).coord());faces.points.push_back(mesh.verts(edge.edgeToVert[1]).coord());faces.connectivity.push_back({offset,offset+1});
            const auto trace=p.flow->normalPolynomial(mesh,*p.domain,face.nativeRow);const double data[]={trace.mean,trace.linear,trace.tangentCenter,trace.tangentScale,trace.area,trace.volumeRate()};
            for(int j=0;j<6;j++)faces.fields[j].values.push_back(data[j]);
        }
        TecplotWriter::write(OutputLayout::fields()/(stem(prefix+"_faces",state.acceptedSteps)+".dat"),faces);
        if(p.twoPhase){
            const auto& components=p.flow->interfaceState().parametricGeometry().boundaryHandle()->components();TecplotWriter::Zone boundary;
            std::ofstream arcs(OutputLayout::fields()/(stem(prefix+"_arcs",state.acceptedSteps)+".csv"));
            arcs<<std::setprecision(17)<<"step,time,component,axisymmetric,begin_kind,end_kind,a,b,p0x,p0y,p1x,p1y,p2x,p2y,p3x,p3y,p4x,p4y\n";
            const auto kind=[](Quantum::Geometry::Parametric::CurveBoundary::Kind k){using K=Quantum::Geometry::Parametric::CurveBoundary::Kind;return k==K::Axis?"axis":k==K::Fixed?"fixed":"periodic";};
            for(const auto& component:components)for(const auto& arc:component.curve->arcs){
                arcs<<state.acceptedSteps<<','<<state.time<<','<<component.id<<','<<p.domain->axisymmetric()<<','<<kind(component.curve->boundary.kind[0])<<','<<kind(component.curve->boundary.kind[1])<<','<<arc.a<<','<<arc.b;
                for(const auto& point:arc.p)arcs<<','<<point.x<<','<<point.y;arcs<<'\n';
            }
            arcs.flush();if(!arcs)throw std::runtime_error("cannot write accepted Hermite polynomial frame");
            boundary.title="Accepted Hermite liquid interface";boundary.element=TecplotWriter::Element::Line;boundary.solutionTime=state.time;
            boundary.auxiliary={{"GeometrySamples","polynomial critical points plus eight uniform intervals per arc"},
                {"CurvatureMeaning",p.domain->axisymmetric()?"one-sided meridional plus azimuthal curvature":"one-sided planar curvature"}};
            boundary.fields={{"curvature_one_sided",TecplotWriter::Location::Nodal,{}},{"component",TecplotWriter::Location::Nodal,{}},{"curvature_model",TecplotWriter::Location::Nodal,{}},
                {"curvature_meridional",TecplotWriter::Location::Nodal,{}},{"curvature_azimuthal",TecplotWriter::Location::Nodal,{}}};
            for(const auto& component:components){const auto& geometry=*component.curve;
            for(const auto& arc:geometry.arcs){const int start=int(boundary.points.size());
                auto samples=Quantum::Geometry::Parametric::surfacePartition(arc,geometry.boundary,p.domain->axisymmetric());for(int k=0;k<=8;k++)samples.push_back(k/8.);
                std::sort(samples.begin(),samples.end());samples.erase(std::unique(samples.begin(),samples.end()),samples.end());
                for(int k=0;k<int(samples.size());k++){const double t=samples[k];
                const auto meridian=Quantum::Geometry::Parametric::sampleSurface(arc,t);
                const auto q=p.domain->axisymmetric()?Quantum::Geometry::Parametric::sampleSurfaceRZ(arc,t,geometry.boundary):meridian;
                boundary.points.push_back({q.point.x,q.point.y});boundary.fields[0].values.push_back(q.curvature);boundary.fields[1].values.push_back(double(component.id));
                boundary.fields[2].values.push_back(p.flow->interfaceState().parametricModelCurvature(component.id,arc.a+(arc.b-arc.a)*t)?1:0);
                boundary.fields[3].values.push_back(meridian.curvature);boundary.fields[4].values.push_back(q.curvature-meridian.curvature);
                if(k)boundary.connectivity.push_back({start+k-1,start+k});}}
            for(int k=0;k<=512;k++){const double t=double(k)/512;const auto x=geometry.sample(t);p.curves<<state.acceptedSteps<<','<<state.time<<','<<t<<','<<x.x<<','<<x.y<<','<<component.id<<','<<p.flow->interfaceState().parametricModelCurvature(component.id,t)<<'\n';}}
            // An absent surface has no line zone. The component/phase history
            // records that state explicitly; never output a dummy droplet.
            if(!components.empty())TecplotWriter::write(OutputLayout::fields()/(stem(prefix+"_interface",state.acceptedSteps)+".dat"),boundary);
            p.curves.flush();if(!p.curves)throw std::runtime_error("cannot append accepted moment interface");
        }
    }
    p.writtenStep=state.acceptedSteps;
}
