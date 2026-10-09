#include "Hydrodynamics/HydroConservativeInterface.h"
#include "Modules/Hydrodynamics/Hydrodynamics.h"
#include "Modules/Hydrodynamics/HydrodynamicsConfig.h"
#include "Modules/Hydrodynamics/HydroFittedVelocity.h"
#include "Modules/Hydrodynamics/HydroInterfaceTransport.h"
#include "Modules/Hydrodynamics/HydroUniformLoad.h"
#include "Modules/Hydrodynamics/IO/HydroGeometryLoader.h"
#include "Modules/Interface/Compatibility/InterfaceNodalConfig.h"
#include "Modules/Interface/Compatibility/InterfaceConservativeGeometry.h"
#include "Modules/Interface/Compatibility/InterfaceTopology.h"
#include "Modules/Interface/Compatibility/InterfaceMotion.h"
#include "Modules/Interface/Compatibility/InterfaceNeckCap.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/UniformLeafGrid.h"
#include "Core/IO/OutputLayout.h"
#include "Core/IO/TecplotWriter.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>

HydroConservativeInterface::HydroConservativeInterface(QuadtreeMesh& mesh,const VertField<double>* phi,
    const HydrodynamicsConfig& flow,const InterfaceNodalConfig& input,const HydroGeometryInput::AxisymmetricColumn* column,
    bool retainHistory):mesh_(mesh),flow_(flow),retainHistory_(retainHistory),neckBreakup_(input.neckBreakup)
{
    const bool rz=flow.coordinateSystem=="axisymmetricRZ";
    if(neckBreakup_ && (!column || !rz || !flow.boundaryPatches.empty()))
        throw std::invalid_argument("neckBreakup requires the conservative RZ reflection-column consumer");
    if(neckBreakup_ && (neckBreakup_->maximumLevel<0 || neckBreakup_->maximumLevel>24
        || !(neckBreakup_->radiusCells>0 && neckBreakup_->radiusCells<1)))
        throw std::invalid_argument("invalid terminal neck model controls");
    FaceSegmentGraph graph;graph.build(mesh);domain_=NativeLeafDomain::build(mesh,graph,rz);
    const UniformLeafGrid grid(mesh);
    Quantum::Geometry::ActiveHermitePatch patch;patch.activeCells.assign(mesh.cellIndexNum(),1);
    patch.nodes.resize(mesh.vertIndexNum());
    std::vector<LevelSetInterfacePanels::SymmetryPlane> sections;
    if(column) {
        if(phi || !rz || !flow.boundaryPatches.empty())throw std::invalid_argument("column initializer requires only the RZ reflection domain");
        columnHalfWaves_=column->halfWaves;
        const double k=columnHalfWaves_*std::acos(-1.)/(mesh.maxCoord.y-mesh.minCoord.y);
        columnRadius_=column->radius;
        const double scale=column->radius/std::sqrt(1+column->relativeAmplitude*column->relativeAmplitude/2);
        if(neckBreakup_ && !(scale*(1-column->relativeAmplitude)>
            neckBreakup_->radiusCells*(mesh.maxCoord.x-mesh.minCoord.x)/mesh.Nx/std::ldexp(1.,neckBreakup_->maximumLevel)))
            throw std::invalid_argument("initial column neck must exceed the terminal cutoff");
        sections={{1,mesh.minCoord.y},{1,mesh.maxCoord.y}};
        for(int v=0;v<mesh.vertIndexNum();++v) {
            const auto p=mesh.verts(mesh.vertIndexToLabel[v]).coord();
            const double theta=k*(p.y-mesh.minCoord.y);
            const double radius=scale*(1+column->relativeAmplitude*std::cos(theta));
            const double derivative=(p.y==mesh.minCoord.y || p.y==mesh.maxCoord.y)?0:
                -scale*column->relativeAmplitude*k*std::sin(theta);
            patch.nodes[v]={p.x*p.x-radius*radius,2*p.x,-2*radius*derivative,0};
        }
    } else {
    if(!phi)throw std::invalid_argument("conservative initialization requires a scalar or analytic column");
    const LevelSetHermiteGeometry::Snapshot initial(mesh,*phi,rz && mesh.minCoord.x==0);
    // One-time import of the existing geometry initializer's cubic field.
    // Nodal legacy geometry is discarded by the caller after this conversion.
    for(int j=0;j<=grid.ny();++j) for(int i=0;i<=grid.nx();++i) {
        const auto jet=initial.scalarJetAtCellLocal(std::min(i,grid.nx()-1),std::min(j,grid.ny()-1),
            {i==grid.nx()?1.:0.,j==grid.ny()?1.:0.});
        patch.nodes[mesh.verts(grid.vertices()[i+(grid.nx()+1)*j]).index()]={jet.value,jet.dx,jet.dy,jet.dxy};
    }
    }
    state_=std::make_unique<InterfaceState>(mesh);
    InterfaceGeometryReconstruction::initialize(mesh,std::move(patch),rz,{},*state_,sections);
    const auto transport=HydroInterfaceTransport::domainView(mesh,domain_,flow.boundaryPatches);
    // Zero-time inventory seeding uses geometric measurement only. This zero
    // diagnostic flux is never used as a physical NS/advection time layer.
    const VEdgeField<double> zeroX(mesh,0);const HEdgeField<double> zeroY(mesh,0);
    const HydroFittedVelocitySampler zero(mesh,domain_,zeroX,zeroY,
        std::vector<HydroFittedVelocitySampler::VertexTrace>(mesh.vertIndexNum()));
    state_->initializePhase(mesh,transport,HydroInterfaceTransport::normalFlux(zero));
    state_->initializeConservativePhase(mesh,domain_);
    if(neckBreakup_) {
        // Use the ordinary synchronized scalar gauge before introducing a
        // terminal region. Splicing an arbitrary analytic r^2 amplitude into
        // a distance-gauge reconstruction would distort the seam next step.
        const InterfaceConservativeGeometry::BoundaryInput boundary{transport,{}, {}};
        auto synchronized=std::make_unique<InterfaceState>(mesh);
        InterfaceConservativeGeometry::resynchronize(mesh,domain_,*state_,*synchronized,
            std::max(4.,flow_.interfaceThickness/std::min(grid.spacing().x,grid.spacing().y)+2),&boundary,
            InterfaceConservativeGeometry::CurvatureDiagnostics::Omit);
        state_=InterfaceNeckCap::prepareTerminalCurvature(mesh_,domain_,std::move(synchronized),
            mesh_.maxCoord.y,neckBreakup_->maximumLevel,true);
    }
    phases_.assign(domain_.nativeFaceCount(),0);
    for(const auto face:domain_.faces(mesh)) if(transport.boundary(face.nativeRow)==InterfaceTransportDomain::Boundary::SinglePhasePort) {
        const auto side=face.horizontal?(face.outwardSign<0?HydroBoundarySide::South:HydroBoundarySide::North)
            :(face.outwardSign<0?HydroBoundarySide::West:HydroBoundarySide::East);
        const double position=face.horizontal?face.center.x:face.center.y;
        for(const auto& p:flow.boundaryPatches) if(p.side==side && position>p.start && position<p.end)
            for(const auto& phase:input.portPhases) if(phase.patch==p.name) phases_[face.nativeRow]=phase.sign;
        const int sign=phases_[face.nativeRow],c=face.negativeCell>=0?face.negativeCell:face.positiveCell;
        const auto& inventory=state_->conservativePhase();
        if(sign==0 || state_->topology().cellPhase()[c]!=sign || inventory.liquidVolumes()[c]!=(sign<0?inventory.capacities()[c]:0))
            throw std::invalid_argument("conservative port phase disagrees with the initial pure geometric/inventory layer");
    }
}
HydroConservativeInterface::MeshCandidate::MeshCandidate()=default;
HydroConservativeInterface::MeshCandidate::~MeshCandidate()=default;
HydroConservativeInterface::MeshCandidate::MeshCandidate(MeshCandidate&&) noexcept=default;

HydroConservativeInterface::HydroConservativeInterface(QuadtreeMesh& mesh,const HydrodynamicsConfig& flow,
    NativeLeafDomain domain,std::unique_ptr<InterfaceState> current,std::unique_ptr<InterfaceState> previous,double radius,
    std::optional<InterfaceNodalConfig::NeckBreakup> breakup,int halfWaves)
    :mesh_(mesh),flow_(flow),domain_(std::move(domain)),state_(std::move(current)),previous_(std::move(previous)),
    retainHistory_(true),phases_(domain_.nativeFaceCount(),0),columnRadius_(radius),columnHalfWaves_(halfWaves),neckBreakup_(breakup)
{}

HydroConservativeInterface::MeshCandidate HydroConservativeInterface::prepareAdaptation(
    Hydrodynamics& hydro,const InterfaceAdaptationControls& controls) const
{
    MeshCandidate candidate;
    // Both physical histories must exist; initialization uses the declared
    // root grid, which must support the first physical step on its own.
    if(!previous_)return candidate;
    if(!flow_.boundaryPatches.empty())throw std::invalid_argument("conservative Hydro AMR requires a closed domain");
    auto policy=controls.policy;
    if(neckBreakup_)policy.terminalAxialSection=mesh_.maxCoord.y;
    const auto spacing=rootSpacing();
    policy.distanceBandCells=std::max(4.,flow_.interfaceThickness/std::min(spacing.x,spacing.y)+2);
    // Review the actual accepted geometries each step. The fitted sampler has
    // no certified speed bound: nodal max-speed*dt is NOT a travel certificate.
    policy.maximumTravel=0;
    candidate.plan=InterfaceResolutionPlan::evaluate(mesh_,domain_,*state_,*previous_,policy,controls.cellBudget);
    if(candidate.plan.saturated)throw std::length_error("Hydro AMR interface resolution request exceeds maximumLevel");
    if(!candidate.plan.changesMesh())return candidate;
    candidate.mesh=InterfaceResolutionPlan::refine(mesh_,candidate.plan,controls.cellBudget);
    auto& target=*candidate.mesh;
    FaceSegmentGraph faces;faces.build(target);
    auto targetDomain=NativeLeafDomain::build(target,faces,domain_.axisymmetric());
    auto current=std::make_unique<InterfaceState>(target),previous=std::make_unique<InterfaceState>(target);
    const auto oldTransport=HydroInterfaceTransport::domainView(mesh_,domain_,flow_.boundaryPatches);
    const auto targetTransport=HydroInterfaceTransport::domainView(target,targetDomain,flow_.boundaryPatches);
    const InterfaceConservativeGeometry::BoundaryInput oldBoundary{oldTransport,{},{}},targetBoundary{targetTransport,{},{}};
    const auto transfer=[&](const InterfaceState& source,InterfaceState& destination) {
        const auto donors=InterfaceConservativeGeometry::reconstructTransportDonors(mesh_,domain_,source,columnRadius_>0?&oldBoundary:nullptr);
        return InterfaceMeshTransfer::prepareConservative(mesh_,domain_,source,target,targetDomain,donors.normals,donors.curved,
            {controls.remapScalarTolerance,controls.remapRelativeGeometryTolerance*source.conservativePhase().report().liquidAfter,
                true,true},destination,columnRadius_>0?&targetBoundary:nullptr);
    };
    candidate.currentGeometry=transfer(*state_,*current);candidate.previousGeometry=transfer(*previous_,*previous);
    candidate.hydro=hydro.prepareConservativeMeshTransfer(mesh_,domain_,*state_,*previous_,target,targetDomain,*current,*previous,
        flow_.densityNegative,flow_.densityPositive,flow_.viscosityNegative,flow_.viscosityPositive,&candidate.flow);
    candidate.adapter=std::unique_ptr<HydroConservativeInterface>(new HydroConservativeInterface(target,flow_,
        std::move(targetDomain),std::move(current),std::move(previous),columnRadius_,neckBreakup_,columnHalfWaves_));
    return candidate;
}

void HydroConservativeInterface::inheritOutput(HydroConservativeInterface& old) noexcept
{
    history_.swap(old.history_);columnHistory_.swap(old.columnHistory_);
    componentHistory_.swap(old.componentHistory_);
    eventHistory_.swap(old.eventHistory_);std::swap(writtenEvents_,old.writtenEvents_);
}

void HydroConservativeInterface::materials(Hydrodynamics& hydro,double time) const
{ materials(hydro,*state_,time); }
void HydroConservativeInterface::materials(Hydrodynamics& hydro,const InterfaceState& state,double time) const
{
    state.validateMesh(mesh_);
    if(state.stage()!=InterfaceState::Stage::Ready || state.progress().endTime!=time)
        throw std::invalid_argument("conservative materials require the supplied accepted/candidate interface time");
    hydro.updateMaterials(mesh_,{HydroStepGeometry::complete(domain_.axisymmetric()),time,HydroMaterialEvolution::TimeDependent,
        HydroConservativeMaterial{&domain_,&state.conservativePhase(),flow_.densityNegative,flow_.densityPositive,
            flow_.viscosityNegative,flow_.viscosityPositive}});
}
Vector2D HydroConservativeInterface::rootSpacing() const
{
    const auto extent=mesh_.maxCoord-mesh_.minCoord;
    return {extent.x/mesh_.Nx,extent.y/mesh_.Ny};
}
HydroFittedVelocitySampler HydroConservativeInterface::velocity(const Hydrodynamics& hydro) const
{
    const auto state=hydro.state();
    auto ux=state.velocityX.copy();auto uy=state.velocityY.copy();
    return HydroInterfaceTransport::completeVelocity(mesh_,domain_,ux,uy,flow_.boundaryPatches);
}
std::unique_ptr<InterfaceState> HydroConservativeInterface::prepare(const Hydrodynamics& hydro,const InterfaceStepTime& time,
    InterfaceConservativeGeometry::AdvanceReport* report,MaterialMotion* material) const
{
    if(material)*material={};
    const auto spacing=rootSpacing();
    const auto transport=HydroInterfaceTransport::domainView(mesh_,domain_,flow_.boundaryPatches);
    const auto frozen=velocity(hydro);const auto flux=HydroInterfaceTransport::normalFlux(frozen);
    const InterfaceMotionVelocityLayers layers{[&](const Vector2D& p){return frozen.sample(p);},
        [&](const Vector2D& p){return frozen.sample(p);},true};
    const InterfaceConservativeGeometry::BoundaryInput boundary{transport,[&](int row){return phases_.at(row);},
        [&](Boundary side,Vector2D p,double) {
            int sign=0;
            for(const auto face:domain_.faces(mesh_)) if(phases_[face.nativeRow]) {
                const auto faceSide=face.horizontal?(face.outwardSign<0?Boundary::south:Boundary::north)
                    :(face.outwardSign<0?Boundary::west:Boundary::east);
                const auto& edge=face.horizontal?mesh_.hEdges(face.edgeLabel):mesh_.vEdges(face.edgeLabel);
                const auto a=mesh_.verts(edge.edgeToVert[0]).coord(),b=mesh_.verts(edge.edgeToVert[1]).coord();
                const double s=face.horizontal?p.x:p.y,lo=face.horizontal?a.x:a.y,hi=face.horizontal?b.x:b.y;
                if(side==faceSide && s>=lo && s<=hi) {
                    if(sign && sign!=phases_[face.nativeRow]) throw std::invalid_argument("ambiguous conservative scalar inlet phase");
                    sign=phases_[face.nativeRow];
                }
            }
            if(!sign) throw std::invalid_argument("no conservative scalar inlet phase");
            return static_cast<double>(sign)*std::max(spacing.x,spacing.y);
        }};
    auto next=std::make_unique<InterfaceState>(mesh_);
    // Terminal provider selection below discards ordinary midpoint diagnostics.
    // Keep the actual force/curvature gates; avoid constructing an unused copy.
    using Diagnostics=InterfaceConservativeGeometry::CurvatureDiagnostics;
    std::optional<InterfaceConservativeGeometry::NeckEventPolicy> neckEvents;
    if(neckBreakup_)neckEvents=InterfaceConservativeGeometry::NeckEventPolicy{
        neckBreakup_->radiusCells,neckBreakup_->maximumLevel,material!=nullptr};
    InterfaceConservativeGeometry::advance(mesh_,domain_,*state_,time,layers,flux,flux,*next,
        std::max(4.,flow_.interfaceThickness/std::min(spacing.x,spacing.y)+2),{},report,!neckBreakup_ || material,&boundary,true,
        neckBreakup_?Diagnostics::Omit:Diagnostics::Record,neckEvents?&*neckEvents:nullptr);
    if(material)material->physical=next->conservativeGeometry().surfaceTransfer;
    if(neckBreakup_) {
        const auto neck=[&](const InterfaceState& state)->std::optional<std::pair<double,int>> {
            if(state.topology().surfaces().size()!=1)return {};
            const auto& chain=state.topology().surfaces()[0].geometry.nativeChain();
            if(chain.closed())return {};
            for(bool end:{false,true})if(const auto* section=chain.endpointSymmetryPlane(end))
                if(section->normalAxis==1 && section->coordinate==mesh_.maxCoord.y) {
                    const int g=end?static_cast<int>(chain.graphs().size())-1:0;
                    return std::make_pair(chain.frameAt(g,end?1:0).point.x,chain.graphs()[g].nativeCellIndex);
                }
            return {};
        };
        const auto old=neck(*state_),current=neck(*next);
        if(current) {
            const auto& cell=mesh_.cells(mesh_.cellIndexToLabel[current->second]);
            if(current->first<=neckBreakup_->radiusCells*cell.size().x) {
                if(cell.level()!=neckBreakup_->maximumLevel)
                    throw std::runtime_error("terminal neck reached before its declared finest level");
                if(!old || !(current->first<old->first))
                    throw std::runtime_error("terminal neck requires a shrinking physical source interval");
                auto event=InterfaceNeckCap::prepare(mesh_,domain_,*next,mesh_.maxCoord.y,neckBreakup_->radiusCells);
                if(material) {
                    material->event=event.materialTransfer();
                    next=std::move(event).takeMaterialCandidate(time,*material->event);
                } else next=std::move(event).takeUnchargedCandidate(time);
            }
        }
        next=InterfaceNeckCap::prepareTerminalCurvature(mesh_,domain_,std::move(next),
            mesh_.maxCoord.y,neckBreakup_->maximumLevel,true);
    }
    return next;
}
std::optional<InterfaceColumnMode> HydroConservativeInterface::columnMode(const InterfaceState& state) const
{
    if(columnRadius_==0)return {};
    state.validateMesh(mesh_);
    return InterfaceColumnMode::evaluate(state,mesh_.minCoord.y,mesh_.maxCoord.y,columnRadius_,columnHalfWaves_);
}
void HydroConservativeInterface::write(const Hydrodynamics& hydro,int step,double time,const std::string& prefix,bool fields,double forceTime,
    const HydroSolveReport* report,const Quantum::Hydrodynamics::ConservativeCapillaryReport* capillary,const InterfaceColumnMode* columnMode)
{
    const auto& phase=state_->conservativePhase();
    const auto material=hydro.state().materials;
    if(!material || material->sourceTime!=time || state_->progress().endTime!=time || state_->progress().completedSteps!=step)
        throw std::logic_error("conservative output requires accepted geometry and material at the same case time");
    if(columnRadius_>0 && (!columnMode || columnMode->source!=state_->geometry()))
        throw std::logic_error("column output requires a diagnostic from the accepted geometry");
    double geometric=0;for(const auto& component:state_->topology().liquids()) geometric+=component.volume;
    if(!history_.is_open()) {
        history_.open(OutputLayout::series()/"interface_conservative.csv");
        history_<<"step,accepted_time,conserved_volume,geometric_volume,inventory_step_closure,geometric_inventory_l1,dual_mass_relative_closure,conservative_momentum,material_time,force_time";
        if(neckBreakup_)history_<<",event_curvature_cells,event_inventory_fraction,event_geometry_l1_fraction,event_count,event_force_faces,max_event_force,event_projection_force_spread,bounded_tail_ridge_faces,max_scaled_tail_force_uncertainty,max_tail_force_uncertainty,surface_integrated_faces,max_normal_moment_defect,event_inventory_iterations,event_inventory_backtracks,event_inventory_cut_cells,event_inventory_initial_cut_error,event_inventory_final_cut_error,event_inventory_pure_tail_error";
        history_<<'\n';
    }
    long double eventVolume=0,eventDifference=0;
    for(int c:state_->eventCurvatureCells()) {
        eventVolume+=phase.liquidVolumes()[c];
        if(state_->hasConservativeGeometry())eventDifference+=std::abs(
            state_->conservativeGeometry().geometricVolumeByCell.at(c)-phase.liquidVolumes()[c]);
    }
    history_<<std::setprecision(17)<<step<<','<<time<<','<<phase.report().liquidAfter<<','<<geometric<<','<<phase.report().closure
        <<','<<(state_->hasConservativeGeometry()?state_->conservativeGeometry().l1Difference:0)
        <<','<<(report?report->conservativeMassRelativeClosure:0)<<','<<(report?report->conservativeMomentumUsed:false)
        <<','<<material->sourceTime<<','<<forceTime;
    if(neckBreakup_)history_<<','<<state_->eventCurvatureCells().size()
        <<','<<eventVolume/phase.report().liquidAfter<<','<<eventDifference/phase.report().liquidAfter
        <<','<<state_->topology().events().size()<<','<<(capillary?capillary->eventModelFaces:0)
        <<','<<(capillary?capillary->maxEventModelForce:0)<<','<<(capillary?capillary->maxEventProjectionForceSpread:0)
        <<','<<(capillary?capillary->boundedTailRidgeFaces:0)<<','<<(capillary?capillary->maxScaledTailForceUncertainty:0)
        <<','<<(capillary?capillary->maxTailForceUncertainty:0)
        <<','<<(capillary?capillary->surfaceIntegratedFaces:0)<<','<<(capillary?capillary->maxNormalMomentDefect:0);
    if(neckBreakup_) {
        const auto* correction=state_->hasConservativeGeometry()?&state_->conservativeGeometry():nullptr;
        history_<<','<<(correction?correction->eventInventoryIterations:0)
            <<','<<(correction?correction->eventInventoryBacktracks:0)
            <<','<<(correction?correction->eventInventoryCutCells:0)
            <<','<<(correction?correction->eventInventoryInitialCutError:0)
            <<','<<(correction?correction->eventInventoryFinalCutError:0)
            <<','<<(correction?correction->eventInventoryPureTailError:0);
    }
    history_<<'\n';
    history_.flush();if(!history_) throw std::runtime_error("cannot write conservative Hydro history");
    if(!componentHistory_.is_open()) {
        componentHistory_.open(OutputLayout::series()/"interface_components.csv");
        componentHistory_<<"step,accepted_time,liquid_id,inventory_volume,geometric_volume,touches_boundary,material_inventory_volume,material_transport_iterations,material_mass_residual,material_row_residual\n";
    }
    // Topology owns the liquid regions. Inventory in cells outside those
    // regions remains a separate row (reserved ID zero), never assigned to
    // the closest drop or silently included in its volume. The final zero-ID
    // row also closes this accepted time group for readers of a live file.
    std::vector<std::uint8_t> assigned(mesh_.cellIndexNum(),0);
    for(const auto& liquid:state_->topology().liquids()) {
        long double inventory=0;
        for(int c:liquid.cells){assigned[c]=1;inventory+=phase.liquidVolumes()[c];}
        const auto& components=state_->componentInventories();
        const auto material=std::find_if(components.begin(),components.end(),[&](const auto& item){return item.liquidId==liquid.id;});
        long double materialVolume=std::numeric_limits<long double>::quiet_NaN();
        if(material!=components.end()) {materialVolume=0;for(int c:domain_.cellIndexByRow())materialVolume+=material->liquidVolumes[c];}
        const auto* diagnostics=state_->hasConservativeGeometry()?&state_->conservativeGeometry():nullptr;
        componentHistory_<<std::setprecision(17)<<step<<','<<time<<','<<liquid.id<<','<<inventory
            <<','<<liquid.volume<<','<<liquid.touchesBoundary<<','<<materialVolume
            <<','<<(diagnostics?diagnostics->componentTransportIterations:0)
            <<','<<(diagnostics?diagnostics->componentMassResidual:0)
            <<','<<(diagnostics?diagnostics->componentRowResidual:0)<<'\n';
    }
    long double unassigned=0;
    for(int c:domain_.cellIndexByRow())if(!assigned[c])unassigned+=phase.liquidVolumes()[c];
    componentHistory_<<std::setprecision(17)<<step<<','<<time<<",0,"<<unassigned<<",0,0,nan,0,0,0\n";
    componentHistory_.flush();if(!componentHistory_)throw std::runtime_error("cannot write accepted liquid component history");
    if(neckBreakup_) {
        if(!eventHistory_.is_open())eventHistory_.open(OutputLayout::series()/"interface_events.jsonl");
        const auto& events=state_->topology().events();
        if(writtenEvents_>events.size())throw std::logic_error("accepted interface event history regressed");
        const auto ids=[&](const std::vector<std::uint64_t>& values) {
            eventHistory_<<'[';bool comma=false;
            for(auto id:values){if(comma)eventHistory_<<',';eventHistory_<<id;comma=true;}
            eventHistory_<<']';
        };
        for(std::size_t i=writtenEvents_;i<events.size();++i) {
            const auto& event=events[i];const char* kind=nullptr;
            switch(event.kind) {
                case InterfaceTopologyEvent::Kind::Split:kind="split";break;
                case InterfaceTopologyEvent::Kind::Merge:kind="merge";break;
                case InterfaceTopologyEvent::Kind::Reconnection:kind="reconnection";break;
                case InterfaceTopologyEvent::Kind::SectionDetachment:kind="section_detachment";break;
            }
            if(!kind)throw std::logic_error("unknown accepted interface event kind");
            eventHistory_<<std::setprecision(17)<<"{\"accepted_step\":"<<step<<",\"accepted_time\":"<<time
                <<",\"event_index\":"<<i<<",\"kind\":\""<<kind<<"\",\"begin\":"<<event.begin<<",\"end\":"<<event.end
                <<",\"parent_liquid_ids\":";ids(event.parents);eventHistory_<<",\"child_liquid_ids\":";ids(event.children);
            eventHistory_<<",\"detached_sections\":[";bool comma=false;
            for(const auto& section:event.detachedSections) {
                if(comma)eventHistory_<<',';comma=true;
                eventHistory_<<"{\"normal_axis\":"<<section.normalAxis<<",\"coordinate\":"<<section.coordinate<<'}';
            }
            eventHistory_<<"],\"interior_necks\":[";comma=false;
            for(const auto& neck:event.interiorNecks) {
                if(comma)eventHistory_<<',';comma=true;
                eventHistory_<<"{\"section\":"<<neck.section<<",\"radial_width\":"<<neck.radialWidth
                    <<",\"axial_width\":"<<neck.axialWidth<<",\"cutoff_cells\":"<<neck.cutoffCells
                    <<",\"source_radius\":"<<neck.sourceRadius<<",\"predicted_radius\":"<<neck.predictedRadius
                    <<",\"gas_value\":"<<neck.gasValue<<",\"cells\":[";
                bool separate=false;for(int c:neck.cells){if(separate)eventHistory_<<',';separate=true;eventHistory_<<c;}
                eventHistory_<<"]}";
            }
            eventHistory_<<"],\"surfaces\":[";comma=false;
            for(const auto& surface:state_->topology().surfaces()) {
                if(comma)eventHistory_<<',';comma=true;
                eventHistory_<<"{\"id\":"<<surface.id<<",\"liquid_id\":"<<surface.liquidId<<",\"parents\":";
                ids(surface.parents);eventHistory_<<'}';
            }
            eventHistory_<<"],\"liquids\":[";comma=false;
            for(const auto& liquid:state_->topology().liquids()) {
                long double inventory=0;for(int c:liquid.cells)inventory+=phase.liquidVolumes()[c];
                if(comma)eventHistory_<<',';comma=true;
                eventHistory_<<"{\"id\":"<<liquid.id<<",\"assigned_inventory_volume\":"<<inventory
                    <<",\"geometric_volume\":"<<liquid.volume<<'}';
            }
            eventHistory_<<"],\"total_inventory_volume\":"<<phase.report().liquidAfter
                <<",\"total_geometric_volume\":"<<geometric<<",\"unassigned_inventory_volume\":"
                <<state_->conservativeGeometry().unassignedInventoryVolume<<"}\n";
        }
        eventHistory_.flush();if(!eventHistory_)throw std::runtime_error("cannot write accepted interface event history");
        writtenEvents_=events.size();
    }
    if(columnRadius_>0) {
        if(!columnHistory_.is_open()) {
            columnHistory_.open(OutputLayout::series()/"column_mode.csv");
            columnHistory_<<"step,accepted_time,single_radius_graph,axial_quadrature_length,mean_radius,cosine_amplitude,relative_amplitude\n";
        }
        columnHistory_<<std::setprecision(17)<<step<<','<<time<<','<<columnMode->singleRadiusGraph<<','<<columnMode->axialLength
            <<','<<columnMode->meanRadius<<','<<columnMode->cosineAmplitude<<','<<columnMode->relativeAmplitude<<'\n';
        columnHistory_.flush();if(!columnHistory_)throw std::runtime_error("cannot write column mode history");
    }
    if(!fields) return;
    writeFields(hydro,step,time,prefix,forceTime,capillary);
}

void HydroConservativeInterface::writeEventSourceFrame(const Hydrodynamics& hydro,const std::string& prefix) const
{
    const auto progress=state_->progress();
    writeFields(hydro,progress.completedSteps,progress.endTime,prefix,progress.endTime,nullptr);
}

void HydroConservativeInterface::writeFields(const Hydrodynamics& hydro,int step,double time,const std::string& prefix,
    double forceTime,const Quantum::Hydrodynamics::ConservativeCapillaryReport* capillary) const
{
    const auto material=hydro.state().materials;
    if(!material || material->sourceTime!=time || state_->progress().endTime!=time || state_->progress().completedSteps!=step)
        throw std::logic_error("interface frame requires an accepted geometry/material clock");
    const auto& phase=state_->conservativePhase();
    using namespace TecplotWriter;
    Zone zone;zone.title="Accepted conservative Hydro state";
    zone.coordinateNames=domain_.axisymmetric()?std::array<std::string,2>{"r_code","z_code"}:std::array<std::string,2>{"x_code","y_code"};
    zone.solutionTime=time;zone.strandId=1;
    zone.auxiliary={{"TimeRole","completed_step_time"},{"MaterialSource","conserved_C"},{"Transport","first_order_frozen_NS"}};
    zone.fields={{"phi",Location::Nodal,{}},{domain_.axisymmetric()?"ur":"ux",Location::Nodal,{}},{domain_.axisymmetric()?"uz":"uy",Location::Nodal,{}},
        {"C",Location::CellCentered,{}},{"rho",Location::CellCentered,{}},{"mu",Location::CellCentered,{}},
        {"p",Location::CellCentered,{}},{"divU",Location::CellCentered,{}},{"C_Gamma",Location::CellCentered,{}}};
    for(const auto& component:state_->componentInventories())
        zone.fields.push_back({"C_liquid_"+std::to_string(component.liquidId),Location::CellCentered,{}});
    const auto solved=hydro.state();const auto frozen=velocity(hydro);
    for(int k=0;k<mesh_.vertIndexNum();++k) {
        const auto p=mesh_.verts(mesh_.vertIndexToLabel[k]).coord(),u=frozen.sample(p);zone.points.push_back(p);
        zone.fields[0].values.push_back(state_->displayValue(k));zone.fields[1].values.push_back(u.x);zone.fields[2].values.push_back(u.y);
    }
    for(int c:domain_.cellIndexByRow()) {
        std::vector<int> corners;for(int label:mesh_.cells(mesh_.cellIndexToLabel[c]).cellToVert) corners.push_back(mesh_.verts(label).index());
        zone.connectivity.push_back(std::move(corners));
        zone.fields[3].values.push_back(static_cast<double>(phase.liquidVolumes()[c]/phase.capacities()[c]));
        zone.fields[4].values.push_back(solved.density[c]);zone.fields[5].values.push_back(solved.viscosity[c]);
        zone.fields[6].values.push_back(solved.physicalPressure[c]);zone.fields[7].values.push_back(solved.divergence[c]);
        const auto volume=state_->hasConservativeGeometry()?state_->conservativeGeometry().geometricVolumeByCell[c]
            :std::visit([&](const auto& measure){return measure.negativeVolumeByCell[c];},state_->phase().volumeMeasure);
        zone.fields[8].values.push_back(static_cast<double>(volume/phase.capacities()[c]));
        for(std::size_t k=0;k<state_->componentInventories().size();++k)
            zone.fields[9+k].values.push_back(static_cast<double>(state_->componentInventories()[k].liquidVolumes[c]/phase.capacities()[c]));
    }
    std::ostringstream name;name<<prefix<<"_step"<<std::setw(4)<<std::setfill('0')<<step<<".dat";
    TecplotWriter::write(OutputLayout::fields()/name.str(),zone);
    // Sample the accepted native zero surface itself. A contour extracted
    // from the displayed nodal phi would be a different, lower-order shape.
    std::ostringstream surfaceName;surfaceName<<prefix<<"_interface_geometry_step"<<std::setw(4)<<std::setfill('0')<<step<<".csv";
    std::ofstream surfaceFile(OutputLayout::fields()/surfaceName.str());
    surfaceFile<<"accepted_step,accepted_time,surface_id,graph,fraction,x,y,normal_x,normal_y,physical_weight\n";
    for(const auto& surface:state_->topology().surfaces()) {
        const auto& chain=surface.geometry.nativeChain();const auto& panels=surface.geometry.panels();
        for(std::size_t g=0;g<chain.graphs().size();++g)for(double fraction:{0.,.5,1.}) {
            const auto frame=chain.frameAt(static_cast<int>(g),fraction);
            const double weight=panels[g].physicalMeasure*(fraction==.5?2./3:1./6);
            surfaceFile<<std::setprecision(17)<<step<<','<<time<<','<<surface.id<<','<<g<<','<<fraction
                <<','<<frame.point.x<<','<<frame.point.y<<','<<frame.normal.x<<','<<frame.normal.y<<','<<weight<<'\n';
        }
    }
    surfaceFile.flush();if(!surfaceFile)throw std::runtime_error("cannot write accepted conservative interface geometry");
    if(state_->hasConservativeGeometry()) {
        std::ostringstream jetName;jetName<<prefix<<"_interface_jets_step"<<std::setw(4)<<std::setfill('0')<<step<<".csv";
        std::ofstream jets(OutputLayout::fields()/jetName.str());
        jets<<"accepted_step,accepted_time,vertex,x,y,phi,dx_before,dy_before,dxy_before,dx_mixed,dy_mixed,dxy_mixed,dx_pure,dy_pure,dxy_pure,source,rejected_blend,rejected_height,rejected_graph,scalar_correction\n";
        for(const auto& jet:state_->conservativeGeometry().nodalJets) {
            const auto p=mesh_.verts(mesh_.vertIndexToLabel[jet.vertex]).coord();
            jets<<std::setprecision(17)<<step<<','<<time<<','<<jet.vertex<<','<<p.x<<','<<p.y<<','<<jet.value;
            for(const auto& stage:{jet.before,jet.mixed,jet.pure}) for(double value:stage) jets<<','<<value;
            jets<<','<<jet.source;for(double value:jet.rejectedValues) {jets<<',';if(std::isfinite(value)) jets<<value;}
            jets<<','<<jet.scalarCorrection;
            jets<<'\n';
        }
        jets.flush();if(!jets) throw std::runtime_error("cannot write conservative interface jet samples");
    }
    if(capillary) {
        std::ostringstream diagnosticName;diagnosticName<<prefix<<"_capillary_source_step"<<std::setw(4)<<std::setfill('0')<<step<<".csv";
        std::ofstream diagnostic(OutputLayout::fields()/diagnosticName.str());
        diagnostic<<"accepted_step,accepted_time,force_time,native_row,face_x,face_y,interface_x,interface_y,gradient_C,curvature,meridional,azimuthal,force,primary_height_function,primary_even_axis,primary_implicit_axis,interpolated,primary_method,secondary_method,secondary_weight,secondary_even_axis,secondary_implicit_axis,unresolved_inventory_tail,projection_distance,max_minority_fraction,same_chart_tail_candidates,regularized_curvature,secondary_curvature_weight,regularization_status";
        if(neckBreakup_)diagnostic<<",event_chart_candidates,event_curvature_spread,bounded_tail_ridge,tail_curvature_spread,scaled_tail_force_uncertainty,tail_alternatives,liquid_id,surface_integrated,surface_samples,normal_moment_defect";
        diagnostic<<'\n';
        for(const auto& face:capillary->faces) {diagnostic<<std::setprecision(17)<<step<<','<<time<<','<<forceTime<<','<<face.nativeRow
            <<','<<face.center.x<<','<<face.center.y<<','<<face.interfacePoint.x<<','<<face.interfacePoint.y
            <<','<<face.gradient<<','<<face.curvature<<','<<face.meridional<<','<<face.azimuthal<<','<<face.force
            <<','<<face.height<<','<<face.evenAxis<<','<<face.implicitAxis
            <<','<<face.interpolated<<','<<face.primaryMethod<<','<<face.secondaryMethod<<','<<face.secondaryWeight
            <<','<<face.secondaryEvenAxis<<','<<face.secondaryImplicitAxis
            <<','<<face.unresolvedInventoryTail<<','<<face.projectionDistance<<','<<face.maxMinorityFraction
            <<','<<face.sameChartTailCandidates<<','<<face.regularizedCurvature<<','<<face.secondaryCurvatureWeight
            <<','<<face.regularizationStatus;
            if(neckBreakup_)diagnostic<<','<<face.eventChartCandidates<<','<<face.eventCurvatureSpread
                <<','<<face.boundedTailRidge<<','<<face.tailCurvatureSpread<<','<<face.scaledTailForceUncertainty<<','<<face.tailAlternatives<<','<<face.liquidId
                <<','<<face.surfaceIntegrated<<','<<face.surfaceSamples<<','<<face.normalMomentDefect;
            diagnostic<<'\n';
        }
        diagnostic.flush();if(!diagnostic) throw std::runtime_error("cannot write conservative capillary source samples");
    }
}
