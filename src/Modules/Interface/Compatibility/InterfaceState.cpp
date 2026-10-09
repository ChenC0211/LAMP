#include "Modules/Interface/Compatibility/Detail/InterfaceStateInternal.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include <utility>

InterfaceState::Data::Data(const QuadtreeMesh& mesh)
    : freeLevelSet(mesh, 0.0), activeVertex(mesh.vertIndexNum(), 0) {}
InterfaceState::InterfaceState(const QuadtreeMesh& mesh) : data_(std::make_unique<Data>(mesh)) {}
InterfaceState::InterfaceState()=default;
InterfaceState::~InterfaceState() = default;
InterfaceState::InterfaceState(InterfaceState&&) noexcept = default;
InterfaceState& InterfaceState::operator=(InterfaceState&&) noexcept = default;
InterfaceState::Stage InterfaceState::stage() const noexcept { return parametric_?Stage::Ready:(data_ ? data_->stage : Stage::Failed); }
bool InterfaceState::hasParametricGeometry() const noexcept {return bool(parametric_);}
const InterfaceParametricProjection::Candidate& InterfaceState::parametricGeometry() const {
    if(!parametric_)throw std::logic_error("Interface state has no parametric geometry");
    return parametric_->candidate;
}
const InterfaceParametricTransport::FaceVolumes& InterfaceState::parametricFluxes() const {
    if(!parametric_ || !parametric_->previousGeometry)throw std::logic_error("Interface state has no parametric transport step");
    return parametric_->fluxes;
}
const std::vector<InterfaceParametricProjection::NeckEvent>& InterfaceState::parametricNeckEvents() const {
    if(!parametric_)throw std::logic_error("Interface state has no parametric event history");
    return parametric_->neckEvents;
}
bool InterfaceState::parametricModelCurvature(std::uint64_t component,double parameter) const {
    return parametricModelEvent(component,parameter).has_value();
}
std::optional<std::uint64_t> InterfaceState::parametricModelEvent(std::uint64_t component,double parameter) const {
    if(!parametric_||!std::isfinite(parameter)||parameter<0||parameter>1)throw std::invalid_argument("invalid parametric material query");
    bool present=false;for(const auto& part:parametric_->candidate.boundaryHandle()->components())present=present||part.id==component;
    if(!present)throw std::invalid_argument("model-curvature query references an absent component");
    for(auto it=parametric_->neckEvents.rbegin();it!=parametric_->neckEvents.rend();++it){const auto& r=it->report;
        for(int child=0;child<2;child++)if(r.children[child]==component){
            const auto cap=r.capIntervals[child];if(parameter>=cap[0]&&parameter<=cap[1])return r.parent;
            const double begin=child?cap[1]:0,end=child?1:cap[0];const auto retained=r.retainedParentIntervals[child];
            parameter=retained[0]+(parameter-begin)/(end-begin)*(retained[1]-retained[0]);component=r.parent;break;
        }
    }
    return {};
}
bool InterfaceState::parametricFollows(const InterfaceState& previous,const InterfaceStepTime& time) const noexcept {
    if(!parametric_ || !previous.parametric_)return false;
    const auto& p=*parametric_;const auto& old=*previous.parametric_;const auto& t=p.stepTime;
    return p.binding==old.binding && p.previousGeometry==old.candidate.boundaryHandle() &&
        t.index==time.index && t.begin==time.begin && t.end==time.end && t.dt==time.dt &&
        t.index==old.progress.completedSteps+1 && t.begin==old.progress.endTime;
}
InterfaceStepProgress InterfaceState::progress() const {
    if(parametric_)return parametric_->progress;
    const auto& data = InterfaceStateAccess::get(*this);
    if (data.stage != Stage::Ready) throw std::logic_error("Interface progress requires ready geometry");
    return data.progress;
}
void InterfaceState::validateMesh(const QuadtreeMesh& mesh) const {
    if(parametric_){if(parametric_->binding->mesh!=&mesh)throw std::logic_error("parametric geometry belongs to another mesh owner");parametric_->binding->identity->validate(mesh);return;}
    const auto& data = InterfaceStateAccess::get(*this);
    if ((data.stage != Stage::Ready && data.stage != Stage::Prepared) || !data.meshIdentity)
        throw std::logic_error("Interface state has no usable mesh binding");
    data.meshIdentity->validate(mesh);
}
bool InterfaceState::hasPhase() const noexcept {
    return data_ && data_->stage == Stage::Ready && data_->phase.has_value();
}
const InterfacePhaseState& InterfaceState::phase() const {
    if (!hasPhase()) throw std::logic_error("Interface state has no prepared phase measures/fluxes");
    return *data_->phase;
}
bool InterfaceState::hasConservativePhase() const noexcept {
    return data_ && data_->stage==Stage::Ready && data_->conservativePhase.has_value();
}
const InterfaceConservativePhase& InterfaceState::conservativePhase() const {
    if (!hasConservativePhase()) throw std::logic_error("Interface has no conservative phase inventory");
    return *data_->conservativePhase;
}
const std::vector<InterfaceComponentInventory>& InterfaceState::componentInventories() const {
    return InterfaceStateAccess::get(*this).componentInventories;
}
bool InterfaceState::hasConservativeGeometry() const noexcept {
    return data_ && data_->stage==Stage::Ready && data_->conservativeGeometry.has_value();
}
const std::vector<int>& InterfaceState::eventCurvatureCells() const {
    return InterfaceStateAccess::get(*this).eventCurvatureCells;
}
const InterfaceConservativeGeometryReport& InterfaceState::conservativeGeometry() const {
    if (!hasConservativeGeometry()) throw std::logic_error("Interface has no conservative geometry diagnostics");
    return *data_->conservativeGeometry;
}
void InterfaceState::initializeConservativePhase(const QuadtreeMesh& mesh,const NativeLeafDomain& domain) {
    auto& data=InterfaceStateAccess::get(*this);
    if (!hasPhase() || data.conservativePhase || data.progress.completedSteps!=0 || data.progress.endTime!=0)
        throw std::logic_error("conservative initialization requires time-zero geometry and no existing inventory");
    validateMesh(mesh);
    InterfacePhaseGeometry::validateLayout(mesh,domain,*data.geometry());
    auto inventory=std::visit([&](const auto& measure) {
        auto volumes=measure.negativeVolumeByCell;
        // Pure cells use their geometric sign certificate and the inventory's
        // own full capacity. Summing remote boundary moments can differ by a
        // few ulps from that capacity even for an entirely filled rectangle.
        // This is initial geometric classification, never clipping a transported
        // value. Mixed cells retain their raw measured volume and must pass the
        // strict bounds check. phase() retains every original integral for audit.
        for (std::size_t row=0;row<domain.cellIndexByRow().size();++row) {
            const int cell=domain.cellIndexByRow()[row];
            const int sign=data.topology ? data.topology->cellPhase()[cell]
                : data.geometry()->closedCellSignExceptVertices(
                    mesh.cells(mesh.cellIndexToLabel[cell]).centerCoord(),15u);
            if (sign==-1) volumes[cell]=domain.cellVolumes()[row];
            else if (sign==1) volumes[cell]=0;
        }
        return InterfaceConservativePhase::initialize(mesh,domain,volumes);
    },data.phase->volumeMeasure);
    data.conservativePhase.emplace(std::move(inventory));
}
void InterfaceState::initializePhase(const QuadtreeMesh& mesh, const InterfaceTransportDomain& transport,
    const InterfaceNormalFluxView& flux) {
    auto& data = InterfaceStateAccess::get(*this);
    if (data.stage != Stage::Ready || data.phase || !data.geometry() || (!data.topology && !data.panels()))
        throw std::logic_error("phase initialization requires ready geometry without an existing phase state");
    validateMesh(mesh);
    const auto& domain = transport.geometry();
    InterfacePhaseGeometry::MappedFacePartitions partitions(mesh, domain, data.geometry());
    InterfacePhaseState prepared;
    if (data.topology) {
        if (data.geometry()->reconstruction()==LevelSetHermiteGeometry::Reconstruction::EnrichedQuartic)
            prepared.volumeMeasure = InterfaceCellMeasure::evaluateQuarticRefined(mesh,domain,partitions,*this);
        else prepared.volumeMeasure = InterfaceCellMeasure::evaluateRefined(mesh,domain,partitions,*this);
    }
    else if (data.geometry()->usesNativeLeaves())
        prepared.volumeMeasure = InterfaceCellMeasure::evaluateQuarticRefined(mesh, domain, partitions, data.panels()->nativeChain());
    else
        prepared.volumeMeasure = InterfaceCellMeasure::evaluateRefined(mesh, domain, partitions, *data.panels());
    prepared.phaseBudget = InterfacePhaseFlux::phaseFluxBudget(mesh, transport, partitions, flux);
    data.phase.emplace(std::move(prepared));
}
int InterfaceState::vertexCount() const { return InterfaceStateAccess::get(*this).freeLevelSet.size(); }
double InterfaceState::displayValue(int vertex) const {
    const auto& field = InterfaceStateAccess::get(*this).freeLevelSet;
    if (vertex < 0 || vertex >= field.size()) throw std::out_of_range("Interface display vertex");
    return field.value[vertex];
}
const std::vector<std::uint8_t>& InterfaceState::activeVertices() const { return InterfaceStateAccess::get(*this).activeVertex; }
const std::vector<int>& InterfaceState::pinnedVertices() const { return InterfaceStateAccess::get(*this).pinnedVertexIndices; }
const std::vector<Vector2D>& InterfaceState::pinnedCoordinates() const { return InterfaceStateAccess::get(*this).pinnedCoordinates; }
const std::vector<InterfaceInitialContactJet>& InterfaceState::initialContactJets() const { return InterfaceStateAccess::get(*this).initialContactJets; }
const InterfaceScalarPatch& InterfaceState::scalarPatch() const { return InterfaceStateAccess::get(*this).scalarPatch; }
const std::vector<Quantum::Geometry::ActiveHermitePatch::NodalDof>& InterfaceState::nodalData() const {
    return std::visit([](const auto& patch) -> const auto& { return patch.nodes; }, scalarPatch());
}
const std::vector<Quantum::Geometry::ActiveHermitePatch::CornerSector>& InterfaceState::cornerSectors() const {
    return std::visit([](const auto& patch) -> const auto& { return patch.cornerSectors; }, scalarPatch());
}
const std::shared_ptr<const LevelSetHermiteGeometry::Snapshot>& InterfaceState::geometry() const { return InterfaceStateAccess::get(*this).source; }
const LevelSetInterfacePanels* InterfaceState::panels() const { return InterfaceStateAccess::get(*this).panels(); }
const LevelSetInterfacePanels* InterfaceState::Data::panels() const {
    if (topology) {
        if (topology->surfaces().size() != 1)
            throw std::logic_error("single-interface consumer requires exactly one surface component");
        const auto& panels = topology->surfaces().front().geometry;
        panels.cellToPanel(); // Reject several arcs per cell for legacy consumers.
        return &panels;
    }
    return panelData.get();
}
bool InterfaceState::hasTopology() const noexcept {
    return data_ && data_->stage == Stage::Ready && static_cast<bool>(data_->topology);
}
const InterfaceTopology& InterfaceState::topology() const {
    if (!hasTopology()) throw std::logic_error("Interface state has no complete topology inventory");
    return *data_->topology;
}
void InterfaceState::Data::retainGeometry(std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> geometry,
    std::unique_ptr<LevelSetInterfacePanels>&& panels) {
    if (!geometry || !panels || !panels->usesSource(*geometry))
        throw std::invalid_argument("interface scalar and panels must retain the same immutable source");
    source = std::move(geometry); panelData = std::move(panels);
    topology.reset();
    phase.reset();
    conservativePhase.reset();
    componentInventories.clear();
    conservativeGeometry.reset();
}
