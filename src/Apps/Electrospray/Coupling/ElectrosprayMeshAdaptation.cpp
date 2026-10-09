#include "Electrospray/Coupling/ElectrosprayMeshAdaptation.h"
#include "Modules/Hydrodynamics/HydroInterfaceTransport.h"
#include "Electrospray/Setup/ElectrosprayModuleInputs.h"

#include "Electrospray/Setup/ElectrosprayElectrostatics.h"
#include "Electrospray/Setup/ElectrosprayHydroDomain.h"
#include "Electrospray/Coupling/ElectrosprayConducting.h"
#include "Core/Geometry/ActiveQuarticHermiteTransfer.h"
#include "Modules/Interface/Compatibility/InterfaceMeshTransfer.h"
#include "Core/Mesh/QuadtreeMeshRefiner.h"
#include "Core/Mesh/NativeMeshIdentity.h"
#include "Modules/Hydrodynamics/Hydrodynamics.h"
#include "Modules/Hydrodynamics/HydroFittedDomain.h"
#include "Modules/Hydrodynamics/HydroFittedVelocity.h"
#include "Modules/Hydrodynamics/HydroFittedRefinement.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <type_traits>

namespace
{
using Patch = Quantum::Geometry::ActiveQuarticHermitePatch;
using Domain = HydroFittedDomain;
void require(bool condition, const char* message)
{
    if (!condition) throw std::runtime_error(std::string("[ElectrosprayMeshAdaptation] ") + message);
}

int sourceIndex(const QuadtreeMesh& mesh, Vector2D point)
{
    const int label = mesh.findCellLabel(point);
    require(label >= 0 && label < mesh.cellLabelCapacity(), "target point is outside source mesh");
    const int index = mesh.cells(label).index();
    require(index >= 0 && index < mesh.cellIndexNum(), "source lookup did not return a leaf");
    return index;
}

std::array<const void*,7> owners(const ElectrosprayAcceptedStateView& s)
{
    return {&s.mesh,&s.domain,&s.binding,&s.interfaceState,&s.hydro,&s.conducting,&s.materials};
}
}

struct ElectrosprayMeshCandidate::Data
{
    explicit Data(const ElectrosprayAcceptedStateView& s)
        : sourceOwners(owners(s)), sourceMesh(s.domain.geometry().meshIdentity()),
          sourceInterface(s.interfaceState.surface.geometry()), sourceInterfaceTime(s.interfaceState.surface.progress()),
          sourceElectrical(s.conducting.electrical), sourceForces(s.conducting.forces.total.source())
    {
        const auto flow=s.hydro.state();
        require(flow.initialStateSet && flow.materials.has_value(),"AMR source lacks initialized Hydro material history");
        sourceFlow=flow.progress;sourceMaterial=*flow.materials;
    }
    std::array<const void*,7> sourceOwners;
    std::shared_ptr<const NativeMeshIdentity> sourceMesh;
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> sourceInterface;
    InterfaceStepProgress sourceInterfaceTime;
    Quantum::Electricity::FittedChargeState sourceElectrical;
    std::shared_ptr<const HydroForceSource> sourceForces;
    HydroStepProgress sourceFlow;
    HydroMaterialState sourceMaterial;
    bool pending=true;
    std::unique_ptr<QuadtreeMesh> mesh;
    std::unique_ptr<HydroFittedDomain> domain;
    std::unique_ptr<ElectrosprayElectricState> binding;
    std::unique_ptr<ElectrosprayInterfaceState> interfaceState;
    std::unique_ptr<Hydrodynamics> hydro;
    std::unique_ptr<ElectrosprayConductingState> conducting;
    ElectrosprayMaterialState materials;
    ElectrosprayMeshReport report;
};
ElectrosprayMeshCandidate::ElectrosprayMeshCandidate(const ElectrosprayAcceptedStateView& source)
    : data_(std::make_unique<Data>(source)) {}
ElectrosprayMeshCandidate::~ElectrosprayMeshCandidate() = default;
bool ElectrosprayMeshCandidate::pending() const noexcept { return data_ && data_->pending; }
const ElectrosprayMeshReport& ElectrosprayMeshCandidate::report() const { return data_->report; }
const ElectrosprayInterfaceState& ElectrosprayMeshCandidate::targetInterfaceState() const {
    require(pending(),"mesh candidate has already been committed");return *data_->interfaceState;
}
const Hydrodynamics& ElectrosprayMeshCandidate::targetHydro() const {
    require(pending(),"mesh candidate has already been committed");return *data_->hydro;
}
const Quantum::Electricity::FittedChargeState& ElectrosprayMeshCandidate::targetElectrical() const {
    require(pending(),"mesh candidate has already been committed");return data_->conducting->electrical;
}
void ElectrosprayMeshCandidate::validateSource(const ElectrosprayAcceptedStateView& source) const
{
    require(pending(),"mesh candidate has already been committed");
    const auto& d=*data_;
    require(owners(source)==d.sourceOwners,"mesh candidate source owners changed after preparation");
    d.sourceMesh->validate(source.mesh);
    const auto flow=source.hydro.state();const auto& p=flow.progress;const auto& old=d.sourceFlow;
    require(flow.initialStateSet && !p.stopped && !p.pendingAcceptance && p.geometry==old.geometry
        && p.completedSteps==old.completedSteps && p.acceptedSteps==old.acceptedSteps
        && p.endTime==old.endTime && p.acceptedEndTime==old.acceptedEndTime && flow.materials
        && flow.materials->revision==d.sourceMaterial.revision && flow.materials->sourceTime==d.sourceMaterial.sourceTime
        && flow.materials->evolution==d.sourceMaterial.evolution,"mesh candidate source Hydro history changed after preparation");
    const auto surface=source.interfaceState.surface.progress();
    require(source.interfaceState.surface.geometry()==d.sourceInterface
        && surface.completedSteps==d.sourceInterfaceTime.completedSteps && surface.endTime==d.sourceInterfaceTime.endTime,
        "mesh candidate source Interface history changed after preparation");
    require(&source.conducting.electrical.electric()==&d.sourceElectrical.electric()
        && source.conducting.forces.total.source()==d.sourceForces,"mesh candidate source electric/force state changed after preparation");
}
void ElectrosprayMeshCandidate::commit(QuadtreeMesh& mesh, std::unique_ptr<HydroFittedDomain>& domain,
    std::unique_ptr<ElectrosprayElectricState>& binding, std::unique_ptr<ElectrosprayInterfaceState>& surface,
    std::unique_ptr<Hydrodynamics>& hydro, std::unique_ptr<ElectrosprayConductingState>& conducting,
    ElectrosprayMaterialState& materials)
{
    require(domain && binding && surface && hydro && conducting,"mesh commit requires every source owner");
    validateSource({mesh,*domain,*binding,*surface,*hydro,*conducting,materials});
    auto& d=*data_;
    // Core checks background compatibility before changing storage. Every
    // following exchange is noexcept, with no solver, callback or allocation.
    static_assert(std::is_nothrow_swappable_v<ElectrosprayMaterialState>);
    mesh.swapStorage(*d.mesh);
    domain.swap(d.domain);binding.swap(d.binding);hydro.swap(d.hydro);
    surface.swap(d.interfaceState);conducting.swap(d.conducting);std::swap(materials,d.materials);
    d.pending=false;
}

std::unique_ptr<ElectrosprayMeshCandidate> prepareElectrosprayMeshAdaptation(
    ElectrosprayMeshAdaptation& policy, const ElectrosprayAcceptedStateView& accepted,
    const ElectrosprayModuleInputs& inputs,
    const UnitsContext& units,
    const ElectrosprayMeshObserver& observer)
{
    const auto& mesh=accepted.mesh;const auto& domain=accepted.domain;const auto& binding=accepted.binding;
    const auto& interfaceState=accepted.interfaceState;const auto& hydro=accepted.hydro;const auto& conducting=accepted.conducting;
    const auto& flowConfig = inputs.flow;
    const auto& electricConfig = inputs.electric;
    const auto& liquid = inputs.materials.phases.at(flowConfig.liquidPhase);
    const auto& gas = inputs.materials.phases.at(flowConfig.gasPhase);
    const double surfaceTensionSI = inputs.materials.surfaceTensionSI.at(std::minmax(flowConfig.liquidPhase, flowConfig.gasPhase));
    const double electricLinearTolerance = inputs.systems.at("ElectricPotential").relativeTolerance;
    const auto started = std::chrono::steady_clock::now();
    domain.validateMeshIdentity(mesh);
    conducting.electrical.validateMesh(mesh);
    const auto progress = hydro.state().progress;
    const double timeAllowance = 64.0 * std::numeric_limits<double>::epsilon()
        * std::max({std::abs(progress.acceptedEndTime), std::abs(conducting.electrical.time()), inputs.controls.dt});
    require(!progress.stopped && !progress.pendingAcceptance && progress.completedSteps > 0
        && progress.acceptedSteps == progress.completedSteps
        && progress.acceptedSteps == conducting.electrical.completedStep()
        && std::abs(progress.acceptedEndTime - conducting.electrical.time()) <= timeAllowance,
        "AMR requires matching accepted Hydro and electric step/time levels");
    require(conducting.forces.total.source() && conducting.forces.total.source()->sourceTime()==conducting.electrical.time(),
        "AMR requires matching accepted electric and force time levels");
    require(!domain.uniformLeaves() && interfaceState.surface.geometry() && interfaceState.surface.geometry()->usesNativeLeaves()
        && std::holds_alternative<Patch>(interfaceState.surface.scalarPatch()) && interfaceState.surface.panels()
        && !interfaceState.surface.pinnedCoordinates().empty() && conducting.electrical.charge() && conducting.electrical.electric().report.success,
        "runtime AMR requires an accepted native-Q4 conducting state");
    const auto& acceptedScalar = interfaceState.surface.geometry();
    const auto interfaceProgress = interfaceState.surface.progress();
    require(interfaceProgress.completedSteps == conducting.electrical.completedStep()
        && interfaceProgress.endTime == conducting.electrical.time(),
        "AMR requires matching accepted Interface and electric step/time levels");
    const auto& acceptedMeasure = interfaceState.surface.phase().volumeMeasure;
    require(conducting.electrical.electric().potential.source() == acceptedScalar
        && conducting.electrical.charge()->geometry().usesSource(*acceptedScalar)
        && conducting.forces.total.source()->interfaceGeometry() == acceptedScalar,
        "AMR requires matching accepted Interface, electric and force geometry");
    require(policy.everySteps > 0 && policy.maxSweeps > 0 && policy.maxSweeps <= 32,
        "invalid runtime adaptation schedule");
    require(binding.phaseByCell.size() == static_cast<std::size_t>(mesh.cellIndexNum()), "stale source binding");
    std::vector<Vector2D> line;
    for (const auto& panel : interfaceState.surface.panels()->panels())
    {
        if (line.empty()) line.push_back(panel.begin);
        for (const auto& node : panel.quadratureNodes) line.push_back(node.frame.point);
        line.push_back(panel.end);
    }
    policy.refiner.replaceBoundaryPolylines(inputs.interfaceConfig.movingBoundaryId, {line});
    require(std::isfinite(flowConfig.capillaryHalfWidthCells) && flowConfig.capillaryHalfWidthCells > 0
        && flowConfig.capillaryHalfWidthCells < std::numeric_limits<int>::max() - 2,
        "invalid interface-force stencil protection width");
    const int halo = static_cast<int>(std::ceil(flowConfig.capillaryHalfWidthCells)) + 2;
    std::vector<std::uint8_t> frozen(mesh.cellIndexNum(), 0);
    InterfaceMeshTransfer::protectScalarSupport(mesh, interfaceState.surface, domain.geometry(), frozen);
    protectHydroFittedBoundarySupport(mesh, domain, frozen);
    mesh.expandCellMask(frozen, halo);
    // Shared scalar constraints are part of the candidate's geometric support,
    // not a reason to first construct and reject a less-protected mesh. Keep
    // the ORIGINAL scalar acceptance mask below; only topology requests use
    // the existing native carrier/vertex-incidence closure.
    const auto frozenTopology = ActiveQuarticHermiteTransfer::topologyProtection(mesh, frozen);
    auto sourceIndicators = policy.flowIndicators;
    if (!sourceIndicators.empty())
    {
        const auto velocity = hydro.fittedVelocity(mesh, domain, HydroVelocityTimeLevel::Current,
            HydroFittedVelocitySampler::Reconstruction::RZFluxCompatible);
        auto fields = evaluateHydroFittedRefinement(mesh, domain, velocity);
        for (auto& field : sourceIndicators)
            field.cellValue = field.name == "Velocity" ? std::move(fields.speed) : std::move(fields.velocityGradient);
    }
    QuadtreeMeshRefiner::CurrentSource refinementSource{mesh, sourceIndicators, frozenTopology, {}};
    if (policy.currentPhases)
        refinementSource.phaseAt = [&](Vector2D point, int source) {
            return (domain.cellRowByIndex()[source] < 0 ? electricConfig.conductorPhaseId
                : interfaceState.surface.geometry()->value(point) <= 0 ? liquid.phaseId : gas.phaseId) - 1;
        };
    // Core projects the same accepted indicators and protection on every
    // sweep, and decides no-op before a candidate clones physical state.
    CellField<double> levels(mesh, 0.0), indicator(mesh, 0.0);
    if (!policy.refiner.buildMappedTargetLevels(mesh, refinementSource, levels, &indicator)) return nullptr;

    // Only a real request may stage the full mesh and physical histories.
    // A request can still be vetoed by the existing adapt() balancing rules;
    // keep the post-adapt actual-topology check and all acceptance gates.
    auto result = std::unique_ptr<ElectrosprayMeshCandidate>(new ElectrosprayMeshCandidate(accepted));
    auto& staged = *result->data_;
    staged.mesh = mesh.clone();
    auto& target = *staged.mesh;
    auto& report = staged.report;
    report.oldCells = mesh.cellIndexNum();
    report.protectedCells = static_cast<int>(std::count(frozenTopology.begin(), frozenTopology.end(), 1));
    HydroFittedMeshTransfer flow(mesh, hydro);
    auto potential=conducting.electrical.electric().potential.copyForMeshTransfer();
    std::vector<VertField<double>> unusedLevelSets;
    std::vector<VertField<double>*> vertices;
    std::vector<HEdgeField<double>*> horizontal;
    std::vector<VEdgeField<double>*> vertical;
    // clone() preserves native label/index ordering. The first request fields
    // already own their values/labels and can be used without recalculation.
    std::vector<CellField<double>*> cells{&levels, &potential};
    flow.appendTransferFields(vertices, horizontal, vertical, cells);
    for (int sweep = 0; sweep < policy.maxSweeps; ++sweep)
    {
        if (sweep > 0 && !policy.refiner.buildMappedTargetLevels(target, refinementSource, levels, &indicator)) break;
        target.adapt(unusedLevelSets, indicator, vertices, horizontal, vertical, cells, &levels);
        ++report.sweeps;
    }
    for (int i = 0; i < mesh.cellIndexNum(); ++i)
    {
        const auto& old = mesh.cells(mesh.cellIndexToLabel[i]);
        const int level = target.cells(target.cellIndexToLabel[sourceIndex(target, old.centerCoord())]).level();
        report.refinedSourceCells += level > old.level();
        report.coarsenedSourceCells += level < old.level();
    }
    if (!report.refinedSourceCells && !report.coarsenedSourceCells) return nullptr;
    report.newCells = target.cellIndexNum();
    if (observer) observer(ElectrosprayMeshStage::Grid, report.newCells);

    std::vector<int> phases(target.cellIndexNum());
    for (int i = 0; i < target.cellIndexNum(); ++i)
        phases[i] = binding.phaseByCell[sourceIndex(mesh, target.cells(target.cellIndexToLabel[i]).centerCoord())];
    staged.binding = std::make_unique<ElectrosprayElectricState>(bindElectrosprayElectrostatics(
        target, policy.definition, electricConfig, liquid, gas, units, phases, &potential));
    staged.domain = std::make_unique<Domain>(buildElectrosprayHydroDomain(
        target, *staged.binding, policy.definition, flowConfig,
        liquid, gas, electricConfig.conductorPhaseId, units));
    require(!staged.domain->uniformLeaves(), "runtime native-Q4 workflow cannot switch to a uniform-Q3 path");
    auto protectedScalar = frozen;
    for (int i = 0; i < mesh.cellIndexNum(); ++i)
        if (domain.cellRowByIndex()[i] < 0) protectedScalar[i] = 0;
    auto state = std::make_unique<ElectrosprayInterfaceState>(target);
    state->freeLineIds = interfaceState.freeLineIds;
    state->fixedLineIds = interfaceState.fixedLineIds;
    state->sourceLevelSetIndex = interfaceState.sourceLevelSetIndex;
    const auto transferred = InterfaceMeshTransfer::preparePinnedToAxis(mesh, interfaceState.surface,
        target, staged.domain->geometry(), protectedScalar, state->surface);
    auto charge = Quantum::Electricity::prepareFittedChargeMeshTransfer(mesh, conducting.electrical,
        target, state->surface.geometry(), transferred.sourceToTargetCell, state->surface.panels());
    staged.interfaceState = std::move(state);
    auto& surface = *staged.interfaceState;
    staged.materials = prepareElectrosprayMaterials(target, *staged.domain, *surface.surface.geometry(),
        flowConfig, electricConfig, liquid, gas, units, *staged.binding);
    const ElectrosprayHydroScales scales(liquid.densitySI, units);
    updateElectrosprayHydroMaterials(target, *staged.domain,
        HydroMaterialFields{flow.density, flow.viscosity, flow.vertexDensity, flow.vertexViscosity},
        staged.materials.phaseByCell, liquid, gas, scales);
    if (observer) observer(ElectrosprayMeshStage::InterfaceAndCharge, report.newCells);
    auto transferredFlow = flow.prepare(target, *staged.domain, electrosprayHydroNumerics(inputs));
    staged.hydro = std::move(transferredFlow.solver);
    report.hydro = transferredFlow.report;
    if (observer) observer(ElectrosprayMeshStage::Flow, report.newCells);

    staged.conducting = std::make_unique<ElectrosprayConductingState>(target);
    auto& next = *staged.conducting;
    const auto& oldQ=conducting.electrical.charge()->transport().panelCharge();
    const auto& newQ=charge.charge()->transport().panelCharge();
    report.chargeBefore=std::accumulate(oldQ.begin(),oldQ.end(),0.0);
    report.chargeAfter=std::accumulate(newQ.begin(),newQ.end(),0.0);
    const Quantum::Electricity::FittedChargeProblem problem{staged.binding->dielectric.graph(),staged.binding->dielectric.boundaries(),
        bindElectrosprayConductingParameters(electricConfig,liquid,gas,electricLinearTolerance,units)};
    next.electrical=Quantum::Electricity::solveFittedChargeMeshTransfer(target,problem,charge,
        surface.surface.geometry(),potential);
    const auto source=next.electrical.source();
    next.forces = buildElectrosprayConductingForces(target, *staged.domain, surface.surface.geometry(),
        *surface.surface.panels(), source, next.electrical.electric(), liquid, gas, surfaceTensionSI, electricConfig, flowConfig.capillaryHalfWidthCells, units,
        conducting.forces.total.source()->sourceTime());
    if (observer) observer(ElectrosprayMeshStage::ElectricAndForces, report.newCells);

    const auto velocity = staged.hydro->fittedVelocity(target, *staged.domain, HydroVelocityTimeLevel::Current,
        HydroFittedVelocitySampler::Reconstruction::RZFluxCompatible);
    surface.surface.initializePhase(target, HydroInterfaceTransport::domainView(*staged.domain), HydroInterfaceTransport::normalFlux(velocity));
    report.volumeBefore = std::visit([](const auto& m) { return m.summedNegativeVolume; }, acceptedMeasure);
    report.volumeAfter = std::visit([](const auto& m) { return m.summedNegativeVolume; }, surface.surface.phase().volumeMeasure);
    // Core already certifies identical zero sets. This separate bound covers
    // floating summation of the same volume over different leaf partitions;
    // it is not a new physical volume-correction tolerance or a volume reset.
    report.volumeRoundoffBound = 8.0L * std::numeric_limits<double>::epsilon()
        * (report.oldCells + report.newCells + 2) * std::abs(report.volumeBefore);
    require(std::isfinite(report.volumeBefore) && report.volumeBefore > 0
        && std::isfinite(report.volumeAfter) && report.volumeAfter > 0
        && std::abs(report.volumeAfter - report.volumeBefore) <= report.volumeRoundoffBound,
        "zero-time mesh transfer changed liquid volume beyond summation roundoff");
    report.liquidFluxBefore = hydro.fittedPhaseOutwardFlux(mesh, domain, *interfaceState.surface.geometry());
    report.liquidFluxAfter = staged.hydro->fittedPhaseOutwardFlux(target, *staged.domain, *surface.surface.geometry());
    const auto& budget = surface.surface.phase().phaseBudget;
    require(std::abs(budget.assemblyClosure) <= flowConfig.phaseFluxBalance * budget.absoluteFaceFlux
        && std::abs(budget.boundaryOutward - report.liquidFluxAfter) <= flowConfig.phaseFluxBalance
            * (std::abs(budget.boundaryOutward) + std::abs(report.liquidFluxAfter)),
        "target phase flux ledger failed original assembly/pressure-port closure");
    report.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
    if (observer) observer(ElectrosprayMeshStage::Volumes, report.newCells);
    return result;
}
