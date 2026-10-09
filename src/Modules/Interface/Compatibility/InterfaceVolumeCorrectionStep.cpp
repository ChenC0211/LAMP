#include "Modules/Interface/Compatibility/InterfaceVolumeCorrection.h"
#include "Modules/Interface/Compatibility/Detail/InterfaceCorrectionRegularization.h"
#include "Modules/Interface/Compatibility/InterfaceTopology.h"
#include "Modules/Interface/Compatibility/InterfaceMotion.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Geometry/InterfaceGeometryResolution.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <tuple>
#include <type_traits>

namespace InterfaceVolumeCorrection
{
namespace
{
using Patch = Quantum::Geometry::ActiveHermitePatch;
enum class SourceLayout { UniformLeaves, NativeLeaves };
void require(bool condition, const std::string& reason)
{
    if (!condition) throw std::runtime_error(reason);
}
struct Sum
{
    long double value = 0.0L, compensation = 0.0L;
    void add(long double term)
    {
        require(std::isfinite(term), "nonfinite residual summand");
        const long double increment = term-compensation, next = value+increment;
        compensation = (next-value)-increment;
        value = next;
        require(std::isfinite(value) && std::isfinite(compensation), "nonfinite residual accumulation");
    }
};
template<class State> struct TrialState
{
    State state;
    std::unique_ptr<InterfacePhaseGeometry::MappedFacePartitions> partitions;
    using TimeResult=std::conditional_t<std::is_same_v<State,QuarticCorrection>,
        InterfaceTemporalPhaseFlux::Result,InterfaceTemporalPhaseFlux::CubicResult>;
    std::unique_ptr<TimeResult> timeFlux;
    Residual residual;
    std::vector<long double> defect;
    std::vector<int> pressureSigns;
    bool targetsInBounds = true;
    long double targetViolation = 0.0L; // Sum of required range corrections; no target is clipped.
    std::string targetFailure;
};

std::string describeCell(int cell, long double oldVolume, long double volume,
    long double transport, long double capacity)
{
    std::ostringstream text;
    text << std::setprecision(17) << "cell=" << cell << ", R=" << (volume-oldVolume)+transport
        << ", Vold=" << oldVolume << ", Vtrial=" << volume << ", transport=" << transport
        << ", capacity=" << capacity;
    return text.str();
}

template<class ScalarPatch>
std::vector<std::uint8_t> validatePatch(const QuadtreeMesh& mesh, const InterfaceTransportDomain& transport,
    const ScalarPatch& patch, const std::vector<Vector2D>& pins, const std::vector<std::uint8_t>& fixed,
    SourceLayout layout)
{
    const auto& domain = transport.geometry();
    require(patch.nodes.size() == static_cast<std::size_t>(mesh.vertIndexNum())
        && fixed.size() == patch.nodes.size()
        && patch.activeCells.size() == static_cast<std::size_t>(mesh.cellIndexNum())
        && domain.cellRowByIndex().size() == patch.activeCells.size(), "native patch/mask indexing mismatch");
    std::vector<std::uint8_t> active(patch.nodes.size(),0);
    for (int cell=0; cell<mesh.cellIndexNum(); ++cell)
    {
        require(patch.activeCells[cell] <= 1
            && static_cast<bool>(patch.activeCells[cell]) == (domain.cellRowByIndex()[cell] >= 0),
            "correction cannot change the actual solid/fluid cell mask");
        if (!patch.activeCells[cell]) continue;
        for (int label : mesh.cells(mesh.cellIndexToLabel[cell]).cellToVert)
        {
            const int vertex = mesh.verts(label).index();
            require(vertex >= 0 && vertex < mesh.vertIndexNum(), "active cell has invalid native vertex");
            active[vertex] = 1;
        }
    }
    for (std::size_t vertex=0; vertex<fixed.size(); ++vertex)
    {
        require(fixed[vertex] <= 1 && (active[vertex] || fixed[vertex]), "invalid fixed-value mask or free inactive vertex");
        if (!active[vertex]) continue; // Inactive payload is never scalar data.
        if (layout == SourceLayout::NativeLeaves
            && mesh.verts(mesh.vertIndexToLabel[vertex]).supCell >= 0)
        {
            // The scalar owns all four coarse-edge constraints. These are
            // not independent columns, even when the caller mask says free.
            require(!fixed[vertex], "a dependent hanging value cannot be prescribed by the fixed-value mask");
            active[vertex] = 0;
            continue;
        }
        const auto& node = patch.nodes[vertex];
        require(std::isfinite(node.value) && std::isfinite(node.dx) && std::isfinite(node.dy)
            && std::isfinite(node.dxy), "nonfinite active Hermite DOF");
    }
    for (const auto& face : domain.faces(mesh))
    {
        if (transport.boundary(face.nativeRow) != InterfaceTransportDomain::Boundary::StationaryWall) continue;
        if (face.negativeCell < 0 && face.positiveCell < 0) continue;
        const auto& edge = face.horizontal ? mesh.hEdges(face.edgeLabel) : mesh.vEdges(face.edgeLabel);
        for (int label : edge.edgeToVert)
            require(fixed[mesh.verts(label).index()] == 1, "actual no-slip wall value is not fixed");
    }
    for (const auto pin : pins)
    {
        require(std::isfinite(pin.x) && std::isfinite(pin.y), "nonfinite pinned coordinate");
        int found = -1;
        for (int vertex=0; vertex<mesh.vertIndexNum(); ++vertex)
        {
            if (!active[vertex]) continue;
            const auto point = mesh.verts(mesh.vertIndexToLabel[vertex]).coord();
            if (point.x != pin.x || point.y != pin.y) continue;
            require(found < 0, "pin has multiple actual native vertex identities");
            found = vertex;
        }
        require(found >= 0 && fixed[found] && patch.nodes[found].value == 0.0,
            "pin must be an actual fixed zero vertex; snapping is unsupported");
    }
    for (const auto& sector : patch.cornerSectors)
    {
        require(sector.vertexIndex >= 0 && sector.vertexIndex < mesh.vertIndexNum()
            && active[sector.vertexIndex] && fixed[sector.vertexIndex]
            && sector.dof.value == 0.0 && patch.nodes[sector.vertexIndex].value == 0.0,
            "corner sector must retain its fixed exact-zero value");
    }
    // Axis values deliberately remain under the caller's mask. Snapshot
    // enforces dx=dxy=0 there, while the axial interface endpoint can move.
    return active;
}


template<class ScalarPatch> auto& quarticDatum(ScalarPatch& patch,const QuarticDof& id)
{
    using Kind=QuarticDof::Kind;
    if(id.kind==Kind::CellCenter) return patch.cellCenters.at(id.index);
    if(id.kind==Kind::EdgeMidpoint)
    {
        auto& edge=patch.edgeMidpoints.at(id.index);
        return id.component==0?edge.value:edge.normalDerivative;
    }
    auto& node=patch.nodes.at(id.index);
    return id.component==0?node.value:id.component==1?node.dx:node.dy;
}

// Eligibility only; the scalar constructor remains the authority for traces,
// axis data, topology and the full-degree sector certificate at every trial.
void validateQuarticData(const QuadtreeMesh& mesh,const InterfaceTransportDomain& transport,
    const Quantum::Geometry::ActiveQuarticHermitePatch& patch,
    const std::vector<std::uint8_t>& active,const std::vector<std::uint8_t>& fixed,
    const std::vector<QuarticDof>& freeData)
{
    const auto& domain = transport.geometry();
    using Kind=QuarticDof::Kind;
    using EdgeKey=std::array<int,2>;
    std::set<EdgeKey> dependentSides;
    std::vector<std::array<Vector2D,2>> walls;
    for(const auto& face:domain.faces(mesh))
    {
        if(face.negativeCell<0 && face.positiveCell<0) continue;
        const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
        const auto& a=mesh.verts(edge.edgeToVert[0]);const auto& b=mesh.verts(edge.edgeToVert[1]);
        if(face.negativeCell>=0 && face.positiveCell>=0
            && mesh.cells(mesh.cellIndexToLabel[face.negativeCell]).level()
                !=mesh.cells(mesh.cellIndexToLabel[face.positiveCell]).level())
            dependentSides.insert({{a.index(),b.index()}});
        if(transport.boundary(face.nativeRow)==InterfaceTransportDomain::Boundary::StationaryWall)
            walls.push_back({{a.coord(),b.coord()}});
    }
    std::set<std::tuple<int,int,int>> declared;
    for(const auto& id:freeData)
    {
        require(id.index>=0 && std::isfinite(id.scale) && id.scale>0,"invalid quartic free-data identity/scale");
        require(declared.emplace(static_cast<int>(id.kind),id.index,id.component).second,
            "duplicate quartic free-data identity");
        if(id.kind==Kind::Vertex)
        {
            require(static_cast<std::size_t>(id.index)<active.size() && active[id.index]
                && id.component>=0 && id.component<=2,"quartic free vertex is inactive, hanging or not a value/first derivative");
            require(id.component!=0 || !fixed[id.index],"quartic correction cannot free a fixed/pinned vertex value");
            require(id.component!=1 || !domain.axisymmetric() || mesh.verts(mesh.vertIndexToLabel[id.index]).coord().x!=0,
                "quartic correction cannot free the radial-axis normal derivative");
        }
        else if(id.kind==Kind::EdgeMidpoint)
        {
            require(static_cast<std::size_t>(id.index)<patch.edgeMidpoints.size()
                && (id.component==0 || id.component==1),"invalid free quartic edge data");
            const auto& input=patch.edgeMidpoints[id.index];
            require(dependentSides.count(input.vertices)==0,"dependent fine-side midpoint cannot be a correction unknown");
            for(int vertex:input.vertices)require(vertex>=0 && vertex<mesh.vertIndexNum(),"invalid quartic edge endpoint");
            const auto a=mesh.verts(mesh.vertIndexToLabel[input.vertices[0]]).coord();
            const auto b=mesh.verts(mesh.vertIndexToLabel[input.vertices[1]]).coord();
            require(id.component!=1 || !domain.axisymmetric() || a.x!=0 || b.x!=0,"quartic correction cannot free the axis-edge normal derivative");
            if(id.component==0)for(const auto& wall:walls)
            {
                const auto c=wall[0],d=wall[1];
                const bool horizontal=a.y==b.y;
                const bool overlaps=horizontal
                    ? c.y==d.y && a.y==c.y && std::max(std::min(a.x,b.x),std::min(c.x,d.x))<std::min(std::max(a.x,b.x),std::max(c.x,d.x))
                    : c.x==d.x && a.x==c.x && std::max(std::min(a.y,b.y),std::min(c.y,d.y))<std::min(std::max(a.y,b.y),std::max(c.y,d.y));
                require(!overlaps,"quartic correction cannot free a no-slip wall midpoint value");
            }
        }
        else
            require(id.kind==Kind::CellCenter && static_cast<std::size_t>(id.index)<patch.activeCells.size()
                && patch.activeCells[id.index] && id.component==0,"invalid free quartic center data");
        require(std::isfinite(quarticDatum(patch,id)),"nonfinite independent quartic correction data");
    }
}

bool balanceConverged(const Residual& residual, const IterationOptions& options)
{
    const bool relative = options.absoluteVolumeTolerance > 0.0L
        ? residual.sumAbsolute <= options.absoluteVolumeTolerance
            + options.relativeBalanceTolerance*residual.transportedMagnitude
        : residual.transportedMagnitude == 0.0L ? residual.sumAbsolute == 0.0L
            : residual.relativeBalance <= options.relativeBalanceTolerance;
    return relative && residual.maxRelativeCell <= options.maxCellRelativeTolerance;
}

std::string describe(const Report& report)
{
    std::ostringstream text;
    text << std::setprecision(17) << "iterations=" << report.iterations << ", backtracks=" << report.backtracks
        << ", rows=" << report.rows << ", columns=" << report.columns << ", usedSuperLU=" << report.usedSuperLU
        << ", maxCorrection=" << report.maxCorrection
        << ", maxGradientCorrection=" << report.maxGradientCorrection
        << ", maxLinearGradientCorrection=" << report.maxLinearGradientCorrection
        << ", maxLinearCorrection=" << report.maxLinearCorrection
        << ", scaledDefectCell=" << report.largestScaledDefectCell
        << ", maxScaledDefect=" << report.largestScaledDefect
        << ", freeRowNorm=" << report.correspondingRowNorm
        << ", initial=(sumAbs=" << report.initial.sumAbsolute << ", transport=" << report.initial.transportedMagnitude
        << ", relative=" << report.initial.relativeBalance << ", maxCell=" << report.initial.maxRelativeCell << ')'
        << ", final=(sumAbs=" << report.final.sumAbsolute << ", transport=" << report.final.transportedMagnitude
        << ", relative=" << report.final.relativeBalance << ", maxCell=" << report.final.maxRelativeCell << ')'
        << ", worst={" << describeCell(report.final.worstCell,report.final.worstOldVolume,
            report.final.worstVolume,report.final.worstTransport,report.final.worstCapacity) << '}';
    return text.str();
}
template<bool Quartic>
auto correctImpl(const QuadtreeMesh& mesh, const InterfaceTransportDomain& transport,
    const std::conditional_t<Quartic,InterfaceCellMeasure::QuarticMeasure,InterfaceCellMeasure::Measure>& oldMeasure,
    const InterfacePhaseFlux::FluxBudget& oldFlux,
    const InterfaceNormalFluxView& currentVelocity, double dt,
    const std::conditional_t<Quartic,Quantum::Geometry::ActiveQuarticHermitePatch,Patch>& advectedPatch,
    const std::vector<Vector2D>& pinnedCoordinates, const std::vector<std::uint8_t>& fixedValue,
    const std::conditional_t<Quartic,QuarticOptions,Options>& options, SourceLayout layout,
    bool completeComponents = false,const InterfaceState* previousState = nullptr,
    const InterfaceNormalFluxView* previousVelocity = nullptr,
    const InterfaceTemporalPhaseFlux::Options* timeFluxOptions = nullptr,
    const InterfaceTemporalPhaseFlux::Trajectory* trajectory = nullptr,
    const TemporalInput<std::conditional_t<Quartic,Quantum::Geometry::ActiveQuarticHermitePatch,Patch>>* temporal = nullptr)
{
    const auto& domain = transport.geometry();
    using LocalPatch=std::conditional_t<Quartic,Quantum::Geometry::ActiveQuarticHermitePatch,Patch>;
    using Trial=TrialState<std::conditional_t<Quartic,QuarticCorrection,Correction>>;
    Report report;
    const auto notify = [&](Progress::Stage stage, int backtrack = 0, double damping = 1.0,
                            const std::string& detail = std::string()) {
        if (options.progress) options.progress({stage,report,backtrack,damping,detail});
    };
    try
    {
        domain.validateMeshIdentity(mesh);
        require(!completeComponents || (previousState && previousVelocity && timeFluxOptions),
            "complete correction requires its old state, velocity trace and time quadrature options");
        require(!completeComponents || !temporal,"two independent time-flux authorities supplied");
        if constexpr(Quartic) {
            require(std::isfinite(options.geometricRegularization) && options.geometricRegularization>=0,
                "geometric regularization strength must be finite and nonnegative");
            require(completeComponents || options.geometricRegularization==0,
                "geometric regularization requires the complete component workflow");
        }
        require(domain.axisymmetric() || completeComponents, "legacy volume correction requires the fitted RZ domain");
        require(std::isfinite(dt) && dt > 0.0, "positive finite correction time step is required");
        require(std::isfinite(options.linearTolerance) && options.linearTolerance > 0.0
            && std::isfinite(options.relativeBalanceTolerance) && options.relativeBalanceTolerance > 0.0L
            && std::isfinite(options.maxCellRelativeTolerance) && options.maxCellRelativeTolerance > 0.0L
            && std::isfinite(options.absoluteVolumeTolerance) && options.absoluteVolumeTolerance >= 0.0L
            && options.maxIterations > 0 && options.maxBacktracks >= 0,
            "invalid correction tolerances or iteration limits");
        const auto active = validatePatch(mesh,transport,advectedPatch,pinnedCoordinates,fixedValue,layout);
        if (layout == SourceLayout::NativeLeaves)
            for (const auto& arc : oldMeasure.arcs)
                require(arc.nativeSourceCell >= 0 && !arc.independentSensitivity.empty(),
                    "native correction requires an old native-leaf volume measure");
        if constexpr(Quartic)
        {
            require(layout==SourceLayout::NativeLeaves,"quartic correction requires native leaves");
            validateQuarticData(mesh,transport,advectedPatch,active,fixedValue,options.freeData);
        }
        else
        {
            std::set<std::pair<int,unsigned>> gradientIdentities;
            for (const auto& dof : options.freeGradients)
            {
                require(dof.vertexIndex >= 0 && static_cast<std::size_t>(dof.vertexIndex) < active.size()
                    && active[dof.vertexIndex] && (dof.component == 1 || dof.component == 2)
                    && std::isfinite(dof.scale) && dof.scale > 0.0,
                    "appended gradient must have a valid active native identity, component and scale");
                require(gradientIdentities.emplace(dof.vertexIndex,dof.component).second,
                    "duplicate appended first-derivative identity");
            }
        }
        const auto cells = static_cast<std::size_t>(mesh.cellIndexNum());
        require(oldMeasure.negativeVolumeByCell.size() == cells
            && oldMeasure.fixedBoundary.capacityByCell.size() == cells && oldFlux.outwardFluxByCell.size() == cells,
            "old native volume/flux indexing mismatch");
        for (int cell : domain.cellIndexByRow())
        {
            const auto capacity = oldMeasure.fixedBoundary.capacityByCell[cell];
            const auto volume = oldMeasure.negativeVolumeByCell[cell];
            require(std::isfinite(capacity) && capacity > 0.0L && std::isfinite(volume)
                && volume >= 0.0L && volume <= capacity && std::isfinite(oldFlux.outwardFluxByCell[cell]),
                "old active volume is nonfinite/out of bounds or old phase flux is nonfinite at cell="+std::to_string(cell));
        }
        std::set<int> pressureCellSet;
        for (const auto& face : domain.faces(mesh))
            if (transport.boundary(face.nativeRow) == InterfaceTransportDomain::Boundary::SinglePhasePort)
            {
                const int owner = face.negativeCell >= 0 ? face.negativeCell : face.positiveCell;
                require(owner >= 0 && owner < mesh.cellIndexNum(), "pressure face lacks an active owner");
                pressureCellSet.insert(owner);
            }
        const std::vector<int> pressureCells(pressureCellSet.begin(),pressureCellSet.end());
        const long double halfDt = 0.5L*dt;
        require(std::isfinite(halfDt) && halfDt > 0.0L, "unrepresentable half time step");
        const auto evaluate = [&](const LocalPatch& patch) {
            using Clock = std::chrono::steady_clock;
            const auto started = Clock::now();
            Trial trial;
            trial.state.patch = patch;
            if (completeComponents)
            {
                trial.state.completeGeometry = std::make_unique<InterfaceState>(mesh);
                InterfaceGeometryReconstruction::initialize(mesh,patch,domain.axisymmetric(),pinnedCoordinates,
                    *trial.state.completeGeometry);
                trial.state.geometry = trial.state.completeGeometry->geometry();
            }
            else if (layout == SourceLayout::NativeLeaves)
            {
                trial.state.geometry = std::make_shared<LevelSetHermiteGeometry::Snapshot>(
                    LevelSetHermiteGeometry::Snapshot::fromNativeLeaves(mesh,patch,true));
                trial.state.nativeChain = std::make_unique<LevelSetInterfacePanels::NativeGraphChain>(
                    LevelSetInterfacePanels::discoverNativeGraphs(mesh,trial.state.geometry,true,pinnedCoordinates));
            }
            else
            {
                if constexpr(!Quartic)
                {
                    trial.state.geometry = std::make_shared<LevelSetHermiteGeometry::Snapshot>(mesh,patch,true,
                        CartesianGridScope::UniformLeaves);
                    trial.state.panels = std::make_unique<LevelSetInterfacePanels>(mesh,trial.state.geometry,true,
                        pinnedCoordinates,CartesianGridScope::UniformLeaves);
                }
                else require(false,"quartic correction has no uniform-leaf fallback");
            }
            const auto geometryReady = Clock::now();
            report.timing.geometry += std::chrono::duration<double>(geometryReady-started).count();
            InterfacePhaseGeometry::MappedFacePartitions partitions(mesh,domain,trial.state.geometry);
            const auto partitionsReady = Clock::now();
            report.timing.facePartitions += std::chrono::duration<double>(partitionsReady-geometryReady).count();
            if constexpr(Quartic)
                trial.state.measure=completeComponents
                    ? InterfaceCellMeasure::evaluateQuarticRefined(mesh,domain,partitions,
                        *trial.state.completeGeometry,options.maxRootRefinements)
                    : InterfaceCellMeasure::evaluateQuarticRefined(mesh,domain,partitions,
                        *trial.state.nativeChain,options.maxRootRefinements);
            else
                trial.state.measure = layout == SourceLayout::NativeLeaves
                    ? InterfaceCellMeasure::evaluateRefined(mesh,domain,partitions,
                        *trial.state.nativeChain,options.maxRootRefinements)
                    : InterfaceCellMeasure::evaluateRefined(mesh,domain,partitions,
                        *trial.state.panels,options.maxRootRefinements);
            const auto volumeReady = Clock::now();
            report.timing.volume += std::chrono::duration<double>(volumeReady-partitionsReady).count();
            // Refinement invalidates earlier integrals: flux is formed only
            // after the complete measure has resolved this common partition.
            trial.state.flux = InterfacePhaseFlux::phaseFluxBudget(mesh,transport,partitions,currentVelocity);
            if constexpr(Quartic) if (completeComponents)
                trial.timeFlux=std::make_unique<InterfaceTemporalPhaseFlux::Result>(InterfaceTemporalPhaseFlux::integrate(
                    mesh,transport,std::get<Quantum::Geometry::ActiveQuarticHermitePatch>(previousState->scalarPatch()),patch,
                    previousState->geometry(),trial.state.geometry,*previousVelocity,currentVelocity,dt,*timeFluxOptions,trajectory));
            if (temporal) trial.timeFlux=std::make_unique<typename Trial::TimeResult>(
                InterfaceTemporalPhaseFlux::integrate(mesh,transport,temporal->previousPatch,patch,
                    temporal->previousGeometry,trial.state.geometry,temporal->previousVelocity,currentVelocity,
                    dt,temporal->options,temporal->trajectory));
            const auto fluxReady = Clock::now();
            report.timing.phaseFlux += std::chrono::duration<double>(fluxReady-volumeReady).count();
            trial.defect.assign(cells,std::numeric_limits<long double>::quiet_NaN());
            Sum absolute,transport,signedResidual,targetViolation;
            for (int cell : domain.cellIndexByRow())
            {
                const auto capacity = trial.state.measure.fixedBoundary.capacityByCell[cell];
                const auto volume = trial.state.measure.negativeVolumeByCell[cell];
                const auto oldVolume = oldMeasure.negativeVolumeByCell[cell];
                const auto oldQ = oldFlux.outwardFluxByCell[cell], newQ = trial.state.flux.outwardFluxByCell[cell];
                require(capacity == oldMeasure.fixedBoundary.capacityByCell[cell], "actual CV capacity changed during correction");
                const long double change = trial.timeFlux ? trial.timeFlux->integral.volume.outwardFluxByCell[cell]
                    : halfDt*(oldQ+newQ);
                const long double residual = (volume-oldVolume)+change;
                const long double target = oldVolume-change;
                require(std::isfinite(volume) && std::isfinite(newQ) && std::isfinite(change)
                    && std::isfinite(residual) && std::isfinite(target),
                    "nonfinite local phase-flux balance; "+describeCell(cell,oldVolume,volume,change,capacity));
                require(volume >= 0.0L && volume <= capacity,
                    "trial raw volume out of bounds; "+describeCell(cell,oldVolume,volume,change,capacity));
                const long double violation=target<0?-target:target>capacity?target-capacity:0;
                const long double uncertainty=trial.timeFlux
                    ?trial.timeFlux->integral.integrationUncertaintyByCell[cell]
                        +2*std::numeric_limits<double>::epsilon()*(std::abs(oldVolume)+std::abs(change))
                    :0;
                trial.residual.maxTargetViolation=std::max(trial.residual.maxTargetViolation,violation);
                trial.residual.maxTargetUncertainty=std::max(trial.residual.maxTargetUncertainty,uncertainty);
                if (violation > uncertainty)
                {
                    trial.targetsInBounds = false;
                    targetViolation.add(violation-uncertainty);
                    if (trial.targetFailure.empty())
                    {
                        std::ostringstream reason;
                        reason << std::setprecision(17) << "phase-flux target out of bounds at cell=" << cell
                            << ", target=" << target << ", numericalUncertainty=" << uncertainty
                            << "; " << describeCell(cell,oldVolume,volume,change,capacity);
                        trial.targetFailure = reason.str();
                    }
                }
                trial.defect[cell] = -residual; // No target-minus-large-volume cancellation.
                absolute.add(std::abs(residual));
                signedResidual.add(residual);
                transport.add(trial.timeFlux ? trial.timeFlux->integral.absoluteTransportByCell[cell]
                    : halfDt*(std::abs(oldQ)+std::abs(newQ)));
                const long double relativeCell = std::abs(residual)/capacity;
                if (trial.residual.worstCell < 0 || relativeCell > trial.residual.maxRelativeCell)
                {
                    trial.residual.maxRelativeCell = relativeCell;
                    trial.residual.worstCell = cell;
                    trial.residual.worstResidual = residual;
                    trial.residual.worstOldVolume = oldVolume;
                    trial.residual.worstVolume = volume;
                    trial.residual.worstTransport = change;
                    trial.residual.worstCapacity = capacity;
                }
            }
            trial.residual.sumAbsolute = absolute.value;
            trial.residual.transportedMagnitude = transport.value;
            trial.residual.signedSum = signedResidual.value;
            trial.targetViolation = targetViolation.value;
            require(trial.targetsInBounds ? trial.targetViolation == 0.0L : trial.targetViolation > 0.0L,
                "target range violation lost its positive measure");
            if (!trial.targetsInBounds)
            {
                std::ostringstream reason;
                reason << std::setprecision(17) << "totalTargetViolation=" << trial.targetViolation
                    << "; " << trial.targetFailure;
                trial.targetFailure = reason.str();
            }
            // Infinity is a diagnostic for a nonzero defect with no transport,
            // not a valid numeric input or a way to pass the zero-budget gate.
            trial.residual.relativeBalance = transport.value > 0.0L ? absolute.value/transport.value
                : absolute.value == 0.0L ? 0.0L : std::numeric_limits<long double>::infinity();
            require(std::isfinite(trial.residual.maxRelativeCell), "nonfinite capacity-normalized residual");
            for (int cell : pressureCells)
            {
                const int sign = trial.state.geometry->closedCellSign(mesh.cells(mesh.cellIndexToLabel[cell]).centerCoord());
                require(sign != 0, "correction reached a non-pure pressure-port cell="+std::to_string(cell));
                trial.pressureSigns.push_back(sign);
            }
            report.timing.balanceChecks += std::chrono::duration<double>(Clock::now()-fluxReady).count();
            ++report.timing.evaluations;
            // Only the native balance Jacobian consumes cuts after evaluate.
            // Uniform/Picard callers retain the original stack lifetime and
            // do not keep a second large face partition during the next trial.
            if (layout == SourceLayout::NativeLeaves)
                trial.partitions = std::make_unique<InterfacePhaseGeometry::MappedFacePartitions>(std::move(partitions));
            return trial;
        };

        Trial current = evaluate(advectedPatch);
        const auto pressureSigns = current.pressureSigns;
        report.initial = report.final = current.residual;
        notify(Progress::Stage::PredictorReady);
        for (int iteration=0; ; ++iteration)
        {
            if (current.targetsInBounds && balanceConverged(current.residual,options))
            {
                current.state.report = report;
                if (current.timeFlux) current.state.stepFlux=std::make_unique<InterfaceTemporalPhaseFlux::Integral>(
                    std::move(current.timeFlux->integral));
                notify(Progress::Stage::Converged);
                return std::move(current.state);
            }
            require(iteration < options.maxIterations, "nonlinear correction iteration limit exhausted"
                +(current.targetsInBounds ? std::string() : "; "+current.targetFailure));
            ++report.iterations;
            const auto linearStarted = std::chrono::steady_clock::now();
            // Trials update explicit shared data only; declared one-sided
            // contact jets remain fixed. Pins/shared traces are still checked.
            const auto linear = [&] {
                if constexpr(Quartic) {
                    if (current.timeFlux) return solveQuarticBalanceDefect(current.state.measure,
                        current.timeFlux->derivative,1.0L,current.defect,options.freeData,options.linearTolerance,
                        completeComponents?geometricRegularization(mesh,*current.state.completeGeometry,options.geometricRegularization)
                            :QuarticRegularization{});
                    return solveQuarticBalanceDefect(current.state.measure,
                        InterfacePhaseFlux::quarticPhaseFluxSensitivity(mesh,transport,*current.partitions,currentVelocity,
                            Quantum::Geometry::CornerSectorVariation::HoldFixed),
                        halfDt,current.defect,options.freeData,options.linearTolerance);
                }
                else {
                if (current.timeFlux) return solveBalanceDefect(current.state.measure,
                    current.timeFlux->derivative,1.0L,current.defect,fixedValue,options.linearTolerance,options.freeGradients);
                return layout==SourceLayout::NativeLeaves
                ? solveBalanceDefect(current.state.measure,
                    InterfacePhaseFlux::phaseFluxSensitivity(mesh,transport,*current.partitions,currentVelocity,
                        Quantum::Geometry::CornerSectorVariation::HoldFixed),
                    halfDt,current.defect,fixedValue,options.linearTolerance,options.freeGradients)
                : solveDefect(current.state.measure,current.defect,fixedValue,options.linearTolerance,options.freeGradients);
                }
            }();
            report.timing.linearSolve += std::chrono::duration<double>(
                std::chrono::steady_clock::now()-linearStarted).count();
            report.rows = linear.rows; report.columns = linear.columns;
            report.largestScaledDefectCell = linear.largestScaledDefectCell;
            report.largestScaledDefect = linear.largestScaledDefect;
            report.correspondingRowNorm = linear.correspondingRowNorm;
            if constexpr(Quartic)
            {
                require(linear.deltaData.size()==options.freeData.size(),"quartic linear correction indexing mismatch");
                for(std::size_t i=0;i<linear.deltaData.size();++i)
                {
                    require(std::isfinite(linear.deltaData[i]),"nonfinite quartic linear increment");
                    auto& maximum=options.freeData[i].component==0?report.maxLinearCorrection:report.maxLinearGradientCorrection;
                    maximum=std::max(maximum,std::abs(linear.deltaData[i]));
                }
                report.usedSuperLU=report.usedSuperLU || linear.usedSuperLU;
            }
            else
            {
                for (const auto increment : linear.deltaValue)
                    report.maxLinearCorrection = std::max(report.maxLinearCorrection, std::abs(increment));
                for (const auto increment : linear.deltaGradient)
                    report.maxLinearGradientCorrection = std::max(report.maxLinearGradientCorrection,std::abs(increment));
                report.usedSuperLU = report.usedSuperLU || linear.usedSuperLU;
                require(linear.deltaValue.size() == advectedPatch.nodes.size(), "linear correction indexing mismatch");
                require(linear.deltaGradient.size() == options.freeGradients.size(), "linear gradient correction indexing mismatch");
                for (const auto increment : linear.deltaGradient)
                    require(std::isfinite(increment), "nonfinite linear first-derivative correction");
                for (std::size_t vertex=0; vertex<linear.deltaValue.size(); ++vertex)
                    require(std::isfinite(linear.deltaValue[vertex])
                        && (!fixedValue[vertex] || linear.deltaValue[vertex] == 0.0), "linear correction changed a fixed value or is nonfinite");
            }
            notify(Progress::Stage::LinearStepReady);

            bool accepted = false;
            std::string lastFailure;
            double damping = 1.0;
            for (int backtrack=0; backtrack<=options.maxBacktracks; ++backtrack)
            {
                if (backtrack != 0) { damping *= 0.5; ++report.backtracks; }
                notify(Progress::Stage::TrialStarted,backtrack,damping);
                try
                {
                    // Change only the declared DOFs. All other derivatives,
                    // sectors and masks stay unchanged; there is no Q4 refit.
                    LocalPatch patch = current.state.patch;
                    bool changed = false;
                    if constexpr(Quartic)
                    {
                        for(std::size_t i=0;i<options.freeData.size();++i)
                        {
                            auto& data=quarticDatum(patch,options.freeData[i]);
                            const double increment=damping*linear.deltaData[i],next=data+increment;
                            require(std::isfinite(increment) && std::isfinite(next),"nonfinite quartic trial data");
                            changed=changed || next!=data;
                            data=next;
                        }
                    }
                    else
                    {
                        for (std::size_t vertex=0; vertex<patch.nodes.size(); ++vertex)
                        {
                            if (!active[vertex] || fixedValue[vertex]) continue;
                            const double increment = damping*linear.deltaValue[vertex];
                            const double value = patch.nodes[vertex].value+increment;
                            require(std::isfinite(increment) && std::isfinite(value), "nonfinite trial nodal value");
                            changed = changed || value != patch.nodes[vertex].value;
                            patch.nodes[vertex].value = value;
                        }
                        for (std::size_t index=0; index<options.freeGradients.size(); ++index)
                        {
                            const auto& dof = options.freeGradients[index];
                            auto& node = patch.nodes[dof.vertexIndex];
                            auto& derivative = dof.component == 1 ? node.dx : node.dy;
                            const double increment = damping*linear.deltaGradient[index];
                            const double next = derivative+increment;
                            require(std::isfinite(increment) && std::isfinite(next), "nonfinite trial first derivative");
                            changed = changed || next != derivative;
                            derivative = next;
                        }
                    }
                    require(changed, "trial correction has no representable nodal change");
                    Trial trial = evaluate(patch);
                    require(trial.pressureSigns == pressureSigns, "trial changed a fixed pressure-port phase");
                    // An internal iterate is not an accepted physical step.
                    // New cut cells must enter the next linearization before
                    // their flux/volume mismatch can be corrected. evaluate()
                    // still requires actual raw volumes in [0,capacity]; the
                    // return gate also requires every flux-derived target to
                    // be in bounds. Never clip either quantity or drop a row.
                    // Honor the complete stopping criterion first. Before the
                    // balance gates pass, decrease the freshly evaluated L1.
                    // Afterwards preserve those SAME gates and decrease the
                    // total target-range violation instead. Do not combine the
                    // two merit rules with OR: that could trade improvements
                    // back and forth without making the final state feasible.
                    const bool trialBalance = balanceConverged(trial.residual,options);
                    const bool currentBalance = balanceConverged(current.residual,options);
                    const bool finalAcceptable = trial.targetsInBounds && trialBalance;
                    const bool progress = currentBalance
                        ? trialBalance && trial.targetViolation < current.targetViolation
                        : trial.residual.sumAbsolute < current.residual.sumAbsolute;
                    if (!finalAcceptable && !progress)
                    {
                        std::ostringstream reason;
                        reason << std::setprecision(17) << "fresh correction trial made no admissible progress"
                            << "; criterion=" << (currentBalance ? "bounded-target feasibility" : "complete L1 balance")
                            << "; currentL1=" << current.residual.sumAbsolute
                            << ", trialL1=" << trial.residual.sumAbsolute
                            << ", currentTargetViolation=" << current.targetViolation
                            << ", trialTargetViolation=" << trial.targetViolation
                            << ", trialRelative=" << trial.residual.relativeBalance
                            << ", trialMaxCell=" << trial.residual.maxRelativeCell
                            << ", trialTargetsInBounds=" << trial.targetsInBounds
                            << "; worst={" << describeCell(trial.residual.worstCell,
                                trial.residual.worstOldVolume,trial.residual.worstVolume,
                                trial.residual.worstTransport,trial.residual.worstCapacity) << '}';
                        if (!trial.targetsInBounds) reason << "; " << trial.targetFailure;
                        throw std::runtime_error(reason.str());
                    }
                    double maxCorrection = 0.0, maxGradientCorrection = 0.0;
                    if constexpr(Quartic)
                    {
                        for(const auto& id:options.freeData)
                        {
                            const double change=quarticDatum(trial.state.patch,id)-quarticDatum(advectedPatch,id);
                            require(std::isfinite(change),"nonfinite total quartic correction");
                            auto& maximum=id.component==0?maxCorrection:maxGradientCorrection;
                            maximum=std::max(maximum,std::abs(change));
                        }
                    }
                    else
                    {
                        for (std::size_t vertex=0; vertex<active.size(); ++vertex)
                            if (active[vertex]) maxCorrection = std::max(maxCorrection,
                                std::abs(trial.state.patch.nodes[vertex].value-advectedPatch.nodes[vertex].value));
                        require(std::isfinite(maxCorrection), "nonfinite total nodal correction");
                        for (const auto& dof : options.freeGradients)
                        {
                            const auto& node = trial.state.patch.nodes[dof.vertexIndex];
                            const auto& initial = advectedPatch.nodes[dof.vertexIndex];
                            const double change = dof.component == 1 ? node.dx-initial.dx : node.dy-initial.dy;
                            require(std::isfinite(change), "nonfinite total first-derivative correction");
                            maxGradientCorrection = std::max(maxGradientCorrection,std::abs(change));
                        }
                    }
                    current = std::move(trial);
                    report.final = current.residual;
                    report.maxCorrection = maxCorrection;
                    report.maxGradientCorrection = maxGradientCorrection;
                    accepted = true;
                }
                catch (const std::exception& error) { lastFailure = error.what(); }
                // Keep observer failures outside the numerical-trial catch:
                // a logging error must not masquerade as rejected geometry.
                if (accepted)
                {
                    notify(Progress::Stage::TrialAccepted,backtrack,damping,
                        current.targetsInBounds ? std::string() : current.targetFailure);
                    break;
                }
                notify(Progress::Stage::TrialRejected,backtrack,damping,lastFailure);
            }
            require(accepted, "binary line search exhausted without admissible balance/feasibility progress; last trial: "+lastFailure);
        }
    }
    catch (const std::exception& error)
    {
        const char* prefix = completeComponents ? "[InterfaceVolumeCorrection::correctQuarticComponents] "
            : Quartic ? "[InterfaceVolumeCorrection::correctQuarticNative] " : layout == SourceLayout::NativeLeaves
            ? "[InterfaceVolumeCorrection::correctNative] " : "[InterfaceVolumeCorrection::correct] ";
        if (const auto* resolution=dynamic_cast<const Quantum::Geometry::InterfaceGeometryResolution*>(&error))
            throw Quantum::Geometry::InterfaceGeometryResolution(resolution->reason(),resolution->cellIndex(),
                prefix+std::string(error.what())+"; "+describe(report));
        throw std::runtime_error(prefix+std::string(error.what())+"; "+describe(report));
    }
}
} // namespace

Correction correct(const QuadtreeMesh& mesh, const InterfaceTransportDomain& transport,
    const InterfaceCellMeasure::Measure& oldMeasure,
    const InterfacePhaseFlux::FluxBudget& oldFlux,
    const InterfaceNormalFluxView& currentVelocity, double dt, const Patch& advectedPatch,
    const std::vector<Vector2D>& pinnedCoordinates, const std::vector<std::uint8_t>& fixedValue,
    const Options& options,const CubicTemporalInput* temporal)
{
    return correctImpl<false>(mesh,transport,oldMeasure,oldFlux,currentVelocity,dt,advectedPatch,
        pinnedCoordinates,fixedValue,options,SourceLayout::UniformLeaves,false,nullptr,nullptr,nullptr,nullptr,temporal);
}

Correction correctNative(const QuadtreeMesh& mesh, const InterfaceTransportDomain& transport,
    const InterfaceCellMeasure::Measure& oldMeasure,
    const InterfacePhaseFlux::FluxBudget& oldFlux,
    const InterfaceNormalFluxView& currentVelocity, double dt, const Patch& advectedPatch,
    const std::vector<Vector2D>& pinnedCoordinates, const std::vector<std::uint8_t>& fixedValue,
    const Options& options)
{
    return correctImpl<false>(mesh,transport,oldMeasure,oldFlux,currentVelocity,dt,advectedPatch,
        pinnedCoordinates,fixedValue,options,SourceLayout::NativeLeaves);
}
QuarticCorrection correctQuarticNative(const QuadtreeMesh& mesh,const InterfaceTransportDomain& transport,
    const InterfaceCellMeasure::QuarticMeasure& oldMeasure,
    const InterfacePhaseFlux::FluxBudget& oldFlux,
    const InterfaceNormalFluxView& currentVelocity,double dt,
    const Quantum::Geometry::ActiveQuarticHermitePatch& advectedPatch,
    const std::vector<Vector2D>& pinnedCoordinates,const std::vector<std::uint8_t>& fixedValue,
    const QuarticOptions& options,const QuarticTemporalInput* temporal)
{
    return correctImpl<true>(mesh,transport,oldMeasure,oldFlux,currentVelocity,dt,advectedPatch,
        pinnedCoordinates,fixedValue,options,SourceLayout::NativeLeaves,false,nullptr,nullptr,nullptr,nullptr,temporal);
}

QuarticCorrection correctQuarticComponents(const QuadtreeMesh& mesh,const InterfaceTransportDomain& transport,
    const InterfaceState& previous,
    const InterfaceNormalFluxView& previousVelocity,const InterfaceNormalFluxView& currentVelocity,double dt,
    const Quantum::Geometry::ActiveQuarticHermitePatch& advectedPatch,
    const std::vector<Vector2D>& boundaryContacts,const std::vector<std::uint8_t>& fixedValue,
    const QuarticOptions& options,const InterfaceTemporalPhaseFlux::Options& timeFluxOptions,
    const InterfaceTemporalPhaseFlux::Trajectory* trajectory)
{
    return correctImpl<true>(mesh,transport,std::get<InterfaceCellMeasure::QuarticMeasure>(previous.phase().volumeMeasure),
        previous.phase().phaseBudget,currentVelocity,dt,advectedPatch,
        boundaryContacts,fixedValue,options,SourceLayout::NativeLeaves,true,&previous,&previousVelocity,&timeFluxOptions,trajectory);
}

}
