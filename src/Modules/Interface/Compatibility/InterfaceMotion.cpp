#include "Modules/Interface/Compatibility/Detail/InterfaceStateInternal.h"
#include "Modules/Interface/Compatibility/InterfaceMotion.h"
#include "Core/Geometry/ActiveLevelSetTransport.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/UniformLeafGrid.h"
#include "Modules/Interface/Compatibility/Detail/InterfaceMotionSupport.h"
#include "Core/Mesh/NativeLeafDomain.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <vector>

struct InterfaceMotionAccess
{
    static void bindMaterialMap(InterfaceMotionResult& result,
        std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source,
        std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> target, const InterfaceStepTime& time,
        std::function<Vector2D(Vector2D)> forward)
    {
        result.materialSource_ = std::move(source);
        result.materialTarget_ = std::move(target);
        result.materialTime_ = time;
        result.forwardMaterialPoint_ = std::move(forward);
    }
};

std::function<Vector2D(Vector2D)> InterfaceMotionResult::materialMapFor(
    const std::shared_ptr<const LevelSetHermiteGeometry::Snapshot>& source,
    const std::shared_ptr<const LevelSetHermiteGeometry::Snapshot>& target, const InterfaceStepTime& time) const
{
    if (!forwardMaterialPoint_ || !source || !target || source != materialSource_
        || target != materialTarget_ || time.index != materialTime_.index || time.begin != materialTime_.begin
        || time.end != materialTime_.end || time.dt != materialTime_.dt)
        throw std::invalid_argument("Interface material trajectory requires its retained source, final target and step times");
    return forwardMaterialPoint_;
}

namespace
{
void require(bool valid, const char *reason)
{
    if (!valid)
        throw std::runtime_error(std::string("[InterfaceMotion] ") + reason);
}
bool same(Vector2D a, Vector2D b)
{
    return a.x == b.x && a.y == b.y;
}
double area(const LevelSetInterfacePanels &panels)
{
    double result = 0.0;
    for (const auto &panel : panels.panels())
    {
        require(std::isfinite(panel.physicalMeasure) && panel.physicalMeasure > 0.0,
                "nonpositive or nonfinite free-panel physical measure");
        result += panel.physicalMeasure;
    }
    require(std::isfinite(result) && result > 0.0, "nonfinite total free-surface area");
    return result;
}
} // namespace

void InterfaceMotion::validatePinnedToAxisChain(const LevelSetInterfacePanels &panels, Vector2D pin)
{
    using Kind = LevelSetInterfacePanels::EndpointKind;
    require(panels.axisymmetric() && !panels.closed() && !panels.panels().empty() &&
                panels.beginEndpointKind() == Kind::DeclaredBoundary &&
                panels.endEndpointKind() == Kind::Axis && panels.boundaryContacts().size() == 1 &&
                same(panels.boundaryContacts().front(), pin) && same(panels.panels().front().begin, pin) &&
                panels.panels().back().end.x == 0.0,
            "free interface must retain one directed pinned-to-axis chain");
}

namespace
{
template<class Patch, class Measure>
InterfaceMotionResult prepareMotion(
    const QuadtreeMesh &mesh, const InterfaceTransportDomain &transport, const InterfaceStateAccess::Data &previous,
    const InterfaceStepTime& time, const LevelSetBoundaryTrace::Velocity &newVelocity,
    const LevelSetBoundaryTrace::Velocity &oldVelocity, double requiredHalfWidthCells,
    const Patch &previousPatch, const Measure &oldMeasure,
    const InterfacePhaseFlux::FluxBudget &oldFlux, const InterfaceNormalFluxView &currentVelocity,
    const InterfaceNormalFluxView* previousFlux,
    double linearTolerance, const InterfaceVolumeProgress &progress,
    const InterfaceGeometryObserver &observeGeometry, bool independentQueries, InterfaceStateAccess::Data& destination)
{
    const double dt = time.dt;
    const auto& domain=transport.geometry();
    const auto motionStarted=std::chrono::steady_clock::now();
    domain.validateMeshIdentity(mesh);
    constexpr bool quartic = std::is_same_v<Patch, Quantum::Geometry::ActiveQuarticHermitePatch>;
    const bool native = !domain.uniformLeaves();
    require(native == quartic, "moving scalar representation differs from its native/uniform binding");
    std::optional<UniformLeafGrid> grid;
    if (!native)
        grid.emplace(mesh, CartesianGridScope::UniformLeaves);
    require(domain.axisymmetric() && mesh.minCoord.x == 0.0 && previous.geometry() && previous.panels(),
            "a retained active RZ scalar and free chain are required");
    require(previous.pinnedCoordinates.size() == 1 && previous.pinnedVertexIndices.size() == 1,
            "this moving binding requires one retained pinned contact");
    require(previousPatch.activeCells.size() == static_cast<std::size_t>(mesh.cellIndexNum()) &&
                previousPatch.nodes.size() == static_cast<std::size_t>(mesh.vertIndexNum()) &&
                previous.activeVertex.size() == static_cast<std::size_t>(mesh.vertIndexNum()) &&
                previous.freeLevelSet.size() == mesh.vertIndexNum() && previous.freeLevelSet.value,
            "retained active scalar storage mismatch");
    const auto &source = *previous.geometry();
    require(source.restrictedDomain() && source.radialAxisSymmetry() && source.usesNativeLeaves() == native &&
                same(source.minCoord(), mesh.minCoord) && same(source.maxCoord(), mesh.maxCoord) &&
                (native ? source.nativeCellCount() == static_cast<std::size_t>(mesh.cellIndexNum())
                        : source.nx() == grid->nx() && source.ny() == grid->ny() &&
                              same(source.cellSize(), grid->spacing())),
            "retained scalar layout differs from the fitted domain");
    const Vector2D pin = previous.pinnedCoordinates.front();
    InterfaceMotion::validatePinnedToAxisChain(*previous.panels(), pin);
    if (observeGeometry)
        observeGeometry(InterfaceMotionStage::Previous, *previous.panels());
    const int pinned = previous.pinnedVertexIndices.front();
    require(pinned >= 0 && pinned < mesh.vertIndexNum() && previous.activeVertex[pinned] &&
                previousPatch.nodes[pinned].value == 0.0 && previous.freeLevelSet.value[pinned] == 0.0,
            "retained pinned scalar was changed or is not active");
    require(std::isfinite(requiredHalfWidthCells) && requiredHalfWidthCells > 0.0,
            "positive finite capillary half width is required");
    const auto minimumSpacing =
        native ? previous.panels()->nativeChain().minimumActiveSpacing() : grid->spacing();
    const double maxh = std::max(minimumSpacing.x, minimumSpacing.y);
    ActiveLevelSetTransport::Options options;
    const InterfaceMotionSupport::Band initialBand(requiredHalfWidthCells, quartic, maxh);
    options.radialAxisSymmetry = true;
    options.fixedZeroVertices = previous.pinnedVertexIndices;

    std::vector<std::array<int, 2>> logical(mesh.cellIndexNum());
    std::vector<std::uint8_t> activeVertices(mesh.vertIndexNum(), 0);
    for (int position = 0; position < mesh.cellIndexNum(); ++position)
    {
        const int k = native ? position : mesh.cells(grid->cells()[position]).index();
        if (!native)
            logical[k] = {{position % grid->nx(), position / grid->nx()}};
        const bool active = domain.cellRowByIndex()[k] >= 0;
        const bool scalarActive =
            native ? source.nativeCellActive(k) : source.cellActive(logical[k][0], logical[k][1]);
        if (native)
        {
            const auto &cell = mesh.cells(mesh.cellIndexToLabel[k]);
            const auto bounds = source.nativeCellBounds(k);
            require(same(bounds[0], mesh.verts(cell.cellToVert[0]).coord()) &&
                        same(bounds[1], mesh.verts(cell.cellToVert[2]).coord()),
                    "retained native leaf bounds changed");
        }
        require(previousPatch.activeCells[k] <= 1 &&
                    static_cast<bool>(previousPatch.activeCells[k]) == active &&
                    scalarActive == active,
                "retained active patch changed the fixed solid mask");
        if (active)
            for (int label : mesh.cells(mesh.cellIndexToLabel[k]).cellToVert)
                activeVertices[mesh.verts(label).index()] = 1;
    }
    require(activeVertices == previous.activeVertex,
            "retained vertex mask differs from actual active-cell corners");
    for (int k = 0; k < mesh.vertIndexNum(); ++k)
        if (previous.activeVertex[k])
            require(previous.activeVertex[k] == 1 &&
                        (native && mesh.verts(mesh.vertIndexToLabel[k]).supCell >= 0
                             ? source.value(mesh.verts(mesh.vertIndexToLabel[k]).coord())
                             : previousPatch.nodes[k].value) == previous.freeLevelSet.value[k],
                    "display scalar differs from authoritative active nodal value");

    InterfaceMotionSupport support(mesh, transport, source, grid ? &*grid : nullptr,
        logical, maxh, initialBand);
    support.prepare(*previous.panels());
    options.fullWidth = support.band().fullWidth;
    options.zeroWidth = support.band().zeroWidth;
    const auto& pressureCells = support.pressureCells();

    // Identify the wet wall by its retained negative trace, not by nozzle
    // coordinates or a fitted contact angle. Its incident cell supplies only
    // the single-phase extension; the other side retains the evolved cap jet.
    const auto normal = previous.panels()->panels().front().beginFrame.normal;
    std::array<int, 2> pinLogical{};
    if (!native)
    {
        const auto vertexPosition =
            std::find(grid->vertices().begin(), grid->vertices().end(), mesh.vertIndexToLabel[pinned]);
        require(vertexPosition != grid->vertices().end(), "pin is absent from the actual leaf layout");
        const auto offset = static_cast<int>(vertexPosition - grid->vertices().begin());
        pinLogical = {{offset % (grid->nx() + 1), offset / (grid->nx() + 1)}};
    }
    std::optional<ActiveLevelSetTransport::FixedWallCornerExtension> wetExtension;
    for (const auto &face : domain.faces(mesh))
    {
        if (transport.boundary(face.nativeRow) != InterfaceTransportDomain::Boundary::StationaryWall)
            continue;
        const auto &edge = face.horizontal ? mesh.hEdges(face.edgeLabel) : mesh.vEdges(face.edgeLabel);
        const int a = mesh.verts(edge.edgeToVert[0]).index(), b = mesh.verts(edge.edgeToVert[1]).index();
        const int other = a == pinned ? b : b == pinned ? a : -1;
        if (other < 0 || previous.freeLevelSet.value[other] >= 0.0)
            continue;
        require(!wetExtension, "this contact binding requires one unique wet wall at the pin");
        const int owner = face.negativeCell >= 0 ? face.negativeCell : face.positiveCell;
        require(owner >= 0 && domain.cellRowByIndex()[owner] >= 0, "wet wall has no active owner");
        unsigned corner = 4;
        const auto &ownerCell = mesh.cells(mesh.cellIndexToLabel[owner]);
        for (unsigned c = 0; c < 4; ++c)
            if (mesh.verts(ownerCell.cellToVert[c]).index() == pinned)
                corner = c;
        require(corner < 4, "wet wall owner is not incident to the pin");
        require(source.closedCellSignExceptVertices(mesh.cells(mesh.cellIndexToLabel[owner]).centerCoord(),
                                                    1u << corner) == -1,
                "wet corner is not certified negative except at the pin; a new contact representation is "
                "required");
        wetExtension = ActiveLevelSetTransport::FixedWallCornerExtension{owner, pinned, -1, face.horizontal};
    }
    require(wetExtension.has_value(), "no declared wet no-slip corner extension was found");
    options.fixedWallCornerExtensions.push_back(*wetExtension);
    const int seamAxis = wetExtension->wallHorizontal ? 0 : 1;
    const auto wetCenter = mesh.cells(mesh.cellIndexToLabel[wetExtension->cellIndex]).centerCoord();
    const int wetSide = native ? ((seamAxis == 0 ? wetCenter.x < pin.x : wetCenter.y < pin.y) ? -1 : 1)
                               : (logical[wetExtension->cellIndex][seamAxis] < pinLogical[seamAxis] ? -1 : 1);
    // Reconstruct increments on each smooth side of the explicit seam.
    // Its direction comes from the wall/cell topology, not the initial or
    // current free-surface angle. On the seam the BASE jet is the cap jet.
    options.oneSidedSamples.assign(mesh.vertIndexNum(), 0);
    options.extensionSideSamples.assign(mesh.vertIndexNum(), 0);
    for (int position = 0; position < mesh.vertIndexNum(); ++position)
    {
        const int k = native ? position : mesh.verts(grid->vertices()[position]).index();
        if (!previous.activeVertex[k])
            continue;
        const auto point = mesh.verts(mesh.vertIndexToLabel[k]).coord();
        const double side =
            native ? (seamAxis == 0 ? point.x - pin.x : point.y - pin.y) * wetSide
                   : ((seamAxis == 0 ? position % (grid->nx() + 1) : position / (grid->nx() + 1)) -
                      pinLogical[seamAxis]) *
                         wetSide;
        options.oneSidedSamples[k] = side <= 0;
        options.extensionSideSamples[k] = side >= 0;
    }
    if constexpr (quartic)
    {
        options.quarticStencils = previous.transportStencils ? previous.transportStencils
            : std::make_shared<ActiveQuarticHermiteReconstruction::StencilCache>();
        // The provider explicitly declares read-only concurrent queries.
        // Prescribed callbacks remain sequential unless they make that promise.
        options.independentVelocityQueries = independentQueries;
    }
    const auto transportStarted = std::chrono::steady_clock::now();
    auto transported = [&]() {
        if constexpr (std::is_same_v<Patch, Quantum::Geometry::ActiveQuarticHermitePatch>)
            return ActiveLevelSetTransport::advanceQuarticNative(mesh, previousPatch, dt, newVelocity, oldVelocity, options);
        else
            return ActiveLevelSetTransport::advance(mesh, previousPatch, dt, newVelocity, oldVelocity, options);
    }();
    std::optional<InterfaceTemporalPhaseFlux::PhaseTrajectory<Patch>> trajectory;
    std::optional<InterfaceVolumeCorrection::TemporalInput<Patch>> temporal;
    if (previousFlux) {
        trajectory.emplace();
        trajectory->samples[0]=previousPatch;
        trajectory->samples[4]=transported.patch;
        for (int k=1;k<4;++k) {
            const double theta=.25*k;
            const LevelSetBoundaryTrace::Velocity velocity=[oldVelocity,newVelocity,theta](const Vector2D& point) {
                return oldVelocity(point)*(1-theta)+newVelocity(point)*theta;
            };
            if constexpr (quartic)
                trajectory->samples[k]=ActiveLevelSetTransport::advanceQuarticNative(
                    mesh,previousPatch,theta*dt,velocity,oldVelocity,options).patch;
            else trajectory->samples[k]=ActiveLevelSetTransport::advance(
                mesh,previousPatch,theta*dt,velocity,oldVelocity,options).patch;
        }
        temporal.emplace(InterfaceVolumeCorrection::TemporalInput<Patch>{
            previousPatch,previous.geometry(),*previousFlux,&*trajectory,{}});
    }
    const auto transportFinished=std::chrono::steady_clock::now();
    const double transportSeconds =
        std::chrono::duration<double>(transportFinished - transportStarted).count();
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> geometry(std::move(transported.geometry));
    auto panels =
        native ? std::make_unique<LevelSetInterfacePanels>(LevelSetInterfacePanels::discoverNativeGraphs(
                     mesh, geometry, true, previous.pinnedCoordinates))
               : std::make_unique<LevelSetInterfacePanels>(mesh, geometry, true, previous.pinnedCoordinates,
                                                           CartesianGridScope::UniformLeaves);
    InterfaceMotion::validatePinnedToAxisChain(*panels, pin);
    if (observeGeometry)
        observeGeometry(InterfaceMotionStage::Transported, *panels);
    support.checkCandidate(*panels, *geometry, InterfaceMotionSupport::Stage::Transported);
    for (const auto &face : domain.faces(mesh))
        if (transport.boundary(face.nativeRow) == InterfaceTransportDomain::Boundary::StationaryWall)
        {
            const auto &edge = face.horizontal ? mesh.hEdges(face.edgeLabel) : mesh.vEdges(face.edgeLabel);
            for (int label : edge.edgeToVert)
            {
                const int k = mesh.verts(label).index();
                require(transported.patch.nodes[k].value == previousPatch.nodes[k].value,
                        "a fixed no-slip wall scalar value moved");
            }
        }

    // Translate declared interface boundaries; the correction kernel receives only
    // native fixed-value identities. Axis values remain free so the pole moves;
    // radial symmetry is enforced by the retained zero dx/dxy, not a pinned phi.
    std::vector<std::uint8_t> fixedValue(mesh.vertIndexNum(), 0);
    for (int k = 0; k < mesh.vertIndexNum(); ++k)
        fixedValue[k] = !previous.activeVertex[k];
    for (const auto &face : domain.faces(mesh))
        if (transport.boundary(face.nativeRow) == InterfaceTransportDomain::Boundary::StationaryWall)
        {
            const auto &edge = face.horizontal ? mesh.hEdges(face.edgeLabel) : mesh.vEdges(face.edgeLabel);
            for (int label : edge.edgeToVert)
                fixedValue[mesh.verts(label).index()] = 1;
        }
    for (int k : pressureCells)
        for (int label : mesh.cells(mesh.cellIndexToLabel[k]).cellToVert)
            fixedValue[mesh.verts(label).index()] = 1;
    for (int k : previous.pinnedVertexIndices)
        fixedValue[k] = 1;

    std::conditional_t<quartic, InterfaceVolumeCorrection::QuarticOptions,
        InterfaceVolumeCorrection::Options> correctionOptions;
    correctionOptions.linearTolerance = linearTolerance;
    // New iteration criteria, not a replacement for the historical global
    // volume CSV or a declaration that a coupled run is physically validated.
    correctionOptions.relativeBalanceTolerance = 1e-3;
    correctionOptions.maxCellRelativeTolerance = 1e-10;
    correctionOptions.progress = progress;
    // Position pinning does not prescribe a contact angle. The wet-wall
    // sector retains its zero tangential derivative; the free cap's base
    // tangential derivative is an independent physical contact-slope DOF.
    // Scale its unknown to the same phi units as the nodal-value unknowns.
    // No other derivative is freed or refitted by this contact capability.
    const unsigned contactComponent = wetExtension->wallHorizontal ? 1u : 2u;
    const auto contactBounds =
        native ? source.nativeCellBounds(wetExtension->cellIndex) : std::array<Vector2D, 2>{};
    const double contactSpacing =
        native ? (wetExtension->wallHorizontal ? contactBounds[1].x - contactBounds[0].x
                                               : contactBounds[1].y - contactBounds[0].y)
               : (wetExtension->wallHorizontal ? grid->spacing().x : grid->spacing().y);
    if constexpr (!quartic)
        correctionOptions.freeGradients.push_back({pinned, contactComponent, 1.0 / contactSpacing});
    std::set<std::array<int, 2>> fixedCoarseSides;
    InterfaceMotionResult result;
    if (native)
    {
        // A fixed dependent boundary value constrains its real coarse-edge
        // endpoint values, not the unused hanging payload. Derivatives remain
        // fixed except the independent contact slope explicitly freed above.
        std::function<void(int)> fixSupport = [&](int k) {
            const auto &vertex = mesh.verts(mesh.vertIndexToLabel[k]);
            if (vertex.supCell < 0)
            {
                fixedValue[k] = 1;
                return;
            }
            const auto &coarse = mesh.cells(vertex.supCell);
            const auto point = vertex.coord();
            int endpoints = 0;
            std::array<int, 2> coarseSide{{-1, -1}};
            const auto lo = mesh.verts(coarse.cellToVert[0]).coord(),
                       hi = mesh.verts(coarse.cellToVert[2]).coord();
            const bool vertical = point.x == lo.x || point.x == hi.x;
            require(vertical || point.y == lo.y || point.y == hi.y,
                    "fixed hanging value is not on a coarse edge");
            for (int label : coarse.cellToVert)
            {
                const auto p = mesh.verts(label).coord();
                if (vertical ? p.x == point.x : p.y == point.y)
                {
                    const int support = mesh.verts(label).index();
                    require(support != pinned,
                            "a dependent fixed boundary value touches the free contact slope");
                    fixSupport(support);
                    require(endpoints < 2, "fixed coarse side has extra endpoints");
                    coarseSide[endpoints++] = support;
                }
            }
            require(endpoints == 2, "fixed hanging value lacks two actual coarse-edge endpoints");
            const auto a = mesh.verts(mesh.vertIndexToLabel[coarseSide[0]]).coord();
            const auto b = mesh.verts(mesh.vertIndexToLabel[coarseSide[1]]).coord();
            if (vertical ? a.y > b.y : a.x > b.x) std::swap(coarseSide[0], coarseSide[1]);
            fixedCoarseSides.insert(coarseSide);
            fixedValue[k] = 0;
        };
        for (int k = 0; k < mesh.vertIndexNum(); ++k)
            if (previous.activeVertex[k] && fixedValue[k] &&
                mesh.verts(mesh.vertIndexToLabel[k]).supCell >= 0)
                fixSupport(k);
    }
    const auto correctionStarted=std::chrono::steady_clock::now();
    auto corrected = [&]() {
        if constexpr (std::is_same_v<Patch, Quantum::Geometry::ActiveQuarticHermitePatch>)
        {
            using Kind = Quantum::Geometry::ActiveQuarticHermiteSensitivity::Kind;
            std::map<std::pair<bool,double>,std::vector<std::pair<double,double>>> fixedWallSpans;
            for (const auto& face : domain.faces(mesh))
                if (transport.boundary(face.nativeRow) == InterfaceTransportDomain::Boundary::StationaryWall)
                {
                    const auto& wall = face.horizontal ? mesh.hEdges(face.edgeLabel) : mesh.vEdges(face.edgeLabel);
                    const auto a = mesh.verts(wall.edgeToVert[0]).coord(), b = mesh.verts(wall.edgeToVert[1]).coord();
                    fixedWallSpans[{face.horizontal,face.horizontal ? a.y : a.x}].push_back(
                        face.horizontal ? std::make_pair(a.x,b.x) : std::make_pair(a.y,b.y));
                }
            // A fixed hanging value also fixes its true coarse-side midpoint.
            // Pressure cells retain ALL their scalar values, not only corners.
            constexpr int pairs[4][2] = {{0,1},{3,2},{0,3},{1,2}};
            for (int c : pressureCells)
            {
                const auto& cell = mesh.cells(mesh.cellIndexToLabel[c]);
                for (const auto& side : pairs)
                    fixedCoarseSides.insert({{mesh.verts(cell.cellToVert[side[0]]).index(),
                                              mesh.verts(cell.cellToVert[side[1]]).index()}});
            }
            for (int k = 0; k < mesh.vertIndexNum(); ++k)
                if (previous.activeVertex[k] && !fixedValue[k] && mesh.verts(mesh.vertIndexToLabel[k]).supCell < 0)
                    correctionOptions.freeData.push_back({Kind::Vertex,k,0,1});
            for (std::size_t i = 0; i < transported.patch.edgeMidpoints.size(); ++i)
            {
                const auto& edge = transported.patch.edgeMidpoints[i];
                if (!std::isfinite(edge.value)) continue; // Declared dependent payload from Core.
                bool fixed = fixedCoarseSides.count(edge.vertices) != 0;
                const auto a = mesh.verts(mesh.vertIndexToLabel[edge.vertices[0]]).coord();
                const auto b = mesh.verts(mesh.vertIndexToLabel[edge.vertices[1]]).coord();
                const bool horizontal = a.y == b.y;
                const auto spans = fixedWallSpans.find({horizontal,horizontal ? a.y : a.x});
                if (spans != fixedWallSpans.end())
                    for (const auto& span : spans->second)
                        fixed |= std::max(horizontal ? a.x : a.y,span.first)
                            < std::min(horizontal ? b.x : b.y,span.second);
                if (!fixed) correctionOptions.freeData.push_back({Kind::EdgeMidpoint,static_cast<int>(i),0,1});
                else if (edge.value != previousPatch.edgeMidpoints[i].value)
                {
                    const auto point = (a + b) * 0.5;
                    std::ostringstream message;
                    message << std::setprecision(17)
                        << "[InterfaceMotion] a fixed no-slip/pressure/coarse-side midpoint moved"
                        << "; edge=" << i << "; vertices=" << edge.vertices[0] << ',' << edge.vertices[1]
                        << "; point=(" << point.x << ',' << point.y << "); old="
                        << previousPatch.edgeMidpoints[i].value << "; transported=" << edge.value;
                    throw std::runtime_error(message.str());
                }
            }
            for (int c : domain.cellIndexByRow())
                if (!pressureCells.count(c)) correctionOptions.freeData.push_back({Kind::CellCenter,c,0,1});
            correctionOptions.freeData.push_back({Kind::Vertex,pinned,static_cast<int>(contactComponent),1.0/contactSpacing});
            return InterfaceVolumeCorrection::correctQuarticNative(mesh, transport, oldMeasure, oldFlux,
                currentVelocity, dt, transported.patch, previous.pinnedCoordinates, fixedValue, correctionOptions,
                temporal?&*temporal:nullptr);
        }
        else
            return InterfaceVolumeCorrection::correct(mesh, transport, oldMeasure, oldFlux, currentVelocity, dt,
                transported.patch, previous.pinnedCoordinates, fixedValue, correctionOptions,temporal?&*temporal:nullptr);
    }();
    const auto correctionFinished=std::chrono::steady_clock::now();
    if (native)
    {
        require(corrected.nativeChain && !corrected.panels,
                "native correction returned mixed geometry ownership");
        corrected.panels = std::make_unique<LevelSetInterfacePanels>(std::move(*corrected.nativeChain));
        corrected.nativeChain.reset();
    }
    InterfaceMotion::validatePinnedToAxisChain(*corrected.panels, pin);
    if (observeGeometry)
        observeGeometry(InterfaceMotionStage::VolumeCorrected, *corrected.panels);
    support.checkCandidate(*corrected.panels, *corrected.geometry, InterfaceMotionSupport::Stage::Corrected);
    for (int k = 0; k < mesh.vertIndexNum(); ++k)
        if (previous.activeVertex[k] && fixedValue[k])
            require(corrected.patch.nodes[k].value == transported.patch.nodes[k].value,
                    "volume correction changed a fixed wall, contact or pressure-port value");
    auto &report = result.report;
    report.transportSeconds = transportSeconds;
    if constexpr (quartic)
    {
        report.transportStencilRequests = transported.reconstruction.stencilRequests;
        report.quarticTimings = transported.timings;
        report.transportStencilBuilds = transported.reconstruction.stencilBuilds;
        report.transportScalarRequests = transported.scalarIncrementRequests;
        report.transportScalarEvaluations = transported.scalarIncrementEvaluations;
    }
    report.activeNodes = transported.activeNodes;
    report.taperedNodes = transported.taperedNodes;
    report.maxTraceRatio = transported.maxTraceRatio;
    for (int k = 0; k < mesh.vertIndexNum(); ++k)
        if (previous.activeVertex[k] && (!native || mesh.verts(mesh.vertIndexToLabel[k]).supCell < 0))
        {
            const double increment = corrected.patch.nodes[k].value - previousPatch.nodes[k].value;
            report.changedValues += increment != 0.0;
            report.stationaryNodes += increment == 0.0;
            report.maxValueIncrement = std::max(report.maxValueIncrement, std::abs(increment));
        }
    report.supportCells = static_cast<int>(support.supportCells().size());
    report.pressureCells = static_cast<int>(pressureCells.size());
    report.dt = dt;
    report.fullWidth = options.fullWidth;
    report.zeroWidth = options.zeroWidth;
    report.areaBefore = area(*previous.panels());
    report.areaAfter = area(*corrected.panels);
    report.pinNormalBefore = normal;
    report.pinNormalAfter = corrected.panels->panels().front().beginFrame.normal;
    report.poleBefore = previous.panels()->panels().back().end;
    report.poleAfter = corrected.panels->panels().back().end;
    auto &state = destination;
    if constexpr (quartic) state.transportStencils = options.quarticStencils;
    state.freeLevelSet = previous.freeLevelSet;
    state.activeVertex = previous.activeVertex;
    state.pinnedVertexIndices = previous.pinnedVertexIndices;
    state.pinnedCoordinates = previous.pinnedCoordinates;
    state.initialContactJets = previous.initialContactJets;
    state.meshIdentity = previous.meshIdentity;
    for (int k = 0; k < mesh.vertIndexNum(); ++k)
        if (state.activeVertex[k])
            state.freeLevelSet.value[k] =
                native && mesh.verts(mesh.vertIndexToLabel[k]).supCell >= 0
                    ? corrected.geometry->value(mesh.verts(mesh.vertIndexToLabel[k]).coord())
                    : corrected.patch.nodes[k].value;
    state.scalarPatch = std::move(corrected.patch);
    state.retainGeometry(std::move(corrected.geometry), std::move(corrected.panels));
    state.phase.emplace(InterfacePhaseState{InterfaceVolumeMeasure(std::move(corrected.measure)), std::move(corrected.flux)});
    if (previousFlux) {
        require(static_cast<bool>(corrected.stepFlux),"pinned correction lost its integrated transport ledger");
        state.phase->transport.emplace(InterfacePhaseState::Transport{
            time,previous.geometry(),state.geometry(),std::move(*corrected.stepFlux)});
    }
    result.conservation = corrected.report;
    require(transported.forwardMaterialPoint && same(transported.forwardMaterialPoint(pin), pin),
            "the common material trajectory did not preserve the fixed no-slip contact");
    InterfaceMotionAccess::bindMaterialMap(result, previous.geometry(), state.geometry(), time,
        std::move(transported.forwardMaterialPoint));
    report.preparationSeconds=std::chrono::duration<double>(transportStarted-motionStarted).count();
    report.correctionPreparationSeconds=std::chrono::duration<double>(correctionStarted-transportFinished).count();
    report.correctionSeconds=std::chrono::duration<double>(correctionFinished-correctionStarted).count();
    report.finalizationSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-correctionFinished).count();
    return result;
}
} // namespace


InterfaceMotionResult InterfaceMotion::preparePinnedToAxis(
    const QuadtreeMesh& mesh, const InterfaceTransportDomain& transport, const InterfaceState& previousState,
    const InterfaceStepTime& time, const InterfaceMotionVelocityLayers& velocities, const InterfaceNormalFluxView& currentFlux,
    double requiredHalfWidthCells, double linearTolerance,
    InterfaceState& destinationState, const InterfaceVolumeProgress& progress,
    const InterfaceGeometryObserver& observeGeometry,const InterfaceNormalFluxView* previousFlux)
{
    require(previousState.stage() == InterfaceState::Stage::Ready
        && destinationState.stage() == InterfaceState::Stage::Empty, "interface motion requires ready source and fresh candidate");
    require(!previousState.hasConservativePhase(),"legacy pinned correction cannot advance conservative inventory");
    previousState.validateMesh(mesh);
    const auto completed = previousState.progress();
    const double timeScale = std::max({std::abs(time.begin), std::abs(time.end), time.dt});
    require(time.index > 0 && completed.completedSteps < std::numeric_limits<int>::max()
        && time.index == completed.completedSteps + 1 && time.begin == completed.endTime
        && std::isfinite(time.begin) && time.begin >= 0.0 && std::isfinite(time.end) && time.end > time.begin
        && std::isfinite(time.dt) && time.dt > 0.0
        && std::abs((time.end-time.begin)-time.dt) <= 16.0*std::numeric_limits<double>::epsilon()*timeScale,
        "interface step index/times do not continue its retained history and numerical dt");
    const auto& oldPhase = previousState.phase();
    const auto& oldMeasure = oldPhase.volumeMeasure;
    const auto& oldFlux = oldPhase.phaseBudget;
    const auto& previous = InterfaceStateAccess::get(previousState);
    auto& destination = InterfaceStateAccess::get(destinationState);
    InterfacePreparation preparation(destination);
    auto result = [&]() {
    require(velocities.current && velocities.previous, "interface motion requires both velocity time levels");
    require(&destination != &previous && !destination.geometry() && !destination.panels()
        && destination.freeLevelSet.size() == mesh.vertIndexNum(),
        "interface motion requires a fresh distinct destination on the current mesh");
    if (transport.geometry().uniformLeaves())
    {
        const auto* patch = std::get_if<Quantum::Geometry::ActiveHermitePatch>(&previous.scalarPatch);
        const auto* measure = std::get_if<InterfaceCellMeasure::Measure>(&oldMeasure);
        require(patch && measure, "uniform motion requires its retained cubic scalar and measure");
        return prepareMotion(mesh,transport,previous,time,velocities.current,velocities.previous,requiredHalfWidthCells,
            *patch,*measure,oldFlux,currentFlux,previousFlux,linearTolerance,progress,observeGeometry,
            velocities.independentQueries,destination);
    }
    const auto* patch = std::get_if<Quantum::Geometry::ActiveQuarticHermitePatch>(&previous.scalarPatch);
    const auto* measure = std::get_if<InterfaceCellMeasure::QuarticMeasure>(&oldMeasure);
    require(patch && measure, "native motion requires its retained quartic scalar and complete measure");
    return prepareMotion(mesh,transport,previous,time,velocities.current,velocities.previous,requiredHalfWidthCells,
        *patch,*measure,oldFlux,currentFlux,previousFlux,linearTolerance,progress,observeGeometry,
        velocities.independentQueries,destination);
    }();
    destination.progress = {time.index, time.end};
    preparation.complete(InterfaceState::Stage::Ready);
    return result;
}
