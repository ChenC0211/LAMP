#include "Modules/Interface/Compatibility/InterfaceInitialization.h"
#include "Modules/Interface/Compatibility/Detail/InterfaceStateInternal.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/CartesianLeafLocator.h"
#include "Core/Mesh/UniformLeafGrid.h"
#include "Core/Geometry/LevelSetCurvature.h"
#include "Core/Geometry/ActiveQuarticHermiteReconstruction.h"
#include "Core/Geometry/LevelSetQueryReconstruction.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <type_traits>

namespace {
void require(bool valid, const char* message) {
    if (!valid) throw std::invalid_argument(std::string("[InterfaceInitialization] ") + message);
}
bool same(Vector2D a, Vector2D b) { return a.x == b.x && a.y == b.y; }
std::vector<std::uint8_t> initialContactSamples(const QuadtreeMesh &mesh, const UniformLeafGrid &grid,
                                                const InterfaceStateAccess::Data &state, Vector2D point,
                                                Vector2D inwardTangent)
{
    // Select the cap-side half-plane from its declared INITIAL primitive.
    // This is a local derivative support, not a clipping/contact-angle law.
    // In particular the endpoint's distance cone and solid payload are not
    // sampled. No circle radius/analytic curvature is inserted into the jet.
    std::vector<std::uint8_t> samples(state.activeVertex.size(), 0);
    for (int j = 0; j <= grid.ny(); ++j)
        for (int i = 0; i <= grid.nx(); ++i)
        {
            const int k =
                mesh.verts(grid.vertices()[static_cast<std::size_t>(j) * (grid.nx() + 1) + i]).index();
            if (!state.activeVertex[k])
                continue;
            const double x = i == grid.nx() ? grid.maxCoord().x : grid.minCoord().x + i * grid.spacing().x;
            const double y = j == grid.ny() ? grid.maxCoord().y : grid.minCoord().y + j * grid.spacing().y;
            const long double side = (static_cast<long double>(x) - point.x) * inwardTangent.x +
                                     (static_cast<long double>(y) - point.y) * inwardTangent.y;
            samples[k] = side >= 0.0L;
        }
    return samples;
}

InterfaceInitialContactJet initialContactJet(const QuadtreeMesh &mesh, const UniformLeafGrid &grid,
                                               const InterfaceStateAccess::Data &state, Vector2D point,
                                               Vector2D inwardTangent)
{
    const auto samples = initialContactSamples(mesh, grid, state, point, inwardTangent);
    const LevelSetCurvature::Reconstruction<QuadtreeMesh> reconstruction(
        mesh, state.freeLevelSet, 4, LevelSetCurvature::Basis::TensorProduct, true,
        CartesianGridScope::UniformLeaves, samples);
    const auto derivatives = reconstruction.evaluate(point);
    require(point.x > 0.0, "a wall contact jet at the axis requires a separate regular axis limit");
    InterfaceInitialContactJet jet;
    jet.inwardTangent = inwardTangent;
    jet.normal = Vector2D(derivatives.nx, derivatives.ny);
    jet.dx = derivatives.dx;
    jet.dy = derivatives.dy;
    jet.dxx = derivatives.dxx;
    jet.dxy = derivatives.dxy;
    jet.dyy = derivatives.dyy;
    jet.meridionalCurvature = -derivatives.curvature;
    jet.axisymmetricCurvature = jet.meridionalCurvature + derivatives.nx / point.x;
    require(std::isfinite(jet.axisymmetricCurvature), "non-finite initial contact RZ curvature");
    jet.admittedSamples = static_cast<int>(std::count(samples.begin(), samples.end(), 1));
    return jet;
}

// Setup binding only: native geometry chooses support; the shared Core Q4
// kernel reconstructs the same INITIAL primitive distance. No evolving field
// is sampled here, and no initial sampler survives into the time loop.
void initializeNativeJets(
    const QuadtreeMesh &mesh, const InterfaceTransportDomain& transport, const InterfaceInitialCurveInput& input,
    InterfaceStateAccess::Data& state,
    const std::vector<Vector2D> &tangents,
    const std::function<InterfaceInitialCurveSample(Vector2D)> &primitive)
{
    const auto& fluid = transport.geometry();
    const auto faces = fluid.faces(mesh);
    require(state.pinnedCoordinates.size() == 1 && tangents.size() == 1,
            "native initial cap requires one declared pin and its initial tangent");
    const CartesianLeafLocator locator(mesh);
    const auto pin = state.pinnedCoordinates.front(), tangent = tangents.front();
    std::vector<std::uint8_t> capSamples(mesh.vertIndexNum(), 0);
    for (int k = 0; k < mesh.vertIndexNum(); ++k)
    {
        const auto point = mesh.verts(mesh.vertIndexToLabel[k]).coord();
        const long double side = (static_cast<long double>(point.x) - pin.x) * tangent.x +
                                 (static_cast<long double>(point.y) - pin.y) * tangent.y;
        capSamples[k] = state.activeVertex[k] && side >= 0;
    }
    const double infinity = std::numeric_limits<double>::infinity();
    std::vector<Vector2D> spacing(mesh.vertIndexNum(), {infinity, infinity});
    std::vector<std::array<Vector2D, 2>> bounds(mesh.cellIndexNum());
    std::vector<std::vector<int>> wallsByCell(mesh.cellIndexNum());
    Quantum::Geometry::ActiveQuarticHermitePatch patch;
    patch.activeCells.assign(mesh.cellIndexNum(), 0);
    for (int c = 0; c < mesh.cellIndexNum(); ++c)
    {
        const auto &cell = mesh.cells(mesh.cellIndexToLabel[c]);
        bounds[c] = {mesh.verts(cell.cellToVert[0]).coord(), mesh.verts(cell.cellToVert[2]).coord()};
        if (fluid.cellRowByIndex()[c] < 0)
            continue;
        patch.activeCells[c] = 1;
        const auto h = bounds[c][1] - bounds[c][0];
        for (int label : cell.cellToVert)
        {
            auto &local = spacing[mesh.verts(label).index()];
            local.x = std::min(local.x, h.x);
            local.y = std::min(local.y, h.y);
        }
    }
    for (int f = 0; f < static_cast<int>(faces.size()); ++f)
        if (transport.boundary(f) == InterfaceTransportDomain::Boundary::StationaryWall)
        {
            const auto &face = faces[f];
            wallsByCell[face.negativeCell >= 0 ? face.negativeCell : face.positiveCell].push_back(f);
        }
    const auto sample = [&](Vector2D point) {
        const auto raw = primitive(point);
        if (same(point, pin))
        {
            require(raw.unsignedDistance == 0, "initial pin sample is not an exact primitive zero");
            return 0.0;
        }
        bool negative = raw.negative;
        int wallPhase = -1;
        const auto incident = locator.incident(point);
        bool active = false;
        for (int n = 0; n < incident.count; ++n)
        {
            const int c = incident.indices[n];
            if (!patch.activeCells[c])
                continue;
            active = true;
            for (int f : wallsByCell[c])
            {
                const auto &face = faces[f];
                const auto &edge =
                    face.horizontal ? mesh.hEdges(face.edgeLabel) : mesh.vEdges(face.edgeLabel);
                const auto a = mesh.verts(edge.edgeToVert[0]).coord(),
                           b = mesh.verts(edge.edgeToVert[1]).coord();
                if (point.x < std::min(a.x, b.x) || point.x > std::max(a.x, b.x) ||
                    point.y < std::min(a.y, b.y) || point.y > std::max(a.y, b.y))
                    continue;
                const int phase = input.phaseByCell[c];
                require(wallPhase < 0 || wallPhase == phase,
                        "initial query meets an undeclared wall contact");
                wallPhase = phase;
            }
        }
        require(active, "initial derivative query has no actual active-fluid cell");
        if (wallPhase >= 0)
        {
            require((wallPhase == input.negativePhase || wallPhase == input.positivePhase) &&
                        raw.unsignedDistance > 0,
                    "invalid initial fluid-side wall trace");
            negative = wallPhase == input.negativePhase;
        }
        return negative ? -raw.unsignedDistance : raw.unsignedDistance;
    };
    ActiveQuarticHermiteReconstruction::Options dataOptions;
    dataOptions.radialAxisSymmetry = true;
    dataOptions.oneSidedSamples = capSamples;
    patch = ActiveQuarticHermiteReconstruction::reconstruct(mesh, patch.activeCells, sample, dataOptions);
    // Keep the historical one-sided 5x5 initial-contact diagnostic and its
    // resolution gate. It is not the new scalar's authority or a dynamic
    // curvature/angle condition; the full quartic data above use Core's
    // existing recursive-extension generator, including edge/center data.
    const LevelSetPolynomial::QueryReconstruction reconstruction(true);
    for (int k = 0; k < mesh.vertIndexNum(); ++k)
    {
        const auto &vertex = mesh.verts(mesh.vertIndexToLabel[k]);
        if (!state.activeVertex[k] || vertex.supCell >= 0)
            continue; // dependent data are not DOFs
        const auto point = vertex.coord();
        require(patch.nodes[k].value == state.freeLevelSet.value[k],
                "native quartic initial sampler differs from the existing nodal distance/sign");
        if (k != state.pinnedVertexIndices.front())
            continue;
        const auto &mask = capSamples[k] ? capSamples : state.activeVertex;
        const auto admitted = [&](Vector2D lo, Vector2D hi) {
            if (!std::isfinite(lo.x) || !std::isfinite(lo.y) || !std::isfinite(hi.x) ||
                !std::isfinite(hi.y) || lo.x < mesh.minCoord.x || lo.y < mesh.minCoord.y ||
                hi.x > mesh.maxCoord.x || hi.y > mesh.maxCoord.y || !(lo.x < hi.x && lo.y < hi.y))
                return false;
            bool covered = false;
            for (int c : locator.intersectingBox(lo, hi))
            {
                if (std::max(lo.x, bounds[c][0].x) >= std::min(hi.x, bounds[c][1].x) ||
                    std::max(lo.y, bounds[c][0].y) >= std::min(hi.y, bounds[c][1].y))
                    continue;
                if (!patch.activeCells[c])
                    return false;
                for (int label : mesh.cells(mesh.cellIndexToLabel[c]).cellToVert)
                    if (!mask[mesh.verts(label).index()])
                        return false;
                covered = true;
            }
            return covered;
        };
        try
        {
            long double variation = 0;
            const auto derivatives =
                reconstruction.evaluateDerivatives(point, spacing[k], admitted, sample, variation);
            const Quantum::Geometry::ActiveHermitePatch::NodalDof node{
                state.freeLevelSet.value[k], derivatives.dx, derivatives.dy, derivatives.dxy};
            if (k == state.pinnedVertexIndices.front())
            {
                const double magnitude = std::hypot(node.dx, node.dy);
                const long double resolutionFloor = 64.0L * std::numeric_limits<double>::epsilon() *
                                                    reconstruction.conditionEstimate() * variation /
                                                    std::min(spacing[k].x, spacing[k].y);
                require(point.x > 0 && std::isfinite(magnitude) && std::isfinite(resolutionFloor) &&
                            magnitude > resolutionFloor && variation > 0,
                        "degenerate native initial contact jet");
                InterfaceInitialContactJet jet;
                jet.inwardTangent = tangent;
                jet.normal = {node.dx / magnitude, node.dy / magnitude};
                jet.dx = node.dx;
                jet.dy = node.dy;
                jet.dxx = derivatives.dxx;
                jet.dxy = derivatives.dxy;
                jet.dyy = derivatives.dyy;
                jet.meridionalCurvature =
                    (jet.normal.y * jet.normal.y * jet.dxx - 2 * jet.normal.x * jet.normal.y * jet.dxy +
                     jet.normal.x * jet.normal.x * jet.dyy) /
                    magnitude;
                jet.axisymmetricCurvature = jet.meridionalCurvature + jet.normal.x / point.x;
                require(std::isfinite(jet.axisymmetricCurvature),
                        "nonfinite native initial contact curvature");
                jet.admittedSamples = static_cast<int>(std::count(capSamples.begin(), capSamples.end(), 1));
                state.initialContactJets.push_back(jet);
            }
        }
        catch (const std::exception &error)
        {
            std::ostringstream message;
            message.precision(17);
            message << "[InterfaceInitialization] native initial vertex=" << k << ", point=(" << point.x << ','
                    << point.y << "): " << error.what();
            throw std::runtime_error(message.str());
        }
    }
    state.scalarPatch = std::move(patch);
}
} // namespace

void InterfaceInitialization::prepareCurve(const QuadtreeMesh& mesh, const InterfaceTransportDomain& transport,
    const InterfaceInitialCurveInput& input, InterfaceState& destination)
{
    const auto& fluid = transport.geometry();
    fluid.validateMeshIdentity(mesh);
    require(destination.stage() == InterfaceState::Stage::Empty && destination.vertexCount() == mesh.vertIndexNum(),
        "initial curve requires a fresh destination on its mesh");
    require(fluid.axisymmetric() && mesh.minCoord.x == 0.0, "initial fitted curve requires the actual RZ axis");
    require(input.phaseByCell.size() == static_cast<std::size_t>(mesh.cellIndexNum())
        && input.negativePhase != input.positivePhase && !input.pinnedCoordinates.empty(), "invalid initial material/contact input");
    for (const auto* field : {&input.unsignedDistance, &input.compositeSign}) {
        require(field->value && field->size() == mesh.vertIndexNum(), "initial curve input storage mismatch");
        for (int i=0; i<field->size(); ++i) require(field->label(i) == mesh.vertIndexToLabel[i], "stale initial curve input labels");
    }
    require(fluid.uniformLeaves() || static_cast<bool>(input.primitive), "native initial curve requires its primitive sampler");
    auto& state = InterfaceStateAccess::get(destination);
    InterfacePreparation preparation(state);
    state.meshIdentity = fluid.meshIdentity();
    std::optional<UniformLeafGrid> grid;
    if (fluid.uniformLeaves()) grid.emplace(mesh, CartesianGridScope::UniformLeaves);
    const auto begin = input.begin, end = input.end;
    for (int cell : fluid.cellIndexByRow())
        for (int label : mesh.cells(mesh.cellIndexToLabel[cell]).cellToVert)
            state.activeVertex[mesh.verts(label).index()] = 1;
    std::vector<std::uint8_t> wallVertex(mesh.vertIndexNum(), 0), pinned(mesh.vertIndexNum(), 0);
    for (const auto &face : fluid.faces(mesh))
        if (transport.boundary(face.nativeRow) == InterfaceTransportDomain::Boundary::StationaryWall)
        {
            const auto &edge = face.horizontal ? mesh.hEdges(face.edgeLabel) : mesh.vEdges(face.edgeLabel);
            for (int label : edge.edgeToVert)
                wallVertex[mesh.verts(label).index()] = 1;
        }
    for (const Vector2D point : input.pinnedCoordinates)
    {
        require(same(point, begin) || same(point, end),
                "pin must be an actual endpoint of the selected free chain");
        const auto address = [](double coordinate, double lower, double upper, double h, int count) {
            const double rounded = std::round((coordinate - lower) / h);
            require(std::isfinite(rounded) && rounded >= 0.0 && rounded <= count,
                    "pinned point is outside the canonical leaf grid");
            const int index = static_cast<int>(rounded);
            require(coordinate == (index == count ? upper : lower + index * h),
                    "pin must be an EXACT canonical grid point; nearest-node snapping is unsupported");
            return index;
        };
        int label = -1;
        if (grid)
        {
            const int x =
                address(point.x, grid->minCoord().x, grid->maxCoord().x, grid->spacing().x, grid->nx());
            const int y =
                address(point.y, grid->minCoord().y, grid->maxCoord().y, grid->spacing().y, grid->ny());
            label = grid->vertices()[static_cast<std::size_t>(y) * (grid->nx() + 1) + x];
        }
        else
        {
            for (int k = 0; k < mesh.vertIndexNum(); ++k)
                if (same(mesh.verts(mesh.vertIndexToLabel[k]).coord(), point))
                {
                    require(label < 0, "native pin has duplicate actual mesh vertices");
                    label = mesh.vertIndexToLabel[k];
                }
            require(label >= 0 && mesh.verts(label).supCell < 0,
                    "pin must be an EXACT independent native vertex; snapping is unsupported");
        }
        const int vertexIndex = mesh.verts(label).index();
        require(state.activeVertex[vertexIndex] && wallVertex[vertexIndex] && !pinned[vertexIndex],
                "pin must identify one unique actual fluid/solid boundary vertex");
        // Logical identity alone is insufficient at a solid contact: a
        // one-ULP difference can put an evolving cap below the actual wall.
        // Fresh meshes use the shared dyadic coordinate generator; reject
        // inconsistent legacy geometry instead of moving or snapping a pin.
        const auto actual = mesh.verts(label).coord();
        require(same(actual, point), "pinned Point differs from the actual mesh vertex; rebuild with "
                                     "consistent Cartesian coordinates");
        pinned[vertexIndex] = 1;
        state.pinnedVertexIndices.push_back(vertexIndex);
        state.pinnedCoordinates.push_back(point);
    }
    for (const Vector2D point : {begin, end})
    {
        bool declared = false;
        for (const auto contact : state.pinnedCoordinates)
            declared |= same(point, contact);
        require(declared || point.x == mesh.minCoord.x,
                "free-chain endpoint is neither pinned wall contact nor axis");
    }
    // Reuse the existing Core geometric primitives, selecting only free lines.
    // Sign on a fixed wall comes from its actual fluid-side material trace, NOT
    // parity on the sign-only solid closure or the old wall's displaced zero.
    const auto& unsignedDistance = input.unsignedDistance;
    for (int i = 0; i < mesh.vertIndexNum(); ++i)
    {
        if (!state.activeVertex[i])
            continue;
        const double distance = unsignedDistance.value[i];
        require(std::isfinite(distance) && distance >= 0.0 && std::isfinite(input.compositeSign.value[i]),
                "nonfinite or negative initial distance/sign input");
        bool negative = input.compositeSign.value[i] < 0.0;
        if (pinned[i])
        {
            require(distance == 0.0,
                    "canonical pin representative has nonzero computed initial curve distance");
            state.freeLevelSet.value[i] = 0.0;
            continue;
        }
        if (wallVertex[i])
        {
            int phase = -1;
            for (int label : mesh.verts(mesh.vertIndexToLabel[i]).vertToCell)
            {
                if (label < 0)
                    continue;
                const int cell = mesh.cells(label).index();
                if (fluid.cellRowByIndex()[cell] < 0)
                    continue;
                const int material = input.phaseByCell[cell];
                require(material == input.negativePhase || material == input.positivePhase,
                        "nonfluid active material");
                require(phase < 0 || phase == material, "undeclared contact on the fixed solid boundary");
                phase = material;
            }
            require(phase >= 0 && distance > 0.0,
                    "fixed wall trace lacks a fluid phase or touches an undeclared free interface");
            negative = phase == input.negativePhase;
        }
        state.freeLevelSet.value[i] = negative ? -distance : distance;
    }
    std::vector<Vector2D> contactTangents;
    for (const auto point : state.pinnedCoordinates)
    {
        const Vector2D tangent = same(point, begin) ? input.inwardEndpointTangents[0] : input.inwardEndpointTangents[1];
        contactTangents.push_back(tangent);
        if (grid)
            state.initialContactJets.push_back(initialContactJet(mesh, *grid, state, point, tangent));
    }
    if (!grid)
        initializeNativeJets(mesh, transport, input, state, contactTangents, input.primitive);
    preparation.complete(InterfaceState::Stage::Prepared);
}

void InterfaceInitialization::publishCurve(const QuadtreeMesh& mesh, const InterfaceTransportDomain& transport,
    InterfaceState& destination, const InterfaceInitialGeometryObserver& observer)
{
    const auto& fluid = transport.geometry();
    require(destination.stage() == InterfaceState::Stage::Prepared, "initial scalar requires prepared curve data");
    destination.validateMesh(mesh);
    auto& state = InterfaceStateAccess::get(destination);
    require(state.meshIdentity == fluid.meshIdentity(), "initial scalar domain incarnation changed");
    InterfacePreparation preparation(state);
    require(!state.geometry() && !state.panels(), "initial active geometry already exists");
    require(state.initialContactJets.size() == 1 && state.pinnedCoordinates.size() == 1,
            "current fitted initial cap binding requires one pinned contact and one axis endpoint");
    fluid.validateMeshIdentity(mesh);
    InterfaceScalarPatch prepared;
    if (fluid.uniformLeaves())
    {
        auto &patch = std::get<Quantum::Geometry::ActiveHermitePatch>(prepared);
        const UniformLeafGrid grid(mesh, CartesianGridScope::UniformLeaves);
        const auto capSamples = initialContactSamples(mesh, grid, state, state.pinnedCoordinates.front(),
                                                      state.initialContactJets.front().inwardTangent);
        using Reconstruction = LevelSetCurvature::Reconstruction<QuadtreeMesh>;
        const Reconstruction fluidDerivatives(mesh, state.freeLevelSet, 4,
                                              LevelSetCurvature::Basis::TensorProduct, true,
                                              CartesianGridScope::UniformLeaves, state.activeVertex);
        const Reconstruction capDerivatives(mesh, state.freeLevelSet, 4,
                                            LevelSetCurvature::Basis::TensorProduct, true,
                                            CartesianGridScope::UniformLeaves, capSamples);
        patch.activeCells.assign(mesh.cellIndexNum(), 0);
        for (int cell : fluid.cellIndexByRow())
            patch.activeCells[cell] = 1;
        patch.nodes.resize(mesh.vertIndexNum());
        for (int j = 0; j <= grid.ny(); ++j)
            for (int i = 0; i <= grid.nx(); ++i)
            {
                const int k =
                    mesh.verts(grid.vertices()[static_cast<std::size_t>(j) * (grid.nx() + 1) + i]).index();
                if (!state.activeVertex[k])
                    continue;
                const Vector2D point(
                    i == grid.nx() ? grid.maxCoord().x : grid.minCoord().x + i * grid.spacing().x,
                    j == grid.ny() ? grid.maxCoord().y : grid.minCoord().y + j * grid.spacing().y);
                const auto jet =
                    (capSamples[k] ? capDerivatives : fluidDerivatives).evaluateDerivatives(point);
                auto &node = patch.nodes[k];
                node.value = state.freeLevelSet.value[k];
                node.dx = jet.dx;
                node.dy = jet.dy;
                node.dxy = jet.dxy;
            }
    }
    else
    {
        require(std::holds_alternative<Quantum::Geometry::ActiveQuarticHermitePatch>(state.scalarPatch),
                "native initial quartic data have not been prepared");
        prepared = state.scalarPatch;
    }
    // INITIAL data compatibility: phi is constant along the declared source
    // curve, so its endpoint tangential derivative is zero. A truncation-error
    // tangent can otherwise create a second zero on a fixed wall arbitrarily
    // close to the pin. Retain the Q4 normal component and mixed derivative;
    // do not insert curvature or prescribe this angle for subsequent motion.
    // The raw one-sided Q4 diagnostic is deliberately left unchanged.
    std::visit([&](auto &patch) {
        for (std::size_t p = 0; p < state.pinnedVertexIndices.size(); ++p)
        {
            auto &node = patch.nodes[state.pinnedVertexIndices[p]];
            const auto tangent = state.initialContactJets[p].inwardTangent;
            const double tangential = node.dx * tangent.x + node.dy * tangent.y;
            const double normal = node.dx * tangent.y - node.dy * tangent.x;
            require(node.value == 0.0 && std::isfinite(normal) && normal != 0.0,
                    "initial pinned curve lacks a zero value or a nondegenerate normal derivative");
            node.dx = normal * tangent.y;
            node.dy = -normal * tangent.x;
            if (observer.initialTangencyAdjusted)
                observer.initialTangencyAdjusted(state.pinnedVertexIndices[p], tangential, normal);
        }
    }, prepared);
    // Retain the actual DOFs even if strict zero-topology construction rejects
    // them; these are the scalar inputs, not a fabricated successful surface.
    if (observer.patchPrepared)
        observer.patchPrepared(prepared);
    auto geometry = std::visit([&](const auto &patch) {
        if constexpr (std::is_same_v<std::decay_t<decltype(patch)>, Quantum::Geometry::ActiveQuarticHermitePatch>)
            return std::make_shared<const LevelSetHermiteGeometry::Snapshot>(
                LevelSetHermiteGeometry::Snapshot::fromNativeLeaves(mesh, patch, true));
        else
            return std::make_shared<const LevelSetHermiteGeometry::Snapshot>(
                mesh, patch, true, CartesianGridScope::UniformLeaves);
    }, prepared);
    if (observer.scalarConstructed)
        observer.scalarConstructed(*geometry);
    auto panels =
        fluid.uniformLeaves()
            ? std::make_unique<LevelSetInterfacePanels>(mesh, geometry, true, state.pinnedCoordinates,
                                                        CartesianGridScope::UniformLeaves)
            : std::make_unique<LevelSetInterfacePanels>(LevelSetInterfacePanels::discoverNativeGraphs(
                  mesh, geometry, true, state.pinnedCoordinates));
    require(!panels->closed() &&
                panels->beginEndpointKind() == LevelSetInterfacePanels::EndpointKind::DeclaredBoundary &&
                panels->endEndpointKind() == LevelSetInterfacePanels::EndpointKind::Axis,
            "active free chain must be liquid-to-gas oriented from the pin to the axis");
    if (observer.panelsConstructed)
        observer.panelsConstructed(*panels);
    state.scalarPatch = std::move(prepared);
    state.retainGeometry(std::move(geometry), std::move(panels));
    if (state.geometry()->usesNativeLeaves())
        for (int k = 0; k < mesh.vertIndexNum(); ++k)
            if (state.activeVertex[k] && mesh.verts(mesh.vertIndexToLabel[k]).supCell >= 0)
                state.freeLevelSet.value[k] =
                    state.geometry()->value(mesh.verts(mesh.vertIndexToLabel[k]).coord());
    preparation.complete(InterfaceState::Stage::Ready);
}
