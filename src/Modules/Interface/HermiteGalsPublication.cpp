#include "Modules/Interface/HermiteGalsPublication.h"

#include "Core/Geometry/ActiveQuarticHermitePatch.h"
#include "Core/Geometry/CellChartPoint.h"
#include "Core/Geometry/LevelSetMappedSegment.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace Quantum::Interface
{
namespace
{
void require(bool ok, const char* message) { if (!ok) throw std::invalid_argument(message); }
bool tolerance(double x) { return std::isfinite(x) && x > 0 && x < 1; }
using Chain = LevelSetInterfacePanels::NativeGraphChain;
bool axisEndpoint(const Chain& chain, std::size_t graph, double fraction)
{
    using Kind = LevelSetInterfacePanels::EndpointKind;
    return (graph == 0 && fraction == 0 && chain.beginEndpointKind() == Kind::Axis)
        || (graph + 1 == chain.graphs().size() && fraction == 1 && chain.endEndpointKind() == Kind::Axis);
}
HermiteGalsSurfaceSample checkedSample(const HermiteGalsState& state,
    LevelSetHermiteGeometry::Frame frame, bool atAxis, const HermiteGalsPublicationConfig& config)
{
    const auto& scalar = state.geometry().snapshot();
    require(frame.nativeCellIndex >= 0 && frame.hasCellChartPoint, "sample lacks native chart provenance");
    const auto bounds = scalar.nativeCellBounds(frame.nativeCellIndex);
    const auto widths = bounds[1] - bounds[0];
    const double h = std::min(widths.x, widths.y);
    require(std::isfinite(frame.levelSetResidual) && frame.levelSetResidual >= 0
        && frame.levelSetResidual <= config.relativeInterfaceResidualTolerance * h, "interface residual exceeds tolerance");
    const double norm = std::hypot(frame.normal.x, frame.normal.y);
    require(std::isfinite(norm) && std::abs(norm - 1) <= config.unitNormalTolerance, "nonunit/nonfinite interface normal");
    require(std::isfinite(frame.curvature), "nonfinite meridional curvature");
    HermiteGalsSurfaceSample result{state.identity(),state.timeLayer(),std::move(frame),0,0};
    if (state.primary().measure == HermiteGalsMeasure::AxisymmetricVolume) {
        const auto coordinate = result.frame.cellChartPoint.x;
        const bool chartAxis = bounds[0].x == 0
            && Quantum::Geometry::compareCellChartCoordinate(coordinate, {0,false}) == 0;
        if (chartAxis) {
            require(atAxis && result.frame.point.x == 0 && result.frame.normal.x == 0,
                "RZ pole requires an exact retained axis endpoint and even normal");
            result.azimuthalCurvature = result.frame.curvature;
        } else {
            require(!atAxis, "axis endpoint/chart mismatch");
            // Keep the original endpoint-relative chart for the metric too;
            // do not recover a small radius by inverting Frame.point.
            const long double width = static_cast<long double>(bounds[1].x) - bounds[0].x;
            const long double radius = coordinate.fromUpper
                ? bounds[1].x - static_cast<long double>(coordinate.offset) * width
                : bounds[0].x + static_cast<long double>(coordinate.offset) * width;
            require(std::isfinite(radius) && radius > 0, "unrepresentable positive radial coordinate");
            result.azimuthalCurvature = static_cast<double>(result.frame.normal.x / radius);
        }
    }
    result.curvature = result.frame.curvature + result.azimuthalCurvature;
    require(std::isfinite(result.azimuthalCurvature) && std::isfinite(result.curvature), "nonfinite physical curvature");
    return result;
}
}
struct HermiteGalsPublication::Data
{
    Data(const QuadtreeMesh& mesh, HermiteGalsState input, NativeLeafDomain leaves,
        LevelSetInterfacePanels::NativeComponents curves, std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> scalar,
        HermiteGalsPublicationConfig controls)
        : state(std::move(input)),domain(std::move(leaves)),topology(std::move(curves)),
          partitions(state.primary().transport==HermiteGalsTransport::FirstJet
              ? InterfacePhaseGeometry::MappedFacePartitions::prepareNativeTrace(mesh,domain,scalar)
              : InterfacePhaseGeometry::MappedFacePartitions(mesh,domain,scalar)),
          volumes(state.primary().transport==HermiteGalsTransport::FirstJet
              ? integrateHermiteGalsNativeTracePhaseVolumes(mesh,state,domain,partitions,topology,controls.volume)
              : integrateHermiteGalsPhaseVolumes(mesh,state,controls.volume)),config(std::move(controls)) {}
    HermiteGalsState state;
    NativeLeafDomain domain;
    LevelSetInterfacePanels::NativeComponents topology;
    InterfacePhaseGeometry::MappedFacePartitions partitions;
    HermiteGalsPhaseVolumeResult volumes;
    HermiteGalsPublicationConfig config;
    std::vector<std::vector<std::vector<HermiteGalsSurfaceQuadrature>>> quadrature;
    HermiteGalsPublicationReport report;
};
HermiteGalsPublication::HermiteGalsPublication(std::shared_ptr<const Data> data) : data_(std::move(data)) {}
const HermiteGalsPublication::Data& HermiteGalsPublication::checked() const
{
    require(static_cast<bool>(data_), "[Interface/HermiteGalsPublication] moved-from result");
    return *data_;
}
const HermiteGalsState& HermiteGalsPublication::state() const { return checked().state; }
const HermiteGalsPhaseVolumeResult& HermiteGalsPublication::volumes() const { return checked().volumes; }
const NativeLeafDomain& HermiteGalsPublication::domain() const { return checked().domain; }
const LevelSetInterfacePanels::NativeComponents& HermiteGalsPublication::topology() const { return checked().topology; }
const InterfacePhaseGeometry::MappedFacePartitions& HermiteGalsPublication::facePartitions() const { return checked().partitions; }
const std::vector<HermiteGalsSurfaceQuadrature>& HermiteGalsPublication::quadrature(std::size_t component, std::size_t graph) const
{ return checked().quadrature.at(component).at(graph); }
HermiteGalsSurfaceSample HermiteGalsPublication::sample(std::size_t component, std::size_t graph, double fraction) const
{
    const auto& d = checked();
    try {
        const auto& chain = d.topology.curves.at(component);
        require(graph < chain.graphs().size(), "graph index out of range");
        return checkedSample(d.state, chain.frameAt(static_cast<int>(graph), fraction), axisEndpoint(chain, graph, fraction), d.config);
    } catch (const std::exception& e) {
        throw std::invalid_argument("[Interface/HermiteGalsPublication geometry=" + std::to_string(d.state.identity().value())
            + " sample component=" + std::to_string(component) + " graph=" + std::to_string(graph) + "] " + e.what());
    }
}
const HermiteGalsPublicationReport& HermiteGalsPublication::report() const { return checked().report; }
void HermiteGalsPublication::validateMesh(const QuadtreeMesh& mesh) const
{
    const auto& d = checked(); d.state.validateMesh(mesh); d.partitions.validate(mesh, d.domain);
}

HermiteGalsPublication publishHermiteGalsGeometry(const QuadtreeMesh& mesh,
    HermiteGalsCandidate&& candidate, const HermiteGalsBubbleResult& bubble, const HermiteGalsPublicationConfig& config)
{ return HermiteGalsPublication::publish(mesh,std::move(candidate),&bubble,nullptr,nullptr,config); }

HermiteGalsPublication publishHermiteGalsFirstJetGeometry(const QuadtreeMesh& mesh,
    HermiteGalsCandidate&& candidate, const HermiteGalsAmrClosureResult& mixed,
    const HermiteGalsCenterPrediction& centers, const HermiteGalsPublicationConfig& config)
{ return HermiteGalsPublication::publish(mesh,std::move(candidate),nullptr,&mixed,&centers,config); }

HermiteGalsPublication HermiteGalsPublication::publish(const QuadtreeMesh& mesh,
    HermiteGalsCandidate&& candidate, const HermiteGalsBubbleResult* bubble,
    const HermiteGalsAmrClosureResult* mixed, const HermiteGalsCenterPrediction* centers,
    const HermiteGalsPublicationConfig& config)
{
    HermiteGalsIdentity id;
    const char* stage = "input";
    bool delegated = false;
    try {
        id = candidate.identity();
        const bool conservative = bubble != nullptr;
        HermiteGalsClosure direct;
        if (conservative) {
            bubble->validateCandidate(mesh, candidate);
            require(candidate.primary().transport == HermiteGalsTransport::ConservativeVolume,
                "conservative publication requires independent inventory");
        } else {
            candidate.validateMesh(mesh);
            require(candidate.primary().transport == HermiteGalsTransport::FirstJet && candidate.primary().phaseVolumes.empty(),
                "FirstJet publication requires state without inventory");
            require(!candidate.hasGeometry(), "FirstJet publication requires unpublished candidate");
            require(mixed && centers && mixed->primaryIdentity() == id && centers->primaryIdentity == id,
                "FirstJet primary revision mismatch");
            mixed->skeleton().validateMesh(mesh);
            const auto time = candidate.timeLayer(), mt = mixed->timeLayer();
            require(time.step == mt.step && time.begin == mt.begin && time.end == mt.end && time.position == mt.position,
                "FirstJet mixed closure time mismatch");
            require(mixed->vertexMixedDerivatives().size() == candidate.primary().vertices.size()
                && centers->values.size() == candidate.primary().geometryCells.size(), "FirstJet closure layout mismatch");
            for (double value : centers->values) require(std::isfinite(value), "FirstJet nonfinite center prediction");
            direct = {mixed->vertexMixedDerivatives(),centers->values};
        }
        require((!conservative || tolerance(config.relativeInventoryTolerance)) && tolerance(config.relativeInterfaceResidualTolerance)
            && tolerance(config.unitNormalTolerance), "invalid publication tolerance");
        const auto& primary = candidate.primary(); const auto& closure = conservative ? bubble->closure() : direct;
        Quantum::Geometry::ActiveQuarticHermitePatch patch;
        const double nan = std::numeric_limits<double>::quiet_NaN();
        patch.activeCells.assign(mesh.cellIndexNum(), 0);
        patch.nodes.assign(mesh.vertIndexNum(), {nan,nan,nan,nan});
        patch.cellCenters.assign(mesh.cellIndexNum(), nan);
        for (std::size_t k = 0; k < primary.vertices.size(); ++k) {
            const auto& v = primary.vertices[k];
            patch.nodes[v.vertex] = {v.phi,v.dx,v.dy,closure.vertexMixedDerivatives.at(k)};
        }
        for (const auto& e : primary.edges) patch.edgeMidpoints.push_back({e.vertices,e.phi,e.dn});
        for (std::size_t k = 0; k < primary.geometryCells.size(); ++k) {
            const int c = primary.geometryCells[k]; patch.activeCells[c] = 1;
            patch.cellCenters[c] = closure.cellBubbleCoefficients.at(k);
        }
        const bool radial = primary.measure == HermiteGalsMeasure::AxisymmetricVolume;
        stage = "private snapshot";
        auto scalar = LevelSetHermiteGeometry::Snapshot::fromNativeLeaves(mesh, patch, radial);
        delegated = true;
        return publishScalar(mesh,std::move(candidate),closure,std::move(scalar),config);
    } catch (const std::exception& e) {
        if (delegated) throw;
        throw std::invalid_argument("[Interface/HermiteGalsPublication geometry=" + std::to_string(id.value())
            + " stage=" + stage + "] " + e.what());
    }
}

HermiteGalsPublication HermiteGalsPublication::publishScalar(const QuadtreeMesh& mesh,
    HermiteGalsCandidate&& candidate, const HermiteGalsClosure& closure,
    LevelSetHermiteGeometry::Snapshot&& snapshot, const HermiteGalsPublicationConfig& config)
{
    HermiteGalsIdentity id;
    const char* stage = "private snapshot";
    try {
        id = candidate.identity(); candidate.validateMesh(mesh);
        const auto& primary = candidate.primary();
        const bool conservative = primary.transport == HermiteGalsTransport::ConservativeVolume;
        const bool radial = primary.measure == HermiteGalsMeasure::AxisymmetricVolume;
        require((!conservative || tolerance(config.relativeInventoryTolerance)) && tolerance(config.relativeInterfaceResidualTolerance)
            && tolerance(config.unitNormalTolerance), "invalid publication tolerance");
        auto prepared = candidate.publicationCopy();
        prepared.bindGeometry(mesh, id, closure, std::move(snapshot));
        auto state = std::move(prepared).finish();
        // This const allocation is the final returned scalar, never a trial
        // clone followed by a different 'equivalent' publication snapshot.
        auto owner = std::make_shared<const HermiteGalsGeometryView>(state.geometry());
        const std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> scalar(owner, &owner->snapshot());
        stage = "geometric inventory";
        stage = "same-source topology and cuts";
        FaceSegmentGraph graph;
        graph.build(mesh, [&](const QuadtreeMesh&, const Cell& c) {
            return std::binary_search(primary.geometryCells.begin(),primary.geometryCells.end(),c.index()); });
        auto domain = NativeLeafDomain::build(mesh, graph, radial);
        auto topology = conservative ? LevelSetInterfacePanels::discoverNativeComponents(mesh, scalar, radial,
            config.volume.boundaryContacts, config.volume.symmetryPlanes)
            : LevelSetInterfacePanels::discoverNativeTraceComponents(mesh, scalar, radial,config.volume.boundaryContacts,config.volume.symmetryPlanes);
        // The cut cache retains the domain's address: bind only after domain
        // reaches its final immutable allocation, never to a movable local.
        auto data = std::make_shared<HermiteGalsPublication::Data>(mesh, state, std::move(domain),
            std::move(topology), scalar, config);
        require(&data->partitions.source() == scalar.get() && &data->volumes.geometry().snapshot() == scalar.get(), "derived geometry changed scalar owner");
        const auto& measured=data->volumes;
        require(!conservative || measured.cells().size() == primary.phaseVolumes.size(), "inventory support mismatch");
        long double maxInventoryResidual = conservative ? 0 : std::numeric_limits<long double>::quiet_NaN();
        for (std::size_t k = 0; conservative && k < measured.cells().size(); ++k) {
            const auto& cell = measured.cells()[k]; const auto target = primary.phaseVolumes[k];
            require(target.cell == cell.cell && target.liquid >= 0 && target.liquid <= cell.capacity, "inventory bounds/owner mismatch");
            const auto residual = std::abs(cell.rawLiquid - target.liquid) / cell.capacity;
            require(std::isfinite(residual) && residual <= config.relativeInventoryTolerance, "geometric inventory residual exceeds tolerance");
            maxInventoryResidual = std::max(maxInventoryResidual, residual);
        }
        stage = "normal/curvature and quadrature";
        std::size_t graphCount = 0, nodeCount = 0, axisCount = 0, halfFaces = 0;
        double maxResidual = 0, maxNormalDefect = 0, maxCurvature = 0;
        const auto record = [&](const HermiteGalsSurfaceSample& sample) {
            const auto& f = sample.frame; const auto b = scalar->nativeCellBounds(f.nativeCellIndex); const auto h = b[1] - b[0];
            maxResidual = std::max(maxResidual, f.levelSetResidual / std::min(h.x,h.y));
            maxNormalDefect = std::max(maxNormalDefect, std::abs(std::hypot(f.normal.x,f.normal.y) - 1));
            maxCurvature = std::max(maxCurvature, std::abs(sample.curvature));
        };
        for (const auto& chain : data->topology.curves) {
            require(chain.usesSource(*scalar) && chain.axisymmetric() == radial, "curve changed source/measure");
            std::vector<std::vector<HermiteGalsSurfaceQuadrature>> component;
            for (std::size_t g = 0; g < chain.graphs().size(); ++g) {
                ++graphCount;
                for (double fraction : {0.,1.}) {
                    const bool axis = axisEndpoint(chain,g,fraction); axisCount += axis;
                    record(checkedSample(state, chain.frameAt(static_cast<int>(g),fraction), axis, config));
                }
                std::vector<HermiteGalsSurfaceQuadrature> nodes;
                for (const auto& q : chain.quadrature(static_cast<int>(g))) {
                    require(std::isfinite(q.weight) && q.weight > 0 && std::isfinite(q.physicalWeight)
                        && q.physicalWeight > 0, "invalid interface quadrature weight");
                    auto sample = checkedSample(state, q.frame, false, config); record(sample);
                    nodes.push_back({std::move(sample),q.weight,q.physicalWeight}); ++nodeCount;
                }
                require(!nodes.empty(), "empty native graph quadrature");
                component.push_back(std::move(nodes));
            }
            data->quadrature.push_back(std::move(component));
        }
        for (const auto f : data->domain.faces(mesh)) if (f.negativeCell >= 0 || f.positiveCell >= 0)
            for (int half = 0; half < 2; ++half) {
                require(conservative ? data->partitions.halfFace(f.nativeRow,half).usesSource(*scalar)
                    : data->partitions.traceHalfFace(f.nativeRow,half).usesSource(*scalar), "half-face cut changed scalar owner");
                ++halfFaces;
            }
        data->report = {data->topology.curves.size(),graphCount,nodeCount,halfFaces,axisCount,
            maxInventoryResidual,maxResidual,maxNormalDefect,maxCurvature};
        HermiteGalsPublication result(std::move(data));
        // Everything that can fail is complete. No caller-visible candidate
        // mutation or acceptance occurred before this no-throw consumption.
        candidate.consumePublication();
        return result;
    } catch (const std::exception& e) {
        throw std::invalid_argument("[Interface/HermiteGalsPublication geometry=" + std::to_string(id.value())
            + " stage=" + stage + "] " + e.what());
    }
}
}
