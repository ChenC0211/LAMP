#include "Modules/Interface/HermiteGalsState.h"

#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Mesh/NativeMeshIdentity.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Geometry/Detail/ExactExpansion.h"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>

namespace LAMP::Interface
{
namespace Detail
{
struct HermiteGalsData
{
    std::shared_ptr<const NativeMeshIdentity> mesh;
    HermiteGalsIdentity identity;
    HermiteGalsTimeLayer time;
    HermiteGalsPrimaryState primary;
    HermiteGalsClosure closure;
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> geometry;
    HermiteGalsData(const QuadtreeMesh& inputMesh, HermiteGalsTimeLayer inputTime,
        const HermiteGalsPrimaryState& inputPrimary)
        : mesh(std::make_shared<const NativeMeshIdentity>(inputMesh)), time(inputTime), primary(inputPrimary) {}
    HermiteGalsData(std::shared_ptr<const NativeMeshIdentity> inputMesh,HermiteGalsTimeLayer inputTime,
        const HermiteGalsPrimaryState& inputPrimary)
        :mesh(std::move(inputMesh)),time(inputTime),primary(inputPrimary) {}
};
}
namespace
{
void require(bool condition, const std::string& message, HermiteGalsIdentity id = {})
{
    if (!condition) throw std::invalid_argument("[Interface/HermiteGalsState id="
        + std::to_string(id.value()) + "] " + message);
}
template<class Pointer> decltype(auto) retained(const Pointer& pointer)
{
    require(static_cast<bool>(pointer), "access: moved-from or consumed handle");
    return *pointer;
}
std::uint64_t nextIdentity()
{
    static std::atomic<std::uint64_t> next{1};
    auto id = next.load(std::memory_order_relaxed);
    do {
        require(id != std::numeric_limits<std::uint64_t>::max(), "identity space exhausted");
    } while (!next.compare_exchange_weak(id, id + 1, std::memory_order_relaxed));
    return id;
}
template<class Point,class Independent>
void validatePrimaryOnTopology(int vertices,int cells,double minimumRadius,const HermiteGalsPrimaryState& p,
    HermiteGalsIdentity identity,const Point& point,const Independent& independent)
{
    const auto check = [&](bool ok, const std::string& message) {
        require(ok, "primary: " + message, identity);
    };
    check(p.measure == HermiteGalsMeasure::CartesianArea ||
        p.measure == HermiteGalsMeasure::AxisymmetricVolume, "unknown measure");
    check(p.measure != HermiteGalsMeasure::AxisymmetricVolume || minimumRadius >= 0,
        "RZ requires nonnegative radius");
    check(p.transport == HermiteGalsTransport::ConservativeVolume || p.transport == HermiteGalsTransport::FirstJet,
        "unknown transport mode");
    check(p.transport != HermiteGalsTransport::FirstJet || p.phaseVolumes.empty(),
        "FirstJet state cannot carry conserved inventory");
    int previous = -1;
    for (const auto& v : p.vertices) {
        check(v.vertex > previous && v.vertex < vertices,
            "master vertex indices must be unique, ascending and native");
        check(independent(v.vertex),
            "hanging vertex cannot be an independent first jet");
        check(std::isfinite(v.phi) && std::isfinite(v.dx) && std::isfinite(v.dy),
            "non-finite master vertex=" + std::to_string(v.vertex));
        previous = v.vertex;
    }
    std::set<std::array<int, 2>> sides;
    for (const auto& e : p.edges) {
        check(e.vertices[0] >= 0 && e.vertices[1] >= 0 &&
            e.vertices[0] < vertices && e.vertices[1] < vertices,
            "edge endpoint is not a native vertex");
        const auto a = point(e.vertices[0]);
        const auto b = point(e.vertices[1]);
        check((a.x == b.x && a.y < b.y) || (a.y == b.y && a.x < b.x),
            "edge endpoints must follow the positive Cartesian tangent");
        check(sides.insert(e.vertices).second, "duplicate master edge");
        check(std::isfinite(e.phi) && std::isfinite(e.dn), "non-finite master edge");
    }
    previous = -1;
    for (const auto& m : p.phaseVolumes) {
        check(m.cell > previous && m.cell < cells,
            "inventory cell indices must be unique, ascending and native");
        check(std::isfinite(m.liquid) && m.liquid >= 0,
            "non-finite or negative phase volume at cell=" + std::to_string(m.cell));
        previous = m.cell;
    }
    previous = -1;
    for (int cell : p.geometryCells) {
        check(cell > previous && cell < cells,
            "geometry cell indices must be unique, ascending and native");
        const auto found = std::lower_bound(p.phaseVolumes.begin(), p.phaseVolumes.end(), cell,
            [](const auto& volume, int index) { return volume.cell < index; });
        check(p.transport == HermiteGalsTransport::FirstJet || (found != p.phaseVolumes.end() && found->cell == cell),
            "geometry cell lacks conserved inventory: cell=" + std::to_string(cell));
        previous = cell;
    }
    // Complete master coverage/restriction and admissible volume bounds are
    // numerical/topological producer contracts (I06-I10), not inferred here.
}
void validatePrimary(const QuadtreeMesh& mesh,const HermiteGalsPrimaryState& p,HermiteGalsIdentity identity)
{
    validatePrimaryOnTopology(mesh.vertIndexNum(),mesh.cellIndexNum(),mesh.minCoord.x,p,identity,
        [&](int v){return mesh.verts(mesh.vertIndexToLabel[v]).coord();},
        [&](int v){return mesh.verts(mesh.vertIndexToLabel[v]).type()==TVertType::non;});
}

void validateMesh(const Detail::HermiteGalsData& data, const QuadtreeMesh& mesh)
{
    try { data.mesh->validate(mesh); }
    catch (const std::invalid_argument& e) {
        require(false, std::string("mesh binding: ") + e.what(), data.identity);
    }
}
}

HermiteGalsGeometryView::HermiteGalsGeometryView(std::shared_ptr<const Detail::HermiteGalsData> data)
    : data_(std::move(data)) {}
HermiteGalsIdentity HermiteGalsGeometryView::identity() const { return retained(data_).identity; }
HermiteGalsTimeLayer HermiteGalsGeometryView::timeLayer() const { return retained(data_).time; }
const LevelSetHermiteGeometry::Snapshot& HermiteGalsGeometryView::snapshot() const
{ return *retained(data_).geometry; }

HermiteGalsState::HermiteGalsState(std::shared_ptr<const Detail::HermiteGalsData> data)
    : data_(std::move(data)) {}
HermiteGalsIdentity HermiteGalsState::identity() const { return retained(data_).identity; }
HermiteGalsTimeLayer HermiteGalsState::timeLayer() const { return retained(data_).time; }
const HermiteGalsPrimaryState& HermiteGalsState::primary() const { return retained(data_).primary; }
const HermiteGalsClosure& HermiteGalsState::closure() const { return retained(data_).closure; }
HermiteGalsGeometryView HermiteGalsState::geometry() const
{ retained(data_); return HermiteGalsGeometryView(data_); }
void HermiteGalsState::validateMesh(const QuadtreeMesh& mesh) const
{ LAMP::Interface::validateMesh(retained(data_), mesh); }

HermiteGalsCandidate::HermiteGalsCandidate(const QuadtreeMesh& mesh, HermiteGalsTimeLayer time,
    const HermiteGalsPrimaryState& primary)
{
    time.time();
    auto data = std::make_unique<Detail::HermiteGalsData>(mesh, time, primary);
    data->identity = HermiteGalsIdentity(nextIdentity());
    validatePrimary(mesh, data->primary, data->identity);
    data_ = std::move(data);
}
HermiteGalsCandidate::~HermiteGalsCandidate() = default;
HermiteGalsCandidate::HermiteGalsCandidate(std::unique_ptr<Detail::HermiteGalsData> data) noexcept
    : data_(std::move(data)) {}
HermiteGalsCandidate HermiteGalsCandidate::publicationCopy() const
{
    const auto& data = retained(data_);
    require(!data.geometry, "publication: candidate already has geometry", data.identity);
    return HermiteGalsCandidate(std::make_unique<Detail::HermiteGalsData>(data));
}
void HermiteGalsCandidate::consumePublication() noexcept { data_.reset(); }
HermiteGalsCandidate::HermiteGalsCandidate(HermiteGalsCandidate&&) noexcept = default;
HermiteGalsCandidate& HermiteGalsCandidate::operator=(HermiteGalsCandidate&&) noexcept = default;
HermiteGalsIdentity HermiteGalsCandidate::identity() const { return retained(data_).identity; }
HermiteGalsTimeLayer HermiteGalsCandidate::timeLayer() const { return retained(data_).time; }
const HermiteGalsPrimaryState& HermiteGalsCandidate::primary() const { return retained(data_).primary; }
bool HermiteGalsCandidate::hasGeometry() const { return static_cast<bool>(retained(data_).geometry); }
void HermiteGalsCandidate::validateMesh(const QuadtreeMesh& mesh) const
{ LAMP::Interface::validateMesh(retained(data_), mesh); }
void HermiteGalsCandidate::replacePrimary(const QuadtreeMesh& mesh, const HermiteGalsPrimaryState& primary)
{
    validateMesh(mesh);
    HermiteGalsCandidate replacement(mesh, data_->time, primary);
    data_.swap(replacement.data_);
}
void HermiteGalsCandidate::bindGeometry(const QuadtreeMesh& mesh, HermiteGalsIdentity source,
    const HermiteGalsClosure& closure, LevelSetHermiteGeometry::Snapshot&& snapshot)
{
    auto& data = retained(data_);
    const auto check = [&](bool ok, const std::string& message) {
        require(ok, "geometry binding: " + message, data.identity);
    };
    check(source == data.identity, "stale/foreign primary identity=" + std::to_string(source.value()));
    check(!data.geometry, "primary revision already has its unique snapshot");
    validateMesh(mesh);
    check(closure.vertexMixedDerivatives.size() == data.primary.vertices.size(), "mixed closure size");
    check(closure.cellBubbleCoefficients.size() == data.primary.geometryCells.size(), "bubble closure size");
    for (double value : closure.vertexMixedDerivatives) check(std::isfinite(value), "non-finite mixed closure");
    for (double value : closure.cellBubbleCoefficients) check(std::isfinite(value), "non-finite bubble closure");
    check(snapshot.usesNativeLeaves() &&
        (snapshot.reconstruction() == LevelSetHermiteGeometry::Reconstruction::EnrichedQuartic
            || snapshot.usesC2Macro()),"requires native quartic geometry");
    if(snapshot.usesC2Macro())check(data.primary.transport==HermiteGalsTransport::FirstJet
        &&data.primary.measure==HermiteGalsMeasure::CartesianArea
        &&data.primary.geometryCells.size()==static_cast<std::size_t>(mesh.cellIndexNum()),
        "C2 macro requires complete XY FirstJet support");
    check(snapshot.nativeCellCount() == static_cast<std::size_t>(mesh.cellIndexNum()), "native cell count mismatch");
    for (int cell = 0; cell < mesh.cellIndexNum(); ++cell) {
        const bool active = std::binary_search(data.primary.geometryCells.begin(), data.primary.geometryCells.end(), cell);
        check(snapshot.nativeCellActive(cell) == active, "activity mismatch at cell=" + std::to_string(cell));
        if (!active) continue;
        const auto& native = mesh.cells(mesh.cellIndexToLabel[cell]);
        // Match the snapshot's native corner coordinates directly. Rebuilding
        // bounds as center +/- size/2 can differ by an ulp on valid meshes.
        const auto lower = mesh.verts(native.cellToVert[0]).coord();
        const auto upper = mesh.verts(native.cellToVert[2]).coord();
        const auto bounds = snapshot.nativeCellBounds(cell);
        check(bounds[0].x == lower.x && bounds[0].y == lower.y &&
            bounds[1].x == upper.x && bounds[1].y == upper.y,
            "native bounds mismatch at cell=" + std::to_string(cell));
    }
    auto copiedClosure = closure;
    auto ownedGeometry = std::make_shared<const LevelSetHermiteGeometry::Snapshot>(std::move(snapshot));
    data.closure = std::move(copiedClosure);
    data.geometry = std::move(ownedGeometry);
}
HermiteGalsState HermiteGalsCandidate::bindPreparedLocal(const HermiteGalsState& source,HermiteGalsTimeLayer time,
    const HermiteGalsPrimaryState& primary,const HermiteGalsClosure& closure,
    LevelSetHermiteGeometry::Snapshot&& snapshot,const std::vector<int>& actualCells,
    const std::vector<Vector2D>& nativePoints,const std::vector<unsigned char>& independent)
{
    namespace Exact=LAMP::Geometry::Detail::ExactExpansion;
    Exact::validateActualEnvironment();time.time();
    const auto& old=retained(source.data_);const auto& meshSource=source.geometry().snapshot();
    require(bool(old.mesh)&&nativePoints.size()==independent.size()
        &&nativePoints.size()==meshSource.vertexCount(),"prepared local mesh certificate mismatch");
    auto data=std::make_unique<Detail::HermiteGalsData>(old.mesh,time,primary);
    data->identity=HermiteGalsIdentity(nextIdentity());
    validatePrimaryOnTopology(static_cast<int>(nativePoints.size()),static_cast<int>(meshSource.nativeCellCount()),
        meshSource.minCoord().x,data->primary,data->identity,
        [&](int v){return nativePoints[v];},[&](int v){return independent[v]!=0;});
    require(primary.transport==HermiteGalsTransport::FirstJet&&primary.measure==HermiteGalsMeasure::CartesianArea
        &&primary.phaseVolumes.empty(),"prepared local requires pure XY FirstJet",data->identity);
    require(primary.geometryCells==actualCells,"prepared local activity/payload mismatch",data->identity);
    require(closure.vertexMixedDerivatives.size()==primary.vertices.size()
        &&closure.cellBubbleCoefficients.size()==primary.geometryCells.size(),"prepared local closure size",data->identity);
    for(double value:closure.vertexMixedDerivatives)require(std::isfinite(value),"non-finite mixed closure",data->identity);
    for(double value:closure.cellBubbleCoefficients)require(std::isfinite(value),"non-finite bubble closure",data->identity);
    require(snapshot.usesNativeLeaves()&&!snapshot.hasCornerSectors()
        &&snapshot.reconstruction()==LevelSetHermiteGeometry::Reconstruction::EnrichedQuartic
        &&snapshot.nativeCellCount()==meshSource.nativeCellCount(),"prepared local scalar support mismatch",data->identity);
    // Only the friend preparation can enter. It constructs this snapshot from
    // the SAME checked compact payload, so no arbitrary active mask is trusted.
    for(int cell:actualCells){
        require(snapshot.nativeCellActive(cell),"prepared local inactive native cell",data->identity);
        const auto a=snapshot.nativeCellBounds(cell),b=meshSource.nativeCellBounds(cell);
        require(a[0].x==b[0].x&&a[0].y==b[0].y&&a[1].x==b[1].x&&a[1].y==b[1].y,
            "prepared local native bounds mismatch",data->identity);
    }
    data->closure=closure;
    data->geometry=std::make_shared<const LevelSetHermiteGeometry::Snapshot>(std::move(snapshot));
    Exact::validateActualEnvironment();
    return HermiteGalsState(std::shared_ptr<const Detail::HermiteGalsData>(std::move(data)));
}
HermiteGalsState HermiteGalsCandidate::finish() &&
{
    const auto& data = retained(data_);
    require(static_cast<bool>(data.geometry), "finish: geometry/closure has not been bound", data.identity);
    return HermiteGalsState(std::shared_ptr<const Detail::HermiteGalsData>(std::move(data_)));
}
} // namespace LAMP::Interface
