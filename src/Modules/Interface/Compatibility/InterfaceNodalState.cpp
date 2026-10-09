#include "Modules/Interface/Compatibility/InterfaceNodalState.h"

#include "Core/Field/FieldChecks.h"
#include "Core/Geometry/ImmersedGeometry.h"
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Mesh/NativeMeshIdentity.h"
#include "Core/Mesh/CartesianLeafMesh.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <cmath>
#include <cstring>
#include <stdexcept>
#include <utility>

namespace
{
struct MotionLease { int step = 0; bool ready = false; };
void requireLease(const std::weak_ptr<const MotionLease>& lease, int step)
{
    const auto current = lease.lock();
    if (!current || !current->ready || current->step != step)
        throw std::logic_error("nodal material map no longer belongs to the current Interface substep");
}

InterfaceBoundaryGeometry::Report inspectNodal(const QuadtreeMesh& mesh,
    const ImmersedGeometry<QuadtreeMesh>& geometry, int interfaceIndex,
    const InterfaceNodalOptions& options, bool checkTrace)
{
    using FieldChecks::checkFinite;
    checkFinite(geometry.levelSets[interfaceIndex], "level set", mesh.vertIndexNum());
    checkFinite(geometry.LS_dx(interfaceIndex), "level-set dx", mesh.vertIndexNum());
    checkFinite(geometry.LS_dy(interfaceIndex), "level-set dy", mesh.vertIndexNum());
    checkFinite(geometry.LS_dxx(interfaceIndex), "level-set dxx", mesh.vertIndexNum());
    checkFinite(geometry.LS_dyy(interfaceIndex), "level-set dyy", mesh.vertIndexNum());
    checkFinite(geometry.LS_dxy(interfaceIndex), "level-set dxy", mesh.vertIndexNum());
    if (options.prescribedBoundary)
    {
        const auto& trace = geometry.traceDiagnostics();
        if (checkTrace && (trace.nonFiniteInterpolationCount != 0 || !std::isfinite(trace.maxTraceRatio)))
            throw std::runtime_error("level-set boundary advection encountered non-finite trace/interpolation values");
        return InterfaceBoundaryGeometry::inspect(mesh, geometry.levelSets[interfaceIndex],
            options.axisymmetric, options.fixedContacts);
    }
    InterfaceBoundaryGeometry::Report result;
    for (int i = 0; i < mesh.vertIndexNum(); ++i)
    {
        const auto point = mesh.verts(mesh.vertIndexToLabel[i]).coord();
        if ((point.x == mesh.minCoord.x || point.x == mesh.maxCoord.x ||
             point.y == mesh.minCoord.y || point.y == mesh.maxCoord.y) &&
            geometry.levelSets[interfaceIndex].value[i] <= 0.0)
        {
            if (options.allowAxisContact && point.x == 0.0 &&
                point.y > mesh.minCoord.y && point.y < mesh.maxCoord.y &&
                geometry.levelSets[interfaceIndex].value[i] < 0.0)
            {
                result.axisContact = true;
                continue;
            }
            throw std::runtime_error("fluid interface reached Domain boundary; contact-line/wall coupling is unsupported");
        }
    }
    const auto& trace = geometry.traceDiagnostics();
    if (checkTrace && (trace.nonFiniteInterpolationCount != 0 || !std::isfinite(trace.maxTraceRatio)))
        throw std::runtime_error("level-set advection encountered non-finite trace/interpolation values");
    result.cutCells = geometry.cutInfoDiagnostics().cutCellCount;
    if (result.cutCells <= 0)
        throw std::runtime_error("the selected fluid interface has no cut cells; it is unresolved or left the domain");
    return result;
}
}

struct InterfaceNodalState::Impl
{
    const QuadtreeMesh* mesh;
    std::shared_ptr<const NativeMeshIdentity> identity;
    std::unique_ptr<ImmersedGeometry<QuadtreeMesh>> geometry;
    int interfaceIndex;
    InterfaceNodalOptions options;
    InterfaceBoundaryGeometry::Report boundary;
    std::shared_ptr<MotionLease> lease = std::make_shared<MotionLease>();
    int completedSteps = 0;
    double dt = 0.0;
    bool failed = false;
    const char* stage = "geometry-initialize";

    void refresh(QuadtreeMesh& mesh, double time)
    {
        if (options.prescribedBoundary) geometry->refreshLevelSetDerivatives(mesh);
        else
        {
            geometry->refreshLevelSetDerivedData(mesh, time, false);
            geometry->checkCutInfo(mesh, true);
        }
    }
};

struct InterfaceNodalMaterialMap::Impl
{
    const QuadtreeMesh* mesh;
    const ImmersedGeometry<QuadtreeMesh>* geometry;
    std::shared_ptr<const NativeMeshIdentity> identity;
    std::weak_ptr<const MotionLease> lease;
    int interfaceIndex, step;
    double sourceTime, dt;
    bool axisContact;
    std::unique_ptr<const LevelSetHermiteGeometry::Snapshot> snapshot;
    std::function<Vector2D(Vector2D)> forward;
};

InterfaceNodalMaterialMap::InterfaceNodalMaterialMap(std::unique_ptr<Impl> impl) : impl_(std::move(impl)) {}
InterfaceNodalMaterialMap::~InterfaceNodalMaterialMap() = default;
InterfaceNodalMaterialMap::Motion InterfaceNodalMaterialMap::motionFor(
    const QuadtreeMesh& mesh, const ImmersedGeometry<QuadtreeMesh>* geometry,
    int interfaceIndex, const std::vector<double>& sourceScalar, int step,
    double sourceTime, double dt, bool axisContact) const
{
    const auto& m = *impl_;
    requireLease(m.lease, m.step);
    if (&mesh != m.mesh || geometry != m.geometry || interfaceIndex != m.interfaceIndex ||
        step != m.step || sourceTime != m.sourceTime || dt != m.dt || axisContact != m.axisContact)
        throw std::invalid_argument("nodal material map source, topology or time layer mismatch");
    m.identity->validate(mesh);
    // updateLevelSets saves phi^n in this existing history before writing
    // phi^(n+1). Redistance changes only the latter. The lease guards access
    // to both the owner and this particular history layer.
    const auto& previous = m.geometry->levelSetsPre[m.interfaceIndex];
    if (sourceScalar.size() != static_cast<std::size_t>(previous.size()) ||
        std::memcmp(sourceScalar.data(), previous.value.get(), sourceScalar.size()*sizeof(double)) != 0)
        throw std::invalid_argument("nodal material map does not match the consumer's source scalar");
    return {*m.snapshot, m.forward};
}

InterfaceNodalState::InterfaceNodalState(QuadtreeMesh& mesh,
    std::unique_ptr<ImmersedGeometry<QuadtreeMesh>> geometry, int interfaceIndex,
    InterfaceNodalOptions options) : impl_(std::make_unique<Impl>())
{
    if (!geometry || geometry->levelSets.size() != 1 || interfaceIndex != 0)
        throw std::invalid_argument("nodal Interface state requires one initialized liquid level set");
    if (options.reinitializeEvery < 0 ||
        (options.prescribedBoundary && options.reinitializeEvery != 0) ||
        (options.axisymmetric && mesh.minCoord.x < 0.0) ||
        (options.allowAxisContact && (!options.axisymmetric || mesh.minCoord.x != 0.0)) ||
        (!options.prescribedBoundary && !options.fixedContacts.empty()))
        throw std::invalid_argument("invalid nodal Interface boundary/redistance options");
    if (validateCartesianLeafMesh(mesh).mixedLevels)
        throw std::invalid_argument("nodal Interface state requires complete uniform leaves; use native Hermite transport for AMR");
    auto& s = *impl_;
    s.mesh = &mesh;
    s.identity = std::make_shared<const NativeMeshIdentity>(mesh);
    s.geometry = std::move(geometry);
    s.interfaceIndex = interfaceIndex;
    s.options = std::move(options);
    s.refresh(mesh, 0.0);
    s.boundary = inspectNodal(mesh, *s.geometry, interfaceIndex, s.options, false);
    s.lease->ready = true;
}
InterfaceNodalState::~InterfaceNodalState() = default;
const ImmersedGeometry<QuadtreeMesh>& InterfaceNodalState::geometry() const
{
    if (impl_->failed) throw std::logic_error("failed nodal Interface state has no usable geometry");
    return *impl_->geometry;
}
const InterfaceBoundaryGeometry::Report& InterfaceNodalState::boundaryGeometry() const
{
    if (impl_->failed) throw std::logic_error("failed nodal Interface state has no usable boundary report");
    return impl_->boundary;
}
const char* InterfaceNodalState::stage() const noexcept { return impl_->stage; }
bool InterfaceNodalState::failed() const noexcept { return impl_->failed; }

std::unique_ptr<InterfaceNodalMaterialMap> InterfaceNodalState::advance(QuadtreeMesh& mesh,
    int step, double sourceTime, double dt, InterfaceNodalVelocityLayers velocity, bool retainMaterialMap)
{
    auto& s = *impl_;
    if (s.failed) throw std::logic_error("failed nodal Interface state cannot advance again");
    s.stage = "geometry-advection";
    if (&mesh != s.mesh || step != s.completedSteps + 1 || !std::isfinite(dt) || !(dt > 0.0) ||
        !std::isfinite(sourceTime) || sourceTime != (step-1)*dt ||
        (s.completedSteps && dt != s.dt) || !s.geometry->levelSets[s.interfaceIndex].moving())
        throw std::invalid_argument("nodal Interface requires its next fixed-dt step on the original mesh");
    s.identity->validate(mesh);
    if (s.options.prescribedBoundary ? (!velocity.boundaryVelocity.current || !velocity.boundaryVelocity.previous)
        : (bool(velocity.boundaryVelocity.current) || bool(velocity.boundaryVelocity.previous)))
        throw std::invalid_argument("nodal Interface velocity boundary policy mismatch");
    if (retainMaterialMap && s.options.prescribedBoundary)
        throw std::invalid_argument("prescribedScalar material maps are not implemented");
    if (!velocity.currentX || !velocity.currentY || !velocity.previousX || !velocity.previousY)
        throw std::invalid_argument("nodal Interface requires both vertex velocity layers");
    using FieldChecks::checkFinite;
    checkFinite(*velocity.currentX, "advection velocity x", mesh.vertIndexNum());
    checkFinite(*velocity.currentY, "advection velocity y", mesh.vertIndexNum());
    checkFinite(*velocity.previousX, "advection history x", mesh.vertIndexNum());
    checkFinite(*velocity.previousY, "advection history y", mesh.vertIndexNum());

    // Core's existing updater mutates histories in place. Fail-stop begins
    // before the first such write; old material callbacks are invalid now.
    s.failed = true;
    s.lease->ready = false;
    s.geometry->updateLevelSets(mesh, *velocity.currentX, *velocity.currentY,
        *velocity.previousX, *velocity.previousY, dt,
        s.boundary.axisContact || (s.options.prescribedBoundary && s.options.axisymmetric),
        (step-1)*dt, s.options.prescribedBoundary ? &velocity.boundaryVelocity : nullptr);
    if (s.options.reinitializeEvery > 0 && step % s.options.reinitializeEvery == 0)
    {
        s.stage = "geometry-redistance";
        s.geometry->exactRedistanceMovingLevelSets(mesh, std::vector<int>{s.interfaceIndex});
    }
    s.stage = "geometry-derived-data";
    s.refresh(mesh, step*dt);
    const auto updated = inspectNodal(mesh, *s.geometry, s.interfaceIndex, s.options, true);
    if (updated.axisContact != s.boundary.axisContact)
        throw std::runtime_error("interface attachment to the radial axis changed; detachment/contact topology changes are unsupported");
    std::unique_ptr<InterfaceNodalMaterialMap> result;
    if (retainMaterialMap)
    {
        auto m = std::make_unique<InterfaceNodalMaterialMap::Impl>();
        m->mesh = &mesh; m->geometry = s.geometry.get(); m->identity = s.identity;
        m->lease = s.lease; m->interfaceIndex = s.interfaceIndex; m->step = step;
        m->sourceTime = sourceTime; m->dt = dt; m->axisContact = updated.axisContact;
        // Query geometry belongs to the completed Interface result. Ordinary
        // EHD formerly constructed this same Hermite snapshot in its adapter.
        // UniformLeaves also preserves this state's equal-level mesh scope.
        m->snapshot = std::make_unique<const LevelSetHermiteGeometry::Snapshot>(mesh,
            s.geometry->levelSets[s.interfaceIndex], s.options.axisymmetric && mesh.minCoord.x == 0.0,
            CartesianGridScope::UniformLeaves);
        const auto fields = std::make_shared<const InterfaceNodalVelocityLayers>(std::move(velocity));
        m->forward = [fields, lease=m->lease, meshPtr=&mesh, geometryPtr=s.geometry.get(),
            step, dt, axis=updated.axisContact](Vector2D point) {
            requireLease(lease, step);
            return geometryPtr->getNextCoord(*meshPtr, point, *fields->currentX, *fields->currentY,
                *fields->previousX, *fields->previousY, dt, axis);
        };
        result.reset(new InterfaceNodalMaterialMap(std::move(m)));
    }
    s.boundary = updated;
    s.completedSteps = step; s.dt = dt;
    s.lease->step = step; s.lease->ready = true;
    s.failed = false;
    return result;
}
