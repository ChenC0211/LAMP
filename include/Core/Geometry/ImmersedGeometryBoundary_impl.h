#pragma once

#include "Core/Geometry/LevelSetBoundaryTrace.h"
#include <iomanip>
#include <sstream>

namespace ImmersedGeometryBoundaryDetail
{
inline bool same(Vector2D a, Vector2D b) { return a.x == b.x && a.y == b.y; }
inline bool finite(Vector2D p) { return std::isfinite(p.x) && std::isfinite(p.y); }
inline void require(bool value, const char* reason)
{ if (!value) throw std::runtime_error(std::string("[ImmersedGeometry boundary] ") + reason); }
inline bool onSide(Vector2D p, Boundary side, Vector2D lower, Vector2D upper)
{
    switch (side)
    {
    case Boundary::west: return p.x == lower.x;
    case Boundary::south: return p.y == lower.y;
    case Boundary::east: return p.x == upper.x;
    case Boundary::north: return p.y == upper.y;
    default: return false;
    }
}
inline double tangentCoordinate(Vector2D p, Boundary side)
{ return side == Boundary::west || side == Boundary::east ? p.y : p.x; }
inline double outwardVelocity(Vector2D u, Boundary side)
{
    switch (side)
    {
    case Boundary::west: return -u.x;
    case Boundary::south: return -u.y;
    case Boundary::east: return u.x;
    case Boundary::north: return u.y;
    default: throw std::runtime_error("[ImmersedGeometry boundary] invalid side");
    }
}
}

template<typename Mesh>
bool ImmersedGeometry<Mesh>::hasLevelSetBoundaryConditions() const
{
    return std::any_of(levelSetBoundaryConditions_.begin(), levelSetBoundaryConditions_.end(),
        [](const auto& boundary) { return !boundary.empty(); });
}

template<typename Mesh>
void ImmersedGeometry<Mesh>::validateBoundaryMesh(const Mesh& mesh) const
{
    using namespace ImmersedGeometryBoundaryDetail;
    require(!periodicGeometry_, "explicit scalar boundaries do not support periodic geometry");
    require(mesh.Nx >= 4 && mesh.Ny >= 4 && finite(mesh.minCoord) && finite(mesh.maxCoord)
        && finite(mesh.cellSize) && mesh.cellSize.x > 0. && mesh.cellSize.y > 0.
        && mesh.minCoord.x < mesh.maxCoord.x && mesh.minCoord.y < mesh.maxCoord.y,
        "invalid fixed Cartesian layout");
    require(static_cast<long long>(mesh.Nx)*mesh.Ny == mesh.cellIndexNum()
        && (static_cast<long long>(mesh.Nx)+1)*(static_cast<long long>(mesh.Ny)+1) == mesh.vertIndexNum()
        && mesh.cellIndexToLabel && mesh.vertIndexToLabel, "only complete fixed uniform meshes are supported");
    for (int k = 0; k < mesh.cellIndexNum(); ++k)
    {
        const int label = mesh.cellIndexToLabel[k];
        require(ImmersedGeometryDetail::validCellLabel(mesh, label), "invalid boundary mesh cell label");
        const auto& c = mesh.cells(label);
        require(c.index() == k && c.label() == label && c.level() == 0 && c.parent == -1
            && std::all_of(c.son, c.son+4, [](int n) { return n == -1; })
            && same(c.size(), mesh.cellSize) && finite(c.centerCoord())
            && ImmersedGeometryDetail::validCellForInterpolation(mesh, c),
            "AMR/split or inconsistent boundary mesh cell");
        if (boundaryMeshLayout_)
        {
            require(boundaryMeshLayout_->cells.size() == static_cast<std::size_t>(mesh.cellIndexNum()),
                "cell count changed after boundary registration");
            const auto& saved = boundaryMeshLayout_->cells[k];
            require(saved.label == label && same(saved.center, c.centerCoord())
                && std::equal(saved.vertices.begin(), saved.vertices.end(), c.cellToVert),
                "cell location/topology changed after boundary registration");
        }
    }
    if (boundaryMeshLayout_)
    {
        const auto& layout = *boundaryMeshLayout_;
        require(layout.nx == mesh.Nx && layout.ny == mesh.Ny && same(layout.lower, mesh.minCoord)
            && same(layout.upper, mesh.maxCoord) && same(layout.spacing, mesh.cellSize)
            && layout.vertices.size() == static_cast<std::size_t>(mesh.vertIndexNum()),
            "mesh changed after scalar boundary registration");
    }
    for (int k = 0; k < mesh.vertIndexNum(); ++k)
    {
        const int label = mesh.vertIndexToLabel[k];
        require(label >= 0 && label < mesh.vertLabelCapacity() && mesh.verts(label).index() == k,
            "invalid boundary mesh vertex index/label");
        const auto p = mesh.verts(label).coord();
        require(finite(p) && p.x >= mesh.minCoord.x && p.x <= mesh.maxCoord.x
            && p.y >= mesh.minCoord.y && p.y <= mesh.maxCoord.y, "invalid boundary mesh vertex coordinate");
        if (boundaryMeshLayout_)
            require(same(p, boundaryMeshLayout_->vertices[k]), "vertex ordering/coordinates changed after boundary registration");
    }
}

template<typename Mesh>
void ImmersedGeometry<Mesh>::setLevelSetBoundaryConditions(const Mesh& mesh, int levelSetI,
    LevelSetBoundaryConditions conditions)
{
    using namespace ImmersedGeometryBoundaryDetail;
    require(levelSetI >= 0 && static_cast<std::size_t>(levelSetI) < levelSets.size(), "invalid boundary level-set index");
    require(!boundaryHistoryInitialized_ && levelSetBoundaryConditions_[levelSetI].empty(),
        "scalar boundary configuration is immutable once registered/initialized");
    require(!conditions.empty(), "omit boundary registration for the legacy path");
    validateBoundaryMesh(mesh);
    for (std::size_t k = 0; k < conditions.inflow.size(); ++k)
    {
        const auto& patch = conditions.inflow[k];
        require(patch.side == Boundary::west || patch.side == Boundary::south
            || patch.side == Boundary::east || patch.side == Boundary::north, "inflow patch has invalid side");
        require(static_cast<bool>(patch.value) && std::isfinite(patch.start) && std::isfinite(patch.end)
            && patch.start < patch.end && patch.start >= tangentCoordinate(mesh.minCoord, patch.side)
            && patch.end <= tangentCoordinate(mesh.maxCoord, patch.side), "invalid inflow interval or missing value callback");
        for (std::size_t j = 0; j < k; ++j)
            require(patch.side != conditions.inflow[j].side ||
                std::max(patch.start, conditions.inflow[j].start) >= std::min(patch.end, conditions.inflow[j].end),
                "inflow patches overlap in their interiors");
    }
    for (std::size_t k = 0; k < conditions.fixedZeroVertices.size(); ++k)
    {
        const auto p = conditions.fixedZeroVertices[k];
        require(finite(p) && (onSide(p, Boundary::west, mesh.minCoord, mesh.maxCoord)
            || onSide(p, Boundary::south, mesh.minCoord, mesh.maxCoord)
            || onSide(p, Boundary::east, mesh.minCoord, mesh.maxCoord)
            || onSide(p, Boundary::north, mesh.minCoord, mesh.maxCoord)), "fixed zero must lie on a physical boundary");
        int matches = 0;
        for (int i = 0; i < mesh.vertIndexNum(); ++i)
            if (same(p, mesh.verts(mesh.vertIndexToLabel[i]).coord())) ++matches;
        require(matches == 1, "fixed zero must identify one exact mesh vertex; snapping is unsupported");
        for (std::size_t j = 0; j < k; ++j)
            require(!same(p, conditions.fixedZeroVertices[j]), "duplicate fixed-zero vertex");
    }
    if (!boundaryMeshLayout_)
    {
        BoundaryMeshLayout layout;
        layout.lower = mesh.minCoord; layout.upper = mesh.maxCoord; layout.spacing = mesh.cellSize;
        layout.nx = mesh.Nx; layout.ny = mesh.Ny;
        for (int k = 0; k < mesh.vertIndexNum(); ++k) layout.vertices.push_back(mesh.verts(mesh.vertIndexToLabel[k]).coord());
        for (int k = 0; k < mesh.cellIndexNum(); ++k)
        {
            const int label = mesh.cellIndexToLabel[k];
            const auto& cell = mesh.cells(label);
            layout.cells.push_back({label, cell.centerCoord(),
                {cell.cellToVert[0], cell.cellToVert[1], cell.cellToVert[2], cell.cellToVert[3]}});
        }
        boundaryMeshLayout_ = std::move(layout);
    }
    levelSetBoundaryConditions_[levelSetI] = std::move(conditions);
}

template<typename Mesh>
bool ImmersedGeometry<Mesh>::prescribedBoundaryVertex(const Mesh& mesh, int levelSetI, Vector2D point) const
{
    using namespace ImmersedGeometryBoundaryDetail;
    const auto& boundary = levelSetBoundaryConditions_[levelSetI];
    for (const auto& p : boundary.fixedZeroVertices) if (same(p, point)) return true;
    for (const auto& patch : boundary.inflow)
        if (onSide(point, patch.side, mesh.minCoord, mesh.maxCoord)
            && tangentCoordinate(point, patch.side) >= patch.start && tangentCoordinate(point, patch.side) <= patch.end) return true;
    return false;
}

template<typename Mesh>
std::optional<double> ImmersedGeometry<Mesh>::prescribedBoundaryValue(const Mesh& mesh, int levelSetI,
    Vector2D point, double time, std::optional<Boundary> side) const
{
    using namespace ImmersedGeometryBoundaryDetail;
    require(std::isfinite(time), "boundary values require explicit finite physical time");
    std::optional<double> result;
    for (const auto& patch : levelSetBoundaryConditions_[levelSetI].inflow)
    {
        if ((side && patch.side != *side) || !onSide(point, patch.side, mesh.minCoord, mesh.maxCoord)
            || tangentCoordinate(point, patch.side) < patch.start || tangentCoordinate(point, patch.side) > patch.end) continue;
        const double value = patch.value(point, time);
        require(std::isfinite(value), "inflow value callback returned a nonfinite scalar");
        require(!result || *result == value, "inflow values disagree at a shared patch endpoint/corner");
        result = value;
    }
    for (const auto& p : levelSetBoundaryConditions_[levelSetI].fixedZeroVertices) if (same(p, point))
    {
        require(!result || *result == 0., "inflow profile conflicts with a fixed-zero vertex");
        result = 0.;
    }
    return result;
}

template<typename Mesh>
void ImmersedGeometry<Mesh>::refreshLevelSetDerivatives(const Mesh& mesh)
{
    if (hasLevelSetBoundaryConditions())
    {
        validateBoundaryMesh(mesh);
        if (boundaryHistoryInitialized_)
            for (std::size_t i = 0; i < levelSets.size(); ++i)
                ImmersedGeometryBoundaryDetail::require(levelSetsPre[i].size() == mesh.vertIndexNum()
                    && levelSetsPrePre[i].size() == mesh.vertIndexNum()
                    && levelSetsPre_dx[i].size() == mesh.vertIndexNum() && levelSetsPre_dy[i].size() == mesh.vertIndexNum()
                    && levelSetsPrePre_dx[i].size() == mesh.vertIndexNum() && levelSetsPrePre_dy[i].size() == mesh.vertIndexNum(),
                    "boundary history resized; refusing to replace old time levels by current values");
    }
    // Validate all scalar inputs before derivative/history storage is touched.
    for (int i = 0; i < static_cast<int>(levelSets.size()); ++i)
    {
        ImmersedGeometryBoundaryDetail::require(levelSets[i].value && levelSets[i].size() == mesh.vertIndexNum(),
            "scalar derivative refresh has invalid field storage");
        for (int k = 0; k < mesh.vertIndexNum(); ++k)
        {
            ImmersedGeometryBoundaryDetail::require(std::isfinite(levelSets[i].value[k]), "scalar derivative input is nonfinite");
            for (const auto& zero : levelSetBoundaryConditions_[i].fixedZeroVertices)
                if (ImmersedGeometryBoundaryDetail::same(zero, mesh.verts(mesh.vertIndexToLabel[k]).coord()))
                    ImmersedGeometryBoundaryDetail::require(levelSets[i].value[k] == 0., "fixed-zero node changed before scalar derivative refresh");
        }
    }
    for (int i = 0; i < static_cast<int>(levelSets.size()); ++i)
    {
        prepareLevelSetDerivedStorage(mesh, i);
        computeSecondDerivatives(mesh, i);
        computeFirstDerivativesAndMixed(mesh, i);
        syncHistoricalLevelSetStorageToCurrent(mesh, i);
        summarizeLevelSetDiagnostics(mesh, i);
    }
    if (hasLevelSetBoundaryConditions()) boundaryHistoryInitialized_ = true;
}

template<typename Mesh>
void ImmersedGeometry<Mesh>::stageBoundaryLevelSet(const Mesh& mesh, int levelSetI,
    const VertField<double>& Vx, const VertField<double>& Vy,
    const VertField<double>& preVx, const VertField<double>& preVy,
    double dt, bool radialAxisSymmetry, double previousTime,
    std::vector<double>& values, std::vector<TraceQueryDiagnostics>& diagnostics,
    const LevelSetAdvectionVelocity* advectionVelocity, const LevelSetAdvectionScalar* advectionScalar) const
{
    using namespace ImmersedGeometryBoundaryDetail;
    require(boundaryHistoryInitialized_, "initialize scalar derivatives/history before boundary advancement");
    require(std::isfinite(dt) && dt >= 0. && std::isfinite(previousTime) && std::isfinite(previousTime+dt)
        && (dt == 0. || previousTime+dt > previousTime), "boundary advancement requires a resolved finite time interval");
    require(levelSets[levelSetI].value && levelSets[levelSetI].size() == mesh.vertIndexNum(), "boundary scalar storage mismatch");
    for (const auto* field : {&Vx, &Vy, &preVx, &preVy})
        require(field->value && field->size() == mesh.vertIndexNum(), "boundary velocity storage mismatch");
    for (const auto* field : {&levelSetsPre[levelSetI], &levelSetsPrePre[levelSetI],
        &levelSetsPre_dx[levelSetI], &levelSetsPre_dy[levelSetI], &levelSetsPrePre_dx[levelSetI], &levelSetsPrePre_dy[levelSetI],
        &LSs_dx_[levelSetI], &LSs_dy_[levelSetI]})
        require(field->value && field->size() == mesh.vertIndexNum(), "boundary history/derivative storage mismatch");
    for (const auto& patch : levelSetBoundaryConditions_[levelSetI].inflow)
        require(!radialAxisSymmetry || patch.side != Boundary::west, "radial symmetry axis cannot be an inflow patch");
    const auto velocity = [&](Vector2D point, const VertField<double>& x, const VertField<double>& y)
    {
        require(finite(point) && point.x >= mesh.minCoord.x && point.x <= mesh.maxCoord.x
            && point.y >= mesh.minCoord.y && point.y <= mesh.maxCoord.y, "boundary trace attempted exterior velocity sampling");
        const int label = mesh.findCellLabel(point);
        require(ImmersedGeometryDetail::validCellLabel(mesh, label)
            && ImmersedGeometryDetail::validCellForInterpolation(mesh, mesh.cells(label)), "boundary trace lacks an interpolation cell");
        return Vector2D{Interpolation::biQuadInterpolation(mesh, x, mesh.cells(label), point),
            Interpolation::biQuadInterpolation(mesh, y, mesh.cells(label), point)};
    };
    const auto newVelocity = [&](Vector2D point) {
        return advectionVelocity ? advectionVelocity->current(point) : velocity(point, Vx, Vy);
    };
    const auto oldVelocity = [&](Vector2D point) {
        return advectionVelocity ? advectionVelocity->previous(point) : velocity(point, preVx, preVy);
    };
    values.resize(mesh.vertIndexNum()); diagnostics.resize(mesh.vertIndexNum());
    const LevelSetBoundaryTrace::Domain domain{mesh.minCoord, mesh.maxCoord, mesh.cellSize, radialAxisSymmetry};
    for (int i = 0; i < mesh.vertIndexNum(); ++i)
    {
        const auto point = mesh.verts(mesh.vertIndexToLabel[i]).coord();
        require(std::isfinite(levelSets[levelSetI].value[i]), "boundary update input scalar is nonfinite");
        for (const auto& zero : levelSetBoundaryConditions_[levelSetI].fixedZeroVertices)
            if (same(point, zero)) require(levelSets[levelSetI].value[i] == 0., "a declared fixed-zero node was externally changed");
        const auto prescribed = prescribedBoundaryValue(mesh, levelSetI, point, previousTime+dt);
        if (prescribed)
        {
            const auto unew = newVelocity(point), uold = oldVelocity(point);
            require(finite(unew) && finite(uold) && std::hypot(unew.x, unew.y)*dt <= .5*std::min(mesh.cellSize.x, mesh.cellSize.y)
                && std::hypot(uold.x, uold.y)*dt <= .5*std::min(mesh.cellSize.x, mesh.cellSize.y),
                "prescribed-boundary velocity is nonfinite or exceeds the half-cell trace limit");
            for (const auto& patch : levelSetBoundaryConditions_[levelSetI].inflow)
                if (onSide(point, patch.side, mesh.minCoord, mesh.maxCoord)
                    && tangentCoordinate(point, patch.side) >= patch.start && tangentCoordinate(point, patch.side) <= patch.end)
                    require(finite(unew) && finite(uold) && outwardVelocity(unew, patch.side) <= 0.
                        && outwardVelocity(uold, patch.side) <= 0., "inflow-value patch has reversed/outflow velocity");
            values[i] = *prescribed;
            diagnostics[i].prescribedBoundaryVertexCount = 1;
            continue;
        }
        const auto trace = LevelSetBoundaryTrace::trace(domain, point, dt, newVelocity, oldVelocity);
        if (trace.boundaryHit)
        {
            const auto value = prescribedBoundaryValue(mesh, levelSetI, trace.point, previousTime+dt-trace.elapsed, trace.side);
            if (!value)
            {
                const auto unew = newVelocity(trace.point), uold = oldVelocity(trace.point);
                std::ostringstream reason;
                reason << std::setprecision(17) << "[ImmersedGeometry boundary] level-set backtrace from ("
                    << point.x << ',' << point.y << ") leaves through an undeclared inflow at ("
                    << trace.point.x << ',' << trace.point.y << "), elapsed=" << trace.elapsed
                    << ", newVelocity=(" << unew.x << ',' << unew.y << "), oldVelocity=("
                    << uold.x << ',' << uold.y << ')';
                throw std::runtime_error(reason.str());
            }
            values[i] = *value;
            diagnostics[i].inflowBoundaryCount = 1;
        }
        else
        {
            if (advectionScalar) values[i] = advectionScalar->previous(levelSetI, trace.point);
            else
            {
                const int label = mesh.findCellLabel(trace.point);
                require(ImmersedGeometryDetail::validCellLabel(mesh, label)
                    && ImmersedGeometryDetail::validCellForInterpolation(mesh, mesh.cells(label)), "boundary scalar trace lacks an interpolation cell");
                values[i] = Interpolation::biQuadInterpolation(mesh, levelSets[levelSetI], mesh.cells(label), trace.point);
            }
        }
        require(std::isfinite(values[i]), "boundary trace produced a nonfinite scalar");
        diagnostics[i].traceRatio = std::hypot(point.x-trace.point.x, point.y-trace.point.y)/std::min(mesh.cellSize.x, mesh.cellSize.y);
    }
}
