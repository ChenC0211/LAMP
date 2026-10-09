#include "Core/Field/StaggeredVelocitySampler.h"
#include "Core/Field/VertField.h"
#include "Core/Field/VEdgeField.h"
#include "Core/Field/HEdgeField.h"
#include "Core/Mesh/QuadtreeMesh.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace
{
void require(bool ok, const char* reason)
{
    if (!ok) throw std::runtime_error(std::string("[StaggeredVelocitySampler] ") + reason);
}
bool finite(Vector2D p) { return std::isfinite(p.x) && std::isfinite(p.y); }
bool metricNear(double actual, double expected)
{
    return std::isfinite(actual) && actual > 0.0
        && std::abs(static_cast<long double>(actual) - expected)
            <= 64.0L * std::numeric_limits<double>::epsilon() * expected;
}

struct Layout
{
    int nx, ny;
    Vector2D lower, upper, h;
    double tolerance;

    explicit Layout(const QuadtreeMesh& mesh)
        : nx(mesh.Nx), ny(mesh.Ny), lower(mesh.minCoord), upper(mesh.maxCoord), h(mesh.cellSize)
    {
        require(nx >= 2 && ny >= 2, "at least 2x2 cells are required for both native bilinear grids");
        require(finite(lower) && finite(upper) && finite(h) && h.x > 0.0 && h.y > 0.0
            && upper.x > lower.x && upper.y > lower.y, "invalid Cartesian domain or spacing");
        const long long count = static_cast<long long>(nx) * ny;
        const long long vertices = (static_cast<long long>(nx) + 1) * (static_cast<long long>(ny) + 1);
        const long long xFaces = (static_cast<long long>(nx) + 1) * ny;
        const long long yFaces = static_cast<long long>(nx) * (static_cast<long long>(ny) + 1);
        require(vertices <= std::numeric_limits<int>::max()
            && count == mesh.cellIndexNum() && vertices == mesh.vertIndexNum()
            && xFaces == mesh.vEdgeIndexNum() && yFaces == mesh.hEdgeIndexNum(),
            "complete fixed uniform Cartesian counts required; AMR is unsupported");
        require(mesh.cellIndexToLabel && mesh.vEdgeIndexToLabel && mesh.hEdgeIndexToLabel
            && mesh.cellLabelCapacity() >= count && mesh.vEdgeLabelCapacity() >= xFaces
            && mesh.hEdgeLabelCapacity() >= yFaces, "missing or undersized mesh maps/storage");
        tolerance = 256.0 * std::numeric_limits<double>::epsilon()
            * std::max({h.x, h.y, std::abs(lower.x), std::abs(lower.y), std::abs(upper.x), std::abs(upper.y)});
        require(std::isfinite(tolerance) && tolerance < std::min(h.x, h.y) / 16.0,
            "Cartesian spacing is unresolved at this coordinate scale");
        require(near(upper.x, static_cast<long double>(lower.x) + nx * static_cast<long double>(h.x))
            && near(upper.y, static_cast<long double>(lower.y) + ny * static_cast<long double>(h.y)),
            "domain extent disagrees with Cartesian count/spacing");
        std::vector<int> cellGrid(static_cast<std::size_t>(count), -1);
        for (int k = 0; k < mesh.cellIndexNum(); ++k)
        {
            const int label = mesh.cellIndexToLabel[k];
            require(label >= 0 && label < mesh.cellLabelCapacity(), "cell label outside storage");
            const auto& cell = mesh.cells(label);
            require(cell.index() == k && cell.label() == label, "cell solver-index/label mismatch");
            require(cell.level() == 0 && cell.parent == -1
                && std::all_of(cell.son, cell.son + 4, [](int child) { return child == -1; }),
                "only level-zero leaf cells are supported, not AMR/history parents");
            require(metricNear(cell.size().x, h.x) && metricNear(cell.size().y, h.y), "nonuniform cell dimensions");
            const auto p = cell.centerCoord();
            const int i = position(p.x, true, .5, nx), j = position(p.y, false, .5, ny);
            int& slot = cellGrid[i + nx * j];
            require(slot == -1, "duplicate Cartesian cell location");
            slot = k;
        }
    }
    bool near(double actual, long double expected) const
    {
        return std::isfinite(actual) && std::isfinite(expected)
            && std::abs(static_cast<long double>(actual) - expected) <= tolerance;
    }
    int position(double coordinate, bool x, double offset, int count) const
    {
        require(std::isfinite(coordinate), "non-finite native face coordinate");
        const double origin = x ? lower.x : lower.y, spacing = x ? h.x : h.y;
        const long double normalized = (static_cast<long double>(coordinate) - origin) / spacing - offset;
        const long double rounded = std::round(normalized);
        require(std::isfinite(normalized) && rounded >= 0.0L && rounded < count,
            "native coordinate outside the Cartesian lattice");
        const int result = static_cast<int>(rounded);
        require(near(coordinate, static_cast<long double>(origin)
            + (static_cast<long double>(result) + offset) * spacing), "nonuniform native face locations");
        return result;
    }
};

struct Component
{
    int nx = 0, ny = 0;
    std::vector<double> x, y, values; // Actual midpoint axes; values in solver-index order.
    std::vector<int> indices; // Cartesian lattice -> solver index.

    Component(const QuadtreeMesh& mesh, const Layout& layout, bool horizontal,
        const Field<double>& oldField, const Field<double>* newField)
    {
        nx = horizontal ? layout.nx : layout.nx + 1;
        ny = horizontal ? layout.ny + 1 : layout.ny;
        const int count = horizontal ? mesh.hEdgeIndexNum() : mesh.vEdgeIndexNum();
        require(oldField.value && oldField.size() == count
            && (!newField || (newField->value && newField->size() == count)),
            "native velocity field has missing storage or wrong size");
        const double missing = std::numeric_limits<double>::quiet_NaN();
        x.assign(nx, missing); y.assign(ny, missing);
        values.resize(count); indices.assign(count, -1);
        const auto coordinate = [](std::vector<double>& axis, int index, double value)
        {
            if (std::isnan(axis[index])) axis[index] = value;
            else require(axis[index] == value, "native face coordinates do not form one Cartesian product grid");
        };
        for (int k = 0; k < count; ++k)
        {
            const int label = horizontal ? mesh.hEdgeIndexToLabel[k] : mesh.vEdgeIndexToLabel[k];
            require(label >= 0 && label < (horizontal ? mesh.hEdgeLabelCapacity() : mesh.vEdgeLabelCapacity()),
                "native face label outside allocated storage");
            const auto& edge = horizontal ? mesh.hEdges(label) : mesh.vEdges(label);
            require(edge.index() == k && edge.label() == label && oldField.label(k) == label
                && (!newField || newField->label(k) == label), "native face/field solver-index/label mismatch");
            require(edge.isHorizontal() == horizontal && edge.parent == -1 && edge.son[0] == -1 && edge.son[1] == -1,
                "wrong face orientation or nonuniform/split native edge");
            require(metricNear(edge.length(), horizontal ? layout.h.x : layout.h.y), "nonuniform native edge length");
            const auto p = edge.midpointCoord();
            const int i = layout.position(p.x, true, horizontal ? .5 : 0.0, nx);
            const int j = layout.position(p.y, false, horizontal ? 0.0 : .5, ny);
            int& slot = indices[i + nx * j];
            require(slot == -1, "duplicate native Cartesian face location");
            slot = k;
            coordinate(x, i, p.x); coordinate(y, j, p.y);
            const double oldValue = oldField.value[k];
            require(std::isfinite(oldValue), "non-finite native velocity");
            if (!newField) values[k] = oldValue;
            else
            {
                const double newValue = newField->value[k];
                require(std::isfinite(newValue), "non-finite new native velocity");
                const double sum = oldValue + newValue;
                values[k] = .5 * sum;
                require(std::isfinite(sum) && std::isfinite(values[k]) && (sum == 0.0 || values[k] != 0.0),
                    "native midpoint velocity overflow or underflow");
            }
        }
        // Unique occupied positions plus matching counts imply a complete
        // native grid. These checks retain actual coordinates for interpolation
        // rather than silently replacing them with an idealized dense-label grid.
        const auto inspectAxis = [](const std::vector<double>& axis)
        {
            for (std::size_t k = 0; k < axis.size(); ++k)
            {
                require(std::isfinite(axis[k]), "incomplete native coordinate axis");
                // position() already checked each coordinate against its
                // lattice location at the coordinate's rounding scale. A
                // second spacing-relative check on their difference loses
                // that scale and rejects valid fine grids far from zero.
                if (k) require(axis[k] > axis[k - 1], "native coordinate axis is not strictly increasing");
            }
        };
        inspectAxis(x); inspectAxis(y);
    }

    void addBoundary(const Layout& layout, bool horizontal,
        const std::vector<double>& low, const std::vector<double>& high)
    {
        const int count = horizontal ? ny : nx;
        require(low.size() == static_cast<std::size_t>(count) && high.size() == low.size(),
            "tangential boundary trace has the wrong size");
        for (std::size_t k = 0; k < low.size(); ++k)
            require(std::isfinite(low[k]) && std::isfinite(high[k]), "non-finite tangential boundary trace");
        require(horizontal ? (y.front() == layout.lower.y && y.back() == layout.upper.y
                && x.front() > layout.lower.x && x.back() < layout.upper.x)
            : (x.front() == layout.lower.x && x.back() == layout.upper.x
                && y.front() > layout.lower.y && y.back() < layout.upper.y),
            "native normal boundary axes must equal the physical domain for explicit wall traces");
        const long long nextNx = static_cast<long long>(nx) + (horizontal ? 2 : 0);
        const long long nextNy = static_cast<long long>(ny) + (horizontal ? 0 : 2);
        require(nextNx * nextNy <= std::numeric_limits<int>::max(), "extended velocity grid exceeds index capacity");
        const int newNx = static_cast<int>(nextNx), newNy = static_cast<int>(nextNy);
        std::vector<int> extended(static_cast<std::size_t>(newNx) * newNy, -1);
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
                extended[i + (horizontal ? 1 : 0) + newNx * (j + (horizontal ? 0 : 1))] = indices[i + nx*j];
        const int lowBegin = static_cast<int>(values.size());
        values.insert(values.end(), low.begin(), low.end());
        const int highBegin = static_cast<int>(values.size());
        values.insert(values.end(), high.begin(), high.end());
        for (int k = 0; k < count; ++k)
        {
            extended[horizontal ? newNx*k : k] = lowBegin+k;
            extended[horizontal ? newNx-1 + newNx*k : k + newNx*(newNy-1)] = highBegin+k;
        }
        auto& axis = horizontal ? x : y;
        axis.insert(axis.begin(), horizontal ? layout.lower.x : layout.lower.y);
        axis.push_back(horizontal ? layout.upper.x : layout.upper.y);
        nx = newNx; ny = newNy;
        indices.swap(extended);
    }

    static int interval(const std::vector<double>& axis, double coordinate)
    {
        require(coordinate >= axis.front() && coordinate <= axis.back(), "point outside native velocity grid");
        if (coordinate == axis.back()) return static_cast<int>(axis.size()) - 2;
        const auto upper = std::upper_bound(axis.begin(), axis.end(), coordinate);
        const int result = static_cast<int>(upper - axis.begin()) - 1;
        require(result >= 0 && result + 1 < static_cast<int>(axis.size()), "native interpolation interval missing");
        return result;
    }
    static double lerp(double a, double b, double weight)
    {
        require(std::isfinite(weight) && weight >= 0.0 && weight <= 1.0, "invalid bilinear interpolation coordinate");
        if (weight == 0.0) return a;
        if (weight == 1.0) return b;
        if (a == b) return a;
        const long double value = (1.0L - weight) * a + static_cast<long double>(weight) * b;
        const double result = static_cast<double>(value);
        require(std::isfinite(value) && std::isfinite(result) && (value == 0.0L || result != 0.0),
            "sampled velocity overflow or underflow");
        return result;
    }
    double sample(Vector2D point) const
    {
        const int i = interval(x, point.x), j = interval(y, point.y);
        const double tx = (point.x - x[i]) / (x[i + 1] - x[i]);
        const double ty = (point.y - y[j]) / (y[j + 1] - y[j]);
        const double low = lerp(values[indices[i + nx * j]], values[indices[i + 1 + nx * j]], tx);
        const double high = lerp(values[indices[i + nx * (j + 1)]], values[indices[i + 1 + nx * (j + 1)]], tx);
        return lerp(low, high, ty);
    }
};

StaggeredVelocityBoundaryTrace retainedBoundary(const Layout& layout,
    const StaggeredVelocityBoundaryTrace& oldBoundary, const StaggeredVelocityBoundaryTrace* newBoundary)
{
    const auto retain = [](const std::vector<double>& oldValues, const std::vector<double>* newValues, std::size_t count)
    {
        require(oldValues.size() == count && (!newValues || newValues->size() == count),
            "tangential boundary trace has the wrong size");
        std::vector<double> result(count);
        for (std::size_t k = 0; k < count; ++k)
        {
            require(std::isfinite(oldValues[k]) && (!newValues || std::isfinite((*newValues)[k])),
                "non-finite tangential boundary trace");
            if (!newValues) result[k] = oldValues[k];
            else
            {
                const double sum = oldValues[k] + (*newValues)[k];
                result[k] = .5 * sum;
                require(std::isfinite(sum) && std::isfinite(result[k]) && (sum == 0.0 || result[k] != 0.0),
                    "boundary midpoint velocity overflow or underflow");
            }
        }
        return result;
    };
    return {retain(oldBoundary.uxSouth, newBoundary ? &newBoundary->uxSouth : nullptr, layout.nx+1),
        retain(oldBoundary.uxNorth, newBoundary ? &newBoundary->uxNorth : nullptr, layout.nx+1),
        retain(oldBoundary.uyWest, newBoundary ? &newBoundary->uyWest : nullptr, layout.ny+1),
        retain(oldBoundary.uyEast, newBoundary ? &newBoundary->uyEast : nullptr, layout.ny+1)};
}
}

struct StaggeredVelocitySampler::Impl
{
    Component ux, uy;
    Vector2D lower, upper;
    bool explicitBoundary = false;

    Impl(const QuadtreeMesh& mesh, const Layout& layout,
        const VEdgeField<double>& uxOld, const HEdgeField<double>& uyOld,
        const VEdgeField<double>* uxNew, const HEdgeField<double>* uyNew)
        : ux(mesh, layout, false, uxOld, uxNew), uy(mesh, layout, true, uyOld, uyNew),
          lower(std::max(ux.x.front(), uy.x.front()), std::max(ux.y.front(), uy.y.front())),
          upper(std::min(ux.x.back(), uy.x.back()), std::min(ux.y.back(), uy.y.back()))
    {
        require(lower.x < upper.x && lower.y < upper.y, "native component grids have no shared two-dimensional interior");
    }

    Impl(const QuadtreeMesh& mesh, const Layout& layout,
        const VEdgeField<double>& uxOld, const HEdgeField<double>& uyOld,
        const VEdgeField<double>* uxNew, const HEdgeField<double>* uyNew,
        const StaggeredVelocityBoundaryTrace& boundaryOld,
        const StaggeredVelocityBoundaryTrace* boundaryNew, bool radialAxisSymmetry)
        : Impl(mesh, layout, uxOld, uyOld, uxNew, uyNew)
    {
        const auto boundary = retainedBoundary(layout, boundaryOld, boundaryNew);
        if (radialAxisSymmetry)
        {
            require(layout.lower.x == 0.0, "radial-axis boundary sampling requires min.r=0");
            for (int j = 0; j < ux.ny; ++j)
            {
                const int k = ux.indices[ux.nx*j];
                require(uxOld.value[k] == 0.0 && (!uxNew || uxNew->value[k] == 0.0),
                    "each native time level must have exactly zero radial velocity on the symmetry axis");
            }
            require(boundaryOld.uxSouth.front() == 0.0 && boundaryOld.uxNorth.front() == 0.0 &&
                (!boundaryNew || (boundaryNew->uxSouth.front() == 0.0 && boundaryNew->uxNorth.front() == 0.0)),
                "each boundary time level must have exactly zero radial velocity at axis corners");
        }
        ux.addBoundary(layout, false, boundary.uxSouth, boundary.uxNorth);
        uy.addBoundary(layout, true, boundary.uyWest, boundary.uyEast);
        lower = layout.lower; upper = layout.upper;
        explicitBoundary = true;
    }
};

std::pair<std::unique_ptr<VertField<double>>, std::unique_ptr<VertField<double>>>
staggeredToVertexVelocity(const QuadtreeMesh& mesh,
    const VEdgeField<double>& ux, const HEdgeField<double>& uy, bool radialAxisSymmetry)
{
    if (!radialAxisSymmetry) return {ux.toVertField(mesh), uy.toVertField(mesh)};
    require(mesh.minCoord.x == 0.0, "radial-axis vertex conversion requires min.r=0");
    const Layout layout(mesh);
    const Component radial(mesh, layout, false, ux, nullptr);
    const Component axial(mesh, layout, true, uy, nullptr);
    for (int j = 0; j < radial.ny; ++j)
        require(radial.values[radial.indices[radial.nx*j]] == 0.0,
            "native radial velocity on the symmetry axis must be exactly zero");

    auto vx = ux.toVertField(mesh), vy = uy.toVertField(mesh);
    for (int k = 0; k < mesh.vertIndexNum(); ++k)
    {
        const auto point = mesh.verts(mesh.vertIndexToLabel[k]).coord();
        if (point.x != 0.0) continue;
        const int j = layout.position(point.y, false, 0.0, axial.ny);
        // u_z is even in r. Interpolate in r^2 using the first two actual
        // horizontal-face columns (h/2 and 3h/2), not a clipped query point.
        const long double a = axial.values[axial.indices[axial.nx*j]];
        const long double b = axial.values[axial.indices[1 + axial.nx*j]];
        const double axisValue = static_cast<double>((9.0L*a - b)/8.0L);
        require(std::isfinite(axisValue), "axis axial-velocity reconstruction overflow");
        vx->value[k] = 0.0;
        vy->value[k] = axisValue;
    }
    return {std::move(vx), std::move(vy)};
}

std::pair<std::unique_ptr<VertField<double>>, std::unique_ptr<VertField<double>>>
staggeredToVertexVelocity(const QuadtreeMesh& mesh,
    const VEdgeField<double>& ux, const HEdgeField<double>& uy,
    const StaggeredVelocityBoundaryTrace& boundary, bool radialAxisSymmetry)
{
    const Layout layout(mesh);
    const StaggeredVelocitySampler sampler(mesh, ux, uy, boundary, radialAxisSymmetry);
    require(mesh.vertIndexToLabel && mesh.vertLabelCapacity() >= mesh.vertIndexNum(), "missing vertex conversion map/storage");
    std::vector<bool> occupied(static_cast<std::size_t>(mesh.vertIndexNum()), false);
    std::vector<Vector2D> values(mesh.vertIndexNum());
    for (int k = 0; k < mesh.vertIndexNum(); ++k)
    {
        const int label = mesh.vertIndexToLabel[k];
        require(label >= 0 && label < mesh.vertLabelCapacity() && mesh.verts(label).index() == k
            && mesh.verts(label).label() == label, "invalid vertex conversion index/label");
        const auto point = mesh.verts(label).coord();
        const int i = layout.position(point.x, true, 0.0, layout.nx+1);
        const int j = layout.position(point.y, false, 0.0, layout.ny+1);
        const std::size_t slot = static_cast<std::size_t>(i) + static_cast<std::size_t>(layout.nx+1)*j;
        require(!occupied[slot], "duplicate Cartesian vertex in velocity conversion");
        occupied[slot] = true;
        values[k] = sampler.sample(point);
    }
    auto vx = std::make_unique<VertField<double>>(), vy = std::make_unique<VertField<double>>();
    vx->initialize(mesh, 0.0); vy->initialize(mesh, 0.0);
    for (int k = 0; k < mesh.vertIndexNum(); ++k)
    { vx->value[k] = values[k].x; vy->value[k] = values[k].y; }
    return {std::move(vx), std::move(vy)};
}

StaggeredVelocitySampler::StaggeredVelocitySampler(const QuadtreeMesh& mesh,
    const VEdgeField<double>& ux, const HEdgeField<double>& uy)
    : impl_(std::make_unique<Impl>(mesh, Layout(mesh), ux, uy, nullptr, nullptr)) {}

StaggeredVelocitySampler::StaggeredVelocitySampler(const QuadtreeMesh& mesh,
    const VEdgeField<double>& uxOld, const HEdgeField<double>& uyOld,
    const VEdgeField<double>& uxNew, const HEdgeField<double>& uyNew)
    : impl_(std::make_unique<Impl>(mesh, Layout(mesh), uxOld, uyOld, &uxNew, &uyNew)) {}

StaggeredVelocitySampler::StaggeredVelocitySampler(const QuadtreeMesh& mesh,
    const VEdgeField<double>& ux, const HEdgeField<double>& uy,
    const StaggeredVelocityBoundaryTrace& boundary, bool radialAxisSymmetry)
    : impl_(std::make_unique<Impl>(mesh, Layout(mesh), ux, uy, nullptr, nullptr, boundary, nullptr, radialAxisSymmetry)) {}

StaggeredVelocitySampler::StaggeredVelocitySampler(const QuadtreeMesh& mesh,
    const VEdgeField<double>& uxOld, const HEdgeField<double>& uyOld,
    const VEdgeField<double>& uxNew, const HEdgeField<double>& uyNew,
    const StaggeredVelocityBoundaryTrace& boundaryOld,
    const StaggeredVelocityBoundaryTrace& boundaryNew, bool radialAxisSymmetry)
    : impl_(std::make_unique<Impl>(mesh, Layout(mesh), uxOld, uyOld, &uxNew, &uyNew,
        boundaryOld, &boundaryNew, radialAxisSymmetry)) {}

StaggeredVelocitySampler::~StaggeredVelocitySampler() = default;
StaggeredVelocitySampler::StaggeredVelocitySampler(StaggeredVelocitySampler&&) noexcept = default;
StaggeredVelocitySampler& StaggeredVelocitySampler::operator=(StaggeredVelocitySampler&&) noexcept = default;

Vector2D StaggeredVelocitySampler::sample(Vector2D point) const
{
    require(static_cast<bool>(impl_), "moved-from sampler");
    require(finite(point), "non-finite velocity query point");
    require(point.x >= impl_->lower.x && point.x <= impl_->upper.x
        && point.y >= impl_->lower.y && point.y <= impl_->upper.y,
        impl_->explicitBoundary ? "point outside the closed physical domain; boundary clipping/extrapolation unsupported"
            : "point outside the two native grids' common domain; clipping/extrapolation unsupported");
    return {impl_->ux.sample(point), impl_->uy.sample(point)};
}
