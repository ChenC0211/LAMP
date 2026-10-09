#include "Modules/Hydrodynamics/HydroBoundaryVelocity.h"
#include "Modules/Hydrodynamics/HydroBoundary.h"
#include "Core/Field/HEdgeField.h"
#include "Core/Field/VEdgeField.h"
#include "Core/Field/VertField.h"
#include "Core/Mesh/QuadtreeMesh.h"

#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
void require(bool ok, const char* reason)
{
    if (!ok) throw std::invalid_argument(std::string("[HydroBoundaryVelocity] ") + reason);
}
}

StaggeredVelocityBoundaryTrace hydroZeroShearBoundaryTrace(const QuadtreeMesh& mesh,
    const VEdgeField<double>& ux, const HEdgeField<double>& uy, bool radialAxisSymmetry)
{
    // Reuse Core's complete uniform mesh, coordinate, label, finite-field and
    // storage validation. The maps below only locate real native values; they
    // do not introduce a second spatial interpolation or assumed label order.
    const StaggeredVelocitySampler nativeValidation(mesh, ux, uy);
    (void)nativeValidation;
    const int nx = mesh.Nx, ny = mesh.Ny;
    const double hx = mesh.cellSize.x, hy = mesh.cellSize.y;
    require(!radialAxisSymmetry || mesh.minCoord.x == 0.0,
        "radial-axis wall conversion requires min.r=0");
    std::vector<int> uIndex(mesh.vEdgeIndexNum(), -1), vIndex(mesh.hEdgeIndexNum(), -1);
    for (int k = 0; k < mesh.vEdgeIndexNum(); ++k)
    {
        const auto p = mesh.vEdges(mesh.vEdgeIndexToLabel[k]).midpointCoord();
        const int i = static_cast<int>(std::llround((static_cast<long double>(p.x)-mesh.minCoord.x)/hx));
        const int j = static_cast<int>(std::llround((static_cast<long double>(p.y)-mesh.minCoord.y)/hy-.5L));
        require(i >= 0 && i <= nx && j >= 0 && j < ny, "radial/x face outside validated native grid");
        require((i != 0 || p.x == mesh.minCoord.x) && (i != nx || p.x == mesh.maxCoord.x),
            "normal x boundary faces are not on the physical walls");
        int& index = uIndex[i+(nx+1)*j];
        require(index == -1, "duplicate radial/x native face");
        index = k;
        if (radialAxisSymmetry && i == 0)
            require(ux.value[k] == 0.0, "native radial velocity on the symmetry axis must be exactly zero");
    }
    for (int k = 0; k < mesh.hEdgeIndexNum(); ++k)
    {
        const auto p = mesh.hEdges(mesh.hEdgeIndexToLabel[k]).midpointCoord();
        const int i = static_cast<int>(std::llround((static_cast<long double>(p.x)-mesh.minCoord.x)/hx-.5L));
        const int j = static_cast<int>(std::llround((static_cast<long double>(p.y)-mesh.minCoord.y)/hy));
        require(i >= 0 && i < nx && j >= 0 && j <= ny, "axial/y face outside validated native grid");
        require((j != 0 || p.y == mesh.minCoord.y) && (j != ny || p.y == mesh.maxCoord.y),
            "normal y boundary faces are not on the physical walls");
        int& index = vIndex[i+nx*j];
        require(index == -1, "duplicate axial/y native face");
        index = k;
    }
    const auto u = [&](int i, int j) { return ux.value[uIndex[i+(nx+1)*j]]; };
    const auto v = [&](int i, int j) { return uy.value[vIndex[i+nx*j]]; };
    StaggeredVelocityBoundaryTrace result;
    result.uxSouth.resize(nx+1); result.uxNorth.resize(nx+1);
    result.uyWest.resize(ny+1); result.uyEast.resize(ny+1);
    const auto wallTrace = [](double interior, double hn, double sign, double difference, double ht, bool corner) {
        // The corner donor in Hydro first evaluates its one-sided derivative.
        // Non-corner tangential rows retain their multiply-then-divide order.
        return corner ? HydroBoundaryDetail::zeroShearTrace(interior, hn, sign, difference/ht)
            : HydroBoundaryDetail::zeroShearTrace(interior, hn, sign, difference, ht);
    };
    for (int i = 0; i <= nx; ++i)
    {
        if (radialAxisSymmetry && i == 0)
        {
            result.uxSouth[i] = result.uxNorth[i] = 0.0;
            continue;
        }
        const bool corner = i == 0 || i == nx;
        const int low = i == 0 ? 0 : i == nx ? nx-2 : i-1;
        const int high = low+1;
        result.uxSouth[i] = wallTrace(u(i, 0), hy, -1.0, v(high, 0)-v(low, 0), hx, corner);
        result.uxNorth[i] = wallTrace(u(i, ny-1), hy, 1.0, v(high, ny)-v(low, ny), hx, corner);
        require(std::isfinite(result.uxSouth[i]) && std::isfinite(result.uxNorth[i]),
            "nonfinite horizontal-wall tangential trace");
    }
    for (int j = 0; j <= ny; ++j)
    {
        const bool corner = j == 0 || j == ny;
        const int low = j == 0 ? 0 : j == ny ? ny-2 : j-1;
        const int high = low+1;
        if (radialAxisSymmetry)
        {
            // Same even-r^2 axial trace as staggeredToVertexVelocity's axis
            // path. A zero-area symmetry axis is not a physical slip wall.
            const long double a = v(0, j), b = v(1, j);
            result.uyWest[j] = static_cast<double>((9.0L*a-b)/8.0L);
        }
        else result.uyWest[j] = wallTrace(v(0, j), hx, -1.0, u(0, high)-u(0, low), hy, corner);
        result.uyEast[j] = wallTrace(v(nx-1, j), hx, 1.0, u(nx, high)-u(nx, low), hy, corner);
        require(std::isfinite(result.uyWest[j]) && std::isfinite(result.uyEast[j]),
            "nonfinite vertical-wall tangential trace");
    }
    return result;
}

std::pair<std::unique_ptr<VertField<double>>, std::unique_ptr<VertField<double>>>
hydroVertexVelocity(const QuadtreeMesh& mesh,
    const VEdgeField<double>& ux, const HEdgeField<double>& uy,
    bool explicitBoundaryPatches, bool radialAxisSymmetry)
{
    if (!explicitBoundaryPatches) return staggeredToVertexVelocity(mesh, ux, uy, radialAxisSymmetry);
    const auto boundary = hydroZeroShearBoundaryTrace(mesh, ux, uy, radialAxisSymmetry);
    return staggeredToVertexVelocity(mesh, ux, uy, boundary, radialAxisSymmetry);
}

std::function<Vector2D(Vector2D)> hydroBoundaryAdvectionVelocity(const QuadtreeMesh& mesh,
    const VEdgeField<double>& ux, const HEdgeField<double>& uy,
    const std::vector<HydroBoundaryPatch>& patches, bool radialAxisSymmetry)
{
    validateHydroBoundaryPatches(patches, mesh.minCoord, mesh.maxCoord, mesh.cellSize, radialAxisSymmetry);
    const auto boundary = hydroZeroShearBoundaryTrace(mesh, ux, uy, radialAxisSymmetry);
    const auto sampler = std::make_shared<StaggeredVelocitySampler>(mesh, ux, uy, boundary, radialAxisSymmetry);
    const auto lower = mesh.minCoord, upper = mesh.maxCoord;
    // Capture real patch data, not caller storage. Missing normal patches are
    // the existing exact zero-normal walls. An outlet returns no constraint.
    const auto wallValue = [patches](HydroBoundarySide side, double coordinate) -> std::optional<double> {
        std::optional<double> normal;
        bool outlet = false;
        for (const auto& patch : patches)
        {
            if (patch.side != side || coordinate < patch.start || coordinate > patch.end) continue;
            if (patch.type == HydroBoundaryType::PressureOutlet) outlet = true;
            else
            {
                require(!normal || *normal == patch.value,
                    "conflicting one-sided normal-velocity traces at a patch endpoint");
                normal = patch.value;
            }
        }
        if (normal)
        {
            const double sign = side == HydroBoundarySide::West || side == HydroBoundarySide::South ? -1.0 : 1.0;
            return sign * *normal;
        }
        if (outlet) return std::nullopt;
        return 0.0;
    };
    const auto matches = [&](HydroBoundarySide side, double coordinate, double actual) {
        const auto expected = wallValue(side, coordinate);
        if (!expected) return;
        const long double scale = std::max(std::abs(static_cast<long double>(*expected)),
            std::abs(static_cast<long double>(actual)));
        require(std::isfinite(actual) && std::abs(static_cast<long double>(actual)-*expected)
            <= 64.0L*std::numeric_limits<double>::epsilon()*scale,
            "native normal velocity disagrees with its specified boundary constraint");
    };
    // Use actual native layer coordinates rather than assuming that a rounded
    // lower+h equals the stored face position. No artificial exterior samples.
    double westInner = upper.x, eastInner = lower.x, southInner = upper.y, northInner = lower.y;
    for (int k = 0; k < mesh.vEdgeIndexNum(); ++k)
    {
        const auto point = mesh.vEdges(mesh.vEdgeIndexToLabel[k]).midpointCoord();
        if (point.x == lower.x) matches(HydroBoundarySide::West, point.y, ux.value[k]);
        else if (point.x == upper.x) matches(HydroBoundarySide::East, point.y, ux.value[k]);
        else
        {
            westInner = std::min(westInner, point.x);
            eastInner = std::max(eastInner, point.x);
        }
    }
    for (int k = 0; k < mesh.hEdgeIndexNum(); ++k)
    {
        const auto point = mesh.hEdges(mesh.hEdgeIndexToLabel[k]).midpointCoord();
        if (point.y == lower.y) matches(HydroBoundarySide::South, point.x, uy.value[k]);
        else if (point.y == upper.y) matches(HydroBoundarySide::North, point.x, uy.value[k]);
        else
        {
            southInner = std::min(southInner, point.y);
            northInner = std::max(northInner, point.y);
        }
    }
    require(westInner > lower.x && westInner < upper.x && eastInner >= westInner && eastInner < upper.x
        && southInner > lower.y && southInner < upper.y && northInner >= southInner && northInner < upper.y,
        "normal boundary interpolation requires an actual interior native layer on every side");
    return [sampler, wallValue, lower, upper, westInner, eastInner, southInner, northInner,
        radialAxisSymmetry](Vector2D point) {
        auto result = sampler->sample(point); // Also validates finite, closed-domain queries.
        const auto blend = [](double wall, double interior, double weight) {
            require(std::isfinite(weight) && weight >= 0.0 && weight <= 1.0,
                "invalid normal boundary-layer interpolation coordinate");
            if (weight == 0.0) return wall;
            if (weight == 1.0 || wall == interior) return interior;
            const long double value = (1.0L-weight)*wall + static_cast<long double>(weight)*interior;
            const double converted = static_cast<double>(value);
            require(std::isfinite(value) && std::isfinite(converted) && (value == 0.0L || converted != 0.0),
                "normal boundary-layer velocity is unrepresentable");
            return converted;
        };
        const auto xBoundary = [&](HydroBoundarySide side, double wall, double interior) {
            const auto value = wallValue(side, point.y);
            if (value) result.x = point.x == wall ? *value
                : blend(*value, sampler->sample({interior, point.y}).x, (point.x-wall)/(interior-wall));
        };
        const auto yBoundary = [&](HydroBoundarySide side, double wall, double interior) {
            const auto value = wallValue(side, point.x);
            if (value) result.y = point.y == wall ? *value
                : blend(*value, sampler->sample({point.x, interior}).y, (point.y-wall)/(interior-wall));
        };
        if (point.x < westInner) xBoundary(HydroBoundarySide::West, lower.x, westInner);
        else if (point.x > eastInner) xBoundary(HydroBoundarySide::East, upper.x, eastInner);
        if (point.y < southInner) yBoundary(HydroBoundarySide::South, lower.y, southInner);
        else if (point.y > northInner) yBoundary(HydroBoundarySide::North, upper.y, northInner);
        // Two physical corner constraints act on independent components. The
        // axis is stricter: each native input was already checked exactly zero.
        if (radialAxisSymmetry && point.x == 0.0) result.x = 0.0;
        return result;
    };
}
