#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "Core/Base/Vector2D.h"
#include "Core/Geometry/LevelSetPhaseFraction.h"
#include "Modules/Hydrodynamics/HydroInterfaceLoadGeometry.h"
#include "Core/Field/CellField.h"
#include "Core/Field/VEdgeField.h"
#include "Core/Field/HEdgeField.h"

// Fixed uniform-grid balanced capillary load, independent of a named case.
// The contour overload uses the original SW--NE PL geometric Hcap. The
// shared-geometry overload consumes its selected PL/regularized indicator.
// Neither is transported VOF or the diffuse material/permittivity tanh field.
// The caller supplies an ordered interface contour and its signed curvature;
// no circle fitting, global curvature averaging or empirical rescaling occurs.
namespace HydroBalancedForce
{
struct Node
{
    Vector2D coord;
    double curvature = 0.0;
};

struct Report
{
    // Sums over interior staggered dual volumes, independently for x and y.
    // RZ multiplies each dual area by 2*pi*face_radius. These are the ACTIVE
    // face loads, distinct from legacy node/Peskin quadrature diagnostics.
    // With declared wall contacts these are still INTERIOR free-row loads,
    // not a mechanical wall reaction or a contact-line support force.
    double forceX = 0.0;
    double forceY = 0.0;
    // Components live at different locations: this is max(|fx|,|fy|), not a
    // cell-centered vector magnitude reconstructed from the staggered fields.
    double maxAbsComponent = 0.0;
};

inline void require(bool condition, const char* reason)
{
    if (!condition) throw std::runtime_error(std::string("[HydroBalancedForce] ") + reason);
}

namespace Detail
{
using LevelSetPhaseFraction::Detail::Vertex;
using LevelSetPhaseFraction::Detail::negativeArea;

inline void validateNodes(const std::vector<Node>& nodes)
{
    require(nodes.size() >= 3, "at least three ordered contour nodes are required");
    for (std::size_t i = 0; i < nodes.size(); ++i)
    {
        const auto& node = nodes[i];
        const auto& next = nodes[(i + 1) % nodes.size()];
        require(std::isfinite(node.coord.x) && std::isfinite(node.coord.y) && std::isfinite(node.curvature),
            "non-finite contour coordinate or curvature");
        const double dx = next.coord.x - node.coord.x, dy = next.coord.y - node.coord.y;
        require(std::isfinite(dx * dx + dy * dy) && dx * dx + dy * dy > 0.0,
            "coincident or unrepresentable contour segment");
    }
}

inline double nearestCurvatureUnchecked(const Vector2D& point, const std::vector<Node>& nodes)
{
    double bestDistance = std::numeric_limits<double>::infinity(), result = 0.0;
    for (std::size_t i = 0; i < nodes.size(); ++i)
    {
        const auto& a = nodes[i];
        const auto& b = nodes[(i + 1) % nodes.size()];
        const double dx = b.coord.x - a.coord.x, dy = b.coord.y - a.coord.y;
        const double lengthSquared = dx * dx + dy * dy;
        const double t = std::clamp(((point.x - a.coord.x) * dx + (point.y - a.coord.y) * dy) / lengthSquared, 0.0, 1.0);
        const double ex = point.x - (a.coord.x + t * dx), ey = point.y - (a.coord.y + t * dy);
        const double distance = ex * ex + ey * ey;
        require(std::isfinite(t) && std::isfinite(distance), "non-finite nearest-segment interpolation");
        if (distance < bestDistance)
        {
            bestDistance = distance;
            result = a.curvature == b.curvature ? a.curvature
                : (1.0 - t) * a.curvature + t * b.curvature;
        }
    }
    require(std::isfinite(bestDistance) && std::isfinite(result), "no finite nearest-segment curvature");
    return result;
}
}

inline double positiveFraction(const std::array<double, 4>& phi)
{
    return LevelSetPhaseFraction::positiveFraction(phi);
}

inline double nearestSegmentCurvature(const Vector2D& point, const std::vector<Node>& nodes)
{
    require(std::isfinite(point.x) && std::isfinite(point.y), "non-finite sampling coordinate");
    Detail::validateNodes(nodes);
    return Detail::nearestCurvatureUnchecked(point, nodes);
}

// Explicit Hermite mode only. The public geometry has kappa=+div(n), whereas
// the legacy contour above supplies -div(n); capillarity is -sigma*kappa*G(H).
// Both standalone Hydro and a coupled provider consume this same snapshot.
// RZ uses surfaceCurvature=kappa_meridional+n_r/r while the frame retains
// meridional curvature for other consumers. Force density is not multiplied
// by r; only the diagnostic integral uses the shared physical dualMeasure.
// Common last-stage force arithmetic for independently validated geometry
// records. Fitted domains supply their actual C(H)/dualVolume gradient and
// physical dual measure; no second capillary formula lives in that adapter.
template<typename Mesh>
Report buildCapillaryFaceForce(const Mesh& mesh, const std::vector<HydroInterfaceLoadGeometry::Face>& faces, double sigma,
    VEdgeField<double>& fx, HEdgeField<double>& fy)
{
    require(std::isfinite(sigma) && sigma >= 0.0, "surface tension must be finite and nonnegative");
    std::vector<double> faceX(mesh.vEdgeIndexNum(), 0.0), faceY(mesh.hEdgeIndexNum(), 0.0);
    std::vector<bool> usedX(faceX.size(), false), usedY(faceY.size(), false);
    Report report;
    long double integralX = 0.0L, integralY = 0.0L;
    for (const auto& face : faces)
    {
        auto& used = face.horizontal ? usedY : usedX;
        require(face.edgeIndex >= 0 && face.edgeIndex < static_cast<int>(used.size()) && !used[face.edgeIndex],
            "invalid or duplicate capillary native face");
        used[face.edgeIndex] = true;
        require(std::isfinite(face.dualMeasure) && face.dualMeasure > 0.0, "invalid capillary physical dual measure");
        const double coefficient = -sigma * face.surfaceCurvature;
        const double force = coefficient * face.gradient;
        require(std::isfinite(coefficient) && std::isfinite(force)
            && (sigma == 0.0 || face.surfaceCurvature == 0.0 || coefficient != 0.0)
            && (coefficient == 0.0 || face.gradient == 0.0 || force != 0.0), "unrepresentable Hermite capillary face force");
        (face.horizontal ? faceY : faceX)[face.edgeIndex] = force;
        const long double contribution = static_cast<long double>(force) * face.dualMeasure;
        require(std::isfinite(contribution) && (force == 0.0 || contribution != 0.0),
            "unrepresentable Hermite capillary dual-measure load");
        (face.horizontal ? integralY : integralX) += contribution;
        report.maxAbsComponent = std::max(report.maxAbsComponent, std::abs(force));
    }
    report.forceX = static_cast<double>(integralX);
    report.forceY = static_cast<double>(integralY);
    require(std::isfinite(report.forceX) && std::isfinite(report.forceY)
        && (integralX == 0.0L || report.forceX != 0.0)
        && (integralY == 0.0L || report.forceY != 0.0), "integrated Hermite capillary force overflow or underflow");
    fx.initialize(mesh, 0.0);
    fy.initialize(mesh, 0.0);
    for (int i = 0; i < mesh.vEdgeIndexNum(); ++i) fx.value[i] = faceX[i];
    for (int i = 0; i < mesh.hEdgeIndexNum(); ++i) fy.value[i] = faceY[i];
    return report;
}

template<typename Mesh, typename LevelSet>
Report buildCapillaryForce(const Mesh& mesh, const LevelSet& levelSet,
    const HydroInterfaceLoadGeometry& geometry, double sigma,
    VEdgeField<double>& fx, HEdgeField<double>& fy)
{
    require(std::isfinite(sigma) && sigma >= 0.0, "surface tension must be finite and nonnegative");
    geometry.validateSource(mesh, levelSet);
    return buildCapillaryFaceForce(mesh, geometry.faces(), sigma, fx, fy);
}

template<typename Mesh, typename LevelSet>
Report buildCapillaryForce(const Mesh& mesh, const LevelSet& levelSet,
    const std::vector<Node>& nodes, double sigma, bool axisymmetric,
    VEdgeField<double>& fx, HEdgeField<double>& fy)
{
    require(std::isfinite(sigma) && sigma >= 0.0, "surface tension must be finite and nonnegative");
    require(mesh.Nx >= 2 && mesh.Ny >= 2 && mesh.cellIndexNum() == mesh.Nx * mesh.Ny,
        "only a complete fixed uniform mesh is supported");
    const double hx = mesh.cellSize.x, hy = mesh.cellSize.y;
    require(std::isfinite(hx) && std::isfinite(hy) && hx > 0.0 && hy > 0.0,
        "invalid uniform mesh spacing");
    require(levelSet.value && levelSet.size() == mesh.vertIndexNum(), "level set must match the mesh vertices");
    require(!axisymmetric || (std::isfinite(mesh.minCoord.x) && mesh.minCoord.x >= 0.0), "RZ requires nonnegative radius");
    Detail::validateNodes(nodes);
    auto validCell = [&](int label)
    {
        return label >= 0 && label < mesh.cellLabelCapacity() && mesh.cells(label).index() >= 0
            && mesh.cells(label).index() < mesh.cellIndexNum()
            && mesh.cellIndexToLabel[mesh.cells(label).index()] == label;
    };
    std::vector<double> fraction(mesh.cellIndexNum());
    const double tolerance = 64.0 * std::numeric_limits<double>::epsilon();
    for (int i = 0; i < mesh.cellIndexNum(); ++i)
    {
        const int label = mesh.cellIndexToLabel[i];
        require(validCell(label), "invalid cell index/label mapping");
        const auto& cell = mesh.cells(label);
        require(std::isfinite(cell.size().x) && std::isfinite(cell.size().y)
            && std::abs(cell.size().x - hx) <= tolerance * hx && std::abs(cell.size().y - hy) <= tolerance * hy,
            "coarse-fine/AMR cells are unsupported");
        for (int side = 0; side < 2; ++side)
            require(cell.cellToVEdge[side][1] == -1 && cell.cellToHEdge[side][1] == -1,
                "split faces are unsupported");
        std::array<double, 4> phi;
        for (int j = 0; j < 4; ++j)
        {
            const int vertex = cell.cellToVert[j];
            require(vertex >= 0 && vertex < mesh.vertLabelCapacity(), "invalid cell vertex label");
            const int index = mesh.verts(vertex).index();
            require(index >= 0 && index < mesh.vertIndexNum() && mesh.vertIndexToLabel[index] == vertex,
                "invalid cell vertex index");
            phi[j] = levelSet.value[index];
        }
        fraction[i] = positiveFraction(phi);
    }
    std::vector<double> faceX(mesh.vEdgeIndexNum(), 0.0), faceY(mesh.hEdgeIndexNum(), 0.0);
    Report report;
    long double integralX = 0.0L, integralY = 0.0L;
    auto valueAtFace = [&](const Edge& edge, bool horizontal, long double& integral)
    {
        if (edge.masterCell == -1 || edge.slaveCell == -1) return 0.0;
        require(validCell(edge.masterCell) && validCell(edge.slaveCell), "invalid interior face cell mapping");
        const auto& master = mesh.cells(edge.masterCell);
        const auto& slave = mesh.cells(edge.slaveCell);
        require(master.level() == slave.level(), "unequal-level face is unsupported");
        // This is precisely the fixed uniform pressureFaceData gradient and
        // dual measure. Do not take a cell gradient then average it to a face.
        const double gradient = (fraction[master.index()] - fraction[slave.index()]) / (horizontal ? hy : hx);
        if (gradient == 0.0 || sigma == 0.0) return 0.0;
        const double curvature = Detail::nearestCurvatureUnchecked(edge.midpointCoord(), nodes);
        // Keep kappa_f * G(H), NOT G(kappa*H), when curvature is nonconstant.
        const double force = sigma * curvature * gradient;
        const double measure = master.area() * (axisymmetric ? 2.0 * std::acos(-1.0) * edge.midpointCoord().x : 1.0);
        require(std::isfinite(force) && std::isfinite(measure) && measure > 0.0,
            "invalid capillary face force or dual volume");
        integral += static_cast<long double>(force) * measure;
        report.maxAbsComponent = std::max(report.maxAbsComponent, std::abs(force));
        return force;
    };
    for (int i = 0; i < mesh.vEdgeIndexNum(); ++i)
    {
        const int label = mesh.vEdgeIndexToLabel[i];
        require(label >= 0 && label < mesh.vEdgeLabelCapacity() && mesh.vEdges(label).index() == i,
            "invalid vertical face mapping");
        faceX[i] = valueAtFace(mesh.vEdges(label), false, integralX);
    }
    for (int i = 0; i < mesh.hEdgeIndexNum(); ++i)
    {
        const int label = mesh.hEdgeIndexToLabel[i];
        require(label >= 0 && label < mesh.hEdgeLabelCapacity() && mesh.hEdges(label).index() == i,
            "invalid horizontal face mapping");
        faceY[i] = valueAtFace(mesh.hEdges(label), true, integralY);
    }
    report.forceX = static_cast<double>(integralX);
    report.forceY = static_cast<double>(integralY);
    require(std::isfinite(report.forceX) && std::isfinite(report.forceY), "integrated capillary face force overflow");
    // Publish only after all geometry/force checks have succeeded.
    fx.initialize(mesh, 0.0);
    fy.initialize(mesh, 0.0);
    for (int i = 0; i < mesh.vEdgeIndexNum(); ++i) fx.value[i] = faceX[i];
    for (int i = 0; i < mesh.hEdgeIndexNum(); ++i) fy.value[i] = faceY[i];
    return report;
}
}
