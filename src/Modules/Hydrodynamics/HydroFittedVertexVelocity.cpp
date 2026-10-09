#include "Modules/Hydrodynamics/Detail/HydroFittedVertexVelocity.h"
#include "Core/Base/Interpolation.h"
#include "Core/Field/HEdgeField.h"
#include "Core/Field/VEdgeField.h"
#include "Core/Mesh/CartesianLeafMesh.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Modules/Hydrodynamics/HydroBoundary.h"
#include "Modules/Hydrodynamics/Detail/HydroFittedPressureTrace.h"
#include "Modules/Hydrodynamics/HydroFittedDomain.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace HydroFittedDetail
{
namespace
{
using Trace = HydroFittedVelocitySampler::VertexTrace;
using Type = HydroFittedDomain::BoundaryType;
using Face = HydroFittedDomain::NativeFace;
constexpr long double twoPi = 6.283185307179586476925286766559L;
void require(bool condition, const char *reason)
{
    if (!condition)
        throw std::invalid_argument(std::string("[HydroFittedVertexVelocity] ") + reason);
}
double checked(long double value)
{
    const double result = static_cast<double>(value);
    require(std::isfinite(value) && std::isfinite(result) && (value == 0.0L || result != 0.0),
            "nonfinite or unrepresentable shared vertex value");
    return result;
}
bool pressure(Type type)
{
    return type == Type::PressureOutlet || type == Type::PressureReservoir;
}

class Reconstruction
{
    const QuadtreeMesh &mesh;
    const HydroFittedDomain &domain;
    std::vector<double> speed;
    std::vector<unsigned char> used, wall;
    std::vector<Trace> traces;
    int nv;

    const Face &face(int row) const
    {
        return domain.faces()[row];
    }
    const Edge &edge(int row) const
    {
        const auto &f = face(row);
        return f.horizontal ? mesh.hEdges(f.edgeLabel) : mesh.vEdges(f.edgeLabel);
    }
    bool active(int row) const
    {
        return row >= 0 && (face(row).negativeCell >= 0 || face(row).positiveCell >= 0);
    }
    int row(bool horizontal, int label) const
    {
        if (label < 0)
            return -1;
        return horizontal ? nv + mesh.hEdges(label).index() : mesh.vEdges(label).index();
    }
    double value(int native) const
    {
        require(active(native) && std::isfinite(speed[native]),
                "velocity stencil consumes a solid-only face");
        return speed[native];
    }
    const Cell &cell(int index) const
    {
        require(index >= 0 && index < mesh.cellIndexNum() && domain.cellRowByIndex()[index] >= 0,
                "velocity stencil leaves the active fluid domain");
        return mesh.cells(mesh.cellIndexToLabel[index]);
    }
    int singleFace(int owner, bool horizontal, int side) const
    {
        const auto &c = cell(owner);
        const int *labels = horizontal ? c.cellToHEdge[side] : c.cellToVEdge[side];
        require(labels[0] >= 0 && labels[1] < 0,
                "physical boundary trace requires a locally unsplit native stencil");
        const int native = row(horizontal, labels[0]);
        require(active(native), "physical boundary trace has no active native face");
        return native;
    }
    double pressureTrace(int native, int tangentSide) const
    {
        // Same fitted incomingTangent/pressure-corner closure. The derivative
        // comes from actual transverse pressure faces, never a solid ghost.
        const auto &f = face(native);
        const auto transverse = [&](int owner) { return singleFace(owner, !f.horizontal, tangentSide); };
        double derivative = 0.0;
        if (f.negativeCell >= 0 && f.positiveCell >= 0)
        {
            const auto trace = pressureTraceGeometry(mesh, domain, native, tangentSide);
            derivative = (value(trace.highNative) - value(trace.lowNative)) / trace.derivativeSpacing;
            return checked(HydroBoundaryDetail::zeroShearTrace(value(native), trace.traceSpacing,
                tangentSide ? 1.0 : -1.0, derivative));
        }
        else
        {
            require(pressure(f.boundary.type), "pressure corner lacks its second pressure side");
            const int owner = f.negativeCell >= 0 ? f.negativeCell : f.positiveCell;
            const int inner = singleFace(owner, f.horizontal, f.outwardSign > 0.0 ? 0 : 1);
            const int innerOwner = f.outwardSign > 0.0 ? face(inner).negativeCell : face(inner).positiveCell;
            const int here = transverse(owner), inside = transverse(innerOwner);
            const auto &a = face(here);
            const auto &b = face(inside);
            require(pressure(a.boundary.type) && pressure(b.boundary.type) &&
                        a.outwardSign == b.outwardSign &&
                        (a.horizontal ? a.center.y == b.center.y : a.center.x == b.center.x),
                    "pressure-corner trace leaves its real boundary plane");
            derivative = f.outwardSign * (value(here) - value(inside)) / f.normalSpacing;
        }
        return checked(HydroBoundaryDetail::zeroShearTrace(value(native), f.tangentSpacing,
                                                           tangentSide ? 1.0 : -1.0, derivative));
    }
    std::vector<int> incidentCells(const Vert &vertex) const
    {
        std::vector<int> result;
        for (int label : vertex.vertToCell)
            if (label >= 0 && domain.cellRowByIndex()[mesh.cells(label).index()] >= 0)
                result.push_back(mesh.cells(label).index());
        std::sort(result.begin(), result.end());
        result.erase(std::unique(result.begin(), result.end()), result.end());
        return result;
    }
    double hangingEstimate(const Vert &vertex, bool horizontal) const
    {
        // Reuse the old first-order native edge-to-vertex least-squares fit,
        // but collect every REAL subface of the active incident leaves. The
        // old fixed slot lists are not valid for a fitted solid mask.
        std::vector<int> rows;
        Vector2D scale{0, 0};
        for (int owner : incidentCells(vertex))
        {
            const auto &c = cell(owner);
            scale.x = std::max(scale.x, c.size().x);
            scale.y = std::max(scale.y, c.size().y);
            for (int side = 0; side < 2; ++side)
                for (int part = 0; part < 2; ++part)
                {
                    const int native =
                        row(horizontal, horizontal ? c.cellToHEdge[side][part] : c.cellToVEdge[side][part]);
                    if (active(native))
                        rows.push_back(native);
                }
        }
        std::sort(rows.begin(), rows.end());
        rows.erase(std::unique(rows.begin(), rows.end()), rows.end());
        require(rows.size() >= 3 && scale.x > 0.0 && scale.y > 0.0,
                "hanging velocity fit lacks native support");
        std::vector<Vector2D> points;
        const auto origin = vertex.coord();
        for (int native : rows)
            points.push_back(
                {(face(native).center.x - origin.x) / scale.x, (face(native).center.y - origin.y) / scale.y});
        // Precheck rank before the old normal-equation inverse; never invoke
        // its legacy exit(0) for a short stencil or replace failure with zero.
        long double maxArea = 0.0L;
        for (std::size_t i = 1; i < points.size(); ++i)
            for (std::size_t j = i + 1; j < points.size(); ++j)
                maxArea = std::max(maxArea, std::abs((static_cast<long double>(points[i].x) - points[0].x) *
                                                         (points[j].y - points[0].y) -
                                                     (static_cast<long double>(points[j].x) - points[0].x) *
                                                         (points[i].y - points[0].y)));
        require(maxArea > 0.0L, "singular hanging velocity stencil");
        const auto weights = Interpolation::leastSquareInterpolationCoeff(points, {0, 0});
        long double estimate = 0.0L, one = 0.0L, x = 0.0L, y = 0.0L;
        for (std::size_t k = 0; k < weights.size(); ++k)
        {
            require(std::isfinite(weights[k]), "nonfinite hanging interpolation coefficient");
            estimate += static_cast<long double>(weights[k]) * value(rows[k]);
            one += weights[k];
            x += weights[k] * points[k].x;
            y += weights[k] * points[k].y;
        }
        require(std::abs(one - 1) <= 1e-10L && std::abs(x) <= 1e-10L && std::abs(y) <= 1e-10L,
                "hanging interpolation lost affine consistency");
        return checked(estimate);
    }
    double component(const Vert &vertex, bool horizontal) const
    {
        const int *labels = horizontal ? vertex.vertToHEdge : vertex.vertToVEdge;
        const int lo = row(horizontal, labels[0]), hi = row(horizontal, labels[1]);
        const auto p = vertex.coord();
        if (active(lo) && active(hi))
        {
            const long double a = horizontal ? face(lo).center.x : face(lo).center.y;
            const long double b = horizontal ? face(hi).center.x : face(hi).center.y;
            const long double t = ((horizontal ? p.x : p.y) - a) / (b - a);
            require(t >= 0.0L && t <= 1.0L, "native vertex is not bracketed by incident face centers");
            if (value(lo) == value(hi))
                return value(lo);
            return checked((1.0L - t) * value(lo) + t * value(hi));
        }
        const int native = active(lo) ? lo : hi;
        require(active(native), "active vertex lacks an incident component face");
        if (face(native).boundary.type == Type::Symmetry)
            return value(native);
        const int missingSide = active(lo) ? 1 : 0;
        const int *transverseLabels = horizontal ? vertex.vertToVEdge : vertex.vertToHEdge;
        int boundary = -1;
        for (int k = 0; k < 2; ++k)
        {
            const int transverse = row(!horizontal, transverseLabels[k]);
            if (!active(transverse) || face(transverse).boundary.type == Type::Interior)
                continue;
            if (boundary >= 0)
                require(face(boundary).boundary.type == face(transverse).boundary.type &&
                            face(boundary).outwardSign == face(transverse).outwardSign,
                        "vertex changes physical boundary type along its trace");
            boundary = transverse;
        }
        if (boundary < 0)
            return hangingEstimate(vertex, horizontal);
        const auto type = face(boundary).boundary.type;
        if (type == Type::NoSlipWall)
            return 0.0;
        if (type == Type::Symmetry)
        {
            if (horizontal && p.x == 0.0)
            {
                const auto &end = mesh.verts(edge(native).edgeToVert[1]);
                const int second = row(true, end.vertToHEdge[1]);
                require(active(second) && face(second).center.y == p.y,
                        "axis even trace lacks two actual inward horizontal samples");
                const long double ra = face(native).center.x, rb = face(second).center.x;
                require(rb > ra && ra > 0.0L, "axis trace has invalid radii");
                return checked(static_cast<long double>(value(native)) -
                               (static_cast<long double>(value(second)) - value(native)) * ra * ra /
                                   ((rb - ra) * (rb + ra)));
            }
            return value(native);
        }
        require(pressure(type), "unsupported physical tangential velocity trace");
        return pressureTrace(native, missingSide);
    }

  public:
    Reconstruction(const QuadtreeMesh &m, const HydroFittedDomain &d, const VEdgeField<double> &ux,
                   const HEdgeField<double> &uy)
        : mesh(m), domain(d), nv(m.vEdgeIndexNum())
    {
        domain.validateMeshIdentity(mesh);
        validateCartesianLeafMesh(mesh);
        require(domain.axisymmetric(), "native fitted vertex reconstruction requires RZ");
        require(ux.size() == nv && uy.size() == mesh.hEdgeIndexNum() && ux.value && uy.value,
                "native field layout mismatch");
        const double nan = std::numeric_limits<double>::quiet_NaN();
        speed.assign(domain.faces().size(), nan);
        used.assign(mesh.vertIndexNum(), 0);
        wall = used;
        traces.assign(mesh.vertIndexNum(), {{nan, nan}, nan});
        for (const auto &f : domain.faces())
        {
            if (!active(f.nativeRow))
                continue;
            const int index = f.horizontal ? f.nativeRow - nv : f.nativeRow;
            require((f.horizontal ? uy.label(index) : ux.label(index)) == f.edgeLabel,
                    "stale active native field label");
            const double v = f.horizontal ? uy.value[index] : ux.value[index];
            require(std::isfinite(v) && (f.velocityRow >= 0 || v == 0.0),
                    "invalid active or constrained native speed");
            speed[f.nativeRow] = v;
            for (int label : edge(f.nativeRow).edgeToVert)
            {
                const int k = mesh.verts(label).index();
                used[k] = 1;
                if (f.boundary.type == Type::NoSlipWall)
                    wall[k] = 1;
            }
        }
    }
    std::vector<Trace> build()
    {
        for (int k = 0; k < mesh.vertIndexNum(); ++k)
            if (used[k])
            {
                const auto &v = mesh.verts(mesh.vertIndexToLabel[k]);
                traces[k].velocity =
                    wall[k] ? Vector2D{0, 0} : Vector2D{component(v, false), component(v, true)};
            }
        std::vector<long double> net(mesh.cellIndexNum(), 0.0L);
        for (const auto &f : domain.faces())
            if (active(f.nativeRow))
            {
                const long double q = static_cast<long double>(f.area) * value(f.nativeRow);
                if (f.negativeCell >= 0)
                    net[f.negativeCell] += q;
                if (f.positiveCell >= 0)
                    net[f.positiveCell] -= q;
            }
        for (int k = 0; k < mesh.vertIndexNum(); ++k)
            if (used[k])
            {
                if (wall[k])
                {
                    traces[k].axialDerivative = 0.0;
                    continue;
                }
                const auto &vertex = mesh.verts(mesh.vertIndexToLabel[k]);
                const auto p = vertex.coord();
                const auto owners = incidentCells(vertex);
                long double sum = 0.0L;
                require(!owners.empty(), "shared derivative has no active incident cell");
                for (int owner : owners)
                {
                    const auto &c = cell(owner);
                    std::array<Vector2D, 4> u;
                    for (int j = 0; j < 4; ++j)
                        u[j] = traces[mesh.verts(c.cellToVert[j]).index()].velocity;
                    const auto low = mesh.verts(c.cellToVert[0]).coord();
                    const auto h = c.size();
                    const long double x = (static_cast<long double>(p.x) - low.x) / h.x,
                                      y = (static_cast<long double>(p.y) - low.y) / h.y;
                    const long double dz = ((1.0L - x) * (u[3].y - u[0].y) + x * (u[2].y - u[1].y)) / h.y;
                    if (p.x == 0.0)
                    {
                        sum += dz;
                        continue;
                    }
                    const long double dr = ((1.0L - y) * (u[1].x - u[0].x) + y * (u[2].x - u[3].x)) / h.x;
                    // Physical form of the original bilinear mixed streamfunction
                    // estimate. Average independent axial and radial estimates,
                    // INCLUDING native FV divergence, then share once per vertex.
                    sum += 0.5L *
                           (dz - dr - traces[k].velocity.x / p.x + net[owner] / (twoPi * p.x * h.x * h.y));
                }
                traces[k].axialDerivative = checked(sum / owners.size());
            }
        return std::move(traces);
    }
};
} // namespace
std::vector<Trace> reconstructVertexVelocity(const QuadtreeMesh &mesh, const HydroFittedDomain &domain,
                                             const VEdgeField<double> &ux, const HEdgeField<double> &uy)
{
    return Reconstruction(mesh, domain, ux, uy).build();
}
} // namespace HydroFittedDetail
