#include "Modules/Hydrodynamics/HydroFittedVelocity.h"
#include "Modules/Hydrodynamics/HydroFittedDomain.h"
#include "Modules/Hydrodynamics/HydroBoundary.h"
#include "Modules/Hydrodynamics/Detail/HydroFittedAdaptiveVelocity.h"
#include "Modules/Hydrodynamics/Detail/HydroFittedVertexVelocity.h"
#include "Modules/Hydrodynamics/Detail/HydroVelocityHermiteBasis.h"
#include "Core/Field/VEdgeField.h"
#include "Core/Field/HEdgeField.h"
#include "Core/Geometry/LevelSetMappedSegment.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/UniformLeafGrid.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
using Type = HydroFittedDomain::BoundaryType;
using Face = HydroFittedDomain::NativeFace;
void require(bool condition, const char* reason)
{
    if (!condition) throw std::invalid_argument(std::string("[HydroFittedVelocitySampler] ") + reason);
}
bool pressure(Type type) { return type == Type::PressureReservoir || type == Type::PressureOutlet; }
double checked(long double value)
{
    const double result = static_cast<double>(value);
    require(std::isfinite(value) && std::isfinite(result) && (value == 0.0L || result != 0.0),
        "unrepresentable interpolated velocity");
    return result;
}
// Same expression and endpoint/constant shortcuts as the Core native sampler.
double lerp(double a, double b, double weight)
{
    require(std::isfinite(weight) && weight >= 0.0 && weight <= 1.0, "invalid interpolation coordinate");
    if (weight == 0.0) return a;
    if (weight == 1.0) return b;
    if (a == b) return a;
    return checked((1.0L - weight) * a + static_cast<long double>(weight) * b);
}


using Quad = std::array<Vector2D, 4>; // Counterclockwise: SW, SE, NE, NW.
constexpr double twoPi = 6.283185307179586476925286766559;
struct Mapping
{
    Vector2D point, xi, eta, mixed;
    double determinant = 0.0;
};
Mapping mapped(const Quad& q, double u, double v)
{
    Mapping map;
    map.point = {lerp(lerp(q[0].x,q[1].x,u),lerp(q[3].x,q[2].x,u),v),
                 lerp(lerp(q[0].y,q[1].y,u),lerp(q[3].y,q[2].y,u),v)};
    // Each retained subquad spans one half of each whole-cell coordinate.
    map.xi = {2.0*lerp(q[1].x-q[0].x,q[2].x-q[3].x,v),
              2.0*lerp(q[1].y-q[0].y,q[2].y-q[3].y,v)};
    map.eta = {2.0*lerp(q[3].x-q[0].x,q[2].x-q[1].x,u),
               2.0*lerp(q[3].y-q[0].y,q[2].y-q[1].y,u)};
    map.mixed = {4.0*((q[2].x-q[3].x)-(q[1].x-q[0].x)),
                 4.0*((q[2].y-q[3].y)-(q[1].y-q[0].y))};
    map.determinant = map.xi.x*map.eta.y-map.xi.y*map.eta.x;
    require(std::isfinite(map.determinant) && map.determinant>0.0,"invalid physical flux map Jacobian");
    return map;
}
using HydroFittedDetail::Basis;
long double cross(Vector2D a, Vector2D b, Vector2D p)
{
    return (static_cast<long double>(b.x) - a.x) * (static_cast<long double>(p.y) - a.y)
        - (static_cast<long double>(b.y) - a.y) * (static_cast<long double>(p.x) - a.x);
}

// Actual shared vertices/midpoints define containment, not a regenerated axis
// or an epsilon-expanded box. Only AFTER this test may inverse-map roundoff be
// bounded to [0,1]; the physical query is never moved or accepted by tolerance.
bool coordinates(const Quad& q, Vector2D point, std::array<double, 2>& result)
{
    for (int k = 0; k < 4; ++k)
    {
        require(cross(q[k], q[(k + 1) % 4], q[(k + 2) % 4]) > 0.0L,
            "unresolved or nonconvex actual interpolation element");
        if (cross(q[k], q[(k + 1) % 4], point) < 0.0L) return false;
    }
    if (q[0].x == q[3].x && q[1].x == q[2].x
        && q[0].y == q[1].y && q[2].y == q[3].y)
    {
        result = {{(point.x - q[0].x) / (q[1].x - q[0].x),
                   (point.y - q[0].y) / (q[3].y - q[0].y)}};
        return true;
    }
    // Isoparametric bilinear coordinates of the retained four physical points.
    // Centered differences retain the small non-Cartesian rounding offsets.
    const long double ax = static_cast<long double>(q[1].x) - q[0].x;
    const long double ay = static_cast<long double>(q[1].y) - q[0].y;
    const long double bx = static_cast<long double>(q[3].x) - q[0].x;
    const long double by = static_cast<long double>(q[3].y) - q[0].y;
    const long double cx = (static_cast<long double>(q[2].x) - q[3].x) - ax;
    const long double cy = (static_cast<long double>(q[2].y) - q[3].y) - ay;
    const long double px = static_cast<long double>(point.x) - q[0].x;
    const long double py = static_cast<long double>(point.y) - q[0].y;
    const long double det = ax * by - ay * bx;
    require(det > 0.0L && std::isfinite(det), "singular actual interpolation element");
    long double u = (px * by - py * bx) / det, v = (ax * py - ay * px) / det;
    for (int iteration = 0; iteration < 8; ++iteration)
    {
        const long double ex = ax * u + bx * v + cx * u * v - px;
        const long double ey = ay * u + by * v + cy * u * v - py;
        const long double ux = ax + cx * v, uy = ay + cy * v;
        const long double vx = bx + cx * u, vy = by + cy * u;
        const long double jacobian = ux * vy - uy * vx;
        require(jacobian > 0.0L && std::isfinite(jacobian), "singular actual bilinear map");
        const long double du = (ex * vy - ey * vx) / jacobian;
        const long double dv = (ux * ey - uy * ex) / jacobian;
        u -= du; v -= dv;
        if (du == 0.0L && dv == 0.0L) break;
    }
    const long double tolerance = 64.0L * std::numeric_limits<double>::epsilon();
    require(std::isfinite(u) && std::isfinite(v) && u >= -tolerance && u <= 1.0L + tolerance
        && v >= -tolerance && v <= 1.0L + tolerance,
        "actual bilinear inverse did not resolve its contained query");
    const long double residual = std::max(std::abs(ax * u + bx * v + cx * u * v - px),
        std::abs(ay * u + by * v + cy * u * v - py));
    require(residual <= tolerance * std::max({std::abs(ax), std::abs(ay), std::abs(bx), std::abs(by)}),
        "actual bilinear inverse residual is unresolved");
    result = {{static_cast<double>(std::max(0.0L, std::min(1.0L, u))),
               static_cast<double>(std::max(0.0L, std::min(1.0L, v)))}};
    return true;
}
}

struct HydroFittedVelocitySampler::Impl
{
    int nx = 0, ny = 0;
    bool axisymmetric = false;
    Vector2D lower, spacing;
    std::vector<Vector2D> vertices, centers; // Every actual point, not first-row/column axes.
    std::vector<int> cells, vertical, horizontal;
    std::vector<std::array<int, 2>> logicalCell;
    std::vector<unsigned char> active, noSlipVertex;
    std::vector<Face> faces;
    std::vector<double> values; // Solid-only entries remain NaN and are never consumed.
    Reconstruction reconstruction;
    struct FluxCell
    {
        double q[4]{}; // W,E,S,N, all oriented in the positive coordinate direction.
        // SW,SE,NW,NE; streamfunction values are zero. On axis cells the
        // first coordinate is s=xi^2: these are Phi_s, Phi_eta, Phi_s_eta.
        // Elsewhere they remain the original psi_xi, psi_eta, psi_xi_eta.
        struct Jet { double dx=0.0,dy=0.0,dxy=0.0; } jet[4];
        unsigned walls=0;
        bool axis=false;
    };
    std::vector<FluxCell> fluxCells;
    std::vector<std::optional<HydroFittedFaceFlux>> nativeFlux;
    std::unique_ptr<const HydroFittedDetail::AdaptiveVelocity> adaptive;

    Impl(const QuadtreeMesh& mesh, const HydroFittedDomain& domain,
        const VEdgeField<double>& ux, const HEdgeField<double>& uy,
        const std::vector<VertexTrace>& traces)
        : reconstruction(Reconstruction::RZFluxCompatible),
          adaptive(std::make_unique<const HydroFittedDetail::AdaptiveVelocity>(mesh,domain,ux,uy,traces)) {}

    Impl(const QuadtreeMesh& mesh, const NativeLeafDomain& domain,
        const VEdgeField<double>& ux, const HEdgeField<double>& uy,
        const std::vector<VertexTrace>& traces)
        : reconstruction(Reconstruction::RZFluxCompatible),
          adaptive(std::make_unique<const HydroFittedDetail::AdaptiveVelocity>(mesh,domain,ux,uy,traces)) {}

    Impl(const QuadtreeMesh& mesh, const HydroFittedDomain& domain,
        const VEdgeField<double>& ux, const HEdgeField<double>& uy, Reconstruction method)
        : reconstruction(method)
    {
        require(method==Reconstruction::StaggeredBilinear || method==Reconstruction::RZFluxCompatible,
            "unknown velocity reconstruction");
        domain.validateMeshIdentity(mesh);
        if (!domain.uniformLeaves() && method == Reconstruction::RZFluxCompatible)
        {
            adaptive = std::make_unique<const HydroFittedDetail::AdaptiveVelocity>(mesh, domain, ux, uy,
                HydroFittedDetail::reconstructVertexVelocity(mesh, domain, ux, uy));
            return;
        }
        const UniformLeafGrid grid(mesh, CartesianGridScope::UniformLeaves);
        nx = grid.nx(); ny = grid.ny(); axisymmetric = domain.axisymmetric();
        require(nx >= 2 && ny >= 2, "at least two actual leaf cells per direction required");
        const int nv = mesh.vEdgeIndexNum(), nh = mesh.hEdgeIndexNum();
        require(ux.size() == nv && uy.size() == nh && ux.value && uy.value,
            "native velocity storage does not match the mesh");
        faces = domain.faces();
        require(faces.size() == static_cast<std::size_t>(nv) + nh, "incomplete fitted native layout");
        lower = grid.minCoord(); spacing = grid.spacing();
        vertices.resize(grid.vertices().size());
        for (int j = 0; j <= ny; ++j) for (int i = 0; i <= nx; ++i)
        {
            const std::size_t k = i + static_cast<std::size_t>(nx + 1) * j;
            const auto point = mesh.verts(grid.vertices()[k]).coord();
            require(std::isfinite(point.x) && std::isfinite(point.y)
                && (!i || point.x > vertices[k - 1].x)
                && (!j || point.y > vertices[k - (nx + 1)].y), "unresolved actual vertex coordinates");
            require(!axisymmetric || point.x >= 0.0, "RZ coordinates must be nonnegative");
            vertices[k] = point;
        }
        cells.resize(grid.cells().size()); active.resize(cells.size()); centers.resize(cells.size());
        logicalCell.resize(mesh.cellIndexNum());
        for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i)
        {
            const int slot = i + nx * j, k = mesh.cells(grid.cells()[slot]).index();
            cells[slot] = k; logicalCell[k] = {{i, j}};
            active[slot] = domain.cellRowByIndex()[k] >= 0;
            centers[slot] = mesh.cells(grid.cells()[slot]).centerCoord();
        }
        vertical.resize(grid.vertical().size()); horizontal.resize(grid.horizontal().size());
        for (std::size_t i = 0; i < vertical.size(); ++i) vertical[i] = mesh.vEdges(grid.vertical()[i]).index();
        for (std::size_t i = 0; i < horizontal.size(); ++i) horizontal[i] = nv + mesh.hEdges(grid.horizontal()[i]).index();
        values.assign(faces.size(), std::numeric_limits<double>::quiet_NaN());
        noSlipVertex.assign(static_cast<std::size_t>(nx + 1) * (ny + 1), 0);
        for (const auto& face : faces)
        {
            if (face.negativeCell < 0 && face.positiveCell < 0) continue;
            const int index = face.horizontal ? face.nativeRow - nv : face.nativeRow;
            require(face.nativeRow >= 0 && index >= 0 && index < (face.horizontal ? nh : nv)
                && (face.horizontal ? uy.label(index) : ux.label(index)) == face.edgeLabel,
                "active native field label/index mismatch");
            const double value = face.horizontal ? uy.value[index] : ux.value[index];
            require(std::isfinite(value), "nonfinite active native velocity");
            if (face.velocityRow < 0)
                require((face.boundary.type == Type::NoSlipWall || face.boundary.type == Type::Symmetry) && value == 0.0,
                    "actual wall/symmetry normal velocity must be exactly zero");
            values[face.nativeRow] = value;
            if (face.boundary.type == Type::NoSlipWall)
            {
                const int owner = face.negativeCell >= 0 ? face.negativeCell : face.positiveCell;
                const auto ij = logicalCell[owner];
                const int i = ij[0] + (!face.horizontal && face.negativeCell >= 0 ? 1 : 0);
                const int j = ij[1] + (face.horizontal && face.negativeCell >= 0 ? 1 : 0);
                noSlipVertex[i + (nx + 1) * j] = 1;
                noSlipVertex[i + (face.horizontal ? 1 : 0) + (nx + 1) * (j + (face.horizontal ? 0 : 1))] = 1;
            }
        }
        if (reconstruction==Reconstruction::RZFluxCompatible)
        {
            validateFluxMeasures(domain);
            initializeFluxReconstruction();
            initializeNativeFluxTraces();
        }
    }

    bool cellActive(int i, int j) const
    { return i >= 0 && i < nx && j >= 0 && j < ny && active[i + nx * j]; }
    int cellFace(int i, int j, bool h, int side) const
    {
        require(cellActive(i, j), "velocity stencil attempts to consume an inactive cell");
        return h ? horizontal[i + nx * (j + side)] : vertical[i + side + (nx + 1) * j];
    }
    double value(int native) const
    {
        require(native >= 0 && native < static_cast<int>(values.size()) && std::isfinite(values[native]),
            "velocity stencil attempts to consume a solid-only face");
        return values[native];
    }
    Vector2D vertex(int i, int j) const { return vertices[i + (nx + 1) * j]; }
    std::array<Quad,4> cellParts(int i, int j) const
    {
        // A shared physical side is vertex -> its actual native midpoint ->
        // vertex. The four subelements therefore meet without first-column
        // approximations, including at the non-binary-spacing lip.
        const auto sw = vertex(i, j), se = vertex(i + 1, j);
        const auto nw = vertex(i, j + 1), ne = vertex(i + 1, j + 1);
        const auto w = faces[cellFace(i, j, false, 0)].center;
        const auto e = faces[cellFace(i, j, false, 1)].center;
        const auto s = faces[cellFace(i, j, true, 0)].center;
        const auto n = faces[cellFace(i, j, true, 1)].center;
        const auto c = centers[i + nx * j];
        return {{Quad{{sw,s,c,w}}, Quad{{s,se,e,c}}, Quad{{w,c,n,nw}}, Quad{{c,e,ne,n}}}};
    }
    Mapping cellMap(int i, int j, double u, double v) const
    {
        const int x=u>=0.5?1:0,y=v>=0.5?1:0;
        return mapped(cellParts(i,j)[x+2*y],2*u-x,2*v-y);
    }
    bool cellCoordinates(int i, int j, Vector2D point, std::array<double, 2>& result) const
    {
        if (!cellActive(i, j)) return false;
        const auto parts=cellParts(i,j);
        for (int part = 0; part < 4; ++part)
        {
            std::array<double, 2> local;
            if (coordinates(parts[part], point, local))
            {
                result = {{0.5 * ((part % 2) + local[0]), 0.5 * ((part / 2) + local[1])}};
                return true;
            }
        }
        return false;
    }
    std::array<int, 2> selectCell(Vector2D point) const
    {
        // Canonical coordinates only narrow the search. All acceptance below
        // uses retained actual geometry, not rounded indices or tolerances.
        const double gx = std::floor((point.x - lower.x) / spacing.x);
        const double gy = std::floor((point.y - lower.y) / spacing.y);
        require(std::isfinite(gx) && std::isfinite(gy) && gx >= -1.0 && gx <= nx
            && gy >= -1.0 && gy <= ny, "velocity query outside actual closed domain");
        const int i = static_cast<int>(gx), j = static_cast<int>(gy);
        std::array<double, 2> local;
        for (int jj : {j, j - 1, j + 1}) for (int ii : {i, i - 1, i + 1})
            if (cellCoordinates(ii, jj, point, local)) return {{ii, jj}};
        throw std::invalid_argument("[HydroFittedVelocitySampler] query has no actual active incident cell");
    }
    bool wallVertex(int i, int j) const { return noSlipVertex[i + (nx + 1) * j] != 0; }

    double pressureTrace(int native, int ownerI, int ownerJ, int tangentSide) const
    {
        const auto& face = faces[native];
        const bool h = face.horizontal;
        const auto transverse = [&](int cell) {
            const auto ij = logicalCell[cell];
            return cellFace(ij[0], ij[1], !h, tangentSide);
        };
        double derivative = 0.0;
        if (face.negativeCell >= 0 && face.positiveCell >= 0)
        {
            const int low = transverse(face.negativeCell), high = transverse(face.positiveCell);
            require(pressure(faces[low].boundary.type) && pressure(faces[high].boundary.type),
                "zero-shear trace requires the same real pressure boundary on both sides");
            derivative = (value(high) - value(low)) / face.normalSpacing;
        }
        else
        {
            const int inward = face.outwardSign > 0.0 ? -1 : 1;
            const int ii = ownerI + (h ? 0 : inward), jj = ownerJ + (h ? inward : 0);
            require(cellActive(ii, jj), "pressure-corner trace lacks its second inward fluid cell");
            const int here = cellFace(ownerI, ownerJ, !h, tangentSide), inside = cellFace(ii, jj, !h, tangentSide);
            require(pressure(faces[here].boundary.type) && pressure(faces[inside].boundary.type)
                && faces[here].horizontal == faces[inside].horizontal
                && faces[here].outwardSign == faces[inside].outwardSign,
                "pressure-corner trace leaves its actual boundary plane");
            derivative = face.outwardSign * (value(here) - value(inside)) / face.normalSpacing;
        }
        // Same pressure-boundary tangential closure as fitted incomingTangent.
        return checked(HydroBoundaryDetail::zeroShearTrace(value(native), face.tangentSpacing,
            tangentSide ? 1.0 : -1.0, derivative));
    }

    double edgeSample(int native, int i, int j, double coordinate) const
    {
        const auto& face = faces[native];
        if (face.boundary.type == Type::NoSlipWall || face.boundary.type == Type::Symmetry) return value(native);
        const bool h = face.horizontal;
        const double center = h ? face.center.x : face.center.y;
        const double current = value(native);
        if (coordinate == center) return current;
        const int side = coordinate > center ? 1 : 0, direction = side ? 1 : -1;
        const int normalSide = native == cellFace(i, j, h, 1) ? 1 : 0;
        const int vi = i + (h ? side : normalSide), vj = j + (h ? normalSide : side);
        const auto actualEndpoint = vertex(vi, vj);
        const double endpoint = h ? actualEndpoint.x : actualEndpoint.y;
        const double endpointWeight = (coordinate - center) / (endpoint - center);
        if (wallVertex(vi, vj)) return lerp(current, 0.0, endpointWeight);
        if (axisymmetric && h && endpoint == 0.0)
        {
            require(i == 0 && side == 0 && cellActive(1, j), "axis even trace lacks two inward active cells");
            const int second = cellFace(1, j, true, normalSide);
            const long double ra = center, rb = faces[second].center.x, r = coordinate;
            require(rb > ra && ra > 0.0, "invalid native axis-even radii");
            return checked(static_cast<long double>(current) + (static_cast<long double>(value(second)) - current)
                * ((r - ra) * (r + ra)) / ((rb - ra) * (rb + ra)));
        }
        const int ni = i + (h ? direction : 0), nj = j + (h ? 0 : direction);
        if (cellActive(ni, nj))
        {
            const int neighbor = cellFace(ni, nj, h, normalSide);
            const double other = h ? faces[neighbor].center.x : faces[neighbor].center.y;
            if (side) return lerp(current, value(neighbor), (coordinate - center) / (other - center));
            return lerp(value(neighbor), current, (coordinate - other) / (center - other));
        }
        const auto& boundary = faces[cellFace(i, j, !h, side)];
        double trace = 0.0;
        if (boundary.boundary.type == Type::NoSlipWall) trace = 0.0;
        else if (boundary.boundary.type == Type::Symmetry) trace = current;
        else
        {
            require(pressure(boundary.boundary.type), "missing physical boundary tangential trace");
            trace = pressureTrace(native, i, j, side);
        }
        return lerp(current, trace, endpointWeight);
    }

    double component(Vector2D point, int i, int j, bool h) const
    {
        const int first = cellFace(i, j, h, 0), second = cellFace(i, j, h, 1);
        const double tangent = h ? point.x : point.y;
        const double center = h ? faces[first].center.x : faces[first].center.y;
        const int preferred = tangent >= center ? 1 : 0;
        for (int side : {preferred, 1 - preferred})
        {
            const int direction = side ? 1 : -1;
            const int ni = i + (h ? direction : 0), nj = j + (h ? 0 : direction);
            const bool cleanVertex = h ? (!wallVertex(i + side, j) && !wallVertex(i + side, j + 1))
                : (!wallVertex(i, j + side) && !wallVertex(i + 1, j + side));
            if (!cellActive(ni, nj) || !cleanVertex) continue;
            const int otherFirst = cellFace(ni, nj, h, 0), otherSecond = cellFace(ni, nj, h, 1);
            if (faces[first].velocityRow >= 0 && faces[second].velocityRow >= 0
                && faces[otherFirst].velocityRow >= 0 && faces[otherSecond].velocityRow >= 0)
            {
                // Complete native element: strictly rectangular inputs retain
                // the Core coordinates and x-then-y arithmetic exactly. Other
                // accepted leaf layouts use all four actual native positions.
                const int a = side ? first : otherFirst, b = side ? otherFirst : first;
                const int c = side ? second : otherSecond, d = side ? otherSecond : second;
                std::array<double, 2> uv;
                if (h)
                {
                    if (coordinates(Quad{{faces[a].center, faces[b].center, faces[d].center, faces[c].center}}, point, uv))
                        return lerp(lerp(value(a), value(b), uv[0]), lerp(value(c), value(d), uv[0]), uv[1]);
                }
                else
                {
                    if (coordinates(Quad{{faces[a].center, faces[c].center, faces[d].center, faces[b].center}}, point, uv))
                        return lerp(lerp(value(a), value(c), uv[0]), lerp(value(b), value(d), uv[0]), uv[1]);
                }
            }
        }
        std::array<double, 2> uv;
        require(cellCoordinates(i, j, point, uv), "query left its actual selected cell");
        const double along = h ? uv[0] : uv[1];
        const auto edgeCoordinate = [&](int native, int normalSide) {
            const auto midpoint = faces[native].center;
            const bool upper = along > 0.5;
            const auto endpoint = vertex(i + (h ? (upper ? 1 : 0) : normalSide),
                j + (h ? normalSide : (upper ? 1 : 0)));
            const double a = h ? midpoint.x : midpoint.y, b = h ? endpoint.x : endpoint.y;
            // This is a point ON this retained half-edge, not a ghost sample
            // or reuse of a different column's endpoint coordinate.
            return lerp(a, b, upper ? 2.0 * along - 1.0 : 1.0 - 2.0 * along);
        };
        const double lo = edgeSample(first, i, j, edgeCoordinate(first, 0));
        const double hi = edgeSample(second, i, j, edgeCoordinate(second, 1));
        return lerp(lo, hi, h ? uv[1] : uv[0]);
    }

    Vector2D legacySample(Vector2D point) const
    {
        require(std::isfinite(point.x) && std::isfinite(point.y), "nonfinite query coordinate");
        require(!axisymmetric || point.x >= 0.0, "negative-radius velocity query");
        const auto ij = selectCell(point);
        const int i = ij[0], j = ij[1];
        for (int vj : {j, j + 1}) for (int vi : {i, i + 1})
            if (point.x == vertex(vi, vj).x && point.y == vertex(vi, vj).y && wallVertex(vi, vj)) return {0.0, 0.0};
        return {component(point, i, j, false), component(point, i, j, true)};
    }

    void validateFluxMeasures(const HydroFittedDomain& domain) const
    {
        require(axisymmetric,"flux-compatible reconstruction currently requires RZ");
        // Stored Cartesian measures and the retained physical subquads need
        // only agree to coordinate roundoff. Neither is overwritten to make
        // them agree. Piola preserves net Q, not a fictitious constant netQ/V.
        const double extent=std::max({1.0,std::abs(lower.x)+nx*spacing.x,std::abs(lower.y)+ny*spacing.y});
        const double tolerance=256.0*std::numeric_limits<double>::epsilon()*extent/std::min(spacing.x,spacing.y);
        const double t[2]={0.5-0.5/std::sqrt(3.0),0.5+0.5/std::sqrt(3.0)};
        for(int j=0;j<ny;++j) for(int i=0;i<nx;++i) if(cellActive(i,j))
        {
            long double volume=0.0L;
            for(const auto& part:cellParts(i,j)) for(double u:t) for(double v:t)
            {
                const auto map=mapped(part,u,v);
                volume+=0.0625L*twoPi*map.point.x*map.determinant;
            }
            const double native=domain.cellVolumes()[domain.cellRowByIndex()[cells[i+nx*j]]];
            require(std::isfinite(volume) && volume>0.0L
                && std::abs(volume-native)<=tolerance*std::max(volume,static_cast<long double>(native)),
                "actual flux map volume differs from fitted FV volume beyond coordinate roundoff");
        }
        for(const auto& face:faces)
        {
            const int owner=face.negativeCell>=0?face.negativeCell:face.positiveCell;
            if(owner<0) continue;
            const auto ij=logicalCell[owner];
            const int side=face.negativeCell>=0?1:0;
            const int i=ij[0]+(!face.horizontal?side:0),j=ij[1]+(face.horizontal?side:0);
            const auto a=vertex(i,j),b=face.center,c=vertex(i+(face.horizontal?1:0),j+(face.horizontal?0:1));
            const auto measure=[](Vector2D x,Vector2D y) {
                return 0.5L*twoPi*(static_cast<long double>(x.x)+y.x)
                    *std::hypot(static_cast<long double>(y.x)-x.x,static_cast<long double>(y.y)-x.y);
            };
            const long double area=measure(a,b)+measure(b,c);
            require(std::isfinite(area) && area>=0.0L
                && std::abs(area-face.area)<=tolerance*std::max(area,static_cast<long double>(face.area)),
                "actual flux map face area differs from native FV area beyond coordinate roundoff");
        }
    }
    void initializeFluxReconstruction()
    {
        require(axisymmetric,"flux-compatible reconstruction currently requires RZ");
        std::vector<Vector2D> trace(vertices.size());
        std::vector<unsigned char> used(vertices.size(),0);
        for(int j=0;j<ny;++j) for(int i=0;i<nx;++i) if(cellActive(i,j))
            for(int y=0;y<2;++y) for(int x=0;x<2;++x)
            {
                const int k=i+x+(nx+1)*(j+y);
                if(!used[k]) { trace[k]=legacySample(vertices[k]); used[k]=1; }
            }
        fluxCells.resize(cells.size());
        for(int j=0;j<ny;++j) for(int i=0;i<nx;++i) if(cellActive(i,j))
        {
            auto& f=fluxCells[i+nx*j];
            f.axis=vertex(i,j).x==0.0 && vertex(i,j+1).x==0.0;
            for(int k=0;k<4;++k)
            {
                const int native=cellFace(i,j,k>=2,k%2);
                f.q[k]=checked(static_cast<long double>(faces[native].area)*value(native));
                if(faces[native].boundary.type==Type::NoSlipWall) f.walls|=1u<<k;
            }
            require(!f.axis || f.q[0]==0.0,"radial axis has nonzero native flux");
            require((f.walls&3u)==0 || (f.walls&12u)==0,
                "perpendicular no-slip walls require a compatible nonconstant divergence lifting");
            Vector2D velocity[4];
            for(int y=0;y<2;++y) for(int x=0;x<2;++x)
                velocity[x+2*y]=trace[i+x+(nx+1)*(j+y)];
            for(int y=0;y<2;++y) for(int x=0;x<2;++x)
            {
                const int k=x+2*y;
                const auto map=cellMap(i,j,x,y);
                const auto u=velocity[k];
                const double r=vertex(i+x,j+y).x;
                const double normalX=map.eta.y*u.x-map.eta.x*u.y;
                const double normalY=-map.xi.y*u.x+map.xi.x*u.y;
                const double baseX=f.axis?f.q[1]*x*x:lerp(f.q[0],f.q[1],x);
                const double baseY=(f.axis?2.0*x:1.0)*lerp(f.q[2],f.q[3],y);
                auto& jet=f.jet[k];
                jet.dx=baseY-twoPi*r*normalY;
                jet.dy=twoPi*r*normalX-baseX;
                // Differentiate a bilinear physical vertex-velocity estimate
                // through this actual Piola cofactor; retain metric derivatives.
                const auto du=velocity[1+2*y]-velocity[2*y];
                const auto dv=velocity[x+2]-velocity[x];
                const double desiredX=twoPi*(map.xi.x*normalX+r*(map.mixed.y*u.x
                    +map.eta.y*du.x-map.mixed.x*u.y-map.eta.x*du.y));
                const double desiredY=twoPi*(map.eta.x*normalY+r*(-map.mixed.y*u.x
                    -map.xi.y*dv.x+map.mixed.x*u.y+map.xi.x*dv.y));
                const double baseDx=f.axis?2.0*f.q[1]*x:f.q[1]-f.q[0];
                const double baseDy=(f.axis?2.0*x:1.0)*(f.q[3]-f.q[2]);
                jet.dxy=0.5*((desiredX-baseDx)+(baseDy-desiredY));
                if((f.walls&(1u<<x)) || (f.axis && x==0)) jet.dxy=baseDy;
                if(f.walls&(1u<<(2+y))) jet.dxy=-baseDx;
                if(f.axis)
                {
                    // psi(xi,eta)=Phi(xi^2,eta). A cubic in xi itself would
                    // introduce an O(r) term in u_z even on a rectangular
                    // mesh. Cubic Hermite in s retains radial evenness there,
                    // while all four zero corner values still preserve Q.
                    if(x==0)
                    {
                        require(r==0.0 && u.x==0.0 && map.eta.x==0.0 && map.xi.x>0.0,
                            "invalid retained radial-axis velocity or metric");
                        // Inverse weighted Piola loses finite u_z at r=0.
                        // Supply its analytic limit instead:
                        // u_z(0)=(Q_face-Phi_s(0))/(pi*r_xi^2).
                        const long double radial=map.xi.x;
                        jet.dx=checked(lerp(f.q[2],f.q[3],y)
                            -0.5L*twoPi*radial*radial*u.y);
                        jet.dy=0.0;
                        // Differentiate the same limit along eta, using the
                        // retained bilinear metric and shared axis trace.
                        const double axialDifference=velocity[2].y-velocity[0].y;
                        jet.dxy=checked(static_cast<long double>(f.q[3])-f.q[2]
                            -twoPi*radial*map.mixed.x*u.y
                            -0.5L*twoPi*radial*radial*axialDifference);
                    }
                    else
                    {
                        // d/dxi=2*xi*d/ds; xi=1 at the east corners.
                        jet.dx*=0.5;
                        jet.dxy*=0.5;
                    }
                    // On a horizontal no-slip wall Phi_eta=-Q_e*s. Its
                    // s-derivative is -Q_e at BOTH radial endpoints, including
                    // the axis where psi_xi_eta itself is identically zero.
                    if(f.walls&(1u<<(2+y))) jet.dxy=-f.q[1];
                }
                for(double a:{jet.dx,jet.dy,jet.dxy}) require(std::isfinite(a),"nonfinite flux streamfunction jet");
            }
        }
        shareTangentialJets();
    }
    void shareTangentialJets()
    {
        // The mixed correction jet does not enter a face's NORMAL trace or
        // integrated flux. Its independent cell estimates nevertheless enter
        // the TANGENTIAL trace and formerly gave two different velocities at
        // the same face. Share the physical axial-flux tangential derivative
        // at each retained vertex, rather than averaging query velocities.
        // For rectangular cells the shared quantity is
        // (baseY_eta - psi_xi_eta)/r_xi. The radial tangential trace then differs
        // only by the retained FV divergence residual; that residual is NOT
        // set to zero. All curl corrections still have identically zero div.
        std::vector<long double> sum(vertices.size(),0.0L);
        std::vector<unsigned> count(vertices.size(),0);
        for(int j=0;j<ny;++j) for(int i=0;i<nx;++i) if(cellActive(i,j))
        {
            const auto& f=fluxCells[i+nx*j];
            for(int y=0;y<2;++y) for(int x=0;x<2;++x)
            {
                const auto map=cellMap(i,j,x,y);
                const int k=i+x+(nx+1)*(j+y);
                if(f.axis && x==0)
                {
                    // At r=0 share the regular coefficient in Phi(s,eta),
                    // not a singular 1/r velocity. qEast+Phi_s_eta controls
                    // the odd radial trace on horizontal faces.
                    require(map.eta.y>0.0,"unresolved axial metric for shared axis jet");
                    sum[k]+=(static_cast<long double>(f.q[1])+f.jet[2*y].dxy)/map.eta.y;
                    ++count[k];
                    continue;
                }
                require(map.xi.x>0.0,"unresolved radial metric for shared tangential jet");
                const double baseDy=(f.axis?2.0*x:1.0)*(f.q[3]-f.q[2]);
                const double mixed=f.axis?2.0*f.jet[x+2*y].dxy:f.jet[x+2*y].dxy;
                sum[k]+=(static_cast<long double>(baseDy)-mixed)/map.xi.x;
                ++count[k];
            }
        }
        for(int j=0;j<ny;++j) for(int i=0;i<nx;++i) if(cellActive(i,j))
        {
            auto& f=fluxCells[i+nx*j];
            for(int y=0;y<2;++y) for(int x=0;x<2;++x)
            {
                // Whole-wall identities remain exact, including their actual
                // nonzero FV residual. Do not replace them by a zero residual.
                if((f.walls&(1u<<x)) || (f.walls&(1u<<(2+y)))) continue;
                const int k=i+x+(nx+1)*(j+y);
                require(count[k]>0,"missing shared tangential vertex jet");
                const long double shared=wallVertex(i+x,j+y)?0.0L:sum[k]/count[k];
                const auto map=cellMap(i,j,x,y);
                if(f.axis && x==0)
                {
                    f.jet[2*y].dxy=checked(map.eta.y*shared-static_cast<long double>(f.q[1]));
                    continue;
                }
                const double baseDy=(f.axis?2.0*x:1.0)*(f.q[3]-f.q[2]);
                const double mixed=checked(static_cast<long double>(baseDy)-map.xi.x*shared);
                f.jet[x+2*y].dxy=f.axis?0.5*mixed:mixed;
            }
        }
    }
    std::array<double,3> correction(const FluxCell& f,double u,double v) const
    {
        if(f.axis)
        {
            const Basis a(u*u),b(v);
            double ds=0.0,dy=0.0,dss=0.0;
            for(int y=0;y<2;++y) for(int x=0;x<2;++x)
            {
                const auto& c=f.jet[x+2*y];
                ds+=c.dx*a.dg[x]*b.h[y]+c.dy*a.dh[x]*b.g[y]+c.dxy*a.dg[x]*b.g[y];
                dy+=c.dx*a.g[x]*b.dh[y]+c.dy*a.h[x]*b.dg[y]+c.dxy*a.g[x]*b.dg[y];
                dss+=c.dx*a.ddg[x]*b.h[y]+c.dy*a.ddh[x]*b.g[y]+c.dxy*a.ddg[x]*b.g[y];
            }
            // Exact chain rule, including the nonzero axis second derivative.
            // With an orthogonal affine physical map this gives u_z even in
            // r and u_r odd. Retained non-Cartesian subquad roundoff is not
            // erased or advertised as exact physical-coordinate parity.
            return {{2.0*u*ds,dy,2.0*ds+4.0*u*u*dss}};
        }
        const Basis a(u),b(v);
        double dx=0.0,dy=0.0,dxx=0.0;
        for(int y=0;y<2;++y) for(int x=0;x<2;++x)
        {
            const auto& c=f.jet[x+2*y];
            dx+=c.dx*a.dg[x]*b.h[y]+c.dy*a.dh[x]*b.g[y]+c.dxy*a.dg[x]*b.g[y];
            dy+=c.dx*a.g[x]*b.dh[y]+c.dy*a.h[x]*b.dg[y]+c.dxy*a.g[x]*b.dg[y];
            dxx+=c.dx*a.ddg[x]*b.h[y]+c.dy*a.ddh[x]*b.g[y]+c.dxy*a.ddg[x]*b.g[y];
        }
        return {{dx,dy,dxx}};
    }
    std::array<int,3> fluxFaceOwner(int native) const
    {
        require(reconstruction==Reconstruction::RZFluxCompatible,
            "face flux charts require the flux-compatible sampler");
        require(native>=0 && static_cast<std::size_t>(native)<faces.size(), "invalid native flux face row");
        const auto& face=faces[native];
        require(face.nativeRow==native, "native flux face identity mismatch");
        // The negative-coordinate active owner is authoritative on an internal
        // face. Both callers therefore integrate the SAME retained trace, not
        // separately rounded reconstructions from opposite incident cells.
        const int owner=face.negativeCell>=0?face.negativeCell:face.positiveCell;
        require(owner>=0 && static_cast<std::size_t>(owner)<logicalCell.size(),
            "solid-only face has no active flux chart");
        const auto ij=logicalCell[owner];
        const int side=face.negativeCell>=0?1:0;
        require(cellActive(ij[0],ij[1]) && cells[ij[0]+nx*ij[1]]==owner
            && cellFace(ij[0],ij[1],face.horizontal,side)==native,
            "flux face owner is not its actual active incident cell");
        return {{ij[0],ij[1],side}};
    }
    FaceChart faceChart(int native) const
    {
        const auto owner=fluxFaceOwner(native);
        const int i=owner[0],j=owner[1],side=owner[2];
        const auto& face=faces[native];
        return {vertex(i+(face.horizontal?0:side),j+(face.horizontal?side:0)),
            face.center,vertex(i+(face.horizontal?1:side),j+(face.horizontal?side:1))};
    }
    void initializeNativeFluxTraces()
    {
        nativeFlux.resize(faces.size());
        for(const auto& face:faces)
        {
            if(face.negativeCell<0 && face.positiveCell<0) continue;
            const auto owner=fluxFaceOwner(face.nativeRow);
            const auto& f=fluxCells[owner[0]+nx*owner[1]];
            const int side=owner[2];
            const long double q=f.q[(face.horizontal?2:0)+side];
            const long double d0=face.horizontal?f.jet[2*side].dx:f.jet[side].dy;
            const long double d1=face.horizontal?f.jet[1+2*side].dx:f.jet[side+2].dy;
            nativeFlux[face.nativeRow].emplace(faceChart(face.nativeRow),
                face.horizontal,face.horizontal && f.axis,q,d0,d1);
        }
    }
    const HydroFittedFaceFlux& fluxTrace(int native) const
    {
        if(adaptive) return adaptive->face(native);
        require(reconstruction==Reconstruction::RZFluxCompatible && native>=0
            && static_cast<std::size_t>(native)<nativeFlux.size() && nativeFlux[native],
            "native face has no retained flux trace");
        return *nativeFlux[native];
    }
    long double integratedNormalFlux(int native,int half,double anchor0,double anchor1,
        double sBegin,double sEnd) const
    { return fluxTrace(native).integrate(half,anchor0,anchor1,sBegin,sEnd); }
    long double integratedNormalFlux(int native,int half,const LevelSetMappedSegment& segment) const
    { return fluxTrace(native).integrate(half,segment); }
    Vector2D fluxSample(Vector2D point, double* divergence) const
    {
        if(adaptive)
        {
            if(divergence) *divergence=adaptive->divergence(point);
            return adaptive->sample(point);
        }
        require(std::isfinite(point.x) && std::isfinite(point.y) && point.x>=0.0,"invalid RZ flux query");
        const auto ij=selectCell(point); const int i=ij[0],j=ij[1];
        std::array<double,2> uv;
        require(cellCoordinates(i,j,point,uv),"flux query left its actual cell");
        const double u=uv[0],v=uv[1];
        const auto& f=fluxCells[i+nx*j]; const auto map=cellMap(i,j,u,v);
        const auto psi=correction(f,u,v);
        const double net=(f.q[1]-f.q[0])+(f.q[3]-f.q[2]);
        if(point.x==0.0)
        {
            require(f.axis && u==0.0 && map.eta.x==0.0,"unresolved radial-axis flux map");
            if(divergence) *divergence=checked(net/(0.5L*twoPi*map.xi.x*map.determinant));
            return {0.0,checked((2.0L*lerp(f.q[2],f.q[3],v)-psi[2])/(twoPi*map.xi.x*map.xi.x))};
        }
        if(divergence) *divergence=checked((f.axis?2.0L*u:1.0L)*net/(twoPi*point.x*map.determinant));
        // These exact boundary zeros follow from the retained zero normal
        // flux and the compatible wall Hermite trace, not post-hoc pinning.
        if(((f.walls&1u)&&u==0.0)||((f.walls&2u)&&u==1.0)
            ||((f.walls&4u)&&v==0.0)||((f.walls&8u)&&v==1.0)) return {0.0,0.0};
        const double fx=(f.axis?f.q[1]*u*u:lerp(f.q[0],f.q[1],u))+psi[1];
        const double fy=(f.axis?2.0*u:1.0)*lerp(f.q[2],f.q[3],v)-psi[0];
        const long double measure=static_cast<long double>(twoPi)*point.x*map.determinant;
        return {checked((static_cast<long double>(map.xi.x)*fx+map.eta.x*fy)/measure),
                checked((static_cast<long double>(map.xi.y)*fx+map.eta.y*fy)/measure)};
    }
    Vector2D sample(Vector2D point) const
    { return reconstruction==Reconstruction::RZFluxCompatible?fluxSample(point,nullptr):legacySample(point); }
};

HydroFittedVelocitySampler::HydroFittedVelocitySampler(const QuadtreeMesh& mesh,
    const HydroFittedDomain& domain, const VEdgeField<double>& ux, const HEdgeField<double>& uy, Reconstruction method)
    : impl_(std::make_unique<Impl>(mesh, domain, ux, uy, method)) {}
HydroFittedVelocitySampler::~HydroFittedVelocitySampler() = default;
HydroFittedVelocitySampler::HydroFittedVelocitySampler(const QuadtreeMesh& mesh,
    const NativeLeafDomain& domain, const VEdgeField<double>& ux, const HEdgeField<double>& uy,
    const std::vector<VertexTrace>& traces)
    : impl_(std::make_unique<Impl>(mesh,domain,ux,uy,traces)) {}
HydroFittedVelocitySampler::HydroFittedVelocitySampler(const QuadtreeMesh& mesh,
    const HydroFittedDomain& domain, const VEdgeField<double>& ux, const HEdgeField<double>& uy,
    const std::vector<VertexTrace>& traces)
    : impl_(std::make_unique<Impl>(mesh,domain,ux,uy,traces)) {}
HydroFittedVelocitySampler::HydroFittedVelocitySampler(HydroFittedVelocitySampler&&) noexcept = default;
HydroFittedVelocitySampler& HydroFittedVelocitySampler::operator=(HydroFittedVelocitySampler&&) noexcept = default;
Vector2D HydroFittedVelocitySampler::sample(Vector2D point) const
{
    require(static_cast<bool>(impl_), "moved-from sampler");
    return impl_->sample(point);
}
double HydroFittedVelocitySampler::divergence(Vector2D point) const
{
    require(impl_ && impl_->reconstruction==Reconstruction::RZFluxCompatible,"divergence requires the flux-compatible sampler");
    double result=0.0; impl_->fluxSample(point,&result); return result;
}
HydroFittedVelocitySampler::FirstJet HydroFittedVelocitySampler::firstJet(Vector2D point) const
{
    require(impl_ && impl_->adaptive,
            "firstJet requires the piecewise native face/streamfunction reconstruction");
    return impl_->adaptive->firstJet(point);
}
HydroFittedVelocitySampler::FirstJet HydroFittedVelocitySampler::volumeFirstJet(Vector2D point) const
{
    require(impl_ && impl_->adaptive,
            "volumeFirstJet requires the piecewise native face/streamfunction reconstruction");
    return impl_->adaptive->volumeFirstJet(point);
}
HydroFittedVelocitySampler::FaceChart HydroFittedVelocitySampler::faceChart(int nativeRow) const
{
    require(static_cast<bool>(impl_), "moved-from sampler");
    return impl_->fluxTrace(nativeRow).chart();
}
long double HydroFittedVelocitySampler::normalFluxDensity(int nativeRow,int half,bool fromEnd,long double offset) const
{
    require(static_cast<bool>(impl_), "moved-from sampler");
    return impl_->fluxTrace(nativeRow).density(half,fromEnd,offset);
}
Quantum::Geometry::UnitPolynomialRoots HydroFittedVelocitySampler::normalFluxRoots(int nativeRow,int half) const
{
    require(static_cast<bool>(impl_),"moved-from sampler");
    return impl_->fluxTrace(nativeRow).densityRoots(half);
}
std::vector<Quantum::Geometry::UnitPolynomialTerm> HydroFittedVelocitySampler::normalFluxPolynomial(int nativeRow,int half) const
{
    require(static_cast<bool>(impl_),"moved-from sampler");
    return impl_->fluxTrace(nativeRow).densityPolynomial(half);
}
long double HydroFittedVelocitySampler::integratedNormalFlux(int nativeRow,int half,double tBegin,double tEnd) const
{
    return integratedNormalFlux(nativeRow,half,0.0,1.0,tBegin,tEnd);
}
long double HydroFittedVelocitySampler::integratedNormalFlux(int nativeRow,int half,double anchor0,double anchor1,
    double sBegin,double sEnd) const
{
    require(static_cast<bool>(impl_), "moved-from sampler");
    return impl_->integratedNormalFlux(nativeRow,half,anchor0,anchor1,sBegin,sEnd);
}

long double HydroFittedVelocitySampler::integratedNormalFlux(int nativeRow,int half,
    const LevelSetMappedSegment& segment) const
{
    require(static_cast<bool>(impl_), "moved-from sampler");
    return impl_->integratedNormalFlux(nativeRow,half,segment);
}
