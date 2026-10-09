#include "Modules/Hydrodynamics/Detail/HydroFittedBoundaryMomentum.h"
#include "Modules/Hydrodynamics/HydroFittedDomain.h"
#include "Core/Base/Interpolation.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/CartesianMeasure.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <string>
#include <tuple>

namespace HydroFittedDetail
{
namespace
{
using Face = HydroFittedDomain::NativeFace;
using Type = HydroFittedDomain::BoundaryType;
void require(bool condition, const char* message)
{
    if (!condition) throw std::invalid_argument(std::string("[HydroFittedBoundaryMomentum] ") + message);
}
double at(Vector2D point, int axis) { return axis == 0 ? point.x : point.y; }
bool pressure(Type type) { return type == Type::PressureReservoir || type == Type::PressureOutlet; }
const Cell& cellAt(const QuadtreeMesh& mesh, int index) { return mesh.cells(mesh.cellIndexToLabel[index]); }
const Edge& edgeAt(const QuadtreeMesh& mesh, const Face& face)
{
    return face.horizontal ? mesh.hEdges(face.edgeLabel) : mesh.vEdges(face.edgeLabel);
}
int nativeAt(const QuadtreeMesh& mesh, int component, int label)
{
    if (label < 0) return -1;
    return component == 0 ? mesh.vEdges(label).index() : mesh.vEdgeIndexNum() + mesh.hEdges(label).index();
}
int singleFace(const QuadtreeMesh& mesh, const Cell& cell, int component, int side)
{
    const auto* labels = component == 0 ? cell.cellToVEdge[side] : cell.cellToHEdge[side];
    return labels[1] < 0 ? nativeAt(mesh, component, labels[0]) : -1;
}
void add(MomentumTrace& target, const MomentumTrace& source, double factor = 1.0)
{
    for (const auto& entry : source)
    {
        auto found = std::find_if(target.begin(), target.end(), [&](const auto& old) { return old.first == entry.first; });
        if (found == target.end()) target.emplace_back(entry.first, factor * entry.second);
        else found->second += factor * entry.second;
    }
}
MomentumTrace difference(int positive, int negative, double spacing)
{
    require(std::isfinite(spacing) && spacing > 0.0, "nonpositive native derivative distance");
    return {{positive, 1.0 / spacing}, {negative, -1.0 / spacing}};
}
std::vector<int> supportCells(const QuadtreeMesh& mesh, const HydroFittedDomain& domain,
    const std::vector<int>& seed)
{
    return nativeTraceSupport(mesh, domain.geometry(), seed);
}
using Fit = NativeVelocityTrace;
Fit fit(const QuadtreeMesh& mesh, const HydroFittedDomain& domain,
    const std::vector<int>& cells, int component, Vector2D point)
{
    return nativeVelocityTrace(mesh, domain.geometry(), cells, component, point);
}
std::vector<int> owners(const Face& a, const Face& b)
{
    return {a.negativeCell, a.positiveCell, b.negativeCell, b.positiveCell};
}
int commonVertex(const QuadtreeMesh& mesh, const Face& a, const Face& b)
{
    for (int x : edgeAt(mesh, a).edgeToVert) for (int y : edgeAt(mesh, b).edgeToVert)
        if (x == y) return mesh.verts(x).index();
    return -1;
}
bool sameLevel(const QuadtreeMesh& mesh, const Face& a, const Face& b)
{
    int level = -1;
    for (int index : owners(a, b)) if (index >= 0)
    {
        if (level < 0) level = cellAt(mesh, index).level();
        else if (level != cellAt(mesh, index).level()) return false;
    }
    return level >= 0;
}
MomentumTrace inwardTrace(const QuadtreeMesh& mesh, const Face& face)
{
    const int component = face.horizontal ? 1 : 0;
    const int owner = face.negativeCell >= 0 ? face.negativeCell : face.positiveCell;
    const auto& cell = cellAt(mesh, owner);
    const int side = face.outwardSign > 0.0 ? 0 : 1;
    const auto* labels = component == 0 ? cell.cellToVEdge[side] : cell.cellToHEdge[side];
    MomentumTrace result;
    for (int part = 0; part < 2; ++part) if (labels[part] >= 0)
    {
        const auto& edge = component == 0 ? mesh.vEdges(labels[part]) : mesh.hEdges(labels[part]);
        result.emplace_back(nativeAt(mesh, component, labels[part]), edge.length() / at(cell.size(), 1 - component));
    }
    require(!result.empty(), "pressure continuation lacks an inward native face");
    return result;
}
MomentumTrace incomingTangent(const QuadtreeMesh& mesh, const HydroFittedDomain& domain,
    const Face& tangential, const Face& boundary)
{
    const int normal = boundary.horizontal ? 1 : 0, tangent = 1 - normal;
    const double plane = at(boundary.center, normal);
    const auto cells = supportCells(mesh, domain, owners(tangential, boundary));
    std::vector<int> rows;
    for (int index : cells) for (int side = 0; side < 2; ++side)
    {
        const int row = singleFace(mesh, cellAt(mesh, index), normal, side);
        if (row < 0) continue;
        const auto& face = domain.faces()[row];
        if (pressure(face.boundary.type) && face.outwardSign == boundary.outwardSign && at(face.center, normal) == plane)
            rows.push_back(row);
    }
    std::sort(rows.begin(), rows.end(), [&](int a, int b) { return at(domain.faces()[a].center, tangent) < at(domain.faces()[b].center, tangent); });
    rows.erase(std::unique(rows.begin(), rows.end()), rows.end());
    require(rows.size() >= 2, "zero-shear inflow trace needs two real pressure-face samples");
    const double target = at(tangential.center, tangent);
    auto upper = std::lower_bound(rows.begin(), rows.end(), target,
        [&](int row, double x) { return at(domain.faces()[row].center, tangent) < x; });
    std::size_t hi = static_cast<std::size_t>(upper - rows.begin());
    hi = std::max<std::size_t>(1, std::min(hi, rows.size() - 1));
    const auto& a = domain.faces()[rows[hi - 1]]; const auto& b = domain.faces()[rows[hi]];
    require(at(mesh.verts(edgeAt(mesh, a).edgeToVert[1]).coord(), tangent)
        == at(mesh.verts(edgeAt(mesh, b).edgeToVert[0]).coord(), tangent),
        "zero-shear inflow trace crosses a gap in the pressure boundary");
    MomentumTrace result{{tangential.nativeRow, 1.0}};
    // u_t(boundary)=u_t(interior)-signed_distance*d_t(u_n), the original
    // zero-shear continuation with true physical sampling distances.
    add(result, difference(b.nativeRow, a.nativeRow, at(b.center, tangent) - at(a.center, tangent)),
        -(plane - at(tangential.center, normal)));
    return result;
}
}

std::shared_ptr<const BoundaryMomentum> BoundaryMomentum::build(const QuadtreeMesh& mesh,
    const HydroFittedDomain& domain, std::vector<unsigned char>& nativeInterior,
    std::vector<unsigned char> boundaryRows)
{
    auto result = std::shared_ptr<BoundaryMomentum>(new BoundaryMomentum);
    result->duals_ = MomentumDuals::build(mesh, domain);
    result->selected_ = std::move(boundaryRows);
    result->prepare(mesh, domain, nativeInterior);
    return result;
}

void BoundaryMomentum::prepare(const QuadtreeMesh& mesh, const HydroFittedDomain& domain,
    std::vector<unsigned char>& nativeInterior)
{
    // Combine contiguous half-pieces belonging to one pair of velocity CVs.
    // A regular full shear side therefore retains the old shared vertex rule.
    std::map<std::tuple<int, int, int, int>, int> byPair;
    for (std::size_t id = 0; id < duals_.segments().size(); ++id)
    {
        const auto& part = duals_.segments()[id];
        const bool exterior = part.negativeNative < 0 || part.positiveNative < 0;
        const auto key = std::make_tuple(part.negativeNative, part.positiveNative, part.normalAxis,
            exterior ? part.primalNative : -1);
        auto inserted = byPair.emplace(key, static_cast<int>(connections_.size()));
        if (inserted.second)
        {
            Connection connection;
            connection.negative = part.negativeNative; connection.positive = part.positiveNative;
            connection.normal = part.normalAxis;
            connection.kind = exterior ? Kind::Exterior : Kind::General;
            connections_.push_back(std::move(connection));
        }
        connections_[inserted.first->second].parts.push_back(static_cast<int>(id));
    }
    for (auto& connection : connections_)
    {
        if (connection.kind == Kind::Exterior) continue;
        const auto& a = domain.faces()[connection.negative]; const auto& b = domain.faces()[connection.positive];
        const int component = a.horizontal ? 1 : 0;
        if (!sameLevel(mesh, a, b)) continue;
        if (at(a.center, 1 - connection.normal) != at(b.center, 1 - connection.normal)) continue;
        if (component == connection.normal)
        {
            const int owner = a.positiveCell;
            if (owner < 0 || owner != b.negativeCell) continue;
            const auto& cell = cellAt(mesh, owner);
            if (singleFace(mesh, cell, component, 0) != a.nativeRow || singleFace(mesh, cell, component, 1) != b.nativeRow) continue;
            connection.kind = Kind::Normal; connection.materialCell = owner;
        }
        else
        {
            const int vertex = commonVertex(mesh, a, b);
            if (vertex < 0) continue;
            if (pressure(a.boundary.type) && pressure(b.boundary.type) && a.outwardSign == b.outwardSign)
            {
                connection.kind = Kind::PressurePair; connection.materialVertex = vertex;
                continue;
            }
            if (a.negativeCell < 0 || a.positiveCell < 0 || b.negativeCell < 0 || b.positiveCell < 0) continue;
            const int low = singleFace(mesh, cellAt(mesh, a.negativeCell), connection.normal, 1);
            const int high = singleFace(mesh, cellAt(mesh, a.positiveCell), connection.normal, 1);
            if (low < 0 || high < 0 || domain.faces()[low].boundary.type != Type::Interior
                || domain.faces()[high].boundary.type != Type::Interior) continue;
            connection.kind = Kind::Shear; connection.materialVertex = vertex;
        }
    }
    // Never splice a physical flux into the middle of an old virtual stencil.
    // Expand only this connected affected component, ending at regular sides.
    bool changed = true;
    while (changed)
    {
        changed = false;
        for (const auto& connection : connections_)
        {
            if (connection.negative < 0 || connection.positive < 0) continue;
            for (const auto pair : {std::make_pair(connection.negative, connection.positive),
                    std::make_pair(connection.positive, connection.negative)})
                if (selected_[pair.first] && !selected_[pair.second] && domain.faces()[pair.second].velocityRow >= 0
                    && (nativeInterior[pair.second] || connection.kind == Kind::General))
                {
                    selected_[pair.second] = 1; changed = true;
                }
        }
    }
    for (std::size_t row = 0; row < selected_.size(); ++row) if (selected_[row]) nativeInterior[row] = 0;
    connections_.erase(std::remove_if(connections_.begin(), connections_.end(), [&](const Connection& connection) {
        return (connection.negative < 0 || !selected_[connection.negative])
            && (connection.positive < 0 || !selected_[connection.positive]);
    }), connections_.end());

    std::map<int, Fit> faceFit;
    const auto derivatives = [&](int row) -> const Fit& {
        auto found = faceFit.find(row);
        if (found == faceFit.end())
        {
            const auto& face = domain.faces()[row];
            found = faceFit.emplace(row, fit(mesh, domain, supportCells(mesh, domain,
                {face.negativeCell, face.positiveCell}), face.horizontal ? 1 : 0, face.center)).first;
        }
        return found->second;
    };
    inward_.resize(selected_.size()); radialSlope_.resize(selected_.size()); ownerVolumes_.resize(selected_.size());
    for (const auto& face : domain.faces()) if (selected_[face.nativeRow])
    {
        const auto& volume = duals_.volumes()[face.nativeRow];
        for (int side = 0; side < 2; ++side)
        {
            const int owner = side == 0 ? face.negativeCell : face.positiveCell;
            if (owner < 0) continue;
            const auto center = cellAt(mesh, owner).centerCoord();
            if (face.horizontal)
                ownerVolumes_[face.nativeRow][side] = CartesianMeasure::twoPi * 0.5 * (volume.upper.x + volume.lower.x)
                    * (volume.upper.x - volume.lower.x) * std::abs(face.center.y - center.y);
            else
                ownerVolumes_[face.nativeRow][side] = 0.5 * CartesianMeasure::twoPi * std::abs(face.center.x - center.x)
                    * (face.center.x + center.x) * (volume.upper.y - volume.lower.y);
        }
        if (pressure(face.boundary.type)) inward_[face.nativeRow] = inwardTrace(mesh, face);
        if (!face.horizontal)
        {
            if (pressure(face.boundary.type))
            {
                const double distance = -face.outwardSign * face.normalSpacing;
                add(radialSlope_[face.nativeRow], inward_[face.nativeRow], 1.0 / distance);
                add(radialSlope_[face.nativeRow], {{face.nativeRow, 1.0}}, -1.0 / distance);
            }
            else radialSlope_[face.nativeRow] = derivatives(face.nativeRow).gradient[0];
        }
    }
    for (auto& connection : connections_)
    {
        const int row = connection.negative >= 0 ? connection.negative : connection.positive;
        const auto& a = domain.faces()[row];
        const int component = a.horizontal ? 1 : 0, normal = connection.normal;
        for (int id : connection.parts) connection.area += duals_.segments()[id].measure;
        require(connection.area >= 0.0 && std::isfinite(connection.area), "invalid shared momentum area");
        if (connection.kind == Kind::Exterior)
        {
            const auto& boundary = domain.faces()[duals_.segments()[connection.parts[0]].primalNative];
            connection.advector = {{boundary.nativeRow, 1.0}};
            if (pressure(boundary.boundary.type))
                connection.incoming = component == normal ? inwardTrace(mesh, boundary) : incomingTangent(mesh, domain, a, boundary);
            else if (boundary.boundary.type == Type::NoSlipWall && component != normal)
            {
                const double distance = boundary.outwardSign * (at(boundary.center, normal) - at(a.center, normal));
                require(distance > 0.0, "wall shear lacks a physical inward sample");
                connection.strain = {{a.nativeRow, -boundary.outwardSign / distance}};
            }
            continue;
        }
        const auto& b = domain.faces()[connection.positive];
        if (connection.kind == Kind::Normal)
        {
            const auto& cell = cellAt(mesh, connection.materialCell);
            connection.area = 2.0 * Pi * (component == 0 ? cell.centerCoord().x : a.center.x)
                * at(cell.size(), 1 - component);
            connection.strain = difference(b.nativeRow, a.nativeRow, 0.5 * at(cell.size(), component));
            connection.advector = {{a.nativeRow, 0.5}, {b.nativeRow, 0.5}};
        }
        else if (connection.kind == Kind::Shear)
        {
            const auto& master = cellAt(mesh, a.positiveCell);
            const auto& vertex = mesh.verts(mesh.vertIndexToLabel[connection.materialVertex]);
            connection.area = component == 0 ? 2.0 * Pi * a.center.x * master.size().x
                : 2.0 * Pi * vertex.coord().x * master.size().y;
            const int low = singleFace(mesh, cellAt(mesh, a.negativeCell), normal, 1);
            const int high = singleFace(mesh, master, normal, 1);
            connection.strain = difference(b.nativeRow, a.nativeRow, a.tangentSpacing);
            connection.currentStrain = difference(b.nativeRow, a.nativeRow, at(b.center, normal) - at(a.center, normal));
            const auto cross = difference(high, low, a.normalSpacing);
            add(connection.strain, cross); add(connection.currentStrain, cross);
            connection.advector = {{high, 0.5}, {low, 0.5}};
        }
        else if (connection.kind == Kind::PressurePair)
        {
            require(connection.parts.size() == 1, "pressure half-duals have a split regular side");
            connection.advector = {{duals_.segments()[connection.parts[0]].primalNative, 1.0}};
        }
        else
        {
            const double distance = at(b.center, normal) - at(a.center, normal);
            const double offset = at(b.center, 1 - normal) - at(a.center, 1 - normal);
            connection.strain = difference(b.nativeRow, a.nativeRow, distance);
            if (offset != 0.0)
            {
                add(connection.strain, derivatives(a.nativeRow).gradient[1 - normal], -0.5 * offset / distance);
                add(connection.strain, derivatives(b.nativeRow).gradient[1 - normal], -0.5 * offset / distance);
            }
            Vector2D center{0.0, 0.0};
            for (int id : connection.parts)
            {
                const auto& part = duals_.segments()[id];
                double r = 0.5 * (part.begin.x + part.end.x);
                if (normal == 1)
                    r = (2.0 / 3.0) * (part.begin.x * part.begin.x + part.begin.x * part.end.x + part.end.x * part.end.x)
                        / (part.begin.x + part.end.x);
                center.x += part.measure * r;
                center.y += part.measure * 0.5 * (part.begin.y + part.end.y);
            }
            require(connection.area > 0.0, "interior momentum connection has zero area");
            center.x /= connection.area; center.y /= connection.area;
            const auto normalFit = fit(mesh, domain, supportCells(mesh, domain, owners(a, b)), normal, center);
            connection.advector = normalFit.value;
            if (component == normal)
                for (auto& coefficient : connection.strain) coefficient.second *= 2.0;
            else add(connection.strain, normalFit.gradient[component]);
        }
    }
}
}
