#include "Core/Geometry/ActiveQuarticHermiteTransfer.h"
#include "Core/Geometry/Detail/LevelSetHermiteQuarticLeafData.h"
#include "Core/Geometry/Detail/NativeHermiteTransferSampling.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/CartesianLeafMesh.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>

namespace ActiveQuarticHermiteTransfer
{
namespace
{
using Patch = Quantum::Geometry::ActiveQuarticHermitePatch;
namespace Scalar = Quantum::Geometry::Detail::LevelSetHermiteScalar;
namespace Exact = Quantum::Geometry::Detail::ExactExpansion;
using Source = Scalar::QuarticLeafData;
using PointKey = std::array<double, 2>;
using BoxKey = std::array<double, 4>;
using EdgeKey = std::array<int, 2>;

void require(bool valid, const std::string& reason)
{
    if (!valid) throw std::runtime_error("[ActiveQuarticHermiteTransfer] " + reason);
}
PointKey key(Vector2D p) { return {{p.x, p.y}}; }
BoxKey key(const std::array<Vector2D, 2>& b) { return {{b[0].x, b[0].y, b[1].x, b[1].y}}; }
std::array<Vector2D, 2> bounds(const QuadtreeMesh& mesh, int c)
{
    const auto& cell = mesh.cells(mesh.cellIndexToLabel[c]);
    return {{mesh.verts(cell.cellToVert[0]).coord(), mesh.verts(cell.cellToVert[2]).coord()}};
}
bool overlaps(const std::array<Vector2D, 2>& a, const std::array<Vector2D, 2>& b)
{
    return std::max(a[0].x, b[0].x) < std::min(a[1].x, b[1].x)
        && std::max(a[0].y, b[0].y) < std::min(a[1].y, b[1].y);
}
void validateMask(const std::vector<std::uint8_t>& mask, int size, const char* name)
{
    require(mask.size() == static_cast<std::size_t>(size), std::string(name) + " has wrong native-cell count");
    for (auto bit : mask) require(bit <= 1, std::string(name) + " must be binary");
}

struct Layout
{
    std::vector<std::uint8_t> usedVertices;
    // Only identities are enumerated. Hanging data are still resolved solely
    // by QuarticLeafData; this is not another constraint/interpolation graph.
    std::map<EdgeKey, bool> sides; // true = horizontal
    Layout(const QuadtreeMesh& mesh, const std::vector<std::uint8_t>& active)
        : usedVertices(mesh.vertIndexNum(), 0)
    {
        constexpr int pairs[4][2] = {{0,1}, {3,2}, {0,3}, {1,2}};
        for (int c = 0; c < mesh.cellIndexNum(); ++c)
        {
            if (!active[c]) continue;
            const auto& cell = mesh.cells(mesh.cellIndexToLabel[c]);
            for (int label : cell.cellToVert) usedVertices[mesh.verts(label).index()] = 1;
            for (int side = 0; side < 4; ++side)
                sides.emplace(EdgeKey{{mesh.verts(cell.cellToVert[pairs[side][0]]).index(),
                    mesh.verts(cell.cellToVert[pairs[side][1]]).index()}}, side < 2);
        }
    }
};

using DatumSampler = Quantum::Geometry::Detail::NativeHermiteTransfer::DatumSampler<Source>;}

ProtectedPolynomialChange::ProtectedPolynomialChange(int sourceCell, int targetCell)
    : std::runtime_error("[ActiveQuarticHermiteTransfer] protected single-phase polynomial changed through shared constraints; source cell="
        + std::to_string(sourceCell) + "; target cell=" + std::to_string(targetCell))
{
}

std::vector<std::uint8_t> topologyProtection(const QuadtreeMesh& mesh,
    const std::vector<std::uint8_t>& preservedCells)
{
    validateMask(preservedCells, mesh.cellIndexNum(), "preserved scalar cells");
    validateCartesianLeafMesh(mesh);
    auto carriers = preservedCells;
    std::vector<int> pending;
    for (int i = 0; i < mesh.cellIndexNum(); ++i) if (carriers[i]) pending.push_back(i);
    for (std::size_t k = 0; k < pending.size(); ++k)
        for (int vertex : mesh.cells(mesh.cellIndexToLabel[pending[k]]).cellToVert)
        {
            const int coarse = mesh.verts(vertex).supCell;
            if (coarse < 0) continue;
            const int index = mesh.cells(coarse).index();
            if (!carriers[index]) { carriers[index] = 1; pending.push_back(index); }
        }
    auto result = carriers;
    // Keep incidence around carrier corners so an independent datum cannot
    // become a hanging trace of a newly coarsened neighbor (or vice versa).
    for (int index : pending)
        for (int vertex : mesh.cells(mesh.cellIndexToLabel[index]).cellToVert)
            for (int label : mesh.verts(vertex).vertToCell)
                if (label >= 0) result[mesh.cells(label).index()] = 1;
    return result;
}

Result prepare(const QuadtreeMesh& sourceMesh, const Patch& sourcePatch,
    const QuadtreeMesh& targetMesh, const std::vector<std::uint8_t>& targetActive,
    bool radial, const std::vector<std::uint8_t>& protectedSource)
{
    Exact::validateEnvironment();
    require(key(sourceMesh.minCoord) == key(targetMesh.minCoord)
        && key(sourceMesh.maxCoord) == key(targetMesh.maxCoord)
        && sourceMesh.Nx == targetMesh.Nx && sourceMesh.Ny == targetMesh.Ny,
        "remeshing requires the same background domain and lattice");
    validateMask(targetActive, targetMesh.cellIndexNum(), "target activity");
    if (!protectedSource.empty()) validateMask(protectedSource, sourceMesh.cellIndexNum(), "protected support");
    const Source source(sourceMesh, sourcePatch, radial); // validates the complete old authority
    const Layout oldLayout(sourceMesh, sourcePatch.activeCells), nextLayout(targetMesh, targetActive);
    Result result;
    result.sourceToTargetCell.assign(sourceMesh.cellIndexNum(), -1);
    result.sourceToTargetVertex.assign(sourceMesh.vertIndexNum(), -1);
    std::map<BoxKey, int> oldCells, nextCells;
    std::vector<int> sourceSigns(sourceMesh.cellIndexNum(), 0), sameBounds(targetMesh.cellIndexNum(), -1);
    std::vector<std::uint8_t> mustPreserve(sourceMesh.cellIndexNum(), 0);
    for (int c = 0; c < sourceMesh.cellIndexNum(); ++c) oldCells.emplace(key(source.leafBounds(c)), c);
    for (int c = 0; c < targetMesh.cellIndexNum(); ++c)
    {
        const auto box = bounds(targetMesh, c);
        nextCells.emplace(key(box), c);
        const auto previous = oldCells.find(key(box));
        if (previous != oldCells.end()) sameBounds[c] = previous->second;
        bool covered = false;
        for (int old : source.leafCandidates(box[0], box[1]))
        {
            if (!overlaps(box, source.leafBounds(old))) continue;
            covered = true;
            require(source.cellActive(old) == (targetActive[c] != 0),
                "remeshing changes fixed active/solid coverage at target cell=" + std::to_string(c));
        }
        require(covered, "target leaf is not covered by the source mesh");
    }
    for (int c = 0; c < sourceMesh.cellIndexNum(); ++c)
    {
        if (!source.cellActive(c)) continue;
        sourceSigns[c] = source.cellSignExceptVertices(c, 0);
        mustPreserve[c] = sourceSigns[c] == 0 || (!protectedSource.empty() && protectedSource[c]);
        if (sourceSigns[c] == 0) ++result.report.protectedInterfaceCells;
        if (mustPreserve[c])
            require(nextCells.count(key(source.leafBounds(c))) != 0,
                "interface/support leaf would be refined or coarsened; enlarge the advance refinement band; source cell="
                    + std::to_string(c));
    }

    std::map<PointKey, int> oldVertices, nextVertices;
    for (int k = 0; k < sourceMesh.vertIndexNum(); ++k)
        if (oldLayout.usedVertices[k]) oldVertices.emplace(key(sourceMesh.verts(sourceMesh.vertIndexToLabel[k]).coord()), k);
    for (int k = 0; k < targetMesh.vertIndexNum(); ++k)
        nextVertices.emplace(key(targetMesh.verts(targetMesh.vertIndexToLabel[k]).coord()), k);
    for (int k = 0; k < sourceMesh.vertIndexNum(); ++k)
    {
        const auto found = nextVertices.find(key(sourceMesh.verts(sourceMesh.vertIndexToLabel[k]).coord()));
        if (found != nextVertices.end()) result.sourceToTargetVertex[k] = found->second;
    }

    const double nan = std::numeric_limits<double>::quiet_NaN();
    auto& patch = result.patch;
    patch.activeCells = targetActive;
    patch.nodes.assign(targetMesh.vertIndexNum(), {nan, nan, nan, nan});
    patch.cellCenters.assign(targetMesh.cellIndexNum(), nan);
    const DatumSampler sampler(source);
    for (int k = 0; k < targetMesh.vertIndexNum(); ++k)
    {
        const auto& vertex = targetMesh.verts(targetMesh.vertIndexToLabel[k]);
        if (!nextLayout.usedVertices[k] || vertex.supCell >= 0) continue;
        const auto found = oldVertices.find(key(vertex.coord()));
        if (found != oldVertices.end())
        {
            const int old = found->second;
            if (sourceMesh.verts(sourceMesh.vertIndexToLabel[old]).supCell < 0)
            {
                patch.nodes[k] = sourcePatch.nodes[old];
                ++result.report.copiedVertices;
            }
            else
            {
                const auto jet = source.nodeJet(old); // old shared coarse trace, not ignored input payload
                patch.nodes[k] = {jet[0], jet[1], jet[2], jet[3]};
            }
        }
        else
        {
            const auto point = sampler.locate(vertex.coord(), vertex.coord());
            patch.nodes[k] = {sampler.value(point), sampler.value(point, 1, 0),
                sampler.value(point, 0, 1), sampler.value(point, 1, 1)};
        }
    }
    std::map<BoxKey, const Patch::EdgeMidpoint*> oldMidpoints;
    for (const auto& edge : sourcePatch.edgeMidpoints)
        if (source.independentMidpoint(edge.vertices))
        {
            const auto a = sourceMesh.verts(sourceMesh.vertIndexToLabel[edge.vertices[0]]).coord();
            const auto b = sourceMesh.verts(sourceMesh.vertIndexToLabel[edge.vertices[1]]).coord();
            oldMidpoints.emplace(BoxKey{{a.x, a.y, b.x, b.y}}, &edge);
        }
    for (const auto& entry : nextLayout.sides)
    {
        const auto a = targetMesh.verts(targetMesh.vertIndexToLabel[entry.first[0]]).coord();
        const auto b = targetMesh.verts(targetMesh.vertIndexToLabel[entry.first[1]]).coord();
        Patch::EdgeMidpoint edge{entry.first, nan, nan};
        const auto old = oldMidpoints.find(BoxKey{{a.x, a.y, b.x, b.y}});
        if (old != oldMidpoints.end())
        {
            edge.value = old->second->value;
            edge.normalDerivative = old->second->normalDerivative;
            ++result.report.copiedMidpoints;
        }
        else
        {
            const auto point = sampler.locate(a, b);
            edge.value = sampler.value(point);
            edge.normalDerivative = sampler.value(point, entry.second ? 0 : 1, entry.second ? 1 : 0);
        }
        patch.edgeMidpoints.push_back(edge); // dependent subedge payload is ignored by the scalar owner
    }
    for (int c = 0; c < targetMesh.cellIndexNum(); ++c)
    {
        if (!targetActive[c]) continue;
        if (sameBounds[c] >= 0)
        {
            patch.cellCenters[c] = sourcePatch.cellCenters[sameBounds[c]];
            ++result.report.copiedCenters;
        }
        else
        {
            const auto box = bounds(targetMesh, c);
            patch.cellCenters[c] = sampler.value(sampler.locate(box[0], box[1]));
        }
    }
    for (const auto& old : sourcePatch.cornerSectors)
    {
        auto sector = old;
        sector.cellIndex = nextCells.at(key(source.leafBounds(old.cellIndex)));
        sector.vertexIndex = result.sourceToTargetVertex.at(old.vertexIndex);
        require(sector.vertexIndex >= 0, "declared corner-sector vertex was removed");
        patch.cornerSectors.push_back(sector);
    }

    // Rebuild the SAME shared constraints and sector/axis certificates. Exact
    // comparison is after constraints, so a new hanging dependency cannot
    // silently change a retained interface polynomial through its neighbors.
    const Source candidate(targetMesh, patch, radial);
    // Keep the same public ownership convention as native reconstruction.
    // Geometry ignores fine-side payload, but transport/correction consumers
    // recognize NaN as dependent data, never as an independent unknown.
    // Ask the existing scalar authority; do not build another constraint graph.
    for (auto& edge : patch.edgeMidpoints)
        if (!candidate.independentMidpoint(edge.vertices))
            edge.value = edge.normalDerivative = nan;
    for (int c = 0; c < targetMesh.cellIndexNum(); ++c)
    {
        if (!targetActive[c]) continue;
        const int old = sameBounds[c];
        if (old >= 0 && source.identicalCell(old, candidate, c))
        {
            result.sourceToTargetCell[old] = c;
            ++result.report.identicalActiveCells;
            continue;
        }
        const int sign = candidate.cellSignExceptVertices(c, 0);
        require(sign != 0, "changed target cell has no strict single-phase certificate; target cell=" + std::to_string(c));
        const auto box = candidate.leafBounds(c);
        for (int previous : source.leafCandidates(box[0], box[1]))
            if (overlaps(box, source.leafBounds(previous)))
            {
                require(sourceSigns[previous] == sign,
                    "changed scalar reaches the interface/protected support or changes phase; target cell=" + std::to_string(c));
                if (mustPreserve[previous]) throw ProtectedPolynomialChange(previous, c);
            }
        ++result.report.certifiedChangedCells;
    }
    for (int c = 0; c < sourceMesh.cellIndexNum(); ++c)
        if (mustPreserve[c]) require(result.sourceToTargetCell[c] >= 0,
            "retained interface/support polynomial changed through shared constraints; source cell=" + std::to_string(c));
    return result;
}
}
