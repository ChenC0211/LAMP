#include "Modules/Hydrodynamics/HydroFittedMomentum.h"
#include "Modules/Hydrodynamics/Detail/HydroFittedBoundaryMomentum.h"
#include "Modules/Hydrodynamics/HydroFittedDomain.h"
#include "Core/Mesh/QuadtreeMesh.h"

#include <stdexcept>

std::shared_ptr<const HydroFittedMomentumLayout> HydroFittedMomentumLayout::build(
    const QuadtreeMesh& mesh, const HydroFittedDomain& domain)
{
    domain.validateMeshIdentity(mesh);
    auto result = std::shared_ptr<HydroFittedMomentumLayout>(new HydroFittedMomentumLayout);
    result->domain_ = std::make_shared<const HydroFittedDomain>(domain);
    const int count = mesh.vEdgeIndexNum() + mesh.hEdgeIndexNum();
    result->nativeInterior_.assign(count, 0);
    if (domain.uniformLeaves()) return result;
    std::vector<unsigned char> boundary(count, 0);
    std::vector<int> visited(mesh.cellIndexNum(), -1), frontier, next;
    bool needsBoundary = false;
    for (int row : domain.nativeByVelocityRow())
    {
        const auto& face = domain.faces()[row];
        const int owner = face.negativeCell >= 0 ? face.negativeCell : face.positiveCell;
        const int level = mesh.cells(mesh.cellIndexToLabel[owner]).level();
        bool same = true;
        for (int index : {face.negativeCell, face.positiveCell}) if (index >= 0)
        {
            const auto& cell = mesh.cells(mesh.cellIndexToLabel[index]);
            same = same && cell.level() == level;
            for (const auto& side : cell.neiCell) for (int label : side)
                if (label >= 0) same = same && mesh.cells(label).level() == level;
            for (int label : cell.cellToVert) same = same && mesh.verts(label).supCell < 0;
        }
        if (same) continue;
        bool safe = face.negativeCell >= 0 && face.positiveCell >= 0;
        frontier.clear();
        for (int index : {face.negativeCell, face.positiveCell}) if (index >= 0)
        {
            frontier.push_back(index); visited[index] = row;
        }
        // Retain the old native virtual-cell read halo verbatim. An unsafe
        // row is routed to physical boundary segments, never fed to old RZ.
        for (int depth = 0; depth < 3 && safe; ++depth)
        {
            next.clear();
            for (int index : frontier)
                for (int vertex : mesh.cells(mesh.cellIndexToLabel[index]).cellToVert)
                    for (int label : mesh.verts(vertex).vertToCell)
                    {
                        if (label < 0) { safe = false; continue; }
                        const int neighbor = mesh.cells(label).index();
                        if (domain.cellRowByIndex()[neighbor] < 0) { safe = false; continue; }
                        if (visited[neighbor] != row)
                        {
                            visited[neighbor] = row; next.push_back(neighbor);
                        }
                    }
            frontier.swap(next);
        }
        if (safe) result->nativeInterior_[row] = 1;
        else { boundary[row] = 1; needsBoundary = true; }
    }
    if (needsBoundary)
        result->boundary_ = HydroFittedDetail::BoundaryMomentum::build(
            mesh, domain, result->nativeInterior_, std::move(boundary));
    return result;
}

void HydroFittedMomentumLayout::validate(const QuadtreeMesh& mesh, const HydroFittedDomain& domain) const
{
    if (!domain_) throw std::invalid_argument("unprepared fitted momentum layout");
    domain_->validateMeshIdentity(mesh);
    domain.validateMeshIdentity(mesh);
    if (domain_->axisymmetric() != domain.axisymmetric()
        || domain_->cellRowByIndex() != domain.cellRowByIndex()
        || domain_->velocityRowByNative() != domain.velocityRowByNative()
        || domain_->faces().size() != domain.faces().size())
        throw std::invalid_argument("fitted momentum layout has a different fluid domain");
    for (std::size_t i = 0; i < domain.faces().size(); ++i)
    {
        const auto& a = domain_->faces()[i]; const auto& b = domain.faces()[i];
        if (a.negativeCell != b.negativeCell || a.positiveCell != b.positiveCell
            || a.dualVolume != b.dualVolume || a.area != b.area
            || a.normalSpacing != b.normalSpacing || a.tangentSpacing != b.tangentSpacing
            || a.boundary.type != b.boundary.type || a.boundary.id != b.boundary.id
            || a.boundary.pressure != b.boundary.pressure
            || a.boundary.incomingDensity != b.boundary.incomingDensity
            || a.boundary.incomingViscosity != b.boundary.incomingViscosity)
            throw std::invalid_argument("fitted momentum layout has stale geometry or boundary data");
    }
}
