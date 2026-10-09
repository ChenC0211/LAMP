#include "Modules/Interface/Compatibility/InterfaceBoundaryGeometry.h"

#include "Core/Geometry/LevelSetInterfacePanels.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <cmath>
#include <stdexcept>

namespace InterfaceBoundaryGeometry
{
double negativeFaceMeasure(const QuadtreeMesh& mesh, const Edge& edge,
    const LevelSetHermiteGeometry::Snapshot& scalar, bool rz)
{
    const auto a = mesh.verts(edge.edgeToVert[0]).coord();
    const auto b = mesh.verts(edge.edgeToVert[1]).coord();
    double result = 0.0;
    const double pi = std::acos(-1.0);
    for (const auto& piece : scalar.partitionSegment(a, b))
        if (piece.negative)
            result += (piece.end-piece.begin).mod() *
                (rz ? pi*(piece.begin.x+piece.end.x) : 1.0);
    if (!std::isfinite(result) || result < 0.0)
        throw std::runtime_error("invalid negative-phase boundary measure");
    return result;
}

Report inspect(const QuadtreeMesh& mesh, const VertField<double>& phi,
    bool rz, const std::vector<Vector2D>& fixedContacts)
{
    LevelSetHermiteGeometry::Snapshot scalar(mesh, phi, rz);
    LevelSetInterfacePanels panels(mesh, scalar, rz, fixedContacts);
    Report result;
    result.axisContact = panels.beginEndpointKind() == LevelSetInterfacePanels::EndpointKind::Axis
        || panels.endEndpointKind() == LevelSetInterfacePanels::EndpointKind::Axis;
    result.cutCells = static_cast<int>(panels.panels().size());
    if (panels.closed() || result.cutCells == 0)
        throw std::runtime_error("prescribedScalar requires a resolved declared-boundary interface chain");
    // div(0,y-ymin)=1 (also under the RZ divergence). Other domain sides
    // have zero flux for this auxiliary vector field, except the north side.
    long double volume = 0.0L;
    for (const auto& panel : panels.panels())
        for (const auto& node : panel.quadratureNodes)
            volume += static_cast<long double>(node.physicalWeight)
                * (node.frame.point.y-mesh.minCoord.y) * node.frame.normal.y;
    for (int i = 0; i < mesh.hEdgeIndexNum(); ++i)
    {
        const auto& edge = mesh.hEdges(mesh.hEdgeIndexToLabel[i]);
        if (edge.boundary() == Boundary::north)
            volume += (mesh.maxCoord.y-mesh.minCoord.y)*negativeFaceMeasure(mesh, edge, scalar, rz);
    }
    const double domainVolume = (mesh.maxCoord.y-mesh.minCoord.y) * (rz
        ? std::acos(-1.0)*(mesh.maxCoord.x*mesh.maxCoord.x-mesh.minCoord.x*mesh.minCoord.x)
        : mesh.maxCoord.x-mesh.minCoord.x);
    if (!std::isfinite(volume) || !(volume > 0.0L) || !(volume < domainVolume))
        throw std::runtime_error("invalid geometric negative-phase volume for prescribedScalar interface");
    result.negativeVolume = static_cast<double>(volume);
    return result;
}
}
