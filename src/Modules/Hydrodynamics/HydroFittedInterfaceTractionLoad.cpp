#include "Modules/Hydrodynamics/HydroFittedInterfaceTractionLoad.h"
#include "Modules/Hydrodynamics/HydroFittedDomain.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/CartesianMeasure.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace Quantum::Hydrodynamics
{
namespace {
void require(bool ok, const char* message)
{
    if (!ok) throw std::invalid_argument(std::string("[HydroFittedInterfaceTractionLoad] ") + message);
}
double checked(long double value, const char* message)
{
    const double result = static_cast<double>(value);
    require(std::isfinite(value) && std::isfinite(result), message);
    return result;
}
}

HydroFittedInterfaceTractionLoad::HydroFittedInterfaceTractionLoad(const QuadtreeMesh& mesh)
    : forceX(mesh, 0.0), forceY(mesh, 0.0) {}

HydroFittedInterfaceTractionLoad buildFittedInterfaceTractionLoad(const QuadtreeMesh& mesh,
    const HydroFittedDomain& domain, const HydroFittedInterfaceLoadGeometry& geometry,
    const std::shared_ptr<const InterfaceTraceGeometry>& traceGeometry,
    const std::vector<Vector2D>& traction)
{
    geometry.validate(mesh, domain, traceGeometry);
    const auto& faces = geometry.faces();
    require(traction.size() == faces.size(), "traction count must match the complete retained load support");
    const bool native = traceGeometry->source()->usesNativeLeaves();
    const auto h = native ? Vector2D{} : domain.spacing();
    const auto metricMatches = [](double supplied, double actual) {
        return std::isfinite(supplied) && supplied > 0.0 && std::isfinite(actual) && actual > 0.0
            && std::abs(supplied - actual) <= 256.0 * std::numeric_limits<double>::epsilon() * std::max(supplied, actual);
    };
    HydroFittedInterfaceTractionLoad result(mesh);
    long double totalX = 0.0L, totalY = 0.0L;
    for (std::size_t p = 0; p < faces.size(); ++p)
    {
        const auto& face = faces[p];
        const int row = face.edgeIndex + (face.horizontal ? mesh.vEdgeIndexNum() : 0);
        require(row >= 0 && row < static_cast<int>(domain.faces().size()), "load face outside native storage");
        const auto& receiving = domain.faces()[row];
        require(receiving.velocityRow >= 0 && receiving.horizontal == face.horizontal
            && receiving.negativeCell == face.lowCellIndex && receiving.positiveCell == face.highCellIndex
            && face.lowCellIndex >= 0 && face.highCellIndex >= 0
            && receiving.center.x == face.point.x && receiving.center.y == face.point.y,
            "load face differs from its actual free active-active receiving row");
        require(std::isfinite(face.gradient) && std::isfinite(face.transverseGradient)
            && (face.gradient != 0.0 || face.transverseGradient != 0.0), "invalid retained gradient union");
        const auto& edge = face.horizontal ? mesh.hEdges(receiving.edgeLabel) : mesh.vEdges(receiving.edgeLabel);
        const auto& lowCell = mesh.cells(edge.slaveCell);
        const auto& highCell = mesh.cells(edge.masterCell);
        auto low = lowCell.centerCoord(), high = highCell.centerCoord();
        if (native && lowCell.level() != highCell.level())
        {
            const auto& fine = lowCell.level() > highCell.level() ? lowCell : highCell;
            const double width = face.horizontal ? fine.size().y : fine.size().x;
            if (lowCell.level() < highCell.level()) {
                low = face.point;
                if (face.horizontal) low.y -= .5 * width;
                else low.x -= .5 * width;
            } else {
                high = face.point;
                if (face.horizontal) high.y += .5 * width;
                else high.x += .5 * width;
            }
        }
        const double distance = face.horizontal ? high.y - low.y : high.x - low.x;
        const double area = CartesianMeasure::face(edge.length(), face.point.x, domain.axisymmetric());
        const double dualMeasure = domain.axisymmetric() && !face.horizontal
            ? 0.5 * CartesianMeasure::twoPi * (high.x + low.x) * distance * edge.length()
            : area * distance;
        require(metricMatches(face.dualArea, native ? distance * edge.length() : h.x * h.y)
            && metricMatches(face.dualMeasure, dualMeasure), "load dual measures disagree with actual native geometry");
        const double tau = traction[p].x, tangent = traction[p].y;
        require(std::isfinite(tau) && std::isfinite(tangent), "non-finite physical traction");
        const double normalLoad = checked(static_cast<long double>(tau) * face.gradient,
            "non-finite normal face load");
        const double tangentLoad = checked(static_cast<long double>(tangent) * face.transverseGradient,
            "non-finite transverse face load");
        require(tau == 0.0 || face.gradient == 0.0 || normalLoad != 0.0, "nonzero normal face load underflow");
        require(tangent == 0.0 || face.transverseGradient == 0.0 || tangentLoad != 0.0,
            "nonzero transverse face load underflow");
        const double force = checked(static_cast<long double>(normalLoad) + (face.horizontal ? tangentLoad : -tangentLoad),
            "non-finite full interface face load");
        const double integral = checked(static_cast<long double>(force) * face.dualMeasure, "dual-face integral overflow");
        require(force == 0.0 || integral != 0.0, "nonzero dual-face integral underflow");
        if (face.horizontal) result.forceY.value[face.edgeIndex] = force;
        else result.forceX.value[face.edgeIndex] = force;
        result.maxForce = std::max(result.maxForce, std::abs(force));
        (face.horizontal ? totalY : totalX) += static_cast<long double>(force) * face.dualMeasure;
    }
    result.totalForceX = checked(totalX, "x-force integral overflow");
    result.totalForceY = checked(totalY, "y-force integral overflow");
    require((totalX == 0.0L || result.totalForceX != 0.0) && (totalY == 0.0L || result.totalForceY != 0.0),
        "nonzero total-force integral underflow");
    return result;
}
}
