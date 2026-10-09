#pragma once

#include "Core/Mesh/QuadtreeMesh.h"
#include "Modules/Hydrodynamics/HydroFittedDomain.h"

#include <cmath>
#include <stdexcept>

namespace HydroFittedDetail
{
// Geometry for the existing tau_nt=0 tangential trace. The two samples are
// actual pressure-normal boundary faces, ordered along the tangential row's
// normal direction. No new boundary data, ghost velocity, or time layer.
struct PressureTraceGeometry
{
    int lowNative = -1, highNative = -1;
    double derivativeSpacing = 0.0;
    double traceSpacing = 0.0; // Twice the distance to the pressure plane.
};

inline PressureTraceGeometry pressureTraceGeometry(const QuadtreeMesh& mesh,
    const HydroFittedDomain& domain, int native, int tangentSide)
{
    using Type = HydroFittedDomain::BoundaryType;
    const auto require = [](bool valid, const char* reason) {
        if (!valid) throw std::invalid_argument(reason);
    };
    require(native >= 0 && native < static_cast<int>(domain.faces().size())
        && (tangentSide == 0 || tangentSide == 1), "invalid fitted pressure trace address");
    const auto& f = domain.faces()[native];
    require(f.negativeCell >= 0 && f.positiveCell >= 0,
        "two-owner pressure trace requires an interior tangential face");
    const auto& lowCell = mesh.cells(mesh.cellIndexToLabel[f.negativeCell]);
    const auto& highCell = mesh.cells(mesh.cellIndexToLabel[f.positiveCell]);
    const auto boundaryRow = [&](const Cell& cell) {
        const auto* labels = f.horizontal ? cell.cellToVEdge[tangentSide] : cell.cellToHEdge[tangentSide];
        require(labels[0] >= 0 && labels[1] < 0,
            "pressure trace requires a complete physical owner side");
        return f.horizontal ? mesh.vEdges(labels[0]).index()
            : mesh.vEdgeIndexNum() + mesh.hEdges(labels[0]).index();
    };
    PressureTraceGeometry result{boundaryRow(lowCell), boundaryRow(highCell),
        f.normalSpacing, f.tangentSpacing};
    const auto& a = domain.faces()[result.lowNative];
    const auto& b = domain.faces()[result.highNative];
    const auto pressure = [](Type type) {
        return type == Type::PressureReservoir || type == Type::PressureOutlet;
    };
    const double plane = f.horizontal ? a.center.x : a.center.y;
    require(pressure(a.boundary.type) && pressure(b.boundary.type)
        && a.outwardSign == b.outwardSign && a.outwardSign == (tangentSide ? 1.0 : -1.0)
        && plane == (f.horizontal ? b.center.x : b.center.y),
        "zero-shear samples must belong to the same physical pressure plane");
    if (lowCell.level() != highCell.level())
    {
        // A 2:1 tangential row's normalSpacing is the VIRTUAL fine dual width.
        // These boundary samples are at real coarse/fine owner centers. Their
        // separation is 3h/2, not h; using h amplifies the cross derivative.
        result.derivativeSpacing = f.horizontal ? b.center.y - a.center.y : b.center.x - a.center.x;
        result.traceSpacing = 2 * a.outwardSign * (plane - (f.horizontal ? f.center.x : f.center.y));
    }
    require(std::isfinite(result.derivativeSpacing) && result.derivativeSpacing > 0
        && std::isfinite(result.traceSpacing) && result.traceSpacing > 0,
        "unresolved fitted pressure trace distances");
    return result; // Same-level rows retain their original nominal arithmetic.
}
} // namespace HydroFittedDetail
