#pragma once

#include <cmath>
#include <initializer_list>
#include <stdexcept>
#include <vector>

#include "Core/Base/Interpolation.h"
#include "Core/Base/Structures.h"
#include "Core/Mesh/Cell.h"
#include "Core/Mesh/Edge.h"

// Shared Cartesian 2:1 virtual-child interpolation and positive-coordinate
// face gradient, extracted from the original Hydro pressure stencil. Callers
// own the active-domain closure; no inactive or solid values are supplied.
namespace QuadtreeFaceGradient
{
template <class Mesh> supCoeff childCenter(const Mesh &mesh, const Cell &cell, int child)
{
    if (child < 0 || child > 3)
        throw std::invalid_argument("native pressure child must be 0..3");
    auto first = [](std::initializer_list<int> labels) {
        for (int label : labels)
            if (label >= 0)
                return label;
        throw std::invalid_argument("native pressure child lacks a directional neighbor");
    };
    const auto &n = cell.neiCell;
    const int x = child == 0   ? first({n[0][0], n[2][0]})
                  : child == 1 ? first({n[2][0], n[0][0]})
                  : child == 2 ? first({n[2][1], n[2][0], n[0][1], n[0][0]})
                               : first({n[0][1], n[0][0], n[2][1], n[2][0]});
    const int y = child == 0   ? first({n[1][0], n[3][0]})
                  : child == 1 ? first({n[1][1], n[1][0], n[3][1], n[3][0]})
                  : child == 2 ? first({n[3][1], n[3][0], n[1][1], n[1][0]})
                               : first({n[3][0], n[1][0]});
    if (x >= mesh.cellLabelCapacity() || y >= mesh.cellLabelCapacity())
        throw std::invalid_argument("native pressure child has an invalid neighbor label");
    supCoeff result;
    result.ind = {cell.index(), mesh.cells(x).index(), mesh.cells(y).index()};
    const auto h = cell.size();
    const auto origin = cell.centerCoord();
    if (!(h.x > 0.0) || !(h.y > 0.0) || !std::isfinite(h.x) || !std::isfinite(h.y) ||
        !std::isfinite(origin.x) || !std::isfinite(origin.y))
        throw std::invalid_argument("native pressure child has invalid local coordinates");
    std::vector<Vector2D> coords(3);
    for (int i = 0; i < 3; ++i)
    {
        const int index = result.ind[i];
        if (index < 0 || index >= mesh.cellIndexNum())
            throw std::invalid_argument("native pressure child references a non-leaf cell");
        coords[i] = mesh.cells(mesh.cellIndexToLabel[index]).centerCoord();
        if (!std::isfinite(coords[i].x) || !std::isfinite(coords[i].y))
            throw std::invalid_argument("native pressure child has nonfinite coordinates");
        coords[i].x = (coords[i].x - origin.x) / h.x;
        coords[i].y = (coords[i].y - origin.y) / h.y;
    }
    const Vector2D aim((child == 1 || child == 2) ? .25 : -.25, child >= 2 ? .25 : -.25);
    result.coe = Interpolation::leastSquareInterpolationCoeff(coords, aim);
    for (double coefficient : result.coe)
        if (!std::isfinite(coefficient))
            throw std::invalid_argument("native pressure child has singular/nonfinite interpolation");
    return result;
}

inline int childAtFace(const Cell &coarse, const Cell &fine, bool horizontal, bool master)
{
    return horizontal ? (coarse.centerCoord().x < fine.centerCoord().x ? (master ? 1 : 2) : (master ? 0 : 3))
                      : (coarse.centerCoord().y < fine.centerCoord().y ? (master ? 3 : 2) : (master ? 0 : 1));
}

// Emit the same gradient used above, for FV pressure assembly. Coarse-cell
// virtual values must appear in BOTH adjacent continuity rows; replacing just
// the fine row creates a pressure operator inconsistent with corrected fluxes.
template <class Mesh, class Centers, class Visit>
void visitGradient(const Mesh &mesh, const Edge &edge, bool horizontal, const Centers &centers,
                   const Visit &visit)
{
    const Cell &master = mesh.cells(edge.masterCell);
    const Cell &slave = mesh.cells(edge.slaveCell);
    const Cell &fine = master.level() >= slave.level() ? master : slave;
    const double spacing = horizontal ? fine.size().y : fine.size().x;
    const auto side = [&](const Cell &cell, bool high) {
        const double sign = high ? 1.0 : -1.0;
        if (cell.level() == fine.level())
            visit(cell.index(), sign / spacing);
        else
        {
            const auto &stencil = centers(cell, childAtFace(cell, fine, horizontal, high));
            for (std::size_t j = 0; j < stencil.ind.size(); ++j)
                visit(stencil.ind[j], sign * stencil.coe[j] / spacing);
        }
    };
    side(master, true);
    side(slave, false);
}
} // namespace QuadtreeFaceGradient
