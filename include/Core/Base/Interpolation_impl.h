#pragma once

#include "Core/Base/Interpolation.h"
using namespace std;

namespace Interpolation
{

template<typename Mesh, typename Field>
double biQuadInterpolation
(
    const Mesh& mesh, 
    Field& field, 
    const Cell& cell,
    const Vector2D& coord
)
{
    if (cell.parent == -1)
    {
        return biLinearInterpolation(mesh, field, cell, coord);
    }

    double vertValue[9];

    const Cell& parentCell = mesh.cells(cell.parent);

    const Cell& sonCell0 = mesh.cells(parentCell.son[0]);
    const Cell& sonCell1 = mesh.cells(parentCell.son[1]);
    const Cell& sonCell2 = mesh.cells(parentCell.son[2]);
    const Cell& sonCell3 = mesh.cells(parentCell.son[3]);

    const Vert& cenVert = mesh.verts(sonCell0.cellToVert[2]);

    vertValue[0] = field.value[mesh.verts(sonCell0.cellToVert[0]).index()];
    vertValue[1] = field.value[mesh.verts(sonCell0.cellToVert[1]).index()];
    vertValue[2] = field.value[mesh.verts(sonCell1.cellToVert[1]).index()];
    vertValue[3] = field.value[mesh.verts(sonCell0.cellToVert[3]).index()];
    vertValue[4] = field.value[mesh.verts(sonCell0.cellToVert[2]).index()];
    vertValue[5] = field.value[mesh.verts(sonCell1.cellToVert[2]).index()];
    vertValue[6] = field.value[mesh.verts(sonCell3.cellToVert[3]).index()];
    vertValue[7] = field.value[mesh.verts(sonCell3.cellToVert[2]).index()];
    vertValue[8] = field.value[mesh.verts(sonCell2.cellToVert[2]).index()];

    double x = (coord.x - cenVert.coord().x) / cell.size().x;
    double y = (coord.y - cenVert.coord().y) / cell.size().y;

    double result = vertValue[0] * 0.25 * x * (x - 1.0) * y * (y - 1.0)
        + vertValue[1] * 0.5 * (1.0 - x * x) * y * (y - 1.0)
        + vertValue[2] * 0.25 * x * (x + 1.0) * y * (y - 1.0)
        + vertValue[3] * 0.5 * x * (x - 1.0) * (1.0 - y * y)
        + vertValue[4] * (1.0 - x * x) * (1.0 - y * y)
        + vertValue[5] * 0.5 * x * (x + 1.0) * (1.0 - y * y)
        + vertValue[6] * 0.25 * x * (x - 1.0) * y * (y + 1.0)
        + vertValue[7] * 0.5 * (1.0 - x * x) * y * (y + 1.0)
        + vertValue[8] * 0.25 * x * (x + 1.0) * y * (y + 1.0);

    return result;
}

template<typename Mesh, typename Field>
double biLinearInterpolation
(
    const Mesh& mesh, 
    Field& field, 
    const Cell& cell, 
    const Vector2D& coord
)
{
    double vertValue[4];
    for (int i = 0; i < 4; i++)
    {
        vertValue[i] = field.value[mesh.verts(cell.cellToVert[i]).index()];
    }

    double x = (coord.x - mesh.verts(cell.cellToVert[0]).coord().x) / cell.size().x;
    double y = (coord.y - mesh.verts(cell.cellToVert[0]).coord().y) / cell.size().y;

    double result = vertValue[0] * (1.0 - x) * (1.0 - y)
        + vertValue[1] * x * (1.0 - y)
        + vertValue[2] * x * y
        + vertValue[3] * (1.0 - x) * y;

    return result;
}

template<typename Mesh, typename Field>
double averageInterpolation
(
    const Mesh& mesh, 
    Field& field, 
    int preI
)
{
    double result = 0.0;
    result = field.value[preI];

    return result;
}

}; // namespace Interpolation