#pragma once

#include "Modules/Hydrodynamics/HydroNativeRZRow_impl.h"
#include "Core/Base/Base.h"
#include "Core/Field/CellField.h"
#include "Core/Field/VertField.h"
#include "Core/Mesh/Cell.h"
#include "Core/Mesh/Edge.h"
#include "Core/Mesh/Vert.h"

// Common physical data for an already selected regular native RZ stencil.
// Domain strategies keep their own active-halo/boundary/coarse-fine gates and
// index maps. In particular, this function never samples a rejected stencil.
// Inline arithmetic keeps the original caller's floating-point boundary.
namespace HydroNativeRZRowDetail
{
template<class VertexViscosityAt>
void prepareRegularData(HydroNativeRZRow::Row& row, bool radial,
    const Edge& edge, const Cell& master, const Cell& slave,
    const Vert& start, const Vert& end,
    const CellField<double>& density, const CellField<double>& viscosity,
    double shearLowDistance, double shearHighDistance, VertexViscosityAt vertexViscosityAt)
{
    auto& g = row.metrics;
    g.volume = 2.0 * Pi * edge.midpointCoord().x * master.area();
    g.dx = master.size().x;
    g.dy = master.size().y;
    g.areaNS = 2.0 * Pi * edge.midpointCoord().x * master.size().x;
    g.areaEast = 2.0 * Pi * (radial ? master.centerCoord().x : end.coord().x) * master.size().y;
    g.areaWest = 2.0 * Pi * (radial ? slave.centerCoord().x : start.coord().x) * master.size().y;
    // Preserve each strategy's actual native midpoint source and subtraction.
    g.shearLowDistance = shearLowDistance;
    g.shearHighDistance = shearHighDistance;
    auto& m = row.materials;
    m.rho = 0.5 * (density.value[master.index()] + density.value[slave.index()]);
    if (radial)
    {
        m.muWest = viscosity.value[slave.index()];
        m.muEast = viscosity.value[master.index()];
        m.muSouth = vertexViscosityAt(start.index());
        m.muNorth = vertexViscosityAt(end.index());
        m.rhoWest = density.value[slave.index()];
        m.rhoEast = density.value[master.index()];
    }
    else
    {
        m.muSouth = viscosity.value[slave.index()];
        m.muNorth = viscosity.value[master.index()];
        m.muWest = vertexViscosityAt(start.index());
        m.muEast = vertexViscosityAt(end.index());
        m.rhoSouth = density.value[slave.index()];
        m.rhoNorth = density.value[master.index()];
    }
}

// Only convection consumes transverse density. Native matrix-only assembly
// does not read it; fitted callers retain their checked active-fluid sampler.
template<class VertexDensityAt>
void prepareTransverseDensity(HydroNativeRZRow::Materials& material, bool radial,
    int start, int end, VertexDensityAt vertexDensityAt)
{
    if (radial)
    {
        material.rhoSouth = vertexDensityAt(start);
        material.rhoNorth = vertexDensityAt(end);
    }
    else
    {
        material.rhoWest = vertexDensityAt(start);
        material.rhoEast = vertexDensityAt(end);
    }
}
}

namespace HydroNativeRZRowDetail
{
// Traverse the low-coordinate side of a validated quadtree cell. A historical
// child can itself be refined: its old edge index is then NOT a velocity DOF.
// Preserve the old equal-half side average, descending only where necessary.
// Callers validate mesh topology before parallel assembly; this traversal
// neither reads field values nor throws across an OpenMP region.
template<typename Mesh, typename Emit>
void visitLeafLowSide(const Mesh& mesh, const Cell& cell, bool radial, double coefficient, const Emit& emit)
{
    if (cell.son[0] >= 0)
    {
        visitLeafLowSide(mesh, mesh.cells(cell.son[0]), radial, 0.5 * coefficient, emit);
        visitLeafLowSide(mesh, mesh.cells(cell.son[radial ? 3 : 1]), radial, 0.5 * coefficient, emit);
        return;
    }
    const int* labels = radial ? cell.cellToVEdge[0] : cell.cellToHEdge[0];
    const double weight = labels[1] >= 0 ? 0.5 * coefficient : coefficient;
    for (int label : {labels[0], labels[1]})
        if (label >= 0)
            emit((radial ? 0 : mesh.vEdgeIndexNum())
                + (radial ? mesh.vEdges(label).index() : mesh.hEdges(label).index()), weight);
}

// This is only a mesh-to-native-data adapter. The shared row arithmetic has
// no mesh/level assumptions. A conservative complete 1-ring gate leaves every
// old boundary, hanging-node and coarse/fine branch untouched.
template<typename Mesh>
bool regularRow(Mesh& mesh, const Edge& edge, bool radial,
    const CellField<double>& density, const CellField<double>& viscosity,
    const VertField<double>& vertexViscosity, HydroNativeRZRow::Row& row)
{
    if (edge.masterCell < 0 || edge.slaveCell < 0) return false;
    const int level = mesh.cells(edge.masterCell).level();
    const auto regularCell = [&](int label) {
        if (label < 0) return false;
        const auto& cell = mesh.cells(label);
        if (cell.level() != level || cell.supCellIndex != -1) return false;
        for (int side = 0; side < 2; ++side)
            if (cell.cellToVEdge[side][0] < 0 || cell.cellToHEdge[side][0] < 0
                || cell.cellToVEdge[side][1] != -1 || cell.cellToHEdge[side][1] != -1)
                return false;
        return true;
    };
    const auto regularEdge = [&](const Edge& candidate) {
        return candidate.index() >= 0 && regularCell(candidate.masterCell)
            && regularCell(candidate.slaveCell);
    };
    if (!regularEdge(edge)) return false;
    const auto& master = mesh.cells(edge.masterCell);
    const auto& slave = mesh.cells(edge.slaveCell);
    const auto& start = mesh.verts(edge.edgeToVert[0]);
    const auto& end = mesh.verts(edge.edgeToVert[1]);
    const int lowLabel = radial ? start.vertToVEdge[0] : start.vertToHEdge[0];
    const int highLabel = radial ? end.vertToVEdge[1] : end.vertToHEdge[1];
    if (start.supCell != -1 || end.supCell != -1 || lowLabel < 0 || highLabel < 0)
        return false;
    const auto& low = radial ? mesh.vEdges(lowLabel) : mesh.hEdges(lowLabel);
    const auto& high = radial ? mesh.vEdges(highLabel) : mesh.hEdges(highLabel);
    const auto& normalLow = radial ? mesh.vEdges(slave.cellToVEdge[0][0])
        : mesh.hEdges(slave.cellToHEdge[0][0]);
    const auto& normalHigh = radial ? mesh.vEdges(master.cellToVEdge[1][0])
        : mesh.hEdges(master.cellToHEdge[1][0]);
    const auto& crossSW = radial ? mesh.hEdges(slave.cellToHEdge[0][0])
        : mesh.vEdges(slave.cellToVEdge[0][0]);
    const auto& crossSE = radial ? mesh.hEdges(master.cellToHEdge[0][0])
        : mesh.vEdges(slave.cellToVEdge[1][0]);
    const auto& crossNW = radial ? mesh.hEdges(slave.cellToHEdge[1][0])
        : mesh.vEdges(master.cellToVEdge[0][0]);
    const auto& crossNE = radial ? mesh.hEdges(master.cellToHEdge[1][0])
        : mesh.vEdges(master.cellToVEdge[1][0]);
    for (const Edge* candidate : {&low, &high, &normalLow, &normalHigh,
        &crossSW, &crossSE, &crossNW, &crossNE})
        if (!regularEdge(*candidate)) return false;

    const int normalOffset = radial ? 0 : mesh.vEdgeIndexNum();
    const int crossOffset = radial ? mesh.vEdgeIndexNum() : 0;
    auto& i = row.indices;
    i.center = normalOffset + edge.index();
    i.west = normalOffset + (radial ? normalLow.index() : low.index());
    i.east = normalOffset + (radial ? normalHigh.index() : high.index());
    i.south = normalOffset + (radial ? low.index() : normalLow.index());
    i.north = normalOffset + (radial ? high.index() : normalHigh.index());
    i.crossSW = crossOffset + crossSW.index();
    i.crossSE = crossOffset + crossSE.index();
    i.crossNW = crossOffset + crossNW.index();
    i.crossNE = crossOffset + crossNE.index();

    prepareRegularData(row, radial, edge, master, slave, start, end, density, viscosity,
        radial ? edge.midpointCoord().y - low.midpointCoord().y : edge.midpointCoord().x - low.midpointCoord().x,
        radial ? high.midpointCoord().y - edge.midpointCoord().y : high.midpointCoord().x - edge.midpointCoord().x,
        [&](int index) { return vertexViscosity.value[index]; });
    return true;
}
}
