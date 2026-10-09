#pragma once

#include "Solvers/Operators/QuadtreeFaceGradient.h"

// Compatibility name for existing Hydro callers. Both momentum projection
// and native electric FV assembly use the SAME virtual-child mesh operator.
namespace HydroNativePressureStencil
{
using QuadtreeFaceGradient::childAtFace;
using QuadtreeFaceGradient::childCenter;
using QuadtreeFaceGradient::visitGradient;

struct FaceData
{
    double gradient, density, dualArea;
};

// centers(cell, child) borrows an already prepared, mesh-owned-generation
// stencil. All labels and 2:1 topology must be checked before a parallel pass.
// No allocation, topology mutation, or exception is introduced in evaluation.
template <class Mesh, class Field, class Centers>
FaceData faceData(const Mesh &mesh, const Edge &edge, bool horizontal, const Field &pressure,
                  const Field &density, const Centers &centers)
{
    const Cell &master = mesh.cells(edge.masterCell);
    const Cell &slave = mesh.cells(edge.slaveCell);
    double pMaster = pressure.value[master.index()], pSlave = pressure.value[slave.index()];
    double rhoMaster = density.value[master.index()], rhoSlave = density.value[slave.index()];
    const Cell &fine = master.level() >= slave.level() ? master : slave;
    if (master.level() != slave.level())
    {
        const bool coarseMaster = master.level() < slave.level();
        const Cell &coarse = coarseMaster ? master : slave;
        const auto &stencil = centers(coarse, childAtFace(coarse, fine, horizontal, coarseMaster));
        double &p = coarseMaster ? pMaster : pSlave;
        double &rho = coarseMaster ? rhoMaster : rhoSlave;
        p = rho = 0.0;
        for (std::size_t j = 0; j < stencil.ind.size(); ++j)
        {
            p += pressure.value[stencil.ind[j]] * stencil.coe[j];
            rho += density.value[stencil.ind[j]] * stencil.coe[j];
        }
    }
    return {(pMaster - pSlave) / (horizontal ? fine.size().y : fine.size().x), .5 * (rhoMaster + rhoSlave),
            fine.area()};
}

} // namespace HydroNativePressureStencil
