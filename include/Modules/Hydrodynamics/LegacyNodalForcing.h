#pragma once

#include "Core/Lagrangian/InterfaceNode.h"
#include "Core/Field/CellField.h"

// Legacy 25-cell IBM arithmetic. Inputs must have valid node indices, complete
// stencils and correctly sized fields. Spreading adds to the caller's fields.
// Keep label traversal and sequential accumulation order for compatibility.
namespace LegacyNodalForcing
{
void initializeCapillary(const std::vector<int>& labels, const std::vector<InterfaceNode>& nodes, double inverseWeber, Field<double>& nodeFx, Field<double>& nodeFy);
double fixedVelocityForce(int count, const Field<double>& nodeUx, const Field<double>& nodeUy, const Vector2D& fixedVel, double alpha, Field<double>& nodeFx, Field<double>& nodeFy);
void interpolate(const std::vector<int>& labels, const std::vector<InterfaceNode>& nodes, const CellField<double>& ux, const CellField<double>& uy, Field<double>& nodeUx, Field<double>& nodeUy);
void interpolateNoSlip(const std::vector<InterfaceNode>& nodes, int count, const CellField<double>& ux, const CellField<double>& uy, Field<double>& nodeUx, Field<double>& nodeUy);

template<typename Mesh>
void spreadSurface(Mesh& mesh, const std::vector<int>& labels, const std::vector<InterfaceNode>& nodes, const Field<double>& nodeFx, const Field<double>& nodeFy, CellField<double>& fx, CellField<double>& fy)
{
    for (int i = 0; i < labels.size(); i++)
    {
        const InterfaceNode& node = nodes[labels[i]];

        for (int j = 0; j < 25; j++)
        {
            int cellI = node.cellStencil[j];
            double D = node.cellStencilCoeff[j];

            Cell& cell = mesh.cells(mesh.cellIndexToLabel[cellI]);

            fx.value[cellI] += nodeFx.value[node.index()] * D * node.ds() / cell.area();
            fy.value[cellI] += nodeFy.value[node.index()] * D * node.ds() / cell.area();
        }
    }
}
void spreadVolume(const std::vector<InterfaceNode>& nodes, int count, const Field<double>& nodeFx, const Field<double>& nodeFy, CellField<double>& fx, CellField<double>& fy);
void spreadEvaporation(const std::vector<int>& labels, const std::vector<InterfaceNode>& nodes, const std::vector<double>& evaporationVelocity, CellField<double>& evapUx, CellField<double>& evapUy);
}
