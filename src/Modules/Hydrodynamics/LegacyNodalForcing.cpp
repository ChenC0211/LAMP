#include "Modules/Hydrodynamics/LegacyNodalForcing.h"

namespace LegacyNodalForcing
{
void initializeCapillary(const std::vector<int>& labels, const std::vector<InterfaceNode>& nodes, double inverseWeber, Field<double>& nodeFx, Field<double>& nodeFy)
{
    const int nNodes = static_cast<int>(labels.size());
    nodeFx.initialize(nNodes, 0.0);
    nodeFy.initialize(nNodes, 0.0);
#pragma omp parallel for num_threads(thread_num)
    for (int i = 0; i < nNodes; i++)
    {
        const int nodeLabel = labels[i];
        const InterfaceNode& node = nodes[nodeLabel];

        nodeFx.value[i] = node.curvature() * inverseWeber * node.normal().x;
        nodeFy.value[i] = node.curvature() * inverseWeber * node.normal().y;
    }
}

double fixedVelocityForce(int count, const Field<double>& nodeUx, const Field<double>& nodeUy, const Vector2D& fixedVel, double alpha, Field<double>& nodeFx, Field<double>& nodeFy)
{
    nodeFx.initialize(count, 0.0);
    nodeFy.initialize(count, 0.0);
    double _residual2 = 0.0;
#pragma omp parallel for num_threads(thread_num) reduction(max:_residual2)
    for (int i = 0; i < count; i++)
    {
        double dUx = fixedVel.x - nodeUx.value[i];
        double dUy = fixedVel.y - nodeUy.value[i];
        nodeFx.value[i] = alpha * dUx;
        nodeFy.value[i] = alpha * dUy;

        _residual2 = std::max(_residual2, dUx * dUx + dUy * dUy);
    }

    return sqrt(_residual2);
}

void interpolate(const std::vector<int>& labels, const std::vector<InterfaceNode>& nodes, const CellField<double>& ux, const CellField<double>& uy, Field<double>& nodeUx, Field<double>& nodeUy)
{
    nodeUx.initialize(labels.size(), 0.0);
    nodeUy.initialize(labels.size(), 0.0);
#pragma omp parallel for num_threads(thread_num)
    for (int i = 0; i < labels.size(); i++)
    {
        const InterfaceNode& node = nodes[labels[i]];

        double Ux_sum = 0.0;
        double Uy_sum = 0.0;
        for (int j = 0; j < 25; j++)
        {
            int cellI = node.cellStencil[j];
            double D = node.cellStencilCoeff[j];
            
            Ux_sum += ux.value[cellI] * D;
            Uy_sum += uy.value[cellI] * D;
        }

        nodeUx.value[i] = Ux_sum;
        nodeUy.value[i] = Uy_sum;
    }
}

void interpolateNoSlip(const std::vector<InterfaceNode>& nodes, int count, const CellField<double>& ux, const CellField<double>& uy, Field<double>& nodeUx, Field<double>& nodeUy)
{
    nodeUx.initialize(count, 0.0);
    nodeUy.initialize(count, 0.0);
#pragma omp parallel for num_threads(thread_num)
    for (int i = 0; i < count; i++)
    {
        const InterfaceNode& node = nodes[i];

        double Ux_sum = 0.0;
        double Uy_sum = 0.0;
        for (int j = 0; j < 25; j++)
        {
            int cellI = node.cellStencil[j];
            double D = node.cellStencilCoeff[j];
            
            Ux_sum += ux.value[cellI] * D;
            Uy_sum += uy.value[cellI] * D;
        }

        nodeUx.value[i] = Ux_sum;
        nodeUy.value[i] = Uy_sum;
    }
}

void spreadVolume(const std::vector<InterfaceNode>& nodes, int count, const Field<double>& nodeFx, const Field<double>& nodeFy, CellField<double>& fx, CellField<double>& fy)
{
    for (int i = 0; i < count; i++)
    {
        const InterfaceNode& node = nodes[i];

        for (int j = 0; j < 25; j++)
        {
            int cellI = node.cellStencil[j];
            double D = node.cellStencilCoeff[j];

            
            fx.value[cellI] += nodeFx.value[node.index()] * D;
            fy.value[cellI] += nodeFy.value[node.index()] * D;
        }
    }
}

void spreadEvaporation(const std::vector<int>& labels, const std::vector<InterfaceNode>& nodes, const std::vector<double>& evaporationVelocity, CellField<double>& evapUx, CellField<double>& evapUy)
{
    for (int i = 0; i < labels.size(); i++)
    {
        const InterfaceNode& node = nodes[labels[i]];

        double _evapVx = -evaporationVelocity[labels[i]] * node.normal().x;
        double _evapVy = -evaporationVelocity[labels[i]] * node.normal().y;

        for (int j = 0; j < 25; j++)
        {
            int cellI = node.cellStencil[j];
            double D = node.cellStencilCoeff[j];
            
            evapUx.value[cellI] += _evapVx * D;
            evapUy.value[cellI] += _evapVy * D;
        }
    }
}

}
