#include "Modules/Electricity/LegacyNodalMaxwell.h"

namespace LegacyNodalMaxwell
{
void add(const std::vector<int>& labels, const std::vector<InterfaceNode>& nodes, const Field<double>& nMaxwellStress, const Field<double>& tMaxwellStress, Field<double>& nodeFx, Field<double>& nodeFy)
{
#pragma omp parallel for num_threads(thread_num)
    for (int i = 0; i < labels.size(); i++)
    {
        const int nodeLabel = labels[i];
        const InterfaceNode& node = nodes[nodeLabel];

        nodeFx.value[i] += nMaxwellStress.value[i] * node.normal().x - tMaxwellStress.value[i] * node.normal().y;
        nodeFy.value[i] += nMaxwellStress.value[i] * node.normal().y + tMaxwellStress.value[i] * node.normal().x;
    }
}
}
