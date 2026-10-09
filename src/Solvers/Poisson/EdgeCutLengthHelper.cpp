#include "Solvers/Poisson/EdgeCutLengthHelper.h"

#include <algorithm>
#include <cmath>

#include "Core/Mesh/QuadtreeMesh.h"

using std::max;
using std::min;

namespace Solvers
{
namespace Poisson
{

static inline double clamp01(double x)
{
    return std::max(0.0, std::min(1.0, x));
}

static void computeForEdge(const QuadtreeMesh& mesh, Edge& e, const VertField<double>& psi, int levelSetI, double clampEps)
{
    const int vLbl0 = e.edgeToVert[0];
    const int vLbl1 = e.edgeToVert[1];
    const Vert& V0 = mesh.verts(vLbl0);
    const Vert& V1 = mesh.verts(vLbl1);
    const double p0 = psi.value[V0.index()];
    const double p1 = psi.value[V1.index()];
    const double len = e.length();

    const bool isCut = (p0 * p1 <= 0.0) && (len > 0.0);
    e.isCut[levelSetI] = isCut;
    e.setCut(isCut);

    if (!isCut)
    {
        e.setCutLength(0.0);
        return;
    }

    // Linear interpolation parameter from v0 to v1.
    double t = 0.5;
    const double denom = (p0 - p1);
    if (std::abs(denom) > 1e-30)
        t = p0 / denom;

    t = clamp01(t);
    double cutLen = len * t;

    // Avoid exactly 0 or len to prevent downstream division / degeneracy.
    const double eps = std::max(0.0, std::min(clampEps, 0.49 * len));
    cutLen = std::max(eps, std::min(len - eps, cutLen));

    e.setCutLength(cutLen);
}

void EdgeCutLengthHelper::computeLinear(QuadtreeMesh& mesh, const VertField<double>& psi, int levelSetI, double clampEps)
{
    // Horizontal edges
    for (int i = 0; i < mesh.hEdgeIndexNum(); ++i)
    {
        const int lbl = mesh.hEdgeIndexToLabel[i];
        Edge& e = mesh.hEdges(lbl);
        computeForEdge(mesh, e, psi, levelSetI, clampEps);
    }

    // Vertical edges
    for (int i = 0; i < mesh.vEdgeIndexNum(); ++i)
    {
        const int lbl = mesh.vEdgeIndexToLabel[i];
        Edge& e = mesh.vEdges(lbl);
        computeForEdge(mesh, e, psi, levelSetI, clampEps);
    }
}

} // namespace Poisson
} // namespace Solvers
