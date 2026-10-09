#pragma once

#include <vector>
#include "Core/Base/Vector2D.h"

class QuadtreeMesh;
class Edge;
template<typename T> class VertField;
namespace LevelSetHermiteGeometry { class Snapshot; }

namespace InterfaceBoundaryGeometry
{
struct Report
{
    bool axisContact = false;
    int cutCells = 0;
    // Actual Hermite zero curve + domain boundary, using the original GL8
    // divergence-theorem integral; not an independently transported VOF.
    double negativeVolume = 0.0;
};
Report inspect(const QuadtreeMesh& mesh, const VertField<double>& phi,
    bool axisymmetric, const std::vector<Vector2D>& fixedContacts);
double negativeFaceMeasure(const QuadtreeMesh& mesh, const Edge& edge,
    const LevelSetHermiteGeometry::Snapshot& scalar, bool axisymmetric);
}
