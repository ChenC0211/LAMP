#pragma once

#include "Core/Base/Vector2D.h"

class QuadtreeMesh;

// Validate the existing complete Cartesian quadtree, including its original
// labels/ancestry, real split faces and hanging-node identities. No resampling,
// geometry mutation, field interpolation or boundary-condition inference.
// The returned spacing is a minimum, never a global mixed-mesh stencil width.
struct CartesianLeafMetrics
{
    Vector2D minimumSpacing{0.0, 0.0};
    bool mixedLevels = false;
};

CartesianLeafMetrics validateCartesianLeafMesh(const QuadtreeMesh &mesh);
