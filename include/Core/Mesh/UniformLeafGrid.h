#pragma once

#include <vector>
#include "Core/Base/Vector2D.h"

class QuadtreeMesh;

enum class CartesianGridScope { RootOnly, UniformLeaves };

// Retained logical layout of the complete, nonperiodic Cartesian leaf mesh.
// Labels, indices and refinement ancestry remain the mesh's actual identities.
// RootOnly additionally requires the original level-zero/root-grid contract;
// UniformLeaves derives spacing and dimensions from actual equal-level leaves.
// This value owns its maps, has no mesh reference, and does not track changes.
// It validates topology/metadata only: no field, material or physical measure.
class UniformLeafGrid
{
public:
    explicit UniformLeafGrid(const QuadtreeMesh& mesh,
        CartesianGridScope scope = CartesianGridScope::RootOnly);

    int nx() const noexcept { return nx_; }
    int ny() const noexcept { return ny_; }
    int level() const noexcept { return level_; }
    Vector2D spacing() const { return spacing_; }
    Vector2D minCoord() const { return lower_; }
    Vector2D maxCoord() const { return upper_; }
    // Actual labels in x-fast Cartesian order; NOT solver indices.
    const std::vector<int>& cells() const noexcept { return cells_; }
    const std::vector<int>& vertices() const noexcept { return vertices_; }
    const std::vector<int>& vertical() const noexcept { return vertical_; }
    const std::vector<int>& horizontal() const noexcept { return horizontal_; }

private:
    int nx_ = 0, ny_ = 0, level_ = 0;
    Vector2D spacing_{0.0, 0.0}, lower_{0.0, 0.0}, upper_{0.0, 0.0};
    std::vector<int> cells_, vertices_, vertical_, horizontal_;
};
