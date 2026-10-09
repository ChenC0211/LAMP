#pragma once

#include "Core/Base/Vector2D.h"
#include <array>
#include <vector>

class QuadtreeMesh;

// Read-only search index of actual rectangular leaves; owns coordinates and
// native indices, not a second mesh or a resampled field. Exact closed-box
// containment: no wall crossing by tolerance and no snapping of the query.
class CartesianLeafLocator
{
  public:
    explicit CartesianLeafLocator(const QuadtreeMesh &mesh);
    struct IncidentCells
    {
        std::array<int, 4> indices{};
        int count = 0;
    };
    IncidentCells incident(Vector2D point) const;
    // Broad-phase native leaf candidates intersecting an exact closed box.
    // The caller applies its segment/curve predicate; no field is resampled.
    std::vector<int> intersectingBox(Vector2D lower, Vector2D upper) const;

  private:
    struct Node
    {
        Vector2D lower, upper;
        int left = -1, right = -1, index = -1;
    };
    std::vector<Node> nodes_;
    int build(std::vector<Node> &leaves, int begin, int end);
    void visit(int node, Vector2D point, IncidentCells &result) const;
    void visitBox(int node, Vector2D lower, Vector2D upper, std::vector<int> &result) const;
};
