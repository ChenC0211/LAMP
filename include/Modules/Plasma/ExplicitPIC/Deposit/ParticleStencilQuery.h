#pragma once

#include <cstddef>
#include <vector>

#include "Core/Base/Vector2D.h"

class BackgroundMesh;
class Cell;

// ----------------------------------------------------------------------------
// ParticleStencilQuery
//
// Enumerate the set of leaf cells whose axis-aligned bounding box (AABB)
//   Enumerate the set of leaf cells whose axis-aligned bounding box (AABB)
//   intersects a particle's CIC support box: [xp-hp, xp+hp] x [yp-hp, yp+hp].
//
// This utility is intentionally independent of PIC species/particles; it only
// needs access to mesh topology (leaf-cell neighbor relations) and geometry
// (cell center + size).
// ----------------------------------------------------------------------------

class ParticleStencilQuery
{
public:
    struct AABB
    {
        double xmin = 0.0;
        double xmax = 0.0;
        double ymin = 0.0;
        double ymax = 0.0;
    };

    static AABB supportAABB(double xp, double yp, double hp);

    // Returns true if the cell's rectangle overlaps the AABB with positive area.
    static bool cellIntersectsAABB(const Cell& cell, const AABB& box);

    // Enumerate intersecting leaf cells using a BFS from hostCell.
    // If hostCell < 0, it will be computed using mesh.findCellLabel({xp,yp}).
    //
    // debugPrint: when true, prints one-line stencil summary:
    //   bbox and number of leaf cells returned.
    // Overload that writes into a caller-provided vector to avoid per-call allocations.
    // The output vector will be cleared and filled with leaf-cell labels.
    static void enumerateLeafCellsIntersectParticleSupport(
        const BackgroundMesh& mesh,
        double xp,
        double yp,
        double hp,
        std::vector<int>* stencilOut,
        bool debugPrint = false,
        int hostCell = -1,
        bool deterministic = true);

    static std::vector<int> enumerateLeafCellsIntersectParticleSupport(
        const BackgroundMesh& mesh,
        double xp,
        double yp,
        double hp,
        bool debugPrint = false,
        int hostCell = -1,
        bool deterministic = true);

private:
    static bool isLeafCell(const Cell& cell);
};
