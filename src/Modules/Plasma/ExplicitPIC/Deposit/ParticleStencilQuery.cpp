#include "Modules/Plasma/ExplicitPIC/Deposit/ParticleStencilQuery.h"

#include <iostream>
#include <algorithm>
#include <vector>

#include "Core/IO/ConsoleControl.h"

#include "Core/Mesh/BackgroundMesh.h"
#include "Core/Mesh/Cell.h"

ParticleStencilQuery::AABB ParticleStencilQuery::supportAABB(
    const double xp, const double yp, const double hp)
{
    AABB box;
    box.xmin = xp - hp;
    box.xmax = xp + hp;
    box.ymin = yp - hp;
    box.ymax = yp + hp;
    return box;
}

bool ParticleStencilQuery::isLeafCell(const Cell& cell)
{
    if (cell.label() == -1)
        return false;
    return (cell.son[0] == -1 && cell.son[1] == -1 && cell.son[2] == -1 && cell.son[3] == -1);
}

bool ParticleStencilQuery::cellIntersectsAABB(const Cell& cell, const AABB& box)
{
    const Vector2D c = cell.centerCoord();
    const Vector2D s = cell.size();
    const double xmin = c.x - 0.5 * s.x;
    const double xmax = c.x + 0.5 * s.x;
    const double ymin = c.y - 0.5 * s.y;
    const double ymax = c.y + 0.5 * s.y;

    // Positive-area overlap test.
    if (xmax <= box.xmin || xmin >= box.xmax)
        return false;
    if (ymax <= box.ymin || ymin >= box.ymax)
        return false;
    return true;
}

void ParticleStencilQuery::enumerateLeafCellsIntersectParticleSupport(
    const BackgroundMesh& mesh,
    const double xp,
    const double yp,
    const double hp,
    std::vector<int>* stencilOut,
    const bool debugPrint,
    const int hostCellIn,
    const bool deterministic)
{
    if (!stencilOut)
        return;

    const AABB box = supportAABB(xp, yp, hp);

    const int hostCell = (hostCellIn >= 0) ? hostCellIn : mesh.findCellLabel(Vector2D(xp, yp));

    // NOTE:
    // This routine is called once per particle during CouplingWeightsCache::build.
    // Avoid per-call dynamic allocations (unordered_set/queue/function recursion).
    // The support stencil is typically small (O(1) leaf cells), so a linear
    // "visited" set is fast and allocation-free.
    struct Scratch
    {
        std::vector<int> visited;
        std::vector<int> queue;     // BFS frontier of *leaf* cells
        std::vector<int> descend;   // stack for descending non-leaf cells
    };

    thread_local Scratch scratch;
    scratch.visited.clear();
    scratch.queue.clear();
    scratch.descend.clear();

    // Conservative reserves: keep capacities across calls (thread-local).
    scratch.visited.reserve(128);
    scratch.queue.reserve(64);
    scratch.descend.reserve(64);

    // Output is written into the caller-provided vector to avoid per-call allocations.
    std::vector<int>& stencil = *stencilOut;
    stencil.clear();
    stencil.reserve(32);

    auto visitedContains = [&](const int cellLabel) -> bool
    {
        return std::find(scratch.visited.begin(), scratch.visited.end(), cellLabel)
            != scratch.visited.end();
    };

    // Push a cell label, descending to leaf children if necessary.
    // Matches legacy behavior:
    //   - mark visited before intersection test
    //   - only enqueue *leaf* cells that intersect
    auto pushCellOrDescend = [&](const int cellLabel)
    {
        if (cellLabel < 0)
            return;
        scratch.descend.push_back(cellLabel);
        while (!scratch.descend.empty())
        {
            const int cur = scratch.descend.back();
            scratch.descend.pop_back();
            if (cur < 0)
                continue;
            if (visitedContains(cur))
                continue;

            scratch.visited.push_back(cur);

            const Cell& cell = mesh.cells(cur);
            if (!cellIntersectsAABB(cell, box))
                continue;

            if (isLeafCell(cell))
            {
                stencil.push_back(cur);
                scratch.queue.push_back(cur);
                continue;
            }

            // Not a leaf: descend into children that may intersect.
            for (int k = 0; k < 4; ++k)
            {
                const int son = cell.son[k];
                if (son != -1)
                    scratch.descend.push_back(son);
            }
        }
    };

    pushCellOrDescend(hostCell);

    // ----------------------------------------------------------------------
    // Fast path for the common case hp <= hHost (CIC with hp tied to host cell size).
    //
    // If hp is no larger than the host cell width, the particle support
    // box extends at most one cell width in each direction. In that situation,
    // the support can only overlap:
    //   - the host leaf cell, and
    //   - its immediate face neighbors, and
    //   - at most one "ring" of corner neighbors.
    //
    // This avoids a full BFS expansion (which is correct but more expensive).
    // We keep the generic BFS as a fallback for larger hp.
    // ----------------------------------------------------------------------
    bool usedFastRing = false;
    if (hostCell >= 0)
    {
        const Cell& host = mesh.cells(hostCell);
        if (isLeafCell(host))
        {
            const Vector2D hs = host.size();
            const double hHost = std::min(hs.x, hs.y);
            // Conservative guard: only use ring enumeration when hp <= hHost.
            if (hp <= hHost * (1.0 + 1e-12))
            {
                usedFastRing = true;

                // Collect immediate face neighbors of the host.
                std::vector<int> faceNbrs;
                faceNbrs.reserve(16);

                for (int dir = 0; dir < 4; ++dir)
                {
                    for (int slot = 0; slot < 2; ++slot)
                    {
                        const int nb = host.neiCell[dir][slot];
                        if (nb != -1)
                        {
                            faceNbrs.push_back(nb);
                            pushCellOrDescend(nb);
                        }
                    }
                }

                // One ring of corners: take perpendicular neighbors of each
                // face neighbor (two directions). This captures diagonal
                // overlaps without a full BFS.
                for (const int nb : faceNbrs)
                {
                    if (nb < 0)
                        continue;
                    const Cell& c = mesh.cells(nb);

                    // For each neighbor, visit left/right and up/down.
                    for (int perpDir = 0; perpDir < 4; ++perpDir)
                    {
                        // Skip the direction that points back to host to reduce
                        // redundant candidates (visited set still guarantees correctness).
                        // We still iterate all dirs for simplicity; the overhead
                        // is tiny compared to the old BFS.
                        for (int slot = 0; slot < 2; ++slot)
                        {
                            const int nb2 = c.neiCell[perpDir][slot];
                            if (nb2 != -1)
                                pushCellOrDescend(nb2);
                        }
                    }
                }
            }
        }
    }

    if (!usedFastRing)
    {
        // Generic BFS from intersecting leaf cells to cover the whole support region.
        std::size_t qHead = 0;
        while (qHead < scratch.queue.size())
        {
            const int cLabel = scratch.queue[qHead++];
            const Cell& c = mesh.cells(cLabel);

            // 4 directions, up to 2 neighbors per direction (coarse-fine interfaces).
            for (int dir = 0; dir < 4; ++dir)
            {
                for (int slot = 0; slot < 2; ++slot)
                {
                    const int nb = c.neiCell[dir][slot];
                    if (nb != -1)
                        pushCellOrDescend(nb);
                }
            }
        }
    }

    // Deterministic ordering is helpful for debugging/reproducibility.
    if (deterministic)
        std::sort(stencil.begin(), stencil.end());

    if (debugPrint && ConsoleControl::enabled(ConsoleControl::Level::kTrace))
    {
        std::cout << "[StencilQuery] xp=" << xp
                  << " yp=" << yp
                  << " hp=" << hp
                  << " bbox=[" << box.xmin << "," << box.xmax
                  << "]x[" << box.ymin << "," << box.ymax << "]"
                  << " stencilCells=" << stencil.size()
                  << (usedFastRing ? " (fast-ring)" : " (bfs)")
                  << std::endl;
    }

}

std::vector<int> ParticleStencilQuery::enumerateLeafCellsIntersectParticleSupport(
    const BackgroundMesh& mesh,
    const double xp,
    const double yp,
    const double hp,
    const bool debugPrint,
    const int hostCellIn,
    const bool deterministic)
{
    std::vector<int> stencil;
    enumerateLeafCellsIntersectParticleSupport(mesh, xp, yp, hp, &stencil, debugPrint, hostCellIn, deterministic);
    return stencil;
}

